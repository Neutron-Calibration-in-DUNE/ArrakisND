"""
"""
import h5py
import numpy as np
from scipy.interpolate import splprep, splev

from arrakis_nd.utils.utils import profiler
from arrakis_nd.utils.utils import fiducialized_vertex
from arrakis_nd.utils.track_utils import fit_track
from arrakis_nd.plugins.plugin import Plugin
from arrakis_nd.arrakis.common import particle_data_type


class ParticlePlugin(Plugin):
    """
    A plugin for creating particle standard record objects.

    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        A particle standard record object is essentially a fully reconstructed
        particle that can be associated back to the particles in the trajectories
        within mc_truth.  Not all information is (or can be) reconstructed, such
        as the process, endprocess, etc., but basic information about the start
        and end points, the starting energy and momentum and whether the particle
        is contained within the TPC are present.

        Constructing these objects is fairly simple, we just need to loop over
        all particles in the trajectories for this event, and see if each particle
        has reconstructed hits or not. If so, then construct the associated
        particle standard record object.

        Since it is in principle impossible to determine the "actual" start and
        end points of particle trajectories, we instead identify the hits
        which are closest to the start and end points.  Even if one could pin point
        the exact starting location for trajectories which begin inside the detector,
        there is no way to determine this for particles outside.  Likewise for the momentum,
        since particles which begin outside cannot have their momentum reconstructed
        faithfully, we instead determine the momentum from fitting a curve to the
        trajectory of the particle and estimating the gradient at the start point.

        Another quantity which seems rather useless is 'tgtA', which is the target
        nucleus that the particle originated from.  This is also impossible to determine
        outside of the detector, so only argon is a reasonable quantity to use, hence
        tgtA is redundant.
        """
        super(ParticlePlugin, self).__init__(config)

        self.input_products = [
            'daughters',
            'track_id_hit_map',
            'track_id_hit_segment_map',
            'track_id_hit_t0_map'
        ]
        self.output_products = ['particle']

    @profiler
    def process_event(
        self,
        event: int,
        flow_file: h5py.File,
        arrakis_file: h5py.File,
        event_indices: dict,
        event_products: dict,
    ):
        """
        """
        interactions = flow_file['mc_truth/interactions/data'][event_indices['interactions']]
        trajectories = flow_file['mc_truth/trajectories/data'][event_indices['trajectories']]
        charge = flow_file['charge/calib_final_hits/data'][event_indices['charge']]
        track_id_hit_map = event_products['track_id_hit_map']
        track_id_hit_t0_map = event_products['track_id_hit_t0_map']

        trajectories_traj_ids = trajectories['traj_id']
        trajectories_vertex_ids = trajectories['vertex_id']
        trajectories_pdg_ids = trajectories['pdg_id']
        trajectories_parent_ids = trajectories['parent_id']
        trajectories_xyz_start = trajectories['xyz_start']
        trajectories_xyz_end = trajectories['xyz_end']
        trajectories_pxyz_start = trajectories['pxyz_start']
        trajectories_E = trajectories['E_start']
        charge_x = charge['x']
        charge_y = charge['y']
        charge_z = charge['z']
        charge_E = charge['E']

        """Iterate over the individual (traj_id, vertex_id, traj_index)"""
        for ii, (particle_id, vertex_id, traj_index) in enumerate(zip(
            trajectories_traj_ids,
            trajectories_vertex_ids,
            event_indices['trajectories']
        )):
            """Get indices in charge_segments arrays"""
            particle_hits = track_id_hit_map[(particle_id, vertex_id)]

            """Check if particle has no hits, go to the next particle"""
            if not any(particle_hits):
                continue

            """Get the associated t0 values"""
            particle_hit_t0s = track_id_hit_t0_map[(particle_id, vertex_id)]

            """Determine if the particle is a primary"""
            if trajectories_parent_ids[ii] == -1:
                particle_primary = 1
            else:
                particle_primary = 0

            """Find the closest start and end hits"""
            particle_charge_xyz = np.array([
                charge_x[particle_hits],
                charge_y[particle_hits],
                charge_z[particle_hits]
            ]).T
            particle_charge_E = charge_E[particle_hits]
            particle_xyz_start = trajectories_xyz_start[ii]
            particle_xyz_end = trajectories_xyz_end[ii]
            particle_pxyz_start = trajectories_pxyz_start[ii]
            particle_E = trajectories_E[ii]

            """
            To do this we find the closest hits to the actual particle beginning
            and ending from the trajectories.  If these points end up being the
            same, we skip this particle, since there's no way to distinguish
            that point as a particle anyways.
            """
            particle_start_distances = np.sqrt(
                ((particle_charge_xyz - particle_xyz_start) ** 2).sum(axis=1)
            )
            particle_end_distances = np.sqrt(
                ((particle_charge_xyz - particle_xyz_end) ** 2).sum(axis=1)
            )
            closest_start_index = np.argmin(particle_start_distances)
            closest_end_index = np.argmin(particle_end_distances)

            closest_particle_xyz_start = particle_charge_xyz[closest_start_index]
            closest_particle_xyz_end = particle_charge_xyz[closest_end_index]

            """Parameterize the trajectory of this particle using t0"""
            track_data = fit_track(particle_hit_t0s, particle_charge_xyz)

            """Determine tgtA"""
            interaction = interactions[(interactions['vertex_id'] == vertex_id)]
            if fiducialized_vertex(particle_xyz_start):
                tgtA = 18
            else:
                if particle_primary:
                    tgtA = interaction['target']
                else:
                    tgtA = 0

            """Let's create the particle object"""
            """
            Particles in the standard record have the following variables:
                particle_data_type = np.dtype([
                    ('event_id', 'i4'),
                    ('particle_id', 'i4'),
                    ('vertex_id', 'i8'),
                    ('primary', 'i4'),
                    ('pdg', 'i4'),
                    ('tgtA', 'i4'),
                    ('score', 'f4'),
                    ('Evis', 'f4'),
                    ('E', 'f4'),
                    ('E_method', 'i4'),
                    ('p', 'f4', (1, 3)),
                    ('dir_hit', 'f4', (1, 3)),
                    ('start', 'f4', (1, 3)),
                    ('end', 'f4', (1, 3)),
                    ('start_hit', 'f4', (1, 3)),
                    ('end_hit', 'f4', (1, 3)),
                    ('contained', 'i4'),
                    ('truth', 'i4', (1, 20)),
                    ('truthOverlap', 'f4', (1, 20)),
                ])
            """

            particle = np.array([(
                event,
                particle_id,
                vertex_id,
                particle_primary,
                trajectories_pdg_ids[ii],
                tgtA,
                1,
                sum(particle_charge_E),
                particle_E,
                0,
                particle_pxyz_start,
                track_data['track_dir'],
                particle_xyz_start,
                particle_xyz_end,
                closest_particle_xyz_start,
                closest_particle_xyz_end,
                fiducialized_vertex(particle_xyz_start) and fiducialized_vertex(particle_xyz_end),
                [[traj_index]+[0]*19],
                [[1]+[0]*19])],
                dtype=particle_data_type
            )

            """Add new particle to the event products"""
            event_products['particle'].append(particle)
