"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
from arrakis_nd.utils.track_utils import fit_track
from arrakis_nd.plugins.plugin import Plugin
from arrakis_nd.dataset.common import (
    Topology, Physics,
    Track
)
from arrakis_nd.arrakis.common import (
    tracklette_data_type,
    track_data_type
)


class TrackPlugin(Plugin):
    """
    A plugin for labeling mips in charge and light data.

    mips are MIPs which deposit energy through MIP ionization (i.e. the
    standard dE/dx mechanism).  This plugin grabs all mips in an event
    and labels their energy deposits according to the following scheme:

        (1) hits coming directly from mip ionization are labeled as
                topology = track
                physics = mip
            mips are defined as muons, tauons, pions and kaons
        (2) hits coming directly from hip ionization are labeled as
                topology = track
                physics = hip
            hips are defined as protons

    Apart from these labels, track beginning and ending points are also
    labelled.

    StandardRecord objects are also generated.
    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        """
        super(TrackPlugin, self).__init__(config)

        self.input_products = [
            'daughters',
            'track_id_hit_map',
            'track_id_hit_segment_map',
            'track_id_hit_t0_map'
        ]
        self.output_products = ['tracklette', 'track']

        self.mip_labels = {
            'topology': Topology.Track.value,
            'physics': Physics.MIP.value,
        }
        self.hip_labels = {
            'topology': Topology.Track.value,
            'physics': Physics.HIP.value,
        }

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
        trajectories = flow_file['mc_truth/trajectories/data'][event_indices['trajectories']]
        charge = flow_file['charge/calib_final_hits/data'][event_indices['charge']]
        arrakis_charge = arrakis_file['charge_segment/calib_final_hits/data'][event_indices['charge']]
        track_id_hit_map = event_products['track_id_hit_map']
        track_id_hit_segment_map = event_products['track_id_hit_segment_map']
        track_id_hit_t0_map = event_products['track_id_hit_t0_map']

        trajectories_traj_ids = trajectories['traj_id']
        trajectories_vertex_ids = trajectories['vertex_id']
        trajectories_pdg_ids = trajectories['pdg_id']
        trajectories_xyz_start = trajectories['xyz_start']
        trajectories_xyz_end = trajectories['xyz_end']
        trajectories_pxyz_start = trajectories['pxyz_start']
        trajectories_pxyz_end = trajectories['pxyz_end']
        trajectories_E = trajectories['E_start']
        charge_x = charge['x']
        charge_y = charge['y']
        charge_z = charge['z']
        charge_E = charge['E']
        charge_io_group = charge['io_group']

        """Grab the mips/hips"""
        particle_mask = (
            (abs(trajectories_pdg_ids) == 13) |
            (abs(trajectories_pdg_ids) == 15) |
            (abs(trajectories_pdg_ids) == 211) |
            (abs(trajectories_pdg_ids) == 321) |
            (abs(trajectories_pdg_ids) == 2212)
        )

        """Iterate over the individual (particle_id, vertex_id, traj_index)"""
        for ii, (particle_id, vertex_id, traj_index) in enumerate(zip(
            trajectories_traj_ids[particle_mask],
            trajectories_vertex_ids[particle_mask],
            event_indices['trajectories'][particle_mask]
        )):
            """Get indices in charge_segments arrays"""
            particle_hits = track_id_hit_map[(particle_id, vertex_id)]

            """Check if mip has no hits, go to the next particle"""
            if not any(particle_hits):
                continue

            particle_segments = track_id_hit_segment_map[(particle_id, vertex_id)]
            particle_hit_segments = (
                particle_hits, particle_segments
            )

            """Get the associated t0 values"""
            particle_hit_t0s = track_id_hit_t0_map[(particle_id, vertex_id)]

            """############################### Set standard labels ###############################"""
            """Set the event_id"""
            arrakis_charge['event_id'][particle_hits] = event

            """Iterate over standard labels"""
            if abs(trajectories_pdg_ids[particle_mask][ii]) == 2212:
                track_type = Track.HIP.value
                for label, value in self.hip_labels.items():
                    arrakis_charge[label][particle_hit_segments] = value
            else:
                track_type = Track.MIP.value
                for label, value in self.mip_labels.items():
                    arrakis_charge[label][particle_hit_segments] = value

            """Set the particle label"""
            arrakis_charge['particle'][particle_hit_segments] = trajectories_pdg_ids[particle_mask][ii]

            """Generate unique label"""
            arrakis_charge['unique_topology'][particle_hit_segments] = traj_index

            """Determine track begin and end points"""
            particle_charge_xyz = np.array([
                charge_x[particle_hits],
                charge_y[particle_hits],
                charge_z[particle_hits]
            ]).T
            particle_charge_E = charge_E[particle_hits]
            particle_xyz_start = trajectories_xyz_start[particle_mask][ii]
            particle_xyz_end = trajectories_xyz_end[particle_mask][ii]
            particle_pxyz_start = trajectories_pxyz_start[particle_mask][ii]
            particle_pxyz_end = trajectories_pxyz_end[particle_mask][ii]
            particle_E = trajectories_E[particle_mask][ii]

            """
            To do this we find the closest hits to the actual track beginning
            and ending from the trajectories.  If these points end up being the
            same, we skip this track, since there's no way to distinguish
            that point as a track anyways.
            """
            particle_start_distances = np.sqrt(
                ((particle_charge_xyz - particle_xyz_start) ** 2).sum(axis=1)
            )
            particle_end_distances = np.sqrt(
                ((particle_charge_xyz - particle_xyz_end) ** 2).sum(axis=1)
            )
            closest_start_index = np.argmin(particle_start_distances)
            closest_end_index = np.argmin(particle_end_distances)

            particle_start_index = (
                particle_hits[closest_start_index],
                particle_segments[closest_start_index]
            )
            particle_end_index = (
                particle_hits[closest_end_index],
                particle_segments[closest_end_index]
            )

            if closest_start_index == closest_end_index:
                continue

            """Set track beginning and ending points"""
            arrakis_charge['tracklette_begin'][particle_start_index] = 1
            arrakis_charge['tracklette_end'][particle_end_index] = 1

            """Parameterize the trajectory of this track using t0"""
            track_data = fit_track(particle_hit_t0s, particle_charge_xyz)

            """############################### Construct Tracklettes ###############################"""

            """
            Now we must determine if the track is broken along its
            trajectory.  We can do this by checking if the beginning and
            ending points are in different TPCs.  If so, then we will have
            to determine the beginning and ending points of each tracklette.

            We can do this by isolating hits according to the io_group,
            which defines the particular TPC region that a hit is in.  For
            each unique io_group associated to the hits of this particle,
            we will follow the same procedure and find the closest
            points to the track beginning and ending to get the tracklette
            begin and end points.
            """
            particle_io_group = charge_io_group[particle_hits]
            unique_particle_io_group = np.unique(particle_io_group)
            particle_hits = np.array(particle_hits)
            particle_segments = np.array(particle_segments)
            particle_hit_t0s = np.array(particle_hit_t0s)
            for io_group in unique_particle_io_group:
                """Isolate hits from this io_group"""
                io_group_mask = (particle_io_group == io_group)
                if not sum(io_group_mask):
                    continue
                io_group_hits = particle_hits[io_group_mask]
                io_group_segments = particle_segments[io_group_mask]
                io_group_t0s = particle_hit_t0s[io_group_mask]

                """Find the closest points to track begin/end for this io_group"""
                io_group_start_distances = particle_start_distances[io_group_mask]
                io_group_end_distances = particle_end_distances[io_group_mask]
                tracklette_start_index = np.argmin(io_group_start_distances)
                tracklette_end_index = np.argmin(io_group_end_distances)

                io_group_start_index = (
                    io_group_hits[tracklette_start_index],
                    io_group_segments[tracklette_start_index]
                )
                io_group_end_index = (
                    io_group_hits[tracklette_end_index],
                    io_group_segments[tracklette_end_index]
                )

                if tracklette_start_index == tracklette_end_index:
                    continue

                arrakis_charge['tracklette_begin'][io_group_start_index] = 1
                arrakis_charge['tracklette_end'][io_group_end_index] = 1

                """Parameterize the trajectory of this tracklette using t0"""
                io_group_xyz = particle_charge_xyz[io_group_mask]
                track_data = fit_track(io_group_t0s, io_group_xyz)

                """Now generate the associated tracklette CAF object"""
                """
                The tracklette data object has the following entries
                that must be filled (see arrakis.common).
                """
                io_group_xyz = particle_charge_xyz[io_group_mask]

                """Create the CAF object for this tracklette"""
                tracklette_data = np.array([(
                    event,
                    io_group,
                    particle_id,
                    vertex_id,
                    track_type,
                    [0, 0, 0],
                    [0, 0, 0],
                    [io_group_xyz[tracklette_start_index]],
                    [io_group_xyz[tracklette_end_index]],
                    [0, 0, 0],
                    [0, 0, 0],
                    [track_data['track_dir']],
                    [track_data['track_enddir']],
                    sum(particle_charge_E[io_group_mask]),
                    1,
                    track_data['track_len_gcm2'],
                    track_data['track_len_cm'],
                    0,
                    [[traj_index]+[0]*19],
                    [[1]+[0]*19])],
                    dtype=tracklette_data_type
                )

                """Add the new tracklette to the CAF objects"""
                event_products['tracklette'].append(tracklette_data)

            """############################### Construct Track ###############################"""
            """Now generate the associated track CAF object"""
            """
            The track data object has the following entries
            that must be filled (see arrakis.common).
            """
            io_group_xyz = particle_charge_xyz[io_group_mask]

            """Create the CAF object for this track"""
            particle_pxyz_start_magnitude = np.linalg.norm(particle_pxyz_start)
            particle_pxyz_end_magnitude = np.linalg.norm(particle_pxyz_end)
            if particle_pxyz_start_magnitude > 0:
                particle_pxyz_start /= particle_pxyz_start_magnitude
            if particle_pxyz_end_magnitude > 0:
                particle_pxyz_end /= particle_pxyz_end_magnitude

            track_data = np.array([(
                event,
                particle_id,
                vertex_id,
                track_type,
                particle_xyz_start,
                particle_xyz_end,
                [particle_charge_xyz[closest_start_index]],
                [particle_charge_xyz[closest_end_index]],
                particle_pxyz_start,
                particle_pxyz_end,
                [track_data['track_dir']],
                [track_data['track_enddir']],
                sum(particle_charge_E),
                1,
                track_data['track_len_gcm2'],
                track_data['track_len_cm'],
                particle_E,
                [[traj_index]+[0]*19],
                [[1]+[0]*19])],
                dtype=track_data_type
            )

            """Add the new track to the CAF objects"""
            event_products['track'].append(track_data)

        """Write changes to arrakis_file"""
        arrakis_file['charge_segment/calib_final_hits/data'][event_indices['charge']] = arrakis_charge
