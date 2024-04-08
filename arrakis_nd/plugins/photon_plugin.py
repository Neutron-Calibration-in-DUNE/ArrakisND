"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
from arrakis_nd.plugins.plugin import Plugin
from arrakis_nd.utils.track_utils import fit_track
from arrakis_nd.dataset.common import (
    SubProcessType,
    Topology, Physics,
    Fragment
)
from arrakis_nd.arrakis.common import (
    fragment_data_type,
    blip_data_type
)


class PhotonPlugin(Plugin):
    """_summary_

    Args:
        Plugin (_type_): _description_
    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        """
        super(PhotonPlugin, self).__init__(config)

        self.input_products = [
            'daughters',
            'track_id_hit_map',
            'track_id_hit_segment_map',
            'track_id_hit_t0_map',
            'parent_pdg_id',
        ]
        self.output_products = ['fragment', 'blip']

        if "shower_size_threshold" not in self.config:
            self.config["shower_size_threshold"] = 10
        self.shower_size_threshold = self.config["shower_size_threshold"]

        self.photon_shower_labels = {
            'topology': Topology.Shower.value,
            'physics': Physics.GammaCompton.value,
        }
        self.photon_blip_labels = {
            'topology': Topology.Blip.value,
            'physics': Physics.GammaCompton.value
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
        ancestor_traj_id_map = event_products['ancestor_traj_id_map']

        trajectories_traj_ids = trajectories['traj_id']
        trajectories_vertex_ids = trajectories['vertex_id']
        trajectories_pdg_ids = trajectories['pdg_id']
        trajectories_xyz_start = trajectories['xyz_start']
        trajectories_xyz_end = trajectories['xyz_end']
        trajectories_pxyz_start = trajectories['pxyz_start']
        trajectories_pxyz_end = trajectories['pxyz_end']
        trajectories_start_subprocess = trajectories['start_subprocess']
        trajectories_E = trajectories['E_start']
        charge_x = charge['x']
        charge_y = charge['y']
        charge_z = charge['z']
        charge_E = charge['E']
        charge_io_group = charge['io_group']

        """Grab the photons"""
        particle_mask = (
            (abs(trajectories_pdg_ids) == 22)
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

            """Get ancestor traj_id"""
            ancestor_id = ancestor_traj_id_map[(particle_id, vertex_id)]

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
            """Determine track begin and end points"""
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

            """############################### Set standard labels ###############################"""
            """Set the event_id"""
            arrakis_charge['event_id'][particle_hits] = event

            """Variable determining if this particle made a fragment or blip"""
            fragment_type = Fragment.Undefined.value
            if (trajectories_start_subprocess[particle_mask][ii] == SubProcessType.Bremsstrahlung.value):
                fragment_type = Fragment.Bremsstrahlung.value
            elif (trajectories_start_subprocess[particle_mask][ii] == SubProcessType.Annihilation.value):
                fragment_type = Fragment.Annihilation.value
            else:
                fragment_type = Fragment.GammaCompton.value
                """This is either neutron capture or hadron inelastic"""

            """############################### Photons ###############################"""
            is_fragment = True
            if len(particle_hits) >= self.shower_size_threshold:
                for label, value in self.photon_shower_labels.items():
                    arrakis_charge[label][particle_hit_segments] = value
            else:
                is_fragment = False
                for label, value in self.photon_blip_labels.items():
                    arrakis_charge[label][particle_hit_segments] = value

            """Set the particle label"""
            arrakis_charge['particle'][particle_hit_segments] = trajectories_pdg_ids[particle_mask][ii]

            """Generate unique label"""
            arrakis_charge['unique_topology'][particle_hit_segments] = traj_index

            """############################### Construct Fragments ###############################"""

            """
            Now we must determine if the fragment is broken along its
            trajectory.  We can do this by checking if the beginning and
            ending points are in different TPCs.  If so, then we will have
            to determine the beginning and ending points of each fragment.

            We can do this by isolating hits according to the io_group,
            which defines the particular TPC region that a hit is in.  For
            each unique io_group associated to the hits of this particle,
            we will follow the same procedure and find the closest
            points to the track beginning and ending to get the fragment
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
                fragment_start_index = np.argmin(io_group_start_distances)
                fragment_end_index = np.argmin(io_group_end_distances)

                io_group_start_index = (
                    io_group_hits[fragment_start_index],
                    io_group_segments[fragment_start_index]
                )
                io_group_end_index = (
                    io_group_hits[fragment_end_index],
                    io_group_segments[fragment_end_index]
                )

                if (fragment_start_index != fragment_end_index) and (is_fragment):
                    arrakis_charge['fragment_begin'][io_group_start_index] = 1
                    arrakis_charge['fragment_end'][io_group_end_index] = 1

                """Parameterize the trajectory of this fragment using t0"""
                io_group_xyz = particle_charge_xyz[io_group_mask]
                track_data = fit_track(io_group_t0s, io_group_xyz)

                """Now generate the associated fragment CAF object"""
                """
                The fragment data object has the following entries
                that must be filled (see arrakis.common).
                """
                io_group_xyz = particle_charge_xyz[io_group_mask]

                if is_fragment:
                    """Create the CAF object for this fragment"""
                    fragment_data = np.array([(
                        event,
                        io_group,
                        particle_id,
                        vertex_id,
                        ancestor_id,
                        fragment_type,
                        particle_xyz_start,
                        particle_xyz_end,
                        [io_group_xyz[fragment_start_index]],
                        [io_group_xyz[fragment_end_index]],
                        particle_pxyz_start,
                        particle_pxyz_end,
                        [track_data['track_dir']],
                        [track_data['track_enddir']],
                        sum(particle_charge_E[io_group_mask]),
                        1,
                        track_data['track_len_gcm2'],
                        track_data['track_len_cm'],
                        particle_E,
                        [[traj_index]+[0]*19],
                        [[1]+[0]*19])],
                        dtype=fragment_data_type
                    )

                    """Add the new fragment to the CAF objects"""
                    event_products['fragment'].append(fragment_data)
                else:
                    """Create the CAF object for this fragment"""
                    blip_data = np.array([(
                        event,
                        io_group,
                        particle_id,
                        vertex_id,
                        ancestor_id,
                        fragment_type,
                        particle_xyz_start,
                        particle_xyz_end,
                        [io_group_xyz[fragment_start_index]],
                        [io_group_xyz[fragment_end_index]],
                        particle_pxyz_start,
                        particle_pxyz_end,
                        [track_data['track_dir']],
                        [track_data['track_enddir']],
                        sum(particle_charge_E[io_group_mask]),
                        1,
                        track_data['track_len_gcm2'],
                        track_data['track_len_cm'],
                        particle_E,
                        [[traj_index]+[0]*19],
                        [[1]+[0]*19])],
                        dtype=blip_data_type
                    )

                    """Add the new blip to the CAF objects"""
                    event_products['blip'].append(blip_data)

        """Write changes to arrakis_file"""
        arrakis_file['charge_segment/calib_final_hits/data'][event_indices['charge']] = arrakis_charge
