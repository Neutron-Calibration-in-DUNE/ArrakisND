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
    Blip
)
from arrakis_nd.arrakis.common import blip_data_type


class NeutronPlugin(Plugin):
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
        super(NeutronPlugin, self).__init__(config)

        self.input_products = [
            'daughters',
            'track_id_hit_map',
            'track_id_hit_segment_map',
            'track_id_hit_t0_map',
            'parent_pdg_id',
        ]
        self.output_products = ['tracklette', 'track', 'blip']

        self.neutron_labels = {
            'topology': Topology.Blip.value,
            'physics': Physics.NuclearRecoil.value
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
        daughters = event_products['daughters']

        trajectories_traj_ids = trajectories['traj_id']
        trajectories_vertex_ids = trajectories['vertex_id']
        trajectories_parent_ids = trajectories['parent_id']
        trajectories_pdg_ids = trajectories['pdg_id']
        trajectories_xyz_start = trajectories['xyz_start']
        trajectories_xyz_end = trajectories['xyz_end']
        trajectories_pxyz_start = trajectories['pxyz_start']
        trajectories_pxyz_end = trajectories['pxyz_end']
        trajectories_E = trajectories['E_start']
        trajectories_start_subprocess = trajectories['start_subprocess']
        trajectories_end_subprocess = trajectories['end_subprocess']
        charge_x = charge['x']
        charge_y = charge['y']
        charge_z = charge['z']
        charge_E = charge['E']
        charge_io_group = charge['io_group']

        """Grab the neutron which result in ionization"""
        particle_mask = (
            (abs(trajectories_pdg_ids) == 2112)
        )

        """Iterate over the individual (particle_id, vertex_id, traj_index)"""
        """
        First, we generate blip objects associated to the neutron itself.
        These will include elastic/inelastic scatters, but not neutron
        captures.  Those will be handled in a second loop.
        """
        for ii, (particle_id, vertex_id, traj_index) in enumerate(zip(
            trajectories_traj_ids[particle_mask],
            trajectories_vertex_ids[particle_mask],
            event_indices['trajectories'][particle_mask]
        )):
            """Get indices in charge_segments arrays"""
            particle_hits = track_id_hit_map[(particle_id, vertex_id)]

            """Check if neutron has no hits, go to the next particle"""
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

            """Iterate over standard labels"""
            for label, value in self.neutron_labels.items():
                arrakis_charge[label][particle_hit_segments] = value

            """Set the particle label"""
            arrakis_charge['particle'][particle_hit_segments] = trajectories_pdg_ids[particle_mask][ii]

            """Generate unique label"""
            arrakis_charge['unique_topology'][particle_hit_segments] = traj_index

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
                io_group_t0s = particle_hit_t0s[io_group_mask]

                """Find the closest points to track begin/end for this io_group"""
                io_group_start_distances = particle_start_distances[io_group_mask]
                io_group_end_distances = particle_end_distances[io_group_mask]
                blip_start_index = np.argmin(io_group_start_distances)
                blip_end_index = np.argmin(io_group_end_distances)

                """Parameterize the trajectory of this blip using t0"""
                io_group_xyz = particle_charge_xyz[io_group_mask]
                track_data = fit_track(io_group_t0s, io_group_xyz)

                """Now generate the associated blip CAF object"""
                """
                The blip data object has the following entries
                that must be filled (see arrakis.common).
                """
                io_group_xyz = particle_charge_xyz[io_group_mask]

                """Set blip type"""
                blip_type = Blip.NuclearRecoil.value

                """Create the CAF object for this blip"""
                blip_data = np.array([(
                    event,
                    io_group,
                    particle_id,
                    vertex_id,
                    ancestor_id,
                    blip_type,
                    particle_xyz_start,
                    particle_xyz_end,
                    [io_group_xyz[blip_start_index]],
                    [io_group_xyz[blip_end_index]],
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

        """Grab the neutrons which result in capture"""
        particle_mask = (
            (abs(trajectories_pdg_ids) == 2112) &
            (
                (abs(trajectories_end_subprocess) == SubProcessType.HadronCapture.value) |
                (abs(trajectories_end_subprocess) == SubProcessType.HadronCaptureAtRest.value)
            )
        )

        """
        Now we iterate over the neutron captures.
        """
        for ii, (particle_id, vertex_id, traj_index) in enumerate(zip(
            trajectories_traj_ids[particle_mask],
            trajectories_vertex_ids[particle_mask],
            event_indices['trajectories'][particle_mask]
        )):
            """Grab gamma daughters"""
            daughter_ids = daughters[particle_mask][ii]
            neutron_daughter_traj_ids = trajectories_traj_ids[daughter_ids]
            neutron_daughter_pdg_ids = trajectories_pdg_ids[daughter_ids]
            neutron_daughter_start_subprocess = trajectories_start_subprocess[daughter_ids]

            """Iterate over gammas"""
            for jj, daughter_id in enumerate(neutron_daughter_traj_ids):
                if abs(neutron_daughter_pdg_ids[jj]) != 22:
                    continue
                if (
                    (abs(neutron_daughter_start_subprocess[jj]) != SubProcessType.HadronCapture.value) &
                    (abs(neutron_daughter_start_subprocess[jj]) != SubProcessType.HadronCaptureAtRest.value)
                ):
                    continue

                """Get indices in charge_segments arrays"""
                gamma_hits = track_id_hit_map[(daughter_id, vertex_id)]

                gamma_daughter_traj_ids = trajectories_traj_ids[
                    (trajectories_parent_ids == daughter_id) & (trajectories_vertex_ids == vertex_id)
                ]
                for kk, gamma_daughter in enumerate(gamma_daughter_traj_ids):
                    gamma_daughter_hits = track_id_hit_map[(gamma_daughter, vertex_id)]
                    gamma_hits += gamma_daughter_hits

