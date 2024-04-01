"""
"""
import h5py
import numpy as np
from scipy.interpolate import splprep, splev
import warnings
from scipy.integrate import quad, IntegrationWarning

from arrakis_nd.utils.utils import profiler, integrand
from arrakis_nd.plugins.plugin import Plugin
from arrakis_nd.dataset.common import (
    ProcessType, SubProcessType,
    Topology, Physics
)


class ElectronPlugin(Plugin):
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
        super(ElectronPlugin, self).__init__(config)

        self.input_products = [
            'daughters',
            'track_id_hit_map',
            'track_id_hit_segment_map',
            'track_id_hit_t0_map',
            'parent_pdg_id',
        ]
        self.output_products = ['fragment', 'blip']

        self.compton_labels = {
            'topology': Topology.Shower.value,
            'physics': Physics.GammaCompton.value,
        }
        self.conversion_labels = {
            'topology': Topology.Shower.value,
            'physics': Physics.GammaConversion.value,
        }
        self.low_energy_labels = {
            'topology': Topology.Blip.value,
            'physics': Physics.ElectronIonization.value,
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
        parent_pdg_ids = event_products['parent_pdg_id']

        trajectories_traj_ids = trajectories['traj_id']
        trajectories_vertex_ids = trajectories['vertex_id']
        trajectories_pdg_ids = trajectories['pdg_id']
        trajectories_parent_ids = trajectories['parent_id']
        trajectories_xyz_start = trajectories['xyz_start']
        trajectories_xyz_end = trajectories['xyz_end']
        trajectories_pxyz_start = trajectories['pxyz_start']
        trajectories_pxyz_end = trajectories['pxyz_end']
        trajectories_E = trajectories['E_start']
        trajectories_start_process = trajectories['start_process']
        trajectories_start_subprocess = trajectories['start_subprocess']
        charge_x = charge['x']
        charge_y = charge['y']
        charge_z = charge['z']
        charge_E = charge['E']
        charge_io_group = charge['io_group']

        """Grab the electrons which are from compton scatters, conversions or electron ionization"""
        particle_mask = (
            (abs(trajectories_pdg_ids) == 11) &
            (
                (
                    (abs(trajectories_start_subprocess) == SubProcessType.ComptonScattering.value) |
                    (abs(trajectories_start_subprocess) == SubProcessType.GammaConversion.value)
                ) |
                (
                    (abs(parent_pdg_ids) == 11) &
                    (abs(trajectories_start_subprocess) == SubProcessType.Ionization.value)
                )
            )
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

            """############################### Electrons ###############################"""
            """
            
            """
            particle_charge_xyz = np.array([
                charge_x[particle_hits],
                charge_y[particle_hits],
                charge_z[particle_hits]
            ]).T
            particle_charge_E = charge_E[particle_hits]
            particle_xyz_start = trajectories_xyz_start[particle_mask][ii]
            particle_xyz_end = trajectories_xyz_end[particle_mask][ii]
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

            if closest_start_index != closest_end_index:
                """Add delta/michel vertex"""
                arrakis_charge['vertex'][particle_start_index] = 1
                arrakis_charge['tracklette_end'][particle_end_index] = 1

            """If this particle is a Delta, add delta labels and vertex"""
            if (
                (abs(trajectories_start_subprocess[particle_mask][ii]) == SubProcessType.ComptonScattering.value)
            ):
                """Iterate over standard labels"""
                for label, value in self.compton_labels.items():
                    arrakis_charge[label][particle_hit_segments] = value
            elif (
                (abs(trajectories_start_subprocess[particle_mask][ii]) == SubProcessType.GammaConversion.value)
            ):
                """If this particle is a Michel electron, add labels and CAF"""
                for label, value in self.conversion_labels.items():
                    arrakis_charge[label][particle_hit_segments] = value
            else:
                """
                If the electrons parent is a MIP/HIP, but has a process and subprocess that
                corresponds to an intermediate particle, such as a photon, then we have to
                be more careful with the logic. Edep-sim can be somewhat difficult to deal
                with since it hides information about intermediate photons. 
                """
                """Iterate over standard labels"""
                for label, value in self.low_energy_labels.items():
                    arrakis_charge[label][particle_hit_segments] = value

            """Set the particle label"""
            arrakis_charge['particle'][particle_hit_segments] = trajectories_pdg_ids[particle_mask][ii]

            """Generate unique label"""
            arrakis_charge['unique_topology'][particle_hit_segments] = traj_index

        """Write changes to arrakis_file"""
        arrakis_file['charge_segment/calib_final_hits/data'][event_indices['charge']] = arrakis_charge
