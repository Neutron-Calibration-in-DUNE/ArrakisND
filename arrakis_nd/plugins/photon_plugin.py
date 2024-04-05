"""
"""
import h5py

from arrakis_nd.utils.utils import profiler
from arrakis_nd.plugins.plugin import Plugin
from arrakis_nd.dataset.common import (
    Topology, Physics
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
        arrakis_charge = arrakis_file['charge_segment/calib_final_hits/data'][event_indices['charge']]
        track_id_hit_map = event_products['track_id_hit_map']
        track_id_hit_segment_map = event_products['track_id_hit_segment_map']

        trajectories_traj_ids = trajectories['traj_id']
        trajectories_vertex_ids = trajectories['vertex_id']
        trajectories_pdg_ids = trajectories['pdg_id']

        """Grab the electrons which are from compton scatters, conversions or electron ionization"""
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

            """############################### Set standard labels ###############################"""
            """Set the event_id"""
            arrakis_charge['event_id'][particle_hits] = event

            """############################### Photons ###############################"""
            if len(particle_hits) >= self.shower_size_threshold:
                for label, value in self.photon_shower_labels.items():
                    arrakis_charge[label][particle_hit_segments] = value
            else:
                for label, value in self.photon_blip_labels.items():
                    arrakis_charge[label][particle_hit_segments] = value

            """Set the particle label"""
            arrakis_charge['particle'][particle_hit_segments] = trajectories_pdg_ids[particle_mask][ii]

            """Generate unique label"""
            arrakis_charge['unique_topology'][particle_hit_segments] = traj_index

        """Write changes to arrakis_file"""
        arrakis_file['charge_segment/calib_final_hits/data'][event_indices['charge']] = arrakis_charge
