"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
from arrakis_nd.plugins.plugin import Plugin


class MuonPlugin(Plugin):
    """
    A plugin for labeling muons in charge and light data.

    Muons are MIPs which deposit energy through MIP ionization (i.e. the
    standard dE/dx mechanism).  This plugin grabs all muons in an event
    and labels their energy deposits according to the following scheme:

        (1) hits coming directly from muon ionization are labeled as
                topology = track
                physics_micro = mip_ionization
                physics_macro = mip
        (2) delta electrons which are daughters of muons are labeled as
                topology = track
                physics_micro = electron_ionization
                physics_macro = delta_electron
        (3) michel electrons which are daughters of muons are labeled as
                topology = track
                physics_micro = electron_ionization
                physics_macro = michel_electron

    Apart from these labels, track beginning and ending points are also
    labelled.  Vertices are placed at points where Michel electrons are
    generated.

    StandardRecord objects are also generated.
    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        """
        super(MuonPlugin, self).__init__(config)

        self.input_products = ['daughters', 'track_id_hit_map', 'track_id_hit_segment_map']
        self.output_products = None

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
        trajctories_pdg_ids = trajectories['pdg_id']

        muon_mask = (abs(trajctories_pdg_ids) == 13)

        for (muon_id, vertex_id) in zip(trajectories_traj_ids[muon_mask], trajectories_vertex_ids[muon_mask]):
            muon_hits = track_id_hit_map[(muon_id, vertex_id)]
            muon_segments = track_id_hit_segment_map[(muon_id, vertex_id)]
            arrakis_charge['topology'][muon_hits, muon_segments] = 1
