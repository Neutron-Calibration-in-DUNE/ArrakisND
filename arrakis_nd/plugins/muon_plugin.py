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

        self.input_products = ['daughters', 'track_id_hit_map']
        self.output_product = None

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
        charge = flow_file['charge/calib_final_hits/data']

        muons = abs(trajectories['pdg_id']) == 13
