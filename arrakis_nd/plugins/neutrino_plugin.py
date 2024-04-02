"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
from arrakis_nd.plugins.plugin import Plugin


class NeutrinoPlugin(Plugin):
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
        super(NeutrinoPlugin, self).__init__(config)

        self.input_products = [
            'daughters',
            'track_id_hit_map',
            'track_id_hit_segment_map',
            'track_id_hit_t0_map',
            'parent_pdg_id',
        ]
        self.output_products = []

        if "segment_influence_cut" not in self.config:
            self.config["segment_influence_cut"] = 1.0
        self.segment_influence_cut = self.config["segment_influence_cut"]

    @profiler
    def fiducialized_vertex(self, vert_pos):
        lar_x = [
            (-63.9273, -3.0652),
            (3.0652, 63.9273)
        ]
        lar_y = [(-62.055, 62.055)]
        lar_z = [
            (-64.51125, -2.48125),
            (2.48125, 64.51125)
        ]
        for i in lar_z:
            if vert_pos[2] < i[0] and vert_pos[2] > i[1]:
                return False
        for i in lar_x:
            if vert_pos[0] < i[0] and vert_pos[0] > i[1]:
                return False
        for i in lar_y:
            if vert_pos[1] < i[0] and vert_pos[1] > i[1]:
                return False
        return True

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
        """Determine neutrino type"""
        argon_neutrino_targets = (interactions['target'] == 18)
        if not sum(argon_neutrino_targets):
            neutrino_type = 0
            fiducialized = 0
        else:
            """We only want events with one neutrino hitting the argon"""
            if sum(argon_neutrino_targets) > 1:
                return
            neutrino = interactions[argon_neutrino_targets]
            fiducialized = self.fiducialized_vertex(neutrino['vertex'][0])
            if neutrino['isCC'] is False:
                neutrino_type = 1
            else:
                if abs(neutrino['nu_pdg']) == 12:
                    neutrino_type = 2
                elif abs(neutrino['nu_pdg']) == 14:
                    neutrino_type = 3
                else:
                    neutrino_type = 4

        neutrino_event_data_type = np.dtype([
            ('event_id', 'i4'),
            ('neutrino_type', 'i4'),
            ('fiducialized', 'i4')
        ])

        neutrino_data = np.array([(
            event,
            neutrino_type,
            fiducialized)],
            dtype=neutrino_event_data_type
        )

        """Add the new track to the CAF objects"""
        event_products['neutrino'].append(neutrino_data)
