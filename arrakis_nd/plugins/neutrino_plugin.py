"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
from arrakis_nd.utils.utils import fiducialized_vertex
from arrakis_nd.plugins.plugin import Plugin
from arrakis_nd.dataset.common import Neutrino


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
            neutrino_type = Neutrino.Undefined.value
            fiducialized = 0
        else:
            """We only want events with one neutrino hitting the argon"""
            if sum(argon_neutrino_targets) > 1:
                return
            neutrino = interactions[argon_neutrino_targets]
            fiducialized = fiducialized_vertex(neutrino['vertex'][0])

            if neutrino['nu_pdg'] == 12:
                if neutrino['isCC'] is False:
                    neutrino_type = Neutrino.NCElectronNeutrino.value
                else:
                    neutrino_type = Neutrino.CCElectronNeutrino.value
            elif neutrino['nu_pdg'] == -12:
                if neutrino['isCC'] is False:
                    neutrino_type = Neutrino.NCAntiElectronNeutrino.value
                else:
                    neutrino_type = Neutrino.CCAntiElectronNeutrino.value
            elif neutrino['nu_pdg'] == 14:
                if neutrino['isCC'] is False:
                    neutrino_type = Neutrino.NCMuonNeutrino.value
                else:
                    neutrino_type = Neutrino.CCMuonNeutrino.value
            elif neutrino['nu_pdg'] == -14:
                if neutrino['isCC'] is False:
                    neutrino_type = Neutrino.NCAntiMuonNeutrino.value
                else:
                    neutrino_type = Neutrino.CCAntiMuonNeutrino.value
            elif neutrino['nu_pdg'] == 16:
                if neutrino['isCC'] is False:
                    neutrino_type = Neutrino.NCTauonNeutrino.value
                else:
                    neutrino_type = Neutrino.CCTauonNeutrino.value
            else:
                if neutrino['isCC'] is False:
                    neutrino_type = Neutrino.NCAntiTauonNeutrino.value
                else:
                    neutrino_type = Neutrino.CCAntiTauonNeutrino.value

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
        event_products['neutrino_event'].append(neutrino_data)
