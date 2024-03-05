"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
from arrakis_nd.plugins.plugin import Plugin


class DaughterPlugin(Plugin):
    """
    A plugin for constructing daughter arrays for event trajectories.
    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        """
        super(DaughterPlugin, self).__init__(config)

        self.input_products = None
        self.output_product = 'daughters'

    @profiler
    def process_event(
        self,
        flow_file: h5py.File,
        arrakis_file: h5py.File,
        event_indices: dict,
        event_products: dict,
    ):
        """
        """
        trajectories = flow_file['mc_truth/trajectories/data'][event_indices['trajectories']]
        daughters = [[] for ii in range(len(trajectories))]

        for ii, traj_id in enumerate(trajectories['traj_id']):
            indices = np.where(trajectories['parent_id'] == traj_id)[0]
            # Store these indices in the dictionary
            daughters[ii] = list(indices)

        event_products['daughters'] = daughters
