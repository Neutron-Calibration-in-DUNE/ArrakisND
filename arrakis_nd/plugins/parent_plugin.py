"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
from arrakis_nd.plugins.plugin import Plugin


class ParentPlugin(Plugin):
    """
    A plugin for constructing an array of the parent particles
    pdg code corresponding to the index of the original particle
    in the trajectories array.
    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        """
        super(ParentPlugin, self).__init__(config)

        self.input_products = None
        self.output_products = 'parent_pdg_id'

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
        parent_pdg_id = [-1 for ii in range(len(trajectories))]

        traj_ids = trajectories['traj_id']
        vertex_ids = trajectories['vertex_id']
        parent_ids = trajectories['parent_id']
        pdg_ids = trajectories['pdg_id']

        for ii, (traj_id, vertex_id, parent_id) in enumerate(zip(
            traj_ids,
            vertex_ids,
            parent_ids
        )):
            if parent_id != -1:
                parent_index = np.where(
                    (traj_ids == parent_id) & (vertex_ids == vertex_id)
                )[0]
                parent_pdg_id[ii] = pdg_ids[parent_index][0]
        event_products['parent_pdg_id'] = np.array(parent_pdg_id)
