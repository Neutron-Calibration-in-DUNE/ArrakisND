"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
from arrakis_nd.plugins.plugin import Plugin


class AncestorPlugin(Plugin):
    """
    A plugin for constructing a map of the ancestor information
    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        """
        super(AncestorPlugin, self).__init__(config)

        self.input_products = None
        self.output_products = [
            'ancestor_traj_index_map',
            'ancestor_traj_id_map',
            'ancestor_pdg_id_map',
            'ancestor_level_map'
        ]

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

        traj_ids = trajectories['traj_id']
        vertex_ids = trajectories['vertex_id']
        parent_ids = trajectories['parent_id']
        pdg_ids = trajectories['pdg_id']

        """Create empty map with (traj_id, vertex_id) as keys"""
        ancestor_traj_index_map = {
            (traj_id, vertex_id): -1
            for (traj_id, vertex_id) in zip(traj_ids, vertex_ids)
        }
        ancestor_traj_id_map = {
            (traj_id, vertex_id): -1
            for (traj_id, vertex_id) in zip(traj_ids, vertex_ids)
        }
        ancestor_pdg_id_map = {
            (traj_id, vertex_id): -1
            for (traj_id, vertex_id) in zip(traj_ids, vertex_ids)
        }
        ancestor_level_map = {
            (traj_id, vertex_id): 0
            for (traj_id, vertex_id) in zip(traj_ids, vertex_ids)
        }

        for ii, (traj_id, vertex_id, parent_id) in enumerate(zip(
            traj_ids,
            vertex_ids,
            parent_ids
        )):
            """If parent_id == -1, then this is a primary"""
            if parent_id == -1:
                ancestor_traj_index_map[(traj_id, vertex_id)] = ii
                ancestor_traj_id_map[(traj_id, vertex_id)] = traj_id
                ancestor_pdg_id_map[(traj_id, vertex_id)] = pdg_ids[ii]
                ancestor_level_map[(traj_id, vertex_id)] = 0
                continue

            """Otherwise, continue to loop back through until we get to a primary"""
            temp_parent_id = parent_id
            parent_index = -1
            ancestor_level = 0
            while temp_parent_id != -1:
                ancestor_level += 1
                parent_index = np.where(
                    (traj_ids == temp_parent_id) & (vertex_ids == vertex_id)
                )[0][0]
                temp_parent_id = parent_ids[parent_index]
                if temp_parent_id != -1:
                    if ancestor_traj_id_map[(temp_parent_id, vertex_id)] != -1:
                        parent_index = ancestor_traj_index_map[(temp_parent_id, vertex_id)]
                        ancestor_level = ancestor_level_map[(temp_parent_id, vertex_id)] + 1
                        temp_parent_id = -1

            """Assign values from index to maps"""
            ancestor_traj_index_map[(traj_id, vertex_id)] = parent_index
            ancestor_traj_id_map[(traj_id, vertex_id)] = traj_ids[parent_index]
            ancestor_pdg_id_map[(traj_id, vertex_id)] = pdg_ids[parent_index]
            ancestor_level_map[(traj_id, vertex_id)] = ancestor_level

        event_products['ancestor_traj_index_map'] = ancestor_traj_index_map
        event_products['ancestor_traj_id_map'] = ancestor_traj_id_map
        event_products['ancestor_pdg_id_map'] = ancestor_pdg_id_map
        event_products['ancestor_level_map'] = ancestor_level_map
