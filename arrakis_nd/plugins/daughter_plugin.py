"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
from arrakis_nd.plugins.plugin import Plugin


class DaughterPlugin(Plugin):
    """
    A plugin for constructing daughter arrays for event trajectories.

    The FLOW files do not contain daughter information, most likely because
    it would need to exist either as a jagged array, which is not the best
    option for h5 files, or as a large array possibly 100's of columns in
    order to account for trajectories which have a large number of daughters.
    This is similar to what happens with the hit -> segment ID associations
    which use a 200 wide array in order to accomodate hits which have a large
    number of contributing segments in truth.  Most of the entries are zero
    however, which is also the case with daughters.  We therefore have to
    construct the list of daughters on the fly, which is time consuming, but
    otherwise necessary for other parts of the logic.

    In MiniRun4.5 there was a change to the traj_id label, which was converted
    back to the label from edep-sim, so that each traj_id no longer resets in a
    file, but resets per interaction.  Because of this, finding associations
    between traj_id, parent_id and now daughters must include a match to
    vertex_id.  The old way to do this would have been:

        traj_ids = trajectories['traj_id']
        parent_ids = trajectories['parent_id']

        for ii, traj_id in enumerate(traj_ids):
            indices = np.where(parent_ids == traj_id)[0]
            # Store these indices in the dictionary
            daughters[ii] = list(indices)

    but is now changed to reflect the new logic.  This is confirmed to work with the
    old <4.5 miniruns as well.

    The above logic is quite slow (repeated calls to np.where).  In order to
    optimize this, we opt for this better solution:

        unique_pairs = set(zip(parent_ids, vertex_ids))
        pair_to_indices = {pair: [] for pair in unique_pairs}

        for index, pair in enumerate(zip(parent_ids, vertex_ids)):
            pair_to_indices[pair].append(index)

        for ii, (traj_id, vertex_id) in enumerate(zip(traj_ids, vertex_ids)):
            # Direct lookup instead of np.where
            daughters[ii] = pair_to_indices.get((traj_id, vertex_id), [])

        event_products['daughters'] = daughters

    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        """
        super(DaughterPlugin, self).__init__(config)

        self.input_products = None
        self.output_products = 'daughters'

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
        daughters = [[] for ii in range(len(trajectories))]

        traj_ids = trajectories['traj_id']
        vertex_ids = trajectories['vertex_id']
        parent_ids = trajectories['parent_id']

        unique_pairs = set(zip(parent_ids, vertex_ids))
        pair_to_indices = {pair: [] for pair in unique_pairs}

        for index, pair in enumerate(zip(parent_ids, vertex_ids)):
            pair_to_indices[pair].append(index)

        for ii, (traj_id, vertex_id) in enumerate(zip(traj_ids, vertex_ids)):
            # Direct lookup instead of np.where
            daughters[ii] = pair_to_indices.get((traj_id, vertex_id), [])

        event_products['daughters'] = np.array(daughters, dtype=object)
