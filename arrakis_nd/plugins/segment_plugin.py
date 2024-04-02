"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
from arrakis_nd.plugins.plugin import Plugin


class SegmentPlugin(Plugin):
    """
    A plugin for ...
    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        This plugin generates ...

        """
        super(SegmentPlugin, self).__init__(config)

        self.input_products = [
            'parent_pdg_id'
        ]
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
        Here we associate (traj_id, vertex_id) pairs to indices in the event_indices['charge']
        mask so that later plugins can easily grab the hits associated to each particle.
        To do this we must first associate (traj_id, vertex_id) to (segment_id) and then
        use the back tracking information in 'mc_truth/calib_final_hit_backtrack/data'
        """
        trajectories = flow_file['mc_truth/trajectories/data'][event_indices['trajectories']]
        segments = flow_file['mc_truth/segments/data'][event_indices['segments']]
        charge = flow_file['charge/calib_final_hits/data'][event_indices['charge']]
        arrakis_charge = arrakis_file['charge_segment/calib_final_hits/data'][event_indices['charge']]
        charge_back_track = flow_file['mc_truth/calib_final_hit_backtrack/data'][event_indices['charge']]

        trajectories_traj_ids = trajectories['traj_id']
        trajectories_vertex_ids = trajectories['vertex_id']
        trajectories_pdg_ids = trajectories['pdg_id']
        trajectories_start_process = trajectories['start_process']
        trajectories_start_subprocess = trajectories['start_subprocess']
        segments_traj_ids = segments['traj_id']
        segments_vertex_ids = segments['vertex_id']
        segments_segment_ids = segments['segment_id']
        segments_x = segments['x']
        segments_y = segments['y']
        segments_z = segments['z']
        charge_x = charge['x']
        charge_y = charge['y']
        charge_z = charge['z']
        charge_segment_ids = charge_back_track['segment_id'].astype(int)
        charge_segment_fraction = charge_back_track['fraction']
        charge_segment_fraction_mask = (charge_segment_fraction == 0)
        charge_segment_ids[charge_segment_fraction_mask] = -1

        """Loop over segment ids"""
        for ii, segment_ids in enumerate(charge_segment_ids):
            """Get segment_ids of segments of each hit"""
            segment_ids = segment_ids[(segment_ids != -1)]
            """Determine where these ids are in the segments array"""
            segment_indices = np.array([
                np.where(segments_segment_ids == segment_id)[0][0]
                if segment_id in segments_segment_ids else -1
                for segment_id in segment_ids
            ])
            """Associate this hit to the (traj_id, vertex_id) pairs"""
            for jj, segment_index in enumerate(segment_indices):
                traj_id = segments_traj_ids[segment_index]
                vertex_id = segments_vertex_ids[segment_index]
                traj_index = np.where(
                    (trajectories_traj_ids == traj_id) &
                    (trajectories_vertex_ids == vertex_id)
                )
                arrakis_charge['particle'][ii][jj] = trajectories_pdg_ids[traj_index]
                arrakis_charge['segment_parent'][ii][jj] = event_products['parent_pdg_id'][traj_index]
                arrakis_charge['segment_start_process'][ii][jj] = trajectories_start_process[traj_index]
                arrakis_charge['segment_start_subprocess'][ii][jj] = trajectories_start_subprocess[traj_index]
                arrakis_charge['segment_distance'][ii][jj] = np.sqrt(
                    (charge_x[ii] - segments_x[segment_index]) ** 2 +
                    (charge_y[ii] - segments_y[segment_index]) ** 2 +
                    (charge_z[ii] - segments_z[segment_index]) ** 2
                )
        arrakis_file['charge_segment/calib_final_hits/data'][event_indices['charge']] = arrakis_charge
