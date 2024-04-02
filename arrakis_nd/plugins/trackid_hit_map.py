"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
from arrakis_nd.plugins.plugin import Plugin


class TrackIDHitMapPlugin(Plugin):
    """
    A plugin for ...
    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        This plugin generates a map from (traj_id, vertex_it) pairs to
        indices in the calib_final_hits array.  This allows other plugins
        to easily find which hits originate from each particle so that
        labels can be assigned in the corresponding ARRAKIS array.

        There are several approaches to associating (traj_id, vertex_id) to hits.  The
        slowest way is to iterate over the (traj_id, vertex_id) pairs in the
        trajectories array and try to see if any hits correspond to those ids.
        This can be done like this:

            for ii, (traj_id, vertex_id) in enumerate(
                zip(trajectories_traj_ids, trajectories_vertex_ids)
            ):
                segment_indices = np.where(
                    (segments_traj_ids == traj_id) & (segments_vertex_ids == vertex_id)
                )[0]
                segment_ids = segments_segment_ids[segment_indices]
                charge_ids = np.any(
                    np.isin(
                        charge_segment_ids, segment_ids
                    ),
                    axis=1,
                )

        This is unfortunately very slow, so instead we opt for the solution
        which is currently implemented in process_event.

        """
        super(TrackIDHitMapPlugin, self).__init__(config)

        self.input_products = None
        self.output_products = [
            'track_id_hit_map',
            'track_id_hit_segment_map',
            'track_id_hit_t0_map'
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
        Here we associate (traj_id, vertex_id) pairs to indices in the event_indices['charge']
        mask so that later plugins can easily grab the hits associated to each particle.
        To do this we must first associate (traj_id, vertex_id) to (segment_id) and then
        use the back tracking information in 'mc_truth/calib_final_hit_backtrack/data'
        """
        trajectories = flow_file['mc_truth/trajectories/data'][event_indices['trajectories']]
        segments = flow_file['mc_truth/segments/data'][event_indices['segments']]
        charge_back_track = flow_file['mc_truth/calib_final_hit_backtrack/data'][event_indices['charge']]

        trajectories_traj_ids = trajectories['traj_id']
        trajectories_vertex_ids = trajectories['vertex_id']
        segments_traj_ids = segments['traj_id']
        segments_vertex_ids = segments['vertex_id']
        segments_segment_ids = segments['segment_id']
        segments_t0s = segments['t0']
        charge_segment_ids = charge_back_track['segment_id'].astype(int)
        charge_segment_fraction = charge_back_track['fraction']
        charge_segment_fraction_mask = (charge_segment_fraction == 0)
        charge_segment_ids[charge_segment_fraction_mask] = -1

        """Create empty map with (traj_id, vertex_id) as keys"""
        track_id_hit_map = {
            (traj_id, vertex_id): []
            for (traj_id, vertex_id) in zip(trajectories_traj_ids, trajectories_vertex_ids)
        }
        track_id_hit_segment_map = {
            (traj_id, vertex_id): []
            for (traj_id, vertex_id) in zip(trajectories_traj_ids, trajectories_vertex_ids)
        }
        track_id_hit_t0_map = {
            (traj_id, vertex_id): []
            for (traj_id, vertex_id) in zip(trajectories_traj_ids, trajectories_vertex_ids)
        }

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
                track_id_hit_map[
                    (segments_traj_ids[segment_index], segments_vertex_ids[segment_index])
                ].append(ii)
                track_id_hit_segment_map[
                    (segments_traj_ids[segment_index], segments_vertex_ids[segment_index])
                ].append(jj)
                track_id_hit_t0_map[
                    (segments_traj_ids[segment_index], segments_vertex_ids[segment_index])
                ].append(segments_t0s[segment_index])

        event_products['track_id_hit_map'] = track_id_hit_map
        event_products['track_id_hit_segment_map'] = track_id_hit_segment_map
        event_products['track_id_hit_t0_map'] = track_id_hit_t0_map
