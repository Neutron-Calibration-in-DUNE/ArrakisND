"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
from arrakis_nd.plugins.plugin import Plugin
from arrakis_nd.dataset.common import (
    Topology, Physics
)
from arrakis_nd.utils.logger import ArrakisError
from arrakis_nd.arrakis.common import undefined_data_type

class ConstraintPlugin(Plugin):
    """
    This plugin takes all the semantic information determine from the other
    plugins and decides how to "officially" label each calib_final_hit.  This is
    done in two steps,

        (1) apply an influence cut based on the distance of each contributing
            segment to the hit.
        (2) make heirarchical choices to the labels based on the following criteria:
            (a) michels take 1st precendence
            (b) deltas take 2nd precendence
            (c) tracks take 3rd precendence
            (d) showers take 4th precendence
            (e) blips take last precendence

    Args:
        Plugin (_type_): _description_
    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        """
        super(ConstraintPlugin, self).__init__(config)

        self.input_products = [
            'daughters',
            'track_id_hit_map',
            'track_id_hit_segment_map',
            'track_id_hit_t0_map',
            'parent_pdg_id',
        ]
        self.output_products = [
            'undefined'
        ]

        if "segment_influence_cut" not in self.config:
            self.config["segment_influence_cut"] = 1.0
        self.segment_influence_cut = self.config["segment_influence_cut"]

        if "constraint_mode" not in self.config:
            self.config["constraint_mode"] = "min_distance"
        if self.config["constraint_mode"] not in ["max_fraction", "min_distance"]:
            raise ArrakisError(f"specified 'constraint_mode' {self.config['constraint_mode']} not allowed!")
        self.constraint_mode = self.config["constraint_mode"]
        
        if "replace_undefined" not in self.config:
            self.config["replace_undefined"] = -1
        self.replace_undefined = self.config["replace_undefined"]

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
        arrakis_charge = arrakis_file['charge/calib_final_hits/data'][event_indices['charge']]
        arrakis_segment_charge = arrakis_file['charge_segment/calib_final_hits/data'][event_indices['charge']]
        charge_back_track = flow_file['mc_truth/calib_final_hit_backtrack/data'][event_indices['charge']]
        trajectories = flow_file['mc_truth/trajectories/data'][event_indices['trajectories']]
        segments = flow_file['mc_truth/segments/data'][event_indices['segments']]

        trajectories_traj_ids = trajectories['traj_id']
        trajectories_vertex_ids = trajectories['vertex_id']
        trajectories_parent_ids = trajectories['parent_id']
        trajectories_pdg_ids = trajectories['pdg_id']
        trajectories_start_process = trajectories['start_process']
        trajectories_start_subprocess = trajectories['start_subprocess']
        trajectories_end_process = trajectories['end_process']
        trajectories_end_subprocess = trajectories['end_subprocess']
        trajectories_parent_pdg_ids = event_products['parent_pdg_id']
        trajectories_ancestor_traj_ids = event_products['ancestor_traj_id_map']
        trajectories_ancestor_pdg_ids = event_products['ancestor_pdg_id_map']
        trajectories_ancestor_levels = event_products['ancestor_level_map']
        
        segments_traj_ids = segments['traj_id']
        segments_vertex_ids = segments['vertex_id']
        segments_segment_ids = segments['segment_id']

        segments_event = arrakis_segment_charge['event_id']
        segments_particle = arrakis_segment_charge['particle']
        segments_distance = arrakis_segment_charge['segment_distance']
        segments_topology = arrakis_segment_charge['topology']
        segments_physics = arrakis_segment_charge['physics']
        segments_unique_topology = arrakis_segment_charge['unique_topology']
        segments_fraction = arrakis_segment_charge['segment_fraction']
        segments_vertex = arrakis_segment_charge['vertex']
        segments_tracklette_begin = arrakis_segment_charge['tracklette_begin']
        segments_tracklette_end = arrakis_segment_charge['tracklette_end']
        segments_fragment_begin = arrakis_segment_charge['fragment_begin']
        segments_fragment_end = arrakis_segment_charge['fragment_end']
        segments_shower_begin = arrakis_segment_charge['shower_begin']

        charge_segment_ids = charge_back_track['segment_id'].astype(int)
        charge_segment_fraction = charge_back_track['fraction']
        charge_segment_fraction_mask = (charge_segment_fraction == 0)
        charge_segment_ids[charge_segment_fraction_mask] = -1

        """Iterate over the individual (hits)"""
        for ii, particle in enumerate(segments_particle):
            segment_event = segments_event[ii]
            segment_particle = segments_particle[ii]
            segment_distance = segments_distance[ii]
            segment_topology = segments_topology[ii]
            segment_physics = segments_physics[ii]
            segment_unique_topology = segments_unique_topology[ii]
            segment_fraction = segments_fraction[ii]
            segment_vertex = segments_vertex[ii]
            segment_tracklette_begin = segments_tracklette_begin[ii]
            segment_tracklette_end = segments_tracklette_end[ii]
            segment_fragment_begin = segments_fragment_begin[ii]
            segment_fragment_end = segments_fragment_end[ii]
            segment_shower_begin = segments_shower_begin[ii]
            """Make a cut on segment_influence distance"""
            """
            We want to limit the influence of segments to labeling by how far away they are
            For example, if we have a delta electron coming from a muon, it is very likely
            that induction from the muon will dominate the fraction of energy contained
            in the delta hits.
            """
            """First check to see if any segments remain after the cut"""
            if sum(segment_distance <= self.segment_influence_cut):
                segment_fraction[(segment_distance > self.segment_influence_cut)] = 0
            
            """Then, give priority to tracks"""
            if Topology.Track.value in segment_topology:
                segment_fraction[(segment_topology != Topology.Track.value)] = 0
            """Then, to showers"""
            if Topology.Shower.value in segment_topology:
                segment_fraction[(segment_topology != Topology.Shower.value)] = 0
            segment_distance[(segment_fraction == 0.0)] = 10e10

            if self.constraint_mode == "max_fraction":
                constraint_mask = np.argmax(segment_fraction)
            else:
                constraint_mask = np.argmin(segment_distance)

            arrakis_charge['event_id'][ii] = segment_event
            arrakis_charge['topology'][ii] = segment_topology[constraint_mask]
            arrakis_charge['particle'][ii] = segment_particle[constraint_mask]
            arrakis_charge['physics'][ii] = segment_physics[constraint_mask]
            arrakis_charge['unique_topology'][ii] = segment_unique_topology[constraint_mask]
            arrakis_charge['vertex'][ii] = int(np.any(segment_vertex))
            arrakis_charge['tracklette_begin'][ii] = int(np.any(segment_tracklette_begin))
            arrakis_charge['tracklette_end'][ii] = int(np.any(segment_tracklette_end))
            arrakis_charge['fragment_begin'][ii] = int(np.any(segment_fragment_begin))
            arrakis_charge['fragment_end'][ii] = int(np.any(segment_fragment_end))
            arrakis_charge['shower_begin'][ii] = int(np.any(segment_shower_begin))
            
            """Check for undefined values"""
            if arrakis_charge['topology'][ii] == -1:
                """Replace undefined values"""
                arrakis_charge['topology'][ii] = self.replace_undefined
                arrakis_charge['physics'][ii] = self.replace_undefined
                
                """Log the undefined information to the arrakis file"""
                segment_index = np.where(
                    segments_segment_ids == charge_segment_ids[ii][constraint_mask]
                )
                traj_id = segments_traj_ids[segment_index]
                vertex_id = segments_vertex_ids[segment_index]
                traj_index = np.where(
                    (trajectories_traj_ids == traj_id) &
                    (trajectories_vertex_ids == vertex_id)
                )
                parent_id = trajectories_parent_ids[traj_index]
                pdg_id = trajectories_pdg_ids[traj_index]
                parent_pdg_id = trajectories_parent_pdg_ids[traj_index]
                ancestor_id = trajectories_ancestor_traj_ids[(traj_id[0], vertex_id[0])]
                ancestor_pdg_id = trajectories_ancestor_pdg_ids[(traj_id[0], vertex_id[0])]
                ancestor_level = trajectories_ancestor_levels[(traj_id[0], vertex_id[0])]
                start_process = trajectories_start_process[traj_index]
                start_subprocess = trajectories_start_subprocess[traj_index]
                end_process = trajectories_end_process[traj_index]
                end_subprocess = trajectories_end_subprocess[traj_index]
                parent_index = np.where(
                    (trajectories_traj_ids == parent_id) &
                    (trajectories_vertex_ids == vertex_id)
                )
                if parent_id != -1:
                    parent_start_process = trajectories_start_process[parent_index]
                    parent_start_subprocess = trajectories_start_subprocess[parent_index]
                    parent_end_process = trajectories_end_process[parent_index]
                    parent_end_subprocess = trajectories_end_subprocess[parent_index]
                else:
                    parent_start_process = -1
                    parent_start_subprocess = -1
                    parent_end_process = -1
                    parent_end_subprocess = -1

                undefined_data = np.array([(
                    event,
                    vertex_id,
                    traj_id,
                    parent_id,
                    pdg_id,
                    parent_pdg_id,
                    ancestor_id,
                    ancestor_pdg_id,
                    ancestor_level,
                    start_process,
                    start_subprocess,
                    end_process,
                    end_subprocess,
                    parent_start_process,
                    parent_start_subprocess,
                    parent_end_process,
                    parent_end_subprocess)],
                    dtype=undefined_data_type
                )
                
                """Add the undefined data to the event products"""
                event_products['undefined'].append(undefined_data)
                
        """Write changes to arrakis_file"""
        arrakis_file['charge/calib_final_hits/data'][event_indices['charge']] = arrakis_charge
