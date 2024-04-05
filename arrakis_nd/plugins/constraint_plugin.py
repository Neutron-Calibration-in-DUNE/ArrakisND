"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
from arrakis_nd.plugins.plugin import Plugin
from arrakis_nd.dataset.common import (
    Topology, Physics
)


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
        self.output_products = []

        if "segment_influence_cut" not in self.config:
            self.config["segment_influence_cut"] = 1.0
        self.segment_influence_cut = self.config["segment_influence_cut"]

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
            if Physics.MichelElectron.value in segment_physics:
                segment_fraction[(segment_physics != Physics.MichelElectron.value)] = 0
            if Physics.DeltaElectron.value in segment_physics:
                segment_fraction[(segment_physics != Physics.DeltaElectron.value)] = 0
            if Topology.Track.value in segment_topology:
                segment_fraction[(segment_topology != Topology.Track.value)] = 0
            if Topology.Shower.value in segment_topology:
                segment_fraction[(segment_topology != Topology.Shower.value)] = 0

            max_fraction = np.argmax(segment_fraction)
            arrakis_charge['event_id'][ii] = segment_event
            arrakis_charge['topology'][ii] = segment_topology[max_fraction]
            arrakis_charge['particle'][ii] = segment_particle[max_fraction]
            arrakis_charge['physics'][ii] = segment_physics[max_fraction]
            arrakis_charge['unique_topology'][ii] = segment_unique_topology[max_fraction]
            arrakis_charge['vertex'][ii] = int(np.any(segment_vertex))
            arrakis_charge['tracklette_begin'][ii] = int(np.any(segment_tracklette_begin))
            arrakis_charge['tracklette_end'][ii] = int(np.any(segment_tracklette_end))
            arrakis_charge['fragment_begin'][ii] = int(np.any(segment_fragment_begin))
            arrakis_charge['fragment_end'][ii] = int(np.any(segment_fragment_end))
            arrakis_charge['shower_begin'][ii] = int(np.any(segment_shower_begin))

        """Write changes to arrakis_file"""
        arrakis_file['charge/calib_final_hits/data'][event_indices['charge']] = arrakis_charge
