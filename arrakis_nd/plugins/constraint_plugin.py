"""
"""
import h5py
import numpy as np
from scipy.interpolate import splprep, splev
import warnings
from scipy.integrate import quad, IntegrationWarning

from arrakis_nd.utils.utils import profiler, integrand
from arrakis_nd.plugins.plugin import Plugin
from arrakis_nd.dataset.common import (
    ProcessType, SubProcessType,
    Topology, Physics
)


class ConstraintPlugin(Plugin):
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
        trajectories = flow_file['mc_truth/trajectories/data'][event_indices['trajectories']]
        charge = flow_file['charge/calib_final_hits/data'][event_indices['charge']]
        arrakis_charge = arrakis_file['charge/calib_final_hits/data'][event_indices['charge']]
        arrakis_segment_charge = arrakis_file['charge_segment/calib_final_hits/data'][event_indices['charge']]

        segments_event = arrakis_segment_charge['event_id']
        segments_particle = arrakis_segment_charge['particle']
        segments_parent = arrakis_segment_charge['segment_parent']
        segments_start_process = arrakis_segment_charge['segment_start_process']
        segments_start_subprocess = arrakis_segment_charge['segment_start_subprocess']
        segments_distance = arrakis_segment_charge['segment_distance']
        segments_topology = arrakis_segment_charge['topology']
        segments_physics = arrakis_segment_charge['physics']
        segments_fraction = arrakis_segment_charge['segment_fraction']

        """Iterate over the individual (hits)"""
        for ii, particle in enumerate(segments_particle):
            segment_event = segments_event[ii]
            segment_particle = segments_particle[ii]
            segment_parent = segments_parent[ii]
            segment_start_process = segments_start_process[ii]
            segment_start_subprocess = segments_start_subprocess[ii]
            segment_distance = segments_distance[ii]
            segment_topology = segments_topology[ii]
            segment_physics = segments_physics[ii]
            segment_fraction = segments_fraction[ii]
            """Make a cut on segment_influence distance"""
            """
            We want to limit the influence of segments to labeling by how far away they are
            For example, if we have a delta electron coming from a muon, it is very likely
            that induction from the muon will dominate the fraction of energy contained
            in the delta hits.
            """
            """First check to see if any segments remain after the cut"""
            if 3 in segment_physics:
                delta_index = np.where(segment_physics == 3)
                print(f"delta distance: {segment_distance[delta_index]}")
            if sum(segment_distance <= self.segment_influence_cut):
                segment_fraction[(segment_distance > self.segment_influence_cut)] = 0

            max_fraction = np.argmax(segment_fraction)
            arrakis_charge['event_id'][ii] = segment_event
            arrakis_charge['topology'][ii] = segment_topology[max_fraction]
            arrakis_charge['physics'][ii] = segment_physics[max_fraction]

        """Write changes to arrakis_file"""
        arrakis_file['charge/calib_final_hits/data'][event_indices['charge']] = arrakis_charge
