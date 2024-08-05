"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
from arrakis_nd.plugins.plugin import Plugin
from arrakis_nd.dataset.common import (
    Topology,
    Physics
)
from arrakis_nd.utils.logger import ArrakisError
from arrakis_nd.arrakis.common import undefined_data_type


class HeatMapPlugin(Plugin):
    """
    This plugin generates a Gaussian heat map around points of
    interest (vertices, track begin/end, fragment begin/end, etc.)

    Args:
        Plugin (_type_): _description_
    """
    def __init__(
        self,
        config: dict = {},
        meta: dict = {}
    ):
        """
        """
        super(HeatMapPlugin, self).__init__(config, meta)

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

        if "vertex_variance" not in self.config:
            self.config["vertex_variance"] = 1.0
        self.vertex_variance = self.config["vertex_variance"]

        if "tracklette_begin_variance" not in self.config:
            self.config["tracklette_begin_variance"] = 1.0
        self.tracklette_begin_variance = self.config["tracklette_begin_variance"]

        if "tracklette_end_variance" not in self.config:
            self.config["tracklette_end_variance"] = 1.0
        self.tracklette_end_variance = self.config["tracklette_end_variance"]

        if "fragment_begin_variance" not in self.config:
            self.config["fragment_begin_variance"] = 1.0
        self.fragment_begin_variance = self.config["fragment_begin_variance"]

        if "fragment_end_variance" not in self.config:
            self.config["fragment_end_variance"] = 1.0
        self.fragment_end_variance = self.config["fragment_end_variance"]

        if "shower_begin_variance" not in self.config:
            self.config["shower_begin_variance"] = 1.0
        self.shower_begin_variance = self.config["shower_begin_variance"]

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
        charge = flow_file[f'charge/calib_{self.meta["hit_type"]}_hits/data'][event_indices['charge']]
        arrakis_charge = arrakis_file[f'charge/calib_{self.meta["hit_type"]}_hits/data'][event_indices['charge']]
        charge_back_track = flow_file[f'mc_truth/calib_{self.meta["hit_type"]}_hit_backtrack/data'][event_indices['charge']]

        charge_x = charge['x']
        charge_y = charge['y']
        charge_z = charge['z']

        charge_points = np.array([
            charge_x,
            charge_y,
            charge_z
        ]).T

        charge_segment_ids = charge_back_track['segment_ids'].astype(int)
        charge_segment_fraction = charge_back_track['fraction']
        charge_segment_fraction_mask = (charge_segment_fraction == 0)
        charge_segment_ids[charge_segment_fraction_mask] = -1

        """Generate heat map for vertices"""
        vertex_indices = np.where(
            (arrakis_charge['vertex'] == 1)
        )[0]
        for jj, index in enumerate(vertex_indices):
            distances = np.sum((charge_points - charge_points[index]) ** 2, axis=1)
            arrakis_charge['vertex_heat_map'] += np.exp(-distances / (2 * self.vertex_variance))
        if len(vertex_indices) > 0:
            if np.max(arrakis_charge['vertex_heat_map']) > 0:
                arrakis_charge['vertex_heat_map'] = (
                    arrakis_charge['vertex_heat_map'] / np.max(arrakis_charge['vertex_heat_map'])
                )

        """Generate heat map for vertices"""
        tracklette_begin_indices = np.where(
            (arrakis_charge['tracklette_begin'] == 1)
        )[0]
        for jj, index in enumerate(tracklette_begin_indices):
            distances = np.sum((charge_points - charge_points[index]) ** 2, axis=1)
            arrakis_charge['tracklette_begin_heat_map'] += np.exp(-distances / (2 * self.tracklette_begin_variance))
        if len(tracklette_begin_indices) > 0:
            if np.max(arrakis_charge['tracklette_begin_heat_map']) > 0:
                arrakis_charge['tracklette_begin_heat_map'] = (
                    arrakis_charge['tracklette_begin_heat_map'] / np.max(arrakis_charge['tracklette_begin_heat_map'])
                )

        """Generate heat map for vertices"""
        tracklette_end_indices = np.where(
            (arrakis_charge['tracklette_end'] == 1)
        )[0]
        for jj, index in enumerate(tracklette_end_indices):
            distances = np.sum((charge_points - charge_points[index]) ** 2, axis=1)
            arrakis_charge['tracklette_end_heat_map'] += np.exp(-distances / (2 * self.tracklette_end_variance))
        if len(tracklette_end_indices) > 0:
            if np.max(arrakis_charge['tracklette_end_heat_map']) > 0:
                arrakis_charge['tracklette_end_heat_map'] = (
                    arrakis_charge['tracklette_end_heat_map'] / np.max(arrakis_charge['tracklette_end_heat_map'])
                )

        """Generate heat map for vertices"""
        fragment_begin_indices = np.where(
            (arrakis_charge['fragment_begin'] == 1)
        )[0]
        for jj, index in enumerate(fragment_begin_indices):
            distances = np.sum((charge_points - charge_points[index]) ** 2, axis=1)
            arrakis_charge['fragment_begin_heat_map'] += np.exp(-distances / (2 * self.fragment_begin_variance))
        if len(fragment_begin_indices) > 0:
            if np.max(arrakis_charge['fragment_begin_heat_map']) > 0:
                arrakis_charge['fragment_begin_heat_map'] = (
                    arrakis_charge['fragment_begin_heat_map'] / np.max(arrakis_charge['fragment_begin_heat_map'])
                )

        """Generate heat map for vertices"""
        fragment_end_indices = np.where(
            (arrakis_charge['fragment_end'] == 1)
        )[0]
        for jj, index in enumerate(fragment_end_indices):
            distances = np.sum((charge_points - charge_points[index]) ** 2, axis=1)
            arrakis_charge['fragment_end_heat_map'] += np.exp(-distances / (2 * self.fragment_end_variance))
        if len(fragment_end_indices) > 0:
            if np.max(arrakis_charge['fragment_end_heat_map']) > 0:
                arrakis_charge['fragment_end_heat_map'] = (
                    arrakis_charge['fragment_end_heat_map'] / np.max(arrakis_charge['fragment_end_heat_map'])
                )

        """Generate heat map for vertices"""
        shower_begin_indices = np.where(
            (arrakis_charge['shower_begin'] == 1)
        )[0]
        for jj, index in enumerate(shower_begin_indices):
            distances = np.sum((charge_points - charge_points[index]) ** 2, axis=1)
            arrakis_charge['shower_begin_heat_map'] += np.exp(-distances / (2 * self.shower_begin_variance))
        if len(shower_begin_indices) > 0:
            if np.max(arrakis_charge['shower_begin_heat_map']) > 0:
                arrakis_charge['shower_begin_heat_map'] = (
                    arrakis_charge['shower_begin_heat_map'] / np.max(arrakis_charge['shower_begin_heat_map'])
                )

        """Write changes to arrakis_file"""
        arrakis_file[f'charge/calib_{self.meta["hit_type"]}_hits/data'][event_indices['charge']] = arrakis_charge
