"""
"""
import numpy as np


class DetectorPointCloud:
    def __init__(
        self,
    ):
        self.clear()

    def clear(self):
        self.data = {
            'event': -1,
            'x': np.array([]),
            'y': np.array([]),
            'z': np.array([]),
            't_drift': np.array([]),
            'ts_pps': np.array([]),
            'Q': np.array([]),
            'E': np.array([]),
            'segment_ids': np.array([]),
            'segment_fractions': np.array([]),

            'topology_label':   np.array([]),
            'particle_label':   np.array([]),
            'physics_micro_label':    np.array([]),
            'physics_meso_label':     np.array([]),
            'physics_macro_label':    np.array([]),

            'unique_topology_label':   np.array([]),
            'unique_particle_label':   np.array([]),
            'unique_physics_micro_label':    np.array([]),
            'unique_physics_meso_label':     np.array([]),
            'unique_physics_macro_label':    np.array([]),

            'topology_labels':   np.array([]),
            'particle_labels':   np.array([]),
            'physics_micro_labels':    np.array([]),
            'physics_meso_labels':     np.array([]),
            'physics_macro_labels':    np.array([]),

            'unique_topology_labels':   np.array([]),
            'unique_particle_labels':   np.array([]),
            'unique_physics_micro_labels':    np.array([]),
            'unique_physics_meso_labels':     np.array([]),
            'unique_physics_macro_labels':    np.array([]),
        }

    def add_event(
        self,
        x:          float,
        y:          float,
        z:          float,
        t_drift:    float,
        ts_pps:     float,
        Q:          float,
        E:          float,
        segment_ids:    list = [],
        segment_fractions:  list = [],
    ):
        # Create a dictionary with the event data
        event_data = {
            'x': x,
            'y': y,
            'z': z,
            't_drift': t_drift,
            'ts_pps': ts_pps,
            'Q': Q,
            'E': E,
            'topology_label': np.full(x.shape, -1),
            'particle_label': np.full(x.shape, -1),
            'physics_micro_label': np.full(x.shape, -1),
            'physics_meso_label': np.full(x.shape, -1),
            'physics_macro_label': np.full(x.shape, -1),
            'unique_topology_label': np.full(x.shape, -1),
            'unique_particle_label': np.full(x.shape, -1),
            'unique_physics_micro_label': np.full(x.shape, -1),
            'unique_physics_meso_label': np.full(x.shape, -1),
            'unique_physics_macro_label': np.full(x.shape, -1),
            'segment_ids': segment_ids,
            'segment_fractions': segment_fractions,
            'topology_labels': np.full((x.shape[0], len(segment_ids)), -1),
            'particle_labels': np.full((x.shape[0], len(segment_ids)), -1),
            'physics_micro_labels': np.full((x.shape[0], len(segment_ids)), -1),
            'physics_meso_labels': np.full((x.shape[0], len(segment_ids)), -1),
            'physics_macro_labels': np.full((x.shape[0], len(segment_ids)), -1),
            'unique_topology_labels': np.full((x.shape[0], len(segment_ids)), -1),
            'unique_particle_labels': np.full((x.shape[0], len(segment_ids)), -1),
            'unique_physics_micro_labels': np.full((x.shape[0], len(segment_ids)), -1),
            'unique_physics_meso_labels': np.full((x.shape[0], len(segment_ids)), -1),
            'unique_physics_macro_labels': np.full((x.shape[0], len(segment_ids)), -1),
        }
        for key, item in event_data.items():
            self.data[key] = item

    def add_point(
        self,
        x:          float,
        y:          float,
        z:          float,
        t_drift:    float,
        ts_pps:     float,
        Q:          float,
        E:          float,
        segment_ids:    list = [],
        segment_fractions:  list = [],
    ):
        for key, value in {
            'x': x,
            'y': y,
            'z': z,
            't_drift': t_drift,
            'ts_pps': ts_pps,
            'Q': Q,
            'E': E,
            'topology_label': -1,
            'particle_label': -1,
            'physics_micro_label': -1,
            'physics_meso_label': -1,
            'physics_macro_label': -1,
            'unique_topology_label': -1,
            'unique_particle_label': -1,
            'unique_physics_micro_label': -1,
            'unique_physics_meso_label': -1,
            'unique_physics_macro_label': -1,
        }.items():
            self.data[key] = np.append(self.data[key], value)

        for key, value in {
            'segment_ids': segment_ids,
            'segment_fractions': segment_fractions,
        }.items():
            if self.data[key].size == 0:
                self.data[key] = np.array(value)
            else:
                self.data[key] = np.vstack((self.data[key], value))

        minus_one_array = np.full((1, len(segment_ids)), -1)        # max number of segments per hit is 200, reality max 10

        for key in [
            'topology_labels',
            'particle_labels',
            'physics_micro_labels',
            'physics_meso_labels',
            'physics_macro_labels',
            'unique_topology_labels',
            'unique_particle_labels',
            'unique_physics_micro_labels',
            'unique_physics_meso_labels',
            'unique_physics_macro_labels',
        ]:
            if self.data[key].size == 0:
                self.data[key] = minus_one_array
            else:
                self.data[key] = np.vstack((self.data[key], minus_one_array))

    def show_point(self):
        for key, value in self.data.items():
            print(f"{key}: {value.shape}")
