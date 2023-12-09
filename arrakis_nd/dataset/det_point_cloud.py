"""
"""
import numpy as np

class DetectorPointCloud:
    def __init__(self):
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

            'source_label': np.array([]),
            'topology_label': np.array([]),
            'particle_label': np.array([]),
            'physics_label': np.array([]),

            'unique_topology': np.array([]),
            'unique_particle': np.array([]),
            'unique_physics': np.array([]),

            'source_labels': np.array([]),
            'topology_labels': np.array([]),
            'particle_labels': np.array([]),
            'physics_labels': np.array([]),

            'unique_topologies': np.array([]),
            'unique_particles': np.array([]),
            'unique_physicses': np.array([]),
        }

    def add_point(self,
        x:      float, 
        y:      float, 
        z:      float, 
        t_drift:    float, 
        ts_pps: float, 
        Q:      float, 
        E:      float, 
        segment_ids:    list=[],
        segment_fractions:  list=[],
    ):
        for key, value in {
            'x': x,
            'y': y,
            'z': z,
            't_drift': t_drift,
            'ts_pps': ts_pps,
            'Q': Q,
            'E': E,
            'source_label': -1,
            'topology_label': -1,
            'particle_label': -1,
            'physics_label': -1,
            'unique_topology': -1,
            'unique_particle': -1,
            'unique_physics': -1,
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

        minus_one_array = np.full((1, len(segment_ids)), -1) # max number of segments per hit is 200, reality max 10

        for key in [
            'source_labels',
            'topology_labels',
            'particle_labels',
            'physics_labels',
            'unique_topologies',
            'unique_particles',
            'unique_physicses',
        ]:
            if self.data[key].size == 0:
                self.data[key] = minus_one_array
            else:
                self.data[key] = np.vstack((self.data[key], minus_one_array))

    def show_point(self):
        for key, value in self.data.items():
            print(f"{key}: {value.shape}")
