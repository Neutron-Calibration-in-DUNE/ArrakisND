"""
"""
import numpy as np
import copy


class DetectorPointCloud:
    def __init__(
        self,
    ):
        self.clear()

    def clear(self):
        self.data = {
            "event": -1,
            # low level reconstructed variables
            "x": np.array([]),
            "y": np.array([]),
            "z": np.array([]),
            "Q": np.array([]),
            "P": np.array([]),
            "ts_pps": np.array([]),
            # low level mc variables
            "E": np.array([]),
            "t_drift": np.array([]),
            "n_photons": np.array([]),
            # corresponding mc segments
            "segment_ids": np.array([]),
            "segment_fractions": np.array([]),
            # high level reconstructed variables
            "topology_label": np.array([]),
            "particle_label": np.array([]),
            "physics_micro_label": np.array([]),
            "physics_meso_label": np.array([]),
            "physics_macro_label": np.array([]),
            # high level clustering variables
            "unique_topology_label": np.array([]),
            "unique_particle_label": np.array([]),
            "unique_physics_micro_label": np.array([]),
            "unique_physics_meso_label": np.array([]),
            "unique_physics_macro_label": np.array([]),
            # high level reconstructed variable lists
            "topology_labels": np.array([]),
            "particle_labels": np.array([]),
            "physics_micro_labels": np.array([]),
            "physics_meso_labels": np.array([]),
            "physics_macro_labels": np.array([]),
            "unique_topology_labels": np.array([]),
            "unique_particle_labels": np.array([]),
            "unique_physics_micro_labels": np.array([]),
            "unique_physics_meso_labels": np.array([]),
            "unique_physics_macro_labels": np.array([]),
            # high level track topology variables
            "vertex": np.array([]),
            "track_begin": np.array([]),
            "track_end": np.array([]),
            "shower_begin": np.array([]),
            "shower_end": np.array([]),
            # high level causality variables
            "causal_parent": np.array([]),
            # high level event scalar variables
        }

    def add_event(
        self,
        x: float,
        y: float,
        z: float,
        t_drift: float,
        ts_pps: float,
        Q: float,
        E: float,
        n_photons: float,
        segment_ids: list = [],
        segment_fractions: list = [],
        segment_track_ids:  list = [],
    ):
        # Create a dictionary with the event data
        # segment_ids = np.array(segment_ids)
        # segment_fractions = np.array(segment_fractions)
        mask = segment_ids != 0
        segment_ids = np.array(
            [segment_ids[ii][m].astype(int) for ii, m in enumerate(mask)], dtype=object
        )
        segment_fractions = np.array(
            [segment_fractions[ii][m] for ii, m in enumerate(mask)], dtype=object
        )
        empty_labels = np.array(
            [np.full_like(subarray, -1) for subarray in segment_ids], dtype=object
        )
        event_data = {
            "x": x,
            "y": y,
            "z": z,
            "Q": Q,
            "ts_pps": ts_pps,
            "E": E,
            "t_drift": t_drift,
            "n_photons": n_photons,
            "segment_ids": segment_ids,
            "segment_fractions": segment_fractions,
            "segment_track_ids": segment_track_ids,
            "topology_label": np.full(x.shape, -1),
            "particle_label": np.full(x.shape, -1),
            "physics_micro_label": np.full(x.shape, -1),
            "physics_meso_label": np.full(x.shape, -1),
            "physics_macro_label": np.full(x.shape, -1),
            "unique_topology_label": np.full(x.shape, -1),
            "unique_particle_label": np.full(x.shape, -1),
            "unique_physics_micro_label": np.full(x.shape, -1),
            "unique_physics_meso_label": np.full(x.shape, -1),
            "unique_physics_macro_label": np.full(x.shape, -1),
            "topology_labels": copy.deepcopy(empty_labels),
            "particle_labels": copy.deepcopy(empty_labels),
            "physics_micro_labels": copy.deepcopy(empty_labels),
            "physics_meso_labels": copy.deepcopy(empty_labels),
            "physics_macro_labels": copy.deepcopy(empty_labels),
            "unique_topology_labels": copy.deepcopy(empty_labels),
            "unique_particle_labels": copy.deepcopy(empty_labels),
            "unique_physics_micro_labels": copy.deepcopy(empty_labels),
            "unique_physics_meso_labels": copy.deepcopy(empty_labels),
            "unique_physics_macro_labels": copy.deepcopy(empty_labels),
            "vertex": np.full(x.shape, 0),
            "track_begin": np.full(x.shape, 0),
            "track_end": np.full(x.shape, 0),
            "delta_begin": np.full(x.shape, 0),
            "delta_end": np.full(x.shape, 0),
            "shower_begin": np.full(x.shape, 0),
            "shower_end": np.full(x.shape, 0),
            "causal_parent": np.full(x.shape, -1),
        }
        for key, item in event_data.items():
            self.data[key] = item

    def add_point(
        self,
        x: float,
        y: float,
        z: float,
        t_drift: float,
        ts_pps: float,
        Q: float,
        E: float,
        n_photons: list = [],
        segment_ids: list = [],
        segment_fractions: list = [],
    ):
        for key, value in {
            "x": x,
            "y": y,
            "z": z,
            "t_drift": t_drift,
            "ts_pps": ts_pps,
            "Q": Q,
            "E": E,
            "topology_label": -1,
            "particle_label": -1,
            "physics_micro_label": -1,
            "physics_meso_label": -1,
            "physics_macro_label": -1,
            "unique_topology_label": -1,
            "unique_particle_label": -1,
            "unique_physics_micro_label": -1,
            "unique_physics_meso_label": -1,
            "unique_physics_macro_label": -1,
        }.items():
            self.data[key] = np.append(self.data[key], value)

        for key, value in {
            "segment_ids": segment_ids,
            "segment_fractions": segment_fractions,
            "n_photons": n_photons,
        }.items():
            if self.data[key].size == 0:
                self.data[key] = np.array(value)
            else:
                self.data[key] = np.vstack((self.data[key], value))

        minus_one_array = np.full(
            (1, len(segment_ids)), -1
        )  # max number of segments per hit is 200, reality max 10

        for key in [
            "topology_labels",
            "particle_labels",
            "physics_micro_labels",
            "physics_meso_labels",
            "physics_macro_labels",
            "unique_topology_labels",
            "unique_particle_labels",
            "unique_physics_micro_labels",
            "unique_physics_meso_labels",
            "unique_physics_macro_labels",
        ]:
            if self.data[key].size == 0:
                self.data[key] = minus_one_array
            else:
                self.data[key] = np.vstack((self.data[key], minus_one_array))

    def show_point(self):
        for key, value in self.data.items():
            print(f"{key}: {value.shape}")
