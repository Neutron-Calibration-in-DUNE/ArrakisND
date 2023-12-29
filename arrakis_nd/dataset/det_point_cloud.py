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
            "event": -1,
            "x": np.array([]),
            "y": np.array([]),
            "z": np.array([]),
            "t_drift": np.array([]),
            "ts_pps": np.array([]),
            "Q": np.array([]),
            "E": np.array([]),
            "n_photons": np.array([]),
            "segment_ids": np.array([]),
            "segment_fractions": np.array([]),
            "topology_label": np.array([]),
            "particle_label": np.array([]),
            "physics_micro_label": np.array([]),
            "physics_meso_label": np.array([]),
            "physics_macro_label": np.array([]),
            "unique_topology_label": np.array([]),
            "unique_particle_label": np.array([]),
            "unique_physics_micro_label": np.array([]),
            "unique_physics_meso_label": np.array([]),
            "unique_physics_macro_label": np.array([]),
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
            "t_drift": t_drift,
            "ts_pps": ts_pps,
            "Q": Q,
            "E": E,
            "n_photons": n_photons,
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
            "segment_ids": segment_ids,
            "segment_fractions": segment_fractions,
            "topology_labels": empty_labels,
            "particle_labels": empty_labels,
            "physics_micro_labels": empty_labels,
            "physics_meso_labels": empty_labels,
            "physics_macro_labels": empty_labels,
            "unique_topology_labels": empty_labels,
            "unique_particle_labels": empty_labels,
            "unique_physics_micro_labels": empty_labels,
            "unique_physics_meso_labels": empty_labels,
            "unique_physics_macro_labels": empty_labels,
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


class OpticalPointCloud:
    def __init__(
        self,
    ):
        self.clear()

    def clear(self):
        self.data = {
            "event": -1,
            "tpc": np.array([]), # there are 8 tpc's
            "channel": np.array([]), # 64? channels per tpc
            "tick": np.array([]), # 1000 ticks per channel, 1 tick is 16 ns
            "max_peak": np.array([]), # the height of the peak in ADC counts
            "fwhm_peak": np.array([]), # the full width at half maximum of the peak in ticks
            "integral_peak": np.array([]), # the integral of the peak in ADC counts
            "n_photons": np.array([]), # the number of photons in the peak (truth)
            "segment_ids": np.array([]), # the segment ids that created the peak (truth)
        }
    
    def add_point(
        self,
        tpc: float,
        channel: float,
        tick: float,
        max_peak: float,
        fwhm_peak: float,
        integral_peak: float,
        n_photons: float,
        segment_ids: list = [],
    ):
        """
        Add a single hit to the point cloud. So, a single peak, caused by possibly 
        multiple segments, registered by a single channel.
        """
        # Create a dictionary with the event data
        event_data = {
            "tpc": tpc,
            "channel": channel,
            "tick": tick,
            "max_peak": max_peak,
            "fwhm_peak": fwhm_peak,
            "integral_peak": integral_peak,
            "n_photons": n_photons,
            "segment_ids": segment_ids,
        }
        for key, item in event_data.items():
            self.data[key] = item

    def add_event(
        self,
        tpc: list = [],
        channel: list = [],
        tick: list = [],
        max_peak: list = [],
        fwhm_peak: list = [],
        integral_peak: list = [],
        n_photons: list = [],
        segment_id: list = [],
    ):
        """
        Add all the peaks of a single event to the point cloud. So, all the peaks,
        registered on any channel.
        """
        for key, value in {
            "tpc": tpc,
            "channel": channel,
            "tick": tick,
            "max_peak": max_peak,
            "fwhm_peak": fwhm_peak,
            "integral_peak": integral_peak,
            "n_photons": n_photons,
            "segment_id": segment_id,
        }.items():
            self.data[key] = np.array(value)


    def show_point(self):
        for key, value in self.data.items():
            print(f"{key}: {value.shape}")
