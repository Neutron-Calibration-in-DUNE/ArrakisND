"""
"""
import numpy as np


class LightPointCloud:
    def __init__(
        self,
    ):
        self.clear()

    def clear(self):
        self.data = {
            "event": -1,
            "tpc": np.array([]),        # there are 8 tpc's
            "channel": np.array([]),    # 64? channels per tpc
            "tick": np.array([]),       # 1000 ticks per channel, 1 tick is 16 ns
            "max_peak": np.array([]),   # the height of the peak in ADC counts
            "fwhm_peak": np.array([]),  # the full width at half maximum of the peak in ticks
            "integral_peak": np.array([]),  # the integral of the peak in ADC counts
            "segment_ids": np.array([]),    # the segment ids that created the peak (truth)
        }

    def add_point(
        self,
        tpc: float,
        channel: float,
        tick: float,
        max_peak: float,
        fwhm_peak: float,
        integral_peak: float,
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
            "segment_id": segment_id,
        }.items():
            self.data[key] = np.array(value)

    def show_point(self):
        for key, value in self.data.items():
            print(f"{key}: {value.shape}")
