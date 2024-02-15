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
            "true_segment_ids": np.array([]),    # the segment ids that created the peak (truth)
            "true_tick": np.array([]),           # the tick of the segment (truth)
            "true_channel": np.array([]),        # the channel of the segment (truth)
            "true_pe": np.array([]),             # the number of photoelectrons of the segment (truth)
        }

    def add_point(
        self,
        tpc: float,
        channel: float,
        tick: float,
        max_peak: float,
        fwhm_peak: float,
        integral_peak: float,
        true_segment_ids: list = [],
        true_tick: list = [],
        true_channel: list = [],
        true_pe: list = [],
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
            "true_segment_ids": true_segment_ids,
            "true_tick": true_tick,
            "true_channel": true_channel,
            "true_pe": true_pe,
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
        true_segment_id: list = [],
        true_tick: list = [],
        true_channel: list = [],
        true_pe: list = [],        
    ):
        """
        Add all the peaks of a single event to the point cloud. So, all the peaks,
        registered on any channel.
        """
        print(true_segment_id)
        print(tick)
        print(channel)
        for key, value in {
            "tpc": tpc,
            "channel": channel,
            "tick": tick,
            "max_peak": max_peak,
            "fwhm_peak": fwhm_peak,
            "integral_peak": integral_peak,
            "true_segment_id": true_segment_id,
            "true_tick": true_tick,
            "true_channel": true_channel,
            "true_pe": true_pe,
        }.items():
            self.data[key] = np.array(value)

    def show_point(self):
        for key, value in self.data.items():
            print(f"{key}: {value.shape}")
