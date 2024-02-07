"""
"""
import numpy as np
import copy


class Track:
    def __init__(
        self,
    ):
        self.clear()

    def clear(self):
        self.mc_track = {
            "event": -1,
        }

    def add_mc_tracks(
        self,
    ):
        mc_track = {
        }
        for key, item in mc_track.items():
            self.mc_track[key] = item