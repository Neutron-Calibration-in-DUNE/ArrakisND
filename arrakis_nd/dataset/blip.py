"""
"""
import numpy as np
import copy


class Blip:
    def __init__(
        self,
    ):
        self.clear()

    def clear(self):
        self.mc_blip = {
            "event": -1,
        }

    def add_mc_blips(
        self,
    ):
        mc_blip = {
        }
        for key, item in mc_blip.items():
            self.mc_blip[key] = item