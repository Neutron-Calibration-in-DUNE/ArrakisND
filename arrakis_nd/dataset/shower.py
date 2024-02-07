"""
"""
import numpy as np
import copy


class Shower:
    def __init__(
        self,
    ):
        self.clear()

    def clear(self):
        self.mc_shower = {
            "event": -1,
        }

    def add_mc_showers(
        self,
    ):
        mc_shower = {
        }
        for key, item in mc_shower.items():
            self.mc_shower[key] = item