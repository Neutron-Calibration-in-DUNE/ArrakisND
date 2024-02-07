"""
"""
import numpy as np
import copy


class Scalar:
    def __init__(
        self,
    ):
        self.clear()

    def clear(self):
        self.mc_scalar = {
            "event": -1,
        }

    def add_mc_scalars(
        self,
    ):
        mc_scalar = {
        }
        for key, item in mc_scalar.items():
            self.mc_scalar[key] = item