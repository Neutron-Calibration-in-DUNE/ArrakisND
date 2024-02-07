"""
"""
import numpy as np
import copy


class Interaction:
    def __init__(
        self,
    ):
        self.clear()

    def clear(self):
        self.mc_interaction = {
            "event": -1,
        }

    def add_mc_interactions(
        self,
    ):
        mc_interaction = {
        }
        for key, item in mc_interaction.items():
            self.mc_interaction[key] = item
