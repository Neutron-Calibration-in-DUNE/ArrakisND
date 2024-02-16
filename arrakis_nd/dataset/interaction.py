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
        self.mc_interactions = {
            "event": -1,
            "vertex_id": [],
            "vertex": [],
            "target": [],
            "reaction": [],
            "Enu": [],
            "nu_4mom": [],
            "nu_pdg": [],
        }

    def add_mc_interactions(
        self,
        vertex_id,
        vertex,
        target,
        reaction,
        Enu,
        nu_4mom,
        nu_pdg
    ):
        mc_interactions = {
            "vertex_id": vertex_id,
            "vertex": vertex,
            "target": target,
            "reaction": reaction,
            "Enu": Enu,
            "nu_4mom": nu_4mom,
            "nu_pdg": nu_pdg
        }
        for key, item in mc_interactions.items():
            self.mc_interactions[key] = item
