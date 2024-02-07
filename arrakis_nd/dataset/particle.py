"""
"""
import numpy as np
import copy


class Particle:
    def __init__(
        self,
    ):
        self.clear()

    def clear(self):
        self.mc_particles = {
            "event": -1,
            "track_id": np.array([]),
            "pdg_code": np.array([]),
            "energy_start":   np.array([]),
            "momentum_start": np.array([]),
            "t_start":        np.array([]),
            "energy_end":     np.array([]),
            "momentum_end":   np.array([]),
            "t_end":          np.array([]),
            "process":  np.array([]),
            "end_process": np.array([]),
            "parent_id": np.array([]),
            "daughters": np.array([]),
            "descendants":  np.array([]),
            "ancestry": np.array([]),
        }

    def add_mc_particles(
        self,
        track_id,
        pdg_code,
        energy_start,
        momentum_start,
        t_start,
        energy_end,
        momentum_end,
        t_end,
        process,
        end_process,
        parent_id,
        daughters,
        descendants,
        ancestry,
    ):
        mc_particles = {
            "track_id": track_id,
            "pdg_code": pdg_code,
            "energy_start":   energy_start,
            "momentum_start": momentum_start,
            "t_start": t_start,
            "energy_end":   energy_end,
            "momentum_end": momentum_end,
            "t_end": t_end,
            "process":  process,
            "end_process":  end_process,
            "parent_id":    parent_id,
            "daughters":    daughters,
            "descendants":  descendants,
            "ancestry": ancestry
        }
        for key, item in mc_particles.items():
            self.mc_particles[key] = item