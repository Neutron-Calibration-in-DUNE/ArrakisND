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
            "vertex_id": np.array([]),
            "traj_id": np.array([]),
            "pdg_id": np.array([]),
            "E_start": np.array([]),
            "xyz_start": np.array([]),
            "pxyz_start": np.array([]),
            "t_start": np.array([]),
            "E_end": np.array([]),
            "xyz_end": np.array([]),
            "pxyz_end": np.array([]),
            "t_end": np.array([]),
            "start_process": np.array([]),
            "start_subprocess": np.array([]),
            "end_process": np.array([]),
            "end_subprocess": np.array([]),
            "parent_id": np.array([]),
            "daughters": np.array([]),
            "descendants": np.array([]),
            "ancestry": np.array([]),
        }

    def add_mc_particles(
        self,
        vertex_id,
        traj_id,
        pdg_id,
        E_start,
        xyz_start,
        pxyz_start,
        t_start,
        E_end,
        xyz_end,
        pxyz_end,
        t_end,
        start_process,
        start_subprocess,
        end_process,
        end_subprocess,
        parent_id,
        daughters,
        descendants,
        ancestry,
    ):
        mc_particles = {
            "vertex_id": vertex_id,
            "traj_id": traj_id,
            "pdg_id": pdg_id,
            "E_start": E_start,
            "xyz_start": xyz_start,
            "pxyz_start": pxyz_start,
            "t_start": t_start,
            "E_end": E_end,
            "xyz_end": xyz_end,
            "pxyz_end": pxyz_end,
            "t_end": t_end,
            "start_process": start_process,
            "start_subprocess": start_subprocess,
            "end_process": end_process,
            "end_subprocess": end_subprocess,
            "parent_id": parent_id,
            "daughters": daughters,
            "descendants": descendants,
            "ancestry": ancestry
        }
        for key, item in mc_particles.items():
            self.mc_particles[key] = item
