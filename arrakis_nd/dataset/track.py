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
        self.mc_tracks = {
            "event": -1,
            "track_id": [],
            "E_start": [],
            "xyz_start": [],
            "pxyz_start": [],
            "t_start": [],
            "hit_start": [],
            "E_end": [],
            "xyz_end": [],
            "pxyz_end": [],
            "t_end": [],
            "hit_end": [],
            "track_length": [],
            "dEdx": [],
            "Q_total": [],
            "hits": [],
        }

    def add_mc_track(
        self,
        track_id,
        E_start,
        xyz_start,
        pxyz_start,
        t_start,
        hit_start,
        E_end,
        xyz_end,
        pxyz_end,
        t_end,
        hit_end,
        track_length,
        dEdx,
        Q_total,
        hits
    ):
        mc_tracks = {
            "track_id": track_id,
            "E_start": E_start,
            "xyz_start": xyz_start,
            "pxyz_start": pxyz_start,
            "t_start": t_start,
            "hit_start": hit_start,
            "E_end": E_end,
            "xyz_end": xyz_end,
            "pxyz_end": pxyz_end,
            "t_end": t_end,
            "hit_end": hit_end,
            "track_length": track_length,
            "dEdx": dEdx,
            "Q_total": Q_total,
            "hits": hits,
        }
        for key, item in mc_tracks.items():
            self.mc_tracks[key].append(item)
