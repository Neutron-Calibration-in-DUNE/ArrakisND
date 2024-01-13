"""

"""
import getpass
import socket
from datetime import datetime
import numpy as np
import copy

from arrakis_nd.utils.logger import Logger
from arrakis_nd.dataset.det_point_cloud import DetectorPointCloud
from arrakis_nd.dataset.light_point_cloud import LightPointCloud
from arrakis_nd.dataset.common import (
    TopologyLabel,
    PhysicsMicroLabel,
    PhysicsMesoLabel,
    PhysicsMacroLabel,
)
from arrakis_nd.dataset.common import classification_labels
from arrakis_nd.wrangler.common import wrangler_modes


class SimulationWrangler:
    """ """

    def __init__(self, name: str = "", config: dict = {}, meta: dict = {}):
        self.name = name + "_simulation_wrangler"
        self.config = config
        self.meta = meta

        if "device" in self.meta:
            self.device = self.meta["device"]
        else:
            self.device = "cpu"
        if meta["verbose"]:
            self.logger = Logger(self.name, output="both", file_mode="w")
        else:
            self.logger = Logger(self.name, level="warning", file_mode="w")

        self.topology_label = 0

        self.parse_config()

    def parse_config(self):
        self.check_config()
        self.parse_point_cloud()

    def check_config(self):
        if "wrangler_mode" not in self.config:
            self.logger.warn('wranger_mode not specified in config! setting to "map"')
            self.config["wrangler_mode"] = "map"
        if self.config["wrangler_mode"] not in wrangler_modes:
            self.logger.error(
                f'specified wrangler_mode {self.config["wrangler_mode"]} not allowed!'
            )
        self.wrangler_mode = self.config["wrangler_mode"]

        if self.wrangler_mode == "map":
            self.clear_event = self.clear_event_maps
            self.get_total_hit_energy = self.get_total_hit_energy_map
        elif self.wrangler_mode == "numpy":
            self.clear_event = self.clear_event_numpy
            self.get_total_hit_energy = self.get_total_hit_energy_numpy

        if "debug" not in self.config.keys():
            self.logger.warn('debug not specified in config! setting to "False"')
            self.config["debug"] = False
        self.debug = self.config["debug"]

    def parse_point_cloud(self):
        self.det_point_cloud = DetectorPointCloud()
        self.det_point_clouds = {}

        self.light_point_cloud = LightPointCloud()
        self.light_point_clouds = {}

        self.vertexid_vertex = {}
        self.vertexid_target = {}
        self.vertexid_reaction = {}
        self.vertexid_pdgcode = {}
        self.vertexid_label = {}

        self.trackid_vertexid = {}
        self.trackid_parentid = {}
        self.trackid_pdgcode = {}
        self.trackid_process = {}
        self.trackid_subprocess = {}
        self.trackid_endprocess = {}
        self.trackid_endsubprocess = {}
        self.trackid_energy = {}
        self.trackid_tstart = {}
        self.trackid_daughters = {}
        self.trackid_progeny = {}
        self.trackid_descendants = {}
        self.trackid_ancestorlevel = {}
        self.trackid_ancestry = {}

        self.trackid_segmentid = {}
        self.segmentid_trackid = {}

        self.trackid_hit = {}
        self.segmentid_hit = {}

        self.track_id_array = None

    def clear_event_maps(self):
        self.det_point_cloud.clear()
        self.light_point_cloud.clear()

        self.vertexid_vertex = {}
        self.vertexid_target = {}
        self.vertexid_reaction = {}
        self.vertexid_pdgcode = {}
        self.vertexid_label = {}

        self.trackid_vertexid = {}
        self.trackid_parentid = {}
        self.trackid_pdgcode = {}
        self.trackid_process = {}
        self.trackid_subprocess = {}
        self.trackid_endprocess = {}
        self.trackid_endsubprocess = {}
        self.trackid_energy = {}
        self.trackid_tstart = {}
        self.trackid_daughters = {}
        self.trackid_progeny = {}
        self.trackid_descendants = {}
        self.trackid_ancestorlevel = {}
        self.trackid_ancestry = {}

        self.trackid_segmentid = {}
        self.segmentid_trackid = {}

        self.trackid_hit = {}
        self.segmentid_hit = {}

    def clear_event_numpy(self):
        self.track_id_array = None

    def clear_point_clouds(self):
        self.det_point_clouds = {}
        self.light_point_clouds = {}

    def print_particle_data(self, particle):
        self.logger.info("## MCParticle #######################################")
        self.logger.info(f"## EventID:            [{'.' * 25}{self.det_point_cloud.data['event']}] ##")
        self.logger.info(f"## TrackID:            [{'.' * 25}{particle}] ##")
        self.logger.info(
            f"## PDG:                [{'.' * 25}{self.trackid_pdgcode[particle]}] ##"
        )
        self.logger.info(
            f"## Energy [MeV]:       [{'.' * 25}{self.trackid_energy[particle]}] ##"
        )
        self.logger.info(
            f"## Process:            [{'.' * 25}{self.trackid_process[particle]}] ##"
        )
        self.logger.info(
            f"## SubProcess:         [{'.' * 25}{self.trackid_subprocess[particle]}] ##"
        )
        self.logger.info(
            f"## EndProcess:         [{'.' * 25}{self.trackid_endprocess[particle]}] ##"
        )
        self.logger.info(
            f"## EndSubProcess:      [{'.' * 25}{self.trackid_endsubprocess[particle]}] ##"
        )
        self.logger.info(
            f"## Parent TrackID:     [{'.' * 25}{self.trackid_parentid[particle]}] ##"
        )
        if self.trackid_parentid[particle] != -1:
            self.logger.info(
                f"## Parent PDG:         [{'.' * 25}{self.trackid_pdgcode[self.trackid_parentid[particle]]}] ##"
            )
            self.logger.info(
                f"## Ancestor TrackID:   [{'.' * 25}{self.trackid_ancestry[particle][-1]}] ##"
            )
            self.logger.info(
                f"## Ancestor PDG:       [{'.' * 25}{self.trackid_pdgcode[self.trackid_ancestry[particle][-1]]}] ##"
            )
            self.logger.info(
                f"## Ancestor level:     [{'.' * 25}{self.trackid_ancestorlevel[particle]}] ##"
            )
        else:
            self.logger.info(f"## Parent PDG:         [{'.' * 25}{-1}] ##")
            self.logger.info(f"## Ancestor TrackID:   [{'.' * 25}{-1}] ##")
            self.logger.info(f"## Ancestor PDG:       [{'.' * 25}{-1}] ##")
            self.logger.info(f"## Ancestor level:     [{'.' * 25}{-1}] ##")
        self.logger.info("## Progeny  [.....level] [...TrackID] [.......PDG] ##")
        progeny = self.trackid_descendants[particle]
        particle_level = self.trackid_ancestorlevel[particle]
        for progeny_track_id in progeny:
            self.logger.info(
                f"##          [{'.' * 10}{self.trackid_ancestorlevel[progeny_track_id] - particle_level}] "
                + f"[{'.' * 10}{progeny_track_id}] [{'.' * 10}{self.trackid_pdgcode[progeny_track_id]}] ##"
            )
        self.logger.info("#####################################################")

    def copy_hit_labels(
        self,
        hit1,
        hit2,
    ):
        for label in [
            "topology_label",
            "physics_micro_label",
            "physics_meso_label",
            "physics_macro_label",
            "unique_topology_label",
            "unique_physics_micro_label",
            "unique_physics_meso_label",
            "unique_physics_macro_label"
        ]:
            self.det_point_cloud.data[label][hit1] = self.det_point_cloud.data[label][hit2]

    def set_hit_labels(
        self,
        hits: list = [],
        segments: list = [],
        trackid: int = 0,
        topology: TopologyLabel = TopologyLabel.Undefined,
        physics_micro: PhysicsMicroLabel = PhysicsMicroLabel.Undefined,
        physics_meso: PhysicsMesoLabel = PhysicsMesoLabel.Undefined,
        unique_topology: int = 0,
        unique_physics_micro: int = 0,
        unique_physics_meso: int = 0,
    ):
        for hit in hits:
            point_cloud_index = list(self.get_index_trackid(hit, trackid))
            vertexid = self.trackid_vertexid[trackid]
            physics_macro = self.vertexid_label[vertexid]
            if len(point_cloud_index) != 0:
                self.det_point_cloud.data["topology_labels"][hit][
                    point_cloud_index
                ] = topology.value
                self.det_point_cloud.data["particle_labels"][hit][
                    point_cloud_index
                ] = self.trackid_pdgcode[trackid]
                self.det_point_cloud.data["physics_micro_labels"][hit][
                    point_cloud_index
                ] = physics_micro.value
                self.det_point_cloud.data["physics_meso_labels"][hit][
                    point_cloud_index
                ] = physics_meso.value
                self.det_point_cloud.data["physics_macro_labels"][hit][
                    point_cloud_index
                ] = physics_macro.value
                self.det_point_cloud.data["unique_topology_labels"][hit][
                    point_cloud_index
                ] = unique_topology
                self.det_point_cloud.data["unique_particle_labels"][hit][
                    point_cloud_index
                ] = trackid
                self.det_point_cloud.data["unique_physics_micro_labels"][hit][
                    point_cloud_index
                ] = unique_physics_micro
                self.det_point_cloud.data["unique_physics_meso_labels"][hit][
                    point_cloud_index
                ] = unique_physics_meso
                self.det_point_cloud.data["unique_physics_macro_labels"][hit][
                    point_cloud_index
                ] = vertexid

                self.det_point_cloud.data["topology_label"][hit] = topology.value
                self.det_point_cloud.data["particle_label"][hit] = self.trackid_pdgcode[
                    trackid
                ]
                self.det_point_cloud.data["physics_micro_label"][hit] = physics_micro.value
                self.det_point_cloud.data["physics_meso_label"][hit] = physics_meso.value
                self.det_point_cloud.data["physics_macro_label"][hit] = physics_macro.value
                self.det_point_cloud.data["unique_topology_label"][hit] = unique_topology
                self.det_point_cloud.data["unique_particle_label"][hit] = trackid
                self.det_point_cloud.data["unique_physics_micro_label"][
                    hit
                ] = unique_physics_micro
                self.det_point_cloud.data["unique_physics_meso_label"][
                    hit
                ] = unique_physics_meso
                self.det_point_cloud.data["unique_physics_macro_label"][
                    hit
                ] = vertexid

    def set_hit_labels_list(
        self,
        hits: list = [[]],
        segments: list = [[]],
        trackid: int = 0,
        topology: TopologyLabel = TopologyLabel.Undefined,
        physics_micro: PhysicsMicroLabel = PhysicsMicroLabel.Undefined,
        physics_meso: PhysicsMesoLabel = PhysicsMesoLabel.Undefined,
        unique_topology: int = 0,
        unique_physics_micro: int = 0,
        unique_physics_meso: int = 0,
    ):
        for ii in range(len(hits)):
            self.set_hit_labels(
                hits[ii],
                segments[ii],
                trackid[ii],
                topology,
                physics_micro,
                physics_meso,
                unique_topology,
                unique_physics_micro,
                unique_physics_meso,
            )

    def set_labels_array(
        self,
        hits: list = [[[]]],
        segments: list = [[[]]],
        trackid: int = 0,
        topology: TopologyLabel = TopologyLabel.Undefined,
        physics_micro: PhysicsMicroLabel = PhysicsMicroLabel.Undefined,
        physics_meso: PhysicsMesoLabel = PhysicsMesoLabel.Undefined,
        unique_topology: int = 0,
        unique_physics_micro: int = 0,
        unique_physics_meso: int = 0,
    ):
        for ii in range(len(hits)):
            self.set_hit_labels_list(
                hits[ii],
                segments[ii],
                trackid[ii],
                topology,
                physics_micro,
                physics_meso,
                unique_topology,
                unique_physics_micro,
                unique_physics_meso,
            )

    def process_event(
        self,
        event_id,
        event_interactions,
        event_trajectories,
        event_segments,
        event_stacks,
        event_hits_back_track,
        event_hits,
        event_light,
    ):
        """
        Go through each of the individual logic functions
        to process mid-level and high-level variables.

        The list here should be exhaustive.  In debug mode, a
        set of check functions are run which gathers statistics
        on any failures, as well as timing and memory usage
        information for each function.
        """
        if self.debug:
            self._process_event_with_timing(
                event_id,
                event_interactions,
                event_trajectories,
                event_segments,
                event_stacks,
                event_hits_back_track,
                event_hits,
                event_light,
            )
        else:
            self._process_event_without_timing(
                event_id,
                event_interactions,
                event_trajectories,
                event_segments,
                event_stacks,
                event_hits_back_track,
                event_hits,
                event_light,
            )

    def _process_event_without_timing(
        self,
        event_id,
        event_interactions,
        event_trajectories,
        event_segments,
        event_stacks,
        event_hits_back_track,
        event_hits,
        event_light,
    ):
        self.clear_event()
        self.det_point_cloud.data["event"] = event_id
        self.light_point_cloud.data["event"] = event_id
        self.process_event_interactions(event_interactions)
        self.process_event_trajectories(event_trajectories)
        self.process_event_stacks(event_stacks)
        self.process_event_segments(event_segments)
        self.process_event_hits(event_hits, event_hits_back_track, event_segments)
        self.process_event_light(event_light)

    def _process_event_with_timing(
        self,
        event_id,
        event_interactions,
        event_trajectories,
        event_segments,
        event_stacks,
        event_hits_back_track,
        event_hits,
        event_light,
    ):
        self.clear_event()
        self.det_point_cloud.data["event"] = event_id
        self.light_point_cloud.data["event"] = event_id

        self.meta["timers"].start("wrangler_process_all")
        self.meta["memory_trackers"].start("wrangler_process_all")

        self.meta["timers"].start("wrangler_process_event_interactions")
        self.meta["memory_trackers"].start("wrangler_process_event_interactions")
        self.process_event_interactions(event_interactions)
        self.meta["timers"].end("wrangler_process_event_interactions")
        self.meta["memory_trackers"].end("wrangler_process_event_interactions")

        self.meta["timers"].start("wrangler_process_event_trajectories")
        self.meta["memory_trackers"].start("wrangler_process_event_trajectories")
        self.process_event_trajectories(event_trajectories)
        self.meta["timers"].end("wrangler_process_event_trajectories")
        self.meta["memory_trackers"].end("wrangler_process_event_trajectories")

        self.meta["timers"].start("wrangler_process_event_stacks")
        self.meta["memory_trackers"].start("wrangler_process_event_stacks")
        self.process_event_stacks(event_stacks)
        self.meta["timers"].end("wrangler_process_event_stacks")
        self.meta["memory_trackers"].end("wrangler_process_event_stacks")

        self.meta["timers"].start("wrangler_process_event_segments")
        self.meta["memory_trackers"].start("wrangler_process_event_segments")
        self.process_event_segments(event_segments)
        self.meta["timers"].end("wrangler_process_event_segments")
        self.meta["memory_trackers"].end("wrangler_process_event_segments")

        self.meta["timers"].start("wrangler_process_event_hits")
        self.meta["memory_trackers"].start("wrangler_process_event_hits")
        self.process_event_hits(event_hits, event_hits_back_track, event_segments)
        self.meta["timers"].end("wrangler_process_event_hits")
        self.meta["memory_trackers"].end("wrangler_process_event_hits")

        self.meta["timers"].start("wrangler_process_event_light")
        self.meta["memory_trackers"].start("wrangler_process_event_light")
        self.process_event_light(event_light)
        self.meta["timers"].end("wrangler_process_event_light")
        self.meta["memory_trackers"].end("wrangler_process_event_light")

        self.meta["timers"].end("wrangler_process_all")
        self.meta["memory_trackers"].end("wrangler_process_all")

    def save_event(self):
        self.det_point_cloud.data["x"] = np.array(self.det_point_cloud.data["x"])
        self.det_point_cloud.data["y"] = np.array(self.det_point_cloud.data["y"])
        self.det_point_cloud.data["z"] = np.array(self.det_point_cloud.data["z"])
        self.det_point_cloud.data["t_drift"] = np.array(
            self.det_point_cloud.data["t_drift"]
        )
        self.det_point_cloud.data["ts_pps"] = np.array(
            self.det_point_cloud.data["ts_pps"]
        )
        self.det_point_cloud.data["Q"] = np.array(self.det_point_cloud.data["Q"])
        self.det_point_cloud.data["E"] = np.array(self.det_point_cloud.data["E"])
        self.det_point_cloud.data["n_photons"] = np.array(
            self.det_point_cloud.data["n_photons"]
        )
        self.det_point_cloud.data["topology_label"] = np.array(
            self.det_point_cloud.data["topology_label"]
        )
        self.det_point_cloud.data["particle_label"] = np.array(
            self.det_point_cloud.data["particle_label"]
        )
        self.det_point_cloud.data["physics_micro_label"] = np.array(
            self.det_point_cloud.data["physics_micro_label"]
        )
        self.det_point_cloud.data["physics_meso_label"] = np.array(
            self.det_point_cloud.data["physics_meso_label"]
        )
        self.det_point_cloud.data["physics_macro_label"] = np.array(
            self.det_point_cloud.data["physics_macro_label"]
        )
        self.det_point_cloud.data["unique_topology_label"] = np.array(
            self.det_point_cloud.data["unique_topology_label"]
        )
        self.det_point_cloud.data["unique_particle_label"] = np.array(
            self.det_point_cloud.data["unique_particle_label"]
        )
        self.det_point_cloud.data["unique_physics_micro_label"] = np.array(
            self.det_point_cloud.data["unique_physics_micro_label"]
        )
        self.det_point_cloud.data["unique_physics_meso_label"] = np.array(
            self.det_point_cloud.data["unique_physics_meso_label"]
        )
        self.det_point_cloud.data["unique_physics_macro_label"] = np.array(
            self.det_point_cloud.data["unique_physics_macro_label"]
        )
        self.det_point_clouds[self.det_point_cloud.data["event"]] = copy.deepcopy(
            self.det_point_cloud
        )

    def save_events(self, simulation_file):
        output_file = simulation_file.replace(".h5", "")
        output_file += ".arrakis_nd.npz"
        meta = {
            "who_created": getpass.getuser(),
            "when_created": datetime.now().strftime("%m-%d-%Y-%H:%M:%S"),
            "where_created": socket.gethostname(),
            "input_file":   simulation_file,
            "det_features": {"x": 0, "y": 1, "z": 2, "Q": 3},
            "mc_features": {"t_drift": 0, "ts_pps": 1, "E": 2, "n_photons": 3},
            "classes": {
                "particle": 0,
                "topology": 1,
                "physics_micro": 2,
                "physics_meso": 3,
                "physics_macro": 4,
            },
            "clusters": {"topology": 0, "particle": 1, "physics": 2},
            "topology_labels": {
                key: value for key, value in classification_labels["topology"].items()
            },
            "particle_labels": {
                key: value for key, value in classification_labels["particle"].items()
            },
            "physics_micro_labels": {
                key: value
                for key, value in classification_labels["physics_micro"].items()
            },
            "physics_meso_labels": {
                key: value
                for key, value in classification_labels["physics_meso"].items()
            },
            "physics_macro_labels": {
                key: value
                for key, value in classification_labels["physics_macro"].items()
            },
            "events":   [int(key) for key in self.det_point_clouds.keys()]

        }
        det_features = np.array(
            [
                np.vstack(
                    (
                        self.det_point_clouds[ii].data["x"],
                        self.det_point_clouds[ii].data["y"],
                        self.det_point_clouds[ii].data["z"],
                        self.det_point_clouds[ii].data["Q"],
                    )
                ).T
                for ii in self.det_point_clouds.keys()
            ],
            dtype=object,
        )
        mc_features = np.array(
            [
                np.vstack(
                    (
                        self.det_point_clouds[ii].data["t_drift"],
                        self.det_point_clouds[ii].data["ts_pps"],
                        self.det_point_clouds[ii].data["E"],
                        self.det_point_clouds[ii].data["n_photons"],
                    )
                ).T
                for ii in self.det_point_clouds.keys()
            ],
            dtype=object,
        )
        classes = np.array(
            [
                np.vstack(
                    (
                        self.det_point_clouds[ii].data["particle_label"],
                        self.det_point_clouds[ii].data["topology_label"],
                        self.det_point_clouds[ii].data["physics_micro_label"],
                        self.det_point_clouds[ii].data["physics_meso_label"],
                        self.det_point_clouds[ii].data["physics_macro_label"],
                    )
                ).T
                for ii in self.det_point_clouds.keys()
            ],
            dtype=object,
        )
        clusters = np.array(
            [
                np.vstack(
                    (
                        self.det_point_clouds[ii].data["unique_particle_label"],
                        self.det_point_clouds[ii].data["unique_topology_label"],
                        self.det_point_clouds[ii].data["unique_physics_micro_label"],
                        self.det_point_clouds[ii].data["unique_physics_meso_label"],
                        self.det_point_clouds[ii].data["unique_physics_macro_label"],
                    )
                ).T
                for ii in self.det_point_clouds.keys()
            ],
            dtype=object,
        )

        np.savez(
            output_file,
            det_features=det_features,
            mc_features=mc_features,
            classes=classes,
            clusters=clusters,
            meta=meta,
        )
        self.clear_point_clouds()

    def process_event_interactions(self, event_interactions):
        for ii, interaction in enumerate(event_interactions):
            vertex_id = interaction['vertex_id']
            self.vertexid_vertex[vertex_id] = interaction['vertex']
            self.vertexid_target[vertex_id] = interaction['target']
            self.vertexid_reaction[vertex_id] = interaction['reaction']
            self.vertexid_pdgcode[vertex_id] = interaction['nu_pdg']
            self.vertexid_label[vertex_id] = PhysicsMacroLabel.Undefined

    def process_event_trajectories(self, event_trajectories):
        if self.wrangler_mode == "map":
            self.process_event_trajectories_map(event_trajectories)
        elif self.wrangler_mode == "numpy":
            self.process_event_trajectories_numpy(event_trajectories)

    def process_event_trajectories_map(self, event_trajectories):
        for ii, particle in enumerate(event_trajectories):
            track_id = particle["traj_id"]
            self.trackid_vertexid[track_id] = particle["vertex_id"]
            self.trackid_parentid[track_id] = particle["parent_id"]
            self.trackid_pdgcode[track_id] = particle["pdg_id"]
            self.trackid_process[track_id] = particle["start_process"]
            self.trackid_subprocess[track_id] = particle["start_subprocess"]
            self.trackid_endprocess[track_id] = particle["end_process"]
            self.trackid_endsubprocess[track_id] = particle["end_subprocess"]
            self.trackid_energy[track_id] = particle["E_end"]  # E_start or E_end?
            self.trackid_tstart[track_id] = particle["t_start"]

            # iterate over daughters
            self.trackid_daughters[track_id] = []
            self.trackid_descendants[track_id] = []

            if particle["parent_id"] != -1:
                self.trackid_daughters[particle["parent_id"]].append(track_id)
                self.trackid_descendants[particle["parent_id"]].append(track_id)
            self.trackid_progeny[track_id] = []
            # iterate over ancestry
            level = 0
            mother = particle["parent_id"]
            temp_track_id = particle["traj_id"]
            ancestry = []
            while mother != -1:
                level += 1
                if track_id not in self.trackid_descendants[mother]:
                    self.trackid_descendants[mother].append(track_id)
                temp_track_id = mother
                ancestry.append(mother)
                mother = self.trackid_parentid[temp_track_id]

                if level > 1 and mother != -1:
                    self.trackid_progeny[mother].append(temp_track_id)

            self.trackid_ancestorlevel[track_id] = level
            self.trackid_ancestry[track_id] = ancestry
            self.trackid_hit[track_id] = []
            self.trackid_segmentid[track_id] = []

    def process_event_trajectories_numpy(self, event_trajectories):
        self.trackid_parentid = event_trajectories["parent_id"]
        self.trackid_pdgcode = event_trajectories["pdg_id"]
        self.trackid_process = event_trajectories["start_process"]
        self.trackid_subprocess = event_trajectories["start_subprocess"]
        self.trackid_endprocess = event_trajectories["end_process"]
        self.trackid_endsubprocess = event_trajectories["end_subprocess"]
        self.trackid_energy = event_trajectories["E_end"]  # E_start or E_end?
        # how to prevent looping for finding daughters?

    def process_event_stacks(self, event_stacks):
        if self.wrangler_mode == "map":
            self.process_event_stacks_map(event_stacks)
        elif self.wrangler_mode == "numpy":
            self.process_event_stacks_numpy(event_stacks)

    def process_event_stacks_map(self, event_stacks):
        for ii, stack in enumerate(event_stacks):
            vertex_id = stack['vertex_id']
            if vertex_id not in self.vertexid_pdgcode.keys():
                self.vertexid_pdgcode[vertex_id] = stack['part_pdg']

    def process_event_stacks_numpy(self, event_stacks):
        pass

    def process_event_segments(self, event_segments):
        if self.wrangler_mode == "map":
            self.process_event_segments_map(event_segments)
        elif self.wrangler_mode == "numpy":
            self.process_event_segments_numpy(event_segments)

    def process_event_segments_map(self, event_segments):
        for ii, segment in enumerate(event_segments):
            self.trackid_segmentid[segment["traj_id"]].append(
                segment["segment_id"]
            )  # segment_ids belonging to a track_id
            self.segmentid_trackid[segment["segment_id"]] = segment[
                "traj_id"
            ]  # track_id belonging to a segment_id
            self.segmentid_hit[
                segment["segment_id"]
            ] = []  # hit_ids belonging to a segment_id

    def process_event_segments_numpy(self, event_segments):
        pass

    def process_event_hits(self, event_hits, event_hits_back_track, event_segments):
        if self.wrangler_mode == "map":
            self.process_event_hits_map(event_hits, event_hits_back_track, event_segments)
        elif self.wrangler_mode == "numpy":
            self.process_event_hits_numpy(event_hits, event_hits_back_track, event_segments)

    def process_event_hits_map(
        self,
        event_hits,
        event_hits_back_track,
        event_segments
    ):
        segment_ids = np.array(event_hits_back_track["segment_id"])
        segment_fractions = np.array(event_hits_back_track["fraction"])
        segment_track_ids = []
        for segment_id in segment_ids:
            segment_id = segment_id[segment_id != 0]
            segment_id_indices = np.where(np.isin(event_segments["segment_id"], segment_id))
            segment_track_ids.append(event_segments["traj_id"][segment_id_indices])
        segment_track_ids = np.array(segment_track_ids, dtype=object)
        # Create a dictionary mapping from segment IDs to number of photons
        segment_id_to_n_photons = {seg['segment_id']: seg['n_photons'] for seg in event_segments}

        # Create a list of number of photons for each hit
        n_photons_per_hit = [
            [segment_id_to_n_photons[seg_id] for seg_id in seg_ids if seg_id in segment_id_to_n_photons]
            for seg_ids in event_hits_back_track["segment_id"]
        ]
        total_photons_per_hit = [np.sum(hit) for hit in n_photons_per_hit]
        self.det_point_cloud.add_event(
            event_hits["x"],
            event_hits["y"],
            event_hits["z"],
            event_hits["t_drift"],
            event_hits["ts_pps"],
            event_hits["Q"],
            event_hits["E"],
            total_photons_per_hit,
            segment_ids,
            segment_fractions,
            segment_track_ids,
        )
        for ii, hit in enumerate(event_hits):
            segment_ids = event_hits_back_track["segment_id"][ii]
            segment_fractions = event_hits_back_track["fraction"][ii]
            for segmentid in segment_ids[(segment_ids != 0)]:
                if segmentid in self.segmentid_hit.keys():
                    self.segmentid_hit[segmentid].append(ii)
                    self.trackid_hit[self.segmentid_trackid[segmentid]].append(ii)

    def process_event_hits_numpy(
        self,
        event_hits,
        event_hits_back_track,
        event_segments
    ):
        segment_ids = event_hits_back_track["segment_id"]
        segment_fractions = event_hits_back_track["fraction"]
        # Create a dictionary mapping from segment IDs to number of photons
        segment_id_to_n_photons = {seg['segment_id']: seg['n_photons'] for seg in event_segments}

        # Create a list of number of photons for each hit
        n_photons_per_hit = [
            [segment_id_to_n_photons[seg_id] for seg_id in seg_ids if seg_id in segment_id_to_n_photons]
            for seg_ids in event_hits_back_track["segment_id"]
        ]
        total_photons_per_hit = [np.sum(hit) for hit in n_photons_per_hit]
        self.det_point_cloud.add_event(
            event_hits["x"],
            event_hits["y"],
            event_hits["z"],
            event_hits["t_drift"],
            event_hits["ts_pps"],
            event_hits["Q"],
            event_hits["E"],
            total_photons_per_hit,
            segment_ids,
            segment_fractions,
        )
        for ii, hit in enumerate(event_hits):
            for segmentid in segment_ids[(segment_ids != 0)]:
                if segmentid in self.segmentid_hit.keys():
                    self.segmentid_hit[segmentid].append(ii)
                    self.trackid_hit[self.segmentid_trackid[segmentid]].append(ii)

    def process_event_light(
            self,
            event_light,
    ):
        self.light_point_cloud.add_event(
            event_light["tpc"],
            event_light["det"],
            event_light["sample_idx"],
            event_light["max"],
            event_light["fwhm_spline"],
            event_light["sum"],
            [],                         # TODO: segment_ids from truth
        )

    def get_total_hit_energy_map(self, hits):
        """Get total energy from a list of hits

        Args:
            hits (_type_): _description_

        Returns:
            _type_: _description_
        """
        energy = 0.0
        for hit in hits:
            energy += self.det_point_cloud.data["E"][hit]
        return energy

    def get_total_hit_energy_numpy(self, hits):
        pass

    def get_primaries_generator_label(self, label):
        pass

    def get_primaries_pdg_code(self, pdg):
        primaries = []
        for track_id, parent_id in self.trackid_parentid.items():
            if parent_id == -1 and self.trackid_pdgcode[track_id] == pdg:
                primaries.append(track_id)
        return primaries

    def get_primaries_abs_pdg_code(self, pdg):
        primaries = []
        for track_id, parent_id in self.trackid_parentid.items():
            if parent_id == -1 and abs(self.trackid_pdgcode[track_id]) == abs(pdg):
                primaries.append(track_id)
        return primaries

    def get_hits_trackid(self, trackids):
        trackids_np = np.array(trackids).astype(int).flatten()
        hits = [self.trackid_hit[track_id] for track_id in trackids_np]
        return hits

    def get_segments_trackid(self, trackids):
        trackids_np = np.array(trackids).astype(int).flatten()
        segments = [self.trackid_segmentid[track_id] for track_id in trackids_np]
        return segments

    def get_trackid_pdg_code(self, pdg):
        trackid = []
        for track_id, pdg_code in self.trackid_pdgcode.items():
            if pdg_code == pdg:
                trackid.append(track_id)
        return trackid

    def get_trackid_abs_pdg_code(self, pdg):
        trackid = []
        for track_id, pdg_code in self.trackid_pdgcode.items():
            if abs(pdg_code) == abs(pdg):
                trackid.append(track_id)
        return trackid

    def get_parentid_trackid(self, trackids):
        parentid = [self.trackid_parentid[track_id] for track_id in trackids]
        return parentid

    def get_daughters_trackid(self, trackids):
        daughters = [self.trackid_daughters[track_id] for track_id in trackids]
        return daughters

    def get_descendants_trackid(self, trackids):
        descendants = [self.trackid_descendants[track_id] for track_id in trackids]
        return descendants

    def get_tstart_trackid(self, trackids):
        tstart = [self.trackid_tstart[track_id] for track_id in trackids]
        return tstart

    def filter_trackid_not_pdg_code(self, trackids, pdg):
        trackid = [
            track_id for track_id in trackids if self.trackid_pdgcode[track_id] != pdg
        ]
        return trackid

    def filter_trackid_pdg_code(self, trackids, pdg):
        trackid = [
            track_id for track_id in trackids if self.trackid_pdgcode[track_id] == pdg
        ]
        return trackid

    def filter_trackid_not_abs_pdg_code(self, trackids, pdg):
        trackid = [
            track_id
            for track_id in trackids
            if abs(self.trackid_pdgcode[track_id]) != abs(pdg)
        ]
        return trackid

    def filter_trackid_abs_pdg_code(self, trackids, pdg):
        trackid = [
            track_id
            for track_id in trackids
            if abs(self.trackid_pdgcode[track_id]) == abs(pdg)
        ]
        return trackid

    def filter_trackid_not_process(self, trackids, process):
        trackid = [
            track_id
            for track_id in trackids
            if self.trackid_process[track_id] != process
        ]
        return trackid

    def filter_trackid_process(self, trackids, process):
        trackid = [
            track_id
            for track_id in trackids
            if self.trackid_process[track_id] == process
        ]
        return trackid

    def filter_trackid_not_subprocess(self, trackids, subprocess):
        trackid = [
            track_id
            for track_id in trackids
            if self.trackid_subprocess[track_id] != subprocess
        ]
        return trackid

    def filter_trackid_not_process_and_subprocess(self, trackids, process, subprocess):
        trackid = [
            track_id
            for track_id in trackids
            if not (
                self.trackid_process[track_id] == process
                and self.trackid_subprocess[track_id] == subprocess
            )
        ]
        return trackid

    def filter_trackid_subprocess(self, trackids, subprocess):
        trackid = [
            track_id
            for track_id in trackids
            if self.trackid_subprocess[track_id] == subprocess
        ]
        return trackid

    def filter_trackid_not_endprocess(self, trackids, endprocess):
        trackid = [
            track_id
            for track_id in trackids
            if self.trackid_endprocess[track_id] != endprocess
        ]
        return trackid

    def filter_trackid_endprocess(self, trackids, endprocess):
        trackid = [
            track_id
            for track_id in trackids
            if self.trackid_endprocess[track_id] == endprocess
        ]
        return trackid

    def filter_trackid_not_endsubprocess(self, trackids, endsubprocess):
        trackid = [
            track_id
            for track_id in trackids
            if self.trackid_endsubprocess[track_id] != endsubprocess
        ]
        return trackid

    def filter_trackid_endsubprocess(self, trackids, endsubprocess):
        trackid = [
            track_id
            for track_id in trackids
            if self.trackid_endsubprocess[track_id] == endsubprocess
        ]
        return trackid

    def get_index_trackid(self, hit, track_id):
        index = np.where(self.det_point_cloud.data["segment_track_ids"][hit] == track_id)
        return index[0]

    def get_hit_distance(self, hit1, hit2):
        x1 = self.det_point_cloud.data['x'][hit1]
        x2 = self.det_point_cloud.data['x'][hit2]

        y1 = self.det_point_cloud.data['y'][hit1]
        y2 = self.det_point_cloud.data['y'][hit2]

        z1 = self.det_point_cloud.data['z'][hit1]
        z2 = self.det_point_cloud.data['z'][hit2]

        dist = np.sqrt(
            (x1 - x2) * (x1 - x2) +
            (y1 - y2) * (y1 - y2) +
            (z1 - z2) * (z1 - z2)
        )
        return dist

    def get_closest_hit(self, hit, hits):
        closest_hit = -1
        closest_hit_distance = 10e10
        for h in hits:
            hit_distance = self.get_hit_distance(hit, h)
            if hit_distance < closest_hit_distance:
                closest_hit = h
                closest_hit_distance = hit_distance
        return closest_hit
