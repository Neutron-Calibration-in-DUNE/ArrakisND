"""
Simluation Labeling Logic

Developers: Nicholas Carrara        [nmcarrara@ucdavis.edu]
            Marjolein van Nuland    [mnuland@nikhef.nl]

ChangeLog:  12/17/2023 - started putting together shower logic.
"""
import numpy as np
from arrakis_nd.utils.logger import Logger
from arrakis_nd.dataset.common import (
    ProcessType,
    SubProcessType,
    ParticleLabel,
    TopologyLabel,
    PhysicsMicroLabel,
    PhysicsMesoLabel,
    PhysicsMacroLabel,
)
from arrakis_nd.utils.utils import ResetableIterator, remove_sublist


class SimulationLabelingLogic:
    """
    Simulation labeling logic conducts mid-level and high-level
    variable construction from low-level information, such as
    pdg code, particle heirarchy, process and subprocess, etc.
    """

    def __init__(
        self,
        name: str = "",
        config: dict = {},
        meta: dict = {},
    ):
        self.name = name + "_simulation_labeling_logic"
        self.config = config
        self.meta = meta

        if "device" in self.meta:
            self.device = self.meta["device"]
            self.gpu = self.meta["gpu"]
        else:
            self.device = "cpu"
            self.gpu = False
        if meta["verbose"]:
            self.logger = Logger(self.name, output="both", file_mode="w")
        else:
            self.logger = Logger(self.name, level="warning", file_mode="w")

        self.unique_topology = ResetableIterator()
        self.unique_physics_micro = ResetableIterator()
        self.unique_physics_meso = ResetableIterator()
        self.unique_physics_macro = ResetableIterator()

        self.parse_config()

    def parse_config(self):
        self.check_config()

    def check_config(self):
        if "simulation_wrangler" not in self.meta.keys():
            self.logger.error("no simulation_wrangler defined in meta!")
        self.simulation_wrangler = self.meta["simulation_wrangler"]
        if "shower_threshold" not in self.config.keys():
            self.logger.warn('no shower_threshold in config! setting to "20"')
            self.config["shower_threshold"] = 20
        self.shower_threshold = self.config["shower_threshold"]
        if "ionization_influence_distance" not in self.config.keys():
            self.logger.warn('no ionization_influence_distance specified in config! setting to 4.')
            self.config["ionization_influence_distance"] = 4
        self.ionization_influence_distance = self.config["ionization_influence_distance"]
        if "debug" not in self.config.keys():
            self.logger.warn('debug not specified in config! setting to "False"')
            self.config["debug"] = False
        self.debug = self.config["debug"]

    def process_event(self):
        """
        Go through each of the individual logic functions
        to process mid-level and high-level variables.

        The list here should be exhaustive.  In debug mode, a
        set of check functions are run which gathers statistics
        on any failures, as well as timing and memory usage
        information for each function.
        """
        if self.debug:
            self._process_event_with_timing()
        else:
            self._process_event_without_timing()

    def process_light_event(self):
        pass

    def _process_event_without_timing(self):
        self.unique_topology.reset
        self.unique_physics_micro.reset
        self.unique_physics_meso.reset
        self.unique_physics_macro.reset
        self.process_mc_truth()
        self.process_electrons()
        self.process_positrons()
        self.process_gammas()
        self.process_muons()
        self.process_pion0s()
        self.process_pions()
        self.process_kaon0s()
        self.process_kaons()
        self.process_protons()
        self.process_neutrons()
        self.process_nuclear_recoils()
        self.process_electron_recoils()
        self.process_baryons()
        self.process_ar39()
        self.process_ar42()
        self.process_kr85()
        self.process_rn222()
        self.process_cosmics()
        self.process_primaries()
        self.clean_up_labels()

    def _process_event_with_timing(self):
        self.unique_topology.reset
        self.unique_physics_micro.reset
        self.unique_physics_meso.reset
        self.unique_physics_macro.reset

        self.meta["timers"].start("logic_process_all")
        self.meta["memory_trackers"].start("logic_process_all")

        self.meta["timers"].start("logic_process_mc_truth")
        self.meta["memory_trackers"].start("logic_process_mc_truth")
        self.process_mc_truth()
        self.meta["timers"].end("logic_process_mc_truth")
        self.meta["memory_trackers"].end("logic_process_mc_truth")

        self.meta["timers"].start("logic_process_electrons")
        self.meta["memory_trackers"].start("logic_process_electrons")
        self.process_electrons()
        self.meta["timers"].end("logic_process_electrons")
        self.meta["memory_trackers"].end("logic_process_electrons")

        self.meta["timers"].start("logic_process_positrons")
        self.meta["memory_trackers"].start("logic_process_positrons")
        self.process_positrons()
        self.meta["timers"].end("logic_process_positrons")
        self.meta["memory_trackers"].end("logic_process_positrons")

        self.meta["timers"].start("logic_process_gammas")
        self.meta["memory_trackers"].start("logic_process_gammas")
        self.process_gammas()
        self.meta["timers"].end("logic_process_gammas")
        self.meta["memory_trackers"].end("logic_process_gammas")

        self.meta["timers"].start("logic_process_muons")
        self.meta["memory_trackers"].start("logic_process_muons")
        self.process_muons()
        self.meta["timers"].end("logic_process_muons")
        self.meta["memory_trackers"].end("logic_process_muons")

        self.meta["timers"].start("logic_process_pion0s")
        self.meta["memory_trackers"].start("logic_process_pion0s")
        self.process_pion0s()
        self.meta["timers"].end("logic_process_pion0s")
        self.meta["memory_trackers"].end("logic_process_pion0s")

        self.meta["timers"].start("logic_process_pions")
        self.meta["memory_trackers"].start("logic_process_pions")
        self.process_pions()
        self.meta["timers"].end("logic_process_pions")
        self.meta["memory_trackers"].end("logic_process_pions")

        self.meta["timers"].start("logic_process_kaon0s")
        self.meta["memory_trackers"].start("logic_process_kaon0s")
        self.process_kaon0s()
        self.meta["timers"].end("logic_process_kaon0s")
        self.meta["memory_trackers"].end("logic_process_kaon0s")

        self.meta["timers"].start("logic_process_kaons")
        self.meta["memory_trackers"].start("logic_process_kaons")
        self.process_kaons()
        self.meta["timers"].end("logic_process_kaons")
        self.meta["memory_trackers"].end("logic_process_kaons")

        self.meta["timers"].start("logic_process_protons")
        self.meta["memory_trackers"].start("logic_process_protons")
        self.process_protons()
        self.meta["timers"].end("logic_process_protons")
        self.meta["memory_trackers"].end("logic_process_protons")

        self.meta["timers"].start("logic_process_neutrons")
        self.meta["memory_trackers"].start("logic_process_neutrons")
        self.process_neutrons()
        self.meta["timers"].end("logic_process_neutrons")
        self.meta["memory_trackers"].end("logic_process_neutrons")

        self.meta["timers"].start("logic_process_nuclear_recoils")
        self.meta["memory_trackers"].start("logic_process_nuclear_recoils")
        self.process_nuclear_recoils()
        self.meta["timers"].end("logic_process_nuclear_recoils")
        self.meta["memory_trackers"].end("logic_process_nuclear_recoils")

        self.meta["timers"].start("logic_process_electron_recoils")
        self.meta["memory_trackers"].start("logic_process_electron_recoils")
        self.process_electron_recoils()
        self.meta["timers"].end("logic_process_electron_recoils")
        self.meta["memory_trackers"].end("logic_process_electron_recoils")

        self.meta["timers"].start("logic_process_baryons")
        self.meta["memory_trackers"].start("logic_process_baryons")
        self.process_baryons()
        self.meta["timers"].end("logic_process_baryons")
        self.meta["memory_trackers"].end("logic_process_baryons")

        self.meta["timers"].start("logic_process_ar39")
        self.meta["memory_trackers"].start("logic_process_ar39")
        self.process_ar39()
        self.meta["timers"].end("logic_process_ar39")
        self.meta["memory_trackers"].end("logic_process_ar39")

        self.meta["timers"].start("logic_process_ar42")
        self.meta["memory_trackers"].start("logic_process_ar42")
        self.process_ar42()
        self.meta["timers"].end("logic_process_ar42")
        self.meta["memory_trackers"].end("logic_process_ar42")

        self.meta["timers"].start("logic_process_kr85")
        self.meta["memory_trackers"].start("logic_process_kr85")
        self.process_kr85()
        self.meta["timers"].end("logic_process_kr85")
        self.meta["memory_trackers"].end("logic_process_kr85")

        self.meta["timers"].start("logic_process_rn222")
        self.meta["memory_trackers"].start("logic_process_rn222")
        self.process_rn222()
        self.meta["timers"].end("logic_process_rn222")
        self.meta["memory_trackers"].end("logic_process_rn222")

        self.meta["timers"].start("logic_process_cosmics")
        self.meta["memory_trackers"].start("logic_process_cosmics")
        self.process_cosmics()
        self.meta["timers"].end("logic_process_cosmics")
        self.meta["memory_trackers"].end("logic_process_cosmics")

        self.meta["timers"].start("logic_clean_up_labels")
        self.meta["memory_trackers"].start("logic_clean_up_labels")
        self.clean_up_labels()
        self.meta["timers"].end("logic_clean_up_labels")
        self.meta["memory_trackers"].end("logic_clean_up_labels")

        self.meta["timers"].start("logic_process_all")
        self.meta["memory_trackers"].end("logic_process_all")

        self.meta["timers"].start("logic_check_labels")
        self.meta["memory_trackers"].start("logic_check_labels")
        self.check_labels()
        self.meta["timers"].end("logic_check_labels")
        self.meta["memory_trackers"].end("logic_check_labels")

    def check_labels(self):
        for particle in self.simulation_wrangler.trackid_hit.keys():
            hits = self.simulation_wrangler.trackid_hit[particle]
            topology_labels = [
                self.simulation_wrangler.det_point_cloud.data["topology_label"][hit]
                for hit in hits
            ]
            if -1 in topology_labels:
                self.logger.debug('###################################################')
                self.logger.debug(f'## Missing hit labels for hits: {hits}')
                self.logger.debug(f'## Topology labels:             {topology_labels}')
                self.simulation_wrangler.print_particle_data(particle)

    def clean_up_labels(self):
        # take care of ions which are not specifically itemized in the particle list
        for particle in self.simulation_wrangler.trackid_hit.keys():
            if self.simulation_wrangler.trackid_pdgcode[particle] > 2212:
                hits = self.simulation_wrangler.trackid_hit[particle]
                topology_labels = np.array([
                    self.simulation_wrangler.det_point_cloud.data["topology_labels"][hit]
                    for hit in hits
                ], dtype=object)
                if len(topology_labels) > 0:
                    topology_labels = np.concatenate(topology_labels)
                if -1 in topology_labels.flatten():
                    ion_hits = self.simulation_wrangler.get_hits_trackid(particle)
                    ion_segments = self.simulation_wrangler.get_segments_trackid(particle)
                    self.simulation_wrangler.set_hit_labels(
                        ion_hits[0],
                        ion_segments[0],
                        particle,
                        TopologyLabel.Blip,
                        PhysicsMicroLabel.HIPIonization,
                        PhysicsMesoLabel.NuclearRecoil,
                        next(self.unique_topology),
                        next(self.unique_physics_micro),
                        next(self.unique_physics_meso),
                    )
                for hit in hits:
                    self.simulation_wrangler.det_point_cloud.data['particle_label'][hit] = ParticleLabel.Ion.value
        # infer ionization from surrounding hits, and reconcile hits which have multiple
        # contributing particles.
        for particle in self.simulation_wrangler.trackid_hit.keys():
            hits = self.simulation_wrangler.trackid_hit[particle]
            for hit in hits:
                if self.simulation_wrangler.det_point_cloud.data["topology_label"][hit] == TopologyLabel.Undefined.value:
                    if self.simulation_wrangler.trackid_subprocess[particle] == SubProcessType.Ionization.value:
                        self.infer_ionization(particle, hit)
                else:
                    if len(self.simulation_wrangler.det_point_cloud.data["topology_labels"][hit]) > 1:
                        self.reconcile_hit(particle, hit)

    def infer_ionization(
        self,
        particle,
        hit
    ):
        """
        This function cleans up points which are generated from electron ionization
        but which are not caught by any of the logic.
        """
        parent = self.simulation_wrangler.trackid_parentid[particle]
        if parent != -1:
            parent_hits = self.simulation_wrangler.get_hits_trackid(parent)
            closest_hit = self.simulation_wrangler.get_closest_hit(hit, parent_hits[0])
            hit_distance = self.simulation_wrangler.get_hit_distance(hit, closest_hit)
            pdg_code = self.simulation_wrangler.trackid_pdgcode[particle]
            self.simulation_wrangler.det_point_cloud.data['particle_label'][hit] = pdg_code
            self.simulation_wrangler.det_point_cloud.data['unique_particle_label'][hit] = particle
            if hit_distance <= self.ionization_influence_distance:
                self.simulation_wrangler.copy_hit_labels(hit, closest_hit)
            else:
                self.simulation_wrangler.set_hit_labels(
                    [hit],
                    [],
                    particle,
                    TopologyLabel.Blip,
                    PhysicsMicroLabel.ElectronIonization,
                    PhysicsMesoLabel.LowEnergyIonization,
                    next(self.unique_topology),
                    next(self.unique_physics_micro),
                    next(self.unique_physics_meso),
                )

    def reconcile_hit(
        self,
        particle,
        hit
    ):
        """
        This function takes care of hits that have multiple contributions from
        different particles.  The heirarchy is that tracks take precendence
        over all other type (blips/showers).  So first check if there are any track
        labels.  Then, adjust the other labels accordingly
        """
        """
        Multiple possible situations here, either:
            1. a single track is crossing another object (switch all labels to the track)
            2. two or more tracks are crossing each other (pick the largest contributor)
            3. a delta is coming off a track (assign labels to the track not the delta)
        """
        segment_fractions = self.simulation_wrangler.det_point_cloud.data['segment_fractions'][hit]
        largest_fraction = np.argmax(segment_fractions)

        # set labels to largest object labels
        topology_label = self.simulation_wrangler.det_point_cloud.data['topology_labels'][hit][largest_fraction]
        particle_label = self.simulation_wrangler.det_point_cloud.data['particle_labels'][hit][largest_fraction]
        physics_micro_label = self.simulation_wrangler.det_point_cloud.data['physics_micro_labels'][hit][largest_fraction]
        physics_meso_label = self.simulation_wrangler.det_point_cloud.data['physics_meso_labels'][hit][largest_fraction]
        physics_macro_label = self.simulation_wrangler.det_point_cloud.data['physics_macro_labels'][hit][largest_fraction]
        unique_topology_label = self.simulation_wrangler.det_point_cloud.data['unique_topology_labels'][hit][largest_fraction]
        unique_particle_label = self.simulation_wrangler.det_point_cloud.data['unique_particle_labels'][hit][largest_fraction]
        unique_physics_micro_label = self.simulation_wrangler.det_point_cloud.data['unique_physics_micro_labels'][hit][
            largest_fraction
        ]
        unique_physics_meso_label = self.simulation_wrangler.det_point_cloud.data['unique_physics_meso_labels'][hit][
            largest_fraction
        ]
        unique_physics_macro_label = self.simulation_wrangler.det_point_cloud.data['unique_physics_macro_labels'][hit][
            largest_fraction
        ]

        self.simulation_wrangler.det_point_cloud.data['topology_label'][hit] = topology_label
        self.simulation_wrangler.det_point_cloud.data['particle_label'][hit] = particle_label
        self.simulation_wrangler.det_point_cloud.data['physics_micro_label'][hit] = physics_micro_label
        self.simulation_wrangler.det_point_cloud.data['physics_meso_label'][hit] = physics_meso_label
        self.simulation_wrangler.det_point_cloud.data['physics_macro_label'][hit] = physics_macro_label
        self.simulation_wrangler.det_point_cloud.data['unique_topology_label'][hit] = unique_topology_label
        self.simulation_wrangler.det_point_cloud.data['unique_particle_label'][hit] = unique_particle_label
        self.simulation_wrangler.det_point_cloud.data['unique_physics_micro_label'][hit] = unique_physics_micro_label
        self.simulation_wrangler.det_point_cloud.data['unique_physics_meso_label'][hit] = unique_physics_meso_label
        self.simulation_wrangler.det_point_cloud.data['unique_physics_macro_label'][hit] = unique_physics_macro_label

    def find_closest_hit(
        self,
        position,
        hits:   list = [],
    ):
        # add in ability to search through all hits if hits list is empty
        closest_hit = hits[0]
        closest_distance = 10e10
        for ii in range(len(hits)):
            hit_x = self.simulation_wrangler.det_point_cloud.data["x"][hits[ii]]
            hit_y = self.simulation_wrangler.det_point_cloud.data["y"][hits[ii]]
            hit_z = self.simulation_wrangler.det_point_cloud.data["z"][hits[ii]]
            hit_distance = np.sqrt(
                (position[0] - hit_x)**2 +
                (position[1] - hit_y)**2 +
                (position[2] - hit_z)**2
            )
            if hit_distance < closest_distance:
                closest_hit = hits[ii]
                closest_distance = hit_distance
        return closest_hit

    def process_mc_truth(self):
        """
        This processes the macro labels for each hit.  We check vertex ids to see
        what kind of interaction is associated with each hit.
        """
        for vertex_id, reaction in self.simulation_wrangler.vertexid_reaction.items():
            if abs(reaction) in [1, 2, 3, 4, 5]:
                if abs(self.simulation_wrangler.vertexid_pdgcode[vertex_id] == 12):
                    self.simulation_wrangler.vertexid_label[vertex_id] = PhysicsMacroLabel.CCNue
                elif abs(self.simulation_wrangler.vertexid_pdgcode[vertex_id] == 14):
                    self.simulation_wrangler.vertexid_label[vertex_id] = PhysicsMacroLabel.CCNuMu
            else:
                self.simulation_wrangler.vertexid_label[vertex_id] = PhysicsMacroLabel.NC

    def process_track(
        self,
        particle: int = 0,
        unique_topology: int = 0
    ):
        """
        """
        particle_hits = self.simulation_wrangler.trackid_hit[particle]
        if len(particle_hits) > 0:
            hit_start = self.find_closest_hit(
                self.simulation_wrangler.trackid_xyz_start[particle],
                particle_hits
            )
            hit_end = self.find_closest_hit(
                self.simulation_wrangler.trackid_xyz_end[particle],
                particle_hits
            )
            self.simulation_wrangler.det_point_cloud.data["track_begin"][hit_start] = 1
            self.simulation_wrangler.det_point_cloud.data["track_end"][hit_end] = 1
            self.simulation_wrangler.track.add_mc_track(
                track_id=particle,
                E_start=self.simulation_wrangler.trackid_energy_start[particle],
                xyz_start=self.simulation_wrangler.trackid_xyz_start[particle],
                pxyz_start=self.simulation_wrangler.trackid_momentum_start[particle],
                t_start=self.simulation_wrangler.trackid_t_start[particle],
                hit_start=hit_start,
                E_end=self.simulation_wrangler.trackid_energy_end[particle],
                xyz_end=self.simulation_wrangler.trackid_xyz_end[particle],
                pxyz_end=self.simulation_wrangler.trackid_momentum_end[particle],
                t_end=self.simulation_wrangler.trackid_t_end[particle],
                hit_end=hit_end,
                track_length=0,
                dEdx=0,
                Q_total=0,
                hits=particle_hits,
            )

    def process_delta(
        self,
        particle: int = 0,
        unique_topology: int = 0
    ):
        """
        """
        particle_hits = self.simulation_wrangler.trackid_hit[particle]
        if len(particle_hits) > 0:
            hit_start = self.find_closest_hit(
                self.simulation_wrangler.trackid_xyz_start[particle],
                particle_hits
            )
            hit_end = self.find_closest_hit(
                self.simulation_wrangler.trackid_xyz_end[particle],
                particle_hits
            )
            self.simulation_wrangler.det_point_cloud.data["delta_begin"][hit_start] = 1
            self.simulation_wrangler.det_point_cloud.data["delta_end"][hit_end] = 1
            self.simulation_wrangler.track.add_mc_track(
                track_id=particle,
                E_start=self.simulation_wrangler.trackid_energy_start[particle],
                xyz_start=self.simulation_wrangler.trackid_xyz_start[particle],
                pxyz_start=self.simulation_wrangler.trackid_momentum_start[particle],
                t_start=self.simulation_wrangler.trackid_t_start[particle],
                hit_start=hit_start,
                E_end=self.simulation_wrangler.trackid_energy_end[particle],
                xyz_end=self.simulation_wrangler.trackid_xyz_end[particle],
                pxyz_end=self.simulation_wrangler.trackid_momentum_end[particle],
                t_end=self.simulation_wrangler.trackid_t_end[particle],
                hit_end=hit_end,
                track_length=0,
                dEdx=0,
                Q_total=0,
                hits=particle_hits,
            )

    def process_showers(self, particle: int = 0, unique_topology: int = 0):
        """
        Showers can come in two flavors, electromagnetic and hadronic.  Only
        dealing with the electromagnetic version for now.  The methodology
        is to assume that each trackid is coming from an electron, positron
        or gamma, but we should check this just in case.  The full shower,
        if there is one, will be associated back to whichever particle
        created the shower.

        Showers currently consist of
            (1) topology label          - TopologyLabel.Shower = 3,
            (3) physics_micro labels    - PhysicsMicroLabel.ElectronIonization = 3
                                          PhysicsMicroLabel.GammaCompton = 4
                                          PhysicsMicroLabel.GammaConversion = 5
            (3) physics_meso labels     - PhysicsMesoLabel.ElectronShower = 5
                                          PhysicsMesoLabel.PositronShower = 6
                                          PhysicsMesoLabel.PhotonShower = 7
            (4) physics_macro labels    - PhysicsMacroLabel.CCNue = 1
                                          PhysicsMacroLabel.CCNuMu = 2
                                          PhysicsMacroLabel.NC = 3
                                          PhysicsMacroLabel.Cosmics = 5

        The physics_micro labels are determined by performing a series of steps.
        First, we look to see if there are any gamma conversion events within
        the hierarchy.  This is denoted as subprocess 14 within edep-sim.
        Each set of descendants from these gammas will get the physics_micro
        label PhysicsMicroLabel.GammaConversion.  We then look for any compton
        scatters or photo-electric effect.  Comptons which have the same parent
        are linked by unique_physics_micro numbers, while photo-electric effect
        have independent ones.  The photo-electric effect gets the labels
        PhysicsMicroLabel.ElectronIonization.

        The physics_meso labels are determined simply by the type of originating
        particle for the shower.  If there is no shower, then the physics_meso
        label ...

        The physics_macro labels are determined by a separate function.

        """
        # grab descendants by type
        particle_subprocess = self.simulation_wrangler.trackid_subprocess[particle]
        particle_descendants = self.simulation_wrangler.trackid_descendants[particle]

        ionization_descendants = self.simulation_wrangler.filter_trackid_subprocess(
            particle_descendants, SubProcessType.Ionization.value
        )
        bremsstrahlung_descendants = self.simulation_wrangler.filter_trackid_subprocess(
            particle_descendants, SubProcessType.Bremsstrahlung.value
        )
        pairprodcharge_descendants = self.simulation_wrangler.filter_trackid_subprocess(
            particle_descendants, SubProcessType.PairProdByCharge.value
        )
        annihilation_descendants = self.simulation_wrangler.filter_trackid_subprocess(
            particle_descendants, SubProcessType.Annihilation.value
        )
        photoelectric_descendants = self.simulation_wrangler.filter_trackid_subprocess(
            particle_descendants, SubProcessType.PhotoElectricEffect.value
        )
        compton_descendants = self.simulation_wrangler.filter_trackid_subprocess(
            particle_descendants, SubProcessType.ComptonScattering.value
        )
        conversion_descendants = self.simulation_wrangler.filter_trackid_subprocess(
            particle_descendants, SubProcessType.GammaConversion.value
        )
        scintillation_descendants = self.simulation_wrangler.filter_trackid_subprocess(
            particle_descendants, SubProcessType.Scintillation.value
        )

        # count number of each subprocess, including the particle itself
        if particle_subprocess == SubProcessType.Ionization.value:
            ionization_descendants.append(particle)
        elif particle_subprocess == SubProcessType.Bremsstrahlung.value:
            bremsstrahlung_descendants.append(particle)
        elif particle_subprocess == SubProcessType.PairProdByCharge.value:
            pairprodcharge_descendants.append(particle)
        elif particle_subprocess == SubProcessType.Annihilation.value:
            annihilation_descendants.append(particle)
        elif particle_subprocess == SubProcessType.PhotoElectricEffect.value:
            photoelectric_descendants.append(particle)
        elif particle_subprocess == SubProcessType.ComptonScattering.value:
            compton_descendants.append(particle)
        elif particle_subprocess == SubProcessType.GammaConversion.value:
            conversion_descendants.append(particle)
        elif particle_subprocess == SubProcessType.Scintillation.value:
            scintillation_descendants.append(particle)

        num_conversion = len(conversion_descendants)
        num_bremsstrahlung = len(bremsstrahlung_descendants)
        # shower or no shower?

        """
        Check times for conversion and bremsstrahlung.  If there are conversions,
        find the earliest conversion, and get it's parents pdg code.  Same
        with bremsstrahlung.  If the particle is a primary, use its pdg code instead
        """
        if num_conversion > 0:
            conversion_times = self.simulation_wrangler.get_t_start_trackid(conversion_descendants)
            earliest_conversion_time = min(conversion_times)
            earliest_conversion = conversion_descendants[conversion_times.index(earliest_conversion_time)]
            if self.simulation_wrangler.trackid_parentid[earliest_conversion] == -1:
                earliest_conversion_pdg_code = self.simulation_wrangler.trackid_pdgcode[earliest_conversion]
            else:
                earliest_conversion_pdg_code = self.simulation_wrangler.trackid_pdgcode[
                    self.simulation_wrangler.trackid_parentid[earliest_conversion]
                ]
        else:
            earliest_conversion_time = 10e10
            earliest_conversion = None
            earliest_conversion_pdg_code = None
        if num_bremsstrahlung > 0:
            bremsstrahlung_times = self.simulation_wrangler.get_t_start_trackid(bremsstrahlung_descendants)
            earliest_bremsstrahlung_time = min(bremsstrahlung_times)
            earliest_bremsstrahlung = bremsstrahlung_descendants[bremsstrahlung_times.index(earliest_bremsstrahlung_time)]
            if self.simulation_wrangler.trackid_parentid[earliest_bremsstrahlung] == -1:
                earliest_bremsstrahlung_pdg_code = self.simulation_wrangler.trackid_pdgcode[earliest_bremsstrahlung]
            else:
                earliest_bremsstrahlung_pdg_code = self.simulation_wrangler.trackid_pdgcode[
                    self.simulation_wrangler.trackid_parentid[earliest_bremsstrahlung]
                ]
        else:
            earliest_bremsstrahlung_time = 10e10
            earliest_bremsstrahlung = None
            earliest_bremsstrahlung_pdg_code = None

        """
        Now check what the pdg codes associated to the earliest of either
        conversion or bremsstrahlung are.  If the conversion is earliest,
        and there is a pi0 associated to it, then label as a pi0decay,
        otherwise we have a photon shower.

        If there are no bremss or conversions, then check to see
        how much energy is associated with the incoming particle,
        and use that to determine whether it should be a shower
        or not.
        """
        if (num_conversion > 0):
            topology = TopologyLabel.Shower
            if earliest_conversion_time <= earliest_bremsstrahlung_time:
                if earliest_conversion_pdg_code == 111:
                    physics_meso = PhysicsMesoLabel.Pi0Decay
                elif earliest_conversion_pdg_code == 22:
                    physics_meso = PhysicsMesoLabel.PhotonShower
                elif earliest_conversion_pdg_code == 11:
                    physics_meso = PhysicsMesoLabel.ElectronShower
                else:
                    physics_meso = PhysicsMesoLabel.PositronShower
            else:
                if earliest_bremsstrahlung_pdg_code == 11:
                    physics_meso = PhysicsMesoLabel.ElectronShower
                else:
                    physics_meso = PhysicsMesoLabel.PositronShower
        elif (num_bremsstrahlung > 0):
            topology = TopologyLabel.Shower
            if earliest_bremsstrahlung_pdg_code == 11:
                physics_meso = PhysicsMesoLabel.ElectronShower
            else:
                physics_meso = PhysicsMesoLabel.PositronShower
        else:
            if self.simulation_wrangler.trackid_energy_start[particle] > self.shower_threshold:
                topology = TopologyLabel.Shower
                if self.simulation_wrangler.trackid_pdgcode[particle] == 11:
                    physics_meso = PhysicsMesoLabel.ElectronShower
                elif self.simulation_wrangler.trackid_pdgcode[particle] == -11:
                    physics_meso = PhysicsMesoLabel.PositronShower
                else:
                    physics_meso = PhysicsMesoLabel.LowEnergyIonization
            else:
                topology = TopologyLabel.Blip
                physics_meso = PhysicsMesoLabel.LowEnergyIonization

        unique_physics_meso = next(self.unique_physics_meso)

        ionization_hits = self.simulation_wrangler.get_hits_trackid(ionization_descendants)
        ionization_segments = self.simulation_wrangler.get_segments_trackid(ionization_descendants)
        self.simulation_wrangler.set_hit_labels_list(
            ionization_hits,
            ionization_segments,
            ionization_descendants,
            topology,
            PhysicsMicroLabel.ElectronIonization,
            physics_meso,
            unique_topology,
            next(self.unique_physics_micro),
            unique_physics_meso,
        )
        bremsstrahlung_hits = self.simulation_wrangler.get_hits_trackid(bremsstrahlung_descendants)
        bremsstrahlung_segments = self.simulation_wrangler.get_segments_trackid(bremsstrahlung_descendants)
        self.simulation_wrangler.set_hit_labels_list(
            bremsstrahlung_hits,
            bremsstrahlung_segments,
            bremsstrahlung_descendants,
            topology,
            PhysicsMicroLabel.Bremsstrahlung,
            physics_meso,
            unique_topology,
            next(self.unique_physics_micro),
            unique_physics_meso,
        )
        pairprodcharge_hits = self.simulation_wrangler.get_hits_trackid(pairprodcharge_descendants)
        pairprodcharge_segments = self.simulation_wrangler.get_segments_trackid(pairprodcharge_descendants)
        self.simulation_wrangler.set_hit_labels_list(
            pairprodcharge_hits,
            pairprodcharge_segments,
            pairprodcharge_descendants,
            topology,
            PhysicsMicroLabel.ElectronIonization,
            PhysicsMesoLabel.LowEnergyIonization,
            unique_topology,
            next(self.unique_physics_micro),
            unique_physics_meso,
        )
        photoelectric_hits = self.simulation_wrangler.get_hits_trackid(photoelectric_descendants)
        photoelectric_segments = self.simulation_wrangler.get_segments_trackid(photoelectric_descendants)
        self.simulation_wrangler.set_hit_labels_list(
            photoelectric_hits,
            photoelectric_segments,
            photoelectric_descendants,
            TopologyLabel.Blip,
            PhysicsMicroLabel.PhotoElectric,
            physics_meso,
            unique_topology,
            next(self.unique_physics_micro),
            unique_physics_meso,
        )

        compton_hits = self.simulation_wrangler.get_hits_trackid(compton_descendants)
        compton_segments = self.simulation_wrangler.get_segments_trackid(compton_descendants)
        self.simulation_wrangler.set_hit_labels_list(
            compton_hits,
            compton_segments,
            compton_descendants,
            TopologyLabel.Blip,
            PhysicsMicroLabel.GammaCompton,
            physics_meso,
            unique_topology,
            next(self.unique_physics_micro),
            unique_physics_meso,
        )

        annihilation_hits = self.simulation_wrangler.get_hits_trackid(annihilation_descendants)
        annihilation_segments = self.simulation_wrangler.get_segments_trackid(annihilation_descendants)
        self.simulation_wrangler.set_hit_labels_list(
            annihilation_hits,
            annihilation_segments,
            annihilation_descendants,
            topology,
            PhysicsMicroLabel.Annihilation,
            physics_meso,
            unique_topology,
            next(self.unique_physics_micro),
            unique_physics_meso,
        )

        scintillation_hits = self.simulation_wrangler.get_hits_trackid(scintillation_descendants)
        scintillation_segments = self.simulation_wrangler.get_segments_trackid(scintillation_descendants)
        self.simulation_wrangler.set_hit_labels_list(
            scintillation_hits,
            scintillation_segments,
            scintillation_descendants,
            TopologyLabel.Blip,
            PhysicsMicroLabel.ElectronIonization,
            physics_meso,
            unique_topology,
            next(self.unique_physics_micro),
            unique_physics_meso,
        )

        conversion_hits = self.simulation_wrangler.get_hits_trackid(conversion_descendants)
        conversion_segments = self.simulation_wrangler.get_segments_trackid(conversion_descendants)
        self.simulation_wrangler.set_hit_labels_list(
            conversion_hits,
            conversion_segments,
            conversion_descendants,
            TopologyLabel.Shower,
            PhysicsMicroLabel.GammaConversion,
            physics_meso,
            unique_topology,
            next(self.unique_physics_micro),
            unique_physics_meso,
        )

        gamma_descendants = self.simulation_wrangler.filter_trackid_abs_pdg_code(particle_descendants, 22)
        hadron_elastic = self.simulation_wrangler.filter_trackid_subprocess(
            gamma_descendants, SubProcessType.HadronElastic.value
        )
        hadron_inelastic = self.simulation_wrangler.filter_trackid_subprocess(
            gamma_descendants, SubProcessType.HadronInelastic.value
        )
        hadron_at_rest = self.simulation_wrangler.filter_trackid_subprocess(
            gamma_descendants, SubProcessType.HadronCaptureAtRest.value
        )
        if (
            self.simulation_wrangler.trackid_pdgcode[particle] == 22 and
            self.simulation_wrangler.trackid_subprocess[particle] == SubProcessType.HadronElastic.value
        ):
            hadron_elastic.append(particle)
        if (
            self.simulation_wrangler.trackid_pdgcode[particle] == 22 and
            self.simulation_wrangler.trackid_subprocess[particle] == SubProcessType.HadronInelastic.value
        ):
            hadron_inelastic.append(particle)
        if (
            self.simulation_wrangler.trackid_pdgcode[particle] == 22 and
            self.simulation_wrangler.trackid_subprocess[particle] == SubProcessType.HadronCaptureAtRest.value
        ):
            hadron_at_rest.append(particle)

        hadron_elastic_hits = self.simulation_wrangler.get_hits_trackid(hadron_elastic)
        hadron_elastic_segments = self.simulation_wrangler.get_segments_trackid(hadron_elastic)
        self.simulation_wrangler.set_hit_labels_list(
            hadron_elastic_hits,
            hadron_elastic_segments,
            hadron_elastic,
            TopologyLabel.Blip,
            PhysicsMicroLabel.HadronElastic,
            PhysicsMesoLabel.NuclearRecoil,
            next(self.unique_topology),
            next(self.unique_physics_micro),
            next(self.unique_physics_meso),
        )

        hadron_inelastic_hits = self.simulation_wrangler.get_hits_trackid(hadron_inelastic)
        hadron_inelastic_segments = self.simulation_wrangler.get_segments_trackid(hadron_inelastic)
        self.simulation_wrangler.set_hit_labels_list(
            hadron_inelastic_hits,
            hadron_inelastic_segments,
            hadron_inelastic,
            TopologyLabel.Blip,
            PhysicsMicroLabel.HadronInelastic,
            PhysicsMesoLabel.NuclearRecoil,
            next(self.unique_topology),
            next(self.unique_physics_micro),
            next(self.unique_physics_meso),
        )

        hadron_at_rest_hits = self.simulation_wrangler.get_hits_trackid(hadron_at_rest)
        hadron_at_rest_segments = self.simulation_wrangler.get_segments_trackid(hadron_at_rest)
        self.simulation_wrangler.set_hit_labels_list(
            hadron_at_rest_hits,
            hadron_at_rest_segments,
            hadron_at_rest,
            TopologyLabel.Blip,
            PhysicsMicroLabel.GammaCompton,
            PhysicsMesoLabel.LowEnergyIonization,
            next(self.unique_topology),
            next(self.unique_physics_micro),
            next(self.unique_physics_meso),
        )

    def process_showers_list(self, particles, topology_label):
        for particle in particles:
            self.process_showers(particle, topology_label)

    def process_showers_array(self, particles):
        for particle in particles:
            self.process_showers_list(particle, next(self.unique_topology))

    """
    Primaries may have to be handled differently due to the physics_macro scale.
    Perhaps all we need is the GENIE information, and perhaps information about
    the radiological generation.  The physics_macro for neutrino and radiological
    will have to be handled separately.
    """

    def process_electrons(self):
        """
        Primary electrons can come from ... cc_nue? nc? what else?
        For now, only sending to showers.  Sometimes electrons will
        have a process that is equal to zero, which means it's a primary.
        In order to label these correctly, we change the subprocess
        to emionization.
        """
        electrons = self.simulation_wrangler.get_primaries_pdg_code(11)
        for electron in electrons:
            if self.simulation_wrangler.trackid_subprocess[electron] == SubProcessType.Primary.value:
                self.simulation_wrangler.trackid_subprocess[electron] = SubProcessType.Ionization.value
            self.process_showers(electron, next(self.unique_topology))

    def process_positrons(self):
        """
        Primary positrons can come from ... cc_nue? nc? what else?
        For now, only sending to showers.
        """
        positrons = self.simulation_wrangler.get_primaries_pdg_code(-11)
        for positron in positrons:
            if self.simulation_wrangler.trackid_subprocess[positron] == SubProcessType.Primary.value:
                self.simulation_wrangler.trackid_subprocess[positron] = SubProcessType.Ionization.value
            self.process_showers(positron, next(self.unique_topology))

    def process_gammas(self):
        """
        Primary gammas can come from ... cc_nue? nc? what else?
        For now, only sending to showers.
        """
        gammas = self.simulation_wrangler.get_primaries_abs_pdg_code(22)
        for gamma in gammas:
            self.process_showers(gamma, next(self.unique_topology))

    def process_muons(self, debug=False):
        """
        Getting all muons, not just primaries, since they all need to
        be labelled the same.  Therefore, these physics_macro labels
        will need to be set somewhere else.  Or maybe that won't work
        with the unique labels.
        """
        muons = self.simulation_wrangler.get_trackid_abs_pdg_code(
            13
        )  # this is fast // XO
        for i, muon in enumerate(muons):
            # process MIP ionization
            muon_topology = next(self.unique_topology)
            muon_hits = self.simulation_wrangler.trackid_hit[muon]
            muon_segments = self.simulation_wrangler.trackid_segmentid[muon]
            self.simulation_wrangler.set_hit_labels(
                muon_hits,
                muon_segments,
                muon,
                TopologyLabel.Track,
                PhysicsMicroLabel.MIPIonization,
                PhysicsMesoLabel.MIP,
                muon_topology,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )
            # add muon track to tracks object
            self.process_track(muon, muon_topology)
            # process daughters (michel or delta, ignore the rest)
            muon_daughters = self.simulation_wrangler.trackid_daughters[muon]
            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(
                muon_daughters, 11
            )

            # process Michel electrons
            decay_daughters = self.simulation_wrangler.filter_trackid_process(
                elec_daughters, ProcessType.Decay.value
            )
            decay_daughters += self.simulation_wrangler.filter_trackid_subprocess(
                elec_daughters, SubProcessType.HadronCaptureAtRest.value
            )
            michel_decay_hits = self.simulation_wrangler.get_hits_trackid(
                decay_daughters
            )
            michel_decay_segments = self.simulation_wrangler.get_segments_trackid(
                decay_daughters
            )
            for ii, decay in enumerate(decay_daughters):
                decay_topology = next(self.unique_topology)
                self.simulation_wrangler.set_hit_labels(
                    michel_decay_hits[ii],
                    michel_decay_segments[ii],
                    decay_daughters[ii],
                    TopologyLabel.Track,
                    PhysicsMicroLabel.ElectronIonization,
                    PhysicsMesoLabel.MichelElectron,
                    decay_topology,
                    next(self.unique_physics_micro),
                    next(self.unique_physics_meso),
                )
                # add michel track to tracks object
                self.process_track(decay, decay_topology)
            michel_descendants = self.simulation_wrangler.get_descendants_trackid(decay_daughters)
            self.process_showers_array(michel_descendants)

            # process deltas
            elec_em_daughters = self.simulation_wrangler.filter_trackid_process(
                elec_daughters, ProcessType.Electromagnetic.value
            )
            delta_daughters = self.simulation_wrangler.filter_trackid_subprocess(
                elec_em_daughters, SubProcessType.Ionization.value
            )
            delta_hits = self.simulation_wrangler.get_hits_trackid(delta_daughters)
            delta_segments = self.simulation_wrangler.get_segments_trackid(
                delta_daughters
            )
            for ii, delta in enumerate(delta_daughters):
                delta_topology = next(self.unique_topology)
                self.simulation_wrangler.set_hit_labels(
                    delta_hits[ii],
                    delta_segments[ii],
                    delta_daughters[ii],
                    TopologyLabel.Track,
                    PhysicsMicroLabel.ElectronIonization,
                    PhysicsMesoLabel.DeltaElectron,
                    delta_topology,
                    next(self.unique_physics_micro),
                    next(self.unique_physics_meso),
                )
                # add deltas to tracks
                self.process_delta(delta, delta_topology)
            delta_descendants = self.simulation_wrangler.get_descendants_trackid(delta_daughters)
            self.process_showers_array(delta_descendants)
            # process other electrons not michel
            other_daughters = remove_sublist(
                muon_daughters, decay_daughters + delta_daughters
            )
            self.process_showers_list(other_daughters, next(self.unique_topology))

    def process_pion0s(self):
        """
        (From Wikipedia: https://en.wikipedia.org/wiki/Pion)
        The dominant π0 decay mode, with a branching ratio of BR2γ = 0.98823,
        is into two photons: π0 → 2γ.

        The decay π0 → 3γ (as well as decays into any odd number of photons) is
        forbidden by the C-symmetry of the electromagnetic interaction: The intrinsic
        C-parity of the π0 is +1, while the C-parity of a system of n photons is (−1)n.

        The second largest π0 decay mode ( BRγee = 0.01174 ) is the Dalitz decay
        (named after Richard Dalitz), which is a two-photon decay with an internal
        photon conversion resulting a photon and an electron-positron pair in the final state:
            π0 → γ + e− + e+.

        The third largest established decay mode ( BR2e2e = 3.34×10−5 ) is the double-Dalitz
        decay, with both photons undergoing internal conversion which leads to further
        suppression of the rate:
            π0 → e− + e+ + e− +	e+.

        The fourth largest established decay mode is the loop-induced and therefore suppressed
        (and additionally helicity-suppressed) leptonic decay mode ( BRee = 6.46×10−8 ):
            π0 → e− + e+.

        The neutral pion has also been observed to decay into positronium with a branching
        fraction on the order of 10−9. No other decay modes have been established experimentally.
        """
        pi0s = self.simulation_wrangler.get_trackid_abs_pdg_code(111)
        for pi0 in pi0s:
            self.process_showers(pi0, next(self.unique_topology))

    def process_pions(self):
        # pi plus can decay into muons or electrons
        pions = self.simulation_wrangler.get_trackid_abs_pdg_code(211)
        for pion in pions:
            # label pion as HIP ionization
            pion_daughters = self.simulation_wrangler.trackid_daughters[pion]
            pion_hits = self.simulation_wrangler.trackid_hit[pion]
            pion_segments = self.simulation_wrangler.trackid_segmentid[pion]
            cluster_label = next(self.unique_topology)
            self.simulation_wrangler.set_hit_labels(
                pion_hits,
                pion_segments,
                pion,
                TopologyLabel.Track,
                PhysicsMicroLabel.HIPIonization,
                PhysicsMesoLabel.HIP,
                cluster_label,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )
            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(
                pion_daughters, 11
            )
            elec_em_daughters = self.simulation_wrangler.filter_trackid_process(
                elec_daughters, ProcessType.Electromagnetic.value
            )
            delta_daughters = self.simulation_wrangler.filter_trackid_subprocess(
                elec_em_daughters, SubProcessType.Ionization.value
            )
            delta_hits = self.simulation_wrangler.get_hits_trackid(delta_daughters)
            delta_segments = self.simulation_wrangler.get_segments_trackid(
                delta_daughters
            )
            self.simulation_wrangler.set_hit_labels_list(
                delta_hits,
                delta_segments,
                delta_daughters,
                TopologyLabel.Track,
                PhysicsMicroLabel.ElectronIonization,
                PhysicsMesoLabel.DeltaElectron,
                cluster_label,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )
            delta_descendants = self.simulation_wrangler.get_descendants_trackid(delta_daughters)
            self.process_showers_array(delta_descendants)
            other_daughters = remove_sublist(
                pion_daughters, delta_daughters
            )
            self.process_showers_list(other_daughters, next(self.unique_topology))

    def process_kaon0s(self):
        ka0s = self.simulation_wrangler.get_trackid_pdg_code(311)
        for ka0 in ka0s:
            self.process_showers(ka0, next(self.unique_topology))

    def process_kaons(self):
        kaons = self.simulation_wrangler.get_trackid_abs_pdg_code(321)
        for kaon in kaons:
            kaon_daughters = self.simulation_wrangler.trackid_daughters[kaon]
            kaon_hits = self.simulation_wrangler.get_hits_trackid(kaon)
            kaon_segments = self.simulation_wrangler.get_segments_trackid(kaon)

            track_label = next(self.unique_topology)
            self.simulation_wrangler.set_hit_labels(
                kaon_hits[0],
                kaon_segments[0],
                kaon,
                TopologyLabel.Track,
                PhysicsMicroLabel.HIPIonization,
                PhysicsMesoLabel.HIP,
                track_label,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )
            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(
                kaon_daughters, 11
            )
            elec_em_daughters = self.simulation_wrangler.filter_trackid_process(
                elec_daughters, ProcessType.Electromagnetic.value
            )
            delta_daughters = self.simulation_wrangler.filter_trackid_subprocess(
                elec_em_daughters, SubProcessType.Ionization.value
            )
            delta_hits = self.simulation_wrangler.get_hits_trackid(delta_daughters)
            delta_segments = self.simulation_wrangler.get_segments_trackid(
                delta_daughters
            )
            self.simulation_wrangler.set_hit_labels_list(
                delta_hits,
                delta_segments,
                delta_daughters,
                TopologyLabel.Track,
                PhysicsMicroLabel.ElectronIonization,
                PhysicsMesoLabel.DeltaElectron,
                track_label,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )
            delta_descendants = self.simulation_wrangler.get_descendants_trackid(delta_daughters)
            self.process_showers_array(delta_descendants)
            other_daughters = remove_sublist(
                kaon_daughters, delta_daughters
            )
            self.process_showers_list(other_daughters, next(self.unique_topology))

    def process_protons(self):
        protons = self.simulation_wrangler.get_trackid_pdg_code(2212)
        for proton in protons:
            # label protons as HIP ionization
            proton_hits = self.simulation_wrangler.get_hits_trackid(proton)
            proton_segments = self.simulation_wrangler.get_segments_trackid(proton)
            cluster_label = next(self.unique_topology)
            self.simulation_wrangler.set_hit_labels(
                proton_hits[0],
                proton_segments[0],
                proton,
                TopologyLabel.Track,
                PhysicsMicroLabel.HIPIonization,
                PhysicsMesoLabel.HIP,
                cluster_label,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )
            proton_daughters = self.simulation_wrangler.trackid_daughters[proton]
            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(
                proton_daughters, 11
            )
            elec_em_daughters = self.simulation_wrangler.filter_trackid_process(
                elec_daughters, ProcessType.Electromagnetic.value
            )
            delta_daughters = self.simulation_wrangler.filter_trackid_subprocess(
                elec_em_daughters, SubProcessType.Ionization.value
            )
            delta_hits = self.simulation_wrangler.get_hits_trackid(delta_daughters)
            delta_segments = self.simulation_wrangler.get_segments_trackid(
                delta_daughters
            )
            self.simulation_wrangler.set_hit_labels_list(
                delta_hits,
                delta_segments,
                delta_daughters,
                TopologyLabel.Track,
                PhysicsMicroLabel.ElectronIonization,
                PhysicsMesoLabel.DeltaElectron,
                cluster_label,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )
            delta_descendants = self.simulation_wrangler.get_descendants_trackid(delta_daughters)
            self.process_showers_array(delta_descendants)
            other_daughters = remove_sublist(
                proton_daughters, delta_daughters
            )
            self.process_showers_list(other_daughters, next(self.unique_topology))

    def process_neutrons(self):
        """
        Process all different types of neutron interactions, which include
        (1) neutron elastic
        (2) neutron inelastic (which can lead to a proton)
        (3) neutron captures on ar40/ar38/ar36
        (4) fission

        TODO: Only doing ar40 captures right now, need to update this.
        """
        neutrons = self.simulation_wrangler.get_trackid_abs_pdg_code(2112)
        elastic_neutrons = self.simulation_wrangler.filter_trackid_subprocess(
            neutrons, SubProcessType.HadronElastic.value
        )
        inelastic_neutrons = self.simulation_wrangler.filter_trackid_subprocess(
            neutrons, SubProcessType.HadronInelastic.value
        )
        inelastic_neutrons += self.simulation_wrangler.filter_trackid_endsubprocess(
            neutrons, SubProcessType.HadronInelastic.value
        )
        hadron_at_rest_neutrons = self.simulation_wrangler.filter_trackid_subprocess(
            neutrons, SubProcessType.HadronCaptureAtRest.value
        )

        elastic_neutrons_hits = self.simulation_wrangler.get_hits_trackid(elastic_neutrons)
        elastic_neutrons_segments = self.simulation_wrangler.get_segments_trackid(elastic_neutrons)
        self.simulation_wrangler.set_hit_labels_list(
            elastic_neutrons_hits,
            elastic_neutrons_segments,
            elastic_neutrons,
            TopologyLabel.Blip,
            PhysicsMicroLabel.HadronElastic,
            PhysicsMesoLabel.NuclearRecoil,
            next(self.unique_topology),
            next(self.unique_physics_micro),
            next(self.unique_physics_meso),
        )

        inelastic_neutrons_hits = self.simulation_wrangler.get_hits_trackid(inelastic_neutrons)
        inelastic_neutrons_segments = self.simulation_wrangler.get_segments_trackid(inelastic_neutrons)
        self.simulation_wrangler.set_hit_labels_list(
            inelastic_neutrons_hits,
            inelastic_neutrons_segments,
            inelastic_neutrons,
            TopologyLabel.Blip,
            PhysicsMicroLabel.HadronInelastic,
            PhysicsMesoLabel.NuclearRecoil,
            next(self.unique_topology),
            next(self.unique_physics_micro),
            next(self.unique_physics_meso),
        )

        hadron_at_rest_neutrons_hits = self.simulation_wrangler.get_hits_trackid(hadron_at_rest_neutrons)
        hadron_at_rest_neutrons_segments = self.simulation_wrangler.get_segments_trackid(hadron_at_rest_neutrons)
        self.simulation_wrangler.set_hit_labels_list(
            hadron_at_rest_neutrons_hits,
            hadron_at_rest_neutrons_segments,
            hadron_at_rest_neutrons,
            TopologyLabel.Blip,
            PhysicsMicroLabel.HadronInelastic,
            PhysicsMesoLabel.NuclearRecoil,
            next(self.unique_topology),
            next(self.unique_physics_micro),
            next(self.unique_physics_meso),
        )

        other_neutrons = remove_sublist(neutrons, elastic_neutrons + inelastic_neutrons + hadron_at_rest_neutrons)
        other_neutrons_hits = self.simulation_wrangler.get_hits_trackid(other_neutrons)
        other_neutrons_segments = self.simulation_wrangler.get_segments_trackid(other_neutrons)
        self.simulation_wrangler.set_hit_labels_list(
            other_neutrons_hits,
            other_neutrons_segments,
            other_neutrons,
            TopologyLabel.Blip,
            PhysicsMicroLabel.HadronElastic,
            PhysicsMesoLabel.NuclearRecoil,
            next(self.unique_topology),
            next(self.unique_physics_micro),
            next(self.unique_physics_meso),
        )

        for neutron in neutrons:
            # process neutron hits

            neutron_daughters = self.simulation_wrangler.trackid_daughters[neutron]
            gamma_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(
                neutron_daughters, 22
            )
            other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(
                neutron_daughters, 22
            )
            capture_daughters = self.simulation_wrangler.filter_trackid_process(
                gamma_daughters, SubProcessType.HadronCapture.value
            )
            capture_daughters += self.simulation_wrangler.filter_trackid_process(
                gamma_daughters, SubProcessType.HadronCaptureAtRest.value
            )
            other_gammas = remove_sublist(gamma_daughters, capture_daughters)

            for capture in capture_daughters:
                for gamma in capture:
                    gamma_energy = self.simulation_wrangler.get_total_hit_energy(
                        gamma, 5
                    )
                    gamma_hits = self.simulation_wrangler.get_hits_trackid(gamma)
                    gamma_segments = self.simulation_wrangler.get_segments_trackid(
                        gamma
                    )

                    physics_meso_label = PhysicsMesoLabel.NeutronCaptureGammaOther
                    if gamma_energy in {0.00474, 0.00475}:
                        physics_meso_label = PhysicsMesoLabel.NeutronCaptureGamma474
                    elif gamma_energy in {0.00336, 0.00337}:
                        physics_meso_label = PhysicsMesoLabel.NeutronCaptureGamma336
                    elif gamma_energy in {0.00256, 0.00257}:
                        physics_meso_label = PhysicsMesoLabel.NeutronCaptureGamma256
                    elif gamma_energy in {0.00118, 0.00119}:
                        physics_meso_label = PhysicsMesoLabel.NeutronCaptureGamma118
                    elif gamma_energy in {0.00083, 0.00084}:
                        physics_meso_label = PhysicsMesoLabel.NeutronCaptureGamma083
                    elif gamma_energy in {0.00051, 0.00052}:
                        physics_meso_label = PhysicsMesoLabel.NeutronCaptureGamma051
                    elif gamma_energy in {0.00016, 0.00017}:
                        physics_meso_label = PhysicsMesoLabel.NeutronCaptureGamma016

                    cluster_label = next(self.unique_topology)
                    self.simulation_wrangler.set_hit_labels(
                        gamma_hits,
                        gamma_segments,
                        gamma,
                        TopologyLabel.Blip,
                        PhysicsMicroLabel.GammaCompton,
                        physics_meso_label,
                        cluster_label,
                        next(self.unique_physics_micro),
                        next(self.unique_physics_meso),
                    )

            self.process_showers_list(other_daughters, next(self.unique_topology))
            self.process_showers_list(other_gammas, next(self.unique_topology))

    def process_nuclear_recoils(self):
        """
        This section concerns recoiling nuclei, such as argon, sulfur, chlorine,
        and others.  These nuclei are usually scattered from nuclear interactions
        such as neutron capture, or scattering from neutrons/neutrinos.
        """
        sulfur = self.simulation_wrangler.get_trackid_pdg_code(1000160320)
        sulfur += self.simulation_wrangler.get_trackid_pdg_code(1000160330)
        sulfur += self.simulation_wrangler.get_trackid_pdg_code(1000160340)
        sulfur += self.simulation_wrangler.get_trackid_pdg_code(1000160350)
        sulfur += self.simulation_wrangler.get_trackid_pdg_code(1000160360)

        chlorine = self.simulation_wrangler.get_trackid_pdg_code(1000170350)
        chlorine += self.simulation_wrangler.get_trackid_pdg_code(1000170360)
        chlorine += self.simulation_wrangler.get_trackid_pdg_code(1000170370)
        chlorine += self.simulation_wrangler.get_trackid_pdg_code(1000170380)
        chlorine += self.simulation_wrangler.get_trackid_pdg_code(1000170390)
        chlorine += self.simulation_wrangler.get_trackid_pdg_code(1000170400)

        argon = self.simulation_wrangler.get_trackid_pdg_code(1000180360)
        argon += self.simulation_wrangler.get_trackid_pdg_code(1000180370)
        argon += self.simulation_wrangler.get_trackid_pdg_code(1000180380)
        argon += self.simulation_wrangler.get_trackid_pdg_code(1000180390)
        argon += self.simulation_wrangler.get_trackid_pdg_code(1000180400)
        argon += self.simulation_wrangler.get_trackid_pdg_code(1000180410)

        for s in sulfur:
            sulfur_hits = self.simulation_wrangler.get_hits_trackid(s)
            sulfur_segments = self.simulation_wrangler.get_segments_trackid(s)
            self.simulation_wrangler.set_hit_labels(
                sulfur_hits[0],
                sulfur_segments[0],
                s,
                TopologyLabel.Blip,
                PhysicsMicroLabel.HadronElastic,
                PhysicsMesoLabel.NuclearRecoil,
                next(self.unique_topology),
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )

        for cl in chlorine:
            chlorine_hits = self.simulation_wrangler.get_hits_trackid(cl)
            chlorine_segments = self.simulation_wrangler.get_segments_trackid(cl)
            self.simulation_wrangler.set_hit_labels(
                chlorine_hits[0],
                chlorine_segments[0],
                cl,
                TopologyLabel.Blip,
                PhysicsMicroLabel.HadronElastic,
                PhysicsMesoLabel.NuclearRecoil,
                next(self.unique_topology),
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )

        for ar in argon:
            argon_hits = self.simulation_wrangler.get_hits_trackid(ar)
            argon_segments = self.simulation_wrangler.get_segments_trackid(ar)
            self.simulation_wrangler.set_hit_labels(
                argon_hits[0],
                argon_segments[0],
                ar,
                TopologyLabel.Blip,
                PhysicsMicroLabel.HadronElastic,
                PhysicsMesoLabel.NuclearRecoil,
                next(self.unique_topology),
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )

    def process_baryons(self):
        """
        This section concerns light and strange baryons, such as the Delta,
        Lambda and Sigma varieties.
        """
        delta_baryons = self.simulation_wrangler.get_trackid_pdg_code(2224)
        delta_baryons += self.simulation_wrangler.get_trackid_pdg_code(2214)
        delta_baryons += self.simulation_wrangler.get_trackid_pdg_code(2114)
        delta_baryons += self.simulation_wrangler.get_trackid_pdg_code(1114)

        lambda_baryons = self.simulation_wrangler.get_trackid_pdg_code(3122)

        sigma_baryons = self.simulation_wrangler.get_trackid_pdg_code(3222)
        sigma_baryons += self.simulation_wrangler.get_trackid_pdg_code(3212)
        sigma_baryons += self.simulation_wrangler.get_trackid_pdg_code(3112)

        for delta in delta_baryons:
            delta_hits = self.simulation_wrangler.get_hits_trackid(delta)
            delta_segments = self.simulation_wrangler.get_segments_trackid(delta)
            self.simulation_wrangler.set_hit_labels(
                delta_hits[0],
                delta_segments[0],
                delta,
                TopologyLabel.Blip,
                PhysicsMicroLabel.HadronElastic,
                PhysicsMesoLabel.NuclearRecoil,
                next(self.unique_topology),
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )

        for lamb in lambda_baryons:
            lamb_hits = self.simulation_wrangler.get_hits_trackid(lamb)
            lamb_segments = self.simulation_wrangler.get_segments_trackid(lamb)
            self.simulation_wrangler.set_hit_labels(
                lamb_hits[0],
                lamb_segments[0],
                lamb,
                TopologyLabel.Blip,
                PhysicsMicroLabel.HadronElastic,
                PhysicsMesoLabel.NuclearRecoil,
                next(self.unique_topology),
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )

        for sigma in sigma_baryons:
            sigma_hits = self.simulation_wrangler.get_hits_trackid(sigma)
            sigma_segments = self.simulation_wrangler.get_segments_trackid(sigma)
            self.simulation_wrangler.set_hit_labels(
                sigma_hits[0],
                sigma_segments[0],
                sigma,
                TopologyLabel.Blip,
                PhysicsMicroLabel.HadronElastic,
                PhysicsMesoLabel.NuclearRecoil,
                next(self.unique_topology),
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )

    def process_electron_recoils(self):
        deuterons = self.simulation_wrangler.get_trackid_pdg_code(1000010020)
        tritons = self.simulation_wrangler.get_trackid_pdg_code(1000010030)
        alphas = self.simulation_wrangler.get_trackid_pdg_code(1000020040)

        for deuteron in deuterons:
            deuteron_hits = self.simulation_wrangler.get_hits_trackid(deuteron)
            deuteron_segments = self.simulation_wrangler.get_segments_trackid(deuteron)

            self.simulation_wrangler.set_hit_labels(
                deuteron_hits[0],
                deuteron_segments[0],
                deuteron,
                TopologyLabel.Blip,
                PhysicsMicroLabel.HadronElastic,
                PhysicsMesoLabel.ElectronRecoil,
                next(self.unique_topology),
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )

        for triton in tritons:
            triton_hits = self.simulation_wrangler.get_hits_trackid(triton)
            triton_segments = self.simulation_wrangler.get_segments_trackid(triton)

            self.simulation_wrangler.set_hit_labels(
                triton_hits[0],
                triton_segments[0],
                triton,
                TopologyLabel.Blip,
                PhysicsMicroLabel.HadronElastic,
                PhysicsMesoLabel.ElectronRecoil,
                next(self.unique_topology),
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )

        for alpha in alphas:
            alpha_hits = self.simulation_wrangler.get_hits_trackid(alpha)
            alpha_segments = self.simulation_wrangler.get_segments_trackid(alpha)

            self.simulation_wrangler.set_hit_labels(
                alpha_hits[0],
                alpha_segments[0],
                alpha,
                TopologyLabel.Blip,
                PhysicsMicroLabel.HadronElastic,
                PhysicsMesoLabel.ElectronRecoil,
                next(self.unique_topology),
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
            )

    def process_ar39(self):
        """
        Argon-39 decays via beta decay into Potassium-39,
        with a Q-value of 565 keV: http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=180039.
        """
        pass

        # mc_data = SimulationWrangler()
        # ar39 = mc_data.get_primaries_generator_label(PhysicsLabel.Ar39)
        # # ar39_daughters = mc_data.get_daughters_trackid(ar39)

        # if ar39 is not None:
        #     for ar in ar39:
        #         ar39_hits = mc_data.get_hits_trackid(ar)
        #         ar39_segments = mc_data.get_segments_trackid(ar)

        #         self.simulation_wrangler.set_hit_labels(
        #             ar39_hits, ar39_segments, ar,
        #             TopologyLabel.Blip, PhysicsLabel.Ar39,
        #             next(self.unique_topology)
        #         )

    def process_ar42(self):
        """
        Argon-42 decays via beta decay into Potassium-42,
        with a Q-value of 599 keV: http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=180042.
        There is a subtlety in the way this decay is simulated. The lifetime of K42 is approximately
        12 hours, which beta decays to Calcium-42 with a 3.5 MeV beta. Calcium-42 is stable.
        See here for some details:
        https://indico.fnal.gov/event/50121/contributions/220205/attachments/145404/185102/20210721_Decay0_Lasorak.pdf

        If the energy of the decay primary is 600 keV or less, we label it as an
        Ar42 decay, otherwise, it's a K42 decay.
        """
        pass
        # mc_data = SimulationWrangler()
        # ar42 = mc_data.get_primaries_generator_label(PhysicsLabel.Ar42)
        # #ar42_betas = mc_data.filter_trackid_pdgcode(ar42, 1000180420)
        # # ar42_daughters = mc_data.get_daughters_trackid(ar42_betas)
        # if ar42 is not None:
        #     for ar in ar42:
        #         ar42_hits = mc_data.get_hits_trackid(ar)
        #         ar42_segments = mc_data.get_segments_trackid(ar)

        #         if mc_data.get_energy_trackid(ar) < 0.600:
        #             self.simulation_wrangler.set_hit_labels(
        #                 ar42_hits, ar42_segments, ar,
        #                 TopologyLabel.Blip, PhysicsLabel.Ar42,
        #                 next(self.unique_topology)
        #             )
        #         else:
        #             self.simulation_wrangler.set_hit_labels(
        #                 ar42_hits, ar42_segments, ar,
        #                 TopologyLabel.Blip, PhysicsLabel.K42,
        #                 next(self.unique_topology)
        #             )

    def process_kr85(self):
        """
        Krypton-85 decays via beta decay into Rubidium 85
        with two prominent betas with energies of 687 keV (99.56 %) and
        173 keV (.43 %): http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=360085.
        """
        pass
        # mc_data = SimulationWrangler()
        # kr85 = mc_data.get_primaries_generator_label(PhysicsLabel.Kr85)
        # # kr85_daughters = mc_data.get_daughters_trackid(kr85)

        # if kr85 is not None:
        #     for kr in kr85:
        #         kr85_hits= mc_data.get_hits_trackid(kr)
        #         kr85_segments = mc_data.get_segments_trackid(kr)

        #         self.simulation_wrangler.set_hit_labels(
        #             kr85_hits, kr85_segments, kr,
        #             TopologyLabel.Blip, PhysicsLabel.Kr85,
        #             next(self.unique_topology)
        # )

    def process_rn222(self):
        """
        Radon-222 decays via alpha decay through a chain that ends in lead
        (https://en.wikipedia.org/wiki/Radon-222).  The alpha has an energy of
        5.5904 MeV, which bounces around locally in Argon, but quickly thermalizes
        due to the short scattering length.  The CSDA range of a 5.5 MeV alpha in
        Argon is about 7.5e-3 g/cm^2.  Using a density of 1.3954 g/cm^3, the
        scattering length is (~0.005 cm) or (~50 um).

        The simulation of radiologicals is done through a wrap of Decay0, which does
        not save information about what decay each simulated particle came from.
        Instead, it saves each decay product (beta, gamma, alpha) as a primary of its
        own.  This can be seen from the following lines of the larsimrad file
        https://github.com/LArSoft/larsimrad/blob/develop/larsimrad/BxDecay0/Decay0Gen_module.cc,
        which has (starting at line 162):
        ...

        We will therefore need some way to figure out which particle came from what decay.
        The possible decays are:
        Rn222 -> Po218 - alpha ~ 5.590 MeV
        Po218 -> Pb214 - alpha ~ 6.115 MeV
        Po218 -> At218 - beta  ~ 0.294 MeV
        At218 -> Bi214 - alpha ~ 6.874 MeV
        At218 -> Rn218 - beta  ~ 2.883 MeV
        Rn218 -> Po214 - alpha ~ 7.263 MeV
        Pb214 -> Bi214 - beta  ~ 1.024 MeV
        Bi214 -> Tl210 - alpha ~ 5.617 MeV
        Bi214 -> Po214 - beta  ~ 3.272 MeV
        Po214 -> Pb210 - alpha ~ 7.833 MeV
        Tl210 -> Pb210 - beta  ~ 5.484 MeV
        Pb210 -> Bi210 - beta  ~ 0.064 MeV
        Pb210 -> Hg206 - alpha ~ 3.792 MeV
        Bi210 -> Po210 - beta  ~ 1.163 MeV
        Bi210 -> Tl206 - alpha ~ 5.037 MeV
        Po210 -> Pb206 - alpha ~ 5.407 MeV

        We therefore have 9 distinct alpha energies
        (7.833, 7.263, 6.874, 6.115, 5.617, 5.590, 5.407, 5.037, 3.792)

        and 7 distinct beta energies
        (5.484, 3.272, 2.883, 1.163, 1.024, 0.294, 0.064)

        We can assign labels then based on whatever energy is closest to each primary.
        """
        pass
        # mRn222Decays = (
        #     PhysicsLabel.Rn222,
        #     PhysicsLabel.Po218a, PhysicsLabel.Po218b,
        #     PhysicsLabel.At218a, PhysicsLabel.At218b,
        #     PhysicsLabel.Rn218,
        #     PhysicsLabel.Pb214,
        #     PhysicsLabel.Bi214a, PhysicsLabel.Bi214b,
        #     PhysicsLabel.Po214,
        #     PhysicsLabel.Tl210,
        #     PhysicsLabel.Pb210a, PhysicsLabel.Pb210b,
        #     PhysicsLabel.Bi210a, PhysicsLabel.Bi210b,
        #     PhysicsLabel.Po210
        # )
        # mRn222PDGs = (
        #     1000020040,
        #     1000020040, 11,
        #     1000020040, 11,
        #     1000020040,
        #     11,
        #     1000020040, 11,
        #     1000020040,
        #     11,
        #     1000020040, 11,
        #     1000020040, 11,
        #     1000020040
        # )
        # mRn222Energies = (
        #     5.590,
        #     6.115, 0.294,
        #     6.874, 2.883,
        #     7.263,
        #     1.024,
        #     5.627, 3.272,
        #     7.833,
        #     5.484,
        #     3.792, 0.064,
        #     5.037, 1.163,
        #     5.407
        # )

        # mc_data = SimulationWrangler()
        # rn222 = mc_data.get_primaries_generator_label(PhysicsLabel.Rn222)
        # if rn222 is not None:
        #     for rn in rn222:
        #         index = 0
        #         energy_diff = 10e10
        #         for ii in range(len(mRn222Energies)):
        #             if mc_data.get_pdg_code_trackid(rn) != mRn222PDGs[ii]:
        #                 continue
        #             if abs(mc_data.get_energy_trackid(rn) - mRn222Energies[ii]) < energy_diff:
        #                 index = ii
        #                 energy_diff = mc_data.get_energy_trackid(rn) - mRn222Energies[ii]

        #         rn222_hits = mc_data.get_hits_trackid(rn)
        #         rn222_segments = mc_data.get_all_segments_trackid(rn)

        #         self.simulation_wrangler.set_hit_labels(
        #             rn222_hits, rn222_segments, rn,
        #             TopologyLabel.Blip, mRn222Decays[index],
        #             next(self.unique_topology)
        #         )

    def process_cosmics(self):
        """ """
        pass
        # mc_data = SimulationWrangler()
        # cosmics = mc_data.get_primaries_generator_label(ParticleLabel.Cosmics) # there's no cosmics label?
        # electrons = mc_data.filter_trackid_pdg_code(cosmics, 11)
        # positrons = mc_data.filter_trackid_pdg_code(cosmics, -11)
        # # gammas = mc_data.filter_trackid_abs_pdg_code(cosmics, 22)
        # # neutrons = mc_data.filter_trackid_pdg_code(cosmics, 2112)
        # # anti_neutrons = mc_data.filter_trackid_pdg_code(cosmics, -2112)
        # if electrons is not None:
        #     for electron in electrons:
        #         electron_daughters = mc_data.get_daughters_trackid(electron)
        #         elec_daughters = mc_data.filter_trackid_abs_pdg_code(electron_daughters, 11)
        #         other_daughters = mc_data.filter_trackid_not_abs_pdg_code(electron_daughters, 11)
        #         elec_hits = mc_data.get_segments_trackid(elec_daughters)
        #         elec_segments = mc_data.get_hits_trackid(elec_daughters)

        #         electron_progeny = mc_data.get_progeny_trackid(electron)
        #         electron_hits = mc_data.get_hits_trackid(electron)
        #         electron_segments= mc_data.get_segments_trackid(electron)

        #         # Set electron detsim labels to Shower::ElectronShower
        #         shower_label = next(self.unique_topology)
        #         self.simulation_wrangler.set_hit_labels(
        #             electron_hits, electron_segments, electron,
        #             TopologyLabel.Shower, PhysicsLabel.ElectronRecoil,
        #             shower_label
        #         )

        #         self.simulation_wrangler.set_hit_labels(
        #             elec_hits, elec_segments, elec_daughters,
        #             TopologyLabel.Shower, PhysicsLabel.ElectronRecoil,
        #             shower_label
        #         )

        #         self.process_showers(electron_progeny, shower_label)
        #         self.process_showers(other_daughters, next(self.unique_topology))
        # if positrons is not None:
        #     for positron in positrons:
        #         positron_daughters = mc_data.get_daughters_trackid(positron)
        #         elec_daughters = mc_data.filter_trackid_abs_pdg_code(positron_daughters, 11)
        #         other_daughters = mc_data.filter_trackid_not_abs_pdg_code(positron_daughters, 11)
        #         elec_hits = mc_data.get_hits_trackid(elec_daughters)
        #         elec_segments = mc_data.get_segments_trackid(elec_daughters)

        #         positron_progeny = mc_data.get_progeny_trackid(positron)
        #         positron_hits = mc_data.get_hits_trackid(positron)
        #         positron_segments = mc_data.get_segments_trackid(positron)

        #         # Set positron detsim labels to Shower::PositronShower
        #         shower_label = next(self.unique_topology)
        #         self.simulation_wrangler.set_hit_labels(
        #             positron_hits, positron_segments, positron,
        #             TopologyLabel.Shower, PhysicsLabel.PositronShower,
        #             shower_label
        #         )

        #         self.simulation_wrangler.set_labels_list(
        #             elec_hits, elec_segments, elec_daughters,
        #             TopologyLabel.Shower, PhysicsLabel.PositronShower,
        #             shower_label
        #         )

        #         self.process_showers(positron_progeny, shower_label)
        #         self.process_showers(other_daughters, next(self.unique_topology))

    def process_primaries(self):
        pass
