"""
Simluation Labeling Logic

Developers: Nicholas Carrara        [nmcarrara@ucdavis.edu]
            Marjolein von Nuland    [mnuland@nikhef.nl]

ChangeLog:  12/17/2023 - started putting together shower logic.
"""
import numpy as np

from arrakis_nd.utils.logger import Logger
from arrakis_nd.dataset.common import TopologyLabel, PhysicsMicroLabel, PhysicsMesoLabel, PhysicsMacroLabel
from arrakis_nd.utils.timing import Timers
from arrakis_nd.utils.memory import MemoryTrackers
from arrakis_nd.utils.utils import ResetableIterator


class SimulationLabelingLogic:
    """
    Simulation labeling logic conducts mid-level and high-level
    variable construction from low-level information, such as
    pdg code, particle heirarchy, process and subprocess, etc.
    
    
    """
    def __init__(
        self,
        name:   str = '',
        config: dict = {},
        meta:   dict = {},
    ):
        self.name = name + '_simulation_labeling_logic'
        self.config = config
        self.meta = meta

        if "device" in self.meta:
            self.device = self.meta['device']
            self.gpu = self.meta['gpu']
        else:
            self.device = 'cpu'
            self.gpu = False
        if meta['verbose']:
            self.logger = Logger(self.name, output="both",   file_mode="w")
        else:
            self.logger = Logger(self.name, level='warning', file_mode="w")

        self.unique_topology = ResetableIterator()
        self.unique_physics_micro = ResetableIterator()
        self.unique_physics_meso = ResetableIterator()
        self.unique_physics_macro = ResetableIterator()

        self.parse_config()

    def parse_config(self):
        self.check_config()
        self.parse_timers()
        self.parse_memory_trackers()

    def check_config(self):
        if 'simulation_wrangler' not in self.meta.keys():
            self.logger.error('no simulation_wrangler defined in meta!')
        self.simulation_wrangler = self.meta['simulation_wrangler']
        if 'shower_threshold' not in self.config.keys():
            self.logger.warn('no shower_threshold in config! setting to "20"')
            self.config["shower_threshold"] = 20
        self.shower_threshold = self.config["shower_threshold"]
        if "debug" not in self.config.keys():
            self.logger.warn('debug not specified in config! setting to "False"')
            self.config['debug'] = False
        self.debug = self.config['debug']

    def parse_timers(self):
        self.timers = Timers(gpu=self.gpu)

    def parse_memory_trackers(self):
        self.memory_trackers = MemoryTrackers(gpu=self.gpu)

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

    def _process_event_without_timing(self):
        self.unique_topology.reset
        self.unique_physics_micro.reset
        self.unique_physics_meso.reset
        self.unique_physics_macro.reset
        self.process_electrons()
        self.process_positrons()
        self.process_gammas()
        self.process_muons()
        self.process_anti_muons()
        self.process_pion0s()
        self.process_pion_plus()
        self.process_pion_minus()
        self.process_kaon0s()
        self.process_kaon_plus()
        self.process_kaon_minus()
        self.process_protons()
        self.process_neutrons()
        self.process_nuclear_recoils()
        self.process_electron_recoils()
        self.process_ar39()
        self.process_ar42()
        self.process_kr85()
        self.process_rn222()
        self.process_cosmics()
        self.process_primaries()

    def _process_event_with_timing(self):
        self.unique_topology.reset
        self.unique_physics_micro.reset
        self.unique_physics_meso.reset
        self.unique_physics_macro.reset

        self.timers.timers['process_all'].start()
        self.memory_trackers.memory_trackers['process_all'].start()

        self.timers.timers['process_electrons'].start()
        self.memory_trackers.memory_trackers['process_electrons'].start()
        self.process_electrons()
        self.timers.timers['process_electrons'].end()
        self.memory_trackers.memory_trackers['process_electrons'].end()

        self.timers.timers['process_positrons'].start()
        self.memory_trackers.memory_trackers['process_positrons'].start()
        self.process_positrons()
        self.timers.timers['process_positrons'].end()
        self.memory_trackers.memory_trackers['process_positrons'].end()

        self.timers.timers['process_gammas'].start()
        self.memory_trackers.memory_trackers['process_gammas'].start()
        self.process_gammas()
        self.timers.timers['process_gammas'].end()
        self.memory_trackers.memory_trackers['process_gammas'].end()

        self.timers.timers['process_muons'].start()
        self.memory_trackers.memory_trackers['process_muons'].start()
        self.process_muons()
        self.timers.timers['process_muons'].end()
        self.memory_trackers.memory_trackers['process_muons'].end()

        self.timers.timers['process_anti_muons'].start()
        self.memory_trackers.memory_trackers['process_anti_muons'].start()
        self.process_anti_muons()
        self.timers.timers['process_anti_muons'].end()
        self.memory_trackers.memory_trackers['process_anti_muons'].end()

        self.timers.timers['process_pion0s'].start()
        self.memory_trackers.memory_trackers['process_pion0s'].start()
        self.process_pion0s()
        self.timers.timers['process_pion0s'].end()
        self.memory_trackers.memory_trackers['process_pion0s'].end()

        self.timers.timers['process_pion_plus'].start()
        self.memory_trackers.memory_trackers['process_pion_plus'].start()
        self.process_pion_plus()
        self.timers.timers['process_pion_plus'].end()
        self.memory_trackers.memory_trackers['process_pion_plus'].end()

        self.timers.timers['process_pion_minus'].start()
        self.memory_trackers.memory_trackers['process_pion_minus'].start()
        self.process_pion_minus()
        self.timers.timers['process_pion_minus'].end()
        self.memory_trackers.memory_trackers['process_pion_minus'].end()

        self.timers.timers['process_kaon0s'].start()
        self.memory_trackers.memory_trackers['process_kaon0s'].start()
        self.process_kaon0s()
        self.timers.timers['process_kaon0s'].end()
        self.memory_trackers.memory_trackers['process_kaon0s'].end()

        self.timers.timers['process_kaon_plus'].start()
        self.memory_trackers.memory_trackers['process_kaon_plus'].start()
        self.process_kaon_plus()
        self.timers.timers['process_kaon_plus'].end()
        self.memory_trackers.memory_trackers['process_kaon_plus'].end()

        self.timers.timers['process_kaon_minus'].start()
        self.memory_trackers.memory_trackers['process_kaon_minus'].start()
        self.process_kaon_minus()
        self.timers.timers['process_kaon_minus'].end()
        self.memory_trackers.memory_trackers['process_kaon_minus'].end()

        self.timers.timers['process_protons'].start()
        self.memory_trackers.memory_trackers['process_protons'].start()
        self.process_protons()
        self.timers.timers['process_protons'].end()
        self.memory_trackers.memory_trackers['process_protons'].end()

        self.timers.timers['process_neutrons'].start()
        self.memory_trackers.memory_trackers['process_neutrons'].start()
        self.process_neutrons()
        self.timers.timers['process_neutrons'].end()
        self.memory_trackers.memory_trackers['process_neutrons'].end()

        self.timers.timers['process_nuclear_recoils'].start()
        self.memory_trackers.memory_trackers['process_nuclear_recoils'].start()
        self.process_nuclear_recoils()
        self.timers.timers['process_nuclear_recoils'].end()
        self.memory_trackers.memory_trackers['process_nuclear_recoils'].end()

        self.timers.timers['process_electron_recoils'].start()
        self.memory_trackers.memory_trackers['process_electron_recoils'].start()
        self.process_electron_recoils()
        self.timers.timers['process_electron_recoils'].end()
        self.memory_trackers.memory_trackers['process_electron_recoils'].end()

        self.timers.timers['process_ar39'].start()
        self.memory_trackers.memory_trackers['process_ar39'].start()
        self.process_ar39()
        self.timers.timers['process_ar39'].end()
        self.memory_trackers.memory_trackers['process_ar39'].end()

        self.timers.timers['process_ar42'].start()
        self.memory_trackers.memory_trackers['process_ar42'].start()
        self.process_ar42()
        self.timers.timers['process_ar42'].end()
        self.memory_trackers.memory_trackers['process_ar42'].end()

        self.timers.timers['process_kr85'].start()
        self.memory_trackers.memory_trackers['process_kr85'].start()
        self.process_kr85()
        self.timers.timers['process_kr85'].end()
        self.memory_trackers.memory_trackers['process_kr85'].end()

        self.timers.timers['process_rn222'].start()
        self.memory_trackers.memory_trackers['process_rn222'].start()
        self.process_rn222()
        self.timers.timers['process_rn222'].end()
        self.memory_trackers.memory_trackers['process_rn222'].end()

        self.timers.timers['process_cosmics'].start()
        self.memory_trackers.memory_trackers['process_cosmics'].start()
        self.process_cosmics()
        self.timers.timers['process_cosmics'].end()
        self.memory_trackers.memory_trackers['process_cosmics'].end()

        self.timers.timers['process_all'].end()
        self.memory_trackers.memory_trackers['process_all'].end()

        self.timers.timers['check_labels'].start()
        self.memory_trackers.memory_trackers['check_labels'].start()
        self.check_labels()
        self.timers.timers['check_labels'].end()
        self.memory_trackers.memory_trackers['check_labels'].end()

    def check_labels(self):
        for particle in self.simulation_wrangler.trackid_hit.keys():
            hits = self.simulation_wrangler.trackid_hit[particle]
            topology_labels = [self.simulation_wrangler.det_point_cloud.data["topology_label"][hit] for hit in hits]
            # xs = [self.simulation_wrangler.det_point_cloud.data["x"][hit] for hit in hits]
            # ys = [self.simulation_wrangler.det_point_cloud.data["y"][hit] for hit in hits]
            # zs = [self.simulation_wrangler.det_point_cloud.data["z"][hit] for hit in hits]
            if -1 in topology_labels:
                self.simulation_wrangler.print_particle_data(particle)
                # pdg = self.simulation_wrangler.trackid_pdgcode[particle]
                # parentid = self.simulation_wrangler.trackid_parentid[particle]
                # parent_pdg = self.simulation_wrangler.trackid_pdgcode[parentid]
                # ancestry = self.simulation_wrangler.trackid_ancestry[particle]
                # ancestry_pdgs = [self.simulation_wrangler.trackid_pdgcode[ancestor] for ancestor in ancestry]
                # self.logger.warn('#############################################')
                # self.logger.warn(f'# Undefined points for trackid: {particle} #')
                # self.logger.warn(f'# pdg:       {pdg}')
                # self.logger.warn(f'# parent:    {parentid}')
                # self.logger.warn(f'# parent_pdg:    {parent_pdg}')
                # self.logger.warn(f'# ancestry:  {ancestry}')
                # self.logger.warn(f'# ancestry_pdg:  {ancestry_pdgs}')
                
            #         f"{particle}: undefined points for particle: {self.simulation_wrangler.trackid_pdgcode[particle]} with process: {self.simulation_wrangler.trackid_process[particle]}:{self.simulation_wrangler.trackid_subprocess[particle]}")
            #     for ii, label in enumerate(topology_labels):
            #         if label == -1:
            #             print(f"(x,y,z): ({xs[ii]},{ys[ii]},{zs[ii]}])")
            # descendants = self.simulation_wrangler.trackid_descendants[particle]
            # for descendant in descendants:
            #     hits = self.simulation_wrangler.trackid_hit[descendant]
            #     topology_labels = [self.simulation_wrangler.det_point_cloud.data["topology_label"][hit] for hit in hits]
            #     xs = [self.simulation_wrangler.det_point_cloud.data["x"][hit] for hit in hits]
            #     ys = [self.simulation_wrangler.det_point_cloud.data["y"][hit] for hit in hits]
            #     zs = [self.simulation_wrangler.det_point_cloud.data["z"][hit] for hit in hits]
            #     if -1 in topology_labels:
            #         print(f"{particle}: ancestor: {self.simulation_wrangler.trackid_pdgcode[particle]} with process: {self.simulation_wrangler.trackid_process[particle]}:{self.simulation_wrangler.trackid_subprocess[particle]}")
            #         print(f"{descendant}:{self.simulation_wrangler.trackid_ancestry[descendant]}:{[self.simulation_wrangler.trackid_pdgcode[anc] for anc in self.simulation_wrangler.trackid_ancestry[descendant]]}: undefined points for particle: {self.simulation_wrangler.trackid_pdgcode[descendant]} with process: {self.simulation_wrangler.trackid_process[descendant]}:{self.simulation_wrangler.trackid_subprocess[descendant]}")
            #         for ii, label in enumerate(topology_labels):
            #             if label == -1:
            #                 print(f"(x,y,z): ({xs[ii]},{ys[ii]},{zs[ii]}])")

    def process_showers(
        self,
        particle:           int = 0,
        unique_topology:    int = 0
    ):
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
        # First determine if this is a shower.
        shower = True
        shower_energy_threshold = 20.0

        # shower or no shower?
        if shower or shower_energy_threshold > 0.0:
            topology = TopologyLabel.Shower
        else:
            topology = TopologyLabel.Blip

        # unique_physics_meso is one number, and
        #   physics_meso is one of the three for everyone.

        # if no shower, then physics_meso is the particle type for each segment.
        #   if the total energy is low enough, the topology type is changed to Blip.

        if self.simulation_wrangler.trackid_pdgcode[particle] == 11:
            physics_meso = PhysicsMesoLabel.ElectronShower
        elif self.simulation_wrangler.trackid_pdgcode[particle] == -11:
            physics_meso = PhysicsMesoLabel.PositronShower
        elif self.simulation_wrangler.trackid_pdgcode[particle] == 22:
            physics_meso = PhysicsMesoLabel.PhotonShower
        else:
            physics_meso = PhysicsMesoLabel.Undefined    # TODO: fix this in the future of course
        unique_physics_meso = next(self.unique_physics_meso)

        particle_descendants = self.simulation_wrangler.trackid_descendants[particle]
        compton_descendants = self.simulation_wrangler.filter_trackid_subprocess(particle_descendants, 13)
        bremmstrahlung_descendants = self.simulation_wrangler.filter_trackid_subprocess(particle_descendants, 3)
        photoelectric_descendants = self.simulation_wrangler.filter_trackid_subprocess(particle_descendants, 12)
        gamma_descendants = self.simulation_wrangler.filter_trackid_abs_pdg_code(particle_descendants, 22)
        conversion_descendants = self.simulation_wrangler.filter_trackid_subprocess(particle_descendants, 14)

        # process gamma conversions
        for conversion in conversion_descendants:
            conversion_daughters = self.simulation_wrangler.trackid_daughters[conversion]
            conversion_hits = self.simulation_wrangler.get_hits_trackid(conversion)
            conversion_segments = self.simulation_wrangler.get_segments_trackid(conversion)
            daughter_hits = self.simulation_wrangler.get_hits_trackid(conversion_daughters)
            daughter_segments = self.simulation_wrangler.get_segments_trackid(conversion_daughters)
            unique_physics_micro = next(self.unique_physics_micro)
            self.simulation_wrangler.set_hit_labels(
                conversion_hits,
                conversion_segments,
                conversion,
                topology,
                PhysicsMicroLabel.GammaConversion,
                physics_meso,
                PhysicsMacroLabel.Undefined,
                unique_topology,
                unique_physics_micro,
                unique_physics_meso,
                0,
            )
            # self.simulation_wrangler.set_hit_labels_list(
            #     daughter_hits,
            #     daughter_segments,
            #     conversion_daughters,
            #     topology,
            #     PhysicsMicroLabel.GammaConversion,
            #     physics_meso,
            #     PhysicsMacroLabel.Undefined,
            #     unique_topology,
            #     next(self.unique_physics_micro),
            #     unique_physics_meso,
            #     0,
            # )

        # # process compton scatters
        # compton_parents = np.array(self.simulation_wrangler.get_parentid_trackid(compton_descendants))
        # compton_descendants = np.array(compton_descendants)
        # unique_compton_parents = np.unique(compton_parents)
        # for unique_compton_parent in unique_compton_parents:
        #     unique_compton_descendants = compton_descendants[(compton_parents == unique_compton_parent)]
        #     # these get labeled with the same unique_physics_micro label
        #     unique_compton_descendants_hits = self.simulation_wrangler.get_hits_trackid(unique_compton_descendants)
        #     unique_compton_descendants_segments = self.simulation_wrangler.get_segments_trackid(unique_compton_descendants)
        #     self.simulation_wrangler.set_hit_labels_list(
        #         unique_compton_descendants_hits,
        #         unique_compton_descendants_segments,
        #         unique_compton_descendants,
        #         topology,
        #         PhysicsMicroLabel.GammaCompton,
        #         physics_meso,
        #         PhysicsMacroLabel.Undefined,
        #         unique_topology,
        #         next(self.unique_physics_micro),
        #         unique_physics_meso,
        #         0
        #     )
        # # process photo-electric effect
        # photoelectric_parents = np.array(self.simulation_wrangler.get_parentid_trackid(photoelectric_descendants))
        # photoelectric_descendants = np.array(photoelectric_descendants)
        # unique_photoelectric_parents = np.unique(photoelectric_parents)
        # for unique_photoelectric_parent in unique_photoelectric_parents:
        #     unique_photoelectric_descendants = photoelectric_descendants[(photoelectric_parents == unique_photoelectric_parent)]
        #     # these get labeled with the same unique_physics_micro label
        #     unique_photoelectric_descendants_hits = self.simulation_wrangler.get_hits_trackid(unique_photoelectric_descendants)
        #     unique_photoelectric_descendants_segments = self.simulation_wrangler.get_segments_trackid(unique_photoelectric_descendants)
        #     self.simulation_wrangler.set_hit_labels_list(
        #         unique_photoelectric_descendants_hits,
        #         unique_photoelectric_descendants_segments,
        #         unique_photoelectric_descendants,
        #         topology,
        #         PhysicsMicroLabel.ElectronIonization,
        #         physics_meso,
        #         PhysicsMacroLabel.Undefined,
        #         unique_topology,
        #         next(self.unique_physics_micro),
        #         unique_physics_meso,
        #         0
        #     )

    def process_showers_list(
        self,
        particles,
        topology_label
    ):
        for particle in particles:
            self.process_showers(particle, topology_label)

    def process_showers_array(
        self,
        particles
    ):
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
        For now, only sending to showers.
        """
        electrons = self.simulation_wrangler.get_primaries_pdg_code(11)
        for electron in electrons:
            self.process_showers(electron, next(self.unique_topology))

    def process_positrons(self):
        """
        Primary positrons can come from ... cc_nue? nc? what else?
        For now, only sending to showers.
        """
        positrons = self.simulation_wrangler.get_primaries_pdg_code(-11)
        for positron in positrons:
            self.process_showers(positron, next(self.unique_topology))

    def process_gammas(self):
        """
        Primary gammas can come from ... cc_nue? nc? what else?
        For now, only sending to showers.
        """
        gammas = self.simulation_wrangler.get_primaries_abs_pdg_code(22)
        for gamma in gammas:
            self.process_showers(gamma, next(self.unique_topology))

    def process_muons(self):
        """
        Getting all muons, not just primaries, since they all need to
        be labelled the same.  Therefore, these physics_macro labels
        will need to be set somewhere else.  Or maybe that won't work
        with the unique labels.
        """
        muons = self.simulation_wrangler.get_trackid_pdg_code(13)   # this is fast // XO
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
                PhysicsMacroLabel.Undefined,
                muon_topology,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
                0
            )
            # process daughters (michel or delta, ignore the rest)
            muon_daughters = self.simulation_wrangler.trackid_daughters[muon]
            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(muon_daughters, 11)

            # process Michel electrons
            decay_daughters = self.simulation_wrangler.filter_trackid_process(elec_daughters, 6)
            decay_daughters += self.simulation_wrangler.filter_trackid_subprocess(elec_daughters, 151)
            michel_decay_hits = self.simulation_wrangler.get_hits_trackid(decay_daughters)
            michel_decay_segments = self.simulation_wrangler.get_segments_trackid(decay_daughters)
            self.simulation_wrangler.set_hit_labels_list(
                michel_decay_hits,
                michel_decay_segments,
                decay_daughters,
                TopologyLabel.Track,
                PhysicsMicroLabel.ElectronIonization,
                PhysicsMesoLabel.MichelElectron,
                PhysicsMacroLabel.Undefined,
                muon_topology,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
                0
            )
            # process deltas
            elec_em_daughters = self.simulation_wrangler.filter_trackid_process(elec_daughters, 2)
            delta_daughters = self.simulation_wrangler.filter_trackid_subprocess(elec_em_daughters, 2)
            delta_hits = self.simulation_wrangler.get_hits_trackid(delta_daughters)
            delta_segments = self.simulation_wrangler.get_segments_trackid(delta_daughters)
            self.simulation_wrangler.set_hit_labels_list(
                delta_hits,
                delta_segments,
                delta_daughters,
                TopologyLabel.Track,
                PhysicsMicroLabel.ElectronIonization,
                PhysicsMesoLabel.DeltaElectron,
                PhysicsMacroLabel.Undefined,
                muon_topology,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
                0
            )
            # process other electrons not delta and not michel
            # these are the not delta daughters
            other_daughters = list(filter(lambda x: x not in decay_daughters, muon_daughters))
            # not_delta_daughters = self.simulation_wrangler.filter_trackid_not_process_and_subprocess(muon_daughters, 2, 2)
            # # these are the not michel electrons of the not delta daughters
            # other_elec_daughters = self.simulation_wrangler.filter_trackid_not_process(not_delta_daughters, 6)

            self.process_showers_list(other_daughters, next(self.unique_topology))

    def process_anti_muons(self):
        muons = self.simulation_wrangler.get_trackid_pdg_code(-13)
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
                PhysicsMacroLabel.Undefined,
                muon_topology,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
                0
            )
            # process daughters (michel or delta, ignore the rest)
            muon_daughters = self.simulation_wrangler.trackid_daughters[muon]
            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(muon_daughters, 11)

            # process Michel electrons
            decay_daughters = self.simulation_wrangler.filter_trackid_process(elec_daughters, 6)
            decay_daughters += self.simulation_wrangler.filter_trackid_subprocess(elec_daughters, 151)
            michel_decay_hits = self.simulation_wrangler.get_hits_trackid(decay_daughters)
            michel_decay_segments = self.simulation_wrangler.get_segments_trackid(decay_daughters)
            self.simulation_wrangler.set_hit_labels_list(
                michel_decay_hits,
                michel_decay_segments,
                decay_daughters,
                TopologyLabel.Track,
                PhysicsMicroLabel.ElectronIonization,
                PhysicsMesoLabel.MichelElectron,
                PhysicsMacroLabel.Undefined,
                muon_topology,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
                0
            )
            # process deltas
            elec_em_daughters = self.simulation_wrangler.filter_trackid_process(elec_daughters, 2)
            delta_daughters = self.simulation_wrangler.filter_trackid_subprocess(elec_em_daughters, 2)
            delta_hits = self.simulation_wrangler.get_hits_trackid(delta_daughters)
            delta_segments = self.simulation_wrangler.get_segments_trackid(delta_daughters)
            self.simulation_wrangler.set_hit_labels_list(
                delta_hits,
                delta_segments,
                delta_daughters,
                TopologyLabel.Track,
                PhysicsMicroLabel.ElectronIonization,
                PhysicsMesoLabel.DeltaElectron,
                PhysicsMacroLabel.Undefined,
                muon_topology,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
                0
            )
            # process other electrons not delta and not michel
            # these are the not delta daughters
            other_daughters = list(filter(lambda x: x not in decay_daughters, muon_daughters))
            # not_delta_daughters = self.simulation_wrangler.filter_trackid_not_process_and_subprocess(muon_daughters, 2, 2)
            # # these are the not michel electrons of the not delta daughters
            # other_elec_daughters = self.simulation_wrangler.filter_trackid_not_process(not_delta_daughters, 6)

            self.process_showers_list(other_daughters, next(self.unique_topology))

    # Start adding from here. Not 100% sure I'm doing it correctly! 
    
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
            # # process pi0 decay into two photons
            # pi0_hits = self.simulation_wrangler.trackid_hit[pi0]
            # pi0_segments = self.simulation_wrangler.trackid_segmentid[pi0]

            # cluster_label = next(self.unique_topology)
            # self.simulation_wrangler.set_hit_labels(
            #     pi0_hits,
            #     pi0_segments,
            #     pi0,
            #     TopologyLabel.Shower,
            #     PhysicsMicroLabel.PhotonShower,
            #     cluster_label
            # )
            # label all pi0 descendants as photon showers
            # pi0_descendants = self.simulation_wrangler.trackid_descendants[pi0]
            # descendants_hits = self.simulation_wrangler.get_hits_trackid(pi0_descendants)
            # descendants_segments = self.simulation_wrangler.get_segments_trackid(pi0_descendants)
            # self.simulation_wrangler.set_hit_labels_list(
            #     descendants_hits, descendants_segments, pi0_descendants,
            #     TopologyLabel.Shower, PhysicsLabel.Pi0Shower,               # TODO: check this
            #     cluster_label
            # )
            # label pi0 descendants as showers
            # this is very inefficient for some reason...
            # pi0_daughters = self.simulation_wrangler.trackid_daughters[pi0]
            # pi0_progeny = self.simulation_wrangler.trackid_progeny[pi0]
            # self.process_showers_list(pi0_daughters, next(self.unique_topology)) # but these will not be
            # # labeled photon showers
            # self.process_showers_list(pi0_progeny, next(self.unique_topology)) # these will

    def process_pion_plus(self):
        # pi plus can decay into muons or electrons
        pipluses = self.simulation_wrangler.get_trackid_pdg_code(211)
        for piplus in pipluses:
            # label piplus as HIP ionization
            piplus_daughters = self.simulation_wrangler.trackid_daughters[piplus]
            # piplus_progeny = self.simulation_wrangler.trackid_progeny[piplus]
            piplus_hits = self.simulation_wrangler.trackid_hit[piplus]
            piplus_segments = self.simulation_wrangler.trackid_segmentid[piplus]

            track_label = next(self.unique_topology)
            self.simulation_wrangler.set_hit_labels(
                piplus_hits,
                piplus_segments,
                piplus,
                TopologyLabel.Track,
                PhysicsMicroLabel.HIPIonization,
                PhysicsMesoLabel.HIP,
                PhysicsMacroLabel.Undefined,
                track_label,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
                0
            )
            # label electron daughters as ..? showers?
            # elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(piplus_daughters, 11)
            # elec_segments = self.simulation_wrangler.get_segments_trackid(elec_daughters)
            # elec_hits = self.simulation_wrangler.get_hits_trackid(elec_daughters)
            # self.simulation_wrangler.set_labels_list(
            #     elec_hits, elec_segments, elec_daughters,
            #     TopologyLabel.Track, PhysicsLabel.ElectronShower, # TODO: check this
            #     track_label
            # )
            # # label muon daughters as ..? tracks?
            # muon_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(piplus_daughters, 13)
            # muon_segments = self.simulation_wrangler.get_segments_trackid(muon_daughters)
            # muon_hits = self.simulation_wrangler.get_hits_trackid(muon_daughters)
            # self.simulation_wrangler.set_labels_list(
            #     muon_hits, muon_segments, muon_daughters,
            #     TopologyLabel.Track, PhysicsLabel.MIPIonization, # TODO: check this
            #     next(self.unique_topology)
            # )
            # label all other daughters and progeny as showers
            # not_elec_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(piplus_daughters, 11)
            # other_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(not_elec_daughters, 13)
            # self.process_showers_list(other_daughters, next(self.unique_topology))
            self.process_showers_list(piplus_daughters, next(self.unique_topology))

    def process_pion_minus(self):
        pimins = self.simulation_wrangler.get_trackid_pdg_code(-211)
        for pimin in pimins:
            # label pimin as HIP ionization
            pimin_daughters = self.simulation_wrangler.trackid_daughters[pimin]
            # pimin_progeny = self.simulation_wrangler.trackid_progeny[pimin]
            pimin_hits = self.simulation_wrangler.trackid_hit[pimin]
            pimin_segments = self.simulation_wrangler.trackid_segmentid[pimin]

            track_label = next(self.unique_topology)
            self.simulation_wrangler.set_hit_labels(
                pimin_hits,
                pimin_segments,
                pimin,
                TopologyLabel.Track,
                PhysicsMicroLabel.HIPIonization,
                PhysicsMesoLabel.HIP,
                PhysicsMacroLabel.Undefined,
                track_label,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
                0
            )
            # # elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(pimin_daughters, 11)
            # # other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(pimin_daughters, 11)

            # # elec_segments = self.simulation_wrangler.get_segments_trackid(elec_daughters)
            # # elec_hits = self.simulation_wrangler.get_hits_trackid(elec_daughters)

            # track_label = next(self.unique_topology)
            # self.simulation_wrangler.set_hit_labels(
            #     pimin_hits, pimin_segments, pimin,
            #     TopologyLabel.Track, PhysicsLabel.HIPIonization,
            #     track_label
            # )

            # # self.simulation_wrangler.set_labels_list(
            # #     elec_hits, elec_segments, elec_daughters,
            # #     TopologyLabel.Track, PhysicsLabel.HIPIonization,
            # #     track_label
            # # )

            # # self.process_showers_list(other_daughters, next(self.unique_topology))
            # # self.process_showers_list(pimin_progeny, next(self.unique_topology))
            self.process_showers_list(pimin_daughters, next(self.unique_topology))

    def process_kaon0s(self):
        ka0s = self.simulation_wrangler.get_trackid_pdg_code(311)
        pass
        # for ka0 in ka0s:
        #     ka0_daughters = self.simulation_wrangler.trackid_daughters[ka0]
        #     ka0_hits = self.simulation_wrangler.trackid_hit[ka0]
        #     ka0_segments = self.simulation_wrangler.trackid_segmentid[ka0]

        #     cluster_label = next(self.unique_topology)
        #     self.simulation_wrangler.set_hit_labels(
        #         ka0_hits, ka0_segments, ka0,
        #         TopologyLabel.Shower, PhysicsLabel.PhotonShower,
        #         cluster_label
        #     )
        #     self.process_showers_list(ka0_daughters, next(self.unique_topology))

    def process_kaon_plus(self):

        kaonpluses = self.simulation_wrangler.get_trackid_pdg_code(321)
        pass
        # for kaonplus in kaonpluses:
        #     kaonplus_daughters = self.simulation_wrangler.trackid_daughters[kaonplus]
        #     # elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(kaonplus_daughters, 11)
        #     # other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(kaonplus_daughters, 11)
        #     # elec_hits = self.simulation_wrangler.get_hits_trackid(elec_daughters)
        #     # elec_segments = self.simulation_wrangler.get_segments_trackid(elec_daughters)

        #     # kaonplus_progeny = self.simulation_wrangler.trackid_progeny[kaonplus]
        #     kaonplus_hits = self.simulation_wrangler.get_hits_trackid(kaonplus)
        #     kaonplus_segments = self.simulation_wrangler.get_segments_trackid(kaonplus)

        #     # Set kaonplus detsim labels to Track:kaononMinus
        #     track_label = next(self.unique_topology)
        #     self.simulation_wrangler.set_hit_labels(
        #         kaonplus_hits, kaonplus_segments, kaonplus,
        #         TopologyLabel.Track, PhysicsLabel.HIPIonization,
        #         track_label
        #     )

        #     # self.simulation_wrangler.set_labels_list(
        #     #     elec_hits, elec_segments, elec_daughters,
        #     #     TopologyLabel.Track, PhysicsLabel.HIPIonization,
        #     #     track_label
        #     # )

        #     # self.process_showers_list(other_daughters, next(self.unique_topology))
        #     # self.process_showers_list(kaonplus_progeny, next(self.unique_topology))
        #     self.process_showers_list(kaonplus_daughters, next(self.unique_topology))

    def process_kaon_minus(self):
        kaonmins = self.simulation_wrangler.get_trackid_pdg_code(-321)
        pass
        # for kaonmin in kaonmins:
        #     kaonmin_daughters = self.simulation_wrangler.trackid_daughters[kaonmin]
        #     # elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(kaonmin_daughters, 11)
        #     # other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(kaonmin_daughters, 11)
        #     # elec_hits = self.simulation_wrangler.get_hits_trackid(elec_daughters)
        #     # elec_segments = self.simulation_wrangler.get_segments_trackid(elec_daughters)

        #     # kaonmin_progeny = self.simulation_wrangler.trackid_progeny[kaonmin]
        #     kaonmin_hits = self.simulation_wrangler.get_hits_trackid(kaonmin)
        #     kaonmin_segments = self.simulation_wrangler.get_segments_trackid(kaonmin)

        #     # Set kaonplus detsim labels to Track:kaononMinus
        #     track_label = next(self.unique_topology)
        #     self.simulation_wrangler.set_hit_labels(
        #         kaonmin_hits, kaonmin_segments, kaonmin,
        #         TopologyLabel.Track, PhysicsLabel.HIPIonization,
        #         track_label
        #     )

        #     # self.simulation_wrangler.set_labels_list(
        #     #     elec_hits, elec_segments, elec_daughters,
        #     #     TopologyLabel.Track, PhysicsLabel.HIPIonization,
        #     #     track_label
        #     # )

        #     # self.process_showers_list(other_daughters, next(self.unique_topology))
        #     # self.process_showers_list(kaonmin_progeny, next(self.unique_topology))
        #     self.process_showers_list(kaonmin_daughters, next(self.unique_topology))

    def process_protons(self):
        protons = self.simulation_wrangler.get_trackid_pdg_code(2212)
        for proton in protons:
            # label protons as HIP ionization
            proton_hits = self.simulation_wrangler.get_hits_trackid(proton)
            proton_segments = self.simulation_wrangler.get_segments_trackid(proton)
            cluster_label = next(self.unique_topology)
            self.simulation_wrangler.set_hit_labels(
                proton_hits,
                proton_segments,
                proton,
                TopologyLabel.Track,
                PhysicsMicroLabel.HIPIonization,
                PhysicsMesoLabel.HIP,
                PhysicsMacroLabel.Undefined,
                cluster_label,
                next(self.unique_physics_micro),
                next(self.unique_physics_meso),
                0
            )
            proton_daughters = self.simulation_wrangler.trackid_daughters[proton]
            # elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(proton_daughters, 11)

            # # process deltas
            # elec_em_daughters = self.simulation_wrangler.filter_trackid_process(elec_daughters, 2)
            # delta_daughters = self.simulation_wrangler.filter_trackid_subprocess(elec_em_daughters, 2)
            # delta_hits = self.simulation_wrangler.get_hits_trackid(delta_daughters)
            # delta_segments = self.simulation_wrangler.get_segments_trackid(delta_daughters)
            # self.simulation_wrangler.set_labels_list(
            #     delta_hits, delta_segments, delta_daughters,
            #     TopologyLabel.Track, PhysicsLabel.DeltaElectron,
            #     next(self.unique_topology)
            # )

            # # process all other daughters as showers too
            # # process all other electrons as showers
            # other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(proton_daughters, 11)
            # other_elec_daughters = self.simulation_wrangler.filter_trackid_not_process_and_subprocess(elec_daughters, 2, 2)

            # self.process_showers_list(other_elec_daughters, next(self.unique_topology))

            # self.process_showers_list(other_daughters, next(self.unique_topology))
            self.process_showers_list(proton_daughters, next(self.unique_topology))
            # we do not need these anymore (I think) # TODO: check this
            # proton_progeny = self.simulation_wrangler.trackid_progeny[proton]
            # self.process_showers_list(proton_progeny, next(self.unique_topology))

    def process_neutrons(self):
        """
        Process all different types of neutron interactions, which include
        (1) neutron elastic
        (2) neutron inelastic (which can lead to a proton)
        (3) neutron captures on ar40/ar38/ar36
        (4) fission

        TODO: Only doing ar40 captures right now, need to update this.
        """
        neutrons = self.simulation_wrangler.get_trackid_pdg_code(2112)
        pass
        # for neutron in neutrons:
        #     neutron_daughters = self.simulation_wrangler.trackid_daughters[neutron]
        #     gamma_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(neutron_daughters, 22)
        #     other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(neutron_daughters, 22)
        #     capture_daughters = self.simulation_wrangler.filter_trackid_process(gamma_daughters, 131)
        #     # 131 = GEANT process capture
        #     other_gammas = self.simulation_wrangler.filter_trackid_not_process(gamma_daughters, 131)

        #     for capture in capture_daughters:
        #         for gamma in capture:
        #             gamma_energy = self.simulation_wrangler.get_total_hit_energy(gamma, 5)
        #             gamma_hits = self.simulation_wrangler.get_hits_trackid(gamma)
        #             gamma_segments = self.simulation_wrangler.get_segments_trackid(gamma)

        #             particle_label = PhysicsLabel.NeutronCaptureGammaOther
        #             if gamma_energy in {0.00474, 0.00475}:
        #                 particle_label = PhysicsLabel.NeutronCaptureGamma474
        #             elif gamma_energy in {0.00336, 0.00337}:
        #                 particle_label = PhysicsLabel.NeutronCaptureGamma336
        #             elif gamma_energy in {0.00256, 0.00257}:
        #                 particle_label = PhysicsLabel.NeutronCaptureGamma256
        #             elif gamma_energy in {0.00118, 0.00119}:
        #                 particle_label = PhysicsLabel.NeutronCaptureGamma118
        #             elif gamma_energy in {0.00083, 0.00084}:
        #                 particle_label = PhysicsLabel.NeutronCaptureGamma083
        #             elif gamma_energy in {0.00051, 0.00052}:
        #                 particle_label = PhysicsLabel.NeutronCaptureGamma051
        #             elif gamma_energy in {0.00016, 0.00017}:
        #                 particle_label = PhysicsLabel.NeutronCaptureGamma016

        #             cluster_label = next(self.unique_topology)
        #             self.simulation_wrangler.set_hit_labels(
        #                 gamma_hits, gamma_segments, gamma,
        #                 TopologyLabel.Blip, particle_label,
        #                 cluster_label
        #             )

        #     self.process_showers_list(other_daughters, next(self.unique_topology))
        #     self.process_showers_list(other_gammas, next(self.unique_topology))

    def process_nuclear_recoils(self):
        pass
        # ar41 = self.simulation_wrangler.get_trackid_pdg_code(1000180410)
        # ar40 = self.simulation_wrangler.get_trackid_pdg_code(1000180400)
        # ar39 = self.simulation_wrangler.get_trackid_pdg_code(1000180390)
        # ar38 = self.simulation_wrangler.get_trackid_pdg_code(1000180380)
        # ar37 = self.simulation_wrangler.get_trackid_pdg_code(1000180370)
        # ar36 = self.simulation_wrangler.get_trackid_pdg_code(1000180360)

        # # ar41_daughters = self.simulation_wrangler.get_daughters_trackid(ar41)
        # # ar40_daughters = self.simulation_wrangler.get_daughters_trackid(ar40)
        # # ar39_daughters = self.simulation_wrangler.get_daughters_trackid(ar39)
        # # ar38_daughters = self.simulation_wrangler.get_daughters_trackid(ar38)
        # # ar37_daughters = self.simulation_wrangler.get_daughters_trackid(ar37)
        # # ar36_daughters = self.simulation_wrangler.get_daughters_trackid(ar36)

        # s33 = self.simulation_wrangler.get_trackid_pdg_code(1000160330)
        # # s33_daughters = self.simulation_wrangler.get_daughters_trackid(s33)
        # s35 = self.simulation_wrangler.get_trackid_pdg_code(1000160350)
        # # s35_daughters = self.simulation_wrangler.get_daughters_trackid(s35)
        # s36 = self.simulation_wrangler.get_trackid_pdg_code(1000160360)
        # # s36_daughters = self.simulation_wrangler.get_daughters_trackid(s36)

        # cl36 = self.simulation_wrangler.get_trackid_pdg_code(1000170360)
        # # cl36_daughters = self.simulation_wrangler.get_daughters_trackid(cl36)
        # cl37 = self.simulation_wrangler.get_trackid_pdg_code(1000170370)
        # # cl37_daughters = self.simulation_wrangler.get_daughters_trackid(cl37)
        # cl39 = self.simulation_wrangler.get_trackid_pdg_code(1000170390)
        # # cl39_daughters = self.simulation_wrangler.get_daughters_trackid(cl39)
        # cl40 = self.simulation_wrangler.get_trackid_pdg_code(1000170400)
        # # cl40_daughters = self.simulation_wrangler.get_daughters_trackid(cl40)

        # for ar in ar41:
        #     ar41_hits = self.simulation_wrangler.get_hits_trackid(ar)
        #     ar41_segments = self.simulation_wrangler.get_segments_trackid(ar)
        #     self.simulation_wrangler.set_hit_labels(
        #         ar41_hits, ar41_segments, ar,
        #         TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
        #         next(self.unique_topology)
        #     )

        # for ar in ar40:
        #     ar40_hits = self.simulation_wrangler.get_hits_trackid(ar)
        #     ar40_segments = self.simulation_wrangler.get_segments_trackid(ar)
        #     self.simulation_wrangler.set_hit_labels(
        #         ar40_hits, ar40_segments, ar,
        #         TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
        #         next(self.unique_topology)
        #     )

        # for ar in ar39:
        #     ar39_hits = self.simulation_wrangler.get_hits_trackid(ar)
        #     ar39_segments = self.simulation_wrangler.get_segments_trackid(ar)
        #     self.simulation_wrangler.set_hit_labels(
        #         ar39_hits, ar39_segments, ar,
        #         TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
        #         next(self.unique_topology)
        #     )

        # for ar in ar38:
        #     ar38_hits = self.simulation_wrangler.get_hits_trackid(ar)
        #     ar38_segments = self.simulation_wrangler.get_segments_trackid(ar)
        #     self.simulation_wrangler.set_hit_labels(
        #         ar38_hits, ar38_segments, ar,
        #         TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
        #         next(self.unique_topology)
        #     )

        # for ar in ar37:
        #     ar37_hits = self.simulation_wrangler.get_hits_trackid(ar)
        #     ar37_segments = self.simulation_wrangler.get_segments_trackid(ar)
        #     self.simulation_wrangler.set_hit_labels(
        #         ar37_hits, ar37_segments, ar,
        #         TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
        #         next(self.unique_topology)
        #     )

        # for ar in ar36:
        #     ar36_hits = self.simulation_wrangler.get_hits_trackid(ar)
        #     ar36_segments = self.simulation_wrangler.get_segments_trackid(ar)
        #     self.simulation_wrangler.set_hit_labels(
        #         ar36_hits, ar36_segments, ar,
        #         TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
        #         next(self.unique_topology)
        #     )

        # for s in s33:
        #     s33_hits = self.simulation_wrangler.get_hits_trackid(s)
        #     s33_segments = self.simulation_wrangler.get_segments_trackid(s)
        #     self.simulation_wrangler.set_hit_labels(
        #         s33_hits, s33_segments, s,
        #         TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
        #         next(self.unique_topology)
        #     )

        # for s in s35:
        #     s35_hits = self.simulation_wrangler.get_hits_trackid(s)
        #     s35_segments = self.simulation_wrangler.get_segments_trackid(s)
        #     self.simulation_wrangler.set_hit_labels(
        #         s35_hits, s35_segments, s,
        #         TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
        #         next(self.unique_topology)
        #     )

        # for s in s36:
        #     s36_hits = self.simulation_wrangler.get_hits_trackid(s)
        #     s36_segments = self.simulation_wrangler.get_segments_trackid(s)
        #     self.simulation_wrangler.set_hit_labels(
        #         s36_hits, s36_segments, s,
        #         TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
        #         next(self.unique_topology)
        #     )

        # for cl in cl36:
        #     cl36_hits = self.simulation_wrangler.get_hits_trackid(cl)
        #     cl36_segments = self.simulation_wrangler.get_segments_trackid(cl)
        #     self.simulation_wrangler.set_hit_labels(
        #         cl36_hits, cl36_segments, cl,
        #         TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
        #         next(self.unique_topology)
        #     )

        # for cl in cl37:
        #     cl37_hits = self.simulation_wrangler.get_hits_trackid(cl)
        #     cl37_segments = self.simulation_wrangler.get_segments_trackid(cl)
        #     self.simulation_wrangler.set_hit_labels(
        #         cl37_hits, cl37_segments, cl,
        #         TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
        #         next(self.unique_topology)
        #     )

        # for cl in cl39:
        #     cl39_hits = self.simulation_wrangler.get_hits_trackid(cl)
        #     cl39_segments = self.simulation_wrangler.get_segments_trackid(cl)
        #     self.simulation_wrangler.set_hit_labels(
        #         cl39_hits, cl39_segments, cl,
        #         TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
        #         next(self.unique_topology)
        #     )

        # for cl in cl40:
        #     cl40_hits = self.simulation_wrangler.get_hits_trackid(cl)
        #     cl40_segments = self.simulation_wrangler.get_segments_trackid(cl)
        #     self.simulation_wrangler.set_hit_labels(
        #         cl40_hits, cl40_segments, cl,
        #         TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
        #         next(self.unique_topology)
        #     )

    def process_electron_recoils(self):
        pass
        # deuterons = self.simulation_wrangler.get_trackid_pdg_code(1000010020)
        # tritons = self.simulation_wrangler.get_trackid_pdg_code(1000010030)
        # alphas = self.simulation_wrangler.get_trackid_pdg_code(1000020040)

        # inelastic_alphas = self.simulation_wrangler.filter_trackid_process(alphas, PhysicsLabel.NuclearRecoil)
        # # TODO: fix this

        # # deuteron_daughters = self.simulation_wrangler.get_daughters_trackid(deuterons)
        # # triton_daughters = self.simulation_wrangler.get_daughters_trackid(tritons)
        # # inelastic_alpha_daughters = self.simulation_wrangler.get_daughters_trackid(inelastic_alphas)

        # for deuteron in deuterons:
        #     deuteron_hits = self.simulation_wrangler.get_hits_trackid(deuteron)
        #     deuteron_segments = self.simulation_wrangler.get_segments_trackid(deuteron)

        #     self.simulation_wrangler.set_hit_labels(
        #         deuteron_hits, deuteron_segments, deuteron,
        #         TopologyLabel.Blip, PhysicsLabel.ElectronRecoil,
        #         next(self.unique_topology)
        #     )

        # for triton in tritons:
        #     triton_hits = self.simulation_wrangler.get_hits_trackid(triton)
        #     triton_segments = self.simulation_wrangler.get_segments_trackid(triton)

        #     self.simulation_wrangler.set_hit_labels(
        #         triton_hits, triton_segments, triton,
        #         TopologyLabel.Blip, PhysicsLabel.ElectronRecoil,
        #         next(self.unique_topology)
        #     )

        # for inelastic_alpha in inelastic_alphas:
        #     inelastic_alpha_hits = self.simulation_wrangler.get_hits_trackid(inelastic_alpha)
        #     inelastic_alpha_segments = self.simulation_wrangler.get_segments_trackid(inelastic_alpha)

        #     self.simulation_wrangler.set_hit_labels(
        #         inelastic_alpha_hits, inelastic_alpha_segments, inelastic_alpha,
        #         TopologyLabel.Blip, PhysicsLabel.ElectronRecoil,
        #         next(self.unique_topology)
        #     )

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
        See here for some details: https://indico.fnal.gov/event/50121/contributions/220205/attachments/145404/185102/20210721_Decay0_Lasorak.pdf

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
        """
        """
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