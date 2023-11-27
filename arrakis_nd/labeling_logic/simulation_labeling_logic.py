"""
"""
from arrakis_nd.wrangler.simulation_wrangler import SimulationWrangler
from arrakis_nd.dataset.common import *

class SimulationLabelingLogic:
    """
    """
    def __init__(self,
        simulation_wrangler:    SimulationWrangler=None
    ):
        self.simulation_wrangler = simulation_wrangler
        self.shower_threshold = 20

        self.topology_label = 0
    
    def iterate_topology_label(self):
        self.topology_label += 1
        return self.topology_label
    
    def set_labels(self,
        hits, segments, trackid,
        topology, physics, unique_topology
    ):
        for hit in hits:
            self.simulation_wrangler.set_hit_labels(
                hit, 
                trackid, 
                topology, 
                trackid, 
                physics, 
                unique_topology
            )
    
    def set_labels_list(self,
        hits, segments, trackid,
        topology, physics, unique_topology
    ):
        for ii in range(len(hits)):
            self.set_labels(
                hits[ii],
                segments[ii], 
                trackid[ii], 
                topology, 
                physics, 
                unique_topology
            )
    
    def set_labels_array(self,
        hits, segments, trackid,
        topology, physics, unique_topology
    ):
        for ii in range(len(hits)):
            self.set_labels_list(
                hits[ii],
                segments[ii], 
                trackid[ii], 
                topology, 
                physics, 
                unique_topology
            )

    def process_event(self):
        
        self.topology_label = 0
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
        self.process_neutron_captures()
        self.process_nuclear_recoils()
        self.process_electron_recoils()
        self.process_ar39()
        self.process_ar42()
        self.process_kr85()
        self.process_rn222()
        self.process_cosmics()

    def process_showers(self,
        particle, topology_label
    ):
        particle_hits = self.simulation_wrangler.trackid_hit[particle]
        particle_segments = self.simulation_wrangler.trackid_segmentid[particle]
        if (
            self.simulation_wrangler.trackid_process[particle] == "primary" or
            self.simulation_wrangler.trackid_process[particle] == "conv" or
            self.simulation_wrangler.trackid_process[particle] == "compt" or
            self.simulation_wrangler.get_total_hit_energy(particle_hits) >= self.shower_threshold
        ):
            if self.simulation_wrangler.trackid_pdgcode[particle] == 11:
                self.set_labels(
                    particle_hits, particle_segments, particle,
                    TopologyLabel.Shower, PhysicsLabel.ElectronShower,
                    topology_label
                )
            elif self.simulation_wrangler.trackid_pdgcode[particle] == -11:
                self.set_labels(
                    particle_hits, particle_segments, particle,
                    TopologyLabel.Shower, PhysicsLabel.PositronShower,
                    topology_label
                )
            elif abs(self.simulation_wrangler.trackid_pdgcode[particle]) == 22:
                self.set_labels(
                    particle_hits, particle_segments, particle,
                    TopologyLabel.Shower, PhysicsLabel.PhotonShower,
                    topology_label
                )
        else:
            self.set_labels(
                particle_hits, particle_segments, particle,
                TopologyLabel.Blip, PhysicsLabel.ElectronRecoil,
                topology_label
            )
        daughters = self.simulation_wrangler.trackid_daughters[particle]
        elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(daughters, 11)
        photon_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(daughters, 22)
        self.process_showers_list(elec_daughters, topology_label)
        self.process_showers_list(photon_daughters, topology_label)
    
    def process_showers_list(self,
        particles, topology_label
    ):
        for particle in particles:
            self.process_showers(particle, topology_label)
    
    def process_showers_array(self,
        particles
    ):
        for particle in particles:
            self.process_showers_list(particle, self.iterate_topology_label())

    def process_electrons(self):
        electrons = self.simulation_wrangler.get_primaries_pdg_code(11)
        for electron in electrons:
            electron_daughters = self.simulation_wrangler.trackid_daughters[electron]
            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(electron_daughters, 11)
            other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(electron_daughters, 11)
            elec_hits = self.simulation_wrangler.get_hits_trackid(elec_daughters)
            elec_segments = self.simulation_wrangler.get_segments_trackid(elec_daughters)

            electron_progeny = self.simulation_wrangler.trackid_progeny[electron]
            electron_hits = self.simulation_wrangler.trackid_hit[electron]
            electron_segment = self.simulation_wrangler.trackid_segmentid[electron]

            shower_label = self.iterate_topology_label()
            self.set_labels(
                electron_hits, electron_segment, electron,
                TopologyLabel.Shower, PhysicsLabel.ElectronShower,
                shower_label
            )
            self.set_labels_list(
                elec_hits, elec_segments, elec_daughters,
                TopologyLabel.Shower, PhysicsLabel.ElectronShower,
                shower_label
            )
            self.process_showers_list(electron_progeny, shower_label)
            self.process_showers_list(other_daughters, self.iterate_topology_label())

    def process_positrons(self):
        positrons = self.simulation_wrangler.get_primaries_pdg_code(-11)
        for positron in positrons:
            positron_daughters = self.simulation_wrangler.trackid_daughters[positron]
            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(positron_daughters, 11)
            other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(positron_daughters, 11)
            elec_hits = self.simulation_wrangler.get_hits_trackid(elec_daughters)
            elec_segments = self.simulation_wrangler.get_segments_trackid(elec_daughters)

            positron_progeny = self.simulation_wrangler.trackid_progeny[positron]
            positron_hits = self.simulation_wrangler.trackid_hit[positron]
            positron_segment = self.simulation_wrangler.trackid_segmentid[positron]

            shower_label = self.iterate_topology_label()
            self.set_labels(
                positron_hits, positron_segment, positron,
                TopologyLabel.Shower, PhysicsLabel.ElectronShower,
                shower_label
            )
            self.set_labels_list(
                elec_hits, elec_segments, elec_daughters,
                TopologyLabel.Shower, PhysicsLabel.ElectronShower,
                shower_label
            )
            self.process_showers_list(positron_progeny, shower_label)
            self.process_showers_list(other_daughters, self.iterate_topology_label())
        
    def process_gammas(self):
        gammas = self.simulation_wrangler.get_primaries_abs_pdg_code(22)
        for gamma in gammas:
            gamma_daughters = self.simulation_wrangler.trackid_daughters[gamma]
            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(gamma_daughters, 11)
            other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(gamma_daughters, 11)
            elec_hits = self.simulation_wrangler.get_hits_trackid(elec_daughters)
            elec_segments = self.simulation_wrangler.get_segments_trackid(elec_daughters)

            gamma_progeny = self.simulation_wrangler.trackid_progeny[gamma]
            gamma_hits = self.simulation_wrangler.trackid_hit[gamma]
            gamma_segment = self.simulation_wrangler.trackid_segmentid[gamma]

            shower_label = self.iterate_topology_label()
            self.set_labels(
                gamma_hits, gamma_segment, gamma,
                TopologyLabel.Shower, PhysicsLabel.ElectronShower,
                shower_label
            )
            self.set_labels_list(
                elec_hits, elec_segments, elec_daughters,
                TopologyLabel.Shower, PhysicsLabel.ElectronShower,
                shower_label
            )
            self.process_showers_list(gamma_progeny, shower_label)

    def process_muons(self):
        muons = self.simulation_wrangler.get_trackid_pdg_code(13)
        for muon in muons:
            muon_daughters = self.simulation_wrangler.trackid_daughters[muon]
            muon_progeny = self.simulation_wrangler.trackid_progeny[muon]
            muon_hits = self.simulation_wrangler.trackid_hit[muon]
            muon_segments = self.simulation_wrangler.trackid_segmentid[muon]
            cluster_label = self.iterate_topology_label()
            self.set_labels(
                muon_hits, muon_segments, muon,
                TopologyLabel.Track, PhysicsLabel.MIPIonization,
                cluster_label
            )
            
            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(muon_daughters, 11)
            other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(muon_daughters, 11)

            decay_daughters = self.simulation_wrangler.filter_trackid_process(elec_daughters, "decay")
            capture_daughters = self.simulation_wrangler.filter_trackid_process(elec_daughters, "muonCaptureAtRest")
            michel_decay_hits = self.simulation_wrangler.get_hits_trackid(decay_daughters)
            michel_capture_hits = self.simulation_wrangler.get_hits_trackid(capture_daughters)
            michel_decay_segments = self.simulation_wrangler.get_segments_trackid(decay_daughters)
            michel_capture_segments = self.simulation_wrangler.get_segments_trackid(capture_daughters)

            delta_daughters = self.simulation_wrangler.filter_trackid_process(elec_daughters, "muIoni")
            delta_hits = self.simulation_wrangler.get_hits_trackid(delta_daughters)
            delta_segments = self.simulation_wrangler.get_segments_trackid(delta_daughters)

            not_decay_elec_daughters = self.simulation_wrangler.filter_trackid_not_process(elec_daughters, "decay")
            not_muon_capture_elec_daughters = self.simulation_wrangler.filter_trackid_not_process(not_decay_elec_daughters, "muonCaptureAtRest")
            other_elec_daughters = self.simulation_wrangler.filter_trackid_not_process(not_muon_capture_elec_daughters, "muIoni")

            self.set_labels_list(
                michel_decay_hits, michel_decay_segments, decay_daughters,
                TopologyLabel.Track, PhysicsLabel.MichelElectron,
                self.iterate_topology_label()
            )
            self.set_labels_list(
                michel_capture_hits, michel_capture_segments, capture_daughters,
                TopologyLabel.Track, PhysicsLabel.MichelElectron,
                self.iterate_topology_label()
            )
            self.set_labels_list(
                delta_hits, delta_segments, delta_daughters,
                TopologyLabel.Track, PhysicsLabel.DeltaElectron,
                self.iterate_topology_label()
            )

            michel_decay_elec_daughters = self.simulation_wrangler.get_daughters_trackid(decay_daughters)
            michel_capture_elec_daughters = self.simulation_wrangler.get_daughters_trackid(capture_daughters)
            delta_elec_daughters = self.simulation_wrangler.get_daughters_trackid(delta_daughters)

            self.process_showers_array(michel_decay_elec_daughters)
            self.process_showers_array(michel_capture_elec_daughters)
            self.process_showers_array(delta_elec_daughters)
            self.process_showers_list(other_elec_daughters, self.iterate_topology_label())
            self.process_showers_list(other_daughters, self.iterate_topology_label())
            self.process_showers_list(muon_progeny, self.iterate_topology_label())
    
    def process_anti_muons(self):
        muons = self.simulation_wrangler.get_trackid_pdg_code(-13)
        for muon in muons:
            muon_daughters = self.simulation_wrangler.trackid_daughters[muon]
            muon_progeny = self.simulation_wrangler.trackid_progeny[muon]
            muon_hits = self.simulation_wrangler.trackid_hit[muon]
            muon_segments = self.simulation_wrangler.trackid_segmentid[muon]

            cluster_label = self.iterate_topology_label()
            self.set_labels(
                muon_hits, muon_segments, muon,
                TopologyLabel.Track, PhysicsLabel.MIPIonization,
                cluster_label
            )
            
            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(muon_daughters, 11)
            other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(muon_daughters, 11)

            decay_daughters = self.simulation_wrangler.filter_trackid_process(elec_daughters, "decay")
            capture_daughters = self.simulation_wrangler.filter_trackid_process(elec_daughters, "muonCaptureAtRest")
            michel_decay_hits = self.simulation_wrangler.get_hits_trackid(decay_daughters)
            michel_capture_hits = self.simulation_wrangler.get_hits_trackid(capture_daughters)
            michel_decay_segments = self.simulation_wrangler.get_segments_trackid(decay_daughters)
            michel_capture_segments = self.simulation_wrangler.get_segments_trackid(capture_daughters)

            delta_daughters = self.simulation_wrangler.filter_trackid_process(elec_daughters, "muIoni")
            delta_hits = self.simulation_wrangler.get_hits_trackid(delta_daughters)
            delta_segments = self.simulation_wrangler.get_segments_trackid(delta_daughters)

            not_decay_elec_daughters = self.simulation_wrangler.filter_trackid_not_process(elec_daughters, "decay")
            not_muon_capture_elec_daughters = self.simulation_wrangler.filter_trackid_not_process(not_decay_elec_daughters, "muonCaptureAtRest")
            other_elec_daughters = self.simulation_wrangler.filter_trackid_not_process(not_muon_capture_elec_daughters, "muIoni")

            self.set_labels_list(
                michel_decay_hits, michel_decay_segments, decay_daughters,
                TopologyLabel.Track, PhysicsLabel.MichelElectron,
                self.iterate_topology_label()
            )
            self.set_labels_list(
                michel_capture_hits, michel_capture_segments, capture_daughters,
                TopologyLabel.Track, PhysicsLabel.MichelElectron,
                self.iterate_topology_label()
            )
            self.set_labels_list(
                delta_hits, delta_segments, delta_daughters,
                TopologyLabel.Track, PhysicsLabel.DeltaElectron,
                self.iterate_topology_label()
            )

            michel_decay_elec_daughters = self.simulation_wrangler.get_daughters_trackid(decay_daughters)
            michel_capture_elec_daughters = self.simulation_wrangler.get_daughters_trackid(capture_daughters)
            delta_elec_daughters = self.simulation_wrangler.get_daughters_trackid(delta_daughters)

            self.process_showers_array(michel_decay_elec_daughters)
            self.process_showers_array(michel_capture_elec_daughters)
            self.process_showers_array(delta_elec_daughters)
            self.process_showers_list(other_elec_daughters, self.iterate_topology_label())
            self.process_showers_list(other_daughters, self.iterate_topology_label())
            self.process_showers_list(muon_progeny, self.iterate_topology_label())

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
        pi0s = self.simulation_wrangler.get_trackid_pdg_code(111)
        for pi0 in pi0s:
            pi0_daughters = self.simulation_wrangler.trackid_daughters[pi0]
            pi0_progeny = self.simulation_wrangler.trackid_progeny[pi0]
            pi0_hits = self.simulation_wrangler.trackid_hit[pi0]
            pi0_segments = self.simulation_wrangler.trackid_segmentid[pi0]

            cluster_label = self.iterate_topology_label()
            self.set_labels(
                pi0_hits, pi0_segments, pi0,
                TopologyLabel.Shower, PhysicsLabel.PhotonShower,
                cluster_label
            )

    def process_pion_plus(self):
        pipluses = self.simulation_wrangler.get_trackid_pdg_code(211)
        for piplus in pipluses:
            piplus_daughters = self.simulation_wrangler.trackid_daughters[piplus]
            piplus_progeny = self.simulation_wrangler.trackid_progeny[piplus]
            piplus_hits = self.simulation_wrangler.trackid_hit[piplus]
            piplus_segments = self.simulation_wrangler.trackid_segmentid[piplus]

            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(piplus_daughters, 11)
            other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(piplus_daughters, 11)
            
            elec_segments = self.simulation_wrangler.get_segments_trackid(elec_daughters)
            elec_hits = self.simulation_wrangler.get_hits_trackid(elec_daughters)

            track_label = self.iterate_topology_label()
            self.set_labels(
                piplus_hits, piplus_segments, piplus,
                TopologyLabel.Track, PhysicsLabel.HIPIonization,
                track_label
            )
            
            self.set_labels(
                elec_hits, elec_segments, elec_daughters,
                TopologyLabel.Track, PhysicsLabel.HIPIonization,
                track_label
            )
            
            self.process_showers_list(other_daughters, self.iterate_topology_label())
            self.process_showers_list(piplus_progeny, self.iterate_topology_label())

    def process_pion_minus(self):
        pimins = self.simulation_wrangler.get_trackid_pdg_code(-211)
        for pimin in pimins:
            pimin_daughters = self.simulation_wrangler.trackid_daughters[pimin]
            pimin_progeny = self.simulation_wrangler.trackid_progeny[pimin]
            pimin_hits = self.simulation_wrangler.trackid_hit[pimin]
            pimin_segments = self.simulation_wrangler.trackid_segmentid[pimin]

            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(pimin_daughters, 11)
            other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(pimin_daughters, 11)
            
            elec_segments = self.simulation_wrangler.get_segments_trackid(elec_daughters)
            elec_hits = self.simulation_wrangler.get_hits_trackid(elec_daughters)

            track_label = self.iterate_topology_label()
            self.set_labels(
                pimin_hits, pimin_segments, pimin,
                TopologyLabel.Track, PhysicsLabel.HIPIonization,
                track_label
            )
            
            self.set_labels(
                elec_hits, elec_segments, elec_daughters,
                TopologyLabel.Track, PhysicsLabel.HIPIonization,
                track_label
            )
            
            self.process_showers_list(other_daughters, self.iterate_topology_label())
            self.process_showers_list(pimin_progeny, self.iterate_topology_label())

    def process_kaon0s(self):
        ka0s = self.simulation_wrangler.get_trackid_pdg_code(311)
        for ka0 in ka0s:
            ka0_daughters = self.simulation_wrangler.trackid_daughters[ka0]
            ka0_hits = self.simulation_wrangler.trackid_hits[ka0]
            ka0_segments = self.simulation_wrangler.trackid_segments[ka0]

            cluster_label = self.iterate_topology_label()
            self.set_labels(
                ka0_hits, ka0_segments, ka0,
                TopologyLabel.Shower, PhysicsLabel.PhotonShower,
                cluster_label
            )


    def process_kaon_plus(self):

        kaonpluses = self.simulation_wrangler.get_trackid_pdg_code(321)
        for kaonplus in kaonpluses:
            kaonplus_daughters = self.simulation_wrangler.trackid_daughters[kaonplus]
            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(kaonplus_daughters, 11)
            other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(kaonplus_daughters, 11)
            elec_hits = self.simulation_wrangler.get_hits_trackid(elec_daughters)
            elec_segments = self.simulation_wrangler.get_segments_trackid(elec_daughters)

            kaonplus_progeny = self.simulation_wrangler.trackid_progeny[kaonplus]
            kaonplus_hits = self.simulation_wrangler.get_hits_trackid(kaonplus)
            kaonplus_segments = self.simulation_wrangler.get_segments_trackid(kaonplus)

            # Set kaonplus detsim labels to Track:kaononMinus
            track_label = self.iterate_topology_label()
            self.set_labels(
                kaonplus_hits, kaonplus_segments, kaonplus,
                TopologyLabel.Track, PhysicsLabel.HIPIonization,
                track_label
            )
            
            self.set_labels(
                elec_hits, elec_segments, elec_daughters,
                TopologyLabel.Track, PhysicsLabel.HIPIonization,
                track_label
            )

            self.process_showers_list(other_daughters, self.iterate_topology_label())
            self.process_showers_list(kaonplus_progeny, self.iterate_topology_label())


    def process_kaon_minus(self):
        kaonmins = self.simulation_wrangler.get_trackid_pdg_code(-321)
        for kaonmin in kaonmins:
            kaonmin_daughters = self.simulation_wrangler.trackid_daughters[kaonmin]
            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(kaonmin_daughters, 11)
            other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(kaonmin_daughters, 11)
            elec_hits = self.simulation_wrangler.get_hits_trackid(elec_daughters)
            elec_segments = self.simulation_wrangler.get_segments_trackid(elec_daughters)

            kaonmin_progeny = self.simulation_wrangler.trackid_progeny[kaonmin]
            kaonmin_hits = self.simulation_wrangler.get_hits_trackid(kaonmin)
            kaonmin_segments = self.simulation_wrangler.get_segments_trackid(kaonmin)

            # Set kaonplus detsim labels to Track:kaononMinus
            track_label = self.iterate_topology_label()
            self.set_labels(
                kaonmin_hits, kaonmin_segments, kaonmin,
                TopologyLabel.Track, PhysicsLabel.HIPIonization,
                track_label
            )
            
            self.set_labels(
                elec_hits, elec_segments, elec_daughters,
                TopologyLabel.Track, PhysicsLabel.HIPIonization,
                track_label
            )

            self.process_showers_list(other_daughters, self.iterate_topology_label())
            self.process_showers_list(kaonmin_progeny, self.iterate_topology_label())

    def process_protons(self):
        protons = self.simulation_wrangler.get_trackid_pdg_code(2212)
        for proton in protons:
            proton_daughters = self.simulation_wrangler.trackid_daughters[proton]
            elec_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(proton_daughters, 11)
            other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(proton_daughters, 11)

            delta_daughters = self.simulation_wrangler.filter_trackid_process(elec_daughters, PhysicsLabel.HadronIonization)
            other_elec_daughters = self.simulation_wrangler.filter_trackid_not_process(elec_daughters, PhysicsLabel.HadronIonization)

            delta_hits = self.simulation_wrangler.get_hits_trackid(delta_daughters)
            other_elec_hits = self.simulation_wrangler.get_hits_trackid(other_elec_daughters)
            delta_segments = self.simulation_wrangler.get_segments_trackid(delta_daughters)
            other_elec_segments = self.simulation_wrangler.get_segments_trackid(other_elec_daughters)

            proton_progeny = self.simulation_wrangler.trackid_progeny[proton]
            proton_hits = self.simulation_wrangler.get_hits_trackid(proton)
            proton_segments = self.simulation_wrangler.get_segments_trackid(proton)

            # Set proton detsim labels to Track:HIPIonization
            track_label = self.iterate_topology_label()
            self.set_labels(
                proton_hits, proton_segments, proton,
                TopologyLabel.Track, PhysicsLabel.HIPIonization,
                track_label
            )
            
            self.set_labels(
                delta_hits, delta_segments, delta_daughters,
                TopologyLabel.Track, PhysicsLabel.DeltaElectron,
                track_label
            )
            
            self.process_showers_list(other_elec_daughters, self.iterate_topology_label())
            self.process_showers_list(other_daughters, self.iterate_topology_label())
            self.process_showers_list(proton_progeny, self.iterate_topology_label())


    def process_neutron_captures(self):
        neutrons = self.simulation_wrangler.get_trackid_pdg_code(2112)
        neutron_daughters = self.simulation_wrangler.get_daughter_trackid_trackid(neutrons)
        gamma_daughters = self.simulation_wrangler.filter_trackid_abs_pdg_code(neutron_daughters, 22)
        other_daughters = self.simulation_wrangler.filter_trackid_not_abs_pdg_code(neutron_daughters, 22)
        capture_daughters = self.simulation_wrangler.filter_trackid_process(gamma_daughters, PhysicsLabel.NeutronCapture)
        other_gammas = self.simulation_wrangler.filter_trackid_not_process(gamma_daughters, PhysicsLabel.NeutronCapture)

        for capture in capture_daughters:
            for gamma in capture:
                gamma_energy = self.simulation_wrangler.get_total_hit_energy(gamma, 5)
                gamma_hits = self.simulation_wrangler.get_hits_trackid(gamma)
                gamma_segments = self.simulation_wrangler.get_segments_trackid(gamma)

                particle_label = PhysicsLabel.NeutronCaptureGammaOther
                if gamma_energy in {0.00474, 0.00475}:
                    particle_label = PhysicsLabel.NeutronCaptureGamma474
                elif gamma_energy in {0.00336, 0.00337}:
                    particle_label = PhysicsLabel.NeutronCaptureGamma336
                elif gamma_energy in {0.00256, 0.00257}:
                    particle_label = PhysicsLabel.NeutronCaptureGamma256
                elif gamma_energy in {0.00118, 0.00119}:
                    particle_label = PhysicsLabel.NeutronCaptureGamma118
                elif gamma_energy in {0.00083, 0.00084}:
                    particle_label = PhysicsLabel.NeutronCaptureGamma083
                elif gamma_energy in {0.00051, 0.00052}:
                    particle_label = PhysicsLabel.NeutronCaptureGamma051
                elif gamma_energy in {0.00016, 0.00017}:
                    particle_label = PhysicsLabel.NeutronCaptureGamma016

                cluster_label = self.iterate_topology_label()
                self.set_labels(
                    gamma_hits, gamma_segments, gamma,
                    TopologyLabel.Blip, particle_label,
                    cluster_label
                )
                
        self.process_showers(other_daughters)
        self.process_showers(other_gammas)


    def process_nuclear_recoils(self):
        ar41 = self.simulation_wrangler.get_trackid_pdg_code(1000180410)
        ar40 = self.simulation_wrangler.get_trackid_pdg_code(1000180400)
        ar39 = self.simulation_wrangler.get_trackid_pdg_code(1000180390)
        ar38 = self.simulation_wrangler.get_trackid_pdg_code(1000180380)
        ar37 = self.simulation_wrangler.get_trackid_pdg_code(1000180370)
        ar36 = self.simulation_wrangler.get_trackid_pdg_code(1000180360)
        
        ar41_daughters = self.simulation_wrangler.get_daughters_trackid(ar41)
        ar40_daughters = self.simulation_wrangler.get_daughters_trackid(ar40)
        ar39_daughters = self.simulation_wrangler.get_daughters_trackid(ar39)
        ar38_daughters = self.simulation_wrangler.get_daughters_trackid(ar38)
        ar37_daughters = self.simulation_wrangler.get_daughters_trackid(ar37)
        ar36_daughters = self.simulation_wrangler.get_daughters_trackid(ar36)

        s33 = self.simulation_wrangler.get_trackid_pdg_code(1000160330)
        s33_daughters = self.simulation_wrangler.get_daughters_trackid(s33)
        s35 = self.simulation_wrangler.get_trackid_pdg_code(1000160350)
        s35_daughters = self.simulation_wrangler.get_daughters_trackid(s35)
        s36 = self.simulation_wrangler.get_trackid_pdg_code(1000160360)
        s36_daughters = self.simulation_wrangler.get_daughters_trackid(s36)

        cl36 = self.simulation_wrangler.get_trackid_pdg_code(1000170360)
        cl36_daughters = self.simulation_wrangler.get_daughters_trackid(cl36)
        cl37 = self.simulation_wrangler.get_trackid_pdg_code(1000170370)
        cl37_daughters = self.simulation_wrangler.get_daughters_trackid(cl37)
        cl39 = self.simulation_wrangler.get_trackid_pdg_code(1000170390)
        cl39_daughters = self.simulation_wrangler.get_daughters_trackid(cl39)
        cl40 = self.simulation_wrangler.get_trackid_pdg_code(1000170400)
        cl40_daughters = self.simulation_wrangler.get_daughters_trackid(cl40)

        for ar in ar41:
            ar41_hits = self.simulation_wrangler.get_hits_trackid(ar)
            ar41_segments = self.simulation_wrangler.get_segments_trackid(ar)
            self.set_labels(
                ar41_hits, ar41_segments, ar,
                TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
                self.iterate_topology_label()
            )

        for ar in ar40:
            ar40_hits = self.simulation_wrangler.get_hits_trackid(ar)
            ar40_segments = self.simulation_wrangler.get_segments_trackid(ar)
            self.set_labels(
                ar40_hits, ar40_segments, ar,
                TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
                self.iterate_topology_label()
            )

        for ar in ar39:
            ar39_hits = self.simulation_wrangler.get_hits_trackid(ar)
            ar39_segments = self.simulation_wrangler.get_segments_trackid(ar)
            self.set_labels(
                ar39_hits, ar39_segments, ar,
                TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
                self.iterate_topology_label()
            )

        for ar in ar38:
            ar38_hits = self.simulation_wrangler.get_hits_trackid(ar)
            ar38_segments = self.simulation_wrangler.get_segments_trackid(ar)
            self.set_labels(
                ar38_hits, ar38_segments, ar,
                TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
                self.iterate_topology_label()
            )

        for ar in ar37:
            ar37_hits = self.simulation_wrangler.get_hits_trackid(ar)
            ar37_segments = self.simulation_wrangler.get_segments_trackid(ar)
            self.set_labels(
                ar37_hits, ar37_segments, ar,
                TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
                self.iterate_topology_label()
            )

        for ar in ar36:
            ar36_hits = self.simulation_wrangler.get_hits_trackid(ar)
            ar36_segments = self.simulation_wrangler.get_segments_trackid(ar)
            self.set_labels(
                ar36_hits, ar36_segments, ar,
                TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
                self.iterate_topology_label()
            )

        for s in s33:
            s33_hits = self.simulation_wrangler.get_hits_trackid(s)
            s33_segments = self.simulation_wrangler.get_segments_trackid(s)
            self.set_labels(
                s33_hits, s33_segments, s,
                TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
                self.iterate_topology_label()
            )

        for s in s35:
            s35_hits = self.simulation_wrangler.get_hits_trackid(s)
            s35_segments = self.simulation_wrangler.get_segments_trackid(s)
            self.set_labels(
                s35_hits, s35_segments, s,
                TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
                self.iterate_topology_label()
            )

        for s in s36:
            s36_hits = self.simulation_wrangler.get_hits_trackid(s)
            s36_segments = self.simulation_wrangler.get_segments_trackid(s)
            self.set_labels(
                s36_hits, s36_segments, s,
                TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
                self.iterate_topology_label()
            )
        
        for cl in cl36:
            cl36_hits = self.simulation_wrangler.get_hits_trackid(cl)
            cl36_segments = self.simulation_wrangler.get_segments_trackid(cl)
            self.set_labels(
                cl36_hits, cl36_segments, cl,
                TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
                self.iterate_topology_label()
            )

        for cl in cl37:
            cl37_hits = self.simulation_wrangler.get_hits_trackid(cl)
            cl37_segments = self.simulation_wrangler.get_segments_trackid(cl)
            self.set_labels(
                cl37_hits, cl37_segments, cl,
                TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
                self.iterate_topology_label()
            )

        for cl in cl39:
            cl39_hits = self.simulation_wrangler.get_hits_trackid(cl)
            cl39_segments = self.simulation_wrangler.get_segments_trackid(cl)
            self.set_labels(
                cl39_hits, cl39_segments, cl,
                TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
                self.iterate_topology_label()
            )

        for cl in cl40:
            cl40_hits = self.simulation_wrangler.get_hits_trackid(cl)
            cl40_segments = self.simulation_wrangler.get_segments_trackid(cl)
            self.set_labels(
                cl40_hits, cl40_segments, cl,
                TopologyLabel.Blip, PhysicsLabel.NuclearRecoil,
                self.iterate_topology_label()
            )


    def process_electron_recoils(self):
        deuterons = self.simulation_wrangler.get_trackid_pdg_code(1000010020)
        tritons = self.simulation_wrangler.get_trackid_pdg_code(1000010030)
        alphas = self.simulation_wrangler.get_trackid_pdg_code(1000020040)
        
        inelastic_alphas = self.simulation_wrangler.filter_trackid_process(alphas, PhysicsLabel.NeutronInelastic)
        
        deuteron_daughters = self.simulation_wrangler.get_daughters_trackid(deuterons)
        triton_daughters = self.simulation_wrangler.get_daughters_trackid(tritons)
        inelastic_alpha_daughters = self.simulation_wrangler.get_daughters_trackid(inelastic_alphas)
        
        for deuteron in deuterons:
            deuteron_hits = self.simulation_wrangler.get_hits_trackid(deuteron)
            deuteron_segments = self.simulation_wrangler.get_segments_trackid(deuteron)
            
            self.set_labels(
                deuteron_hits, deuteron_segments, deuteron,
                TopologyLabel.Blip, PhysicsLabel.ElectronRecoil,
                self.iterate_topology_label()
            )
                   
        for triton in tritons:
            triton_hits = self.simulation_wrangler.get_hits_trackid(triton)
            triton_segments = self.simulation_wrangler.get_segments_trackid(triton)
            
            self.set_labels(
                triton_hits, triton_segments, triton,
                TopologyLabel.Blip, PhysicsLabel.ElectronRecoil,
                self.iterate_topology_label()
            )
        
        for inelastic_alpha in inelastic_alphas:
            inelastic_alpha_hits = self.simulation_wrangler.get_hits_trackid(inelastic_alpha)
            inelastic_alpha_segments = self.simulation_wrangler.get_segments_trackid(inelastic_alpha)
            
            self.set_labels(
                inelastic_alpha_hits, inelastic_alpha_segments, inelastic_alpha,
                TopologyLabel.Blip, PhysicsLabel.ElectronRecoil,
                self.iterate_topology_label()
            )
            

    def process_ar39(self):
        pass

    def process_ar42(self):
        pass

    def process_kr85(self):
        pass

    def process_rn222(self):
        pass

    def process_cosmics(self):
        pass
