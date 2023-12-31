"""
"""
from enum import Enum


class ParticleLabel(Enum):
    Undefined = -1
    Noise = 0
    Electron = 11
    Positron = -11
    ElectronNeutrino = 12
    AntiElectronNeutrino = -12
    Muon = 13
    AntiMuon = -13
    MuonNeutrino = 14
    AntiMuonNeutrino = -14
    Tauon = 15
    AntiTauon = -15
    TauonNeutrino = 16
    AntiTauonNeutrino = -16
    Gamma = 22
    Pion0 = 111
    PionPlus = 211
    PionMinus = -211
    Kaon0 = 311
    KaonPlus = 321
    KaonMinus = -321
    Neutron = 2112
    AntiNeutron = -2112
    Proton = 2212
    AntiProton = -2212
    Deuteron = 1000010020
    Triton = 1000010030
    Alpha = 1000020040


class TopologyLabel(Enum):
    """
    High-level description of the shape of certain
    event types.
    """

    Undefined = -1
    Noise = 0
    Blip = 1
    Track = 2
    Shower = 3


class PhysicsMicroLabel(Enum):
    """
    micro-level descriptions of topological types.
    This further breaks down the blip/track/shower
    topology labels into micro-level physics.
    """

    Undefined = -1
    Noise = 0
    MIPIonization = 1
    HIPIonization = 2
    ElectronIonization = 3
    Bremsstrahlung = 4
    Annihilation = 5
    PhotoElectric = 6
    GammaCompton = 7
    GammaConversion = 8
    HadronElastic = 9
    HadronInelastic = 10


class PhysicsMesoLabel(Enum):
    """
    meso-level descriptions of topological types.
    """

    Undefined = -1
    Noise = 0
    MIP = 1
    HIP = 2
    DeltaElectron = 3
    MichelElectron = 4
    ElectronShower = 5
    PositronShower = 6
    PhotonShower = 7
    LowEnergyIonization = 8
    NeutronCaptureGamma474 = 9
    NeutronCaptureGamma336 = 10
    NeutronCaptureGamma256 = 11
    NeutronCaptureGamma118 = 12
    NeutronCaptureGamma083 = 13
    NeutronCaptureGamma051 = 14
    NeutronCaptureGamma016 = 15
    NeutronCaptureGammaOther = 16
    Pi0Decay = 17
    AlphaDecay = 18
    BetaDecay = 19
    GammaDecay = 20
    NuclearRecoil = 21
    ElectronRecoil = 22


class PhysicsMacroLabel(Enum):
    """
    macro-level descriptions of topological types
    """

    Undefined = -1
    Noise = 0

    # Neutrino interactions
    CCNue = 1
    CCNuMu = 2
    NC = 3

    Cosmics = 4

    # Radiological interactions
    Ar39 = 5
    Ar42 = 6
    K42 = 7
    Kr85 = 8
    Rn222 = 9
    Po218a = 10
    Po218b = 11
    At218a = 12
    At218b = 13
    Rn218 = 14
    Pb214 = 15
    Bi214a = 16
    Bi214b = 17
    Po214 = 18
    Tl210 = 19
    Pb210a = 20
    Pb210b = 21
    Bi210a = 22
    Bi210b = 23
    Po210 = 24


classification_labels = {
    "particle": {
        -1: "undefined",
        0: "noise",
        11: "electron",
        -11: "positron",
        12: "electron_neutrino",
        -12: "anti-electron_neutrino",
        13: "muon",
        -13: "anti-muon",
        14: "muon_neutrino",
        -14: "anti-muon_neutrino",
        15: "tauon",
        -15: "anti-tauon",
        16: "tauon_neutrino",
        -16: "anti-tauon_neutrino",
        22: "gamma",
        111: "pion0",
        211: "pion_plus",
        -211: "pion_minus",
        311: "kaon0",
        321: "kaon_plus",
        -321: "kaon_minus",
        2112: "neutron",
        -2112: "anti-neutron",
        2212: "proton",
        -2212: "anti-proton",
        1000010020: "deuteron",
        1000010030: "triton",
        1000020040: "alpha",
        10000160330: "sulfur_33",
        10000160340: "sulfur_34",
        10000160350: "sulfur_35",
        10000160360: "sulfur_36",
        10000170360: "chlorine_36",
        10000170370: "chlorine_37",
        10000170380: "chlorine_38",
        10000170390: "chlorine_39",
        10000170400: "chlorine_40",
        10000180360: "argon_36",
        10000180370: "argon_37",
        10000180380: "argon_38",
        10000180390: "argon_39",
        10000180400: "argon_40",
        10000180410: "argon_41",
    },
    "topology": {
        -1: "undefined",
        0: "noise",
        1: "blip",
        2: "track",
        3: "shower",
    },
    "physics_micro": {
        -1: "undefined",
        0: "noise",
        1: "mip_ionization",
        2: "hip_ionization",
        3: "electron_ionization",
        4: "bremsstrahlung",
        5: "annihilation",
        6: "photo_electric",
        7: "gamma_compton",
        8: "gamma_conversion",
        9: "hadron_elastic",
        10: "hadron_inelastic",
    },
    "physics_meso": {
        -1: "undefined",
        0: "noise",
        1: "mip",
        2: "hip",
        3: "delta_electron",
        4: "michel_electron",
        5: "electron_shower",
        6: "positron_shower",
        7: "photon_shower",
        8: "low_energy_ionization",
        9: "neutron_capture_gamma_474",
        10: "neutron_capture_gamma_336",
        11: "neutron_capture_gamma_256",
        12: "neutron_capture_gamma_118",
        13: "neutron_capture_gamma_083",
        14: "neutron_capture_gamma_051",
        15: "neutron_capture_gamma_016",
        16: "neutron_capture_gamma_other",
        17: "pi0_decay",
        18: "alpha_decay",
        19: "beta_decay",
        20: "gamma_decay",
        21: "nuclear_recoil",
        22: "electron_recoil"
    },
    "physics_macro": {
        -1: "undefined",
        0: "noise",
        1: "cc_nu_e",
        2: "cc_nu_mu",
        3: "nc",
        4: "cosmics",
        5: "ar39",
        6: "ar42",
        7: "k42",
        8: "kr85",
        9: "rn222",
        10: "po218a",
        11: "po218b",
        12: "at218a",
        13: "at218b",
        14: "rn218",
        15: "pb214",
        16: "bi214a",
        17: "bi214b",
        18: "po214",
        19: "tl210",
        20: "pb210a",
        21: "pb210b",
        22: "bi210a",
        23: "bi210b",
        24: "po210",
    },
    "hit": {
        0: "induction",
        1: "hit",
    },
}
