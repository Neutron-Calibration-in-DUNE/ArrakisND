"""
"""
from enum import Enum


class ProcessType(Enum):
    """
    Geant4 process identifications.
    """
    Undefined = -1
    NotDefined = 0
    Transportation = 1
    Electromagnetic = 2
    Optical = 3
    Hadronic = 4
    PhotoLeptonHadron = 5
    Decay = 6
    General = 7
    Parameterization = 8
    UserDefined = 9
    Parallel = 10
    Phonon = 11
    UCN = 12


class SubProcessType(Enum):
    """
    Geant4 subprocess identifications.
    """
    Undefined = -1
    Primary = 0
    CoulombScattering = 1
    Ionization = 2
    Bremsstrahlung = 3
    PairProdByCharge = 4
    Annihilation = 5
    AnnihilationToMuMu = 6
    AnnihilationToHadrons = 7
    NuclearStopping = 8
    MultipleScattering = 9
    Msc = 10
    Rayleigh = 11
    PhotoElectricEffect = 12
    ComptonScattering = 13
    GammaConversion = 14
    GammaConversionToMuMu = 15
    Cerenkov = 21
    Scintillation = 22
    SynchrotronRadiation = 23
    TransitionRadiation = 24
    OpticalAbsorption = 31
    OpticalBoundary = 32
    OpticalRayleigh = 33
    OpticalWLS = 34
    OpticalMieHG = 35
    UCNLoss = 41
    UCNAbsorption = 42
    UCNBoundary = 43
    UCNMultiScattering = 44
    LowEnergyElastic = 51
    LowEnergyExcitation = 52
    LowEnergyIonization = 53
    LowEnergyVibrationalExcitation = 54
    LowEnergyAttachment = 55
    LowEnergyChargeDecrease = 56
    LowEnergyChargeIncrease = 57
    LowEnergyElectronSolvation = 58
    LowEnergyMolecularDecay = 59
    LowEnergyTransportation = 60
    LowEnergyBrownianTransportation = 61
    LowEnergyDoubleIonization = 62
    LowEnergyDoubleCap = 63
    LowEnergyIoniTransfer = 64
    LowEnergyStaticMol = 65
    LowEnergyScavenger = 66
    Transportation = 91
    CoupledTransportation = 92
    HadronElastic = 111
    NeutronGeneral = 116
    HadronInelastic = 121
    HadronCapture = 131
    MuAtomicCapture = 132
    HadronFission = 141
    HadronCaptureAtRest = 151
    LeptonCaptureAtRest = 152
    ChargeExchange = 161
    NuOscillation = 165
    NuElectron = 166
    NuNucleus = 167
    Decay = 201
    DecayRadioactive = 210
    DecayUnknown = 211
    DecayMuAtom = 221
    DecayExternal = 231
    FastSimManagerProcess = 301
    EMDissociation = 310
    StepLimiter = 401
    UserSpecialCuts = 402
    NeutronKiller = 403
    ParallelWorld = 491
    DNAUnknownModel = 11000
    Ritchie1994eSolvation = 11001
    Terrisol1990eSolvation = 11002
    Meesungnoen2002eSolvation = 11003
    Kreipl2009eSolvation = 11004
    Meesungnoensolid2002eSolvation = 11005


class Topology(Enum):
    """
    Topological decsriptions of energy
    deposits.
    """
    Undefined = -1
    Track = 0
    Shower = 1
    Blip = 2


class Physics(Enum):
    """
    micro-level descriptions of topological types.
    This further breaks down the blip/track/shower
    topology labels into micro-level physics.
    """
    Undefined = -1
    MIP = 0
    HIP = 1
    ElectronIonization = 2
    DeltaElectron = 3
    MichelElectron = 4
    GammaCompton = 5
    GammaConversion = 6
    NuclearRecoil = 7
    ElectronRecoil = 8


class Tracklette(Enum):
    """
    Macro-level descriptions of tracklette types

    Args:
        Enum (_type_): _description_
    """
    Undefined = -1
    MIP = 0
    HIP = 1
    Delta = 2
    Michel = 3


class Track(Enum):
    """
    Macro-level descriptions of track types

    Args:
        Enum (_type_): _description_
    """
    Undefined = -1
    MIP = 0
    HIP = 1
    Delta = 2
    Michel = 3


class Fragment(Enum):
    """
    Macro-level descriptions of fragment types

    Args:
        Enum (_type_): _description_
    """
    Undefined = -1
    GammaCompton = 0
    GammaConversion = 1
    PairProductionByCharge = 2
    Bremsstrahlung = 3
    PhotoelectricEffect = 4
    ElectronIonization = 5
    Annihilation = 6


class Shower(Enum):
    """
    Macro-level descriptions of shower types
    Args:
        Enum (_type_): _description_
    """
    Undefined = -1
    Electromagnetic = 0
    Pion0Decay = 1
    PiPlusDecay = 2
    PiMinusDecay = 3
    Kaon0Decay = 4
    KaonShortDecay = 5
    KaonLongDecay = 6
    KaonPlusDecay = 7
    KaonMinusDecay = 8
    D0Decay = 9
    DPlusDecay = 10
    DMinusDecay = 11
    LambdaDecay = 12
    Sigma0Decay = 13
    SigmaPlusDecay = 14
    SigmaMinusDecay = 15


class Blip(Enum):
    """
    Macro-level descriptions of blip types.
    """
    Undefined = -1
    GammaCompton = 0
    GammaConversion = 1
    PairProductionByCharge = 2
    Bremsstrahlung = 3
    PhotoelectricEffect = 4
    ElectronIonization = 5
    Annihilation = 6
    NeutronCaptureGamma474 = 7
    NeutronCaptureGamma336 = 8
    NeutronCaptureGamma256 = 9
    NeutronCaptureGamma118 = 10
    NeutronCaptureGamma083 = 11
    NeutronCaptureGamma051 = 12
    NeutronCaptureGamma016 = 13
    NeutronCaptureGammaOther = 14
    AlphaDecay = 15
    BetaDecay = 16
    GammaDecay = 17
    NuclearRecoil = 18
    ElectronRecoil = 19


class Interaction(Enum):
    """
    Macro-level descriptions of final state interaction types
    Args:
        Enum (_type_): _description_
    """
    Undefined = -1
    Noise = 0


class Neutrino(Enum):
    """
    Macro-level descriptions of neutrino types.
    """
    Undefined = -1
    NCElectronNeutrino = 0
    CCElectronNeutrino = 1
    NCAntiElectronNeutrino = 2
    CCAntiElectronNeutrino = 3
    NCMuonNeutrino = 4
    CCMuonNeutrino = 5
    NCAntiMuonNeutrino = 6
    CCAntiMuonNeutrino = 7
    NCTauonNeutrino = 8
    CCTauonNeutrino = 9
    NCAntiTauonNeutrino = 10
    CCAntiTauonNeutrino = 11


process_type_dict = {
    -1: 'undefined',
    0: 'not_defined',
    1: 'transportation',
    2: 'electromagnetic',
    3: 'optical',
    4: 'hadronic',
    5: 'photo_lepton_hadron',
    6: 'decay',
    7: 'general',
    8: 'parameterization',
    9: 'user_defined',
    10: 'parallel',
    11: 'phonon',
    12: 'ucn'
}

sub_process_type_dict = {
    -1: 'undefined',
    0: 'primary',
    1: 'coulomb_scattering',
    2: 'ionization',
    3: 'bremsstrahlung',
    4: 'pair_prod_by_charge',
    5: 'annihilation',
    6: 'annihilation_to_mu_mu',
    7: 'annihilation_to_hadrons',
    8: 'nuclear_stopping',
    9: 'multiple_scattering',
    10: 'msc',
    11: 'rayleigh',
    12: 'photo_electric_effect',
    13: 'compton_scattering',
    14: 'gamma_conversion',
    15: 'gamma_conversion_to_mu_mu',
    21: 'cerenkov',
    22: 'scintillation',
    23: 'synchrotron_radiation',
    24: 'transition_radiation',
    31: 'optical_absorption',
    32: 'optical_boundary',
    33: 'optical_rayleigh',
    34: 'optical_wls',
    35: 'optical_mie_h_g',
    41: 'ucn_loss',
    42: 'ucn_absorption',
    43: 'ucn_boundary',
    44: 'ucn_multi_scattering',
    51: 'low_energy_elastic',
    52: 'low_energy_excitation',
    53: 'low_energy_ionization',
    54: 'low_energy_vibrational_excitation',
    55: 'low_energy_attachment',
    56: 'low_energy_charge_decrease',
    57: 'low_energy_charge_increase',
    58: 'low_energy_electron_solvation',
    59: 'low_energy_molecular_decay',
    60: 'low_energy_transportation',
    61: 'low_energy_brownian_transportation',
    62: 'low_energy_double_ionization',
    63: 'low_energy_double_cap',
    64: 'low_energy_ioni_transfer',
    65: 'low_energy_static_mol',
    66: 'low_energy_scavenger',
    91: 'transportation',
    92: 'coupled_transportation',
    111: 'hadron_elastic',
    116: 'neutron_general',
    121: 'hadron_inelastic',
    131: 'hadron_capture',
    132: 'mu_atomic_capture',
    141: 'hadron_fission',
    151: 'hadron_capture_at_rest',
    152: 'lepton_capture_at_rest',
    161: 'charge_exchange',
    165: 'nu_oscillation',
    166: 'nu_electron',
    167: 'nu_nucleus',
    201: 'decay',
    210: 'decay_radioactive',
    211: 'decay_unknown',
    221: 'decay_mu_atom',
    231: 'decay_external',
    301: 'fast_sim_manager_process',
    310: 'em_dissociation',
    401: 'step_limiter',
    402: 'user_special_cuts',
    403: 'neutron_killer',
    491: 'parallel_world',
    11000: 'dna_unknown_model',
    11001: 'ritchie1994e_solvation',
    11002: 'terrisol1990e_solvation',
    11003: 'meesungnoen2002e_solvation',
    11004: 'kreipl2009e_solvation',
    11005: 'meesungnoen_solid2002e_solvation'
}

reverse_process_type_dict = {
    'undefined': -1,
    'not_defined': 0,
    'transportation': 1,
    'electromagnetic': 2,
    'optical': 3,
    'hadronic': 4,
    'photo_lepton_hadron': 5,
    'decay': 6,
    'general': 7,
    'parameterization': 8,
    'user_defined': 9,
    'parallel': 10,
    'phonon': 11,
    'ucn': 12
}

reverse_sub_process_type_dict = {
    'undefined': -1,
    'primary': 0,
    'coulomb_scattering': 1,
    'ionization': 2,
    'bremsstrahlung': 3,
    'pair_prod_by_charge': 4,
    'annihilation': 5,
    'annihilation_to_mu_mu': 6,
    'annihilation_to_hadrons': 7,
    'nuclear_stopping': 8,
    'multiple_scattering': 9,
    'msc': 10,
    'rayleigh': 11,
    'photo_electric_effect': 12,
    'compton_scattering': 13,
    'gamma_conversion': 14,
    'gamma_conversion_to_mu_mu': 15,
    'cerenkov': 21,
    'scintillation': 22,
    'synchrotron_radiation': 23,
    'transition_radiation': 24,
    'optical_absorption': 31,
    'optical_boundary': 32,
    'optical_rayleigh': 33,
    'optical_wls': 34,
    'optical_mie_h_g': 35,
    'ucn_loss': 41,
    'ucn_absorption': 42,
    'ucn_boundary': 43,
    'ucn_multi_scattering': 44,
    'low_energy_elastic': 51,
    'low_energy_excitation': 52,
    'low_energy_ionization': 53,
    'low_energy_vibrational_excitation': 54,
    'low_energy_attachment': 55,
    'low_energy_charge_decrease': 56,
    'low_energy_charge_increase': 57,
    'low_energy_electron_solvation': 58,
    'low_energy_molecular_decay': 59,
    'low_energy_transportation': 60,
    'low_energy_brownian_transportation': 61,
    'low_energy_double_ionization': 62,
    'low_energy_double_cap': 63,
    'low_energy_ioni_transfer': 64,
    'low_energy_static_mol': 65,
    'low_energy_scavenger': 66,
    'transportation': 91,
    'coupled_transportation': 92,
    'hadron_elastic': 111,
    'neutron_general': 116,
    'hadron_inelastic': 121,
    'hadron_capture': 131,
    'mu_atomic_capture': 132,
    'hadron_fission': 141,
    'hadron_capture_at_rest': 151,
    'lepton_capture_at_rest': 152,
    'charge_exchange': 161,
    'nu_oscillation': 165,
    'nu_electron': 166,
    'nu_nucleus': 167,
    'decay': 201,
    'decay_radioactive': 210,
    'decay_unknown': 211,
    'decay_mu_atom': 221,
    'decay_external': 231,
    'fast_sim_manager_process': 301,
    'em_dissociation': 310,
    'step_limiter': 401,
    'user_special_cuts': 402,
    'neutron_killer': 403,
    'parallel_world': 491,
    'dna_unknown_model': 11000,
    'ritchie1994e_solvation': 11001,
    'terrisol1990e_solvation': 11002,
    'meesungnoen2002e_solvation': 11003,
    'kreipl2009e_solvation': 11004,
    'meesungnoen_solid2002e_solvation': 11005
}
