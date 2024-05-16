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
