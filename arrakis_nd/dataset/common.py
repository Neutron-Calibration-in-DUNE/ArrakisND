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
    Noise = 0
    Track = 1
    Shower = 2
    Blip = 3


class Physics(Enum):
    """
    micro-level descriptions of topological types.
    This further breaks down the blip/track/shower
    topology labels into micro-level physics.
    """
    Undefined = -1
    Noise = 0
    MIP = 1
    HIP = 2
    ElectronIonization = 3
    DeltaElectron = 4
    MichelElectron = 5
    GammaCompton = 6
    GammaConversion = 7
    NuclearRecoil = 8
    ElectronRecoil = 9


class Tracklette(Enum):
    """
    Macro-level descriptions of tracklette types

    Args:
        Enum (_type_): _description_
    """
    Undefined = -1
    Noise = 0
    MIP = 1
    HIP = 2
    Delta = 3
    Michel = 4


class Track(Enum):
    """
    Macro-level descriptions of track types

    Args:
        Enum (_type_): _description_
    """
    Undefined = -1
    Noise = 0
    MIP = 1
    HIP = 2
    Delta = 3
    Michel = 4


class Fragment(Enum):
    """
    Macro-level descriptions of fragment types

    Args:
        Enum (_type_): _description_
    """
    Undefined = -1
    Noise = 0
    GammaCompton = 1
    GammaConversion = 2
    PairProductionByCharge = 3
    Bremsstrahlung = 4
    PhotoelectricEffect = 5
    ElectronIonization = 6
    Annihilation = 7


class Shower(Enum):
    """
    Macro-level descriptions of shower types
    Args:
        Enum (_type_): _description_
    """
    Undefined = -1
    Noise = 0
    Electromagnetic = 1
    Pion0Decay = 2
    PiPlusDecay = 3
    PiMinusDecay = 4
    Kaon0Decay = 5
    KaonShortDecay = 6
    KaonLongDecay = 7
    KaonPlusDecay = 8
    KaonMinusdecay = 9
    D0Decay = 10
    DPlusDecay = 11
    DMinusDecay = 12
    LambdaDecay = 13
    Sigma0Decay = 14
    SigmaPlusDecay = 15
    SigmaMinusDecay = 16


class Blip(Enum):
    """
    Macro-level descriptions of blip types.
    """
    Undefined = -1
    Noise = 0
    GammaCompton = 1
    GammaConversion = 2
    PairProductionByCharge = 3
    Bremsstrahlung = 4
    PhotoelectricEffect = 5
    ElectronIonization = 6
    Annihilation = 7
    NeutronCaptureGamma474 = 8
    NeutronCaptureGamma336 = 9
    NeutronCaptureGamma256 = 10
    NeutronCaptureGamma118 = 11
    NeutronCaptureGamma083 = 12
    NeutronCaptureGamma051 = 13
    NeutronCaptureGamma016 = 14
    NeutronCaptureGammaOther = 15
    AlphaDecay = 16
    BetaDecay = 17
    GammaDecay = 18
    NuclearRecoil = 19
    ElectronRecoil = 20


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
    Noise = 0
    NCElectronNeutrino = 1
    CCElectronNeutrino = 2
    NCAntiElectronNeutrino = 3
    CCAntiElectronNeutrino = 4
    NCMuonNeutrino = 5
    CCMuonNeutrino = 6
    NCAntiMuonNeutrino = 7
    CCAntiMuonNeutrino = 8
    NCTauonNeutrino = 9
    CCTauonNeutrino = 10
    NCAntiTauonNeutrino = 11
    CCAntiTauonNeutrino = 12
