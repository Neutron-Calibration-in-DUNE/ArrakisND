"""
"""
from enum import Enum


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


class PhysicsMicro(Enum):
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


class PhysicsMacro(Enum):
    """
    Macro-lebel decsriptions of topological types.
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
