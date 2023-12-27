"""
ArrakisND data
"""

wrangler_modes = ["map", "numpy"]

process_type = {
    "not_defined":  0,
    "transportation":   1,
    "electromagentic":  2,
    "optical":          3,
    "hadronic":         4,
    "photo_lepton_hadron":  5,
    "decay":            6,
    "general":          7,
    "parameterization": 8,
    "user_defined":     9,
}

process_subtype = {
    "coulomb_scattering":   1,
    "em_ionization":        2,
    "em_bremsstrahlung":    3,
    "em_pair_production_by_charge":   4,
    "em_nuclear_stopping":      8,
    "em_multiple_scattering":   10,
    "em_photoelectric":         12,
    "em_compton_scattering":    13,
    "em_gamma_conversion":      14,

    "hadron_elastic":   111,
    "hadron_inelastic": 121,
    "hadron_capture":   131,
    "hadron_charge_exchange":   161
}
