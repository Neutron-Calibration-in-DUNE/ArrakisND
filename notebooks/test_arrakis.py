import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

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


arr = np.load("/home/ncarrara/workspace/MiniRun4_1E19_RHC.flow.00002.FLOW.arrakis_nd.npz", allow_pickle=True)
event = 2
det_features = arr["det_features"][event]
mc_features = arr["mc_features"][event]
classes = arr["classes"][event]
clusters = arr["clusters"][event]

x = det_features[:, 0]
y = det_features[:, 1]
z = det_features[:, 2]
q = det_features[:, 3]

E = mc_features[:, 2]

num_photons = mc_features[:, 3]
num_photons /= np.mean(num_photons)

particle_labels = classes[:, 0]
topology_labels = classes[:, 1]
physics_micro_labels = classes[:, 2]
physics_meso_labels = classes[:, 3]
physics_macro_labels = classes[:, 4]

unique_particle_labels = clusters[:, 0]
unique_topology_labels = clusters[:, 1]
unique_physics_micro_labels = clusters[:, 2]
unique_physics_meso_labels = clusters[:, 3]
unique_physics_macro_labels = clusters[:, 4]

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
scatter = {}

for unique_label in np.unique(physics_meso_labels):
    mask = (physics_meso_labels == unique_label)
    if unique_label not in classification_labels['physics_meso']:
        continue
    scatter[unique_label] = ax.scatter(
        z[mask], x[mask], y[mask],
        label=classification_labels['physics_meso'][unique_label],
        s=num_photons[mask]*20
    )
ax.set_xlabel("z [mm]")
ax.set_ylabel("x [mm]")

plt.legend()


def update(frame):
    ax.view_init(elev=10, azim=frame)  # Change the azimuthal angle for rotation
    return scatter.values(),


# Create the animation
animation = FuncAnimation(fig, update, frames=np.arange(0, 360, 1), interval=50)

# Save the animation as a GIF
animation.save('rotating_plot.gif', writer='pillow')

plt.show()
