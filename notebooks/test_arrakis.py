import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

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
