import numpy as np
import matplotlib.pyplot as plt


arr = np.load("/home/ncarrara/workspace/MiniRun4_1E19_RHC.flow.00002.FLOW.arrakis_nd.npz", allow_pickle=True)
event = 1
det_features = arr["det_features"][event]
mc_features = arr["mc_features"][event]
classes = arr["classes"][event]
clusters = arr["clusters"][event]

x = det_features[:, 0]
y = det_features[:, 1]
z = det_features[:, 2]
q = det_features[:, 3]

E = mc_features[:, 2]

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
ax = fig.add_subplot(111, projection="3d")
ax.scatter(
    x, y, z, 
    c=topology_labels, cmap='tab20',
    s=q
)
plt.legend()
plt.show()
