"""
"""
import numpy as np
import h5py

class DetectorPointCloud:
    """
    A point cloud structure for an event, which consists of a three dimensional
    set of points (x,y,z) together with the reconstructed charge (Q).  Extra variables
    from simulation are the kinetic energy (E), the drift time (t_drift), the
    ... (ts_pps), and the various labels.
    """
    def __init__(self):
        self.clear()
    
    def clear(self):
        """
        
        """
        self.x = []
        self.y = []
        self.z = []
        self.t_drift = []
        self.ts_pps = []
        self.Q = []
        self.E = []
        self.segment_ids = []
        self.segment_fractions = []

        self.source_label = []
        self.topology_label = []
        self.particle_label = []
        self.physics_label = []

        self.unique_topology = []
        self.unique_particle = []
        self.unique_physics = []

        self.source_labels = []
        self.topology_labels = []
        self.particle_labels = []
        self.physics_labels = []

        self.unique_topologies = []
        self.unique_particles = []
        self.unique_physicses = []

    def add_points(self, x, y, z, t_drift, ts_pps, Q, E, segment_ids, segment_fractions):
        for i in range(len(x)):
            self.add_point(
                x[i], y[i], z[i], t_drift[i], ts_pps[i], Q[i], E[i],
                segment_ids[i], segment_fractions[i]
            )

    def add_point(self,
            x:      float, 
            y:      float, 
            z:      float, 
            t_drift:    float, 
            ts_pps: float, 
            Q:      float, 
            E:      float, 
            segment_ids:    list=[],
            segment_fractions:  list=[],
        ):

        self.x.extend([x])
        self.y.extend([y])
        self.z.extend([z])
        self.t_drift.extend([t_drift])
        self.ts_pps.extend([ts_pps])
        self.Q.extend([Q])
        self.E.extend([E])
        self.segment_ids.extend([segment_ids])
        self.segment_fractions.extend([segment_fractions])

        self.source_label.extend([-1])
        self.topology_label.extend([-1])
        self.particle_label.extend([-1])
        self.physics_label.extend([-1])

        self.unique_topology.extend([-1])
        self.unique_particle.extend([-1])
        self.unique_physics.extend([-1])

        self.source_labels.extend([-1 for _ in range(len([segment_ids]))])
        self.topology_labels.extend([-1 for _ in range(len([segment_ids]))])
        self.particle_labels.extend([-1 for _ in range(len([segment_ids]))])
        self.physics_labels.extend([-1 for _ in range(len([segment_ids]))])

        self.unique_topologies.extend([-1 for _ in range(len([segment_ids]))])
        self.unique_particles.extend([-1 for _ in range(len([segment_ids]))])
        self.unique_physicses.extend([-1 for _ in range(len([segment_ids]))])



        