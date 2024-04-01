"""
"""
import h5py
import numpy as np
from scipy.interpolate import splprep, splev
import warnings
from scipy.integrate import quad, IntegrationWarning
import fastcluster

from arrakis_nd.utils.utils import profiler, integrand
from arrakis_nd.plugins.plugin import Plugin
from arrakis_nd.dataset.common import (
    ProcessType, SubProcessType,
    Topology, Physics
)


class ProtonPlugin(Plugin):
    """

    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        """
        super(ProtonPlugin, self).__init__(config)

        self.input_products = [
            'daughters',
            'track_id_hit_map',
            'track_id_hit_segment_map',
            'track_id_hit_t0_map'
        ]
        self.output_products = []
    
    @profiler
    def fit_track(
        self,
        track_t0s,
        track_xyz
    ):
        combined = sorted(zip(track_t0s, track_xyz), key=lambda x: x[0])
        sorted_t0, sorted_xyz = zip(*combined)
        sorted_t0 = np.array(sorted_t0)
        sorted_xyz = np.array(sorted_xyz)

        try:
            """Try to fit a spline curve to the xyz data"""
            tck, u = splprep(
                [sorted_xyz[:, 0], sorted_xyz[:, 1], sorted_xyz[:, 2]],
                s=0
            )

            """Try to integrate this curve"""
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", IntegrationWarning)
                try:
                    curve_length, _ = quad(integrand, 0, 1, args=(tck,))
                except Exception:
                    curve_length = 0

            """Get the derivatives at the beginning and ending points"""
            dxdt_start, dydt_start, dzdt_start = splev(0, tck, der=1)
            dxdt_end, dydt_end, dzdt_end = splev(1, tck, der=1)
            dx_start = np.array([dxdt_start, dydt_start, dzdt_start])
            dx_end = np.array([dxdt_end, dydt_end, dzdt_end])
            dx_start_magnitude = np.linalg.norm(dx_start)
            dx_end_magnitude = np.linalg.norm(dx_end)
            
            if dx_start_magnitude > 0:
                track_dir = [dxdt_start, dydt_start, dzdt_start] / dx_start_magnitude
            else:
                track_dir = [dxdt_start, dydt_start, dzdt_start]
            
            if dx_end_magnitude > 0:
                track_enddir = [dxdt_end, dydt_end, dzdt_end] / dx_end_magnitude
            else:
                track_enddir = [dxdt_end, dydt_end, dzdt_end]

            """Fill the values"""
            track_data = {
                'track_dir': track_dir,
                'track_enddir': track_enddir,
                'track_len_gcm2': 0,
                'track_len_cm': curve_length
            }
        except Exception:
            """Otherwise, these values are undefined"""
            track_data = {
                'track_dir': [0, 0, 0],
                'track_enddir': [0, 0, 0],
                'track_len_gcm2': 0,
                'track_len_cm': 0
            }
        return track_data

    def tpc_bounds(self, i):
        active_tpc_widths=[30.6, 130., 64.] # cm
        tpcs_relative_to_module=[[-15.7,0.,0.],[15.7,0.,0.]]
        modules_relative_to_2x2=[[-33.5,0.,-33.5],
                                [33.5,0.,-33.5],
                                [-33.5,0.,33.5],
                                [33.5,0.,33.5]]
        detector_center=[0.,52.25,0.]
        tpc_bounds=np.array([-active_tpc_widths[i]/2., active_tpc_widths[i]/2.])
        tpc_bounds_relative_to_2x2=[]
        for tpc in tpcs_relative_to_module:
            tpc_bound_relative_to_module = tpc_bounds + tpc[i]
            for module in modules_relative_to_2x2:
                bound = tpc_bound_relative_to_module + module[i]
                tpc_bounds_relative_to_2x2.append(bound)
        bounds_relative_to_NDhall = np.array(tpc_bounds_relative_to_2x2) + detector_center[i]
        return np.unique(bounds_relative_to_NDhall, axis=0)

    def fiducialized_vertex(self, vert_pos):
        lar_x=[(-63.9273,-3.0652),(3.0652,63.9273)]
        lar_y=[(-62.055,62.055)]
        lar_z=[(-64.51125,-2.48125),(2.48125,64.51125)]
        flagX, flagY, flagZ=False, False, False
        
        for i in lar_z:
            if vert_pos[2]>i[0] and vert_pos[2]<i[1]: flagZ=True
        if flagZ==False: return False

        for i in lar_x:
            if vert_pos[0]>i[0] and vert_pos[0]<i[1]: flagX=True
        if flagX==False: return False
        
        for i in lar_y:
            if vert_pos[1]>i[0] and vert_pos[1]<i[1]: flagY=True
        if flagY==False: return False
        
        return True

    @profiler
    def process_event(
        self,
        event: int,
        flow_file: h5py.File,
        arrakis_file: h5py.File,
        event_indices: dict,
        event_products: dict,
    ):
        """
        """
        trajectories = flow_file['mc_truth/trajectories/data'][event_indices['trajectories']]
        interactions = flow_file['mc_truth/interactions/data'][event_indices['interactions']]
        charge = flow_file['charge/calib_final_hits/data'][event_indices['charge']]
        charge_backtrack = flow_file['mc_truth/calib_final_hit_backtrack/data'][event_indices['charge']]
        track_id_hit_map = event_products['track_id_hit_map']
        track_id_hit_segment_map = event_products['track_id_hit_segment_map']
        track_id_hit_t0_map = event_products['track_id_hit_t0_map']
        charge_back_track_segments = charge_backtrack['segment_id']

        trajectories_traj_ids = trajectories['traj_id']
        trajectories_vertex_ids = trajectories['vertex_id']
        trajectories_parent_ids = trajectories['parent_id']
        trajectories_pdg_ids = trajectories['pdg_id']
        trajectories_xyz_start = trajectories['xyz_start']
        trajectories_xyz_end = trajectories['xyz_end']
        trajectories_pxyz_start = trajectories['pxyz_start']
        trajectories_pxyz_end = trajectories['pxyz_end']
        trajectories_E = trajectories['E_start']
        trajectories_t_start = trajectories['t_start']
        trajectories_start_process = trajectories['start_process']
        trajectories_start_subprocess = trajectories['start_subprocess']
        charge_x = charge['x']
        charge_y = charge['y']
        charge_z = charge['z']
        charge_E = charge['E']
        charge_io_group = charge['io_group']

        """Grab the mips/hips"""
        particle_mask = (
            (abs(trajectories_pdg_ids) == 13) |
            (abs(trajectories_pdg_ids) == 15) |
            (abs(trajectories_pdg_ids) == 211) |
            (abs(trajectories_pdg_ids) == 321) |
            (trajectories_pdg_ids == 2212)
        )
        proton_mask = (
            (trajectories_pdg_ids == 2212) &
            (trajectories_start_subprocess == 121)
        )
        particle_hits = []
        for (particle_id, vertex_id) in zip(
            trajectories_traj_ids[particle_mask],
            trajectories_vertex_ids[particle_mask]
        ):
            particle_hits.append(np.array(track_id_hit_map[(particle_id, vertex_id)], dtype=int))
        
        particle_hits = np.concatenate(particle_hits, dtype=int)
        if not any(particle_hits):
            return
        
        particle_charge_xyz = np.array([
            charge_x[np.unique(particle_hits)],
            charge_y[np.unique(particle_hits)],
            charge_z[np.unique(particle_hits)]
        ]).T
           

        """Iterate over the individual (particle_id, vertex_id, traj_index)"""
        for ii, (proton_id, vertex_id, traj_index) in enumerate(zip(
            trajectories_traj_ids[proton_mask],
            trajectories_vertex_ids[proton_mask],
            event_indices['trajectories'][proton_mask]
        )):
            """Get indices in charge_segments arrays"""
            proton_hits = track_id_hit_map[(proton_id, vertex_id)]

            """Check if mip has no hits, go to the next proton"""
            if not any(proton_hits):
                continue

            proton_segments = track_id_hit_segment_map[(proton_id, vertex_id)]
            proton_hit_segments = (
                proton_hits, proton_segments
            )

            """Get the associated t0 values"""
            proton_hit_t0s = track_id_hit_t0_map[(proton_id, vertex_id)]
            proton_t_start = trajectories_t_start[proton_mask][ii]
            
            """Get neutrino vertex and check whether it's fiducialized"""
            vert_mask = interactions['vertex_id'] == vertex_id
            nu_vert = interactions[vert_mask]
            vert_loc = nu_vert['vertex'][0]
            neutrino_vertex_fiducialized = self.fiducialized_vertex(nu_vert['vertex'][0])
                        
            """Get proton parent info"""
            parent_pdg = 0
            parent_E = 0.0
            parent_distance = 0.0
            parent_pxyz_start = [[0, 0, 0]]
            parent_pxyz_end = [[0, 0, 0]]
            proton_parent_dt = 0
            proton_parent_id = trajectories_parent_ids[proton_mask][ii]

            if proton_parent_id != -1:
                proton_parent_index = np.where(
                    (trajectories_traj_ids == proton_parent_id) & (trajectories_vertex_ids == vertex_id)
                )
                proton_parent = trajectories[proton_parent_index]
                parent_pdg = proton_parent['pdg_id']
                parent_E = proton_parent['E_start']
                parent_distance = np.sqrt(
                    (proton_parent['xyz_start'][0][0] - proton_parent['xyz_end'][0][0]) ** 2 + 
                    (proton_parent['xyz_start'][0][1] - proton_parent['xyz_end'][0][1]) ** 2 +
                    (proton_parent['xyz_start'][0][2] - proton_parent['xyz_end'][0][2]) ** 2
                )
                parent_pxyz_start = proton_parent['pxyz_start']
                parent_pxyz_end = proton_parent['pxyz_end']
                proton_parent_dt = (proton_t_start - vert_loc[3])
            
            """Get grandparent info"""
            grandparent_pdg = 0
            if proton_parent_id != -1:
                proton_grandparent_id = proton_parent['parent_id']
                if proton_grandparent_id != -1:
                    proton_grandparent_index = np.where(
                        (trajectories_traj_ids == proton_grandparent_id) & (trajectories_vertex_ids == vertex_id)
                    )
                    grandparent_pdg = trajectories[proton_grandparent_index]['pdg_id']

            """############################### Generate track data ###############################"""
            proton_charge_xyz = np.array([
                charge_x[proton_hits],
                charge_y[proton_hits],
                charge_z[proton_hits]
            ]).T
            proton_charge_E = charge_E[proton_hits]
            proton_xyz_start = trajectories_xyz_start[proton_mask][ii]
            proton_xyz_end = trajectories_xyz_end[proton_mask][ii]
            proton_pxyz_start = trajectories_pxyz_start[proton_mask][ii]
            proton_pxyz_end = trajectories_pxyz_end[proton_mask][ii]
            proton_E = trajectories_E[proton_mask][ii]
            
            """Check whether proton is fiducialized"""
            proton_xyz_start_fiducialized = self.fiducialized_vertex(proton_xyz_start)
            proton_xyz_end_fiducialized = self.fiducialized_vertex(proton_xyz_end)
            
            """
            To do this we find the closest hits to the actual track beginning
            and ending from the trajectories.  If these points end up being the
            same, we skip this track, since there's no way to distinguish
            that point as a track anyways.
            """
            proton_start_distances = np.sqrt(
                ((proton_charge_xyz - proton_xyz_start) ** 2).sum(axis=1)
            )
            proton_end_distances = np.sqrt(
                ((proton_charge_xyz - proton_xyz_end) ** 2).sum(axis=1)
            )
            closest_start_index = np.argmin(proton_start_distances)
            closest_end_index = np.argmin(proton_end_distances)

            proton_start_index = (
                proton_hits[closest_start_index],
                proton_segments[closest_start_index]
            )
            proton_end_index = (
                proton_hits[closest_end_index],
                proton_segments[closest_end_index]
            )
            
            """Generate track fit data"""
            track_data = self.fit_track(proton_hit_t0s, proton_charge_xyz)
            
            """Determine amount by which hits are shared from things other than the proton"""
            proton_hit_segment_ids = charge_back_track_segments[np.unique(proton_hits)]
            true_proton_hit_purity = float(len(proton_hits)) / float(np.count_nonzero(proton_hit_segment_ids))
            
            """############################### Generate clusters ###############################"""
            linked = fastcluster.linkage_vector(particle_charge_xyz, method='ward')
            
            n_clusters = len(linked)
            cluster_mapping = {i: [i] for i in range(len(particle_charge_xyz))}

            # Initialize a dictionary to hold lists of indices for each cluster
            clusters = {i: [] for i in range(-1, n_clusters)}
            
            for ii, row in enumerate(linked):
                cluster_1, cluster_2 = int(row[0]), int(row[1])
                new_cluster_index = len(particle_charge_xyz) + ii  # New index for the newly formed cluster

                # Combine the points from the clusters being merged and update the mapping
                combined_points = cluster_mapping[cluster_1] + cluster_mapping[cluster_2]
                cluster_mapping[new_cluster_index] = combined_points
                clusters[new_cluster_index] = combined_points
            
            # Iterate over the linkage matrix
            lowest_completeness_cluster = 0
            lowest_completeness_label = 0
            lowest_completeness = 0.0
            lowest_completeness_purity = 0.0
            
            """Set up list of proton points"""
            proton_points = np.zeros(len(np.unique(particle_hits)), dtype=bool)
            particle_proton_points = np.where(
                np.isin(np.unique(proton_hits), np.unique(particle_hits))
            )[0]
            proton_points[particle_proton_points] = True
            
            for ii, (label, combined_points) in enumerate(clusters.items()):
                if label == -1:
                    continue
                particle_points = np.zeros(len(np.unique(particle_hits)), dtype=bool)
                particle_points[combined_points] = True
                    
                if not sum(particle_points):
                    continue
                
                completeness = sum(particle_points & proton_points) / sum(proton_points)
                purity = sum(particle_points & proton_points) / sum(particle_points)
                
                if (completeness >= lowest_completeness) and (lowest_completeness < 1.0):
                    lowest_completeness = completeness
                    lowest_completeness_label = label
                    lowest_completeness_cluster = ii
                    lowest_completeness_purity = purity

            """Construct output CAF"""
            nar_inelastic_data_type = np.dtype([
                ('event_id', 'i4'),
                ('vertex_id', 'i4'),
                ('proton_id', 'i4'),
                ('proton_xyz_start', 'f4', (1, 3)),
                ('nu_vertex', 'f4', (1, 3)),
                ('proton_total_energy', 'f4'),
                ('proton_vis_energy', 'f4'),
                ('proton_length', 'f4'),
                ('proton_pxyz_start', 'f4', (1, 3)),
                ('proton_pxyz_end', 'f4', (1, 3)),
                ('parent_total_energy', 'f4'),
                ('parent_length', 'f4'),
                ('parent_pxyz_start', 'f4', (1, 3)),
                ('parent_pxyz_end', 'f4', (1, 3)),
                ('nu_proton_dt', 'f4'),
                ('nu_proton_distance', 'f4'),
                ('parent_pdg', 'i4'),
                ('grandparent_pdg', 'i4'),
                ('primary_pdg', 'i4'),
                ('primary_length', 'f4'),
                ('neutrino_vertex_fiducialized', 'i4'),
                ('proton_xyz_start_fiducialized', 'i4'),
                ('proton_xyz_end_fiducialized', 'i4'),
                ('truth_segment_overlap', 'f4'),
                ('best_completeness_cluster', 'i4'),
                ('best_completeness', 'f4'),
                ('best_purity', 'f4')
            ])
            nar_inelastic_data = np.array([(
                event,
                vertex_id,
                proton_id,
                proton_xyz_start,
                vert_loc[:3],
                proton_E,
                sum(proton_charge_E),
                np.sqrt(
                    (proton_xyz_start[0] - proton_xyz_end[0]) ** 2 + 
                    (proton_xyz_start[1] - proton_xyz_end[1]) ** 2 +
                    (proton_xyz_start[2] - proton_xyz_end[2]) ** 2
                ),
                proton_pxyz_start,
                proton_pxyz_end,
                parent_E,
                parent_distance,
                parent_pxyz_start,
                parent_pxyz_end,
                proton_parent_dt,
                np.sqrt(
                    (proton_xyz_start[0] - vert_loc[0]) ** 2 + 
                    (proton_xyz_start[1] - vert_loc[1]) ** 2 +
                    (proton_xyz_start[2] - vert_loc[2]) ** 2
                ),
                parent_pdg,
                grandparent_pdg,
                0,
                0,
                neutrino_vertex_fiducialized,
                proton_xyz_start_fiducialized,
                proton_xyz_end_fiducialized,
                true_proton_hit_purity,
                lowest_completeness_cluster,
                lowest_completeness,
                lowest_completeness_purity)],
                dtype=nar_inelastic_data_type
            )
            event_products['nar_inelastic'].append(nar_inelastic_data)
                