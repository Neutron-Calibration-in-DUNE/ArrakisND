"""
"""
import h5py
import numpy as np
from scipy.interpolate import splprep, splev
import warnings
from scipy.integrate import quad, IntegrationWarning

from arrakis_nd.utils.utils import profiler, integrand
from arrakis_nd.plugins.plugin import Plugin
from arrakis_nd.dataset.common import (
    ProcessType, SubProcessType,
    Topology, PhysicsMicro, PhysicsMacro
)


class MuonPlugin(Plugin):
    """
    A plugin for labeling muons in charge and light data.

    Muons are MIPs which deposit energy through MIP ionization (i.e. the
    standard dE/dx mechanism).  This plugin grabs all muons in an event
    and labels their energy deposits according to the following scheme:

        (1) hits coming directly from muon ionization are labeled as
                topology = track
                physics_micro = mip_ionization
                physics_macro = mip
        (2) delta electrons which are daughters of muons are labeled as
                topology = track
                physics_micro = electron_ionization
                physics_macro = delta_electron
        (3) michel electrons which are daughters of muons are labeled as
                topology = track
                physics_micro = electron_ionization
                physics_macro = michel_electron

    Apart from these labels, track beginning and ending points are also
    labelled.  Vertices are placed at points where Michel electrons are
    generated.

    StandardRecord objects are also generated.
    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        """
        super(MuonPlugin, self).__init__(config)

        self.input_products = [
            'daughters', 
            'track_id_hit_map', 
            'track_id_hit_segment_map',
            'track_id_hit_t0_map'
        ]
        self.output_products = ['tracklette', 'track']

        self.muon_labels = {
            'topology': Topology.Track.value,
            'physics_micro': PhysicsMicro.MIPIonization.value,
            'physics_macro': PhysicsMacro.MIP.value
        }

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
        charge = flow_file['charge/calib_final_hits/data'][event_indices['charge']]
        arrakis_charge = arrakis_file['charge_segment/calib_final_hits/data'][event_indices['charge']]
        track_id_hit_map = event_products['track_id_hit_map']
        track_id_hit_segment_map = event_products['track_id_hit_segment_map']
        track_id_hit_t0_map = event_products['track_id_hit_t0_map']

        trajectories_traj_ids = trajectories['traj_id']
        trajectories_vertex_ids = trajectories['vertex_id']
        trajectories_pdg_ids = trajectories['pdg_id']
        trajectories_xyz_start = trajectories['xyz_start']
        trajectories_xyz_end = trajectories['xyz_end']
        trajectories_E = trajectories['E_start']
        charge_x = charge['x']
        charge_y = charge['y']
        charge_z = charge['z']
        charge_Q = charge['Q']
        charge_E = charge['E']
        charge_io_group = charge['io_group']
        

        """Grab the muons"""
        muon_mask = (abs(trajectories_pdg_ids) == 13)

        """Iterate over the individual (traj_id, vertex_id, traj_index)"""
        for ii, (muon_id, vertex_id, traj_index) in enumerate(zip(
            trajectories_traj_ids[muon_mask],
            trajectories_vertex_ids[muon_mask],
            event_indices['trajectories'][muon_mask]
        )):
            """Get indices in charge_segments arrays"""
            muon_hits = track_id_hit_map[(muon_id, vertex_id)]

            """Check if muon has no hits, go to the next particle"""
            if not any(muon_hits):
                continue

            muon_segments = track_id_hit_segment_map[(muon_id, vertex_id)]
            muon_hit_segments = (
                muon_hits, muon_segments
            )
            
            """Get the associated t0 values"""
            muon_hit_t0s = track_id_hit_t0_map[(muon_id, vertex_id)]

            """Set the event_id"""
            arrakis_charge['event_id'][muon_hits] = event

            """Iterate over standard labels"""
            for label, value in self.muon_labels.items():
                arrakis_charge[label][muon_hit_segments] = value

            """Generate unique label"""
            arrakis_charge['unique_topology'][muon_hit_segments] = traj_index

            """Determine track begin and end points"""
            muon_charge_xyz = np.array([
                charge_x[muon_hits],
                charge_y[muon_hits],
                charge_z[muon_hits]
            ]).T
            muon_charge_Q = charge_Q[muon_hits]
            muon_charge_E = charge_E[muon_hits]
            muon_xyz_start = trajectories_xyz_start[muon_mask][ii]
            muon_xyz_end = trajectories_xyz_end[muon_mask][ii]
            muon_E = trajectories_E[muon_mask][ii]

            """
            To do this we find the closest hits to the actual track beginning
            and ending from the trajectories.  If these points end up being the
            same, we skip this track, since there's no way to distinguish
            that point as a track anyways.
            """
            muon_start_distances = np.sqrt(
                ((muon_charge_xyz - muon_xyz_start) ** 2).sum(axis=1)
            )
            muon_end_distances = np.sqrt(
                ((muon_charge_xyz - muon_xyz_end) ** 2).sum(axis=1)
            )
            closest_start_index = np.argmin(muon_start_distances)
            closest_end_index = np.argmin(muon_end_distances)

            muon_start_index = (
                muon_hits[closest_start_index],
                muon_segments[closest_start_index]
            )
            muon_end_index = (
                muon_hits[closest_end_index],
                muon_segments[closest_end_index]
            )
            
            if closest_start_index == closest_end_index:
                continue
            
            """Set track beginning and ending points"""
            arrakis_charge['tracklette_begin'][muon_start_index] = 1
            arrakis_charge['tracklette_end'][muon_end_index] = 1

            """
            Now we must determine if the track is broken along its
            trajectory.  We can do this by checking if the beginning and
            ending points are in different TPCs.  If so, then we will have
            to determine the beginning and ending points of each tracklette.
            
            We can do this by isolating hits according to the io_group,
            which defines the particular TPC region that a hit is in.  For
            each unique io_group associated to the hits of this particle,
            we will follow the same procedure and find the closest
            points to the track beginning and ending to get the tracklette 
            begin and end points.
            """
            muon_io_group = charge_io_group[muon_hits]
            unique_muon_io_group = np.unique(muon_io_group)
            muon_hits = np.array(muon_hits)
            muon_segments = np.array(muon_segments)
            muon_hit_t0s = np.array(muon_hit_t0s)
            muon_tracklette_ids = []
            for io_group in unique_muon_io_group:
                """Isolate hits from this io_group"""
                io_group_mask = (muon_io_group == io_group)
                if not sum(io_group_mask):
                    continue
                io_group_hits = muon_hits[io_group_mask]
                io_group_segments = muon_segments[io_group_mask]
                io_group_t0s = muon_hit_t0s[io_group_mask]
                
                """Find the closest points to track begin/end for this io_group"""
                io_group_start_distances = muon_start_distances[io_group_mask]
                io_group_end_distances = muon_end_distances[io_group_mask]
                closest_start_index = np.argmin(io_group_start_distances)
                closest_end_index = np.argmin(io_group_end_distances)

                io_group_start_index = (
                    io_group_hits[closest_start_index],
                    io_group_segments[closest_start_index]
                )
                io_group_end_index = (
                    io_group_hits[closest_end_index],
                    io_group_segments[closest_end_index]
                )
                
                if closest_start_index == closest_end_index:
                    continue
                
                arrakis_charge['tracklette_begin'][io_group_start_index] = 1
                arrakis_charge['tracklette_end'][io_group_end_index] = 1
                
                """Parameterize the trajectory of this tracklette using t0"""
                io_group_xyz = muon_charge_xyz[io_group_mask]
                combined = sorted(zip(io_group_t0s, io_group_xyz), key=lambda x: x[0])
                sorted_t0, sorted_xyz = zip(*combined)
                sorted_t0 = np.array(sorted_t0)
                sorted_xyz = np.array(sorted_xyz)

                try:
                    """Try to fit a spline curve to the xyz data"""
                    t_param = np.arange(len(sorted_xyz))
                    tck, u = splprep(
                        [sorted_xyz[:,0], sorted_xyz[:,1], sorted_xyz[:,2]], 
                        s=0
                    )
                    
                    """Try to integrate this curve"""
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", IntegrationWarning)
                        try:
                            curve_length, _ = quad(integrand, 0, 1, args=(tck,))
                        except Exception as e:
                            curve_length = 0
                        
                    """Get the derivatives at the beginning and ending points"""
                    dxdt_start, dydt_start, dzdt_start = splev(0, tck, der=1)
                    dxdt_end, dydt_end, dzdt_end = splev(1, tck, der=1)
                    dx_start = np.array([dxdt_start, dydt_start, dzdt_start])
                    dx_end = np.array([dxdt_end, dydt_end, dzdt_end])
                    dx_start_magnitude = np.linalg.norm(dx_start)
                    dx_end_magnitude = np.linalg.norm(dx_end)
                    
                    """Fill the values"""
                    tracklette_dir = [dxdt_start, dydt_start, dzdt_start] / dx_start_magnitude
                    tracklette_enddir = [dxdt_end, dydt_end, dzdt_end] / dx_end_magnitude
                    tracklette_len_gcm2 = 0
                    tracklette_len_cm = 0
                except Exception as e:
                    """Otherwise, these values are undefined"""
                    tracklette_dir = [0, 0, 0]
                    tracklette_enddir = [0, 0, 0]
                    tracklette_len_gcm2 = 0
                    tracklette_len_cm = 0
                
                """Now generate the associated tracklette CAF object"""
                """
                The tracklette data object has the following entries
                that must be filled.
                
                    tracklette_data_type = np.dtype([
                        ('event_id', 'i4'),
                        ('tracklette_id', 'i4'),
                        ('start', 'f4', (1, 3)),
                        ('end', 'f4', (1, 3)),
                        ('dir', 'f4', (1, 3)),
                        ('enddir', 'f4', (1, 3)),
                        ('Evis', 'f4'),
                        ('qual', 'f4'),
                        ('len_gcm2', 'f4'),
                        ('len_cm', 'f4'),
                        ('E', 'f4'),
                        ('truth', 'i4', (1, 20)),
                        ('truthOverlap', 'f4', (1, 20)),
                    ])
                """
                io_group_xyz = muon_charge_xyz[io_group_mask]
                
                tracklette_data_type = np.dtype([
                    ('event_id', 'i4'),
                    ('tracklette_id', 'i4'),
                    ('start', 'f4', (1, 3)),
                    ('end', 'f4', (1, 3)),
                    ('dir', 'f4', (1, 3)),
                    ('enddir', 'f4', (1, 3)),
                    ('Evis', 'f4'),
                    ('qual', 'f4'),
                    ('len_gcm2', 'f4'),
                    ('len_cm', 'f4'),
                    ('E', 'f4'),
                    ('truth', 'i4', (1, 20)),
                    ('truthOverlap', 'f4', (1, 20)),
                ])
                
                """Assign this tracklette to the muon track"""
                muon_tracklette_ids.append(io_group_end_index[0])
                
                """Create the CAF object for this tracklette"""
                tracklette_data = np.array([(
                    event, 
                    io_group_end_index[0], 
                    [io_group_xyz[closest_start_index]], 
                    [io_group_xyz[closest_end_index]], 
                    [tracklette_dir], 
                    [tracklette_enddir], 
                    sum(muon_charge_E[io_group_mask]), 
                    1, 
                    tracklette_len_gcm2, 
                    tracklette_len_cm, 
                    0, 
                    [[traj_index]+[0]*19], 
                    [[1]+[0]*19])], 
                    dtype=tracklette_data_type
                )

                """Add the new tracklette to the CAF objects"""
                event_products['tracklette'].append(tracklette_data)
                
        """Write changes to arrakis_file"""
        arrakis_file['charge_segment/calib_final_hits/data'][event_indices['charge']] = arrakis_charge
