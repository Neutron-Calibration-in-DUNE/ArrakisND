"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
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

        self.input_products = ['daughters', 'track_id_hit_map', 'track_id_hit_segment_map']
        self.output_products = None

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

        trajectories_traj_ids = trajectories['traj_id']
        trajectories_vertex_ids = trajectories['vertex_id']
        trajectories_pdg_ids = trajectories['pdg_id']
        trajectories_xyz_start = trajectories['xyz_start']
        trajectories_xyz_end = trajectories['xyz_end']
        charge_x = charge['x']
        charge_y = charge['y']
        charge_z = charge['z']

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

            muon_xyz_start = trajectories_xyz_start[muon_mask][ii]
            muon_xyz_end = trajectories_xyz_end[muon_mask][ii]

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

            """
            Now we must determine if the track is broken along its
            trajectory.  We can do this by checking if the beginning and
            ending points are in different TPCs.  If so, then we will have
            to determine the beginning and ending points of each tracklette.
            """


            """Set track beginning and ending points"""
            arrakis_charge['tracklette_begin'][muon_start_index] = 1
            arrakis_charge['tracklette_end'][muon_end_index] = 1

        """Write changes to arrakis_file"""
        arrakis_file['charge_segment/calib_final_hits/data'][event_indices['charge']] = arrakis_charge
