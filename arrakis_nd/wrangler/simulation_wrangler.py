"""

"""
import numpy as np
from matplotlib import pyplot as plt
import copy
import h5py

from arrakis_nd.dataset.det_point_cloud import DetectorPointCloud

class SimulationWrangler:
    """
    """
    def __init__(self):

        self.det_point_cloud = DetectorPointCloud()
        self.det_point_clouds = {}

        self.trackid_parentid = {}
        self.trackid_pdgcode = {}
        self.trackid_process = {}
        self.trackid_subprocess = {}
        self.trackid_endprocess = {}
        self.trackid_endsubprocess = {}
        self.trackid_energy = {}
        self.trackid_daughters = {}
        self.trackid_progeny = {}
        self.trackid_descendants = {}
        self.trackid_ancestorlevel = {}
        self.trackid_ancestry = {} 

        self.trackid_segmentid = {}
        self.segmentid_trackid = {}

        self.trackid_hit = {}
        self.segmentid_hit = {}


    def clear_event(self):

        self.det_point_cloud.clear()

        self.trackid_parentid = {}
        self.trackid_pdgcode = {}
        self.trackid_process = {}
        self.trackid_subprocess = {}
        self.trackid_endprocess = {}
        self.trackid_endsubprocess = {}
        self.trackid_energy = {}
        self.trackid_daughters = {}
        self.trackid_progeny = {}
        self.trackid_descendants = {}
        self.trackid_ancestorlevel = {}
        self.trackid_ancestry = {}

        self.trackid_segmentid = {}
        self.segmentid_trackid = {}

        self.trackid_hit = {}
        self.segmentid_hit = {}
    
    def clear_point_clouds(self):
        self.det_point_clouds = {}
    
    def set_hit_labels(self,
        hit, trackid, topology, 
        particle, physics, unique_topology
    ):
        track_index = self.get_index_trackid(hit, trackid)
        if track_index != -1:
            self.det_point_cloud.topology_labels[hit][track_index] = topology.value
            self.det_point_cloud.particle_labels[hit][track_index] = self.trackid_pdgcode[particle]
            self.det_point_cloud.physics_labels[hit][track_index] = physics.value
            self.det_point_cloud.unique_topologies[hit][track_index] = unique_topology
            self.det_point_cloud.unique_particles[hit][track_index] = trackid
        self.det_point_cloud.topology_label[hit] = topology.value
        self.det_point_cloud.particle_label[hit] = self.trackid_pdgcode[particle]
        self.det_point_cloud.physics_label[hit] = physics.value
        self.det_point_cloud.unique_topology[hit] = unique_topology
        self.det_point_cloud.unique_particle[hit] = trackid

    def process_event(self,
        event_id,
        event_trajectories,
        event_segments,
        event_stacks,
        hits_back_track,
        hits
    ):
        self.clear_event()
        self.det_point_cloud.event = event_id
        self.process_event_trajectories(event_trajectories)
        self.process_event_stacks(event_stacks)
        self.process_event_segments(event_segments)
        self.process_event_hits(hits, hits_back_track)
    
    def save_event(self):
        self.det_point_cloud.x = np.array(self.det_point_cloud.x)
        self.det_point_cloud.y = np.array(self.det_point_cloud.y)
        self.det_point_cloud.z = np.array(self.det_point_cloud.z)
        self.det_point_cloud.t_drift = np.array(self.det_point_cloud.t_drift)
        self.det_point_cloud.ts_pps = np.array(self.det_point_cloud.ts_pps)
        self.det_point_cloud.Q = np.array(self.det_point_cloud.Q)
        self.det_point_cloud.E = np.array(self.det_point_cloud.E)
        self.det_point_cloud.source_label = np.array(self.det_point_cloud.source_label)
        self.det_point_cloud.topology_label = np.array(self.det_point_cloud.topology_label)
        self.det_point_cloud.particle_label = np.array(self.det_point_cloud.particle_label)
        self.det_point_cloud.physics_label = np.array(self.det_point_cloud.physics_label)
        self.det_point_cloud.unique_topology = np.array(self.det_point_cloud.unique_topology)
        self.det_point_cloud.unique_particle = np.array(self.det_point_cloud.unique_particle)
        self.det_point_cloud.unique_physics = np.array(self.det_point_cloud.unique_physics)
        self.det_point_clouds[self.det_point_cloud.event] = copy.deepcopy(self.det_point_cloud)

    def process_event_trajectories(self,
        event_trajectories
    ):
        for ii, particle in enumerate(event_trajectories):
            if abs(particle['pdg_id']) == 13:
                print(ii, particle['pdg_id'], particle['traj_id'], particle['start_process'], particle['parent_id'])
            track_id = particle['traj_id']                          
            self.trackid_parentid[track_id] = particle['parent_id']
            self.trackid_pdgcode[track_id] = particle['pdg_id']
            self.trackid_process[track_id] = particle['start_process']
            self.trackid_subprocess[track_id] = particle['start_subprocess']
            self.trackid_endprocess[track_id] = particle['end_process']
            self.trackid_endsubprocess[track_id] = particle['end_subprocess']
            self.trackid_energy[track_id] = particle['E_end'] # E_start or E_end?
            
            # iterate over daughters
            self.trackid_daughters[track_id] = []
            self.trackid_descendants[track_id] = []
            if particle['parent_id'] != -1:
                self.trackid_daughters[particle['parent_id']].append(track_id)
                self.trackid_descendants[particle['parent_id']].append(track_id)
            self.trackid_progeny[track_id] = []
            # iterate over ancestry
            level = 0
            mother = particle['parent_id']
            temp_track_id = particle['traj_id']
            ancestry = []
            while mother != -1:
                level += 1
                temp_track_id = mother
                ancestry.append(mother)
                mother = self.trackid_parentid[temp_track_id]
                
                if level > 1 and mother != -1:
                    self.trackid_progeny[mother].append(temp_track_id)
                    self.trackid_descendants[mother].append(temp_track_id)

            self.trackid_ancestorlevel[track_id] = level
            self.trackid_ancestry[track_id] = ancestry
            self.trackid_hit[track_id] = []
            self.trackid_segmentid[track_id] = []

    def process_event_stacks(self,
        event_stacks
    ):
        pass
        
    def process_event_segments(self,
        event_segments
    ):
        for ii, segment in enumerate(event_segments):
            self.trackid_segmentid[segment['traj_id']].append(segment['segment_id'])
            self.segmentid_trackid[segment['segment_id']] = segment['traj_id']
            self.segmentid_hit[segment['segment_id']] = []
    
    def process_event_hits(self,
        event_hits,
        event_hits_back_track
    ):
        for ii, hit in enumerate(event_hits):
            segment_ids = event_hits_back_track['segment_id'][ii]
            segment_fractions = event_hits_back_track['fraction'][ii]
            self.det_point_cloud.add_point(
                hit['x'], 
                hit['y'],
                hit['z'],
                hit['t_drift'], 
                hit['ts_pps'], 
                hit['Q'], 
                hit['E'], 
                segment_ids[(segment_ids != 0)],
                #segment_fractions[(segment_ids != 0)]
            )
            for segmentid in segment_ids[(segment_ids != 0)]:
                if segmentid in self.segmentid_hit.keys():
                    self.segmentid_hit[segmentid].append(ii)
                    self.trackid_hit[self.segmentid_trackid[segmentid]].append(ii)
    
    def get_total_hit_energy(self, 
        hits
    ):
        energy = 0.0
        for hit in hits:
            energy += self.det_point_cloud.E[hit]
        return energy

    def get_primaries_generator_label(self,
        label
    ):
        pass
    
    def get_primaries_pdg_code(self,
        pdg
    ):
        primaries = []
        for track_id, parent_id in self.trackid_parentid.items():
            if parent_id == 0 and self.trackid_pdgcode[track_id] == pdg:
                primaries.append(track_id)
        return primaries

    def get_primaries_abs_pdg_code(self,
        pdg
    ):
        primaries = []
        for track_id, parent_id in self.trackid_parentid.items():
            if parent_id == 0 and abs(self.trackid_pdgcode[track_id]) == abs(pdg):
                primaries.append(track_id)
        return primaries
    
    def get_hits_trackid(self,
        trackids
    ):
        hits = [
            self.trackid_hit[track_id]
            for track_id in trackids
        ]
        return hits
    
    def get_segments_trackid(self,
        trackids
    ):
        segments = [
            self.trackid_segmentid[track_id]
            for track_id in trackids
        ]
        return segments

    def get_trackid_pdg_code(self,
        pdg
    ):
        trackid = []
        for track_id, pdg_code in self.trackid_pdgcode.items():
            if pdg_code == pdg:
                trackid.append(track_id)
        return trackid   
    
    def get_daughters_trackid(self,
        trackids
    ):
        daughters = [
            self.trackid_daughters[track_id]
            for track_id in trackids
        ]
        return daughters
    
    def filter_trackid_not_pdg_code(self,
        trackids, pdg
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_pdgcode[track_id] != pdg
        ]
        return trackid

    def filter_trackid_pdg_code(self,
        trackids, pdg
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_pdgcode[track_id] == pdg
        ]
        return trackid

    def filter_trackid_not_abs_pdg_code(self,
        trackids, pdg
    ):
        trackid = [
            track_id for track_id in trackids
            if abs(self.trackid_pdgcode[track_id]) != abs(pdg)
        ]
        return trackid

    def filter_trackid_abs_pdg_code(self,
        trackids, pdg
    ):
        trackid = [
            track_id for track_id in trackids
            if abs(self.trackid_pdgcode[track_id]) == abs(pdg)
        ]
        return trackid

    def filter_trackid_not_process(self,
        trackids, process
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_process[track_id] != process
        ]
        return trackid

    def filter_trackid_process(self,
        trackids, process
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_process[track_id] == process
        ]
        return trackid
    
    def filter_trackid_not_subprocess(self,
        trackids, subprocess
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_subprocess[track_id] != subprocess
        ]
        return trackid

    def filter_trackid_subprocess(self,
        trackids, subprocess
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_subprocess[track_id] == subprocess
        ]
        return trackid

    def filter_trackid_not_endprocess(self,
        trackids, endprocess
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_endprocess[track_id] != endprocess
        ]
        return trackid

    def filter_trackid_endprocess(self,
        trackids, endprocess
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_endprocess[track_id] == endprocess
        ]
        return trackid
    
    def filter_trackid_not_endsubprocess(self,
        trackids, endsubprocess
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_endsubprocess[track_id] != endsubprocess
        ]
        return trackid

    def filter_trackid_endsubprocess(self,
        trackids, endsubprocess
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_endsubprocess[track_id] == endsubprocess
        ]
        return trackid
    
    def get_index_trackid(self,
        hit, trackid
    ):
        for ii, particle in enumerate(self.det_point_cloud.particle_label[hit]):
            if particle == trackid:
                return ii
        return -1