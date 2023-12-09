"""

"""
import uproot,os,getpass,socket
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import copy
import h5py

from arrakis_nd.dataset.det_point_cloud import DetectorPointCloud
from arrakis_nd.dataset.common import *

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

    #@profile
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
    
    #@profile
    def clear_point_clouds(self):
        self.det_point_clouds = {}
    
    #@profile
    def set_hit_labels(self,
        hit, segment_id, topology, 
        particle, physics, unique_topology
    ):
        point_cloud_index = self.get_index_trackid(hit, segment_id)
        if point_cloud_index != -1:
            self.det_point_cloud.data["topology_labels"][hit][point_cloud_index] = topology.value
            self.det_point_cloud.data["particle_labels"][hit][point_cloud_index] = self.trackid_pdgcode[particle]
            self.det_point_cloud.data["physics_labels"][hit][point_cloud_index] = physics.value
            self.det_point_cloud.data["unique_topologies"][hit][point_cloud_index] = unique_topology
            self.det_point_cloud.data["unique_particles"][hit][point_cloud_index] = segment_id
        self.det_point_cloud.data["topology_label"][hit] = topology.value
        self.det_point_cloud.data["particle_label"][hit] = self.trackid_pdgcode[particle]
        try:
            self.det_point_cloud.data["physics_label"][hit] = physics.value
        except:
            self.det_point_cloud.data["physics_label"][hit] = physics.value[0]
        self.det_point_cloud.data["unique_topology"][hit] = unique_topology
        self.det_point_cloud.data["unique_particle"][hit] = segment_id

    #@profile
    def process_event(self,
        event_id,
        event_trajectories,
        event_segments,
        event_stacks,
        hits_back_track,
        hits
    ):
        self.clear_event()
        self.det_point_cloud.data['event'] = event_id
        
        print("Processing trajectories")
        self.process_event_trajectories(event_trajectories)

        print("Processing stacks")
        self.process_event_stacks(event_stacks)

        print("Processing segments")
        self.process_event_segments(event_segments)

        print("Processing hits")
        self.process_event_hits(hits, hits_back_track)
        
        print("Done processing event")
    
    #@profile
    def save_event(self):
        self.det_point_cloud.data['x'] = np.array(self.det_point_cloud.data['x'])
        self.det_point_cloud.data['y'] = np.array(self.det_point_cloud.data['y'])
        self.det_point_cloud.data['z'] = np.array(self.det_point_cloud.data['z'])
        self.det_point_cloud.data['t_drift'] = np.array(self.det_point_cloud.data['t_drift'])
        self.det_point_cloud.data['ts_pps'] = np.array(self.det_point_cloud.data['ts_pps'])
        self.det_point_cloud.data['Q'] = np.array(self.det_point_cloud.data['Q'])
        self.det_point_cloud.data['E'] = np.array(self.det_point_cloud.data['E'])
        self.det_point_cloud.data['source_label'] = np.array(self.det_point_cloud.data['source_label'])
        self.det_point_cloud.data['topology_label'] = np.array(self.det_point_cloud.data['topology_label'])
        self.det_point_cloud.data['particle_label'] = np.array(self.det_point_cloud.data['particle_label'])
        self.det_point_cloud.data['physics_label'] = np.array(self.det_point_cloud.data['physics_label'])
        self.det_point_cloud.data['unique_topology'] = np.array(self.det_point_cloud.data['unique_topology'])
        self.det_point_cloud.data['unique_particle'] = np.array(self.det_point_cloud.data['unique_particle'])
        self.det_point_cloud.data['unique_physics'] = np.array(self.det_point_cloud.data['unique_physics'])
        self.det_point_clouds[self.det_point_cloud.data['event']] = copy.deepcopy(self.det_point_cloud)
 
        
    #@profile
    def save_events(self,
        simulation_file
    ):
        output_file = simulation_file.replace('.h5', '')
        output_file += '.arrakis_nd.npz'
        meta = {
            "who_created":      getpass.getuser(),
            "when_created":     datetime.now().strftime("%m-%d-%Y-%H:%M:%S"),
            "where_created":    socket.gethostname(),
            "det_features": {
                "x": 0, "y": 1, "z": 2, "Q": 3
            },
            "mc_features": {
                "t_drift": 0, "ts_pps": 1, "E": 2
            },
            "classes": {
                "source": 0, "topology": 1, "particle": 2, "physics": 3
            },
            "clusters": {
                "topology":  0, "particle": 1, "physics": 2
            },
            "source_labels": {
                key: value
                for key, value in classification_labels["source"].items()
            },
            "topology_labels": {
                key: value
                for key, value in classification_labels["topology"].items()
            },
            "particle_labels": {
                key: value
                for key, value in classification_labels["particle"].items()
            },      
            "physics_labels": {
                key: value
                for key, value in classification_labels["physics"].items()
            },  
            "hit_labels": {
                key: value
                for key, value in classification_labels["hit"].items()
            },    
        }
        np.savez(
            output_file,
            data=self.det_point_clouds,
            meta=meta,
        )

    #@profile
    def process_event_trajectories(self,
        event_trajectories
    ):
        for ii, particle in enumerate(event_trajectories):
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

    #@profile
    def process_event_stacks(self,
        event_stacks
    ):
        pass
        
    #@profile
    def process_event_segments(self,
        event_segments
    ):
        for ii, segment in enumerate(event_segments):
            self.trackid_segmentid[segment['traj_id']].append(segment['segment_id']) # segment_ids belonging to a track_id
            self.segmentid_trackid[segment['segment_id']] = segment['traj_id'] # track_id belonging to a segment_id
            self.segmentid_hit[segment['segment_id']] = [] # hit_ids belonging to a segment_id
    
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
                segment_ids,
                segment_fractions
            )
            for segmentid in segment_ids[(segment_ids != 0)]:
                if segmentid in self.segmentid_hit.keys():
                    self.segmentid_hit[segmentid].append(ii)
                    self.trackid_hit[self.segmentid_trackid[segmentid]].append(ii)
    
    #@profile
    def get_total_hit_energy(self, 
        hits
    ):
        energy = 0.0
        for hit in hits:
            try:
                energy += self.det_point_cloud.data['E'][hit]
            except:
                energy += 0.0
                print("Warning in get_total_hit_energy")
        return energy

    #@profile
    def get_primaries_generator_label(self,
        label
    ):
        pass
    
    #@profile
    def get_primaries_pdg_code(self,
        pdg
    ):
        primaries = []
        for track_id, parent_id in self.trackid_parentid.items():
            if parent_id == -1 and self.trackid_pdgcode[track_id] == pdg:
                primaries.append(track_id)
        return primaries

    #@profile
    def get_primaries_abs_pdg_code(self,
        pdg
    ):
        primaries = []
        for track_id, parent_id in self.trackid_parentid.items():
            if parent_id == -1 and abs(self.trackid_pdgcode[track_id]) == abs(pdg):
                primaries.append(track_id)
        return primaries
    
    #@profile
    def get_hits_trackid(self,
        trackids
    ):
        trackids_np = np.array(trackids).astype(int).flatten()
        hits = [
            self.trackid_hit[track_id]
            for track_id in trackids_np
        ]
        return hits
    
    #@profile
    def get_segments_trackid(self,
        trackids
    ):
        trackids_np = np.array(trackids).astype(int).flatten()
        segments = [
            self.trackid_segmentid[track_id]
            for track_id in trackids_np
        ]
        return segments

    #@profile
    def get_trackid_pdg_code(self,
        pdg
    ):
        trackid = []
        for track_id, pdg_code in self.trackid_pdgcode.items():
            if pdg_code == pdg:
                trackid.append(track_id)
        return trackid   
    
    #@profile
    def get_daughters_trackid(self,
        trackids
    ):
        daughters = [
            self.trackid_daughters[track_id]
            for track_id in trackids
        ]
        return daughters
    
    #@profile
    def filter_trackid_not_pdg_code(self,
        trackids, pdg
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_pdgcode[track_id] != pdg
        ]
        return trackid
    

    #@profile
    def filter_trackid_pdg_code(self,
        trackids, pdg
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_pdgcode[track_id] == pdg
        ]
        return trackid
    
    #@profile
    def filter_trackid_not_abs_pdg_code(self,
        trackids, pdg
    ):
        try:
            trackid = [
                track_id for track_id in trackids
                if abs(self.trackid_pdgcode[track_id]) != abs(pdg)
            ]
            return trackid
        except:
            print("Warning in filter_trackid_not_abs_pdg_code")
            print("Trackids", trackids)
            print("Pdg", pdg)
            return []

    #@profile
    def filter_trackid_abs_pdg_code(self,
        trackids, pdg
    ):
        try:
            trackid = [
                track_id for track_id in trackids
                if abs(self.trackid_pdgcode[track_id]) == abs(pdg)
            ]
            return trackid
        except:
            print("Warning in filter_trackid_abs_pdg_code")
            print("Trackids", trackids)
            print("Pdg", pdg)
            return []
    
    #@profile
    def filter_trackid_not_process(self,
        trackids, process
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_process[track_id] != process
        ]
        return trackid

    #@profile
    def filter_trackid_process(self,
        trackids, process
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_process[track_id] == process
        ]
        return trackid

    #@profile
    def filter_trackid_not_subprocess(self,
        trackids, subprocess
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_subprocess[track_id] != subprocess
        ]
        return trackid

    #@profile
    def filter_trackid_subprocess(self,
        trackids, subprocess
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_subprocess[track_id] == subprocess
        ]
        return trackid

    #@profile
    def filter_trackid_not_endprocess(self,
        trackids, endprocess
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_endprocess[track_id] != endprocess
        ]
        return trackid

    #@profile
    def filter_trackid_endprocess(self,
        trackids, endprocess
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_endprocess[track_id] == endprocess
        ]
        return trackid
    
    #@profile
    def filter_trackid_not_endsubprocess(self,
        trackids, endsubprocess
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_endsubprocess[track_id] != endsubprocess
        ]
        return trackid

    #@profile
    def filter_trackid_endsubprocess(self,
        trackids, endsubprocess
    ):
        trackid = [
            track_id for track_id in trackids
            if self.trackid_endsubprocess[track_id] == endsubprocess
        ]
        return trackid
    
    #@profile
    def get_index_trackid(self,
        hit, segment_id
    ):
        index = np.where(
            self.det_point_cloud.data['segment_ids'][hit] == segment_id
        )
        return index
