"""
"""
import numpy as np
from matplotlib import pyplot as plt
import os
import sys
import h5flow
import h5py
import numpy.lib.recfunctions as rfn
from collections import defaultdict
import json
from tqdm import tqdm

from h5flow.core import H5FlowStage, H5FlowDataManager, resources

from arrakis_nd.utils.logger import Logger, default_logger
from arrakis_nd.utils.config import ConfigParser
from arrakis_nd.wrangler.simulation_wrangler import SimulationWrangler
from arrakis_nd.labeling_logic.simulation_labeling_logic import SimulationLabelingLogic

class Arrakis(H5FlowStage):
    class_version = '0.0.0' # keep track of a version number for each class

    default_custom_param = None
    default_obj_name = 'obj0'
    """
    The module class helps to organize meta data and objects related to different tasks
    and execute those tasks based on a configuration file.  The spirit of the 'Module' class
    is to mimic some of the functionality of LArSoft, e.g. where you can specify a chain
    of tasks to be completed, the ability to have nested config files where default parameters
    can be overwritten.

    We'd like to put this into H5FlowStage format so that it can be used in larnd-sim.

    larnd-sim files have the following sets of arrays

    trajectories (mc truth):
        These are the true particle trajectories (or paths) through the detector for all particles, 
        both neutral and charged, excluding the incident neutrino. Each true particle may have multiple 
        trajectories if the trajectory was split/broken by edep-sim with each having their own unique track ID.

        event_id:       unique ID for an interesting window of time; for beam events this corresponds to a spill
        vertex_id:      the vertex ID number, corresponds to an individual generator interaction
        traj_id:        the monotonic trajectory (track) ID, guaranteed to be unique within a file
        local_traj_id:  the original edep-sim trajectory (track) ID, may not be unique
        parent_id:      the trajectory (track) ID of the parent trajectory, if the trajectory is a primary particle the ID is -1
        E_start:        the total energy in [MeV] at the start of the trajectory
        pxyz_start:     the momentum 3-vector (px, py, pz) in [MeV] at the start of the trajectory
        xyz_start:      the start position 3-vector (x, y, z) in [cm] of the trajectory (specifically the position of the first trajectory point)
        t_start:        the start time of the trajectory in [us]
        E_end:          the total energy in [MeV] at the end of the trajectory
        pxyz_end:       the momentum 3-vector (px, py, pz) in [MeV] at the end of the trajectory
        xyz_end:        the end position 3-vector (x, y, z) in [cm] of the trajectory (specifically the position of the last trajectory point)
        t_end:          the end time of the trajectory in [us]
        pdg_id:         the PDG code of the particle
        start_process:  physics process for the start of the trajectory as defined by GEANT4
        start_subprocess: physics subprocess for the start of the trajectory as defined by GEANT4
        end_process:    physics process for the end of the trajectory as defined by GEANT4
        end_subprocess: physics subprocess for the end of the trajectory as defined by GEANT4

    segments (energy depositions):
        These are the true energy deposits (or energy segments) for active parts of the detector from edep-sim. 
        Each segment corresponds to some amount of energy deposited over some distance. Some variables are filled 
        during the larndsim stage of processing.

        event_id:       unique ID for an interesting window of time; for beam events this corresponds to a spill
        vertex_id:      the vertex ID number, corresponds to an individual generator interaction
        segment_id:     the segment ID number
        traj_id:        the trajectory (track) ID of the edep-sim trajectory that created this energy deposit
        x_start:        the x start position [cm]
        y_start:        the y start position [cm]
        z_start:        the z start position [cm]
        t0_start:       the start time [us]
        x_end:          the x end position [cm]
        y_end:          the y end position [cm]
        z_end:          the z end position [cm]
        t0_end:         the start time [us]
        x:              the x mid-point of the segment [cm] -> (x_start + x_end) / 2
        y:              the y mid-point of the segment [cm] -> (y_start + y_end) / 2
        z:              the z mid-point of the segment [cm] -> (z_start + z_end) / 2
        t0:             the time mid-point [us] -> (t0_start + t0_end) / 2
        pdg_id:         PDG code of the particle that created this energy deposit
        dE:             the energy deposited in this segment [MeV]
        dx:             the length of this segment [cm]
        dEdx:           the calculated energy per length [MeV/cm]
        tran_diff:      (ADD INFO)
        long_diff:      (ADD INFO)
        n_electrons:    (ADD INFO)
        n_photons:      (ADD INFO)
        pixel_plane:    (ADD INFO)
        t/t_start/t_end: (ADD INFO)
    
    flow files have the following sets of arrays

    charge:

    mc_truth:
        calib_final_hit_backtrack:
            fraction:   fraction of the segment associated to the hit
            segment_id: segment id associated to the hit
        interactions:
        light:
            segment_id: segment id associated to the hit
            n_photons_det:
            t0_det:
        packet_fraction:
        segments:
        stack:
        trajectories:

    Associations between calib_final_hits and particles/edeps can be made with the 'calib_final_hit_backtrack' 
    array inside of the mc_truth dataset in the flow files.  Each calib_final_hit has an associated segment id and a
    fraction of the edep that corresponds to the hit.

    First, we'll have to collect information using an event mask, and then arange the mc_truth info for particles
    and edeps.  Then, we will construct the 3d charge points and apply the labeling logic.
    
    H5FlowStage
    '''
        Base class for loop stage. Provides the following attributes:

         - ``name``: instance name of stage (declared in configuration file)
         - ``classname``: stage class
         - ``class_version``: a ``str`` version number (``'major.minor.fix'``, default = ``'0.0.0'``)
         - ``data_manager``: an ``H5FlowDataManager`` instance used to access the output file
         - ``requires``: a list of dataset names to load when calling ``H5FlowStage.load()``
         - ``comm``: MPI world communicator (if needed, else ``None``)
         - ``rank``: MPI group rank
         - ``size``: MPI group size

         To build a custom stage, inherit from this base class and implement
         the ``init()`` and the ``run()`` methods.

         Example::

            class ExampleStage(H5FlowStage):
                class_version = '0.0.0' # keep track of a version number for each class

                default_custom_param = None
                default_obj_name = 'obj0'

                def __init__(**params):
                    super(ExampleStage,self).__init__(**params)

                    # grab parameters from configuration file here, e.g.
                    self.custom_param = params.get('custom_param', self.default_custom_param)
                    self.obj_name = self.name + '/' + params.get('obj_name', self.default_obj_name)

                def init(self, source_name):
                    # declare any new datasets and set dataset metadata, e.g.

                    self.data_manager.set_attrs(self.obj_name,
                        classname=self.classname,
                        class_version=self.class_version,
                        custom_param=self.custom_param,
                        )
                    self.data_manager.create_dset(self.obj_name)

                def run(self, source_name, source_slice):
                    # load, process, and save new data objects

                    data = self.load(source_name, source_slice)
    """
    def __init__(self,
        name:   str='default',
        config: dict={},
        meta:   dict={},
        classname:      str='none',
        data_manager:   H5FlowDataManager=None,
        requires:       list=None,
        **params
    ):
        super(Arrakis, self).__init__(
            name, classname, data_manager, requires, **params
        )
        self.simulation_wrangler = SimulationWrangler()
        self.simulation_labeling_logic = SimulationLabelingLogic(self.simulation_wrangler)

        self.name = name + "_arrakis_nd"
        self.config = config
        self.meta = meta
        self.logger = Logger(self.name, output="both", file_mode='w')

        self.parse_config()
    
    def init(self, source_name):
        """
        This method is run if arrakis is being used in an H5Flow stage.
        The labels for the dataset are then added as an extra data product
        for the various arrays (hits).
        """
        pass

    def set_config(self,
        config_file:    str
    ):
        self.logger.info(f"parsing config file: {config_file}.")
        self.config_file = config_file
        self.config = ConfigParser(self.config_file).data
        self.parse_config()
    
    def parse_config(self):
        """
        """
        self.check_config()
        self.parse_input_files()
        self.run_arrakis_nd()

    def check_config(self):
        pass

    def parse_input_files(self):
        # default to what's in the configuration file. May decide to deprecate in the future
        if ("simulation_folder" in self.config):
            self.simulation_folder = self.config["simulation_folder"]
            self.logger.info(
                f"Set simulation file folder from configuration. " +
                f" simulation_folder : {self.simulation_folder}"
            )
        elif ('ARRAKIS_ND_SIMULATION_PATH' in os.environ ):
            self.logger.debug(f'Found ARRAKIS_ND_SIMULATION_PATH in environment')
            self.simulation_folder = os.environ['ARRAKIS_ND_SIMULATION_PATH']
            self.logger.info(
                f"Setting simulation path from Enviroment." +
                f" ARRAKIS_ND_SIMULATION_PATH = {self.simulation_folder}"
            )
        else:
            self.logger.error(f'No simluation_folder specified in environment or configuration file!')
        if not os.path.isdir(self.simulation_folder):
            self.logger.error(f'Specified simulation folder "{self.simulation_folder}" does not exist!')
        if ('simulation_files' not in self.config):
            self.logger.error(f'No simulation files specified in configuration file!')
        self.simulation_files = self.config['simulation_files']
        for ii, simulation_file in enumerate(self.simulation_files):
            if not os.path.isfile(self.simulation_folder + '/' + simulation_file):
                self.logger.error(f'Specified file "{self.simulation_folder}/{simulation_file}" does not exist!')

    def run(self, source_name, source_slice):
        # load, process, and save new data objects
        pass

    def run_arrakis_nd(self):
        for ii, simulation_file in enumerate(self.simulation_files):
            #flow_file = h5py.File(self.simulation_folder + '/' + simulation_file, 'r')
            flow_file = h5flow.data.H5FlowDataManager(self.simulation_folder + '/' + simulation_file, "r")
            try:
                charge = flow_file['charge']
                combined = flow_file['combined']
                geometry_info = flow_file['geometry_info']
                lar_info = flow_file['lar_info']
                light = flow_file['light']
                mc_truth = flow_file['mc_truth']
                run_info = flow_file['run_info']
            except:
                self.logger.error(f'there was a problem processing flow file {simulation_file}')
            
            truth = flow_file["/mc_truth/stack", "/mc_truth/trajectories", "/mc_truth/segments"]
            hits = flow_file["/mc_truth/calib_final_hit_backtrack", "/charge/calib_final_hits", "/charge/calib_prompt_hits", "/charge/packets", "/mc_truth/segments"]

            #trajectories = mc_truth['trajectories']['data']
            #segments = mc_truth['segments']['data']
            #stacks = mc_truth['stack']['data']
            #hits_back_track = mc_truth['calib_final_hit_backtrack']['data']
            #hits = charge['calib_final_hits']['data']

            trajectory_events = truth[['event_id','traj_id']]
            segment_events = truth[['event_id','segment_id']]
            stack_events = truth['event_id']
            hit_events = hits["event_id", "segment_id"]

            #unique_events = np.unique(segment_events['segment_id']) # gives the unique segment ids
            unique_events = np.unique(stack_events) # I think we want the event_id instead?
            
            event_loop = tqdm(
                enumerate(unique_events, 0), 
                total=len(unique_events), 
                leave=True,
                position=0,
                colour='green'
            )
            for ii, event in event_loop:
                trajectory_event_mask = (trajectory_events["event_id"] == event)
                segment_event_mask = (segment_events["event_id"] == event)
                stack_event_mask = (stack_events == event)
                
                
                hits_back_track_mask = np.any(
                    np.isin(hit_events['event_id'], truth[segment_event_mask]['event_id']), # I think this is still wrong
                    axis=1
                )
                # Work In Progress
                self.simulation_wrangler.process_event(
                    event,
                    event_trajectories=trajectory_events[trajectory_event_mask],
                    event_segments=segment_events[segment_event_mask],
                    event_stacks=stack_events[stack_event_mask],
                    hits_back_track=hit_events[hits_back_track_mask],
                    hits=hit_events[hits_back_track_mask] # this is now a duplicate?
                )
                self.simulation_labeling_logic.process_event()
                self.simulation_wrangler.save_event()
                event_loop.set_description(f"Running ArrakisND - Event: [{ii+1}/{len(unique_events)}]")


