"""
"""

import os
import sys
import glob
import h5flow
import h5py
from tqdm import tqdm
import numpy as np
from matplotlib import pyplot as plt
from h5flow.core import H5FlowStage, H5FlowDataManager, resources

from arrakis_nd.arrakis.common import process_types
from arrakis_nd.utils.logger import Logger
from arrakis_nd.utils.config import ConfigParser
from arrakis_nd.wrangler.simulation_wrangler import SimulationWrangler
from arrakis_nd.labeling_logic.simulation_labeling_logic import SimulationLabelingLogic


class Arrakis(H5FlowStage):
    class_version = '0.0.1'

    default_custom_param = None
    default_obj_name = 'obj0'
    """

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
    def __init__(
        self,
        name:   str = 'default',
        config: dict = {},
        meta:   dict = {},
        classname:      str = 'none',
        data_manager:   H5FlowDataManager = None,
        requires:       list = None,
        **params
    ):
        self.name = name + "_arrakis_nd"
        self.config = config
        self.meta = meta

        if "device" in self.meta:
            self.device = self.meta['device']
        else:
            self.device = 'cpu'
        if meta['verbose']:
            self.logger = Logger(self.name, output="both",   file_mode="w")
        else:
            self.logger = Logger(self.name, level='warning', file_mode="w")

        super(Arrakis, self).__init__(
            name, classname, data_manager, requires, **params
        )

        self.parse_config()

    def init(
        self,
        source_name
    ):
        """
        This method is run if arrakis is being used in an H5Flow stage.
        The labels for the dataset are then added as an extra data product
        for the various arrays (hits).
        """
        pass

    def set_config(
        self,
        config_file
    ):
        self.logger.info(f"parsing config file: {config_file}.")
        self.config_file = config_file
        self.config = ConfigParser(self.config_file).data
        self.parse_config()

    def parse_config(self):
        """
        """
        self.check_config()
        self.parse_dataset_folder()
        self.parse_dataset_files()
        self.parse_simulation_wrangler()
        self.parse_simulation_labeling_logic()

    def check_config(self):
        if "process_type" not in self.config:
            self.logger.warn('process type not specified in config! setting to "npz"')
            self.config['process_type'] = 'npz'
        else:
            if self.config['process_type'] not in process_types:
                self.logger.error(f'specified process_type {self.config["process_type"]} not allowed!')

    def parse_dataset_folder(self):
        # default to what's in the configuration file. May decide to deprecate in the future
        if ("simulation_folder" in self.config):
            self.simulation_folder = self.config["simulation_folder"]
            self.logger.info(
                "Set simulation file folder from configuration. " +
                f" simulation_folder : {self.simulation_folder}"
            )
        elif ('ARRAKIS_ND_SIMULATION_PATH' in os.environ):
            self.logger.debug('Found ARRAKIS_ND_SIMULATION_PATH in environment')
            self.simulation_folder = os.environ['ARRAKIS_ND_SIMULATION_PATH']
            self.logger.info(
                "Setting simulation path from Enviroment." +
                f" ARRAKIS_ND_SIMULATION_PATH = {self.simulation_folder}"
            )
        else:
            self.logger.error('No simulation_folder specified in environment or configuration file!')
        if not os.path.isdir(self.simulation_folder):
            self.logger.error(f'Specified simulation folder "{self.simulation_folder}" does not exist!')

    def parse_dataset_files(self):
        if ('simulation_files' not in self.config):
            self.logger.error('no simulation files specified in configuration file!')
        if isinstance(self.config["simulation_files"], list):
            self.simulation_files = [
                self.simulation_folder + input_file
                for input_file in self.config["simulation_files"]
            ]
        elif isinstance(self.config["simulation_files"], str):
            if self.config["simulation_files"] == "all":
                self.logger.info(f'searching {self.simulation_folder} recursively for all .npz files.')
                self.simulation_files = glob.glob(self.simulation_folder + '**/*.npz', recursive=True)
            else:
                try:
                    self.logger.info(
                        f'searching {self.simulation_folder} recursively for all {self.config["simulation_files"]} files.')
                    self.simulation_files = glob.glob(
                        self.simulation_folder + f'**/{self.config["simulation_files"]}',
                        recursive=True
                    )
                except:
                    self.logger.error(
                        f'specified "simulation_files" parameter: {self.config["simulation_files"]} incompatible!'
                    )
        else:
            self.logger.error(f'specified "simulation_files" parameter: {self.config["simulation_files"]} incompatible!')
        self.simulation_files = self.config['simulation_files']
        for ii, simulation_file in enumerate(self.simulation_files):
            if not os.path.isfile(self.simulation_folder + '/' + simulation_file):
                self.logger.error(f'specified file "{self.simulation_folder}/{simulation_file}" does not exist!')

    def parse_simulation_wrangler(self):
        self.simulation_wrangler = SimulationWrangler(
            self.name,
            self.config,
            self.meta
        )
        self.meta['simulation_wrangler'] = self.simulation_wrangler

    def parse_simulation_labeling_logic(self):
        self.simulation_labeling_logic = SimulationLabelingLogic(
            self.name,
            self.config,
            self.meta
        )
        self.meta['simulation_labeling_logic'] = self.simulation_labeling_logic

    def run(
        self,
        source_name,
        source_slice
    ):
        # load, process, and save new data objects
        pass

    def run_arrakis_nd(self):
        if self.config["process_type"] == "npz":
            self.run_arrakis_nd_npz()
        elif self.config["process_type"] == "flow":
            self.run_arrakis_nd_flow()
        else:
            self.logger.error(f'specified process_type {self.config["process_type"]} not allowed!')
    
    # TODO: break this up into smaller functions
    def run_arrakis_nd_npz(self):
        self.logger.info('running arrakis_nd in npz mode')
        for ii, simulation_file in enumerate(self.simulation_files):
            try:
                flow_file = h5flow.data.H5FlowDataManager(self.simulation_folder + '/' + simulation_file, "r")
            except:
                self.logger.error(f'there was a problem processing flow file {simulation_file}')
            
            trajectories = flow_file['mc_truth/trajectories/data']["event_id", "traj_id", "parent_id", "pdg_id", "start_process", "start_subprocess", "end_process", "end_subprocess", "E_end"]
            segments = flow_file['mc_truth/segments/data']['event_id', 'segment_id', 'traj_id']
            stacks = flow_file["mc_truth/stack/data"]['event_id']
            hits_back_track = flow_file["mc_truth/calib_final_hit_backtrack/data"]
            hits = flow_file["charge/calib_final_hits/data"]

            event_ids = np.unique(flow_file['mc_truth/segments/data']['event_id'])
            event_loop = tqdm(
                enumerate(event_ids, 0),
                total=len(event_ids),
                leave=True,
                position=1,
                colour='green'
            )
            for jj, event_id in event_loop:
                if event_id == 0:           # skip the first event for now, it is very large and takes forever to process
                    continue
                if event_id > 4:
                    self.simulation_labeling_logic.timers.evaluate_run()
                    self.simulation_labeling_logic.memory_trackers.evaluate_run()
                    break
                event_trajectories = trajectories[trajectories['event_id'] == event_id]
                event_segments = segments[segments['event_id'] == event_id]
                event_stacks = stacks[stacks == event_id]
                hits_back_track_mask = np.any(
                    np.isin(hits_back_track['segment_id'], event_segments['segment_id']),
                    axis=1
                )
                event_back_track_hits = hits_back_track[hits_back_track_mask]
                event_hits = hits[hits_back_track_mask]

                self.simulation_wrangler.process_event(
                    event_id,
                    event_trajectories,
                    event_segments,
                    event_stacks,
                    event_back_track_hits,
                    event_hits
                )
                self.simulation_labeling_logic.process_event()
                self.simulation_wrangler.save_event()
                event_loop.set_description(f"File: [{ii+1}/{len(self.simulation_files)}]")
                # event_loop.set_postfix_str(f"num_process={:.2e}")
            self.simulation_wrangler.save_events(
                "test.npz"
            )
        if self.simulation_labeling_logic.debug:
            self.simulation_labeling_logic.timers.evaluate_run()
            self.simulation_labeling_logic.memory_trackers.evaluate_run()

    def run_arrakis_nd_flow(self):
        pass
