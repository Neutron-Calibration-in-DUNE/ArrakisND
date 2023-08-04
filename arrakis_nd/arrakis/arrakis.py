"""
"""
import numpy as np
from matplotlib import pyplot as plt
import os
import sys
import h5py

from arrakis_nd.utils.logger import Logger, default_logger
from arrakis_nd.utils.config import ConfigParser
from arrakis_nd.wrangler.simulation_wrangler import SimulationWrangler

class Arrakis:
    """
    The module class helps to organize meta data and objects related to different tasks
    and execute those tasks based on a configuration file.  The spirit of the 'Module' class
    is to mimic some of the functionality of LArSoft, e.g. where you can specify a chain
    of tasks to be completed, the ability to have nested config files where default parameters
    can be overwritten.
    """
    def __init__(self,
        config_file:    str="",
    ):
        self.simulation_wrangler = SimulationWrangler()
        self.simulation_labeling_logic = None

        self.config_file = config_file
        if config_file != "":
            self.config = ConfigParser(config_file).data
            self.name = "arrakis_nd"
            self.logger = Logger(self.name, output="both", file_mode='w')
            system_info = self.logger.get_system_info()
            for key, value in system_info.items():
                self.logger.info(f"system_info - {key}: {value}")
            self.logger.info(f"configuring arrakis_nd.")
        else:
            default_logger.error(f"no config file specified for arrakis_nd at constructor.")
        
        self.logger.info(f"parsing config file: {config_file}.")
        self.meta = {}
        self.parse_config()

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
        if ("simulation_folder" in self.config['arrakis_nd']):
            self.simulation_folder = self.config['arrakis_nd']["simulation_folder"]
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
        if ('simulation_files' not in self.config['arrakis_nd']):
            self.logger.error(f'No simulation files specified in configuration file!')
        self.simulation_files = self.config['arrakis_nd']['simulation_files']
        for ii, simulation_file in enumerate(self.simulation_files):
            if not os.path.isfile(self.simulation_folder + '/' + simulation_file):
                self.logger.error(f'Specified file "{self.simulation_folder}/{simulation_file}" does not exist!')

    def run_arrakis_nd(self):
        for ii, simulation_file in enumerate(self.simulation_files):
            temp_file = h5py.File(self.simulation_folder + '/' + simulation_file, 'r')
            trajectory_events = temp_file['mc_truth']['trajectories']['data']['eventID']
            track_events = temp_file['mc_truth']['tracks']['data']['eventID']
            charge_events = temp_file['charge']['events']['data']
            print(charge_events)
            unique_trajectory_events = np.unique(trajectory_events)
            for event in unique_trajectory_events:
                trajectory_event_mask = (trajectory_events == event)
                track_event_mask = (track_events == event)
                charge_event_mask = (charge_events == event)
                self.simulation_wrangler.process_event(
                    event_trajectories=temp_file['mc_truth']['trajectories']['data'][trajectory_event_mask],
                    event_tracks=temp_file['mc_truth']['tracks']['data'][track_event_mask],
                    event_charge=temp_file['charge']['calib_final_hits']['data'][charge_event_mask]
                )


