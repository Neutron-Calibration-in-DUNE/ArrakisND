"""
"""
from mpi4py import MPI
import h5py
import os
import importlib.util
import sys
import inspect
import glob
from tqdm import tqdm
import torch
import numpy as np
from datetime import datetime
import traceback
from matplotlib import pyplot as plt
from importlib.metadata import version

from arrakis_nd.utils.logger import Logger
from arrakis_nd.utils.utils import (
    profiler,
    timing_manager,
    memory_manager
)
from arrakis_nd.arrakis.common import (
    tracklette_data_type,
    track_data_type,
    fragment_data_type,
    shower_data_type,
    blip_data_type,
    particle_data_type,
    interaction_data_type,
    neutrino_data_type,
    interaction_event_data_type,
    neutrino_event_data_type,
    nar_inelastic_data_type,
    undefined_data_type,
    output_error_data_type
)
from arrakis_nd.plugins.plugin import Plugin


class Arrakis:
    """
    Main Arrakis class for running jobs. Arrakis is designed
    to work with MPI and H5.  The user must run Arrakis using mpirun
    with at least two processes (one master and N-1 workers).
    """
    @profiler
    def __init__(
        self,
        config: dict = {},
        meta:   dict = {},
    ):
        """_summary_

        Args:
            config (dict): config file for running Arrakis.
            meta (dict): dictionary of meta information to
            be shared across all nodes.
        """

        """Get mpi communication parameters"""
        try:
            self.comm = MPI.COMM_WORLD
            self.rank = self.comm.Get_rank()
            self.size = self.comm.Get_size()
        except Exception as exception:
            raise RuntimeError(f"unable to obtain MPI parameters: {exception}")

        """Set input parameters and set up loggers"""
        try:
            if self.rank == 0:
                self.logger = Logger(f"master: {self.rank}")
            else:
                self.logger = Logger(f"worker: {self.rank}")
        except Exception as exception:
            raise RuntimeError(f"unable to set up logging system: {exception}")

        self.config = config
        self.meta = meta

        """Setting error status for this node"""
        self.error_status = None
        self.exc_type = None
        self.exc_value = None
        self.exc_traceback = None
        self.line_number = None
        self.file_name = None
        self.tb_str = None
        self.traceback_details = None
        self.event_errors = []
        self.plugin_errors = []
        self.event_plugin_errors = []
        self.event_plugin_exc_types = []
        self.event_plugin_exc_values = []
        self.event_plugin_exc_tracebacks = []
        self.event_plugin_line_numbers = []
        self.event_plugin_file_names = []
        self.event_plugin_tb_strs = []
        self.event_plugin_traceback_details = []

        """Parse config"""
        try:
            self.parse_config()
        except Exception as exception:
            self.report_error(exception=exception)
        self.barrier()

        """Construct plugins"""
        self.plugins = {}
        try:
            self.construct_plugins()
        except Exception as exception:
            self.report_error(exception=exception)
        self.barrier()

        self.flow_files = []
        """Distributed events and indices for various arrays"""
        self.distributed_events = {}
        self.distributed_interactions_indices = {}
        self.distributed_segments_indices = {}
        self.distributed_stack_indices = {}
        self.distributed_trajectories_indices = {}
        self.distributed_charge_indices = {}
        self.distributed_light_indices = {}

    @profiler
    def clear_indices(self):
        self.distributed_events.clear()
        self.distributed_interactions_indices.clear()
        self.distributed_segments_indices.clear()
        self.distributed_stack_indices.clear()
        self.distributed_trajectories_indices.clear()
        self.distributed_charge_indices.clear()
        self.distributed_light_indices.clear()

    @profiler
    def report_error(
        self,
        exception: Exception = None
    ):
        self.error_status = exception
        self.exc_type, self.exc_value, self.exc_traceback = sys.exc_info()
        # Extracting the line number from the traceback
        self.line_number = self.exc_traceback.tb_lineno
        self.file_name = self.exc_traceback.tb_frame.f_code.co_filename
        # Optionally, use traceback to format a string of the entire traceback
        self.tb_str = traceback.format_exception(self.exc_type, self.exc_value, self.exc_traceback)
        self.traceback_details = "".join(self.tb_str)

    @profiler
    def barrier(self):
        errors = self.comm.allgather(self.error_status)
        exc_types = self.comm.allgather(self.exc_type)
        line_numbers = self.comm.allgather(self.line_number)
        file_names = self.comm.allgather(self.file_name)
        if any(errors):
            if self.rank == 0:
                errors_count = sum(1 for error in errors if error is not None)
                errors_indices = [index for index, error in enumerate(errors) if error is not None]
                self.logger.critical(f"{errors_count} errors encountered in Arrakis program")
                for index in errors_indices:
                    self.logger.critical(f"error encountered in worker {index}: ")
                    self.logger.critical(f"exception:   {errors[index]}")
                    self.logger.critical(f"exc_type:    {exc_types[index]}")
                    self.logger.critical(f"line_number:     {line_numbers[index]}")
                    self.logger.critical(f"file_name:       {file_names[index]}")
                self.comm.Abort(1)
        else:
            self.comm.Barrier()

    @profiler
    def parse_config(self):
        """
        Set up config parameters from input config file
        """
        if self.rank == 0:
            try:
                """Startup main arrakis program"""
                self.logger.info(
                    f'############################ ARRAKIS  v. [{version("arrakis-nd")}] ############################'
                )

                """Check for main arrakis_nd parameters"""
                if "arrakis_nd" not in self.config.keys():
                    self.logger.error("arrakis_nd section not in config!")
                arrakis_nd_config = self.config["arrakis_nd"]

                """Try to grab system info and display to the logger"""
                system_info = self.logger.get_system_info()
                time = datetime.now()
                now = f"{time.hour}:{time.minute}:{time.second} [{time.day}/{time.month}/{time.year}]"
                self.logger.info(f'system_info - local time: {now}')
                for key, value in system_info.items():
                    self.logger.info(f"system_info - {key}: {value}")

                """Set the verbosity of the Arrakis program"""
                if "verbose" in arrakis_nd_config:
                    if not isinstance(arrakis_nd_config["verbose"], bool):
                        self.logger.error(
                            f'"arrakis_nd:verbose" must be of type bool, but got {type(arrakis_nd_config["verbose"])}!'
                        )
                    self.meta["verbose"] = arrakis_nd_config["verbose"]
                else:
                    self.meta["verbose"] = False

                """See if any CUDA devices are available on the system"""
                if torch.cuda.is_available():
                    self.logger.info("CUDA is available with devices:")
                    for ii in range(torch.cuda.device_count()):
                        device_properties = torch.cuda.get_device_properties(ii)
                        self.logger.info(f" -- device: {ii}")
                        self.logger.info(f" {' ':<{5}} name: {device_properties.name}")
                        self.logger.info(f" {' ':<{5}} compute: {device_properties.major}.{device_properties.minor}")
                        self.logger.info(f" {' ':<{5}} memory: {round(device_properties.total_memory / (1024**3))} GB")

                """Check for GPU devices and configure them"""
                if "gpu" not in arrakis_nd_config.keys():
                    self.logger.info('"arrakis_nd:gpu" not specified in config')
                    gpu = None
                else:
                    gpu = arrakis_nd_config["gpu"]
                if "gpu_device" not in arrakis_nd_config.keys():
                    self.logger.info('"arrakis_nd:gpu_device" not specified in config')
                    gpu_device = None
                else:
                    gpu_device = arrakis_nd_config["gpu_device"]

                """Set the GPU parameters"""
                if gpu:
                    if torch.cuda.is_available():
                        if gpu_device >= torch.cuda.device_count() or gpu_device < 0:
                            self.logger.warn(
                                f"desired gpu_device '{gpu_device}' not available, using device '0'"
                            )
                            gpu_device = 0
                        self.meta["device"] = torch.device(f"cuda:{gpu_device}")
                        self.logger.info(
                            f"CUDA is available, using device {gpu_device}"
                            + f": {torch.cuda.get_device_name(gpu_device)}"
                        )
                    else:
                        gpu is False
                        self.logger.warn("CUDA not available! Using the cpu")
                        self.meta["device"] = torch.device("cpu")
                else:
                    self.logger.info("using cpu as device")
                    self.meta["device"] = torch.device("cpu")
                self.meta["gpu"] = gpu

                """Send new device meta data to workers"""
                for ii in range(1, self.comm.size):
                    self.comm.send(self.meta, dest=ii, tag=0)
            except Exception as exception:
                self.report_error(exception=exception)
            self.barrier()
        else:
            self.barrier()
            try:
                self.meta = self.comm.recv(source=0, tag=0)
            except Exception as exception:
                self.report_error(exception=exception)

    @profiler
    def construct_plugins(self):
        """
        Construct all plugin objects.  First, we look for all associated
        classes which inherit from plugins.  These include both the standard
        plugins which are merged into the plugins/ folder, as well as
        any which are specified by the user at run time.
        """
        """First collect all available plugins"""
        self.collect_plugins()

        """Check that plugins section exists and has some actual plugins in it"""
        if 'plugins' not in self.config.keys():
            self.logger.error('no plugins section specified in config!')
        plugins_config = self.config['plugins']
        if plugins_config is None:
            self.logger.error('no plugins specified in config!')

        """Iterate over specified plugins and check that they exist"""
        for plugin, config in plugins_config.items():
            if plugin not in self.available_plugins.keys():
                self.logger.error(f'specified plugin "{plugin}" not available!')
            try:
                self.plugins[plugin] = self.available_plugins[plugin](config)
            except Exception as exception:
                self.logger.error(f'failed to construct plugin {plugin} with config {config}: {exception}')

        """Check that intermediate products are set up correctly"""
        running_output_products = []
        for ii, plugin_name in enumerate(self.plugins.keys()):
            if ii == 0:
                if self.plugins[plugin_name].input_products is not None:
                    self.logger.error(
                        f'first specfied plugin "{plugin_name}" should have input_products == None, \
                        but has input_products == "{self.plugins[plugin_name].input_products}"!'
                    )
            else:
                if (
                    (self.plugins[plugin_name].input_products is None) or
                    (self.plugins[plugin_name].input_products == [])
                ):
                    pass
                else:
                    if isinstance(self.plugins[plugin_name].input_products, str):
                        if self.plugins[plugin_name].input_products not in running_output_products:
                            self.logger.error(
                                f'specified plugin "{plugin_name}" at location "{ii}" expects input_products \
                                "{self.plugins[plugin_name].input_products}" but current products only contain \
                                {running_output_products}!'
                            )
                    elif isinstance(self.plugins[plugin_name].input_products, list):
                        for input_product in self.plugins[plugin_name].input_products:
                            if input_product not in running_output_products:
                                self.logger.error(
                                    f'specified plugin "{plugin_name}" at location "{ii}" expects input_products \
                                    "{self.plugins[plugin_name].input_products}" but current products only contain \
                                    {running_output_products}!'
                                )
            if isinstance(self.plugins[plugin_name].output_products, str):
                running_output_products.append(self.plugins[plugin_name].output_products)
            elif isinstance(self.plugins[plugin_name].output_products, list):
                running_output_products += self.plugins[plugin_name].output_products
            else:
                if self.plugins[plugin_name].output_products is not None:
                    self.logger.error(
                        f'specified plugin "{plugin_name}" at location "{ii}" should have output_products == None, \
                        or output_products == str, or output_products == list(str), but has type \
                        "{type(self.plugins[plugin_name].output_products)}"!'
                    )

    @profiler
    def collect_plugins(self):
        """
        This function collects all available plugins from the ArrakisND
        plugins folder, as well as any plugin locations specified in the
        config file.
        """
        self.available_plugins = {}
        self.plugin_files = [
            os.path.dirname(__file__) + '/../plugins/' + file
            for file in os.listdir(path=os.path.dirname(__file__) + '/../plugins/')
        ]
        if 'plugin_locations' in self.config.keys():
            if (
                (not isinstance(self.config['plugin_locations'], list)) and
                (self.config['plugin_locations'] is not None)
            ):
                self.logger.error(
                    f"plugin_locations specified in config is of type {type(self.config['plugin_locations'])}, \
                    but must be a list or none!"
                )
            self.plugin_files.extend(self.config['plugin_locations'])

        """Iterate over found plugin files and attempt to load them"""
        for plugin_file in self.plugin_files:
            if (
                ("__init__.py" in plugin_file) or
                ("__pycache__.py" in plugin_file) or
                ("__pycache__" in plugin_file) or
                (".py" not in plugin_file)
            ):
                continue
            try:
                self.load_plugin_function(plugin_file)
            except Exception as exception:
                self.logger.warn(f'problem loading plugin from file: {plugin_file}: {exception}')

    @profiler
    def load_plugin_function(
        self,
        plugin_file: str
    ):
        """
        This function attempts to load a plugin which inherits from
        the Plugin class.

        Args:
            plugin_file (str): _description_
        """
        spec = importlib.util.spec_from_file_location(
            f'{plugin_file.removesuffix(".py")}.name',
            plugin_file
        )
        custom_plugin_file = importlib.util.module_from_spec(spec)
        sys.modules[f'{plugin_file.removesuffix(".py")}.name'] = custom_plugin_file
        spec.loader.exec_module(custom_plugin_file)
        for name, obj in inspect.getmembers(sys.modules[f'{plugin_file.removesuffix(".py")}.name']):
            if inspect.isclass(obj):
                custom_class = getattr(custom_plugin_file, name)
                if issubclass(custom_class, Plugin):
                    self.available_plugins[name] = custom_class

    @profiler
    def set_up_input_files(self):
        """
        Iterate over flow_folder directory and determine
        if listed input files exist, or if "all" is selected for
        input files, construct the list of all .h5 files in the
        flow_folder.

        We create the data members
            self.flow_folder - location of the flow files specified in config
            self.flow_files - names of all the flow files to process
            self.arrakis_folder - location to put the arrakis files
        """
        if self.rank == 0:
            try:
                arrakis_dict = self.config['arrakis_nd']

                """Check for parameters"""
                if 'flow_folder' not in arrakis_dict.keys():
                    self.logger.error('flow_folder not specified in config!')
                if 'flow_files' not in arrakis_dict.keys():
                    self.logger.error('flow_files not specified in config!')

                flow_folder = arrakis_dict['flow_folder']
                flow_files = arrakis_dict["flow_files"]

                """Check for arrakis folder"""
                if 'arrakis_folder' not in arrakis_dict.keys():
                    self.logger.warn('arrakis_folder not specified in config! setting to "/local_scratch"')
                    arrakis_dict['arrakis_folder'] = '/local_scratch'
                arrakis_folder = arrakis_dict['arrakis_folder'].replace('flow', 'arrakis').replace('FLOW', 'ARRAKIS')

                """Check that flow folder exists"""
                if not os.path.isdir(flow_folder):
                    self.logger.error(f'specified flow_folder {flow_folder} does not exist!')

                """Check that flow folder has a '/' at the end"""
                if flow_folder[-1] != '/':
                    flow_folder += '/'

                """Check that arrakis folder has a '/' at the end"""
                if arrakis_folder[-1] != '/':
                    arrakis_folder += '/'

                """Check that arrakis folder exists"""
                if not os.path.isdir(arrakis_folder):
                    os.makedirs(arrakis_folder)

                if isinstance(arrakis_dict["flow_files"], list):
                    """
                    If the flow_files parameter is a list, look through
                    the list and make sure each specified file actually exists
                    in the flow_folder.
                    """
                    flow_files = [
                        input_file for input_file in arrakis_dict["flow_files"]
                        if input_file not in arrakis_dict["skip_files"]
                    ]
                    for flow_file in flow_files:
                        if not os.path.isfile(flow_folder + flow_file):
                            self.logger.error(
                                f"specified file {flow_file} does not exist in directory {flow_folder}!"
                            )
                elif isinstance(arrakis_dict["flow_files"], str):
                    """
                    If the flow_files parameter is a string, check if its
                    the phrase 'all', and if so, recursively grab all h5
                    files in the flow_folder.

                    Otherwise, assume that the flow_files parameter is a
                    file extension, and search recursively for all files
                    with that extension.
                    """
                    if arrakis_dict["flow_files"] == "all":
                        self.logger.info(
                            f"searching {flow_folder} recursively for all .h5 FLOW files."
                        )
                        flow_files = [
                            os.path.basename(input_file) for input_file in glob.glob(
                                f"{flow_folder}*.h5", recursive=True
                            )
                            if 'FLOW' in input_file and input_file not in arrakis_dict["skip_files"]
                        ]
                    else:
                        try:
                            self.logger.info(
                                f'searching {flow_folder} recursively for all {arrakis_dict["flow_files"]} files.'
                            )
                            flow_files = [
                                os.path.basename(input_file) for input_file in glob.glob(
                                    f'{flow_folder}/*.{arrakis_dict["flow_files"]}',
                                    recursive=True,
                                )
                                if input_file not in arrakis_dict["skip_files"]
                            ]
                        except Exception as exception:
                            self.logger.error(
                                f'specified "flow_files" parameter: {arrakis_dict["flow_files"]} incompatible!'
                                + f" exception: {exception}"
                            )
                else:
                    self.logger.error(
                        f'specified "flow_files" parameter: {arrakis_dict["flow_files"]} incompatible!'
                    )
                self.flow_folder = flow_folder
                self.flow_files = flow_files
                self.arrakis_folder = arrakis_folder
                self.logger.info(f'setting arrakis_folder to {self.arrakis_folder}.')
                for ii in range(1, self.comm.size):
                    self.comm.send(self.flow_folder, dest=ii, tag=80)
                    self.comm.send(self.flow_files, dest=ii, tag=81)
                    self.comm.send(self.arrakis_folder, dest=ii, tag=82)
                self.logger.info(f'found {len(self.flow_files)} flow files for processing.')
            except Exception as exception:
                self.report_error(exception=exception)
            self.barrier()
        else:
            self.barrier()
            try:
                self.flow_folder = self.comm.recv(source=0, tag=80)
                self.flow_files = self.comm.recv(source=0, tag=81)
                self.arrakis_folder = self.comm.recv(source=0, tag=82)
            except Exception as exception:
                self.report_error(exception=exception)

    @profiler
    def reset_event_errors(self):
        self.event_errors = []
        self.plugin_errors = []
        self.event_plugin_errors = []
        self.event_plugin_exc_types = []
        self.event_plugin_exc_values = []
        self.event_plugin_exc_tracebacks = []
        self.event_plugin_line_numbers = []
        self.event_plugin_file_names = []

    @profiler
    def report_event_errors(
        self,
        file_name: str
    ):
        """
        Report any event/plugin errors that have occurred during
        this file.  These errors, if there are any, are reported
        and saved to the corresponding ARRAKIS file in the meta
        information.
        """
        event_errors = np.concatenate(self.comm.allgather(self.event_errors))
        plugin_errors = np.concatenate(self.comm.allgather(self.plugin_errors))
        event_plugin_errors = np.concatenate(self.comm.allgather(self.event_plugin_errors))
        event_plugin_line_numbers = np.concatenate(self.comm.allgather(self.event_plugin_line_numbers))
        event_plugin_file_names = np.concatenate(self.comm.allgather(self.event_plugin_file_names))
        if self.rank == 0:
            if len(event_errors) > 0:
                """Report event errors to log"""
                self.logger.warn(f"event/plugin errors occurred when processing file {file_name}")
                
                """Display error information"""
                for ii, event_error in enumerate(event_errors):
                    self.logger.warn(f"event: {event_error} / plugin: {plugin_errors[ii]}")
                    self.logger.warn(f" {' ':<{5}} exception: {event_plugin_errors[ii]}")
                    self.logger.warn(f" {' ':<{5}} line_number: {event_plugin_line_numbers[ii]}")
                    self.logger.warn(f" {' ':<{5}} file_name: {event_plugin_file_names[ii]}")
                errors = np.empty(len(event_errors), dtype=output_error_data_type)
                errors['event_id'] = event_errors.astype('i4')
                errors['plugin_name'] = plugin_errors.astype('S50')
                errors['plugin_error'] = event_plugin_errors.astype('S500')
                errors['plugin_line_number'] = event_plugin_line_numbers.astype('S50')
                errors['plugin_file_name'] = event_plugin_file_names.astype('S50')
                
                """Report error information to ARRAKIS file"""
                arrakis_file_name = file_name.replace('FLOW', 'ARRAKIS').replace('flow', 'arrakis')
                with h5py.File(self.arrakis_folder + arrakis_file_name, 'a') as arrakis_file:
                    original_size = arrakis_file['arrakis/errors'].shape[0]
                    additional_size = len(errors)
                    new_size = original_size + additional_size
                    """Insert the objects into the array"""
                    arrakis_file['arrakis/errors'].resize((new_size,))
                    arrakis_file['arrakis/errors'][original_size:new_size] = errors
        else:
            pass

    @profiler
    def run_begin_of_file(
        self,
        file_name: str
    ):
        """
        Plugins to be run at the beginning of the file,
        which act on the entire file.

        Args:
            file_name (_str_): name of the input flow file
        """
        pass

    @profiler
    def set_up_output_arrays(
        self,
        file_name: str
    ):
        """
        Construct output arrays for various data objects.
        This will take the flow file and create a new .h5 file
        replacing the words "FLOW" with "ARRAKIS" and "flow"
        with "arrakis".

        Args:
            file_name (str): name of the input flow file
        """

        """Determine the arrakis file name"""
        arrakis_file_name = file_name.replace('FLOW', 'ARRAKIS').replace('flow', 'arrakis')

        """Open the flow file and determine the output array shapes"""
        with h5py.File(self.flow_folder + file_name, 'r') as flow_file, \
             h5py.File(self.arrakis_folder + arrakis_file_name, 'a') as arrakis_file:
            """
            First we make the labels for the charge dataset.  These consist of six
            main labels that we wish to generate with various plugins:
                (1) topology - a descriptor of basic shapes
                (2) particle - the unique pdg of the particle
                (3) physics - a descriptor of local physics
                (4) unique_topology - unique labels for instances of topology
                (5) vertex - binary variable denoting whether a vertex is at this hit
                (6) tracklette_begin - binary variable denoting whether a track beginning is at this hit
                (7) tracklette_end - binary variable denoting whether a track end is at this hit
                (8) fragment_begin - binary variable denoting whether a fragment beginning is at this hit
                (9) fragment_end - binary variable denoting whether a fragment end is at this hit
                (10) shower_begin - binary variable denoting whether a shower beginning is at this hit

            These labels are assigned to each reconstructed charge hit and
            written in the corresponding ARRAKIS file.
            """
            charge_name = 'charge/calib_final_hits/data'
            charge_segment_name = 'charge_segment/calib_final_hits/data'
            charge_segments = flow_file['mc_truth/calib_final_hit_backtrack/data']['segment_id']
            charge_segments_fraction = flow_file['mc_truth/calib_final_hit_backtrack/data']['fraction']
            non_zero_charge_segments = [row[row != 0] for row in charge_segments]
            max_length = len(max(non_zero_charge_segments, key=len))
            num_charge = len(flow_file[charge_name]['x'])

            new_charge_segment_data_type = np.dtype([
                ('event_id', 'i4'),
                ('segment_id', 'i4', (max_length,)),
                ('segment_fraction', 'f4', (max_length,)),
                ('segment_distance', 'f4', (max_length,)),
                ('segment_parent', 'i4', (max_length,)),
                ('segment_start_process', 'i4', (max_length,)),
                ('segment_start_subprocess', 'i4', (max_length,)),
                ('topology', 'i4', (max_length,)),
                ('particle', 'i4', (max_length,)),
                ('physics', 'i4', (max_length,)),
                ('unique_topology', 'i4', (max_length,)),
                ('vertex', 'i4', (max_length,)),
                ('tracklette_begin', 'i4', (max_length,)),
                ('tracklette_end', 'i4', (max_length,)),
                ('fragment_begin', 'i4', (max_length,)),
                ('fragment_end', 'i4', (max_length,)),
                ('shower_begin', 'i4', (max_length,))
            ])

            new_charge_segment_data = np.full(num_charge, -1, dtype=new_charge_segment_data_type)
            new_charge_segment_data['segment_id'] = charge_segments[:, :max_length]
            new_charge_segment_data['segment_fraction'] = charge_segments_fraction[:, :max_length]
            new_charge_segment_data['segment_distance'][:, :max_length] = np.nan
            new_charge_segment_data['vertex'][:, :max_length] = 0
            new_charge_segment_data['tracklette_begin'][:, :max_length] = 0
            new_charge_segment_data['tracklette_end'][:, :max_length] = 0
            new_charge_segment_data['fragment_begin'][:, :max_length] = 0
            new_charge_segment_data['fragment_end'][:, :max_length] = 0
            new_charge_segment_data['shower_begin'][:, :max_length] = 0

            if charge_segment_name in arrakis_file:
                del arrakis_file[charge_segment_name]

            arrakis_file.create_dataset(charge_segment_name, data=new_charge_segment_data)

            new_charge_data_type = np.dtype([
                ('event_id', 'i4'),
                ('topology', 'i4'),
                ('particle', 'i4'),
                ('physics', 'i4'),
                ('unique_topology', 'i4'),
                ('vertex', 'i4'),
                ('tracklette_begin', 'i4'),
                ('tracklette_end', 'i4'),
                ('fragment_begin', 'i4'),
                ('fragment_end', 'i4'),
                ('shower_begin', 'i4')
            ])

            new_charge_data = np.full(num_charge, -1, dtype=new_charge_data_type)
            new_charge_data['vertex'][:] = 0
            new_charge_data['tracklette_begin'][:] = 0
            new_charge_data['tracklette_end'][:] = 0
            new_charge_data['fragment_begin'][:] = 0
            new_charge_data['fragment_end'][:] = 0
            new_charge_data['shower_begin'][:] = 0

            if charge_name in arrakis_file:
                del arrakis_file[charge_name]

            arrakis_file.create_dataset(charge_name, data=new_charge_data)

            """
            Second, we set up arrays for CAF-like reco objects which consist
            of the following types:
            """

            """Construct track data types"""
            tracklette_name = 'standard_record/tracklette'
            track_name = 'standard_record/track'
            if tracklette_name in arrakis_file:
                del arrakis_file[tracklette_name]
            if track_name in arrakis_file:
                del arrakis_file[track_name]

            arrakis_file.create_dataset(tracklette_name, shape=(0,), maxshape=(None,), dtype=tracklette_data_type)
            arrakis_file.create_dataset(track_name, shape=(0,), maxshape=(None,), dtype=track_data_type)

            """Construct shower data types"""
            fragment_name = 'standard_record/fragment'
            shower_name = 'standard_record/shower'
            if fragment_name in arrakis_file:
                del arrakis_file[fragment_name]
            if shower_name in arrakis_file:
                del arrakis_file[shower_name]

            arrakis_file.create_dataset(fragment_name, shape=(0,), maxshape=(None,), dtype=fragment_data_type)
            arrakis_file.create_dataset(shower_name, shape=(0,), maxshape=(None,), dtype=shower_data_type)

            """Construct blip data types"""
            blip_name = 'standard_record/blip'
            if blip_name in arrakis_file:
                del arrakis_file[blip_name]

            arrakis_file.create_dataset(blip_name, shape=(0,), maxshape=(None,), dtype=blip_data_type)

            """Construct particle data types"""
            particle_name = 'standard_record/particle'
            if particle_name in arrakis_file:
                del arrakis_file[particle_name]

            arrakis_file.create_dataset(particle_name, shape=(0,), maxshape=(None,), dtype=particle_data_type)

            """Construct interaction data types"""
            interaction_name = 'standard_record/interaction'
            if interaction_name in arrakis_file:
                del arrakis_file[interaction_name]

            arrakis_file.create_dataset(interaction_name, shape=(0,), maxshape=(None,), dtype=interaction_data_type)

            """Construct neutrino data types"""
            neutrino_name = 'standard_record/neutrino'
            if neutrino_name in arrakis_file:
                del arrakis_file[neutrino_name]

            arrakis_file.create_dataset(neutrino_name, shape=(0,), maxshape=(None,), dtype=neutrino_data_type)

            """Construct interaction event label types"""
            interaction_event_name = 'standard_record/interaction_event'
            if interaction_event_name in arrakis_file:
                del arrakis_file[interaction_event_name]

            arrakis_file.create_dataset(
                interaction_event_name, shape=(0,), maxshape=(None,), dtype=interaction_event_data_type
            )

            """Construct neutrino event label types"""
            neutrino_event_name = 'standard_record/neutrino_event'
            if neutrino_event_name in arrakis_file:
                del arrakis_file[neutrino_event_name]

            arrakis_file.create_dataset(neutrino_event_name, shape=(0,), maxshape=(None,), dtype=neutrino_event_data_type)

            """Construct output error data types"""
            output_error_name = 'arrakis/errors'
            if output_error_name in arrakis_file:
                del arrakis_file[output_error_name]

            arrakis_file.create_dataset(output_error_name, shape=(0,), maxshape=(None,), dtype=output_error_data_type)

            """Construct n-Ar inelastic data types"""
            nar_inelastic_name = 'standard_record/nar_inelastic'
            if nar_inelastic_name in arrakis_file:
                del arrakis_file[nar_inelastic_name]

            arrakis_file.create_dataset(nar_inelastic_name, shape=(0,), maxshape=(None,), dtype=nar_inelastic_data_type)
            
            """Construct undefined objects"""
            undefined_name = 'standard_record/undefined'
            if undefined_name in arrakis_file:
                del arrakis_file[undefined_name]
            
            arrakis_file.create_dataset(undefined_name, shape=(0,), maxshape=(None,), dtype=undefined_data_type)

    @profiler
    def distribute_tasks(
        self,
        file_name: str
    ):
        """
        Determine the indices which correspond to unique
        events in the FLOW output files and distribute
        those indices among the workers.

        The slowest part of this code is backtracking through hits->segments,
        which must be done in order to create the reverse map from segments->hits.
        This must be done since there is no way to determine what hits go
        with what events with only the calib_final_hits data alone.

        Args:
            file_name (_str_): _description_
        """
        with h5py.File(self.flow_folder + file_name, 'r', driver='mpio', comm=self.comm) as file:
            if self.rank == 0:
                """Collect event ids from mc_truth and charge/light data"""
                interactions_events = file['mc_truth/interactions/data']['event_id']
                segments_events = file['mc_truth/segments/data']['event_id']
                stack_events = file['mc_truth/stack/data']['event_id']
                trajectories_events = file['mc_truth/trajectories/data']['event_id']
                charge_segments = file['mc_truth/calib_final_hit_backtrack/data']['segment_id'].astype(int)
                charge_fraction = file['mc_truth/calib_final_hit_backtrack/data']['fraction']
                charge_fraction_mask = (charge_fraction == 0)
                charge_segments[charge_fraction_mask] = -1
                non_zero_charge_segments = [row[row != 0] for row in charge_fraction]
                max_length = len(max(non_zero_charge_segments, key=len))
                segments_ids = file['mc_truth/segments/data']['segment_id']

                """Determine unique events and distribute among workers"""
                unique_events = np.unique(trajectories_events)
                self.num_events = len(unique_events)
                self.distributed_events.clear()
                self.distributed_interactions_indices.clear()
                self.distributed_segments_indices.clear()
                self.distributed_stack_indices.clear()
                self.distributed_trajectories_indices.clear()
                self.distributed_charge_indices.clear()
                self.distributed_light_indices.clear()

                self.distributed_events = {i: [] for i in range(1, self.size)}
                self.distributed_interactions_indices = {i: [] for i in range(1, self.size)}
                self.distributed_segments_indices = {i: [] for i in range(1, self.size)}
                self.distributed_stack_indices = {i: [] for i in range(1, self.size)}
                self.distributed_trajectories_indices = {i: [] for i in range(1, self.size)}
                self.distributed_charge_indices = {i: [] for i in range(1, self.size)}
                self.distributed_light_indices = {i: [] for i in range(1, self.size)}

                """Determine indices for mc_truth and charge/light data"""
                for ii, event_id in enumerate(unique_events):
                    worker_rank = 1 + ii % (self.size - 1)  # Distribute round-robin among wor
                    self.distributed_events[worker_rank].append(event_id)
                    self.distributed_interactions_indices[worker_rank].append(
                        np.where(interactions_events == event_id)[0]
                    )
                    self.distributed_segments_indices[worker_rank].append(
                        np.where(segments_events == event_id)[0]
                    )
                    self.distributed_stack_indices[worker_rank].append(
                        np.where(stack_events == event_id)[0]
                    )
                    self.distributed_trajectories_indices[worker_rank].append(
                        np.where(trajectories_events == event_id)[0]
                    )
                    """For charge data we must backtrack through segments"""
                    hits_to_segments = np.any(
                        np.isin(
                            charge_segments[:, :max_length], segments_ids[(segments_events == event_id)]
                        ),
                        axis=1,
                    )
                    self.distributed_charge_indices[worker_rank].append(
                        hits_to_segments
                    )
                    """Likewise for light data, we must backtrack through segments"""
                    self.distributed_light_indices[worker_rank].append([])

                """Distribute indices for events among workers"""
                for ii in range(1, self.comm.size):
                    self.comm.send(self.distributed_events[ii], dest=ii, tag=2)
                    self.comm.send(self.distributed_interactions_indices[ii], dest=ii, tag=3)
                    self.comm.send(self.distributed_segments_indices[ii], dest=ii, tag=4)
                    self.comm.send(self.distributed_stack_indices[ii], dest=ii, tag=5)
                    self.comm.send(self.distributed_trajectories_indices[ii], dest=ii, tag=6)
                    self.comm.send(self.distributed_charge_indices[ii], dest=ii, tag=7)
                    self.comm.send(self.distributed_light_indices[ii], dest=ii, tag=8)
            else:
                self.distributed_events[self.rank] = self.comm.recv(source=0, tag=2)
                self.distributed_interactions_indices[self.rank] = self.comm.recv(source=0, tag=3)
                self.distributed_segments_indices[self.rank] = self.comm.recv(source=0, tag=4)
                self.distributed_stack_indices[self.rank] = self.comm.recv(source=0, tag=5)
                self.distributed_trajectories_indices[self.rank] = self.comm.recv(source=0, tag=6)
                self.distributed_charge_indices[self.rank] = self.comm.recv(source=0, tag=7)
                self.distributed_light_indices[self.rank] = self.comm.recv(source=0, tag=8)

    @profiler
    def process_events_master(
        self,
        flow_file: h5py.File,
        arrakis_file: h5py.File,
    ):
        """
        Special code run by the master node on an event,
        which occurs before any of the workers.  Any work done
        here should be by plugins which create data that is
        needed by the other plugins.

        Args:
            flow_file (_h5py.File_): input flow_file
        """
        """Clear event indices so that file closes properly!"""
        self.standard_record_objects = {
            'tracklette': [],
            'track': [],
            'fragment': [],
            'shower': [],
            'blip': [],
            'particle': [],
            'interaction': [],
            'neutrino': [],
            'interaction_event': [],
            'neutrino_event': [],
            'nar_inelastic': [],
            'undefined': [],
        }
        self.clear_indices()

    @profiler
    def process_events_worker(
        self,
        flow_file: h5py.File,
        arrakis_file: h5py.File,
    ):
        """
        This function loops over all plugins and hands off the flow_file,
        arrakis_file and the associated indices for each that correspond
        to this workers event.

        Args:
            flow_file (_h5py.File_): input flow_file
        """
        self.standard_record_objects = {
            'tracklette': [],
            'track': [],
            'fragment': [],
            'shower': [],
            'blip': [],
            'particle': [],
            'interaction': [],
            'neutrino': [],
            'interaction_event': [],
            'neutrino_event': [],
            'nar_inelastic': [],
            'undefined': [],
        }
        for ii, event in enumerate(self.distributed_events[self.rank]):
            """Grab event index information"""
            event_indices = {
                'interactions': self.distributed_interactions_indices[self.rank][ii],
                'segments': self.distributed_segments_indices[self.rank][ii],
                'stack': self.distributed_stack_indices[self.rank][ii],
                'trajectories': self.distributed_trajectories_indices[self.rank][ii],
                'charge': self.distributed_charge_indices[self.rank][ii],
                'light': self.distributed_light_indices[self.rank][ii]
            }
            self.event_products = {
                'tracklette': [],
                'track': [],
                'fragment': [],
                'shower': [],
                'blip': [],
                'particle': [],
                'interaction': [],
                'neutrino': [],
                'interaction_event': [],
                'neutrino_event': [],
                'nar_inelastic': [],
                'undefined': [],
            }
            """Iterate over plugins"""
            for plugin_name, plugin in self.plugins.items():
                try:
                    plugin.process_event(
                        event=event,
                        flow_file=flow_file,
                        arrakis_file=arrakis_file,
                        event_indices=event_indices,
                        event_products=self.event_products,
                    )
                except Exception as exception:
                    self.event_errors.append(event)
                    self.plugin_errors.append(plugin_name)
                    self.event_plugin_errors.append(exception)
                    exc_type, exc_value, exc_traceback = sys.exc_info()
                    if exc_traceback.tb_next:
                        exc_traceback = exc_traceback.tb_next
                    if exc_traceback.tb_next:
                        exc_traceback = exc_traceback.tb_next
                    self.event_plugin_exc_types.append(exc_type)
                    self.event_plugin_exc_values.append(exc_value)
                    self.event_plugin_exc_tracebacks.append(exc_traceback)
                    self.event_plugin_line_numbers.append(exc_traceback.tb_lineno)
                    self.event_plugin_file_names.append(exc_traceback.tb_frame.f_code.co_filename)
            for object in self.standard_record_objects.keys():
                self.standard_record_objects[object] += self.event_products[object]
        """Clear event indices so that file closes properly!"""
        self.clear_indices()

    @profiler
    def run_end_of_file(
        self,
        file_name: str
    ):
        """
        Plugins to be run on the entire file after all
        events have been evaluated within the file.

        One required task is for the master node to gather up all
        reco object information and append it to the arrays in the
        ARRAKIS file.

        Args:
            file_name (_type_): _description_
        """
        standard_record_objects = self.comm.allgather(self.standard_record_objects)

        if self.rank == 0:
            """Add those standard record objects to the ARRAKIS file"""
            arrakis_file_name = file_name.replace('FLOW', 'ARRAKIS').replace('flow', 'arrakis')

            """Open the flow file and determine the output array shapes"""
            with h5py.File(self.arrakis_folder + arrakis_file_name, 'a') as arrakis_file:
                for object in self.standard_record_objects.keys():
                    """Gather up all standard record objects"""
                    objects = [
                        standard_record_objects[ii][object]
                        for ii in range(self.size)
                        if len(standard_record_objects[ii][object]) != 0
                    ]
                    if len(objects):
                        """Flatten the objects array"""
                        flat_objects_list = [
                            item for sublist in objects for item in sublist
                        ]
                        objects = np.concatenate(flat_objects_list)
                        original_size = arrakis_file[f'standard_record/{object}'].shape[0]
                        additional_size = len(objects)
                        new_size = original_size + additional_size
                        """Insert the objects into the array"""
                        arrakis_file[f'standard_record/{object}'].resize((new_size,))
                        arrakis_file[f'standard_record/{object}'][original_size:new_size] = objects
        else:
            pass
    
    @profiler
    def run_undefined_check(
        self,
        file_name: str
    ):
        """
        We check for undefined labels in the arrakis output file and report 
        information about them.

        Args:
            file_name (str): _description_
        """
        arrakis_file_name = file_name.replace('FLOW', 'ARRAKIS').replace('flow', 'arrakis')
        with h5py.File(self.flow_folder + file_name, 'r') as flow_file, \
             h5py.File(self.arrakis_folder + arrakis_file_name, 'a') as arrakis_file:
            arrakis_charge = arrakis_file['charge/calib_final_hits/data']
            arrakis_topology = arrakis_charge['topology']
            arrakis_physics = arrakis_charge['physics']
            undefined_topology = np.where(arrakis_topology == -1)
            undefined_physics = np.where(arrakis_physics == -1)
            

    @profiler
    def run_end_of_arrakis(self):
        """
        Set of functions to be run at the end of the entire
        Arrakis job.  Some default operations are to create
        profiling plots.
        """
        if len(self.flow_files) != 0:
            self.generate_timing_and_memory_plots()
        if self.rank == 0:
            try:
                self.logger.info("Arrakis program ran successfully. Closing out.")
            except Exception as exception:
                self.report_error(exception=exception)
        else:
            pass

    @profiler
    def generate_timing_and_memory_plots(self):
        """Send over timing and memory information from each worker"""
        if self.rank == 0:
            try:
                collected_timings = {ii: self.comm.recv(source=ii, tag=100) for ii in range(1, self.size)}
                collected_memory = {ii: self.comm.recv(source=ii, tag=101) for ii in range(1, self.size)}
            except Exception as exception:
                self.report_error(exception=exception)
            self.barrier()
        else:
            try:
                self.comm.send(timing_manager.timings, dest=0, tag=100)
                self.comm.send(memory_manager.memory, dest=0, tag=101)
            except Exception as exception:
                self.report_error(exception=exception)
            self.barrier()

        """Generate timing and memory plots"""
        if self.rank == 0:
            try:
                collected_timings[0] = timing_manager.timings
                collected_memory[0] = memory_manager.memory

                timings = {}
                memory = {}
                for key in collected_timings.keys():
                    for name, value in collected_timings[key].items():
                        if name not in timings:
                            timings[name] = value
                            memory[name] = collected_memory[key][name]
                        else:
                            timings[name] += value
                            memory[name] += collected_memory[key][name]

                timing_averages = {}
                timing_stds = {}
                for item in timings.keys():
                    timing_averages[item] = np.mean(timings[item])
                    timing_stds[item] = np.std(timings[item])

                """Separate plugin timings from Arrakis timings"""
                arrakis_fig, arrakis_axs = plt.subplots(figsize=(15, 10))
                arrakis_box_values = []
                arrakis_labels = []
                for item in timings.keys():
                    if 'Plugin' in item:
                        continue
                    else:
                        arrakis_box_values.append(timings[item])
                        arrakis_axs.plot(
                            [],
                            [],
                            marker="",
                            linestyle="-",
                            label=f'{item}\n({timing_averages[item]:.2f} +/- {timing_stds[item]:.2f})',
                        )
                        arrakis_labels.append(item)
                arrakis_axs.boxplot(arrakis_box_values, vert=True, patch_artist=True, labels=arrakis_labels)
                arrakis_axs.set_ylabel(r"$\langle\Delta t\rangle$ (ms)")
                arrakis_axs.set_xticklabels(arrakis_labels, rotation=45, ha="right")
                arrakis_axs.set_yscale("log")
                plt.title(r"Arrakis $\langle\Delta t\rangle$ (ms) vs. function")
                plt.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
                plt.tight_layout()
                plt.savefig("/local_scratch/arrakis_nd_timing_avg.png")

                plugin_fig, plugin_axs = plt.subplots(figsize=(15, 10))
                plugin_box_values = []
                plugin_labels = []
                for item in timings.keys():
                    if 'Plugin' in item:
                        plugin_box_values.append(timings[item])
                        plugin_axs.plot(
                            [],
                            [],
                            marker="",
                            linestyle="-",
                            label=f'{item}\n({timing_averages[item]:.2f} +/- {timing_stds[item]:.2f})',
                        )
                        plugin_labels.append(item)
                    else:
                        continue
                plugin_axs.boxplot(plugin_box_values, vert=True, patch_artist=True, labels=plugin_labels)
                plugin_axs.set_ylabel(r"$\langle\Delta t\rangle$ (ms)")
                plugin_axs.set_xticklabels(plugin_labels, rotation=45, ha="right")
                plugin_axs.set_yscale("log")
                plt.title(r"Plugin $\langle\Delta t\rangle$ (ms) vs. function")
                plt.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
                plt.tight_layout()
                plt.savefig("/local_scratch/arrakis_nd_plugin_timing_avg.png")

                """Separate plugin memory from Arrakis memory"""
                memory_averages = {}
                memory_stds = {}
                for item in memory.keys():
                    memory_averages[item] = np.mean(memory[item])
                    memory_stds[item] = np.std(memory[item])

                arrakis_fig, arrakis_axs = plt.subplots(figsize=(15, 10))
                arrakis_box_values = []
                arrakis_labels = []
                for item in memory.keys():
                    if 'Plugin' in item:
                        continue
                    else:
                        arrakis_box_values.append(memory[item])
                        arrakis_axs.plot(
                            [],
                            [],
                            marker="",
                            linestyle="-",
                            label=f'{item}\n({memory_averages[item]:.2f} +/- {memory_stds[item]:.2f})',
                        )
                        arrakis_labels.append(item)
                arrakis_axs.boxplot(arrakis_box_values, vert=True, patch_artist=True, labels=arrakis_labels)
                arrakis_axs.set_ylabel(r"$\langle\Delta m\rangle$ (Mb)")
                arrakis_axs.set_xticklabels(arrakis_labels, rotation=45, ha="right")
                arrakis_axs.set_yscale("log")
                plt.title(r"Arrakis $\langle\Delta m\rangle$ (Mb) vs. function")
                plt.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
                plt.tight_layout()
                plt.savefig("/local_scratch/arrakis_nd_memory_avg.png")

                plugin_fig, plugin_axs = plt.subplots(figsize=(15, 10))
                plugin_box_values = []
                plugin_labels = []
                for item in memory.keys():
                    if 'Plugin' in item:
                        plugin_box_values.append(memory[item])
                        plugin_axs.plot(
                            [],
                            [],
                            marker="",
                            linestyle="-",
                            label=f'{item}\n({memory_averages[item]:.2f} +/- {memory_stds[item]:.2f})',
                        )
                        plugin_labels.append(item)
                    else:
                        continue
                plugin_axs.boxplot(plugin_box_values, vert=True, patch_artist=True, labels=plugin_labels)
                plugin_axs.set_ylabel(r"$\langle\Delta m\rangle$ (Mb)")
                plugin_axs.set_xticklabels(plugin_labels, rotation=45, ha="right")
                plugin_axs.set_yscale("log")
                plt.title(r"Plugin $\langle\Delta m\rangle$ (Mb) vs. function")
                plt.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
                plt.tight_layout()
                plt.savefig("/local_scratch/arrakis_nd_plugin_memory_avg.png")
            except Exception as exception:
                self.report_error(exception=exception)
            self.barrier()
        else:
            self.barrier()

    @profiler
    def run_arrakis_nd(self):
        """
        Main Arrakis program loop.

        Every function call here should be encapsulated in a try/except
        block to prevent hangups with comm.Barrier.

        I. Set up input files for run
        ~ Loop over each file
            a. Load file to determine unique events and indices for each array type
            b. Run whole file over file begin plugins
            c. Set up output arrays in output h5 file
            d. Send index, output information to each worker
            ~ Loop over each event in each worker
                i. Grab event related information and pass to plugins in order
                ii. Collect output information and add to output files
            e. Run whole file over file end plugins
            f. Run end of file functions
        II. Run end of program functions
        """
        self.barrier()

        """Set up input files"""
        try:
            self.set_up_input_files()
        except Exception as exception:
            self.report_error(exception=exception)
        self.barrier()

        """Set up progress bar"""
        if self.rank == 0:
            self.logger.info(f'running Arrakis with {self.size} workers.')
            self.progress_bar = tqdm(
                total=len(self.flow_files),
                ncols=100,
                colour='MAGENTA',
                leave=True,
            )
        self.barrier()

        """Loop over files and call master/worker methods for each."""
        for ii, file_name in enumerate(self.flow_files):
            """First run begin of file"""
            if self.rank == 0:
                try:
                    self.progress_bar.set_description_str(f'File [{ii+1}/{len(self.flow_files)}]')
                    self.run_begin_of_file(file_name)
                except Exception as exception:
                    self.report_error(exception=exception)
            self.barrier()

            """Reset event/plugin errors"""
            self.reset_event_errors()
            self.barrier()

            """Set up output arrays in ARRAKIS file"""
            if self.rank == 0:
                try:
                    self.set_up_output_arrays(file_name)
                except Exception as exception:
                    self.report_error(exception=exception)
            self.barrier()

            """Prepare indices for workers"""
            try:
                self.distribute_tasks(file_name)
            except Exception as exception:
                self.report_error(exception=exception)
            self.barrier()

            """Now process the file"""
            try:
                arrakis_file_name = file_name.replace('FLOW', 'ARRAKIS').replace('flow', 'arrakis')
                with h5py.File(self.flow_folder + file_name, 'r', driver='mpio', comm=self.comm) as flow_file, \
                     h5py.File(self.arrakis_folder + arrakis_file_name, 'r+', driver='mpio', comm=self.comm) as arrakis_file:
                    self.barrier()
                    if self.rank == 0:
                        """Process event in master node"""
                        try:
                            self.progress_bar.set_postfix_str(f'# Events: [{self.num_events}]')
                            self.process_events_master(flow_file, arrakis_file)
                        except Exception as exception:
                            self.report_error(exception=exception)
                        self.barrier()
                    else:
                        """Process event in worker node"""
                        self.barrier()
                        try:
                            self.process_events_worker(flow_file, arrakis_file)
                        except Exception as exception:
                            self.report_error(exception=exception)
                    """Ensure that changes are pushed to the arrakis file"""
                    self.barrier()
                    arrakis_file.flush()
            except Exception as exception:
                self.report_error(exception=exception)
            self.barrier()

            """Run end of file plugins"""
            try:
                self.run_end_of_file(file_name)
            except Exception as exception:
                self.report_error(exception=exception)
            self.barrier()
            
            """Run undefined check"""
            if self.rank == 0:
                try:
                    self.run_undefined_check(file_name)
                except Exception as exception:
                    self.report_error(exception=exception)
            self.barrier()
            
            """Update progress bar"""
            if self.rank == 0:
                try:
                    self.progress_bar.update(1)
                except Exception as exception:
                    self.report_error(exception=exception)
            self.barrier()

            """Report any event errors"""
            try:
                self.report_event_errors(file_name)
            except Exception as exception:
                self.report_error(exception=exception)
            self.barrier()

        """Run end of program functions"""
        if self.rank == 0:
            try:
                self.progress_bar.close()
            except Exception as exception:
                self.report_error(exception=exception)
        self.barrier()
        self.run_end_of_arrakis()
