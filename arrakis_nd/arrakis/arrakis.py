"""
"""
from mpi4py import MPI
import h5py
import os
import importlib.util
import traceback
import sys
import inspect
import glob
from tqdm import tqdm
import torch
import numpy as np
from datetime import datetime
from matplotlib import pyplot as plt
from importlib.metadata import version

from arrakis_nd.utils.logger import Logger
from arrakis_nd.utils.utils import (
    profiler,
    timing_manager,
    memory_manager
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
        except Exception as e:
            raise RuntimeError(f"unable to obtain MPI parameters: {e}")

        """Set input parameters and set up loggers"""
        try:
            if self.rank == 0:
                self.logger = Logger(f"master: {self.rank}", output="both")
            else:
                self.logger = Logger(f"worker: {self.rank}", output="both")
        except Exception as e:
            raise RuntimeError(f"unable to set up logging system: {e}")

        self.config = config
        self.meta = meta

        """Setting error status for this node"""
        self.error_status = None

        """Parse config"""
        try:
            self.parse_config()
        except Exception as e:
            self.error_status = e
        self.barrier()

        """Construct plugins"""
        self.plugins = {}
        try:
            self.construct_plugins()
        except Exception as e:
            self.error_status = e
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
    def barrier(self):
        errors = self.comm.allgather(self.error_status)
        if any(errors):
            if self.rank == 0:
                errors_count = sum(1 for error in errors if error is not None)
                errors_indices = [index for index, error in enumerate(errors) if error is not None]
                self.logger.critical(f"{errors_count} errors encountered in Arrakis program")
                for index in errors_indices:
                    self.logger.critical(f"error encountered in worker {index}: ")
                    self.logger.critical(errors[index])
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
                self.logger.info(f'############################ ARRAKIS  v. [{version("arrakis-nd")}] ############################')
                
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
            except Exception as e:
                self.error_status = e
            self.barrier()
        else:
            self.barrier()
            try:
                self.meta = self.comm.recv(source=0, tag=0)
            except Exception as e:
                self.error_status = e

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
            except Exception as e:
                self.logger.error(f'failed to construct plugin {plugin} with config {config}: {e}')

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
            running_output_products.append(self.plugins[plugin_name].output_product)

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
            except Exception as e:
                self.logger.warn(f'problem loading plugin from file: {plugin_file}: {e}')

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
        """
        if self.rank == 0:
            try:
                arrakis_dict = self.config['arrakis_nd']
                
                """Check for parameters"""
                if 'flow_folder' not in arrakis_dict.keys():
                    self.logger.error(f'flow_folder not specified in config!')
                if 'flow_files' not in arrakis_dict.keys():
                    self.logger.error(f'flow_files not specified in config!')
                    
                flow_folder = arrakis_dict['flow_folder']
                flow_files = arrakis_dict["flow_files"]
                
                """Check for arrakis folder"""
                if 'arrakis_folder' not in arrakis_dict.keys():
                    self.logger.warn(f'arrakis_folder not specified in config! setting to "/local_scratch"')
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
                        if not os.path.isfile(flow_file):
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
                                    f'{flow_folder}/{arrakis_dict["flow_files"]}',
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
            except Exception as e:
                self.error_status = e
            self.barrier()
        else:
            self.barrier()
            try:
                self.flow_folder = self.comm.recv(source=0, tag=80)
                self.flow_files = self.comm.recv(source=0, tag=81)
                self.arrakis_folder = self.comm.recv(source=0, tag=82)
            except Exception as e:
                self.error_status = e

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
                (2) physics_micro - a descriptor of local physics
                (3) physics_macro - a descriptor of larger scale physics objects
                (4) unique_topology - unique labels for instances of topology
                (5) unique_physics_micro - unique labels for instances of physics_micro
                (6) unique_physics_macro - unique labels for instances of physics_macro
                (7) vertex - binary variable denoting whether a vertex is at this hit
                (8) tracklette_begin - binary variable denoting whether a track beginning is at this hit
                (9) tracklette_end - binary variable denoting whether a track end is at this hit
                (10) fragment_begin - binary variable denoting whether a fragment beginning is at this hit
                (11) fragment_end - binary variable denoting whether a fragment end is at this hit
                (12) shower_begin - binary variable denoting whether a shower beginning is at this hit

            These labels are assigned to each reconstructed charge hit and
            written in the corresponding ARRAKIS file.
            """
            charge_name = 'charge/calib_final_hits/data'

            num_charge = len(flow_file[charge_name]['x'])

            new_charge_data_type = np.dtype([
                ('event_id', 'i4'),
                ('topology', 'i4'),
                ('physics_micro', 'i4'),
                ('physics_macro', 'i4'),
                ('unique_topology', 'i4'),
                ('unique_physics_micro', 'i4'),
                ('unique_physics_macro', 'i4'),
                ('vertex', 'bool'),
                ('tracklette_begin', 'bool'),
                ('tracklette_end', 'bool'),
                ('fragment_begin', 'bool'),
                ('fragment_end', 'bool'),
                ('shower_begin', 'bool')
            ])

            new_charge_data = np.full(num_charge, -1, dtype=new_charge_data_type)

            if charge_name in arrakis_file:
                del arrakis_file[charge_name]

            arrakis_file.create_dataset(charge_name, data=new_charge_data)

            """
            Second, we set up arrays for CAF-like reco objects which consist
            of the following types:
                (1) tracklettes - pieces of contiguous track like objects.
                    (a) start (x,y,z) - track start point
                    (b) end (x,y,z) - track end point
                    (c) dir (u_x,u_y,u_z) - estimate of track direction taken from start point
                    (d) enddir (u_x,u_y,u_z) - estimate of track direction takend from end point
                    (e) Evis - visible energy in voxels corresponding to this track
                    (f) qual - reco specific quality metric
                    (g) len_gcm2 - track length in g/cm2
                    (h) len_cm - track length in centimeter
                    (i) E - track energy estimate in MeV
                    (j) truth - associated true particle in flow file (if relevant) (i.e. track_id)
                    (k) truthOverlap - fractional overlap between this track and true particle
                (2) tracks - collections of tracklettes which define a complete track.
                (3) fragments - pieces of contiguous shower like objects.
                    (a) start (x,y,z) - shower start point
                    (b) dir (u_x,u_y,u_z) - shower direction
                    (c) Evis - visible energy in voxels corresponding to this shower
                    (d) qual - reco specific quality metric
                    (e) len_gcm2 - track length in g/cm2
                    (f) len_cm - track length in centimeter
                    (g) E - track energy estimate in MeV
                    (h) truth - associated true particle in flow file (if relevant) (i.e. track_id)
                    (i) truthOverlap - fractional overlap between this shower and true particle
                (4) showers - collections of fragments which define a complete shower.
                (5) blips - collections of points associated to blip-like activity.
                    (a) start (x,y,z) - start position of this blip object
                    (b) Evis - visible energy in voxels corresponding to this blip
                    (c) E - reconstructed energy (GeV)
                    (d) bliphyp - hypothesis for this blip's identity
                    (e) truth - associated true particle in flow file (if relevant) (i.e track_id)
                    (f) truthOverlap - fractional overlap between this blip and true particle
                (6) particles - particles within an event that are associated to tracks/showers/blips.
                    (a) primary - is this reco particle a primary one (i.e. eminates directly from vertex)?
                    (b) pdg - pdg code inferred for this particle
                    (c) tgtA - atomic number of nucleus this particle was reconstructed in
                    (d) score - PID score for this particle
                    (e) E - reconstructed energy (GeV)
                    (f) E_method - method used to determine energy for the particle
                    (g) p (p_x,p_y,p_z) - reconstructed momentum for this particle
                    (h) start (x,y,z) - reconstructed start point of this particle
                    (i) end (x,y,z) - reconstructed end point of this particle
                    (j) contained - contained in LAr TPC?
                    (k) truth - associated true particle in flow file (if relevant) (i.e. track_id)
                    (l) truthOverlap - fractional overlap between this reco particle and true particle
                (7) interactions - collections of tracks/showers/blips which define interactions of interest.
                    (a) vtx (x,y,z) - reconstructed vertex location
                    (b) dir (u_x,u_y,u_z) - hypothesis for this interaction's parent particle direction
                        (b.i) lngtrk - direction using longest track
                        (b.ii) heshw - direction using highest energy shower
                    (c) nuhyp - hypothesis for this interaction's neutrino identity
                        (c.i) isnubar
                        (c.ii) nue
                        (c.iii) numu
                        (c.iv) nutau
                        (c.v) nc
                        (c.vi) protons0
                        (c.vii) protons1
                        (c.viii) protons2
                        (c.ix) protonsN
                        (c.x) chgpi0
                        (c.xi) chgpi1
                        (c.xii) chgpi2
                        (c.xiii) chgpiN
                        (c.xiv) pizero0
                        (c.xv) pizero1
                        (c.xvi) pizero2
                        (c.xvii) pizeroN
                        (c.xviii) neutron0
                        (c.xix) neutron1
                        (c.xx) neutron2
                        (c.xxi) neutronN
                    (d) Enu - hypothesis for this interaction's neutrino energy
                    (e) part - collections of reconstructed particles
                    (f) truth - indicies of true_interactions in flow file (if relevant)
                    (g) truthOverlap - fractional overlap between this reco interaction and each true interaction
            """

            """Construct track data types"""
            tracklette_name = 'standard_record/tracklette'
            track_name = 'standard_record/track'
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
            track_data_type = np.dtype([
                ('event_id', 'i4'),
                ('track_id', 'i4'),
                ('tracklette_ids', 'i4', (1, 20)),
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
            if tracklette_name in arrakis_file:
                del arrakis_file[tracklette_name]
            if track_name in arrakis_file:
                del arrakis_file[track_name]

            arrakis_file.create_dataset(tracklette_name, shape=(0,), maxshape=(None,), dtype=tracklette_data_type)
            arrakis_file.create_dataset(track_name, shape=(0,), maxshape=(None,), dtype=track_data_type)

            """Construct shower data types"""
            fragment_name = 'standard_record/fragment'
            shower_name = 'standard_record/shower'
            fragment_data_type = np.dtype([
                ('event_id', 'i4'),
                ('fragment_id', 'i4'),
                ('start', 'f4', (1, 3)),
                ('dir', 'f4', (1, 3)),
                ('Evis', 'f4'),
                ('qual', 'f4'),
                ('len_gcm2', 'f4'),
                ('len_cm', 'f4'),
                ('E', 'f4'),
                ('truth', 'i4', (1, 20)),
                ('truthOverlap', 'f4', (1, 20)),
            ])
            shower_data_type = np.dtype([
                ('event_id', 'i4'),
                ('shower_id', 'i4'),
                ('fragment_ids', 'i4', (1, 20)),
                ('start', 'f4', (1, 3)),
                ('dir', 'f4', (1, 3)),
                ('Evis', 'f4'),
                ('qual', 'f4'),
                ('len_gcm2', 'f4'),
                ('len_cm', 'f4'),
                ('E', 'f4'),
                ('truth', 'i4', (1, 20)),
                ('truthOverlap', 'f4', (1, 20)),
            ])
            if fragment_name in arrakis_file:
                del arrakis_file[fragment_name]
            if shower_name in arrakis_file:
                del arrakis_file[shower_name]

            arrakis_file.create_dataset(fragment_name, shape=(0,), maxshape=(None,), dtype=fragment_data_type)
            arrakis_file.create_dataset(shower_name, shape=(0,), maxshape=(None,), dtype=shower_data_type)

            """Construct blip data types"""
            blip_name = 'standard_record/blip'
            blip_data_type = np.dtype([
                ('event_id', 'i4'),
                ('blip_id', 'i4'),
                ('start', 'f4', (1, 3)),
                ('Evis', 'f4'),
                ('E', 'f4'),
                ('bliphyp', 'i4'),
                ('truth', 'i4', (1, 20)),
                ('truthOverlap', 'f4', (1, 20)),
            ])
            if blip_name in arrakis_file:
                del arrakis_file[blip_name]

            arrakis_file.create_dataset(blip_name, shape=(0,), maxshape=(None,), dtype=blip_data_type)

            """Construct particle data types"""
            particle_name = 'standard_record/particle'
            particle_data_type = np.dtype([
                ('event_id', 'i4'),
                ('particle_id', 'i4'),
                ('primary', 'bool'),
                ('pdg', 'i4'),
                ('tgtA', 'i4'),
                ('score', 'f4'),
                ('E', 'f4'),
                ('E_method', 'i4'),
                ('p', 'f4', (1, 3)),
                ('start', 'f4', (1, 3)),
                ('end', 'f4', (1, 3)),
                ('contained', 'bool'),
                ('truth', 'i4', (1, 20)),
                ('truthOverlap', 'f4', (1, 20)),
            ])
            if particle_name in arrakis_file:
                del arrakis_file[particle_name]

            arrakis_file.create_dataset(particle_name, shape=(0,), maxshape=(None,), dtype=particle_data_type)

            """Construct interaction data types"""
            interaction_name = 'standard_record/interaction'
            interaction_data_type = np.dtype([
                ('event_id', 'i4'),
                ('interaction_id', 'i4'),
                ('vtx', 'f4', (1, 3)),
                ('dir_lngtrk', 'f4', (1, 3)),
                ('dir_heshw', 'f4', (1, 3)),
                ('nuhyp', 'f4', (1, 20)),
                ('Enu', 'f4'),
                ('E_method', 'i4'),
                ('part', 'i4', (1, 20)),
                ('truth', 'i4', (1, 20)),
                ('truthOverlap', 'f4', (1, 20)),
            ])
            if interaction_name in arrakis_file:
                del arrakis_file[interaction_name]

            arrakis_file.create_dataset(interaction_name, shape=(0,), maxshape=(None,), dtype=interaction_data_type)

    @profiler
    def distribute_tasks(
        self,
        file_name: str
    ):
        """
        Determine the indices which correspond to unique
        events in the FLOW output files and distribute
        those indices among the workers.

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
                charge_segments = file['mc_truth/calib_final_hit_backtrack/data']['segment_id']
                non_zero_charge_segments = [row[row != 0] for row in charge_segments]
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
                for i, event_id in enumerate(unique_events):
                    worker_rank = 1 + i % (self.size - 1)  # Distribute round-robin among wor
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
                    self.distributed_charge_indices[worker_rank].append(
                        np.any(
                            np.isin(
                                charge_segments[:, :max_length], segments_ids[(segments_events == event_id)]
                            ),
                            axis=1,
                        )
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
        pass

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
            event_products = {}
            """Iterate over plugins"""
            for plugin_name, plugin in self.plugins.items():
                plugin.process_event(
                    event=event,
                    flow_file=flow_file,
                    arrakis_file=arrakis_file,
                    event_indices=event_indices,
                    event_products=event_products,
                )

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
        pass

    @profiler
    def run_end_of_arrakis(self):
        """
        Set of functions to be run at the end of the entire
        Arrakis job.  Some default operations are to create
        profiling plots.
        """
        self.generate_timing_and_memory_plots()
        if self.rank == 0:
            try:
                self.logger.info("Arrakis program ran successfully. Closing out.")
            except Exception as e:
                self.error_status = e
        else:
            pass

    @profiler
    def generate_timing_and_memory_plots(self):
        """Send over timing and memory information from each worker"""
        if self.rank == 0:
            try:
                collected_timings = {ii: self.comm.recv(source=ii, tag=100) for ii in range(1, self.size)}
                collected_memory = {ii: self.comm.recv(source=ii, tag=101) for ii in range(1, self.size)}
            except Exception as e:
                self.error_status = e
            self.barrier()
        else:
            try:
                self.comm.send(timing_manager.timings, dest=0, tag=100)
                self.comm.send(memory_manager.memory, dest=0, tag=101)
            except Exception as e:
                self.error_status = e
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
            except Exception as e:
                self.error_status = e
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
        except Exception as e:
            self.error_status = e
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
                except Exception as e:
                    self.error_status = e
            self.barrier()

            """Set up output arrays in ARRAKIS file"""
            if self.rank == 0:
                try:
                    self.set_up_output_arrays(file_name)
                except Exception as e:
                    self.error_status = e
            self.barrier()

            """Prepare indices for workers"""
            try:
                self.distribute_tasks(file_name)
            except Exception as e:
                self.error_status = e
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
                        except Exception as e:
                            self.error_status = e
                        self.barrier()
                    else:
                        """Process event in worker node"""
                        self.barrier()
                        try:
                            self.process_events_worker(flow_file, arrakis_file)
                        except Exception as e:
                            self.error_status = e
            except Exception as e:
                self.error_status = e
            self.barrier()

            """Run end of file plugins"""
            if self.rank == 0:
                try:
                    self.run_end_of_file(file_name)
                    self.progress_bar.update(1)
                except Exception as e:
                    self.error_status = e
            self.barrier()

        """Run end of program functions"""
        if self.rank == 0:
            self.progress_bar.close()
        self.run_end_of_arrakis()
