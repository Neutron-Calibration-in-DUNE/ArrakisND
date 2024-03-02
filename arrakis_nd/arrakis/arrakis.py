"""
"""
from mpi4py import MPI
import h5py
import os
import glob
from tqdm import tqdm
import torch
import numpy as np

from arrakis_nd.utils.logger import Logger
from arrakis_nd.utils.config import ConfigParser
from arrakis_nd.utils.utils import (
    profiler,
    timing_manager,
    memory_manager
)


class Arrakis:
    """
    """
    @profiler
    def __init__(
        self,
        config: dict = {},
        meta:   dict = {},
        number_of_files:    int = -1,
    ):
        """_summary_

        Args:
            config (dict): _description_
            meta (dict): _description_
            number_of_files (int): _description_
        """
        """
        Parse config parameters
        Configure plugins
        """

        """Get mpi communication parameters"""
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

        """Set input parameters and set up loggers"""
        if self.rank == 0:
            self.logger = Logger(f"master: {self.rank}", output="both")
        else:
            self.logger = Logger(f"worker: {self.rank}", output="both")
        self.config = config
        self.meta = meta

        """Parse config"""
        self.parse_config()
        self.comm.Barrier()

        """Construct plugins"""
        self.construct_plugins()
        self.comm.Barrier()

        self.flow_files = []
        """Distributed events and indices for various arrays"""
        self.distributed_events = {}
        self.distributed_interactions_indices = {}
        self.distributed_segments_indices = {}
        self.distributed_stack_indices = {}
        self.distributed_trajectories_indices = {}
        self.distributed_charge_indices = {}

    @profiler
    def parse_config(self):
        """_summary_
        """
        if self.rank == 0:
            if "arrakis_nd" not in self.config.keys():
                self.logger.error("arrakis_nd section not in config!")
            arrakis_nd_config = self.config["arrakis_nd"]

            system_info = self.logger.get_system_info()
            for key, value in system_info.items():
                self.logger.info(f"system_info - {key}: {value}")

            if "verbose" in arrakis_nd_config:
                if not isinstance(arrakis_nd_config["verbose"], bool):
                    self.logger.error(
                        f'"arrakis_nd:verbose" must be of type bool, but got {type(arrakis_nd_config["verbose"])}!'
                    )
                self.meta["verbose"] = arrakis_nd_config["verbose"]
            else:
                self.meta["verbose"] = False

            # Eventually we will want to check that the order of the arrakis_nds makes sense,
            # and that the data products are compatible and available for the different modes.

            # check for devices
            if "gpu" not in arrakis_nd_config.keys():
                self.logger.warn('"arrakis_nd:gpu" not specified in config!')
                gpu = None
            else:
                gpu = arrakis_nd_config["gpu"]
            if "gpu_device" not in arrakis_nd_config.keys():
                self.logger.warn('"arrakis_nd:gpu_device" not specified in config!')
                gpu_device = None
            else:
                gpu_device = arrakis_nd_config["gpu_device"]

            if torch.cuda.is_available():
                self.logger.info("CUDA is available with devices:")
                for ii in range(torch.cuda.device_count()):
                    device_properties = torch.cuda.get_device_properties(ii)
                    cuda_stats = f"name: {device_properties.name}, "
                    cuda_stats += (
                        f"compute: {device_properties.major}.{device_properties.minor}, "
                    )
                    cuda_stats += f"memory: {device_properties.total_memory}"
                    self.logger.info(f" -- device: {ii} - " + cuda_stats)

            # set gpu settings
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
            for ii in range(1, self.comm.size):
                self.comm.send(self.meta, dest=ii, tag=0)
        else:
            self.meta = self.comm.recv(source=0, tag=0)

    @profiler
    def set_up_input_files(self):
        """_summary_
        """
        if self.rank == 0:
            arrakis_dict = self.config['arrakis_nd']
            flow_folder = arrakis_dict['flow_folder']
            flow_files = arrakis_dict["flow_files"]
            if isinstance(arrakis_dict["flow_files"], list):
                flow_files = [
                    input_file for input_file in arrakis_dict["flow_files"]
                    if input_file not in arrakis_dict["skip_files"]
                ]
            elif isinstance(arrakis_dict["flow_files"], str):
                if arrakis_dict["flow_files"] == "all":
                    self.logger.info(
                        f"searching {flow_folder} recursively for all .h5 files."
                    )
                    flow_files = [
                        os.path.basename(input_file) for input_file in glob.glob(
                            f"{flow_folder}*.h5", recursive=True
                        )
                        if input_file not in arrakis_dict["skip_files"]
                    ]
                else:
                    try:
                        self.logger.info(
                            f'searching {flow_folder} recursively for all {arrakis_dict["flow_files"]} files.'
                        )
                        flow_files = [
                            os.path.basename(input_file) for input_file in glob.glob(
                                f'{flow_folder}{arrakis_dict["flow_files"]}',
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
            self.flow_files = [flow_folder + flow_file for flow_file in flow_files]
            for ii in range(1, self.comm.size):
                self.comm.send(self.flow_files, dest=ii, tag=1)
        else:
            self.flow_files = self.comm.recv(source=0, tag=1)

    @profiler
    def construct_plugins(self):
        """_summary_
        """
        pass

    @profiler
    def run_begin_of_file(
        self,
        file_name:  str
    ):
        """_summary_

        Args:
            file_name (_type_): _description_
        """
        pass

    @profiler
    def set_up_output_arrays(
        self,
        file_name:  str
    ):
        """_summary_

        Args:
            file_name (str): _description_
        """
        with h5py.File(file_name, 'r+') as f:
            charge = f['charge/calib_final_hits/data']
            # Define the compound data type with named fields
            data_type = np.dtype([
                ('topology', 'i4'),  # Example: 4-byte integer
                ('physics_micro', 'i4'),  # Example: 8-byte float
                ('physics_meso', 'i4')  # Example: String of up to 10 characters
            ])
            N = len(charge['x'])
            self.data = np.full(N, -1, dtype=data_type)
            if 'classes' in f:
                del f['classes']  # Delete the existing dataset
            f.create_dataset('classes', data=self.data)

    @profiler
    def distribute_tasks(
        self,
        file_name:  str
    ):
        """_summary_

        Args:
            file_name (_type_): _description_
        """
        with h5py.File(file_name, 'r+', driver='mpio', comm=self.comm) as file:
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
                self.distributed_events.clear()
                self.distributed_interactions_indices.clear()
                self.distributed_segments_indices.clear()
                self.distributed_stack_indices.clear()
                self.distributed_trajectories_indices.clear()
                self.distributed_charge_indices.clear()

                self.distributed_events = {i: [] for i in range(1, self.size)}
                self.distributed_interactions_indices = {i: [] for i in range(1, self.size)}
                self.distributed_segments_indices = {i: [] for i in range(1, self.size)}
                self.distributed_stack_indices = {i: [] for i in range(1, self.size)}
                self.distributed_trajectories_indices = {i: [] for i in range(1, self.size)}
                self.distributed_charge_indices = {i: [] for i in range(1, self.size)}

                """Determine indices for mc_truth and charge/light data"""
                for i, event_id in enumerate(unique_events):
                    worker_rank = 1 + i % (self.size - 1)  # Distribute round-robin among workers
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

                """Distribute indices for events among workers"""
                for ii in range(1, self.comm.size):
                    self.comm.send(self.distributed_events[ii], dest=ii, tag=2)
                    self.comm.send(self.distributed_interactions_indices[ii], dest=ii, tag=3)
                    self.comm.send(self.distributed_segments_indices[ii], dest=ii, tag=4)
                    self.comm.send(self.distributed_stack_indices[ii], dest=ii, tag=5)
                    self.comm.send(self.distributed_trajectories_indices[ii], dest=ii, tag=6)
                    self.comm.send(self.distributed_charge_indices[ii], dest=ii, tag=7)
            else:
                self.distributed_events[self.rank] = self.comm.recv(source=0, tag=2)
                self.distributed_interactions_indices[self.rank] = self.comm.recv(source=0, tag=3)
                self.distributed_segments_indices[self.rank] = self.comm.recv(source=0, tag=4)
                self.distributed_stack_indices[self.rank] = self.comm.recv(source=0, tag=5)
                self.distributed_trajectories_indices[self.rank] = self.comm.recv(source=0, tag=6)
                self.distributed_charge_indices[self.rank] = self.comm.recv(source=0, tag=7)

    @profiler
    def process_events_master(
        self,
        file:   h5py.File
    ):
        """_summary_

        Args:
            file (_type_): _description_
        """
        pass

    @profiler
    def process_events_worker(
        self,
        file:   h5py.File
    ):
        """_summary_

        Args:
            file (_type_): _description_
        """
        progress_bar = tqdm(
            enumerate(self.distributed_events[self.rank], 0),
            total=len(self.distributed_events[self.rank]),
        )
        for ii, event in progress_bar:
            trajectories = file['mc_truth/trajectories/data'][self.distributed_trajectories_indices[self.rank][ii]]
            topology = file['classes']['topology']
            progress_bar.set_description(f"Worker {self.rank}:")
            progress_bar.set_postfix_str(f"")

    @profiler
    def run_end_of_file(
        self,
        file_name:  str
    ):
        """_summary_

        Args:
            file_name (_type_): _description_
        """
        pass

    @profiler
    def run_end_of_arrakis(self):
        """_summary_
        """
        for key in timing_manager.timings:
            self.logger.info(
                f"func: {key} - " +
                f"avg_time: {np.mean(timing_manager.timings[key])} - " +
                f"avg_memory: {np.mean(memory_manager.memory[key])}"
            )

    @profiler
    def run_arrakis_nd(self):
        """_summary_
        """
        """
        Loop over each file
            Load file to determine unique events and indices for each array type
            Run whole file over file begin plugins
            Set up output arrays in output h5 file
            Send index, output information to each worker
            Loop over each event in each worker
                Grab event related information and pass to plugins in order
                Collect output information and add to output files
            Run whole file over file end plugins
            Run end of file functions
        Run end of program functions
        """
        self.comm.Barrier()
        """Set up input files"""
        self.set_up_input_files()
        self.comm.Barrier()
        """Loop over files and call master/worker methods for each."""
        for file_name in self.flow_files:
            """First set up output data and prepare indices for workers"""
            if self.rank == 0:
                self.run_begin_of_file(file_name)
                self.set_up_output_arrays(file_name)
            self.distribute_tasks(file_name)
            self.comm.Barrier()
            """Now process the file"""
            with h5py.File(file_name, 'r+', driver='mpio', comm=self.comm) as file:
                self.comm.Barrier()
                if self.rank == 0:
                    self.process_events_master(file)
                    self.comm.Barrier()
                else:
                    self.comm.Barrier()
                    self.process_events_worker(file)
            self.comm.Barrier()
            """Run end of file plugins"""
            if self.rank == 0:
                self.run_end_of_file(file_name)
            self.comm.Barrier()
        """Run end of program functions"""
        self.run_end_of_arrakis()
