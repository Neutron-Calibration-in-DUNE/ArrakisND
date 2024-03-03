"""
"""
from mpi4py import MPI
import h5py
import os
import glob
from tqdm import tqdm
import torch
import numpy as np
from matplotlib import pyplot as plt

from arrakis_nd.utils.logger import Logger
from arrakis_nd.utils.utils import (
    profiler,
    timing_manager,
    memory_manager
)


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
            self.error_status = str(e)
        self.barrier()

        """Construct plugins"""
        try:
            self.construct_plugins()
        except Exception as e:
            self.error_status = str(e)
        self.barrier()

        self.flow_files = []
        """Distributed events and indices for various arrays"""
        self.distributed_events = {}
        self.distributed_interactions_indices = {}
        self.distributed_segments_indices = {}
        self.distributed_stack_indices = {}
        self.distributed_trajectories_indices = {}
        self.distributed_charge_indices = {}

    @profiler
    def barrier(self):
        errors = self.comm.allgather(self.error_status)
        if any(errors):
            self.logger.error(f"errors encountered in Arrakis run: {errors}")
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
                """Check for main arrakis_nd parameters"""
                if "arrakis_nd" not in self.config.keys():
                    self.logger.error("arrakis_nd section not in config!")
                arrakis_nd_config = self.config["arrakis_nd"]

                """Try to grab system info and display to the logger"""
                system_info = self.logger.get_system_info()
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

                """Check for GPU devices and configure them"""
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

                """See if any CUDA devices are available on the system"""
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
                self.error_status = str(e)
            self.barrier()
        else:
            self.barrier()
            try:
                self.meta = self.comm.recv(source=0, tag=0)
            except Exception as e:
                self.error_status = str(e)

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
                flow_folder = arrakis_dict['flow_folder']
                flow_files = arrakis_dict["flow_files"]

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
            except Exception as e:
                self.error_status = str(e)
            self.barrier()
        else:
            self.barrier()
            try:
                self.flow_files = self.comm.recv(source=0, tag=1)
            except Exception as e:
                self.error_status = str(e)

    @profiler
    def construct_plugins(self):
        """
        Construct all plugin objects.
        """
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
            file_name (_str_): _description_
        """
        pass

    @profiler
    def set_up_output_arrays(
        self,
        file_name: str
    ):
        """
        Construct output arrays for various data objects.

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
        file_name: str
    ):
        """
        Determine the indices which correspond to unique
        events in the FLOW output files and distribute
        those indices among the workers.

        Args:
            file_name (_str_): _description_
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
        file: h5py.File
    ):
        """
        Special code run by the master node on an event,
        which occurs before any of the workers.  Any work done
        here should be by plugins which create data that is
        needed by the other plugins.

        Args:
            file (_h5py.File_): _description_
        """
        pass

    @profiler
    def process_events_worker(
        self,
        file: h5py.File
    ):
        """_summary_

        Args:
            file (_h5py.File_): _description_
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
        file_name: str
    ):
        """
        Plugins to be run on the entire file after all 
        events have been evaluated.

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
                self.error_status = str(e)
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
                self.error_status = str(e)
            self.barrier()
        else:
            try:
                self.comm.send(timing_manager.timings, dest=0, tag=100)
                self.comm.send(memory_manager.memory, dest=0, tag=101)
            except Exception as e:
                self.error_status = str(e)
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

                fig, axs = plt.subplots(figsize=(15, 10))
                box_values = []
                labels = []
                for item in timings.keys():
                    box_values.append(timings[item])
                    axs.plot(
                        [],
                        [],
                        marker="",
                        linestyle="-",
                        label=f'{item}\n({timing_averages[item]:.2f} +/- {timing_stds[item]:.2f})',
                    )
                    labels.append(item)
                axs.boxplot(box_values, vert=True, patch_artist=True, labels=labels)
                axs.set_ylabel(r"$\langle\Delta t\rangle$ (ms)")
                axs.set_xticklabels(labels, rotation=45, ha="right")
                axs.set_yscale("log")
                plt.title(r"$\langle\Delta t\rangle$ (ms) vs. function")
                plt.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
                plt.tight_layout()
                plt.savefig("/local_scratch/arrakis_nd_timing_avg.png")

                memory_averages = {}
                memory_stds = {}
                for item in memory.keys():
                    memory_averages[item] = np.mean(memory[item])
                    memory_stds[item] = np.std(memory[item])

                fig, axs = plt.subplots(figsize=(15, 10))
                box_values = []
                labels = []
                for item in memory.keys():
                    box_values.append(memory[item])
                    axs.plot(
                        [],
                        [],
                        marker="",
                        linestyle="-",
                        label=f'{item}\n({memory_averages[item]:.2f} +/- {memory_stds[item]:.2f})',
                    )
                    labels.append(item)
                axs.boxplot(box_values, vert=True, patch_artist=True, labels=labels)
                axs.set_ylabel(r"$\langle\Delta m\rangle$ (Mb)")
                axs.set_xticklabels(labels, rotation=45, ha="right")
                axs.set_yscale("log")
                plt.title(r"$\langle\Delta m\rangle$ (Mb) vs. function")
                plt.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
                plt.tight_layout()
                plt.savefig("/local_scratch/arrakis_nd_memory_avg.png")
            except Exception as e:
                self.error_status = str(e)
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
            self.error_status = str(e)
        self.barrier()

        """Loop over files and call master/worker methods for each."""
        for file_name in self.flow_files:
            """First set up output data and prepare indices for workers"""
            if self.rank == 0:
                try:
                    self.run_begin_of_file(file_name)
                except Exception as e:
                    self.error_status = str(e)
                try:
                    self.set_up_output_arrays(file_name)
                except Exception as e:
                    self.error_status = str(e)

            """Prepare indices for workers"""
            try:
                self.distribute_tasks(file_name)
            except Exception as e:
                self.error_status = str(e)
            self.barrier()

            """Now process the file"""
            try:
                with h5py.File(file_name, 'r+', driver='mpio', comm=self.comm) as file:
                    self.barrier()
                    if self.rank == 0:
                        """Process event in master node"""
                        try:
                            self.process_events_master(file)
                        except Exception as e:
                            self.error_status = str(e)
                        self.barrier()
                    else:
                        """Process event in worker node"""
                        self.barrier()
                        try:
                            self.process_events_worker(file)
                        except Exception as e:
                            self.error_status = str(e)
            except Exception as e:
                self.error_status = str(e)
            self.barrier()

            """Run end of file plugins"""
            if self.rank == 0:
                try:
                    self.run_end_of_file(file_name)
                except Exception as e:
                    self.error_status = str(e)
            self.barrier()

        """Run end of program functions"""
        self.run_end_of_arrakis()
