"""
Classes for storing ML memory information.
"""
import torch
import psutil
from matplotlib import pyplot as plt

from arrakis_nd.utils.common import process_functions


class MemoryTracker:
    """
    Internal class for recording memory information.
    """

    def __init__(
        self,
        name: str,
        gpu: bool = True,
    ):
        self.name = name
        self.gpu = gpu
        # initialized tensor
        self.gpu_memory_values = torch.empty(size=(0, 1), dtype=torch.float)
        self.cpu_memory_values = torch.empty(size=(0, 1), dtype=torch.float)
        if self.gpu:
            self.gpu_memory_values.cuda()
        self.gpu_memory_start = 0.0
        self.gpu_memory_end = 0.0
        self.cpu_memory_start = 0.0
        self.cpu_memory_end = 0.0

    def synchronize(self):
        torch.cuda.synchronize()

    def start(self):
        if self.gpu:
            self.gpu_memory_start = torch.cuda.memory_stats()[
                "allocated_bytes.all.allocated"
            ]
        self.cpu_memory_start = psutil.virtual_memory().used

    def end(self):
        if self.gpu:
            self.gpu_memory_end = torch.cuda.memory_stats()[
                "allocated_bytes.all.allocated"
            ]
            torch.cuda.synchronize()
            self.gpu_memory_values = torch.cat(
                (
                    self.gpu_memory_values,
                    torch.tensor([[self.gpu_memory_end - self.gpu_memory_start]]),
                )
            )
        self.cpu_memory_end = psutil.virtual_memory().used
        self.cpu_memory_values = torch.cat(
            (
                self.cpu_memory_values,
                torch.tensor([[self.cpu_memory_end - self.cpu_memory_start]]),
            )
        )


class MemoryTrackers:
    """
    Collection of memory_trackers for ML tasks.
    """

    def __init__(
        self,
        gpu: bool = True,
    ):
        self.gpu = gpu
        self.reset_trackers()

    def reset_trackers(self):
        self.memory_trackers = {
            f"{process}": MemoryTracker(f"{process}", gpu=self.gpu)
            for process in process_functions
        }

    def start(self, function):
        self.memory_trackers[function].start()

    def end(self, function):
        self.memory_trackers[function].end()

    def synchronize(self):
        torch.cuda.synchronize()

    def evaluate_run(self):
        cpu_averages = {}
        cpu_stds = {}
        for item in self.memory_trackers.keys():
            if len(self.memory_trackers[item].cpu_memory_values) == 0:
                continue
            if sum(self.memory_trackers[item].cpu_memory_values) == 0:
                continue
            temp_times = self.memory_trackers[item].cpu_memory_values.squeeze()
            cpu_averages[item] = temp_times.mean()
            cpu_stds[item] = temp_times.std()

        fig, axs = plt.subplots(figsize=(15, 10))
        box_values = torch.empty(
            size=(0, len(self.memory_trackers[item].cpu_memory_values.squeeze()))
        )
        labels = []
        for item in self.memory_trackers.keys():
            if len(self.memory_trackers[item].cpu_memory_values) == 0:
                continue
            if sum(self.memory_trackers[item].cpu_memory_values) == 0:
                continue
            temp_times = self.memory_trackers[item].cpu_memory_values.squeeze()
            box_values = torch.cat((box_values, temp_times.unsqueeze(0)), dim=0)
            axs.plot(
                [],
                [],
                marker="",
                linestyle="-",
                label=f'{item.replace("process_","")}\n({cpu_averages[item]:.2f} +/- {cpu_stds[item]:.2f})',
            )
            labels.append(f'{item.replace("process_","")}')
        if len(labels) != 0:
            axs.boxplot(box_values, vert=True, patch_artist=True, labels=labels)
            axs.set_ylabel(r"CPU $\langle\Delta m\rangle$ (bytes)")
            axs.set_xticklabels(labels, rotation=45, ha="right")
            axs.set_yscale("log")
            plt.title(r"CPU $\langle\Delta m\rangle$ (bytes) vs. function")
            plt.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
            plt.tight_layout()
            plt.savefig("/local_scratch/arrakis_nd_cpu_memory_avg.png")

        gpu_averages = {}
        gpu_stds = {}
        if not self.gpu:
            return
        for item in self.memory_trackers.keys():
            if len(self.memory_trackers[item].gpu_memory_values) == 0:
                continue
            if sum(self.memory_trackers[item].gpu_memory_values) == 0:
                continue
            temp_times = self.memory_trackers[item].gpu_memory_values.squeeze()
            gpu_averages[item] = temp_times.mean()
            gpu_stds[item] = temp_times.std()

        fig, axs = plt.subplots(figsize=(15, 10))
        box_values = torch.empty(
            size=(0, len(self.memory_trackers[item].gpu_memory_values.squeeze()))
        )
        labels = []
        for item in self.memory_trackers.keys():
            if len(self.memory_trackers[item].gpu_memory_values) == 0:
                continue
            if sum(self.memory_trackers[item].gpu_memory_values) == 0:
                continue
            temp_times = self.memory_trackers[item].gpu_memory_values.squeeze()
            box_values = torch.cat((box_values, temp_times.unsqueeze(0)), dim=0)
            axs.plot(
                [],
                [],
                marker="",
                linestyle="-",
                label=f'{item.replace("process_","")}\n({gpu_averages[item]:.2f} +/- {gpu_stds[item]:.2f})',
            )
            labels.append(f'{item.replace("process_","")}')
        if len(labels) != 0:
            axs.boxplot(box_values, vert=True, patch_artist=True, labels=labels)
            axs.set_ylabel(r"GPU $\langle\Delta m\rangle$ (bytes)")
            axs.set_xticklabels(labels, rotation=45, ha="right")
            axs.set_yscale("log")
            plt.title(r"GPU $\langle\Delta m\rangle$ (bytes) vs. function")
            plt.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
            plt.tight_layout()
            plt.savefig("/local_scratch/arrakis_nd_gpu_memory_avg.png")
