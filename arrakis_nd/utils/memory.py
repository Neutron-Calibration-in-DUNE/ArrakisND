"""
Classes for storing ML memory information.
"""
import torch
from matplotlib import pyplot as plt

from arrakis_nd.utils.common import process_functions


class MemoryTracker:
    """
    Internal class for recording memory information.
    """
    def __init__(
        self,
        name:   str,
        gpu:    bool = True,
    ):
        self.name = name
        self.gpu = gpu
        # initialized tensor
        self.memory_values = torch.empty(size=(0, 1), dtype=torch.float)
        if self.gpu:
            self.memory_values.cuda()
            self.memory_start = 0.0
            self.memory_end = 0.0
            self.start = self._start_cuda
            self.end = self._end_cuda
        else:
            self.memory_start = 0.0
            self.memory_end = 0.0
            self.start = self._start_cpu
            self.end = self._end_cpu

    def synchronize(self):
        torch.cuda.synchronize()

    def _start_cuda(self):
        self.memory_start = torch.cuda.memory_stats()['allocated_bytes.all.allocated']

    def _start_cpu(self):
        self.memory_start = 0.0

    def _end_cuda(self):
        self.memory_end = torch.cuda.memory_stats()['allocated_bytes.all.allocated']
        torch.cuda.synchronize()
        self.memory_values = torch.cat(
            (self.memory_values, torch.tensor([[self.memory_end - self.memory_start]]))
        )

    def _end_cpu(self):
        self.memory_end = 0.0


class MemoryTrackers:
    """
    Collection of memory_trackers for ML tasks.
    """
    def __init__(
        self,
        gpu:    bool = True,
    ):
        self.gpu = gpu
        self.reset_trackers()

    def reset_trackers(self):
        self.memory_trackers = {
            f'process_{process}': MemoryTracker(f'process_{process}', gpu=self.gpu) for process in process_functions
        }

    def synchronize(self):
        torch.cuda.synchronize()

    def evaluate_run(self):
        averages = {}
        stds = {}
        for item in self.memory_trackers.keys():
            if len(self.memory_trackers[item].memory_values) == 0:
                continue
            temp_times = self.memory_trackers[item].memory_values.squeeze()
            averages[item] = temp_times.mean()
            stds[item] = temp_times.std()

        fig, axs = plt.subplots(figsize=(15, 10))
        box_values = torch.empty(size=(0, len(self.memory_trackers[item].memory_values.squeeze())))
        labels = []
        for item in self.memory_trackers.keys():
            if len(self.memory_trackers[item].memory_values) == 0:
                continue
            temp_times = self.memory_trackers[item].memory_values.squeeze()
            box_values = torch.cat((box_values, temp_times.unsqueeze(0)), dim=0)
            axs.plot(
                [], [],
                marker='', linestyle='-',
                label=f'{item.replace("process_","")}\n({averages[item]:.2f} +/- {stds[item]:.2f})'
            )
            labels.append(f'{item.replace("process_","")}')
        axs.boxplot(
            box_values,
            vert=True,
            patch_artist=True,
            labels=labels
        )
        axs.set_ylabel(r"$\langle\Delta m\rangle$ (bytes)")
        axs.set_xticklabels(labels, rotation=45, ha='right')
        axs.set_yscale('log')
        plt.title(r"$\langle\Delta m\rangle$ (bytes) vs. function")
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.tight_layout()
        plt.savefig("/local_scratch/arrakis_nd_memory_avg.png")

