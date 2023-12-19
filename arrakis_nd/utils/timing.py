"""
Classes for storing ML timing information.
"""
import torch
import time
from matplotlib import pyplot as plt

from arrakis_nd.utils.common import process_functions


class Timer:
    """
    Internal class for recording timing information.
    """
    def __init__(
        self,
        name:   str,
        gpu:    bool = True,
    ):
        self.name = name
        self.gpu = gpu
        # initialized tensor
        self.timer_values = torch.empty(size=(0, 1), dtype=torch.float)
        if self.gpu:
            self.timer_values.cuda()
            self.timer_start = torch.cuda.Event(enable_timing=True)
            self.timer_end = torch.cuda.Event(enable_timing=True)
            self.start = self._start_cuda
            self.end = self._end_cuda
        else:
            self.timer_start = 0
            self.timer_end = 0
            self.start = self._start_cpu
            self.end = self._end_cpu

    def synchronize(self):
        torch.cuda.synchronize()

    def _start_cuda(self):
        self.timer_start.record()

    def _start_cpu(self):
        self.timer_start = time.time()

    def _end_cuda(self):
        self.timer_end.record()
        torch.cuda.synchronize()
        self.timer_values = torch.cat(
            (self.timer_values, torch.tensor([[self.timer_start.elapsed_time(self.timer_end)]]))
        )

    def _end_cpu(self):
        self.timer_end = time.time()


class Timers:
    """
    Collection of timers for ML tasks.
    """
    def __init__(
        self,
        gpu:    bool = True,
    ):
        self.gpu = gpu
        self.reset_timers()

    def reset_timers(self):
        self.timers = {
            f'{process}': Timer(f'{process}', gpu=self.gpu) for process in process_functions
        }

    def synchronize(self):
        torch.cuda.synchronize()

    def evaluate_run(self):
        averages = {}
        stds = {}
        for item in self.timers.keys():
            if len(self.timers[item].timer_values) == 0:
                continue
            temp_times = self.timers[item].timer_values.squeeze()
            averages[item] = temp_times.mean()
            stds[item] = temp_times.std()

        fig, axs = plt.subplots(figsize=(15, 10))
        box_values = torch.empty(size=(0, len(self.timers[item].timer_values.squeeze())))
        labels = []
        for item in self.timers.keys():
            if len(self.timers[item].timer_values) == 0:
                continue
            temp_times = self.timers[item].timer_values.squeeze()
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
        axs.set_ylabel(r"$\langle\Delta t\rangle$ (ms)")
        axs.set_xticklabels(labels, rotation=45, ha='right')
        axs.set_yscale('log')
        plt.title(r"$\langle\Delta t\rangle$ (ms) vs. function")
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.tight_layout()
        plt.savefig("/local_scratch/arrakis_nd_timing_avg.png")
