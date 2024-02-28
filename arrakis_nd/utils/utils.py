"""
Utilities for ArrakisND
"""
import numpy as np
from collections import defaultdict
import functools
import psutil
import time
import os


class ResetableIterator:
    def __init__(self):
        self._value = 0

    def __iter__(self):
        return self

    def __next__(self):
        self._value += 1
        return self._value

    def __call__(self):
        return self.__next__()

    def reset(self):
        self._value = 0


def remove_sublist(original_list, remove_list):
    return list(filter(lambda x: x not in remove_list, original_list))


# A class to store and manage timings
class TimingManager:
    def __init__(self):
        self.timings = defaultdict(list)

    def record_timing(self, func_name, elapsed):
        self.timings[func_name].append(elapsed)


# A class to store and manage timings
class MemoryManager:
    def __init__(self):
        self.memory = defaultdict(list)

    def record_memory(self, func_name, elapsed):
        self.memory[func_name].append(elapsed)


# Global instance to store timings
timing_manager = TimingManager()
# Global instance to store profiles
memory_manager = MemoryManager()


def profiler(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Redirect stdout to a StringIO object to capture the memory_profiler output
        process = psutil.Process(os.getpid())
        mem_before = process.memory_info().rss / 1024 ** 2  # Convert bytes to MB

        start_time = time.time()
        result = func(*args, **kwargs)

        end_time = time.time()
        elapsed_time = end_time - start_time

        # Get the process after the function execution
        mem_after = process.memory_info().rss / 1024 ** 2  # Convert bytes to MB
        mem_used = mem_after - mem_before

        # Capture the memory profile output
        try:
            memory_manager.record_memory(func.__qualname__, mem_used)
            timing_manager.record_timing(func.__qualname__, elapsed_time)
        finally:
            memory_manager.record_memory(func.__name__, mem_used)
            timing_manager.record_timing(func.__name__, elapsed_time)

        return result
    return wrapper
