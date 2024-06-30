"""
Utilities for ArrakisND
"""
from collections import defaultdict
import functools
from datetime import datetime
import psutil
import time
import os
import numpy as np
from scipy.interpolate import splev


def get_datetime():
    """
    Returns a datetime in the format of:
    YYYY.MM.DD.HH.MM.SS

    Returns:
        _type_: _description_
    """
    time = datetime.now()
    now = f"{time.year}.{time.month}.{time.day}.{time.hour}.{time.minute}.{time.second}"
    return now


def integrand(t, tck):
    dxdt, dydt, dzdt = splev(t, tck, der=1)
    return np.sqrt(dxdt**2 + dydt**2 + dzdt**2)


def fiducialized_vertex(vert_pos):
    lar_x = [
        (-63.9273, -3.0652),
        (3.0652, 63.9273)
    ]
    lar_y = [(-62.055, 62.055)]
    lar_z = [
        (-64.51125, -2.48125),
        (2.48125, 64.51125)
    ]
    for i in lar_z:
        if vert_pos[2] < i[0] and vert_pos[2] > i[1]:
            return False
    for i in lar_x:
        if vert_pos[0] < i[0] and vert_pos[0] > i[1]:
            return False
    for i in lar_y:
        if vert_pos[1] < i[0] and vert_pos[1] > i[1]:
            return False
    return True


class ResetableIterator:
    """_summary_
    """
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
    """_summary_

    Args:
        original_list (_type_): _description_
        remove_list (_type_): _description_

    Returns:
        _type_: _description_
    """
    return list(filter(lambda x: x not in remove_list, original_list))


class TimingManager:
    """_summary_
    """
    def __init__(self):
        self.timings = defaultdict(list)

    def record_timing(self, func_name, elapsed):
        self.timings[func_name].append(elapsed)


class MemoryManager:
    """_summary_
    """
    def __init__(self):
        self.memory = defaultdict(list)

    def record_memory(self, func_name, elapsed):
        self.memory[func_name].append(elapsed)


# Global instance to store timings
timing_manager = TimingManager()
# Global instance to store profiles
memory_manager = MemoryManager()


def profiler(func):
    """_summary_

    Args:
        func (_type_): _description_

    Returns:
        _type_: _description_
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        process = psutil.Process(os.getpid())
        mem_before = process.memory_info().rss / 1024 ** 2  # Convert bytes to MB

        start_time = time.time()
        result = func(*args, **kwargs)

        end_time = time.time()
        elapsed_time = end_time - start_time

        mem_after = process.memory_info().rss / 1024 ** 2  # Convert bytes to MB
        mem_used = mem_after - mem_before

        # Capture the memory profile output
        try:
            memory_manager.record_memory(func.__qualname__, mem_used)
            timing_manager.record_timing(func.__qualname__, elapsed_time)
        except Exception:
            memory_manager.record_memory(func.__name__, mem_used)
            timing_manager.record_timing(func.__name__, elapsed_time)

        return result
    return wrapper
