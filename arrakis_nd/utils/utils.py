"""
Utilities for ArrakisND
"""
import numpy as np
import h5py


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


@np.vectorize
def fill_with_minus_one(arr):
    return np.full_like(arr, -1)


def read_hdf5(path):
    keys = []
    with h5py.File(path, 'r') as f:     # open file
        f.visit(keys.append)            # append all keys to list
        for key in keys:
            print(f[key].name)
