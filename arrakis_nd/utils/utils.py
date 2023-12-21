"""
Utilities for ArrakisND
"""


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


def remove_sublist(
    original_list,
    remove_list
):
    return list(filter(lambda x: x not in remove_list, original_list))
