from typing import Union
import numpy as np
import pandas as pd
import multiprocessing as mp
import re

from ._chemivec import _set_option, _get_option

SUPPORTED_OPTIONS = {
    "num_cores": int
}

def _set_str_option(name: str, value: str):
    if not isinstance(value, str):
        raise TypeError(f"'{name}' value '{value}' is not a string")
    _set_option(name, value)


def _check_num_cores(value: int) -> str:
    if isinstance(value, float):
        raise TypeError(f"float type not allowed, int or string expected")
    value = int(value)
    if value < 0:
        raise ValueError(f"Negative 'num_cores' not allowed")
    elif value == 0 or value > mp.cpu_count():
        value = mp.cpu_count()
    return str(value)



def set_option(name: str, value: Union[str, int]):
    """
    Set global option in Chemivec module.
    :param name: (str) option name
    :param value: option value
    :return:
    """
    if not name in SUPPORTED_OPTIONS:
        raise ValueError(f"Option `{name}` not supported, must be one of : {SUPPORTED_OPTIONS.keys()}")

    # num_cores
    if name == "num_cores":
        value = _check_num_cores(value)
        _set_option(name, value)


def get_option(name: str):
    """
    Get global option from Chemivec module by name
    :param name: option name
    :return:
    """
    if not name in SUPPORTED_OPTIONS:
        raise ValueError(f"Option `{name}` not supported, must be one of : {SUPPORTED_OPTIONS.keys()}")
    if SUPPORTED_OPTIONS[name] == int:
        return int(_get_option(name))
    elif SUPPORTED_OPTIONS[name] == str:
        return str(_get_option(name))
