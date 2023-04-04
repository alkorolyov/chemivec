from typing import Union
import numpy as np
import pandas as pd
import multiprocessing as mp
import re

from ._chemivec import _set_option, _get_option

SUPPORTED_OPTIONS = {
    "n_jobs": int
}

# _set_option(name: str, value)
"""
Low level function to set option.
Option names passed as Unicode Python strings. Option values are passed
as Python object, with correct type. Type checking is done on the python side 
"""
# _get_option(name: str)
"""
Low level function to get option on module level.
Options names are always passed as strings. 
"""



def _set_str_option(name: str, value: str):
    if not isinstance(value, str):
        raise TypeError(f"'{name}' value '{value}' is not a string")
    _set_option(name, value)


def _process_n_jobs(value) -> int:
    if isinstance(value, float):
        raise TypeError(f"float type not allowed, int or string expected")
    value = int(value)
    if value < 0:
        raise ValueError(f"Negative 'n_jobs' not allowed")
    elif value == 0 or value > mp.cpu_count():
        value = mp.cpu_count()
    return value


def set_option(name: str, value: Union[str, int]):
    """
    Set global option in Chemivec module.
    :param name: (str) option name
    :param value: option value
    :return:
    """
    if not name in SUPPORTED_OPTIONS:
        raise ValueError(f"Option `{name}` not supported, must be one of : {SUPPORTED_OPTIONS.keys()}")

    # n_jobs
    if name == "n_jobs":
        value = _process_n_jobs(value)

    _set_option(name, value)


def get_option(name: str):
    """
    Get global option from Chemivec module by name
    :param name: (str) option name
    :return: option value
    """
    if not name in SUPPORTED_OPTIONS:
        raise ValueError(f"Option `{name}` not supported, must be one of : {SUPPORTED_OPTIONS.keys()}")
    if SUPPORTED_OPTIONS[name] == int:
        return int(_get_option(name))
    elif SUPPORTED_OPTIONS[name] == str:
        return str(_get_option(name))
