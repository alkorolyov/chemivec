from typing import Union
import numpy as np
import pandas as pd
import re

from ._chemivec import _set_option, _get_option

OPTION_TYPE = {
    "num_cores": int
}

def _set_int_option(name: str, value: Union[str, int]):
    # integer values
    if isinstance(value, str):
        if not re.match("^[0-9]+$", value):
            raise ValueError(f"'{name}' value '{value}' is not a number")
        _set_option(name, value)
    elif isinstance(value, int):
        _set_option(name, str(value))


def _set_str_option(name: str, value: str):
    if not isinstance(value, str):
        raise ValueError(f"'{name}' value '{value}' is not a string")
    _set_option(name, value)


def set_option(name: str, value: Union[str, int]):
    """
    Set global option in Chemivec module.
    :param name: (str) option name
    :param value: option value
    :return:
    """
    if not name in OPTION_TYPE:
        raise ValueError(f"Option `{name}` not allowed")
    if OPTION_TYPE[name] == int:
        _set_int_option(name, value)
    if OPTION_TYPE[name] == str:
        _set_str_option(name, value)

def get_option(name: str):
    """
    Get global option from Chemivec module by name
    :param name: option name
    :return:
    """
    if not name in OPTION_TYPE:
        raise ValueError(f"Option `{name}` not found")
    if OPTION_TYPE[name] == int:
        return int(_get_option(name))
    elif OPTION_TYPE[name] == str:
        return str(_get_option(name))
