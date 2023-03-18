from typing import Union
import numpy as np
import pandas as pd
import re
from ._chemivec import _rxn_match, _set_option, _get_option


def _convert_to_numpy(arr: Union[np.ndarray, pd.DataFrame, pd.Series, list]) -> np.ndarray:
    # Check the array type and convert everything to numpy
    if isinstance(arr, pd.DataFrame):
        if arr.shape[1] > 1:
            raise ValueError("Input dataframe has more than one column, "
                             "please use Series, 1D Numpy array or single column Dataframe")
        arr = arr.squeeze().to_numpy()
    elif isinstance(arr, pd.Series):
        arr = arr.to_numpy()
    elif isinstance(arr, list):
        arr = np.array(arr, dtype=object)
    elif isinstance(arr, np.ndarray):
        pass
    else:
        raise ValueError("Input array can be from the following types: list, np.ndrray, pd.Series or pd.Dataframe,"
                         f"got {type(arr)} type instead")
    return arr

def rxn_match(arr: Union[np.ndarray, pd.DataFrame, pd.Series, list], query_smarts: str = None,
              aam_mode: str = "DAYLIGHT-AAM") -> np.ndarray:
    """
    Vectorized reaction substructure search. Input SMILES array and query SMARTS. Both should
    be reactions, e.g. contains ">>" sign. By default uses daylight atom-to-atom mapping rules:
    https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html (Section 4.6 Reaction Queries)
    If no atom mapping found in query - atom mappings are ignored
    Example:
        rxn_match([ '[C:1]=O>>[C:1]O', 'C=O>>CO' ],
                  query_smarts = '[C:1]=O>>[C:1]O'
                  )
        output: array([ True, False])

        rxn_match([ '[C:1]=O>>[C:1]O', 'C=O>>CO' ],
                  query_smarts='C=O>>CO'
                  )
        output: array([ True, True])

    :param arr: input array of reaction SMILES, supported inputs: np.ndarray, pd.DataFrame, pd.Series, list
    :param query_smarts: (str) reaction SMARTS
    :param aam_mode: (str) by defaylt "DAYLIGHT-AAM"
    :return: (np.ndarray[bool]) boolean result as numpy array
    """
    # query smarts
    if query_smarts is None or not query_smarts:
        raise ValueError(f"query_smarts could not be empty or None, should be a SMARTS string")

    arr = _convert_to_numpy(arr)

    num_cores = 1

    # check item type
    # first check 'np.str_' because it is subclass of 'str'
    if isinstance(arr[0], np.str_):
        return _rxn_match(arr.astype(object), query_smarts, aam_mode, num_cores)
    elif isinstance(arr[0], str):
        return _rxn_match(arr.astype(object), query_smarts, aam_mode, num_cores)

    raise ValueError(f"Input should be array of python or numpy strings, instead got array of {type(arr[0])}")

""" Options setters and getters """

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


