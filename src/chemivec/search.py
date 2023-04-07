from typing import Union
import numpy as np
import pandas as pd
import multiprocessing as mp

from ._chemivec import _rxn_subsearch, _rxn_smarts_isok, _mol_subsearch
from .options import get_option, set_option, _process_n_jobs

def _convert_to_numpy(arr: Union[np.ndarray, pd.DataFrame, pd.Series, list]) -> np.ndarray:
    # Check the array type and convert everything to numpy
    if isinstance(arr, pd.DataFrame):
        if arr.shape[1] > 1:
            raise ValueError("Input dataframe has more than one column, "
                             "please use Series, 1D Numpy array or single column Dataframe")
        np_array = arr.squeeze().to_numpy()
    elif isinstance(arr, pd.Series):
        np_array = arr.to_numpy()
    elif isinstance(arr, list):
        np_array = np.array(arr, dtype=object)
    elif isinstance(arr, np.ndarray):
        np_array = arr
    else:
        raise ValueError("Input array can be from the following types: list, np.ndrray, pd.Series or pd.Dataframe,"
                         f"got {type(arr)} type instead")

    # now as arr is numpy ndarray, convert to C-CONTIGUOUS
    if not np_array.flags.c_contiguous:
        np_array = np.ascontiguousarray(arr)

    return np_array

def _convert_items_to_str(input: np.ndarray) -> np.ndarray:
    # check item type
    # first check 'np.str_' because it is subclass of 'str'
    out = input
    if isinstance(input[0], np.str_):
        out = input.astype(object)
    if not isinstance(input[0], str):
        raise ValueError(f"Input should be array of python or numpy strings, instead got array of {type(input[0])}")
    return out

def _validate_shape(arr: np.ndarray):
    # check array dims
    if arr.ndim != 1:
        raise ValueError(f"Multidimensional input arrays not allowed")
    if arr.shape[0] == 0:
        raise ValueError(f"Input array cannot be empty")
        # return np.array([], dtype=bool)

def _validate_input_arr(arr) -> np.ndarray:
    np_arr = _convert_to_numpy(arr)
    _validate_shape(np_arr)
    return _convert_items_to_str(np_arr)

def _validate_n_jobs(n_jobs):
    if n_jobs:
        return _process_n_jobs(n_jobs)
    else:
        return get_option("n_jobs")

def _validate_rxn_mode(mode: str):
    if mode != "DAYLIGHT-AAM":
        if not isinstance(mode, str):
            raise TypeError(f"mode expected to be str, instead {type(mode)} type received")


def _validate_mol_mode(mode: str):
    if not isinstance(mode, str):
        raise TypeError(f"mode expected to be str, instead {type(mode)} type received")


def _validate_rxn_query(query: str):
    # query smarts
    if not isinstance(query, str):
        raise TypeError(f"Query must be of string type, instead {type(query)} type received")
    if query is None or not query:
        raise ValueError(f"Query could not be empty or None")
    if not _rxn_smarts_isok(query):
        raise ValueError(f"Invalid reaction SMARTS:\n{query}")

def _validate_mol_query(query: str):
    # query smarts
    if not isinstance(query, str):
        raise TypeError(f"Query must be of string type, instead {type(query)} type received")
    if query is None or not query:
        raise ValueError(f"Query could not be empty or None")
    #TODO implement check mol query


def rxn_subsearch(arr: Union[np.ndarray, pd.DataFrame, pd.Series, list],
                  query: str = None,
                  n_jobs: Union[int, None] = None,
                  mode: str = "DAYLIGHT-AAM"
                  ) -> np.ndarray:
    """
    Vectorized reaction substructure search. Input SMILES array and query SMARTS. Both should
    be reactions, e.g. contains ">>" sign. By default, uses daylight atom-to-atom mapping rules:
    https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html (Section 4.6 Reaction Queries)
    If no atom mapping found in query - atom mappings are ignored. By default, uses all available cores
    for parallel computation. This number can be set globally `chemivec.set_option('n_jobs', 12)`

    Example:
        rxn_subsearch([ '[C:1]=O>>[C:1]O', 'C=O>>CO' ],
                  query = '[C:1]=O>>[C:1]O'
                  )
        output: array([ True, False])

        rxn_subsearch([ '[C:1]=O>>[C:1]O', 'C=O>>CO' ],
                  query='C=O>>CO'
                  )
        output: array([ True, True])

    :param arr: input array of reaction SMILES, supported inputs: np.ndarray, pd.DataFrame, pd.Series, list
    :param query: (str) reaction SMARTS
    :param mode: (str) by defaylt "DAYLIGHT-AAM"
    :param n_jobs: (int) number of threads or parallel computation, max by default
    :return: (np.ndarray[bool]) boolean result as numpy array
    """
    _validate_rxn_query(query)
    _validate_rxn_mode(mode)
    n_jobs_val = _validate_n_jobs(n_jobs)
    np_arr = _validate_input_arr(arr)

    return _rxn_subsearch(np_arr, query, mode, n_jobs_val)


def mol_subsearch(arr: Union[np.ndarray, pd.DataFrame, pd.Series, list],
                  query: str = None,
                  n_jobs: Union[int, None] = None,
                  mode: str = ""
                  ) -> np.ndarray:
    # query smarts
    _validate_mol_query(query)
    _validate_mol_mode(mode)
    n_jobs_val = _validate_n_jobs(n_jobs)
    np_arr = _validate_input_arr(arr)
    return _mol_subsearch(np_arr, query, mode, n_jobs_val)
