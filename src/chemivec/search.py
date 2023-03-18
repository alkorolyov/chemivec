from typing import Union
import numpy as np
import pandas as pd
import multiprocessing as mp

from ._chemivec import _rxn_match
from .options import get_option, set_option, _check_num_cores

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


def rxn_match(arr: Union[np.ndarray, pd.DataFrame, pd.Series, list],
              query_smarts: str = None,
              aam_mode: str = "DAYLIGHT-AAM",
              num_cores: Union[int, None] = None) -> np.ndarray:
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

    :param num_cores:
    :param arr: input array of reaction SMILES, supported inputs: np.ndarray, pd.DataFrame, pd.Series, list
    :param query_smarts: (str) reaction SMARTS
    :param aam_mode: (str) by defaylt "DAYLIGHT-AAM"
    :return: (np.ndarray[bool]) boolean result as numpy array
    """
    # query smarts
    if query_smarts is None or not query_smarts:
        raise ValueError(f"query_smarts could not be empty or None, should be a SMARTS string")
    # TODO check smarts

    # num_cores
    if num_cores:
        num_cores = _check_num_cores(int(num_cores))
    else:
        num_cores = get_option("num_cores")

    arr = _convert_to_numpy(arr)

    # check item type
    # first check 'np.str_' because it is subclass of 'str'
    if isinstance(arr[0], np.str_):
        return _rxn_match(arr.astype(object), query_smarts, aam_mode, num_cores)
    elif isinstance(arr[0], str):
        return _rxn_match(arr, query_smarts, aam_mode, num_cores)

    raise ValueError(f"Input should be array of python or numpy strings, instead got array of {type(arr[0])}")



