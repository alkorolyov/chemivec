from typing import Union
import numpy as np
import pandas as pd
from ._chemivec import _rxn_match


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

    # query smarts
    if query_smarts is None or not query_smarts:
        raise ValueError(f"query_smarts could not be empty or None, should be a SMARTS string")

    arr = _convert_to_numpy(arr)

    # check item type
    # first check 'np.str_' because it is subclass of 'str'
    if isinstance(arr[0], np.str_):
        return _rxn_match(arr.astype(object), query_smarts, aam_mode)
    elif isinstance(arr[0], str):
        return _rxn_match(arr.astype(object), query_smarts, aam_mode)

    raise ValueError(f"Input should be array of python or numpy strings, instead got array of {type(arr[0])}")
