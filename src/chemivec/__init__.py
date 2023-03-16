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

def rxn_match(arr: Union[np.ndarray, pd.DataFrame, pd.Series, list],
            query_smarts: str = None,
            use_aam: bool = False) -> np.ndarray:

    # query smarts
    if query_smarts is None or not query_smarts:
        raise ValueError(f"query_smarts could not be empty or None, should be a SMARTS string")

    # atom-to-atom mapping mode
    aam_mode = ""
    if use_aam:
        aam_mode = "DAYLIGHT-AAM"

    arr = _convert_to_numpy(arr)

    # check item type
    # using direct type check because 'np.str_' is also instance of 'str'
    if type(arr[0] == str):
        return _rxn_match(arr, query_smarts, aam_mode)
    elif type(arr[0] == np.str_):
        return _rxn_match(arr.astype(object), query_smarts, aam_mode)
    else:
        raise ValueError(f"Input should be array of python or numpy strings, instead got array of {type(arr[0])}")
