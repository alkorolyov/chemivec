import numpy as np
import pandas as pd
import pytest
import os
import sys
import multiprocessing as mp

try:
    # skbuild
    # package available through pip install
    from chemivec import mol_subsearch, set_option, get_option
except ModuleNotFoundError:
    # clion build
    # run command `pytest ./tests` from root project folder
    sys.path.append(os.getcwd())
    from src.chemivec import mol_subsearch, set_option, get_option

MAX_CPU_COUNT = mp.cpu_count()

""" TEST FIRST ARGUMENT - array """

def test_numpy_npstr():
    arr = np.array(['C=O', 'CO'])
    query = "C=O"
    res = mol_subsearch(arr, query)
    assert res[0]
    assert not res[1]
    assert res.dtype == np.bool_


def test_numpy_pystr():
    arr = np.array(['C=O', 'CO'], dtype=object)
    query = "C=O"
    res = mol_subsearch(arr, query)
    assert res[0]
    assert not res[1]


def test_pylist():
    arr = ['C=O', 'CO']
    query = "C=O"
    res = mol_subsearch(arr, query)
    assert res[0]
    assert not res[1]


def test_pandas_df():
    arr = pd.DataFrame(['C=O', 'CO'])
    query = "C=O"
    res = mol_subsearch(arr, query)
    assert res[0]
    assert not res[1]


def test_pandas_series():
    arr = pd.Series(['C=O', 'CO'])
    query = "C=O"
    res = mol_subsearch(arr, query)
    assert res[0]
    assert not res[1]


def test_multi_dim():
    arr = np.array([[1, 2], [2, 3]])
    query = "C=O"
    with pytest.raises(ValueError, match="Multidimensional input arrays not allowed"):
        mol_subsearch(arr, query)


def test_empty_array():
    arr = np.array([])
    query = "C=O"
    res = mol_subsearch(arr, query)
    assert res.shape[0] == 0
    assert res.dtype == np.bool_


# def test_bad_reaction_smiles():
#     arr = np.array(['C]>>'])
#     query = "C>>"
#     res = mol_subsearch(arr, query)
#     assert not res[0]
#
#
# def test_long_reaction_smiles():
#     arr = np.array(['C'+']>>'])
#     query = "C>>"
#     res = mol_subsearch(arr, query)
#     assert not res[0]
#
#
# """ TEST SECOND ARGUMENT - query """
#
# def test_bad_query():
#     arr = np.array(['C'])
#     query = "[C>>"
#     with pytest.raises(ValueError, match="Invalid reaction SMARTS"):
#         mol_subsearch(arr, query)
#
#
def test_bad_query_type():
    arr = np.array(['C'])
    query = 0
    with pytest.raises(TypeError, match="Query must be of string type"):
        mol_subsearch(arr, query)


def test_empty_query():
    arr = np.array(['C'])
    query = ""
    with pytest.raises(ValueError, match="Query could not be empty or None"):
        mol_subsearch(arr, query)


""" TEST THIRD ARGUMENT - mode """


""" TEST THIRD ARGUMENT - n_jobs """

def test_add_option_n_jobs():
    arr = np.array(['C'])
    query = "C"
    res = mol_subsearch(arr, query=query, n_jobs=1)
    assert res[0]