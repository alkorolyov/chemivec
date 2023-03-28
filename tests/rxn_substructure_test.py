import numpy as np
import pandas as pd
import pytest
import os
import sys
import multiprocessing as mp

try:
    # skbuild
    # package available through pip install
    from chemivec import rxn_subsearch, set_option, get_option
except ModuleNotFoundError:
    # clion build
    # run command `pytest ./tests` from root project folder
    sys.path.append(os.getcwd())
    from src.chemivec import rxn_subsearch, set_option, get_option


MAX_CPU_COUNT = mp.cpu_count()

def test_numpy_npstr():
    arr = np.array(['[C:1]=O>>[C:1]O', 'C=O>>CO'])
    query = "[C:1]=O>>[C:1]O"
    res = rxn_subsearch(arr, query_smarts=query)
    assert res[0]
    assert not res[1]
    assert res.dtype == np.bool_


def test_numpy_pystr():
    arr = np.array(['[C:1]=O>>[C:1]O', 'C=O>>CO'], dtype=object)
    query = "[C:1]=O>>[C:1]O"
    res = rxn_subsearch(arr, query_smarts=query)
    assert res[0]
    assert not res[1]


def test_pylist():
    arr = ['[C:1]=O>>[C:1]O', 'C=O>>CO']
    query = "[C:1]=O>>[C:1]O"
    res = rxn_subsearch(arr, query_smarts=query)
    assert res[0]
    assert not res[1]


def test_pandas_df():
    arr = pd.DataFrame(['[C:1]=O>>[C:1]O', 'C=O>>CO'])
    query = "[C:1]=O>>[C:1]O"
    res = rxn_subsearch(arr, query_smarts=query)
    assert res[0]
    assert not res[1]


def test_pandas_series():
    arr = pd.Series(['[C:1]=O>>[C:1]O', 'C=O>>CO'])
    query = "[C:1]=O>>[C:1]O"
    res = rxn_subsearch(arr, query_smarts=query)
    assert res[0]
    assert not res[1]


def test_multi_dim():
    arr = np.array([[1, 2], [2, 3]])
    query = "[C:1]=O>>[C:1]O"
    with pytest.raises(ValueError, match="Multidimensional input arrays not allowed"):
        rxn_subsearch(arr, query_smarts=query)


def test_empty_array():
    arr = np.array([])
    query = "[C:1]=O>>[C:1]O"
    res = rxn_subsearch(arr, query_smarts=query)
    assert res.shape[0] == 0
    assert res.dtype == np.bool_


def test_bad_reaction_smiles():
    arr = np.array(['C]>>'])
    query = "C>>"
    res = rxn_subsearch(arr, query_smarts=query)
    assert not res[0]


def test_bad_query():
    arr = np.array(['C>>'])
    query = "[C>>"
    with pytest.raises(ValueError, match="Invalid reaction SMARTS"):
        rxn_subsearch(arr, query_smarts=query)


def test_aam_mode():
    arr = np.array(['[C:1]=O>>[C:1]O',
                    'C=O>>[C:1]O',
                    '[C:1]=O>>CO',
                    '[C:1]=O>>C[O:1]',
                    'C=O>>CO'
                    ])
    query = "[C:1]=O>>[C:1]O"
    res = rxn_subsearch(arr, query_smarts=query)
    assert res[0]
    assert not res[1]
    assert not res[2]
    assert not res[3]
    assert not res[4]


def test_no_aam_query():
    arr = np.array(['[C:1]=O>>[C:1]O',
                    'C=O>>[C:1]O',
                    '[C:1]=O>>CO',
                    '[C:1]=O>>C[O:1]',
                    'C=O>>CO'
                    ])
    query = "C=O>>[C:1]O"
    res = rxn_subsearch(arr, query_smarts=query)
    assert res[0]
    assert res[1]
    assert res[2]
    assert res[3]
    assert res[4]


def test_get_default_n_jobs():
    assert get_option("n_jobs") == MAX_CPU_COUNT

def test_set_option_n_jobs_str():
    set_option("n_jobs", "1")
    assert get_option("n_jobs") == 1

def test_set_option_n_jobs_int():
    set_option("n_jobs", MAX_CPU_COUNT // 2)
    assert get_option("n_jobs") == MAX_CPU_COUNT // 2

def test_set_float_n_jobs():
    with pytest.raises(TypeError, match="float type not allowed, int or string expected"):
        set_option("n_jobs", 1.1)


def test_set_negative_n_jobs():
    with pytest.raises(ValueError, match="Negative 'n_jobs' not allowed"):
        set_option("n_jobs", -1)

def test_set_bad_str_n_jobs():
    with pytest.raises(ValueError):
        set_option("n_jobs", "1.1")
    with pytest.raises(ValueError):
        set_option("n_jobs", "1a")


def test_set_zero_n_jobs():
    set_option("n_jobs", 0)
    assert get_option("n_jobs") == MAX_CPU_COUNT

def test_set_big_n_jobs():
    set_option("n_jobs", MAX_CPU_COUNT * 2)
    assert get_option("n_jobs") == MAX_CPU_COUNT




