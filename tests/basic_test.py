import numpy as np
import pandas as pd
import pytest
import re
import io
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


def test_pandas_pd():
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


def test_bad_reaction_smiles(capfd):
    arr = np.array(['C]>>'])
    query = "C>>"

    res = rxn_subsearch(arr, query_smarts=query)
    captured = capfd.readouterr()

    assert not res[0]
    # assert captured.out == "[23 0] Invalid SMILES: C]>>"


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


def test_get_default_num_cores():
    assert get_option("num_cores") == MAX_CPU_COUNT

def test_set_option_num_cores_str():
    set_option("num_cores", "1")
    assert get_option("num_cores") == 1

def test_set_option_num_cores_int():
    set_option("num_cores", MAX_CPU_COUNT // 2)
    assert get_option("num_cores") == MAX_CPU_COUNT // 2

def test_set_float_num_cores():
    with pytest.raises(TypeError, match="float type not allowed, int or string expected"):
        set_option("num_cores", 1.1)


def test_set_negative_num_cores():
    with pytest.raises(ValueError, match="Negative 'num_cores' not allowed"):
        set_option("num_cores", -1)

def test_set_bad_str_num_cores():
    with pytest.raises(ValueError):
        set_option("num_cores", "1.1")
    with pytest.raises(ValueError):
        set_option("num_cores", "1a")


def test_set_zero_num_cores():
    set_option("num_cores", 0)
    assert get_option("num_cores") == MAX_CPU_COUNT

def test_set_big_num_cores():
    set_option("num_cores", MAX_CPU_COUNT * 2)
    assert get_option("num_cores") == MAX_CPU_COUNT




