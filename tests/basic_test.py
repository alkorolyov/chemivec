import numpy as np
import pandas as pd
import pytest
import os
import sys


try:
    # skbuild
    # package available through pip install
    from chemivec import rxn_match
except ModuleNotFoundError:
    # clion build
    # run command `pytest ./tests` from root project folder
    sys.path.append(os.getcwd())
    from src.chemivec import rxn_match

def test_numpy_npstr():
    arr = np.array(['[C:1]=O>>[C:1]O', 'C=O>>CO'])
    query = "[C:1]=O>>[C:1]O"
    res = rxn_match(arr.astype(object), query_smarts=query, use_aam=True)
    assert res[0]
    assert not res[1]

def test_numpy_pystr():
    arr = np.array(['[C:1]=O>>[C:1]O', 'C=O>>CO'], dtype=object)
    query = "[C:1]=O>>[C:1]O"
    res = rxn_match(arr, query_smarts=query, use_aam=True)
    assert res[0]
    assert not res[1]

def test_pylist():
    arr = ['[C:1]=O>>[C:1]O', 'C=O>>CO']
    query = "[C:1]=O>>[C:1]O"
    res = rxn_match(arr, query_smarts=query, use_aam=True)
    assert res[0]
    assert not res[1]

def test_pandas_pd():
    arr = pd.DataFrame(['[C:1]=O>>[C:1]O', 'C=O>>CO'])
    query = "[C:1]=O>>[C:1]O"
    res = rxn_match(arr, query_smarts=query, use_aam=True)
    assert res[0]
    assert not res[1]

def test_pandas_series():
    arr = pd.Series(['[C:1]=O>>[C:1]O', 'C=O>>CO'])
    query = "[C:1]=O>>[C:1]O"
    res = rxn_match(arr, query_smarts=query, use_aam=True)
    assert res[0]
    assert not res[1]

def test_incorrect_reaction_smiles():
    arr = np.array(['C]>>'])
    query = "C>>"
    res = rxn_match(arr, query_smarts=query, use_aam=True)
    assert not res[0]

def test_incorrect_query():
    arr = np.array(['C>>'])
    query = "[C>>"
    with pytest.raises(ValueError, match="Invalid SMARTS"):
        rxn_match(arr, query_smarts=query, use_aam=True)