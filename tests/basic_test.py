import numpy as np
import pandas as pd
from chemivec import rxn_match

def test_import_module():
    import chemivec._chemivec


def test_import_func():
    from chemivec import rxn_match


def test_numpy_npstr():
    arr = np.array(['[C:1]=O>>[C:1]O', 'C=O>>CO'])
    query = "[C:1]=[O]>>[C:1]-[O]"
    res = rxn_match(arr.astype(object), query_smarts=query, use_aam=True)

    correct_res = np.array([True, False])
    for i in range(len(res)):
        assert res[i] == correct_res[i]


def test_numpy_pystr():
    arr = np.array(['[C:1]=O>>[C:1]O', 'C=O>>CO'], dtype=object)
    query = "[C:1]=[O]>>[C:1]-[O]"
    res = rxn_match(arr, query_smarts=query, use_aam=True)

    correct_res = np.array([True, False])
    for i in range(len(res)):
        assert res[i] == correct_res[i]


def test_pylist():
    arr = ['[C:1]=O>>[C:1]O', 'C=O>>CO']
    query = "[C:1]=[O]>>[C:1]-[O]"
    res = rxn_match(arr, query_smarts=query, use_aam=True)

    correct_res = np.array([True, False])
    for i in range(len(res)):
        assert res[i] == correct_res[i]


def test_pandas_pd():
    arr = pd.DataFrame(['[C:1]=O>>[C:1]O', 'C=O>>CO'])
    query = "[C:1]=[O]>>[C:1]-[O]"
    res = rxn_match(arr, query_smarts=query, use_aam=True)

    correct_res = np.array([True, False])
    for i in range(len(res)):
        assert res[i] == correct_res[i]

def test_pandas_series():
    arr = pd.Series(['[C:1]=O>>[C:1]O', 'C=O>>CO'])
    query = "[C:1]=[O]>>[C:1]-[O]"
    res = rxn_match(arr, query_smarts=query, use_aam=True)

    correct_res = np.array([True, False])
    for i in range(len(res)):
        assert res[i] == correct_res[i]


def test_incorrect_smi():
    arr = ['C:1=O>>[C:1]O']
    query = "[C:1]=[O]>>[C:1]-[O]"
    res = rxn_match(arr, query_smarts=query, use_aam=True)
    assert not res[0]



# def test_1000_timeit():
#     arr = np.load("./tests/test_1000.npy", allow_pickle=True)
#     query = "[B;X3,4]-[C,c:1].[C,c:2]-[Cl,Br,I,$([O]-S)]>>[C,c:1]-[C,c:2]"
#     res = rxn_match(arr, query, "DAYLIGHT-AAM")