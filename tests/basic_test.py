import numpy as np
import pandas as pd

from chemivec import rxn_match

def try_import_module(name: str):
    print(f"Importing {name} ... ", end="")
    try:
        exec(f"import {name}")
        print(f"OK")
        return True
    except Exception as e:
        print(f"FAILED {e}")
        return False


def test_import_module():
    import chemivec._chemivec

def test_import_func():
    from chemivec import rxn_match

def test_numpy():
    arr = np.array(["[C:1](=O)C>>[C:1](O)C",
                    "C=O>>CO",
                    "[C:2]=O>>[C:2]O",
                    "[C:1](=O)C>>C(O)[C:1]",
                    "[C:2]=O>>[C:2]O",
                    "[C:1](=O)C>>C(O)[C:1]"], dtype=object)
    query = "[C:1]=[O]>>[C:1]-[OX2]"
    res = rxn_match(arr, query_smarts=query, use_aam=True)
    correct_res = np.array([True, False, True, False, True, False])
    for i in range(len(res)):
        assert res[i] == correct_res[i]

def test_pandas_pd():
    arr = pd.DataFrame(["[C:1](=O)C>>[C:1](O)C",
                    "C=O>>CO",
                    "[C:2]=O>>[C:2]O",
                    "[C:1](=O)C>>C(O)[C:1]",
                    "[C:2]=O>>[C:2]O",
                    "[C:1](=O)C>>C(O)[C:1]"], dtype=object)
    query = "[C:1]=[O]>>[C:1]-[OX2]"
    res = rxn_match(arr, query_smarts=query, use_aam=True)
    correct_res = np.array([True, False, True, False, True, False])
    for i in range(len(res)):
        assert res[i] == correct_res[i]

def test_pandas_series():
    arr = pd.Series(["[C:1](=O)C>>[C:1](O)C",
                        "C=O>>CO",
                        "[C:2]=O>>[C:2]O",
                        "[C:1](=O)C>>C(O)[C:1]",
                        "[C:2]=O>>[C:2]O",
                        "[C:1](=O)C>>C(O)[C:1]"], dtype=object)
    query = "[C:1]=[O]>>[C:1]-[OX2]"
    res = rxn_match(arr, query_smarts=query, use_aam=True)
    correct_res = np.array([True, False, True, False, True, False])
    for i in range(len(res)):
        assert res[i] == correct_res[i]


def test_1000_timeit():
    arr = np.load("tests/test_1000.npy", allow_pickle=True)
    query = "[B;X3,4]-[C,c:1].[C,c:2]-[Cl,Br,I,$([O]-S)]>>[C,c:1]-[C,c:2]"
    res = rxn_match(arr, query, "DAYLIGHT-AAM")