import numpy as np
from time import time


def try_import_module(name: str):
    print(f"Importing {name} ... ", end="")
    try:
        exec(f"import {name}")
        print(f"OK")
        return True
    except Exception as e:
        print(f"FAILED {e}")
        return False


def test_import():
    try_import_module("chemivec")


def test_basic():
    from chemivec import _rxn_match
    arr = np.array(["[C:1](=O)C>>[C:1](O)C",
                    "C=O>>CO",
                    "[C:2]=O>>[C:2]O",
                    "[C:1](=O)C>>C(O)[C:1]",
                    "[C:2]=O>>[C:2]O",
                    "[C:1](=O)C>>C(O)[C:1]"], dtype=object)
    query = "[C:1]=[O]>>[C:1]-[OX2]"
    assert isinstance(arr, np.ndarray)
    assert isinstance(arr[0], str)
    print(f"Running _rxn_match ... ", end="")
    try:
        res = _rxn_match(arr, query_smarts=query, aam_mode="DAYLIGHT-AAM")
        print(f"OK")
    except Exception as e:
        print(f"FAILED\n{e}")
        return
    correct_res = np.array([True, False, True, False, True, False])
    for i in range(len(res)):
        assert res[i] == correct_res[i]


def test_timeit():
    from chemivec import _rxn_match
    arr = np.load("test_10000.npy", allow_pickle=True)
    query = "[B;X3,4]-[C,c:1].[C,c:2]-[Cl,Br,I,$([O]-S)]>>[C,c:1]-[C,c:2]"

    start = time()
    print("Timing 10000 samples ... ", end="")
    try:
        res = _rxn_match(arr, query, "DAYLIGHT-AAM")
        print(f"OK {time() - start:.3f}s")
    except Exception as e:
        print(f"FAILED\n{e}")


test_import()
test_basic()
test_timeit()
