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
    try_import_module("_chemivec")

def test_import_package():
    try_import_module("chemivec")


def test_basic():
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
    arr = np.load("tests/test_10000.npy", allow_pickle=True)
    query = "[B;X3,4]-[C,c:1].[C,c:2]-[Cl,Br,I,$([O]-S)]>>[C,c:1]-[C,c:2]"

    start = time()
    print("Timing 10000 samples ... ", end="")
    try:
        res = _rxn_match(arr, query, "DAYLIGHT-AAM")
        print(f"OK {time() - start:.3f}s")
    except Exception as e:
        print(f"FAILED\n{e}")

def test_empty_str():
    arr = np.array(["[C:1](=O)C>>[C:1](O)C",
                    ""], dtype=object)
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
    correct_res = np.array([True, False])
    for i in range(len(res)):
        assert res[i] == correct_res[i]


def test_numpy_npstr():
    arr = np.array(['[C:1]=O>>[C:1]O', 'C=O>>CO'])

    query = "[C:1]=[O]>>[C:1]-[O]"
    print("type(arr[0])", type(arr[0]))
    print("type(arr[0]) == str", type(arr[0]) == str)
    print("type(arr[0]) == np.str_", type(arr[0]) == np.str_)
    print("isinstance(arr[0], str)", isinstance(arr[0], str))
    print("isinstance(arr[0], np.str_)", isinstance(arr[0], np.str_))

    print(f"Running npstr _rxn_match ... ", end="")
    print(np.version.version)
    if isinstance(arr[0], str):
        print("Pystring")
        return

    try:
        res = rxn_match(arr, query_smarts=query, use_aam=True)
        print(f"OK")
    except Exception as e:
        print(f"FAILED\n{e}")
        return

    correct_res = np.array([True, False])
    for i in range(len(res)):
        assert res[i] == correct_res[i]


# test_import()
# test_import_package()

from src.chemivec._chemivec import _rxn_match
from src.chemivec import rxn_match

# test_basic()
# test_timeit()
# test_empty_str()
test_numpy_npstr()




