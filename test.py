# import sys
# print(sys.path)

def try_import_module(name: str):
    print(f"Importing {name} ... ")
    try:
        exec(f"import {name}")
        print(f"    import successful")
        return True
    except Exception as e:
        print(f"    {e}")
        return False


if try_import_module("chemivec"):
    from chemivec import _rxn_match
    import numpy as np
    arr = np.array(["[C:1](=O)C>>[C:1](O)C",
                    "C=O>>CO",
                    "[C:2]=O>>[C:2]O",
                    "[C:1](=O)C>>C(O)[C:1]",
                    "[C:2]=O>>[C:2]O",
                    "[C:1](=O)C>>C(O)[C:1]"], dtype=object)
    query = "[C:1]=[O]>>[C:1]-[OX2]"
    assert isinstance(arr, np.ndarray)
    assert isinstance(arr[0], str)

    try:
        res = _rxn_match(arr, query_smarts=query, aam_mode="DAYLIGHT-AAM")
        print(f"rxn_match successful\n{res}")
    except Exception as e:
        print(f"rxn_match failed\n{e}")


