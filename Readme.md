# Chemivec

Vectorized Cheminformatics Python library, based on EPAM Indigo toolkit C-API
and using NumPy for input/output.

### Supported operations:
```
rxn_subsearch(input, query) - reaction substructure match
    input : reaction SMILES array (numpy, pandas and python list supported)
    query : reaction query SMARTS, ex "C=C>>C-C"
```

### Example usage:

```python
import numpy as np
import chemivec

arr = np.array(['[C:1]=O>>[C:1]O', 'C=O>>CO'])
query = "[C:1]=O>>[C:1]O"
res = chemivec.rxn_subsearch(arr, query=query)
print(res)

# Output: array([ True, False]) 
```

### Multithreading

Multithreading realized by OpenMP library. By default, tries to use maximum available number of cores.
Number of cores can be specified as a global option or passed as a parameter.

```python
import chemivec

chemivec.rxn_subsearch(arr, query=query)   # default max available cores
chemivec.set_option("n_jobs", 12)                 # change defaults
chemivec.rxn_subsearch(arr, query=query, n_jobs=8)
```

### Atom-to-atom matching (AAM) 
If atom mapping is present in the query, ex `[C:1]>>[C:1]` chemivec follows the standard DAYLIGHT SMARTS rules
declared here https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html (Section 4.6 Reaction Queries)


### Install

Download from pip

`pip install chemivec`

### Build from sources

`python3 -m twine check wheelhouse/*`


### Misc
To check dependencies of your `*.pyd` library
dumpbin should be run from developer command prompt of VS 2022

`dumpbin mylib_c_ext.pyd /DEPENDENTS`
