# Chemivec

Vectorized Cheminformatics Python library, based on EPAM Indigo toolkit C-API
and using Pandas as NumPy for input/output.

### Supported operations:
rxn_match - reaction substructure match

### Example usage:
```python
import numpy as np
from chemivec import rxn_match

arr = np.array(['[C:1]=O>>[C:1]O', 'C=O>>CO'])
query = "[C:1]=[O]>>[C:1]-[O]"
res = rxn_match(arr, query_smarts=query, use_aam=True)
print(res)

# Output: array([ True, False]) 
```


### Build

To build extension:

`python setup.py build_ext`

To build distribution wheel:

`python setup.py bdist_wheel`

To build and install as a pip package:

`pip install .`


### general notes on configuring packages
You may follow the example of https://github.com/scikit-build/scikit-build-sample-projects/tree/master/projects/hello-cpp

If you want to change the library of your package for example to `src/mylib`
```
project/
src/
├─ mylib/
│  ├─ __init__.py
├─ mylib_c_ext.c
tests/
CMakeLists.txt
pyproject.toml
setup.py
```

Then you need to modify `setup.py` to include the path to packages location.

```python
from skbuild import setup
setup(
    name="mylib",
    version="1.0.0",
    packages=['mylib'],
    package_dir={"": "src"}
)
```

Right now your package would be installed into ` ... \cmake-install\src\build`
So you need to modify your target with C extension to match this folder, so
setuptools would be aware of it.

```cmake
add_library(mylib_c_ext MODULE mylib_c_ext.c)
install(TARGETS mylib_c_ext LIBRARY DESTINATION src/mylib)
```



### to check dependencies of your *.pyd library
### dumpbin should be run from developer command prompt of VS 2022
`dumpbin mylib_c_ext.pyd /DEPENDENTS`