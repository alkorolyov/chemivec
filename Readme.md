# Chemivec

Vectorized Cheminformatics Python library, based on EPAM Indigo toolkit C-API
and using Pandas as NumPy for pinput/poutput.

### Supported operations:
rxn_match - reaction substructure match

### Example usage:
```python
import numpy as np
from chemivec import rxn_match

arr = np.array(['[C:1]=O>>[C:1]O', 'C=O>>CO'])
query = "[C:1]=[O]>>[C:1]-[O]"
res = rxn_match(arr, query_smarts=query)
print(res)

# Output: array([ True, False]) 
```

`Atom-to-atom matching (use_aam)` follows the standard DAYLIGHT SMARTS rules
declared here https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html (Section 4.6 Reaction Queries).
If set to `use_aam=False` all atom-mappings are ignored.


### Build

To build extension:

`python setup.py build_ext`

To build distribution wheel:

`python setup.py bdist_wheel`

To build and install as a pip package:

`pip install .`


### General notes on configuring skbuild packages
You may follow the example of https://github.com/scikit-build/scikit-build-sample-projects/tree/master/projects/hello-cpp

Or if you want to change the directory of your package, for example to `src/mylib`
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

Then you need to modify `setup.py` to include package location.

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
You need to modify your target with C extension to match this folder, so
`setuptools` would be aware of it and will include it into build and wheel.

```cmake
python3_add_library(mylib_c_ext MODULE src/mylib_c_ext.c)
install(TARGETS mylib_c_ext LIBRARY DESTINATION src/mylib)
```

### Using cibuildwheels to create distro
cibuildwheels --windows


### Misc
To check dependencies of your `*.pyd` library
dumpbin should be run from developer command prompt of VS 2022

`dumpbin mylib_c_ext.pyd /DEPENDENTS`