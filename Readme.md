# Chemivec

Vectorized Cheminformatics Python library, based on EPAM Indigo toolkit C-API
and using NumPy for input/output.

### Supported operations:
`rxn_subsearch` - reaction substructure match

    `input`         : reaction SMILES array (numpy, pandas and python list supported)
    `query_smarts`  : reaction query SMARTS

### Example usage:

```python
import numpy as np
import chemivec

arr = np.array(['[C:1]=O>>[C:1]O', 'C=O>>CO'])
query = "[C:1]=O>>[C:1]O"
res = chemivec.rxn_subsearch(arr, query_smarts=query)
print(res)

# Output: array([ True, False]) 
```

### Multithreading

Multithreading realized by OpenMP library. By default, tries to use maximum available number of cores.
Number of cores can be specified as a global option or passed as a parameter.

```python
import chemivec

chemivec.rxn_subsearch(arr, query_smarts=query)   # default max available cores
chemivec.set_option("n_jobs", 12)                 # change defaults
chemivec.rxn_subsearch(arr, query_smarts=query, n_jobs=8)
```

`Atom-to-atom matching` follows the standard DAYLIGHT SMARTS rules
declared here https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html (Section 4.6 Reaction Queries).


### Install

Download from pip

`pip install chemivec`

### Build from sources

#### Ubuntu
sudo apt install build-essential ninja-build mc wget git libcairo2-dev zlib1g-dev -y
git clone https://github.com/alkorolyov/chemivec

wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh;chmod +x Mambaforge-Linux-x86_64.sh;bash Mambaforge-Linux-x86_64.sh;export MAMBA_NO_BANNER=1
# if conda still not seen then ~/.bashrc is not sourced when you log in using SSH.
# You need to source it in your ~/.bash_profile like this:
echo "if [ -f ~/.bashrc ]; then
. ~/.bashrc
fi" >> ~/.bash_profile
# restart shell
conda config --set auto_activate_base false
mamba create -n cp310 python=3.10
mamba activate cp310
mamba install pip pytest -y
pip install .

# (optional) to build in cibuildwheel
pip install cibuildwheel
sudo apt-get install docker.io -y; sudo groupadd docker; sudo usermod -aG docker $USER
sudo reboot now
cd chemivec
cibuildwheel --platform linux

#### Windows

# mingw64 on windows
# download stable mingw64 release, extract and add to %Path%
https://github.com/brechtsanders/winlibs_mingw/releases/download/11.2.0-10.0.0-msvcrt-r1/winlibs-x86_64-posix-seh-gcc-11.2.0-mingw-w64msvcrt-10.0.0-r1.zip
# download ninja and also add to %Path%
https://github.com/ninja-build/ninja/releases/download/v1.11.1/ninja-win.zip
cmake -B build -G "Ninja" -D CMAKE_C_COMPILER=gcc.exe -D CMAKE_CXX_COMPILER=g++.exe .
cmake --build build --target _chemivec

### Misc
To check dependencies of your `*.pyd` library
dumpbin should be run from developer command prompt of VS 2022

`dumpbin mylib_c_ext.pyd /DEPENDENTS`
