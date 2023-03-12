# python with numpy and scikit-build should be installed
pip install numpy
pip install scikit-build
pip install setuptools


# to run PowerShell in CLion with conda environment
Set-ExecutionPolicy -ExecutionPolicy Unrestricted -Scope CurrentUser

# to check dependencies of your *.pyd library
# dumpbin should be run from developer command prompt of VS 2022
dumpbin <library name>.pyd /DEPENDENTS
