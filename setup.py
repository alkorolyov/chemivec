from skbuild import setup

setup(
    name="chemivec",
    version="1.0.0",
    description="Vectorized cheminformatics library",
    license="MIT",
    python_requires=">=3.7",
    # packages=['chemivec'],
    cmake_install_target="chemivec",
    # cmake_with_sdist=True
)