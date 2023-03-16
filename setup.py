from skbuild import setup

setup(
    name="chemivec",
    version="0.1.0",
    description="Vectorized cheminformatics library leveraging EPAM Indigo Toolkit",
    license="MIT",
    packages=['chemivec'],
    package_dir={"": "src"},
    # cmake_install_target="_chemivec",
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.7",
        "pandas"
    ]
)