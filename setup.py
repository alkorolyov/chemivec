from skbuild import setup

setup(
    name="chemivec",
    version="1.0.0",
    description="Vectorized cheminformatics library leveraging EPAM Indigo Toolkit",
    license="MIT",
    python_requires=">=3.5",
    packages=['chemivec'],
    package_dir={"": "src"},

)