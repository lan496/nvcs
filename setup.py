from setuptools import find_packages, setup

setup(
    name="nvc",
    version="0.0.1",
    license="MIT",
    author="Kohei Shinohara",
    author_email="kohei19950508@gmail.com",
    description="nglview wrapper for crystal structure",
    packages=find_packages("nvc"),
    python_requires=">=3.7",
    install_requires=[
        "setuptools",
        "wheel",
        "notebook",
        "nglview",
        "ase",
        "pymatgen",
        "scipy",
        "numpy",
    ],
    extra_require={
        "dev": [
            "pre-commit",
            "black",
            "flake8",
            "mypy",
            "isort",
        ],
    },
)
