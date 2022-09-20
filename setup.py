from setuptools import find_packages, setup

setup(
    name="nvc",
    version="0.0.2",
    license="MIT",
    author="Kohei Shinohara",
    author_email="kshinohara0508@gmail.com",
    description="nglview wrapper for crystal structure",
    package_dir={"": "src"},
    packages=find_packages(where="src", include=["nvc"]),
    package_data={},
    include_package_data=True,
    zip_safe=False,
    python_requires=">=3.8",
    install_requires=[
        "setuptools",
        "wheel",
        "ipykernel",
        "notebook>=4.2",
        # nglview does not work with ipywidgets>=8, https://github.com/nglviewer/nglview/issues/1032
        "ipywidgets>=7.0, <8",
        "nglview>=3.0.3",
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
            "pyupgrade",
        ],
    },
)
