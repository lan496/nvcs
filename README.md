# nvcs
[![testing](https://github.com/lan496/nvcs/actions/workflows/testing.yml/badge.svg)](https://github.com/lan496/nvcs/actions/workflows/testing.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/nvcs.svg)](https://badge.fury.io/py/nvcs)
![PyPI - Downloads](https://img.shields.io/pypi/dm/nvcs)

nglview wrapper for crystal structure

- GitHub: https://github.com/lan496/nvcs/
- PyPI: https://pypi.org/project/nvcs/

## Gallery

<p>
    <img src="examples/rutile.png" alt="rutile.png" width=200>
    <img src="examples/NiS.png" alt="NiS.png" width=200>
    <img src="examples/hi_quartz.png" alt="hi_quartz.png" width=200>
</p>

## Installation

```
pip install nvcs
```

## Example

See [examples/example.ipynb](examples/example.ipynb) for usage

```python
from nvcs import viewer
from pymatgen.core import Structure

structure = Structure.from_file('some_structure.cif')

view = viewer(structure)
view
```
