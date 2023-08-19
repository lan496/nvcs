# nvcs
nglview wrapper for crystal structure

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
