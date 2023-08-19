from __future__ import annotations

import numpy as np
import pytest
from pymatgen.core import Structure


@pytest.fixture
def rutile() -> Structure:
    a = 4.603
    c = 2.969
    x_4f = 0.3046
    lattice = np.diag([a, a, c])
    frac_coords = np.array(
        [  # Fractional coordinates
            [0, 0, 0],  # Mn(2a)
            [0.5, 0.5, 0.5],  # Mn(2a)
            [x_4f, x_4f, 0],  # F(4f)
            [-x_4f, -x_4f, 0],  # F(4f)
            [-x_4f + 0.5, x_4f + 0.5, 0.5],  # F(4f)
            [x_4f + 0.5, -x_4f + 0.5, 0.5],  # F(4f)
        ]
    )
    species = ["Ti"] * 2 + ["O"] * 4
    return Structure(lattice, species, frac_coords)


@pytest.fixture
def afm_nio() -> Structure:
    a = 5.687

    lattice = a * np.eye(3)
    species = ["Ni"] * 4 + ["O"] * 4
    frac_coords = np.array(
        [
            [0.0, 0.0, 0.0],  # Ni
            [0.0, 0.5, 0.5],  # Ni
            [0.5, 0.0, 0.5],  # Ni
            [0.5, 0.5, 0.0],  # Ni
            [0.5, 0.5, 0.5],  # O
            [0.5, 0.0, 0.0],  # O
            [0.0, 0.5, 0.0],  # O
            [0.0, 0.0, 0.5],  # O
        ]
    )
    mx = 0.569
    magmom = np.array(
        [
            [mx, mx, mx],
            [-mx, mx, -mx],
            [-mx, -mx, mx],
            [mx, -mx, -mx],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]
    )
    structure = Structure(lattice, species, frac_coords, site_properties={"magmom": magmom})

    return structure
