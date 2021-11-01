from itertools import product

from nglview import NGLWidget, show_pymatgen
from pymatgen.core import Structure
from pymatgen.core.sites import PeriodicSite
import numpy as np


def viewer(structure: Structure, spacefill=True, show_axes=True) -> NGLWidget:
    eps = 1e-8
    sites = []
    for site in structure:
        species = site.species
        frac_coords = np.remainder(site.frac_coords, 1)
        for jimage in product([0, 1 - eps], repeat=3):
            new_frac_coords = frac_coords + np.array(jimage)
            if np.all(new_frac_coords < 1 + eps):
                new_site = PeriodicSite(species=species, coords=new_frac_coords, lattice=structure.lattice)
                sites.append(new_site)
    structure_display = Structure.from_sites(sites)

    view = show_pymatgen(structure_display)
    view.add_unitcell()

    if spacefill:
        view.add_spacefill(radius_type='vdw', radius=0.5, color_scheme='element')
        view.remove_ball_and_stick()
    else:
        view.add_ball_and_stick()

    if show_axes:
        view.shape.add_arrow([-4, -4, -4], [0, -4, -4], [1, 0, 0], 0.5, "x-axis")
        view.shape.add_arrow([-4, -4, -4], [-4, 0, -4], [0, 1, 0], 0.5, "y-axis")
        view.shape.add_arrow([-4, -4, -4], [-4, -4, 0], [0, 0, 1], 0.5, "z-axis")

    view.camera = "perspective"
    return view
