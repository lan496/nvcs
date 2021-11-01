from itertools import product

from nglview import NGLWidget, show_pymatgen
from pymatgen.core import Structure
from pymatgen.core.sites import PeriodicSite
import numpy as np


def viewer(
    structure: Structure,
    show_unitcell: bool =True,
    show_axes: bool =True,
) -> NGLWidget:
    """
    Args:
        structure:
        show_unitcell: show frame of unit cell iff true
        show_axes: show a, b, and c axes iff true
    """
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
    view.center()
    if show_unitcell:
        view.add_unitcell()

    view.add_spacefill(radius=0.5, color_scheme='element')
    view.remove_ball_and_stick()

    # Ref: https://github.com/pyiron/pyiron_atomistics/blob/c5df5e87745d7b575463f7b2a0b588e18007dc40/pyiron_atomistics/atomistics/structure/_visualize.py#L388-L403
    if show_axes:
        axes_start = -np.ones(3)
        arrow_radius = 0.1
        text_size = 1
        text_color = [0, 0, 0]  # black
        arrow_colors = [
            [1.0, 0.0, 0.0],  # Red
            [0.0, 1.0, 0.0],  # Green
            [0.0, 0.0, 1.0],  # Blue
        ]
        arrow_names = ['a', 'b', 'c']

        for i, (arrow_color, arrow_name) in enumerate(zip(arrow_colors, arrow_names)):
            start = list(axes_start)
            basis_i = structure.lattice.matrix[i]
            shift = basis_i / np.linalg.norm(basis_i)
            end = list(axes_start + shift)

            view.shape.add_arrow(start, end, arrow_color, arrow_radius, f"{arrow_name}-axis")
            view.shape.add_text(end, text_color, text_size, arrow_name)

    view.camera = "perspective"

    return view
