from __future__ import annotations
from itertools import product
from typing import Optional, List, Tuple
from dataclasses import dataclass

from nglview import NGLWidget, show_pymatgen
from pymatgen.core import Structure
from pymatgen.core.sites import PeriodicSite
from pymatgen.analysis.local_env import NearNeighbors, CrystalNN
from pymatgen.analysis.graphs import StructureGraph
import numpy as np


@dataclass
class PeriodicSiteImage:
    site: PeriodicSite
    site_index: int
    jimage: Tuple[int, int, int]

    def __lt__(self, rhs: PeriodicSiteImage):
        if self.site_index != rhs.site_index:
            return (self.site_index < rhs.site_index)
        else:
            return (self.jimage < rhs.jimage)

    def __gt__(self, rhs: PeriodicSiteImage):
        return rhs.__lt__(self)


def viewer(
    structure: Structure,
    show_unitcell: bool = True,
    show_bonds: bool = True,
    show_outside_bonds: bool = False,
    show_axes: bool = True,
    local_env_strategy: Optional[NearNeighbors] = None,
) -> NGLWidget:
    """
    Args:
        structure:
        show_unitcell: show frame of unit cell iff true
        show_bonds: show chemical bonds with `local_env_strategy` iff true
        show_outside_bonds:
        show_axes: show a, b, and c axes iff true
    """
    # wrap frac_coords in [0, 1)
    wrapped_sites = []
    for site in structure:
        frac_coords = np.remainder(site.frac_coords, 1)
        wrapped_sites.append(
            PeriodicSite(
                species=site.species,
                coords=frac_coords,
                lattice=structure.lattice
            )
        )
    wrapped_structure = Structure.from_sites(wrapped_sites)

    # show image atoms near unitcell
    displayed = _get_displayed(wrapped_structure)
    structure_display = Structure.from_sites([si.site for si in displayed])

    view = show_pymatgen(structure_display)
    view.center()

    view.add_spacefill(radius=0.5, color_scheme='element')
    view.remove_ball_and_stick()

    if local_env_strategy is None:
        local_env_strategy = CrystalNN()
    sg = StructureGraph.with_local_env_strategy(wrapped_structure, local_env_strategy)

    if show_bonds:
        view = _add_bonds(view, sg, displayed, show_outside_bonds)

    if show_unitcell:
        view.add_unitcell()

    if show_axes:
        view = _add_axes(view, structure.lattice.matrix)

    view.camera = "perspective"

    return view


def _get_displayed(structure: Structure, eps: float = 1e-8) -> List[PeriodicSiteImage]:
    ghosts = []
    for site_index, site in enumerate(structure):
        for jimage in product([0, 1 - eps], repeat=3):
            # skip original site
            if np.allclose(jimage, 0):
                continue

            new_frac_coords = site.frac_coords + np.array(jimage)
            if np.all(new_frac_coords < 1 + eps):
                new_site = PeriodicSite(
                    species=site.species,
                    coords=new_frac_coords,
                    lattice=structure.lattice
                )
                ghosts.append(
                    PeriodicSiteImage(
                        site=new_site,
                        site_index=site_index,
                        jimage=tuple(map(int, np.around(jimage).tolist()))
                    )
                )

    displayed = []
    for site_index, site in enumerate(structure):
        displayed.append(PeriodicSiteImage(site=site, site_index=site_index, jimage=(0, 0, 0)))
    displayed.extend(ghosts)

    return displayed


def _add_bonds(
    view: NGLWidget,
    sg: StructureGraph,
    displayed: List[PeriodicSiteImage],
    show_outside_bonds: bool,
) -> NGLWidget:
    bonds = []
    for from_si in displayed:
        for connected_site in sg.get_connected_sites(from_si.site_index, from_si.jimage):
            to_si = PeriodicSiteImage(
                site=connected_site.site,
                site_index=connected_site.index,
                jimage=connected_site.jimage,
            )

            if (not show_outside_bonds) and (to_si not in displayed):
                continue

            # TODO: avoid double counting
            bonds.append((from_si, to_si))

    for from_si, to_si in bonds:
        # Ref: https://github.com/nglviewer/nglview/issues/912
        color = [0, 0, 0]
        cylinder_radius = 0.1
        view.shape.add_cylinder(
            from_si.site.coords.tolist(),  # position1
            to_si.site.coords.tolist(),  # position2
            color,  # color
            cylinder_radius,  # radius
        )

    return view


def _add_axes(view: NGLWidget, matrix: np.ndarray) -> NGLWidget:
    # Ref: https://github.com/pyiron/pyiron_atomistics/blob/c5df5e87745d7b575463f7b2a0b588e18007dc40/pyiron_atomistics/atomistics/structure/_visualize.py#L388-L403
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
        basis_i = matrix[i]
        shift = basis_i / np.linalg.norm(basis_i)
        end = list(axes_start + shift)

        view.shape.add_arrow(start, end, arrow_color, arrow_radius, f"{arrow_name}-axis")
        view.shape.add_text(end, text_color, text_size, arrow_name)

    return view
