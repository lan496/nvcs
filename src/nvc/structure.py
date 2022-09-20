from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import product

import numpy as np
from nglview import NGLWidget, show_pymatgen
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import CrystalNN, NearNeighbors
from pymatgen.core import Structure
from pymatgen.core.sites import PeriodicSite
from pymatgen.electronic_structure.core import Magmom
from scipy.spatial import ConvexHull

from nvc.color import ColorScheme


@dataclass
class PeriodicSiteImage:
    site: PeriodicSite
    site_index: int
    jimage: tuple[int, int, int]


def viewer(
    structure: Structure,
    show_unitcell: bool = True,
    show_bonds: bool = True,
    show_outside_bonds: bool = False,
    show_polyhedrons: bool = True,
    local_env_strategy: NearNeighbors | None = None,
    show_magmom: bool = False,
    magmom_scale: float = 2.0,
    arrow_radius: float = 0.1,
    show_axes: bool = True,
    width: int | None = None,
    height: int | None = None,
) -> NGLWidget:
    """
    Args:
        structure: pymatgen's Structure object
        show_unitcell: Iff true, show frame of unit cell
        show_bonds: Iff true, show bonds with `local_env_strategy`
        show_outside_bonds: Iff true, show bonds and polyhedrons connected to outside sites
        show_polyhedrons: Iff true, show coordination polyhedrons
        local_env_strategy: pymatgen's NearNeighbors object to detect connections of sites. If not specified, CrystalNN is used.
        show_magmom: Iff true, show magnetic moments by arrows
        magmon_scale: length scale of arrows to represent magmom
        show_axes: Iff true, show a, b, and c axes
        width: in pixel
        height: in pixel

    Returns:
        view: NGLWidget
    """
    if not isinstance(structure, Structure):
        raise ValueError("Only support pymatgen.core.Structure.")

    # wrap frac_coords in [0, 1)
    wrapped_sites = []
    for site in structure:
        frac_coords = np.remainder(site.frac_coords, 1)
        wrapped_sites.append(
            PeriodicSite(
                species=site.species,
                coords=frac_coords,
                lattice=structure.lattice,
                properties=site.properties,
            )
        )
    wrapped_structure = Structure.from_sites(wrapped_sites)

    # show image atoms near unitcell
    displayed = _get_displayed(wrapped_structure)
    structure_display = Structure.from_sites([si.site for si in displayed])

    view = show_pymatgen(structure_display)
    view.clear()
    view.center()

    cc = ColorScheme(scheme="jmol")

    view = _add_sites(view, cc, displayed, show_magmom, magmom_scale, arrow_radius)

    if local_env_strategy is None:
        local_env_strategy = CrystalNN()
    sg = StructureGraph.with_local_env_strategy(wrapped_structure, local_env_strategy)

    if show_bonds or show_polyhedrons:
        view = _add_connections(
            view, cc, sg, displayed, show_bonds, show_outside_bonds, show_polyhedrons
        )

    if show_unitcell:
        view.add_unitcell()

    if show_axes:
        view = _add_axes(view, structure.lattice.matrix)

    # ref: https://github.com/nglviewer/nglview/issues/900
    view.control.spin([1, 0, 0], -math.pi / 2)
    view.control.spin([0, 0, 1], math.pi * 0.45)
    view.control.spin([0, 1, 0], math.pi * 0.1)
    view.camera = "perspective"

    if (width is not None) and (height is not None):
        view._set_size(w=f"{width}px", h=f"{height}px")

    return view


def _get_displayed(structure: Structure, eps: float = 1e-8) -> list[PeriodicSiteImage]:
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
                    lattice=structure.lattice,
                    properties=site.properties,
                )
                ghosts.append(
                    PeriodicSiteImage(
                        site=new_site,
                        site_index=site_index,
                        jimage=tuple(map(int, np.around(jimage).tolist())),  # type: ignore
                    )
                )

    displayed = []
    for site_index, site in enumerate(structure):
        displayed.append(PeriodicSiteImage(site=site, site_index=site_index, jimage=(0, 0, 0)))
    displayed.extend(ghosts)

    return displayed


def _add_sites(
    view: NGLWidget,
    cc: ColorScheme,
    displayed: list[PeriodicSiteImage],
    show_magmom: bool,
    magmom_scale: float,
    arrow_radius: float = 0.1,
) -> NGLWidget:
    for i, si in enumerate(displayed):
        # ref: https://github.com/nglviewer/nglview/issues/913
        # selection=[i] is equivalent to f"@{i}" in "selection language".
        # See https://nglviewer.org/ngl/api/manual/usage/selection-language.html
        hex_color = cc.get_hex_color(si.site)
        view.add_spacefill(
            radius=0.5,
            selection=[i],
            color=hex_color,
        )

    if show_magmom:
        positions1 = []
        positions2 = []
        for si in displayed:
            magmom = Magmom(si.site.properties.get("magmom", np.zeros(3)))
            vector = magmom.get_moment() * magmom_scale
            if np.allclose(vector, 0):
                continue
            start = si.site.coords - 0.5 * vector
            end = si.site.coords + 0.5 * vector
            positions1.append(start)
            positions2.append(end)

        color = [1, 0, 0]  # red
        colors = [color for _ in range(len(positions1))]
        radii = [arrow_radius for _ in range(len(positions1))]
        view.shape.add_buffer(
            "arrow",
            position1=np.array(positions1).flatten().tolist(),
            position2=np.array(positions2).flatten().tolist(),
            color=np.array(colors).flatten().tolist(),
            radius=radii,
        )

    return view


def _add_connections(
    view: NGLWidget,
    cc: ColorScheme,
    sg: StructureGraph,
    displayed: list[PeriodicSiteImage],
    show_bonds: bool,
    show_outside_bonds: bool,
    show_polyhedrons: bool,
) -> NGLWidget:
    bonds = []
    polyhedrons = []
    for from_si in displayed:
        connected_sites = sg.get_connected_sites(from_si.site_index, from_si.jimage)
        draw_polyhedron = True
        for connected_site in connected_sites:
            to_si = PeriodicSiteImage(
                site=connected_site.site,
                site_index=connected_site.index,
                jimage=connected_site.jimage,
            )

            if (not show_outside_bonds) and (to_si not in displayed):
                draw_polyhedron = False
                continue

            # We use the same strategy with Crystal Toolkit for drawing coordination polyhedrons
            # ref: https://github.com/materialsproject/crystaltoolkit/blob/main/crystal_toolkit/renderables/site.py
            if (from_si.site.specie > to_si.site.specie) or (
                from_si.site.specie == to_si.site.specie
            ):
                draw_polyhedron = False

            bonds.append((from_si, to_si))

        if draw_polyhedron and (len(connected_sites) > 3):
            polyhedrons.append((from_si.site, connected_sites))

    if show_bonds:
        positions1 = []
        positions2 = []
        colors = []
        for from_si, to_si in bonds:
            # Ref: https://github.com/nglviewer/nglview/issues/912
            color = cc.get_color(from_si.site)
            start = from_si.site.coords
            end = (from_si.site.coords + to_si.site.coords) / 2
            positions1.append(start)
            positions2.append(end)
            colors.append(color)

        cylinder_radius = 0.1
        radii = [cylinder_radius for _ in range(len(colors))]
        view.shape.add_buffer(
            "cylinder",
            position1=np.array(positions1).flatten().tolist(),
            position2=np.array(positions2).flatten().tolist(),
            color=np.array(colors).flatten().tolist(),
            radius=radii,
        )

    if show_polyhedrons:
        for center_site, vertices in polyhedrons:
            # ref: https://github.com/nglviewer/ngl/issues/754
            positions = np.array([list(csite.site.coords) for csite in vertices])
            indices = _get_mesh(positions)

            color = cc.get_color(center_site)
            colors = [color for _ in range(len(positions))]

            # Currently, add_buffer(name='mesh') is not supported
            view.shape.add_mesh(
                positions.flatten().tolist(),  # position
                np.array(colors).flatten().tolist(),  # color
                indices.flatten().tolist(),  # index
            )
            # TODO: specify opacity here via params: BufferParameters
            # ref: http://nglviewer.org/ngl/api/class/src/buffer/mesh-buffer.js~MeshBuffer.html

            # ref: https://github.com/nglviewer/nglview/issues/785#issuecomment-486786411
            opacity = 0.5
            view.update_representation(component=len(view._ngl_component_ids), opacity=opacity)

    return view


def _get_mesh(positions: np.ndarray) -> np.ndarray:
    """
    Args:
        positions: (num, 3)

    Returns:
        indices: (nfacet, 3), indices of triangles composing of convex hull
    """
    hull = ConvexHull(positions)
    center = np.mean(positions, axis=0)

    indices = []
    for triangle in hull.simplices:
        # sort indices in counter clockwise
        d0 = positions[triangle[0]] - center
        d1 = positions[triangle[1]] - center
        d2 = positions[triangle[2]] - center
        if np.dot(np.cross(d0, d1), d2) > 0:
            indices.append([triangle[0], triangle[1], triangle[2]])
        else:
            indices.append([triangle[0], triangle[2], triangle[1]])

    return np.array(indices)


def _add_axes(view: NGLWidget, matrix: np.ndarray) -> NGLWidget:
    # Ref: https://github.com/pyiron/pyiron_atomistics/blob/c5df5e87745d7b575463f7b2a0b588e18007dc40/pyiron_atomistics/atomistics/structure/_visualize.py#L388-L403
    axes_start = matrix[0] + np.array([0, -2, 0])
    arrow_radius = 0.1
    text_size = 1
    text_color = [0, 0, 0]  # black
    arrow_colors = [
        [1.0, 0.0, 0.0],  # Red
        [0.0, 1.0, 0.0],  # Green
        [0.0, 0.0, 1.0],  # Blue
    ]
    arrow_names = ["a", "b", "c"]

    for i, (arrow_color, arrow_name) in enumerate(zip(arrow_colors, arrow_names)):
        start = list(axes_start)
        basis_i = matrix[i]
        shift = basis_i / np.linalg.norm(basis_i)
        end = list(axes_start + shift)

        view.shape.add_arrow(start, end, arrow_color, arrow_radius, f"{arrow_name}-axis")
        view.shape.add_text(end, text_color, text_size, arrow_name)

    return view
