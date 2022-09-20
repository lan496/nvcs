from __future__ import annotations

import json
import os

from pymatgen.core.periodic_table import Element, Species
from pymatgen.core.sites import PeriodicSite


class ColorScheme:
    def __init__(self, scheme="jmol") -> None:
        # TODO: Add more color_scheme
        with open(os.path.join(os.path.dirname(__file__), "color_scheme.json")) as f:
            color_scheme = json.load(f)
        self.color_scheme = color_scheme[scheme]

    def get_color(self, site: PeriodicSite) -> list[float]:
        key = self._get_key(site.specie)
        color = [c / 256 for c in self.color_scheme[key]]
        return color

    def get_hex_color(self, site: PeriodicSite) -> str:
        key = self._get_key(site.specie)
        color = self.color_scheme[key]
        hex = f"#{color[0]:02x}{color[1]:02x}{color[2]:02x}"
        return hex

    def _get_key(self, specie: Element | Species) -> str:
        if isinstance(specie, Element):
            return str(specie)
        else:
            # with oxidation state
            return str(specie.element)
