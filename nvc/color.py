from typing import Tuple
import os
import json


class ColorScheme:

    def __init__(self, scheme='jmol') -> None:
        # TODO: Add more color_scheme
        with open(os.path.join(os.path.dirname(__file__), 'color_scheme.json')) as f:
            color_scheme = json.load(f)
        self.color_scheme = color_scheme[scheme]

    def get_color(self, key) -> Tuple[float, float, float]:
        color = [c / 256 for c in self.color_scheme[key]]
        return color

    def get_hex_color(self, key) -> str:
        color = self.color_scheme[key]
        hex = f'#{color[0]:02x}{color[1]:02x}{color[2]:02x}'
        return hex
