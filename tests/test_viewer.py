from __future__ import annotations

import pytest
from pymatgen.analysis.local_env import VoronoiNN

from nvc import viewer


@pytest.mark.parametrize("structure_fixture", ["rutile", "afm_nio"])
@pytest.mark.parametrize("local_env_strategy", [None, VoronoiNN()])
def test_viewer(request, structure_fixture, local_env_strategy):
    structure = request.getfixturevalue(structure_fixture)
    viewer(structure, local_env_strategy=local_env_strategy)
