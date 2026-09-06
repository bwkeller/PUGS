"""
Shared fixtures.

The tests run against the NUGS128 test simulation.  Point ``PUGS_SIMULATION``
at it (``source test_vars``) or place it at ``~/data/NUGS128``.
"""

import os
from pathlib import Path

import pytest

from pugs import export
from pugs.io import read_catalog
from pugs.simulation import find_simulation

SIMULATION_NAME = "NUGS128"


def simulation_folder() -> Path:
    folder = Path(os.environ.get("PUGS_SIMULATION", Path.home() / "data" / SIMULATION_NAME))
    if not folder.is_dir():
        pytest.skip(f"test simulation not found at {folder}; set PUGS_SIMULATION")
    return folder


@pytest.fixture(scope="session")
def folder() -> Path:
    return simulation_folder()


@pytest.fixture(scope="session")
def simulation(folder):
    return find_simulation(folder, SIMULATION_NAME)


@pytest.fixture(scope="session")
def final_snapshot(simulation):
    return simulation.final


@pytest.fixture(scope="session")
def snapshot_data(final_snapshot):
    """The final snapshot loaded once, in physical units, with its membership."""
    from pugs import ahf, halo_properties

    sim = final_snapshot.load()
    membership = ahf.read_particle_membership(final_snapshot.ahf("particles"))
    return sim, membership, halo_properties.critical_density(sim)


@pytest.fixture(scope="session")
def catalog_dir(folder, tmp_path_factory):
    """A full catalog, built once for the whole session."""
    out = tmp_path_factory.mktemp("catalog")
    export.export_simulation(folder, out, name=SIMULATION_NAME)
    return out


@pytest.fixture(scope="session")
def catalog(catalog_dir):
    return read_catalog(catalog_dir)
