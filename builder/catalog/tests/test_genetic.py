"""Tests for the zoom-in IC helpers."""

import numpy as np
import pytest

from pugs import ahf, genetic
from pugs.simulation import halo_id


@pytest.fixture(scope="module")
def first_halo(simulation):
    """The catalog halo_id of AHF halo 0 in the final snapshot."""
    return halo_id(simulation.final.index, 0)


def test_resolve_halo_finds_the_snapshot(simulation, first_halo):
    snapshot, finder_id = genetic.resolve_halo(simulation, first_halo)
    assert snapshot.extension == simulation.final.extension
    assert finder_id == 0


def test_resolve_halo_rejects_an_unknown_snapshot(simulation):
    with pytest.raises(ValueError, match="snapshot"):
        genetic.resolve_halo(simulation, halo_id(len(simulation) + 5, 0))


def test_particle_ids_default_to_the_finder_membership(folder, simulation, first_halo):
    ids = genetic.particle_ids(folder, first_halo)
    membership = ahf.read_particle_membership(simulation.final.ahf("particles"))
    assert np.array_equal(np.sort(ids), np.sort(membership.particles(0)))


def test_radius_factor_expands_the_region(folder, first_halo):
    """
    A zoom region is usually larger than the halo the finder identified, so the
    ids are selected from the snapshot rather than read from AHF.
    """
    members = genetic.particle_ids(folder, first_halo)
    wide = genetic.particle_ids(folder, first_halo, radius_factor=2)
    assert len(wide) > len(members)


def test_regions_are_nested_by_radius(folder, first_halo):
    """
    Ids come back sorted by radius, so a smaller region is a prefix of a larger
    one. That is what lets a caller trim a written id file instead of
    re-selecting the region.
    """
    small = genetic.particle_ids(folder, first_halo, radius_factor=2)
    large = genetic.particle_ids(folder, first_halo, radius_factor=8)
    assert len(large) > len(small)
    assert np.array_equal(large[: len(small)], small)


def test_write_particle_ids_round_trips(folder, first_halo, tmp_path):
    target = tmp_path / "id_file.txt"
    written = genetic.write_particle_ids(folder, first_halo, str(target), radius_factor=2)
    recovered = np.loadtxt(target, dtype=np.int64)
    assert written == len(recovered)
    assert np.array_equal(recovered, genetic.particle_ids(folder, first_halo, radius_factor=2))


def test_build_param_file_fills_the_template(tmp_path):
    template = tmp_path / "template.txt"
    template.write_text(
        "outname {outname}\nbase {base_grid}\nzoom {zoom_grid}\n" "pad {autopad}\nsub {subsample}\n"
    )
    out = tmp_path / "genetIC_zoom.txt"
    genetic.build_param_file(
        filename=str(out),
        outname="halo1",
        base_grid=2048,
        zoom_grid=[10, 2048],
        autopad=1,
        subsample=8,
        template=str(template),
    )
    text = out.read_text()
    assert "outname halo1" in text
    assert "base 2048" in text
    assert "{" not in text, "every placeholder must be substituted"
