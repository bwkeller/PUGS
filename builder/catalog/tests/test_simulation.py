"""Tests for snapshot discovery."""

import numpy as np
import pytest

from pugs.simulation import HALO_ID_STRIDE, find_simulation, halo_id, split_halo_id


def test_finds_every_snapshot(simulation):
    assert len(simulation) == 16
    assert simulation.name == "NUGS128"


def test_snapshots_are_in_time_order(simulation):
    redshifts = [snapshot.redshift for snapshot in simulation]
    assert redshifts == sorted(redshifts, reverse=True)
    assert simulation.final.redshift == pytest.approx(0.0, abs=1e-6)
    assert [s.index for s in simulation] == list(range(len(simulation)))


def test_header_redshift_refines_the_filename(simulation):
    """
    The AHF filename rounds the redshift to three decimals; the snapshot header
    is exact. They must still agree to the filename's precision.
    """
    for snapshot in simulation:
        assert snapshot.redshift == pytest.approx(snapshot.catalog_redshift, abs=5e-4)
        assert snapshot.scalefactor == pytest.approx(1.0 / (1.0 + snapshot.redshift))


def test_ahf_files_resolve(simulation):
    for snapshot in simulation:
        assert snapshot.has_ahf("halos")
        assert snapshot.has_ahf("particles")
        assert snapshot.ahf("halos").name.endswith(".AHF_halos")
    # the first snapshot has no previous one, so no merger tree
    assert not simulation[0].has_ahf("croco")
    assert all(snapshot.has_ahf("croco") for snapshot in simulation[1:])


def test_halo_id_round_trips():
    for snapshot_index, finder_id in [(0, 0), (3, 17), (63, 30_000_000)]:
        value = halo_id(snapshot_index, finder_id)
        assert split_halo_id(value) == (snapshot_index, finder_id)


def test_halo_id_is_vectorised_and_ordered():
    finder = np.arange(5, dtype=np.int64)
    assert np.array_equal(halo_id(2, finder), 2 * HALO_ID_STRIDE + finder)
    # ids from a later snapshot always sort after an earlier one
    assert halo_id(3, 0) > halo_id(2, HALO_ID_STRIDE - 1)


def test_missing_directory_is_an_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        find_simulation(tmp_path)
