"""
Tests for the merger forest built from AHF's tree files.
"""

import numpy as np
import pytest

from pugs import ahf, merger_forest
from pugs.simulation import halo_id, split_halo_id


@pytest.fixture(scope="module")
def forest(simulation):
    return merger_forest.build_forest(simulation)


def test_forest_holds_every_ahf_halo(simulation, forest):
    """
    The forest is built from all AHF halos, not just those above the catalog's
    particle cut: a halo whose two progenitors both fall below the cut would
    otherwise look like it had no merger at all.
    """
    total = sum(ahf.read_halo_table(snapshot.ahf("halos")).num_rows for snapshot in simulation)
    assert len(forest) == total
    assert np.all(np.diff(forest.halo_id) > 0), "halo_id must be sorted for lookup"


def test_halo_ids_decode_to_their_snapshot(forest):
    snapshot_index, finder_id = split_halo_id(forest.halo_id)
    assert np.array_equal(snapshot_index, forest.snapshot_index)
    assert np.array_equal(finder_id, forest.finder_id)


def test_index_of_round_trips(forest):
    sample = forest.halo_id[:: max(len(forest) // 200, 1)]
    assert np.array_equal(forest.halo_id[forest.index_of(sample)], sample)
    assert forest.index_of([-1])[0] == -1


def test_links_point_one_snapshot_at_a_time(forest):
    """Progenitors are in the previous snapshot, descendants in the next."""
    has_progenitor = forest.main_progenitor >= 0
    assert np.all(
        forest.snapshot_index[forest.main_progenitor[has_progenitor]]
        == forest.snapshot_index[has_progenitor] - 1
    )
    has_descendant = forest.descendant >= 0
    assert np.all(
        forest.snapshot_index[forest.descendant[has_descendant]]
        == forest.snapshot_index[has_descendant] + 1
    )


def test_pruned_graph_is_a_forest(simulation):
    """
    Requiring more than half of a progenitor to end up in the descendant makes
    a second descendant arithmetically impossible. The single-sweep merger
    accumulation depends on that, so check the raw links honour it.
    """
    assert merger_forest.MIN_SHARED_FRACTION >= 0.5
    for snapshot in simulation:
        if not snapshot.has_ahf("croco"):
            continue
        links = ahf.read_merger_links(snapshot.ahf("croco"))
        fraction = links.shared / np.maximum(links.n_progenitor, 1)
        kept = links.progenitor[fraction > merger_forest.MIN_SHARED_FRACTION]
        _, counts = np.unique(kept, return_counts=True)
        assert counts.max(initial=0) <= 1, f"{snapshot.extension}: a progenitor has two descendants"


def test_first_snapshot_has_no_progenitors(simulation, forest):
    first = forest.snapshot_index == 0
    assert np.all(forest.main_progenitor[first] == -1)
    assert np.all(forest.n_progenitors[first] == 0)


def test_final_snapshot_has_no_descendants(simulation, forest):
    final = forest.snapshot_index == simulation.final.index
    assert np.all(forest.descendant[final] == -1)


def test_merger_counts_are_consistent(forest):
    assert np.all(forest.n_mm >= 0)
    # the sentinel and the count must agree in both directions
    assert np.all(forest.z_lmm[forest.n_mm == 0] == -1.0)
    assert np.all(forest.z_lmm[forest.n_mm > 0] > 0.0)


def test_merger_counts_never_decrease_along_the_main_branch(forest):
    """A halo's progenitor tree is contained in its descendant's."""
    has_descendant = np.flatnonzero(forest.descendant >= 0)
    assert np.all(forest.n_mm[forest.descendant[has_descendant]] >= forest.n_mm[has_descendant])


def test_main_progenitor_is_the_highest_merit_link(simulation, forest):
    for snapshot in simulation[1:]:
        links = ahf.read_merger_links(snapshot.ahf("croco"))
        fraction = links.shared / np.maximum(links.n_progenitor, 1)
        keep = fraction > merger_forest.MIN_SHARED_FRACTION
        if not keep.any():
            continue
        best: dict[int, tuple[float, int]] = {}
        for descendant, progenitor, merit in zip(
            links.descendant[keep], links.progenitor[keep], links.merit[keep]
        ):
            if merit > best.get(int(descendant), (-1.0, -1))[0]:
                best[int(descendant)] = (float(merit), int(progenitor))
        descendants = np.array(sorted(best))
        rows = forest.index_of(halo_id(snapshot.index, descendants))
        for row, descendant in zip(rows, descendants):
            if row < 0:
                continue  # a croco descendant with no row in AHF_halos
            expected = halo_id(snapshot.index - 1, best[int(descendant)][1])
            assert forest.halo_id[forest.main_progenitor[row]] == expected


def test_assembly_redshifts_are_ordered(simulation, forest):
    final = np.flatnonzero(forest.snapshot_index == simulation.final.index)
    z25, z50, z75 = forest.mass_percentile_redshifts(final, (0.25, 0.50, 0.75))
    assert np.all(z25 >= z50)
    assert np.all(z50 >= z75)
    # a halo with any progenitor must have assembled at some positive redshift
    with_history = forest.main_progenitor[final] >= 0
    assert np.all(z25[with_history] > 0)


def test_assembly_redshifts_bracket_the_main_branch(simulation, forest):
    """z50 must be a redshift the main branch actually passed through."""
    final = np.flatnonzero(forest.snapshot_index == simulation.final.index)
    z50 = forest.mass_percentile_redshifts(final, (0.5,))[0]
    available = set(np.round(forest.redshift, 9)) | {0.0}
    assert set(np.round(z50, 9)) <= available


def test_particle_cut_shrinks_the_forest(simulation, forest):
    """The cut is available for comparison, but is not what the catalog uses."""
    cut = merger_forest.build_forest(simulation, min_particles=500)
    assert len(cut) < len(forest)
    assert cut.n_mm.sum() < forest.n_mm.sum()


def test_tree_files_may_name_halos_the_catalog_does_not_have(simulation):
    """
    AHF's croco files reference at least one descendant that has no row in the
    matching AHF_halos file. Such links have to be dropped rather than indexed,
    or they silently point at whatever halo sits at index -1.
    """
    dangling = 0
    for snapshot in simulation[1:]:
        links = ahf.read_merger_links(snapshot.ahf("croco"))
        present = set(ahf.read_halo_table(snapshot.ahf("halos")).column("ID").to_pylist())
        dangling += len(set(links.descendant.tolist()) - present)
    assert dangling > 0, "expected the test data to exercise this path"

    forest = merger_forest.build_forest(simulation)
    linked = forest.main_progenitor >= 0
    assert np.all(forest.main_progenitor[linked] < len(forest))
    assert np.all(forest.descendant[forest.descendant >= 0] < len(forest))
