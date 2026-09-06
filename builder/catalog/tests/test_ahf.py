"""Tests for the AHF file readers."""

import numpy as np
import pytest

from pugs import ahf


def test_halo_table_has_ahf_column_names(final_snapshot):
    table = ahf.read_halo_table(final_snapshot.ahf("halos"))
    assert table.column_names[:5] == ["ID", "hostHalo", "numSubStruct", "Mhalo", "npart"]
    # the "(N)" suffixes in the header must be stripped, not carried into names
    assert not any("(" in name for name in table.column_names)
    assert table.num_rows > 0


def test_halo_table_types_are_inferred(final_snapshot):
    """Integer columns must stay integers; the leading whitespace in AHF's
    fixed-width output otherwise turns them into strings or floats."""
    import pyarrow as pa

    table = ahf.read_halo_table(final_snapshot.ahf("halos"))
    for name in ("ID", "hostHalo", "numSubStruct", "npart"):
        assert pa.types.is_integer(table.schema.field(name).type), name


def test_halos_have_substructure_information(final_snapshot):
    """
    hostHalo is genuinely populated in the AHF output. The TANGOS pipeline
    imported nothing for it, so this is data the catalog previously lacked.
    """
    table = ahf.read_halo_table(final_snapshot.ahf("halos"))
    host = np.asarray(table.column("hostHalo").to_numpy())
    assert (host >= 0).sum() > 0, "expected some subhalos"
    assert (host == -1).sum() > 0, "expected some field halos"


def test_particle_membership_matches_the_halo_table(final_snapshot):
    table = ahf.read_halo_table(final_snapshot.ahf("halos"))
    membership = ahf.read_particle_membership(final_snapshot.ahf("particles"))

    assert len(membership) == table.num_rows
    npart = dict(zip(table.column("ID").to_pylist(), table.column("npart").to_pylist()))
    for halo in list(npart)[:50]:
        assert len(membership.particles(halo)) == npart[halo]


def test_particle_membership_offsets_are_consistent(final_snapshot):
    membership = ahf.read_particle_membership(final_snapshot.ahf("particles"))
    assert membership.offsets[0] == 0
    assert membership.offsets[-1] == len(membership.indices)
    assert np.all(np.diff(membership.offsets) >= 0)
    assert np.all(np.diff(membership.halo_id) > 0), "halo_id must be sorted for lookup"


def test_particle_membership_rejects_unknown_halo(final_snapshot):
    membership = ahf.read_particle_membership(final_snapshot.ahf("particles"))
    with pytest.raises(KeyError):
        membership.particles(10**9)


def test_particle_ids_are_snapshot_indices(final_snapshot):
    """
    The tipsy files carry no iord array, so AHF's particle values index the
    snapshot directly. If that stopped being true every measured property
    would silently be computed on the wrong particles.
    """
    membership = ahf.read_particle_membership(final_snapshot.ahf("particles"))
    sim = final_snapshot.load()
    assert membership.indices.min() >= 0
    assert membership.indices.max() < len(sim)


def test_substructure_agrees_with_the_halo_table(final_snapshot):
    table = ahf.read_halo_table(final_snapshot.ahf("halos"))
    substructure = ahf.read_substructure(final_snapshot.ahf("substructure"))
    declared = dict(zip(table.column("ID").to_pylist(), table.column("numSubStruct").to_pylist()))
    for halo, children in substructure.items():
        assert len(children) == declared[halo], halo


def test_merger_links_carry_real_shared_particle_counts(final_snapshot):
    """
    AHF_croco records SharedPart and merit for every link. AHF_mtree holds the
    same edges with those numbers stripped out, leaving only an ordering; the
    merit is what makes a physically meaningful link cut possible.
    """
    links = ahf.read_merger_links(final_snapshot.ahf("croco"))
    assert len(links) > 0
    assert np.all(links.shared > 0)
    assert np.all(links.merit > 0)
    # merit is shared**2 / (n_descendant * n_progenitor).  atol matters: AHF
    # prints merit to six decimal places, so the smallest merits carry a
    # relative rounding error of order 1e-4 all on their own.
    expected = links.shared.astype(float) ** 2 / (links.n_descendant * links.n_progenitor)
    np.testing.assert_allclose(links.merit, expected, rtol=1e-3, atol=1e-6)


def test_merger_link_ranks_follow_ahf_ordering(final_snapshot):
    """rank 0 is AHF's own best progenitor for each descendant."""
    links = ahf.read_merger_links(final_snapshot.ahf("croco"))
    starts = np.flatnonzero(np.r_[True, links.descendant[1:] != links.descendant[:-1]])
    assert np.all(links.rank[starts] == 0)
    # a progenitor can never share more particles than it has
    assert np.all(links.shared <= links.n_progenitor)
