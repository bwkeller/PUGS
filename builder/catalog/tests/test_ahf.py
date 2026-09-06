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


# ------------------------------------------------------ delimiter handling ---

_HEADER_FIELDS = ["ID(1)", "hostHalo(2)", "numSubStruct(3)", "Mhalo(4)", "npart(5)"]
_ROWS = [[0, -1, 6, 1.36606e14, 26156], [1, 0, 0, 5.5e12, 4337]]


def _write_ahf_halos(path, separator, pad=False):
    def field(value):
        text = f"{value:g}" if isinstance(value, float) else str(value)
        return text.rjust(12) if pad else text

    lines = ["#" + separator.join(_HEADER_FIELDS)]
    lines += [separator.join(field(v) for v in row) for row in _ROWS]
    path.write_text("\n".join(lines) + "\n")
    return path


def test_delimiter_is_detected_from_the_header():
    assert ahf._delimiter("#ID(1)\thostHalo(2)\n") == "\t"
    assert ahf._delimiter("#ID(1) hostHalo(2)\n") == " "


@pytest.mark.parametrize(
    "separator,pad",
    [("\t", True), ("\t", False), (" ", False)],
    ids=["tab-padded", "tab-plain", "space"],
)
def test_reads_both_ahf_separator_conventions(tmp_path, separator, pad):
    """
    AHF is not consistent between runs: the NUGS128 catalogs are tab-separated
    with the columns space-padded, the NUGS2048 ones single-space with no tabs
    anywhere. Assuming either one makes the reader return a single column of
    text on the other, which then fails much later as a missing column.
    """
    path = _write_ahf_halos(tmp_path / "sim.z0.000.AHF_halos", separator, pad)
    table = ahf.read_halo_table(path)

    assert table.column_names == ["ID", "hostHalo", "numSubStruct", "Mhalo", "npart"]
    assert table.num_rows == len(_ROWS)
    assert table.column("ID").to_pylist() == [0, 1]
    assert table.column("npart").to_pylist() == [26156, 4337]


def test_header_without_hash_is_rejected(tmp_path):
    path = tmp_path / "sim.z0.000.AHF_halos"
    path.write_text("ID(1) npart(5)\n0 10\n")
    with pytest.raises(ValueError, match="header line"):
        ahf.read_halo_table(path)


def test_no_host_sentinel_is_not_assumed(final_snapshot):
    """
    The 'no host' value depends on where the run numbers halo ids from: -1
    where they start at 0 (NUGS128), 0 where they start at 1 (NUGS2048).
    Anything keying on hostHalo must compare against the smallest id present.
    """
    table = ahf.read_halo_table(final_snapshot.ahf("halos"))
    ids = np.asarray(table.column("ID").to_numpy())
    host = np.asarray(table.column("hostHalo").to_numpy())
    assert (host < ids.min()).sum() > 0, "expected some halos with no host"
    assert set(host[host >= ids.min()]) <= set(ids.tolist())


# ------------------------------------------------- fpos-indexed membership ---


def test_fpos_has_one_offset_per_halo(final_snapshot):
    offsets = ahf.read_fpos(final_snapshot.ahf("fpos"))
    table = ahf.read_halo_table(final_snapshot.ahf("halos"))
    assert len(offsets) == table.num_rows
    assert np.all(np.diff(offsets) > 0), "offsets must be strictly increasing"


def test_fpos_rejects_a_non_monotonic_index(tmp_path):
    path = tmp_path / "sim.z0.000.AHF_fpos"
    path.write_text("100\n50\n200\n")
    with pytest.raises(ValueError, match="increasing"):
        ahf.read_fpos(path)


def test_indexed_membership_matches_the_full_parse(final_snapshot):
    """
    The index is only useful if it lands on exactly the same particles. Checked
    for every halo, not a sample: an off-by-one in the offsets would give a
    plausible-looking but wrong particle set.
    """
    table = ahf.read_halo_table(final_snapshot.ahf("halos"))
    reference = ahf.read_particle_membership(final_snapshot.ahf("particles"))

    with ahf.open_membership(final_snapshot, table) as indexed:
        assert isinstance(indexed, ahf.IndexedMembership)
        for halo in indexed.halo_id:
            got = indexed.particles(int(halo))
            expected = reference.particles(int(halo))
            assert np.array_equal(np.sort(got), np.sort(expected)), halo


def test_indexed_membership_respects_declared_counts(final_snapshot):
    table = ahf.read_halo_table(final_snapshot.ahf("halos"))
    declared = dict(zip(table.column("ID").to_pylist(), table.column("npart").to_pylist()))
    with ahf.open_membership(final_snapshot, table) as indexed:
        for halo in list(declared)[:100]:
            assert len(indexed.particles(halo)) == declared[halo]


def test_block_delimiter_is_detected_per_file(final_snapshot):
    """
    Inside a block the separator is uniform, but it is not the same across
    catalogs: NUGS128 uses tabs between the particle columns while NUGS2048
    uses spaces. The block headers use yet another convention, which is why
    seeking past them matters.
    """
    table = ahf.read_halo_table(final_snapshot.ahf("halos"))
    with ahf.open_membership(final_snapshot, table) as indexed:
        assert indexed.delimiter in ("\t", " ")


def test_indexed_membership_rejects_unknown_halo(final_snapshot):
    table = ahf.read_halo_table(final_snapshot.ahf("halos"))
    with ahf.open_membership(final_snapshot, table) as indexed:
        with pytest.raises(KeyError):
            indexed.particles(10**9)


def test_missing_fpos_falls_back_and_warns(final_snapshot, tmp_path, caplog):
    """Without the index there is no choice but to parse the whole file."""
    import shutil

    stem = tmp_path / "DM.z0.000"
    for suffix in ("halos", "particles"):
        shutil.copy(final_snapshot.ahf(suffix), f"{stem}.AHF_{suffix}")

    from pugs.simulation import Snapshot

    snapshot = Snapshot(
        index=0,
        extension="DM",
        path=final_snapshot.path,
        ahf_stem=stem,
        catalog_redshift=0.0,
    )
    assert not snapshot.has_ahf("fpos")

    table = ahf.read_halo_table(snapshot.ahf("halos"))
    with caplog.at_level("WARNING"):
        membership = ahf.open_membership(snapshot, table)
    assert isinstance(membership, ahf.ParticleMembership)
    assert "AHF_fpos" in caplog.text


def test_mismatched_fpos_length_is_an_error(final_snapshot, tmp_path):
    table = ahf.read_halo_table(final_snapshot.ahf("halos"))
    short = tmp_path / "short.AHF_fpos"
    short.write_text("\n".join(str(v) for v in [10, 20, 30]) + "\n")
    with pytest.raises(ValueError, match="offsets for"):
        ahf.open_indexed_membership(
            final_snapshot.ahf("particles"),
            short,
            table.column("ID").to_numpy(),
            table.column("npart").to_numpy(),
        )
