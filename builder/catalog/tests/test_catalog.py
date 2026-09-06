"""
Tests for the Parquet halo catalog: the export, and the metadata contract that
makes an exported catalog self-describing.
"""

import json

import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from numpy.testing import assert_allclose

from pugs import ahf, export, io
from pugs.simulation import split_halo_id

SIM = "NUGS128"


def populated(simulation, min_particles=export.MIN_PARTICLES):
    """Snapshots that keep at least one halo after the particle cut."""
    kept = {}
    for snapshot in simulation:
        npart = ahf.read_halo_table(snapshot.ahf("halos")).column("npart").to_numpy()
        kept[snapshot.extension] = int((np.asarray(npart) >= min_particles).sum())
    return kept


# ------------------------------------------------------------ data dictionary --


def test_every_documented_column_is_complete():
    for name, entry in io.dictionary()["columns"].items():
        assert entry.keys() >= {"dtype", "unit", "system", "description"}, name
        assert entry["description"].strip(), name
        assert entry["system"] in {"physical", "finder", "pugs"}, name


def test_exported_columns_are_all_documented(catalog):
    assert set(catalog.columns) <= set(io.documented_columns())


# ----------------------------------------------------------------- structure --


def test_one_file_per_snapshot(catalog_dir, simulation):
    assert len(io.catalog_files(catalog_dir)) == len(simulation)


def test_schema_identical_across_snapshots(catalog_dir):
    """Every file carries the full column set, so the directory unions without
    `union_by_name`."""
    schemas = [pq.read_schema(f) for f in io.catalog_files(catalog_dir)]
    for schema in schemas[1:]:
        assert schema.names == schemas[0].names
        assert schema.types == schemas[0].types


def test_row_counts_match_the_particle_cut(catalog_dir, simulation):
    expected = populated(simulation)
    for snapshot in simulation:
        table = pq.read_table(catalog_dir / f"halos_{snapshot.extension}.parquet")
        assert table.num_rows == expected[snapshot.extension], snapshot.extension


def test_empty_snapshots_still_get_a_file(catalog_dir, simulation):
    expected = populated(simulation)
    empty = [name for name, count in expected.items() if count == 0]
    assert empty, "expected the particle cut to leave some snapshots empty"
    for name in empty:
        table = pq.read_table(catalog_dir / f"halos_{name}.parquet")
        assert table.num_rows == 0
        assert "M200" in table.column_names


def test_halo_ids_are_unique_and_decodable(catalog):
    halo_id = catalog.table.column("halo_id").to_numpy()
    assert len(np.unique(halo_id)) == len(halo_id)
    snapshot_index, finder_id = split_halo_id(halo_id)
    assert np.array_equal(snapshot_index, catalog.table.column("snapshot_index").to_numpy())
    assert np.array_equal(finder_id, catalog.table.column("finder_id").to_numpy())


def test_snapshot_and_redshift_are_columns(catalog, simulation):
    """Per-snapshot values must be columns: pyarrow takes footer metadata from
    a single file of a multi-file dataset."""
    table = catalog.table
    assert table.column("snapshot").null_count == 0
    z = table.column("z").to_numpy()
    a = table.column("a").to_numpy()
    assert_allclose(a, 1.0 / (1.0 + z))
    assert len(np.unique(z)) == sum(1 for c in populated(simulation).values() if c)


# -------------------------------------------------------------------- values --


def test_ahf_columns_round_trip(catalog_dir, simulation):
    """AHF columns must reach the catalog unchanged, matched by halo id."""
    for snapshot in simulation:
        source = ahf.read_halo_table(snapshot.ahf("halos"))
        table = pq.read_table(catalog_dir / f"halos_{snapshot.extension}.parquet")
        if table.num_rows == 0:
            continue
        by_id = dict(zip(source.column("ID").to_pylist(), range(source.num_rows)))
        rows = [by_id[f] for f in table.column("finder_id").to_pylist()]
        for name in ("Mhalo", "npart", "hostHalo", "Rhalo", "Vmax", "cNFW", "Xc"):
            expected = np.asarray(source.column(name).to_numpy())[rows]
            got = np.asarray(table.column(name).to_numpy(zero_copy_only=False))
            assert_allclose(got.astype(float), expected.astype(float), rtol=1e-12, err_msg=name)


def test_substructure_survives_the_export(catalog):
    """
    hostHalo is real data in the AHF output. The TANGOS pipeline imported
    nothing for it, so an all-null column here would be a regression to that.
    """
    host = catalog.table.column("hostHalo")
    assert host.null_count < len(host)
    values = np.asarray(host.to_numpy(zero_copy_only=False))
    assert (values >= 0).sum() > 0, "expected subhalos"
    assert (values == -1).sum() > 0, "expected field halos"


def test_measured_columns_recompute(catalog_dir, simulation, snapshot_data):
    """Spot-check the particle-derived columns against a direct measurement."""
    from pugs import halo_properties

    sim, membership, rho = snapshot_data
    snapshot = simulation.final
    table = pq.read_table(catalog_dir / f"halos_{snapshot.extension}.parquet")
    finder = table.column("finder_id").to_pylist()[:5]
    for row, halo in enumerate(finder):
        measured = halo_properties.measure_halo(sim, membership.particles(halo), rho).as_dict()
        for name in ("finder_mass", "max_radius", "R200", "R500", "M200", "M500"):
            got = float(table.column(name)[row].as_py())
            assert_allclose(got, measured[name], rtol=1e-10, err_msg=f"{halo} {name}")


def test_virial_ratio_follows_the_energies(catalog):
    table = catalog.table
    ekin = np.asarray(table.column("Ekin").to_numpy(zero_copy_only=False), dtype=float)
    epot = np.asarray(table.column("Epot").to_numpy(zero_copy_only=False), dtype=float)
    ratio = np.asarray(table.column("virial_ratio").to_numpy(zero_copy_only=False), dtype=float)
    good = np.isfinite(ratio) & (epot != 0)
    assert_allclose(ratio[good], (2 * ekin / epot)[good], rtol=1e-12)


# ------------------------------------------------------- merger tree columns --


def test_tree_pointers_reference_real_halos(catalog):
    """
    A descendant may sit below the particle cut and so be absent from the
    catalog, but its snapshot must still be the next one along.
    """
    table = catalog.table
    halo_id = table.column("halo_id").to_numpy()
    snapshot_index = table.column("snapshot_index").to_numpy()

    for column, offset in (("main_progenitor_id", -1), ("descendant_id", +1)):
        pointer = table.column(column)
        present = ~np.asarray(pointer.is_null())
        target = np.asarray(pointer.to_numpy(zero_copy_only=False))[present]
        target_snapshot, _ = split_halo_id(target.astype(np.int64))
        assert np.array_equal(target_snapshot, snapshot_index[present] + offset), column
        assert present.sum() > 0, column
    assert len(np.unique(halo_id)) == len(halo_id)


def test_assembly_history_is_written_for_every_snapshot(catalog_dir, simulation):
    """
    Unlike the TANGOS pipeline, which only wrote these at z=0, the forest is
    built once for the whole simulation so they cost nothing extra.
    """
    for snapshot in simulation:
        table = pq.read_table(
            catalog_dir / f"halos_{snapshot.extension}.parquet",
            columns=["N_mm", "z_lmm", "z50_mass"],
        )
        if table.num_rows == 0:
            continue
        for name in ("N_mm", "z_lmm", "z50_mass"):
            assert table.column(name).null_count == 0, f"{snapshot.extension} {name}"


def test_assembly_redshifts_are_ordered(catalog):
    table = catalog.table
    z25, z50, z75 = (
        np.asarray(table.column(n).to_numpy(zero_copy_only=False), dtype=float)
        for n in ("z25_mass", "z50_mass", "z75_mass")
    )
    assert np.all(z25 >= z50)
    assert np.all(z50 >= z75)


def test_merger_sentinels_agree(catalog):
    table = catalog.table
    n_mm = np.asarray(table.column("N_mm").to_numpy(zero_copy_only=False), dtype=float)
    z_lmm = np.asarray(table.column("z_lmm").to_numpy(zero_copy_only=False), dtype=float)
    assert np.all(z_lmm[n_mm == 0] == -1.0)
    assert np.all(z_lmm[n_mm > 0] > 0.0)


# ------------------------------------------------------------------ metadata --


def test_units_travel_with_the_data(catalog):
    entries = io.dictionary()["columns"]
    for name in catalog.columns:
        assert catalog.unit(name) == entries[name]["unit"], name
        assert catalog.system(name) == entries[name]["system"], name
    assert catalog.unit("M200") == "Msol"
    assert catalog.unit("Mhalo") == "Msol / littleh"
    assert catalog.system("M200") == "physical"
    assert catalog.system("Mhalo") == "finder"


def test_non_arrow_readers_can_see_units(catalog_dir):
    """Arrow hides field metadata inside a base64 ARROW:schema blob that DuckDB
    cannot read, so units are mirrored into a file-level JSON entry."""
    path = str(io.catalog_files(catalog_dir)[0])
    rows = duckdb.sql(f"select key, value from parquet_kv_metadata('{path}')").fetchall()
    footer = {k.decode(): v.decode() for k, v in rows}
    assert footer["pugs_schema_version"] == str(io.schema_version())
    mirrored = json.loads(footer["columns"])
    assert mirrored["M200"]["unit"] == "Msol"
    assert mirrored["Mhalo"]["unit"] == "Msol / littleh"


def test_duckdb_can_query_the_whole_catalog(catalog_dir, simulation):
    glob = str(catalog_dir / io.CATALOG_GLOB)
    expected = populated(simulation)
    n_rows, n_snaps = duckdb.sql(
        f"select count(*), count(distinct snapshot) from '{glob}'"
    ).fetchone()
    assert n_rows == sum(expected.values())
    assert n_snaps == sum(1 for count in expected.values() if count)
    assert duckdb.sql(f"select count(*) from '{glob}' where M200 > 1e12").fetchone()[0] > 0


def test_provenance_sidecar(catalog_dir, catalog, simulation):
    sidecar = json.loads((catalog_dir / io.PROVENANCE_FILENAME).read_text())
    assert sidecar["simulation"] == SIM
    assert sidecar["pugs_schema_version"] == io.schema_version()
    assert sidecar["n_snapshots"] == len(simulation)
    assert sidecar["min_particles"] == export.MIN_PARTICLES
    assert sidecar["particle_properties_measured"] is True
    assert sidecar["versions"]["pynbody"]
    assert "tangos" not in sidecar["versions"]
    assert sidecar["forest"]["n_halos_all"] > sidecar["n_halos"]
    assert sidecar["cosmology"]
    assert catalog.provenance == sidecar


def test_no_halo_failed_to_measure(catalog_dir):
    sidecar = json.loads((catalog_dir / io.PROVENANCE_FILENAME).read_text())
    assert sidecar["unmeasured_halos"] == {}


def test_column_subset_reads_only_what_was_asked(catalog_dir):
    subset = io.read_catalog(catalog_dir, columns=["halo_id", "M200"])
    assert subset.columns == ["halo_id", "M200"]
    assert subset.unit("M200") == "Msol"


# --------------------------------------------------------- failure behaviour --


def test_undocumented_column_is_rejected(tmp_path):
    table = pa.table({"halo_id": pa.array([1, 2], pa.int64()), "not_a_column": [1.0, 2.0]})
    with pytest.raises(io.UndocumentedColumnError, match="not_a_column"):
        io.write_snapshot(table, tmp_path, "DM128.00512")


def test_version_mismatch_is_an_error(tmp_path, monkeypatch):
    table = pa.table({"halo_id": pa.array([1, 2], pa.int64())})
    io.write_snapshot(table, tmp_path, "DM128.00512")
    io.write_provenance(tmp_path, io.build_provenance(simulation=SIM))
    io.read_catalog(tmp_path)  # matching version is fine

    monkeypatch.setattr(io, "schema_version", lambda: 999)
    with pytest.raises(io.SchemaVersionError, match="schema version"):
        io.read_catalog(tmp_path)


def test_unstamped_file_is_rejected(tmp_path):
    pq.write_table(pa.table({"halo_id": pa.array([1], pa.int64())}), tmp_path / "halos_x.parquet")
    with pytest.raises(io.SchemaVersionError, match="not written by pugs.io"):
        io.read_catalog(tmp_path)


# ------------------------------------------------------- dependency revisions --


def test_provenance_records_dependency_revisions(catalog):
    """
    A version number cannot identify a build: the pynbody fork carrying the
    deterministic-units fix reports the same version as the release it was
    branched from. The commit is what distinguishes them.
    """
    revisions = catalog.provenance["revisions"]
    assert "pynbody" in revisions, revisions
    entry = revisions["pynbody"]
    assert len(entry["commit"]) == 40
    assert entry["resolved_from"] in {"install record", "working tree"}


def test_revisions_do_not_borrow_the_enclosing_repo(catalog):
    """
    A virtualenv usually sits inside a source checkout, so walking up from
    site-packages finds *that project's* git repo. Every dependency would then
    be stamped with the PUGS commit, which looks plausible and is wrong.
    """
    revisions = catalog.provenance["revisions"]
    pugs_sha = catalog.provenance["git_sha"]
    if pugs_sha:
        borrowed = [n for n, r in revisions.items() if r.get("commit") == pugs_sha]
        assert not borrowed, f"{borrowed} were attributed the PUGS commit"


def test_installed_copies_are_not_treated_as_checkouts():
    from pathlib import Path

    assert io._is_installed_copy(Path("/x/.venv/lib/python3.12/site-packages/numpy"))
    assert io._is_installed_copy(Path("/x/lib/python3/dist-packages/numpy"))
    assert not io._is_installed_copy(Path("/home/me/code/pynbody/pynbody"))


def test_vcs_install_record_is_preferred(tmp_path, monkeypatch):
    """When pip installed from a git URL it records the exact commit (PEP 610)."""

    class FakeDistribution:
        def read_text(self, name):
            assert name == "direct_url.json"
            return json.dumps(
                {
                    "url": "https://github.com/bwkeller/pynbody.git",
                    "vcs_info": {"vcs": "git", "commit_id": "a" * 40, "requested_revision": "pugs"},
                }
            )

        def locate_file(self, name):
            return tmp_path / name

    monkeypatch.setattr(io.importlib_metadata, "distribution", lambda name: FakeDistribution())
    revision = io._package_revision("pynbody")
    assert revision["commit"] == "a" * 40
    assert revision["resolved_from"] == "install record"
    assert revision["url"].endswith("pynbody.git")


def test_missing_package_yields_no_revision(monkeypatch):
    def missing(name):
        raise io.importlib_metadata.PackageNotFoundError(name)

    monkeypatch.setattr(io.importlib_metadata, "distribution", missing)
    assert io._package_revision("nope") is None
