"""
Tests for the stored particle-id shells.

The shells are checked against the snapshot rather than against each other:
membership is recomputed from the particle positions, so a bug in the shell
assignment cannot agree with itself.
"""

import numpy as np
import pyarrow.parquet as pq
import pytest
from numpy.testing import assert_allclose

from pugs import halo_properties, io, particle_ids


@pytest.fixture(scope="module")
def shells(shells_dir):
    return io.read_shells(shells_dir)


@pytest.fixture(scope="module")
def by_halo(shells):
    """{halo_id: {shell: ids}}."""
    table = shells.table
    grouped: dict[int, dict[int, np.ndarray]] = {}
    for halo, shell, ids in zip(
        table.column("halo_id").to_pylist(),
        table.column("shell").to_pylist(),
        table.column("particle_ids").to_pylist(),
    ):
        grouped.setdefault(halo, {})[shell] = np.asarray(ids, dtype=np.int64)
    return grouped


# ------------------------------------------------------------- shell grid ----


def test_shell_grid_is_half_rvir_out_to_five():
    assert particle_ids.SHELL_STEP == 0.5
    assert particle_ids.SHELL_MAX == 5.0
    assert particle_ids.N_SHELLS == 10
    assert_allclose(particle_ids.SHELL_EDGES, np.arange(0.0, 5.5, 0.5))


def test_shell_of_assigns_edges_correctly():
    """The innermost shell is a sphere; anything at or beyond 5 R_vir is out."""
    rvir = 100.0
    radii = np.array([0.0, 49.9, 50.0, 99.9, 100.0, 499.9, 500.0, 1e4])
    expected = [0, 0, 1, 1, 2, 9, particle_ids.N_SHELLS, particle_ids.N_SHELLS]
    assert list(particle_ids.shell_of(radii, rvir)) == expected


def test_every_halo_gets_the_full_shell_grid(shells, by_halo):
    for halo, grouped in by_halo.items():
        assert sorted(grouped) == list(range(particle_ids.N_SHELLS)), halo


def test_edges_match_the_grid(shells):
    table = shells.table
    pairs = set(
        zip(table.column("r_inner_rvir").to_pylist(), table.column("r_outer_rvir").to_pylist())
    )
    expected = {(round(e, 6), round(e + 0.5, 6)) for e in np.arange(0.0, 5.0, 0.5)}
    assert pairs == expected


def test_physical_edges_scale_with_the_reference_radius(shells):
    table = shells.table
    reference = np.asarray(table.column("reference_radius").to_numpy(), dtype=float)
    for edge, physical in (("r_inner_rvir", "r_inner"), ("r_outer_rvir", "r_outer")):
        fraction = np.asarray(table.column(edge).to_numpy(), dtype=float)
        assert_allclose(
            np.asarray(table.column(physical).to_numpy(), dtype=float),
            fraction * reference,
            rtol=1e-12,
        )


# ------------------------------------------------- the delta-encoding premise --


def test_ids_are_ascending_within_every_shell(shells):
    """
    The whole point of bucketing by radius instead of sorting by it: the ids
    stay in ascending order, so consecutive values differ by very little and
    DELTA_BINARY_PACKED has something to compress. Radial sorting would
    scramble them.
    """
    for ids in shells.table.column("particle_ids").to_pylist():
        values = np.asarray(ids, dtype=np.int64)
        if values.size > 1:
            assert np.all(np.diff(values) > 0)


def test_ids_are_delta_encoded_on_disk(shells_dir):
    """The encoding is pinned in the writer; confirm Parquet actually used it."""
    path = io.catalog_files(shells_dir, io.SHELL_GLOB)[0]
    metadata = pq.ParquetFile(path).metadata
    index = metadata.schema.names.index("element")
    encodings = set()
    for group in range(metadata.num_row_groups):
        encodings.update(metadata.row_group(group).column(index).encodings)
    assert "DELTA_BINARY_PACKED" in encodings, encodings


def test_storage_beats_raw_int64(shells, shells_dir):
    """A weak bound, so the test fails if the encoding silently stops working."""
    n_ids = int(np.sum(shells.table.column("n_particles").to_numpy()))
    on_disk = sum(f.stat().st_size for f in io.catalog_files(shells_dir, io.SHELL_GLOB))
    assert n_ids * 8 / on_disk > 3.0


# --------------------------------------------------------------- membership --


def test_n_particles_matches_the_id_lists(shells):
    table = shells.table
    counts = np.asarray(table.column("n_particles").to_numpy())
    lengths = np.array([len(v) for v in table.column("particle_ids").to_pylist()])
    assert np.array_equal(counts, lengths)


def test_shells_partition_the_region(by_halo, catalog_dir, simulation, snapshot_data):
    """
    The ten shells must be disjoint and must union to exactly the 5 R_vir
    sphere -- checked against a fresh region selection from the snapshot.
    """
    import pynbody

    sim, _, _ = snapshot_data
    halos = io.read_catalog(
        catalog_dir / f"halos_{simulation.final.extension}.parquet",
        columns=[
            "halo_id",
            particle_ids.REFERENCE_RADIUS,
            "shrink_center_x",
            "shrink_center_y",
            "shrink_center_z",
        ],
    ).table

    for row in range(0, halos.num_rows, 40):
        halo = halos.column("halo_id")[row].as_py()
        radius = halos.column(particle_ids.REFERENCE_RADIUS)[row].as_py()
        centre = np.array([halos.column(f"shrink_center_{a}")[row].as_py() for a in "xyz"])

        stored = np.concatenate([by_halo[halo][s] for s in range(particle_ids.N_SHELLS)])
        assert len(stored) == len(np.unique(stored)), f"halo {halo}: shells overlap"

        region = sim[pynbody.filt.Sphere(particle_ids.SHELL_MAX * radius, centre)]
        expected = np.asarray(region.get_index_list(sim.ancestor))
        with halo_properties.recentred(region, centre):
            radii = np.sqrt((np.asarray(region["pos"], dtype=np.float64) ** 2).sum(axis=1))
        expected = expected[radii < particle_ids.SHELL_MAX * radius]
        assert np.array_equal(np.sort(stored), np.sort(expected)), f"halo {halo}"


def test_particles_lie_inside_their_shell(by_halo, catalog_dir, simulation, snapshot_data):
    sim, _, _ = snapshot_data
    halos = io.read_catalog(
        catalog_dir / f"halos_{simulation.final.extension}.parquet",
        columns=[
            "halo_id",
            particle_ids.REFERENCE_RADIUS,
            "shrink_center_x",
            "shrink_center_y",
            "shrink_center_z",
        ],
    ).table

    for row in range(0, halos.num_rows, 60):
        halo = halos.column("halo_id")[row].as_py()
        radius = halos.column(particle_ids.REFERENCE_RADIUS)[row].as_py()
        centre = np.array([halos.column(f"shrink_center_{a}")[row].as_py() for a in "xyz"])
        for shell, ids in by_halo[halo].items():
            if ids.size == 0:
                continue
            subset = sim[ids]
            with halo_properties.recentred(subset, centre):
                scaled = (
                    np.sqrt((np.asarray(subset["pos"], dtype=np.float64) ** 2).sum(axis=1)) / radius
                )
            assert scaled.min() >= shell * 0.5 - 1e-6
            assert scaled.max() < (shell + 1) * 0.5 + 1e-6


def test_only_the_final_snapshot_is_written(shells_dir, simulation):
    files = io.catalog_files(shells_dir, io.SHELL_GLOB)
    assert [f.name for f in files] == [f"shells_{simulation.final.extension}.parquet"]


# ------------------------------------------------------------- metadata -----


def test_shell_columns_carry_units(shells):
    assert shells.unit("reference_radius") == "kpc"
    assert shells.unit("r_outer") == "kpc"
    assert shells.system("r_outer") == "physical"
    assert shells.unit("r_outer_rvir") == "dimensionless"


def test_provenance_records_the_shell_settings(shells):
    recorded = shells.provenance["particle_ids"]
    assert recorded["shell_step_rvir"] == particle_ids.SHELL_STEP
    assert recorded["shell_max_rvir"] == particle_ids.SHELL_MAX
    assert recorded["n_shells"] == particle_ids.N_SHELLS
    assert recorded["reference_radius"] == particle_ids.REFERENCE_RADIUS
    assert recorded["snapshots"]


def test_reading_geometry_does_not_need_the_id_column(shells_dir):
    """n_particles exists so shell sizes can be inspected cheaply."""
    lean = io.read_shells(shells_dir, columns=["halo_id", "shell", "n_particles"])
    assert "particle_ids" not in lean.columns
    assert lean.table.num_rows > 0
