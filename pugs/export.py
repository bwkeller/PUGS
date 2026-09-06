"""
Build the Parquet halo catalog directly from AHF output and the snapshots.

This is the whole pipeline that used to be ``build_tangos.sh``: discover the
snapshots, read AHF's halo tables, build the merger forest from AHF's own tree
files, measure the particle-derived properties with pynbody, and write one
Parquet file per snapshot.  Nothing is stored in a database along the way.

    pugs-export /path/to/NUGS128 --output NUGS128_catalog
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pyarrow as pa

from . import ahf, halo_properties, io, merger_forest
from .simulation import Simulation, find_simulation, halo_id

logger = logging.getLogger(__name__)

#: Minimum particles for a halo to enter the catalog.  AHF finds a long tail of
#: marginal objects; the merger forest still uses all of them, so trees are not
#: biased by this cut.
MIN_PARTICLES = 500

#: Columns copied verbatim from AHF_halos.  AHF's _gas and _star column
#: families are omitted: this is a collisionless volume and they are all zero.
AHF_COLUMNS: tuple[str, ...] = (
    "hostHalo",
    "numSubStruct",
    "Mhalo",
    "npart",
    "Xc",
    "Yc",
    "Zc",
    "VXc",
    "VYc",
    "VZc",
    "Rhalo",
    "Rmax",
    "r2",
    "mbp_offset",
    "com_offset",
    "Vmax",
    "v_esc",
    "sigV",
    "lambda",
    "lambdaE",
    "Lx",
    "Ly",
    "Lz",
    "b",
    "c",
    "Eax",
    "Eay",
    "Eaz",
    "Ebx",
    "Eby",
    "Ebz",
    "Ecx",
    "Ecy",
    "Ecz",
    "ovdens",
    "nbins",
    "fMhires",
    "Ekin",
    "Epot",
    "SurfP",
    "Phi0",
    "cNFW",
)

#: Identity columns, written for every halo.
IDENTITY_COLUMNS: tuple[str, ...] = (
    "halo_id",
    "finder_id",
    "snapshot",
    "snapshot_index",
    "z",
    "a",
)

#: Merger-tree and assembly-history columns, derived from the forest.  Unlike
#: the TANGOS pipeline these are written for every snapshot, not only z=0: the
#: forest is built once for the whole simulation, so they cost nothing extra.
TREE_COLUMNS: tuple[str, ...] = (
    "main_progenitor_id",
    "descendant_id",
    "n_progenitors",
    "N_mm",
    "z_lmm",
    "z25_mass",
    "z50_mass",
    "z75_mass",
)

#: Columns measured from particle data, one pynbody pass per halo.
MEASURED_COLUMNS: tuple[str, ...] = (
    "finder_mass",
    "max_radius",
    "shrink_center_x",
    "shrink_center_y",
    "shrink_center_z",
    "R200",
    "R500",
    "M200",
    "M500",
)

#: Columns derived arithmetically from AHF columns.
DERIVED_COLUMNS: tuple[str, ...] = ("virial_ratio",)

MASS_FRACTIONS = (0.25, 0.50, 0.75)


def export_simulation(
    folder: Path | str,
    output: Path | str,
    name: str | None = None,
    min_particles: int = MIN_PARTICLES,
    measure: bool = True,
) -> Path:
    """Build a catalog for the simulation in ``folder``.

    ``measure=False`` skips the particle-derived columns, leaving them null.
    That turns a run that is dominated by pynbody into one that takes seconds,
    which is useful when only the AHF columns or the trees are of interest.
    """
    output = Path(output)
    simulation = find_simulation(folder, name)
    logger.info("%s: %d snapshots", simulation.name, len(simulation))

    tables = {
        snapshot.extension: ahf.read_halo_table(snapshot.ahf("halos")) for snapshot in simulation
    }
    forest = merger_forest.build_forest(simulation, tables)
    logger.info(
        "merger forest: %d halos, %d links", len(forest), int((forest.descendant >= 0).sum())
    )

    schema = io.schema_for(_catalog_columns())
    rows_per_snapshot: dict[str, int] = {}
    failures: dict[str, int] = {}

    for snapshot in simulation:
        table, failed = _snapshot_table(
            simulation,
            snapshot,
            tables[snapshot.extension],
            forest,
            schema,
            min_particles=min_particles,
            measure=measure,
        )
        io.write_snapshot(
            table, output, snapshot=snapshot.extension, pugs_simulation=simulation.name
        )
        rows_per_snapshot[snapshot.extension] = table.num_rows
        if failed:
            failures[snapshot.extension] = failed
        logger.info(
            "%s: %d halos%s",
            snapshot.extension,
            table.num_rows,
            f" ({failed} unmeasured)" if failed else "",
        )

    io.write_provenance(
        output,
        io.build_provenance(
            simulation=simulation.name,
            source_directory=str(Path(folder).resolve()),
            n_snapshots=len(simulation),
            n_halos=sum(rows_per_snapshot.values()),
            rows_per_snapshot=rows_per_snapshot,
            final_snapshot=simulation.final.extension,
            min_particles=min_particles,
            particle_properties_measured=measure,
            unmeasured_halos=failures,
            forest=_forest_provenance(forest),
            cosmology=_cosmology(simulation),
            columns=sorted(schema.names),
        ),
    )
    return output


def _catalog_columns() -> list[str]:
    return list(IDENTITY_COLUMNS + MEASURED_COLUMNS + TREE_COLUMNS + AHF_COLUMNS + DERIVED_COLUMNS)


def _snapshot_table(
    simulation,
    snapshot,
    halo_table,
    forest,
    schema,
    min_particles,
    measure,
) -> tuple[pa.Table, int]:
    """Assemble one snapshot's halos into a table matching ``schema``."""
    finder = np.asarray(halo_table.column("ID").to_numpy(), dtype=np.int64)
    npart = np.asarray(halo_table.column("npart").to_numpy(), dtype=np.int64)

    selected = np.flatnonzero(npart >= min_particles)
    order = selected[np.argsort(finder[selected], kind="stable")]
    finder = finder[order]
    n = finder.size

    columns: dict[str, pa.Array] = {
        "halo_id": pa.array(halo_id(snapshot.index, finder), type=pa.int64()),
        "finder_id": pa.array(finder, type=pa.int64()),
        "snapshot": pa.array([snapshot.extension] * n, type=pa.string()),
        "snapshot_index": pa.array(np.full(n, snapshot.index, dtype=np.int64), type=pa.int64()),
        "z": pa.array(np.full(n, snapshot.redshift), type=pa.float64()),
        "a": pa.array(np.full(n, snapshot.scalefactor), type=pa.float64()),
    }

    for name in AHF_COLUMNS:
        if name not in halo_table.column_names:
            logger.warning("%s: AHF_halos has no column %s", snapshot.extension, name)
            continue
        values = np.asarray(halo_table.column(name).to_numpy())[order]
        columns[name] = pa.array(values, type=schema.field(name).type)

    columns.update(_tree_columns(forest, snapshot, finder, schema))
    columns.update(_derived_columns(columns, n))

    failed = 0
    if measure and n:
        measured, failed = _measure_snapshot(snapshot, halo_table, finder, schema)
        columns.update(measured)

    for field in schema:
        if field.name not in columns:
            columns[field.name] = pa.nulls(n, type=field.type)

    return pa.table({f.name: columns[f.name] for f in schema}, schema=schema), failed


def _tree_columns(forest, snapshot, finder, schema) -> dict[str, pa.Array]:
    """Merger-tree pointers and assembly history for the selected halos."""
    index = forest.index_of(halo_id(snapshot.index, finder))
    known = index >= 0
    safe = np.where(known, index, 0)

    def ids_of(pointer: np.ndarray) -> pa.Array:
        target = pointer[safe]
        present = known & (target >= 0)
        values = np.where(present, forest.halo_id[np.where(present, target, 0)], 0)
        return pa.array(values, type=pa.int64(), mask=~present)

    percentiles = forest.mass_percentile_redshifts(safe, MASS_FRACTIONS)
    missing = ~known

    columns = {
        "main_progenitor_id": ids_of(forest.main_progenitor),
        "descendant_id": ids_of(forest.descendant),
        "n_progenitors": pa.array(forest.n_progenitors[safe], type=pa.int64(), mask=missing),
        "N_mm": pa.array(forest.n_mm[safe], type=pa.int64(), mask=missing),
        "z_lmm": pa.array(forest.z_lmm[safe], type=pa.float64(), mask=missing),
    }
    for row, name in enumerate(("z25_mass", "z50_mass", "z75_mass")):
        columns[name] = pa.array(percentiles[row], type=pa.float64(), mask=missing)
    return columns


def _derived_columns(columns: dict[str, pa.Array], n: int) -> dict[str, pa.Array]:
    """Columns that are arithmetic on the AHF columns."""
    derived: dict[str, pa.Array] = {}
    if "Ekin" in columns and "Epot" in columns:
        kinetic = np.asarray(columns["Ekin"].to_numpy(zero_copy_only=False), dtype=np.float64)
        potential = np.asarray(columns["Epot"].to_numpy(zero_copy_only=False), dtype=np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = 2.0 * kinetic / potential
        bad = ~np.isfinite(ratio)
        derived["virial_ratio"] = pa.array(np.where(bad, 0.0, ratio), type=pa.float64(), mask=bad)
    return derived


def _measure_snapshot(snapshot, halo_table, finder, schema) -> tuple[dict[str, pa.Array], int]:
    """Measure the particle-derived properties for every selected halo.

    Membership comes through :func:`pugs.ahf.open_membership`, which seeks to
    each halo's block via the ``AHF_fpos`` index rather than parsing the whole
    particles file.  Only halos above the cut are ever read, and nothing larger
    than one halo is held at a time.
    """
    sim = snapshot.load()
    rho_crit = halo_properties.critical_density(sim)

    n = finder.size
    values = {name: np.zeros(n, dtype=np.float64) for name in MEASURED_COLUMNS}
    missing = np.ones(n, dtype=bool)
    failed = 0

    with ahf.open_membership(snapshot, halo_table) as membership:
        for row, halo in enumerate(finder):
            try:
                indices = membership.particles(int(halo))
                measurement = halo_properties.measure_halo(sim, indices, rho_crit).as_dict()
            except Exception as exc:
                failed += 1
                logger.warning("%s halo %d: %s", snapshot.extension, halo, exc)
                continue
            for name in MEASURED_COLUMNS:
                values[name][row] = measurement[name]
            missing[row] = False

    return {
        name: pa.array(values[name], type=schema.field(name).type, mask=missing)
        for name in MEASURED_COLUMNS
    }, failed


def _forest_provenance(forest) -> dict[str, Any]:
    return {
        "n_halos_all": int(len(forest)),
        "n_with_descendant": int((forest.descendant >= 0).sum()),
        "n_with_progenitor": int((forest.main_progenitor >= 0).sum()),
        "min_shared_fraction": merger_forest.MIN_SHARED_FRACTION,
        "major_merger_ratio": merger_forest.MAJOR_MERGER_RATIO,
        "mass_column": merger_forest.MASS_COLUMN,
    }


def _cosmology(simulation: Simulation) -> dict[str, Any]:
    """Cosmology from the final snapshot's header, recorded verbatim.

    The tipsy headers carry no omegaM0/omegaL0, so pynbody assumes defaults for
    them; ``inputs/planck_2018_CAMB.ini`` remains the authority for this
    volume's cosmology, not anything derived from a snapshot.
    """
    properties: dict[str, Any] = {}
    try:
        header = simulation.final._header
    except Exception as exc:
        logger.debug("could not read snapshot header: %s", exc)
        return properties
    for key, value in header.items():
        if isinstance(value, (str, bool, int, float)):
            properties[key] = value
        elif value is not None:
            properties[key] = str(value)
    return properties


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("folder", help="directory holding the snapshots and AHF output")
    parser.add_argument(
        "-o", "--output", default=None, help="catalog directory (default: ./<name>_catalog)"
    )
    parser.add_argument(
        "-n", "--name", default=None, help="simulation name (default: the folder's name)"
    )
    parser.add_argument(
        "--min-particles",
        type=int,
        default=MIN_PARTICLES,
        help=f"minimum particles per halo (default: {MIN_PARTICLES})",
    )
    parser.add_argument(
        "--no-measure",
        action="store_true",
        help="skip the particle-derived columns, leaving them null",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")
    name = args.name or Path(args.folder).resolve().name
    output = Path(args.output or f"{name}_catalog")
    export_simulation(
        args.folder,
        output,
        name=name,
        min_particles=args.min_particles,
        measure=not args.no_measure,
    )
    logger.info("wrote catalog to %s", output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
