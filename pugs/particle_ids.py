"""
Particle-id regions for the halos of the final snapshot.

For each halo the particles within ``SHELL_MAX`` virial radii are split into
concentric shells of ``SHELL_STEP`` R_vir each, and their snapshot indices are
stored one row per (halo, shell).  The innermost shell is a sphere rather than
a shell, since its inner edge is zero.

A zoom region out to any multiple of R_vir is then the concatenation of the
shells below it -- no decompression of a larger region, and no per-halo blob to
unpack.

Why shells rather than one radially sorted list
-----------------------------------------------
The obvious layout is a single list per halo sorted by radius, which makes any
sub-region a prefix.  But radial order scrambles the ids, and an int64 column
in random order compresses badly: there is nothing for a delta encoder to bite
on.

Bucketing by radius instead of sorting by it keeps the ids of each shell in
their natural ascending order -- the region selection already returns them
sorted -- so consecutive values differ by very little.  In the NUGS128 z=0
catalog the median gap is 1-3 inside R_vir and around 100 in the outer shells,
against ids spanning millions.  Parquet's DELTA_BINARY_PACKED encoding then
stores the differences rather than the values.

Measured over the 40 most massive z=0 halos of NUGS128 (534k ids):

    raw int64                    4.27 MB   1.0x
    one radially sorted blob     0.88 MB   4.8x   (the old zlib scheme)
    shells, DELTA + zstd         0.55 MB   7.8x

The ids stay ordinary Parquet int64 values, readable by anything, with the
compression handled by the format rather than by a custom codec.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pyarrow as pa

from . import io
from .simulation import find_simulation

logger = logging.getLogger(__name__)

#: Width of each shell, in units of the reference radius.
SHELL_STEP = 0.5

#: Outer edge of the outermost shell, in units of the reference radius.
SHELL_MAX = 5.0

#: Shell edges: 0, 0.5, ... 5.0 R_vir.
SHELL_EDGES = np.round(np.arange(0.0, SHELL_MAX + SHELL_STEP / 2, SHELL_STEP), 6)

#: Number of shells (one fewer than the number of edges).
N_SHELLS = len(SHELL_EDGES) - 1

#: Catalog column used as the reference radius, "R_vir".
#:
#: R200 is the halo's actual virial radius.  `max_radius` -- the distance to
#: the outermost halo-finder particle -- is the other defensible choice and is
#: what the older TANGOS-based pipeline used; the two agree to a few parts in a
#: thousand for the NUGS128 halos.  Pass ``reference="max_radius"`` to switch.
REFERENCE_RADIUS = "R200"

SHELL_COLUMNS: tuple[str, ...] = (
    "halo_id",
    "snapshot",
    "shell",
    "reference_radius",
    "r_inner_rvir",
    "r_outer_rvir",
    "r_inner",
    "r_outer",
    "n_particles",
    "particle_ids",
)

#: Storing the ids as deltas is the entire point of the shell layout, so the
#: encoding is pinned rather than left to pyarrow's defaults.  Dictionary
#: encoding has to be off: it would replace the values with dictionary indices
#: and there would be nothing left for the delta encoder to compress.
ID_COLUMN_ENCODING = {"particle_ids.list.element": "DELTA_BINARY_PACKED"}


def shell_of(radii: np.ndarray, reference_radius: float) -> np.ndarray:
    """Shell index of each radius, ``N_SHELLS`` for anything beyond the last."""
    with np.errstate(divide="ignore", invalid="ignore"):
        index = np.floor(radii / (SHELL_STEP * reference_radius)).astype(np.int64)
    return np.where((index < 0) | (radii >= SHELL_MAX * reference_radius), N_SHELLS, index)


def halo_shells(sim, centre, reference_radius: float) -> list[np.ndarray]:
    """Snapshot indices of each shell around one halo, innermost first.

    The returned arrays are in ascending index order.  That is not enforced by
    sorting -- ``Sphere`` already yields sorted indices, and bucketing by
    radius preserves the order within each bucket.  It is what makes the ids
    delta-encodable.
    """
    import pynbody

    if not np.isfinite(reference_radius) or reference_radius <= 0:
        raise ValueError(f"reference radius must be positive and finite, got {reference_radius}")

    # Sphere is periodic-aware, so a halo near a box face still gets its whole
    # region rather than a truncated one.
    region = sim[pynbody.filt.Sphere(SHELL_MAX * reference_radius, centre)]
    indices = np.asarray(region.get_index_list(sim.ancestor), dtype=np.int64)
    if indices.size == 0:
        return [np.empty(0, dtype=np.int64) for _ in range(N_SHELLS)]

    with io_recentred(region, centre) as centred:
        radii = np.sqrt((np.asarray(centred["pos"], dtype=np.float64) ** 2).sum(axis=1))

    membership = shell_of(radii, reference_radius)
    return [indices[membership == shell] for shell in range(N_SHELLS)]


def io_recentred(region, centre):
    """Local alias so the import stays inside :mod:`pugs.halo_properties`."""
    from .halo_properties import recentred

    return recentred(region, centre)


def build_shell_table(snapshot, halos, sim, reference: str, schema: pa.Schema):
    """Shell rows for every halo of one snapshot.

    ``halos`` is the catalog table for that snapshot, which supplies the centre
    and reference radius rather than having them recomputed here.
    """
    halo_id = np.asarray(halos.column("halo_id").to_numpy(), dtype=np.int64)
    radius = np.asarray(halos.column(reference).to_numpy(zero_copy_only=False), dtype=np.float64)
    centres = np.column_stack(
        [
            np.asarray(halos.column(f"shrink_center_{axis}").to_numpy(zero_copy_only=False))
            for axis in ("x", "y", "z")
        ]
    )

    out_halo, out_shell, out_reference, out_counts, out_ids = [], [], [], [], []
    skipped = 0

    for row in range(len(halo_id)):
        if not np.isfinite(radius[row]) or radius[row] <= 0:
            skipped += 1
            logger.warning(
                "%s halo %d: %s is %r; skipping",
                snapshot.extension,
                halo_id[row],
                reference,
                radius[row],
            )
            continue
        try:
            shells = halo_shells(sim, centres[row], float(radius[row]))
        except Exception as exc:
            skipped += 1
            logger.warning("%s halo %d: %s", snapshot.extension, halo_id[row], exc)
            continue

        for shell, ids in enumerate(shells):
            out_halo.append(int(halo_id[row]))
            out_shell.append(shell)
            out_reference.append(float(radius[row]))
            out_counts.append(len(ids))
            out_ids.append(ids)

    inner = np.array([SHELL_EDGES[s] for s in out_shell], dtype=np.float64)
    outer = np.array([SHELL_EDGES[s + 1] for s in out_shell], dtype=np.float64)
    reference_radius = np.asarray(out_reference, dtype=np.float64)

    columns = {
        "halo_id": pa.array(out_halo, type=pa.int64()),
        "snapshot": pa.array([snapshot.extension] * len(out_halo), type=pa.string()),
        "shell": pa.array(out_shell, type=pa.int64()),
        "reference_radius": pa.array(reference_radius, type=pa.float64()),
        "r_inner_rvir": pa.array(inner, type=pa.float64()),
        "r_outer_rvir": pa.array(outer, type=pa.float64()),
        "r_inner": pa.array(inner * reference_radius, type=pa.float64()),
        "r_outer": pa.array(outer * reference_radius, type=pa.float64()),
        "n_particles": pa.array(out_counts, type=pa.int64()),
        "particle_ids": pa.array(out_ids, type=pa.list_(pa.int64())),
    }
    table = pa.table({field.name: columns[field.name] for field in schema}, schema=schema)
    return table, skipped


def export_particle_ids(
    folder: Path | str,
    catalog: Path | str,
    name: str | None = None,
    reference: str = REFERENCE_RADIUS,
    snapshots: Sequence[str] | None = None,
) -> Path:
    """Write the shell files into an existing catalog directory.

    Only the final snapshot is done by default: the regions are large and the
    zoom-in machinery only ever asks about z=0 halos.  ``snapshots`` overrides
    that with an explicit list of snapshot extensions.
    """
    catalog = Path(catalog)
    simulation = find_simulation(folder, name)
    schema = io.schema_for(SHELL_COLUMNS)

    wanted = list(snapshots) if snapshots else [simulation.final.extension]
    by_extension = {snapshot.extension: snapshot for snapshot in simulation}
    written: dict[str, dict[str, int]] = {}

    for extension in wanted:
        if extension not in by_extension:
            raise KeyError(f"{extension} is not a snapshot of {simulation.name}")
        snapshot = by_extension[extension]

        halos = io.read_catalog(
            catalog / f"halos_{extension}.parquet",
            columns=["halo_id", reference, "shrink_center_x", "shrink_center_y", "shrink_center_z"],
        ).table
        logger.info("%s: %d halos", extension, halos.num_rows)

        sim = snapshot.load()
        table, skipped = build_shell_table(snapshot, halos, sim, reference, schema)
        io.write_snapshot(
            table,
            catalog,
            snapshot=extension,
            prefix=io.SHELL_PREFIX,
            column_encoding=ID_COLUMN_ENCODING,
            pugs_simulation=simulation.name,
            pugs_shell_reference=reference,
        )
        written[extension] = {
            "halos": halos.num_rows - skipped,
            "skipped": skipped,
            "rows": table.num_rows,
            "particle_ids": int(np.sum(table.column("n_particles").to_numpy())),
        }
        logger.info(
            "%s: %d shell rows, %d ids%s",
            extension,
            table.num_rows,
            written[extension]["particle_ids"],
            f" ({skipped} halos skipped)" if skipped else "",
        )

    _record_provenance(catalog, reference, written)
    return catalog


def _record_provenance(catalog: Path, reference: str, written: dict[str, Any]) -> None:
    """Add the shell settings to the catalog's provenance sidecar."""
    provenance = io.read_provenance(catalog)
    provenance["particle_ids"] = {
        "reference_radius": reference,
        "shell_step_rvir": SHELL_STEP,
        "shell_max_rvir": SHELL_MAX,
        "n_shells": N_SHELLS,
        "edges_rvir": [float(edge) for edge in SHELL_EDGES],
        "encoding": "DELTA_BINARY_PACKED + zstd",
        "snapshots": written,
    }
    io.write_provenance(catalog, provenance)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("folder", help="directory holding the snapshots and AHF output")
    parser.add_argument("catalog", help="existing catalog directory to write the shells into")
    parser.add_argument("-n", "--name", default=None, help="simulation name")
    parser.add_argument(
        "--reference",
        default=REFERENCE_RADIUS,
        choices=("R200", "R500", "max_radius"),
        help=f"catalog column used as R_vir (default: {REFERENCE_RADIUS})",
    )
    parser.add_argument(
        "--snapshots",
        nargs="*",
        default=None,
        help="snapshot extensions to process (default: the final snapshot only)",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")
    export_particle_ids(
        args.folder,
        args.catalog,
        name=args.name,
        reference=args.reference,
        snapshots=args.snapshots,
    )
    logger.info("wrote particle-id shells into %s", args.catalog)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
