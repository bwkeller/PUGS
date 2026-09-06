"""
Helpers for generating zoom-in initial conditions for individual halos.

There are two ways to get a zoom region:

:func:`particle_ids_from_catalog`
    Reads the precomputed shells written by :mod:`pugs.particle_ids`.  Needs
    only the catalog, not the snapshots, so a zoom IC can be set up from the
    distributed catalog alone.  Limited to the shell grid -- radii in steps of
    ``SHELL_STEP`` out to ``SHELL_MAX`` R_vir.

:func:`particle_ids`
    Selects from the snapshot directly.  Slower, and needs the simulation on
    disk, but works at any radius and for halos at any snapshot.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import ahf, halo_properties, io
from .particle_ids import REFERENCE_RADIUS, SHELL_MAX, SHELL_STEP
from .simulation import Simulation, find_simulation, split_halo_id


def resolve_halo(simulation: Simulation, halo_id: int):
    """Find the snapshot and AHF id a catalog ``halo_id`` refers to."""
    snapshot_index, finder_id = split_halo_id(int(halo_id))
    if not 0 <= snapshot_index < len(simulation):
        raise ValueError(
            f"halo_id {halo_id} names snapshot {snapshot_index}, "
            f"but {simulation.name} has {len(simulation)} snapshots"
        )
    return simulation[snapshot_index], int(finder_id)


def particle_ids(
    folder: Path | str,
    halo_id: int,
    radius_factor: float | None = None,
    reference: str = REFERENCE_RADIUS,
) -> np.ndarray:
    """Snapshot indices of the particles making up a halo's zoom region.

    With no ``radius_factor`` these are exactly the particles AHF assigned to
    the halo.  Given one, the region is every particle within
    ``radius_factor`` times the halo's ``reference`` radius of its centre.

    ``reference`` defaults to the same radius the stored shells use, so this
    and :func:`particle_ids_from_catalog` return the same region for the same
    factor.  ``max_radius`` is the other defensible choice -- it is what the
    older TANGOS-based pipeline used -- and differs by a few parts in a
    thousand.

    The returned indices are sorted by radius, so a smaller region is a prefix
    of a larger one.
    """
    import pynbody

    simulation = find_simulation(folder)
    snapshot, finder_id = resolve_halo(simulation, halo_id)
    membership = ahf.read_particle_membership(snapshot.ahf("particles"))
    members = membership.particles(finder_id)

    if radius_factor is None:
        return members

    sim = snapshot.load()
    measurement = halo_properties.measure_halo(sim, members, halo_properties.critical_density(sim))
    centre = measurement.shrink_center
    radius = (
        measurement.max_radius if reference == "max_radius" else measurement.overdensity[reference]
    )

    region = sim[pynbody.filt.Sphere(radius_factor * radius, centre)]
    indices = region.get_index_list(sim.ancestor)

    with halo_properties.recentred(region, centre):
        radii = np.sqrt((np.asarray(region["pos"], dtype=np.float64) ** 2).sum(axis=1))
    return np.asarray(indices)[np.argsort(radii)]


def particle_ids_from_catalog(
    catalog: Path | str, halo_id: int, radius_factor: float = SHELL_MAX
) -> np.ndarray:
    """Zoom region for a halo, read from the catalog's stored shells.

    No snapshot is needed: the region is assembled from the shells written by
    :mod:`pugs.particle_ids`.  ``radius_factor`` must land on the shell grid,
    i.e. be a multiple of ``SHELL_STEP`` and no larger than ``SHELL_MAX``.

    Ids come back innermost shell first, each shell in ascending order.  They
    are not globally sorted; GenetIC treats the file as a set, and sorting
    would discard the radial grouping for no benefit.
    """
    shells_needed = radius_factor / SHELL_STEP
    if (
        radius_factor <= 0
        or radius_factor > SHELL_MAX
        or abs(shells_needed - round(shells_needed)) > 1e-9
    ):
        raise ValueError(
            f"radius_factor must be a positive multiple of {SHELL_STEP} up to "
            f"{SHELL_MAX}; got {radius_factor}. Use particle_ids() to select "
            "an arbitrary radius from the snapshot instead."
        )

    table = io.read_shells(catalog, columns=["halo_id", "shell", "particle_ids"]).table
    wanted = int(round(shells_needed))

    ids, found = [], False
    rows = zip(
        table.column("halo_id").to_pylist(),
        table.column("shell").to_pylist(),
        table.column("particle_ids").to_pylist(),
    )
    for row_halo, shell, values in sorted(rows, key=lambda row: row[1]):
        if row_halo != int(halo_id):
            continue
        found = True
        if shell < wanted:
            ids.append(np.asarray(values, dtype=np.int64))
    if not found:
        raise KeyError(f"halo {halo_id} has no stored shells in {catalog}")
    return np.concatenate(ids) if ids else np.empty(0, dtype=np.int64)


def write_particle_ids(
    folder: Path | str,
    halo_id: int,
    filename: str = "id_file.txt",
    radius_factor: float | None = None,
) -> int:
    """Write a halo's particle ids to the plain-text file GenetIC expects.

    Selects from the snapshot; see :func:`write_particle_ids_from_catalog` for
    the version that needs only a catalog.

    :param folder: directory holding the snapshots and AHF output
    :param halo_id: catalog ``halo_id`` of the halo to zoom on
    :param filename: output path, one integer per line
    :param radius_factor: expand the region to this multiple of ``max_radius``
    :returns: the number of ids written
    """
    ids = particle_ids(folder, halo_id, radius_factor=radius_factor)
    np.savetxt(filename, ids, fmt="%d")
    return len(ids)


def write_particle_ids_from_catalog(
    catalog: Path | str,
    halo_id: int,
    filename: str = "id_file.txt",
    radius_factor: float = SHELL_MAX,
) -> int:
    """Write a zoom region from the catalog's stored shells.

    :returns: the number of ids written
    """
    ids = particle_ids_from_catalog(catalog, halo_id, radius_factor=radius_factor)
    np.savetxt(filename, ids, fmt="%d")
    return len(ids)


def build_param_file(
    filename="genetIC_zoom.txt",
    outname="PUGS_zoom",
    base_grid=2048,
    zoom_grid=[10, 2048],
    autopad=1,
    subsample=8,
    template="../inputs/zoom_template.txt",
):
    """Fill the GenetIC zoom template with per-halo parameters."""
    with open(template, "r") as handle:
        contents = handle.read()
    with open(filename, "w") as out:
        out.write(
            contents.format(
                outname=outname,
                base_grid=base_grid,
                zoom_grid=zoom_grid,
                autopad=autopad,
                subsample=subsample,
            )
        )
