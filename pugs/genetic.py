"""
Helpers for generating zoom-in initial conditions for individual halos.

Particle ids are read straight from AHF's membership file, or -- when the zoom
region needs to be larger than the halo the finder identified -- selected from
the snapshot within a multiple of the halo's radius.  The PUGS catalog no
longer stores precomputed id lists; the snapshot is already required to build
one, so computing the region on demand avoids keeping a compressed copy of it
in the catalog.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import ahf, halo_properties
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


def particle_ids(folder: Path | str, halo_id: int, radius_factor: float | None = None):
    """Snapshot indices of the particles making up a halo's zoom region.

    With no ``radius_factor`` these are exactly the particles AHF assigned to
    the halo.  Given one, the region is instead every particle within
    ``radius_factor`` times the halo's ``max_radius`` of its centre, which is
    how a zoom region larger than the halo itself is built.  The returned
    indices are sorted by radius, so a smaller region is a prefix of a larger
    one.
    """
    simulation = find_simulation(folder)
    snapshot, finder_id = resolve_halo(simulation, halo_id)
    membership = ahf.read_particle_membership(snapshot.ahf("particles"))
    members = membership.particles(finder_id)

    if radius_factor is None:
        return members

    import pynbody

    sim = snapshot.load()
    centre, max_radius = halo_properties.centre_and_extent(sim[members])
    region = sim[pynbody.filt.Sphere(radius_factor * max_radius, centre)]
    indices = region.get_index_list(sim.ancestor)

    with halo_properties.recentred(region, centre):
        radii = np.sqrt((np.asarray(region["pos"], dtype=np.float64) ** 2).sum(axis=1))
    return np.asarray(indices)[np.argsort(radii)]


def write_particle_ids(
    folder: Path | str,
    halo_id: int,
    filename: str = "id_file.txt",
    radius_factor: float | None = None,
) -> int:
    """Write a halo's particle ids to the plain-text file GenetIC expects.

    :param folder: directory holding the snapshots and AHF output
    :param halo_id: catalog ``halo_id`` of the halo to zoom on
    :param filename: output path, one integer per line
    :param radius_factor: expand the region to this multiple of ``max_radius``
    :returns: the number of ids written
    """
    ids = particle_ids(folder, halo_id, radius_factor=radius_factor)
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
