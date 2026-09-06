"""
Discovering the snapshots of a simulation and their AHF catalogs.

This replaces what ``tangos add`` used to do: walk a simulation directory, pair
each snapshot with its AHF output, and put them in time order.  Nothing is
written to a database -- a :class:`Simulation` is just an ordered list of
:class:`Snapshot` records built by looking at the filesystem.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

#: e.g. "DM128.08192.z0.000.AHF_halos" -> snapshot "DM128.08192", z 0.000
_AHF_HALOS = re.compile(r"^(?P<snapshot>.+)\.z(?P<redshift>-?\d+\.\d+)\.AHF_halos$")

#: Halo ids are ``snapshot_index * HALO_ID_STRIDE + finder_id``.  A stride of
#: 2**32 keeps the two parts separable by eye in hex and leaves room for far
#: more halos per snapshot than any catalog will hold, while staying inside
#: int64 for any plausible number of snapshots.
HALO_ID_STRIDE = 1 << 32


def halo_id(snapshot_index: int, finder_id):
    """Globally unique halo id from a snapshot index and an AHF halo id.

    Deterministic, so a catalog can be rebuilt and still refer to the same
    halos -- unlike a database row id, which depends on insertion order.
    """
    return snapshot_index * HALO_ID_STRIDE + finder_id


def split_halo_id(value):
    """Inverse of :func:`halo_id`: ``(snapshot_index, finder_id)``."""
    return value // HALO_ID_STRIDE, value % HALO_ID_STRIDE


@dataclass(frozen=True)
class Snapshot:
    """One snapshot and the AHF catalog that goes with it."""

    index: int
    extension: str
    path: Path
    ahf_stem: Path
    catalog_redshift: float

    def ahf(self, suffix: str) -> Path:
        """Path to an AHF file, e.g. ``snapshot.ahf("halos")``."""
        return self.ahf_stem.with_name(f"{self.ahf_stem.name}.AHF_{suffix}")

    def has_ahf(self, suffix: str) -> bool:
        return self.ahf(suffix).is_file()

    @cached_property
    def _header(self) -> dict:
        """Cosmology from the snapshot header, read without loading arrays."""
        import pynbody

        return dict(pynbody.load(str(self.path)).properties)

    @property
    def scalefactor(self) -> float:
        return float(self._header["a"])

    @property
    def redshift(self) -> float:
        """Redshift from the snapshot header.

        Preferred over the value in the AHF filename, which is rounded to three
        decimal places.
        """
        return 1.0 / self.scalefactor - 1.0

    def load(self, physical_units: bool = True):
        """Load the snapshot with pynbody."""
        import pynbody

        sim = pynbody.load(str(self.path))
        if physical_units:
            sim.physical_units()
        return sim


@dataclass(frozen=True)
class Simulation:
    """An ordered set of snapshots sharing a directory."""

    name: str
    path: Path
    snapshots: list[Snapshot]

    def __len__(self) -> int:
        return len(self.snapshots)

    def __iter__(self):
        return iter(self.snapshots)

    def __getitem__(self, index):
        return self.snapshots[index]

    @property
    def final(self) -> Snapshot:
        return self.snapshots[-1]


def find_simulation(folder: Path | str, name: str | None = None) -> Simulation:
    """Build a :class:`Simulation` from a directory of snapshots and AHF output.

    Snapshots are ordered by redshift, latest last, so ``simulation.final`` is
    the z=0 snapshot.  A snapshot is included only when both the snapshot file
    and its ``AHF_halos`` file are present.
    """
    folder = Path(folder)
    if not folder.is_dir():
        raise NotADirectoryError(folder)

    found: list[tuple[float, str, Path, Path]] = []
    for halos_file in sorted(folder.glob("*.AHF_halos")):
        match = _AHF_HALOS.match(halos_file.name)
        if match is None:
            continue
        extension = match.group("snapshot")
        snapshot_path = folder / extension
        if not snapshot_path.is_file():
            continue
        stem = halos_file.with_name(halos_file.name[: -len(".AHF_halos")])
        found.append((float(match.group("redshift")), extension, snapshot_path, stem))

    if not found:
        raise FileNotFoundError(f"no snapshot/AHF_halos pairs found in {folder}")

    found.sort(key=lambda entry: -entry[0])  # decreasing redshift == increasing time
    snapshots = [
        Snapshot(
            index=index,
            extension=extension,
            path=path,
            ahf_stem=stem,
            catalog_redshift=redshift,
        )
        for index, (redshift, extension, path, stem) in enumerate(found)
    ]
    return Simulation(name=name or folder.name, path=folder, snapshots=snapshots)
