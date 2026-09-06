"""
Readers for AHF halo-finder output.

AHF writes one set of plain-text files per snapshot, all sharing the stem
``<snapshot>.z<redshift>.AHF_``:

``AHF_halos``
    One row per halo, 83 whitespace-separated columns, names in a ``#`` header.
``AHF_particles``
    Halo membership: a ``<npart> <halo_id>`` line followed by that many
    ``<particle> <type>`` lines.  The particle values are *indices into the
    snapshot*, not ``iord`` values -- the tipsy files carry no ``iord`` array.
``AHF_substructure``
    A ``<halo_id> <n_sub>`` line followed by a line listing the subhalo ids.
``AHF_croco``
    The merger tree with real shared-particle counts and AHF's merit function.
``AHF_mtree``
    The same tree with the shared-particle counts stripped out.

Halo ids are per-snapshot and, in these catalogs, run from 0.  Progenitor ids
in the tree files refer to the *previous* snapshot's numbering.

A note on the tree files: PUGS reads ``AHF_croco``, which carries
``SharedPart`` and ``merit`` for every link.  ``AHF_mtree`` holds the same
edges with those numbers removed, leaving only the progenitor ordering; code
that reads it has to invent a weight from the rank instead.  Since both files
are written side by side there is no reason to use the weaker one.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.csv as pv

#: Column header of an AHF_halos file: "#ID(1)\thostHalo(2)\t..."
_HEADER_COLUMN = re.compile(r"^(.*?)\((\d+)\)$")


def read_halo_table(path: Path | str) -> pa.Table:
    """Read an ``AHF_halos`` file into an Arrow table with AHF's column names.

    Column types are inferred by pyarrow, which gets the integer columns
    (``ID``, ``hostHalo``, ``numSubStruct``, ``npart``) right and everything
    else as double.  pyarrow's CSV reader is used rather than ``np.loadtxt``
    because the production catalogs have millions of rows.
    """
    path = Path(path)
    names = _halo_column_names(path)
    table = pv.read_csv(
        path,
        read_options=pv.ReadOptions(skip_rows=1, autogenerate_column_names=True),
        parse_options=pv.ParseOptions(delimiter="\t", ignore_empty_lines=True),
    )
    if table.num_columns < len(names):
        raise ValueError(
            f"{path}: header names {len(names)} columns but the data has {table.num_columns}"
        )
    # A trailing delimiter on the header line can leave an unnamed extra column.
    table = table.select(list(range(len(names))))
    return table.rename_columns(names)


def _halo_column_names(path: Path) -> list[str]:
    """Column names from the ``#ID(1)\\thostHalo(2)...`` header line."""
    with open(path) as handle:
        header = handle.readline()
    if not header.startswith("#"):
        raise ValueError(f"{path}: expected a '#' header line, got {header[:40]!r}")

    names = []
    for field in header[1:].strip().split("\t"):
        field = field.strip()
        if not field:
            continue
        match = _HEADER_COLUMN.match(field)
        names.append(match.group(1) if match else field)
    if not names:
        raise ValueError(f"{path}: could not parse any column names from the header")
    return names


@dataclass(frozen=True)
class ParticleMembership:
    """Halo membership in compressed-sparse-row form.

    ``indices[offsets[i]:offsets[i + 1]]`` are the snapshot indices belonging to
    the halo whose id is ``halo_id[i]``.
    """

    halo_id: np.ndarray
    offsets: np.ndarray
    indices: np.ndarray

    def __len__(self) -> int:
        return len(self.halo_id)

    def particles(self, halo_id: int) -> np.ndarray:
        """Snapshot indices belonging to one halo."""
        row = np.searchsorted(self.halo_id, halo_id)
        if row >= len(self.halo_id) or self.halo_id[row] != halo_id:
            raise KeyError(f"halo {halo_id} is not in this membership file")
        return self.indices[self.offsets[row] : self.offsets[row + 1]]


def read_particle_membership(path: Path | str) -> ParticleMembership:
    """Read an ``AHF_particles`` file.

    The file interleaves per-halo headers with particle rows, and both are two
    integers wide, so the blocks can only be separated by walking the counts.
    The bulk read is still vectorised; only the block walk is a Python loop, of
    one iteration per halo rather than per particle.
    """
    path = Path(path)
    raw = np.loadtxt(path, skiprows=1, dtype=np.int64, ndmin=2)
    if raw.size == 0:
        empty = np.empty(0, dtype=np.int64)
        return ParticleMembership(empty, np.zeros(1, dtype=np.int64), empty)

    halo_ids: list[int] = []
    starts: list[int] = []
    counts: list[int] = []
    row = 0
    while row < len(raw):
        count, halo_id = int(raw[row, 0]), int(raw[row, 1])
        halo_ids.append(halo_id)
        starts.append(row + 1)
        counts.append(count)
        row += 1 + count
    if row != len(raw):
        raise ValueError(f"{path}: block structure overruns the file at row {row}")

    keep = np.concatenate(
        [np.arange(s, s + c, dtype=np.int64) for s, c in zip(starts, counts)]
        or [np.empty(0, dtype=np.int64)]
    )
    indices = raw[keep, 0]
    offsets = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
    halo_id = np.asarray(halo_ids, dtype=np.int64)

    order = np.argsort(halo_id, kind="stable")
    if not np.array_equal(order, np.arange(len(halo_id))):
        # keep halo_id sorted so `particles` can binary-search
        regrouped = [indices[offsets[i] : offsets[i + 1]] for i in order]
        halo_id = halo_id[order]
        offsets = np.concatenate([[0], np.cumsum([len(g) for g in regrouped])]).astype(np.int64)
        indices = np.concatenate(regrouped) if regrouped else indices
    return ParticleMembership(halo_id, offsets, indices)


def read_substructure(path: Path | str) -> dict[int, np.ndarray]:
    """Read an ``AHF_substructure`` file into ``{halo_id: subhalo_ids}``.

    Only halos that have substructure appear in the file.
    """
    path = Path(path)
    result: dict[int, np.ndarray] = {}
    with open(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            header = header.strip()
            if not header or header.startswith("#"):
                continue
            halo_id, n_sub = (int(field) for field in header.split()[:2])
            children = handle.readline().split()
            if len(children) != n_sub:
                raise ValueError(
                    f"{path}: halo {halo_id} declares {n_sub} subhalos, found {len(children)}"
                )
            result[halo_id] = np.asarray(children, dtype=np.int64)
    return result


@dataclass(frozen=True)
class MergerLinks:
    """Tree edges between one snapshot and the one before it.

    Each entry links a halo in the current snapshot (``descendant``) to one in
    the previous snapshot (``progenitor``).  ``merit`` is AHF's figure of
    merit, ``shared**2 / (n_descendant * n_progenitor)``; ``rank`` is the
    position in AHF's own ordering, so ``rank == 0`` is the best progenitor.
    """

    descendant: np.ndarray
    progenitor: np.ndarray
    shared: np.ndarray
    merit: np.ndarray
    rank: np.ndarray
    n_descendant: np.ndarray
    n_progenitor: np.ndarray

    def __len__(self) -> int:
        return len(self.descendant)


def read_merger_links(path: Path | str) -> MergerLinks:
    """Read an ``AHF_croco`` merger-tree file.

    Format, after three ``#`` header lines, is a block per descendant::

        <halo_id>  <npart>  <n_progenitors>
          <shared>  <prog_id>  <prog_npart>  <merit>     x n_progenitors
    """
    path = Path(path)
    descendant, progenitor = [], []
    shared, merit, rank = [], [], []
    n_descendant, n_progenitor = [], []

    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            halo_id, halo_npart, n_prog = int(fields[0]), int(fields[1]), int(fields[2])
            for position in range(n_prog):
                entry = handle.readline().split()
                if len(entry) < 4:
                    raise ValueError(
                        f"{path}: halo {halo_id} progenitor {position} has {len(entry)} fields"
                    )
                shared.append(int(entry[0]))
                progenitor.append(int(entry[1]))
                n_progenitor.append(int(entry[2]))
                merit.append(float(entry[3]))
                descendant.append(halo_id)
                n_descendant.append(halo_npart)
                rank.append(position)

    as_int = lambda values: np.asarray(values, dtype=np.int64)  # noqa: E731
    return MergerLinks(
        descendant=as_int(descendant),
        progenitor=as_int(progenitor),
        shared=as_int(shared),
        merit=np.asarray(merit, dtype=np.float64),
        rank=as_int(rank),
        n_descendant=as_int(n_descendant),
        n_progenitor=as_int(n_progenitor),
    )
