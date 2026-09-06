"""
Readers for AHF halo-finder output.

AHF writes one set of plain-text files per snapshot, all sharing the stem
``<snapshot>.z<redshift>.AHF_``:

``AHF_halos``
    One row per halo, 83 whitespace-separated columns, names in a ``#`` header.
    The separator is not consistent between AHF runs -- the NUGS128 catalogs
    are tab-separated with space padding, the NUGS2048 ones single-space -- so
    it is detected from the header rather than assumed.
``AHF_particles``
    Halo membership: a ``<npart> <halo_id>`` line followed by that many
    ``<particle> <type>`` lines.  The particle values are *indices into the
    snapshot*, not ``iord`` values -- the tipsy files carry no ``iord`` array.
``AHF_substructure``
    A ``<halo_id> <n_sub>`` line followed by a line listing the subhalo ids.
``AHF_fpos``
    A byte offset into ``AHF_particles`` for each halo, in the same order as
    ``AHF_halos``.  It turns membership lookup from a whole-file parse into a
    seek, which is the difference between minutes and microseconds on the
    production catalogs.
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

import logging
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.csv as pv

logger = logging.getLogger(__name__)

#: Column header of an AHF_halos file: "#ID(1)\thostHalo(2)\t..."
_HEADER_COLUMN = re.compile(r"^(.*?)\((\d+)\)$")


def read_halo_table(path: Path | str) -> pa.Table:
    """Read an ``AHF_halos`` file into an Arrow table with AHF's column names.

    Column types are inferred by pyarrow, which gets the integer columns
    (``ID``, ``hostHalo``, ``numSubStruct``, ``npart``) right and everything
    else as double.  pyarrow's CSV reader is used rather than ``np.loadtxt``
    because the production tables have millions of rows and are gigabytes
    apiece.
    """
    path = Path(path)
    header = _header_line(path)
    delimiter = _delimiter(header)
    names = _halo_column_names(header, delimiter, path)

    table = pv.read_csv(
        path,
        read_options=pv.ReadOptions(skip_rows=1, autogenerate_column_names=True),
        parse_options=pv.ParseOptions(delimiter=delimiter, ignore_empty_lines=True),
    )
    if table.num_columns < len(names):
        raise ValueError(
            f"{path}: header names {len(names)} columns but the data parsed as "
            f"{table.num_columns} with delimiter {delimiter!r}"
        )
    # A trailing delimiter on the header line can leave an unnamed extra column.
    table = table.select(list(range(len(names))))
    return table.rename_columns(names)


def _header_line(path: Path) -> str:
    with open(path) as handle:
        header = handle.readline()
    if not header.startswith("#"):
        raise ValueError(f"{path}: expected a '#' header line, got {header[:40]!r}")
    return header


def _delimiter(header: str) -> str:
    """Field separator of an AHF table, taken from its header line.

    AHF is not consistent about this between runs: the NUGS128 catalogs are
    tab-separated with the columns space-padded to width, while the NUGS2048
    ones are single-space separated with no tabs anywhere.  Assuming either one
    makes the reader return a single column of text on the other, which then
    fails much later with a confusing missing-column error.
    """
    return "\t" if "\t" in header else " "


def _halo_column_names(header: str, delimiter: str, path: Path) -> list[str]:
    """Column names from the ``#ID(1) hostHalo(2) ...`` header line."""
    names = []
    for field in header[1:].strip().split(delimiter):
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

    def __enter__(self) -> "ParticleMembership":
        return self

    def __exit__(self, *exc) -> None:
        return None

    def close(self) -> None:
        """No-op; the whole file is already in memory."""

    def particles(self, halo_id: int) -> np.ndarray:
        """Snapshot indices belonging to one halo."""
        row = np.searchsorted(self.halo_id, halo_id)
        if row >= len(self.halo_id) or self.halo_id[row] != halo_id:
            raise KeyError(f"halo {halo_id} is not in this membership file")
        return self.indices[self.offsets[row] : self.offsets[row + 1]]


def read_particle_membership(path: Path | str) -> ParticleMembership:
    """Read a whole ``AHF_particles`` file into memory.

    This is the fallback for snapshots with no ``AHF_fpos`` index; prefer
    :func:`open_membership`, which seeks to individual halos instead.  The
    production files run to tens of gigabytes apiece, and ``np.loadtxt`` reads
    them at roughly 80 MB/s.

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


def read_fpos(path: Path | str) -> np.ndarray:
    """Read an ``AHF_fpos`` file: one byte offset per halo, in ``AHF_halos`` order.

    Each offset points at the *first particle row* of that halo's block in
    ``AHF_particles``, i.e. just past the ``<npart> <halo_id>`` header line.
    """
    offsets = np.loadtxt(path, dtype=np.int64, ndmin=1)
    if offsets.size and np.any(np.diff(offsets) <= 0):
        raise ValueError(f"{path}: offsets are not strictly increasing")
    return offsets


class IndexedMembership:
    """Random access into ``AHF_particles``, using the ``AHF_fpos`` index.

    Reads one halo's block on demand rather than parsing the whole file.  On
    the NUGS2048 z=0 catalog that file is 60 GB across ~4.6 billion rows, so
    the difference is a seek against a full parse.

    It also sidesteps a parsing problem: the block header rows and the particle
    rows do not always use the same separator -- in the NUGS128 catalogs the
    headers are space-padded while the particle rows are tab-separated, which
    no single-delimiter CSV parser can read.  Seeking past the headers leaves
    each block internally uniform.
    """

    def __init__(self, path, halo_id, offset, end, count, delimiter):
        order = np.argsort(halo_id, kind="stable")
        self.path = Path(path)
        self.halo_id = np.asarray(halo_id, dtype=np.int64)[order]
        self.offset = np.asarray(offset, dtype=np.int64)[order]
        self.end = np.asarray(end, dtype=np.int64)[order]
        self.count = np.asarray(count, dtype=np.int64)[order]
        self.delimiter = delimiter
        self._handle = None

    def __len__(self) -> int:
        return len(self.halo_id)

    def __enter__(self) -> "IndexedMembership":
        self._handle = open(self.path, "rb")
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def particles(self, halo_id: int) -> np.ndarray:
        """Snapshot indices belonging to one halo."""
        row = np.searchsorted(self.halo_id, halo_id)
        if row >= len(self.halo_id) or self.halo_id[row] != halo_id:
            raise KeyError(f"halo {halo_id} is not in {self.path}")
        count = int(self.count[row])
        if count == 0:
            return np.empty(0, dtype=np.int64)

        handle = self._handle or open(self.path, "rb")
        try:
            handle.seek(int(self.offset[row]))
            raw = handle.read(int(self.end[row]) - int(self.offset[row]))
        finally:
            if self._handle is None:
                handle.close()

        return _parse_block(raw, count, self.delimiter, self.path, halo_id)


def _parse_block(raw: bytes, count: int, delimiter: str, path: Path, halo_id: int) -> np.ndarray:
    """Parse the first ``count`` rows of a block, returning column 0.

    The byte range runs to the *next* halo's offset, so it also contains that
    halo's header line; the rows are trimmed to ``count`` before parsing.
    """
    newlines = np.flatnonzero(np.frombuffer(raw, dtype=np.uint8) == 0x0A)
    if newlines.size < count:
        raise ValueError(
            f"{path}: halo {halo_id} declares {count} particles but its block "
            f"holds {newlines.size} rows"
        )
    table = pv.read_csv(
        pa.BufferReader(raw[: newlines[count - 1] + 1]),
        read_options=pv.ReadOptions(autogenerate_column_names=True),
        parse_options=pv.ParseOptions(delimiter=delimiter),
        convert_options=pv.ConvertOptions(include_columns=["f0"]),
    )
    return np.asarray(table.column(0).to_numpy(), dtype=np.int64)


def open_membership(snapshot, halo_table):
    """Fastest available membership access for a snapshot.

    Uses the ``AHF_fpos`` byte index when it is present, which reads only the
    halos actually asked for.  Without it there is no choice but to parse the
    entire particles file.
    """
    if snapshot.has_ahf("fpos"):
        return open_indexed_membership(
            snapshot.ahf("particles"),
            snapshot.ahf("fpos"),
            halo_table.column("ID").to_numpy(),
            halo_table.column("npart").to_numpy(),
        )
    logger.warning(
        "%s has no AHF_fpos; falling back to parsing the whole particles file (%.1f GB)",
        snapshot.extension,
        snapshot.ahf("particles").stat().st_size / 1e9,
    )
    return read_particle_membership(snapshot.ahf("particles"))


def open_indexed_membership(
    particles: Path | str, fpos: Path | str, halo_id, npart
) -> IndexedMembership:
    """Build an :class:`IndexedMembership` from the AHF files of one snapshot.

    ``halo_id`` and ``npart`` come from ``AHF_halos`` and must be in file
    order, since ``AHF_fpos`` is written in that same order.
    """
    particles = Path(particles)
    offsets = read_fpos(fpos)
    halo_id = np.asarray(halo_id, dtype=np.int64)
    npart = np.asarray(npart, dtype=np.int64)
    if len(offsets) != len(halo_id):
        raise ValueError(
            f"{fpos}: {len(offsets)} offsets for {len(halo_id)} halos in the halo table"
        )

    end = np.empty_like(offsets)
    end[:-1] = offsets[1:]
    end[-1] = particles.stat().st_size
    return IndexedMembership(
        particles, halo_id, offsets, end, npart, _block_delimiter(particles, offsets)
    )


def _block_delimiter(particles: Path, offsets: np.ndarray) -> str:
    """Separator used *inside* a block, read from its first particle row."""
    if offsets.size == 0:
        return " "
    with open(particles, "rb") as handle:
        handle.seek(int(offsets[0]))
        row = handle.readline()
    return "\t" if b"\t" in row else " "


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
