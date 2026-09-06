"""
Reading and writing the PUGS halo catalog.

The catalog is a directory of Parquet files, one per snapshot, plus a
``provenance.json`` sidecar::

    NUGS128/
        halos_DM128.00512.parquet
        halos_DM128.01024.parquet
        ...
        provenance.json

Every file carries the full column set, so the schema is identical across the
whole directory and readers can union the files without ``union_by_name``.
Columns that only exist at the final snapshot (the assembly history) are null
everywhere else.

Metadata lives in three places, deliberately:

* ``columns.toml``, shipped in the package, is the hand-maintained data
  dictionary: dtype, unit, unit system and description for each column.
* ``provenance.json``, written next to the data, records what produced *this*
  catalog -- code versions, git SHA, cosmology, row counts.
* The Parquet footer of each file carries the units mirrored as JSON plus a
  small anchor (``pugs_schema_version``, ``pugs_git_sha``).  The anchor is what
  lets :func:`read_catalog` refuse a catalog whose sidecar has drifted, and it
  keeps a file that has been separated from its directory self-identifying.

Anything that varies per snapshot -- redshift, scale factor, snapshot name --
is a *column*, never file metadata.  pyarrow takes the footer metadata of a
multi-file dataset from a single file without warning, so per-snapshot values
stored there would be silently wrong for every other file.
"""

from __future__ import annotations

import json
import subprocess
import tomllib
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from importlib import resources
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import pyarrow as pa
import pyarrow.parquet as pq

#: Files making up a catalog; chosen so the JSON sidecar is not swept up by a glob.
CATALOG_GLOB = "halos_*.parquet"
PROVENANCE_FILENAME = "provenance.json"

#: zstd beats gzip on both ratio and speed for this data, and every reader we
#: care about (pyarrow, duckdb, polars, Spark) supports it.
COMPRESSION = "zstd"
COMPRESSION_LEVEL = 9

_DICTIONARY: dict[str, Any] | None = None


class SchemaVersionError(RuntimeError):
    """A catalog was written by an incompatible version of the data dictionary."""


class UndocumentedColumnError(KeyError):
    """A column was written that ``columns.toml`` does not describe."""


def dictionary() -> Mapping[str, Any]:
    """The parsed contents of the packaged ``columns.toml``."""
    global _DICTIONARY
    if _DICTIONARY is None:
        with resources.files("pugs").joinpath("columns.toml").open("rb") as handle:
            _DICTIONARY = tomllib.load(handle)
    return _DICTIONARY


def schema_version() -> int:
    """Schema version of the packaged data dictionary."""
    return int(dictionary()["schema_version"])


def describe(column: str) -> Mapping[str, str]:
    """Unit, unit system and description for a single column."""
    try:
        return dictionary()["columns"][column]
    except KeyError:
        raise UndocumentedColumnError(f"{column!r} is not described in pugs/columns.toml") from None


def documented_columns() -> tuple[str, ...]:
    """Every column name the data dictionary knows about."""
    return tuple(dictionary()["columns"])


#: Arrow types the dictionary is allowed to declare.
_ARROW_TYPES: dict[str, pa.DataType] = {
    "int64": pa.int64(),
    "float64": pa.float64(),
    "string": pa.string(),
    "bool": pa.bool_(),
}


def arrow_type(column: str) -> pa.DataType:
    """Declared Arrow type of a column."""
    declared = describe(column)["dtype"]
    try:
        return _ARROW_TYPES[declared]
    except KeyError:
        raise ValueError(
            f"{column!r} declares dtype {declared!r}; "
            f"columns.toml may only use {sorted(_ARROW_TYPES)}"
        ) from None


def schema_for(columns: Iterable[str]) -> pa.Schema:
    """Arrow schema for ``columns``, taken entirely from the data dictionary.

    Deriving the schema from the dictionary rather than from whatever the
    database happens to hold is what keeps every snapshot file identical: a
    property that is missing at one snapshot yields a null column, not a
    different schema.
    """
    return pa.schema([pa.field(name, arrow_type(name)) for name in columns])


# --------------------------------------------------------------------- write --


def apply_schema_metadata(table: pa.Table, **footer: str) -> pa.Table:
    """Attach the data dictionary to ``table`` as Arrow metadata.

    Per-field metadata (``unit``, ``system``, ``description``) round-trips
    through pyarrow but is invisible to non-Arrow readers, which only see an
    opaque base64 ``ARROW:schema`` footer entry.  The same information is
    therefore mirrored into a file-level ``columns`` JSON blob, which DuckDB and
    Spark *can* read.  It costs a few hundred bytes.

    Raises :class:`UndocumentedColumnError` if any column is missing from
    ``columns.toml``, so the dictionary cannot silently fall behind the writer.
    """
    entries = dictionary()["columns"]
    undocumented = [name for name in table.column_names if name not in entries]
    if undocumented:
        raise UndocumentedColumnError(
            f"columns absent from pugs/columns.toml: {sorted(undocumented)}. "
            "Add them to the dictionary before writing."
        )

    # Cast to the declared types as well as attaching metadata. pyarrow casts
    # safely by default, so a value that will not survive the declared dtype
    # raises here rather than being quietly truncated.
    fields = [
        pa.field(name, arrow_type(name)).with_metadata(
            {key: str(value) for key, value in entries[name].items()}
        )
        for name in table.column_names
    ]
    mirror = {name: entries[name] for name in table.column_names}
    file_metadata = {
        "pugs_schema_version": str(schema_version()),
        "columns": json.dumps(mirror, sort_keys=True),
        **{key: str(value) for key, value in footer.items()},
    }
    return table.cast(pa.schema(fields, metadata=file_metadata))


def write_snapshot(table: pa.Table, directory: Path | str, snapshot: str, **footer: str) -> Path:
    """Write one snapshot's halos to ``<directory>/halos_<snapshot>.parquet``."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"halos_{snapshot}.parquet"
    pq.write_table(
        apply_schema_metadata(table, pugs_snapshot=snapshot, **footer),
        path,
        compression=COMPRESSION,
        compression_level=COMPRESSION_LEVEL,
    )
    return path


def build_provenance(**extra: Any) -> dict[str, Any]:
    """Record of what produced a catalog: code versions, git SHA, timestamp."""
    provenance: dict[str, Any] = {
        "pugs_schema_version": schema_version(),
        "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_sha": _git_sha(),
        "versions": _versions(),
    }
    provenance.update(extra)
    return provenance


def write_provenance(directory: Path | str, provenance: Mapping[str, Any]) -> Path:
    """Write the ``provenance.json`` sidecar."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / PROVENANCE_FILENAME
    path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    return path


def _git_sha() -> str | None:
    """Current PUGS commit, or None outside a git checkout (e.g. an installed wheel)."""
    package_root = Path(__file__).resolve().parent
    try:
        result = subprocess.run(
            ["git", "-C", str(package_root), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout.strip() or None if result.returncode == 0 else None


def _versions() -> dict[str, str | None]:
    packages = ("PUGS", "pynbody", "pyarrow", "numpy")
    versions: dict[str, str | None] = {}
    for name in packages:
        try:
            versions[name] = importlib_metadata.version(name)
        except importlib_metadata.PackageNotFoundError:
            versions[name] = None
    return versions


# ---------------------------------------------------------------------- read --


@dataclass(frozen=True)
class Catalog:
    """A halo catalog: the table, plus the metadata needed to interpret it."""

    table: pa.Table
    provenance: dict[str, Any]
    path: Path

    def __len__(self) -> int:
        return self.table.num_rows

    @property
    def columns(self) -> list[str]:
        return self.table.column_names

    def unit(self, column: str) -> str:
        """Unit string for ``column``, taken from the file it was read from."""
        return self._field_metadata(column).get("unit", "dimensionless")

    def system(self, column: str) -> str:
        """Unit system for ``column``: 'physical', 'finder' or 'pugs'.

        'finder' columns come straight from AHF and are comoving and h-scaled;
        'physical' columns are measured by PUGS with pynbody.  Mixing them
        without converting is the most likely way to get a wrong answer.
        """
        return self._field_metadata(column).get("system", "unknown")

    def description(self, column: str) -> str:
        return self._field_metadata(column).get("description", "")

    def _field_metadata(self, column: str) -> dict[str, str]:
        field = self.table.schema.field(column)
        raw = field.metadata or {}
        return {key.decode(): value.decode() for key, value in raw.items()}

    def to_pandas(self):
        """Convert to a DataFrame.

        Note that pandas discards Arrow field metadata, so units are lost on the
        way out -- query them from the :class:`Catalog` before converting.
        """
        return self.table.to_pandas()


def catalog_files(path: Path | str) -> list[Path]:
    """The Parquet files making up a catalog, sorted by snapshot."""
    path = Path(path)
    if path.is_file():
        return [path]
    return sorted(path.glob(CATALOG_GLOB))


def read_catalog(
    path: Path | str,
    columns: Sequence[str] | None = None,
    *,
    check_version: bool = True,
) -> Catalog:
    """Read a catalog directory (or a single Parquet file).

    ``columns`` restricts the read to a subset, which is the whole point of a
    columnar format -- pulling three columns out of a fifty column catalog reads
    only those three from disk.

    With ``check_version`` (the default) the schema version stamped in each
    file's footer must match the packaged data dictionary.  This is what turns a
    stale sidecar or an out-of-date PUGS install from silently wrong units into
    an exception.
    """
    path = Path(path)
    files = catalog_files(path)
    if not files:
        raise FileNotFoundError(f"no {CATALOG_GLOB} files under {path}")

    if check_version:
        for file in files:
            _check_version(file)

    table = pq.read_table(files, columns=list(columns) if columns else None)
    return Catalog(table=table, provenance=read_provenance(path), path=path)


def read_provenance(path: Path | str) -> dict[str, Any]:
    """Read the provenance sidecar; empty dict if the catalog has none."""
    path = Path(path)
    directory = path.parent if path.is_file() else path
    sidecar = directory / PROVENANCE_FILENAME
    if not sidecar.is_file():
        return {}
    return json.loads(sidecar.read_text())


def _check_version(file: Path) -> None:
    metadata = pq.read_schema(file).metadata or {}
    raw = metadata.get(b"pugs_schema_version")
    if raw is None:
        raise SchemaVersionError(
            f"{file} carries no pugs_schema_version; it was not written by pugs.io"
        )
    found = int(raw)
    expected = schema_version()
    if found != expected:
        raise SchemaVersionError(
            f"{file} was written with schema version {found}, but this PUGS "
            f"provides version {expected}. Column units may have changed; "
            "install a matching PUGS or re-export the catalog."
        )
