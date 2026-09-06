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
from urllib.parse import unquote, urlparse

import pyarrow as pa
import pyarrow.parquet as pq

#: Files making up a catalog; chosen so the JSON sidecar is not swept up by a glob.
CATALOG_PREFIX = "halos"
CATALOG_GLOB = f"{CATALOG_PREFIX}_*.parquet"

#: Particle-id shell files, written alongside the halo tables.
SHELL_PREFIX = "shells"
SHELL_GLOB = f"{SHELL_PREFIX}_*.parquet"

PROVENANCE_FILENAME = "provenance.json"

#: zstd beats gzip on both ratio and speed for this data, and every reader we
#: care about (pyarrow, duckdb, polars, Spark) supports it.
COMPRESSION = "zstd"
COMPRESSION_LEVEL = 9

#: Dependencies whose source revision is recorded alongside their version.
_REVISION_PACKAGES = ("pynbody", "pyarrow", "numpy")

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
    "list<int64>": pa.list_(pa.int64()),
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


def write_snapshot(
    table: pa.Table,
    directory: Path | str,
    snapshot: str,
    prefix: str = CATALOG_PREFIX,
    column_encoding: dict[str, str] | None = None,
    **footer: str,
) -> Path:
    """Write one snapshot's table to ``<directory>/<prefix>_<snapshot>.parquet``.

    ``column_encoding`` pins the Parquet encoding of particular columns, given
    as ``{"column.list.element": "DELTA_BINARY_PACKED"}`` for a list column's
    values.  Dictionary encoding is switched off whenever it is used: pyarrow
    will not apply both, and a dictionary would replace the values with
    dictionary indices, leaving the delta encoder nothing to compress.  The
    cost is negligible here -- the only column that dictionary-encodes well is
    the constant ``snapshot`` string, which zstd handles anyway.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{prefix}_{snapshot}.parquet"
    pq.write_table(
        apply_schema_metadata(table, pugs_snapshot=snapshot, **footer),
        path,
        compression=COMPRESSION,
        compression_level=COMPRESSION_LEVEL,
        use_dictionary=column_encoding is None,
        column_encoding=column_encoding,
    )
    return path


def build_provenance(**extra: Any) -> dict[str, Any]:
    """Record of what produced a catalog: code versions, revisions, timestamp.

    ``versions`` alone cannot identify a build -- a fork carrying a patch
    reports the same version as the release it branched from -- so
    ``revisions`` records the source commit of each dependency where one can be
    determined.
    """
    provenance: dict[str, Any] = {
        "pugs_schema_version": schema_version(),
        "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_sha": _git_sha(),
        "versions": _versions(),
        "revisions": _revisions(),
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


def _git_sha(path: Path | None = None) -> str | None:
    """HEAD of the git checkout containing ``path``, or None if there isn't one.

    Defaults to the PUGS package itself.  Returns None for an installed wheel,
    which is not in a checkout.
    """
    target = Path(path) if path is not None else Path(__file__).resolve().parent
    return _git(target, "rev-parse", "HEAD")


def _git_is_dirty(path: Path) -> bool | None:
    """Whether the checkout containing ``path`` has uncommitted *tracked* changes.

    Untracked files are ignored: a checkout used for development accumulates
    scratch output that has no bearing on the installed code, and counting it
    would leave the flag permanently true and therefore useless.
    """
    status = _git(Path(path), "status", "--porcelain", "--untracked-files=no")
    return None if status is None else bool(status)


#: Directories that hold installed copies rather than checkouts.
_INSTALL_DIRS = ("site-packages", "dist-packages")


def _is_installed_copy(path: Path) -> bool:
    """Whether ``path`` is an installed copy rather than a working tree.

    This matters because a virtualenv often sits *inside* a source checkout.
    Walking up from ``.venv/lib/python3.12/site-packages/numpy`` finds the
    enclosing project's git repository and would otherwise attribute that
    project's commit to numpy.
    """
    return any(part in _INSTALL_DIRS for part in Path(path).parts)


def _git(path: Path, *args: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(path), *args],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip()


def _revisions() -> dict[str, dict[str, Any]]:
    """Source revision of each dependency, where one can be determined.

    A version number is not enough to identify a build: a fork carrying a patch
    reports the same version as the release it was branched from.  This records
    the commit as well, so a catalog can be traced to the exact code that
    produced it.
    """
    found = {}
    for name in _REVISION_PACKAGES:
        revision = _package_revision(name)
        if revision:
            found[name] = revision
    return found


def _package_revision(name: str) -> dict[str, Any] | None:
    """Commit and source URL of an installed distribution, if discoverable.

    Three install shapes, in order of reliability:

    1. Installed from a VCS URL -- pip records the exact commit in
       ``direct_url.json`` (PEP 610).  Authoritative.
    2. Installed from a local directory -- the URL is recorded but no commit,
       so the checkout is queried.  Best effort: it reports where that checkout
       stands *now*, which may have moved since the install, hence
       ``resolved_from``.
    3. Installed in place (editable) -- the package directory is itself inside
       a checkout.  Installed *copies* are excluded from this last case: a
       virtualenv commonly sits inside a source checkout, so walking up from
       ``site-packages`` would attribute the enclosing project's commit to
       every dependency.
    """
    try:
        distribution = importlib_metadata.distribution(name)
    except importlib_metadata.PackageNotFoundError:
        return None

    revision: dict[str, Any] = {}
    raw = None
    try:
        raw = distribution.read_text("direct_url.json")
    except Exception:  # nosec - metadata may be absent or unreadable
        raw = None

    if raw:
        try:
            direct = json.loads(raw)
        except json.JSONDecodeError:
            direct = {}
        if direct.get("url"):
            revision["url"] = direct["url"]
        commit = (direct.get("vcs_info") or {}).get("commit_id")
        if commit:
            revision["commit"] = commit
            revision["resolved_from"] = "install record"
        elif str(direct.get("url", "")).startswith("file://"):
            source = Path(unquote(urlparse(direct["url"]).path))
            sha = _git_sha(source)
            if sha:
                revision.update(
                    commit=sha, resolved_from="working tree", dirty=_git_is_dirty(source)
                )

    if "commit" not in revision:
        installed = Path(distribution.locate_file(name))
        if not _is_installed_copy(installed):
            sha = _git_sha(installed)
            if sha:
                revision.update(
                    commit=sha, resolved_from="working tree", dirty=_git_is_dirty(installed)
                )
    return revision or None


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


def catalog_files(path: Path | str, glob: str = CATALOG_GLOB) -> list[Path]:
    """The Parquet files making up a catalog, sorted by snapshot."""
    path = Path(path)
    if path.is_file():
        return [path]
    return sorted(path.glob(glob))


def read_catalog(
    path: Path | str,
    columns: Sequence[str] | None = None,
    *,
    check_version: bool = True,
    glob: str = CATALOG_GLOB,
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
    files = catalog_files(path, glob)
    if not files:
        raise FileNotFoundError(f"no {glob} files under {path}")

    if check_version:
        for file in files:
            _check_version(file)

    table = pq.read_table(files, columns=list(columns) if columns else None)
    return Catalog(table=table, provenance=read_provenance(path), path=path)


def read_shells(
    path: Path | str,
    columns: Sequence[str] | None = None,
    *,
    check_version: bool = True,
) -> Catalog:
    """Read the particle-id shell files of a catalog.

    The ``particle_ids`` column is by far the largest thing in a PUGS catalog,
    so pass ``columns`` to inspect the shell geometry without reading it --
    ``["halo_id", "shell", "n_particles"]`` costs almost nothing.
    """
    return read_catalog(path, columns, check_version=check_version, glob=SHELL_GLOB)


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
