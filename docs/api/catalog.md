---
title: pugs.io / pugs.export
---

# The Parquet halo catalog

`pugs.export` builds the catalog directly from the AHF output and the
snapshots; `pugs.io` reads it back with its units and provenance attached.
There is no database anywhere in the pipeline.

## Layout on disk

```
NUGS128_catalog/
    halos_DM128.00512.parquet
    halos_DM128.01024.parquet
    ...
    halos_DM128.08192.parquet
    provenance.json
```

One file per snapshot, all sharing an identical schema, so the whole directory
unions without `union_by_name`. Snapshots left empty by the particle cut still
get a zero-row file.

## Building one

```bash
cd builder/catalog
source test_vars          # or prod_vars
./build_catalog.sh test_vars
```

or from Python:

```python
from pugs.export import export_simulation

export_simulation("/path/to/NUGS128", "NUGS128_catalog", name="NUGS128")
```

Pass `measure=False` (or `--no-measure`) to skip the particle-derived columns
and leave them null; the run drops from minutes to seconds.

## Reading one

```python
from pugs.io import read_catalog

catalog = read_catalog("NUGS128_catalog")
print(len(catalog))                      # halos across all snapshots
print(catalog.unit("M200"))              # 'Msol'
print(catalog.system("M200"))            # 'physical'
print(catalog.description("z_lmm"))      # 'Redshift of the last major merger; ...'
print(catalog.provenance["versions"])    # code versions that produced it
```

Pass `columns` to read only what you need — the whole point of a columnar
format is that three columns out of fifty cost three columns of I/O:

```python
masses = read_catalog("NUGS128_catalog", columns=["halo_id", "z", "M200"])
```

### With DuckDB

The catalog is plain Parquet, so no PUGS install is required to query it:

```sql
SELECT snapshot, count(*), max(M200)
FROM 'NUGS128_catalog/halos_*.parquet'
WHERE M200 > 1e10
GROUP BY snapshot
ORDER BY snapshot;
```

## Two unit systems

This is the single most likely source of a wrong answer, so every column
declares which system it is in, queryable via `catalog.system(name)`:

| `system` | Meaning |
|---|---|
| `physical` | Measured by PUGS from the snapshot with pynbody. Physical units, `h` divided out — `Msol`, `kpc`. |
| `finder` | Read verbatim from AHF. Comoving lengths, `h`-scaled masses and lengths — `Msol / littleh`, `kpc / littleh`. |
| `pugs` | Bookkeeping and derived columns: `halo_id`, `snapshot`, `z`, `a`, the tree pointers, and the assembly history. |

`M200` and `Mhalo` are both halo masses, but they are **not** directly
comparable: `M200` is physical `Msol`, `Mhalo` is `Msol / littleh`.

## Where metadata lives

Three places, each for a reason:

| Location | Holds | Why there |
|---|---|---|
| `pugs/columns.toml` | dtype, unit, unit system, description | Hand-maintained, versioned in git, shipped with the package. TOML because YAML 1.1 parses `1e10` and `2.775e11` as *strings* — unusable for a catalog of masses and densities. |
| `provenance.json` | code versions, git SHA, cosmology, row counts, timestamp | Machine-written per export; never hand-edited, so JSON's lack of comments costs nothing. |
| Parquet footer | units mirrored as JSON, plus `pugs_schema_version` and `pugs_git_sha` | Travels *inside* each file, so a file separated from its directory is still self-identifying. |

Arrow field metadata carries the units for pyarrow readers, but it is stashed
inside an opaque base64 `ARROW:schema` footer entry that DuckDB and Spark
cannot read. The units are therefore mirrored into a plain file-level JSON
entry that any reader can reach:

```sql
SELECT key, value FROM parquet_kv_metadata('NUGS128_catalog/halos_DM128.08192.parquet');
```

### Schema versioning

`columns.toml` carries a `schema_version`, stamped into every file's footer.
`read_catalog` compares the two and raises `SchemaVersionError` on a mismatch,
which converts a stale sidecar or an out-of-date PUGS install from *silently
wrong units* into an exception. Bump it whenever a column is renamed, removed,
or has its unit changed; adding a column is backwards compatible and needs no
bump.

Pass `check_version=False` to read a catalog anyway.

## Adding a column

1. Add a `[columns.<name>]` block to `pugs/columns.toml` with `dtype`, `unit`,
   `system` and `description`.
2. Add the name to `AHF_COLUMNS`, `MEASURED_COLUMNS`, `TREE_COLUMNS` or
   `DERIVED_COLUMNS` in `pugs/export.py`.
3. A vector-valued measurement is written as one scalar column per component,
   as `shrink_center_x/y/z`; document each separately.

The writer raises `UndocumentedColumnError` if you skip step 1, so the
dictionary cannot fall behind the data.

## Identity

`halo_id` is `snapshot_index * 2**32 + finder_id`, so it is deterministic: a
rebuilt catalog refers to the same halos, and `pugs.simulation.split_halo_id`
recovers the two parts by arithmetic rather than a lookup.

`finder_id` is AHF's own halo id within its snapshot. (The TANGOS-derived
catalogs had a `halo_number` that was offset from it by one.)

## Merger tree columns

`main_progenitor_id` and `descendant_id` are `halo_id` values, so the tree can
be walked inside the catalog with a self-join and no separate links table:

```sql
SELECT h.halo_id, h.M200, p.M200 AS progenitor_M200
FROM 'cat/halos_*.parquet' h
JOIN 'cat/halos_*.parquet' p ON p.halo_id = h.main_progenitor_id;
```

A `descendant_id` may name a halo below the particle cut, which is therefore
absent from the catalog. The forest is built from every halo AHF found, so the
pointer is still correct; it simply has no row.

`N_mm`, `z_lmm` and `z25/z50/z75_mass` are written for **every** snapshot, not
only z=0: the forest is built once for the whole simulation, so they cost
nothing extra.

## Nulls

Halos missing a property get a null, not a sentinel. The one deliberate
sentinel is `z_lmm = -1`, meaning "this halo has never had a major merger",
which is distinct from "unknown".
