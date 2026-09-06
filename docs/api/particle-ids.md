---
title: pugs.particle_ids
---

# `pugs.particle_ids`

Stores the particle ids around each z=0 halo, so a zoom-in region can be built
from the distributed catalog without the snapshots.

For every halo, the particles within 5 R_vir are split into concentric shells
of 0.5 R_vir and written one row per (halo, shell). The innermost shell is a
sphere rather than a shell, since its inner edge is zero.

| | |
|---|---|
| Shell width | 0.5 R_vir |
| Outermost edge | 5.0 R_vir |
| Shells per halo | 10 |
| Reference radius | `R200` by default |
| File | `shells_<snapshot>.parquet`, beside the halo tables |

---

## Why shells instead of one sorted list

The obvious layout is one list per halo sorted by radius, so that any
sub-region is a prefix. That is what the older TANGOS pipeline stored, as a
zlib-compressed blob.

But radial order scrambles the ids, and an int64 column in random order barely
compresses: there is nothing for a delta encoder to work with.

Bucketing by radius instead of sorting by it keeps each shell's ids in their
natural ascending order — the region selection already returns them sorted, so
**no sort is applied anywhere** — and consecutive values then differ by very
little. In the NUGS128 z=0 catalog the median gap is 1–3 inside R_vir and
around 100 in the outer shells, against ids spanning millions. Parquet's
`DELTA_BINARY_PACKED` encoding stores those differences instead of the values.

Measured on the 40 most massive z=0 halos of NUGS128 (534k ids):

| Layout | Size | vs raw |
|---|---|---|
| raw int64 | 4.27 MB | 1.0× |
| Parquet, zstd only | 0.68 MB | 6.3× |
| one radially sorted blob + zlib (the old scheme) | 0.88 MB | 4.8× |
| **shells, DELTA + zstd** | **0.55 MB** | **7.8×** |

The ids stay ordinary Parquet int64 values readable by anything; the
compression is the format's job, not a custom codec's.

---

## Building

```bash
pugs-particle-ids /path/to/NUGS128 NUGS128_catalog
```

or as part of the pipeline, by setting `PUGS_PARTICLE_IDS=1` before
`build_catalog.sh`. It is opt-in because the shells are much the largest part
of a catalog.

The catalog must already exist: the halo centres and `R200` are read from it
rather than recomputed.

| Flag | Meaning |
|---|---|
| `--reference` | Catalog column used as R_vir: `R200` (default), `R500`, `max_radius` |
| `--snapshots` | Snapshot extensions to process (default: the final snapshot only) |

### On the reference radius

`R200` is the halo's actual virial radius. `max_radius` — the distance to the
outermost halo-finder particle — is the other defensible choice and is what the
TANGOS pipeline used; the two agree to a few parts in a thousand. Whichever was
used is recorded in `provenance.json` under `particle_ids.reference_radius`,
and per row in the `reference_radius` column.

---

## Reading

```python
from pugs.io import read_shells

shells = read_shells("NUGS128_catalog")
print(shells.unit("r_outer"))          # 'kpc'
```

`particle_ids` is by far the largest column in a catalog, so read the geometry
without it when that is all you need:

```python
geometry = read_shells("NUGS128_catalog", columns=["halo_id", "shell", "n_particles"])
```

`n_particles` exists precisely so shell sizes can be inspected without touching
the ids.

With DuckDB, no PUGS import:

```sql
SELECT halo_id, sum(n_particles) AS n_within_2rvir
FROM 'NUGS128_catalog/shells_*.parquet'
WHERE r_outer_rvir <= 2.0
GROUP BY halo_id;
```

---

## Building a zoom region

```python
from pugs.genetic import write_particle_ids_from_catalog

n = write_particle_ids_from_catalog("NUGS128_catalog", halo_id, "id_file.txt", radius_factor=3.0)
```

The region is the union of the shells below `radius_factor`, so the factor must
land on the shell grid — a multiple of 0.5 up to 5.0. For an arbitrary radius,
or a halo at a snapshot with no stored shells, use
{func}`pugs.genetic.particle_ids`, which selects from the snapshot instead.
Both return the same particles for the same factor.

---

## Scale

Per-halo cost is one periodic sphere query at 5 R_vir plus a radius sort into
buckets. NUGS128's 219 z=0 halos take about 30 seconds and produce 1.16 M ids
in a 1.5 MB file.

Projected to NUGS2048's 461,416 z=0 halos above the 500-particle cut, which
hold 4.28 billion finder particles between them, and scaling by the 3.0 stored
ids per finder particle measured on NUGS128:

| | |
|---|---|
| stored ids | ~13 billion |
| raw int64 | ~104 GB |
| **on disk** | **~16 GB** (at 1.25 bytes/id) |

Halos overlap, so a particle in a crowded region is stored once per halo whose
5 R_vir sphere contains it. That duplication is inherent to per-halo regions,
not an artefact of this layout.
