---
title: Halo Catalog
---

# Stage 2 — The halo catalog

**Directory:** `builder/catalog/`

Takes the N-body snapshots produced by running the volume IC through your
N-body code, together with the AHF halo-finder output beside them, and builds
a Parquet halo catalog. There is no intermediate database: every column is
derived directly from the AHF files and the snapshots.

```
snapshots + AHF output
        │
        ├── AHF_halos       ──►  43 finder columns, verbatim
        ├── AHF_particles   ──►  halo membership (snapshot indices)
        ├── AHF_croco       ──►  merger tree, with real shared-particle merits
        │                            │
        │                            ▼
        │                    pugs.merger_forest
        │                            │
        │                            ▼
        │              main progenitor / descendant pointers,
        │              N_mm, z_lmm, z25/z50/z75_mass
        │
        └── snapshot particles ──►  pugs.halo_properties
                                        │
                                        ▼
                          shrink_center, max_radius,
                          R200, R500, M200, M500
                                        │
                                        ▼
                            halos_<snapshot>.parquet  (one per snapshot)
                                 + provenance.json
```

## Running it

```bash
cd builder/catalog
source test_vars          # or prod_vars
./build_catalog.sh test_vars
```

or directly:

```bash
python -m pugs.export /path/to/NUGS128 --output NUGS128_catalog
```

| Variable | Meaning |
|---|---|
| `PUGS_SIMULATION` | Directory holding the snapshots and AHF output |
| `PUGS_NAME` | Simulation name recorded in the catalog |
| `PUGS_CATALOG` | Output catalog directory |
| `PUGS_MIN_PARTICLES` | Minimum particles for a halo to be kept (default 500) |

Useful flags:

- `--no-measure` skips the particle-derived columns, leaving them null. The run
  drops from minutes to seconds, which is what you want when only the finder
  columns or the merger trees matter.
- `--min-particles N` changes the cut. It applies to the *catalog*, not to the
  merger forest, which always uses every halo AHF found.

## What comes from where

**`AHF_halos`** is a plain tab-separated table with a `#` header. All 43
non-gas, non-star columns are copied through unchanged, keeping AHF's own
comoving, h-scaled units. This includes `hostHalo` and `numSubStruct`, which
describe the substructure hierarchy.

**`AHF_particles`** lists the particles of each halo. In these catalogs the
values are *indices into the snapshot*, not `iord` values — the tipsy files
carry no `iord` array. This is what makes the particle-derived measurements
possible without any halo-finder re-run.

**`AHF_croco`** is the merger tree, carrying `SharedPart` and AHF's merit
function `shared² / (n_descendant × n_progenitor)` for every link. PUGS reads
this rather than `AHF_mtree`, which holds the same edges with those numbers
stripped out — leaving only an ordering, from which a weight can only be
invented.

## The particle-derived columns

For every halo above the cut, PUGS loads the snapshot particles AHF assigned to
it and measures:

| Column | How |
|---|---|
| `finder_mass` | Sum of the member particle masses |
| `shrink_center` | Shrinking-sphere centre, shrink factor 0.8 |
| `max_radius` | Distance to the outermost member particle |
| `R200`, `R500` | Bisection on enclosed density over *all* particles in the region |
| `M200`, `M500` | `contrast × ρ_crit × 4/3 π R³` |

`max_radius` is **not** an overdensity radius. It measures the extent of AHF's
particle set and is close to, but not equal to, `R200`; it is the radius that
zoom regions are defined in multiples of.

### Periodicity

Halos straddling a box boundary are handled in two different ways, on purpose.
Centring wraps explicitly, because the shrinking-sphere algorithm knows nothing
about the box: without the wrap, such a halo's particles span the full box
width and the centre lands in the empty space between the two pieces, after
which the overdensity bisection fails outright. Region selection instead relies
on `pynbody.filt.Sphere` being periodic-aware, which avoids wrapping the whole
snapshot once per halo.

## The merger forest

A link is kept when **strictly more than half** of a progenitor's particles end
up in the descendant.

Thresholding on the merit instead would be wrong for merger counting, because
the merit penalises a size mismatch: a small halo swallowed whole by a much
larger one has a merit near zero while having contributed every one of its
particles. In the NUGS128 z=0 catalog, halo 0's second progenitor has merit
0.0087 and yet gave up 227 of its 227 particles — a real merger that a merit
cut would discard.

Requiring strictly more than half also makes the pruned graph a forest by
construction, since a progenitor cannot give more than half of itself to two
different descendants. That is what lets the merger counts be accumulated in a
single sweep over the snapshots.

The forest is built from **every** halo AHF found, not just those above the
catalog's particle cut. Merger counts and assembly histories are properties of
the tree, and dropping small progenitors before building it biases them: a halo
whose two progenitors both fall below the cut would look like it had no merger
at all.

## Output layout

```
NUGS128_catalog/
    halos_DM128.00512.parquet
    ...
    halos_DM128.08192.parquet
    provenance.json
```

One file per snapshot, all sharing an identical schema, so the directory unions
without `union_by_name`. Snapshots left empty by the particle cut still get a
(zero-row) file.

See [the catalog API](../api/catalog.md) for reading it back, and the
[metadata contract](../api/catalog.md#where-metadata-lives).
