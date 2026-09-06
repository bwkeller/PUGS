---
title: API Reference
---

# API Reference

The `pugs` Python package builds and reads the Parquet halo catalog. It has no
database layer: every module works from the AHF output files and the snapshots
directly.

---

## Building a catalog

### `pugs.export`

The whole pipeline. Discovers snapshots, reads AHF's halo tables, builds the
merger forest, measures the particle-derived properties, and writes one Parquet
file per snapshot.

See [the catalog reference](catalog.md).

### `pugs.ahf`

Readers for AHF's output: `AHF_halos`, `AHF_particles`, `AHF_substructure` and
the `AHF_croco` merger tree.

See [pugs.ahf](ahf.md).

### `pugs.simulation`

Snapshot discovery and ordering, and the deterministic `halo_id` scheme.

### `pugs.halo_properties`

Per-halo physics measured from the snapshot with pynbody: shrinking-sphere
centre, `max_radius`, and the spherical-overdensity radii and masses.

### `pugs.merger_forest`

The merger forest, built from AHF's own tree files. Major-merger counts, the
redshift of the last major merger, and mass-assembly redshifts.

---

## Reading a catalog

### `pugs.io`

The format layer: the schema, the units, and the provenance. `read_catalog()`
returns a `Catalog` exposing `.table`, `.unit()`, `.system()`, `.description()`
and `.provenance`.

See [the catalog reference](catalog.md).

---

## Zoom-in initial conditions

### `pugs.genetic`

Helpers for turning a halo in the catalog into GenetIC input.

See [pugs.genetic](genetic.md).

---

## Quick example

```{toctree}
:hidden:

catalog
ahf
genetic
```

```python
from pugs.io import read_catalog

catalog = read_catalog("NUGS128_catalog")

print(len(catalog))                    # halos across all snapshots
print(catalog.unit("M200"))            # 'Msol'
print(catalog.system("Mhalo"))         # 'finder' -- comoving, h-scaled

# Read only the columns you need
masses = read_catalog("NUGS128_catalog", columns=["halo_id", "z", "M200", "N_mm"])
```

Or with no PUGS import at all:

```sql
SELECT snapshot, count(*), max(M200)
FROM 'NUGS128_catalog/halos_*.parquet'
WHERE M200 > 1e12
GROUP BY snapshot
ORDER BY snapshot;
```

Generating a zoom IC for a halo:

```python
from pugs.genetic import build_param_file, write_particle_ids

n = write_particle_ids("/path/to/NUGS2048", halo_id, "id_file.txt", radius_factor=8)
build_param_file(filename="genetIC_zoom.txt", outname="halo1", base_grid=2048)
```
