---
title: pugs.genetic
---

# `pugs.genetic`

Helpers for generating zoom-in initial conditions for individual halos.

There are two ways to get a zoom region:

| | Needs | Radii |
|---|---|---|
| `particle_ids_from_catalog` | the catalog only | the shell grid: multiples of 0.5 up to 5 R_vir |
| `particle_ids` | the snapshots on disk | any radius, any snapshot |

Both return the same particles for the same factor, so the stored shells are a
true substitute for re-reading the snapshot. See
[pugs.particle_ids](particle-ids.md) for how the shells are stored.

The typical workflow is:

1. Pick a halo from the catalog and note its `halo_id`.
2. Call `write_particle_ids` to export the region to disk.
3. Call `build_param_file` to create the GenetIC parameter file.
4. Run GenetIC.

---

## `particle_ids`

```python
pugs.genetic.particle_ids(folder, halo_id, radius_factor=None) -> numpy.ndarray
```

Snapshot indices of the particles making up a halo's zoom region.

With no `radius_factor` these are exactly the particles AHF assigned to the
halo. Given one, the region is instead every particle within `radius_factor`
times the halo's `max_radius` of its centre. The returned indices are sorted by
radius, so a smaller region is a prefix of a larger one.

| Name | Type | Default | Description |
|---|---|---|---|
| `folder` | `str \| Path` | — | Directory holding the snapshots and AHF output |
| `halo_id` | `int` | — | Catalog `halo_id` of the halo |
| `radius_factor` | `float \| None` | `None` | Expand the region to this multiple of `max_radius` |

---

## `write_particle_ids`

```python
pugs.genetic.write_particle_ids(folder, halo_id, filename="id_file.txt", radius_factor=None) -> int
```

Write a halo's particle ids to the plain-text file GenetIC expects, one integer
per line. Returns the number of ids written.

**Example**

```python
from pugs.genetic import write_particle_ids
from pugs.io import read_catalog

catalog = read_catalog("NUGS2048_catalog", columns=["halo_id", "snapshot", "M200"])
# ... pick a halo_id from the table ...

n = write_particle_ids("/data/sims/NUGS2048", halo_id, "halo1_ids.txt", radius_factor=8)
print(f"wrote {n} particle ids")
```

---

## `resolve_halo`

```python
pugs.genetic.resolve_halo(simulation, halo_id) -> (Snapshot, int)
```

Find the snapshot and AHF id a catalog `halo_id` refers to. Since `halo_id` is
`snapshot_index * 2**32 + finder_id`, this is arithmetic rather than a lookup.

---

## `build_param_file`

```python
pugs.genetic.build_param_file(
    filename="genetIC_zoom.txt",
    outname="PUGS_zoom",
    base_grid=2048,
    zoom_grid=[10, 2048],
    autopad=1,
    subsample=8,
    template="../inputs/zoom_template.txt",
)
```

Fill the `inputs/zoom_template.txt` template with per-halo parameters.

| Name | Type | Default | Description |
|---|---|---|---|
| `filename` | `str` | `"genetIC_zoom.txt"` | Output parameter file |
| `outname` | `str` | `"PUGS_zoom"` | Prefix for GenetIC output files |
| `base_grid` | `int` | `2048` | Base grid resolution (match the volume IC) |
| `zoom_grid` | `list[int, int]` | `[10, 2048]` | `[high_res_cells, coarse_cells]` |
| `autopad` | `int` | `1` | Cells of padding around the zoom region |
| `subsample` | `int` | `8` | Subsampling of the coarse outer region |
| `template` | `str` | `"../inputs/zoom_template.txt"` | Template to fill |

Then run GenetIC:

```bash
docker run --rm -v `pwd`:/w/ apontzen/genetic:1.5.0 /w/genetIC_zoom.txt
```

---

## Full zoom workflow

```python
import os

from pugs.genetic import build_param_file, write_particle_ids
from pugs.io import read_catalog

SIM = "/data/sims/NUGS2048"

# 1. Pick the third most massive halo at z=0
catalog = read_catalog("NUGS2048_catalog", columns=["halo_id", "snapshot", "M200", "z50_mass"])
table = catalog.table  # filter it however you like (pandas, polars, DuckDB, ...)

# 2. Write its particle ids, out to 8 virial radii
os.makedirs("zoom_halo3", exist_ok=True)
write_particle_ids(SIM, halo_id, "zoom_halo3/id_file.txt", radius_factor=8)

# 3. Write the GenetIC parameter file
build_param_file(
    filename="zoom_halo3/genetIC_zoom.txt",
    outname="halo3_zoom",
    base_grid=2048,
    zoom_grid=[10, 2048],
)

# 4. Run GenetIC
os.system(
    "docker run --rm -v $(pwd)/zoom_halo3:/w/ "
    "--user $(id -u):$(id -g) apontzen/genetic:1.5.0 /w/genetIC_zoom.txt"
)
```
