---
title: pugs.genetic
---

# `pugs.genetic`

Helper functions for generating zoom-in initial conditions for individual
halos identified in the PUGS TANGOS database.

The typical workflow is:

1. Identify the halo ID from TANGOS.
2. Call `write_particle_ids` to export particle indices to disk.
3. Call `build_param_file` to create the GenetIC parameter file.
4. Run GenetIC with the produced files to generate the zoom IC.

---

## Functions

### `write_particle_ids`

```python
pugs.genetic.write_particle_ids(halo_id, filename="id_file.txt")
```

Export the particle IDs belonging to a halo to a plain-text file suitable for
use as GenetIC's `id_file` input.

The function queries the `ids()` live property from TANGOS (which
decompresses `zlib_ids` on the fly) and writes one integer per line.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `halo_id` | `int` | — | TANGOS database ID of the halo |
| `filename` | `str` | `"id_file.txt"` | Path for the output particle ID file |

**Returns:** `None`

**Example**

```python
import tangos
import pugs.genetic as genetic

# Find the halo you want
halo = tangos.get_halo("NUGS128/016/%1")

# Write its particle IDs
genetic.write_particle_ids(halo.id, filename="halo1_ids.txt")
```

---

### `build_param_file`

```python
pugs.genetic.build_param_file(
    filename="genetIC_zoom.txt",
    outname="PUGS_zoom",
    base_grid=2048,
    zoom_grid=[10, 2048],
    autopad=1,
    subsample=8,
)
```

Generate a GenetIC parameter file for a zoom-in simulation by filling the
`inputs/zoom_template.txt` template with per-halo parameters.

The template contains Python `str.format()` placeholders; this function
substitutes them and writes the result to `filename`.

**Parameters**

| Name | Type | Default | Description |
|---|---|---|---|
| `filename` | `str` | `"genetIC_zoom.txt"` | Path for the output parameter file |
| `outname` | `str` | `"PUGS_zoom"` | Prefix for GenetIC output files |
| `base_grid` | `int` | `2048` | Base grid resolution (should match the volume IC) |
| `zoom_grid` | `list[int, int]` | `[10, 2048]` | `[high_res_cells, coarse_cells]` for the nested zoom region |
| `autopad` | `int` | `1` | Number of cells of padding around the zoom region |
| `subsample` | `int` | `8` | Subsampling factor applied to the coarse outer region |

**Returns:** `None`

**Example**

```python
import pugs.genetic as genetic

genetic.build_param_file(
    filename="halo1_zoom.txt",
    outname="halo1",
    base_grid=2048,
    zoom_grid=[10, 2048],
    autopad=1,
    subsample=8,
)
```

Then run GenetIC:

```bash
docker run --rm -v `pwd`:/w/ apontzen/genetic:1.5.0 /w/halo1_zoom.txt
```

---

## Parameter template (`inputs/zoom_template.txt`)

The template sets the same cosmology and transfer function as the volume IC,
and additionally specifies:

- `mapper_relative_to genetIC_volume.txt` — links the zoom IC to the original
  volume by reusing the same random phase realisation.
- `id_file id_file.txt` — the particle ID file written by `write_particle_ids`.
- `centre_output` — centres the output coordinates on the zoom region.

The placeholders replaced by `build_param_file` are:

| Placeholder | Parameter |
|---|---|
| `{outname}` | `outname` |
| `{base_grid}` | `base_grid` |
| `{zoom_grid}` | `zoom_grid` |
| `{autopad}` | `autopad` |
| `{subsample}` | `subsample` |

---

## Full zoom workflow example

```python
import os
import tangos
import pugs.genetic as genetic

# 1. Identify target halo
halo = tangos.get_halo("NUGS128/016/%3")   # third most-massive halo at z=0
print(f"M200 = {halo['M200']:.2e} Msol/h")
print(f"Assembly: z50 = {halo['z50_mass']:.2f}, N_mm = {halo['N_mm']}")

# 2. Write particle IDs
os.makedirs("zoom_halo3", exist_ok=True)
genetic.write_particle_ids(halo.id, filename="zoom_halo3/id_file.txt")

# 3. Write GenetIC parameter file
genetic.build_param_file(
    filename="zoom_halo3/genetIC_zoom.txt",
    outname="halo3_zoom",
    base_grid=2048,
    zoom_grid=[10, 2048],
    autopad=1,
    subsample=8,
)

# 4. Run GenetIC (in zoom_halo3 directory)
os.system(
    "docker run --rm -v $(pwd)/zoom_halo3:/w/ "
    "--user $(id -u):$(id -g) "
    "apontzen/genetic:1.5.0 /w/genetIC_zoom.txt"
)
```
