---
title: Getting Started
---

# Getting Started

## Requirements

- **Python 3.9+** (tested through 3.13)
- **Docker** (tested with 28.5.1) — required for running GenetIC when building
  the volume IC
- A Unix-like shell (Linux or macOS)

For the full production volume build you also need at least **320 GB RAM** and
roughly **96 CPU cores** for a reasonable wall time (~90 minutes).

---

## Installation

Clone the repository and install the package with its development dependencies:

```bash
git clone https://github.com/bwkeller/PUGS.git
cd PUGS
pip install .[dev]
```

`pip install .` (without `[dev]`) is sufficient for building and reading halo
catalogs; the `[dev]` extras add the build tools (CAMB, testing, linters).

Reading an existing catalog needs no PUGS install at all — it is plain Parquet,
readable by DuckDB, polars, pandas or any other Arrow-aware tool. PUGS is only
required to *build* one, or to get the units and provenance attached
automatically.

---

## Quickstart: test volume

The fastest way to verify your installation is to build the small 128³-particle
test volume:

```bash
cd builder/volume_ic
./build.sh -t
```

This runs CAMB (cosmological transfer function) and GenetIC (via Docker) in
test mode. The build takes a few minutes and produces files in
`builder/volume_ic/outputs/`. You can verify bitwise reproducibility against
the checked-in MD5 checksums:

```bash
md5sum -c .checksums.txt
```

---

## Quickstart: halo catalog

The catalog build requires a directory holding N-body snapshots together with
the AHF halo-finder output beside them (`.AHF_halos`, `.AHF_particles` and
`.AHF_croco` per snapshot).

Edit `builder/catalog/test_vars` to match your paths, then:

```bash
cd builder/catalog
source test_vars
./build_catalog.sh test_vars
```

This writes one Parquet file per snapshot plus a `provenance.json` sidecar.
Add `--no-measure` to `pugs.export` to skip the particle-derived columns, which
turns a run of minutes into one of seconds.

Once the catalog is built, run the test suite to confirm correctness:

```bash
pytest
```

---

## Using the Python API

Read a catalog, with units and provenance attached:

```python
from pugs.io import read_catalog

catalog = read_catalog("NUGS128_catalog")
print(len(catalog), "halos")
print(catalog.unit("M200"), catalog.system("M200"))   # 'Msol' 'physical'
print(catalog.unit("Mhalo"), catalog.system("Mhalo")) # 'Msol / littleh' 'finder'
```

Or query it with DuckDB, no PUGS import needed:

```python
import duckdb

duckdb.sql("SELECT snapshot, count(*) FROM 'NUGS128_catalog/halos_*.parquet' GROUP BY snapshot")
```

Generate zoom-in IC inputs for a halo:

```python
import pugs.genetic as genetic

# Export particle IDs out to 8 virial radii (needed by GenetIC for zoom ICs)
genetic.write_particle_ids("/path/to/NUGS128", halo_id, "halo1_ids.txt", radius_factor=8)

# Generate a GenetIC parameter file for the zoom simulation
genetic.build_param_file(
    filename="genetIC_zoom.txt",
    outname="halo1_zoom",
    base_grid=2048,
    zoom_grid=[10, 2048],
)
```

See the [API reference](api/index.md) for full details.

---

## Environment variables

| Variable | Description | Example |
|---|---|---|
| `PUGS_SIMULATION` | Directory holding the snapshots and AHF output | `$HOME/data/NUGS128` |
| `PUGS_NAME` | Simulation name recorded in the catalog | `NUGS128` |
| `PUGS_CATALOG` | Output catalog directory | `NUGS128_catalog` |
| `PUGS_MIN_PARTICLES` | Minimum particles for a halo to be kept | `500` |

These are set by sourcing `builder/catalog/test_vars` (or `prod_vars`). Each
defers to a value already exported, so you can override any of them from the
environment.
