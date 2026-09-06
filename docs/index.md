# PUGS: Portable Universal Galaxy Sampler

[![build_volume_ic](https://github.com/bwkeller/PUGS/actions/workflows/build_volume_ic.yml/badge.svg)](https://github.com/bwkeller/PUGS/actions/workflows/build_volume_ic.yml)
[![Catalog Build Test](https://github.com/bwkeller/PUGS/actions/workflows/build_catalog.yml/badge.svg)](https://github.com/bwkeller/PUGS/actions/workflows/build_catalog.yml)
[![Container Build Test](https://github.com/bwkeller/PUGS/actions/workflows/build_container.yml/badge.svg)](https://github.com/bwkeller/PUGS/actions/workflows/build_container.yml)
[![linter](https://github.com/bwkeller/PUGS/actions/workflows/linter.yml/badge.svg)](https://github.com/bwkeller/PUGS/actions/workflows/linter.yml)

PUGS is a Python package for generating and analyzing cosmological N-body simulation data.
It provides two integrated pipelines:

1. **Volume IC generation** — produces a 50 $h^{-1}$ Mpc, $2048^3$-particle collisionless
   dark matter volume using [CAMB](https://camb.info) and [GenetIC](https://github.com/pynbody/genetIC),
   seeded with Planck 2018 cosmology.

2. **Halo catalog** — turns N-body snapshots and their AHF halo-finder output
   into a Parquet catalog with halo properties, merger trees, and assembly
   histories. Every column is derived directly from the AHF files and the
   snapshots; there is no database in the pipeline, and reading the result
   needs nothing but a Parquet reader.

The `pugs` Python package also exposes helpers for generating **zoom-in initial
conditions** for individual halos identified in the volume.

---

## Quick navigation

::::{grid} 1 2 2 2

:::{grid-item-card} Getting Started
:link: getting-started
:link-type: doc
Install PUGS and run your first build.
:::

:::{grid-item-card} The Pipeline
:link: pipeline/index
:link-type: doc
Understand the data flow from cosmological parameters to halo catalogs.
:::

:::{grid-item-card} API Reference
:link: api/index
:link-type: doc
Full reference for building and reading the halo catalog.
:::

:::{grid-item-card} Cosmological Background
:link: cosmology
:link-type: doc
Planck 2018 parameters, power spectrum, and transfer function details.
:::

::::

---

## Package layout

```{toctree}
:hidden:
:maxdepth: 2

getting-started
pipeline/index
api/index
cosmology
NUGS2048
development
```

```
pugs/
├── ahf.py              # Readers for AHF halo-finder output
├── simulation.py       # Snapshot discovery, ordering, halo ids
├── halo_properties.py  # Per-halo physics measured with pynbody
├── merger_forest.py    # Merger trees from AHF's own tree files
├── export.py           # Builds the Parquet catalog
├── io.py               # Catalog format: schema, units, provenance
├── columns.toml        # The data dictionary
└── genetic.py          # Zoom-in IC generation helpers

builder/
├── volume_ic/      # Stage 1 — Volume initial conditions
├── catalog/        # Stage 2 — Parquet halo catalog
└── container/      # Stage 3 — Apptainer/Singularity packaging

inputs/
├── planck_2018_CAMB.ini         # CAMB cosmology configuration
├── genetIC_volume.txt           # GenetIC volume parameters
└── zoom_template.txt            # Zoom IC parameter template
```
