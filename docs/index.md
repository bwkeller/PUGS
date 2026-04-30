# PUGS: Portable Universal Galaxy Sampler

[![build_volume_ic](https://github.com/bwkeller/PUGS/actions/workflows/build_volume_ic.yml/badge.svg)](https://github.com/bwkeller/PUGS/actions/workflows/build_volume_ic.yml)
[![Tangos DB Build Test](https://github.com/bwkeller/PUGS/actions/workflows/build_tangos_db.yml/badge.svg)](https://github.com/bwkeller/PUGS/actions/workflows/build_tangos_db.yml)
[![linter](https://github.com/bwkeller/PUGS/actions/workflows/linter.yml/badge.svg)](https://github.com/bwkeller/PUGS/actions/workflows/linter.yml)

PUGS is a Python package for generating and analyzing cosmological N-body simulation data.
It provides two integrated pipelines:

1. **Volume IC generation** — produces a 50 $h^{-1}$ Mpc, $2048^3$-particle collisionless
   dark matter volume using [CAMB](https://camb.info) and [GenetIC](https://github.com/pynbody/genetIC),
   seeded with Planck 2018 cosmology.

2. **TANGOS halo catalog** — ingests N-body snapshots into a
   [TANGOS](https://tangos.readthedocs.io) database and computes a rich set of halo
   properties, merger trees, and assembly histories.

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

:::{grid-item-card} The Two-Stage Pipeline
:link: pipeline/index
:link-type: doc
Understand the data flow from cosmological parameters to halo catalogs.
:::

:::{grid-item-card} API Reference
:link: api/index
:link-type: doc
Full reference for `pugs.properties` and `pugs.genetic`.
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
├── properties.py   # TANGOS halo property calculators
└── genetic.py      # Zoom-in IC generation helpers

builder/
├── volume_ic/      # Stage 1 — Volume initial conditions
└── tangos_db/      # Stage 2 — TANGOS database construction

inputs/
├── planck_2018_CAMB.ini         # CAMB cosmology configuration
├── genetIC_volume.txt           # GenetIC volume parameters
└── zoom_template.txt            # Zoom IC parameter template
```
