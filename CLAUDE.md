# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PUGS (Portable Universal Galaxy Sampler) is a Python package for generating and processing cosmological N-body simulation data. It provides two main pipelines: (1) generating a large collisionless dark matter volume simulation, and (2) building a TANGOS halo catalog database from simulation snapshots.

## Commands

### Install
```bash
pip install .[dev]
```

### Linting
```bash
black --check .
flake8 .
isort --check-only .
```
Line length is 100 chars. `black` and `isort` (Black profile) are enforced via pre-commit hooks.

### Tests
Tests require the NUGS128 test simulation data and a configured TANGOS database. In CI, the test data is downloaded from nbody.shop; locally, you need `config_vars` set first:
```bash
cd builder/tangos_db
source config_vars
pytest
```
To run a single test file:
```bash
pytest builder/tangos_db/tests/test_properties.py
```

### Build Volume IC (test mode)
```bash
cd builder/volume_ic
./build.sh -t    # 128³ particles, ~4 OMP threads
```
Production (`./build.sh` without `-t`) requires 320 GB RAM and ~90 minutes.

### Build TANGOS Database
```bash
cd builder/tangos_db
./build_tangos.sh    # Requires TANGOS_DB_CONNECTION and TANGOS_SIMULATION_FOLDER env vars
```

## Architecture

### Two-Stage Pipeline

**Stage 1 — Volume IC** (`builder/volume_ic/`): CAMB (transfer function from `inputs/planck_2018_CAMB.ini`) → GenetIC (`inputs/genetIC_volume.txt`, 2048³ particles in a 50 h⁻¹ Mpc box at z=200) → `fix_header.py` (patches 32-bit tipsy header integer overflow) → `DM_volume.tipsy`

**Stage 2 — TANGOS DB** (`builder/tangos_db/`): N-body snapshots → TANGOS load simulation → add halos → import properties (from `pugs.properties`) → build merger trees

### `pugs` Package (`pugs/`)

**`pugs.properties`** — Seven TANGOS property classes registered as entry points under `tangos.property_modules`. Each extends a TANGOS base class and defines a `region_specification()` (sphere) and optional `requires_property()` for prerequisites.

- `Radius500c` / `Mass500c` / `Mass200c`: Virial radius and mass at 500× and 200× critical density
- `StoreIords` / `GetIords`: Store/retrieve particle IDs within 8 virial radii as zlib-compressed bytes
- `StoreIordIndices`: Integer index boundaries for particles at each integer multiple of the virial radius (1–7 R_vir), enabling fast slicing without decompression
- `MassPercentileRedshifts`: Redshifts when the halo assembled 25%, 50%, and 75% of its final mass
- `MergerHistory`: Number of major mergers (`N_mm`) and redshift of the last major merger (`z_lmm`)

Expensive properties (`StoreIords`, `StoreIordIndices`, `MergerHistory`) are only written at z=0 (final timestep), controlled by checking `halo.timestep.extension` against the simulation's last snapshot.

**`pugs.genetic`** — Helpers for zoom-in IC generation:
- `write_particle_ids(halo_id, filename)`: Export particle IDs from TANGOS to a file for GenetIC input
- `build_param_file(...)`: Fill `inputs/zoom_template.txt` with per-halo parameters

### Entry Points

`pugs.properties` classes are auto-discovered by TANGOS via the `[project.entry-points."tangos.property_modules"]` section in `pyproject.toml`. Adding a new property class requires registering it there.

### Test Structure

`builder/tangos_db/tests/` contains two test modules:
- `test_counts.py`: Structural validation (number of simulations, snapshots, halos)
- `test_properties.py`: Physics validation — compares TANGOS-stored values against direct pynbody calculations for radii, masses, particle ID round-trips, assembly history, and merger counts
