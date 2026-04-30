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

`pip install .` (without `[dev]`) is sufficient for using the `pugs` Python
package to query an existing TANGOS database; the `[dev]` extras add the build
tools (CAMB, testing, linters).

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

## Quickstart: TANGOS database

The TANGOS database build requires:

- A set of N-body simulation snapshots stored under `$TANGOS_SIMULATION_FOLDER`
- The environment variable `TANGOS_DB_CONNECTION` pointing to a SQLite file path

Edit `builder/tangos_db/config_vars` to match your paths, then:

```bash
cd builder/tangos_db
source config_vars
./build_tangos.sh
```

Once the database is built, run the test suite to confirm correctness:

```bash
pytest
```

---

## Using the Python API

After installing the package and sourcing `config_vars`, you can query halo
data directly:

```python
import tangos
import pugs.genetic as genetic

# Export particle IDs for halo 1 to disk (needed by GenetIC for zoom ICs)
genetic.write_particle_ids(1, filename="halo1_ids.txt")

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
| `TANGOS_DB_CONNECTION` | Path to the SQLite database file | `pugs.db` |
| `TANGOS_SIMULATION_FOLDER` | Parent directory containing snapshot folders | `$HOME/data` |
| `SIM` | Simulation folder name inside `TANGOS_SIMULATION_FOLDER` | `NUGS128` |

These are set by sourcing `builder/tangos_db/config_vars`.
