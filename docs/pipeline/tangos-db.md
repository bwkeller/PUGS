---
title: TANGOS Database
---

# Stage 2: TANGOS Halo Catalog

**Directory:** `builder/tangos_db/`

This stage reads N-body simulation snapshots and builds a TANGOS SQLite
database containing halos with a full set of precomputed physical properties.

---

## Prerequisites

- N-body snapshots in a directory accessible as `$TANGOS_SIMULATION_FOLDER/$SIM`
- `TANGOS_DB_CONNECTION` and `TANGOS_SIMULATION_FOLDER` environment variables set
- `pugs` package installed (so its entry points are registered with TANGOS)

Set all three by sourcing:

```bash
source builder/tangos_db/config_vars
```

Default values in `config_vars`:

```bash
export TANGOS_DB_CONNECTION=pugs.db
export TANGOS_SIMULATION_FOLDER=$HOME/data
export SIM=NUGS128
```

---

## Build script (`build_tangos.sh`)

Run from `builder/tangos_db/`:

```bash
./build_tangos.sh
```

The script runs the following TANGOS commands in sequence.

### Step 1 — Load simulation

```bash
tangos add --min-particles 500 $SIM
```

Registers the simulation and all its snapshots in the database. The
`--min-particles 500` threshold prevents spurious very small halos from being
cataloged.

### Step 2 — Import properties

```bash
tangos import-properties
```

Discovers all registered TANGOS property modules (including `pugs.properties`
via the entry point) and computes every property for every halo at every
snapshot. Prerequisites are handled automatically.

### Step 3 — Import merger trees

```bash
tangos import-ahf-trees
```

Reads the AHF-format merger tree files alongside the snapshots and links halos
across timesteps.

### Step 4 — Remove duplicates

```bash
tangos remove-duplicates
```

Cleans up any halo entries that were registered more than once.

### Step 5 — Write fast-access properties (all timesteps)

```bash
tangos write finder_mass shrink_center max_radius R500 M200 M500 \
    --for $SIM --with-prerequisites
```

Pre-computes and caches the lightweight per-halo properties at every snapshot.

### Step 6 — Write expensive properties (z=0 only)

```bash
tangos write zlib_ids Rvir_indices N_mm z_lmm z25_mass z50_mass z75_mass \
    --with-prerequisites --for $SIM --latest
```

Computes the storage-intensive and merger-tree-dependent properties only at
the final (z=0) timestep.

---

## Computed properties

### Properties computed at every timestep

| TANGOS name | Class | Description |
|---|---|---|
| `finder_mass` | (TANGOS built-in) | Total halo mass from halo finder |
| `shrink_center` | (TANGOS built-in) | Shrinking-sphere center of mass |
| `max_radius` | (TANGOS built-in) | Outer radius of halo |
| `R500` | `Radius500c` | Radius enclosing 500× critical density |
| `M200` | `Mass200c` | Mass within $R_{200c}$ (= `finder_mass`) |
| `M500` | `Mass500c` | Mass within $R_{500c}$ |

### Properties computed at z=0 only

| TANGOS name | Class | Description |
|---|---|---|
| `zlib_ids` | `StoreIords` | zlib-compressed particle IDs sorted by radius, out to 8 $R_\mathrm{vir}$ |
| `Rvir_indices` | `StoreIordIndices` | Index boundaries at 1–7 $R_\mathrm{vir}$ in the sorted ID array |
| `ids` | `GetIords` | Live decompressor for `zlib_ids` (not stored) |
| `z25_mass`, `z50_mass`, `z75_mass` | `MassPercentileRedshifts` | Redshifts of 25%, 50%, 75% mass assembly |
| `N_mm` | `MergerHistory` | Number of major mergers (mass ratio ≥ 1:4) |
| `z_lmm` | `MergerHistory` | Redshift of last major merger (−1 if none) |

---

## Particle ID storage design

`StoreIords` compresses all particle IDs within 8 $R_\mathrm{vir}$ to a single
zlib-compressed blob. Particles are sorted by radius before compression, so
any radial shell $[0, n \cdot R_\mathrm{vir}]$ can be extracted by slicing the
decompressed array at the index stored in `Rvir_indices[n-1]`.

This design avoids storing per-particle float coordinates and allows fast
construction of zoom regions without re-running the N-body code:

```python
import zlib
import numpy as np
import tangos

halo = tangos.get_halo("NUGS128/016/%1")

# Decompress all IDs (out to 8 Rvir)
all_ids = np.frombuffer(
    zlib.decompress(halo["zlib_ids"].tobytes()), dtype=np.int64
)

# Slice to 3 Rvir using the pre-computed index
idx_3rvir = halo["Rvir_indices"][2]   # index for the 3rd shell boundary
ids_within_3rvir = all_ids[:idx_3rvir]
```

Or use the live `ids()` property to get all IDs decompressed automatically:

```python
all_ids = halo.calculate("ids()")
```

---

## Test suite

```bash
cd builder/tangos_db
source config_vars
pytest
```

### `test_counts.py`

Verifies the database was built correctly:

| Test | Assertion |
|---|---|
| `test_sim_count` | Exactly 1 simulation in the database |
| `test_snap_count` | 16 timesteps for NUGS128 |
| `test_halo_count` | 219 halos at the final snapshot |

### `test_properties.py`

Validates computed property values against independent pynbody calculations:

| Test | What it checks |
|---|---|
| `test_property_counts` | All expected properties are present at the right timesteps |
| `test_id_compression` | Decompressed IDs are a superset of direct pynbody particle indices |
| `test_rvir_idx` | `Rvir_indices` boundaries fall in the correct radial shells |
| `test_virial_radii` | `R500 < R200`; `R500` matches pynbody `virial_radius(..., overden=500)` |
| `test_virial_mass` | `M500 < M200`; both within 0.1% of pynbody direct sums |
| `test_mass_percentiles` | `z25_mass ≥ z50_mass ≥ z75_mass > 0` for all halos |
| `test_merger_history` | `z_lmm = −1` iff `N_mm = 0`; `z_lmm > 0` when `N_mm > 0` |
