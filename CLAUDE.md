# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PUGS (Portable Universal Galaxy Sampler) is a Python package for generating and processing cosmological N-body simulation data. It provides two pipelines: (1) generating a large collisionless dark matter volume simulation, and (2) building a Parquet halo catalog from simulation snapshots and their AHF halo-finder output.

There is no database anywhere in the pipeline. PUGS used to build a TANGOS SQLite database; every column is now derived directly from the AHF files and the snapshots. Reading a finished catalog requires nothing but a Parquet reader.

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
Line length is 100 chars. `black` and `isort` (Black profile) are enforced via pre-commit hooks. Shell scripts are checked with `shellcheck` and `shfmt`.

### Tests
Tests require the NUGS128 test simulation. In CI it is downloaded from nbody.shop; locally, point `PUGS_SIMULATION` at it (or place it at `~/data/NUGS128`). No catalog needs to exist first — the fixtures build one into a temporary directory.
```bash
cd builder/catalog
source test_vars
pytest
```
To run a single test file:
```bash
pytest builder/catalog/tests/test_halo_properties.py
```
The full suite takes ~2.5 minutes, dominated by building a catalog once as a session fixture.

### Build Volume IC (test mode)
```bash
cd builder/volume_ic
./build.sh -t    # 128³ particles, ~4 OMP threads
```
Production (`./build.sh` without `-t`) requires 320 GB RAM and ~90 minutes.

### Build the halo catalog
```bash
cd builder/catalog
source test_vars          # or prod_vars
./build_catalog.sh test_vars
```
or directly:
```bash
python -m pugs.export /path/to/NUGS128 --output NUGS128_catalog
```
`--no-measure` skips the particle-derived columns (leaving them null), turning a run of minutes into one of seconds. Env vars: `PUGS_SIMULATION`, `PUGS_NAME`, `PUGS_CATALOG`, `PUGS_MIN_PARTICLES`. Each vars file defers to an already-exported value.

## Architecture

### Three-Stage Pipeline

**Stage 1 — Volume IC** (`builder/volume_ic/`): CAMB (transfer function from `inputs/planck_2018_CAMB.ini`) → GenetIC (`inputs/genetIC_volume.txt`, 2048³ particles in a 50 h⁻¹ Mpc box at z=200) → `fix_header.py` (patches 32-bit tipsy header integer overflow) → `DM_volume.tipsy`

**Stage 2 — Halo catalog** (`builder/catalog/`): snapshots + AHF output → `pugs.export` → one Parquet file per snapshot plus a `provenance.json` sidecar.

**Stage 3 — Container** (`builder/container/`): bundles GenetIC, the Python stack, and a catalog into an Apptainer SIF.

### `pugs` Package (`pugs/`)

**`pugs.ahf`** — Readers for AHF output. `read_halo_table` (the 83-column TSV, via pyarrow's CSV reader so it scales to millions of rows), `read_particle_membership` (CSR-form halo membership), `read_substructure`, and `read_merger_links`.

Two facts the whole pipeline rests on:
- AHF's particle values are **indices into the snapshot**, not `iord` values — the tipsy files carry no `iord` array.
- PUGS reads `AHF_croco`, not `AHF_mtree`. Both hold the same tree edges, but `croco` carries the shared-particle counts and AHF's merit `shared²/(n_desc·n_prog)`, while `mtree` has them stripped out, leaving only an ordering from which a weight can only be invented.

**`pugs.simulation`** — Snapshot discovery and time ordering. `halo_id` is `snapshot_index * 2**32 + finder_id`: deterministic, so a rebuilt catalog refers to the same halos, and decodable by arithmetic via `split_halo_id`.

**`pugs.halo_properties`** — Per-halo physics with pynbody: shrinking-sphere centre, `max_radius`, `R200`/`R500` by bisection on enclosed density, `M200`/`M500`.

`max_radius` is **not** an overdensity radius — it is the distance to the outermost halo-finder particle, and it is what zoom regions are defined in multiples of. It is close to `R200` but not bounded by it, since `R200` is measured over every particle in the region while `max_radius` covers only the particles AHF left bound.

Periodicity is handled two different ways on purpose. Centring wraps explicitly, because `shrink_sphere_center` knows nothing about the box; without it a boundary-straddling halo's centre lands in empty space and the bisection then fails outright. Region selection relies on `pynbody.filt.Sphere` being periodic-aware, avoiding an O(N) wrap per halo.

**`pugs.merger_forest`** — The merger forest from AHF's tree files. A link is kept when **strictly more than** `MIN_SHARED_FRACTION` (0.5) of the progenitor's particles end up in the descendant. Thresholding on merit would be wrong: merit penalises size mismatch, so a small halo swallowed whole scores near zero while having contributed all of itself. The strict >0.5 rule also makes the graph a forest by construction, which is what lets merger counts accumulate in one sweep. The forest uses **every** AHF halo, not just those above the catalog's particle cut.

**`pugs.export`** — Orchestration. `AHF_COLUMNS` / `MEASURED_COLUMNS` / `TREE_COLUMNS` / `DERIVED_COLUMNS` define what is written. Assembly history is written for every snapshot, not only z=0.

**`pugs.particle_ids`** — Stores the particles around each z=0 halo as 10 concentric shells of 0.5 R_vir out to 5 R_vir (the innermost is a sphere), one row per (halo, shell), in `shells_<snapshot>.parquet`.

Bucketing by radius rather than sorting by it is the whole point: the region selection already returns ascending indices, so **no sort is applied** and each shell's ids stay nearly consecutive, which lets Parquet's `DELTA_BINARY_PACKED` encoding store differences instead of values. Measured on NUGS128, that is 7.8x smaller than raw int64 against 4.8x for the old radially-sorted zlib blobs. Dictionary encoding must stay off or there is nothing left to delta.

Opt-in via `PUGS_PARTICLE_IDS=1`; reads centres and `R200` from the existing catalog rather than recomputing.

**`pugs.io`** — The catalog format. Owns the schema (built from `pugs/columns.toml`, never inferred from the data), attaches units as Arrow field metadata plus a file-level JSON mirror for non-Arrow readers, and refuses to read a catalog whose `schema_version` does not match. `read_catalog()` returns a `Catalog` with `.table`, `.unit()`, `.system()`, `.description()`, `.provenance`.

**`pugs/columns.toml`** — The data dictionary: `dtype`, `unit`, `system`, `description` per column, plus `schema_version`. TOML rather than YAML because YAML 1.1 parses `1e10` and `2.775e11` as *strings*. Adding a column to the export requires adding it here first — the writer raises `UndocumentedColumnError` otherwise.

**`pugs.genetic`** — Zoom-in IC helpers. `particle_ids_from_catalog` assembles a region from the stored shells (catalog only, radii on the 0.5 R_vir grid); `particle_ids` selects from the snapshot (any radius). Both use `R200` as the reference and return identical particles for the same factor. `build_param_file` fills `inputs/zoom_template.txt`.

### Two unit systems

Every column declares a `system`, queryable via `catalog.system(name)`:
- `physical` — measured by PUGS with pynbody; `Msol`, `kpc`, h divided out.
- `finder` — read verbatim from AHF; comoving, h-scaled (`Msol / littleh`, `kpc / littleh`).
- `pugs` — bookkeeping and derived columns.

`M200` and `Mhalo` are both halo masses and are **not** directly comparable.

### Test Structure

`builder/catalog/tests/`:
- `test_ahf.py`: The AHF readers and the assumptions they rest on
- `test_simulation.py`: Snapshot discovery, ordering, the `halo_id` scheme
- `test_halo_properties.py`: Physics, recomputed from particle data rather than compared against another stored value
- `test_merger_forest.py`: Tree structure, merger counts, assembly redshifts
- `test_particle_ids.py`: Shell geometry, that the ids really are delta-encoded on disk, and that the shells partition the 5 R_vir sphere
- `test_genetic.py`: Zoom-region selection, and that the catalog and snapshot paths agree
- `test_catalog.py`: Parquet output, schema uniformity, units reaching Arrow and non-Arrow (DuckDB) readers, the schema-version guard

Two production hazards have regression tests because the 128³ test data cannot expose them on its own: AHF's separator differs between runs (tab-padded in NUGS128, single-space in NUGS2048), and the `hostHalo` "no host" sentinel is -1 or 0 depending on whether the run numbers ids from 0 or 1.
