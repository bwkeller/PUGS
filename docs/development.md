---
title: Development Guide
---

# Development Guide

## Setup

```bash
git clone https://github.com/bwkeller/PUGS.git
cd PUGS
pip install .[dev]
pre-commit install
```

---

## Code style

PUGS enforces consistent style via pre-commit hooks that run automatically on
`git commit`:

| Tool | Purpose | Config |
|---|---|---|
| [Black](https://black.readthedocs.io) | Code formatter | `pyproject.toml` |
| [isort](https://pycqa.github.io/isort/) | Import sorter (Black profile) | `pyproject.toml` |
| [flake8](https://flake8.pycqa.org) | Linter | `.flake8` |

Line length is **100 characters** throughout.

Run the checks manually:

```bash
black --check .
flake8 .
isort --check-only .
```

---

## Tests

Tests require the NUGS128 test simulation. In CI the data is downloaded from
`nbody.shop`; locally, point `PUGS_SIMULATION` at it (or place it at
`~/data/NUGS128`). No catalog needs to exist beforehand — the fixtures build
one into a temporary directory.

```bash
cd builder/catalog
source test_vars
pytest                                                   # All tests
pytest tests/test_halo_properties.py                     # Single file
pytest tests/test_halo_properties.py::test_finder_mass_is_the_member_particle_mass
```

| File | What it covers |
|---|---|
| `test_ahf.py` | The AHF file readers, and the assumptions they rest on |
| `test_simulation.py` | Snapshot discovery, ordering, and the `halo_id` scheme |
| `test_halo_properties.py` | Physics, checked against direct particle computations |
| `test_merger_forest.py` | Tree structure, merger counts, assembly redshifts |
| `test_catalog.py` | Parquet output, schema, units, provenance |

The physics tests deliberately recompute each quantity from the particle data
rather than comparing against another stored value, so a bug in the measurement
code cannot agree with itself.

### Test data

The tests use the NUGS128 simulation, a 128³-particle version of the PUGS
volume intended for fast CI runs. The snapshots can be downloaded from:

```
https://nbody.shop/NUGS128.tar.gz
```

After extracting, set `PUGS_SIMULATION` to the simulation directory itself
(e.g. `$HOME/data/NUGS128`).

---

## Adding a new column

1. Add a `[columns.<name>]` block to `pugs/columns.toml` giving its `dtype`,
   `unit`, `system` and `description`. The writer raises
   `UndocumentedColumnError` if you skip this, so the data dictionary cannot
   fall behind the data.

2. Add the name to the right tuple in `pugs/export.py`:

   | Tuple | For |
   |---|---|
   | `AHF_COLUMNS` | A column copied verbatim from `AHF_halos` |
   | `MEASURED_COLUMNS` | Something measured from particle data |
   | `TREE_COLUMNS` | Something derived from the merger forest |
   | `DERIVED_COLUMNS` | Arithmetic on columns already present |

3. For a measured column, implement it in `pugs/halo_properties.py` and return
   it from `HaloMeasurement.as_dict()`. A vector quantity is written as one
   scalar column per component (as `shrink_center_x/y/z`), each documented
   separately.

4. Bump `schema_version` in `columns.toml` **only** if you renamed or removed a
   column, or changed a unit. Adding a column is backwards compatible.

5. Add tests to `builder/catalog/tests/`.

### Cost

Columns read from AHF or derived from the forest are effectively free — the
files are read once per snapshot regardless. Columns measured from particle
data cost one pynbody pass per halo and dominate the runtime; `--no-measure`
skips them entirely when they are not needed.

---

## CI/CD

Three GitHub Actions workflows run on every push and pull request:

| Workflow | File | What it tests |
|---|---|---|
| Linter | `linter.yml` | Black, flake8, isort, shellcheck |
| Volume IC | `build_volume_ic.yml` | 128³ IC build + MD5 checksums |
| Catalog | `build_catalog.yml` | Full catalog build + pytest |
| Container | `build_container.yml` | SIF build + smoke test |

The catalog workflow downloads NUGS128 from `nbody.shop`, builds the catalog,
and runs the test suite in a single job.

---

## Versioning

The package version is derived dynamically from the git tag (via
`setuptools-scm`). Tag releases following [Semantic Versioning](https://semver.org):

```bash
git tag v1.2.3
git push origin v1.2.3
```

---

## Dependencies

Core runtime dependencies (pinned in `pyproject.toml`):

| Package | Purpose |
|---|---|
| `pynbody` | Snapshot reading and N-body particle analysis |
| `pyarrow` | Parquet output and the Arrow schema/metadata |
| `numpy` | Array operations |

Reading a finished catalog requires none of these — it is plain Parquet.

Dev/build extras:

| Package | Version | Purpose |
|---|---|---|
| `camb` | 1.5.5 | Transfer function generation |
| `black` | 25.9.0 | Code formatting |
| `flake8` | 7.3.0 | Linting |
| `isort` | 7.0.0 | Import sorting |
| `pytest` | 9.0.2 | Testing |
| `pytest-timeout` | 2.4.0 | Fails a hung test with a traceback instead of stalling the run |
| `duckdb` | >=1.0 | Verifying the catalog through a non-Arrow reader |
| `requests` | 2.32.5 | HTTP downloads (CI) |
| `shellcheck-py` | latest | Shell script linting |
| `pre-commit` | latest | Git hooks |

### Test timeouts

`pyproject.toml` sets a per-test timeout of 3600 s using
`timeout_method = "thread"`. Both values are deliberate:

- The limit is generous because nearly all of the suite's cost is the
  session-scoped catalog fixture, which lands on whichever test triggers it
  first — around 2.5 minutes on an idle machine, but closer to 20 on a shared
  one under load. A tighter limit would fail on contention rather than on a
  real hang.
- The `thread` method is used rather than the default `signal` because a
  signal-based timeout cannot interrupt a C call: it only fires once Python
  regains control. Measured against a single long BLAS call with a 5 s limit,
  `signal` reported at 9.7 s (when the call finished on its own) while
  `thread` reported at 5.2 s. With a genuine hang inside pynbody or pyarrow,
  `signal` would never fire at all.

Note that `thread` aborts the whole session rather than failing one test and
continuing, which is the right trade when the run is already doomed.

GenetIC itself is not a Python dependency — it is run via the official Docker
image `apontzen/genetic:1.5.0`.
