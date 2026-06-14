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

Tests require the NUGS128 test simulation and a configured TANGOS database.
In CI the data is downloaded from `nbody.shop`; locally you need to have run
the database build first.

```bash
cd builder/tangos_db
source config_vars
pytest                                          # All tests
pytest tests/test_properties.py                 # Single file
pytest tests/test_properties.py::test_virial_radii  # Single test
```

### Test data

The tests use the NUGS128 simulation, a 128³-particle version of the PUGS
volume intended for fast CI runs. The snapshots can be downloaded from:

```
https://nbody.shop/NUGS128.tar.gz
```

After extracting, set `TANGOS_SIMULATION_FOLDER` to the parent directory
and `SIM=NUGS128`.

---

## Adding a new property

1. Add a class to `pugs/properties.py` that inherits from the appropriate
   TANGOS base class.

2. Set `names` to the TANGOS property name (or a tuple of names if the
   class returns multiple values).

3. Implement `calculate(self, particle_data, existing_properties)`.

4. Optionally implement:
   - `region_specification(existing_properties)` — returns a pynbody filter
     defining which particles to load.
   - `requires_property(self)` — returns a list of prerequisite property names.
   - `preloop(cls, sim, db_timestep)` — classmethod called once per timestep
     to cache data shared across halos (avoids repeated TANGOS queries).

5. Register the class in `pyproject.toml`:

   ```toml
   [project.entry-points."tangos.property_modules"]
   properties = "pugs.properties"
   ```

   Because the entire module is registered, any new class in
   `pugs/properties.py` is automatically discovered — no additional entry
   needed.

6. Add corresponding tests to `builder/tangos_db/tests/test_properties.py`.

### Expensive properties (z=0 only)

Properties that traverse the full merger tree or store large arrays should
be written only at the final snapshot. The convention is to add them to the
`tangos write ... --latest` command in `build_tangos.sh`.

---

## CI/CD

Three GitHub Actions workflows run on every push and pull request:

| Workflow | File | What it tests |
|---|---|---|
| Linter | `linter.yml` | Black, flake8, isort |
| Volume IC | `build_volume_ic.yml` | 128³ IC build + MD5 checksums |
| TANGOS DB | `build_tangos_db.yml` | Full DB build + pytest |

The TANGOS DB workflow downloads NUGS128 from `nbody.shop`, builds the
database, and runs the test suite in a single job.

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

| Package | Version | Purpose |
|---|---|---|
| `tangos` | 1.10.0 | Halo database and merger trees |
| `pynbody` | 2.3.2 | N-body particle analysis |
| `numpy` | (transitive) | Array operations |

Dev/build extras:

| Package | Version | Purpose |
|---|---|---|
| `camb` | 1.5.5 | Transfer function generation |
| `black` | 25.9.0 | Code formatting |
| `flake8` | 7.3.0 | Linting |
| `isort` | 7.0.0 | Import sorting |
| `pytest` | 9.0.2 | Testing |
| `requests` | 2.32.5 | HTTP downloads (CI) |
| `pre-commit` | latest | Git hooks |

GenetIC itself is not a Python dependency — it is run via the official Docker
image `apontzen/genetic:1.5.0`.
