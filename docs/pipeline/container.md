---
title: Container Builder
---

# Stage 3: Containerization

**Directory:** `builder/container/`

Bundles the full PUGS pipeline — GenetIC, the scientific Python stack, the
`pugs` package, and the pre-computed Parquet halo catalog — into a single
Apptainer/Singularity image that can be moved to any HPC site without
rebuilding the dependency tree.

---

## Overview

| Property | Value |
|---|---|
| Format | Apptainer/Singularity SIF (squashfs) |
| Base image | `python:3.13-bookworm` |
| GenetIC version | `v1.5.0` (kept in lockstep with `builder/volume_ic/build.sh`) |
| Bundled catalog path | `/opt/PUGS/data/catalog` |
| Code root | `/opt/PUGS/` |
| Default output | `builder/container/pugs.sif` |

---

## Files in `builder/container/`

| File | Purpose |
|---|---|
| `pugs.def` | Apptainer/Singularity definition file |
| `build_container.sh` | Wrapper that stages the catalog and invokes the builder |
| `catalog.placeholder/` | Inert placeholder bundled when no real catalog is supplied |
| `.gitignore` | Excludes built `*.sif` images and the staged catalog |

---

## Build flow

```
  builder/container/pugs.def
            +
   pugs/, inputs/, builder/*, ...
            +
  user-supplied catalog dir (or catalog.placeholder/)
                       │
                       ▼
        build_container.sh -c /path/to/catalog
                       │
                       │   (stages catalog at _catalog.staged)
                       ▼
              apptainer build pugs.sif
                       │
                       ▼
                  pugs.sif
       (portable, content-addressed SIF)
```

---

## Build script

```{code-block} bash
./build_container.sh                        # Bundle the placeholder catalog
./build_container.sh -c /path/to/catalog    # Bundle the real catalog
./build_container.sh -o /tmp/pugs.sif       # Custom output path
./build_container.sh -h                     # Print help
```

### What the script does

1. **Pick a builder** — auto-detects `apptainer` or `singularity` on
   `PATH`, preferring `apptainer`.
2. **Stage the catalog** — copies the user-supplied catalog directory (or
   `catalog.placeholder/`) to `_catalog.staged`. The `.def` file's `%files`
   section references this fixed name, so the recipe itself stays free of
   user-specific paths.
3. **Warn on placeholder** — when no `-c` was given, prints a warning that
   the bundled catalog holds no halo data.
4. **Run the builder** — `apptainer build --force <output> pugs.def`.
5. **Clean up** — removes `_catalog.staged` via an `EXIT` trap, leaving
   only `catalog.placeholder/` in the build directory.

---

## Catalog staging

The halo catalog is baked into the SIF at build time, not mounted at
runtime. Two reasons:

- HPC sites often have inconsistent paths for shared resources; baking
  the catalog in sidesteps mount-path coordination.
- The Apptainer SIF format is content-addressed, so two SIFs built from
  different catalogs are bitwise distinct — a free reproducibility check.

### Placeholder mode

If no `-c` argument is given, the `catalog.placeholder/` directory (holding
only a README) is bundled in place of the real catalog. This lets the build
pipeline (and the CI workflow) succeed end-to-end before the $2048^3$ catalog
is ready.

`read_catalog` against the placeholder raises `FileNotFoundError` because it
contains no `halos_*.parquet` files — by design, as a loud failure mode.

### Replacing with the real catalog

Once the $2048^3$ catalog has been built (see [Halo Catalog](catalog.md)):

```bash
cd builder/container
./build_container.sh -c /path/to/NUGS2048_catalog
```

The resulting `pugs.sif` will have the full halo catalog at
`/opt/PUGS/data/catalog`.

---

## Container contents

### System packages

Installed via `apt-get` in the `%post` section:
`build-essential`, `cmake`, `gfortran`, `git`, `libfftw3-dev`,
`libgsl-dev`, `libhdf5-dev`, `pkg-config`, `wget`, `ca-certificates`.

### GenetIC

Cloned from [pynbody/genetIC](https://github.com/pynbody/genetIC) at the
pinned tag `v1.5.0` and built from source against the system FFTW and GSL.
The resulting binary is installed to `/usr/local/bin/genetIC`.

```{warning}
The pinned GenetIC version must match the Docker tag used in
`builder/volume_ic/build.sh` (currently `apontzen/genetic:1.5.0`).
Bump `GENETIC_VERSION` in `pugs.def` in lockstep with any change to
that Docker tag, otherwise zoom ICs generated inside the container
will not be bitwise-compatible with the bundled volume.
```

### Python environment

The `pugs` package is installed via `pip install /opt/PUGS[dev]`, which
pulls in the full `dev` dependency set:

- Runtime: `pynbody`, `pyarrow`, `numpy`
- Tooling: `camb==1.5.5` (for transfer-function regeneration),
  `pytest==9.0.2`, `requests`, plus the project linters

### Repository content bundled

The `%files` section copies a curated subset of the repo:

- `pugs/` — the `pugs` Python package
- `inputs/` — CAMB and GenetIC parameter files
- `pyproject.toml`, `README.md`
- `builder/volume_ic/` — `build.sh`, `fix_header.py`, `.checksums.txt`
- `builder/catalog/` — `build_catalog.sh`, `test_vars`

`builder/container/` is intentionally not copied (it is itself the
build context).

---

## Runtime usage

### Environment variables (pre-set inside the SIF)

| Variable | Value |
|---|---|
| `PUGS_CATALOG` | `/opt/PUGS/data/catalog` |
| `PUGS_SIMULATION` | `/data/sims` |
| `PUGS_HOME` | `/opt/PUGS` |
| `PYTHONUNBUFFERED` | `1` |

### Interactive shell

```bash
apptainer run pugs.sif
```

Drops into `bash` with the environment above already set.

### One-shot commands

```bash
# Summarise the bundled catalog
apptainer exec pugs.sif python -c \
    "import os; from pugs.io import read_catalog; \
     c = read_catalog(os.environ['PUGS_CATALOG']); \
     print(len(c), 'halos from', c.provenance['simulation'])"

# Or query it with any Parquet reader
apptainer exec pugs.sif python -c \
    "import duckdb, os; duckdb.sql(f\"select count(*) from '{os.environ['PUGS_CATALOG']}/halos_*.parquet'\").show()"

apptainer exec pugs.sif genetIC /path/to/params.txt
```

### Bind-mounting simulation data

Snapshot data lives outside the SIF — it is far too large to bundle. Mount
it at the path `PUGS_SIMULATION` points to:

```bash
apptainer exec -B /scratch/sims:/data/sims pugs.sif \
    python -m pugs.export /data/sims/NUGS2048 -o /data/catalog
```

Multiple binds compose normally:

```bash
apptainer exec \
    -B /scratch/sims:/data/sims \
    -B /scratch/zoom_ics:/data/zoom_ics \
    pugs.sif python my_zoom_pipeline.py
```

---

## Continuous integration

`.github/workflows/build_container.yml` builds the container with the
placeholder catalog on every push and pull request, then smoke-tests the
result by:

- Importing `pynbody` and `pugs`
- Confirming `tangos` is *not* importable, so the dependency cannot creep back
- Confirming `genetIC` is on `PATH`
- Confirming the bundled catalog is present at `/opt/PUGS/data/catalog`
- Confirming each environment variable in the table above matches its
  expected value

The build itself runs ~12–18 minutes on `ubuntu-latest` — most of that is
compiling GenetIC and CAMB's Fortran extension.
