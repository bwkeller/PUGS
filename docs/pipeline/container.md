---
title: Container Builder
---

# Stage 3: Containerization

**Directory:** `builder/container/`

Bundles the full PUGS pipeline — GenetIC, the scientific Python stack, the
`pugs` package, and the pre-computed TANGOS halo catalog — into a single
Apptainer/Singularity image that can be moved to any HPC site without
rebuilding the dependency tree.

---

## Overview

| Property | Value |
|---|---|
| Format | Apptainer/Singularity SIF (squashfs) |
| Base image | `python:3.13-bookworm` |
| GenetIC version | `v1.5.0` (kept in lockstep with `builder/volume_ic/build.sh`) |
| Bundled catalog path | `/opt/PUGS/data/pugs.db` |
| Code root | `/opt/PUGS/` |
| Default output | `builder/container/pugs.sif` |

---

## Files in `builder/container/`

| File | Purpose |
|---|---|
| `pugs.def` | Apptainer/Singularity definition file |
| `build_container.sh` | Wrapper that stages the catalog and invokes the builder |
| `pugs.db.placeholder` | Inert placeholder bundled when no real catalog is supplied |
| `.gitignore` | Excludes built `*.sif` images and the staged catalog |

---

## Build flow

```
  builder/container/pugs.def
            +
   pugs/, inputs/, builder/*, ...
            +
  user-supplied pugs.db (or pugs.db.placeholder)
                       │
                       ▼
            build_container.sh -d pugs.db
                       │
                       │   (stages catalog at _pugs.db.staged)
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
./build_container.sh                       # Bundle the placeholder catalog
./build_container.sh -d /path/to/pugs.db   # Bundle the real catalog
./build_container.sh -o /tmp/pugs.sif      # Custom output path
./build_container.sh -h                    # Print help
```

### What the script does

1. **Pick a builder** — auto-detects `apptainer` or `singularity` on
   `PATH`, preferring `apptainer`.
2. **Stage the catalog** — copies the user-supplied database (or
   `pugs.db.placeholder`) to `_pugs.db.staged`. The `.def` file's `%files`
   section references this fixed name, so the recipe itself stays free of
   user-specific paths.
3. **Warn on placeholder** — when no `-d` was given, prints a warning that
   the bundled catalog is not a valid SQLite file.
4. **Run the builder** — `apptainer build --force <output> pugs.def`.
5. **Clean up** — removes `_pugs.db.staged` via an `EXIT` trap, leaving
   only `pugs.db.placeholder` in the build directory.

---

## Catalog staging

The TANGOS catalog is baked into the SIF at build time, not mounted at
runtime. Two reasons:

- HPC sites often have inconsistent paths for shared resources; baking
  the catalog in sidesteps mount-path coordination.
- The Apptainer SIF format is content-addressed, so two SIFs built from
  different catalogs are bitwise distinct — a free reproducibility check.

### Placeholder mode

If no `-d` argument is given, `pugs.db.placeholder` (a plain text file
with a notice inside) is bundled in place of the SQLite catalog. This
lets the build pipeline (and the CI workflow) succeed end-to-end before
the $2048^3$ catalog is ready.

Any TANGOS query against the bundled `pugs.db` in placeholder mode will
fail with an SQLite "file is not a database" error — by design, as a
loud failure mode.

### Replacing with the real catalog

Once the $2048^3$ TANGOS database has been built (see
[TANGOS Database](tangos-db.md)):

```bash
cd builder/container
./build_container.sh -d /path/to/pugs.db
```

The resulting `pugs.sif` will have the full halo catalog at
`/opt/PUGS/data/pugs.db`.

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

- Runtime: `tangos==1.10.0`, `pynbody==2.3.2`, `pip==25.2`
- Tooling: `camb==1.5.5` (for transfer-function regeneration),
  `pytest==9.0.2`, `requests`, plus the project linters

### Repository content bundled

The `%files` section copies a curated subset of the repo:

- `pugs/` — the `pugs` Python package
- `inputs/` — CAMB and GenetIC parameter files
- `pyproject.toml`, `README.md`
- `builder/volume_ic/` — `build.sh`, `fix_header.py`, `.checksums.txt`
- `builder/tangos_db/` — `build_tangos.sh`, `config_vars`, `tests/`

`builder/container/` is intentionally not copied (it is itself the
build context).

---

## Runtime usage

### Environment variables (pre-set inside the SIF)

| Variable | Value |
|---|---|
| `TANGOS_DB_CONNECTION` | `/opt/PUGS/data/pugs.db` |
| `TANGOS_SIMULATION_FOLDER` | `/data/sims` |
| `PUGS_HOME` | `/opt/PUGS` |
| `PYTHONUNBUFFERED` | `1` |

### Interactive shell

```bash
apptainer run pugs.sif
```

Drops into `bash` with the environment above already set.

### One-shot commands

```bash
apptainer exec pugs.sif tangos list-simulations
apptainer exec pugs.sif python -c "import tangos, pugs"
apptainer exec pugs.sif genetIC /path/to/params.txt
```

### Bind-mounting simulation data

Snapshot data lives outside the SIF — it is far too large to bundle. Mount
it at the path TANGOS expects:

```bash
apptainer exec -B /scratch/sims:/data/sims pugs.sif tangos write ...
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

- Importing `tangos`, `pynbody`, and `pugs`
- Confirming `genetIC` is on `PATH`
- Confirming the bundled catalog is present at `/opt/PUGS/data/pugs.db`
- Confirming each environment variable in the table above matches its
  expected value

The build itself runs ~12–18 minutes on `ubuntu-latest` — most of that is
compiling GenetIC and CAMB's Fortran extension.
