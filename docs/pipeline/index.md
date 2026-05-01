---
title: The Two-Stage Pipeline
---

# The Two-Stage Pipeline

PUGS processes cosmological data in two sequential stages. Each stage has its
own builder directory and can be run independently once its prerequisites are
met.

---

## Stage 1 — Volume initial conditions

**Directory:** `builder/volume_ic/`

Takes Planck 2018 cosmological parameters and produces a single large
tipsy-format initial-condition file ready for an N-body code.

```
inputs/planck_2018_CAMB.ini
        │
        ▼
      CAMB  ──► planck_2018_transfer_out.dat
                        │
        inputs/genetIC_volume.txt
                        │
                        ▼
                    GenetIC (Docker)  ──► DM_volume.tipsy
                                                │
                                         fix_header.py  (production only)
                                                │
                                                ▼
                                      outputs/DM_volume.tipsy
```

**Key outputs:** `outputs/DM_volume.tipsy` — the $2048^3$-particle,
50 $h^{-1}$ Mpc IC file.

See [Volume IC](volume-ic.md) for the full description.

---

## Stage 2 — TANGOS halo catalog

**Directory:** `builder/tangos_db/`

Takes the N-body snapshots produced by running the IC through your chosen
N-body code, and builds a fully-featured TANGOS database containing halos
with precomputed physical properties.

```
N-body snapshots
        │
        ▼
  tangos add  ──► halos registered
        │
        ▼
  tangos import-properties
        │   (pugs.properties entry point)
        ▼
  tangos import-ahf-trees  ──► merger trees built
        │
        ▼
  tangos write (per-halo properties)
        │
        ▼
  pugs.db  (SQLite TANGOS database)
```

**Key outputs:** `pugs.db` — the SQLite database of halo properties.

See [TANGOS Database](tangos-db.md) for the full description.

---

## Stage 3 — Containerization

**Directory:** `builder/container/`

Bundles the Python stack, GenetIC, and the TANGOS halo catalog into a single
Apptainer/Singularity image so the pipeline can be moved to any HPC site
without rebuilding the dependency tree.

```
builder/container/pugs.def
            +
   pugs/, inputs/, builder/* ...
            +
        pugs.db
              │
              ▼
    build_container.sh -d pugs.db
              │
              ▼
        pugs.sif  (portable SIF)
```

**Key output:** `pugs.sif` — the portable, content-addressed container image.

See [Container Builder](container.md) for the full description.

---

## Zoom-in ICs

After the TANGOS database is built, individual halos can be used to generate
zoom-in initial conditions for higher-resolution re-simulations:

```
pugs.db  ──► pugs.genetic.write_particle_ids()  ──► id_file.txt
                                                          │
          inputs/zoom_template.txt                        │
                     │                                    │
                     ▼                                    │
       pugs.genetic.build_param_file()                    │
                     │                                    │
                     ▼                                    ▼
              genetIC_zoom.txt  ──►  GenetIC  ──►  zoom IC
```

See [pugs.genetic](../api/genetic.md) for the API.

```{toctree}
:hidden:

volume-ic
tangos-db
container
```
