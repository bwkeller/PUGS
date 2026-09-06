---
title: The Pipeline
---

# The Pipeline

PUGS processes cosmological data in sequential stages. Each stage has its own
builder directory and can be run independently once its prerequisites are met.

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

## Stage 2 — Halo catalog

**Directory:** `builder/catalog/`

Takes the N-body snapshots produced by running the IC through your chosen
N-body code, together with the AHF halo-finder output beside them, and builds
a Parquet halo catalog. There is no intermediate database.

```
N-body snapshots + AHF output
        │
        ├── AHF_halos      ──► 43 finder columns, verbatim
        ├── AHF_particles  ──► halo membership (snapshot indices)
        ├── AHF_croco      ──► merger tree with shared-particle merits
        │                          │
        │                          ▼
        │                   pugs.merger_forest
        │                          │
        │              N_mm, z_lmm, z25/z50/z75_mass,
        │              main progenitor / descendant pointers
        │
        └── snapshot particles ──► pugs.halo_properties
                                       │
                        shrink_center, max_radius, R200/R500, M200/M500
                                       │
                                       ▼
                        halos_<snapshot>.parquet, one per snapshot
                             + provenance.json
```

**Key outputs:** a catalog directory — one Parquet file per snapshot plus a
`provenance.json` sidecar, readable without installing PUGS.

See [Halo Catalog](catalog.md) for the full description.

---

## Stage 3 — Containerization

**Directory:** `builder/container/`

Bundles the Python stack, GenetIC, and the Parquet halo catalog into a single
Apptainer/Singularity image so the pipeline can be moved to any HPC site
without rebuilding the dependency tree.

```
builder/container/pugs.def
            +
   pugs/, inputs/, builder/* ...
            +
     the catalog directory
              │
              ▼
    build_container.sh -c /path/to/catalog
              │
              ▼
        pugs.sif  (portable SIF)
```

**Key output:** `pugs.sif` — the portable, content-addressed container image.

See [Container Builder](container.md) for the full description.

---

## Zoom-in ICs

After the catalog is built, individual halos can be used to generate zoom-in
initial conditions for higher-resolution re-simulations. The particle ids are
read from the AHF membership, or selected from the snapshot within a multiple
of the halo's radius:

```
halo_id + snapshots  ──► pugs.genetic.write_particle_ids()  ──► id_file.txt
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
catalog
container
```
