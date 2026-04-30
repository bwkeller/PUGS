---
title: API Reference
---

# API Reference

The `pugs` Python package provides two public submodules.

---

## `pugs.properties`

Seven TANGOS property classes that compute halo physics from N-body particle
data. These classes are auto-discovered by TANGOS via the
`tangos.property_modules` entry point registered in `pyproject.toml`.

You do not normally call these classes directly — TANGOS invokes them when you
run `tangos write` or `tangos import-properties`. Their outputs are then
accessible via the TANGOS query interface.

See [pugs.properties](properties.md) for the full class reference.

---

## `pugs.genetic`

Helper functions for constructing zoom-in initial conditions for individual
halos. These are thin wrappers around TANGOS queries and template file
manipulation, intended to be called from a script or notebook.

See [pugs.genetic](genetic.md) for the full function reference.

---

## Quick example

```{toctree}
:hidden:

properties
genetic
```

```python
import tangos
import pugs.genetic as genetic

# Point TANGOS at the database
# (or set TANGOS_DB_CONNECTION before importing tangos)

# Get the first halo at the final snapshot
halo = tangos.get_halo("NUGS128/016/%1")

# Query precomputed properties
r200 = halo["max_radius"]          # kpc/h
m200 = halo["M200"]                # Msol/h
r500 = halo["R500"]                # kpc/h
m500 = halo["M500"]                # Msol/h
z50  = halo["z50_mass"]            # redshift of 50% mass assembly
n_mm = halo["N_mm"]                # number of major mergers
z_lmm = halo["z_lmm"]             # redshift of last major merger

# Decompress and retrieve the particle IDs
ids = halo.calculate("ids()")      # numpy int64 array

# Generate zoom-in IC files for this halo
genetic.write_particle_ids(halo.id, filename="halo1_ids.txt")
genetic.build_param_file(
    filename="genetIC_zoom.txt",
    outname="halo1_zoom",
    base_grid=2048,
    zoom_grid=[10, 2048],
)
```
