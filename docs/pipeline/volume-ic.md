---
title: Volume Initial Conditions
---

# Stage 1: Volume Initial Conditions

**Directory:** `builder/volume_ic/`

This stage produces the large collisionless dark matter volume that all PUGS
zoom simulations are drawn from.

---

## Overview

| Property | Value |
|---|---|
| Box size | 50 $h^{-1}$ Mpc |
| Particle count (production) | $2048^3 \approx 8.6 \times 10^9$ |
| Initial redshift | $z = 200$ |
| Cosmology | Planck 2018 LCDM |
| Output format | tipsy |
| RAM required (production) | ≥ 320 GB |
| Wall time (production, 96 cores) | ~90 minutes |

---

## Build script

```{code-block} bash
./build.sh         # Production: 2048³, 96 OMP threads
./build.sh -t      # Test mode: 128³, 4 OMP threads
./build.sh -h      # Print help
```

### What the script does

1. **Create `outputs/`** — all outputs are written here.

2. **Run CAMB**
   ```bash
   camb ../../inputs/planck_2018_CAMB.ini
   ```
   Produces `planck_2018_transfer_out.dat` — the linear matter transfer
   function evaluated at $z = 0$.

3. **Prepare the GenetIC parameter file** — `inputs/genetIC_volume.txt` is
   copied into the working directory. In test mode every occurrence of `2048`
   is replaced with `128`.

4. **Run GenetIC via Docker**
   ```bash
   docker run -e OMP_NUM_THREADS=$NUM_THREADS --rm \
       -v `pwd`:/w/ --user $(id -u):$(id -g) \
       apontzen/genetic:1.5.0 /w/genetIC_volume.txt
   ```
   GenetIC reads the transfer function and produces `DM_volume.tipsy`.

5. **Fix 32-bit header (production only)** — calls `./fix_header.py` (see
   below).

6. **Move outputs** — `planck_2018_transfer_out.dat`, `tipsy.param`, and
   `IC_output.params` are relocated to `outputs/`.

---

## GenetIC parameters (`inputs/genetIC_volume.txt`)

The key physical parameters are:

| Parameter | Value | Description |
|---|---|---|
| `Om` | 0.315823 | Matter density $\Omega_m$ |
| `Ol` | 0.684097 | Dark energy density $\Omega_\Lambda$ |
| `Ob` | 0 | Baryon density (DM-only run) |
| `s8` | 0.8119 | $\sigma_8$ normalization |
| `ns` | 0.9660499 | Scalar spectral index $n_s$ |
| `hubble` | 0.6732117 | Hubble parameter $h$ |
| `z` | 200 | Starting redshift |
| `random_seed` | 8896131 | Fixed seed for reproducibility |
| `base_grid` | 50 Mpc/h, 2048 | Box size and grid resolution |

The transfer function is loaded from `planck_2018_transfer_out.dat` (generated
by CAMB in the previous step).

---

## The 32-bit header fix (`fix_header.py`)

The tipsy file format stores particle counts in a 32-bit unsigned integer.
$2048^3 = 4{,}294{,}967{,}296 = 2^{32}$, which is exactly one bit beyond the
32-bit range.

`fix_header.py` resolves the overflow by splitting the high bit into the
tipsy header's `pad` field:

```python
n = 2048**3
pad = ((n >> 32) & 0xFF) + ((n >> 16) & 0xFF0000)
n &= 0xFFFFFFFF
```

The header is then rewritten as:
`(time, n_truncated, ndim=3, ng=0, nd=n_truncated, ns=0, pad=overflow_bits)`

This scheme is specific to how TANGOS parses the tipsy format; it is not
needed for the $128^3$ test volume.

---

## Reproducibility

The IC is bitwise reproducible **only when run with the same number of OMP
threads**. GenetIC's FFT-based convolutions are thread-count dependent. To
match the checked-in checksums, build with exactly 96 threads (or use the
test mode with 4 threads for the 128³ checksums in `.checksums.txt`).

Verify the test build:

```bash
md5sum -c .checksums.txt
```

Expected checksums (128³ test volume):

| File | MD5 |
|---|---|
| `outputs/DM_volume.tipsy` | `993da5b41eb021a6026ea03055536aa8` |
| `outputs/planck_2018_transfer_out.dat` | `0c6c7d44132b7ab6ab5c2e2145648a5b` |
| `outputs/tipsy.param` | `5783d988ce4380a5eb7577b08011c53a` |
| `outputs/camb.log` | `8e83f6ea9ba24f383dd6ccdbf7ad6d19` |

---

## $\sigma_8$ note

The default GenetIC transfer function produces $\sigma_8 = 0.8117$, slightly
below the Planck 2018 Plik best-fit value of $0.8120$. This discrepancy arises
because the GenetIC default disables scalar CMB spectra (`get_scalar_cls = F`).
PUGS uses a CAMB parameter file that enables the full scalar calculation,
yielding the correct $\sigma_8 = 0.8119$.
