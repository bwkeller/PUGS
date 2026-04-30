---
title: Cosmological Background
---

# Cosmological Background

PUGS uses the Planck 2018 LCDM best-fit parameters throughout. This page
documents those parameters, the transfer function computation, and the
physical volume properties.

---

## Planck 2018 LCDM parameters

| Symbol | CAMB name | Value | Description |
|---|---|---|---|
| $\Omega_m$ | `Om` | 0.315823 | Total matter density |
| $\Omega_\Lambda$ | `Ol` | 0.684097 | Dark energy density |
| $\Omega_b$ | `Ob` | 0 (DM-only) | Baryon density |
| $h$ | `hubble` | 0.6732117 | Hubble parameter ($H_0 = 100h$ km/s/Mpc) |
| $\sigma_8$ | `s8` | 0.8119 | Linear variance at 8 $h^{-1}$ Mpc |
| $n_s$ | `ns` | 0.9660499 | Scalar spectral index |
| $A_s$ | `scalar_amp` | $2.100549 \times 10^{-9}$ | Scalar amplitude |
| $\tau$ | `re_optical_depth` | 0.05430842 | Reionization optical depth |

These are the Plik best-fit values from [Planck 2018 Results VI](https://arxiv.org/abs/1807.06209).

---

## Transfer function

The linear matter transfer function is computed using
[CAMB](https://camb.info) from `inputs/planck_2018_CAMB.ini`.

Key CAMB settings:

```ini
get_transfer = T
do_nonlinear = 1          # HALOFIT for non-linear P(k)
transfer_kmax = 1000      # Maximum wavenumber (h/Mpc)
transfer_num_redshifts = 1
transfer_redshift(1) = 0  # Evaluate at z=0
```

The output `planck_2018_transfer_out.dat` contains 13 columns:

```
k/h [Mpc⁻¹]  CDM  baryon  photon  nu  mass_nu  total
no_nu  total_de  Weyl  v_CDM  v_b  v_b-v_c
```

GenetIC uses this file to seed the initial power spectrum at $z = 200$.

### $\sigma_8$ note

The default GenetIC transfer function (generated with `get_scalar_cls = F`)
gives $\sigma_8 = 0.8117$, slightly below the Planck best-fit value of
$0.8120$. PUGS enables the full scalar calculation, which yields
$\sigma_8 = 0.8119$ — within 0.001 of the target.

---

## Volume properties

| Property | Value |
|---|---|
| Box side length | 50 $h^{-1}$ Mpc |
| Initial redshift | $z = 200$ |
| Number of particles | $2048^3 \approx 8.59 \times 10^9$ |
| Particle mass | $\approx 8.6 \times 10^6$ M$_\odot/h$ |
| Mean inter-particle separation | $\approx 24$ $h^{-1}$ kpc |
| Force softening (typical) | $\sim 1$ $h^{-1}$ kpc |

The particle mass is:

$$m_p = \frac{\Omega_m \, \rho_\mathrm{crit} \, L^3}{N} \approx 8.6 \times 10^6 \, M_\odot/h$$

where $L = 50 \, h^{-1}$ Mpc and $N = 2048^3$.

---

## Halo finding

The PUGS TANGOS database uses the AHF (Amiga Halo Finder) halo catalog,
with halos defined by an overdensity of $\Delta = 200$ relative to the
critical density ($M_{200c}$). Only halos with at least 500 particles
($M_\mathrm{min} \approx 4 \times 10^9 \, M_\odot/h$) are included.

---

## Critical density and virial radii

The critical density at redshift $z$ is:

$$\rho_\mathrm{crit}(z) = \frac{3 H(z)^2}{8\pi G}$$

PUGS stores two virial definitions:

| Property | Overdensity | Method |
|---|---|---|
| $R_{200c}$, $M_{200c}$ | $200\,\rho_\mathrm{crit}$ | AHF finder mass |
| $R_{500c}$, $M_{500c}$ | $500\,\rho_\mathrm{crit}$ | Iterative radius solver (TANGOS `Radius` base class) |
