"""Compare AHF halo catalog HMF to the Fernandez-Garcia et al. 2026 fit."""

import os

import matplotlib.pyplot as plt
import numpy as np

# colossus 1.3.x uses np.trapz which was removed in NumPy 2.0
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid

from colossus.cosmology import cosmology
from colossus.lss import mass_function

AHF_FILE = "/home/kellerbw/data/NUGS2048/DM2048.08192.z0.000.AHF_halos"
BOX_SIZE_MPCH = 50.0  # h^-1 Mpc
Z_SNAP = 0.0
MDEF = "200m"  # AHF ovdens ~200x mean background density
N_BINS = 25
OUTPUT_PNG = os.path.join(os.path.dirname(__file__), "hmf_comparison.png")

cosmology.setCosmology("planck18")

# ---------------------------------------------------------------------------
# Load halo catalog
# ---------------------------------------------------------------------------
# Columns 1 = hostHalo, 3 = Mhalo (0-indexed)
data = np.loadtxt(AHF_FILE, comments="#", usecols=(1, 3))
host_mask = data[:, 0] == 0
hosts = data[host_mask, 1]  # M_sun / h
log10_m = np.log10(hosts)

# ---------------------------------------------------------------------------
# Measured HMF: dn / d log10 M  [(Mpc/h)^-3]
# ---------------------------------------------------------------------------
volume = BOX_SIZE_MPCH**3  # (Mpc/h)^3

m_min = log10_m.min()
m_max = log10_m.max()
bin_edges = np.linspace(m_min, m_max, N_BINS + 1)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
dlog10m = bin_edges[1] - bin_edges[0]

counts, _ = np.histogram(log10_m, bins=bin_edges)
dn_dlog10m = counts / (volume * dlog10m)

# Only plot bins with at least one halo
mask = counts > 0
bin_centers = bin_centers[mask]
dn_dlog10m = dn_dlog10m[mask]
counts = counts[mask]

# Poisson error bars: sigma_N / (V * dlog10m)
err = np.sqrt(counts) / (volume * dlog10m)

# ---------------------------------------------------------------------------
# Theoretical HMF from Fernandez-Garcia et al. 2026
# ---------------------------------------------------------------------------
m_theory = np.logspace(8, 15.5, 300)  # M_sun / h
mf_theory = mass_function.massFunction(
    m_theory, Z_SNAP, mdef=MDEF, model="fernandezgarcia26", q_out="dndlnM"
)
# dndlnM -> dn/dlog10M: multiply by ln(10)
mf_theory_log10 = mf_theory * np.log(10)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))

ax.errorbar(
    10**bin_centers,
    dn_dlog10m,
    yerr=err,
    fmt="o",
    ms=4,
    color="steelblue",
    elinewidth=1.2,
    capsize=3,
    label="Simulation",
    zorder=5,
)

ax.plot(
    m_theory,
    mf_theory_log10,
    color="tomato",
    lw=2,
    label="Fernandez-Garcia et al. 2026",
)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$M\ [h^{-1}\,M_\odot]$", fontsize=13)
ax.set_ylabel(r"$dn/d\log_{10}M\ [(h^{-1}\,\mathrm{Mpc})^{-3}]$", fontsize=13)
ax.set_title(
    r"NUGS2048 Halo Mass Function at $z=0$",
    fontsize=12,
)
ax.legend(fontsize=11)
ax.set_xlim(10**8, 10**15)
ax.set_ylim(1e-6, 1e2)
ax.tick_params(which="both", direction="in", top=True, right=True)

fig.tight_layout()
fig.savefig(OUTPUT_PNG, dpi=150)
print(f"Saved {OUTPUT_PNG}")
