"""Hexbin grid of halo properties vs. Mhalo from the AHF halo catalog."""

import os

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

AHF_FILE = "/home/kellerbw/data/NUGS2048/DM2048.08192.z0.000.AHF_halos"
OUTPUT_PNG = os.path.join(os.path.dirname(__file__), "halo_properties.png")

# 0-indexed column positions
COL_HOST = 1
COL_MHALO = 3
COL_VMAX = 16
COL_SIGV = 18
COL_LAMBDA = 19
COL_C_SHAPE = 25
COL_EKIN = 38
COL_EPOT = 39
COL_CNFW = 42

# ---------------------------------------------------------------------------
# Load columns
# ---------------------------------------------------------------------------
usecols = (
    COL_HOST,
    COL_MHALO,
    COL_VMAX,
    COL_SIGV,
    COL_LAMBDA,
    COL_C_SHAPE,
    COL_EKIN,
    COL_EPOT,
    COL_CNFW,
)
data = np.loadtxt(AHF_FILE, comments="#", usecols=usecols)

host_mask = data[:, 0] == 0
data = data[host_mask]

mhalo = data[:, 1]  # M_sun / h
vmax = data[:, 2]  # km/s
sigv = data[:, 3]  # km/s
lam = data[:, 4]  # dimensionless spin
c_shape = data[:, 5]  # minor/major axis ratio c
ekin = data[:, 6]
epot = data[:, 7]
cnfw = data[:, 8]

# Virial parameter: 2*Ekin / |Epot|  (= 1 for virialized)
good_epot = epot < 0
virial = np.full(len(mhalo), np.nan)
virial[good_epot] = 2.0 * ekin[good_epot] / np.abs(epot[good_epot])

log10_m = np.log10(mhalo)

# ---------------------------------------------------------------------------
# Panel definitions: (y-data, y-label, log-y, y-limits)
# ---------------------------------------------------------------------------
panels = [
    (vmax, r"$V_\mathrm{max}\ [\mathrm{km\,s^{-1}}]$", True, (5, 2000)),
    (lam, r"Spin Parameter $\lambda$", True, (3e-4, 0.3)),
    (sigv, r"$\sigma_V\ [\mathrm{km\,s^{-1}}]$", True, (5, 2000)),
    (cnfw, r"Concentration $c_\mathrm{NFW}$", True, (2, 100)),
    (c_shape, r"Sphericity $c/a$", False, (0, 1.05)),
    (virial, r"$2E_\mathrm{kin}/|E_\mathrm{pot}|$", False, (0.3, 2.5)),
]

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(12, 8))
axes = axes.flatten()

for ax, (y, ylabel, log_y, ylim) in zip(axes, panels):
    # Drop NaN / non-positive values where log scale requires positivity
    finite = np.isfinite(y)
    if log_y:
        finite &= y > 0
    x_plot = log10_m[finite]
    y_plot = np.log10(y[finite]) if log_y else y[finite]
    y_lim_plot = (np.log10(ylim[0]), np.log10(ylim[1])) if log_y else ylim

    hb = ax.hexbin(
        x_plot,
        y_plot,
        gridsize=50,
        mincnt=1,
        bins="log",
        cmap="viridis",
        linewidths=0.2,
    )
    cb = fig.colorbar(hb, ax=ax, pad=0.02)

    ax.set_xlim(8, 14.5)
    ax.set_ylim(y_lim_plot)
    ax.set_xlabel(r"$\log_{10}(M_\mathrm{halo}\ [h^{-1}\,M_\odot])$", fontsize=10)
    ax.tick_params(which="both", direction="in", top=True, right=True)

    ax.set_ylabel(ylabel, fontsize=11)
    if log_y:
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"$10^{{{int(v)}}}$"))
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

fig.suptitle("NUGS2048 Halo Properties at $z=0$ (host halos only)", fontsize=13, y=1.01)
fig.tight_layout()
fig.savefig(OUTPUT_PNG, dpi=150, bbox_inches="tight")
print(f"Saved {OUTPUT_PNG}")
