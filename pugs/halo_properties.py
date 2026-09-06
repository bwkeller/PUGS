"""
Per-halo physics measured directly from the snapshot.

These are the quantities TANGOS used to compute through its property-calculator
framework, written here as plain functions over a pynbody snapshot and a set of
particle indices.  Nothing here needs a database: give it a snapshot and the
AHF membership for a halo and it returns numbers.

Two conventions worth stating, because the names invite confusion:

``max_radius``
    The distance from the halo centre to the *outermost halo-finder particle*.
    It is a measure of the extent of AHF's particle set, not an overdensity
    radius.  It is close to AHF's own ``Rhalo`` but not equal to it, and it is
    the radius the particle-ID regions are defined in multiples of.

``R200`` / ``R500``
    Genuine spherical-overdensity radii, found by bisection on the enclosed
    density using every particle in the region, not only the halo's own.

Periodicity is handled in two different ways, deliberately.  Centring wraps
explicitly, because ``shrink_sphere_center`` knows nothing about the box and a
halo straddling a boundary would otherwise have its particles spread across the
full box width and its centre land in the empty space between the two pieces --
in the NUGS128 test data the third most massive halo does exactly this.  Region
selection instead relies on ``pynbody.filt.Sphere`` being periodic-aware, which
avoids wrapping the whole snapshot once per halo.
"""

from __future__ import annotations

import contextlib
import logging
from dataclasses import dataclass, field

import numpy as np

logger = logging.getLogger(__name__)

#: Radius searched for the overdensity bisection, in units of ``max_radius``.
#: The overdensity radius of a halo lies well inside its finder particle set,
#: so three times its extent is a generous bracket.
SEARCH_FACTOR = 3.0

#: Shrink factor for the shrinking-sphere centre.
SHRINK_FACTOR = 0.8

#: Overdensity contrasts measured for every halo, as ``(radius, mass, contrast)``.
OVERDENSITIES = (("R200", "M200", 200), ("R500", "M500", 500))


@dataclass
class HaloMeasurement:
    """Everything measured from particle data for a single halo."""

    finder_mass: float
    max_radius: float
    shrink_center: np.ndarray
    overdensity: dict[str, float] = field(default_factory=dict)

    def as_dict(self) -> dict[str, float]:
        centre = np.asarray(self.shrink_center, dtype=np.float64)
        values = {
            "finder_mass": self.finder_mass,
            "max_radius": self.max_radius,
            "shrink_center_x": float(centre[0]),
            "shrink_center_y": float(centre[1]),
            "shrink_center_z": float(centre[2]),
        }
        values.update(self.overdensity)
        return values


def _boxsize(sim) -> float | None:
    """Box side length in the snapshot's position units, or None if unknown.

    The header stores a comoving, h-scaled length, so converting it to the
    position units needs the snapshot's conversion context; without that
    pynbody raises "Not convertible" rather than returning anything useful.
    """
    boxsize = sim.properties.get("boxsize")
    if boxsize is None:
        return None
    positions = sim["pos"]
    try:
        return float(boxsize.in_units(positions.units, **positions.conversion_context()))
    except Exception as exc:
        # Not fatal -- centring still wraps via pynbody -- but the returned
        # centre will not be folded back into the box, so say so rather than
        # silently disabling it.
        logger.warning("could not determine boxsize (%s); centres will not be wrapped", exc)
        return None


def wrap_position(position: np.ndarray, boxsize: float | None) -> np.ndarray:
    """Fold a position back into the box, which is centred on the origin."""
    if not boxsize:
        return position
    return ((position + boxsize / 2.0) % boxsize) - boxsize / 2.0


@contextlib.contextmanager
def recentred(view, centre):
    """Temporarily move ``view`` to ``centre`` and wrap it around the box.

    Positions are restored on exit; the array is shared with the parent
    snapshot, so leaving it shifted would corrupt every later halo.
    """
    original = np.array(view["pos"])
    try:
        view["pos"] -= centre
        view.wrap()
        yield view
    finally:
        view["pos"] = original


def centre_and_extent(particles) -> tuple[np.ndarray, float]:
    """Shrinking-sphere centre of a halo, and its outermost particle radius.

    The halo is temporarily recentred on one of its own particles and wrapped
    first.  Without the wrap, a halo straddling the periodic boundary has its
    particles spread across the full box width and the centre lands in the
    empty space between the two pieces.
    """
    import pynbody

    if len(particles) == 0:
        raise ValueError("cannot centre a halo with no particles")

    boxsize = _boxsize(particles.ancestor)
    anchor = np.array(particles["pos"][0], dtype=np.float64)

    with recentred(particles, anchor):
        centre = pynbody.analysis.halo.shrink_sphere_center(
            particles, shrink_factor=SHRINK_FACTOR, particles_for_velocity=0
        )
        centre = np.asarray(centre, dtype=np.float64)
        if not np.all(np.isfinite(centre)):
            raise RuntimeError("shrinking-sphere centre did not converge to a finite position")

        offset = np.asarray(particles["pos"], dtype=np.float64) - centre
        max_radius = float(np.sqrt((offset**2).sum(axis=1)).max())

    return wrap_position(centre + anchor, boxsize), max_radius


def overdensity_radius(region, contrast: int, rho_def: str = "critical") -> float:
    """Radius enclosing ``contrast`` times the reference density.

    ``region`` must already be centred on the halo.  Recentring before the call
    rather than passing a centre to pynbody matters: pynbody derives its
    bisection bracket from the coordinate range of the region, so leaving the
    halo at its raw box coordinates changes the bracket and shifts the root.
    """
    import pynbody

    return float(pynbody.analysis.halo.virial_radius(region, overden=contrast, rho_def=rho_def))


def overdensity_mass(radius: float, contrast: int, rho_crit: float) -> float:
    """Mass implied by an overdensity radius.

    Exact by the definition of the radius -- ``M(<R) / (4/3 pi R^3)`` is
    ``contrast * rho_crit`` -- so this agrees with summing the enclosed particle
    masses to within the bisection tolerance, without a second particle pass.
    """
    return contrast * (4.0 / 3.0) * np.pi * rho_crit * radius**3


def critical_density(sim) -> float:
    """Critical density of the snapshot in Msol kpc^-3."""
    from pynbody.analysis import cosmology

    return float(cosmology.rho_crit(sim, unit="Msol kpc**-3"))


def measure_halo(
    sim,
    indices: np.ndarray,
    rho_crit: float,
    search_factor: float = SEARCH_FACTOR,
    overdensities=OVERDENSITIES,
) -> HaloMeasurement:
    """Measure one halo, given the snapshot indices AHF assigned to it.

    ``sim`` must be in physical units.  ``rho_crit`` is passed in rather than
    recomputed because it depends only on the snapshot.
    """
    import pynbody

    particles = sim[indices]
    centre, max_radius = centre_and_extent(particles)

    measurement = HaloMeasurement(
        finder_mass=float(particles["mass"].sum()),
        max_radius=max_radius,
        shrink_center=centre,
    )

    # Sphere is periodic-aware, so this picks up the far side of a halo that
    # straddles a boundary without wrapping the whole snapshot.
    region = sim[pynbody.filt.Sphere(search_factor * max_radius, centre)]
    with recentred(region, centre):
        for radius_name, mass_name, contrast in overdensities:
            radius = overdensity_radius(region, contrast)
            measurement.overdensity[radius_name] = radius
            measurement.overdensity[mass_name] = overdensity_mass(radius, contrast, rho_crit)

    return measurement
