"""
Physics validation.

Every quantity is checked against an independent computation from the particle
data rather than against another stored value, so a bug in the measurement code
cannot agree with itself.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pugs import ahf, halo_properties


@pytest.fixture(scope="module")
def sample(final_snapshot):
    """The most massive few halos, by AHF's own particle count."""
    table = ahf.read_halo_table(final_snapshot.ahf("halos"))
    ids = np.asarray(table.column("ID").to_numpy())
    npart = np.asarray(table.column("npart").to_numpy())
    return ids[np.argsort(-npart)][:5]


def test_finder_mass_is_the_member_particle_mass(snapshot_data, sample):
    sim, membership, rho = snapshot_data
    for halo in sample:
        indices = membership.particles(int(halo))
        measured = halo_properties.measure_halo(sim, indices, rho)
        assert_allclose(measured.finder_mass, float(sim[indices]["mass"].sum()), rtol=1e-12)


def test_max_radius_is_the_outermost_member(snapshot_data, sample):
    """max_radius is the extent of AHF's particle set, not an overdensity radius."""
    sim, membership, rho = snapshot_data
    for halo in sample:
        indices = membership.particles(int(halo))
        measured = halo_properties.measure_halo(sim, indices, rho)
        with halo_properties.recentred(sim[indices], measured.shrink_center) as centred:
            radii = np.sqrt((np.asarray(centred["pos"], dtype=np.float64) ** 2).sum(axis=1))
        assert_allclose(measured.max_radius, radii.max(), rtol=1e-10)
        # and it really is an outer bound on the finder's particles
        assert radii.max() <= measured.max_radius * (1 + 1e-10)


def test_overdensity_radii_are_ordered(snapshot_data, sample):
    sim, membership, rho = snapshot_data
    for halo in sample:
        measured = halo_properties.measure_halo(sim, membership.particles(int(halo)), rho)
        assert measured.overdensity["R500"] < measured.overdensity["R200"]
        assert measured.overdensity["M500"] < measured.overdensity["M200"]


@pytest.mark.parametrize("radius_name,contrast", [("R200", 200), ("R500", 500)])
def test_enclosed_density_matches_the_contrast(snapshot_data, sample, radius_name, contrast):
    """
    The defining property: the mean density inside R is `contrast` times the
    critical density. Checked by counting particles, independently of the
    bisection that produced R.
    """
    import pynbody

    sim, membership, rho = snapshot_data
    for halo in sample:
        measured = halo_properties.measure_halo(sim, membership.particles(int(halo)), rho)
        radius = measured.overdensity[radius_name]
        region = sim[pynbody.filt.Sphere(radius, measured.shrink_center)]
        enclosed = float(region["mass"].sum())
        mean_density = enclosed / (4.0 / 3.0 * np.pi * radius**3)
        assert_allclose(mean_density, contrast * rho, rtol=2e-3)


@pytest.mark.parametrize("mass_name,radius_name", [("M200", "R200"), ("M500", "R500")])
def test_overdensity_mass_matches_enclosed_particles(snapshot_data, sample, mass_name, radius_name):
    import pynbody

    sim, membership, rho = snapshot_data
    for halo in sample:
        measured = halo_properties.measure_halo(sim, membership.particles(int(halo)), rho)
        region = sim[pynbody.filt.Sphere(measured.overdensity[radius_name], measured.shrink_center)]
        assert_allclose(measured.overdensity[mass_name], float(region["mass"].sum()), rtol=2e-3)


def test_centre_lies_inside_the_box(snapshot_data, final_snapshot):
    sim, membership, rho = snapshot_data
    box = halo_properties._boxsize(sim)
    assert box is not None, "boxsize must resolve, or centres are never wrapped"
    for halo in membership.halo_id[:40]:
        centre, _ = halo_properties.centre_and_extent(sim[membership.particles(int(halo))])
        assert np.all(np.abs(centre) <= box / 2 + 1e-6)


def test_halo_straddling_the_boundary_is_centred_correctly(snapshot_data):
    """
    A halo spanning the periodic boundary has particles at both ends of the
    box. Without wrapping, the shrinking sphere converges to empty space
    between the two pieces and the overdensity bisection then fails outright.
    """
    sim, membership, rho = snapshot_data
    box = halo_properties._boxsize(sim)

    straddling = []
    for halo in membership.halo_id:
        indices = membership.particles(int(halo))
        if len(indices) < 500:
            continue
        span = np.ptp(np.asarray(sim[indices]["pos"], dtype=np.float64), axis=0)
        if span.max() > 0.5 * box:
            straddling.append(int(halo))
    assert straddling, "expected at least one halo across the boundary in NUGS128"

    for halo in straddling:
        measured = halo_properties.measure_halo(sim, membership.particles(int(halo)), rho)
        # the halo must come out compact, not box-sized
        assert measured.max_radius < 0.1 * box
        assert measured.overdensity["R500"] < measured.overdensity["R200"] < 0.1 * box
        assert np.all(np.abs(measured.shrink_center) <= box / 2 + 1e-6)
        # R200 tracks the finder's extent closely but is not bounded by it:
        # it is measured over every particle in the region, while max_radius
        # only covers the particles AHF left bound to the halo.
        assert 0.5 < measured.overdensity["R200"] / measured.max_radius < 2.0


def test_empty_halo_is_rejected(snapshot_data):
    sim, _, _ = snapshot_data
    with pytest.raises(ValueError):
        halo_properties.centre_and_extent(sim[np.empty(0, dtype=np.int64)])


def test_recentred_restores_positions(snapshot_data):
    """
    Positions are shared with the parent snapshot, so a leaked shift would
    corrupt every halo measured afterwards.
    """
    sim, membership, _ = snapshot_data
    indices = membership.particles(int(membership.halo_id[0]))
    before = np.array(sim[indices]["pos"])
    with halo_properties.recentred(sim[indices], np.array([10.0, 20.0, 30.0])):
        pass
    assert_allclose(np.asarray(sim[indices]["pos"]), before, rtol=0, atol=0)
