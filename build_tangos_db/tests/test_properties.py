import contextlib

import numpy as np
import pynbody as pyn
import pytest
import tangos
from numpy.testing import assert_allclose, assert_equal


@pytest.fixture(scope="function")
def last_halo():
    snap = tangos.get_simulation("NUGS128").timesteps[-1]
    sim = pyn.load(snap.filename)
    sim.physical_units()
    h = sim.halos()
    return sim, h[0], snap.halos[0]


def test_id_compression(last_halo):
    """
    Ensure that the compressed IDs decompress to what is obtained from pynbody
    """
    psim, ph, th = last_halo
    tangos_ids = th.calculate("ids()")
    snap_ids = ph.get_index_list(psim)
    assert np.all(np.isin(snap_ids, tangos_ids))


def test_virial_radii(last_halo):
    """
    Ensure that the virial radii R200c and R500c are correct.
    """
    psim, ph, th = last_halo
    filt = pyn.filt.Sphere(th["max_radius"] * 3, th["shrink_center"])
    snap_r200c = pyn.analysis.halo.virial_radius(
        psim[filt], cen=th["shrink_center"], overden=200, rho_def="critical"
    )
    snap_r500c = pyn.analysis.halo.virial_radius(
        psim[filt], cen=th["shrink_center"], overden=500, rho_def="critical"
    )
    tangos_r200c = th["r200c"]
    tangos_r500c = th["r500c"]
    assert tangos_r500c < tangos_r200c
    assert_equal(tangos_r200c, snap_r200c)
    assert_equal(tangos_r500c, snap_r500c)


def test_virial_mass(last_halo):
    """
    Ensure that the virial masses R200c and R500c are correct.
    """
    psim, ph, th = last_halo
    filt = pyn.filt.Sphere(th["max_radius"] * 3, th["shrink_center"])
    snap_r200c = pyn.analysis.halo.virial_radius(
        psim[filt], cen=th["shrink_center"], overden=200, rho_def="critical"
    )
    snap_r500c = pyn.analysis.halo.virial_radius(
        psim[filt], cen=th["shrink_center"], overden=500, rho_def="critical"
    )
    tangos_m200c = th["m200c"]
    tangos_m500c = th["m500c"]
    snap_m200c = psim[pyn.filt.Sphere(snap_r200c, th["shrink_center"])]["mass"].sum()
    snap_m500c = psim[pyn.filt.Sphere(snap_r500c, th["shrink_center"])]["mass"].sum()
    assert tangos_m500c < tangos_m200c
    assert_allclose(tangos_m200c, snap_m200c, rtol=1e-3)
    assert_allclose(tangos_m500c, snap_m500c, rtol=1e-3)


def test_mass_percentiles(last_halo):
    """
    Check that the halo mass assembly percentiles are sane.
    """
    psim, ph, th = last_halo
    assert th.calculate("z25_mass()") > th.calculate("z50_mass()")
    assert th.calculate("z50_mass()") > th.calculate("z75_mass()")
