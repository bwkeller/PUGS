import numpy as np
import pynbody as pyn
import pytest
import tangos
from numpy.testing import assert_allclose, assert_array_equal, assert_array_less


@pytest.fixture(scope="function")
def last_halo():
    snap = tangos.get_simulation("NUGS128").timesteps[-1]
    sim = pyn.load(snap.filename)
    sim.physical_units()
    h = sim.halos()
    return sim, h[0], snap.halos[0]


def test_property_counts():
    for ts in tangos.get_simulation("NUGS128").timesteps:
        halo_count = ts.halos.count()
        assert ts.calculate_all("shrink_center")[0].size == halo_count * 3  # x,y,z
        assert ts.calculate_all("max_radius")[0].size == halo_count
        assert ts.calculate_all("finder_mass")[0].size == halo_count
        assert ts.calculate_all("M200")[0].size == halo_count
        assert ts.calculate_all("M500")[0].size == halo_count
        assert ts.calculate_all("R500")[0].size == halo_count
        if ts == tangos.get_simulation("NUGS128").timesteps[-1]:
            assert ts.calculate_all("zlib_ids")[0].size == halo_count
            assert ts.calculate_all("N_mm")[0].size == halo_count
            assert ts.calculate_all("z_lmm")[0].size == halo_count
            assert ts.calculate_all("z25_mass")[0].size == halo_count
            assert ts.calculate_all("z50_mass")[0].size == halo_count
            assert ts.calculate_all("z75_mass")[0].size == halo_count


def test_id_compression(last_halo):
    """
    Ensure that the compressed IDs decompress to what is obtained from pynbody
    """
    psim, ph, th = last_halo
    tangos_ids = th.calculate("ids()")
    snap_ids = ph.get_index_list(psim)
    assert np.all(np.isin(snap_ids, tangos_ids))


def test_rvir_idx(last_halo):
    """
    Ensure that the index of particles at 1,2,...7 Rvir are correct.
    """
    psim, ph, th = last_halo
    indices = th["Rvir_indices"]
    tangos_ids = th.calculate("ids()")
    assert len(indices) == 7
    rvir_ids = tangos_ids[indices]
    with psim.translate(-th["shrink_center"]):
        particle_r = psim["r"]
        for i in range(7):
            assert particle_r[rvir_ids[i]] > (i + 1) * th["max_radius"] * pyn.units.kpc
            assert particle_r[rvir_ids[i]] < (i + 2) * th["max_radius"] * pyn.units.kpc


def test_virial_radii(last_halo):
    """
    Ensure that the virial radii R200c and R500c are correct.
    """
    psim, ph, th = last_halo
    filt = pyn.filt.Sphere(th["max_radius"] * 3, th["shrink_center"])
    snap_r500c = pyn.analysis.halo.virial_radius(
        psim[filt], cen=th["shrink_center"], overden=500, rho_def="critical"
    )
    tangos_r200c = th["max_radius"]
    tangos_r500c = th["R500"]
    assert tangos_r500c < tangos_r200c
    assert_array_equal(tangos_r500c, snap_r500c)


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
    tangos_m200c = th["M200"]
    tangos_m500c = th["M500"]
    snap_m200c = psim[pyn.filt.Sphere(snap_r200c, th["shrink_center"])]["mass"].sum()
    snap_m500c = psim[pyn.filt.Sphere(snap_r500c, th["shrink_center"])]["mass"].sum()
    assert tangos_m500c < tangos_m200c
    assert_allclose(tangos_m200c, snap_m200c, rtol=1e-3)
    assert_allclose(tangos_m500c, snap_m500c, rtol=1e-3)


def test_mass_percentiles():
    """
    Check that the halo mass assembly percentiles are sane.
    """
    ts = tangos.get_simulation("NUGS128").timesteps[-1]
    z25_mass, z50_mass, z75_mass = ts.calculate_all("z25_mass", "z50_mass", "z75_mass")
    assert np.all(z25_mass >= z50_mass)
    assert np.all(z50_mass >= z75_mass)
    assert_array_less(0, z25_mass)
    assert_array_less(0, z50_mass)
    assert_array_less(0, z75_mass)


def test_merger_history():
    """
    Check that the merger history quantities are sane
    """
    ts = tangos.get_simulation("NUGS128").timesteps[-1]
    z_lmm, N_mm = ts.calculate_all("z_lmm", "N_mm")
    assert z_lmm.min() >= -1.0
    assert N_mm.min() >= 0
    assert_array_less(0, z_lmm[N_mm > 0])
    assert_array_equal(-1.0, z_lmm[N_mm == 0])
    assert_array_less(0, N_mm[N_mm != 0])
