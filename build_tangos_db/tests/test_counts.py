import tangos


def test_sim_count():
    """
    Make sure that the correct number of simulations are loaded.
    """
    assert len(tangos.all_simulations()) == 1


def test_snap_count():
    """
    Make sure that the correct number of snapshots are loaded.
    """
    assert len(tangos.get_simulation("NUGS128").timesteps) == 16


def test_halo_count():
    """
    Make sure we have the correct number of halos.
    """
    assert tangos.get_simulation("NUGS128").timesteps[-1].halos.count() == 219
