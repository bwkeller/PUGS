"""
Tests for the bulk merger-forest calculation that backs N_mm, z_lmm and the
mass-assembly redshifts.

These build their own small databases rather than using NUGS128, so the merger
tree can be shaped deliberately.
"""

import gc
import os

import numpy as np
import pytest
import sqlalchemy
import tangos
from tangos import core, relation_finding

from pugs import merger_forest


def _new_db(path):
    core.init_db("sqlite:///%s" % path)
    from tangos.testing.simulation_generator import SimulationGeneratorForTests

    generator = SimulationGeneratorForTests("sim", max_steps=64)
    session = core.get_default_session()
    relation = core.dictionary.get_or_create_dictionary_item(
        session, merger_forest.LINK_RELATIONS[0]
    )
    session.commit()
    return generator, session, relation


def _link(session, links, early, late, relation, weight):
    """Store a link in both directions, as the AHF tree importer does."""
    links.append(core.halo_data.HaloLink(early, late, relation, weight))
    links.append(core.halo_data.HaloLink(late, early, relation, weight))


@pytest.fixture(scope="module")
def restore_db():
    """Put the default connection back, so the NUGS128 tests are unaffected."""
    previous = os.environ.get("TANGOS_DB_CONNECTION")
    yield
    core.close_db()
    if previous is not None:
        core.init_db(previous)


@pytest.fixture(scope="module")
def synthetic_sim(tmp_path_factory, restore_db):
    """A randomly wired merger tree with the structure of an AHF tree: each
    halo has exactly one descendant, linked above the pruning threshold."""
    generator, session, relation = _new_db(
        tmp_path_factory.mktemp("merger_forest") / "synthetic.db"
    )
    n_steps, n_halos = 12, 30
    rng = np.random.default_rng(12345)

    steps = []
    for _ in range(n_steps):
        generator.add_timestep()
        generator.add_objects_to_timestep(
            n_halos, NDM=[int(x) for x in rng.integers(500, 200000, n_halos)]
        )
        generator.add_properties_to_halos(Mhalo=lambda _i: float(rng.uniform(1e10, 1e14)))
        steps.append(generator.sim.timesteps[-1])

    links = []
    for early_step, late_step in zip(steps[:-1], steps[1:]):
        late = late_step.objects.all()
        for halo in early_step.objects.all():
            target = late[int(rng.integers(len(late)))]
            _link(session, links, halo, target, relation, float(rng.uniform(0.11, 1.0)))
    session.add_all(links)
    session.commit()

    return tangos.get_simulation("sim")


def _reference_merger_history(paths, halo_number):
    """The original per-halo implementation, kept here as the reference."""
    n_mm = 0
    z_lmm = -1.0
    names, masses, z = tangos.get_halo(paths[halo_number - 1]).calculate_for_descendants(
        "path()",
        "Mhalo",
        "z()",
        strategy=relation_finding.MultiHopAllProgenitorsStrategy,
    )
    for z_i in np.unique(z)[::-1]:
        this_z = z == z_i
        if this_z.sum() < 2:
            continue
        descends = np.array([tangos.get_halo(n).next for n in names[this_z]])
        seen = set()
        merged = {x for x in descends if x in seen or seen.add(x)}
        for d in merged:
            merge_mass = masses[this_z][descends == d]
            if merge_mass.max() * 0.25 > np.sort(merge_mass)[-2]:
                continue
            n_mm += 1
            z_lmm = z_i
    return n_mm, z_lmm


def _reference_percentiles(paths, halo_number, final_mass):
    m, z = tangos.get_halo(paths[halo_number - 1]).calculate_for_progenitors("Mhalo", "z()")
    return (
        z[m > 0.25 * final_mass][-1],
        z[m > 0.5 * final_mass][-1],
        z[m > 0.75 * final_mass][-1],
    )


def test_merger_history_matches_reference(synthetic_sim):
    ts = synthetic_sim.timesteps[-1]
    forest = merger_forest.MergerForest(synthetic_sim)
    ids, numbers = ts.calculate_all("dbid()", "halo_number()", object_type="halo")
    paths = ts.calculate_all("path()", object_type="halo")[0]
    indices = forest.index_of(ids)
    assert (indices >= 0).all()

    for k, halo_number in enumerate(numbers):
        expected = _reference_merger_history(paths, int(halo_number))
        got = (int(forest.n_mm[indices[k]]), float(forest.z_lmm[indices[k]]))
        assert got == pytest.approx(expected), f"halo {halo_number}"


def test_mass_percentiles_match_reference(synthetic_sim):
    ts = synthetic_sim.timesteps[-1]
    forest = merger_forest.MergerForest(synthetic_sim)
    ids, numbers, masses = ts.calculate_all("dbid()", "halo_number()", "Mhalo", object_type="halo")
    paths = ts.calculate_all("path()", object_type="halo")[0]
    indices = forest.index_of(ids)
    got = forest.mass_percentile_redshifts(indices, (0.25, 0.5, 0.75))

    for k, halo_number in enumerate(numbers):
        expected = _reference_percentiles(paths, int(halo_number), masses[k])
        assert got[:, k] == pytest.approx(expected), f"halo {halo_number}"


def test_multihop_queries_do_not_accumulate_mapper_state(synthetic_sim):
    """Each multi-hop query builds a throwaway ORM class; tangos must let go of
    it again, or a per-halo tree walk leaks without bound."""
    halo = synthetic_sim.timesteps[-1].halos[0]
    halo_mapper = sqlalchemy.inspect(core.halo.SimulationObjectBase)

    def walk():
        halo.calculate_for_descendants(
            "Mhalo", strategy=relation_finding.MultiHopAllProgenitorsStrategy
        )

    walk()
    gc.collect()
    before = (
        len(core.Base.registry._managers),
        len(core.Base.registry._class_registry),
        len(halo_mapper._dependency_processors),
    )

    for _ in range(20):
        walk()
    gc.collect()
    after = (
        len(core.Base.registry._managers),
        len(core.Base.registry._class_registry),
        len(halo_mapper._dependency_processors),
    )

    assert after == before


@pytest.fixture(scope="module")
def shadowed_route_sim(tmp_path_factory, restore_db):
    """The smallest tree in which TANGOS' multi-hop walk loses a progenitor.

    ``A`` is a genuine progenitor of ``R``, reachable via ``D``.  It is also
    reachable via ``E``, on a route with a higher aggregated weight but whose
    reverse link is under the 0.1 threshold.  MultiHopStrategy thins each hop
    down to the strongest route *before* applying that threshold, so the ``E``
    route displaces the ``D`` route and is then rejected, and ``A`` is dropped.
    """
    generator, session, relation = _new_db(
        tmp_path_factory.mktemp("merger_forest_shadow") / "shadow.db"
    )

    generator.add_timestep()
    (a,) = generator.add_objects_to_timestep(1, NDM=[10000])
    generator.add_properties_to_halos(Mhalo=lambda _i: 1e12)

    generator.add_timestep()
    d, e = generator.add_objects_to_timestep(2, NDM=[20000, 20000])
    generator.add_properties_to_halos(Mhalo=lambda _i: 2e12)

    generator.add_timestep()
    (r,) = generator.add_objects_to_timestep(1, NDM=[40000])
    generator.add_properties_to_halos(Mhalo=lambda _i: 4e12)

    links = []
    _link(session, links, a, d, relation, 0.5)  # A's only above-threshold descendant
    _link(session, links, a, e, relation, 0.09)  # below threshold, but a stronger route
    _link(session, links, d, r, relation, 0.11)
    _link(session, links, e, r, relation, 0.9)
    session.add_all(links)
    session.commit()

    return tangos.get_simulation("sim"), {"a": a.id, "d": d.id, "e": e.id, "r": r.id}


def test_forest_keeps_progenitor_that_multihop_drops(shadowed_route_sim):
    sim, ids = shadowed_route_sim
    root = tangos.get_halo(ids["r"])

    walked = root.calculate_for_descendants(
        "dbid()", strategy=relation_finding.MultiHopAllProgenitorsStrategy
    )[0]
    assert set(walked) == {ids["r"], ids["d"], ids["e"]}, "expected TANGOS to drop A"

    forest = merger_forest.MergerForest(sim)
    index = forest.index_of(np.array([ids["a"], ids["d"], ids["e"], ids["r"]]))
    a_i, d_i, e_i, r_i = (int(i) for i in index)

    # A is in the tree, below D, so R's progenitors are A, D and E
    assert int(forest.next_idx[a_i]) == d_i
    assert int(forest.next_idx[d_i]) == r_i
    assert int(forest.next_idx[e_i]) == r_i

    # D and E merge into R and are equal in mass, so that is a major merger;
    # the walk above sees it too. The point is that A is not silently lost.
    assert forest.n_mm[r_i] == 1
