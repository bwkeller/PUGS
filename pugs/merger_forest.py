"""
Bulk merger-tree analysis for a whole TANGOS simulation.

The per-halo alternative -- calling ``calculate_for_progenitors`` /
``calculate_for_descendants`` once for every halo -- costs O(N_progenitors)
individual SQL round trips per halo, and every ``MultiHopStrategy`` execution
also leaks SQLAlchemy schema and mapper state (a mapper, two dependency
processors on the ``SimulationObjectBase`` mapper, and a few thousand tracked
Python objects per call).  On a 461k-halo z=0 snapshot that makes the run both
extremely slow and progressively slower.

This module instead loads the whole merger forest once, with a handful of
table-scale queries, and derives every quantity with linear numpy passes.
"""

import warnings

import numpy as np
from tangos import core

#: Dictionary items under which merger-tree links are stored.  PUGS databases
#: use "ahf_tree_link", written by
#: ``tangos.tools.ahf_merger_tree_importer.AHFTreeImporter.create_links``, which
#: stores each link in both directions with the same weight; "ptcls_in_common"
#: is what TANGOS' other tree importers use.  Whichever exist are accepted, so
#: that links from unrelated relations (e.g. crosslink's "match") are ignored.
LINK_RELATIONS = ("ahf_tree_link", "ptcls_in_common")

#: Minimum shared-particle fraction for a link to count as a tree edge.  This
#: mirrors ``min_onehop_reverse_weight`` in
#: ``tangos.relation_finding.MultiHopAllProgenitorsStrategy`` (which
#: ``MultiHopMajorProgenitorsStrategy`` inherits), i.e. the threshold TANGOS
#: applies when walking progenitors.  Because the AHF importer stores symmetric
#: weights, the forward and reverse thresholds coincide.
MIN_LINK_WEIGHT = 0.1

#: A merger is "major" when the second most massive progenitor is at least this
#: fraction of the most massive one.
MAJOR_MERGER_RATIO = 0.25

_FETCH_CHUNK = 1_000_000


class MergerForest:
    """The merger forest of a single simulation, held in numpy arrays.

    Halos are referred to by *index*: ``halo_id`` is sorted, so
    :meth:`index_of` maps a TANGOS database id to its position in every other
    array here.
    """

    def __init__(self, db_simulation):
        self._load_halos(db_simulation)
        self._load_links(db_simulation)
        self._accumulate_merger_history()

    # ------------------------------------------------------------------
    # loading
    # ------------------------------------------------------------------

    def _load_halos(self, db_simulation):
        timesteps = list(db_simulation.timesteps)  # ordered by time_gyr
        self.redshift = np.array([ts.redshift for ts in timesteps], dtype=np.float64)

        ids, masses, numbers, ts_index = [], [], [], []
        for i, ts in enumerate(timesteps):
            ts_ids, ts_mass, ts_number = ts.calculate_all(
                "dbid()", "Mhalo", "halo_number()", object_type="halo", sanitize=False
            )
            ts_ids = np.asarray(ts_ids, dtype=np.int64)
            ids.append(ts_ids)
            masses.append(_as_float(ts_mass))
            numbers.append(np.asarray(ts_number, dtype=np.int64))
            ts_index.append(np.full(ts_ids.size, i, dtype=np.int32))

        halo_id = np.concatenate(ids)
        order = np.argsort(halo_id, kind="stable")

        self.halo_id = halo_id[order]
        self.mass = np.concatenate(masses)[order]
        self.halo_number = np.concatenate(numbers)[order]
        self.ts_index = np.concatenate(ts_index)[order]

        # halos grouped by timestep, so each level is a contiguous slice
        self._by_ts = np.argsort(self.ts_index, kind="stable")
        self._ts_bounds = np.searchsorted(self.ts_index[self._by_ts], np.arange(len(timesteps) + 1))

    def _load_links(self, db_simulation):
        """Load the merger-tree links and reduce them to major descendant and
        major progenitor pointers."""
        session = core.get_default_session()
        relations = [
            rid
            for rid in (
                core.dictionary.get_dict_id(name, None, session=session) for name in LINK_RELATIONS
            )
            if rid is not None
        ]
        if not relations:
            raise RuntimeError(
                "No merger-tree links found; expected one of %s. "
                "Has `tangos import-ahf-trees` been run?" % (LINK_RELATIONS,)
            )
        timestep_ids = [ts.id for ts in db_simulation.timesteps]

        link = core.halo_data.HaloLink.__table__
        halo = core.halo.SimulationObjectBase.__table__
        stmt = (
            link.select()
            .with_only_columns(link.c.halo_from_id, link.c.halo_to_id, link.c.weight)
            .select_from(link.join(halo, halo.c.id == link.c.halo_from_id))
            .where(link.c.relation_id.in_(relations))
            .where(halo.c.timestep_id.in_(timestep_ids))
            .where(link.c.weight > MIN_LINK_WEIGHT)
        )

        from_idx, to_idx, weight = self._fetch_links(session, stmt)

        t_from = self.ts_index[from_idx]
        t_to = self.ts_index[to_idx]

        # A link pointing forwards in time is a candidate descendant link, one
        # pointing backwards a candidate progenitor link.  Links within a
        # timestep (there should be none in an AHF tree) are ignored.
        forwards = t_to > t_from
        backwards = t_to < t_from

        # major descendant: earliest later timestep, then heaviest link, then
        # lowest halo number -- matching HopMajorDescendantStrategy's ordering
        self.next_idx = self._best_neighbour(
            from_idx[forwards], to_idx[forwards], weight[forwards], t_to[forwards]
        )
        # major progenitor: latest earlier timestep, same tie-breaks -- matching
        # MultiHopMajorProgenitorsStrategy
        self.main_prog_idx = self._best_neighbour(
            from_idx[backwards],
            to_idx[backwards],
            weight[backwards],
            -t_to[backwards].astype(np.int64),
        )

        self._check_is_forest(from_idx[forwards])

    @staticmethod
    def _check_is_forest(link_sources):
        """Warn if a halo has more than one descendant above the threshold.

        Aggregating over the tree in a single sweep relies on the pruned links
        forming a forest, which is what makes a halo's progenitor set exactly
        its subtree here.  AHF trees satisfy this -- in the NUGS2048 database
        all 29,984,548 halos with a tree link above the threshold have exactly
        one -- but if a future catalogue does not, the counts would start to
        drift from a per-halo `MultiHopAllProgenitorsStrategy` walk, so say so
        rather than failing quietly.
        """
        if link_sources.size == 0:
            return
        _, counts = np.unique(link_sources, return_counts=True)
        branching = int((counts > 1).sum())
        if branching:
            warnings.warn(
                "%d halos have more than one descendant linked above a shared-particle "
                "fraction of %g. The merger tree is not a forest, so merger counts may "
                "differ from a per-halo progenitor walk." % (branching, MIN_LINK_WEIGHT),
                RuntimeWarning,
            )

    def _fetch_links(self, session, stmt):
        """Stream the link query into index arrays, dropping links that point
        outside the halos we loaded."""
        conn = session.connection().execution_options(stream_results=True)
        result = conn.execute(stmt)

        from_chunks, to_chunks, weight_chunks = [], [], []
        while True:
            rows = result.fetchmany(_FETCH_CHUNK)
            if not rows:
                break
            n = len(rows)
            from_chunks.append(np.fromiter((r[0] for r in rows), dtype=np.int64, count=n))
            to_chunks.append(np.fromiter((r[1] for r in rows), dtype=np.int64, count=n))
            weight_chunks.append(np.fromiter((r[2] for r in rows), dtype=np.float64, count=n))

        if not from_chunks:
            empty_i = np.empty(0, dtype=np.int64)
            return empty_i, empty_i, np.empty(0, dtype=np.float64)

        from_id = np.concatenate(from_chunks)
        to_id = np.concatenate(to_chunks)
        weight = np.concatenate(weight_chunks)

        from_idx = self.index_of(from_id)
        to_idx = self.index_of(to_id)
        keep = (from_idx >= 0) & (to_idx >= 0)
        return from_idx[keep], to_idx[keep], weight[keep]

    def _best_neighbour(self, source, target, weight, time_key):
        """For each ``source``, pick the ``target`` minimising ``time_key``,
        then maximising ``weight``, then minimising halo number.

        Returns a full-length array of target indices, -1 where a source has no
        candidate.
        """
        best = np.full(self.halo_id.size, -1, dtype=np.int64)
        if source.size == 0:
            return best

        # np.lexsort applies the last key first, so this sorts by source, then
        # time_key ascending, then weight descending, then halo number ascending
        order = np.lexsort((self.halo_number[target], -weight, time_key, source))
        source_sorted = source[order]
        starts = np.flatnonzero(np.r_[True, source_sorted[1:] != source_sorted[:-1]])
        best[source_sorted[starts]] = target[order][starts]
        return best

    # ------------------------------------------------------------------
    # derived quantities
    # ------------------------------------------------------------------

    def _accumulate_merger_history(self):
        """Count major mergers over each halo's progenitor tree.

        Every halo has at most one major descendant, so the links form a
        forest and the counts can be accumulated in a single sweep from the
        earliest timestep to the latest.  A halo's progenitors are grouped by
        their own timestep, matching the per-redshift grouping of the original
        per-halo implementation.
        """
        n_halos = self.halo_id.size
        self.n_mm = np.zeros(n_halos, dtype=np.int64)
        self.z_lmm = np.full(n_halos, -1.0, dtype=np.float64)

        for level in range(self.redshift.size):
            members = self._by_ts[self._ts_bounds[level] : self._ts_bounds[level + 1]]
            progenitors = members[self.next_idx[members] >= 0]
            if progenitors.size == 0:
                continue

            descendants = self.next_idx[progenitors]
            # group by descendant, heaviest progenitor first within each group
            order = np.lexsort((-self.mass[progenitors], descendants))
            progenitors = progenitors[order]
            descendants = descendants[order]

            starts = np.flatnonzero(np.r_[True, descendants[1:] != descendants[:-1]])
            counts = np.diff(np.r_[starts, descendants.size])
            group = descendants[starts]

            # a major merger needs at least two progenitors, with the second
            # most massive at least MAJOR_MERGER_RATIO of the most massive
            event = np.zeros(starts.size, dtype=bool)
            multiple = counts >= 2
            mass = self.mass[progenitors]
            event[multiple] = (
                mass[starts[multiple] + 1] >= MAJOR_MERGER_RATIO * mass[starts[multiple]]
            )

            self.n_mm[group] += np.add.reduceat(self.n_mm[progenitors], starts) + event

            # z_lmm is the lowest redshift of any major merger in the subtree
            inherited = np.minimum.reduceat(_no_merger_to_inf(self.z_lmm[progenitors]), starts)
            here = np.where(event, self.redshift[level], np.inf)
            merged = np.minimum(np.minimum(inherited, here), _no_merger_to_inf(self.z_lmm[group]))
            self.z_lmm[group] = np.where(np.isfinite(merged), merged, -1.0)

    def mass_percentile_redshifts(self, indices, fractions):
        """Redshifts at which each halo's main branch last exceeded each of
        ``fractions`` times its final mass.

        This reproduces ``z[m > f * M][-1]`` over the results of
        ``calculate_for_progenitors``, whose rows run from low to high
        redshift: the answer is the *highest* redshift at which the main branch
        was above the threshold.
        """
        indices = np.asarray(indices, dtype=np.int64)
        current = indices.copy()
        final_mass = self.mass[indices]
        out = np.zeros((len(fractions), indices.size), dtype=np.float64)

        while True:
            live = current >= 0
            if not live.any():
                break
            live_at = np.flatnonzero(live)
            here = current[live]
            redshift = self.redshift[self.ts_index[here]]
            mass = self.mass[here]
            # walking backwards in time means redshift only increases, so each
            # assignment overwrites with a higher redshift
            for row, fraction in enumerate(fractions):
                above = mass > fraction * final_mass[live]
                out[row, live_at[above]] = redshift[above]
            current[live] = self.main_prog_idx[here]

        return out

    # ------------------------------------------------------------------

    def index_of(self, halo_id):
        """Map TANGOS database ids to indices into this forest, -1 if absent."""
        halo_id = np.asarray(halo_id, dtype=np.int64)
        idx = np.searchsorted(self.halo_id, halo_id)
        np.clip(idx, 0, self.halo_id.size - 1, out=idx)
        return np.where(self.halo_id[idx] == halo_id, idx, -1)


def _as_float(values):
    """Convert a possibly object-dtype array containing None into float64."""
    values = np.asarray(values, dtype=object)
    return np.array([np.nan if v is None else float(v) for v in values], dtype=np.float64)


def _no_merger_to_inf(z):
    """Map the "no merger" sentinel of -1 onto +inf so it loses any minimum."""
    return np.where(z < 0, np.inf, z)


_cache = {}


def get_merger_forest(db_simulation):
    """Return the (cached) :class:`MergerForest` for a simulation.

    Only one forest is kept at a time; it is several GB for a 2048^3 run.
    """
    forest = _cache.get(db_simulation.id)
    if forest is None:
        forest = MergerForest(db_simulation)
        _cache.clear()
        _cache[db_simulation.id] = forest
    return forest
