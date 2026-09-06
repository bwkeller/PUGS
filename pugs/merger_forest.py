"""
Merger trees built directly from AHF's own tree output.

For each snapshot AHF writes an ``AHF_croco`` file listing, for every halo, the
halos in the previous snapshot that contributed particles to it, together with
the number of shared particles and AHF's merit function
``shared**2 / (n_descendant * n_progenitor)``.  That is everything needed to
build the tree; no database and no re-derivation from particle lists.

Which links to keep
-------------------
A link is kept when at least :data:`MIN_SHARED_FRACTION` of the *progenitor's*
particles ended up in the descendant -- that is, when most of the progenitor
really did become part of that halo.

Thresholding on the merit instead would be wrong for merger counting, because
the merit penalises a size mismatch: a small halo that is swallowed whole by a
much larger one has a merit near zero while having contributed every one of its
particles.  In the NUGS128 z=0 catalog, halo 0's second progenitor has merit
0.0087 and yet gave up 227 of its 227 particles.  A merit cut would discard it;
it is a real merger.

Requiring *strictly* more than half also makes the pruned graph a forest by
construction -- a progenitor cannot give more than half of itself to two
different descendants -- which is what lets the merger counts below be
accumulated in a single sweep.  The comparison has to be strict: with a
non-strict one, a progenitor split exactly evenly between two descendants would
be linked to both and the guarantee would fail.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import ahf
from .simulation import Simulation, halo_id

#: A link is kept when strictly more than this fraction of the progenitor's
#: particles end up in the descendant.  Must be at least 0.5, and the
#: comparison must be strict, for the pruned graph to be a forest.
MIN_SHARED_FRACTION = 0.5

#: A merger is "major" when the second most massive progenitor is at least this
#: fraction of the most massive one.
MAJOR_MERGER_RATIO = 0.25

#: Column of the AHF halo table used for merger ratios and assembly history.
#: Masses only ever enter as ratios, so AHF's h-scaled mass is fine and costs
#: no particle work.
MASS_COLUMN = "Mhalo"


@dataclass
class MergerForest:
    """The merger forest of a simulation, in flat numpy arrays.

    Halos are addressed by *index* into these arrays.  ``halo_id`` is sorted,
    so :meth:`index_of` maps a catalog halo id to its row.
    """

    halo_id: np.ndarray
    snapshot_index: np.ndarray
    finder_id: np.ndarray
    mass: np.ndarray
    npart: np.ndarray
    redshift: np.ndarray
    #: index of the best progenitor of each halo, -1 if none
    main_progenitor: np.ndarray
    #: index of the halo each halo merges into, -1 if none
    descendant: np.ndarray
    #: number of progenitors above the threshold
    n_progenitors: np.ndarray
    n_mm: np.ndarray
    z_lmm: np.ndarray

    def __len__(self) -> int:
        return len(self.halo_id)

    def index_of(self, values) -> np.ndarray:
        """Map catalog halo ids to indices into this forest, -1 if absent."""
        values = np.asarray(values, dtype=np.int64)
        position = np.searchsorted(self.halo_id, values)
        np.clip(position, 0, max(len(self.halo_id) - 1, 0), out=position)
        if len(self.halo_id) == 0:
            return np.full(values.shape, -1, dtype=np.int64)
        return np.where(self.halo_id[position] == values, position, -1)

    def mass_percentile_redshifts(self, indices, fractions) -> np.ndarray:
        """Redshifts at which each halo's main branch last held ``fractions``
        of its final mass.

        Walking backwards in time, redshift only increases, so the answer is
        the highest redshift at which the main branch was still above the
        threshold.  Returned as ``(len(fractions), len(indices))``.
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
            redshift = self.redshift[self.snapshot_index[here]]
            mass = self.mass[here]
            for row, fraction in enumerate(fractions):
                above = mass > fraction * final_mass[live]
                out[row, live_at[above]] = redshift[above]
            current[live] = self.main_progenitor[here]

        return out


def build_forest(
    simulation: Simulation,
    halo_tables: dict[str, "object"] | None = None,
    min_shared_fraction: float = MIN_SHARED_FRACTION,
    min_particles: int = 0,
) -> MergerForest:
    """Build the merger forest of ``simulation`` from its AHF tree files.

    By default the forest includes *every* halo AHF found, not only those above
    the catalog's particle cut.  Merger counts and assembly histories are
    properties of the tree, and dropping small progenitors before building it
    biases them: a halo whose two progenitors both fall below the cut looks
    like it had no merger at all.  ``min_particles`` is provided to reproduce
    that narrower tree for comparison, not because it is the better choice.
    """
    halo_id_parts, snapshot_parts, finder_parts = [], [], []
    mass_parts, npart_parts = [], []
    finder_by_snapshot: list[np.ndarray] = []

    for snapshot in simulation:
        table = (halo_tables or {}).get(snapshot.extension)
        if table is None:
            table = ahf.read_halo_table(snapshot.ahf("halos"))
        finder = np.asarray(table.column("ID").to_numpy(), dtype=np.int64)
        counts = np.asarray(table.column("npart").to_numpy(), dtype=np.int64)
        order = np.argsort(finder, kind="stable")
        if min_particles:
            order = order[counts[order] >= min_particles]
        finder = finder[order]

        finder_by_snapshot.append(finder)
        finder_parts.append(finder)
        halo_id_parts.append(halo_id(snapshot.index, finder))
        snapshot_parts.append(np.full(finder.size, snapshot.index, dtype=np.int64))
        mass_parts.append(np.asarray(table.column(MASS_COLUMN).to_numpy(), dtype=np.float64)[order])
        npart_parts.append(np.asarray(table.column("npart").to_numpy(), dtype=np.int64)[order])

    forest = MergerForest(
        halo_id=_concat(halo_id_parts, np.int64),
        snapshot_index=_concat(snapshot_parts, np.int64),
        finder_id=_concat(finder_parts, np.int64),
        mass=_concat(mass_parts, np.float64),
        npart=_concat(npart_parts, np.int64),
        redshift=np.array([snapshot.redshift for snapshot in simulation], dtype=np.float64),
        main_progenitor=np.empty(0, dtype=np.int64),
        descendant=np.empty(0, dtype=np.int64),
        n_progenitors=np.empty(0, dtype=np.int64),
        n_mm=np.empty(0, dtype=np.int64),
        z_lmm=np.empty(0, dtype=np.float64),
    )

    offsets = np.concatenate([[0], np.cumsum([f.size for f in finder_by_snapshot])]).astype(
        np.int64
    )
    _link_forest(forest, simulation, finder_by_snapshot, offsets, min_shared_fraction)
    _accumulate_merger_history(forest)
    return forest


def _concat(parts, dtype) -> np.ndarray:
    return np.concatenate(parts).astype(dtype) if parts else np.empty(0, dtype=dtype)


def _link_forest(forest, simulation, finder_by_snapshot, offsets, min_shared_fraction) -> None:
    """Fill in ``main_progenitor``, ``descendant`` and ``n_progenitors``."""
    n_halos = len(forest.halo_id)
    forest.main_progenitor = np.full(n_halos, -1, dtype=np.int64)
    forest.descendant = np.full(n_halos, -1, dtype=np.int64)
    forest.n_progenitors = np.zeros(n_halos, dtype=np.int64)

    # a progenitor keeps the descendant it gave the largest fraction of itself to
    best_fraction = np.zeros(n_halos, dtype=np.float64)
    best_merit = np.zeros(n_halos, dtype=np.float64)

    for snapshot in simulation:
        if snapshot.index == 0 or not snapshot.has_ahf("croco"):
            continue
        links = ahf.read_merger_links(snapshot.ahf("croco"))
        if len(links) == 0:
            continue

        fraction = np.where(
            links.n_progenitor > 0, links.shared / np.maximum(links.n_progenitor, 1), 0.0
        )
        keep = fraction > min_shared_fraction
        if not keep.any():
            continue

        descendant_idx = _locate(
            links.descendant[keep], finder_by_snapshot[snapshot.index], offsets[snapshot.index]
        )
        progenitor_idx = _locate(
            links.progenitor[keep],
            finder_by_snapshot[snapshot.index - 1],
            offsets[snapshot.index - 1],
        )
        merit = links.merit[keep]
        fraction = fraction[keep]

        valid = (descendant_idx >= 0) & (progenitor_idx >= 0)
        descendant_idx, progenitor_idx = descendant_idx[valid], progenitor_idx[valid]
        merit, fraction = merit[valid], fraction[valid]
        if descendant_idx.size == 0:
            continue

        np.add.at(forest.n_progenitors, descendant_idx, 1)

        # main progenitor of each descendant: highest merit, AHF's own ranking
        _assign_best(forest.main_progenitor, best_merit, descendant_idx, progenitor_idx, merit)
        # descendant of each progenitor: wherever most of it ended up
        _assign_best(forest.descendant, best_fraction, progenitor_idx, descendant_idx, fraction)


def _assign_best(target, best_score, keys, values, scores) -> None:
    """For each key, keep the value with the highest score."""
    order = np.lexsort((-scores, keys))
    keys, values, scores = keys[order], values[order], scores[order]
    starts = np.flatnonzero(np.r_[True, keys[1:] != keys[:-1]])
    winners, winning_score, winning_value = keys[starts], scores[starts], values[starts]
    better = winning_score > best_score[winners]
    target[winners[better]] = winning_value[better]
    best_score[winners[better]] = winning_score[better]


def _locate(finder_ids, sorted_finder, offset) -> np.ndarray:
    """Global forest indices for AHF ids within one snapshot, -1 if unknown."""
    if sorted_finder.size == 0:
        return np.full(len(finder_ids), -1, dtype=np.int64)
    position = np.searchsorted(sorted_finder, finder_ids)
    np.clip(position, 0, sorted_finder.size - 1, out=position)
    return np.where(sorted_finder[position] == finder_ids, position + offset, -1)


def _accumulate_merger_history(forest: MergerForest) -> None:
    """Count major mergers over each halo's progenitor tree.

    The threshold on shared particles makes the links a forest, so a halo's
    progenitor set is exactly its subtree and the counts accumulate in one
    sweep from the earliest snapshot to the latest.
    """
    n_halos = len(forest.halo_id)
    forest.n_mm = np.zeros(n_halos, dtype=np.int64)
    forest.z_lmm = np.full(n_halos, -1.0, dtype=np.float64)

    order = np.argsort(forest.snapshot_index, kind="stable")
    bounds = np.searchsorted(forest.snapshot_index[order], np.arange(forest.redshift.size + 1))

    for level in range(forest.redshift.size):
        members = order[bounds[level] : bounds[level + 1]]
        progenitors = members[forest.descendant[members] >= 0]
        if progenitors.size == 0:
            continue

        descendants = forest.descendant[progenitors]
        # group by descendant, heaviest progenitor first within each group
        grouped = np.lexsort((-forest.mass[progenitors], descendants))
        progenitors, descendants = progenitors[grouped], descendants[grouped]

        starts = np.flatnonzero(np.r_[True, descendants[1:] != descendants[:-1]])
        counts = np.diff(np.r_[starts, descendants.size])
        group = descendants[starts]

        mass = forest.mass[progenitors]
        event = np.zeros(starts.size, dtype=bool)
        multiple = counts >= 2
        event[multiple] = mass[starts[multiple] + 1] >= MAJOR_MERGER_RATIO * mass[starts[multiple]]

        forest.n_mm[group] += np.add.reduceat(forest.n_mm[progenitors], starts) + event

        inherited = np.minimum.reduceat(_no_merger_to_inf(forest.z_lmm[progenitors]), starts)
        here = np.where(event, forest.redshift[level], np.inf)
        merged = np.minimum(np.minimum(inherited, here), _no_merger_to_inf(forest.z_lmm[group]))
        forest.z_lmm[group] = np.where(np.isfinite(merged), merged, -1.0)


def _no_merger_to_inf(z: np.ndarray) -> np.ndarray:
    """Map the "no merger" sentinel of -1 onto +inf so it loses any minimum."""
    return np.where(z < 0, np.inf, z)
