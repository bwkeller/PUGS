---
title: pugs.ahf
---

# `pugs.ahf`

Readers for AHF halo-finder output. AHF writes one set of plain-text files per
snapshot, sharing the stem `<snapshot>.z<redshift>.AHF_`.

Halo ids are per-snapshot and, in the PUGS catalogs, run from 0. Progenitor ids
in the tree files refer to the *previous* snapshot's numbering.

---

## `read_halo_table`

```python
pugs.ahf.read_halo_table(path) -> pyarrow.Table
```

Read an `AHF_halos` file into an Arrow table with AHF's own column names, taken
from the `#ID(1)\thostHalo(2)...` header line (the `(N)` suffixes are stripped).

pyarrow's CSV reader is used rather than `np.loadtxt` because production
catalogs have millions of rows. It also types the integer columns (`ID`,
`hostHalo`, `numSubStruct`, `npart`) correctly despite AHF's whitespace
padding.

---

## `read_particle_membership`

```python
pugs.ahf.read_particle_membership(path) -> ParticleMembership
```

Read an `AHF_particles` file into compressed-sparse-row form:

```python
membership.halo_id      # sorted halo ids
membership.offsets      # halo i occupies indices[offsets[i]:offsets[i+1]]
membership.indices      # snapshot indices, concatenated
membership.particles(7) # the indices belonging to halo 7
```

**The values are indices into the snapshot, not `iord` values.** The tipsy
files carry no `iord` array, so AHF numbers particles by read order. Every
particle-derived measurement depends on this.

The file interleaves per-halo headers with particle rows, and both are two
integers wide, so blocks can only be separated by walking the counts. The bulk
read is vectorised; only the block walk is a Python loop, of one iteration per
halo rather than per particle.

---

## `read_substructure`

```python
pugs.ahf.read_substructure(path) -> dict[int, numpy.ndarray]
```

Read an `AHF_substructure` file into `{halo_id: subhalo_ids}`. Only halos that
have substructure appear in the file.

---

## `read_merger_links`

```python
pugs.ahf.read_merger_links(path) -> MergerLinks
```

Read an `AHF_croco` merger-tree file. Each entry links a halo in the current
snapshot to one in the previous snapshot:

| Field | Meaning |
|---|---|
| `descendant` | halo id in this snapshot |
| `progenitor` | halo id in the previous snapshot |
| `shared` | number of particles the two have in common |
| `merit` | AHF's figure of merit, `shared² / (n_descendant × n_progenitor)` |
| `rank` | position in AHF's ordering; `rank == 0` is its best progenitor |
| `n_descendant`, `n_progenitor` | particle counts of each |

### Why `AHF_croco` and not `AHF_mtree`

Both files describe the same edges. `AHF_croco` carries the shared-particle
counts and merits; `AHF_mtree` has them stripped out, leaving only the
progenitor ordering, from which a weight can only be invented from the rank.
Since AHF writes both side by side, there is no reason to use the weaker one.

Note that AHF's `croco` files can reference a descendant that has no row in the
matching `AHF_halos` file. Such links are dropped when the forest is built.
