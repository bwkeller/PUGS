---
title: pugs.properties
---

# `pugs.properties`

TANGOS property calculator classes. All classes are registered as entry points
under `tangos.property_modules` and are discovered automatically by TANGOS.

---

## Inheritance overview

```
tangos.properties.PropertyCalculation
├── Mass200c
├── MassPercentileRedshifts
└── MergerHistory

tangos.properties.pynbody.PynbodyPropertyCalculation
├── Mass500c
├── StoreIords
└── StoreIordIndices

tangos.properties.pynbody.radius.Radius
└── Radius500c

tangos.properties.LivePropertyCalculation
└── GetIords
```

---

## Virial radius and mass

### `Radius500c`

**TANGOS name:** `R500`

Computes the radius enclosing a mean density of 500 times the critical
density.

Inherits from the TANGOS `Radius` base class and overrides two methods:

```python
class Radius500c(Radius):
    names = "R500"

    @staticmethod
    def _get_overdensity_contrast():
        return 500

    @staticmethod
    def _get_reference_definition():
        return "critical"
```

**Units:** kpc/$h$ (pynbody physical units)

---

### `Mass500c`

**TANGOS name:** `M500`

Computes the total mass within $R_{500c}$ analytically from the critical
density and `R500`:

$$M_{500c} = 500 \cdot \rho_\mathrm{crit} \cdot \frac{4}{3}\pi \, R_{500c}^3$$

**Prerequisites:** `R500`

**Units:** M$_\odot$/$h$

---

### `Mass200c`

**TANGOS name:** `M200`

Retrieves the halo finder mass (`finder_mass`), which is defined as the mass
within $R_{200c}$ in the AHF halo finder.

```python
def calculate(self, particle_data, existing_properties):
    return existing_properties["finder_mass"]
```

**Prerequisites:** `finder_mass`

**Units:** M$_\odot$/$h$

---

## Particle ID storage

### `StoreIords`

**TANGOS name:** `zlib_ids`

Stores the iorders (global particle indices) of all particles within
8 $R_\mathrm{vir}$ of the halo center, sorted by radius, as a
zlib-compressed byte array.

**Calculation:**
1. Select all particles within a sphere of radius $8 \times R_\mathrm{vir}$
   centered on `shrink_center`.
2. Sort particle indices by distance from center.
3. Compress the sorted index array with `zlib.compress`.
4. Return as a `numpy.int8` array (the raw compressed bytes).

**Prerequisites:** `shrink_center`, `max_radius`

**Written:** z=0 only (final snapshot)

**Storage size:** typically 5–15 MB per halo after compression.

---

### `GetIords`

**TANGOS name:** `ids`

Live property (not stored in the database) that decompresses `zlib_ids` on
demand:

```python
def calculate(self, _, halo):
    return np.frombuffer(
        zlib.decompress(halo["zlib_ids"].tobytes()), dtype=np.int64
    )
```

Call via `halo.calculate("ids()")` or `tangos.get_halo(...).calculate("ids()")`.

**Prerequisites:** `zlib_ids`

---

### `StoreIordIndices`

**TANGOS name:** `Rvir_indices`

Stores 7 integer indices marking the positions in the sorted `zlib_ids` array
where the cumulative particle list crosses $n \times R_\mathrm{vir}$ for
$n = 1, 2, \ldots, 7$.

**Purpose:** Enables $O(1)$ extraction of particles within any integer
multiple of $R_\mathrm{vir}$ without decompressing and re-sorting the full
ID array:

```python
ids_within_3rvir = all_ids[:halo["Rvir_indices"][2]]
```

**Calculation:**
```python
r_values = np.arange(1, 8) * existing_properties["max_radius"]
radii = np.sort(particle_data["r"])
return [np.argwhere(radii > r)[0] for r in r_values]
```

**Prerequisites:** `shrink_center`, `max_radius`

**Written:** z=0 only (final snapshot)

---

## Assembly history

### `MassPercentileRedshifts`

**TANGOS names:** `z25_mass`, `z50_mass`, `z75_mass`

Computes the redshifts at which the halo first assembled 25%, 50%, and 75%
of its present-day $M_{200c}$.

**Algorithm:**

1. Walk the merger tree backwards via `calculate_for_progenitors("M200", "z()")`.
2. At each progenitor snapshot, sum the progenitor masses.
3. Find the highest redshift at which the running mass exceeds each threshold.

```python
m, z = halo.calculate_for_progenitors("M200", "z()")
z25 = z[m > 0.25 * halo["M200"]][-1]
z50 = z[m > 0.50 * halo["M200"]][-1]
z75 = z[m > 0.75 * halo["M200"]][-1]
```

**Prerequisites:** `M200`

**Written:** z=0 only

**Invariant:** `z25_mass ≥ z50_mass ≥ z75_mass > 0`

---

### `MergerHistory`

**TANGOS names:** `N_mm` (int), `z_lmm` (float)

Counts major mergers and records the redshift of the most recent one.

**Definition:** A merger is *major* when the mass ratio of the two merging
halos satisfies

$$\frac{M_\mathrm{minor}}{M_\mathrm{major}} \geq \frac{1}{4}$$

**Algorithm:**

1. Retrieve the full progenitor tree with `MultiHopAllProgenitorsStrategy`.
2. For each redshift slice (from high $z$ to low $z$):
   - Identify halos that share the same descendant at the next timestep
     (i.e., halos that merge).
   - For each merging group, compare the largest and second-largest
     progenitor masses.
   - If the ratio exceeds 1:4, increment `N_mm` and update `z_lmm`.

**Return values:**

| Value | Meaning |
|---|---|
| `N_mm ≥ 1`, `z_lmm > 0` | Halo has experienced major mergers |
| `N_mm = 0`, `z_lmm = −1.0` | No major mergers in the tree |

**Prerequisites:** `M200`

**Written:** z=0 only
