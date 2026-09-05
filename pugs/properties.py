import zlib

import numpy as np
import pynbody
from tangos.properties import LivePropertyCalculation, PropertyCalculation
from tangos.properties.pynbody import PynbodyPropertyCalculation
from tangos.properties.pynbody.centring import centred_calculation
from tangos.properties.pynbody.radius import Radius

from . import merger_forest


class Radius200c(Radius):
    names = "R200"

    @staticmethod
    def _get_overdensity_contrast():
        return 200

    @staticmethod
    def _get_reference_definition():
        return "critical"


class Radius500c(Radius):
    names = "R500"

    @staticmethod
    def _get_overdensity_contrast():
        return 500

    @staticmethod
    def _get_reference_definition():
        return "critical"


class StoreIordIndices(PynbodyPropertyCalculation):
    """
    This property stores the index position in the zlib_ids array of the
    particle at between 1 and 7 Virial radii, allowing you to select how large
    the zoom region is.  The full set is out to 8 Rvir.
    """

    names = "Rvir_indices"

    @centred_calculation
    def calculate(self, particle_data, existing_properties):
        r_values = np.arange(1, 8) * existing_properties["max_radius"]
        radii = np.sort(particle_data["r"])
        return [np.argwhere(radii > r)[0] for r in r_values]

    def region_specification(self, existing_properties):
        return pynbody.filt.Sphere(
            existing_properties["max_radius"] * 8, existing_properties["shrink_center"]
        )

    def requires_property(self):
        return ["shrink_center", "max_radius"]


class StoreIords(PynbodyPropertyCalculation):
    """
    This property stores the zlib compressed iorders of the particles within 2
    virial radii of the center of mass of the halo.
    """

    names = "zlib_ids"

    @centred_calculation
    def calculate(self, particle_data, existing_properties):
        radii = particle_data["r"]
        indices = particle_data.get_index_list(particle_data.ancestor)[np.argsort(radii)]
        iords = np.frombuffer(
            zlib.compress(indices.tobytes()),
            dtype=np.int8,
        )
        return iords

    def region_specification(self, existing_properties):
        return pynbody.filt.Sphere(
            existing_properties["max_radius"] * 8, existing_properties["shrink_center"]
        )

    def requires_property(self):
        return ["shrink_center", "max_radius"]


class VirialRatio(PropertyCalculation):
    """
    This property calculates the Virial Ratio (2T/U) of the halo.
    """

    names = "virial_ratio"

    def calculate(self, _, halo):
        return 2 * halo["Ekin"] / halo["Epot"]

    def requires_property(self):
        return ["Ekin", "Epot"]


class GetIords(LivePropertyCalculation):
    """
    This property loads and decompresses the zlib compressed iorders of the
    particles stored by the StoreIords property.
    """

    names = "ids"

    def calculate(self, _, halo):
        iords = np.frombuffer(zlib.decompress(halo["zlib_ids"].tobytes()), dtype=np.int64)
        return iords

    def requires_property(self):
        return ["zlib_ids"]


class Mass500c(PynbodyPropertyCalculation):
    names = "M500"

    def calculate(self, particle_data, existing_properties):
        from pynbody.analysis import cosmology

        rho_crit = cosmology.rho_crit(particle_data, unit="Msol kpc**-3")
        m_crit = 4 / 3 * np.pi * rho_crit
        return 500 * m_crit * existing_properties["R500"] ** 3

    def requires_property(self):
        return ["R500"]


class Mass200c(PynbodyPropertyCalculation):
    names = "M200"

    def calculate(self, particle_data, existing_properties):
        from pynbody.analysis import cosmology

        rho_crit = cosmology.rho_crit(particle_data, unit="Msol kpc**-3")
        m_crit = 4 / 3 * np.pi * rho_crit
        return 200 * m_crit * existing_properties["R200"] ** 3

    def requires_property(self):
        return ["R200"]


class MassPercentileRedshifts(LivePropertyCalculation):
    """
    Redshifts at which the halo's main branch last held 25%, 50% and 75% of its
    final mass.  The whole timestep is computed in `preloop` from the bulk
    merger forest; walking the tree per halo costs thousands of SQL round trips
    each and leaks SQLAlchemy state on every call.
    """

    names = "z25_mass", "z50_mass", "z75_mass"

    fractions = (0.25, 0.5, 0.75)

    @classmethod
    def preloop(cls, sim, db_timestep):
        forest, indices, cls.positions, cls.known = _forest_for_timestep(db_timestep)
        cls.redshifts = forest.mass_percentile_redshifts(
            np.where(indices < 0, 0, indices), cls.fractions
        )

    def calculate(self, _, halo):
        i = self.positions[halo.id]
        if not self.known[i]:
            return None, None, None
        return tuple(float(z) for z in self.redshifts[:, i])


class MergerHistory(LivePropertyCalculation):
    """
    Number of major mergers in the halo's progenitor tree, and the redshift of
    the most recent one (-1 if there was none).  As for
    `MassPercentileRedshifts`, this is computed for the whole timestep at once.
    """

    names = "N_mm", "z_lmm"

    @classmethod
    def preloop(cls, sim, db_timestep):
        forest, indices, cls.positions, cls.known = _forest_for_timestep(db_timestep)
        safe = np.where(indices < 0, 0, indices)
        cls.n_mm = forest.n_mm[safe]
        cls.z_lmm = forest.z_lmm[safe]

    def calculate(self, _, halo):
        i = self.positions[halo.id]
        if not self.known[i]:
            return None, None
        return int(self.n_mm[i]), float(self.z_lmm[i])


def _forest_for_timestep(db_timestep):
    """Bulk merger forest for a timestep's simulation, plus the lookups needed
    to go from a halo to its row in the per-timestep result arrays.

    Returns the forest; the forest index of each halo in the timestep; a map
    from database id to row; and a mask of the halos we can answer for, which
    excludes any halo the forest does not know about or has no mass for.
    """
    forest = merger_forest.get_merger_forest(db_timestep.simulation)
    ids = db_timestep.calculate_all("dbid()", object_type="halo")[0]
    positions = {int(halo_id): i for i, halo_id in enumerate(ids)}
    indices = forest.index_of(ids)
    known = (indices >= 0) & np.isfinite(forest.mass[np.where(indices < 0, 0, indices)])
    return forest, indices, positions, known
