import zlib

import numpy as np
import pynbody
import tangos
from tangos.properties import LivePropertyCalculation, PropertyCalculation
from tangos.properties.pynbody import PynbodyPropertyCalculation
from tangos.properties.pynbody.radius import Radius


class Radius500c(Radius):
    names = "R500"

    @staticmethod
    def _get_overdensity_contrast():
        return 500

    @staticmethod
    def _get_reference_definition():
        return "critical"


class StoreIords(PynbodyPropertyCalculation):
    """
    This property stores the zlib compressed iorders of the particles within 2
    virial radii of the center of mass of the halo.
    """

    names = "zlib_ids"

    def calculate(self, particle_data, existing_properties):
        iords = np.frombuffer(
            zlib.compress(particle_data.get_index_list(particle_data.ancestor).tobytes()),
            dtype=np.int8,
        )
        return iords

    def region_specification(self, existing_properties):
        return pynbody.filt.Sphere(
            existing_properties["max_radius"] * 2, existing_properties["shrink_center"]
        )

    def requires_property(self):
        return ["shrink_center", "max_radius"]


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


class Mass200c(PropertyCalculation):
    names = "M200"

    def calculate(self, particle_data, existing_properties):
        return existing_properties["finder_mass"]

    def requires_property(self):
        return ["finder_mass"]


class MassPercentileRedshifts(PropertyCalculation):
    names = "z25_mass", "z50_mass", "z75_mass"

    @classmethod
    def preloop(self, sim, db_timestep):
        self.paths = db_timestep.calculate_all("path()")[0]

    def calculate(self, _, halo):
        m, z = tangos.get_halo(self.paths[halo.halo_number - 1]).calculate_for_progenitors(
            "M200", "z()"
        )
        return (
            z[m > 0.25 * halo["M200"]][-1],
            z[m > 0.5 * halo["M200"]][-1],
            z[m > 0.75 * halo["M200"]][-1],
        )

    def requires_property(self):
        return ["M200"]


class MergerHistory(PropertyCalculation):
    names = "N_mm", "z_lmm"

    @classmethod
    def preloop(self, sim, db_timestep):
        self.paths = db_timestep.calculate_all("path()")[0]

    def calculate(self, _, halo):
        N_mm = 0
        z_lmm = -1.0
        names, masses, z = tangos.get_halo(
            self.paths[halo.halo_number - 1]
        ).calculate_for_descendants(
            "path()",
            "M200",
            "z()",
            strategy=tangos.relation_finding.MultiHopAllProgenitorsStrategy,
        )
        all_z = np.unique(z)
        # Iterate through the redshifts in reverse
        for z_i in all_z[-1::-1]:
            this_z = z == z_i
            if (this_z).sum() < 2:
                continue
            descends = np.array([tangos.get_halo(n).next for n in names[this_z]])
            seen = set()
            merged = {x for x in descends if x in seen or seen.add(x)}
            if len(merged) == 0:
                continue
            for d in merged:
                merge_mass = masses[this_z][descends == d]
                if merge_mass.max() * 0.25 > np.sort(merge_mass)[-2]:
                    continue
                N_mm += 1
                z_lmm = z_i
        return N_mm, z_lmm

    def requires_property(self):
        return ["M200"]
