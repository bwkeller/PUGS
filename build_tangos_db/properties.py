import zlib

import numpy as np
import pynbody
import tangos
from tangos.properties import LivePropertyCalculation
from tangos.properties.pynbody import PynbodyPropertyCalculation
from tangos.properties.pynbody.centring import centred_calculation
from tangos.properties.pynbody.radius import Radius


class Radius500c(Radius):
    names = "r500c"

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
            existing_properties["r200c"] * 2, existing_properties["shrink_center"]
        )

    def requires_property(self):
        return ["shrink_center", "r200c"]


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


class VirialMasses(PynbodyPropertyCalculation):
    names = "m200c", "m500c"

    @centred_calculation
    def calculate(self, particle_data, existing_properties):
        from pynbody.analysis import cosmology

        rho_crit = cosmology.rho_crit(particle_data, unit="Msol kpc**-3")
        m_crit = 4 / 3 * np.pi * rho_crit
        return (
            200 * m_crit * existing_properties["r200c"] ** 3,
            500 * m_crit * existing_properties["r500c"] ** 3,
        )

    def requires_property(self):
        return ["r200c", "r500c"]


class MassPercentileRedshifts(LivePropertyCalculation):
    names = "z25_mass", "z50_mass", "z75_mass"

    def calculate(self, _, halo):
        m,z = halo.calculate_for_progenitors('m200c', 'z()')
        return z[m > 0.25*halo['m200c']][-1], \
               z[m > 0.5*halo['m200c']][-1], \
               z[m > 0.75*halo['m200c']][-1] 

    def requires_property(self):
        return ["m200c"]
