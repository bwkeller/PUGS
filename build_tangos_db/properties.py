from tangos.properties.pynbody import PynbodyPropertyCalculation
from tangos.properties.pynbody.profile import HaloDensityProfile
from tangos.properties.pynbody.BH import BHAccHistogram
from tangos.properties.pynbody.centring import centred_calculation
from tangos.properties import PropertyCalculation, LivePropertyCalculation, TimeChunkedProperty
from tangos.live_calculation import NoResultsError
import pynbody
import numpy as np
import zlib
import scipy

class StoreIords(PynbodyPropertyCalculation):
    names = "zlib_ids"

    def calculate(self, particle_data, existing_properties):
        iords =  np.frombuffer(zlib.compress(particle_data.get_index_list(particle_data.ancestor).tobytes()), dtype=np.int8)
        return iords

    def region_specification(self, existing_properties):
        return pynbody.filt.Sphere(existing_properties['max_radius']*2,
                                   existing_properties['shrink_center'])

    def requires_property(self):
            return ["shrink_center", "max_radius"]

class GetIords(LivePropertyCalculation):
    names = "ids"

    def calculate(self, particle_data, existing_properties):
        iords = np.frombuffer(zlib.decompress(existing_properties["zlib_ids"].tobytes()), dtype=np.int64)
        return iords

    def requires_property(self):
            return ["zlib_ids"]
