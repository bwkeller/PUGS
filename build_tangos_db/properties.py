import zlib

import numpy as np
import pynbody
from tangos.properties import LivePropertyCalculation
from tangos.properties.pynbody import PynbodyPropertyCalculation


class StoreIords(PynbodyPropertyCalculation):
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
    names = "ids"

    def calculate(self, particle_data, existing_properties):
        iords = np.frombuffer(
            zlib.decompress(existing_properties["zlib_ids"].tobytes()), dtype=np.int64
        )
        return iords

    def requires_property(self):
        return ["zlib_ids"]
