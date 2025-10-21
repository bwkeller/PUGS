import tangos
import pynbody
import numpy as np
from tangos.properties import PropertyCalculation
from tangos.properties.pynbody import PynbodyPropertyCalculation

class HalfMassRedShiftLive(PropertyCalculation):
    names = "z_half_live"

    def calculate(self, particle_data, existing_properties):
        mass_now = existing_properties['finder_mass']
        z = existing_properties.timestep.redshift
        descendant = existing_properties.previous
        while descendant is not None:
            if descendant['finder_mass'] < 0.5*mass_now:
                break
            z = descendant.timestep.redshift
            descendant = descendant.previous
        return z

    def requires_property(self):
        return ["finder_mass"]

class HalfMassRedShift(PropertyCalculation):
    names = "z_half"

    def calculate(self, particle_data, existing_properties):
        return self.z[existing_properties['finder_id']]

    def preloop(self, particle_data, timestep_object):
        self.z, = timestep_object.calculate_all("z_half_live()")

    def requires_property(self):
        return ["finder_mass"]

class VirialRadius(PynbodyPropertyCalculation):
    names = "r200"
    
    def calculate(self, particle_data, existing_properties):
        with pynbody.transformation.translate(particle_data, -existing_properties['shrink_center']):
            return pynbody.analysis.halo.virial_radius(particle_data, overden=200)

    def region_specification(self, existing_properties):
        return pynbody.filt.Sphere(existing_properties['max_radius']*4,
                                   existing_properties['shrink_center'])

    def requires_property(self):
        return ["shrink_center", "max_radius"]

