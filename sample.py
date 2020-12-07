import numpy as np
from constants import *
from util import *
import xraylib
import periodictable


class Composite(object):

    def __init__(self, matrix_compound, matrix_density, feature_compound, feature_density, matrix_eloss, feature_eloss,
                 matrix_thickness, feature_thickness, variable='t_b', feature_alias=None, matrix_alias=None):

        self.variable = variable
        self.type = 'composite'
        self.matrix = matrix_compound
        self.matrix_alias = matrix_alias if matrix_alias is not None else matrix_compound
        self.matrix_den = matrix_density
        self.feature = feature_compound
        self.feature_alias = feature_alias if feature_alias is not None else feature_compound
        self.feature_den = feature_density
        self.matrix_eloss = matrix_eloss
        self.feature_eloss = feature_eloss
        self.matrix_elements, self.matrix_stoic, self.matrix_mw = parser(self.matrix)
        self.feature_elements, self.feature_stoic, self.feature_mw = parser(self.feature)
        # Provide max value if variable
        self.t_b = matrix_thickness
        self.t_f = feature_thickness

    def get_ri_delta_beta(self, beam):

        delta_f = 1 - xraylib.Refractive_Index_Re(self.feature, beam.energy, self.feature_den)
        delta_b = 1 - xraylib.Refractive_Index_Re(self.matrix, beam.energy, self.matrix_den)
        beta_f = xraylib.Refractive_Index_Im(self.feature, beam.energy, self.feature_den)
        beta_b = xraylib.Refractive_Index_Im(self.matrix, beam.energy, self.matrix_den)
        return (delta_f, beta_f), (delta_b, beta_b)

    def get_ri_delta_beta_henke(self, beam):

        energy = beam.energy
        lmbda = 1.24 / energy * 1.e-9
        r0 = 2.818e-15
        delta_b = 0.
        beta_b = 0.
        delta_f = 0.
        beta_f = 0.

        pt_formula = periodictable.formula(self.matrix)
        parse_dict = pt_formula.atoms
        mw = self.matrix_mw
        density = self.matrix_den * 1.e6
        n = 1. / (mw / density) * Nav
        prefactor = n * r0 * lmbda ** 2 / (2 * np.pi)
        for i, element in enumerate(parse_dict.keys()):
            f1, f2 = element.xray.scattering_factors(energy=energy)
            s = parse_dict.values()[i]
            delta_b += prefactor * f1 * s
            beta_b += prefactor * f2 * s

        pt_formula = periodictable.formula(self.feature)
        parse_dict = pt_formula.atoms
        mw = self.feature_mw
        density = self.feature_den * 1.e6
        n = 1. / (mw / density) * Nav
        prefactor = n * r0 * lmbda ** 2 / (2 * np.pi)
        for i, element in enumerate(parse_dict.keys()):
            f1, f2 = element.xray.scattering_factors(energy=energy)
            s = parse_dict.values()[i]
            delta_f += prefactor * f1 * s
            beta_f += prefactor * f2 * s

        return (delta_f, beta_f), (delta_b, beta_b)



class SingleMaterial(object):

    def __init__(self, compound, density, thickness, eloss=None):

        self.type = 'homogeneous'
        self.compound = compound
        self.density = density
        if eloss is not None:
            self.eloss = eloss
        self.elements, self.stoic, self.mw = parser(self.compound)
        self.thickness = thickness

