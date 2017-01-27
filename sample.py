import numpy as np
from constants import *
from util import *


class Composite(object):

    def __init__(self, matrix_compound, matrix_density, feature_compound, feature_density, matrix_eloss, feature_eloss,
                 thickness, feature_thickness):

        self.type = 'composite'
        self.matrix = matrix_compound
        self.matrix_den = matrix_density
        self.feature = feature_compound
        self.feature_den = feature_density
        self.matrix_eloss = matrix_eloss
        self.feature_eloss = feature_eloss
        self.matrix_elements, self.matrix_stoic, self.matrix_mw = parser(self.matrix)
        self.feature_elements, self.feature_stoic, self.feature_mw = parser(self.feature)
        self.thickness = thickness
        self.t_f = feature_thickness


class SingleMaterial(object):

    def __init__(self, compound, density, thickness, eloss=None):

        self.type = 'homogeneous'
        self.compound = compound
        self.density = density
        if eloss is not None:
            self.eloss = eloss
        self.elements, self.stoic, self.mw = parser(self.compound)
        self.thickness = thickness
