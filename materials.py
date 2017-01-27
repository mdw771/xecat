from util import *

class single_material(object):

    def __init__(self, formula, density, ineleloss=None):

        self.formula = formula
        self.density = density
        if ineleloss is not None:
            self.ineleloss = ineleloss
        self.matrix_elements, self.matrix_stoic, self.matrix_mw = parser(formula)
