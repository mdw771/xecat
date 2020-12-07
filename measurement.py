import numpy as np
from util import *
from constants import *
from beam import *


class Measurement(object):

    def __init__(self, pixel_size, n_ccd, working_distance):

        self.pixel = pixel_size
        self.n_ccd = n_ccd
        self.wd = working_distance

    def get_numerical_aperture(self, beam):

        assert isinstance(beam, XrayBeam)
        return (0.5 * beam.wavelength) / (self.pixel * 1.e-6)