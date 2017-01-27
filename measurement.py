import numpy as np
from util import *
from constants import *


class Measurement(object):

    def __init__(self, pixel_size, n_ccd, working_distance):

        self.pixel = pixel_size
        self.n_ccd = n_ccd
        self.wd = working_distance