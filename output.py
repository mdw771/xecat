import numpy
from util import *
from constants import *


class Output(object):

    def __init__(self, sample, step=None, overflow_limit=None):

        if step is None:
            step = sample.thickness / 1000.
        self.step = step
        if overflow_limit is None:
            self.t = np.arange(step, sample.thickness + step, step)
        else:
            self.t = np.arange(step, overflow_limit + step, step)
        if sample.type == 'composite':
            self.t_b = self.t - sample.t_f


