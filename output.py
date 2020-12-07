import numpy
from util import *
from constants import *


class Output(object):

    def __init__(self, sample, step=None, overflow_limit=None, resolution_range=(1e-3, 2)):

        if sample.type == 'composite':
            if sample.variable == 't_b':
                self.t_f = sample.t_f
                if step is None:
                    step = sample.t_b / 1000.
                self.step = step
                if step != 0:
                    if overflow_limit is None:
                        self.t_b = np.arange(step, sample.t_b + step, step)
                    else:
                        self.t_b = np.arange(step, overflow_limit + step, step)
                else:
                    self.t_b = np.array([sample.t_b])
            elif sample.variable == 't_f':
                self.t_b = sample.t_b
                if step is None:
                    step = sample.t_f / 1000.
                if overflow_limit is None:
                    self.t_f = np.arange(step, sample.t_f + step, step)
                else:
                    self.t_f = np.arange(step, overflow_limit + step, step)
            elif sample.variable == 'both':
                if step is None:
                    step = sample.t_f / 1000.
                if overflow_limit is None:
                    self.t_f = np.arange(step, sample.t_f + step, step)
                    self.t_b = np.arange(step, sample.t_b + step, step)
                else:
                    self.t_f = np.arange(step, overflow_limit + step, step)
                    self.t_b = np.arange(step, overflow_limit + step, step)
            elif sample.variable == 'none':
                self.t_f = sample.t_f
                self.t_b = sample.t_b
                self.step = None
            else:
                raise ValueError
            self.t = self.t_b + self.t_f
        else:
            if step is None:
                step = sample.thickness / 1000.
            self.step = step
            if overflow_limit is None:
                self.t = np.arange(step, sample.thickness + step, step)
            else:
                self.t = np.arange(step, overflow_limit + step, step)
        self.delta = np.arange(resolution_range[0], resolution_range[1], 1e-3)
