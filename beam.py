import numpy as np
from constants import *
from util import *


class XrayBeam(object):

    def __init__(self, energy):

        self.energy = energy
        self.wavelength = PlanckConst * SpeedOfLight / (energy * 1000 * ECharge)


class ElectronBeam(object):

    def __init__(self, energy):

        self.energy = energy