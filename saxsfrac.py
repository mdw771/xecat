### this module returns frac1, the fraction in FORWARD SCATTERED (detectable) photons that enter the aperture and
### contribute to PCI signal

import sys
from math import *
from scipy.integrate import quad
import numpy as np


def guinier(s):
    return exp(-1.396 * s ** 2)


def higha(s):
    return 0.0168 * s ** -4 - 0.1514 * s ** -3 + 0.4397 * s ** -2 - 0.2441 * s ** -1 + 0.0399


def inna_frac(wavelen, NA):

    # s is the momentum transfer in nm-1
    s_max = 4 * pi * sin(pi / 4) / (wavelen * 1e9)
    s_glim = 1.3 / 4.1888
    s_na = 4 * pi * sin(NA / 2) / (wavelen * 1e9)
    if s_na > s_max:
        print 'WARNING: Check NA.'
        sys.exit()

    full_guinier, err = quad(guinier, 0, s_glim)
    full_ha, err = quad(higha, s_glim, 3)
    if s_na < s_glim:
        na_int, err = quad(guinier, 0, s_na)
    else:
        na_int, err = quad(higha, s_glim, s_na)
        na_int = na_int + full_guinier
    if s_max > 3:
        frac1 = na_int / (full_guinier + full_ha)
    elif s_max < s_glim:
        max_int, err = quad(guinier, 0, s_max)
        frac1 = na_int / max_int
    else:
        max_int, err = quad(higha, s_glim, s_max)
        max_int += full_guinier
        frac1 = na_int / max_int

    return frac1
