import xraylib
import numpy as np
from math import *
from scipy.integrate import simps
from scipy.integrate import quad
from constants import *


def parser(compounds):

    comp = xraylib.CompoundParser(compounds)
    stoic_allowed = range(10)
    stoic_allowed = map(str, stoic_allowed)
    stoic_allowed.append('.')
    _compounds = list(compounds)
    elements = comp['Elements']
    stoic = []
    for i in elements:
        symb = xraylib.AtomicNumberToSymbol(i)
        ind = compounds.find(symb) + len(symb)
        if ind >= len(_compounds):
            stoic.append(1.0)
        else:
            if _compounds[ind] not in stoic_allowed:
                stoic.append(1.0)
            else:
                temp = []
                while _compounds[ind] in stoic_allowed:
                    temp.append(_compounds[ind])
                    ind += 1
                    if ind >= len(_compounds):
                        break
                temp = ''.join(temp)
                stoic.append(float(temp))
    mw = np.sum(np.asarray(stoic) * np.asarray(map(xraylib.AtomicWeight, elements)))

    return elements, stoic, mw


def dsigma_th(theta):

    return sin(theta) * (1 + cos(theta) ** 2)


def dsigma_kn(theta, energy):

    k = energy / 510.999
    return sin(theta) * (1 + k * (1 - cos(theta))) ** (-2) * (
        1 + cos(theta) ** 2 + (k ** 2 * (1 - cos(theta)) ** 2) / (1 + k * (1 - cos(theta))))


def dsigma_el(theta, energy, mw, elements, stoic):

    dsigma = 0
    for i in range(0, len(elements)):
        q = xraylib.MomentTransf(energy, theta)
        w = xraylib.AtomicWeight(elements[i]) * stoic[i] / mw
        temp = dsigma_th(theta) * xraylib.FF_Rayl(elements[i], q) ** 2
        temp = temp / xraylib.AtomicWeight(elements[i]) * w
        dsigma += temp
    dsigma *= mw
    return dsigma


def dsigma_inel(theta, energy, mw, elements, stoic):

    dsigma = 0
    for i in range(0, len(elements)):
        q = xraylib.MomentTransf(energy, theta)
        w = xraylib.AtomicWeight(elements[i]) * stoic[i] / mw
        temp = dsigma_kn(theta, energy) * xraylib.SF_Compt(elements[i], q)
        temp = temp / xraylib.AtomicWeight(elements[i]) * w
        dsigma += temp
    dsigma *= mw
    return dsigma


def scat_efficiencies(energy, mw, elements, stoic):

    # eta is reversely scattered fraction
    theta = np.arange(1, 181, 1)
    theta = theta / 180. * pi
    theta_in = np.arange(1, 91, 1)
    theta_in = theta_in / 180. * pi

    dsigma_el_ls = [dsigma_el(thx, energy, mw, elements, stoic) for thx in theta]
    dsigma_inel_ls = [dsigma_inel(thx, energy, mw, elements, stoic) for thx in theta]
    dsigma_el_in = [dsigma_el(thx, energy, mw, elements, stoic) for thx in theta_in]
    dsigma_inel_in = [dsigma_inel(thx, energy, mw, elements, stoic) for thx in theta_in]
    dsigma_el_full = simps(dsigma_el_ls, theta)
    dsigma_inel_full = simps(dsigma_inel_ls, theta)

    eta_el = simps(dsigma_el_in, theta_in) / dsigma_el_full
    eta_inel = simps(dsigma_inel_in, theta_in) / dsigma_inel_full

    return eta_el, eta_inel


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


def cse_inelastic(eloss, energy, elements, stoic):
    """
    Calculate electron inelastic scattering cross section in cm2.
    """
    i = 0
    cs_inel = 0
    gamma = float(1.+energy/ElectronMass)
    beta = sqrt(1.-1./(gamma*gamma))
    theta = 1.0e-3*eloss/(beta*beta*(energy+ElectronMass))
    for Z in elements:
        Z = float(Z)
        if abs(Z-1) < 0.01:
            gamma100 = float(1.+100./ElectronMass)
            beta100sq = 1. - 1./gamma100**2
            const = 8.8/(1.5*log(2./theta))*beta100sq
            Z_local = const**2
            cs_inel += 1.0e-14*1.5e-6*sqrt(Z_local)*log(2./theta)/(beta*beta)*stoic[i]
        else:
            cs_inel += 1.0e-14*1.5e-6*sqrt(Z)*log(2./theta)/(beta*beta)*stoic[i]
        i += 1

    return cs_inel

def cse_elastic(energy, elements, stoic):
    """
    Calculate electron elastic scattering cross section in cm2.
    """
    i = 0
    cs_el = 0
    gamma = float(1.+energy/ElectronMass)
    beta = sqrt(1.-1./(gamma*gamma))
    for Z in elements:
        Z = float(Z)
        cs_el += 1.0e-14*1.4e-6*(Z**1.5)*(1.-(0.26*Z/(137.*beta)))/(beta*beta)*stoic[i]
        i += 1

    return cs_el