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
    return exp(-7.256 * s ** 2)

def higha(s):
    return 0.0168 * s ** -4 - 0.1514 * s ** -3 + 0.4397 * s ** -2 - 0.2441 * s ** -1 + 0.0399

def saxs_polynomial(s):
    return 10 ** (0.007466 * s ** 5 - 0.1068 * s ** 4 + 0.5305 * s ** 3 - 0.8888 * s ** 2 - 0.6426 * s + 0.08601)

def inna_frac(wavelen, NA):

    # s is the momentum transfer in nm-1
    s_max = 4 * pi * sin(pi / 4) / (wavelen * 1e9)
    # s_glim = 1.3 / 4.1888
    s_na = 4 * pi * sin(NA / 2) / (wavelen * 1e9)
    if s_na > s_max:
        raise ValueError('NA is too large. Booooooooooooooom')
    int_top, _ = quad(saxs_polynomial, 0, min([s_na, 3]))
    int_full, _ = quad(saxs_polynomial, 0, min([s_max, 3]))
    frac1 = int_top / int_full

    return frac1


    # full_guinier, _ = quad(guinier, 0, s_glim)
    # full_ha, _ = quad(higha, s_glim, 3)
    # if s_na < s_glim:
    #     na_int, _ = quad(guinier, 0, s_na)
    # else:
    #     na_int, _ = quad(higha, s_glim, s_na)
    #     na_int = na_int + full_guinier
    # if s_max > 3:
    #     frac1 = na_int / (full_guinier + full_ha)
    # elif s_max < s_glim:
    #     max_int, _ = quad(guinier, 0, s_max)
    #     frac1 = na_int / max_int
    # else:
    #     max_int, _ = quad(higha, s_glim, s_max)
    #     max_int += full_guinier
    #     frac1 = na_int / max_int


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
    beta = sqrt(1. - 1./(gamma*gamma))
    for Z in elements:
        Z = float(Z)
        cs_el += 1.0e-14*1.4e-6*(Z**1.5)*(1.-(0.26*Z/(137.*beta)))/(beta*beta)*stoic[i]
        i += 1

    return cs_el


def rudolph_theta(mu_b, t_b, mu_f, t_f, delta_b, delta_f, wavelength_um, output_dichotomy=False):

    etabt = 2 * np.pi * delta_b / wavelength_um * t_f
    etaft = 2 * np.pi * delta_f / wavelength_um * t_f
    mubt = mu_b * t_f
    muft = mu_f * t_f
    mubttf = mu_b * t_b

    # i_f = 2 * np.exp(-mu_b * t) + np.exp(-mu_f * t) + 2 * np.exp(-mu_f / 2 * t) * np.exp(-mu_b / 2 * t) * \
    #                                                       ((eta_f * t - eta_b * t) - 1)
    # i_b = np.exp(-mu_b * t)

    i_f = np.exp(-mubttf) * (2. * np.exp(-mubt) + np.exp(-muft) + 2. * np.exp(-0.5 * muft - 0.5 * mubt)
                             * np.cos(etaft - etabt - 0.5*np.pi) - 2. * np.exp(-0.5 * muft - 0.5 * mubt) * np.cos(etaft - etabt)
                             - 2. * np.exp(-mubt) * np.cos(0.5*np.pi))
    i_b = np.exp(-mubttf) * np.exp(-mubt)

    theta = np.abs(i_f - i_b) / np.sqrt(i_f + i_b)

    if output_dichotomy:
        return theta, i_f, i_b
    else:
        return theta

      # i_f_plusphase(0:(enum-1),i_t) = exp(-mu2ttf(0:(enum-1),i_t)) * $
      #   ( (1.+exp(-mupt_plus))*exp(-mu2t) + exp(-mu1t) + $
      #     2.*exp(-0.5*mu1t-0.5*mu2t-0.5*mupt_plus)* $
      #     cos(eta1t-eta2t-etapt_plus) - $
      #     2.*exp(-0.5*mu1t-0.5*mu2t)*cos(eta1t-eta2t) - $
      #     2.*exp(-mu2t-0.5*mupt_plus)*cos(etapt_plus) )
      # i_b_plusphase(0:(enum-1),i_t) = exp(-mu2ttf(0:(enum-1),i_t)) * $
      #   exp(-mu2t-mupt_plus)