### the modules returns the X-ray elastic/inelastic efficiency eta_el and eta_inel
import numpy as np
import xraylib
from math import *
from scipy.integrate import simps
import matplotlib.pyplot as plt


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

    # plt.figure()
    # plt.plot(theta, dsigma_inel_ls)
    # plt.show()

    dsigma_el_full = simps(dsigma_el_ls, theta)
    dsigma_inel_full = simps(dsigma_inel_ls, theta)

    eta_el = simps(dsigma_el_in, theta_in) / dsigma_el_full
    eta_inel = simps(dsigma_inel_in, theta_in) / dsigma_inel_full

    return (eta_el, eta_inel)
