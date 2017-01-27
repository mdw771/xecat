from math import sqrt,log
import numpy as np
from constants import *

def cse_inelastic(eloss,energy,elements,stoic):
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

def cse_elastic(energy,elements,stoic):
    i = 0
    cs_el = 0
    gamma = float(1.+energy/ElectronMass)
    beta = sqrt(1.-1./(gamma*gamma))
    for Z in elements:
        Z = float(Z)
        cs_el += 1.0e-14*1.4e-6*(Z**1.5)*(1.-(0.26*Z/(137.*beta)))/(beta*beta)*stoic[i]
        i += 1
    return cs_el