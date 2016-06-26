import sys
import xraylib
from math import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.integrate import quad
import csv

# define constants
Nav = 6.02214179e+23
PlanckConst = 6.62607004e-34
SpeedOfLight = 299792458.
ECharge = 1.60217662e-19
pi = 3.14159

# define pixel size in um and instrumental constants: working distance (a) in um, imaging system dimension
pixel = 1
a = 1e4
N = 1024

# beam energy in keV
energyls = [5, 10, 20, 40]
nenergy = len(energyls)

# max thickness for all energies in integer nm
thickls = [25000000, 25000000, 25000000, 25000000]
stepls = [10000, 10000, 10000, 10000]

# matrix compound information
matrix = "H48.6C32.9N8.9O8.9S0.6"
matrix_den = 1.35
elements = [1, 6, 7, 8, 16]
stoic = [48.6, 32.9, 8.9, 8.9, 0.6]
mw = 0.
j = 0
for i in elements:
    mw += xraylib.AtomicWeight(i)*stoic[j]
    j += 1

fig_index = 1
fig = plt.figure(figsize=(10,50))
matplotlib.rcParams['pdf.fonttype'] = 'truetype'
fontProperties = {'family' : 'serif', 'serif' : ['Times New Roman'], 'weight' : 'normal', 'size' : 12}
plt.rc('font', **fontProperties)

f = open('unimatrix_data.csv','wb')


for energy in energyls:
    # eta is reversely scattered fraction
    # NA based on typical resolution of 1 um
    # Abbe criterion applied
    eta = 0.5
    wavelen = PlanckConst*SpeedOfLight/(energy*1000*ECharge)
    NA = (0.5*wavelen)/(pixel*1e-6)
    print NA

    # frac1 is the fraction in FORWARD SCATTERED (detectable) photons that enter the aperture and contribute to PCI signal
    # s in nm^-1
    s_max = 4*pi*sin(pi/4)/(wavelen*1e9)
    s_glim = 1.3/4.1888
    s_na = 4*pi*sin(NA/2)/(wavelen*1e9)
    if s_max < 3:
        print 'WARNING: Check s_max.'
        sys.exit()
    def guinier(s):
        return exp(-1.396*s**2)
    def higha(s):
        return 0.0168*s**-4 - 0.1514*s**-3 + 0.4397*s**-2 - 0.2441*s**-1 + 0.0399
    full_guinier, err = quad(guinier, 0, s_glim)
    full_ha, err = quad(higha, s_glim, 3)
    if s_na < s_glim:
        na_int, err = quad(guinier, 0, s_na)
    else:
        na_int, err = quad(higha, s_glim, s_na)
        na_int = na_int + full_guinier
    frac1 = na_int/(full_guinier+full_ha)

    # establish table of thickness in nm (must be integer)
    thickness = thickls[fig_index-1]
    step = stepls[fig_index-1]
    t = range(step, thickness + 1, step)

    # frac2 is the fraction of FORWARD scattered inelastic scattered or plural scattered photons that enter the aperture
    # and contribute to PCI background
    frac2 = NA**2/2

    # retrieve cross sections in cm2/g
    cs_inel = xraylib.CS_Compt_CP(matrix, energy)
    cs_el = xraylib.CS_Rayl_CP(matrix, energy)
    cs_pi = xraylib.CS_Photo_CP(matrix, energy)
    cs_tot = xraylib.CS_Total_CP(matrix, energy)

    # probability per thickness
    # eta is the fraction of photons scattered out of detectable range (backscattered)
    k_el = cs_el*matrix_den
    k_inel = cs_inel*matrix_den
    k_elin = cs_el*(1-eta)*matrix_den
    k_inelin = cs_inel*(1-eta)*matrix_den
    k_elinnopc = cs_el*(1-eta)*(1-frac1)*matrix_den
    k_out = cs_el*eta*matrix_den + cs_inel*eta*matrix_den
    k_pi = cs_pi*matrix_den
    k_tot = k_inel + k_el + k_pi

    # intensity fractions relative to primary beam
    # t is converted to cm
    i_noscat = [exp(-k_tot*t[i/step-1]*(1e-7)) for i in t]
    i_1el = [k_elin*t[i/step-1]*(1e-7)*i_noscat[i/step-1] for i in t]
    i_1elpc = [frac1*i_1el[i/step-1] for i in t]
    i_pc = [sqrt(i_noscat[i/step-1]*i_1elpc[i/step-1]) for i in t]
    i_elpl = [exp(-(k_inelin + k_out + k_pi)*t[i/step-1]*(1e-7)) - i_noscat[i/step-1] - i_1el[i/step-1] for i in t]
    i_elplpc = [i_elpl[i/step-1]*frac2 for i in t]
    i_out = [(k_out/(k_out+k_pi))*(1 - exp(-(k_out+k_pi)*t[i/step-1]*(1e-7))) for i in t]
    i_pi = [(k_pi/(k_out+k_pi))*(1 - exp(-(k_out+k_pi)*t[i/step-1]*(1e-7))) for i in t]
    i_inel = [exp(-(k_out + k_pi)*t[i/step-1]*(1e-7)) - exp(-(k_inelin + k_out + k_pi)*t[i/step-1]*(1e-7)) for i in t]
    i_inelpc = [i_inel[i/step-1]*frac2 for i in t]

    # output data to csv
    writer = csv.writer(f)
    writer.writerows([["Energy (keV)", energy], t, i_pc, i_elplpc, i_out, i_pi, i_inelpc, []])

    # run report
    print "---------------------"
    print "Energy: %.2f keV" % energy
    print "Thickness limit: %.2f nm (%.2f cm)" % (thickness, thickness/1e7)
    print "Thickness resolution: %.2f nm" % step
    print "---------------------"

    # plot intensities
    plt.subplot(nenergy+2,1,fig_index)

    pl_pc, = plt.plot(t, i_pc, label = 'PCI signal')
    pl_pc.set_dashes([1, 1])
    pl_elplpc, = plt.plot(t, i_elplpc, label = 'Plural elastically scattered in PCI background')
    pl_elplpc.set_dashes([3, 1])
    pl_inelpc, = plt.plot(t, i_inelpc, label = 'Inelastically scattered in PCI background')
    pl_inelpc.set_dashes([1, 3, 5, 3])
    pl_inel, = plt.plot(t, i_inel, label = 'Inelastically scattered')
    pl_inel.set_dashes([1, 2, 4, 2])
    pl_out, = plt.plot(t, i_out, label = 'Scattered out')
    pl_out.set_dashes([1, 2, 1, 2, 1, 2, 5, 2])
    pl_pi, = plt.plot(t, i_pi, label = 'Absorbed')
    pl_pi.set_dashes([1, 2, 1, 2, 4, 2, 4, 2])
    plt.semilogy()
    plt.xlabel('Thickness (nm)')
    plt.ylabel('Fraction')
    fig_index = fig_index + 1
    if energy == energyls[-1]:
        plt.legend(bbox_to_anchor = (0,-0.8), loc = 2)

plt.savefig('unimatrix_fig.eps', format = 'eps')

plt.show()
sys.exit()

