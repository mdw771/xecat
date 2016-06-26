import sys
import xraylib
from math import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.integrate import dblquad
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
fig = plt.figure(figsize=(14,28))
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

    # frac1 is the fraction in FORWARD SCATTERED (detectable) photons that enter the aperture and contribute to PCI signal
    def fsquare(phi, theta)
        return (1-sin(atan(tan(phi)*sin(theta)))**2*cos(acos(sin(phi)*cos(theta)))**2)*sin(phi)


    # establish table of thickness in nm (must be integer)
    thickness = thickls[fig_index-1]
    step = stepls[fig_index-1]
    t = range(step, thickness + 1, step)

    # frac2 is the fraction of FORWARD scattered inelastic scattered or plural scattered photons that enter the aperture
    # and contribute to PCI background
    # for each t calculate frac2
    # all units in um except t
    zstep = 250./1000
    NAa = tan(asin(NA))*a
    frac2 = [0 for i in t]
    index = 0
    def tri_a1():
        return 2*NAa
    def tri_b1(x, y, z):
        return sqrt((sqrt(x**2+y**2) + NAa)**2 + (a+z)**2)
    def tri_c1(x, y, z):
        return sqrt((sqrt(x**2+y**2) - NAa)**2 + (a+z)**2)
    def tri_a2():
        return 2*NAa
    def tri_b2(x, y, z):
        return sqrt(x**2 + y**2 + (a+z)**2 + NAa**2)
    def tri_c2(x, y, z):
        return sqrt(x**2 + y**2 + (a+z)**2 + NAa**2)
    def theta(x, y, z):
        return (acos((tri_b1(x, y, z)**2 + tri_c1(x, y, z)**2 - tri_a1()**2)/(2*tri_b1(x, y, z)*tri_c1(x, y, z))) +
                acos((tri_b2(x, y, z)**2 + tri_c2(x, y, z)**2 - tri_a2()**2)/(2*tri_b2(x, y, z)*tri_c2(x, y, z))))/2
    for i in t:
        print "%i / %i\n" % (i, thickness)
        ti = i/1000.
        Nt = ti/zstep
        volume = (pixel*N)**2*ti
        def integrand(x, y, z):
            return 1/volume*(1-cos(theta(x,y,z)))
        result = tplquad(integrand, 0, ti, lambda z: -N*pixel/2, lambda z: N*pixel/2, lambda z, y: -N*pixel/2, lambda
            z, y: N*pixel/2)
        frac2[index] = result[0]
        index = index + 1

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
    ###i_1elnopc = [k_elinnopc*t[i/step-1]*(1e-7)*i_noscat[i/step-1] for i in t]
    i_elpl = [exp(-(k_inelin + k_out + k_pi)*t[i/step-1]*(1e-7)) - i_noscat[i/step-1] - i_1el[i/step-1] for i in t]
    i_elplpc = [i_elpl[i/step-1]*frac2[i/step-1] for i in t]
    ###i_elplnopc = [i_elpl[i/step-1]*(1 - frac2[i/step-1]) for i in t]
    i_out = [(k_out/(k_out+k_pi))*(1 - exp(-(k_out+k_pi)*t[i/step-1]*(1e-7))) for i in t]
    i_pi = [(k_pi/(k_out+k_pi))*(1 - exp(-(k_out+k_pi)*t[i/step-1]*(1e-7))) for i in t]
    i_inel = [exp(-(k_out + k_pi)*t[i/step-1]*(1e-7)) - exp(-(k_inelin + k_out + k_pi)*t[i/step-1]*(1e-7)) for i in t]
    i_inelpc = [i_inel[i/step-1]*frac2[i/step-1] for i in t]

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
    # pl_1elnopc, = plt.plot(t, i_1elnopc, label = 'Single elastically scattered not contributing to PCI')
    # pl_1elnopc.set_dashes([1, 2])
    pl_elplpc, = plt.plot(t, i_elplpc, label = 'Plural elastically scattered in PCI background')
    pl_elplpc.set_dashes([3, 1])
    # pl_elplnopc, = plt.plot(t, i_elplnopc, label = 'Plural elastically scattered not contributing to PCI')
    # pl_elplnopc.set_dashes([3, 2])
    pl_out, = plt.plot(t, i_out, label = 'Scattered out')
    pl_out.set_dashes([1, 2, 1, 2, 1, 2, 5, 2])
    pl_pi, = plt.plot(t, i_pi, label = 'Absorbed')
    pl_pi.set_dashes([1, 2, 1, 2, 4, 2, 4, 2])
    pl_inel, = plt.plot(t, i_inel, label = 'Inelastically scattered')
    pl_inel.set_dashes([1, 2, 4, 2])
    pl_inelpc, = plt.plot(t, i_inelpc, label = 'Inelastically scattered in PCI background')
    pl_inelpc.set_dashes([1, 3, 5, 3])
    plt.semilogy()
    plt.xlabel('Thickness (nm)')
    plt.ylabel('Fraction')
    fig_index = fig_index + 1
    if energy == energyls[-1]:
        plt.legend(bbox_to_anchor = (0,-0.8), loc = 2)

plt.savefig('unimatrix_fig.pdf', format = 'pdf')
plt.show()
sys.exit()

