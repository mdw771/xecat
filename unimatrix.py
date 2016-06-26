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
thickls = [int(1e7), int(100e6), int(500e6), int(1000e6)]
stepls = [10000, 20000, 40000, 80000]

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

# prepare plots
fig_index = 1
fig_lab = ['(a)','(b)','(c)','(d)']
fig = plt.figure(figsize=(12,10))
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
    k_el = cs_el*matrix_den
    k_inel = cs_inel*matrix_den
    k_elin = cs_el*(1-eta)*matrix_den
    k_inelin = cs_inel*(1-eta)*matrix_den
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
    #writer = csv.writer(f)
    #writer.writerows([["Energy (keV)", energy], t, i_pc, i_elplpc, i_out, i_pi, i_inelpc, []])

    # run report
    print "---------------------"
    print "Energy: %.2f keV" % energy
    print "Thickness limit: %.2f nm (%.2f cm)" % (thickness, thickness/1e7)
    print "Thickness resolution: %.2f nm" % step
    print "---------------------"

    # plot intensities
    plt.subplot(ceil(nenergy/2),2,fig_index)

    t_cm = [t[i/step-1]*1e-7 for i in t]

    label_pos = int(thickness/step/3)
    plt.title(fig_lab[fig_index-1]+' Protein at '+str(energy)+' keV')

    pl_pc, = plt.plot(t_cm, i_pc, label = 'Phase contrast image (PCI) signal', color='blue',linewidth=2.0)
    if fig_index == 1:
        plt.text(0.41,1.8e-13,'Phase contrast\nimage (PCI) signal', fontsize=9, color='blue',rotation=-40)
    if fig_index == 2:
        plt.text(3.4,1.8e-13,'Phase contrast\nimage (PCI) signal', fontsize=9, color='blue',rotation=-50.2)
    if fig_index == 3:
        plt.text(22.5,1.8e-13,'Phase contrast\nimage (PCI) signal', fontsize=9, color='blue',rotation=-44)
    if fig_index == 4:
        plt.text(40,6e-9,'Phase contrast image (PCI) signal', fontsize=9, color='blue',rotation=-34.2)

    pl_elplpc, = plt.plot(t_cm, i_elplpc, label = 'Plural elastically scattered in PCI background', color='green',linewidth=2.0)
    if fig_index == 1:
        plt.text(0.02,3e-16,'Plural elastically\nscattered in PCI\nbackground', fontsize=9, color='green',rotation=-35)
    if fig_index == 2:
        plt.text(0.4,3e-16,'Plural elastically\nscattered in PCI\nbackground', fontsize=9, color='green',rotation=-45)
    if fig_index == 3:
        plt.text(1.8,1.5e-14,'Plural elastically scattered\nin PCI background', fontsize=9, color='green',rotation=-38)
    if fig_index == 4:
        plt.text(9,2e-15,'Plural elastically scattered\nin PCI background', fontsize=9, color='green',rotation=-27)

    pl_inelpc, = plt.plot(t_cm, i_inelpc, label = 'Inelastically scattered in PCI background', color='red',linewidth=2.0)
    if fig_index == 1:
        plt.text(0.12,5e-13,'Inelastically scattered\nin PCI background', fontsize=9, color='red',rotation=-40)
    if fig_index == 2:
        plt.text(1,5e-13,'Inelastically scattered\nin PCI background', fontsize=9, color='red',rotation=-47)
    if fig_index == 3:
        plt.text(3,5e-11,'Inelastically scattered in PCI background', fontsize=9, color='red',rotation=-38)
    if fig_index == 4:
        plt.text(14,1e-11,'Inelastically scattered in PCI background', fontsize=9, color='red',rotation=-22)

    pl_noscat, = plt.plot(t_cm, i_noscat, label = 'Unscattered', color='salmon')
    if fig_index == 1:
        plt.text(0.51,1e-10,'Unscattered', fontsize=9, color='salmon',rotation=-43)
    if fig_index == 2:
        plt.text(1.15,1e-3,'Unscattered', fontsize=9, color='salmon',rotation=-50)
    if fig_index == 3:
        plt.text(3.6,3e-2,'Unscattered', fontsize=9, color='salmon',rotation=0)
    if fig_index == 4:
        plt.text(21,1e-3,'Unscattered', fontsize=9, color='salmon',rotation=-35)

    pl_1el, = plt.plot(t_cm, i_1el, label = 'Single elastically scattered', color='black')
    if fig_index == 1:
        plt.text(0.09,1e-3,'Single elastically scattered', fontsize=9, color='black',rotation=-42)
    if fig_index == 2:
        plt.text(t_cm[label_pos]+0.18,i_1el[label_pos],'Single elastically scattered', fontsize=9, color='black',rotation=-50)
    if fig_index == 3:
        plt.text(t_cm[label_pos]+0.4,i_1el[label_pos],'Single elastically scattered', fontsize=9, color='black',rotation=-43)
    if fig_index == 4:
        plt.text(60,3e-9,'Single elastically scattered', fontsize=9, color='black',rotation=-33)

    pl_elpl, = plt.plot(t_cm, i_elpl, label = 'Plural elastically scattered', color='darkcyan')
    if fig_index == 1:
        plt.text(0.1,4e-8,'Plural elastically\nscattered', fontsize=9, color='darkcyan',rotation=-40)
    if fig_index == 2:
        plt.text(0.1,8e-7,'Plural\nelastically\nscattered', fontsize=9, color='darkcyan',rotation=-36)
    if fig_index == 3:
        plt.text(32,1.5e-13,'Plural elastically scattered', fontsize=9, color='darkcyan',rotation=-43)
    if fig_index == 4:
        plt.text(3,8e-6,'Plural\nelastically\nscattered', fontsize=9, color='darkcyan',rotation=-11)

    pl_inel, = plt.plot(t_cm, i_inel, label = 'Inelastically scattered', color='magenta')
    if fig_index == 1:
        plt.text(0.04,1e-4,'Inelastically scattered', fontsize=9, color='magenta',rotation=-39)
    if fig_index == 2:
        plt.text(6,3e-16,'Inelastically\nscattered', fontsize=9, color='magenta',rotation=-52)
    if fig_index == 3:
        plt.text(t_cm[label_pos]+0.4,i_inel[label_pos],'Inelastically scattered', fontsize=9, color='magenta',rotation=-39)
    if fig_index == 4:
        plt.text(t_cm[label_pos]+2,i_inel[label_pos],'Inelastically scattered', fontsize=9, color='magenta',rotation=-22)

    pl_out, = plt.plot(t_cm, i_out, label = 'Scattered out', color='orange')
    if fig_index == 1:
        plt.text(t_cm[label_pos],i_out[label_pos]+4e-3,'Scattered out', fontsize=9, color='orange',rotation=0)
    if fig_index == 2:
        plt.text(t_cm[label_pos],i_out[label_pos]+2e-2,'Scattered out', fontsize=9, color='orange',rotation=0)
    if fig_index == 3:
        plt.text(t_cm[label_pos],i_out[label_pos]-0.15,'Scattered out', fontsize=9, color='orange',rotation=0)
    if fig_index == 4:
        plt.text(t_cm[label_pos],i_out[label_pos]-0.6,'Scattered out', fontsize=9, color='orange',rotation=0)

    pl_pi, = plt.plot(t_cm, i_pi, label = 'Absorbed', color='grey')
    if fig_index == 1:
        plt.text(t_cm[label_pos]+0.48,i_pi[label_pos]-0.74,'Absorbed', fontsize=9, color='grey',rotation=0)
    if fig_index == 2:
        plt.text(t_cm[label_pos]+3.4,i_pi[label_pos]-0.66,'Absorbed', fontsize=9, color='grey',rotation=0)
    if fig_index == 3:
        plt.text(t_cm[label_pos],i_pi[label_pos]-0.5,'Absorbed', fontsize=9, color='grey',rotation=0)
    if fig_index == 4:
        plt.text(t_cm[label_pos]+30,i_pi[label_pos]-0.26,'Absorbed', fontsize=9, color='grey',rotation=0)

    plt.semilogy()
    plt.ylim(1e-18,1.1)
    plt.xlabel('Thickness (cm)')
    plt.ylabel('Fraction')
    fig_index = fig_index + 1

fig.savefig('unimatrix_fig.pdf', format = 'pdf')

plt.show()

