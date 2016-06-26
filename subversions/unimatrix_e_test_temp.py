import sys
from math import *
from xraylib import AtomicWeight
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from scipy.integrate import quad
import csv
from matplotlib.ticker import AutoMinorLocator


# define constants
Nav = 6.02214179e+23
PlanckConst = 6.62607004e-34
SpeedOfLight = 299792458. # m/s
ECharge = 1.60217662e-19 # C
EMass = 9.10938356e-31 # kg
pi = 3.14159

# define pixel size in um and instrumental constants: working distance (a) in um, imaging system dimension
pixel = 1
a = 1e4
N = 1024

# beam energy in keV
energyls = [300, 1000]
nenergy = len(energyls)

# max thickness for all energies in integer nm
thickls = [int(3000), int(3000)]
stepls = [5,15]

# matrix compound information
matrix_den = 0.92
elements = [1,8]
stoic = [2,1]
mw = 0.
j = 0
for i in elements:
    mw += AtomicWeight(i)*stoic[j]
    j += 1

# prepare plots
fig_index = 1
fig_lab = ['(a)','(b)']
fig = plt.figure(figsize=(12,5))
matplotlib.rcParams['pdf.fonttype'] = 'truetype'
fontProperties = {'family' : 'serif', 'serif' : ['Times New Roman'], 'weight' : 'normal', 'size' : 12}
plt.rc('font', **fontProperties)

#f = open('unimatrix_data.csv','wb')

for energy in energyls:
    # eta is the out-scattered fraction
    eta = 1 - 4.12/10

    # establish table of thickness in nm (must be integer)
    thickness = thickls[fig_index-1]
    step = stepls[fig_index-1]
    t = range(step, thickness + 1, step)

    # compute cross sections in cm2 (see Langmore 1992)
    i = 0
    cs_el = 0
    cs_inel = 0
    eloss = 39.3 # eV
    electron_mass = 511.003414
    gamma = float(1.+energy/electron_mass)
    beta = sqrt(1.-1./(gamma*gamma))
    theta = 1.0e-3*eloss/(beta*beta*(energy+electron_mass))
    for Z in elements:
        Z = float(Z)
        cs_el += 1.0e-14*1.4e-6*(Z**1.5)*(1.-(0.26*Z/(137.*beta)))/(beta*beta)*stoic[i]
        if abs(Z-1) < 0.01:
            gamma100 = float(1.+100./electron_mass)
            beta100sq = 1. - 1./gamma100**2
            const = 8.8/(1.5*log(2./theta))*beta100sq
            Z_local = const**2
            cs_inel += 1.0e-14*1.5e-6*sqrt(Z_local)*log(2./theta)/(beta*beta)*stoic[i]
        else:
            cs_inel += 1.0e-14*1.5e-6*sqrt(Z)*log(2./theta)/(beta*beta)*stoic[i]
        i += 1

    # probability per thickness
    # eta is the fraction of photons scattered out of detectable range (backscattered)
    # delta is the number density of molecules (cm-3)
    delta = matrix_den/mw*Nav
    k_el = cs_el*delta
    k_inel = cs_inel*delta
    k_elin = cs_el*(1-eta)*delta
    k_out = cs_el*eta*delta
    k_tot = k_inel + k_elin + k_out

    # intensity fractions relative to primary beam
    # t is converted to cm
    i_noscat = [exp(-k_tot*t[i/step-1]*(1e-7)) for i in t]
    i_1el = [k_elin*t[i/step-1]*(1e-7)*i_noscat[i/step-1] for i in t]
    i_elpl = [exp(-(k_inel + k_out)*t[i/step-1]*(1e-7)) - i_noscat[i/step-1] - i_1el[i/step-1] for i in t]
    i_out = [1 - exp(-k_out*t[i/step-1]*(1e-7)) for i in t]
    i_inel = [exp(-k_out*t[i/step-1]*(1e-7)) - exp(-(k_inel + k_out)*t[i/step-1]*(1e-7)) for i in t]

    # output data to csv
    # writer = csv.writer(f)
    # writer.writerows([["Energy (keV)", energy], t, i_pc, i_elplpc, i_out, i_pi, i_inelpc, []])

    # run report
    print "---------------------"
    print "Energy: %.2f keV" % energy
    print "Thickness limit: %.2f nm" % (thickness)
    print "Thickness resolution: %.2f nm" % step
    print "---------------------"

    # plot intensities
    minorLocator = AutoMinorLocator()
    ax = plt.subplot(ceil(nenergy/2),2,fig_index)
    ax.xaxis.set_minor_locator(minorLocator)
    label_pos = int(thickness/step/3)
    plt.title(fig_lab[fig_index-1]+' Ice at '+str(energy)+' keV')

    pl_noscat, = plt.plot(t, i_noscat, label = 'Unscattered', color='salmon')
    if fig_index == 1:
        plt.text(0.51,1e-10,'Unscattered', fontsize=9, color='salmon',rotation=-43)
    if fig_index == 2:
        plt.text(1.15,1e-3,'Unscattered', fontsize=9, color='salmon',rotation=-50)
    if fig_index == 3:
        plt.text(3.6,3e-2,'Unscattered', fontsize=9, color='salmon',rotation=0)
    if fig_index == 4:
        plt.text(21,1e-3,'Unscattered', fontsize=9, color='salmon',rotation=-35)

    pl_1el, = plt.plot(t, i_1el, label = 'Single elastically scattered', color='black')
    if fig_index == 1:
        plt.text(0.09,1e-3,'Single elastically scattered', fontsize=9, color='black',rotation=-42)
    if fig_index == 2:
        plt.text(t[label_pos]+0.18,i_1el[label_pos],'Single elastically scattered', fontsize=9, color='black',rotation=-50)
    if fig_index == 3:
        plt.text(t[label_pos]+0.4,i_1el[label_pos],'Single elastically scattered', fontsize=9, color='black',rotation=-43)
    if fig_index == 4:
        plt.text(60,3e-9,'Single elastically scattered', fontsize=9, color='black',rotation=-33)

    pl_elpl, = plt.plot(t, i_elpl, label = 'Plural elastically scattered', color='darkcyan')
    if fig_index == 1:
        plt.text(0.1,4e-8,'Plural elastically\nscattered', fontsize=9, color='darkcyan',rotation=-40)
    if fig_index == 2:
        plt.text(0.1,8e-7,'Plural\nelastically\nscattered', fontsize=9, color='darkcyan',rotation=-36)
    if fig_index == 3:
        plt.text(32,1.5e-13,'Plural elastically scattered', fontsize=9, color='darkcyan',rotation=-43)
    if fig_index == 4:
        plt.text(3,8e-6,'Plural\nelastically\nscattered', fontsize=9, color='darkcyan',rotation=-11)

    pl_inel, = plt.plot(t, i_inel, label = 'Inelastically scattered', color='magenta')
    if fig_index == 1:
        plt.text(0.04,1e-4,'Inelastically scattered', fontsize=9, color='magenta',rotation=-39)
    if fig_index == 2:
        plt.text(6,3e-16,'Inelastically\nscattered', fontsize=9, color='magenta',rotation=-52)
    if fig_index == 3:
        plt.text(t[label_pos]+0.4,i_inel[label_pos],'Inelastically scattered', fontsize=9, color='magenta',rotation=-39)
    if fig_index == 4:
        plt.text(t[label_pos]+2,i_inel[label_pos],'Inelastically scattered', fontsize=9, color='magenta',rotation=-22)

    pl_out, = plt.plot(t, i_out, label = 'Scattered out', color='orange')
    if fig_index == 1:
        plt.text(t[label_pos],i_out[label_pos]+4e-3,'Scattered out', fontsize=9, color='orange',rotation=0)
    if fig_index == 2:
        plt.text(t[label_pos],i_out[label_pos]+2e-2,'Scattered out', fontsize=9, color='orange',rotation=0)
    if fig_index == 3:
        plt.text(t[label_pos],i_out[label_pos]-0.15,'Scattered out', fontsize=9, color='orange',rotation=0)
    if fig_index == 4:
        plt.text(t[label_pos],i_out[label_pos]-0.6,'Scattered out', fontsize=9, color='orange',rotation=0)

    plt.semilogy()
    plt.ylim(1e-3,1.1)
    plt.xlabel('Thickness (nm)')
    plt.ylabel('Fraction')
    fig_index = fig_index + 1

fig.savefig('unimatrix_fig_e.pdf', format = 'pdf')

plt.show()
sys.exit()

