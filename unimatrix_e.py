import sys
from math import *
from xraylib import AtomicWeight
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from matplotlib.ticker import AutoMinorLocator
import e_xsections

# define constants
Nav = 6.02214179e+23
PlanckConst = 6.62607004e-34
SpeedOfLight = 299792458.  # m/s
ECharge = 1.60217662e-19  # C
EMass = 9.10938356e-31  # kg
pi = 3.14159

# define pixel size in um and instrumental constants: working distance (a) in um, imaging system dimension
pixel = 1
a = 1e4
N = 1024

# beam energy in keV
energyls = [300, 1000]
nenergy = len(energyls)

# max thickness for all energies in integer nm
thickls = [int(2000), int(3000)]  # sample specific
stepls = [5, 5]

# matrix compound information
matrix = "H48.6C32.9N8.9O8.9S0.6"
matrix_den = 1.35
elements = np.array([1, 6, 7, 8, 16])
stoic = [48.6, 32.9, 8.9, 8.9, 0.6]
ineleloss = 37.5  # sample specific
mw = 0.
j = 0
for i in elements:
    mw += AtomicWeight(i) * stoic[j]
    j += 1

# prepare plots
matplotlib.rcParams['pdf.fonttype'] = 'truetype'
fontProperties = {'family': 'serif', 'serif': ['Times New Roman'], 'weight': 'normal', 'size': 12}
plt.rc('font', **fontProperties)
fig_index = 1
fig_lab = ['(a)', '(b)']
fig = plt.figure(figsize=(12, 5))

# f = open('unimatrix_data.csv','wb')

for energy in energyls:
    # eta is the out-scattered fraction (Langmore, 1992)
    eta = 1 - 4.12 / 10

    # establish table of thickness in nm (must be integer)
    thickness = thickls[fig_index - 1]
    step = stepls[fig_index - 1]
    t = range(step, thickness + 1, step)

    # compute cross sections in cm2 (see Langmore 1992)
    cs_el = e_xsections.cse_elastic(energy, elements, stoic)
    cs_inel = e_xsections.cse_inelastic(ineleloss, energy, elements, stoic)

    # probability per thickness in nm-1
    # delta is the number density of molecules (cm-3)
    delta = matrix_den / mw * Nav
    k_el = cs_el * delta * 1e-7
    k_inel = cs_inel * delta * 1e-7
    k_elin = cs_el * (1 - eta) * delta * 1e-7
    k_out = cs_el * eta * delta * 1e-7
    k_tot = k_inel + k_el

    # intensity fractions relative to primary beam
    # t is converted to cm
    i_noscat = [exp(-k_tot * t[i / step - 1]) for i in t]
    i_1el = [k_elin * t[i / step - 1] * i_noscat[i / step - 1] for i in t]
    i_pc = [sqrt(i_noscat[i / step - 1] * i_1el[i / step - 1]) for i in t]
    i_elpl = [exp(-(k_inel + k_out) * t[i / step - 1]) - i_noscat[i / step - 1] - i_1el[i / step - 1] for i in t]
    i_out = [1 - exp(-k_out * t[i / step - 1]) for i in t]
    i_inel = [exp(-k_out * t[i / step - 1]) - exp(-(k_inel + k_out) * t[i / step - 1]) for i in t]

    # output data to csv
    # writer = csv.writer(f)
    # writer.writerows([["Energy (keV)", energy], t, i_pc, i_elplpc, i_out, i_pi, i_inelpc, []])

    # run report
    print "---------------------"
    print "Energy: %.2f keV" % energy
    print "Thickness limit: %.2f nm (%.2f cm)" % (thickness, thickness / 1e7)
    print "Thickness resolution: %.2f nm" % step
    print "---------------------"


    # acquire slope arctan at specified point of a function curve


    def getatan(x, y, xdim, ydim, xsize, ysize, i0, leng):
        logy = [log10(i) for i in y]
        slope = (logy[i0 + leng] - logy[i0]) / (x[i0 + leng] - x[i0]) * xdim / ydim * ysize / xsize
        return atan(slope) / pi * 180


    def offset(x, y, xdim, ydim, xsize, ysize, i0, dist):
        logy = [log10(i) for i in y]
        aslant = getatan(x, y, xdim, ydim, xsize, ysize, i0, 80)

        class coords:
            pass

        ncoords = coords()
        ncoords.xcoord = x[i0] + dist * -sin(aslant) * xdim / xsize
        ncoords.ycoord = 10 ** (logy[i0] + dist * cos(aslant) * ydim / ysize)
        return ncoords


    # plot intensities
    minorLocator = AutoMinorLocator()
    ax = plt.subplot(ceil(nenergy / 2), 2, fig_index)
    ax.xaxis.set_minor_locator(minorLocator)
    label_pos = int(thickness / step / 3)
    npts = thickness / step
    plt.title(fig_lab[fig_index - 1] + ' Protein at ' + str(energy) + ' keV')
    xdim = thickness
    ydim = 3
    xsize = 677
    ysize = 641

    pl_pc, = plt.plot(t, i_pc, label='Phase contrast image (PCI) signal', color='blue', linewidth=2.0)
    if fig_index == 1:
        labpt = int(npts / 10 * 6)
        rot = getatan(t, i_pc, xdim, ydim, xsize, ysize, labpt, 80)
        plt.text(t[labpt], i_pc[labpt], 'Phase contrast\nimage (PCI) signal', fontsize=9, color='blue', rotation=rot)
    if fig_index == 2:
        labpt = int(npts / 10 * 5.5)
        rot = getatan(t, i_pc, xdim, ydim, xsize, ysize, labpt, 80)
        plt.text(t[labpt], i_pc[labpt], 'Phase contrast\nimage (PCI) signal', fontsize=9, color='blue', rotation=rot)

    pl_noscat, = plt.plot(t, i_noscat, label='Unscattered', color='salmon')
    if fig_index == 1:
        labpt = int(npts / 10 * 1.2)
        rot = getatan(t, i_noscat, xdim, ydim, xsize, ysize, labpt, 80)
        plt.text(t[labpt], i_noscat[labpt], 'Unscattered', fontsize=9, color='salmon', rotation=rot)
    if fig_index == 2:
        labpt = int(npts / 10 * 1.1)
        rot = getatan(t, i_noscat, xdim, ydim, xsize, ysize, labpt, 80)
        plt.text(t[labpt], i_noscat[labpt], 'Unscattered', fontsize=9, color='salmon', rotation=rot)

    pl_1el, = plt.plot(t, i_1el, label='Single elastically scattered', color='black')
    if fig_index == 1:
        labpt = int(npts / 10 * 2.9)
        rot = getatan(t, i_1el, xdim, ydim, xsize, ysize, labpt, 80)
        ofst = offset(t, i_1el, xdim, ydim, xsize, ysize, labpt, -25)
        plt.text(ofst.xcoord, ofst.ycoord, 'Single elastically scattered', fontsize=9, color='black', rotation=rot)
    if fig_index == 2:
        labpt = int(npts / 10 * 2.7)
        rot = getatan(t, i_1el, xdim, ydim, xsize, ysize, labpt, 80)
        ofst = offset(t, i_1el, xdim, ydim, xsize, ysize, labpt, -55)
        plt.text(ofst.xcoord, ofst.ycoord, 'Single elastically scattered', fontsize=9, color='black', rotation=rot)
        plt.text(ofst.xcoord, ofst.ycoord, 'Single elastically scattered', fontsize=9, color='black', rotation=rot)

    pl_elpl, = plt.plot(t, i_elpl, label='Plural elastically scattered', color='darkcyan')
    if fig_index == 1:
        labpt = int(npts / 10 * 2.3)
        rot = getatan(t, i_elpl, xdim, ydim, xsize, ysize, labpt, 80)
        ofst = offset(t, i_elpl, xdim, ydim, xsize, ysize, labpt, -45)
        plt.text(ofst.xcoord, ofst.ycoord, 'Plural elastically\nscattered', fontsize=9, color='darkcyan', rotation=rot)
    if fig_index == 2:
        labpt = int(npts / 10 * 2)
        rot = getatan(t, i_elpl, xdim, ydim, xsize, ysize, labpt, 80)
        ofst = offset(t, i_elpl, xdim, ydim, xsize, ysize, labpt, -55)
        plt.text(ofst.xcoord, ofst.ycoord, 'Plural elastically\nscattered', fontsize=9, color='darkcyan', rotation=rot)

    pl_inel, = plt.plot(t, i_inel, label='Inelastically scattered', color='magenta')
    if fig_index == 1:
        labpt = int(npts / 10 * 5.9)
        rot = getatan(t, i_inel, xdim, ydim, xsize, ysize, labpt, 80)
        ofst = offset(t, i_inel, xdim, ydim, xsize, ysize, labpt, -30)
        plt.text(ofst.xcoord, ofst.ycoord, 'Inelastically scattered', fontsize=9, color='magenta', rotation=rot)
    if fig_index == 2:
        labpt = int(npts / 10 * 6)
        rot = getatan(t, i_inel, xdim, ydim, xsize, ysize, labpt, 80)
        ofst = offset(t, i_inel, xdim, ydim, xsize, ysize, labpt, 5)
        plt.text(ofst.xcoord, ofst.ycoord, 'Inelastically scattered', fontsize=9, color='magenta', rotation=rot)

    pl_out, = plt.plot(t, i_out, label='Scattered out', color='orange')
    if fig_index == 1:
        labpt = int(npts / 10 * 6)
        rot = getatan(t, i_out, xdim, ydim, xsize, ysize, labpt, 80)
        ofst = offset(t, i_out, xdim, ydim, xsize, ysize, labpt, -60)
        plt.text(ofst.xcoord, ofst.ycoord, 'Scattered out', fontsize=9, color='orange', rotation=rot)
    if fig_index == 2:
        labpt = int(npts / 10 * 6)
        rot = getatan(t, i_out, xdim, ydim, xsize, ysize, labpt, 80)
        ofst = offset(t, i_out, xdim, ydim, xsize, ysize, labpt, -60)
        plt.text(ofst.xcoord, ofst.ycoord, 'Scattered out', fontsize=9, color='orange', rotation=rot)

    plt.semilogy()
    plt.ylim(10 ** (-ydim), 1)
    plt.xlabel('Thickness (nm)')
    plt.ylabel('Fraction')
    fig_index = fig_index + 1

fig.savefig('unimatrix_fig_e.pdf', format='pdf')
plt.show()

sys.exit()
