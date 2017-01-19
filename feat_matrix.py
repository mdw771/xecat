import sys
import xraylib
import numpy as np
from math import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
import e_xsections
import numpy as np
import saxsfrac


# acquire slope arctan at specified point of a function curve
def getatan(x, y, xdim, ydim, xsize, ysize, i0, leng):
    logy = [log10(i) for i in y]
    logx = [log10(i) for i in x]
    slope = (logy[i0 + leng] - logy[i0]) / (logx[i0 + leng] - logx[i0]) * xdim / ydim * ysize / xsize
    return atan(slope) / pi * 180


def offset(x, y, xdim, ydim, xsize, ysize, i0, dist):
    logy = [log10(i) for i in y]
    logx = [log10(i) for i in x]
    aslant = getatan(x, y, xdim, ydim, xsize, ysize, i0, 80)

    class coords:
        pass

    ncoords = coords()
    ncoords.xcoord = 10 ** (logx[i0] + dist * -sin(aslant) * xdim / xsize)
    ncoords.ycoord = 10 ** (logy[i0] + dist * cos(aslant) * ydim / ysize)
    return ncoords


Nav = 6.02214179e23
PlanckConst = 6.62607004e-34
SpeedOfLight = 299792458.  # m/s
ECharge = 1.60217662e-19  # C

### sample/instrument specifics
t_f = 10  # feature thickness in integer nm
matrix = 'H2O'
feat = 'H48.6C32.9N8.9O8.9S0.6'
matrix_den = 0.92
feat_den = 1.35
elements_b = [1, 8]
stoic_b = [2, 1]
elements_f = [1, 6, 7, 8, 16]
stoic_f = [48.6, 32.9, 8.9, 8.9, 0.6]
eloss_b = 46
eloss_f = 37.5
mw_b = 0.
mw_f = 0.
j = 0
for i in elements_b:
    mw_b += xraylib.AtomicWeight(i) * stoic_b[j]
    j += 1
j = 0
for i in elements_f:
    mw_f += xraylib.AtomicWeight(i) * stoic_f[j]
    j += 1
pixel = 0.01  # pixel size in um

### prepare plot
fig = plt.figure(figsize=(6, 5))
matplotlib.rcParams['pdf.fonttype'] = 'truetype'
fontProperties = {'family': 'serif', 'serif': ['Times New Roman'], 'weight': 'normal', 'size': 12}
plt.rc('font', **fontProperties)
ax = fig.add_subplot(1, 1, 1)
ax.set_yscale('log')
ax.set_xscale('log')
ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%g'))
ax.xaxis.set_major_locator(plt.FixedLocator([0.05, 0.1, 0.3, 1, 3, 10, 30, 100, 200]))

### prepare X-ray data

# pixel size in um

# beam energy in keV
energyls = [0.5, 10]
nenergy = len(energyls)

# establish table of thickness in um (t must be integer)
t_max = 200.e3  # max total thickness in integer nm
step = 10.  # thickness step length in nm
t = np.arange(10, t_max + 1, step).astype(float)
t *= 1.e-3  # convert to um from nm
t_f *= 1.e-3
t_b = t - t_f

for energy in energyls:

    wavelen = PlanckConst * SpeedOfLight / (energy * 1000 * ECharge)  # wavelength in m

    cs_tot_b = xraylib.CS_Total_CP(matrix, energy)
    cs_pi_f = xraylib.CS_Photo_CP(feat, energy)
    cs_tot_f = xraylib.CS_Total_CP(feat, energy)

    k_tot_b = cs_tot_b * matrix_den * 1e-4
    k_pi_f = cs_pi_f * feat_den * 1e-4
    k_tot_f = cs_tot_f * feat_den * 1e-4

    # calculate contrast parameter
    ridelta_f = 1 - xraylib.Refractive_Index_Re(feat, energy, feat_den)
    ridelta_b = 1 - xraylib.Refractive_Index_Re(matrix, energy, matrix_den)


    #-----------------------------
    if energy == energyls[0]:
        ridelta_f = 0.00107625
        ridelta_b = 0.000569931
        k_tot_f = 4 * np.pi * 0.000302952 / wavelen * 1e-6
        k_tot_b = 4 * np.pi * 2.2036e-5 / wavelen * 1e-6
    else:
        ridelta_f = 2.9945e-6
        ridelta_b = 2.1281e-6
        k_tot_f = 4 * np.pi * 5.3916e-9 / wavelen * 1e-6
        k_tot_b = 4 * np.pi * 4.47286e-9 / wavelen * 1e-6
    print 'kb', k_tot_b
    print 'kf', k_tot_f
    #-----------------------------


    theta_x = 2 * np.sqrt(2) * np.pi * t_f / (wavelen * 1e6) * np.abs(ridelta_f - ridelta_b) * np.exp(-k_tot_b * t_b / 2)

    # ============================
    # 500 ev agreement test
    # beta_f = 0.000302952
    # beta_b = 2.2036e-5
    # delta_f = 0.00107625
    # delta_b = 0.000569931
    # k = 2 * np.pi / wavelen * 1e-9
    # thicknesses_background = t_b * 1e3
    # t_feature = t_f * 1e3
    # kbptp = k*beta_z*t_phase_nm
    # kbbttf = k*beta_b*(thicknesses_background-t_feature)
    # kdft = k*delta_f*t_feature
    # kdbt = k*delta_b*t_feature
    # kbft = k*beta_f*t_feature
    # kbbt = k*beta_b*t_feature
    # etafbt = k*(delta_f-delta_b)*t_feature
    # if_zernike = exp(-kbbt)*(1+exp(-kbptp))+exp(-kbft)+ \
    #              2.*exp(-0.5*kbft)*exp(-0.5*kbbt)*exp(-0.5*kbptp)*sin(etafbt)- \
    #              2.*exp(-0.5*kbft)*exp(-0.5*kbbt)*cos(etafbt)
    # if_zernike = if_zernike*exp(-kbbttf)
    # ib_zernike = exp(-kbbt)*exp(-kbptp)
    # ib_zernike = ib_zernike*exp(-kbbttf)
    # theta_phase = abs(if_zernike-ib_zernike)/sqrt(if_zernike+ib_zernike)
    # ============================

    # calculate dose in Gray (J/kg)
    nprobe_x = 25. / theta_x ** 2
    # --------
    print 'nprobe_x at t = 100 nm: ' + str(nprobe_x[10])
    if energy is energyls[0]:
        nprobe_x[10] = 7.6843e4
    else:
        nprobe_x[10] = 6.48207e+07
    # --------

    i_abs_f = np.exp(-k_tot_b * t_b / 2.) * (1. - np.exp(-k_tot_f * t_f))
    print 'i_abs_f', i_abs_f[10]
    # affected feature mass in kg
    feat_mass = pixel ** 2. * feat_den * t_f * 1.e-15
    dose_x = nprobe_x * energy * ECharge * 1.e3 * i_abs_f / feat_mass
    print 'inv_feat_mass', str(1/feat_mass)
    #dose_x = nprobe_x * energy * 1000 * k_tot_f / (pixel ** 2 * 1e-15 * feat_den) * ECharge
    print dose_x[10]
    print '=============='

    # plot
    labpt = int(len(t) / 10 * 0.02)
    rot = getatan(t, dose_x, log10(200) - log10(0.05), 7, 6, 5, labpt, 50)
    if energy == energyls[0]:
        ofst = offset(t, dose_x, log10(200) - log10(0.05), 7, 6, 5, labpt, -0.2)
        plt.plot(t, dose_x, color='red')
        plt.text(ofst.xcoord, ofst.ycoord, str(energy) + ' keV soft X-ray', fontsize=9, color='red', rotation=rot)
    else:
        ofst = offset(t, dose_x, log10(200) - log10(0.05), 7, 6, 5, labpt, 0.05)
        plt.plot(t, dose_x, color='green')
        plt.text(ofst.xcoord, ofst.ycoord, str(energy) + ' keV soft X-ray', fontsize=9, color='green', rotation=rot)
    plt.hold(True)


### prepare electron data
# eta is the out-scattered fraction
eta = 1 - 4.12 / 10

# electron energy in keV
energy = 300

# establish table of thickness in um (t must be integer)
t_max = int(10e3)  # max total thickness in integer nm
step = 10.  # thickness step length in nm
t = np.arange(step, t_max + 1, step)
t = t.astype(float)
t *= 1.e-3  # convert nm to um
t_b = t - t_f

# compute cross sections in cm2 (see Langmore 1992)
cs_el_b = e_xsections.cse_elastic(energy, elements_b, stoic_b)
cs_inel_b = e_xsections.cse_inelastic(eloss_b, energy, elements_b, stoic_b)
cs_el_f = e_xsections.cse_elastic(energy, elements_f, stoic_f)
cs_inel_f = e_xsections.cse_inelastic(eloss_f, energy, elements_f, stoic_f)

# probability per thickness in micron-1
delta_b = matrix_den / mw_b * Nav
delta_f = feat_den / mw_f * Nav
k_el_b = cs_el_b * delta_b * 1.e-4
k_inel_b = cs_inel_b * delta_b * 1.e-4
k_elin_b = cs_el_b * (1 - eta) * delta_b * 1.e-4
k_out_b = cs_el_b * eta * delta_b * 1.e-4
k_tot_b = k_inel_b + k_el_b
k_el_f = cs_el_f * delta_f * 1.e-4
k_inel_f = cs_inel_f * delta_f * 1.e-4
k_elin_f = cs_el_f * (1 - eta) * delta_f * 1.e-4
k_out_f = cs_el_f * eta * delta_f * 1.e-4
k_tot_f = k_inel_f + k_el_f

# intensity fractions for matrix relative to primary beam
i_noscat_b = np.exp(-k_tot_b * t)
i_1el_b = k_elin_b * t * i_noscat_b
i_innoinel_b = np.exp(-(k_inel_b + k_out_b) * t)
i_elpl_b = i_innoinel_b - i_noscat_b - i_1el_b
i_in_b = np.exp(-k_out_b * t)
i_inel_b = i_in_b - i_innoinel_b

# intensity fractions for feature relative to primary beam
i_noscat_f = np.exp(-k_tot_b * t_b) * np.exp(-k_tot_f * t_f)
i_1el_f = (k_elin_b * t_b + k_elin_f * t_f) * i_noscat_f
i_1elf_f = k_elin_f * t_f * i_noscat_f
i_innoinel_f = np.exp(-(k_out_b + k_inel_b) * t_b) * np.exp(-(k_out_f + k_inel_f) * t_f)
i_elpl_f = i_innoinel_f - i_noscat_f - i_1el_f
i_in_f = np.exp(-k_out_b * t_b) * np.exp(-k_out_f * t_f)
i_inel_f = i_in_f - i_innoinel_f

# calculate contrast parameters
theta_e = np.abs(i_innoinel_f - i_innoinel_b) + 2 * np.sqrt(i_noscat_f * i_1elf_f)
theta_e = theta_e / np.sqrt(i_in_f + i_in_b)
theta_ef = np.abs(i_innoinel_f - i_innoinel_b) + 2 * np.sqrt(i_noscat_f * i_1elf_f)
theta_ef = theta_ef / np.sqrt(i_innoinel_f + i_innoinel_b)

# calculate dose
nprobe_e = 25. / theta_e ** 2
nprobe_ef = 25. / theta_ef ** 2
dose_e = nprobe_e * eloss_f * k_inel_f / (pixel ** 2 * 1e-12 * feat_den) * ECharge * 1e3
dose_ef = nprobe_ef * eloss_f * k_inel_f / (pixel ** 2 * 1e-12 * feat_den) * ECharge * 1e3

# plot
# no energy filter
labpt = int(len(t) / 10 * 1.5)
rot = getatan(t, dose_e, log10(10) - log10(0.05), 7, 6, 5, labpt, 100)
ofst = offset(t, dose_e, log10(10) - log10(0.05), 7, 6, 5, labpt, 0)
plt.plot(t, dose_e, color='gray')
plt.text(ofst.xcoord, ofst.ycoord, '(no energy filter)', fontsize=9, color='grey', rotation=rot)
# energy filter
labpt = int(len(t) / 10 * 1)
rot = getatan(t, dose_ef, log10(10) - log10(0.05), 7, 6, 5, labpt, 100)
ofst = offset(t, dose_ef, log10(10) - log10(0.05), 7, 6, 5, labpt, 0)
plt.plot(t, dose_ef, color='black')
plt.text(ofst.xcoord, ofst.ycoord, '(energy filter)', fontsize=9, color='black', rotation=rot)

plt.hold(False)
plt.axis([0.05, 200, 1e5, 1e12])
plt.xlabel('Thickness ($\mu$m)')
plt.ylabel('Dose (Gray)')

plt.show()
fig.savefig('dose.pdf', type='pdf')
