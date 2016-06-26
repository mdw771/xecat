import sys
import xraylib
from math import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.integrate import quad
import e_xsections

Nav = 6.02214179e+23
PlanckConst = 6.62607004e-34
SpeedOfLight = 299792458.  # m/s
ECharge = 1.60217662e-19  # C
EMass = 9.10938356e-31  # kg
pi = 3.14159

### sample specifics
t_f = 10  # feature thickness in integer nm
matrix = 'H2O'
feat = 'H48.6C32.9N8.9O8.9S0.6'
matrix_den = 0.92
feat_den = 1.35
elements_b = [1, 8]
stoic_b = [2, 1]
elements_f = [1, 6, 7, 8, 16]
stoic_f = [48.6, 32.9, 8.9, 8.9, 0.6]
eloss_b = 39.3
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

### prepare plot
fig = plt.figure(figsize=(6, 5))
matplotlib.rcParams['pdf.fonttype'] = 'truetype'
fontProperties = {'family': 'serif', 'serif': ['Times New Roman'], 'weight': 'normal', 'size': 12}
plt.rc('font', **fontProperties)
ax = fig.add_subplot(1,1,1)
ax.set_yscale('log')
ax.set_xscale('log')
ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%g'))
ax.xaxis.set_major_locator(plt.FixedLocator([0.05, 0.1, 0.3, 1, 3, 10, 30, 100, 200]))

### prepare X-ray data

# define pixel size in um and instrumental constants: working distance (a) in um, imaging system dimension
pixel = 1
a = 1e4
N = 1024

# beam energy in keV
energyls = [0.5, 10]
nenergy = len(energyls)

# establish table of thickness in um (t must be integer)
t_max = int(200e3)  # max total thickness in integer nm
step = 100  # thickness step length in nm
t = range(5, t_max + 1, step)
t_index = range(0, len(t), 1)
t = [t[i] * 1e-3 for i in t_index]  # convert to um from nm
t_f *= 1e-3
t_b = [t[i] - t_f for i in t_index]


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


for energy in energyls:

    # eta is reversely scattered fraction
    # NA based on typical resolution of 1 um
    # Abbe criterion applied
    eta = 0.5
    wavelen = PlanckConst * SpeedOfLight / (energy * 1000 * ECharge)
    NA = (0.5 * wavelen) / (pixel * 1e-6)

    # frac1 is the fraction in FORWARD SCATTERED (detectable) photons that enter the aperture and contribute to PCI signal
    # s in nm^-1
    s_max = 4 * pi * sin(pi / 4) / (wavelen * 1e9)
    s_glim = 1.3 / 4.1888
    s_na = 4 * pi * sin(NA / 2) / (wavelen * 1e9)
    if s_na > s_max:
        print 'WARNING: Check NA.'
        sys.exit()


    def guinier(s):
        return exp(-1.396 * s ** 2)


    def higha(s):
        return 0.0168 * s ** -4 - 0.1514 * s ** -3 + 0.4397 * s ** -2 - 0.2441 * s ** -1 + 0.0399


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

    # frac1m is the analog of frac1 for small molecule matrix
    frac1m = NA / 1.414

    # frac2 is the fraction of FORWARD scattered inelastic scattered or plural scattered photons that enter the aperture
    # and contribute to PCI background
    frac2 = NA ** 2 / 2

    # retrieve cross sections in cm2/g
    cs_inel_b = xraylib.CS_Compt_CP(matrix, energy)
    cs_el_b = xraylib.CS_Rayl_CP(matrix, energy)
    cs_pi_b = xraylib.CS_Photo_CP(matrix, energy)
    cs_tot_b = xraylib.CS_Total_CP(matrix, energy)
    cs_inel_f = xraylib.CS_Compt_CP(feat, energy)
    cs_el_f = xraylib.CS_Rayl_CP(feat, energy)
    cs_pi_f = xraylib.CS_Photo_CP(feat, energy)
    cs_tot_f = xraylib.CS_Total_CP(feat, energy)

    # probability per thickness in micron-1
    k_el_b = cs_el_b * matrix_den * 1e-4
    k_inel_b = cs_inel_b * matrix_den * 1e-4
    k_elin_b = cs_el_b * (1 - eta) * matrix_den * 1e-4
    k_inelin_b = cs_inel_b * (1 - eta) * matrix_den * 1e-4
    k_out_b = cs_el_b * eta * matrix_den * 1e-4 + cs_inel_b * eta * matrix_den * 1e-4
    k_pi_b = cs_pi_b * matrix_den * 1e-4
    k_tot_b = k_inel_b + k_el_b + k_pi_b
    k_el_f = cs_el_f * feat_den * 1e-4
    k_inel_f = cs_inel_f * feat_den * 1e-4
    k_elin_f = cs_el_f * (1 - eta) * feat_den * 1e-4
    k_inelin_f = cs_inel_f * (1 - eta) * feat_den * 1e-4
    k_out_f = cs_el_f * eta * feat_den * 1e-4 + cs_inel_f * eta * feat_den * 1e-4
    k_pi_f = cs_pi_f * feat_den * 1e-4
    k_tot_f = k_inel_f + k_el_f + k_pi_f

    # intensity fractions for matrix relative to primary beam
    i_noscat_b = [exp(-k_tot_b * t[i]) for i in t_index]
    i_1el_b = [k_elin_b * t[i] * i_noscat_b[i] for i in t_index]
    i_1elpc_b = [frac1m * i_1el_b[i] for i in t_index]
    i_pc_b = [sqrt(i_noscat_b[i] * i_1elpc_b[i]) for i in t_index]
    i_innoinel_b = [exp(-(k_inelin_b + k_out_b + k_pi_b) * t[i]) for i in t_index]
    i_elpl_b = [i_innoinel_b[i] - i_noscat_b[i] - i_1el_b[i] for i in t_index]
    i_elplpc_b = [i_elpl_b[i] * frac2 for i in t_index]
    i_out_b = [(k_out_b / (k_out_b + k_pi_b)) * (1 - exp(-(k_out_b + k_pi_b) * t[i])) for i in t_index]
    i_pi_b = [(k_pi_b / (k_out_b + k_pi_b)) * (1 - exp(-(k_out_b + k_pi_b) * t[i])) for i in t_index]
    i_in_b = [exp(-(k_out_b + k_pi_b)) for i in t_index]
    i_inel_b = [1 - i_out_b[i] - i_pi_b[i] - i_innoinel_b[i] for i in t_index]
    i_inelpc_b = [i_inel_b[i] * frac2 for i in t_index]
    i_innoinelpc_b = [i_noscat_b[i] + i_1elpc_b[i] + i_elplpc_b[i] for i in t_index]
    i_inpc_b = [i_innoinelpc_b[i] + i_inelpc_b[i] for i in t_index]

    # intensity fractions for feature relative to primary beam
    i_noscat_f = [exp(-k_tot_b * t_b[i]) * exp(-k_tot_f * t_f) for i in t_index]
    i_1el_f = [(k_elin_b * t_b[i] + k_elin_f * t_f) * i_noscat_f[i] for i in t_index]
    i_1elpc_f = [frac1 * i_1el_f[i] for i in t_index]
    i_1elf_f = [k_elin_f * t_f * i_noscat_f[i] for i in t_index]
    i_1elfpc_f = [frac1 * i_1elf_f[i] for i in t_index]
    i_innoinel_f = [exp(-(k_out_b + k_inelin_b + k_pi_b) * t_b[i]) * exp(-(k_out_f + k_inelin_f + k_pi_f) * t_f) for i
                    in t_index]
    i_elpl_f = [i_innoinel_f[i] - i_noscat_f[i] - i_1el_f[i] for i in t_index]
    i_elplpc_f = [i_elpl_f[i] * frac2 for i in t_index]
    i_out_f = [(1 - i_out_b[i] - i_pi_b[i]) * (k_out_f / (k_out_f + k_pi_f)) * (1 - exp(-(k_out_f + k_pi_f) * t_f)) for
               i in t_index]
    i_pi_f = [(1 - i_out_b[i] - i_pi_b[i]) * (k_pi_f / (k_out_f + k_pi_f)) * (1 - exp(-(k_out_f + k_pi_f) * t_f)) for i
              in t_index]
    i_in_f = [exp(-(k_out_b + k_pi_b) * t_b[i]) * exp(-(k_out_f + k_pi_f) * t_f) for i in t_index]
    i_inel_f = [1 - i_out_f[i] - i_pi_f[i] - i_innoinel_f[i] for i in t_index]
    i_inelpc_f = [i_inel_f[i] * frac2 for i in t_index]
    i_innoinelpc_f = [i_noscat_f[i] + i_1elpc_f[i] + i_elplpc_f[i] for i in t_index]
    i_inpc_f = [i_innoinelpc_f[i] + i_inelpc_f[i] for i in t_index]

    # calculate contrast parameter
    # theta_x = [abs(i_innoinelpc_f[i] - i_innoinelpc_b[i]) + 2 * sqrt(i_noscat_f[i] * i_1elfpc_f[i]) for i in t_index]
    theta_x = [abs(i_innoinelpc_f[i] - i_innoinelpc_b[i]) + 2 * sqrt(i_noscat_f[i] * i_1elf_f[i]) for i in t_index]
    theta_x = [theta_x[i] / sqrt(i_inpc_f[i] + i_inpc_b[i]) for i in t_index]

    # calculate dose in Gray (J/kg)
    nprobe_x = [25 / theta_x[i] ** 2 for i in t_index]
    dose_x = [nprobe_x[i] * energy * k_pi_f / (pixel ** 2 * 1e-12 * feat_den) * ECharge * 1e6 for i in t_index]

    # plot
    labpt = int(len(t) / 10 * 0.02)
    rot = getatan(t, dose_x, log10(200) - log10(0.05), 7, 6, 5, labpt, 50)
    print rot
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
step = 5  # thickness step length in nm
t = range(step, t_max + 1, step)
t_index = range(0, len(t), 1)
t = [t[i] * 1e-3 for i in t_index]  # convert to um from nm
t_f *= 1e-3
t_b = [t[i] - t_f for i in t_index]

# compute cross sections in cm2 (see Langmore 1992)
cs_el_b = e_xsections.cse_elastic(energy, elements_b, stoic_b)
cs_inel_b = e_xsections.cse_inelastic(eloss_b, energy, elements_b, stoic_b)
cs_el_f = e_xsections.cse_elastic(energy, elements_f, stoic_f)
cs_inel_f = e_xsections.cse_inelastic(eloss_f, energy, elements_f, stoic_f)

# probability per thickness in micron-1
delta_b = matrix_den / mw_b * Nav
delta_f = feat_den / mw_f * Nav
k_el_b = cs_el_b * delta_b * 1e-4
k_inel_b = cs_inel_b * delta_b * 1e-4
k_elin_b = cs_el_b * (1 - eta) * delta_b * 1e-4
k_inelin_b = cs_inel_b * (1 - eta) * delta_b * 1e-4
k_out_b = cs_el_b * eta * delta_b * 1e-4 + cs_inel_b * eta * delta_b * 1e-4
k_tot_b = k_inel_b + k_el_b
k_el_f = cs_el_f * delta_f * 1e-4
k_inel_f = cs_inel_f * delta_f * 1e-4
k_elin_f = cs_el_f * (1 - eta) * delta_f * 1e-4
k_inelin_f = cs_inel_f * (1 - eta) * delta_f * 1e-4
k_out_f = cs_el_f * eta * delta_f * 1e-4 + cs_inel_f * eta * delta_f * 1e-4
k_tot_f = k_inel_f + k_el_f

# intensity fractions for matrix relative to primary beam
i_noscat_b = [exp(-k_tot_b * t[i]) for i in t_index]
i_1el_b = [k_elin_b * t[i] * i_noscat_b[i] for i in t_index]
i_innoinel_b = [exp(-(k_inelin_b + k_out_b) * t[i]) for i in t_index]
i_elpl_b = [i_innoinel_b[i] - i_noscat_b[i] - i_1el_b[i] for i in t_index]
i_in_b = [exp(-k_out_b * t[i]) for i in t_index]
i_inel_b = [i_in_b[i] - i_innoinel_b[i] for i in t_index]

# intensity fractions for feature relative to primary beam
i_noscat_f = [exp(-k_tot_b * t_b[i]) * exp(-k_tot_f * t_f) for i in t_index]
i_1el_f = [(k_elin_b * t_b[i] + k_elin_f * t_f) * i_noscat_f[i] for i in t_index]
i_1elf_f = [k_elin_f * t_f * i_noscat_f[i] for i in t_index]
i_innoinel_f = [exp(-(k_out_b + k_inelin_b) * t_b[i]) * exp(-(k_out_f + k_inelin_f) * t_f) for i in t_index]
i_elpl_f = [i_innoinel_f[i] - i_noscat_f[i] - i_1el_f[i] for i in t_index]
i_in_f = [exp(-k_out_b * t_b[i]) * exp(-k_out_f * t_f) for i in t_index]
i_inel_f = [i_in_f[i] - i_innoinel_f[i] for i in t_index]

# calculate contrast parameters
theta_e = [(abs(i_innoinel_f[i] - i_innoinel_b[i]) + 2 * sqrt(i_noscat_f[i] * i_1elf_f[i])) for i in t_index]
theta_e = [theta_e[i] / sqrt(i_in_f[i] + i_in_b[i]) for i in t_index]
theta_ef = [(abs(i_innoinel_f[i] - i_innoinel_b[i]) + 2 * sqrt(i_noscat_f[i] * i_1elf_f[i])) for i in t_index]
theta_ef = [theta_ef[i] / sqrt(i_innoinel_f[i] + i_innoinel_b[i]) for i in t_index]

# calculate dose
nprobe_e = [25 / theta_e[i] ** 2 for i in t_index]
nprobe_ef = [25 / theta_ef[i] ** 2 for i in t_index]
dose_e = [nprobe_e[i] * eloss_f * k_inel_f / (pixel ** 2 * 1e-12 * feat_den) * ECharge * 1e3 for i in t_index]
dose_ef = [nprobe_ef[i] * eloss_f * k_inel_f / (pixel ** 2 * 1e-12 * feat_den) * ECharge * 1e3 for i in t_index]

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
plt.axis([0.05, 200, 1e3, 1e10])
plt.xlabel('Thickness ($\mu$m)')
plt.ylabel('Dose (Gray)')

plt.show()
fig.savefig('dose.pdf', type='pdf')
