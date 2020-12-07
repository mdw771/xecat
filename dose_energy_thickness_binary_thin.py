#####################################
# Input units expected:             #
# Thickness: um                     #
# Energy keV                        #
# Density: g/cm3                    #
# Intermediate quantities may use   #
# different unit systems.           #
#####################################

import sys
import xraylib
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
from constants import *
from util import *
from results import *
from sample import *
from beam import *
from measurement import *
from output import *
from feat_matrix import CompositeSimulator
import h5py


energy = np.linspace(np.log10(0.1), np.log10(20), 500)
energy = 10 ** energy

t_b = np.linspace(np.log10(1), np.log10(1e4), 500)
t_b = 10 ** t_b

t_f = 0.01
# protein in ice ===============================================================
feature_den = 1.35
sample = Composite(matrix_compound='H2O', matrix_density=0.92, feature_compound='H48.6C32.9N8.9O8.9S0.6',
                   feature_density=feature_den, matrix_eloss=39.3, feature_eloss=37.5, matrix_thickness=50,
                   feature_thickness=t_f, variable='t_b', feature_alias='protein', matrix_alias='ice')

# protein in tissue ============================================================
# feature_den = 1.35
# mass_frac_water = 0.708
# mass_frac_protein = 1 - mass_frac_water
# molar_mass_water = 18.
# molar_mass_protein = 729.6
# stoic_water = 1.
# stoic_protein = molar_mass_water / mass_frac_water * mass_frac_protein / molar_mass_protein
# formula_mixture =  'H{:.5f}C{:.5f}N{:.5f}O{:.5f}S{:.5f}'.format(2 * stoic_water + 48.6 * stoic_protein,
#                                             0 * stoic_water + 32.9 * stoic_protein,
#                                             0 * stoic_water + 8.9 * stoic_protein,
#                                             1 * stoic_water + 8.9 * stoic_protein,
#                                             0 * stoic_water + 0.6 * stoic_protein,)
# print(formula_mixture)
# sample = Composite(matrix_compound=formula_mixture, matrix_density=0.92 * 0.708 + feature_den * 0.292, feature_compound='H48.6C32.9N8.9O8.9S0.6',
#                    feature_density=feature_den, matrix_eloss=39.3, feature_eloss=37.5, matrix_thickness=50,
#                    feature_thickness=t_f, variable='t_b', feature_alias='protein', matrix_alias='tissue')

# copper in silicon ==============================================================
# feature_den = 8.96
# matrix_den = 2.32
# sample = Composite(matrix_compound='Si', matrix_density=matrix_den, feature_compound='Cu',
#                    feature_density=feature_den, matrix_eloss=None, feature_eloss=None, matrix_thickness=50,
#                    feature_thickness=t_f, variable='t_b', feature_alias='Cu', matrix_alias='Si')

# ================================================================================

measure = Measurement(pixel_size=t_f, n_ccd=1024, working_distance=1e4)
simulator_thickness = CompositeSimulator(sample)
output = Output(sample, step=0.001)
output.t_b = t_b
output.t = t_b + t_f
simulator = CompositeSimulator(sample)

arr_zpc = np.zeros([t_b.size, energy.size])
arr_abs = np.zeros([t_b.size, energy.size])
n_zpc = np.zeros([t_b.size, energy.size])
n_abs = np.zeros([t_b.size, energy.size])

for i, kev in enumerate(energy):
    print('energy = {}'.format(kev))
    x_beam = XrayBeam(kev)
    i_matrix, i_feature, i_ftf, i_btb, i_btf = simulator.get_xray_categories(x_beam, measure, output, return_aux=True)
    theta_zpc_thin = simulator.get_xray_theta_thin(i_matrix, i_feature, contrast_mode='zpc')
    theta_abs_thin = simulator.get_xray_theta_thin(i_matrix, i_feature, contrast_mode='abs')
    nprobe_x_zpc = 25. / theta_zpc_thin ** 2
    nprobe_x_abs = 25. / theta_abs_thin ** 2
    # k_pi_f = i_feature.constants['k_pi_b']
    k_pi_f = i_feature.constants['k_pi_f']
    mfp = 1. / k_pi_f
    # pen_depth = mfp
    pen_depth = np.min([t_f, mfp])
    # feat_mass = t_f ** 2. * matrix_den * pen_depth * 1.e-15
    feat_mass = t_f ** 2. * feature_den * pen_depth * 1.e-15
    i_abs_f = i_ftf.i_pi.data
    dose_x_zpc = nprobe_x_zpc * kev * ECharge * 1.e3 * i_abs_f / feat_mass
    dose_x_abs = nprobe_x_abs * kev * ECharge * 1.e3 * i_abs_f / feat_mass
    arr_zpc[:, i] = dose_x_zpc
    arr_abs[:, i] = dose_x_abs
    n_zpc[:, i] = nprobe_x_zpc
    n_abs[:, i] = nprobe_x_abs

f = h5py.File('DEt_bin_thin_{}_in_{}_r{}um.h5'.format(sample.feature_alias, sample.matrix_alias, t_f), 'w')
f.create_dataset('energy_kev', dtype='float32', data=energy)
f.create_dataset('t_b_um', dtype='float32', data=t_b)
f.create_dataset('dose_zpc', dtype='float32', data=arr_zpc)
f.create_dataset('dose_abs', dtype='float32', data=arr_abs)
f.create_dataset('n_zpc', dtype='float32', data=n_zpc)
f.create_dataset('n_abs', dtype='float32', data=n_abs)

f.close()