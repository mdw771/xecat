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


energy = np.linspace(np.log10(5), np.log10(100), 500)
energy = 10 ** energy

t_b = np.linspace(np.log10(1000), np.log10(100e3), 500)
t_b = 10 ** t_b

t_f = 0.1

feature_den = 1.35

sample = Composite(matrix_compound='H2O', matrix_density=0.92, feature_compound='H48.6C32.9N8.9O8.9S0.6',
                   feature_density=feature_den, matrix_eloss=39.3, feature_eloss=37.5, matrix_thickness=50,
                   feature_thickness=t_f, variable='t_b')
measure = Measurement(pixel_size=0.1, n_ccd=1024, working_distance=1e4)
simulator_thickness = CompositeSimulator(sample)
output = Output(sample, step=0.001)
output.t_b = t_b
output.t = t_b + t_f
simulator = CompositeSimulator(sample)

arr_zpc = np.zeros([t_b.size, energy.size])
arr_abs = np.zeros([t_b.size, energy.size])

for i, kev in enumerate(energy):
    print('energy = {}'.format(kev))
    x_beam = XrayBeam(kev)
    i_matrix, i_feature, i_ftf, i_btb, i_btf = simulator.get_xray_categories(x_beam, measure, output, return_aux=True)
    theta_zpc_complete = simulator.get_xray_theta_complete(i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='zpc')
    theta_abs_complete = simulator.get_xray_theta_complete(i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='abs')
    nprobe_x_zpc = 25. / theta_zpc_complete ** 2
    nprobe_x_abs = 25. / theta_abs_complete ** 2
    k_pi_f = i_feature.constants['k_pi_f']
    mfp = 1. / k_pi_f
    pen_depth = np.min([t_f, mfp])
    feat_mass = t_f ** 2. * feature_den * pen_depth * 1.e-15
    i_abs_f = i_feature.i_pi.data
    dose_x_zpc = nprobe_x_zpc * kev * ECharge * 1.e3 * i_abs_f / feat_mass
    dose_x_abs = nprobe_x_abs * kev * ECharge * 1.e3 * i_abs_f / feat_mass
    arr_zpc[:, i] = dose_x_zpc
    arr_abs[:, i] = dose_x_abs

f = h5py.File('dose_energy_thickness_limit.h5', 'w')
f.create_dataset('energy_kev', dtype='float32', data=energy)
f.create_dataset('t_b_um', dtype='float32', data=t_b)
f.create_dataset('dose_zpc', dtype='float32', data=arr_zpc)
f.create_dataset('dose_abs', dtype='float32', data=arr_abs)

import dxchange
dxchange.write_tiff(np.log10(arr_zpc), 'tmp/arr_zpc', dtype='float32', overwrite=True)

f.close()