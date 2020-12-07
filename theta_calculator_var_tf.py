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


measure = Measurement(pixel_size=1, n_ccd=1024, working_distance=1e4)

matplotlib.rcParams['pdf.fonttype'] = 'truetype'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams['font.weight'] = 'normal'
plt.rcParams['font.size'] = 12
# fontProperties = {'family': 'sans-serif', 'sans-serif': ['Helvetica'], 'weight': 'normal', 'size': 12}
# plt.rc('font', **fontProperties)

fig = plt.figure(figsize=(12, 5))

# fix t_b, vary t_f, zpc
ax = fig.add_subplot(1, 2, 1)

# =============================================================================================
energy = 15
t_f_max = 1000
t_b = 0
step = 0.1
sample = Composite(matrix_compound='H2O', matrix_density=0.92, feature_compound='H48.6C32.9N8.9O8.9S0.6',
                   feature_density=1.35, matrix_eloss=39.3, feature_eloss=37.5, matrix_thickness=t_b,
                   feature_thickness=t_f_max, variable='t_f', feature_alias='protein', matrix_alias='ice')
# =============================================================================================

simulator = CompositeSimulator(sample)
output = Output(sample, step=step)
x_beam = XrayBeam(energy)
i_matrix, i_feature, i_ftf, i_btb, i_btf = simulator.get_xray_categories(x_beam, measure, output, return_aux=True)
theta_zpc_complete = simulator.get_xray_theta_complete(i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='zpc', numerical=False)
theta_zpc_simple = simulator.get_xray_theta_simple(i_matrix, i_feature, contrast_mode='zpc')
theta_zpc_thin = simulator.get_xray_theta_thin(i_matrix, i_feature, contrast_mode='zpc')

ax.plot(output.t_f, theta_zpc_complete, label='ZPC complete model')
ax.plot(output.t_f, theta_zpc_simple, label='ZPC conventional model', linestyle='--', linewidth=2)
ax.plot(output.t_f, theta_zpc_thin, label='ZPC pure-phase model')
plt.xlabel('$t_b$ ($\mu$m)')
plt.ylabel('Contrast parameter $\Theta$')
plt.legend(frameon=False)

plt.xlim([0, t_f_max])
plt.ylim([0, 7])

# ax = fig.add_subplot(1, 2, 2)
#
# theta_abs_complete = simulator.get_xray_theta_complete(i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='abs')
# theta_abs_simple = simulator.get_xray_theta_simple(i_matrix, i_feature, contrast_mode='abs')
# theta_abs_thin = simulator.get_xray_theta_thin(i_matrix, i_feature, contrast_mode='abs')
#
# ax.plot(output.t_b, theta_abs_complete, label='Absorption complete model')
# ax.plot(output.t_b, theta_abs_simple, label='Absorption conventional model')
# plt.xlabel('$t_b$ ($\mu$m)')
# plt.ylabel('Contrast parameter $\Theta$')
# plt.xlim([0, t_f_max])
# # plt.ylim([0, 0.02])
plt.legend(frameon=False)

plt.savefig('zpc_tf{}_tb{}_{}kev_var_tf.pdf'.format(t_f_max, t_b, energy), format='pdf')
# plt.show()