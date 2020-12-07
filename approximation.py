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

energy = 0.5
t_f = 0.02
t_b_max = 40
step = 0.01
sample = Composite(matrix_compound='H2O', matrix_density=0.92,
                   feature_compound='H48.6C32.9N8.9O8.9S0.6', feature_density=1.35,
                   matrix_eloss=39.3, feature_eloss=37.5,
                   matrix_thickness=t_b_max, feature_thickness=t_f,
                   variable='t_b')
simulator = CompositeSimulator(sample)
output = Output(sample, step=step)
x_beam = XrayBeam(energy)
i_matrix, i_feature, i_ftf, i_btb, i_btf = simulator.get_xray_categories(x_beam, measure, output, return_aux=True)
theta_zpc_complete = simulator.get_xray_theta_complete(i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='zpc')
theta_zpc_simple = simulator.get_xray_theta_simple(i_matrix, i_feature, contrast_mode='zpc')
theta_zpc_thin = simulator.get_xray_theta_thin(i_matrix, i_feature, contrast_mode='zpc')

ax.plot(output.t_b, theta_zpc_complete, label='ZPC complete expression')
ax.plot(output.t_b, theta_zpc_simple, label='ZPC simple expression')
ax.plot(output.t_b, theta_zpc_thin, label='ZPC pure-phase thin sample')
plt.xlabel('$t_b$ ($\mu$m)')
plt.ylabel('Contrast parameter $\Theta$')
plt.xlim([0, t_b_max])
plt.ylim([0, 0.05])
plt.legend(frameon=False)

ax = fig.add_subplot(1, 2, 2)

theta_abs_complete = simulator.get_xray_theta_complete(i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='abs')
theta_abs_simple = simulator.get_xray_theta_simple(i_matrix, i_feature, contrast_mode='abs')
theta_abs_thin = simulator.get_xray_theta_thin(i_matrix, i_feature, contrast_mode='abs')

ax.plot(output.t_b, theta_abs_complete, label='Absorption complete expression')
ax.plot(output.t_b, theta_abs_simple, label='Absorption simple expression')
plt.xlabel('$t_b$ ($\mu$m)')
plt.ylabel('Contrast parameter $\Theta$')
plt.xlim([0, t_b_max])
plt.ylim([0, 0.02])
plt.legend(frameon=False)

plt.savefig('zpc_abs_tf0.02_tb40v_0.5kev.pdf', format='pdf')

# ============================================================================

fig = plt.figure(figsize=(12, 5))
ax = fig.add_subplot(1, 2, 1)

energy = 15
t_b_max = 1000
t_f = 0.02
step = 0.01
sample = Composite(matrix_compound='H2O', matrix_density=0.92,
                   feature_compound='H48.6C32.9N8.9O8.9S0.6', feature_density=1.35,
                   matrix_eloss=39.3, feature_eloss=37.5,
                   matrix_thickness=t_b_max, feature_thickness=t_f,
                   variable='t_b')
simulator = CompositeSimulator(sample)
output = Output(sample, step=step)
x_beam = XrayBeam(energy)
i_matrix, i_feature, i_ftf, i_btb, i_btf = simulator.get_xray_categories(x_beam, measure, output, return_aux=True)
theta_zpc_complete = simulator.get_xray_theta_complete(i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='zpc')
theta_zpc_simple = simulator.get_xray_theta_simple(i_matrix, i_feature, contrast_mode='zpc')
# theta_zpc_thin = simulator.get_xray_theta_thin(i_matrix, i_feature, contrast_mode='zpc')
ax.plot(output.t_b, theta_zpc_complete, label='ZPC complete expression')
ax.plot(output.t_b, theta_zpc_simple, label='ZPC simple expression')
# ax.plot(output.t_b, theta_zpc_thin, label='ZPC pure-phase thin sample')
plt.xlabel('$t_b$ ($\mu$m)')
plt.ylabel('$Contrast parameter \Theta$')
plt.xlim([0, t_b_max])
plt.ylim([0.00074, 0.00084])
plt.legend(frameon=False)

ax = fig.add_subplot(1, 2, 2)
t_f = 1000
sample = Composite(matrix_compound='H2O', matrix_density=0.92,
                   feature_compound='H48.6C32.9N8.9O8.9S0.6', feature_density=1.35,
                   matrix_eloss=39.3, feature_eloss=37.5,
                   matrix_thickness=t_b_max, feature_thickness=t_f,
                   variable='t_b')
simulator = CompositeSimulator(sample)
output = Output(sample, step=step)
x_beam = XrayBeam(energy)
i_matrix, i_feature, i_ftf, i_btb, i_btf = simulator.get_xray_categories(x_beam, measure, output, return_aux=True)
theta_zpc_complete = simulator.get_xray_theta_complete(i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='zpc')
theta_zpc_simple = simulator.get_xray_theta_simple(i_matrix, i_feature, contrast_mode='zpc')
# theta_zpc_thin = simulator.get_xray_theta_thin(i_matrix, i_feature, contrast_mode='zpc')
ax.plot(output.t_b, theta_zpc_complete, label='ZPC complete expression')
ax.plot(output.t_b, theta_zpc_simple, label='ZPC simple expression')
# ax.plot(output.t_b, theta_zpc_thin, label='ZPC pure-phase thin sample')
plt.xlabel('$t_b$ ($\mu$m)')
plt.ylabel('$\Theta$')
plt.xlim([0, t_b_max])
plt.ylim([6.4, 7.1])
plt.legend(frameon=False)

plt.savefig('zpc_tf0.02_tf1000_tb1000v_15kev.pdf', format='pdf')

# ============================================================================

fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(1, 1, 1)

energy = 15
t_b_max = 1000
t_f_max = 1000
step = 0.01
sample = Composite(matrix_compound='H2O', matrix_density=0.92,
                   feature_compound='H48.6C32.9N8.9O8.9S0.6', feature_density=1.35,
                   matrix_eloss=39.3, feature_eloss=37.5,
                   matrix_thickness=t_b_max, feature_thickness=t_f_max,
                   variable='both')
simulator = CompositeSimulator(sample)
output = Output(sample, step=step)
x_beam = XrayBeam(energy)

i_matrix, i_feature, i_ftf, i_btb, i_btf = simulator.get_xray_categories(x_beam, measure, output, return_aux=True)
theta_zpc_complete = simulator.get_xray_theta_complete(i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='zpc')
theta_zpc_simple = simulator.get_xray_theta_simple(i_matrix, i_feature, contrast_mode='zpc')
theta_zpc_thin = simulator.get_xray_theta_thin(i_matrix, i_feature, contrast_mode='zpc')
ax.plot(output.t_b, theta_zpc_complete, label='ZPC complete expression')
ax.plot(output.t_b, theta_zpc_simple, label='ZPC simple expression')
ax.plot(output.t_b, theta_zpc_thin, label='ZPC pure-phase thin sample')
plt.xlabel('$t_b$ or $t_f$ ($\mu$m)')
plt.ylabel('$Contrast parameter \Theta$')
plt.xlim([0, t_b_max])
plt.ylim([0, 7])
plt.legend(frameon=False)

plt.savefig('zpc_tf1000v_tb1000v_15kev.pdf', format='pdf')

