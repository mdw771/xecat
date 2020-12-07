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


def get_xray_theta_complete(sample_protein, sample_ice, i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='zpc'):
    (delta_f, _), (delta_bb, _) = sample_protein.get_ri_delta_beta(i_matrix.beam)
    (delta_b, _), (_, _) = sample_ice.get_ri_delta_beta(i_matrix.beam)
    # (delta_f, _), (delta_b, _) = self.sample.get_ri_delta_beta_henke(i_matrix.beam)
    # mu_f, mu_b = 4 * np.pi * np.array([beta_f, beta_b]) / (i_matrix.beam.wavelength * 1e6)
    mu_bb = i_matrix.constants['k_pi_b']
    mu_f = i_feature.constants['k_pi_f']
    mu_b = i_matrix.constants['k_pi_f']
    eta_f, eta_b, eta_bb = 2 * np.pi * np.array([delta_f, delta_b, delta_bb]) / (i_matrix.beam.wavelength * 1e6)
    t_f = i_matrix.output.t_f
    t_b = i_matrix.output.t_b

    i_noscat_ftf = i_ftf.i_noscat.data
    i_1el_ftf = i_ftf.i_1el.data
    i_elpl_ftf = i_ftf.i_elpl.data
    i_inel_ftf = i_ftf.i_inel.data
    i_out_ftf = i_ftf.i_out.data
    i_abs_ftf = i_ftf.i_pi.data
    fftf = np.sqrt((1 - i_abs_ftf - i_elpl_ftf - i_inel_ftf - i_out_ftf) / (1 - i_abs_ftf))

    i_noscat_btb = i_btb.i_noscat.data
    i_1el_btb = i_btb.i_1el.data
    i_elpl_btb = i_btb.i_elpl.data
    i_inel_btb = i_btb.i_inel.data
    i_out_btb = i_btb.i_out.data
    i_abs_btb = i_btb.i_pi.data
    # fbtb = np.sqrt((1 - i_abs_btb - i_elpl_btb - i_inel_btb - i_out_btb) / (1 - i_abs_btb))
    fbtb = np.sqrt((i_noscat_btb + i_1el_btb) / (1 - i_abs_btb))

    i_noscat_btf = i_btf.i_noscat.data
    i_1el_btf = i_btf.i_1el.data
    i_elpl_btf = i_btf.i_elpl.data
    i_inel_btf = i_btf.i_inel.data
    i_out_btf = i_btf.i_out.data
    i_abs_btf = i_btf.i_pi.data
    fbtf = np.sqrt((1 - i_abs_btf - i_elpl_btf - i_inel_btf - i_out_btf) / (1 - i_abs_btf))

    i_signal_f = (1 + 2 * (eta_f - eta_b) * t_f) * np.exp(-mu_b * t_b)
    i_signal_b = np.exp(-mu_b * t_b)
    phi_f = eta_f * i_matrix.output.t_f / 2 + eta_b * i_matrix.output.t_b / 2
    i_noscat_f = i_feature.i_noscat.data
    i_1el_f = i_feature.i_1el.data
    i_elpl_f = i_feature.i_elpl.data
    i_inel_f = i_feature.i_inel.data
    i_out_f = i_feature.i_out.data
    i_abs_f = i_feature.i_pi.data

    phi_b = eta_b * i_matrix.output.t / 2
    i_noscat_b = i_matrix.i_noscat.data
    i_1el_b = i_matrix.i_1el.data
    i_elpl_b = i_matrix.i_elpl.data
    i_inel_b = i_matrix.i_inel.data
    i_out_b = i_matrix.i_out.data
    i_abs_b = i_matrix.i_pi.data

    print(mu_f, mu_b, mu_bb, eta_f, eta_b, eta_bb)

    if contrast_mode == 'zpc':
        # _, i_signal_f, i_signal_b = rudolph_theta(mu_b, t_b, mu_f, t_f, delta_b, delta_f,
        #                                           i_matrix.beam.wavelength*1e6, output_dichotomy=True)
        # print fftf, 'fftf'
        # print fbtf
        # print fbtb
        # i_f = np.exp(-mu_b * t_b) * fbtb
        # i_f = i_f * (2 * np.exp(-mu_b * t_f) * fbtf  + np.exp(-mu_f * t_f) * fftf  +
        #              2 * (-1 + (eta_f - eta_b) * t_f) * np.exp(-(mu_f + mu_b) * t_f / 2) * np.sqrt(fbtf * fftf))
        # i_b = np.exp(-mu_b * t_b) * fbtb  * np.exp(-mu_b * t_f) * fbtf
        print(fbtb, fbtf, fftf)
        i_f, i_b = unified_intensities(t_b, t_f, mu_b, mu_f, mu_bb, eta_b, eta_f, eta_bb, fbtb, fbtf, fftf, phi=np.pi / 2)
        print('if, ib: ', i_f[0], i_b[0])
        i_noise_f = i_f + i_inel_f + i_elpl_f
        i_noise_b = i_b + i_inel_b + i_elpl_b
        print('inf, inb:', i_noise_f[0], i_noise_b[0])
        print('i_inel_b, i_elpl_b:', i_inel_b[0], i_elpl_b[0])
        theta = np.abs(i_f - i_b) / np.sqrt(i_noise_f + i_noise_b)

    elif contrast_mode == 'abs':
        # i_f = i_noscat_f + i_1el_f
        # i_b = i_noscat_b + i_1el_b
        i_f, i_b = unified_intensities(t_b, t_f, mu_b, mu_f, mu_bb, eta_b, eta_f, eta_bb, fbtb, fbtf, fftf, phi=0)
        i_noise_f = i_f + i_inel_f + i_elpl_f
        i_noise_b = i_b + i_inel_b + i_elpl_b
        theta = np.abs(i_f - i_b) / np.sqrt(i_noise_f + i_noise_b)
    else:
        raise ValueError('Invalid contrast mode.')

    return theta


def unified_intensities(t_b, t_f, mu_b, mu_f, mu_bb, eta_b, eta_f, eta_bb, cbtb, cbtf, cftf, phi=np.pi/2):

    a_f = np.exp(-mu_bb / 2 * t_b) * np.exp(1j * eta_bb * t_b) * np.exp(-mu_f / 2 * t_f) * np.exp(1j * eta_f * t_f) * np.sqrt(cbtb * cftf)
    a_b = np.exp(-mu_bb / 2 * t_b) * np.exp(1j * eta_bb * t_b) * np.exp(-mu_b / 2 * t_f) * np.exp(1j * eta_b * t_f) * np.sqrt(cbtb * cbtf)

    t_p = np.exp(1j * phi)

    a_b_ = a_b * t_p
    a_d = a_f - a_b
    a_f_ = a_b_ + a_d

    i_f = abs(a_f_) ** 2
    i_b = abs(a_b_) ** 2

    return i_f, i_b


# -----------------------------------------
# |                                       |
# |             ice + protein             |
# -----------------------------------------
# |     ice        |protein|    ice       |
# -----------------------------------------
# |             ice + protein             |
# |                                       |
# -----------------------------------------

# protein in ice ===========================================================
energy = np.linspace(np.log10(0.1), np.log10(30), 500)
# energy = np.linspace(np.log10(0.1), np.log10(30), 50)
energy = 10 ** energy
t_b = np.linspace(np.log10(1), np.log10(1e4), 500)
t_b = 10 ** t_b
t_f = 0.02

feature_den = 1.35
mass_frac_water = 0.708
mass_frac_protein = 1 - mass_frac_water
molar_mass_water = 18.
molar_mass_protein = 729.6
stoic_water = 1.
stoic_protein = molar_mass_water / mass_frac_water * mass_frac_protein / molar_mass_protein
formula_mixture =  'H{:.5f}C{:.5f}N{:.5f}O{:.5f}S{:.5f}'.format(2 * stoic_water + 48.6 * stoic_protein,
                                            0 * stoic_water + 32.9 * stoic_protein,
                                            0 * stoic_water + 8.9 * stoic_protein,
                                            1 * stoic_water + 8.9 * stoic_protein,
                                            0 * stoic_water + 0.6 * stoic_protein,)
density_mixture = 0.92 * 0.708 + feature_den * 0.292
print(formula_mixture)

sample_protein = Composite(matrix_compound=formula_mixture, matrix_density=density_mixture, feature_compound='H48.6C32.9N8.9O8.9S0.6',
                     feature_density=feature_den, matrix_eloss=39.3, feature_eloss=37.5, matrix_thickness=50,
                     feature_thickness=t_f, variable='t_b', feature_alias='protein', matrix_alias='tissue')
sample_ice = Composite(matrix_compound=formula_mixture, matrix_density=density_mixture, feature_compound='H2O',
                     feature_density=0.92, matrix_eloss=39.3, feature_eloss=37.5, matrix_thickness=50,
                     feature_thickness=t_f, variable='t_b', feature_alias='ice', matrix_alias='tissue')

# ==========================================================================

# Complete model:
# -----------------------------------------
# |                                       |
# |             ice + protein             |
# -----------------------------------------
# |     ice        |protein|    ice       |
# -----------------------------------------
# |             ice + protein             |
# |                                       |
# -----------------------------------------
# Protein part:
# -----------------------------------------
# |                                       |
# |             ice + protein             |
# -----------------------------------------
# | ice + protein |protein| ice + protein |
# -----------------------------------------
# |             ice + protein             |
# |                                       |
# -----------------------------------------
# Ice part:
# -----------------------------------------
# |                                       |
# |             ice + protein             |
# -----------------------------------------
# | ice + protein |  ice  | ice + protein |
# -----------------------------------------
# |             ice + protein             |
# |                                       |
# -----------------------------------------

# sample_protein = Composite(matrix_compound='H2O', matrix_density=0.92, feature_compound='H48.6C32.9N8.9O8.9S0.6',
#                      feature_density=feature_den, matrix_eloss=39.3, feature_eloss=37.5, matrix_thickness=50,
#                      feature_thickness=t_f, variable='t_b')
# sample_ice = Composite(matrix_compound='H2O', matrix_density=0.92, feature_compound='H2O',
#                      feature_density=0.92, matrix_eloss=39.3, feature_eloss=37.5, matrix_thickness=50,
#                      feature_thickness=t_f, variable='t_b')

measure = Measurement(pixel_size=t_f, n_ccd=1024, working_distance=1e4)

simulator_thickness_protein = CompositeSimulator(sample_protein)
output_protein = Output(sample_protein, step=0.001)
output_protein.t_b = t_b
output_protein.t = t_b + t_f
simulator_protein = CompositeSimulator(sample_protein)

simulator_thickness_ice = CompositeSimulator(sample_ice)
output_ice = Output(sample_ice, step=0.001)
output_ice.t_b = t_b
output_ice.t = t_b + t_f
simulator_ice = CompositeSimulator(sample_ice)

arr_zpc = np.zeros([t_b.size, energy.size])
arr_abs = np.zeros([t_b.size, energy.size])
n_zpc = np.zeros([t_b.size, energy.size])
n_abs = np.zeros([t_b.size, energy.size])

for i, kev in enumerate(energy):
    print('energy = {}'.format(kev))
    x_beam = XrayBeam(kev)

    i_matrix, i_protein, i_ftf_protein, i_btb_protein, i_btf_protein = simulator_protein.get_xray_categories(x_beam, measure, output_protein, return_aux=True)
    _,        i_ice,     i_ftf_ice,     i_btb_ice,     i_btf_ice     = simulator_ice.get_xray_categories(x_beam, measure, output_ice, return_aux=True)

    theta_zpc_complete = get_xray_theta_complete(sample_protein, sample_ice, i_ice,    i_protein, i_ftf_protein, i_btb_protein, i_ftf_ice, contrast_mode='zpc')
    # theta_zpc_complete = get_xray_theta_complete(sample_protein, sample_ice,   i_ice,    i_protein, i_ftf_protein, i_btb_protein, i_ftf_ice, contrast_mode='zpc') # replace args 1 by 1
    # theta_zpc_complete = get_xray_theta_complete(sample_protein, sample_ice, i_matrix, i_protein, i_ftf_protein, i_btb_protein, i_btf_protein, contrast_mode='zpc') # all args from simulator_protein
    # theta_zpc_complete = simulator.get_xray_theta_complete(i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='zpc')
    theta_abs_complete = get_xray_theta_complete(sample_protein, sample_ice, i_ice, i_protein, i_ftf_protein, i_btb_protein, i_ftf_ice, contrast_mode='abs')
    # theta_abs_complete = simulator.get_xray_theta_complete(i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='abs')
    nprobe_x_zpc = 25. / theta_zpc_complete ** 2
    nprobe_x_abs = 25. / theta_abs_complete ** 2
    # k_pi_f = i_matrix.constants['k_pi_b']
    k_pi_f = i_protein.constants['k_pi_f']
    k_pi_b = i_matrix.constants['k_pi_b']
    # mfp = 1. / k_pi_f
    # pen_depth = mfp
    # pen_depth = np.min([t_f, mfp])
    # print(t_f, mfp)
    # feat_mass = t_f ** 2. * density_mixture * pen_depth * 1.e-15
    # feat_mass = t_f ** 2. * feature_den * pen_depth * 1.e-15
    # # i_abs_f = i_protein.i_pi.data
    # i_abs_f = i_ftf_protein.i_pi.data
    # dose_x_zpc = nprobe_x_zpc * kev * ECharge * 1.e3 * i_abs_f / feat_mass
    # dose_x_abs = nprobe_x_abs * kev * ECharge * 1.e3 * i_abs_f / feat_mass

    dose_x_zpc = nprobe_x_zpc * kev * ECharge * 1.e3 * k_pi_f * np.exp(-k_pi_b * t_b / 2) / (feature_den * t_f ** 2 * 1e-15)
    dose_x_abs = nprobe_x_abs * kev * ECharge * 1.e3 * k_pi_f * np.exp(-k_pi_b * t_b / 2) / (feature_den * t_f ** 2 * 1e-15)

    arr_zpc[:, i] = dose_x_zpc
    arr_abs[:, i] = dose_x_abs
    n_zpc[:, i] = nprobe_x_zpc
    n_abs[:, i] = nprobe_x_abs

f = h5py.File('DEt_ter_{}_w_{}_in_{}_r{}um.h5'.format(sample_protein.feature_alias, sample_ice.feature_alias, sample_ice.matrix_alias, t_f), 'w')
f.create_dataset('energy_kev', dtype='float32', data=energy)
f.create_dataset('t_b_um', dtype='float32', data=t_b)
f.create_dataset('dose_zpc', dtype='float32', data=arr_zpc)
f.create_dataset('dose_abs', dtype='float32', data=arr_abs)
f.create_dataset('n_zpc', dtype='float32', data=n_zpc)
f.create_dataset('n_abs', dtype='float32', data=n_abs)

print(n_zpc[0, :])
print(arr_zpc[0, :])

f.close()