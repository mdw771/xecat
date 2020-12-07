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


class CompositeSimulator():
    def __init__(self, sample):
        """
        Initialize object.
        :param pixel: Pixel size in um
        :param wd: working distance in um
        :param N: number of pixels along the side of CCD array
        :param matrix_compound: formula of matrix material
        :param matrix_density: density of matrix material in g/cm3
        """

        self.results = []
        self.doses_x_thickness = []
        self.doses_x_resolution = []
        self.doses_e = []
        self.sample = sample

    def get_xray_categories(self, xray_beam, measurement, output, return_aux=False):

        energy = xray_beam.energy

        # Abbe criterion applied
        wavelen = xray_beam.wavelength
        num_aperture = measurement.get_numerical_aperture(xray_beam)
        (eta_el_f, eta_inel_f) = scat_efficiencies(energy, self.sample.feature_mw, self.sample.feature_elements,
                                                   self.sample.feature_stoic)
        (eta_el_b, eta_inel_b) = scat_efficiencies(energy, self.sample.matrix_mw, self.sample.matrix_elements,
                                                   self.sample.matrix_stoic)

        # frac1 is the fraction in forward scattered (detectable) photons that finally enter the aperture
        frac1 = inna_frac(wavelen, num_aperture)

        # frac2 is the fraction of forward scattered inelastic scattered or plural scattered photons that enter the
        # aperture and contribute to PCI background
        frac2 = num_aperture ** 2 / 2

        cs_inel_b = xraylib.CS_Compt_CP(self.sample.matrix, energy)
        cs_el_b = xraylib.CS_Rayl_CP(self.sample.matrix, energy)
        cs_pi_b = xraylib.CS_Photo_CP(self.sample.matrix, energy)
        cs_tot_b = xraylib.CS_Total_CP(self.sample.matrix, energy)

        cs_inel_f = xraylib.CS_Compt_CP(self.sample.feature, energy)
        cs_el_f = xraylib.CS_Rayl_CP(self.sample.feature, energy)
        cs_pi_f = xraylib.CS_Photo_CP(self.sample.feature, energy)
        cs_tot_f = xraylib.CS_Total_CP(self.sample.feature, energy)

        # probability per length in um-1
        k_el_b = cs_el_b * self.sample.matrix_den * 1e-4
        k_inel_b = cs_inel_b * self.sample.matrix_den * 1e-4
        k_elin_b = cs_el_b * (1 - eta_el_b) * self.sample.matrix_den * 1e-4
        k_inelin_b = cs_inel_b * (1 - eta_inel_b) * self.sample.matrix_den * 1e-4
        k_out_b = (cs_el_b * eta_el_b * self.sample.matrix_den + cs_inel_b * eta_inel_b * self.sample.matrix_den) * 1e-4
        k_pi_b = cs_pi_b * self.sample.matrix_den * 1e-4
        k_tot_b = k_inel_b + k_el_b + k_pi_b

        # probability per length in um-1
        k_el_f = cs_el_f * self.sample.feature_den * 1e-4
        k_inel_f = cs_inel_f * self.sample.feature_den * 1e-4
        k_elin_f = cs_el_f * (1 - eta_el_f) * self.sample.feature_den * 1e-4
        k_inelin_f = cs_inel_f * (1 - eta_inel_f) * self.sample.feature_den * 1e-4
        k_out_f = (
                  cs_el_f * eta_el_f * self.sample.feature_den + cs_inel_f * eta_inel_f * self.sample.feature_den) * 1e-4
        k_pi_f = cs_pi_f * self.sample.feature_den * 1e-4
        k_tot_f = k_inel_f + k_el_f + k_pi_f

        t_b = output.t_b
        t_f = output.t_f
        t = output.t

        # intensity fractions for matrix relative to primary beam
        i_noscat_b_func = lambda t: np.exp(-k_tot_b * t)
        i_1el_b_func = lambda t: k_elin_b * t * i_noscat_b_func(t)
        i_1elpc_b_func = lambda t: frac1 * i_1el_b_func(t)
        i_pc_b_func = lambda t: i_noscat_b_func(t) + i_1elpc_b_func(t)
        i_innoinel_b_func = lambda t: np.exp(-(k_inelin_b + k_out_b + k_pi_b) * t)
        i_elpl_b_func = lambda t: i_innoinel_b_func(t) - i_noscat_b_func(t) - i_1el_b_func(t)
        i_elplpc_b_func = lambda t: i_elpl_b_func(t) * frac2
        i_out_b_func = lambda t: (k_out_b / (k_out_b + k_pi_b)) * (1 - np.exp(-(k_out_b + k_pi_b) * t))
        i_pi_b_func = lambda t: (k_pi_b / (k_out_b + k_pi_b)) * (1 - np.exp(-(k_out_b + k_pi_b) * t))
        i_inel_b_func = lambda t: np.exp(-(k_out_b + k_pi_b) * t_b) - i_innoinel_b_func(t)
        i_inelpc_b_func = lambda t: i_inel_b_func(t) * frac2

        i_noscat_b = np.exp(-k_tot_b * t)
        i_1el_b = k_elin_b * t * i_noscat_b
        i_1elpc_b = frac1 * i_1el_b
        i_pc_b = i_noscat_b + i_1elpc_b
        i_innoinel_b = np.exp(-(k_inelin_b + k_out_b + k_pi_b) * t)
        i_elpl_b = i_innoinel_b - i_noscat_b - i_1el_b
        i_elplpc_b = i_elpl_b * frac2
        i_out_b = (k_out_b / (k_out_b + k_pi_b)) * (1 - np.exp(-(k_out_b + k_pi_b) * t))
        i_pi_b = (k_pi_b / (k_out_b + k_pi_b)) * (1 - np.exp(-(k_out_b + k_pi_b) * t))
        i_inel_b = np.exp(-(k_out_b + k_pi_b) * t_b) - i_innoinel_b
        i_inelpc_b = i_inel_b * frac2

        i_noscat_ftf = np.exp(-k_tot_f * t_f)
        i_1el_ftf = k_elin_f * t_f * i_noscat_ftf
        i_1elpc_ftf = frac1 * i_1el_ftf
        i_pc_ftf = i_noscat_ftf + i_1elpc_ftf
        i_innoinel_ftf = np.exp(-(k_inelin_f + k_out_f + k_pi_f) * t_f)
        i_elpl_ftf = i_innoinel_ftf - i_noscat_ftf - i_1el_ftf
        i_elplpc_ftf = i_elpl_ftf * frac2
        i_out_ftf = (k_out_f / (k_out_f + k_pi_f)) * (1 - np.exp(-(k_out_f + k_pi_f) * t_f))
        i_pi_ftf = (k_pi_f / (k_out_f + k_pi_f)) * (1 - np.exp(-(k_out_f + k_pi_f) * t_f))
        i_inel_ftf = np.exp(-(k_out_f + k_pi_f) * t_f) - i_innoinel_ftf
        i_inelpc_ftf = i_inel_ftf * frac2

        i_noscat_btf = np.exp(-k_tot_b * t_f)
        i_1el_btf = k_elin_b * t_f * i_noscat_btf
        i_innoinel_btf = np.exp(-(k_inelin_b + k_out_b + k_pi_b) * t_f)
        i_elpl_btf = i_innoinel_btf - i_noscat_btf - i_1el_btf
        i_out_btf = (k_out_b / (k_out_b + k_pi_b)) * (1 - np.exp(-(k_out_b + k_pi_b) * t_f))
        i_pi_btf = (k_pi_b / (k_out_b + k_pi_b)) * (1 - np.exp(-(k_out_b + k_pi_b) * t_f))
        i_inel_btf = np.exp(-(k_out_b + k_pi_b) * t_f) - i_innoinel_btf

        i_noscat_btb = np.exp(-k_tot_b * t_b)
        i_1el_btb = k_elin_b * t_b * i_noscat_btb
        i_innoinel_btb = np.exp(-(k_inelin_b + k_out_b + k_pi_b) * t_b)
        i_elpl_btb = i_innoinel_btb - i_noscat_btb - i_1el_btb
        i_out_btb = (k_out_b / (k_out_b + k_pi_b)) * (1 - np.exp(-(k_out_b + k_pi_b) * t_b))
        i_pi_btb = (k_pi_b / (k_out_b + k_pi_b)) * (1 - np.exp(-(k_out_b + k_pi_b) * t_b))
        i_inel_btb = np.exp(-(k_out_b + k_pi_b) * t_b) - i_innoinel_btb

        i_noscat_f = np.exp(-k_tot_b * t_b) * np.exp(-k_tot_f * t_f)
        i_1el_f = (k_elin_b * t_b + k_elin_f * t_f) * i_noscat_f
        i_innoinel_f = np.exp(-(k_out_b + k_inelin_b + k_pi_b) * t_b) * np.exp(-(k_out_f + k_inelin_f + k_pi_f) * t_f)
        i_elpl_f = i_innoinel_f - i_noscat_f - i_1el_f

        i_out_ff = (1 - i_out_b_func(t_b / 2) - i_pi_b_func(t_b / 2)) * (k_out_f / (k_out_f + k_pi_f)) * (1 - np.exp(-(k_out_f + k_pi_f) * t_f))
        i_pi_ff = (1 - i_out_b_func(t_b / 2) - i_pi_b_func(t_b / 2)) * (k_pi_f / (k_out_f + k_pi_f)) * (1 - np.exp(-(k_out_f + k_pi_f) * t_f))


        i_out_f = i_out_b_func(t_b / 2) + i_out_ff + (1 - i_out_b_func(t_b / 2) - i_pi_b_func(t_b / 2) - i_out_ff - i_pi_ff) * i_out_b_func(t_b / 2)
        # i_out_f = (1 - i_out_b_func(t_b / 2) - i_pi_b_func(t_b / 2)) * (k_out_f / (k_out_f + k_pi_f)) * (1 - np.exp(-(k_out_f + k_pi_f) * t_f)) + i_out_b_func(t_b)
        i_pi_f = i_pi_b_func(t_b / 2) + i_pi_ff + (1 - i_out_b_func(t_b / 2) - i_pi_b_func(t_b / 2) - i_out_ff - i_pi_ff) * i_pi_b_func(t_b / 2)
        # i_pi_f = (1 - i_out_b_func(t_b / 2) - i_pi_b_func(t_b / 2)) * (k_out_f / (k_out_f + k_pi_f)) * (1 - np.exp(-(k_out_f + k_pi_f) * t_f)) + i_pi_b_func(t_b)

        # i_inel_f = np.exp(-(k_out_f + k_pi_f) * t_f) - i_innoinel_f
        i_inel_f = 1 - i_out_f - i_pi_f - i_innoinel_f

        # print('inside i_inel', i_inel_f[0], i_inel_b[0])
        # print('inside i_innoinel', i_innoinel_f[0], i_innoinel_b[0])
        # print('inside i_out', i_out_f[0], i_out_b[0])
        # print('inside i_pi', i_pi_f[0], i_pi_b[0])
        print('i_out_ff', i_out_ff)
        print('i_pi_ff', i_pi_ff)
        print('i_out_f', i_out_f)
        print('i_pi_f', i_pi_f)
        print('i_out_b', i_out_f)
        print('i_pi_b', i_pi_b)

        x_matrix = result_holder(xray_beam, output, measurement,
                                 constants={'k_tot_b': k_tot_b, 'k_pi_f': k_pi_f, 'k_pi_b': k_pi_b, 'k_el_b': k_el_b},
                                 i_noscat=i_noscat_b, i_1el=i_1el_b,
                                 i_innoinel=i_innoinel_b, i_elpl=i_elpl_b, i_out=i_out_b,
                                 i_pi=i_pi_b, i_inel=i_inel_b)
        x_feature = result_holder(xray_beam, output, measurement,
                                  constants={'k_tot_f': k_tot_f, 'k_pi_f': k_pi_f, 'k_pi_b': k_pi_b, 'k_el_f': k_el_f},
                                  i_noscat=i_noscat_f, i_1el=i_1el_f,
                                  i_innoinel=i_innoinel_f, i_elpl=i_elpl_f, i_out=i_out_f,
                                  i_pi=i_pi_f, i_inel=i_inel_f)
        x_ftf = result_holder(xray_beam, output, measurement,
                                 constants={'k_tot_b': k_tot_b, 'k_pi_f': k_pi_f, 'k_pi_b': k_pi_b, 'k_el_b': k_el_b},
                                 i_noscat=i_noscat_ftf, i_1el=i_1el_ftf, i_1elpc=i_1elpc_ftf, i_pc=i_pc_ftf,
                                 i_innoinel=i_innoinel_ftf, i_elpl=i_elpl_ftf, i_elplpc=i_elplpc_ftf, i_out=i_out_ftf,
                                 i_pi=i_pi_ftf, i_inel=i_inel_ftf, i_inelpc=i_inelpc_ftf, i_pi_ff=i_pi_ff)
        x_btf = result_holder(xray_beam, output, measurement,
                                 constants={'k_tot_b': k_tot_b, 'k_pi_f': k_pi_f, 'k_pi_b': k_pi_b, 'k_el_b': k_el_b},
                                 i_noscat=i_noscat_btf, i_1el=i_1el_btf,
                                 i_innoinel=i_innoinel_btf, i_elpl=i_elpl_btf, i_out=i_out_btf,
                                 i_pi=i_pi_btf, i_inel=i_inel_btf)
        x_btb = result_holder(xray_beam, output, measurement,
                                 constants={'k_tot_b': k_tot_b, 'k_pi_f': k_pi_f, 'k_pi_b': k_pi_b, 'k_el_b': k_el_b},
                                 i_noscat=i_noscat_btb, i_1el=i_1el_btb,
                                 i_innoinel=i_innoinel_btb, i_elpl=i_elpl_btb, i_out=i_out_btb,
                                 i_pi=i_pi_btb, i_inel=i_inel_btb)

        if return_aux:
            return x_matrix, x_feature, x_ftf, x_btb, x_btf
        else:
            return x_matrix, x_feature

    def get_e_categories(self, e_beam, measurement, output):

        eta = 1 - 4.12 / 10

        # compute cross sections in cm2 (see Langmore 1992)
        cs_el_b = cse_elastic(e_beam.energy, self.sample.matrix_elements, self.sample.matrix_stoic)
        cs_inel_b = cse_inelastic(self.sample.matrix_eloss, e_beam.energy, self.sample.matrix_elements,
                                  self.sample.matrix_stoic)
        cs_el_f = cse_elastic(e_beam.energy, self.sample.feature_elements, self.sample.feature_stoic)
        cs_inel_f = cse_inelastic(self.sample.feature_eloss, e_beam.energy, self.sample.feature_elements,
                                  self.sample.feature_stoic)

        # probability per thickness in um-1
        delta_b = self.sample.matrix_den / self.sample.matrix_mw * Nav
        delta_f = self.sample.feature_den / self.sample.feature_mw * Nav
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
        self.k_inel_f = k_inel_f
        self.k_inel_b = k_inel_b
        self.k_el_b = k_el_b
        self.k_el_f = k_el_f

        t_f = output.t_f
        t_b = output.t_b
        t = output.t

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

        e_matrix = result_holder(e_beam, output, measurement,
                                 constants={'k_inel_b': k_inel_b, 'k_tot_b': k_tot_b, 'k_elin_b': k_elin_b,
                                            'k_out_b': k_out_b, 'k_el_b': k_el_b},
                                 i_noscat=i_noscat_b, i_1el=i_1el_b, i_innoinel=i_innoinel_b, i_elpl=i_elpl_b,
                                 i_in=i_in_b, i_inel=i_inel_b)
        e_feature = result_holder(e_beam, output, measurement,
                                  constants={'k_inel_f': k_inel_f, 'k_tot_f': k_tot_f, 'k_elin_f': k_elin_f,
                                             'k_out_f': k_out_f, 'k_el_f': k_el_f},
                                  i_noscat=i_noscat_f, i_1el=i_1el_f, i_innoinel=i_innoinel_f, i_elpl=i_elpl_f,
                                  i_in=i_in_f, i_inel=i_inel_f, i_1elf_f=i_1elf_f)

        return e_matrix, e_feature

    def get_xray_dose_thickness(self, i_matrix, i_feature, i_ftf, i_btb, i_btf, hold_thickness=10):

        assert isinstance(i_matrix, result_holder)
        assert isinstance(i_feature, result_holder)

        theta_x_zpc = self.get_xray_theta_complete(i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='zpc')
        theta_x_abs = self.get_xray_theta_complete(i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='abs')

        # calculate contrast parameter
        ridelta_f = 1 - xraylib.Refractive_Index_Re(self.sample.feature, i_feature.beam.energy, self.sample.feature_den)
        ridelta_b = 1 - xraylib.Refractive_Index_Re(self.sample.matrix, i_matrix.beam.energy, self.sample.matrix_den)
        k_pi_f = i_feature.constants['k_pi_f']
        k_pi_b = i_matrix.constants['k_pi_b']
        i_abs_f = i_ftf.i_pi_ff.data

        # calculate dose wrt thickness
        t_f = i_matrix.output.t_f
        t_b = i_matrix.output.t_b
        wavelen = i_feature.beam.wavelength
        pixel = i_feature.measurement.pixel
        k_tot_b_x = i_matrix.constants['k_tot_b']
        feature_den = self.sample.feature_den
        energy = i_feature.beam.energy

        nprobe_x_zpc = 25. / theta_x_zpc ** 2
        nprobe_x_abs = 25. / theta_x_abs ** 2
        # feat_mass = pixel ** 2. * feature_den * t_f * 1.e-15
        # dose_x_zpc = nprobe_x_zpc * energy * ECharge * 1.e3 * i_abs_f / feat_mass
        # dose_x_abs = nprobe_x_abs * energy * ECharge * 1.e3 * i_abs_f / feat_mass

        dose_x_zpc = nprobe_x_zpc * energy * ECharge * 1.e3 * k_pi_f * np.exp(-k_pi_b * t_b / 2) / (feature_den * t_f ** 2 * 1e-15)
        dose_x_abs = nprobe_x_abs * energy * ECharge * 1.e3 * k_pi_f * np.exp(-k_pi_b * t_b / 2) / (feature_den * t_f ** 2 * 1e-15)

        print('===== Complete model ======')
        print('theta_x_zpc', theta_x_zpc)
        print('nprobe_x_zpc', nprobe_x_zpc)
        print('dose_x_zpc', dose_x_zpc)
        print('theta_x_abs', theta_x_abs)
        print('nprobe_x_abs', nprobe_x_abs)
        print('dose_x_abs', dose_x_abs)
        print('===========================')

        self.doses_x_thickness.append((dose_x_zpc, dose_x_abs, i_matrix, i_feature))

    def get_xray_dose_resolution(self, i_matrix, i_feature, i_ftf, i_btb, i_btf):

        assert isinstance(i_matrix, result_holder)
        assert isinstance(i_feature, result_holder)

        theta_x_zpc = self.get_xray_theta_complete(i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='zpc')
        theta_x_abs = self.get_xray_theta_complete(i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='abs')

        # calculate dose wrt resolution
        t_f = i_feature.output.t_f
        delta_ls = t_f
        t_b = i_matrix.output.t_b
        wavelen = i_feature.beam.wavelength
        pixel = i_feature.measurement.pixel
        feature_den = self.sample.feature_den
        energy = i_feature.beam.energy

        # calculate contrast parameter
        k_pi_f = i_feature.constants['k_pi_f']
        k_pi_b = i_matrix.constants['k_pi_b']
        i_abs_f = i_ftf.i_pi_ff.data

        nprobe_x_zpc = 25. / theta_x_zpc ** 2
        nprobe_x_abs = 25. / theta_x_abs ** 2

        # take mean free path or feature thickness, whichever is smaller
        # mu_arr = [1. / k_pi_f] * t_f.size
        # pen = np.squeeze(np.dstack([mu_arr, t_f]))
        # print(mu_arr)
        # pen_depth = np.min(pen, axis=1)
        # feat_mass = delta_ls ** 2. * feature_den * pen_depth * 1.e-15
        # dose_x_zpc = nprobe_x_zpc * energy * ECharge * 1.e3 * i_abs_f / feat_mass
        # dose_x_abs = nprobe_x_abs * energy * ECharge * 1.e3 * i_abs_f / feat_mass

        dose_x_zpc = nprobe_x_zpc * energy * ECharge * 1.e3 * k_pi_f * np.exp(-k_pi_b * t_b / 2) / (feature_den * delta_ls ** 2 * 1e-15)
        dose_x_abs = nprobe_x_abs * energy * ECharge * 1.e3 * k_pi_f * np.exp(-k_pi_b * t_b / 2) / (feature_den * delta_ls ** 2 * 1e-15)

        self.doses_x_resolution.append((dose_x_zpc, dose_x_abs, i_matrix, i_feature))

    def get_e_dose(self, i_matrix, i_feature, hold_thickness=10):

        assert isinstance(i_matrix, result_holder) and isinstance(i_feature, result_holder)

        # calculate dose wrt thickness

        i_innoinel_f = i_feature.i_innoinel.data
        i_noscat_f = i_feature.i_noscat.data
        i_1elf_f = i_feature.i_1el.data
        i_1el_f = i_feature.i_1el.data
        i_in_f = i_feature.i_in.data
        i_innoinel_b = i_matrix.i_innoinel.data
        i_in_b = i_matrix.i_in.data
        i_noscat_b = i_matrix.i_noscat.data
        i_1el_b = i_matrix.i_1el.data
        # calculate contrast parameters
        theta_e = np.abs(i_noscat_f - i_noscat_b + i_1el_f - i_1el_b) + 2 * np.sqrt(i_noscat_f * i_1elf_f)
        theta_e = theta_e / np.sqrt(i_in_f + i_in_b)
        theta_ef = np.abs(i_noscat_f - i_noscat_b + i_1el_f - i_1el_b) + 2 * np.sqrt(i_noscat_f * i_1elf_f)
        theta_ef = theta_ef / np.sqrt(i_innoinel_f + i_innoinel_b)
        eloss_f = self.sample.feature_eloss
        eloss_b = self.sample.matrix_eloss
        pixel = i_matrix.measurement.pixel
        feature_den = self.sample.feature_den
        nprobe_e = 25. / theta_e ** 2
        nprobe_ef = 25. / theta_ef ** 2
        dose_e = nprobe_e * eloss_f * self.k_inel_f / (pixel ** 2 * 1e-12 * feature_den) * ECharge * 1e3
        dose_ef = nprobe_ef * eloss_f * self.k_inel_f / (pixel ** 2 * 1e-12 * feature_den) * ECharge * 1e3

        pixel_area = i_matrix.measurement.pixel ** 2
        fluence_e = nprobe_e / pixel_area * 1e-6
        fluence_ef = nprobe_ef / pixel_area * 1e-6
        print('Electron fluence: {} (no filter); {} (filtered); in nm^-2'.format(fluence_e, fluence_ef))

        print('n_e_filtered', nprobe_ef)
        print('n_e_unfiltered', nprobe_e)
        print('dose_e_filtered', dose_ef)
        print('dopse_e_unfiltered', dose_e)

        # calculate dose wrt resolution

        t_f = i_matrix.output.delta
        t_b = hold_thickness
        k_tot_b = i_matrix.constants['k_tot_b']
        k_elin_b = i_matrix.constants['k_elin_b']
        k_out_b = i_matrix.constants['k_out_b']
        k_inel_b = i_matrix.constants['k_inel_b']
        k_tot_f = i_feature.constants['k_tot_f']
        k_inel_f = i_feature.constants['k_inel_f']
        k_elin_f = i_feature.constants['k_elin_f']
        k_out_f = i_feature.constants['k_out_f']

        i_noscat_f = np.exp(-k_tot_b * t_b) * np.exp(-k_tot_f * t_f)
        i_1el_f = (k_elin_b * t_b + k_elin_f * t_f) * i_noscat_f
        i_1elf_f = k_elin_f * t_f * i_noscat_f
        i_innoinel_f = np.exp(-(k_out_b + k_inel_b) * t_b) * np.exp(-(k_out_f + k_inel_f) * t_f)
        i_in_f = np.exp(-k_out_b * t_b) * np.exp(-k_out_f * t_f)
        i_innoinel_b = np.exp(-(k_inel_b + k_out_b) * (t_b + t_f))
        i_in_b = np.exp(-k_out_b * (t_b + t_f))

        theta_e = np.abs(i_innoinel_f - i_innoinel_b) + 2 * np.sqrt(i_noscat_f * i_1elf_f)
        theta_e = theta_e / np.sqrt(i_in_f + i_in_b)
        theta_ef = np.abs(i_innoinel_f - i_innoinel_b) + 2 * np.sqrt(i_noscat_f * i_1elf_f)
        theta_ef = theta_ef / np.sqrt(i_innoinel_f + i_innoinel_b)

        delta_ls = i_matrix.output.delta
        feature_den = self.sample.feature_den
        nprobe_e = 25. / theta_e ** 2
        nprobe_ef = 25. / theta_ef ** 2

        # nprobe_ef = 1.e6
        # dose_ef_res = nprobe_ef * eloss_f * self.k_inel_f / (feature_den * 1.e-12) * ECharge * 1e3
        # print dose_ef_res

        dose_e_res = nprobe_e * eloss_f * self.k_inel_f / (delta_ls ** 2 * 1e-12 * feature_den) * ECharge * 1e3
        dose_ef_res = nprobe_ef * eloss_f * self.k_inel_f / (delta_ls ** 2 * 1e-12 * feature_den) * ECharge * 1e3

        self.doses_e.append((dose_e, dose_ef, i_matrix, i_feature))

    def get_xray_theta_complete(self, i_matrix, i_feature, i_ftf, i_btb, i_btf, contrast_mode='zpc', numerical=False):

        (delta_f, beta_f), (delta_b, beta_b) = self.sample.get_ri_delta_beta(i_matrix.beam)
        # (delta_f, _), (delta_b, _) = self.sample.get_ri_delta_beta_henke(i_matrix.beam)
        # mu_f, mu_b = 4 * np.pi * np.array([beta_f, beta_b]) / (i_matrix.beam.wavelength * 1e6)
        mu_b = i_matrix.constants['k_pi_b']
        mu_f = i_feature.constants['k_pi_f']
        eta_f, eta_b = 2 * np.pi * np.array([delta_f, delta_b]) / (i_matrix.beam.wavelength * 1e6)
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
            print('delta_f, delta_b:', delta_f, delta_b)
            print('beta_f, beta_b:', beta_f, beta_b)
            print('C_b(tb), C_b(tf), C_f(tf), ', fbtb, fbtf, fftf)
            i_f, i_b = unified_intensities(t_b, t_f, mu_b, mu_f, eta_b, eta_f, fbtb, fbtf, fftf, phi=np.pi/2, numerical=numerical)
            i_noise_f = i_f + i_inel_f + i_elpl_f
            i_noise_b = i_b + i_inel_b + i_elpl_b
            print('I_f, I_b: ', i_f[0], i_b[0])
            print('I_noise_f, I_noise_b:', i_noise_f[0], i_noise_b[0])
            print('I_inel_b, I_elpl_b:', i_inel_b[0], i_elpl_b[0])
            theta = np.abs(i_f - i_b) / np.sqrt(i_noise_f + i_noise_b)
            print(theta)


        elif contrast_mode == 'abs':
            # i_f = i_noscat_f + i_1el_f
            # i_b = i_noscat_b + i_1el_b
            i_f, i_b = unified_intensities(t_b, t_f, mu_b, mu_f, eta_b, eta_f, fbtb, fbtf, fftf, phi=0)
            i_noise_f = i_f + i_inel_f + i_elpl_f
            i_noise_b = i_b + i_inel_b + i_elpl_b
            theta = np.abs(i_f - i_b) / np.sqrt(i_noise_f + i_noise_b)
        else:
            raise ValueError('Invalid contrast mode.')

        return theta

    def get_xray_theta_simple(self, i_matrix, i_feature, contrast_mode='zpc'):

        (delta_f, _), (delta_b, _) = self.sample.get_ri_delta_beta(i_matrix.beam)


        # mu_f, mu_b = 4 * np.pi * np.array([beta_f, beta_b]) / (i_matrix.beam.wavelength * 1e6)
        mu_b = i_matrix.constants['k_pi_b']
        mu_f = i_feature.constants['k_pi_f']
        eta_f, eta_b = 2 * np.pi * np.array([delta_f, delta_b]) / (i_matrix.beam.wavelength * 1e6)

        t_b = i_matrix.output.t_b
        t_f = i_matrix.output.t_f

        # i_noscat_f = i_feature.i_noscat.data
        # i_1el_f = i_feature.i_1el.data
        # i_elpl_f = i_feature.i_elpl.data
        # i_inel_f = i_feature.i_inel.data
        # i_out_f = i_feature.i_out.data
        # i_abs_f = i_feature.i_pi.data
        #
        # i_noscat_b = i_matrix.i_noscat.data
        # i_1el_b = i_matrix.i_1el.data
        # i_elpl_b = i_matrix.i_elpl.data
        # i_inel_b = i_matrix.i_inel.data
        # i_out_b = i_matrix.i_out.data
        # i_abs_b = i_matrix.i_pi.data
        if contrast_mode == 'zpc':
            fbtb = fftf = fbtf = 1
            i_f = np.exp(-mu_b * t_b) * fbtb ** 2
            i_f = i_f * (2 * np.exp(-mu_b * t_f) * fbtf ** 2 + np.exp(-mu_f * t_f) * fftf ** 2 +
                         2 * (-1 + (eta_f - eta_b) * t_f) * np.exp(-(mu_f + mu_b) * t_f / 2) * fbtf * fftf)
            i_b = np.exp(-mu_b * t_b) * fbtb ** 2 * np.exp(-mu_b * t_f) * fbtf ** 2
            # i_f = np.exp(-mu_b * (t_b - t_f)) * (
            #     2 * np.exp(-mu_b * t_f) + np.exp(-mu_f * t_f) +
            #     2 * np.exp(-mu_f * t_f / 2 - mu_b * t_f / 2) * np.cos(eta_f * t_f - eta_b * t_f - np.pi / 2) -
            #     2 * np.exp(-mu_f * t_f / 2 - mu_b * t_f / 2) * np.cos(eta_f * t_f - eta_b * t_f)
            # )
            # i_b = np.exp(-mu_b * (t_b - t_f)) * np.exp(-mu_b * t_f)

            theta = np.abs(i_f - i_b) / np.sqrt(i_f + i_b)

        elif contrast_mode == 'abs':
            i_f = np.exp(-mu_f * t_f) * np.exp(-mu_b * t_b)
            i_b = np.exp(-mu_b * (t_f + t_b))
            # i_f = i_noscat_f + i_1el_f
            # i_b = i_noscat_b + i_1el_b
            theta = np.abs(i_f - i_b) / np.sqrt(i_f + i_b)
        else:
            raise ValueError('Invalid contrast mode.')

        return theta

    def get_xray_dose_simple(self, theta_x_zpc, theta_x_abs, i_feature, i_matrix):

        # calculate dose wrt resolution
        t_f = i_feature.output.t_f
        delta_ls = t_f
        t_b = i_matrix.output.t_b
        wavelen = i_feature.beam.wavelength
        pixel = i_feature.measurement.pixel
        feature_den = self.sample.feature_den
        energy = i_feature.beam.energy

        # calculate contrast parameter
        k_pi_f = i_feature.constants['k_pi_f']
        k_pi_b = i_matrix.constants['k_pi_b']

        nprobe_x_zpc = 25. / theta_x_zpc ** 2
        nprobe_x_abs = 25. / theta_x_abs ** 2

        dose_x_zpc = nprobe_x_zpc * energy * ECharge * 1.e3 * k_pi_f * np.exp(-k_pi_b * t_b / 2) / (feature_den * delta_ls ** 2 * 1e-15)
        dose_x_abs = nprobe_x_abs * energy * ECharge * 1.e3 * k_pi_f * np.exp(-k_pi_b * t_b / 2) / (feature_den * delta_ls ** 2 * 1e-15)

        return dose_x_zpc, dose_x_abs


    def get_xray_theta_thin(self, i_matrix, i_feature, contrast_mode='zpc'):

        (delta_f, _), (delta_b, _) = self.sample.get_ri_delta_beta(i_matrix.beam)
        # mu_f, mu_b = 4 * np.pi * np.array([beta_f, beta_b]) / (i_matrix.beam.wavelength * 1e6)
        mu_b = i_matrix.constants['k_pi_b']
        mu_f = i_feature.constants['k_pi_f']
        eta_f, eta_b = 2 * np.pi * np.array([delta_f, delta_b]) / (i_matrix.beam.wavelength * 1e6)
        t_b = i_matrix.output.t_b
        t_f = i_matrix.output.t_f
        if contrast_mode == 'zpc':
            theta = np.sqrt(2) * t_f * np.abs(eta_f - eta_b) * np.exp(-mu_b * t_b / 2)
        elif contrast_mode == 'abs':
            theta = t_f / np.sqrt(2) * np.abs(mu_f - mu_b) * np.exp(-mu_b * t_b / 2)
        else:
            raise ValueError('Invalid contrast mode.')

        return theta


def unified_intensities(t_b, t_f, mu_b, mu_f, eta_b, eta_f, cbtb, cbtf, cftf, phi=np.pi/2, numerical=False):

    if not numerical:
        if phi == np.pi/2:
            term = 2 * ((eta_f - eta_b) * t_f - 1)

        else:
            term = 2 * (np.cos(phi + (eta_b - eta_f) * t_f) - np.cos((eta_b - eta_f) * t_f))

        i_f = np.exp(-mu_b * t_b) * cbtb * (term *
              np.exp(-(mu_b + mu_f) * t_f / 2) * np.sqrt(cbtf * cftf) +
              (2 - 2 * cos(phi)) * np.exp(-mu_b * t_f) * cbtf + np.exp(-mu_f * t_f) * cftf)
        i_b = np.exp(-mu_b * (t_b + t_f)) * cbtb * cbtf
    else:
        a_f = np.exp(-mu_b / 2 * t_b) * np.exp(1j * eta_b * t_b) * np.exp(-mu_f / 2 * t_f) * np.exp(1j * eta_f * t_f) * np.sqrt(cbtb * cftf)
        a_b = np.exp(-mu_b / 2 * t_b) * np.exp(1j * eta_b * t_b) * np.exp(-mu_b / 2 * t_f) * np.exp(1j * eta_b * t_f) * np.sqrt(cbtb * cbtf)

        t_p = np.exp(1j * phi)

        a_b_ = a_b * t_p
        a_d = a_f - a_b
        a_f_ = a_b_ + a_d

        i_f = abs(a_f_) ** 2
        i_b = abs(a_b_) ** 2

    return i_f, i_b


if __name__ == '__main__':

    sample = Composite(matrix_compound='H2O', matrix_density=0.92, feature_compound='H48.6C32.9N8.9O8.9S0.6',
                       feature_density=1.35, matrix_eloss=39.3, feature_eloss=37.5, matrix_thickness=0,
                       feature_thickness=0.01, variable='t_b')
    measure = Measurement(pixel_size=0.01, n_ccd=1024, working_distance=1e4)
    simulator_thickness = CompositeSimulator(sample)
    output_x = Output(sample, step=0)
    output_e = Output(sample, step=0)

    # xray with thickness
    for energy in [0.5, 10]:
        x_beam = XrayBeam(energy)
        i_matrix, i_feature, i_ftf, i_btb, i_btf = simulator_thickness.get_xray_categories(x_beam, measure, output_x,
                                                                                           return_aux=True)
        # =================================
        simulator_thickness.get_xray_dose_thickness(i_matrix, i_feature, i_ftf, i_btb, i_btf)
        # =================================
        theta_x_rudolph_zpc = simulator_thickness.get_xray_theta_simple(i_matrix, i_feature, 'zpc')
        theta_x_rudolph_abs = simulator_thickness.get_xray_theta_simple(i_matrix, i_feature, 'abs')
        print('theta_x_rudolph_zpc', theta_x_rudolph_zpc)
        print('theta_x_rudolph_abs', theta_x_rudolph_abs)
        print('n_x_rudolph_zpc', 25. / theta_x_rudolph_zpc ** 2)
        print('n_x_rudolph_abs', 25. / theta_x_rudolph_abs ** 2)
        dose_x_rudolph_zpc, dose_x_rudolph_abs = simulator_thickness.get_xray_dose_simple(theta_x_rudolph_zpc, theta_x_rudolph_abs, i_feature, i_matrix)
        print('dose_x_rudolph_zpc', dose_x_rudolph_zpc)
        print('dose_x_rudolph_abs', dose_x_rudolph_abs)
        # =================================


    # e- with thickness
    for energy in [100, 300]:
        e_beam = ElectronBeam(energy)
        i_matrix, i_feature = simulator_thickness.get_e_categories(e_beam, measure, output_e)
        simulator_thickness.get_e_dose(i_matrix, i_feature)

    # xray with resolution
    sample = Composite(matrix_compound='H2O', matrix_density=0.92, feature_compound='H48.6C32.9N8.9O8.9S0.6',
                       feature_density=1.35, matrix_eloss=39.3, feature_eloss=37.5, matrix_thickness=10,
                       feature_thickness=2, variable='t_f')
    measure = Measurement(pixel_size=0.01, n_ccd=1024, working_distance=1e4)
    simulator_resolution = CompositeSimulator(sample)
    output_x = Output(sample, step=0.001)
    for energy in [0.5, 10]:
        x_beam = XrayBeam(energy)
        i_matrix, i_feature, i_ftf, i_btb, i_btf = simulator_resolution.get_xray_categories(x_beam, measure, output_x,
                                                                                 return_aux=True)
        simulator_resolution.get_xray_dose_resolution(i_matrix, i_feature, i_ftf, i_btb, i_btf)


    sys.exit()
