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
        self.doses_x = []
        self.doses_e = []
        self.sample = sample

    def get_xray_categories(self, xray_beam, measurement, output):

        cs_tot_b = xraylib.CS_Total_CP(self.sample.matrix, xray_beam.energy)
        cs_pi_f = xraylib.CS_Photo_CP(self.sample.feature, xray_beam.energy)
        cs_tot_f = xraylib.CS_Total_CP(self.sample.feature, xray_beam.energy)

        # probability per length in um-1
        k_tot_b = cs_tot_b * self.sample.matrix_den * 1e-4
        k_pi_f = cs_pi_f * self.sample.feature_den * 1e-4
        k_tot_f = cs_tot_f * self.sample.feature_den * 1e-4

        t_b = output.t_b
        t_f = self.sample.t_f

        i_abs_f = np.exp(-k_tot_b * t_b / 2.) * (1. - np.exp(-k_tot_f * t_f))
        i_xray = result_holder(xray_beam, output, measurement, constants={'k_tot_b_x': k_tot_b}, i_pi=i_abs_f)

        return i_xray

    def get_e_categories(self, e_beam, measurement, output):

        eta = 1 - 4.12 / 10

        # compute cross sections in cm2 (see Langmore 1992)
        cs_el_b = cse_elastic(e_beam.energy, self.sample.matrix_elements, self.sample.matrix_stoic)
        cs_inel_b = cse_inelastic(self.sample.matrix_eloss, e_beam.energy, self.sample.matrix_elements, self.sample.matrix_stoic)
        cs_el_f = cse_elastic(e_beam.energy, self.sample.feature_elements, self.sample.feature_stoic)
        cs_inel_f = cse_inelastic(self.sample.feature_eloss, e_beam.energy, self.sample.feature_elements, self.sample.feature_stoic)

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

        t = output.t
        t_b = output.t_b
        t_f = self.sample.t_f

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

        e_matrix = result_holder(e_beam, output, measurement, constants={'k_inel_f': k_inel_f}, i_noscat=i_noscat_b,
                                 i_1el=i_1el_b, i_innoinel=i_innoinel_b, i_elpl=i_elpl_b, i_in=i_in_b, i_inel=i_inel_b)
        e_feature = result_holder(e_beam, output, measurement, constants={'k_inel_f': k_inel_f}, i_noscat=i_noscat_f,
                                  i_1el=i_1el_f, i_innoinel=i_innoinel_f, i_elpl=i_elpl_f, i_in=i_in_f, i_inel=i_inel_f,
                                  i_1elf_f=i_1elf_f)

        return e_matrix, e_feature

    def get_xray_dose(self, i_xray):

        assert isinstance(i_xray, result_holder)

        i_abs_f = i_xray.i_pi.data

        # calculate contrast parameter
        ridelta_f = 1 - xraylib.Refractive_Index_Re(self.sample.feature, i_xray.beam.energy, self.sample.feature_den)
        ridelta_b = 1 - xraylib.Refractive_Index_Re(self.sample.matrix, i_xray.beam.energy, self.sample.matrix_den)

        t_f = self.sample.t_f
        t_b = i_xray.output.t_b
        wavelen = i_xray.beam.wavelength
        pixel = i_xray.measurement.pixel
        k_tot_b_x = i_xray.constants['k_tot_b_x']
        feature_den = self.sample.feature_den
        energy = i_xray.beam.energy
        theta_x = 2 * np.sqrt(2) * np.pi * t_f / (wavelen * 1e6) * np.abs(ridelta_f - ridelta_b) * np.exp(-k_tot_b_x * t_b / 2)
        nprobe_x = 25. / theta_x ** 2
        feat_mass = pixel ** 2. * feature_den * t_f * 1.e-15
        dose_x = nprobe_x * energy * ECharge * 1.e3 * i_abs_f / feat_mass

        self.doses_x.append((dose_x, i_xray))

    def get_e_dose(self, i_matrix, i_feature):

        assert isinstance(i_matrix, result_holder) and isinstance(i_feature, result_holder)

        i_innoinel_f = i_feature.i_innoinel.data
        i_noscat_f = i_feature.i_noscat.data
        i_1elf_f = i_feature.i_1elf_f.data
        i_in_f = i_feature.i_in.data
        i_innoinel_b = i_matrix.i_innoinel.data
        i_in_b = i_matrix.i_in.data

        # calculate contrast parameters
        theta_e = np.abs(i_innoinel_f - i_innoinel_b) + 2 * np.sqrt(i_noscat_f * i_1elf_f)
        theta_e = theta_e / np.sqrt(i_in_f + i_in_b)
        theta_ef = np.abs(i_innoinel_f - i_innoinel_b) + 2 * np.sqrt(i_noscat_f * i_1elf_f)
        theta_ef = theta_ef / np.sqrt(i_innoinel_f + i_innoinel_b)

        # calculate dose
        eloss_f = self.sample.feature_eloss
        eloss_b = self.sample.matrix_eloss
        pixel = i_matrix.measurement.pixel
        feature_den = self.sample.feature_den
        nprobe_e = 25. / theta_e ** 2
        nprobe_ef = 25. / theta_ef ** 2
        dose_e = nprobe_e * eloss_f * self.k_inel_f / (pixel ** 2 * 1e-12 * feature_den) * ECharge * 1e3
        dose_ef = nprobe_ef * eloss_f * self.k_inel_f / (pixel ** 2 * 1e-12 * feature_den) * ECharge * 1e3

        self.doses_e.append((dose_e, dose_ef, i_matrix, i_feature))

    def plot_dose(self):

        fig = plt.figure(figsize=(6, 5))
        matplotlib.rcParams['pdf.fonttype'] = 'truetype'
        fontProperties = {'family': 'serif', 'serif': ['Times New Roman'], 'weight': 'normal', 'size': 12}
        plt.rc('font', **fontProperties)
        ax = fig.add_subplot(1, 1, 1)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%g'))
        ax.xaxis.set_major_locator(plt.FixedLocator([0.05, 0.1, 0.3, 1, 3, 10, 30, 100, 200]))

        t = self.doses_x[0][1].output.t
        labpt = int(t.size / 10 * 0.02)
        for (dose_x, i_xray) in self.doses_x:
            energy = i_xray.beam.energy
            plt.plot(t, dose_x, color='red')
            plt.text(t[labpt], dose_x[labpt], str(energy) + ' keV soft X-ray', fontsize=9, color='red')
        plt.hold(True)

        t = self.doses_e[0][2].output.t
        labpt = int(t.size / 10 * 1.5)
        for (dose_e, dose_ef, i_matrix, i_feature) in self.doses_e:
            plt.plot(t, dose_e, color='gray')
            plt.text(t[labpt], dose_e[labpt], '(no energy filter)', fontsize=9, color='grey')
            plt.plot(t, dose_ef, color='black')
            plt.text(t[[labpt]], dose_ef[labpt], '(energy filter)', fontsize=9, color='black')

        plt.hold(False)
        plt.axis([0.05, 200, 1e5, 1e12])
        plt.xlabel('Thickness ($\mu$m)')
        plt.ylabel('Dose (Gray)')

        plt.show()
        fig.savefig('dose.pdf', type='pdf')


thickness = 200
step = 0.01
feature_thickness = 0.01
sample = Composite(matrix_compound='H2O', matrix_density=0.92, feature_compound='H48.6C32.9N8.9O8.9S0.6',
                   feature_density=1.35, matrix_eloss=46, feature_eloss=37.5, thickness=200, feature_thickness=0.01)
measure = Measurement(pixel_size=0.01, n_ccd=1024, working_distance=1e4)
simulator = CompositeSimulator(sample)
output_x = Output(sample, step=0.01)
output_e = Output(sample, step=0.01, overflow_limit=10)

for energy in [0.5, 10]:
    x_beam = XrayBeam(energy)
    i_xray = simulator.get_xray_categories(x_beam, measure, output_x)
    simulator.get_xray_dose(i_xray)

for energy in [300]:
    e_beam = ElectronBeam(energy)
    i_matrix, i_feature = simulator.get_e_categories(e_beam, measure, output_e)
    simulator.get_e_dose(i_matrix, i_feature)

simulator.plot_dose()


