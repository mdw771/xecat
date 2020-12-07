#####################################
# Input units expected:             #
# Thickness: um                     #
# Energy keV                        #
# Density: g/cm3                    #
# Intermediate quantities may use   #
# different unit systems.           #
#####################################

import xraylib
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from constants import *
from util import *
from results import *
from sample import SingleMaterial
from beam import *
from output import *
from measurement import *
from itertools import izip
import os


class XraySingleMatSimulator(object):

    def __init__(self, sample=None):
        """
        Initialize object.
        """
        self.results = []
        if sample is not None:
            self.sample = sample

    def get_xray_categories(self, xray_beam, measurement, output, new_sample=None):
        """
        :param energy: beam energy in keV
        :param thickness: sample thickness in cm
        :param step: step length in cm
        :return:
        """
        if new_sample is not None:
            self.sample = sample
        thickness = self.sample.thickness
        step = output.step
        energy = xray_beam.energy

        print("-----------------------------------------------")
        print("Energy: %.2f keV" % energy)
        print("Thickness limit: %.2f um" % thickness)
        print("Thickness resolution: %.2f um" % (step))

        # Abbe criterion applied
        wavelen = xray_beam.wavelength
        num_aperture = measurement.get_numerical_aperture(xray_beam)
        (eta_el, eta_inel) = scat_efficiencies(energy, self.sample.mw, self.sample.elements, self.sample.stoic)

        # frac1 is the fraction in forward scattered (detectable) photons that finally enter the aperture
        frac1 = inna_frac(wavelen, num_aperture)

        # frac2 is the fraction of forward scattered inelastic scattered or plural scattered photons that enter the
        # aperture and contribute to PCI background
        frac2 = num_aperture ** 2 / 2

        # retrieve cross sections in cm2/g
        cs_inel = xraylib.CS_Compt_CP(self.sample.compound, energy)
        cs_el = xraylib.CS_Rayl_CP(self.sample.compound, energy)
        cs_pi = xraylib.CS_Photo_CP(self.sample.compound, energy)
        cs_tot = xraylib.CS_Total_CP(self.sample.compound, energy)

        # probability per thickness
        k_el = cs_el * self.sample.density * 1e-4
        k_inel = cs_inel * self.sample.density * 1e-4
        k_elin = cs_el * (1 - eta_el) * self.sample.density * 1e-4
        k_inelin = cs_inel * (1 - eta_inel) * self.sample.density * 1e-4
        k_out = (cs_el * eta_el * self.sample.density + cs_inel * eta_inel * self.sample.density) * 1e-4
        k_pi = cs_pi * self.sample.density * 1e-4
        k_tot = k_inel + k_el + k_pi

        # intensity fractions relative to primary beam
        t = output.t
        i_noscat = np.exp(-k_tot * t)
        i_1el = k_elin * t * i_noscat
        i_1elpc = frac1 * i_1el
        i_pc = i_noscat + i_1elpc
        i_elpl = np.exp(-(k_inelin + k_out + k_pi) * t) - i_noscat - i_1el
        i_elplpc = i_elpl * frac2
        i_out = (k_out / (k_out + k_pi)) * (1 - np.exp(-(k_out + k_pi) * t))
        i_pi = (k_pi / (k_out + k_pi)) * (1 - np.exp(-(k_out + k_pi) * t))
        i_inel = np.exp(-(k_out + k_pi) * t) - np.exp(-(k_inelin + k_out + k_pi) * t)
        i_inelpc = i_inel * frac2

        # res = result_holder(xray_beam, output, measurement, i_noscat=i_noscat, i_1el=i_1el, i_1elpc=None, i_pc=i_pc,
        #                     i_elpl=i_elpl, i_elplpc=i_elplpc, i_out=i_out, i_pi=i_pi, i_inel=i_inel, i_inelpc=i_inelpc,
        #                     i_df=i_1elpc)
        res = result_holder(xray_beam, output, measurement, i_noscat=i_noscat, i_1el=i_1el,
                            i_elpl=i_elpl, i_out=i_out, i_pi=i_pi, i_inel=i_inel)
        self.results.append(res)

        return

    def plot_all(self, dest_folder='.', dest_fname='unimatrix_xray_fig.eps', show=False):
        """
        Plot all results
        """
        n_figs = len(self.results)
        matplotlib.rcParams['pdf.fonttype'] = 'truetype'
        fontProperties = {'family': 'serif', 'serif': ['Times New Roman'], 'weight': 'normal', 'size': 12}
        plt.rc('font', **fontProperties)

        fig, axes = plt.subplots(nrows=int(n_figs/3), ncols=3)
        fig.set_figheight(6*int(n_figs/3))
        fig.set_figwidth(20)

        for i_fig in range(n_figs):

            axes = axes.flatten()
            ax = axes[i_fig]
            res = self.results[i_fig]
            energy = res.beam.energy
            t_cm = res.output.t * 1e-4
            thickness = t_cm[-1]
            label_pos = int(t_cm.size / 3)

            for i in res.categories:
                ax.plot(t_cm, i.data, linestyle=i.style, label=i.label, color=i.color)
                ax.text(t_cm[label_pos], i.data[label_pos], i.label, fontsize=9, color=i.color)

            ax.set_yscale('log')
            ax.set_ylim(1e-10, 1.1)
            ax.set_xlim(0, thickness)
            ax.set_xlabel('Thickness (cm)')
            ax.set_ylabel('Fraction')
            ax.set_title('(%d) EPON at %d keV' % (i_fig+1, energy))

        fig.savefig(os.path.join(dest_folder, dest_fname), format='pdf')

        if show:
            plt.show()


energyls = [5, 15, 45]
thickls = [0.1e4, 2e4, 10e4]
stepls = [1, 1, 1]

measurement = Measurement(pixel_size=1, n_ccd=1024, working_distance=1e4)
simulator = XraySingleMatSimulator()

# EPON formula http://www.sigmaaldrich.com/catalog/product/aldrich/181196?lang=en&region=US
for (energy, thickness, step) in izip(energyls, thickls, stepls):

    sample = SingleMaterial(compound='H2O', density=0.94, thickness=thickness)
    xray_beam = XrayBeam(energy)
    output = Output(sample, step=step)
    simulator.get_xray_categories(xray_beam, measurement, output, new_sample=sample)

simulator.plot_all(show=False, dest_fname='unimatrix_fig_x_ice.pdf')

