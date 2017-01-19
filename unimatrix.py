import xraylib
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from constants import *
from util import *
from results import *
import os


class x_ray_beam(object):

    def __init__(self, pixel, wd, n_ccd, matrix_compound, matrix_density):
        """
        Initialize object.
        :param pixel: Pixel size in um
        :param wd: working distance in um
        :param N: number of pixels along the side of CCD array
        :param matrix_compound: formula of matrix material
        :param matrix_density: density of matrix material in g/cm3
        """
        self.pixel = pixel
        self.wd = wd
        self.n = n_ccd
        self.results = []
        self.matrix = matrix_compound
        self.matrix_den = matrix_density
        self.matrix_elements, self.matrix_stoic, self.matrix_mw = parser(self.matrix)

    def get_xray_categories(self, energy, thickness, step):
        """
        :param energy: beam energy in keV
        :param thickness: sample thickness in cm
        :param step: step length in cm
        :return:
        """
        print("-----------------------------------------------")
        print("Energy: %.2f keV" % energy)
        print("Thickness limit: %.2f cm" % thickness)
        print("Thickness resolution: %.2f um" % (step*1e4))

        # Abbe criterion applied
        wavelen = PlanckConst * SpeedOfLight / (energy * 1000 * ECharge)
        num_aperture = (0.5 * wavelen) / (self.pixel * 1e-6)
        (eta_el, eta_inel) = scat_efficiencies(energy, self.matrix_mw, self.matrix_elements, self.matrix_stoic)

        # frac1 is the fraction in forward scattered (detectable) photons that finally enter the aperture
        frac1 = inna_frac(wavelen, num_aperture)

        # establish table of thickness in nm (must be integer)
        t = np.arange(step, thickness + 1, step)

        # frac2 is the fraction of forward scattered inelastic scattered or plural scattered photons that enter the
        # aperture and contribute to PCI background
        frac2 = num_aperture ** 2 / 2

        # retrieve cross sections in cm2/g
        cs_inel = xraylib.CS_Compt_CP(self.matrix, energy)
        cs_el = xraylib.CS_Rayl_CP(self.matrix, energy)
        cs_pi = xraylib.CS_Photo_CP(self.matrix, energy)
        cs_tot = xraylib.CS_Total_CP(self.matrix, energy)

        # probability per thickness
        k_el = cs_el * self.matrix_den
        k_inel = cs_inel * self.matrix_den
        k_elin = cs_el * (1 - eta_el) * self.matrix_den
        k_inelin = cs_inel * (1 - eta_inel) * self.matrix_den
        k_out = cs_el * eta_el * self.matrix_den + cs_inel * eta_inel * self.matrix_den
        k_pi = cs_pi * self.matrix_den
        k_tot = k_inel + k_el + k_pi

        # intensity fractions relative to primary beam
        # t is converted to cm
        i_noscat = np.exp(-k_tot * t)
        i_1el = k_elin * t * i_noscat
        i_1elpc = frac1 * i_1el
        i_pc = np.sqrt(i_noscat * i_1elpc)
        i_elpl = np.exp(-(k_inelin + k_out + k_pi) * t) - i_noscat - i_1el
        i_elplpc = i_elpl * frac2
        i_out = (k_out / (k_out + k_pi)) * (1 - np.exp(-(k_out + k_pi) * t))
        i_pi = (k_pi / (k_out + k_pi)) * (1 - np.exp(-(k_out + k_pi) * t))
        i_inel = np.exp(-(k_out + k_pi) * t) - np.exp(-(k_inelin + k_out + k_pi) * t)
        i_inelpc = i_inel * frac2

        res = result_holder(energy, thickness, step, t, i_noscat, i_1el, i_1elpc, i_pc, i_elpl, i_elplpc, i_out, i_pi, i_inel, i_inelpc)
        self.results.append(res)

        return

    def plot_all(self, dest_folder='.', dest_fname='unimatrix_xray_fig.pdf', show=False):
        """
        Plot all results
        """
        n_figs = len(self.results)
        matplotlib.rcParams['pdf.fonttype'] = 'truetype'
        fontProperties = {'family': 'serif', 'serif': ['Times New Roman'], 'weight': 'normal', 'size': 12}
        plt.rc('font', **fontProperties)

        fig, axes = plt.subplots(nrows=int(n_figs/2), ncols=2)
        fig.set_figheight(6*int(n_figs/2))
        fig.set_figwidth(14)

        for i_fig in range(n_figs):

            axes = axes.flatten()
            ax = axes[i_fig]
            res = self.results[i_fig]
            energy = res.energy
            thickness = res.thickness
            step = res.step
            t = res.t
            label_pos = int(thickness / step / 3)

            pl_pc, = ax.plot(t, res.i_pc, '--', label='Phase contrast image (PCI) signal', color='blue')
            ax.text(t[label_pos], res.i_pc[label_pos], 'Phase contrast image (PCI) signal', fontsize=9, color='blue')

            pl_elplpc, = ax.plot(t, res.i_elplpc, '--', label='Plural elastically scattered in PCI background', color='green')
            ax.text(t[label_pos], res.i_elplpc[label_pos], 'Plural elastically scattered in PCI background', fontsize=9, color='green')

            pl_inelpc, = ax.plot(t, res.i_inelpc, '--', label='Inelastically scattered in PCI background', color='red')
            ax.text(t[label_pos], res.i_inelpc[label_pos], 'Inelastically scattered in PCI background', fontsize=9, color='red')

            pl_noscat, = ax.plot(t, res.i_noscat, label='Unscattered', color='salmon')
            ax.text(t[label_pos], res.i_noscat[label_pos], 'Unscattered', fontsize=9, color='salmon')

            pl_1el, = ax.plot(t, res.i_1el, label='Single elastically scattered', color='black')
            ax.text(t[label_pos], res.i_1el[label_pos], 'Single elastically scattered', fontsize=9, color='black')

            pl_elpl, = ax.plot(t, res.i_elpl, label='Plural elastically scattered', color='darkcyan')
            ax.text(t[label_pos], res.i_elpl[label_pos], 'Plural elastically scattered', fontsize=9, color='darkcyan')

            pl_inel, = ax.plot(t, res.i_inel, label='Inelastically scattered', color='magenta')
            ax.text(t[label_pos], res.i_inel[label_pos], 'Inelastically scattered', fontsize=9, color='magenta')

            pl_out, = ax.plot(t, res.i_out, label='Scattered out', color='orange')
            ax.text(t[label_pos], res.i_out[label_pos], 'Scattered out', fontsize=9, color='orange')

            pl_pi, = ax.plot(t, res.i_pi, label='Absorbed', color='grey')
            ax.text(t[label_pos], res.i_pi[label_pos], 'Absorbed', fontsize=9, color='grey')

            ax.set_yscale('log')
            ax.set_ylim(1e-18, 1.1)
            ax.set_xlim(0, thickness)
            ax.set_xlabel('Thickness (cm)')
            ax.set_ylabel('Fraction')
            ax.set_title('(%d) Protein at %d keV' % (i_fig, energy))

        fig.savefig(os.path.join(dest_folder, dest_fname), format='pdf')

        if show:
            plt.show()

