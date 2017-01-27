#####################################
# Input units expected:             #
# Thickness: um                     #
# Energy keV                        #
# Density: g/cm3                    #
# Intermediate quantities may use   #
# different unit systems.           #
#####################################

import sys
from math import *
from xraylib import AtomicWeight
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from matplotlib.ticker import AutoMinorLocator
import e_xsections
from constants import *
from util import *
from results import *
from unimatrix_xray import x_ray_beam
from itertools import izip
import os


class electron_beam(x_ray_beam):

    def get_xray_categories(self, energy, thickness, step):

        raise AttributeError('Invalid mathod for electron beam class.')

    def get_e_categories(self, energy, thickness, step):

        print("-----------------------------------------------")
        print("Energy: %.2f keV" % energy)
        print("Thickness limit: %.2f um" % thickness)
        print("Thickness resolution: %.2f um" % step)

        eta = 1 - 4.12 / 10

        t = np.arange(step, thickness + 1, step)

        # compute cross sections in cm2 (see Langmore 1992)
        cs_el = e_xsections.cse_elastic(energy, self.matrix_elements, self.matrix_stoic)
        cs_inel = e_xsections.cse_inelastic(self.ineleloss, energy, self.matrix_elements, self.matrix_stoic)

        # probability per thickness in um-1
        # delta is the number density of molecules (cm-3)
        delta = self.matrix_den / self.matrix_mw * Nav
        k_el = cs_el * delta * 1e-4
        k_inel = cs_inel * delta * 1e-4
        k_elin = cs_el * (1 - eta) * delta * 1e-4
        k_out = cs_el * eta * delta * 1e-4
        k_tot = k_inel + k_el

        # intensity fractions relative to primary beam
        # t is converted to cm
        i_noscat = np.exp(-k_tot * t)
        i_1el = k_elin * t * i_noscat
        i_pc = np.sqrt(i_noscat * i_1el)
        i_elpl = np.exp(-(k_inel + k_out) * t) - i_noscat - i_1el
        i_out = 1 - np.exp(-k_out * t)
        i_inel = np.exp(-k_out * t) - np.exp(-(k_inel + k_out) * t)

        res = result_holder(energy, thickness, step, t, i_noscat=i_noscat, i_1el=i_1el, i_pc=i_pc, i_elpl=i_elpl,
                            i_out=i_out, i_inel=i_inel)
        self.results.append(res)

        return

    def plot_all(self, dest_folder='.', dest_fname='unimatrix_e_fig.pdf', show=False):
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

            for i in res.categories:
                ax.plot(t, i.data, linestyle=i.style, label=i.label, color=i.color)
                ax.text(t[label_pos], i.data[label_pos], i.label, fontsize=9, color=i.color)

            ax.set_yscale('log')
            ax.set_ylim(1e-3, 1)
            ax.set_xlim(0, thickness)
            ax.set_xlabel('Thickness ($\mu$m)')
            ax.set_ylabel('Fraction')
            ax.set_title('(%d) Protein at %d keV' % (i_fig, energy))

        fig.savefig(os.path.join(dest_folder, dest_fname), format='pdf')

        if show:
            plt.show()


unimatrix_e = electron_beam(pixel=1, wd=1e4, n_ccd=1024, matrix_compound='H48.6C32.9N8.9O8.9S0.6', matrix_density=1.35,
                            ineleloss=37.5)

energyls = [300, 1000]
thickls = [2, 3]
stepls = [5e-3, 5e-3]

for energy, thickness, step in izip(energyls, thickls, stepls):
    unimatrix_e.get_e_categories(energy, thickness, step)

unimatrix_e.plot_all()