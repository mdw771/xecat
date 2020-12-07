#####################################
# Input units expected:             #
# Thickness: um                     #
# Energy keV                        #
# Density: g/cm3                    #
# Intermediate quantities may use   #
# different unit systems.           #
#####################################

import sys
from xraylib import AtomicWeight
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from matplotlib.ticker import AutoMinorLocator
from constants import *
from util import *
from results import *
from measurement import *
from output import *
from sample import *
from beam import *
from itertools import izip
import os


class ElectronSingleMatSimulator(object):

    def __init__(self, sample=None):
        """
        Initialize object.
        """
        self.results = []
        if sample is not None:
            self.sample = sample

    def get_xray_categories(self, xray_beam, measurement, output, new_sample=None):

        raise AttributeError('Invalid mathod for electron class.')

    def get_e_categories(self, e_beam, measurement, output, new_sample=None):

        assert isinstance(e_beam, ElectronBeam)
        energy = e_beam.energy

        if new_sample is not None:
            self.sample = new_sample

        print("-----------------------------------------------")
        print("Energy: %.2f keV" % energy)
        print("Thickness limit: %.2f um" % thickness)
        print("Thickness resolution: %.2f um" % step)

        eta = 1 - 4.12 / 10

        # compute cross sections in cm2 (see Langmore 1992)
        cs_el = cse_elastic(energy, self.sample.elements, self.sample.stoic)
        cs_inel = cse_inelastic(self.sample.eloss, energy, self.sample.elements, self.sample.stoic)

        # probability per thickness in um-1
        delta = self.sample.density / self.sample.mw * Nav
        k_el = cs_el * delta * 1e-4
        k_inel = cs_inel * delta * 1e-4
        k_elin = cs_el * (1 - eta) * delta * 1e-4
        k_out = cs_el * eta * delta * 1e-4
        k_tot = k_inel + k_el

        # intensity fractions relative to primary beam
        t = output.t
        i_noscat = np.exp(-k_tot * t)
        i_1el = k_elin * t * i_noscat
        i_pc = np.sqrt(i_noscat * i_1el)
        i_elpl = np.exp(-(k_inel + k_out) * t) - i_noscat - i_1el
        i_out = 1 - np.exp(-k_out * t)
        i_inel = np.exp(-k_out * t) - np.exp(-(k_inel + k_out) * t)

        res = result_holder(e_beam, output, measurement, i_noscat=i_noscat, i_1el=i_1el, i_pc=i_pc, i_elpl=i_elpl,
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
            energy = res.beam.energy
            t = res.output.t
            thickness = t[-1]
            label_pos = int(t.size / 3)

            for i in res.categories:
                ax.plot(t, i.data, linestyle=i.style, label=i.label, color=i.color)
                ax.text(t[label_pos], i.data[label_pos], i.label, fontsize=9, color=i.color)

            ax.set_yscale('log')
            ax.set_ylim(1e-3, 1)
            ax.set_xlim(0, thickness)
            ax.set_xlabel('Thickness ($\mu$m)')
            ax.set_ylabel('Fraction')
            ax.set_title('(%d) EPON at %d keV' % (i_fig+1, energy))
            minorLocator = AutoMinorLocator()
            ax.xaxis.set_minor_locator(minorLocator)

        fig.savefig(os.path.join(dest_folder, dest_fname), format='pdf')

        if show:
            plt.show()


energyls = [100, 300]
thickls = [2, 2]
stepls = [5e-3, 5e-3]

measurement = Measurement(pixel_size=1, n_ccd=1024, working_distance=1e4)
simulator = ElectronSingleMatSimulator()

for energy, thickness, step in izip(energyls, thickls, stepls):
    e_beam = ElectronBeam(energy)
    # sample = SingleMaterial(compound='H48.6C32.9N8.9O8.9S0.6', density=1.35, thickness=thickness, eloss=38.7)
    sample = SingleMaterial(compound='C18H21O3Cl', density=1.20, thickness=thickness, eloss=38.7)
    output = Output(sample, step=step)
    simulator.get_e_categories(e_beam, measurement, output, new_sample=sample)

simulator.plot_all(dest_fname='unimatrix_fig_e_epon.pdf')