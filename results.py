class dset(object):

    def __init__(self, data, label, color, style=None):
        self.data = data
        self.label = label
        self.color = color
        self.style=style

class result_holder(object):

    def __init__(self, beam, output, measurement, constants=None, i_noscat=None, i_1el=None, i_1elpc=None, i_pc=None, i_elpl=None,
                 i_elplpc=None, i_out=None, i_pi=None, i_inel=None, i_inelpc=None, i_innoinel=None, i_in=None,
                 i_1elf_f=None):
        self.beam = beam
        self.measurement = measurement
        self.output = output
        self.constants = constants
        self.categories = []
        if i_noscat is not None:
            self.i_noscat = dset(i_noscat, 'Unscattered', 'salmon')
            self.categories.append(self.i_noscat)
        if i_1el is not None:
            self.i_1el = dset(i_1el, 'Single elastically scattered', 'black')
            self.categories.append(self.i_1el)
        if i_1elpc is not None:
            self.i_1elpc = dset(i_1elpc, 'Single elastically scattered in imaging signal', 'black', style='--')
            self.categories.append(self.i_1elpc)
        if i_pc is not None:
            self.i_pc = dset(i_pc, 'Phase contrasti imaging (PCI) signal', 'blue', style='--')
            self.categories.append(self.i_pc)
        if i_elpl is not None:
            self.i_elpl = dset(i_elpl, 'Plural elastically scattered', 'darkcyan')
            self.categories.append(self.i_elpl)
        if i_elplpc is not None:
            self.i_elplpc = dset(i_elplpc, 'Plural elastically scattered in PCI background', 'green', style='--')
            self.categories.append(self.i_elplpc)
        if i_out is not None:
            self.i_out = dset(i_out, 'Scattered out', 'orange')
            self.categories.append(self.i_out)
        if i_pi is not None:
            self.i_pi = dset(i_pi, 'Absorbed', 'grey')
            self.categories.append(self.i_pi)
        if i_inel is not None:
            self.i_inel = dset(i_inel, 'Inelastically scattered', 'magenta')
            self.categories.append(self.i_inel)
        if i_inelpc is not None:
            self.i_inelpc = dset(i_inelpc, 'Inelastically scattered in PCI background', 'red', style='--')
            self.categories.append(self.i_inelpc)
        if i_innoinel is not None:
            self.i_innoinel = dset(i_innoinel, 'Scattered in without inelastic scattering', 'red')
            self.categories.append(self.i_innoinel)
        if i_in is not None:
            self.i_in = dset(i_in, 'Scattered in', 'red')
            self.categories.append(self.i_in)
        if i_1elf_f is not None:
            self.i_1elf_f = dset(i_1elf_f, 'Single elastic scattered by feature', 'red')
            self.categories.append(self.i_1elf_f)