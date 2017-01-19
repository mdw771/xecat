class result_holder(object):

    def __init__(self, energy, thickness, step, t, i_noscat, i_1el, i_1elpc, i_pc, i_elpl, i_elplpc, i_out, i_pi, i_inel, i_inelpc):
        self.energy = energy
        self.thickness = thickness
        self.step = step
        self.t = t
        self.i_noscat = i_noscat
        self.i_1el = i_1el
        self.i_1elpc = i_1elpc
        self.i_pc = i_pc
        self.i_elpl = i_elpl
        self.i_elplpc = i_elplpc
        self.i_out = i_out
        self.i_pi = i_pi
        self.i_inel = i_inel
        self.i_inelpc = i_inelpc