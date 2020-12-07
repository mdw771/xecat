import numpy as np

from feat_matrix import CompositeSimulator
from sample import Composite
from measurement import *
from output import *


# grimm_jmic_1996 (1)
sample = Composite(matrix_compound='H2O', matrix_density=0.92, feature_compound='H48.6C32.9N8.9O8.9S0.6',
                   feature_density=1.35, matrix_eloss=39.3, feature_eloss=37.5,
                   matrix_thickness=0.1,
                   feature_thickness=0.004,
                   variable='none')
measure = Measurement(pixel_size=0.004, n_ccd=1024, working_distance=1e4)
simulator = CompositeSimulator(sample)
output_e = Output(sample, step=None, overflow_limit=11)
e_beam = ElectronBeam(120)
i_matrix, i_feature = simulator.get_e_categories(e_beam, measure, output_e)
simulator.get_e_dose(i_matrix, i_feature)

# grimm_jmic_1996 (2)
sample = Composite(matrix_compound='H2O', matrix_density=0.92, feature_compound='H48.6C32.9N8.9O8.9S0.6',
                   feature_density=1.35, matrix_eloss=39.3, feature_eloss=37.5,
                   matrix_thickness=0.1,
                   feature_thickness=0.0007,
                   variable='none')
measure = Measurement(pixel_size=0.0007, n_ccd=1024, working_distance=1e4)
simulator = CompositeSimulator(sample)
output_e = Output(sample, step=None, overflow_limit=11)
e_beam = ElectronBeam(120)
i_matrix, i_feature = simulator.get_e_categories(e_beam, measure, output_e)
simulator.get_e_dose(i_matrix, i_feature)

# grimm_jmic_1996 (3)
sample = Composite(matrix_compound='H2O', matrix_density=0.92, feature_compound='H48.6C32.9N8.9O8.9S0.6',
                   feature_density=1.35, matrix_eloss=39.3, feature_eloss=37.5,
                   matrix_thickness=0.1,
                   feature_thickness=0.0004,
                   variable='none')
measure = Measurement(pixel_size=0.0004, n_ccd=1024, working_distance=1e4)
simulator = CompositeSimulator(sample)
output_e = Output(sample, step=None, overflow_limit=11)
e_beam = ElectronBeam(120)
i_matrix, i_feature = simulator.get_e_categories(e_beam, measure, output_e)
simulator.get_e_dose(i_matrix, i_feature)

# schroder_jsb_1990
sample = Composite(matrix_compound='H2O', matrix_density=0.92, feature_compound='H48.6C32.9N8.9O8.9S0.6',
                   feature_density=1.35, matrix_eloss=39.3, feature_eloss=37.5,
                   matrix_thickness=0.05,
                   feature_thickness=0.002,
                   variable='none')
measure = Measurement(pixel_size=0.002, n_ccd=1024, working_distance=1e4)
simulator = CompositeSimulator(sample)
output_e = Output(sample, step=None, overflow_limit=11)
e_beam = ElectronBeam(80)
i_matrix, i_feature = simulator.get_e_categories(e_beam, measure, output_e)
simulator.get_e_dose(i_matrix, i_feature)

# nicastro_jsb_2000
sample = Composite(matrix_compound='H2O', matrix_density=0.92, feature_compound='H48.6C32.9N8.9O8.9S0.6',
                   feature_density=1.35, matrix_eloss=39.3, feature_eloss=37.5,
                   matrix_thickness=0.6,
                   feature_thickness=0.007,
                   variable='none')
measure = Measurement(pixel_size=0.007, n_ccd=1024, working_distance=1e4)
simulator = CompositeSimulator(sample)
output_e = Output(sample, step=None, overflow_limit=11)
e_beam = ElectronBeam(120)
i_matrix, i_feature = simulator.get_e_categories(e_beam, measure, output_e)
simulator.get_e_dose(i_matrix, i_feature)