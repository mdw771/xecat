from unimatrix import *
from itertools import izip
import constants


unimatrix_xray = x_ray_beam(pixel=1, wd=1e4, n_ccd=1024, matrix_compound='H48.6C32.9N8.9O8.9S0.6', matrix_density=1.35)

energyls = [5, 10, 20, 40]
thickls = [1, 10, 50, 100]
stepls = [1e-3, 2e-3, 4e-3, 8e-3]

for energy, thickness, step in izip(energyls, thickls, stepls):
    unimatrix_xray.get_xray_categories(energy, thickness, step)

unimatrix_xray.plot_all()