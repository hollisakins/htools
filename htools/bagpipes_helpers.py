import os
import numpy as np

import h5py

def get_spectrum(h5path):
    f = dd.io.load(h5path)



delayed = {}
delayed["massformed"] = (7.5, 11.5)   # Log_10 total stellar mass formed: M_Solar
delayed["metallicity"] = (0., 2.5)  # Metallicity: Z_sol = 0.02
delayed["age"] = (0., 14)           # Time since SF began: Gyr
delayed["tau"] = (0., 10)           # Timescale of decrease: Gyr
nebular = {}
nebular["logU"] = (-3,-1)          # Log_10 of the ionization parameter.
dust = {}
dust["type"] = 'Calzetti'   # Attenuation law: "Calzetti", "Cardelli", "CF00" or "Salim"
dust["Av"] = (0., 4.)     # Absolute attenuation in the V band: magnitudes
fit_instructions = {}
fit_instructions["delayed"] = delayed
fit_instructions["nebular"] = nebular
fit_instructions["dust"] = dust
fit_instructions["redshift"] = (0., 15.)  # Vary observed redshift from 0 to 15


galaxy = bp.model_galaxy('CWeb-Mz8', load_data, filt_list=filters, photometry_exists=True, phot_units='mujy', spectrum_exists=False)