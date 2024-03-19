import astropy.units as u 
import prospect.io.read_results as reader
from sedpy.observate import load_filters
import numpy as np

def get_prospector_result(dir, objname, which):
    if which=='latest':
        import glob, os
        files = glob.glob(f'{dir}/{objname}-prospector_*_result.h5')
        hfile = max(files, key=os.path.getctime)
        print(f'Using latest prospector output {hfile}')
    else:
        hfile = f'{dir}/{objname}-prospector_{which}_result.h5'
    result, obs, _ = reader.results_from(hfile, dangerous=False)
    obs['filters'] = load_filters(obs['filters'], directory='/Users/hba423/Drive/Research/Codes/sedpy_filters/')
    return result, obs


def get_prospector_model(dir, objname, which):    
    result, obs = get_prospector_result(dir, objname, which)
    model = reader.get_model(result)
    sps = reader.get_sps(result)
    run_params = result['run_params']

    imax = np.argmax(result['lnprobability'])
    i, j = np.unravel_index(imax, result['lnprobability'].shape)
    theta_max = result['chain'][i, j, :].copy()

    phot_lam = obs['phot_lam'] * u.angstrom
    spec_lam = sps.wavelengths * (1 + run_params['object_redshift']) * u.angstrom
    phot_lam = phot_lam.to(u.um)
    spec_lam = spec_lam.to(u.um)
    spec_maggies, phot_maggies, _ = model.mean_model(theta_max, obs, sps=sps)
    spec_fnu = spec_maggies * 3631 * u.Jy
    phot_fnu_mod = phot_maggies * 3631 * u.Jy

    phot_maggies_obs = obs['maggies']
    phot_fnu_obs = phot_maggies_obs * 3631 * u.Jy

    return spec_lam, spec_fnu, phot_lam, phot_fnu_mod, phot_fnu_obs

def get_prospector_kde(dir, objname, which, param, logspace=True):
    from scipy.stats import gaussian_kde
    result, obs = get_prospector_result(dir, objname, which)

    assert param in result['theta_labels'] 
    i = result['theta_labels'].index(param)
    samples = result['chain'][:,:,i].flatten()

    kde = gaussian_kde(samples)
    
    if logspace:
        x = np.logspace(np.log10(np.min(samples)), np.log10(np.max(samples)), 1000)
    else:
        x = np.linspace(np.min(samples), np.max(samples), 100)

    return x, kde

