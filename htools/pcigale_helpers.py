import numpy as np

def get_cigale_prob(dir, object, param, normalize=False, downsample=False):
    fname = dir + '/' + object + '_' + param + '_chi2-block-0.npy'
    
    # Read the cigale output chi file
    chi_data = np.memmap(fname, dtype=np.float64)
    # Separate chi values from parameter values
    chi_data = chi_data.reshape( (2, int(len(chi_data)/2)) )
    # Filter out nan models
    chi_data = chi_data[:, np.isfinite(chi_data[1])]
    # Reset chi=nan to chi=1e99
    chi_data[0, np.isnan(chi_data[0])] = 1e99

    if downsample:
        chi_data = chi_data[:,::downsample]

    # Get all the unique values
    unq_vls, unq_idxs = np.unique( chi_data[1] , return_inverse=True)
    # Count the probability of each unique value
    N = len(unq_vls)
    ps = np.zeros(N)
    np.add.at(ps, unq_idxs, np.exp(-chi_data[0]/2))
    # Normalize the probability
    if normalize:
        ps = np.array(ps) / np.sum(ps)

    return unq_vls, ps

def get_cigale_prob_2d(d, obj, param1, param2, normalize=False, downsample=False):
    fname1 = d + '/' + obj + '_' + param1 + '_chi2-block-0.npy'
    fname2 = d + '/' + obj + '_' + param2 + '_chi2-block-0.npy'
    x_chi_data = np.memmap(fname1, dtype=np.float64)
    y_chi_data = np.memmap(fname2, dtype=np.float64)

    x_chi_data = x_chi_data.reshape( (2, int(len(x_chi_data)/2)) )
    y_chi_data = y_chi_data.reshape( (2, int(len(y_chi_data)/2)) )

    x_chi_data = x_chi_data[:, np.isfinite(x_chi_data[1])]
    y_chi_data = y_chi_data[:, np.isfinite(y_chi_data[1])]

    x_chi_data[0, np.isnan(x_chi_data[0])] = 1e99
    y_chi_data[0, np.isnan(y_chi_data[0])] = 1e99

    #vls = np.array([tuple(sorted([m, n])) for m, n in zip(x_chi_data[1], y_chi_data[1])])
    
    x, y = x_chi_data[1], y_chi_data[1]

    if downsample is not False:
        x = x[::downsample]
        y = y[::downsample]
        ps = np.exp(-x_chi_data[0][::downsample]/2)
    else:
        ps = np.exp(-x_chi_data[0]/2)

    if normalize:
        ps = ps/np.sum(ps)

    return x, y, ps


def get_cigale_kde(dir, object, param, logspace=True, N=1000, bw_method='scott', return_samples=True):
    #import time
    #t1 = time.time()
    x, P = get_cigale_prob(dir, object, param, normalize=True)
    if logspace:
        x = np.log10(x)

    # bins = np.arange(np.min(x), np.max(x), sep)

    # bc = (bins[1:]+bins[:-1])/2
    from scipy.stats import gaussian_kde
    # y, bins, binnumber = binned_statistic(x, P, bins=bins, statistic='sum')
    # y = y/np.sum(y)

    np.random.seed(123)
    resamples = np.random.choice(x, p=P, size=N)
    kde = gaussian_kde(resamples, bw_method=bw_method)
    #t2 = time.time()
    #print(f'getting KDE took {(t2-t1)/60:.1f} minutes with N={N}')
    if return_samples:
        return x, kde, resamples
    return x, kde


def get_cigale_kde_2d(d, obj, param1, param2, logspace1=True, logspace2=True, sep1=0.05, sep2=0.05, N=100):
    #import time
    #t1 = time.time()
    x, y, P = get_cigale_prob_2d(d, obj, param1, param2, normalize=False)
    if logspace1:
        x = np.log10(x)
    if logspace2:
        y = np.log10(y)

    xbins = np.arange(np.min(x), np.max(x), sep1)
    ybins = np.arange(np.min(y), np.max(y), sep2)

    from scipy.stats import binned_statistic_2d, gaussian_kde
    z, _, _, _ = binned_statistic_2d(y, x, P, bins=[ybins, xbins], statistic='sum')
    z = z/np.sum(z)

    # bc = (bins[1:]+bins[:-1])/2
    # resamples = np.random.choice(bc, p=y, size=N)

    

    # kde = gaussian_kde(resamples)
    #t2 = time.time()
    #print(f'getting KDE took {(t2-t1)/60:.1f} minutes with N={N}')
    return x, y, z


def get_cigale_model_sed(dir, object):
    from astropy.io import fits
    import astropy.units as u

    fname = dir + '/' + object + '_best_model.fits'
    
    with fits.open(fname) as hdu:
        data = hdu[1].data
        zbest = float(hdu[1].header['universe.redshift'])

    wavelength =(data['wavelength'] * 1e-3 * u.micron).value
    Fnu = (data['Fnu']*u.mJy).value

    return wavelength, Fnu, zbest



def get_cigale_model_fluxes(dir, object):
    from astropy.io import fits
    import astropy.units as u
    from pcigale.data import SimpleDatabase as Database

    fname1 = dir + '/results.fits'
    fname2 = dir + '/observations.fits'

    with fits.open(fname1) as hdu:
        mod = hdu[1].data
    with fits.open(fname2) as hdu:
        obs = hdu[1].data

    mod = mod[mod['id'] == object]
    obs = obs[obs['id'] == object]

    names = [f for f in obs.names if not (f.endswith('_err') or f in ['id','redshift'])]
    model_fluxes = np.zeros(len(names))
    observed_fluxes = np.zeros(len(names))

    for i, f in enumerate(names):
        model_fluxes[i] = mod['best.'+f]
        observed_fluxes[i] = obs[f]

    model_fluxes *= u.mJy
    observed_fluxes *= u.mJy
    with Database("filters") as db:
        filters = {name: db.get(name=name) for name in names}

    filters_wl = np.array([filt.pivot for filt in filters.values()]) * 1e-3 * u.micron

    return filters_wl, model_fluxes, observed_fluxes
