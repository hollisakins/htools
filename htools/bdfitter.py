from astropy.io import fits
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from htools.utils import get_fnu, get_lameff
import glob, sys, os
from dotmap import DotMap
from scipy.stats import binned_statistic

sonora_path = '/Users/hba423/Drive/Research/Proposals/JWST/Cycle3/ERO_AGN_spectra/sonora-spectra/'

def compute_sonora_photometry():
    filepaths = glob.glob(sonora_path+'*')
    filepaths = [f for f in filepaths if not 'parameters' in f]
    filepaths = [f for f in filepaths if not 'photometry' in f]
    names = [f.split('/')[-1].split('_')[1] for f in filepaths]
    temperatures = np.array([f.split('t')[1].split('g')[0] for f in names],dtype=int)
    gravities = np.array([f.split('g')[1].split('n')[0] for f in names],dtype=int)
    filters = ['hst_acs_f606w','hst_acs_f814w','jwst_nircam_f090w','jwst_nircam_f115w','jwst_nircam_f150w','jwst_nircam_f200w','jwst_nircam_f277w','jwst_nircam_f356w','jwst_nircam_f410m','jwst_nircam_f444w','jwst_miri_f770w','jwst_miri_f1800w']
    import tqdm
    phot_mod = np.zeros((len(filepaths), len(filters)))
    for i in tqdm.tqdm(range(len(filepaths))):
        filepath = filepaths[i]
        lam, fnu = np.loadtxt(filepath, skiprows=2).T
        phot_mod[i,:] = get_fnu(lam, fnu, filters)
        
        
    names = ['sp_'+name for name in names]
    print(names)
    columns = []
    columns.append(fits.Column(name='T', format='K', unit='K', array=temperatures))
    columns.append(fits.Column(name='g', format='K', unit='cgs', array=gravities))
    columns.append(fits.Column(name='name', format='A20', array=names))
    for i,filt in enumerate(filters):
        columns.append(fits.Column(name=filt, format='D', unit='uJy', array=phot_mod[:,i]))

    t = fits.BinTableHDU.from_columns(fits.ColDefs(columns))
    t.writeto(f'{sonora_path}/sonora_photometry.fits', overwrite=True)


    #### The flux received at Earth is that given in the table scaled by (R/D)2  where R is the radius of the object (given in the companion evolution tables) and D its distance. 



# from astropy.io import fits

# ID = 35666
# input_catalog = '/Users/hba423/Drive/Data/PRIMER/catalogs/primer-cosmos-grizli-v0.2.fits'
# catalog = fits.open(input_catalog)[1].data
# cat = catalog[catalog['ID']==ID]
# phot = np.array([cat[f'FLUX_KRON_{b}'][0] for b in ['F606W','F814W','F090W','F115W','F150W','F200W','F277W','F356W','F410M','F444W','F770W','F1800W']])
# phot_err = np.array([cat[f'FLUX_ERR_KRON_{b}'][0] for b in ['F606W','F814W','F090W','F115W','F150W','F200W','F277W','F356W','F410M','F444W','F770W','F1800W']])
# filters = ['hst_acs_f606w','hst_acs_f814w','jwst_nircam_f090w','jwst_nircam_f115w','jwst_nircam_f150w','jwst_nircam_f200w','jwst_nircam_f277w','jwst_nircam_f356w','jwst_nircam_f410m','jwst_nircam_f444w','jwst_miri_f770w','jwst_miri_f1800w']
# from htools.utils import get_fnu, get_lameff
# phot_lam = get_lameff(filters)

def fit_sonora(phot, phot_err, filters, Nbins=1000):
    try:
        sonora_phot = fits.getdata(f'{sonora_path}/sonora_photometry.fits')
    except:
        print(f'File "sonora_photometry.fits" not found at {sonora_path}')
        raise

    phot_mod = np.array([sonora_phot[filt] for filt in filters]).T 
    C, chi2 = np.zeros(len(sonora_phot)), np.zeros(len(sonora_phot))
    for i in range(len(sonora_phot)):
        C[i] = np.sum(phot_mod[i,:] * phot/phot_err**2) / np.sum(phot_mod[i,:]**2/phot_err**2)
        chi2[i] = np.sum(((phot-phot_mod[i,:]*C[i])/phot_err)**2)

    i = np.argmin(chi2)
    temperatures = sonora_phot['T']
    gravities = sonora_phot['g']
    
    filepath = f"{sonora_path}/{sonora_phot['name'][i]}_m0.0"
    lam, fnu = np.loadtxt(filepath, skiprows=2).T
    lam_bins = np.logspace(np.log10(np.min(lam)),np.log10(np.max(lam)), Nbins)
    fnu, _, _ = binned_statistic(lam, fnu, bins=lam_bins, statistic='mean')
    lam = 0.5*(lam_bins[1:]+lam_bins[:-1])
    return DotMap(lam=lam, fnu=fnu*C[i], phot=phot_mod[i,:]*C[i], chi2=chi2[i], T=temperatures[i], g=gravities[i])

# from htools.utils import get_fnu, get_lameff
# phot_lam = get_lameff(filters)



# plt.figure()
# plt.errorbar(phot_lam, phot, yerr=phot_err, linewidth=0, marker='s', ms=6, 
#             mfc='none', mec='k', elinewidth=1, ecolor='k', capthick=1, capsize=2, zorder=1000,
#             # label=f'COS-{id_classic} photometry')
#             label=f'COS-{ID} photometry')
# plt.step(lam, fnu*C[i], where='mid', color='tab:red', linewidth=0.5, alpha=0.5)#, label=fr'Sonora model, $T={fit1.T} K$, $\log g = {np.log10(fit1.g*100):.2f} cgs$, $\chi^2 = {fit1.chi2:.1f}$')
# plt.scatter(phot_lam, phot_mod[i,:]*C[i], color='tab:red', s=3)
# # plt.stairs(fit2.fnu, fit2.lam, color='royalblue', linewidth=0.5, alpha=0.5, label=fr'Sonora model, $T={fit2.T} K$, $\log g = {np.log10(fit2.g*100):.2f} cgs$, $\chi^2 = {fit2.chi2:.1f}$')
# # plt.scatter(phot_lam, fit2.phot, color='royalblue', s=3)
# plt.semilogx()
# plt.xlim(0.5, 20)
# # plt.xlim(1.05, 1.062)
# plt.ylim(-0.1, 2.5)
# plt.legend(loc='upper left')
# plt.show()

# # C[i] = np.sum(phot_mod * phot / phot_err**2) / np.sum(phot_mod**2 / phot_err**2)
# # chi2[i] = np.sum(((phot - phot_mod*C[i])/phot_err)**2)


# # i = np.argmin(chi2)
# # lam, fnu = np.loadtxt(filepaths[i], skiprows=3).T
# # lam, fnu = np.flip(lam), np.flip(fnu)
# # phot_mod = get_fnu(lam, fnu, filters)
# # # fnu *= 1e29 # convert to uJy
# # from scipy.stats import binned_statistic
# # lam_bins = np.logspace(np.log10(np.min(lam)),np.log10(np.max(lam)), Nbins)
# # fnu_binned, _, _ = binned_statistic(lam, fnu, bins=lam_bins, statistic='mean')

# # # return DotMap(lam=lam_bins, fnu=fnu_binned*C[i], phot=phot_mod*C[i], chi2=chi2[i], T=temperatures[i], g=gravities[i])




# # fit1 = fit_sonora(phot, phot_err, filters)
