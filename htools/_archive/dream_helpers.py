import htools
import os
import numpy as np
from fsps import StellarPopulation
from astropy.io import fits
import astropy.units as u


def load_intrinsic(path_to_intrinsic, which='all'):
    cat = fits.open(path_to_intrinsic)[1].data
    if which != 'all':
        cat = cat[cat['ID']==which]
    return cat

def load_true_spectrum(which, path_to_intrinsic='/Users/hba423/Drive/Research/COSMOS-Web/DREaM/data/DREaM_intrinsic.fits'):
    '''Generates the intrinsic galaxy spectrum with fsps. 
       Returns observed-frame wavelength in microns and flux in Jy'''
    cat = load_intrinsic(path_to_intrinsic, which=which)
    ztrue = cat['redshift']
    M_gal = cat['M_gal']
    SF = cat['SF']
    t_start = cat['t_start']
    tau = cat['tau']
    logZ = cat['logZ']
    dust = cat['dust']
    logUS = cat['logUS']

    sp = StellarPopulation(imf_type=1, zcontinuous=1, add_igm_absorption=True, add_neb_emission=True, zred=ztrue, 
                        logzsol=logZ, gas_logu=logUS, gas_logz=logZ, sfh=4, sf_start=t_start, add_stellar_remnants=False,
                        add_agb_dust_model=True, tau=tau, tburst=14, dust_type=2, dust1=0., dust2=dust)

    wave, spec = sp.get_spectrum(tage=htools.cosmo.age(ztrue).value)
    spec = spec * (1+ztrue) * u.Lsun/u.Hz
    spec = spec/(4*np.pi*htools.cosmo.luminosity_distance(ztrue)**2)
    spec = spec/sp.stellar_mass * M_gal 
    spec = spec.to(u.Jy).value
    wave = wave/1e4*(1+ztrue)

    return wave, spec


def load_photometry(filtlist, which='all', path_to_photo1='/Users/hba423/Drive/Research/COSMOS-Web/DREaM/data/DREaM_photo.fits', path_to_photo2='/Users/hba423/Drive/Research/COSMOS-Web/DREaM/data/DREaM_photo2.fits'):
    '''Loads mock photometry from the two separate photometry catalogs. 
       Returns an array with magnitudes for each filter in filtlist. 
       First column is the galaxy ID number. 
       If which=='all', each galaxy gets its own row. '''

    cat1 = fits.open(path_to_photo1)[1].data
    cat2 = fits.open(path_to_photo2)[1].data

    if which != 'all':
        cat1 = cat1[cat1['ID']==which]
        cat2 = cat2[cat2['ID']==which]

    outmags = np.zeros((len(cat1),len(filtlist)+1))
    outmags[:,0] = cat1['ID']
    for j,f in enumerate(filtlist):
        if f in cat1.dtype.fields.keys():
            fmag = cat1[f]
        elif f in cat2.dtype.fields.keys():
            fmag = cat2[f]
        outmags[:,j+1] = fmag

    return outmags
        




