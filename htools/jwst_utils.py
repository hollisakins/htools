from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.stats import sigma_clipped_stats
from photutils.aperture import CircularAperture, aperture_photometry
from astropy.cosmology import Planck18 as cosmo
import os, sys, tqdm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from copy import copy, deepcopy
import warnings
from astropy.wcs import FITSFixedWarning
from astropy.nddata.utils import NoOverlapError
warnings.simplefilter('ignore', FITSFixedWarning)
warnings.simplefilter('ignore')
import glob
from .imaging import gen_rgb_image
from .utils import get_lameff, get_lamrange#, gen_cutout
from .eazy_helpers import get_obs_sed,get_tem_sed, get_pz, load_eazy_catalog

config = {
    'bands': {
        'cosmos-web':['f814w','f115w','f150w','f277w','f444w','f770w'],
        'primer-cosmos':['f606w','f814w','f090w','f115w','f150w','f200w','f277w','f356w','f410m','f444w','f770w','f1800w'],
        'ceers':['f606w','f814w','f115w','f150w','f200w','f277w','f356w','f410m','f444w']
    },
    'filters': {
        'cosmos-web':['hst_acs_f814w','jwst_nircam_f115w','jwst_nircam_f150w','jwst_nircam_f277w','jwst_nircam_f444w','jwst_miri_f770w'],
        'primer-cosmos':['hst_acs_f606w','hst_acs_f814w','jwst_nircam_f090w','jwst_nircam_f115w','jwst_nircam_f150w','jwst_nircam_f200w','jwst_nircam_f277w','jwst_nircam_f356w','jwst_nircam_f410m','jwst_nircam_f444w','jwst_miri_f770w','jwst_miri_f1800w']
    },
    'colors': {
        'cosmos-web':['darkmagenta','#0088e7', '#03a1a1', '#83b505','#ab0202','#e74001'],
        'primer-cosmos':['indigo', 'darkmagenta','#0c00e7', '#0088e7', '#03a1a1', '#009a49', '#83b505', '#c2a206','darkorange','#ab0202','#e74001','#630606']
    }
}


def get_filepath(field, band, ext, tile=None):
    '''
    general helper function to load a given JWST (or other) image
    some fields are broken into tiles to save memory; for these, specify e.g. tile='A5'
    `field` must be one of 'cosmos-web', 'primer-cosmos', 'ceers' (for now)
    `band` must be in the above config
    '''

    ############################ some code to change the default path based on the hostname
    import socket
    hostname = socket.gethostname()
    if hostname == 'cns-s-pmaa65432': # my desktop
<<<<<<< HEAD
        simmons_prefix = '/Users/hba423/simmons//'
=======
        simmons_prefix = '/V/simmons/'
>>>>>>> 55572f3 (test)
    else: # otherwise, on hard drive
        simmons_prefix = '/V/simmons/'
    ############################

    ############################ for ground-based images
    if field in ['cosmos-web','primer-cosmos']:
<<<<<<< HEAD
        if band in ['Y','J','H','Ks']:
            if not ext.startswith('sci') or ext.startswith('wht'): 
                raise Exception('we dont have ERR maps for UVISTA! only SCI, WHT')
            if tile.startswith('A'):
                filepath, hdu_index = f'{simmons_prefix}/data/cosmos/mosaics/cutout-HSC-{band.upper()}-9813-pdr3_dud_rev_{ext}_zp-28.09_{tile}.fits', 0
            elif tile.startswith('B'):
                filepath, hdu_index = f'{simmons_prefix}/data/cosmos/mosaics/{tile}--cutout-HSC-{band.upper()}-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
        if band in ['g','r','i','z','y']:
            if not ext.startswith('sci'):
                raise Exception('we dont have ERR or WHT maps for HSC! only SCI')
            if tile.startswith('A'):
                filepath, hdu_index = f'{simmons_prefix}/data/cosmos/mosaics/cutout-HSC-{band.upper()}-9813-pdr3_dud_rev_{ext}_zp-28.09_{tile}.fits', 0
            elif tile.startswith('B'):
                filepath, hdu_index = f'{simmons_prefix}/data/cosmos/mosaics/{tile}--cutout-HSC-{band.upper()}-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
=======
        if band.startswith('IRAC'):
            if ext == 'sci':
                filepath, hdu_index = f"{simmons_prefix}/cosmos-web/mosaics_ground/irac.{band.split('IRAC')[-1]}.mosaic.Fconv_resamp015_zp-28.09_{tile}.fits", 0
            elif ext == 'wht':
                filepath, hdu_index = f"{simmons_prefix}/cosmos-web/mosaics_ground/irac.{band.split('IRAC')[-1]}.mosaic.Fconv_resamp015_weight_zp-28.09_{tile}.fits", 0
            else:
                raise Exception('we dont have ERR maps for IRAC! only SCI,WHT')

        if band in ['NB118','Y','J','H','Ks']:
            if ext == 'sci': 
                filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/UVISTA_{band}_12_01_24_allpaw_skysub_015_dr6_rc_v1_zp-28.09_{tile}.fits', 0
            elif ext == 'wht': 
                filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/UVISTA_{band}_12_01_24_allpaw_skysub_015_dr6_rc_v1.weight_zp-28.09_{tile}.fits', 0
            else:
                raise Exception('we dont have ERR maps for UVISTA! only SCI, WHT')
            
        if band in ['g','r','i','z','y']:
            if ext in ['sci','wht']:
                if tile in ['A4', 'A5', 'A9', 'A10']:
                    if band=='g': filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/cutout-HSC-G-9813-pdr3_dud_rev-230413-130357_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='r': filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/cutout-HSC-R-9813-pdr3_dud_rev-230413-130346_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='i': filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/cutout-HSC-I-9813-pdr3_dud_rev-230413-130351_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='z': filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/cutout-HSC-Z-9813-pdr3_dud_rev-230413-130355_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='y': filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/cutout-HSC-Y-9813-pdr3_dud_rev-230413-130357_{ext}_zp-28.09_{tile}.fits', 0
                elif tile in ['A1', 'A2', 'A3', 'A8', 'A7', 'A6']:
                    if band=='g': filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/cutout-HSC-G-9813-pdr3_dud_rev-230412-135737_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='r': filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/cutout-HSC-R-9813-pdr3_dud_rev-230413-121613_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='i': filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/cutout-HSC-I-9813-pdr3_dud_rev-230413-121625_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='z': filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/cutout-HSC-Z-9813-pdr3_dud_rev-230413-121629_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='y': filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/cutout-HSC-Y-9813-pdr3_dud_rev-230413-121631_{ext}_zp-28.09_{tile}.fits', 0
                else:
                    if band=='g': filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/{tile}--cutout-HSC-G-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
                    if band=='r': filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/{tile}--cutout-HSC-R-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
                    if band=='i': filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/{tile}--cutout-HSC-I-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
                    if band=='z': filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/{tile}--cutout-HSC-Z-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
                    if band=='y': filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics_ground/{tile}--cutout-HSC-Y-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
            else:
                raise Exception('we dont have ERR maps for HSC! only SCI,WHT')
>>>>>>> 55572f3 (test)

    # if len(suffix)>0:
    #     typ += f'_{suffix}'

    if field=='cosmos-web':
        if band=='f814w':
            if ext.startswith('sci'): ext = ext.replace('sci','drz')
        if tile.startswith('B'):
            if band=='f814w':
<<<<<<< HEAD
                filepath, hdu_index = f'{simmons_prefix}/data/cosmos-web/mosaics/mosaic_cosmos_web_2024jan_30mas_tile_{tile}_hst_acs_wfc_f814w_{ext}.fits', 0
            elif band=='f770w':
                filepath, hdu_index = f'{simmons_prefix}/data/cosmos-web/mosaics/mosaic_miri_f770w_COSMOS-Web_30mas_{tile}_{ext}.fits', 1
            elif band in ['f115w','f150w','f277w','f444w']:
                filepath, hdu_index = f'{simmons_prefix}/data/cosmos-web/mosaics/mosaic_nircam_{band}_COSMOS-Web_30mas_{tile}_{ext}.fits', 0
                # if 'psfMatched' in ext: hdu_index=0
        elif tile.startswith('A'):
            if band=='f814w':
                filepath, hdu_index = f'{simmons_prefix}/data/cosmos-web/mosaics/mosaic_cosmos_web_2023apr_30mas_tile_{tile}_hst_acs_wfc_f814w_{ext}.fits', 0
            elif band=='f770w':
                filepath, hdu_index = f'{simmons_prefix}/data/cosmos-web/mosaics/mosaic_miri_f770w_COSMOS-Web_30mas_{tile}_{ext}.fits', 1           
                if tile in ['A7','A9','A10']:
                    hdu_index = 0
            elif band in ['f115w','f150w','f277w','f444w']:
                filepath, hdu_index = f'{simmons_prefix}/data/cosmos-web/mosaics/mosaic_nircam_{band}_COSMOS-Web_30mas_{tile}_{ext}.fits', 0
=======
                filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics/mosaic_cosmos_web_2024jan_30mas_tile_{tile}_hst_acs_wfc_f814w_{ext}.fits', 0
            elif band=='f770w':
                if tile in ['B1','B2','B3','B4','B5','B6','B7','B8','B9']:
                    filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics/mosaic_miri_f770w_COSMOS-Web_30mas_{tile}_v0_7_{ext}.fits', 1
                if tile in ['B10']:
                    filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics/mosaic_miri_f770w_COSMOS-Web_30mas_{tile}_v0_6_{ext}.fits', 1
            elif band in ['f115w','f150w','f277w','f444w']:
                filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics/mosaic_nircam_{band}_COSMOS-Web_30mas_{tile}_v0_8_{ext}.fits', 0
                # if 'psfMatched' in ext: hdu_index=0
        elif tile.startswith('A'):
            if band=='f814w':
                filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics/mosaic_cosmos_web_2023apr_30mas_tile_{tile}_hst_acs_wfc_f814w_{ext}.fits', 0
            elif band=='f770w':
                if tile in ['A1']:
                    filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics/mosaic_miri_f770w_COSMOS-Web_30mas_{tile}_v0_7_{ext}.fits', 1           
                if tile in ['A2','A3','A4','A5','A6','A8']:
                    filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics/mosaic_miri_f770w_COSMOS-Web_30mas_{tile}_v0_6_{ext}.fits', 1           
                if tile in ['A7','A9','A10']:
                    filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics/mosaic_miri_f770w_COSMOS-Web_30mas_{tile}_v0_3_{ext}.fits', 0           
            elif band in ['f115w','f150w','f277w','f444w']:
                filepath, hdu_index = f'{simmons_prefix}/cosmos-web/mosaics/mosaic_nircam_{band}_COSMOS-Web_30mas_{tile}_v0_8_{ext}.fits', 0
>>>>>>> 55572f3 (test)
                # if 'psfMatched' in ext: hdu_index=0

        else:
            print(f'tile "{tile}" not understood')
            return
                
    
    elif field=='primer-cosmos':
        if band in ['f606w','f814w']:
            if ext.startswith('sci'): ext = ext.replace('sci','drz')
<<<<<<< HEAD
            filepath, hdu_index = f'{simmons_prefix}/data/primer/mosaics/mosaic_cosmos_primer_30mas_acs_wfc_{band}_{ext}.fits', 1
        elif band in ['f125w','f160w']:
            if ext.startswith('sci'): ext = ext.replace('sci','drz')
            filepath, hdu_index = f'{simmons_prefix}/data/primer/mosaics/mosaic_cosmos_primer_30mas_wfc3_ir_{band}_{ext}.fits', 1
        elif band in ['f770w','f1800w']:
            filepath, hdu_index = f'{simmons_prefix}/data/primer/mosaics/mosaic_miri_{band}_PRIMER-COSMOS_60mas_{ext}.fits', 1
        elif band in ['f090w','f115w','f150w','f200w','f277w','f356w','f410m','f444w']:
            filepath, hdu_index = f'{simmons_prefix}/data/primer/mosaics/mosaic_nircam_{band}_PRIMER-COSMOS_30mas_{ext}.fits', 1
=======
            filepath, hdu_index = f'{simmons_prefix}/primer/mosaics/mosaic_cosmos_primer_30mas_acs_wfc_{band}_{ext}.fits', 0
        elif band in ['f125w','f160w']:
            if ext.startswith('sci'): ext = ext.replace('sci','drz')
            filepath, hdu_index = f'{simmons_prefix}/primer/mosaics/mosaic_cosmos_primer_30mas_wfc3_ir_{band}_{ext}.fits', 1
        elif band in ['f770w','f1800w']:
            filepath, hdu_index = f'{simmons_prefix}/primer/mosaics/mosaic_miri_{band}_PRIMER-COSMOS_30mas_{ext}.fits', 0
        elif band in ['f090w','f115w','f150w','f200w','f277w','f356w','f410m','f444w']:
            filepath, hdu_index = f'{simmons_prefix}/primer/mosaics/mosaic_nircam_{band}_PRIMER-COSMOS_30mas_{ext}.fits', 0
>>>>>>> 55572f3 (test)

    elif field=='ceers':
        if band in ['f606w','f814w','f105w','f125w','f140w','f160w']:
            if ext=='wht': raise Exception("WHT map don't exist for CEERS HST")
            if ext=='sci_pfMatched': raise Exception("I don't yet have PSF matched mosaics for CEERS")
            
            if band in ['f606w','f814w']:
<<<<<<< HEAD
                if ext=='sci': filepath, hdu_index = f'{simmons_prefix}/data/ceers/mosaics/{tile}/egs_all_acs_wfc_{band}_030mas_v1.9_{tile}_mbkgsub1.fits', 0
                if ext=='err': filepath, hdu_index = f'{simmons_prefix}/data/ceers/mosaics/{tile}/egs_all_acs_wfc_{band}_030mas_v1.9_{tile}_rms.fits', 0
            elif band in ['f105w','f125w','f140w','f160w']:
                if ext=='sci': filepath, hdu_index = f'{simmons_prefix}/data/ceers/mosaics/{tile}/egs_all_wfc3_ir_{band}_030mas_v1.9_{tile}_mbkgsub1.fits', 0
                if ext=='err': filepath, hdu_index = f'{simmons_prefix}/data/ceers/mosaics/{tile}/egs_all_wfc3_ir_{band}_030mas_v1.9.1_{tile}_rms.fits', 0
=======
                if ext=='sci': filepath, hdu_index = f'{simmons_prefix}/ceers/mosaics/{tile}/egs_all_acs_wfc_{band}_030mas_v1.9_{tile}_mbkgsub1.fits', 0
                if ext=='err': filepath, hdu_index = f'{simmons_prefix}/ceers/mosaics/{tile}/egs_all_acs_wfc_{band}_030mas_v1.9_{tile}_rms.fits', 0
            elif band in ['f105w','f125w','f140w','f160w']:
                if ext=='sci': filepath, hdu_index = f'{simmons_prefix}/ceers/mosaics/{tile}/egs_all_wfc3_ir_{band}_030mas_v1.9_{tile}_mbkgsub1.fits', 0
                if ext=='err': filepath, hdu_index = f'{simmons_prefix}/ceers/mosaics/{tile}/egs_all_wfc3_ir_{band}_030mas_v1.9.1_{tile}_rms.fits', 0
>>>>>>> 55572f3 (test)
        elif band in ['f115w','f150w','f200w','f277w','f356w','f410m','f444w']:
            if ext=='sci': hdu_index=1
            if ext=='err': hdu_index=3
            if ext=='wht': hdu_index=5
<<<<<<< HEAD
            if ext=='sci_pfMatched': raise Exception("I don't yet have PSF matched mosaics for CEERS")            
            filepath = f'{simmons_prefix}/data/ceers/mosaics/{tile}/ceers_{tile}_{band}_v1_mbkgsub1.fits'
=======
            if ext=='sci_psfMatched': raise Exception("I don't yet have PSF matched mosaics for CEERS")            
            filepath = f'{simmons_prefix}/ceers/mosaics/{tile}/ceers_{tile}_{band}_v1_mbkgsub1.fits'
>>>>>>> 55572f3 (test)
        elif band in ['f560w','f770w']:
            if ext=='sci': hdu_index=1
            if ext=='err': hdu_index=2
            if ext=='wht': hdu_index=4
            if ext=='sci_pfMatched': raise Exception("I don't yet have PSF matched mosaics for CEERS")            
<<<<<<< HEAD
            if hostname == 'cns-s-pmaa65432':
                filepath = f'{simmons_prefix}/data/ceers/mosaics/{tile}/ceers_{tile}_{band}_i2d.fits'
            else:
                filepath = f'/data/ceers/mosaics/{tile}/ceers_{tile}_{band}_i2d.fits'
=======
            filepath = f'{simmons_prefix}/ceers/mosaics/{tile}/ceers_{tile}_{band}_i2d.fits'
>>>>>>> 55572f3 (test)


    else:
        print(f'field "{field}" not understood')
        return
    return filepath, hdu_index
    
            


def load_image(field, band, ext, tile=None, convert_units=False):
    '''
        'tile' == 'A1'--'A10' or 'B1'--'B10' or 'primer-cosmos'
        'band' == one of 'f814w', 'f115w', 'f150w', 'f277w', 'f444w', 'f770w'
        'typ' == one of 'sci', 'err', 'wht'
    '''
    filepath, hdu_index = get_filepath(field, band, ext, tile=tile)  
    f = fits.open(filepath)
    # f = fits.open(filepath)
    # if np.ndim(f[0].data)==0:
    #     data, header = f[1].data, f[1].header
    #     if 'NSCI_MU' in f[0].header:
    #         f[1].header['NSCI_MU'] = f[0].header['NSCI_MU']
    #     if 'NSCI_SIG' in f[0].header:
    #         f[1].header['NSCI_SIG'] = f[0].header['NSCI_SIG']
    # else:
        # data, header = f[0].data, f[0].header
    data, header = f[hdu_index].data, f[hdu_index].header
    del f


    if convert_units:
        if 'sci' in ext or 'err' in ext:
            if band=='f814w':
                photflam, photplam = 6.99715969242424E-20, 8047.468423484849
                conversion = 3.33564e10 * (photplam)**2 * photflam
            elif band=='f606w':
                photflam, photplam = 7.86211131043771E-20, 5921.891224579125
                conversion = 3.33564e10 * (photplam)**2 * photflam
            elif band=='f125w':
                photflam, photplam = 2.2446000E-20, 1.2486060E+04
                conversion = 3.33564e10 * (photplam)**2 * photflam
            elif band=='f160w':
                photflam, photplam = 1.9429001E-20,  15369.176
                conversion = 3.33564e10 * (photplam)**2 * photflam
            else:
                conversion = 1e12 * np.pi**2/(180**2) * (1/(3600**2)) * 0.03**2
            data *= conversion
        else:
            print('not converting units for some reason')

    return data, header

<<<<<<< HEAD
def load_hdr(field, band, tile=None):
    return load_image(field, band, 'sci', tile)[1]
=======
def load_hdr(field, band, tile=None, suffix=None):
    if suffix:
        return load_image(field, band, f'sci_{suffix}', tile)[1]
    else:
        return load_image(field, band, 'sci', tile)[1]
>>>>>>> 55572f3 (test)
def load_sci(field, band, tile=None, convert_units=False, suffix=None):
    if suffix:
        return load_image(field, band, f'sci_{suffix}', tile, convert_units=convert_units)[0]
    else:
        return load_image(field, band, 'sci', tile, convert_units=convert_units)[0]
def load_wcs(field, band, tile=None):
    return WCS(load_hdr(field, band, tile))
def load_err(field, band, tile=None, convert_units=False):
    return load_image(field, band, 'err', tile, convert_units=convert_units)[0]
def load_wht(field, band, tile=None):
    return load_image(field, band, 'wht', tile)[0]
# def load_vrp(field, band, tile=tile, suffix=''):
#     return load_image(field, tile, band, 'vrp', suffix=suffix)[0]
# def load_exp(field, band, tile=tile, suffix=''):
#     d, h = load_image(field, tile, band, 'exp', suffix=suffix)
#     return d, WCS(h)

def gen_cutout(field, band, coord, width, suffix='', tile=None, err_method='simple', verbose=False):
    '''
    coord = astropy.coordinates.SkyCoord object or tuple of (ra,dec) in degrees
    width = astropy.units.quantity.Quantity object or float assumed in arcsec
    '''
    if verbose: print(f'Generating cutout for {band}...')
    import warnings
    from astropy.wcs import FITSFixedWarning
    warnings.filterwarnings("ignore", category=FITSFixedWarning)
    from astropy.nddata.utils import Cutout2D
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    if type(coord)==tuple:
        coord = SkyCoord(*coord, unit='deg')
    if not type(width)==u.quantity.Quantity:
        width = width * u.arcsec
    
    if band=='detec':
        return gen_detec_cutout(field, tile, coord, width)
    if band=='segm':
        return gen_segm_cutout(field, tile, coord, width)

    if len(suffix)>0:
        sci, hdr = load_image(field, band, f'sci_{suffix}', tile=tile)
    else:
        sci, hdr = load_image(field, band, 'sci', tile=tile)
    wcs = WCS(hdr)
    cutout = Cutout2D(sci, coord, width, wcs=wcs)
    cutout_wcs = cutout.wcs
    sci = cutout.data

    if err_method != 'none':
        from astropy.stats import sigma_clipped_stats 
        mean, median, std = sigma_clipped_stats(sci)
        sci -= median
        if err_method=='simple': # per pixel rms 
            err = np.ones_like(sci)*std
        else:
            try:
                err = load_err(field, tile, band)
                err = Cutout2D(err, coord, width, wcs=wcs).data
            except:
                if verbose: print('No ERR map found, using WHT')
                wht  = load_wht(field, tile, band)
                err = 1/np.sqrt(Cutout2D(wht, coord, width, wcs=wcs).data)
    
    if band=='f814w':
        photflam, photplam = 6.99715969242424E-20, 8047.468423484849
        conversion = 3.33564e10 * (photplam)**2 * photflam
    elif band=='f606w':
        photflam, photplam = 7.86211131043771E-20, 5921.891224579125
        conversion = 3.33564e10 * (photplam)**2 * photflam
    elif band=='f125w':
        photflam, photplam = 2.2446000E-20, 1.2486060E+04
        conversion = 3.33564e10 * (photplam)**2 * photflam
    elif band=='f160w':
        photflam, photplam = 1.9429001E-20,  15369.176
        conversion = 3.33564e10 * (photplam)**2 * photflam
    else:
        conversion = 1e12 * np.pi**2/(180**2) * (1/(3600**2)) * 0.03**2
    sci *= conversion
    if err_method != 'none':
        err *= conversion
        setattr(cutout, 'err', err)
        setattr(cutout, 'snr', sci/err)


    cutout.data = sci
    ps = wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value
    size = np.shape(sci)[0]
    extent = [-size*ps/2, size*ps/2, -size*ps/2, size*ps/2]

    setattr(cutout, 'sci', sci)
    setattr(cutout, 'extent', extent)

    return cutout
    
def gen_detec_cutout(field, tile, coord, width):
    import warnings
    from astropy.wcs import FITSFixedWarning
    warnings.filterwarnings("ignore", category=FITSFixedWarning)
    from astropy.nddata.utils import Cutout2D
    import astropy.units as u
    import socket 
    hostname = socket.gethostname()
    if hostname=='cns-s-pmaa65432':
<<<<<<< HEAD
        simmons_prefix = '/Users/hba423/simmons//'
    else:
        simmons_prefix = '/V/simmons/'
    if field=='primer-cosmos':
        sci = fits.getdata(f'{simmons_prefix}/data/primer/mosaics/detection_chi2pos_SWLW_primer-cosmos.fits')
        wcs = WCS(fits.getheader(f'{simmons_prefix}/data/primer/mosaics/detection_chi2pos_SWLW_primer-cosmos.fits'))
    elif field=='cosmos-web':
        sci = fits.getdata(f'{simmons_prefix}/data/cosmos-web/mosaics/detection_chi2pos_SWLW_{tile}.fits')
        wcs = WCS(fits.getheader(f'{simmons_prefix}/data/cosmos-web/mosaics/detection_chi2pos_SWLW_{tile}.fits'))
=======
        simmons_prefix = '/V/simmons/'
    else:
        simmons_prefix = '/V/simmons/'
    if field=='primer-cosmos':
        sci = fits.getdata(f'{simmons_prefix}/primer/mosaics/detection_chi2pos_SWLW_primer-cosmos.fits')
        wcs = WCS(fits.getheader(f'{simmons_prefix}/primer/mosaics/detection_chi2pos_SWLW_primer-cosmos.fits'))
    elif field=='cosmos-web':
        sci = fits.getdata(f'{simmons_prefix}/cosmos-web/mosaics/detection_chi2pos_SWLW_{tile}.fits')
        wcs = WCS(fits.getheader(f'{simmons_prefix}/cosmos-web/mosaics/detection_chi2pos_SWLW_{tile}.fits'))
>>>>>>> 55572f3 (test)
    cutout = Cutout2D(sci, coord, width, wcs=wcs)
    ps = cutout.wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value
    size = np.shape(cutout.data)[0]
    extent = [-size*ps/2, size*ps/2, -size*ps/2, size*ps/2]
    setattr(cutout, 'sci', cutout.data)
    setattr(cutout, 'extent', extent)
    return cutout

def gen_segm_cutout(field, tile, coord, width):
    import warnings
    from astropy.wcs import FITSFixedWarning
    warnings.filterwarnings("ignore", category=FITSFixedWarning)
    from astropy.nddata.utils import Cutout2D
    import astropy.units as u
    import socket 
    hostname = socket.gethostname()
    if hostname=='cns-s-pmaa65432':
<<<<<<< HEAD
        simmons_prefix = '/Users/hba423/simmons//'
    else:
        simmons_prefix = '/V/simmons/'
    if field=='primer-cosmos':
        segm = np.load(f'{simmons_prefix}/data/primer/mosaics/detection_chi2pos_SWLW_primer-cosmos_{tile}_segmap.npy')
        detec = fits.getdata(f'{simmons_prefix}/data/primer/mosaics/detection_chi2pos_SWLW_primer-cosmos.fits')
        wcs = WCS(fits.getheader(f'{simmons_prefix}/data/primer/mosaics/detection_chi2pos_SWLW_primer-cosmos.fits'))
=======
        simmons_prefix = '/V/simmons/'
    else:
        simmons_prefix = '/V/simmons/'
    if field=='primer-cosmos':
        segm = np.load(f'{simmons_prefix}/primer/mosaics/detection_chi2pos_SWLW_primer-cosmos_{tile}_segmap.npy')
        detec = fits.getdata(f'{simmons_prefix}/primer/mosaics/detection_chi2pos_SWLW_primer-cosmos.fits')
        wcs = WCS(fits.getheader(f'{simmons_prefix}/primer/mosaics/detection_chi2pos_SWLW_primer-cosmos.fits'))
>>>>>>> 55572f3 (test)
        if tile == 'south':
            c0 = wcs.pixel_to_world(12480, 12500)
            cutout = Cutout2D(detec, c0, size=[11.0*u.arcmin,9.8*u.arcmin], wcs=wcs)
        elif tile == 'north': 
            c0 = wcs.pixel_to_world(12480, 33000)
            cutout = Cutout2D(detec, c0, size=[11.0*u.arcmin,9.8*u.arcmin], wcs=wcs)
        wcs = cutout.wcs
    elif field=='cosmos-web':
<<<<<<< HEAD
        segm = np.load(f'{simmons_prefix}/data/cosmos-web/mosaics/detection_chi2pos_SWLW_{tile}_segmap.npy')
        wcs = WCS(fits.getheader(f'{simmons_prefix}/data/cosmos-web/mosaics/detection_chi2pos_SWLW_{tile}.fits'))
=======
        segm = fits.getdata(f'{simmons_prefix}/cosmos-web/mosaics/detection_chi2pos_SWLW_{tile}_segmap_cleaned_v1.2.fits')
        wcs = WCS(fits.getheader(f'{simmons_prefix}/cosmos-web/mosaics/detection_chi2pos_SWLW_{tile}.fits'))
>>>>>>> 55572f3 (test)

    cutout = Cutout2D(segm, coord, width, wcs=wcs)
    ps = cutout.wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value
    size = np.shape(cutout.data)[0]
    extent = [-size*ps/2, size*ps/2, -size*ps/2, size*ps/2]
    sci = cutout.data

    for i,unq_val in enumerate(np.sort(np.unique(sci))):
        sci[sci==unq_val] = i

    cutout.data = sci
    setattr(cutout, 'sci', sci)
    setattr(cutout, 'extent', extent)
    from photutils.utils.colormaps import make_random_cmap
    from matplotlib.colors import to_rgba
    cmap = make_random_cmap(len(np.unique(sci)))
    cmap.colors[0] = to_rgba('k')
    setattr(cutout, 'cmap', cmap)
    return cutout




def make_inspection_plot(ID, 
                         nickname=None,
                         field = 'cosmos-web', # or primer-cosmos
                         catalog_path = '/data/COSMOS-Web/catalogs/',
                         catalog = 'COSMOS-Web_aperture_catalog_v0.6.fits',
                         catalog_name = 'hba_v0.6',
                         outdir='inspection_plots/', 
                         phot_columns=['auto'],
                         cutout_width = 15*u.arcsec,
                         display_width = 'auto', 
                         vmax = 'auto',
                         sepp_matching_radius = 1.0,
                         overwrite=True,
                         eazy_run=None, 
                         eazy_outdir=None, 
                         bqf_run=None, 
                         bagpipes_run=None, 
                         show_bd=False,
                         suffix=None,
                         rgb_params=dict(noiselum=0.1,noisesig=1,satpercent=0.5),
                         plot_msa_slits=False,
                         plot_msa_wavelength_coverage=False,
                         msa_target_info=None,
                         show=False, 
                         verbose=True):
                         
    with plt.style.context(('hba_sans')):
        input_catalog = f'{catalog_path}/{catalog}'
        catalog = fits.getdata(input_catalog)
        aper_diams = [0.2,0.3,0.5,0.75,1.0]

        cat = catalog[catalog['id']==ID]
        ra = cat['ra'][0]
        dec = cat['dec'][0]
        coord = SkyCoord(ra=ra, dec=dec, unit=u.deg)
        if field=='cosmos-web':
            tile = cat['tile'][0]
        elif field=='primer-cosmos':
            tile = cat['tile'][0]

        if display_width == 'auto':
            display_width = 8*np.max([cat['kron1_a'][0],0.2]) * u.arcsec

        if display_width > cutout_width:
            cutout_width = 2*display_width
            if verbose: print(f'Increasing cutout_width to {cutout_width.value:.2f}')

        if field in ['cosmos-web','primer-cosmos']:    
            secondary_catalog = f'/data/COSMOS-Web/catalogs/COSMOSWeb_JAN24_v3.0.0_assoc_cold+hot_sersic_cgs.fits'
            cat_sepp = fits.getdata(secondary_catalog)
            sepplus_coords = SkyCoord(ra=cat_sepp['RA_DETEC'], dec=cat_sepp['DEC_DETEC'], unit='deg')
            sep = coord.separation(sepplus_coords).to(u.arcsec).value
            i_sepplus = np.argmin(sep)
            if sep[i_sepplus] > sepp_matching_radius:
                if verbose: print('No match in SE++ catalog')
                i_sepplus = -99



        if verbose: print(f'{field} galaxy {ID}, tile {tile}, RA={coord.ra.value:.5f}, DEC={coord.dec.value:.5f}')
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if suffix is not None: 
            path = f"{outdir}/{field}_inspec_{catalog_name.replace('-','_')}_{ID}_{suffix}.pdf"
        else:
            path = f"{outdir}/{field}_inspec_{catalog_name.replace('-','_')}_{ID}.pdf"

        if type(phot_columns) == str:
            phot_columns = [phot_columns]
        if type(phot_columns) == tuple:
            phot_columns = list(phot_columns)

        if not overwrite: 
            if os.path.exists(path):
                if verbose: print(f'\t Skipping, plot already exists at {path}.')
                return

        if field=='cosmos-web':
            fig = plt.figure(figsize=(8.0, 9.2), constrained_layout=False)
            gs = mpl.gridspec.GridSpec(ncols=12, nrows=6, height_ratios=[2.2,0.0,1,4.5,0.5,3.0],figure=fig)
            gs.update(hspace=0.08, wspace=0.04, left=0.03, bottom=0.06, right=0.97, top=0.88)
            ax_f814 = plt.subplot(gs[0,0:2])
            ax_f115 = plt.subplot(gs[0,2:4])
            ax_f150 = plt.subplot(gs[0,4:6])
            ax_f277 = plt.subplot(gs[0,6:8])
            ax_f444 = plt.subplot(gs[0,8:10])
            ax_f770 = plt.subplot(gs[0,10:12])
            cutout_axes = [ax_f814,ax_f115,ax_f150,ax_f277,ax_f444,ax_f770]
        
            ax_segm = plt.subplot(gs[-1,3:6])
            ax_detc = plt.subplot(gs[-1,0:3])
            ax_rgb1 = plt.subplot(gs[-1,6:9])
            ax_rgb2 = plt.subplot(gs[-1,9:12])

            ax_g = plt.subplot(gs[2, 0])
            ax_r = plt.subplot(gs[2, 1])
            ax_i = plt.subplot(gs[2, 2])
            ax_z = plt.subplot(gs[2, 3])
            ax_y = plt.subplot(gs[2, 4])
            ax_Y = plt.subplot(gs[2, 5])
            ax_J = plt.subplot(gs[2, 6])
            ax_H = plt.subplot(gs[2, 7])
            ax_K = plt.subplot(gs[2, 8])
            ax_irac1 = plt.subplot(gs[2, 9])
            ax_irac3 = plt.subplot(gs[2, 10])
            ax_irac4 = plt.subplot(gs[2, 11])

            ax_sed = plt.subplot(gs[3, 1:-1])
            ax_pz = ax_sed.inset_axes([0.03,0.7,0.3,0.27])

        elif field=='primer-cosmos':
            fig = plt.figure(figsize=(8.0, 10.5), constrained_layout=False)
            gs = mpl.gridspec.GridSpec(ncols=12, nrows=7, height_ratios=[2.2,2.2,0.0,1,4.5,0.5,3.0],figure=fig)
            gs.update(hspace=0.08, wspace=0.04, left=0.03, bottom=0.05, right=0.97, top=0.88)
            ax_f606 = plt.subplot(gs[0,0:2])
            ax_f814 = plt.subplot(gs[0,2:4])
            ax_f090 = plt.subplot(gs[0,4:6])
            ax_f115 = plt.subplot(gs[0,6:8])
            ax_f150 = plt.subplot(gs[0,8:10])
            ax_f200 = plt.subplot(gs[0,10:12])
            ax_f277 = plt.subplot(gs[1,0:2])
            ax_f356 = plt.subplot(gs[1,2:4])
            ax_f410 = plt.subplot(gs[1,4:6])
            ax_f444 = plt.subplot(gs[1,6:8])
            ax_f770 = plt.subplot(gs[1,8:10])
            ax_f1800 = plt.subplot(gs[1,10:12])
            cutout_axes = [ax_f606,ax_f814,ax_f090,ax_f115,ax_f150,ax_f200,ax_f277,ax_f356,ax_f410,ax_f444,ax_f770,ax_f1800]
        
            ax_segm = plt.subplot(gs[-1,3:6])
            ax_detc = plt.subplot(gs[-1,0:3])
            ax_rgb1 = plt.subplot(gs[-1,6:9])
            ax_rgb2 = plt.subplot(gs[-1,9:12])

            ax_g = plt.subplot(gs[3, 0])
            ax_r = plt.subplot(gs[3, 1])
            ax_i = plt.subplot(gs[3, 2])
            ax_z = plt.subplot(gs[3, 3])
            ax_y = plt.subplot(gs[3, 4])
            ax_Y = plt.subplot(gs[3, 5])
            ax_J = plt.subplot(gs[3, 6])
            ax_H = plt.subplot(gs[3, 7])
            ax_K = plt.subplot(gs[3, 8])
            ax_irac1 = plt.subplot(gs[3, 9])
            ax_irac3 = plt.subplot(gs[3, 10])
            ax_irac4 = plt.subplot(gs[3, 11])

            ax_sed = plt.subplot(gs[4, 1:-1])
            ax_pz = ax_sed.inset_axes([0.03,0.7,0.3,0.27])

        ax_g.set_title('HSC $g$', fontsize=7, pad=0.5)
        ax_r.set_title('HSC $r$', fontsize=7, pad=0.5)
        ax_i.set_title('HSC $i$', fontsize=7, pad=0.5)
        ax_z.set_title('HSC $z$', fontsize=7, pad=0.5)
        ax_y.set_title('HSC $y$', fontsize=7, pad=0.5)
        ax_Y.set_title('UVISTA $Y$', fontsize=7, pad=0.5)
        ax_J.set_title('UVISTA $J$', fontsize=7, pad=0.5)
        ax_H.set_title('UVISTA $H$', fontsize=7, pad=0.5)
        ax_K.set_title('UVISTA $K$', fontsize=7, pad=0.5)
        ax_irac1.set_title('IRAC1', fontsize=7, pad=0.5)
        ax_irac3.set_title('IRAC3', fontsize=7, pad=0.5)
        ax_irac3.set_title('IRAC3', fontsize=7, pad=0.5)
        ax_irac4.set_title('IRAC4', fontsize=7, pad=0.5)

    

        cmap1 = plt.colormaps['Greys_r']
        cmap1.set_bad('0.7')
        cmap2 = plt.colormaps['Reds']
        cmap2 = cmap2(np.arange(cmap2.N))
        cmap2[:len(cmap2)//10,-1] = np.linspace(0, 1, len(cmap2)//10)
        cmap2 = mpl.colors.ListedColormap(cmap2)

        ### make cutouts
        if verbose: print('\t Making cutouts...')
        cutout_names = config['bands'][field]
        colors = config['colors'][field]
        N = len(cutout_names)
        cutouts = []
        vm = []
        for i in range(N):
            try:
                cutout = gen_cutout(field, cutout_names[i], coord, width=cutout_width, tile=tile)
                cutouts.append(cutout)
                vm.append(np.nanpercentile(cutout.snr, 99))
            except:
                cutouts.append(None)
                vm.append(0)
        
        if vmax=='auto':
            vmax = np.nanmax(vm)
            vmin = -0.15*vmax
        else:
            assert (type(vmax)==int) or (type(vmax)==float)
            vmin = -3

        for i in range(N):
            ############################## plot cutouts ##############################
            name = cutout_names[i]
            ax = cutout_axes[i]
            cutout = cutouts[i]
            if cutout is not None:
<<<<<<< HEAD
                if np.all(np.isnan(cutout.snr)):
                    ax.imshow(np.zeros_like(cutout.snr), vmin=vmin,vmax=vmax, cmap='Greys', origin='lower', extent=cutout.extent)
                else:
                    ax.imshow(cutout.snr, vmin=vmin, vmax=vmax, cmap='Greys', origin='lower', extent=cutout.extent)
                a, b, theta = cat['kron1_a'], cat['kron1_b'], np.degrees(cat['kron_pa'])
=======
                if np.all(np.isnan(cutout.snr)) or np.all(cutout.snr==0):
                    ax.imshow(np.zeros_like(cutout.snr), vmin=vmin,vmax=vmax, cmap='Greys', origin='lower', extent=cutout.extent)
                else:
                    ax.imshow(cutout.snr, vmin=vmin, vmax=vmax, cmap='Greys', origin='lower', extent=cutout.extent)
                a, b, theta = cat['kron1_a'], cat['kron1_b'], np.degrees(cat['theta_image'])
>>>>>>> 55572f3 (test)
                ax.add_patch(mpl.patches.Ellipse((0,0), width=2*a, height=2*b, angle=theta, facecolor='none', edgecolor='salmon', linestyle='-', linewidth=0.5))
                ax.add_patch(mpl.patches.Circle((0,0), radius=0.15, facecolor='none', edgecolor='w', linestyle='--', linewidth=0.5))
                if name.startswith('f'):
                    ax.set_xlim(-0.5*display_width.to(u.arcsec).value, 0.5*display_width.to(u.arcsec).value)
                    ax.set_ylim(-0.5*display_width.to(u.arcsec).value, 0.5*display_width.to(u.arcsec).value)
                else:
                    ax.set_xlim(-display_width.to(u.arcsec).value, display_width.to(u.arcsec).value)
                    ax.set_ylim(-display_width.to(u.arcsec).value, display_width.to(u.arcsec).value)
            else:
                ax.imshow(np.zeros((100,100)), vmin=vmin,vmax=vmax, cmap='Greys')
            # ax.imshow(dist, extent=cutout.extent)
            ax.axis('off')

        
        cutout_names2 = ['g','r','i','z','y', 'Y','J','H','Ks','IRAC1','IRAC3','IRAC4']
        cutout_axes2 = [ax_g,ax_r,ax_i,ax_z,ax_y,ax_Y,ax_J,ax_H,ax_K,ax_irac1,ax_irac3,ax_irac4]
        N = len(cutout_names2)
        for i in range(N):
            ############################## plot cutouts ##############################
            name = cutout_names2[i]
            ax = cutout_axes2[i]
            try:
                cutout = gen_cutout(field, name, coord, cutout_width, err_method='simple', tile=tile)
                ax.imshow(cutout.snr, vmin=-3, vmax=8, cmap='Greys', origin='lower', extent=cutout.extent)
                ax.add_patch(mpl.patches.Circle((0,0), radius=0.15, facecolor='none', edgecolor='w', linestyle='--', linewidth=0.5))
            except:
                ax.imshow(np.zeros((100,100)), vmin=-3, vmax=8, cmap='Greys', origin='lower', extent=[-display_width.to(u.arcsec).value, display_width.to(u.arcsec).value, -display_width.to(u.arcsec).value, display_width.to(u.arcsec).value])
                pass
            ax.set_xlim(-2*display_width.to(u.arcsec).value, 2*display_width.to(u.arcsec).value)
            ax.set_ylim(-2*display_width.to(u.arcsec).value, 2*display_width.to(u.arcsec).value)
            ax.axis('off')
            ax.set_aspect('equal')


        if verbose: print('\t Making RGB image...')
        if field=='cosmos-web':
            f115w_cutout_pm = gen_cutout(field, 'f115w', coord, cutout_width, tile=tile, suffix='psfMatched')
            f150w_cutout_pm = gen_cutout(field, 'f150w', coord, cutout_width, tile=tile, suffix='psfMatched')
            f277w_cutout_pm = gen_cutout(field, 'f277w', coord, cutout_width, tile=tile, suffix='psfMatched')
            f444w_cutout_pm = gen_cutout(field, 'f444w', coord, cutout_width, tile=tile)
            pm115 = f115w_cutout_pm.sci
            pm150 = f150w_cutout_pm.sci
            pm277 = f277w_cutout_pm.sci
            pm444 = f444w_cutout_pm.sci
            pm115[np.isnan(pm115)] = 0
            pm150[np.isnan(pm150)] = 0
            pm277[np.isnan(pm277)] = 0
            pm444[np.isnan(pm444)] = 0

            input_dict = {}
            input_dict['f115w'] = {'colors':np.array([0.1, 0.3, 1.0]), 'data':pm115}
            input_dict['f150w'] = {'colors':np.array([0, 1.0, 0]),     'data':pm150}
            input_dict['f277w'] = {'colors':np.array([0.8, 0.2, 0]),   'data':pm277}
            input_dict['f444w'] = {'colors':np.array([1.0, 0, 0]),     'data':pm444}
        
        elif field=='primer-cosmos':
            f090w_cutout_pm = gen_cutout(field, 'f090w', coord, cutout_width, tile=tile, suffix='psfMatched')
            f115w_cutout_pm = gen_cutout(field, 'f115w', coord, cutout_width, tile=tile, suffix='psfMatched')
            f150w_cutout_pm = gen_cutout(field, 'f150w', coord, cutout_width, tile=tile, suffix='psfMatched')
            f200w_cutout_pm = gen_cutout(field, 'f200w', coord, cutout_width, tile=tile, suffix='psfMatched')
            f277w_cutout_pm = gen_cutout(field, 'f277w', coord, cutout_width, tile=tile, suffix='psfMatched')
            f356w_cutout_pm = gen_cutout(field, 'f356w', coord, cutout_width, tile=tile, suffix='psfMatched')
            f410m_cutout_pm = gen_cutout(field, 'f410m', coord, cutout_width, tile=tile, suffix='psfMatched')
            f444w_cutout_pm = gen_cutout(field, 'f444w', coord, cutout_width, tile=tile)
            pm090 = f090w_cutout_pm.sci
            pm115 = f115w_cutout_pm.sci
            pm150 = f150w_cutout_pm.sci
            pm200 = f200w_cutout_pm.sci
            pm277 = f277w_cutout_pm.sci
            pm356 = f356w_cutout_pm.sci
            pm410 = f410m_cutout_pm.sci
            pm444 = f444w_cutout_pm.sci
            pm090[np.isnan(pm090)] = 0
            pm115[np.isnan(pm115)] = 0
            pm150[np.isnan(pm150)] = 0
            pm200[np.isnan(pm200)] = 0
            pm277[np.isnan(pm277)] = 0
            pm356[np.isnan(pm356)] = 0
            pm410[np.isnan(pm410)] = 0
            pm444[np.isnan(pm444)] = 0

            input_dict = {}
            input_dict['f090w'] = {'colors':np.array([0.2, 0.3, 1.0]), 'data':f090w_cutout_pm.sci}
            input_dict['f115w'] = {'colors':np.array([0.2, 0.3, 1.0]), 'data':f115w_cutout_pm.sci}
            input_dict['f150w'] = {'colors':np.array([0, 1.0, 0]),     'data':f150w_cutout_pm.sci}
            input_dict['f200w'] = {'colors':np.array([0, 1.0, 0]),     'data':f200w_cutout_pm.sci}
            input_dict['f277w'] = {'colors':np.array([0.95, 0.55, 0]), 'data':f277w_cutout_pm.sci}
            input_dict['f356w'] = {'colors':np.array([0.95, 0.55, 0]), 'data':f356w_cutout_pm.sci}
            input_dict['f410m'] = {'colors':np.array([1.0, 0, 0]),     'data':f410m_cutout_pm.sci}
            input_dict['f444w'] = {'colors':np.array([1.0, 0, 0]),     'data':f444w_cutout_pm.sci}

        imrgb = gen_rgb_image(input_dict, **rgb_params)
        ax_rgb1.imshow(imrgb, extent=f444w_cutout_pm.extent)
        ax_rgb1.tick_params(labelleft=False,labelbottom=False,left=False,right=False,bottom=False,top=False,which='both')
        ax_rgb1.spines['left'].set_visible(False)    
        ax_rgb1.spines['bottom'].set_visible(False)    
        ax_rgb1.spines['right'].set_visible(False)    
        ax_rgb1.spines['top'].set_visible(False)    
        ax_rgb1.set_xlabel('NIRCAM RGB')
        ax_rgb1.set_xlim(-display_width.to(u.arcsec).value/2, display_width.to(u.arcsec).value/2)
        ax_rgb1.set_ylim(-display_width.to(u.arcsec).value/2, display_width.to(u.arcsec).value/2)


        if plot_msa_slits:
            from astropy.table import Table
            assert msa_target_info is not None
            if type(msa_target_info)==str:
                msa_target_info = [msa_target_info]
            shutter_size_x = 0.26785
            open_size_x = 0.20
            shutter_size_y = 0.5232
            open_size_y = 0.46
            for msa_target_info_file in msa_target_info:
                tab = Table.read(msa_target_info_file)
                tab = tab[tab['ID']==ID]
                if len(tab)==0: continue
                pa = float(tab['Aperture PA (Degrees)'])
                for addl_y_offset in [-1,0,1]:
                    offset_x = float(tab['Offset (x)'])
                    offset_y = float(tab['Offset (y)'])
                    offset_x = offset_x*shutter_size_x
                    offset_y = offset_y*shutter_size_y + addl_y_offset*shutter_size_y
                    rot_offset_x1 =  np.cos(np.radians(pa))*offset_x - np.sin(np.radians(pa))*offset_y 
                    rot_offset_y1 =  np.sin(np.radians(pa))*offset_x + np.cos(np.radians(pa))*offset_y 
                    rot_offset_x2 =  np.cos(np.radians(pa))*(offset_x-0.035) - np.sin(np.radians(pa))*(offset_y-0.035)
                    rot_offset_y2 =  np.sin(np.radians(pa))*(offset_x-0.035) + np.cos(np.radians(pa))*(offset_y-0.035)

                    p = mpl.patches.Rectangle((rot_offset_x1, rot_offset_y1), -shutter_size_x, -shutter_size_y, angle=pa, rotation_point='xy', facecolor='none', edgecolor='w', linestyle='--', alpha=0.4, linewidth=0.5)
                    ax_rgb1.add_patch(p)
                    p = mpl.patches.Rectangle((rot_offset_x2, rot_offset_y2), -open_size_x, -open_size_y, angle=pa, rotation_point='xy', facecolor='none', edgecolor='w', alpha=0.4, linewidth=0.5)
                    ax_rgb1.add_patch(p)



        ax_rgb2.imshow(imrgb, extent=f444w_cutout_pm.extent)
        ax_rgb2.axis('off')
        ax_rgb2.set_xlim(-cutout_width.to(u.arcsec).value/2, cutout_width.to(u.arcsec).value/2)
        ax_rgb2.set_ylim(-cutout_width.to(u.arcsec).value/2, cutout_width.to(u.arcsec).value/2)
        d = display_width.to(u.arcsec).value/2
        ax_rgb2.plot([-d]*2,[-d,d],color='w',linewidth=0.5,linestyle=':',alpha=0.5)
        ax_rgb2.plot([d]*2,[-d,d],color='w',linewidth=0.5,linestyle=':',alpha=0.5)
        ax_rgb2.plot([-d,d],[-d]*2,color='w',linewidth=0.5,linestyle=':',alpha=0.5)
        ax_rgb2.plot([-d,d],[d]*2,color='w',linewidth=0.5,linestyle=':',alpha=0.5)

        if field in ['primer-cosmos','cosmos-web']: 
            # ancillary data only available in COSMOS, for now
            if verbose: print('\t Plotting ancillary data...')
            from reproject import reproject_interp

            mips = fits.open('/data/COSMOS/mosaics/mips_24_GO3_sci_10.fits')[0]
            mips_cutout = Cutout2D(mips.data, coord, size=cutout_width*5, wcs=WCS(mips.header,naxis=2))
            mips_err = sigma_clipped_stats(mips_cutout.data)[2]
            mips_data, _ = reproject_interp((mips_cutout.data/mips_err, mips_cutout.wcs), f444w_cutout_pm.wcs, shape_out=np.shape(f444w_cutout_pm))
            if np.max(mips_data)>7:
                levels = np.linspace(3,np.max(mips_data),3)
            else:
                levels = [3,5]
            p1 = ax_rgb2.contour(mips_data, extent=f444w_cutout_pm.extent, linewidths=0.8, colors='gold', levels=levels)

            vla = fits.open('/data/COSMOS/mosaics/vla_3ghz_msmf.fits')[0]
            vla_cutout = Cutout2D(vla.data[0][0], coord, size=cutout_width*5, wcs=WCS(vla.header,naxis=2))
            vla_err = sigma_clipped_stats(vla_cutout.data)[2]
            vla_data, _ = reproject_interp((vla_cutout.data/vla_err, vla_cutout.wcs), f444w_cutout_pm.wcs, shape_out=np.shape(f444w_cutout_pm))
            if np.max(vla_data)>7:
                levels = np.linspace(3,np.max(vla_data),3)
            else:
                levels = [3,5]
            p2 = ax_rgb2.contour(vla_data, extent=f444w_cutout_pm.extent, linewidths=0.8, colors='lime', levels=levels)


            scuba2 = fits.open('/data/COSMOS/mosaics/S2COSMOS_20180927_850_snr_mf_crop.fits')[0]
            scuba2_cutout = Cutout2D(scuba2.data[0], coord, size=cutout_width*5, wcs=WCS(scuba2.header,naxis=2))
            scuba2_data, _ = reproject_interp((scuba2_cutout.data, scuba2_cutout.wcs), f444w_cutout_pm.wcs, shape_out=np.shape(f444w_cutout_pm))
            if np.max(scuba2_data)>7:
                levels = np.linspace(3.5,np.max(scuba2_data),5)
            else:
                levels = [3.5,5]
            p3 = ax_rgb2.contour(scuba2_data, extent=f444w_cutout_pm.extent, linewidths=0.8, colors='w', levels=levels)
            ax_rgb2.contour(scuba2_data, extent=f444w_cutout_pm.extent, linewidths=0.8, colors='w', linestyles=':', levels=[3])

            try:
                exmora = fits.open('/data/COSMOS/mosaics/exmora_feb23.snr.fits')[0]
                exmora_cutout = Cutout2D(exmora.data, coord, size=cutout_width*5, wcs=WCS(exmora.header,naxis=2))
                exmora_data, _ = reproject_interp((exmora_cutout.data, exmora_cutout.wcs), f444w_cutout_pm.wcs, shape_out=np.shape(f444w_cutout_pm))
                p4 = ax_rgb2.contour(exmora_data, extent=f444w_cutout_pm.extent, linewidths=0.8, colors='violet', levels=[5,7,10])
                ax_rgb2.contour(exmora_data, extent=f444w_cutout_pm.extent, linewidths=0.8, colors='violet', linestyles=':', levels=[3,3.5,4])
            except NoOverlapError:
                pass

            p1 = plt.Line2D((0, 1), (0, 0), linewidth=0.8, color='gold')
            p2 = plt.Line2D((0, 1), (0, 0), linewidth=0.8, color='lime')
            p3 = plt.Line2D((0, 1), (0, 0), linewidth=0.8, color='w')
            p4 = plt.Line2D((0, 1), (0, 0), linewidth=0.8, color='violet')

            ax_rgb2.legend((p1,p2,p3,p4),('MIPS 24µm','VLA 3GHz','SCUBA-2 850µm','ALMA 2mm'),loc=('lower center'),bbox_to_anchor=(0.5,-0.2),ncol=2,facecolor='k',edgecolor='0.3',labelcolor='w', frameon=True,fontsize=7)

        cutout = gen_cutout(field, 'detec', coord,cutout_width/2, tile=tile)
        d = np.sqrt(cutout.data)
        cen = np.shape(d)[0]//2
        ax_detc.imshow(d, extent=cutout.extent, vmin=1, vmax=np.nanpercentile(d[cen-20:cen+20,cen-20:cen+20],99), cmap='Greys_r')
        ax_detc.tick_params(labelleft=False,labelbottom=False,left=False,right=False,bottom=False,top=False,which='both')
        ax_detc.spines['left'].set_visible(False)    
        ax_detc.spines['bottom'].set_visible(False)    
        ax_detc.spines['right'].set_visible(False)    
        ax_detc.spines['top'].set_visible(False)    
        ax_detc.set_xlabel(r'DETECTION')
        ax_detc.set_xlim(-display_width.to(u.arcsec).value/2, display_width.to(u.arcsec).value/2)
        ax_detc.set_ylim(-display_width.to(u.arcsec).value/2, display_width.to(u.arcsec).value/2)
        a, b, theta = cat['kron1_a'], cat['kron1_b'], np.degrees(cat['theta_image'])
        ax_detc.add_patch(mpl.patches.Ellipse((0,0), width=2*a, height=2*b, angle=theta, facecolor='none', edgecolor='salmon', linestyle='-', linewidth=0.5))
        ax_detc.add_patch(mpl.patches.Circle((0,0), radius=0.15, facecolor='none', edgecolor='w', linestyle='--', linewidth=0.5))

        cutout = gen_cutout(field, 'segm',coord,cutout_width, tile=tile)
        ax_segm.imshow(cutout.data, extent=cutout.extent, cmap=cutout.cmap)
        ax_segm.tick_params(labelleft=False,labelbottom=False,left=False,right=False,bottom=False,top=False,which='both')
        ax_segm.spines['left'].set_visible(False)    
        ax_segm.spines['bottom'].set_visible(False)    
        ax_segm.spines['right'].set_visible(False)    
        ax_segm.spines['top'].set_visible(False)    
        ax_segm.set_xlabel('SEGMENTATION')
        ax_segm.set_xlim(-display_width.to(u.arcsec).value/2, display_width.to(u.arcsec).value/2)
        ax_segm.set_ylim(-display_width.to(u.arcsec).value/2, display_width.to(u.arcsec).value/2)

        def plot_data(ax, wav, wav_min, wav_max, flux, flux_err, colors, zorder=10, annotate=True, label=None):
            colors = np.array(colors)
            wav = wav[flux_err>0]
            wav_min = wav_min[flux_err>0]
            wav_max = wav_max[flux_err>0]
            colors = colors[flux_err>0]
            flux = flux[flux_err>0]
            flux_err = flux_err[flux_err>0]
            for w, w1, w2, f, f_err, c in zip(wav, wav_min, wav_max, flux, flux_err, colors):
                if verbose: print(f'\t {w:.2f}, {f:.2f}, {f_err:.2f}, {f/f_err:.1f}')
                s = f/f_err
                if s > 1:
                    ax.errorbar(w, f, yerr=f_err, xerr=[[w-w1],[w2-w]], linewidth=0, marker='s', ms=6, 
                                mfc='none', mec=c, elinewidth=1, ecolor=c, capthick=1, capsize=2, zorder=zorder)
                    if annotate: ax.annotate(fr'${s:.1f}\sigma$', (w, 1.15*(f+f_err)), ha='center', va='bottom', color=c, fontsize=6, bbox=dict(facecolor='w', edgecolor='none', pad=0.01, alpha=0.7), zorder=1000)
                else:
                    ax.errorbar(w, 2*f_err, yerr=0.5*f_err, xerr=[[w-w1],[w2-w]], uplims=True, linewidth=0,
                                mfc='none', mec=c, elinewidth=1, ecolor=c, capthick=1, capsize=2, zorder=zorder)
                    if annotate: ax.annotate(fr'${s:.1f}\sigma$', (w, 1.15*2*f_err), ha='center', va='bottom', color=c, fontsize=6, bbox=dict(facecolor='w', edgecolor='none', pad=0.01, alpha=0.7), zorder=1000)
            if label is not None:
                ax.errorbar(100, 1, yerr=1, xerr=1, linewidth=0, marker='s', ms=6, 
                            mfc='none', mec=c, elinewidth=1, ecolor=c, capthick=1, capsize=2, zorder=zorder, label=label)


        if verbose: print('\t Plotting photometry...')
        in_cat = np.array([f'f_auto_{b}' in cat.dtype.names for b in config['bands'][field]])
        if np.any(~in_cat):
            print('Missing filters in catalog', np.array(config['bands'][field])[~in_cat])
        filters = np.array(config['filters'][field])[in_cat]
        bands = np.array(config['bands'][field])[in_cat]
        wav = get_lameff(filters)
        wav_min, wav_max = get_lamrange(filters)
        colors = np.array(colors)[in_cat]

        maxflux = 0
        if phot_columns[0] == 'auto':
            flux = np.array([cat[f'f_auto_{b}'][0] for b in bands])
            flux_err = np.array([cat[f'e_auto_{b}'][0] for b in bands])
            maxflux = np.nanmax([maxflux,np.nanmax(flux)])
            plot_data(ax_sed, wav, wav_min, wav_max, flux, flux_err, colors, label='AUTO')
        elif 'auto' in phot_columns:
            flux = np.array([cat[f'f_auto_{b}'][0] for b in bands])
            flux_err = np.array([cat[f'e_auto_{b}'][0] for b in bands])
            maxflux = np.nanmax([maxflux,np.nanmax(flux)])
            plot_data(ax_sed, wav, wav_min, wav_max, flux, flux_err, ['0.7']*len(flux), annotate=False, label='AUTO')
        if phot_columns[0].startswith('aper'):
            which_aper = int(phot_columns[0][-1])-1
            flux = np.array([cat[f'f_aper_{b}'][0,which_aper] for b in bands])
            flux_err = np.array([cat[f'e_aper_{b}'][0,which_aper] for b in bands])
            maxflux = np.nanmax([maxflux,np.nanmax(flux)])
            plot_data(ax_sed, wav, wav_min, wav_max, flux, flux_err, colors, label=f"{aper_diams[which_aper]}'' APER")
        elif ('aper1' in phot_columns or 'aper2' in phot_columns or 'aper3' in phot_columns or 'aper4' in phot_columns or 'aper5' in phot_columns):
            which_aper = int(phot_columns[1][-1])-1
            flux = np.array([cat[f'f_aper_{b}'][0,which_aper] for b in bands])
            flux_err = np.array([cat[f'e_aper_{b}'][0,which_aper] for b in bands])
            maxflux = np.nanmax([maxflux,np.nanmax(flux)])
            plot_data(ax_sed, wav, wav_min, wav_max, flux, flux_err, ['0.7']*len(flux), annotate=False, label=f"{aper_diams[which_aper]}'' APER")

        if field in ['cosmos-web','primer-cosmos']:
            if i_sepplus>=0:
                filters = ['subaru_hsc_g','subaru_hsc_r','subaru_hsc_i','subaru_hsc_z','subaru_hsc_y','uvista_Y', 'uvista_J', 'uvista_H', 'uvista_Ks','spitzer_irac_ch1','spitzer_irac_ch3','spitzer_irac_ch4','hst_acs_f814w','jwst_nircam_f115w','jwst_nircam_f150w','jwst_nircam_f277w','jwst_nircam_f444w','jwst_miri_f770w']
                wav = get_lameff(filters)
                wav_min, wav_max = get_lamrange(filters)
                flux = np.array([cat_sepp[f'FLUX_MODEL_{b}'][i_sepplus] for b in ['HSC-g','HSC-r','HSC-i','HSC-z','HSC-y', 'UVISTA-Y', 'UVISTA-J', 'UVISTA-H', 'UVISTA-Ks', 'IRAC-ch1', 'IRAC-ch2', 'IRAC-ch3', 'HST-F814W', 'F115W', 'F150W', 'F277W', 'F444W', 'F770W']])*1e29
                flux_err = np.array([cat_sepp[f'FLUX_ERR_MODEL_{b}'][i_sepplus] for b in ['HSC-g','HSC-r','HSC-i','HSC-z','HSC-y', 'UVISTA-Y', 'UVISTA-J', 'UVISTA-H', 'UVISTA-Ks', 'IRAC-ch1', 'IRAC-ch2', 'IRAC-ch3', 'HST-F814W', 'F115W', 'F150W', 'F277W', 'F444W', 'F770W']])*1e29
                # depths = np.array([28.1, 27.8, 27.6, 27.2, 26.5, 26.02, 25.85, 25.48, 25.14, ])        
                # flux_err = np.sqrt(flux_err**2 + (3631e6*np.power(10., -0.4*depths)/5)**2)
                c2 = np.array(['lightblue']*len(filters))
                ID_sepplus = cat_sepp['ID'][i_sepplus]
                plot_data(ax_sed, wav, wav_min, wav_max, flux, flux_err, c2, zorder=-1000, annotate=False, label=f'SE++ {ID_sepplus}')

        if field=='cosmos-web' and 'super' in input_catalog:
            in_primer = bool(cat['in_primer'][0])
            if in_primer:
                ID_primer = int(cat['id_primer'][0])
                filters = ['hst_acs_f606w','hst_acs_f814w','jwst_nircam_f090w','jwst_nircam_f115w',
                           'jwst_nircam_f150w','jwst_nircam_f200w','jwst_nircam_f277w','jwst_nircam_f356w',
                           'jwst_nircam_f410m','jwst_nircam_f444w']
                wav = get_lameff(filters)
                wav_min, wav_max = get_lamrange(filters)
                flux = np.array([cat[f'f_primer_auto_{b}'][0] for b in ['f606w','f814w','f090w','f115w','f150w','f200w','f277w','f356w','f410m','f444w']])
                flux_err = np.array([cat[f'e_primer_auto_{b}'][0] for b in ['f606w','f814w','f090w','f115w','f150w','f200w','f277w','f356w','f410m','f444w']])
                c2 = np.array(['darkorange']*len(filters))
                plot_data(ax_sed, wav, wav_min, wav_max, flux, flux_err, c2, zorder=-1000, annotate=False, label=f'PRIMER {ID_primer}')

<<<<<<< HEAD
        if eazy_run is not None:
            if type(eazy_run)==str:
                ezrs = [eazy_run]
            else:
                ezrs = eazy_run
            
            ez_c = ['slategrey','mediumslateblue']
            for k,ezr in enumerate(ezrs):
                if verbose: print('\t Plotting EAzY results...')
                try:
                    suffix = ezr.split('_')[-1]
                    twav, tflux, zbest = get_tem_sed(main=ezr, outdir=eazy_outdir, which=ID)
                    f = load_eazy_catalog(main=ezr, outdir=eazy_outdir, ext=1)
                    chi2 = f['chi2'][f['ID']==ID][0]
                    ax_sed.plot(twav, tflux, label=f'EAzY {suffix}, $z={zbest:.2f}$ ($\chi^2={chi2:.1f}$)', color=ez_c[k], linewidth=1, zorder=5)

=======
        if field=='cosmos-web' and 'supercatalog_pm' in input_catalog:
            has_pm = ~np.isnan(cat['f_primer_f1800w'][0])
            if has_pm:
                filters = ['jwst_miri_f770w','jwst_miri_f1800w']
                wav = get_lameff(filters)
                wav_min, wav_max = get_lamrange(filters)
                flux = np.array([cat[f'f_primer_{b}'][0] for b in ['f770w','f1800w']])
                flux_err = np.array([cat[f'e_primer_{b}'][0] for b in ['f770w','f1800w']])
                c2 = np.array(['darkorange']*len(filters))
                plot_data(ax_sed, wav, wav_min, wav_max, flux, flux_err, c2, zorder=-1000, annotate=False)

        if eazy_run is not None:
            if type(eazy_run)==str:
                ezrs = [eazy_run]
            else:
                ezrs = eazy_run
            
            ez_c = ['slategrey','mediumslateblue']
            for k,ezr in enumerate(ezrs):
                if verbose: print('\t Plotting EAzY results...')
                try:
                    suffix = ezr.split('_')[-1]
                    twav, tflux, zbest = get_tem_sed(main=ezr, outdir=eazy_outdir, which=ID)
                    f = load_eazy_catalog(main=ezr, outdir=eazy_outdir, ext=1)
                    chi2 = f['chi2'][f['ID']==ID][0]
                    ax_sed.plot(twav, tflux, label=f'EAzY {suffix}, $z={zbest:.2f}$ ($\chi^2={chi2:.1f}$)', color=ez_c[k], linewidth=1, zorder=5)

>>>>>>> 55572f3 (test)
                    ### P(z)
                    z, Pz = get_pz(main=ezr, outdir=eazy_outdir, which=ID)
                    Pz /= np.nanmax(Pz)
                    ax_pz.fill_between(z, Pz, color=ez_c[k], alpha=0.1, zorder=5)
                    ax_pz.plot(z, Pz, color=ez_c[k], linewidth=1, zorder=10)
                except:
                    if verbose: print(f'\t ! EAzY run {ezr} not found')
                    pass

        ################################################################################################################
        ################################################################################################################
        ################################################################################################################
            
        if show_bd:
            if verbose: print('\t Plotting brown dwarf results...')
            bd = fits.getdata('/data/COSMOS-Web/catalogs/COSMOS-Web_aperture_catalog_v0.7_brown_dwarf_results.fits')
            bd = bd[bd['id']==ID]
            family, chi2, temp, grav, zmet, dist = [bd[f][0] for f in ['model_family', 'chi2', 'temp', 'grav', 'zmet', 'dist']]
            lam, fnu = bd['lam'][0], bd['fnu'][0]
            ax_sed.plot(lam, fnu, color='0.7', linewidth=0.5, label=f'{family} $T={int(temp)}$ K, $D={int(dist)}$ pc ($\chi^2 = {chi2:.1f}$)')
            filters = ['hst_acs_f814w','jwst_nircam_f115w','jwst_nircam_f150w','jwst_nircam_f277w','jwst_nircam_f444w','jwst_miri_f770w']
            flux = [bd[f'flux_mod_{filt}'][0] for filt in filters]
            wav = get_lameff(filters)
            ax_sed.scatter(wav, flux, color='0.7')
        
        ################################################################################################################
        ################################################################################################################
        if bqf_run is not None:
            try:
                if verbose: print('\t Plotting QSO results...')
                hdr = fits.getheader(f'bqf/posterior/{bqf_run}/{ID}_bqf_results.fits', ext=1)
                chi2, zbest = hdr['chi2'], hdr['zbest']
                Muv = np.median(fits.getdata(f'bqf/posterior/{bqf_run}/{ID}_bqf_results.fits', ext=4)['Muv'])
                bqf = fits.getdata(f'bqf/posterior/{bqf_run}/{ID}_bqf_results.fits', ext=1)
                ax_sed.plot(bqf['wav_obs'], bqf['f_nu'], color='royalblue', linewidth=0.5, label=f'BQF $z={zbest:.2f}$ ($\chi^2={chi2:.1f}$' + r', $M_{\rm UV,int} = ' + f'{Muv:.1f}$)')

                bqf = fits.getdata(f'bqf/posterior/{bqf_run}/{ID}_bqf_results.fits', ext=2)
                ax_sed.scatter(bqf['wav_obs'], bqf['f_nu_mod'], color='royalblue')

                bqf = fits.getdata(f'bqf/posterior/{bqf_run}/{ID}_bqf_results.fits', ext=3)
                z, Pz = bqf['zgrid'], bqf['Pz']
                Pz /= np.max(Pz)
                ax_pz.fill_between(z, Pz, color='royalblue', alpha=0.1, zorder=15)
                ax_pz.plot(z, Pz, color='royalblue', linewidth=1, zorder=20)
            except:
                if verbose: print(f'file bqf/posterior/{bqf_run}/{ID}_bqf_results.fits not found!')

        if bagpipes_run is not None:
            if type(bagpipes_run)==str:
                bagpipes_run = [bagpipes_run]
            for bp_run in bagpipes_run:
                try:
                    if verbose: print('\t Plotting BAGPIPES results...')
                    hdr = fits.getheader(f'pipes/posterior/{bp_run}/{ID}_bagpipes_results.fits', ext=1)
                    chi2, zbest = hdr['chi2'], hdr['zbest']
                    massformed = np.median(fits.getdata(f'pipes/posterior/{bp_run}/{ID}_bagpipes_results.fits', ext=4)['massformed'])
                    pipes = fits.getdata(f'pipes/posterior/{bp_run}/{ID}_bagpipes_results.fits', ext=1)
                    ax_sed.plot(pipes['wav_obs'], pipes['f_nu'], color='purple', linewidth=0.5, label=f'BAGPIPES $z={zbest:.2f}$ ($\chi^2={chi2:.1f}$' + r', $\log M_{\star}/M_\odot = ' + f'{massformed:.1f}$)')

                    pipes = fits.getdata(f'pipes/posterior/{bp_run}/{ID}_bagpipes_results.fits', ext=2)
                    ax_sed.scatter(pipes['wav_obs'], pipes['f_nu_mod'], color='purple')

                    pipes = fits.getdata(f'pipes/posterior/{bp_run}/{ID}_bagpipes_results.fits', ext=3)
                    z, Pz = pipes['zgrid'], pipes['Pz']
                    Pz /= np.max(Pz)
                    ax_pz.fill_between(z, Pz, color='purple', alpha=0.1, zorder=15)
                    ax_pz.plot(z, Pz, color='purple', linewidth=1, zorder=20)
                except:
                    if verbose: print(f'file pipes/posterior/{bp_run}/{ID}_bagpipes_results.fits not found!')

        if plot_msa_wavelength_coverage:
            assert msa_target_info is not None
            from astropy.table import Table
            if type(msa_target_info)==str:
                msa_target_info = [msa_target_info]
            for msa_target_info_file in msa_target_info:
                tab = Table.read(msa_target_info_file)
                tab = tab[tab['ID']==ID]
                if len(tab)==0: continue

                prism_min_wav, prism_max_wav = 0.6, 5.3
                print(float(tab['NRS1 Min Wave']), float(tab['NRS1 Max Wave']), float(tab['NRS2 Min Wave']), float(tab['NRS2 Max Wave']))
                if tab['NRS1 Max Wave'] == -1: # spectrum only in NRS2
                    if tab['NRS2 Min Wave'] == -1:
                        min_wav = prism_min_wav
                    else:
                        min_wav = float(tab['NRS2 Min Wave'])
                    if tab['NRS2 Max Wave'] == -2:
                        max_wav = prism_max_wav
                    else:
                        max_wav = float(tab['NRS2 Max Wave'])
                    break_wav1, break_wav2 = min_wav, min_wav
                elif tab['NRS2 Min Wave'] == -2: # spectrum only in NRS1
                    if tab['NRS1 Min Wave'] == -1:
                        min_wav = prism_min_wav
                    else:
                        min_wav = float(tab['NRS1 Min Wave'])
                    if tab['NRS1 Max Wave'] == -2:
                        max_wav = prism_max_wav
                    else:
                        max_wav = float(tab['NRS2 Max Wave'])
                    break_wav1, break_wav2 = min_wav, min_wav
                else:
                    min_wav, max_wav = prism_min_wav, prism_max_wav
                    break_wav1, break_wav2 = float(tab['NRS1 Max Wave']), float(tab['NRS2 Min Wave'])

                ax_sed.fill_between([min_wav, break_wav1],[1e-10]*2,[1e10]*2,edgecolor='none', facecolor='k',alpha=0.04,zorder=-999)
                ax_sed.fill_between([break_wav2, max_wav],[1e-10]*2,[1e10]*2,edgecolor='none', facecolor='k',alpha=0.04,zorder=-999)




        
        ax_sed.legend(loc='lower right',  markerfirst=False, fontsize=7)#, title=legend_title)
        ax_sed.set_xlabel('Observed Wavelength [$\mu$m]')
        ax_sed.set_ylabel('Flux Density [$\mu$Jy]')
        ax_sed.set_xlim(0.25, 20)
        ax_sed.set_ylim(1e-3, 1e2)
        if maxflux > 1e2:
            ax_sed.set_ylim(1e-1, 1e3*maxflux)
        ax_sed.loglog()

        # ax_pz.set_ylabel('$P(z)$')
        ax_pz.set_xlabel('Redshift', size=7, labelpad=0.2)
        ax_pz.set_ylim(0, 1.05)
        ax_pz.set_xlim(0, 16)


        fig.text(0.03, 0.970, f'{field}', va='top', ha='left',fontsize=12)
        fig.text(0.03, 0.945, f'catalog: {catalog_name}', va='top', ha='left',fontsize=12)
        if nickname is None:
            fig.text(0.03, 0.920, f'ID: {ID}', va='top', ha='left',fontsize=12)
        else:
            fig.text(0.03, 0.920, f'ID: {ID}, aka: {nickname}', va='top', ha='left',fontsize=12)
        
        coordstring = coord.to_string('hmsdms', precision=2).split(' ')
        fig.text(0.97, 0.970, f'RA, Dec: ({coordstring[0]}, {coordstring[1]})', va='top', ha='right',fontsize=12)
        fig.text(0.97, 0.945, f'({coord.ra.value:.7f}, {coord.dec.value:.6f})', va='top', ha='right',fontsize=12)
        fig.text(0.97, 0.920, f'Tile: {tile}', va='top', ha='right',fontsize=12)
        
        snrs = [cat[f'snr_{b}'][0] for b in bands]
        for ax, snr, band, color in zip(cutout_axes, snrs, bands, colors):
            ax.annotate(band.upper(), (0.48, 1.03), xycoords='axes fraction', color=color, va='bottom', ha='right', fontsize=9)
            ax.annotate(fr'({snr:.1f}$\sigma$)', (0.52, 1.03), xycoords='axes fraction', color=color, va='bottom', ha='left', fontsize=9)
        
        ax_f444.annotate('$C = '+str(round(cat['C_f444w'][0],2))+'$', (-0.01, -0.01+0.02), color='0.9', xycoords='axes fraction', ha='left', va='bottom', fontsize=8, bbox=dict(boxstyle='round,pad=0.3', facecolor='0.9', edgecolor='k', linewidth=0.5))
        ax_f444.annotate('$C = '+str(round(cat['C_f444w'][0],2))+'$', (-0.01, -0.02+0.02), color='k',  xycoords='axes fraction', ha='left', va='bottom', fontsize=8)

        ax_pz.tick_params(labelleft=False,left=False,right=False,top=False,which='both')
        ax_pz.tick_params(direction='inout', which='both', labelsize=7, pad=1.5)
        ax_pz.spines['left'].set_visible(False)
        ax_pz.spines['right'].set_visible(False)
        ax_pz.spines['top'].set_visible(False)
        ax_sed.tick_params(right=False, which='both')
        ax_sed.set_xticks([0.3,0.4,0.6,1.0,1.5,2.0,3.0,4.0,5,7,10,15,20],['0.3','0.4','0.6','1','1.5','2','3','4','5','7','10','15','20'])

        def flux2mag(x):
            return -2.5*np.log10(x/3631e6)
        def mag2flux(x):
            return 3631e6*np.power(10., -0.4*x)

        secax1 = ax_sed.secondary_yaxis('right', functions=(flux2mag, mag2flux))
        secax1.set_ylabel('AB mag')
        secax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(2))
        secax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
        secax1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
        secax1.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
        
        plt.savefig(path, dpi=500)

        if show:
            plt.show()
        else:
            plt.close()

# # from astropy.coordinates import SkyCoord
# # def pc(i):
# #     c = SkyCoord(ra=f[f['id']==i]['RA_DETEC'][0], dec=f[f['id']==i]['DEC_DETEC'][0], unit='deg')
# #     print(c.to_string('hmsdms', precision=2))

def find_nearest(c, matching_radius=0.5, catalog='/data/COSMOS-Web/catalogs/COSMOS-Web_hot+cold_aperture_catalog_v1.2.fits'):
    f = fits.getdata(catalog)
    cs = SkyCoord(f['ra'], f['dec'], unit='deg')
    sep = cs.separation(c).to(u.arcsec).value
    i = np.argmin(sep)
    if sep[i]<matching_radius:
        print(f"Found match {f['id'][i]} located {sep[i]:.2f}'' away in tile {f['tile'][i]}")
        return f['id'][i]
    else:
        print(f"No match within {matching_radius:.2f}''")
        return None


    # cosmos_inspection_plot(ID, 
    #     outdir=f'inspection_plots/', 
    #     bqf_run='cweb_run_8', 
    #     bagpipes_run='cweb_run_18_burst', suffix='burst',
    #     plot_eazy=False, show=False, overwrite=True, verbose=False)
    