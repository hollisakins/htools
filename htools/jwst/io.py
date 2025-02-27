
# from htools.jwst.io import load_image

# from .. import host
from . import config
from . import paths

from astropy.io import fits
from astropy.wcs import WCS

def load_image(field, band, ext, tile=None, convert_units=False, ps='30mas'):
    '''
    Load an image from a given JWST field.

    Args: 
        'field': JWST field, one of config.all_jwst_fields
        'band': JWST band, one of config.all_jwst_bands or config.all_hst_bands
        'tile': 
        'ext': 'sci', 'err', or 'wht', or something like that
        'ps': pixel scale, '30mas' or '60mas'
        'convert_units': If True, convert the units of the image to microJy/pixel.
            Uses the header information to do this. 
            
    '''

    assert field in config.all_fields, f"Field {field} not recognized. Options are: {config.all_fields}"
    assert band in config.all_hst_bands or band in config.all_jwst_bands or band in config.misc_bands, f"Band {band} not recognized. Options are: {config.all_hst_bands + config.all_jwst_bands + config.misc_bands}"

    if field == 'cosmos-web':
        assert tile is not None, "For COSMOS-Web, you must specify a tile"
        assert band in config.bands['cosmos-web'] or band in config.bands['cosmos'], f"Band {band} not available for field {field}. Options are: {config.bands['cosmos-web'] + config.bands['cosmos']}"
        if band in config.bands['cosmos']:
            field = 'cosmos'
    else:
        assert band in config.bands[field], f"Band {band} not available for field {field}. Options are: {config.bands[field]}"

    filepath, hdu_index = paths.get_filepath(field, band, ext, tile=tile, ps=ps)  
    f = fits.open(filepath)
    
    data, header = f[hdu_index].data, f[hdu_index].header
    del f

    if convert_units:
        # import numpy as np
        # if 'sci' in ext or 'err' in ext:
        #     if band in config.all_hst_bands: 
        #         photflam, photplam = header['PHOTFLAM'], header['PHOTPLAM']
        #         conversion = 3.33564e10 * (photplam)**2 * photflam
            
        #     elif band in config.all_jwst_bands: # JWST
        #         if header['BUNIT'] == 'MJy/sr':
        #             wcs = WCS(header)
        #             pa = wcs.proj_plane_pixel_area().to('sr').value
        #             conversion = 1e12 * pa
        #         conversion = 1e12 * np.pi**2/(180**2) * (1/(3600**2)) * 0.03**2

        #     else:
        #         raise Exception(f"Not sure how to convert units for {band}.")
        #     data *= conversion
        # else:
        #     raise Exception(f"Not converting units for extension {ext}")
        raise Exception(f"Not yet configured to convert units.")

    return data, header


def load_sci(field, band, tile=None, ext='sci', convert_units=False, ps='30mas'):
    return load_image(field, band, ext, tile=tile, ps=ps, convert_units=convert_units)[0]

def load_err(field, band, tile=None, ext='err', convert_units=False, ps='30mas'):
    return load_image(field, band, ext, tile=tile, ps=ps, convert_units=convert_units)[0]

def load_wht(field, band, tile=None, ext='wht', ps='30mas'):
    return load_image(field, band, ext, tile=tile, ps=ps, convert_units=False)[0]

def load_hdr(field, band, tile=None, ext='sci', ps='30mas'):
    return load_image(field, band, ext, tile=tile, ps=ps)[1]

def load_wcs(field, band, tile=None, ext='sci', ps='30mas'):
    return WCS(load_hdr(field, band, ext=ext, tile=tile, ps=ps))

def load_cosmos_web_detec(tile):
    filepath, hdu_index = paths.get_cosmos_web_detec_filepath(tile)
    return fits.open(filepath)[hdu_index]

def load_cosmos_web_segm(tile, catalog_version='v1.3'):
    filepath, hdu_index = paths.get_cosmos_web_segm_filepath(tile, catalog_version=catalog_version)
    return fits.open(filepath)[hdu_index]