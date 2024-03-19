from rich.console import Console
from rich.theme import Theme
from rich.traceback import install

custom_theme = Theme(
    {
        "data": "dark_sea_green4",
        "info": "yellow4",
        "warning": "red",
        "error": "bold red",
        "repr.number": "bold bright_blue",
        "rule.line": "bright_yellow",
        "panel": "yellow4",
    }
)
console = Console(theme=custom_theme)
INFO = "[[info]INFO[/info]]"
WARNING = "[[warning]WARNING[/warning]]"
ERROR = "[[error]ERROR[/error]]"

install()


##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
from scipy.special import gammainc
import numpy as np 
from scipy.interpolate import interp1d

# def normal_asymm(loc=0, scale_up=1, scale_lo=1, size=1):
#     scale = np.random.choice([scale_up, scale_lo], size=1)
#     norm1 = np.random.normal(loc=loc, scale=scale_lo, size=size)
#     norm2 = np.random.normal(loc=loc, scale=scale_up, size=size)


def get_filter(filt, filterpath='/Users/hba423/Dropbox/research/COSMOS-Web/DREaM/bagpipes/bagpipes_filters/'):
    '''Returns filter transmission curve (w,t) = (wavelength in angstroms, transmission).'''
    f = np.loadtxt(filterpath+filt+'.dat')
    w = f[:,0]
    t = f[:,1]
    return w,t

def get_fnu(x, y, filt, filterpath='/Users/hba423/Dropbox/research/COSMOS-Web/DREaM/bagpipes/bagpipes_filters/'):
    '''Compute fnu in a particular filter by convolving the filter curve with an SED
    Assumed the filter file has wavelengths in angstroms, and the input SED has wavelengths in microns.'''
    s = interp1d(x, y, bounds_error=False, fill_value=0)
    if type(filt)==str:
        f = np.loadtxt(filterpath+filt+'.dat')
        lam = f[:,0]/1e4 # in micron
        nu = 299800/lam # in GHz
        fi = f[:,1]
        fnu = np.trapz(s(lam)*fi/nu, x=nu)/np.trapz(fi/nu, x=nu)

    else:
        fnu = np.zeros(len(filt))
        for i in range(len(filt)):
            f = np.loadtxt(filterpath+filt[i]+'.dat')
            lam = f[:,0]/1e4 # in micron
            nu = 299800/lam # in GHz
            fi = f[:,1]
            fnu[i] = np.trapz(s(lam)*fi/nu, x=nu)/np.trapz(fi/nu, x=nu)

    return fnu



def get_lameff(filt, filterpath='/Users/hba423/Dropbox/research/COSMOS-Web/DREaM/bagpipes/bagpipes_filters/'):
    if type(filt)==str:
        f = np.loadtxt(filterpath+filt+'.dat')
        fi = f[:,1]
        lam = f[:,0]/1e4 # in micron
        lam_mean = np.trapz(lam*fi, x=lam)/np.trapz(fi, x=lam)
    else:
        lam_mean = np.zeros(len(filt))
        for i in range(len(filt)):
            f = np.loadtxt(filterpath+filt[i]+'.dat')
            fi = f[:,1]
            lam = f[:,0]/1e4 # in micron
            lam_mean[i] = np.trapz(lam*fi, x=lam)/np.trapz(fi, x=lam)

    return lam_mean


def get_lamrange(filt, filterpath='/Users/hba423/Dropbox/research/COSMOS-Web/DREaM/bagpipes/bagpipes_filters/'):
    if type(filt)==str:
        f = np.loadtxt(filterpath+filt+'.dat')
        wi = f[:,0]/1e4
        fi = f[:,1]/np.max(f[:,1])
        w_min = np.min(wi[fi>0.5])
        w_max = np.max(wi[fi>0.5])
    else:
        w_min, w_max = np.zeros(len(filt)),np.zeros(len(filt))
        for i in range(len(filt)):
            f = np.loadtxt(filterpath+filt[i]+'.dat')
            wi = f[:,0]/1e4
            fi = f[:,1]/np.max(f[:,1])
            w_min[i] = np.min(wi[fi>0.5])
            w_max[i] = np.max(wi[fi>0.5])

    return w_min, w_max
    

def redChiSq(ydata,ymod,std,deg):
    z = (ydata-ymod)/std
    chisq = np.sum(z**2)  
    nu = len(ydata)-deg  
    prob = 1-gammainc(0.5*nu, 0.5*chisq)
    return chisq/nu, prob



def get_tile(coord):
    from astropy.io import fits
    import warnings
    warnings.simplefilter('ignore')
    from astropy.wcs import WCS
    wcs = WCS(fits.open('/Users/hba423/data/COSMOS-Web/mosaics/Apr23/mask_f115w_tile_A10_60mas.fits')[0].header)
    from regions import Regions
    tiles = np.array(['A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','B1','B2','B3','B4','B5','B6','B7','B8','B9','B10'])
    in_tile = np.array(np.zeros(len(tiles)),dtype=bool)
    for i in range(len(tiles)):
        t = Regions.read(f'/Users/hba423/data/COSMOS-Web/regionfiles/tile_{tiles[i]}.reg', format='ds9')
        if t[0].contains(coord, wcs):
            in_tile[i] = True
        
    if not any(in_tile):
        raise Exception
    else:
        return tiles[in_tile][0]




from astropy.io import fits
from astropy.wcs import WCS

def get_filepath(field, tile, band, typ, suffix=''):
    if field in ['cosmos-web','primer-cosmos']:
        if band in ['Y','J','H','Ks']:
            if typ!='sci': 
                raise Exception('we dont have ERR or WHT maps for UVISTA! only SCI')
            return f'/data/COSMOS/mosaics/UVISTA_{band}_12_07_22_allpaw_skysub_015_dr5_rc_v1.fits'
        if band in ['g','r','i','z','y']:
            if typ!='sci':
                raise Exception('we dont have ERR or WHT maps for HSC! only SCI')
            return f'/data/COSMOS/mosaics/cutout-HSC-{band.upper()}-9813-pdr3_dud_rev_sci_zp-28.09_{tile}.fits'

    if band=='f814w':
        if typ.startswith('sci'): typ = typ.replace('sci','drz')
    if len(suffix)>0:
        typ += f'_{suffix}'

    if field=='cosmos-web':
        if tile.startswith('B'):
            if band=='f814w':
                return f'/V/simmons/data/cosmos-web/mosaics/mosaic_cosmos_web_2024jan_30mas_tile_{tile}_hst_acs_wfc_f814w_{typ}.fits'
            elif band=='f770w':
                return f'/V/simmons/data/cosmos-web/mosaics/mosaic_miri_f770w_COSMOS-Web_30mas_{tile}_{typ}.fits'
            else:
                return f'/V/simmons/data/cosmos-web/mosaics/mosaic_nircam_{band}_COSMOS-Web_30mas_{tile}_{typ}.fits'
        elif tile.startswith('A'):
            if band=='f814w':
                return f'/V/simmons/data/cosmos-web/mosaics/mosaic_cosmos_web_2023apr_30mas_tile_{tile}_hst_acs_wfc_f814w_{typ}.fits'
            elif band=='f770w':
                return f'/V/simmons/data/cosmos-web/mosaics/mosaic_miri_f770w_COSMOS-Web_30mas_{tile}_{typ}.fits'
            else:
                return f'/V/simmons/data/cosmos-web/mosaics/mosaic_nircam_{band}_COSMOS-Web_30mas_{tile}_{typ}.fits'
        else:
            print(f'tile "{tile}" not understood')
            return
                
    
    elif field=='primer-cosmos':
        if band in ['f606w','f814w']:
            if typ.startswith('sci'): typ = typ.replace('sci','drz')
            return f'/V/simmons/data/primer/mosaics/mosaic_cosmos_primer_30mas_acs_wfc_{band}_{typ}.fits'
        elif band in ['f125w','f160w']:
            if typ.startswith('sci'): typ = typ.replace('sci','drz')
            return f'/V/simmons/data/primer/mosaics/mosaic_cosmos_primer_30mas_wfc3_ir_{band}_{typ}.fits'
        elif band in ['f770w','f1800w']:
            return f'/V/simmons/data/primer/mosaics/mosaic_miri_{band}_PRIMER-COSMOS_epoch1_60mas_v0_2_{typ}.fits'
        else:
            return f'/V/simmons/data/primer/mosaics/mosaic_nircam_{band}_PRIMER-COSMOS_30mas_{typ}.fits'
    else:
        print(f'field "{field}" not understood')
        return
            


def load_image(field, tile, band, typ, suffix='', convert_units=False):
    '''
        'tile' == 'A1'--'A10' or 'B1'--'B10' or 'primer-cosmos'
        'band' == one of 'f814w', 'f115w', 'f150w', 'f277w', 'f444w', 'f770w'
        'typ' == one of 'sci', 'err', 'wht'
    '''
    filepath = get_filepath(field, tile, band, typ, suffix=suffix)  
    f = fits.open(filepath)
    if np.ndim(f[0].data)==0:
        data, header = f[1].data, f[1].header
        if 'NSCI_MU' in f[0].header:
            f[1].header['NSCI_MU'] = f[0].header['NSCI_MU']
        if 'NSCI_SIG' in f[0].header:
            f[1].header['NSCI_SIG'] = f[0].header['NSCI_SIG']
    else:
        data, header = f[0].data, f[0].header
    del f


    if convert_units:
        if typ=='sci' or typ=='err' or typ=='drz': # convert from e/s to uJy
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

    return data, header

def load_hdr(field, tile, band, typ, suffix=''):
    return load_image(field, tile, band, typ, suffix=suffix)[1]
def load_sci(field, tile, band, suffix='', convert_units=False):
    return load_image(field, tile, band, 'sci', suffix=suffix, convert_units=convert_units)[0]
def load_wcs(field, tile, band, suffix=''):
    return WCS(load_hdr(field, tile, band, 'sci', suffix=suffix))
def load_err(field, tile, band, suffix='', convert_units=False):
    return load_image(field, tile, band, 'err', suffix=suffix, convert_units=convert_units)[0]
def load_rms(field, tile, band, suffix='', convert_units=False):
    return load_image(field, tile, band, 'rms', suffix=suffix, convert_units=convert_units)[0]
def load_wht(field, tile, band, suffix=''):
    return load_image(field, tile, band, 'wht', suffix=suffix)[0]
def load_vrp(field, tile, band, suffix=''):
    return load_image(field, tile, band, 'vrp', suffix=suffix)[0]
def load_exp(field, tile, band, suffix=''):
    d, h = load_image(field, tile, band, 'exp', suffix=suffix)
    return d, WCS(h)


def gen_cutout(field, tile, band, coord, width, suffix='', err_method='simple', verbose=False):
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

    sci = load_sci(field, tile, band, suffix=suffix)
    wcs = load_wcs(field, tile, band, suffix=suffix)
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
    if field=='primer-cosmos':
        sci = fits.getdata(f'/V/simmons/data/primer/mosaics/detection_chi2pos_SWLW_primer-cosmos.fits')
        wcs = WCS(fits.getheader(f'/V/simmons/data/primer/mosaics/detection_chi2pos_SWLW_primer-cosmos.fits'))
    elif field=='cosmos-web':
        sci = fits.getdata(f'/V/simmons/data/cosmos-web/mosaics/detection_chi2pos_SWLW_{tile}.fits')
        wcs = WCS(fits.getheader(f'/V/simmons/data/cosmos-web/mosaics/detection_chi2pos_SWLW_{tile}.fits'))
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
    if field=='primer-cosmos':
        segm = np.load(f'/V/simmons/data/primer/mosaics/detection_chi2pos_SWLW_primer-cosmos_{tile}_segmap.npy')
        detec = fits.getdata(f'/V/simmons/data/primer/mosaics/detection_chi2pos_SWLW_primer-cosmos.fits')
        wcs = WCS(fits.getheader(f'/V/simmons/data/primer/mosaics/detection_chi2pos_SWLW_primer-cosmos.fits'))
        if tile == 'south':
            c0 = wcs.pixel_to_world(12480, 12500)
            cutout = Cutout2D(detec, c0, size=[11.0*u.arcmin,9.8*u.arcmin], wcs=wcs)
        elif tile == 'north': 
            c0 = wcs.pixel_to_world(12480, 33000)
            cutout = Cutout2D(detec, c0, size=[11.0*u.arcmin,9.8*u.arcmin], wcs=wcs)
        wcs = cutout.wcs
    elif field=='cosmos-web':
        segm = np.load(f'/V/simmons/data/cosmos-web/mosaics/detection_chi2pos_SWLW_{tile}_segmap.npy')
        wcs = WCS(fits.getheader(f'/V/simmons/data/cosmos-web/mosaics/detection_chi2pos_SWLW_{tile}.fits'))

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
