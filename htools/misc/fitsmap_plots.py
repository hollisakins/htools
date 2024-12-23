import numpy as np
import os, sys, tqdm, pickle
from PIL import Image
from rich.console import Console
console = Console()
log = console.log

from .. import host
from ..plotting import set_style
from ..utils.filters import Filters
from ..jwst import paths, io

import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.nddata.utils import NoOverlapError
from astropy.stats import sigma_clipped_stats

# plotting library
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib import transforms
import matplotlib as mpl
set_style('sans')

import warnings
from astropy.wcs import FITSFixedWarning
warnings.simplefilter('ignore', FITSFixedWarning)
warnings.simplefilter('ignore')




err_table = {
        'f814w':            {'pixel_scale':0.030, 'std1':76.888,  'alpha':0.346, 'beta':1.582},
        'f115w':            {'pixel_scale':0.030, 'std1':17.939,  'alpha':0.753, 'beta':1.296},
        'f150w':            {'pixel_scale':0.030, 'std1':19.057,  'alpha':0.643, 'beta':1.342},
        'f277w':            {'pixel_scale':0.030, 'std1':17.296,  'alpha':0.918, 'beta':1.399},
        'f444w':            {'pixel_scale':0.030, 'std1':24.260,  'alpha':1.077, 'beta':1.312},
        'f770w':            {'pixel_scale':0.030, 'std1':61.707,  'alpha':2.107, 'beta':1.264},
        'F814W_psfMatched': {'pixel_scale':0.030, 'std1':29.784,  'alpha':0.668, 'beta':1.647},
        'F115W_psfMatched': {'pixel_scale':0.030, 'std1':4.7430,  'alpha':2.230, 'beta':1.347},
        'F150W_psfMatched': {'pixel_scale':0.030, 'std1':4.8720,  'alpha':1.979, 'beta':1.393},
        'F277W_psfMatched': {'pixel_scale':0.030, 'std1':8.5610,  'alpha':1.523, 'beta':1.443},
        'CFHT-u':           {'pixel_scale':0.150, 'std1':0.579,   'alpha':0.453, 'beta':1.689},
        'g':            {'pixel_scale':0.167, 'std1':158.455, 'alpha':0.350, 'beta':1.731},
        'r':            {'pixel_scale':0.167, 'std1':160.421, 'alpha':0.368, 'beta':1.698},
        'i':            {'pixel_scale':0.167, 'std1':159.004, 'alpha':0.331, 'beta':1.780},
        'z':            {'pixel_scale':0.167, 'std1':158.813, 'alpha':0.439, 'beta':1.587},
        'y':            {'pixel_scale':0.167, 'std1':159.362, 'alpha':0.515, 'beta':1.483},
        'HSC-NB0816':       {'pixel_scale':0.167, 'std1':158.962, 'alpha':0.643, 'beta':1.322},
        'HSC-NB0921':       {'pixel_scale':0.167, 'std1':158.512, 'alpha':0.627, 'beta':1.340},
        'HSC-NB1010':       {'pixel_scale':0.167, 'std1':157.914, 'alpha':0.770, 'beta':1.185},
        'Y':         {'pixel_scale':0.150, 'std1':0.652,   'alpha':0.761, 'beta':1.843},
        'J':         {'pixel_scale':0.150, 'std1':0.604,   'alpha':0.779, 'beta':1.744},
        'H':         {'pixel_scale':0.150, 'std1':0.579,   'alpha':0.989, 'beta':1.587},
        'Ks':        {'pixel_scale':0.150, 'std1':0.558,   'alpha':1.213, 'beta':1.434},
        'UVISTA-NB118':     {'pixel_scale':0.150, 'std1':0.533,   'alpha':1.537, 'beta':1.301},
        'SC-IA484':         {'pixel_scale':0.150, 'std1':0.019,   'alpha':1.703, 'beta':1.277},
        'SC-IA527':         {'pixel_scale':0.150, 'std1':0.018,   'alpha':1.807, 'beta':1.303},
        'SC-IA624':         {'pixel_scale':0.150, 'std1':0.015,   'alpha':1.659, 'beta':1.358},
        'SC-IA679':         {'pixel_scale':0.150, 'std1':0.086,   'alpha':1.935, 'beta':1.000}, # fixed to beta=1
        'SC-IA738':         {'pixel_scale':0.150, 'std1':0.016,   'alpha':1.806, 'beta':1.308},
        'SC-IA767':         {'pixel_scale':0.150, 'std1':0.053,   'alpha':1.887, 'beta':1.000}, # fixed to beta=1
        'SC-IB427':         {'pixel_scale':0.150, 'std1':0.056,   'alpha':1.948, 'beta':1.000}, # fixed to beta=1
        'SC-IB505':         {'pixel_scale':0.150, 'std1':0.030,   'alpha':1.559, 'beta':1.147},
        'SC-IB574':         {'pixel_scale':0.150, 'std1':0.079,   'alpha':2.141, 'beta':1.000}, # fixed to beta=1
        'SC-IB709':         {'pixel_scale':0.150, 'std1':0.050,   'alpha':1.955, 'beta':1.000}, # fixed to beta=1
        'SC-IB827':         {'pixel_scale':0.150, 'std1':0.072,   'alpha':2.008, 'beta':1.000}, # fixed to beta=1
        'SC-NB711':         {'pixel_scale':0.150, 'std1':0.009,   'alpha':2.028, 'beta':1.419},
        'SC-NB816':         {'pixel_scale':0.150, 'std1':0.081,   'alpha':1.798, 'beta':1.000}, # fixed to beta=1
        'IRAC1':            {'pixel_scale':0.150, 'std1':581.215, 'alpha':0.989, 'beta':1.935},
        'IRAC2':            {'pixel_scale':0.150, 'std1':326.814, 'alpha':1.152, 'beta':1.818},
        'IRAC3':            {'pixel_scale':0.150, 'std1':209.957, 'alpha':2.273, 'beta':1.392},
        'IRAC4':            {'pixel_scale':0.150, 'std1':195.262, 'alpha':2.480, 'beta':1.341}, 
        }

filters = Filters(['f814w','f115w','f150w','f277w','f444w','f770w',
                    'cfht_u', 'hsc_g', 'hsc_r', 'hsc_i', 'hsc_z', 'hsc_y', 
                    'uvista_Y', 'uvista_J', 'uvista_H', 'uvista_Ks',
                    'IB427','IB484','IB505','IB527','IB574','IB624',
                    'IB679','IB709','IB738','IB767','IB827','NB711',
                    'NB816','HSC-NB816','HSC-NB921','HSC-NB1010','NB118',
                    'irac_ch1', 'irac_ch3', 'irac_ch4'])
wav, wav_min, wav_max = filters.wav.to(u.micron).value, filters.wav_min.to(u.micron).value, filters.wav_max.to(u.micron).value
# hues = np.interp(wav, [0.3, 4], [0.8, 0.05], left=0.8, right=0.05)
hues = np.interp(np.log10(wav), [-0.5, 0.6], [0.8, 0.05], left=0.8, right=0.05)
saturations = [1]*16 + [0.3]*20
values = [0.6]*16 + [0.9]*20
colors = []
for i in range(len(hues)):
    colors.append(mpl.colors.hsv_to_rgb([hues[i], saturations[i], values[i]]))
short_names = ['f814w','f115w','f150w','f277w','f444w','f770w','u','g','r','i','z','y','Y','J','H','Ks','IB427','IB484','IB505','IB527','IB574','IB624',
                    'IB679','IB709','IB738','IB767','IB827','NB711','NB816','HSC-NB816','HSC-NB921','HSC-NB1010','NB118','IRAC1','IRAC3','IRAC4']
colors_dict = dict(zip(short_names, colors))

def main(IDs, 
         logo = None,
         catalog_path = None,
         catalog_filename = None,
         catalog_shortname = None,
         outdir='inspection_plots/',
         out_format='pdf',
         dpi=None,
         overwrite=True,
         display_width = 4*u.arcsec,
         vmax_min = 10,
         cmap = 'Greys',
         lephare_spec_path=None,
         verbose=False,
    ):

    tiles = list(IDs.keys())

    catalog = fits.getdata(os.path.join(catalog_path, catalog_filename))

    log('Loading full RGB image...')
    rgb_path = paths.get_cosmos_web_rgb_filepath()
    with open(rgb_path, 'rb') as f:
        rgb_full = pickle.load(f)
    log('Done')


    for tile in tiles:
        log(f'Processing tile {tile}...')
        
        sci_dict, hdr_dict, wht_dict = {}, {}, {}
        for band in ['f814w','f115w','f150w','f277w','f444w','f770w','g','r','i','z','y','Y','J','H','Ks','IRAC1','IRAC3','IRAC4']:
            log(f'Reading in band {band}')
            sci_dict[band] = io.load_sci('cosmos-web', band, tile=tile) 
            hdr_dict[band] = io.load_hdr('cosmos-web', band, tile=tile) 
            wht_dict[band] = io.load_wht('cosmos-web', band, tile=tile) 
            # sci_filepath, hdu_index = get_filepath('cosmos-web', band, 'sci', tile=tile) # fits.open(sci_filepath)[hdu_index]
            # wht_filepath, hdu_index = get_filepath('cosmos-web', band, 'wht', tile=tile) #fits.open(wht_filepath)[hdu_index]

        model_path, hdu_index = paths.get_filepath('cosmos-web', 'f277w', 'mod', tile=tile, ps='60mas')
        model = fits.open(model_path)[hdu_index]
        detec = io.load_cosmos_web_detec(tile)
        segm = io.load_cosmos_web_segm(tile, catalog_version='v1.3')

        for i in (pbar := tqdm.tqdm(range(len(IDs[tile])))):
            ID = IDs[tile][i]
            pbar.set_description(str(ID))

            outfilename = f"cosmos-web_sed_{catalog_shortname.replace('-','_')}_{ID}"
            outpath = os.path.join(outdir, outfilename)
            if out_format == 'png':
                outpath += '.png'
            elif out_format == 'pdf':
                outpath += '.pdf'

            if not overwrite: 
                if os.path.exists(outpath):
                    if verbose: print(f'\t Skipping, plot already exists at {outpath}.')
                    continue

            ######################################################################################################################################
            cat = catalog[catalog['ID_SE++']==ID]
            ra = cat['RA_MODEL'][0]
            dec = cat['DEC_MODEL'][0]
            coord = SkyCoord(ra=ra, dec=dec, unit=u.deg)
            
            Reff = cat['RADIUS'][0] * 3600
            Reff_err = cat['RADIUS_err'][0] * 3600
            n_sersic = cat['SERSIC'][0]
            n_sersic_err = cat['SERSIC_err'][0]

            display_width_options = [4, 5, 6, 8, 10, 15]
            scalebar_size_options = [2, 3, 3, 4, 4, 5]
            i = np.argmin(np.abs(np.array(display_width_options)-display_width.to(u.arcsec).value))
            display_width = display_width_options[i]*u.arcsec
            scalebar_size = scalebar_size_options[i]*u.arcsec
            i = 0
            while Reff*10 > display_width.to(u.arcsec).value:
                try:
                    display_width = display_width_options[i]*u.arcsec
                    scalebar_size = scalebar_size_options[i]*u.arcsec
                    i += 1
                except:
                    display_width = display_width_options[-1]*u.arcsec
                    scalebar_size = scalebar_size_options[-1]*u.arcsec
                    break
            cutout_width = display_width*1.5



            ######################################################################################################################################
            fig = plt.figure(figsize=(14.4, 8.2), constrained_layout=False)
            gs = mpl.gridspec.GridSpec(ncols=54, nrows=30, width_ratios=[1.01]*9 + [1]*36 + [1.01]*9, figure=fig)
            gs.update(hspace=0.08, wspace=0.15, left=0.015, bottom=0.025, right=1-0.015, top=1-0.025)

            if logo is not None:
                logo_img = Image.open(logo)
                logo_img = logo_img.transpose(Image.FLIP_TOP_BOTTOM)

                # Create an inset axes
                axins = fig.add_axes([0.843, 0.883, 0.142, 0.11])
                axins.imshow(logo_img, alpha=1, aspect='auto')
                # axins.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labeltop=False, labelleft=False, labelright=False)
                axins.axis('off')


            fig.text(0.015, 0.980, f'ID ({catalog_shortname}): {ID}', va='top', ha='left',fontsize=12, weight='bold')
            coordstring = coord.to_string('hmsdms', precision=2).split(' ')
            fig.text(0.015, 0.955, f'Tile: {tile}', va='top', ha='left',fontsize=12)
            fig.text(0.015, 0.930, f'RA, Dec: ({coordstring[0]}, {coordstring[1]})', va='top', ha='left',fontsize=12)
            fig.text(0.069 , 0.905, f'({coord.ra.value:.7f}, {coord.dec.value:.6f})', va='top', ha='left',fontsize=12)


            ax_rgb = plt.subplot(gs[3:12,0:9])
            ax_detc = plt.subplot(gs[12:21,0:9])
            ax_segm = plt.subplot(gs[21:30,0:9])
            ax_mod = plt.subplot(gs[3:12,-9:])

            ax_f814 =  plt.subplot(gs[3:9,9:15])
            ax_f115 =  plt.subplot(gs[3:9,15:21])
            ax_f150 =  plt.subplot(gs[3:9,21:27])
            ax_f277 =  plt.subplot(gs[3:9,27:33])
            ax_f444 =  plt.subplot(gs[3:9,33:39])
            ax_f770 =  plt.subplot(gs[3:9,39:45])

            ax_g =     plt.subplot(gs[9:12,9:12])
            ax_r =     plt.subplot(gs[9:12,12:15])
            ax_i =     plt.subplot(gs[9:12,15:18])
            ax_z =     plt.subplot(gs[9:12,18:21])
            ax_y =     plt.subplot(gs[9:12,21:24])
            ax_Y =     plt.subplot(gs[9:12,24:27])
            ax_J =     plt.subplot(gs[9:12,27:30])
            ax_H =     plt.subplot(gs[9:12,30:33])
            ax_K =     plt.subplot(gs[9:12,33:36])
            ax_irac1 = plt.subplot(gs[9:12,36:39])
            ax_irac3 = plt.subplot(gs[9:12,39:42])
            ax_irac4 = plt.subplot(gs[9:12,42:45])

            ax_f814.text(0.05, 0.95, 'F814W', transform=ax_f814.transAxes, color='k', va='top', ha='left', path_effects=[pe.withStroke(linewidth=2, foreground='0.9')], weight='bold', size=12)
            ax_f115.text(0.05, 0.95, 'F115W', transform=ax_f115.transAxes, color='k', va='top', ha='left', path_effects=[pe.withStroke(linewidth=2, foreground='0.9')], weight='bold', size=12)
            ax_f150.text(0.05, 0.95, 'F150W', transform=ax_f150.transAxes, color='k', va='top', ha='left', path_effects=[pe.withStroke(linewidth=2, foreground='0.9')], weight='bold', size=12)
            ax_f277.text(0.05, 0.95, 'F277W', transform=ax_f277.transAxes, color='k', va='top', ha='left', path_effects=[pe.withStroke(linewidth=2, foreground='0.9')], weight='bold', size=12)
            ax_f444.text(0.05, 0.95, 'F444W', transform=ax_f444.transAxes, color='k', va='top', ha='left', path_effects=[pe.withStroke(linewidth=2, foreground='0.9')], weight='bold', size=12)
            ax_f770.text(0.05, 0.95, 'F770W', transform=ax_f770.transAxes, color='k', va='top', ha='left', path_effects=[pe.withStroke(linewidth=2, foreground='0.9')], weight='bold', size=12)

            ax_g.text(0.05, 0.95, 'HSC $g$',    transform=ax_g.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=1.3, foreground='k')], size=8)
            ax_r.text(0.05, 0.95, 'HSC $r$',    transform=ax_r.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=1.3, foreground='k')], size=8)
            ax_i.text(0.05, 0.95, 'HSC $i$',    transform=ax_i.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=1.3, foreground='k')], size=8)
            ax_z.text(0.05, 0.95, 'HSC $z$',    transform=ax_z.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=1.3, foreground='k')], size=8)
            ax_y.text(0.05, 0.95, 'HSC $y$',    transform=ax_y.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=1.3, foreground='k')], size=8)
            ax_Y.text(0.05, 0.95, 'UVISTA $Y$', transform=ax_Y.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=1.3, foreground='k')], size=8)
            ax_J.text(0.05, 0.95, 'UVISTA $J$', transform=ax_J.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=1.3, foreground='k')], size=8)
            ax_H.text(0.05, 0.95, 'UVISTA $H$', transform=ax_H.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=1.3, foreground='k')], size=8)
            ax_K.text(0.05, 0.95, 'UVISTA $K$', transform=ax_K.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=1.3, foreground='k')], size=8)
            ax_irac1.text(0.05, 0.95, 'IRAC1',  transform=ax_irac1.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=1.3, foreground='k')], size=8)
            ax_irac3.text(0.05, 0.95, 'IRAC3',  transform=ax_irac3.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=1.3, foreground='k')], size=8)
            ax_irac3.text(0.05, 0.95, 'IRAC3',  transform=ax_irac3.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=1.3, foreground='k')], size=8)
            ax_irac4.text(0.05, 0.95, 'IRAC4',  transform=ax_irac4.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=1.3, foreground='k')], size=8)

            ax_detc.text(0.04, 0.96, r'$\chi^2_{+}$ detection',  transform=ax_detc.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=2, foreground='k')], size=12)
            ax_segm.text(0.04, 0.96, r'segmentation',  transform=ax_segm.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=2, foreground='k')], size=12)
            ax_rgb.text(0.04, 0.96, r'NIRCam RGB',  transform=ax_rgb.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=2, foreground='k')], size=12)
            ax_mod.text(0.04, 0.96, 'Model (F277W)', transform=ax_mod.transAxes, color='k', va='top', ha='left', path_effects=[pe.withStroke(linewidth=2, foreground='0.9')], weight='bold', size=12)

            if Reff < 0.1:
                rstr = r"$R_{\rm eff} = " + fr"{Reff*1000:.1f} \pm {Reff_err*1000:.1f}$ mas" + '\n' + r'$n = ' + fr"{n_sersic:.2f} \pm {n_sersic_err:.2f}$"
            elif Reff < 0.5:
                rstr = r"$R_{\rm eff} = " + fr"{Reff:.3f} \pm {Reff_err:.3f}$ arcsec" + '\n' + r'$n = ' + fr"{n_sersic:.2f} \pm {n_sersic_err:.2f}$"
            else:
                rstr = r"$R_{\rm eff} = " + fr"{Reff:.2f} \pm {Reff_err:.2f}$ arcsec" + '\n' + r'$n = ' + fr"{n_sersic:.2f} \pm {n_sersic_err:.2f}$"
            ax_mod.text(0.04, 0.04, rstr, transform=ax_mod.transAxes, color='k', va='bottom', ha='left', path_effects=[pe.withStroke(linewidth=1, foreground='0.9')], size=10)

            # ax_sed = plt.subplot(gs[3, 1:-1])
            # ax_pz = ax_sed.inset_axes([0.03,0.7,0.3,0.27])
            img_axes = [ax_f814, ax_f115, ax_f150, ax_f277, ax_f444, ax_f770, ax_g, ax_r, ax_i, ax_z, ax_y, ax_Y, ax_J, ax_H, ax_K, ax_irac1, ax_irac3, ax_irac4, ax_rgb, ax_detc, ax_segm, ax_mod]
            for ax in img_axes:
                ax.set_aspect('equal')
                ax.tick_params(labelleft=False,labelbottom=False,left=False,right=False,top=False,bottom=False,which='both')

            ax_sed = plt.subplot(gs[13:-2,11:-13])
            ax_sed.set_xlabel('Observed Wavelength [µm]')
            ax_sed.set_ylabel('AB mag')

            ax_sed.tick_params(right=False, which='both')
            ax_sed.semilogx()
            ax_sed.set_xlim(0.3, 11)
            ax_sed.set_ylim(-30.5, -21.5)
            ax_sed.set_xticks([0.3,0.4,0.6,1.0,1.5,2.0,3.0,4.0,5,7,10],['0.3','0.4','0.6','1','1.5','2','3','4','5','7','10'])

            # def flux2mag(x):
            #     return -2.5*np.log10(x/3631e6)
            # def mag2flux(x):
            #     return 3631e6*np.power(10., -0.4*x)

            # secax1 = ax_sed.secondary_yaxis('right', functions=(flux2mag, mag2flux))
            # secax1.set_ylabel('AB mag')
            # secax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(2))
            # secax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
            # secax1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
            # secax1.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

            ax_pz = plt.subplot(gs[13:20,-12:])
            # ax_pz = ax_sed.inset_axes([0.03,0.7,0.38,0.27])
            ax_pz.set_xlabel('Redshift')#, size=9, labelpad=1)
            ax_pz.set_ylim(0, 1.1)
            ax_pz.set_xlim(0, 20)
            ax_pz.tick_params(labelleft=False,left=False,right=False,top=False,which='both')
            ax_pz.tick_params(direction='inout', which='both')#, labelsize=8, pad=1.5)
            ax_pz.tick_params(axis='x', which='major', length=5)
            ax_pz.tick_params(axis='x', which='minor', length=3)
            ax_pz.spines['left'].set_visible(False)
            ax_pz.spines['right'].set_visible(False)
            ax_pz.spines['top'].set_visible(False)

            ax_leg1 = plt.subplot(gs[21:25,-13:])
            ax_leg1.axis('off')
            ax_leg2 = plt.subplot(gs[25:,-13:])
            ax_leg2.axis('off')
            legend_handles, legend_labels = [], []


            ######################################################################################################################################
            ######################################################################################################################################
            ######################################################################################################################################

            # make cutouts
            if verbose: print('\t Making HST/JWST cutouts...')
            cutout_axes = [ax_f814, ax_f115, ax_f150, ax_f277, ax_f444, ax_f770]
            cutout_names = ['f814w','f115w','f150w','f277w','f444w','f770w']
            N = len(cutout_names)
            cutouts = []
            vm = []

            photflam, photplam = 6.99715969242424E-20, 8047.468423484849
            conversion1 = 3.33564e13 * (photplam)**2 * photflam
            conversion2 = 1e15*((0.03*u.arcsec)**2).to(u.sr).value # MJy/sr to nJy/pix
            conversion = [conversion1] + [conversion2]*5

            for i in range(N):
                try:
                    sci_cutout = Cutout2D(sci_dict[cutout_names[i]], coord, size=cutout_width, wcs=WCS(hdr_dict[cutout_names[i]]))
                    wht_cutout = Cutout2D(wht_dict[cutout_names[i]], coord, size=cutout_width, wcs=WCS(hdr_dict[cutout_names[i]]))
                    sci = sci_cutout.data * conversion[i]
                    err = err_table[cutout_names[i]]['std1'] / np.sqrt(wht_cutout.data)
                    snr = sci/err
                    wcs = sci_cutout.wcs

                    ps = wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value
                    size = np.shape(sci)[0]
                    extent = [-size*ps/2, size*ps/2, -size*ps/2, size*ps/2]
                    
                    cutouts.append((snr, wcs, extent))
                    vm.append(np.nanpercentile(snr[size//2-50:size//2+50,size//2-50:size//2+50], 95.0))
                except:
                    print(cutout_names[i], 'missing')
                    cutouts.append(None)
                    vm.append(0)
                    raise

            vmax = np.nanmax(vm)
            if vmax < vmax_min: vmax = vmax_min
            elif vmax > 50: vmax = 50
            vmin = -0.15*vmax
            tr = transforms.Affine2D().rotate_deg_around(0, 0, 20)
            for i in range(N):
                ############################## plot cutouts ##############################
                name = cutout_names[i]
                ax = cutout_axes[i]
                cutout = cutouts[i]
                if cutout is not None:
                    snr, wcs, extent = cutout
                    if np.all(np.isnan(snr)) or np.all(snr==0):
                        ax.imshow(np.zeros_like(snr), vmin=vmin,vmax=vmax, cmap=cmap, origin='lower', extent=extent)
                    else:
                        ax.imshow(snr, vmin=vmin, vmax=vmax, cmap=cmap, origin='lower', extent=extent, transform=tr + ax.transData)


                    # a, b, theta = cat['kron1_a'], cat['kron1_b'], np.degrees(cat['theta_image'])
                    # ax.add_patch(mpl.patches.Ellipse((0,0), width=2*a, height=2*b, angle=theta, facecolor='none', edgecolor='salmon', linestyle='-', linewidth=0.5))
                    # ax.add_patch(mpl.patches.Circle((0,0), radius=0.15, facecolor='none', edgecolor='w', linestyle='--', linewidth=0.5))
                    ax.set_xlim(-0.5*display_width.to(u.arcsec).value, 0.5*display_width.to(u.arcsec).value)
                    ax.set_ylim(-0.5*display_width.to(u.arcsec).value, 0.5*display_width.to(u.arcsec).value)
                    ax.spines[:].set_color(colors_dict[name])
                else:
                    ax.imshow(np.zeros((100,100)), vmin=vmin,vmax=vmax, cmap='Greys')
                
            
            if verbose: log('\t Making ground/IRAC cutouts...')
            cutout_names2 = ['g','r','i','z','y', 'Y','J','H','Ks','IRAC1','IRAC3','IRAC4']
            cutout_axes2 = [ax_g,ax_r,ax_i,ax_z,ax_y,ax_Y,ax_J,ax_H,ax_K,ax_irac1,ax_irac3,ax_irac4]
            N = len(cutout_names2)
            zp = 28.09
            conversion = 3631e9*np.power(10., -0.4*zp) 

            for i in range(N):
                name = cutout_names2[i]
                ax = cutout_axes2[i]

                sci_cutout = Cutout2D(sci_dict[cutout_names2[i]], coord, size=cutout_width, wcs=WCS(hdr_dict[cutout_names2[i]]))
                wht_cutout = Cutout2D(wht_dict[cutout_names2[i]], coord, size=cutout_width, wcs=WCS(hdr_dict[cutout_names2[i]]))
                sci = sci_cutout.data * conversion
                err = err_table[cutout_names2[i]]['std1'] / np.sqrt(wht_cutout.data)
                snr = sci/err
                wcs = sci_cutout.wcs

                ps = wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value
                size = np.shape(sci)[0]
                extent = [-size*ps/2, size*ps/2, -size*ps/2, size*ps/2]
                
                ax.imshow(snr, vmin=vmin, vmax=vmax, cmap='Greys', origin='lower', extent=extent)
                ax.set_xlim(-0.5*display_width.to(u.arcsec).value, 0.5*display_width.to(u.arcsec).value)
                ax.set_ylim(-0.5*display_width.to(u.arcsec).value, 0.5*display_width.to(u.arcsec).value)
                ax.spines[:].set_color(colors_dict[name])
                

            if verbose: log('\t Making detection cutout...')
            detec_cutout = Cutout2D(detec.data, coord, size=cutout_width*3/2, wcs=WCS(detec.header))
            wcs = detec_cutout.wcs
            ps = wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value
            size = np.shape(detec_cutout.data)[0]
            extent = [-size*ps/2, size*ps/2, -size*ps/2, size*ps/2]

            d = np.sqrt(detec_cutout.data)
            cen = np.shape(d)[0]//2
            vmax = np.nanpercentile(d[cen-20:cen+20,cen-20:cen+20],95)
            if vmax < 12: 
                vmax = 12
            ax_detc.imshow(d, extent=extent, vmin=1, vmax=vmax, cmap='Greys_r', transform=tr + ax_detc.transData)
            ax_detc.set_xlim(-display_width.to(u.arcsec).value*3/4, display_width.to(u.arcsec).value*3/4)
            ax_detc.set_ylim(-display_width.to(u.arcsec).value*3/4, display_width.to(u.arcsec).value*3/4)

            if verbose: log('\t Making segmentation cutout...')
            segm_cutout = Cutout2D(segm.data.astype(int), coord, size=cutout_width*3/2, wcs=WCS(segm.header))
            d = segm_cutout.data
            for i,unq_val in enumerate(np.sort(np.unique(d))):
                d[d==unq_val] = i

            from photutils.utils.colormaps import make_random_cmap
            from matplotlib.colors import to_rgba
            segm_cmap = make_random_cmap(len(np.unique(d)))
            segm_cmap.colors[0] = to_rgba('k')

            ax_segm.imshow(d, extent=extent, cmap=segm_cmap, transform=tr + ax_segm.transData, interpolation='none')
            ax_segm.set_xlim(-display_width.to(u.arcsec).value*3/4, display_width.to(u.arcsec).value*3/4)
            ax_segm.set_ylim(-display_width.to(u.arcsec).value*3/4, display_width.to(u.arcsec).value*3/4)


            if verbose: log('\t Making RGB cutout...')
            wcs = WCS(fits.getheader(paths.get_cosmos_web_full_f444w_filepath()))
            b = Cutout2D(np.flip(rgb_full[:,:,0],axis=0), coord, size=cutout_width*3/2, wcs=wcs)
            g = Cutout2D(np.flip(rgb_full[:,:,1],axis=0), coord, size=cutout_width*3/2, wcs=wcs)
            r = Cutout2D(np.flip(rgb_full[:,:,2],axis=0), coord, size=cutout_width*3/2, wcs=wcs)
            imrgb = np.dstack((r.data,g.data,b.data))

            wcs = r.wcs
            ps = wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value
            size = np.shape(r)[0]
            extent = [-size*ps/2, size*ps/2, -size*ps/2, size*ps/2]

            ax_rgb.imshow(imrgb, extent=extent)
            ax_rgb.set_xlim(-display_width.to(u.arcsec).value*3/4, display_width.to(u.arcsec).value*3/4)
            ax_rgb.set_ylim(-display_width.to(u.arcsec).value*3/4, display_width.to(u.arcsec).value*3/4)
            x0, y0 = -display_width.to(u.arcsec).value*3/4, -display_width.to(u.arcsec).value*3/4
            x0 += 0.05*display_width.to(u.arcsec).value
            y0 += 0.07*display_width.to(u.arcsec).value
            x0 += scalebar_size.to(u.arcsec).value/2
            ax_rgb.errorbar(x0, y0, xerr=scalebar_size.to(u.arcsec).value/2, yerr=0, linewidth=1, color='w', capsize=3, capthick=1, marker='none')
            y0 += 0.025*display_width.to(u.arcsec).value
            ax_rgb.annotate(f'{scalebar_size.to(u.arcsec).value:.0f}"', (x0, y0), ha='center', va='bottom', color='w', fontsize=8, path_effects=[pe.withStroke(linewidth=0.5, foreground='k')])

            if verbose: log('\t Making F277W model cutout...')
            model_cutout = Cutout2D(model.data, coord, size=cutout_width*3/2, wcs=WCS(model.header))
            wcs = model_cutout.wcs
            ps = wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value
            size = np.shape(model_cutout.data)[0]
            extent = [-size*ps/2, size*ps/2, -size*ps/2, size*ps/2]

            d = model_cutout.data
            cen = np.shape(d)[0]//2
            vmax = np.nanpercentile(d[cen-10:cen+10,cen-10:cen+10],95)
            if vmax < 0: vmax = np.nanpercentile(d[cen-20:cen+20,cen-20:cen+20],99)
            if vmax < 0: vmax = np.nanmax(d)
            ax_mod.imshow(d, extent=extent, vmin=-vmax/6, vmax=vmax, cmap='Greys', transform=tr + ax_mod.transData)
            ax_mod.set_xlim(-display_width.to(u.arcsec).value*3/4, display_width.to(u.arcsec).value*3/4)
            ax_mod.set_ylim(-display_width.to(u.arcsec).value*3/4, display_width.to(u.arcsec).value*3/4)


            def plot_data(ax, wav, wav_min, wav_max, flux, flux_err, colors, sizes, zorders, 
                          plot_xerr=True, 
                          annotate=True, 
                          label=None):
                min_mag, max_mag = 0, -30

                colors = np.array(colors)
                wav = wav[flux_err>0]
                wav_min = wav_min[flux_err>0]
                wav_max = wav_max[flux_err>0]
                colors = colors[flux_err>0]
                flux = flux[flux_err>0]
                flux_err = flux_err[flux_err>0]
                
                snrs = flux/flux_err
                flux_uplim = 2.5*np.log10((flux+2*flux_err)/3631e6)
                flux = 2.5*np.log10(flux/3631e6)
                flux_upper_err = 2.5*np.log10((flux+flux_err)/3631e6) - flux
                flux_lower_err = flux - 2.5*np.log10((flux-flux_err)/3631e6)
                for w, w1, w2, f, f_up_err, f_lo_err, f_up_lim, snr, c, s, z in zip(wav, wav_min, wav_max, flux, flux_upper_err, flux_lower_err, flux_uplim, snrs, colors, sizes, zorders):
                    # print(f'{w:.2f} {f:.2f} {snr:.2f}')
                    if snr > 1.5:
                        ax.errorbar(w, f, yerr=[[f_lo_err],[f_up_err]], linewidth=0, marker='o', ms=s, 
                                    mfc=c, mec=c, elinewidth=1, ecolor=c, mew=1, capthick=0, capsize=0, zorder=z)
                        if plot_xerr:
                            ax.errorbar(w, f, xerr=[[w-w1],[w2-w]], linewidth=0, marker='none',
                                        elinewidth=1, ecolor=c, mew=1, capthick=0, capsize=0, zorder=z)
                        if annotate: 
                            ax.annotate(fr'${snr:.1f}\sigma$', (w, 1.15*(f+f_err)), ha='center', va='bottom', color=c, fontsize=6, bbox=dict(facecolor='w', edgecolor='none', pad=0.01, alpha=0.7), zorder=z)
                        if not np.isinf(f):
                            min_mag = np.nanmin([min_mag, f])
                            max_mag = np.nanmax([max_mag, f])
                    else:
                        ax.errorbar(w, f_up_lim, yerr=0.3, uplims=True, linewidth=0,
                                    mfc='none', mec=c, elinewidth=1, ecolor=c, mew=1, capthick=1, capsize=s/2, zorder=z)
                        if plot_xerr:
                            ax.errorbar(w, f_up_lim, xerr=[[w-w1],[w2-w]], linewidth=0, marker='none', 
                                        elinewidth=1, ecolor=c, mew=1, capthick=0, capsize=0, zorder=z)
                        if annotate: ax.annotate(fr'${snr:.1f}\sigma$', (w, 1.15*2*f_err), ha='center', va='bottom', color=c, fontsize=6, bbox=dict(facecolor='w', edgecolor='none', pad=0.01, alpha=0.7), zorder=z)
                        if not np.isinf(f_up_lim):
                            min_mag = np.nanmin([min_mag, f_up_lim])
                            max_mag = np.nanmax([max_mag, f_up_lim])
                return min_mag, max_mag


            if verbose: log('\t Plotting photometry...')
            wav, wav_min, wav_max = filters.wav.to(u.micron).value, filters.wav_min.to(u.micron).value, filters.wav_max.to(u.micron).value
            flux_cols = ['FLUX_MODEL_HST-F814W','FLUX_MODEL_F115W','FLUX_MODEL_F150W','FLUX_MODEL_F277W','FLUX_MODEL_F444W','FLUX_MODEL_F770W',
                        'FLUX_MODEL_CFHT-u','FLUX_MODEL_HSC-g','FLUX_MODEL_HSC-r','FLUX_MODEL_HSC-i','FLUX_MODEL_HSC-z','FLUX_MODEL_HSC-y',
                        'FLUX_MODEL_UVISTA-Y','FLUX_MODEL_UVISTA-J','FLUX_MODEL_UVISTA-H','FLUX_MODEL_UVISTA-Ks',
                        'FLUX_MODEL_SC-IB427', 'FLUX_MODEL_SC-IA484', 'FLUX_MODEL_SC-IB505', 'FLUX_MODEL_SC-IA527', 
                        'FLUX_MODEL_SC-IB574', 'FLUX_MODEL_SC-IA624', 'FLUX_MODEL_SC-IA679', 'FLUX_MODEL_SC-IB709', 
                        'FLUX_MODEL_SC-IA738', 'FLUX_MODEL_SC-IA767', 'FLUX_MODEL_SC-IB827', 'FLUX_MODEL_SC-NB711', 
                        'FLUX_MODEL_SC-NB816', 'FLUX_MODEL_HSC-NB0816', 'FLUX_MODEL_HSC-NB0921', 'FLUX_MODEL_HSC-NB1010', 
                        'FLUX_MODEL_UVISTA-NB118', 
                        'FLUX_MODEL_IRAC-ch1','FLUX_MODEL_IRAC-ch3','FLUX_MODEL_IRAC-ch4']


            flux_err_cols = [f.replace('FLUX_','FLUX_ERR-CAL_') for f in flux_cols]
            flux = np.array([cat[f][0] for f in flux_cols])*1e29
            flux_err = np.array([cat[f][0] for f in flux_err_cols])*1e29


            labels = ['F814W','F115W','F150W','F277W','F444W','F770W','$u$','$g$','$r$','$i$','$z$','$y$','$Y$','$J$','$H$',r'$K_s$']
            labels += ['']*17
            labels += ['[3.6]', '[5.8]', '[8.0]']
            yoff = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0,0,0,0,0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0.1,-0.1]
            xoff = [1.01,1,1,1,1,1,1,1,1,0.95,1.03,1,1,1.03,1.03,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1,1.2]
            sizes = [7]*6 + [4]*10 + [3]*20
            zorders = 1000-np.arange(len(colors))

            
            # maxflux = 0
            
            # maxflux = np.nanmax([maxflux,np.nanmax(flux)*2])
            min_mag, max_mag = plot_data(ax_sed, wav, wav_min, wav_max, flux, flux_err, colors, sizes, zorders, annotate=False, plot_xerr=True)
            
            ymin, ymax = ax_sed.get_ylim()
            if max_mag > -24: 
                ymax = max_mag + 3
            
            if min_mag > -28:
                ymin = min_mag - 2

            ax_sed.set_ylim(ymin, ymax)

                

            # if maxflux > 1e2:
            #     ax_sed.set_ylim(2e-3, 3*maxflux)

            ymin, ymax = ax_sed.get_ylim()
            for i in range(len(filters)):
                w, t = filters[i].T
                scale = 0.6*sizes[i]/7
                t = scale*t/np.max(t) + ymin
                ax_sed.fill_between(w/1e4, ymin, t, edgecolor='none', facecolor=colors[i], alpha=0.08*sizes[i]/4)
                ax_sed.plot(w/1e4, t, color=colors[i], linewidth=0.5, alpha=0.4*sizes[i]/4)
                ax_sed.annotate(labels[i], (wav[i]*xoff[i], np.max(t)*yoff[i]), ha='center', va='bottom', color=colors[i], size=6)


            if lephare_spec_path is not None:
                lph_color1 = 'darkorange'
                lph_color2 = 'steelblue'
                lph_color3 = 'darkmagenta'
                lph_color4 = '0.7'

                from .lephare_helpers import LephareResult
                lph = LephareResult.read(os.path.join(lephare_spec_path, f'Id{ID}.0.spec'))
                lph.pz.normalize()
                ax_pz.plot(lph.pz.zgrid, lph.pz.Pz, color='k', linewidth=0.8)
                ax_pz.fill_between(lph.pz.zgrid, lph.pz.Pz, edgecolor='none', facecolor='k', alpha=0.25)
                # ax_pz.plot(lph.pz.zgrid, lph.pz.Pz_bayesian, color=lph_color, linewidth=0.5, linestyle='--')

                cdf = np.cumsum(lph.pz.Pz)/np.sum(lph.pz.Pz)
                zmin = np.interp(0.01, cdf, lph.pz.zgrid)
                zmax = np.interp(0.99, cdf, lph.pz.zgrid)
                z16 = np.interp(0.16, cdf, lph.pz.zgrid)
                z50 = np.interp(0.50, cdf, lph.pz.zgrid)
                z84 = np.interp(0.84, cdf, lph.pz.zgrid)
                if any(np.isnan(zmin,zmax,z16,z50,z84)) or any(np.isinf(zmin,zmax,z16,z50,z84)):
                    ax_pz.set_xlim(0, 20)
                elif zmin > 0.75:
                    ax_pz.set_xlim(zmin-0.3, zmax+0.3)
                else:
                    ax_pz.set_xlim(0, zmax+0.5)
                # if zmax < 3: 
                #     ax_pz.set_xlim(0, 3)
                # if zmax < 4: 
                #     ax_pz.set_xlim(0, 4)
                # if zmax < 5: 
                #     ax_pz.set_xlim(0, 5)
                # elif zmax < 6: 
                #     ax_pz.set_xlim(0, 6)
                # elif zmax < 7: 
                #     ax_pz.set_xlim(0, 7)
                # elif zmax < 8: 
                #     ax_pz.set_xlim(0, 8)
                # elif zmax < 9: 
                #     ax_pz.set_xlim(0, 9)
                # elif zmax < 10: 
                #     ax_pz.set_xlim(0, 10)


                ax_pz.errorbar(z50, 1.05, xerr=[[z50-z16],[z84-z50]], linewidth=0, marker='o', ms=6, ecolor='k', mfc='k', mec='k', elinewidth=1, mew=1, capthick=1, capsize=3, zorder=1000)



                phot = lph.phot.get_filters(filters.names[:-3])
                w = wav[:-3]
                s = sizes[:-3]
                for i in range(len(phot)):
                    si = s[i] + 2
                    e = ax_sed.errorbar(w[i], 2.5*np.log10(phot.model[i]/3631e6), marker='s', mew=1, mec=lph_color1, mfc='none', ms=si, zorder=100)

                p, = ax_sed.plot(lph.models['GAL-1']['wav_obs']/1e4, 2.5*np.log10(lph.models['GAL-1']['fnu']/3631e6), color=lph_color1, alpha=0.5, zorder=-999, linewidth=0.8)
                chi2 = lph.models['GAL-1']['chi2']
                zbest = lph.models['GAL-1']['zphot']
                ax_pz.axvline(zbest, color=lph_color1, linestyle='--', linewidth=1)
                legend_labels.append(fr'Galaxy best-fit ($z={zbest:.2f}$, $\chi^2 = {chi2:.1f}$)')
                legend_handles.append((e, p, ))

                if lph.models['GAL-2']['zphot'] > 0:
                    p, = ax_sed.plot(lph.models['GAL-2']['wav_obs']/1e4, 2.5*np.log10(lph.models['GAL-2']['fnu']/3631e6), color=lph_color2, alpha=0.4, zorder=-999, linewidth=0.5)
                    chi2 = lph.models['GAL-2']['chi2']
                    zbest = lph.models['GAL-2']['zphot']
                    ax_pz.axvline(zbest, color=lph_color2, linestyle='--', linewidth=1)
                    legend_labels.append(fr'Secondary ($z={zbest:.2f}$, $\chi^2 = {chi2:.1f}$)')
                    legend_handles.append(p)

                if lph.models['QSO']['zphot'] > 0:
                    ax_sed.plot(lph.models['QSO']['wav_obs']/1e4, 2.5*np.log10(lph.models['QSO']['fnu']/3631e6), color=lph_color3, alpha=0.4, zorder=-999, linewidth=0.5)
                    p, = ax_sed.plot([1e-5,2e-5], [-50,-50], color=lph_color3, alpha=0.6, zorder=-999, linewidth=0.8)
                    chi2 = lph.models['QSO']['chi2']
                    zbest = lph.models['QSO']['zphot']
                    ax_pz.axvline(zbest, color=lph_color3, linestyle='--', linewidth=1)
                    legend_labels.append(fr'QSO ($z={zbest:.2f}$, $\chi^2 = {chi2:.1f}$)')
                    legend_handles.append(p)

                ax_sed.plot(lph.models['STAR']['wav_obs']/1e4, 2.5*np.log10(lph.models['STAR']['fnu']/3631e6), color=lph_color4, alpha=0.4, zorder=-999, linewidth=0.5)
                p, = ax_sed.plot([1e-5,2e-5], [-50,-50], color=lph_color4, alpha=0.7, zorder=-999, linewidth=0.8)
                chi2 = lph.models['STAR']['chi2']
                legend_labels.append(fr'Star ($\chi^2 = {chi2:.1f}$)')
                legend_handles.append(p)

                if np.min([z50-z16, z84-z50]) < 0.1:
                    s = r'LePhare $z_{\rm phot} = %.2f^{+%.2f}_{-%.2f}$' % (z50, z84-z50, z50-z16)
                else:
                    s = r'LePhare $z_{\rm phot} = %.1f^{+%.1f}_{-%.1f}$' % (z50, z84-z50, z50-z16)

                s += '\n'
                m16, m50, m84 = cat['LP_mass_l68_PDF'][0], cat['LP_mass_med_PDF'][0], cat['LP_mass_u68_PDF'][0]
                s16, s50, s84 = cat['LP_ssfr_l68_PDF'][0], cat['LP_ssfr_med_PDF'][0], cat['LP_ssfr_u68_PDF'][0]
                if np.min([m50-m16, m84-m50]) < 0.1:
                    s += r'$\log M_\star/M_\odot = %.2f^{+%.2f}_{-%.2f}$' % (m50, m84-m50, m50-m16)
                else:
                    s += r'$\log M_\star/M_\odot = %.1f^{+%.1f}_{-%.1f}$' % (m50, m84-m50, m50-m16)
                s += '\n'
                if np.min([s50-s16, s84-s50]) < 0.1:
                    s += r'$\log {\rm sSFR}/M_\odot\,{\rm yr}^{-1} = %.2f^{+%.2f}_{-%.2f}$' % (s50, s84-s50, s50-s16)
                else:
                    s += r'$\log {\rm sSFR}/M_\odot\,{\rm yr}^{-1} = %.1f^{+%.1f}_{-%.1f}$' % (s50, s84-s50, s50-s16)
                ax_leg1.annotate(s, (0.03, 0.), ha='left', va='bottom', color='k', xycoords='axes fraction', fontsize=12)

                # ax_leg1.annotate(s, (0.03, 0.2), ha='left', va='center', color='k', xycoords='axes fraction', fontsize=12)


            leg1 = ax_leg2.legend(handles=tuple(legend_handles), labels=tuple(legend_labels), frameon=True, loc='upper left')


            def labels(x, pos):
                return f'{int(round(-x,0))}'
            ax_sed.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(labels))

            if out_format == 'pdf':
                plt.savefig(outpath)
            elif out_format == 'png':
                assert dpi is not None, 'dpi must be specified for png output'
                plt.savefig(outpath, dpi=300)
            else:
                raise Exception("out_format must be 'pdf' or 'png'")

            plt.close()





if __name__ == '__main__':
    catalog_filename = 'COSMOSWeb_master_v3.1.0-sersic-cgs_err-calib_LePhare+CIGALE.fits'
    catalog_shortname = 'v3.1.0'

    if host.hostname == 'patrick':
        logo = '/Users/hba423/Library/CloudStorage/Dropbox/presentations/assets/cosmos_web_logos/cosmosweb_black_horz_transparent.png'
        catalog_path = '/data/COSMOS-Web/catalogs/'
        outdir = '/research/COSMOS-Web/inspection_plots/'
        lephare_spec_path = '/Users/hba423/Downloads/'

    elif host.hostname == 'candide':
        logo = '/n23data2/hakins/exchg/COSMOS-Web/cwlogo_black.png'
        catalog_path = '/n17data/shuntov/COSMOS-Web/Catalogs/'
        outdir = '/n23data2/hakins/exchg/COSMOS-Web/inspec_plots/'
        lephare_spec_path = '/home/ilbert/n07data/COSMOS-Web/photoz_MASTER_v3.1.0/PHOTOZ_BC03/SPEC_v3.1.0/'
            

    which_ids = sys.argv[1]
    if which_ids == 'all':
        log('Starting plot generation for all IDs in catalog...')
        
        # figure out which tile each ID is in
        f = fits.getdata(os.path.join(catalog_path, catalog_filename))

        IDs_all = np.array(f['ID_SE++'], dtype=int)
        tiles_all = np.array(f['TILE'], dtype=str)

        # split IDs_all array into dict of arrays, one for each tile
        IDs = {}
        for tile in np.unique(tiles_all):
            IDs[tile] = IDs_all[tiles_all==tile]

    elif which_ids.startswith('A') or which_ids.startswith('B'):
        tile = which_ids
        log(f'Starting plot generation for all IDs in tile {tile}...')
        f = fits.getdata(os.path.join(catalog_path, catalog_filename))
        IDs_all = np.array(f[f['TILE']==tile]['ID_SE++'], dtype=int)
        IDs = {tile: IDs_all}

    else:
        f = fits.getdata(os.path.join(catalog_path, catalog_filename))
        # raise Exception('Invalid input for which_ids')
        IDs_all = np.array(which_ids.split(','), dtype=int)
        tiles_all = []
        for ID in IDs_all:
            fi = f[f['ID_SE++']==ID]
            tiles_all.append(fi['TILE'][0])
        tiles_all = np.array(tiles_all)
        
        IDs = {}
        for tile in np.unique(tiles_all):
            IDs[tile] = IDs_all[tiles_all==tile]
    
    main(IDs, 
         logo=logo, 
         catalog_path=catalog_path, 
         catalog_filename=catalog_filename, 
         catalog_shortname=catalog_shortname, 
         outdir=outdir, 
         out_format='png', dpi=400, 
         overwrite=False, 
         display_width=4*u.arcsec, 
         vmax_min=8,
         cmap='Greys', 
         lephare_spec_path = lephare_spec_path, 
         verbose=False)