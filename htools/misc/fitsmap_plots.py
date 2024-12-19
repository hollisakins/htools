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


def main(IDs, 
         tiles,
         logo = None,
         catalog_path = None,
         catalog_filename = None,
         catalog_shortname = None,
         outdir='inspection_plots/',
         out_format='pdf',
         dpi=None,
         overwrite=True,
         display_width = 4*u.arcsec,
         cutout_width = 6*u.arcsec,
         vmax_in = 10,
         cmap = 'Greys',
         lephare_spec_path=None,
         verbose=False,
    ):

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
            print(f'Reading in band {band}')
            sci_dict[band] = io.load_sci('cosmos-web', band, tile=tile) 
            hdr_dict[band] = io.load_hdr('cosmos-web', band, tile=tile) 
            wht_dict[band] = io.load_wht('cosmos-web', band, tile=tile) 
            # sci_filepath, hdu_index = get_filepath('cosmos-web', band, 'sci', tile=tile) # fits.open(sci_filepath)[hdu_index]
            # wht_filepath, hdu_index = get_filepath('cosmos-web', band, 'wht', tile=tile) #fits.open(wht_filepath)[hdu_index]

        detec = io.load_cosmos_web_detec(tile)
        segm = io.load_cosmos_web_segm(tile, catalog_version='v1.3')

        colors = {
            'f814w': 'darkmagenta',
            'f115w': '#0088e7',
            'f150w': '#03a1a1',
            'f277w': '#83b505',
            'f444w': '#ab0202',
            'f770w': '#e74001'
        }

        for i in tqdm.tqdm(range(len(IDs[tile]))):
            ID = IDs[tile][i]

            ######################################################################################################################################
            cat = catalog[catalog['ID_SE++']==ID]
            ra = cat['RA_MODEL'][0]
            dec = cat['DEC_MODEL'][0]
            coord = SkyCoord(ra=ra, dec=dec, unit=u.deg)

            outfilename = f"cosmos-web_sed_{catalog_shortname.replace('-','_')}_{ID}"
            outpath = os.path.join(outdir, outfilename)
            if out_format == 'png':
                outpath += '.png'
            elif out_format == 'pdf':
                outpath += '.pdf'

            if not overwrite: 
                if os.path.exists(outpath):
                    if verbose: print(f'\t Skipping, plot already exists at {outpath}.')


            ######################################################################################################################################
            fig = plt.figure(figsize=(12, 8.2), constrained_layout=False)
            gs = mpl.gridspec.GridSpec(ncols=15, nrows=30, width_ratios=[1.01]*3 + [1]*12, figure=fig)
            gs.update(hspace=0.08, wspace=0.05, left=0.03, bottom=0.03, right=1-0.03, top=1-0.03)

            if logo is not None:
                logo_img = Image.open(logo)
                logo_img = logo_img.transpose(Image.FLIP_TOP_BOTTOM)

                # Create an inset axes
                axins = fig.add_axes([0.03, 0.88, 0.17, 0.11])
                axins.imshow(logo_img, alpha=1, aspect='auto')
                # axins.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labeltop=False, labelleft=False, labelright=False)
                axins.axis('off')


            fig.text(0.97, 0.977, f'ID ({catalog_shortname}): {ID}', va='top', ha='right',fontsize=12, weight='bold')
            coordstring = coord.to_string('hmsdms', precision=2).split(' ')
            fig.text(0.97, 0.952, f'COSMOS-Web tile: {tile}', va='top', ha='right',fontsize=12)
            fig.text(0.97, 0.927, f'RA, Dec: ({coordstring[0]}, {coordstring[1]})', va='top', ha='right',fontsize=12)
            fig.text(0.97, 0.902, f'({coord.ra.value:.7f}, {coord.dec.value:.6f})', va='top', ha='right',fontsize=12)


            ax_rgb = plt.subplot(gs[3:12,0:3])
            ax_detc = plt.subplot(gs[12:21,0:3])
            ax_segm = plt.subplot(gs[21:30,0:3])

            ax_f814 =  plt.subplot(gs[3:9,3:5])
            ax_f115 =  plt.subplot(gs[3:9,5:7])
            ax_f150 =  plt.subplot(gs[3:9,7:9])
            ax_f277 =  plt.subplot(gs[3:9,9:11])
            ax_f444 =  plt.subplot(gs[3:9,11:13])
            ax_f770 =  plt.subplot(gs[3:9,13:15])

            ax_g =     plt.subplot(gs[9:12,3])
            ax_r =     plt.subplot(gs[9:12,4])
            ax_i =     plt.subplot(gs[9:12,5])
            ax_z =     plt.subplot(gs[9:12,6])
            ax_y =     plt.subplot(gs[9:12,7])
            ax_Y =     plt.subplot(gs[9:12,8])
            ax_J =     plt.subplot(gs[9:12,9])
            ax_H =     plt.subplot(gs[9:12,10])
            ax_K =     plt.subplot(gs[9:12,11])
            ax_irac1 = plt.subplot(gs[9:12,12])
            ax_irac3 = plt.subplot(gs[9:12,13])
            ax_irac4 = plt.subplot(gs[9:12,14])

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

            ax_detc.text(0.05, 0.95, r'$\chi^2_{+}$ detection',  transform=ax_detc.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=2, foreground='k')], size=12)
            ax_segm.text(0.05, 0.95, r'segmentation',  transform=ax_segm.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=2, foreground='k')], size=12)
            ax_rgb.text(0.05, 0.95, r'NIRCam RGB',  transform=ax_rgb.transAxes, color='white', va='top', ha='left', path_effects=[pe.withStroke(linewidth=2, foreground='k')], size=12)

            # ax_sed = plt.subplot(gs[3, 1:-1])
            # ax_pz = ax_sed.inset_axes([0.03,0.7,0.3,0.27])
            img_axes = [ax_f814, ax_f115, ax_f150, ax_f277, ax_f444, ax_f770, ax_g, ax_r, ax_i, ax_z, ax_y, ax_Y, ax_J, ax_H, ax_K, ax_irac1, ax_irac3, ax_irac4, ax_rgb, ax_detc, ax_segm]
            for ax in img_axes:
                ax.set_aspect('equal')
                ax.tick_params(labelleft=False,labelbottom=False,left=False,right=False,top=False,bottom=False,which='both')

            ax_sed = plt.subplot(gs[13:-2,4:-1])
            ax_sed.set_xlabel('Observed Wavelength [µm]')
            ax_sed.set_ylabel('Flux Density [µJy]')

            ax_sed.tick_params(right=False, which='both')
            ax_sed.set_xticks([0.3,0.4,0.6,1.0,1.5,2.0,3.0,4.0,5,7,10,15,20],['0.3','0.4','0.6','1','1.5','2','3','4','5','7','10','15','20'])
            ax_sed.set_xlim(0.25, 14)
            ax_sed.set_ylim(2e-3, 1e2)
            ax_sed.loglog()

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

            ax_pz = ax_sed.inset_axes([0.03,0.7,0.38,0.27])
            ax_pz.set_xlabel('Redshift', size=9, labelpad=1)
            ax_pz.set_ylim(0, 1.05)
            ax_pz.set_xlim(0, 20)
            ax_pz.tick_params(labelleft=False,left=False,right=False,top=False,which='both')
            ax_pz.tick_params(direction='inout', which='both', labelsize=8, pad=1.5)
            ax_pz.spines['left'].set_visible(False)
            ax_pz.spines['right'].set_visible(False)
            ax_pz.spines['top'].set_visible(False)



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
                    vm.append(np.nanpercentile(snr[size//2-50:size//2+50,size//2-50:size//2+50], 99.0))
                except:
                    print(cutout_names[i], 'missing')
                    cutouts.append(None)
                    vm.append(0)

            if vmax_in=='auto':
                vmax = np.nanmax(vm)
                if vmax < 8: vmax = 8
                if vmax > 30: vmax = 30
                vmin = -0.15*vmax
            else:
                assert (type(vmax_in)==int) or (type(vmax_in)==float)
                vmax = vmax_in
                vmin = -3
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
                    ax.spines[:].set_color(colors[name])
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
                
                ax.imshow(snr, vmin=-3, vmax=8, cmap='Greys', origin='lower', extent=extent)
                ax.set_xlim(-0.5*display_width.to(u.arcsec).value, 0.5*display_width.to(u.arcsec).value)
                ax.set_ylim(-0.5*display_width.to(u.arcsec).value, 0.5*display_width.to(u.arcsec).value)
                

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
            wcs = WCS(fits.getheader('/Users/hba423/fitsmap/data3/CW_f444w_60mas_tot_v8.fits'))
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

            def plot_data(ax, wav, wav_min, wav_max, flux, flux_err, colors, zorder=10, annotate=True, label=None):
                colors = np.array(colors)
                wav = wav[flux_err>0]
                wav_min = wav_min[flux_err>0]
                wav_max = wav_max[flux_err>0]
                colors = colors[flux_err>0]
                flux = flux[flux_err>0]
                flux_err = flux_err[flux_err>0]
                for w, w1, w2, f, f_err, c, z in zip(wav, wav_min, wav_max, flux, flux_err, colors, zorder):
                    # if verbose: print(f'\t {w:.2f}, {f:.2f}, {f_err:.2f}, {f/f_err:.1f}')
                    s = f/f_err
                    if s > 1.5:
                        ax.errorbar(w, f, yerr=f_err, xerr=[[w-w1],[w2-w]], linewidth=0, marker='o', ms=6, 
                                    mfc=c, mec=c, elinewidth=1, ecolor=c, capthick=1, capsize=2, zorder=z)
                        if annotate: ax.annotate(fr'${s:.1f}\sigma$', (w, 1.15*(f+f_err)), ha='center', va='bottom', color=c, fontsize=6, bbox=dict(facecolor='w', edgecolor='none', pad=0.01, alpha=0.7), zorder=z)
                    else:
                        ax.errorbar(w, 2*f_err, yerr=0.5*f_err, xerr=[[w-w1],[w2-w]], uplims=True, linewidth=0,
                                    mfc='none', mec=c, elinewidth=1, ecolor=c, capthick=1, capsize=2, zorder=z)
                        if annotate: ax.annotate(fr'${s:.1f}\sigma$', (w, 1.15*2*f_err), ha='center', va='bottom', color=c, fontsize=6, bbox=dict(facecolor='w', edgecolor='none', pad=0.01, alpha=0.7), zorder=z)
                if label is not None:
                    ax.errorbar(100, 1, yerr=1, xerr=1, linewidth=0, marker='s', ms=6, 
                                mfc='none', mec=c, elinewidth=1, ecolor=c, capthick=1, capsize=2, zorder=zorder, label=label)


            if verbose: log('\t Plotting photometry...')
            filters = Filters(['f814w','f115w','f150w','f277w','f444w','f770w',
                            'cfht_u', 'hsc_g', 'hsc_r', 'hsc_i', 'hsc_z', 'hsc_y', 
                            'uvista_Y', 'uvista_J', 'uvista_H', 'uvista_Ks',
                            'irac_ch1', 'irac_ch3', 'irac_ch4'])
            wav, wav_min, wav_max = filters.wav.to(u.micron).value, filters.wav_min.to(u.micron).value, filters.wav_max.to(u.micron).value
            flux_cols = ['FLUX_MODEL_HST-F814W','FLUX_MODEL_F115W','FLUX_MODEL_F150W','FLUX_MODEL_F277W','FLUX_MODEL_F444W','FLUX_MODEL_F770W',
                        'FLUX_MODEL_CFHT-u','FLUX_MODEL_HSC-g','FLUX_MODEL_HSC-r','FLUX_MODEL_HSC-i','FLUX_MODEL_HSC-z','FLUX_MODEL_HSC-y',
                        'FLUX_MODEL_UVISTA-Y','FLUX_MODEL_UVISTA-J','FLUX_MODEL_UVISTA-H','FLUX_MODEL_UVISTA-Ks',
                        'FLUX_MODEL_IRAC-ch1','FLUX_MODEL_IRAC-ch3','FLUX_MODEL_IRAC-ch4']
            flux_err_cols = [f.replace('FLUX_','FLUX_ERR-CAL_') for f in flux_cols]
            flux = np.array([cat[f][0] for f in flux_cols])*1e29
            flux_err = np.array([cat[f][0] for f in flux_err_cols])*1e29

            cs = [colors['f814w'], colors['f115w'], colors['f150w'], colors['f277w'], colors['f444w'], colors['f770w']] + ['k'] * 10 + ['0.8'] * 3
            labels = ['F814W','F115W','F150W','F277W','F444W','F770W','$u$','$g$','$r$','$i$','$z$','$y$','$Y$','$J$','$H$',r'$K_s$', '[3.6]', '[5.8]', '[8.0]']
            yoff = [1.05,1.05,1.05,1.05,1.05,1.05,1,1,1,1,1,1,1,1,1,1,1,1.05,0.9]
            xoff = [1.01,1,1,1,1,1,1,1,1,0.95,1.03,1,1,1.03,1.03,1,1,1,1.2]
            zorders = 100-np.arange(len(cs))

            maxflux = 0
            
            maxflux = np.nanmax([maxflux,np.nanmax(flux)*2])
            plot_data(ax_sed, wav, wav_min, wav_max, flux, flux_err, cs, zorder=zorders, annotate=False)

            if maxflux > 1e2:
                ax_sed.set_ylim(2e-3, 3*maxflux)

            ymin, ymax = ax_sed.get_ylim()
            for i in range(len(filters)):
                w, t = filters[i].T
                t = np.power(10., 0.3*t/np.max(t)+np.log10(ymin))
                ax_sed.fill_between(w/1e4, ymin, t, edgecolor='none', facecolor=cs[i], alpha=0.05)
                ax_sed.plot(w/1e4, t, color=cs[i], linewidth=0.5, alpha=0.4)
                ax_sed.annotate(labels[i], (wav[i]*xoff[i], np.max(t)*yoff[i]), ha='center', va='bottom', color=cs[i], size=6)


            if lephare_spec_path is not None:
                lph_color = 'steelblue'
                from .lephare_helpers import LephareResult
                lph = LephareResult.read(os.path.join(lephare_spec_path, f'ID{ID}.0.spec'))
                lph.pz.normalize()
                ax_pz.plot(lph.pz.zgrid, lph.pz.Pz, color=lph_color, linewidth=1)
                ax_pz.fill_between(lph.pz.zgrid, lph.pz.Pz, edgecolor='none', facecolor=lph_color, alpha=0.1)
                ax_pz.plot(lph.pz.zgrid, lph.pz.Pz_bayesian, color=lph_color, linewidth=0.5, linestyle='--')

                phot = lph.phot.get_filters(filters.names[:-3])
                w = wav[:-3]
                for i in range(len(phot)):
                    ax_sed.errorbar(w[i], phot.model[i], marker='o', mew=2, mec=lph_color, mfc='none', ms=10, zorder=999)

                ax_sed.plot(lph.models['GAL-1']['wav_obs']/1e4, lph.models['GAL-1']['fnu'], color=lph_color, alpha=0.4, zorder=-999, linewidth=0.5)


            # if brisket_run is not None:
            #     import pickle
            #     try:
            #         with open(f'{brisket_run}/{nickname}_sed.pickle', 'rb') as f:
            #             sed = pickle.load(f)
            #     except:
            #         with open(f'{brisket_run}/{ID}_sed.pickle', 'rb') as f:
            #             sed = pickle.load(f)
                

            #     from getdist import MCSamples
            #     samples = np.array([sed['z']]).T
            #     settings = {
            #                 "contours":[0.68, 0.95, 0.99], 
            #                 "range_ND_contour":1, 
            #                 "range_confidence":0.001,
            #                 "fine_bins":200,
            #                 "fine_bins_2d":80,
            #                 "smooth_scale_1D":0.1,
            #                 "smooth_scale_2D":0.5,
            #                 "tight_gap_fraction":0.15
            #                 }
            #     samples = MCSamples(samples=samples, names=['z'], settings=settings)
            #     density = samples.get1DDensity('z')
            #     x = np.linspace(0, 20, 1000)
            #     y = density(x)
            #     y[y<0] = 0
            #     y = y/np.max(y)
            #     ax_pz.fill_between(x, y, facecolor='royalblue', edgecolor='none', alpha=0.2)
            #     ax_pz.plot(x, y, color='b')
            #     z0 = x[np.argmax(y)]
            #     if z0<1:
            #         z0 = np.median(sed['z'])
            #     ax_pz.axvline(z0, color='b', linestyle='--')
            #     # t = Table({'z': x, 'Pz': y})
            #     # t.write(f'{out_sed_dir}/{ID}_brisket_pz.txt', format='ascii', overwrite=True)
            #     # zmin = np.percentile(np.random.choice(x, p=y/np.sum(y), size=1000), 16)
            #     # zmax = np.percentile(np.random.choice(x, p=y/np.sum(y), size=1000), 84)

            #     # y16, y50, y84 = samples.confidence('z', 0.16), samples.confidence('z', 0.50), samples.confidence('z', 0.84)
            #     # ax.set_title(r'$z =' + f'{y50:.2f}' + '^{+' + f'{y84-y50:.2f}' + '}_{-' + f'{y50-y16:.2f}' + '}$', fontsize=8)
            #     wav0 = np.logspace(-1, 1.5, 3000)
            #     y = np.interp(wav0, sed['wav_rest'][0]*(1+z0), np.median(sed['f_nu'],axis=0))
            #     ax_sed.plot(wav0, y, color='royalblue', linewidth=1)

            #     sed_samples = np.zeros((len(sed['f_nu']), len(wav0)))
            #     for i in range(len(sed_samples)):
            #         z = sed['z'][i]
            #         sed_samples[i] = np.interp(wav0, sed['wav_rest'][0]*(1+z), sed['f_nu'][i])

            #     yl95, yl68, yu68, yu95 = np.percentile(sed_samples, [2.5,16,84,97.5], axis=0)
            #     ax_sed.fill_between(wav0, yl68, yu68, facecolor='royalblue', edgecolor='none', alpha=0.15, zorder=-999)
            #     ax_sed.fill_between(wav0, yl95, yu95, facecolor='royalblue', edgecolor='none', alpha=0.1, zorder=-1000)
            #     # t = Table({'wav_obs': wav0,'fnu_l95': yl95, 'fnu_l68': yl68, 'fnu_med':y, 'fnu_u68':yu68, 'fnu_u95':yu95})
            #     # t.write(f'{out_sed_dir}/{ID}_brisket_sed.txt', format='ascii', overwrite=True)


            if out_format == 'pdf':
                plt.savefig(outpath)
            elif out_format == 'png':
                assert dpi is not None, 'dpi must be specified for png output'
                plt.savefig(outpath, dpi=300)
            else:
                raise Exception("out_format must be 'pdf' or 'png'")

            plt.close()





if __name__ == '__main__':


    IDs_all = np.array([27244])
    catalog_filename = 'COSMOSWeb_master_v3.1.0_assoc_cold+hot_sersic_cgs_err-calib.fits'
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
            

    # figure out which tile each ID is in
    f = fits.getdata(os.path.join(catalog_path, catalog_filename))
    tiles_all = []
    for ID in IDs_all:
        fi = f[f['ID_SE++']==ID]
        tiles_all.append(fi['TILE'][0])
    tiles_all = np.array(tiles_all)

    # split IDs_all array into dict of arrays, one for each tile
    IDs = {}
    for tile in np.unique(tiles_all):
        IDs[tile] = IDs_all[tiles_all==tile]
    tiles = np.unique(tiles_all)

    
    main(IDs, 
         tiles,
         logo=logo, 
         catalog_path=catalog_path, 
         catalog_filename=catalog_filename, 
         catalog_shortname=catalog_shortname, 
         outdir=outdir, 
         out_format='png', dpi=600, 
         overwrite=True, 
         display_width=4*u.arcsec, 
         cutout_width=6*u.arcsec, 
         vmax_in=10,
         cmap='Greys', 
         lephare_spec_path = lephare_spec_path, 
        #  zphot_cat = zphot_cat, 
        #  zpdf_cat = zpdf_cat,
        #  eazy_run=None, 
        #  eazy_outdir=None, 
         verbose=False)