'''
This module contains the paths to various JWST datasets.
'''

from .. import host
from . import config

def get_cosmos_web_rgb_filepath():
    if host.hostname == 'patrick':
        return '/Users/hba423/fitsmap/RGB_60mas_tot_edited_gimp_compressed.pickle'
    elif host.hostname == 'candide':
        return '/n23data2/hakins/COSMOS-Web/RGB_60mas_tot_edited_gimp_compressed.pickle'

def get_cosmos_web_full_f444w_filepath():
    if host.hostname == 'patrick':
        return '/Users/hba423/fitsmap/data3/CW_f444w_60mas_tot_v8.fits'
    elif host.hostname == 'candide':
        return '/n23data2/hakins/exchg/segmap_60mas_tot.fits'

def get_cosmos_web_detec_filepath(tile):
    if host.hostname == 'patrick':
        filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics/detection_chi2pos_SWLW_{tile}.fits', 0
    elif host.hostname == 'candide':
        filepath, hdu_index = f'/n23data2/hakins/COSMOS-Web/detection_images/detection_chi2pos_SWLW_{tile}.fits', 0
    return filepath, hdu_index

def get_cosmos_web_segm_filepath(tile, catalog_version='v1.3'):
    if host.hostname == 'patrick':
        filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics/detection_chi2pos_SWLW_{tile}_segmap_{catalog_version}.fits', 0
    elif host.hostname == 'candide':
        filepath, hdu_index = f'/n23data2/hakins/COSMOS-Web/segmaps/{catalog_version}/detection_chi2pos_SWLW_{tile}_segmap_{catalog_version}.fits', 0
    return filepath, hdu_index

                
def get_filepath(field, band, ext, ps='30mas', tile=None):
    '''
    general helper function to load a given JWST (or other) image
    some fields are broken into tiles to save memory; for these, specify e.g. tile='A5'
    `field` must be one of 'cosmos-web', 'primer-cosmos', 'ceers' (for now)
    `band` must be in the above config
    '''

    if field not in config.all_fields:
        raise ValueError(f"{field} not understood")
    
    if field == 'cosmos-web':
        assert tile is not None, "For COSMOS-Web, you must specify a tile"
        if band in config.bands['cosmos-web']:
            return get_cosmos_web_filepath(band, ext, tile, ps)
        elif band in config.bands['cosmos']:
            return get_cosmos_filepath(band, ext, tile)
    
    elif field == 'primer-cosmos':
        return get_primer_cosmos_filepath(band, ext, ps)
    
    elif field == 'ceers':
        return get_ceers_filepath(band, ext, tile)
    
    elif field == 'cosmos':
        return get_cosmos_filepath(band, ext, tile)
        


def get_cosmos_web_filepath(band, ext, tile, ps='30mas'):
    if band not in config.bands['cosmos-web']:
        raise ValueError(f"{band} not available for COSMOS-Web")
    
    if ext not in ['sci','wht','err','mod']:
        raise ValueError(f"{ext} not understood")
    
    if ps not in ['30mas','60mas']:
        raise ValueError(f"{ps} not understood")
    
    if tile not in ['A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','B1','B2','B3','B4','B5','B6','B7','B8','B9','B10']:
        raise ValueError(f"{tile} not understood")

    if band=='f814w':
        if ext.startswith('sci'): ext = ext.replace('sci','drz')

    if host.hostname == 'patrick':
        if tile.startswith('B'):
            if band=='f814w':
                filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics/mosaic_cosmos_web_2024jan_{ps}_tile_{tile}_hst_acs_wfc_f814w_{ext}.fits', 0
            elif band=='f770w':
                if tile in ['B1','B2','B3','B4','B5','B6','B7','B8','B9']:
                    filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics/mosaic_miri_f770w_COSMOS-Web_{ps}_{tile}_v0_7_{ext}.fits', 1
                if tile in ['B10']:
                    filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics/mosaic_miri_f770w_COSMOS-Web_{ps}_{tile}_v0_6_{ext}.fits', 1
            elif band in ['f115w','f150w','f277w','f444w']:
                filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics/mosaic_nircam_{band}_COSMOS-Web_{ps}_{tile}_v0_8_{ext}.fits', 0

        elif tile.startswith('A'):
            if band=='f814w':
                filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics/mosaic_cosmos_web_2023apr_{ps}_tile_{tile}_hst_acs_wfc_f814w_{ext}.fits', 0
            elif band=='f770w':
                if tile in ['A1']:
                    filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics/mosaic_miri_f770w_COSMOS-Web_{ps}_{tile}_v0_7_{ext}.fits', 1           
                if tile in ['A2','A3','A4','A5','A6','A8']:
                    filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics/mosaic_miri_f770w_COSMOS-Web_{ps}_{tile}_v0_6_{ext}.fits', 1           
                if tile in ['A7','A9','A10']:
                    filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics/mosaic_miri_f770w_COSMOS-Web_{ps}_{tile}_v0_3_{ext}.fits', 0           
            elif band in ['f115w','f150w','f277w','f444w']:
                filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics/mosaic_nircam_{band}_COSMOS-Web_{ps}_{tile}_v0_8_{ext}.fits', 0
        
        if ext == 'mod':
            filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_model/v3.1.0/CheckImg_COSMOSWeb__model_mosaic_nircam_{band}_COSMOS-Web_60mas_{tile}_v0_8_sci_1.fits', 0

    elif host.hostname == 'candide':
        if tile.startswith('B'):
            if band=='f814w':
                filepath, hdu_index = f'/n23data2/cosmosweb/COSMOS-Web_Jan24/ACS/mosaic_cosmos_web_2024jan_{ps}_tile_{tile}_hst_acs_wfc_f814w_{ext}.fits', 0 
            elif band=='f770w':
                if tile in ['B1','B2','B3','B4','B5','B6','B7','B8','B9']:
                    filepath, hdu_index = f'/n23data2/hakins/COSMOS-Web/mosaics/mosaic_miri_f770w_COSMOS-Web_{ps}_{tile}_v0_7_{ext}.fits', 1
                if tile in ['B10']:
                    filepath, hdu_index = f'/n23data2/hakins/COSMOS-Web/mosaics/mosaic_miri_f770w_COSMOS-Web_{ps}_{tile}_v0_6_{ext}.fits', 1
            elif band in ['f115w','f150w','f277w','f444w']:
                filepath, hdu_index = f'/n23data2/hakins/COSMOS-Web/mosaics/mosaic_nircam_{band}_COSMOS-Web_{ps}_{tile}_v0_8_{ext}.fits', 0
        elif tile.startswith('A'):
            if band=='f814w':
                filepath, hdu_index = f'/n23data2/cosmosweb/COSMOS-Web_Apr23/ACS/mosaic_cosmos_web_2023apr_{ps}_tile_{tile}_hst_acs_wfc_f814w_{ext}.fits', 0
            elif band=='f770w':
                if tile in ['A1']:
                    filepath, hdu_index = f'/n23data2/hakins/COSMOS-Web/mosaics/mosaic_miri_f770w_COSMOS-Web_{ps}_{tile}_v0_7_{ext}.fits', 1 
                if tile in ['A2','A3','A4','A5','A6','A8']:
                    filepath, hdu_index = f'/n23data2/hakins/COSMOS-Web/mosaics/mosaic_miri_f770w_COSMOS-Web_{ps}_{tile}_v0_6_{ext}.fits', 1
                if tile in ['A7','A9','A10']:
                    filepath, hdu_index = f'/n23data2/hakins/COSMOS-Web/mosaics/mosaic_miri_f770w_COSMOS-Web_{ps}_{tile}_v0_3_{ext}.fits', 0           
            elif band in ['f115w','f150w','f277w','f444w']:
                filepath, hdu_index = f'/n23data2/hakins/COSMOS-Web/mosaics/mosaic_nircam_{band}_COSMOS-Web_{ps}_{tile}_v0_8_{ext}.fits', 0
        
        if ext == 'mod':
            filepath, hdu_index = f'/n17data/shuntov/COSMOS-Web/CheckImages/JAN24-{tile}_v3.1.0-ASC/CheckImg_COSMOSWeb__model_mosaic_nircam_{band}_COSMOS-Web_60mas_{tile}_v0_8_sci_1.fits', 0
    
    return filepath, hdu_index

def get_ceers_filepath(band, ext, tile):
    if band not in config.bands['ceers']:
        raise ValueError(f"{band} not available for CEERS")

    if tile not in ['miri1', 'miri2', 'miri3', 'miri6', 'miri7', 'miri9', 'nircam5']:
        raise ValueError(f"{tile} not understood")
    
    if hostname.host == 'patrick':
        if band in ['f606w','f814w','f105w','f125w','f140w','f160w']:
            if ext=='wht': raise Exception("WHT map don't exist for CEERS HST")
            if ext=='sci_pfMatched': raise Exception("I don't yet have PSF matched mosaics for CEERS")
            
            if band in ['f606w','f814w']:
                if ext=='sci': filepath, hdu_index = f'/V/simmons/ceers/mosaics/{tile}/egs_all_acs_wfc_{band}_030mas_v1.9_{tile}_mbkgsub1.fits', 0
                if ext=='err': filepath, hdu_index = f'/V/simmons/ceers/mosaics/{tile}/egs_all_acs_wfc_{band}_030mas_v1.9_{tile}_rms.fits', 0
            elif band in ['f105w','f125w','f140w','f160w']:
                if ext=='sci': filepath, hdu_index = f'/V/simmons/ceers/mosaics/{tile}/egs_all_wfc3_ir_{band}_030mas_v1.9_{tile}_mbkgsub1.fits', 0
                if ext=='err': filepath, hdu_index = f'/V/simmons/ceers/mosaics/{tile}/egs_all_wfc3_ir_{band}_030mas_v1.9.1_{tile}_rms.fits', 0
        elif band in ['f115w','f150w','f200w','f277w','f356w','f410m','f444w']:
            if ext=='sci': hdu_index=1
            if ext=='err': hdu_index=3
            if ext=='wht': hdu_index=5
            if ext=='sci_psfMatched': raise Exception("I don't yet have PSF matched mosaics for CEERS")            
            filepath = f'/V/simmons/ceers/mosaics/{tile}/ceers_{tile}_{band}_v1_mbkgsub1.fits'
        elif band in ['f560w','f770w']:
            if ext=='sci': hdu_index=1
            if ext=='err': hdu_index=2
            if ext=='wht': hdu_index=4
            if ext=='sci_pfMatched': raise Exception("I don't yet have PSF matched mosaics for CEERS")            
            filepath = f'/V/simmons/ceers/mosaics/{tile}/ceers_{tile}_{band}_i2d.fits'
    else:
        raise Exception("I don't yet have CEERS mosaics on candide")

    return filepath, hdu_index
    

def get_primer_cosmos_filepath(band, ext, ps='30mas'):
    if band not in config.bands['primer-cosmos']:
        raise ValueError(f"{band} not available for PRIMER-COSMOS")
    
    if ext not in ['sci','wht','err']:
        raise ValueError(f"{ext} not understood")
     
    if ps not in ['30mas','60mas']:
        raise ValueError(f"{ps} not understood")

    if hostname.host == 'patrick':
        if band in ['f606w','f814w']:
            if ext.startswith('sci'): ext = ext.replace('sci','drz')
            filepath, hdu_index = f'/V/simmons/primer/mosaics/mosaic_cosmos_primer_{ps}_acs_wfc_{band}_{ext}.fits', 0
        elif band in ['f125w','f160w']:
            if ext.startswith('sci'): ext = ext.replace('sci','drz')
            filepath, hdu_index = f'/V/simmons/primer/mosaics/mosaic_cosmos_primer_{ps}_wfc3_ir_{band}_{ext}.fits', 1
        elif band in ['f770w','f1800w']:
            filepath, hdu_index = f'/V/simmons/primer/mosaics/mosaic_miri_{band}_PRIMER-COSMOS_{ps}_{ext}_reproj.fits', 0
        elif band in ['f090w','f115w','f150w','f200w','f277w','f356w','f410m','f444w']:
            filepath, hdu_index = f'/V/simmons/primer/mosaics/mosaic_nircam_{band}_PRIMER-COSMOS_{ps}_{ext}.fits', 0
    elif hostname.host == 'candide':
        raise Exception("I don't yet have PRIMER-COSMOS mosaics on candide")
    
    return filepath, hdu_index

            
def get_cosmos_filepath(band, ext, tile=None):
    '''Filepaths for COSMOS ground-based and IRAC data'''

    if band not in config.bands['cosmos']:
        if band.startswith('f') and band.endswith('w'):
            raise ValueError(f"{band} not available for COSMOS. Did you mean to use COSMOS-Web?")
        raise ValueError(f"{band} not available for COSMOS")

    if ext not in ['sci','wht']:
        raise ValueError(f"{ext} not understood")
    if host.hostname == 'patrick':
        if band.startswith('IRAC'):
            if tile is not None:
                if ext == 'sci':
                    filepath, hdu_index = f"/V/simmons/cosmos-web/mosaics_ground/irac.{band.split('IRAC')[-1]}.mosaic.Fconv_resamp015_zp-28.09_{tile}.fits", 0
                else:
                    filepath, hdu_index = f"/V/simmons/cosmos-web/mosaics_ground/irac.{band.split('IRAC')[-1]}.mosaic.Fconv_resamp015_weight_zp-28.09_{tile}.fits", 0
            else:
                # load full area IRAC mosaics
                raise NotImplementedError("IRAC mosaics for the full COSMOS area not implemented ")

        if band in ['NB118','Y','J','H','Ks']:
            if tile is not None:
                if ext == 'sci': 
                    filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/UVISTA_{band}_12_01_24_allpaw_skysub_015_dr6_rc_v1_zp-28.09_{tile}.fits', 0
                else: 
                    filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/UVISTA_{band}_12_01_24_allpaw_skysub_015_dr6_rc_v1.weight_zp-28.09_{tile}.fits', 0
            else:
                # load full area UVISTA mosaics
                raise NotImplementedError("UVISTA mosaics for the full COSMOS area not implemented ")
            
        if band in ['g','r','i','z','y']:
            if tile is not None:
                if tile in ['A4', 'A5', 'A9', 'A10']:
                    if band=='g': filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/cutout-HSC-G-9813-pdr3_dud_rev-230413-130357_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='r': filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/cutout-HSC-R-9813-pdr3_dud_rev-230413-130346_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='i': filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/cutout-HSC-I-9813-pdr3_dud_rev-230413-130351_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='z': filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/cutout-HSC-Z-9813-pdr3_dud_rev-230413-130355_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='y': filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/cutout-HSC-Y-9813-pdr3_dud_rev-230413-130357_{ext}_zp-28.09_{tile}.fits', 0
                elif tile in ['A1', 'A2', 'A3', 'A8', 'A7', 'A6']:
                    if band=='g': filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/cutout-HSC-G-9813-pdr3_dud_rev-230412-135737_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='r': filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/cutout-HSC-R-9813-pdr3_dud_rev-230413-121613_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='i': filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/cutout-HSC-I-9813-pdr3_dud_rev-230413-121625_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='z': filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/cutout-HSC-Z-9813-pdr3_dud_rev-230413-121629_{ext}_zp-28.09_{tile}.fits', 0
                    if band=='y': filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/cutout-HSC-Y-9813-pdr3_dud_rev-230413-121631_{ext}_zp-28.09_{tile}.fits', 0
                else:
                    if band=='g': filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/{tile}--cutout-HSC-G-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
                    if band=='r': filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/{tile}--cutout-HSC-R-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
                    if band=='i': filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/{tile}--cutout-HSC-I-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
                    if band=='z': filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/{tile}--cutout-HSC-Z-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
                    if band=='y': filepath, hdu_index = f'/V/simmons/cosmos-web/mosaics_ground/{tile}--cutout-HSC-Y-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
            else:
                raise NotImplementedError("HSC mosaics for the full COSMOS area not implemented ")
                # load full area HSC mosaics
    
    elif host.hostname == 'candide':
        if band.startswith('IRAC'):
            if tile is not None:
                if ext == 'sci':
                    filepath, hdu_index = f"/n23data2/hakins/exchg/grounddata/irac.{band.split('IRAC')[-1]}.mosaic.Fconv_resamp015_zp-28.09_{tile}.fits", 0
                else:
                    filepath, hdu_index = f"/n23data2/hakins/exchg/grounddata/irac.{band.split('IRAC')[-1]}.mosaic.Fconv_resamp015_weight_zp-28.09_{tile}.fits", 0
            else:
                # load full area IRAC mosaics
                raise NotImplementedError("IRAC mosaics for the full COSMOS area not implemented ")

        elif band in ['NB118','Y','J','H','Ks']:
            if tile is not None:
                if ext == 'sci': 
                    filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/UVISTA_{band}_12_01_24_allpaw_skysub_015_dr6_rc_v1_zp-28.09_{tile}.fits', 0
                else: 
                    filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/UVISTA_{band}_12_01_24_allpaw_skysub_015_dr6_rc_v1.weight_zp-28.09_{tile}.fits', 0
            else:
                # load full area UVISTA mosaics
                raise NotImplementedError("UVISTA mosaics for the full COSMOS area not implemented ")
            
        elif band in ['g','r','i','z','y']:
            if tile is not None:
                if tile in ['A4', 'A5', 'A9', 'A10']:
                    if band=='g': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-HSC-G-9813-pdr3_dud_rev-230413-130357_{ext}_zp-28.09_{tile}.fits', 0
                    elif band=='r': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-HSC-R-9813-pdr3_dud_rev-230413-130346_{ext}_zp-28.09_{tile}.fits', 0
                    elif band=='i': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-HSC-I-9813-pdr3_dud_rev-230413-130351_{ext}_zp-28.09_{tile}.fits', 0
                    elif band=='z': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-HSC-Z-9813-pdr3_dud_rev-230413-130355_{ext}_zp-28.09_{tile}.fits', 0
                    elif band=='y': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-HSC-Y-9813-pdr3_dud_rev-230413-130357_{ext}_zp-28.09_{tile}.fits', 0
                elif tile in ['A1', 'A2', 'A3', 'A8', 'A7', 'A6']:
                    if band=='g': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-HSC-G-9813-pdr3_dud_rev-230412-135737_{ext}_zp-28.09_{tile}.fits', 0
                    elif band=='r': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-HSC-R-9813-pdr3_dud_rev-230413-121613_{ext}_zp-28.09_{tile}.fits', 0
                    elif band=='i': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-HSC-I-9813-pdr3_dud_rev-230413-121625_{ext}_zp-28.09_{tile}.fits', 0
                    elif band=='z': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-HSC-Z-9813-pdr3_dud_rev-230413-121629_{ext}_zp-28.09_{tile}.fits', 0
                    elif band=='y': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-HSC-Y-9813-pdr3_dud_rev-230413-121631_{ext}_zp-28.09_{tile}.fits', 0
                    else: raise ValueError(f"{band} not understood")
                else:
                    if band=='g': filepath, hdu_index = f'/n23data2/hakins/exchg/hscdata/{tile}--cutout-HSC-G-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
                    elif band=='r': filepath, hdu_index = f'/n23data2/hakins/exchg/hscdata/{tile}--cutout-HSC-R-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
                    elif band=='i': filepath, hdu_index = f'/n23data2/hakins/exchg/hscdata/{tile}--cutout-HSC-I-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
                    elif band=='z': filepath, hdu_index = f'/n23data2/hakins/exchg/hscdata/{tile}--cutout-HSC-Z-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
                    elif band=='y': filepath, hdu_index = f'/n23data2/hakins/exchg/hscdata/{tile}--cutout-HSC-Y-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
                    else: raise ValueError(f"{band} not understood")
            else:
                raise NotImplementedError("HSC mosaics for the full COSMOS area not implemented ")
                # load full area HSC mosaics

        elif band in ['NB0816', 'NB0921', 'NB1010']:
            if tile in ['A4', 'A5', 'A9', 'A10']:
                if band == 'NB0816': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-NB0816-9813-pdr3_dud_rev-230413-130357_{ext}_zp-28.09_{tile}.fits', 0
                elif band == 'NB0921': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-NB0921-9813-pdr3_dud_rev-230413-130357_{ext}_zp-28.09_{tile}.fits', 0
                elif band == 'NB1010': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-NB1010-9813-pdr3_dud_rev-230413-130357_{ext}_zp-28.09_{tile}.fits', 0
                else: raise ValueError(f"{band} not understood")
            elif tile in ['A1', 'A2', 'A3', 'A8', 'A7', 'A6']:
                if band == 'NB0816': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-NB0816-9813-pdr3_dud_rev-230413-121622_{ext}_zp-28.09_{tile}.fits', 0
                elif band == 'NB0921': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-NB0921-9813-pdr3_dud_rev-230413-121626_{ext}_zp-28.09_{tile}.fits', 0
                elif band == 'NB1010': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/cutout-NB1010-9813-pdr3_dud_rev-230413-121845_{ext}_zp-28.09_{tile}.fits', 0
                else: raise ValueError(f"{band} not understood")
            else:
                if band == 'NB0816': filepath, hdu_index = f'/n23data2/hakins/exchg/hscdata/{tile}--cutout-NB0816-9813-pdr3_dud_rev_{ext}_zp-28.09.fits', 0
                elif band == 'NB0921': filepath, hdu_index = f'/n23data2/hakins/exchg/hscdata/{tile}--cutout-NB0921-9813-pdr3_dud_rev_{ext}_zp-28.09.fits.fits', 0
                elif band == 'NB1010': filepath, hdu_index = f'/n23data2/hakins/exchg/hscdata/{tile}--cutout-NB1010-9813-pdr3_dud_rev_{ext}_zp-28.09.fits.fits', 0
                else: raise ValueError(f"{band} not understood")
                
        else:
            if ext == 'sci':
                if band == 'IB427': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L427_20-09-29a_cosmos_zp-28.09_{tile}.fits'
                elif band == 'IB464': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L464_20-09-29a_cosmos_zp-28.09_{tile}.fits'
                elif band == 'IA484': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L484_20-09-29a_cosmos_zp-28.09_{tile}.fits'
                elif band == 'IB505': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L505_20-09-29a_cosmos_zp-28.09_{tile}.fits'
                elif band == 'IA527': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L527_20-09-29a_cosmos_zp-28.09_{tile}.fits'
                elif band == 'IB574': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L574_20-09-29a_cosmos_zp-28.09_{tile}.fits'
                elif band == 'IA624': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L624_20-09-29a_cosmos_zp-28.09_{tile}.fits'
                elif band == 'IA679': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L679_20-09-29a_cosmos_zp-28.09_{tile}.fits'
                elif band == 'IB709': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L709_20-09-29a_cosmos_zp-28.09_{tile}.fits'
                elif band == 'NB711': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L711_20-09-29a_cosmos_zp-28.09_{tile}.fits'
                elif band == 'IA738': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L738_20-09-29a_cosmos_zp-28.09_{tile}.fits'
                elif band == 'IA767': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L767_20-09-29a_cosmos_zp-28.09_{tile}.fits'
                elif band == 'NB816': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L816_20-09-29a_cosmos_zp-28.09_{tile}.fits'
                elif band == 'IB827': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L827_20-09-29a_cosmos_zp-28.09_{tile}.fits'
                else: raise ValueError(f"{band} not understood")
            else:
                if band == 'IB427': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L427_20-09-29a_cosmos.weight_zp-28.09_{tile}.fits'
                elif band == 'IB464': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L464_20-09-29a_cosmos.weight_zp-28.09_{tile}.fits'
                elif band == 'IA484': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L484_20-09-29a_cosmos.weight_zp-28.09_{tile}.fits'
                elif band == 'IB505': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L505_20-09-29a_cosmos.weight_zp-28.09_{tile}.fits'
                elif band == 'IA527': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L527_20-09-29a_cosmos.weight_zp-28.09_{tile}.fits'
                elif band == 'IB574': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L574_20-09-29a_cosmos.weight_zp-28.09_{tile}.fits'
                elif band == 'IA624': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L624_20-09-29a_cosmos.weight_zp-28.09_{tile}.fits'
                elif band == 'IA679': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L679_20-09-29a_cosmos.weight_zp-28.09_{tile}.fits'
                elif band == 'IB709': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L709_20-09-29a_cosmos.weight_zp-28.09_{tile}.fits'
                elif band == 'NB711': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L711_20-09-29a_cosmos.weight_zp-28.09_{tile}.fits'
                elif band == 'IA738': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L738_20-09-29a_cosmos.weight_zp-28.09_{tile}.fits'
                elif band == 'IA767': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L767_20-09-29a_cosmos.weight_zp-28.09_{tile}.fits'
                elif band == 'NB816': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L816_20-09-29a_cosmos.weight_zp-28.09_{tile}.fits'
                elif band == 'IB827': filepath, hdu_index = f'/n23data2/hakins/exchg/grounddata/SPC_L827_20-09-29a_cosmos.weight_zp-28.09_{tile}.fits'
                else: raise ValueError(f"{band} not understood")

    return filepath, hdu_index
