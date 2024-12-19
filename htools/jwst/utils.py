

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



