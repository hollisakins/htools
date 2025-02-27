import matplotlib as mpl
import numpy as np
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from photutils.aperture import CircularAperture, RectangularAperture, CircularAnnulus, aperture_photometry

# imagesRGB = {"R": ['../multiple_bands/HSC/J114444.8+001347.0-HSC-I-pdr3_wide.fits[1]'], \
#              "G": ['../multiple_bands/HSC/J114444.8+001347.0-HSC-R-pdr3_wide.fits[1]'], \
#              "B": ['../multiple_bands/HSC/J114444.8+001347.0-HSC-G-pdr3_wide.fits[1]']}
# noiselums = {'R': 0.5, 'G': 0.5, 'B': 0.5}
      
# trilogy.Trilogy(infile = None, samplesize = 20000, stampsize = 20000, maxstampsize = 20000, \
#                 deletetests = 1, deletefilters = 1, testfirst = 0, showwith = "PIL", \
#                 mode = 'RGB', imagesorder = 'RGB', imagesRGB = imagesRGB, noiselums = noiselums, images = None, \
#                 outname = 'multiple_band_image', satpercent = 0.0009, noiselum = 0.5, noisesig = 50, \
#                 noisesig0 = 10, correctbias = 0, colorsatfac = 1, combine = 'sum', show = True).run()


def gen_rgb_image(input_dict, noisesig=1, noiselum=0.1, satpercent=0.5, save=False):
    # from scipy.optimize import golden
    # def da(k):
    #     a1 = k * (x1 - x0) + 1
    #     a2 = k * (x2 - x0) + 1
    #     a1n = a1**n
    #     a1n = np.abs(a1n)  # Don't want the solutions where a1 & a2 are both negative!
    #     da1 = a1n - a2
    #     k = np.abs(k)
    #     if k == 0:
    #         return da(1e-10)
    #     else:
    #         da1 = da1 / k  # To avoid solution k = 0!
    #     return abs(da1)
    # def imscale2(data, levels, y1):
    #     # x0, x1, x2  YIELD  0, y1, 1,  RESPECTIVELY
    #     # y1 = noiselum
    #     global n, x0, x1, x2  # So that golden can use them
    #     x0, x1, x2 = levels
    #     if y1 == 0.5:
    #         k = (x2 - 2 * x1 + x0) / float(x1 - x0) ** 2
    #     else:
    #         n = 1 / y1
    #         k = np.abs(golden(da))
    #     r1 = np.log10( k * (x2 - x0) + 1)
    #     v = np.ravel(data)
    #     v = clip2(v, 0, None)
    #     d = k * (v - x0) + 1
    #     d = clip2(d, 1e-30, None)
    #     z = np.log10(d) / r1
    #     z = np.clip(z, 0, 1)
    #     z.shape = data.shape
    #     z = z * 255
    #     z = z.astype(np.uint8)
    #     return z
    # def clip2(m, m_min=None, m_max=None):
    #     # nanmin and nanmax important to ignore nan values
    #     # otherwise you'll get all 0's
    #     if m_min == None:
    #         m_min = np.nanmin(m)
    #     if m_max == None:
    #         m_max = np.nanmax(m)
    #     return np.clip(m, m_min, m_max)
    # # # PREVIOUSLY in colorimage.py
    # def set_levels(data, pp, stripneg=False, sortedalready=False):
    #     if sortedalready:
    #         vs = data
    #     else:
    #         print('sorting...')
    #         vs = np.sort(data.flat)
    #     if stripneg:  
    #         i = np.searchsorted(vs, 0)
    #         vs = vs[i+1:]
    #     else:  # Clip negative values to zero
    #         vs = clip2(vs, 0, None)
    #     ii = np.array(pp) * len(vs)
    #     ii = ii.astype(int)
    #     ii = np.clip(ii, 0, len(vs)-1)
    #     levels = vs.take(ii)
    #     return levels
    # def determine_scaling(data, unsatpercent, noisesig=1, correctbias=True, noisefloorsig=2):
        # """Determines data values (x0,x1,x2) which will be scaled to (0,noiselum,1)"""
        # Robust mean & standard deviation
        # datasorted = data + 0
        # datasorted[np.isnan(datasorted)]=0  # set all nan values to zero
        # datasorted = np.sort(datasorted.flat)
        # if datasorted[0] == datasorted[-1]:  # data is all one value
        #     levels = 0, 1, 100  # whatever
        # else:
        #     data_mean, data_median, data_stddev = sigma_clipped_stats(datasorted)
        #     m = data_mean
        #     r = data_stddev
        #     #print('%g +/- %g' % (m, r))

        #     if correctbias:
        #         x0 = m - noisefloorsig * r
        #     else:
        #         x0 = 0
        #     x1 = m + noisesig * r
        #     x2 = set_levels(datasorted, np.array([unsatpercent]), sortedalready=True)[0]
        #     levels = x0, x1, x2
        # return levels

    # filters = list(input_dict.keys())
    # filter_colors = {f:input_dict[f]['colors'] for f in filters}
    # image_data_dict = {f:input_dict[f]['data'] for f in filters}

    # rgb_lum_sum = np.zeros(3)
    # for i, filt in enumerate(filters):
    #     rgb_lum_sum += np.array(filter_colors[filt])

    # unsatpercent = 1 - 0.01 * satpercent
    # # # Scaling parameters
    # # noiselum   = 0.2
    # # satpercent = 0.
    # # noisesig = 4
    # # correctbias = False  ## yes because need to dip below 0 by noisefloorsig-sigma
    # # noisefloorsig = 1  #  set black to e.g., 2-sigma below
    # # upscale = 0.7*np.sqrt(stack_snr)

    # scaled_images = {}
    # levels_all = {}
    # for filt in filters:
    #     levels = determine_scaling(image_data_dict[filt].ravel(), unsatpercent, noisesig, correctbias, noisefloorsig)
    #     scaled = imscale2(image_data_dict[filt], levels, noiselum)
    #     levels_all[filt] = levels
    #     scaled_images[filt] = scaled

    # rgb_total = 0
    # for filt in filters:
    #     rgb = r, g, b = filter_colors[filt][:, np.newaxis, np.newaxis] * scaled_images[filt]
    #     rgb_total += rgb

    # r, g, b = rgb_average = rgb_total / rgb_lum_sum[:, np.newaxis, np.newaxis]
    # if upscale is not None:
    #     r = np.where(r*upscale>255,255,r*upscale)
    #     g = np.where(g*upscale>255,255,g*upscale)
    #     b = np.where(b*upscale>255,255,b*upscale)
        

    filters = list(input_dict.keys())
    filter_colors = {f:input_dict[f]['colors'] for f in filters}
    image_data_dict = {f:input_dict[f]['data'] for f in filters}
    rgb_lum_sum = np.zeros(3)
    for i, filt in enumerate(filters):
        rgb_lum_sum += np.array(filter_colors[filt])
    unsatpercent = 1 - 0.01 * satpercent

    from astropy.stats import sigma_clipped_stats
    rgb_total = 0
    for filt in filters:
        rgb = r, g, b = filter_colors[filt][:, np.newaxis, np.newaxis] * image_data_dict[filt]
        rgb_total += rgb
    r, g, b = rgb_average = rgb_total / rgb_lum_sum[:, np.newaxis, np.newaxis]

    stds = [sigma_clipped_stats(i)[2] for i in [r,g,b]]
    blackpoint = noisesig*np.max(stds)
    whitepoint = np.nanpercentile([r,g,b], 100*unsatpercent)
    r = (np.log10(r)-np.log10(blackpoint))/(np.log10(whitepoint)-np.log10(blackpoint))
    r = r * (255*(1-noiselum)) + 255*noiselum
    g = (np.log10(g)-np.log10(blackpoint))/(np.log10(whitepoint)-np.log10(blackpoint))
    g = g * (255*(1-noiselum)) + 255*noiselum
    b = (np.log10(b)-np.log10(blackpoint))/(np.log10(whitepoint)-np.log10(blackpoint))
    b = b * (255*(1-noiselum)) + 255*noiselum
    r = np.where(r>255,255,r)
    g = np.where(g>255,255,g)
    b = np.where(b>255,255,b)
    r = np.where(np.isnan(r)|(r<0),0,r)
    g = np.where(np.isnan(g)|(g<0),0,g)
    b = np.where(np.isnan(b)|(b<0),0,b)

    imrgb = np.array([r, g, b]).transpose((1,2,0)).astype(np.uint8)
    if save:
        import matplotlib as mpl
        mpl.image.imsave(save, imrgb)
    return imrgb



####################################################################################################
#################################### Helpful plotting functions ####################################
####################################################################################################

def plotBeam(ax, hdu, xy=None, **kwargs):
    # kwargs['fc'] = 'w'
    # kwargs['ec'] = 'k'
    # kwargs['alpha'] = 1
    # kwargs['linewidth'] = 0.5
    hdr = hdu.header
    
    Bmaj = dict(hdr)['BMAJ']*60*60
    Bmin = dict(hdr)['BMIN']*60*60
    Bpa = dict(hdr)['BPA']

    if xy==None:
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        xrange = xlim[1]-xlim[0]
        yrange = ylim[1]-ylim[0]
        xy = (xlim[0]+xrange*0.1, ylim[0]+yrange*0.1)
    
    e = mpl.patches.Ellipse(xy, height=Bmaj, width=Bmin, angle=-Bpa, **kwargs)
    ax.add_patch(e)
    
# setattr(mpl.axes.Axes, "plotBeam", plotBeam)


####################################################################################################
###################################### ALMA Image (Map) Class ######################################
####################################################################################################

splatalogue = {'CII158' : 1900.53690000, 
               'OIII88' : 3393.00624400,
               'OIII52' : 5785.87958900,
               'NII122' : 2459.38010085}


class fitsimage:
    def __init__(self, filepath, telescope=None, **kwargs):# crop=False, dx=120, dy=180, x0=None, y0=None):
        self.filepath = filepath
        self.kwargs = kwargs

        assert filepath.endswith('.fits')
        with fits.open(filepath) as h:
            self.header = h[0].header
            self.data = h[0].data


        if np.ndim(self.data) == 4: # ALMA data cubes & collapsed cubes i.e. continuum maps
            self.data = self.data[0] # remove stokes channel
            self.nchans = np.shape(self.data)[0]
            self.telescope = 'ALMA'
            
            with fits.open(filepath.replace('.fits','_pb.fits')) as h:
                self.pb = h[0].data[0]

            if self.nchans==1: # continuum maps have 1 velocity channel
                self.type = 'continuum'
                self.data = self.data[0]
                self.pb = self.pb[0]
                self.size = np.shape(self.data)[0]
            else: # cubes have more than one velocity channel
                self.type = 'cube'
                self.size= np.shape(self.data)[1]
        if np.ndim(self.data) == 3: # misc data... e.g. SCUBA-2 map
            self.nchans = np.shape(self.data)[0]
            self.telescope = '?'
            if self.nchans==1: # continuum maps have 1 velocity channel
                self.type = 'continuum'
                self.data = self.data[0]
                self.size = np.shape(self.data)[0]
            else: # cubes have more than one velocity channel
                self.type = 'cube'
                self.size= np.shape(self.data)[1]

        elif np.ndim(self.data) == 2: # other data, e.g. JWST or HST images
            self.telescope = '?'
            self.type = 'continuum'
            self.nchans = 1
            self.size = np.shape(self.data)[0]

                
        if telescope is not None:
            self.telescope = telescope
                        
        self._x0 = self.size//2
        self._y0 = self.size//2

            
        # if crop:
        #     if img_code.startswith('cube'):
        #         self.data = self.data[:,(y0-dy):(y0+dy),(x0-dx):(x0+dx)]
        #     elif img_code.startswith('cont'):
        #         self.data = self.data[(y0-dy):(y0+dy),(x0-dx):(x0+dx)]
        #     #self.std_map = self.std_map[(y0-dy):(y0+dy),(x0-dx):(x0+dx)]
        #     #self.error = self.error[(y0-dy):(y0+dy),(x0-dx):(x0+dx)]
        #     self.x0 = dx
        #     self.y0 = dy
           
        ### construct error map
        if self.telescope=='ALMA':
            if self.type == 'cube':
                mean, median, std = sigma_clipped_stats(self.data*self.pb, sigma=3, mask_value=np.nan)
                self.std = std
                self.std_map = self.std/np.median(self.pb, axis=0)
                self.error = self.std_map * np.sqrt(self.NPixPerBeam)
            elif self.type == 'continuum':
                mean, median, std = sigma_clipped_stats(self.data*self.pb, sigma=3)
                self.std = std
                self.std_map = self.std/self.pb
                self.error = self.std_map * np.sqrt(self.NPixPerBeam)
        else:
            if self.type == 'cube':
                mean, median, std = sigma_clipped_stats(self.data, sigma=3, mask_value=np.nan)
                self.std = std
                self.std_map = np.zeros(shape=np.shape(self.data)[1:]) + self.std
            elif self.type == 'continuum':
                mean, median, std = sigma_clipped_stats(self.data, sigma=3)
                self.std = std
                self.std_map = np.zeros(shape=np.shape(self.data)) + self.std
  

        
        # ### for line mfs maps, convert to moment 0 map
        # if obstype=='linemfs':
        #     if line=='OIII':
        #         nu_min = 416.6 # GHz
        #         nu_max = 417.8 
        #         nu_rest = 417.2 
        #     elif line=='CII':
        #         nu_min = 233.4 # GHz
        #         nu_max = 234.0 
        #         nu_rest = 233.7
                
        #     from astropy.constants import c
        #     c = c.to(u.m/u.s).value
        #     v_max = c*(nu_rest - nu_min)/nu_rest
        #     v_min = c*(nu_rest - nu_max)/nu_rest
        #     self.delta_v = (v_max - v_min)/1000
        #     self.data *= self.delta_v
        #     self.residual *= self.delta_v
        #     self.std *= self.delta_v
        #     self.std_map *= self.delta_v
        #     self.error *= self.delta_v
        
        ### convert from Jy/beam to Jy/pix
        if self.telescope=='ALMA':
            self.data /= self.NPixPerBeam
            self.std /= self.NPixPerBeam
            self.std_map /= self.NPixPerBeam
            self.error /= self.NPixPerBeam
        
    @property
    def x0(self):
        return self._x0

    @x0.setter
    def x0(self, value):
        self._x0 = value

    @property
    def y0(self):
        return self._y0

    @y0.setter
    def y0(self, value):
        self._y0 = value

    @property
    def extent(self):
        '''Pass to matplotlib imshow or contour functions to set image extent in units of arcsec from center'''
        left = self.x0*self.cell
        right = -left
        bottom = -self.y0*self.cell
        top = -bottom
        return (left,right,bottom,top)

    @property
    def dists(self):
        x, y = np.arange(0, self.size, 1), np.arange(0, self.size, 1)
        x, y = np.meshgrid(x,y)
        return np.sqrt((x-self.x0)**2 + (y-self.y0)**2)*self.cell
    
    @property
    def cell(self):
        assert ('RA' in self.header['CTYPE1']) or ('DEC' in self.header['CTYPE1'])
        try:
            return np.abs(self.header['CDELT1']*60*60)
        except: 
            from astropy.wcs.utils import proj_plane_pixel_scales
            wcs = WCS(self.header, naxis=2)
            return proj_plane_pixel_scales(wcs)[0]*3600



    
    @property
    def BeamArea(self):
        assert self.telescope == 'ALMA', 'BeamArea is only defined for ALMA images'
        Bmaj = self.header['BMAJ']*60*60
        Bmin = self.header['BMIN']*60*60
        return np.pi/(4*np.log(2))*Bmaj*Bmin

    @property
    def NPixPerBeam(self):
        assert self.telescope == 'ALMA', 'NPixPerBeam is only defined for ALMA images'
        return self.BeamArea/(self.cell**2)

    @property
    def restfreq(self):
        # check if z, line, or restfreq provided in kwargs
        assert 'z' in self.kwargs, 'please provide redshift (e.g. `z=6`) in kwargs for initial call'
        self.z = self.kwargs['z']
        assert 'line' in self.kwargs, "please provide line (e.g. `line='CII158'`) in kwargs for initial call"
        self.line = self.kwargs['line']

        self.line_freq_rest = splatalogue[self.line]
        self.line_freq_obs = self.line_freq_rest/(1+self.z)
        return self.line_freq_obs

    @restfreq.setter
    def restfreq(self, value):
        self.line_freq_obs = value
            

    @property
    def freq(self):
        assert self.type=='cube', 'Frequency array is only provided for datacubes'
        freq = (np.arange(0,self.nchans,1)-1)*self.header['CDELT3'] + self.header['CRVAL3']
        freq = freq / 1e9
        return freq
    
    @property
    def vel(self):
        assert self.type=='cube', 'Velocity array is only provided for datacubes'
        restfreq = self.restfreq
        vel = 2.998e5 * (restfreq - self.freq)/restfreq
        return vel

    def moment(self, i, vrange):
        assert type(vrange) == tuple, 'please provide vrange as a tuple of (min, max) e.g. (-500,500) [km/s]'
        d = copy(self.data)
        v = copy(self.vel)
        dv = np.abs(np.mean(v[1:]-v[:-1]))
        d = d[(v > vrange[0]) & (v < vrange[1])]
        v = v[(v > vrange[0]) & (v < vrange[1])]

        mom0 = np.sum(d,axis=0)*dv # Jy km/s
        mean, median, mom0_std = sigma_clipped_stats(mom0, sigma=3)
        self.mom0 = mom0
        self.mom0_std = mom0_std
        self.mom0_err = np.zeros(shape=np.shape(mom0)) + mom0_std * np.sqrt(self.NPixPerBeam)
        if i==0:
            return mom0

        if i >= 1:
            mom1 = np.sum(np.array([v[i]*d[i] for i in range(np.shape(d)[0])]),axis=0)/mom0*dv
            if i==1:
                return mom1

            mom2 = np.sqrt(np.sum(np.array([d[i]*(v[i]-mom1)**2 for i in range(np.shape(d)[0])]),axis=0)/mom0*dv)
            return mom2
        

    def RadialProfile(self, normalized=True, bins=np.arange(0.01, 4, 0.2), cutoff=True):
        '''Test'''
        from photutils.aperture import CircularAnnulus, aperture_photometry

        x0, y0 = self.x0, self.y0
            
        ### construct aperture
        apertures = [CircularAnnulus([self.x0,self.y0], r_in=r_in/self.cell, r_out=r_out/self.cell) for r_in,r_out in zip(bins[:-1],bins[1:])]
        self.apertures = apertures
            
        ### setup data and error
        data = copy(self.data)
        error = copy(self.error)
        
        ### perform photometry 
        phot_table = aperture_photometry(data, apertures, error=self.error, mask=np.isnan(data))
        self.phot_table = phot_table
        area = np.array([a.area for a in apertures])
        sb = np.array([self.phot_table[f'aperture_sum_{i}'][0] for i in range(len(bins)-1)])/area
        sb_err = np.array([self.phot_table[f'aperture_sum_err_{i}'][0] for i in range(len(bins)-1)])/area
        
        if normalized:
            sb_err /= np.max(sb)
            sb /= np.max(sb)
            
        if cutoff:
            if any(sb < 0):
                i = np.min(np.arange(len(sb))[sb < 0])
                sb[i:] = 0
            
        bc = 0.5*(bins[1:]+bins[:-1])
        
        return bc, sb, sb_err

    
    
    def BeamProfile(self, normalized=True, bins=np.arange(0.01, 16, 2.25)):
        from photutils.aperture import EllipticalAnnulus, aperture_photometry

        if self.plane == 'image':
            raise Exception('Beam profile is not implemented for image plane maps')
        elif self.plane == 'source':
            x0, y0 = self.y0-1.8, self.y0-1.8
            
        a = self.header['BMAJ']*60*60/self.cell # in pix
        b = self.header['BMIN']*60*60/self.cell # in pix
        theta = np.pi/2-self.header['BPA']*np.pi/180
        
        delta = np.mean(bins[1:]-bins[:-1])
        print(f'Beam size: {a*self.cell_kpc:.2f} x {b*self.cell_kpc:.2f} kpc, bin separation is {delta:.2f} kpc')

        blist = bins/self.cell_kpc
        alist = blist/b*a

        apertures = [EllipticalAnnulus((self.x0-1.8, self.y0-1.8), a_in=alist[i-1], a_out=alist[i], b_out=blist[i], theta=theta) for i in range(1,len(bins))]
        self.apertures = apertures

        error = self.error
        data = self.data

        data[np.isnan(data)] = 0 # set nan values to 0 (since astropy doesn't like them) but so they won't contribute to flux
        phot_table = aperture_photometry(data, apertures)
        self.phot_table = phot_table
        
        area = np.array([a.area for a in apertures])
        aperture_sum = np.array([self.phot_table[f'aperture_sum_{i}'][0] for i in range(len(bins)-1)])/area
        
        if self.obstype=='HST':
            aperture_sum_err = np.array([self.phot_table[f'aperture_sum_err_{i}'][0] for i in range(len(bins)-1)])/area
        else:
            aperture_sum_err = error / np.sqrt(area/self.NPixPerBeam)
        
        #aperture_sum_err = np.array([self.phot_table[f'aperture_sum_err_{i}'][0] for i in range(len(bins)-1)])/area
        
        if normalized:
            aperture_sum_err /= np.max(aperture_sum)
            aperture_sum /= np.max(aperture_sum)
            
        
        bc = 0.5*(blist[1:]+blist[:-1])*self.cell_kpc
        ac = 0.5*(alist[1:]+alist[:-1])*self.cell_kpc
        
        
        
        return bc, ac, aperture_sum, aperture_sum_err
        

        
    def AxisProfile(self, axis, FWHM=None, boxwidth=3, theta=35, theta_source=-12, normalized=True, average=False, N=20, kpc=True):
        print(f"Beam size: {self.header['BMAJ']*60*60:.2f}'' x {self.header['BMIN']*60*60:.2f}''")
        if self.plane=='image':
            x0,y0 = self.x0, self.y0
        elif self.plane=='source':
            x0,y0 = self.x0-1.8, self.y0-1.8
            theta = theta_source
            
        if FWHM == None:
            if self.plane=='source':
                FWHM = self.header['BMAJ']*60*60
                print(f"Using box height of {FWHM:.2f}''")
            else:
                FWHM = np.mean([self.header['BMAJ'],self.header['BMIN']])*60*60
                print(f"Using box height of {FWHM:.2f}''")
        
        theta *= np.pi/180 # convert theta to radians
                
        self.major_slope = np.cos(theta)/np.sin(theta)
        self.minor_slope = np.cos(theta-np.pi/2)/np.sin(theta-np.pi/2)
        
        boxheight = FWHM/self.cell # must be close to the beam size or 1/2 beam size! 
        boxwidth = boxwidth/self.cell

        if axis=='major':
            centers = np.array([(x0-np.sin(theta)*boxheight*i,y0+np.cos(theta)*boxheight*i) for i in np.arange(-N,N+1,1)])
            posangle = theta
        if axis=='minor':
            centers = np.array([(x0+np.cos(theta)*boxheight*i,y0+np.sin(theta)*boxheight*i) for i in np.arange(-N,N+1,1)])
            posangle = theta+np.pi/2
            
        from photutils.aperture import RectangularAperture, aperture_photometry
        self.data[np.isnan(self.data)] = 0 # set nan values to 0 (since astropy doesn't like them) but so they won't contribute to flux
        apertures = RectangularAperture(centers, h=boxheight, w=boxwidth, theta=posangle)
        self.aperture = apertures
        phot_table = aperture_photometry(self.data, apertures, error=self.error)
        
        delta_x = x0 - np.array(phot_table['xcenter'])
        delta_y = np.array(phot_table['ycenter']) - y0
        
        delta_RA = delta_x * self.cell
        delta_Dec = delta_y * self.cell
        d = np.sqrt(delta_RA**2 + delta_Dec**2) * np.sign(delta_x)
        sb = np.array(phot_table['aperture_sum'])
        sb_err = np.array(phot_table['aperture_sum_err'])

        if normalized:
            sb_err /= sb[np.argmin(np.abs(d))]
            sb /= sb[np.argmin(np.abs(d))]
            
        if average:
            sb = 0.5*(sb[d >= 0] + np.flip(sb[d <= 0]))
            sb_err = 0.5*np.sqrt(sb_err[d >= 0]**2 + np.flip(sb_err[d <= 0])**2)
            d = d[d >= 0]
            
        if self.plane=='source' and kpc:
            d = d/self.cell * self.cell_kpc
            
        return d, sb, sb_err
    
    def GrowthCurve(self, radii=np.arange(0.1, 6, 0.1)):
        from photutils.aperture import CircularAperture, aperture_photometry

        apertures = [CircularAperture([self.x0,self.y0], r/self.cell) for r in radii]

        error = self.total_error
        phot_table = aperture_photometry(self.data, apertures, error=error)
        self.phot_table = phot_table
        # error = np.sqrt(sum_{all pixels in aperture} sigma_{pixel}^2)
        # for HST data, this error calculation should be just fine, as we have pixel-level errors for the entire array
        # for ALMA data, we have a single sigma for the whole map but we need to scale by the beam size
        aperture_sum = np.array([self.phot_table[f'aperture_sum_{i}'][0] for i in range(len(radii))])
        aperture_sum_err = np.array([self.phot_table[f'aperture_sum_err_{i}'][0] for i in range(len(radii))])/np.sqrt(self.NPixPerBeam)
        
        aperture_sum_err /= aperture_sum[-1]
        aperture_sum /= aperture_sum[-1]
        
        return radii, aperture_sum, aperture_sum_err
    
    
    def Spectrum(self, aperture, restfreq=None, aperture_units='arcsec'):
        '''Returns the object's spectrum (frequency [GHz], flux [mJy]) in an specified aperture.
           Specify aperture as (x0,y0,R) where x0 = central right ascension in arcsec from center, 
           y0 = central declination in arcsec from center, R = radius in arcsec'''
        
        # if aperture is a tuple specifying (x0,y0,R)
        if type(aperture)==tuple:
            if len(aperture)==3: # circular aperture
                x0, y0, R = aperture
                x0 = -x0/self.cell + self.x0
                y0 = y0/self.cell + self.y0
                R = R/self.cell
                aperture = CircularAperture([x0,y0], R)
                self.aperture_patch = mpl.patches.Circle((-(x0-self.x0)*self.cell, (y0-self.y0)*self.cell),
                                                         radius=R*self.cell, fc='none', ec='w', lw=0.8, zorder=2000)
            elif len(aperture)==5: # elliptical aperture
                from photutils.aperture import EllipticalAperture
                x0, y0, a, b, theta = aperture
                x0 = -x0/self.cell + self.x0
                y0 = y0/self.cell + self.y0
                a = a/self.cell
                b = b/self.cell
                theta = theta*np.pi/180
                aperture = EllipticalAperture([x0,y0],a,b,theta)
                #self.aperture_patch = mpl.patches.Circle((x0, y0), radius=R, fc='none', ec='w', lw=0.8, zorder=2000)
                    
            self.aperture = aperture

            flux, flux_err = np.zeros(shape=self.nchans),np.zeros(shape=self.nchans)
            for v in range(self.nchans):
                im = self.data[v,:,:]
                phot_table = aperture_photometry(im, aperture, error=self.error)
                flux[v] = np.array(phot_table['aperture_sum'])[0]*1000
                flux_err[v] = np.array(phot_table['aperture_sum_err'])[0]*1000
         
        # if aperture is a boolean array (a mask)
        elif type(aperture)==np.ndarray:
            flux = np.zeros(shape=self.nchans)
            flux_err = np.zeros(shape=self.nchans)
            for v in range(self.nchans):
                im = self.data[v,:,:]
                flux[v] = np.sum(im[aperture])*1000
                flux_err[v] = np.sqrt(np.sum(np.power(self.error[aperture],2)))*1000
        
        # if aperture is a photutils.aperture object
        else: 
            self.aperture = aperture

            flux = np.zeros(shape=self.nchans)
            for v in range(self.nchans):
                im = self.data[v,:,:]
                phot_table = aperture_photometry(im, aperture)
                flux[v] = np.array(phot_table['aperture_sum'])[0]*1000

        if restfreq==None:
            x = self.freq
        else:
            self.restfreq = restfreq
            x = self.vel
        
        return x, flux, flux_err
    
    
    def Reconstruct(self, cell=0.035, beam=False, f=np.nanmean):
        '''Performs source plane reconstruction and replaces image with source-plane map.'''
        
        if not beam: 
            print('Performing source-plane reconstruction. Image-plane properties will be overwritten with source-plane properties.')
        dx, dy, mu = open_lens_model()
        x, y = np.arange(0, self.size, 1), np.arange(0, self.size, 1)
        x, y = np.meshgrid(x, y)
        x = (x - self.header['CRPIX1'])*self.header['CDELT1'] + self.header['CRVAL1']
        y = (y - self.header['CRPIX2'])*self.header['CDELT2'] + self.header['CRVAL2']
        
        x_source = x - dx
        y_source = y - dy
        
        crpix1 = self.size//2
        crpix2 = self.size//2
        crval1 = 197.87615818272903
        crval2 = -1.3371902006535314
        cdelt1 = cell/60/60
        cdelt2 = cell/60/60
        
        xbins = np.arange(0, self.size+1, 1)
        ybins = np.arange(0, self.size+1, 1)
        xbins = (xbins-crpix1)*cdelt1 + crval1
        ybins = (ybins-crpix2)*cdelt2 + crval2
        
        old_cell = self.cell
        
        self.data *= self.NPixPerBeam
        self.std *= self.NPixPerBeam
        self.error *= self.NPixPerBeam
        if self.obstype != 'HST':
            self.std_map *= self.NPixPerBeam
            self.residual *= self.NPixPerBeam
        
        from scipy.stats import binned_statistic_2d
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if self.obstype=='cube':
                x, y, z = x_source.flatten(), y_source.flatten(), self.mom0.flatten()
                values, xbins1, ybins1, binnumber = binned_statistic_2d(y, x, z, bins=(ybins, xbins), statistic=f)
                self.mom0 = values
            else:   
                x, y, z = x_source.flatten(), y_source.flatten(), copy(self.data).flatten()
                values, xbins1, ybins1, binnumber = binned_statistic_2d(y, x, z, bins=(ybins, xbins), statistic=f)
                self.data = values
                
                if not beam and self.obstype != 'HST':
                    z = copy(self.residual).flatten()
                    values, xbins1, ybins1, binnumber = binned_statistic_2d(y, x, z, bins=(ybins, xbins), statistic=f)
                    self.residual = values
                    
                    z = copy(self.std_map).flatten()
                    values, xbins1, ybins1, binnumber = binned_statistic_2d(y, x, z, bins=(ybins, xbins), statistic=f)
                    self.std_map = values

        
        w = WCS(naxis=2)
        w.wcs.crpix = [-crpix1,crpix2]
        w.wcs.cdelt = np.array([-cdelt1, cdelt2])
        w.wcs.crval = [crval1, crval2]
        w.wcs.ctype = ["RA---SIN", "DEC--SIN"]
        self.header.update(w.to_header())
        
        self.x0, self.y0 = self.size//2, self.size//2

        self.plane = 'source'
        
        cell_factor = (old_cell**2)/(cell**2)
        
        if not beam:
            self.update_beam_size()
            
            self.data /= self.NPixPerBeam
            self.data *= cell_factor
            self.std /= self.NPixPerBeam
            self.std *= cell_factor
            
            if self.obstype != 'HST':
                self.std_map *= cell_factor
                self.error = self.std_map * np.sqrt(self.NPixPerBeam) 
                self.std_map /= self.NPixPerBeam
                self.residual *= cell_factor
                self.residual /= self.NPixPerBeam
            self.error /= self.NPixPerBeam
            
            
        
        
    def update_beam_size(self):
        # compute new beam size for reconstructed map
        Bmaj = self.header['Bmaj']*60*60/self.cell # in pixels
        Bmin = self.header['Bmin']*60*60/self.cell # in pixels
        theta = self.header['BPA']*np.pi/180 # in radians
        #print(f"Beam Size: {Bmaj*im.cell:.2f}'' x {Bmin*im.cell:.2f}'', BPA = {theta/np.pi*180:.2f} degrees")

        sigma_maj = Bmaj / (2*np.sqrt(2*np.log(2)))
        sigma_min = Bmin / (2*np.sqrt(2*np.log(2)))

        # produce gaussian model for the beam
        from astropy.convolution import Gaussian2DKernel
        beam_model = Gaussian2DKernel(sigma_maj,sigma_min,theta*np.pi/180-np.pi/2, x_size=self.size, y_size=self.size)
        
        try:
            im = image(self.line, self.obstype, self.weighting)
        except AttributeError: # for HST images
            im = image(self.psfmatch[0], self.psfmatch[1], self.psfmatch[2])
            
        
        im.data = beam_model.array # set psf array to the gaussian model for the psf (simpler!)
        im.Reconstruct(beam=True) # use beam=True to tell Reconstruct to not run update_beam_size again
        im.data[np.isnan(im.data)] = 0 # remove nan values from resulting reconstructed image
        
        # fit a Gaussian to the source-plane beam
        from astropy.modeling import models, fitting

        p_init = models.Gaussian2D(amplitude=1, x_mean=im.x0, y_mean=im.y0, x_stddev=0.5/im.cell, y_stddev=0.1/im.cell, theta=90/180*np.pi)
        fit_p = fitting.LevMarLSQFitter()
        x, y = np.arange(0, im.size, 1), np.arange(0, im.size, 1)
        x, y = np.meshgrid(x, y)
        p = fit_p(p_init, x, y, im.data, maxiter=1000)
        
        self.header['Bmaj'] = p.y_stddev.value * im.cell * 2.355 / 60 /60
        self.header['Bmin'] = p.x_stddev.value * im.cell * 2.355 / 60 /60
        self.header['BPA'] = p.theta.value*180/np.pi - 180
        
        self.header['Bmaj_kpc'] = p.y_stddev.value * im.cell_kpc * 2.355 
        self.header['Bmin_kpc'] = p.x_stddev.value * im.cell_kpc * 2.355 
            