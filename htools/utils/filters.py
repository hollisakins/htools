import numpy as np
from astropy.units import Unit
import toml, os
from ..config import filter_directory

if hasattr(np, 'trapezoid'):
    trapz = np.trapezoid
else:
    trapz = np.trapz

class Filters(object):
    """Class for loading and manipulating sets of filter curves. 

    Parameters
    ----------

    filter_names : list
        List of names of filters, as defined in filter_directory.toml
    """

    def __init__(self, names, verbose=False):
        self.wavelengths = None
        self.verbose = verbose
        self.names = names
        self._load_filter_curves()
        self._calculate_effective_wavelengths()

    def _load_filter_curves(self):
        """ Loads filter files for the specified filter_names and truncates
        any zeros from either of their edges. """
        
        all_nicknames = {}
        
        self.filt_db = toml.load(filter_directory)
        for key in self.filt_db:
            for n in self.filt_db[key]['nicknames']:
                all_nicknames[n] = key
        
        self.filt_dict = {}
        self.nicknames = []
        
        if self.verbose: print("Loading filters")
        if self.verbose: print("-" * 80)
        l = np.max([len(i) for i in self.names])+3
        if self.verbose: print(f"Nickname".rjust(l) + ' -> ' + "Filter ID")
        for filt in self.names:
            if filt in all_nicknames:
                self.nicknames.append(all_nicknames[filt])
                self.filt_dict[filt] = np.loadtxt(os.path.join(os.path.split(filter_directory)[0], f'{all_nicknames[filt]}'))
                if self.verbose: print(f"{filt}".rjust(l) + ' -> ' +  f"{all_nicknames[filt]}") #+ f"{self.filt_db[all_nicknames[filt]]['description']}".ljust(32))
            else:
                raise ValueError(f"""Failed to match {filt} to any filter curve or nickname in database. Make sure it is named properly, or add it to the filter database.""")
                

        if self.verbose: print("-" * 80)

    def _calculate_effective_wavelengths(self):
        """ Calculates effective wavelengths for each filter curve. """

        self.wav = np.zeros(len(self))
        self.wav_min = np.zeros(len(self))
        self.wav_max = np.zeros(len(self))

        for i in range(len(self)):
            filt = self.names[i]

            w, T = self.filt_dict[filt][:, 0], self.filt_dict[filt][:, 1]
            self.wav[i] = trapz(w*T, x=w)/trapz(T, x=w)
            self.wav_min[i] = np.min(w[T/np.max(T) > 0.01])
            self.wav_max[i] = np.max(w[T/np.max(T) > 0.01])

        self.wav *= Unit('angstrom')
        self.wav_min *= Unit('angstrom')
        self.wav_max *= Unit('angstrom')

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.filt_dict[key]
        elif isinstance(key, int):
            return self.filt_dict[self.names[key]]
    
    def get_photometry(self, wav_obs: np.ndarray, f_nu: np.ndarray):
        """ Calculates photometric fluxes. The filters are first re-
        sampled onto the same wavelength grid with transmission values
        blueshifted by (1+z). This is followed by an integration over
        the observed spectrum in the rest frame:

        flux = integrate[(f_lambda*lambda*T(lambda*(1+z))*dlambda)]
        norm = integrate[(lambda*T(lambda*(1+z))*dlambda))]
        photometry = flux/norm

        lambda:            rest-frame wavelength array
        f_lambda:          observed spectrum
        T(lambda*(1+z)):   transmission of blueshifted filters
        dlambda:           width of each wavelength bin

        The integrals over all filters are done in one array operation
        to improve the speed of the code.
        """

        # Array containing filter profiles on new wavelength sampling
        filt_array = np.zeros((len(wav_obs), len(self.names)))

        for i in range(len(self.names)):
            filt = self.names[i]
            filt_array[:, i] = np.interp(wav_obs,
                                         self.filt_dict[filt][:, 0],
                                         self.filt_dict[filt][:, 1],
                                         left=0, right=0)

        print(filt_array)
        return
        
        # Calculate numerator of expression
        flux = np.expand_dims(getattr(sed, which)*self.widths*self.wavelengths, axis=1)
        flux = np.sum(flux*filters_z, axis=0)

        # Calculate denominator of expression
        norm = filters_z*np.expand_dims(self.widths*self.wavelengths, axis=1)
        norm = np.sum(norm, axis=0)

        photometry = np.squeeze(flux/norm)

        # # This is a little dodgy as pointed out by Ivo, it should depend
        # # on the spectral shape however only currently used for UVJ mags
        # if unit_conv == "cgs_to_mujy":
        #     photometry /= (10**-29*2.9979*10**18/self.eff_wavs**2)

        return photometry


    def __len__(self):
        return len(self.names)

    def __repr__(self):
        return f"Filters({', '.join(self.nicknames)})"





# def normal_asymm(loc=0, scale_up=1, scale_lo=1, size=1):
#     scale = np.random.choice([scale_up, scale_lo], size=1)
#     norm1 = np.random.normal(loc=loc, scale=scale_lo, size=size)
#     norm2 = np.random.normal(loc=loc, scale=scale_up, size=size)


# def get_filter(filt, filterpath='/Users/hba423/Dropbox/research/COSMOS-Web/archive/DREaM/bagpipes/bagpipes_filters/'):
#     '''Returns filter transmission curve (w,t) = (wavelength in angstroms, transmission).'''
#     f = np.loadtxt(filterpath+filt+'.dat')
#     w = f[:,0]
#     t = f[:,1]
#     return w,t

# def get_fnu(x, y, filt):
#     '''Compute fnu in a particular filter by convolving the filter curve with an SED
#     Assumed the filter file has wavelengths in angstroms, and the input SED has wavelengths in microns.'''
#     s = interp1d(x, y, bounds_error=False, fill_value=0)
#     if type(filt)==str:
#         w, fi = get_filter(filt)
#         lam = w/1e4 # in micron
#         nu = 299800/lam # in GHz
#         fnu = np.trapz(s(lam)*fi/nu, x=nu)/np.trapz(fi/nu, x=nu)

#     else:
#         fnu = np.zeros(len(filt))
#         for i in range(len(filt)):
#             w, fi = get_filter(filt[i])
#             lam = w/1e4 # in micron
#             nu = 299800/lam # in GHz
#             fnu[i] = np.trapz(s(lam)*fi/nu, x=nu)/np.trapz(fi/nu, x=nu)

#     return fnu

# def get_flam(x, y, filt):
#     '''Compute flam in a particular filter by convolving the filter curve with an SED
#     Assumed the filter file has wavelengths in angstroms, and the input SED has wavelengths in microns.'''
#     s = interp1d(x, y, bounds_error=False, fill_value=0)
#     if type(filt)==str:
#         w, fi = get_filter(filt)
#         lam = w/1e4 # in micron
#         flam = np.trapz(s(lam)*fi/lam, x=lam)/np.trapz(fi/lam, x=lam)

#     else:
#         flam = np.zeros(len(filt))
#         for i in range(len(filt)):
#             w, fi = get_filter(filt[i])
#             lam = w/1e4 # in micron
#             flam[i] = np.trapz(s(lam)*fi/lam, x=lam)/np.trapz(fi/lam, x=lam)

#     return flam



# def get_lameff(filt):
#     if type(filt)==str:
#         w, fi = get_filter(filt)
#         lam = w/1e4 # in micron
#         lam_mean = np.trapz(lam*fi, x=lam)/np.trapz(fi, x=lam)
#     else:
#         lam_mean = np.zeros(len(filt))
#         for i in range(len(filt)):
#             w, fi = get_filter(filt[i])
#             lam = w/1e4 # in micron
#             lam_mean[i] = np.trapz(lam*fi, x=lam)/np.trapz(fi, x=lam)

#     return lam_mean


# def get_lamrange(filt):
#     if type(filt)==str:
#         wi, fi = get_filter(filt)
#         wi /= 1e4
#         fi /= np.max(fi)
#         w_min = np.min(wi[fi>0.5])
#         w_max = np.max(wi[fi>0.5])
#     else:
#         w_min, w_max = np.zeros(len(filt)),np.zeros(len(filt))
#         for i in range(len(filt)):
#             wi, fi = get_filter(filt[i])
#             wi /= 1e4
#             fi /= np.max(fi)
#             w_min[i] = np.min(wi[fi>0.5])
#             w_max[i] = np.max(wi[fi>0.5])

#     return w_min, w_max
    

