### Far Infrared SED Fitting Procedures
from htools import *
import numpy as np
from scipy.stats import logistic
from astropy import units as u
from astropy.cosmology import Planck18 as cosmo
from numpy.random import Generator, SFC64
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde
from scipy import interpolate
import os
from dotmap import DotMap


os.environ["OMP_NUM_THREADS"] = "1"

T_CMB_z0 = 2.73

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################

from astropy.constants import c, h, k_B
hck = (h*c/k_B).to(u.K*u.micron).value
c = c.to(u.micron/u.s).value



def BB(lam, Nbb, T, beta): 
    '''Basic blackbody spectrum given by Planck's Law.'''
    return np.power(10, Nbb) * np.power(c/lam, 3) / (np.exp(hck/(lam*T))-1)

def MBB(lam, Nbb, T, beta, lam0):
    '''
    Modified blackbody spectrum.
        * lam0 : wavelength at which tau = 1 
          for optically thin case, set lam0 = None
          for general opacity model, start with lam0 = 200 * u.micron
    '''
    if lam0 == 'ot': # optically thin 
        return np.power(c/lam, beta) * BB(lam, Nbb, T, beta)
    else: # general opacity
        return (1-np.exp(-(lam0/lam)**beta)) * BB(lam, Nbb, T, beta)

def derivativeLogMBB(T, beta, lam0):
    '''Estimate the derivative of the modified blackbody'''
    lam_fine = np.logspace(0.1, 2.5, 1000)
    log_MBB = np.log10(MBB(lam_fine, 1, T, beta, lam0))
    delta_y = np.diff(log_MBB)
    delta_x = np.diff(np.log10(lam_fine))
    deriv = delta_y / delta_x
    return deriv

def lam_intersect(alpha, T, beta, lam0):
    '''Compute the wavelength where the derivative of the log of MBB equals the slope of the power law'''
    MBB_deriv = np.flip(derivativeLogMBB(T, beta, lam0))
    lam_fine = np.logspace(0.1, 2.5, 1000)
    return lam_fine[np.searchsorted(MBB_deriv, alpha)]

def powerLaw(lam, Npl, alpha):
    """Equation of the power law portion of SED"""
    return Npl * lam**alpha

def S(lam, Nbb, T, beta, alpha, lam0, z):
    '''Combined modified blackbody and power law.'''
    lam_int = lam_intersect(alpha, T, beta, lam0)
    mbb = MBB(lam, Nbb, T, beta, lam0)
    Npl = MBB(lam_int, Nbb, T, beta, lam0) * lam_int**(-alpha)
    pl = powerLaw(lam, Npl, alpha)
    return np.where(lam < lam_int, pl, mbb)
    
def CMB(lam,T,beta, z):
    '''Correction for CMB. Ratio of f_observed / f_intrinsic from da Cunha et al. 2013'''
    T_CMB = T_CMB_z0*(1+z)
    return 1 - (BB(lam, 1, T_CMB, beta)/BB(lam, 1, Tdust_z(T, beta, z), beta))

def Tdust_z(T, beta, z):
    '''Correction for CMB. From da Cunha et al. 2013 but adapted from Novak et al. 2019'''
    return (T**(beta+4) + T_CMB_z0**(beta+4) * ((1+z)**(beta+4)-1))**(1/(beta+4))

def S_CMB(lam, Nbb, T, beta, alpha, lam0, z):
    '''Combined modified blackbody and power law, corrected for the CMB heating following da Cunha et al. (2013).'''
    return S(lam, Nbb, Tdust_z(T,beta,z), beta, alpha, lam0, z) * CMB(lam, Tdust_z(T,beta,z),beta, z)


def MBB_CMB(lam, Nbb, T, beta, lam0, z):
    '''Combined modified blackbody and power law, corrected for the CMB heating following da Cunha et al. (2013).'''
    return MBB(lam, Nbb, Tdust_z(T,beta,z), beta, lam0) * CMB(lam, Tdust_z(T,beta,z),beta, z)
        

def fourPiLumDistSquared(z):
    '''4 pi luminosity distance squared'''
    return 4*np.pi*(cosmo.luminosity_distance(z).to(u.Mpc).value)**2 

def logLIR(lam, Nbb, T, beta, alpha, lam0, z):
    """Calculate LIR"""
    x_lam = np.linspace(8, 1000, 10000)
    x_nu = c/x_lam
    return np.log10(-np.trapz(S(x_lam, Nbb, T, beta, alpha, lam0, z), x=x_nu) / (1+z) * fourPiLumDistSquared(z) * 2.4873056783618645e-11)  # mJy Hz Mpc^2 to Lsun

def lambdaPeak(Nbb, T, beta, alpha, lam0, z):
    """Calculate Peak Wavelength"""
    x_lam = np.logspace(np.log10(8), 3, 10000)
    return x_lam[np.argmax(S(x_lam, Nbb, T, beta, alpha, lam0, z))]


# def L_IR(Nbbrand, Trand, betarand, alpha=2, lam0=200, fullout=False):
#     N1, N2, N3 = len(Nbbrand), len(Trand), len(betarand)
#     assert N1==N2
#     assert N2==N3
#     N = N1

#     FIR_rand = np.zeros(N)
#     x_lam = np.linspace(8, 1000, 1000) * u.micron
#     x_nu = c/x_lam
#     for i in tqdm.tqdm(range(N)):
#         Nbb = Nbbrand[i]
#         T = Tdust_rand[i]
#         beta = beta_rand[i]
#         Snu = S(x_lam, Nbb, T, beta, alpha, lam0)
#         FIR_rand[i] = np.trapz(Snu, x=x_nu)

#     L_IR_rand = (FIR_rand / (1+z) * fourPiLumDistSquared(z)).to(u.Lsun)

#     if fullout:
#         return L_IR_rand
#     else:
#         print('returning L_IR, -L_IR_err, +L_IR_err')
#         return np.median(L_IR_rand), np.median(L_IR_rand)-np.percentile(L_IR_rand, 16), np.percentile(L_IR_rand, 84)-np.median(L_IR_rand)


def SFR_IR(L_IR):
    '''Relation from Murphy et al. 2011'''
    return (L_IR.to(u.erg/u.s)).value * 3.88e-44 * u.Msun/u.yr
    

##########################################################################################################################
################################################### FITTING PROCEDURES ###################################################
##########################################################################################################################


from scipy.optimize import fmin_slsqp

def model(theta, func, lam, z, fix_alpha, fix_lam0):
    if fix_alpha is not None and fix_lam0 is not None: # both parameters fixed
        Nbb, T, beta = theta
        alpha, lam0 = fix_alpha, fix_lam0
    elif fix_alpha is not None and fix_lam0 is None:
        Nbb, T, beta, lam0 = theta
        alpha = fix_alpha
    elif fix_alpha is None and fix_lam0 is not None:
        Nbb, T, beta, alpha = theta
        lam0 = fix_lam0
    else:
        Nbb, T, beta, alpha, lam0 = theta
    return func(lam, Nbb, T, beta, alpha, lam0, z)

def lnlike(theta, func, lam, f, f_err, z, fix_alpha, fix_lam0):
    f_model = model(theta, func, lam, z, fix_alpha, fix_lam0)
    return -0.5 * np.sum(np.power((f-f_model)/f_err,2))

def lnprior(theta, bounds):
    prior_cond = [(t>=b[0]) and (t<=b[1]) for t, b in zip(theta, bounds)]
    if all(prior_cond):
        return 0.
    else:
        return -np.inf

def lnprob(theta, func, lam, f, f_err, bounds, z, fix_alpha, fix_lam0):
    lp = lnprior(theta, bounds)
    if np.isinf(lp):
        return -np.inf
    return lp + lnlike(theta, func, lam, f, f_err, z, fix_alpha, fix_lam0)


def mcfirsedfit(lam, f, f_err, func=S, p0=None, bounds=(), fix_alpha=None, fix_lam0=None, z=0, nwalkers=500, niter=1000, parallel=False, writedir=None, overwrite=None):#, uplim_thresh=3):
    '''Perform MCMC fitting to the FIR SED. 
       Options include: 
        * func: function to fit to, should be one of S, S_ot, or S_CMB above
            in order to fix certain parameters (e.g. alpha or lam0), input e.g.:
            `func = lambda lam, Nbb, T, beta : S(lam, Nbb, T, beta 4., 200.)`
        * bounds : bounds for flat priors in the MCMC fitting
        * p0 : initial guesses for fitting parameters
        * nwalkers : number of MCMC random walkers (default 500)
        * niter : number of MCMC random walkers (default 1000)
        * parallel : whether to employ pool.map multiprocessing to speed up runtime
        * uplim_thresh: sigma threshold below which to consider a 3-sigma upper limit rather than an actual measurement 
    '''
    import emcee
    import os
    if writedir is not None:
        writepath = writedir + 'mcfirsed-latest.pickle'
        if os.path.exists(writepath):
            print(f'mcfirsed output already saved at {writepath}')
            if overwrite is None:
                yn = input('would you like to run and overwrite it? [y/n] ')
                if yn!='y':
                    print('returning saved output...')
                    import pickle
                    with open(writepath, 'rb') as f:
                        out = pickle.load(f)
                        return out
            if overwrite is False:
                print('returning saved output...')
                import pickle
                with open(writepath, 'rb') as f:
                    out = pickle.load(f)
                    return out




    if hasattr(lam, 'unit'):
        lam = np.array(lam.to(u.micron).value)
    if hasattr(f,'unit'):
        f = np.array(f.to(u.mJy).value)
    if hasattr(f_err,'unit'):
        f_err = np.array(f_err.to(u.mJy).value)
    
    

    if p0 is None:
        from scipy.optimize import curve_fit
        popt, pcov = curve_fit(func, lam, f, sigma=f_err, bounds=np.array(bounds).T)
        p0 = popt
    else:
        p0 = np.array(p0)

    print('Initial guesses:')
    print(p0)
    ndim = len(p0)
    
    steps = np.random.normal(loc=p0, scale=np.abs(0.001*p0), size=(nwalkers, ndim))

    import htools.istarmap
    from multiprocessing import Pool, cpu_count
    if parallel:
        ncores = cpu_count()
        print(f'Entering pool.istarmap multiprocessing with {ncores} cores')
    
    from contextlib import nullcontext
    with Pool() if parallel else nullcontext() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(func, lam, f, f_err, bounds, z, fix_alpha, fix_lam0), pool=pool)

        print('Running burn-in...')
        steps, _, _ = sampler.run_mcmc(steps, niter//5, progress=True)
        sampler.reset()

        print('Running production...')
        pos, prob, state = sampler.run_mcmc(steps, niter, progress=True)
    
    out = DotMap()
    out.samples = sampler.flatchain
    out.pmax = out.samples[np.argmax(sampler.flatlnprobability)]
    out.p16 = np.percentile(out.samples, 16, axis=0)
    out.p50 = np.percentile(out.samples, 50, axis=0)
    out.p84 = np.percentile(out.samples, 84, axis=0)

    out.x_lam = np.logspace(0.9, 3, 100000)
    n = len(out.samples)//100
    print(f'Computing median and 16th-84th percentile SEDs using {n} samples')
    f_rand = np.zeros((len(out.samples),len(out.x_lam)))
    # for i in tqdm.tqdm(range(len(out.samples))):
    #     f_rand[i,:] = model(out.samples[i], func, out.x_lam, z, fix_alpha, fix_lam0)

    iterable = [[out.samples[i], func, out.x_lam, z, fix_alpha, fix_lam0] for i in np.random.choice(range(len(out.samples)), size=n, replace=False)]
    if parallel:
        with Pool() as pool:
            f_rand = list(tqdm.tqdm(pool.istarmap(model, iterable), total=len(iterable)))
    else:
        f_rand = np.zeros((len(iterable),len(out.x_lam)))
        for i in tqdm.tqdm(range(len(iterable))):
            f_rand[i,:] = model(*tuple(iterable[i]))
        # f_rand = list(tqdm.tqdm(map(model, iterable[0], iterable[1], iterable[2], iterable[3], iterable[4], iterable[5]), total=len(iterable)))

    cs = np.percentile(f_rand, [16, 50, 84], axis=0)
    out.f16 = cs[0,:]
    out.f50 = cs[1,:]
    out.f84 = cs[2,:]

    print('Converting Nbb to LIR for plotting purposes')
    samples_LIR = copy(out.samples)
    iterable = [[out.samples[i], logLIR, lam, z, fix_alpha, fix_lam0] for i in range(len(out.samples))]
    
    if parallel:
        with Pool() as pool:
            logLIRs = list(tqdm.tqdm(pool.istarmap(model, iterable), total=len(iterable)))
    else:
        logLIRs = np.zeros((len(iterable)))
        for i in tqdm.tqdm(range(len(iterable))):
            logLIRs[i] = model(*tuple(iterable[i]))
    samples_LIR[:,0] = np.power(10, logLIRs)

    out.samples_LIR = samples_LIR

    # for i, var in enumerate(func.__code__.co_varnames[1:]):
        # out[var] = out.p50[i]

    out.LIR50 = np.percentile(np.power(10, logLIRs), 50)
    out.LIR16 = np.percentile(np.power(10, logLIRs), 16)
    out.LIR84 = np.percentile(np.power(10, logLIRs), 84)

    if writedir is not None:
        import pickle
        print(f'writing mcfirsed output to {writepath}')
        with open(writepath, 'wb') as f:
            pickle.dump(out, f)


    return out





        

    # '''Perform MCMC fitting to the FIR SED. 
    #    Options include: 
    #     * uplim_thresh: sigma threshold below which to consider a 3-sigma upper limit rather than an actual measurement 
    #     * func: function to fit to, should be one of S, S_ot above
    #     * coeff0 : initial guesses for fitting parameters
    #     * bounds : bounds for the fitting parameters, optional 
    #     * N : number of MCMC iterations (default 10,000)
    # '''
    # def chisq(coeff):
    #     '''chi-squared ptimization function. equals 1 when coeff is optimized'''
    #     r = newf - func(lam, *coeff)
    #     return np.sum((r/std)**2)
    
    # # set up constraints
    # # first determine which values are measurements vs. upper limits
    # uplims = np.where(f/f_err < 3)[0]
    # ieqcons = []
    # for l, std in zip(lam[uplims], f_err[uplims]):
    #     # constrains fitting such that all the values of the returned array will be >= 0 at the end
    #     # in this case, this requires that func(lam, *coeff) >= 3*std in the optimized solution 
    #     ieqcons.append(lambda coeff : 3*std - func(l, *coeff)) 
    #     # additionally we require that the flux is >=0 at the uplim position
    #     ieqcons.append(lambda coeff : func(l, *coeff)) 

    # if coeff0 == ():
    #     raise Exception('You need to specify initital guesses for the fitting parameters!')

    # coeff, fx, its, imode, smode = fmin_slsqp(chisq, coeff0, iprint=0, full_output=True, ieqcons=ieqcons)#, bounds=[(-100,-20),(10,80),(0.5,3)])

    # newcoeffs = np.zeros(shape=(N,len(coeff)))
    # np.random.seed(123)
    # for i in tqdm.tqdm(range(N)):
    #     with warnings.catch_warnings():
    #         warnings.simplefilter('ignore')
    #         newf = np.random.normal(loc=f, scale=std)
    #         newcoeffs[i,:] = fmin_slsqp(chisq, coeff, 
    #                                     ieqcons=ieqcons,
    #                                     #bounds=[(-100,-20),(10,80),(0.5,3)],
    #                                     iprint = 0)

    # yvalues = [func(x_lam, newcoeffs[i,0], newcoeffs[i,1],newcoeffs[i,2]) for i in range(N)]
    # yvalues = np.array(yvalues)
    # ymeans = np.mean(yvalues,axis=0)
    # ystds = np.std(yvalues,axis=0)

    # coeff = np.zeros(shape=np.shape(coeff0))
    # coeff_err_upper = np.zeros(shape=np.shape(coeff0))
    # coeff_err_lower = np.zeros(shape=np.shape(coeff0))
    # for i in range(len(coeff0)):
    #     cs = np.percentile(newcoeffs[:,i], [16, 50, 84])
    #     q = np.diff(cs)
    #     coeff_err_upper[i] = q[1]
    #     coeff_err_lower[i] = q[0]
    #     coeff[i] = cs[1]

    # class mcfirsedfit_output(self):
    #     def __init__(coeff, coeff_err_lower, coeff_err_upper, ymeans, ystds, x_lam):
    #         self.coeff = coeff
    #         self.coeff_err_lower = coeff_err_lower
    #         self.coeff_err_upper = coeff_err_upper
    #         self.ymeans = ymeans
    #         self.ystds = ystds
    #         self.x_lam = x_lam

    # return mcfitsedfit_output(coeff, coeff_err_lower, coeff_err_upper, ymeans, ystds)


##########################################################################################################################
################################################### PLOTTING FUNCTIONS ###################################################
##########################################################################################################################

def confRegion(ax,x,y,levels=[0.68,0.95],fill=False,smooth=None, bins=20, range=None, **kwargs):
    if range is None:
        range = [[x.min(), x.max()], [y.min(), y.max()]]
    H, X, Y = np.histogram2d(x.flatten(),y.flatten(),bins=bins,range=list(map(np.sort, range)))

    if smooth is not None:
        from scipy.ndimage import gaussian_filter
        H = gaussian_filter(H, smooth)

    Hflat = H.flatten()
    inds = np.argsort(Hflat)[::-1]
    Hflat = Hflat[inds]
    sm = np.cumsum(Hflat)
    sm /= sm[-1]
    V = np.empty(len(levels))
    for i, v0 in enumerate(levels):
        try:
            V[i] = Hflat[sm <= v0][-1]
        except IndexError:
            V[i] = Hflat[0]
    V.sort()
    m = np.diff(V) == 0
    while np.any(m):
        V[np.where(m)[0][0]] *= 1.0 - 1e-4
        m = np.diff(V) == 0
    V.sort()

    # Compute the bin centers.
    X1, Y1 = 0.5 * (X[1:] + X[:-1]), 0.5 * (Y[1:] + Y[:-1])

    # Extend the array for the sake of the contours at the plot edges.
    H2 = H.min() + np.zeros((H.shape[0] + 4, H.shape[1] + 4))
    H2[2:-2, 2:-2] = H
    H2[2:-2, 1] = H[:, 0]
    H2[2:-2, -2] = H[:, -1]
    H2[1, 2:-2] = H[0]
    H2[-2, 2:-2] = H[-1]
    H2[1, 1] = H[0, 0]
    H2[1, -2] = H[0, -1]
    H2[-2, 1] = H[-1, 0]
    H2[-2, -2] = H[-1, -1]
    X2 = np.concatenate(
        [
            X1[0] + np.array([-2, -1]) * np.diff(X1[:2]),
            X1,
            X1[-1] + np.array([1, 2]) * np.diff(X1[-2:]),
        ]
    )
    Y2 = np.concatenate(
        [
            Y1[0] + np.array([-2, -1]) * np.diff(Y1[:2]),
            Y1,
            Y1[-1] + np.array([1, 2]) * np.diff(Y1[-2:]),
        ]
    )


    if fill:
        cs = ax.contourf(X2, Y2, H2.T, np.concatenate([V, [H.max() * (1 + 1e-4)]]), **kwargs)
    else:
        cs = ax.contour(X2, Y2, H2.T, V, **kwargs)
        
    return cs

setattr(mpl.axes.Axes, "confRegion", confRegion)


def firsedplot(lam, f, f_err, func, coeff_opt, yunit=u.uJy, figsize=(4,3.5)):
    lam = lam.to(u.micron).value
    f = f.to(yunit).value
    f_err = f.to(yunit).value

    fig, ax = plt.subplots(1,1,figsize=figsize, constrained_layout=False)

    ax.errorbar(lam,f, yerr=f_err, mec='r', mfc='none', mew=0.75, ecolor='r', elinewidth=0.75, capsize=2, capthick=0.75, 
                marker='s', linestyle='none', label='Data')


    ax.plot(x_lam, func(x_lam, *coeff_opt).to(yunit).value, 
            color='k', linewidth=1)

    ax.fill_between(x_lam, 1e6*(ymeans-ystds), 1e6*(ymeans+ystds), ec='none', fc='0.7', alpha=0.5, zorder=-1000)

    ax.errorbar([CO_lam0], 3*CO_std*1e6, yerr=CO_std*1e6, uplims=True, 
                mec='r', mfc='none', mew=0.75, ecolor='r', elinewidth=0.75, capsize=2, capthick=0.75, 
                marker='s', linestyle='none')

    ax.loglog()
    ax.set_ylabel(r'Flux~$S\,\mu^{-1}$ [$\mu$Jy]')
    ax.set_xlabel(r'Rest Wavelength [$\mu$m]')
    ax.set_xlim(20, 1400)
    ax.set_ylim(1.1 ,2500)
    ax.legend(loc='upper left')
    ax.tick_params(direction='in', which='both')

    iax = ax.inset_axes([0.6,0.6,0.35,0.35])

    from scipy.stats import gaussian_kde
    xmin, xmax, ymin, ymax = 10, 100, 0.5, 3.5
    bins = (np.linspace(xmin,xmax,40),np.linspace(ymin,ymax,40))
    extent = [xmin,xmax,ymin,ymax]
    sigmas = np.array([1])
    levels = [1 - np.exp(-(x)**2/2) for x in sigmas]
    x = np.array(newcoeffs[:,1])
    y = np.array(newcoeffs[:,2])
    iax.confRegion(x, y, smooth=1.2, fill=True, levels=levels, bins=bins, colors='0.85', zorder=2, extent=extent, origin='lower')
    cs = iax.confRegion(x, y, smooth=1.2, fill=False, linewidths=0.5, levels=levels, bins=bins, colors='k', zorder=2, extent=extent, origin='lower')

    # iax.contourf(Z, levels=[0.68,0.95,1], colors=['0.8','0.6'], zorder=1, extent=[xmin,xmax,ymin,ymax], origin='lower')
    # cs = iax.contour(Z, levels=[0.68,0.95,1], linewidths=0.5, colors='k', zorder=2, extent=[xmin,xmax,ymin,ymax], origin='lower')
    key = {cs.levels[i]:np.flip(sigmas)[i] for i in range(len(levels))}
    def fmt(x):
        global key
        x = key[x]
        #s = f"{int(x*100)}"
        return fr"\bf {int(x)}$\sigma$"

    iax.clabel(cs, cs.levels, inline=True, fmt=fmt, fontsize=6.5, inline_spacing=0.5)
    iax.scatter([coeff[1]],[coeff[2]],marker='x',c='r',s=10, zorder=3)
    iax.set_xlabel(r'$T_{\rm dust}$ [K]')
    iax.set_ylabel(r'$\beta$')
    iax.set_xlim(25,60)
    iax.set_ylim(1,2.75)
    iax.tick_params(direction='in',which='both')
    iax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))
    iax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.25))

    ## change at will :)
    R = 2.8

    ax = plt.subplot(gs[0,2])
    im = image('OIII','continuum',weighting)
    i = ax.imshow(im.data, cmap='Greys_r', origin='lower', extent=im.extent, vmin=-3*im.std, vmax=5*im.std)
    ax.contour(im.data, extent=im.extent, colors='r', levels=np.array([2,3,5,8,12,16,20,25])*im.std, linewidths=0.5, zorder=999)
    ax.contour(im.data, extent=im.extent, colors='w', levels=np.array([-3,-2])*im.std, linewidths=0.5, zorder=999)
    ax.plotBeam(im, (0.75*R,-0.75*R))
    ax.annotate(r"90~$\mu$m", (0.05,0.95), xycoords='axes fraction', va='top', color='w', zorder=1002)
    ax.tick_params(labelbottom=False, labelleft=False, left=False, right=False, top=False, bottom=False)
    ax.set_xlim(R,-R)
    ax.set_ylim(-R,R)

    center = ((im.x0-x0)*im.cell, (y0-im.y0)*im.cell)

    el = mpl.patches.Ellipse(center, width=width, height=height, angle=angle, fc='none', ec='gold', linestyle='--', zorder=1003)
    ax.add_patch(el)

    ax.fill_between([0.95*R,0.35*R],[0.7*R,0.7*R],[0.95*R,0.95*R],color='k', zorder=1001)


    ax = plt.subplot(gs[0,1])
    im = image('OIII52','continuum',weighting)
    i = ax.imshow(im.data, cmap='Greys_r', origin='lower', extent=im.extent, vmin=-3*im.std, vmax=5*im.std)
    ax.contour(im.data, extent=im.extent, colors='r', levels=np.array([2,3,5,8,10])*im.std, linewidths=0.5, zorder=999)
    ax.contour(im.data, extent=im.extent, colors='w', levels=np.array([-3,-2])*im.std, linewidths=0.5, zorder=999)
    ax.plotBeam(im, (0.75*R,-0.75*R))
    ax.annotate(r"53~$\mu$m", (0.05,0.95), xycoords='axes fraction', va='top', color='w', zorder=1002)
    ax.tick_params(labelbottom=False, labelleft=False, left=False, right=False, top=False, bottom=False)
    ax.set_xlim(R,-R)
    ax.set_ylim(-R,R)

    el = mpl.patches.Ellipse(center, width=width, height=height, angle=angle, fc='none', ec='gold', linestyle='--', zorder=1003)
    ax.add_patch(el)

    ax.fill_between([0.95*R,0.35*R],[0.7*R,0.7*R],[0.95*R,0.95*R],color='k', zorder=1001)


    ax = plt.subplot(gs[1,1])
    im = image('NII','continuum',weighting)
    i = ax.imshow(im.data, cmap='Greys_r', origin='lower', extent=im.extent, vmin=-3*im.std, vmax=5*im.std)
    ax.contour(im.data, extent=im.extent, colors='r', levels=np.array([2,3,5,8,12,16,20,25])*im.std, linewidths=0.5, zorder=999)
    ax.contour(im.data, extent=im.extent, colors='w', levels=np.array([-3,-2])*im.std, linewidths=0.5, zorder=999)
    ax.plotBeam(im, (0.75*R,-0.75*R))
    ax.annotate(r"107~$\mu$m", (0.05,0.95), xycoords='axes fraction', va='top', color='w', zorder=1002)
    ax.tick_params(labelbottom=False, labelleft=False, left=False, right=False, top=False, bottom=False)
    ax.set_xlim(R,-R)
    ax.set_ylim(-R,R)

    el = mpl.patches.Ellipse(center, width=width, height=height, angle=angle, fc='none', ec='gold', linestyle='--', zorder=1003)
    ax.add_patch(el)

    ax.fill_between([0.95*R,0.25*R],[0.7*R,0.7*R],[0.95*R,0.95*R],color='k', zorder=1001)

    ax.plot([-0.85*R+2,-0.85*R],[-0.85*R,-0.85*R],color='w', linewidth=1.5, zorder=1000)
    ax.annotate("2.0''", (np.mean([-0.85*R+2,-0.85*R]),-0.75*R), color='w', ha='center', zorder=2000)


    ax = plt.subplot(gs[1,2])
    im = image('CII','continuum','natural_uv0.7')
    i = ax.imshow(im.data, cmap='Greys_r', origin='lower', extent=im.extent, vmin=-3*im.std, vmax=5*im.std)
    ax.contour(im.data, extent=im.extent, colors='r', levels=np.array([2,3,5,8,12,16,20,25])*im.std, linewidths=0.5, zorder=999)
    ax.contour(im.data, extent=im.extent, colors='w', levels=np.array([-3,-2])*im.std, linewidths=0.5, zorder=999)
    ax.plotBeam(im, (0.75*R,-0.75*R))
    ax.annotate(r"163~$\mu$m", (0.05,0.95), xycoords='axes fraction', va='top', color='w', zorder=1002)
    ax.tick_params(labelbottom=False, labelleft=False, left=False, right=False, top=False, bottom=False)
    ax.set_xlim(R,-R)
    ax.set_ylim(-R,R)

    el = mpl.patches.Ellipse(center, width=width, height=height, angle=angle, fc='none', ec='gold', linestyle='--', zorder=1003)
    ax.add_patch(el)
    ax.fill_between([0.95*R,0.25*R],[0.7*R,0.7*R],[0.95*R,0.95*R],color='k', zorder=1001)


    # plt.savefig('../Plots/FIR_SED_4.pdf')
    plt.show()