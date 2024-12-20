
import numpy as np

class LephareResult:

    def __init__(self, ID, zspec, zphot, phot, pz, models):
        self.ID = ID
        self.zspec = zspec
        self.zphot = zphot
        self.phot = phot
        self.pz = pz
        self.models = models

    def read(filepath):
        ID, zspec, zphot = np.loadtxt(filepath, skiprows=1, max_rows=1, usecols=(0, 1, 2)).T
        ID = int(ID)
        if zspec < 0:
            zspec = None
        else:
            zspec = float(zspec)
        if zphot < 0: 
            zphot = None
        else:
            zphot = float(zphot)

        nfilters = int(np.loadtxt(filepath, skiprows=3, max_rows=1, usecols=1, dtype=int))
        npz = int(np.loadtxt(filepath, skiprows=5, max_rows=1, usecols=1, dtype=int))
        
        i = 13

        mag_obs, mag_obs_err, mag_model = np.loadtxt(filepath, skiprows=i, max_rows=nfilters, usecols=(0, 1, 8)).T
        fnu_obs, fnu_model = 3631e6*np.power(10., -0.4*mag_obs), 3631e6*np.power(10., -0.4*mag_model)
        fnu_obs_err = mag_obs_err * np.log(10)/2.5 * fnu_obs
        phot = LepharePhot(fnu_obs, fnu_obs_err, fnu_model)
        i += nfilters

        zgrid, Pz, PzB = np.loadtxt(filepath, skiprows=i, max_rows=npz, usecols=(0, 1, 2)).T
        pz = LepharePz(zgrid, Pz, PzB)
        i += npz

        # Type Nline Model Library Nband  Zphot Zinf Zsup Chi2  PDF  Extlaw EB-V Lir Age  Mass SFR SSFR
        mod_type = np.loadtxt(filepath, skiprows=7, max_rows=6, dtype=str, usecols=0)
        Nrows = np.loadtxt(filepath, skiprows=7, max_rows=6, dtype=int, usecols=1)
        zphots = np.loadtxt(filepath, skiprows=7, max_rows=6, dtype=float, usecols=5)
        chi2 = np.loadtxt(filepath, skiprows=7, max_rows=6, dtype=float, usecols=8)


        models = {}
        for t,n,c,z in zip(mod_type,Nrows,chi2,zphots):
            x, y = np.loadtxt(filepath, skiprows=i, max_rows=n, usecols=(0, 1)).T
            models[t] = {'wav_obs':x, 'fnu':3631e6*np.power(10., -0.4*y), 'chi2':c, 'zphot':z}
            i += n
        
        return LephareResult(ID, zspec, zphot, phot, pz, models)



class LepharePz:
    def __init__(self, zgrid, Pz, Pz_bayesian):
        self.zgrid = zgrid
        self.Pz = Pz
        self.Pz_bayesian = Pz_bayesian

        self.normalized = False

    def normalize(self):
        if self.normalized:
            return
        self.Pz /= np.trapz(self.Pz, x=self.zgrid)
        self.Pz_bayesian /= np.trapz(self.Pz_bayesian, x=self.zgrid)
        self.Pz /= np.max(self.Pz)
        self.Pz_bayesian /= np.max(self.Pz)
        self.normalized = True



class LepharePhot:
    def __init__(self, obs, obs_err, model, filter_list = None):
        self.obs = obs
        self.obs_err = obs_err
        self.model = model
        if filter_list is None:
            self.filter_list = ['cfht_u','hsc_g','hsc_r','hsc_i','hsc_z','hsc_y','uvista_Y','uvista_J','uvista_H','uvista_Ks','f814w','f115w','f150w','f277w','f444w','f770w','IB427','IB484','IB505','IB527','IB574','IB624','IB679','IB709','IB738','IB767','IB827','NB711','NB816','HSC-NB816','HSC-NB921','HSC-NB1010','NB118','FUV','NUV','U','V']
        else:
            self.filter_list = filter_list

    def get_filters(self, filters):
        mo, moe, mm = [], [], []
        for f in filters:
            i = np.where(np.array(self.filter_list) == f)[0][0]
            mo.append(self.obs[i])
            moe.append(self.obs_err[i])
            mm.append(self.model[i])
        return LepharePhot(mo, moe, mm, filter_list=filters)
        
    def __len__(self):
        return len(self.filter_list)
