import os, sys
import numpy as np
from astropy.io import fits

def readEazyBinary(main='photz', outdir='./output'):
    """
    Read Eazy binary output files
    """
        
    root = os.path.join(outdir, main) 
    cache_file = root+'.tempfilt'
    
    if os.path.exists(cache_file) is False:
        print(('File, %s, not found.' %(cache_file)))
        return -1,-1,-1,-1
    
    f = open(cache_file,'rb')
    
    s = np.fromfile(file=f,dtype=np.int32, count=4)
    NFILT=s[0]
    NTEMP=s[1]
    NZ=s[2]
    NOBJ=s[3]
    tempfilt = np.fromfile(file=f,dtype=np.double,count=NFILT*NTEMP*NZ).reshape((NZ,NTEMP,NFILT)).transpose()
    lc = np.fromfile(file=f,dtype=np.double,count=NFILT)
    zgrid = np.fromfile(file=f,dtype=np.double,count=NZ)
    fnu = np.fromfile(file=f,dtype=np.double,count=NFILT*NOBJ).reshape((NOBJ,NFILT)).transpose()
    efnu = np.fromfile(file=f,dtype=np.double,count=NFILT*NOBJ).reshape((NOBJ,NFILT)).transpose()
    
    f.close()
    
    tempfilt  = {'NFILT':NFILT,'NTEMP':NTEMP,'NZ':NZ,'NOBJ':NOBJ,\
                 'tempfilt':tempfilt,'lc':lc,'zgrid':zgrid,'fnu':fnu,'efnu':efnu}
    
    ###### .coeff
    f = open(root+'.coeff','rb')
    
    s = np.fromfile(file=f,dtype=np.int32, count=4)
    NFILT=s[0]
    NTEMP=s[1]
    NZ=s[2]
    NOBJ=s[3]
    coeffs = np.fromfile(file=f,dtype=np.double,count=NTEMP*NOBJ).reshape((NOBJ,NTEMP)).transpose()
    izbest = np.fromfile(file=f,dtype=np.int32,count=NOBJ)
    tnorm = np.fromfile(file=f,dtype=np.double,count=NTEMP)
    
    f.close()
    
    coeffs = {'NFILT':NFILT,'NTEMP':NTEMP,'NZ':NZ,'NOBJ':NOBJ,\
              'coeffs':coeffs,'izbest':izbest,'tnorm':tnorm}
              
    ###### .temp_sed
    f = open(root+'.temp_sed','rb')
    s = np.fromfile(file=f,dtype=np.int32, count=3)
    NTEMP=s[0]
    NTEMPL=s[1]
    NZ=s[2]
    templam = np.fromfile(file=f,dtype=np.double,count=NTEMPL)
    temp_seds = np.fromfile(file=f,dtype=np.double,count=NTEMPL*NTEMP).reshape((NTEMP,NTEMPL)).transpose()
    da = np.fromfile(file=f,dtype=np.double,count=NZ)
    db = np.fromfile(file=f,dtype=np.double,count=NZ)
    
    f.close()
    
    temp_sed = {'NTEMP':NTEMP,'NTEMPL':NTEMPL,'NZ':NZ,\
              'templam':templam,'temp_seds':temp_seds,'da':da,'db':db}
                
    ###### .pz
    if os.path.exists(root+'.pz'):
        f = open(root+'.pz','rb')
        s = np.fromfile(file=f,dtype=np.int32,count=2)
        NZ=s[0]
        NOBJ=s[1]
        chi2fit = np.fromfile(file=f,dtype=np.double,count=NZ*NOBJ).reshape((NOBJ,NZ)).transpose()
        p = np.exp(-chi2fit)
        p /= np.sum(p)
        pz = {'NZ':NZ,'NOBJ':NOBJ, 'chi2fit':chi2fit,'zgrid':zgrid, 'pz':p }
        
        f.close()

    return tempfilt, coeffs, temp_sed, pz


def write_eazy_catalog(main='photz', outdir='./output'):
    print(f'Loading EAzY binaries for {main}')
    tempfilt, coeffs, temp_seds, pz = readEazyBinary(main, outdir)

    print(f'Performing IGM correction...')
    temp_sed = (np.dot(temp_seds['temp_seds'],coeffs['coeffs']).T*(temp_seds['templam']/5500.)**2).T
    lim1 = np.where(temp_seds['templam'] < 912)
    temp_sed[lim1,:] *= 0
    lim2 = np.where((temp_seds['templam'] >= 912) & (temp_seds['templam'] < 1026))
    db = 1.-temp_seds['db'][coeffs['izbest']]
    temp_sed[lim2,:] *= db
    lim3 = np.where((temp_seds['templam'] >= 1026) & (temp_seds['templam'] < 1216))
    da = 1.-temp_seds['da'][coeffs['izbest']]
    temp_sed[lim3,:] *= da
    
    # zbest = tempfilt['zgrid'][coeffs['izbest']]
    temp_sed = temp_sed.T
    temp_lam = temp_seds['templam']/1e4

    #temp_lam = temp_lam*(1+zbest[:,np.newaxis])

    print(f'Generating FITS PROPERTIES catalog...')
    ID, z_spec, z_a, z_m1, chi2, zl68, zu68, zl95, zu95 = np.loadtxt(f'{outdir}/{main}.zout', skiprows=2, usecols=(0,1,2,3,4,5,6,7,8)).T
    columns = []
    columns.append(fits.Column(name=f'ID', array=ID, format='K'))
    columns.append(fits.Column(name=f'z_spec', array=z_spec, format='E'))
    columns.append(fits.Column(name=f'z_a', array=z_a, format='E'))
    columns.append(fits.Column(name=f'chi2', array=chi2, format='E'))
    columns.append(fits.Column(name=f'z_m1', array=z_m1, format='E'))
    columns.append(fits.Column(name=f'zl68', array=zl68, format='E'))
    columns.append(fits.Column(name=f'zu68', array=zu68, format='E'))
    columns.append(fits.Column(name=f'zl95', array=zl95, format='E'))
    columns.append(fits.Column(name=f'zu95', array=zu95, format='E'))
    t1 = fits.BinTableHDU.from_columns(fits.ColDefs(columns), header=fits.Header({'EXTNAME':'PROPERTIES'}))

    print(f'Generating FITS TEMP_SED catalog...')
    ID = np.append(-1, ID)
    spec = np.concatenate((temp_lam[np.newaxis,:], temp_sed), axis=0)
    columns = []
    columns.append(fits.Column(name='ID', array=ID, format='K'))
    columns.append(fits.Column(name='spec', array=spec, format=f'{np.shape(spec)[1]}E', unit='uJy'))
    t2 = fits.BinTableHDU.from_columns(fits.ColDefs(columns), header=fits.Header({'EXTNAME':'TEMP_SED'}))

    print(f'Generating FITS Pz catalog...')
    z = pz['zgrid']
    pz = pz['pz'].T
    Pz = np.concatenate((z[np.newaxis,:], pz), axis=0)
    columns = []
    columns.append(fits.Column(name='ID', array=ID, format='K'))
    columns.append(fits.Column(name='Pz', array=Pz, format=f'{np.shape(Pz)[1]}E'))
    t3 = fits.BinTableHDU.from_columns(fits.ColDefs(columns), header=fits.Header({'EXTNAME':'PZ'}))

    print(f'Generating FITS OBS_SED catalog...')
    obs_sed = tempfilt['fnu'].T
    obs_err = tempfilt['efnu'].T
    obs_lam = tempfilt['lc']/1e4
    Nbands = len(obs_lam)

    phot = np.concatenate((obs_lam[np.newaxis,:], obs_sed), axis=0)
    phot_err = np.concatenate((obs_lam[np.newaxis,:], obs_err), axis=0)

    columns = []
    columns.append(fits.Column(name='ID', array=ID, format='K'))
    columns.append(fits.Column(name='phot', array=phot, format=f'{Nbands}D', unit='uJy'))
    columns.append(fits.Column(name='phot_err', array=phot_err, format=f'{Nbands}D', unit='uJy'))
    t4 = fits.BinTableHDU.from_columns(fits.ColDefs(columns), header=fits.Header({'EXTNAME':'OBS_SED'}))

    print(f'Writing to {outdir}/{main}.fits')
    t0 = fits.PrimaryHDU(header=fits.Header({'EXTEND':'T'}))
    t = fits.HDUList([t0,t1,t2,t3,t4])
    t.writeto(f'{outdir}/{main}.fits')

def load_eazy_catalog(main='photz', outdir='./output', ext=0):
    if not os.path.exists(f'{outdir}/{main}.fits'):
        print(f"Looks like you haven't yet generated a FITS summary catalog ({main}.fits)")
        decision = input(f"Would you like to now? [y/n] ")
        if decision=='y':
            write_eazy_catalog(main, outdir)
            return fits.getdata(f'{outdir}/{main}.fits', ext=ext)
        else:
            print('Exiting...')
            sys.exit()
    else:
        return fits.getdata(f'{outdir}/{main}.fits', ext=ext)


def get_obs_sed(main='photz', outdir='./output', which='all'):
    cat = load_eazy_catalog(main, outdir, ext=4)
    wav = cat['phot'][0,:]
    phot = cat['phot'][1:,:]
    phot_err = cat['phot_err'][1:,:]
    ids = cat['ID'][1:]

    if which=='all':
        return wav, phot, phot_err
    else:
        i = np.argwhere(ids==which)[0][0]
        return wav, phot[i,:], phot_err[i,:]


def get_tem_sed(main='photz', outdir='./output', which='all', obs_frame=True):
    cat = load_eazy_catalog(main, outdir, ext=2)
    wav = cat['spec'][0,:]
    fnu = cat['spec'][1:,:]
    ids = cat['ID'][1:]
    z = load_eazy_catalog(main, outdir, ext=1)['z_a']

    if which=='all':
        if obs_frame:
            return wav*(1+z[:,np.newaxis]), fnu, z
        else:
            return wav, fnu, z
    else:
        i = np.argwhere(ids==which)[0][0]
        if obs_frame:
            return wav*(1+z[i]), fnu[i,:], z[i]
        else:
            return wav, fnu[i,:], z[i]

def get_pz(main='photz', outdir='./output', which='all'):
    cat = load_eazy_catalog(main, outdir, ext=3)
    z = cat['Pz'][0,:]
    Pz = cat['Pz'][1:,:]
    ids = cat['ID'][1:]
    if which=='all':
        return z, Pz
    else:
        i = np.argwhere(ids==which)[0][0]
        return z, Pz[i,:]


# def get_chi2(main='photz', outdir='./output', which='all'):
#     _, _, _, pz = readEazyBinary(main, outdir)
#     z = pz['zgrid']
#     chi2 = pz['chi2fit'].T
#     if which=='all':
#         return z, chi2
#     else:
#         f = np.loadtxt(outdir+'/'+main+'.zout')
#         idx = f[:,0]
#         cond = np.where(idx==which)[0]
#         return z, chi2[cond][0]


