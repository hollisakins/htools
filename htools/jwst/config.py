

all_fields = ['cosmos-web', 'primer-cosmos', 'ceers', 'cosmos']

all_hst_bands = ['f435w', 'f606w', 'f814w', 'f105w', 'f125w', 'f140w', 'f160w']
all_jwst_bands = ['f070w', 'f090w', 'f115w', 'f150w', 'f200w', 'f277w', 'f356w', 'f444w',
                  'f140m', 'f162m', 'f182m', 'f210m', 'f250m', 'f300m', 'f335m', 'f360m', 'f410m', 'f430m', 'f460m', 'f480m', 
                  'f560w', 'f770w', 'f1000w', 'f1130w', 'f1280w', 'f1500w', 'f1800w', 'f2100w', 'f2550w']
misc_bands = ['u','g','r','i','z','y', 'NB118','Y','J','H','Ks', 'IRAC1','IRAC2','IRAC3','IRAC4']

bands = {
    'cosmos-web':['f814w','f115w','f150w','f277w','f444w','f770w'],
    'primer-cosmos':['f606w','f814w','f090w','f115w','f150w','f200w','f277w','f356w','f410m','f444w','f770w','f1800w'],
    'ceers':['f606w','f814w','f115w','f150w','f200w','f277w','f356w','f410m','f444w'],
    'cosmos': ['u','g','r','i','z','y', 'NB118','Y','J','H','Ks', 'IRAC1','IRAC2','IRAC3','IRAC4'],
    }

colors = {
    'cosmos-web':['darkmagenta','#0088e7', '#03a1a1', '#83b505','#ab0202','#e74001'],
    'primer-cosmos':['indigo', 'darkmagenta','#0c00e7', '#0088e7', '#03a1a1', '#009a49', '#83b505', '#c2a206','darkorange','#ab0202','#e74001','#630606']
    }

