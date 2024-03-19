from htools.imaging import gen_rgb_image
from htools.utils import gen_cutout
import numpy as np

f606w = gen_cutout('primer-cosmos', '', 'f606w', (150.090,2.25), width=60, suffix='psfMatched')
f814w = gen_cutout('primer-cosmos', '', 'f814w', (150.090,2.25), width=60, suffix='psfMatched')
f090w = gen_cutout('primer-cosmos', '', 'f090w', (150.090,2.25), width=60, suffix='psfMatched')
f115w = gen_cutout('primer-cosmos', '', 'f115w', (150.090,2.25), width=60, suffix='psfMatched')
f150w = gen_cutout('primer-cosmos', '', 'f150w', (150.090,2.25), width=60, suffix='psfMatched')
f200w = gen_cutout('primer-cosmos', '', 'f200w', (150.090,2.25), width=60, suffix='psfMatched')
f277w = gen_cutout('primer-cosmos', '', 'f277w', (150.090,2.25), width=60, suffix='psfMatched')
f356w = gen_cutout('primer-cosmos', '', 'f356w', (150.090,2.25), width=60, suffix='psfMatched')
f410m = gen_cutout('primer-cosmos', '', 'f410m', (150.090,2.25), width=60, suffix='psfMatched')
f444w = gen_cutout('primer-cosmos', '', 'f444w', (150.090,2.25), width=60)
f150w.sci[np.isnan(f150w.sci)] = 0
f277w.sci[np.isnan(f277w.sci)] = 0
f444w.sci[np.isnan(f444w.sci)] = 0

input_dict = {}
input_dict['f606w'] = {'colors':np.array([0.0, 0.0, 1.0]), 'data':f606w.sci}
input_dict['f814w'] = {'colors':np.array([0.0, 0.0, 1.0]), 'data':f814w.sci}
input_dict['f090w'] = {'colors':np.array([0.0, 0.0, 1.0]), 'data':f090w.sci}
input_dict['f115w'] = {'colors':np.array([0.0, 0.5, 0.5]), 'data':f115w.sci}
input_dict['f150w'] = {'colors':np.array([0.0, 1.0, 0.0]), 'data':f150w.sci}
input_dict['f200w'] = {'colors':np.array([0.0, 1.0, 0.0]), 'data':f200w.sci}
input_dict['f277w'] = {'colors':np.array([0.3, 0.7, 0.0]), 'data':f277w.sci}
input_dict['f356w'] = {'colors':np.array([0.7, 0.3, 0.0]), 'data':f356w.sci}
input_dict['f410m'] = {'colors':np.array([1.0, 0.0, 0.0]), 'data':f410m.sci}
input_dict['f444w'] = {'colors':np.array([1.0, 0.0, 0.0]), 'data':f444w.sci}

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(4,4))

im = gen_rgb_image(input_dict,
                   noisesig=1, # map 1sigma flux to...
                   noiselum=0.1, # 10% grey
                   satpercent=0.8, # top 0.8% of pixels are saturated
                   save=False)
ax.imshow(im)
ax.axis('off')
plt.show()


