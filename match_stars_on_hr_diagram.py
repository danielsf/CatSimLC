import numpy as np
from scipy.spatial import KDTree

chisq_cutoff = 700.0
span_cutoff = 300.0

kepler_dtype = np.dtype([('id', int), ('ra', float), ('dec', float), ('u', float),
                         ('g', float), ('r', float), ('i', float), ('z', float),
                         ('dist', float), ('teff_kic', float), ('teff_stellar', float),
                         ('Av', float), ('ebv', float)])
kep_file = 'KIC/kic_data.txt'

kep_data = np.genfromtxt(kep_file, dtype=kepler_dtype)
valid = np.where(kep_data['dist']>0.0)
kep_data = kep_data[valid]

rms_dtype = np.dtype([('name', str, 100), ('rms', float), ('chisq', float),
                      ('span', float), ('min_flux', float)])

rms_data = np.genfromtxt('rms_variability_lookup.txt', dtype=rms_dtype)
valid = np.where(np.logical_and(rms_data['min_flux']>0.0,
                 np.logical_and(rms_data['chisq']<chisq_cutoff,
                                rms_data['span']>span_cutoff)))

rms_data = rms_data[valid]
rms_dict = dict([(int(nn.split('_')[0][4:]),rr) for nn, rr in zip(rms_data['name'], rms_data['rms'])])

del rms_data

valid = []
for ix in range(len(kep_data)):
    if kep_data['id'][ix] in rms_dict:
        valid.append(ix)

valid = np.array(valid)
kep_data = kep_data[valid]
print len(kep_data)
exit()

# eqns 1 and 2 of
# Pisonneault et al 2012
# ApJS 199:30

new_mags = np.copy(kep_data)

new_mags['g'] = np.where(np.logical_and(kep_data['g']>0.0, kep_data['r']>0.0),
                 kep_data['g'] + 0.0921*(kep_data['g']-kep_data['r']) - 0.0985,
                 -999.0)

new_mags['r'] = np.where(np.logical_and(kep_data['r']>0.0, kep_data['i']>0.0),
                 kep_data['r'] + 0.0548*(kep_data['r']-kep_data['i']) - 0.0383,
                 -999.0)

new_mags['i'] = np.where(np.logical_and(kep_data['i']>0.0, kep_data['z']>0.0),
                 kep_data['i'] +0.0696*(kep_data['r']-kep_data['i']) - 0.0583,
                 -999.0)

new_mags['z'] = np.where(np.logical_and(kep_data['i']>0.0, kep_data['z']>0.0),
                 kep_data['z']+0.1587*(kep_data['i']-kep_data['z']) - 0.0597,
                 -999.0)

kep_abs_r = new_mags['r'] - 5.0*np.log10(new_mags['dist']/10.0)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

star_dtype = np.dtype([('id', long), ('u', float), ('g', float),
                       ('r', float), ('i', float), ('z', float),
                       ('parallax', float)])

star_data = np.genfromtxt(os.path.join('data','stars_to_match.txt'), dtype=star_dtype)

from lsst.sims.utils import radiansFromArcsec
_au_to_parsec = 1.0/206265.0
star_distance = _au_to_parsec/radiansFromArcsec(star_data['parallax']*0.001)
star_abs_r = star_data['r']-5.0*np.log10(star_distance/10.0)

plt.figsize = (30,30)

i_fig = 0

for mag1 in ('g', 'r', 'i'):
    for mag2 in ('r', 'i', 'z'):
        kep_color = new_mags[mag1]-new_mags[mag2]
        kep_params = np.array([kep_abs_r, kep_color]).transpose()
        kd_tree = KDTree(kep_params, leafsize=20)
        star_color = star_data[mag1]-star_data[mag2]
        star_params = np.array([star_abs_r, star_color]).transpose()
        match_dist, match_dex = kd_tree.query(star_params, k=1)
        i_fig += 1
        print mag1,mag2,len(np.unique(match_dex))
