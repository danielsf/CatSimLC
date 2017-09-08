import numpy as np
import os
from scipy.spatial import KDTree
import gc

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
rms_dict = dict([(int(nn.split('_')[0][4:]),ix) for ix, nn in enumerate(rms_data['name'])])

valid_kep = []
valid_rms = []
for ix in range(len(kep_data)):
    if kep_data['id'][ix] in rms_dict:
        valid_kep.append(ix)
        valid_rms.append(rms_dict[kep_data['id'][ix]])

valid_kep = np.array(valid_kep)
kep_data = kep_data[valid_kep]

valid_rms = np.array(valid_rms)
rms_data = rms_data[valid_rms]

print len(kep_data), len(rms_data)

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

header_list = []
label_list = []

import sys
sys.setrecursionlimit(100000)

for mag1, mag2 in zip(('g', 'r', 'i'), ('r', 'i', 'z')):
    with open(os.path.join('data', 'stars_to_match_%s_%s_association.txt' % (mag1,mag2)), 'w') as out_file:
        out_file.write('# catsim_id r_abs %s %s lc_name dex rms\n' % (mag1, mag2))
        kep_color = new_mags[mag1]-new_mags[mag2]

        kep_params = np.array([kep_abs_r, kep_color]).transpose()
        kd_tree = KDTree(kep_params, leafsize=1)
        star_color = star_data[mag1]-star_data[mag2]
        star_params = np.array([star_abs_r, star_color]).transpose()
        print 'doing search'
        match_dist, match_dex = kd_tree.query(star_params, k=1)

        for i_star in range(len(star_params)):
            out_file.write('%d %e %e %e %s %d %e\n' % (star_data['id'][i_star],star_abs_r[i_star],
                                                       star_data[mag1][i_star], star_data[mag2][i_star],
                                                       rms_data['name'][match_dex[i_star]],
                                                       match_dex[i_star],
                                                       rms_data['rms'][match_dex[i_star]]))

        i_fig += 1
        unq, unq_cts = np.unique(match_dex, return_counts=True)
        n_unq = len(unq)
        sorted_cts = np.sort(unq_cts)
        sum_cts = 0
        first_quartile = None
        third_quartile = None
        second_quartile = None
        for val in sorted_cts:
            sum_cts += val
            if sum_cts > 3*len(star_data)/4 and third_quartile is None:
                third_quartile = val
            if sum_cts > len(star_data)/2 and second_quartile is None:
                second_quartile = val
            if sum_cts > len(star_data)/4 and first_quartile is None:
                first_quartile = val
        print mag1,mag2,' -- ',n_unq,first_quartile,second_quartile,third_quartile
        actual_rms = rms_data['rms'][match_dex]

        sorted_rms = np.sort(actual_rms)
        cumulative_distribution = np.array([float(len(sorted_rms)-ix)/len(sorted_rms)
                                            for ix in range(len(sorted_rms))])

        hh, = plt.plot(sorted_rms*1000.0, cumulative_distribution)

        header_list.append(hh)
        label_list.append('%s-%s; %d %d %d %d' %
                         (mag1, mag2, n_unq, first_quartile, second_quartile, third_quartile))



plt.xscale('log')
plt.yscale('log')
plt.axvline(1.0,color='r',linestyle='--')
plt.axvline(10.0,color='r',linestyle='--')
plt.axvline(100.0,color='r',linestyle='--')
plt.axhline(0.1,color='r',linestyle='--')
plt.axhline(0.01,color='r',linestyle='--')
plt.ylim(1.0e-6,1.0)
plt.xlim(0.1,1000)
plt.xlabel('rms variability (mmag)')
plt.ylabel('cumulative distribution')
plt.legend(header_list, label_list, loc=0)
plt.title('%d stars' % len(star_data))
plt.savefig('variability_distribution.png')
plt.close()
