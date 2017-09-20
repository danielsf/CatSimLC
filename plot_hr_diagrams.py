import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

catsim_dtype = np.dtype([('sedname', str, 300), ('kep_abs', float),
                         ('u', float), ('g', float), ('r', float), ('i', float),
                         ('z', float),
                         ('dist', float), ('teff', float), ('r_abs', float)])
catsim_file = 'catsim_star_data_same_pointing_sdss.txt'
catsim_data = np.genfromtxt(catsim_file, dtype=catsim_dtype)

kepler_dtype = np.dtype([('id', int), ('ra', float), ('dec', float), ('u', float),
                         ('g', float), ('r', float), ('i', float), ('z', float),
                         ('dist', float), ('teff_kic', float), ('teff_stellar', float),
                         ('Av', float), ('ebv', float)])
kep_file = 'KIC/kic_data.txt'

kep_data = np.genfromtxt(kep_file, dtype=kepler_dtype)
valid = np.where(kep_data['dist']>0.0)
kep_data = kep_data[valid]

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

catsim_abs_r = catsim_data['r'] - 5.0*np.log10(catsim_data['dist']/10.0)
kep_abs_r = new_mags['r'] - 5.0*np.log10(new_mags['dist']/10.0)

mag_list = ['g', 'r', 'i', 'z']
plt.figsize=(30,30)
i_fig = 0
for i_mag_1 in range(len(mag_list)):
    for i_mag_2 in range(len(mag_list)):
        if i_mag_2<=i_mag_1:
            continue
        i_fig += 1
        mag1 = mag_list[i_mag_1]
        mag2 = mag_list[i_mag_2]

        valid = np.where(np.logical_and(new_mags['r']>0.0,
                         np.logical_and(new_mags[mag1]>0.0,
                                        new_mags[mag2]>0.0)))

        print 'r ',new_mags['r'].min(),new_mags['r'].max()
        print mag1,new_mags[mag1].min(),new_mags[mag2].max()
        print mag2,new_mags[mag2].min(),new_mags[mag2].max()

        kep_r = new_mags['r'][valid] - 5.0*np.log10(kep_data['dist'][valid]/10.0)
        kep_color = new_mags[mag1][valid] - new_mags[mag2][valid]

        print 'color ',kep_color.min(),kep_color.max()
        print 'dist ',kep_data['dist'][valid].min(),kep_data['dist'][valid].max()

        n_kep = len(kep_r)
        trim = 30
        kep_r_sorted = np.sort(kep_r)
        r_min = kep_r_sorted[n_kep/trim]
        r_max = kep_r_sorted[(trim-1)*n_kep/trim]

        kep_color_sorted = np.sort(kep_color)
        color_min = kep_color_sorted[n_kep/trim]
        color_max = kep_color_sorted[(trim-1)*n_kep/trim]

        catsim_r = catsim_data['r'] - 5.0*np.log10(catsim_data['dist']/10.0)
        catsim_color = catsim_data[mag1] - catsim_data[mag2]

        plt.subplot(3,2,i_fig)
        plt.xlim((color_min, color_max))
        plt.ylim((r_min, r_max))

        if i_fig == 1:
            plt.title('Blue is Kepler; Red is CatSim')

        counts, xbins, ybins = np.histogram2d(catsim_color, catsim_r, bins=100)
        catsim = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                             colors='r', alpha=0.5, zorder=3)

        counts,xbins,ybins = np.histogram2d(kep_color, kep_r, bins=200)
        kep = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                          colors='blue', alpha=0.5, zorder=1)

        plt.ylabel('r')
        plt.xlabel('%s-%s' % (mag1,mag2))
        plt.gca().invert_yaxis()

plt.tight_layout()
plt.savefig('hr_diagaram_color_sdss.png')
plt.close()
