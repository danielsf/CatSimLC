import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

catsim_dtype = np.dtype([('sedname', str, 300), ('kep_abs', float),
                         ('u', float), ('g', float), ('r', float), ('i', float), ('z', float),
                         ('dist', float), ('teff', float), ('r_abs', float)])
catsim_file = 'catsim_star_data_same_pointing.txt'
catsim_data = np.genfromtxt(catsim_file, dtype=catsim_dtype)

kepler_dtype = np.dtype([('id', int), ('ra', float), ('dec', float), ('u', float),
                         ('g', float), ('r', float), ('i', float), ('z', float),
                         ('dist', float), ('teff_kic', float), ('teff_stellar', float),
                         ('Av', float), ('ebv', float)])
kep_file = 'KIC/kic_data.txt'

kep_data = np.genfromtxt(kep_file, dtype=kepler_dtype)

valid = np.where(kep_data['Av']>-990.0)
kep_data = kep_data[valid]

valid = np.where(kep_data['dist']>0.0)
kep_data = kep_data[valid]

# eqns 1 and 2 of
# Pisonneault et al 2012
# ApJS 199:30

new_mags = {}

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

# from An et al 2009
# ApJ 700, 523
ag_factor = 1.196
ar_factor = 0.874
ai_factor = 0.672
az_factor = 0.488

new_mags['g'] = np.where(new_mags['g']>0.0, new_mags['g'] - ag_factor*kep_data['Av'], -999.0)
new_mags['r'] = np.where(new_mags['r']>0.0, new_mags['r'] - ar_factor*kep_data['Av'], -999.0)
new_mags['i'] = np.where(new_mags['i']>0.0, new_mags['i'] - ai_factor*kep_data['Av'], -999.0)
new_mags['z'] = np.where(new_mags['z']>0.0, new_mags['z'] - az_factor*kep_data['Av'], -999.0)

kep_data['g'] = np.where(kep_data['g']>0.0, kep_data['g'] - ag_factor*kep_data['Av'], -999.0)
kep_data['r'] = np.where(kep_data['r']>0.0, kep_data['r'] - ar_factor*kep_data['Av'], -999.0)
kep_data['i'] = np.where(kep_data['i']>0.0, kep_data['i'] - ai_factor*kep_data['Av'], -999.0)
kep_data['z'] = np.where(kep_data['z']>0.0, kep_data['z'] - az_factor*kep_data['Av'], -999.0)


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

        kep_r = new_mags['r'][valid] - 5.0*np.log10(kep_data['dist'][valid]/10.0)
        kep_color = new_mags[mag1][valid] - new_mags[mag2][valid]

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

        valid = np.where(np.logical_and(kep_data['r']>0.0,
                         np.logical_and(kep_data[mag1]>0.0,
                                        kep_data[mag2]>0.0)))

        orig_r = kep_data['r'][valid] - 5.0*np.log10(kep_data['dist'][valid]/10.0)
        orig_color = kep_data[mag1][valid]-kep_data[mag2][valid]

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

        counts,xbins,ybins = np.histogram2d(orig_color, orig_r, bins=200)
        kep = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                          colors='green', alpha=0.5, zorder=2)

        plt.ylabel('r')
        plt.xlabel('%s-%s' % (mag1,mag2))
        plt.gca().invert_yaxis()

plt.tight_layout()
plt.savefig('hr_diagaram_dered.png')
plt.close()

mag_list = ['g', 'r', 'i', 'z']
teff_list = ['teff_kic', 'teff_stellar']

catsim_teff = catsim_data['teff']

for teff in teff_list:
    plt.figsize = (30, 30)
    i_fig = 0
    for mag in mag_list:
        i_fig += 1
        plt.subplot(2,2,i_fig)

        catsim_mag = catsim_data[mag] - 5.0*np.log10(catsim_data['dist']/10.0)

        valid = np.where(np.logical_and(new_mags[mag]>0.0, kep_data[teff]>0.0))

        kep_mag = new_mags[mag][valid] - 5.0*np.log10(kep_data['dist'][valid]/10.0)
        kep_teff = kep_data[teff][valid]

        trim = 30
        n_kep = len(kep_mag)
        mag_sorted = np.sort(kep_mag)
        mag_min = mag_sorted[n_kep/trim]
        mag_max = mag_sorted[(trim-1)*n_kep/trim]
        teff_sorted = np.sort(kep_teff)
        teff_min = teff_sorted[n_kep/trim]
        teff_max = teff_sorted[(trim-1)*n_kep/trim]

        if i_fig == 1:
            plt.title('Blue is Kepler; Red is CatSim')

        counts, xbins, ybins = np.histogram2d(catsim_teff, catsim_mag, bins=100)
        catsim = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                             colors='r', alpha=0.5, zorder=3)

        counts,xbins,ybins = np.histogram2d(kep_teff, kep_mag, bins=200)
        kep = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                          colors='blue', alpha=0.5, zorder=1)

        plt.xlim((teff_min, teff_max))
        plt.ylim((mag_min, mag_max))

        plt.ylabel(mag)
        plt.xlabel(teff)
        plt.gca().invert_yaxis()
        plt.gca().invert_xaxis()

    plt.tight_layout()
    plt.savefig('hr_diagaram_%s.png' % teff)
    plt.close()
