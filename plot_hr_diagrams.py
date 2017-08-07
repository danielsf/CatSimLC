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

# from An et al 2009
# ApJ 700, 523
ag_factor = 1.196
ar_factor = 0.874
ai_factor = 0.672
az_factor = 0.488

un_dereddened = np.copy(kep_data)
un_dereddened['g'] = new_mags['g']
un_dereddened['r'] = new_mags['r']
un_dereddened['i'] = new_mags['i']
un_dereddened['z'] = new_mags['z']

new_mags['g'] = np.where(new_mags['g']>0.0, new_mags['g'] - ag_factor*kep_data['Av'], -999.0)
new_mags['r'] = np.where(new_mags['r']>0.0, new_mags['r'] - ar_factor*kep_data['Av'], -999.0)
new_mags['i'] = np.where(new_mags['i']>0.0, new_mags['i'] - ai_factor*kep_data['Av'], -999.0)
new_mags['z'] = np.where(new_mags['z']>0.0, new_mags['z'] - az_factor*kep_data['Av'], -999.0)

kep_data['g'] = np.where(kep_data['g']>0.0, kep_data['g'] - ag_factor*kep_data['Av'], -999.0)
kep_data['r'] = np.where(kep_data['r']>0.0, kep_data['r'] - ar_factor*kep_data['Av'], -999.0)
kep_data['i'] = np.where(kep_data['i']>0.0, kep_data['i'] - ai_factor*kep_data['Av'], -999.0)
kep_data['z'] = np.where(kep_data['z']>0.0, kep_data['z'] - az_factor*kep_data['Av'], -999.0)


catsim_abs_r = catsim_data['r'] - 5.0*np.log10(catsim_data['dist']/10.0)
kep_abs_r = new_mags['r'] - 5.0*np.log10(new_mags['dist']/10.0)
un_dered_abs_r = un_dereddened['r'] - 5.0*np.log10(un_dereddened['dist']/10.0)

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

        valid = np.where(np.logical_and(un_dereddened['r']>0.0,
                         np.logical_and(un_dereddened[mag1]>0.0,
                                        un_dereddened[mag2]>0.0)))

        orig_r = un_dereddened['r'][valid] - 5.0*np.log10(un_dereddened['dist'][valid]/10.0)
        orig_color = un_dereddened[mag1][valid]-un_dereddened[mag2][valid]

        plt.subplot(3,2,i_fig)
        plt.xlim((color_min, color_max))
        plt.ylim((r_min, r_max))

        if i_fig == 1:
            plt.title('Blue is Kepler; Red is CatSim\nGreen is Kepler without dereddening')

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
plt.savefig('hr_diagaram_color.png')
plt.close()

plt.figsize = (30,30)
i_fig = 0
trim = 30
for imag1 in range(len(mag_list)):
    for imag2 in range(len(mag_list)):
        if imag2<=imag1:
            continue
        i_fig += 1
        plt.subplot(3,2,i_fig)
        mag1 = mag_list[imag1]
        mag2 = mag_list[imag2]

        catsim_color = catsim_data[mag1]-catsim_data[mag2]

        valid = np.where(np.logical_and(new_mags[mag1]>0.0,
                                        new_mags[mag2]>0.0))

        kep_color = new_mags[mag1][valid]-new_mags[mag2][valid]

        kep_color_sorted = np.sort(kep_color)
        n_kep = len(kep_color)
        xmin = kep_color_sorted[n_kep/trim]
        xmax = kep_color_sorted[(trim-1)*n_kep/trim]
        plt.xlim(xmin,xmax)

        valid = np.where(np.logical_and(un_dereddened[mag1]>0.0,
                                        un_dereddened[mag2]>0.0))
        un_dered_color = un_dereddened[mag1][valid]-un_dereddened[mag2][valid]

        if i_fig == 1:
            plt.title('Blue is Kepler; Red is CatSim;\nGreen is Kepler without dereddening',
                      fontsize=10)

        plt.hist(catsim_color, bins=1000, color='r', edgecolor='r', normed=True, zorder=1)
        plt.hist(kep_color, bins=1000, color='b', edgecolor='b', normed=True, zorder=2,
                 alpha=0.25)
        plt.hist(un_dered_color, bins=1000, color='g', edgecolor='g', normed=True, zorder=3,
                 alpha=0.25)
        plt.xlabel('%s-%s' % (mag1, mag2), fontsize=10)

plt.tight_layout()
plt.savefig('color_histograms.png')
plt.close()

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

        valid = np.where(np.logical_and(un_dereddened[mag]>0.0, un_dereddened[teff]>0.0))
        und_mag = un_dereddened[mag][valid] - 5.0*np.log10(un_dereddened['dist'][valid]/10.0)
        und_teff = un_dereddened[teff][valid]

        trim = 30
        n_kep = len(kep_mag)
        mag_sorted = np.sort(kep_mag)
        mag_min = mag_sorted[n_kep/trim]
        mag_max = mag_sorted[(trim-1)*n_kep/trim]
        teff_sorted = np.sort(kep_teff)
        teff_min = teff_sorted[n_kep/trim]
        teff_max = teff_sorted[(trim-1)*n_kep/trim]

        if i_fig == 1:
            plt.title('Blue is Kepler; Red is CatSim\nGreen is Kepler without dereddening')

        counts, xbins, ybins = np.histogram2d(catsim_teff, catsim_mag, bins=100)
        catsim = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                             colors='r', alpha=0.5, zorder=3)

        counts,xbins,ybins = np.histogram2d(kep_teff, kep_mag, bins=200)
        kep = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                          colors='blue', alpha=0.5, zorder=1)

        counts,xbins,ybins = np.histogram2d(und_teff, und_mag, bins=200)
        kep = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                          colors='green', alpha=0.5, zorder=2)

        plt.xlim((teff_min, teff_max))
        plt.ylim((mag_min, mag_max))

        plt.ylabel(mag)
        plt.xlabel(teff)
        plt.gca().invert_yaxis()
        plt.gca().invert_xaxis()

    plt.tight_layout()
    plt.savefig('hr_diagaram_%s.png' % teff)
    plt.close()


s2_coeffs = {'u':-0.249, 'g':0.794, 'r':-0.555, 'offset':0.234}
w2_coeffs = {'g':-0.227, 'r':0.792, 'i':-0.567, 'offset':0.05}
x2_coeffs = {'g':0.707, 'r':-0.707, 'offset':-0.988}
y2_coeffs = {'r':-0.270, 'i':0.8, 'z':-0.534, 'offset':0.054}

p2_coeff_dict = {'s':s2_coeffs, 'w':w2_coeffs, 'x':x2_coeffs, 'y':y2_coeffs}

s1_coeffs = {'u':0.91, 'g':-0.495, 'r':-0.415, 'offset':-1.28}
w1_coeffs = {'g':0.928, 'r':-0.556, 'i':-0.372, 'offset':-0.425}
x1_coeffs = {'r':1.0, 'i':-1.0, 'offset':0.0}
y1_coeffs = {'r':0.895, 'i':-0.448, 'z':-0.447, 'offset':-0.6}

p1_coeff_dict = {'s':s1_coeffs, 'w':w1_coeffs, 'x':x1_coeffs, 'y':y1_coeffs}

plt.figsize = (30,30)
trim = 30
for i_fig, color_name in enumerate(['w', 'x', 'y']):
    plt.subplot(2,2,i_fig+1)
    p2_coeffs = p2_coeff_dict[color_name]
    p1_coeffs = p1_coeff_dict[color_name]
    print 'kep_data ',len(kep_data)
    kep_valid = new_mags
    for tag in p2_coeffs.keys():
        if tag != 'offset':
            valid_dex = np.where(kep_valid[tag]>0.0)
            kep_valid = kep_valid[valid_dex]

    for tag in p1_coeffs.keys():
        if tag != 'offset':
            valid_dex = np.where(kep_valid[tag]>0.0)
            kep_valid = kep_valid[valid_dex]

    un_dered_valid = un_dereddened
    for tag in p2_coeffs.keys():
        if tag != 'offset':
            valid_dex = np.where(un_dered_valid[tag]>0.0)
            un_dered_valid = un_dered_valid[valid_dex]

    for tag in p1_coeffs.keys():
        if tag != 'offset':
            valid_dex = np.where(un_dered_valid[tag]>0.0)
            un_dered_valid = un_dered_valid[valid_dex]

    print 'becomes ',len(kep_valid)

    kep_p2 = np.ones(len(kep_valid))*p2_coeffs['offset']
    catsim_p2 = np.ones(len(catsim_data))*p2_coeffs['offset']
    un_dered_p2 = np.ones(len(un_dered_valid))*p2_coeffs['offset']
    for tag in p2_coeffs.keys():
        if tag == 'offset':
            continue
        kep_p2 += p2_coeffs[tag]*kep_valid[tag]
        catsim_p2 += p2_coeffs[tag]*catsim_data[tag]
        un_dered_p2 += p2_coeffs[tag]*un_dered_valid[tag]

    kep_p1 = np.ones(len(kep_valid))*p1_coeffs['offset']
    catsim_p1 = np.ones(len(catsim_data))*p1_coeffs['offset']
    un_dered_p1 = np.ones(len(un_dered_valid))*p1_coeffs['offset']
    for tag in p1_coeffs.keys():
        if tag == 'offset':
            continue
        kep_p1 += p1_coeffs[tag]*kep_valid[tag]
        catsim_p1 += p1_coeffs[tag]*catsim_data[tag]
        un_dered_p1 += p1_coeffs[tag]*un_dered_valid[tag]

    counts, xbins, ybins = np.histogram2d(catsim_p1, catsim_p2, bins=100)
    catsim = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                         colors='r',zorder=3)


    counts,xbins,ybins = np.histogram2d(kep_p1, kep_p2, bins=200)
    kep = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                      colors='blue',zorder=1)

    counts,xbins,ybins = np.histogram2d(un_dered_p1, un_dered_p2, bins=200)
    kep = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                      colors='green',zorder=2)


    if i_fig ==0:
        plt.title('blue is Kepler; red is CatSim;\ngreen is Kepler without dereddening', fontsize=12)

    plt.ylabel(color_name)
    plt.xlabel('P1(%s)' % color_name)
    #plt.xlim((min(kep_color.min(),catsim_color.min()),max(kep_color.max(),catsim_color.max())))

    kep_p1_sorted = np.sort(kep_p1)
    catsim_p1_sorted = np.sort(catsim_p1)

    xmin = min(kep_p1_sorted[len(kep_p1)/trim], catsim_p1_sorted[len(catsim_p1)/trim])
    xmax = max(kep_p1_sorted[(trim-1)*len(kep_p1)/trim], catsim_p1_sorted[(trim-1)*len(catsim_p1)/trim])

    kep_p2_sorted = np.sort(kep_p2)
    catsim_p2_sorted = np.sort(catsim_p2)

    ymin = min(kep_p2_sorted[len(kep_p2)/trim], catsim_p2_sorted[len(catsim_p2)/trim])
    ymax = max(kep_p2_sorted[(trim-1)*len(kep_p2)/trim], catsim_p2_sorted[(trim-1)*len(catsim_p2)/trim])

    if color_name == 'w':
        ymax = 0.06
        xmax = 0.8
    elif color_name == 'y':
        xmax = -0.1
    elif color_name == 'x':
        xmax = 0.4
        ymax = -0.3

    plt.ylim((ymin,ymax))
    plt.xlim((xmin,xmax))
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

plt.tight_layout()
plt.savefig('principal_colors.png')
plt.close()

