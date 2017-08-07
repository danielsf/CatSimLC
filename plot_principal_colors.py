import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np


kepler_dtype = np.dtype([('id', int), ('ra', float), ('dec', float),
                         ('u', float), ('g', float), ('r', float),
                         ('i', float), ('z', float), ('dist', float),
                         ('teff', float), ('teff_stellar', float),
                         ('Av', float), ('ebv', float)])

transformed_file = 'KIC/kic_data_transformed.txt'
kic_file = 'KIC/kic_data.txt'

transformed_data = np.genfromtxt(transformed_file, dtype=kepler_dtype)
kic_data = np.genfromtxt(kic_file, dtype=kepler_dtype)


ag_factor = 1.196
ar_factor = 0.874
ai_factor = 0.672
az_factor = 0.488

trans_valid = np.where(transformed_data['Av']>-990.0)
transformed_data = transformed_data[trans_valid]


print 'max_g ',np.abs(ag_factor*transformed_data['Av']).max()
print 'max_r ',np.abs(ar_factor*transformed_data['Av']).max()
print 'extreme Av', transformed_data['Av'].min(),transformed_data['Av'].max()


transformed_data['g'] -= ag_factor*transformed_data['Av']
transformed_data['r'] -= ar_factor*transformed_data['Av']
transformed_data['i'] -= ai_factor*transformed_data['Av']
transformed_data['z'] -= az_factor*transformed_data['Av']

catsim_dtype = np.dtype([('sedname', str, 300), ('teff', float), ('feh', float), ('logg', float), ('m', float),
                   ('r_abs', float), ('norm', float),
                   ('u', float), ('g', float), ('r', float), ('i', float), ('z', float)])

catsim_file = 'catsim_star_data_same_pointing_cutoff.txt'

catsim_data = np.genfromtxt(catsim_file, dtype=catsim_dtype)


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
trim = 50
for i_fig, color_name in enumerate(['s', 'w', 'x', 'y']):
    plt.subplot(2,2,i_fig+1)
    p2_coeffs = p2_coeff_dict[color_name]
    p1_coeffs = p1_coeff_dict[color_name]
    print 'kic_data ',len(kic_data)
    kic_valid = kic_data
    for tag in p2_coeffs.keys():
        if tag != 'offset':
            valid_dex = np.where(kic_valid[tag]>0.0)
            kic_valid = kic_valid[valid_dex]

    for tag in p1_coeffs.keys():
        if tag != 'offset':
            valid_dex = np.where(kic_valid[tag]>0.0)
            kic_valid = kic_valid[valid_dex]

    print 'becomes ',len(kic_valid)

    trans_valid = transformed_data
    for tag in p2_coeffs.keys():
        if tag != 'offset':
            valid_dex = np.where(trans_valid[tag]>0.0)
            trans_valid = trans_valid[valid_dex]
    for tag in p1_coeffs.keys():
        if tag != 'offset':
            valid_dex = np.where(trans_valid[tag]>0.0)
            trans_valid = trans_valid[valid_dex]

    trans_p2 = np.ones(len(trans_valid))*p2_coeffs['offset']
    kep_p2 = np.ones(len(kic_valid))*p2_coeffs['offset']
    catsim_p2 = np.ones(len(catsim_data))*p2_coeffs['offset']
    for tag in p2_coeffs.keys():
        if tag == 'offset':
            continue
        trans_p2 += p2_coeffs[tag]*trans_valid[tag]
        kep_p2 += p2_coeffs[tag]*kic_valid[tag]
        catsim_p2 += p2_coeffs[tag]*catsim_data[tag]

    trans_p1 = np.ones(len(trans_valid))*p1_coeffs['offset']
    kep_p1 = np.ones(len(kic_valid))*p1_coeffs['offset']
    catsim_p1 = np.ones(len(catsim_data))*p1_coeffs['offset']
    for tag in p1_coeffs.keys():
        if tag == 'offset':
            continue
        trans_p1 += p1_coeffs[tag]*trans_valid[tag]
        kep_p1 += p1_coeffs[tag]*kic_valid[tag]
        catsim_p1 += p1_coeffs[tag]*catsim_data[tag]

    counts, xbins, ybins = np.histogram2d(catsim_p1, catsim_p2, bins=100)
    catsim = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                         colors='r',zorder=3, alpha=0.5)


    counts,xbins,ybins = np.histogram2d(kep_p1, kep_p2, bins=200)
    kep = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                      colors='blue',zorder=1)

    counts,xbins,ybins = np.histogram2d(trans_p1, trans_p2, bins=200)
    trans = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                      colors='green',zorder=2)

    if i_fig ==0:
        plt.title('blue is KIC; green is Pinsonneault; is CatSim')

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
plt.savefig('kepler_principal_colors_dereddened_2d.png')
plt.close()

