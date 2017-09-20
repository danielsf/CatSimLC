import numpy as np
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


#stellar_temp_dict = {}
#with open(os.path.join('data', 'teff_to_type.txt'), 'r') as input_file:
#    for line in input_file:
#        if line[0] == '#':
#            continue
#        params = line.strip().split()
#        stellar_temp_dict[params[0]] = float(params[1])


print "built stellar temp dict"

dtype = np.dtype([('id', long), ('r_abs', float), ('g', float), ('r', float),
                   ('lc', str, 100), ('dex', int), ('rms', float),
                   ('sed', str, 200), ('teff', float), ('logg', float)])

star_data = np.genfromtxt(os.path.join('data', 'stars_to_match_g_r_association.txt'),
                          dtype=dtype)


#star_data = np.genfromtxt(os.path.join('data', 'junk.txt'),
#                          dtype=dtype)


valid = np.where(np.logical_and(star_data['teff']>-998.0, star_data['logg']>-998.0))

star_data = star_data[valid]

teff_arr = star_data['teff']
logg_arr = star_data['logg']
rms_arr = star_data['rms']

del star_data
import gc
gc.collect()
print "read in data ",len(teff_arr),len(logg_arr),len(rms_arr)

stellar_type_list = []

stellar_type_list.append(('K2', 'K5', 3.75, -0.5, -1.1, None, None, None, None, 0.4, 1.5, 4410.0, 4830.0))
stellar_type_list.append(('K0', 'K2', 4.0, -0.75, -1.2, None, None, None, None, 0.55, 1.5, 4830.0, 5150.0))
stellar_type_list.append(('G8', 'K0', 4.1,-0.6, -0.9, None, None, None, None, 0.55, 1.5, 5150.0, 5310.0))
stellar_type_list.append(('G5', 'G8', 4.2, -0.1, -1.15, None, None, -0.2, -1.4, None, None, 5310.0, 5560.0))
stellar_type_list.append(('G2', 'G5', 4.3, -0.65, -0.9, None, None, None, None, 1.5, 1.45, 5560.0, 5790.0))
stellar_type_list.append(('G0', 'G2', 4.3, -0.75, -0.95, None, None, -0.65, -1.2, None, None, 5790.0, 5940.0))
stellar_type_list.append(('F8', 'G0', 4.3, -0.75, -0.85, None, None, -0.75, -1.2, None, None, 5940.0, 6250.0))
stellar_type_list.append(('F5', 'F8', 4.25, -0.7, -0.8, None, None, -0.75, -1.0, None, None, 6250.0, 6650.0))
stellar_type_list.append(('F2', 'F5', 4.3, -0.5, -0.7, None, None, -0.6, -0.9, None, None, 6650.0, 7000.0))

empirical_sigma = np.arange(0.1,200.0,0.01)

plt.figsize = (30,30)
sorted_rms = np.sort(rms_arr)*1000.0
cumulative_dist = np.array([float(len(sorted_rms)-ix)/len(sorted_rms)
                            for ix in range(len(sorted_rms))])
hh, = plt.plot(sorted_rms, cumulative_dist, color='k')
plt.xscale('log')
plt.yscale('log')
plt.axvline(1.0,color='b',linestyle=':')
plt.axvline(10.0,color='b',linestyle=':')
plt.axvline(100.0,color='b',linestyle=':')
plt.axhline(0.1,color='b',linestyle=':')
plt.axhline(0.01,color='b',linestyle=':')
plt.axhline(0.001,color='b',linestyle=':')
plt.axhline(0.0001,color='b',linestyle=':')
plt.axhline(0.00001,color='b',linestyle=':')
plt.ylim(1.0e-6,1.0)
plt.xlim(0.1,1000)
plt.xlabel('rms variability (mmag)')
plt.ylabel('cumulative distribution')
plt.savefig('variability_distribution_total.png')
plt.close()

for plot_params in stellar_type_list:
    t1 = plot_params[0]
    t2 = plot_params[1]
    g_div = plot_params[2]
    t_max = plot_params[12]
    t_min = plot_params[11]
    valid = np.where(np.logical_and(teff_arr<t_max,
                                    teff_arr>=t_min))

    teff_t_valid = teff_arr[valid]
    logg_t_valid = logg_arr[valid]
    rms_t_valid = rms_arr[valid]

    valid = np.where(logg_t_valid<g_div)
    logg_g_valid = logg_t_valid[valid]
    teff_g_valid = teff_t_valid[valid]
    rms_g_valid = rms_t_valid[valid]

    sorted_rms = np.sort(rms_g_valid)*1000.0
    cumulative_dist = np.array([float(len(sorted_rms)-ix)/len(sorted_rms)
                                for ix in range(len(sorted_rms))])

    if plot_params[3] is None:
        empirical = 1.0/(1.0+plot_params[5]*np.power(empirical_sigma, plot_params[6]))
    else:
        empirical = np.power(10.0, plot_params[3])*np.power(empirical_sigma, plot_params[4])


    #emp_sum = empirical.sum()
    #cumulative_empirical=np.array([empirical[ix:].sum()/emp_sum
    #                               for ix in range(len(empirical))])

    header_list = []
    label_list = []
    plt.figsize=(30,30)
    hh, = plt.plot(sorted_rms, cumulative_dist, color='k')
    plt.plot(empirical_sigma, empirical, color='k', linestyle='--')
    header_list.append(hh)
    label_list.append('log(g)<%.2f' % g_div)

    valid = np.where(logg_t_valid>=g_div)
    logg_g_valid = logg_t_valid[valid]
    teff_g_valid = teff_t_valid[valid]
    rms_g_valid = rms_t_valid[valid]

    sorted_rms = np.sort(rms_g_valid)*1000.0
    cumulative_dist = np.array([float(len(sorted_rms)-ix)/len(sorted_rms)
                                for ix in range(len(sorted_rms))])

    if plot_params[7] is None:
        empirical = 1.0/(1.0+plot_params[9]*np.power(empirical_sigma, plot_params[10]))
    else:
        empirical = np.power(10.0, plot_params[7])*np.power(empirical_sigma, plot_params[8])

    #emp_sum = empirical.sum()
    #cumulative_empirical=np.array([empirical[ix:].sum()/emp_sum
    #                               for ix in range(len(empirical))])

    hh, = plt.plot(sorted_rms, cumulative_dist, color='r')
    plt.plot(empirical_sigma, empirical, color='r', linestyle='--')
    header_list.append(hh)
    label_list.append('log(g)>=%.2f' % g_div)

    plt.title('%s-%s stars; ct is %d' % (t1, t2, len(teff_t_valid)))
    plt.xscale('log')
    plt.yscale('log')
    plt.axvline(1.0,color='b',linestyle=':')
    plt.axvline(10.0,color='b',linestyle=':')
    plt.axvline(100.0,color='b',linestyle=':')
    plt.axhline(0.1,color='b',linestyle=':')
    plt.axhline(0.01,color='b',linestyle=':')
    plt.ylim(0.001,1.0)
    plt.xlim(0.1,100)
    plt.xlabel('rms variability (mmag)')
    plt.ylabel('cumulative distribution')
    plt.legend(header_list, label_list, loc=0)
    plt.savefig('variability_distribution_%s_%s.png' % (t1, t2))
    plt.close()


