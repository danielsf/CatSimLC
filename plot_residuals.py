import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import os

def make_2d_histogram(xx, yy, dx, dy):
    """
    returns indices and counts of unique points on the map
    """
    i_color1 = np.round(xx/dx).astype(int)
    i_color2 = np.round(yy/dy).astype(int)
    dex_reverse = np.array([i_color1, i_color2])
    dex_arr = dex_reverse.transpose()
    # see http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
    dex_raw = np.ascontiguousarray(dex_arr).view(np.dtype((np.void, dex_arr.dtype.itemsize*dex_arr.shape[1])))
    _, unique_rows, unique_counts = np.unique(dex_raw, return_index=True, return_counts=True)

    return unique_rows, unique_counts


def plot_color(xx, yy, dx, dy):
    dexes, cts = make_2d_histogram(xx, yy, dx, dy)
    sorted_dex = np.argsort(cts)
    dexes = dexes[sorted_dex]
    cts = cts[sorted_dex]
    plt.scatter(xx[dexes], yy[dexes], c=cts, s=5,
                cmap=plt.cm.gist_ncar, edgecolor='')

    plt.colorbar()


kep_dtype = np.dtype([('kepid', int),
                      ('tm_designation', str, 200),
                      ('teff', float),
                      ('teff_err1', float),
                      ('teff_err2', float),
                      ('logg', float),
                      ('logg_err1', float),
                      ('logg_err2', float),
                      ('feh', float),
                      ('feh_err1', float),
                      ('feh_err2', float),
                      ('mass', float),
                      ('mass_err1', float),
                      ('mass_err2', float),
                      ('st_radius', float),
                      ('radius_err1', float),
                      ('radius_err2', float),
                      ('dens', float),
                      ('dens_err1', float),
                      ('dens_err2', float),
                      ('prov_sec', str, 200),
                      ('kepmag', float),
                      ('dist', float),
                      ('dist_err1', float),
                      ('dist_err2', float),
                      ('nconfp', int),
                      ('nkoi', int),
                      ('ntce', int),
                      ('datalink_dvr', str, 200),
                      ('st_delivname', str, 200),
                      ('st_vet_date_str', str, 200),
                      ('degree_ra', float),
                      ('degree_dec', float),
                      ('st_quarters', int),
                      ('teff_prov', str, 200),
                      ('logg_prov', str, 200),
                      ('feh_prov', str, 200),
                      ('jmag', float),
                      ('jmag_err', float),
                      ('hmag', float),
                      ('hmag_err', float),
                      ('kmag', float),
                      ('kmag_err', float),
                      ('dutycycle', float),
                      ('dataspan', float),
                      ('mesthres01p5', float),
                      ('mesthres02p0', float),
                      ('mesthres02p5', float),
                      ('mesthres03p0', float),
                      ('mesthres03p5', float),
                      ('mesthres04p5', float),
                      ('mesthres05p0', float),
                      ('mesthres06p0', float),
                      ('mesthres07p5', float),
                      ('mesthres09p0', float),
                      ('mesthres10p5', float),
                      ('mesthres12p0', float),
                      ('mesthres12p5', float),
                      ('mesthres15p0', float),
                      ('rrmscdpp01p5', float),
                      ('rrmscdpp02p0', float),
                      ('rrmscdpp02p5', float),
                      ('rrmscdpp03p0', float),
                      ('rrmscdpp03p5', float),
                      ('rrmscdpp04p5', float),
                      ('rrmscdpp05p0', float),
                      ('rrmscdpp06p0', float),
                      ('rrmscdpp07p5', float),
                      ('rrmscdpp09p0', float),
                      ('rrmscdpp10p5', float),
                      ('rrmscdpp12p0', float),
                      ('rrmscdpp12p5', float),
                      ('rrmscdpp15p0', float),
                      ('av', float),
                      ('av_err1', float),
                      ('av_err2', float)])

catalog_name = 'data/kepler_stellar17.csv'

if not os.path.exists(catalog_name):
   raise RuntimeError('Need to download Kepler catalog data using\n'
                      'wget -nH --cut-dirs=3 https://archive.stsci.edu/pub/kepler/catalogs/kepler_stellar17.csv.gz\n'
                      'and place it in the data/ directory')


kep_data = np.genfromtxt(catalog_name, dtype=kep_dtype,
                         delimiter='|', skip_header=1)

valid_teff = np.where(np.logical_and(np.logical_not(np.isnan(kep_data['teff'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['teff_err1'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['teff_err2'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['logg'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['logg_err1'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['logg_err2'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['kepmag'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['feh'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['feh_err1'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['feh_err2'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['dist'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['dist_err1'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['dist_err2'])),
                                     kep_data['dist']>1.0e-3))))))))))))))

kep_data = kep_data[valid_teff]

kep_abs_mag = kep_data['kepmag']-5.0*np.log10(kep_data['dist']/10.0)


dtype = np.dtype([('dex', int), ('kep_teff', float), ('catsim_teff', float),
                  ('kep_feh', float), ('catsim_feh', float),
                  ('kep_mag', float), ('catsim_mag', float)])

for kk in (1, 10):
    data_file = 'test_star_fits_k%d.txt' % kk
    data = np.genfromtxt(data_file, dtype=dtype)
    print 'read in data'

    used_dexes = np.unique(data['dex'])
    unused_dexes = [ix for ix in range(len(kep_data)) if ix not in used_dexes]
    unused_dexes = np.array(unused_dexes)
    print ('used ',len(used_dexes),' unused ',len(unused_dexes),' available ',
           len(kep_data))

    plt.figsize = (30, 30)
    plt.subplot(3,2,1)
    plot_color(kep_data['teff'][unused_dexes], kep_abs_mag[unused_dexes],
               100.0, 0.1)
    t_min = kep_data['teff'].min()
    t_max = kep_data['teff'].max()
    xticks = np.arange(np.round(t_min/1000.0)*1000.0,
                       np.round(t_max/1000.0)*1000.0,
                       2000.0)

    xlabels = ['%d' % xx for xx in xticks]
    plt.xlabel('Teff', fontsize=7)
    plt.ylabel('Absolute Kepler Magnitude', fontsize=7)
    plt.xticks(xticks, xlabels, fontsize=7)
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()

    plt.tight_layout()
    plt.savefig('fit_plots_nn%d.png' % kk)
    plt.close()
