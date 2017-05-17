from __future__ import with_statement
import numpy as np
import os

def get_kurucz_phys(sed_name):
    """
    Read in the name of a kurucz SED file.  Return it's
    T_eff, metallicity, log(surface gravity)
    """
    if sed_name[1] == 'p':
        metallicity_sgn = 1.0
    elif sed_name[1] == 'm':
        metallicity_sgn = -1.0
    else:
        raise RuntimeError('Cannot parse metallicity sign of %s' % sed_name)

    new_name = sed_name.replace('.','_').split('_')

    metallicity = 0.1*metallicity_sgn*float(new_name[0][2:])

    teff = float(new_name[-2])

    logg = 0.1*np.float(new_name[3][1:])

    return teff, metallicity, logg


def get_wd_phys(sed_name):
    """
    Read in the name of a white dwarf SED,
    return its T_eff, metallicity (which we don't actually have),
    and log(surface gravity)
    """
    new_name = sed_name.replace('.','_').split('_')
    teff = float(new_name[-2])
    if new_name[1]!='He':
        logg = 0.1*float(new_name[2])
    else:
        logg = 0.1*float(new_name[3])

    return teff, -999.0, logg


def get_mlt_phys(sed_name):
    """
    Read in the name of an M/L/T dwarf SED and return
    its T_eff, metallicity, and log(surface gravity)
    """

    new_name = sed_name.replace('+','-').replace('a','-').split('-')

    logg_sgn_dex = len(new_name[0])

    if sed_name[logg_sgn_dex] == '-':
        logg_sgn = 1.0
    elif sed_name[logg_sgn_dex] == '+':
        logg_sgn = -1.0
    else:
        raise RuntimeError('Cannot get logg_sgn for %s' % sed_name)

    metallicity_sgn_dex = len(new_name[0]) + len(new_name[1]) + 1

    if sed_name[metallicity_sgn_dex] == '-':
        metallicity_sgn = -1.0
    elif sed_name[metallicity_sgn_dex] == '+':
        metallicity_sgn = 1.0
    else:
        raise RuntimeError('Cannot get metallicity_sgn for %s' % sed_name)

    teff = 100.0*float(new_name[0][3:])
    metallicity = metallicity_sgn*float(new_name[2])
    logg = logg_sgn*float(new_name[1])

    return teff, metallicity, logg


def get_physical_characteristics(sed_name):
    """
    Read in the name of an SED file.
    Return (in this order) Teff, metallicity (FeH), log(g)
    """
    sed_name = sed_name.strip()

    if not hasattr(get_physical_characteristics, 'teff_dict'):
        get_physical_characteristics.teff_dict = {}
        get_physical_characteristics.logg_dict = {}
        get_physical_characteristics.metal_dict = {}

    if sed_name in get_physical_characteristics.teff_dict:
        return (get_physical_characteristics.teff_dict[sed_name],
                get_physical_characteristics.metal_dict[sed_name],
                get_physical_characteristics.logg_dict[sed_name])

    if sed_name.startswith('bergeron'):
        sub_dir = 'wDs'
    elif sed_name.startswith('lte'):
        sub_dir = 'mlt'
    elif sed_name[0] == 'k':
        sub_dir = 'kurucz'
    else:
        raise RuntimeError("Do not understand name %s" % sed_name)



    if 'kurucz' in sub_dir:
        tt,mm,gg =  get_kurucz_phys(sed_name)
    elif sub_dir == 'wDs':
        tt,mm,gg = get_wd_phys(sed_name)
    elif sub_dir == 'mlt':
        tt,mm,gg = get_mlt_phys(sed_name)
    else:
        raise RuntimeError('Do not know how to get '
                           'physical characteristics for '
                           'sub_dir %s' % sub_dir)

    get_physical_characteristics.teff_dict[sed_name] = tt
    get_physical_characteristics.metal_dict[sed_name] = mm
    get_physical_characteristics.logg_dict[sed_name] = gg

    return tt, mm, gg

#need to create lookup table of temp, fe/h, logg
#and see how those things correlate when we just
#match on colors

from lsst.utils import getPackageDir

sed_dir = getPackageDir('sims_sed_library')
kurucz_dir = os.path.join(sed_dir, 'starSED', 'kurucz')

list_of_sed_files = os.listdir(kurucz_dir)

logg_dict = {}
teff_dict = {}
metallicity_dict = {}

for file_name in list_of_sed_files:
    tt,ff,ll = get_physical_characteristics(file_name)
    key_name = file_name.replace('.txt','').replace('.gz','')
    logg_dict[key_name] = ll
    teff_dict[key_name] = tt
    metallicity_dict[key_name] = ff

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
                      np.logical_and(np.logical_not(np.isnan(kep_data['jmag'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['hmag'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['kmag'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['dist'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['dist_err1'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['dist_err2'])),
                                     kep_data['dist']>1.0e-3))))))))))))))

kep_data = kep_data[valid_teff]

kep_j_h = kep_data['jmag']-kep_data['hmag']
kep_h_k = kep_data['hmag']-kep_data['kmag']

kep_colors = np.array([kep_j_h, kep_h_k]).transpose()


teff_up = kep_data['teff'] + kep_data['teff_err1']
teff_down = kep_data['teff'] + kep_data['teff_err2']

assert kep_data['teff_err2'].max()<0.0

dtemp = np.median(teff_up-teff_down)
assert dtemp>0.0
assert (teff_up-teff_down).min()>0.0

logg_up = kep_data['logg'] + kep_data['logg_err1']
logg_down = kep_data['logg'] + kep_data['logg_err2']

assert kep_data['logg_err2'].max()<0.0

dg = np.median(logg_up-logg_down)
assert dg>0.0
assert (logg_up-logg_down).min()>0.0

abs_mag = kep_data['kepmag']-5.0*np.log10(kep_data['dist']/10.0)

try:
    assert kep_data['dist_err2'].max()<0.0
except:
    print kep_data['dist_err2'].min(),kep_data['dist_err2'].max()
    raise

abs_up = kep_data['kepmag']-5.0*np.log10((kep_data['dist']+kep_data['dist_err2'])/10.0)
abs_down = kep_data['kepmag']-5.0*np.log10((kep_data['dist']+kep_data['dist_err1'])/10.0)

dmag = np.median(abs_up-abs_down)
assert dmag>0.0
assert (abs_up-abs_down).min()>0.0


with open('abs_mag_check.txt', 'w') as output_file:
    for ix in range(len(kep_data)):
        if kep_data['dist'][ix]>1.0e-4:
            output_file.write('%e %e %e %e %e\n' % (kep_data['kepmag'][ix], abs_mag[ix],
                                          kep_data['dist'][ix],
                                          kep_data['teff'][ix],kep_data['logg'][ix]))

from scipy.spatial import KDTree

print 'dmag ',dmag,abs_mag.max(),abs_mag.min(),kep_data['kepmag'].max(),kep_data['kepmag'].min()
print 'n_kep ',len(kep_data)

kep_params = np.array([kep_data['teff']/dtemp, kep_data['logg']/dg, abs_mag/dmag]).transpose()

print kep_params

import sys

sys.setrecursionlimit(10000)
kep_param_kdtree = KDTree(kep_params)


from lsst.sims.photUtils import cache_LSST_seds

sed_wavelen_min = 290.0
sed_wavelen_max = 2600.0

cache_LSST_seds(wavelen_min=sed_wavelen_min, wavelen_max=sed_wavelen_max)

from lsst.sims.photUtils import SedList, BandpassDict, Bandpass

throughput_dir = getPackageDir('throughputs')
twomass_dir = os.path.join(throughput_dir, '2MASS')

bp_list = []
bp_name_list = []

dtype = np.dtype([('w', float), ('s', float)])

bp_wavelen_max = 2500.0
bp_wavelen_min = 300.0

for file_name, bp_name in zip((os.path.join(twomass_dir,'2MASS_J.dat'),
                               os.path.join(twomass_dir,'2MASS_H.dat'),
                               os.path.join(twomass_dir,'2MASS_Ks.dat'),
                               'data/kepler_throughput.txt'),
                              ('j', 'h', 'k', 'kep')):

    data = np.genfromtxt(file_name, dtype=dtype)
    w_pad = np.arange(bp_wavelen_min, data['w'][0], 1.0)
    s_pad = np.zeros(len(w_pad))

    sb = np.append(s_pad, data['s'])
    wavelen = np.append(w_pad, data['w'])

    w_pad = np.arange(wavelen[-1], bp_wavelen_max, 1.0)
    s_pad = np.zeros(len(w_pad))

    sb = np.append(sb, s_pad)
    wavelen = np.append(wavelen, w_pad)

    bp_wavelen = np.arange(bp_wavelen_min, bp_wavelen_max, 1.0)
    bp_sb = np.interp(bp_wavelen, wavelen, sb)

    bp = Bandpass(wavelen=bp_wavelen, sb=bp_sb)

    bp_list.append(bp)
    bp_name_list.append(bp_name)



bp_dict = BandpassDict(bp_list, bp_name_list)

from lsst.sims.catalogs.db import DBObject

db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
              port=1433, driver='mssql+pymssql')

table_name = 'stars_obafgk_part_0870'
query = 'SELECT TOP 3000 sedfilename, flux_scale, ebv, parallax FROM %s' % table_name

query_dtype = np.dtype([('sedfilename', str, 200), ('scale', float),
                        ('ebv', float), ('parallax', float)])

star_iter = db.get_arbitrary_chunk_iterator(query, dtype=query_dtype,
                                            chunk_size=1000)


sed_list = None

out_name = 'test_star_fits.txt'

import time

with open(out_name, 'w') as output_file:
    output_file.write('# teff, logg, absmag, J-H, H-K\n')


from lsst.sims.utils import radiansFromArcsec

t_start = time.time()
ct=0
_au_to_parsec = 1.0/206265.0
for chunk in star_iter:
    ct += len(chunk)
    av = 3.1*chunk['ebv']
    magNorm = -2.5*np.log(chunk['scale'])/np.log(10.0) - 18.402732642

    if sed_list is None:
        sed_list = SedList(chunk['sedfilename'], magNorm,
                           galacticAvList=av, wavelenMatch=bp_dict.wavelenMatch)
    else:
        sed_list.flush()
        sed_list.loadSedsFromList(chunk['sedfilename'], magNorm,
                                  galacticAvList=av)


    mag_list = bp_dict.magListForSedList(sed_list)
    mag_list = mag_list.transpose()
    color_list = np.array([mag_list[0]-mag_list[1], mag_list[1]-mag_list[2]])

    #color_dist, color_dex = kep_kdtree.query(color_list)

    teff = []
    logg = []
    for name in chunk['sedfilename']:
        name = name.strip().replace('.txt', '').replace('.gz','')
        teff.append(teff_dict[name])
        logg.append(logg_dict[name])
    teff = np.array(teff)
    logg = np.array(logg)

    catsim_dist = _au_to_parsec/radiansFromArcsec(0.001*chunk['parallax'])

    catsim_abs_mag = mag_list[3]-5.0*np.log10(catsim_dist/10.0)

    pts = np.array([teff/dtemp, logg/dg, catsim_abs_mag/dmag]).transpose()
    param_dist, param_dex = kep_param_kdtree.query(pts)

    print '    mean dist ',np.mean(param_dist),' median dist ',np.median(param_dist),len(np.unique(param_dex))

    with open(out_name, 'a') as output_file:
        for ix, (name, dx) in enumerate(zip(chunk['sedfilename'], param_dex)):
            name = name.strip().replace('.txt','').replace('.gz','')
            output_file.write('%e %e %e %e %e %e %e %e %e %e\n' %
                              (kep_data['teff'][dx],teff_dict[name],
                               kep_data['logg'][dx], logg_dict[name],
                               abs_mag[dx], catsim_abs_mag[ix],
                               kep_j_h[dx], color_list[0][ix],
                               kep_h_k[dx], color_list[1][ix]))


    print 'did %d in %e ' % (ct, time.time()-t_start)

print 'data points ',len(kep_data)
print 'dtemp ',dtemp
print 'dg ',dg
print 'dmag ',dmag
