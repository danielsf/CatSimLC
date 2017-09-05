import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import os

from lsst.sims.utils import radiansFromArcsec

def make_grid(xx, yy, dx, dy):
    xdex = np.round(xx/dx).astype(int)
    ydex = np.round(yy/dy).astype(int)
    x_grid = []
    y_grid = []
    pairs = []
    for xd, yd in zip(xdex, ydex):
        rep = '%d_%d' % (xd,yd)
        if rep not in pairs:
            x_grid.append(xd)
            y_grid.append(yd)
            pairs.append(rep)

    return np.array(x_grid)*dx, np.array(y_grid)*dy

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

    teff = float(new_name[-1])

    logg = 0.1*np.float(new_name[3][1:])

    return teff, metallicity, logg


def get_wd_phys(sed_name):
    """
    Read in the name of a white dwarf SED,
    return its T_eff, metallicity (which we don't actually have),
    and log(surface gravity)
    """
    new_name = sed_name.replace('.','_').split('_')
    teff = float(new_name[-1])
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

kep_data = np.genfromtxt('data/kepler_stellar17.csv', dtype=kep_dtype,
                         delimiter='|', skip_header=1)

valid_teff = np.where(np.logical_and(np.logical_not(np.isnan(kep_data['teff'])),
                      np.logical_and(np.logical_not(np.isnan(kep_data['logg'])),
                                     np.logical_not(np.isnan(kep_data['kepmag'])))))

kep_data = kep_data[valid_teff]


kep_j_h = kep_data['jmag']-kep_data['hmag']
kep_h_k = kep_data['hmag']-kep_data['kmag']

print 'ra ',kep_data['degree_ra'].min(),kep_data['degree_ra'].max()
print 'dec ',kep_data['degree_dec'].min(),kep_data['degree_dec'].max()

ra_min = kep_data['degree_ra'].min()
ra_max = kep_data['degree_ra'].max()
dec_min = kep_data['degree_dec'].min()
dec_max = kep_data['degree_dec'].max()
mag_max = kep_data['kepmag'].max()

teff_catsim = []
logg_catsim = []

mag_grid_dict = {}
with open('data/lsst_color_to_kepler_grid.txt', 'r') as input_file:
    input_lines = input_file.readlines()
    for line in input_lines:
        if line[0] == '#':
            continue
        line = line.split()
        name = line[0].replace('.txt','').replace('.gz','')
        if name not in mag_grid_dict:
            mag_grid_dict[name] = {}
            mag_grid_dict[name]['magnorm'] = []
            mag_grid_dict[name]['kep'] = []
        mag_grid_dict[name]['magnorm'].append(float(line[1]))
        mag_grid_dict[name]['kep'].append(float(line[2]))

for sed_name in mag_grid_dict:
    mag_grid_dict[sed_name]['magnorm'] = np.array(mag_grid_dict[sed_name]['magnorm'])
    mag_grid_dict[sed_name]['kep'] = np.array(mag_grid_dict[sed_name]['kep'])

print 'built grid'



from lsst.sims.photUtils import Sed

dummy_sed = Sed()

ct=0
sed_list = None
import time
t_start = time.time()

percent_per_table = 0.05

_au_to_parsec = 1.0/206265.0

colnames = ['sedfilename', 'parallax', 'ebv', 'varParamStr', 'flux_scale']
for name in ('u', 'g', 'r', 'i', 'z'):
    colnames.append('sdss%s' % name)


out_name = 'catsim_star_data_same_pointing_sdss.txt'

with open(out_name, 'w') as out_file:
    out_file.write('# name abs_kep_mag ugriz ugriz_noatm dist Teff sdss_r_abs\n')

from lsst.sims.catUtils.baseCatalogModels import StarObj
db = StarObj(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
             port=1433, driver='mssql+pymssql')

from lsst.sims.utils import ObservationMetaData
obs = ObservationMetaData(pointingRA=0.5*(ra_min+ra_max),
                          pointingDec=0.5*(dec_min+dec_max),
                          boundType='circle',boundLength=max(0.5*(ra_max-ra_min),0.5*(dec_max-dec_min)))
                          #boundLength=(0.5*(ra_max-ra_min),
                          #             0.5*(dec_max-dec_min)))

star_iter = db.query_columns(colnames=colnames, obs_metadata=obs, chunk_size=10000,
                             constraint='sdssr<20.298')

for chunk in star_iter:

    valid = np.where(np.char.find(chunk['varParamStr'], 'None')==0)
    if len(valid[0])==0:
        print 'No valid rows'
        print chunk['varParamStr']
        continue

    chunk = chunk[valid]

    ct += len(chunk)
    magnorm = -2.5*np.log(chunk['flux_scale'])/np.log(10.0)-18.402732642

    catsim_kep_mag = []
    for name, mm in zip(chunk['sedfilename'], magnorm):
        kep_mag = np.interp(mm, mag_grid_dict[name]['magnorm'], mag_grid_dict[name]['kep'])
        catsim_kep_mag.append(kep_mag)
    catsim_kep_mag = np.array(catsim_kep_mag)

    #valid_mag = np.where(catsim_kep_mag<mag_max)
    #print 'valid_mag ',len(valid_mag[0])
    #chunk = chunk[valid_mag]
    #catsim_kep_mag = catsim_kep_mag[valid_mag]
    #magnorm_list = magnorm[valid_mag]


    catsim_dist = _au_to_parsec/chunk['parallax']
    catsim_abs_mag = catsim_kep_mag-5.0*np.log10(catsim_dist/10.0)
    r_abs_mag = chunk['sdssr']-5.0*np.log10(catsim_dist/10.0)

    print 'chunk size ',len(chunk),r_abs_mag.min(),r_abs_mag.max()


    with open(out_name, 'a') as out_file:
        for i_obj, (name, kep, umag, gmag, rmag, imag, zmag, dist, rr) in \
        enumerate(zip(chunk['sedfilename'], catsim_abs_mag,
                      chunk['sdssu'], chunk['sdssg'], chunk['sdssr'],
                      chunk['sdssi'], chunk['sdssz'], catsim_dist, r_abs_mag)):
            tt, feh, logg = get_physical_characteristics(name)
            out_file.write('%s %e %e %e %e %e %e %e %e %e\n' %
            (name, kep,
            umag, gmag, rmag, imag, zmag,
            dist,tt,rr))

