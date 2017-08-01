from __future__ import with_statement
import numpy as np
import os
import gzip
import time

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

kep_stellar_data = np.genfromtxt('../data/kepler_stellar17.csv',
                                 dtype=kep_dtype,
                                 delimiter='|', skip_header=1)

print 'dec range ',kep_stellar_data['degree_dec'].min(),kep_stellar_data['degree_dec'].max()

list_of_files = os.listdir('.')

to_find = np.copy(kep_stellar_data['kepid'])
stellar_dex_dict = {}
for i_line, val in enumerate(kep_stellar_data['kepid']):
    stellar_dex_dict[val] = i_line
been_found = []

id_dex = 0
pos_dex = 1
sdss_dex = 3
param_dex = 14

from lsst.sims.utils import angularSeparation

distance_max = -1.0

t_start = time.time()
with open('kic_data.txt', 'w') as out_file:
    out_file.write('# kepid sdss_g sdss_r distance(in parsecs) Teff(from KIC)\n')
    for file_name in list_of_files:
        if '.dat' in file_name:
            if file_name[0] == 'n':
                dec_sgn = 1.0
            else:
                dec_sgn = -1.0
            dec_val = float(file_name[1:3])*dec_sgn
            possible_sources = np.where(np.logical_and(kep_stellar_data['degree_dec']>=dec_val,
                                                       kep_stellar_data['degree_dec']<=dec_val+1.0))
            if len(possible_sources[0]) == 0:
                continue
            print file_name,len(possible_sources[0]),len(been_found)
            elapsed = time.time()-t_start
            elapsed = elapsed/3600.0
            n_found = len(been_found)
            print '%d in %.2e hours; should take %.2e hours' % \
            (n_found,elapsed,len(to_find)*elapsed/max(n_found,1))
            if file_name.endswith('.gz'):
                open_cmd = gzip.open
            else:
                open_cmd = open
            with open_cmd(file_name, 'r') as in_file:
                in_cat_lines = in_file.readlines()

            id_dict = {}
            for i_line, line in enumerate(in_cat_lines):
                if line[0] == '#':
                    continue
                data_v = line.split('|')
                id_dict[int(data_v[0])] = i_line

            sources_to_be_found = kep_stellar_data['kepid'][possible_sources]
            for source_id in sources_to_be_found:
                if source_id not in id_dict:
                    continue
                line_dex = id_dict[source_id]
                line = in_cat_lines[line_dex]
                data_v = line.strip().split('|')
                data_id = int(data_v[0])
                assert source_id == data_id
                if data_id in been_found:
                    raise RuntimeError('found %d again in %s' % data_id,file_name)
                if data_id in to_find:
                    been_found.append(data_id)
                    stellar_dex = stellar_dex_dict[data_id]
                    try:
                        pos_row = data_v[pos_dex]
                        sdss_row = data_v[sdss_dex]
                        if len(data_v)>param_dex:
                            param_row = data_v[param_dex]
                        else:
                            param_row = None
                    except:
                        print 'offending row'
                        print line
                        raise
                    ra = float(pos_row[:10])
                    dec = float(pos_row[10:])
                    distance = angularSeparation(ra, dec,
                                                 kep_stellar_data['degree_ra'][stellar_dex],
                                                 kep_stellar_data['degree_dec'][stellar_dex])
                    distance_arcsec = distance*3600.0
                    if distance_arcsec > distance_max:
                        distance_max = distance_arcsec
                        print '    %d distance_max %.2e in arcsec' % (data_id, distance_max)
                        print '    %e %e %e %e' % (ra, dec,
                               kep_stellar_data['degree_ra'][stellar_dex],
                               kep_stellar_data['degree_dec'][stellar_dex])

                    try:
                        sdss_g = float(sdss_row[7:14])
                        sdss_r = float(sdss_row[14:21])
                    except:
                        continue
                    if param_row is not None:
                        teff = float(param_row[:6])
                    else:
                        teff = -999.0
                    if np.isnan(kep_stellar_data['dist'][stellar_dex]):
                        dd = -999.0
                    else:
                        dd = kep_stellar_data['dist'][stellar_dex]
                    out_file.write('%d %le %le %le %le\n' %
                    (data_id, sdss_g, sdss_r, dd, teff))
