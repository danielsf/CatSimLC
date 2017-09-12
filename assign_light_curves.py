import numpy as np
import os
from scipy.spatial import KDTree
import gc
import sys
import multiprocessing as mproc

import argparse

class LightCurveFitter(object):

    def __init__(self, kep_data, rms_data):
        """
        kep_data = a numpy recarray with the valid data from the Kepler catalog
        rms_data = a numpy recarray with data from the light curve catalog
        """

        if len(kep_data) != len(rms_data):
            raise RuntimeError("len of kep_data does not match len of rms_data")

        # absolute r-band magnitude
        kep_abs_r = new_mags['r'] - 5.0*np.log10(new_mags['dist']/10.0)

        # build the KD tree we will use to associate CatSim
        # stars to Kepler light curves
        kep_color = kep_data['g'] - kep_data['r']
        kep_params = np.array([kep_abs_r, kep_color]).transpose()
        print('creating KDTree')
        self.kep_kd_tree = KDTree(kep_params, leafsize=1)
        self.kep_data = kep_data

    def fit_lc(self, catsim_data=None, out_name=None, lock=None, t_off=None):

        catsim_color = catsim_data['g'] - catsim_data['r']
        # CatSim parallaxes are in milliarcseconds
        catsim_distance = _au_to_parsec/radiansFromArcsec(catsim_data['parallax']*0.001)
        catsim_abs_r = catsim_data['r'] - 5.0*np.log10(catsim_distance/10.0)

        catsim_params = np.array([catsim_abs_r, catsim_color]).transpose()

        match_dist, out_dexes = self.kep_kd_tree.query(catsim_params, k=1)

        ids_in_use = self.kep_data['id'][out_dexes]

        lock.acquire()
        with open(out_name, 'a') as out_file:
            for i_star in range(len(ids_in_use)):
                paramStr = '{"m":"kplr", "p":{"lc":%d, "t0":%.3f}}' % (ids_in_use[i_star],
                                                                       t_off[i_star])

                out_file.write('%d;%d;%s\n' %
                               (catsim_data['id'][i_star],
                                catsim_data['htmid'][i_star],
                                paramStr))

        lock.release()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--kepler_catalog', type=str,
                        default='KIC/kic_data.txt',
                        help='file containing the input catalog of kepler sources')

    parser.add_argument('--lc_file', type=str,
                        default=None,
                        help='The file containing light curve chisquared/dof and '
                        'time span data (produced by get_rms_lookup.py')

    parser.add_argument('--output', type=str,
                        default=None,
                        help='file where we will write the final light curve '
                        'assignments')

    parser.add_argument('--seed', type=int,
                        default=None,
                        help='seed for the random number generator that will '
                        'be used to create the phase offsets for light curves')

    parser.add_argument('--chisq_cutoff', type=float,
                        default=700.0,
                        help='The maximum chisquared/dof value that will be '
                        'allowed for valid light curves')

    parser.add_argument('--span_cutoff', type=float,
                        default=300.0,
                        help='The minimum time span taht will be allowed '
                        'for valid light curves')

    parser.add_argument('--catsim_table', type=str,
                        default=None,
                        help='The CatSim table to fit')

    parser.add_argument('--n_proc', type=int,
                        default=1,
                        help='Number of processes to spawn')

    args = parser.parse_args()
    if args.seed is None:
        raise RuntimeError('must specify a seed')

    if args.output is None:
        raise RuntimeError('must specify a file for output')

    if args.lc_file is None:
        raise RuntimeError('must specify a file to read '
                           'light curve parameters from')

    if args.catsim_table is None:
        raise RuntimeError('must specify a CatSim table to fit '
                           'to the Kepler light curves')

    sys.setrecursionlimit(100000)
    rng = np.random.RandomState(args.seed)

    from lsst.sims.catalogs.db import DBObject

    db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                  port=1433, driver='mssql+pymssql')

    # only accept Kepler stars for which there is a valid distance
    kepler_dtype = np.dtype([('id', int), ('ra', float), ('dec', float), ('u', float),
                             ('g', float), ('r', float), ('i', float), ('z', float),
                             ('dist', float), ('teff_kic', float), ('teff_stellar', float),
                             ('Av', float), ('ebv', float)])
    kep_data = np.genfromtxt(args.kepler_catalog, dtype=kepler_dtype)
    valid = np.where(kep_data['dist']>0.0)
    kep_data = kep_data[valid]

    # only accept light curves that have non-negative flux,
    # valid chisquared/dof and valid time spans of data
    rms_dtype = np.dtype([('name', str, 100), ('rms', float), ('chisq', float),
                          ('span', float), ('min_flux', float)])
    rms_data = np.genfromtxt(args.lc_file, dtype=rms_dtype)
    valid = np.where(np.logical_and(rms_data['min_flux']>0.0,
                     np.logical_and(rms_data['chisq']<args.chisq_cutoff,
                                    rms_data['span']>args.span_cutoff)))

    rms_data = rms_data[valid]

    # Throw out light curves that do not appear in both the catalog and
    # the list of valid light curves

    rms_dict = dict([(int(nn.split('_')[0][4:]),ix)
                     for ix, nn in enumerate(rms_data['name'])])

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

    # correct the Kepler catalog photometry according to
    # eqns 1 and 2 of
    # Pinsonneault et al 2012
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

    del kep_data
    del rms_dict
    gc.collect()

    fitter_list = []
    for i_process in range(args.n_proc):
        fitter_list.append(LightCurveFitter(new_mags, rms_data))

    if sys.version_info.major == 2:
        id_type = long
    else:
        id_type = int

    catsim_dtype = np.dtype([('id', id_type), ('htmid', id_type), ('g', float),
                             ('r', float), ('parallax', float)])

    from lsst.sims.utils import radiansFromArcsec
    _au_to_parsec = 1.0/206265.0

    query = 'SELECT simobjid, htmid, sdssg, sdssr, parallax '
    query += 'FROM %s ' % args.catsim_table

    chunk_iterator = db.get_arbitrary_chunk_iterator(query, chunk_size=1000000,
                                                     dtype=catsim_dtype)

    with open(args.output, 'w') as out_file:
        out_file.write('# simobjid htmid varParamStr\n')

    print("starting search")
    import time
    total = 0
    t_start = time.time()

    for chunk in chunk_iterator:

        process_list = []
        lock = mproc.Lock()
        t_offset_chunk = rng.random_sample(len(chunk))*3652.5
        for i_process in range(args.n_proc):
            if args.n_proc==1:
                catsim_data = chunk
            else:
                if i_process==args.n_proc-1:
                    i_start = (args.n_proc-1)*len(chunk)//args.n_proc
                    catsim_data = chunk[i_start:]
                    t_offset = t_offset_chunk[i_start:]
                else:
                    i_start = i_process*len(chunk)//args.n_proc
                    i_end = (i_process+1)*len(chunk)//args.n_proc
                    catsim_data = chunk[i_start:i_end]
                    t_offset = t_offset_chunk[i_start:i_end]

            p = mproc.Process(target=fitter_list[i_process].fit_lc,
                              kwargs={'catsim_data': catsim_data,
                                      'lock': lock,
                                      't_off': t_offset,
                                      'out_name': args.output})
            p.start()
            process_list.append(p)

        for p in process_list:
            p.join()

        total += len(chunk)
        elapsed = time.time()-t_start
        project_billion = 1.0e9*elapsed/total
        elapsed = elapsed/3600.0
        project_billion = project_billion/3600.0
        print('did %d in %e hours; could do a billion in %e hours' %
              (total, elapsed, project_billion))

    print('max_len of paramStr %d' % max_len_paramStr)
