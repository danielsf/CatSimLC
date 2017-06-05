from __future__ import with_statement
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--obj_list', type=str, default=None)
parser.add_argument('--out_dir', type=str, default=None)
parser.add_argument('--write_every', type=int, default=104)
parser.add_argument('--in_dir_root', type=str, default='lc_nostitch')
parser.add_argument('--n_in_dir', type=int, default=6)

args = parser.parse_args()
if args.out_dir is None:
    raise RuntimeError('must specify out_dir')

if args.obj_list is None:
    raise RuntimeERror('must specify obj_list')

dtype = np.dtype([('lc', str, 300)])
obj_list = np.genfromtxt(args.obj_list, dtype=dtype)

out_dir = args.out_dir
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

dtype = np.dtype([('t', float), ('f', float), ('s', float)])

dir_list = []
for ix in range(1,args.n_in_dir+1):
    dir_list.append('%s_%d' % (args.in_dir_root, ix))



files_in_dir = {}
for in_dir in dir_list:
    files_in_dir[in_dir] = os.listdir(in_dir)

import time
t_start = time.time()

data_dict = {}

for file_name in obj_list['lc']:
    print 'trying to create ',file_name
    flux = None
    mjd = None
    sig = None
    last_t = -1.0
    for in_dir in dir_list:
        if not file_name in files_in_dir[in_dir]:
            continue

        in_name = os.path.join(in_dir, file_name)
        data = np.genfromtxt(in_name, dtype=dtype)
        if data['t'].min() < last_t:
            continue

        last_t = data['t'].max()

        if flux is None:
            flux = data['f']
            mjd = data['t']
            sig = data['s']
        else:
            flux = np.append(flux, data['f'] - data['f'][0] + flux[-1])
            mjd = np.append(mjd, data['t'])
            sig = np.append(sig, data['s'])

    data_dict[file_name] = (mjd, flux, sig)

    if (len(data_dict) >= args.write_every or
        file_name == obj_list['lc'][-1]):

        for file_name in data_dict.keys():
            out_name = os.path.join(out_dir, file_name)
            with open(out_name, 'w') as out_file:
               obj_data = data_dict[file_name]
               for (tt, ff, ss) in zip(obj_data[0], obj_data[1], obj_data[2]):
                   out_file.write('%.12e %e %e\n' % (tt, ff, ss))

        del data_dict
        data_dict = {}
