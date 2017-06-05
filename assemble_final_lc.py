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
header_dict = {}

for file_name in obj_list['lc']:
    try:
        flux = None
        mjd = None
        sig = None
        last_t = -1.0
        for in_dir in dir_list:
            if not file_name in files_in_dir[in_dir]:
                continue

            in_name = os.path.join(in_dir, file_name)
            data = np.genfromtxt(in_name, dtype=dtype)

            if flux is None:
                flux = data['f']
                mjd = data['t']
                sig = data['s']
            else:
                flux = np.append(flux, data['f'])
                mjd = np.append(mjd, data['t'])
                sig = np.append(sig, data['s'])

            # get the "header" information encoding
            # where segments begin and end
            with open(in_name, 'r') as in_file:
                in_lines = in_file.readlines()
            in_lines = np.array(in_lines)
            comment_dexes = np.where(np.char.find(in_lines, '#')==0)
            if file_name not in header_dict:
                header_dict[file_name] = []
            for ix in comment_dexes[0][1:]:
                end_pts = in_lines[ix].strip().split()
                header_dict[file_name].append((float(end_pts[1]), float(end_pts[2])))

        sorted_dexes = np.argsort(mjd)
        mjd = mjd[sorted_dexes]
        flux = flux[sorted_dexes]
        sig = sig[sorted_dexes]

        data_dict[file_name] = (mjd, flux, sig)

        if (len(data_dict) >= args.write_every or
            file_name == obj_list['lc'][-1]):

            for file_name in data_dict.keys():
                out_name = os.path.join(out_dir, file_name)
                with open(out_name, 'w') as out_file:
                    for end_pts in header_dict[file_name]:
                        out_file.write('# %.12e %.12e\n' % (end_pts[0], end_pts[1]))
                    obj_data = data_dict[file_name]
                    for (tt, ff, ss) in zip(obj_data[0], obj_data[1], obj_data[2]):
                        out_file.write('%.12e %e %e\n' % (tt, ff, ss))

            del data_dict
            data_dict = {}
            del header_dict
            header_dict = {}
    except:
        print('failed on %s' % file_name)
        raise
