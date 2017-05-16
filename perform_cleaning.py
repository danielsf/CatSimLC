from __future__ import with_statement

from PressRybicki import get_clean_spectrum_PressRybicki

import argparse
import os
import numpy as np
import time

parser = argparse.ArgumentParser()

parser.add_argument('--lc_list', type=str, default=None)
parser.add_argument('--out_file', type=str, default=None)
parser.add_argument('--lc_dir', type=str, default=None)
parser.add_argument('--components', type=int, default=5)
parser.add_argument('--flush_every', type=int, default=1000)

args = parser.parse_args()
if args.lc_list is None:
    raise RuntimeError('must specify list of light curve files')
if args.out_file is None:
    raise RuntimeError('must specify out_file')

lc_file_list = []
with open(args.lc_list, 'r') as input_file:
    input_lines = input_file.readlines()
    for line in input_lines:
        lc_file_list.append(line.strip())


dtype = np.dtype([('time', float), ('flux', float), ('sigma', float)])

aa_dict = {}
bb_dict = {}
cc_dict = {}
omega_dict = {}
tau_dict = {}
omega_dict = {}
span_dict = {}
nt_dict = {}


with open(args.out_file, 'w') as output_file:
    output_file.write('# lc_filename, number of components, timespan, n time steps\n')
    output_file.write('# A, B, C, tau, omega (f = A*cos(omega*(t-tau)) + B*sin(omega*(t-tau)) + C)\n')

t_start = time.time()
for file_name in lc_file_list:
    full_name = os.path.join(args.lc_dir, file_name)
    data = np.genfromtxt(full_name, dtype=dtype)
    sorted_dex = np.argsort(data['time'])
    dt = np.diff(data['time'][sorted_dex])
    dt_min = dt.min()
    (aa, bb, cc,
     omega, tau,
     freq) = get_clean_spectrum_PressRybicki(data['time'],
                                             data['flux'],
                                             data['sigma'],
                                             0.1*dt_min,
                                             max_components=args.components)

    aa_dict[file_name] = aa
    bb_dict[file_name] = bb
    cc_dict[file_name] = cc
    tau_dict[file_name] = tau
    omega_dict[file_name] = omega
    span_dict[file_name] = data['time'].max()-data['time'].min()
    nt_dict[file_name] = len(data['time'])

    #print('done with %d after %e' % (len(aa_dict), time.time()-t_start))

    if len(aa_dict) >= args.flush_every or file_name == lc_file_list[-1]:
        with open(args.out_file, 'a') as output_file:
            for file_name in aa_dict.keys():

                aa = aa_dict[file_name]
                bb = bb_dict[file_name]
                cc = cc_dict[file_name]
                tau = tau_dict[file_name]
                omega = omega_dict[file_name]

                output_file.write('%s %d %.3e %d ' %
                                  (file_name, len(aa),
                                   span_dict[file_name],
                                   nt_dict[file_name]))

                for ix in range(len(aa)):
                    output_file.write('%.6e %.6e %.6e %.6e %.6e '
                                      % (aa[ix], bb[ix], cc[ix],
                                         tau[ix], omega[ix]))
                output_file.write('\n')

        aa_dict = {}
        bb_dict = {}
        cc_dict = {}
        tau_dict = {}
        omega_dict = {}
        span_dict = {}
        nt_dict = {}

print('generating clean light curves took %e' % (time.time()-t_start))
