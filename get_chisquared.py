from __future__ import with_statement
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--lcdir', type=str, default=None)
parser.add_argument('--params', type=str, default=None)

args = parser.parse_args()

param_dict = {}
with open(args.params, 'r') as input_file:
    for line in input_file:
        if line[0] == '#':
            continue
        line = line.strip().split()
        name = line[0]
        param_dict[name] = []
        n_components = int(line[1])
        for ix in range(min(3,n_components)):
            base_dex = 4+5*ix
            aa = float(line[base_dex])
            bb = float(line[base_dex+1])
            cc = float(line[base_dex+2])
            tau = float(line[base_dex+3])
            omega = float(line[base_dex+4])
            param_dict[name].append((aa, bb, cc, tau, omega))

dtype = np.dtype([('t', float), ('f', float), ('s', float)])
chisquared_dict = {}
ct = 0
import time
t_start = time.time()
with open('chisquared_of_models.txt', 'w') as out_file:
    for file_name in param_dict:
        out_file.write('%s %d ' % (file_name, len(param_dict[file_name])))

        full_name = os.path.join(args.lcdir, file_name)
        data = np.genfromtxt(full_name, dtype=dtype)
        model = np.zeros(len(data['t']))
        for ix in range(min(len(param_dict[name]), 5)):
            params = param_dict[name][ix]
            model += params[2]
            arg = params[4]*(data['t']-data['t'].min()-params[3])
            model += params[0]*np.cos(arg)
            model += params[1]*np.sin(arg)
            chisq = np.power((data['f']-model)/data['s'],2).sum()
            out_file.write('%e ' % chisq)
        if ix<4:
            while ix<5:
                out_file.write('%e ' % chisq)
                ix += 1
        ct += 1
        if ct%1000 == 0:
            print '%d took %e' % (ct, (time.time()-t_start)/60.0)
        out_file.write('\n')

