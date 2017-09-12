from __future__ import with_statement
from __future__ import print_function
import os
import numpy as np

import argparse
import time

parser = argparse.ArgumentParser()

parser.add_argument('--cutoff', type=float, default=700.0,
                    help='chisq/doff cutoff')

parser.add_argument('--out', type=str, default='rms_variability_lookup.txt',
                    help='file associating LC name with rms variability')

parser.add_argument('--workdir', type=str,
                    default=os.path.join('workspace/validation_170824'),
                    help='root dir for workspace')

parser.add_argument('--params', type=str,
                    default='master_lc_params_170824.txt',
                    help='file containing parameters')


args = parser.parse_args()

param_name = os.path.join(args.workdir, args.params)

full_model_dict = {}
name_arr = []
chisq_arr = []
t_span_arr = []

with open(param_name, 'r') as input_file:
    for line in input_file:
        if line[0] == '#':
            continue
        params = line.strip().split()
        name = params[0]
        t_span = float(params[2])
        n_c = int(params[3])
        median = float(params[4+n_c])
        if median < 0.0:
            continue
        if n_c < 2:
            continue
        dof = int(params[1])
        chisq = float(params[3+n_c])

        local_aa = []
        local_bb = []
        local_cc = []
        local_omega = []
        local_tau = []

        for i_c in range(n_c):
            base_dex = 5+n_c+i_c*5
            local_aa.append(float(params[base_dex]))
            local_bb.append(float(params[base_dex+1]))
            local_cc.append(float(params[base_dex+2]))
            local_omega.append(float(params[base_dex+3]))
            local_tau.append(float(params[base_dex+4]))
        local_aa = np.array(local_aa)
        local_bb = np.array(local_bb)
        local_cc = np.array(local_cc)
        local_tau = np.array(local_tau)
        local_omega = np.array(local_omega)
        local_power = np.power(local_aa,2)+np.power(local_bb,2)
        max_dex = np.argmax(local_power)

        name_arr.append(name)
        t_span_arr.append(t_span)
        flux_power = np.sqrt(local_power[max_dex])
        mag_power = 2.5*np.log10(1.0 + flux_power/median)
        chisq_arr.append(chisq/dof)

        full_model_dict[name] = {}
        full_model_dict[name]['a'] = local_aa
        full_model_dict[name]['b'] = local_bb
        full_model_dict[name]['c'] = local_cc
        full_model_dict[name]['median'] = median
        full_model_dict[name]['tau'] = local_tau
        full_model_dict[name]['omega'] = local_omega


print('valid models %d' % len(t_span_arr))

t_model = np.arange(0.0, 5.0*365.25, 30.0/(24.0*60.0))
f_model = np.ones(len(t_model), dtype=float)

t_start = time.time()
ct = 0
with open(args.out, 'w') as out_file:
    for i_model, name in enumerate(name_arr):
        if i_model>0 and i_model%10==0:
            print("%d took %.2e hours; expect total %.2e -- %d" %
                  (i_model,(time.time()-t_start)/3600.0,
                   len(name_arr)*(time.time()-t_start)/(3600.0*i_model),
                   ct))
        full_model = full_model_dict[name]
        quiescent_flux = full_model['median']
        for cc in full_model['c']:
            quiescent_flux += cc

        f_model *= 0.0
        f_model += quiescent_flux
        for aa, bb, omega, tau in \
        zip(full_model['a'], full_model['b'], full_model['omega'], full_model['tau']):
            arg = omega*(t_model-tau)
            f_model += aa*np.cos(arg)
            f_model += bb*np.sin(arg)

        min_flux = f_model.min()
        quiescent_flux = np.median(f_model)
        if quiescent_flux<0.0:
            continue

        quiescent_mag = -2.5*np.log10(quiescent_flux)

        if f_model.min()<0.0:
            invalid = np.where(f_model<0.0)
            f_model[invalid] = 1.0e-20
        else:
            ct += 1

        m_model = -2.5*np.log10(f_model)

        rms_var = np.sqrt(np.power(m_model-quiescent_mag,2).sum()/(len(m_model)-1.0))

        out_file.write('%s %e %e %e %e\n' % (name, rms_var, chisq_arr[i_model],
                                             t_span_arr[i_model],min_flux))

print("reasonable models %d" % ct)
