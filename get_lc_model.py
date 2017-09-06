from __future__ import with_statement
import os
import numpy as np

import subprocess

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--id', type=str, nargs='+', default = None)
args = parser.parse_args()

if not isinstance(args.id, list):
    id_list = [args.id]
else:
    id_list = args.id

work_dir = os.path.join('workspace', 'validation_170824')
stitch_dir = os.path.join(work_dir, 'stitched')
full_dir = os.path.join(work_dir, 'full_models')

lc_name = os.path.join(work_dir, 'master_lc_params_170824.txt')

name_arr = []
period_arr = []
power_arr = []
chisq_arr = []

full_model_dict = {}

with open(lc_name, 'r') as input_file:
    for line in input_file:
        if line[0] == '#':
            continue
        params = line.strip().split()
        name = params[0]
        n_c = int(params[3])
        median = float(params[4+n_c])
        if median < 0.0:
            continue
        if n_c < 2:
            continue
        dof = int(params[1])
        chisq = float(params[3+n_c])
        if name == 'kplr009472174_lc.txt':
            print 'target chisq dof ',chisq/dof

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
        flux_power = np.sqrt(local_power[max_dex])
        mag_power = 2.5*np.log10(1.0 + flux_power/median)
        power_arr.append(mag_power)
        period_arr.append(2.0*np.pi/local_omega[max_dex])
        chisq_arr.append(chisq/dof)

        full_model_dict[name] = {}
        full_model_dict[name]['a'] = local_aa
        full_model_dict[name]['b'] = local_bb
        full_model_dict[name]['c'] = local_cc
        full_model_dict[name]['median'] = median
        full_model_dict[name]['tau'] = local_tau
        full_model_dict[name]['omega'] = local_omega
        full_model_dict[name]['chisq/dof'] = chisq/dof

print 'valid models %d' % len(name_arr)

stitch_dtype = np.dtype([('t', float), ('f', float), ('s', float)])

for lc_id in id_list:
    stitch_name = os.path.join(stitch_dir, 'kplr%s_lc_stitched.txt' % lc_id)
    if not os.path.exists(stitch_name):
        subprocess.call(["bash", "get_lc.sh", "%s" % lc_id])

    stitch_data = np.genfromtxt(stitch_name, dtype=stitch_dtype)
    name = 'kplr%s_lc.txt' % lc_id
    full_model = full_model_dict[name]
    d_t = np.diff(stitch_data['t']).min()*0.1
    t_model = np.arange(stitch_data['t'].min(),stitch_data['t'].max(),d_t)
    f_model =np.ones(len(t_model))*full_model['median']
    f_data = np.ones(len(stitch_data['t']))*full_model['median']
    for i_c in range(len(full_model['a'])):
        f_model += full_model['c'][i_c]
        arg = full_model['omega'][i_c]*(t_model-stitch_data['t'].min()-full_model['tau'][i_c])
        f_model += full_model['a'][i_c]*np.cos(arg)
        f_model += full_model['b'][i_c]*np.sin(arg)

        f_data += full_model['c'][i_c]
        arg = full_model['omega'][i_c]*(stitch_data['t']-stitch_data['t'].min()-full_model['tau'][i_c])
        f_data += full_model['a'][i_c]*np.cos(arg)
        f_data += full_model['b'][i_c]*np.sin(arg)

    print name,full_model['chisq/dof']

    with open(os.path.join(full_dir, 'kplr%s_lc_model.txt' % lc_id), 'w') as out_model:
        for i_t in range(len(t_model)):
            out_model.write('%e %e\n' % (t_model[i_t], f_model[i_t]))
