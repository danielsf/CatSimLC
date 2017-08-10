from __future__ import with_statement

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np

to_plot = ['kplr008240861_lc.txt', 'kplr007740566_lc.txt',
           'kplr002719436_lc.txt', 'kplr003549994_lc.txt',
            'kplr010447902_lc.txt', 'kplr009202969_lc.txt',
            'kplr011710139_lc.txt', 'kplr011401845_lc.txt',
            'kplr009388303_lc.txt', 'kplr009470175_lc.txt',
             'kplr007119876_lc.txt', 'kplr009345163_lc.txt']

def fold_light_curve(lc_name, period):
    dtype = np.dtype([('t', float), ('f', float), ('s', float)])

    data = np.genfromtxt(lc_name, dtype=dtype)
    t_folded = data['t'] % period
    dtavg = 0.01*(t_folded.max()-t_folded.min())
    t_avg_raw = np.arange(t_folded.min(),t_folded.max(),dtavg)

    f_avg = []
    s_avg = []
    t_avg = []
    for tt in t_avg_raw:
        valid_dex = np.where(np.abs(t_folded-tt)<0.5*dtavg)
        if len(valid_dex[0])==0:
            continue
        wgts = 1.0/np.power(data['s'][valid_dex],2)
        wgt_sum = wgts.sum()
        wgts = wgts/wgt_sum
        mean = data['f'][valid_dex]*wgts
        mean = mean.sum()
        f_avg.append(mean)
        var = np.power(wgts*data['s'][valid_dex],2)
        var = var.sum()
        s_avg.append(np.sqrt(var))
        t_avg.append(tt)
    f_avg = np.array(f_avg)
    s_avg= np.array(s_avg)
    t_avg = np.array(t_avg)
    n_periods = (data['t'].max()-data['t'].min())/period

    return t_avg, f_avg, s_avg, data['t'].min(),n_periods

rng = np.random.RandomState(65234)
work_dir = os.path.join('workspace', 'validation_170705')
fig_dir = os.path.join(work_dir, 'figs')
stitch_dir = os.path.join(work_dir, 'stitched')

param_file_name = os.path.join(work_dir, 'lc_params_master.txt')

# read in all models

chisq_arr = []
nt_arr = []
nc_arr = []
model_dict = {}
name_arr = []
with open(param_file_name, 'r') as in_file:
    for line in in_file:
        if line[0] == '#':
            continue
        v_line = line.strip().split()
        nt = int(v_line[1])
        nc = int(v_line[3])
        chisq = float(v_line[3+nc])
        nt_arr.append(nt)
        chisq_arr.append(chisq)
        nc_arr.append(nc)
        name_arr.append(v_line[0])
        model_dict[v_line[0]] = v_line

chisq_arr = np.array(chisq_arr)
nt_arr = np.array(nt_arr)
nc_arr = np.array(nc_arr)
name_arr = np.array(name_arr)

# limit to models with more than 0 Fourier components

nc_greater_than_0 = np.where(nc_arr>0)
chisq_arr = chisq_arr[nc_greater_than_0]
nc_arr = nc_arr[nc_greater_than_0]
nt_arr = nt_arr[nc_greater_than_0]
name_arr = name_arr[nc_greater_than_0]

power_arr= []
omega_arr = []
aa_arr = []
bb_arr = []
cc_arr = []
tau_arr = []
median_arr = []
full_models = {}
for name, nc in zip(name_arr, nc_arr):
    model = model_dict[name]
    local_power = []
    local_aa = []
    local_bb = []
    local_cc = []
    local_omega = []
    local_tau = []
    median = float(model[4+nc])
    median_arr.append(median)
    base_dex = 5+nc
    for ix in range(nc):
        aa = float(model[base_dex+5*ix])
        bb = float(model[base_dex+5*ix+1])
        cc = float(model[base_dex+5*ix+2])
        omega = float(model[base_dex+5*ix+3])
        tau = float(model[base_dex+5*ix+4])
        local_power.append(np.sqrt(aa*aa+bb*bb)/median)
        local_aa.append(aa)
        local_bb.append(bb)
        local_cc.append(cc)
        local_omega.append(omega)
        local_tau.append(tau)
    full_models[name] = {}
    full_models[name]['aa'] = local_aa
    full_models[name]['bb'] = local_bb
    full_models[name]['cc'] = local_cc
    full_models[name]['tau'] = local_tau
    full_models[name]['omega'] = local_omega
    full_models[name]['median'] = median
    local_power=np.array(local_power)
    max_dex = np.argmax(local_power)
    power_arr.append(local_power[max_dex])
    aa_arr.append(local_aa[max_dex])
    bb_arr.append(local_bb[max_dex])
    cc_arr.append(local_cc[max_dex])
    tau_arr.append(local_tau[max_dex])
    omega_arr.append(local_omega[max_dex])

#low_power_cases = np.where(np.log10(power_arr)<=-5.0)
#models_to_plot = rng.choice(low_power_cases[0], size=12, replace=False)

fig_dir = os.path.join('workspace', 'figs_170810')

dtype = np.dtype([('t', float), ('f', float), ('s', float)])

for lc_name in to_plot:
    model = full_models[lc_name]
    stitch_name = lc_name.replace('lc.txt', 'lc_stitched.txt')
    stitch_name = os.path.join(stitch_dir, stitch_name)
    stitch_data = np.genfromtxt(stitch_name, dtype=dtype)

    plt.figsize = (30,30)

    t_step = 0.1*np.diff(stitch_data['t']).min()
    print lc_name,t_step

    t0 = stitch_data['t'].min()
    t_model = np.arange(t0, stitch_data['t'].max(), t_step)
    f_model = np.ones(len(t_model))*model['median']
    for ix in range(len(model['aa'])):
        arg = model['omega'][ix]*(t_model-t0-model['tau'][ix])
        f_model += model['aa'][ix]*np.cos(arg)
        f_model += model['bb'][ix]*np.sin(arg)
        f_model += model['cc'][ix]

    fig_name = lc_name.replace("lc.txt", "fig.png")
    dt = (stitch_data['t'].max()-t0)/3.0
    dt = 5.0
    print lc_name,stitch_data['t'].min(),stitch_data['t'].max()
    t_min = t0
    for i_fig in range(3):
        plt.subplot(3,1,i_fig+1)
        if i_fig>0:
            t_min += 0.25*(stitch_data['t'].max()-stitch_data['t'].min())

        valid = []
        t_max = t_min + dt
        valid = np.where(np.logical_and(stitch_data['t']>=t_min,
                                        stitch_data['t']<=t_max))

        while len(valid[0])<50 and t_max<stitch_data['t'].max():
            t_min += dt
            t_max = t_min+dt
            valid = np.where(np.logical_and(stitch_data['t']>=t_min,
                                            stitch_data['t']<=t_max))

        plt.scatter(stitch_data['t'], stitch_data['f'], color='b', zorder=1)
        plt.plot(t_model, f_model, color='r', zorder=2)
        plt.xlim(t_min, t_max)
        valid = np.where(np.logical_and(stitch_data['t']<=t_max,
                                        stitch_data['t']>=t_min))
        if len(valid[0])==0:
            continue
        f_window = stitch_data['f'][valid]
        f_window = np.sort(f_window)
        plt.ylim(f_model.min(),f_model.max())
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, fig_name))
    plt.close()
