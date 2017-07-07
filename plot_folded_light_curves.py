from __future__ import with_statement

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np

def fold_light_curve(lc_name, period):
    dtype = np.dtype([('t', float), ('f', float), ('s', float)])

    data = np.genfromtxt(lc_name, dtype=dtype)
    t_folded = data['t'] % period
    dtavg = 0.01*(t_folded.max()-t_folded.min())
    t_avg = np.arange(t_folded.min(),t_folded.max(),dtavg)

    f_avg = []
    s_avg = []
    for tt in t_avg:
        valid_dex = np.where(np.abs(t_folded-tt)<0.5*dtavg)
        wgts = 1.0/np.power(data['s'][valid_dex],2)
        wgt_sum = wgts.sum()
        wgts = wgts/wgt_sum
        mean = data['f'][valid_dex]*wgts
        mean = mean.sum()
        f_avg.append(mean)
        var = np.power(wgts*data['s'][valid_dex],2)
        var = var.sum()
        s_avg.append(np.sqrt(var))
    f_avg = np.array(f_avg)
    s_avg= np.array(s_avg)
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

# take only the top 95 percent in chisq/dof

chisq_arr_dof = chisq_arr/nt_arr
chisq_arr_dof_sorted = np.sort(chisq_arr_dof)
cutoff_95 = chisq_arr_dof_sorted[95*len(chisq_arr_dof_sorted)/100]
print 'chisq_dof cutoff ',cutoff_95

percentile_95 = np.where(chisq_arr_dof<=cutoff_95)
nc_arr = nc_arr[percentile_95]
chisq_arr = chisq_arr[percentile_95]
chisq_arr_dof = chisq_arr_dof[percentile_95]
nt_arr = nt_arr[percentile_95]
name_arr = name_arr[percentile_95]

# for each light curve, find the maximum power component

power_arr= []
omega_arr = []
aa_arr = []
bb_arr = []
cc_arr = []
tau_arr = []
median_arr = []
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
    local_power=np.array(local_power)
    max_dex = np.argmax(local_power)
    power_arr.append(local_power[max_dex])
    aa_arr.append(local_aa[max_dex])
    bb_arr.append(local_bb[max_dex])
    cc_arr.append(local_cc[max_dex])
    tau_arr.append(local_tau[max_dex])
    omega_arr.append(local_omega[max_dex])

low_power_cases = np.where(np.log10(power_arr)<=-5.0)
models_to_plot = rng.choice(low_power_cases[0], size=12, replace=False)

for name in name_arr[models_to_plot]:
    stitch_name = os.path.join(stitch_dir,name.replace('.txt','_stitched.txt'))
    if not os.path.exists(stitch_name):
        print 'need to get %s' % name_arr[models_to_plot]
        exit()

i_fig = 0
i_plot=0
plt.figsize = (30,30)
for ix in models_to_plot:
    name = name_arr[ix]
    power = power_arr[ix]
    aa = aa_arr[ix]
    bb = bb_arr[ix]
    cc = cc_arr[ix]
    omega = omega_arr[ix]
    tau = tau_arr[ix]
    median = median_arr[ix]

    period = 2.0*np.pi/omega
    stitch_name = os.path.join(stitch_dir, name.replace('.txt', '_stitched.txt'))
    t_fold, f_fold, s_fold, t0, n_periods = fold_light_curve(stitch_name, period)

    t_model = np.arange(t_fold.min(), t_fold.min()+period, period*0.001)
    f_model = np.ones(len(t_model))*median
    f_model += cc
    f_model += aa*np.cos(omega*(t_model-t0-tau))
    f_model += bb*np.sin(omega*(t_model-t0-tau))

    i_fig += 1
    plt.subplot(3,2,i_fig)
    plt.errorbar(t_fold,f_fold,yerr=s_fold,zorder=1,color='b', fmt='o')
    plt.plot(t_model,f_model,color='r', zorder=2)
    title = '%s\n$\chi^2$/dof=%.2e;\npower=%.2e\nN periods %.2e' % \
    (name, chisq_arr_dof[ix], power_arr[ix],n_periods)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    ylim0 = plt.ylim()
    plt.text(t_model.max()-0.5*(t_model.max()-t_model.min()),
             ylim0[1]+0.05*(ylim0[1]-ylim0[0]),
             title,fontsize=7)

    plt.ylim((ylim0[0],ylim0[1]+0.5*(ylim0[1]-ylim0[0])))

    if i_fig==6:
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir,'low_power_examples_%d.eps' % i_plot))
        plt.close()
        i_fig=0
        i_plot+=1

