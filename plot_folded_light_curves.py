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

print 'max chisq_arr_doff ',chisq_arr_dof.max()

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

"""
models_to_plot = []
monday_dir = os.path.join(work_dir, 'monday_stitched')
stitch_dir = monday_dir
monday_files = os.listdir(monday_dir)
for file_name in monday_files:
    key_name = file_name.replace('_stitched','')
    for ix in range(len(name_arr)):
        if name_arr[ix] == key_name:
            models_to_plot.append(ix)
models_to_plot = np.array(models_to_plot)
"""

sorted_dex = np.argsort(chisq_arr_dof)
models_to_plot = []
models_to_plot.append(sorted_dex[-1])
models_to_plot.append(sorted_dex[-2])
models_to_plot.append(sorted_dex[-3])
models_to_plot.append(sorted_dex[-4])

ii = len(sorted_dex)/2
models_to_plot.append(sorted_dex[ii])
models_to_plot.append(sorted_dex[ii+1])
models_to_plot.append(sorted_dex[ii+2])
models_to_plot.append(sorted_dex[ii+3])

ii=len(sorted_dex)/10
models_to_plot.append(sorted_dex[ii])
models_to_plot.append(sorted_dex[ii+1])
models_to_plot.append(sorted_dex[ii+2])
models_to_plot.append(sorted_dex[ii+3])

for ii in models_to_plot:
    print 'going to plot ',chisq_arr_dof[ii]

not_there = False
for name in name_arr[models_to_plot]:
    stitch_name = os.path.join(stitch_dir,name.replace('.txt','_stitched.txt'))
    if not os.path.exists(stitch_name):
        print 'need to get %s' % stitch_name
        not_there = True
if not_there:
    exit()

i_fig = 0
i_plot=0

fold_dict = {}

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
    plt.subplot(2,2,i_fig)
    plt.errorbar(t_fold,f_fold,yerr=s_fold,zorder=1,color='b', fmt='o')
    plt.plot(t_model,f_model,color='r', zorder=2)
    fold_dict[ix] = {}
    fold_dict[ix]['t_fold'] = t_fold
    fold_dict[ix]['f_fold'] = f_fold
    fold_dict[ix]['s_fold'] = s_fold
    fold_dict[ix]['t_model'] = t_model
    fold_dict[ix]['f_model'] = f_model
    #plt.plot(t_model,f_full_model,color='g',zorder=2)
    title = '%s\n$\chi^2$/dof=%.2e;\npower=%.2e\nN periods %.2e\nN components %d' % \
    (name, chisq_arr_dof[ix], power_arr[ix],n_periods,len(full_models[name]['aa']))

    print 'plotting ',chisq_arr_dof[ix]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    ylim0 = plt.ylim()
    plt.text(t_model.max()-0.5*(t_model.max()-t_model.min()),
             ylim0[1]+0.05*(ylim0[1]-ylim0[0]),
             title,fontsize=7)

    plt.ylim((ylim0[0],ylim0[1]+0.5*(ylim0[1]-ylim0[0])))

    if i_fig==1:
        plt.xlabel('MJD')
        plt.ylabel('flux')

    if i_fig==4 or ix == models_to_plot[-1]:
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir,'monday_talk_examples_again_%d.eps' % i_plot))
        plt.close()
        plt.figsize=(30,30)
        i_fig=0
        i_plot+=1


dtype = np.dtype([('t', float), ('f', float), ('s', float)])
for i_place, ix in enumerate(models_to_plot):
    name = name_arr[ix]
    period = 2.0*np.pi/omega_arr[ix]
    stitch_name = os.path.join(stitch_dir, name.replace('.txt', '_stitched.txt'))
    data = np.genfromtxt(stitch_name, dtype=dtype)
    #t_model = np.arange(data['t'].min(), data['t'].min()+0.05*(data['t'].max()-data['t'].min()), 0.01)
    t_range = max(5.0*period, 50.0)
    t_model = np.arange(210.0, 210.0+t_range, 0.01)
    median = median_arr[ix]
    t0 = data['t'].min()
    f_full_model = np.ones(len(t_model))*median
    for ic in range(len(full_models[name]['aa'])):
        sin_arg = full_models[name]['omega'][ic]*(t_model-t0-full_models[name]['tau'][ic])
        f_full_model += full_models[name]['aa'][ic]*np.cos(sin_arg)
        f_full_model += full_models[name]['bb'][ic]*np.sin(sin_arg)
        f_full_model += full_models[name]['cc'][ic]

    f_sorted = np.sort(data['f'])
    f_min = f_sorted[len(f_sorted)/20]
    f_max = f_sorted[19*len(f_sorted)/20]
    f_valid_dex = np.where(np.logical_and(data['t']>t_model.min(), data['t']<t_model.min()+2.0*t_range/3.0))
    f_valid = data['f'][f_valid_dex]
    f_valid_sorted = np.sort(f_valid)
    f_min = f_valid_sorted[5]
    f_max = f_valid_sorted[-6]

    plt.figsize = (30,30)

    plt.subplot(3,1,3)
    plt.errorbar(fold_dict[ix]['t_fold'],fold_dict[ix]['f_fold'],yerr=fold_dict[ix]['s_fold'],
                 zorder=1,color='b', fmt='o')
    plt.plot(fold_dict[ix]['t_model'],fold_dict[ix]['f_model'],color='r', zorder=2)
    plt.xlabel('MJD')
    plt.ylabel('flux')
    plt.title('folded light curve')

    for i_fig in range(2):
        t_min = t_model.min() + i_fig*t_range/3.0
        if i_fig>0:
            t_min -= t_range/15.0
        t_max = t_model.min() + (i_fig+1)*t_range/4.0
        t_max += t_range/15.0

        plt.subplot(3,1,i_fig+1)
        plt.errorbar(data['t'],data['f'],data['s'],zorder=1,color='b',linestyle='')
        plt.plot(t_model,f_full_model,color='r',zorder=2)
        plt.xlim(t_min,t_max)
        plt.ylim(f_min,f_max)
        if i_fig==0:
            plt.xlabel('MJD')
            plt.ylabel('flux')
            plt.title('$\chi^2$ per dof = %.2e' % chisq_arr_dof[ix])
    plt.tight_layout()
    fig_name = os.path.join(fig_dir,name.replace('.txt','.png'))
    if i_place>7:
        fig_name = fig_name.replace('figs/','figs/good_')
    plt.savefig(fig_name)
    plt.close()
