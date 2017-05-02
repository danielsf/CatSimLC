import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

def brute_sums(time_arr, f_arr, freq_arr):

    sh = np.array([(f_arr*np.sin(2.0*np.pi*nu*time_arr)).sum() for nu in freq_arr])
    ch = np.array([(f_arr*np.cos(2.0*np.pi*nu*time_arr)).sum() for nu in freq_arr])

    return ch, sh


def f_of_t(time_arr):
    y = 0.1*np.sin(0.98*time_arr)
    y += 0.35*np.cos(1.6*time_arr)
    y += 0.77*np.sin(2.6*time_arr)
    y += 0.45*np.cos(0.7*time_arr)
    return y

import time
from fft import extirp_sums

rng = np.random.RandomState(45)
time_arr = rng.random_sample(20)*56.9
time_arr = np.sort(time_arr)
fn_arr = f_of_t(time_arr)

delta = (time_arr.max()-0.0)/(4*2048.0)

t_start = time.time()
cos_fft_sum, sin_fft_sum, ttk, hhk = extirp_sums(time_arr, fn_arr, delta)
t_extirp = time.time()-t_start

#hhk_control = f_of_t(ttk)

#print('hhk residual %e ' % np.abs(hhk-hhk_control).max())
#exit()

freq_arr = np.array([k/(delta*len(cos_fft_sum)) for k in range(len(cos_fft_sum))])
t_start = time.time()
cos_brute, sin_brute = brute_sums(time_arr, fn_arr, freq_arr)
t_brute = time.time()-t_start

print 't_brute ',t_brute
print 't_extirp ',t_extirp
print 'len ',len(cos_fft_sum)


rat_arr = []
k_arr = []
nu_arr = []
for nu, cos_f, cos_b, sin_f, sin_b in zip(freq_arr, cos_fft_sum, cos_brute,
                                          sin_fft_sum, sin_brute):


    k = nu*delta*len(ttk)
    nu_arr.append(nu)
    k_arr.append(k)
    if np.abs(cos_f/cos_b-1.0) > np.abs(sin_f/sin_b-1.0):
        rat_arr.append(cos_f/cos_b)
    else:
        rat_arr.append(sin_f/sin_b)


plt.figsize = (30,30)
plt.subplot(2,1,1)
plt.plot(nu_arr,rat_arr)
plt.xlabel('$\\nu$')
plt.yscale('log')
plt.xscale('log')
plt.subplot(2,1,2)
plt.plot(k_arr,rat_arr)
plt.xlabel('k')
plt.yscale('log')
plt.xscale('log')
plt.tight_layout()
plt.savefig('sum_plot.png')

#print 'cos residuals ', np.abs(cos_compare-cos_brute).max()
#print 'sin residuals ', np.abs(sin_compare-sin_brute).max()
#print 'n_freq ',len(freq_arr),delta
#for ff, cc1, cc2 in zip(freq_arr, cos_brute, cos_compare):
#    print ff,cc1,cc2

