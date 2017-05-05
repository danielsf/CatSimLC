import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from clean import get_ls
from PressRybicki import get_ls_PressRybicki

import numpy as np
import time

def ff_lc(time):
    f_arr = 0.4*np.cos(9.0*(time-91.3))
    f_arr += 0.6*np.sin(25.0*(time+26.0))
    f_arr += 0.4*np.cos(76.0*(time-15.0))
    #f_arr = 2.0*np.cos(5.0*(time-15.6))
    return f_arr


rng = np.random.RandomState(88)

time_arr = np.array([],dtype=float)
for mn, mx in zip(np.arange(0.0, 90.0, 10.0), np.arange(5.0, 95.0, 10.0)):
    sub_sample = rng.random_sample(10)*5.0+mn
    time_arr = np.append(time_arr, sub_sample)

time_arr = np.sort(time_arr)
f_arr = ff_lc(time_arr)

print('time arr max %e ' % time_arr.max())

time_control_arr = np.arange(0.0,100.0, 0.1)
f_control_arr = ff_lc(time_control_arr)

#sigma_arr = rng.random_sample(len(time_arr))*0.5+0.25
sigma_arr = rng.random_sample(len(time_arr))*0.3 + 0.1
noise_arr = np.array([rng.normal(0.0, ss) for ss in sigma_arr])

f_noisy_arr = f_arr + noise_arr


freq_arr = np.arange(0.001,4000.0,0.1)


t_start = time.time()
(pp_test, freq_arr, tau_test,
aa_test, bb_test, cc_test) = get_ls_PressRybicki(time_arr, f_noisy_arr,
                                                 sigma_arr)

t_PR = time.time()-t_start


t_start = time.time()
pp,tau,aa,bb,cc = get_ls(time_arr-time_arr.min(), f_noisy_arr, sigma_arr, 2.0*np.pi*freq_arr)
t_brute = time.time()-t_start

plt.figsize = (30,30)
plt.subplot(3,1,1)
plt.title('brute force')
plt.plot(freq_arr, pp)
plt.xlabel('$\\nu$')
plt.subplot(3,1,2)
plt.title('Press and Rybicki')
plt.plot(freq_arr, pp_test)
plt.xlabel('$\\nu$')
plt.subplot(3,1,3)
plt.title("residual")
plt.plot(1.0/freq_arr, np.abs(pp_test/pp-1.0))
plt.xlabel('$\\tau$')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.savefig('periodogram_comparison.png')
plt.close()

plt.figsize=(30,30)
for ix, (title, control, test) in enumerate(zip(('c', 's', 'cs'),
                                      (get_ls._c_tau, get_ls._s_tau,
                                       get_ls._cs),
                                      (get_ls_PressRybicki.c,
                                       get_ls_PressRybicki.s,
                                       get_ls_PressRybicki.cs))):

    plt.subplot(3,1,ix+1)
    plt.title(title)
    plt.plot(1.0/freq_arr,np.abs(test/control-1.0))
    plt.yscale('log')
    plt.xscale('log')
    print title
    dexes = np.where(np.logical_and(np.logical_not(np.isnan(test)),
                                    np.logical_not(np.isnan(control))))
    test = test[dexes]
    control = control[dexes]
    print test.min(),test.max()
    print control.min(),control.max()
    print np.median(np.abs(test-control))
    """
    header_list = []
    label_list = []
    plt.title(title)
    hh, = plt.plot(freq_arr, control)
    header_list.append(hh)
    label_list.append('control')
    hh, = plt.plot(freq_arr, test)
    header_list.append(hh)
    label_list.append('PR')
    plt.legend(header_list, label_list, loc=0)
    """

plt.tight_layout()
plt.savefig('periodogram_components_1.png')
plt.close()


for ix, (title, control, test) in enumerate(zip(('cc', 'ss', 'd'),
                                      (get_ls._cc, get_ls._ss,
                                       get_ls._d),
                                      (get_ls_PressRybicki.cc,
                                       get_ls_PressRybicki.ss,
                                       get_ls_PressRybicki.d))):

    plt.subplot(3,1,ix+1)
    plt.title(title)
    plt.plot(1.0/freq_arr,np.abs(test/control-1.0))
    plt.yscale('log')
    plt.xscale('log')
    print title
    dexes = np.where(np.logical_and(np.logical_not(np.isnan(test)),
                                    np.logical_not(np.isnan(control))))
    test = test[dexes]
    control = control[dexes]
    print test.min(),test.max()
    print control.min(),control.max()
    print np.median(np.abs(test-control))
    """
    header_list = []
    label_list = []
    plt.title(title)
    hh, = plt.plot(freq_arr, control)
    header_list.append(hh)
    label_list.append('control')
    hh, = plt.plot(freq_arr, test)
    header_list.append(hh)
    label_list.append('PR')
    plt.legend(header_list, label_list, loc=0)
    """

plt.tight_layout()
plt.savefig('periodogram_components_2.png')

print 'brute ',t_brute
print 'PR ',t_PR
dex = np.argmin(np.abs(freq_arr*2.0*np.pi-25.0))
print 'closest freq ',freq_arr[dex]
sorted_dexes = np.argsort(pp)
print('max classical freq %.4f %.4f %.4f %.4f' %
      (freq_arr[sorted_dexes[-1]]*2.0*np.pi,
         freq_arr[sorted_dexes[-2]]*2.0*np.pi,
         freq_arr[sorted_dexes[-3]]*2.0*np.pi,
         freq_arr[sorted_dexes[-4]]*2.0*np.pi))

sorted_dexes = np.argsort(pp_test)
print ('PR freqs: %.4f %.4f %.4f %.4f' %
        (freq_arr[sorted_dexes[-1]]*2.0*np.pi,
         freq_arr[sorted_dexes[-2]]*2.0*np.pi,
         freq_arr[sorted_dexes[-3]]*2.0*np.pi,
         freq_arr[sorted_dexes[-4]]*2.0*np.pi))
