import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from clean import get_ls, get_window_function

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


freq_arr = np.arange(0.001,4000.0,0.01)

window = get_window_function(time_arr, freq_arr)
pp,tau,aa,bb,cc = get_ls(time_arr, f_noisy_arr, sigma_arr, freq_arr)

dex = np.argmax(pp)
print pp[dex]
print freq_arr[dex]
print tau[dex]

dex_list = np.argsort(pp)
print 'next freq ',freq_arr[dex_list[-2]],pp[dex_list[-2]]

from clean import get_clean_spectrum
t_start = time.time()
(aa, bb, cc,
 omega, tau) = get_clean_spectrum(time_arr, f_noisy_arr,
                                  sigma_arr, freq_arr)

print 'cleaning took ',time.time()-t_start

model = np.zeros(len(time_control_arr))
for a, b, c, o, t in zip(aa, bb, cc, omega, tau):
    print('model omega %e %e %e %e %e' % (o,np.sqrt(a*a+b*b),a,b,c))
    model += c
    model += a*np.cos(o*(time_control_arr-t))
    model += b*np.sin(o*(time_control_arr-t))


plt.figsize=(30,30)
plt.subplot(2,1,1)
plt.plot(freq_arr, pp)
plt.subplot(2,1,2)
h_control, = plt.plot(time_control_arr,f_control_arr,zorder=1)
h_model, = plt.plot(time_control_arr, model, zorder=2)
plt.scatter(time_arr,f_noisy_arr,color='r',edgecolor='',zorder=3)
plt.legend([h_control, h_model], ['truth', 'model'], loc=0)
plt.xlim(0,15)
plt.savefig('test_clean_periodogram.png')
print('time arr max %e ' % time_arr.max())
