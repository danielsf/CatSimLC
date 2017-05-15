import numpy as np

def brute_force_ft(time_arr, f_arr):
    delta = time_arr[1]-time_arr[0]
    ft_re = delta*np.array([(f_arr*np.cos(2.0*np.pi*k*time_arr/(delta*len(f_arr)))).sum()
                            for k in range(len(f_arr))])
    ft_im = delta*np.array([(f_arr*np.sin(2.0*np.pi*k*time_arr/(delta*len(f_arr)))).sum()
                            for k in range(len(f_arr))])

    return ft_re, ft_im

from fft import _bit_reverse
import copy

def fft_four(time_arr, ff_arr):
    if len(ff_arr)!=4:
        raise RuntimeError("need to give me a 4d array")
    delta = time_arr[1]-time_arr[0]
    fft_re = np.zeros(len(ff_arr))
    fft_im = np.zeros(len(ff_arr))
    ff_rev = np.zeros(len(ff_arr))
    for ii in range(len(ff_arr)):
        opp = _bit_reverse(ii, 2)
        print 'mapping ',ii,opp
        ff_rev[ii] = ff_arr[opp]

    fft2_re = copy.deepcopy(ff_rev)
    fft2_im = np.zeros(len(ff_rev))
    fft_temp_re = copy.deepcopy(fft2_re)
    fft_temp_im = copy.deepcopy(fft2_im)

    w_re = 1.0
    w_im = 0.0
    fft2_re[0] = fft_temp_re[0]+fft_temp_re[1]*w_re-fft_temp_im[1]*w_im
    fft2_im[0] = fft_temp_im[0]+fft_temp_re[1]*w_im+fft_temp_im[1]*w_re

    fft2_re[2] = fft_temp_re[2]+fft_temp_re[3]*w_re-fft_temp_im[3]*w_im
    fft2_im[2] = fft_temp_im[2]+fft_temp_re[3]*w_im+fft_temp_im[3]*w_re

    w_re = np.cos(2.0*np.pi/2)
    w_im = np.sin(2.0*np.pi/2)
    fft2_re[1] = fft_temp_re[0]+fft_temp_re[1]*w_re-fft_temp_im[1]*w_im
    fft2_im[1] = fft_temp_im[0]+fft_temp_re[1]*w_im+fft_temp_im[1]*w_re

    fft2_re[3] = fft_temp_re[2]+fft_temp_re[3]*w_re-fft_temp_im[3]*w_im
    fft2_im[3] = fft_temp_im[2]+fft_temp_re[3]*w_im+fft_temp_im[3]*w_re

    temp_re = copy.deepcopy(fft2_re)
    temp_im = copy.deepcopy(fft2_im)

    w_re = np.cos(0.0)
    w_im = np.sin(0.0)
    odd_re = temp_re[2]*w_re-temp_im[2]*w_im
    odd_im = temp_re[2]*w_im+temp_im[2]*w_re
    fft_re[0] = temp_re[0] + odd_re
    fft_im[0] = temp_im[0] + odd_im

    fft_re[2] = temp_re[0] - odd_re
    fft_im[2] = temp_im[0] - odd_im

    w_re = np.cos(2.0*np.pi/4)
    w_im = np.sin(2.0*np.pi/4)
    odd_re = temp_re[3]*w_re-temp_im[3]*w_im
    odd_im = temp_re[3]*w_im+temp_im[3]*w_re
    fft_re[1] = temp_re[1] + odd_re
    fft_im[1] = temp_im[1] + odd_im

    fft_re[3] = temp_re[1] - odd_re
    fft_im[3] = temp_im[1] - odd_im

    return fft_re*delta, fft_im*delta

def f_of_t(time_arr):
    y = np.cos(2.4*time_arr)
    y += 0.3*np.sin(8.7*time_arr)
    y += -5.5*np.sin(1.3*time_arr)
    y+= 0.3*np.cos(11.3*time_arr)
    y += -0.9*np.sin(0.462*time_arr)
    return y



dt = 0.037
time_arr = np.arange(0.0, 32*1024*dt, dt)
y_arr = f_of_t(time_arr)

from fft import FFTransformer
import time

ffter = FFTransformer()

t_start = time.time()
f_re, f_im = ffter.fft_real(time_arr, y_arr)
t_fft = time.time()-t_start

t_start = time.time()
brute_re, brute_im = brute_force_ft(time_arr, y_arr)
t_brute = time.time()-t_start

#for k in range(len(f_re)):
#    print k/(len(time_arr)*(time_arr[1]-time_arr[0])),f_re[k],brute_re[k],f_im[k],brute_im[k]
print 'max residuals: ',np.abs(f_re-brute_re).max(),np.abs(f_im-brute_im).max()
print 'median amplitudes: ',np.median(np.abs(brute_re)),np.median(np.abs(brute_im))
#print np.abs(brute_re).min(),np.abs(brute_im).min()
print 'time_fft  ',t_fft,'time_brute ',t_brute
print 'len ',len(f_re)
"""
print '\n'
print 'four'

re4, im4 = fft_four(time_arr, y_arr)

for k in range(len(f_re)):
    print k/(len(time_arr)*(time_arr[1]-time_arr[0])),re4[k],brute_re[k],im4[k],brute_im[k]

print np.abs(brute_re-re4).max(),np.abs(brute_im-im4).max()
"""
