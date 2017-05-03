"""
Implement a Fast Fourier Transform as described in section 12.2.1 of

Numerical Recipes, Third Edition
by Press, Teukolsky, Vetterling, and Flannery
2007, Cambridge University Press
"""

import numpy as np
import copy

def _bit_reverse(in_val, num_bits):
    """
    Return the bit-reversed equivalent of the integer ii
    """
    out = 0
    for i_bit in range(num_bits):
        active_bit = num_bits-1-i_bit
        if in_val & 2**active_bit != 0:
            out += 2**i_bit
    return out

import time
def fft_real(time_arr, f_arr):
    """
    Fast Fourier Transform a real function as described in section 12.2.1
    of Numerical Recipes, Third Edition

    Parameters
    ----------
    time_arr is an evenly-spaced numpy array of time values

    f_arr_in is the corresponding numpy array of function values

    Returns
    -------
    """

    t_start = time.time()
    fft_re = np.zeros(len(f_arr))
    fft_im = np.zeros(len(f_arr))

    if len(time_arr) & (len(time_arr)-1) != 0:
        raise RuntimeError("FFT input arrays must have a power of 2 "
                           "number of elements; you gave %d" % len(time_arr))

    n_bits = int(np.log(len(time_arr))/np.log(2.0))
    delta = time_arr[1]-time_arr[0]
    f_arr_reversed = np.zeros(len(f_arr))
    for ii in range(len(f_arr)):
        opp = _bit_reverse(ii, n_bits)
        fft_re[opp] = f_arr[ii]

    print 'prep took ',time.time()-t_start
    t_start = time.time()
    n_pts = 1
    for i_bit in range(n_bits):
       n_pts *= 2
       n_strides = len(f_arr)//n_pts
       for k in range(n_pts//2):
          w_re = np.cos(2.0*np.pi*k/n_pts)
          w_im = np.sin(2.0*np.pi*k/n_pts)
          for i_stride in range(n_strides):
              i_even = i_stride*n_pts + k
              i_odd = i_even + n_pts//2
              temp_re_even = fft_re[i_even]
              temp_im_even = fft_im[i_even]
              temp_re_odd = fft_re[i_odd]*w_re - fft_im[i_odd]*w_im
              temp_im_odd = fft_im[i_odd]*w_re + fft_re[i_odd]*w_im
              fft_re[i_even] += temp_re_odd
              fft_im[i_even] += temp_im_odd
              fft_re[i_odd] = temp_re_even - temp_re_odd
              fft_im[i_odd] = temp_im_even - temp_im_odd

    print 'work took ',time.time()-t_start
    return fft_re*delta, fft_im*delta
