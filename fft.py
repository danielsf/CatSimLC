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
    active_val = 2**(num_bits-1)
    additive_val = 1
    for i_bit in range(num_bits):
        if in_val & active_val != 0:
            out += additive_val
        additive_val *= 2
        active_val = active_val//2
    return out


def _bit_reverse_vector(in_val, num_bits):

    out = np.zeros(len(in_val), dtype=int)
    active_val = 2**(num_bits-1)
    additive_val = 1
    for i_bit in range(num_bits):
        dexes = np.where(in_val & active_val != 0)
        out[dexes] += additive_val
        additive_val *= 2
        active_val = active_val//2
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
    forward_dexes = np.arange(len(f_arr), dtype=int)
    rev_dexes = _bit_reverse_vector(forward_dexes, n_bits)
    fft_re[rev_dexes] = f_arr

    n_pts = 1
    n_strides = len(f_arr)
    tot_pts = len(time_arr)
    cache_dexes = np.arange(tot_pts//2)

    if (not hasattr(fft_real, 'cos_cache') or
        not np.array_equal(cache_dexes, fft_real.cache_dexes)):

        fft_real.cache_dexes = copy.deepcopy(cache_dexes)

        fft_real.cos_cache = np.cos(2.0*np.pi*cache_dexes/tot_pts)
        fft_real.sin_cache = np.sin(2.0*np.pi*cache_dexes/tot_pts)

    print 'prep took ',time.time()-t_start
    t_start = time.time()

    for i_bit in range(n_bits):
        n_pts *= 2
        n_strides = n_strides//2
        base_dexes = np.arange(0, n_strides*n_pts, n_pts)
        even_dexes = np.array([base_dexes + k for k in range(n_pts//2)]).flatten()
        odd_dexes = even_dexes + n_pts//2
        cache_dexes = np.array([[k*tot_pts//n_pts]*len(base_dexes) for k in range(n_pts//2)]).flatten()
        w_re = fft_real.cos_cache[cache_dexes]
        w_im = fft_real.sin_cache[cache_dexes]
        temp_re_even = fft_re[even_dexes]
        temp_im_even = fft_im[even_dexes]
        temp_re_odd = fft_re[odd_dexes]*w_re - fft_im[odd_dexes]*w_im
        temp_im_odd = fft_im[odd_dexes]*w_re + fft_re[odd_dexes]*w_im
        fft_re[even_dexes] += temp_re_odd
        fft_im[even_dexes] += temp_im_odd
        fft_re[odd_dexes] = temp_re_even - temp_re_odd
        fft_im[odd_dexes] = temp_im_even - temp_im_odd

    print 'work took ',time.time()-t_start
    return fft_re*delta, fft_im*delta
