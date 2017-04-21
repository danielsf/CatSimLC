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

def _f_k_real(f_arr, k, f_arr_reversed=None):

    n_bits = int(np.log(len(f_arr))/np.log(2.0))
    n_f = len(f_arr)

    if f_arr_reversed is None:
        f_arr_reversed = np.zeros(len(f_arr))
        for ii in range(len(f_arr)):
            opp = _bit_reverse(ii, n_bits)
            f_arr_reversed[opp] = f_arr[ii]

    w_real = np.ones(len(f_arr))
    w_im = np.zeros(len(f_arr))
    n_exp = 2*len(f_arr)
    for i_bit in range(n_bits):
        n_exp = n_exp/2
        base = 2**i_bit
        term_re = np.cos(2.0*np.pi*k/n_exp)
        term_im = np.sin(2.0*np.pi*k/n_exp)
        for ii in range(len(f_arr)):
            if ii & base != 0:
                old_re = w_real[ii]
                old_im = w_im[ii]
                w_real[ii] = old_re*term_re - old_im*term_im
                w_im[ii] = old_re*term_im + old_im*term_re

    f_k_real = 0.0
    f_k_im = 0.0
    for ii in range(len(f_arr)):
        f_k_real += f_arr[ii]*w_real[ii]
        f_k_im += f_arr[ii]*w_im[ii]

    return f_k_real, f_k_im, f_arr_reversed


def fft_real(time_arr, f_arr_in):
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

    f_arr = copy.deepcopy(f_arr_in)

    if len(time_arr) & (len(time_arr)-1) != 0:
        raise RuntimeError("FFT input arrays must have a power of 2 "
                           "number of elements")

    n_bits = int(np.log(len(time_arr))/np.log(2.0))+1
    delta = time_arr[1]-time_arr[0]
    f_k_real = np.zeros(len(time_arr))
    f_k_im = np.zeros(len(time_arr))
    f_arr_reversed = None
    for k in range(len(time_arr)):
        (f_k_real[k],
         f_k_im[k],
         f_arr_reversed) = _f_k_real(f_arr, k,
                                     f_arr_reversed=f_arr_reversed)

    return f_k_real*delta, f_k_im*delta
