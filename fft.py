"""
Implement a Fast Fourier Transform as described in section 12.2.1 of

Numerical Recipes, Third Edition
by Press, Teukolsky, Vetterling, and Flannery
2007, Cambridge University Press
"""

import numpy as np
import copy
import time


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


class FFTransformer(object):

    def __init__(self):
        pass

    def _initialize_fft_real(self, n_bits, n_strides, calc_dexes):

        self.cos_cache = np.cos(2.0*np.pi*calc_dexes/self.tot_pts)
        self.sin_cache = np.sin(2.0*np.pi*calc_dexes/self.tot_pts)

        n_pts = 1

        self.even_dexes = []
        self.odd_dexes = []
        self.cache_dexes = []
        for i_bit in range(n_bits):
            n_pts *= 2
            n_strides = n_strides//2
            base_dexes = np.arange(0, n_strides*n_pts, n_pts)
            even_dexes = np.array([base_dexes + k for k in range(n_pts//2)]).flatten()
            self.even_dexes.append(even_dexes)
            self.odd_dexes.append(even_dexes + n_pts//2)
            self.cache_dexes.append(np.array([[k*self.tot_pts//n_pts]*len(base_dexes)
                                               for k in range(n_pts//2)]).flatten())


    def fft_real(self, time_arr, f_arr):
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

        tot_pts = len(time_arr)

        if (not hasattr(self, 'cos_cache') or
            n_bits != self.n_bits or
            tot_pts != self.tot_pts):

            calc_dexes = np.arange(tot_pts//2)
            n_strides = len(f_arr)

            self.n_bits = n_bits
            self.tot_pts = tot_pts
            self._initialize_fft_real(n_bits, n_strides, calc_dexes)

        #print 'prep took ',time.time()-t_start
        t_start = time.time()

        for i_bit in range(n_bits):
            cache_dexes = self.cache_dexes[i_bit]
            odd_dexes = self.odd_dexes[i_bit]
            even_dexes = self.even_dexes[i_bit]
            w_re = self.cos_cache[cache_dexes]
            w_im = self.sin_cache[cache_dexes]
            temp_re_even = fft_re[even_dexes]
            temp_im_even = fft_im[even_dexes]
            temp_re_odd = fft_re[odd_dexes]*w_re - fft_im[odd_dexes]*w_im
            temp_im_odd = fft_im[odd_dexes]*w_re + fft_re[odd_dexes]*w_im
            fft_re[even_dexes] += temp_re_odd
            fft_im[even_dexes] += temp_im_odd
            fft_re[odd_dexes] = temp_re_even - temp_re_odd
            fft_im[odd_dexes] = temp_im_even - temp_im_odd

        #print 'work took ',time.time()-t_start
        return fft_re*delta, fft_im*delta
