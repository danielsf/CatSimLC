"""
This script will implement the CLEAN algorithm as presented in
Roberts et al. 1987 (AJ 93, 968)
"""

import numpy as np
import copy
import time


def get_window_function(time_arr, freq_arr):
    """
    Parameters
    ----------
    time_arr is a numpy array containing the times at which
    the data was sampled

    freq_arr is a numpy array containing the grid of frequencies
    (inverse units as time_arr) being considered

    Returns
    -------
    A numpy array representing the Fourier transform of the window
    function, W(nu) from equation (8) of

    Robert et al. 1987 (AJ 93, 968)

    The first row will be the real part; the second row will be
    the imaginary part.
    """

    npts = len(time_arr)
    real_part = np.array([np.cos(2.0*np.pi*nu*time_arr).sum()
                          for nu in freq_arr])

    im_part = np.array([-1.0*np.sin(2.0*np.pi*nu*time_arr).sum()
                        for nu in freq_arr])

    return np.array([real_part/npts, im_part/npts])


def get_ls(time_arr, f_arr, sigma_arr, freq_arr):
    """
    Evaluate the generalized Lomb-Scargle periodogram for a
    ligth curve, as in

    Zechmeister and Kurster 2009 (A&A 496, 577)

    Parameters
    ----------
    time_arr is a numpy array of the times at which the light curve
    is sampled

    f_arr is a numpy array of the values of the light curve

    freq_arr is a numpy array of the angular frequencies at which to
    evaluate the periodogram

    Returns
    -------
    a numpy array of the periodogram (equation 20 of Zechmeister and
    Kurster 2009)

    a numpy array of the time offsets tau (equation 19 of
    Zechmeister and Kurster 2009)
    """

    if (not hasattr(get_ls, '_freq_cache') or
        not hasattr(get_ls, '_time_cache') or
        not hasattr(get_ls, '_sigma_cache') or
        not np.array_equal(freq_arr, get_ls._freq_cache) or
        not np.array_equal(time_arr, get_ls._time_cache) or
        not np.array_equal(sigma_arr, get_ls._sigma_cache)):

        get_ls._time_cache = copy.deepcopy(time_arr)
        get_ls._freq_cache = copy.deepcopy(freq_arr)
        get_ls._sigma_cache = copy.deepcopy(sigma_arr)

        _w = (1.0/np.power(sigma_arr, 2)).sum()
        get_ls._wgt_cache = 1.0/(_w*np.power(sigma_arr, 2))

        _c = np.dot(np.array([np.cos(omega*time_arr)
                              for omega in freq_arr]),
                    get_ls._wgt_cache)


        _s = np.dot(np.array([np.sin(omega*time_arr)
                              for omega in freq_arr]),
                    get_ls._wgt_cache)

        _omega_tau = np.arctan2(2*_c*_s, _c*_c-_s*_s)
        _tau = _omega_tau/(2.0*freq_arr)

        del _s
        del _c

        get_ls._cos_tau = np.array([np.cos(omega*(time_arr-tt))
                                    for omega, tt in zip(freq_arr, _tau)])

        get_ls._sin_tau = np.array([np.sin(omega*(time_arr-tt))
                                    for omega, tt in zip(freq_arr, _tau)])

        _cchat = np.dot(get_ls._cos_tau*get_ls._cos_tau,
                        get_ls._wgt_cache)
        _sshat = np.dot(get_ls._sin_tau*get_ls._sin_tau,
                        get_ls._wgt_cache)
        _cshat = np.dot(get_ls._cos_tau*get_ls._sin_tau,
                        get_ls._wgt_cache)

        get_ls._c_tau = np.dot(get_ls._cos_tau, get_ls._wgt_cache)
        get_ls._s_tau = np.dot(get_ls._sin_tau, get_ls._wgt_cache)

        get_ls._cc = _cchat - get_ls._c_tau*get_ls._c_tau
        get_ls._ss = _sshat - get_ls._s_tau*get_ls._s_tau

        get_ls._cs = _cshat - get_ls._c_tau*get_ls._s_tau
        get_ls._d = get_ls._cc*get_ls._ss - get_ls._cs*get_ls._cs

        get_ls._tau = _tau

        assert len(get_ls._cc) == len(freq_arr)
        assert len(get_ls._ss) == len(freq_arr)

    t_start = time.time()
    y = (get_ls._wgt_cache*f_arr).sum()
    yy = (get_ls._wgt_cache*np.power(f_arr-y,2)).sum()
    yc = np.dot(get_ls._cos_tau, get_ls._wgt_cache*(f_arr-y))
    ys = np.dot(get_ls._sin_tau, get_ls._wgt_cache*(f_arr-y))

    aa = (yc*get_ls._ss-ys*get_ls._cs)/get_ls._d
    bb = (ys*get_ls._cc-yc*get_ls._cs)/get_ls._d

    cc = y-aa*get_ls._c_tau-bb*get_ls._s_tau

    pgram = ((yc*yc/get_ls._cc) + (ys*ys/get_ls._ss))/yy
    print('repetitive part took %e' % (time.time()-t_start))
    return pgram, get_ls._tau, aa, bb, cc


def get_clean_spectrum(time_arr, f_arr, sigma_arr, freq_arr):
    """
    Clean a time series according to the algorithm presented in
    Roberts et al. 1987 (AJ 93, 968) (though this works in real
    space)

    Will return parameters needed to reconstruct a clean version
    of the light curve as

    \sum_i a_i cos(omega_i (t-tau_i)) + b_i sin(omega_i (t-tau_i)) + c_i

    Parameters
    ----------
    time_arr is a numpy array of when the light curve was sampled

    f_arr is a numpy array of light curve flux/magnitude values

    sigma_arr is a numpy array of uncertainties on f_arr

    freq_arr is a numpy array of the angular frequencies to be
    meaured

    Returns
    -------
    aa a numpy array of a_i parameters from the model

    bb a numpy array of b_i parameters from the model

    cc a numpy array of c_i parameters from the model

    omega a numpy array of omega_i parameters from the model

    tau a numpy array of tau_i parameters from the model
    """

    iteration = 10
    gain = 1.0

    window = get_window_function(time_arr, freq_arr)

    residual_arr = copy.deepcopy(f_arr)

    pspec, tau, aa, bb, cc = get_ls(time_arr, residual_arr,
                                    sigma_arr, freq_arr)

    aa_list = []
    bb_list = []
    cc_list = []
    tau_list = []
    omega_list = []

    for it in range(1, iteration+1):
        max_dex = np.argmax(pspec)
        freq_max = freq_arr[max_dex]
        tau_max = tau[max_dex]
        aa_max = aa[max_dex]*gain
        bb_max = bb[max_dex]*gain
        cc_max = cc[max_dex]*gain

        aa_list.append(aa_max)
        bb_list.append(bb_max)
        cc_list.append(cc_max)
        tau_list.append(tau_max)
        omega_list.append(freq_max)

        model = aa_max*np.cos(freq_max*(time_arr-tau_max))
        model += bb_max*np.sin(freq_max*(time_arr-tau_max))
        model += cc_max

        residual_arr -= model
        if it<iteration:
            pspec, tau, aa, bb, cc = get_ls(time_arr, residual_arr,
                                            sigma_arr, freq_arr)

    return (np.array(aa_list), np.array(bb_list),
            np.array(cc_list), np.array(omega_list),
            np.array(tau_list))
