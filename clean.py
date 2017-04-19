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
        not np.array_equal(freq_arr, get_ls._freq_cache) or
        not np.array_equal(time_arr, get_ls._time_cache)):

        get_ls._time_cache = copy.deepcopy(time_arr)
        get_ls._freq_cache = copy.deepcopy(freq_arr)
        get_ls._cos_cache = np.array([np.cos(omega*time_arr)
                                      for omega in freq_arr])

        get_ls._sin_cache = np.array([np.sin(omega*time_arr)
                                      for omega in freq_arr])

    if (not hasattr(get_ls, '_sigma_cache') or
        not np.array_equal(sigma_arr, get_ls._sigma_cache)):

        get_ls._sigma_cache = copy.deepcopy(sigma_arr)
        _w = (1.0/np.power(sigma_arr, 2)).sum()
        get_ls._wgt_cache = 1.0/(_w*np.power(sigma_arr, 2))
        _c = np.dot(get_ls._cos_cache, get_ls._wgt_cache)
        _s = np.dot(get_ls._sin_cache, get_ls._wgt_cache)

        _omega_tau = np.arctan2(2*_c*_s, _c*_c-_s*_s)
        _tau = _omega_tau/(2.0*freq_arr)

        get_ls._cos_tau = np.array([np.cos(omega*(time_arr-tt))
                                    for omega, tt in zip(freq_arr, _tau)])

        get_ls._sin_tau = np.array([np.sin(omega*(time_arr-tt))
                                    for omega, tt in zip(freq_arr, _tau)])

        get_ls._sin2 = get_ls._sin_tau*get_ls._sin_tau
        get_ls._cos2 = get_ls._cos_tau*get_ls._cos_tau

        _cchat = np.dot(get_ls._cos_tau*get_ls._cos_tau,
                        get_ls._wgt_cache)
        _sshat = np.dot(get_ls._sin_tau*get_ls._sin_tau,
                        get_ls._wgt_cache)
        _cshat = np.dot(get_ls._cos_tau*get_ls._sin_tau,
                        get_ls._wgt_cache)

        _c_tau = np.dot(get_ls._cos_tau, get_ls._wgt_cache)
        _s_tau = np.dot(get_ls._sin_tau, get_ls._wgt_cache)

        get_ls._cc = _cchat - _c_tau*_c_tau
        get_ls._ss = _sshat - _s_tau*_s_tau
        get_ls._tau = _tau

    assert len(get_ls._cc) == len(freq_arr)
    assert len(get_ls._ss) == len(freq_arr)

    t_start = time.time()
    y = (get_ls._wgt_cache*f_arr).sum()
    yy = (get_ls._wgt_cache*np.power(f_arr-y,2)).sum()
    yc = np.dot(get_ls._cos_tau, get_ls._wgt_cache*(f_arr-y))
    ys = np.dot(get_ls._sin_tau, get_ls._wgt_cache*(f_arr-y))

    pgram = ((yc*yc/get_ls._cc) + (ys*ys/get_ls._ss))/yy
    print('repetitive part took %e' % (time.time()-t_start))
    return pgram, get_ls._tau
