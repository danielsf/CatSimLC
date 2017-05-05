import numpy as np
import hashlib
from fft import fft_real


def extirp_sums(tt_arr, ff_arr, delta):
    """
    Take an arbitrary function sampled irregularly and return the
    sums in equation (5) of Press and Rybicki 1988 (ApJ 338, 277)
    using the FFT code in this module.
    """

    ttk = np.arange(0.0,
                    tt_arr.max(),
                    delta)

    print('actual len(ttk) %d' % len(ttk))
    time_dexes = np.round((tt_arr-ttk.min())/delta).astype(int)
    hk = np.zeros(len(ttk))
    for ij, (tt, ff) in enumerate(zip(tt_arr, ff_arr)):
        tj=time_dexes[ij]
        if tj==0:
            dexes = [0, 1, 2]
        elif tj>=len(ttk)-1:
            dexes = [len(ttk)-1, len(ttk)-2, len(ttk)-3]
        else:
            dexes = [tj-1, tj, tj+1]

        hk[dexes[0]] += ff*(tt-ttk[dexes[1]])*(tt-ttk[dexes[2]])/((ttk[dexes[0]]-ttk[dexes[1]])*(ttk[dexes[0]]-ttk[dexes[2]]))
        hk[dexes[1]] += ff*(tt-ttk[dexes[0]])*(tt-ttk[dexes[2]])/((ttk[dexes[1]]-ttk[dexes[0]])*(ttk[dexes[1]]-ttk[dexes[2]]))
        hk[dexes[2]] += ff*(tt-ttk[dexes[0]])*(tt-ttk[dexes[1]])/((ttk[dexes[2]]-ttk[dexes[0]])*(ttk[dexes[2]]-ttk[dexes[1]]))

    print 'max hk ',np.abs(hk).max(),ff_arr.max()
    ft_re, ft_im = fft_real(ttk, hk)

    return (ft_re/delta, ft_im/delta, ttk, hk)


def _initialize_PressRybicki(time_arr, sigma_arr):
    """
    Initialize arrays for the Press and Rybicki periodogram
    that only depend on t and sigma

    Parameters
    ----------
    time_arr is a numpy array of the times of the measurement

    sigma_arr is a numpy array of the uncertainties of the
    measurement

    Returns
    -------

    wgt = the array of weights (w_i from Zechmeister and Kurster 2009)
    delta = the delta_t in used to construct the frequency array
    freq_arr = the array of frequencies (not angular frequencies)
    tau = the array of time offsets tau from Zechmeister and Kurster 2009 eqn (19)
    cos_omega_tau = an array of cos(omega*tau)
    sin_omega_tau = an array of sin(omega*tau)
    an array of \sum w_i cos(omega (t_i-tau))
    an array of \sum w_i sin(omega (t_i-tau))
    an array of \sum w_i cos(omega (t_i-tau)) sin(omega (t_i-tau))
    an array of CC from Zechmeister and Kurster eqn 13
    an array of SS from Zechmeister and Kurster eqn 14
    an array of CS from Zechmeister and Kurster eqn 15
    an array of D from Zechmeister and Kurster eqn 6
    """

    delta_guess = 0.1*np.diff(time_arr).min()
    delta = time_arr.max()/1024.0
    n_t = 1024
    while delta>delta_guess:
        delta *= 0.5
        n_t *= 2
    print('n_t %d\ndelta %e' % (int(n_t), delta))

    freq_arr = np.array([k/(delta*n_t) for k in range(n_t)])

    w = (1.0/np.power(sigma_arr, 2)).sum()
    wgt_fn = 1.0/(w*np.power(sigma_arr, 2))

    c, s, tk, hk = extirp_sums(time_arr, wgt_fn, delta)
    del tk
    del hk

    omega_tau = np.arctan2(2*c*s, c*c-s*s)
    tau = omega_tau/(4.0*np.pi*freq_arr)

    cos_omega_tau = np.cos(2.0*np.pi*freq_arr*tau)
    sin_omega_tau = np.sin(2.0*np.pi*freq_arr*tau)
    cos_tau = c*cos_omega_tau + s*sin_omega_tau
    sin_tau = s*cos_omega_tau - c*sin_omega_tau

    del c
    del s

    c2_raw, s2_raw, tk, hk = extirp_sums(2.0*time_arr, wgt_fn, delta)
    del tk
    del hk
    dexes = range(0,len(c2_raw), 2)
    c2 = c2_raw[dexes]
    del c2_raw
    s2 = s2_raw[dexes]
    del s2_raw
    cos_2omega_tau = np.cos(4.0*np.pi*freq_arr*tau)
    sin_2omega_tau = np.sin(4.0*np.pi*freq_arr*tau)
    w_sum = wgt_fn.sum()

    csq = 0.5*w_sum + 0.5*c2*cos_2omega_tau + 0.5*s2*sin_2omega_tau
    ssq = 0.5*w_sum - 0.5*c2*cos_2omega_tau - 0.5*s2*sin_2omega_tau
    csomega = 0.5*(s2*cos_2omega_tau - c2*sin_2omega_tau)  # cos(theta)*sin(theta) = 0.5*sin(2*theta)

    del s2
    del c2

    cs = csomega - cos_tau*sin_tau
    ss = ssq - sin_tau*sin_tau
    cc = csq - cos_tau*cos_tau
    d = cc*ss - cs*cs

    # return:
    # wgt
    # delta
    # freq_arr
    # tau
    # cos(omega*tau)
    # sin(omega*tau)
    # \sum w_i cos(omega (t_i-tau))
    # \sum w_i sin(omega (t_i-tau))
    # \sum w_i cos(omega (t_i-tau)) sin(omega (t_i-tau))
    # CC from Zechmeister and Kurster eqn 13
    # SS from Zechmeister and Kurster eqn 14
    # CS from Zechmeister and Kurster eqn 15
    # D from Zechmeister and Kurster eqn 6
    return (wgt_fn, delta, freq_arr, tau,
            cos_omega_tau, sin_omega_tau,
            cos_tau, sin_tau, cc, ss, cs, d)

def get_ls_PressRybicki(time_arr_in, f_arr_in, sigma_arr_in):
    """
    Evaluate the generalized Lomb-Scargle periodogram for a
    ligth curve, as in

    Zechmeister and Kurster 2009 (A&A 496, 577)

    Using the FFT trick from Press and Rybicki 1989 (ApJ 338, 277)
    to quickly evaluate the sums over sinusoids.

    Parameters
    ----------
    time_arr_in is a numpy array of the times at which the light curve
    is sampled

    f_arr_in is a numpy array of the values of the light curve

    freq_arr_in is a numpy array of the angular frequencies at which to
    evaluate the periodogram

    Returns
    -------
    a numpy array of the periodogram (equation 20 of Zechmeister and
    Kurster 2009)

    a numpy array of the time offsets tau (equation 19 of
    Zechmeister and Kurster 2009)
    """

    time_offset = time_arr_in.min()
    time_arr = time_arr_in - time_offset
    sorted_dex = np.argsort(time_arr)
    time_arr = time_arr[sorted_dex]
    f_arr = f_arr_in[sorted_dex]
    sigma_arr = sigma_arr_in[sorted_dex]

    if hasattr(get_ls_PressRybicki, 'initialized'):
        if get_ls_PressRybicki.initialized:
            local_hasher = hashlib.sha1()
            local_hasher.update(time_arr_in)
            if local_hasher.hexdigest() != get_ls_PressRybicki.time_hash.hexdigest():
                get_ls_PressRybicki.initialized = False
            local_hasher.update(f_arr_in)
            if local_hasher.hexdigest() != get_ls_PressRybicki.f_hash.hexdigest():
                get_ls_PressRybicki.initialized = False
            local_hasher.update(sigma_arr_in)
            if local_hasher.hexdigest() != get_ls_PressRybicki.sigma_hash.hexdigest():
                get_ls_PressRybicki.initialized = False

    if (not hasattr(get_ls_PressRybicki, 'initialized') or
        not get_ls_PressRybicki.initialized):

        get_ls_PressRybicki.initialized = True
        get_ls_PressRybicki.time_hash = hashlib.sha1()
        get_ls_PressRybicki.time_hash.update(time_arr_in)
        get_ls_PressRybicki.f_hash = hashlib.sha1()
        get_ls_PressRybicki.f_hash.update(f_arr_in)
        get_ls_PressRybicki.sigma_hash = hashlib.sha1()
        get_ls_PressRybicki.sigma_hash.update(sigma_arr_in)

        (get_ls_PressRybicki.w,
         get_ls_PressRybicki.delta,
         get_ls_PressRybicki.freq_arr,
         get_ls_PressRybicki.tau,
         get_ls_PressRybicki.cos_omega_tau,
         get_ls_PressRybicki.sin_omega_tau,
         get_ls_PressRybicki.c,
         get_ls_PressRybicki.s,
         get_ls_PressRybicki.cc,
         get_ls_PressRybicki.ss,
         get_ls_PressRybicki.cs,
         get_ls_PressRybicki.d) = _initialize_PressRybicki(time_arr, sigma_arr)

    y_bar = (f_arr*get_ls_PressRybicki.w).sum()
    yy = (get_ls_PressRybicki.w*np.power(f_arr-y_bar,2)).sum()
    y_fn = get_ls_PressRybicki.w*(f_arr-y_bar)
    y_c_raw, y_s_raw, tk, hk = extirp_sums(time_arr, y_fn, get_ls_PressRybicki.delta)
    del tk
    del hk

    y_c = (y_c_raw*get_ls_PressRybicki.cos_omega_tau +
           y_s_raw*get_ls_PressRybicki.sin_omega_tau)

    y_s = (y_s_raw*get_ls_PressRybicki.cos_omega_tau -
           y_c_raw*get_ls_PressRybicki.sin_omega_tau)

    del y_s_raw
    del y_c_raw

    aa = (y_c*get_ls_PressRybicki.ss - y_s*get_ls_PressRybicki.cs)/get_ls_PressRybicki.d
    bb = (y_s*get_ls_PressRybicki.cc - y_c*get_ls_PressRybicki.cs)/get_ls_PressRybicki.d
    cc = y_bar - aa*get_ls_PressRybicki.c - bb*get_ls_PressRybicki.s

    pgram = ((y_c*y_c/get_ls_PressRybicki.cc) + (y_s*y_s/get_ls_PressRybicki.ss))/yy

    return pgram, get_ls_PressRybicki.freq_arr, get_ls_PressRybicki.tau, aa, bb, cc
