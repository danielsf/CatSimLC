import numpy as np
import hashlib
import copy
from fft import fft_real


def _do_extirpation(hk, ff, tt, ttk, dexes):

    for ix in range(len(dexes)):
        other_dexes = dexes[np.where(dexes!=dexes[ix])]
        term = np.product(tt-ttk[other_dexes])
        term /= np.product(ttk[dexes[ix]]-ttk[other_dexes])
        term *= ff
        hk[dexes[ix]] += term


def extirp_sums(tt_arr, ff_arr, delta, n_t):
    """
    Take an arbitrary function sampled irregularly and return the
    sums in equation (5) of Press and Rybicki 1988 (ApJ 338, 277)
    using the FFT code in this module.
    """

    ttk = np.arange(0.0,
                    n_t*delta,
                    delta)

    print('actual len(ttk) %d' % len(ttk))
    time_dexes = np.round((tt_arr-ttk.min())/delta).astype(int)
    hk = np.zeros(len(ttk))
    dexes = np.zeros(25, dtype=int)
    half_dexes = len(dexes)//2
    for ij, (tt, ff) in enumerate(zip(tt_arr, ff_arr)):
        tj=time_dexes[ij]
        if tj<=half_dexes+1:
            for ix in range(len(dexes)):
                dexes[ix] = ix
        elif tj>=len(ttk)-(half_dexes+1):
            for ix in range(len(dexes)):
                dexes[ix] = len(ttk)-ix-1
        else:
            for ix in range(len(dexes)):
                dexes[ix] = tj - half_dexes + ix

        _do_extirpation(hk, ff, tt, ttk, dexes)

    print 'max hk ',np.abs(hk).max(),ff_arr.max()
    ft_re, ft_im = fft_real(ttk, hk)

    return (ft_re/delta, ft_im/delta, ttk, hk)


def _initialize_PressRybicki(time_arr, sigma_arr, delta):
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
    n_t = the number of steps in the frequency_array used for FFT
    freq_arr = the array of frequencies (not angular frequencies)
    cut_off_freq = frequency at which Press and Rybicki approximation breaks down
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

    n_t_init = time_arr.max()/delta
    n_t = 2
    while n_t < n_t_init:
        n_t *= 2
    print('n_t %d\ndelta %e' % (int(n_t), delta))

    freq_arr = np.array([k/(delta*n_t) for k in range(n_t)])

    w = (1.0/np.power(sigma_arr, 2)).sum()
    wgt_fn = 1.0/(w*np.power(sigma_arr, 2))

    c, s, tk, hk = extirp_sums(time_arr, wgt_fn, delta, n_t)
    del tk
    del hk

    print 'testing validity of PR trick'
    cos_test = np.dot(np.array([np.cos(2.0*np.pi*nu*time_arr)
                                for nu in freq_arr]),
                      wgt_fn)

    threshold = 0.01
    offending_dex = np.where(np.logical_or(np.logical_and(np.abs(cos_test)<0.001,
                                                          np.abs(c-cos_test)>threshold),
                                           np.abs(c/cos_test-1.0)>threshold))

    worst_cos = np.min(offending_dex[0])
    del cos_test

    sin_test = np.dot(np.array([np.sin(2.0*np.pi*nu*time_arr)
                                for nu in freq_arr]),
                      wgt_fn)


    offending_dex = np.where(np.logical_or(np.logical_and(np.abs(sin_test)<0.001,
                                                          np.abs(s-sin_test)>threshold),
                                           np.abs(s/sin_test-1.0)>threshold))

    worst_sin = np.min(offending_dex[0])
    del sin_test

    if worst_cos<worst_sin:
        cut_off_freq = 0.5*freq_arr[worst_cos]
    else:
        cut_off_freq = 0.5*freq_arr[worst_sin]

    print 'done assessing validity'

    c2_raw, s2_raw, tk, hk = extirp_sums(2.0*time_arr, wgt_fn, delta, n_t*2)
    del tk
    del hk
    dexes = range(0,len(c2_raw), 2)
    c2 = c2_raw[dexes]
    del c2_raw
    s2 = s2_raw[dexes]
    del s2_raw

    omega_tau = np.arctan2(s2-2*c*s, c2-c*c+s*s)
    tau = omega_tau/(4.0*np.pi*freq_arr)

    cos_omega_tau = np.cos(2.0*np.pi*freq_arr*tau)
    sin_omega_tau = np.sin(2.0*np.pi*freq_arr*tau)
    cos_tau = c*cos_omega_tau + s*sin_omega_tau
    sin_tau = s*cos_omega_tau - c*sin_omega_tau

    del c
    del s

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
    # n_t
    # freq_arr
    # cut_off_freq
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
    return (wgt_fn, delta, n_t, freq_arr, cut_off_freq, tau,
            cos_omega_tau, sin_omega_tau,
            cos_tau, sin_tau, cc, ss, cs, d)

def get_ls_PressRybicki(time_arr_in, f_arr_in, sigma_arr_in, delta):
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

    delta is a float; maximum frequency considered will be 1/delta

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

    """
    if hasattr(get_ls_PressRybicki, 'initialized'):
        if get_ls_PressRybicki.initialized:
            local_hasher = hashlib.sha1()
            local_hasher.update(time_arr_in)
            if local_hasher.hexdigest() != get_ls_PressRybicki.time_hash.hexdigest():
                get_ls_PressRybicki.initialized = False

            local_hasher = hashlib.sha1()
            local_hasher.update(sigma_arr_in)
            if local_hasher.hexdigest() != get_ls_PressRybicki.sigma_hash.hexdigest():
                get_ls_PressRybicki.initialized = False

            if delta != get_ls_PressRybicki.delta:
                get_ls_PressRybicki.initialized = False
    """

    if (not hasattr(get_ls_PressRybicki, 'initialized') or
        not get_ls_PressRybicki.initialized):

        print '\n\ninitializing periodogram\n\n'

        get_ls_PressRybicki.initialized = True

        #get_ls_PressRybicki.time_hash = hashlib.sha1()
        #get_ls_PressRybicki.time_hash.update(time_arr_in)

        #get_ls_PressRybicki.sigma_hash = hashlib.sha1()
        #get_ls_PressRybicki.sigma_hash.update(sigma_arr_in)

        (get_ls_PressRybicki.w,
         get_ls_PressRybicki.delta,
         get_ls_PressRybicki.n_t,
         get_ls_PressRybicki.freq_arr,
         get_ls_PressRybicki.cut_off_freq,
         get_ls_PressRybicki.tau,
         get_ls_PressRybicki.cos_omega_tau,
         get_ls_PressRybicki.sin_omega_tau,
         get_ls_PressRybicki.c,
         get_ls_PressRybicki.s,
         get_ls_PressRybicki.cc,
         get_ls_PressRybicki.ss,
         get_ls_PressRybicki.cs,
         get_ls_PressRybicki.d) = _initialize_PressRybicki(time_arr, sigma_arr, delta)

    y_bar = (f_arr*get_ls_PressRybicki.w).sum()
    yy = (get_ls_PressRybicki.w*np.power(f_arr-y_bar,2)).sum()
    y_fn = get_ls_PressRybicki.w*(f_arr-y_bar)
    y_c_raw, y_s_raw, tk, hk = extirp_sums(time_arr, y_fn,
                                           get_ls_PressRybicki.delta,
                                           get_ls_PressRybicki.n_t)
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


def _is_significant(aa, bb, cc, omega, tau,
                    aa_test, bb_test, cc_test, omega_test, tau_test,
                    time_arr, f_arr, sig_arr):

    if (not hasattr(_is_significant, 'model') or
        _is_significant.model is None):

        model = np.zeros(len(time_arr))
        for a,b,c,o,t in zip(aa, bb, cc, omega, tau):
            model += c
            model += a*np.cos(o*(time_arr-time_arr.min()-t))
            model += b*np.sin(o*(time_arr-time_arr.min()-t))

        chi_0 = np.power((f_arr-model)/sig_arr,2).sum()
        _is_significant.bic_0 = 3.0*len(aa)*np.log(len(time_arr)) + chi_0
        _is_significant.model = model

    _is_significant.model += cc_test
    _is_significant.model += aa_test*np.cos(omega_test*(time_arr-time_arr.min()-tau_test))
    _is_significant.model += bb_test*np.sin(omega_test*(time_arr-time_arr.min()-tau_test))

    chi_1 = np.power((f_arr-_is_significant.model)/sig_arr,2).sum()
    bic_1 = 3.0*(len(aa)+1)*np.log(len(time_arr)) + chi_1

    print bic_1,_is_significant.bic_0

    if bic_1 < _is_significant.bic_0:
        _is_significant.bic_0 = bic_1
        return True
    return False

def get_clean_spectrum_PressRybicki(time_arr, f_arr, sigma_arr, delta):
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

    delta is a float; maximum frequency considered will be 1/delta

    Returns
    -------
    aa a numpy array of a_i parameters from the model

    bb a numpy array of b_i parameters from the model

    cc a numpy array of c_i parameters from the model

    omega a numpy array of omega_i parameters from the model

    tau a numpy array of tau_i parameters from the model
    """

    _is_significant.model = None

    iteration = 30
    gain = 1.0

    residual_arr = copy.deepcopy(f_arr)

    (pspec, freq_arr,
     tau, aa, bb, cc) = get_ls_PressRybicki(time_arr, residual_arr, sigma_arr, delta)

    aa_list = []
    bb_list = []
    cc_list = []
    tau_list = []
    omega_list = []

    for it in range(1, iteration+1):
        valid = np.where(np.logical_and(np.logical_not(np.isnan(pspec)),
                                        freq_arr<get_ls_PressRybicki.cut_off_freq))
        pspec = pspec[valid]
        aa = aa[valid]
        bb = bb[valid]
        cc = cc[valid]
        valid_freq = freq_arr[valid]
        tau = tau[valid]
        max_dex = np.argmax(pspec)

        freq_max = valid_freq[max_dex]
        tau_max = tau[max_dex]
        aa_max = aa[max_dex]*gain
        bb_max = bb[max_dex]*gain
        cc_max = cc[max_dex]*gain
        omega_max = freq_max*2.0*np.pi

        print 'a ',aa_max
        print 'b ',bb_max
        print 'c ',cc_max
        print 'omega ',freq_max*2.0*np.pi

        if _is_significant(aa_list, bb_list, cc_list, omega_list, tau_list,
                           aa_max, bb_max, cc_max, omega_max, tau_max,
                           time_arr, f_arr, sigma_arr):

            aa_list.append(aa_max)
            bb_list.append(bb_max)
            cc_list.append(cc_max)
            tau_list.append(tau_max)
            omega_list.append(omega_max)

        else:
            print 'is not significant'
            break

        model = np.array([cc_max]*len(time_arr))
        model += aa_max*np.cos(omega_max*(time_arr-time_arr.min()-tau_max))
        model += bb_max*np.sin(omega_max*(time_arr-time_arr.min()-tau_max))

        residual_arr -= model
        if it<iteration:
            (pspec, freq_arr,
             tau, aa, bb, cc) = get_ls_PressRybicki(time_arr, residual_arr, sigma_arr, delta)

    print'cut off is ',get_ls_PressRybicki.cut_off_freq*2.0*np.pi,get_ls_PressRybicki.cut_off_freq

    return (np.array(aa_list), np.array(bb_list),
            np.array(cc_list), np.array(omega_list),
            np.array(tau_list), freq_arr)
