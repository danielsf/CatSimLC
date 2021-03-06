import gc
import numpy as np
import copy
from fft import FFTransformer


def extirp_sums(tt_arr, ff_arr, delta, n_t, ffter, cache_fft=False):
    """
    Take an arbitrary function sampled irregularly and return the
    sums in equation (5) of Press and Rybicki 1988 (ApJ 338, 277)
    using the FFT code in this module.

    ffter is an instantiation of FFTransformer

    cache_fft is a boolean.  If true, use memory intensive caches
    to speed up the FFT.
    """
    ttk = np.arange(0.0,
                    n_t*delta,
                    delta)
    hk = np.zeros(len(ttk))

    if (not hasattr(extirp_sums, '_n_t_cache') or
        not extirp_sums._n_t_cache == n_t or
        not extirp_sums._delta_cache == delta or
        not np.array_equal(tt_arr, extirp_sums._tt_cache)):

        extirp_sums._n_t_cache = n_t
        extirp_sums._delta_cache = delta
        extirp_sums._tt_cache = copy.deepcopy(tt_arr)

        n_extirp_terms = 25

        dexes = np.zeros(n_extirp_terms, dtype=int)
        dex_range = np.arange(len(dexes),dtype=int)
        half_dexes = len(dexes)//2
        half_dex_range = np.arange(-half_dexes, half_dexes+1, 1, dtype=int)
        assert len(half_dex_range) == n_extirp_terms

        time_dexes = np.round((tt_arr-ttk.min())/delta).astype(int)
        dex_arr = np.array([tj + half_dex_range for tj in time_dexes])

        n_neg = 1
        while n_neg>0:
            negative_dexes = np.where(dex_arr<0)
            n_neg = len(negative_dexes[0])
            dex_arr[negative_dexes[0],:] += 1
        n_pos = 1
        while n_pos>0:
            positive_dexes = np.where(dex_arr>=len(ttk))
            n_pos = len(positive_dexes[0])
            dex_arr[positive_dexes[0],:] -= 1

        extirp_sums.dex_arr = dex_arr

        col_range = np.arange(dex_arr.shape[1], dtype=int)
        col_range_matrix = np.array([np.where(col_range != i_col)[0]
                                     for i_col in range(n_extirp_terms)])

        other_times_list = []
        for i_col in range(dex_arr.shape[1]):
            col_dexes = dex_arr[:,col_range_matrix[i_col]]
            other_times_list.append(np.array([ttk[cc] for cc in col_dexes]).transpose())


        extirp_sums.coeff_cache = []

        for i_col in range(extirp_sums.dex_arr.shape[1]):
            target_dexes = extirp_sums.dex_arr[:,i_col]
            other_times = other_times_list[i_col]
            num = np.product((tt_arr - other_times), axis=0)
            denom = np.product((ttk[target_dexes] - other_times), axis=0)
            extirp_sums.coeff_cache.append(num/denom)

    for i_col in range(extirp_sums.dex_arr.shape[1]):
        target_dexes = extirp_sums.dex_arr[:,i_col]
        term = ff_arr*extirp_sums.coeff_cache[i_col]
        unq_targets, unq_dexes, ct = np.unique(target_dexes, return_counts=True,
                                               return_index=True)
        duplicates = np.where(ct>1)
        for dd, ix in zip(unq_targets[duplicates], unq_dexes[duplicates]):
            sum_dexes = np.where(target_dexes==dd)
            term[ix] = term[sum_dexes].sum()

        hk[unq_targets] += term[unq_dexes]

    ft_re, ft_im = ffter.fft_real(ttk, hk, cache=cache_fft)

    return (ft_re/delta, ft_im/delta)

class LombScargle_PressRybicki(object):

    def __init__(self):
        pass

    def _initialize_PressRybicki(self, time_arr, sigma_arr, delta, ffter, ffter2,
                                 cache_fft=False):
        """
        Initialize arrays for the Press and Rybicki periodogram
        that only depend on t and sigma

        Parameters
        ----------
        time_arr is a numpy array of the times of the measurement

        sigma_arr is a numpy array of the uncertainties of the
        measurement

        ffter and ffter2 are instantiations of FFTransformer.  One is
        for extirpolating sums over omega*t; one is for extirpolating
        sums over 2*omega*t

        cache_fft is a boolean.  If true, use memory intensive caches
        to speed up the FFT.

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
        #print('n_t %d\ndelta %e' % (int(n_t), delta))

        freq_arr = np.array([k/(delta*n_t) for k in range(n_t)])

        w = (1.0/np.power(sigma_arr, 2)).sum()
        wgt_fn = 1.0/(w*np.power(sigma_arr, 2))

        c2_raw, s2_raw = extirp_sums(2.0*time_arr, wgt_fn,
                                    delta, n_t*2, ffter2,
                                    cache_fft=cache_fft)
        dexes = range(0,len(c2_raw), 2)
        c2 = c2_raw[dexes]
        del c2_raw
        s2 = s2_raw[dexes]
        del s2_raw
        gc.collect()

        c, s = extirp_sums(time_arr, wgt_fn, delta, n_t, ffter,
                           cache_fft=cache_fft)

        cut_off_freq = np.exp(-1.3093286772)*np.power(delta, -0.97075831145)
        cut_off_freq *=0.5

        tau = np.arctan2(s2-2*c*s, c2-c*c+s*s)/(4.0*np.pi*freq_arr)

        cos_omega_tau = np.cos(2.0*np.pi*freq_arr*tau)
        sin_omega_tau = np.sin(2.0*np.pi*freq_arr*tau)
        cos_tau = c*cos_omega_tau + s*sin_omega_tau
        sin_tau = s*cos_omega_tau - c*sin_omega_tau

        del c
        del s
        gc.collect()

        cos_2omega_tau = np.cos(4.0*np.pi*freq_arr*tau)
        sin_2omega_tau = np.sin(4.0*np.pi*freq_arr*tau)
        w_sum = wgt_fn.sum()

        csq = 0.5*w_sum + 0.5*c2*cos_2omega_tau + 0.5*s2*sin_2omega_tau
        ssq = 0.5*w_sum - 0.5*c2*cos_2omega_tau - 0.5*s2*sin_2omega_tau
        csomega = 0.5*(s2*cos_2omega_tau - c2*sin_2omega_tau)  # cos(theta)*sin(theta) = 0.5*sin(2*theta)

        del s2
        del c2
        gc.collect()

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

    def get_ls_PressRybicki(self, time_arr_in, f_arr_in, sigma_arr_in, delta,
                            cache_fft=False):
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

        cache_fft is a boolean.  If true, use memory intensive caches
        to speed up the FFT.

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

        if hasattr(self, 'initialized'):
            if not np.array_equal(sigma_arr_in, self.sig_cache):
                self.initialized = False
            if not np.array_equal(time_arr_in, self.time_cache):
                self.initialized = False

        if not hasattr(self, 'ffter'):
            self.ffter = FFTransformer()
            self.ffter2 = FFTransformer()

        if (not hasattr(self, 'initialized') or
            not self.initialized):

            self.initialized = True
            self.sig_cache = copy.deepcopy(sigma_arr_in)
            self.time_cache = copy.deepcopy(time_arr_in)

            (self.w,
             self.delta,
             self.n_t,
             self.freq_arr,
             self.cut_off_freq,
             self.tau,
             self.cos_omega_tau,
             self.sin_omega_tau,
             self.c,
             self.s,
             self.cc,
             self.ss,
             self.cs,
             self.d) = self._initialize_PressRybicki(time_arr, sigma_arr, delta,
                                                     self.ffter, self.ffter2,
                                                     cache_fft=cache_fft)

        y_bar = (f_arr*self.w).sum()
        yy = (self.w*np.power(f_arr-y_bar,2)).sum()
        y_fn = self.w*(f_arr-y_bar)
        y_c_raw, y_s_raw = extirp_sums(time_arr, y_fn,
                                       self.delta,
                                       self.n_t,
                                       self.ffter,
                                       cache_fft=cache_fft)

        y_c = (y_c_raw*self.cos_omega_tau +
               y_s_raw*self.sin_omega_tau)

        y_s = (y_s_raw*self.cos_omega_tau -
               y_c_raw*self.sin_omega_tau)

        del y_s_raw
        del y_c_raw
        gc.collect()

        aa = (y_c*self.ss - y_s*self.cs)/self.d
        bb = (y_s*self.cc - y_c*self.cs)/self.d
        cc = y_bar - aa*self.c - bb*self.s

        pgram = ((y_c*y_c/self.cc) + (y_s*y_s/self.ss))/yy

        return pgram, self.freq_arr, self.tau, aa, bb, cc

    def get_clean_spectrum_PressRybicki(self, time_arr, f_arr,
                                        sigma_arr, delta,
                                        max_components=None,
                                        min_components=None,
                                        cut_off_omega=None,
                                        cache_fft=False):
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

        max_components is an integer indicating the maximum number of
        components to find

        cut_off_omega is an optional maximum allowed angular frequency
        for components

        cache_fft is a boolean.  If true, use memory intensive caches
        to speed up the FFT.

        Returns
        -------
        aa a numpy array of a_i parameters from the model

        bb a numpy array of b_i parameters from the model

        cc a numpy array of c_i parameters from the model

        omega a numpy array of omega_i parameters from the model

        tau a numpy array of tau_i parameters from the model
        """

        residual_arr = copy.deepcopy(f_arr)

        aa_list = []
        bb_list = []
        cc_list = []
        tau_list = []
        omega_list = []
        chisq_arr = []

        model = np.zeros(len(time_arr))
        median_flux = np.median(f_arr)

        residual_arr -= median_flux
        model += median_flux
        chisq = np.power((model-f_arr)/sigma_arr,2).sum()
        ln_n_data = np.log(len(f_arr))
        bic_1 = chisq + ln_n_data
        bic_0 = np.power(f_arr/sigma_arr,2).sum()

        (pspec, freq_arr,
         tau, aa, bb, cc) = self.get_ls_PressRybicki(time_arr, residual_arr,
                                                     sigma_arr, delta,
                                                     cache_fft=cache_fft)

        significance = False
        if bic_1<bic_0:
            significance = True
        bic_0=bic_1

        if cut_off_omega is not None:
            if cut_off_omega/(2.0*np.pi) < self.cut_off_freq:
                cut_off_freq = cut_off_omega/(2.0*np.pi)
            else:
                cut_off_freq = self.cut_off_freq
        else:
            cut_off_freq = self.cut_off_freq

        while significance:

            if max_components is not None and len(aa_list)>=max_components:
                break

            valid = np.where(np.logical_and(np.logical_not(np.isnan(pspec)),
                                            freq_arr<cut_off_freq))
            pspec = pspec[valid]
            valid_freq = freq_arr[valid]
            max_dex = np.argmax(pspec)

            poss = np.where(pspec>=0.75*pspec[max_dex])
            freq_poss = valid_freq[poss]
            best_dex = np.argmin(freq_poss)
            freq_best = freq_poss[best_dex]

            best_dex = poss[0][best_dex]
            best_dex = valid[0][best_dex]

            tau_best = tau[best_dex]
            aa_best = aa[best_dex]
            bb_best = bb[best_dex]
            cc_best = cc[best_dex]
            omega_best = freq_best*2.0*np.pi

            t_arg = omega_best*(time_arr-time_arr.min()-tau_best)
            local_model = np.array([cc_best]*len(time_arr))
            local_model += aa_best*np.cos(t_arg)
            local_model += bb_best*np.sin(t_arg)

            residual_arr -= local_model

            model += local_model

            chisq=np.power((model-f_arr)/sigma_arr,2).sum()
            bic_1 = chisq + (1.0+4.0*len(aa_list))*ln_n_data

            if bic_1>bic_0 and (min_components is None or
                                len(aa_list)>=min_components):
                break

            aa_list.append(aa_best)
            bb_list.append(bb_best)
            cc_list.append(cc_best)
            tau_list.append(tau_best)
            omega_list.append(omega_best)
            chisq_arr.append(chisq)
            #print "%d components %e ; bic %e %e" % (len(aa_list), omega_best, bic_0, bic_1)

            bic_0 = bic_1

            (pspec, freq_arr,
             tau, aa, bb, cc) = self.get_ls_PressRybicki(time_arr,
                                                         residual_arr,
                                                         sigma_arr, delta,
                                                         cache_fft=cache_fft)

        return (median_flux,np.array(aa_list), np.array(bb_list),
                np.array(cc_list), np.array(omega_list),
                np.array(tau_list), np.array(chisq_arr))
