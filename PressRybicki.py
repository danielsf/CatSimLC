import numpy as np
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
