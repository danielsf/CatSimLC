from __future__ import with_statement
import os
import numpy as np
import argparse
import time

from PressRybicki import LombScargle_PressRybicki

__all__ = ['clean_spectra']

def _fit_and_offset(PRobj,
                    time_to_fit, flux_to_fit, sigma_to_fit,
                    time_to_offset, flux_to_offset, sigma_to_offset,
                    cache_fft=False, dt_factor=0.1):

        dt = dt_factor*np.diff(np.unique(time_to_fit)).min()

        (median_flux, aa, bb, cc,
         omega, tau, chisq_arr) = PRobj.get_clean_spectrum_PressRybicki(time_to_fit,
                                                                        flux_to_fit,
                                                                        sigma_to_fit, dt,
                                                                        min_components=3,
                                                                        max_components=3,
                                                                        cache_fft=cache_fft,
                                                                        cut_off_omega=200.0)

        model = np.array([median_flux]*len(time_to_offset))
        for ix in range(len(aa)):
            model += cc[ix]
            t_arg = omega[ix]*(time_to_offset - time_to_fit.min() - tau[ix])
            model += aa[ix]*np.cos(t_arg)
            model += bb[ix]*np.sin(t_arg)

        offset_num = ((flux_to_offset-model)/np.power(sigma_to_offset,2)).sum()
        offset_denom = (1.0/np.power(sigma_to_offset,2)).sum()
        offset = offset_num/offset_denom
        chisq = np.power((model-flux_to_offset-offset)/sigma_to_offset,2).sum()

        return offset, chisq


def re_calibrate_lc(PRobj, time_arr, flux_arr, sigma_arr, segments, cache_fft=False,
                    dt_factor=0.1):
    """
    Stitch together the differently calibrated segments of a light curve.

    Parameters
    ----------
    PRobj -- an instantiation of LombScargle_PressRybicki to do the work

    time_arr -- a numpy array containing the time axis of the light curve

    flux_arr -- a numpy array containing the flux of the light curve

    sigma_arr -- a numpy array containing the uncertainty of the flux

    segments -- a list of tuples; each tuple contains the start time and end
    time of each light curve segment.

    Return
    ------
    stitch_time -- a numpy array of the stitched-together light curve time axis

    stitch_flux -- a numpy array of the stitched-together light curve flux axis

    stitch_sigma -- a numpy array of the stitched-together light curve
    flux uncertainty
    """

    sorted_dex = np.argsort(time_arr)
    time = time_arr[sorted_dex]
    flux = flux_arr[sorted_dex]
    sigma = sigma_arr[sorted_dex]

    start_times = []
    end_times = []
    for seg in segments:
        start_times.append(seg[0])
        end_times.append(seg[1])
    start_times = np.array(start_times)
    end_times = np.array(end_times)

    sorted_dex = np.argsort(start_times)
    start_times = start_times[sorted_dex]
    end_times = end_times[sorted_dex]

    tol = 1.0e-6

    # In some cases, the first Kepler observing quarter contains two
    # light curves (with identical time axes) for the same source.
    # I did not know this when I started processing the data.  Rather than
    # separate these light curves from each other and try to pick one, I am
    # just going to discard any light curve segments that have double reporting.
    first_segment = 0
    if len(start_times)>1:
        is_valid = False
        while not is_valid:
            is_valid = True
            if first_segment < len(start_times)-1:
                if (start_times[first_segment] > start_times[first_segment+1]-tol and
                    start_times[first_segment] < end_times[first_segment+1]+tol):

                    is_valid = False

                if (end_times[first_segment] > start_times[first_segment+1]-tol and
                    end_times[first_segment] < end_times[first_segment+1]+tol):

                    is_valid = False

            if first_segment>0:
               if (start_times[first_segment] > start_times[first_segment-1]-tol and
                   start_times[first_segment] < end_times[first_segment-1]-tol):

                   is_valid = False

               if (end_times[first_segment] > start_times[first_segment-1]-tol and
                   end_times[first_segment] < end_times[first_segment-1]-tol):

                   is_valid = False

            if not is_valid:
                first_segment += 1

    #print '    first_segment ',first_segment

    first_dexes = np.where(np.logical_and(time>=start_times[first_segment]-tol,
                                          time<=end_times[first_segment]+tol))

    time_out = np.zeros(len(time), dtype=float)
    flux_out = np.zeros(len(time), dtype=float)
    sigma_out = np.zeros(len(time), dtype=float)

    time_buffer = np.zeros(len(time), dtype=float)
    flux_buffer = np.zeros(len(time), dtype=float)
    sigma_buffer = np.zeros(len(time), dtype=float)

    n_out = len(first_dexes[0])

    time_out[:n_out] = time[first_dexes]
    flux_out[:n_out] = flux[first_dexes]
    sigma_out[:n_out] = sigma[first_dexes]

    for i_seg in range(first_segment+1, len(segments)):

        # again: discard light curve segments with double reporting
        use_segment = True
        if i_seg < len(segments)-1:
            if ((start_times[i_seg] > start_times[i_seg+1]-tol and
                 start_times[i_seg] < end_times[i_seg+1]+tol) or
                (end_times[i_seg] > start_times[i_seg+1]-tol and
                 end_times[i_seg] < end_times[i_seg+1]+tol)):

                use_segment = False

        if ((start_times[i_seg] > start_times[i_seg-1]-tol and
             start_times[i_seg] < end_times[i_seg-1]+tol) or
            (end_times[i_seg] > start_times[i_seg-1]-tol and
             end_times[i_seg] < end_times[i_seg-1]+tol)):

            use_segment = False

        if not use_segment:
            continue


        valid_dex = np.where(np.logical_and(time>=start_times[i_seg]-tol,
                                            time<=end_times[i_seg]+tol))

        next_time = time[valid_dex]
        next_flux = flux[valid_dex]
        next_sigma = sigma[valid_dex]

        if n_out < len(next_time):
            time_to_fit_master = next_time
            flux_to_fit_master = next_flux
            sigma_to_fit_master = next_sigma
            time_to_offset = time_out[:n_out]
            flux_to_offset = flux_out[:n_out]
            sigma_to_offset = sigma_out[:n_out]
        else:
            time_to_fit_master = time_out[:n_out]
            flux_to_fit_master = flux_out[:n_out]
            sigma_to_fit_master = sigma_out[:n_out]
            time_to_offset = next_time
            flux_to_offset = next_flux
            sigma_to_offset = next_sigma

        if 2*len(time_to_offset) < len(time_to_fit_master):
            n_to_fit = 2*len(time_to_offset)
        else:
            n_to_fit = len(time_to_fit_master)

        #print 'segment %e %e %s %e %e' % \
        #(next_time.min(),next_time.max(),(n_to_fit==len(time_to_fit_master)),
        # time_to_fit_master.min(), time_to_fit_master.max())

        time_to_fit = time_to_fit_master[-n_to_fit:]
        flux_to_fit = flux_to_fit_master[-n_to_fit:]
        sigma_to_fit = sigma_to_fit_master[-n_to_fit:]

        offset, chisq = _fit_and_offset(PRobj,
                                        time_to_fit, flux_to_fit, sigma_to_fit,
                                        time_to_offset, flux_to_offset, sigma_to_offset,
                                        cache_fft=cache_fft,
                                        dt_factor=dt_factor)

        med_fit = np.median(flux_to_fit)
        med_offset = np.median(flux_to_offset-offset)
        stdev_fit = np.sqrt(np.power(flux_to_fit-med_fit,2).sum()/(len(flux_to_fit)+1))
        stdev_offset = np.sqrt(np.power(flux_to_offset-offset-med_offset,2).sum()/(len(flux_to_offset)+1))

        if (np.abs(med_fit-med_offset) > min(stdev_fit, stdev_offset) and
            n_to_fit<len(time_to_fit_master)):
            time_to_fit = time_to_fit_master
            flux_to_fit = flux_to_fit_master
            sigma_to_fit = sigma_to_fit_master

            offset, chisq = _fit_and_offset(PRobj,
                                            time_to_fit,
                                            flux_to_fit,
                                            sigma_to_fit,
                                            time_to_offset,
                                            flux_to_offset,
                                            sigma_to_offset,
                                            cache_fft=cache_fft,
                                            dt_factor=dt_factor)

        if time_to_offset[-1] < time_to_fit[0]:
            n_first = len(time_to_offset)
            n_out = n_first + len(time_to_fit_master)
            time_buffer[:n_first] = time_to_offset
            flux_buffer[:n_first] = flux_to_offset-offset
            sigma_buffer[:n_first] = sigma_to_offset

            time_buffer[n_first:n_out] = time_to_fit_master
            flux_buffer[n_first:n_out] = flux_to_fit_master
            sigma_buffer[n_first:n_out] = sigma_to_fit_master

        else:
            n_first = len(time_to_fit_master)
            n_out = n_first + len(time_to_offset)

            time_buffer[:n_first] = time_to_fit_master
            flux_buffer[:n_first] = flux_to_fit_master
            sigma_buffer[:n_first] = sigma_to_fit_master

            time_buffer[n_first:n_out] = time_to_offset
            flux_buffer[n_first:n_out] = flux_to_offset-offset
            sigma_buffer[n_first:n_out] = sigma_to_offset

        time_out[:n_out] = time_buffer[:n_out]
        flux_out[:n_out] = flux_buffer[:n_out]
        sigma_out[:n_out] = sigma_buffer[:n_out]

        #print '    calculated offset %e %d of %d %e -- dt %e %d -- %d' % \
        #(offset,i_seg,len(segments),chisq,dt,len(time_out), len(aa))

    return (time_out[:n_out], flux_out[:n_out], sigma_out[:n_out])


def clean_spectra(list_of_lc, out_file_name, in_dir=None,
                  do_stitch=False, stitch_dir=None, fig_dir=None,
                  log_file=None, write_every=100,
                  max_components=51, dt_factor=0.1, cache_fft=False):
    """
    Clean a list of time series

    Parameters
    ----------
    list_of_lc is a list of filenames, each file storing a light curve

    out_file_name is the name of the file where the cleaning paramters
    will be saved

    in_dir is the directory where the light curve files live

    log_file is the name of a file where we will write logging
    information

    write_every is an integer; write output this offten

    max_components is the maximum number of Fourier components
    to keep for each light curve

    dt_factor controls the size of the time step fed to the FFT

    cache_fft is a boolean controlling whether or not to use
    memory intensive caching in the FFT

    do_stitch is a boolean controlling whether or not we need to
    stitch together the light curves

    stitch_dir is the directory where to which we will write the
    stitched light curves

    fig_dir is the directory to which we will write plots of our
    light curve data and smoothings
    """

    PRobj = LombScargle_PressRybicki()

    if fig_dir is not None:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)


    if os.path.exists(log_file):
        os.unlink(log_file)

    if do_stitch:
        if stitch_dir is None:
            raise RuntimeError("If DO_STITCH is True, must specify STITCH_DIR")

    if do_stitch:
        if not os.path.isdir(stitch_dir):
            os.mkdir(stitch_dir)

    dtype = np.dtype([('t', float), ('f', float), ('s', float)])

    with open(out_file_name, 'w') as out_file:
        out_file.write('# lc_name n_t_steps t_span n_components chisquared median_flux ')
        out_file.write('aa bb cc omega tau ')
        out_file.write('{f = cc + aa*cos(omega*(t-tmin-tau)) + bb*sin(omega*(t-tmin-tau))}\n')

    ct = 0
    t_start = time.time()

    data_dict = {}
    segment_dict = {}
    output_dict = {}
    stitch_dict = {}

    for lc_name_global in list_of_lc:

        if in_dir is not None:
            full_name = os.path.join(in_dir, lc_name_global)
        else:
            full_name = os.path.join(lc_name_global)

        data = np.genfromtxt(full_name, dtype=dtype)

        segments = []
        with open(full_name, 'r') as in_file:
            for line in in_file:
                if line[0] != '#':
                    break
                line = line.strip().split()
                segments.append((float(line[1]), float(line[2])))

        data_dict[lc_name_global] = data
        segment_dict[lc_name_global] = segments
        #print 'read in ',len(data_dict),time.time()-t_start

        if len(data_dict) >= write_every or lc_name_global==list_of_lc[-1]:
            for lc_name in data_dict:
                data = data_dict[lc_name]
                segments = segment_dict[lc_name]

                if do_stitch:
                    #print data.shape
                    time_arr, flux_arr, sigma_arr = re_calibrate_lc(PRobj,
                                                                    data['t'], data['f'],
                                                                    data['s'], segments,
                                                                    cache_fft=cache_fft,
                                                                    dt_factor=dt_factor)

                    stitch_name = lc_name.replace('.txt','')
                    stitch_name = os.path.join(stitch_dir, stitch_name+'_stitched.tx')
                    stitch_dict[lc_name] = (time_arr, flux_arr, sigma_arr)
                else:
                    time_arr = data['t']
                    flux_arr = data['f']
                    sigma_arr = data['s']

                try:
                    assert len(time_arr) == len(np.unique(time_arr))
                except:
                    print 'unique(time) failed on ',lc_name
                    continue

                dt = dt_factor*np.diff(np.unique(time_arr)).min()

                (median_flux,
                 aa, bb,
                 cc, omega,
                 tau, chisq_arr) = PRobj.get_clean_spectrum_PressRybicki(time_arr, flux_arr,
                                                                         sigma_arr, dt,
                                                                         max_components=max_components,
                                                                         cut_off_omega=200.0,
                                                                         cache_fft=cache_fft)

                output_dict[lc_name] = {}
                output_dict[lc_name]['span'] = time_arr.max() - time_arr.min()
                output_dict[lc_name]['tsteps'] = len(time_arr)
                output_dict[lc_name]['chisq'] = chisq_arr
                output_dict[lc_name]['median'] = median_flux
                output_dict[lc_name]['aa'] = aa
                output_dict[lc_name]['bb'] = bb
                output_dict[lc_name]['cc'] = cc
                output_dict[lc_name]['omega'] = omega
                output_dict[lc_name]['tau'] = tau

                ct += 1

                #print 'done with %d in %e' % (ct, time.time()-t_start)
                if ct%2 == 0 or ct==len(list_of_lc):
                    with open(log_file, 'a') as out_file:
                        out_file.write('finished %d in %e sec; should take %e days\n' %\
                        (ct, time.time()-t_start, (len(list_of_lc)*(time.time()-t_start)/ct)/(24.0*3600.0)))

            data_dict = {}
            segment_dict = {}

        if len(output_dict) >= write_every or lc_name_global == list_of_lc[-1]:
            with open(out_file_name, 'a') as out_file:
                for lc_name in output_dict:
                    out_file.write('%s %d %e %d ' % (lc_name,
                                                     output_dict[lc_name]['tsteps'],
                                                     output_dict[lc_name]['span'],
                                                     len(output_dict[lc_name]['aa'])))

                    assert len(output_dict[lc_name]['chisq']) == len(output_dict[lc_name]['aa'])

                    for ix in range(len(output_dict[lc_name]['chisq'])):
                        out_file.write('%e ' % output_dict[lc_name]['chisq'][ix])

                    out_file.write('%e ' % output_dict[lc_name]['median'])

                    for ix in range(len(output_dict[lc_name]['aa'])):
                         out_file.write('%e %e %e %e %e ' % (output_dict[lc_name]['aa'][ix],
                                                             output_dict[lc_name]['bb'][ix],
                                                             output_dict[lc_name]['cc'][ix],
                                                             output_dict[lc_name]['omega'][ix],
                                                             output_dict[lc_name]['tau'][ix]))
                    out_file.write('\n')

            if do_stitch:
                for lc_name in stitch_dict:
                    out_name = os.path.join(stitch_dir,
                                            lc_name.replace('.txt','')+'_stitched.txt')
                    with open(out_name, 'w') as out_file:
                        for tt, ff, ss in zip(stitch_dict[lc_name][0],
                                              stitch_dict[lc_name][1],
                                              stitch_dict[lc_name][2]):

                            out_file.write('%.12e %e %e\n' % (tt, ff, ss))

                if fig_dir is not None:
                    for lc_name in output_dict:
                        tt = stitch_dict[lc_name][0]
                        ff = stitch_dict[lc_name][1]
                        f_sorted = np.sort(ff)
                        y_min = f_sorted[len(ff)/10]
                        y_max = f_sorted[len(ff)*9/10]
                        dt = 0.1*np.min(np.diff(np.unique(tt)))
                        model_t = np.arange(tt.min(), tt.max(), dt)
                        model = np.zeros(len(model_t))
                        model += output_dict[lc_name]['median']
                        aa = output_dict[lc_name]['aa']
                        bb = output_dict[lc_name]['bb']
                        cc = output_dict[lc_name]['cc']
                        omega = output_dict[lc_name]['omega']
                        tau = output_dict[lc_name]['tau']
                        for ix in range(len(aa)):
                            model += cc[ix]
                            t_arg = omega[ix]*(model_t-tt.min()-tau[ix])
                            model += aa[ix]*np.cos(t_arg)
                            model += bb[ix]*np.sin(t_arg)
                        with open(os.path.join(fig_dir, lc_name.replace('.txt', '_model.txt')), 'w') as out_file:
                            for t_val, f_val in zip(model_t, model):
                                out_file.write('%e %e\n' % (t_val, f_val))

                        plt.figsize = (30,30)
                        plt.subplot(3,1,1)
                        t_dex = np.where(tt<tt.min()+(tt.max()-tt.min())/3.0)
                        plt.scatter(tt[t_dex], ff[t_dex], s=5, color='k', edgecolor='',
                                    zorder=1)

                        plt.ylim((y_min,y_max))
                        plt.xticks(fontsize=10)
                        plt.yticks(fontsize=10)
                        t_dex = np.where(model_t<tt.min()+(tt.max()-tt.min())/3.0)
                        plt.plot(model_t[t_dex], model[t_dex], color='r', linewidth=1,
                                 zorder=2)

                        plt.subplot(3,1,2)
                        t_dex = np.where(np.logical_and(tt<tt.min()+2.0*(tt.max()-tt.min())/3.0,
                                                        tt>tt.min()+(tt.max()-tt.min())/4.0))
                        plt.scatter(tt[t_dex], ff[t_dex], s=5, color='k', edgecolor='',
                                    zorder=1)

                        plt.ylim((y_min,y_max))
                        plt.xticks(fontsize=10)
                        plt.yticks(fontsize=10)
                        t_dex = np.where(np.logical_and(model_t<tt.min()+2.0*(tt.max()-tt.min())/3.0,
                                                        model_t>tt.min()+(tt.max()-tt.min())/4.0))
                        plt.plot(model_t[t_dex], model[t_dex], color='r', linewidth=1,
                                 zorder=2)

                        plt.subplot(3,1,3)
                        t_dex = np.where(tt>tt.min()+2.0*(tt.max()-tt.min())/3.0-0.1*(tt.max()-tt.min()))
                        plt.scatter(tt[t_dex], ff[t_dex], s=5, color='k', edgecolor='',
                                    zorder=1)

                        plt.ylim((y_min,y_max))
                        plt.xticks(fontsize=10)
                        plt.yticks(fontsize=10)
                        t_dex = np.where(model_t>tt.min()+2.0*(tt.max()-tt.min())/3.0-0.1*(tt.max()-tt.min()))
                        plt.plot(model_t[t_dex], model[t_dex], color='r', linewidth=1,
                                 zorder=2)

                        plt.tight_layout()
                        plt.savefig(os.path.join(fig_dir, lc_name.replace('.txt','.png')))
                        plt.close()

            output_dict = {}
            stitch_dict = {}

    #print 'that took ',time.time()-t_start


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--list', type=str, default=None,
                        help='text file containing a list of the names of the '
                             'light curves to be processed')

    parser.add_argument('--in_dir', type=str, default=None,
                        help='directory where the light curve files referenced in '
                             'LIST reside')

    parser.add_argument('--do_stitch', type=str, default='True',
                        help = 'do we need to stitch these light curves together '
                               '(i.e. do we need to deal with the variable Kepler '
                               'calibration).  True/False.')

    parser.add_argument('--stitch_dir', type=str, default=None,
                        help='if DO_STITCH is True, in what directory should we '
                             'save the stitched together light curves?')

    parser.add_argument('--out_file', type=str, default=None,
                        help='in what file should we save the parametrization of '
                             'the light curves')

    parser.add_argument('--max_components', type=int, default=51,
                        help='maximum number of components to use when smoothing '
                             'light curves')

    parser.add_argument('--log_file', type=str, default='lc_timing_log.txt',
                        help='log file where timing information is written')

    parser.add_argument('--dt', type=float, default=0.1,
                        help='what fraction of minimum data dt do we use as '
                             'dt input to Lomb-Scargle periodogram')

    parser.add_argument('--write_every', type=int, default=100,
                        help='how often do we write to output '
                             '(in units of completed light curves)')

    parser.add_argument('--fig_dir', type=str, default=None,
                        help='directory in which to output plots')

    parser.add_argument('--cache_fft', type=str, default='False',
                        help='should use the memory-intensive index cache '
                             'when FFTing the data')

    args = parser.parse_args()

    if args.cache_fft.lower()[0] == 'f':
        cache_fft = False
    else:
        cache_fft = True

    if args.do_stitch.lower()[0] == 't':
        do_stitch = True
    else:
        do_stitch = False

    if args.list is None:
        raise RuntimeError("Must specify LIST")

    if args.out_file is None:
        raise RuntimeError("Must specify OUT_FILE")


    list_of_lc = []
    with open(args.list, 'r') as in_file:
        for line in in_file:
            list_of_lc.append(line.strip())

    #print 'list ',list_of_lc

    clean_spectra(list_of_lc, args.out_file, cache_fft=cache_fft,
                  do_stitch=do_stitch, stitch_dir=args.stitch_dir,
                  log_file=args.log_file, in_dir=args.in_dir,
                  dt_factor=args.dt, max_components=args.max_components,
                  write_every=args.write_every, fig_dir=args.fig_dir)
