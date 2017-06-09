from __future__ import with_statement
import os
import numpy as np
import argparse
import time

from PressRybicki import get_clean_spectrum_PressRybicki

def re_calibrate_lc(time_arr, flux_arr, sigma_arr, segments):
    """
    Stitch together the differently calibrated segments of a light curve.

    Parameters
    ----------
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

    time_out = time[first_dexes]
    flux_out = flux[first_dexes]
    sigma_out = sigma[first_dexes]

    seg_ct = 1

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

        if len(time_out) < len(next_time):
            time_to_fit_master = next_time
            flux_to_fit_master = next_flux
            sigma_to_fit_master = next_sigma
            time_to_offset = time_out
            flux_to_offset = flux_out
            sigma_to_offset = sigma_out
        else:
            time_to_fit_master = time_out
            flux_to_fit_master = flux_out
            sigma_to_fit_master = sigma_out
            time_to_offset = next_time
            flux_to_offset = next_flux
            sigma_to_offset = next_sigma

        if 2*len(time_to_offset) < len(time_to_fit_master):
            n_to_fit = 2*len(time_to_offset)
        else:
            n_to_fit = len(time_to_fit_master)

        seg_ct += 1

        time_to_fit = time_to_fit_master[-n_to_fit:]
        flux_to_fit = flux_to_fit_master[-n_to_fit:]
        sigma_to_fit = sigma_to_fit_master[-n_to_fit:]

        dt = 0.1*np.diff(np.unique(time_to_fit)).min()

        (median_flux, aa, bb, cc,
         omega, tau, chisq_arr) = get_clean_spectrum_PressRybicki(time_to_fit,
                                                                  flux_to_fit,
                                                                  sigma_to_fit, dt,
                                                                  min_components=3,
                                                                  max_components=3)

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

        if time_to_offset[-1] < time_to_fit[0]:
            time_out = np.append(time_to_offset, time_to_fit_master)
            flux_out = np.append(flux_to_offset-offset, flux_to_fit_master)
            sigma_out = np.append(sigma_to_offset, sigma_to_fit_master)
        else:
            time_out = np.append(time_to_fit_master, time_to_offset)
            flux_out = np.append(flux_to_fit_master, flux_to_offset-offset)
            sigma_out = np.append(sigma_to_fit_master, sigma_to_offset)

        #print '    calculated offset %e %d of %d %e -- dt %e %d -- %d' % \
        #(offset,i_seg,len(segments),chisq,dt,len(time_out), len(aa))

    return (time_out, flux_out, sigma_out)


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

args = parser.parse_args()

if args.do_stitch.lower()[0] == 't':
    do_stitch = True
else:
    do_stitch = False

if args.list is None:
    raise RuntimeError("Must specify LIST")

if args.out_file is None:
    raise RuntimeError("Must specify OUT_FILE")

if do_stitch:
    if args.stitch_dir is None:
        raise RuntimeError("If DO_STITCH is True, must specify STITCH_DIR")

if do_stitch:
    if not os.path.isdir(args.stitch_dir):
        os.mkdir(args.stitch_dir)

list_of_lc = []
with open(args.list, 'r') as in_file:
    for line in in_file:
        list_of_lc.append(line.strip())

dtype = np.dtype([('t', float), ('f', float), ('s', float)])

write_every = 500

with open(args.out_file, 'w') as out_file:
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

    full_name = os.path.join(args.in_dir, lc_name_global)
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
    print 'read in ',len(data_dict),time.time()-t_start

    if len(data_dict) >= write_every:
        for lc_name in data_dict:
            data = data_dict[lc_name]
            segments = segment_dict[lc_name]

            if do_stitch:
                time_arr, flux_arr, sigma_arr = re_calibrate_lc(data['t'], data['f'],
                                                                data['s'], segments)

                stitch_name = lc_name.replace('.txt','')
                stitch_name = os.path.join(args.stitch_dir, stitch_name+'_stitched.tx')
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

            dt = 0.1*np.diff(np.unique(time_arr)).min()

            (median_flux,
             aa, bb,
             cc, omega,
             tau, chisq_arr) = get_clean_spectrum_PressRybicki(time_arr, flux_arr,
                                                               sigma_arr, dt,
                                                               max_components=args.max_components)

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
            if ct%10 == 0:
                with open(args.log_file, 'a') as out_file:
                    out_file.write('finished %d in %e sec; should take %e days\n' %\
                    (ct, time.time()-t_start, (len(list_of_lc)*(time.time()-t_start)/ct)/(24.0*3600.0)))

        data_dict = {}
        segment_dict = {}

    if len(output_dict) >= write_every or lc_name == list_of_lc[-1]:
        with open(args.out_file, 'a') as out_file:
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
                out_name = os.path.join(args.stitch_dir,
                                        lc_name.replace('.txt','')+'_stitched.txt')
                with open(out_name, 'w') as out_file:
                    for tt, ff, ss in zip(stitch_dict[lc_name][0],
                                          stitch_dict[lc_name][1],
                                          stitch_dict[lc_name][2]):

                        out_file.write('%.12e %e %e\n' % (tt, ff, ss))

        output_dict = {}
        stitch_dict = {}

print 'that took ',time.time()-t_start
