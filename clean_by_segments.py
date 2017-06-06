from __future__ import with_statement
import os
import numpy as np
import argparse

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
    first_dexes = np.where(np.logical_and(time>=start_times[0]-tol,
                                          time<=end_times[0]+tol))

    time_out = time[first_dexes]
    flux_out = flux[first_dexes]
    sigma_out = sigma[first_dexes]

    for i_seg in range(1, len(segments)):

        dt = 0.1*np.diff(np.unique(time_out)).min()

        (median_flux, aa, bb, cc,
         omega, tau, freq) = get_clean_spectrum_PressRybicki(time_out, flux_out,
                                                             sigma_out, dt,
                                                             min_components=3,
                                                             max_components=3)

        valid_dex = np.where(np.logical_and(time>=start_times[i_seg]-tol,
                                            time<=end_times[i_seg]+tol))

        local_time = time[valid_dex]
        local_flux = flux[valid_dex]
        local_sigma = sigma[valid_dex]

        model = np.array([median_flux]*len(local_time))
        for ix in range(len(aa)):
            model += cc[ix]
            t_arg = omega[ix]*(local_time - local_time.min() - tau[ix])
            model += aa[ix]*np.cos(t_arg)
            model += bb[ix]*np.sin(t_arg)

        offset_num = ((local_flux-model)/np.power(local_sigma,2)).sum()
        offset_denom = (1.0/np.power(local_sigma,2)).sum()
        offset = offset_num/offset_denom

        time_out = np.append(time_out, local_time)
        flux_out = np.append(flux_out, local_flux-offset)
        sigma_out = np.append(sigma_out, local_sigma)

    return (time_out, flux_out, sigma_out)


parser = argparse.ArgumentParser()

parser.add_argument('--list', type=str, default=None,
                    help='text file containing a list of the names of the '
                         'light curves to be processed')

parser.add_argument('--in_dir', type=str, default=None,
                    help='directory where the light curve files referenced in '
                         'LIST reside')

parser.add_argument('--do_stitch', type=bool, default=True,
                    help = 'do we need to stitch these light curves together '
                           '(i.e. do we need to deal with the variable Kepler '
                           'calibration)')

parser.add_argument('--stitch_dir', type=str, default=None,
                    help='if DO_STITCH is True, in what directory should we '
                         'save the stitched together light curves?')

parser.add_argument('--out_file', type=str, default=None,
                    help='in what file should we save the parametrization of '
                         'the light curves')

parser.add_argument('--max_components', type=int, default=51,
                    help='maximum number of components to use when smoothing '
                         'light curves')

args = parser.parse_args()

if args.list is None:
    raise RuntimeError("Must specify LIST")

if args.out_file is None:
    raise RuntimeError("Must specify OUT_FILE")

if args.do_stitch:
    if args.stitch_dir is None:
        raise RuntimeError("If DO_STITCH is True, must specify STITCH_DIR")

list_of_lc = []
with open(args.list, 'r') as in_file):
    for line in in_file:
        list_of_lc.append(line.strip())

dtype = np.dtype([('t', float), ('f', float), ('s', float)])

write_every = 500

output_dict = {}
stitch_dict = {}

with open(args.out_file, 'w') as out_file:
    out_file.write('# lc_name n_t_steps t_span chisquared n_components median_flux '
    out_file.write('aa bb cc omega tau '
    out_file.write('{f = cc + aa*cos(omega*(t-tmin-tau)) + bb*sin(omega*(t-tmin-tau))}\n')

for lc_name in list_of_lc:
    full_name = os.path.join(args.in_dir, lc_name)
    data = np.genfromtxt(full_name, dtype=dtype)

    segments = []
    with open(full_name, 'r') as in_file:
        for line in in_file:
            if line[0] != '#':
                break
            line = line.strip().split()
            segments.append((float(line[1]), float(line[2])))

    if args.do_stitch:
        time_arr, flux_arr, sigma_arr = re_calibrate_lc(data['t'], data['f'],
                                                        data['s'], segments)

        stitch_name = lc_name.replace('.txt','')
        stitch_name = os.path.join(args.stitch_dir, stitch_name+'_stitched.tx')
        stitch_dict[lc_name] = (time_arr, flux_arr, sigma_arr)
    else:
        time_arr = data['t']
        flux_arr = data['f']
        sigma_arr = data['s']

    dt = 0.1*np.diff(np.unique(time_arr)).min()

    (median_flux,
     aa, bb,
     cc, omega,
     tau, freq) = get_clean_spectrum_PressRybicki(time_arr, flux_arr,
                                                  sigma_arr, dt,
                                                  max_components=args.max_components)

    model = np.array([median_flux]*len(time_arr))
    for ix in range(len(aa)):
        model += cc[ix]
        t_arg = omega[ix]*(time_arr-time_arr.min()-tau[ix])
        model += aa[ix]*np.cos(t_arg)
        model += bb[ix]*np.sin(t_arg)

    chisq = np.power((model-flux_arr)/sigma_arr,2).sum()
    output_dict[lc_name] = {}
    output_dict[lc_name]['span'] = time_arr.max() - time_arr.min()
    output_dict[lc_name]['tsteps'] = len(time_arr)
    output_dict[lc_name]['chisq'] = chisq
    output_dict[lc_name]['median'] = median_flux
    output_dict[lc_name]['aa'] = aa
    output_dict[lc_name]['bb'] = bb
    output_dict[lc_name]['cc'] = cc
    output_dict[lc_name]['omega'] = omega
    output_dict[lc_name]['tau'] = tau

    if len(output_dict) >= write_every or lc_name == list_of_lc[-1]:
        with open(args.out_file, 'w') as out_file:
            out_file.write('%s %d %e %e %d %e' % (lc_name,
                                                  output_dict[lc_name]['tsteps'],
                                                  output_dict[lc_name]['span'],
                                                  output_dict[lc_name]['chisq'],
                                                  len(output_dict[lc_name]['aa']),
                                                  output_dict[lc_name]['median']))

            for ix in range(len(output_dict[lc_name]['aa'])):
                 out_file.write('%e %e %e %e %e\n' % (output_dict[lc_name]['aa'][ix],
                                                      output_dict[lc_name]['bb'][ix],
                                                      output_dict[lc_name]['cc'][ix],
                                                      output_dict[lc_name]['omega'][ix],
                                                      output_dict[lc_name]['tau'][ix]))
