from __future__ import with_statement
import argparse
import time

from multiprocessing import Process

from clean_by_segments import clean_spectra

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

    parser.add_argument('--n_p', type=int, default=32,
                        help='number of processes to spawn')

    parser.add_argument('--sleep', type=int, default=60,
                        help='time to wait between processes')

    args = parser.parse_args()

    lc_lists = []
    for ix in range(args.n_p):
        lc_lists.append([])

    i_list = 0
    with open(args.list, 'r') as input_file:
        for line in input_file:
            print 'reading ',line
            lc_lists[i_list].append(line.strip())
            i_list += 1
            if i_list >= len(lc_lists):
                i_list = 0

    if args.do_stitch.lower()[0] == 't':
        do_stitch = True
    else:
        do_stitch = False

    if args.cache_fft.lower()[0] == 't':
        cache_fft = True
    else:
        cache_fft = False

    print 'lc_lists ',lc_lists

    t_start = time.time()
    for ix in range(len(lc_lists)):
        p = Process(target=clean_spectra,
                    args=(lc_lists[ix],
                          '%s_%d.txt' % (args.out_file, ix)),
                    kwargs={'in_dir': args.in_dir,
                            'do_stitch': do_stitch,
                            'stitch_dir': args.stitch_dir,
                            'max_components': args.max_components,
                            'log_file': '%s_%d.txt' % (args.log_file, ix),
                            'dt_factor': args.dt,
                            'write_every': args.write_every,
                            'fig_dir': args.fig_dir,
                            'cache_fft': cache_fft})

        p.start()
        time.sleep(args.sleep)
        print 'started %d at %e ' % (ix, time.time()-t_start)
