from __future__ import with_statement
import os
import numpy as np
import astropy.io.fits as fits
import tarfile

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--in_dirs', type=str, default=None, nargs='+')
parser.add_argument('--out_dir', type=str, default=None)

args = parser.parse_args()
if args.out_dir is None:
    raise RuntimeError('must specify out_dir')

if args.in_dirs is None:
    raise RuntimeError('must specify in_dirs')

if not isinstance(args.in_dirs, list):
    in_dirs = [args.in_dirs]
else:
    in_dirs = args.in_dirs

out_dir = args.out_dir
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

has_written = []

end_points = {}

for lc_dir in in_dirs:
    tar_file_list = os.listdir(lc_dir)
    for tar_file_name in tar_file_list:
        if (not tar_file_name.endswith('.tgz') or
            not 'long' in tar_file_name):

            continue

        print 'running on ',tar_file_name
        tar_file = tarfile.open(os.path.join(lc_dir, tar_file_name), mode='r')
        member_list = tar_file.getmembers()
        for member in member_list:
            fits_name = member.name
            if not fits_name.endswith('fits'):
                continue

            try:
                obj_name = fits_name.split('/')[-1].split('-')[0]
                out_name = os.path.join(out_dir, '%s_lc.txt' % obj_name)
                file_obj = tar_file.extractfile(member)
                fits_file = fits.open(file_obj)
                mjd = fits_file[1].data['TIME']
                flux = fits_file[1].data['PDCSAP_FLUX']
                sig = fits_file[1].data['PDCSAP_FLUX_ERR']
                valid_dexes = np.where(np.logical_and(np.logical_not(np.isnan(mjd)),
                                       np.logical_and(np.logical_not(np.isnan(flux)),
                                                      np.logical_not(np.isnan(sig)))))

                mjd = mjd[valid_dexes]
                flux = flux[valid_dexes]
                sig = sig[valid_dexes]

                if out_name not in end_points:
                    end_points[out_name] = []
                end_points[out_name].append((mjd[0], mjd[-1]))

                if out_name in has_written:
                    mode = 'a'
                else:
                    mode = 'w'
                    has_written.append(out_name)

                with open(out_name, mode) as output_file:
                    for tt, ff, ss in zip(mjd,flux,sig):
                        output_file.write('%.12e %e %e\n' % (tt, ff, ss))
            except:
                print 'failed on: %s' % fits_name
                raise


        print 'processed in %s' % tar_file_name

for out_name in end_points:
    with open(out_name, 'a') as output_file:
        output_file.write('# observing quarter end points\n')
        for span in end_points[out_name]:
            output_file.write('# %.12e %.12e\n' % (span[0], span[1]))
