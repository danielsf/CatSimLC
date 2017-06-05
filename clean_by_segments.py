from __future__ import with_statement
import os
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--objid', type=str, default=None)
parser.add_argument('--outfile', type=str, default=None)

args = parser.parse_args()
if args.objid is None:
    raise RuntimeError('Must specify objid')
if args.outfile is None:
    raise RuntimeError('Must specify outfile')

data_dir = os.path.join('workspace', 'data')
list_of_files = os.listdir(data_dir)
list_of_files.sort()

import astropy.io.fits as fits

segment_list = []
for file_name in list_of_files:
    if args.objid in file_name:
        data = fits.open(os.path.join(data_dir,file_name))
        d = data[1].data
        valid_dexes = np.where(np.logical_and(np.logical_not(np.isnan(d['TIME'])),
                               np.logical_and(np.logical_not(np.isnan(d['PDCSAP_FLUX'])),
                                              np.logical_not(np.isnan(d['PDCSAP_FLUX_ERR'])))))
        segment_list.append((data[1].data['TIME'][valid_dexes],
                             data[1].data['PDCSAP_FLUX'][valid_dexes],
                             data[1].data['PDCSAP_FLUX_ERR'][valid_dexes]))

time_arr = None
flux_arr = None
sig_arr = None

from PressRybicki import get_clean_spectrum_PressRybicki
n_steps = 1048576
#n_steps= 1024
components = 3

import time

t_start = time.time()
for i_seg, segment in enumerate(segment_list):
    if time_arr is None:
        time_arr = segment[0]
        flux_arr = segment[1]
        sig_arr = segment[2]
        continue

    #dt = (time_arr.max()-time_arr.min())/n_steps
    dt = 0.1*np.diff(np.unique(time_arr)).min()

    (aa, bb, cc,
     omega, tau, freq) = get_clean_spectrum_PressRybicki(time_arr,
                                                         flux_arr,
                                                         sig_arr,
                                                         dt,
                                                         max_components=components)


    model = np.zeros(len(segment[0]))
    for ix in range(min(components, len(aa))):
        model += cc[ix]
        arg = omega[ix]*(segment[0]-time_arr.min()-tau[ix])
        model += aa[ix]*np.cos(arg)
        model += bb[ix]*np.sin(arg)

    offset_num = ((segment[1]-model)/np.power(segment[2],2)).sum()
    offset_denom = (1.0/np.power(segment[2],2)).sum()
    offset = offset_num/offset_denom

    #print 'did segment time %e n_seg %d off %e' % \
    #(time.time()-t_start, len(segment_list), offset)

    #print '%e %e' % (omega[0], cc[0])
    #print np.median(time_arr)

    with open('kplr%s_segment_%d.txt' % (args.objid, i_seg), 'w') as out_file:
        for tt, ff, ss in zip(segment[0], segment[1]-offset, segment[2]):
            out_file.write('%.12e %e %e\n' % (tt, ff, ss))

    time_arr = np.append(time_arr, segment[0])
    flux_arr = np.append(flux_arr, segment[1]-offset)
    sig_arr = np.append(sig_arr, segment[2])


with open('kplr%s_stitched.txt' % args.objid, 'w') as out_file:
    for tt, ff, ss in zip(time_arr, flux_arr, sig_arr):
        out_file.write('%.12e %e %e\n' % (tt, ff, ss))


dt = 0.1*np.diff(np.unique(time_arr)).min()

(aa, bb, cc,
 omega, tau, freq) = get_clean_spectrum_PressRybicki(time_arr,
                                                     flux_arr,
                                                     sig_arr,
                                                     dt, gain=0.2,
                                                     max_components=10)

model3 = np.zeros(len(time_arr))
model_max = np.zeros(len(time_arr))
print 'final component count %d' % len(aa)
for ix in range(min(len(aa), 50)):
    xx = omega[ix]*(time_arr-time_arr.min()-tau[ix])
    if ix<3:
        model3 += cc[ix]
        model3 += aa[ix]*np.cos(xx)
        model_max += bb[ix]*np.sin(xx)
    model_max += cc[ix]
    model_max += aa[ix]*np.cos(xx)
    model_max += bb[ix]*np.sin(xx)

with open(args.outfile, 'w') as out_file:
    for tt, ff, ss, m3, m_max in zip(time_arr, flux_arr, sig_arr, model3, model_max):
        out_file.write('%.12e %e %e %e %e\n' % (tt, ff, ss, m3, m_max))


print 'smoothing data took ',time.time()-t_start

