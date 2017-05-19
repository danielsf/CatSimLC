"""
This script will create a look up table that will map LSST colors to
the ratio of Kepler flux to the sum of g, r, i, z fluxes.  This should
allow us to approximate the Kepler magnitudes of objects in CatSim without
having to read in their SEDs and integrate over the Kepler bandpass directly.
"""

from __future__ import with_statement
import numpy as np
import os

sed_wavelen_min = 250.0
sed_wavelen_max = 1500.0

bp_wavelen_min = 300.0
bp_wavelen_max = 1150.0

from lsst.sims.photUtils import cache_LSST_seds

cache_LSST_seds(wavelen_min=sed_wavelen_min, wavelen_max=sed_wavelen_max)

from lsst.sims.photUtils import Bandpass, BandpassDict, Sed, SedList

from lsst.utils import getPackageDir
throughput_dir = getPackageDir('throughputs')
bp_list = []
bp_name_list = []

for name in ('u', 'g', 'r', 'i', 'z', 'y'):
    file_name = os.path.join(throughput_dir, 'baseline', 'total_%s.dat' % name)
    bp = Bandpass()
    bp.readThroughput(file_name)
    bp_list.append(bp)
    bp_name_list.append(name)

file_name =  'data/kepler_throughput.txt'
dtype = np.dtype([('w', float), ('sb', float)])
data = np.genfromtxt(file_name, dtype=dtype)

w_pad = np.arange(bp_wavelen_min, data['w'][0]-1.0, 1.0)
s_pad = np.zeros(len(w_pad))

wav = np.append(w_pad, data['w'])
sb = np.append(s_pad, data['sb'])

w_pad = np.arange(wav[-1]+1.0, bp_wavelen_max, 1.0)
s_pad = np.zeros(len(w_pad))

wav = np.append(wav, w_pad)
sb = np.append(sb, s_pad)

print len(wav),len(sb)
bp = Bandpass(wavelen=wav, sb=sb)

bp_list.append(bp)
bp_name_list.append('kepler')

bp_dict = BandpassDict(bp_list, bp_name_list)

sed_dir = getPackageDir('sims_sed_library')
star_dir = os.path.join(sed_dir, 'starSED')

sed_list = None

converter = Sed()

out_name = 'data/lsst_color_to_kepler_grid.txt'

with open(out_name, 'w') as output_file:
    output_file.write('# ebv u-g g-r r-i i-z z-y kep_flux/griz_flux\n')

for dir_name in ('kurucz', 'mlt', 'wDs'):
    print 'processing ',dir_name

    list_of_files = os.listdir(os.path.join(star_dir, dir_name))
    list_of_files = np.array(list_of_files)
    magnorm_list = 15.0*np.ones(len(list_of_files))
    for ebv in np.arange(0.1, 7.0, 0.1):
        print '    ebv ',ebv,dir_name
        av_list = 3.1*np.ones(len(list_of_files))*ebv
        if sed_list is None:
           sed_list = SedList(list_of_files, magnorm_list,
                              galacticAvList=av_list,
                              wavelenMatch=bp_dict.wavelenMatch)
        else:
            sed_list.flush()
            sed_list.loadSedsFromList(list_of_files, magnorm_list,
                                      galacticAvList=av_list)

        flux_list = bp_dict.fluxListForSedList(sed_list).transpose()
        mag_list = converter.magFromFlux(flux_list)

        u_g = mag_list[0]-mag_list[1]
        g_r = mag_list[1]-mag_list[2]
        r_i = mag_list[2]-mag_list[3]
        i_z = mag_list[3]-mag_list[4]
        z_y = mag_list[4]-mag_list[5]
        griz_flux = flux_list[1]+flux_list[2]+flux_list[3]+flux_list[4]
        flux_rat = flux_list[6]/griz_flux

        with open(out_name, 'a') as output_file:
            for ix in range(len(flux_rat)):
                output_file.write('%e %e %e %e %e %e %e\n'
                                  % (ebv, u_g[ix], g_r[ix], r_i[ix],
                                     i_z[ix], z_y[ix], flux_rat[ix]))
