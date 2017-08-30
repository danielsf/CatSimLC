from __future__ import with_statement
import os

hrs = 15

for ix in range(30):
    #if ix==0 or ix==20:
    #    continue
        # because I already ran these two datasets when testing Cori

    script_name = 'cori_batch_script_%d.sl' % ix
    with open(script_name, 'w') as out_file:
        out_file.write('#!/bin/bash -l\n')
        out_file.write('\n')
        out_file.write('#SBATCH -p regular\n')
        out_file.write('#SBATCH -N 1\n')
        out_file.write('#SBATCH -t %d:00:00\n' % hrs)
        out_file.write('#SBATCH -C haswell\n')
        out_file.write('#SBATCH -e cori_batch_%d_stderr.txt\n' % ix)
        out_file.write('#SBATCH -o cori_batch_%d_stdout.txt\n' % ix)
        out_file.write('\n')
        out_file.write('#SBATCH -A m1727\n')
        out_file.write('\n')
        out_file.write('cd $HOME/CatSimLC\n')
        out_file.write('\n')
        out_file.write("echo 'starting now'\n")
        out_file.write("date\n")
        out_file.write("\n")
        out_file.write("module load python/2.7-anaconda\n")
        out_file.write("\n")
        out_file.write("python clean_multiprocess.py \\\n")
        out_file.write("--list $SCRATCH/kepler_lightcurves/lc_lists/lc_list_%d.txt \\\n" % ix)
        out_file.write("--in_dir $SCRATCH/kepler_lightcurves/lc_master_3/ \\\n")
        out_file.write("--do_stitch True --stitch_dir $SCRATCH/kepler_lightcurves/stitch_170824/ \\\n")
        out_file.write("--out_file lc_param_dir_170824/lc_params_batch_%d \\\n" % ix)
        out_file.write("--log_file log_dir/cori_log_batch_%d \\\n" % ix)
        out_file.write("--dt 0.1 --write_every 100 --cache_fft True --n_p 32 --sleep 10\n")
        out_file.write("\n")
        out_file.write("echo 'ending now'\n")
        out_file.write("date\n")
