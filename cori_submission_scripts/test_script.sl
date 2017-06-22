#!/bin/bash -l

#SBATCH -p regular
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH -C haswell
#SBATCH -e test_lc_stderr.txt
#SBATCH -o test_lc_stdout.txt

#SBATCH -A m1727

cd $HOME/CatSimLC

echo 'starting now'
date

module load python/2.7-anaconda

python clean_multiprocess.py --list lc_full_test_list.txt \
--in_dir $SCRATCH/kepler_lightcurves/lc_master_3/ \
--do_stitch True --stitch_dir $SCRATCH/kepler_lightcurves/stitch_170619/ \
--out_file lc_params_170619 --log_file log_dir/lc_log_170619 --dt 0.1 \
--write_every 100 --cache_fft True --n_p 32 --sleep 10

echo 'ending now'
date
