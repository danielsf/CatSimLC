#!/bin/bash -l

#SBATCH -p regular
#SBATCH -N 1
#SBATCH -t 15:00:00
#SBATCH -C haswell
#SBATCH -e cori_batch_20_stderr.txt
#SBATCH -o cori_batch_20_stdout.txt

#SBATCH -A m1727

cd $HOME/CatSimLC

echo 'starting now'
date

module load python/2.7-anaconda

python clean_multiprocess.py \
--list $SCRATCH/kepler_lightcurves/lc_lists/lc_list_20.txt \
--in_dir $SCRATCH/kepler_lightcurves/lc_master_3/ \
--do_stitch True --stitch_dir $SCRATCH/kepler_lightcurves/stitch_170811/ \
--out_file lc_param_dir_170811/lc_params_batch_20 \
--log_file log_dir/cori_log_batch_20 \
--dt 0.1 --write_every 100 --cache_fft True --n_p 32 --sleep 10

echo 'ending now'
date
