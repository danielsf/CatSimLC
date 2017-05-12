# This script will use wget to download all of the Kepler light
# curves from stsci.  It will put the resulting FITS files in
# directories named lc_N where N refers to the observing
# quarter and ranges from 0 to 17.

for ii in {0..17};
do
    mkdir lc_${ii}
    cd lc_${ii}
    nohup nice wget -q -nH --cut-dirs=6 -r -l0 -c -N -np \
-R 'index*' -erobots=off \
http://archive.stsci.edu/pub/kepler/lightcurves/tarfiles/Q${ii}_public/ \
>& wget_stdout_${ii}.txt

    cd ..
done
