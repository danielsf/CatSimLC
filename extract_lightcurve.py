from __future__ import with_statement
import os
import numpy as np
import astropy.io.fits as fits

obj_id = '009449503'
#obj_id = '008097275'

data_dir = os.path.join('workspace', 'data')
list_of_files = os.listdir(data_dir)
list_of_files.sort()

obj_name = None
ra = None
dec = None
tunit = None

mjd = None
flux = None
sig = None

target_flux = None

for file_name in list_of_files:
    if obj_id in file_name and file_name.endswith('fits'):
        data = fits.open(os.path.join(data_dir, file_name))
        if obj_name is None:
            obj_name = data[0].header['OBJECT']
            ra = data[1].header['RA_OBJ']
            dec = data[1].header['DEC_OBJ']
            tunit = data[1].header['TUNIT1']
        else:
            try:
                assert obj_name == data[0].header['OBJECT']
                assert tunit == data[1].header['TUNIT1']
                assert ra == data[1].header['RA_OBJ']
                assert dec == data[1].header['DEC_OBJ']
            except:
                print file_name
                print obj_name, data[0].header['OBJECT']
                print tunit, data[1].header['TUNIT1']
                print ra, data[1].header['RA_OBJ']
                print dec, data[1].header['DEC_OBJ']
                raise

        valid_dexes = np.where(np.logical_not(np.isnan(data[1].data['PDCSAP_FLUX'])))
        valid_flux = data[1].data['PDCSAP_FLUX'][valid_dexes]
        valid_mjd = data[1].data['TIME'][valid_dexes]
        valid_sig = data[1].data['PDCSAP_FLUX_ERR'][valid_dexes]

        print valid_mjd.min(),valid_mjd.max()

        if mjd is None:
            mjd = valid_mjd
            flux = valid_flux
            sig = valid_sig
            target_flux = valid_flux[-1]
        else:
            mjd = np.append(mjd, valid_mjd)
            #flux = np.append(flux, valid_flux)
            flux = np.append(flux, valid_flux+target_flux-valid_flux[0])
            sig = np.append(sig, valid_sig)
            target_flux = flux[-1]

print mjd
print flux
print sig

valid = np.where(np.logical_not(np.isnan(mjd)))
mjd= mjd[valid]
flux = flux[valid]
sig = sig[valid]

valid = np.where(np.logical_not(np.isnan(flux)))
mjd=mjd[valid]
flux=flux[valid]
sig=sig[valid]

valid = np.where(np.logical_not(np.isnan(sig)))
mjd=mjd[valid]
flux=flux[valid]
sig=sig[valid]
dexes = np.argsort(mjd)
delta = np.diff(mjd[dexes]).min()
print 'delta %.4e' % delta
#print type(delta),type(mjd.max()),type(mjd.min())
#print (mjd.max()-mjd.min())/delta
print 'n %.4e' % ((mjd.max()-mjd.min())/delta)
with open('kplr%s_matched_lc.txt' % obj_id, 'w') as output_file:
    for tt, ff, ss in zip(mjd[dexes], flux[dexes], sig[dexes]):
        output_file.write('%.12e %.5e %.5e\n' % (tt, ff, ss))
