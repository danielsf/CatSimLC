from __future__ import with_statement
import os
import numpy as np
import astropy.io.fits as fits

def linear_bridge(left_x, left_y, left_s,
                  right_x, right_y, right_s):

    theta_1=(left_x/np.power(left_s,2)).sum()/(1.0/np.power(left_s,2)).sum()
    theta_2=(right_x/np.power(right_s,2)).sum()/(1.0/np.power(right_s,2)).sum()
    mm_denom = ((np.power(left_x,2)-left_x*theta_1)/np.power(left_s,2)).sum()
    mm_denom += ((np.power(right_x,2)-right_x*theta_2)/np.power(right_s,2)).sum()

    theta_1=(left_y/np.power(left_s,2)).sum()/(1.0/np.power(left_s,2)).sum()
    theta_2=(right_y/np.power(right_s,2)).sum()/(1.0/np.power(right_s,2)).sum()
    mm_num = ((left_y*left_x-left_x*theta_1)/np.power(left_s,2)).sum()
    mm_num += ((right_y*right_x-right_x*theta_2)/np.power(right_s,2)).sum()

    mm = mm_num/mm_denom

    bb = ((left_y-mm*left_x)/np.power(left_s,2)).sum()
    bb = bb/((1.0/np.power(left_s,2)).sum())

    aa = ((right_y-mm*right_x)/np.power(right_s,2)).sum()
    aa = aa/((1.0/np.power(right_s,2)).sum())
    aa -= bb

    return mm, bb, aa

x=np.arange(10.0, 20.0, 1.0)
y=x*2.3-7.2
x2 = np.arange(40.0, 50.0, 10.0)
y2 = x2*2.3-7.2-9.3

mm, bb, aa = linear_bridge(x,y,np.ones(len(x)),
                           x2,y2,np.ones(len(x)))

print mm,bb,aa
#mm=2.3
#bb=-7.2
#aa=-9.3
chisq = np.power(y-mm*x-bb,2).sum()
chisq += np.power(y2-mm*x2-bb-aa,2).sum()
print 'chisq ',chisq


overhang = 100

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

line_ct = 0

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
        else:
            mm, bb, aa = linear_bridge(mjd[-overhang:],
                                       flux[-overhang:],
                                       sig[-overhang:],
                                       valid_mjd[:overhang],
                                       valid_flux[:overhang],
                                       valid_sig[:overhang])

            with open('line_file_%d.txt' % line_ct, 'w') as out_file:
                for t in mjd[-overhang:]:
                    out_file.write('%e %e\n' % (t,mm*t+bb))
                for t in valid_mjd[:overhang]:
                    out_file.write('%e %e\n' % (t,mm*t+bb))
            line_ct += 1
            mjd = np.append(mjd, valid_mjd)
            flux = np.append(flux, valid_flux-aa)
            sig = np.append(sig, valid_sig)
            target_flux = flux[-1]

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
