from __future__ import with_statement

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import os

def make_2d_histogram(xx, yy, dx, dy):
    """
    returns indices and counts of unique points on the map
    """
    i_color1 = np.round(xx/dx).astype(int)
    i_color2 = np.round(yy/dy).astype(int)
    dex_reverse = np.array([i_color1, i_color2])
    dex_arr = dex_reverse.transpose()
    # see http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
    dex_raw = np.ascontiguousarray(dex_arr).view(np.dtype((np.void, dex_arr.dtype.itemsize*dex_arr.shape[1])))
    _, unique_rows, unique_counts = np.unique(dex_raw, return_index=True, return_counts=True)

    return unique_rows, unique_counts


def plot_color(xx, yy, dx, dy):
    dexes, cts = make_2d_histogram(xx, yy, dx, dy)
    sorted_dex = np.argsort(cts)
    dexes = dexes[sorted_dex]
    cts = cts[sorted_dex]
    plt.scatter(xx[dexes], yy[dexes], c=cts, s=5,
                cmap=plt.cm.gist_ncar, edgecolor='')

    plt.colorbar()


def plot_color_mesh(xx, yy, dx, dy, vmin=None, vmax=None):
    i_x_arr = np.round((xx-xx.min())/dx).astype(int)
    i_y_arr = np.round((yy-yy.min())/dy).astype(int)
    new_x = i_x_arr*dx
    new_y = i_y_arr*dy
    dex_list, ct_list = make_2d_histogram(new_x, new_y, dx, dy)

    if i_x_arr.min()<0 or i_y_arr.min()<0:
        raise RuntimeError('negative dex %e %d %e %d' %
                           (xx.min(), i_x_arr.min(), yy.min(), i_y_arr.min()))

    x_mesh=np.arange(xx.min(),xx.max()+0.1,dx)
    y_mesh=np.arange(yy.min(),yy.max()+0.1,dy)
    x_mesh,y_mesh = np.meshgrid(x_mesh,y_mesh,indexing='xy')
    print 'x mesh shape ',x_mesh.shape
    z_mesh = np.zeros(shape=x_mesh.shape, dtype=int)
    ct_1000b = 0

    for dex, ct in zip(dex_list, ct_list):
        ix = i_x_arr[dex]
        iy = i_y_arr[dex]
        z_mesh[iy][ix] += ct

    z_mesh = np.ma.masked_where(z_mesh==0,z_mesh)
    plt.pcolormesh(x_mesh,y_mesh,z_mesh, vmin=vmin, vmax=vmax)
                   #norm=matplotlib.colors.LogNorm(vmin=1.0,
                   #                               vmax=1.2e6))
    plt.colorbar(label='sources per pixel')

    ct_1000 = 0
    big_min = np.round((2.8-xx.min())/dx).astype(int)
    big_max = x_mesh.shape[0]

    print 'big x ',xx.min(),big_min,big_max,ct_1000,ct_1000b


def plot_color_mesh_set_color(xx, yy, color, dx, dy, vmin=None, vmax=None,
                              color_label=None):
    i_x_arr = np.round((xx-xx.min())/dx).astype(int)
    i_y_arr = np.round((yy-yy.min())/dy).astype(int)
    new_x = i_x_arr*dx
    new_y = i_y_arr*dy
    dex_list, ct_list = make_2d_histogram(new_x, new_y, dx, dy)

    if i_x_arr.min()<0 or i_y_arr.min()<0:
        raise RuntimeError('negative dex %e %d %e %d' %
                           (xx.min(), i_x_arr.min(), yy.min(), i_y_arr.min()))

    x_mesh=np.arange(xx.min(),xx.max()+0.1,dx)
    y_mesh=np.arange(yy.min(),yy.max()+0.1,dy)
    x_mesh,y_mesh = np.meshgrid(x_mesh,y_mesh,indexing='xy')
    z_mesh = np.ones(shape=x_mesh.shape, dtype=float)*(-999.0)

    for dex in dex_list:
        ix = i_x_arr[dex]
        iy = i_y_arr[dex]

        xmin = xx.min()+ix*dx-0.5*dx
        xmax = xx.min()+ix*dx+0.5*dx
        ymin = yy.min()+iy*dy-0.5*dy
        ymax = yy.min()+iy*dy+0.5*dy

        valid_pts = np.where(np.logical_and(xx>=xmin,
                             np.logical_and(xx<=xmax,
                             np.logical_and(yy>=ymin,yy<=ymax))))

        color_valid = color[valid_pts]

        if len(color_valid)>0:
            val = np.median(color_valid)
            z_mesh[iy][ix] = val

    z_mesh = np.ma.masked_where(z_mesh<-990.0,z_mesh)
    plt.pcolormesh(x_mesh,y_mesh,z_mesh, vmin=vmin, vmax=vmax)
                   #norm=matplotlib.colors.LogNorm(vmin=1.0,
                   #                               vmax=1.2e6))
    plt.colorbar(label=color_label)


mag_dict = {}
with open('data/kepler_stellar17.csv', 'r') as in_file:
    for line in in_file:
        params = line.strip().split('|')
        try:
            mag = float(params[41])
            name_dex = params[0]
            name = 'kplr%9d_lc.txt' % int(name_dex)
            name = name.replace(' ','0')
            #print name
            mag_dict[name] = mag
        except:
            pass

lc_file = os.path.join('workspace', 'validation_170705', 'lc_params_master.txt')

chisq_dof = []
full_models = []
omega_max = []
amp_max = []

with open(lc_file, 'r') as input_file:
    for line in input_file:
        if line[0] == '#':
            continue
        params = line.strip().split()
        n_c = int(params[3])
        if n_c == 0:
            continue
        if params[0] not in mag_dict:
            continue
        n_t = int(params[1])
        chisq = float(params[3+n_c])
        chisq_dof.append(chisq/n_t)
        aa = []
        bb = []
        cc = []
        omega = []
        tau = []
        median = float(params[4+n_c])
        base_dex = 5+n_c
        for i_c in range(n_c):
            local_base = base_dex+5*i_c
            aa.append(float(params[local_base]))
            bb.append(float(params[local_base+1]))
            cc.append(float(params[local_base+2]))
            omega.append(float(params[local_base+3]))
            tau.append(float(params[local_base+4]))

        aa = np.array(aa)
        bb = np.array(bb)
        cc = np.array(cc)
        omega = np.array(omega)
        tau = np.array(tau)
        model = {}
        model['median'] = median
        model['aa'] = aa
        model['bb'] = bb
        model['cc'] = cc
        model['omega'] = omega
        model['tau'] = tau
        model['name'] = params[0]
        full_models.append(model)

t        amp = np.sqrt(aa*aa+bb*bb)
        """
        amp_dexes = np.argsort(amp)
        if len(amp)>5:
            max_dex = amp_dexes[-6]
        else:
            max_dex = amp_dexes[0]
        """
        max_dex = np.argmax(amp)
        amp_max.append(amp[max_dex])
        omega_max.append(omega[max_dex])

amp_max = np.array(amp_max)
omega_max = np.array(omega_max)
full_models = np.array(full_models)
chisq_dof = np.array(chisq_dof)

chisq_sorted = np.sort(chisq_dof)

chisq_cut = chisq_sorted[95*len(chisq_sorted)/100]
chisq_cut = 700.0
#chisq_cut = 3.0

print 'chisq_cut is ',chisq_cut

chisq_valid = np.where(chisq_dof<=chisq_cut)
amp_max = amp_max[chisq_valid]
omega_max = omega_max[chisq_valid]
full_models = full_models[chisq_valid]
chisq_dof = chisq_dof[chisq_valid]

period_max = 2.0*np.pi/omega_max

period_cut = 1000.0
period_valid = np.where(period_max<=period_cut)

amp_max = amp_max[period_valid]
omega_max = omega_max[period_valid]
full_models = full_models[period_valid]
chisq_dof = chisq_dof[period_valid]
period_max = period_max[period_valid]

print '%d lc after chisq and period cut' % len(full_models)

median_arr = []
for model in full_models:
    median_arr.append(model['median'])
median_arr = np.array(median_arr)

print 'min amp_max ',amp_max.min()
print 'median_arr min ',median_arr.min()

median_valid = np.where(median_arr>0.0)
amp_max = amp_max[median_valid]
omega_max = omega_max[median_valid]
full_models = full_models[median_valid]
chisq_dof = chisq_dof[median_valid]
period_max = period_max[median_valid]
median_arr = median_arr[median_valid]

print 'after positive median cut ',len(median_arr)

mag_meas = []
mag_th = []
with open('mag_to_flux.txt', 'w') as out_file:
    out_file.write('# mag flux\n')
    for model in full_models:
        name = model['name']
        median = model['median']
        if name in mag_dict:
            mag = mag_dict[name]
            m_th = -2.5*np.log10(median)
            mag_meas.append(mag_dict[name])
            mag_th.append(m_th)
            out_file.write('%e %e %e\n' % (mag, median, m_th-mag))

mag_meas = np.array(mag_meas)
mag_th = np.array(mag_th)

mag_offset = np.mean(mag_meas-mag_th)

median_th = []
for model in full_models:
    name = model['name']
    mag = mag_dict[name]
    flux = np.power(10.0,-0.4*(mag-mag_offset))
    median_th.append(flux)

median_th = np.array(median_th)

print 'find full amplitudes'
print amp_max.min(),amp_max.max()

mag_amp = 2.5*np.log10(1.0+amp_max/median_th)

already_calculated = []

if os.path.exists('name_to_amp.txt'):
    with open('name_to_amp.txt', 'r') as in_file:
        for line in in_file:
            params = line.strip().split()
            already_calculated.append(params[0])


# make the name_to_amp file
if len(already_calculated) != len(full_models):
    lines_to_write = []
    amp_th = np.zeros(len(amp_max))
    mag_amp_th = np.ones(len(amp_max))*(-999.0)
    t_model = np.arange(0.0,17.0*90.0,0.1)
    for i_model, model in enumerate(full_models):
        if model['name'] in already_calculated:
            continue
        f_model = np.zeros(len(t_model))
        for ix in range(len(model['aa'])):
            arg = model['omega'][ix]*(t_model-model['tau'][ix])
            f_model += model['aa'][ix]*np.cos(arg)
            f_model += model['bb'][ix]*np.sin(arg)

        rms = np.sqrt(np.mean((f_model-np.median(f_model))**2))
        amp_th[i_model] = rms
        mag_amp_th[i_model] = 2.5*np.log10(1.0+amp_th[i_model]/median_th[i_model])

        lines_to_write.append('%s %e %e %e\n' % (model['name'], amp_th[i_model], mag_amp_th[i_model], chisq_dof[i_model]))

        if i_model%1000 == 0:
            greater = np.where(mag_amp_th[:i_model+1]>0.1)
            greater_old = np.where(mag_amp[:i_model+1]>0.1)
            ratio = amp_th[:i_model+1]/amp_max[:i_model+1]
            print('%d med ratio %.2e max ratio %.2e med amp %.2e %.2e  max amp %.2e %.2e %d %d' %
                                (i_model,
                                np.median(ratio),
                                ratio.max(),
                                np.median(mag_amp_th[:i_model+1]),
                                np.median(mag_amp[:i_model+1]),
                                mag_amp_th[:i_model+1].max(),
                                mag_amp[:i_model+1].max(),
                                len(greater[0]),len(greater_old[0])))


    if len(lines_to_write)>0:
        with open("name_to_amp.txt", "a") as amp_file:
            for line in lines_to_write:
                amp_file.write(line)


mag_amp_dict = {}
flux_amp_dict = {}
with open('name_to_amp.txt', 'r') as input_file:
    for line in input_file:
        if line[0]=='#':
            continue
        params = line.strip().split()
        flux_amp_dict[params[0]] = float(params[1])
        mag_amp_dict[params[0]] = float(params[2])

mag_amp_th = []
for model in full_models:
    mag_amp_th.append(mag_amp_dict[model['name']])
mag_amp_th = np.array(mag_amp_th)

plt.figsize = (30,30)
plt.subplot(2,1,1)
plt.title('amplitude is sqrt(a^2+b^2) of first component', fontsize=10)
plot_color_mesh(np.log10(period_max), np.log10(mag_amp), 0.05, 0.05,
                vmin=0.0,vmax=1000.0)
plt.xlabel('log10(period in days)')
plt.ylabel('log10(amplitude in mags)')
plt.axvline(np.log10(90.0), linestyle='--', color='r')
plt.axvline(np.log10(45.0), linestyle='--', color='r')
#plt.axvline(np.log10(22.5), linestyle='--', color='r')
#plt.axvline(np.log10(1.0/48.0), linestyle='--',color='r')
plt.axvline(-0.6, linestyle='--', color='r')
plt.axhline(-2, linestyle='--', color='r')
plt.axhline(-3, linestyle='--', color='r')

plt.ylim(-6,-1)
plt.xlim(-2,3)

plt.subplot(2,1,2)
plt.title("amplitude is rms of sum of all components", fontsize=10)
plot_color_mesh(np.log10(period_max), np.log10(mag_amp_th), 0.05, 0.05,
                vmin=0.0,vmax=1000.0)
plt.xlabel('log10(period in days)')
plt.ylabel('log10(amplitude in mags)')
plt.axvline(np.log10(90.0), linestyle='--', color='r')
plt.axvline(np.log10(45.0), linestyle='--', color='r')
#plt.axvline(np.log10(22.5), linestyle='--', color='r')
#plt.axvline(np.log10(1.0/48.0), linestyle='--',color='r')
plt.axvline(-0.6, linestyle='--', color='r')
plt.axhline(-2, linestyle='--', color='r')
plt.axhline(-3, linestyle='--', color='r')

plt.ylim(-6,-1)
plt.xlim(-2,3)

plt.tight_layout()

plt.savefig('period_amp_mag.png')
plt.close()

"""
curious = np.where(np.logical_and(np.abs(np.log10(period_max)-np.log10(45.0))<0.1,
                                         np.log10(mag_amp)<-4.0))

print 'curious ',curious

curious_models = full_models[curious]
for model in curious_models:
    print model['name']

"""

#kic_dtype = np.dtype([('id', int), ('ra', float), ('dec', float), ('u', float),
#                         ('g', float), ('r', float), ('i', float), ('z', float),
#                         ('dist', float), ('teff_kic', float), ('teff_stellar', float),
#                         ('Av', float), ('ebv', float)])
kic_file = 'KIC/kic_data.txt'

mag_dict = {}
with open(kic_file, 'r') as in_file:
    for line in in_file:
        if line[0] == '#':
            continue
        params = line.strip().split()
        name = 'kplr%9d_lc.txt' % int(params[0])
        name = name.replace(' ','0')
        mag_list = [float(params[3]), float(params[4]), float(params[5]),
                    float(params[6]), float(params[7])]
        if mag_list[1]>0.0 and mag_list[2]>0.0:
            mag_dict[name] = mag_list

plt.figsize = (30,30)

period_to_plot = []
color_to_plot = []
amp_mag_to_plot = []
for ix, model in enumerate(full_models):
    if model['name'] in mag_dict:
        period_to_plot.append(period_max[ix])
        color_to_plot.append(mag_dict[model['name']][1]-mag_dict[model['name']][2])
        amp_mag_to_plot.append(mag_amp[ix])

period_to_plot = np.array(period_to_plot)
color_to_plot = np.array(color_to_plot)
amp_mag_to_plot = np.array(amp_mag_to_plot)

print 'plotting ',len(period_to_plot),len(mag_dict)

plt.subplot(2,1,1)
plot_color_mesh(color_to_plot, np.log10(period_to_plot), 0.05, 0.05)
plt.ylabel('log10(period in days')
plt.xlabel('g-r')
plt.xlim(-0.5,2)

print 'making plot of color vs period shaded by amp'
plt.subplot(2,1,2)
plot_color_mesh_set_color(color_to_plot, np.log10(period_to_plot), amp_mag_to_plot,
                          0.05, 0.05, color_label='median amplitude (mag)')
plt.ylabel('log10(period in days')
plt.xlabel('g-r')
plt.xlim(-0.5,2)
plt.ylim(-1.5, 0.25)


plt.tight_layout()
plt.savefig('color_vs_period.png')
plt.close()

big_amp = np.where(mag_amp>0.1)
chisq_big_amp = chisq_dof[big_amp]
model_big_amp = full_models[big_amp]
sorted_dex = np.argsort(chisq_big_amp)

for ix in sorted_dex[-6:]:
    print model_big_amp[ix]['name'],chisq_big_amp[ix]
