import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

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
        raise RuntimeError('negative dex')

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

catsim_dtype = np.dtype([('sedname', str, 300), ('teff', float), ('feh', float), ('logg', float), ('m', float),
                   ('r_abs', float), ('norm', float),
                   ('u', float), ('g', float), ('r', float), ('i', float), ('z', float)])

kepler_dtype = np.dtype([('id', int), ('ra', float), ('dec', float),
                         ('u', float), ('g', float), ('r', float),
                         ('i', float), ('z', float), ('dist', float),
                         ('teff', float), ('teff_stellar', float)])

kep_file = 'KIC/kic_data.txt'
catsim_file = 'catsim_star_data_same_pointing_cutoff.txt'
suffix = 'kic'

raw_kep_data = np.genfromtxt(kep_file, dtype=kepler_dtype)

valid_dex = np.where(np.logical_and(raw_kep_data['dist']>0.0, raw_kep_data['teff']>0.0))
kep_data = raw_kep_data[valid_dex]
print kep_data['r'].min(),kep_data['g'].min()
valid_dex = np.where(np.logical_and(kep_data['r']>0.0, kep_data['g']>0.0))
kep_data = kep_data[valid_dex]
kep_r_abs = kep_data['r']-5.0*np.log10(kep_data['dist']/10.0)
kep_color = kep_data['g']-kep_data['r']

catsim_data = np.genfromtxt(catsim_file, dtype=catsim_dtype)

#cut_dex = np.where(catsim_data['r_abs']>4.0)
#catsim_data= catsim_data[cut_dex]

catsim_color = catsim_data['g']-catsim_data['r']

color_min = min(kep_color.min(), catsim_color.max())
color_max = max(kep_color.max(), catsim_color.max())
m_min = min(kep_r_abs.min(), catsim_data['r_abs'].min())
m_max = max(kep_r_abs.max(), catsim_data['r_abs'].max())

m_min = 0.0
m_max = 8.0
color_min=-0.5
color_max=1.5

#t_ticks = np.arange(np.round(color_min/1000.0)*1000.0, np.round(color_max/1000.0)*1000.0, 2000.0)
#t_labels = ['%d' % tt for tt in t_ticks]

dm = 0.02
dfeh = 0.1
dlogg = 0.1

plt.figsize = (30, 30)
plt.subplot(2,2,1)

#plot_color(catsim_color, catsim_data['r_abs'], dm, dm)
plot_color_mesh(catsim_color, catsim_data['r_abs'], dm, dm)
plt.title('CatSim')
plt.xlabel('g-r')
plt.ylabel('r')
plt.xlim(color_min, color_max)
#plt.xticks(t_ticks, t_labels, fontsize=7)
plt.ylim(m_min, m_max)
#plt.gca().invert_xaxis()
plt.gca().invert_yaxis()

print 'made CatSim HR diagram'

plt.subplot(2,2,2)
plt.title('Kepler')
plt.xlabel('g-r')
plt.ylabel('r')
#plot_color(kep_color, kep_r_abs, dm, dm)
plot_color_mesh(kep_color, kep_r_abs, dm, dm)
plt.xlim(color_min, color_max)
#plt.xticks(t_ticks, t_labels, fontsize=7)
plt.ylim(m_min, m_max)
#plt.gca().invert_xaxis()
plt.gca().invert_yaxis()

print 'made Kepler HR diagram'

plt.subplot(2,2,3)
m_min=2.0
m_max=8.0
color_min=-0.5
color_max=1.5
counts, xbins, ybins = np.histogram2d(catsim_color, catsim_data['r_abs'], bins=100)
catsim = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
            colors='r', alpha=0.5)

print 'catsim counts shape ',counts.shape

counts,xbins,ybins = np.histogram2d(kep_color, kep_r_abs, bins=200)
kep = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
            colors='blue', alpha=0.5)

print 'kepler counts shape ',counts.shape

plt.xlabel('g-r')
plt.ylabel('r')

#m_min=-5
#m_max=10

#t_min=3000.0
#t_max=7000.0
#t_ticks = np.arange(np.round(t_min/1000.0)*1000.0, np.round(t_max/1000.0)*1000.0, 1000.0)
#t_labels = ['%d' % tt for tt in t_ticks]


plt.xlim(color_min, color_max)
#plt.xticks(t_ticks, t_labels, fontsize=7)
plt.ylim(m_min, m_max)

plt.text(color_max-0.9*(color_max-color_min),m_min+0.25*(m_max-m_min),
         'Blue is Kepler; Red is CatSim')

#plt.gca().invert_xaxis()
plt.gca().invert_yaxis()

"""
plt.subplot(2,2,4)
valid_dex = np.where(catsim_data['feh']>-900.0)
plt.hist(kep_data['feh'], bins=1000, color='b', zorder=1, edgecolor='b', normed=True)
plt.hist(catsim_data['feh'][valid_dex], bins=1000, color='r', zorder=2, edgecolor='r',
         alpha=0.5, normed=True)

plt.xlabel('FeH', fontsize=10)
"""

plt.tight_layout()
plt.savefig('kepler_hr_diagram_%s.png' % suffix)
plt.close()


plt.figsize = (30,30)

plt.subplot(1,2,1)
plt.hist(kep_r_abs, bins=1000, color='b', zorder=1, edgecolor='b', normed=True)
plt.hist(catsim_data['r_abs'], bins=1000, color='r', zorder=2, edgecolor='r',
         alpha=0.1, normed=True)

plt.xlabel('r')
plt.title('blue is Kepler; red is CatSim')

plt.subplot(1,2,2)
plt.hist(kep_color, bins=1000, color='b', zorder=1, edgecolor='b', normed=True)
plt.hist(catsim_color, bins=1000, color='r', zorder=2, edgecolor='r',
         alpha=0.1, normed=True)
plt.xlabel('g-r')

plt.tight_layout()
plt.savefig('kepler_1d_dist_%s.png' % suffix)
plt.close()

print 'mean CatSim color %.3e' % np.mean(catsim_color)
print 'mean Kepler color %.3e' % np.mean(kep_color)
print 'mean CatSim r %.3e' % np.mean(catsim_data['r_abs'])
print 'mean Kepler r %3e' % np.mean(kep_r_abs)

#print 'catsim feh range ',catsim_data['feh'][valid_dex].min(),catsim_data['feh'].max()

# now plot the distributions in principal colors

s2_coeffs = {'u':-0.249, 'g':0.794, 'r':-0.555, 'offset':0.234}
w2_coeffs = {'g':-0.227, 'r':0.792, 'i':-0.567, 'offset':0.05}
x2_coeffs = {'g':0.707, 'r':-0.707, 'offset':-0.988}
y2_coeffs = {'r':-0.270, 'i':0.8, 'z':-0.534, 'offset':0.054}

p2_coeff_dict = {'s':s2_coeffs, 'w':w2_coeffs, 'x':x2_coeffs, 'y':y2_coeffs}

s1_coeffs = {'u':0.91, 'g':-0.495, 'r':-0.415, 'offset':-1.28}
w1_coeffs = {'g':0.928, 'r':-0.556, 'i':-0.372, 'offset':-0.425}
x1_coeffs = {'r':1.0, 'i':-1.0, 'offset':0.0}
y1_coeffs = {'r':0.895, 'i':-0.448, 'z':-0.447, 'offset':-0.6}

p1_coeff_dict = {'s':s1_coeffs, 'w':w1_coeffs, 'x':x1_coeffs, 'y':y1_coeffs}

plt.figsize = (30,30)
trim = 50
for i_fig, color_name in enumerate(['s', 'w', 'x', 'y']):
    plt.subplot(2,2,i_fig+1)
    coeffs = p2_coeff_dict[color_name]
    print 'raw_kep_data ',len(raw_kep_data)
    kep_data = raw_kep_data
    for tag in coeffs.keys():
        if tag != 'offset':
            valid_dex = np.where(kep_data[tag]>0.0)
            kep_data = kep_data[valid_dex]

    print 'becomes ',len(kep_data)
    kep_color = np.ones(len(kep_data))*coeffs['offset']
    catsim_color = np.ones(len(catsim_data))*coeffs['offset']
    for tag in coeffs.keys():
        if tag == 'offset':
            continue
        kep_color += coeffs[tag]*kep_data[tag]
        catsim_color += coeffs[tag]*catsim_data[tag]

    if i_fig ==0:
        plt.title('blue is Kepler; red is CatSim')
    plt.hist(kep_color, bins=1000, color='b', zorder=1, edgecolor='b', normed=True)
    plt.hist(catsim_color, bins=1000, color='r', zorder=2, edgecolor='r',
             alpha=0.1, normed=True)
    plt.xlabel(color_name)
    #plt.xlim((min(kep_color.min(),catsim_color.min()),max(kep_color.max(),catsim_color.max())))

    kep_color_sorted = np.sort(kep_color)
    catsim_color_sorted = np.sort(catsim_color)

    xmin = min(kep_color_sorted[len(kep_color)/trim], catsim_color_sorted[len(catsim_color)/trim])
    xmax = max(kep_color_sorted[(trim-1)*len(kep_color)/trim], catsim_color_sorted[(trim-1)*len(catsim_color)/trim])

    plt.xlim((xmin,xmax))
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

plt.tight_layout()
plt.savefig('kepler_principal_colors_%s.png' % suffix)
plt.close()

# 2d plots

plt.figsize = (30,30)
trim_2d = trim
for i_fig, color_name in enumerate(['s', 'w', 'x', 'y']):
    plt.subplot(2,2,i_fig+1)
    p2_coeffs = p2_coeff_dict[color_name]
    p1_coeffs = p1_coeff_dict[color_name]
    print 'raw_kep_data ',len(raw_kep_data)
    kep_data = raw_kep_data
    for tag in p2_coeffs.keys():
        if tag != 'offset':
            valid_dex = np.where(kep_data[tag]>0.0)
            kep_data = kep_data[valid_dex]

    for tag in p1_coeffs.keys():
        if tag != 'offset':
            valid_dex = np.where(kep_data[tag]>0.0)
            kep_data = kep_data[valid_dex]

    print 'becomes ',len(kep_data)
    kep_p2 = np.ones(len(kep_data))*p2_coeffs['offset']
    catsim_p2 = np.ones(len(catsim_data))*p2_coeffs['offset']
    for tag in p2_coeffs.keys():
        if tag == 'offset':
            continue
        kep_p2 += p2_coeffs[tag]*kep_data[tag]
        catsim_p2 += p2_coeffs[tag]*catsim_data[tag]

    kep_p1 = np.ones(len(kep_data))*p1_coeffs['offset']
    catsim_p1 = np.ones(len(catsim_data))*p1_coeffs['offset']
    for tag in p1_coeffs.keys():
        if tag == 'offset':
            continue
        kep_p1 += p1_coeffs[tag]*kep_data[tag]
        catsim_p1 += p1_coeffs[tag]*catsim_data[tag]

    counts, xbins, ybins = np.histogram2d(catsim_p1, catsim_p2, bins=100)
    catsim = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                         colors='r', alpha=0.5)


    counts,xbins,ybins = np.histogram2d(kep_p1, kep_p2, bins=200)
    kep = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
                      colors='blue', alpha=0.5)

    if i_fig ==0:
        plt.title('blue is Kepler; red is CatSim')

    plt.ylabel(color_name)
    plt.xlabel('P1(%s)' % color_name)
    #plt.xlim((min(kep_color.min(),catsim_color.min()),max(kep_color.max(),catsim_color.max())))

    kep_p1_sorted = np.sort(kep_p1)
    catsim_p1_sorted = np.sort(catsim_p1)

    xmin = min(kep_p1_sorted[len(kep_p1)/trim_2d], catsim_p1_sorted[len(catsim_p1)/trim_2d])
    xmax = max(kep_p1_sorted[(trim_2d-1)*len(kep_p1)/trim_2d], catsim_p1_sorted[(trim_2d-1)*len(catsim_p1)/trim_2d])

    plt.xlim((xmin,xmax))


    kep_p2_sorted = np.sort(kep_p2)
    catsim_p2_sorted = np.sort(catsim_p2)

    ymin = min(kep_p2_sorted[len(kep_p2)/trim_2d], catsim_p2_sorted[len(catsim_p2)/trim_2d])
    ymax = max(kep_p2_sorted[(trim_2d-1)*len(kep_p2)/trim_2d], catsim_p2_sorted[(trim_2d-1)*len(catsim_p2)/trim_2d])

    plt.ylim((ymin,ymax))
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

plt.tight_layout()
plt.savefig('kepler_principal_colors_2d_%s.png' % suffix)
plt.close()


# now apply p1 restrictions

rmax_dict = {'s': 19.0, 'w':20.0, 'x':19.0, 'y':19.5}
p1_min_dict = {'s':-0.2, 'w':-0.2, 'x':0.8, 'y':0.1}
p1_max_dict = {'s':0.8, 'w':0.6, 'x':1.6, 'y':1.2}

plt.figsize = (30,30)
for i_fig, color_name in enumerate(['s', 'w', 'x', 'y']):
    plt.subplot(2,2,i_fig+1)
    p2_coeffs = p2_coeff_dict[color_name]
    p1_coeffs = p1_coeff_dict[color_name]
    print 'raw_kep_data ',len(raw_kep_data)
    kep_data = raw_kep_data
    for tag in p2_coeffs.keys():
        if tag != 'offset':
            valid_dex = np.where(kep_data[tag]>0.0)
            kep_data = kep_data[valid_dex]

    for tag in p1_coeffs.keys():
        if tag != 'offset':
            valid_dex = np.where(kep_data[tag]>0.0)
            kep_data = kep_data[valid_dex]

    print 'becomes ',len(kep_data)
    kep_color = np.ones(len(kep_data))*p2_coeffs['offset']
    catsim_color = np.ones(len(catsim_data))*p2_coeffs['offset']
    for tag in p2_coeffs.keys():
        if tag == 'offset':
            continue
        kep_color += p2_coeffs[tag]*kep_data[tag]
        catsim_color += p2_coeffs[tag]*catsim_data[tag]

    kep_p1 = np.ones(len(kep_data))*p1_coeffs['offset']
    catsim_p1 = np.ones(len(catsim_data))*p1_coeffs['offset']
    for tag in p1_coeffs.keys():
        if tag == 'offset':
            continue
        kep_p1 += p1_coeffs[tag]*kep_data[tag]
        catsim_p1 += p1_coeffs[tag]*catsim_data[tag]

    rmax = rmax_dict[color_name]
    p1max = p1_max_dict[color_name]
    p1min = p1_min_dict[color_name]

    catsim_valid = np.where(np.logical_and(catsim_data['r']<rmax,
                            np.logical_and(catsim_p1<p1max, catsim_p1>p1min)))

    kep_valid = np.where(np.logical_and(kep_data['r']<rmax,
                         np.logical_and(kep_p1<p1max, kep_p1>p1min)))

    print '%s catsim_valid %d %e %e' % (color_name, len(catsim_valid[0]),catsim_p1.min(),
                                        catsim_p1.max())
    print '%s kep_valid %d %e %e' % (color_name, len(kep_valid[0]), kep_p1.min(),
                                     kep_p1.max())

    catsim_color = catsim_color[catsim_valid]
    kep_color = kep_color[kep_valid]


    if i_fig ==0:
        plt.title('blue is Kepler; red is CatSim')
    if len(kep_color)>0:
        plt.hist(kep_color, bins=1000, color='b', zorder=1, edgecolor='b', normed=True)
    if len(catsim_color)>0:
        plt.hist(catsim_color, bins=1000, color='r', zorder=2, edgecolor='r',
                 alpha=0.1, normed=True)
    plt.xlabel(color_name)
    #plt.xlim((min(kep_color.min(),catsim_color.min()),max(kep_color.max(),catsim_color.max())))

    kep_color_sorted = np.sort(kep_color)
    catsim_color_sorted = np.sort(catsim_color)

    if len(kep_color>0) and len(catsim_color)>0:
        xmin = min(kep_color_sorted[len(kep_color)/trim], catsim_color_sorted[len(catsim_color)/trim])
        xmax = max(kep_color_sorted[(trim-1)*len(kep_color)/trim], catsim_color_sorted[(trim-1)*len(catsim_color)/trim])

        plt.xlim((xmin,xmax))
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

plt.tight_layout()
plt.savefig('kepler_principal_colors_restricted_%s.png' % suffix)
plt.close()
