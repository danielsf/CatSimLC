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

catsim_dtype = np.dtype([('sedname', str, 300), ('teff', float), ('feh', float), ('logg', float), ('m', float),
                   ('mr', float), ('norm', float),
                   ('u', float), ('g', float), ('r', float), ('i', float), ('z', float)])

kepler_dtype = np.dtype([('teff', float), ('feh', float), ('logg', float), ('m', float)])


kep_file = 'kep_star_data.txt'
catsim_file = 'catsim_star_data_same_pointing.txt'

kep_data = np.genfromtxt(kep_file, dtype=kepler_dtype)
catsim_data = np.genfromtxt(catsim_file, dtype=catsim_dtype)


t_min = min(kep_data['teff'].min(), catsim_data['teff'].max())
t_max = max(kep_data['teff'].max(), catsim_data['teff'].max())
m_min = min(kep_data['m'].min(), catsim_data['m'].min())
m_max = max(kep_data['m'].max(), catsim_data['m'].max())

t_min = 3000.0
t_max = 10000.0
m_min = -10.0
m_max = 15.0

t_ticks = np.arange(np.round(t_min/1000.0)*1000.0, np.round(t_max/1000.0)*1000.0, 2000.0)
t_labels = ['%d' % tt for tt in t_ticks]

dt = 100
dm = 0.1
dfeh = 0.1
dlogg = 0.1

plt.figsize = (30, 30)
plt.subplot(2,2,1)

plot_color(catsim_data['teff'], catsim_data['m'], dt, dm)
plt.title('CatSim')
plt.xlabel('Teff')
plt.ylabel('Absolute Kepler magnitude')
plt.xlim(t_min, t_max)
plt.xticks(t_ticks, t_labels, fontsize=7)
plt.ylim(m_min, m_max)
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()

plt.subplot(2,2,2)
plt.title('Kepler')
plt.xlabel('Teff')
plt.ylabel('Absolute Kepler magnitude')
plot_color(kep_data['teff'], kep_data['m'], dt, dm)
plt.xlim(t_min, t_max)
plt.xticks(t_ticks, t_labels, fontsize=7)
plt.ylim(m_min, m_max)
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()


plt.subplot(2,2,3)
counts, xbins, ybins = np.histogram2d(catsim_data['teff'], catsim_data['m'], bins=100)
catsim = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
            colors='r', alpha=0.5)

print 'catsim counts shape ',counts.shape

counts,xbins,ybins = np.histogram2d(kep_data['teff'], kep_data['m'], bins=100)
kep = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
            colors='blue', alpha=0.5)

print 'kepler counts shape ',counts.shape

plt.xlabel('Teff')
plt.ylabel('Absolute Kepler magnitude')

m_min=-5
m_max=10

t_min=3000.0
t_max=7000.0
t_ticks = np.arange(np.round(t_min/1000.0)*1000.0, np.round(t_max/1000.0)*1000.0, 1000.0)
t_labels = ['%d' % tt for tt in t_ticks]


plt.xlim(t_min, t_max)
plt.xticks(t_ticks, t_labels, fontsize=7)
plt.ylim(m_min, m_max)

plt.text(t_max-0.1*(t_max-t_min),m_min+0.25*(m_max-m_min),
         'Blue is Kepler; Red is CatSim')

plt.gca().invert_xaxis()
plt.gca().invert_yaxis()


plt.subplot(2,2,4)
valid_dex = np.where(catsim_data['feh']>-900.0)
plt.hist(kep_data['feh'], bins=1000, color='b', zorder=1, edgecolor='b', normed=True)
plt.hist(catsim_data['feh'][valid_dex], bins=1000, color='r', zorder=2, edgecolor='r',
         alpha=0.5, normed=True)

plt.xlabel('FeH', fontsize=10)
plt.tight_layout()
plt.savefig('kepler_hr_diagram_same_pointing.eps')
plt.close()

print 'catsim feh range ',catsim_data['feh'][valid_dex].min(),catsim_data['feh'].max()

# now with 0.05 mag offset

plt.figsize = (30,30)
counts, xbins, ybins = np.histogram2d(catsim_data['teff'], catsim_data['m'], bins=100)
catsim = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
            colors='r', alpha=0.5)

print 'catsim counts shape ',counts.shape

counts,xbins,ybins = np.histogram2d(kep_data['teff'], kep_data['m']-0.05, bins=100)
kep = plt.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],
            colors='blue', alpha=0.5)

print 'kepler counts shape ',counts.shape

plt.xlabel('Teff')
plt.ylabel('Absolute Kepler magnitude')
plt.title('subtract 0.05 from magnitude of Kepler stars')

m_min=-5
m_max=10

t_min=3000.0
t_max=7000.0
t_ticks = np.arange(np.round(t_min/1000.0)*1000.0, np.round(t_max/1000.0)*1000.0, 1000.0)
t_labels = ['%d' % tt for tt in t_ticks]


plt.xlim(t_min, t_max)
plt.xticks(t_ticks, t_labels, fontsize=7)
plt.ylim(m_min, m_max)

plt.text(t_max-0.1*(t_max-t_min),m_min+0.25*(m_max-m_min),
         'Blue is Kepler; Red is CatSim')

plt.gca().invert_xaxis()
plt.gca().invert_yaxis()

plt.savefig('kepler_hr_diagram_same_pointing_offset.eps')
plt.close()

