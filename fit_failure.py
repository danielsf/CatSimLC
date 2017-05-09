import numpy as np

dtype = np.dtype([('n', int), ('delta', float), ('freq', float)])

data = np.genfromtxt('failure_data.txt', dtype=dtype)

lnfreq = np.log(data['freq'])
lndelta = np.log(data['delta'])

sum_x = lndelta.sum()
sum_xx = np.power(lndelta,2).sum()
sum_xy = (lnfreq*lndelta).sum()
sum_y = lnfreq.sum()

bb = sum_y-sum_x*sum_xy/sum_xx
bb = bb/(len(lndelta)-sum_x*sum_x/sum_xx)

mm = (sum_xy - bb*sum_x)/sum_xx

print mm,bb

rng = np.random.RandomState(18352)

x= rng.random_sample(50)*100.0
x=np.sort(x)
y=1.2*x*x-9.34*x+12.3

delta = 0.01
from PressRybicki import extirp_sums

n_t_guess = x.max()/delta
n_t = 2
while n_t<n_t_guess:
    n_t*=2

c_t, s_t, tk, hk = extirp_sums(x,y,delta,n_t)
freq_arr = np.array([k/(n_t*delta) for k in range(n_t)])

cut_off = np.exp(bb)*np.power(delta,mm)

valid = np.where(freq_arr<cut_off)

sin_truth = np.dot([np.sin(2.0*np.pi*nu*x) for nu in freq_arr[valid]],
                   y)


cos_truth = np.dot([np.cos(2.0*np.pi*nu*x) for nu in freq_arr[valid]],
                   y)


print np.abs(c_t[valid]/cos_truth-1.0).max()
print np.nanmax(np.abs(s_t[valid]/sin_truth-1.0))

c_t = c_t[valid]
s_t = s_t[valid]

small = np.where(np.abs(sin_truth)<0.001)
print np.max(np.abs(s_t[small]-sin_truth[small]))
