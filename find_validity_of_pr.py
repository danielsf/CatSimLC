from __future__ import with_statement

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from PressRybicki import extirp_sums

import numpy as np

rng = np.random.RandomState(76)

threshold = 0.01
max_iterations = 20

delta_arr = []
freq_failures = []
n_arr = []

for iteration in range(max_iterations):
    time_arr = rng.random_sample(400)*200.0
    time_arr = np.sort(time_arr)
    fn = np.ones(len(time_arr))*rng.random_sample(1)*7.0-3.5

    for ix in range(4):
        aa = rng.random_sample(1)*9.0-4.5
        omega = rng.random_sample(1)*100.0-50.0
        tau = rng.random_sample(1)*100.0-50.0
        fn += aa*np.cos(omega*(time_arr-tau))

        aa = rng.random_sample(1)*9.0-4.5
        omega = rng.random_sample(1)*100.0-50.0
        tau = rng.random_sample(1)*100.0-50.0
        fn += aa*np.sin(omega*(time_arr-tau))

    for log_delta in rng.random_sample(5)*(-2.6)-1.0:
        delta = np.power(10.0, log_delta)
        n_t_guess = time_arr.max()/delta
        n_t = 2
        while n_t<n_t_guess:
            n_t *= 2

        freq_arr = np.array([k/(n_t*delta) for k in range(n_t)])

        cos_test, sin_test, tk, hk = extirp_sums(time_arr, fn, delta, n_t)

        del tk
        del hk

        sin_truth = np.dot(np.array([np.sin(2.0*np.pi*nu*(time_arr))
                                     for nu in freq_arr]),
                           fn)

        offending_dex = np.where(np.logical_or(np.logical_and(np.abs(sin_truth<0.001),
                                                             np.abs(sin_test-sin_truth)>threshold),
                                              np.abs(sin_test/sin_truth-1.0)>threshold))

        del sin_truth

        sin_worst = np.min(offending_dex[0])


        cos_truth = np.dot(np.array([np.cos(2.0*np.pi*nu*(time_arr))
                                     for nu in freq_arr]),
                           fn)


        offending_dex = np.where(np.logical_or(np.logical_and(np.abs(cos_truth<0.001),
                                                             np.abs(cos_test-cos_truth)>threshold),
                                              np.abs(cos_test/cos_truth-1.0)>threshold))

        del cos_truth

        cos_worst = np.min(offending_dex[0])
        if cos_worst<sin_worst:
            worst_dex = cos_worst
        else:
            worst_dex = sin_worst

        freq_worst = freq_arr[worst_dex]
        delta_arr.append(delta)
        freq_failures.append(freq_worst)
        n_arr.append(n_t)

n_arr = np.array(n_arr)
delta_arr = np.array(delta_arr)
freq_failures = np.array(freq_failures)
sorted_dex = np.argsort(delta_arr)
plt.figsize = (30, 30)
plt.subplot(3,1,1)
plt.scatter(delta_arr[sorted_dex], freq_failures[sorted_dex]/delta_arr[sorted_dex])
plt.xlabel('delta')
plt.ylabel('freq/delta')
plt.xscale('log')
plt.yscale('log')

plt.subplot(3,1,2)
plt.scatter(delta_arr[sorted_dex], freq_failures[sorted_dex])
plt.xlabel('delta')
plt.ylabel('freq')
plt.xscale('log')
plt.yscale('log')

plt.subplot(3,1,3)
plt.scatter(n_arr, freq_failures)
plt.xlabel('n_t')
plt.ylabel('freq')
plt.xscale('log')
plt.yscale('log')
plt.savefig('failure_scale.png')
plt.close()

with open('failure_data.txt', 'w') as output_file:
    output_file.write('# n_t delta freq\n')
    for nn, dd, ff in zip(n_arr, delta_arr, freq_failures):
        output_file.write('%d %.6e %.6e\n' % (nn, dd, ff))
