"""Plot a histogram of the posterior distribution of the number of events for each terrace offset

Jonathan Griffin, Geoscience Australia
"""

import numpy as np
import matplotlib
from matplotlib import pyplot as plt

#datafiles = ['outputs/hyde_alpha_norm_1_0.0625_mu_norm_10_0.0004_tpe_lnorm_0.66_16.5/df_posterior_1_hyde.csv']
datafiles = ['outputs/dunstan_alpha_norm_1_0.0625_mu_norm_10_0.0004_tpe_lnorm_0.58_4.48/df_posterior_1_dunstan.csv',
              'outputs/dunstan_alpha_norm_1_0.0625_mu_norm_10_0.0004_tpe_lnorm_0.58_4.48/df_posterior_2_dunstan.csv',
              'outputs/dunstan_alpha_norm_1_0.0625_mu_norm_10_0.0004_tpe_lnorm_0.58_4.48/df_posterior_3_dunstan.csv',
              'outputs/dunstan_alpha_norm_1_0.0625_mu_norm_10_0.0004_tpe_lnorm_0.58_4.48/df_posterior_4_dunstan.csv']
name = datafiles[0].split('/')[1].split('_')[0]
print('name', name)
if name == 'hyde':
    max_val = 12.5
else:
    max_val = 35.5
column_names = ['"n_events[1]"', '"n_events[2]"', '"n_events[3]"','"n_events[4]"']

#Get header
f_in = open(datafiles[0], 'r')
header = f_in.readline().split(',')
print(header)


for i, datafile in enumerate(datafiles):
    data_tmp = np.genfromtxt(datafile, skip_header=1, delimiter = ',')
    if i == 0:
        data = data_tmp
    else:
        data = np.append(data, data_tmp, axis=1)

for i, column_name in enumerate(column_names):
    plt.clf()
    try:
        index  = header.index(column_name)
    except ValueError:
        continue
    ns = data[:,index]
    plt.hist(ns, bins=np.arange(-0.5, max_val, 1), density=True,
             facecolor='0.75', edgecolor='0.3')
    plt.xlabel('$n_i$', fontsize=16)
    plt.ylabel('Density', fontsize=16)
    figname = 'plots/n_hist_%s_%i.png' % (name, i)
    plt.savefig(figname)
