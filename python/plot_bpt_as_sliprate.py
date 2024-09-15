"""Plots tAo illustrate links between brownian oscilator and geological observations
of fault slip.
Jonathan Griffin
"""

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.stats import lognorm, norm
#import random

np.random.seed(7) #7

offset = 2
datafile = './brownian_oscillators_sigma_0.75_0.5.csv'
#datafile = './brownian_oscillators_sigma_0.75_1.csv'
figfilename = datafile.rstrip('csv') + 'png'
data = np.genfromtxt(datafile, delimiter=',', skip_header=1)[:,1]
data = data[np.where(data>0)]
offsets = np.zeros(len(data))
#random_offsets = lognorm(0.7, loc=8.4)  #dlnorm(0.7, 8.4)
r = lognorm.rvs(s=0.5, loc=0.7, size=len(data))

for i in range(len(data)):
    if i == 0:
        offsets[i] = r[i]
    else:
        offsets[i] = offsets[i-1] + r[i] 
# Now reverse to get correct way forward for time
offsets = offsets[::-1]

# Repeat to plot as steps
data = np.repeat(data,2)
#data = np.delete(data, 0)
#data = np.append(data, 0)
offsets = np.repeat(offsets,2)
offsets = np.delete(offsets, 0)
offsets = np.append(offsets, 0)

# Now subsample to make 'observations'
if datafile == './brownian_oscillators_sigma_0.75_1.csv':
    indices = [20, 42, 54, 59, 64, 66, 68]
else:
    indices = [10, 14, 22, 30, 32, 34]
data_sub = data[indices]
offsets_sub = offsets[indices]

# Add errors
xerrors = norm.rvs(scale=0.5, loc=0, size=len(data_sub))
xerrors = 0.1*(offsets_sub)*xerrors
xerrors[0] = -2
yerrors = norm.rvs(scale=1.0, loc=0, size=len(data_sub))
yerrors = 0.1*(data_sub)*yerrors
yerrors[0] = 2
print('xerrors', xerrors)
data_sub = data_sub + 0.5*xerrors
offsets_sub = offsets_sub+yerrors

# Plot 'true' slip rate
plt.plot(data, offsets, c='0.5', linestyle='-', linewidth=0.5)
plt.plot(data_sub, offsets_sub, c='k')
plt.scatter(data_sub[-3:], offsets_sub[-3:], c='red', marker='s')
plt.scatter(data_sub[:-3], offsets_sub[:-3], c='darkorange', marker='s') 
# Add error bars
plt.errorbar(data_sub[-3:], offsets_sub[-3:], yerr=1.5*yerrors[-3:], xerr = 1.5*xerrors[-3:], ecolor='red') #0.5*data_sub*yerrors, xerr=0.3*offsets_sub*xerrors)
plt.errorbar(data_sub[:-3], offsets_sub[:-3], yerr=1.5*yerrors[:-3], xerr = 1.5*xerrors[:-3], ecolor='darkorange')
#plt.errorbar(data_sub, offsets_sub, xerr=2*xerrors)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.xlabel('Age', fontsize=16)
plt.ylabel('Displacement', fontsize=16)
ax = plt.gca()
ax.set_aspect(0.25)
plt.xlim(0, 28)
plt.ylim(0, 50)
#plt.tight_layout()
#figfilename = 'brownian_oscillator_slip_rate.png'
plt.savefig(figfilename)
