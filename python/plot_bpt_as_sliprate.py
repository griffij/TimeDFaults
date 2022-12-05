"""Plots to illustrate links between brownian oscilator and geological observations
of fault slip.
Jonathan Griffin
"""

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.stats import lognorm

offset = 2
datafile = './brownian_oscillators_sigma_0.75_1.csv'
data = np.genfromtxt(datafile, delimiter=',', skip_header=1)[:,1]
print(data)
data = data[np.where(data>0)]
print(data)
offsets = np.zeros(len(data))
random_offsets = lognorm(0.7, loc=8.4)  #dlnorm(0.7, 8.4)
r = lognorm.rvs(s=0.5, loc=0.7, size=len(data))
print('r', r)
#print(len(data))
for i in range(len(data)):
##    print(i)
#    offsets[i] = offset*(i+1)
    if i == 0:
        offsets[i] = r[i]
    else:
        offsets[i] = offsets[i-1] + r[i] 
# Now reverse to get correct way forward for time
#offsets = reversed(offsets)
offsets = offsets[::-1]
print(offsets)

# Repeat to plot as steps
data = np.repeat(data,2)
data = np.delete(data, 0)
data = np.append(data, 0)
offsets = np.repeat(offsets,2)
print(data)
print(offsets)
#print(len(offsets))
# Now subsample to make 'observations'
indices = [1, 29, 45, 55, 59, 65, 67, 69]
data_sub = data[indices]
offsets_sub = offsets[indices]

# Plot 'true' slip rate
plt.plot(data, offsets, c='0.5', linestyle='-', linewidth=0.5)
plt.plot(data_sub, offsets_sub, c='k')
plt.scatter(data_sub, offsets_sub)
plt.xlabel('Age')
plt.ylabel('Vertical displacement')
plt.xlim(0, 35)
plt.ylim(0, 72)
figfilename = 'brownian_oscillator_slip_rate.png'
plt.savefig(figfilename)
