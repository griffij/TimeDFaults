""""Simple calculation of distribution of inter-event times from chronology samples
"""

import numpy as np

datafile = '../data/chronologies/Hyde4event_10000_chronologies.csv'

data = np.genfromtxt(datafile, delimiter=',')
print(data)
ie_times = np.diff(data)
print(ie_times)

print('Mean', np.mean(ie_times))
print('2.5p', np.percentile(ie_times, 2.5))
print('97.5p', np.percentile(ie_times, 97.5))  
