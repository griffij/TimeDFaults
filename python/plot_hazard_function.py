"""Plot hazard function against aperiodicity value for BPT distribution
"""

import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from BPT import bpt_pdf, bpt_cdf, bpt_hazard_function # In Matthews et al parameterisation

mu = 10000 # Mean inter-event time, years
alphas = np.arange(0.01, 10, 0.01)

# Calculate asymptotic hazard function from Matthews et al (2002 BSSA)
# HFinf = 1/(2*mu*alpha**2)
hfs = 1/(2*mu*alphas**2)

plt.semilogy(alphas, hfs, c='0.7')
plt.xlabel(r'$\alpha$')
plt.ylabel('Asymptotic hazard rate')
label = r'$\mu$ = %i' % mu
plt.annotate(label, (0.8, 0.9), xycoords = 'axes fraction', fontsize = 12)
plt.savefig('BPT_asymptotic_hazard_function.png')

# Now plot CDF 
times = np.arange(1, 30000, 100)
alpha = 1.5
cdf = bpt_cdf(mu, alpha, times)
plt.clf()
plt.plot(times, cdf)
plt.savefig('BPT_CDF.png')

# Now plot hazard function
hf = bpt_hazard_function(mu, alpha, times)
plt.clf()
plt.plot(times, hf)
plt.savefig('BPT_hazard_function.png')
