"""Plot hazard function against aperiodicity value for BPT distribution
"""

import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import gamma, expon, weibull_min
from BPT import bpt_pdf, bpt_cdf, bpt_hazard_function # In Matthews et al parameterisation

mu = 10000 # Mean inter-event time, years
alphas = np.arange(0.01, 10, 0.01)

# Calculate asymptotic hazard function from Matthews et al (2002 BSSA)
# HFinf = 1/(2*mu*alpha**2)
hfs = 1/(2*mu*alphas**2)

plt.clf()
plt.semilogy(alphas, hfs, c='0.7')
plt.xlabel(r'$\alpha$')
plt.ylabel('Asymptotic hazard rate')
mean_label = r'$\mu$ = %i' % mu
plt.annotate(mean_label, (0.8, 0.9), xycoords = 'axes fraction', fontsize = 12)
plt.savefig('BPT_asymptotic_hazard_function.png')

# Now plot CDF
plt.clf()
times = np.arange(1, 30000, 100)
alpha_vals = [0.5, 1.0, 1.5, 2.0, 3.]
linestyles = ['solid', 'dotted', 'dashed', 'dashdot', (0, (5, 10))]
for i,alpha in enumerate(alpha_vals):
    cdf = bpt_cdf(mu, alpha, times)
    plt.plot(times, cdf, linestyle=linestyles[i], c = '0.3', label=(r'$\alpha$ = ' + str(alpha)))
plt.xlabel('Time (years)')
plt.ylabel('Cumulative density')
plt.legend(fontsize=12, handlelength=3)
plt.annotate(mean_label, (0.4, 0.1), xycoords = 'axes fraction', fontsize = 12) 
plt.savefig('BPT_CDF.png')

# Now plot hazard function
plt.clf()
for i,alpha in enumerate(alpha_vals):
    hf = bpt_hazard_function(mu, alpha, times)
    plt.plot(times, hf, linestyle=linestyles[i], c = '0.3', label=(r'$\alpha$ = ' + str(alpha)))
# Add exponential function for comparision
exp_dist = expon(scale=mu)
hf_expon = exp_dist.pdf(times) / (1 - exp_dist.cdf(times))
plt.plot(times, hf_expon, c = '0.7', linestyle=linestyles[0], label='Exponential')
plt.xlabel('Time since last event (years)')
plt.ylabel('Hazard rate')
plt.annotate(mean_label, (0.4, 0.9), xycoords = 'axes fraction', fontsize = 12)   
plt.legend(fontsize=12, handlelength=3)
plt.savefig('BPT_hazard_function.png')



# Now we do for gamma distribution
# Now we do for gamma distribution
plt.clf()
gamma_alpha = [2.2, 1.0, 0.9]
for i, alpha in enumerate(gamma_alpha):
    gam = gamma(alpha, scale=mu)
    # Get hazard function at time times
    gam_cdf = gam.cdf(times)
    plt.plot(times, gam_cdf, linestyle=linestyles[i], c = '0.3', label=(r'$\alpha$ = ' + str(alpha)))
plt.xlabel('Time since last event (years)')
plt.ylabel('Cumulative density')
plt.annotate(mean_label, (0.4, 0.9), xycoords = 'axes fraction', fontsize = 12)
plt.legend(fontsize=12, handlelength=3)
plt.savefig('gamma_CDF.png')


plt.clf()
for i, alpha in enumerate(gamma_alpha):
    gam = gamma(alpha, scale=mu)
    # Get hazard function at time times
    gam_hf = gam.pdf(times) / (1 - gam.cdf(times))
    plt.plot(times, gam_hf, linestyle=linestyles[i], c = '0.3', label=(r'$\alpha$ = ' + str(alpha)))
plt.xlabel('Time since last event (years)')
plt.ylabel('Hazard rate')
plt.annotate(mean_label, (0.4, 0.9), xycoords = 'axes fraction', fontsize = 12)
plt.legend(fontsize=12, handlelength=3)
plt.savefig('gamma_hazard_function.png')
