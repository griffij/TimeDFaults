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
plt.semilogy(alphas, hfs, c='0.5')
# Add line showing crossing point
yvals = [1e-4, 1e-4, 1e-7]
xvals = [0, 1/np.sqrt(2), 1/np.sqrt(2)]
extraticks = [1/np.sqrt(2)]
extraticklabels = [r'$\frac{1}{\sqrt{2}}$']
#print(plt.xticks())
#print(plt.xticks()[1])
#print(plt.xticks()[1] + extraticklabels)
ax = plt.gca()
#ax.set_xticks(list(plt.xticks()[0]) + extraticks)
#ax.set_xticklabels(plt.xticks()[1])# + extraticklabels)
plt.xticks(list(plt.xticks()[0]) + extraticks, list(plt.xticks()[0]) + extraticklabels)
plt.semilogy(xvals, yvals, linestyle='dashed', c='0.7')
plt.ylim(1e-7, 1)
plt.xlim(0, 10)
plt.xlabel(r'$\alpha$')
plt.ylabel('Asymptotic hazard rate')
mean_label = r'$\mu$ = %i' % mu
plt.annotate(mean_label, (0.8, 0.9), xycoords = 'axes fraction', fontsize = 12)
plt.tight_layout()
plt.savefig('BPT_asymptotic_hazard_function.png')

# Now plot CDF
plt.clf()
times = np.arange(1, 50000, 100)
alpha_vals = [4.3, 3.9, 2.3] # Dunstan
mu_vals = [33.1, 26.3, 11.2] # Dunstan
fault_name = 'Dunstan'
#alpha_vals = [1.3, 0.8, 0.3]
#mu_vals = [19.1, 13.3, 11.8]
#fault_name = 'Hyde' 
linestyles = ['solid', 'dotted', 'dashed', 'dashdot', (0, (5, 10))]
labels = ['Mean', 'Median', 'Mode']
for i,alpha in enumerate(alpha_vals):
    cdf = bpt_cdf(mu_vals[i]*1000, alpha, times)
    plt.plot(times/1000, cdf, linestyle=linestyles[i], c = '0.3',
             label=(labels[i] + ': ' + r'$\alpha$ = ' + str(alpha) + ', $\mu$ = ' + str(mu_vals[i])))
plt.xlabel('Inter-event time (kyr)', fontsize=12)
plt.ylabel('Cumulative density', fontsize=12)
plt.xlim(0,50)
plt.ylim(0,1)
plt.legend(fontsize=12, handlelength=3)
#plt.annotate(mean_label, (0.4, 0.1), xycoords = 'axes fraction', fontsize = 12) 
plt.savefig('BPT_CDF_%s.png' % fault_name)
print('Mode prob < 2000 years', bpt_cdf(mu_vals[2]*1000, alpha_vals[2], 2000))
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

# Now plot BPT pdf
plt.clf()
for i,alpha in enumerate(alpha_vals):
    pdf = bpt_pdf(mu_vals[i]*1000, alpha, times)
    plt.plot(times/1000, pdf*1000, linestyle=linestyles[i], c = '0.3',
#             linewidth=3,
             label=(labels[i] + ': ' + r'$\alpha$ = ' + str(alpha) + ', $\mu$ = ' + str(mu_vals[i])))
# Add exponential function for comparision
#exp_dist = expon(scale=mu)
#pdf_expon = exp_dist.pdf(times)
#plt.plot(times, pdf_expon, c = '0.7', linestyle=linestyles[0], label='Exponential')
plt.xlabel('Inter-event time (kyr)', fontsize=12)
plt.ylabel('Density', fontsize=12)
#plt.annotate(mean_label, (0.4, 0.1), xycoords = 'axes fraction', fontsize = 12) 
plt.legend(fontsize=12, handlelength=3)
#ax = plt.gca()
#ax.axes.xaxis.set_visible(False)
#ax.axes.yaxis.set_visible(False)
plt.xlim(0,50)
plt.ylim(0,)
plt.tight_layout()
plt.savefig('BPT_pdf_%s.png' % fault_name)

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
