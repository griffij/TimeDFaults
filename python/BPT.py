"""Define likelihoods of Brownian Passage Time (Inverse Gaussian)
distribution
"""

import os, sys
import numpy as np
from scipy.stats import invgauss, norm, uniform
import matplotlib
import matplotlib.pyplot as plt

alphas = np.arange(0.1, 5.0, 0.1) # Shape parameter; Adjust later to real alpha def
means = np.arange(1, 2000, 20)
#means = np.arange(1, 30000, 1000)

#data = [-3.414500000000000000e+03,-2.804500000000000000e+03,-2.714500000000000000e+03,-2.034500000000000000e+03,
#        -1.329500000000000000e+03,-7.695000000000000000e+02,-6.945000000000000000e+02,7.855000000000000000e+02,
#        1.535500000000000000e+03]
datafile = '../data/testdata/chronologies100.csv'
n_samples = 100 # Number of Monte Carlo samples of earthquake chronology in datafile
data = np.genfromtxt(datafile, delimiter=',')
#Cadell
#data = [[32400, 38500, 45000, 55000, 62500, 70000, 1000000, 2000000, 4500000]]
# Dunstan
#data = [[-5000, 13600, 15100, 16700, 23000, 50000]]
# Altyn Tagh
#data = [[450,  170,  515,  785, 490, 170, 1530,  925]]
# Define a prior distribution on the mean (e.g. based on slip rate and single event displacements
#prior_mean = norm(loc=10000, scale=5000)
prior_mean = uniform(min(means), max(means)) 
prior_mean = norm(loc=630, scale=1000) 
prior_mean_vals = prior_mean.pdf(means)

prior_alpha = uniform(min(alphas), max(alphas))
prior_alpha_vals = prior_alpha.pdf(alphas)
prior_matrix = np.outer(prior_mean_vals, prior_alpha_vals)
ie_times = np.diff(data)
#ie_times = np.array(data) # If input ie times directly
print('ie_times', ie_times)
print('mean ie time', np.mean(ie_times))

def bpt_pdf(mu, alpha, x_vals):
    """Probability distributino function using
    parameterisation given in Matthews et al (2002) BSSA.
    Used to generate likelihoods
    """
    pdf = np.power((mu/(2*np.pi*(alpha**2)*np.power(x_vals, 3))), 1/2) * \
        np.exp(-1*np.power((x_vals - mu), 2)/(2*mu*x_vals*np.power(alpha, 2)))
#    print(pdf)
    return(pdf)

def bpt_cdf(mu, alpha, x_vals):
    """Evalue cumulative distribution function at x_vals
    given parameters mu and alpha"""
    print('alpha', alpha)
    u1 = (1/alpha) * (np.power(x_vals, 1/2)*np.power(mu, -1/2) - \
                      np.power(x_vals, -1/2)*np.power(mu, 1/2))
    u2 = (1/alpha) * (np.power(x_vals, 1/2)*np.power(mu, -1/2) + \
                      np.power(x_vals, -1/2)*np.power(mu, 1/2))   
    cdf = norm.cdf(u1) + np.exp(2/np.power(alpha,2))*norm.cdf((-1*u2))
    return cdf

def bpt_hazard_function(mu, alpha, x_vals):
    """Evaluate hazard function of BPT distribution at x_vals 
    given mu and alpha"""
    hf = bpt_pdf(mu, alpha, x_vals) / (1 - bpt_cdf(mu, alpha, x_vals))
    return hf

def bpt_likelihoods(means, alphas, ie_times):
    """Function for calculating likelihoods for parameter values of
    BPT distribution given observed data"""
    likelihoods = []
    likelihoods_matt = [] 
    for b in means:
        for a in alphas:
            bpt = invgauss(a, scale=b)
            bpt_l = bpt.pdf(ie_times)
            bpt_matt = bpt_pdf(b, a, ie_times)
            likelihood = np.cumprod(bpt_l)[-1]
            likelihood_matt = np.cumprod(bpt_matt)[-1]   
#            print(b, a, likelihood)
            likelihoods.append(likelihood)
            likelihoods_matt.append(likelihood_matt)    
    likelihoods = np.array(likelihoods)
    likelihoods_matt = np.array(likelihoods_matt) 
    max_llh = np.max(likelihoods)
    max_llh_idx = np.argmax(likelihoods)
    return likelihoods, max_llh, max_llh_idx, likelihoods_matt

# Get matrix of priors

all_llhs = [] # Store all likelihoods
all_llhs_matt = [] 
max_llhs = []
max_llh_idxs = []
mle2s = []
all_posterior = []
for i in range(1):
    print(ie_times[i,:])
    likelihoods, max_llh, max_llh_idx, likelihoods_matt = \
        bpt_likelihoods(means, alphas, ie_times[i,:])
    max_llhs.append(max_llh)
    max_llh_idxs.append(max_llh_idx)
    mle2 = invgauss.fit(ie_times[i,:], loc=0)
    mle2s.append(mle2)
    # Normalise likelihoods (i.e. assume uniform prior distribution over search space)
#    likelihoods = likelihoods/sum(likelihoods)
    # Do Bayesian calculation
    num = likelihoods_matt * prior_matrix.flatten()
    denom = np.sum(num)
    posterior = num / denom
    all_llhs.append(likelihoods)
#    likelihoods_matt = likelihoods_matt/sum(likelihoods_matt)
    all_llhs_matt.append(likelihoods_matt)
    all_posterior.append(posterior)

    
all_llhs = np.array(all_llhs)
print('all_lhs', all_llhs)
all_llhs_matt = np.array(all_llhs_matt)
all_posterior = np.array(all_posterior)
mean_llhs = np.mean(all_llhs, axis=0)
mean_llhs_matt = np.mean(all_llhs_matt, axis=0) 
print(mean_llhs, mean_llhs.size, type(mean_llhs))
#mean_llhs = np.array(mean_llhs)
#likelihoods = np.array(likelihoods) 
#print(max(likelihoods))
xi, yi = np.meshgrid(alphas, means, sparse=False)
print(xi)
print(yi)
#print(likelihoods.reshape(xi.shape))
plt.pcolormesh(yi, xi, mean_llhs.reshape(xi.shape), shading='gouraud',
               cmap=plt.cm.Greys)
plt.xlabel('Mean')
plt.ylabel('Shape parameter')

# Add maximum likelihood estimate(s)
plt.scatter(yi.flatten()[max_llh_idxs], xi.flatten()[max_llh_idxs])
print('Maximum likelihood estimate mu and shape', yi.flatten()[max_llh_idxs], xi.flatten()[max_llh_idxs])
print('mle2s', mle2s)

#plt.pcolormesh(yi, xi, likelihoods.reshape(xi.shape), shading='gouraud',
#               cmap=plt.cm.Greys)
#plt.pcolormesh(xi.flatten(), yi.flatten(), likelihoods, shading='gouraud',
#               cmap=plt.cm.Greys)
plt.savefig('BPT_likelihoods.png')

# Plot Matthews version
plt.clf()
plt.pcolormesh(yi, xi, mean_llhs_matt.reshape(xi.shape), shading='gouraud',
               cmap=plt.cm.Greys)
plt.xlabel('Mean')
plt.ylabel('Alpha')
plt.savefig('BPT_likelihoods_matthews.png')

# Now plot posterior
plt.clf()
# Update later with mean of posteriors
plt.pcolormesh(yi, xi, all_posterior.reshape(xi.shape), shading='gouraud',
               cmap=plt.cm.Greys)
# Add MLE estimate - Need to convert to same parameterisation
#alpha_mle, beta_mle, mu_mle = mle2s[0]
#plt.scatter(mu_mle, alpha_mle)
plt.xlabel('Mean')
plt.ylabel('Alpha')
plt.savefig('BPT_posterior_matthews.png')
# Now plot marginal distributions
#print(len(means))
#print(len(alphas))
llhs_x = mean_llhs.reshape(xi.shape)
#print(llhs_x, llhs_x.shape, xi.shape)
llhs_x = np.sum(llhs_x, axis=1)
#print(llhs_x, llhs_x.shape)

llhs_y = mean_llhs.reshape(xi.shape)
#print(llhs_y, llhs_y.shape, xi.shape)
llhs_y = np.sum(llhs_y.T, axis=1)
#print(llhs_y, llhs_y.shape) 

plt.clf()
plt.plot(means, llhs_x)
plt.savefig('BPT_mean_margninal_dist.png')

plt.clf()
plt.plot(alphas, llhs_y)
plt.savefig('BPT_shape_margninal_dist.png')

# Plot distribution with mle estimates
plt.clf()
bpt_plot = invgauss(xi.flatten()[max_llh_idxs], scale = yi.flatten()[max_llh_idxs])
print('mean bpt', bpt_plot.mean())
xvals = np.arange(0, 2*max(means))
plt.plot(xvals, bpt_plot.pdf(xvals), c='b')
# With second estimate
a,l,s = mle2s[0]
bpt_plot = invgauss(a, scale=s)
plt.plot(xvals, bpt_plot.pdf(xvals), c='r')
print('mean bpt', bpt_plot.mean()) 
plt.savefig('BPT_pdf.png')
