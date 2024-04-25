import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import gamma, expon, weibull_min, lognorm, norm

x = np.arange(0.01, 10, 0.1)

# Desired mean and standard deviation
#Hyde
#mu_x = 2
#sig_x = 0.5
# Dunstan
mu_x = 2
sig_x = 1.5
#Convert to parameters
mu = np.log(mu_x**2/np.sqrt(mu_x**2 + sig_x**2))
sig2 = np.log(1  + sig_x**2/mu_x**2)
sig = np.sqrt(sig2)
            
print('mu', mu)
print('sigma', sig)
print('1/sigma', 1/sig)
tau = 1/(sig2)
print('tau', tau)
#lnorm = lognorm(s=1, scale=np.exp(0.1), loc=0.5)
lnorm = lognorm(s=sig, scale=np.exp(mu))#, loc=0)#, loc=0.5)
y =lnorm.pdf(x)
norm2 = norm()#loc=2, scale=1)
y2 = np.exp(norm2.pdf(x))
plt.plot(x,y)
#plt.plot(x, y2, c='r')
figname = 'lognorm_%.2f_%.2f.png' % (mu_x, sig_x)
plt.savefig(figname)
