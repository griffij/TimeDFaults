import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import gamma, expon, weibull_min, lognorm, norm

x = np.arange(0.01, 10, 0.01)

# Desired mean and standard deviation
#Hyde
#mu_x = 2
#sig_x = 0.5
# Dunstan
#median = 2.5
mu_xs = [2.0, 2.4]
sig_xs = [0.5, 0.75]
#Convert to parameters
#mu = np.log(median)

# Data from each fault to add to plot
hyde_throws_rc = [2.0, 1.8+0.1] # Add a bit to express uncertainty in second obs]
hyde_throws_ge = [True, False] # Flag for whether measurement is a minimum
dunstan_throws_devonshire = [1, 19./4] # min and max

for i, mu_x in enumerate(mu_xs):
    sig_x = sig_xs[i]
    mu = np.log(mu_x**2/np.sqrt(mu_x**2 + sig_x**2))
    sig2 = np.log(1  + sig_x**2/mu_x**2)
    sig = np.sqrt(sig2)
    #sig = np.log(sig_x)
    
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

# Add measured throws
plt.bar(hyde_throws_rc, [1.25, 0.833], width=[0.2, 0.3],
        facecolor='cornflowerblue', edgecolor='mediumblue',
        linewidth=1, alpha=0.5)
dunstan_density_value = 1./(dunstan_throws_devonshire[1] - dunstan_throws_devonshire[0])
plt.bar(np.mean(dunstan_throws_devonshire), dunstan_density_value,
        width=dunstan_throws_devonshire[1] - dunstan_throws_devonshire[0],
        facecolor='sandybrown', edgecolor='darkorange', linewidth=1, alpha=0.5)
plt.xlabel('Single-event vertical displacement (m)')
plt.ylabel('Density')
figname = 'lognorm_%.2f_%.2f.png' % (mu_x, sig_x)
#figname = 'lognorm_median%.2f_%.2f.png' % (median, sig_x)
plt.savefig(figname, dpi=300)
