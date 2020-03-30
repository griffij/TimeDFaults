# Simulate a Brownian motion oscillator
# Based on Matthews et al. 2002. A Brownian Model for Recurrent Earthquakes. BSSA
# Note that to reproduce their results in Figure 2 of the paper we must assume that
# these authors are giving values of the variance (sigma^2, i.e. the perturbation
# rate paramaeter) rather than the standard deviation (sigma) as reported.

# Jonathan Griffin
# University of Otago, March 2020


# Define simulation parameters
t = 4 # Total simulation time
dt = 0.01 # Timestep
x0 = 0 # Value immediatly after failure
xf = 1 # Failure threshold value
delta = xf-x0
mu = 0 # Mean value of normally distributed white noise, set to zero
var = 1/2. # Perturbation rate parameter for Brownian oscillator, variance of
      	   # normal distribution
sigma = sqrt(var) # Standard deviation
lambda = 1 # Mean loading rate (i.e. simulates constant tectonic loading)

# File for saving figures
fig_filename = paste('brownian_oscillator', 'sig', sigma, 'lambda', lambda,
	     'xf', xf, 'figures.pdf', sep='_')
pdf(fig_filename)


# Vector of timesteps  
n=rep(dt, t/dt)
n = cumsum(n)
# Simulate Brownian motion
gauss_n = rnorm(t/dt, mu, sigma)
bro = cumsum(gauss_n)*sqrt(dt) # Check this scaling for timestep
# Plot Brownian motion
title = bquote('Brownian Motion in 1D ,' ~  sigma == .(sigma))
plot(n, bro, type = 'l', main = title, xlab = 'Time',
	  ylab = 'Displacement')

# Calculate random variable X(t)
X = lambda*(n) + sigma*bro
# Plot X(t)
plot(n, X, type = 'l', main = bquote('X(t),' ~ sigma == .(sigma)), xlab = 'Time',                                                              
          ylab = 'X(t)')

# Calculate state as defined by random variable Y(t)
xs = x0
Xt = 0
xf0 = xf # Original failure value
xff = 0.8*xf0 # Value of failure state immediately after failure
hc = 0
print(X)
for (i in seq_along(X)) {
    if (hc > 0) {
       xf = xf0 - log10(hc)*0.05
       hc = hc - 1
       }
    else { xf = xf0}
    print(hc)
    print(xf)
    if (i==1) {
       y = xs + X[i] - Xt
       Y = y
       }
    else {xs
    	 y = xs + X[i] - Xt
	 Y = c(Y,y)
	 }
    if (Y[i] >= xf) {
       Y[i] = xf
       xs = x0
       Xt = X[i]
       xf = xff
       hc = 50} # Healing counter
    else {xs = x0}
#    if (Y[i] < x0) {
#       Y[i] = x0
#       }
    }
#Y = x0 + X  - X[-1]
plot(n, Y, type = 'l', main = bquote('Brownian Oscillator,' ~ sigma == .(sigma)), xlab = 'Time', ylab = 'State (Y(t)), failure at Y=1', asp=0.5, ylim=c(-2,2))
dev.off()