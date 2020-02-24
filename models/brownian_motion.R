# Simple simulator of Brownian motion

# Number of simulations
t = 4 # Total simulation time
dt = 0.01 # Timestep
x0 = 0
xf = 1
delta = xf-x0
mu = 0
sigma = 1. # Scale parameter for Brownian motion
sigma_b = sigma 
lambda = 1

# File for saving figures
fig_filename = paste('brownian_oscillator', 'sig', sigma, 'lambda', lambda,
	     'xf', xf, 'figures.pdf', sep='_')
pdf(fig_filename)


# Vector of timesteps  
n=rep(dt, t/dt)
n = cumsum(n)
# Simulate Brownian motion
gauss_n = rnorm(t/dt, mu, sigma_b)
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
    else {
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