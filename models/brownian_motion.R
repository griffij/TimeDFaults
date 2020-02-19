# Simple simulator of Brownian motion

# Number of simulations
t = 4 # Total simulation time
dt = 0.01 # Timestep
x0 = 0
xf = 1
delta = xf-x0
mu = 0
sigma = 1.1 # Scale parameter for Brownian motion
sigma_b = sigma 
lambda = 1
# Vector of timesteps  
n=rep(dt, t/dt)
n = cumsum(n)
# Simulate Brownian motion
gauss_n = rnorm(t/dt, mu, sigma_b)
bro = cumsum(gauss_n)*sqrt(dt) # Check this scaling for timestep
# Plot Brownian motion
plot(n, bro, type = 'l', main = 'Brownian Motion in 1D', xlab = 'Time',
	  ylab = 'Displacement')

# Calculate random variable X(t)
X = lambda*(n) + sigma*bro
# Plot X(t)
plot(n, X, type = 'l', main = 'X(t)', xlab = 'Time',                                                                                  
          ylab = 'X(t)')

# Calculate state as defined by random variable Y(t)
xs = x0
Xt = 0
print(X)
for (i in seq_along(X)) {
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
       Xt = X[i]}
    else {xs = x0}
#    if (Y[i] < x0) {
#       Y[i] = x0
#       }
    }
#Y = x0 + X  - X[-1]
plot(n, Y, type = 'l', main = 'Y(t)', xlab = 'Time', ylab = 'Y(t)', asp=0.5, ylim=c(-2,2))
dev.off()