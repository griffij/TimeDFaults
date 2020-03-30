# Script for running simulations of a Brownian oscillator

source('BPT.R') # Call module containing brownian oscillator function

# Define simulation parameters
t = 8 # Total simulation time
dt = 0.01 # Timestep
x0 = 0 # Value immediatly after failure
xf = 1 # Failure threshold value
mu = 0 # Mean value of normally distributed white noise, set to zero
sigma = c(1/4, 1/2, 3/4, 0.9, 1, 1.1, 1.25) # Standard deviation  
var = sigma^2 # Perturbation rate parameter for Brownian oscillator, variance of
      	   # normal distribution
lambda = 1 # Mean loading rate (i.e. simulates constant tectonic loading)
rseed = 2 # Fix random seed for repeatability

fig_filename = 'brownian_oscillators.pdf'
pdf(fig_filename)
par(mfrow=c(4, 1))
for (i in seq_along(sigma)){
    oscillator = brownian_oscillator(lambda, t, sigma[i], mu, dt,
    			      	     x0, xf, plot=TRUE,
    	      			     healing=FALSE, rseed=rseed)
    plot(oscillator$realisation$n, oscillator$realisation$Y, type = 'l',
    	 main = bquote('Brownian Oscillator,' ~ sigma == .(sigma[i])),
    	 xlab = 'Time',   ylab = 'State (Y(t)), failure at Y=1', ylim=c(-2,2))
	 #asp=0.5)
    for (i in seq_along(oscillator$event_times)){
    	 lines(c(oscillator$event_times[i], oscillator$event_times[i]), c(-10,10),
	       lty=3)
	 }
    }
dev.off()