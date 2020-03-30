# Script for running simulations of a Brownian oscillator

source('BPT.R') # Call module containing brownian oscillator function

# Define simulation parameters
t = 4 # Total simulation time
dt = 0.01 # Timestep
x0 = 0 # Value immediatly after failure
xf = 1 # Failure threshold value
mu = 0 # Mean value of normally distributed white noise, set to zero
var = 1/2. # Perturbation rate parameter for Brownian oscillator, variance of
      	   # normal distribution
sigma = sqrt(var) # Standard deviation
lambda = 1 # Mean loading rate (i.e. simulates constant tectonic loading)

oscillator = brownian_oscillator(lambda, t, sigma, mu, dt, x0, xf, plot=TRUE,
	      			  healing=FALSE)