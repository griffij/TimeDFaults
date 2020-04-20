# Script for running simulations of a Brownian oscillator
source('BPT.R') # Call module containing brownian oscillator function

# Define simulation parameters
t = 8 # Total simulation time
dt = 0.01 # Timestep
x0 = 0 # Value immediatly after failure
xf = 1 # Failure threshold value
mu = 0 # Mean value of normally distributed white noise, set to zero
#sigma = c(0.1, 1/4, 1/2, 3/4, 0.9, 1, 1.1, 1.25) # Standard deviation  
sigma = 0.75
var = sigma^2 # Perturbation rate parameter for Brownian oscillator, variance of
      	   # normal distribution
#lambda = 1 # Mean loading rate (i.e. simulates constant tectonic loading)
lambda = c(8, 4, 2, 1) #, 1, 1/2)
rseed = 3 # Fix random seed for repeatability

fig_filename = 'brownian_oscillators.pdf'
pdf(fig_filename)
par(mfrow=c(4, 1))
#for (i in seq_along(sigma)){
for (i in seq_along(lambda)){
     oscillator = brownian_oscillator(lambda[i], t, sigma, mu, dt,
    			      	     x0, xf, plot=TRUE,
    	      			     healing=FALSE, rseed=rseed)
    interevent_times = numeric(length(oscillator$event_times))
    plot(oscillator$realisation$n, oscillator$realisation$Y, type = 'l',
    	 main = bquote('Brownian Oscillator,' ~ lambda == .(lambda[i])),
    	 xlab = 'Time',   ylab = 'State (Y(t)), failure at Y=1', ylim=c(-2,2))
	 #asp=0.5)
    for (i in seq_along(oscillator$event_times)){
    	 lines(c(oscillator$event_times[i], oscillator$event_times[i]), c(-10,10),
	       lty=3)
	 if (i==1){
	    interevent_times[i] = oscillator$event_times[i]
	    }
	 else{
	    interevent_times[i] = oscillator$event_times[i] - oscillator$event_times[i-1]
	    }  
	 }
#    print(interevent_times)
    mean_ie_time = mean(interevent_times)
    print(mean_ie_time)
    std_ie_time = sd(interevent_times)
    print(std_ie_time)
    cov_ie_time = std_ie_time/mean_ie_time
    print(cov_ie_time)
    }
dev.off()