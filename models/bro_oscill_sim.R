# Script for running simulations of a Brownian oscillator
source('BPT.R') # Call module containing brownian oscillator function

# Define simulation parameters
t = 100 # Total simulation time
dt = 0.01 # Timestep
x0 = 0 # Value immediatly after failure
xf = 1 # Failure threshold value
mu = 0 # Mean value of normally distributed white noise, set to zero
#sigma = c(0.1, 1/4, 1/2, 3/4, 0.9, 1, 1.1, 1.25) # Standard deviation  
sigma = 0.8
var = sigma^2 # Perturbation rate parameter for Brownian oscillator, variance of
      	   # normal distribution
#lambda = 1 # Mean loading rate (i.e. simulates constant tectonic loading)
#lambda = c(4, 2, 1, 1/2, 0.3, 0.25, 0.2, 0.1)
lambda = c(4, 1, 0.25)
offset = 0 #8 # Offset plot to show good explanatory behaviour
rseed = 5#3 # Fix random seed for repeatability. Seed of 5 gives good explanatory behaviour

fig_filename = 'brownian_oscillators.pdf'
pdf(fig_filename)
par(mfrow=c(4, 1), mar=c(1.1,4.2,4.1,1.1))
#for (i in seq_along(sigma)){
for (i in seq_along(lambda)){
     oscillator = brownian_oscillator(lambda[i], t, sigma, mu, dt,
    			      	     x0, xf, plot=FALSE,
    	      			     healing=FALSE, rseed=rseed)
    interevent_times = numeric(length(oscillator$event_times))
    plot((oscillator$realisation$n-offset), oscillator$realisation$Y, type = 'l',
#    	 main = bquote('Brownian Oscillator,' ~ lambda == .(lambda[i])),
    	 xlab =  '', ylab = 'State', xlim=c(0,8), xaxs='i',
	 ylim=c(-1.15,1.15), cex.lab = 1.2, cex.main=1.2)
	 #asp=0.5)
    if (i==1){
	 title(main = c(bquote('Brownian Oscillator,' ~ sigma == .(sigma[i])),
	 	    bquote(lambda == .(lambda[i]))))
#         title(main = bquote(paste('Brownian Oscillator,' ~ sigma == .(sigma[i]),
#	 	    ~ lambda == .(lambda[i]))))
             }
    if (i==length(lambda)){
       title(xlab='Time', cex.lab=1.5)
       }
    for (j in seq_along(oscillator$event_times)){
    	 lines(c(oscillator$event_times[j]-offset, oscillator$event_times[j]-offset), c(-10,10),
	       lty=3)
	 if (j==1){
	    interevent_times[j] = oscillator$event_times[j]
	    }
	 else{
	    interevent_times[j] = oscillator$event_times[j] - oscillator$event_times[j-1]
	    }  
	 }
#    print(interevent_times)
    mean_ie_time = mean(interevent_times)
    print(mean_ie_time)
    std_ie_time = sd(interevent_times)
    print(std_ie_time)
    cov_ie_time = std_ie_time/mean_ie_time
    print(cov_ie_time)
    # Add lambda and COV to plots
    mtext(bquote(lambda == .(lambda[i]) ~ ',' ~ COV == .(round(cov_ie_time, 2))), cex=0.8)
    }
dev.off()