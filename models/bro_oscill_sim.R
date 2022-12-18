#OA Script for running simulations of a Brownian oscillator
source('BPT.R') # Call module containing brownian oscillator function

# Define simulation parameters
t = 100#8 # Total simulation time
tplot = 35 # Don't plot whole sequence
dt = 0.001 # Timestep
x0 = 0 # Value immediatly after failure
xf = 1 # Failure threshold value
mu = 0 # Mean value of normally distributed white noise, set to zero
#sigma = c(0.1, 1/4, 1/2, 3/4, 0.9, 1, 1.1, 1.25) # Standard deviation  
sigma = 0.8 #0.7 #0.8
var = sigma^2 # Perturbation rate parameter for Brownian oscillator, variance of
      	   # normal distribution
#lambda = 1 # Mean loading rate (i.e. simulates constant tectonic loading)
#lambda = c(4, 2, 1, 1/2, 0.3, 0.25, 0.2, 0.1)
lambda = c(4, 2,  1, 0.5, 0.4, 0.33, 0.25)
offset = 0 #8 # Offset plot to show good explanatory behaviour
rseed = 8#7ok#2ok#3#3 # Fix random seed for repeatability. Seed of 5 gives good explanatory behaviour

#fig_filename = 'brownian_oscillators.pdf'
fig_filename = paste0('brownian_oscillators_sigma_', sigma, '.pdf') 
pdf(fig_filename, width=8, height=4)
#par(mfrow=c(4, 1), mar=c(1.1,4.2,4.1,1.1))
#dev.new(width=8, height=4, unit="in") 
#for (i in seq_along(sigma)){
for (i in seq_along(lambda)){
#     fig_filename = paste0('brownian_oscillators_', lambda[i], '.pdf')
#     print(pdf(fig_filename))
     oscillator = brownian_oscillator(lambda[i], t, sigma, mu, dt,
    			      	     x0, xf, plot=TRUE,
    	      			     healing=FALSE, rseed=rseed)
    interevent_times = numeric(length(oscillator$event_times))
#    fig_filename = paste0('brownian_oscillators_', lambda[i], '_', sigma, '.pdf')
#    print(pdf(fig_filename))
#    dev.new(width=8, height=4, unit="in")
    print(plot(((-1*(oscillator$realisation$n-offset)+tplot)), oscillator$realisation$Y, type = 'l',
#    	 main = bquote('Brownian Oscillator,' ~ lambda == .(lambda[i])),
    	 xlab =  'Age', ylab = 'State', xlim=c(0,tplot-10), xaxs='i',
	 ylim=c(-1.05,1.15), cex.lab = 1.2, cex.main=1.2))
#	 fig.dim=c(8,1)))
	 #asp=0.5)
    if (i==1){
	 print(title(main = c(bquote('Brownian Oscillator,' ~ sigma == .(sigma[i])),
	 	    bquote(lambda == .(lambda[i])))))
#         title(main = bquote(paste('Brownian Oscillator,' ~ sigma == .(sigma[i]),
#	 	    ~ lambda == .(lambda[i]))))
             }
#    if (i==length(lambda)){
#       print(title(xlab='Time', cex.lab=1.5))
#       }
    for (j in seq_along(oscillator$event_times)){
#         print(lines(c(oscillator$event_times[j]-offset, oscillator$event_times[j]-offset), c(-10,10), lty=3))
	 print(lines(c(-1*(oscillator$event_times[j]-offset) +tplot, -1*(oscillator$event_times[j]-offset) + tplot),
	 						     c(-10,10), lty=3))  
#    	 print(lines(c((-1*oscillator$event_times[j]-offset) + tplot), oscillator$event_times[j]-offset), c(-10,10), lty=3)
	 if (j==1){
	    interevent_times[j] = -1*oscillator$event_times[j]-offset
	    }
	 else{
	    interevent_times[j] = -1*oscillator$event_times[j] - -1*oscillator$event_times[j-1]
	    }  
	 }
    print('interevent_times')	 
    print(interevent_times)
    mean_ie_time = abs(mean(interevent_times))
    print(mean_ie_time)
    std_ie_time = sd(interevent_times)
    print(std_ie_time)
    cov_ie_time = std_ie_time/mean_ie_time
    print(cov_ie_time)
    # Add lambda and COV to plots
    print(mtext(bquote(lambda == .(lambda[i]) ~ ',' ~ COV == .(round(cov_ie_time, 1))), cex=0.8))
#    print(mtext(bquote(lambda == .(lambda[i])), cex=0.8))
#    print(dev.off())
#    dev.off()
    # Write event times to file
    event_times_filename = paste0('brownian_oscillators_sigma_', sigma, '_', lambda[i], '.csv')
    write.csv(-1*(oscillator$event_times-offset)+tplot, event_times_filename)
    }
dev.off()