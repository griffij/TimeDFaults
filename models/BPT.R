# Contains functions for modelling Brownian motion, including simulation
# of a Brownian oscillator.

# Jonathan Griffin
# University of Otago, March 2020

brownian_oscillator <- function(lambda, t, sigma, mu=0, dt=0.01, x0=0, xf=1,
    plot=TRUE, healing=FALSE, fsc = 0.8, rseed=NULL){
    
    #' Simulate a Brownian motion oscillator
    #' Based on Matthews et al. 2002. A Brownian Model for Recurrent Earthquakes. BSSA
    #' Note that to reproduce their results in Figure 2 of the paper we must assume that
    #' these authors are giving values of the variance (sigma^2, i.e. the perturbation
    #' rate paramaeter) rather than the standard deviation (sigma) as reported.

    #' The relaxation oscillator is defined by the random variable Y(t) as:
    #' Y(t) = x0 + X(t) - X(t-dt)
    #' X(t) is a random variable defined as:
    #' X(t) = lambda*t + sigma*bro,
    #' where bro is Brownian motion defined by:
    #' N(mu, sigma) where N is the normal distribution and
    #' bro = cumulative_sum(N)*sqrt(dt). The sqrt(dt) term provides
    #' scaling for the timestep.
    #' X(t) increases in time with mean rate lambda, with variations due to
    #' the applied Gaussian noise. If the value of X(t) > xf, where xf is a defined failure
    #' threshold, an event is said to have occured, and X(t) is reset to an initial value
    #' X(t) = x0, where x0 defaults to zero. 

    #' @param lambda = mean loading rate (i.e. simulates constant tectonic loading)
    #' @param t = total simulation time
    #' @param sigma = standard deviation of Gaussian white noise added to oscillator;
    #'     var = sigm^2 is the perturbation rate parameter.
    #' @param mu = mean value of normally distributed white noise, default to zero
    #' @param dt = model timestep. This value must also be used to scale the white noise
    #' @param x0 = oscillator value immediately post-failure
    #' @param xf = failure threshold value
    #' @param plot = TRUE. Do some basic plots
    #' @param healing = FALSE. If true, reudce failure threshold immediately following an
    #'     event to simulate fault healing. Failure threshold then gradually increases back to
    #'     intial value xf to simulate fault healing. This is intended to promote clusters.
    #' @param fsc = failure strenght scale. Scale inital value of failure threshold xf by fsc
    #'     to model weakning fault, i.e. post-failure threshold xff = fsc*xf.
    #' @param seed = random seed, which can be fixed to ensure repeatability of simulations
    #' @returns oscillator = dataframe (t, Y(t), event_times), where event times is a vector
    #'     of t|(Y(t)>xf

    # Set random seed to be fixed
    if (!(is.null(rseed))){
        print('Fixing random seed')
	set.seed(rseed)
	}
    delta=xf-x0

    # Vector of timesteps  
    n=rep(dt, t/dt)
    n = cumsum(n)
    # Simulate Brownian motion
    gauss_n = rnorm(t/dt, mu, sigma)
    bro = cumsum(gauss_n)*sqrt(dt) # Check this scaling for timestep

    # Plot Brownian motion
    if (plot){
        # File for saving figures
	fig_filename = paste('brownian_oscillator', 'sig', sigma, 'lambda', lambda,
		     'xf', xf, 'figures.pdf', sep='_')
  	pdf(fig_filename) 
        title = bquote('Brownian Motion in 1D ,' ~  sigma == .(sigma))
    	plot(n, bro, type = 'l', main = title, xlab = 'Time',
            ylab = 'Displacement')
	}
	
    # Calculate random variable X(t)
    X = lambda*(n) + sigma*bro

    if (plot){
        # Plot X(t)
    	plot(n, X, type = 'l', main = bquote('X(t),' ~ sigma == .(sigma)), xlab = 'Time',
    	    ylab = 'X(t)')
	}
	
    # Calculate state as defined by random variable Y(t)
    xs = x0
    Xt = 0
    xf0 = xf # Original failure value
    xff = fsc*xf0 # Value of failure state immediately after failure, if fault
    	  	  # healing is being modelled.
    hc = 0
    event_times = vector()
#    print(X)
    for (i in seq_along(X)) {
    	if (hc > 0) {
       	    xf = xf0 - log10(hc)*0.05
	    hc = hc - 1
       	    }
    	else { xf = xf0}
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
	    event_times = c(event_times, n[i])
	    if (healing){
	       xf = xff
       	       hc = 50} # Healing counter
	    }
    	else {xs = x0}
        }

   if (plot){
       plot(n, Y, type = 'l', main = bquote('Brownian Oscillator,' ~ sigma == .(sigma)),
       	       xlab = 'Time',	ylab = 'State (Y(t)), failure at Y=1',
	       asp=0.5, ylim=c(-2,2))
       }
   dev.off()

   # Return data frame of timesteps and values of Y(t), and vector of event_times
   realisation = data.frame(n, Y)
   output = list()
   output$realisation = realisation
   output$event_times = event_times
   return(output)
   }