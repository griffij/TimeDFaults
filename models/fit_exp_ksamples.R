# Function for creating an R2jags model for an exponential distribution
# where we have k samples of that distribution.
# Y = data
# N = number of interevent_times
# k number of samples (e.g. considering k samples of earthquake chronologies
#   from the output of an OxCal model

library(R2jags)
library(lattice)

fit_exp_k_samples <- function(Y, N, k, plot_lambda=FALSE){
	library(lattice)
	library(R2jags)	
	sim.data.jags <- list("Y", "N", "k")
	
	# Define the parameters whose posterior distributions we want to calculate
	if (plot_lambda){
	   bayes.mod.params <- c( "mu", "lambda_k")
	   }
	else{
	   bayes.mod.params <- c( "mu") #Don't plot out many lambdas
	   }

	#Define starting values
	bayes.mod.inits <- function(){
			#list("mu"=1/0.01)
			list("lambda_k[i]"=0.01)
			}

	bayes.mod.fit <- jags(data = sim.data.jags, inits = bayes.mod.inits,
		            parameters.to.save = bayes.mod.params, n.chains = 3,
			    n.iter = 9000, n.burnin = 1000, model.file = 'exp_ksamples.jags')

	print(bayes.mod.fit)
	plot(bayes.mod.fit)
	traceplot(bayes.mod.fit)

	# Convert to an MCMC object
	bayes.mod.fit.mcmc <- as.mcmc(bayes.mod.fit)
	summary(bayes.mod.fit.mcmc)
	return(bayes.mod.fit.mcmc)
	}