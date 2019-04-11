# Fit exponential distirbution to data for earthquake inter-event times
# i.e. we are assuming the earthquakes occur as a Poisson process

library(R2jags)
library(lattice)
# Fix random seed
set.seed(23)

setwd('.')

datafile = '/Users/PCUser/Documents/PhD/modelling/QuakeRates/dataman/chronologies100.csv'
data = read.csv(datafile, header=FALSE)#, delimiter=',')
#print(data)
k=100

# Convert data to inter-event times
m = data.matrix(data)
inter_event_m = t(diff(t(m))) # Transpose, take difference and transpose back
inter_event_v = as.vector(inter_event_m) # Flatten into single vector of inter_event times
#N = length(inter_event_v)

#print(inter_event_v)
#n.sim <- 1000
#lam <- 1./100 # Define as inverse of rate
#mean_time <- 1./lam
#y <- rexp(n=n.sim, rate=lam) # Exponential DV

#sim.dat <- data.frame(y)#, x)
#Y <- sim.dat$y
#print(Y)
#N <- nrow(sim.dat)

Y = inter_event_v
sim.data.jags <- list("Y", "N", "k")

# Define the parameters whose posterior distributions we want to calculate
bayes.mod.params <- c("mu", "lambda")

#Define starting values
bayes.mod.inits <- function(){
		list("lambda"=1/500.)
			     }

bayes.mod.fit <- jags(data = sim.data.jags, #inits = bayes.mod.inits,
	      parameters.to.save = bayes.mod.params, n.chains = 3,
	      n.iter = 9000, n.burnin = 1000, model.file = 'exp_ksamples.bug')

print(bayes.mod.fit)
plot(bayes.mod.fit)
traceplot(bayes.mod.fit)

# Convert to an MCMC object
bayes.mod.fit.mcmc <- as.mcmc(bayes.mod.fit)
summary(bayes.mod.fit.mcmc)

# Somore more plots
xyplot(bayes.mod.fit.mcmc, layout=c(2,2), aspect="fill")

# Density plot
densityplot(bayes.mod.fit.mcmc, layout=c(2,2), aspect="fill")

#Auto-correlation plot
autocorr.plot(bayes.mod.fit.mcmc)

dev.off()

#summary(glm(sim.dat$y, family=poisson))  