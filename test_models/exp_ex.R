# Generate example exponentially distributed data and run test model to fit
library(R2jags)
library(lattice)

# Fix random seed
set.seed(23)

setwd('.')

n.sim <- 1000
lam <- 1./100 # Define as inverse of rate
mean_time <- 1./lam
y <- rexp(n=n.sim, rate=lam) # Exponential DV

sim.dat <- data.frame(y)#, x)
Y <- sim.dat$y
print(Y)
N <- nrow(sim.dat)
sim.data.jags <- list("Y", "N")

# Define the parameters whose posterior distributions we want to calculate
bayes.mod.params <- c("mu", "lambda")

#Define starting values
bayes.mod.inits <- function(){
		list("lambda"=1/100.)
			     }

bayes.mod.fit <- jags(data = sim.data.jags, #inits = bayes.mod.inits,
	      parameters.to.save = bayes.mod.params, n.chains = 3,
	      n.iter = 9000, n.burnin = 1000, model.file = 'exp.bug')

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