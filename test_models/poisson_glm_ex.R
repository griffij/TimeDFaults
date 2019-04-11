# Generate example Poisson data and run test model to fit
library(R2jags)
library(lattice)

# Fix random seed
set.seed(23)

setwd('.')

n.sim <- 100
beta0 <- 1 # intercept
beta1 <- 1 # slope
x <- rnorm(n=n.sim) # Standard normal predictor
mu <- beta0*1 + beta1*x
lambda <- exp(mu)
y <- rpois(n=n.sim, lambda=lambda) # Poisson DV

sim.dat <- data.frame(y, x)
y <- sim.dat$y
X <- cbind(1,sim.dat$x)
N <- nrow(sim.dat)
mu.beta=rep(0,2)
tau.beta=diag(.0001,2)
sim.data.jags <- list("y", "X", "N", "mu.beta", "tau.beta")

# Define the parameters whose posterior distributions we want to calculate
bayes.mod.params <- c("beta")

#Define starting values
bayes.mod.inits <- function(){
		list("beta"=1)
			     }
bayes.mod.fit <- jags(data = sim.data.jags, #inits = bayes.mod.inits,
	      parameters.to.save = bayes.mod.params, n.chains = 3,
	      n.iter = 9000, n.burnin = 1000, model.file = 'poisson.bug')

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

summary(glm(sim.dat$y~sim.dat$x, family=poisson))  