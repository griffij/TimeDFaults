# Fit gamma distirbution to data for earthquake inter-event times
# i.e. we are assuming the earthquakes occur as a Poisson process

library(R2jags)
library(lattice)
# Fix random seed
#set.seed(23)

setwd('.')


###########
# Real data

datafile = '../../dataman/testdata/chronologies1000.csv'
#datafile = '../../data/Akatore4eventBdy_output_10000_chronologies.csv'
#datafile = 'chronologies100.csv'
data = read.csv(datafile, header=FALSE)#, delimiter=',')

k=1000

# Name of figure file
pdf('gamma_fit.pdf')

# Convert data to inter-event times
m = data.matrix(data)
inter_event_m = t(diff(t(m))) # Transpose, take difference and transpose back
Y = inter_event_m
#print(Y)

##############
# Test data
#n.sim <- 100
#lam <- 1./100 # Define as inverse of rate
#mean_time <- 1./lam
#y_init <- rexp(n=n.sim, rate=lam) # Exponential DV
#k = 10 #2 # 10
##y = cbind(y_init, y_init) # if k=2
#y = cbind(y_init, y_init, y_init, y_init, y_init, y_init, y_init, y_init, y_init, y_init) # k=10
#y = data.matrix(y)
##print(y)
#print(t(y))
#Y = t(y)
#################
N <- length(Y[1,])
print(N)
print(Y[1,1])
print(Y[1,2])
print(mean(Y[1,]))
print(mean(Y[2,]))
print(mean(Y[3,]))
print(mean(Y[4,]))
print(mean(Y))
###############

sim.data.jags <- list("Y", "N", "k")

# Define the parameters whose posterior distributions we want to calculate
bayes.mod.params <- c( "alpha", "theta") 

#Define starting values
#bayes.mod.inits <- function(){
#		list("mu"=1/0.01)
#			     }

bayes.mod.fit <- jags(data = sim.data.jags, #inits = bayes.mod.inits,
	      parameters.to.save = bayes.mod.params, n.chains = 3,
	      n.iter = 9000, n.burnin = 1000, model.file = 'gamma_ksamples.jags')

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
