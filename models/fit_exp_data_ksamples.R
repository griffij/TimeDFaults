# Fit exponential distirbution to data for earthquake inter-event times
# i.e. we are assuming the earthquakes occur as a Poisson process

library(R2jags)
library(lattice)
source('fit_exp_ksamples.R')
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
pdf('exponential_fit.pdf')

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

bayes.mod.fit.mcmc <- fit_exp_k_samples(Y, N, k)

# Some more plots
xyplot(bayes.mod.fit.mcmc, layout=c(2,2), aspect="fill")

# Density plot
print('Densityplot')
densityplot(bayes.mod.fit.mcmc, layout=c(2,2), aspect="fill")

#Auto-correlation plot
autocorr.plot(bayes.mod.fit.mcmc)

dev.off()