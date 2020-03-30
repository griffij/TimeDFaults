# Run some dummy data through fit_exp_data_ksamples to test that the model
# is giving sensible results

# Clean-up any old variables in the environment
rm(list = ls())
library(R2jags)
library(lattice)
source('fit_exp_ksamples.R')
# Fix random seed
set.seed(1.2)

setwd('.')

# Name of figure file
pdf('test_exponential_fit.pdf')

# Test data
k=10
n.sim <- 10
lam <- 1./100 # Define rate
mean_time <- 1./lam
for (i in 1:k){
    y_init <- rexp(n=n.sim, rate=lam) # Exponential DV
    if (i==1){
       y = y_init
       }
    else{y=cbind(y, y_init)}
    }   
#y = cbind(y_init, y_init) # if k=2
#y = cbind(y_init, y_init, y_init, y_init, y_init, y_init, y_init, y_init, y_init, y_init) # k=10
y = data.matrix(y)
##print(y)
#print(t(y))
Y = t(y)
#################
N <- length(Y[1,])
print(N)
print(Y[1,1])
print(Y[1,2])
print(mean(Y[1,]))
print(mean(Y[2,]))
print(mean(Y[3,]))
#print(mean(Y[4,]))
print(mean(Y[,]))
###############

bayes.mod.fit.mcmc <- fit_exp_k_samples(Y, N, k, plot_lambda=TRUE)

# Some more plots
xyplot(bayes.mod.fit.mcmc, layout=c(2,2), aspect="fill")

# Density plot
print('Densityplot')
densityplot(bayes.mod.fit.mcmc, layout=c(2,2), aspect="fill")

#Auto-correlation plot
autocorr.plot(bayes.mod.fit.mcmc)

dev.off()