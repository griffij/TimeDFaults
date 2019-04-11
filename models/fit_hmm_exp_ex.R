library(R2jags)
library(lattice)


# Test data                                                                    
n.sim <- 100                                                                   
lam <- 1./100 # Define as inverse of rate                                      
lam2 <-1./800
mean_time <- 1./lam                                                            
y_init <- rexp(n=n.sim, rate=lam) # Exponential DV 
print(y_init)
y_init2 <- rexp(n=n.sim, rate=lam2)
print(y_init2)
Y = cbind(y_init, y_init2)
Y = data.matrix(Y)
Y = as.vector(Y)
print(Y)
N = n.sim*2

HMM <- function(){
  for(i in 1:m){
    delta[i] <- 1/m
    v[i] <- 1}
  s[1] ~ dcat(delta[])
  for(i in 2:200){
    s[i] ~ dcat(Gamma[s[i-1],])}
  states[1] ~ dcat(Gamma[s[200],])
  x[1]~dexp(lambda[states[1]])
  for(i in 2:n){
    states[i]~dcat(Gamma[states[i-1],])
    x[i]~dexp(lambda[states[i]])}
  for(i in 1:m){
    tau[i]~dgamma(1,0.08)
    Gamma[i,1:m]~ddirch(v[])}
  lambda[1]<-tau[1]
  mu[1] <- 1/lambda[1]
  for(i in 2:m){
    lambda[i]<-lambda[i-1]+tau[i]
    mu[i] <- 1/lambda[i]}}
    
#x = dat[,2]
#n = dim(dat)[1]
m = 2
x=Y
n=N
bayes.mod.fit = jags(data = list("x", "n", "m" ), inits = NULL,
    parameters.to.save = c("lambda","Gamma", "mu"),
    model.file = HMM, n.iter = 9000, n.burnin=1000, n.chains = 1)

print(bayes.mod.fit)
plot(bayes.mod.fit)
traceplot(bayes.mod.fit)

densityplot(bayes.mod.fit.mcmc, layout=c(2,2), aspect="fill")

#Auto-correlation plot
autocorr.plot(bayes.mod.fit.mcmc)

dev.off()

#summary(glm(sim.dat$y, family=poisson)) 