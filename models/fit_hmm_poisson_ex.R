library(R2jags)
library(lattice)

dat <- read.table("http://www.hmms-for-time-series.de/second/data/earthquakes.txt")

HMM <- function(){
  for(i in 1:m){
    delta[i] <- 1/m
    v[i] <- 1}
  s[1] ~ dcat(delta[])
  for(i in 2:100){
    s[i] ~ dcat(Gamma[s[i-1],])}
  states[1] ~ dcat(Gamma[s[100],])
  x[1]~dpois(lambda[states[1]])
  for(i in 2:n){
    states[i]~dcat(Gamma[states[i-1],])
    x[i]~dpois(lambda[states[i]])}
  for(i in 1:m){
    tau[i]~dgamma(1,0.08)
    Gamma[i,1:m]~ddirch(v[])}
  lambda[1]<-tau[1]
  for(i in 2:m){
    lambda[i]<-lambda[i-1]+tau[i]}}

x = dat[,2]
n = dim(dat)[1]
m = 2

bayes.mod.fit = jags(data = list("x", "n", "m" ), inits = NULL,
    parameters.to.save = c("lambda","Gamma"),
    model.file = HMM, n.iter = 10000, n.chains = 1)

print(bayes.mod.fit)
plot(bayes.mod.fit)
traceplot(bayes.mod.fit)

densityplot(bayes.mod.fit.mcmc, layout=c(2,2), aspect="fill")

#Auto-correlation plot
autocorr.plot(bayes.mod.fit.mcmc)

dev.off()

#summary(glm(sim.dat$y, family=poisson)) 