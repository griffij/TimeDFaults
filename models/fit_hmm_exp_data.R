library(R2jags)
library(lattice)
gc()
set.seed(23)

# Test data                                                                    
#n.sim <- 100                                                                   
#lam <- 1./100 # Define as inverse of rate                                      
#lam2 <-1./800
#mean_time <- 1./lam                                                            
#y_init <- rexp(n=n.sim, rate=lam) # Exponential DV 
#print(y_init)
#y_init2 <- rexp(n=n.sim, rate=lam2)
#print(y_init2)
#Y = cbind(y_init, y_init2)
#Y = data.matrix(Y)
#Y = as.vector(Y)
#print(Y)
#N = n.sim*2

#dates = c(-125000, -1000, 850, 1150) # Akatore
dates <- c(-4500000, -2000000, -1000000, -70000, -62500, -55000, -45000, -38500, -32000) # Cadell
#dates  <- rev(c(-15000, 15100, 17000,19000, 21000, 22500, 24000, 55000))*-1 # Dunstan 6 events. add dummy date for 1500 years in future
#dates <- c(1314, 1350, 1388, 1569, 1597, 1613, 1631, 1658, 1703, 1797, 1833, 2007, 2010) #Sumatra (Mentawai Segment) from Philibosian et al 2017, Figures 20)
#dates <- c(365, 526, 530, 601, 692, 733, 764, 794, 824, 840, 1208, 1237, 1314, 1350, 1388, 1569, 1597, 1613, 1631, 1658, 1703, 1797, 1833, 2007, 2010) #Sumatra (Mentawai Segment) from Philibosian et al 2017, Figures 19,20)
interevent_times = diff(dates)
#interevent_times = (dates[1] - dates)*(-1) 
print(interevent_times)
Y = interevent_times
N = length(interevent_times)

HMM <- function(){
  for(i in 1:m){
    delta[i] <- 1/m
    v[i] <- 1}
  s[1] ~ dcat(delta[])
  for(i in 2:100){
    s[i] ~ dcat(Gamma[s[i-1],])}
  states[1] ~ dcat(Gamma[s[100],])
  x[1]~dexp(lambda[states[1]])
  for(i in 2:n){
    states[i]~dcat(Gamma[states[i-1],])
    x[i]~dexp(lambda[states[i]])}
  for(i in 1:m){
#    tau[i]~dgamma(1,0.08)
#     tau[i]~dgamma(0.01,0.01) # Changing this makes a big difference, need to check
      tau[i]~dunif(0,0.01)
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
print(x)
print(n)
bayes.mod.fit = jags(data = list("x", "n", "m" ), inits = NULL,
    parameters.to.save = c("lambda","Gamma", "mu", "states"),
    model.file = HMM, n.iter = 10000, n.burnin=2000, n.chains = 1)

print(bayes.mod.fit)
plot(bayes.mod.fit)
traceplot(bayes.mod.fit)

# Convert to an MCMC object                                                    
bayes.mod.fit.mcmc <- as.mcmc(bayes.mod.fit)                                   
summary(bayes.mod.fit.mcmc) 

densityplot(bayes.mod.fit.mcmc, layout=c(2,2), aspect="fill")

#Auto-correlation plot
autocorr.plot(bayes.mod.fit.mcmc)

# Calculate present day rate
#print(Gamma)


dev.off()

#summary(glm(sim.dat$y, family=poisson)) 