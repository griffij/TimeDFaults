# Fit inverse Gaussian (BPTAA) distirbution to data for earthquake inter-event times
# when sampling multiple chronologies. Fit each Monte Carlo sample of the earthquake
# chronology individually, and then combine all chains into single MCMC object
# to calculate final statistics

library(R2jags)
library(runjags)
library(lattice)
# Fix random seed
#set.seed(23)

setwd('.')


###########
# Real data

datafile = '../data/chronologies/Dunstan_10_chronologies.csv'
data = read.csv(datafile, header=FALSE)#, delimiter=',')
print(data)
datalist = data
# Dunstan - Use NA as placeholder for open intervals
#data1 = cbind(NA, 13600, 15100, 16700, 23000, NA)
#data2 = cbind(NA, 13200, 15000, 17500, 24000, NA)
#datalist = list(data1, data2)
print(nrow(datalist))
# Name of figure file 
pdf('invgauss_fit_censored_combined.pdf')

# Initialise list for storing each MCMC object
mcmclist = vector("list", 3*length(datalist))
for (i in 1:nrow(datalist)){
    data = list(cbind(NA, datalist[i,], NA))#[1,]
    data = data[[1]]
    print('data')
    print(data)
    print(typeof(data))
    print(data[4])
    print(length(datalist)-1)
    sl = (length(data) - 1)[1]
    print(sl)
    censorLimitVec = cbind(abs(data[2]), 1, 1, 1, 50000-abs(data[sl]))
    # Convert data to inter-event times
    m = data.matrix(data)
    inter_event_m = t(diff(t(m))) # Transpose, take difference and transpose back
    isCensored = (inter_event_m < 0)
    print(isCensored)
    isCensored[1] = TRUE
    isCensored[length(isCensored)] = TRUE
    print(isCensored)
    print(censorLimitVec)
    Y = inter_event_m
    Y = Y[1,]
    print(Y)

    N <- length(Y)
    print(mean(Y))	

    isCensored = as.numeric(isCensored)
    Y[1] = censorLimitVec[1]
    Y[length(Y)] = censorLimitVec[length(Y)]
    print(Y)
    censorLimitVec = as.numeric(censorLimitVec)
    print(censorLimitVec)

    # Define data
    sim.data.jags <- list("Y", "N", "censorLimitVec", "isCensored")

    # Define the parameters whose posterior distributions we want to calculate
    bayes.mod.params <- c("alpha", "mu")

    alphaInit = 1.0
    muInit = 1000 # Rough estimate of mean(Y)

    print(isCensored)
    print(censorLimitVec)
    print(censorLimitVec[isCensored])

    # Define starting values
    bayes.mod.inits <- function(){list("mu"=muInit, "alpha"=alphaInit)}

    # The model
    bayes.mod.fit <- jags(data = sim.data.jags, inits = bayes.mod.inits,
    	      parameters.to.save = bayes.mod.params, n.chains = 3,
	      n.iter = 12000, n.burnin = 1000, model.file = 'invgauss_censored.jags')
	      
    print(bayes.mod.fit)
#    plot(bayes.mod.fit)
#    traceplot(bayes.mod.fit)

    # Convert to an MCMC object
    bayes.mod.fit.mcmc <- as.mcmc(bayes.mod.fit)
    # Density plot   
#    densityplot(bayes.mod.fit.mcmc, layout=c(2,2), aspect="fill")
    # Add to list of all MCMC objects
    print(typeof(bayes.mod.fit.mcmc)) # list
    # Get all three chains
    mcmclist[3*i-2] = bayes.mod.fit.mcmc[1]
    mcmclist[3*i-1] = bayes.mod.fit.mcmc[2]
    mcmclist[3*i] = bayes.mod.fit.mcmc[3] 
#    mcmclist = cbind(mcmclist, bayes.mod.fit.mcmc)
    summary(bayes.mod.fit.mcmc)
#    print(bayes.mod.fit.mcmc)
#    # Some more plots - only plot summary plots now
#    xyplot(bayes.mod.fit.mcmc, layout=c(2,2), aspect="fill")
#
#    #Auto-correlation plot
#    autocorr.plot(bayes.mod.fit.mcmc)
    }
#print(mcmclist)
# Convert all MCMC models into one MCMC object

bayes.mcmc.combined = as.mcmc.list(mcmclist) # Use this to plot all individual chains
#bayes.mcmc.combined = combine.mcmc(mcmclist) # Use this to plot one combined density
# Print summary from combined mcmc objects
#bayes.mcmc.combined = unlist(bayes.mcmc.combined)
#print(typeof(bayes.mcmc.combined))
#print(bayes.mcmc.combined[,][,1])
#print(bayes.mcmc.combined[1])
summary(bayes.mcmc.combined) 
#print(bayes.mcmc.combined)
# Now plot combined results
xyplot(bayes.mcmc.combined, layout=c(2,2), aspect="fill")
# Density plot
densityplot(bayes.mcmc.combined, layout=c(2,2), aspect="fill")
#Auto-correlation plot
#autocorr.plot(bayes.mcmc.combined)

dev.off()
