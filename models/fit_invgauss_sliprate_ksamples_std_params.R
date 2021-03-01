# FOAit inverse Gaussian (BPTAA) distirbution to data for earthquake inter-event times
# when sampling multiple chronologies. Fit each Monte Carlo sample of the earthquake
# chronology individually, and then combine all chains into single MCMC object
# to calculate final statistics

#library(ggsci)
library(R2jags)
library(runjags)
library(bayesplot)
library(ggplot2)
library(coda)
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
throws = cbind(15, 18, 20, 30)
slip_times = cbind(50000, 90000, 130000, 150000)#, 55000, 1000000)
isSlipCensored = (slip_times < 0)
isSlipCensored = as.numeric(isSlipCensored)
print(isSlipCensored)
print(nrow(datalist))
# Name of figure file 
pdf('invgauss_fit_sliprate_std_param.pdf')

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
#    print(isCensored)
    isCensored[1] = TRUE
    isCensored[length(isCensored)] = TRUE
#    print(isCensored)
#    print(censorLimitVec)
    Y = inter_event_m
    Y = Y[1,]
#    print(Y)
    V = throws[1,]
    T = slip_times[1,]
    N <- length(Y)
    M = length(V)
    isCensored = as.numeric(isCensored)
    Y[1] = censorLimitVec[1]
    Y[length(Y)] = censorLimitVec[length(Y)]
    print(Y)
    print('T')
    print(T)
    print(V)
    censorLimitVec = as.numeric(censorLimitVec)
#    print(censorLimitVec)

    # Define data
    sim.data.jags <- list("Y", "N", "V", "T", "M",  "censorLimitVec", "isCensored",
    		  "isSlipCensored")

    # Define the parameters whose posterior distributions we want to calculate
    bayes.mod.params <- c("lambda", "mu", "alpha", "n_events")

    lambdaInit = 1.0
    muInit = 1000 # Rough estimate of mean(Y)

#    print(isCensored)
#    print(censorLimitVec)
#    print(censorLimitVec[isCensored])

    # Define starting values
    bayes.mod.inits <- function(){list("mu"=muInit, "lambda"=lambdaInit)}

    # The model
    bayes.mod.fit <- jags(data = sim.data.jags, inits = bayes.mod.inits,
    	      parameters.to.save = bayes.mod.params, n.chains = 3,
	      n.iter = 12000, n.burnin = 1000, model.file = 'invgauss_sliprate_std_param.jags')
	      
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
bayes.mcmc.combinedlist = as.mcmc.list(mcmclist)
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
# Do separate density plot of all combined
bayes.mcmc.combined = combine.mcmc(mcmclist)
densityplot(bayes.mcmc.combined, layout=c(2,2), aspect="fill")

# Do some 2D plots
posterior <- as.array(bayes.mcmc.combined)
#print(posterior)
dim(posterior)
#color_scheme_set('gray')
mcmc_scatter(posterior, pars = c("mu", "alpha"),
			size = 1.5, alpha = 0.5)
h = mcmc_hex(posterior, pars = c("mu", "alpha"))
#h + plot_bg(fill = "gray95") + panel_bg(fill = "gray70")
#h + stat_binhex(aes(colour = ..density.., fill = ..density..))
h = h + geom_hex(aes_(color = ~ scales::rescale(..density..)))
h + scale_color_gradientn("Density", colors = unlist(color_scheme_get()),
  breaks = c(.1, .9), labels = c("low", "high"))
df_post = do.call(rbind.data.frame, bayes.mcmc.combinedlist)
#print(df_post)
ggplot(df_post, aes(x=mu, y=alpha)) +
	     stat_density_2d(aes(fill = ..level..), geom = "polygon")
ggplot(df_post, aes(x=mu, y=alpha)) +
		stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
		scale_fill_distiller(direction=1) +
		labs(colour = "Density") +
		scale_x_continuous(expand = c(0, 0)) +
		scale_y_continuous(expand = c(0, 0))
#		theme(
#    		legend.position='none'
#  		)
dev.off()
