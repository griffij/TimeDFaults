# Fit inverse Gaussian (BPT) distirbution to data for earthquake inter-event times
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
library(ks)
# Fix random seed
#set.seed(23)

setwd('.')


###########
# Real data

datafile = '../data/chronologies/Dunstan_10_chronologies.csv'
data = read.csv(datafile, header=FALSE)#, delimiter=',')
print(data)
# reverse order
#datalist = data[,order(ncol(data):1)]
datalist = data*-1 # Make ages positive for years before present
print(datalist)
#datalist = data
# Dunstan - Use NA as placeholder for open intervals
#data1 = cbind(NA, 13600, 15100, 16700, 23000, NA)
#data2 = cbind(NA, 13200, 15000, 17500, 24000, NA)
#datalist = list(data1, data2)
throws = cbind(28, 13, 7)#, 3)# Vertical offsets in meters  
#throws = cbind(7, 13, 28)
#V_sigma = cbind(2, 2, 2)
V_sigma = cbind(2, 2, 2)#, 1) # Uncertainty on throw (metres)
V_tau = 1/(V_sigma**2)[1,]
slip_times = cbind(340000, 200000, 100000)#, 31000)
T_sigma = cbind(20000, 10000, 10000)#, 2500)
#slip_times = cbind(100000, 200000, 340000)
#T_sigma = cbind(10000, 10000, 20000)
T_tau = 1/(T_sigma**2)[1,]
isSlipCensored = (slip_times < 0)
isSlipCensored = as.numeric(isSlipCensored)
print(isSlipCensored)
print(nrow(datalist))
# Times for evaluating hazard function
hf_times = seq(from = 1, to = 20000, by = 1000)

# Name of figure file 
pdf('plots/invgauss_fit_sliprate_std_param.pdf')

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
##    censorLimitVec = cbind(abs(data[2]), 1, 1, 1, 50000-abs(data[sl]))
    # Most recent event in year before present defines minimum length of
    # present open interval
    MRE = abs(data[sl])
    print('MRE')
    print(MRE)
    censorLimitVec = cbind(50000-abs(data[2]), 1, 1, 1, MRE)
    # Convert data to inter-event times
    m = data.matrix(data)
    inter_event_m = t(abs(diff(t(m)))) # Transpose, take difference and transpose back
    isCensored = (inter_event_m < 0)
    isCensored[1] = TRUE
    isCensored[length(isCensored)] = TRUE
    Y = inter_event_m
    print("Y")
    print(Y)
    Y = Y[1,]
    print(Y)
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
    print('censorLimitVec')
    print(censorLimitVec)
    print(isCensored)
    print(isSlipCensored)
    # Lets deal with V and T as uncertain observations
    V_obs = V
    T_obs = T
    # Define data
    sim.data.jags <- list("Y", "N"
    		  ,"V_obs", "V_tau", "T_obs", "T_tau"
    		  ,"M", "isSlipCensored"
		  ,"censorLimitVec", "isCensored",
		  "hf_times"
		  )

    # Define the parameters whose posterior distributions we want to calculate
    bayes.mod.params <- c("lambda", "mu", "alpha", "n_events", "hf"
    		     ,"V", "T", "Y", "y[1]", "V_sum", "T_sum", "V_obs", "T_obs","future_event_prob", "MRE"
		     ) 
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
	      n.iter = 10000, n.burnin = 1000, model.file = 'invgauss_sliprate_std_param.jags')
	      
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
	     stat_density_2d(aes(fill = ..level..), geom = "polygon") +
	     xlab(expression(mu)) +
	     ylab(expression(alpha))
dev.off()
png(file='plots/invgauss_fit_sliprate_std_param.png', units="in", width=5, height=5, res=300)
# Estimate density value containing 95% of posterior ditribution
df_param = data.frame(df_post$mu, df_post$alpha)
kd <- ks::kde(df_param, compute.cont=TRUE)
contour_95 =with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]], 
	   z=estimate, levels=cont["5%"])[[1]])
contour_95 = data.frame(contour_95)
contour_65 =with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
	   z=estimate, levels=cont["35%"])[[1]])
contour_65 = data.frame(contour_65)   
ggplot(df_post, aes(x=mu, y=alpha)) +
		stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
		geom_path(aes(x=x, y=y), data=contour_95) +
		geom_path(aes(x=x, y=y), data=contour_65) +
		scale_fill_distiller(direction=1) +
		labs(colour = "Density") +
		xlab(expression(mu)) +
		ylab(expression(alpha)) +
		scale_x_continuous(expand = c(0, 0), limits = c(0, 150000)) +
		scale_y_continuous(expand = c(0, 0), limits = c(0, 10))
#		theme(
#    		legend.position='none'
#  		)
#pp + ggplot(geom_path(aes(x=x, y=y), data=contour_95))
dev.off()
