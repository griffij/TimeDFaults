q# Fit inverse GausOAsian (BPT) distirbution to data for earthquake inter-event times
# when sampling multiple chronologies. Fit each Monte Carlo sample of the earthquake
# chronology individually, and then combine all chains into single MCMC object
# to calculate final statistics

# Clear any previously stored variables
rm(list = ls())

#library(ggsci)
library(R2jags)
library(runjags)
library(bayesplot)
library(ggplot2)
library(rlist)
library(RColorBrewer)
library(coda)
library(lattice)
library(ks)
# Fix random seed
#set.seed(23)

setwd('.')


###########
# Real data

#datafile = '../data/chronologies/Dunstan4event_100_chronologies.csv'
#datafiles = c('../data/chronologies/Dunstan4event_10_chronologies.csv',
#	  '../data/chronologies/Dunstan5event_10_chronologies.csv',
#	  '../data/chronologies/Dunstan6event_10_chronologies.csv')
datafiles = c('../data/chronologies/Dunstan4eventOxcal_100_chronologies.csv',
	  '../data/chronologies/Dunstan5eventOxcal_100_chronologies.csv',
	  '../data/chronologies/Dunstan5eventOxcalv2_100_chronologies.csv',
	  '../data/chronologies/Dunstan6eventOxcal_100_chronologies.csv')
print(datafiles)
for (i in 1:length(datafiles)){
    print(i)
    data = read.csv(datafiles[i], header=FALSE)#, delimiter=',')
    print(data)
    # reverse order
    #datalist = data[,order(ncol(data):1)]
    dl = data*-1 # Make ages positive for years before present
    if (i==1){
       datalists = list(dl)
       }else{
       datalists = list.append(datalists, dl)
       }
    }
print(datalists)
#datalist = data
# Dunstan - Use NA as placeholder for open intervals
#data1 = cbind(NA, 13600, 15100, 16700, 23000, NA)
#data2 = cbind(NA, 13200, 15000, 17500, 24000, NA)
#datalist = list(data1, data2)
throws = cbind(27.5, 17.5, 12.5)#, 3)# Vertical offsets in meters  
#throws = cbind(7, 13, 28)
#V_sigma = cbind(2, 2, 2)
V_sigma = cbind(4, 3, 2)#, 1) # Uncertainty on throw (metres)
V_tau = 1/(V_sigma**2)[1,]
slip_times = cbind(321100, 182400, 92000)#, 31000)
T_sigma = cbind(8600, 13600, 5600)#, 2500)
#slip_times = cbind(100000, 200000, 340000)
#T_sigma = cbind(10000, 10000, 20000)
T_tau = 1/(T_sigma**2)[1,]
isSlipCensored = (slip_times < 0)
isSlipCensored = as.numeric(isSlipCensored)
print(isSlipCensored)
print(ncol(datalists))

j=1
for (datalist in datalists){
    # Initialise list for storing each MCMC object
    mcmclist = vector("list", 3*length(datalist))
    index = 1
    figure_filename = paste0('plots/invgauss_fit_sliprate_std_param_', j, '.pdf')
    pdf(figure_filename)

	for (i in 1:nrow(datalist)){
	    data = list(cbind(NA, datalist[i,], NA))#[1,]
	    data = data[[1]]
	    sl = (length(data) - 1)[1]
	    # Most recent event in year before present defines minimum length of
	    # present open interval
	    MRE = abs(data[sl])
#	    print('MRE')
#	    print(MRE)
	    censorLimitVec = rep(1,length(data)-1)
	    censorLimitVec[1] = 30000-abs(data[2])
	    censorLimitVec[length(censorLimitVec)] = MRE
	    # Convert data to inter-event times
	    m = data.matrix(data)
	    inter_event_m = t(abs(diff(t(m)))) # Transpose, take difference and transpose back
	    isCensored = (inter_event_m < 0)
	    isCensored[1] = TRUE
	    isCensored[length(isCensored)] = TRUE
	    Y = inter_event_m
	    Y = Y[1,]
	    V = throws[1,]
	    T = slip_times[1,]
	    N <- length(Y)
	    M = length(V)
	    isCensored = as.numeric(isCensored)
	    Y[1] = censorLimitVec[1]
	    Y[length(Y)] = censorLimitVec[length(Y)]
	    censorLimitVec = as.numeric(censorLimitVec)
	    # Lets deal with V and T as uncertain observations
	    V_obs = V
	    T_obs = T
	    # Define data
	    sim.data.jags <- list("Y", "N"
	    		  ,"V_obs", "V_tau", "T_obs", "T_tau"
	    		  ,"M", "isSlipCensored"
			  ,"censorLimitVec", "isCensored", "MRE"
			  )
	
	    # Define the parameters whose posterior distributions we want to calculate
	    bayes.mod.params <- c("lambda", "mu", "alpha", "n_events"
	    		     ,"V", "T", "Y", "y[1]", "V_sum", "T_sum", "V_obs", "T_obs", "mre"
			     ) 
	    lambdaInit = 1.0
	    muInit = 1000 # Rough estimate of mean(Y)
	
	    # Define starting values
	    bayes.mod.inits <- function(){list("mu"=muInit, "lambda"=lambdaInit)}
	
	    # The model
	    bayes.mod.fit <- jags(data = sim.data.jags, inits = bayes.mod.inits,
	    	      parameters.to.save = bayes.mod.params, n.chains = 3,
		      n.iter = 6000, n.burnin = 1000, n.thin=20,
		      model.file = 'invgauss_sliprate_std_param.jags')
		      
	    print(bayes.mod.fit)
	
	    # Convert to an MCMC object
	    bayes.mod.fit.mcmc <- as.mcmc(bayes.mod.fit)
	    # Add to list of all MCMC objects
	    # Get all three chains
	    mcmclist[index] = bayes.mod.fit.mcmc[1]
	    mcmclist[index+1] = bayes.mod.fit.mcmc[2]
	    mcmclist[index+2] = bayes.mod.fit.mcmc[3]
	    index = index+3
	    summary(bayes.mod.fit.mcmc)
	    }
	    
	# Convert all MCMC models into one MCMC object
	bayes.mcmc.combined = as.mcmc.list(mcmclist) # Use this to plot all individual chains
	bayes.mcmc.combinedlist = as.mcmc.list(mcmclist)
	summary(bayes.mcmc.combined) 
	
	# Now plot combined results
	print(xyplot(bayes.mcmc.combined, layout=c(2,2), aspect="fill"))
	# Density plot
	print(densityplot(bayes.mcmc.combined, layout=c(2,2), aspect="fill"))
	#Auto-correlation plot
	#autocorr.plot(bayes.mcmc.combined)
	# Do separate density plot of all combined
	bayes.mcmc.combined = combine.mcmc(mcmclist)
	print(densityplot(bayes.mcmc.combined, layout=c(2,2), aspect="fill"))
	
	# Do some 2D plots
	posterior <- as.array(bayes.mcmc.combined)
	dim(posterior)
	print(mcmc_scatter(posterior, pars = c("mu", "alpha"),
				size = 1.5, alpha = 0.5))
	h = mcmc_hex(posterior, pars = c("mu", "alpha"))
	h = h + geom_hex(aes_(color = ~ scales::rescale(..density..)))
	h + scale_color_gradientn("Density", colors = unlist(color_scheme_get()),
	  breaks = c(.1, .9), labels = c("low", "high"))
	df_post = do.call(rbind.data.frame, bayes.mcmc.combinedlist)
	
	print(ggplot(df_post, aes(x=mu, y=alpha)) +
		     stat_density_2d(aes(fill = ..level..), geom = "polygon") +
		     xlab(expression(mu)) +
		     ylab(expression(alpha)))
	
	# Dump data to file
	filename = paste0('outputs/df_posterior_', j, '_dunstan.csv')
	print(filename)
	write.csv(df_post, filename, row.names=FALSE)
	print(typeof(df_post$mu))
	if (j==1){
	   mu_list = df_post$mu
	   alpha_list = df_post$alpha
	}else{
	   mu_list = c(mu_list, df_post$mu)
	   alpha_list[j] = c(alpha_list, df_post$alpha)
	   }
	print("dev.off()")
	dev.off() 
	j = j+1
}

png(file='plots/invgauss_fit_sliprate_std_param.png', units="in", width=5, height=5, res=300)
# Estimate density value containing 95% of posterior ditribution
df_param = data.frame(mu_list, alpha_list)
#df_param = data.frame(df_post$mu, df_post$alpha)
kd <- ks::kde(df_param, gridsize=rep(1401,2) , compute.cont=TRUE)
contour_95 =with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]], 
	   z=estimate, levels=cont["5%"])[[1]])
contour_95 = data.frame(contour_95)
contour_68 =with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
	   z=estimate, levels=cont["32%"])[[1]])
contour_68 = data.frame(contour_68)   
ggplot(df_param, aes(x=mu_list, y=alpha_list)) +
		stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
		geom_path(aes(x=x, y=y), data=contour_95, lwd=0.25) +
		geom_path(aes(x=x, y=y), data=contour_68, lwd=0.5) +
		scale_fill_distiller(palette="Greys", direction=1) +
		labs(colour = "Density") +
		xlab(expression(mu)) +
		ylab(expression(alpha)) +
		scale_x_continuous(expand = c(0, 0), limits = c(0, 150000)) +
		scale_y_continuous(expand = c(0, 0), limits = c(0, 10))
#		theme(
#    		legend.position='none'
#  		)

dev.off()
