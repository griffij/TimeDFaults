q# Fit inverse GausOAsian (BPT) distirbution to data for earthquake inter-event times
# when sampling multiple chronologies. In this version pass list of Monte Carlo samples
# of earthquake chronologies and sample from these within jags
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

datafiles = c('../data/chronologies/Hyde_10000_chronologies.csv')

print(datafiles)
for (i in 1:length(datafiles)){
    data = read.csv(datafiles[i], header=FALSE)#, delimiter=',')
    # Convert to ka rather than years to look at sensitivity to priors
    data = data/1000
    # reverse order
    #datalist = data[,order(ncol(data):1)]
    dl = data*-1 # Make ages positive for years before present
    if (i==1){
       datalists = list(dl)
       }else{
       datalists = list.append(datalists, dl)
       }
    }
#print(datalists)

throws = cbind(25.0, 14, 4)# Vertical offsets in meters  
total_throw = throws[1]
V_sigma = cbind(1.5, 1.0, 0.5) # Uncertainty on throw (metres)
total_throw_sigma = V_sigma[1]
V_tau = 1/(V_sigma**2)[1,]
slip_times = cbind(112800, 97800, 34200)
T_sigma = cbind(4400, 8100, 800)
T_tau = 1/(T_sigma**2)[1,]

# Convert to ka rather than years to look at sensitivity to priors
slip_times = slip_times/1000
T_sigma = T_sigma/1000
T_tau = 1/(T_sigma**2)[1,]

total_time = slip_times[1]
total_time_sigma = T_sigma[1]
isSlipCensored = (slip_times < 0)
isSlipCensored = as.numeric(isSlipCensored)
print("isSlipCensored")
print(isSlipCensored)



j=1
for (datalist in datalists){
    #print("datalist")
    #print(datalist)
    # Initialise list for storing each MCMC object
    mcmclist = vector("list", 3)#*length(datalist))
    index = 1
    figure_filename = paste0('plots/invgauss_fit_sliprate_std_param_eqsample_hyde', j, '.pdf')
    pdf(figure_filename)
    # Get number of events
    sl = length(datalist[0,])#[1]

    # Most recent event in year before present defines minimum length of
    # present open interval
    MRE = abs(datalist[,sl])
    censorLimitVec = matrix(1, length(datalist[,1]), length(datalist[1,])-2)
    censorLimitVec[,length(censorLimitVec[1,])] = MRE

    # Convert data to inter-event times
    m = data.matrix(datalist)
    inter_event_m = t(abs(diff(t(m)))) # Transpose, take difference and transpose back
    isCensored = matrix(FALSE, length(datalist[,1]), length(datalist[1,]))
    isCensored[,length(datalist)] = TRUE
    print(isCensored)
    Y = matrix(1, length(datalist[,1]), length(datalist[1,])) 
    for (p in 1:length(datalist[1,])-1){
    	Y[,p] = inter_event_m[,p]
	}
    Y[,length(datalist[1,])] = censorLimitVec[,length(censorLimitVec[1,])]
    print(Y)
    V = throws[1,]
    T = slip_times[1,]
    N <- length(Y[1,])
    print(N) # Number of inter-event times (including censored)
    N_MC = length(Y[,1]) # Number of Monte Carlo samples of eq chronology
    M = length(V)
    isCensored = as.numeric(isCensored)
    print(Y)
    print(isCensored)
    # Lets deal with V and T as uncertain observations
    V_obs = V
    T_obs = T
    Y_obs = Y
    # Define data
    sim.data.jags <- list("Y_obs", "N"#, "N_MC"
    		  ,"V_obs", "V_tau", "T_obs", "T_tau"
	  	  ,"M", "isSlipCensored"
		  , "isCensored", "MRE"
		  )
    # Define the parameters whose posterior distributions we want to calculate
    bayes.mod.params <- c("lambda", "mu", "alpha", "n_events", "n_events_cont",
        "V", "T", "V_sum", "y", "T_sum", "V_obs", "T_obs", "y_ind", "mre"
	)
    alphaInit = 1.0
    #lambdaInit = 1.0
    muInit = 10 #10000 # Rough estimate of mean(Y)
	
    # Define starting values
    bayes.mod.inits <- function(){list("alpha"=alphaInit, "mu"=muInit)}

    # The model
    bayes.mod.fit <- jags(data = sim.data.jags, inits = bayes.mod.inits,
    		  parameters.to.save = bayes.mod.params, n.chains = 3,
		  n.iter = 1000000, n.burnin = 20000, n.thin = 20,
		  model.file = 'invgauss_sliprate_std_param_eqsample_hyde.jags')
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
    print(mcmclist)
    # Convert all MCMC models into one MCMC object
    #	bayes.mcmc.combined = as.mcmc.list(mcmclist) # Use this to plot all individual chains
    bayes.mcmc.combinedlist = as.mcmc.list(mcmclist)
    # Dump data to file
    df_post = do.call(rbind.data.frame, bayes.mcmc.combinedlist)
    filename = paste0('outputs/df_posterior_', j, '_hyde.csv')
    print(filename)
    write.csv(df_post, filename, row.names=FALSE)
    bayes.mcmc.combined = as.mcmc.list(mcmclist) # Use this to plot all individual chains
    summary(bayes.mcmc.combined) 

    # Gelman-Rubin convergence diagnostic
#    gd = gelman.diag(bayes.mod.fit.mcmc, transform=TRUE)
    g <- matrix(NA, nrow=nvar(bayes.mcmc.combined), ncol=2)
    for (v in 1:nvar(bayes.mcmc.combined)) {
    	g[v,] <- gelman.diag(bayes.mcmc.combined[,v])$psrf
	}
    print('Gelman-Rubin diagnostics, values near 1 indicate convergence, NaN may be data variables so check')
    print(g)

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
#	df_post = do.call(rbind.data.frame, bayes.mcmc.combinedlist
	
	print(ggplot(df_post, aes(x=mu, y=alpha)) +
		     stat_density_2d(aes(fill = ..level..), geom = "polygon") +
		     xlab(expression(mu)) +
		     ylab(expression(alpha)))
	
	#print(typeof(df_post$mu))
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


png(file='plots/invgauss_fit_sliprate_std_param_hyde.png', units="in", width=5, height=5, res=300)
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
		scale_x_continuous(expand = c(0, 0), limits = c(0, 150)) +
		scale_y_continuous(expand = c(0, 0), limits = c(0, 10))
#		theme(
#    		legend.position='none'
#  		)

dev.off()
