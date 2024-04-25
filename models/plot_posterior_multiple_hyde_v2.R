#Some separate plotting functions
# These ones for combining outputs from multiple MCMC runs

library(bayesplot)
library(ggplot2)
library(coda)
library(lattice)
library(ks)

# Read in posterior dataset
posterior_files = c('outputs/hyde_alpha_norm_1_0.0625_mu_norm_10_0.0004_tpe_lnorm_0.66_16.5/df_posterior_1_hyde.csv')

i=1
for (filename in posterior_files){
    df_post = read.csv(filename)
    print(head(df_post))
    if (i==1){
       mu = df_post$mu
       lambda = df_post$lambda
       y = df_post$mre
       alpha = df_post$alpha      
    }else{
       mu = c(mu, df_post$mu)
       lambda = c(lambda, df_post$lambda)
       y = c(y, df_post$mre)
       alpha = c(alpha, df_post$alpha)     
       }
    i = i+1
    }

#print(mu)
#print(lambda)
#print(y)

figname = 'plots/posterior_mu_hyde.png'  
png(figname, units="in", width=6, height=6, res=300)     
print("Mean mu")
print(mean(mu))
mu_percentiles = quantile(mu, probs = c(0.025, 0.26, 0.5, 0.84, 0.975), na.rm=TRUE)
print("0.025, 0.26, 0.5, 0.84, 0.975")
print(mu_percentiles)
d = density(mu)
plot(d)
dev.off()

figname = 'plots/posterior_alpha_hyde.png'
png(figname, units="in", width=6, height=6, res=300)
print("mean alpha")
print(mean(alpha))
alpha_percentiles = quantile(alpha, probs = c(0.025, 0.26, 0.5, 0.84, 0.975), na.rm=TRUE)
print("0.025, 0.26, 0.5, 0.84, 0.975")
print(alpha_percentiles)
da = density(alpha)
plot(da)
dev.off()

figname = 'plots/posterior_hazard_rate1_hyde.png'
png(figname, units="in", width=6, height=6, res=300)

# Now calculate hazard function directly from posterior
# Probably more efficient to do this way
#seq1 = seq(1, 5000, 20)
#seq2 = seq(5500, 30000, 500)
seq1 = seq(0.01, 0.55, 0.01)
seq2 = seq(0.5500, 3, 0.05)
hf_times = c(seq1, seq2)
#hf_times = hf_times/1e4 # Convert to ka
#print(hf_times)
#hf_times = seq(1, 30000, 100)
u1_f = matrix(, nrow = length(mu), ncol=length(hf_times))
u2_f = matrix(, nrow = length(mu), ncol=length(hf_times))    
cdf_f = matrix(, nrow = length(mu), ncol=length(hf_times))    
li_f = matrix(, nrow = length(mu), ncol=length(hf_times))    
pdf_f = matrix(, nrow = length(mu), ncol=length(hf_times))    
hf = matrix(, nrow = length(mu), ncol=length(hf_times))    

for (t in 1:length(hf_times)){
    u1_f[,t] = ((lambda/(hf_times[t]))^(1/2)) * ((hf_times[t])/mu - 1)  
    u2_f[,t] = -1*((lambda/(hf_times[t]))^(1/2))*((hf_times[t])/mu + 1)  
    cdf_f[,t] = pnorm(u1_f[,t], 0, 1) + exp((2*lambda)/mu)*pnorm(u2_f[,t], 0, 1)  
    # Hazard function - first calculate pdf
    li_f[,t] = 0.5*(log(lambda) - log(2*3.1416) - 3*log(hf_times[t])) -
    	    0.5*lambda*((hf_times[t] - mu)^2)/((mu^2)*hf_times[t])
    pdf_f[,t] = exp(li_f[,t])
    hf[,t] = pdf_f[,t] / (1 - cdf_f[,t])
    #print(hf[,t])
    if (t==1){
       xvals = rep(hf_times[t], length(hf[,t]))
       yvals = hf[,t]
       mean_hf = mean(hf[,t])
       # Get percentiles
       xval_percentiles = rep(hf_times[t], 5)
       yval_percentiles = quantile(hf[,t], probs = c(0.025, 0.26, 0.5, 0.84, 0.975), na.rm=TRUE)
       }else{
       xvals = c(xvals, rep(hf_times[t], length(hf[,t])))
       yvals = c(yvals, hf[,t])
       #print(yvals)
       mean_hf = c(mean_hf, mean(hf[,t]))
       # Get percentiles
       xval_percentiles = cbind(xval_percentiles, rep(hf_times[t], 5))
       yval_percentiles = cbind(yval_percentiles, quantile(hf[,t], probs = c(0.025, 0.26, 0.5, 0.84, 0.975), na.rm=TRUE))
       }
    }

xval_percentiles = do.call(rbind, as.list(xval_percentiles))
yval_percentiles = do.call(rbind, as.list(yval_percentiles))
xval_percentiles = matrix(xval_percentiles, nrow=5)*1e4
yval_percentiles = matrix(yval_percentiles, nrow=5)/1e4
print(yval_percentiles)
df = data.frame(xvals, yvals)
df2 = data.frame(xval_percentiles, yval_percentiles)
linestyles = c(3,2,4,2,3)
for (i in 1:5){
    if (i==1){
       plot(xval_percentiles[i,], yval_percentiles[i,], type='l', lty=linestyles[i],
       				  ylim = c(0, max(yval_percentiles, na.rm=T)), xlim = c(0, max(hf_times)*1e4), 
#       				  ylim = c(0, max(yval_percentiles)), xlim = c(0, max(hf_times)),
				  xlab = 'Time elapsed since most recent event (years)',
				  ylab = 'Hazard rate')
       }else{
    lines(xval_percentiles[i,], yval_percentiles[i,], lty=linestyles[i])
       }
    }
# Add mean curve
lines(xval_percentiles[1,], mean_hf/1e4, lty=1, lwd=2)
# Add legend
legend(17000, 6.1e-4, legend=c('Mean', 'Median', '68% bounds', '95% bounds'), lty=c(1,4,2,3), lwd=c(2,1,1,1))

dev.off()
figname = 'plots/posterior_hazard_rate2_hyde.png'
png(figname, units="in", width=6, height=6, res=300)

# Now redo but calculate hazard function for next 500 years, ie
# taking into account uncertainty on the time of the most recent event,
# and hence length of the current open interval
conditional_time = 500
conditional_step = 10
percentiles = c(0.025, 0.16, 0.5, 0.84, 0.975)
hf_times = seq(0, conditional_time, conditional_step) # These will be added to length of current open interval
hf_times = hf_times/1e4
cdf_f = matrix(, nrow = length(mu), ncol=length(hf_times))    
li_f = matrix(, nrow = length(mu), ncol=length(hf_times))    
pdf_f = matrix(, nrow = length(mu), ncol=length(hf_times))    
hf = matrix(, nrow = length(mu), ncol=length(hf_times))
for (t in 1:length(hf_times)){
    u1_f[,t] = ((lambda/(hf_times[t] + y))^(1/2)) * ((hf_times[t]+y)/mu - 1)  
    u2_f[,t] = -1*((lambda/(hf_times[t] + y))^(1/2))*((hf_times[t] + y)/mu + 1)  
    cdf_f[,t] = pnorm(u1_f[,t], 0, 1) + exp((2*lambda)/mu)*pnorm(u2_f[,t], 0, 1) 
    # Hazard function - first calculate pdf
    li_f[,t] = 0.5*(log(lambda) - log(2*3.1416) - 3*log(hf_times[t]+y)) -
            0.5*lambda*((hf_times[t] + y - mu)^2)/((mu^2)*(hf_times[t]+y))
    pdf_f[,t] = exp(li_f[,t])
    hf[,t] = pdf_f[,t] / (1 - cdf_f[,t])
    if (t==1){
       xvals = rep(hf_times[t], length(hf[,t]))
       yvals = hf[,t]
       mean_hf = mean(hf[,t])
       # Get percentiles
       xval_percentiles = rep(hf_times[t], 5)
       yval_percentiles = quantile(hf[,t], probs = percentiles, na.rm=TRUE)
       }else{
       xvals = c(xvals, rep(hf_times[t], length(hf[,t])))
       yvals = c(yvals, hf[,t])
       mean_hf = c(mean_hf, mean(hf[,t]))
       # Get percentiles
       xval_percentiles = cbind(xval_percentiles, rep(hf_times[t], 5))
       yval_percentiles = cbind(yval_percentiles, quantile(hf[,t], probs = c(0.025, 0.16, 0.5, 0.84, 0.975), na.rm=TRUE))
       }
    }

xval_percentiles = do.call(rbind, as.list(xval_percentiles))
yval_percentiles = do.call(rbind, as.list(yval_percentiles))
xval_percentiles = matrix(xval_percentiles, nrow=5)*1e4
yval_percentiles = matrix(yval_percentiles, nrow=5)/1e4
df = data.frame(xvals, yvals)
df2 = data.frame(xval_percentiles, yval_percentiles)
linestyles = c(3,2,4,2,3)
# Get mean curve
#mean_hf = mean(hf[,t])
print(yval_percentiles)
for (i in 1:5){
    if (i==1){
       plot(xval_percentiles[i,], yval_percentiles[i,], type='l', lty=linestyles[i],
                                  ylim = c(0, max(yval_percentiles ,na.rm=T)), xlim = c(0, max(hf_times)*1e4),
                                  xlab = 'Time elapsed since 2020 (years)',
                                  ylab = 'Hazard rate')
       }else{
    lines(xval_percentiles[i,], yval_percentiles[i,], lty=linestyles[i], lwd=1)
       }
    }
lines(xval_percentiles[1,], mean_hf, lty=1, lwd=2)
# Add legend
legend(300, 1.2e-4, legend=c('Mean', 'Median', '68% bounds', '95% bounds'), lty=c(1,4,2,3), lwd=c(2,1,1,1))   

# Now calculate conditional probability
# Based on Rhoades et al 1994
# Approximate integration by summation
for (i in 1:5){
    hf_int = sum(yval_percentiles[i,])*conditional_step
    conditional_prob = 1 - exp(-1*hf_int)
    cat('Conditional probability', percentiles[i], ' = ', conditional_prob, '\n')
    }
# Now get mean conditional prob
hf_int = sum(mean_hf)*conditional_step
conditional_prob = 1 - exp(-1*hf_int)
cat('Mean conditional probability', ' = ', conditional_prob, '\n')

#ggplot(df2, aes(xval_percentiles, yval_percentiles)) +
#	   geom_point()   
#ggplot(df, aes(xvals, yvals)) +
#	   stat_density_2d(aes(fill = ..level..), geom = "polygon")


dev.off()