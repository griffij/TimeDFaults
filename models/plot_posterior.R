#Some separate plotting functions

library(bayesplot)
library(ggplot2)
library(coda)
library(lattice)
library(ks)

# Read in posterior dataset
df_post = read.csv('outputs/df_posterior.csv')
figname = 'plots/posterior_figs.pdf'
pdf(figname)
head(df_post)
#hf_times = unique(df_post$hf_times)
str1 = "lambda"
#for (t in 1:4){
#    str1 = paste0("hf_times.", t, ".")
#    print(str1)
#    print(df_post[c(str)]) 
#    print(df_post[c("hf_times.1.")])
#    print(df_post$hf_times.1.)
#    i=df_post[c(str1)][1,1]
#    print("i")
#    print(i)
#    str2 = paste0("hf.", i, ".")
##    print(str2)
#    print("df_post[c(str2)]")
#    print(df_post[c(str2)])
#    if(t==1){
#    	  xvals = df_post[,c(str1)]
#	  yvals = df_post[,c(str2)]
#	  }
#    else{
#	  xvals = c(xvals, df_post[,c(str1)])
#	  yvals = c(yvals, df_post[,c(str2)])
#	}
#    plot(df_post[,c(str1)], df_post[,c(str2)])
#    ggplot(df_post, aes(df_post[,c(str1)], df_post[,c(str2)])) +
#     		     geom_point()
#    }
#plot(xvals, yvals)
#df = data.frame(xvals, yvals)
#ggplot(df, aes(xvals, yvals)) +
#	   geom_point()

# Now calculate hazard function directly from posterior
# Probably more efficient to do this way
seq1 = seq(1, 500, 10)
seq2 = seq(1000, 30000, 500)
hf_times = c(seq1, seq2)
print(hf_times)
#hf_times = seq(1, 30000, 100)
u1_f = matrix(, nrow = length(df_post[,c(str1)]), ncol=length(hf_times))
u2_f = matrix(, nrow = length(df_post[,c(str1)]), ncol=length(hf_times))    
cdf_f = matrix(, nrow = length(df_post[,c(str1)]), ncol=length(hf_times))    
li_f = matrix(, nrow = length(df_post[,c(str1)]), ncol=length(hf_times))    
pdf_f = matrix(, nrow = length(df_post[,c(str1)]), ncol=length(hf_times))    
hf = matrix(, nrow = length(df_post[,c(str1)]), ncol=length(hf_times))    

for (t in 1:length(hf_times)){
    print(t[1])
    u1_f[,t] = ((df_post$lambda/(hf_times[t] + df_post$Y.5.))^(1/2)) * ((hf_times[t]+df_post$Y.5.)/df_post$mu - 1)  
    u2_f[,t] = -1*((df_post$lambda/(hf_times[t] + df_post$Y.5.))^(1/2))*((hf_times[t] + df_post$Y.5.)/df_post$mu + 1)  
    cdf_f[,t] = pnorm(u1_f[,t], 0, 1) + exp((2*df_post$lambda)/df_post$mu)*pnorm(u2_f[,t], 0, 1)  
    # Hazard function - first calculate pdf
    li_f[,t] = 0.5*(log(df_post$lambda) - log(2*3.1416) - 3*log(hf_times[t])) -
    	    0.5*df_post$lambda*((hf_times[t] - df_post$mu)^2)/((df_post$mu^2)*hf_times[t])
    pdf_f[,t] = exp(li_f[,t])
    hf[,t] = pdf_f[,t] / (1 - cdf_f[,t])
    if (t==1){
       xvals = rep(hf_times[t], length(hf[,t]))
       yvals = hf[,t]
       # Get percentiles
       xval_percentiles = rep(hf_times[t], 5)
       yval_percentiles = quantile(hf[,t], probs = c(0.025, 0.26, 0.5, 0.84, 0.975))
       print(xval_percentiles)
       print(yval_percentiles)
#       xval_50 = hf_times[t]
#       yval_50 = median(hf[,t])
       }else{
       xvals = c(xvals, rep(hf_times[t], length(hf[,t])))
       yvals = c(yvals, hf[,t])
       # Get percentiles
       xval_percentiles = cbind(xval_percentiles, rep(hf_times[t], 5))
       yval_percentiles = cbind(yval_percentiles, quantile(hf[,t], probs = c(0.025, 0.26, 0.5, 0.84, 0.975)))
       }
    }
xval_percentiles = do.call(rbind, as.list(xval_percentiles))
#xval_percentiles = sapply(xval_percentiles) 
print(xval_percentiles)
yval_percentiles = do.call(rbind, as.list(yval_percentiles))
print(yval_percentiles)
xval_percentiles = matrix(xval_percentiles, nrow=5)
yval_percentiles = matrix(yval_percentiles, nrow=5)
print(xval_percentiles[1,])
print(yval_percentiles[1,]) 
df = data.frame(xvals, yvals)
df2 = data.frame(xval_percentiles, yval_percentiles)
#df2 = as.data.frame(t(df2))
print(df2)
for (i in 1:5){
    if (i==1){
       plot(xval_percentiles[i,], yval_percentiles[i,], type='l',
       				  ylim = c(0, max(yval_percentiles)),
				  xlab = 'Time elapsed since most recent event (years)',
				  ylab = 'Hazard rate')
       }else{
    lines(xval_percentiles[i,], yval_percentiles[i,])
       }
    }
#ggplot(df2, aes(xval_percentiles, yval_percentiles)) +
#	   geom_point()   
#ggplot(df, aes(xvals, yvals)) +
#	   stat_density_2d(aes(fill = ..level..), geom = "polygon")


dev.off()