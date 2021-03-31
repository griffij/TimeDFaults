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
for (t in 1:4){
    str1 = paste0("hf_times.", t, ".")
    print(str1)
#    print(df_post[c(str)]) 
#    print(df_post[c("hf_times.1.")])
#    print(df_post$hf_times.1.)
    i=df_post[c(str1)][1,1]
    print("i")
    print(i)
    str2 = paste0("hf.", i, ".")
#    print(str2)
#    print("df_post[c(str2)]")
    print(df_post[c(str2)])
    if(t==1){
    	  xvals = df_post[,c(str1)]
	  yvals = df_post[,c(str2)]
	  }
    else{
	  xvals = c(xvals, df_post[,c(str1)])
	  yvals = c(yvals, df_post[,c(str2)])
	}
#    plot(df_post[,c(str1)], df_post[,c(str2)])
#    ggplot(df_post, aes(df_post[,c(str1)], df_post[,c(str2)])) +
#     		     geom_point()
    }
#plot(xvals, yvals)
df = data.frame(xvals, yvals)
ggplot(df, aes(xvals, yvals)) +
	   geom_point()

# Now calculate hazard function directly from posterior
# Probably more efficient to do this way
hf_times = seq(1, 20000, 1000)
u1_f = matrix(, nrow = length(df_post[,c(str1)]), ncol=length(hf_times))
u2_f = matrix(, nrow = length(df_post[,c(str1)]), ncol=length(hf_times))    
cdf_f = matrix(, nrow = length(df_post[,c(str1)]), ncol=length(hf_times))    
li_f = matrix(, nrow = length(df_post[,c(str1)]), ncol=length(hf_times))    
pdf_f = matrix(, nrow = length(df_post[,c(str1)]), ncol=length(hf_times))    
hf = matrix(, nrow = length(df_post[,c(str1)]), ncol=length(hf_times))    

for (t in length(hf_times)){
    print(t[1])
    u1_f[,t] = (df_post$lambda/(t + df_post$Y.5.)^(1/2)) * ((t+df_post$Y.5.)/df_post$mu - 1)  
    u2_f[,t] = -1*(df_post$lambda/(t + df_post$Y.5.)^(1/2))*((t + df_post$Y.5.)/df_post$mu + 1)  
    cdf_f[,t] = pnorm(u1_f[,t], 0, 1) + exp((2*df_post$lambda)/df_post$mu)*pnorm(u2_f[,t], 0, 1)  
    # Hazard function - first calculate pdf
    li_f[,t] = 0.5*(log(df_post$lambda) - log(2*3.1416) - 3*log(t)) -
    	    0.5*df_post$lambda*((t - df_post$mu)^2)/((df_post$mu^2)*t)
    pdf_f[,t] = exp(li_f[,t])
    hf[,t] = pdf_f[,t] / (1 - cdf_f[,t])
    if (t==1){
       xvals = rep(hf_times[t], length(hf[,t]))
       yvals = hf[,t]
       }else{
       xvals = c(xvals, rep(hf_times[t], length(hf[,t])))
       yvals = c(yvals, hf[,t])
       }
    }
df = data.frame(xvals, yvals)
ggplot(df, aes(xvals, yvals)) +
	   geom_point()

dev.off()