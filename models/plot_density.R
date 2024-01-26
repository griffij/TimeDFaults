# Plot density of multiple runs on one figures

library(ggplot2)
library(ks)
library(lattice)
library(grid)
library(gridExtra)
library(hdrcde)

# Read in posterior dataset

#posterior_files = c('outputs/df_posterior_1_eq_only_dunstan.csv',
#              'outputs/df_posterior_2_eq_only_dunstan.csv',
#              'outputs/df_posterior_3_eq_only_dunstan.csv',
#               'outputs/df_posterior_4_eq_only_dunstan.csv')
#posterior_files = c('outputs/dunstan_alpha_unif_0_10_mu_unif_0_150_tpe_lnorm_0.2_0.8/df_posterior_1_dunstan.csv',
#              'outputs/dunstan_alpha_unif_0_10_mu_unif_0_150_tpe_lnorm_0.2_0.8/df_posterior_2_dunstan.csv',
#              'outputs/dunstan_alpha_unif_0_10_mu_unif_0_150_tpe_lnorm_0.2_0.8/df_posterior_3_dunstan.csv',
#              'outputs/dunstan_alpha_unif_0_10_mu_unif_0_150_tpe_lnorm_0.2_0.8/df_posterior_4_dunstan.csv')
posterior_files = c('outputs/dunstan_alpha_norm_1_0.005_mu_norm_0_0.0001_tpe_lnorm_0.2_0.8/df_posterior_1_dunstan.csv',
		'outputs/dunstan_alpha_norm_1_0.005_mu_norm_0_0.0001_tpe_lnorm_0.2_0.8/df_posterior_2_dunstan.csv',
		'outputs/dunstan_alpha_norm_1_0.005_mu_norm_0_0.0001_tpe_lnorm_0.2_0.8/df_posterior_3_dunstan.csv',
		'outputs/dunstan_alpha_norm_1_0.005_mu_norm_0_0.0001_tpe_lnorm_0.2_0.8/df_posterior_4_dunstan.csv')
#MRE_position = c(5,6,6,7)             

plot_posterior_2d <-function(mu, alpha, fig_lab, lab_x=15, lab_y=9.5){
    # Estimate density value containing 95% of posterior ditribution
    # Get mean values
    mean_mu = mean(mu)
    mean_alpha = mean(alpha)
    # Get median
    median_mu = median(mu)
    median_alpha = median(alpha)
    # Get mode
    mode_vals = hdr.2d(mu, alpha, prob=c(1),
           den=NULL, kde.package = c("ks"),
           h = NULL)
#    print(mode_vals)
    mode = mode_vals$mode
    print('Mode')
    print(mode)
    print('Median')
    cat(median_mu, median_alpha, '\n')
    print('Mean')
    cat(mean_mu, mean_alpha)
    # Here let's add the prior on mu
    xvals = seq(0,150,by=1)
    mu_prior = dnorm(xvals, 10, sqrt(10000))*333
#    mu_prior = dnorm(xvals, 0, sqrt(1000))*33
#    mu_prior = dunif(xvals, 0, 150)*333
    # And on alpha
    yvals = seq(0,10,by=0.01)
    #alpha_prior = dnorm(yvals, 1, sqrt(10))*333
    alpha_prior = dnorm(yvals, 1, sqrt(20))*333
#    alpha_prior = dunif(yvals, 0, 10)*333  
    print(yvals)
    print(alpha_prior)
#    print(mu_prior)
    df_mu_prior = data.frame(xvals, mu_prior)
    df_alpha_prior = data.frame(yvals, alpha_prior)

    df_param = data.frame(mu, alpha)
    kd <- ks::kde(df_param, compute.cont=TRUE, positive=TRUE)
    contour_90 =with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]], 
               z=estimate, levels=cont["10%"])[[1]])
    contour_90 = data.frame(contour_90)
    # Need to get second polygon as first is tiny and disappears somewhere
#    contour_90_2 =with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
#    	       z=estimate, levels=cont["10%"])[[2]])
#    contour_90_2 = data.frame(contour_90_2)
    contour_50 =with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
               z=estimate, levels=cont["50%"])[[1]])
    contour_50 = data.frame(contour_50) 
    contour_25 =with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
    	       z=estimate, levels=cont["75%"])[[1]])
    contour_25 = data.frame(contour_25)   
    p1 = ggplot(df_param, aes(x=mu, y=alpha)) +
               stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
               geom_path(aes(x=x, y=y), data=contour_90, lwd=0.25) +
#	       geom_path(aes(x=x, y=y), data=contour_90_2, lwd=0.25) +
	       geom_path(aes(x=x, y=y), data=contour_50, lwd=0.5) + 
               geom_path(aes(x=x, y=y), data=contour_25, lwd=0.75) +
               scale_fill_distiller(palette="Greys", direction=1) +
               labs(colour = "Density") +
               xlab(expression("Mean ("*mu*")")) +
               ylab(expression("Aperiodicity ("*alpha*")")) +
               scale_x_continuous(expand = c(0, 0), limits = c(0, 120),
	       				 breaks=c(0, 50, 100), labels=c("0", "50", "100")) +
               scale_y_continuous(expand = c(0, 0), limits = c(0, 10), breaks = c(0,2,4,6,8,10),
	       				 labels=c("0","2","4","6","8","10")) +
	       theme(
	           legend.position='none'
  	       )+
#	       geom_point(x = mean_mu, y=mean_alpha, shape=4, colour='red') +
#	       geom_text(label=fig_lab, x=lab_x, y=lab_y)
	       geom_path(data=df_mu_prior, aes(x = xvals, y = mu_prior), colour='black') +
               geom_path(data=df_alpha_prior, aes(x = alpha_prior, y = yvals), colour='black') +
               geom_point(x = mean_mu, y=mean_alpha, shape=0, colour='red') +
               geom_point(x = median_mu, y=median_alpha, shape=1, colour='blue') +
               geom_point(x = mode[1], y=mode[2], shape=2, colour='green') +
              geom_text(label=fig_lab, x=lab_x, y=lab_y)
    return(p1)
    }

labels = c('a', 'b', 'c', 'd', 'e')
i=1
plot_list <- vector("list", length(posterior_files)+1)
for (filename in posterior_files){
    df_post = read.csv(filename)
    print(head(df_post))
    l = paste0(labels[[i]], ')')
    if (i==1){
       mu = df_post$mu
       lambda = df_post$lambda
       y = df_post$mre
       alpha = df_post$alpha
       pl = plot_posterior_2d(mu=mu, alpha=alpha, fig_lab=l)
       plot_list[[i]] = pl
    }else {
#       mu = c(mu, df_post$mu)
       mu = df_post$mu 
       lambda = c(lambda, df_post$lambda)
       y = c(y, df_post$mre)
#       alpha = c(alpha, df_post$alpha)
       alpha = df_post$alpha
       pt = plot_posterior_2d(mu=mu, alpha=alpha, fig_lab=l)
       plot_list[[i]] = pt 
       }
    i = i+1
    }

# Now do combined plot
l = paste0(labels[[i]], ')') 
p_combined = plot_posterior_2d(mu, alpha, fig_lab=l, lab_x=7, lab_y=9.8)# +
#	   geom_text(label=l, x=7, y=9.8)
plot_list[[i]] = p_combined
figname = 'plots/posterior_density_dunstan.png'
#figname = 'plots/posterior_density_eq_only_dunstan.png'
png(file=figname, units="in", width=7, height=5, res=300)
#pl = list(pl, pt)
print(plot_list)
grid.arrange(#pl, pt, nrow=2)
    grobs = plot_list,
    widths=c(0.1,1,1,2,0.1),
    layout_matrix = rbind(c(NA, 1, 2, 5, NA),
    		    	  c(NA, 3, 4, 5, NA))
    )
print(warnings())
warnings()
dev.off()