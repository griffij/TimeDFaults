# Plot density of multiple runs on one figures

library(ggplot2)
library(ks)
library(lattice)
library(grid)
library(gridExtra)

# Read in posterior dataset
posterior_files = c('outputs/df_posterior_1_dunstan.csv',
              'outputs/df_posterior_2_dunstan.csv',
              'outputs/df_posterior_3_dunstan.csv',
              'outputs/df_posterior_4_dunstan.csv')
#MRE_position = c(5,6,6,7)             

plot_posterior_2d <-function(mu, alpha, fig_lab, lab_x=15000, lab_y=9.5){
    # Estimate density value containing 95% of posterior ditribution
    # Get mean values
    mean_mu = median(mu)
    mean_alpha = median(alpha)
    df_param = data.frame(mu, alpha)
    kd <- ks::kde(df_param, gridsize=rep(1401,2) , compute.cont=TRUE)
    contour_95 =with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]], 
               z=estimate, levels=cont["5%"])[[1]])
    contour_95 = data.frame(contour_95)
    contour_68 =with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
    	       z=estimate, levels=cont["32%"])[[1]])
    contour_68 = data.frame(contour_68)   
    p1 = ggplot(df_param, aes(x=mu, y=alpha)) +
               stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
               geom_path(aes(x=x, y=y), data=contour_95, lwd=0.25) +
               geom_path(aes(x=x, y=y), data=contour_68, lwd=0.5) +
               scale_fill_distiller(palette="Greys", direction=1) +
               labs(colour = "Density") +
               xlab(expression("Mean ("*mu*")")) +
               ylab(expression("Aperiodicity ("*alpha*")")) +
               scale_x_continuous(expand = c(0, 0), limits = c(0, 150000),
	       				 breaks=c(0, 50000, 100000), labels=c("0", "50000", "100000")) +
               scale_y_continuous(expand = c(0, 0), limits = c(0, 10), breaks = c(0,2,4,6,8,10),
	       				 labels=c("0","2","4","6","8","10")) +
	       theme(
	           legend.position='none'
  	       )+
	       geom_point(x = mean_mu, y=mean_alpha, shape=4, colour='red') +
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
       mu = c(mu, df_post$mu)
       lambda = c(lambda, df_post$lambda)
       y = c(y, df_post$mre)
       alpha = c(alpha, df_post$alpha)
       pt = plot_posterior_2d(mu=df_post$mu, alpha=df_post$alpha, fig_lab=l)
       plot_list[[i]] = pt 
       }
    i = i+1
    }

# Now do combined plot
l = paste0(labels[[i]], ')') 
p_combined = plot_posterior_2d(mu, alpha, fig_lab=l, lab_x=7000, lab_y=9.8)# +
#	   geom_text(label=l, x=7000, y=9.8)
plot_list[[i]] = p_combined
figname = 'plots/posterior_density_dunstan.png'
png(file=figname, units="in", width=7, height=5, res=300)
#pl = list(pl, pt)
print(plot_list)
grid.arrange(#pl, pt, nrow=2)
    grobs = plot_list,
    widths=c(0.1,1,1,2,0.1),
    layout_matrix = rbind(c(NA, 1, 2, 5, NA),
    		    	  c(NA, 3, 4, 5, NA))
    )

dev.off()