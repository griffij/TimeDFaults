# Plot pdf of gamm distribution
alpha = 2.063
theta = 3.284e-03

pdf('gamma_density.pdf')
x = cumsum(rep(1, 2000))
fn = dgamma(x, alpha, scale=1/theta)
#plot(density(fn))
plot(x, fn, type='l')
dev.off()