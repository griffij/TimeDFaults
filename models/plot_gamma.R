# Plot pdf of gamma distribution
alphas = list(2.001, 2.001, 2.063, 2.129, 2.129)
thetas = list(3.182e-03, 3.404e-03, 3.284e-03, 3.182e-03,3.404e-03)

pdf('gamma_density.pdf')
x = cumsum(rep(1, 2000))
#fn = dgamma(x, alpha, scale=1/theta)
#plot(density(fn))
xrange = range(0,2000)
yrange = range(0, 0.0015)
plot(xrange, yrange, type='n', xlab = 'Interevent time (years)',
	ylab = 'Probability')
for (i in seq_along(alphas)) {
    alpha = as.numeric(alphas[i])
    theta = as.numeric(thetas[i])
    print(alpha)
    print(theta[1])
    fn = dgamma(x, alpha, scale=1/theta)
    lines(x, fn)
    }

dev.off()