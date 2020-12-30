library("msm") # Truncated Multivariate Normal Distribution

scLGM_cellParam = function(x, Y, nMCMC = 3000, bMCMC = 2000) {
	n = ncol(Y)
	q = nrow(Y)
	nq = n * q
	kappaMCMC = tauMCMC = matrix(NA, n, nMCMC+bMCMC)
	# Initial Estimation
	alpha = beta = kappa = tau = rep(NA, n)
	for (i in 1:n) {
		y = Y[,i]
		non0id = (y != 0)
		fitkt = glm( c(non0id*1) ~ log(x + 1e-7),
					family=binomial(link="probit"))
		kappa[i] = coefficients(fitkt)[1]
		tau[i] = coefficients(fitkt)[2]
		fitab = lm(log(y[non0id]) ~ log(x[non0id] + 1e-7))
		alpha[i] = coefficients(fitab)[1]
		beta[i] = coefficients(fitab)[2]
	}
	# MCMC
	v = log(x + 1e-7)
	W = matrix(NA, n, q)
	kappamean = mean(kappa)
	kappavar = var(kappa)
	taumean = mean(tau)
	tauvar = var(tau)
	loglambda = alpha + outer(beta, v)
	lambda = exp(loglambda)
	expneglambda = exp(-lambda)
	for (iMCMC in 1:(nMCMC+bMCMC)) {
		tau_v = outer(tau, v)
		psi = kappa + tau_v
		Wpos = matrix(rtnorm(nq, psi, 1, 0, Inf), n, q)
		Wneg = matrix(rtnorm(nq, psi, 1, -Inf, 0), n, q)
		Wpneg = matrix(pnorm(0, psi, 1), n, q)
		Wpneg = Wpneg / (Wpneg + (1-Wpneg)*expneglambda)
		Windicator = matrix(rbinom(nq, 1, Wpneg), n, q)
		W = Wpos
		negid = which(t(Y==0) & (Windicator==1), arr.ind=T)
		W[negid] = Wneg[negid]
		sigma2 = 1 / (q + 1/kappavar)
		mu = sigma2 * (rowSums(W-tau_v) + kappamean/kappavar)
		kappa = kappaMCMC[,iMCMC] = rnorm(n, mu, sqrt(sigma2))
		sigma2 = 1 / (sum(v^2) + 1/tauvar)
		mu = sigma2 * (c((W-kappa) %*% v) + taumean/tauvar)
		tau = tauMCMC[,iMCMC] = rnorm(n, mu, sqrt(sigma2))
	}
	return( data.frame(
		alpha = alpha, beta = beta,
		kappa = rowMeans( kappaMCMC[,-(1:bMCMC)] ),
		tau = rowMeans( tauMCMC[ , -(1:bMCMC) ] ) 
	) )
}
