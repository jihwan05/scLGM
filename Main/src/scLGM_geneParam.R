library("MASS")
library("msm") # Truncated Multivariate Normal Distribution
library("glasso")
library("fastclime")
library("flare")

algR = function(Vmean, VvarD, Wprob,
		Y, N, alpha, beta) {
	n = nrow(Y)
	p = ncol(Y)
	loglambda = alpha + beta * Vmean
	rho = beta^2 * VvarD + (loglambda - log(N))^2
	rho = sqrt(abs(Wprob * rho))
	Rmean = matrix(N / 4, n, p)
	non0id = which(rho > 1e-5, arr.ind=T)
	Rmean[non0id] =	(N / 2 / rho * tanh(rho/2))[non0id]
	return(list(rho = rho, mean = Rmean))
}

algW = function(Vmean, VvarD, Rmean, Y, N,
		alpha, beta, kappa, tau) {
	n = nrow(Y)
	p = ncol(Y)
	psi = kappa + tau * Vmean
	loglambda = alpha + beta * Vmean
	nu = N*log(sqrt(N)/2) - N/2*loglambda - Rmean/2*
			(beta^2 * VvarD + (loglambda - log(N))^2)
	expnu = exp(nu)
	PP = matrix(pnorm(-psi), n, p)
	DD = matrix(dnorm(-psi), n, p)
	Wprobneg = PP
	Wprobpos = (1 - PP) * expnu
	Wmean = (PP * psi - DD) +
			expnu * ((1 - PP) * psi + DD)
	Wmean = Wmean / (Wprobneg + Wprobpos)
	Wvar = (expnu + (1 - expnu) * PP) * (psi^2 + 1) -
			(1 - expnu) * psi * DD
	Wvar = Wvar / (Wprobneg + Wprobpos) - Wmean^2
	Wprob = Wprobpos / (Wprobneg + Wprobpos)
	non0id = which(Y != 0, arr.ind=T)
	Wmean[non0id] = (psi + DD / (1 - PP))[non0id]
	Wvar[non0id] = (1 - psi*DD/(1-PP) - (DD/(1-PP))^2)[non0id]
	Wprob[non0id] = 1
	return(list(nu = nu, theta = psi,
			mean = Wmean, var=Wvar, prob = Wprob))
}

algV = function(Wmean, Wprob, Rmean, Y, N, mu, Omega,
		alpha, beta, kappa, tau) {
	n = nrow(Y)
	p = ncol(Y)
	Vmean = matrix(NA, n, p)
	Vvar = array(NA, c(n, p, p))
	VvarD = matrix(NA, n, p)
	Omegamu = c(Omega %*% mu)
	for (i in 1:n) {
		Vvar[i,,] = Omega + diag(tau[i]^2 +
				Wprob[i,] * Rmean[i,] * beta[i]^2)
		Vvar[i,,] = chol2inv(chol(Vvar[i,,]))
		Vmean[i,] = as.vector( Omegamu + tau[i] * (Wmean[i,] - kappa[i]) +
			Wprob[i,] * beta[i] * (
				as.numeric(Y[i,]) - N/2 - Rmean[i,] * (alpha[i] - log(N))
			) )
		Vmean[i,] = as.vector( Vvar[i,,] %*% as.vector( Vmean[i,] ) )
		VvarD[i,] = diag(Vvar[i,,])
	}
	return(list(mean = Vmean, var = Vvar, varD = VvarD))
}

checkConv = function(tp, mu, Omega, fitV, fitW, Y,
		alpha, beta, kappa, tau) {
	n = nrow(Y)
	p = ncol(Y)
	psi = kappa + tau * fitV$mean
	loglambda = alpha + beta * fitV$mean
	lambda = exp(loglambda)
	result = - tp * sum(abs(Omega[upper.tri(Omega)]))
	result = result + n / 2 * sum(log(
			Re(eigen(Omega, only.values=T)$values)))
	temp = 0; for (i in 1:n) {
		temp = temp + fitV$var[i,,]
	}
	temp = temp / n + cov(fitV$mean) * ((n-1) / n)
	result = result - n / 2 * sum(Omega * temp)
	result = result - sum(fitW$var) / 2
	result = result - sum(tau^2 * fitV$varD) / 2
	result = result - sum((fitW$mean - psi)^2) / 2
	result = result + sum(fitW$prob * (-lfactorial(Y)
			-lambda * exp(beta^2 / 2 * fitV$varD)
			+ Y * loglambda))
	return(result)			
}

algMainSmall = function(Y, tp, alpha, beta, kappa, tau,
		Mstep = "glasso", niter=20, printll=F) {
	n = nrow(Y)
	p = ncol(Y)
	N = 100 * max(Y)
	# Initialization: Latent
	non0id = which(Y != 0, arr.ind=T)
	Vmean = matrix(-kappa/tau, n, p)
	Vmean[non0id] = ((log(Y) - alpha) / beta)[non0id]
	VvarD = matrix(1, n, p)
	fitVOld = list(mean = Vmean, varD = VvarD)
	Wprob = (Y != 0) * 1
	fitWOld = list(prob = Wprob)
	rm(non0id, Vmean, VvarD, Wprob)
	# Initialization: Gene Parameter
	muOld = colMeans(fitVOld$mean)
	OmegaOld = diag(1/diag(cov(fitVOld$mean)))
	# Algorithm
	llOld = -Inf
	for (it in 1:niter) {
		# E-step
		fitRNew = algR(fitVOld$mean, fitVOld$varD, fitWOld$prob,
				Y, N, alpha, beta)
		fitWNew = algW(fitVOld$mean, fitVOld$varD, fitRNew$mean,
				Y, N, alpha, beta, kappa, tau)
		fitVNew = algV(fitWNew$mean, fitWNew$prob, fitRNew$mean,
				Y, N, muOld, OmegaOld, alpha, beta, kappa, tau)
		# M-step
		muNew = colMeans(fitVNew$mean)
		S = cov(fitVNew$mean) * (n - 1) / n
		for (i in 1:n) {
			S = S + fitVNew$var[i,,]/n
		}
#		if (Mstep == "fastclime") {
#			OmegaNew = fastclime(S, tp)
#			OmegaNew = OmegaNew$icovlist[[1]]
#		} else
		if (Mstep == "flare") {
			OmegaNew = sugm(S, method="clime", lambda=tp)
			OmegaNew = (OmegaNew$icov)[[1]]
		} else {
			#OmegaNew = glasso(S, tp, penalize.diagonal=T)
			ttp = matrix( tp, nrow(S), ncol(S) )
			diag(ttp) = 1e-7
			OmegaNew = glasso(S, ttp, penalize.diagonal=T)
			OmegaNew = OmegaNew$wi
		}
		# Convergence
		llNew = checkConv(tp, muNew, OmegaNew,
				fitVNew, fitWNew, Y,
				alpha, beta, kappa, tau)
		if (printll == T) {
			print(c(it, llOld, llNew, llNew-llOld))
		}
		if (it >= 5) {
			if (abs(llNew - llOld) < 1) {
				print(paste("It took", it,
						"iterations to converge.", sep=" "))
				break
			} else if (abs(llNew - llOld) < 0.01*abs(llOld)) {
				print(paste("It took", it,
					"iterations to converge.", sep=" "))
				break
			}
		}
		# Re-iteration
		llOld = llNew
		muOld = muNew
		OmegaOld = OmegaNew
		fitROld = fitRNew
		fitWOld = fitWNew
		fitVOld = fitVNew
	}
	return(list(mu = muOld, Omega = OmegaOld,
			V = fitVOld, W = fitWOld, R = fitROld))
}

logLike = function(mu, Omega, fitV, fitW, Y,
		alpha, beta, kappa, tau) {
	n = nrow(Y)
	p = ncol(Y)
	psi = kappa + tau * fitV$mean
	loglambda = alpha + beta * fitV$mean
	lambda = exp(loglambda)
	result = - n * p * log(2 * pi)
	result = result + n / 2 * sum(log(
			Re(eigen(Omega, only.values=T)$values)))
	temp = 0; for (i in 1:n) {
		temp = temp + fitV$var[i,,]
	}; temp = temp / n + cov(fitV$mean) * (n-1) / n
	result = result - n / 2 * sum(Omega * temp)
	result = result - sum(fitW$var) / 2
	result = result - sum(tau^2 * fitV$varD) / 2
	result = result - sum((fitW$mean - psi)^2) / 2
	result = result + sum(fitW$prob * (-lfactorial(Y)
			-lambda * exp(beta^2 / 2 * fitV$varD)
			+ Y * loglambda))
	return(result)			
}

BIC = function(mu, Omega, fitV, fitW, Y,
		alpha, beta, kappa, tau) {
	n = nrow(Y)
	p = ncol(Y)
	result = -2 * logLike(mu, Omega, fitV, fitW, Y,
			alpha, beta, kappa, tau) + log(n) *
			(p + p + sum(Omega[upper.tri(Omega)] != 0))
	return(result)
}

algMain = function(Y, tp, alpha, beta, kappa, tau,
		Mstep = "glasso", niter=20, printll=F) {
	n = nrow(Y)
	p = ncol(Y)
	N = 100 * max(Y)
	# Initialization: Latent
	non0id = which(Y != 0, arr.ind=T)
	Vmean = matrix(-kappa/tau, n, p)
	Vmean[non0id] = ((log(Y) - alpha) / beta)[non0id]
	VvarD = matrix(1, n, p)
	fitVOld = list(mean = Vmean, varD = VvarD)
	Wprob = (Y != 0) * 1
	fitWOld = list(prob = Wprob)
	rm(non0id, Vmean, VvarD, Wprob)
	# Initialization: Gene Parameter
	muOld = colMeans(fitVOld$mean)
	OmegaOld = diag(1/diag(cov(fitVOld$mean)))
	# Algorithm
	llOld = -Inf
	for (it in 1:niter) {
		# E-step
		fitRNew = algR(fitVOld$mean, fitVOld$varD, fitWOld$prob,
				Y, N, alpha, beta)
		fitWNew = algW(fitVOld$mean, fitVOld$varD, fitRNew$mean,
				Y, N, alpha, beta, kappa, tau)
		fitVNew = algV(fitWNew$mean, fitWNew$prob, fitRNew$mean,
				Y, N, muOld, OmegaOld, alpha, beta, kappa, tau)
		# M-step
		muNew = colMeans(fitVNew$mean)
		S = cov(fitVNew$mean) * (n - 1) / n
		for (i in 1:n) {
			S = S + fitVNew$var[i,,]/n
		}
#		if (Mstep == "fastclime") {
#			OmegaNew = fastclime(S, tp)
#			OmegaNew = OmegaNew$icovlist[[1]]
#		} else
		if (Mstep == "flare") {
			OmegaNew = sugm(S, method="clime", lambda=tp)
			OmegaNew = (OmegaNew$icov)[[1]]
		} else {
			OmegaNew = glasso(S, tp)
			OmegaNew = OmegaNew$wi
		}
		# Convergence
		llNew = checkConv(tp, muNew, OmegaNew,
				fitVNew, fitWNew, Y,
				alpha, beta, kappa, tau)
		if (printll == T) {
			print(c(it, llOld, llNew, llNew-llOld))
		}
		if (it >= 5) {
			if (abs(llNew - llOld) < 1) {
				print(paste("It took", it,
						"iterations to converge.", sep=" "))
				break
			} else if (abs(llNew - llOld) < 0.01*abs(llOld)) {
				print(paste("It took", it,
					"iterations to converge.", sep=" "))
				break
			}
		}
		# Re-iteration
		llOld = llNew
		muOld = muNew
		OmegaOld = OmegaNew
		fitROld = fitRNew
		fitWOld = fitWNew
		fitVOld = fitVNew
	}
	return(list(mu = muOld, Omega = OmegaOld,
			V = fitVOld, W = fitWOld, R = fitROld))
}
