#############
### Setup ###
#############

## Preparation
inNames = outNames = list()

## Packages and Libraries
library("MASS")

## Command Line Arguments
#source("./src/_commArgs.R")
isce = 77

## Packages and Libraries
source("./simul_A1a_data.R")

## Randomness
set.seed(2020)


##############################
# Generating Gene Parameters #
##############################

## Input : gene parameter estimate
inNames$bestOursGlasso = paste(
	"../", "Real", "/", "A", "7",
	"/", "x", "_", "best", "Ours", "Glasso",
	".rds", sep=""
)
bestOursGlasso = readRDS(inNames$bestOursGlasso)

##
mu = bestOursGlasso$mu
Omega = bestOursGlasso$Omega

##
OmegaNew = matrix(0, nrow(Omega), ncol(Omega))
OmegaNewUt = Omega[upper.tri(Omega)]
OmegaNewUtNz = (Omega[upper.tri(Omega)] != 0)
OmegaNewUt[OmegaNewUtNz] = sin(2*pi*runif(sum(OmegaNewUtNz))) / 5
OmegaNew[upper.tri(Omega)] = OmegaNewUt
OmegaNew = OmegaNew + t(OmegaNew)
diag(OmegaNew) = abs(min(eigen(OmegaNew)$values))+0.01
OmegaNew = diag(1/sqrt(diag(OmegaNew))) %*%
		OmegaNew %*% diag(1/sqrt(diag(OmegaNew)))
Omega = diag(sqrt(diag(Omega))) %*%
		OmegaNew %*% diag(sqrt(diag(Omega)))
Sigma = chol2inv(chol(Omega))

##
GeneParam = list(mu = mu, Omega = Omega)
fileName = paste(
	"../", "Simul", "/", "A", "1",
	"/", "a", "_", "GeneParam",
	"_", "sce", isce,
	".rds", sep=""
)
saveRDS(GeneParam, fileName)
rm(fileName, GeneParam)


##############################
# Generating Cell Parameters #
##############################

##
fileName = "../Real/A4/a_cellParam_CA1Pyr2.rds"
CellParam = readRDS(fileName)
rm(fileName)

##
alpha = CellParam$alpha
beta = CellParam$beta
kappa = CellParam$kappa
tau = CellParam$tau
CellParam = data.frame(
	alpha=alpha, beta=beta, kappa=kappa, tau=tau
)

##
fileName = paste(
	"../", "Simul", "/", "A", "1",
	"/", "a", "_", "CellParam",
	"_", "sce", isce,
	".rds", sep=""
)
saveRDS(CellParam, fileName)
rm(CellParam, fileName)


########################
# Generating Spike-ins #
########################

##
fileName = "../Real/A0/b_ospikeData.rds"
SpikeIns = readRDS(fileName)
rm(fileName)

##
SpikeInsInput = SpikeIns$mpcNew
rm(SpikeIns)
q = length(SpikeInsInput)
psi = kappa + outer(tau, log(SpikeInsInput + 1e-7))
loglambda = alpha + outer(beta, log(SpikeInsInput + 1e-7))
lambda = exp(loglambda)
SpikeInsOutput = NULL
for (irep in 1:nrep) {
	Wsi = matrix(rbinom(n*q, 1, pnorm(psi)), n, q)
	Ysi = matrix(rpois(n*q, lambda), n, q)
	Ysi[which(Wsi==0, arr.ind=T)] = 0
	SpikeInsOutput[[irep]] = t(Ysi)
}; rm(irep, Wsi, Ysi)
rm(q, psi, loglambda, lambda)

##
SpikeIns = list(input=SpikeInsInput, output=SpikeInsOutput)
fileName = paste(
	"../", "Simul", "/", "A", "1",
	"/", "a", "_", "SpikeIns",
	"_", "sce", isce,
	".rds", sep=""
)
saveRDS(SpikeIns, fileName)
rm(fileName, SpikeIns)
rm(SpikeInsInput, SpikeInsOutput)


##########################
# Generating mRNA Counts #
##########################

##
Data = NULL
for (irep in 1:nrep) {
	V = mvrnorm(n, mu, Sigma)
	psi = kappa + tau * V
	lambda = exp(alpha + beta * V)
	W = matrix(rnorm(n*p, psi, 1), n, p)
	negID = which(W < 0, arr.ind=T)
	Y = matrix(rpois(n*p, lambda), n, p)
	Y[negID] = 0
	Data[[irep]] = list(V=V, W=W, Y=Y)
}; rm(irep, psi, lambda, negID, V, W, Y)

##
fileName = paste(
	"../", "Simul", "/", "A", "1",
	"/", "a", "_", "Data",
	"_", "sce", isce,
	".rds", sep=""
)
saveRDS(Data, fileName)
rm(fileName, Data)
