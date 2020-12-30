#############
### Setup ###
#############

## Preparation
inNames = outNames = list()

## Packages and Libraries
library("MASS")
library("spacejam")

## Command Line Arguments
source("./src/_commArgs.R")
# isce = 41

## Packages and Libraries
source("./simul_A1a_data.R")

## Randomness
set.seed(2020)


##############################
# Generating Gene Parameters #
##############################

##
mu = rnorm(p, 4)

##
G = rdag(p, p)
amG = as.matrix(get.adjacency(moralize(G)))

##
GeneParam = list(mu = mu, Omega = amG)
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
alphabeta = mvrnorm( n, c( -1.136, 0.900 ), matrix(
	c( 0.077, -0.006, -0.006, 0.002 ), 2, 2
) )
kappatau = mvrnorm( n, c( -0.520, 0.869 ), matrix(
	c( 0.058, 0.001, 0.001, 0.016 ), 2, 2
) )
alpha = alphabeta[,1]
beta = alphabeta[,2]
kappa = kappatau[,1]
tau = kappatau[,2]
rm(alphabeta, kappatau)

##
CellParam = data.frame(
	alpha=alpha, beta=beta, kappa=kappa, tau=tau
)
fileName = paste(
	"../", "Simul", "/", "A", "1",
	"/", "a", "_", "CellParam",
	"_", "sce", isce,
	".rds", sep=""
)
saveRDS(CellParam, fileName)
rm(fileName, CellParam)


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
VV = generate.dag.data(G, n * nrep, bfuns=function(x){cbind(x)})$X
VV = VV + outer(rep(1, n * nrep), mu)
for (irep in 1:nrep) {
	V = VV[(irep - 1) * n + 1:n,]
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


###############
# Termination #
###############
warnings()
