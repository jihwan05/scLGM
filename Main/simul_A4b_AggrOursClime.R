#############
### Setup ###
#############

## Preparation
inNames = outNames = list()
setwd("~/project/Jihwan/scLGM/Code")

## Command Line Arguments
source("./src/_commArgs.R")
# isce = 11

## Packages and Libraries
source("./src/scLGM_geneParam.R")
source("./src/ROC.R")
source("./simul_A1a_data.R")

## Randomness
set.seed(2020)


#############
### Input ###
#############

## Gene Parameters
fileName = paste(
	"../", "Simul", "/", "A", "1",
	"/", "a", "_", "GeneParam",
	"_", "sce", isce,
".rds", sep=""
)
GeneParam = readRDS(fileName)
rm(fileName)


##################
### OursGlasso ###
##################

##
AdjMats = array(NA, c(p, p, nrep, nreg))
BICs = matrix(NA, nrep, nreg)

##
for (irep in 1:nrep) {

	## Data
	fileName = paste(
		"../", "Simul", "/", "A", "1",
		"/", "a", "_", "Data",
		"_", "sce", isce,
		".rds", sep=""
	)
	YObs = readRDS(fileName)[[irep]]$Y
	rm(fileName)

	## Cell Parameters
	fileName = paste(
		"../", "Simul", "/", "A", "2",
		"/", "a", "/", "Ours", "Cell",
		"_", "sce", isce, "rep", irep,
	".rds", sep=""
	)
	CellParam = readRDS(fileName)
	rm(fileName)

	## Gene Parameters
	for (ireg in 1:nreg) {
		fileName = paste(
			"../", "Simul", "/", "A", "2",
			"/", "c", "/", "Ours", "Gene", "Clime",
			"_", "sce", isce, "rep", irep, "reg", ireg,
			".rds", sep=""
		)
		Fit = readRDS(fileName)
		rm(fileName)
		AdjMats[,,irep,ireg] = (Fit$Omega != 0)
		BICs[irep, ireg] = BIC(
			Fit$mu, Fit$Omega, Fit$V, Fit$W, YObs,
			CellParam$alpha, CellParam$beta, CellParam$kappa, CellParam$tau
		)
	}
}

##
AdjMatTrue = ( GeneParam$Omega != 0 )
ROCs = ROCs(AdjMatTrue, AdjMats, regsOurs_Glasso)
AdjMatBest = array( NA, c(p, p, nrep) )
for (irep in 1:nrep) {
	AdjMatBest[,,irep] = AdjMats[
		, , irep, which.min( BICs[irep,] )
	]
}
ROC1 = ROC1(AdjMatTrue, AdjMatBest)


##############
### Output ###
##############

ROCs = list(curve=ROCs, bests=ROC1)
fileName = paste(
	"../", "Simul", "/", "A", "4",
	"/", "b", "/", "Aggr", "Ours", "Clime",
	"_", "sce", isce,
	".rds", sep=""
)
saveRDS(ROCs, fileName)
rm(fileName)
