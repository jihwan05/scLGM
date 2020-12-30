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
library("LAM")
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

	##
	YObs = log( YObs + 1e-7 )
	M = colMeans(YObs)
	S = cov(YObs)

	## Gene Parameters
	for (ireg in 1:nreg) {
		fileName = paste(
			"../", "Simul", "/", "A", "3",
			"/", "a", "/", "Comp", "Glasso",
			"_", "sce", isce, "rep", irep, "reg", ireg,
			".rds", sep=""
		)
		Fit = readRDS(fileName)
		rm(fileName)
		AdjMats[,,irep,ireg] = Fit$wi
		BICs[irep,ireg] = sum( ( Fit$wi != 0 ) & ( col(S) > row(S) ) )
		BICs[irep,ireg] = -2 * ( loglike_mvnorm(
			M, S, rep( 0, ncol(YObs) ), Fit$w, nrow(YObs)
		) ) + BICs[irep,ireg] * log( nrow(YObs) )
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
	"/", "c", "/", "Aggr", "Comp", "Glasso",
	"_", "sce", isce,
	".rds", sep=""
)
saveRDS(ROCs, fileName)
rm(fileName)
