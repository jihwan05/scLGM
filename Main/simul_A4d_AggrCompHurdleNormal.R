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
if (isce %% 10 == 1) { nreg = 84
} else if (isce %% 10 == 2) { nreg=100
} else if (isce %% 10 == 3) { nreg=123
} else { nreg=114 }

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
	YObs = scale( log( YObs + 1e-7 ) )
	M = colMeans(YObs)
	S = cov(YObs)

	## Gene Parameters
	fileName = paste(
		"../", "Simul", "/", "A", "3",
		"/", "b", "/", "Comp", "HurdleNormal",
		"_", "sce", isce, "rep", irep,
		".rds", sep=""
	)
	Fits = readRDS(fileName)
	rm(fileName)
	for (ireg in 1:nreg) {
		AdjMats[,,irep,ireg] = as.matrix( Fits$adjMat[[ireg]] )
	}
	BICs[irep,] = Fits$BIC

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
	"/", "d", "/", "Aggr", "Comp", "HurdleNormal",
	"_", "sce", isce,
	".rds", sep=""
)
saveRDS(ROCs, fileName)
rm(fileName)
