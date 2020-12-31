#############
### Setup ###
#############

## Preparation
inNames = outNames = list()
setwd("~/project/Jihwan/scLGM/Code")

## Packages and Libraries
source("./src/ROC.R")

## Command Line Arguments
source("./src/_commArgs.R")
cId = 1
sId = 1
mId = 3


#############
### Input ###
#############

## Input : list of special classes : big three
inNames$specClass = paste("../Real/B1/b_specClass.rds")
specClass = readRDS(inNames$specClass)
oneClass = specClass[cId]


##############
### Result ###
##############

## Preparation
minOursGlasso = rep(NA, 5)
rocOursGlasso = NULL
minOursClime = rep(NA, 5)
rocOursClime = NULL
minCompGlasso = rep(NA, 5)
rocCompGlasso = NULL
minCompHurdleNormal = rep(NA, 5)
rocCompHurdleNormal = NULL

## Result Aggregation
for ( gId in 1:5 ) {

	## Input : graph information from KEGG database
	inNames$keggGraph = paste(
		"../", "Real", "/", "B", "3",
		"/", "c", "/", "keggGraph", "Equal",
		"_", oneClass, "_", sId, "_", gId,
		".rds", sep=""
	)
	keggGraph = readRDS(inNames$keggGraph)
	
	## Input : gene parameter estimate : oursGlasso
	inNames$aggrOursGlasso = paste(
		"../", "Real", "/", "B", "6",
		"/", "a", "/", "aggr", "_", "ours", "Glasso",
		"_", oneClass, "_", sId, "_", mId, "_", gId,
		".rds", sep=""
	)
	aggrOursGlasso = readRDS(inNames$aggrOursGlasso)

	## Input : gene parameter estimate : oursClime
	inNames$aggrOursClime = paste(
		"../", "Real", "/", "B", "6",
		"/", "b", "/", "aggr", "_", "ours", "Clime",
		"_", oneClass, "_", sId, "_", mId, "_", gId,
		".rds", sep=""
	)
	aggrOursClime = readRDS(inNames$aggrOursClime)

	## Input : gene parameter estimate : compGlasso
	inNames$aggrCompGlasso = paste(
		"../", "Real", "/", "B", "6",
		"/", "c", "/", "aggr", "_", "comp", "Glasso",
		"_", oneClass, "_", sId, "_", mId, "_", gId,
		".rds", sep=""
	)
	aggrCompGlasso = readRDS(inNames$aggrCompGlasso)

	## Input : gene parameter estimate : compHurdleNormal
	inNames$compHurdleNormal = paste(
		"../", "Real", "/", "B", "5",
		"/", "b", "/", "comp", "HurdleNormal",
		"_", oneClass, "_", sId, "_", mId, "_", gId,
		".rds", sep=""
	)
	aggrCompHurdleNormal = readRDS(inNames$compHurdleNormal)

	## Best : gene parameter estimate : oursGlasso
	minOursGlasso[gId] = which.min(aggrOursGlasso)[1]
	inNames$bestOursGlasso = paste(
		"../", "Real", "/", "B", "4",
		"/", "b", "/", "oursGene", "Glasso",
		"_", oneClass, "_", sId, "_", mId, "_", gId, "_", minOursGlasso[gId],
		".rds", sep=""
	)
	bestOursGlasso = readRDS(inNames$bestOursGlasso)

	## Best : gene parameter estimate : oursClime
	minOursClime[gId] = which.min(aggrOursClime)[1]
	inNames$bestOursClime = paste(
		"../", "Real", "/", "B", "4",
		"/", "c", "/", "oursGene", "Clime",
		"_", oneClass, "_", sId, "_", mId, "_", gId, "_", minOursClime[gId],
		".rds", sep=""
	)
	bestOursClime = readRDS(inNames$bestOursClime)

	## Best : gene parameter estimate : compGlasso
	minCompGlasso[gId] = which.min(aggrCompGlasso)[1]
	inNames$bestCompGlasso = paste(
		"../", "Real", "/", "B", "5",
		"/", "a", "/", "comp", "Glasso",
		"_", oneClass, "_", sId, "_", mId, "_", gId, "_", minCompGlasso[gId],
		".rds", sep=""
	)
	bestCompGlasso = readRDS(inNames$bestCompGlasso)

	## Best : gene parameter estimate : compHurdleNormal
	minCompHurdleNormal[gId] = which.min(aggrCompHurdleNormal$BIC)[1]
	bestCompHurdleNormal = as.matrix( aggrCompHurdleNormal$adjMat[[ minCompHurdleNormal[gId] ]] )
	
	## ROC
	graphGs = matrix(F, 200, 200)
	for(j in 1:nrow(keggGraph$E)) {
		graphGs[ keggGraph$E[j,1], keggGraph$E[j,2] ] = T
	}
	rocOursGlasso = rbind(
		rocOursGlasso,
		ROCnew( graphGs, bestOursGlasso$Omega )
	)
	rocOursClime = rbind(
		rocOursClime,
		ROCnew( graphGs, bestOursClime$Omega )
	)
	rocCompGlasso = rbind(
		rocCompGlasso,
		ROCnew( graphGs, bestCompGlasso$wi )
	)
	rocCompHurdleNormal = rbind(
		rocCompHurdleNormal,
		ROCnew( graphGs, bestCompHurdleNormal )
	)
	
}

##
ROCs = rbind(
	colMeans(rocOursGlasso),
	colMeans(rocOursClime),
	colMeans(rocCompHurdleNormal),
	colMeans(rocCompGlasso)
)
ROCfin = ROCs[,c(1,2,2,5)]
ROCfin[,2] = ROCfin[,3] + ROCfin[,4]
print( round( ROCfin, 3) )
