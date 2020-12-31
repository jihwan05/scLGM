#############
### Setup ###
#############

## Preparation
inNames = outNames = list()
setwd("~/project/Jihwan/scLGM/Code")

## Packages and Libraries
library("VennDiagram")
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
inNames$specClass = paste("../Real/A1/b_specClass.rds")
specClass = readRDS(inNames$specClass)
oneClass = specClass[cId]


##############
### Result ###
##############

## Preparation
minOursGlasso = rep(NA, 5)
minOursClime = rep(NA, 5)
minCompGlasso = rep(NA, 5)
minCompHurdleNormal = rep(NA, 5)

## Result Aggregation
overlapAggr = NULL
for ( gId in 1:5 ) {

	## Input : graph information from KEGG database
	inNames$keggGraph = paste(
		"../", "Real", "/", "A", "3",
		"/", "c", "/", "keggGraph", "Equal",
		"_", oneClass, "_", sId, "_", gId,
		".rds", sep=""
	)
	keggGraph = readRDS(inNames$keggGraph)
	
	## Input : gene parameter estimate : oursGlasso
	inNames$aggrOursGlasso = paste(
		"../", "Real", "/", "A", "6",
		"/", "a", "/", "aggr", "_", "ours", "Glasso",
		"_", oneClass, "_", sId, "_", mId, "_", gId,
		".rds", sep=""
	)
	aggrOursGlasso = readRDS(inNames$aggrOursGlasso)

	## Input : gene parameter estimate : oursClime
	inNames$aggrOursClime = paste(
		"../", "Real", "/", "A", "6",
		"/", "b", "/", "aggr", "_", "ours", "Clime",
		"_", oneClass, "_", sId, "_", mId, "_", gId,
		".rds", sep=""
	)
	aggrOursClime = readRDS(inNames$aggrOursClime)

	## Input : gene parameter estimate : compGlasso
	inNames$aggrCompGlasso = paste(
		"../", "Real", "/", "A", "6",
		"/", "c", "/", "aggr", "_", "comp", "Glasso",
		"_", oneClass, "_", sId, "_", mId, "_", gId,
		".rds", sep=""
	)
	aggrCompGlasso = readRDS(inNames$aggrCompGlasso)

	## Input : gene parameter estimate : compHurdleNormal
	inNames$compHurdleNormal = paste(
		"../", "Real", "/", "A", "5",
		"/", "b", "/", "comp", "HurdleNormal",
		"_", oneClass, "_", sId, "_", mId, "_", gId,
		".rds", sep=""
	)
	aggrCompHurdleNormal = readRDS(inNames$compHurdleNormal)

	## Best : gene parameter estimate : oursGlasso
	minOursGlasso[gId] = which.min(aggrOursGlasso)[1]
	inNames$bestOursGlasso = paste(
		"../", "Real", "/", "A", "4",
		"/", "b", "/", "oursGene", "Glasso",
		"_", oneClass, "_", sId, "_", mId, "_", gId, "_", minOursGlasso[gId],
		".rds", sep=""
	)
	bestOursGlasso = readRDS(inNames$bestOursGlasso)

	## Best : gene parameter estimate : oursClime
	minOursClime[gId] = which.min(aggrOursClime)[1]
	inNames$bestOursClime = paste(
		"../", "Real", "/", "A", "4",
		"/", "c", "/", "oursGene", "Clime",
		"_", oneClass, "_", sId, "_", mId, "_", gId, "_", minOursClime[gId],
		".rds", sep=""
	)
	bestOursClime = readRDS(inNames$bestOursClime)

	## Best : gene parameter estimate : compGlasso
	minCompGlasso[gId] = which.min(aggrCompGlasso)[1]
	inNames$bestCompGlasso = paste(
		"../", "Real", "/", "A", "5",
		"/", "a", "/", "comp", "Glasso",
		"_", oneClass, "_", sId, "_", mId, "_", gId, "_", minCompGlasso[gId],
		".rds", sep=""
	)
	bestCompGlasso = readRDS(inNames$bestCompGlasso)

	## Best : gene parameter estimate : compHurdleNormal
	minCompHurdleNormal[gId] = which.min(aggrCompHurdleNormal$BIC)[1]
	bestCompHurdleNormal = as.matrix(
		aggrCompHurdleNormal$adjMat[[ minCompHurdleNormal[gId] ]]
	)
	
	## Overlap
	overlapAggr[[gId]] = cbind(
		OursGlasso = bestOursGlasso$Omega[ upper.tri( bestOursGlasso$Omega ) ],
		OursClime = bestOursClime$Omega[ upper.tri( bestOursClime$Omega ) ],
		Glasso = bestCompGlasso$wi[ upper.tri( bestCompGlasso$wi ) ],
		HurdleNormal = bestCompHurdleNormal[ upper.tri( bestCompHurdleNormal ) ]
	)
	
}

##
outNames$overlapAggr = paste(
	"../", "Real", "/", "A", "8",
	"/", "c", "_", "overlap", "Equal",
	".rds", sep=""
)
saveRDS(overlapAggr, outNames$overlapAggr)
