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

## Input : class data
inNames$classData = paste("../Real/A0/a_classData.rds")
classData = readRDS(inNames$classData)

## Input : mRNA UMI count data
inNames$mrnaData = paste("../Real/A0/a_mrnaData.rds")
mrnaData = readRDS(inNames$mrnaData)

## Input : list of special classes : big three
inNames$specClass = paste("../Real/A1/b_specClass.rds")
specClass = readRDS(inNames$specClass)
oneClass = specClass[cId]

## Input : selections of 200 genes
inNames$geneSelec = paste("../Real/A2/a_geneSelecTop.rds")
geneSelec = readRDS(inNames$geneSelec)


##############
### Output ###
##############

## Input : gene parameter estimate : oursGlasso
inNames$aggrOursGlasso = paste(
	"../", "Real", "/", "A", "6",
	"/", "a", "/", "aggr", "_", "ours", "Glasso",
	"_", oneClass, "_", sId, "_", mId, "_", 1,
	".rds", sep=""
)
aggrOursGlasso = readRDS(inNames$aggrOursGlasso)

## Best : gene parameter estimate : oursGlasso
minOursGlasso = which.min(aggrOursGlasso)[1]
inNames$bestOursGlasso = paste(
	"../", "Real", "/", "A", "4",
	"/", "b", "/", "oursGene", "Glasso",
	"_", oneClass, "_", sId, "_", mId, "_", 1, "_", minOursGlasso,
	".rds", sep=""
)
bestOursGlasso = readRDS(inNames$bestOursGlasso)

## Output : gene parameter estimate
outNames$bestOursGlasso = paste(
	"../", "Real", "/", "A", "7",
	"/", "x", "_", "best", "Ours", "Glasso",
	".rds", sep=""
)
saveRDS(bestOursGlasso, outNames$bestOursGlasso)

