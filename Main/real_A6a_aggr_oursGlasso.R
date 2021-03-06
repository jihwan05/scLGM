#############
### Setup ###
#############

## Preparation
inNames = outNames = list()
setwd("~/project/Jihwan/scLGM/Code")

## Packages and Libraries
source("./src/scLGM_geneParam.R")
source("./src/ROC.R")

## Command Line Arguments
source("./src/_commArgs.R")
# cId = 1
# sId = 1
# mId = 1
# gId = 1


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
oneDir = c("a", "b", "c", "d")[mId]
oneMethod = c("Top", "Rand", "Equal", "Quart")[mId]
inNames$geneSelec = paste(
	"../", "Real", "/", "A", "2",
	"/", oneDir, "_", "geneSelec", oneMethod,
	".rds", sep=""
)
geneSelec = readRDS(inNames$geneSelec)

## Input : cell parameter estimate
inNames$cellParam = paste(
	"../", "Real", "/", "A", "4",
	"/", "a", "_", "cellParam",
	"_", oneClass,
	".rds", sep=""
)
cellParam = readRDS(inNames$cellParam)


##############
### Result ###
##############

## Data Processing
classIds = which( classData$level2class == oneClass )
geneIds = geneSelec[,oneClass,sId,gId]
oneDataMrna = mrnaData[classIds,geneIds]

## Result Aggregation
aggrOursGlasso = rep(NA, 25)
for (tId in 1:25) {

	## Input : gene parameter estimate
	inNames$geneParam = paste(
		"../", "Real", "/", "A", "4",
		"/", "b", "/", "oursGene", "Glasso",
		"_", oneClass, "_", sId, "_", mId, "_", gId, "_", tId,
		".rds", sep=""
	)
	geneParam = readRDS(inNames$geneParam)

	## Quality of Tuning
	aggrOursGlasso[tId] = BIC(
		geneParam$mu, geneParam$Omega, geneParam$V, geneParam$W,
		oneDataMrna,
		cellParam$alpha, cellParam$beta, cellParam$kappa, cellParam$tau
	)

}


##############
### Output ###
##############

## Output : gene parameter estimate
outNames$aggrOursGlasso = paste(
	"../", "Real", "/", "A", "6",
	"/", "a", "/", "aggr", "_", "ours", "Glasso",
	"_", oneClass, "_", sId, "_", mId, "_", gId,
	".rds", sep=""
)
saveRDS(aggrOursGlasso, outNames$aggrOursGlasso)
