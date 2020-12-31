#############
### Setup ###
#############

## Preparation
inNames = outNames = list()
setwd("~/project/Jihwan/scLGM/Code")

## Packages and Libraries
library("LAM")
source("./src/ROC.R")
#source("./src/scLGM_geneParam.R")

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


##############
### Result ###
##############

## Data Processing
classIds = which( classData$level2class == oneClass )
geneIds = geneSelec[,oneClass,sId,gId]
oneDataMrna = mrnaData[classIds,geneIds]

## Result Aggregation
aggrCompGlasso = rep(NA, 25)
oneData = log( oneDataMrna + 1e-7 )
M = colMeans(oneData)
S = cov(oneData)
for (tId in 1:25) {

	## Input : gene parameter estimate
	inNames$geneParam = paste(
		"../", "Real", "/", "A", "5",
		"/", "a", "/", "comp", "Glasso",
		"_", oneClass, "_", sId, "_", mId, "_", gId, "_", tId,
		".rds", sep=""
	)
	geneParam = readRDS(inNames$geneParam)

	## Quality of Tuning
	
	aggrCompGlasso[tId] = sum( ( geneParam$wi != 0 ) & ( col(S) > row(S) ) )
	aggrCompGlasso[tId] = -2 * (
		loglike_mvnorm( M, S, rep( 0, ncol(oneDataMrna) ), geneParam$w, nrow(oneDataMrna) )
	) + aggrCompGlasso[tId] * log( nrow(oneDataMrna) )

}


##############
### Output ###
##############

## Output : gene parameter estimate
outNames$aggrCompGlasso = paste(
	"../", "Real", "/", "A", "6",
	"/", "c", "/", "aggr", "_", "comp", "Glasso",
	"_", oneClass, "_", sId, "_", mId, "_", gId,
	".rds", sep=""
)
saveRDS(aggrCompGlasso, outNames$aggrCompGlasso)
