#############
### Setup ###
#############

## Preparation
inNames = outNames = list()
setwd("~/project/Jihwan/scLGM/Code")

## Packages and Libraries
library("HurdleNormal")

## Command Line Arguments
source("./src/_commArgs.R")
# cId = 1
# sId = 1
# mId = 1
# gId = 1

## Randomness
set.seed(2020)


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



###############################
### Gene Parameter Estimate ###
###############################

## Data Processing
classIds = which( classData$level2class == oneClass )
geneIds = geneSelec[,oneClass,sId,gId]
oneData = mrnaData[classIds,geneIds]

## Estimation
geneParam = fitHurdle( log(oneData+1), parallel=F )


##############
### Output ###
##############

## Output : gene parameter estimate
outNames$geneParam = paste(
	"../", "Real", "/", "A", "5",
	"/", "b", "/", "comp", "HurdleNormal",
	"_", oneClass, "_", sId, "_", mId, "_", gId,
	".rds", sep=""
)
saveRDS(geneParam, outNames$geneParam)
