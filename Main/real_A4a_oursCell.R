#############
### Setup ###
#############

## Preparation
inNames = outNames = list()
setwd("~/project/Jihwan/scLGM/Code")

## Packages and Libraries
#suppressMessages( library("org.Mm.eg.db") )
source("./src/scLGM_cellParam.R")

## Command Line Arguments
source("./src/_commArgs.R")
# cId = 1

## Randomness
set.seed(2020)


###############################
### Cell Parameter Estimate ###
###############################

## Input : class data
inNames$classData = paste("../Real/A0/a_classData.rds")
classData = readRDS(inNames$classData)

## Input : detected ERCC spike-in UMI count data
inNames$dspikeData = paste("../Real/A0/a_dspikeData.rds")
dspikeData = readRDS(inNames$dspikeData)

## Input : original ERCC spike-in molecules per cell data
inNames$ospikeData = paste("../Real/A0/b_ospikeData.rds")
ospikeData = readRDS(inNames$ospikeData)

## Input : list of special classes : big three
inNames$specClass = paste("../Real/A1/b_specClass.rds")
specClass = readRDS(inNames$specClass)

## Name Matching
oneClass = specClass[cId]
classIds = which( classData$level2class == oneClass )
oneData = dspikeData[classIds,]
cns = substr( colnames(oneData), 6, 15 )
ids = match( cns, ospikeData$id )
ospikeSub = ospikeData[ids,]
print( sum( cns != ospikeSub$id ) )
print( sum( cns == ospikeSub$id ) )

## MCMC
cellParam = scLGM_cellParam(ospikeSub$mpcNew, t(oneData), 3000, 2000)
#rownames(estim) = colnames(Output)


##############
### Output ###
##############

## Output : cell parameter estimate
outNames$cellParam = paste(
	"../", "Real", "/", "A", "4",
	"/", "a", "_", "cellParam",
	"_", oneClass,
	".rds", sep=""
)
saveRDS(cellParam, outNames$cellParam)
