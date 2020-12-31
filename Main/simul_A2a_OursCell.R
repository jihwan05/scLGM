#############
### Setup ###
#############

## Preparation
inNames = outNames = list()
setwd("~/project/Jihwan/scLGM/Code")

## Command Line Arguments
source("./src/_commArgs.R")
# isce = 11
# irep = 1

## Packages and Libraries
source("./src/scLGM_cellParam.R")

## Randomness
set.seed(2020)


##############################
# Estimating Cell Parameters #
##############################

##
fileName = paste(
	"../", "Simul", "/", "A", "1",
	"/", "a", "_", "SpikeIns",
	"_", "sce", isce,
	".rds", sep=""
)
SpikeIns = readRDS(fileName)
rm(fileName)

##
CellParam = scLGM_cellParam(
	SpikeIns$input,
	SpikeIns$output[[irep]],
	3000, 2000
)


###################################
# Saving Cell Parameter Estimates #
###################################

fileName = paste(
	"../", "Simul", "/", "A", "2",
	"/", "a", "/", "Ours", "Cell",
	"_", "sce", isce, "rep", irep,
	".rds", sep=""
)
saveRDS(CellParam, fileName)
rm(fileName)
