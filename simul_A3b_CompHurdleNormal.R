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
# ireg = 10

## Packages and Libraries
library("HurdleNormal")
source("./simul_A1a_data.R")

## Randomness
set.seed(2020)


#####################
# Reading True Data #
#####################

fileName = paste(
	"../", "Simul", "/", "A", "1",
	"/", "a", "_", "Data",
	"_", "sce", isce,
	".rds", sep=""
)
YObs = readRDS(fileName)[[irep]]$Y
rm(fileName)


##############################
# Estimating Gene Parameters #
##############################

##
Fits = fitHurdle( log(YObs+1), parallel=F )


####################
# Saving Estimates #
####################

fileName = paste(
	"../", "Simul", "/", "A", "3",
	"/", "b", "/", "Comp", "HurdleNormal",
	"_", "sce", isce, "rep", irep,
	".rds", sep=""
)
saveRDS(Fits, fileName)
rm(fileName)
