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
source("./src/scLGM_geneParam.R")
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


###########################
# Reading Cell Parameters #
###########################

fileName = paste(
	"../", "Simul", "/", "A", "2",
	"/", "a", "/", "Ours", "Cell",
	"_", "sce", isce, "rep", irep,
	".rds", sep=""
)
CellParam = readRDS(fileName)
rm(fileName)


##############################
# Estimating Gene Parameters #
##############################

##
Fit = algMain(
	YObs, regsOurs_Glasso[ireg],
	CellParam$alpha, CellParam $beta,
	CellParam$kappa, CellParam $tau,
	Mstep="glasso", printll=T
)


####################
# Saving Estimates #
####################
fileName = paste(
	"../", "Simul", "/", "A", "2",
	"/", "b", "/", "Ours", "Gene", "Glasso",
	"_", "sce", isce, "rep", irep, "reg", ireg,
	".rds", sep=""
)
saveRDS(Fit, fileName)
rm(fileName)
