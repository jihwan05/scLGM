#############
### Setup ###
#############

## Preparation
inNames = outNames = list()
#setwd("~/project/Jihwan/scLGM/Code")

## Packages and Libraries
suppressMessages( library("org.Mm.eg.db") )
source("./src/KEGG.R")

## Command Line Arguments
source("./src/_commArgs.R")
# cId = 1
# tId = 1
# k = 1


#############
### Input ###
#############

## Input : list of special classes : big three
inNames$specClass = paste("../Real/A1/b_specClass.rds")
specClass = readRDS(inNames$specClass)

## Input : selections of 200 genes
inNames$geneSelec = paste(
	"../", "Real", "/", "A", "2",
	"/", "c", "_", "geneSelec", "Equal",
	".rds", sep=""
)
geneSelec = readRDS(inNames$geneSelec)


############################################
### Graph Information from KEGG Database ###
############################################

oneClass = specClass[cId]
names = geneSelec[,oneClass,tId,k]
eIds = AnnotationDbi::select(
	org.Mm.eg.db, keys=names,
	columns='ENTREZID', keytype='SYMBOL'
)
keggGraph = ENT2Graph( eIds[,2], edge=TRUE )


##############
### Output ###
##############

## Output : graph information from KEGG database
outNames$keggGraph = paste(
	"../", "Real", "/", "A", "3",
	"/", "c", "/", "keggGraph", "Equal",
	"_", oneClass, "_", tId, "_", k,
	".rds", sep=""
)
saveRDS(keggGraph, outNames$keggGraph)
