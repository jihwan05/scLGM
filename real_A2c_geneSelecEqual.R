#############
### Setup ###
#############

## Preparation
inNames = outNames = list()

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


######################
### Gene Selection ###
######################

## Selection : top 100 + random 100 genes
geneSelec = array( NA, c(200, length(specClass), 2, 100) )
for ( cId in 1:length(specClass) ) { # cId = 1
	##
	oneClass = specClass[cId]
	classIds = which( classData$level2class == oneClass )
	oneDataMrna = mrnaData[classIds,]
	##
	cm1 = colMeans(oneDataMrna)
	cm0 = colMeans( oneDataMrna > 0 )
	cm1s = cm0s = cbind(cm1, cm0)
	##
	cm1s = cm1s[ order( cm1s[,2], decreasing=T ), ]
	cm1s = cm1s[ order( cm1s[,1], decreasing=T ), ]
	##
	cm0s = cm0s[ order( cm0s[,1], decreasing=T ), ]
	cm0s = cm0s[ order( cm0s[,2], decreasing=T ), ]
	##
	cm1s = cm1s[ which( cm1s[,2] >= 0.5 ), ]
	cm0s = cm0s[ which( cm0s[,2] >= 0.5 ), ]
	##
	for (k in 1:100) {
		geneSelec[1:100,cId,1,k] = rownames(cm1s)[1:100]
		geneSelec[1:100,cId,2,k] = rownames(cm0s)[1:100]
	}
	##
	cm1s = cm1s[ -(1:100), ]
	cm0s = cm0s[ -(1:100), ]
	for (k in 1:100) {
		geneSelec[101:200,cId,1,k] = sample( rownames(cm1s), 100 )
		geneSelec[101:200,cId,2,k] = sample( rownames(cm0s), 100 )
	}
}
dimnames(geneSelec)[[2]] = specClass


##############
### Output ###
##############

## Output : selections of 200 genes
outNames$geneSelec = paste(
	"../", "Real", "/", "A", "2",
	"/", "c", "_", "geneSelec", "Equal",
	".rds", sep=""
)
saveRDS(geneSelec, outNames$geneSelec)
