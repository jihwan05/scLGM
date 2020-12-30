#############
### Setup ###
#############

## Preparation
inNames = outNames = list()

## Packages and Libraries
library(xtable)


#############
### Input ###
#############

## Input : class data
inNames$classData = paste("../Real/A0/a_classData.rds")
classData = readRDS(inNames$classData)


#######################
### Data Processing ###
#######################

## Table : cell counts for each class
classTable = table(classData)
classTable = classTable[ names( sort( table(classData$level1class), T ) ), ]
classTable = classTable[ , names( sort( table(classData$level2class), T ) ) ]
classTable = cbind( classTable, total=rowSums(classTable) )
classTable = rbind( classTable, total=colSums(classTable) )


##############
### Output ###
##############

## Output : cell counts for each class
print( t(classTable) )
xtable( t( classTable ), digits=0 )
