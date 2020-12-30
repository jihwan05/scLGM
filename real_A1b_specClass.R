#############
### Setup ###
#############

## Preparation
inNames = outNames = list()


##############
### Output ###
##############

## Output : list of special classes : big three
specClass = c("CA1Pyr2", "CA1Pyr1", "Oligo6", "Oligo5", "Oligo4", "Vend2")
outNames$specClass = paste("../Real/A1/b_specClass.rds")
saveRDS(specClass, outNames$specClass)
