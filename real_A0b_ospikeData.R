#############
### Setup ###
#############

## Preparation
inNames = outNames = list()


#############
### Input ###
#############

## Input : observed ERCC spike-in UMI count data
inNames$Spike = paste("../jia2017accounting/ERCC_controls.csv")
SpikeRaw = read.csv(inNames$Spike)


#######################
### Data Processing ###
#######################

ospikeData = data.frame(
	id = substr( SpikeRaw$ERCC.ID, 6, 20 ),
	mpcOrig = SpikeRaw$MoleculesPerCell,
	mpcNew = SpikeRaw$CorrectMPC
)


##############
### Output ###
##############

## Output : original ERCC spike-in molecules per cell data
outNames$ospikeData = paste("../Real/A0/b_ospikeData.rds")
saveRDS(ospikeData, outNames$ospikeData)
