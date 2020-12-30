#############
### Setup ###
#############

## Preparation
inNames = outNames = list()


#############
### Input ###
#############

## Input : given mRNA matrix
inNames$Mrna = paste("../zeisel2015cell/expression_mRNA_17-Aug-2014.txt")
MrnaRaw = read.table(inNames$Mrna, sep="\t", stringsAsFactors=F, skip=7)
Mrna = MrnaRaw
Mrna = Mrna[-4,]
Mrna[-(1:3),2] = Mrna[-(1:3),1]
Mrna = Mrna[,-1]
rownames(Mrna) = Mrna[,1]
colnames(Mrna) = Mrna[1,]
Mrna = t( Mrna[-1,-1] )

## Input : ERCC Spike-Ins
inNames$Spike = paste("../zeisel2015cell/expression_spikes_17-Aug-2014.txt")
SpikeRaw = read.table(inNames$Spike, sep="\t", stringsAsFactors=F, skip=7)
Spike = SpikeRaw
Spike = Spike[-4,]
Spike[-(1:3),2] = Spike[-(1:3),1]
Spike = Spike[,-1]
rownames(Spike) = Spike[,1]
colnames(Spike) = Spike[1,]
Spike = t( Spike[-1,-1] )


#######################
### Data Processing ###
#######################

## Data : all class information
check = cbind( ( Mrna[,1:2] != Spike[,1:2] ), ( Mrna[,1:2] == Spike[,1:2] ) )
print( colSums(check) )
classData = data.frame( Mrna[,1:2] )

## Data : all UMI counts
mrnaData = matrix( NA, nrow(Mrna), ncol(Mrna)-2 )
for ( i in 1:nrow(mrnaData) ) {
	mrnaData[i,] = as.numeric( Mrna[ i, -(1:2) ] )
}
rownames(mrnaData) = rownames(Mrna)
colnames(mrnaData) = colnames(Mrna)[-(1:2)]
mrnaData = data.frame(mrnaData)

## Data : all ERCC spike-in counts
dspikeData = matrix( NA, nrow(Spike), ncol(Spike)-2 )
for ( i in 1:nrow(dspikeData) ) {
	dspikeData[i,] = as.numeric( Spike[ i, -(1:2) ] )
}
rownames(dspikeData) = rownames(Spike)
colnames(dspikeData) = colnames(Spike)[-(1:2)]
dspikeData = data.frame(dspikeData)


##############
### Output ###
##############

## Output : raw data
rawData = list( mrna=MrnaRaw, spike=SpikeRaw )
outNames$rawData = paste("../Real/A0/a_rawData.rds" )
saveRDS(rawData, outNames$rawData)

## Output : class data
outNames$classData = paste("../Real/A0/a_classData.rds")
saveRDS(classData, outNames$classData)

## Output : mRNA UMI count data
outNames$mrnaData = paste("../Real/A0/a_mrnaData.rds")
saveRDS(mrnaData, outNames$mrnaData)

## Output : detected ERCC spike-in UMI count data
outNames$dspikeData = paste("../Real/A0/a_dspikeData.rds")
saveRDS(dspikeData, outNames$dspikeData)
