#############
### Setup ###
#############

## Preparation
inNames = outNames = list()


##############
### Output ###
##############

## Output : list of special classes : big three
specClass = c(
	"pyramidal CA1", "oligodendrocytes", "pyramidal SS",
	"interneurons", "endothelial-mural",
	"astrocytes_ependymal", "microglia"
)
outNames$specClass = paste(
	"../", "Real", "/", "B", "1",
	"/", "b", "_", "specClass",
	".rds", sep=""
)
saveRDS(specClass, outNames$specClass)
