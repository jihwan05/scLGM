###############
# Preparation #
###############

# isce = c(11, 12, 13, 21, 22, 23, 31, 32, 33, 41, 42, 43, 77)
#isce = 77
 

###################
# Reading Results #
###################

##
fileName = paste(
	"../", "Simul", "/", "A", "4",
	"/", "a", "/", "Aggr", "Ours", "Glasso",
	"_", "sce", isce,
	".rds", sep=""
)
ROC_Ours_Glasso = readRDS(fileName)
rm(fileName)
Summar = ROC_Ours_Glasso$bests[,c(3,6,8,9,10)]
Summar[,1] = Summar[,1] + Summar[,2]
Summar = Summar[,-2]
round(colMeans(Summar), 3)
round(sqrt(diag(cov(Summar))), 3)

##
fileName = paste(
	"../", "Simul", "/", "A", "4",
	"/", "b", "/", "Aggr", "Ours", "Clime",
	"_", "sce", isce,
	".rds", sep=""
)
ROC_Ours_Clime = readRDS(fileName)
rm(fileName)
Summar = ROC_Ours_Clime$bests[,c(3,6,8,9,10)]
Summar[,1] = Summar[,1] + Summar[,2]
Summar = Summar[,-2]
round(colMeans(Summar), 3)
round(sqrt(diag(cov(Summar))), 3)

##
fileName = paste(
	"../", "Simul", "/", "A", "4",
	"/", "c", "/", "Aggr", "Comp", "Glasso",
	"_", "sce", isce,
	".rds", sep=""
)
ROC_HurdleNormal = readRDS(fileName)
rm(fileName)
Summar = ROC_HurdleNormal$bests[,c(3,6,8,9,10)]
Summar[,1] = Summar[,1] + Summar[,2]
Summar = Summar[,-2]
round(colMeans(Summar), 3)
round(sqrt(diag(cov(Summar))), 3)

##
fileName = paste(
	"../", "Simul", "/", "A", "4",
	"/", "d", "/", "Aggr", "Comp", "HurdleNormal",
	"_", "sce", isce,
	".rds", sep=""
)
ROC_Glasso = readRDS(fileName)
rm(fileName)
Summar = ROC_Glasso$bests[,c(3,6,8,9,10)]
Summar[,1] = Summar[,1] + Summar[,2]
Summar = Summar[,-2]
round(colMeans(Summar), 3)
round(sqrt(diag(cov(Summar))), 3)

##
RocOursGlasso = ROC_Ours_Glasso$bests[ , c(3,6,9,8,10) ]
RocOursClime = ROC_Ours_Clime$bests[ , c(3,6,9,8,10) ]
RocCompHurdleNormal = ROC_HurdleNormal$bests[ , c(3,6,9,8,10) ]
RocCompGlasso = ROC_Glasso$bests[ , c(3,6,9,8,10) ]

##
RocOursGlasso[,1] = RocOursGlasso[,1] + RocOursGlasso[,2]
RocOursGlasso = RocOursGlasso[,-2]
RocOursClime[,1] = RocOursClime[,1] + RocOursClime[,2]
RocOursClime = RocOursClime[,-2]
RocCompHurdleNormal[,1] = RocCompHurdleNormal[,1] + RocCompHurdleNormal[,2]
RocCompHurdleNormal = RocCompHurdleNormal[,-2]
RocCompGlasso[,1] = RocCompGlasso[,1] + RocCompGlasso[,2]
RocCompGlasso = RocCompGlasso[,-2]

##
Roc1 = rbind(
	scLGM_Glasso = colMeans( RocOursGlasso ),
	scLGM_Clime = colMeans( RocOursClime ),
	HurdlrNormal = colMeans( RocCompHurdleNormal ),
	Glasso = colMeans( RocCompGlasso )
)
Roc2 = rbind(
	scLGM_Glasso = diag( var( RocOursGlasso ) ),
	scLGM_Clime = diag( var( RocOursClime ) ),
	HurdlrNormal = diag( var( RocCompHurdleNormal ) ),
	Glasso = diag( var( RocCompGlasso ) )
)
Roc = round( cbind( Roc1, sqrt(Roc2) )[ , c(1,5,2,6,3,7,4,8) ], 3 )
print(Roc)

##
fileName = paste(
	"../", "Simul", "/", "A", "5",
	"/", "a", "/", "ROC",
	"_", "sce", isce,
	".rds", sep=""
)
saveRDS(Roc, fileName)
rm(fileName)


######################
# Drawing ROC Curves #
######################

##
ROC_Glasso$curve[ nrow(ROC_Glasso$curve):1, ] = ROC_Glasso$curve
if ( (isce%%10) !=3 ) {
	ROC_Ours_Glasso$curve[ nrow(ROC_Ours_Glasso$curve)+1, ] = 1
	ROC_Ours_Clime$curve[ nrow(ROC_Ours_Clime$curve)+1, ] = 1
	ROC_HurdleNormal$curve[ nrow(ROC_HurdleNormal$curve)+1, ] = 1
	ROC_Glasso$curve[ nrow(ROC_Glasso$curve)+1, ] = 1
}

##
quartz(width=3.5, height=3.5)
par(mar=0.25+c(3,3,0,0))
par(mgp=c(2,1,0))

##
plot(0,0, type="n",
		xlim=c(0,1), ylim=c(0,1), xlab="FPR", ylab="TPR")
lines(c(0,1), c(0,1), lty=1, lwd=0.25, col="lightgrey")

## Ours_Glasso
lines(ROC_Ours_Glasso$curve[,c(9,8)], lwd=0.75, lty=1, col="blue")
#points(ROCs$Ours_Glasso[,c(4,2)], cex=0.25, pch=1, col="blue")
points(x=colMeans(ROC_Ours_Glasso$bests)[9],
		y=colMeans(ROC_Ours_Glasso$bests)[8],
		cex=1, pch=1, col="blue")

## Ours_Clime
lines(ROC_Ours_Clime$curve[,c(9,8)], lwd=0.75, lty=1, col="red")
#points(ROCs$Ours_Clime[,c(4,2)], cex=0.25, pch=2, col="red")
points(x=colMeans(ROC_Ours_Clime$bests)[9],
		y=colMeans(ROC_Ours_Clime$bests)[8],
		cex=1, pch=2, col="red")

## HurdleNormal
lines(ROC_HurdleNormal$curve[,c(9,8)], lwd=0.75, lty=2, col="black")
#points(ROCs$HurdleNormal[,c(4,2)], cex=0.25, pch=4, col="black")
points(x=colMeans(ROC_HurdleNormal$bests)[9],
		y=colMeans(ROC_HurdleNormal$bests)[8],
		cex=1, pch=3, col="black")

## Glasso
lines(ROC_Glasso$curve[,c(9,8)], lwd=0.75, lty=3, col="darkgreen")
#points(ROCs$Glasso[,c(4,2)], cex=0.25, pch=3, col="darkgreen")
points(
	x=colMeans(ROC_Glasso$bests)[9],
	y=colMeans(ROC_Glasso$bests)[8],
	cex=1, pch=4, col="darkgreen"
)

## Legends
if ( isce==43 | isce==77 )  {
	legend(
		"bottomright", bty="n",
		c(
			"scLGM (Glasso)", "scLGM (Clime)",
			"HurdleNormal", "Glasso"
		),
		lwd=1, lty=c(1,1,2,3), cex=1, pch=1:4,
		col=c("blue", "red", "black", "darkgreen")
	)
}

################
# Saving Graph #
################

##
fileName = paste(
	"../", "Simul", "/", "A", "5",
	"/", "a", "/", "ROC",
	"_", "sce", isce,
	".pdf", sep=""
)
quartz.save(fileName, type="pdf")
rm(fileName)
