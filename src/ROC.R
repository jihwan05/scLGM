ROCs = function(TrueGraph, EstmGraphs, regs, nrep=NULL) {
	if (dim(TrueGraph)[1] == dim(TrueGraph)[2]) {
		p = dim(TrueGraph)[1]
	} else {
		stop("unmatching p")
	}
	if (is.null(nrep)) {
		nrep = dim(EstmGraphs)[3]
	}
	nreg = dim(EstmGraphs)[4]
	errorMat = matrix(0, nreg, 9)
	colnames(errorMat) = c("regs",
			"nRP", "nTP", "nFN", "nRN", "nFP", "nTN",
			"TPR", "FPR")
	errorMat = data.frame(errorMat)
	TG = TrueGraph
	errorMat$nRP = sum(TG[upper.tri(TG)] != 0) * nrep
	errorMat$nRN = sum(TG[upper.tri(TG)] == 0) * nrep
	for (irep in 1:nrep) { for (ireg in 1:nreg) {
		EG = EstmGraphs[,,irep,ireg]
		errorMat$nTP[ireg] = errorMat$nTP[ireg] +
				sum((TG[upper.tri(TG)] != 0) & (EG[upper.tri(EG)] != 0))
		errorMat$nFN[ireg] = errorMat$nFN[ireg] +
				sum((TG[upper.tri(TG)] != 0) & (EG[upper.tri(EG)] == 0))
		errorMat$nFP[ireg] = errorMat$nFP[ireg] +
				sum((TG[upper.tri(TG)] == 0) & (EG[upper.tri(EG)] != 0))
		errorMat$nTN[ireg] = errorMat$nTN[ireg] +
				sum((TG[upper.tri(TG)] == 0) & (EG[upper.tri(EG)] == 0))
	}}
	errorMat = errorMat / nrep
	errorMat$TPR = errorMat$nTP / errorMat$nRP
	errorMat$FPR = errorMat$nFP / errorMat$nRN
	return(errorMat)
}

ROC1 = function(TrueGraph, EstmGraphs) {
	if (nrow(TrueGraph) == ncol(TrueGraph)) {
		p = nrow(TrueGraph)[1]
	} else {
		stop("unmatching p")
	}
	nrep = dim(EstmGraphs)[3]
	errorMat = matrix(0, nrep, 10)
	colnames(errorMat) = c("regs",
			"nRP", "nTP", "nFN", "nRN", "nFP", "nTN",
			"TPR", "FPR", "MCC")
	errorMat = data.frame(errorMat)
	TG = TrueGraph
	errorMat$nRP = sum(TG[upper.tri(TG)] != 0)
	errorMat$nRN = sum(TG[upper.tri(TG)] == 0)
	for (irep in 1:nrep) {
		EG = EstmGraphs[,,irep]
		errorMat$nTP[irep] = errorMat$nTP[irep] +
				sum((TG[upper.tri(TG)] != 0) & (EG[upper.tri(EG)] != 0))
		errorMat$nFN[irep] = errorMat$nFN[irep] +
				sum((TG[upper.tri(TG)] != 0) & (EG[upper.tri(EG)] == 0))
		errorMat$nFP[irep] = errorMat$nFP[irep] +
				sum((TG[upper.tri(TG)] == 0) & (EG[upper.tri(EG)] != 0))
		errorMat$nTN[irep] = errorMat$nTN[irep] +
				sum((TG[upper.tri(TG)] == 0) & (EG[upper.tri(EG)] == 0))
	}
	errorMat$TPR = errorMat$nTP / errorMat$nRP
	errorMat$FPR = errorMat$nFP / errorMat$nRN
	errorMat$MCC = errorMat$nTP*errorMat$nTN - errorMat$nFP*errorMat$nFN
	for (irep in 1:nrep) {
		if ( (errorMat$nTP[irep]+errorMat$nFP[irep]) == 0 ||
				(errorMat$nTP[irep]+errorMat$nFN[irep]) == 0 ||
				(errorMat$nTN[irep]+errorMat$nFP[irep]) == 0 ||
				(errorMat$nTN[irep]+errorMat$nFN[irep]) == 0 ) {
			errorMat$MCC[irep] = 0
		} else {
			errorMat$MCC[irep] = errorMat$MCC[irep] /
				sqrt((errorMat$nTP[irep]+errorMat$nFP[irep])*
						(errorMat$nTP[irep]+errorMat$nFN[irep])*
						(errorMat$nTN[irep]+errorMat$nFP[irep])*
						(errorMat$nTN[irep]+errorMat$nFN[irep]))
		}
	}
	return(errorMat)
}

ROCnew = function(TrueGraph, EstmGraph) {
	results = data.frame(
		nRP=0, nTP=0, nFN=0, nRN=0, nFP=0, nTN=0,
		TPR=0, FPR=0, MCC=0
	)
	TG = TrueGraph
	EG = EstmGraph
	results$nRP = sum( TG[upper.tri(TG)] != 0 )
	results$nRN = sum( TG[upper.tri(TG)] == 0 )
	results$nTP = sum(
		( TG[ upper.tri(TG) ] != 0 )
		& ( EG[ upper.tri(EG) ] != 0 )
	)
	results$nFN = sum(
		( TG[ upper.tri(TG) ] != 0 )
		& ( EG[ upper.tri(EG) ] == 0 )
	)
	results$nFP = sum(
		( TG[ upper.tri(TG) ] == 0 )
		& ( EG[ upper.tri(EG) ] != 0 )
	)
	results$nTN = sum(
		( TG[ upper.tri(TG) ] == 0 )
		& ( EG[ upper.tri(EG) ] == 0 )
	)
	results$TPR = results$nTP / results$nRP
	results$FPR = results$nFP / results$nRN
	results$MCC = results$nTP*results$nTN - results$nFP*results$nFN
	if (
		( results$nTP + results$nFP ) == 0 ||
		( results$nTP + results$nFN ) == 0 ||
		( results$nTN + results$nFP ) == 0 ||
		( results$nTN + results$nFN ) == 0
	) {
		results$MCC = 0
	} else {
		results$MCC = results$MCC / (
			sqrt( results$nTP + results$nFP ) *
			sqrt( results$nTP + results$nFN ) *
			sqrt( results$nTN + results$nFP ) *
			sqrt( results$nTN + results$nFN )
		)
	}
	return(results)
}