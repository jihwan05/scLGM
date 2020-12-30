###################
# Setting: Common #
###################

if ( (isce %/% 10) == 7 ) {
	n=447
	p=194
} else {
	n = 100
	if ( (isce %% 10) == 0 ) {
		p = 5
	} else if ( (isce %% 10) == 1 ) {
		p = 50
	} else if ( (isce %% 10) == 2 ) {
		p = 100
	} else if ( (isce %% 10) == 3 ) {
		p = 300
	}
}


###################
# Setting: Common #
###################

nrep = 100
nreg = 25
if ( (isce %/% 10) == 7 ) {
	regsOurs_Glasso = exp(seq(log(5), log(0.01), length.out=nreg))
	regsOurs_Clime = exp(seq(log(0.5), log(0.01), length.out=nreg))
	regsGlasso = exp(seq(log(10), log(0.05), length.out=nreg))
} else {
	regsOurs_Glasso = exp(seq(log(5), log(0.01), length.out=nreg))
	regsOurs_Clime = exp(seq(log(0.5), log(0.01), length.out=nreg))
	regsGlasso = exp(seq(log(20), log(0.1), length.out=nreg))
}
