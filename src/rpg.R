rpg.changgee <- function(r, n, z, trunc=2000) {
	# r : number of observations
	# (n, z) : two parameters
	gam = matrix(rgamma(r*trunc,n),r)
	co = t(matrix((1:trunc-0.5)^2,trunc,r))*2*pi^2 + z^2/2
	return(apply(gam/co,1,sum))
}
