argVec = commandArgs(T)

if ( length(argVec) == 0 ) {
	warning("No arguments has been supplied.")
} else {
	print( paste(
		"Total",
		length(argVec),
		"arguments are supplied.",
		sep = " "
	) )
	for ( eachArg in argVec ) {
		print( paste("    ", eachArg, sep="") )
		eval( parse( text = eachArg ) )
	}
	rm(eachArg)
}
rm(argVec)