library("KEGGgraph")
library("KEGGREST")
library("graph")

##
KEGG2Group <- function(KEGG.ID)
{
	nKegg = length(KEGG.ID)
	PATH.ID = c()
	for ( i in 0:(nKegg%/%100) )
	{
		print(i)
		PATH.ID = c(PATH.ID,keggLink("pathway",KEGG.ID[100*i+1:100]))
	}
	PATH.ID = unique(PATH.ID)
	nPath = length(PATH.ID)
	Pathway = matrix(FALSE,nKegg,nPath)
	for ( i in 0:(nKegg%/%100) )
	{
		print(i)
		PATH.ID = c(PATH.ID,keggLink("pathway",KEGG.ID[100*i+1:100]))
	}
}


##
KEGG2Graph <- function(KEGG.ID,edge=TRUE)
{
	nKegg = length(KEGG.ID)
	PATH.ID = c()
	for ( i in 0:(nKegg%/%100) )
	{
		print(i)
		PATH.ID = c(PATH.ID,keggLink("pathway",KEGG.ID[100*i+1:100]))
	}
	PATH.ID = unique(PATH.ID)
	nPath = length(PATH.ID)
	print(nPath)
	KEGG.PATH = matrix(FALSE,nKegg,nPath)
	KEGG.EDGE = matrix(0,0,2)
	for ( i in 1:nPath )
	{
		print(i)
		kgml = keggGet(PATH.ID[i],"kgml")
		
		graph = parseKGML2Graph(kgml)
		nds = match(nodes(graph),KEGG.ID)
		nds = nds[!is.na(nds)]
		KEGG.PATH[nds,i] = TRUE
		if ( edge )
			for ( node in nds )
			{
				adju = unique(adj(graph,KEGG.ID[node])[[1]])
				adjs = match(adju,KEGG.ID)
				adjs = adjs[!is.na(adjs)]
				if ( length(adjs) > 0 )
					KEGG.EDGE = rbind(KEGG.EDGE,cbind(node,adjs))
			}
	}
	KEGG.PATH = matrix(KEGG.PATH[,apply(KEGG.PATH,2,sum)>1],nKegg)
	ret = list(PATH=KEGG.PATH)
	if ( edge )
	{
		KEGG.EDGE.D = unique(KEGG.EDGE[KEGG.EDGE[,1]!=KEGG.EDGE[,2],])
		ret = c(ret,list(E=unique(rbind(KEGG.EDGE.D,KEGG.EDGE.D[,c(2,1)]))))
	}
	ret
}

##
ENT2KEGG <- function(ENTREZ.ID)
{
	nEnt = length(ENTREZ.ID)
	
	nKegg = 0
	KEGG.ID = c()
	
	KEGG2ENT = c()
	ENT2KEGG = rep(0,nEnt)
	
	print(nEnt)
	for ( i in 1:nEnt )
	{
		print(i)
		kid = keggConv("genes",paste("ncbi-geneid",ENTREZ.ID[i],sep=":"))
		if ( length(kid) > 0 )
		{
			nKegg = nKegg + 1
			KEGG.ID[nKegg] = kid
			KEGG2ENT[nKegg] = i
			ENT2KEGG[i] = nKegg
		}
	}
	list(ENTREZ.ID=ENTREZ.ID,KEGG.ID=KEGG.ID,KEGG2ENT=KEGG2ENT,ENT2KEGG=ENT2KEGG)
}

##
ENT2Graph <- function(ENTREZ.ID,edge=TRUE)
{
	map = ENT2KEGG(ENTREZ.ID)
	KEGG.Graph = KEGG2Graph(map$KEGG.ID,edge)
	
	nPATH = ncol(KEGG.Graph$PATH)
	print(nPATH)
	nEnt = length(ENTREZ.ID)
	PATH = matrix(FALSE,nEnt,nPATH)
	if ( nPATH > 0 )
		for ( i in 1:nPATH )
			PATH[map$KEGG2ENT[which(KEGG.Graph$PATH[,i])],i] = TRUE
	ret = list(PATH=PATH)
	if ( edge )
	{
		E = matrix(map$KEGG2ENT[KEGG.Graph$E],ncol=2)
		ret = c(ret,list(E=E))
	}
	ret
}

