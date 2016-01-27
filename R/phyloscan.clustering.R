#' @export
hivc.clu.brdist.stats.subtree<- function(node, tree, distmat, eval.dist.btw="leaf", stat.fun=max)
{
	require(ape)
	require(igraph)
	require(geiger)	
	if(eval.dist.btw=="leaf")
	{
		nlist	<- tips(tree,node)
		foo 	<- distmat[nlist,nlist] 		
	}
	else if(eval.dist.btw=="all")
	{
		nlist	<- tips(tree,node)
		elist 	<- tree$edge[which.edge(tree,nlist),2]
		foo 	<- distmat[elist,elist] 	
	}
	else	
		stop("unrecognized eval.dist.btw")
	return( stat.fun(foo[upper.tri(foo,diag=FALSE)]) )
}
######################################################################################
#' Compute a statistic of patristic distances for each subtree in the tree
#' @export
#' @import igraph geiger ape phytools
hivc.clu.brdist.stats<- function(tree, distmat=NULL, eval.dist.btw="leaf", stat.fun=max)
{
	require(phytools)
	if(is.null(distmat))
	{
		if(eval.dist.btw=='leaf')
			distmat	<- cophenetic.phylo(tree)
		else if(eval.dist.btw=='all')
			distmat <- dist.nodes(tree)
		else	stop("unrecognized eval.dist.btw")
	}
	ntips			<- Ntip(tree)
	nint 			<- tree$Nnode 		## number of internal nodes
	return(sapply(seq.int(ntips+1,ntips+nint), function(x) 		hivc.clu.brdist.stats.subtree(x,tree,distmat,eval.dist.btw=eval.dist.btw, stat.fun=stat.fun)		))	
}
######################################################################################
#' @useDynLib phyloscan
#' @export
hivc.clu.min.transmission.cascade<- function(brlmat)
{
	if(is.matrix(brlmat))
		ans	<- .Call("hivc_clu_mintransmissioncascade", brlmat[upper.tri(brlmat)])  						
	else
		ans	<- .Call("hivc_clu_mintransmissioncascade", brlmat)
	ans		
}
######################################################################################
#' @export
hivc.clu.clusterbythresh<- function(ph,thresh.brl=NULL,dist.brl=NULL,thresh.nodesupport=NULL,nodesupport=NULL,retval="all")
{
	require(ape)
	require(igraph)
	require(geiger)	
	if(all( is.null(c(thresh.brl, thresh.nodesupport))) ) 		stop("all threshold criteria NULL - provide at least one")
	if(!is.null(thresh.nodesupport) && is.null(nodesupport))	stop("node support threshold set but no node support values provided")
	if(!is.null(thresh.brl) && is.null(dist.brl))				stop("branch length threshold set but no branch length distances provided")
	if(!is.null(nodesupport) && any(nodesupport>1+EPS))			warning("Found nodesupport values above 1")
	if(is.null(thresh.brl))
	{		
		dist.brl			<- rep(1,Nnode(ph))
		thresh.brl			<- 1
	}	
	if(is.null(thresh.nodesupport))
	{
		nodesupport			<- rep(1,Nnode(ph))
		thresh.nodesupport	<- 1		
	}
#print(nodesupport[1:10]); print(thresh.nodesupport)
	## set up clustering
	ntips		<- Ntip(ph)	
	clu.i 		<- 0 ## cluster number
	clu.mem 	<- rep(NA,ntips+ph$Nnode) ## cluster member assignment
	clu.idx		<- rep(NA,ntips+ph$Nnode) ## cluster index assignment
	igraph.ph	<- graph.edgelist(ph$edge) ## ph in igraph form
	dfs 		<- graph.dfs(igraph.ph,root=ntips+1,neimode='out',order=TRUE,dist=TRUE)
#print(c(length(dfs$order),Ntip(ph),Nnode(ph)))
	## travese the ph in depth first order
	for(i in 1:length(dfs$order))
	{
		node <- dfs$order[i]
#print( c(node,node-ntips, is.na(clu.mem[node]), thresh.brl, dist.brl[node-ntips], dist.brl[node-ntips]<=thresh.brl, thresh.nodesupport, nodesupport[node-ntips], nodesupport[node-ntips]<thresh.nodesupport) )
		if(	node > ntips	&&											## skip leaves
				is.na(clu.mem[node]) &&										## only consider unassigned nodes
				nodesupport[node-ntips]>=thresh.nodesupport-EPS &&
				dist.brl[node-ntips]<=thresh.brl+EPS
				)	 
		{				
#print(nodesupport[node-ntips]);print(c(node,i))			
			clu.i 			<- clu.i+1
			subtree 		<- graph.dfs(igraph.ph,node,neimode='out',unreachable=FALSE)$order
			subtree 		<- subtree[! is.nan(subtree)]
			clu.mem[subtree]<- clu.i
			clu.idx[node]	<- clu.i
		}
	}
	ans <- list( clu.mem=clu.mem, clu.idx=which(!is.na(clu.idx)), size.all=table(clu.mem), size.tips=table(clu.mem[1:ntips]), ntips=ntips, thresh.brl=thresh.brl, thresh.nodesupport=thresh.nodesupport)
	if(retval=="tips")
		ans$clu.mem 		<- ans$clu.mem[1:ntips]
	return(ans)
}