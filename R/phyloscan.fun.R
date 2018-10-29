#' @export
#' @import data.table
#' @title Determine stage 1 phyloscanner runs
#' @description This function groups individuals for phyloscanner analyses, so that phylogenetic linkage between every pair of individuals is assessed at least once.
#'   Specifically, individuals are grouped into batches of specified size, and then, all possible pairs of batches are formed. 
#'   Each of these pairs of batches defines a group of individuals between whom phylogenetic linkages are assessed in one phyloscanner run.
#'   The number of individuals in each group is twice the batch size.
#' @param x Character vector of individual identifiers. 
#' @param batch.size Batch size. Default is 50.
#' @return data.table with rows 'IND' (individual identifiers), 'PTY_RUN' group for phyloscanner analysis, and 'BATCH' batch of individuals (not used further, but there should be two batches of individuals in each phyloscanner analysis).
phsc.define.stage1.analyses<- function(x, batch.size=50)	
{
	dind		<- data.table(IND=x)
	#	assign batches
	set(dind, NULL, 'BATCH', dind[, 1L+floor((seq_len(nrow(dind))-1)/batch.size)])
	batches		<- dind[, unique(BATCH)]
	#	assign phyloscanner runs
	pty.runs	<- as.data.table(t(combn(batches,2)))
	setnames(pty.runs, c('V1','V2'), c('BATCH','BATCH2'))	
	pty.runs[, PTY_RUN:= seq_len(nrow(pty.runs))]
	#	merge individuals to phyloscanner runs based on batches
	tmp			<- merge(pty.runs, dind, by='BATCH', allow.cartesian=TRUE)
	setnames(dind, c('IND','BATCH'), c('IND2','BATCH2'))
	tmp2		<- merge(pty.runs, dind, by='BATCH2', allow.cartesian=TRUE)
	setnames(dind, c('IND2','BATCH2'), c('IND','BATCH'))
	setnames(tmp2, 'IND2', 'IND')
	pty.runs	<- rbind(tmp, tmp2)
	pty.runs
}

#' @export
#' @import data.table
#' @importFrom igraph graph.data.frame clusters
#' @title Find phylogenetic transmission networks
#' @param rtp Pairs of individuals between whom linkage is not excluded, stored as data.table.
#' @param rplkl Summary of phylogenetic relationship counts for each pair, stored as data.table.
#' @param conf.cut Threshold on the proportion of deep-sequence phylogenies with distant/disconnected subgraphs above which pairs are considered phylogenetically unlinked. Default: 0.6
#' @param neff.cut Threshold on the minimum number of deep-sequence phylogenies with sufficient reads from two individuals to make any phylogenetic inferences. Default: 3.
#' @param verbose Flag to switch on/off verbose mode. Default: TRUE. 
#' @return list of two R objects 'transmission.networks', 'most.likely.transmission.chains'. See description.
#' @description This function reconstructs phylogenetic transmission networks from pairs of individuals between whom linkage is not excluded. 
#' Two R objects are generated: 
#' 	  'transmission.networks' is a data.table describing transmission networks, with information of phylogenetic support for transmission in direction 12, direction 21, linkage with direction unclear, and no linkage. 
#' 	  'most.likely.transmission.chains' is a data.table describing the most likely transmission chain for each transmission network.
phsc.find.transmission.networks.from.linked.pairs<- function(rtp, rplkl, conf.cut=conf.cut, neff.cut=neff.cut, verbose=TRUE)
{
	#	internal variables
	linked.group	<- 'TYPE_CHAIN_TODI'
	linked.no		<- 'distant'
	linked.yes		<- 'chain'
	scores.group	<- 'TYPE_ADJ_NETWORK_SCORES'
	scores.no		<- c('ambiguous','not close/disconnected') 
	dir.group		<- 'TYPE_ADJ_DIR_TODI2'
	
	#
	#	construct max prob network among all possible pairs
	rtn		<- subset(rtp, GROUP==linked.group & TYPE==linked.yes & NEFF>neff.cut)
	
	#	define potential transmission networks
	if(verbose) cat('\nReconstruct transmission networks among linked pairs, n=',nrow(rtn))
	tmp		<- subset(rtn, select=c(ID1, ID2))			
	tmp		<- graph.data.frame(tmp, directed=FALSE, vertices=NULL)
	rtc		<- data.table(ID=V(tmp)$name, CLU=clusters(tmp, mode="weak")$membership)	
	tmp2	<- rtc[, list(CLU_SIZE=length(ID)), by='CLU']
	setkey(tmp2, CLU_SIZE)
	tmp2[, IDCLU:=rev(seq_len(nrow(tmp2)))]
	rtc			<- subset( merge(rtc, tmp2, by='CLU'))
	rtc[, CLU:=NULL]
	setkey(rtc, IDCLU)	
	#	add info on edges: linked vs not linked
	setnames(rtc, c('ID'), c('ID1'))
	rtn		<- merge(rtn, rtc, by='ID1')
	rtc[, CLU_SIZE:=NULL]
	setnames(rtc, c('ID1'), c('ID2'))
	rtn		<- merge(rtn, rtc, by=c('ID2','IDCLU'))
	#	add info on edges: direction 12, direction 21, direction ambiguous, unlinked 
	tmp		<- subset(rtn, select=c(ID1, ID2, PTY_RUN, IDCLU, CLU_SIZE))	
	tmp2	<- c('ID1','ID2','PTY_RUN')
	tmp		<- merge(tmp, subset(rplkl, GROUP==scores.group), by=tmp2)
	tmp		<- merge(tmp, tmp[, list(TYPE=TYPE, POSTERIOR_SCORE=(POSTERIOR_ALPHA-1)/sum(POSTERIOR_ALPHA-1) ), by=c('ID1','ID2','PTY_RUN')], by=c('ID1','ID2','PTY_RUN','TYPE'))
	rtn		<- rbind(tmp, rtn)
	if(verbose) cat('\nFound transmission networks, n=',rtn[, length(unique(IDCLU))], '. Number of links (either direction and ambiguous)=', nrow(subset(rtn, GROUP==scores.group & TYPE=='not close/disconnected')), '. Number of individuals=', length(unique(c(rtn$ID1, rtn$ID2))),'.')
	
	#
	#	generate most likely transmission chains
	if(verbose) cat('\nReconstruct most likely transmission chains...')
	rtn		<- subset(rtn, GROUP==scores.group)
	rtnn	<- Phyloscanner.R.utilities:::phsc.get.most.likely.transmission.chains(rtn, verbose=0)	
	#	for TYPE=='ambiguous', this has the cols:
	#	POSTERIOR_SCORE 	posterior prob direction ambiguous before self-consistence
	#	MX_PROB_12			total posterior prob supporting  1 to 2 including 50% ambiguous AFTER self-consistence
	#	MX_PROB_21			total posterior prob supporting  2 to 1 including 50% ambiguous AFTER self-consistence
	#	MX_KEFF_21 			total KEFF supporting  2 to 1 including 50% ambiguous before self-consistence 
	#	MX_KEFF_12 			total KEFF supporting  1 to 2 including 50% ambiguous before self-consistence  
	#	LINK_12 			if there is a directed edge from 1 to 2 in max edge credibility network
	#	LINK_21				if there is a directed edge from 2 to 1 in max edge credibility network
	#	where self-consistence means that 12 xor 21 are set to zero	
	rtnn	<- subset(rtnn, TYPE=='ambiguous', select=c(ID1, ID2, PTY_RUN, IDCLU, POSTERIOR_SCORE, MX_PROB_12, MX_PROB_21, MX_KEFF_21, MX_KEFF_12, LINK_12, LINK_21))
	#	work out prob for linkage in max prob network, when 'inconsistent direction' is ignored
	rtnn[, POSTERIOR_SCORE_LINKED_MECN:= pmax(MX_PROB_12,MX_PROB_21) + 0.5*POSTERIOR_SCORE]
	set(rtnn, NULL, c('POSTERIOR_SCORE','MX_PROB_12','MX_PROB_21'), NULL)
	#	merge POSTERIOR_SCORE_LINKED on max prob network
	#	rationale: this describes prob of linkage. here, any 'inconsistent direction' is still considered as prob for linkage 	
	tmp		<- subset(rplkl, GROUP==linked.group & TYPE==linked.yes, c(ID1,ID2,PTY_RUN,POSTERIOR_ALPHA,POSTERIOR_BETA,N_TYPE))
	tmp[, POSTERIOR_SCORE_LINKED:= (POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE)]
	set(tmp, NULL, c('POSTERIOR_ALPHA','POSTERIOR_BETA','N_TYPE'), NULL)
	rtnn	<- merge(rtnn, tmp, by=c('ID1','ID2','PTY_RUN'), all.x=TRUE)
	#	merge POSTERIOR_SCORE_12 POSTERIOR_SCORE_21 (direction) on max prob network
	#	this is considering in denominator 12 + 21 before reducing probs to achieve self-consistency
	#	rationale: decide on evidence for direction based on comparing only the flows in either direction, 12 vs 21
	tmp		<- subset(rplkl, GROUP==dir.group, c(ID1,ID2,PTY_RUN,TYPE,POSTERIOR_ALPHA,POSTERIOR_BETA,N_TYPE))
	tmp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE)]
	set(tmp, NULL, c('POSTERIOR_ALPHA','POSTERIOR_BETA','N_TYPE'), NULL)
	set(tmp, NULL, 'TYPE', tmp[, paste0('POSTERIOR_SCORE_',TYPE)])
	tmp		<- dcast.data.table(tmp, ID1+ID2+PTY_RUN~TYPE, value.var='POSTERIOR_SCORE')
	rtnn	<- merge(rtnn, tmp, by=c('ID1','ID2','PTY_RUN'), all.x=TRUE)
	#	merge NETWORK_SCORE_12 NETWORK_SCORE_21 on max prob network
	#	this is considering in denominator 12 + 21 + unclear reducing probs to achieve self-consistency
	#	same as MX_PROB_12, MX_PROB_21, after the final step below that sets one of the two probs to zero
	tmp		<- subset(rplkl, GROUP==scores.group, c(ID1,ID2,PTY_RUN,TYPE,POSTERIOR_ALPHA))
	tmp		<- tmp[, list(TYPE=TYPE, POSTERIOR_SCORE=(POSTERIOR_ALPHA-1)/sum(POSTERIOR_ALPHA-1)), by=c('ID1','ID2','PTY_RUN')]
	tmp		<- subset(tmp, !TYPE%in%scores.no)
	set(tmp, NULL, 'TYPE', tmp[, paste0('NETWORK_SCORE_',TYPE)])	
	tmp		<- dcast.data.table(tmp, ID1+ID2+PTY_RUN~TYPE, value.var='POSTERIOR_SCORE')
	rtnn	<- merge(rtnn, tmp, by=c('ID1','ID2','PTY_RUN'), all.x=TRUE)
	#	ensure DIR scores and NETWORK_SCORE scores are compatible with self-consistency in maxprobnetwork
	tmp		<- rtnn[, which(LINK_12==0 & LINK_21==1 & POSTERIOR_SCORE_12>POSTERIOR_SCORE_21)]
	set(rtnn, tmp, c('POSTERIOR_SCORE_12','NETWORK_SCORE_12'), 0)
	tmp		<- rtnn[, which(LINK_12==1 & LINK_21==0 & POSTERIOR_SCORE_21>POSTERIOR_SCORE_12)]
	set(rtnn, tmp, c('POSTERIOR_SCORE_21','NETWORK_SCORE_21'), 0)	
	rtnn	<- subset(rtnn, LINK_12==1 | LINK_21==1)	
	if(verbose) cat('\nFound most likely transmission chains, n=',rtnn[, length(unique(IDCLU))], '. Number of links=', nrow(rtnn), '. Number of individuals=', length(unique(c(rtnn$ID1, rtnn$ID2))),'.')
	if(verbose) cat('\nDone.')
	# return
	list(transmission.networks=rtn, most.likely.transmission.chains=rtnn)
}


#' @export
#' @import data.table 
#' @title Find pairs of individuals between whom linkage is not excluded phylogenetically
#' @param indir Full directory path to output of phyloscanner runs
#' @param batch.regex Regular expression that identifies the batch ID of multiple phyloscanner analyses. Default: '^ptyr([0-9]+)_.*'.
#' @param conf.cut Threshold on the proportion of deep-sequence phylogenies with distant/disconnected subgraphs above which pairs are considered phylogenetically unlinked. Default: 0.6
#' @param neff.cut Threshold on the minimum number of deep-sequence phylogenies with sufficient reads from two individuals to make any phylogenetic inferences. Default: 3.
#' @param verbose Flag to switch on/off verbose mode. Default: TRUE. 
#' @param dmeta Optional individual-level meta-data that is to be added to output. Can be NULL.
#' @return list of three R objects 'linked.pairs', 'relationship.counts', 'windows'. See description.
#' @description This function identifies pairs of individuals between whom linkage is not excluded phylogenetically in a large number of phyloscanner analyses, and provides detailed information on them.
#' Three R objects are generated: 
#' 	  'linked.pairs' is a data.table that describes pairs of individuals between whom linkage is not excluded phylogenetically.
#' 	  'relationship.counts' is a data.table that summarises the phylogenetic relationship counts for each pair. 
#' 	  'windows' is a data.table that describes the basic phyloscanner statistics (distance, adjacency, paths between subgraphs) in each deep-sequence phylogeny for each pair. 
phsc.find.linked.pairs<- function(indir, batch.regex='^ptyr([0-9]+)_.*', conf.cut=0.6, neff.cut=3, verbose=TRUE, dmeta=NULL)
{
	#	internal variables
	linked.group	<- 'TYPE_CHAIN_TODI'
	linked.no		<- 'distant'
	linked.yes		<- 'chain'
	#
	#	from every phyloscanner run,  
	#	select pairs that are not predominantly unlinked by distance + topology
	infiles	<- data.table(F=list.files(indir, pattern='pairwise_relationships.rda', full.names=TRUE))
	infiles[, PTY_RUN:= as.integer(gsub(batch.regex,'\\1',basename(F)))]
	setkey(infiles, PTY_RUN)	
	if(verbose) cat('\nFound phylogenetic relationship files, n=', nrow(infiles))
	if(verbose) cat('\nProcessing...')
	rtp.todi2<- infiles[, {				
				load(F)				
				rtp		<- subset(rplkl, 	GROUP==linked.group & 
								TYPE==linked.no &
								((POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE))<conf.cut,
						c(ID1, ID2))
				ans		<- merge(rtp, subset(rplkl, GROUP==linked.group & TYPE==linked.yes), by=c('ID1','ID2'), all.x=1)								
				ans[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE)]						
				ans				
			}, by=c('PTY_RUN')]		
	rtp.todi2	<- subset(rtp.todi2, POSTERIOR_SCORE>0)
	set(rtp.todi2, NULL, 'GROUP', rtp.todi2[, as.character(GROUP)])	
	if(verbose) cat('\nFound (potentially duplicate) pairs between whom linkage is not excluded phylogenetically, n=', nrow(rtp.todi2))
	
	#
	#	gather phylogenetic relationships for each deep-sequence tree
	if(verbose) cat('\nCollect phylogenetic relationship counts for each pair...')
	rplkl	<- infiles[, {
				load(F)
				rplkl			
			}, by='PTY_RUN']
	#	gather summaries of phylogenetic relationships
	if(verbose) cat('\nCollect basic phyloscanner statistics (distance, adjacency, paths between subgraphs) for each pair...')
	rpw		<- infiles[, {
				load(F)
				dwin			
			}, by='PTY_RUN']
	rpw			<- melt(rpw, variable.name='GROUP', value.name='TYPE', measure.vars=c("TYPE_RAW","TYPE_BASIC","TYPE_PAIR_DI2","TYPE_PAIR_TO","TYPE_PAIR_TODI2x2","TYPE_PAIR_TODI2","TYPE_DIR_TODI2","TYPE_NETWORK_SCORES","TYPE_CHAIN_TODI","TYPE_ADJ_DIR_TODI2","TYPE_ADJ_NETWORK_SCORES"))	
	set(rpw, NULL, 'ID_R_MAX', rpw[, pmax(ID1_R,ID2_R)])
	set(rpw, NULL, 'ID_R_MIN', rpw[, pmin(ID1_R,ID2_R)])
	
	#	make sure IDs are characters
	set(rpw, NULL, 'ID1', as.character(rpw$ID1))
	set(rpw, NULL, 'ID2', as.character(rpw$ID2))
	set(rplkl, NULL, 'ID1', as.character(rplkl$ID1))
	set(rplkl, NULL, 'ID2', as.character(rplkl$ID2))
	set(rtp.todi2, NULL, 'ID1', as.character(rtp.todi2$ID1))
	set(rtp.todi2, NULL, 'ID2', as.character(rtp.todi2$ID2))
	
	#
	#	re-arrange pairs so that ID1<ID2
	if(verbose) cat('\nRe-arrange pairs so that ID1<ID2...')
	tmp			<- subset(rplkl, ID1>ID2)
	setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
	set(tmp, NULL, 'TYPE', tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',TYPE)))])
	rplkl		<- rbind(subset(rplkl, !(ID1>ID2)), tmp)
	tmp			<- subset(rpw, ID1>ID2)
	setnames(tmp, c('ID1','ID2','ID1_L','ID1_R','ID2_L','ID2_R','PATHS_12','PATHS_21'), c('ID2','ID1','ID2_L','ID2_R','ID1_L','ID1_R','PATHS_21','PATHS_12'))
	set(tmp, NULL, 'TYPE', tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',TYPE)))])
	rpw			<- rbind(subset(rpw, !(ID1>ID2)), tmp)
	tmp			<- subset(rtp.todi2, ID1>ID2)
	setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
	set(tmp, NULL, 'TYPE', tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',TYPE)))])
	rtp.todi2	<- rbind(subset(rtp.todi2, !(ID1>ID2)), tmp)
	
	#
	#	select analysis in which each pair has highest neff
	if(verbose) cat('\nIf pairs are in several batches, select batch with most deep-sequence phylogenies...')
	setkey(rplkl, ID1, ID2, PTY_RUN, GROUP, TYPE)
	tmp			<- rplkl[ GROUP==linked.group & TYPE==linked.yes, list(PTY_RUN=PTY_RUN[which.max(NEFF)]), by=c('ID1','ID2')]
	rplkl		<- merge(rplkl, tmp, by=c('ID1','ID2','PTY_RUN'))
	rpw			<- merge(rpw, tmp, by=c('ID1','ID2','PTY_RUN'))
	rtp.todi2	<- merge(rtp.todi2, tmp, by=c('ID1','ID2','PTY_RUN'))
	if(verbose) cat('\nLeft with pairs between whom linkage is not excluded phylogenetically, n=', nrow(rtp.todi2))
	
	#	save phylogenetic relationship info for all pairs 
	#	between whom linkage cannot be ruled out 	
	tmp			<- unique( subset(rtp.todi2, select=c(ID1, ID2, PTY_RUN)) )
	rplkl		<- merge(rplkl, tmp, by=c('ID1','ID2','PTY_RUN'))
	rpw			<- merge(rpw, tmp, by=c('ID1','ID2','PTY_RUN'))
	
	#
	#	add meta-data if provided
	if(!is.null(dmeta))
	{		
		stopifnot( 'ID'%in%colnames(dmeta) )
		set(dmeta, NULL, 'ID', as.character(dmeta$ID))
		if(verbose) cat('\nAdd meta-data...')
		tmp			<- unique(dmeta, by='ID')
		setnames(tmp, colnames(tmp), gsub('ID1_ID','ID1',paste0('ID1_',colnames(tmp))))			
		rplkl		<- merge(rplkl, tmp, by=c('ID1'))
		rpw			<- merge(rpw, tmp, by=c('ID1'))
		rtp.todi2	<- merge(rtp.todi2, tmp, by=c('ID1'))
		setnames(tmp, colnames(tmp), gsub('ID1','ID2',colnames(tmp)))		
		rplkl		<- merge(rplkl, tmp, by=c('ID2'))
		rpw			<- merge(rpw, tmp, by=c('ID2'))	
		rtp.todi2	<- merge(rtp.todi2, tmp, by=c('ID2'))		
	}
	
	if(verbose) cat('\nDone. Found pairs, n=', nrow(rtp.todi2), '. Found relationship counts, n=', nrow(rplkl), '. Found phyloscanner statistics, n=', nrow(rpw), '.')
	#	return
	list(linked.pairs=rtp.todi2, relationship.counts=rplkl, windows=rpw)
}




#' @import Rsamtools
#' @import data.table
#' @export
#' @title Calculate position and length of merged reads
#' @description This function calculates the position and length of the two sequenced segments from a single RNA template, potentially after merging when both segments overlap.    
#' @param bam.file.name full path name to bam file.
#' @return data.table with columns QNAME (template query ID), POS (leftmost position of read), LEN (length of read)
phsc.bam.get.length.and.pos.of.mergedreads<- function(bam.file.name, error.strict=TRUE)
{
	dlen	<- scanBam(bam.file.name, param=ScanBamParam(what=c('qname','qwidth','pos','rname','isize','strand')))[[1]]
	dlen	<- as.data.table(dlen)
	setnames(dlen, colnames(dlen), toupper(colnames(dlen)))
	#	check we have at most two segments per template
	tmp		<- dlen[, list(N_SEGMENTS=length(STRAND)), by='QNAME']
	if(error.strict)
		stopifnot( tmp[, all(N_SEGMENTS<3)] )
	if(!error.strict)
	{
		warning('ignoring QNAMES with more than 2 segments, n=',subset(tmp, N_SEGMENTS>2)[, length(QNAME)])
		tmp		<- subset(tmp, N_SEGMENTS<3)
		dlen 	<- merge(dlen, tmp, by='QNAME')
	}
	dlen[, END:= POS+QWIDTH-1L]
	#	determine if segments overlap
	tmp		<- dlen[, list(OVERLAP= as.numeric(max(POS)<=min(END)), LEN= max(END)-min(POS)+1L ), by='QNAME']
	dlen	<- merge(dlen, tmp, by='QNAME')
	#	set LEN for segments that don t overlap
	tmp		<- dlen[, which(OVERLAP==0)]
	set(dlen, tmp, 'LEN', dlen[tmp, QWIDTH])		
	#	get segments that don t overlap
	tmp		<- subset(dlen, OVERLAP==0)
	set(tmp, NULL, 'QNAME', tmp[, paste0(QNAME,':',STRAND)])
	tmp		<- subset(tmp, select=c(QNAME, POS, LEN))
	#	get segments that overlap
	dlen	<- dlen[, list(POS=min(POS), LEN=LEN[1]), by='QNAME']
	dlen	<- rbind(dlen, tmp)
	dlen
}
