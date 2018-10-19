
#' @export
#' @import data.table igraph sna
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
	rtnn	<- phsc.get.most.likely.transmission.chains(rtn, verbose=0)	
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

#' @export
#' @title Combine output from multiple phyloscanner runs
#' @param in.dir Full directory path to output of phyloscanner runs
#' @param save.file If not NA, the combined output is saved to file.
#' @param postfix.trees Pattern at end of file names to identify short read tree output generated by the phyloscanner tools. By default 'trees.rda'.
#' @param postfix.trmwindowstats Pattern at end of file names to identify transmission statistics output per window, generated by the phyloscanner tools. By default 'trmStatsPerWindow.rda'.
#' @param regex.ind Regular expression that identifies the ID of query individuals
#' @param trmw.min.reads Minimum number of reads per individuals in order to make a transmission assignment that involves this individual
#' @return list of three R objects 'phs', 'dtrees', 'dtrms'. See description.
#' @description This function generates three R objects. 'phs' is a list of all short read trees in ape format.
#' 	  'dtrees' is a data.table that provides info on each short read tree. Columns are 'PTY_RUN' (phyloscanner run id), 'W_FROM' (start of window), 'W_TO' (end of window), 'IDX' (index of tree in phs).
#' 	  'dtrms' is a data.table that provides info on transmission assignments. Columns are 'PTY_RUN' (phyloscanner run id), 'ID1' (identifier of first individual), 'ID2' (identifier of second individual), 'TYPE' (window assignment), 
#'	  'WIN_OF_TYPE' (number of windows with that assignment), 'WIN_TOTAL' (Number of windows where reads from both individuals are present), 'PAIR_ID' (unique ID of each combination of PTY_RUN,ID1, ID2)
phsc.combine.phyloscanner.output<- function(in.dir, save.file=NA, postfix.trees='trees.rda', postfix.trmwindowstats='trmStatsPerWindow.rda', regex.ind="^[0-9]+_[0-9]_[0-9]+", trmw.min.reads=1, trmw.min.tips=1, trmw.close.brl=Inf, trmw.distant.brl=Inf, norm.file.name=NA)
{
	#	read trees
	cat('\nread trees')
	ptyr	<- data.table(FILE=list.files(in.dir, pattern=postfix.trees, full.names=TRUE))	
	phs		<- lapply(seq_len(nrow(ptyr)), function(i)
			{
				load(ptyr[i, FILE])
				phs
			})
	phs		<- do.call('c', phs)
	#	read tree info
	dtrees	<- lapply(seq_len(nrow(ptyr)), function(i)
			{
				load(ptyr[i, FILE])
				dfr
			})
	dtrees	<- do.call('rbind', dtrees)
	dtrees[, IDX:=seq_len(nrow(dtrees))]
	#
	#	read transmission window stats
	#	
	cat('\nread transmission window stats')
	ptyr	<- data.table(FILE=list.files(in.dir, pattern=postfix.trmwindowstats, full.names=TRUE))
	ptyr[, PTY_RUN:=ptyr[, as.integer(gsub('ptyr','',regmatches(FILE, regexpr('ptyr[0-9]+',FILE))))]]
	dwin		<- lapply(seq_len(nrow(ptyr)), function(i)
			{
				load(ptyr[i, FILE])
				tt[, PTY_RUN:= ptyr[i,PTY_RUN]]
				tt
			})
	dwin	<- do.call('rbind', dwin)
	setnames(dwin, 	c('PAT.1','PAT.2','PAT.1_TIPS','PAT.2_TIPS','PAT.1_READS','PAT.2_READS','PATHS.12','PATHS.21'),
					c('ID1','ID2','ID1_L','ID2_L','ID1_R','ID2_R','PATHS_12','PATHS_21'))
	set(dwin, NULL, 'ID1', dwin[, regmatches(ID1, regexpr(regex.ind, ID1))])
	set(dwin, NULL, 'ID2', dwin[, regmatches(ID2, regexpr(regex.ind, ID2))])
	set(dwin, NULL, 'PATRISTIC_DISTANCE', dwin[, as.numeric(PATRISTIC_DISTANCE)])
	#	create W_FROM W_TO from SUFFIX
	set(dwin, NULL, 'W_FROM', dwin[, as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\1', SUFFIX))])
	set(dwin, NULL, 'W_TO', dwin[, as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\2', SUFFIX))])
	set(dwin, NULL, 'SUFFIX', NULL)
	#	normalise if desired
	if(!is.na(norm.file.name))
	{
		tmp			<- load(norm.file.name)
		if(length(tmp)!=1)	stop("Expected one R data.table in file",norm.file.name)
		eval(parse(text=paste("norm.table<- ",tmp,sep='')))
		stopifnot( c('W_FROM','W_TO','MEDIAN_PWD')%in%colnames(norm.table) )	 
		setnames(norm.table, 'MEDIAN_PWD', 'NORM_CONST')
		norm.table[, W_MID:= (W_FROM+W_TO)/2]
		#	standardize to 1 on gag+pol ( prot + first part of RT in total 1300bp )
		#	790 - 3385
		cat('\nstandardise normalising constants to 1 on the gag+pol ( prot + first part of RT in total 1300bp pol ) region')
		tmp		<- subset(norm.table, W_MID>=790L & W_MID<=3385L)
		stopifnot( nrow(tmp)>0 )	# norm.table must contain gag+pol region	
		tmp		<- tmp[, mean(NORM_CONST)]
		stopifnot( is.finite(tmp) )
		set(norm.table, NULL, 'NORM_CONST', norm.table[, NORM_CONST/tmp])
		norm.table	<- subset(norm.table, select=c(W_MID, NORM_CONST))
		
		tmp		<- unique(subset(dwin, select=c(W_FROM, W_TO)))
		setkey(tmp, W_FROM)
		tmp		<- tmp[, {
					list( NORM_CONST=subset(norm.table, W_MID>=W_FROM & W_MID<=W_TO)[, mean(NORM_CONST)]	)
				}, by=c('W_FROM','W_TO')]
		dwin	<- merge(dwin, tmp, by=c('W_FROM','W_TO'))
		set(dwin, NULL, 'PATRISTIC_DISTANCE', dwin[, PATRISTIC_DISTANCE/NORM_CONST])		
	}
	#
	#	build transmission summary stats based on selection criteria
	#	
	cat('\nreduce transmission window stats to windows with at least',trmw.min.reads,'reads and at least',trmw.min.tips,'tips')
	dwin	<- subset(dwin, ID1_R>=trmw.min.reads & ID2_R>=trmw.min.reads & ID1_L>=trmw.min.tips & ID2_L>=trmw.min.tips)
	cat('\ntotal number of windows with trm assignments is',nrow(dwin))		
	#
	#	merge assignments (keeping as much detail as needed for full interpretation)
	#
	dwin[, TYPE_RAW:= TYPE]
	#	chains with no intermediate
	tmp		<- dwin[, which(TYPE=="anc_12" & ADJACENT)]
	cat('\nFound adjacent anc_12, n=', length(tmp),'--> chain with no intermediate')
	set(dwin, tmp, 'TYPE', 'chain_12_nointermediate')
	tmp		<- dwin[, which(TYPE=="multi_anc_12" & ADJACENT)]
	cat('\nFound adjacent multi_anc_12, n=', length(tmp),'--> chain with no intermediate')
	set(dwin, tmp, 'TYPE', 'chain_12_nointermediate')
	#	chains with no intermediate
	tmp		<- dwin[, which(TYPE=="anc_21" & ADJACENT)]
	cat('\nFound adjacent anc_21, n=', length(tmp),'--> chain with no intermediate')
	set(dwin, tmp, 'TYPE', 'chain_21_nointermediate')
	tmp		<- dwin[, which(TYPE=="multi_anc_21" & ADJACENT)]
	cat('\nFound adjacent multi_anc_21, n=', length(tmp),'--> chain with no intermediate')
	set(dwin, tmp, 'TYPE', 'chain_21_nointermediate')	
	#
	#	chains with intermediate
	tmp		<- dwin[, which(TYPE=="anc_12" & !ADJACENT)]
	cat('\nFound non-adjacent anc_12, n=', length(tmp),'--> chain with intermediate')
	set(dwin, tmp, 'TYPE', 'chain_12_withintermediate')
	tmp		<- dwin[, which(TYPE=="multi_anc_12" & !ADJACENT)]
	cat('\nFound non-adjacent multi_anc_12, n=', length(tmp),'--> chain with intermediate')
	set(dwin, tmp, 'TYPE', 'chain_12_withintermediate')
	#	chains with intermediate
	tmp		<- dwin[, which(TYPE=="anc_21" & !ADJACENT)]
	cat('\nFound non-adjacent anc_21, n=', length(tmp),'--> chain with intermediate')
	set(dwin, tmp, 'TYPE', 'chain_21_withintermediate')
	tmp		<- dwin[, which(TYPE=="multi_anc_21" & !ADJACENT)]
	cat('\nFound non-adjacent multi_anc_21, n=', length(tmp),'--> chain with intermediate')
	set(dwin, tmp, 'TYPE', 'chain_21_withintermediate')
	#
	#	intermingled with no intermediate
	tmp		<- dwin[, which(TYPE=="conflict" & ADJACENT)]
	cat('\nFound adjacent conflict, n=', length(tmp),'--> intermingled with no intermediate')
	set(dwin, tmp, 'TYPE', 'intermingled_nointermediate')
	#	intermingled with intermediate
	tmp		<- dwin[, which(TYPE=="conflict" & !ADJACENT)]
	cat('\nFound non-adjacent conflict, n=', length(tmp),'--> intermingled with intermediate')
	set(dwin, tmp, 'TYPE', 'intermingled_withintermediate')
	#
	#	other	
	tmp		<- dwin[, which(ADJACENT & PATHS_12==0 & PATHS_21==0)]
	cat('\nFound adjacent with no paths, n=', length(tmp),'--> other')
	set(dwin, tmp, 'TYPE', 'other_nointermediate')
	tmp		<- dwin[, which(!ADJACENT & PATHS_12==0 & PATHS_21==0)]
	cat('\nFound non-adjacent with no assignment, n=', length(tmp),'--> other')
	set(dwin, tmp, 'TYPE', 'other_withintermediate')
	#	check
	stopifnot( !nrow(subset(dwin, TYPE=='none'))	)
	#
	#	add distance as second dimension
	#
	if(!is.na(trmw.close.brl) & is.finite(trmw.close.brl))
	{
		cat('\nidentifying close pairwise assignments using distance=',trmw.close.brl)
		tmp		<- dwin[, which(PATRISTIC_DISTANCE<trmw.close.brl)]
		cat('\nFound close, n=', length(tmp))
		set(dwin, tmp, 'TYPE', dwin[tmp, paste0(TYPE,'_close')])		
	}
	if(!is.na(trmw.distant.brl) & is.finite(trmw.distant.brl))
	{
		cat('\nidentifying distant pairwise assignments using distance=',trmw.distant.brl)
		tmp		<- dwin[, which(PATRISTIC_DISTANCE>=trmw.distant.brl)]
		cat('\nFound distant, n=', length(tmp))
		set(dwin, tmp, 'TYPE', dwin[tmp, paste0(TYPE,'_distant')])	
	}
	#	
	#	summarise transmission stats
	#	
	dtrms	<- dwin[, list(WIN_OF_TYPE=length(W_FROM), ID1_R=ID1_R[1], ID1_L=ID1_L[1], ID2_R=ID2_R[1], ID2_L=ID2_L[1]), by=c('PTY_RUN','ID1','ID2','TYPE')]
	tmp		<- dtrms[, list(WIN_TOTAL=sum(WIN_OF_TYPE)), by=c('PTY_RUN','ID1','ID2')]
	dtrms	<- merge(dtrms, tmp, by=c('PTY_RUN','ID1','ID2'))
	#	set pair id	
	tmp		<- dtrms[, list(SCORE=sum(WIN_OF_TYPE[grepl('close|chain_12_nointermediate|chain_21_nointermediate|intermingled_nointermediate',TYPE)])), by=c('ID1','ID2','PTY_RUN')]
	dtrms	<- merge(dtrms, tmp, by=c('ID1','ID2','PTY_RUN'))
	#	give every pair an ID
	tmp		<- unique(dtrms, by=c('ID1','ID2'))
	tmp		<- tmp[order(-SCORE),]
	tmp[, PAIR_ID:= seq_len(nrow(tmp))]	
	tmp2	<- unique(dtrms, by=c('ID1','ID2','PTY_RUN'))
	tmp		<- merge(tmp2, subset(tmp, select=c(ID1,ID2,PAIR_ID)), by=c('ID1','ID2'))
	tmp		<- tmp[, list(PAIR_ID= paste(PAIR_ID,'-',seq_along(PAIR_ID),sep=''), PTY_RUN=PTY_RUN), by=c('ID1','ID2')]	
	dtrms	<- merge(dtrms, tmp, by=c('ID1','ID2', 'PTY_RUN'))	
	setkey(dtrms, PAIR_ID)
	#	save to file
	if(!is.na(save.file))
	{
		cat('\nwrite to file', save.file)
		save(phs, dtrees, dtrms, dwin, file=save.file)		
	}
	list(phs=phs, dtrees=dtrees, dtrms=dtrms, dwin=dwin)
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

#' @export    
#' @title Calculate basic pairwise relationships
#' @description This function calculates basic pairwise relationships of two individuals in any window. Several different relationship groups can be calculated, for example just using pairwise distance, or using both pairwise distance and topology to define likely pairs.
#' @param df data.table to which basic pairwise relationships will be added. Must contain columns with name 'ADJACENT','TYPE_RAW','PATRISTIC_DISTANCE','PATHS_12','PATHS_21'. 
#' @param trmw.close.brl  Maximum patristic distance between any two read trees from both individuals in a window to classify the individuals as phylogenetically close.
#' @param trmw.distant.brl   Minimum patristic distance between any two read trees from both individuals in a window to classify the individuals as phylogenetically distant.
#' @return input data.table with new relationship column TYPE_BASIC. 
phsc.get.basic.pairwise.relationships<- function(df, trmw.close.brl, trmw.distant.brl, verbose=TRUE)
{
	df[, TYPE_BASIC:= TYPE_RAW]
	
	stopifnot(c('ADJACENT','CONTIGUOUS','TYPE_BASIC','PATRISTIC_DISTANCE','PATHS_12','PATHS_21')%in%colnames(df))
	#
	#	chains_12 
	#
	tmp		<- df[, which(PATHS_12>0 & PATHS_21==0 & ADJACENT & CONTIGUOUS)]
	if(verbose) cat('\nFound PATHS_12>0 & PATHS_21==0 & CONTIGUOUS, n=', length(tmp),'--> chain_12 with no intermediate')
	set(df, tmp, 'TYPE_BASIC', 'chain_12_nointermediate')
	tmp		<- df[, which(PATHS_12>0 & PATHS_21==0 & ADJACENT & !CONTIGUOUS)]
	if(verbose) cat('\nFound PATHS_12>0 & PATHS_21==0 & !CONTIGUOUS, n=', length(tmp),'--> chain_12 with intermediate')
	set(df, tmp, 'TYPE_BASIC', 'chain_12_withintermediate')	
	#
	#	chains_21
	#
	tmp		<- df[, which(PATHS_12==0 & PATHS_21>0 & ADJACENT & CONTIGUOUS)]
	if(verbose) cat('\nFound PATHS_12==0 & PATHS_21>0 & CONTIGUOUS, n=', length(tmp),'--> chain_21 with no intermediate')
	set(df, tmp, 'TYPE_BASIC', 'chain_21_nointermediate')	
	tmp		<- df[, which(PATHS_12==0 & PATHS_21>0 & ADJACENT & !CONTIGUOUS)]
	if(verbose) cat('\nFound PATHS_12==0 & PATHS_21>0 & !CONTIGUOUS, n=', length(tmp),'--> chain_21 with intermediate')
	set(df, tmp, 'TYPE_BASIC', 'chain_21_withintermediate')
	#
	#	intermingled 
	#
	tmp		<- df[, which(PATHS_12>0 & PATHS_21>0 & ADJACENT & CONTIGUOUS)]
	if(verbose) cat('\nFound PATHS_12==0 & PATHS_21==0 & CONTIGUOUS, n=', length(tmp),'--> intermingled with no intermediate')
	set(df, tmp, 'TYPE_BASIC', 'intermingled_nointermediate')
	tmp		<- df[, which(PATHS_12>0 & PATHS_21>0 & ADJACENT & !CONTIGUOUS)]
	if(verbose) cat('\nFound PATHS_12==0 & PATHS_21==0 & !CONTIGUOUS, n=', length(tmp),'--> intermingled with intermediate')
	set(df, tmp, 'TYPE_BASIC', 'intermingled_withintermediate')
	#
	#	sibling
	#
	tmp		<- df[, which(PATHS_12==0 & PATHS_21==0 & ADJACENT & CONTIGUOUS)]
	if(verbose) cat('\nFound PATHS_12==0 & PATHS_21==0 & ADJACENT & CONTIGUOUS, n=', length(tmp),'--> sibling with no intermediate')
	set(df, tmp, 'TYPE_BASIC', 'sibling_nointermediate')
	tmp		<- df[, which(PATHS_12==0 & PATHS_21==0 & ADJACENT & !CONTIGUOUS)]
	if(verbose) cat('\nFound PATHS_12==0 & PATHS_21==0 & ADJACENT & !CONTIGUOUS, n=', length(tmp),'--> sibling with intermediate')
	set(df, tmp, 'TYPE_BASIC', 'sibling_withintermediate')
	#
	#	other
	#	
	tmp		<- df[, which(PATHS_12==0 & PATHS_21==0 & !ADJACENT )]
	if(verbose) cat('\nFound PATHS_12==0 & PATHS_21==0 & !ADJACENT, n=', length(tmp),'--> other')
	set(df, tmp, 'TYPE_BASIC', 'other')
	tmp		<- df[, which(PATHS_12>0 & PATHS_21==0 & !ADJACENT )]
	if(verbose) cat('\nFound PATHS_12>0 & PATHS_21==0 & !ADJACENT, n=', length(tmp),'--> other')
	set(df, tmp, 'TYPE_BASIC', 'other')	
	tmp		<- df[, which(PATHS_12==0 & PATHS_21>0 & !ADJACENT )]
	if(verbose) cat('\nFound PATHS_12==0 & PATHS_21>0 & !ADJACENT, n=', length(tmp),'--> other')
	set(df, tmp, 'TYPE_BASIC', 'other')	
	tmp		<- df[, which(PATHS_12>0 & PATHS_21>0 & !ADJACENT )]
	if(verbose) cat('\nFound PATHS_12>0 & PATHS_21>0 & !ADJACENT, n=', length(tmp),'--> other')
	set(df, tmp, 'TYPE_BASIC', 'other')	
	
	
	#	check
	stopifnot( !nrow(subset(df, TYPE_BASIC=='none'))	)
	#
	#	add distance as second dimension
	#
	if(!is.na(trmw.close.brl) & is.finite(trmw.close.brl))
	{
		if(verbose) cat('\nidentifying close pairwise assignments using distance=',trmw.close.brl)
		tmp		<- df[, which(PATRISTIC_DISTANCE<trmw.close.brl)]
		if(verbose) cat('\nFound close, n=', length(tmp))
		set(df, tmp, 'TYPE_BASIC', df[tmp, paste0(TYPE_BASIC,'_close')])		
	}
	if(!is.na(trmw.distant.brl) & is.finite(trmw.distant.brl))
	{
		if(verbose) cat('\nidentifying distant pairwise assignments using distance=',trmw.distant.brl)
		tmp		<- df[, which(PATRISTIC_DISTANCE>=trmw.distant.brl)]
		if(verbose) cat('\nFound distant, n=', length(tmp))
		set(df, tmp, 'TYPE_BASIC', df[tmp, paste0(TYPE_BASIC,'_distant')])	
	}
	df
}

#' @export
#' @title Wrapper script to calculate pairwise relationships and likelihoods
phsc.get.pairwise.relationships.likelihoods<- function(dwin, trmw.min.reads, trmw.min.tips, trmw.close.brl, trmw.distant.brl, prior.keff, prior.neff, prior.calibrated.prob, relationship.types, prior.keff.dir=prior.keff, prior.neff.dir=prior.neff, verbose=TRUE)
{
	setnames(dwin, 	c('PAT.1','PAT.2','PAT.1_TIPS','PAT.2_TIPS','PAT.1_READS','PAT.2_READS','PATHS.12','PATHS.21'),
					c('ID1','ID2','ID1_L','ID2_L','ID1_R','ID2_R','PATHS_12','PATHS_21'))
	set(dwin, NULL, 'PATRISTIC_DISTANCE', dwin[, as.numeric(PATRISTIC_DISTANCE)])
	#
	#	selection windows
	#	
	if(verbose) cat('\nReducing transmission window stats to windows with at least',trmw.min.reads,'reads and at least',trmw.min.tips,'tips ...')
	dwin	<- subset(dwin, ID1_R>=trmw.min.reads & ID2_R>=trmw.min.reads & ID1_L>=trmw.min.tips & ID2_L>=trmw.min.tips)
	if(verbose) cat('\nTotal number of windows with trm assignments is',nrow(dwin),'...')		
	#
	#	define basic relationship types
	#
	setnames(dwin, 'TYPE', 'TYPE_RAW')
	if(verbose) cat('\nCalculate basic pairwise relationships for windows n=',nrow(dwin),'...')
	dwin	<- phsc.get.basic.pairwise.relationships(dwin, trmw.close.brl, trmw.distant.brl)
	#
	#	derive other relationship types
	#
	if(verbose) cat('\nCalculate derived pairwise relationships for windows n=',nrow(dwin),'...')
	setnames(dwin, 'TYPE_BASIC', 'TYPE_DIR_TODI7x3')	#for backwards compatibility
	dwin	<- phsc.get.pairwise.relationships(dwin, get.groups=relationship.types, make.pretty.labels=FALSE)
	setnames(dwin, 'TYPE_DIR_TODI7x3', 'TYPE_BASIC')
	#
	#	calculate effective K. 
	#	this is based on windows and contiguous chunks of windows
	#
	#	guess W_FROM W_TO from SUFFIX
	if(verbose) cat('\nCalculate KEFF and NEFF for windows n=',nrow(dwin),'...')
	set(dwin, NULL, 'W_FROM', dwin[, as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\1', SUFFIX))])
	set(dwin, NULL, 'W_TO', dwin[, as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\2', SUFFIX))])
	rplkl	<- phsc.get.pairwise.relationships.keff.and.neff(dwin, relationship.types)
	#
	#	calculate marginal posterior for each pairwise relationship state 
	#	this needs a prior which is calibrated as desired.   
	#	the default calibration is KEFF=2 out of NEFF=3 windows are of type t so that the pair is selected to be of type t with posterior probability=50%
	#
	if(verbose) cat('\nCalculate posterior state probabilities for pairs and relationship groups n=',nrow(rplkl),'...')
	rplkl	<- phsc.get.pairwise.relationships.posterior(rplkl, n.type=prior.keff, n.obs=prior.neff, n.type.dir=prior.keff.dir, n.obs.dir=prior.neff.dir, confidence.cut=prior.calibrated.prob)
	#
	#	make TYBE_BASIC labels nice
	#
	tmp		<- rplkl[, which(GROUP=='TYPE_BASIC')]
	set(rplkl, tmp, 'TYPE', rplkl[tmp, gsub('other_withintermediate_distant','other_distant',gsub('other_withintermediate_close','other_close',gsub('other_withintermediate$','other',gsub('other_nointermediate$','other',gsub('other_nointermediate_distant','other_distant',TYPE)))))])	
	set(rplkl, tmp, 'TYPE', rplkl[tmp, gsub('other_no','other\nno',gsub('([ho])intermediate','\\1 intermediate',gsub('intermediate_','intermediate\n',gsub('intermingled_','intermingled\n',gsub('(chain_[fm][mf])_','\\1\n',gsub('(chain_[12][21])_','\\1\n',TYPE))))))])
	set(rplkl, tmp, 'TYPE', rplkl[tmp, gsub('_',' ',TYPE)])
	set(dwin, NULL, 'TYPE_BASIC', dwin[, gsub('other_withintermediate_distant','other_distant',gsub('other_withintermediate_close','other_close',gsub('other_withintermediate$','other',gsub('other_nointermediate$','other',gsub('other_nointermediate_distant','other_distant',TYPE_BASIC)))))])	
	set(dwin, NULL, 'TYPE_BASIC', dwin[, gsub('other_no','other\nno',gsub('([ho])intermediate','\\1 intermediate',gsub('intermediate_','intermediate\n',gsub('intermingled_','intermingled\n',gsub('(chain_[fm][mf])_','\\1\n',gsub('(chain_[12][21])_','\\1\n',TYPE_BASIC))))))])
	set(dwin, NULL, 'TYPE_BASIC', dwin[, gsub('_',' ',TYPE_BASIC)])
	
	list(dwin=dwin, rplkl=rplkl)
}

#' @title Calculate pairwise relationships
#' @description This function calculates pairwise relationships of two individuals in any window. Several different relationship groups can be calculated, for example just using pairwise distance, or using both pairwise distance and topology to define likely pairs.
#' @export    
#' @param df data.table to which new columns of relationship groups will be added. Must contain a column with name TYPE_DIR_TODI7x3. This column contains fundamental relationship states for every window, from which other relationships are derived. 
#' @param get.groups names of relationship groups  
#' @param make.pretty.labels Logical   
#' @return input data.table with new columns. Each new column defines relationship states for a specific relationship group. 
phsc.get.pairwise.relationships<- function(df, get.groups=c('TYPE_PAIR_DI2','TYPE_PAIR_TO','TYPE_PAIR_TODI2x2','TYPE_PAIR_TODI2','TYPE_DIR_TODI2','TYPE_NETWORK_SCORES','TYPE_CHAIN_TODI','TYPE_ADJ_NETWORK_SCORES','TYPE_ADJ_DIR_TODI2'), make.pretty.labels=TRUE)
{
	#	
	#	group to define likely pair just based on distance
	#
	#if('TYPE_PAIR_DI'%in%get.groups)
	#{
	#	df[, TYPE_PAIR_DI:= 'intermediate\ndistance']
	#	set(df, df[, which(grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_DI', 'close')
	#	set(df, df[, which(grepl('distant',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_DI', 'distant')			
	#}
	#if('TYPE_PAIRSCORE_DI'%in%get.groups)
	#{
	#	df[, TYPE_PAIRSCORE_DI:= NA_character_]
	#	set(df, df[, which(grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_PAIRSCORE_DI', 'close')
	#	set(df, df[, which(grepl('distant',TYPE_DIR_TODI7x3))], 'TYPE_PAIRSCORE_DI', 'distant')			
	#}	
	if('TYPE_PAIR_DI2_NOAMB'%in%get.groups)
	{
		df[, TYPE_PAIR_DI2:= 'not close']
		set(df, df[, which(grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_DI2', 'close')			
	}
	if('TYPE_PAIR_DI2'%in%get.groups)
	{
		df[, TYPE_PAIR_DI2:= 'ambiguous']
		set(df, df[, which(grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_DI2', 'close')
		set(df, df[, which(grepl('distant',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_DI2', 'distant')
	}		
	#	
	#	group to define likely pair just based on topology
	#	
	if('TYPE_PAIR_TO'%in%get.groups)
	{
		df[, TYPE_PAIR_TO:= 'other']
		set(df, df[, which(grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('chain',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_TO', 'ancestral/\nintermingled')
		set(df, df[, which(grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('intermingled',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_TO', 'ancestral/\nintermingled')
	}
	#
	#	group for full 2x2 table of distance and topology
	#
	if('TYPE_PAIR_TODI2x2'%in%get.groups)
	{		
		df[, TYPE_PAIR_TODI2x2:= 'not close other']
		set(df, df[, which(grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_TODI2x2', 'close ancestral/\nintermingled/\nsibling')	
		set(df, df[, which(grepl('nointermediate',TYPE_DIR_TODI7x3) & !grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_TODI2x2', 'not close ancestral/\nintermingled/\nsibling')	
		set(df, df[, which(!grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_TODI2x2', 'close other')
	}
	#
	#	group to determine linkage - only 2 states, linked / unlinked
	#
	if('TYPE_PAIR_TODI2_NOAMB'%in%get.groups)
	{
		df[, TYPE_PAIR_TODI2:= 'unlinked']			 		
		set(df, df[, which(grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_TODI2', 'linked')
	}
	if('TYPE_PAIR_TODI2'%in%get.groups)
	{
		df[, TYPE_PAIR_TODI2:= 'unlinked']	
		set(df, df[, which(grepl('nointermediate',TYPE_DIR_TODI7x3) & !grepl('close',TYPE_DIR_TODI7x3) & !grepl('distant',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_TODI2', 'ambiguous')
		set(df, df[, which(grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_TODI2', 'linked')
	}		
	#
	#	group to determine likely pairs
	#
	#if('TYPE_PAIR_TODI'%in%get.groups)
	#{
	#	df[, TYPE_PAIR_TODI:= 'distant/disconnected']	
	#	set(df, df[, which(!grepl('distant',TYPE_DIR_TODI7x3) & !grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_TODI', 'intermediate distance')
	#	set(df, df[, which(grepl('withintermediate',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_TODI', 'distant/disconnected')
	#	set(df, df[, which(grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_PAIR_TODI', 'likely pair')		
	#}
	#	same as TYPE_PAIR_TODI but do not count ambiguous distance
	#if('TYPE_PAIRSCORE_TODI'%in%get.groups)
	#{
	#	df[, TYPE_PAIRSCORE_TODI:= 'distant/disconnected']
	#	set(df, df[, which(grepl('withintermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_PAIRSCORE_TODI', 'distant/disconnected')
	#	set(df, df[, which(grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_PAIRSCORE_TODI', 'likely pair')	
	#	set(df, df[, which(!grepl('distant',TYPE_DIR_TODI7x3) & !grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_PAIRSCORE_TODI', NA_character_)
	#}	
	#
	#	group to determine direction of transmission in likely pairs
	#
	#if('TYPE_DIR_TODI3'%in%get.groups)
	#{				
	#	df[, TYPE_DIR_TODI3:= NA_character_]	# non-NA: all relationships that are used for likely pair	
	#	set(df, df[, which(grepl('other',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'ambiguous')
	#	set(df, df[, which(grepl('intermingled',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'ambiguous')
	#	set(df, df[, which(grepl('chain_',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'ambiguous')		
	#	set(df, df[, which(grepl('chain_fm',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'fm')
	#	set(df, df[, which(grepl('chain_mf',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'mf')
	#	set(df, df[, which(grepl('chain_12',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', '12')
	#	set(df, df[, which(grepl('chain_21',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', '21')
	#}
	#
	#	group to determine direction of transmission in likely pairs
	#	based on contiguous linkage
	#
	if('TYPE_DIR_TODI2'%in%get.groups)
	{
		df[, TYPE_DIR_TODI2:= NA_character_]	# non-NA: all relationships of a likely pair that have a direction assigned	
		set(df, df[, which(grepl('chain_fm',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI2', 'fm')
		set(df, df[, which(grepl('chain_mf',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI2', 'mf')
		set(df, df[, which(grepl('chain_12',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI2', '12')
		set(df, df[, which(grepl('chain_21',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI2', '21')
	}
	#
	#	group to determine direction of transmission in likely pairs
	#	based on adjacent linkage
	#
	if('TYPE_ADJ_DIR_TODI2'%in%get.groups)
	{
		df[, TYPE_ADJ_DIR_TODI2:= NA_character_]	# non-NA: all relationships of a likely pair that have a direction assigned	
		set(df, df[, which(grepl('chain_fm',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_ADJ_DIR_TODI2', 'fm')
		set(df, df[, which(grepl('chain_mf',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_ADJ_DIR_TODI2', 'mf')
		set(df, df[, which(grepl('chain_12',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_ADJ_DIR_TODI2', '12')
		set(df, df[, which(grepl('chain_21',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_ADJ_DIR_TODI2', '21')
	}
	#
	#	group to determine probabilities for transmission networks
	#	based on contiguous linkage
	#
	if('TYPE_NETWORK_SCORES'%in%get.groups)
	{				
		df[, TYPE_NETWORK_SCORES:= 'not close/disconnected']	# non-NA: all relationships that are used for likely pair	
		set(df, df[, which(grepl('sibling',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_NETWORK_SCORES', 'ambiguous')
		set(df, df[, which(grepl('intermingled',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_NETWORK_SCORES', 'ambiguous')				
		set(df, df[, which(grepl('chain_fm',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_NETWORK_SCORES', 'fm')
		set(df, df[, which(grepl('chain_mf',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_NETWORK_SCORES', 'mf')
		set(df, df[, which(grepl('chain_12',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_NETWORK_SCORES', '12')
		set(df, df[, which(grepl('chain_21',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_NETWORK_SCORES', '21')
	}
	#
	#	group to determine probabilities for transmission networks
	#	based on adjacent linkage
	#
	if('TYPE_ADJ_NETWORK_SCORES'%in%get.groups)
	{				
		df[, TYPE_ADJ_NETWORK_SCORES:= 'not close/disconnected']	# non-NA: all relationships that are used for likely pair	
		set(df, df[, which(grepl('sibling',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_ADJ_NETWORK_SCORES', 'ambiguous')
		set(df, df[, which(grepl('intermingled',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_ADJ_NETWORK_SCORES', 'ambiguous')				
		set(df, df[, which(grepl('chain_fm',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_ADJ_NETWORK_SCORES', 'fm')
		set(df, df[, which(grepl('chain_mf',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_ADJ_NETWORK_SCORES', 'mf')
		set(df, df[, which(grepl('chain_12',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_ADJ_NETWORK_SCORES', '12')
		set(df, df[, which(grepl('chain_21',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_ADJ_NETWORK_SCORES', '21')
	}		
	#
	#	group to transmission networks
	#	
	if('TYPE_CHAIN_TODI_NOAMB'%in%get.groups)
	{
		df[, TYPE_CHAIN_TODI:= 'distant']
		set(df, df[, which(grepl('withintermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_CHAIN_TODI', 'chain')
		set(df, df[, which(grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_CHAIN_TODI', 'chain')		
	}
	if('TYPE_CHAIN_TODI'%in%get.groups)
	{
		df[, TYPE_CHAIN_TODI:= 'distant']
		set(df, df[, which(grepl('withintermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_CHAIN_TODI', 'chain')
		set(df, df[, which(grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_CHAIN_TODI', 'chain')
		set(df, df[, which(grepl('withintermediate',TYPE_DIR_TODI7x3) & !grepl('close',TYPE_DIR_TODI7x3) & !grepl('distant',TYPE_DIR_TODI7x3))], 'TYPE_CHAIN_TODI', 'ambiguous')
		set(df, df[, which(grepl('nointermediate',TYPE_DIR_TODI7x3) & !grepl('close',TYPE_DIR_TODI7x3) & !grepl('distant',TYPE_DIR_TODI7x3))], 'TYPE_CHAIN_TODI', 'ambiguous')		
	}
	#
	#	TYPE_DIR_TODI7x3 make pretty labels
	#
	if(make.pretty.labels)
	{
		set(df, NULL, 'TYPE_DIR_TODI7x3', df[, gsub('other_close','other',gsub('other_distant','other',TYPE_DIR_TODI7x3))])	
		set(df, NULL, 'TYPE_DIR_TODI7x3', df[, gsub('([ho])intermediate','\\1 intermediate',gsub('intermediate_','intermediate\n',gsub('intermingled_','intermingled\n',gsub('(chain_[fm][mf])_','\\1\n',gsub('(chain_[12][21])_','\\1\n',TYPE_DIR_TODI7x3)))))])		
		#set(df, NULL, 'TYPE_DIR_TODI7x3', df[, gsub('other_withintermediate_distant','other_distant',gsub('other_withintermediate_close','other_close',gsub('other_withintermediate$','other',gsub('other_nointermediate$','other',gsub('other_nointermediate_distant','other_distant',TYPE_DIR_TODI7x3)))))])	
		#set(df, NULL, 'TYPE_DIR_TODI7x3', df[, gsub('other_no','other\nno',gsub('([ho])intermediate','\\1 intermediate',gsub('intermediate_','intermediate\n',gsub('intermingled_','intermingled\n',gsub('(chain_[fm][mf])_','\\1\n',gsub('(chain_[12][21])_','\\1\n',TYPE_DIR_TODI7x3))))))])		
		for(x in c('TYPE_DIR_TODI7x3',get.groups))
			set(df, NULL, x, gsub('_',' ',df[[x]]))
	}	
	df
}

#' @export
#' @import sna igraph
#' @title Reconstruct most likely transmission chains
#' @description This function reconstructs most likely transmission chains 
#' from the scores associated with directed and undirected edges.
#' @param rtn data.table with network scores for all individuals that could form a network. Must contain columns 'ID1','ID2','IDCLU','GROUP','TYPE','POSTERIOR_SCORE','KEFF'.   
#' @return new data.table with added columns LINK_12 LINK_21 (either 1 or 0), and MX_PROB_12 MX_PROB_21 (associated posterior probabilities)  
phsc.get.most.likely.transmission.chains<- function(rtnn, verbose=0, method='Edmonds')
{
	if(method=='greedy')
		return(phsc.get.most.likely.transmission.chains.greedy(rtnn, verbose=verbose))
	if(method=='Edmonds')
		return(phsc.get.most.likely.transmission.chains.RBGLedmonds(rtnn, verbose=verbose))
}
		

#' @import sna igraph RBGL
#' @title Construct maximum probability transmission network
#' @description This function reconstructs a maximum probility transmission 
#' network from the scores associated with directed and undirected edges.
#' The algorithm starts by keeping the edge with highest score.
#' It then removes the competitor in the opposite direction, and any conflicting edges that would result in indegrees larger than one.
#' By construction, all removed edges have lower probability.
#' The algorithm proceeds until all edges have been processed.
#' @param rtn data.table with network scores for all individuals that could form a network. Must contain columns 'ID1','ID2','IDCLU','GROUP','TYPE','POSTERIOR_SCORE','KEFF'.   
#' @return new data.table with added columns LINK_12 LINK_21 (either 1 or 0), and MX_PROB_12 MX_PROB_21 (associated posterior probabilities)  
phsc.get.most.likely.transmission.chains.RBGLedmonds<- function(rtnn, verbose=0)
{
	require(igraph)
	require(RBGL)
	
	stopifnot(c('ID1','ID2','IDCLU','TYPE','POSTERIOR_SCORE','KEFF')%in%colnames(rtnn))		
	stopifnot( !length(setdiff(c('ambiguous','21','12'), rtnn[, unique(TYPE)])) )
	
	rtnn[, ID1_IN_WEIGHT:=0]
	set(rtnn, rtnn[, which(TYPE=='ambiguous')],'ID1_IN_WEIGHT', 0.5)
	set(rtnn, rtnn[, which(TYPE=='21')],'ID1_IN_WEIGHT', 1)
	rtnn[, ID2_IN_WEIGHT:=0]
	set(rtnn, rtnn[, which(TYPE=='ambiguous')],'ID2_IN_WEIGHT', 0.5)
	set(rtnn, rtnn[, which(TYPE=='12')],'ID2_IN_WEIGHT', 1)
	rtm		<- rtnn[, list(	PROB_21= sum(POSTERIOR_SCORE*ID1_IN_WEIGHT), 
					KEFF_21= sum(KEFF*ID1_IN_WEIGHT),
					PROB_12= sum(POSTERIOR_SCORE*ID2_IN_WEIGHT),
					KEFF_12= sum(KEFF*ID2_IN_WEIGHT)), by=c('IDCLU','CLU_SIZE','ID1','ID2')]
	#
	#	handle networks of size 2 - this is easy
	#
	ans		<- subset(rtm, CLU_SIZE==2)
	set(ans, NULL, c('LINK_12','LINK_21'), 0L)
	set(ans, ans[, which(PROB_12>PROB_21)], 'LINK_12', 1L)
	set(ans, ans[, which(PROB_21>PROB_12)], 'LINK_21', 1L)
		
	#
	#	handle networks of size >2 - use Edmonds algorithm 
	#
	rtm		<- subset(rtm, CLU_SIZE>2)
	rtmm	<- lapply(rtm[, unique(IDCLU)], function(x) subset(rtm, IDCLU==x))
	for(i in seq_along(rtmm))
	{
		#	
		#i	<- 13 
		if(verbose)
			cat('\nIDCLU ',i)
		adj	<- rtmm[[i]]
		adj2<- melt(adj, id.vars=c('ID1','ID2'), measure.vars=c('PROB_12','PROB_21'), value.name='weight')
		tmp	<- subset(adj2, variable=='PROB_21')
		setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
		adj2<- rbind(subset(adj2, variable=='PROB_12'), tmp)
		adj2<- subset(adj2, weight>0)
		#	maximisation is over sum of weights, so transform to log prob
		#	since edmonds cannot handle negative values or zeros, add some constant
		set(adj2, NULL, 'weight', adj2[, log(weight) - min(log(weight)) + 1])
		adj2[, variable:=NULL]		
		g	<- graph.data.frame(adj2)	
		g2	<- igraph.to.graphNEL(g)
		g3	<- edmondsOptimumBranching(g2)
		g3	<- as.data.table(t(g3$edgeList))
		setnames(g3, c('from','to'), c('ID1','ID2'))
		g3[, LINK_12:= 1L]
		g3[, LINK_21:= 0L]
		tmp	<- copy(g3)
		setnames(tmp, c('ID1','ID2', 'LINK_12', 'LINK_21'), c('ID2','ID1', 'LINK_21', 'LINK_12'))
		tmp	<- rbind(g3, tmp)
		adj	<- merge(adj, tmp, by=c('ID1','ID2'))
		rtmm[[i]] <- adj		
	}
	rtm		<- do.call('rbind',rtmm)
	ans		<- rbind(ans, rtm)
	set(ans, ans[, which(LINK_12==1)], 'PROB_21', 0)
	set(ans, ans[, which(LINK_21==1)], 'PROB_12', 0)	
	setnames(ans, c('PROB_21','PROB_12','KEFF_21','KEFF_12'), c('MX_PROB_21','MX_PROB_12','MX_KEFF_21','MX_KEFF_12'))	
	ans		<- merge(rtnn, ans, by=c('ID1','ID2','IDCLU','CLU_SIZE'))	
	ans
}

#' @export   
#' @title Count observed relationship states
#' @description This function counts for each pair of individuals their relationship states across all windows (KEFF), and the total number of windows (NEFF). Since windows can be overlapping, the count is a real value. 
#' @param df data.table with basic relationship types for paired individuals across windows. Must contain columns 'ID1','ID2','W_FROM','W_TO','TYPE_BASIC'. 
#' @param get.groups names of relationship groups  
#' @return new data.table with columns ID1 ID2 GROUP TYPE K KEFF N NEFF. 
phsc.get.pairwise.relationships.keff.and.neff<- function(df, get.groups, w.slide=NA)
{
	stopifnot(c('ID1','ID2','W_FROM','W_TO','TYPE_BASIC')%in%colnames(df))
	#
	#	identify chunks of contiguous windows
	#	
	setkey(df, ID1, ID2, W_FROM)
	if(is.na(w.slide))
	{
		w.slide	<- df[, {
					ans	<- NA_integer_
					tmp	<- diff(W_FROM)
					if(length(tmp))
						ans	<- min(tmp)
					list(W_SLIDE=ans)
				}, by=c('ID1','ID2')]
		w.slide	<- subset(w.slide, !is.na(W_SLIDE))
		w.slide	<- ifelse(nrow(w.slide), w.slide[, min(W_SLIDE)], 1L)		
	}
	#	define chunks
	setkey(df, ID1, ID2, W_FROM)
	tmp		<- df[, {
				tmp<- as.integer( c(TRUE,(W_FROM[-length(W_FROM)]+w.slide)!=W_FROM[-1]) )
				list(W_FROM=W_FROM, W_TO=W_TO, CHUNK=cumsum(tmp))
			}, by=c('ID1','ID2')]
	df		<- merge(df,tmp,by=c('ID1','ID2','W_FROM','W_TO'))
	#	define chunk length in terms of non-overlapping windows	& number of windows in chunk
	tmp		<- df[, {
				list(W_FROM=W_FROM, W_TO=W_TO, CHUNK_L=(max(W_TO+1L)-min(W_FROM))/(W_TO[1]+1L-W_FROM[1]), CHUNK_N=length(W_FROM))
			}, by=c('ID1','ID2','CHUNK')]
	df		<- merge(df,tmp,by=c('ID1','ID2','CHUNK','W_FROM','W_TO'))	
	#	for each chunk, count: windows by type and effective length of chunk
	#	then sum chunks
	rplkl	<- df[, list(	K= length(W_FROM), KEFF= length(W_FROM)/CHUNK_N[1] * CHUNK_L[1]), by=c('ID1','ID2','CHUNK','TYPE_BASIC')]	
	rplkl	<- rplkl[, list(STAT=c('K','KEFF'), V=c(sum(K),sum(KEFF))), by=c('ID1','ID2','TYPE_BASIC')]
	#
	#	add relationship types
	#
	setnames(rplkl, 'TYPE_BASIC', 'TYPE_DIR_TODI7x3')	#for backwards compatibility
	rplkl	<- phsc.get.pairwise.relationships(rplkl, get.groups=get.groups, make.pretty.labels=FALSE)
	#for(x in get.groups)
	#	set(rplkl, NULL, x, gsub('_',' ',rplkl[[x]]))
	setnames(rplkl, 'TYPE_DIR_TODI7x3', 'TYPE_BASIC')
	#	melt relationship groups
	rplkl	<- melt(rplkl, measure.vars=c(get.groups,'TYPE_BASIC'), variable.name='GROUP', value.name='TYPE')
	rplkl	<- subset(rplkl, !is.na(TYPE))
	#	sum K and KEFF of same relationship state
	rplkl	<- rplkl[, list(V=sum(V)), by=c('ID1','ID2','GROUP','TYPE','STAT')]
	#	add zero-count relationship states (change to wide table and set NA's to zero's)
	tmp		<- unique(subset(rplkl, select=c(GROUP,TYPE)))
	tmp2	<- unique(subset(rplkl, select=c(ID1,ID2, STAT)))
	tmp[, DUMMY:=1L]
	tmp2[, DUMMY:=1L]
	tmp		<- merge(tmp, tmp2, by='DUMMY',allow.cartesian=TRUE)
	set(tmp, NULL, 'DUMMY', NULL)
	rplkl	<- merge(tmp, rplkl, all.x=1, by=c('ID1','ID2','GROUP','TYPE','STAT'))
	set(rplkl, rplkl[,which(is.na(V))], 'V', 0)	
	#	expand KEFF and K columns now that everything is done
	rplkl	<- dcast.data.table(rplkl, ID1+ID2+GROUP+TYPE~STAT, value.var='V')	
	#	calculate N and NEFF
	tmp		<- rplkl[, list(N= sum(K), NEFF= sum(KEFF)), by=c('ID1','ID2','GROUP')]	
	rplkl	<- merge(rplkl, tmp, by=c('ID1','ID2','GROUP'))
	rplkl
}

#' @title Calculate marginal posterior probability for two individuals being in a particular relationship state
#' @description This function calculates the parameters that specify the marginal posterior probability for two individuals being in a particular relationship state. The marginal posterior is Beta distributed and this function calculates the ALPHA and BETA parameters.
#' @export  
#' @param df Input data.table   
#' @param n.type Calibration parameter for the prior: minimum number of windows of state to select a pair of individuals with confidence of at least at least confidence.cut, if the total number of windows is n.obs
#' @param n.obs Calibration parameter for the prior: total number of windows. 
#' @param confidence.cut Calibration parameter for the prior: confidence cut off.  
#' @return Input data.table with two additional columns POSTERIOR_ALPHA and POSTERIOR_BETA
phsc.get.pairwise.relationships.posterior<- function(df, n.type=2, n.obs=3, n.type.dir=n.type, n.obs.dir=n.obs, confidence.cut=0.5)
{
	stopifnot(c('GROUP')%in%colnames(df))
	tmp		<- phsc.get.pairwise.relationships.numbers()	
	tmp		<- tmp[, {
				#if(!grepl('_DIR',GROUP))
				#	z<- phsc.get.prior.parameter.n0(N_TYPE, keff=n.type, neff=n.obs, confidence.cut=confidence.cut)
				#if(grepl('_DIR',GROUP))
				#	z<- phsc.get.prior.parameter.n0(N_TYPE, keff=n.type.dir, neff=n.obs.dir, confidence.cut=confidence.cut)
				#list(PAR_PRIOR=z)
				list(PAR_PRIOR=N_TYPE)
			}, by=c('GROUP','N_TYPE')]
	df		<- merge(df, tmp, by=c('GROUP'))
	df[, POSTERIOR_ALPHA:= PAR_PRIOR/N_TYPE+KEFF]
	df[, POSTERIOR_BETA:= PAR_PRIOR*(1-1/N_TYPE)+NEFF-KEFF]	
	df
}