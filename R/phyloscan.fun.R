
#' @export
#' @title Find phylogenetic transmission networks
#' @param rtp Pairs of individuals between whom linkage is not excluded, stored as data.table.
#' @param rplkl Summary of phylogenetic relationship counts for each pair, stored as data.table.
#' @param conf.cut Threshold on the proportion of deep-sequence phylogenies with distant/disconnected subgraphs above which pairs are considered phylogenetically unlinked. Default: 0.6
#' @param neff.cut Threshold on the minimum number of deep-sequence phylogenies with sufficient reads from two individuals to make any phylogenetic inferences. Default: 3.
#' @param verbose Flag to switch on/off verbose mode. Default: TRUE. 
#' @return list of two R objects 'transmission.networks', 'most.likely.transmission.chains'. See description.
#' @description This function reconstructs phylogenetic transmission networks from pairs of individuals between whom linkage is not excluded.
#' @import data.table, igraph, sna 
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
	if(verbose) cat('\nFound transmission networks, n=',rtn[, length(unique(IDCLU))], '. Number of links (either direction and ambiguous)=', nrow(subset(rtn, GROUP==scores.group & TYPE=='not close/disconnected')), '. Number of individuals=', length(unique(c(rtn$ID1, rtnn$ID2))),'.')
	
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
#' @title Find pairs of individuals between whom linkage is not excluded phylogenetically
#' @param indir Full directory path to output of phyloscanner runs
#' @param batch.regex Regular expression that identifies the batch ID of multiple phyloscanner analyses. Default: '^ptyr([0-9]+)_.*'.
#' @param conf.cut Threshold on the proportion of deep-sequence phylogenies with distant/disconnected subgraphs above which pairs are considered phylogenetically unlinked. Default: 0.6
#' @param neff.cut Threshold on the minimum number of deep-sequence phylogenies with sufficient reads from two individuals to make any phylogenetic inferences. Default: 3.
#' @param verbose Flag to switch on/off verbose mode. Default: TRUE. 
#' @param dmeta Optional individual-level meta-data that is to be added to output. Can be NULL.
#' @return list of three R objects 'linked.pairs', 'relationship.counts', 'windows'. See description.
#' @description This function identifies pairs of individuals between whom linkage is not excluded phylogenetically in a large number of phyloscanner analyses, and provides detailed information on them.
#' @import data.table 
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

#' @export
#' @title Generate bash command to resume a single phyloscanner run
#' @param prefix.infiles File name that points phyloscanner output.
#' @param pty.args List of phyloscanner input variables. See examples.
#' @return Character string of phyloscanner commands.
#' @description This function generates bash commands to resume a single phyloscanner run, from the point where all read alignments and read phylogenies were created. The bash script can be called via 'system' in R, or written to file to run on a UNIX system.
#' @example example/ex.cmd.phyloscanner.one.resume.R      
phsc.cmd.phyloscanner.one.resume<- function(prefix.infiles, pty.args)
{	
	#	create local tmp dir
	cmd			<- paste("CWD=$(pwd)\n",sep='\n')
	cmd			<- paste(cmd,"echo $CWD\n",sep='')
	tmpdir		<- paste('pty','_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')	
	tmpdir		<- paste("$CWD/",tmpdir,sep='')
	cmd			<- paste(cmd,'mkdir -p "',tmpdir,'"\n',sep='')
	#	copy required files to local tmp dir	
	file.patient<- list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*patients.txt',sep=''), full.names=TRUE)
	stopifnot(length(file.patient)==1)	
	cmd		<- paste(cmd,'cp "',file.patient,'" "',tmpdir,'"\n',sep='')	
	tmp		<- list.files(dirname(prefix.infiles), pattern=paste('^',basename(prefix.infiles),'.*fasta.zip',sep=''), full.names=TRUE)
	stopifnot(length(tmp)==1)
	cmd		<- paste(cmd,'unzip "',tmp,'" -d "',tmpdir,'"\n',sep='')
	tmp		<- list.files(dirname(prefix.infiles), pattern=paste('^',basename(prefix.infiles),'.*newick.zip',sep=''), full.names=TRUE)
	stopifnot(length(tmp)==1)
	cmd		<- paste(cmd,'unzip "',tmp,'" -d "',tmpdir,'"\n',sep='')
	#	cd to tmp dir
	cmd		<- paste(cmd, 'cd "',tmpdir,'"\n', sep='')	
	#	add all toolkit commands according to pty.args
	cmd		<- paste(cmd, phsc.cmd.process.phyloscanner.output.in.directory(tmpdir, file.patient, pty.args), collapse='\n',sep='')
	#	move all files starting with current run ID
	run.id	<- gsub('_patients.txt','',basename(file.patient))
	cmd		<- paste(cmd, '\nmv ',run.id,'* "',pty.args$out.dir,'"\n',sep='')	
	#	zip up everything else	
	tmp	<- paste(run.id,'_otherstuff.zip',sep='')
	cmd		<- paste(cmd, 'for file in *; do\n\tzip -ur9XTj ',tmp,' "$file"\ndone\n',sep='')
	cmd		<- paste(cmd, 'mv ',tmp,' "',pty.args$out.dir,'"\n',sep='')
	#	clean up
	cmd		<- paste(cmd,'cd $CWD\nrm -r "',tmpdir,'"\n',sep='')		
	cmd
}

#' @export
#' @title Generate bash command to resume a single window of a phyloscanner run
#' @param prefix.infiles File name that points phyloscanner output.
#' @param pty.args List of phyloscanner input variables. See examples.
#' @return Character string of phyloscanner commands.
#' @description This function generates bash commands to resume a single phyloscanner run, from the point where all read alignments and read phylogenies were created. The bash script can be called via 'system' in R, or written to file to run on a UNIX system.
#' @example example/ex.cmd.phyloscanner.one.resume.R      
phsc.cmd.phyloscanner.one.resume.onewindow<- function(prefix.infiles, pty.args)
{	
	stopifnot(!is.na(pty.args$process.window))
	#	create local tmp dir
	cmd			<- paste("CWD=$(pwd)\n",sep='\n')
	cmd			<- paste(cmd,"echo $CWD\n",sep='')
	tmpdir		<- paste('pty','_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')	
	tmpdir		<- paste("$CWD/",tmpdir,sep='')
	cmd			<- paste(cmd,'mkdir -p "',tmpdir,'"\n',sep='')
	#	copy required files to local tmp dir	
	file.patient<- list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*patients.txt',sep=''), full.names=TRUE)
	stopifnot(length(file.patient)==1)	
	cmd		<- paste(cmd,'cp "',file.patient,'" "',tmpdir,'"\n',sep='')	
	tmp		<- list.files(dirname(prefix.infiles), pattern=paste('^',basename(prefix.infiles),'.*fasta.zip',sep=''), full.names=TRUE)
	stopifnot(length(tmp)==1)	
	cmd	<- paste(cmd,'unzip -j "',tmp,'" "*',pty.args$process.window,'*" -d "',tmpdir,'"\n',sep='')	
	tmp		<- list.files(dirname(prefix.infiles), pattern=paste('^',basename(prefix.infiles),'.*newick.zip',sep=''), full.names=TRUE)
	stopifnot(length(tmp)==1)	
	cmd	<- paste(cmd,'unzip -j "',tmp,'" "*',pty.args$process.window,'*" -d "',tmpdir,'"\n',sep='')
	#	cd to tmp dir
	cmd		<- paste(cmd, 'cd "',tmpdir,'"\n', sep='')	
	#	add all toolkit commands according to pty.args
	cmd		<- paste(cmd, phsc.cmd.process.phyloscanner.output.in.directory(tmpdir, file.patient, pty.args), collapse='\n',sep='')
	#	zip up window output
	run.id	<- gsub('_patients.txt','',basename(file.patient))
	tmp		<- paste(run.id,'_output_Window',pty.args$process.window,'.zip',sep='')
	cmd		<- paste(cmd, '\nfor file in ',run.id,'*; do\n\tzip -ur9XTj ',tmp,' "$file"\ndone\n',sep='')
	cmd		<- paste(cmd, 'mv ',tmp,' "',pty.args$out.dir,'"\n',sep='')
	#	zip up everything else
	tmp		<- paste(run.id,'_otherstuff_Window',pty.args$process.window,'.zip',sep='')
	cmd		<- paste(cmd, 'for file in *; do\n\tzip -ur9XTj ',tmp,' "$file"\ndone\n',sep='')
	cmd		<- paste(cmd, 'mv ',tmp,' "',pty.args$out.dir,'"\n',sep='')
	#	clean up
	cmd		<- paste(cmd,'cd $CWD\nrm -r "',tmpdir,'"\n',sep='')		
	cmd
}

#' @export
#' @title Generate bash command to combine single window resumes of a phyloscanner run
#' @param prefix.infiles File name that points phyloscanner output.
#' @param pty.args List of phyloscanner input variables. See examples.
#' @return Character string of phyloscanner commands.
#' @description This function generates bash commands to resume a single phyloscanner run, from the point where all read alignments and read phylogenies were created. The bash script can be called via 'system' in R, or written to file to run on a UNIX system.
#' @example example/ex.cmd.phyloscanner.one.resume.R      
phsc.cmd.phyloscanner.one.resume.combinewindows<- function(prefix.infiles, pty.args)
{	
	stopifnot(pty.args$combine.processed.windows==1)
	#	create local tmp dir
	cmd			<- paste("CWD=$(pwd)\n",sep='\n')
	cmd			<- paste(cmd,"echo $CWD\n",sep='')
	tmpdir		<- paste('pty','_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')	
	tmpdir		<- paste("$CWD/",tmpdir,sep='')
	cmd			<- paste(cmd,'mkdir -p "',tmpdir,'"\n',sep='')
	#	copy required files to local tmp dir	
	file.patient<- list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*patients.txt',sep=''), full.names=TRUE)
	stopifnot(length(file.patient)==1)	
	cmd		<- paste(cmd,'cp "',file.patient,'" "',tmpdir,'"\n',sep='')	
	tmp		<- list.files(dirname(prefix.infiles), pattern=paste('^',basename(prefix.infiles),'output_Window.*',sep=''), full.names=TRUE)
	for(i in seq_along(tmp))
		cmd	<- paste(cmd,'unzip -n "',tmp[i],'" -d "',tmpdir,'"\n',sep='')
	tmp		<- list.files(dirname(prefix.infiles), pattern=paste('^',basename(prefix.infiles),'otherstuff_Window.*',sep=''), full.names=TRUE)
	for(i in seq_along(tmp))
		cmd	<- paste(cmd,'unzip -n "',tmp[i],'" -d "',tmpdir,'"\n',sep='')	
	#	cd to tmp dir
	cmd		<- paste(cmd, 'cd "',tmpdir,'"\n', sep='')	
	#	add all toolkit commands according to pty.args
	cmd		<- paste(cmd, phsc.cmd.process.phyloscanner.output.in.directory(tmpdir, file.patient, pty.args), collapse='\n',sep='')
	#	move all files starting with current run ID
	run.id	<- gsub('_patients.txt','',basename(file.patient))
	cmd		<- paste(cmd, '\nmv ',run.id,'* "',pty.args$out.dir,'"\n',sep='')	
	#	zip up everything else	
	tmp	<- paste(run.id,'_otherstuff.zip',sep='')
	cmd		<- paste(cmd, 'for file in *; do\n\tzip -ur9XTj ',tmp,' "$file"\ndone\n',sep='')
	cmd		<- paste(cmd, 'mv ',tmp,' "',pty.args$out.dir,'"\n',sep='')
	#	clean up
	cmd		<- paste(cmd,'cd $CWD\nrm -r "',tmpdir,'"\n',sep='')		
	cmd	
}

#' @export
#' @title Generate bash command for a single phyloscanner run
#' @param pty.args List of phyloscanner input variables. See examples.
#' @param file.input File name of the file that contains the list of bam files, reference files, and potentially aliases
#' @param file.patient File name of the file that contains the list of unique individuals/units to which the bam files correspond. Multiple bam files are allowed per individual/unit. 
#' @return Character string of phyloscanner commands.
#' @description This function generates bash commands for a single phyloscanner run, that can be called via 'system' in R, or written to file to run on a UNIX system.
#' @example example/ex.cmd.phyloscanner.one.R    
phsc.cmd.phyloscanner.one<- function(pty.args, file.input, file.patient)
{	
	stopifnot(is.character(file.input),is.character(file.patient))
	#	copy input variables into namespace	 
	attach(pty.args)	
	#	sense checks
	stopifnot(is.logical(all.bootstrap.trees))	
	#	define window coordinates
	window.coord			<- integer(0)
	if(!nchar(pty.args[['window.automatic']]))
	{
		stopifnot(length(pty.args[['win']])==4)
		window.coord		<- seq(pty.args[['win']][1], pty.args[['win']][2]-max(pty.args[['win']][3:4]),pty.args[['win']][3])
		window.coord		<- as.vector(rbind( window.coord,window.coord-1+pty.args[['win']][4] ))											
	}	
	#	
	if(!nchar(window.automatic))	stopifnot( is.numeric(window.coord), !length(window.coord)%%2)
	if(nchar(window.automatic))		stopifnot( !length(window.coord) )
	merge.paired.reads			<- ifelse(!is.na(merge.paired.reads) & merge.paired.reads, '--merge-paired-reads', NA_character_)
	keep.overhangs				<- ifelse(!is.na(keep.overhangs) & keep.overhangs, '--keep-overhangs', NA_character_)
	no.trees					<- ifelse(!is.na(no.trees) & no.trees, '--no-trees', NA_character_)
	num.bootstraps				<- ifelse(is.na(no.trees) & !is.na(num.bootstraps) & is.numeric(num.bootstraps) & num.bootstraps>1, as.integer(num.bootstraps), NA_integer_)
	dont.check.duplicates		<- ifelse(!is.na(dont.check.duplicates) & dont.check.duplicates, '--dont-check-duplicates', NA_character_)
	dont.check.recombination	<- ifelse(!is.na(dont.check.recombination) & dont.check.recombination==FALSE, '--check-recombination', NA_character_)	
	#	create local tmp dir
	cmd		<- paste("CWD=$(pwd)\n",sep='\n')
	cmd		<- paste(cmd,"echo $CWD\n",sep='')
	tmpdir	<- paste('pty','_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')	
	tmpdir	<- paste("$CWD/",tmpdir,sep='')
	cmd		<- paste(cmd,'mkdir -p "',tmpdir,'"\n',sep='')
	#	copy files to local tmp dir
	cmd		<- paste(cmd,'cp "',file.input,'" "',tmpdir,'"\n',sep='')	
	cmd		<- paste(cmd,'cp "',file.patient,'" "',tmpdir,'"\n',sep='')
	#	cd to tmp dir
	cmd		<- paste(cmd, 'cd "',tmpdir,'"\n', sep='')	
	cmd		<- paste(cmd, prog.pty,' "',basename(file.input),'" ',sep='')	
	cmd		<- paste(cmd, '--merging-threshold-a', merge.threshold,'--min-read-count',min.read.count,'--quality-trim-ends', quality.trim.ends, '--min-internal-quality',min.internal.quality,'--keep-output-together')
	if(!is.na(merge.paired.reads))
		cmd	<- paste(cmd,' ',merge.paired.reads,sep='')	
	if(!is.na(dont.check.duplicates))
		cmd	<- paste(cmd,' ',dont.check.duplicates,sep='')
	if(!is.na(dont.check.recombination))
		cmd	<- paste(cmd,' ',dont.check.recombination,sep='')	
	if(!is.na(alignments.pairwise.to))		
		cmd	<- paste(cmd,' --pairwise-align-to ',alignments.pairwise.to,sep='')
	if(nchar(window.automatic))
		cmd	<- paste(cmd,' --auto-window-params ', window.automatic,sep='')
	if(!nchar(window.automatic))
		cmd	<- paste(cmd,' --windows ', paste(as.character(window.coord), collapse=','),sep='')
	if(!is.na(alignments.file))
		cmd	<- paste(cmd,' --alignment-of-other-refs "',alignments.file,'"',sep='')	
	if(!is.na(no.trees))		
		cmd	<- paste(cmd, no.trees)
	if(!is.na(num.bootstraps))
		cmd	<- paste(cmd, ' --num-bootstraps ', num.bootstraps, sep='')
	if(!is.na(keep.overhangs))
		cmd	<- paste(cmd, keep.overhangs)
	cmd		<- paste(cmd, '--x-mafft',prog.mafft)
	if(is.na(no.trees))
		cmd	<- paste(cmd, '--x-raxml',prog.raxml)	
	cmd		<- paste(cmd, '\n')
	run.id	<- gsub('_input.csv','',basename(file.input))
	#	process RAxML files
	if(is.na(no.trees) & (is.na(num.bootstraps) | (!is.na(num.bootstraps) & all.bootstrap.trees)))
		cmd	<- paste(cmd, 'for file in RAxML_bestTree*.tree; do\n\tmv "$file" "${file//RAxML_bestTree\\./',run.id,'_}"\ndone\n',sep='')
	if(is.na(no.trees) & !is.na(num.bootstraps) & !all.bootstrap.trees)
		cmd	<- paste(cmd, 'for file in RAxML_bipartitions.MLtreeWbootstraps*.tree; do\n\tmv "$file" "${file//RAxML_bipartitions.MLtreeWbootstraps/',run.id,'_}"\ndone\n',sep='')	
	#cmd	<- paste(cmd, "for file in AlignedReads*.fasta; do\n\tsed 's/<unknown description>//' \"$file\" > \"$file\".sed\n\tmv \"$file\".sed \"$file\"\ndone\n",sep='')		
	if(!is.na(alignments.file) & !is.na(keep.overhangs))
	{
		cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\tcat "$file" | awk \'{if (substr($0,1,4) == ">REF") censor=1; else if (substr($0,1,1) == ">") censor=0; if (censor==0) print $0}\' > NoRef$file\ndone\n', sep='')		
		cmd	<- paste(cmd, 'for file in NoRefAlignedReads*.fasta; do\n\t',phsc.cmd.mafft.add(alignments.file,'"$file"','Ref"$file"', options='--keeplength --memsave --parttree --retree 1'),'\ndone\n',sep='')		
		cmd	<- paste(cmd, 'for file in RefNoRefAlignedReads*.fasta; do\n\t','mv "$file" "${file//RefNoRefAlignedReads/',run.id,'_}"\ndone\n',sep='')		
	}
	if(is.na(alignments.file) || is.na(keep.overhangs))
	{
		cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\tmv "$file" "${file//AlignedReads/',run.id,'_}"\ndone\n',sep='')	
	}	
	#	move Duplicate Read Counts - only for backward compatibility
	#cmd	<- paste(cmd, 'for file in DuplicateReadCountsProcessed_*.csv; do\n\tmv "$file" "${file//DuplicateReadCountsProcessed_/',run.id,'_DuplicateReadCounts_}"\ndone',sep='')
	#	run phyloscanner tools and compress output
	if(is.na(no.trees))
		cmd	<- paste(cmd, phsc.cmd.process.phyloscanner.output.in.directory(tmpdir, file.patient, pty.args), sep='\n')
	#
	out.dir2<- out.dir	
	if(!is.na(no.trees))
	{
		out.dir2<- file.path(out.dir,paste0(run.id,'_trees'))
		cmd		<- paste(cmd,'\nmkdir -p ',out.dir2)		
	}
	cmd		<- paste(cmd, '\nmv ',run.id,'* "',out.dir2,'"\n',sep='')	
	#	zip up everything else
	tmp		<- ''
	if(length(window.coord)==2)
		tmp	<- window.coord[1]
	if(is.null(mem.save) || is.na(mem.save) || mem.save==0)
	{
		cmd		<- paste(cmd, 'for file in *; do\n\tzip -ur9XTjq ',paste(run.id,'_otherstuff',tmp,'.zip',sep=''),' "$file"\ndone\n',sep='')
		cmd		<- paste(cmd, 'mv ',paste(run.id,'_otherstuff',tmp,'.zip',sep=''),' "',out.dir2,'"\n',sep='')		
	}
	#	clean up
	cmd		<- paste(cmd,'cd $CWD\nrm -r "',tmpdir,'"\n',sep='')
	detach(pty.args)
	cmd
}

#' @export
#' @title Generate bash commands for a multiple phyloscanner runs
#' @param pty.runs Data.table of individual assignments to phyloscanner runs, with columns 'PTY_RUN' (run id), 'SAMPLE_ID' (ID of individuals that are assigned to that run). Optional columns: 'RENAME_ID' (new ID for each bam file in phyloscanner output).
#' @param pty.args List of phyloscanner input variables. See examples.
#' @return Data.table with columns 'PTY_RUN' (run id) and 'CMD' (bash commands for that run). 
#' @description This function generates bash commands for multiple phyloscanner runs, that can be called via 'system' in R, or written to file to run on a UNIX system.
#' @example example/ex.cmd.phyloscanner.multi.R  
phsc.cmd.phyloscanner.multi <- function(pty.runs, pty.args, regex.ref='_ref.fasta$', postfix.sample.id='\\.bam|_ref\\.fasta') 		
{
	#
	#	associate BAM and REF files with each scheduled phylotype run
	#	
	set(pty.runs, NULL, 'SAMPLE_ID',  pty.runs[, gsub('\\.bam$','',SAMPLE_ID)])
	#	get available Bam files
	ptyd		<- data.table(FILE=list.files(pty.args[['data.dir']], full.names=TRUE))
	ptyd[, TYPE:=NA_character_]
	set(ptyd, ptyd[, which(grepl('.bam$',FILE))], 'TYPE', 'BAM')
	#	get Reference files
	#	if reference files that were used in assembly are not specified in pty.runs, search for SAMPLE_ID+'_ref.fasta'
	if(!any(colnames(pty.runs)=='REFERENCE_ID'))
	{		
		set(ptyd, ptyd[, which(grepl(regex.ref,FILE))], 'TYPE', 'REF')		
		ptyd		<- subset(ptyd, !is.na(TYPE))
		ptyd[, SAMPLE_ID:= gsub(postfix.sample.id,'',basename(FILE))]
		ptyd		<- dcast.data.table(ptyd, SAMPLE_ID~TYPE, value.var='FILE')		
	}
	#	if reference files that were used in assembly are specified in pty.runs, use these
	if(any(colnames(pty.runs)=='REFERENCE_ID'))
	{
		tmp			<- subset(ptyd, is.na(TYPE))
		tmp[, REFERENCE_ID:= gsub('\\.bam','',basename(FILE))]
		tmp			<- merge(subset(pty.runs, select=c(SAMPLE_ID, REFERENCE_ID)), subset(tmp, select=c(REFERENCE_ID, FILE)), by='REFERENCE_ID')
		tmp[, TYPE:='REF']
		tmp[, REFERENCE_ID:=NULL]
		ptyd		<- subset(ptyd, !is.na(TYPE))
		ptyd[, SAMPLE_ID:= gsub('\\.bam$','',basename(FILE))]
		ptyd		<- rbind(ptyd, tmp)
		ptyd		<- dcast.data.table(ptyd, SAMPLE_ID~TYPE, value.var='FILE')
	}
	#	merge
	ptyd	<- merge(pty.runs, ptyd, by='SAMPLE_ID', all.x=1)	
	if(!any(is.na(pty.args[['select']])))
		ptyd<- subset(ptyd, PTY_RUN%in%pty.args[['select']])			
	#	if background alignment is specified in pty.runs, use it
	if(any(colnames(ptyd)=='BACKGROUND_ID'))
	{
		tmp	<- unique(data.table(BACKGROUND_ID=ptyd[, BACKGROUND_ID], FILE=file.path(pty.args[['data.dir']], ptyd[, BACKGROUND_ID])))
		#	check if files exist
		tmp	<- tmp[, list(EXISTS=file.exists(FILE)), by=c('BACKGROUND_ID','FILE')]
		if(any(tmp[, !EXISTS]))
			warning('\nCould not find location of BACKGROUND files for all runs in pty.runs, n=', tmp[, length(which(!EXISTS))],'\nRuns with missing background alignments are ignored. Please check.')
		ptyd <- merge(ptyd, subset(tmp, EXISTS, c(BACKGROUND_ID, FILE)), by='BACKGROUND_ID')
		set(ptyd, NULL, 'BACKGROUND_ID', NULL)
		setnames(ptyd, 'FILE', 'BACKGROUND_ID')
	}
	if(ptyd[,any(is.na(BAM))])
		warning('\nCould not find location of BAM files for all individuals in pty.runs, n=', ptyd[, length(which(is.na(BAM)))],'\nMissing individuals are ignored. Please check.')	
	if(ptyd[,any(is.na(REF))])
		warning('\nCould not find location of reference files for all individuals in pty.runs, n=', ptyd[, length(which(is.na(REF)))],'\nMissing individuals are ignored. Please check.')	
	setkey(ptyd, PTY_RUN)
	#
	#	write pty.run files and get pty command lines
	#
	pty.c		<- ptyd[, {
				#	PTY_RUN<- z <- 1; BAM<- subset(ptyd, PTY_RUN==z)[, BAM]; REF<- subset(ptyd, PTY_RUN==z)[, REF]
				#	SAMPLE_ID<- subset(ptyd, PTY_RUN==z)[, SAMPLE_ID]; RENAME_ID<- subset(ptyd, PTY_RUN==z)[, RENAME_ID]; BACKGROUND_ID<- subset(ptyd, PTY_RUN==z)[, BACKGROUND_ID]
				file.input		<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_input.csv',sep=''))
				tmp				<- cbind(BAM[!is.na(BAM)&!is.na(REF)], REF[!is.na(BAM)&!is.na(REF)])
				if(exists('RENAME_ID'))
					tmp			<- cbind(tmp, RENAME_ID[!is.na(BAM)&!is.na(REF)])
				write.table(tmp, file=file.input, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')
				file.patient	<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_patients.txt',sep=''))
				tmp				<- paste0(SAMPLE_ID[!is.na(BAM)&!is.na(REF)],'.bam')
				if(exists('RENAME_ID'))
					tmp			<- unique(RENAME_ID[!is.na(BAM)&!is.na(REF)])
				if(exists('UNIT_ID'))
					tmp			<- unique(UNIT_ID[!is.na(BAM)&!is.na(REF)])				
				cat( paste(tmp,collapse='\n'), file= file.patient	)
				if(exists('BACKGROUND_ID'))
				{
					if(length(unique(BACKGROUND_ID))!=1)
						stop('\nrun',PTY_RUN,' Expected one background alignment file, found: ',paste(unique(BACKGROUND_ID), collapse=''),'. Please check.')
					pty.args$alignments.file<- unique(BACKGROUND_ID)
				}				
				cmd			<- phsc.cmd.phyloscanner.one(	pty.args, 
															file.input, 
															file.patient)
				#cmd			<- paste(cmd, pty.cmd.evaluate.fasta(pty.args[['out.dir']], strip.max.len=pty.args[['strip.max.len']], select=paste('^ptyr',PTY_RUN,'_In',sep=''), min.ureads.individual=pty.args[['min.ureads.individual']]), sep='')
				#cat(cmd)
				list(CMD= cmd)				
			},by='PTY_RUN']
	pty.c
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
#' @import data.table grid ggtree
#' @title Plot all short read phylogenies including two individuals
#' @description This function plots short read phylogenies and highlights the clades of two individuals in red and blue.  
#' @param phs List of trees in ape format
#' @param dfs data.table with mandatory column 'IDX' and optional column 'TITLE'. IDX is the index of all phylogenies in 'phs' that are to be plotted. TITLE is a title for each sub-plot, for example specifying a window.
#' @param id1	Regular expression that identifies the first individual, to be plotted in red.  
#' @param id2	Regular expression that identifies the first individual, to be plotted in blue.
#' @param pdf.h	Height of the pdf file in inches.
#' @param pdf.rw Relative width of the pdf file, internally multiplied by the number of phylogenies to give the total width in inches.
#' @param pdf.ntrees Number of trees per pdf.
#' @param pdf.title.size Size of pdf title in inches.
#' @param plot.file If not missing, the phylogenies will be printed to file.	
#' @return List of ggtree objects, ready for printing.
phsc.plot.phy.selected.pairs<- function(phs, dfs, id1, id2, plot.file=NA, pdf.h=50, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40)
{
	#	determine which phylogenies contain both individuals
	tmp		<- copy(dfs)
	tmp		<- merge(tmp, tmp[, {
				ph	<- phs[[ IDX ]]
				list(HAS_BOTH_IND= any(grepl(id1, attr(ph, "INDIVIDUAL"))) & any(grepl(id2, attr(ph, "INDIVIDUAL"))))  				
			}, by='IDX'], by='IDX')
	tmp		<- subset(tmp, HAS_BOTH_IND)
	phps	<- lapply(seq_len(nrow(tmp)), function(i){
				ph.title	<- NULL
				if('TITLE'%in%colnames(tmp))
					ph.title	<- tmp[i, TITLE]										
				ph			<- phs[[ tmp[i, IDX] ]]
				col			<- rep('grey50', length(attr(ph, "INDIVIDUAL"))) 
				col[ grepl(id1, attr(ph, "INDIVIDUAL")) ]	<- 'red'
				col[ grepl(id2, attr(ph, "INDIVIDUAL")) ]	<- 'blue'
				attr(ph, 'COLOUR')	<- col								
				#						
				p 			<- ggtree(ph, aes(color=I(COLOUR))) +
						geom_point2(shape = 16, size=3, aes(subset=NODE_SHAPES)) +
						scale_fill_hue(na.value="black") +								
						theme(legend.position="none") +
						geom_tiplab(aes(col=I(COLOUR))) +
						theme_tree2() +
						theme(legend.position="bottom", plot.title = element_text(size=pdf.title.size)) + 
						ggplot2::xlim(0, max(node.depth.edgelength(ph)[1:Ntip(ph)])*1.3) +
						labs(x='subst/site', title=ph.title)						
				p
			})
	#
	#	single page plot
	#		
	if(!is.na(plot.file))					
	{
		if(length(phps)<=pdf.ntrees)
		{
			cat('Plotting to file', plot.file,'...\n')
			pdf(file=plot.file, w=pdf.rw*length(phps), h=pdf.h)
			grid.newpage()
			pushViewport(viewport(layout=grid.layout(1, length(phps))))
			for(i in seq_along(phps))
				print(phps[[i]], vp = viewport(layout.pos.row=1, layout.pos.col=i))
			dev.off()
		}
		if(length(phps)>pdf.ntrees)
		{
			pi	<- data.table(IDX=seq_along(phps))
			pi[, PLOT:= ceiling(IDX/pdf.ntrees)]
			pi[, PLOT_IDX:= (IDX-1)%%pdf.ntrees+1]
			pi[,{
						cat('Plotting to file', gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file),'...\n')
						pdf(file=gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file), w=pdf.rw*pdf.ntrees, h=pdf.h)
						grid.newpage()
						pushViewport(viewport(layout=grid.layout(1, pdf.ntrees)))
						for(i in seq_along(IDX))
							print(phps[[IDX[i]]], vp = viewport(layout.pos.row=1, layout.pos.col=PLOT_IDX[i]))
						dev.off()
					}, by='PLOT']
		}
	}
	phps	
}


phsc.plot.windowscan.for.pairs<- function(rpw2, plot.file, plot.w=10, plot.h=10, id.cols=c('ID1','ID2'), ylim=NULL, cols.typet=NULL)
{		
	#	make manual plot to show intermingled
	if(is.null(ylim))
		ylim	<- c(1e-3,0.4)
	if(is.null(cols.typet))
		cols.typet			<- c(	"ancestral 1->2"=brewer.pal(11, 'PiYG')[2],
				'ancestral m->f'='steelblue2',
				"ancestral 2->1"=brewer.pal(11, 'PiYG')[10],
				'ancestral f->m'='hotpink2',
				"intermingled"=brewer.pal(11, 'PuOr')[4], 
				'sibling'=rev(brewer.pal(11, 'PuOr'))[c(3)], 
				"disconnected"=rev(brewer.pal(11, 'RdGy'))[4])		
	group		<- 'TYPE_BASIC'
	stopifnot(group%in%unique(rpw2$GROUP))
	stopifnot(id.cols%in%colnames(rpw2))
	#
	rpw3		<- subset(rpw2, GROUP==group)
	rpw3[, TYPE_TO:= 'disconnected']
	rpw3[, Y:=1e-3]
	set(rpw3, rpw3[, which(PATRISTIC_DISTANCE<1e-3)],'PATRISTIC_DISTANCE',1.1e-3)
	set(rpw3, rpw3[,which(grepl('chain 12', TYPE))], 'TYPE_TO', 'ancestral 1->2')
	set(rpw3, rpw3[,which(grepl('chain 21', TYPE))], 'TYPE_TO', 'ancestral 2->1')
	set(rpw3, rpw3[,which(grepl('chain mf', TYPE))], 'TYPE_TO', 'ancestral m->f')
	set(rpw3, rpw3[,which(grepl('chain fm', TYPE))], 'TYPE_TO', 'ancestral f->m')	
	set(rpw3, rpw3[,which(grepl('intermingled', TYPE))], 'TYPE_TO', 'intermingled')
	set(rpw3, rpw3[,which(grepl('sibling', TYPE))], 'TYPE_TO', 'sibling')	
	set(rpw3, NULL, 'TYPE_TO', rpw3[, factor(TYPE_TO, levels=c('ancestral 1->2','ancestral m->f','ancestral 2->1','ancestral f->m','intermingled','sibling','disconnected'))])	
	setnames(rpw3, id.cols, c('PRIVATECOL_ID1','PRIVATECOL_ID2'))	
	#	
	ggplot(rpw3, aes(x=W_FROM)) +
			geom_hline(yintercept=0.025, colour='grey50') +
			geom_bar(aes(y=Y, fill=TYPE_TO), colour='transparent', stat='identity', width=25) +
			geom_point(aes(y=PATRISTIC_DISTANCE), size=1) +				
			labs(x='\ngenomic position\n(relative to HXB2)', y='subgraph distance\n(subst/site)\n',fill='topological subgraph\nrelationship') +
			scale_x_continuous(breaks=seq(0,1e4,500), minor_breaks=seq(0,1e4,100), limits=c(rpw3[, min(W_FROM)]-diff(rpw3[1:2, W_FROM]), rpw3[, max(W_FROM)]+diff(rpw3[1:2, W_FROM]))) +
			scale_y_log10(labels=percent, expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			coord_cartesian(ylim=ylim) +
			scale_fill_manual(values=cols.typet) +
			theme_bw() + 
			theme(legend.position='bottom', panel.spacing = unit(1, "lines")) +
			facet_grid(PRIVATECOL_ID1+PRIVATECOL_ID2~.)	
	ggsave(file=plot.file, w=plot.w, h=plot.h, useDingbats=FALSE, limitsize=FALSE)
}	

#' @export
#' @import data.table grid ggtree colorspace
#' @title Plot short read phylogenies and highlight individuals
#' @description This function plots short read phylogenies and highlights the clades of two individuals in red and blue.  
#' @param phs List of trees in ape format
#' @param dfs data.table with mandatory column 'IDX' and optional column 'TITLE'. IDX is the index of all phylogenies in 'phs' that are to be plotted. TITLE is a title for each sub-plot, for example specifying a window.
#' @param ids Vector of regular expressions that identify individuals to be highlighted in colour.
#' @param plot.cols Vector of colours for each individual
#' @param group.redo Logical, indicating if the colour groups should be recalculated from ids.
#' @param drop.blacklisted Logical, indicating if all blacklisted taxa should be dropped prior to plotting.  
#' @param pdf.h	Height of the pdf file in inches.
#' @param pdf.rw Relative width of the pdf file, internally multiplied by the number of phylogenies to give the total width in inches.
#' @param pdf.ntrees Number of trees per pdf.
#' @param pdf.title.size Size of pdf title in inches.
#' @param plot.file If not missing, the phylogenies will be printed to file.	
#' @return List of ggtree objects, ready for printing.
phsc.plot.phycollapsed.selected.individuals<- function(phs, dfs, ids, plot.cols=rainbow(length(ids)), plot.cols.background=function(n){ rainbow_hcl(n, start = 60, end = 240, c=30, l=80) }, drop.blacklisted=TRUE, tip.regex='(.*)_read_([0-9]+)_count_([0-9]+)', plot.file=NA, drop.less.than.n.ids=2, pdf.h=50, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40)
{	
	#	determine which phylogenies contain at least one of the requested individuals
	tmp		<- copy(dfs)
	tmp		<- merge(tmp, tmp[, {
						ph	<- phs[[ IDX ]]
						z	<- ph$tip.label
						if(drop.blacklisted)
							z	<- as.character(attr(ph, "INDIVIDUAL"))
						list(HAS_N_IND=length(which(sapply(ids, function(id) any(grepl(id, z)) ))))						  				
					}, by='IDX'], by='IDX')
	tmp		<- subset(tmp, HAS_N_IND>=drop.less.than.n.ids)
	#	define colours so that they don t change across windows
	cols	<- tmp[, {
				ph	<- phs[[ IDX ]]
				phb			<- data.table(	TAXA=ph$tip.label,
											IDX=seq_along(ph$tip.label), 											 
											ID= gsub(tip.regex,'\\1',ph$tip.label),
											COUNT=sub(tip.regex,'\\3',ph$tip.label)		)
				set(phb, phb[, which(ID==COUNT)], 'COUNT', NA_character_)
				set(phb, NULL, 'COUNT', phb[, as.integer(COUNT)])
				set(phb, phb[, which(is.na(COUNT))], 'ID','REFERENCE')	
				list(ID= phb[, unique(ID)])				
			}, by='IDX']
	cols	<- unique(subset(cols, select=ID))
	cols[,COL:= plot.cols.background(nrow(cols))]	
	set(cols, cols[,which(grepl('REFERENCE',ID))], 'COL', 'grey50')
	for(i in seq_along(ids))
		set(cols, cols[, which(grepl(ids[i], ID))], 'COL', plot.cols[i])
	#	
	phps	<- lapply(seq_len(nrow(tmp)), function(i){
				#cat(i,'\n')
				ph.title	<- NULL
				if('TITLE'%in%colnames(tmp))
					ph.title	<- tmp[i, TITLE]										
				ph			<- phs[[ tmp[i, IDX] ]]
				ph			<- phsc.phy.collapse.monophyletic.clades(ph, drop.blacklisted=TRUE, tip.regex=tip.regex)
				#
				#	define IDs, find COUNTS, find references
				phb			<- data.table(	TAXA=ph$tip.label,
											IDX=seq_along(ph$tip.label), 											 
											ID= gsub(tip.regex,'\\1',ph$tip.label),
											COUNT=sub(tip.regex,'\\3',ph$tip.label)		)
				set(phb, phb[, which(ID==COUNT)], 'COUNT', NA_character_)
				set(phb, NULL, 'COUNT', phb[, as.integer(COUNT)])
				set(phb, phb[, which(is.na(COUNT))], 'ID','REFERENCE')	
				set(phb, phb[, which(is.na(COUNT))], 'COUNT', 1L)
				#	define cols for this tree
				tmp					<- data.table(IDX= seq_along(attr(ph, "INDIVIDUAL")), ID= attr(ph, "INDIVIDUAL"))
				set(tmp, tmp[, which(is.na(ID))], 'ID', 'REFERENCE')
				tmp					<- merge(tmp, cols, by='ID')
				attr(ph, 'COLOUR')	<- tmp[order(IDX),][,COL]	
				#						
				p 			<- ggtree(ph, aes(color=I(COLOUR))) %<+% phb +
						geom_tippoint(aes(size=COUNT), shape=18) +
						#geom_point2(shape=16, aes(size=NODE_S, subset=NODE_S)) +
						scale_fill_hue(na.value="black") +								
						theme(legend.position="none") +
						geom_tiplab(align=T, aes(col=I(COLOUR)), size=1, linetype=NA, linesize=NA) +
						guides(size='none') +
						theme_tree2() +
						theme(	axis.line.x=element_line(),
								panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3),
								legend.position="bottom",								
								plot.title = element_text(hjust = 0.5, size=pdf.title.size)) + 
						ggplot2::xlim(0, max(node.depth.edgelength(ph)[1:Ntip(ph)])*1.15) +
						labs(x='subst/site', title=ph.title)
				#ggsave(plot.file, w=10, h=200, limitsize = FALSE)
				p
			})
	#
	#	single page plot
	#		
	if(!is.na(plot.file))					
	{
		if(length(phps)<=pdf.ntrees)
		{
			cat('Plotting to file', plot.file,'...\n')
			pdf(file=plot.file, w=pdf.rw*length(phps), h=pdf.h)
			grid.newpage()
			pushViewport(viewport(layout=grid.layout(1, length(phps))))
			for(i in seq_along(phps))
				print(phps[[i]], vp = viewport(layout.pos.row=1, layout.pos.col=i))
			dev.off()
		}
		if(length(phps)>pdf.ntrees)
		{
			pi	<- data.table(IDX=seq_along(phps))
			pi[, PLOT:= ceiling(IDX/pdf.ntrees)]
			pi[, PLOT_IDX:= (IDX-1)%%pdf.ntrees+1]
			pi[,{
						cat('Plotting to file', gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file),'...\n')
						pdf(file=gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file), w=pdf.rw*pdf.ntrees, h=pdf.h)
						grid.newpage()
						pushViewport(viewport(layout=grid.layout(1, pdf.ntrees)))
						for(i in seq_along(IDX))
							print(phps[[IDX[i]]], vp = viewport(layout.pos.row=1, layout.pos.col=PLOT_IDX[i]))
						dev.off()
					}, by='PLOT']
		}
	}
	phps	
}

#' @export
#' @import data.table grid ggtree ggnet
#' @title Plot maximum probability network
#' @description This function plots a maximum probability network.  
#' @param df data.table with the following columns  "IDCLU","ID1", "ID2", "TYPE","KEFF","LKL_MAX","POSTERIOR_SCORE" 
#' @param di data.table with meta-data to customize the plot with columns  "ID", node.shape, node.label, node.fill 
#' @param point.size size of the individual points
#' @param point.sizec.couple size of the outer ring around individuals in couples
#' @param edge.gap value to adjust start / end points of edges
#' @param edge.size multiplier by which the size of edges is shrunk/magnified
#' @param curvature curvature of directed edges  
#' @param arrow type of arrow to be plotted
#' @param curv.shift offset to place the label for directed edges
#' @param label.size size of label
#' @param node.label Text displayed on top of each node 
#' @param node.shape column name in di by which the shape of each node is drawn 
#' @param node.fill column name in di by which each node is coloured
#' @param node.shape.values named vector associating shapes to the values in the node.shape column
#' @param node.fill.values named vector associating colours to the values in the node.fill column
#' @param threshold.linked treshold value between 0 and 1. Edges with weight above this treshold are shown in black.
#' @return ggplot object
phsc.plot.max.probability.network<- function(df, di, point.size=10, edge.gap=0.04, edge.size=0.4, curvature= -0.2, arrow=arrow(length=unit(0.04, "npc"), type="open"), curv.shift=0.08, label.size=3, node.label='ID', node.shape=NA_character_, node.fill=NA_character_, node.shape.values=NA_integer_, node.fill.values=NA_character_, threshold.linked=NA_real_, layout=NULL)
{
	#point.size=10; point.size.couple=14; edge.gap=0.04; edge.size=0.4; curvature= -0.2; arrow=arrow(length=unit(0.04, "npc"), type="open"); curv.shift=0.08; label.size=3
	#node.label='ID'; node.shape='IN_COUPLE'; node.fill='SEX'
	#node.fill.values=c('F'='hotpink2', 'M'='steelblue2')
	#node.shape.values=c('not in long-term\nrelationship'=18,'in long-term\nrelationship'=16)	
	if(is.na(node.label))
	{
		node.label<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.label, NA_character_)
	}
	if(is.na(node.shape))
	{
		node.shape<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.shape, 'NA')
	}
	if(is.na(node.fill))
	{
		node.fill<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.fill, 'NA')
	}
	if(any(is.na(node.fill.values)))
	{
		z						<- unique(di[[node.fill]])
		node.fill.values		<- heat.colors(length(z))
		names(node.fill.values)	<- z
	}
	if(any(is.na(node.shape.values)))
	{
		z						<- unique(di[[node.shape]])
		node.shape.values		<- seq_along(z)
		names(node.shape.values)<- z
	}
	setnames(di, c(node.label, node.shape, node.fill), c('NODE_LABEL','NODE_SHAPE','NODE_FILL'))
	tmp	<- c('NODE_LABEL','NODE_SHAPE','NODE_FILL')[which(c(node.label, node.shape, node.fill)=='ID')]
	if(length(tmp))
		set(di, NULL, 'ID', di[[tmp]])
	di	<- subset(di, select=c(ID, NODE_LABEL, NODE_SHAPE, NODE_FILL))
	if(is.null(layout))
	{
		layout	<- as.data.table(ggnet2(network(unique(subset(df, select=c(ID1,ID2))), directed=FALSE, matrix.type="edgelist"))$data[,c("label", "x", "y")])		
	}		
	if(any(grepl('label', colnames(layout))))
	{		
		setnames(layout, c('label','x','y'), c('ID1','ID1_X','ID1_Y'))
	}
	if(any(colnames(layout)=='ID'))
	{
		setnames(layout, c('ID','X','Y'), c('ID1','ID1_X','ID1_Y'))		
	}	
	df		<- merge(df, layout, by='ID1')
	setnames(layout, c('ID1','ID1_X','ID1_Y'), c('ID2','ID2_X','ID2_Y'))
	df		<- merge(df, layout, by='ID2')
	setnames(layout, c('ID2','ID2_X','ID2_Y'),  c('ID','X','Y'))	
	layout	<- merge(layout,di, by='ID')		
	df[, EDGETEXT_X:= (ID1_X+ID2_X)/2]
	df[, EDGETEXT_Y:= (ID1_Y+ID2_Y)/2]
	df[, EDGE_LABEL:= paste0('D',round(100*pmax(NETWORK_SCORE_12,NETWORK_SCORE_21),d=1),'%',' // ','L',round(100*POSTERIOR_SCORE_LINKED,d=1),'%' ) ]
	df[, EDGE_COL:= as.character(factor(POSTERIOR_SCORE_LINKED>threshold.linked, levels=c(TRUE,FALSE),labels=c('edge_col_2','edge_col_1')))]	
	#	for edges, move the start and end points on the line between X and Y
	#	define unit gradient
	df[, MX:= (ID2_X - ID1_X)]	
	df[, MY:= (ID2_Y - ID1_Y)]
	tmp		<- df[, sqrt(MX*MX+MY*MY)]
	set(df, NULL, 'MX', df[, MX/tmp])
	set(df, NULL, 'MY', df[, MY/tmp])	
	set(df, NULL, 'ID1_X', df[, ID1_X + MX*edge.gap])
	set(df, NULL, 'ID1_Y', df[, ID1_Y + MY*edge.gap])
	set(df, NULL, 'ID2_X', df[, ID2_X - MX*edge.gap])
	set(df, NULL, 'ID2_Y', df[, ID2_Y - MY*edge.gap])	
	p		<- ggplot() +			
			geom_point(data=layout, aes(x=X, y=Y, colour=NODE_FILL, pch=NODE_SHAPE), size=point.size) +
			geom_segment(data=subset(df, LINK_12==1), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*MX_KEFF_12, colour=EDGE_COL), arrow=arrow, lineend="butt") +			
			geom_segment(data=subset(df, LINK_21==1), aes(x=ID2_X, xend=ID1_X, y=ID2_Y, yend=ID1_Y, size=edge.size*MX_KEFF_21, colour=EDGE_COL), arrow=arrow, lineend="butt") +						
			scale_colour_manual(values=c(node.fill.values, 'edge_col_1'='grey80', 'edge_col_2'='grey40','NA'='grey50')) +
			scale_shape_manual(values=c(node.shape.values, 'NA'=16)) +
			scale_fill_manual(values=c(node.fill.values, 'NA'='grey50')) +
			scale_size_identity() +
			geom_text(data=df, aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL), size=label.size) +
			geom_text(data=layout, aes(x=X, y=Y, label=NODE_LABEL)) +
			theme_void() +
			guides(colour='none', fill='none',size='none', pch='none')
	layout		<- subset(layout, select=c(ID,X,Y))
	setnames(layout, c('ID','X','Y'), c('label','x','y'))
	p$layout	<- layout
	p	
}

#' @export
#' @import data.table grid ggtree ggnet
#' @title Plot probability network
#' @description This function plots a probability network.  
#' @param df data.table with the following columns  "IDCLU","ID1", "ID2", "TYPE","KEFF","LKL_MAX","POSTERIOR_SCORE" 
#' @param di data.table with meta-data to customize the plot with columns  "ID", node.shape, node.label, node.fill 
#' @param point.size size of the individual points
#' @param point.sizec.couple size of the outer ring around individuals in couples
#' @param edge.gap value to adjust start / end points of edges
#' @param edge.size multiplier by which the size of edges is shrunk/magnified
#' @param curvature curvature of directed edges  
#' @param arrow type of arrow to be plotted
#' @param curv.shift offset to place the label for directed edges
#' @param label.size size of label
#' @param node.label Text displayed on top of each node 
#' @param node.shape column name in di by which the shape of each node is drawn 
#' @param node.fill column name in di by which each node is coloured
#' @param node.shape.values named vector associating shapes to the values in the node.shape column
#' @param node.fill.values named vector associating colours to the values in the node.fill column
#' @param threshold.linked treshold value between 0 and 1. Edges with weight above this treshold are shown in black.
#' @return ggplot object
phsc.plot.probability.network<- function(df, di, point.size=10, point.size.couple=point.size*1.4, edge.gap=0.04, edge.size=0.4, curvature= -0.2, arrow=arrow(length=unit(0.04, "npc"), type="open"), curv.shift=0.08, label.size=3, node.label='ID', node.shape=NA_character_, node.fill=NA_character_, node.shape.values=NA_integer_, node.fill.values=NA_character_, threshold.linked=NA_real_)
{	
	#point.size=10; point.size.couple=14; edge.gap=0.04; edge.size=0.4; curvature= -0.2; arrow=arrow(length=unit(0.04, "npc"), type="open"); curv.shift=0.08; label.size=3
	#node.label='ID'; threshold.linked=0.6; node.shape=NA_character_; node.fill=NA_character_; node.shape.values=NA_integer_; node.fill.values=NA_character_
	#node.shape='IN_COUPLE'; node.fill='SEX'
	#node.fill.values=c('F'='hotpink2', 'M'='steelblue2')
	#node.shape.values=c('not in long-term\nrelationship'=18,'in long-term\nrelationship'=16)
	if(is.na(node.label))
	{
		node.label<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.label, NA_character_)
	}
	if(is.na(node.shape))
	{
		node.shape<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.shape, 'NA')
	}
	if(is.na(node.fill))
	{
		node.fill<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.fill, 'NA')
	}
	if(any(is.na(node.fill.values)))
	{
		z						<- unique(di[[node.fill]])
		node.fill.values		<- heat.colors(length(z))
		names(node.fill.values)	<- z
	}
	if(any(is.na(node.shape.values)))
	{
		z						<- unique(di[[node.shape]])
		node.shape.values		<- seq_along(z)
		names(node.shape.values)<- z
	}
	setnames(di, c(node.label, node.shape, node.fill), c('NODE_LABEL','NODE_SHAPE','NODE_FILL'))
	tmp	<- c('NODE_LABEL','NODE_SHAPE','NODE_FILL')[which(c(node.label, node.shape, node.fill)=='ID')]
	if(length(tmp))
		set(di, NULL, 'ID', di[[tmp]])
	di	<- subset(di, select=c(ID, NODE_LABEL, NODE_SHAPE, NODE_FILL))
	
	layout	<- as.data.table(ggnet2(network(unique(subset(df, select=c(ID1,ID2))), directed=FALSE, matrix.type="edgelist"))$data[,c("label", "x", "y")])
	setnames(layout, c('label','x','y'), c('ID1','ID1_X','ID1_Y'))
	df		<- merge(df, layout, by='ID1')
	setnames(layout, c('ID1','ID1_X','ID1_Y'), c('ID2','ID2_X','ID2_Y'))
	df		<- merge(df, layout, by='ID2')
	setnames(layout, c('ID2','ID2_X','ID2_Y'),  c('ID','X','Y'))	
	layout	<- merge(layout,di, by='ID')	
	
	df[, EDGETEXT_X:= (ID1_X+ID2_X)/2]
	df[, EDGETEXT_Y:= (ID1_Y+ID2_Y)/2]
	#
	#	calculate score for linked
	if(is.na(threshold.linked))
	{
		df	<- merge(df,df[, 	{
									z<- rep('edge_col_1', length(TYPE))
									z[which.max(POSTERIOR_SCORE)]	<- 'edge_col_2'
									list(EDGE_COL=z, TYPE=TYPE)	
								}, by=c('ID1','ID2')], by=c('ID1','ID2','TYPE'))		
	}
	if(!is.na(threshold.linked))
	{
		tmp	<- subset(df, TYPE!='not close/disconnected')[, list( EDGE_COL=as.character(factor(sum(POSTERIOR_SCORE)>=threshold.linked, levels=c(TRUE, FALSE), labels=c('edge_col_2','edge_col_1'))) ), by=c('ID1','ID2')]
		df	<- merge(df, tmp, by=c('ID1','ID2'))		
	}	
	#	for edges, move the start and end points on the line between X and Y
	#	define unit gradient
	df[, MX:= (ID2_X - ID1_X)]	
	df[, MY:= (ID2_Y - ID1_Y)]
	tmp		<- df[, sqrt(MX*MX+MY*MY)]
	set(df, NULL, 'MX', df[, MX/tmp])
	set(df, NULL, 'MY', df[, MY/tmp])	
	set(df, NULL, 'ID1_X', df[, ID1_X + MX*edge.gap])
	set(df, NULL, 'ID1_Y', df[, ID1_Y + MY*edge.gap])
	set(df, NULL, 'ID2_X', df[, ID2_X - MX*edge.gap])
	set(df, NULL, 'ID2_Y', df[, ID2_Y - MY*edge.gap])	
	#	label could just be move on the tangent vector to the line
	#	define unit tangent
	df[, TX:= -MY]
	df[, TY:= MX]
	tmp		<- df[, which(TYPE=='12')]
	set(df, tmp, 'EDGETEXT_X', df[tmp, EDGETEXT_X + TX*curv.shift])
	set(df, tmp, 'EDGETEXT_Y', df[tmp, EDGETEXT_Y + TY*curv.shift])
	tmp		<- df[, which(TYPE=='21')]
	set(df, tmp, 'EDGETEXT_X', df[tmp, EDGETEXT_X - TX*curv.shift])
	set(df, tmp, 'EDGETEXT_Y', df[tmp, EDGETEXT_Y - TY*curv.shift])
	#
	p		<- ggplot() +			
			geom_point(data=layout, aes(x=X, y=Y, colour=NODE_FILL, pch=NODE_SHAPE), size=point.size) +
			geom_segment(data=subset(df, EDGE_COL=='edge_col_1' & TYPE=='ambiguous' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), lineend="butt") +
			geom_curve(data=subset(df, EDGE_COL=='edge_col_1' & TYPE=='12' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +
			geom_curve(data=subset(df, EDGE_COL=='edge_col_1' & TYPE=='21' & KEFF>0), aes(x=ID2_X, xend=ID1_X, y=ID2_Y, yend=ID1_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +
			geom_segment(data=subset(df, EDGE_COL=='edge_col_2' & TYPE=='ambiguous' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), lineend="butt") +
			geom_curve(data=subset(df, EDGE_COL=='edge_col_2' & TYPE=='12' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +
			geom_curve(data=subset(df, EDGE_COL=='edge_col_2' & TYPE=='21' & KEFF>0), aes(x=ID2_X, xend=ID1_X, y=ID2_Y, yend=ID1_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +									
			scale_colour_manual(values=c(node.fill.values, 'edge_col_1'='grey80', 'edge_col_2'='grey40','NA'='grey50')) +
			scale_shape_manual(values=c(node.shape.values, 'NA'=16)) +
			scale_fill_manual(values=c(node.fill.values, 'NA'='grey50')) +
			scale_size_identity() +
			geom_text(data=subset(df, TYPE!='not close/disconnected' & KEFF>0), aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=paste0(round(100*POSTERIOR_SCORE,d=1),'%')), size=label.size) +
			geom_text(data=layout, aes(x=X, y=Y, label=NODE_LABEL)) +
			theme_void() +
			guides(colour='none', fill='none',size='none', pch='none') 
	layout		<- subset(layout, select=c(ID,X,Y))
	setnames(layout, c('ID','X','Y'), c('label','x','y'))	
	p$layout	<- layout
	p
}


#' @export
#' @import data.table grid ggtree ggnet
#' @title Plot probability network with most likely edges
#' @description This function plots the network showing the most likely edge types.  
#' @param df data.table with the following columns  "IDCLU","ID1", "ID2", "TYPE","KEFF","LKL_MAX","POSTERIOR_SCORE" 
#' @param di data.table with meta-data to customize the plot with columns  "ID", node.shape, node.label, node.fill 
#' @param point.size size of the individual points
#' @param point.sizec.couple size of the outer ring around individuals in couples
#' @param edge.gap value to adjust start / end points of edges
#' @param edge.size multiplier by which the size of edges is shrunk/magnified
#' @param curvature curvature of directed edges  
#' @param arrow type of arrow to be plotted
#' @param curv.shift offset to place the label for directed edges
#' @param label.size size of label
#' @param node.label Text displayed on top of each node 
#' @param node.shape column name in di by which the shape of each node is drawn 
#' @param node.fill column name in di by which each node is coloured
#' @param node.shape.values named vector associating shapes to the values in the node.shape column
#' @param node.fill.values named vector associating colours to the values in the node.fill column
#' @param threshold.linked treshold value between 0 and 1. Edges with weight above this treshold are shown in black.
#' @return ggplot object
phsc.plot.maxedge.network<- function(df, di, point.size=10, point.size.couple=point.size*1.4, edge.gap=0.04, edge.size=0.4, curvature= -0.2, arrow=arrow(length=unit(0.04, "npc"), type="open"), curv.shift=0.08, label.size=3, node.label='ID', node.shape=NA_character_, node.fill=NA_character_, node.shape.values=NA_integer_, node.fill.values=NA_character_, edge.label=NA_character_, edge.label.values=NA_character_, threshold.linked=NA_real_)
{	
	#point.size=10; point.size.couple=14; edge.gap=0.04; edge.size=0.4; curvature= -0.2; arrow=arrow(length=unit(0.04, "npc"), type="open"); curv.shift=0.08; label.size=3
	#node.label='ID'; threshold.linked=0.6; node.shape=NA_character_; node.fill=NA_character_; node.shape.values=NA_integer_; node.fill.values=NA_character_
	#node.shape='IN_COUPLE'; node.fill='SEX'
	#node.fill.values=c('F'='hotpink2', 'M'='steelblue2')
	#node.shape.values=c('not in long-term\nrelationship'=18,'in long-term\nrelationship'=16)
	if(is.na(node.label))
	{
		node.label<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.label, NA_character_)
	}
	if(is.na(edge.label))
	{
		edge.label<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(df, NULL, edge.label, df[, TYPE])
		edge.label.values	<- c('12'='red','21'='red','ambiguous'='grey50')
	}
	if(is.na(node.shape))
	{
		node.shape<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.shape, 'NA')
	}
	if(is.na(node.fill))
	{
		node.fill<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.fill, 'NA')
	}
	if(any(is.na(node.fill.values)))
	{
		z						<- unique(di[[node.fill]])
		node.fill.values		<- heat.colors(length(z))
		names(node.fill.values)	<- z
	}
	if(any(is.na(node.shape.values)))
	{
		z						<- unique(di[[node.shape]])
		node.shape.values		<- seq_along(z)
		names(node.shape.values)<- z
	}
	setnames(di, c(node.label, node.shape, node.fill), c('NODE_LABEL','NODE_SHAPE','NODE_FILL'))
	tmp	<- c('NODE_LABEL','NODE_SHAPE','NODE_FILL')[which(c(node.label, node.shape, node.fill)=='ID')]
	if(length(tmp))
		set(di, NULL, 'ID', di[[tmp]])
	di	<- subset(di, select=c(ID, NODE_LABEL, NODE_SHAPE, NODE_FILL))
	setnames(df, c(edge.label), c('EDGE_LABEL'))
	
	layout	<- as.data.table(ggnet2(network(unique(subset(df, select=c(ID1,ID2))), directed=FALSE, matrix.type="edgelist"))$data[,c("label", "x", "y")])
	setnames(layout, c('label','x','y'), c('ID1','ID1_X','ID1_Y'))
	df		<- merge(df, layout, by='ID1')
	setnames(layout, c('ID1','ID1_X','ID1_Y'), c('ID2','ID2_X','ID2_Y'))
	df		<- merge(df, layout, by='ID2')
	setnames(layout, c('ID2','ID2_X','ID2_Y'),  c('ID','X','Y'))	
	layout	<- merge(layout,di, by='ID')	
	
	df[, EDGETEXT_X:= (ID1_X+ID2_X)/2]
	df[, EDGETEXT_Y:= (ID1_Y+ID2_Y)/2]
	#
	#	calculate score for linked
	if(is.na(threshold.linked))
	{
		df	<- merge(df,df[, 	{
							z<- rep('alpha_1', length(TYPE))
							z[which.max(POSTERIOR_SCORE)]	<- 'alpha_2'
							list(ALPHA=z, TYPE=TYPE)	
						}, by=c('ID1','ID2')], by=c('ID1','ID2','TYPE'))		
	}
	if(!is.na(threshold.linked))
	{
		tmp	<- subset(df, TYPE!='not close/disconnected')[, list( ALPHA=as.character(factor(sum(POSTERIOR_SCORE)>=threshold.linked, levels=c(TRUE, FALSE), labels=c('alpha_2','alpha_1'))) ), by=c('ID1','ID2')]
		df	<- merge(df, tmp, by=c('ID1','ID2'))		
	}	
	#	for edges, move the start and end points on the line between X and Y
	#	define unit gradient
	df[, MX:= (ID2_X - ID1_X)]	
	df[, MY:= (ID2_Y - ID1_Y)]
	tmp		<- df[, sqrt(MX*MX+MY*MY)]
	set(df, NULL, 'MX', df[, MX/tmp])
	set(df, NULL, 'MY', df[, MY/tmp])	
	set(df, NULL, 'ID1_X', df[, ID1_X + MX*edge.gap])
	set(df, NULL, 'ID1_Y', df[, ID1_Y + MY*edge.gap])
	set(df, NULL, 'ID2_X', df[, ID2_X - MX*edge.gap])
	set(df, NULL, 'ID2_Y', df[, ID2_Y - MY*edge.gap])	
	#	label could just be move on the tangent vector to the line
	#	define unit tangent
	df[, TX:= -MY]
	df[, TY:= MX]
	tmp		<- df[, which(TYPE=='12')]
	set(df, tmp, 'EDGETEXT_X', df[tmp, EDGETEXT_X + TX*curv.shift])
	set(df, tmp, 'EDGETEXT_Y', df[tmp, EDGETEXT_Y + TY*curv.shift])
	tmp		<- df[, which(TYPE=='21')]
	set(df, tmp, 'EDGETEXT_X', df[tmp, EDGETEXT_X - TX*curv.shift])
	set(df, tmp, 'EDGETEXT_Y', df[tmp, EDGETEXT_Y - TY*curv.shift])
	
	tmp		<- df[, list(TYPE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
	df		<- merge(df, tmp, by=c('ID1','ID2','TYPE'))
	#
	p		<- ggplot() +			
			geom_point(data=layout, aes(x=X, y=Y, fill=NODE_FILL, pch=NODE_SHAPE), size=point.size) +
			geom_segment(data=subset(df, TYPE=='ambiguous' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, alpha=ALPHA, colour=EDGE_LABEL), lineend="butt") +
			geom_curve(data=subset(df, TYPE=='12' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, alpha=ALPHA, colour=EDGE_LABEL), curvature=curvature, arrow=arrow, lineend="butt") +
			geom_curve(data=subset(df, TYPE=='21' & KEFF>0), aes(x=ID2_X, xend=ID1_X, y=ID2_Y, yend=ID1_Y, size=edge.size*KEFF, alpha=ALPHA, colour=EDGE_LABEL), curvature=curvature, arrow=arrow, lineend="butt") +
			scale_shape_manual(values=c(node.shape.values, 'NA'=16)) +
			scale_fill_manual(values=c(node.fill.values, 'NA'='grey50')) +
			scale_alpha_manual(values=c('alpha_1'=0.5,'alpha_2'=1, 'NA'=0)) +
			scale_colour_manual(values=c(edge.label.values, 'NA'='grey50')) +
			scale_size_identity() +
			geom_text(data=subset(df, TYPE!='not close/disconnected' & KEFF>0), aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=paste0(round(100*POSTERIOR_SCORE,d=1),'%')), size=label.size) +
			geom_text(data=layout, aes(x=X, y=Y, label=NODE_LABEL)) +
			theme_void() +
			guides(colour='none', fill='none',size='none', pch='none') 
	layout		<- subset(layout, select=c(ID,X,Y))
	setnames(layout, c('ID','X','Y'), c('label','x','y'))	
	p$layout	<- layout
	p
}

#' @export
#' @import data.table grid ggtree
#' @title Plot short read phylogenies and highlight individuals
#' @description This function plots short read phylogenies and highlights the clades of two individuals in red and blue.  
#' @param phs List of trees in ape format
#' @param dfs data.table with mandatory column 'IDX' and optional column 'TITLE'. IDX is the index of all phylogenies in 'phs' that are to be plotted. TITLE is a title for each sub-plot, for example specifying a window.
#' @param ids Vector of regular expressions that identify individuals to be highlighted in colour.
#' @param plot.cols Vector of colours for each individual
#' @param group.redo Logical, indicating if the colour groups should be recalculated from ids.
#' @param drop.blacklisted Logical, indicating if all blacklisted taxa should be dropped prior to plotting.  
#' @param pdf.h	Height of the pdf file in inches.
#' @param pdf.rw Relative width of the pdf file, internally multiplied by the number of phylogenies to give the total width in inches.
#' @param pdf.ntrees Number of trees per pdf.
#' @param pdf.title.size Size of pdf title in inches.
#' @param plot.file If not missing, the phylogenies will be printed to file.	
#' @return List of ggtree objects, ready for printing.
phsc.plot.phy.selected.individuals<- function(phs, dfs, ids, plot.cols=rainbow(length(ids)), plot.file=NA, group.redo=FALSE, drop.less.than.n.ids=2, drop.blacklisted=FALSE, pdf.h=50, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40)
{	
	#	determine which phylogenies contain at least one of the requested individuals
	tmp		<- copy(dfs)
	tmp		<- merge(tmp, tmp[, {
						ph	<- phs[[ IDX ]]
						z	<- ph$tip.label
						if(drop.blacklisted)
							z	<- as.character(attr(ph, "INDIVIDUAL"))
						list(HAS_N_IND=length(which(sapply(ids, function(id) any(grepl(id, z)) ))))						  				
					}, by='IDX'], by='IDX')
	tmp		<- subset(tmp, HAS_N_IND>=drop.less.than.n.ids)
	phps	<- lapply(seq_len(nrow(tmp)), function(i){
				ph.title	<- NULL
				if('TITLE'%in%colnames(tmp))
					ph.title	<- tmp[i, TITLE]										
				ph			<- phs[[ tmp[i, IDX] ]]
				#
				#
				if(drop.blacklisted)
				{
					z						<- data.table(FROM=ph$edge[,1],TO=ph$edge[,2])
					z						<- merge(z, data.table(TO= seq_len(length(attr(ph,'INDIVIDUAL'))), GROUP= attr(ph,'INDIVIDUAL')), by='TO')
					z						<- subset(z, is.na(GROUP) & TO<=Ntip(ph))[, TO]
					ph						<- drop.tip(ph, z)
					attr(ph,'NODE_SHAPES')	<- rep(FALSE, Nnode(ph, internal.only=FALSE))
					attr(ph,'INDIVIDUAL')	<- NULL
					attr(ph,'SPLIT')		<- NULL
				}
				#
				#
				if(group.redo||drop.blacklisted)
				{
					phb			<- data.table(IDX=seq_along(ph$tip.label), TAXA=ph$tip.label, ID='none')
					for(id in ids)
						set(phb, phb[, which(grepl(id,TAXA))], 'ID', id)
					tmp				<- lapply( phb[, unique(ID)], function(x)	subset(phb, ID==x)[, TAXA]	)
					names(tmp)		<- phb[, unique(ID)]
					ph				<- groupOTU(ph, tmp, group='INDIVIDUAL')
					#z	<- merge(data.table(FROM=ph$edge[,1],IDX=ph$edge[,2]), phb, by='IDX', all=1)
					#z[, GROUP:= attr(ph,'INDIVIDUAL')[1:nrow(ph$edge)]]
					#z	<- unique(subset(z, !is.na(ID), select=c(ID, GROUP)))
					#attr(ph,'INDIVIDUAL')	<- factor(attr(ph,'INDIVIDUAL'), levels=c(0,z[,as.character(GROUP)]), labels=c('not characterized',z[,FILE_ID]))
					#ph									
				}				
				#
				#
				cols		<- rep('grey50', length(attr(ph, "INDIVIDUAL"))) 
				for(i in seq_along(ids))
					cols[ grepl(ids[i], attr(ph, "INDIVIDUAL")) ]<- plot.cols[i]
				attr(ph, 'COLOUR')	<- cols								
				#						
				p 			<- ggtree(ph, aes(color=I(COLOUR))) +
						#geom_point2(shape = 16, size=3, aes(subset=NODE_SHAPES)) +
						scale_fill_hue(na.value="black") +								
						theme(legend.position="none") +
						geom_tiplab(aes(col=I(COLOUR))) +
						theme_tree2() +
						theme(legend.position="bottom", plot.title = element_text(size=pdf.title.size)) + 
						ggplot2::xlim(0, max(node.depth.edgelength(ph)[1:Ntip(ph)])*1.3) +
						labs(x='subst/site', title=ph.title)
				#ggsave(plot.file, w=10, h=200, limitsize = FALSE)
				p
			})
	#
	#	single page plot
	#		
	if(!is.na(plot.file))					
	{
		if(length(phps)<=pdf.ntrees)
		{
			cat('Plotting to file', plot.file,'...\n')
			pdf(file=plot.file, w=pdf.rw*length(phps), h=pdf.h)
			grid.newpage()
			pushViewport(viewport(layout=grid.layout(1, length(phps))))
			for(i in seq_along(phps))
				print(phps[[i]], vp = viewport(layout.pos.row=1, layout.pos.col=i))
			dev.off()
		}
		if(length(phps)>pdf.ntrees)
		{
			pi	<- data.table(IDX=seq_along(phps))
			pi[, PLOT:= ceiling(IDX/pdf.ntrees)]
			pi[, PLOT_IDX:= (IDX-1)%%pdf.ntrees+1]
			pi[,{
						cat('Plotting to file', gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file),'...\n')
						pdf(file=gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file), w=pdf.rw*pdf.ntrees, h=pdf.h)
						grid.newpage()
						pushViewport(viewport(layout=grid.layout(1, pdf.ntrees)))
						for(i in seq_along(IDX))
							print(phps[[IDX[i]]], vp = viewport(layout.pos.row=1, layout.pos.col=PLOT_IDX[i]))
						dev.off()
					}, by='PLOT']
		}
	}
	phps	
}

#' @title Calculate basic pairwise relationships
#' @description This function calculates basic pairwise relationships of two individuals in any window. Several different relationship groups can be calculated, for example just using pairwise distance, or using both pairwise distance and topology to define likely pairs.
#' @export    
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
#' @title Reconstruct most likely transmission chains
#' @description This function reconstructs most likely transmission chains 
#' from the scores associated with directed and undirected edges.
#' @import sna igraph
#' @param rtn data.table with network scores for all individuals that could form a network. Must contain columns 'ID1','ID2','IDCLU','GROUP','TYPE','POSTERIOR_SCORE','KEFF'.   
#' @return new data.table with added columns LINK_12 LINK_21 (either 1 or 0), and MX_PROB_12 MX_PROB_21 (associated posterior probabilities)  
phsc.get.most.likely.transmission.chains<- function(rtnn, verbose=0, method='Edmonds')
{
	if(method=='greedy')
		return(phsc.get.most.likely.transmission.chains.greedy(rtnn, verbose=verbose))
	if(method=='Edmonds')
		return(phsc.get.most.likely.transmission.chains.RBGLedmonds(rtnn, verbose=verbose))
}
		


#' @title Construct maximum probability transmission network
#' @description This function reconstructs a maximum probility transmission 
#' network from the scores associated with directed and undirected edges.
#' The algorithm starts by keeping the edge with highest score.
#' It then removes the competitor in the opposite direction, and any conflicting edges that would result in indegrees larger than one.
#' By construction, all removed edges have lower probability.
#' The algorithm proceeds until all edges have been processed.
#' @import sna igraph RBGL
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