Tasp.180513.process.phyloscanneroutput<- function()
{
	require(data.table)	
	require(igraph)
	require(sna)
	
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Q20_results/output_q20w150'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Q20_results/output_q20w160'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Q20_results/output_q20w170'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Q20_results/output_q20w180'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Q20_results/output_q20w190'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Q20_results/output_q20w200'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Q23_results/output_q23w150'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Q23_results/output_q23w160'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Q23_results/output_q23w170'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Q23_results/output_q23w180'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Q23_results/output_q23w190'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Q23_results/output_q23w200'
	outfile.base <- gsub('output','phcs_180513',indir)
	outfile	<- paste0(outfile.base,'.rda')
	
	neff.cut	<- 3
	conf.cut	<- 0.6 
	if(0)
	{
		linked.group	<- 'TYPE_PAIR_TODI2'
		linked.type.yes	<- 'linked'
		linked.type.no	<- 'unlinked'
		scores.group	<- 'TYPE_NETWORK_SCORES'
		scores.type.no	<- c('ambiguous','not close/disconnected') 
		dir.group		<- 'TYPE_DIR_TODI2'		
	}
	if(1)
	{
		close.group		<- 'TYPE_PAIR_DI2'
		close.type.yes	<- 'close'
		close.type.no	<- 'distant'
		linked.group	<- 'TYPE_CHAIN_TODI'
		linked.type.yes	<- 'chain'
		linked.type.no	<- 'distant'	
		scores.group	<- 'TYPE_ADJ_NETWORK_SCORES'
		scores.type.no	<- c('ambiguous','not close/disconnected') 
		dir.group		<- 'TYPE_ADJ_DIR_TODI2'
	}
	
	
	
	#
	#	from every phyloscanner run, select pairs that are most frequently linked 
	#	by: distance, distance + topology
	infiles	<- data.table(F=list.files(indir, pattern='pairwise_relationships.rda', full.names=TRUE))
	infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(F)))]
	setkey(infiles, PTY_RUN)
	rtp.todi2<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p10_d50_stagetwo_rerun23_min30_adj_chain_mean/ptyr40_pairwise_relationships.rda'
				#cat(PTY_RUN,'\n')
				load(F)
				#	all pairs that are not decisively unlinked based on "close.group"
				rtp		<- subset(rplkl, 	GROUP==close.group & 
								TYPE==close.type.no &
								((POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE))<conf.cut,
						c(ID1, ID2))
				ans		<- merge(rtp, subset(rplkl, GROUP==close.group & TYPE==close.type.yes), by=c('ID1','ID2'), all.x=1)
				#	all pairs that are not decisively unlinked based on ML likely transmission pairs by distance + topology
				rtp		<- subset(rplkl, 	GROUP==linked.group & 
								TYPE==linked.type.no &
								((POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE))<conf.cut,
						c(ID1, ID2))
				rtp		<- merge(rtp, subset(rplkl, GROUP==linked.group & TYPE==linked.type.yes), by=c('ID1','ID2'), all.x=1)
				ans		<- rbind(ans, rtp)				
				ans[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE)]						
				ans				
			}, by=c('PTY_RUN')]		
	rtp.todi2	<- subset(rtp.todi2, POSTERIOR_SCORE>0)
	set(rtp.todi2, NULL, 'GROUP', rtp.todi2[, as.character(GROUP)])	
	#
	#	prepare all dwin and rplkl and save separately
	rplkl	<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr1_pairwise_relationships.rda'
				#cat(PTY_RUN,'\n')
				load(F)
				rplkl			
			}, by='PTY_RUN']
	rpw		<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr1_pairwise_relationships.rda'
				#cat(PTY_RUN,'\n')
				load(F)
				dwin			
			}, by='PTY_RUN']
	#	melt rpw
	rpw			<- melt(rpw, variable.name='GROUP', value.name='TYPE', measure.vars=c("TYPE_RAW","TYPE_BASIC","TYPE_PAIR_DI2","TYPE_PAIR_TO","TYPE_PAIR_TODI2x2","TYPE_PAIR_TODI2","TYPE_DIR_TODI2","TYPE_NETWORK_SCORES","TYPE_CHAIN_TODI","TYPE_ADJ_DIR_TODI2","TYPE_ADJ_NETWORK_SCORES"))	
	set(rpw, NULL, 'ID_R_MAX', rpw[, pmax(ID1_R,ID2_R)])
	set(rpw, NULL, 'ID_R_MIN', rpw[, pmin(ID1_R,ID2_R)])		
	#	re-arrange so that ID1<ID2 
	tmp			<- subset(rplkl, ID1>ID2)
	setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
	set(tmp, NULL, 'TYPE', tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',TYPE)))])
	rplkl		<- rbind(subset(rplkl, ID1<=ID2), tmp)
	tmp			<- subset(rpw, ID1>ID2)
	setnames(tmp, c('ID1','ID2','ID1_L','ID1_R','ID2_L','ID2_R','PATHS_12','PATHS_21'), c('ID2','ID1','ID2_L','ID2_R','ID1_L','ID1_R','PATHS_21','PATHS_12'))
	set(tmp, NULL, 'TYPE', tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',TYPE)))])
	rpw			<- rbind(subset(rpw, ID1<=ID2), tmp)
	tmp			<- subset(rtp.todi2, ID1>ID2)
	setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
	set(tmp, NULL, 'TYPE', tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',TYPE)))])
	rtp.todi2	<- rbind(subset(rtp.todi2, ID1<=ID2), tmp)	
	#
	#	select runs with highest neff
	setkey(rplkl, ID1, ID2, PTY_RUN, GROUP, TYPE)	# sort to make sure same run selected in case of tie
	tmp			<- rplkl[ GROUP==linked.group & TYPE==linked.type.yes, list(PTY_RUN=PTY_RUN[which.max(NEFF)]), by=c('ID1','ID2')]
	rplkl		<- merge(rplkl, tmp, by=c('ID1','ID2','PTY_RUN'))
	rpw			<- merge(rpw, tmp, by=c('ID1','ID2','PTY_RUN'))
	rtp.todi2	<- merge(rtp.todi2, tmp, by=c('ID1','ID2','PTY_RUN'))
	#	save
	save(rtp.todi2, rplkl, rpw, file=gsub('\\.rda','_allwindows.rda',outfile))
	
	#
	#	prepare just the dwin and rplkl that we need for further linkage analysis of the pairs 
	#	this is all pairs for whom unlinked is not decisive
	#
	tmp			<- unique( subset(rtp.todi2, select=c(ID1, ID2, PTY_RUN)) )
	rplkl		<- merge(rplkl, tmp, by=c('ID1','ID2','PTY_RUN'))
	rpw			<- merge(rpw, tmp, by=c('ID1','ID2','PTY_RUN'))
	save(rp, rd, rh, ra, rs, rtp.todi2, rplkl, rpw, file=outfile)
	
	
	#
	#	construct max prob network among all possible pairs regardless of gender
	#	above NEFF cut
	rtn		<- subset(rtp.todi2, GROUP==linked.group & TYPE==linked.type.yes & NEFF>neff.cut)
	#	get transmission chains with igraph, to speed up calculations later
	tmp		<- subset(rtn, select=c(ID1, ID2))			
	tmp		<- graph.data.frame(tmp, directed=FALSE, vertices=NULL)
	rtc		<- data.table(ID=V(tmp)$name, CLU=clusters(tmp, mode="weak")$membership)	
	tmp2	<- rtc[, list(CLU_SIZE=length(ID)), by='CLU']
	setkey(tmp2, CLU_SIZE)
	tmp2[, IDCLU:=rev(seq_len(nrow(tmp2)))]
	rtc			<- subset( merge(rtc, tmp2, by='CLU'))
	rtc[, CLU:=NULL]
	setkey(rtc, IDCLU)
	#	add info on edges
	setnames(rtc, c('ID'), c('ID1'))
	rtn		<- merge(rtn, rtc, by='ID1')
	rtc[, CLU_SIZE:=NULL]
	setnames(rtc, c('ID1'), c('ID2'))
	rtn		<- merge(rtn, rtc, by=c('ID2','IDCLU'))
	tmp		<- subset(rtn, select=c(ID1, ID2, PTY_RUN, IDCLU, CLU_SIZE))
	#	add posterior mean for network types	 
	tmp		<- merge(tmp, subset(rplkl, GROUP==scores.group), by=c('ID1','ID2','PTY_RUN'))
	tmp		<- merge(tmp, tmp[, list(TYPE=TYPE, POSTERIOR_SCORE=(POSTERIOR_ALPHA-1)/sum(POSTERIOR_ALPHA-1) ), by=c('ID1','ID2','PTY_RUN')], by=c('ID1','ID2','PTY_RUN','TYPE'))
	rtn		<- rbind(tmp, rtn)
	#	generate maximum branch transmission network
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
	#
	#	work out prob for linkage in max prob network, when 'inconsistent direction' is ignored
	rtnn[, POSTERIOR_SCORE_LINKED_MECN:= pmax(MX_PROB_12,MX_PROB_21) + 0.5*POSTERIOR_SCORE]
	set(rtnn, NULL, c('POSTERIOR_SCORE','MX_PROB_12','MX_PROB_21'), NULL)
	#
	#	merge POSTERIOR_SCORE_LINKED on max prob network
	#	rationale: this describes prob of linkage. here, any 'inconsistent direction' is still considered as prob for linkage 	
	tmp		<- subset(rplkl, GROUP==linked.group & TYPE==linked.type.yes, c(ID1,ID2,PTY_RUN,POSTERIOR_ALPHA,POSTERIOR_BETA,N_TYPE))
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
	tmp		<- subset(tmp, !TYPE%in%scores.type.no)
	set(tmp, NULL, 'TYPE', tmp[, paste0('NETWORK_SCORE_',TYPE)])	
	tmp		<- dcast.data.table(tmp, ID1+ID2+PTY_RUN~TYPE, value.var='POSTERIOR_SCORE')
	rtnn	<- merge(rtnn, tmp, by=c('ID1','ID2','PTY_RUN'), all.x=TRUE)
	#	ensure DIR scores and NETWORK_SCORE scores are compatible with self-consistency in maxprobnetwork
	tmp		<- rtnn[, which(LINK_12==0 & LINK_21==1 & POSTERIOR_SCORE_12>POSTERIOR_SCORE_21)]
	set(rtnn, tmp, c('POSTERIOR_SCORE_12','NETWORK_SCORE_12'), 0)
	tmp		<- rtnn[, which(LINK_12==1 & LINK_21==0 & POSTERIOR_SCORE_21>POSTERIOR_SCORE_12)]
	set(rtnn, tmp, c('POSTERIOR_SCORE_21','NETWORK_SCORE_21'), 0)	
	rtnn	<- subset(rtnn, LINK_12==1 | LINK_21==1)		
	#
	
	#
	#	define the phylogenetic relationship between those pairs for whom linkage was not excluded 
	#
	rtp			<- unique(subset(rtn, select=c(ID1,ID2,PTY_RUN)))
	#	merge linkage and direction probabilities
	rtp			<- merge(rtp, rtnn, by=c('ID1','ID2','PTY_RUN')) #merge data on unlinked for everyone
	tmp			<- subset(rplkl, GROUP==linked.group & TYPE==linked.type.no)	
	tmp[, POSTERIOR_SCORE_UNLINKED:= (POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE)]
	rtp			<- merge(rtp, subset(tmp, select=c(ID1, ID2, PTY_RUN, POSTERIOR_SCORE_UNLINKED)), all.x=1, by=c('ID1','ID2','PTY_RUN'))	
	#	define SELECT
	rtp[, SELECT:= NA_character_]
	set(rtp, rtp[, which(is.na(PTY_RUN))], 'SELECT', 'insufficient deep sequence data for at least one partner of couple')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED>conf.cut)], 'SELECT', 'couple most likely not a pair')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED<=conf.cut & is.na(LINK_12) & is.na(LINK_21))], 'SELECT', 'couple ambiguous if pair or not pair')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED<=conf.cut & !is.na(LINK_12) & !is.na(LINK_21) & POSTERIOR_SCORE_LINKED<=conf.cut)], 'SELECT', 'couple ambiguous if pair or not pair')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED<=conf.cut & !is.na(LINK_12) & !is.na(LINK_21) & POSTERIOR_SCORE_LINKED>conf.cut)], 'SELECT', 'couple most likely a pair direction not resolved')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED<=conf.cut & !is.na(LINK_12) & !is.na(LINK_21) & POSTERIOR_SCORE_LINKED>conf.cut & POSTERIOR_SCORE_12>conf.cut)], 'SELECT', 'couple most likely a pair direction resolved to 12')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED<=conf.cut & !is.na(LINK_12) & !is.na(LINK_21) & POSTERIOR_SCORE_LINKED>conf.cut & POSTERIOR_SCORE_21>conf.cut)], 'SELECT', 'couple most likely a pair direction resolved to 21')	
	setkey(rtp, ID1, ID2)
	
	save(rtp.todi2, rplkl, rpw, rtn, rtnn, rtp, file=gsub('\\.rda','_networksallpairs.rda',outfile))
	
	if(0)
	{
		#	plot transmission networks
		df	<- copy(rtn)	
		di	<- melt(subset(df, select=c(ID1,ID2,PTY_RUN)), id.vars='PTY_RUN', value.name='ID')
		set(di, NULL, 'variable', NULL)
		di	<- unique(di)
		p	<- phsc.plot.probability.network(df, di, point.size=10, 
				edge.gap=0.015, 
				edge.size=0.4, 
				curvature= -0.2, 
				arrow=arrow(length=unit(0.04, "npc"), type="open"), 
				curv.shift=0.02, 
				label.size=3, 
				threshold.linked=0.6)	
		pdf(file=paste0(outfile.base,'transmissionnetwork.pdf'), w=16, h=16)		
		print(p)
		dev.off()
	}
	
	if(0)
	{
		#
		#	plot most likely transmission chain
		df	<- copy(rtnn)	
		di	<- melt(subset(df, select=c(ID1,ID2,PTY_RUN)), id.vars='PTY_RUN', value.name='ID')
		set(di, NULL, 'variable', NULL)
		di	<- unique(di) 
		q		<- phsc.plot.max.probability.network(df, di, point.size=10, 
				edge.gap=0.015, 
				edge.size=0.4, 
				curvature= -0.2, 
				arrow=arrow(length=unit(0.04, "npc"), type="open"), 
				curv.shift=0.02, 
				label.size=3, 
				threshold.linked=0.6,
				layout=p$layout)	
		pdf(file=paste0(outfile.base,'mostlikelychain.pdf'), w=16, h=16)	
		print(q)
		dev.off()
	}
	
	
	#
	#	window scan for all pairs in rtn
	#
	if(1)
	{
		plot.file	<- paste0(outfile.base,'scanplot.pdf')
		rpw2		<- merge(rpw, subset(rtp, select=c(ID1, ID2, PTY_RUN)), by=c('ID1','ID2','PTY_RUN'))
		phsc.plot.windowscan.for.pairs(rpw2, plot.file, plot.w=10, plot.h=35, id.cols=c('ID1','ID2'), ylim=NULL, cols.typet=NULL)		
	}	
}

TasP.calculate.bam.readlength.coverage<- function()
{
	require(ggplot2)
	require(data.table)
	require(Rsamtools)		
	
	pty.data.dir	<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/data'			
	outfile			<- file.path(pty.data.dir,'bam_stats_170928.rda')			
	
	bfiles			<- data.table(FILE=list.files(pty.data.dir, pattern='bam$'))	
	bfiles			<- subset(bfiles, !grepl('Contam',FILE) & !grepl('merged',FILE))
	bfiles[, FILE_ID:=gsub('.bam','',FILE)]	
	#FILE<- "14939_1_80.bam"
	#cat('\nreading',FILE)
	#bam<- scanBam(file.path(pty.data.dir,FILE))[[1]]	#scan entire file
	bam.cov	<- bfiles[,{
				cat('\nCOV reading',FILE)
				z	<- scanBam(file.path(pty.data.dir,FILE), param=ScanBamParam(what=c('qwidth','pos','rname')))[[1]]
				tmp	<- IRanges(start=z$pos, width=z$qwidth)
				tmp <- coverage(tmp)
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths, REF=z$rname[1])
			}, by='FILE_ID']
	#	get lengths of all reads in quality trimmed bam file
	#z	<- scanBam(file.path(pty.data.dir,FILE), param=ScanBamParam(what=c('qwidth','qual')))
	bam.len			<- bfiles[,{
				#FILE<- "14939_1_80.bam"
				cat('\nLEN reading',FILE)
				z	<- scanBam(file.path(pty.data.dir,FILE), param=ScanBamParam(what=c('qwidth')))
				list(QU=seq(0,320,20), CDF=ecdf(z[[1]][['qwidth']])(seq(0,320,20)))
				#list(PR=seq(0.01,0.99,0.01), QU=quantile(z[[1]][['qwidth']], p=seq(0.01,0.99,0.01)))				
			}, by='FILE']
	#
	save(bam.len, bam.cov, file= outfile)
}

TasP.evaluate.window.len<- function()
{	
	tip.regex	<- '(.*)_read_([0-9]+)_count_([0-9]+)'
	min.count	<- 30
	indir		<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch'
	plot.file	<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/171001_evaluate_winlen.pdf'
	infiles		<- data.table(F=list.files(indir, recursive=TRUE, full.names=TRUE, pattern='subgraphs_s_ptyr1_InWindow_.*.rda$'))
	infiles[, WLEN:= as.integer(gsub('.*_winlen([0-9]+).*','\\1',F))]
	infiles[, WFROM:= as.integer(gsub('.*_InWindow_([0-9]+).*','\\1',F))]
	infiles[, PTY_RUN:= as.integer(gsub('.*ptyr([0-9]+)_.*','\\1',F))]
	
	
	dw			<- infiles[, {
				#F			<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/ptyr1_winlen180/subgraphs_s_ptyr1_InWindow_2300_to_2479.rda'
				load(F)
				tmp			<- data.table(TAXA=tree$tip.label)
				tmp			<- subset(tmp, grepl(tip.regex, TAXA))
				tmp[, ID:= gsub(tip.regex,'\\1',TAXA)]
				tmp[, READ:= as.integer(gsub(tip.regex,'\\2',TAXA))]
				tmp[, COUNT:= as.integer(gsub(tip.regex,'\\3',TAXA))]				
				tmp			<- tmp[, list(READ=length(unique(READ)), COUNT=sum(COUNT)),by='ID']
				tmp
			}, by=c('PTY_RUN','WLEN','WFROM')]
	
	tmp			<- subset(dw, COUNT>=min.count)[, list(IDN=length(unique(ID))), by=c('PTY_RUN','WLEN','WFROM')]
	tmp			<- tmp[, list(IDN_L=min(IDN), IDN_M=mean(IDN), IDN_U=max(IDN)), by=c('PTY_RUN','WLEN')]
	
	ggplot(tmp, aes(x=WLEN, y=IDN_M, ymin=IDN_L, ymax=IDN_U)) + 
			geom_point() + 
			geom_errorbar() + 
			scale_y_continuous(lim=c(0,75), breaks=seq(0,75,5), expand=c(0,0)) +
			labs(x='window length', y='number of individuals with min 30 reads')
	ggsave(file=plot.file, w=5, h=6)	
}

TasP.evaluate.bam.len<- function()
{	
	infile.bam		<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/output/bam_stats_170928.rda'
	infile.runs		<- '/Users/Oliver/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/input/ptyRuns_2batchesTest_OR.csv'
	plot.file		<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/output/bam_len_170928.pdf'
	load(infile.bam)
	
	pty.runs		<- as.data.table(read.csv( infile.runs ))
	
	setnames(bam.len, 'FILE', 'SAMPLE_ID')
	set(bam.len,NULL,'SAMPLE_ID',bam.len[, gsub('\\.bam','',SAMPLE_ID)])
	#bam.len			<- bam.len[, list(X=paste0(QU[-length(QU)],'-',QU[-1]-1), PDF=diff(CDF)), by='SAMPLE_ID']
	
	pty.runs		<- unique(pty.runs, by=c('PTY_RUN','SAMPLE_ID'))
	bam.len			<- merge(bam.len, pty.runs, by='SAMPLE_ID')
	setkey(bam.len, PTY_RUN, SAMPLE_ID, QU)
	
	ggplot(subset(bam.len, QU>=40 & QU<320), aes(y=CDF, x=factor(QU), fill=factor(PTY_RUN))) + 
			geom_boxplot() + 
			scale_fill_brewer(palette='Set1') +
			theme_bw() + 
			labs(y='cumulative frequency\nin one individual\n', x='\nlength of quality-trimmed short reads\n(nt)', fill='sequence run') +
			theme(legend.position='bottom')
	ggsave(file=plot.file, w=8, h=6)
}

TasP.evaluate.bam.coverage<- function()
{	
	infile.bam		<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/output/bam_stats_170928.rda'
	infile.runs		<- '/Users/Oliver/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/input/ptyRuns_2batchesTest_OR.csv'
	plot.file		<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/output/bam_coverage_170928.pdf'
	load(infile.bam)
	
	pty.runs		<- as.data.table(read.csv( infile.runs ))
	
	setnames(bam.cov, 'FILE_ID', 'SAMPLE_ID')
	pty.runs		<- unique(pty.runs, by=c('PTY_RUN','SAMPLE_ID'))
	bam.cov			<- merge(bam.cov, pty.runs, by='SAMPLE_ID')
	setkey(bam.cov, PTY_RUN, SAMPLE_ID, POS)
	
	ps				<- lapply(bam.cov[, unique(PTY_RUN)], function(ptyr)
			{
				cat('\nprocess run', ptyr)
				tmp				<- subset(bam.cov, PTY_RUN==ptyr)
				tmp				<- tmp[, {
							z	<- rep(COV,REP)
							list(COV=z, POS=seq_along(z), REF=REF[1])
						}, by=c('SAMPLE_ID','PTY_RUN')]
				p	<- ggplot(tmp, aes(x=POS, y=COV, colour=SAMPLE_ID)) + geom_step() + theme_bw() +
						scale_y_log10() +
						facet_grid(~PTY_RUN) +
						labs(x='genome position\n(relative to reference)', y='read coverage', colour='patient')
				p
			})
	#pdf(file=gsub('\\.rda','_coverage_UG.pdf',infile), w=12, h=5)
	pdf(file=plot.file, w=12, h=5)
	for(i in seq_along(ps))
		print(ps[[i]])	
	dev.off()
}

Tasp.pipeline.testbatch.170925.stage1<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	#
	#	INPUT ARGS PLATFORM
	#	
	if(0)
	{
		#fix Tiago's file
		HOME				<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch'
		in.dir				<- file.path(HOME,"input")
		pty.runs			<- as.data.table(read.csv( file.path(in.dir, 'ptyRuns_2batchesTest.csv') ))
		set(pty.runs, NULL, c('PTY_RUN'), NULL)		
		setnames(pty.runs, c('BATCH'), c('PTY_RUN'))
		set(pty.runs, NULL, 'REFERENCE_ID', pty.runs[, paste0(REFERENCE_ID,'.fasta')])		
		set(pty.runs, NULL, 'UNIT_ID', pty.runs[, paste0('TasP',UNIT_ID)])
		write.csv(pty.runs, row.names=FALSE, file=file.path(in.dir, 'ptyRuns_2batchesTest_OR.csv'))		
	}
	if(1)
	{	
		#olli's paths
		#HOME				<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch'
		in.dir				<- file.path(HOME,"input")
		work.dir			<- file.path(HOME,"tmp")		
		out.dir				<- file.path(HOME,"output")				
		pty.runs			<- as.data.table(read.csv( file.path(in.dir, 'ptyRuns_2batchesTest_OR.csv'), header=TRUE, stringsAsFactors=FALSE ))
		hpc.load			<- "module load intel-suite/2015.1 mpi R/3.3.2 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"
		#hpc.nproc			<- 4	
		hpc.nproc			<- 1
		#prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner_make_trees.py'
		#pty.data.dir		<- '/work/or105/PANGEA_mapout/data'
		prog.pty			<- '/Users/Oliver/git/phylotypes/phyloscanner_make_trees.py'
		pty.data.dir		<- file.path(HOME,"data")
		#prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-SSE3 -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-SSE3 -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))
		prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-AVX -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-AVX -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))		
		pty.select			<- 1			
	}					
	#
	#	INPUT ARGS PHYLOSCANNER RUN
	#	
	if(1)
	{				
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft='mafft', 
				prog.raxml=prog.raxml, 
				data.dir=pty.data.dir, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=NA,	#to be specified as part of pty.runs --> BACKGROUND_ID
				alignments.root='B.FR.83.HXB2.K03455', 
				alignments.pairwise.to='B.FR.83.HXB2.K03455',
				bl.normalising.reference.file=system.file(package="Phyloscanner.R.utilities", "data", "hiv.hxb2.norm.constants.rda"),
				bl.normalising.reference.var='MEDIAN_PWD',
				tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$',		
				window.automatic= '', 
				merge.threshold=2, 
				min.read.count=1, 				
				merge.paired.reads=TRUE, 
				#no.trees=TRUE,
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				dont.check.recombination=TRUE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				#win=c(2300,2800,125,250),
				#win=c(2300,2800,75,150),
				#win=c(2300,2800,100,200),
				#win=c(2300,2800,100,180),
				#win=c(2300,2800,95,190),
				win=c(2300,2800,85,170),			
				keep.overhangs=FALSE,	
				use.blacklisters=c('ParsimonyBasedBlacklister','DownsampleReads'),				
				roguesubtree.kParam=20,
				roguesubtree.prop.threshold=0,
				roguesubtree.read.threshold=20,
				dwns.maxReadsPerPatient=50,	
				multifurcation.threshold=1e-5,
				split.pruneBlacklist=FALSE,
				trms.allowMultiTrans=TRUE,				
				split.rule='s',
				split.kParam=20,
				split.proximityThreshold=0.035,	
				split.readCountsMatterOnZeroBranches=TRUE,
				pw.trmw.min.reads=20,									
				pw.trmw.min.tips=1,
				pw.trmw.close.brl=0.035,
				pw.trmw.distant.brl=0.08,
				pw.prior.keff=2,
				pw.prior.neff=3,
				pw.prior.keff.dir=2,
				pw.prior.neff.dir=3,				
				pw.prior.calibrated.prob=0.66,
				mem.save=0,
				verbose=TRUE,
				select=pty.select
		)		 
	}	
	#
	#	RUN PHYLOSCANNER
	#
	if(1)
	{
		pty.c				<- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)		
		#pty.c[1,cat(CMD)]		
		invisible(pty.c[,	{
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=171, hpc.q="pqeelab", hpc.mem="23600mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=171, hpc.q="pqeelab", hpc.mem="5600mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=24, hpc.q=NA, hpc.mem="1890mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=199, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile		<- paste("scRA2",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
							#stop()
						}, by='PTY_RUN'])
		quit('no')
	}	
	if(0)	#check failing runs
	{
		df	<- data.table(F=readLines('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/all_files.txt'))
		df[, PTY_RUN:= as.integer(gsub('ptyr([0-9]+)_.*','\\1',F))]
		tmp	<- df[, list(DONE=any(grepl('pairwise',F))), by='PTY_RUN']
		tmp	<- subset(tmp, !DONE)
		merge(tmp, subset(pty.runs, SID=='15080_1_28', c('PTY_RUN','SID')), by='PTY_RUN', all=TRUE)
		
		subset(tmp, !DONE)[, cat(paste(sort(PTY_RUN), collapse='","'))]
	}
}