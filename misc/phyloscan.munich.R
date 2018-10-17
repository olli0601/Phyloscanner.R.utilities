Munich.networks.plot.180924<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(igraph)
	require(ggnet)
	require(Phyloscanner.R.utilities)
	
	infile				<- '~/Box Sync/OR_Work/2018/2018_MunichCluster/180924_analysis/180924_MunichCluster_w250_cl25_d50_min30_cut60_networksallpairs.rda'
	infile				<- "/Users/Oliver/Box Sync/OR_Work/2018/2018_MunichCluster/180924_analysis/180924_MunichCluster_w250_cl25_d50_min30_cut60_strongsupport_networksallpairs.rda"
	outfile.base		<- gsub('networksallpairs.rda','',infile)		
	
	#	load preprocessed phyloscanner output
	load(infile)
	#	IDs of networks
	idclus	<- sort(unique(rtn$IDCLU))
	
	
	#	get proportion of windows with GD<1/250 for all pairs in network
	tmp		<- unique(subset(rtn, select=c(ID1,ID2)))
	rpw2	<- merge(tmp, rpw, by=c('ID1','ID2'))
	dvcl	<- rpw2[, list(PVCL= mean(PATRISTIC_DISTANCE<1/250), N=length(PATRISTIC_DISTANCE)), by=c('ID1','ID2')]
	#ggplot(dvcl, aes(x=PVCL)) + geom_histogram(binwidth=0.05)
	dvcl[, PVCLC:= cut(PVCL, breaks=c(-Inf,0.85,0.95,1.05),labels=c('<85%','85-95%','>95%'))]
	
	#
	#	plot max edge network with edges coloured by %identical reads
	pns		<- lapply(seq_along(idclus), function(i)
			{
				idclu	<- idclus[i]								
				df		<- subset(rtn, IDCLU==idclu)
				df		<- merge(df, subset(dvcl, select=c(ID1,ID2,PVCLC)),by=c('ID1','ID2'))
				di		<- unique(melt(subset(df, select=c(ID1, ID2)), measure.vars=c('ID1','ID2'), value.name='ID'), by='ID')
				di[, variable:=NULL]
				di[, TYPE:= 'IDU']
				di[, TYPE2:= 'IDU']
				p		<- phsc.plot.maxedge.network(df, di, point.size=10, 
						edge.gap=0.02, 
						edge.size=0.4, 
						curvature= -0.2, 
						arrow=arrow(length=unit(0.02, "npc"), type="open"), 
						curv.shift=0.08, 
						label.size=3, 	
						edge.label='PVCLC',
						edge.label.values=c('<85%'='grey50','85-95%'='darkorange','>95%'='red'),
						node.label='ID',
						node.shape='TYPE',
						node.fill='TYPE2',  
						node.fill.values=c('IDU'='steelblue2'),
						node.shape.values=c('IDU'=21),
						threshold.linked=0.6)
			})
	pdf(file=paste0(outfile.base,'maxedgenetworks_colourednearidentical.pdf'), w=28, h=28)
	for(i in seq_along(pns))	
		print(pns[[i]])
	dev.off()
	#
	#	plot probability network
	pns		<- lapply(seq_along(idclus), function(i)
			{
				idclu	<- idclus[i]				
				df		<- subset(rtn, IDCLU==idclu)
				di		<- unique(melt(subset(df, select=c(ID1, ID2)), measure.vars=c('ID1','ID2'), value.name='ID'), by='ID')
				di[, variable:=NULL]
				di[, TYPE:= 'IDU']
				di[, TYPE2:= 'IDU']
				p		<- phsc.plot.probability.network(df, di, point.size=10, 
						edge.gap=0.04, 
						edge.size=0.4, 
						curvature= -0.2, 
						arrow=arrow(length=unit(0.04, "npc"), type="open"), 
						curv.shift=0.08, 
						label.size=3, 						 
						node.label='ID',
						node.shape='TYPE',
						node.fill='TYPE2',  
						node.fill.values=c('IDU'='steelblue2'),
						node.shape.values=c('IDU'=16),
						threshold.linked=0.6)
				p	
			})
	pdf(file=paste0(outfile.base,'probabilitynetworks.pdf'), w=28, h=28)
	for(i in seq_along(pns))	
		print(pns[[i]])
	dev.off()
	
	
	#
	#	plot max edge network
	pns		<- lapply(seq_along(idclus), function(i)
			{
				idclu	<- idclus[i]				
				df		<- subset(rtn, IDCLU==idclu)
				di		<- unique(melt(subset(df, select=c(ID1, ID2)), measure.vars=c('ID1','ID2'), value.name='ID'), by='ID')
				di[, variable:=NULL]
				di[, TYPE:= 'IDU']
				di[, TYPE2:= 'IDU']
				p		<- phsc.plot.maxedge.network(df, di, point.size=10, 
						edge.gap=0.02, 
						edge.size=0.4, 
						curvature= -0.2, 
						arrow=arrow(length=unit(0.02, "npc"), type="open"), 
						curv.shift=0.08, 
						label.size=3, 						 
						node.label='ID',
						node.shape='TYPE',
						node.fill='TYPE2',  
						node.fill.values=c('IDU'='steelblue2'),
						node.shape.values=c('IDU'=16),
						threshold.linked=0.6)
				p	
			})
	pdf(file=paste0(outfile.base,'maxedgenetworks.pdf'), w=28, h=28)
	for(i in seq_along(pns))	
		print(pns[[i]])
	dev.off()
	
	#
	#	plot most likely transmission network	
	qns		<- lapply(seq_along(idclus), function(i)
			{				
				idclu	<- idclus[i]
				layout	<- pns[[i]]$layout 
				print(layout)
				df		<- subset(rtnn, IDCLU==idclu)
				di		<- unique(melt(subset(df, select=c(ID1, ID2)), measure.vars=c('ID1','ID2'), value.name='ID'), by='ID')
				di[, variable:=NULL]
				di[, TYPE:= 'IDU']
				di[, TYPE2:= 'IDU']				
				p		<- phsc.plot.max.probability.network(df, di, point.size=10, 
						edge.gap=0.03, 
						edge.size=0.4, 
						curvature= -0.2, 
						arrow=arrow(length=unit(0.04, "npc"), type="open"), 
						curv.shift=0.08, 
						label.size=3, 
						node.label='ID', 
						node.shape='TYPE', 
						node.fill='TYPE2', 
						node.shape.values=c('IDU'=16), 
						node.fill.values=c('IDU'='steelblue2'),
						threshold.linked=0.6,
						layout=layout)
				p		
			})	
	pdf(file=paste0(outfile.base,'maxprobabilitynetworks.pdf'), w=20, h=20)	
	for(i in seq_along(qns))	
		print(qns[[i]])
	dev.off()
}

Munich.phyloscan.plots.180924<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	
	neff.cut				<- 3	
	if(0)
	{
		linked.group	<- 'TYPE_PAIR_TODI2'
		linked.type.yes	<- 'linked'		
		linked.type.no	<- 'unlinked'
		dir.group		<- 'TYPE_DIR_TODI2'		
	}
	if(1)
	{
		linked.group	<- 'TYPE_CHAIN_TODI'
		linked.type.yes	<- 'chain'
		linked.type.no	<- 'distant'
		dir.group		<- 'TYPE_ADJ_DIR_TODI2'		
	}
	
	infile.trmpairs.todi	<- '~/Box Sync/OR_Work/2018/2018_MunichCluster/180924_analysis/180924_MunichCluster_w250_cl25_d50_min30_cut60_networksallpairs.rda'
	infile.trmpairs.todi	<- "/Users/Oliver/Box Sync/OR_Work/2018/2018_MunichCluster/180924_analysis/180924_MunichCluster_w250_cl25_d50_min30_cut60_strongsupport_networksallpairs.rda"
	outfile.base			<- gsub('networksallpairs.rda','',infile.trmpairs.todi)
	
	#	read and update confidence cut
	confidence.cut			<- 0.6 
	if(grepl('_cut[0-9]+',infile.trmpairs.todi))
		confidence.cut		<- as.integer(gsub('^.*_cut([0-9]+).*$','\\1',infile.trmpairs.todi))/100
	
	#	load output from preprocessing step
	load(infile.trmpairs.todi)
	
	#	merge linkage and direction probabilities in max prob network
	rtp			<- unique(subset(rtn, select=c(ID1,ID2,PTY_RUN)))				
	#	merge data on unlinked for everyone
	tmp			<- subset(rplkl, GROUP==linked.group & TYPE==linked.type.no)	
	tmp[, POSTERIOR_SCORE_UNLINKED:= (POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE)]
	rtp			<- merge(rtp, subset(tmp, select=c(ID1, ID2, PTY_RUN, POSTERIOR_SCORE_UNLINKED)), all.x=1, by=c('ID1','ID2','PTY_RUN'))	
	#	merge data on linked for everyone
	tmp			<- subset(rplkl, GROUP==linked.group & TYPE==linked.type.yes)	
	tmp[, POSTERIOR_SCORE_LINKED:= (POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE)]
	rtp			<- merge(rtp, subset(tmp, select=c(ID1, ID2, PTY_RUN, POSTERIOR_SCORE_LINKED)), all.x=1, by=c('ID1','ID2','PTY_RUN'))
	#	merge data on dir 12 for everyone
	tmp			<- subset(rplkl, GROUP==dir.group & TYPE=='12')	
	tmp[, POSTERIOR_SCORE_12:= (POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE)]
	rtp			<- merge(rtp, subset(tmp, select=c(ID1, ID2, PTY_RUN, POSTERIOR_SCORE_12)), all.x=1, by=c('ID1','ID2','PTY_RUN'))
	#	merge data on dir 21 for everyone
	tmp			<- subset(rplkl, GROUP==dir.group & TYPE=='21')	
	tmp[, POSTERIOR_SCORE_21:= (POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE)]
	rtp			<- merge(rtp, subset(tmp, select=c(ID1, ID2, PTY_RUN, POSTERIOR_SCORE_21)), all.x=1, by=c('ID1','ID2','PTY_RUN'))
	#	set POSTERIOR_SCORE for direction to 0 if most likely edge type is ambiguous 
	tmp			<- rtn[, list(TYPE_MAX=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
	rtp			<- merge(rtp, tmp, by=c('ID1','ID2') )
	set(rtp, rtp[, which(TYPE_MAX=='ambiguous')], c('POSTERIOR_SCORE_21','POSTERIOR_SCORE_12'), 0)
	
	#	classify linkages and direction
	rtp[, SELECT:= NA_character_]
	set(rtp, rtp[, which(is.na(PTY_RUN))], 'SELECT', 'insufficient deep sequence data for at least one partner of couple')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED>confidence.cut)], 'SELECT', 'couple most likely not a pair')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED<=confidence.cut & POSTERIOR_SCORE_LINKED<=confidence.cut)], 'SELECT', 'couple ambiguous if pair or not pair')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED<=confidence.cut & POSTERIOR_SCORE_LINKED>confidence.cut)], 'SELECT', 'couple most likely a pair direction not resolved')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED<=confidence.cut & POSTERIOR_SCORE_LINKED>confidence.cut & POSTERIOR_SCORE_12>confidence.cut)], 'SELECT', 'couple most likely a pair direction resolved to 12')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED<=confidence.cut & POSTERIOR_SCORE_LINKED>confidence.cut & POSTERIOR_SCORE_21>confidence.cut)], 'SELECT', 'couple most likely a pair direction resolved to 21')	
	setkey(rtp, ID1, ID2)
	#  SELECT
	#  couple most likely a pair direction not resolved 	couple most likely a pair direction resolved to 12	couple most likely a pair direction resolved to 21 
	#  157                                                  11                                        			21  
	
	
	#
	#	make phyloscan plots	
	rps		<- subset(rtp, !grepl('ambiguous|couple most likely a pair direction not resolved', SELECT), select=c(ID1, ID2, PTY_RUN, SELECT, POSTERIOR_SCORE_UNLINKED, POSTERIOR_SCORE_LINKED, POSTERIOR_SCORE_12, POSTERIOR_SCORE_21))
	setkey(rps, SELECT, POSTERIOR_SCORE_LINKED)
	write.csv(rps, file=paste0(outfile.base,'_summary_lklpairs_withdir.csv'))
	rps[, DUMMY:=seq_len(nrow(rps))]
	rps[, LABEL:=rps[, factor(DUMMY, levels=DUMMY, labels=paste0('id1 ',ID1,' id2 ', ID2,'\n', SELECT,'\nunlinked: ',round(POSTERIOR_SCORE_UNLINKED, d=3), ' linked: ', round(POSTERIOR_SCORE_LINKED, d=3), ' 12: ', round(POSTERIOR_SCORE_12, d=3), ' 21: ', round(POSTERIOR_SCORE_21, d=3),'\nrun',PTY_RUN))]]	
	rpw2		<- merge(rpw, unique(subset(rps, select=c(ID1,ID2,PTY_RUN))), by=c('ID1', 'ID2','PTY_RUN'))
	plot.file	<- paste0(outfile.base,'_phyloscans_lklpairs_withdir.pdf')
	phsc.plot.windowscan.for.pairs(rpw2, plot.file, plot.w=10, plot.h=40, id.cols=c('ID1','ID2'), ylim=NULL, cols.typet=NULL)
	
	rps		<- subset(rtp, grepl('couple most likely a pair direction not resolved', SELECT), select=c(ID1, ID2, PTY_RUN, SELECT, POSTERIOR_SCORE_UNLINKED, POSTERIOR_SCORE_LINKED, POSTERIOR_SCORE_12, POSTERIOR_SCORE_21))
	setkey(rps, SELECT, POSTERIOR_SCORE_LINKED)
	write.csv(rps, file=paste0(outfile.base,'_summary_lklpairs_nodir.csv'))
	rps[, DUMMY:=seq_len(nrow(rps))]
	rps[, LABEL:=rps[, factor(DUMMY, levels=DUMMY, labels=paste0('id1 ',ID1,' id2 ', ID2,'\n', SELECT,'\nunlinked: ',round(POSTERIOR_SCORE_UNLINKED, d=3), ' linked: ', round(POSTERIOR_SCORE_LINKED, d=3), ' 12: ', round(POSTERIOR_SCORE_12, d=3), ' 21: ', round(POSTERIOR_SCORE_21, d=3),'\nrun',PTY_RUN))]]	
	rpw2		<- merge(rpw, unique(subset(rps, select=c(ID1,ID2,PTY_RUN))), by=c('ID1', 'ID2','PTY_RUN'))
	plot.file	<- paste0(outfile.base,'_phyloscans_lklpairs_nodir.pdf')
	phsc.plot.windowscan.for.pairs(rpw2, plot.file, plot.w=10, plot.h=300, id.cols=c('ID1','ID2'), ylim=NULL, cols.typet=NULL)
	
	
	if(0)
	{
		require(colorspace)
		#for(ii in seq_len(nrow(rtpdm))[-1])
		for(ii in 111:252)
		{		
			if(rtpdm[ii, PTY_RUN]!=1)
			{
				indir		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun'
				# load dfr and phs
				load( file.path(indir, paste0('ptyr',rtpdm[ii,PTY_RUN],'_trees.rda')) )
				# setup plotting
				ids			<- c(rtpdm[ii, MALE_RID],rtpdm[ii, FEMALE_RID])
				dfs			<- subset(dfr, select=c(W_FROM, W_TO, IDX))
				dfs[, MALE_RID:=ids[1]]
				dfs[, FEMALE_RID:=ids[2]]
				dfs			<- merge(dfs, subset(rpw, GROUP=='TYPE_RAW'), by=c('MALE_RID','FEMALE_RID','W_FROM','W_TO'), all.x=1)	
				dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', PTY_RUN, '\nwindow ', W_FROM,'-', W_TO,'\npmf',PATHS_MF,' pfm',PATHS_FM, ' ',ADJACENT,' ',CONTIGUOUS,' ',TYPE, '\n', round(PATRISTIC_DISTANCE, d=5), sep='')]]
				plot.file	<- paste0(outfile.base, 'trees/todi_pairs_170516_run_', rtpdm[ii, PTY_RUN],'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
				invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10, tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'))			
			}				
		}
	}				
}

Munich.phyloscan.plots.on.MLE.network.180924<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	
	neff.cut				<- 3	
	if(0)
	{
		linked.group	<- 'TYPE_PAIR_TODI2'
		linked.type.yes	<- 'linked'		
		linked.type.no	<- 'unlinked'
		dir.group		<- 'TYPE_DIR_TODI2'		
	}
	if(1)
	{
		linked.group	<- 'TYPE_CHAIN_TODI'
		linked.type.yes	<- 'chain'
		linked.type.no	<- 'distant'
		dir.group		<- 'TYPE_ADJ_DIR_TODI2'		
	}
	
	infile.trmpairs.todi	<- '~/Box Sync/OR_Work/2018/2018_MunichCluster/180924_analysis/180924_MunichCluster_w250_cl25_d50_min30_cut60_networksallpairs.rda'
	infile.trmpairs.todi	<- "/Users/Oliver/Box Sync/OR_Work/2018/2018_MunichCluster/180924_analysis/180924_MunichCluster_w250_cl25_d50_min30_cut60_strongsupport_networksallpairs.rda"
	outfile.base			<- gsub('networksallpairs.rda','',infile.trmpairs.todi)
	
	#	read and update confidence cut
	confidence.cut			<- 0.6 
	if(grepl('_cut[0-9]+',infile.trmpairs.todi))
		confidence.cut		<- as.integer(gsub('^.*_cut([0-9]+).*$','\\1',infile.trmpairs.todi))/100
	
	#	load output from preprocessing step
	load(infile.trmpairs.todi)
	
	#	merge linkage and direction probabilities in max prob network
	rtp			<- unique(subset(rtn, select=c(ID1,ID2,PTY_RUN)))	
	rtp			<- merge(rtp, rtnn, by=c('ID1','ID2','PTY_RUN'), all.x=TRUE)	
	
	#	merge data on unlinked for everyone
	tmp			<- subset(rplkl, GROUP==linked.group & TYPE==linked.type.no)	
	tmp[, POSTERIOR_SCORE_UNLINKED:= (POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE)]
	rtp			<- merge(rtp, subset(tmp, select=c(ID1, ID2, PTY_RUN, POSTERIOR_SCORE_UNLINKED)), all.x=1, by=c('ID1','ID2','PTY_RUN'))	
	
	#	classify linkages and direction
	rtp[, SELECT:= NA_character_]
	set(rtp, rtp[, which(is.na(PTY_RUN))], 'SELECT', 'insufficient deep sequence data for at least one partner of couple')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED>confidence.cut)], 'SELECT', 'couple most likely not a pair')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED<=confidence.cut & is.na(LINK_12) & is.na(LINK_21))], 'SELECT', 'couple ambiguous if pair or not pair')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED<=confidence.cut & !is.na(LINK_12) & !is.na(LINK_21) & POSTERIOR_SCORE_LINKED<=confidence.cut)], 'SELECT', 'couple ambiguous if pair or not pair')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED<=confidence.cut & !is.na(LINK_12) & !is.na(LINK_21) & POSTERIOR_SCORE_LINKED>confidence.cut)], 'SELECT', 'couple most likely a pair direction not resolved')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED<=confidence.cut & !is.na(LINK_12) & !is.na(LINK_21) & POSTERIOR_SCORE_LINKED>confidence.cut & POSTERIOR_SCORE_12>confidence.cut)], 'SELECT', 'couple most likely a pair direction resolved to 12')
	set(rtp, rtp[, which(!is.na(PTY_RUN) & POSTERIOR_SCORE_UNLINKED<=confidence.cut & !is.na(LINK_12) & !is.na(LINK_21) & POSTERIOR_SCORE_LINKED>confidence.cut & POSTERIOR_SCORE_21>confidence.cut)], 'SELECT', 'couple most likely a pair direction resolved to 21')	
	setkey(rtp, ID1, ID2)
	
	
	#
	#	make phyloscan plots	
	rps		<- subset(rtp, !grepl('ambiguous', SELECT), select=c(ID1, ID2, PTY_RUN, SELECT, POSTERIOR_SCORE_UNLINKED, POSTERIOR_SCORE_LINKED, POSTERIOR_SCORE_12, POSTERIOR_SCORE_21))
	setkey(rps, SELECT, POSTERIOR_SCORE_LINKED)
	write.csv(rps, file=paste0(outfile.base,'_summary_lklpairs.csv'))
	rps[, DUMMY:=seq_len(nrow(rps))]
	rps[, LABEL:=rps[, factor(DUMMY, levels=DUMMY, labels=paste0('id1 ',ID1,' id2 ', ID2,'\n', SELECT,'\nunlinked: ',round(POSTERIOR_SCORE_UNLINKED, d=3), ' linked: ', round(POSTERIOR_SCORE_LINKED, d=3), ' 12: ', round(POSTERIOR_SCORE_12, d=3), ' 21: ', round(POSTERIOR_SCORE_21, d=3),'\nrun',PTY_RUN))]]	
	rpw2		<- merge(rpw, unique(subset(rps, select=c(ID1,ID2,PTY_RUN))), by=c('ID1', 'ID2','PTY_RUN'))
	plot.file	<- paste0(outfile.base,'_phyloscans_lklpairs.pdf')
	phsc.plot.windowscan.for.pairs(rpw2, plot.file, plot.w=10, plot.h=40, id.cols=c('ID1','ID2'), ylim=NULL, cols.typet=NULL)
	
	rps		<- subset(rtp, grepl('ambiguous', SELECT), select=c(ID1, ID2, PTY_RUN, SELECT, POSTERIOR_SCORE_UNLINKED, POSTERIOR_SCORE_LINKED, POSTERIOR_SCORE_12, POSTERIOR_SCORE_21))
	setkey(rps, SELECT, POSTERIOR_SCORE_LINKED)
	write.csv(rps, file=paste0(outfile.base,'_summary_ambiguouspairs.csv'))
	rps[, DUMMY:=seq_len(nrow(rps))]
	rps[, LABEL:=rps[, factor(DUMMY, levels=DUMMY, labels=paste0('id1 ',ID1,' id2 ', ID2,'\n', SELECT,'\nunlinked: ',round(POSTERIOR_SCORE_UNLINKED, d=3), ' linked: ', round(POSTERIOR_SCORE_LINKED, d=3), ' 12: ', round(POSTERIOR_SCORE_12, d=3), ' 21: ', round(POSTERIOR_SCORE_21, d=3),'\nrun',PTY_RUN))]]	
	rpw2		<- merge(rpw, unique(subset(rps, select=c(ID1,ID2,PTY_RUN))), by=c('ID1', 'ID2','PTY_RUN'))
	plot.file	<- paste0(outfile.base,'_phyloscans_ambiguouspairs.pdf')
	phsc.plot.windowscan.for.pairs(rpw2, plot.file, plot.w=10, plot.h=350, id.cols=c('ID1','ID2'), ylim=NULL, cols.typet=NULL)
	
	
	if(0)
	{
		require(colorspace)
		#for(ii in seq_len(nrow(rtpdm))[-1])
		for(ii in 111:252)
		{		
			if(rtpdm[ii, PTY_RUN]!=1)
			{
				indir		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun'
				# load dfr and phs
				load( file.path(indir, paste0('ptyr',rtpdm[ii,PTY_RUN],'_trees.rda')) )
				# setup plotting
				ids			<- c(rtpdm[ii, MALE_RID],rtpdm[ii, FEMALE_RID])
				dfs			<- subset(dfr, select=c(W_FROM, W_TO, IDX))
				dfs[, MALE_RID:=ids[1]]
				dfs[, FEMALE_RID:=ids[2]]
				dfs			<- merge(dfs, subset(rpw, GROUP=='TYPE_RAW'), by=c('MALE_RID','FEMALE_RID','W_FROM','W_TO'), all.x=1)	
				dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', PTY_RUN, '\nwindow ', W_FROM,'-', W_TO,'\npmf',PATHS_MF,' pfm',PATHS_FM, ' ',ADJACENT,' ',CONTIGUOUS,' ',TYPE, '\n', round(PATRISTIC_DISTANCE, d=5), sep='')]]
				plot.file	<- paste0(outfile.base, 'trees/todi_pairs_170516_run_', rtpdm[ii, PTY_RUN],'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
				invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10, tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'))			
			}				
		}
	}				
}

Munich.preprocess.phyloscanneroutput.based.on.exclude.180924<- function()
{
	require(data.table)	
	require(igraph)
	require(sna)
	
	indir	<- '/Users/Oliver/Box Sync/OR_Work/2018/2018_MunichCluster/180924_phscoutput'
	outfile	<- '/Users/Oliver/Box Sync/OR_Work/2018/2018_MunichCluster/180924_analysis/180924_MunichCluster_w250_cl25_d50_min30_cut60.rda'
	
	conf.cut	<- 0.6
	neff.cut	<- 3
	
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
				#F<- '/Users/Oliver/Box Sync/OR_Work/2018/2018_MunichCluster/180924_phscoutput/ptyr1_pairwise_relationships.rda'
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
	#
	#	re-arrange to ID1<ID2
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
	tmp2	<- c('ID1','ID2','PTY_RUN')
	tmp		<- merge(tmp, subset(rplkl, GROUP==scores.group), by=tmp2)
	tmp		<- merge(tmp, tmp[, list(TYPE=TYPE, POSTERIOR_SCORE=(POSTERIOR_ALPHA-1)/sum(POSTERIOR_ALPHA-1) ), by=c('ID1','ID2','PTY_RUN')], by=c('ID1','ID2','PTY_RUN','TYPE'))
	rtn		<- rbind(tmp, rtn)
	#	generate maximum branch transmission network
	rtn		<- subset(rtn, GROUP==scores.group)
	rtnn	<- phsc.get.maximum.probability.transmission.network(rtn, verbose=0)	
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
	#rtnn	<- merge(rtnn, unique(subset(rtp.todi2, select=c(ID1,ID2,ID1_SEX,ID2_SEX))), by=c('ID1','ID2'))
	#
	save(rp, rd, rh, ra, rs, rtp.todi2, rplkl, rpw, rtn, rtnn, file=gsub('\\.rda','_networksallpairs.rda',outfile))
}

Munich.preprocess.phyloscanneroutput.based.on.strongsupport.180924<- function()
{
	require(data.table)	
	require(igraph)
	require(sna)
	
	indir	<- '/Users/Oliver/Box Sync/OR_Work/2018/2018_MunichCluster/180924_phscoutput'
	outfile	<- '/Users/Oliver/Box Sync/OR_Work/2018/2018_MunichCluster/180924_analysis/180924_MunichCluster_w250_cl25_d50_min30_cut60_strongsupport.rda'
	
	conf.cut	<- 0.6
	neff.cut	<- 3
	
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
				#F<- '/Users/Oliver/Box Sync/OR_Work/2018/2018_MunichCluster/180924_phscoutput/ptyr1_pairwise_relationships.rda'
				#cat(PTY_RUN,'\n')
				load(F)
				#	all pairs that are not decisively unlinked based on "close.group"
				rtp		<- subset(rplkl, 	GROUP==close.group & 
								TYPE==close.type.yes &
								((POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE))>=conf.cut,
						c(ID1, ID2))
				ans		<- merge(rtp, subset(rplkl, GROUP==close.group & TYPE==close.type.yes), by=c('ID1','ID2'), all.x=1)
				#	all pairs that are not decisively unlinked based on ML likely transmission pairs by distance + topology
				rtp		<- subset(rplkl, 	GROUP==linked.group & 
								TYPE==linked.type.yes &
								((POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-N_TYPE))>=conf.cut,
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
	#
	#	re-arrange to ID1<ID2
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
	tmp2	<- c('ID1','ID2','PTY_RUN')
	tmp		<- merge(tmp, subset(rplkl, GROUP==scores.group), by=tmp2)
	tmp		<- merge(tmp, tmp[, list(TYPE=TYPE, POSTERIOR_SCORE=(POSTERIOR_ALPHA-1)/sum(POSTERIOR_ALPHA-1) ), by=c('ID1','ID2','PTY_RUN')], by=c('ID1','ID2','PTY_RUN','TYPE'))
	rtn		<- rbind(tmp, rtn)
	#	generate maximum branch transmission network
	rtn		<- subset(rtn, GROUP==scores.group)
	rtnn	<- phsc.get.maximum.probability.transmission.network(rtn, verbose=0)	
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
	#rtnn	<- merge(rtnn, unique(subset(rtp.todi2, select=c(ID1,ID2,ID1_SEX,ID2_SEX))), by=c('ID1','ID2'))
	#
	save(rp, rd, rh, ra, rs, rtp.todi2, rplkl, rpw, rtn, rtnn, file=gsub('\\.rda','_networksallpairs.rda',outfile))
}