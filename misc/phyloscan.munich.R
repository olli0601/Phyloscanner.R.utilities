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
				p		<- phsc.plot.transmission.network(df, di, point.size=10, 
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
				p		<- phsc.plot.most.likely.transmission.chain(df, di, point.size=10, 
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
	phsc.plot.phyloscan(rpw2, plot.file, plot.w=10, plot.h=40, id.cols=c('ID1','ID2'), ylim=NULL, cols.typet=NULL)
	
	rps		<- subset(rtp, grepl('couple most likely a pair direction not resolved', SELECT), select=c(ID1, ID2, PTY_RUN, SELECT, POSTERIOR_SCORE_UNLINKED, POSTERIOR_SCORE_LINKED, POSTERIOR_SCORE_12, POSTERIOR_SCORE_21))
	setkey(rps, SELECT, POSTERIOR_SCORE_LINKED)
	write.csv(rps, file=paste0(outfile.base,'_summary_lklpairs_nodir.csv'))
	rps[, DUMMY:=seq_len(nrow(rps))]
	rps[, LABEL:=rps[, factor(DUMMY, levels=DUMMY, labels=paste0('id1 ',ID1,' id2 ', ID2,'\n', SELECT,'\nunlinked: ',round(POSTERIOR_SCORE_UNLINKED, d=3), ' linked: ', round(POSTERIOR_SCORE_LINKED, d=3), ' 12: ', round(POSTERIOR_SCORE_12, d=3), ' 21: ', round(POSTERIOR_SCORE_21, d=3),'\nrun',PTY_RUN))]]	
	rpw2		<- merge(rpw, unique(subset(rps, select=c(ID1,ID2,PTY_RUN))), by=c('ID1', 'ID2','PTY_RUN'))
	plot.file	<- paste0(outfile.base,'_phyloscans_lklpairs_nodir.pdf')
	phsc.plot.phyloscan(rpw2, plot.file, plot.w=10, plot.h=300, id.cols=c('ID1','ID2'), ylim=NULL, cols.typet=NULL)
	
	
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
	phsc.plot.phyloscan(rpw2, plot.file, plot.w=10, plot.h=40, id.cols=c('ID1','ID2'), ylim=NULL, cols.typet=NULL)
	
	rps		<- subset(rtp, grepl('ambiguous', SELECT), select=c(ID1, ID2, PTY_RUN, SELECT, POSTERIOR_SCORE_UNLINKED, POSTERIOR_SCORE_LINKED, POSTERIOR_SCORE_12, POSTERIOR_SCORE_21))
	setkey(rps, SELECT, POSTERIOR_SCORE_LINKED)
	write.csv(rps, file=paste0(outfile.base,'_summary_ambiguouspairs.csv'))
	rps[, DUMMY:=seq_len(nrow(rps))]
	rps[, LABEL:=rps[, factor(DUMMY, levels=DUMMY, labels=paste0('id1 ',ID1,' id2 ', ID2,'\n', SELECT,'\nunlinked: ',round(POSTERIOR_SCORE_UNLINKED, d=3), ' linked: ', round(POSTERIOR_SCORE_LINKED, d=3), ' 12: ', round(POSTERIOR_SCORE_12, d=3), ' 21: ', round(POSTERIOR_SCORE_21, d=3),'\nrun',PTY_RUN))]]	
	rpw2		<- merge(rpw, unique(subset(rps, select=c(ID1,ID2,PTY_RUN))), by=c('ID1', 'ID2','PTY_RUN'))
	plot.file	<- paste0(outfile.base,'_phyloscans_ambiguouspairs.pdf')
	phsc.plot.phyloscan(rpw2, plot.file, plot.w=10, plot.h=350, id.cols=c('ID1','ID2'), ylim=NULL, cols.typet=NULL)
	
	
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

Munich.phyloscanner.190715.make.trees<- function()
{
	require(data.table)
	require(Phyloscanner.R.utilities)
	
	#	set up working environment	
	if(1)
	{
		prog.pty <- '/rds/general/user/or105/home/phyloscanner/phyloscanner_make_trees.py'
		HOME <<- '/rds/general/user/or105/home/WORK/MUNICH'
		data.dir <- '/rds/general/user/or105/home/WORK/MUNICH/Data_190130'		
	}
	in.dir			<- file.path(HOME,'M190715_phsc_output')		
	work.dir		<- file.path(HOME,"M190715_phsc_work")
	out.dir			<- file.path(HOME,"M190715_phsc_output")	
	
	
	#
	#	generate trees	
	infiles	<- data.table(FI=list.files(in.dir, pattern='fasta$', full.names=TRUE, recursive=TRUE))
	infiles[, FO:= gsub('fasta$','tree',FI)]
	infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FI)))]
	infiles[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FI)))]
	#	check which (if any) trees have already been processed, and remove from TODO list
	tmp		<- data.table(FT=list.files(out.dir, pattern='^ptyr.*tree$', full.names=TRUE, recursive=TRUE))
	tmp[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FT)))]
	tmp[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FT)))]
	infiles	<- merge(infiles, tmp, by=c('PTY_RUN','W_FROM'), all.x=1)
	infiles	<- subset(infiles, is.na(FT))	
	setkey(infiles, PTY_RUN, W_FROM)			
	
	#
	#	define pbs header, max 10,000 trees (re-start a few times)
	#	create header first because RAXML call depends on single-core / multi-core
	infiles[, CASE_ID2:= seq_len(nrow(infiles))]
	infiles[, CASE_ID:= ceiling(CASE_ID2/1)]
	
		
	hpc.load			<- "module load intel-suite/2015.1 mpi raxml/8.2.9"	# make third party requirements available	 
	hpc.select			<- 1						# number of nodes
	hpc.nproc			<- 1						# number of processors on node
	hpc.walltime		<- 123						# walltime
	if(0)		
	{
		hpc.q			<- NA						# PBS queue
		hpc.mem			<- "2gb" 					# RAM
		raxml.pr		<- ifelse(hpc.nproc==1, 'raxmlHPC-SSE3', 'raxmlHPC-PTHREADS-SSE3')	#on older machines without AVX instructions
	}
	#		or run this block to submit a job array to Oliver's machines
	if(1)
	{
		hpc.q			<- "pqeelab"				# PBS queue
		hpc.mem			<- "6gb" 					# RAM
		raxml.pr		<- ifelse(hpc.nproc==1, 'raxmlHPC-AVX','raxmlHPC-PTHREADS-AVX')		#on newer machines with AVX instructions
	}
	#
	#	
	hpc.array			<- infiles[, max(CASE_ID)]	# number of runs for job array	
	#	define PBS header for job scheduler. this will depend on your job scheduler.
	pbshead		<- "#!/bin/sh"
	tmp			<- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
	pbshead		<- paste(pbshead, tmp, sep = "\n")
	tmp			<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
	pbshead 	<- paste(pbshead, tmp, sep = "\n")
	pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")	
	if(!is.na(hpc.array))
		pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')	
	if(!is.na(hpc.q)) 
		pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
	pbshead 	<- paste(pbshead, hpc.load, sep = "\n")	
	cat(pbshead)
	
	#
	#	create UNIX bash script to generate trees with RAxML
	#	
	raxml.args	<- ifelse(hpc.nproc==1, '-m GTRCAT --HKY85 -p 42 -o REF_B_K03455', paste0('-m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42 -o REF_B_K03455'))
	pty.c	<- infiles[, list(CMD=raxml.cmd(FI, outfile=FO, pr=raxml.pr, pr.args=raxml.args)), by=c('CASE_ID','CASE_ID2')]
	pty.c	<- pty.c[, list(CMD=paste(CMD, collapse='\n')), by='CASE_ID']
	pty.c[1,cat(CMD)]
	
	#
	#	make array job
	cmd		<- pty.c[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
	cmd		<- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]			
	cmd		<- paste(pbshead,cmd,sep='\n')	
	#	submit job
	outfile		<- gsub(':','',paste("trs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
	outfile		<- file.path(work.dir, outfile)
	cat(cmd, file=outfile)
	cmd 		<- paste("qsub", outfile)
	cat(cmd)
	cat(system(cmd, intern= TRUE))	
}


Munich.phyloscanner.190715.make.alignments<- function()
{
	#
	#	set up working environment
	require(Phyloscanner.R.utilities)
	if(1)
	{
		prog.pty <- '/rds/general/user/or105/home/phyloscanner/phyloscanner_make_trees.py'
		HOME <<- '/rds/general/user/or105/home/WORK/MUNICH'
		data.dir <- '/rds/general/user/or105/home/WORK/MUNICH/Data_190130'
	}
	in.dir			<- file.path(HOME,'M190715_phsc_input')		
	work.dir		<- file.path(HOME,"M190715_phsc_work")
	out.dir			<- file.path(HOME,"M190715_phsc_output")	
	dir.create(in.dir)
	dir.create(work.dir)
	dir.create(out.dir)
	
	#
	#	load runs
	if(1)
	{
		load(file.path(in.dir, "MunichCluster_190715.rda"))
		# loads pty.runs
	}	
	#
	#	search for bam files and references and merge with runs		
	tmp			<- phsc.find.bam.and.references(data.dir, 
						regex.person='^([A-Za-z0-9]+-[0-9]+)_.*$',
						regex.bam="^(.*)\\.bam$",
						regex.ref="^(.*)\\.fasta$")	
	setnames(tmp, c('IND','SAMPLE'), c('UNIT_ID','SAMPLE_ID'))
	tmp2 <- which(!grepl('Control',tmp$UNIT_ID))
	set(tmp, tmp2, 'UNIT_ID', tmp[tmp2, paste0('MC-',UNIT_ID)])
	set(tmp, NULL, 'UNIT_ID', tmp[, gsub('Control','CNTRL-', UNIT_ID)])	
	pty.runs	<- merge(pty.runs, tmp, by=c('UNIT_ID','SAMPLE_ID'))
	
	
	#
	#	define phyloscanner input args to generate read alignments 
	#	for each window and each run
	ptyi		<- seq(2000,5500,25)		
	pty.c		<- lapply(seq_along(ptyi), function(i)
			{
				pty.args			<- list(	prog.pty=prog.pty, 
						prog.mafft='mafft', 						 
						data.dir=data.dir, 
						work.dir=work.dir, 
						out.dir=out.dir, 
						alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_B.fasta"),
						alignments.root='REF_B_K03455', 
						alignments.pairwise.to='REF_B_K03455',
						window.automatic= '', 
						merge.threshold=1, 
						min.read.count=1, 
						quality.trim.ends=23, 
						min.internal.quality=23, 
						merge.paired.reads=TRUE, 
						no.trees=TRUE, 
						dont.check.duplicates=FALSE,
						dont.check.recombination=TRUE,
						num.bootstraps=1,
						all.bootstrap.trees=TRUE,
						strip.max.len=350, 
						min.ureads.individual=NA, 
						win=c(ptyi[i],ptyi[i]+250,25,250),				 				
						keep.overhangs=FALSE,
						mem.save=0,
						verbose=TRUE,					
						select=NA	
				)											
				pty.c <- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)
				pty.c[, W_FROM:= ptyi[i]]
				pty.c[, PTY_RUN:= as.integer(sub('.*ptyr([0-9])_.*','\\1',CMD))]
				pty.c
			})
	pty.c	<- do.call('rbind', pty.c)	
	tmp		<- data.table(FO=list.files(out.dir, pattern='ptyr.*fasta$', recursive=TRUE, full.names=TRUE))
	tmp[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FO)))]
	tmp[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FO)))]
	pty.c	<- merge(pty.c, tmp, by=c('PTY_RUN','W_FROM'), all.x=1)		
	pty.c	<- subset(pty.c, is.na(FO))		
	#pty.c	<- subset(pty.c, W_FROM==2350)
	#print(pty.c, n=1e3)

	#	define PBS variables
	hpc.load			<- "module load intel-suite/2015.1 mpi raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"	# make third party requirements available	 
	hpc.select			<- 1						# number of nodes
	hpc.nproc			<- 1						# number of processors on node
	hpc.walltime		<- 23						# walltime
	hpc.q				<- "pqeelab"				# PBS queue
	hpc.mem				<- "6gb" 					# RAM	
	hpc.array			<- pty.c[, max(CASE_ID)]	# number of runs for job array	
	#	define PBS header for job scheduler. this will depend on your job scheduler.
	pbshead		<- "#!/bin/sh"
	tmp			<- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
	pbshead		<- paste(pbshead, tmp, sep = "\n")
	tmp			<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
	pbshead 	<- paste(pbshead, tmp, sep = "\n")
	pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")	
	if(!is.na(hpc.array))
		pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')	
	if(!is.na(hpc.q)) 
		pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
	pbshead 	<- paste(pbshead, hpc.load, sep = "\n")	
	cat(pbshead)
	
	#	create PBS job array
	cmd		<- pty.c[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
	cmd		<- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]			
	cmd		<- paste(pbshead,cmd,sep='\n')	
	#	submit job
	outfile		<- gsub(':','',paste("readali",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
	outfile		<- file.path(work.dir, outfile)
	cat(cmd, file=outfile)
	cmd 		<- paste("qsub", outfile)
	cat(cmd)
	cat(system(cmd, intern= TRUE))		
}


Munich.phyloscanner.190715.make.ptyruns<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	
	
	# set up pty.runs file -- version 5
	# with cols	'SAMPLE_ID','RENAME_ID','UNIT_ID'
	# 39 individuals is too many, set up multiple analyses
	if(1)	
	{
		pty.runs		<- data.table(SAMPLE_ID=list.files('~/duke/2018_MunichCluster/Data_190130', pattern='bam$'))
		set(pty.runs, NULL, 'SAMPLE_ID', pty.runs[, gsub('\\.bam','',SAMPLE_ID)])
		#	manually define Doppelmeldung
		#	16-02593 (don t have)
		#	16-01790 (don t have)
		#	we have both 18-01721, 17-02604 but only PRRT from 17-02604. So just use 18-01721.  
		pty.runs <- subset(pty.runs, !grepl('17-02604',SAMPLE_ID))
		tmp	<- pty.runs[, which(!grepl('Control',SAMPLE_ID))]
		set(pty.runs, tmp, 'RENAME_ID', pty.runs[tmp,paste0('MC-',SAMPLE_ID)])
		tmp	<- pty.runs[, which(grepl('Control',SAMPLE_ID))]
		set(pty.runs, tmp, 'RENAME_ID', pty.runs[tmp,gsub('Control','CNTRL-',SAMPLE_ID)])
		pty.runs[, UNIT_ID:= gsub('_INT|_PRRT','',RENAME_ID)]
		#	define runs for individuals in periphery 
		load('~/Box Sync/OR_Work/2018/2018_MunichCluster/180924_analysis/180924_MunichCluster_w250_cl25_d50_min30_cut60_strongsupport_networksallpairs.rda')
		
		conf.cut <- 0.5
		pinds <- c('MC-16-03807', 'MC-16-03259', 'MC-16-03499', 'MC-15-01402', 'MC-16-04475', 'MC-16-03644', 'MC-16-03516', 'MC-15-04719')
		pty.runs2 <- lapply(seq_along(pinds), function(i)
				{
					#pind<- 'MC-16-03807'	
					pind <- pinds[i]
					tmp <- subset(rplkl, (ID1==pind | ID2==pind) & GROUP=="TYPE_CHAIN_TODI" & TYPE=='chain' & KEFF/NEFF>conf.cut)
					data.table( PTY_RUN=i, UNIT_ID= unique(c(tmp$ID1, tmp$ID2)) )			
				})
		pty.runs2<- do.call('rbind',pty.runs2)		
		tmp <- setdiff( subset(pty.runs, !grepl('CNTRL',UNIT_ID))[, unique(UNIT_ID)], pinds )
		tmp <- data.table( PTY_RUN= max(pty.runs2$PTY_RUN)+1L,  UNIT_ID=tmp)
		pty.runs2<- rbind(pty.runs2, tmp)
		#	run #4 is smallest, add all controls
		tmp <- data.table(PTY_RUN= 4, UNIT_ID=subset(pty.runs, grepl('CNTRL',UNIT_ID))[, unique(UNIT_ID)])
		pty.runs2<- rbind(pty.runs2,tmp)		
		#	for all other runs except #4, add two standard controls 
		tmp <- as.data.table(expand.grid(PTY_RUN= setdiff(1:max(pty.runs2$PTY_RUN),4), UNIT_ID=c('CNTRL-16-01405','CNTRL-16-01277')))
		pty.runs2<- rbind(pty.runs2,tmp)
		
		#	check numbers
		pty.runs2[, list(N=length(unique(UNIT_ID))), by='PTY_RUN']
				
		
		#	make overall run file
		pty.runs <- merge(pty.runs2, pty.runs, by='UNIT_ID', allow.cartesian=TRUE)
		outfile	<- '~/duke/2018_MunichCluster/Data_190130/MunichCluster_190715.rda'
		save(pty.runs, file=outfile)
	}
}

Munich.phyloscanner.180605<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	
	# set up pty.runs file
	# with cols	'SAMPLE_ID','RENAME_ID','UNIT_ID'
	if(0)	
	{
		pty.runs		<- data.table(SAMPLE_ID=list.files('~/duke/2018_MunichCluster/Data', pattern='bam$'))
		set(pty.runs, NULL, 'SAMPLE_ID', pty.runs[, gsub('\\.bam','',SAMPLE_ID)])
		pty.runs[, RENAME_ID:= paste0('MC-',SAMPLE_ID)]
		pty.runs[, UNIT_ID:= paste0('MC-',SAMPLE_ID)]
		pty.runs[, PTY_RUN:=1]		
		outfile	<- '~/duke/2018_MunichCluster/Data/MunichCluster_180605.rda'
		save(pty.runs, file=outfile)
	}
	# set up pty.runs file -- version 2
	# with cols	'SAMPLE_ID','RENAME_ID','UNIT_ID'
	if(0)	
	{
		pty.runs		<- data.table(SAMPLE_ID=list.files('~/duke/2018_MunichCluster/Data_180607', pattern='bam$'))
		set(pty.runs, NULL, 'SAMPLE_ID', pty.runs[, gsub('\\.bam','',SAMPLE_ID)])
		pty.runs[, RENAME_ID:= paste0('MC-',SAMPLE_ID)]
		pty.runs[, UNIT_ID:= paste0('MC-',gsub('_INT|_PRRT','',SAMPLE_ID))]		
		pty.runs[, PTY_RUN:=1]		
		outfile	<- '~/duke/2018_MunichCluster/Data/MunichCluster_180607.rda'
		save(pty.runs, file=outfile)
	}
	# set up pty.runs file -- version 3
	# with cols	'SAMPLE_ID','RENAME_ID','UNIT_ID'
	if(0)	
	{
		pty.runs		<- data.table(SAMPLE_ID=list.files('~/duke/2018_MunichCluster/Data_180618', pattern='bam$'))
		set(pty.runs, NULL, 'SAMPLE_ID', pty.runs[, gsub('\\.bam','',SAMPLE_ID)])
		pty.runs[, RENAME_ID:= paste0('MC-',SAMPLE_ID)]
		pty.runs[, UNIT_ID:= paste0('MC-',gsub('_INT|_PRRT','',SAMPLE_ID))]		
		pty.runs[, PTY_RUN:=1]		
		outfile	<- '~/duke/2018_MunichCluster/Data_180618/MunichCluster_180618.rda'
		save(pty.runs, file=outfile)
	}
	# set up pty.runs file -- version 4
	# with cols	'SAMPLE_ID','RENAME_ID','UNIT_ID'
	# 39 individuals is too many
	if(0)	
	{
		pty.runs		<- data.table(SAMPLE_ID=list.files('~/duke/2018_MunichCluster/Data_180815', pattern='bam$'))
		set(pty.runs, NULL, 'SAMPLE_ID', pty.runs[, gsub('\\.bam','',SAMPLE_ID)])
		tmp	<- pty.runs[, which(!grepl('Control',SAMPLE_ID))]
		set(pty.runs, tmp, 'RENAME_ID', pty.runs[tmp,paste0('MC-',SAMPLE_ID)])
		tmp	<- pty.runs[, which(grepl('Control',SAMPLE_ID))]
		set(pty.runs, tmp, 'RENAME_ID', pty.runs[tmp,gsub('Control','CNTRL-',SAMPLE_ID)])
		pty.runs[, UNIT_ID:= gsub('_INT|_PRRT','',RENAME_ID)]	
		pty.runs[, PTY_RUN:=1]		
		outfile	<- '~/duke/2018_MunichCluster/Data_180815/MunichCluster_180815.rda'
		save(pty.runs, file=outfile)
	}
	#
	#	INPUT ARGS
	if(1)
	{	
		#HOME				<<- '/Users/Oliver/duke/2018_MunichCluster'
		HOME				<<- '/work/or105/MUNICH'
		in.dir				<- file.path(HOME,"MunichCluster_180815_in")
		work.dir			<- file.path(HOME,"MunichCluster_180815_work")
		out.dir				<- file.path(HOME,"MunichCluster_180815_out")		
		load( file.path(in.dir, 'MunichCluster_180815.rda') )
		print(pty.runs)
		hpc.load			<- "module load intel-suite/2015.1 mpi R/3.3.3 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"
		hpc.nproc			<- 1
		#prog.pty			<- '/Users/Oliver/git/phylotypes/phyloscanner_make_trees.py'
		#pty.data.dir		<- '~/duke/2018_MunichCluster/Data'
		#prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-SSE3 -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-SSE3 -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))
		prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner_make_trees.py'
		#prog.pty			<- '/work/or105/libs/phyloscanner/phyloscanner_make_trees.py'
		#pty.data.dir		<- '/work/or105/MUNICH/data'
		pty.data.dir		<- '/work/or105/MUNICH/Data_180815'
		prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-AVX -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-AVX -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))
		pty.select			<- 1		
	}	
	
	#
	# generate read alignments
	if(0)
	{		
		#ptyi		<- seq(800,9150,25)
		ptyi		<- seq(2000,5500,25)		
		#ptyi		<- seq(2475,2525,25)
		pty.c		<- lapply(seq_along(ptyi), function(i)
				{
					pty.args			<- list(	prog.pty=prog.pty, 
							prog.mafft='mafft', 
							prog.raxml=prog.raxml, 
							data.dir=pty.data.dir, 
							work.dir=work.dir, 
							out.dir=out.dir, 
							alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_B.fasta"),
							alignments.root='REF_B_K03455', 
							alignments.pairwise.to='REF_B_K03455',
							window.automatic= '', 
							merge.threshold=1, 
							min.read.count=1, 
							quality.trim.ends=23, 
							min.internal.quality=23, 
							merge.paired.reads=TRUE, 
							no.trees=TRUE, 
							dont.check.duplicates=FALSE,
							dont.check.recombination=TRUE,
							num.bootstraps=1,
							all.bootstrap.trees=TRUE,
							strip.max.len=350, 
							min.ureads.individual=NA, 
							win=c(ptyi[i],ptyi[i]+250,25,250),				 				
							keep.overhangs=FALSE,
							mem.save=0,
							verbose=TRUE,					
							select=pty.select	#of 240
					)											
					pty.c				<- phsc.cmd.phyloscanner.multi(pty.runs, pty.args,regex.ref='\\.fasta$', postfix.sample.id='\\.bam|\\.fasta')
					#cat(pty.c$CMD)
					#stop()
					pty.c[, W_FROM:= ptyi[i]]
					pty.c
				})
		pty.c	<- do.call('rbind', pty.c)	
		tmp		<- data.table(FO=list.files(out.dir, pattern='ptyr.*fasta$', recursive=TRUE, full.names=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FO)))]
		tmp[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FO)))]
		pty.c	<- merge(pty.c, tmp, by=c('PTY_RUN','W_FROM'), all.x=1)		
		pty.c	<- subset(pty.c, is.na(FO))
		#pty.c	<- subset(pty.c, W_FROM==2350)
		print(pty.c, n=1e3)
		#	make array job
		pty.c[, CASE_ID:= seq_len(nrow(pty.c))]
		tmp		<- pty.c[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
		tmp		<- tmp[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]		
		cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=18, hpc.q="pqeelab", hpc.mem="6gb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load, hpc.array=pty.c[, max(CASE_ID)])
		cmd		<- paste(cmd,tmp,sep='\n')
		outfile	<- paste("mc",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(work.dir, outfile, cmd)
		
		#print(pty.c)
		#stop()		 
		#invisible(pty.c[,	{
		#					cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=18, hpc.q="pqeelab", hpc.mem="6gb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
		#cmd	<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=71, hpc.q=NA, hpc.mem="1850mb",  hpc.nproc=1, hpc.load=hpc.load)
		#cmd	<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=71, hpc.q=NA, hpc.mem="63800mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
		#					cmd		<- paste(cmd,CMD,sep='\n')
		#					cat(cmd)					
		#					outfile	<- paste("mc",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		#					cmd.hpccaller(work.dir, outfile, cmd)
		#stop()
		#				}, by=c('PTY_RUN','W_FROM')])
		quit('no')
	}
	#
	# generate trees
	if(0)
	{
		#HOME		<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'	
		hpc.load	<- "module load intel-suite/2015.1 mpi raxml/8.2.9"
		hpc.select	<- 1; hpc.nproc<- 1; 	hpc.walltime<- 71; hpc.mem<-"6gb"; hpc.q<- "pqeelab"
		#hpc.select	<- 1; hpc.nproc<- 16; 	hpc.walltime<- 23; hpc.mem<-"124gb"; hpc.q<- NA
		#raxml.pr	<- ifelse(hpc.nproc==1, 'raxmlHPC-SSE3', 'raxmlHPC-PTHREADS-SSE3')	
		raxml.pr	<- ifelse(hpc.nproc==1, 'raxmlHPC-AVX','raxmlHPC-PTHREADS-AVX')
		raxml.args	<- ifelse(hpc.nproc==1, '-m GTRCAT --HKY85 -p 42 -o REF_B_K03455', paste0('-m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42 -o REF_B_K03455'))
		in.dir		<- file.path(HOME,"MunichCluster_180815_out")
		out.dir		<- in.dir
		work.dir	<- file.path(HOME,"MunichCluster_180815_work")
		
		infiles	<- data.table(FI=list.files(in.dir, pattern='fasta$', full.names=TRUE, recursive=TRUE))
		infiles[, FO:= gsub('fasta$','tree',FI)]
		infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FI)))]
		infiles[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FI)))]		
		tmp		<- data.table(FT=list.files(out.dir, pattern='^ptyr.*tree$', full.names=TRUE, recursive=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FT)))]
		tmp[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FT)))]
		infiles	<- merge(infiles, tmp, by=c('PTY_RUN','W_FROM'), all.x=1)
		infiles	<- subset(infiles, is.na(FT))	
		setkey(infiles, PTY_RUN, W_FROM)			
		#infiles	<- subset(infiles, PTY_RUN>=240 & PTY_RUN<300)
		#infiles	<- subset(infiles, PTY_RUN>=186 & PTY_RUN<239)
		#print(infiles)	 		
		
		#	make raxml run
		pty.c	<- infiles[, list(CMD=cmd.raxml(FI, outfile=FO, pr=raxml.pr, pr.args=raxml.args)), by=c('PTY_RUN','W_FROM')]		
		print(pty.c, n=1e3)
		
		#	make array job
		pty.c[, CASE_ID:= seq_len(nrow(pty.c))]
		tmp		<- pty.c[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
		tmp		<- tmp[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]		
		cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=hpc.load, hpc.array=pty.c[, max(CASE_ID)])
		cmd		<- paste(cmd,tmp,sep='\n')
		outfile	<- paste("mct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(work.dir, outfile, cmd)				
	}
	#
	#	combine all the data	
	if(0)
	{		
		indirs 	<- '/Users/Oliver/duke/tmp/ptyr143_trees'
		indirs	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/BEEHIVE_67_180302_out'
		indirs	<- '/work/or105/MUNICH/MunichCluster_180815_out'
		#
		indirs	<- list.files(indirs, pattern='^ptyr[0-9]+_trees$', full.names=TRUE)
		allwin	<- data.table(W_FROM=seq(2050,5075,25))
		#allwin	<- data.table(W_FROM=seq(800,9150,25))		
		#allwin	<- data.table(W_FROM=seq(800,9050,125))
		#indirs	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo/ptyr97_trees'
		for(i in seq_along(indirs))
		{
			indir	<- indirs[i]
			pty.run	<- as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(indir)))
			#	check if we have all fasta and tree files
			infiles	<- data.table(F=list.files(indir,pattern='ptyr.*fasta$',full.names=TRUE))
			infiles[, W_FROM:= as.integer(gsub('.*_InWindow_([0-9]+)_.*','\\1',basename(F)))]
			infiles	<- merge(allwin, infiles, by='W_FROM', all.x=1)			 					
			missfs	<- subset(infiles, is.na(F))[, W_FROM]
			if(length(missfs))
				cat('\nIn',indir,'Found missing fasta files for',paste(missfs,collapse=', '))
			infiles	<- data.table(F=list.files(indir,pattern='ptyr.*tree$',full.names=TRUE))
			infiles[, W_FROM:= as.integer(gsub('.*_InWindow_([0-9]+)_.*','\\1',basename(F)))]
			infiles	<- merge(allwin, infiles, by='W_FROM', all.x=1)
			misstrs	<- subset(infiles, is.na(F))[, W_FROM]
			if(length(misstrs))
				cat('\nIn',indir,'Found missing tree files for',paste(misstrs,collapse=', '))
			zipit	<- 0
			if(!length(missfs) & !length(misstrs))
			{
				cat('\nIn',indir,'Found all fasta and tree files')
				zipit	<- 1
			}				
			if(!length(setdiff(misstrs,missfs)))
			{ 
				cat('\nIn',indir,'Found all tree files for which there is a fasta file')
				zipit	<- 1
			}	
			#zipit	<- 1
			#
			if(zipit)
			{			
				cat('\nProcess',indir)
				#	first combine all zip files into ptyrXXX_otherstuff.zip
				infiles	<- data.table(F=list.files(indir,pattern='zip$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				stopifnot(!nrow(subset(infiles, is.na(PTY_RUN))))		
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_otherstuff.zip'))]
				cat('\nZip to file', tmp,'...\n')
				suppressWarnings(invisible( infiles[, list(RTN= unzip(F, overwrite=FALSE, exdir=file.path(indir,'tmp42'))), by='F'] ))
				invisible( infiles[, list(RTN= file.remove(F)), by='F'] )
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	now zip fasta files
				infiles	<- data.table(F=list.files(indir,pattern='ptyr.*fasta$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				invisible( infiles[, file.rename(F, file.path(indir,'tmp42',basename(F))), by='F'] )
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_trees_fasta.zip'))]
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	now zip tree files
				infiles	<- data.table(F=list.files(indir,pattern='ptyr.*tree$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				invisible( infiles[, file.rename(F, file.path(indir,'tmp42',basename(F))), by='F'] )
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_trees_newick.zip'))]		
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	remove tmp dir
				invisible( file.remove( file.path(indir,'tmp42') ) )
				#	move one level down
				infiles	<- data.table(F=list.files(indir, full.names=TRUE))
				invisible( infiles[, file.rename(F, file.path(dirname(indir),basename(F))), by='F'] )
				cat('\nDone',indir)
			}
			#if(!length(misstrs))
			if(zipit)
				invisible(unlink(indir, recursive=TRUE))
			#	expand again if asked to
			#if(length(misstrs))
			if(0)
			{
				cat('\nExtract',file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_fasta.zip')))
				unzip(file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_fasta.zip')), junkpaths=TRUE, exdir=indir)
				cat('\nExtract',file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_newick.zip')))
				unzip(file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_newick.zip')), junkpaths=TRUE, exdir=indir)
			}
		}					
	}
	#
	#	process all trees in one go
	if(0)	
	{	
		hpc.select	<- 1; hpc.nproc<- 1; 	hpc.walltime<- 471; hpc.mem<-"6gb"; hpc.q<- "pqeelab"
		#HOME				<<- '/Users/Oliver/duke/2018_MunichCluster'
		HOME				<<- '/work/or105/MUNICH'
		in.dir				<- file.path(HOME,"MunichCluster_180815_out")
		out.dir				<- in.dir
		work.dir			<- file.path(HOME,"MunichCluster_180815_work") 			
		dir.create(out.dir, showWarnings=FALSE)
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft=NA, 
				prog.raxml=NA, 
				data.dir=NA, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_B.fasta"),
				alignments.root='REF_B_K03455', 
				alignments.pairwise.to='REF_B_K03455',				
				bl.normalising.reference.file=system.file(package="Phyloscanner.R.utilities", "data", "hiv.hxb2.norm.constants.rda"),
				bl.normalising.reference.var='MEDIAN_PWD',														
				window.automatic= '', 
				merge.threshold=1, 
				min.read.count=1, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				merge.paired.reads=TRUE, 
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				dont.check.recombination=TRUE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				win=c(2000,5750,25,250), 				
				keep.overhangs=FALSE,
				use.blacklisters=c('ParsimonyBasedBlacklister'), #,'DownsampleReads'),
				tip.regex='^(.*)_[A-Z]+_read_([0-9]+)_count_([0-9]+)$',
				roguesubtree.kParam=20,
				roguesubtree.prop.threshold=0,
				roguesubtree.read.threshold=20,
				#dwns.maxReadsPerPatient=1000,	
				multifurcation.threshold=1e-5,
				split.rule='s',
				split.kParam=20,
				split.proximityThreshold=0,
				split.readCountsMatterOnZeroBranches=TRUE,
				split.pruneBlacklist=FALSE,
				trms.allowMultiTrans=TRUE,
				pw.trmw.min.reads=30,									
				pw.trmw.min.tips=1,
				pw.trmw.close.brl=0.025,
				pw.trmw.distant.brl=0.05,
				pw.prior.keff=2,
				pw.prior.neff=3,
				pw.prior.keff.dir=2,
				pw.prior.neff.dir=3,				
				pw.prior.calibrated.prob=0.6,
				mem.save=0,
				verbose=TRUE,				
				select=NA 
		)	
		save(pty.args, file=file.path(out.dir, 'pty.args.rda'))
		
		pty.c	<- data.table(FILE_BAM=list.files(in.dir, pattern='_patients.txt', full.names=TRUE))
		pty.c[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_patients.txt','',basename(FILE_BAM))))]				
		tmp		<- data.table(FILE_TRMW=list.files(out.dir, pattern='_pairwise_relationships.rda', full.names=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_pairwise_relationships.rda','',basename(FILE_TRMW))))]
		pty.c	<- merge(pty.c, tmp, by='PTY_RUN', all.x=1)
		pty.c	<- subset(pty.c, is.na(FILE_TRMW))
		setkey(pty.c, PTY_RUN)		
		pty.c	<- pty.c[, { 
					#FILE_BAM<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr1_bam.txt'
					#FILE_BAM<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr1_bam.txt'
					#cat('\n',FILE_BAM)
					pty.args$process.window<- 2500
					prefix.infiles	<- gsub('patients.txt','',FILE_BAM)
					print(prefix.infiles)
					cmd				<- Phyloscanner.R.utilities:::phsc.cmd.phyloscanner.one.resume.onewindow(prefix.infiles, pty.args)
					cat(cmd)
					stop()
					list(CMD=cmd) 
				}, by='PTY_RUN']		
		pty.c[1,cat(CMD)]		
		#stop()		
		invisible(pty.c[,	{
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=hpc.load)							
							cmd			<- paste(cmd,'cd $TMPDIR',sep='\n')
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile		<- paste("scRAr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')														
							cmd.hpccaller(work.dir, outfile, cmd)														
						}, by='PTY_RUN'])		
		quit('no')		
	}
	#
	#	process windows one by one, part 1
	if(0)	
	{	
		#hpc.select	<- 1; hpc.nproc<- 1; 	hpc.walltime<- 471; hpc.mem<-"6gb"; hpc.q<- "pqeelab"
		hpc.select	<- 1; hpc.nproc<- 1; 	hpc.walltime<- 71; hpc.mem<-"6gb"; hpc.q<- "pqeelab"
		#HOME				<<- '/Users/Oliver/duke/2018_MunichCluster'
		HOME				<<- '/work/or105/MUNICH'
		in.dir				<- file.path(HOME,"MunichCluster_180815_out")
		out.dir				<- in.dir
		work.dir			<- file.path(HOME,"MunichCluster_180815_work") 			
		dir.create(out.dir, showWarnings=FALSE)
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft=NA, 
				prog.raxml=NA, 
				data.dir=NA, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_B.fasta"),
				alignments.root='REF_B_K03455', 
				alignments.pairwise.to='REF_B_K03455',				
				bl.normalising.reference.file=system.file(package="Phyloscanner.R.utilities", "data", "hiv.hxb2.norm.constants.rda"),
				bl.normalising.reference.var='MEDIAN_PWD',														
				window.automatic= '', 
				merge.threshold=1, 
				min.read.count=1, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				merge.paired.reads=TRUE, 
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				dont.check.recombination=TRUE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				win=c(2000,5750,25,250), 				
				keep.overhangs=FALSE,
				use.blacklisters=c('ParsimonyBasedBlacklister'), #,'DownsampleReads'),
				tip.regex='^(.*)_[A-Z]+_read_([0-9]+)_count_([0-9]+)$',
				roguesubtree.kParam=20,
				roguesubtree.prop.threshold=0,
				roguesubtree.read.threshold=20,
				#dwns.maxReadsPerPatient=1000,	
				multifurcation.threshold=1e-5,
				split.rule='s',
				split.kParam=20,
				split.proximityThreshold=0,
				split.readCountsMatterOnZeroBranches=TRUE,
				split.pruneBlacklist=FALSE,
				trms.allowMultiTrans=TRUE,
				pw.trmw.min.reads=30,									
				pw.trmw.min.tips=1,
				pw.trmw.close.brl=0.025,
				pw.trmw.distant.brl=0.05,
				pw.prior.keff=2,
				pw.prior.neff=3,
				pw.prior.keff.dir=2,
				pw.prior.neff.dir=3,				
				pw.prior.calibrated.prob=0.6,
				mem.save=0,
				verbose=TRUE,				
				select=NA 
		)	
		save(pty.args, file=file.path(out.dir, 'pty.args.rda'))
		
		pty.c	<- data.table(FILE_BAM=list.files(in.dir, pattern='_patients.txt', full.names=TRUE))
		pty.c[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_patients.txt','',basename(FILE_BAM))))]				
		tmp		<- data.table(FILE_TRMW=list.files(out.dir, pattern='_pairwise_relationships.rda', full.names=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_pairwise_relationships.rda','',basename(FILE_TRMW))))]
		pty.c	<- merge(pty.c, tmp, by='PTY_RUN', all.x=1)
		pty.c	<- subset(pty.c, is.na(FILE_TRMW))
		setkey(pty.c, PTY_RUN)
		pty.c[, DUMMY:=1]
		tmp		<- data.table(DUMMY=1, W_FROM=seq(2050,5075,25))
		pty.c	<- merge(pty.c, tmp, by='DUMMY')
		pty.c	<- pty.c[, { 
					pty.args$process.window	<- W_FROM
					prefix.infiles			<- gsub('patients.txt','',FILE_BAM)					
					cmd						<- Phyloscanner.R.utilities:::phsc.cmd.phyloscanner.one.resume.onewindow(prefix.infiles, pty.args)
					list(CMD=cmd) 
				}, by=c('PTY_RUN','W_FROM')]		
		pty.c[1,cat(CMD)]	
		
		#	make array job
		pty.c[, CASE_ID:= seq_len(nrow(pty.c))]
		tmp		<- pty.c[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
		tmp		<- tmp[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]		
		cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=hpc.load, hpc.array=pty.c[, max(CASE_ID)])
		cmd		<- paste(cmd,tmp,sep='\n')
		outfile	<- paste("mct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(work.dir, outfile, cmd)
		
		quit('no')		
	}
	#
	#	process windows one by one, part 2
	if(1)	
	{	
		
		hpc.select	<- 1; hpc.nproc<- 1; 	hpc.walltime<- 167; hpc.mem<-"6gb"; hpc.q<- "pqeelab"
		#HOME				<<- '/Users/Oliver/duke/2018_MunichCluster'
		HOME				<<- '/work/or105/MUNICH'
		in.dir				<- file.path(HOME,"MunichCluster_180815_out")
		out.dir				<- in.dir
		work.dir			<- file.path(HOME,"MunichCluster_180815_work") 			
		dir.create(out.dir, showWarnings=FALSE)
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft=NA, 
				prog.raxml=NA, 
				data.dir=NA, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_B.fasta"),
				alignments.root='REF_B_K03455', 
				alignments.pairwise.to='REF_B_K03455',				
				bl.normalising.reference.file=system.file(package="Phyloscanner.R.utilities", "data", "hiv.hxb2.norm.constants.rda"),
				bl.normalising.reference.var='MEDIAN_PWD',														
				window.automatic= '', 
				merge.threshold=1, 
				min.read.count=1, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				merge.paired.reads=TRUE, 
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				dont.check.recombination=TRUE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				win=c(2000,5750,25,250), 				
				keep.overhangs=FALSE,
				use.blacklisters=c('ParsimonyBasedBlacklister'), #,'DownsampleReads'),
				tip.regex='^(.*)_[A-Z]+_read_([0-9]+)_count_([0-9]+)$',
				roguesubtree.kParam=20,
				roguesubtree.prop.threshold=0,
				roguesubtree.read.threshold=20,
				#dwns.maxReadsPerPatient=1000,	
				multifurcation.threshold=1e-5,
				split.rule='s',
				split.kParam=20,
				split.proximityThreshold=0,
				split.readCountsMatterOnZeroBranches=TRUE,
				split.pruneBlacklist=FALSE,
				trms.allowMultiTrans=TRUE,
				pw.trmw.min.reads=30,									
				pw.trmw.min.tips=1,
				pw.trmw.close.brl=0.025,
				pw.trmw.distant.brl=0.05,
				pw.prior.keff=2,
				pw.prior.neff=3,
				pw.prior.keff.dir=2,
				pw.prior.neff.dir=3,				
				pw.prior.calibrated.prob=0.6,
				mem.save=0,
				verbose=TRUE,				
				select=NA 
		)	
		save(pty.args, file=file.path(out.dir, 'pty.args.rda'))
		
		pty.c	<- data.table(FILE_BAM=list.files(in.dir, pattern='_patients.txt', full.names=TRUE))
		pty.c[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_patients.txt','',basename(FILE_BAM))))]				
		tmp		<- data.table(FILE_TRMW=list.files(out.dir, pattern='_pairwise_relationships.rda', full.names=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_pairwise_relationships.rda','',basename(FILE_TRMW))))]
		pty.c	<- merge(pty.c, tmp, by='PTY_RUN', all.x=1)
		pty.c	<- subset(pty.c, is.na(FILE_TRMW))
		setkey(pty.c, PTY_RUN)
		pty.args$process.window				<- NULL
		pty.args$combine.processed.windows	<- 1
		pty.c	<- pty.c[, { 					
					prefix.infiles			<- gsub('patients.txt','',FILE_BAM)					
					cmd						<- Phyloscanner.R.utilities:::phsc.cmd.phyloscanner.one.resume.combinewindows(prefix.infiles, pty.args)
					list(CMD=cmd) 
				}, by=c('PTY_RUN')]		
		pty.c[1,cat(CMD)]	
		
		invisible(pty.c[,	{
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=hpc.load)							
							cmd			<- paste(cmd,'cd $TMPDIR',sep='\n')
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile		<- paste("mcc",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')														
							cmd.hpccaller(work.dir, outfile, cmd)														
						}, by='PTY_RUN'])		
		quit('no')				
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

Munich.readlengths.180605<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	infile.base	<- '~/Box Sync/OR_Work/2018/2018_MunichCluster/mergedfragmentlen__'
	infile.bams	<- paste0(infile.base,c('cov175.rda','cov200.rda','cov225.rda','cov250.rda','cov275.rda','cov300.rda','cov325.rda','cov350.rda','cov375.rda','lendist.rda'))		
	outfile.base<- '~/Box Sync/OR_Work/2018/2018_MunichCluster/MunichCluster_bamstats_'						
	
	for(i in seq_along(infile.bams))		
		load(infile.bams[i])	
	setnames(bam.cov175, 'FILE_ID', 'SID')
	bam.cov175	<- subset(bam.cov175, !is.na(COV))
	setnames(bam.cov200, 'FILE_ID', 'SID')
	bam.cov200	<- subset(bam.cov200, !is.na(COV))
	setnames(bam.cov225, 'FILE_ID', 'SID')
	bam.cov225	<- subset(bam.cov225, !is.na(COV))	
	setnames(bam.cov250, 'FILE_ID', 'SID')
	bam.cov250	<- subset(bam.cov250, !is.na(COV))
	setnames(bam.cov275, 'FILE_ID', 'SID')
	bam.cov275	<- subset(bam.cov275, !is.na(COV))
	setnames(bam.cov300, 'FILE_ID', 'SID')
	bam.cov300	<- subset(bam.cov300, !is.na(COV))
	setnames(bam.cov325, 'FILE_ID', 'SID')
	bam.cov325	<- subset(bam.cov325, !is.na(COV))	
	setnames(bam.cov350, 'FILE_ID', 'SID')
	bam.cov350	<- subset(bam.cov350, !is.na(COV))
	setnames(bam.cov375, 'FILE_ID', 'SID')
	bam.cov375	<- subset(bam.cov375, !is.na(COV))	
	setnames(bam.len, 'FILE_ID', 'SID')
	bam.len[, LONGER:= factor(CDF>1-1e-4, levels=c(TRUE,FALSE), labels=c('no','yes'))]
	
	
	ggplot(subset(bam.len,QU>=40 & QU<320), aes(y=1-CDF, x=factor(QU))) + 
			geom_boxplot() +
			scale_y_continuous(lab=scales:::percent, expand=c(0,0)) +
			theme_bw() + labs(y='proportion of reads longer than x in one individual\n(boxplot over all individuals)', x='\nx=length of quality-trimmed short reads\n(nt)', fill='sequence run') +
			theme(legend.position='bottom')
	ggsave(file=paste0(outfile.base,'bamlen_longer_than_xnt.pdf'), w=6, h=6)
	#	comparing to beehive, perhaps a good threshold is 250	
	
	
	#
	df			<- copy(bam.cov175)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=175L]
	ans			<- copy(tmp)		
	#
	df			<- copy(bam.cov200)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=200L]
	ans			<- rbind(ans, tmp)	
	#
	df			<- copy(bam.cov225)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=225L]
	ans			<- rbind(ans, tmp)	
	#
	df			<- copy(bam.cov250)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=250L]
	ans			<- rbind(ans, tmp)		
	#
	df			<- copy(bam.cov275)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=275L]
	ans			<- rbind(ans, tmp)	
	#
	df			<- copy(bam.cov300)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=300L]
	ans			<- rbind(ans, tmp)	
	#
	df			<- copy(bam.cov325)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=325L]
	ans			<- rbind(ans, tmp)	
	#
	df			<- copy(bam.cov350)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=350L]
	ans			<- rbind(ans, tmp)	
	#
	df			<- copy(bam.cov375)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{						
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=375L]
	ans			<- rbind(ans, tmp)	
	#	make table
	ans[, XCOV_LABEL:= paste0(round(XCOV_MEAN,d=0),'x ( ', round(XCOV_QL, d=1),'x - ',round(XCOV_QU, d=0),'x )')]
	ans[, HCOV_LABEL:= paste0(round(100*HCOV_MEAN,d=1),'% (', round(100*HCOV_QL, d=1),'% - ',round(100*HCOV_QU, d=1),'%)')]
	ans[, P:= paste0( round(100*N / ans[COV_MIN==1 & READL==1, N], d=1), '%')]
	set(ans, NULL, 'COV_MIN', ans[, factor(COV_MIN, levels=c(1,10,20,30,50,100), labels=c('1X','10X','20X','30X','50X','100X'))])
	ans			<- subset(ans, select=c(COV_MIN, READL, N, P, HCOV_LABEL, HCOV_MEAN))
	setkey(ans, COV_MIN, READL)
	write.csv(ans, row.names=FALSE, file=paste0(outfile.base,'NGSoutput_info.csv'))
	
	ggplot(subset(ans, READL>1), aes(x=READL, y=N/max(N), group=COV_MIN, pch=COV_MIN)) + 
			geom_line(colour='grey50') +
			geom_point() + 			 
			scale_y_continuous(label=scales:::percent) +
			theme_bw() + theme(legend.position='bottom') +
			labs(x='\nmin read length', y='subjects retained\n',pch='min sequencing\ndepth')
	ggsave(file=paste0(outfile.base,'subjectsretained_by_minreads.pdf'), w=3.5, h=6)
	ggplot(subset(ans, READL>1), aes(x=READL, y=HCOV_MEAN*9719, group=COV_MIN, pch=COV_MIN)) + 
			geom_line(colour='grey50') +
			geom_point() + 			 
			scale_y_continuous() +
			theme_bw() + theme(legend.position='bottom') +
			labs(x='\nmin read length', y='average coverage of HIV-1 genome\n(nt)',pch='min sequencing\ndepth')
	ggsave(file=paste0(outfile.base,'coverage_by_minreads.pdf'), w=3.5, h=6)
	
	
	#	roughly where are these 250 bp reads?
	bam.cov250.30	<- subset(bam.cov250, COV>=30)
	bam.cov250.30[, RID:= SID]	
	tmp				<- bam.cov250.30[, list(SUM_REP=sum(REP)), by=c('SID','RID')]
	tmp				<- subset(tmp, SUM_REP>700)
	tmp				<- tmp[, list(SID=SID[which.max(SUM_REP)]), by='RID']
	bam.cov250.30	<- merge(tmp, bam.cov250.30, by=c('RID','SID'))
	setkey(bam.cov250.30, RID, POS)
	tmp				<- bam.cov250.30[, 	{
				z	<- rep(COV,REP)
				list(COV=z, POS2=POS+seq_along(z)-1L)
			}, by=c('SID','RID','POS')]
	tmp			<- tmp[, list(N=length(RID)), by='POS2']
	ggplot(tmp, aes(x=POS2, y=N)) + geom_area() +
			theme_bw() +
			scale_x_continuous(expand=c(0,0)) +
			scale_y_continuous(expand=c(0,0)) +
			#coord_cartesian(ylim=c(0,NA)) +
			labs(x='\nposition on HIV-1 genome', y='individuals\nwith NGS reads at position\n') +
			theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
	ggsave(file=paste0(outfile.base,'horizontal_coverage_histogram_minNGSoutput_reads250.pdf'), w=8, h=2)
	
	#	who are we losing due to some drop off in quality?
	#
	df			<- copy(bam.cov350)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	dcast.data.table(bam.covm, SID~COV_MIN, value.var='HCOV')
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