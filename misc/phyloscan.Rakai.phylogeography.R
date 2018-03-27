
RakaiFull.phylogeography.170421<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170428_withmetadata.rda"		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170428_"	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_withmetadata.rda"		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_"	
	zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	load(infile)	
	#set(rtpdm, NULL, c('MALE_FIRSTPOSVIS.y','MALE_FIRSTPOSDATE.y','FEMALE_FIRSTPOSVIS.y','FEMALE_FIRSTPOSDATE.y'), NULL)
	#setnames(rtpdm, c('MALE_FIRSTPOSVIS.x','MALE_FIRSTPOSDATE.x','FEMALE_FIRSTPOSVIS.x','FEMALE_FIRSTPOSDATE.x'), c('MALE_FIRSTPOSVIS','MALE_FIRSTPOSDATE','FEMALE_FIRSTPOSVIS','FEMALE_FIRSTPOSDATE'))
	nrow(rtpdm)
	#	stage 1: 307 transmissions with direction resolved to 307 recipients with unique transmitter
	#	stage 2: 252 transmissions with unique transmitters	
	
	rtpdm[, AGEDIFF:= rtpdm[, FEMALE_BIRTHDATE-MALE_BIRTHDATE]]
	set(rtpdm, NULL, 'PAIR_ID', rtpdm[, paste0(MALE_RID,'-',FEMALE_RID)])
	set(rtpdm, NULL, 'MALE_SEX', 'M')
	set(rtpdm, NULL, 'FEMALE_SEX', 'F')	
	
	subset(rsm, MIN_PNG_OUTPUT>0)[, quantile(MIN_PNG_OUTPUT/HIV, p=c(0,0.5,1))]
	
	#
	#	some helper data.tables
	rmf		<- subset(rtpdm, TYPE=='mf')
	rfm		<- subset(rtpdm, TYPE=='fm')
	rtr2	<- copy(rmf)
	setnames(rtr2,colnames(rtr2),gsub('FEMALE','REC',colnames(rtr2)))
	setnames(rtr2,colnames(rtr2),gsub('MALE','TR',colnames(rtr2)))
	tmp		<- copy(rfm)
	setnames(tmp,colnames(tmp),gsub('FEMALE','TR',colnames(tmp)))
	setnames(tmp,colnames(tmp),gsub('MALE','REC',colnames(tmp)))
	rtr2	<- rbind(rtr2,tmp)
	
	#
	#	Bayesian source attriution model
	#	adjust for incomplete sampling
	#
	dc	<- rtr2[, list(TR_OBS=length(PAIR_ID)), by=c('TR_COMM_NUM_A','REC_COMM_NUM_A')]	
	#	Bayesian model: add uniform prior
	if(0)
	{
		#	(which is Dirichlet 1 among all communities pairs that have a connection either way
		tmp	<- subset(dc, select=c(REC_COMM_NUM_A, TR_COMM_NUM_A))	
		setnames(tmp, c('REC_COMM_NUM_A','TR_COMM_NUM_A'), c('TR_COMM_NUM_A','REC_COMM_NUM_A'))
		tmp	<- merge(tmp, dc, all.x=1)
		tmp	<- subset(tmp, is.na(TR_OBS))
		set(tmp, NULL, 'TR_OBS', 0)
		dc	<- rbind(dc, tmp)		
	}
	if(1)
	{
		#	This is a bit non-standard, I just don t want the prior to have a large impact, so I chose a sparse one. 
		#	(which is Dirichlet 1 among all communities pairs that are closest)
		#	always add self if not present
		tmp	<- subset(dc[, list(UNOBSERVED_SELF=!any(REC_COMM_NUM_A==TR_COMM_NUM_A)), by='TR_COMM_NUM_A'],UNOBSERVED_SELF, TR_COMM_NUM_A)
		tmp[, REC_COMM_NUM_A:=TR_COMM_NUM_A]
		tmp[, TR_OBS:=0]
		dc	<- rbind(dc, tmp)
		#	ensure each community has at least 1 non-self community, if not add closest other community
		#	I really want to keep this sparse, so do not consider non-self to all communities
		tmp	<- unique(zc, by='COMM_NUM_A')
		tmp	<- as.data.table(t(sapply(seq_len(nrow(tmp)), function(i)
								{
									z<- sort( sqrt( (tmp[,longitude]-tmp[i,longitude])^2+(tmp[,latitude]-tmp[i,latitude])^2 ), index.return=TRUE)$ix
									c('TR_COMM_NUM_A'=tmp[i, COMM_NUM_A], 'REC_COMM_NUM_A'=tmp[z[2],COMM_NUM_A])
								})))	
		z	<- subset(dc[, list(REC_N_OBS=length(REC_COMM_NUM_A)), by='TR_COMM_NUM_A'], REC_N_OBS==1)
		tmp	<- merge(tmp, z, by='TR_COMM_NUM_A')
		tmp	<- merge(subset(tmp, select=c(TR_COMM_NUM_A, REC_COMM_NUM_A)), dc, all.x=1, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A'))
		tmp	<- subset(tmp, is.na(TR_OBS))
		set(tmp, NULL, 'TR_OBS', 0)
		dc	<- rbind(dc, tmp)
	}
	#rsm[, list(ELIGIBLE_AVG=sum(ELIGIBLE_AVG)), by='COMM_TYPE']
	dc[, TR_PRIOR:= 0.5]
	#
	#	Bayesian model first hierarchy: define Beta posterior for sampling probabilities (all alpha and betas)
	#
	tmp	<- subset(rsm, select=c(COMM_NUM_A, ELIGIBLE_AVG, PARTICIPATED_AVG, HIV, MIN_PNG_OUTPUT))
	tmp[, P_PART_EMP:= PARTICIPATED_AVG/ELIGIBLE_AVG]
	tmp[, P_PART_ALPHA:= round(PARTICIPATED_AVG)+1]
	tmp[, P_PART_BETA:= round(ELIGIBLE_AVG-PARTICIPATED_AVG)+1]
	tmp[, P_SEQ_EMP:= MIN_PNG_OUTPUT/HIV]
	tmp[, P_SEQ_ALPHA:= round(MIN_PNG_OUTPUT)+1]
	tmp[, P_SEQ_BETA:= round(HIV-MIN_PNG_OUTPUT)+1]	
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc	<- merge(dc, tmp, by='TR_COMM_NUM_A')
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc	<- merge(dc, tmp, by='REC_COMM_NUM_A')
	#
	#	Bayesian model second hierarchy: draw unobserved data to augment likelihood
	#
	mc.it	<- 1e4
	dcb		<- dc[, {
				tmp	<- 	rbeta(mc.it, TR_P_PART_ALPHA, TR_P_PART_BETA)*
						rbeta(mc.it, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA)*
						rbeta(mc.it, REC_P_PART_ALPHA, REC_P_PART_BETA)*
						rbeta(mc.it, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA)
				#print(tmp)
				tmp	<- rnbinom(mc.it, TR_OBS+TR_PRIOR, tmp)
				#print(tmp)
				list(MONTE_CARLO_IT=seq_len(mc.it), TR_PRIOR=TR_PRIOR, TR_OBS=TR_OBS, TR_MISS= tmp)
			}, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A')]	
	#
	#	Bayesian model second hierarchy: Dirichlet posterior for transmission from community i to j, pi_ij with pi_ij summing to 1
	#
	tmp		<- dcb[, list(	REC_COMM_NUM_A= REC_COMM_NUM_A, 
					TR_COMM_NUM_A= TR_COMM_NUM_A, 
					PI_IJ_ALPHA= TR_OBS+TR_MISS+TR_PRIOR				
			), by='MONTE_CARLO_IT']
	dcb		<- merge(dcb, tmp, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A','MONTE_CARLO_IT'))
	#
	#	this is the end of the source attribution inference on the WAIFM matrix
	#
	
	#
	#	transmission hubs
	#	proportion of transmitters in community i that have a recipient outside community
	#	this is a summary on the estimated posterior density of the WAIFM matrix
	#	TODO missing community 94
	#
	#	get parameters of posterior under augmented likelihood
	z		<- dcb[, {
				z<- which(REC_COMM_NUM_A==TR_COMM_NUM_A)
				list(P_RECOUTSIDE_ALPHA= sum(PI_IJ_ALPHA[-z]), P_RECOUTSIDE_BETA= sum(PI_IJ_ALPHA[z]))
			}, by=c('TR_COMM_NUM_A','MONTE_CARLO_IT')]	
	#	aggregate and get quantiles
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- quantile(rbeta(length(P_RECOUTSIDE_ALPHA)*mc.it, P_RECOUTSIDE_ALPHA, P_RECOUTSIDE_BETA), p=seq(0,1,0.01))
				list(P=seq(0,1,0.01), Q=unname(tmp))
			}, by='TR_COMM_NUM_A']
	#	subset to main quantities of interest
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_NUM_A~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('REC_P_OUTSIDE_CL','REC_P_OUTSIDE_IL','REC_P_OUTSIDE_M','REC_P_OUTSIDE_IU','REC_P_OUTSIDE_CU'))
	#	add TR_OBS from TR_COMM_NUM_A
	z		<- merge(z, dc[, list(REC_N_OBS=sum(TR_OBS)), by='TR_COMM_NUM_A'], by='TR_COMM_NUM_A')
	setnames(z, 'TR_COMM_NUM_A', 'COMM_NUM_A')	
	#	now after analysis, remove 94 with unknown long / lat
	#z		<- subset(z, COMM_NUM!='94')
	z	<- merge(z, unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE))), by='COMM_NUM_A')
	ggplot(z, aes(x=COMM_NUM_A, middle=REC_P_OUTSIDE_M, min=REC_P_OUTSIDE_CL, max=REC_P_OUTSIDE_CU, lower=REC_P_OUTSIDE_IL, upper=REC_P_OUTSIDE_IU, fill=REC_P_OUTSIDE_M)) + 
			geom_boxplot(outlier.shape=NA, stat='identity') +
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +
			#geom_point(data=subset(z, BS==0), pch=2, colour='red') +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.2)) +
			facet_grid(~COMM_TYPE, scales='free_x', space='free_x') +
			labs(x='\ncommunity', y="proportion of transmitters from the community\nthat have a recipient outside the community") +
			theme_bw() + guides(fill='none')
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_adjusted_confidenceintervals_anonmyized.pdf'), w=9, h=5)		
	
	#
	#	transmission hubs - central estimates
	#	among transmitters how many have a recipient outside community
	#	adjusted for sampling
	#	TODO missing community 94
	z[, REC_P_OUTSIDE_PREC:= 1 - REC_P_OUTSIDE_IU + REC_P_OUTSIDE_IL]	
	z	<- merge(unique(subset(zc, select=c(COMM_NUM_A, longitude, latitude, LONG_A, LAT_A)), by='COMM_NUM_A'), z, by='COMM_NUM_A')	
	ggmap(zm) +
			geom_point(data=z, aes(x=longitude, y=latitude, pch=COMM_TYPE, size=REC_N_OBS, fill=100*REC_P_OUTSIDE_M), stroke=1.5, alpha=0.8) +
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +			
			scale_size(breaks=c(5,10,20,40,80), range=c(5,20))+
			scale_shape_manual(values=c('agrarian'=21, 'trading'=22, 'fisherfolk'=23)) +
			#scale_colour_brewer(palette='Dark2') +
			theme(legend.position='bottom', legend.box = "vertical") +			
			labs(	size="transmissions with transmitters from community",
					pch="community type",
					fill="proportion of transmitters from the community\nthat have a recipient outside the community")
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_adjusted.pdf'), w=10, h=10)
	ggmap(zm) +
			geom_point(data=z, aes(x=LONG_A, y=LAT_A, pch=COMM_TYPE, size=REC_N_OBS, fill=100*REC_P_OUTSIDE_M), stroke=1.5, alpha=0.8) +
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +			
			scale_size(breaks=c(5,10,20,40,80), range=c(5,20))+
			scale_shape_manual(values=c('agrarian'=21, 'trading'=22, 'fisherfolk'=23)) +
			#scale_colour_brewer(palette='Dark2') +
			theme(legend.position='bottom', legend.box = "vertical") +			
			labs(	size="transmissions with transmitters from community",
					pch="community type",
					fill="proportion of transmitters from the community\nthat have a recipient outside the community")
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_adjusted_anonymized.pdf'), w=10, h=10)
	
	#
	#	recipient hubs - size is 1-IQR
	ggmap(zm) +
			geom_point(data=z, aes(x=longitude, y=latitude, pch=COMM_TYPE, fill=100*REC_P_OUTSIDE_M, size=REC_P_OUTSIDE_PREC), stroke=1.5, alpha=0.8) +
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +
			#scale_fill_distiller(palette = "RdBu") +
			scale_size(breaks=c(.6,.7,.8,.9,.95), range=c(1,10))+
			scale_shape_manual(values=c('agrarian'=21, 'trading'=22, 'fisherfolk'=23)) +			
			theme(legend.position='bottom', legend.box = "vertical") +
			labs(	size="1 - interquantile range",
					pch="community type",
					fill="proportion of transmitters from the community\nthat have a recipient outside the community")
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_adjusted_sizeIsUnctertainty.pdf'), w=10, h=10)
	ggmap(zm) +
			geom_point(data=z, aes(x=LONG_A, y=LAT_A, pch=COMM_TYPE, fill=100*REC_P_OUTSIDE_M, size=REC_P_OUTSIDE_PREC), stroke=1.5, alpha=0.8) +
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +
			#scale_fill_distiller(palette = "RdBu") +
			scale_size(breaks=c(.6,.7,.8,.9,.95), range=c(1,10))+
			scale_shape_manual(values=c('agrarian'=21, 'trading'=22, 'fisherfolk'=23)) +			
			theme(legend.position='bottom', legend.box = "vertical") +
			labs(	size="1 - interquantile range",
					pch="community type",
					fill="proportion of transmitters from the community\nthat have a recipient outside the community")
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_adjusted_sizeIsUnctertainty_anonymized.pdf'), w=10, h=10)
	
	#
	#	transmission hubs
	#	crude estimates		
	z	<- copy(rtr2)
	z	<- z[, list(REC_OUTSIDE_COMM=length(which(REC_COMM_NUM_A!=TR_COMM_NUM_A)), REC_IN_COMM=length(which(REC_COMM_NUM_A==TR_COMM_NUM_A)) ), by=c('TR_COMM_NUM_A','TR_COMM_TYPE')]
	z[, REC_N:=REC_OUTSIDE_COMM+REC_IN_COMM]
	setnames(z, 'TR_COMM_NUM_A', 'COMM_NUM_A')
	z	<- merge(unique(zc, by='COMM_NUM_A'), z, by='COMM_NUM_A')
	ggmap(zm) +
			geom_point(data=z, aes(x=longitude, y=latitude, pch=COMM_TYPE, colour=COMM_TYPE, size=REC_N, fill=100*REC_OUTSIDE_COMM/REC_N), stroke=1.5, alpha=0.8) + 
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +
			scale_size(breaks=c(5,10,20,40,80), range=c(5,20))+
			scale_shape_manual(values=c('agrarian'=21, 'trading'=22, 'fisherfolk'=23)) +
			scale_colour_brewer(palette='Dark2') +			
			theme(legend.position='bottom', legend.box = "vertical") +
			labs(	size="transmissions with transmitters from community", 
					fill="proportion of transmitters from the community\nthat have a recipient outside the community")
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_crude.pdf'), w=10, h=10)
	
	
	#
	#	geography transmitter flows into agrarian/trading/fisherolk recipient communities
	#	crude	
	tmp		<- rtr2[,list(N=length(unique(PAIR_ID))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]
	tmp[, P_CELL:= N/sum(N)]
	tmp		<- merge(tmp, tmp[, list(P_REC= N/sum(N), REC_COMM_TYPE=REC_COMM_TYPE), by='TR_COMM_TYPE'], by=c('TR_COMM_TYPE','REC_COMM_TYPE'))
	tmp		<- merge(tmp, tmp[, list(P_TR= N/sum(N), TR_COMM_TYPE=TR_COMM_TYPE), by='REC_COMM_TYPE'], by=c('TR_COMM_TYPE','REC_COMM_TYPE'))
	#tmp[, LABEL:= paste0(N, ' (',round(P_CELL,d=2)*100,'%)\ntransmitters: ',round(P_TR,d=2)*100,'%\nrecipients: ',round(P_REC,d=2)*100,'%')]
	tmp[, LABEL:= paste0(N, '\n(',round(P_TR,d=2)*100,'%)')]
	ggplot(tmp, aes(x=REC_COMM_TYPE, y=P_TR, fill=TR_COMM_TYPE)) + 
			geom_bar(stat='identity', position='dodge') +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + 
			scale_y_continuous(labels=scales::percent, limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.2)) +			
			labs(x='\nlocation likely recipient',y='location likely transmitter\n',fill='transmitter from') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_sources_crude_prop.pdf'), w=6, h=5)
	ggplot(tmp, aes(x=REC_COMM_TYPE, y=N, fill=TR_COMM_TYPE)) + 
			geom_bar(stat='identity', position='dodge') +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + 
			scale_y_continuous() +			
			labs(x='\nlocation likely recipient',y='location likely transmitter\n',fill='transmitter from') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_sources_crude_count.pdf'), w=6, h=5)	
	ggplot(tmp, aes(x=TR_COMM_TYPE, y=P_REC, fill=REC_COMM_TYPE)) + 
			geom_bar(stat='identity', position='dodge') +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + 
			scale_y_continuous(labels=scales::percent, limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.2)) +			
			labs(x='\nlocation likely transmitters',y='location likely recipients\n',fill='recipient in') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_destinations_crude_prop.pdf'), w=6, h=5)
	ggplot(tmp, aes(x=TR_COMM_TYPE, y=N, fill=REC_COMM_TYPE)) + 
			geom_bar(stat='identity', position='dodge') +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + 
			scale_y_continuous() +			
			labs(x='\nlocation likely transmitter',y='location likely recipient\n',fill='recipient in') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_destinations_crude_count.pdf'), w=6, h=5)	
	
	#
	#	geography transmitter flows into agrarian/trading/fisherolk recipient communities
	#	adjusted		
	groups	<- c('agrarian','trading','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- subset(zc, COMM_TYPE==group)$COMM_NUM_A		
				z		<- subset(dcb, REC_COMM_NUM_A%in%tmp)[, list(PI_ITYPE_ALPHA= sum(PI_IJ_ALPHA)), by=c('TR_COMM_NUM_A','MONTE_CARLO_IT')]
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				z		<- z[, list(PI_TYPETYPE_ALPHA=sum(PI_ITYPE_ALPHA)), by=c('TR_COMM_TYPE','MONTE_CARLO_IT')]	
				#	aggregate and get quantiles
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_TYPETYPE_ALPHA)
							colnames(tmp)	<- TR_COMM_TYPE
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				z		<- melt(z, id.vars='MONTE_CARLO_IT', variable.name='TR_COMM_TYPE', value.name='PI_TYPETYPE')
				z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(PI_TYPETYPE, p=seq(0,1,0.01)))), by='TR_COMM_TYPE']
				#	subset to main quantities of interest
				z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_TYPE~P, value.var='Q')
				z[, REC_COMM_TYPE:=group]
				z
			})
	z		<- do.call('rbind',z)
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('TR_CL','TR_IL','TR_M','TR_IU','TR_CU'))
	ggplot(z, aes(x=REC_COMM_TYPE, y=TR_M, ymin=TR_CL, lower=TR_IL, upper=TR_IU, ymax=TR_CU, fill=TR_COMM_TYPE)) + 
			#geom_boxplot(stat='identity', position=position_dodge(width=0.9)) +
			geom_bar(position='dodge', stat='identity') +
			geom_errorbar(width=0.2, position=position_dodge(width=0.9)) +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous(labels=scales::percent, limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.2)) +			
			labs(x='\nlocation likely recipient',y='sources of transmissions\n',fill='transmitter from') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_sources_adjusted.pdf'), w=6, h=5)	
	#
	#	geography transmitter flows from agrarian/trading/fisherolk transmitter communities
	#	adjusted		
	groups	<- c('agrarian','trading','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- subset(zc, COMM_TYPE==group)$COMM_NUM_A		
				z		<- subset(dcb, TR_COMM_NUM_A%in%tmp)[, list(PI_ITYPE_ALPHA= sum(PI_IJ_ALPHA)), by=c('REC_COMM_NUM_A','MONTE_CARLO_IT')]
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(z, tmp, by='REC_COMM_NUM_A')
				z		<- z[, list(PI_TYPETYPE_ALPHA=sum(PI_ITYPE_ALPHA)), by=c('REC_COMM_TYPE','MONTE_CARLO_IT')]	
				#	aggregate and get quantiles
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_TYPETYPE_ALPHA)
							colnames(tmp)	<- REC_COMM_TYPE
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				z		<- melt(z, id.vars='MONTE_CARLO_IT', variable.name='REC_COMM_TYPE', value.name='PI_TYPETYPE')
				z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(PI_TYPETYPE, p=seq(0,1,0.01)))), by='REC_COMM_TYPE']
				#	subset to main quantities of interest
				z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), REC_COMM_TYPE~P, value.var='Q')
				z[, TR_COMM_TYPE:=group]
				z
			})
	z		<- do.call('rbind',z)
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('REC_CL','REC_IL','REC_M','REC_IU','REC_CU'))
	ggplot(z, aes(x=TR_COMM_TYPE, y=REC_M, ymin=REC_CL, lower=REC_IL, upper=REC_IU, ymax=REC_CU, fill=REC_COMM_TYPE)) + 
			#geom_boxplot(stat='identity', position=position_dodge(width=0.9)) +
			geom_bar(position='dodge', stat='identity') +
			geom_errorbar(width=0.2, position=position_dodge(width=0.9)) +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous(labels=scales::percent, limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.2)) +			
			labs(x='\nlocation likely transmitter',y='destination of transmissions\n',fill='recipients in') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_destinations_barplot_adjusted.pdf'), w=6, h=5)	
	
	
	
	#
	#	geography flows out > flows in for agrarian/trading/fisherolk communities
	#	adjusted		
	groups	<- c('agrarian','trading','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
				set(tmp, tmp[, which(COMM_TYPE!=group)], 'COMM_TYPE', paste0('non-',group))
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
				z[, FLOW:=paste0('from ',TR_COMM_TYPE,' to ',REC_COMM_TYPE)]
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- FLOW
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	transform as desired
				z[, STAT:= paste0('from ',group,' to non-',group,'\n/\nfrom non-',group,' to ',group)]				
				setnames(z, c(paste0('from ',group,' to non-',group),paste0('from non-',group,' to ',group)), c('A','D'))
				#	get quantiles
				z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(A/D, p=seq(0,1,0.01)))), by='STAT']
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), STAT~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('FLOW_CL','FLOW_IL','FLOW_M','FLOW_IU','FLOW_CU'))
	z[, COMM_TYPE:= gsub('^from ([a-z]+).*$', '\\1', STAT)]
	ggplot(z, aes(x=COMM_TYPE, middle=FLOW_M, min=FLOW_CL, lower=FLOW_IL, upper=FLOW_IU, max=FLOW_CU, fill=STAT)) +
			geom_hline(yintercept=1, colour='grey50', size=1.5) +
			geom_boxplot(stat='identity') +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + theme(legend.position='bottom') +
			scale_y_log10(expand=c(0,0), breaks=c(1/4,1/3,1/2,2/3,1,3/2,2,3,4), labels=c('1/4','1/3','1/2','2/3','1','3/2','2','3','4'), limits=c(1/4,5)) +			
			coord_flip() +
			labs(x='community type\n', y='\nflows out / flows in', fill='flow') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_flowsoutin_adjusted.pdf'), w=6, h=4)
	
	
	#
	#	geography who infects whom matrix  
	#	crude
	tmp		<- rtr2[,list(N=length(unique(PAIR_ID))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]
	tmp[, P_CELL:= N/sum(N)]
	tmp		<- merge(tmp, tmp[, list(P_REC= N/sum(N), REC_COMM_TYPE=REC_COMM_TYPE), by='TR_COMM_TYPE'], by=c('TR_COMM_TYPE','REC_COMM_TYPE'))
	tmp		<- merge(tmp, tmp[, list(P_TR= N/sum(N), TR_COMM_TYPE=TR_COMM_TYPE), by='REC_COMM_TYPE'], by=c('TR_COMM_TYPE','REC_COMM_TYPE'))
	#tmp[, LABEL:= paste0(N, ' (',round(P_CELL,d=2)*100,'%)\ntransmitters: ',round(P_TR,d=2)*100,'%\nrecipients: ',round(P_REC,d=2)*100,'%')]
	tmp[, LABEL:= paste0(N, '\n(',round(P_TR,d=2)*100,'%)')]
	#tmp[, LABEL:= paste0(N, '\n(',round(P_CELL,d=2)*100,'%)')]
	ggplot(tmp, aes(x=factor(REC_COMM_TYPE),y=factor(TR_COMM_TYPE))) + 
			geom_point(aes(size=N), colour='grey80') +
			geom_text(aes(label=LABEL), nudge_x=0, nudge_y=0, size=3, colour='black') +			
			theme_bw() + 
			scale_size(range = c(5, 50)) +
			labs(x='\nlocation likely recipient',y='location likely transmitter\n') +
			guides(size='none')
	ggsave(file=paste0(outfile.base,'_commtype_3x3.pdf'), w=5, h=5)
	
	set(tmp, NULL, 'REC_COMM_TYPE', tmp[,paste0('to_',REC_COMM_TYPE)])
	set(tmp, NULL, 'TR_COMM_TYPE', tmp[,paste0('from_',TR_COMM_TYPE)])
	tmp		<- suppressWarnings(melt(tmp, id.vars=c('REC_COMM_TYPE','TR_COMM_TYPE'), measure.vars=c('N','P_CELL')))
	tmp		<- dcast.data.table(tmp, variable+TR_COMM_TYPE~REC_COMM_TYPE, value.var='value')
	write.csv(tmp, row.names=FALSE, paste0(outfile.base,'_commtype_3x3_raw.csv'))
	#
	#	geography who infects whom matrix  
	#	adjusted P
	tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	z		<- z[, list(TR_OBS=sum(TR_OBS), TR_MISS=sum(TR_MISS), PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0('from_',TR_COMM_TYPE,' to_',REC_COMM_TYPE)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z[, TR_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\1',variable)]
	z[, REC_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\2',variable)]	
	z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(value, p=seq(0,1,0.01)))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]	
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('PADJ_CL','PADJ_IL','PADJ_M','PADJ_IU','PADJ_CU'))		
	ans		<- melt(z, id.vars=c('TR_COMM_TYPE','REC_COMM_TYPE'), measure.vars=c('PADJ_M','PADJ_CL','PADJ_CU'))
	ans		<- dcast.data.table(ans, variable+REC_COMM_TYPE~TR_COMM_TYPE, value.var='value')
	#	adjusted N
	tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	z		<- z[, list(TR_ADJ=sum(TR_OBS)+sum(TR_MISS)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	set(z, NULL, 'REC_COMM_TYPE', z[,paste0('to_',REC_COMM_TYPE)])
	set(z, NULL, 'TR_COMM_TYPE', z[,paste0('from_',TR_COMM_TYPE)])
	z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(TR_ADJ, p=seq(0,1,0.01), type=1))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]	
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('NADJ_CL','NADJ_IL','NADJ_M','NADJ_IU','NADJ_CU'))
	z		<- melt(z, id.vars=c('TR_COMM_TYPE','REC_COMM_TYPE'), measure.vars=c('NADJ_M','NADJ_CL','NADJ_CU'))
	z		<- dcast.data.table(z, variable+REC_COMM_TYPE~TR_COMM_TYPE, value.var='value')
	ans		<- rbind(z, ans)
	write.csv(ans, row.names=FALSE, paste0(outfile.base,'_commtype_3x3_adjusted.csv'))
	
	
	#
	#	geography transmitters from outside community
	z	<- copy(rtr2)
	z[, FROM_OUTSIDE:= factor(REC_COMM_NUM_A!=TR_COMM_NUM_A, levels=c(TRUE,FALSE), labels=c('outside community','same community'))]
	set(z, NULL, 'REC_COMM_TYPE', z[, factor(REC_COMM_TYPE)])
	ggplot(z, aes(x=REC_COMM_TYPE, fill=FROM_OUTSIDE)) + geom_bar() + 
			theme_bw() + 
			labs(x='\ncommunity type of recipient', fill='transmitter from')
	ggsave(file=paste0(outfile.base,'_extra_community.pdf'), w=5, h=5)
	ggplot(z, aes(x=REC_COMM_TYPE, fill=FROM_OUTSIDE)) + geom_bar(position='fill') + 
			theme_bw() + 
			scale_y_continuous(labels=scales::percent, expand=c(0,0), breaks=seq(0,1,0.2)) +
			labs(x='\ncommunity type of recipient', y='transmitters\n', fill='transmitter from')
	ggsave(file=paste0(outfile.base,'_extra_community_props.pdf'), w=5, h=5)
	#
	#	odds for unequal transmission flows agrarian/fisherfolk
	z	<- z[, list(N=length(REC_RID)), by=c('REC_COMM_TYPE','FROM_OUTSIDE')]
	tmp	<- dcast.data.table(subset(z, REC_COMM_TYPE!='trading'), FROM_OUTSIDE~REC_COMM_TYPE, value.var='N')
	zz	<- as.matrix(tmp[, 2:3, with=FALSE])
	rownames(zz)	<- tmp[[1]]
	fisher.test(zz)
	#STAGE 1
	#data:  zz
	#p-value = 0.0113
	#alternative hypothesis: true odds ratio is not equal to 1
	#95 percent confidence interval:
	#1.137942 3.390227
	#sample estimates:
	#odds ratio 
	#1.963628
	
	#STAGE 2
	#p-value = 0.0002936
	#alternative hypothesis: true odds ratio is not equal to 1
	#95 percent confidence interval:
	#1.698675 7.052786
	#sample estimates:
	#odds ratio 
	#3.427573 
	
	#
	#	geography recipients in different community
	z	<- copy(rtr2)
	z[, FROM_OUTSIDE:= factor(REC_COMM_NUM_A!=TR_COMM_NUM_A, levels=c(TRUE,FALSE), labels=c('outside community','same community'))]
	set(z, NULL, 'TR_COMM_TYPE', z[, factor(TR_COMM_TYPE)])
	ggplot(z, aes(x=TR_COMM_TYPE, fill=FROM_OUTSIDE)) + geom_bar() + 
			theme_bw() + 
			labs(x='\ncommunity type of transmitters', fill='recipient')
	ggsave(file=paste0(outfile.base,'_extra_community_of_recipients.pdf'), w=5, h=5)	
	ggplot(z, aes(x=TR_COMM_TYPE, fill=FROM_OUTSIDE)) + geom_bar(position='fill') + 
			theme_bw() + 
			scale_y_continuous(labels=scales::percent, expand=c(0,0), breaks=seq(0,1,0.2)) +
			labs(x='\ncommunity type of transmitters', y='recipients\n', fill='recipient')
	ggsave(file=paste0(outfile.base,'_extra_community_of_recipients_props.pdf'), w=5, h=5)	
	z	<- z[, list(N=length(TR_RID)), by=c('TR_COMM_TYPE','FROM_OUTSIDE')]
	tmp	<- dcast.data.table(subset(z, TR_COMM_TYPE!='trading'), FROM_OUTSIDE~TR_COMM_TYPE, value.var='N')
	zz	<- as.matrix(tmp[, 2:3, with=FALSE])
	rownames(zz)	<- tmp[[1]]
	fisher.test(zz)
	#STAGE 1
	#p-value = 0.05827
	#alternative hypothesis: true odds ratio is not equal to 1
	#95 percent confidence interval:
	#		0.9429887 3.1200202
	#sample estimates:
	#		odds ratio 
	#1.719531
	
	#STAGE 2
	#p-value = 0.003303
	#alternative hypothesis: true odds ratio is not equal to 1
	#95 percent confidence interval:
	# 1.346652 6.087743
	#sample estimates:
	#odds ratio 
	#  2.844896 
	
	
	#
	#	geography transmitters from outside community by gender of recipients
	z	<- copy(rtr2)
	z[, FROM_OUTSIDE:= factor(REC_COMM_NUM_A!=TR_COMM_NUM_A, levels=c(TRUE,FALSE), labels=c('outside community','same community'))]
	set(z, NULL, 'REC_COMM_TYPE', z[, factor(REC_COMM_TYPE)])
	set(z, NULL, 'REC_SEX', z[, paste0('gender recipient: ',REC_SEX)])	
	ggplot(z, aes(x=REC_COMM_TYPE, fill=FROM_OUTSIDE)) + geom_bar(position='fill') + 
			theme_bw() + 
			scale_y_continuous(labels=scales::percent, expand=c(0,0), breaks=seq(0,1,0.2)) +
			labs(x='\ncommunity type of recipient', y='transmitters\n', fill='transmitter from') +
			facet_grid(~REC_SEX)
	ggsave(file=paste0(outfile.base,'_extra_community_bygender_props.pdf'), w=7, h=5)	
	
	#
	#	community locations
	#		
	zc	<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/community_spatialLoc.csv'))
	set(zc, NULL, 'COMM_NUM_A', zc[,gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',as.character(COMM_NUM)))))])
	zc	<- merge(zc, unique(subset(rh, select=c(COMM_NUM, COMM_TYPE))), by='COMM_NUM')	
	ggmap(zm) +
			geom_point(data=unique(zc, by='COMM_NUM'), aes(x=longitude, y=latitude, pch=COMM_TYPE, colour=COMM_TYPE), size=8, alpha=0.8) +
			geom_text(data=unique(zc, by='COMM_NUM'), aes(x=longitude, y=latitude, label=COMM_NUM), nudge_x=0, nudge_y=0, size=3, colour='black')
	ggsave(file=paste0(outfile.base,'_hubs_comm_locations.pdf'), w=10, h=10)
	
	
}

RakaiFull.phylogeography.171122.samplinglocations<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl35_prior23_min30_withmetadata.rda"
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_prior23_min30_withmetadata.rda"
	outfile.base			<- gsub('_withmetadata.rda','',infile)
	
	zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	load(infile)	
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	nrow(rtpdm)
	#	284	with 0.035
	#	238	with 0.025
	
	#
	#	all locations
	ggmap(zm) +
			geom_point(data=unique(zc, by='COMM_NUM'), aes(x=longitude, y=latitude, pch=COMM_TYPE, colour=COMM_TYPE), size=8, alpha=0.8) +
			geom_text(data=unique(zc, by='COMM_NUM'), aes(x=longitude, y=latitude, label=COMM_NUM), nudge_x=0, nudge_y=0, size=3, colour='black')
	ggsave(file=paste0(outfile.base,'_comm_all_locations.pdf'), w=10, h=10)
	
	#
	#	locations with number sequence samples
	ds	<- subset(rsm, MIN_PNG_OUTPUT>0, select=c(COMM_NUM, COMM_TYPE, COMM_NUM_A, ELIGIBLE_AVG, PARTICIPATED_AVG, HIV, MIN_PNG_OUTPUT))
	ds[, P_PART_EMP:= PARTICIPATED_AVG/ELIGIBLE_AVG]
	ds[, P_SEQ_EMP:= MIN_PNG_OUTPUT/HIV]
	ds[, P_SEQCOV:= P_PART_EMP*P_SEQ_EMP]
	ds[, P_SEQCOV_C:= cut(P_SEQCOV, breaks=seq(0,1,0.1), labels=c('<10%','10-19%','20-29%','30-39%','40-49%','50-59%','60-69%','70-79%','80-89%','90%-100%'))]
	tmp	<- unique(subset(zc, select=c(COMM_NUM, longitude, latitude)), by='COMM_NUM')
	ds	<- merge(ds, tmp, by='COMM_NUM')
	ggmap(zm) +
			geom_point(data=ds, aes(x=longitude, y=latitude, pch=COMM_TYPE, colour=P_SEQCOV_C, size=MIN_PNG_OUTPUT), alpha=0.8) +
			geom_text(data=ds, aes(x=longitude, y=latitude, label=MIN_PNG_OUTPUT), nudge_x=0, nudge_y=0, size=3, colour='black') +			
			scale_colour_brewer(palette="YlOrRd") +
			scale_size(breaks=c(0,10,50,100,200,400,1000), range=c(5,15))+
			#scale_shape_manual(values=c('agrarian'=21, 'trading'=22, 'fisherfolk'=23)) +			
			theme(legend.position='bottom', legend.box = "vertical") +
			guides(size='none', colour=guide_legend(override.aes=list(size=7)), pch=guide_legend(override.aes=list(size=7))) +
			labs(	x='', y='', pch="community\ntype",
					colour="estimated\nsequence coverage\nof HIV infected population")
	ggsave(file=paste0(outfile.base,'_comm_sequencenumbers_on_map.pdf'), w=7, h=7)
	
	
	#	sequence coverage
	tmp	<- rsm[,  list(STAT=c('Median','CL','CU'), V=as.numeric(binconf(MIN_PNG_OUTPUT, round(HIV*6205/4928)))) , by=c('COMM_NUM_A','COMM_TYPE','MIN_PNG_OUTPUT')]
	tmp	<- dcast.data.table(tmp, COMM_NUM_A+COMM_TYPE+MIN_PNG_OUTPUT~STAT, value.var='V')
	ggplot(tmp, aes(x=COMM_NUM_A, y=Median, ymin=CL, ymax=CU)) +
			geom_point(aes(size=MIN_PNG_OUTPUT)) +
			#geom_errorbar(width=0.3) +
			theme_bw() +
			scale_y_continuous(labels=scales:::percent, lim=c(0.1,0.8)) +
			labs(x='\ncommunity', y='sequence coverage\nHIV infected who have a sequence\n', size='NGS output') +
			facet_grid(~COMM_TYPE, space='free', scales='free')
	ggsave(file=paste0(outfile.base,'_comm_seqcov.pdf'), w=10, h=5)
}


RakaiFull.phylogeography.171122.participation.bias<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/rakai_elibility.rda"
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	load(infile)
	
	#	subset to data of interest
	de	<- as.data.table(eldat)	
	de	<- subset(de, status%in%c('_Participated','Away','Blood refusal','Missing data','Other','Refused','urine sample'))
	de	<- subset(de, visit<17)
	de	<- subset(de, select=c(PERM_ID, CURR_ID, visit, RESIDENT, MOBILITY, REGION, COMM_NUM, HH_NUM, SEX, STUDY_ID, status, date, AGEYRS))
	setnames(de, colnames(de), toupper(colnames(de)))
	#	define PARTICIPATED as "participated and gave blood"
	de[, PARTICIPATED:= as.integer(STATUS%in%c('_Participated','Missing data'))]
	
	#	merge communities that are very close / identical
	setnames(de, 'COMM_NUM', 'COMM_NUM_RAW')
	set(de, NULL, 'COMM_NUM', de[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',COMM_NUM_RAW))))])
	
	#	fixup missing permanent IDs with Study IDs
	stopifnot( nrow(subset(de, is.na(PERM_ID) & is.na(STUDY_ID)))==0 )
	tmp		<- de[, which(is.na(PERM_ID))]
	set(de, tmp, 'PERM_ID', de[tmp, STUDY_ID])
	#	find study participants with multiple permanent IDs and fixup
	tmp		<- subset(de, !is.na(STUDY_ID))[, list(N_PERM=length(unique(PERM_ID))), by='STUDY_ID']
	tmp		<- subset(tmp, N_PERM>1, STUDY_ID)
	tmp[, DUMMY:=1]
	de		<- merge(de, tmp, by='STUDY_ID', all.x=1)
	tmp		<- de[, which(DUMMY==1)]
	set(de, tmp, 'PERM_ID', de[tmp, STUDY_ID])
	de[, DUMMY:=NULL]
	#	calculate number participated yes/no by gender and community
	#	we care about wether individuals participated in any of the visits 15, 15.1, 16
	#	the best we can do is look at permanent ids, even though they are not unique
	tmp		<- de[, list(PARTICIPATED_ANY_VISIT=as.integer(any(PARTICIPATED==1))), by='PERM_ID']
	de		<- merge(de, tmp, by='PERM_ID')
	#	to those of ever participated, find study id
	tmp		<- subset(de, PARTICIPATED_ANY_VISIT==1)
	stopifnot( nrow(subset(tmp, PARTICIPATED==1 & is.na(STUDY_ID)))==0 )
	tmp		<- tmp[, list(DUMMY= STUDY_ID[PARTICIPATED==1][1]), by='PERM_ID']
	de		<- merge(de, tmp, by='PERM_ID', all.x=TRUE)
	tmp		<- de[, which(is.na(STUDY_ID) & !is.na(DUMMY))]
	set(de, tmp, 'STUDY_ID', de[tmp,DUMMY])
	de[, DUMMY:=NULL]
	setnames(de, 'STUDY_ID','RID')
	#	there are a few individuals with Study ID who did not "participate" in the sense above
	#subset(de, PARTICIPATED_ANY_VISIT==0 & !is.na(STUDY_ID))
	
	#	prepare HIV status
	tmp		<- RakaiCirc.epi.get.info.170208()
	rd		<- tmp$rd
	rneuro	<- tmp$rn
	rd		<- unique(subset(rd, !is.na(FIRSTPOSDATE), select=c(RID, SEX, VISIT, DATE, FIRSTPOSDATE, PANGEA, COMM_NUM, COMM_NUM_A, BIRTHYR, EST_DATEDIED)))
	#	select community etc closest to first pos date
	tmp		<- rd[, list(VISIT=VISIT[which.min(abs(DATE-FIRSTPOSDATE))]), by='RID']
	tmp2	<- rd[, list(PANGEA=as.integer(any(PANGEA==1))), by='RID']
	rd		<- merge(rd, tmp, by=c('RID','VISIT'))
	rd[, HIV:= 1L]
	rd[, PANGEA:=NULL]
	rd		<- merge(rd, tmp2, by='RID')
	#rd		<- subset(rd, BIRTHYR>2010-50 & (is.na(EST_DATEDIED) | (!is.na(EST_DATEDIED) & EST_DATEDIED>2010)))
		
	#	add HIV status
	de		<- merge(de, subset(rd, select=c(RID, HIV)), by='RID', all.x=1)	
	
	#	get individuals with at least 750nt overlap with another individual at 20X
	infile	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170704_cl35_withmedian.rda"
	load(infile)
	rn	<- subset(rn, NEFF>3)	
	rn	<- data.table(RID=rn[, unique(c(ID1, ID2))],MIN_PNG_OUTPUT=1)		
	#	add to rn individuals with any sequence data
	load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75_part2.rda')
	tmp		<- unique(subset(dc, !is.na(SID), select=RID))
	tmp[, BAM_OUTPUT:=1L]
	rn		<- merge(tmp, rn, all=1, by='RID')
	set(rn, rn[, which(is.na(MIN_PNG_OUTPUT))],'MIN_PNG_OUTPUT',0L)
	#	add
	de		<- merge(de, rn, by='RID', all.x=1)
	rd		<- merge(rd, rn, by='RID', all.x=1)
	rneuro	<- merge(rneuro, rn, by='RID', all.x=1)
		
	#	for every individual that ever participated, keep first record of when participated
	#	for every individual that never participated, keep last record
	tmp		<- subset(de, PARTICIPATED_ANY_VISIT==0)[, list(VISIT=max(VISIT)), by='PERM_ID']
	tmp		<- merge(de, tmp, by=c('PERM_ID','VISIT'))
	tmp2	<- subset(de, PARTICIPATED_ANY_VISIT==1)[, list(VISIT=min(VISIT[PARTICIPATED==1])), by='PERM_ID']
	de		<- merge(de, tmp2, by=c('PERM_ID','VISIT'))
	de		<- rbind(de, tmp)
	
	#	fill in missing HIV etc
	tmp		<- de[, which(PARTICIPATED_ANY_VISIT==1 & is.na(HIV))]
	set(de, tmp, 'HIV', 0L)
	tmp		<- de[, which(PARTICIPATED_ANY_VISIT==1 & is.na(BAM_OUTPUT))]
	set(de, tmp, c('BAM_OUTPUT','MIN_PNG_OUTPUT'), 0L)
	tmp		<- rd[, which(is.na(BAM_OUTPUT))]
	set(rd, tmp, c('BAM_OUTPUT','MIN_PNG_OUTPUT'), 0L)
	
	rneuro[, sum(MIN_PNG_OUTPUT)]			# 224	
	rd[, sum(MIN_PNG_OUTPUT)]				# 2708
	de[, sum(MIN_PNG_OUTPUT, na.rm=TRUE)]	# 2652
	
	#	calculate number for whom we have deep seq output
	tmp		<- subset(rd, MIN_PNG_OUTPUT==1)[, list(DEEP_SEQ=length(RID)), by=c('COMM_NUM')]	
	tmp		<- subset(tmp, DEEP_SEQ>0)
	tmp[, DEEP_SEQ:=NULL]	
	tmp[, DUMMY:= 1]
	de		<- merge(de, tmp, by='COMM_NUM', all.x=1)
	rd		<- merge(rd, tmp, by='COMM_NUM', all.x=1)
		
	#	subset to only those communities with seq data
	#	4 fishing, 34 inland
	def		<- copy(de)
	de		<- subset(de, !is.na(DUMMY))
	de[, DUMMY:=NULL]
	rd		<- subset(rd, !is.na(DUMMY))
	rd[, DUMMY:=NULL]	
	tmp		<- unique(subset(rd, select=c(COMM_NUM,COMM_NUM_A)))
	de		<- merge(de, tmp, by='COMM_NUM')
	
	#	add age category at midpoint of observation period, 2013.25
	set(de, de[, which(is.na(DATE))], 'DATE', as.Date("2011-09-02"))	
	de[, AGE_AT_MID:= 2013.25 - (hivc.db.Date2numeric(DATE)-AGEYRS)]
	de[, AGE_AT_MID_C:= cut(AGE_AT_MID, breaks=c(10,25,35,65), labels=c('15-24','25-34','35+'), right=FALSE)]
	stopifnot( nrow(subset(de, is.na(AGE_AT_MID_C)))==0 )
	
	#	now calculate participation by community and gender
	des		<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','SEX','PARTICIPATED_ANY_VISIT')]
	set(des, NULL, 'PARTICIPATED_ANY_VISIT', des[, factor(PARTICIPATED_ANY_VISIT, levels=c('1','0'), labels=c('PART_EVER','PART_NEVER'))])
	des		<- dcast.data.table(des, COMM_NUM+COMM_NUM_A+SEX~PARTICIPATED_ANY_VISIT, value.var='N')
	
	#	now calculate HIV pos by end of round 16, and sequenced by end of round 16, by community and gender
	#	among those that participated in rounds 15 - 16
	tmp		<- subset(de, PARTICIPATED_ANY_VISIT==1)
	tmp		<- tmp[, list(HIV_1516_YES=length(which(HIV==1)), HIV_1516_NO=length(which(HIV==0)), DEEP_SEQ_1516=length(which(MIN_PNG_OUTPUT==1)) ), by=c('COMM_NUM','SEX')]
	tmp		<- subset(tmp, DEEP_SEQ_1516>0)	
	des		<- merge(des, tmp, by=c('COMM_NUM','SEX'), all=TRUE)
	
	#	add community type
	des[, COMM_TYPE:= as.character(factor(substr(COMM_NUM_A,1,1)=='f',levels=c(TRUE,FALSE),labels=c('fisherfolk','inland')))]	
	
	#	now do the same by community and gender and age category
	desa	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','SEX','AGE_AT_MID_C','PARTICIPATED_ANY_VISIT')]
	set(desa, NULL, 'PARTICIPATED_ANY_VISIT', desa[, factor(PARTICIPATED_ANY_VISIT, levels=c('1','0'), labels=c('PART_EVER','PART_NEVER'))])
	desa	<- dcast.data.table(desa, COMM_NUM+COMM_NUM_A+SEX+AGE_AT_MID_C~PARTICIPATED_ANY_VISIT, value.var='N')
	tmp		<- subset(de, PARTICIPATED_ANY_VISIT==1)
	tmp		<- tmp[, list(HIV_1516_YES=length(which(HIV==1)), HIV_1516_NO=length(which(HIV==0)), DEEP_SEQ_1516=length(which(MIN_PNG_OUTPUT==1)) ), by=c('COMM_NUM','SEX','AGE_AT_MID_C')]
	tmp		<- subset(tmp, DEEP_SEQ_1516>0)	
	desa	<- merge(desa, tmp, by=c('COMM_NUM','SEX','AGE_AT_MID_C'), all=TRUE)
	desa[, COMM_TYPE:= as.character(factor(substr(COMM_NUM_A,1,1)=='f',levels=c(TRUE,FALSE),labels=c('fisherfolk','inland')))]
	
	#	now calculate HIV pos by end of round 16, and sequenced by end of round 16, by community and gender
	#	among those that ever participated
	#	excluding those that died before 2010, or that reached age 60 before 2010
	tmp		<- subset(rd, is.na(EST_DATEDIED) | (!is.na(EST_DATEDIED) & EST_DATEDIED>=2010))	
	tmp		<- subset(tmp, (2010-BIRTHYR)<60)
	rds		<- tmp[, list(HIV_EVER_YES=length(which(HIV==1)), DEEP_SEQ_EVER=length(which(MIN_PNG_OUTPUT==1)) ), by=c('COMM_NUM','SEX')]
	rds		<- subset(rds, DEEP_SEQ_EVER>0)
	rds[, sum(DEEP_SEQ_EVER)]	# 2665 (among those in rd who are not too old and did not die before 2010)
	
		
	#	there s just 1 person ever sequenced in comms 401 and 55, ignore
	#	this means we are just working with des, not rds
	des		<- merge(des, rds, by=c('COMM_NUM','SEX'), all.x=TRUE)
	des		<- subset(des, !COMM_NUM%in%c('55','401'))
	desa	<- subset(desa, !COMM_NUM%in%c('55','401'))
	#	check for missing and replace by 0		
	for(x in c('PART_EVER','PART_NEVER','HIV_1516_YES','HIV_1516_NO','DEEP_SEQ_1516'))
		set(desa, which(is.na(desa[[x]])), x, 0)
	
	length( setdiff( subset(tmp, !is.na(RID) & HIV==1 & MIN_PNG_OUTPUT==1)[, RID], subset(de, !is.na(RID) & HIV==1 & MIN_PNG_OUTPUT==1)[, RID] ) )
	#	57 individuals that have a sequence who are not in de but in rd
	
	save(des, desa, rds, de, rd, file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender.rda')
	
	#	overall sequence coverage 
	des[, sum(DEEP_SEQ_1516)]	#2652
	des[, sum(HIV_1516_YES)]	#4454
	des[, sum(HIV_EVER_YES)]	#5397
	des[, sum(DEEP_SEQ_1516)/sum(HIV_EVER_YES)]		# 0.4913841
	des[, sum(DEEP_SEQ_1516)/sum(HIV_1516_YES)]		# 0.5954198
	des[, sum(DEEP_SEQ_1516)/( sum(HIV_EVER_YES * (PART_EVER+PART_NEVER)/PART_EVER) )]	#0.3435715
	des[, sum(DEEP_SEQ_1516)/( sum(HIV_1516_YES * (PART_EVER+PART_NEVER)/PART_EVER) )]	#0.4195271
	
	#	make supplementary plots of participation and sequencing by gender
	des[, P_PART_RAW:= PART_EVER/(PART_EVER+PART_NEVER)]
	des[, P_SEQ_RAW_EVER:= DEEP_SEQ_1516/(HIV_EVER_YES)]
	des[, P_SEQ_RAW_1516:= DEEP_SEQ_1516/(HIV_1516_YES)]
	
	
	tmp		<- melt(des, id.vars=c('COMM_TYPE','COMM_NUM_A','SEX'), measure.vars=c('P_PART_RAW','P_SEQ_RAW_EVER','P_SEQ_RAW_1516'))	
	set(tmp, NULL, 'COMM_TYPE', tmp[, factor(COMM_TYPE, levels=c('fisherfolk','inland'), labels=c('fishing site','inland community'))])
	set(tmp, NULL, 'SEX', tmp[, factor(SEX, levels=c('M','F'), labels=c('male','female'))])
	ggplot(tmp, aes(x=COMM_NUM_A, fill=SEX, y=value)) +
			geom_bar(position='dodge', stat='identity') +
			scale_fill_manual(values=c('female'='hotpink2', 'male'='steelblue2')) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			facet_grid(variable~COMM_TYPE, space='free', scales='free') +
			labs(x='\ncommunity', fill='gender') +
			theme_bw()
	ggsave(file=paste0(outfile.base,'sampling_differences_180322.pdf'), w=10, h=10)
	#
	
	#	number communities with some sequences
	des[, length(unique(COMM_NUM))]
	#	36

	#	make table 1 for paper
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_withmetadata.rda"
	outfile.base			<- gsub('_withmetadata.rda','',infile)
	load(infile)
	rtp[, SELECT2:= 1L]
	set(rtp, rtp[, which(grepl('most likely a pair',SELECT))], 'SELECT2', 2L)
	set(rtp, rtp[, which(grepl('mf|fm',SELECT))], 'SELECT2', 3L)
	
	tmp	<- subset(rtp, LINK_MF==1 | LINK_FM==1, select=c(MALE_RID, MALE_COMM_TYPE, MALE_BIRTHDATE, SELECT2))
	setnames(tmp, c('MALE_RID','MALE_COMM_TYPE','MALE_BIRTHDATE'), c('RID','COMM_TYPE','BIRTHDATE'))
	tmp[, SEX:='M']
	tmp2<- subset(rtp, select=c(FEMALE_RID, FEMALE_COMM_TYPE, FEMALE_BIRTHDATE, SELECT2))
	setnames(tmp2, c('FEMALE_RID','FEMALE_COMM_TYPE','FEMALE_BIRTHDATE'), c('RID','COMM_TYPE','BIRTHDATE'))
	tmp2[, SEX:='F']
	tmp	<- rbind(tmp, tmp2)
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')
	set(tmp, tmp[, which(is.na(BIRTHDATE))],'BIRTHDATE',tmp[, mean(BIRTHDATE, na.rm=TRUE)])
	dr	<- tmp[, 	{
						z<- which.max(SELECT2)
						list(SEX=SEX[z], COMM_TYPE=COMM_TYPE[z], SELECT2=SELECT2[z], BIRTHDATE=BIRTHDATE[z])
					}, by='RID']
	dr[, AGE_AT_MID:= 2013.25-BIRTHDATE]
	dr[, AGE_AT_MID_C:= cut(AGE_AT_MID, breaks=c(10,25,35,65), labels=c('15-24','25-34','35+'), right=FALSE)]
			
	ans		<- dr[, list(IN_TRM_CHAIN= length(RID[SELECT2>=1])), by=c('COMM_TYPE','SEX')]
	ans		<- merge(ans, dr[, list(LINKED= length(RID[SELECT2>=2])), by=c('COMM_TYPE','SEX')], by=c('COMM_TYPE','SEX'))
	ans		<- merge(ans, dr[, list(LINKED_DIR= length(RID[SELECT2>=3])), by=c('COMM_TYPE','SEX')], by=c('COMM_TYPE','SEX'))	
	tmp		<- des[, list(	ELIGIBLE_1516= sum(PART_EVER+PART_NEVER), 
							PART_1516= sum(PART_EVER),
							HIV_1516= sum(HIV_1516_YES),
							DEEP_SEQ_1516= sum(DEEP_SEQ_1516)
							), by=c('COMM_TYPE','SEX')]
	ans		<- merge(tmp, ans, by=c('COMM_TYPE','SEX'))	
	
	
	ansa	<- dr[, list(IN_TRM_CHAIN= length(RID[SELECT2>=1])), by=c('COMM_TYPE','SEX','AGE_AT_MID_C')]
	ansa	<- merge(ansa, dr[, list(LINKED= length(RID[SELECT2>=2])), by=c('COMM_TYPE','SEX','AGE_AT_MID_C')], by=c('COMM_TYPE','SEX','AGE_AT_MID_C'))
	ansa	<- merge(ansa, dr[, list(LINKED_DIR= length(RID[SELECT2>=3])), by=c('COMM_TYPE','SEX','AGE_AT_MID_C')], by=c('COMM_TYPE','SEX','AGE_AT_MID_C'))		
	tmp		<- desa[, list(	ELIGIBLE_1516= sum(PART_EVER+PART_NEVER), 
							PART_1516= sum(PART_EVER),
							HIV_1516= sum(HIV_1516_YES),
							DEEP_SEQ_1516= sum(DEEP_SEQ_1516)
							), by=c('COMM_TYPE','SEX','AGE_AT_MID_C')]
	ansa	<- merge(tmp, ansa, by=c('COMM_TYPE','SEX','AGE_AT_MID_C'))
	# add total to ansa
	ansa	<- melt(ansa, id.vars=c('COMM_TYPE','SEX','AGE_AT_MID_C'))
	tmp		<- ansa[, list(SEX=SEX, AGE_AT_MID_C=AGE_AT_MID_C, prop=value/sum(value), LEGEND=paste(value,' (',100*round(value/sum(value),d=2),'%)',sep='')), by=c('COMM_TYPE','variable')]
	ansa	<- merge(ansa, tmp, by=c('variable','COMM_TYPE','SEX','AGE_AT_MID_C'))	
	tmp		<- ansa[, list(SEX='Any', AGE_AT_MID_C='Any', value=sum(value), LEGEND=as.character(sum(value))), by=c('COMM_TYPE','variable')]
	tmp2	<- ansa[, list(COMM_TYPE='Any',SEX='Any', AGE_AT_MID_C='Any', value=sum(value), LEGEND=as.character(sum(value))), by=c('variable')]
	ansa	<- rbind(ansa,tmp,tmp2, fill=TRUE)
	ansa[, LABEL:= paste(COMM_TYPE, SEX, AGE_AT_MID_C, sep='-')]
	ansa	<- dcast.data.table(ansa, LABEL+COMM_TYPE+SEX+AGE_AT_MID_C~variable, value.var='LEGEND')
	
	write.csv(ansa, row.names=FALSE, file=paste0(outfile.base, '_table1.csv'))	
}


RakaiFull.phylogeography.171122.flows.on.map<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl35_prior23_min30_withmetadata.rda"
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_prior23_min30_withmetadata.rda"
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_withmetadata.rda"
	outfile.base			<- gsub('_withmetadata.rda','',infile)
	
	style	<- "feature:road|color:0x17202A&style=feature:water|color:0x677996&style=feature:landscape.natural|color:0xedecda&style=feature:administrative|visibility=off"
	zm		<- get_googlemap(c(lon=31.65, lat=-0.66), scale=2, size=c(550,550), zoom=10, maptype="road", style=style)
	zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	load(infile)	
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	set(rtpdm, NULL, 'MALE_COMM_TYPE', rtpdm[, factor(MALE_COMM_TYPE=='fisherfolk',levels=c(TRUE,FALSE),labels=c('fishing site','inland community') )])
	set(rtpdm, NULL, 'FEMALE_COMM_TYPE', rtpdm[, factor(FEMALE_COMM_TYPE=='fisherfolk',levels=c(TRUE,FALSE),labels=c('fishing site','inland community') )])
	nrow(rtpdm)
	#	284	with 0.035
	#	238	with 0.025
	
	#	locations with number sequence samples
	ds	<- subset(rsm, MIN_PNG_OUTPUT>0, select=c(COMM_NUM, COMM_TYPE, COMM_NUM_A, ELIGIBLE_AVG, PARTICIPATED_AVG, HIV, MIN_PNG_OUTPUT))
	set(ds, NULL, 'COMM_TYPE', ds[, factor(COMM_TYPE=='fisherfolk',levels=c(TRUE,FALSE),labels=c('fishing site','inland community') )])
	ds[, P_PART_EMP:= PARTICIPATED_AVG/ELIGIBLE_AVG]
	ds[, P_SEQ_EMP:= MIN_PNG_OUTPUT/HIV]
	ds[, P_SEQCOV:= P_PART_EMP*P_SEQ_EMP]
	ds[, P_SEQCOV_C:= cut(P_SEQCOV, breaks=seq(0,1,0.1), labels=c('<10%','10-19%','20-29%','30-39%','40-49%','50-59%','60-69%','70-79%','80-89%','90%-100%'))]
	tmp	<- unique(subset(zc, select=c(COMM_NUM, longitude, latitude)), by='COMM_NUM')
	ds	<- merge(ds, tmp, by='COMM_NUM')
	
	#	show location of transmitters, recipients
	rtpdm[, AGEDIFF:= rtpdm[, FEMALE_BIRTHDATE-MALE_BIRTHDATE]]	
	set(rtpdm, NULL, 'MALE_SEX', 'M')
	set(rtpdm, NULL, 'FEMALE_SEX', 'F')	
	rmf		<- subset(rtpdm, grepl('mf',SELECT) )
	rfm		<- subset(rtpdm, grepl('fm',SELECT) )
	rtr2	<- copy(rmf)
	setnames(rtr2,colnames(rtr2),gsub('FEMALE','REC',colnames(rtr2)))
	setnames(rtr2,colnames(rtr2),gsub('MALE','TR',colnames(rtr2)))
	tmp		<- copy(rfm)
	setnames(tmp,colnames(tmp),gsub('FEMALE','TR',colnames(tmp)))
	setnames(tmp,colnames(tmp),gsub('MALE','REC',colnames(tmp)))
	rtr2	<- rbind(rtr2,tmp)	
	dc.tr	<- rtr2[, list(TR_OBS=length(PAIRID)), by=c('TR_COMM_NUM_A')]
	setnames(dc.tr, 'TR_COMM_NUM_A', 'COMM_NUM_A')
	dc.tr	<- merge(dc.tr, ds, by='COMM_NUM_A')
	dc.tr[, TR_OBS_C:= cut(TR_OBS, breaks=c(0,1,2,5,10,20,50,100))]
	dc.rec	<- rtr2[, list(REC_OBS=length(PAIRID)), by=c('REC_COMM_NUM_A')]
	setnames(dc.rec, 'REC_COMM_NUM_A', 'COMM_NUM_A')
	dc.rec	<- merge(dc.rec, ds, by='COMM_NUM_A')
	dc.rec[, REC_OBS_C:= cut(REC_OBS, breaks=c(0,1,2,5,10,20,50,100))]
	
	rtr2[, table(TR_COMM_TYPE, REC_COMM_TYPE)]
	#                  REC_COMM_TYPE
	#TR_COMM_TYPE       fishing site inland community
  	#fishing site              170                9
  	#inland community           14              100
	#
	#	plot number of observed recipients and observed transmitters
	ggmap(zm) +
			geom_point(data=dc.rec, aes(x=longitude, y=latitude, pch=COMM_TYPE, colour=REC_OBS_C, size=REC_OBS), alpha=1) +
			geom_text(data=dc.rec, aes(x=longitude, y=latitude, label=REC_OBS), nudge_x=0, nudge_y=0, size=4, colour='black') +			
			scale_colour_brewer(palette="YlOrRd") +
			scale_shape_manual(values=c('fishing site'=18, 'inland community'=19)) +
			scale_size(breaks=c(1,2,5,10,20,50), range=1.3*c(5,15))+						
			theme(legend.position='bottom', legend.box = "vertical") +
			guides(size='none', colour='none', pch=guide_legend(override.aes=list(size=7))) +
			labs(	x='', y='', pch="community\ntype", colour="phylogenetically observed recipients")
	ggsave(file=paste0(outfile.base,'_comm_recipientnumbers_on_map.pdf'), w=7, h=7)
	ggmap(zm) +
			geom_point(data=dc.tr, aes(x=longitude, y=latitude, pch=COMM_TYPE, colour=TR_OBS_C, size=TR_OBS), alpha=1) +
			geom_text(data=dc.tr, aes(x=longitude, y=latitude, label=TR_OBS), nudge_x=0, nudge_y=0, size=4, colour='black') +			
			scale_colour_brewer(palette="YlOrRd") +
			scale_size(breaks=c(1,2,5,10,20,50), range=1.3*c(5,15))+						
			theme(legend.position='bottom', legend.box = "vertical") +
			guides(size='none', colour='none', pch=guide_legend(override.aes=list(size=7))) +
			labs(	x='', y='', pch="community\ntype", colour="phylogenetically observed transmitters")
	ggsave(file=paste0(outfile.base,'_comm_transmitternumbers_on_map.pdf'), w=7, h=7)
	
	#
	#	plot number of observed transmission flows	
	edge.gap	<- 0.04
	node.size	<- 0.2
	edge.size	<- 0.4
	curvature	<- -0.2
	label.size	<- 5
	arrow		<- arrow(length=unit(0.04, "npc"))
	arrow2		<- arrow(length=unit(0.02, "npc"))
	curv.shift	<- 0.08
	dc		<- rtr2[, list(FLOW=length(PAIRID)), by=c('TR_COMM_NUM_A','REC_COMM_NUM_A')]
	set(dc, NULL, 'TR_COMM_NUM_A', dc[, as.character(TR_COMM_NUM_A)])
	set(dc, NULL, 'REC_COMM_NUM_A', dc[, as.character(REC_COMM_NUM_A)])
	layout	<- subset(ds, select=c(COMM_NUM_A, longitude, latitude))
	set(layout, NULL, 'COMM_NUM_A', layout[, as.character(COMM_NUM_A)])
	setnames(layout, c('COMM_NUM_A','longitude','latitude'), c('TR_COMM_NUM_A','TR_X','TR_Y'))
	dc		<- merge(dc, layout, by='TR_COMM_NUM_A')
	setnames(layout, c('TR_COMM_NUM_A','TR_X','TR_Y'), c('REC_COMM_NUM_A','REC_X','REC_Y'))
	dc		<- merge(dc, layout, by='REC_COMM_NUM_A')
	setnames(layout, c('REC_COMM_NUM_A','REC_X','REC_Y'),  c('COMM_NUM_A','longitude','latitude'))			
	dc[, EDGETEXT_X:= (TR_X+REC_X)/2]
	dc[, EDGETEXT_Y:= (TR_Y+REC_Y)/2]
	dc[, EDGE_LABEL:= FLOW ]
	#	for edges, move the start and end points on the line between X and Y
	#	define unit gradient
	dc[, MX:= (REC_X - TR_X)]	
	dc[, MY:= (REC_Y - TR_Y)]
	#	label could just be move on the tangent vector to the line
	#	define unit tangent
	dc[, TX:= -MY]
	dc[, TY:= MX]
	set(dc, NULL, 'EDGETEXT_X', dc[, EDGETEXT_X + TX*curv.shift])
	set(dc, NULL, 'EDGETEXT_Y', dc[, EDGETEXT_Y + TY*curv.shift])	
	ggmap(zm) +		
			#geom_point(data=subset(dc,  REC_COMM_NUM_A==TR_COMM_NUM_A), aes(x=TR_X, y=TR_Y, size=node.size*EDGE_LABEL)) +
			geom_curve(data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & REC_COMM_NUM_A!=TR_COMM_NUM_A & FLOW>0), aes(x=TR_X, xend=REC_X, y=TR_Y, yend=REC_Y, size=edge.size*FLOW), colour='darkorange1', curvature=curvature, arrow=arrow, lineend="butt") +
			geom_text(data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & FLOW>0), aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL), colour='darkorange4', size=label.size) +
			coord_cartesian() +
			scale_size(breaks=c(1,2,4), range=1*c(0.6,1.2)) +
			labs(x='', y='') +
			theme(legend.position='bottom', legend.box = "vertical")
	ggsave(file=paste0(outfile.base,'_flows_intofishing_on_map.pdf'), w=7, h=7)
	
	subset(dc,  substr(TR_COMM_NUM_A,1,1)=='f')
	edge.size	<- 1
	ggmap(zm) +					
			#geom_point(data=subset(dc,  substr(TR_COMM_NUM_A,1,1)=='f'), aes(x=TR_X, y=TR_Y), size=3, colour='firebrick1') +
			#geom_point(data=subset(dc,  substr(TR_COMM_NUM_A,1,1)=='f'), aes(x=REC_X, y=REC_Y), size=3, colour='firebrick1') +
			geom_curve(data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A!=TR_COMM_NUM_A & FLOW>0), aes(x=TR_X, xend=REC_X, y=TR_Y, yend=REC_Y, size=edge.size*FLOW), colour='firebrick1', curvature=curvature, arrow=arrow, lineend="butt") +
			geom_text(data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f'  & FLOW>0), aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL), colour='firebrick4', size=label.size) +
			coord_cartesian() +
			scale_size(breaks=c(1,2,4), range=1*c(0.6,1.2))+
			labs(x='', y='') +
			theme(legend.position='bottom', legend.box = "vertical")
	ggsave(file=paste0(outfile.base,'_flows_fromfishing_on_map.pdf'), w=7, h=7)
	
	ggmap(zm) +		
			#geom_point(data=subset(dc,  REC_COMM_NUM_A==TR_COMM_NUM_A), aes(x=TR_X, y=TR_Y, size=node.size*EDGE_LABEL)) +
			geom_curve(data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)!='f' & REC_COMM_NUM_A!=TR_COMM_NUM_A & FLOW>0), aes(x=TR_X, xend=REC_X, y=TR_Y, yend=REC_Y, size=edge.size*FLOW), colour='darkorange1', curvature=curvature, arrow=arrow) +			
			geom_curve(data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)!='f' & REC_COMM_NUM_A==TR_COMM_NUM_A & FLOW>0), aes(x=TR_X+0.005, xend=REC_X+0.03, y=TR_Y-0.005, yend=REC_Y, size=edge.size*FLOW), colour='darkorange1', curvature=1) +
			geom_curve(data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)!='f' & REC_COMM_NUM_A==TR_COMM_NUM_A & FLOW>0), aes(x=TR_X+0.03, xend=REC_X+0.005, y=TR_Y, yend=REC_Y+0.005, size=edge.size*FLOW), colour='darkorange1', curvature=1, arrow=arrow2) +
			geom_text(data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)!='f' & FLOW>0), aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL), colour='darkorange4', size=label.size) +
			coord_cartesian() +
			scale_size(breaks=c(1,2,4), range=1*c(0.6,1.2)) +
			labs(x='', y='') +
			theme(legend.position='bottom', legend.box = "vertical")
	ggsave(file=paste0(outfile.base,'_flows_inland_on_map.pdf'), w=7, h=7)
	ggmap(zm) +		
			#geom_point(data=subset(dc,  REC_COMM_NUM_A==TR_COMM_NUM_A), aes(x=TR_X, y=TR_Y, size=node.size*EDGE_LABEL)) +
			geom_curve(data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A!=TR_COMM_NUM_A & FLOW>0), aes(x=TR_X, xend=REC_X, y=TR_Y, yend=REC_Y, size=edge.size*FLOW), colour='firebrick1', curvature=curvature, arrow=arrow) +			
			geom_curve(data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A==TR_COMM_NUM_A & FLOW>0), aes(x=TR_X+0.005, xend=REC_X+0.03, y=TR_Y-0.005, yend=REC_Y, size=edge.size*FLOW), colour='firebrick1', curvature=1) +
			geom_curve(data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A==TR_COMM_NUM_A & FLOW>0), aes(x=TR_X+0.03, xend=REC_X+0.005, y=TR_Y, yend=REC_Y+0.005, size=edge.size*FLOW), colour='firebrick1', curvature=1, arrow=arrow2) +
			geom_text(data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)=='f' & FLOW>0), aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL), colour='firebrick4', size=label.size) +
			coord_cartesian() +
			scale_size(breaks=c(1,2,4), range=1*c(0.6,1.2)) +
			labs(x='', y='') +
			theme(legend.position='bottom', legend.box = "vertical")
	ggsave(file=paste0(outfile.base,'_flows_fishing_on_map.pdf'), w=7, h=7)
	
	ggmap(zm) +					
			geom_curve(data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)!='f' & REC_COMM_NUM_A!=TR_COMM_NUM_A & FLOW>0), aes(x=TR_X, xend=REC_X, y=TR_Y, yend=REC_Y, size=edge.size*FLOW), colour='firebrick1', curvature=curvature, arrow=arrow) +			
			geom_text(data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)!='f' & FLOW>0), aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL), colour='firebrick4', size=label.size) +
			coord_cartesian() +
			scale_size(breaks=c(1,2,4), range=1*c(0.6,1.2)) +
			labs(x='', y='') +
			theme(legend.position='bottom', legend.box = "vertical")
	ggsave(file=paste0(outfile.base,'_flows_fromfishing_on_map.pdf'), w=7, h=7)
	
	ggmap(zm) +		
			geom_curve(data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A!=TR_COMM_NUM_A & FLOW>0), aes(x=TR_X, xend=REC_X, y=TR_Y, yend=REC_Y, size=edge.size*FLOW), colour='darkorange1', curvature=curvature, arrow=arrow) +
			geom_text(data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)=='f' & FLOW>0), aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL), colour='darkorange4', size=label.size) +
			coord_cartesian() +
			scale_size(breaks=c(1,2,4), range=1*c(0.6,1.2)) +
			labs(x='', y='') +
			theme(legend.position='bottom', legend.box = "vertical")
	ggsave(file=paste0(outfile.base,'_flows_frominland_on_map.pdf'), w=7, h=7)
}


RakaiFull.phylogeography.171122.core.inference<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl35_prior23_min30_withmetadata.rda"
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_prior23_min30_withmetadata.rda"
	outfile.base			<- gsub('_withmetadata.rda','',infile)
	
	zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	load(infile)
	
	ds	<- subset(rsm, MIN_PNG_OUTPUT>0, select=c(COMM_NUM, COMM_TYPE, COMM_NUM_A, ELIGIBLE_AVG, PARTICIPATED_AVG, HIV, MIN_PNG_OUTPUT))
	ds[, P_PART_EMP:= PARTICIPATED_AVG/ELIGIBLE_AVG]
	ds[, P_PART_ALPHA:= round(PARTICIPATED_AVG)+1]
	ds[, P_PART_BETA:= round(ELIGIBLE_AVG-PARTICIPATED_AVG)+1]
	ds[, P_SEQ_EMP:= MIN_PNG_OUTPUT/HIV]
	ds[, P_SEQ_ALPHA:= round(MIN_PNG_OUTPUT)+1]
	ds[, P_SEQ_BETA:= round(HIV-MIN_PNG_OUTPUT)+1]		
	tmp	<- unique(subset(zc, select=c(COMM_NUM, longitude, latitude)), by='COMM_NUM')
	ds	<- merge(ds, tmp, by='COMM_NUM')
	set(ds, NULL, 'COMM_NUM_A', ds[, as.character(COMM_NUM_A)])
	set(ds, NULL, 'COMM_TYPE', ds[, as.character(COMM_TYPE)])
	
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	nrow(rtpdm)
	#	284	with 0.035
	#	238	with 0.025
	
	rtpdm[, AGEDIFF:= rtpdm[, FEMALE_BIRTHDATE-MALE_BIRTHDATE]]	
	set(rtpdm, NULL, 'MALE_SEX', 'M')
	set(rtpdm, NULL, 'FEMALE_SEX', 'F')	
	
	#	subset(rsm, MIN_PNG_OUTPUT>0)[, quantile(MIN_PNG_OUTPUT/HIV, p=c(0,0.5,1))]
	
	#
	#	some helper data.tables
	rmf		<- subset(rtpdm, grepl('mf',SELECT) )
	rfm		<- subset(rtpdm, grepl('fm',SELECT) )
	rtr2	<- copy(rmf)
	setnames(rtr2,colnames(rtr2),gsub('FEMALE','REC',colnames(rtr2)))
	setnames(rtr2,colnames(rtr2),gsub('MALE','TR',colnames(rtr2)))
	tmp		<- copy(rfm)
	setnames(tmp,colnames(tmp),gsub('FEMALE','TR',colnames(tmp)))
	setnames(tmp,colnames(tmp),gsub('MALE','REC',colnames(tmp)))
	rtr2	<- rbind(rtr2,tmp)
	
	#
	#	Bayesian source attriution model
	#	adjust for incomplete sampling
	#
	dc	<- rtr2[, list(TR_OBS=length(PAIRID)), by=c('TR_COMM_NUM_A','REC_COMM_NUM_A')]
	set(dc, NULL, 'TR_COMM_NUM_A', dc[, as.character(TR_COMM_NUM_A)])
	set(dc, NULL, 'REC_COMM_NUM_A', dc[, as.character(REC_COMM_NUM_A)])	
	setkey(dc, TR_COMM_NUM_A, REC_COMM_NUM_A)
	#	Bayesian model: add uniform prior
	if(0)
	{
		#	(which is Dirichlet 1 among all communities pairs that have a connection either way
		tmp	<- subset(dc, select=c(REC_COMM_NUM_A, TR_COMM_NUM_A))	
		setnames(tmp, c('REC_COMM_NUM_A','TR_COMM_NUM_A'), c('TR_COMM_NUM_A','REC_COMM_NUM_A'))
		tmp	<- merge(tmp, dc, all.x=1)
		tmp	<- subset(tmp, is.na(TR_OBS))
		set(tmp, NULL, 'TR_OBS', 0)
		dc	<- rbind(dc, tmp)		
	}
	if(1)
	{
		#	This is a bit non-standard, I just don t want the prior to have a large impact, so I chose a sparse one. 
		#	(which is Dirichlet 1 among all communities pairs that are closest)
		#	always add self if not present
		tmp	<- subset(dc[, list(UNOBSERVED_SELF=!any(REC_COMM_NUM_A==TR_COMM_NUM_A)), by='TR_COMM_NUM_A'],UNOBSERVED_SELF, TR_COMM_NUM_A)
		tmp[, REC_COMM_NUM_A:=TR_COMM_NUM_A]
		tmp[, TR_OBS:=0]
		dc	<- rbind(dc, tmp)
		#	there are ZERO transmissions from other commtype to trading, add one
		#dc	<- rbind(dc, data.table(REC_COMM_NUM_A=c('thd','ttp','tpq','tpq'), TR_COMM_NUM_A=c('ald','fpt','fpt','adi'), TR_OBS=0))
		dc	<- rbind(dc, data.table(REC_COMM_NUM_A=c('thd','ttp','tpq'), TR_COMM_NUM_A=c('ald','fpt','adi'), TR_OBS=1))
		#	ensure each community has at least 1 non-self community, if not add closest other community
		#	I really want to keep this sparse, so do not consider non-self to all communities
		tmp	<- unique(ds, by='COMM_NUM_A')
		clo	<- as.data.table(t(sapply(seq_len(nrow(tmp)), function(i)
								{
									z<- sort( sqrt( (tmp[,longitude]-tmp[i,longitude])^2+(tmp[,latitude]-tmp[i,latitude])^2 ), index.return=TRUE)$ix
									c('TR_COMM_NUM_A'=tmp[i, COMM_NUM_A], 'REC_COMM_NUM_A'=tmp[z[2],COMM_NUM_A])
								})))	
		z	<- subset(dc[, list(REC_N_OBS=length(REC_COMM_NUM_A)), by='TR_COMM_NUM_A'], REC_N_OBS==1)
		tmp	<- merge(clo, z, by='TR_COMM_NUM_A')
		tmp	<- merge(subset(tmp, select=c(TR_COMM_NUM_A, REC_COMM_NUM_A)), dc, all.x=1, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A'))
		tmp	<- subset(tmp, is.na(TR_OBS))
		#z	<- copy(tmp)
		#setnames(z, c('REC_COMM_NUM_A','TR_COMM_NUM_A'), c('TR_COMM_NUM_A','REC_COMM_NUM_A'))
		#tmp	<- rbind(tmp,z)
		set(tmp, NULL, 'TR_OBS', 0)
		dc	<- rbind(dc, tmp)
		tmp	<- unique(ds, by='COMM_NUM_A')
		clo	<- as.data.table(t(sapply(seq_len(nrow(tmp)), function(i)
								{
									z<- sort( sqrt( (tmp[,longitude]-tmp[i,longitude])^2+(tmp[,latitude]-tmp[i,latitude])^2 ), index.return=TRUE)$ix
									c('REC_COMM_NUM_A'=tmp[i, COMM_NUM_A], 'TR_COMM_NUM_A'=tmp[z[2],COMM_NUM_A])
								})))			
		z	<- subset(dc[, list(TR_N_OBS=length(TR_COMM_NUM_A)), by='REC_COMM_NUM_A'], TR_N_OBS==1)
		tmp	<- merge(clo, z, by='REC_COMM_NUM_A')
		tmp	<- merge(subset(tmp, select=c(TR_COMM_NUM_A, REC_COMM_NUM_A)), dc, all.x=1, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A'))
		tmp	<- subset(tmp, is.na(TR_OBS))
		#z	<- copy(tmp)
		#setnames(z, c('REC_COMM_NUM_A','TR_COMM_NUM_A'), c('TR_COMM_NUM_A','REC_COMM_NUM_A'))
		#tmp	<- rbind(tmp,z)		
		set(tmp, NULL, 'TR_OBS', 0)
		dc	<- rbind(dc, tmp)		
		
	}
	#rsm[, list(ELIGIBLE_AVG=sum(ELIGIBLE_AVG)), by='COMM_TYPE']
	dc[, TR_PRIOR:= 0.5]
	#
	#	Bayesian model first hierarchy: define Beta posterior for sampling probabilities (all alpha and betas)
	#
	tmp	<- copy(ds)
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc	<- merge(dc, tmp, by='TR_COMM_NUM_A')
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc	<- merge(dc, tmp, by='REC_COMM_NUM_A')
	#
	#	Bayesian model second hierarchy: draw unobserved data to augment likelihood
	#
	mc.it	<- 1e4
	dcb		<- dc[, {
				tmp	<- 	rbeta(mc.it, TR_P_PART_ALPHA, TR_P_PART_BETA)*
						rbeta(mc.it, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA)*
						rbeta(mc.it, REC_P_PART_ALPHA, REC_P_PART_BETA)*
						rbeta(mc.it, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA)
				#print(tmp)
				tmp	<- rnbinom(mc.it, TR_OBS+TR_PRIOR, tmp)
				#print(tmp)
				list(MONTE_CARLO_IT=seq_len(mc.it), TR_PRIOR=TR_PRIOR, TR_OBS=TR_OBS, TR_MISS= tmp)
			}, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A')]	
	#
	#	Bayesian model second hierarchy: Dirichlet posterior for transmission from community i to j, pi_ij with pi_ij summing to 1
	#
	tmp		<- dcb[, list(	REC_COMM_NUM_A= REC_COMM_NUM_A, 
					TR_COMM_NUM_A= TR_COMM_NUM_A, 
					PI_IJ_ALPHA= TR_OBS+TR_MISS+TR_PRIOR				
			), by='MONTE_CARLO_IT']
	dcb		<- merge(dcb, tmp, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A','MONTE_CARLO_IT'))
	#
	#	this is the end of the source attribution inference on the WAIFM matrix
	#
	save(zm, zc, rtpdm, rfm, rmf, rtr2, ds, dc, dcb, file=paste0(outfile.base,'_inference.rda'))	
}

RakaiFull.phylogeography.180322.core.inference.gender<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	#	load des which contains participation and seq counts by comm and gender
	infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender.rda'
	load(infile)
	#	determine posterior parameters for Binomial models of sampling and participiation 
	ds		<- subset(des, DEEP_SEQ_EVER>0)
	set(ds, NULL, c('HIV_1516_NO','HIV_EVER_YES','DEEP_SEQ_EVER'), NULL)
	ds[, P_PART_EMP:= PART_EVER/(PART_EVER+PART_NEVER)]
	ds[, P_PART_ALPHA:= PART_EVER+1]
	ds[, P_PART_BETA:= PART_NEVER+1]
	ds[, P_SEQ_EMP:= DEEP_SEQ_1516/HIV_1516_YES]
	ds[, P_SEQ_ALPHA:= DEEP_SEQ_1516+1]
	ds[, P_SEQ_BETA:= HIV_1516_YES+1]
	#	add long lat
	zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	tmp		<- unique(subset(zc, select=c(COMM_NUM, longitude, latitude)), by='COMM_NUM')
	setnames(tmp, c('longitude','latitude'),c('LONG','LAT'))
	ds		<- merge(ds, tmp, by='COMM_NUM')
	set(ds, NULL, 'COMM_NUM_A', ds[, as.character(COMM_NUM_A)])
	set(ds, NULL, 'COMM_TYPE', ds[, as.character(COMM_TYPE)])
	
	#	get map
	style	<- "feature:road|color:0x17202A&style=feature:water|color:0x2874A6&style=feature:administrative|visibility=off"
	zm		<- get_googlemap(c(lon=31.65, lat=-0.66), scale=2, size=c(550,550), zoom=10, maptype="road", style=style)	
			
	#	load transmission events
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_withmetadata.rda"
	outfile.base			<- gsub('_withmetadata.rda','',infile)
	load(infile)
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	nrow(rtpdm)
	#	293
	
	#	add new variables
	rtpdm[, AGEDIFF:= rtpdm[, FEMALE_BIRTHDATE-MALE_BIRTHDATE]]	
	set(rtpdm, NULL, 'MALE_SEX', 'M')
	set(rtpdm, NULL, 'FEMALE_SEX', 'F')			
	
	#	cast MALE FEMALE to TRM REC
	rmf		<- subset(rtpdm, grepl('mf',SELECT) )
	rfm		<- subset(rtpdm, grepl('fm',SELECT) )
	rtr2	<- copy(rmf)
	setnames(rtr2,colnames(rtr2),gsub('FEMALE','REC',colnames(rtr2)))
	setnames(rtr2,colnames(rtr2),gsub('MALE','TR',colnames(rtr2)))
	tmp		<- copy(rfm)
	setnames(tmp,colnames(tmp),gsub('FEMALE','TR',colnames(tmp)))
	setnames(tmp,colnames(tmp),gsub('MALE','REC',colnames(tmp)))
	rtr2	<- rbind(rtr2,tmp)
	
	#	sum transmission events by community and gender
	dc	<- rtr2[, list(TR_OBS=length(PAIRID)), by=c('TR_COMM_NUM_A','REC_COMM_NUM_A','TR_SEX','REC_SEX')]
	set(dc, NULL, 'TR_COMM_NUM_A', dc[, as.character(TR_COMM_NUM_A)])
	set(dc, NULL, 'REC_COMM_NUM_A', dc[, as.character(REC_COMM_NUM_A)])	
	setkey(dc, TR_COMM_NUM_A, REC_COMM_NUM_A, TR_SEX)
	#	ensure we have always have same community flows
	tmp	<- as.data.table(expand.grid(TR_COMM_NUM_A=unique(c(dc$TR_COMM_NUM_A, dc$REC_COMM_NUM_A)), TR_SEX=c('M','F')))
	tmp[, REC_COMM_NUM_A:= TR_COMM_NUM_A]
	tmp[, REC_SEX:= as.character(factor(TR_SEX=='M', levels=c(TRUE,FALSE),labels=c('F','M')))]
	dc	<- merge(tmp, dc, by=c('TR_COMM_NUM_A','REC_COMM_NUM_A','TR_SEX','REC_SEX'), all=TRUE)
	set(dc, dc[, which(is.na(TR_OBS))], 'TR_OBS', 0L)
	
	#	Bayesian model: add all connections and set to 0 (for uniform prior)
	if(0)
	{
		#	(which is Dirichlet 1 among all communities pairs that have a connection either way
		tmp	<- subset(dc, select=c(REC_COMM_NUM_A, TR_COMM_NUM_A))	
		setnames(tmp, c('REC_COMM_NUM_A','TR_COMM_NUM_A'), c('TR_COMM_NUM_A','REC_COMM_NUM_A'))
		tmp	<- merge(tmp, dc, all.x=1)
		tmp	<- subset(tmp, is.na(TR_OBS))
		set(tmp, NULL, 'TR_OBS', 0)
		dc	<- rbind(dc, tmp)		
	}
	#	Bayesian model: add "sparse" connections and set to 0 (for "sparse" prior)
	if(1)
	{
		#	This is non-standard, I just don t want the prior to have a large impact, so I chose a sparse one. 
		#	(which is Dirichlet 1 among all communities pairs that are closest)	
		#	ensure each community has at least 1 non-self community, if not add closest other community
		#	I really want to keep this sparse, so do not consider non-self to all communities
		tmp	<- unique(ds, by='COMM_NUM_A')
		#	find geographically closest community
		clo	<- as.data.table(t(sapply(seq_len(nrow(tmp)), function(i)
								{
									z<- sort( sqrt( (tmp[,LONG]-tmp[i,LONG])^2+(tmp[,LAT]-tmp[i,LAT])^2 ), index.return=TRUE)$ix
									c('TR_COMM_NUM_A'=tmp[i, COMM_NUM_A], 'REC_COMM_NUM_A'=tmp[z[2],COMM_NUM_A])
								})))	
		#	find transmitter communities with no non-self recipient
		z	<- subset(dc[, list(REC_N_OBS=length(which(REC_COMM_NUM_A!=TR_COMM_NUM_A))), by='TR_COMM_NUM_A'], REC_N_OBS==0)
		z	<- merge(z, data.table(REC_N_OBS=0, TR_SEX=c('M','F')), by='REC_N_OBS',allow.cartesian=TRUE)
		set(z, NULL, 'REC_N_OBS', NULL)
		z	<- merge(clo, z, by='TR_COMM_NUM_A')
		z[, REC_SEX:= as.character(factor(TR_SEX=='M', levels=c(TRUE,FALSE),labels=c('F','M')))]		
		dc	<- rbind(dc, z, fill=TRUE)
		#	find recipient communities with no non-self recipient
		tmp	<- unique(ds, by='COMM_NUM_A')
		clo	<- as.data.table(t(sapply(seq_len(nrow(tmp)), function(i)
								{
									z<- sort( sqrt( (tmp[,LONG]-tmp[i,LONG])^2+(tmp[,LAT]-tmp[i,LAT])^2 ), index.return=TRUE)$ix
									c('REC_COMM_NUM_A'=tmp[i, COMM_NUM_A], 'TR_COMM_NUM_A'=tmp[z[2],COMM_NUM_A])
								})))			
		z	<- subset(dc[, list(TR_N_OBS=length(which(TR_COMM_NUM_A!=REC_COMM_NUM_A))), by='REC_COMM_NUM_A'], TR_N_OBS==0)
		z	<- merge(z, data.table(TR_N_OBS=0, REC_SEX=c('M','F')), by='TR_N_OBS',allow.cartesian=TRUE)
		set(z, NULL, 'TR_N_OBS', NULL)
		z	<- merge(clo, z, by='REC_COMM_NUM_A')
		z[, TR_SEX:= as.character(factor(REC_SEX=='M', levels=c(TRUE,FALSE),labels=c('F','M')))]
		dc	<- rbind(dc, z, fill=TRUE)
		#	0 observations in each of these connections
		set(dc, dc[, which(is.na(TR_OBS))],'TR_OBS',0L)		
	}
	
	# 	add prior observations
	dc[, TR_PRIOR:= 0.5]
	
	#
	#	Bayesian model first hierarchy: define Beta posterior for sampling probabilities (all alpha and betas)
	#
	tmp	<- subset(ds, select=c(COMM_NUM_A, SEX, P_PART_ALPHA, P_PART_BETA, P_SEQ_ALPHA, P_SEQ_BETA))
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc	<- merge(dc, tmp, by=c('TR_COMM_NUM_A','TR_SEX'))
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc	<- merge(dc, tmp, by=c('REC_COMM_NUM_A','REC_SEX'))
	#
	#	Bayesian model second hierarchy: draw unobserved data to augment likelihood
	#
	mc.it	<- 5e4
	dcb		<- dc[, {
				#print( c(TR_P_PART_ALPHA, TR_P_PART_BETA, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA, REC_P_PART_ALPHA, REC_P_PART_BETA, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA) )
				tmp	<- 	rbeta(mc.it, TR_P_PART_ALPHA, TR_P_PART_BETA)*
						rbeta(mc.it, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA)*
						rbeta(mc.it, REC_P_PART_ALPHA, REC_P_PART_BETA)*
						rbeta(mc.it, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA)
				#print(mean(tmp))
				#print(TR_P_PART_ALPHA/(TR_P_PART_ALPHA+TR_P_PART_BETA) * TR_P_SEQ_ALPHA/(TR_P_SEQ_ALPHA+TR_P_SEQ_BETA) * REC_P_PART_ALPHA/(REC_P_PART_ALPHA+REC_P_PART_BETA) * REC_P_SEQ_ALPHA/(REC_P_SEQ_ALPHA+REC_P_SEQ_BETA) )
				#print(tmp)
				tmp	<- rnbinom(mc.it, TR_OBS+TR_PRIOR, tmp)
				#print(tmp)
				#stop()
				list(	MONTE_CARLO_IT=seq_len(mc.it), 
						TR_PRIOR=TR_PRIOR, 
						TR_OBS=TR_OBS, 
						TR_MISS= tmp,
						TR_MISS_P= TR_P_PART_ALPHA/(TR_P_PART_ALPHA+TR_P_PART_BETA) * TR_P_SEQ_ALPHA/(TR_P_SEQ_ALPHA+TR_P_SEQ_BETA) * REC_P_PART_ALPHA/(REC_P_PART_ALPHA+REC_P_PART_BETA) * REC_P_SEQ_ALPHA/(REC_P_SEQ_ALPHA+REC_P_SEQ_BETA)   )
			}, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A','REC_SEX','TR_SEX')]
	#	check
	#tmp	<- dcb[, list(	AUG_E= (TR_PRIOR[1]+TR_OBS[1])/TR_MISS_P[1],
	#					AUG_MEAN= mean(TR_MISS)+TR_PRIOR[1]+TR_OBS[1]
	#					), by=c('REC_COMM_NUM_A','TR_COMM_NUM_A','REC_SEX','TR_SEX')]
	#tmp[, summary(AUG_E-AUG_MEAN)]
	#
	#	Bayesian model second hierarchy: Dirichlet posterior for transmission from community i to j, pi_ij with pi_ij summing to 1
	#
	tmp		<- dcb[, list(	REC_COMM_NUM_A= REC_COMM_NUM_A, 
							TR_COMM_NUM_A= TR_COMM_NUM_A, 
							REC_SEX=REC_SEX,
							TR_SEX=TR_SEX,
							PI_IJ_ALPHA= TR_OBS+TR_MISS+TR_PRIOR				
					), by='MONTE_CARLO_IT']
	dcb		<- merge(dcb, tmp, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A','REC_SEX','TR_SEX','MONTE_CARLO_IT'))
	#
	#	this is the end of the source attribution inference on the WAIFM matrix
	#
	save(zc, rtpdm, rfm, rmf, rtr2, ds, dc, dcb, file=paste0(outfile.base,'_phylogeography_core_inference.rda'))	
}

RakaiFull.phylogeography.180322.core.inference.age<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	#	load des which contains participation and seq counts by comm and gender
	infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender.rda'
	load(infile)
	#	determine posterior parameters for Binomial models of sampling and participiation
	tmp		<- melt(desa,id.vars=c('COMM_TYPE','COMM_NUM','COMM_NUM_A','AGE_AT_MID_C','SEX'))
	tmp		<- tmp[, list(value=sum(value)), by=c('COMM_TYPE','COMM_NUM','COMM_NUM_A','AGE_AT_MID_C','variable')]
	ds		<- dcast.data.table(tmp, COMM_TYPE+COMM_NUM+COMM_NUM_A+AGE_AT_MID_C~variable, value.var='value')	
	ds		<- subset(ds, DEEP_SEQ_1516>0)

	ds[, P_PART_EMP:= PART_EVER/(PART_EVER+PART_NEVER)]
	ds[, P_PART_ALPHA:= PART_EVER+1]
	ds[, P_PART_BETA:= PART_NEVER+1]
	ds[, P_SEQ_EMP:= DEEP_SEQ_1516/HIV_1516_YES]
	ds[, P_SEQ_ALPHA:= DEEP_SEQ_1516+1]
	ds[, P_SEQ_BETA:= HIV_1516_YES+1]
	
	ggplot(ds, aes(x=COMM_NUM_A, colour=AGE_AT_MID_C, y=P_SEQ_EMP)) + geom_point()
	ggsave(file=paste0(outfile.base,'_phylogeographyage_seqdifferences.pdf'), w=12, h=6)
	ggplot(ds, aes(x=COMM_NUM_A, colour=AGE_AT_MID_C, y=P_PART_EMP)) + geom_point()
	ggsave(file=paste0(outfile.base,'_phylogeographyage_partdifferences.pdf'), w=12, h=6)
	
	#	add long lat
	zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	tmp		<- unique(subset(zc, select=c(COMM_NUM, longitude, latitude)), by='COMM_NUM')
	setnames(tmp, c('longitude','latitude'),c('LONG','LAT'))
	ds		<- merge(ds, tmp, by='COMM_NUM')
	set(ds, NULL, 'COMM_NUM_A', ds[, as.character(COMM_NUM_A)])
	set(ds, NULL, 'COMM_TYPE', ds[, as.character(COMM_TYPE)])
	
	#	get map
	style	<- "feature:road|color:0x17202A&style=feature:water|color:0x2874A6&style=feature:administrative|visibility=off"
	zm		<- get_googlemap(c(lon=31.65, lat=-0.66), scale=2, size=c(550,550), zoom=10, maptype="road", style=style)	
	
	#	load transmission events
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_withmetadata.rda"
	outfile.base			<- gsub('_withmetadata.rda','',infile)
	load(infile)
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	nrow(rtpdm)
	#	293
	
	#	add new variables
	rtpdm[, AGEDIFF:= rtpdm[, FEMALE_BIRTHDATE-MALE_BIRTHDATE]]
	rtpdm[, MALE_AGE_AT_MID_C:= cut(2013.25-MALE_BIRTHDATE, breaks=c(10,25,35,65), labels=c('15-24','25-34','35+'), right=FALSE)]
	rtpdm[, FEMALE_AGE_AT_MID_C:= cut(2013.25-FEMALE_BIRTHDATE, breaks=c(10,25,35,65), labels=c('15-24','25-34','35+'), right=FALSE)]
	set(rtpdm, NULL, 'MALE_SEX', 'M')
	set(rtpdm, NULL, 'FEMALE_SEX', 'F')			
	
	#	cast MALE FEMALE to TRM REC
	rmf		<- subset(rtpdm, grepl('mf',SELECT) )
	rfm		<- subset(rtpdm, grepl('fm',SELECT) )
	rtr2	<- copy(rmf)
	setnames(rtr2,colnames(rtr2),gsub('FEMALE','REC',colnames(rtr2)))
	setnames(rtr2,colnames(rtr2),gsub('MALE','TR',colnames(rtr2)))
	tmp		<- copy(rfm)
	setnames(tmp,colnames(tmp),gsub('FEMALE','TR',colnames(tmp)))
	setnames(tmp,colnames(tmp),gsub('MALE','REC',colnames(tmp)))
	rtr2	<- rbind(rtr2,tmp)
	
	#	sum transmission events by community and age
	dc	<- rtr2[, list(TR_OBS=length(PAIRID)), by=c('TR_COMM_NUM_A','REC_COMM_NUM_A','TR_AGE_AT_MID_C','REC_AGE_AT_MID_C')]
	set(dc, NULL, 'TR_COMM_NUM_A', dc[, as.character(TR_COMM_NUM_A)])
	set(dc, NULL, 'REC_COMM_NUM_A', dc[, as.character(REC_COMM_NUM_A)])	
	setkey(dc, TR_COMM_NUM_A, REC_COMM_NUM_A, TR_AGE_AT_MID_C)
	#	ensure we have always have same community flows
	tmp	<- as.data.table(expand.grid(TR_COMM_NUM_A=unique(c(dc$TR_COMM_NUM_A, dc$REC_COMM_NUM_A)), TR_AGE_AT_MID_C=c('15-24','25-34','35+'), REC_AGE_AT_MID_C=c('15-24','25-34','35+') ))
	tmp[, REC_COMM_NUM_A:= TR_COMM_NUM_A]	
	dc	<- merge(tmp, dc, by=c('TR_COMM_NUM_A','REC_COMM_NUM_A','TR_AGE_AT_MID_C','REC_AGE_AT_MID_C'), all=TRUE)
	set(dc, dc[, which(is.na(TR_OBS))], 'TR_OBS', 0L)
	#	ensure for each TR, REC combination we always have all age combinations
	tmp	<- unique(subset(dc, select=c(TR_COMM_NUM_A, REC_COMM_NUM_A)))
	tmp[, DUMMY:=1]
	tmp2	<- as.data.table(expand.grid(TR_AGE_AT_MID_C=c('15-24','25-34','35+'), REC_AGE_AT_MID_C=c('15-24','25-34','35+')))
	tmp2[, DUMMY:=1]		
	tmp	<- merge(tmp, tmp2, by='DUMMY',allow.cartesian=TRUE)
	dc	<- merge(tmp, dc, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A','TR_AGE_AT_MID_C','REC_AGE_AT_MID_C'), all.x=TRUE)
	set(dc, dc[, which(is.na(TR_OBS))], 'TR_OBS', 0L)
	dc[, DUMMY:=NULL]
	#	Bayesian model: add all connections and set to 0 (for uniform prior)
	if(0)
	{
		#	(which is Dirichlet 1 among all communities pairs that have a connection either way
		tmp	<- subset(dc, select=c(REC_COMM_NUM_A, TR_COMM_NUM_A))	
		setnames(tmp, c('REC_COMM_NUM_A','TR_COMM_NUM_A'), c('TR_COMM_NUM_A','REC_COMM_NUM_A'))
		tmp	<- merge(tmp, dc, all.x=1)
		tmp	<- subset(tmp, is.na(TR_OBS))
		set(tmp, NULL, 'TR_OBS', 0)
		dc	<- rbind(dc, tmp)		
	}
	#	Bayesian model: add "sparse" connections and set to 0 (for "sparse" prior)
	if(1)
	{
		#	This is non-standard, I just don t want the prior to have a large impact, so I chose a sparse one. 
		#	(which is Dirichlet 1 among all communities pairs that are closest)	
		#	ensure each community has at least 1 non-self community, if not add closest other community
		#	I really want to keep this sparse, so do not consider non-self to all communities
		tmp	<- unique(ds, by='COMM_NUM_A')
		#	find geographically closest community
		clo	<- as.data.table(t(sapply(seq_len(nrow(tmp)), function(i)
								{
									z<- sort( sqrt( (tmp[,LONG]-tmp[i,LONG])^2+(tmp[,LAT]-tmp[i,LAT])^2 ), index.return=TRUE)$ix
									c('TR_COMM_NUM_A'=tmp[i, COMM_NUM_A], 'REC_COMM_NUM_A'=tmp[z[2],COMM_NUM_A])
								})))	
		#	find transmitter communities with no non-self recipient
		z	<- subset(dc[, list(REC_N_OBS=length(which(REC_COMM_NUM_A!=TR_COMM_NUM_A))), by='TR_COMM_NUM_A'], REC_N_OBS==0)		
		tmp	<- as.data.table(expand.grid(TR_AGE_AT_MID_C=c('15-24','25-34','35+'), REC_AGE_AT_MID_C=c('15-24','25-34','35+')))
		tmp[, REC_N_OBS:=0]		
		z	<- merge(z, tmp, by='REC_N_OBS',allow.cartesian=TRUE)
		set(z, NULL, 'REC_N_OBS', NULL)
		z	<- merge(clo, z, by='TR_COMM_NUM_A')				
		dc	<- rbind(dc, z, fill=TRUE)
		#	find recipient communities with no non-self recipient
		tmp	<- unique(ds, by='COMM_NUM_A')
		clo	<- as.data.table(t(sapply(seq_len(nrow(tmp)), function(i)
								{
									z<- sort( sqrt( (tmp[,LONG]-tmp[i,LONG])^2+(tmp[,LAT]-tmp[i,LAT])^2 ), index.return=TRUE)$ix
									c('REC_COMM_NUM_A'=tmp[i, COMM_NUM_A], 'TR_COMM_NUM_A'=tmp[z[2],COMM_NUM_A])
								})))			
		z	<- subset(dc[, list(TR_N_OBS=length(which(TR_COMM_NUM_A!=REC_COMM_NUM_A))), by='REC_COMM_NUM_A'], TR_N_OBS==0)
		tmp	<- as.data.table(expand.grid(TR_AGE_AT_MID_C=c('15-24','25-34','35+'), REC_AGE_AT_MID_C=c('15-24','25-34','35+')))
		tmp[, TR_N_OBS:=0]		
		z	<- merge(z, tmp, by='TR_N_OBS',allow.cartesian=TRUE)
		set(z, NULL, 'TR_N_OBS', NULL)
		z	<- merge(clo, z, by='REC_COMM_NUM_A')		
		dc	<- rbind(dc, z, fill=TRUE)
		#	0 observations in each of these connections
		set(dc, dc[, which(is.na(TR_OBS))],'TR_OBS',0L)		
	}
	
	# 	add prior observations
	dc[, TR_PRIOR:= 0.5*1/3]
	
	#
	#	Bayesian model first hierarchy: define Beta posterior for sampling probabilities (all alpha and betas)
	#
	tmp	<- subset(ds, select=c(COMM_NUM_A, AGE_AT_MID_C, P_PART_ALPHA, P_PART_BETA, P_SEQ_ALPHA, P_SEQ_BETA))
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc	<- merge(dc, tmp, by=c('TR_COMM_NUM_A','TR_AGE_AT_MID_C'))
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc	<- merge(dc, tmp, by=c('REC_COMM_NUM_A','REC_AGE_AT_MID_C'))
	#
	#	Bayesian model second hierarchy: draw unobserved data to augment likelihood
	#
	mc.it	<- 1e4
	dcb		<- dc[, {
				#print( c(TR_P_PART_ALPHA, TR_P_PART_BETA, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA, REC_P_PART_ALPHA, REC_P_PART_BETA, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA) )
				tmp	<- 	rbeta(mc.it, TR_P_PART_ALPHA, TR_P_PART_BETA)*
						rbeta(mc.it, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA)*
						rbeta(mc.it, REC_P_PART_ALPHA, REC_P_PART_BETA)*
						rbeta(mc.it, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA)
				#print(mean(tmp))
				#print(TR_P_PART_ALPHA/(TR_P_PART_ALPHA+TR_P_PART_BETA) * TR_P_SEQ_ALPHA/(TR_P_SEQ_ALPHA+TR_P_SEQ_BETA) * REC_P_PART_ALPHA/(REC_P_PART_ALPHA+REC_P_PART_BETA) * REC_P_SEQ_ALPHA/(REC_P_SEQ_ALPHA+REC_P_SEQ_BETA) )
				#print(tmp)
				tmp	<- rnbinom(mc.it, TR_OBS+TR_PRIOR, tmp)
				#print(tmp)
				#stop()
				list(	MONTE_CARLO_IT=seq_len(mc.it), 
						TR_PRIOR=TR_PRIOR, 
						TR_OBS=TR_OBS, 
						TR_MISS= tmp,
						TR_MISS_P= TR_P_PART_ALPHA/(TR_P_PART_ALPHA+TR_P_PART_BETA) * TR_P_SEQ_ALPHA/(TR_P_SEQ_ALPHA+TR_P_SEQ_BETA) * REC_P_PART_ALPHA/(REC_P_PART_ALPHA+REC_P_PART_BETA) * REC_P_SEQ_ALPHA/(REC_P_SEQ_ALPHA+REC_P_SEQ_BETA)   )
			}, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A','REC_AGE_AT_MID_C','TR_AGE_AT_MID_C')]
	#	check
	#tmp	<- dcb[, list(	AUG_E= (TR_PRIOR[1]+TR_OBS[1])/TR_MISS_P[1],
	#					AUG_MEAN= mean(TR_MISS)+TR_PRIOR[1]+TR_OBS[1]
	#					), by=c('REC_COMM_NUM_A','TR_COMM_NUM_A','REC_SEX','TR_SEX')]
	#tmp[, summary(AUG_E-AUG_MEAN)]
	#
	#	Bayesian model second hierarchy: Dirichlet posterior for transmission from community i to j, pi_ij with pi_ij summing to 1
	#
	tmp		<- dcb[, list(	REC_COMM_NUM_A= REC_COMM_NUM_A, 
					TR_COMM_NUM_A= TR_COMM_NUM_A, 
					REC_AGE_AT_MID_C=REC_AGE_AT_MID_C,
					TR_AGE_AT_MID_C=TR_AGE_AT_MID_C,
					PI_IJ_ALPHA= TR_OBS+TR_MISS+TR_PRIOR				
			), by='MONTE_CARLO_IT']
	dcb		<- merge(dcb, tmp, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A','REC_AGE_AT_MID_C','TR_AGE_AT_MID_C','MONTE_CARLO_IT'))
	#
	#	this is the end of the source attribution inference on the WAIFM matrix
	#
	save(zc, zm, rtpdm, rfm, rmf, rtr2, ds, dc, dcb, file=paste0(outfile.base,'_phylogeographyage_core_inference.rda'))	
}


RakaiFull.phylogeography.171122.transmission.hubs<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_prior23_min30_inference.rda"
	outfile.base		<- gsub('_inference.rda','',infile.inference)
	load(infile.inference)
	
	#
	#	transmission hubs crude estimates		
	z	<- copy(rtr2)
	z	<- z[, list(REC_OUTSIDE_COMM=length(which(REC_COMM_NUM_A!=TR_COMM_NUM_A)), REC_IN_COMM=length(which(REC_COMM_NUM_A==TR_COMM_NUM_A)) ), by=c('TR_COMM_NUM_A','TR_COMM_TYPE')]
	z[, REC_N:=REC_OUTSIDE_COMM+REC_IN_COMM]
	setnames(z, 'TR_COMM_NUM_A', 'COMM_NUM_A')
	z	<- merge(unique(zc, by='COMM_NUM_A'), z, by='COMM_NUM_A')
	ggmap(zm) +
			geom_point(data=z, aes(x=longitude, y=latitude, pch=COMM_TYPE, size=REC_N, fill=100*REC_OUTSIDE_COMM/REC_N), stroke=1.5, alpha=0.8) + 
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +
			scale_size(breaks=c(5,10,20,40,80), range=c(5,20))+
			scale_shape_manual(values=c('agrarian'=21, 'trading'=22, 'fisherfolk'=23)) +
			scale_colour_brewer(palette='Dark2') +			
			theme(legend.position='bottom', legend.box = "vertical") +
			guides(	pch=guide_legend(override.aes=list(size=7))) +
			labs(	x='', y='', 
					size="phylogenetically inferred transmissions\nwith transmitters from community",
					pch="community type",
					fill="recipients outside community")
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_crude.pdf'), w=7, h=7)
	
	
	#
	#	transmission hubs
	#	proportion of transmitters in community i that have a recipient outside community
	#	this is a summary on the estimated posterior density of the WAIFM matrix
	#	TODO missing community 94
	#
	#	get parameters of posterior under augmented likelihood
	z		<- dcb[, {
				z<- which(REC_COMM_NUM_A==TR_COMM_NUM_A)
				list(P_RECOUTSIDE_ALPHA= sum(PI_IJ_ALPHA[-z]), P_RECOUTSIDE_BETA= sum(PI_IJ_ALPHA[z]))
			}, by=c('TR_COMM_NUM_A','MONTE_CARLO_IT')]	
	#	aggregate and get quantiles
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- quantile(rbeta(length(P_RECOUTSIDE_ALPHA)*mc.it, P_RECOUTSIDE_ALPHA, P_RECOUTSIDE_BETA), p=seq(0,1,0.01))
				list(P=seq(0,1,0.01), Q=unname(tmp))
			}, by='TR_COMM_NUM_A']
	#	subset to main quantities of interest
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_NUM_A~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('REC_P_OUTSIDE_CL','REC_P_OUTSIDE_IL','REC_P_OUTSIDE_M','REC_P_OUTSIDE_IU','REC_P_OUTSIDE_CU'))
	#	add TR_OBS from TR_COMM_NUM_A
	z		<- merge(z, dc[, list(REC_N_OBS=sum(TR_OBS)), by='TR_COMM_NUM_A'], by='TR_COMM_NUM_A')
	setnames(z, 'TR_COMM_NUM_A', 'COMM_NUM_A')	
	#	now after analysis, remove 94 with unknown long / lat
	#z		<- subset(z, COMM_NUM!='94')
	z	<- merge(z, unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE))), by='COMM_NUM_A')
	ggplot(z, aes(x=COMM_NUM_A, middle=REC_P_OUTSIDE_M, min=REC_P_OUTSIDE_CL, max=REC_P_OUTSIDE_CU, lower=REC_P_OUTSIDE_IL, upper=REC_P_OUTSIDE_IU, fill=REC_P_OUTSIDE_M)) + 
			geom_boxplot(outlier.shape=NA, stat='identity') +
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +
			#geom_point(data=subset(z, BS==0), pch=2, colour='red') +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.2)) +
			facet_grid(~COMM_TYPE, scales='free_x', space='free_x') +
			labs(x='\ncommunity', y="recipients outside community") +
			theme_bw() + guides(fill='none')
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_adjusted_confidenceintervals_anonmyized.pdf'), w=9, h=5)		
	
	#
	#	transmission hubs - central estimates
	#	among transmitters how many have a recipient outside community
	#	adjusted for sampling
	#	TODO missing community 94
	z[, REC_P_OUTSIDE_PREC:= 1 - REC_P_OUTSIDE_IU + REC_P_OUTSIDE_IL]	
	z	<- merge(unique(subset(zc, select=c(COMM_NUM_A, longitude, latitude, LONG_A, LAT_A)), by='COMM_NUM_A'), z, by='COMM_NUM_A')	
	ggmap(zm) +
			geom_point(data=z, aes(x=longitude, y=latitude, pch=COMM_TYPE, size=REC_N_OBS, fill=100*REC_P_OUTSIDE_M), stroke=1.5, alpha=0.8) +
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +			
			scale_size(breaks=c(5,10,20,40,80), range=c(5,20))+
			scale_shape_manual(values=c('agrarian'=21, 'trading'=22, 'fisherfolk'=23)) +
			theme(	legend.position='bottom', legend.box = "vertical") +	
			guides(	pch=guide_legend(override.aes=list(size=7))) +
			labs(	x='', y='', size="phylogenetically inferred transmissions\nwith transmitters from community",
					pch="community type",
					fill="proportion of transmitters from community\nwho have recipient outside community")
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_adjusted.pdf'), w=7, h=7)
	ggmap(zm) +
			geom_point(data=z, aes(x=LONG_A, y=LAT_A, pch=COMM_TYPE, size=REC_N_OBS, fill=100*REC_P_OUTSIDE_M), stroke=1.5, alpha=0.8) +
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +			
			scale_size(breaks=c(5,10,20,40,80), range=c(5,20))+
			scale_shape_manual(values=c('agrarian'=21, 'trading'=22, 'fisherfolk'=23)) +
			theme(legend.position='bottom', legend.box = "vertical") +
			guides(	pch=guide_legend(override.aes=list(size=7))) +
			labs(	x='', y='', 
					size="phylogenetically inferred transmissions\nwith transmitters from community",
					pch="community type",
					fill="proportion of transmitters from community\nwho have recipient outside community")
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_adjusted_anonymized.pdf'), w=7, h=7)		
}

RakaiFull.phylogeography.171122.recipient.hubs<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_prior23_min30_inference.rda"
	outfile.base		<- gsub('_inference.rda','',infile.inference)
	load(infile.inference)
	
	#
	#	transmission hubs crude estimates		
	z	<- copy(rtr2)
	z	<- z[, list(TR_OUTSIDE_COMM=length(which(REC_COMM_NUM_A!=TR_COMM_NUM_A)), TR_IN_COMM=length(which(REC_COMM_NUM_A==TR_COMM_NUM_A)) ), by=c('REC_COMM_NUM_A','REC_COMM_TYPE')]
	z[, TR_N:=TR_OUTSIDE_COMM+TR_IN_COMM]
	setnames(z, 'REC_COMM_NUM_A', 'COMM_NUM_A')
	z	<- merge(unique(zc, by='COMM_NUM_A'), z, by='COMM_NUM_A')
	ggmap(zm) +
			geom_point(data=z, aes(x=longitude, y=latitude, pch=COMM_TYPE, size=TR_N, fill=100*TR_OUTSIDE_COMM/TR_N), stroke=1.5, alpha=0.8) + 
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +
			scale_size(breaks=c(5,10,20,40,80), range=c(5,20))+
			scale_shape_manual(values=c('agrarian'=21, 'trading'=22, 'fisherfolk'=23)) +
			scale_colour_brewer(palette='Dark2') +			
			theme(legend.position='bottom', legend.box = "vertical") +
			guides(	pch=guide_legend(override.aes=list(size=7))) +
			labs(	x='', y='', 
					size="phylogenetically inferred transmissions\nwith recipients from community",
					pch="community type",
					fill="proportion of recipients from community\nwho have transmitter outside community")
	ggsave(file=paste0(outfile.base,'_hubs_recipients_crude.pdf'), w=7, h=7)
	
	
	#
	#	transmission hubs
	#	proportion of transmitters in community i that have a recipient outside community
	#	this is a summary on the estimated posterior density of the WAIFM matrix
	#	TODO missing community 94
	#
	#	get parameters of posterior under augmented likelihood
	z		<- dcb[, {
				z<- which(REC_COMM_NUM_A==TR_COMM_NUM_A)
				list(P_TROUTSIDE_ALPHA= sum(PI_IJ_ALPHA[-z]), P_TROUTSIDE_BETA= sum(PI_IJ_ALPHA[z]))
			}, by=c('REC_COMM_NUM_A','MONTE_CARLO_IT')]	
	#	aggregate and get quantiles
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- quantile(rbeta(length(P_TROUTSIDE_ALPHA)*mc.it, P_TROUTSIDE_ALPHA, P_TROUTSIDE_BETA), p=seq(0,1,0.01))
				list(P=seq(0,1,0.01), Q=unname(tmp))
			}, by='REC_COMM_NUM_A']
	#	subset to main quantities of interest
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), REC_COMM_NUM_A~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('TR_P_OUTSIDE_CL','TR_P_OUTSIDE_IL','TR_P_OUTSIDE_M','TR_P_OUTSIDE_IU','TR_P_OUTSIDE_CU'))
	#	add TR_OBS from TR_COMM_NUM_A
	z		<- merge(z, dc[, list(TR_N_OBS=sum(TR_OBS)), by='REC_COMM_NUM_A'], by='REC_COMM_NUM_A')
	setnames(z, 'REC_COMM_NUM_A', 'COMM_NUM_A')	
	#	now after analysis, remove 94 with unknown long / lat
	#z		<- subset(z, COMM_NUM!='94')
	z	<- merge(z, unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE))), by='COMM_NUM_A')
	ggplot(z, aes(x=COMM_NUM_A, middle=TR_P_OUTSIDE_M, min=TR_P_OUTSIDE_CL, max=TR_P_OUTSIDE_CU, lower=TR_P_OUTSIDE_IL, upper=TR_P_OUTSIDE_IU, fill=TR_P_OUTSIDE_M)) + 
			geom_boxplot(outlier.shape=NA, stat='identity') +
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +
			#geom_point(data=subset(z, BS==0), pch=2, colour='red') +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.2)) +
			facet_grid(~COMM_TYPE, scales='free_x', space='free_x') +
			labs(x='\ncommunity', y="transmitters outside community") +
			theme_bw() + guides(fill='none')
	ggsave(file=paste0(outfile.base,'_hubs_recipients_adjusted_confidenceintervals_anonmyized.pdf'), w=9, h=5)					
}

RakaiFull.phylogeography.171122.outinflows.commtype<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_prior23_min30_inference.rda"
	#infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl35_prior23_min30_inference.rda"
	outfile.base		<- gsub('_inference.rda','',infile.inference)
	load(infile.inference)
	#
	#	geography flows out > flows in for agrarian/trading/fisherolk communities
	#	adjusted		
	groups	<- c('agrarian','trading','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
				set(tmp, tmp[, which(COMM_TYPE!=group)], 'COMM_TYPE', paste0('non-',group))
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
				z[, FLOW:=paste0('from ',TR_COMM_TYPE,' to ',REC_COMM_TYPE)]
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- FLOW
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	transform as desired
				z[, STAT:= paste0('from ',group,' to non-',group,'\n/\nfrom non-',group,' to ',group)]				
				setnames(z, c(paste0('from ',group,' to non-',group),paste0('from non-',group,' to ',group)), c('A','D'))
				#	get quantiles
				#	just compare the diagonals, ie total flow.
				z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(A/D, p=seq(0,1,0.01)))), by='STAT']
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), STAT~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('FLOW_CL','FLOW_IL','FLOW_M','FLOW_IU','FLOW_CU'))
	z[, COMM_TYPE:= gsub('^from ([a-z]+).*$', '\\1', STAT)]
	ggplot(z, aes(x=COMM_TYPE, middle=FLOW_M, min=FLOW_CL, lower=FLOW_IL, upper=FLOW_IU, max=FLOW_CU, fill=STAT)) +
			geom_hline(yintercept=1, colour='grey50', size=1.5) +
			geom_boxplot(stat='identity') +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + theme(legend.position='bottom') +
			scale_y_log10(expand=c(0,0), breaks=c(1/4,1/3,1/2,2/3,1,3/2,2,3,4), labels=c('1/4','1/3','1/2','2/3','1','3/2','2','3','4')) +
			coord_flip() +
			labs(x='community type\n', y='\nflows out / flows in') +
			guides(fill='none')
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_flowsoutin_adjusted.pdf'), w=6, h=4)	
}

RakaiFull.phylogeography.180322.fishinland.netflows<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_phylogeography_core_inference.rda"	
	outfile.base		<- gsub('_inference.rda','',infile.inference)
	load(infile.inference)
	#
	#	geography flows out > flows in for agrarian/trading/fisherolk communities
	#	adjusted		
	groups	<- c('inland','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))				
				set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')				
				set(tmp, tmp[, which(COMM_TYPE!=group)], 'COMM_TYPE', paste0('non-',group))
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
				z[, FLOW:=paste0('from ',TR_COMM_TYPE,' to ',REC_COMM_TYPE)]
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- FLOW
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	transform as desired
				z[, STAT:= paste0('from ',group,' to non-',group,'\n/\nfrom non-',group,' to ',group)]				
				setnames(z, c(paste0('from ',group,' to non-',group),paste0('from non-',group,' to ',group)), c('A','D'))
				#	get quantiles
				#	just compare the diagonals, ie total flow.
				qs		<- c(0.025,0.25,0.5,0.75,0.975)
				qsn		<- c('FLOW_CL','FLOW_IL','FLOW_M','FLOW_IU','FLOW_CU')				
				z		<- z[, list(P=qsn, Q=unname(quantile(A/D, p=qs))), by='STAT']
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, STAT~P, value.var='Q')	
	z[, COMM_TYPE:= gsub('^from ([a-z]+).*$', '\\1', STAT)]
	z[, LABEL2:= paste0(round(FLOW_M, d=2), ' (',round(FLOW_CL,d=2),'-',round(FLOW_CU,d=2),')')]
	write.csv(z, row.names=FALSE, paste0(outfile.base,'_sourcesink_fishing_inland_adjusted.csv'))
	#
	ggplot(z, aes(x=COMM_TYPE, middle=FLOW_M, min=FLOW_CL, lower=FLOW_IL, upper=FLOW_IU, max=FLOW_CU, fill=COMM_TYPE)) +
			geom_hline(yintercept=1, colour='grey50', size=1.5) +
			geom_boxplot(stat='identity') +
			theme_bw() + theme(legend.position='bottom') +
			scale_y_log10(expand=c(0,0), breaks=c(1/4,1/3,1/2,2/3,1,3/2,2,3,4), labels=c('1/4','1/3','1/2','2/3','1','3/2','2','3','4')) +
			coord_flip() +
			scale_fill_manual(values=c('inland'='darkorange1','fisherfolk'='firebrick1')) +
			labs(x='community type\n', y='\nflows out / flows in') +
			guides(fill='none')
	ggsave(file=paste0(outfile.base,'_sourcesink_fishing_inland.pdf'), w=6, h=3)	
}

RakaiFull.phylogeography.180322.agefishinland.flows.odds<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_phylogeographyage_core_inference.rda"
	outfile.base		<- gsub('_inference.rda','',infile.inference)
	load(infile.inference)
	dcb[, TR_AGE_AT_MID_C:= factor(as.character(TR_AGE_AT_MID_C), levels=c('15-24','25-34','35+'), labels=c('15to24','25to34','35to50'))]
	dcb[, REC_AGE_AT_MID_C:= factor(as.character(REC_AGE_AT_MID_C), levels=c('15-24','25-34','35+'), labels=c('15to24','25to34','35to50'))]
	
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')		
	#
	#	geography who infects whom matrix  between fisherfolk, inland and gender
	#	adjusted P
	tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
							TR_MISS=sum(TR_MISS), 
							PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
					by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_AGE_AT_MID_C','REC_AGE_AT_MID_C','MONTE_CARLO_IT')]
	z[, FLOW:=paste0('from_',TR_COMM_TYPE,'_',TR_AGE_AT_MID_C,' to_',REC_COMM_TYPE,'_',REC_AGE_AT_MID_C)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]	
	z[, TR_COMM_TYPE:= gsub('from_([a-z]+)_([a-z0-9]+) to_([a-z]+)_([a-z0-9]+)','\\1',variable)]
	z[, TR_AGE_AT_MID_C:= gsub('from_([a-z]+)_([a-z0-9]+) to_([a-z]+)_([a-z0-9]+)','\\2',variable)]
	z[, REC_COMM_TYPE:= gsub('from_([a-z]+)_([a-z0-9]+) to_([a-z]+)_([a-z0-9]+)','\\3',variable)]	
	z[, REC_AGE_AT_MID_C:= gsub('from_([a-z]+)_([a-z0-9]+) to_([a-z]+)_([a-z0-9]+)','\\4',variable)]	
	z		<- dcast.data.table(z, TR_COMM_TYPE+TR_AGE_AT_MID_C+REC_COMM_TYPE+REC_AGE_AT_MID_C~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]	
	setkey(z, REC_COMM_TYPE, REC_AGE_AT_MID_C, TR_COMM_TYPE, TR_AGE_AT_MID_C)
	z[, STAT:='joint']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- copy(z)	
	dt		<- dcast.data.table(z, REC_COMM_TYPE+REC_AGE_AT_MID_C~TR_COMM_TYPE+TR_AGE_AT_MID_C, value.var='LABEL2')
	write.csv(dt, row.names=FALSE, paste0(outfile.base,'_fishing_inland_age_joint_adjusted.csv'))
	
	#
	#	WAIFM
	#
	groups	<- c('inland_15to24','inland_25to34','inland_35to50','fisherfolk_15to24','fisherfolk_25to34','fisherfolk_35to50')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))				
				set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')				
				z		<- subset(z, TR_COMM_TYPE==gsub('^([a-z]+)_([a-z0-9]+)$','\\1',group) & TR_AGE_AT_MID_C==gsub('^([a-z]+)_([a-z0-9]+)$','\\2',group))				
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('REC_COMM_TYPE','REC_AGE_AT_MID_C','MONTE_CARLO_IT')]				
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- paste0('to ',REC_COMM_TYPE,'_',REC_AGE_AT_MID_C)
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	get quantiles
				z		<- melt(z, id.vars='MONTE_CARLO_IT')
				z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]				
				z[, REC_COMM_TYPE:= gsub('to ([a-z]+)_([a-z0-9]+)','\\1',variable)]	
				z[, REC_AGE_AT_MID_C:= gsub('to ([a-z]+)_([a-z0-9]+)','\\2',variable)]	
				z[, TR_COMM_TYPE:= gsub('^([a-z]+)_([a-z0-9]+)$','\\1',group)]
				z[, TR_AGE_AT_MID_C:= gsub('^([a-z]+)_([a-z0-9]+)$','\\2',group)]
				z
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, REC_COMM_TYPE+REC_AGE_AT_MID_C+TR_COMM_TYPE+TR_AGE_AT_MID_C~P, value.var='Q')		
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, TR_AGE_AT_MID_C, REC_COMM_TYPE, REC_AGE_AT_MID_C)
	z[, STAT:='waifm']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans,z)
	dt		<- dcast.data.table(z, REC_COMM_TYPE+REC_AGE_AT_MID_C~TR_COMM_TYPE+TR_AGE_AT_MID_C, value.var='LABEL2')
	write.csv(dt, row.names=FALSE, paste0(outfile.base,'_fishing_inland_age_waifm_adjusted.csv'))
	
	#
	#	sources
	#
	groups	<- c('inland_15to24','inland_25to34','inland_35to50','fisherfolk_15to24','fisherfolk_25to34','fisherfolk_35to50')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))				
				set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')				
				z		<- subset(z, REC_COMM_TYPE==gsub('^([a-z]+)_([a-z0-9]+)$','\\1',group) & REC_AGE_AT_MID_C==gsub('^([a-z]+)_([a-z0-9]+)$','\\2',group))				
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('TR_COMM_TYPE','TR_AGE_AT_MID_C','MONTE_CARLO_IT')]				
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- paste0('from ',TR_COMM_TYPE,'_',TR_AGE_AT_MID_C)
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	get quantiles
				z		<- melt(z, id.vars='MONTE_CARLO_IT')
				z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]				
				z[, TR_COMM_TYPE:= gsub('from ([a-z]+)_([a-z0-9]+)','\\1',variable)]	
				z[, TR_AGE_AT_MID_C:= gsub('from ([a-z]+)_([a-z0-9]+)','\\2',variable)]	
				z[, REC_COMM_TYPE:= gsub('^([a-z]+)_([a-z0-9]+)$','\\1',group)]
				z[, REC_AGE_AT_MID_C:= gsub('^([a-z]+)_([a-z0-9]+)$','\\2',group)]
				z
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, REC_COMM_TYPE+REC_AGE_AT_MID_C+TR_COMM_TYPE+TR_AGE_AT_MID_C~P, value.var='Q')		
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_COMM_TYPE, REC_AGE_AT_MID_C, TR_COMM_TYPE, TR_AGE_AT_MID_C)
	z[, STAT:='sources']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans,z)
	dt		<- dcast.data.table(z, TR_COMM_TYPE+TR_AGE_AT_MID_C~REC_COMM_TYPE+REC_AGE_AT_MID_C, value.var='LABEL2')
	write.csv(dt, row.names=FALSE, paste0(outfile.base,'_fishing_inland_age_sources_adjusted.csv'))
	
	#
	#	odds male/female
	#
	tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
							TR_MISS=sum(TR_MISS), 
							PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
					by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_AGE_AT_MID_C','MONTE_CARLO_IT')]
	z[, FLOW:=paste0('from_',TR_COMM_TYPE,'_',TR_AGE_AT_MID_C,'_to_',REC_COMM_TYPE)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)			
			}, by=c('MONTE_CARLO_IT')]	
	z[, oddsYO_from_inland_to_inland:=  from_inland_15to24_to_inland / from_inland_35to50_to_inland]
	z[, oddsYO_from_fisherfolk_to_inland:=  from_fisherfolk_15to24_to_inland / from_fisherfolk_35to50_to_inland]	
	z[, oddsYO_from_inland_to_fisherfolk:=  from_inland_15to24_to_fisherfolk / from_inland_35to50_to_fisherfolk]
	z[, oddsYO_from_fisherfolk_to_fisherfolk:=  from_fisherfolk_15to24_to_fisherfolk / from_fisherfolk_35to50_to_fisherfolk]		
	z[, oddsMO_from_inland_to_inland:=  from_inland_25to34_to_inland / from_inland_35to50_to_inland]
	z[, oddsMO_from_fisherfolk_to_inland:=  from_fisherfolk_25to34_to_inland / from_fisherfolk_35to50_to_inland]	
	z[, oddsMO_from_inland_to_fisherfolk:=  from_inland_25to34_to_fisherfolk / from_inland_35to50_to_fisherfolk]
	z[, oddsMO_from_fisherfolk_to_fisherfolk:=  from_fisherfolk_25to34_to_fisherfolk / from_fisherfolk_35to50_to_fisherfolk]	
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z		<- subset(z, grepl('odds', variable))
	z[, STAT:= gsub('(odds[A-Z]+)_from_([a-z]+)_to_([a-z]+)','\\1',variable)]
	z[, TR_COMM_TYPE:= gsub('(odds[A-Z]+)_from_([a-z]+)_to_([a-z]+)','\\2',variable)]	
	z[, REC_COMM_TYPE:= gsub('(odds[A-Z]+)_from_([a-z]+)_to_([a-z]+)','\\3',variable)]
	
	z		<- dcast.data.table(z, STAT+TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),' - ',round(CU,d=2),']')]
	z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'-',round(CU,d=2),')')]	
	setkey(z, STAT, REC_COMM_TYPE, TR_COMM_TYPE)
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans,z, fill=TRUE)
	dt		<- dcast.data.table(z, STAT+REC_COMM_TYPE~TR_COMM_TYPE, value.var='LABEL2')
	write.csv(dt, row.names=FALSE, paste0(outfile.base,'_fishing_inland_age_oddMF_adjusted.csv'))
		
	save(ans, file=paste0(outfile.base,'_fishing_inland_age_results.rda'))
	
	tmp		<- subset(ans, STAT=='sources')
	ggplot(tmp, aes(x=paste0(TR_COMM_TYPE,TR_AGE_AT_MID_C), y=paste0(REC_COMM_TYPE,REC_AGE_AT_MID_C), fill=M)) + 
		geom_tile() + 
		scale_fill_gradient(low = "white", high = "steelblue", labels=scales:::percent) +
		geom_text(aes(label=LABEL), size=2) +
		theme_grey(base_size = 9) + 
		labs(x = "from\n",  y = "to\n", fill='sources of\nHIV-1 infection') + 
		scale_x_discrete(expand = c(0, 0), position = "top") +
		scale_y_discrete(expand = c(0, 0)) + 		
		theme(axis.text.x= element_text(size = 9 * 0.8, angle=30, hjust=0 ))
	ggsave(file=paste0(outfile.base,'_fishing_inland_age_sources_adjusted.pdf'), w=7, h=6)
}


RakaiFull.phylogeography.180322.genderfishinland.results<- function()
{
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_phylogeography_core_fishing_inland_gender_results.rda")
}
	
RakaiFull.phylogeography.180322.genderfishinland.flows.odds<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_phylogeography_core_inference.rda"
	outfile.base		<- gsub('_inference.rda','',infile.inference)
	load(infile.inference)
	
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')	
	
	#
	#	geography who infects whom matrix  between fisherfolk, inland and gender
	#	adjusted P
	tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_SEX','REC_SEX','MONTE_CARLO_IT')]
	z[, FLOW:=paste0('from_',TR_COMM_TYPE,'_',TR_SEX,' to_',REC_COMM_TYPE,'_',REC_SEX)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]	
	z[, TR_COMM_TYPE:= gsub('from_([a-z]+)_([A-Z]+) to_([a-z]+)_([A-Z]+)','\\1',variable)]
	z[, TR_SEX:= gsub('from_([a-z]+)_([A-Z]+) to_([a-z]+)_([A-Z]+)','\\2',variable)]
	z[, REC_COMM_TYPE:= gsub('from_([a-z]+)_([A-Z]+) to_([a-z]+)_([A-Z]+)','\\3',variable)]	
	z[, REC_SEX:= gsub('from_([a-z]+)_([A-Z]+) to_([a-z]+)_([A-Z]+)','\\4',variable)]	
	z		<- dcast.data.table(z, TR_COMM_TYPE+TR_SEX+REC_COMM_TYPE+REC_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]	
	setkey(z, REC_COMM_TYPE, TR_COMM_TYPE, TR_SEX)
	z[, STAT:='joint']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- copy(z)	
	dt		<- dcast.data.table(z, REC_COMM_TYPE~TR_COMM_TYPE+TR_SEX, value.var='LABEL2')
	write.csv(dt, row.names=FALSE, paste0(outfile.base,'_fishing_inland_gender_2x2_adjusted.csv'))
	
	#
	#	WAIFM
	#
	groups	<- c('inland_M','inland_F','fisherfolk_M','fisherfolk_F')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))				
				set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')				
				z		<- subset(z, TR_COMM_TYPE==gsub('^([a-z]+)_([A-Z]+)$','\\1',group) & TR_SEX==gsub('^([a-z]+)_([A-Z]+)$','\\2',group))				
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('REC_COMM_TYPE','REC_SEX','MONTE_CARLO_IT')]				
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- paste0('to ',REC_COMM_TYPE,'_',REC_SEX)
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	get quantiles
				z		<- melt(z, id.vars='MONTE_CARLO_IT')
				z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]				
				z[, REC_COMM_TYPE:= gsub('to ([a-z]+)_([A-Z]+)','\\1',variable)]	
				z[, REC_SEX:= gsub('to ([a-z]+)_([A-Z]+)','\\2',variable)]	
				z[, TR_COMM_TYPE:= gsub('^([a-z]+)_([A-Z]+)$','\\1',group)]
				z[, TR_SEX:= gsub('^([a-z]+)_([A-Z]+)$','\\2',group)]
				z
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, REC_COMM_TYPE+REC_SEX+TR_COMM_TYPE+TR_SEX~P, value.var='Q')		
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, TR_SEX, REC_COMM_TYPE)
	z[, STAT:='waifm']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans,z)
	dt		<- dcast.data.table(z, REC_COMM_TYPE~TR_COMM_TYPE+TR_SEX, value.var='LABEL2')
	write.csv(dt, row.names=FALSE, paste0(outfile.base,'_fishing_inland_gender_waifm_adjusted.csv'))
	
	#
	#	sources
	#
	groups	<- c('inland_M','inland_F','fisherfolk_M','fisherfolk_F')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))				
				set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')				
				z		<- subset(z, REC_COMM_TYPE==gsub('^([a-z]+)_([A-Z]+)$','\\1',group) & REC_SEX==gsub('^([a-z]+)_([A-Z]+)$','\\2',group))				
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('TR_COMM_TYPE','TR_SEX','MONTE_CARLO_IT')]				
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- paste0('from ',TR_COMM_TYPE,'_',TR_SEX)
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	get quantiles
				z		<- melt(z, id.vars='MONTE_CARLO_IT')
				z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]				
				z[, TR_COMM_TYPE:= gsub('from ([a-z]+)_([A-Z]+)','\\1',variable)]	
				z[, TR_SEX:= gsub('from ([a-z]+)_([A-Z]+)','\\2',variable)]	
				z[, REC_COMM_TYPE:= gsub('^([a-z]+)_([A-Z]+)$','\\1',group)]
				z[, REC_SEX:= gsub('^([a-z]+)_([A-Z]+)$','\\2',group)]
				z
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, REC_COMM_TYPE+REC_SEX+TR_COMM_TYPE+TR_SEX~P, value.var='Q')		
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_COMM_TYPE, REC_SEX, TR_COMM_TYPE)
	z[, STAT:='sources']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans,z)
	dt		<- dcast.data.table(z, TR_COMM_TYPE~REC_COMM_TYPE+REC_SEX, value.var='LABEL2')
	write.csv(dt, row.names=FALSE, paste0(outfile.base,'_fishing_inland_gender_sources_adjusted.csv'))
	
	#
	#	odds male/female
	#
	tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_SEX','REC_SEX','MONTE_CARLO_IT')]
	z[, FLOW:=paste0('from_',TR_COMM_TYPE,'_',TR_SEX,'_to_',REC_COMM_TYPE,'_',REC_SEX)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)			
			}, by=c('MONTE_CARLO_IT')]
	
	z[, oddsMF_from_inland_to_inland:=  from_inland_M_to_inland_F / from_inland_F_to_inland_M]
	z[, oddsMF_from_fisherfolk_to_inland:=  from_fisherfolk_M_to_inland_F / from_fisherfolk_F_to_inland_M]	
	z[, oddsMF_from_inland_to_fisherfolk:=  from_inland_M_to_fisherfolk_F / from_inland_F_to_fisherfolk_M]
	z[, oddsMF_from_fisherfolk_to_fisherfolk:=  from_fisherfolk_M_to_fisherfolk_F / from_fisherfolk_F_to_fisherfolk_M]	
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z		<- subset(z, grepl('odds', variable))
	z[, TR_COMM_TYPE:= gsub('oddsMF_from_([a-z]+)_to_([a-z]+)','\\1',variable)]	
	z[, REC_COMM_TYPE:= gsub('oddsMF_from_([a-z]+)_to_([a-z]+)','\\2',variable)]			
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),'% - ',round(CU,d=2),'%]')]
	z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'%-',round(CU,d=2),'%)')]	
	setkey(z, REC_COMM_TYPE, TR_COMM_TYPE)
	z[, STAT:='oddsMF']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans,z, fill=TRUE)
	dt		<- dcast.data.table(z, REC_COMM_TYPE~TR_COMM_TYPE, value.var='LABEL2')
	write.csv(dt, row.names=FALSE, paste0(outfile.base,'_fishing_inland_gender_oddMF_adjusted.csv'))
	
	save(ans, file=paste0(outfile.base,'_fishing_inland_gender_results.rda'))
}


RakaiFull.phylogeography.180322.fishinland.flows<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_phylogeography_core_inference.rda"
	outfile.base		<- gsub('_inference.rda','',infile.inference)
	load(infile.inference)
	
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')
	
	#
	#	geography who infects whom matrix  between fisherfolk and others
	#	adjusted P
	tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0('from_',TR_COMM_TYPE,' to_',REC_COMM_TYPE)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]
	z[, TR_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\1',variable)]
	z[, REC_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\2',variable)]		
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_COMM_TYPE, TR_COMM_TYPE, TR_SEX)
	z[, STAT:='joint']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- copy(z)	
	dt		<- dcast.data.table(z, REC_COMM_TYPE~TR_COMM_TYPE, value.var='LABEL2')
	write.csv(dt, row.names=FALSE, paste0(outfile.base,'_fishing_inland_2x2_adjusted.csv'))
	
	ggplot(z, aes(x=factor(REC_COMM_TYPE, levels=c('to_fisherfolk','to_inland'), labels=c('fishing\nsites','inland\ncommunities')), y=factor(TR_COMM_TYPE, levels=c('from_inland','from_fisherfolk'), labels=c('inland\ncommunities','fishing\nsites')))) + 
			geom_point(aes(size=PADJ_M), colour='grey80') +
			geom_text(aes(label=LABEL), nudge_x=0, nudge_y=0, size=3, colour='black') +			
			theme_bw() + 
			scale_size(range=c(5, 50)) +
			labs(x='\nto',y='from\n') +
			guides(size='none')
	ggsave(file=paste0(outfile.base,'_flows_fishing2x2_adjusted.pdf'), w=4, h=4)
	
	#
	#	WAIFM
	#
	groups	<- c('inland','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))				
				set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				z		<- subset(z, TR_COMM_TYPE==group)				
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('REC_COMM_TYPE','MONTE_CARLO_IT')]				
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- REC_COMM_TYPE
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	get quantiles
				z		<- melt(z, id.vars=c('MONTE_CARLO_IT'), variable.name='REC_COMM_TYPE')								
				z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by='REC_COMM_TYPE']
				z[, TR_COMM_TYPE:= group]
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE)
	z[, STAT:='waifm']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans,z)	
	dt		<- dcast.data.table(z, REC_COMM_TYPE~TR_COMM_TYPE, value.var='LABEL2')
	write.csv(dt, row.names=FALSE, paste0(outfile.base,'_fishing_inland_waifm_adjusted.csv'))
	
	#
	#	sources
	#
	groups	<- c('inland','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))				
				set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				z		<- subset(z, REC_COMM_TYPE==group)				
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('TR_COMM_TYPE','MONTE_CARLO_IT')]				
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- TR_COMM_TYPE
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	get quantiles
				z		<- melt(z, id.vars=c('MONTE_CARLO_IT'), variable.name='TR_COMM_TYPE')								
				z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by='TR_COMM_TYPE']
				z[, REC_COMM_TYPE:= group]
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_COMM_TYPE, TR_COMM_TYPE)
	z[, STAT:='sources']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans,z)	
	dt		<- dcast.data.table(z, TR_COMM_TYPE~REC_COMM_TYPE, value.var='LABEL2')
	write.csv(dt, row.names=FALSE, paste0(outfile.base,'_fishing_inland_sources_adjusted.csv'))
	
	save(ans, file=paste0(outfile.base,'_fishing_inland_results.rda'))
}

RakaiFull.phylogeography.171122.flowsbetweencommunities<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_prior23_min30_inference.rda"
	outfile.base		<- gsub('_inference.rda','',infile.inference)
	load(infile.inference)
	
	#
	#	geography who infects whom matrix  
	#	crude
	tmp		<- rtr2[,list(N=length(unique(PAIRID))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]
	tmp[, P_CELL:= N/sum(N)]
	tmp		<- merge(tmp, tmp[, list(P_REC= N/sum(N), REC_COMM_TYPE=REC_COMM_TYPE), by='TR_COMM_TYPE'], by=c('TR_COMM_TYPE','REC_COMM_TYPE'))
	tmp		<- merge(tmp, tmp[, list(P_TR= N/sum(N), TR_COMM_TYPE=TR_COMM_TYPE), by='REC_COMM_TYPE'], by=c('TR_COMM_TYPE','REC_COMM_TYPE'))
	#tmp[, LABEL:= paste0(N, ' (',round(P_CELL,d=2)*100,'%)\ntransmitters: ',round(P_TR,d=2)*100,'%\nrecipients: ',round(P_REC,d=2)*100,'%')]
	tmp[, LABEL:= paste0(N, '\n(',round(P_TR,d=2)*100,'%)')]
	#tmp[, LABEL:= paste0(N, '\n(',round(P_CELL,d=2)*100,'%)')]
	ggplot(tmp, aes(x=factor(REC_COMM_TYPE),y=factor(TR_COMM_TYPE))) + 
			geom_point(aes(size=N), colour='grey80') +
			geom_text(aes(label=LABEL), nudge_x=0, nudge_y=0, size=3, colour='black') +			
			theme_bw() + 
			scale_size(range = c(5, 50)) +
			labs(x='\nlocation likely recipient',y='location likely transmitter\n') +
			guides(size='none')
	ggsave(file=paste0(outfile.base,'_commtype_3x3.pdf'), w=5, h=5)
	
	set(tmp, NULL, 'REC_COMM_TYPE', tmp[,paste0('to_',REC_COMM_TYPE)])
	set(tmp, NULL, 'TR_COMM_TYPE', tmp[,paste0('from_',TR_COMM_TYPE)])
	tmp		<- suppressWarnings(melt(tmp, id.vars=c('REC_COMM_TYPE','TR_COMM_TYPE'), measure.vars=c('N','P_CELL')))
	tmp		<- dcast.data.table(tmp, variable+TR_COMM_TYPE~REC_COMM_TYPE, value.var='value')
	write.csv(tmp, row.names=FALSE, paste0(outfile.base,'_commtype_3x3_raw.csv'))
	
	#
	#	geography who infects whom matrix  
	#	adjusted P
	tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	z		<- z[, list(TR_OBS=sum(TR_OBS), TR_MISS=sum(TR_MISS), PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0('from_',TR_COMM_TYPE,' to_',REC_COMM_TYPE)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z[, TR_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\1',variable)]
	z[, REC_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\2',variable)]	
	z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(value, p=seq(0,1,0.01)))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]	
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('PADJ_CL','PADJ_IL','PADJ_M','PADJ_IU','PADJ_CU'))		
	ans		<- melt(z, id.vars=c('TR_COMM_TYPE','REC_COMM_TYPE'), measure.vars=c('PADJ_M','PADJ_CL','PADJ_CU'))
	ans		<- dcast.data.table(ans, variable+REC_COMM_TYPE~TR_COMM_TYPE, value.var='value')
	#	adjusted N
	tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	z		<- z[, list(TR_ADJ=sum(TR_OBS)+sum(TR_MISS)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	set(z, NULL, 'REC_COMM_TYPE', z[,paste0('to_',REC_COMM_TYPE)])
	set(z, NULL, 'TR_COMM_TYPE', z[,paste0('from_',TR_COMM_TYPE)])
	z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(TR_ADJ, p=seq(0,1,0.01), type=1))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]	
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('NADJ_CL','NADJ_IL','NADJ_M','NADJ_IU','NADJ_CU'))
	z		<- melt(z, id.vars=c('TR_COMM_TYPE','REC_COMM_TYPE'), measure.vars=c('NADJ_M','NADJ_CL','NADJ_CU'))
	z		<- dcast.data.table(z, variable+REC_COMM_TYPE~TR_COMM_TYPE, value.var='value')
	ans		<- rbind(z, ans)
	write.csv(ans, row.names=FALSE, paste0(outfile.base,'_commtype_3x3_adjusted.csv'))
	
	tmp		<- subset(ans, grepl('PADJ', variable))
	tmp		<- melt(tmp, id.vars=c('variable','REC_COMM_TYPE'), variable.name='TR_COMM_TYPE')
	tmp		<- dcast.data.table(tmp, REC_COMM_TYPE+TR_COMM_TYPE~variable, value.var='value')	
	tmp[, LABEL:= paste0(round(PADJ_M*100, d=1), '%\n[',round(PADJ_CL*100,d=1),'% - ',round(PADJ_CU*100,d=1),'%]')]
	
	ggplot(tmp, aes(x=factor(REC_COMM_TYPE), y=factor(TR_COMM_TYPE))) + 
			geom_point(aes(size=PADJ_M), colour='grey80') +
			geom_text(aes(label=LABEL), nudge_x=0, nudge_y=0, size=3, colour='black') +			
			theme_bw() + 
			scale_size(range=c(5, 50)) +
			labs(x='\nlocation likely recipient',y='location likely transmitter\n') +
			guides(size='none')
	ggsave(file=paste0(outfile.base,'_flows_commtype3x3_adjusted.pdf'), w=5, h=5)
	
	
	#
	#	geography who infects whom matrix  between fisherfolk and others
	#	adjusted P
	tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'other')
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	z		<- z[, list(TR_OBS=sum(TR_OBS), TR_MISS=sum(TR_MISS), PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0('from_',TR_COMM_TYPE,' to_',REC_COMM_TYPE)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z[, TR_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\1',variable)]
	z[, REC_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\2',variable)]	
	z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(value, p=seq(0,1,0.01)))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]	
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('PADJ_CL','PADJ_IL','PADJ_M','PADJ_IU','PADJ_CU'))		
	ans		<- melt(z, id.vars=c('TR_COMM_TYPE','REC_COMM_TYPE'), measure.vars=c('PADJ_M','PADJ_CL','PADJ_CU'))
	ans		<- dcast.data.table(ans, variable+REC_COMM_TYPE~TR_COMM_TYPE, value.var='value')	
	#	adjusted N
	tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	z		<- z[, list(TR_ADJ=sum(TR_OBS)+sum(TR_MISS)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	set(z, NULL, 'REC_COMM_TYPE', z[,paste0('to_',REC_COMM_TYPE)])
	set(z, NULL, 'TR_COMM_TYPE', z[,paste0('from_',TR_COMM_TYPE)])
	z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(TR_ADJ, p=seq(0,1,0.01), type=1))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]	
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('NADJ_CL','NADJ_IL','NADJ_M','NADJ_IU','NADJ_CU'))
	z		<- melt(z, id.vars=c('TR_COMM_TYPE','REC_COMM_TYPE'), measure.vars=c('NADJ_M','NADJ_CL','NADJ_CU'))
	z		<- dcast.data.table(z, variable+REC_COMM_TYPE~TR_COMM_TYPE, value.var='value')
	ans		<- rbind(z, ans)
	write.csv(ans, row.names=FALSE, paste0(outfile.base,'_commtype_3x3_adjusted.csv'))
	
	tmp		<- subset(ans, grepl('PADJ', variable))
	tmp		<- melt(tmp, id.vars=c('variable','REC_COMM_TYPE'), variable.name='TR_COMM_TYPE')
	tmp		<- dcast.data.table(tmp, REC_COMM_TYPE+TR_COMM_TYPE~variable, value.var='value')	
	tmp[, LABEL:= paste0(round(PADJ_M*100, d=1), '%\n[',round(PADJ_CL*100,d=1),'% - ',round(PADJ_CU*100,d=1),'%]')]
	ggplot(tmp, aes(x=factor(REC_COMM_TYPE, levels=c('to_fisherfolk','to_other'), labels=c('fishing\ncommunities','other\ncommunities')), y=factor(TR_COMM_TYPE, levels=c('from_other','from_fisherfolk'), labels=c('other\ncommunities','fishing\ncommunities')))) + 
			geom_point(aes(size=PADJ_M), colour='grey80') +
			geom_text(aes(label=LABEL), nudge_x=0, nudge_y=0, size=3, colour='black') +			
			theme_bw() + 
			scale_size(range=c(5, 50)) +
			labs(x='\nto',y='from\n') +
			guides(size='none')
	ggsave(file=paste0(outfile.base,'_flows_fishing2x2_adjusted.pdf'), w=4, h=4)
	
}

RakaiFull.phylogeography.171122<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	
	#
	#	geography transmitter flows into agrarian/trading/fisherolk recipient communities
	#	crude	
	tmp		<- rtr2[,list(N=length(unique(PAIRID))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]
	tmp[, P_CELL:= N/sum(N)]
	tmp		<- merge(tmp, tmp[, list(P_REC= N/sum(N), REC_COMM_TYPE=REC_COMM_TYPE), by='TR_COMM_TYPE'], by=c('TR_COMM_TYPE','REC_COMM_TYPE'))
	tmp		<- merge(tmp, tmp[, list(P_TR= N/sum(N), TR_COMM_TYPE=TR_COMM_TYPE), by='REC_COMM_TYPE'], by=c('TR_COMM_TYPE','REC_COMM_TYPE'))
	#tmp[, LABEL:= paste0(N, ' (',round(P_CELL,d=2)*100,'%)\ntransmitters: ',round(P_TR,d=2)*100,'%\nrecipients: ',round(P_REC,d=2)*100,'%')]
	tmp[, LABEL:= paste0(N, '\n(',round(P_TR,d=2)*100,'%)')]
	ggplot(tmp, aes(x=REC_COMM_TYPE, y=P_TR, fill=TR_COMM_TYPE)) + 
			geom_bar(stat='identity', position='dodge') +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + 
			scale_y_continuous(labels=scales::percent, limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.2)) +			
			labs(x='\nlocation likely recipient',y='location likely transmitter\n',fill='transmitter from') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_sources_crude_prop.pdf'), w=6, h=5)
	ggplot(tmp, aes(x=REC_COMM_TYPE, y=N, fill=TR_COMM_TYPE)) + 
			geom_bar(stat='identity', position='dodge') +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + 
			scale_y_continuous() +			
			labs(x='\nlocation likely recipient',y='location likely transmitter\n',fill='transmitter from') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_sources_crude_count.pdf'), w=6, h=5)	
	ggplot(tmp, aes(x=TR_COMM_TYPE, y=P_REC, fill=REC_COMM_TYPE)) + 
			geom_bar(stat='identity', position='dodge') +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + 
			scale_y_continuous(labels=scales::percent, limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.2)) +			
			labs(x='\nlocation likely transmitters',y='location likely recipients\n',fill='recipient in') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_destinations_crude_prop.pdf'), w=6, h=5)
	ggplot(tmp, aes(x=TR_COMM_TYPE, y=N, fill=REC_COMM_TYPE)) + 
			geom_bar(stat='identity', position='dodge') +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + 
			scale_y_continuous() +			
			labs(x='\nlocation likely transmitter',y='location likely recipient\n',fill='recipient in') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_destinations_crude_count.pdf'), w=6, h=5)	
	
	#
	#	geography transmitter flows into agrarian/trading/fisherolk recipient communities
	#	adjusted		
	groups	<- c('agrarian','trading','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- subset(zc, COMM_TYPE==group)$COMM_NUM_A		
				z		<- subset(dcb, REC_COMM_NUM_A%in%tmp)[, list(PI_ITYPE_ALPHA= sum(PI_IJ_ALPHA)), by=c('TR_COMM_NUM_A','MONTE_CARLO_IT')]
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				z		<- z[, list(PI_TYPETYPE_ALPHA=sum(PI_ITYPE_ALPHA)), by=c('TR_COMM_TYPE','MONTE_CARLO_IT')]	
				#	aggregate and get quantiles
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_TYPETYPE_ALPHA)
							colnames(tmp)	<- TR_COMM_TYPE
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				z		<- melt(z, id.vars='MONTE_CARLO_IT', variable.name='TR_COMM_TYPE', value.name='PI_TYPETYPE')
				z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(PI_TYPETYPE, p=seq(0,1,0.01)))), by='TR_COMM_TYPE']
				#	subset to main quantities of interest
				z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_TYPE~P, value.var='Q')
				z[, REC_COMM_TYPE:=group]
				z
			})
	z		<- do.call('rbind',z)
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('TR_CL','TR_IL','TR_M','TR_IU','TR_CU'))
	ggplot(z, aes(x=REC_COMM_TYPE, y=TR_M, ymin=TR_CL, lower=TR_IL, upper=TR_IU, ymax=TR_CU, fill=TR_COMM_TYPE)) + 
			#geom_boxplot(stat='identity', position=position_dodge(width=0.9)) +
			geom_bar(position='dodge', stat='identity') +
			geom_errorbar(width=0.2, position=position_dodge(width=0.9)) +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous(labels=scales::percent, limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.2)) +			
			labs(x='\nlocation likely recipient',y='sources of transmissions\n',fill='transmitter from') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_sources_adjusted.pdf'), w=6, h=5)	
	#
	#	geography transmitter flows from agrarian/trading/fisherolk transmitter communities
	#	adjusted		
	groups	<- c('agrarian','trading','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- subset(zc, COMM_TYPE==group)$COMM_NUM_A		
				z		<- subset(dcb, TR_COMM_NUM_A%in%tmp)[, list(PI_ITYPE_ALPHA= sum(PI_IJ_ALPHA)), by=c('REC_COMM_NUM_A','MONTE_CARLO_IT')]
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(z, tmp, by='REC_COMM_NUM_A')
				z		<- z[, list(PI_TYPETYPE_ALPHA=sum(PI_ITYPE_ALPHA)), by=c('REC_COMM_TYPE','MONTE_CARLO_IT')]	
				#	aggregate and get quantiles
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_TYPETYPE_ALPHA)
							colnames(tmp)	<- REC_COMM_TYPE
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				z		<- melt(z, id.vars='MONTE_CARLO_IT', variable.name='REC_COMM_TYPE', value.name='PI_TYPETYPE')
				z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(PI_TYPETYPE, p=seq(0,1,0.01)))), by='REC_COMM_TYPE']
				#	subset to main quantities of interest
				z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), REC_COMM_TYPE~P, value.var='Q')
				z[, TR_COMM_TYPE:=group]
				z
			})
	z		<- do.call('rbind',z)
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('REC_CL','REC_IL','REC_M','REC_IU','REC_CU'))
	ggplot(z, aes(x=TR_COMM_TYPE, y=REC_M, ymin=REC_CL, lower=REC_IL, upper=REC_IU, ymax=REC_CU, fill=REC_COMM_TYPE)) + 
			#geom_boxplot(stat='identity', position=position_dodge(width=0.9)) +
			geom_bar(position='dodge', stat='identity') +
			geom_errorbar(width=0.2, position=position_dodge(width=0.9)) +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous(labels=scales::percent, limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.2)) +			
			labs(x='\nlocation likely transmitter',y='destination of transmissions\n',fill='recipients in') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_destinations_barplot_adjusted.pdf'), w=6, h=5)	
	
	
	
	
	
	
	
	#
	#	geography transmitters from outside community
	z	<- copy(rtr2)
	z[, FROM_OUTSIDE:= factor(REC_COMM_NUM_A!=TR_COMM_NUM_A, levels=c(TRUE,FALSE), labels=c('outside community','same community'))]
	set(z, NULL, 'REC_COMM_TYPE', z[, factor(REC_COMM_TYPE)])
	ggplot(z, aes(x=REC_COMM_TYPE, fill=FROM_OUTSIDE)) + geom_bar() + 
			theme_bw() + 
			labs(x='\ncommunity type of recipient', fill='transmitter from')
	ggsave(file=paste0(outfile.base,'_extra_community.pdf'), w=5, h=5)
	ggplot(z, aes(x=REC_COMM_TYPE, fill=FROM_OUTSIDE)) + geom_bar(position='fill') + 
			theme_bw() + 
			scale_y_continuous(labels=scales::percent, expand=c(0,0), breaks=seq(0,1,0.2)) +
			labs(x='\ncommunity type of recipient', y='transmitters\n', fill='transmitter from')
	ggsave(file=paste0(outfile.base,'_extra_community_props.pdf'), w=5, h=5)
	#
	#	odds for unequal transmission flows agrarian/fisherfolk
	z	<- z[, list(N=length(REC_RID)), by=c('REC_COMM_TYPE','FROM_OUTSIDE')]
	tmp	<- dcast.data.table(subset(z, REC_COMM_TYPE!='trading'), FROM_OUTSIDE~REC_COMM_TYPE, value.var='N')
	zz	<- as.matrix(tmp[, 2:3, with=FALSE])
	rownames(zz)	<- tmp[[1]]
	fisher.test(zz)
	#STAGE 1
	#data:  zz
	#p-value = 0.0113
	#alternative hypothesis: true odds ratio is not equal to 1
	#95 percent confidence interval:
	#1.137942 3.390227
	#sample estimates:
	#odds ratio 
	#1.963628
	
	#STAGE 2
	#p-value = 0.0002936
	#alternative hypothesis: true odds ratio is not equal to 1
	#95 percent confidence interval:
	#1.698675 7.052786
	#sample estimates:
	#odds ratio 
	#3.427573 
	
	#
	#	geography recipients in different community
	z	<- copy(rtr2)
	z[, FROM_OUTSIDE:= factor(REC_COMM_NUM_A!=TR_COMM_NUM_A, levels=c(TRUE,FALSE), labels=c('outside community','same community'))]
	set(z, NULL, 'TR_COMM_TYPE', z[, factor(TR_COMM_TYPE)])
	ggplot(z, aes(x=TR_COMM_TYPE, fill=FROM_OUTSIDE)) + geom_bar() + 
			theme_bw() + 
			labs(x='\ncommunity type of transmitters', fill='recipient')
	ggsave(file=paste0(outfile.base,'_extra_community_of_recipients.pdf'), w=5, h=5)	
	ggplot(z, aes(x=TR_COMM_TYPE, fill=FROM_OUTSIDE)) + geom_bar(position='fill') + 
			theme_bw() + 
			scale_y_continuous(labels=scales::percent, expand=c(0,0), breaks=seq(0,1,0.2)) +
			labs(x='\ncommunity type of transmitters', y='recipients\n', fill='recipient')
	ggsave(file=paste0(outfile.base,'_extra_community_of_recipients_props.pdf'), w=5, h=5)	
	z	<- z[, list(N=length(TR_RID)), by=c('TR_COMM_TYPE','FROM_OUTSIDE')]
	tmp	<- dcast.data.table(subset(z, TR_COMM_TYPE!='trading'), FROM_OUTSIDE~TR_COMM_TYPE, value.var='N')
	zz	<- as.matrix(tmp[, 2:3, with=FALSE])
	rownames(zz)	<- tmp[[1]]
	fisher.test(zz)
	#STAGE 1
	#p-value = 0.05827
	#alternative hypothesis: true odds ratio is not equal to 1
	#95 percent confidence interval:
	#		0.9429887 3.1200202
	#sample estimates:
	#		odds ratio 
	#1.719531
	
	#STAGE 2
	#p-value = 0.003303
	#alternative hypothesis: true odds ratio is not equal to 1
	#95 percent confidence interval:
	# 1.346652 6.087743
	#sample estimates:
	#odds ratio 
	#  2.844896 
	
	
	#
	#	geography transmitters from outside community by gender of recipients
	z	<- copy(rtr2)
	z[, FROM_OUTSIDE:= factor(REC_COMM_NUM_A!=TR_COMM_NUM_A, levels=c(TRUE,FALSE), labels=c('outside community','same community'))]
	set(z, NULL, 'REC_COMM_TYPE', z[, factor(REC_COMM_TYPE)])
	set(z, NULL, 'REC_SEX', z[, paste0('gender recipient: ',REC_SEX)])	
	ggplot(z, aes(x=REC_COMM_TYPE, fill=FROM_OUTSIDE)) + geom_bar(position='fill') + 
			theme_bw() + 
			scale_y_continuous(labels=scales::percent, expand=c(0,0), breaks=seq(0,1,0.2)) +
			labs(x='\ncommunity type of recipient', y='transmitters\n', fill='transmitter from') +
			facet_grid(~REC_SEX)
	ggsave(file=paste0(outfile.base,'_extra_community_bygender_props.pdf'), w=7, h=5)	
}


RakaiFull.phylogeography.170522.community.anonymize<- function()
{
	#
	#	anonymize IDs
	#
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv'	
	dc		<- data.table(	COMM_NUM_RAW=	c("1","2","3","4","5","6","7","8","9","14","15","16","18","19","22","23","24","25","29","32","33","34","35","36","38","40","44","45","46","51","52","53","54","55", "56","57","58","59","60","61","62","65","67","74","77","81","84","89","94","95","103","106","107","108","109","120","177", "183", "256", "370","391","401","451", "468","602", "754", "755", "760", "770","771","772","773","774","776"),
			COMM_TYPE=		c("T","A","A","T","A","A","A","A","A", "A", "A", "T", "A", "A", "T", "A", "T", "A", "A", "A", "T", "A", "A", "A", "F", "A", "A", "A", "A", "T", "A", "A", "A", "A",  "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",  "A","A",   "A",  "T",  "A",  "A",  "A",  "A",   "A",   "A",   "A",  "A", "A",  "A",    "A",  "A",  "A",    "A",   "A", "F",  "F",  "A",  "A",  "F",   "T"))
	set(dc, NULL, 'COMM_TYPE', dc[, as.character(factor(COMM_TYPE, levels=c('A','T','F'), labels=c('agrarian','trading','fisherfolk')))])		
	set(dc, NULL, 'COMM_NUM', dc[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',COMM_NUM_RAW))))])
	
	seed	<- 42L 
	set.seed(seed)
	tmp		<- unique(dc, by='COMM_NUM')
	tmp		<- tmp[, list( 	COMM_NUM=COMM_NUM,
					COMM_NUM_A=paste0(substring(COMM_TYPE,1,1), letters[round(runif(length(COMM_NUM), min=1, max=24), d=0)], letters[round(runif(length(COMM_NUM), min=1, max=24), d=0)])
			), by='COMM_TYPE']
	dc		<- merge(dc, tmp, by=c('COMM_TYPE','COMM_NUM'))
	setkey(dc, COMM_NUM_A)
	
	set(dc, dc[, which(COMM_NUM=='89')], 'COMM_NUM_A', 'aop')
	set(dc, dc[, which(COMM_NUM=='45')], 'COMM_NUM_A', 'arb')
	set(dc, dc[, which(COMM_NUM=='755')], 'COMM_NUM_A', 'ate')
	set(dc, dc[, which(COMM_NUM=='84')], 'COMM_NUM_A', 'awb')
	stopifnot( !nrow(unique(dc, by='COMM_NUM')[which(duplicated(unique(dc, by='COMM_NUM'),by='COMM_NUM_A')),]) )
	write.csv(dc, row.names=FALSE, file=outfile)	
	#
	#	anonymize locations
	#
	require(ggmap)
	zm					<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	infile.comms		<- '~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv'
	infile.loc.a		<- '~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/community_spatialLoc.csv'
	infile.loc			<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/community_geography.rda"
	outfile.base		<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_"
	#	original locations	
	dc		<- as.data.table(read.csv(infile.comms))
	zc		<- as.data.table(read.csv(infile.loc.a))
	load(infile.loc)
	comgps	<- as.data.table(comgps)
	setnames(comgps, 'COMM_NUM', 'COMM_NUM_RAW')
	set(comgps, NULL, 'COMM_NUM_RAW', comgps[, as.integer(as.character(COMM_NUM_RAW))])
	dc		<- merge(dc, comgps, by='COMM_NUM_RAW',all.x=1)
	set.seed(42L)
	tmp	<- unique(subset(dc, select=c(COMM_NUM, longitude, latitude)), by='COMM_NUM')
	tmp[, LONG_A:= rnorm(nrow(tmp), mean=longitude, sd=0.1)]
	tmp[, LAT_A:= rnorm(nrow(tmp), mean=latitude, sd=0.1)]
	set(tmp, tmp[, which(COMM_NUM=='770')], 'LONG_A', 31.8)
	set(tmp, tmp[, which(COMM_NUM=='770')], 'LAT_A', -0.59)
	set(tmp, tmp[, which(COMM_NUM=='771')], 'LONG_A', 31.75)
	set(tmp, tmp[, which(COMM_NUM=='771')], 'LAT_A', -0.76)	
	set(tmp, tmp[, which(COMM_NUM=='38')], 'LONG_A', 31.8)
	set(tmp, tmp[, which(COMM_NUM=='38')], 'LAT_A', -1.05)
	set(tmp, tmp[, which(COMM_NUM=='772')], 'LONG_A', 31.5)
	set(tmp, tmp[, which(COMM_NUM=='772')], 'LAT_A', -1)	
	set(tmp, tmp[, which(COMM_NUM=='370')], 'LONG_A', 31.59)
	set(tmp, tmp[, which(COMM_NUM=='370')], 'LAT_A', -1.04)	
	set(tmp, tmp[, which(COMM_NUM=='23')], 'LONG_A', 31.65)
	set(tmp, tmp[, which(COMM_NUM=='23')], 'LAT_A', -0.9)
	set(tmp, tmp[, which(COMM_NUM=='16m')], 'LONG_A', 31.4)
	set(tmp, tmp[, which(COMM_NUM=='16m')], 'LAT_A', -1.1)
	set(tmp, tmp[, which(COMM_NUM=='51m')], 'LONG_A', 31.66)
	set(tmp, tmp[, which(COMM_NUM=='51m')], 'LAT_A', -0.43)
	set(tmp, tmp[, which(COMM_NUM=='40')], 'LONG_A', 31.51)
	set(tmp, tmp[, which(COMM_NUM=='40')], 'LAT_A', -0.77)
	set(tmp, tmp[, which(COMM_NUM=='56')], 'LONG_A', 31.48)
	set(tmp, tmp[, which(COMM_NUM=='56')], 'LAT_A', -0.85)
	set(tmp, tmp[, which(COMM_NUM=='108')], 'LONG_A', 31.3)
	set(tmp, tmp[, which(COMM_NUM=='108')], 'LAT_A', -0.8)
	set(tmp, tmp[, which(COMM_NUM=='34')], 'LONG_A', 31.4)
	set(tmp, tmp[, which(COMM_NUM=='34')], 'LAT_A', -0.55)
	set(tmp, tmp[, which(COMM_NUM=='57')], 'LONG_A', 31.45)
	set(tmp, tmp[, which(COMM_NUM=='57')], 'LAT_A', -0.39)
	set(tmp, tmp[, which(COMM_NUM=='58')], 'LONG_A', 31.52)
	set(tmp, tmp[, which(COMM_NUM=='58')], 'LAT_A', -0.74)
	set(tmp, tmp[, which(COMM_NUM=='74')], 'LONG_A', 31.68)
	set(tmp, tmp[, which(COMM_NUM=='74')], 'LAT_A', -0.40)
	set(tmp, tmp[, which(COMM_NUM=='754')], 'LONG_A', 31.5)
	set(tmp, tmp[, which(COMM_NUM=='754')], 'LAT_A', -0.8)
	set(tmp, tmp[, which(COMM_NUM=='5')], 'LONG_A', 31.63)
	set(tmp, tmp[, which(COMM_NUM=='5')], 'LAT_A', -0.5)
	set(tmp, tmp[, which(COMM_NUM=='6')], 'LONG_A', 31.41)
	set(tmp, tmp[, which(COMM_NUM=='6')], 'LAT_A', -0.48)
	set(tmp, tmp[, which(COMM_NUM=='62')], 'LONG_A', 31.51)
	set(tmp, tmp[, which(COMM_NUM=='62')], 'LAT_A', -0.46)
	set(tmp, tmp[, which(COMM_NUM=='89')], 'LONG_A', 31.55)
	set(tmp, tmp[, which(COMM_NUM=='89')], 'LAT_A', -0.62)
	set(tmp, tmp[, which(COMM_NUM=='391')], 'LONG_A', 31.53)
	set(tmp, tmp[, which(COMM_NUM=='391')], 'LAT_A', -0.65)	
	dc	<- merge(dc, subset(tmp, select=c(COMM_NUM, LONG_A, LAT_A)), by='COMM_NUM')
	write.csv(dc, row.names=FALSE, file=infile.comms)
	ggmap(zm) +
			geom_point(data=unique(dc, by='COMM_NUM'), aes(x=LONG_A, y=LAT_A, pch=COMM_TYPE, colour=COMM_TYPE), size=8, alpha=0.8) +
			geom_text(data=unique(dc, by='COMM_NUM'), aes(x=LONG_A, y=LAT_A, label=COMM_NUM_A), nudge_x=0, nudge_y=0, size=3, colour='black')
	ggsave(file=paste0(outfile.base,'_hubs_comm_locations_anonymized.pdf'), w=10, h=10)
	
}