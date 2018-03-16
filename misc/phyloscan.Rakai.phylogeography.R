
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
	
	#	locations with number sequence samples
	ds	<- subset(rsm, MIN_PNG_OUTPUT>0, select=c(COMM_NUM, COMM_TYPE, COMM_NUM_A, ELIGIBLE_AVG, PARTICIPATED_AVG, HIV, MIN_PNG_OUTPUT))
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
	
	#
	#	plot number of observed recipients and observed transmitters
	ggmap(zm) +
			geom_point(data=dc.rec, aes(x=longitude, y=latitude, pch=COMM_TYPE, colour=REC_OBS_C, size=REC_OBS), alpha=0.8) +
			geom_text(data=dc.rec, aes(x=longitude, y=latitude, label=REC_OBS), nudge_x=0, nudge_y=0, size=3, colour='black') +			
			scale_colour_brewer(palette="YlOrRd") +
			scale_size(breaks=c(1,2,5,10,20,50), range=c(5,15))+						
			theme(legend.position='bottom', legend.box = "vertical") +
			guides(size='none', colour='none', pch=guide_legend(override.aes=list(size=7))) +
			labs(	x='', y='', pch="community\ntype", colour="phylogenetically observed recipients")
	ggsave(file=paste0(outfile.base,'_comm_recipientnumbers_on_map.pdf'), w=7, h=7)
	ggmap(zm) +
			geom_point(data=dc.tr, aes(x=longitude, y=latitude, pch=COMM_TYPE, colour=TR_OBS_C, size=TR_OBS), alpha=0.8) +
			geom_text(data=dc.tr, aes(x=longitude, y=latitude, label=TR_OBS), nudge_x=0, nudge_y=0, size=3, colour='black') +			
			scale_colour_brewer(palette="YlOrRd") +
			scale_size(breaks=c(1,2,5,10,20,50), range=c(5,15))+						
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
	arrow		<- arrow(length=unit(0.04, "npc"), type="closed")
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
			geom_point(data=subset(dc,  REC_COMM_NUM_A==TR_COMM_NUM_A), aes(x=TR_X, y=TR_Y, size=node.size*EDGE_LABEL)) +
			geom_curve(data=subset(dc, substr(REC_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A!=TR_COMM_NUM_A & FLOW>0), aes(x=TR_X, xend=REC_X, y=TR_Y, yend=REC_Y, size=edge.size*FLOW), curvature=curvature, arrow=arrow, lineend="butt") +
			geom_text(data=subset(dc, substr(REC_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A!=TR_COMM_NUM_A & FLOW>0), aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL), size=label.size) +
			coord_cartesian() +
			scale_size_identity() 
	ggsave(file=paste0(outfile.base,'_flows_intofishing_on_map.pdf'), w=7, h=7)		
	ggmap(zm) +		
			geom_point(data=subset(dc,  substr(TR_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A==TR_COMM_NUM_A), aes(x=TR_X, y=TR_Y, size=node.size*EDGE_LABEL), colour='DarkOrange') +
			geom_curve(data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A!=TR_COMM_NUM_A & FLOW>0), aes(x=TR_X, xend=REC_X, y=TR_Y, yend=REC_Y, size=edge.size*FLOW), colour='DarkOrange', curvature=curvature, arrow=arrow, lineend="butt") +
			#geom_text(data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A!=TR_COMM_NUM_A & FLOW>0), aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL), size=label.size) +
			coord_cartesian() +
			scale_size_identity() +
			labs(x='', y='') +
			theme(legend.position='bottom', legend.box = "vertical")
	ggsave(file=paste0(outfile.base,'_flows_fromfishing_on_map.pdf'), w=7, h=7)					
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