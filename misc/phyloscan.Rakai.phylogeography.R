RakaiFull.phylogeography.180521.prevalence.gender<- function()
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
	require(rethinking)	# STAN wrapper
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_withmetadata.rda"
	outfile.base			<- gsub('_withmetadata.rda','',infile)	
	load(infile)	
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	rtpdm[, PAIR_COMM_TYPE:= FEMALE_COMM_TYPE]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_TYPE!=MALE_COMM_TYPE)], 'PAIR_COMM_TYPE', 'mixed')
	set(rtpdm, rtpdm[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(rtpdm, rtpdm[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')
	rtpdm[, COUPLE2:= factor(COUPLE=='no couple', levels=c(TRUE,FALSE), labels=c('no couple','couple'))]
	rtpdm[, SAMEHH:= factor(FEMALE_HH_NUM==MALE_HH_NUM, levels=c(TRUE,FALSE), labels=c('same hh','different hh'))]	
	rtpdm[, PAIR_COMM:= MALE_COMM_NUM_A]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_NUM_A!=MALE_COMM_NUM_A)], 'PAIR_COMM', 'mixed')
	rtpdm[, MALE_SEXP1OUT2:= factor(MALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(MALE_SEXP1OUT=='Unknown')], 'MALE_SEXP1OUT2', 'Unknown')
	rtpdm[, FEMALE_SEXP1OUT2:= factor(FEMALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(FEMALE_SEXP1OUT=='Unknown')], 'FEMALE_SEXP1OUT2', 'Unknown')
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/community_hivprev_byGenderStudyVisit.rda"
	load(infile)
	df	<- as.data.table(melt(female.negative, varnames=c('COMM_NUM', 'VISIT'), value.name='FEMALE_NEG'))	
	tmp	<- as.data.table(melt(male.negative, varnames=c('COMM_NUM', 'VISIT'), value.name='MALE_NEG'))	
	df	<- merge(df, tmp, by=c('COMM_NUM', 'VISIT'), all=TRUE)	
	tmp	<- as.data.table(melt(female.positive, varnames=c('COMM_NUM', 'VISIT'), value.name='FEMALE_POS'))
	df	<- merge(df, tmp, by=c('COMM_NUM', 'VISIT'), all=TRUE)
	tmp	<- as.data.table(melt(male.positive, varnames=c('COMM_NUM', 'VISIT'), value.name='MALE_POS'))
	df	<- merge(df, tmp, by=c('COMM_NUM', 'VISIT'), all=TRUE)
	#	add geo-locations
	load("~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/community_geography.rda")
	comgps	<- as.data.table(comgps)
	set(comgps, NULL, 'COMM_NUM', comgps[, as.integer(as.character(COMM_NUM))])
	df	<- merge(df, comgps, by='COMM_NUM', all.x=TRUE)	
	set(df, NULL, 'COMM_NUM', df[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',as.character(COMM_NUM)))))])
	for(x in c('FEMALE_NEG','MALE_NEG','FEMALE_POS','MALE_POS'))
		set(df, which(is.na(df[[x]])), x, 0L)
	df	<- subset(df, FEMALE_NEG!=0 | MALE_NEG!=0 | FEMALE_POS!=0 | MALE_POS!=0)
	#	sum by merged communities that are essentially the same
	df	<- df[, list(FEMALE_NEG=sum(FEMALE_NEG), MALE_NEG=sum(MALE_NEG), FEMALE_POS=sum(FEMALE_POS), MALE_POS=sum(MALE_POS), longitude=mean(longitude), latitude=mean(latitude)), by= c('COMM_NUM','VISIT')]
	#	add anonymized ID
	tmp	<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	tmp	<- unique(subset(tmp, select=c(COMM_NUM, COMM_NUM_A, COMM_TYPE)))
	df	<- merge(df, tmp, by='COMM_NUM', all.x=TRUE)
	#	subset(df, is.na(COMM_NUM_A))	#comms with missing anonymized IDs are some historic ones, not relevant
	df	<- subset(df, !is.na(COMM_NUM_A))
	
	#	select relevant communities
	if(0)
	{
		tmp	<- unique(subset(rtpdm, select=c(MALE_COMM_NUM,MALE_COMM_TYPE)))
		setnames(tmp, c('MALE_COMM_NUM','MALE_COMM_TYPE'), c('COMM_NUM','COMM_TYPE'))
		tmp2<- unique(subset(rtpdm, select=c(FEMALE_COMM_NUM,FEMALE_COMM_TYPE)))
		setnames(tmp2, c('FEMALE_COMM_NUM','FEMALE_COMM_TYPE'), c('COMM_NUM','COMM_TYPE'))
		tmp	<- unique(rbind(tmp, tmp2))
		df	<- merge(tmp, df, by='COMM_NUM')		
	}
	df[, FEMALE:= FEMALE_NEG+FEMALE_POS]
	df[, MALE:= MALE_NEG+MALE_POS]
	df[, POS:= MALE_POS+FEMALE_POS]
	df[, NEG:= MALE_NEG+FEMALE_NEG]
	#	new clean community numbers
	dg	<- subset(df, VISIT%in%c(15,15.1,16))
	tmp	<- unique(subset(dg, select=c(COMM_NUM,COMM_TYPE)))
	tmp[, COMM_NUM2:= seq_len(nrow(tmp))]
	tmp[, COMM_TYPE2:= as.integer(as.character(factor(COMM_TYPE,levels=c('agrarian','trading','fisherfolk'), labels=c('1','2','3'))))]
	dg	<- merge(dg, tmp, by=c('COMM_NUM','COMM_TYPE'))
	
	if(0)
	{
		#	any change in time in gender specific prevalence?
		#	not really
		dg	<- subset(df, VISIT%in%c(14,15,15.1,16,17))
		ggplot(dg, aes(x=factor(VISIT), group=COMM_NUM)) +
				geom_line(aes(y= MALE_POS/MALE), colour='blue') +
				geom_line(aes(y= FEMALE_POS/FEMALE), colour='hotpink2') +
				geom_point(aes(y= MALE_POS/MALE), colour='blue') +
				geom_point(aes(y= FEMALE_POS/FEMALE), colour='hotpink2') +
				theme_bw() +
				facet_wrap(~COMM_TYPE+COMM_NUM, ncol=4) +
				labs(x='\nvisit', y='gender specific HIV prevalence estimate\n')
		ggsave(file=paste0(outfile.base,'_trmMF_prevalenceratios.pdf'), w=9, h=9)		
	}
	
	#	estimate prevalences and prevalence ratio for each community		
	mpr.1 	<- map2stan(
			alist(
					FEMALE_POS ~ dbinom(FEMALE, seropos_f),
					MALE_POS ~ dbinom(MALE, seropos_m),
					logit(seropos_f) <- afc[COMM_NUM2],
					logit(seropos_m) <- amc[COMM_NUM2],					
					afc[COMM_NUM2] ~ dnorm(0, 10),
					amc[COMM_NUM2] ~ dnorm(0, 10)										
			),
			data=as.data.frame(dg), 
			start=list(	afc=rep(0,length(unique(dg$COMM_NUM2))), amc=rep(0,length(unique(dg$COMM_NUM2)))),
			warmup=5e2, iter=5e3, chains=1, cores=4
	)	
	post	<- extract.samples(mpr.1)			
	dgg		<- unique(subset(dg, select=c(COMM_NUM, COMM_TYPE, COMM_NUM2, COMM_NUM_A, longitude, latitude)))
	tmp		<- dgg[, 
			list(	STAT=c('M','CL','IL','IU','CU'),
					PF= as.numeric(quantile(logistic(post$afc[, COMM_NUM2]), prob=c(0.5,0.025,0.25,0.75,0.975))),
					PM= as.numeric(quantile(logistic(post$amc[, COMM_NUM2]), prob=c(0.5,0.025,0.25,0.75,0.975))),
					RFM= as.numeric(quantile(logistic(post$afc[, COMM_NUM2])/logistic(post$amc[, COMM_NUM2]), prob=c(0.5,0.025,0.25,0.75,0.975)))
			), 
			by='COMM_NUM2']
	tmp		<- melt(tmp, id.vars=c('COMM_NUM2','STAT'))	
	tmp		<- dcast.data.table(tmp, variable+COMM_NUM2~STAT, value.var='value')
	dgg		<- merge(dgg, tmp, by='COMM_NUM2')
	#
	#	 plot
	#
	tmp		<- subset(dgg, variable!='RFM')
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland communities')
	set(tmp, tmp[, which(COMM_TYPE=='fisherfolk')], 'COMM_TYPE', 'fishing sites')
	set(tmp, tmp[, which(variable=='PM')], 'variable', 'men')
	set(tmp, tmp[, which(variable=='PF')], 'variable', 'women')
	tmp2	<- subset(tmp, variable=='men', c(COMM_NUM, COMM_NUM_A, COMM_TYPE, M))	
	tmp2	<- tmp2[order(COMM_TYPE, M),]
	tmp2[, DUMMY:= seq_len(nrow(tmp2))]
	set(tmp2, NULL, 'COMM_NUM_A', tmp2[, factor(DUMMY, levels=DUMMY, labels=COMM_NUM_A)])
	tmp2	<- subset(tmp2, select=c(COMM_NUM,COMM_NUM_A))
	tmp[, COMM_NUM_A:=NULL]
	tmp		<- merge(tmp, tmp2, by='COMM_NUM')
	ggplot(tmp, aes(x=COMM_NUM_A)) +
			geom_boxplot(aes(middle=M, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=variable), stat='identity') +
			theme_bw() + 
			facet_grid(~COMM_TYPE, space='free', scale='free') +
			scale_fill_manual(values=c('women'='hotpink2', 'men'='deepskyblue')) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,0.6)) +
			labs(x='\ncommunity', y='HIV-1 prevalence\namong study participants\n', fill='gender')
	ggsave(file=paste0(outfile.base,'_trmMF_estimatedprevalence.pdf'), w=12, h=6)
	
	ggplot(subset(dgg, variable!='RFM'), aes(x=COMM_NUM3)) +
			geom_boxplot(aes(middle=M, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=variable), stat='identity') +
			theme_bw() + 
			#facet_grid(~COMM_TYPE, space='free', scale='free') +
			scale_fill_manual(values=c('PF'='hotpink2', 'PM'='deepskyblue')) +
			scale_y_continuous(labels=scales:::percent) +
			labs(x='\ncommunity', y='HIV prevalence\nvisits 15, 15.1, 16\n', fill='gender')
	ggsave(file=paste0(outfile.base,'_trmMF_estimatedprevalence_orderedbyFMratio.pdf'), w=12, h=6)
	save(dgg, file=paste0(outfile.base,'_trmMF_estimatedFMprevalence.rda'))
	#
	#	plot on map
	#	
	require(ggmap)
	#zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	style	<- "feature:road|color:0x17202A&style=feature:water|color:0x677996&style=feature:landscape.natural|color:0xedecda&style=feature:administrative|visibility=off"
	zm		<- get_googlemap(c(lon=31.65, lat=-0.66), scale=2, size=c(550,550), zoom=10, maptype="road", style=style)
	#	plot number of observed recipients and observed transmitters
	ggmap(zm) +
			geom_point(data=subset(dgg, variable=='RFM'), aes(x=longitude, y=latitude, pch=COMM_TYPE, colour=M), size=7) +
			geom_text(data=subset(dgg, variable=='RFM'), aes(x=longitude, y=latitude, label=COMM_NUM), nudge_x=0, nudge_y=0, size=3, colour='black') + 
			scale_colour_gradient2(trans='log', breaks=c(1,1.25,1.5,2,2.5,3), low="deepskyblue", mid="orange", high='red', midpoint=log(2.2)) +
			labs(x='\nlongitude',y='latitude\n',colour='female to male\nprevalence ratio', pch='community\ntype')
	ggsave(file=paste0(outfile.base,'trmMF_prevalenceratio_on_map.pdf'), w=7, h=7)	
	ggmap(zm) +
			geom_point(data=subset(dgg, variable=='PF'), aes(x=longitude, y=latitude, pch=factor(COMM_TYPE=='fisherfolk', levels=c(TRUE,FALSE),labels=c('fishing\nsite','inland\ncommunity')), colour=M), size=9) +
			scale_colour_gradient2(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6), lim=c(0.04,0.56), labels=scales:::percent, low="cyan", mid='darkorchid2', high='firebrick1', midpoint=0.35) +
			scale_shape_manual(values=c('fishing\nsite'=17,'inland\ncommunity'=19)) +
			labs(x='\nlongitude',y='latitude\n',colour='female\nHIV-1 prevalence', pch='')
	ggsave(file=paste0(outfile.base,'trmMF_prevalenceF_on_map.pdf'), w=7, h=7, useDingbats=FALSE)
	ggmap(zm) +
			geom_point(data=subset(dgg, variable=='PM'), aes(x=longitude, y=latitude, pch=factor(COMM_TYPE=='fisherfolk', levels=c(TRUE,FALSE),labels=c('fishing\nsite','inland\ncommunity')), colour=M), size=9) +			 
			scale_colour_gradient2(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6), lim=c(0.04,0.56), labels=scales:::percent, low="cyan", mid='darkorchid2', high='firebrick1', midpoint=0.35) +
			scale_shape_manual(values=c('fishing\nsite'=17,'inland\ncommunity'=19)) +
			labs(x='\nlongitude',y='latitude\n',colour='male\nHIV-1 prevalence', pch='')
	ggsave(file=paste0(outfile.base,'trmMF_prevalenceM_on_map.pdf'), w=7, h=7, useDingbats=FALSE)	
}



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
	
	set.seed(123)
	long <- rnorm(50, sd=100)
	lat <- rnorm(50, sd=50)
	d <- data.frame(long=long, lat=lat)
	d <- with(d, d[abs(long) < 150 & abs(lat) < 70,])
	n <- nrow(d)
	d$region <- factor(1:n)
	d$A <- abs(rnorm(n, sd=1))
	d$B <- abs(rnorm(n, sd=2))
	d$C <- abs(rnorm(n, sd=3))
	d$D <- abs(rnorm(n, sd=4))
	d$radius <- 6 * abs(rnorm(n))
	head(d)
	library(ggplot2)
	library(scatterpie)
	
	world <- map_data('world')
	p <- ggplot(world, aes(long, lat)) +
			geom_map(map=world, aes(map_id=region), fill=NA, color="black") +
			coord_quickmap()
	p + geom_scatterpie(aes(x=long, y=lat, group=region, r=radius),
					data=d, cols=LETTERS[1:4], color=NA, alpha=.8) +
			geom_scatterpie_legend(d$radius, x=-160, y=-55)
	
	infile	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/popviraemia_data.rda"
	load(infile)
	df[, SLART_YES:= SLART_YES_M+SLART_YES_F]
	
}

RakaiFull.phylogeography.180618.samplinglocations<- function()
{	
	library(ggplot2)
	library(scatterpie)
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
	#	get map
	#
	style	<- "feature:road|color:0x17202A&style=feature:water|color:0x677996&style=feature:landscape.natural|color:0xedecda&style=feature:administrative|visibility=off"
	zm		<- get_googlemap(c(lon=31.65, lat=-0.66), scale=2, size=c(550,550), zoom=10, maptype="road", style=style)
	zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))			
	zc		<- unique(zc, by='COMM_NUM')
	#
	#	get self reported ART, POS
	#
	infile			<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender_age.rda'		
	load(infile)
	#	get those infected before round 17 and who participated in at least one round 15-16
	df		<- subset(de, HIV_1517==1, select=c(COMM_NUM, RID, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT))
	df		<- df[, list(MIN_PNG_OUTPUT=max(MIN_PNG_OUTPUT)), by=c('COMM_NUM','RID','HIV_1517','SELFREPORTART_AT_FIRST_VISIT')]	
	df		<- df[, list(	POS=sum(HIV_1517),
							SLART_YES=sum(SELFREPORTART_AT_FIRST_VISIT),
							DEEP_SEQ_1516=sum(MIN_PNG_OUTPUT)
					), by='COMM_NUM']
	df[, SLART_NO_NOSEQ:= POS-DEEP_SEQ_1516-SLART_YES]
	df		<- merge(df, zc, by='COMM_NUM')		
							
	df[, POS2:= POS/1000]
	df[, POS3:= POS^(1/3)/300]
	ggmap(zm) +
		geom_scatterpie(aes(x=longitude, y=latitude, group=COMM_NUM, r=POS3),
						data=df, cols=c('SLART_YES','SLART_NO_NOSEQ','DEEP_SEQ_1516')) +		
		scale_fill_manual(values=c('SLART_NO_NOSEQ'='turquoise2','DEEP_SEQ_1516'='firebrick2','SLART_YES'='grey50')) +
		geom_scatterpie_legend(df$POS3, x=31.4, y=-0.4, n=5, labeller=function(radius) round((radius*400)^3) )
	ggsave(file=paste0(outfile.base,'_comm_seqcov_map.pdf'), w=7, h=5)
}

RakaiCirc.epi.get.info.170208<- function()
{
	hivc.db.Date2numeric<- function( x )
	{
		if(!class(x)%in%c('Date','character'))	return( x )
		x	<- as.POSIXlt(x)
		tmp	<- x$year + 1900
		x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
		x	
	}
	#
	infile				<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/alldat_r15or17_witharv.rda"
	load(infile)
	ra		<- as.data.table(alldat)
	ra		<- subset(ra, select=c(RCCS_studyid,REGION,COMM_NUM,HH_NUM,SEX,AGEYRS,visdate,visit,lastNegDate,hiv,firstPosDate,eversex, evermarr, currmarr, polymar, sexpever, sexp1yr, sexp1out, sexgift, sexyear, religion, educate, educyrs, edcat, occ, occ2, arvmed))
	#	a bit of clean up 	
	setnames(ra, colnames(ra), gsub('\\.','_',toupper(colnames(ra))))	
	set(ra, NULL, 'VISDATE', ra[, as.Date(VISDATE, format='%d/%m/%Y')])	
	for(x in colnames(ra))
		if(class(ra[[x]])=='Date')
			set(ra, NULL, x, hivc.db.Date2numeric(ra[[x]]))
	#	
	wdir				<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/circumcision"	
	infile				<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData_v2.rda"
	load(infile)
	#	a bit of clean up 
	rd		<- as.data.table(rccsData)
	setnames(rd, colnames(rd), gsub('\\.','_',toupper(colnames(rd))))	
	for(x in colnames(rd))
		if(class(rd[[x]])=='Date')
			set(rd, NULL, x, hivc.db.Date2numeric(rd[[x]]))
	rh		<- as.data.table(rccsHistory)
	setnames(rh, colnames(rh), gsub('\\.','_',toupper(colnames(rh))))
	rn		<- as.data.table(neuroData)
	setnames(rn, colnames(rn), gsub('\\.','_',toupper(colnames(rn))))
	for(x in colnames(rn))
		if(class(rn[[x]])=='Date')
			set(rn, NULL, x, hivc.db.Date2numeric(rn[[x]]))
	#rd[, table(VISIT)]
	#	make shorter
	setnames(ra, 'RCCS_STUDYID', 'RID')
	setnames(rd, 'RCCS_STUDYID', 'RID')
	setnames(rd, 'PANGEA_ID', 'PID')
	setnames(rd, 'STUDYID', 'SID')
	setnames(rh, 'RCCS_STUDYID', 'RID')	
	setnames(rn, 'STUDYID', 'RID')
	setnames(rn, 'PANGEA_ID', 'PID')
	#	get neuro into format of rd
	setnames(rn, c('SAMPLEDATE','GENDER','CD4','VIRALLOAD'), c('DATE','SEX','RECENTCD4','RECENTVL'))
	set(rn, NULL, c('RECENTVLDATE','RECENTCD4DATE'), rn[, DATE])
	set(rn, NULL, c('TIMESINCEVL','TIMESINCECD4'), 0)
	#	data checks
	setkey(rh, VISIT, RID)
	stopifnot(nrow(rh)==nrow(unique(rh)))
	#	define circumcision	
	set(rh, rh[, which(!CIRCUM%in%c(1,2))], 'CIRCUM', NA_integer_)
	set(rh, NULL, 'CIRCUM', rh[, factor(CIRCUM, levels=c(1,2), labels=c('Y','N'))])
	#	define sexual relationships
	set(rh, rh[, which(is.na(RLTN1))],'RLTN1',99L)
	set(rh, rh[, which(is.na(RLTN2))],'RLTN2',99L)
	set(rh, rh[, which(is.na(RLTN3))],'RLTN3',99L)
	set(rh, rh[, which(is.na(RLTN4))],'RLTN4',99L)
	setnames(rh, c('RLTN1','RLTN2','RLTN3','RLTN4'), c('RLTN1_','RLTN2_','RLTN3_','RLTN4_'))
	warning('undocumented RLTN codes 15 16 88 -> set to Unknown')
	tmp		<- as.data.table(matrix(data=c(	'1','Current wife (at the time)',
							'2','Current consensual partner (at the time)',
							'3','Former wife/consensual partner',
							'4','Girlfriend',
							'5','Occasional or casual friend',
							'6','Visitor (incl. wedding/funeral)',
							'7','Stranger',
							'8','Workmate',
							'9','Boss/work supervisor',
							'10','Employee',
							'11','Fellow student',
							'12','Sugar mummy',
							'13','Relative other than spouse',
							'14','Other non relative',
							'15','Unknown',
							'16','Unknown',
							'88','Unknown',
							'0','Unknown',
							'97','Unknown',
							'98','Unknown',
							'99','Unknown'
					), ncol=2, byrow=TRUE))
	setnames(tmp, colnames(tmp), c('RLTN1_','RLTN1'))
	set(tmp, NULL,'RLTN1_',tmp[, as.integer(RLTN1_)])
	rh		<- merge(rh,tmp,by='RLTN1_')
	setnames(tmp, c('RLTN1_','RLTN1'), c('RLTN2_','RLTN2'))	
	rh		<- merge(rh,tmp,by='RLTN2_')
	setnames(tmp, c('RLTN2_','RLTN2'), c('RLTN3_','RLTN3'))	
	rh		<- merge(rh,tmp,by='RLTN3_')
	setnames(tmp, c('RLTN3_','RLTN3'), c('RLTN4_','RLTN4'))	
	rh		<- merge(rh,tmp,by='RLTN4_')
	set(rh, NULL, c('RLTN1_','RLTN2_','RLTN3_','RLTN4_'), NULL)	
	tmp		<- rh[, which(RLTN1=='Unknown' & RLTN2!='Unknown')]
	set(rh, tmp, 'RLTN1', rh[tmp, RLTN2])
	set(rh, tmp, 'RLTN2', 'Unknown')
	tmp		<- rh[, which(RLTN2=='Unknown' & RLTN3!='Unknown')]
	set(rh, tmp, 'RLTN2', rh[tmp, RLTN3])
	set(rh, tmp, 'RLTN3', 'Unknown')
	tmp		<- rh[, which(RLTN3=='Unknown' & RLTN4!='Unknown')]
	set(rh, tmp, 'RLTN3', rh[tmp, RLTN4])
	set(rh, tmp, 'RLTN4', 'Unknown')
	#	define number named sexual relations last year
	tmp		<- melt(rh, id.vars=c('RID','VISIT'), measure.vars=c('RLTN1','RLTN2','RLTN3','RLTN4'))
	tmp		<- tmp[, list(RLTN_NAMED= length(which(value!='Unknown'))), by=c('RID','VISIT')]
	rh		<- merge(rh, tmp, by=c('RID','VISIT'))
	#	define occupation
	#	set back OCCUP codes according to patient sheet
	setnames(rh, c('OCCUP1','OCCUP2','OCCUP3','OCCUP4'), c('OCCUP11','OCCUP12','OCCUP13','OCCUP14'))
	setnames(rh, c('OCC','OCC2'), c('OCCUP1','OCCUP2'))	
	set(rh, rh[, which(is.na(OCCUP1))],'OCCUP1',99L)
	set(rh, rh[, which(is.na(OCCUP2))],'OCCUP2',99L)
	set(rh, NULL, 'OCCUP1', rh[, as.integer(OCCUP1)])	
	setnames(rh, c('OCCUP1','OCCUP2'), c('OCCUP1_','OCCUP2_'))
	#	handle ra
	setnames(ra, c('OCC','OCC2'), c('OCCUP1','OCCUP2'))	
	set(ra, ra[, which(is.na(OCCUP1))],'OCCUP1',99L)
	set(ra, ra[, which(is.na(OCCUP2))],'OCCUP2',99L)
	set(ra, NULL, 'OCCUP1', ra[, as.integer(OCCUP1)])	
	setnames(ra, c('OCCUP1','OCCUP2'), c('OCCUP1_','OCCUP2_'))	
	tmp		<- as.data.table(matrix(data=c(	'1','Agriculture for home use/barter',
							'2','Agriculture for selling',
							'3','Housework in your own home',
							'4','Housekeeper',
							'5','Home brewing',
							'6','Government/clerical/teaching',
							'7','Fishing',
							'8','Student',
							'9','Military/police',
							'10','Shopkeeper',
							'11','Trading/vending',
							'12','Bar worker or owner',
							'13','Trucker',
							'14','Unemployed',
							'15','Other',
							'88','No additional occupation',
							'16','Medical worker',
							'17','Casual laborer',
							'18','Waitress/Waiter/restaurant owner',
							'19','Hair dresser/Salon owner',
							'20','Construction (brick maker, builder, porter, painter, roofing)',
							'21','Mechanic (automobiles, bicycles, electronics)',
							'22','Boda Boda',
							'23','Client/Sex worker',
							'0','Unknown',
							'98','Unknown',
							'99','Unknown'), ncol=2, byrow=TRUE))
	setnames(tmp, colnames(tmp), c('OCCUP1_','OCCUP1'))
	set(tmp, NULL,'OCCUP1_',tmp[, as.integer(OCCUP1_)])
	rh		<- merge(rh,tmp,by='OCCUP1_')
	setnames(tmp, c('OCCUP1_','OCCUP1'), c('OCCUP2_','OCCUP2') )
	rh		<- merge(rh,tmp,by='OCCUP2_')	
	set(rh, NULL, c('OCCUP1_','OCCUP2_'), NULL)
	setnames(tmp, c('OCCUP2_','OCCUP2'), c('OCCUP1_','OCCUP1') )
	#	handle ra
	ra		<- merge(ra,tmp,by='OCCUP1_')
	setnames(tmp, c('OCCUP1_','OCCUP1'), c('OCCUP2_','OCCUP2') )
	ra		<- merge(ra,tmp,by='OCCUP2_')	
	set(ra, NULL, c('OCCUP1_','OCCUP2_'), NULL)
	tmp		<- rh[, which(OCCUP1=='Unknown' & OCCUP2!='Unknown')]
	set(rh, tmp, 'OCCUP1', rh[tmp, OCCUP2])
	set(rh, tmp, 'OCCUP2', 'Unknown')
	tmp		<- ra[, which(OCCUP1=='Unknown' & OCCUP2!='Unknown')]
	set(ra, tmp, 'OCCUP1', ra[tmp, OCCUP2])
	set(ra, tmp, 'OCCUP2', 'Unknown')	
	#	condense OCCUP1 and OCCUP2
	for(x in c('OCCUP1','OCCUP2'))
	{
		set(rh, which(rh[[x]]%in%c('Boda Boda','Trucker')), x, 'Boda/Trucking')
		set(rh, which(rh[[x]]%in%c('Government/clerical/teaching','Military/police','Medical worker')), x, 'Government/clerical/related')
		set(rh, which(rh[[x]]%in%c('Trading/vending','Shopkeeper','Hair dresser/Salon owner')), x, 'Trading/Shopkeeper/Hair')
		set(rh, which(rh[[x]]%in%c('Agriculture for home use/barter','Agriculture for selling','Housekeeper','Housework in your own home','Home brewing')), x, 'Agro/House')
		set(rh, which(rh[[x]]%in%c('Waitress/Waiter/restaurant owner','Bar worker or owner')), x, 'Bar/waitress')
		set(rh, which(rh[[x]]%in%c('Casual laborer','Unemployed')), x, 'Casual laborer/unemployed')
		set(rh, which(rh[[x]]%in%c('Construction (brick maker, builder, porter, painter, roofing)','Mechanic (automobiles, bicycles, electronics)')), x, 'Construction/Mechanic')
		
		set(ra, which(ra[[x]]%in%c('Boda Boda','Trucker')), x, 'Boda/Trucking')
		set(ra, which(ra[[x]]%in%c('Government/clerical/teaching','Military/police','Medical worker')), x, 'Government/clerical/related')
		set(ra, which(ra[[x]]%in%c('Trading/vending','Shopkeeper','Hair dresser/Salon owner')), x, 'Trading/Shopkeeper/Hair')
		set(ra, which(ra[[x]]%in%c('Agriculture for home use/barter','Agriculture for selling','Housekeeper','Housework in your own home','Home brewing')), x, 'Agro/House')
		set(ra, which(ra[[x]]%in%c('Waitress/Waiter/restaurant owner','Bar worker or owner')), x, 'Bar/waitress')
		set(ra, which(ra[[x]]%in%c('Casual laborer','Unemployed')), x, 'Casual laborer/unemployed')
		set(ra, which(ra[[x]]%in%c('Construction (brick maker, builder, porter, painter, roofing)','Mechanic (automobiles, bicycles, electronics)')), x, 'Construction/Mechanic')
	}
	set(rh, rh[, which(OCCUP1=='No additional occupation' & OCCUP2=='Unknown')], 'OCCUP1', 'Unknown')
	#	refine OCCUP1 OCCUP2 when OCAT==Student
	set(rh, rh[, which(OCAT%in%c('Student'))], 'OCCUP1', 'Student')
	set(rh, rh[, which(OCAT%in%c('Student'))], 'OCCUP2', 'Student')
	#	define own OCCUP at time of RCCS visit
	rh[, OCCUP_OLLI:= OCCUP1]
	for(x in c('Student','Boda/Trucking','Bar/waitress','Client/Sex worker'))
		set(rh, rh[, which(OCCUP1==x | OCCUP2==x )],'OCCUP_OLLI', x)
	#	handle ra
	set(ra, ra[, which(OCCUP1=='No additional occupation' & OCCUP2=='Unknown')], 'OCCUP1', 'Unknown')
	ra[, OCCUP_OLLI:= OCCUP1]
	for(x in c('Student','Boda/Trucking','Bar/waitress','Client/Sex worker'))
		set(ra, ra[, which(OCCUP1==x | OCCUP2==x )],'OCCUP_OLLI', x)
	#	OK this is getting complicated: Kate is using all OCC codes to override OCCUP1 
	#	as deemed sensible
	#	just OCAT for simplicity
	#subset(rh, OCAT=='Bar/waitress' & OCCUP1!='Bar/waitress')
	#subset(rh, RID=='G030852')
	#	define SEXWORK
	set(rh, NULL, 'SEXWORK', rh[, as.character(factor(SEXWORK,levels=c(0,1),labels=c('N','Y')))])
	#	define SEXBAR (either work in bar or have sexual contact with barworker)
	set(rh, NULL, 'SEXBAR', rh[, as.character(factor(SEXBAR,levels=c(0,1),labels=c('N','Y')))])	
	#	extend MARSTAT for rh
	tmp		<- rh[, which((is.na(MARSTAT)|MARSTAT=='Never Married') & (RLTN1=='Current wife (at the time)'|RLTN2=='Current wife (at the time)'|RLTN3=='Current wife (at the time)'|RLTN4=='Current wife (at the time)'))]
	set(rh, tmp, 'EVERMARR', 1L)
	set(rh, tmp, 'CURRMARR', 1L)
	set(rh, tmp, 'MARSTAT', 'Monogamous')
	tmp		<- rh[, which(is.na(MARSTAT) & (RLTN1=='Former wife/consensual partner'|RLTN2=='Former wife/consensual partner'|RLTN3=='Former wife/consensual partner'|RLTN4=='Former wife/consensual partner'))]
	set(rh, tmp, 'EVERMARR', 1L)
	set(rh, tmp, 'PREVMAR', 1L)
	set(rh, tmp, 'MARSTAT', 'Previously Married')
	tmp		<- rh[, which(is.na(MARSTAT) & (!RLTN1%in%c('Unknown','Current consensual partner (at the time)')))]
	set(rh, tmp, 'EVERMARR', 0L)
	set(rh, tmp, 'PREVMAR', 0L)
	set(rh, tmp, 'CURRMARR', 0L)
	set(rh, tmp, 'MARSTAT', 'Never Married')
	set(rh, rh[, which(is.na(MARSTAT))], 'MARSTAT', 'Unknown')
	#	define MARSTAT for ra
	ra[, MARSTAT:='Unknown']
	set(ra, ra[, which(EVERMARR==2)], 'MARSTAT', 'Never Married')
	set(ra, ra[, which(EVERMARR==1 & CURRMARR==2)], 'MARSTAT', 'Previously Married')
	set(ra, ra[, which(EVERMARR==1 & CURRMARR==1)], 'MARSTAT', 'Monogamous')
	set(ra, ra[, which(POLYMAR>1 & !POLYMAR%in%c(97,98))], 'MARSTAT', 'Polygamous')
	#	define SEXYEAR
	set(rh, rh[, which(is.na(SEXYEAR))],'SEXYEAR', 99)
	set(rh, NULL, 'SEXYEAR', rh[, gsub('_[0-9]$','',as.character(factor(SEXYEAR, levels=c(0,1,2,8,99), labels=c('Unknown_1','Y','N','Unknown_2','Unknown_3'))))])
	stopifnot( !nrow(subset(rh, is.na(SEXYEAR))) )
	set(ra, ra[, which(is.na(SEXYEAR))],'SEXYEAR', 99)
	set(ra, NULL, 'SEXYEAR', ra[, gsub('_[0-9]$','',as.character(factor(SEXYEAR, levels=c(0,1,2,8,9,99), labels=c('Unknown_1','Y','N','Unknown_2','Unknown_3','Unknown_4'))))])
	stopifnot( !nrow(subset(ra, is.na(SEXYEAR))) )		
	#	define SEXP1YR 
	setnames(rh, c('SEXP1YR'), c('SEXP1YR_'))	
	rh[, SEXP1YR:= as.character(SEXP1YR_)]
	set(rh, rh[, which(SEXP1YR_==92)],'SEXP1YR','<3')
	set(rh, rh[, which(SEXP1YR_==93)],'SEXP1YR','3+')
	set(rh, rh[, which(SEXP1YR_%in%c(97,98,99) | is.na(SEXP1YR_))],'SEXP1YR','Unknown')
	set(rh, NULL, 'SEXP1YR_', NULL)
	setnames(ra, c('SEXP1YR'), c('SEXP1YR_'))	
	ra[, SEXP1YR:= as.character(SEXP1YR_)]
	set(ra, ra[, which(SEXP1YR_==92)],'SEXP1YR','<3')
	set(ra, ra[, which(SEXP1YR_==93)],'SEXP1YR','3+')
	set(ra, ra[, which(SEXP1YR_%in%c(97,98,99) | is.na(SEXP1YR_))],'SEXP1YR','Unknown')
	set(ra, NULL, 'SEXP1YR_', NULL)	
	#	revisit SEXP1YR based on relationships
	set(rh, rh[, which(SEXP1YR%in%c('Unknown') & RLTN_NAMED==1)], 'SEXP1YR', '1')
	set(rh, rh[, which(SEXP1YR%in%c('Unknown') & RLTN_NAMED==2)], 'SEXP1YR', '2')
	set(rh, rh[, which(SEXP1YR%in%c('Unknown') & RLTN_NAMED==3)], 'SEXP1YR', '3')
	set(rh, rh[, which(SEXP1YR%in%c('Unknown') & RLTN_NAMED==4)], 'SEXP1YR', '3+')
	tmp		<- rh[, which(SEXP1YR%in%c('0') & RLTN_NAMED>0)]
	if( nrow(rh[tmp,]) )
		warning("found SEXP1YR%in%c('0') & RLTN_NAMED>0  --> set to had sex last year, n=", length(tmp))
	set(rh, rh[, which(SEXP1YR%in%c('0') & RLTN_NAMED==1)], 'SEXP1YR', '1')
	set(rh, rh[, which(SEXP1YR%in%c('0') & RLTN_NAMED==2)], 'SEXP1YR', '2')
	set(rh, rh[, which(SEXP1YR%in%c('0') & RLTN_NAMED==3)], 'SEXP1YR', '3')
	set(rh, rh[, which(SEXP1YR%in%c('0') & RLTN_NAMED==4)], 'SEXP1YR', '3+')
	#	revisit SEXYEAR based on SEXP1YR
	set(rh, rh[, which(SEXYEAR%in%c('Unknown') & !SEXP1YR%in%c('0','Unknown'))], 'SEXYEAR', 'Y')
	tmp		<- rh[, which(SEXYEAR%in%c('N') & !SEXP1YR%in%c('0','Unknown'))]
	if( nrow(rh[tmp,]) )
		warning("found SEXYEAR%in%c('N') & !SEXP1YR%in%c('0','Unknown') --> set to had sex last year, n=", length(tmp))	
	set(rh, tmp, 'SEXYEAR', 'Y')		
	set(ra, ra[, which(SEXYEAR%in%c('Unknown') & !SEXP1YR%in%c('0','Unknown'))], 'SEXYEAR', 'Y')
	tmp		<- ra[, which(SEXYEAR%in%c('N') & !SEXP1YR%in%c('0','Unknown'))]
	if( nrow(ra[tmp,]) )
		warning("found SEXYEAR%in%c('N') & !SEXP1YR%in%c('0','Unknown') --> set to had sex last year, n=", length(tmp))	
	set(ra, tmp, 'SEXYEAR', 'Y')		
	#	define SEXP1OUT 
	setnames(rh, c('SEXP1OUT'), c('SEXP1OUT_'))	
	rh[, SEXP1OUT:= as.character(SEXP1OUT_)]
	set(rh, rh[, which(SEXP1OUT_==92)],'SEXP1OUT','<3')
	set(rh, rh[, which(SEXP1OUT_==93)],'SEXP1OUT','3+')
	set(rh, rh[, which(SEXP1OUT_%in%c(97,98,99) | is.na(SEXP1OUT_))],'SEXP1OUT','Unknown')
	set(rh, NULL, 'SEXP1OUT_', NULL)	
	setnames(ra, c('SEXP1OUT'), c('SEXP1OUT_'))	
	ra[, SEXP1OUT:= as.character(SEXP1OUT_)]
	set(ra, ra[, which(SEXP1OUT_==92)],'SEXP1OUT','<3')
	set(ra, ra[, which(SEXP1OUT_==93)],'SEXP1OUT','3+')
	set(ra, ra[, which(SEXP1OUT_%in%c(97,98,99) | is.na(SEXP1OUT_))],'SEXP1OUT','Unknown')
	set(ra, NULL, 'SEXP1OUT_', NULL)	
	#	revisit SEXYEAR based on SEXP1OUT
	set(rh, rh[, which(SEXYEAR%in%c('Unknown') & !SEXP1OUT%in%c('0','Unknown'))], 'SEXYEAR', 'Y')
	tmp		<- rh[, which(SEXYEAR%in%c('N') & !SEXP1OUT%in%c('0','Unknown'))]
	if( nrow(rh[tmp,]) )
		warning("found SEXYEAR%in%c('N') & !SEXP1OUT%in%c('0','Unknown') --> set to had sex last year, n=", length(tmp))
	set(rh, tmp, 'SEXYEAR', 'Y')
	set(ra, ra[, which(SEXYEAR%in%c('Unknown') & !SEXP1OUT%in%c('0','Unknown'))], 'SEXYEAR', 'Y')
	tmp		<- ra[, which(SEXYEAR%in%c('N') & !SEXP1OUT%in%c('0','Unknown'))]
	if( nrow(ra[tmp,]) )
		warning("found SEXYEAR%in%c('N') & !SEXP1OUT%in%c('0','Unknown') --> set to had sex last year, n=", length(tmp))
	set(ra, tmp, 'SEXYEAR', 'Y')	
	#	revisit SEXP1YR based on SEXP1OUT
	set(rh, rh[, which(SEXP1OUT=='3+' & SEXP1YR%in%c('0','1','2','<3','Unknown'))], 'SEXP1YR', '3+')
	set(rh, rh[, which(SEXP1OUT=='<3' & SEXP1YR%in%c('0','Unknown'))], 'SEXP1YR', '<3')
	rh[, DUMMY:= seq_len(nrow(rh))]
	warning("set(rh, rh[, which(RID=='G013746' & VISIT==14)], 'SEXP1OUT', 1) --> think this is typo")
	set(rh, rh[, which(RID=='G013746' & VISIT==14)], 'SEXP1OUT', '1')	
	tmp		<- rh[, which(	!SEXP1OUT%in%c('3+','<3','Unknown') & !SEXP1YR%in%c('3+','<3','Unknown'))]  
	tmp		<- subset(rh[tmp, ], as.numeric(SEXP1OUT)>as.numeric(SEXP1YR))[, DUMMY]
	warning("as.numeric(SEXP1OUT)>as.numeric(SEXP1YR), set to SEXP1OUT, n=", length(tmp))
	set(rh, tmp, 'SEXP1YR', rh[tmp, SEXP1OUT])	
	set(ra, ra[, which(SEXP1OUT=='3+' & SEXP1YR%in%c('0','1','2','<3','Unknown'))], 'SEXP1YR', '3+')
	set(ra, ra[, which(SEXP1OUT=='<3' & SEXP1YR%in%c('0','Unknown'))], 'SEXP1YR', '<3')
	ra[, DUMMY:= seq_len(nrow(ra))]
	tmp		<- ra[, which(	!SEXP1OUT%in%c('3+','<3','Unknown') & !SEXP1YR%in%c('3+','<3','Unknown'))]  
	tmp		<- subset(ra[tmp, ], as.numeric(SEXP1OUT)>as.numeric(SEXP1YR))[, DUMMY]
	warning("as.numeric(SEXP1OUT)>as.numeric(SEXP1YR), set to SEXP1OUT, n=", length(tmp))
	set(ra, tmp, 'SEXP1YR', ra[tmp, SEXP1OUT])	
	#	revisit SEXP1YR based on SEXYEAR
	set(rh, rh[, which(SEXYEAR=='Unknown' & SEXP1YR=='0')], 'SEXP1YR', 'Unknown')
	set(rh, rh[, which(SEXYEAR=='Y' & SEXP1YR=='0')], 'SEXP1YR', 'Unknown')
	set(ra, ra[, which(SEXYEAR=='Unknown' & SEXP1YR=='0')], 'SEXP1YR', 'Unknown')
	set(ra, ra[, which(SEXYEAR=='Y' & SEXP1YR=='0')], 'SEXP1YR', 'Unknown')	
	#	define SEXPEVER
	setnames(rh, c('SEXPEVER'), c('SEXPEVER_'))	
	rh[, SEXPEVER:= as.character(SEXPEVER_)]
	set(rh, rh[, which(SEXPEVER_==92)],'SEXPEVER','<3')
	set(rh, rh[, which(SEXPEVER_==93)],'SEXPEVER','3+')
	set(rh, rh[, which(SEXPEVER_%in%c(97,98,99) | is.na(SEXPEVER_))],'SEXPEVER','Unknown')
	set(rh, NULL, 'SEXPEVER_', NULL)
	stopifnot( !nrow(subset(rh, is.na(SEXPEVER))) )	
	setnames(ra, c('SEXPEVER'), c('SEXPEVER_'))	
	ra[, SEXPEVER:= as.character(SEXPEVER_)]
	set(ra, ra[, which(SEXPEVER_==92)],'SEXPEVER','<3')
	set(ra, ra[, which(SEXPEVER_==93)],'SEXPEVER','3+')
	set(ra, ra[, which(SEXPEVER_%in%c(97,98,99) | is.na(SEXPEVER_))],'SEXPEVER','Unknown')
	set(ra, NULL, 'SEXPEVER_', NULL)
	stopifnot( !nrow(subset(ra, is.na(SEXPEVER))) )	
	#	revisit SEXPEVER based on SEXP1YR
	set(rh, rh[, which(SEXP1YR=='3+' & SEXPEVER%in%c('0','1','2','<3','Unknown'))], 'SEXPEVER', '3+')
	set(rh, rh[, which(SEXP1YR=='<3' & SEXPEVER%in%c('0','Unknown'))], 'SEXPEVER', '<3')	
	tmp		<- rh[, which(	!SEXP1YR%in%c('3+','<3','Unknown') & !SEXPEVER%in%c('3+','<3','Unknown'))]  
	tmp		<- subset(rh[tmp, ], as.numeric(SEXP1YR)>as.numeric(SEXPEVER))[, DUMMY]
	warning("as.numeric(SEXP1YR)>as.numeric(SEXPEVER), set to SEXP1YR, n=", length(tmp))
	set(rh, tmp, 'SEXPEVER', rh[tmp, SEXP1YR])	
	set(ra, ra[, which(SEXP1YR=='3+' & SEXPEVER%in%c('0','1','2','<3','Unknown'))], 'SEXPEVER', '3+')
	set(ra, ra[, which(SEXP1YR=='<3' & SEXPEVER%in%c('0','Unknown'))], 'SEXPEVER', '<3')	
	tmp		<- ra[, which(	!SEXP1YR%in%c('3+','<3','Unknown') & !SEXPEVER%in%c('3+','<3','Unknown'))]  
	tmp		<- subset(ra[tmp, ], as.numeric(SEXP1YR)>as.numeric(SEXPEVER))[, DUMMY]
	warning("as.numeric(SEXP1YR)>as.numeric(SEXPEVER), set to SEXP1YR, n=", length(tmp))
	set(ra, tmp, 'SEXPEVER', ra[tmp, SEXP1YR])	
	#	define EVERSEX 
	set(rh, rh[, which(is.na(EVERSEX))],'EVERSEX', 99)
	set(rh, NULL, 'EVERSEX', rh[, gsub('_[0-9]$','',as.character(factor(EVERSEX, levels=c(0,1,2,3,8,99), labels=c('Unknown_1','Y','N','Unknown_2','Unknown_3','Unknown_4'))))])
	stopifnot( !nrow(subset(rh, is.na(EVERSEX))) )	
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & SEXYEAR=='Y')], 'EVERSEX', 'Y')
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & CURRMARR==1)], 'EVERSEX', 'Y')
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & EVERMARR==1)], 'EVERSEX', 'Y')
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & RLTN_NAMED>0)], 'EVERSEX', 'Y')	
	set(ra, ra[, which(is.na(EVERSEX))],'EVERSEX', 99)
	set(ra, NULL, 'EVERSEX', ra[, gsub('_[0-9]$','',as.character(factor(EVERSEX, levels=c(0,1,2,3,8,99), labels=c('Unknown_1','Y','N','Unknown_2','Unknown_3','Unknown_4'))))])
	stopifnot( !nrow(subset(ra, is.na(EVERSEX))) )	
	set(ra, ra[, which(EVERSEX%in%c('Unknown') & SEXYEAR=='Y')], 'EVERSEX', 'Y')
	set(ra, ra[, which(EVERSEX%in%c('Unknown') & CURRMARR==1)], 'EVERSEX', 'Y')
	set(ra, ra[, which(EVERSEX%in%c('Unknown') & EVERMARR==1)], 'EVERSEX', 'Y')		
	#	revisit EVERSEX based on SEXPEVER
	set(rh, rh[, which(EVERSEX%in%c('Unknown') & !SEXPEVER%in%c('0','Unknown'))], 'EVERSEX', 'Y')
	tmp		<- rh[, which(EVERSEX%in%c('N') & !SEXPEVER%in%c('0','Unknown'))]
	if( nrow(rh[tmp,]) )
		warning("found EVERSEX%in%c('N') & !SEXPEVER%in%c('0','Unknown') --> set to had sex ever, n=", length(tmp))
	set(rh, tmp, 'EVERSEX', 'Y')
	set(ra, ra[, which(EVERSEX%in%c('Unknown') & !SEXPEVER%in%c('0','Unknown'))], 'EVERSEX', 'Y')
	tmp		<- ra[, which(EVERSEX%in%c('N') & !SEXPEVER%in%c('0','Unknown'))]
	if( nrow(ra[tmp,]) )
		warning("found EVERSEX%in%c('N') & !SEXPEVER%in%c('0','Unknown') --> set to had sex ever, n=", length(tmp))
	set(ra, tmp, 'EVERSEX', 'Y')	
	#	revisit SEXPEVER based on EVERSEX
	set(rh, rh[, which(EVERSEX=='Unknown' & SEXPEVER=='0')], 'SEXP1YR', 'Unknown')
	set(rh, rh[, which(EVERSEX=='Y' & SEXPEVER=='0')], 'SEXP1YR', 'Unknown')
	set(ra, ra[, which(EVERSEX=='Unknown' & SEXPEVER=='0')], 'SEXP1YR', 'Unknown')
	set(ra, ra[, which(EVERSEX=='Y' & SEXPEVER=='0')], 'SEXP1YR', 'Unknown')		
	#	check SEXPEVER	
	tmp		<- subset(rh, SEXPEVER=='0' & SEXACTIVE==1)
	if( nrow(tmp) )
		warning("found SEXPEVER=='0' & SEXACTIVE==1 --> report only, n=", nrow(tmp))	
	#	add extra-marital partner to MARSTAT
	tmp	<- rh[, which(MULTIPART>0 & MARSTAT!='Previously Married')]
	set(rh, tmp, 'MARSTAT', rh[tmp, paste0(MARSTAT,' + casual partner')])
	tmp	<- rh[, which(MARSTAT%in%c('Monogamous') & !SEXP1YR%in%c('1','Unknown'))]
	set(rh, tmp, 'MARSTAT', rh[tmp, paste0(MARSTAT,' + casual partner')])	
	tmp	<- rh[, which(MARSTAT%in%c('Polygamous') & !SEXP1YR%in%c('1','Unknown') & POLYMAR<as.numeric(gsub('Unknown','0',gsub('+','',SEXP1YR,fixed=1))) ) ]
	set(rh, tmp, 'MARSTAT', rh[tmp, paste0(MARSTAT,' + casual partner')])
	tmp	<- ra[, which(MARSTAT%in%c('Monogamous') & !SEXP1YR%in%c('1','Unknown'))]
	set(ra, tmp, 'MARSTAT', ra[tmp, paste0(MARSTAT,' + casual partner')])	
	tmp	<- ra[, which(MARSTAT%in%c('Polygamous') & !SEXP1YR%in%c('1','Unknown') & POLYMAR<as.numeric(gsub('Unknown','0',gsub('+','',SEXP1YR,fixed=1))) ) ]
	set(ra, tmp, 'MARSTAT', ra[tmp, paste0(MARSTAT,' + casual partner')])
	#	define ever alcohol use during sex
	set(rh, NULL, 'ALC', rh[, as.character(factor(ALC, levels=c(0,1), labels=c('N','Y')))])
	#	redefine SEXP1YR
	set(rh, rh[, which(!SEXP1YR%in%c('0','1','2','Unknown'))], 'SEXP1YR','3+')
	set(ra, ra[, which(!SEXP1YR%in%c('0','1','2','Unknown'))], 'SEXP1YR','3+')
	#	redefine SEXP1OUT
	set(rh, rh[, which(!SEXP1OUT%in%c('0','1','2','Unknown'))], 'SEXP1OUT','3+')
	set(ra, ra[, which(!SEXP1OUT%in%c('0','1','2','Unknown'))], 'SEXP1OUT','3+')
	#	check circumcision
	tmp		<- rh[, which(CIRCUM=='Y' & SEX=='F')]
	cat('\nWarning: found female circumcised --> set to NA' ,rh[tmp, paste(RID, collapse=' ')])	
	set(rh, tmp, 'CIRCUM', NA_integer_)
	#	define COMM_TYPE
	set(rh, NULL, 'COMM_NUM', rh[, as.character(COMM_NUM)])
	set(ra, NULL, 'COMM_NUM', ra[, as.character(COMM_NUM)])
	tmp		<- data.table(	COMM_NUM=	c("1","2","3","4","5","6","7","8","9","14","15","16","18","19","22","23","24","25","29","32","33","34","35","36","38","40","44","45","46","51","52","53","54","55", "56","57","58","59","60","61","62","65","67","74","77","81","84","89","94","95","103","106","107","108","109","120","177", "183", "256", "370","391","401","451", "468","602", "754", "755", "760", "770","771","772","773","774","776"),
			COMM_TYPE=	c("T","A","A","T","A","A","A","A","A", "A", "A", "T", "A", "A", "T", "A", "T", "A", "A", "A", "T", "A", "A", "A", "F", "A", "A", "A", "A", "T", "A", "A", "A", "A",  "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",  "A","A",   "A",  "T",  "A",  "A",  "A",  "A",   "A",   "A",   "A",  "A", "A",  "A",    "A",  "A",  "A",    "A",   "A", "F",  "F",  "A",  "A",  "F",   "T"))
	set(tmp, NULL, 'COMM_TYPE', tmp[, as.character(factor(COMM_TYPE, levels=c('A','T','F'), labels=c('agrarian','trading','fisherfolk')))])		
	stopifnot(!length(setdiff( rh[, COMM_NUM], tmp[, COMM_NUM] )))
	stopifnot(!length(setdiff( ra[, COMM_NUM], tmp[, COMM_NUM] )))
	rh		<- merge(rh, tmp, by='COMM_NUM')
	ra		<- merge(ra, tmp, by='COMM_NUM')	
	#	merge community numbers for same community
	set(rh, NULL, 'COMM_NUM', rh[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',COMM_NUM))))])
	set(ra, NULL, 'COMM_NUM', ra[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',COMM_NUM))))])
	set(rd, NULL, 'COMM_NUM', rd[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',COMM_NUM))))])
	#	add anonymized IDs
	dc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv'))
	dc		<- unique(subset(dc, select=c('COMM_NUM','COMM_NUM_A')))	
	rh		<- merge(rh, dc, by='COMM_NUM')
	ra		<- merge(ra, dc, by='COMM_NUM')
	rd		<- merge(rd, dc, by='COMM_NUM')
	#	set to NULL
	set(rh, NULL, c('VDEX','EVERMARR','CURRMARR','RELIGION','POLYMAR','DUMMY','RLTN1','RLTN2','RLTN3','RLTN4','RLTN_NAMED'), NULL)
	set(rh, NULL, c('BVDEX','EVERSEX','SEXGIFT','SEXYEAR','EDUCATE','EDUCYRS','EDCAT','ARVMED','CNDEVER1','RNYRCON1','CNDEVER2','RNYRCON2','CNDEVER3','RNYRCON3','CNDEVER4','RNYRCON4','RLTNCON1'),NULL)
	set(rh, NULL, c('RLTNCON2','RLTNCON3','RLTNCON4','ALC1B','ALC2B','ALC3B','ALC4B','ALC1F','ALC2F','ALC3F','ALC4F','OCCUP11','OCCUP12','OCCUP13','OCCUP14','OCCUP21','OCCUP22','OCCUP23','OCCUP24','SEXHIGH','SEXOUT'),NULL)
	set(rh, NULL, c('SEXCAT','PREVMAR','AGECAT','AGECAT2','HIVPREV2','UNDER25','AGE15TO19','AGE20TO24','AGE25TO29','AGE30TO34','AGE35TO39','AGE40TO44','AGE45TO49','OCCLAG1','SUM_ALC'),NULL)
	set(rh, NULL, c('SEXACTIVE','MULTIPART','CAS','SUMCON','CONCON','NEVERSEX','OCCUP1','OCCUP2','SEXPEVER'),NULL)
	set(rn, NULL, c('SAMPLEREASON'), NULL)
	rd[, COHORT:= 'RCCS']
	rh[, COHORT:= 'RCCS']
	ra[, COHORT:= 'RCCS']
	#
	list(rd=rd, rh=rh, ra=ra, rn=rn)
}

RakaiFull.phylogeography.180322.get.data.eligibility.participation.sequenced<- function()
{
	require(data.table)
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/rakai_elibility.rda"
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	load(infile)
	
	#	subset to data of interest
	de	<- as.data.table(eldat)	
	de	<- subset(de, status%in%c('_Participated','Away','Blood refusal','Missing data','Other','Refused','urine sample'))
	de	<- subset(de, visit<17)
	de	<- subset(de, select=c(PERM_ID, CURR_ID, visit, RESIDENT, MOBILITY, REGION, COMM_NUM, HH_NUM, SEX, STUDY_ID, status, date, AGEYRS, inmigrant))
	setnames(de, colnames(de), toupper(colnames(de)))
	#	define PARTICIPATED as "participated, missing data ok"
	#	TODO take out missing data
	de[, PARTICIPATED:= as.integer(STATUS%in%c('_Participated'))]
	de[, DATE:= hivc.db.Date2numeric(DATE)]	
	
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
	#	there are a few individuals with Study ID who did not "participate" when missing not included
	#subset(de, PARTICIPATED_ANY_VISIT==0 & !is.na(STUDY_ID))
	
	#	prepare HIV status
	tmp		<- RakaiCirc.epi.get.info.170208()
	rd		<- tmp$rd
	rneuro	<- tmp$rn
	ra		<- tmp$ra
	#	prepare self report ART
	#	find those who report to be on ART at their first visit in 15-17
	rart	<- merge(ra, ra[VISIT>=14 & VISIT<17, list(VISIT=min(VISIT)), by='RID'], by=c('RID','VISIT'))	
	set(rart, rart[, which(is.na(ARVMED))], 'ARVMED', 2L)
	setnames(rart, 'ARVMED', 'SELFREPORTART_AT_FIRST_VISIT')
	rart	<- subset(rart, select=c(RID, SELFREPORTART_AT_FIRST_VISIT))
	#set(rart, NULL,  'SELFREPORTART', rart[, factor(SELFREPORTART, levels=c(0,1,2), labels=c('no','yes','unknown'))])
	#	select meta-data closest to first pos date	
	tmp		<- rd[, list(VISIT=VISIT[which.min(abs(DATE-FIRSTPOSDATE))]), by='RID']
	tmp2	<- rd[, list(PANGEA=as.integer(any(PANGEA==1))), by='RID']
	rd		<- merge(rd, tmp, by=c('RID','VISIT'))
	rd[, HIV:= 1L]
	rd[, PANGEA:=NULL]
	rd		<- merge(rd, tmp2, by='RID')
	#
	ra		<- unique(subset(ra, !is.na(FIRSTPOSDATE), select=c(RID, SEX, VISIT, VISDATE, FIRSTPOSDATE, ARVMED, COMM_NUM, COMM_NUM_A)))
	tmp		<- ra[, list(VISIT=VISIT[which.min(abs(VISDATE-FIRSTPOSDATE))]), by='RID']	
	ra		<- merge(ra, tmp, by=c('RID','VISIT'))
	ra[, HIV_1517:= 1L]
	#rd		<- subset(rd, BIRTHYR>2010-50 & (is.na(EST_DATEDIED) | (!is.na(EST_DATEDIED) & EST_DATEDIED>2010)))
		
	#	add HIV status
	tmp		<- rd[FIRSTPOSDATE<2015+2/12, list(HIV=max(HIV), PANGEA=max(PANGEA)), by='RID'] 			
	de		<- merge(de, tmp, by='RID', all.x=1)	
	tmp		<- unique(subset(ra, FIRSTPOSDATE<2015+2/12, select=c(RID, HIV_1517, FIRSTPOSDATE)))
	de		<- merge(de, tmp, by='RID', all.x=1)
	
	#	add ART status
	de		<- merge(de, rart, by='RID', all.x=1)
	
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
		
	#	tmp above has 4074 individuals
	#	of those, 3758 are in 'de', ie participated in rounds 15-17
	
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
	tmp		<- de[, which(PARTICIPATED_ANY_VISIT==1 & is.na(HIV_1517))]
	set(de, tmp, 'HIV_1517', 0L)
	tmp		<- de[, which(PARTICIPATED_ANY_VISIT==1 & is.na(BAM_OUTPUT))]
	set(de, tmp, c('BAM_OUTPUT','MIN_PNG_OUTPUT'), 0L)
	tmp		<- rd[, which(is.na(BAM_OUTPUT))]
	set(rd, tmp, c('BAM_OUTPUT','MIN_PNG_OUTPUT'), 0L)
	
	de[, table(HIV, HIV_1517)]	
	# 672 with HIV_1517==1 and HIV==0
	#	0 with HIV_1517==0 and HIV==1
	# use HIV_1517 below
		
	rneuro[, sum(MIN_PNG_OUTPUT)]			# 224	
	rd[, sum(MIN_PNG_OUTPUT)]				# 2746
	de[, sum(MIN_PNG_OUTPUT, na.rm=TRUE)]	# 2689
	
	#	calculate number for whom we have deep seq output
	tmp		<- subset(rd, MIN_PNG_OUTPUT==1)[, list(DEEP_SEQ=length(RID)), by=c('COMM_NUM')]	
	tmp		<- subset(tmp, DEEP_SEQ>0)
	tmp[, DEEP_SEQ:=NULL]	
	tmp[, COMM_ANY_MIN_PNG_OUTPUT:= 1]
	de		<- merge(de, tmp, by='COMM_NUM', all.x=1)
	rd		<- merge(rd, tmp, by='COMM_NUM', all.x=1)
	set(de, de[, which(is.na(COMM_ANY_MIN_PNG_OUTPUT))], 'COMM_ANY_MIN_PNG_OUTPUT', 0L)
	set(rd, rd[, which(is.na(COMM_ANY_MIN_PNG_OUTPUT))], 'COMM_ANY_MIN_PNG_OUTPUT', 0L)
	tmp		<- unique(subset(rd, select=c(COMM_NUM,COMM_NUM_A)))
	de		<- merge(de, tmp, by='COMM_NUM')
	
	#	subset to only those communities with seq data
	#	4 fishing, 34 inland
	#def		<- copy(de)
	#de		<- subset(de, !is.na(DUMMY))
	#de[, DUMMY:=NULL]
	#rd		<- subset(rd, !is.na(DUMMY))
	#rd[, DUMMY:=NULL]	
	#

	#	add age category at midpoint of observation period, 2013.25
	set(de, de[, which(is.na(DATE))], 'DATE', hivc.db.Date2numeric(as.Date("2011-09-02")))	
	de[, AGE_AT_MID:= 2013.25 - (hivc.db.Date2numeric(DATE)-AGEYRS)]
	de[, AGE_AT_MID_C:= cut(AGE_AT_MID, breaks=c(10,25,35,65), labels=c('15-24','25-34','35+'), right=FALSE)]
	stopifnot( nrow(subset(de, is.na(AGE_AT_MID_C)))==0 )
	
	
	
	#	prepare inmigrant -- identify inmigrants from fishing communities and from external
	infile		<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData_v2.rda"
	load(infile)
	inmigrant	<- as.data.table(inmigrant)	
	#	plot fisherfolk	to figure out how much of a radius we need
	if(0)
	{
		zf		<- data.table(longitude=c(31.763,31.7968,31.754,31.838), latitude=c(-0.915, -0.6518, -0.703, -0.497), ID= c('Kasensero','Bukyanju','NearBwende','Fish4'))
		make_circles <- function(centers, radius, nPoints = 100){
			# centers: the data frame of centers with ID
			# radius: radius measured in kilometer
			#
			meanLat <- mean(centers$latitude)
			# length per longitude changes with lattitude, so need correction
			radiusLon <- radius /111 / cos(meanLat/57.3) 
			radiusLat <- radius / 111
			circleDF <- data.frame(ID = rep(centers$ID, each = nPoints))
			angle <- seq(0,2*pi,length.out = nPoints)
			
			circleDF$lon <- unlist(lapply(centers$longitude, function(x) x + radius * cos(angle)))
			circleDF$lat <- unlist(lapply(centers$latitude, function(x) x + radius * sin(angle)))
			return(circleDF)
		}
		zc <- make_circles(zf, 0.01)
		ggmap(zm) +			
				geom_point(data=zf, aes(x=longitude, y=latitude, pch=ID), stroke=1.5, alpha=0.8) +
				geom_polygon(data=zc, aes(lon, lat, group = ID), color = "red", alpha = 0)
		#	radius of length 0.01 should catch					
		tmp		<- inmigrant[, list( 	DIST_KASENSERO= sqrt( (inmig_lon- 31.763)^2 + (inmig_lat - (-0.915))^2),
						DIST_BUKYANJU= sqrt( (inmig_lon- 31.7968)^2 + (inmig_lat - (-0.6518))^2),
						DIST_NEARBWENDE= sqrt( (inmig_lon- 31.754)^2 + (inmig_lat - (-0.703))^2),
						DIST_FISH4= sqrt( (inmig_lon- 31.838)^2 + (inmig_lat - (-0.497))^2)
				), by=c('RCCS_studyid','visit')]
		tmp		<- melt(tmp, id.vars=c('RCCS_studyid','visit'))
		ggplot(subset(tmp, value<0.3), aes(x=value)) +
				geom_histogram(binwidth=0.01) +
				facet_grid(variable~.)
		#	1 looks good
		zfd		<- merge(inmigrant, subset(tmp, value<0.01, c(RCCS_studyid, visit)), by=c('RCCS_studyid','visit'))
		#	so fishing sites are MALEMBO DIMU KASENSERO NAMIREMBE 
		#	but there are spelling mistakes		
	}
	#
	#	clean up inmigrant	
	#
	#	inmigrant[, unique(sort(inmig_place))]
	set(inmigrant, NULL, 'inmig_place', inmigrant[, gsub('DDIMO|DDIMU|DIMO|DIMU','DIMU',inmig_place)])
	set(inmigrant, inmigrant[, which(grepl('MALEMBO',inmig_place))], 'inmig_place', 'MALEMBO')
	set(inmigrant, NULL, 'inmig_place', inmigrant[, gsub("KASEMSERO","KASENSERO",inmig_place)])
	set(inmigrant, inmigrant[, which(grepl('KASENSERO',inmig_place))], 'inmig_place', 'KASENSERO')
	
	#	define from_fishing and from_outside and from_inland
	inmigrant[, INMIG_LOC:= 'inland' ]
	set(inmigrant, inmigrant[, which(grepl('MALEMBO|DIMU|KASENSERO|NAMIREMBE',inmig_place))], 'INMIG_LOC','fisherfolk')
	set(inmigrant, inmigrant[, which(inmig_admin0!='Uganda')], 'INMIG_LOC','external')
	set(inmigrant, inmigrant[, which(inmig_admin1!='Rakai')], 'INMIG_LOC','external')
	set(inmigrant, inmigrant[, which(is.na(inmig_admin1))], 'INMIG_LOC','origin_unknown')
	inmigrant[, table(INMIG_LOC)]
	#	  external fisherfolk     inland    unknown 
    #   	762         40       1571        381 
	inmigrant	<- subset(inmigrant, select=c(RCCS_studyid, date, INMIG_LOC))	
	inmigrant[, date:= hivc.db.Date2numeric(date)]	
	setnames(inmigrant, colnames(inmigrant), gsub('DATE','INMIGRATE_DATE',gsub('RCCS_STUDYID','RID',toupper(colnames(inmigrant)))))
	inmigrant	<- merge(inmigrant, de, by='RID', all.x=TRUE)
	#	ignore inmigrants not seen in 15-16
	inmigrant	<- subset(inmigrant, !is.na(STATUS))
	#	only inmigration date closest to first visit and before the first visit
	tmp			<- inmigrant[, list(FIRSTVISITDATE=min(DATE)), by='RID']
	inmigrant	<- merge(inmigrant, tmp, by='RID')	
	tmp			<- inmigrant[INMIGRATE_DATE<=FIRSTVISITDATE, list(INMIGRATE_DATE= INMIGRATE_DATE[which.min(FIRSTVISITDATE-INMIGRATE_DATE)]), by='RID']
	inmigrant	<- merge(inmigrant, tmp, by=c('RID','INMIGRATE_DATE'))
	inmigrant	<- unique(inmigrant, by=c('RID','VISIT'))
	#	select those inmigrated in the last 2 years
	tmp			<- subset(inmigrant, (FIRSTVISITDATE-INMIGRATE_DATE)<=2, c(RID, INMIG_LOC))
	setnames(tmp, 'INMIG_LOC', 'INMIG_2YRS_LOC')
	tmp[, INMIG_2YRS:=1]
	de			<- merge(de, tmp, by='RID', all.x=TRUE)
	set(de, de[, which(is.na(INMIG_2YRS))], 'INMIG_2YRS', 0)
	set(de, de[, which(is.na(INMIG_2YRS_LOC))], 'INMIG_2YRS_LOC', 'resident')	
	
	#	some fixup
	set(de, de[, which(PARTICIPATED_ANY_VISIT==0 & MIN_PNG_OUTPUT==1)], 'PARTICIPATED_ANY_VISIT', 1L)
	set(de, de[, which(SELFREPORTART_AT_FIRST_VISIT==1 & MIN_PNG_OUTPUT==1)], 'SELFREPORTART_AT_FIRST_VISIT', 0L)
	
	
	de[, table(HIV_1517, SELFREPORTART_AT_FIRST_VISIT, useNA='if')]
	#	        SELFREPORTART
	#	HIV_1517     0         1     2  <NA>
    #		    0    20715     8     7    10
    #			1     3884  1256     2     0
    #			<NA>     0     0     0 11763
		
	#	now calculate participation by community and gender
	des		<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','COMM_ANY_MIN_PNG_OUTPUT','SEX','PARTICIPATED_ANY_VISIT')]
	set(des, NULL, 'PARTICIPATED_ANY_VISIT', des[, factor(PARTICIPATED_ANY_VISIT, levels=c('1','0'), labels=c('PART_EVER','PART_NEVER'))])
	des		<- dcast.data.table(des, COMM_NUM+COMM_NUM_A+COMM_ANY_MIN_PNG_OUTPUT+SEX~PARTICIPATED_ANY_VISIT, value.var='N')
	
	#	now calculate HIV pos by end of round 16, and sequenced by end of round 16, by community and gender
	#	among those that participated in rounds 15 - 16
	tmp		<- subset(de, PARTICIPATED_ANY_VISIT==1)
	tmp		<- tmp[, list(	HIV_1516_YES=length(which(HIV_1517==1)), 
							HIV_1516_NO=length(which(HIV_1517==0)),
							SLART_AT_FIRST_VISIT=length(which(SELFREPORTART_AT_FIRST_VISIT==1)),
							DEEP_SEQ_1516=length(which(MIN_PNG_OUTPUT==1)) ), by=c('COMM_NUM','SEX')]	
	des		<- merge(des, tmp, by=c('COMM_NUM','SEX'), all=TRUE)
	
	#	add community type
	des[, COMM_TYPE:= as.character(factor(substr(COMM_NUM_A,1,1)=='f',levels=c(TRUE,FALSE),labels=c('fisherfolk','inland')))]	
	
	#	now do the same by community and gender and age category
	desa	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','COMM_ANY_MIN_PNG_OUTPUT','SEX','AGE_AT_MID_C','PARTICIPATED_ANY_VISIT')]
	set(desa, NULL, 'PARTICIPATED_ANY_VISIT', desa[, factor(PARTICIPATED_ANY_VISIT, levels=c('1','0'), labels=c('PART_EVER','PART_NEVER'))])
	desa	<- dcast.data.table(desa, COMM_NUM+COMM_NUM_A+COMM_ANY_MIN_PNG_OUTPUT+SEX+AGE_AT_MID_C~PARTICIPATED_ANY_VISIT, value.var='N')
	tmp		<- subset(de, PARTICIPATED_ANY_VISIT==1)
	tmp		<- tmp[, list(	HIV_1516_YES=length(which(HIV_1517==1)), 
							HIV_1516_NO=length(which(HIV_1517==0)), 
							SLART_AT_FIRST_VISIT=length(which(SELFREPORTART_AT_FIRST_VISIT==1)),
							DEEP_SEQ_1516=length(which(MIN_PNG_OUTPUT==1)) ), by=c('COMM_NUM','SEX','AGE_AT_MID_C')]	
	desa	<- merge(desa, tmp, by=c('COMM_NUM','SEX','AGE_AT_MID_C'), all=TRUE)
	desa[, COMM_TYPE:= as.character(factor(substr(COMM_NUM_A,1,1)=='f',levels=c(TRUE,FALSE),labels=c('fisherfolk','inland')))]
	#	check for missing and replace by 0		
	for(x in c('PART_EVER','PART_NEVER','HIV_1516_YES','HIV_1516_NO','SLART_AT_FIRST_VISIT','DEEP_SEQ_1516'))
		set(desa, which(is.na(desa[[x]])), x, 0)
	
	#	now do the same by community and gender and migration category
	desm	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','COMM_ANY_MIN_PNG_OUTPUT','SEX','INMIG_2YRS_LOC','PARTICIPATED_ANY_VISIT')]
	set(desm, NULL, 'PARTICIPATED_ANY_VISIT', desm[, factor(PARTICIPATED_ANY_VISIT, levels=c('1','0'), labels=c('PART_EVER','PART_NEVER'))])
	desm	<- dcast.data.table(desm, COMM_NUM+COMM_NUM_A+COMM_ANY_MIN_PNG_OUTPUT+SEX+INMIG_2YRS_LOC~PARTICIPATED_ANY_VISIT, value.var='N')
	tmp		<- subset(de, PARTICIPATED_ANY_VISIT==1)
	tmp		<- tmp[, list(	HIV_1516_YES=length(which(HIV_1517==1)), 
							HIV_1516_NO=length(which(HIV_1517==0)),
							SLART_AT_FIRST_VISIT=length(which(SELFREPORTART_AT_FIRST_VISIT==1)),
							DEEP_SEQ_1516=length(which(MIN_PNG_OUTPUT==1)) ), by=c('COMM_NUM','SEX','INMIG_2YRS_LOC')]	
	desm	<- merge(desm, tmp, by=c('COMM_NUM','SEX','INMIG_2YRS_LOC'), all=TRUE)
	desm[, COMM_TYPE:= as.character(factor(substr(COMM_NUM_A,1,1)=='f',levels=c(TRUE,FALSE),labels=c('fisherfolk','inland')))]
	#	check for missing and replace by 0		
	for(x in c('PART_EVER','PART_NEVER','HIV_1516_YES','HIV_1516_NO','SLART_AT_FIRST_VISIT','DEEP_SEQ_1516'))
		set(desm, which(is.na(desm[[x]])), x, 0)
	
	
	
	#	now calculate HIV pos by end of round 16, and sequenced by end of round 16, by community and gender
	#	among those that ever participated
	#	excluding those that died before 2010, or that reached age 60 before 2010
	tmp		<- subset(rd, is.na(EST_DATEDIED) | (!is.na(EST_DATEDIED) & EST_DATEDIED>=2010))	
	tmp		<- subset(tmp, (2010-BIRTHYR)<60)
	rds		<- tmp[, list(HIV_EVER_YES=length(which(HIV==1)), DEEP_SEQ_EVER=length(which(MIN_PNG_OUTPUT==1)) ), by=c('COMM_NUM','SEX')]	
	rds[, sum(DEEP_SEQ_EVER)]	
	# 2665 (among those in rd who are not too old and did not die before 2010)
	# now up to 2700

	save(des, desa, desm, rds, de, df, rd, file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender_age.rda')
}

RakaiFull.phylogeography.180322.table1<- function()
{
	require(data.table)
	require(Hmisc)	
	
	infile.allpairs				<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_networksallpairs.rda'	
	infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender_age.rda'
	outfile.base	<- gsub('_networksallpairs.rda','',infile.allpairs)
	load(infile)
	
	#	there s just 1 person ever sequenced in comms 401 and 55, ignore
	#	this means we are just working with des, not rds
	des		<- merge(des, rds, by=c('COMM_NUM','SEX'), all.x=TRUE)
	subset(des[, list(DEEP_SEQ_1516=sum(DEEP_SEQ_1516)), by='COMM_NUM'], DEEP_SEQ_1516==0)
	#	36, 55, 401
	#des		<- subset(des, !COMM_NUM%in%c('36','55','401'))
	#desa	<- subset(desa, !COMM_NUM%in%c('36','55','401'))
	
	length( setdiff( subset(tmp, !is.na(RID) & HIV==1 & MIN_PNG_OUTPUT==1)[, RID], subset(de, !is.na(RID) & HIV==1 & MIN_PNG_OUTPUT==1)[, RID] ) )
	#	57 individuals that have a sequence who are not in de but in rd
			
	#	overall sequence coverage 
	des[, sum(DEEP_SEQ_1516)]	#2652
	des[, sum(SLART_AT_FIRST_VISIT)]	#1264
	des[, sum(HIV_1516_YES)]	#5142
	des[, sum(HIV_EVER_YES)]	#5620
	des[, sum(DEEP_SEQ_1516)/sum(HIV_EVER_YES)]							# 	0.4761135
	des[, sum(DEEP_SEQ_1516)/sum(HIV_1516_YES)]							# 	0.5156584
	des[, sum(DEEP_SEQ_1516)/sum(HIV_1516_YES-SLART_AT_FIRST_VISIT)]	#	0.6838577
	des[, sum(DEEP_SEQ_1516)/( sum(HIV_EVER_YES * (PART_EVER+PART_NEVER)/PART_EVER) )]	# 0.3295489
	des[, sum(DEEP_SEQ_1516)/( sum(HIV_1516_YES * (PART_EVER+PART_NEVER)/PART_EVER) )]	# 0.3618695
	des[, sum(DEEP_SEQ_1516)/( sum((HIV_1516_YES-SLART_AT_FIRST_VISIT) * (PART_EVER+PART_NEVER)/PART_EVER) )]	# 0.481056
	
	#	make supplementary plots of participation and sequencing by gender
	des[, P_PART_RAW:= PART_EVER/(PART_EVER+PART_NEVER)]
	des[, P_SEQ_RAW_1516:= DEEP_SEQ_1516/(HIV_1516_YES-SLART_AT_FIRST_VISIT)]
		
	tmp		<- melt(des, id.vars=c('COMM_TYPE','COMM_NUM_A','SEX'), measure.vars=c('P_PART_RAW','P_SEQ_RAW_1516'))	
	set(tmp, NULL, 'COMM_TYPE', tmp[, factor(COMM_TYPE, levels=c('fisherfolk','inland'), labels=c('fishing site','inland community'))])
	set(tmp, NULL, 'SEX', tmp[, factor(SEX, levels=c('M','F'), labels=c('male','female'))])
	ggplot(tmp, aes(x=COMM_NUM_A, fill=SEX, y=value)) +
			geom_bar(position='dodge', stat='identity') +
			scale_fill_manual(values=c('female'='hotpink2', 'male'='steelblue2')) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			facet_grid(variable~COMM_TYPE, space='free', scales='free') +
			labs(x='\ncommunity', fill='gender') +
			theme_bw()
	ggsave(file=paste0(outfile.base,'sampling_differences_180322.pdf'), w=12, h=6)
	#
	
	#	number communities with some sequences
	des[, length(unique(COMM_NUM))]
	#	36


	#	make table 1 for paper
	#	get all transmission events part of transmission chains regardless of gender combinations		
	confidence.cut		<- 0.6
	load(infile.allpairs)		
	rtnn[, SXO:= paste0(ID1_SEX,ID2_SEX)]
	rtnn[, SELECT:= NA_character_]
	set(rtnn, rtnn[, which(is.na(PTY_RUN))], 'SELECT', 'insufficient deep sequence data for at least one partner of couple')
	set(rtnn, rtnn[, which(!is.na(PTY_RUN) & is.na(LINK_12) & is.na(LINK_21))], 'SELECT', 'couple most likely not a pair')
	set(rtnn, rtnn[, which(!is.na(POSTERIOR_SCORE_LINKED) & POSTERIOR_SCORE_LINKED<=confidence.cut)], 'SELECT', 'couple ambiguous if pair or not pair')
	set(rtnn, rtnn[, which(!is.na(POSTERIOR_SCORE_LINKED) & POSTERIOR_SCORE_LINKED>confidence.cut)], 'SELECT', 'couple most likely a pair direction not resolved')
	set(rtnn, rtnn[, which(!is.na(POSTERIOR_SCORE_LINKED) & POSTERIOR_SCORE_LINKED>confidence.cut & POSTERIOR_SCORE_12>confidence.cut)], 'SELECT', 'couple most likely a pair direction resolved to 12')
	set(rtnn, rtnn[, which(!is.na(POSTERIOR_SCORE_LINKED) & POSTERIOR_SCORE_LINKED>confidence.cut & POSTERIOR_SCORE_21>confidence.cut)], 'SELECT', 'couple most likely a pair direction resolved to 21')	
	rtnn	<- subset(rtnn, select=c(ID1, ID2, ID1_SEX, ID2_SEX, LINK_12, LINK_21, SXO, SELECT))
	#	determine first concordant pos visit; add metadata at first concordant pos visit
	tmp		<- unique(subset(rd, select=c(RID, FIRSTPOSVIS, FIRSTPOSDATE)))
	setnames(tmp, colnames(tmp), gsub('ID1_RID','ID1',paste0('ID1_',colnames(tmp))))
	rtnn	<- merge(rtnn, tmp, by='ID1')	
	setnames(tmp, colnames(tmp), gsub('ID1','ID2',colnames(tmp)))
	rtnn	<- merge(rtnn, tmp, by='ID2')			
	rtnn[, VISIT_FIRSTCONCPOS:= rtnn[, pmax(ID1_FIRSTPOSVIS,ID2_FIRSTPOSVIS)]]
	set(rtnn, NULL, c('ID1_FIRSTPOSVIS','ID1_FIRSTPOSDATE','ID2_FIRSTPOSVIS','ID2_FIRSTPOSDATE'), NULL)
	#	add in date of birth, comm type from nearest visit we have data on
	tmp		<- unique(subset(rd, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'ID1')
	tmp		<- unique(merge(rtnn, tmp, by='ID1')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('ID1','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'ID1', 'RID')
	tmp2	<- unique(subset(rd, select=c('RID','VISIT', "COMM_NUM", "COMM_NUM_A","BIRTHDATE")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), gsub('ID1_VISIT_FIRSTCONCPOS','VISIT_FIRSTCONCPOS',gsub('ID1_RID','ID1',paste0('ID1_',colnames(tmp)))))
	set(tmp, NULL, 'ID1_VISIT', NULL)	
	rtnn	<- merge(rtnn, tmp, by=c('ID1','VISIT_FIRSTCONCPOS'), all.x=1)	
	tmp		<- unique(subset(rd, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'ID2')
	tmp		<- unique(merge(rtnn, tmp, by='ID2')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('ID2','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'ID2', 'RID')
	tmp2	<- unique(subset(rd, select=c('RID','VISIT', "COMM_NUM", "COMM_NUM_A","BIRTHDATE")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), gsub('ID2_VISIT_FIRSTCONCPOS','VISIT_FIRSTCONCPOS',gsub('ID2_RID','ID2',paste0('ID2_',colnames(tmp)))))
	set(tmp, NULL, 'ID2_VISIT', NULL)	
	rtnn	<- merge(rtnn, tmp, by=c('ID2','VISIT_FIRSTCONCPOS'), all.x=1)
	#	get list of unique individuals in these events with metadata
	rtnn[, SELECT2:= 1L]
	set(rtnn, rtnn[, which(grepl('most likely a pair',SELECT) & SXO=='MF')], 'SELECT2', 2L)
	set(rtnn, rtnn[, which(grepl('mf|fm|12|21',SELECT) & SXO=='MF')], 'SELECT2', 3L)
	
	tmp	<- subset(rtnn, LINK_12==1 | LINK_21==1, select=c(ID1, ID1_SEX, ID1_COMM_NUM_A, ID1_BIRTHDATE, SELECT2))
	setnames(tmp, c('ID1','ID1_SEX','ID1_COMM_NUM_A','ID1_BIRTHDATE'), c('RID','SEX','COMM_NUM_A','BIRTHDATE'))	
	tmp2<- subset(rtnn, LINK_12==1 | LINK_21==1, select=c(ID2, ID2_SEX, ID2_COMM_NUM_A, ID2_BIRTHDATE, SELECT2))
	setnames(tmp2, c('ID2','ID2_SEX','ID2_COMM_NUM_A','ID2_BIRTHDATE'), c('RID','SEX','COMM_NUM_A','BIRTHDATE'))	
	tmp	<- rbind(tmp, tmp2)
	set(tmp, NULL, 'COMM_TYPE', tmp[, as.character(factor(substr(COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fishing','inland')))])
	set(tmp, tmp[, which(is.na(BIRTHDATE))],'BIRTHDATE',tmp[, mean(BIRTHDATE, na.rm=TRUE)])
	#	add if person is recent inmigrant
	tmp	<- merge(tmp, subset(de, !is.na(RID), c(RID, INMIG_2YRS_LOC, INMIG_2YRS)), by='RID', all.x=TRUE)
	set(tmp, tmp[, which(is.na(INMIG_2YRS_LOC))], 'INMIG_2YRS_LOC', 'resident')
	dr	<- tmp[, 	{
						z<- which.max(SELECT2)
						list(SEX=SEX[z], COMM_TYPE=COMM_TYPE[z], SELECT2=SELECT2[z], BIRTHDATE=BIRTHDATE[z], INMIG_2YRS_LOC=INMIG_2YRS_LOC[z])
					}, by='RID']
	dr[, AGE_AT_MID:= 2013.25-BIRTHDATE]
	dr[, AGE_AT_MID_C:= cut(AGE_AT_MID, breaks=c(10,25,35,65), labels=c('15-24','25-34','35+'), right=FALSE)]
	
	
	
	#
	# table by community type, sex
	#	
	ans		<- dr[, list(IN_TRM_CHAIN= length(RID[SELECT2>=1])), by=c('COMM_TYPE','SEX')]
	ans		<- merge(ans, dr[, list(LINKED= length(RID[SELECT2>=2])), by=c('COMM_TYPE','SEX')], by=c('COMM_TYPE','SEX'))
	ans		<- merge(ans, dr[, list(LINKED_DIR= length(RID[SELECT2>=3])), by=c('COMM_TYPE','SEX')], by=c('COMM_TYPE','SEX'))
	set(ans, ans[, which(grepl('fish',COMM_TYPE))], 'COMM_TYPE', 'fisherfolk')
	tmp		<- des[, list(	ELIGIBLE_1516= sum(PART_EVER+PART_NEVER), 
							PART_1516= sum(PART_EVER),
							HIV_1516= sum(HIV_1516_YES),
							ARTNAIVE= sum(HIV_1516_YES-SLART_AT_FIRST_VISIT),
							DEEP_SEQ_1516= sum(DEEP_SEQ_1516)
							), by=c('COMM_TYPE','SEX')]
	ans		<- merge(tmp, ans, by=c('COMM_TYPE','SEX'))	
	
	
	#
	# table by community type, sex, age
	#
	ansa	<- dr[, list(IN_TRM_CHAIN= length(RID[SELECT2>=1])), by=c('COMM_TYPE','SEX','AGE_AT_MID_C')]
	ansa	<- merge(ansa, dr[, list(LINKED= length(RID[SELECT2>=2])), by=c('COMM_TYPE','SEX','AGE_AT_MID_C')], by=c('COMM_TYPE','SEX','AGE_AT_MID_C'))
	ansa	<- merge(ansa, dr[, list(LINKED_DIR= length(RID[SELECT2>=3])), by=c('COMM_TYPE','SEX','AGE_AT_MID_C')], by=c('COMM_TYPE','SEX','AGE_AT_MID_C'))
	set(ansa, ansa[, which(grepl('fish',COMM_TYPE))], 'COMM_TYPE', 'fisherfolk')
	tmp		<- desa[, list(	ELIGIBLE_1516= sum(PART_EVER+PART_NEVER), 
							PART_1516= sum(PART_EVER),
							HIV_1516= sum(HIV_1516_YES),
							ARTNAIVE= sum(HIV_1516_YES-SLART_AT_FIRST_VISIT),
							DEEP_SEQ_1516= sum(DEEP_SEQ_1516)
							), by=c('COMM_TYPE','SEX','AGE_AT_MID_C')]
	ansa	<- merge(tmp, ansa, by=c('COMM_TYPE','SEX','AGE_AT_MID_C'), all=TRUE)
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
	
	#
	# table by community type, sex, migration status
	#
	ansm	<- dr[, list(IN_TRM_CHAIN= length(RID[SELECT2>=1])), by=c('COMM_TYPE','SEX','INMIG_2YRS_LOC')]
	ansm	<- merge(ansm, dr[, list(LINKED= length(RID[SELECT2>=2])), by=c('COMM_TYPE','SEX','INMIG_2YRS_LOC')], by=c('COMM_TYPE','SEX','INMIG_2YRS_LOC'))	
	set(ansm, ansm[, which(grepl('fish',COMM_TYPE))], 'COMM_TYPE', 'fisherfolk')
	tmp		<- desm[, list(	#ELIGIBLE_1516= sum(PART_EVER+PART_NEVER), 
							#PART_1516= sum(PART_EVER),
							HIV_1516= sum(HIV_1516_YES),
							ARTNAIVE= sum(HIV_1516_YES-SLART_AT_FIRST_VISIT),
							DEEP_SEQ_1516= sum(DEEP_SEQ_1516)
			), by=c('COMM_TYPE','SEX','INMIG_2YRS_LOC')]
	ansm	<- merge(tmp, ansm, by=c('COMM_TYPE','SEX','INMIG_2YRS_LOC'), all=TRUE)
	# add total to ansm
	ansm	<- melt(ansm, id.vars=c('COMM_TYPE','SEX','INMIG_2YRS_LOC'))
	set(ansm, ansm[, which(is.na(value))], 'value', 0L)
	tmp		<- ansm[, list(SEX=SEX, INMIG_2YRS_LOC=INMIG_2YRS_LOC, prop=value/sum(value), LEGEND=paste(value,' (',100*round(value/sum(value),d=2),'%)',sep='')), by=c('COMM_TYPE','variable')]
	ansm	<- merge(ansm, tmp, by=c('variable','COMM_TYPE','SEX','INMIG_2YRS_LOC'))	
	tmp		<- ansm[, list(SEX='Any', INMIG_2YRS_LOC='Any', value=sum(value), LEGEND=as.character(sum(value))), by=c('COMM_TYPE','variable')]
	tmp2	<- ansm[, list(COMM_TYPE='Any',SEX='Any', INMIG_2YRS_LOC='Any', value=sum(value), LEGEND=as.character(sum(value))), by=c('variable')]
	ansm	<- rbind(ansm,tmp,tmp2, fill=TRUE)
	ansm[, LABEL:= paste(COMM_TYPE, SEX, INMIG_2YRS_LOC, sep='-')]
	set(ansm, NULL, 'INMIG_2YRS_LOC', ansm[, factor(INMIG_2YRS_LOC, levels=c('Any','resident','inland','fisherfolk','external','origin_unknown'))])
	ansm	<- dcast.data.table(ansm, LABEL+COMM_TYPE+SEX+INMIG_2YRS_LOC~variable, value.var='LEGEND')
	setkey(ansm, COMM_TYPE, SEX, INMIG_2YRS_LOC)	
	write.csv(ansm, row.names=FALSE, file=paste0(outfile.base, '_table1_commtype_sex_migration.csv'))	
	
	
	#
	#	table by sex, age
	#
	ansa	<- dr[, list(IN_TRM_CHAIN= length(RID[SELECT2>=1])), by=c('SEX','AGE_AT_MID_C')]
	ansa	<- merge(ansa, dr[, list(LINKED= length(RID[SELECT2>=2])), by=c('SEX','AGE_AT_MID_C')], by=c('SEX','AGE_AT_MID_C'))
	ansa	<- merge(ansa, dr[, list(LINKED_DIR= length(RID[SELECT2>=3])), by=c('SEX','AGE_AT_MID_C')], by=c('SEX','AGE_AT_MID_C'))		
	tmp		<- desa[, list(	ELIGIBLE_1516= sum(PART_EVER+PART_NEVER), 
							PART_1516= sum(PART_EVER),
							HIV_1516= sum(HIV_1516_YES),
							ARTNAIVE= sum(HIV_1516_YES-SLART_AT_FIRST_VISIT),
							DEEP_SEQ_1516= sum(DEEP_SEQ_1516)
					), by=c('SEX','AGE_AT_MID_C')]
	ansa	<- merge(tmp, ansa, by=c('SEX','AGE_AT_MID_C'))
	# add proportions to ansa
	ansa	<- melt(ansa, id.vars=c('SEX','AGE_AT_MID_C'))
	tmp		<- ansa[, list(SEX=SEX, AGE_AT_MID_C=AGE_AT_MID_C, prop=value/sum(value), LEGEND=paste(value,' (',100*round(value/sum(value),d=2),'%)',sep='')), by=c('variable')]
	ansa	<- merge(ansa, tmp, by=c('variable','SEX','AGE_AT_MID_C'))
	# add total to ansa
	tmp		<- ansa[, list(SEX='Any', AGE_AT_MID_C='Any', value=sum(value), LEGEND=as.character(sum(value))), by=c('variable')]
	tmp2	<- ansa[, list(AGE_AT_MID_C='Any', value=sum(value), LEGEND=as.character(sum(value))), by=c('variable','SEX')]
	ansa	<- rbind(ansa,tmp,tmp2, fill=TRUE)
	ansa[, LABEL:= paste(SEX, AGE_AT_MID_C, sep='-')]
	ansa	<- dcast.data.table(ansa, LABEL+SEX+AGE_AT_MID_C~variable, value.var='LEGEND')
	set(ansa, NULL, 'LABEL', ansa[, factor(LABEL, levels=c("Any-Any","F-Any","F-15-24","F-25-34","F-35+","M-Any","M-15-24","M-25-34","M-35+"))])
	setkey(ansa, LABEL)
	write.csv(ansa, row.names=FALSE, file=paste0(outfile.base, '_table1_nocommtype.csv'))	
}


RakaiFull.phylogeography.171122.recpients.on.map<- function()
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

	
	subset(rtr2, TR_INMIGRANT==1)	
}

RakaiFull.phylogeography.180618.figure.flows.on.map<- function()
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
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda"		
	outfile.base			<- gsub('_data_with_inmigrants.rda','',infile)
	load(infile)
	setnames(rtr2, c('TR_INMIGRATE_2YR','REC_INMIGRATE_2YR'), c('TR_INMIGRATE','REC_INMIGRATE'))
	
	#	define location from where inmigrants arrived 
	rtr3	<- subset(rtr2, grepl('inmigrant',TR_INMIGRATE) & !grepl('external',TR_INMIGRATE))
	

	#	update inmigrant with new resolved locations in inmigrant2
	infile		<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData_v2.rda"
	load(infile)
	inmigrant	<- as.data.table(inmigrant)
	infile		<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/migrants_withMissingGPS.csv"
	inmigrant2	<- as.data.table(read.csv(infile))
	inmigrant	<- merge(inmigrant, subset(inmigrant2, select=c('RCCS_studyid','visit','X')), by=c('RCCS_studyid','visit'), all.x=TRUE)
	tmp			<- subset(inmigrant, !is.na(X), c(RCCS_studyid, visit, date, inmigrant))
	inmigrant2	<- merge(inmigrant2, tmp, by=c('RCCS_studyid','visit'))
	inmigrant	<- subset(inmigrant, is.na(X))
	set(inmigrant, NULL, 'X', NULL)
	set(inmigrant2, NULL, c('X','inmig_place_original','recode_city'), NULL)
	inmigrant	<- rbind(inmigrant, inmigrant2)
	#	prepare inmigrant -- identify inmigrants from fishing communities and from external
	set(inmigrant, NULL, 'inmig_place', inmigrant[, gsub('DDIMO|DDIMU|DIMO|DIMU','DIMU',inmig_place)])
	set(inmigrant, inmigrant[, which(grepl('MALEMBO',inmig_place))], 'inmig_place', 'MALEMBO')
	set(inmigrant, NULL, 'inmig_place', inmigrant[, gsub("KASEMSERO","KASENSERO",inmig_place)])
	set(inmigrant, inmigrant[, which(grepl('KASENSERO',inmig_place))], 'inmig_place', 'KASENSERO')
	set(inmigrant, NULL, 'date', inmigrant[, hivc.db.Date2numeric(date)])	
	#	define from_fishing and from_outside and from_inland
	inmigrant[, INMIG_LOC:= 'inland' ]
	set(inmigrant, inmigrant[, which(grepl('MALEMBO|DIMU|KASENSERO|NAMIREMBE',inmig_place))], 'INMIG_LOC','fisherfolk')
	set(inmigrant, inmigrant[, which(inmig_admin0!='Uganda')], 'INMIG_LOC','external')
	set(inmigrant, inmigrant[, which(inmig_admin1!='Rakai')], 'INMIG_LOC','external')
	set(inmigrant, inmigrant[, which(is.na(inmig_admin1))], 'INMIG_LOC','unknown')	
	setnames(inmigrant, c('RCCS_studyid','date'), c('TR_RID','TR_VISIT_DATE'))	
	rtr3	<- merge(unique(subset(rtr3, select=c(TR_RID, REC_RID, VISIT_FIRSTCONCPOS, DATE_FIRSTCONCPOS))), inmigrant, by='TR_RID', all.x=TRUE)
	setnames(rtr3, c('inmig_lon','inmig_lat','INMIG_LOC'), c('TR_INMIG_LON','TR_INMIG_LAT','TR_INMIG_LOC'))
	#	find inmigration event of transmitter that is at or before first time both were recorded infected
	rtr3	<- subset(rtr3, DATE_FIRSTCONCPOS >= TR_VISIT_DATE)
	tmp		<- rtr3[,  list(TR_VISIT_DATE=TR_VISIT_DATE[which.min(DATE_FIRSTCONCPOS-TR_VISIT_DATE)]), by=c('TR_RID','REC_RID','DATE_FIRSTCONCPOS')]
	rtr3	<- merge(rtr3,tmp,by=c('TR_RID','REC_RID','DATE_FIRSTCONCPOS','TR_VISIT_DATE'))	
	rtr3	<- subset(rtr3, select=c(TR_RID, REC_RID, VISIT_FIRSTCONCPOS, TR_INMIG_LON, TR_INMIG_LAT, TR_INMIG_LOC))
	#	add to rtr2
	rtr2	<- merge(rtr2,rtr3,by=c('TR_RID','REC_RID','VISIT_FIRSTCONCPOS'),all.x=TRUE)
	
	#	add coordinates to communities
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, LONG, LAT)))
	setnames(tmp, c('COMM_NUM_A','LONG','LAT'), c('TR_COMM_NUM_A','TR_X','TR_Y'))
	rtr2	<- merge(rtr2, tmp, by='TR_COMM_NUM_A')
	setnames(tmp, c('TR_COMM_NUM_A','TR_X','TR_Y'), c('REC_COMM_NUM_A','REC_X','REC_Y'))
	rtr2	<- merge(rtr2, tmp, by='REC_COMM_NUM_A')
	
	#
	#	get map
	#
	style	<- "feature:road|color:0x17202A&style=feature:water|color:0x677996&style=feature:landscape.natural|color:0xedecda&style=feature:administrative|visibility=off"
	zm		<- get_googlemap(c(lon=31.65, lat=-0.66), scale=2, size=c(550,550), zoom=10, maptype="road", style=style)
	zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))			
	
	#
	#	to plot flows from migrants, set up fake communities from migrant locations
	#
	#	set up dummy location in corner for inmigration from external 
	#tmp		<- rtr2[, c(min(c(REC_X, TR_X)), max(c(REC_Y, TR_Y)))]
	tmp		<- c(31.3,-0.34)
	tmp2	<- rtr2[, which(grepl('external',TR_INMIGRATE))]
	set(rtr2, tmp2, 'TR_INMIG_LON', tmp[1])
	set(rtr2, tmp2, 'TR_INMIG_LAT', tmp[2])
	#	set up dummy location in corner for inmigration from unknown
	#tmp		<- rtr2[, c(max(c(REC_X, TR_X)), min(c(REC_Y, TR_Y)))]
	tmp		<- c(32.0,-1.0)
	tmp2	<- rtr2[, which(TR_INMIG_LOC=='unknown')]
	set(rtr2, tmp2, 'TR_INMIG_LON', tmp[1])
	set(rtr2, tmp2, 'TR_INMIG_LAT', tmp[2])
	#	define fake community IDs
	tmp		<- rtr2[, which(grepl('external',TR_INMIGRATE))]
	set(rtr2, tmp, 'TR_INMIG_LOC', 'external')	
	tmp		<- unique(subset(rtr2, !grepl('resident',TR_INMIGRATE), c(TR_INMIG_LON, TR_INMIG_LAT, TR_INMIG_LOC)))
	tmp[, TR_COMM_NUM_A_MIG:= paste0(substr(TR_INMIG_LOC,1,1),'mig',seq_len(nrow(tmp)))]
	rtr2	<- merge(rtr2, tmp, by=c('TR_INMIG_LON','TR_INMIG_LAT','TR_INMIG_LOC'), all.x=TRUE)
	
	#
	#	change locations of transmitter if transmitter is inmigrant and flow is not between same areas
	#
	tmp		<- rtr2[, which(TR_INMIGRATE!='resident')]	
	set(rtr2, tmp, 'TR_COMM_NUM_A', rtr2[tmp,TR_COMM_NUM_A_MIG])
	set(rtr2, tmp, 'TR_X', rtr2[tmp,TR_INMIG_LON])
	set(rtr2, tmp, 'TR_Y', rtr2[tmp,TR_INMIG_LAT])
	 	
	
	#
	#	plot number of observed transmission flows
	#
	
	#	overall parameters for plotting
	tr.edge.gap	<- 0
	rec.edge.gap<- 0.06
	node.size	<- 0.2
	edge.size	<- 1
	curvature	<- -0.2
	label.size	<- 5
	arrow		<- arrow(length=unit(0.04, "npc"), type="closed", angle=15)
	arrow2		<- arrow(length=unit(0.02, "npc"), type="closed", angle=15)
	curv.shift	<- 0.08
	#
	#	define flows
	#
	dc		<- rtr2[, list(FLOW=length(PAIRID)), by=c('TR_COMM_NUM_A','TR_INMIGRATE','REC_COMM_NUM_A','TR_X','TR_Y','REC_X','REC_Y')]
	tmp		<- dc[, which(substr(TR_COMM_NUM_A,1,1)=='u')]
	set(dc, tmp, 'TR_INMIGRATE', dc[tmp, gsub('inland|fisherfolk','unknownloc',TR_INMIGRATE)])
	set(dc, dc[, which(TR_INMIGRATE=='resident' & substr(TR_COMM_NUM_A,1,1)=='f')], 'TR_INMIGRATE', 'resident_fish')
	set(dc, dc[, which(TR_INMIGRATE=='resident' & substr(TR_COMM_NUM_A,1,1)!='f')], 'TR_INMIGRATE', 'resident_inland')
	#	count flows
	subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)=='f')[, sum(FLOW)]
	subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(TR_COMM_NUM_A,1,1)!='e' & substr(TR_COMM_NUM_A,1,1)!='u' & substr(REC_COMM_NUM_A,1,1)!='f')[, sum(FLOW)]
	subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)!='f')
	subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(TR_COMM_NUM_A,1,1)!='e' & substr(TR_COMM_NUM_A,1,1)!='u' & substr(REC_COMM_NUM_A,1,1)=='f')
	
	#
	#	define stuff for plotting curved flows
	#	
	#dc[, which(substr(TR_COMM_NUM_A,1,1)=='a' & substr(REC_COMM_NUM_A,1,1)=='a' & TR_INMIGRATE!='resident')]
	set(dc, NULL, 'TR_COMM_NUM_A', dc[, as.character(TR_COMM_NUM_A)])
	set(dc, NULL, 'REC_COMM_NUM_A', dc[, as.character(REC_COMM_NUM_A)])
	dc[, EDGETEXT_X:= (TR_X+REC_X)/2]
	dc[, EDGETEXT_Y:= (TR_Y+REC_Y)/2]
	dc[, EDGE_LABEL:= FLOW ]
	dc[, MX:= (REC_X - TR_X)]	
	dc[, MY:= (REC_Y - TR_Y)]	
	set(dc, NULL, 'TR_X_EDGE', dc[, TR_X + MX*tr.edge.gap])
	set(dc, NULL, 'TR_Y_EDGE', dc[, TR_Y + MY*tr.edge.gap])
	set(dc, NULL, 'REC_X_EDGE', dc[, REC_X - MX*rec.edge.gap])
	set(dc, NULL, 'REC_Y_EDGE', dc[, REC_Y - MY*rec.edge.gap])		
	dc[, TX:= -MY]
	dc[, TY:= MX]
	set(dc, NULL, 'EDGETEXT_X', dc[, EDGETEXT_X + TX*curv.shift])
	set(dc, NULL, 'EDGETEXT_Y', dc[, EDGETEXT_Y + TY*curv.shift])
	#
	#	flows inland -> fishing
	#
	ggmap(zm) +		
			geom_point(	data=unique(subset(dc,  substr(REC_COMM_NUM_A,1,1)=='f' & FLOW>0, c(REC_X, REC_Y))), 
						aes(x=REC_X, y=REC_Y), 
						colour='firebrick1', size=4) +
			geom_curve(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)=='f' & FLOW>0), 
						aes(x=TR_X_EDGE, xend=REC_X_EDGE, y=TR_Y_EDGE, yend=REC_Y_EDGE, colour=TR_INMIGRATE), 
						curvature=curvature, arrow=arrow, lineend="butt", linejoin='mitre') +
			geom_text(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)=='f' & FLOW>0), 
						aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL, colour=TR_INMIGRATE), 
						size=label.size) +
			coord_cartesian() +
			scale_colour_manual(values=c('resident_inland'="#66BD63", 'inmigrant_from_inland'="#00441B",'inmigrant_from_unknownloc'='black', 'inmigrant_from_external'='deepskyblue1')) +
			#scale_size(breaks=c(1,2,4), range=c(0.6,1.2)) +
			labs(x='', y='') +
			theme(legend.position='bottom', legend.box = "vertical") + 
			guides(size='none', colour='none') 
	ggsave(file=paste0(outfile.base,'_flows_intofishing_on_map.pdf'), w=7, h=7, useDingbats=FALSE)
	#
	#	flows fishing -> inland
	#		
	ggmap(zm) +					
			geom_point(	data=unique(subset(dc, substr(REC_COMM_NUM_A,1,1)!='f' & FLOW>0, c(REC_X, REC_Y))), 
						aes(x=REC_X, y=REC_Y), 
						colour="#00441B", size=4) +
			geom_curve(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)!='f' & FLOW>0), 
						aes(x=TR_X_EDGE, xend=REC_X_EDGE, y=TR_Y_EDGE, yend=REC_Y_EDGE, colour=TR_INMIGRATE), 
						curvature=curvature, arrow=arrow, lineend="butt", linejoin='mitre') +
			geom_text(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f'  & substr(REC_COMM_NUM_A,1,1)!='f' & FLOW>0), 
						aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL, colour=TR_INMIGRATE), 
						size=label.size) +
			coord_cartesian() +
			#scale_size(breaks=c(1,2,4), range=c(0.6,1.2))+
			scale_colour_manual(values=c('resident_fish'="firebrick1", 'inmigrant_from_fisherfolk'="darkred",'inmigrant_from_unknownloc'='black', 'inmigrant_from_external'='deepskyblue1')) +
			labs(x='', y='') +
			theme(legend.position='bottom', legend.box = "vertical") + 
			guides(size='none', colour='none')
	ggsave(file=paste0(outfile.base,'_flows_fromfishing_on_map.pdf'), w=7, h=7, useDingbats=FALSE)
	#
	#	flows fishing -> fishing and inland -> inland
	#
	ggmap(zm) +		
			geom_point(	data=unique(subset(dc,  substr(REC_COMM_NUM_A,1,1)=='f' & FLOW>0, c(REC_X, REC_Y))), 
					aes(x=REC_X, y=REC_Y), 
					colour='firebrick1', size=4) +
			geom_curve(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A!=TR_COMM_NUM_A & FLOW>0), 
					aes(x=TR_X_EDGE, xend=REC_X_EDGE, y=TR_Y_EDGE, yend=REC_Y_EDGE, colour=TR_INMIGRATE), 
					curvature=curvature, arrow=arrow, lineend="butt", linejoin='mitre') +			
			geom_curve(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A==TR_COMM_NUM_A & FLOW>0), 
					aes(x=TR_X+0.005, xend=REC_X+0.03, y=TR_Y-0.005, yend=REC_Y, colour=TR_INMIGRATE), 
					curvature=1) +
			geom_curve(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A==TR_COMM_NUM_A & FLOW>0), 
					aes(x=TR_X+0.03, xend=REC_X+0.005, y=TR_Y, yend=REC_Y+0.005, colour=TR_INMIGRATE), 
					curvature=1, arrow=arrow2, lineend="butt", linejoin='mitre') +
			geom_text(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)=='f' & FLOW>0), 
					aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL, colour=TR_INMIGRATE), 
					size=label.size) +
			coord_cartesian() +
			scale_colour_manual(values=c('resident_inland'="#66BD63",'resident_fish'="firebrick1", 'inmigrant_from_inland'="#00441B",'inmigrant_from_fisherfolk'="darkred",'inmigrant_from_unknownloc'='black', 'inmigrant_from_external'='deepskyblue1')) +
			labs(x='', y='') +
			theme(legend.position='bottom', legend.box = "vertical") + 
			guides(size='none', colour='none')
	ggsave(file=paste0(outfile.base,'_flows_fishing_on_map.pdf'), w=7, h=7, useDingbats=FALSE)
	ggmap(zm) +		
			geom_point(	data=unique(subset(dc,  substr(REC_COMM_NUM_A,1,1)!='f' & FLOW>0, c(REC_X, REC_Y))), 
					aes(x=REC_X, y=REC_Y), 
					colour="#00441B", size=4) +
			geom_curve(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)!='f' & REC_COMM_NUM_A!=TR_COMM_NUM_A & FLOW>0), 
					aes(x=TR_X_EDGE, xend=REC_X_EDGE, y=TR_Y_EDGE, yend=REC_Y_EDGE, colour=TR_INMIGRATE), 
					curvature=curvature, arrow=arrow, lineend="butt", linejoin='mitre') +			
			geom_curve(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)!='f' & REC_COMM_NUM_A==TR_COMM_NUM_A & FLOW>0), 
					aes(x=TR_X+0.005, xend=REC_X+0.03, y=TR_Y-0.005, yend=REC_Y, colour=TR_INMIGRATE), 
					curvature=1) +
			geom_curve(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)!='f' & REC_COMM_NUM_A==TR_COMM_NUM_A & FLOW>0), 
					aes(x=TR_X+0.03, xend=REC_X+0.005, y=TR_Y, yend=REC_Y+0.005, colour=TR_INMIGRATE), 
					curvature=1, arrow=arrow2, lineend="butt", linejoin='mitre') +
			geom_text(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)!='f' & FLOW>0), 
					aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL, colour=TR_INMIGRATE), 
					size=label.size) +			
			coord_cartesian() +
			scale_colour_manual(values=c('resident_inland'="#66BD63",'resident_fish'="firebrick1", 'inmigrant_from_inland'="#00441B",'inmigrant_from_fisherfolk'="darkred",'inmigrant_from_unknownloc'='black', 'inmigrant_from_external'='deepskyblue1')) +
			labs(x='', y='') +
			theme(legend.position='bottom', legend.box = "vertical") + 
			guides(size='none', colour='none')
	ggsave(file=paste0(outfile.base,'_flows_inland_on_map.pdf'), w=7, h=7, useDingbats=FALSE)
	ggmap(zm) +		
			geom_point(	data=unique(subset(dc,  substr(REC_COMM_NUM_A,1,1)=='f' & FLOW>0, c(REC_X, REC_Y))), 
						aes(x=REC_X, y=REC_Y), 
						colour='firebrick1', size=4) +
			geom_curve(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A!=TR_COMM_NUM_A & FLOW>0), 
						aes(x=TR_X_EDGE, xend=REC_X_EDGE, y=TR_Y_EDGE, yend=REC_Y_EDGE, colour=TR_INMIGRATE), 
						curvature=curvature, arrow=arrow, lineend="butt", linejoin='mitre') +			
			geom_curve(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A==TR_COMM_NUM_A & FLOW>0), 
						aes(x=TR_X+0.005, xend=REC_X+0.03, y=TR_Y-0.005, yend=REC_Y, colour=TR_INMIGRATE), 
						curvature=1) +
			geom_curve(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)=='f' & REC_COMM_NUM_A==TR_COMM_NUM_A & FLOW>0), 
						aes(x=TR_X+0.03, xend=REC_X+0.005, y=TR_Y, yend=REC_Y+0.005, colour=TR_INMIGRATE), 
						curvature=1, arrow=arrow2, lineend="butt", linejoin='mitre') +
			geom_text(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)=='f' & substr(REC_COMM_NUM_A,1,1)=='f' & FLOW>0), 
						aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL, colour=TR_INMIGRATE), 
						size=label.size) +
			geom_point(	data=unique(subset(dc,  substr(REC_COMM_NUM_A,1,1)!='f' & FLOW>0, c(REC_X, REC_Y))), 
						aes(x=REC_X, y=REC_Y), 
						colour="#00441B", size=4) +
			geom_curve(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)!='f' & REC_COMM_NUM_A!=TR_COMM_NUM_A & FLOW>0), 
						aes(x=TR_X_EDGE, xend=REC_X_EDGE, y=TR_Y_EDGE, yend=REC_Y_EDGE, colour=TR_INMIGRATE), 
						curvature=curvature, arrow=arrow, lineend="butt", linejoin='mitre') +			
			geom_curve(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)!='f' & REC_COMM_NUM_A==TR_COMM_NUM_A & FLOW>0), 
						aes(x=TR_X+0.005, xend=REC_X+0.03, y=TR_Y-0.005, yend=REC_Y, colour=TR_INMIGRATE), 
						curvature=1) +
			geom_curve(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)!='f' & REC_COMM_NUM_A==TR_COMM_NUM_A & FLOW>0), 
						aes(x=TR_X+0.03, xend=REC_X+0.005, y=TR_Y, yend=REC_Y+0.005, colour=TR_INMIGRATE), 
						curvature=1, arrow=arrow2, lineend="butt", linejoin='mitre') +
			geom_text(	data=subset(dc, substr(TR_COMM_NUM_A,1,1)!='f' & substr(REC_COMM_NUM_A,1,1)!='f' & FLOW>0), 
						aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL, colour=TR_INMIGRATE), 
						size=label.size) +			
			coord_cartesian() +
			scale_colour_manual(values=c('resident_inland'="#66BD63",'resident_fish'="firebrick1", 'inmigrant_from_inland'="#00441B",'inmigrant_from_fisherfolk'="darkred",'inmigrant_from_unknownloc'='black', 'inmigrant_from_external'='deepskyblue1')) +
			labs(x='', y='') +
			theme(legend.position='bottom', legend.box = "vertical") + 
			guides(size='none', colour='none')
	ggsave(file=paste0(outfile.base,'_flows_fishinginland_on_map.pdf'), w=7, h=7, useDingbats=FALSE)
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

RakaiFull.phylogeography.180322.core.inference.stan.models.noseqsampling<- function()
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
	require(rethinking)
	
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
	rtpdm[, MALE_FISHCOMM:= rtpdm[, as.character(factor(MALE_COMM_TYPE=='fisherfolk', levels=c(TRUE,FALSE), labels=c('Fishing','Inland')))]]
	rtpdm[, FEMALE_FISHCOMM:= rtpdm[, as.character(factor(FEMALE_COMM_TYPE=='fisherfolk', levels=c(TRUE,FALSE), labels=c('Fishing','Inland')))]]
	
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
	
	#	rtr2 is a line list of transmitters -> recipients ( no time )
	#	sum transmission events by community and gender
	dc1	<- rtr2[, list(TR_OBS=length(PAIRID)), by=c('TR_FISHCOMM','REC_FISHCOMM')]
	dc1[, FROM_TO_ID:= seq_len(nrow(dc1))]
	
	#	write stan model 1 between community types
	m1		<- list()
	m1$code	<- "data {
	  int<lower=0> K; 	// number of community combinations 
	  int TR_OBS[K]; 	// transmission counts	   
	}
	transformed data {
	  vector<lower=0>[K] dirprior_comm_flow= rep_vector(1, K); // community flows, positive
	}
	parameters {
	  simplex[K] prop_comm_flow;	//prior on community flows, sum to 1  	  
	}
	model {
	  //dirprior_comm_flow ~ lognormal(0, 5);
	  prop_comm_flow ~ dirichlet(dirprior_comm_flow);	  
	  TR_OBS ~ multinomial(prop_comm_flow);	  
	}"
	m1$code		<- gsub('\t','',m1$code)
	m1$codefile	<- paste0(outfile.base,'_stanm1.stan')
	#write(m1$code, file=m1$codefile)	
	m1$data			<- as.list(dc1)
	m1$data[['K']] 	<- nrow(dc1)	
	m1$fit			<- stan(model_code=m1$code, data=m1$data, warmup=5e2, iter=1e4, chains=1)
	#m1.fit			<- stan(file=m1$codefile, data=m1$data, iter = 1000, chains = 4)
	print(m1.fit)
	plot(m1.fit)	
	traceplot(m1$fit)
	pairs(m1$fit)
	
	#	write stan model 2 between communities with distance based prior
	dc5	<- rtr2[, list(TR_OBS=length(PAIRID)), by=c('TR_FISHCOMM','REC_FISHCOMM','TR_COMM_NUM_A','REC_COMM_NUM_A')]
	tmp	<- unique(subset(ds, select=c(COMM_NUM_A, LONG, LAT)))
	tmp	<- merge(tmp, data.table(COMM_NUM_A= sort(unique(tmp$COMM_NUM_A)), COMM_NUM2= seq_along(unique(tmp$COMM_NUM_A))), by='COMM_NUM_A')
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc5	<- merge(dc5, tmp, by='TR_COMM_NUM_A')
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc5	<- merge(dc5, tmp, by='REC_COMM_NUM_A')
	dc5[, EUCLIDEAN_DIST:= sqrt( (TR_LONG-REC_LONG)^2+(TR_LAT-REC_LAT)^2 )]
	dc5[, FROM_TO_ID:= seq_len(nrow(dc5))]
	
	#m3 has convergence issues: sig_random_pair_effect. coeff_space and base_flow are also extremely highly correlated.
	m2		<- list()
	m2$code	<- "data {
	  int<lower=1> N; 						// number of observations
	  int<lower=1> N_PAIRS;					// number of **observed** community pairs 
	  int TR_OBS[N]; 						// transmission counts
	  vector<lower=0>[N] EUCLIDEAN_DIST; 	// distances for fixed spatial effect
	  int<lower=1> FROM_TO_ID[N]; 			// community pair for ith observation	  
	}
	parameters {
	  simplex[N_PAIRS] prop_comm_flow;  	// community flows, sum to 1
	  vector[N_PAIRS] random_pair_effect;	// random effect on flows between two communities 	  
	  real base_flow;						// overall grand mean
	  real coeff_space;						// strength of spatial effect
	  real<lower=0> rhosq;	 				// rate of decline of spatical effect
	  real<lower=0> sig_random_pair_effect; // hyperprior for random pair effect 	  
	}
	model {	  	  
	  vector[N] fixed_spatial_effect;	// fixed spatial effect
	  vector[N] log_lambda;				// linear predictor
	  base_flow ~ normal(0, 100);		// prior on log mean flow between all communities
	  coeff_space ~ normal(0, 10 );		// prior on strength of spatial effect
	  rhosq ~ exponential( 1 );			// prior on rate of decline of spatical effect
	  random_pair_effect ~ normal(0, sig_random_pair_effect);
	  sig_random_pair_effect ~ exponential( 1 );
	  // build linear predictor
	  for(i in 1:N)	    
	  {
	     fixed_spatial_effect[i]= coeff_space*exp(-rhosq*pow(EUCLIDEAN_DIST[i],2));
	     log_lambda[i]= base_flow + random_pair_effect[FROM_TO_ID[i]] + fixed_spatial_effect[i];	      
	  }	  
	  prop_comm_flow ~ dirichlet_log(log_lambda);	//dirichlet prior to satisfy sum-to-1 constraint
	  TR_OBS ~ multinomial(prop_comm_flow);	  
	}"
	m2$code		<- gsub('\t','',m2$code)
	m2$codefile	<- paste0(outfile.base,'_stanm2.stan')
	#write(m1$code, file=m1$codefile)	
	m2$data					<- as.list(dc5)
	m2$data[['N']] 			<- nrow(dc5)
	m2$data[['N_PAIRS']] 	<- max(dc5$FROM_TO_ID)
	m2$fit		<- stan(model_code=m2$code, 
						data=m2$data, 
						warmup=5e2, iter=1e3, chains=1,
						init=list(list(base_flow=exp(1), coeff_space=0, rhosq=1, sig_random_pair_effect=1, random_pair_effect=rep(0,m2$data[['N_PAIRS']]))))				
	precis(m2$fit, depth=2, prob=0.95)
	post <- as.data.frame(m2$fit)
	pdf(file=paste0(outfile.base,'_stanm2_traces.pdf'), w=15, h=100)
	traceplot(m2$fit, pars=colnames(post), ncol=2)
	dev.off()
	pairs(m2$fit, pars=c('base_flow','coeff_space','rhosq'))
	
	#m3 has convergence issues: sig_random_pair_effect
	m3		<- list()
	m3$code	<- "data {
	  int<lower=1> N; 					// number of observations
	  int<lower=1> N_PAIRS;				// number of **observed** community pairs 
	  int TR_OBS[N]; 					// transmission counts
	  int<lower=1> FROM_TO_ID[N]; 		// community pair for ith observation	  
	}
	parameters {
	  simplex[N_PAIRS] prop_comm_flow;  	// community flows, sum to 1
	  vector[N_PAIRS] random_pair_effect;	// random effect on flows between two communities 	  
	  real base_flow;						// overall grand mean	  
	  real<lower=0> sig_random_pair_effect; // hyperprior for random pair effect	   	  
	}
	model {	  	  
	  vector[N] log_lambda;				// linear predictor
	  base_flow ~ normal(0, 100);		// prior on log mean flow between all communities
	  random_pair_effect ~ normal(0, sig_random_pair_effect);
	  sig_random_pair_effect ~ exponential( 1 );
	  // build linear predictor
	  for(i in 1:N)	    
	  {	  
	  	log_lambda[i]= base_flow + random_pair_effect[FROM_TO_ID[i]];	      
	  }	  
	  prop_comm_flow ~ dirichlet_log(log_lambda);	//dirichlet prior to satisfy sum-to-1 constraint
	  TR_OBS ~ multinomial(prop_comm_flow);	  
	}"
	m3$code		<- gsub('\t','',m3$code)
	m3$codefile	<- paste0(outfile.base,'_stanm3.stan')
	#write(m3$code, file=m1$codefile)	
	m3$data					<- as.list(dc5)
	m3$data[['N']] 			<- nrow(dc5)
	m3$data[['N_PAIRS']] 	<- max(dc5$FROM_TO_ID)
	m3$fit		<- stan(model_code=m3$code, 
			data=m2$data, 
			warmup=5e2, iter=1e3, chains=1,
			init=list(list(base_flow=exp(1), sig_random_pair_effect=1, random_pair_effect=rep(0,m2$data[['N_PAIRS']]))))				
	precis(m3$fit, depth=2, prob=0.95)
	post <- as.data.frame(m3$fit)
	pdf(file=paste0(outfile.base,'_stanm3_traces.pdf'), w=7, h=150)
	traceplot(m3$fit, pars=colnames(post), ncol=1)
	dev.off()
	
	#	this works!
	m4		<- list()
	m4$code	<- "data {
	  int<lower=1> N; 						// number of observations
	  int<lower=1> N_PAIRS;					// number of **observed** community pairs 
	  int TR_OBS[N]; 						// transmission counts	  	  
	}
	transformed data {
	  vector<lower=0>[N_PAIRS] dirprior_comm_flow= rep_vector(1, N_PAIRS); // community flows, positive
	}	
	parameters {
	  simplex[N_PAIRS] prop_comm_flow;  	// community flows, sum to 1	   	  
	}
	model {	  	  
	  prop_comm_flow ~ dirichlet_log(dirprior_comm_flow);	//dirichlet prior to satisfy sum-to-1 constraint
	  TR_OBS ~ multinomial(prop_comm_flow);	  
	}"			
	m4$code		<- gsub('\t','',m4$code)
	m4$codefile	<- paste0(outfile.base,'_stanm4.stan')
	#write(m1$code, file=m1$codefile)	
	m4$data					<- as.list(dc5)
	m4$data[['N']] 			<- nrow(dc5)
	m4$data[['N_PAIRS']] 	<- max(dc5$FROM_TO_ID)
	m4$fit		<- stan(model_code=m4$code, 
			data=m4$data, 
			warmup=5e2, iter=1e4, chains=1,
			init=list(list(base_flow=exp(1), coeff_space=0, rhosq=1, sig_random_pair_effect=1, random_pair_effect=rep(0,m4$data[['N_PAIRS']]))))				
	precis(m4$fit, depth=2, prob=0.95)
	post <- as.data.frame(m4$fit)
	pdf(file=paste0(outfile.base,'_stanm4_traces.pdf'), w=10, h=100)
	traceplot(m4$fit, pars=colnames(post), ncol=1)
	dev.off()
	
	
	#	this works!
	m5		<- list()
	m5$code	<- "data {
	  int<lower=1> N; 						// number of observations
	  int<lower=1> N_PAIRS;					// number of **observed** community pairs 
	  int TR_OBS[N]; 						// transmission counts	  	  
	}
	transformed data {
	  real<lower=0> objective_prior_weight;
	  vector<lower=0>[N_PAIRS] dirprior_comm_flow;
	  objective_prior_weight= N_PAIRS;
	  objective_prior_weight= 0.8/objective_prior_weight; //this is the Berger objective prior with minimal loss compared to marginal Beta reference prior	  
	  dirprior_comm_flow= rep_vector(objective_prior_weight, N_PAIRS); // community flows, positive
	}	
	parameters {
	  simplex[N_PAIRS] prop_comm_flow;  	// community flows, sum to 1	   	  
	}
	model {	  	  
	  prop_comm_flow ~ dirichlet_log(dirprior_comm_flow);	//dirichlet prior to satisfy sum-to-1 constraint
	  TR_OBS ~ multinomial(prop_comm_flow);	  
	}"			
	m5$code		<- gsub('\t','',m5$code)
	m5$codefile	<- paste0(outfile.base,'_stanm5.stan')
	#write(m1$code, file=m1$codefile)	
	m5$data					<- as.list(dc5)
	m5$data[['N']] 			<- nrow(dc5)
	m5$data[['N_PAIRS']] 	<- max(dc5$FROM_TO_ID)
	m5$fit		<- stan(model_code=m5$code, 
			data=m5$data, 
			warmup=5e2, iter=1e4, chains=1,
			init=list(list(base_flow=exp(1), coeff_space=0, rhosq=1, sig_random_pair_effect=1, random_pair_effect=rep(0,m5$data[['N_PAIRS']]))))				
	precis(m5$fit, depth=2, prob=0.95)
	post <- as.data.frame(m5$fit)
	pdf(file=paste0(outfile.base,'_stanm5_traces.pdf'), w=10, h=100)
	traceplot(m5$fit, pars=colnames(post), ncol=1)
	dev.off()
	
	#	save model fits
	save(dc1, dc5, m1, m4, m5, file=paste0(outfile.base,'_stanm1m4m5.rda'))

	#	question: how different are estimates from m2 m4 m5?
	post	<- extract.samples(m1$fit)
	tmp		<- data.table(FROM_TO_ID=1:max(dc1$FROM_TO_ID))
	df1		<- tmp[, list(MCIT=seq_len(nrow(post$prop_comm_flow)), PI=post$prop_comm_flow[,FROM_TO_ID]), by='FROM_TO_ID']
	df1		<- merge(dc1, df1, by='FROM_TO_ID')	
	df1[, MODEL:='area-based']
		
	post	<- extract.samples(m4$fit)
	tmp		<- data.table(FROM_TO_ID=1:max(dc5$FROM_TO_ID))
	df4		<- tmp[, list(MCIT=seq_len(nrow(post$prop_comm_flow)), PI=post$prop_comm_flow[,FROM_TO_ID]), by='FROM_TO_ID']
	df4		<- merge(dc5, df4, by='FROM_TO_ID')
	df4		<- df4[, list(PI=sum(PI)), by=c('MCIT','TR_FISHCOMM','REC_FISHCOMM')]
	df4[, MODEL:='community-based prior Dir(1)']
	
	post	<- extract.samples(m5$fit)
	tmp		<- data.table(FROM_TO_ID=1:max(dc5$FROM_TO_ID))
	df5		<- tmp[, list(MCIT=seq_len(nrow(post$prop_comm_flow)), PI=post$prop_comm_flow[,FROM_TO_ID]), by='FROM_TO_ID']
	df5		<- merge(dc5, df5, by='FROM_TO_ID')
	df5		<- df5[, list(PI=sum(PI)), by=c('MCIT','TR_FISHCOMM','REC_FISHCOMM')]
	df5[, MODEL:='community-based prior Dir(0.8/NPAIR)']
		
	tmp		<- rbind(df1, df4, df5, fill=TRUE)
	ggplot(tmp, aes(x=paste0(TR_FISHCOMM,'-',REC_FISHCOMM), y=PI, fill=MODEL)) +
			geom_boxplot() +
			theme_bw()
	ggsave(file=paste0(outfile.base,'_stan_compare_m1_m4_m5.pdf'), w=7, h=7)
	
	#
	#	GP model
	#
	data(islandsDistMatrix)
	data(Kline2)
	d<- Kline2
	d$society<- 1:10
	m13.7 <- map2stan(
			alist(
					total_tools ~ dpois(lambda),
					log(lambda) <- a + g[society] + bp*logpop,
					g[society] ~ GPL2( Dmat, etasq, rhosq, 0.01),
					a ~ dnorm(0,10),
					bp ~ dnorm(0,1),
					etasq ~ dcauchy(0,1),
					rhosq ~ dcauchy(0,1)
					),
			data=list(total_tools=d$total_tools, logpop=d$logpop, society=d$society, Dmat=islandsDistMatrix),
			warmup=2e3, iter=3e3, chains=1
			)		
}

RakaiFull.phylogeography.180322.core.inference.stan.models.withseqsampling<- function()
{
	require(data.table)
	require(rstan)
	require(scales)
	require(ggplot2)
	require(RColorBrewer)
	
	infile.trms			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_withmetadata.rda"
	infile.sam			<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender.rda'
	outfile.base		<- gsub('_withmetadata.rda','',infile.trms)
	
	#	load previous model fits
	load( paste0(outfile.base,'_stanm1m4m5.rda') )
	
	#	load des which contains participation and seq counts by comm and gender	
	load(infile.sam)
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
	#	add new "community indices" 1...N by community and sex
	tmp	<- unique(subset(ds, select=c(COMM_NUM_A, SEX)))
	setkey(tmp, COMM_NUM_A, SEX)
	tmp[, COMM_NUM2:= seq_len(nrow(tmp))]	
	ds	<- merge(ds, tmp, by=c('COMM_NUM_A','SEX'))
	
	
	
	#	load transmission events
	load(infile.trms)
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	nrow(rtpdm)
	#	293
	
	#	add new variables
	rtpdm[, AGEDIFF:= rtpdm[, FEMALE_BIRTHDATE-MALE_BIRTHDATE]]	
	set(rtpdm, NULL, 'MALE_SEX', 'M')
	set(rtpdm, NULL, 'FEMALE_SEX', 'F')			
	rtpdm[, MALE_FISHCOMM:= rtpdm[, as.character(factor(MALE_COMM_TYPE=='fisherfolk', levels=c(TRUE,FALSE), labels=c('Fishing','Inland')))]]
	rtpdm[, FEMALE_FISHCOMM:= rtpdm[, as.character(factor(FEMALE_COMM_TYPE=='fisherfolk', levels=c(TRUE,FALSE), labels=c('Fishing','Inland')))]]
	
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
	
	#	count transmissions
	dc6	<- rtr2[, list(TR_OBS=length(PAIRID)), by=c('TR_FISHCOMM','REC_FISHCOMM','TR_COMM_NUM_A','REC_COMM_NUM_A','TR_SEX','REC_SEX')]
	#	add new indices and posterior alpha beta
	tmp	<- unique(subset(ds, select=c(COMM_NUM_A, SEX, COMM_NUM2)))	
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc6	<- merge(dc6, tmp, by=c('TR_COMM_NUM_A','TR_SEX'))
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc6	<- merge(dc6, tmp, by=c('REC_COMM_NUM_A','REC_SEX'))	
	dc6[, FROM_TO_ID:= seq_len(nrow(dc6))]	

	m6		<- list()
	m6$code	<- "data {
		int<lower=1> N; 					// number of observations
		int<lower=1> N_PAIRS;				// number of **observed** community pairs
		int TR_OBS[N]; 						// transmission counts	  	  
	}
	transformed data {
		real<lower=0> objective_prior_weight;
		vector<lower=0>[N_PAIRS] dirprior_comm_flow;
		objective_prior_weight= N_PAIRS;
		objective_prior_weight= 0.8/objective_prior_weight; //this is the Berger objective prior with minimal loss compared to marginal Beta reference prior	  
		dirprior_comm_flow= rep_vector(objective_prior_weight, N_PAIRS); // community flows, positive
	}	
	parameters {
		simplex[N_PAIRS] prop_comm_flow;  	// community flows, sum to 1	   	  
	}
	model {	  	  
		prop_comm_flow ~ dirichlet_log(dirprior_comm_flow);	//dirichlet prior to satisfy sum-to-1 constraint
		TR_OBS ~ multinomial(prop_comm_flow);	  
	}"			
	m6$code		<- gsub('\t','',m6$code)
	m6$codefile	<- paste0(outfile.base,'_stanm6.stan')
	#write(m1$code, file=m1$codefile)	
	m6$data					<- as.list(dc6)
	m6$data[['N']] 			<- nrow(dc6)
	m6$data[['N_PAIRS']] 	<- nrow(dc6)
	m6$fit		<- stan(model_code=m6$code, 
			data=m6$data, 
			warmup=5e2, iter=1e4, chains=1,
			init=list(list(base_flow=exp(1), coeff_space=0, rhosq=1, sig_random_pair_effect=1, random_pair_effect=rep(0,m6$data[['N_PAIRS']]))))				
	precis(m6$fit, depth=2, prob=0.95)
	post <- as.data.frame(m6$fit)
	pdf(file=paste0(outfile.base,'_stanm6_traces.pdf'), w=10, h=100)
	traceplot(m6$fit, pars=colnames(post), ncol=1)
	dev.off()
	
	#	setup transmissions
	dc7	<- rtr2[, list(TR_OBS=length(PAIRID)), by=c('TR_FISHCOMM','REC_FISHCOMM','TR_COMM_NUM_A','REC_COMM_NUM_A','TR_SEX','REC_SEX')]
	#	add new indices and posterior alpha beta
	tmp	<- unique(subset(ds, select=c(COMM_NUM_A, COMM_NUM2, SEX, P_PART_EMP, P_SEQ_EMP)))	
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc7	<- merge(dc7, tmp, by=c('TR_COMM_NUM_A','TR_SEX'))
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc7	<- merge(dc7, tmp, by=c('REC_COMM_NUM_A','REC_SEX'))	
	dc7[, FROM_TO_ID:= seq_len(nrow(dc7))]
	#	work out maximum numer missing
	dc7[, P:= TR_P_PART_EMP*TR_P_SEQ_EMP*REC_P_PART_EMP*REC_P_SEQ_EMP]
	dc7[, MAXMISS:= qnbinom(0.05, size=TR_OBS, prob=P, lower.tail=FALSE)]
	
	dc7[, prod(MAXMISS)]
	#	setup sampling
	ds7	<- subset(ds, select=c(COMM_NUM2, P_PART_ALPHA, P_PART_BETA, P_SEQ_ALPHA, P_SEQ_BETA))
	       
	
	m7		<- list()
	m7$code	<- "data {
		// transmissions
		int<lower=1> N; 					// number of observations		
		int<lower=1> N_PAIRS;				// number of **observed** community pairs
		int<lower=1> MXMISS;				// cut off on the maximum number of missing transmissions	
		int<lower=1> TR_COMM_NUM2[N];		// index of transmitting community
		int<lower=1> REC_COMM_NUM2[N];		// index of recipient community
		int TR_OBS[N]; 						// transmission counts
		// sampling
		int<lower=1> N_COM;					// number of communities * 2 (for male female)
		int<lower=1> COMM_NUM2[N_COM];		// index of community-gender pair
		real<lower=0> P_PART_ALPHA[N_COM];	// alpha param of posterior Beta distribution for participation in comm x and gender y
		real<lower=0> P_PART_BETA[N_COM];	// beta param of posterior Beta distribution for participation in comm x and gender y
		real<lower=0> P_SEQ_ALPHA[N_COM];	// alpha param of posterior Beta distribution for seq sampling in comm x and gender y
		real<lower=0> P_SEQ_BETA[N_COM];	// beta param of posterior Beta distribution for seq sampling in comm x and gender y						  	  
	}
	transformed data {
		real<lower=0> objective_prior_weight;		
		vector<lower=0>[N_PAIRS] dirprior_comm_flow;
		objective_prior_weight= N_PAIRS;
		objective_prior_weight= 0.8/objective_prior_weight; //this is the Berger objective prior with minimal loss compared to marginal Beta reference prior			  
		dirprior_comm_flow= rep_vector(objective_prior_weight, N_PAIRS); // community flows, positive
	}	
	parameters {
		simplex[N_PAIRS] prop_comm_flow;  			// community flows, sum to 1
		vector<lower=0, upper=1>[N_COM] prob_part;	// participation probability in comm x and gender y 
		vector<lower=0, upper=1>[N_COM] prob_seq;   // sequencing probability in comm x and gender y
		int TR_MISS[N]; 							// missed transmission counts under NB model 
	}
	transformed parameters {		
		vector<lower=0, upper=1>[N] prob_obs;		//  prob_part * prob_seq (from comm) * prob_part * prob_seq (to comm)
		vector<lower=0, upper=1>[N] nb2_mu;			//  transformation to NegBin2 mu parameter
		prob_obs= prob_part[TR_COMM_NUM2] .* prob_seq[TR_COMM_NUM2] .* prob_part[REC_COMM_NUM2] .* prob_seq[REC_COMM_NUM2];
		nb2_mu= TR_OBS .* (1-prob_obs) ./ prob_obs;
	}
	model {	  	
		prob_part ~ beta(P_PART_ALPHA, P_PART_BETA); 		// Empirical Bayes prior
		prob_seq ~ beta(P_SEQ_ALPHA, P_SEQ_BETA); 	 		// Empirical Bayes prior
		prop_comm_flow ~ dirichlet_log(dirprior_comm_flow);	//dirichlet prior to satisfy sum-to-1 constraint
		TR_MISS ~ neg_binomial_2_rng(nb2_mu, TR_OBS);
		TR_OBS ~ multinomial(prop_comm_flow);	  
	}"			
	m7$code		<- gsub('\t','',m7$code)
	m7$codefile	<- paste0(outfile.base,'_stanm7.stan')
	#write(m1$code, file=m1$codefile)	
	m7$data						<- as.list(ds7)
	m7$data[['N']] 				<- nrow(dc7)
	m7$data[['N_COM']] 			<- nrow(ds7)
	m7$data[['N_PAIRS']] 		<- nrow(dc7)
	m7$data[['MXMISS']] 		<- dc7[, max(MAXMISS)]
	m7$data[['TR_COMM_NUM2']]	<- dc7[,TR_COMM_NUM2]
	m7$data[['REC_COMM_NUM2']]	<- dc7[,REC_COMM_NUM2]
	m7$data[['TR_OBS']]			<- dc7[,TR_OBS]
	m7$fit		<- stan(model_code=m7$code, 
			data=m7$data, 
			warmup=5e2, iter=1e3, chains=1,
			init=list(list(prop_comm_flow=rep(1/nrow(dc7),nrow(dc7)))))				
	precis(m7$fit, depth=2, prob=0.95)
	post <- as.data.frame(m7$fit)
	pdf(file=paste0(outfile.base,'_stanm7_traces.pdf'), w=10, h=100)
	traceplot(m7$fit, pars=colnames(post), ncol=1)
	dev.off()
	
	
	m8$code	<- "
		data {
			// counts
			int<lower=1> N; 					// number of observations									
			int COUNTS[N]; 						// observed counts
			// sampling
			real<lower=0> P_ALPHA[N];			// sampling probabilities prior hyperparameter alpha 
			real<lower=0> P_BETA[N];			// sampling probabilities prior hyperparameter beta									  	  
		}
		transformed data {
			// dirichlet hyperparameters depending on size of data
			real<lower=0> prior_weight;		
			vector<lower=0>[N] dirprior;
			prior_weight= N;
			prior_weight= 0.8/prior_weight; 	// poor man casting	  
			dirprior= rep_vector(prior_weight, N); 
		}	
		parameters {
			simplex[N_PAIRS] props;  					// target parameter
			vector<lower=0, upper=1>[N] prob_obs;		// sampling parameter
			int MISSING[N];								// latent variable, to be integrated out	 
		}
		transformed parameters {		
			vector<lower=0, upper=1>[N] nb2_mu;			//  transformation to NegBin2 mu parameter			
			vector<lower=0, upper=1>[N] nb2_phi;		//  transformation to NegBin2 phi parameter
			int AUX_COUNTS[N];							
			nb2_mu= COUNTS .* (1-prob_obs) ./ prob_obs;
			nb2_phi= COUNTS;
			AUX_COUNTS= COUNTS + MISSING;
		}
		model {	  	
			prob_obs ~ beta(P_ALPHA, P_BETA); 			// prior
			props ~ dirichlet_log(dirprior);			// Dirichlet prior 
			MISSING ~ neg_binomial_2_rng(nb2_mu, nb2_phi);
			AUX_COUNTS ~ multinomial(props);	  		// augmented likelihood
		}"	
	
	#	save model fits
	save(dc1, dc5, m1, m4, m5, file=paste0(outfile.base,'_stanm1m4m5.rda'))
	
	
	
	
	#	question: how different are estimates from m2 m4 m5 m6?
	post	<- extract.samples(m1$fit)
	tmp		<- data.table(FROM_TO_ID=1:max(dc1$FROM_TO_ID))
	df1		<- tmp[, list(MCIT=seq_len(nrow(post$prop_comm_flow)), PI=post$prop_comm_flow[,FROM_TO_ID]), by='FROM_TO_ID']
	df1		<- merge(dc1, df1, by='FROM_TO_ID')	
	df1[, MODEL:='area-based']
	
	post	<- extract.samples(m4$fit)
	tmp		<- data.table(FROM_TO_ID=1:max(dc5$FROM_TO_ID))
	df4		<- tmp[, list(MCIT=seq_len(nrow(post$prop_comm_flow)), PI=post$prop_comm_flow[,FROM_TO_ID]), by='FROM_TO_ID']
	df4		<- merge(dc5, df4, by='FROM_TO_ID')
	df4		<- df4[, list(PI=sum(PI)), by=c('MCIT','TR_FISHCOMM','REC_FISHCOMM')]
	df4[, MODEL:='community-based prior Dir(1)']
	
	post	<- extract.samples(m5$fit)
	tmp		<- data.table(FROM_TO_ID=1:max(dc5$FROM_TO_ID))
	df5		<- tmp[, list(MCIT=seq_len(nrow(post$prop_comm_flow)), PI=post$prop_comm_flow[,FROM_TO_ID]), by='FROM_TO_ID']
	df5		<- merge(dc5, df5, by='FROM_TO_ID')
	df5		<- df5[, list(PI=sum(PI)), by=c('MCIT','TR_FISHCOMM','REC_FISHCOMM')]
	df5[, MODEL:='community-based prior Dir(0.8/NPAIR)']
	
	post	<- extract.samples(m6$fit)
	tmp		<- data.table(FROM_TO_ID=1:nrow(dc6))
	df6		<- tmp[, list(MCIT=seq_len(nrow(post$prop_comm_flow)), PI=post$prop_comm_flow[,FROM_TO_ID]), by='FROM_TO_ID']
	df6		<- merge(dc6, df6, by='FROM_TO_ID')
	df6		<- df6[, list(PI=sum(PI)), by=c('MCIT','TR_FISHCOMM','REC_FISHCOMM')]
	df6[, MODEL:='community-based prior Dir(0.8/NPAIR) M and F sep']
	
	
	tmp		<- rbind(df1, df4, df5, df6, fill=TRUE)
	ggplot(tmp, aes(x=paste0(TR_FISHCOMM,'-',REC_FISHCOMM), y=PI, fill=MODEL)) +
			geom_boxplot() +
			theme_bw()
	ggsave(file=paste0(outfile.base,'_stan_compare_m1_m4_m5.pdf'), w=7, h=7)	
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

RakaiFull.phylogeography.180521.gender.mobility.data<- function(infile.inference)
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
	
	opt.set.missing.migloc.to.inland	<- 0
	opt.set.missing.migloc.to.fishing	<- 1
	
	hivc.db.Date2numeric<- function( x )
	{
		if(!class(x)%in%c('Date','character'))	return( x )
		x	<- as.POSIXlt(x)
		tmp	<- x$year + 1900
		x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
		x	
	}	
	
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
	ds[, P_SEQ_BETA:= HIV_1516_YES-DEEP_SEQ_1516+1]
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
	
	
	#	prepare inmigrant -- identify inmigrants from fishing communities and from external
	infile		<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData_v2.rda"
	load(infile)
	inmigrant	<- as.data.table(inmigrant)
	infile		<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/migrants_withMissingGPS.csv"
	inmigrant2	<- as.data.table(read.csv(infile))
	
	
	#	plot fisherfolk	to figure out how much of a radius we need
	if(0)
	{
		zf		<- data.table(longitude=c(31.763,31.7968,31.754,31.838), latitude=c(-0.915, -0.6518, -0.703, -0.497), ID= c('Kasensero','Bukyanju','NearBwende','Fish4'))
		make_circles <- function(centers, radius, nPoints = 100){
			# centers: the data frame of centers with ID
			# radius: radius measured in kilometer
			#
			meanLat <- mean(centers$latitude)
			# length per longitude changes with lattitude, so need correction
			radiusLon <- radius /111 / cos(meanLat/57.3) 
			radiusLat <- radius / 111
			circleDF <- data.frame(ID = rep(centers$ID, each = nPoints))
			angle <- seq(0,2*pi,length.out = nPoints)
			
			circleDF$lon <- unlist(lapply(centers$longitude, function(x) x + radius * cos(angle)))
			circleDF$lat <- unlist(lapply(centers$latitude, function(x) x + radius * sin(angle)))
			return(circleDF)
		}
		zc <- make_circles(zf, 0.01)
		ggmap(zm) +			
				geom_point(data=zf, aes(x=longitude, y=latitude, pch=ID), stroke=1.5, alpha=0.8) +
				geom_polygon(data=zc, aes(lon, lat, group = ID), color = "red", alpha = 0)
		#	radius of length 0.01 should catch					
		tmp		<- inmigrant[, list( 	DIST_KASENSERO= sqrt( (inmig_lon- 31.763)^2 + (inmig_lat - (-0.915))^2),
										DIST_BUKYANJU= sqrt( (inmig_lon- 31.7968)^2 + (inmig_lat - (-0.6518))^2),
										DIST_NEARBWENDE= sqrt( (inmig_lon- 31.754)^2 + (inmig_lat - (-0.703))^2),
										DIST_FISH4= sqrt( (inmig_lon- 31.838)^2 + (inmig_lat - (-0.497))^2)
									), by=c('RCCS_studyid','visit')]
		tmp		<- melt(tmp, id.vars=c('RCCS_studyid','visit'))
		ggplot(subset(tmp, value<0.3), aes(x=value)) +
				geom_histogram(binwidth=0.01) +
				facet_grid(variable~.)		
		zfd		<- merge(inmigrant, subset(tmp, value<0.01, c(RCCS_studyid, visit)), by=c('RCCS_studyid','visit'))		
	}
	#	so fishing sites are MALEMBO DIMU KASENSERO NAMIREMBE but there are spelling mistakes
	#	clean up inmigrant	
	#
	#	inmigrant[, unique(sort(inmig_place))]
	set(inmigrant, NULL, 'inmig_place', inmigrant[, gsub('DDIMO|DDIMU|DIMO|DIMU','DIMU',inmig_place)])
	set(inmigrant, inmigrant[, which(grepl('MALEMBO',inmig_place))], 'inmig_place', 'MALEMBO')
	set(inmigrant, NULL, 'inmig_place', inmigrant[, gsub("KASEMSERO","KASENSERO",inmig_place)])
	set(inmigrant, inmigrant[, which(grepl('KASENSERO',inmig_place))], 'inmig_place', 'KASENSERO')
	set(inmigrant2, inmigrant2[, which(grepl('MALEMBO',inmig_place))], 'inmig_place', 'MALEMBO')
	
	#	define from_fishing and from_outside and from_inland
	inmigrant[, INMIG_LOC:= 'inland' ]
	set(inmigrant, inmigrant[, which(grepl('MALEMBO|DIMU|KASENSERO|NAMIREMBE',inmig_place))], 'INMIG_LOC','fisherfolk')
	set(inmigrant, inmigrant[, which(inmig_admin0!='Uganda')], 'INMIG_LOC','external')
	set(inmigrant, inmigrant[, which(inmig_admin1!='Rakai')], 'INMIG_LOC','external')
	set(inmigrant, inmigrant[, which(is.na(inmig_admin1))], 'INMIG_LOC','unknown')	
	inmigrant2[, INMIG_LOC:= 'inland' ]
	set(inmigrant2, inmigrant2[, which(grepl('MALEMBO|DIMU|KASENSERO|NAMIREMBE',inmig_place))], 'INMIG_LOC','fisherfolk')
	set(inmigrant2, inmigrant2[, which(inmig_admin0!='Uganda')], 'INMIG_LOC','external')
	set(inmigrant2, inmigrant2[, which(inmig_admin1!='Rakai')], 'INMIG_LOC','external')
	set(inmigrant2, inmigrant2[, which(is.na(inmig_admin1))], 'INMIG_LOC','unknown')
	inmigrant2[, table(INMIG_LOC)]
	if(opt.set.missing.migloc.to.inland)
	{
		set(inmigrant2, inmigrant2[, which(INMIG_LOC=='unknown')], 'INMIG_LOC', 'inland')
	}
	if(opt.set.missing.migloc.to.fishing)
	{
		set(inmigrant2, inmigrant2[, which(INMIG_LOC=='unknown')], 'INMIG_LOC', 'fisherfolk')
	}
	setnames(inmigrant2, 'INMIG_LOC', 'INMIG_LOC2')
	inmigrant2	<- subset(inmigrant2, select=c('RCCS_studyid','visit','INMIG_LOC2'))
	inmigrant	<- merge(inmigrant, inmigrant2, by=c('RCCS_studyid','visit'), all.x=TRUE)
	tmp			<- inmigrant[, which(!is.na(INMIG_LOC2))]
	set(inmigrant, tmp, 'INMIG_LOC', inmigrant[tmp, INMIG_LOC2])
	
	inmigrant[, table(INMIG_LOC)]
	#	  external fisherfolk     inland    unknown 
    #	   763         45       1577        369 
	inmigrant	<- subset(inmigrant, select=c(RCCS_studyid, date, INMIG_LOC))	
	inmigrant[, date:= hivc.db.Date2numeric(date)]	
	setnames(inmigrant, colnames(inmigrant), gsub('DATE','INMIGRATE_DATE',gsub('RCCS_STUDYID','RID',toupper(colnames(inmigrant)))))


	#	load transmission events
	if(is.null(infile.inference))
	{
		infile.inference					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_withmetadata.rda"
		infile.inference					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_withmetadata.rda"
		infile.inference					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min10_withmetadata.rda"
		infile.inference					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min20_withmetadata.rda"
		infile.inference					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min50_withmetadata.rda"		
	}
	outfile.base			<- gsub('171122','180522',gsub('_withmetadata.rda','',infile.inference))
	load(infile.inference)
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	#stopifnot(293==nrow(rtpdm))	
	#	add new variables
	rtpdm[, AGEDIFF:= rtpdm[, FEMALE_BIRTHDATE-MALE_BIRTHDATE]]	
	set(rtpdm, NULL, 'MALE_SEX', 'M')
	set(rtpdm, NULL, 'FEMALE_SEX', 'F')		
	#	add immigrant: to compute if individual was inmigrant within 1 year of first conc pos, we need date first conc pos	
	tmp			<- subset(rtpdm, select=c(PAIRID, MALE_RID, FEMALE_RID, VISIT_FIRSTCONCPOS))
	tmp2		<- unique(subset(rd, select=c(RID, VISIT, DATE)))
	setnames(tmp2, colnames(tmp2), gsub('MALE_VISIT','VISIT_FIRSTCONCPOS',paste0('MALE_',colnames(tmp2))))
	tmp			<- merge(tmp, tmp2, by=c('MALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	setnames(tmp2, colnames(tmp2), gsub('MALE_','FEMALE_', colnames(tmp2)))
	tmp			<- merge(tmp, tmp2, by=c('FEMALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	set(tmp, tmp[, which(is.na(MALE_DATE))], 'MALE_DATE', -1)
	set(tmp, tmp[, which(is.na(FEMALE_DATE))], 'FEMALE_DATE', -1)
	#	argh there are entries where both have missing date in rd ..
	tmp2		<- subset(tmp, MALE_DATE==-1 & FEMALE_DATE==-1, select=c(MALE_RID, FEMALE_RID, VISIT_FIRSTCONCPOS, PAIRID))
	tmp2		<- merge(tmp2, subset(rtpdm, select=c(MALE_RID, FEMALE_RID, VISIT_FIRSTCONCPOS, PAIRID, MALE_DATE, FEMALE_DATE)), by=c('MALE_RID','FEMALE_RID','VISIT_FIRSTCONCPOS','PAIRID'))
	tmp			<- rbind(subset(tmp, !(MALE_DATE==-1 & FEMALE_DATE==-1)), tmp2)	
	set(tmp, NULL, 'DATE_FIRSTCONCPOS', tmp[, pmax(MALE_DATE, FEMALE_DATE)])
	stopifnot( nrow(subset(tmp, DATE_FIRSTCONCPOS==-1))==0 )
	set(tmp, NULL, c('VISIT_FIRSTCONCPOS','MALE_DATE','FEMALE_DATE'), NULL)
	#	find inmigrants among transmission pairs
	tmp2		<- copy(inmigrant)
	setnames(tmp2, colnames(tmp2), paste0('MALE_',colnames(tmp2)))	
	tmp			<- merge(tmp, tmp2, by='MALE_RID', all.x=1)
	setnames(tmp2, colnames(tmp2), gsub('MALE_','FEMALE_', colnames(tmp2)))
	tmp			<- merge(tmp, tmp2, by='FEMALE_RID', all.x=1)
	set(tmp, tmp[, which(is.na(MALE_INMIGRATE_DATE))], 'MALE_INMIGRATE_DATE', -1)
	set(tmp, tmp[, which(is.na(FEMALE_INMIGRATE_DATE))], 'FEMALE_INMIGRATE_DATE', -1)
	#	don t consider inmigrations after transmission occurred
	set(tmp, tmp[, which(DATE_FIRSTCONCPOS<MALE_INMIGRATE_DATE)], c('MALE_INMIG_PLACE','MALE_INMIG_ADMIN0','MALE_INMIG_ADMIN1','MALE_INMIG_ADMIN2'), NA_character_)
	set(tmp, tmp[, which(DATE_FIRSTCONCPOS<MALE_INMIGRATE_DATE)], 'MALE_INMIGRATE_DATE', -1)
	set(tmp, tmp[, which(DATE_FIRSTCONCPOS<FEMALE_INMIGRATE_DATE)], c('FEMALE_INMIG_PLACE','FEMALE_INMIG_ADMIN0','FEMALE_INMIG_ADMIN1','FEMALE_INMIG_ADMIN2'), NA_character_)
	set(tmp, tmp[, which(DATE_FIRSTCONCPOS<FEMALE_INMIGRATE_DATE)], 'FEMALE_INMIGRATE_DATE', -1)
	tmp[, MALE_TIME_INMIGRATED:= DATE_FIRSTCONCPOS-MALE_INMIGRATE_DATE]
	tmp[, FEMALE_TIME_INMIGRATED:= DATE_FIRSTCONCPOS-FEMALE_INMIGRATE_DATE]	
	#	select most recent inmigration events
	tmp[, DUMMY:= seq_len(nrow(tmp))]
	tmp2	<- tmp[, list(DUMMY=DUMMY[min(FEMALE_TIME_INMIGRATED)==FEMALE_TIME_INMIGRATED]), by=c('PAIRID')]
	tmp		<- merge(tmp, tmp2, by=c('PAIRID','DUMMY'))
	tmp2	<- tmp[, list(DUMMY=DUMMY[min(MALE_TIME_INMIGRATED)==MALE_TIME_INMIGRATED]), by=c('PAIRID')]
	tmp		<- merge(tmp, tmp2, by=c('PAIRID','DUMMY'))
	tmp[, DUMMY:=NULL]
	#	define inmigrated in last year, and in last two years
	tmp[, MALE_INMIGRATE_1YR:= as.character(as.integer(MALE_TIME_INMIGRATED<=1))]
	tmp[, MALE_INMIGRATE_2YR:= as.character(as.integer(MALE_TIME_INMIGRATED<=2))]
	tmp[, MALE_INMIGRATE_3YR:= as.character(as.integer(MALE_TIME_INMIGRATED<=3))]
	tmp[, FEMALE_INMIGRATE_1YR:= as.character(as.integer(FEMALE_TIME_INMIGRATED<=1))]
	tmp[, FEMALE_INMIGRATE_2YR:= as.character(as.integer(FEMALE_TIME_INMIGRATED<=2))]	
	tmp[, FEMALE_INMIGRATE_3YR:= as.character(as.integer(FEMALE_TIME_INMIGRATED<=3))]
	tmp		<- subset(tmp, select=c(PAIRID, DATE_FIRSTCONCPOS, MALE_INMIG_LOC, MALE_INMIGRATE_1YR, MALE_INMIGRATE_2YR, MALE_INMIGRATE_3YR, FEMALE_INMIG_LOC, FEMALE_INMIGRATE_1YR, FEMALE_INMIGRATE_2YR, FEMALE_INMIGRATE_3YR))
	#	make compound variable
	tmp		<- melt(tmp, id.vars=c('PAIRID','DATE_FIRSTCONCPOS','MALE_INMIG_LOC','FEMALE_INMIG_LOC'))
	set(tmp, tmp[, which(value=='0')], 'value', 'resident')
	tmp2	<- tmp[, which(value=='1' & grepl('^FEMALE_',variable))]
	set(tmp, tmp2, 'value', tmp[tmp2, paste0('inmigrant_from_',FEMALE_INMIG_LOC)])
	tmp2	<- tmp[, which(value=='1' & grepl('^MALE_',variable))]
	set(tmp, tmp2, 'value', tmp[tmp2, paste0('inmigrant_from_',MALE_INMIG_LOC)])
	tmp		<- dcast.data.table(tmp, PAIRID+DATE_FIRSTCONCPOS~variable, value.var='value')
	
	#	tmp[, table(MALE_INMIGRATE_2YR,  MALE_INMIGRATE_3YR)]
	#	only one more migration event when we allow for 3 years
	#	ignore this
	set(tmp, NULL, c('FEMALE_INMIGRATE_3YR','MALE_INMIGRATE_3YR'), NULL)

	rtpdm		<- merge(rtpdm, tmp, by=c('PAIRID'))
	rtpdm[, table(MALE_INMIGRATE_2YR, FEMALE_INMIGRATE_2YR)]
	#                           	FEMALE_INMIGRATE_2YR
	#	MALE_INMIGRATE_2YR          inmigrant_from_external inmigrant_from_fisherfolk inmigrant_from_inland resident
 	# inmigrant_from_external                         5                         0                     5        8
  	# inmigrant_from_fisherfolk                       0                         0                     0        1
  	# inmigrant_from_inland                           3                         0                    10       22
  	# resident                                       22                         6                    37      174
				 
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
	
	#	save
	save(rtpdm, rtr2, zm, ds, file=paste0(outfile.base,'_phylogeography_data_with_inmigrants.rda'))
}


RakaiFull.phylogeography.180521.flows.fishinlandmigrationgender.netflows<- function()
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
	
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_core_inference.rda"
	outfile.base		<- gsub('_inference.rda','',infile.inference)
	load(infile.inference)
	
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')
	
	
	#
	#	geography flows inland->fish / flows fish->inland
	#	start by gender	
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))		
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset to simple
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	set(z, z[, which(TR_INMIGRATE=='resident')], 'TR_INMIGRATE', 'resident/outmigrant')	
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_SEX','REC_SEX','TR_INMIGRATE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',TR_INMIGRATE,' ',TR_SEX,' ',REC_COMM_TYPE,' ',REC_SEX)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	set(z, NULL, c('inland resident/outmigrant F inland M','external inmigrant_from_external F inland M','external inmigrant_from_external F fisherfolk M','fisherfolk resident/outmigrant F fisherfolk M'), NULL)
	set(z, NULL, c('inland resident/outmigrant M inland F','external inmigrant_from_external M inland F','external inmigrant_from_external M fisherfolk F','fisherfolk resident/outmigrant M fisherfolk F'), NULL)
	z[ , inlanddivfisherfolk_M:= z[['inland resident/outmigrant M fisherfolk F']] / z[['fisherfolk resident/outmigrant M inland F']] ]
	z[ , inlanddivfisherfolk_F:= z[['inland resident/outmigrant F fisherfolk M']] / z[['fisherfolk resident/outmigrant F inland M']] ]				
	z		<- melt(z, id.vars='MONTE_CARLO_IT', measure.vars=c('inlanddivfisherfolk_M','inlanddivfisherfolk_F'))
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, SINK_TYPE:= gsub('([^_]+)_([^_]+)','\\1',variable)]
	z[, SINK_SEX:= gsub('([^_]+)_([^_]+)','\\2',variable)]						
	z		<- dcast.data.table(z, SINK_TYPE+SINK_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),' - ',round(CU,d=2),']')]
	z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'-',round(CU,d=2),')')]
	setkey(z, SINK_TYPE, SINK_SEX )
	ans		<- copy(z)
	
	
	
	#
	#	geography flows inland->fish / flows fish->inland
	#	add by migration status	
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))		
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')	
	#	reset to complex migrant
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_INMIGRATE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',TR_INMIGRATE,' ',REC_COMM_TYPE)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	set(z, NULL, c('inland resident inland','inland outmigrant inland', 'external inmigrant_from_external inland','external inmigrant_from_external fisherfolk','fisherfolk resident fisherfolk'), NULL)	
	z[ , inlanddivfisherfolk_resident:= z[['inland resident fisherfolk']] / z[['fisherfolk resident inland']] ]
	z[ , inlanddivfisherfolk_outmigrant:= z[['inland outmigrant fisherfolk']] / z[['fisherfolk outmigrant inland']] ]				
	z		<- melt(z, id.vars='MONTE_CARLO_IT', measure.vars=c('inlanddivfisherfolk_resident','inlanddivfisherfolk_outmigrant'))
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, SINK_TYPE:= gsub('([^_]+)_([^_]+)','\\1',variable)]
	z[, SINK_MIGRATIONSTATUS:= gsub('([^_]+)_([^_]+)','\\2',variable)]						
	z		<- dcast.data.table(z, SINK_TYPE+SINK_MIGRATIONSTATUS~P, value.var='Q')
	z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),' - ',round(CU,d=2),']')]
	z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'-',round(CU,d=2),')')]
	z[, SINK_SEX:= 'Any']	
	setkey(z, SINK_TYPE, SINK_MIGRATIONSTATUS )	
	ans		<- rbind(z, ans, fill=TRUE)
	
	#
	#	add overall
	#
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))		
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset to simple
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	set(z, z[, which(TR_INMIGRATE=='resident')], 'TR_INMIGRATE', 'resident/outmigrant')	
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_INMIGRATE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',TR_INMIGRATE,' ',REC_COMM_TYPE)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	set(z, NULL, c('inland resident/outmigrant inland','external inmigrant_from_external inland','external inmigrant_from_external fisherfolk','fisherfolk resident/outmigrant fisherfolk'), NULL)
	z[ , inlanddivfisherfolk:= z[['inland resident/outmigrant fisherfolk']] / z[['fisherfolk resident/outmigrant inland']] ]
	z		<- melt(z, id.vars='MONTE_CARLO_IT', measure.vars=c('inlanddivfisherfolk'))
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, SINK_TYPE:= gsub('([^_]+)_([^_]+)','\\1',variable)]
	z[, SINK_SEX:= 'Any']						
	z		<- dcast.data.table(z, SINK_TYPE+SINK_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),' - ',round(CU,d=2),']')]
	z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'-',round(CU,d=2),')')]
	setkey(z, SINK_TYPE, SINK_SEX )
	ans		<- rbind(z, ans)
	
	#	save
	save(ans, file=paste0(outfile.base,'_flows_sinkfishinlandgender.rda'))
	
	
	#
	ggplot(ans, aes(x=SINK_SEX, middle=M, min=CL, lower=IL, upper=IU, max=CU)) +
			geom_hline(yintercept=1, colour='grey50', size=1.5) +
			geom_boxplot(stat='identity', fill='grey50') +
			theme_bw() + theme(legend.position='bottom') +
			scale_y_log10(expand=c(0,0), breaks=c(1/4,1/3,1/2,2/3,1,3/2,2,3,4,8,16,32), labels=c('1/4','1/3','1/2','2/3','1','3/2','2','3','4','8','16','32')) +
			coord_flip(ylim=c(2/3,40)) +			
			labs(x='Gender\n', y='\nflows inland->fishing / flows fishing->inland') +
			guides(fill='none')
	ggsave(file=paste0(outfile.base,'_flows_sinkfishinlandgender.pdf'), w=6, h=3)	
}

RakaiFull.phylogeography.180521.flows.fishinlandgender<- function()
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
	
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_core_inference.rda"
	outfile.base		<- gsub('_inference.rda','',infile.inference)
	load(infile.inference)
	
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')
	
	#
	#	geography who infects whom matrix  between fisherfolk and others account for migration (complex)
	#	adjusted P
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
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
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',TR_SEX,' ',REC_COMM_TYPE, ' ', REC_SEX)]
	#unique(z$FLOW)
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+) ([^ ]+) ([^ ]+)','\\1',variable)]
	z[, TR_SEX:= gsub('([^ ]+) ([^ ]+) ([^ ]+) ([^ ]+)','\\2',variable)]
	z[, REC_COMM_TYPE:= gsub('([^ ]+) ([^ ]+) ([^ ]+) ([^ ]+)','\\3',variable)]
	z[, REC_SEX:= gsub('([^ ]+) ([^ ]+) ([^ ]+) ([^ ]+)','\\4',variable)]
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_SEX+REC_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_SEX, REC_SEX )
	z[, STAT:='joint_sex']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- copy(z)	
	
	#	make plot
	tmp		<- copy(z)	
	set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','External',gsub('inland','Inland',gsub('from_fisherfolk','Fishing',TR_COMM_TYPE)))])
	set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','Inland',gsub('fisherfolk','Fishing',REC_COMM_TYPE))])
	set(tmp, NULL, 'TR_SEX', tmp[, gsub('M','male',TR_SEX)])
	set(tmp, NULL, 'TR_SEX', tmp[, gsub('F','female',TR_SEX)])
	set(tmp, NULL, 'REC_SEX', tmp[, gsub('M','male',REC_SEX)])
	set(tmp, NULL, 'REC_SEX', tmp[, gsub('F','female',REC_SEX)])
	tmp[, MODEL:= 'MonteCarloAdjustSampling']	
	tmp[, X:= paste0(TR_SEX, ' ', TR_COMM_TYPE,'->',REC_SEX, ' ',REC_COMM_TYPE)]
	ggplot(tmp, aes(x=X, fill=MODEL)) +
			geom_boxplot(aes(ymin=CL, lower=IL, upper=IU, ymax=CU, middle=M), stat='identity') +
			theme_bw() + 
			coord_flip() +
			scale_y_continuous(label=scales:::percent) +
			labs(x='Transmissions from -> to\n', y='\nProportion of transmissions')
	ggsave(file=paste0(outfile.base,'_flows_fishinlandgender.pdf'), w=8, h=5)
	
	
	#
	#	WAIFM matrix (simple model)
	#
	groups	<- data.table(TR_COMM_TYPE=c('inland','inland','fisherfolk','fisherfolk'), TR_SEX=c('M','F','M','F'))
	z		<- lapply(1:nrow(groups), function(ii)
			{				
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
				#set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				#
				z		<- subset(z, TR_COMM_TYPE==groups$TR_COMM_TYPE[ii] & TR_SEX==groups$TR_SEX[ii])				
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('REC_COMM_TYPE','REC_SEX','MONTE_CARLO_IT')]
				z[, FLOW:=paste0(REC_COMM_TYPE, ' ', REC_SEX)]
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- FLOW
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	get quantiles
				z		<- melt(z, id.vars=c('MONTE_CARLO_IT'))								
				z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by='variable']
				z[, REC_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\1',variable)]
				z[, REC_SEX:= gsub('([^ ]+) ([^ ]+)','\\2',variable)]				
				z[, TR_COMM_TYPE:= groups$TR_COMM_TYPE[ii] ]
				z[, TR_SEX:= groups$TR_SEX[ii] ] 
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_SEX+REC_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_SEX, REC_SEX)
	z[, STAT:='waifm_sex']	
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	
	#
	#	sources (simple model)
	#	
	groups	<- data.table(REC_COMM_TYPE=c('inland','inland','fisherfolk','fisherfolk'), REC_SEX=c('M','F','M','F'))
	z		<- lapply(1:nrow(groups), function(ii)
			{		
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
				#set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				#
				z		<- subset(z, REC_COMM_TYPE==groups$REC_COMM_TYPE[ii] & REC_SEX==groups$REC_SEX[ii])				
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('TR_COMM_TYPE','TR_SEX','MONTE_CARLO_IT')]
				z[, FLOW:=paste0(TR_COMM_TYPE, ' ', TR_SEX)]
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- FLOW
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	get quantiles
				z		<- melt(z, id.vars=c('MONTE_CARLO_IT'))								
				z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by='variable']
				z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\1',variable)]
				z[, TR_SEX:= gsub('([^ ]+) ([^ ]+)','\\2',variable)]				
				z[, REC_COMM_TYPE:= groups$REC_COMM_TYPE[ii] ]
				z[, REC_SEX:= groups$REC_SEX[ii] ] 							 
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_SEX+REC_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_SEX, REC_SEX)
	z[, STAT:='sources_sex']	
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans, z)		
	
	
	#	write table
	df	<- copy(ans)
	df	<- dcast.data.table(df, TR_COMM_TYPE+TR_SEX+REC_COMM_TYPE+REC_SEX~STAT, value.var='LABEL2')
	setkey(df, TR_COMM_TYPE, TR_SEX, REC_COMM_TYPE, REC_SEX)
	write.csv(df, row.names=FALSE, file=paste0(outfile.base, '_table_fishinlandgender.csv'))
	
	#	save
	save(ans, file=paste0(outfile.base,'_flows_fishinlandgender.rda'))
	
	
	#	plot
	tmp		<- copy(z)	
	set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','External',gsub('inland','Inland',gsub('fisherfolk','Fishing',TR_COMM_TYPE)))])
	set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','Inland',gsub('fisherfolk','Fishing',REC_COMM_TYPE))])
	set(tmp, NULL, 'TR_SEX', tmp[, gsub('M','male',TR_SEX)])
	set(tmp, NULL, 'TR_SEX', tmp[, gsub('F','female',TR_SEX)])
	set(tmp, NULL, 'REC_SEX', tmp[, gsub('M','male',REC_SEX)])
	set(tmp, NULL, 'REC_SEX', tmp[, gsub('F','female',REC_SEX)])	
	tmp[, MODEL:= 'MonteCarloAdjustSampling']	
	tmp[, X:= paste0(TR_SEX, ' ', TR_COMM_TYPE,'->',REC_SEX, ' ',REC_COMM_TYPE)]
	ggplot(tmp, aes(x=X, fill=MODEL)) +
			geom_boxplot(aes(ymin=CL, lower=IL, upper=IU, ymax=CU, middle=M), stat='identity') +
			theme_bw() + 
			coord_flip() +
			scale_y_continuous(label=scales:::percent) +
			labs(x='Transmissions from -> to\n', y='\nProportion of transmissions by source population')
	ggsave(file=paste0(outfile.base,'_waifm_fishinlandgender.pdf'), w=8, h=5)
	
}

RakaiFull.phylogeography.180521.flows.fishinlandmigrant<- function()
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
	
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_core_inference.rda"
	outfile.base		<- gsub('_inference.rda','',infile.inference)
	load(infile.inference)
	
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')
	
	#
	#	geography who infects whom matrix  between fisherfolk and others account for migration (complex)
	#	adjusted P
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset some TR_INMIGRATE
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_INMIGRATE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',TR_INMIGRATE,' ',REC_COMM_TYPE)]
	#unique(z$FLOW)
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+) ([^ ]+)','\\1',variable)]
	z[, TR_INMIGRATE:= gsub('([^ ]+) ([^ ]+) ([^ ]+)','\\2',variable)]
	z[, REC_COMM_TYPE:= gsub('([^ ]+) ([^ ]+) ([^ ]+)','\\3',variable)]		
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE )
	z[, STAT:='joint']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- copy(z)	
	
	#	make plot
	tmp		<- copy(z)	
	set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','External',gsub('inland','Inland',gsub('fisherfolk','Fishing',TR_COMM_TYPE)))])
	set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','Inland',gsub('fisherfolk','Fishing',REC_COMM_TYPE))])
	set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('_',' ',TR_INMIGRATE)])
	set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('inmigrant from external','',TR_INMIGRATE)])
	tmp[, MODEL:= 'MonteCarloAdjustSampling']	
	tmp[, X:= paste0(TR_COMM_TYPE,'->',REC_COMM_TYPE, ' (', TR_INMIGRATE, ')')]
	ggplot(tmp, aes(x=X, fill=MODEL)) +
			geom_boxplot(aes(ymin=CL, lower=IL, upper=IU, ymax=CU, middle=M), stat='identity') +
			theme_bw() + 
			coord_flip() +
			scale_y_continuous(label=scales:::percent) +
			labs(x='Transmissions from -> to\n', y='\nProportion of transmissions')
	ggsave(file=paste0(outfile.base,'_flows_fishinlandmigration.pdf'), w=8, h=5)
	

	#
	#	WAIFM matrix (complex model)
	#
	groups	<- data.table(TR_COMM_TYPE=c('inland','inland','fisherfolk','fisherfolk','external'), TR_INMIGRATE=c('resident','outmigrant','resident','outmigrant','inmigrant_from_external'))
	z		<- lapply(1:nrow(groups), function(ii)
			{				
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
				#set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				#	reset to complex migrant
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
				set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
				#
				z		<- subset(z, TR_COMM_TYPE==groups$TR_COMM_TYPE[ii] & TR_INMIGRATE==groups$TR_INMIGRATE[ii])				
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
				z[, TR_COMM_TYPE:= groups$TR_COMM_TYPE[ii] ]
				z[, TR_INMIGRATE:= groups$TR_INMIGRATE[ii] ] 
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE)
	z[, STAT:='waifm']	
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	#	plot
	tmp		<- copy(z)	
	set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','External',gsub('inland','Inland',gsub('fisherfolk','Fishing',TR_COMM_TYPE)))])
	set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','Inland',gsub('fisherfolk','Fishing',REC_COMM_TYPE))])
	set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('_',' ',TR_INMIGRATE)])
	set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('inmigrant from external','',TR_INMIGRATE)])
	tmp[, MODEL:= 'MonteCarloAdjustSampling']	
	tmp[, X:= paste0(TR_COMM_TYPE,'->',REC_COMM_TYPE, ' (', TR_INMIGRATE, ')')]
	ggplot(tmp, aes(x=X, fill=MODEL)) +
			geom_boxplot(aes(ymin=CL, lower=IL, upper=IU, ymax=CU, middle=M), stat='identity') +
			theme_bw() + 
			coord_flip() +
			scale_y_continuous(label=scales:::percent) +
			labs(x='Transmissions from -> to\n', y='\nProportion of transmissions by source population')
	ggsave(file=paste0(outfile.base,'_waifm_fishinlandmigration.pdf'), w=8, h=5)
	
	#
	#	sources (complex model)
	#
	groups	<- data.table(REC_COMM_TYPE=c('inland','fisherfolk'))
	z		<- lapply(1:nrow(groups), function(ii)
			{				
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
				#set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				#	reset to complex
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
				set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
				#
				z		<- subset(z, REC_COMM_TYPE==groups$REC_COMM_TYPE[ii])				
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('TR_COMM_TYPE','TR_INMIGRATE','MONTE_CARLO_IT')]
				z[, FLOW:= paste0(TR_COMM_TYPE, ' ',TR_INMIGRATE)]
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- FLOW
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	get quantiles
				z		<- melt(z, id.vars=c('MONTE_CARLO_IT'))								
				z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by='variable']
				z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\1',variable)]
				z[, TR_INMIGRATE:= gsub('([^ ]+) ([^ ]+)','\\2',variable)]				
				z[, REC_COMM_TYPE:= groups$REC_COMM_TYPE[ii] ]				 
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_COMM_TYPE, TR_COMM_TYPE, TR_INMIGRATE)
	z[, STAT:='sources']	
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	
	
	
	#
	#	geography who infects whom matrix  between fisherfolk and others account for migration (simple)
	#	adjusted P
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset 
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	set(z, z[, which(TR_INMIGRATE=='resident')], 'TR_INMIGRATE', 'resident/outmigrant')	
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_INMIGRATE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',TR_INMIGRATE,' ',REC_COMM_TYPE)]
	#unique(z$FLOW)
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+) ([^ ]+)','\\1',variable)]
	z[, TR_INMIGRATE:= gsub('([^ ]+) ([^ ]+) ([^ ]+)','\\2',variable)]
	z[, REC_COMM_TYPE:= gsub('([^ ]+) ([^ ]+) ([^ ]+)','\\3',variable)]		
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE )
	z[, STAT:='joint2']
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	#	plot
	tmp		<- copy(z)	
	set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','External',gsub('inland','Inland',gsub('fisherfolk','Fishing',TR_COMM_TYPE)))])
	set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','Inland',gsub('fisherfolk','Fishing',REC_COMM_TYPE))])
	set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('_',' ',TR_INMIGRATE)])
	set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('^resident$','resident/outmigrant',gsub('inmigrant from external','',TR_INMIGRATE))])
	tmp[, MODEL:= 'MonteCarloAdjustSampling']	
	tmp[, X:= paste0(TR_COMM_TYPE,'->',REC_COMM_TYPE, ' (', TR_INMIGRATE, ')')]
	ggplot(tmp, aes(x=X, fill=MODEL)) +
			geom_boxplot(aes(ymin=CL, lower=IL, upper=IU, ymax=CU, middle=M), stat='identity') +
			theme_bw() + 
			coord_flip() +
			scale_y_continuous(label=scales:::percent) +
			labs(x='Transmissions from -> to\n', y='\nProportion of transmissions')
	ggsave(file=paste0(outfile.base,'_flows_fishinlandmigration2.pdf'), w=8, h=5)
	
	
	#
	#	WAIFM matrix (simple model)
	#
	groups	<- data.table(TR_COMM_TYPE=c('inland','fisherfolk','external'), TR_INMIGRATE=c('resident/outmigrant','resident/outmigrant','inmigrant_from_external'))
	z		<- lapply(1:nrow(groups), function(ii)
			{				
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
				#set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				#	reset to simple
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')	
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
				set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
				set(z, z[, which(TR_INMIGRATE=='resident')], 'TR_INMIGRATE', 'resident/outmigrant')	
				#
				z		<- subset(z, TR_COMM_TYPE==groups$TR_COMM_TYPE[ii] & TR_INMIGRATE==groups$TR_INMIGRATE[ii])				
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
				z[, TR_COMM_TYPE:= groups$TR_COMM_TYPE[ii] ]
				z[, TR_INMIGRATE:= groups$TR_INMIGRATE[ii] ] 
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE)
	z[, STAT:='waifm2']	
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	#	plot
	tmp		<- copy(z)	
	set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','External',gsub('inland','Inland',gsub('fisherfolk','Fishing',TR_COMM_TYPE)))])
	set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','Inland',gsub('fisherfolk','Fishing',REC_COMM_TYPE))])
	set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('_',' ',TR_INMIGRATE)])
	set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('inmigrant from external','',TR_INMIGRATE)])
	tmp[, MODEL:= 'MonteCarloAdjustSampling']	
	tmp[, X:= paste0(TR_COMM_TYPE,'->',REC_COMM_TYPE, ' (', TR_INMIGRATE, ')')]
	ggplot(tmp, aes(x=X, fill=MODEL)) +
			geom_boxplot(aes(ymin=CL, lower=IL, upper=IU, ymax=CU, middle=M), stat='identity') +
			theme_bw() + 
			coord_flip() +
			scale_y_continuous(label=scales:::percent) +
			labs(x='Transmissions from -> to\n', y='\nProportion of transmissions by source population')
	ggsave(file=paste0(outfile.base,'_waifm_fishinlandmigration2.pdf'), w=8, h=5)
	
	
	#
	#	sources (simple model)
	#
	groups	<- data.table(REC_COMM_TYPE=c('inland','fisherfolk'))
	z		<- lapply(1:nrow(groups), function(ii)
			{				
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
				#set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				#	reset to simple
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')	
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
				set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
				set(z, z[, which(TR_INMIGRATE=='resident')], 'TR_INMIGRATE', 'resident/outmigrant')	
				#
				z		<- subset(z, REC_COMM_TYPE==groups$REC_COMM_TYPE[ii])				
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('TR_COMM_TYPE','TR_INMIGRATE','MONTE_CARLO_IT')]
				z[, FLOW:= paste0(TR_COMM_TYPE, ' ',TR_INMIGRATE)]
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- FLOW
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	get quantiles
				z		<- melt(z, id.vars=c('MONTE_CARLO_IT'))								
				z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by='variable']
				z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\1',variable)]
				z[, TR_INMIGRATE:= gsub('([^ ]+) ([^ ]+)','\\2',variable)]				
				z[, REC_COMM_TYPE:= groups$REC_COMM_TYPE[ii] ]				 
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_COMM_TYPE, TR_COMM_TYPE, TR_INMIGRATE)
	z[, STAT:='sources2']	
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	
	#	make table
	df	<- subset(ans, !grepl('sources',STAT))
	df[, STAT2:= gsub('joint2','joint',STAT)]
	set(df, NULL, 'TR_COMM_TYPE', df[, factor(TR_COMM_TYPE, levels=c('inland','fisherfolk','external'))])
	set(df, NULL, 'REC_COMM_TYPE', df[, factor(REC_COMM_TYPE, levels=c('inland','fisherfolk'))])
	set(df, NULL, 'TR_INMIGRATE', df[, factor(TR_INMIGRATE, levels=c('resident/outmigrant','resident','outmigrant','inmigrant_from_external'))])
	df	<- subset(df, !(TR_COMM_TYPE=='external' & STAT=='joint2'))
	df	<- dcast.data.table(df, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~STAT2, value.var='LABEL2')
	setkey(df, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE)
	write.csv(df, row.names=FALSE, file=paste0(outfile.base, '_table_fishinlandmigration.csv'))
	
	#	save
	save(ans, file=paste0(outfile.base,'_flows_fishinlandmigration.rda'))

	#	make alluvial sources plot
	require(ggalluvial)
	df	<- subset(ans, STAT=='joint2', select=c(TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE, STAT, M))
	tmp	<- subset(ans, STAT=='sources2', select=c(TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE, STAT, M, LABEL2))
	ggplot(df, aes(weight= M, axis1= TR_COMM_TYPE, axis2 = REC_COMM_TYPE)) +
			geom_alluvium(aes(fill = M, width = 2/12), alpha=0.8) +
			geom_stratum(width = 2/12, fill='grey', color="black") +
			geom_label(stat = "stratum", label.strata = TRUE) +
			geom_label(stat = "flow", label=rep(tmp$LABEL2,2)) +
			scale_x_continuous(breaks = 1:2, labels = c("Source population", "Recipient population")) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0)) +
			scale_fill_gradient2(trans='log', breaks=c(0.05, 0.1, 0.2, 0.4), low="deepskyblue", mid="orange", high='red', midpoint=log(0.15)) +
			#scale_fill_manual(values=c('external'='grey50','fisherfolk'='firebrick4','inland'='darkgreen')) +
			theme_minimal()
	ggsave(file=paste0(outfile.base,'_alluvial_fishinlandmigration2.pdf'), w=5, h=8)
	#	https://stackoverflow.com/questions/43053375/weighted-sankey-alluvial-diagram-for-visualizing-discrete-and-continuous-panel
	#dff	<- melt(df, id.vars=c('TR_INMIGRATE'))
	#ggplot(data=d, aes(x = timeperiod, stratum = discretechoice, alluvium = individual, weight = continuouschoice)) +
	#			geom_stratum(aes(fill = discretechoice)) +
	#			geom_flow()
	
	#	make boxplot
	tmp	<- subset(ans, STAT=='joint2' | STAT=='waifm2', select=c(TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE, STAT, M, CL, CU, IL, IU))
	set(tmp, NULL, 'STAT', tmp[, gsub('joint2','Transmission flow\noverall',gsub('waifm2','Transmission flow\nby source',STAT))])
	set(tmp, NULL, 'STAT', tmp[, factor(STAT, levels=c('Transmission flow\noverall','Transmission flow\nby source'))])
	set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','From\nExternal',gsub('inland','From\nInland',gsub('fisherfolk','From\nFishing',TR_COMM_TYPE)))])
	set(tmp, NULL, 'TR_COMM_TYPE', tmp[, factor(TR_COMM_TYPE, levels=c('From\nInland','From\nFishing','From\nExternal'))])
	set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','To Inland',gsub('fisherfolk','To Fishing',REC_COMM_TYPE))])
	set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('_',' ',TR_INMIGRATE)])
	set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('inmigrant from external','',TR_INMIGRATE)])
	#tmp[, X:= paste0(TR_COMM_TYPE,'->',REC_COMM_TYPE, ' (', TR_INMIGRATE, ')')]
	tmp[, X:= REC_COMM_TYPE]
	ggplot(tmp, aes(x=X)) +
			geom_boxplot(aes(ymin=CL, lower=IL, upper=IU, ymax=CU, middle=M, fill=TR_COMM_TYPE), stat='identity') +
			theme_bw() + 			
			scale_y_continuous(label=scales:::percent) +
			scale_fill_manual(values=c('From\nExternal'='grey50', 'From\nInland'='DarkGreen', 'From\nFishing'='firebrick1')) +
			coord_flip() +
			labs(x='', y='\nHIV transmissions between RCCS communities and external locations') +
			facet_grid(TR_COMM_TYPE~STAT, scales='free', switch='y') +
			theme(strip.text.y = element_text(angle=180), strip.placement = "outside") +
			guides(fill='none')
	ggsave(file=paste0(outfile.base,'_alluvial_fishinlandmigration_jointwaifm.pdf'), w=8, h=5)
	
	#	make barplot
	tmp	<- subset(ans, STAT=='joint2' | STAT=='waifm2', select=c(TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE, STAT, M, CL, CU, IL, IU))
	set(tmp, NULL, 'STAT', tmp[, gsub('joint2','Transmission flow\noverall',gsub('waifm2','Transmission flow\nby source',STAT))])
	set(tmp, NULL, 'STAT', tmp[, factor(STAT, levels=c('Transmission flow\noverall','Transmission flow\nby source'))])
	set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','From\nExternal',gsub('inland','From\nInland',gsub('fisherfolk','From\nFishing',TR_COMM_TYPE)))])
	set(tmp, NULL, 'TR_COMM_TYPE', tmp[, factor(TR_COMM_TYPE, levels=c('From\nInland','From\nFishing','From\nExternal'))])
	set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','To\nInland',gsub('fisherfolk','To\nFishing',REC_COMM_TYPE))])
	set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('_',' ',TR_INMIGRATE)])
	set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('inmigrant from external','',TR_INMIGRATE)])
	#tmp[, X:= paste0(TR_COMM_TYPE,'->',REC_COMM_TYPE, ' (', TR_INMIGRATE, ')')]
	tmp[, X:= REC_COMM_TYPE]
	ggplot(tmp, aes(x=X)) +
			geom_bar(aes(y=M, fill=TR_COMM_TYPE), stat='identity') +
			geom_errorbar(aes(ymin=CL, ymax=CU), width=0.5) +
			theme_bw() + 			
			scale_y_continuous(label=scales:::percent) +
			scale_fill_manual(values=c('From\nExternal'='grey50', 'From\nInland'='DarkGreen', 'From\nFishing'='firebrick1')) +			
			labs(x='', y='HIV transmissions\nbetween RCCS communities and external locations\n') +
			facet_grid(STAT~TR_COMM_TYPE, scales='free') +			
			guides(fill='none')
	ggsave(file=paste0(outfile.base,'_alluvial_fishinlandmigration_jointwaifm_bar.pdf'), w=5, h=8)
	
}	
	

RakaiFull.phylogeography.180521.flows.fishinland<- function()
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
	
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_core_inference.rda"
	outfile.base		<- gsub('_inference.rda','',infile.inference)
	load(infile.inference)
	
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')
	
	#
	#	geography who infects whom matrix  between fisherfolk and others
	#	adjusted P
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
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
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, TR_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\1',variable)]
	z[, REC_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\2',variable)]		
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_COMM_TYPE, TR_COMM_TYPE)
	z[, STAT:='joint']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- copy(z)	
	
	require(rstan)
	require(rethinking)
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_stanm1m4m5.rda")
	#	question: how different are estimates from m2 m4 m5?
	post	<- extract.samples(m1$fit)
	tmp		<- data.table(FROM_TO_ID=1:max(dc1$FROM_TO_ID))
	df1		<- tmp[, list(MCIT=seq_len(nrow(post$prop_comm_flow)), PI=post$prop_comm_flow[,FROM_TO_ID]), by='FROM_TO_ID']
	df1		<- df1[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by='FROM_TO_ID']
	df1		<- merge(dc1, df1, by='FROM_TO_ID')	
	df1[, MODEL:='area-based']
	df1 	<- dcast.data.table(df1, MODEL+TR_FISHCOMM+REC_FISHCOMM~P, value.var='Q')
	
	post	<- extract.samples(m5$fit)
	tmp		<- data.table(FROM_TO_ID=1:max(dc5$FROM_TO_ID))
	df5		<- tmp[, list(MCIT=seq_len(nrow(post$prop_comm_flow)), PI=post$prop_comm_flow[,FROM_TO_ID]), by='FROM_TO_ID']	
	df5		<- merge(dc5, df5, by='FROM_TO_ID')
	df5		<- df5[, list(PI=sum(PI)), by=c('MCIT','TR_FISHCOMM','REC_FISHCOMM')]
	df5		<- df5[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by=c('TR_FISHCOMM','REC_FISHCOMM')]
	df5[, MODEL:='community-based prior Dir(0.8/NPAIR)']
	df5 	<- dcast.data.table(df5, MODEL+TR_FISHCOMM+REC_FISHCOMM~P, value.var='Q')
	
	tmp		<- copy(z)
	setnames(tmp, c('TR_COMM_TYPE','REC_COMM_TYPE'), c('TR_FISHCOMM','REC_FISHCOMM'))
	set(tmp, NULL, 'TR_FISHCOMM', tmp[, gsub('from_inland','Inland',gsub('from_fisherfolk','Fishing',TR_FISHCOMM))])
	set(tmp, NULL, 'REC_FISHCOMM', tmp[, gsub('to_inland','Inland',gsub('to_fisherfolk','Fishing',REC_FISHCOMM))])
	tmp[, MODEL:= 'MonteCarloAdjustSampling']
	tmp		<- rbind(tmp, df1, df5, fill=TRUE)
	
	ggplot(tmp, aes(x=paste0(TR_FISHCOMM,'->',REC_FISHCOMM), fill=MODEL)) +
			geom_boxplot(aes(ymin=CL, lower=IL, upper=IU, ymax=CU, middle=M), stat='identity') +
			theme_bw() + 
			coord_flip() +
			scale_y_continuous(label=scales:::percent) +
			labs(x='Transmissions from -> to\n', y='\nProportion of transmissions')
	ggsave(file=paste0(outfile.base,'_compare_m1_m5_sampling.pdf'), w=8, h=5)
		
	
	#
	#	WAIFM
	#
	groups	<- c('inland','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))				
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
	
	
	
	#
	#	sources
	#
	groups	<- c('inland','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))				
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
	
	save(ans, file=paste0(outfile.base,'_fishing_inland_results.rda'))
}

RakaiFull.phylogeography.180521.gender.mobility.core.inference<- function()
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
	
	infile			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda"
	outfile.base	<- gsub('data_with_inmigrants.rda','',infile)
	load(infile)
	
	#	select which inmigration data to work with
	setnames(rtr2, c('TR_INMIGRATE_2YR','REC_INMIGRATE_2YR'), c('TR_INMIGRATE','REC_INMIGRATE'))
	
	#	sum transmission events by community and gender and inmigration status
	dc	<- rtr2[, list(TR_OBS=length(PAIRID)), by=c('TR_COMM_NUM_A','REC_COMM_NUM_A','TR_SEX','REC_SEX', 'TR_INMIGRATE')]
	set(dc, NULL, 'TR_COMM_NUM_A', dc[, as.character(TR_COMM_NUM_A)])
	set(dc, NULL, 'REC_COMM_NUM_A', dc[, as.character(REC_COMM_NUM_A)])	
	setkey(dc, TR_COMM_NUM_A, REC_COMM_NUM_A, TR_SEX, TR_INMIGRATE)
	# 	use the Berger objective prior with minimal loss compared to marginal Beta reference prior
	#	(https://projecteuclid.org/euclid.ba/1422556416)
	dc[, TR_PRIOR:= 0.8/nrow(dc)]
	
	
	#
	#	Bayesian model sampling prior: 
	#	define Beta prior for sampling probabilities as Empirical Bayes prior from cohort data (all alpha and betas)
	#	currently sampling probs differ by community and sex
	#
	tmp	<- subset(ds, select=c(COMM_NUM_A, SEX, P_PART_ALPHA, P_PART_BETA, P_SEQ_ALPHA, P_SEQ_BETA))
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc	<- merge(dc, tmp, by=c('TR_COMM_NUM_A','TR_SEX'))
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc	<- merge(dc, tmp, by=c('REC_COMM_NUM_A','REC_SEX'))
	
	
	#
	#	Bayesian model missing data: draw unobserved data to augment likelihood
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
				tmp	<- rnbinom(mc.it, TR_OBS, tmp)
				#print(tmp)
				#stop()
				list(	MONTE_CARLO_IT=seq_len(mc.it), 
						TR_PRIOR=TR_PRIOR, 
						TR_OBS=TR_OBS, 
						TR_MISS= tmp,
						TR_MISS_P= TR_P_PART_ALPHA/(TR_P_PART_ALPHA+TR_P_PART_BETA) * TR_P_SEQ_ALPHA/(TR_P_SEQ_ALPHA+TR_P_SEQ_BETA) * REC_P_PART_ALPHA/(REC_P_PART_ALPHA+REC_P_PART_BETA) * REC_P_SEQ_ALPHA/(REC_P_SEQ_ALPHA+REC_P_SEQ_BETA)   )
			}, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A','REC_SEX','TR_SEX','TR_INMIGRATE')]
	#	check
	#tmp	<- dcb[, list(	AUG_E= (TR_PRIOR[1]+TR_OBS[1])/TR_MISS_P[1],
	#					AUG_MEAN= mean(TR_MISS)+TR_PRIOR[1]+TR_OBS[1]
	#					), by=c('REC_COMM_NUM_A','TR_COMM_NUM_A','REC_SEX','TR_SEX')]
	#tmp[, summary(AUG_E-AUG_MEAN)]
	#
	#	Bayesian model augmented posterior: 
	#	the augmented likelihood is multinomial, so the augmented posterior is 
	#	analytically tractable, and a Dirichlet posterior
	#
	tmp		<- dcb[, list(	REC_COMM_NUM_A= REC_COMM_NUM_A, 
					TR_COMM_NUM_A= TR_COMM_NUM_A, 
					REC_SEX=REC_SEX,
					TR_SEX=TR_SEX,
					TR_INMIGRATE=TR_INMIGRATE,
					PI_IJ_ALPHA= TR_OBS+TR_MISS+TR_PRIOR				
			), by='MONTE_CARLO_IT']
	dcb		<- merge(dcb, tmp, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A','REC_SEX','TR_SEX','TR_INMIGRATE','MONTE_CARLO_IT'))
	#
	#	this is the end of the source attribution inference on the WAIFM matrix
	#
	save(zm, rtpdm, rtr2, ds, dc, dcb, file=paste0(outfile.base,'core_inference.rda'))	
}

mcmc.n.z.pi.eachcout.ok<- function(zm, rtpdm, rtr2, ds, dc, outfile.base)
{
	
	# set up mcmc objects
	mc				<- list()
	mc$n			<- 1e6
	mc$pars			<- list() 	
	mc$pars$LAMBDA	<- matrix(NA_real_, ncol=length(unique(dc$COUNT_ID)), nrow=1)		#prior for proportions
	mc$pars$S		<- matrix(NA_real_, ncol=length(unique(dc$COUNT_ID)), nrow=mc$n)	#sampling probabilities
	mc$pars$Z		<- matrix(NA_integer_, ncol=length(unique(dc$COUNT_ID)), nrow=mc$n) #augmented data
	mc$pars$NU		<- NA_real_															#prior for N
	mc$pars$N		<- matrix(NA_integer_, ncol=1, nrow=mc$n)							#total number of counts on augmented data
	mc$pars$PI		<- matrix(NA_real_, ncol=length(unique(dc$COUNT_ID)), nrow=mc$n)	#proportions	
	mc$it.info		<- data.table(	IT= seq.int(1,mc$n),
			BLOCK= rep(NA_integer_, mc$n),
			MHRATIO= rep(NA_real_, mc$n),
			ACCEPT=rep(NA_integer_, mc$n), 
			LOG_LKL=rep(NA_real_, mc$n),
			LOG_PRIOR=rep(NA_real_, mc$n))
	
	# initialise MCMC
	setkey(dc, COUNT_ID)
	mc$verbose			<- 0L
	mc$curr.it			<- 1L
	mc$seed				<- 42L
	set(mc$it.info, mc$curr.it, 'BLOCK', 'INIT')
	set.seed( mc$seed )
	#	prior lambda: use the Berger objective prior with minimal loss compared to marginal Beta reference prior
	#	(https://projecteuclid.org/euclid.ba/1422556416)	
	mc$pars$LAMBDA[1,]	<- 0.8/nrow(dc)
	# 	sampling: set to mean sampling probability
	mc$pars$S[1,]		<- dc$S						
	#	augmented data: proposal draw under mean sampling probability
	mc$pars$Z[1,]		<- dc[, TR_OBS + rnbinom(nrow(dc), TR_OBS, S)]	# dc[, round(TR_OBS/S)]
	#	prior nu: set Poisson rate to the expected augmented counts, under average sampling probability
	mc$pars$NU			<- sum(dc$TR_OBS) / mean(dc$S)
	#	total count: that s just the sum of Z
	mc$pars$N[1,]		<- sum(mc$pars$Z[1,])
	#	proportions: draw from full conditional
	mc$pars$PI[1,]		<- rdirichlet(1, mc$pars$Z[1,] + mc$pars$LAMBDA[1,])[1,]			
	#	store log likelihood
	tmp	<- sum( dbinom(dc$TR_OBS, size=mc$pars$Z[1,], prob=dc$S, log=TRUE) ) +
			dmultinom(mc$pars$Z[1,], size=mc$pars$N[1,], prob=mc$pars$PI[1,], log=TRUE)
	set(mc$it.info, 1L, 'LOG_LKL', tmp)
	# 	store log prior		
	tmp	<- dpois(mc$pars$N[1,], lambda=mc$pars$NU, log=TRUE) +
			sum( log( ddirichlet(mc$pars$PI[1,], alpha=mc$pars$LAMBDA[1,]) ) )
	set(mc$it.info, 1L, 'LOG_PRIOR', tmp)
	
	#	
	# run mcmc
	mc$blocks		<- ncol(mc$pars$Z)
	for(i in 1L:(mc$n-1L))
	{
		mc$curr.it		<- i		
		# determine source-recipient combination that will be updated in this iteration
		update.count	<- i %% mc$blocks
		
		# update Z, N, PI for the source-recipient combination 'update.count'
		#	propose
		Z.prop					<- mc$pars$Z[mc$curr.it,]
		Z.prop[update.count]	<- dc[COUNT_ID==update.count, TR_OBS + rnbinom(1, TR_OBS, S)]
		N.prop					<- sum(Z.prop)
		PI.prop					<- mc$pars$PI[mc$curr.it,]
		tmp						<- c(Z.prop[update.count]+mc$pars$LAMBDA[1,update.count], sum(Z.prop + mc$pars$LAMBDA[1,]))
		PI.prop[update.count]	<- rbeta(1, tmp[1], tmp[2]-tmp[1])							
		
		#	calculate MH ratio
		log.prop.ratio	<- sum(dnbinom(mc$pars$Z[mc$curr.it,update.count], size=dc$TR_OBS[update.count], prob=dc$S[update.count], log=TRUE)) - 
				sum(dnbinom(Z.prop[update.count], size=dc$TR_OBS[update.count], prob=dc$S[update.count], log=TRUE))
		log.fc			<- sum(dbinom(dc$TR_OBS[update.count], size=mc$pars$Z[mc$curr.it,update.count], prob=dc$S[update.count], log=TRUE)) +								
				dpois(mc$pars$N[mc$curr.it,], lambda=mc$pars$NU, log=TRUE)
		log.fc.prop		<- sum(dbinom(dc$TR_OBS[update.count], size=Z.prop[update.count], prob=dc$S[update.count], log=TRUE)) +								
				dpois(N.prop, lambda=mc$pars$NU, log=TRUE)
		log.mh.ratio	<- log.fc.prop - log.fc + log.prop.ratio
		mh.ratio		<- min(1,exp(log.mh.ratio))
		
		#	update
		mc$curr.it				<- mc$curr.it+1L
		set(mc$it.info, mc$curr.it, 'BLOCK', update.count)
		set(mc$it.info, mc$curr.it, 'MHRATIO', mh.ratio)
		set(mc$it.info, mc$curr.it, 'ACCEPT', as.integer(runif(1) < mh.ratio))							
		mc$pars$S[mc$curr.it,]	<- mc$pars$S[mc$curr.it-1L,]
		if(mc$verbose & mc$it.info[mc$curr.it, ACCEPT])
		{
			print(paste0('it ',mc$curr.it,' ACCEPT block ',update.count))
		}
		if(mc$it.info[mc$curr.it, ACCEPT])
		{
			mc$pars$Z[mc$curr.it,]	<- Z.prop
			mc$pars$N[mc$curr.it,]	<- N.prop		
			mc$pars$PI[mc$curr.it,]	<- PI.prop
		}
		if(mc$it.info[mc$curr.it, !ACCEPT])
		{
			mc$pars$Z[mc$curr.it,]	<- mc$pars$Z[mc$curr.it-1L,]
			mc$pars$N[mc$curr.it,]	<- mc$pars$N[mc$curr.it-1L,]
			mc$pars$PI[mc$curr.it,]	<- mc$pars$PI[mc$curr.it-1L,]
		}			
		
		# store log likelihood
		tmp	<- sum( dbinom(dc$TR_OBS, size=mc$pars$Z[mc$curr.it,], prob=dc$S, log=TRUE) ) +
				dmultinom(mc$pars$Z[mc$curr.it,], size=mc$pars$N[mc$curr.it,], prob=mc$pars$PI[mc$curr.it,], log=TRUE)
		set(mc$it.info, mc$curr.it, 'LOG_LKL', tmp)
		# store log prior		
		tmp	<- dpois(mc$pars$N[mc$curr.it,], lambda=mc$pars$NU, log=TRUE) +
				sum( log( ddirichlet(mc$pars$PI[mc$curr.it,], alpha=mc$pars$LAMBDA[1,]) ) )
		set(mc$it.info, mc$curr.it, 'LOG_PRIOR', tmp)		
	}
	save(zm, rtpdm, rtr2, ds, dc, mc, file=paste0(outfile.base,'core_inference_mcmcEachCount.rda'))
	
	
	#	to manage all further calculations, thin right at the start
	tmp			<- seq.int(2, nrow(mc$pars$Z), 100)
	mc$pars$S	<- mc$pars$S[tmp,,drop=FALSE]
	mc$pars$Z	<- mc$pars$Z[tmp,,drop=FALSE]
	mc$pars$PI	<- mc$pars$PI[tmp,,drop=FALSE]
	mc$pars$N	<- mc$pars$N[tmp,,drop=FALSE]
	gc()
	colnames(mc$pars$S)	<- paste0('S-',1:ncol(mc$pars$S))
	colnames(mc$pars$Z)	<- paste0('Z-',1:ncol(mc$pars$Z))
	colnames(mc$pars$PI)<- paste0('PI-',1:ncol(mc$pars$PI))
	colnames(mc$pars$N)	<- 'N'
	#	convert to coda
	mcc			<- mcmc( cbind( mc$pars$N, mc$pars$Z, mc$pars$PI, mc$pars$S) )
	mcc.eff		<- effectiveSize(mcc)
	#	plot trace and autocorrelations for variable with worst 6 effective sizes
	tmp			<- mcc.eff[mcc.eff>0] 
	tmp			<- names(tmp[sort(tmp, index.return=TRUE)$ix[1:6]])	
	pdf(file=paste0(outfile.base,'core_inference_mcmcEachCount_6worstchains_trace.pdf'), w=10, h=7)
	plot( mcc[,tmp] )
	dev.off()
	pdf(file=paste0(outfile.base,'core_inference_mcmcEachCount_6worstchains_autocor.pdf'), w=5, h=7)
	autocorr.plot( mcc[,tmp] )
	dev.off()
}	

mcmc.n.z.pi.joint.low.acceptance.rate<- function(zm, rtpdm, rtr2, ds, dc, outfile.base)
{
	# set up mcmc objects
	mc				<- list()
	mc$n			<- 1e5
	mc$pars			<- list() 	
	mc$pars$LAMBDA	<- matrix(NA_real_, ncol=length(unique(dc$COUNT_ID)), nrow=1)		#prior for proportions
	mc$pars$S		<- matrix(NA_real_, ncol=length(unique(dc$COUNT_ID)), nrow=mc$n)	#sampling probabilities
	mc$pars$Z		<- matrix(NA_integer_, ncol=length(unique(dc$COUNT_ID)), nrow=mc$n) #augmented data
	mc$pars$NU		<- NA_real_															#prior for N
	mc$pars$N		<- matrix(NA_integer_, ncol=1, nrow=mc$n)							#total number of counts on augmented data
	mc$pars$PI		<- matrix(NA_real_, ncol=length(unique(dc$COUNT_ID)), nrow=mc$n)	#proportions	
	mc$it.info		<- data.table(	IT= seq.int(1,mc$n),
			BLOCK= rep(NA_character_, mc$n),
			MHRATIO= rep(NA_real_, mc$n),
			ACCEPT=rep(NA_integer_, mc$n), 
			LOG_LKL=rep(NA_real_, mc$n),
			LOG_PRIOR=rep(NA_real_, mc$n))
	
	# initialise MCMC
	setkey(dc, COUNT_ID)
	mc$verbose			<- 1L
	mc$curr.it			<- 1L
	mc$seed				<- 42L
	set(mc$it.info, mc$curr.it, 'BLOCK', 'INIT')
	set.seed( mc$seed )
	#	prior lambda: use the Berger objective prior with minimal loss compared to marginal Beta reference prior
	#	(https://projecteuclid.org/euclid.ba/1422556416)	
	mc$pars$LAMBDA[1,]	<- 0.8/nrow(dc)
	# 	sampling: set to mean sampling probability
	mc$pars$S[1,]		<- dc$S						
	#	augmented data: proposal draw under mean sampling probability
	mc$pars$Z[1,]		<- dc[, TR_OBS + rnbinom(nrow(dc), TR_OBS, S)]	# dc[, round(TR_OBS/S)]
	#	prior nu: set Poisson rate to the expected augmented counts, under average sampling probability
	mc$pars$NU			<- sum(dc$TR_OBS) / mean(dc$S)
	#	total count: that s just the sum of Z
	mc$pars$N[1,]		<- sum(mc$pars$Z[1,])
	#	proportions: draw from full conditional
	mc$pars$PI[1,]		<- rdirichlet(1, mc$pars$Z[1,] + mc$pars$LAMBDA[1,])[1,]			
	#	store log likelihood
	tmp	<- sum( dbinom(dc$TR_OBS, size=mc$pars$Z[1,], prob=dc$S, log=TRUE) ) +
			dmultinom(mc$pars$Z[1,], size=mc$pars$N[1,], prob=mc$pars$PI[1,], log=TRUE)
	set(mc$it.info, 1L, 'LOG_LKL', tmp)
	# 	store log prior		
	tmp	<- dpois(mc$pars$N[1,], lambda=mc$pars$NU, log=TRUE) +
			sum( log( ddirichlet(mc$pars$PI[1,], alpha=mc$pars$LAMBDA[1,]) ) )
	set(mc$it.info, 1L, 'LOG_PRIOR', tmp)
	
	#	
	# run mcmc
	mc$blocks		<- 1
	for(i in 1L:(mc$n-1L))
	{
		mc$curr.it		<- i		
		#
		# block 0: update Z, N, PI jointly
		if(i %% mc$blocks==0)
		{			
			#	propose
			Z.prop			<- dc[, TR_OBS + rnbinom(nrow(dc), TR_OBS, S)]
			N.prop			<- sum(Z.prop)
			PI.prop			<- rdirichlet(1, Z.prop + mc$pars$LAMBDA[1,])[1,]
			#	calculate MH ratio
			log.prop.ratio	<- sum(dnbinom(mc$pars$Z[mc$curr.it,], size=dc$TR_OBS, prob=dc$S, log=TRUE)) - sum(dnbinom(Z.prop, size=dc$TR_OBS, prob=dc$S, log=TRUE))
			log.fc			<- sum(dbinom(dc$TR_OBS, size=mc$pars$Z[mc$curr.it,], prob=dc$S, log=TRUE)) +								
					dpois(mc$pars$N[mc$curr.it,], lambda=mc$pars$NU, log=TRUE)
			log.fc.prop		<- sum(dbinom(dc$TR_OBS, size=Z.prop, prob=dc$S, log=TRUE)) +								
					dpois(N.prop, lambda=mc$pars$NU, log=TRUE)
			log.mh.ratio	<- log.fc.prop - log.fc + log.prop.ratio
			mh.ratio		<- min(1,exp(log.mh.ratio))
			#	update
			mc$curr.it				<- mc$curr.it+1L
			set(mc$it.info, mc$curr.it, 'BLOCK', 'Z-N-PI')
			set(mc$it.info, mc$curr.it, 'MHRATIO', mh.ratio)
			set(mc$it.info, mc$curr.it, 'ACCEPT', as.integer(runif(1) < mh.ratio))							
			mc$pars$S[mc$curr.it,]	<- mc$pars$S[mc$curr.it-1L,]
			if(mc$verbose & mc$it.info[mc$curr.it, ACCEPT])
			{
				print(paste0('it ',mc$curr.it,' ACCEPT block Z-N-PI'))
			}
			if(mc$it.info[mc$curr.it, ACCEPT])
			{
				mc$pars$Z[mc$curr.it,]	<- Z.prop
				mc$pars$N[mc$curr.it,]	<- N.prop		
				mc$pars$PI[mc$curr.it,]	<- PI.prop
			}
			if(mc$it.info[mc$curr.it, !ACCEPT])
			{
				mc$pars$Z[mc$curr.it,]	<- mc$pars$Z[mc$curr.it-1L,]
				mc$pars$N[mc$curr.it,]	<- mc$pars$N[mc$curr.it-1L,]
				mc$pars$PI[mc$curr.it,]	<- mc$pars$PI[mc$curr.it-1L,]
			}						
		}		
		# store log likelihood
		tmp	<- sum( dbinom(dc$TR_OBS, size=mc$pars$Z[mc$curr.it,], prob=dc$S, log=TRUE) ) +
				dmultinom(mc$pars$Z[mc$curr.it,], size=mc$pars$N[mc$curr.it,], prob=mc$pars$PI[mc$curr.it,], log=TRUE)
		set(mc$it.info, mc$curr.it, 'LOG_LKL', tmp)
		# store log prior		
		tmp	<- dpois(mc$pars$N[mc$curr.it,], lambda=mc$pars$NU, log=TRUE) +
				sum( log( ddirichlet(mc$pars$PI[mc$curr.it,], alpha=mc$pars$LAMBDA[1,]) ) )
		set(mc$it.info, mc$curr.it, 'LOG_PRIOR', tmp)		
	}
}

RakaiFull.phylogeography.181006.gender.mobility.data<- function(infile.inference)
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)	
	
	opt.set.missing.migloc.to.inland	<- 0
	opt.set.missing.migloc.to.fishing	<- 1
	
	hivc.db.Date2numeric<- function( x )
	{
		if(!class(x)%in%c('Date','character'))	return( x )
		x	<- as.POSIXlt(x)
		tmp	<- x$year + 1900
		x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
		x	
	}	
	
	#	load des which contains participation and seq counts by comm and gender
	infile	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender_age.rda"
	load(infile)
	setnames(desm, 'INMIG_2YRS_LOC', 'INMIG_LOC')
	#	handle missing data on migration
	if(opt.set.missing.migloc.to.inland)
	{
		set(desm, desm[, which(INMIG_LOC=='origin_unknown')], 'INMIG_LOC', 'inland')
	}
	if(opt.set.missing.migloc.to.fishing)
	{
		set(desm, desm[, which(INMIG_LOC=='origin_unknown')], 'INMIG_LOC', 'fisherfolk')
	}
	ds		<- desm[, list(PART_EVER=sum(PART_EVER), PART_NEVER=sum(PART_NEVER), HIV_1516_YES=sum(HIV_1516_YES), 
							HIV_1516_NO=sum(HIV_1516_NO), SLART_AT_FIRST_VISIT=sum(SLART_AT_FIRST_VISIT), 
							DEEP_SEQ_1516=sum(DEEP_SEQ_1516)), by=c('COMM_NUM','SEX','INMIG_LOC','COMM_NUM_A','COMM_TYPE')]
	#	determine posterior parameters for Binomial models of sampling and participiation 
	ds[, P_PART_ALPHA:= PART_EVER+1]
	ds[, P_PART_BETA:= PART_NEVER+1]
	ds[, P_SEQ_ALPHA:= DEEP_SEQ_1516+1]
	ds[, P_SEQ_BETA:= HIV_1516_YES-DEEP_SEQ_1516+1]
	set(ds, NULL, c('PART_EVER','PART_NEVER','HIV_1516_YES','HIV_1516_NO','SLART_AT_FIRST_VISIT','DEEP_SEQ_1516'), NULL)
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
	
	
	#	prepare inmigrant -- identify inmigrants from fishing communities and from external
	infile		<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData_v2.rda"
	load(infile)
	inmigrant	<- as.data.table(inmigrant)
	infile		<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/migrants_withMissingGPS.csv"
	inmigrant2	<- as.data.table(read.csv(infile))
	
	
	#	plot fisherfolk	to figure out how much of a radius we need
	if(0)
	{
		zf		<- data.table(longitude=c(31.763,31.7968,31.754,31.838), latitude=c(-0.915, -0.6518, -0.703, -0.497), ID= c('Kasensero','Bukyanju','NearBwende','Fish4'))
		make_circles <- function(centers, radius, nPoints = 100){
			# centers: the data frame of centers with ID
			# radius: radius measured in kilometer
			#
			meanLat <- mean(centers$latitude)
			# length per longitude changes with lattitude, so need correction
			radiusLon <- radius /111 / cos(meanLat/57.3) 
			radiusLat <- radius / 111
			circleDF <- data.frame(ID = rep(centers$ID, each = nPoints))
			angle <- seq(0,2*pi,length.out = nPoints)
			
			circleDF$lon <- unlist(lapply(centers$longitude, function(x) x + radius * cos(angle)))
			circleDF$lat <- unlist(lapply(centers$latitude, function(x) x + radius * sin(angle)))
			return(circleDF)
		}
		zc <- make_circles(zf, 0.01)
		ggmap(zm) +			
				geom_point(data=zf, aes(x=longitude, y=latitude, pch=ID), stroke=1.5, alpha=0.8) +
				geom_polygon(data=zc, aes(lon, lat, group = ID), color = "red", alpha = 0)
		#	radius of length 0.01 should catch					
		tmp		<- inmigrant[, list( 	DIST_KASENSERO= sqrt( (inmig_lon- 31.763)^2 + (inmig_lat - (-0.915))^2),
						DIST_BUKYANJU= sqrt( (inmig_lon- 31.7968)^2 + (inmig_lat - (-0.6518))^2),
						DIST_NEARBWENDE= sqrt( (inmig_lon- 31.754)^2 + (inmig_lat - (-0.703))^2),
						DIST_FISH4= sqrt( (inmig_lon- 31.838)^2 + (inmig_lat - (-0.497))^2)
				), by=c('RCCS_studyid','visit')]
		tmp		<- melt(tmp, id.vars=c('RCCS_studyid','visit'))
		ggplot(subset(tmp, value<0.3), aes(x=value)) +
				geom_histogram(binwidth=0.01) +
				facet_grid(variable~.)		
		zfd		<- merge(inmigrant, subset(tmp, value<0.01, c(RCCS_studyid, visit)), by=c('RCCS_studyid','visit'))		
	}
	#	so fishing sites are MALEMBO DIMU KASENSERO NAMIREMBE but there are spelling mistakes
	#	clean up inmigrant	
	#
	#	inmigrant[, unique(sort(inmig_place))]
	set(inmigrant, NULL, 'inmig_place', inmigrant[, gsub('DDIMO|DDIMU|DIMO|DIMU','DIMU',inmig_place)])
	set(inmigrant, inmigrant[, which(grepl('MALEMBO',inmig_place))], 'inmig_place', 'MALEMBO')
	set(inmigrant, NULL, 'inmig_place', inmigrant[, gsub("KASEMSERO","KASENSERO",inmig_place)])
	set(inmigrant, inmigrant[, which(grepl('KASENSERO',inmig_place))], 'inmig_place', 'KASENSERO')
	set(inmigrant2, inmigrant2[, which(grepl('MALEMBO',inmig_place))], 'inmig_place', 'MALEMBO')
	
	#	define from_fishing and from_outside and from_inland
	inmigrant[, INMIG_LOC:= 'inland' ]
	set(inmigrant, inmigrant[, which(grepl('MALEMBO|DIMU|KASENSERO|NAMIREMBE',inmig_place))], 'INMIG_LOC','fisherfolk')
	set(inmigrant, inmigrant[, which(inmig_admin0!='Uganda')], 'INMIG_LOC','external')
	set(inmigrant, inmigrant[, which(inmig_admin1!='Rakai')], 'INMIG_LOC','external')
	set(inmigrant, inmigrant[, which(is.na(inmig_admin1))], 'INMIG_LOC','unknown')	
	inmigrant2[, INMIG_LOC:= 'inland' ]
	set(inmigrant2, inmigrant2[, which(grepl('MALEMBO|DIMU|KASENSERO|NAMIREMBE',inmig_place))], 'INMIG_LOC','fisherfolk')
	set(inmigrant2, inmigrant2[, which(inmig_admin0!='Uganda')], 'INMIG_LOC','external')
	set(inmigrant2, inmigrant2[, which(inmig_admin1!='Rakai')], 'INMIG_LOC','external')
	set(inmigrant2, inmigrant2[, which(is.na(inmig_admin1))], 'INMIG_LOC','unknown')
	inmigrant2[, table(INMIG_LOC)]
	if(opt.set.missing.migloc.to.inland)
	{
		set(inmigrant2, inmigrant2[, which(INMIG_LOC=='unknown')], 'INMIG_LOC', 'inland')
	}
	if(opt.set.missing.migloc.to.fishing)
	{
		set(inmigrant2, inmigrant2[, which(INMIG_LOC=='unknown')], 'INMIG_LOC', 'fisherfolk')
	}
	setnames(inmigrant2, 'INMIG_LOC', 'INMIG_LOC2')
	inmigrant2	<- subset(inmigrant2, select=c('RCCS_studyid','visit','INMIG_LOC2'))
	inmigrant	<- merge(inmigrant, inmigrant2, by=c('RCCS_studyid','visit'), all.x=TRUE)
	tmp			<- inmigrant[, which(!is.na(INMIG_LOC2))]
	set(inmigrant, tmp, 'INMIG_LOC', inmigrant[tmp, INMIG_LOC2])
	
	inmigrant[, table(INMIG_LOC)]
	#	  external fisherfolk     inland    unknown 
	#	   763         45       1577        369 
	inmigrant	<- subset(inmigrant, select=c(RCCS_studyid, date, INMIG_LOC))	
	inmigrant[, date:= hivc.db.Date2numeric(date)]	
	setnames(inmigrant, colnames(inmigrant), gsub('DATE','INMIGRATE_DATE',gsub('RCCS_STUDYID','RID',toupper(colnames(inmigrant)))))
	
	
	#	load transmission events
	if(is.null(infile.inference))
	{
		infile.inference					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_withmetadata.rda"
		infile.inference					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_withmetadata.rda"
		infile.inference					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min10_withmetadata.rda"
		infile.inference					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min20_withmetadata.rda"
		infile.inference					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min50_withmetadata.rda"		
	}
	outfile.base			<- gsub('171122','180522',gsub('_withmetadata.rda','',infile.inference))
	load(infile.inference)
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	#stopifnot(293==nrow(rtpdm))	
	#	add new variables
	rtpdm[, AGEDIFF:= rtpdm[, FEMALE_BIRTHDATE-MALE_BIRTHDATE]]	
	set(rtpdm, NULL, 'MALE_SEX', 'M')
	set(rtpdm, NULL, 'FEMALE_SEX', 'F')		
	#	add immigrant: to compute if individual was inmigrant within 1 year of first conc pos, we need date first conc pos	
	tmp			<- subset(rtpdm, select=c(PAIRID, MALE_RID, FEMALE_RID, VISIT_FIRSTCONCPOS))
	tmp2		<- unique(subset(rd, select=c(RID, VISIT, DATE)))
	setnames(tmp2, colnames(tmp2), gsub('MALE_VISIT','VISIT_FIRSTCONCPOS',paste0('MALE_',colnames(tmp2))))
	tmp			<- merge(tmp, tmp2, by=c('MALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	setnames(tmp2, colnames(tmp2), gsub('MALE_','FEMALE_', colnames(tmp2)))
	tmp			<- merge(tmp, tmp2, by=c('FEMALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	set(tmp, tmp[, which(is.na(MALE_DATE))], 'MALE_DATE', -1)
	set(tmp, tmp[, which(is.na(FEMALE_DATE))], 'FEMALE_DATE', -1)
	#	argh there are entries where both have missing date in rd ..
	tmp2		<- subset(tmp, MALE_DATE==-1 & FEMALE_DATE==-1, select=c(MALE_RID, FEMALE_RID, VISIT_FIRSTCONCPOS, PAIRID))
	tmp2		<- merge(tmp2, subset(rtpdm, select=c(MALE_RID, FEMALE_RID, VISIT_FIRSTCONCPOS, PAIRID, MALE_DATE, FEMALE_DATE)), by=c('MALE_RID','FEMALE_RID','VISIT_FIRSTCONCPOS','PAIRID'))
	tmp			<- rbind(subset(tmp, !(MALE_DATE==-1 & FEMALE_DATE==-1)), tmp2)	
	set(tmp, NULL, 'DATE_FIRSTCONCPOS', tmp[, pmax(MALE_DATE, FEMALE_DATE)])
	stopifnot( nrow(subset(tmp, DATE_FIRSTCONCPOS==-1))==0 )
	set(tmp, NULL, c('VISIT_FIRSTCONCPOS','MALE_DATE','FEMALE_DATE'), NULL)
	#	find inmigrants among transmission pairs
	tmp2		<- copy(inmigrant)
	setnames(tmp2, colnames(tmp2), paste0('MALE_',colnames(tmp2)))	
	tmp			<- merge(tmp, tmp2, by='MALE_RID', all.x=1)
	setnames(tmp2, colnames(tmp2), gsub('MALE_','FEMALE_', colnames(tmp2)))
	tmp			<- merge(tmp, tmp2, by='FEMALE_RID', all.x=1)
	set(tmp, tmp[, which(is.na(MALE_INMIGRATE_DATE))], 'MALE_INMIGRATE_DATE', -1)
	set(tmp, tmp[, which(is.na(FEMALE_INMIGRATE_DATE))], 'FEMALE_INMIGRATE_DATE', -1)
	#	don t consider inmigrations after transmission occurred
	set(tmp, tmp[, which(DATE_FIRSTCONCPOS<MALE_INMIGRATE_DATE)], c('MALE_INMIG_PLACE','MALE_INMIG_ADMIN0','MALE_INMIG_ADMIN1','MALE_INMIG_ADMIN2'), NA_character_)
	set(tmp, tmp[, which(DATE_FIRSTCONCPOS<MALE_INMIGRATE_DATE)], 'MALE_INMIGRATE_DATE', -1)
	set(tmp, tmp[, which(DATE_FIRSTCONCPOS<FEMALE_INMIGRATE_DATE)], c('FEMALE_INMIG_PLACE','FEMALE_INMIG_ADMIN0','FEMALE_INMIG_ADMIN1','FEMALE_INMIG_ADMIN2'), NA_character_)
	set(tmp, tmp[, which(DATE_FIRSTCONCPOS<FEMALE_INMIGRATE_DATE)], 'FEMALE_INMIGRATE_DATE', -1)
	tmp[, MALE_TIME_INMIGRATED:= DATE_FIRSTCONCPOS-MALE_INMIGRATE_DATE]
	tmp[, FEMALE_TIME_INMIGRATED:= DATE_FIRSTCONCPOS-FEMALE_INMIGRATE_DATE]	
	#	select most recent inmigration events
	tmp[, DUMMY:= seq_len(nrow(tmp))]
	tmp2	<- tmp[, list(DUMMY=DUMMY[min(FEMALE_TIME_INMIGRATED)==FEMALE_TIME_INMIGRATED]), by=c('PAIRID')]
	tmp		<- merge(tmp, tmp2, by=c('PAIRID','DUMMY'))
	tmp2	<- tmp[, list(DUMMY=DUMMY[min(MALE_TIME_INMIGRATED)==MALE_TIME_INMIGRATED]), by=c('PAIRID')]
	tmp		<- merge(tmp, tmp2, by=c('PAIRID','DUMMY'))
	tmp[, DUMMY:=NULL]
	#	define inmigrated in last year, and in last two years
	tmp[, MALE_INMIGRATE_1YR:= as.character(as.integer(MALE_TIME_INMIGRATED<=1))]
	tmp[, MALE_INMIGRATE_2YR:= as.character(as.integer(MALE_TIME_INMIGRATED<=2))]
	tmp[, MALE_INMIGRATE_3YR:= as.character(as.integer(MALE_TIME_INMIGRATED<=3))]
	tmp[, FEMALE_INMIGRATE_1YR:= as.character(as.integer(FEMALE_TIME_INMIGRATED<=1))]
	tmp[, FEMALE_INMIGRATE_2YR:= as.character(as.integer(FEMALE_TIME_INMIGRATED<=2))]	
	tmp[, FEMALE_INMIGRATE_3YR:= as.character(as.integer(FEMALE_TIME_INMIGRATED<=3))]
	tmp		<- subset(tmp, select=c(PAIRID, DATE_FIRSTCONCPOS, MALE_INMIG_LOC, MALE_INMIGRATE_1YR, MALE_INMIGRATE_2YR, MALE_INMIGRATE_3YR, FEMALE_INMIG_LOC, FEMALE_INMIGRATE_1YR, FEMALE_INMIGRATE_2YR, FEMALE_INMIGRATE_3YR))
	#	make compound variable
	tmp		<- melt(tmp, id.vars=c('PAIRID','DATE_FIRSTCONCPOS','MALE_INMIG_LOC','FEMALE_INMIG_LOC'))
	set(tmp, tmp[, which(value=='0')], 'value', 'resident')
	tmp2	<- tmp[, which(value=='1' & grepl('^FEMALE_',variable))]
	set(tmp, tmp2, 'value', tmp[tmp2, paste0('inmigrant_from_',FEMALE_INMIG_LOC)])
	tmp2	<- tmp[, which(value=='1' & grepl('^MALE_',variable))]
	set(tmp, tmp2, 'value', tmp[tmp2, paste0('inmigrant_from_',MALE_INMIG_LOC)])
	tmp		<- dcast.data.table(tmp, PAIRID+DATE_FIRSTCONCPOS~variable, value.var='value')
	
	#	tmp[, table(MALE_INMIGRATE_2YR,  MALE_INMIGRATE_3YR)]
	#	only one more migration event when we allow for 3 years
	#	ignore this
	set(tmp, NULL, c('FEMALE_INMIGRATE_3YR','MALE_INMIGRATE_3YR'), NULL)
	
	rtpdm		<- merge(rtpdm, tmp, by=c('PAIRID'))
	rtpdm[, table(MALE_INMIGRATE_2YR, FEMALE_INMIGRATE_2YR)]
	#                           	FEMALE_INMIGRATE_2YR
	#	MALE_INMIGRATE_2YR          inmigrant_from_external inmigrant_from_fisherfolk inmigrant_from_inland resident
	# inmigrant_from_external                         5                         0                     5        8
	# inmigrant_from_fisherfolk                       0                         0                     0        1
	# inmigrant_from_inland                           3                         0                    10       22
	# resident                                       22                         6                    37      174
	
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
	
	#	save
	save(rtpdm, rtr2, zm, ds, file=gsub('180522','181006',paste0(outfile.base,'_phylogeography_data_with_inmigrants.rda')))
}

RakaiFull.phylogeography.181006.flows.wrapper<- function()
{

	#
	#	make data set
	#
	if(0)
	{
		infiles		<- c("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min10_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min20_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min50_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf50_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf55_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf65_withmetadata.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf70_withmetadata.rda"
		)
		for(ii in seq_along(infiles))
		{
			infile.inference	<- infiles[[ii]]
			RakaiFull.phylogeography.180521.gender.mobility.data(infile.inference)		
		}
	}
	
	#
	#	core inference
	#
	if(1)
	{
		indir		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run'
		indir		<- '/work/or105/Gates_2014/Rakai'
		opt									<- list()
		opt$adjust.sequencing.bias			<- 1
		opt$adjust.participation.bias		<- 1
		opt$exclude.onART.from.denominator	<- 1
		opt$set.missing.migloc.to.inland	<- 0
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland
		infiles		<- c( "RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min10_phylogeography_data_with_inmigrants.rda"
						, "RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min20_phylogeography_data_with_inmigrants.rda"
						, "RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min50_phylogeography_data_with_inmigrants.rda"
						, "RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf50_phylogeography_data_with_inmigrants.rda"
						, "RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf55_phylogeography_data_with_inmigrants.rda"
						, "RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda"
						, "RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf65_phylogeography_data_with_inmigrants.rda"
						, "RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf70_phylogeography_data_with_inmigrants.rda"
						, "todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda"
						)
		for(ii in seq_along(infiles))
		{
			infile.inference	<- file.path(indir, infiles[[ii]])
			RakaiFull.phylogeography.181006.gender.mobility.core.inference(infile.inference=infile.inference, opt=opt)		
		}
		#
		opt									<- list()
		opt$adjust.sequencing.bias			<- 0
		opt$adjust.participation.bias		<- 1
		opt$exclude.onART.from.denominator	<- 1
		opt$set.missing.migloc.to.inland	<- 0
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland
		infile.inference	<- file.path(indir, "todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda")
		RakaiFull.phylogeography.181006.gender.mobility.core.inference(infile.inference=infile.inference, opt=opt)		
		opt									<- list()
		opt$adjust.sequencing.bias			<- 1
		opt$adjust.participation.bias		<- 0
		opt$exclude.onART.from.denominator	<- 1
		opt$set.missing.migloc.to.inland	<- 0
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland
		infile.inference	<- file.path(indir, "todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda")
		RakaiFull.phylogeography.181006.gender.mobility.core.inference(infile.inference=infile.inference, opt=opt)		
		opt									<- list()
		opt$adjust.sequencing.bias			<- 0
		opt$adjust.participation.bias		<- 0
		opt$exclude.onART.from.denominator	<- 1
		opt$set.missing.migloc.to.inland	<- 0
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland
		infile.inference	<- file.path(indir, "todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda")
		RakaiFull.phylogeography.181006.gender.mobility.core.inference(infile.inference=infile.inference, opt=opt)		
		opt									<- list()
		opt$adjust.sequencing.bias			<- 1
		opt$adjust.participation.bias		<- 1
		opt$exclude.onART.from.denominator	<- 0
		opt$set.missing.migloc.to.inland	<- 0
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland
		infile.inference	<- file.path(indir, "todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda")
		RakaiFull.phylogeography.181006.gender.mobility.core.inference(infile.inference=infile.inference, opt=opt)		
		opt									<- list()
		opt$adjust.sequencing.bias			<- 1
		opt$adjust.participation.bias		<- 1
		opt$exclude.onART.from.denominator	<- 1
		opt$set.missing.migloc.to.inland	<- 1
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland
		infile.inference	<- file.path(indir, "todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda")
		RakaiFull.phylogeography.181006.gender.mobility.core.inference(infile.inference=infile.inference, opt=opt)		
	}
		
	if(0)
	{
		#
		#	extract target variables of interest
		#
		infiles		<- c(#"~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_inference.rda"
				#, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_inference_seqExclART.rda"
				#, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min10_phylogeography_core_inference.rda"
				#, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min20_phylogeography_core_inference.rda"
				#, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min50_phylogeography_core_inference.rda"
				"~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf50_phylogeography_core_inference.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf55_phylogeography_core_inference.rda"
		#, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_core_inference.rda"
		#, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf65_phylogeography_core_inference.rda"
		#, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf70_phylogeography_core_inference.rda"
		)
		for(ii in seq_along(infiles))
		{
			infile.inference	<- infiles[[ii]]
			#RakaiFull.phylogeography.180618.flows.fishinlandgender(infile.inference)
			RakaiFull.phylogeography.180618.flows.fishinlandmigrant(infile.inference)
			RakaiFull.phylogeography.180618.flows.fishinlandmigrationgender.netflows(infile.inference)
		}
	}
}

RakaiFull.phylogeography.181006.gender.mobility.core.inference<- function(infile.inference=NULL, opt=NULL)
{
	require(data.table)	
	require(Boom)	
	#require(gtools)	
	
	logistic<- function(x) 1/(1+e(-x))
		
	if(is.null(opt))
	{
		opt									<- list()
		opt$adjust.sequencing.bias			<- 1
		opt$adjust.participation.bias		<- 1
		opt$exclude.onART.from.denominator	<- 1
		opt$set.missing.migloc.to.inland	<- 0
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland		
	}
	if(!is.null(infile.inference))
	{
		indir				<- dirname(infile.inference)
	}
	if(is.null(infile.inference))
	{
		indir				<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run"
		infile.inference	<- file.path(indir,"todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda")		
	}
	
	outfile.base	<- gsub('170704','181006',gsub('data_with_inmigrants.rda','',infile.inference))
	load(infile.inference)
		
	#	select which inmigration data to work with
	setnames(rtr2, c('TR_INMIGRATE_2YR','REC_INMIGRATE_2YR'), c('TR_INMIGRATE','REC_INMIGRATE'))
	#	double check that TR_INMIGRATE does not have any unknown origin
	tmp	<- rtr2[, which(TR_INMIGRATE=='inmigrant_from_unknown')]
	if(opt$set.missing.migloc.to.inland)
		set(rtr2, tmp, 'TR_INMIGRATE', 'inmigrant_from_inland')
	if(opt$set.missing.migloc.to.fishing)
		set(rtr2, tmp, 'TR_INMIGRATE', 'inmigrant_from_fisherfolk')
	#	select variables that we need, define missing ones
	rtr2	<- subset(rtr2, select=c(PAIRID, TR_RID, REC_RID, TR_COMM_NUM_A, REC_COMM_NUM_A, TR_SEX, REC_SEX, TR_INMIGRATE, REC_INMIGRATE, TR_BIRTHDATE, REC_BIRTHDATE, TR_COMM_TYPE, REC_COMM_TYPE))
	set(rtr2, NULL, 'TR_COMM_NUM_A', rtr2[, as.character(TR_COMM_NUM_A)])
	set(rtr2, NULL, 'REC_COMM_NUM_A', rtr2[, as.character(REC_COMM_NUM_A)])		
	#	impute missing birth dates	
	set(rtr2, rtr2[, which(is.na(TR_BIRTHDATE))], 'TR_BIRTHDATE', rtr2[, mean(TR_BIRTHDATE, na.rm=TRUE)])
	set(rtr2, rtr2[, which(is.na(REC_BIRTHDATE))], 'REC_BIRTHDATE', rtr2[, mean(REC_BIRTHDATE, na.rm=TRUE)])
	#	set age at midpoint of study period
	rtr2[, TR_AGE_AT_MID:= 2013.25 - TR_BIRTHDATE]
	rtr2[, REC_AGE_AT_MID:= 2013.25 - REC_BIRTHDATE]
	#	stratify age
	rtr2[, TR_AGE_AT_MID_C:= as.character(cut(TR_AGE_AT_MID, breaks=c(10,25,35,65), labels=c('15-24','25-34','35+'), right=FALSE))]
	rtr2[, REC_AGE_AT_MID_C:= as.character(cut(REC_AGE_AT_MID, breaks=c(10,25,35,65), labels=c('15-24','25-34','35+'), right=FALSE))]
	stopifnot( nrow(subset(rtr2, is.na(TR_AGE_AT_MID_C)))==0 )
	stopifnot( nrow(subset(rtr2, is.na(REC_AGE_AT_MID_C)))==0 )
	#	delete what we don t need
	set(rtr2, NULL, c('TR_BIRTHDATE','REC_BIRTHDATE','TR_AGE_AT_MID','REC_AGE_AT_MID'), NULL)
	
	#
	#	Bayesian model sampling prior: 
	#	define prior for sampling probabilities from previously fitted participation and sequencing models	
	#	sampling probs differ by: age community and sex
	#
	infile.participation	<- file.path(indir,"participation_differences_180322_logisticmodels.rda")	
	infile.sequencing		<- file.path(indir,"sequencing_differences_180322_logisticmodels.rda")
	if(opt$exclude.onART.from.denominator)
		infile.sequencing	<- file.path(indir,"sequencing_differences_180322_exclART_logisticmodels.rda")
	load(infile.participation)
	mp1 <- mp2 <- mp4 <- NULL
	#	binarize covariates
	dc	<- copy(rtr2)
	dc[, TR_AGE_YOUNG:= as.integer(TR_AGE_AT_MID_C=='15-24')]
	dc[, TR_AGE_MID:= as.integer(TR_AGE_AT_MID_C=='25-34')]
	dc[, TR_MALE:= as.integer(TR_SEX=='M')]
	dc[, TR_COMM_TYPE_F:= as.integer(substr(TR_COMM_NUM_A,1,1)=='f')]
	dc[, TR_COMM_TYPE_T:= as.integer(substr(TR_COMM_NUM_A,1,1)=='t')]
	dc[, TR_INMIGRANT:= as.integer(grepl('inmigrant',TR_INMIGRATE))]	
	dc[, REC_AGE_YOUNG:= as.integer(REC_AGE_AT_MID_C=='15-24')]
	dc[, REC_AGE_MID:= as.integer(REC_AGE_AT_MID_C=='25-34')]
	dc[, REC_MALE:= as.integer(REC_SEX=='M')]
	dc[, REC_COMM_TYPE_F:= as.integer(substr(REC_COMM_NUM_A,1,1)=='f')]
	dc[, REC_COMM_TYPE_T:= as.integer(substr(REC_COMM_NUM_A,1,1)=='t')]
	dc[, REC_INMIGRANT:= as.integer(grepl('inmigrant',REC_INMIGRATE))]
	#	define community number to match STAN output for participation model
	tmp		<- unique(subset(dg, select=c(COMM_NUM_A,COMM_NUM_B)))
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc		<- merge(dc, tmp, by='TR_COMM_NUM_A')
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc		<- merge(dc, tmp, by='REC_COMM_NUM_A')
	#	define community number to match STAN output for sequencing model
	load(infile.sequencing)
	ms2 <- ms3 <- ms4 <- NULL
	tmp		<- unique(subset(dg, select=c(COMM_NUM_A,COMM_NUM_B)))
	setnames(tmp, 'COMM_NUM_B', 'COMM_NUM_B2')
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc		<- merge(dc, tmp, by='TR_COMM_NUM_A')
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc		<- merge(dc, tmp, by='REC_COMM_NUM_A')
	#	extract Monte Carlo samples from best WAIC participation and best WAIC sequencing models	
	mps			<- extract.samples(mp3)
	mss			<- extract.samples(ms1)
	mc.it		<- 5e2
	set.seed(42L)
	gc()
	#	for every source-recipient pair, get mc.it samples from their sampling probabilities
	dc		<- dc[, {
				#	get Monte Carlo samples from posterior distribution of logit(participation probs)
				p.tr.part	<- with(mps, a + comm[, TR_COMM_NUM_B] + trading*TR_COMM_TYPE_T + fishing*TR_COMM_TYPE_F + 
								inmigrant*TR_INMIGRANT + inmigrant_young*TR_INMIGRANT*TR_AGE_YOUNG + 
								male*TR_MALE + 
								young_male*TR_AGE_YOUNG*TR_MALE + young_female*TR_AGE_YOUNG*(1-TR_MALE) +
								midage*TR_AGE_MID)
				p.rec.part	<- with(mps, a + comm[, REC_COMM_NUM_B] + trading*REC_COMM_TYPE_T + fishing*REC_COMM_TYPE_F + 
								inmigrant*REC_INMIGRANT + inmigrant_young*REC_INMIGRANT*REC_AGE_YOUNG + 
								male*REC_MALE + 
								young_male*REC_AGE_YOUNG*REC_MALE + young_female*REC_AGE_YOUNG*(1-REC_MALE) +
								midage*REC_AGE_MID)
				#	get Monte Carlo samples from posterior distribution of logit(sequencing probs)			
				p.tr.seq	<- with(mss, a + comm[, TR_COMM_NUM_B2] + trading*TR_COMM_TYPE_T + fishing*TR_COMM_TYPE_F + 
								inmigrant*TR_INMIGRANT + male*TR_MALE + young*TR_AGE_YOUNG + midage*TR_AGE_MID)
				p.rec.seq	<- with(mss, a + comm[, REC_COMM_NUM_B2] + trading*REC_COMM_TYPE_T + fishing*REC_COMM_TYPE_F + 
								inmigrant*REC_INMIGRANT + male*REC_MALE + young*REC_AGE_YOUNG + midage*REC_AGE_MID)
				#	sensitivity analyses			
				if(!opt$adjust.participation.bias)
				{
					p.tr.part <- p.rec.part <- rep(Inf, length(p.tr.part))
				}					
				if(!opt$adjust.sequencing.bias)
				{
					p.tr.seq <- p.rec.seq <- rep(Inf, length(p.tr.seq))
				}				
				#	sample mc.it many, and calculate product of sequencing and participation probs
				mc.idx		<- sample.int(length(p.tr.part), mc.it, replace=TRUE)
				#tmp			<- log1p(exp(-p.tr.part[mc.idx])) + log1p(exp(-p.rec.part[mc.idx]))
				#mc.idx		<- sample.int(length(p.tr.seq), mc.it, replace=TRUE)
				#tmp			<- tmp + log1p(exp(-p.tr.seq[mc.idx])) + log1p(exp(-p.rec.seq[mc.idx]))
				#list( IT=seq_len(mc.it), S=as.vector(exp(-tmp)) )
				list( IT=seq_len(mc.it), 
					  TR_P_PART= as.vector(logistic(p.tr.part[mc.idx])),
					  REC_P_PART= as.vector(logistic(p.rec.part[mc.idx])),
					  TR_P_SEQ= as.vector(logistic(p.tr.seq[mc.idx])),
					  REC_P_SEQ= as.vector(logistic(p.rec.seq[mc.idx]))	)				
		}, by=c('PAIRID','TR_COMM_NUM_A','REC_COMM_NUM_A','TR_SEX','REC_SEX','TR_INMIGRATE','REC_INMIGRATE')]
	dc[, S:= TR_P_PART*REC_P_PART*TR_P_SEQ*REC_P_SEQ]
	#
	#	for each combination of transmitter, recipient covariates:
	#	determine best fitting beta distribution	
	tmp	<- dc[, {
				#	mean and sd
				mu		<- mean(S)
				sd		<- sd(S)				
				#	convert https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
				alpha	<- mu*mu*(1-mu)/(sd*sd)-mu				
				list( BETA_OK=sd*sd <= mu*(1-mu), S_MU=mu, S_SD=sd, S_ALPHA=alpha, S_BETA= alpha*(1/mu-1))
			}, by=c('TR_COMM_NUM_A','REC_COMM_NUM_A','TR_SEX','REC_SEX','TR_INMIGRATE','REC_INMIGRATE')]
	stopifnot(tmp[, all(BETA_OK)])
	
	#	calculate observed number of transmissions
	dc	<- dc[, list( TR_OBS=length(unique(PAIRID))), by=c('TR_COMM_NUM_A','REC_COMM_NUM_A','TR_SEX','REC_SEX','TR_INMIGRATE','REC_INMIGRATE')]
	dc	<- merge(dc, tmp, by=c('TR_COMM_NUM_A','REC_COMM_NUM_A','TR_SEX','REC_SEX','TR_INMIGRATE','REC_INMIGRATE'))	
	set(dc, NULL, 'TR_COMM_NUM_A', dc[, as.character(TR_COMM_NUM_A)])
	set(dc, NULL, 'REC_COMM_NUM_A', dc[, as.character(REC_COMM_NUM_A)])	
	setkey(dc, TR_COMM_NUM_A, REC_COMM_NUM_A, TR_SEX, TR_INMIGRATE, REC_INMIGRATE)
	dc[, COUNT_ID:= seq_len(nrow(dc))]
	#	156 non-zero counts
	
	# calculate mean sampling probability
	dc[, S:=   S_ALPHA/(S_ALPHA+S_BETA)]	
	
	# set up mcmc objects
	mc				<- list()
	mc$n			<- 1e7
	mc$pars			<- list() 	
	mc$pars$LAMBDA	<- matrix(NA_real_, ncol=length(unique(dc$COUNT_ID)), nrow=1)		#prior for proportions
	mc$pars$S		<- matrix(NA_real_, ncol=length(unique(dc$COUNT_ID)), nrow=mc$n)	#sampling probabilities (product)
	mc$pars$Z		<- matrix(NA_integer_, ncol=length(unique(dc$COUNT_ID)), nrow=mc$n) #augmented data
	mc$pars$NU		<- NA_real_															#prior for N
	mc$pars$N		<- matrix(NA_integer_, ncol=1, nrow=mc$n)							#total number of counts on augmented data
	mc$pars$PI		<- matrix(NA_real_, ncol=length(unique(dc$COUNT_ID)), nrow=mc$n)	#proportions	
	mc$it.info		<- data.table(	IT= seq.int(1,mc$n),
			PAR_ID= rep(NA_integer_, mc$n),
			BLOCK= rep(NA_character_, mc$n),
			MHRATIO= rep(NA_real_, mc$n),
			ACCEPT=rep(NA_integer_, mc$n), 
			LOG_LKL=rep(NA_real_, mc$n),
			LOG_PRIOR=rep(NA_real_, mc$n))
	
	#
	# define helper functions
	lddirichlet_vector	<- function(x, nu){
		ans	<- sum((nu - 1) * log(x)) + sum(lgamma(nu)) - lgamma(sum(nu))
		stopifnot(is.finite(ans))
		ans
	}
	
	# initialise MCMC
	setkey(dc, COUNT_ID)
	mc$verbose			<- 0L
	mc$curr.it			<- 1L
	mc$seed				<- 42L
	set(mc$it.info, mc$curr.it, 'BLOCK', 'INIT')
	set(mc$it.info, mc$curr.it, 'PAR_ID', 0L)
	set.seed( mc$seed )
	#	prior lambda: use the Berger objective prior with minimal loss compared to marginal Beta reference prior
	#	(https://projecteuclid.org/euclid.ba/1422556416)	
	mc$pars$LAMBDA[1,]	<- 0.8/nrow(dc)
	# 	sampling: set to draw from prior
	tmp					<- dc[, list(S= rbeta(1L, S_ALPHA, S_BETA)), by=c('COUNT_ID')]	
	setkey(tmp, COUNT_ID)			
	mc$pars$S[1,]		<- tmp$S
	#	augmented data: proposal draw under sampling probability
	mc$pars$Z[1,]		<- dc[, TR_OBS + rnbinom(nrow(dc), TR_OBS, mc$pars$S[1,])]	
	#	prior nu: set Poisson rate to the expected augmented counts, under average sampling probability
	mc$pars$NU			<- sum(dc$TR_OBS) / mean(dc$S)
	#	total count: that s just the sum of Z
	mc$pars$N[1,]		<- sum(mc$pars$Z[1,])
	#	proportions: draw from full conditional
	mc$pars$PI[1,]		<- rdirichlet(1, mc$pars$Z[1,] + mc$pars$LAMBDA[1,])			
	#	store log likelihood
	tmp	<- sum( dbinom(dc$TR_OBS, size=mc$pars$Z[1,], prob=mc$pars$S[1,], log=TRUE) ) +
			dmultinom(mc$pars$Z[1,], size=mc$pars$N[1,], prob=mc$pars$PI[1,], log=TRUE)
	set(mc$it.info, 1L, 'LOG_LKL', tmp)
	# 	store log prior		
	tmp	<- dpois(mc$pars$N[1,], lambda=mc$pars$NU, log=TRUE) +
			lddirichlet_vector(mc$pars$PI[1,], nu=mc$pars$LAMBDA[1,]) +
			sum( dbeta( mc$pars$S[1,], dc[, S_ALPHA], dc[, S_BETA], log=TRUE ) )	
	set(mc$it.info, 1L, 'LOG_PRIOR', tmp)
	
	#	
	# run mcmc
	options(warn=2)
	mc$sweep			<- ncol(mc$pars$Z) + 1L
	for(i in 1L:(mc$n-1L))
	{
		mc$curr.it		<- i		
		# determine source-recipient combination that will be updated in this iteration
		update.count	<- (i-1L) %% mc$sweep + 1L
		# update S, Z, N for the source-recipient combination 'update.count'
		if(update.count<mc$sweep)
		{
			#	propose  
			S.prop					<- mc$pars$S[mc$curr.it,]
			S.prop[update.count]	<- dc[COUNT_ID==update.count, rbeta(1L, S_ALPHA, S_BETA)]		
			Z.prop					<- mc$pars$Z[mc$curr.it,]
			Z.prop[update.count]	<- dc[COUNT_ID==update.count, TR_OBS + rnbinom(1, TR_OBS, S.prop[update.count])]
			N.prop					<- sum(Z.prop)
			#	calculate MH ratio
			log.prop.ratio			<- sum(dnbinom(mc$pars$Z[mc$curr.it,update.count], size=dc$TR_OBS[update.count], prob=mc$pars$S[mc$curr.it,update.count], log=TRUE)) - 
										sum(dnbinom(Z.prop[update.count], size=dc$TR_OBS[update.count], prob=S.prop[update.count], log=TRUE))
			log.fc					<- sum(dbinom(dc$TR_OBS[update.count], size=mc$pars$Z[mc$curr.it,update.count], prob=mc$pars$S[mc$curr.it,update.count], log=TRUE)) +
										dmultinom(mc$pars$Z[mc$curr.it,], prob=mc$pars$PI[mc$curr.it,], log=TRUE) +
										dpois(mc$pars$N[mc$curr.it,], lambda=mc$pars$NU, log=TRUE)
			log.fc.prop				<- sum(dbinom(dc$TR_OBS[update.count], size=Z.prop[update.count], prob=S.prop[update.count], log=TRUE)) +
										dmultinom(Z.prop, prob=mc$pars$PI[mc$curr.it,], log=TRUE) +
										dpois(N.prop, lambda=mc$pars$NU, log=TRUE)
			log.mh.ratio			<- log.fc.prop - log.fc + log.prop.ratio
			mh.ratio				<- min(1,exp(log.mh.ratio))	
			#	update
			mc$curr.it				<- mc$curr.it+1L
			set(mc$it.info, mc$curr.it, 'BLOCK', 'S-Z-N')
			set(mc$it.info, mc$curr.it, 'PAR_ID', update.count)
			set(mc$it.info, mc$curr.it, 'MHRATIO', mh.ratio)
			set(mc$it.info, mc$curr.it, 'ACCEPT', as.integer(runif(1) < mh.ratio))
			mc$pars$PI[mc$curr.it,]	<- mc$pars$PI[mc$curr.it-1L,]
			if(mc$verbose & mc$it.info[mc$curr.it, ACCEPT])
			{
				print(paste0('it ',mc$curr.it,' ACCEPT S-Z-N block ',update.count))
			}
			if(mc$it.info[mc$curr.it, ACCEPT])
			{
				mc$pars$Z[mc$curr.it,]	<- Z.prop
				mc$pars$N[mc$curr.it,]	<- N.prop		
				mc$pars$S[mc$curr.it,]	<- S.prop				
			}
			if(mc$it.info[mc$curr.it, !ACCEPT])
			{
				mc$pars$Z[mc$curr.it,]	<- mc$pars$Z[mc$curr.it-1L,]
				mc$pars$N[mc$curr.it,]	<- mc$pars$N[mc$curr.it-1L,]				
				mc$pars$S[mc$curr.it,]	<- mc$pars$S[mc$curr.it-1L,]
			}	
		}
		# update PI
		if(update.count==mc$sweep)
		{
			#	propose 
			PI.prop					<- rdirichlet(1L, nu= mc$pars$Z[mc$curr.it,]+mc$pars$LAMBDA[1,])
			#	this is the full conditional of PI given S, N, Z
			#	always accept
			#	update
			mc$curr.it				<- mc$curr.it+1L
			set(mc$it.info, mc$curr.it, 'BLOCK', 'PI')
			set(mc$it.info, mc$curr.it, 'PAR_ID', NA_integer_)
			set(mc$it.info, mc$curr.it, 'MHRATIO', 1L)
			set(mc$it.info, mc$curr.it, 'ACCEPT', 1L)			
			mc$pars$S[mc$curr.it,]	<- mc$pars$S[mc$curr.it-1L,]
			mc$pars$Z[mc$curr.it,]	<- mc$pars$Z[mc$curr.it-1L,]
			mc$pars$N[mc$curr.it,]	<- mc$pars$N[mc$curr.it-1L,]				
			mc$pars$PI[mc$curr.it,]	<- PI.prop			
		}		
		# store log likelihood
		tmp	<- sum( dbinom(dc$TR_OBS, size=mc$pars$Z[mc$curr.it,], prob=mc$pars$S[mc$curr.it,], log=TRUE) ) +
				dmultinom(mc$pars$Z[mc$curr.it,], size=mc$pars$N[mc$curr.it,], prob=mc$pars$PI[mc$curr.it,], log=TRUE)
		set(mc$it.info, mc$curr.it, 'LOG_LKL', tmp)
		# store log prior	
		tmp	<- dpois(mc$pars$N[mc$curr.it,], lambda=mc$pars$NU, log=TRUE) +
				lddirichlet_vector(mc$pars$PI[1,], nu=mc$pars$LAMBDA[1,]) +
				sum(dbeta(mc$pars$S[mc$curr.it,], dc[, S_ALPHA], dc[, S_BETA], log=TRUE))
		set(mc$it.info, mc$curr.it, 'LOG_PRIOR', tmp)		
	}		
	save(zm, rtpdm, rtr2, ds, dc, mc, file=paste0(outfile.base,'core_inference_mcmc_',paste0(opt, collapse=''),'.rda'))
}


RakaiFull.phylogeography.180928.gender.mobility.core.inference<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(Boom)	
	
	indir								<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run"	
	if(is.null(infile.inference))
	{		
		infile.inference	<- file.path(indir,"todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda")		
	}
	
	outfile.base	<- gsub('180522','180928',gsub('data_with_inmigrants.rda','',infile.inference))
	load(infile.inference)
	
	infile			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda"
	outfile.base	<- gsub('data_with_inmigrants.rda','',infile)
	load(infile)
	
	#	select which inmigration data to work with
	setnames(rtr2, c('TR_INMIGRATE_2YR','REC_INMIGRATE_2YR'), c('TR_INMIGRATE','REC_INMIGRATE'))
	
	#	sum transmission events by community and gender
	dc	<- rtr2[, list(TR_OBS=length(PAIRID)), by=c('TR_COMM_NUM_A','REC_COMM_NUM_A','TR_SEX','REC_SEX')]
	set(dc, NULL, 'TR_COMM_NUM_A', dc[, as.character(TR_COMM_NUM_A)])
	set(dc, NULL, 'REC_COMM_NUM_A', dc[, as.character(REC_COMM_NUM_A)])	
	setkey(dc, TR_COMM_NUM_A, REC_COMM_NUM_A, TR_SEX)
	dc[, COUNT_ID:= seq_len(nrow(dc))]
	#	101 non-zero counts
	
	#
	#	sampling prior: 
	#	define Beta prior for sampling probabilities as Empirical Bayes prior from cohort data (all alpha and betas)
	#	currently sampling probs differ by community and sex
	#
	tmp	<- subset(ds, select=c(COMM_NUM_A, SEX, P_PART_ALPHA, P_PART_BETA, P_SEQ_ALPHA, P_SEQ_BETA))
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc	<- merge(dc, tmp, by=c('TR_COMM_NUM_A','TR_SEX'))
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc	<- merge(dc, tmp, by=c('REC_COMM_NUM_A','REC_SEX'))
	
	
	# for simplicity, calculate mean for now
	dc[, S:=   TR_P_PART_ALPHA/(TR_P_PART_ALPHA+TR_P_PART_BETA) * 
					TR_P_SEQ_ALPHA/(TR_P_SEQ_ALPHA+TR_P_SEQ_BETA) * 
					REC_P_PART_ALPHA/(REC_P_PART_ALPHA+REC_P_PART_BETA) * 
					REC_P_SEQ_ALPHA/(REC_P_SEQ_ALPHA+REC_P_SEQ_BETA)]	
	
	# set up mcmc objects
	mc				<- list()
	mc$n			<- 1e6
	mc$pars			<- list() 	
	mc$pars$LAMBDA	<- matrix(NA_real_, ncol=length(unique(dc$COUNT_ID)), nrow=1)		#prior for proportions
	mc$pars$S		<- matrix(NA_real_, ncol=length(unique(dc$COUNT_ID)), nrow=mc$n)	#sampling probabilities (product)
	mc$pars$S_DTL	<- data.table()														#sampling probabilities (detailed terms that make up product)
	mc$pars$Z		<- matrix(NA_integer_, ncol=length(unique(dc$COUNT_ID)), nrow=mc$n) #augmented data
	mc$pars$NU		<- NA_real_															#prior for N
	mc$pars$N		<- matrix(NA_integer_, ncol=1, nrow=mc$n)							#total number of counts on augmented data
	mc$pars$PI		<- matrix(NA_real_, ncol=length(unique(dc$COUNT_ID)), nrow=mc$n)	#proportions	
	mc$it.info		<- data.table(	IT= seq.int(1,mc$n),
									PAR_ID= rep(NA_integer_, mc$n),
									BLOCK= rep(NA_character_, mc$n),
									MHRATIO= rep(NA_real_, mc$n),
									ACCEPT=rep(NA_integer_, mc$n), 
									LOG_LKL=rep(NA_real_, mc$n),
									LOG_PRIOR=rep(NA_real_, mc$n))
	
	#
	# define helper functions
	lddirichlet_vector	<- function(x, nu){
		ans	<- sum((nu - 1) * log(x)) + sum(lgamma(nu)) - lgamma(sum(nu))
		stopifnot(is.finite(ans))
		ans
	}
							
	# initialise MCMC
	setkey(dc, COUNT_ID)
	mc$verbose			<- 0L
	mc$curr.it			<- 1L
	mc$seed				<- 42L
	set(mc$it.info, mc$curr.it, 'BLOCK', 'INIT')
	set(mc$it.info, mc$curr.it, 'PAR_ID', 0L)
	set.seed( mc$seed )
	#	prior lambda: use the Berger objective prior with minimal loss compared to marginal Beta reference prior
	#	(https://projecteuclid.org/euclid.ba/1422556416)	
	mc$pars$LAMBDA[1,]	<- 0.8/nrow(dc)
	# 	sampling: set to draw from prior
	tmp					<- dc[, list(	TR_PART_P= rbeta(1L, TR_P_PART_ALPHA, TR_P_PART_BETA),
										TR_SEQ_P= rbeta(1L, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA),
										REC_PART_P= rbeta(1L, REC_P_PART_ALPHA, REC_P_PART_BETA),
										REC_SEQ_P= rbeta(1L, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA)
										),			
								by=c('COUNT_ID','TR_P_PART_ALPHA','TR_P_PART_BETA','REC_P_PART_ALPHA','REC_P_PART_BETA','TR_P_SEQ_ALPHA','TR_P_SEQ_BETA','REC_P_SEQ_ALPHA','REC_P_SEQ_BETA')]
	tmp[, TR_PART_LOGD:= dbeta(TR_PART_P, TR_P_PART_ALPHA, TR_P_PART_BETA, log=TRUE)]
	tmp[, TR_SEQ_LOGD:= dbeta(TR_SEQ_P, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA, log=TRUE)]
	tmp[, REC_PART_LOGD:= dbeta(REC_PART_P, REC_P_PART_ALPHA, REC_P_PART_BETA, log=TRUE)]
	tmp[, REC_SEQ_LOGD:= dbeta(REC_SEQ_P, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA, log=TRUE)]
	mc$pars$S_DTL	<- copy(tmp)
	tmp				<- mc$pars$S_DTL[, list(S=TR_PART_P*TR_SEQ_P*REC_PART_P*REC_SEQ_P), by='COUNT_ID']
	setkey(tmp, COUNT_ID)			
	mc$pars$S[1,]		<- tmp$S
	#	augmented data: proposal draw under sampling probability
	mc$pars$Z[1,]		<- dc[, TR_OBS + rnbinom(nrow(dc), TR_OBS, mc$pars$S[1,])]	
	#	prior nu: set Poisson rate to the expected augmented counts, under average sampling probability
	mc$pars$NU			<- sum(dc$TR_OBS) / mean(dc$S)
	#	total count: that s just the sum of Z
	mc$pars$N[1,]		<- sum(mc$pars$Z[1,])
	#	proportions: draw from full conditional
	mc$pars$PI[1,]		<- rdirichlet(1, mc$pars$Z[1,] + mc$pars$LAMBDA[1,])			
	#	store log likelihood
	tmp	<- sum( dbinom(dc$TR_OBS, size=mc$pars$Z[1,], prob=mc$pars$S[1,], log=TRUE) ) +
			dmultinom(mc$pars$Z[1,], size=mc$pars$N[1,], prob=mc$pars$PI[1,], log=TRUE)
	set(mc$it.info, 1L, 'LOG_LKL', tmp)
	# 	store log prior		
	tmp	<- dpois(mc$pars$N[1,], lambda=mc$pars$NU, log=TRUE) +
			lddirichlet_vector(mc$pars$PI[1,], nu=mc$pars$LAMBDA[1,]) +
			mc$pars$S_DTL[,sum(TR_PART_LOGD) + sum(TR_SEQ_LOGD) + sum(REC_PART_LOGD) + sum(REC_SEQ_LOGD)]	
	set(mc$it.info, 1L, 'LOG_PRIOR', tmp)
		
	#	
	# run mcmc
	options(warn=2)
	mc$sweep			<- ncol(mc$pars$Z) + 1L
	for(i in 1L:(mc$n-1L))
	{
		mc$curr.it		<- i		
		# determine source-recipient combination that will be updated in this iteration
		update.count	<- (i-1L) %% mc$sweep + 1L		
		# update S, Z, N for the source-recipient combination 'update.count'
		if(update.count<mc$sweep)
		{
			#	propose  
			S.prop					<- mc$pars$S[mc$curr.it,]
			S_DTL.prop				<- copy(mc$pars$S_DTL)
			tmp						<- dc[COUNT_ID==update.count, list(	TR_PART_P= rbeta(1L, TR_P_PART_ALPHA, TR_P_PART_BETA),
																		TR_SEQ_P= rbeta(1L, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA),
																		REC_PART_P= rbeta(1L, REC_P_PART_ALPHA, REC_P_PART_BETA),
																		REC_SEQ_P= rbeta(1L, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA)),			
											by=c('COUNT_ID','TR_P_PART_ALPHA','TR_P_PART_BETA','REC_P_PART_ALPHA','REC_P_PART_BETA','TR_P_SEQ_ALPHA','TR_P_SEQ_BETA','REC_P_SEQ_ALPHA','REC_P_SEQ_BETA')]
			tmp[, TR_PART_LOGD:= dbeta(TR_PART_P, TR_P_PART_ALPHA, TR_P_PART_BETA, log=TRUE)]
			tmp[, TR_SEQ_LOGD:= dbeta(TR_SEQ_P, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA, log=TRUE)]
			tmp[, REC_PART_LOGD:= dbeta(REC_PART_P, REC_P_PART_ALPHA, REC_P_PART_BETA, log=TRUE)]
			tmp[, REC_SEQ_LOGD:= dbeta(REC_SEQ_P, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA, log=TRUE)]
			S_DTL.prop				<- rbind(subset(S_DTL.prop, COUNT_ID!=update.count), tmp)
			setkey(S_DTL.prop, COUNT_ID)
			S.prop[update.count]	<- tmp[, TR_PART_P*TR_SEQ_P*REC_PART_P*REC_SEQ_P]
			Z.prop					<- mc$pars$Z[mc$curr.it,]
			Z.prop[update.count]	<- dc[COUNT_ID==update.count, TR_OBS + rnbinom(1, TR_OBS, S.prop)]
			N.prop					<- sum(Z.prop)			
			#	calculate MH ratio
			log.prop.ratio			<- sum(dnbinom(mc$pars$Z[mc$curr.it,update.count], size=dc$TR_OBS[update.count], prob=mc$pars$S[mc$curr.it,update.count], log=TRUE)) - 
										sum(dnbinom(Z.prop[update.count], size=dc$TR_OBS[update.count], prob=S.prop[update.count], log=TRUE))
			log.fc					<- sum(dbinom(dc$TR_OBS[update.count], size=mc$pars$Z[mc$curr.it,update.count], prob=mc$pars$S[mc$curr.it,update.count], log=TRUE)) +
										dmultinom(mc$pars$Z[mc$curr.it,], prob=mc$pars$PI[mc$curr.it,], log=TRUE) +
										dpois(mc$pars$N[mc$curr.it,], lambda=mc$pars$NU, log=TRUE)
			log.fc.prop				<- sum(dbinom(dc$TR_OBS[update.count], size=Z.prop[update.count], prob=S.prop[update.count], log=TRUE)) +
										dmultinom(Z.prop, prob=mc$pars$PI[mc$curr.it,], log=TRUE) +
										dpois(N.prop, lambda=mc$pars$NU, log=TRUE)
			log.mh.ratio			<- log.fc.prop - log.fc + log.prop.ratio
			mh.ratio				<- min(1,exp(log.mh.ratio))	
			#	update
			mc$curr.it				<- mc$curr.it+1L
			set(mc$it.info, mc$curr.it, 'BLOCK', 'S-Z-N')
			set(mc$it.info, mc$curr.it, 'PAR_ID', update.count)
			set(mc$it.info, mc$curr.it, 'MHRATIO', mh.ratio)
			set(mc$it.info, mc$curr.it, 'ACCEPT', as.integer(runif(1) < mh.ratio))
			mc$pars$PI[mc$curr.it,]	<- mc$pars$PI[mc$curr.it-1L,]
			if(mc$verbose & mc$it.info[mc$curr.it, ACCEPT])
			{
				print(paste0('it ',mc$curr.it,' ACCEPT S-Z-N block ',update.count))
			}
			if(mc$it.info[mc$curr.it, ACCEPT])
			{
				mc$pars$Z[mc$curr.it,]	<- Z.prop
				mc$pars$N[mc$curr.it,]	<- N.prop		
				mc$pars$S[mc$curr.it,]	<- S.prop
				mc$pars$S_DTL			<- copy(S_DTL.prop)
			}
			if(mc$it.info[mc$curr.it, !ACCEPT])
			{
				mc$pars$Z[mc$curr.it,]	<- mc$pars$Z[mc$curr.it-1L,]
				mc$pars$N[mc$curr.it,]	<- mc$pars$N[mc$curr.it-1L,]				
				mc$pars$S[mc$curr.it,]	<- mc$pars$S[mc$curr.it-1L,]
			}	
		}
		# update PI
		if(update.count==mc$sweep)
		{
			#	propose 
			PI.prop					<- rdirichlet(1L, nu= mc$pars$Z[mc$curr.it,]+mc$pars$LAMBDA[1,])
			#	this is the full conditional of PI given S, N, Z
			#	always accept
			#	update
			mc$curr.it				<- mc$curr.it+1L
			set(mc$it.info, mc$curr.it, 'BLOCK', 'PI')
			set(mc$it.info, mc$curr.it, 'PAR_ID', update.count)
			set(mc$it.info, mc$curr.it, 'MHRATIO', 1L)
			set(mc$it.info, mc$curr.it, 'ACCEPT', 1L)			
			mc$pars$S[mc$curr.it,]	<- mc$pars$S[mc$curr.it-1L,]
			mc$pars$Z[mc$curr.it,]	<- mc$pars$Z[mc$curr.it-1L,]
			mc$pars$N[mc$curr.it,]	<- mc$pars$N[mc$curr.it-1L,]				
			mc$pars$PI[mc$curr.it,]	<- PI.prop			
		}		
		# store log likelihood
		tmp	<- sum( dbinom(dc$TR_OBS, size=mc$pars$Z[mc$curr.it,], prob=mc$pars$S[mc$curr.it,], log=TRUE) ) +
				dmultinom(mc$pars$Z[mc$curr.it,], size=mc$pars$N[mc$curr.it,], prob=mc$pars$PI[mc$curr.it,], log=TRUE)
		set(mc$it.info, mc$curr.it, 'LOG_LKL', tmp)
		# store log prior	
		tmp	<- dpois(mc$pars$N[mc$curr.it,], lambda=mc$pars$NU, log=TRUE) +
				lddirichlet_vector(mc$pars$PI[1,], nu=mc$pars$LAMBDA[1,]) +
				mc$pars$S_DTL[,sum(TR_PART_LOGD) + sum(TR_SEQ_LOGD) + sum(REC_PART_LOGD) + sum(REC_SEQ_LOGD)]
		set(mc$it.info, mc$curr.it, 'LOG_PRIOR', tmp)		
	}
	save(zm, rtpdm, rtr2, ds, dc, mc, file=paste0(outfile.base,'core_inference_SNPIZ_mcmcEachCount.rda'))
}

RakaiFull.phylogeography.180928.gender.mobility.core.inference.old<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(Boom)	
	
	indir								<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run"	
	if(is.null(infile.inference))
	{		
		infile.inference	<- file.path(indir,"todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda")		
	}
	
	outfile.base	<- gsub('180522','180928',gsub('data_with_inmigrants.rda','',infile.inference))
	load(infile.inference)
	
	infile			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda"
	outfile.base	<- gsub('data_with_inmigrants.rda','',infile)
	load(infile)
	
	#	select which inmigration data to work with
	setnames(rtr2, c('TR_INMIGRATE_2YR','REC_INMIGRATE_2YR'), c('TR_INMIGRATE','REC_INMIGRATE'))
	
	#	sum transmission events by community and gender
	dc	<- rtr2[, list(TR_OBS=length(PAIRID)), by=c('TR_COMM_NUM_A','REC_COMM_NUM_A','TR_SEX','REC_SEX')]
	set(dc, NULL, 'TR_COMM_NUM_A', dc[, as.character(TR_COMM_NUM_A)])
	set(dc, NULL, 'REC_COMM_NUM_A', dc[, as.character(REC_COMM_NUM_A)])	
	setkey(dc, TR_COMM_NUM_A, REC_COMM_NUM_A, TR_SEX)
	dc[, COUNT_ID:= seq_len(nrow(dc))]
	#	101 non-zero counts
	
	#
	#	sampling prior: 
	#	define Beta prior for sampling probabilities as Empirical Bayes prior from cohort data (all alpha and betas)
	#	currently sampling probs differ by community and sex
	#
	tmp	<- subset(ds, select=c(COMM_NUM_A, SEX, P_PART_ALPHA, P_PART_BETA, P_SEQ_ALPHA, P_SEQ_BETA))
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc	<- merge(dc, tmp, by=c('TR_COMM_NUM_A','TR_SEX'))
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc	<- merge(dc, tmp, by=c('REC_COMM_NUM_A','REC_SEX'))
	
	
	# for simplicity, calculate mean for now
	dc[, S:=   TR_P_PART_ALPHA/(TR_P_PART_ALPHA+TR_P_PART_BETA) * 
					TR_P_SEQ_ALPHA/(TR_P_SEQ_ALPHA+TR_P_SEQ_BETA) * 
					REC_P_PART_ALPHA/(REC_P_PART_ALPHA+REC_P_PART_BETA) * 
					REC_P_SEQ_ALPHA/(REC_P_SEQ_ALPHA+REC_P_SEQ_BETA)]	
	
	# set up mcmc objects
	mc				<- list()
	mc$n			<- 1e6
	mc$pars			<- list() 	
	mc$pars$LAMBDA	<- matrix(NA_real_, ncol=length(unique(dc$COUNT_ID)), nrow=1)		#prior for proportions
	mc$pars$S		<- matrix(NA_real_, ncol=length(unique(dc$COUNT_ID)), nrow=mc$n)	#sampling probabilities (product)
	mc$pars$S_DTL	<- data.table()														#sampling probabilities (detailed terms that make up product)
	mc$pars$Z		<- matrix(NA_integer_, ncol=length(unique(dc$COUNT_ID)), nrow=mc$n) #augmented data
	mc$pars$NU		<- NA_real_															#prior for N
	mc$pars$N		<- matrix(NA_integer_, ncol=1, nrow=mc$n)							#total number of counts on augmented data
	mc$pars$PI		<- matrix(NA_real_, ncol=length(unique(dc$COUNT_ID)), nrow=mc$n)	#proportions	
	mc$it.info		<- data.table(	IT= seq.int(1,mc$n),
			BLOCK= rep(NA_integer_, mc$n),
			MHRATIO= rep(NA_real_, mc$n),
			ACCEPT=rep(NA_integer_, mc$n), 
			LOG_LKL=rep(NA_real_, mc$n),
			LOG_PRIOR=rep(NA_real_, mc$n))
	
	#
	# define helper functions
	lddirichlet_vector	<- function(x, nu){
		ans	<- sum((nu - 1) * log(x)) + sum(lgamma(nu)) - lgamma(sum(nu))
		stopifnot(is.finite(ans))
		ans
	}
	
	# initialise MCMC
	setkey(dc, COUNT_ID)
	mc$verbose			<- 0L
	mc$curr.it			<- 1L
	mc$seed				<- 42L
	set(mc$it.info, mc$curr.it, 'BLOCK', 0L)
	set.seed( mc$seed )
	#	prior lambda: use the Berger objective prior with minimal loss compared to marginal Beta reference prior
	#	(https://projecteuclid.org/euclid.ba/1422556416)	
	mc$pars$LAMBDA[1,]	<- 0.8/nrow(dc)
	# 	sampling: set to draw from prior
	tmp					<- dc[, list(	TR_PART_P= rbeta(1L, TR_P_PART_ALPHA, TR_P_PART_BETA),
					TR_SEQ_P= rbeta(1L, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA),
					REC_PART_P= rbeta(1L, REC_P_PART_ALPHA, REC_P_PART_BETA),
					REC_SEQ_P= rbeta(1L, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA)
			),			
			by=c('COUNT_ID','TR_P_PART_ALPHA','TR_P_PART_BETA','REC_P_PART_ALPHA','REC_P_PART_BETA','TR_P_SEQ_ALPHA','TR_P_SEQ_BETA','REC_P_SEQ_ALPHA','REC_P_SEQ_BETA')]
	tmp[, TR_PART_LOGD:= dbeta(TR_PART_P, TR_P_PART_ALPHA, TR_P_PART_BETA, log=TRUE)]
	tmp[, TR_SEQ_LOGD:= dbeta(TR_SEQ_P, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA, log=TRUE)]
	tmp[, REC_PART_LOGD:= dbeta(REC_PART_P, REC_P_PART_ALPHA, REC_P_PART_BETA, log=TRUE)]
	tmp[, REC_SEQ_LOGD:= dbeta(REC_SEQ_P, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA, log=TRUE)]
	mc$pars$S_DTL	<- copy(tmp)
	tmp				<- mc$pars$S_DTL[, list(S=TR_PART_P*TR_SEQ_P*REC_PART_P*REC_SEQ_P), by='COUNT_ID']
	setkey(tmp, COUNT_ID)			
	mc$pars$S[1,]		<- tmp$S
	#	augmented data: proposal draw under sampling probability
	mc$pars$Z[1,]		<- dc[, TR_OBS + rnbinom(nrow(dc), TR_OBS, mc$pars$S[1,])]	
	#	prior nu: set Poisson rate to the expected augmented counts, under average sampling probability
	mc$pars$NU			<- sum(dc$TR_OBS) / mean(dc$S)
	#	total count: that s just the sum of Z
	mc$pars$N[1,]		<- sum(mc$pars$Z[1,])
	#	proportions: draw from full conditional
	mc$pars$PI[1,]		<- rdirichlet(1, mc$pars$Z[1,] + mc$pars$LAMBDA[1,])			
	#	store log likelihood
	tmp	<- sum( dbinom(dc$TR_OBS, size=mc$pars$Z[1,], prob=mc$pars$S[1,], log=TRUE) ) +
			dmultinom(mc$pars$Z[1,], size=mc$pars$N[1,], prob=mc$pars$PI[1,], log=TRUE)
	set(mc$it.info, 1L, 'LOG_LKL', tmp)
	# 	store log prior		
	tmp	<- dpois(mc$pars$N[1,], lambda=mc$pars$NU, log=TRUE) +
			lddirichlet_vector(mc$pars$PI[1,], nu=mc$pars$LAMBDA[1,]) +
			mc$pars$S_DTL[,sum(TR_PART_LOGD) + sum(TR_SEQ_LOGD) + sum(REC_PART_LOGD) + sum(REC_SEQ_LOGD)]	
	set(mc$it.info, 1L, 'LOG_PRIOR', tmp)
	
	#	
	# run mcmc
	options(warn=2)
	mc$blocks		<- ncol(mc$pars$Z)
	for(i in 1L:(mc$n-1L))
	{
		mc$curr.it		<- i		
		# determine source-recipient combination that will be updated in this iteration
		update.count	<- (i-1L) %% mc$blocks + 1L
		
		# update S, Z, N, PI for the source-recipient combination 'update.count'
		#	propose	
		S.prop					<- mc$pars$S[mc$curr.it,]
		S_DTL.prop				<- copy(mc$pars$S_DTL)
		tmp						<- dc[COUNT_ID==update.count, list(	TR_PART_P= rbeta(1L, TR_P_PART_ALPHA, TR_P_PART_BETA),
						TR_SEQ_P= rbeta(1L, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA),
						REC_PART_P= rbeta(1L, REC_P_PART_ALPHA, REC_P_PART_BETA),
						REC_SEQ_P= rbeta(1L, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA)),			
				by=c('COUNT_ID','TR_P_PART_ALPHA','TR_P_PART_BETA','REC_P_PART_ALPHA','REC_P_PART_BETA','TR_P_SEQ_ALPHA','TR_P_SEQ_BETA','REC_P_SEQ_ALPHA','REC_P_SEQ_BETA')]
		tmp[, TR_PART_LOGD:= dbeta(TR_PART_P, TR_P_PART_ALPHA, TR_P_PART_BETA, log=TRUE)]
		tmp[, TR_SEQ_LOGD:= dbeta(TR_SEQ_P, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA, log=TRUE)]
		tmp[, REC_PART_LOGD:= dbeta(REC_PART_P, REC_P_PART_ALPHA, REC_P_PART_BETA, log=TRUE)]
		tmp[, REC_SEQ_LOGD:= dbeta(REC_SEQ_P, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA, log=TRUE)]
		S_DTL.prop				<- rbind(subset(S_DTL.prop, COUNT_ID!=update.count), tmp)
		setkey(S_DTL.prop, COUNT_ID)
		S.prop[update.count]	<- tmp[, TR_PART_P*TR_SEQ_P*REC_PART_P*REC_SEQ_P]
		Z.prop					<- mc$pars$Z[mc$curr.it,]
		Z.prop[update.count]	<- dc[COUNT_ID==update.count, TR_OBS + rnbinom(1, TR_OBS, S.prop)]
		N.prop					<- sum(Z.prop)
		PI.prop					<- mc$pars$PI[mc$curr.it,]
		tmp						<- c(Z.prop[update.count]+mc$pars$LAMBDA[1,update.count], sum(Z.prop + mc$pars$LAMBDA[1,]))
		PI.prop[update.count]	<- rbeta(1, tmp[1], tmp[2]-tmp[1])							
		
		#	calculate MH ratio
		log.prop.ratio	<- sum(dnbinom(mc$pars$Z[mc$curr.it,update.count], size=dc$TR_OBS[update.count], prob=mc$pars$S[mc$curr.it,update.count], log=TRUE)) - 
				sum(dnbinom(Z.prop[update.count], size=dc$TR_OBS[update.count], prob=S.prop[update.count], log=TRUE))
		log.fc			<- sum(dbinom(dc$TR_OBS[update.count], size=mc$pars$Z[mc$curr.it,update.count], prob=mc$pars$S[mc$curr.it,update.count], log=TRUE)) +								
				dpois(mc$pars$N[mc$curr.it,], lambda=mc$pars$NU, log=TRUE)
		log.fc.prop		<- sum(dbinom(dc$TR_OBS[update.count], size=Z.prop[update.count], prob=S.prop[update.count], log=TRUE)) +								
				dpois(N.prop, lambda=mc$pars$NU, log=TRUE)
		log.mh.ratio	<- log.fc.prop - log.fc + log.prop.ratio
		mh.ratio		<- min(1,exp(log.mh.ratio))
		
		#	update
		mc$curr.it				<- mc$curr.it+1L
		set(mc$it.info, mc$curr.it, 'BLOCK', update.count)
		set(mc$it.info, mc$curr.it, 'MHRATIO', mh.ratio)
		set(mc$it.info, mc$curr.it, 'ACCEPT', as.integer(runif(1) < mh.ratio))							
		if(mc$verbose & mc$it.info[mc$curr.it, ACCEPT])
		{
			print(paste0('it ',mc$curr.it,' ACCEPT block ',update.count))
		}
		if(mc$it.info[mc$curr.it, ACCEPT])
		{
			mc$pars$Z[mc$curr.it,]	<- Z.prop
			mc$pars$N[mc$curr.it,]	<- N.prop		
			mc$pars$PI[mc$curr.it,]	<- PI.prop
			mc$pars$S[mc$curr.it,]	<- S.prop
			mc$pars$S_DTL			<- copy(S_DTL.prop)
		}
		if(mc$it.info[mc$curr.it, !ACCEPT])
		{
			mc$pars$Z[mc$curr.it,]	<- mc$pars$Z[mc$curr.it-1L,]
			mc$pars$N[mc$curr.it,]	<- mc$pars$N[mc$curr.it-1L,]
			mc$pars$PI[mc$curr.it,]	<- mc$pars$PI[mc$curr.it-1L,]
			mc$pars$S[mc$curr.it,]	<- mc$pars$S[mc$curr.it-1L,]
		}			
		
		# store log likelihood
		tmp	<- sum( dbinom(dc$TR_OBS, size=mc$pars$Z[mc$curr.it,], prob=mc$pars$S[mc$curr.it,], log=TRUE) ) +
				dmultinom(mc$pars$Z[mc$curr.it,], size=mc$pars$N[mc$curr.it,], prob=mc$pars$PI[mc$curr.it,], log=TRUE)
		set(mc$it.info, mc$curr.it, 'LOG_LKL', tmp)
		# store log prior	
		tmp	<- dpois(mc$pars$N[mc$curr.it,], lambda=mc$pars$NU, log=TRUE) +
				lddirichlet_vector(mc$pars$PI[1,], nu=mc$pars$LAMBDA[1,]) +
				mc$pars$S_DTL[,sum(TR_PART_LOGD) + sum(TR_SEQ_LOGD) + sum(REC_PART_LOGD) + sum(REC_SEQ_LOGD)]
		set(mc$it.info, mc$curr.it, 'LOG_PRIOR', tmp)		
	}
	save(zm, rtpdm, rtr2, ds, dc, mc, file=paste0(outfile.base,'core_inference_SNPIZ_mcmcEachCount.rda'))
}

RakaiFull.phylogeography.181006.mcmc.assess<- function(infile.inference=NULL)
{
	require(coda)
	if(is.null(infile.inference))
	{
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_core_inference_mcmc.rda"
	}			
	outfile.base		<- gsub('.rda$','',gsub('_core_inference','',infile.inference))
	load(infile.inference)
	
	#	thin MCMC output
	tmp			<- seq.int(2, nrow(mc$pars$Z), mc$sweep)
	mc$pars$S	<- mc$pars$S[tmp,,drop=FALSE]
	mc$pars$Z	<- mc$pars$Z[tmp,,drop=FALSE]
	mc$pars$PI	<- mc$pars$PI[tmp,,drop=FALSE]
	mc$pars$N	<- mc$pars$N[tmp,,drop=FALSE]
	gc()
	colnames(mc$pars$S)	<- paste0('S-',1:ncol(mc$pars$S))
	colnames(mc$pars$Z)	<- paste0('Z-',1:ncol(mc$pars$Z))
	colnames(mc$pars$PI)<- paste0('PI-',1:ncol(mc$pars$PI))
	colnames(mc$pars$N)	<- 'N'
	#	convert to coda
	mcc			<- mcmc( cbind( mc$pars$N, mc$pars$Z, mc$pars$PI, mc$pars$S) )
	mcc.eff		<- effectiveSize(mcc)
	#	plot trace and autocorrelations for variable with worst 6 effective sizes
	tmp			<- mcc.eff[mcc.eff>0] 
	tmp			<- names(tmp[sort(tmp, index.return=TRUE)$ix[1:6]])	
	pdf(file=paste0(outfile.base,'_6worstchains_trace.pdf'), w=10, h=7)
	plot( mcc[,tmp] )
	dev.off()
	pdf(file=paste0(outfile.base,'_6worstchains_autocor.pdf'), w=5, h=7)
	autocorr.plot( mcc[,tmp] )
	dev.off()
	#	print effective sample size to file
	mcc.eff		<- data.table(PAR=names(mcc.eff), NEFF=as.numeric(mcc.eff))
	setkey(mcc.eff, NEFF)
}

RakaiFull.phylogeography.181006.flows.fishinland<- function()
{
	require(data.table)	
	require(Hmisc)	
	
	if(is.null(infile.inference))
	{
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_core_inference_mcmc.rda"
	}			
	outfile.base		<- gsub('.rda$','',gsub('_core_inference','',infile.inference))
	load(infile.inference)
	
	#	thin MCMC output
	tmp			<- seq.int(2, nrow(mc$pars$Z), mc$sweep)
	mc$pars$S	<- mc$pars$S[tmp,,drop=FALSE]
	mc$pars$Z	<- mc$pars$Z[tmp,,drop=FALSE]
	mc$pars$PI	<- mc$pars$PI[tmp,,drop=FALSE]
	mc$pars$N	<- mc$pars$N[tmp,,drop=FALSE]
	gc()
	colnames(mc$pars$S)	<- paste0('S-',1:ncol(mc$pars$S))
	colnames(mc$pars$Z)	<- paste0('Z-',1:ncol(mc$pars$Z))
	colnames(mc$pars$PI)<- paste0('PI-',1:ncol(mc$pars$PI))
	colnames(mc$pars$N)	<- 'N'
	
	#
	#	prepare data.table of proportions
	dc[, REC_COMM_TYPE:= as.character(factor(substr(REC_COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fisherfolk','inland')))]
	dc[, TR_COMM_TYPE:= as.character(factor(substr(TR_COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fisherfolk','inland')))]
	mcpi	<- as.data.table(mc$pars$PI)
	mcpi[, IT:= seq.int(2, by=100, length.out=nrow(mcpi))]
	mcpi	<- melt(mcpi, id.vars='IT', variable.name='COUNT_ID', value.name='PI')
	set(mcpi, NULL, 'COUNT_ID', mcpi[, gsub('([A-Z]+)-([0-9]+)','\\2',COUNT_ID)])
	set(mcpi, NULL, 'COUNT_ID', mcpi[, as.integer(COUNT_ID)])
	tmp		<- subset(dc, select=c(REC_COMM_NUM_A, REC_SEX, REC_COMM_TYPE, TR_COMM_NUM_A, TR_SEX, TR_COMM_TYPE, TR_OBS, COUNT_ID))
	mcpi	<- merge(tmp, mcpi, by='COUNT_ID')
	
	#
	#	calculate overall transmission flows fishing-inland
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')
	
	#
	#	geography who infects whom matrix  between fisherfolk and others
	#	adjusted P
	z		<- mcpi[, list(PI=sum(PI)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','IT')]
	z		<- z[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by=c('REC_COMM_TYPE','TR_COMM_TYPE')]
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_COMM_TYPE, TR_COMM_TYPE)
	z[, STAT:='joint']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- copy(z)	

	#
	#	WAIFM
	#
	groups	<- c('inland','fisherfolk')
	z		<- lapply(groups, function(group)
			{												
				z		<- subset(mcpi, TR_COMM_TYPE==group)
				z		<- z[, list(PI=sum(PI)), by=c('REC_COMM_TYPE','IT')]
				z		<- z[, list(REC_COMM_TYPE=REC_COMM_TYPE, PI=PI/sum(PI)), by=c('IT')]				
				z		<- z[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by='REC_COMM_TYPE']
				z[, TR_COMM_TYPE:= group]
				z
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE)
	z[, STAT:='waifm']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans,z)	

	#
	#	sources
	#
	groups	<- c('inland','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				z		<- subset(mcpi , REC_COMM_TYPE==group)
				z		<- z[, list(PI=sum(PI)), by=c('TR_COMM_TYPE','IT')]
				z		<- z[, list(TR_COMM_TYPE=TR_COMM_TYPE, PI=PI/sum(PI)), by=c('IT')]
				z		<- z[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by='TR_COMM_TYPE']
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

	save(ans, file=paste0(outfile.base,'_fishing_inland_results.rda'))
}

RakaiFull.phylogeography.181006.flows.fishinlandgender<- function()
{	
	require(data.table)	
	require(Hmisc)	
	
	if(is.null(infile.inference))
	{
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_core_inference_mcmc.rda"
	}			
	outfile.base		<- gsub('.rda$','',gsub('_core_inference','',infile.inference))
	load(infile.inference)
	
	#	thin MCMC output
	tmp			<- seq.int(2, nrow(mc$pars$Z), mc$sweep)
	mc$pars$S	<- mc$pars$S[tmp,,drop=FALSE]
	mc$pars$Z	<- mc$pars$Z[tmp,,drop=FALSE]
	mc$pars$PI	<- mc$pars$PI[tmp,,drop=FALSE]
	mc$pars$N	<- mc$pars$N[tmp,,drop=FALSE]
	gc()
	colnames(mc$pars$S)	<- paste0('S-',1:ncol(mc$pars$S))
	colnames(mc$pars$Z)	<- paste0('Z-',1:ncol(mc$pars$Z))
	colnames(mc$pars$PI)<- paste0('PI-',1:ncol(mc$pars$PI))
	colnames(mc$pars$N)	<- 'N'
	
	#
	#	prepare data.table of proportions
	dc[, REC_COMM_TYPE:= as.character(factor(substr(REC_COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fisherfolk','inland')))]
	dc[, TR_COMM_TYPE:= as.character(factor(substr(TR_COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fisherfolk','inland')))]
	mcpi	<- as.data.table(mc$pars$PI)
	mcpi[, IT:= seq.int(2, by=100, length.out=nrow(mcpi))]
	mcpi	<- melt(mcpi, id.vars='IT', variable.name='COUNT_ID', value.name='PI')
	set(mcpi, NULL, 'COUNT_ID', mcpi[, gsub('([A-Z]+)-([0-9]+)','\\2',COUNT_ID)])
	set(mcpi, NULL, 'COUNT_ID', mcpi[, as.integer(COUNT_ID)])
	tmp		<- subset(dc, select=c(REC_COMM_NUM_A, REC_SEX, REC_COMM_TYPE, TR_COMM_NUM_A, TR_SEX, TR_COMM_TYPE, TR_OBS, COUNT_ID))
	mcpi	<- merge(tmp, mcpi, by='COUNT_ID')
	
	#
	#	calculate overall transmission flows fishing-inland
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')
	
	#
	#	geography who infects whom matrix  between fisherfolk and others
	#	adjusted P
	z		<- mcpi[, list(PI=sum(PI)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_SEX','REC_SEX','IT')]
	z		<- z[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by=c('TR_COMM_TYPE','TR_SEX','REC_COMM_TYPE','REC_SEX')]
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_SEX+REC_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_SEX, REC_SEX )
	z[, STAT:='joint']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- copy(z)	
	z		<- NULL	
	gc()
	
	#
	#	WAIFM
	#
	groups	<- data.table(TR_COMM_TYPE=c('inland','inland','fisherfolk','fisherfolk'), TR_SEX=c('M','F','M','F'))
	z		<- lapply(1:nrow(groups), function(ii)
			{												
				z		<- subset(mcpi, TR_COMM_TYPE==groups$TR_COMM_TYPE[ii] & TR_SEX==groups$TR_SEX[ii])
				z		<- z[, list(PI=sum(PI)), by=c('REC_COMM_TYPE','REC_SEX','IT')]
				z		<- z[, list(REC_COMM_TYPE=REC_COMM_TYPE, REC_SEX=REC_SEX, PI=PI/sum(PI)), by=c('IT')]				
				z		<- z[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by=c('REC_COMM_TYPE','REC_SEX')]
				z[, TR_COMM_TYPE:= groups$TR_COMM_TYPE[ii] ]
				z[, TR_SEX:= groups$TR_SEX[ii] ]
				z
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_SEX+REC_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, TR_SEX, REC_COMM_TYPE, REC_SEX)
	z[, STAT:='waifm']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans,z)	
	z		<- NULL	
	gc()
	
	#
	#	sources
	#
	groups	<- data.table(REC_COMM_TYPE=c('inland','inland','fisherfolk','fisherfolk'), REC_SEX=c('M','F','M','F'))
	z		<- lapply(1:nrow(groups), function(ii)
			{				
				z		<- subset(mcpi, REC_COMM_TYPE==groups$REC_COMM_TYPE[ii] & REC_SEX==groups$REC_SEX[ii])
				z		<- z[, list(PI=sum(PI)), by=c('TR_COMM_TYPE','TR_SEX','IT')]
				z		<- z[, list(TR_COMM_TYPE=TR_COMM_TYPE, TR_SEX=TR_SEX, PI=PI/sum(PI)), by=c('IT')]
				z		<- z[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by=c('TR_COMM_TYPE','TR_SEX')]
				z[, REC_COMM_TYPE:= groups$REC_COMM_TYPE[ii] ]
				z[, REC_SEX:= groups$REC_SEX[ii] ]
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_SEX+REC_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_COMM_TYPE, REC_SEX, TR_COMM_TYPE, TR_SEX)
	z[, STAT:='sources']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans,z)	
	z		<- NULL	
	gc()
		
	save(ans, file=paste0(outfile.base,'_fishing_inland_results.rda'))
}

RakaiFull.phylogeography.181006.flows.netflows<- function(infile.inference=NULL)
{
	require(data.table)	
	require(Hmisc)	
	
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_core_inference_mcmc.rda"	
	outfile.base		<- gsub('.rda$','',gsub('_core_inference','',infile.inference))
	load(infile.inference)
	
	#	thin MCMC output
	tmp			<- seq.int(2, nrow(mc$pars$Z), mc$sweep)
	mc$pars$S	<- mc$pars$S[tmp,,drop=FALSE]
	mc$pars$Z	<- mc$pars$Z[tmp,,drop=FALSE]
	mc$pars$PI	<- mc$pars$PI[tmp,,drop=FALSE]
	mc$pars$N	<- mc$pars$N[tmp,,drop=FALSE]
	gc()
	colnames(mc$pars$S)	<- paste0('S-',1:ncol(mc$pars$S))
	colnames(mc$pars$Z)	<- paste0('Z-',1:ncol(mc$pars$Z))
	colnames(mc$pars$PI)<- paste0('PI-',1:ncol(mc$pars$PI))
	colnames(mc$pars$N)	<- 'N'
	
	#
	#	prepare data.table of proportions
	dc[, REC_COMM_TYPE:= as.character(factor(substr(REC_COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fisherfolk','inland')))]
	dc[, TR_COMM_TYPE:= as.character(factor(substr(TR_COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fisherfolk','inland')))]
	mcpi	<- as.data.table(mc$pars$PI)
	mcpi[, IT:= seq.int(2, by=100, length.out=nrow(mcpi))]
	mcpi	<- melt(mcpi, id.vars='IT', variable.name='COUNT_ID', value.name='PI')
	set(mcpi, NULL, 'COUNT_ID', mcpi[, gsub('([A-Z]+)-([0-9]+)','\\2',COUNT_ID)])
	set(mcpi, NULL, 'COUNT_ID', mcpi[, as.integer(COUNT_ID)])
	tmp		<- subset(dc, select=c(REC_COMM_NUM_A, REC_SEX, REC_COMM_TYPE, REC_INMIGRATE, TR_COMM_NUM_A, TR_SEX, TR_INMIGRATE, TR_COMM_TYPE, TR_OBS, COUNT_ID))
	mcpi	<- merge(tmp, mcpi, by='COUNT_ID')
	
	
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')
	
	#
	#	geography flows inland->fish / flows fish->inland
	#	by gender and migration status
	z		<- copy(mcpi)
	#	set simple migration status
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	set(z, z[, which(TR_INMIGRATE=='resident')], 'TR_INMIGRATE', 'resident/outmigrant')	
	z		<- z[, list(	TR_OBS=sum(TR_OBS), PI=sum(PI)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_SEX','REC_SEX','TR_INMIGRATE','IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',TR_INMIGRATE,' ',TR_SEX,' ',REC_COMM_TYPE,' ',REC_SEX)]
	z		<- dcast.data.table(z, IT~FLOW, value.var='PI')	
	set(z, NULL, c('inland resident/outmigrant F inland M','external inmigrant_from_external F inland M','external inmigrant_from_external F fisherfolk M','fisherfolk resident/outmigrant F fisherfolk M'), NULL)
	set(z, NULL, c('inland resident/outmigrant M inland F','external inmigrant_from_external M inland F','external inmigrant_from_external M fisherfolk F','fisherfolk resident/outmigrant M fisherfolk F'), NULL)
	z[ , inlanddivfisherfolk_M:= z[['inland resident/outmigrant M fisherfolk F']] / z[['fisherfolk resident/outmigrant M inland F']] ]
	z[ , inlanddivfisherfolk_F:= z[['inland resident/outmigrant F fisherfolk M']] / z[['fisherfolk resident/outmigrant F inland M']] ]				
	z		<- melt(z, id.vars='IT', measure.vars=c('inlanddivfisherfolk_M','inlanddivfisherfolk_F'))
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, SINK_TYPE:= gsub('([^_]+)_([^_]+)','\\1',variable)]
	z[, SINK_SEX:= gsub('([^_]+)_([^_]+)','\\2',variable)]	
	z		<- dcast.data.table(z, SINK_TYPE+SINK_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),' - ',round(CU,d=2),']')]
	z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'-',round(CU,d=2),')')]
	z[, STAT:='by gender']
	setkey(z, SINK_TYPE, SINK_SEX )
	ans		<- copy(z)
			
	#
	#	geography flows inland->fish / flows fish->inland
	#	by complex migration status	
	z		<- copy(mcpi)
	#	reset to complex migrant
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')	
	z		<- z[, list(PI=sum(PI)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_INMIGRATE','IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',TR_INMIGRATE,' ',REC_COMM_TYPE)]	
	z		<- dcast.data.table(z, IT~FLOW, value.var='PI')	
	set(z, NULL, c('inland resident inland','inland outmigrant inland', 'external inmigrant_from_external inland','external inmigrant_from_external fisherfolk','fisherfolk resident fisherfolk'), NULL)	
	z[ , inlanddivfisherfolk_resident:= z[['inland resident fisherfolk']] / z[['fisherfolk resident inland']] ]
	z[ , inlanddivfisherfolk_outmigrant:= z[['inland outmigrant fisherfolk']] / z[['fisherfolk outmigrant inland']] ]				
	z		<- melt(z, id.vars='IT', measure.vars=c('inlanddivfisherfolk_resident','inlanddivfisherfolk_outmigrant'))
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, SINK_TYPE:= gsub('([^_]+)_([^_]+)','\\1',variable)]
	z[, SINK_MIGRATIONSTATUS:= gsub('([^_]+)_([^_]+)','\\2',variable)]						
	z		<- dcast.data.table(z, SINK_TYPE+SINK_MIGRATIONSTATUS~P, value.var='Q')
	z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),' - ',round(CU,d=2),']')]
	z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'-',round(CU,d=2),')')]
	z[, SINK_SEX:= 'Any']	
	z[, STAT:='by migration status']
	setkey(z, SINK_TYPE, SINK_MIGRATIONSTATUS )	
	ans		<- rbind(z, ans, fill=TRUE)
			
	#
	#	add overall
	z		<- copy(mcpi)
	#	reset to simple
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	set(z, z[, which(TR_INMIGRATE=='resident')], 'TR_INMIGRATE', 'resident/outmigrant')	
	z		<- z[, list(PI=sum(PI)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_INMIGRATE','IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',TR_INMIGRATE,' ',REC_COMM_TYPE)]
	z		<- dcast.data.table(z, IT~FLOW, value.var='PI')
	set(z, NULL, c('inland resident/outmigrant inland','external inmigrant_from_external inland','external inmigrant_from_external fisherfolk','fisherfolk resident/outmigrant fisherfolk'), NULL)
	z[ , inlanddivfisherfolk:= z[['inland resident/outmigrant fisherfolk']] / z[['fisherfolk resident/outmigrant inland']] ]
	z		<- melt(z, id.vars='IT', measure.vars=c('inlanddivfisherfolk'))
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, SINK_TYPE:= gsub('([^_]+)_([^_]+)','\\1',variable)]
	z[, SINK_SEX:= 'Any']	
	z		<- dcast.data.table(z, SINK_TYPE+SINK_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),' - ',round(CU,d=2),']')]
	z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'-',round(CU,d=2),')')]
	z[, STAT:='overall']
	setkey(z, SINK_TYPE, SINK_SEX )	
	ans		<- rbind(z, ans, fill=TRUE)
	
	#	save
	save(ans, file=paste0(outfile.base,'_flows_netflows.rda'))
		
	#
	if(0)
	{
		ggplot(ans, aes(x=SINK_SEX, middle=M, min=CL, lower=IL, upper=IU, max=CU)) +
				geom_hline(yintercept=1, colour='grey50', size=1.5) +
				geom_boxplot(stat='identity', fill='grey50') +
				theme_bw() + theme(legend.position='bottom') +
				scale_y_log10(expand=c(0,0), breaks=c(1/4,1/3,1/2,2/3,1,3/2,2,3,4,8,16,32), labels=c('1/4','1/3','1/2','2/3','1','3/2','2','3','4','8','16','32')) +
				coord_flip(ylim=c(2/3,40)) +			
				labs(x='Gender\n', y='\nflows inland->fishing / flows fishing->inland') +
				guides(fill='none')
		ggsave(file=paste0(outfile.base,'_flows_sinkfishinlandgender.pdf'), w=6, h=3)	
	}	
}

RakaiFull.phylogeography.181006.flows.fishinlandmigrant<- function(infile.inference=NULL)
{
	require(data.table)	
	require(Hmisc)	
		
	if(is.null(infile.inference))
	{
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_core_inference_mcmc.rda"
	}			
	outfile.base		<- gsub('.rda$','',gsub('_core_inference','',infile.inference))
	load(infile.inference)
	
	#	thin MCMC output
	tmp			<- seq.int(2, nrow(mc$pars$Z), mc$sweep)
	mc$pars$S	<- mc$pars$S[tmp,,drop=FALSE]
	mc$pars$Z	<- mc$pars$Z[tmp,,drop=FALSE]
	mc$pars$PI	<- mc$pars$PI[tmp,,drop=FALSE]
	mc$pars$N	<- mc$pars$N[tmp,,drop=FALSE]
	gc()
	colnames(mc$pars$S)	<- paste0('S-',1:ncol(mc$pars$S))
	colnames(mc$pars$Z)	<- paste0('Z-',1:ncol(mc$pars$Z))
	colnames(mc$pars$PI)<- paste0('PI-',1:ncol(mc$pars$PI))
	colnames(mc$pars$N)	<- 'N'
	
	#
	#	prepare data.table of proportions
	dc[, REC_COMM_TYPE:= as.character(factor(substr(REC_COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fisherfolk','inland')))]
	dc[, TR_COMM_TYPE:= as.character(factor(substr(TR_COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fisherfolk','inland')))]
	mcpi	<- as.data.table(mc$pars$PI)
	mcpi[, IT:= seq.int(2, by=100, length.out=nrow(mcpi))]
	mcpi	<- melt(mcpi, id.vars='IT', variable.name='COUNT_ID', value.name='PI')
	set(mcpi, NULL, 'COUNT_ID', mcpi[, gsub('([A-Z]+)-([0-9]+)','\\2',COUNT_ID)])
	set(mcpi, NULL, 'COUNT_ID', mcpi[, as.integer(COUNT_ID)])
	tmp		<- subset(dc, select=c(REC_COMM_NUM_A, REC_SEX, REC_COMM_TYPE, REC_INMIGRATE, TR_COMM_NUM_A, TR_SEX, TR_INMIGRATE, TR_COMM_TYPE, TR_OBS, COUNT_ID))
	mcpi	<- merge(tmp, mcpi, by='COUNT_ID')
	
	
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')
	
	#
	#	by migration status (complex) and community status
	#
	
	#	flows
	z		<- copy(mcpi)
	#	reset to complex migration status
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')					
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	z		<- z[, list(PI=sum(PI)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_INMIGRATE','IT')]
	z		<- z[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by=c('TR_COMM_TYPE','TR_INMIGRATE','REC_COMM_TYPE')]
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE )
	z[, STAT:='joint']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- copy(z)	
	z		<- NULL	
	gc()
		
	#	WAIFM matrix
	groups	<- data.table(TR_COMM_TYPE=c('inland','inland','fisherfolk','fisherfolk','external'), TR_INMIGRATE=c('resident','outmigrant','resident','outmigrant','inmigrant_from_external'))
	z		<- lapply(1:nrow(groups), function(ii)
			{				
				z		<- copy(mcpi)				
				#	reset to complex migration status
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')					
				set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
				set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
				z		<- subset(z, TR_COMM_TYPE==groups$TR_COMM_TYPE[ii] & TR_INMIGRATE==groups$TR_INMIGRATE[ii])
				z		<- z[, list(PI=sum(PI)), by=c('REC_COMM_TYPE','IT')]
				z		<- z[, list(REC_COMM_TYPE=REC_COMM_TYPE, PI=PI/sum(PI)), by=c('IT')]				
				z		<- z[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by=c('REC_COMM_TYPE')]
				z[, TR_COMM_TYPE:= groups$TR_COMM_TYPE[ii] ]
				z[, TR_INMIGRATE:= groups$TR_INMIGRATE[ii] ]
				z
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE)
	z[, STAT:='waifm']	
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans, z)
	z		<- NULL
	gc()
	
	#	sources
	groups	<- data.table(REC_COMM_TYPE=c('inland','fisherfolk'))
	z		<- lapply(1:nrow(groups), function(ii)
			{	
				z		<- copy(mcpi)				
				#	reset to complex migration status
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')					
				set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
				set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')				
				z		<- subset(z, REC_COMM_TYPE==groups$REC_COMM_TYPE[ii])								
				z		<- z[, list(PI=sum(PI)), by=c('TR_COMM_TYPE','TR_INMIGRATE','IT')]
				z		<- z[, list(TR_COMM_TYPE=TR_COMM_TYPE, TR_INMIGRATE=TR_INMIGRATE, PI=PI/sum(PI)), by=c('IT')]
				z		<- z[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by=c('TR_COMM_TYPE','TR_INMIGRATE')]
				z[, REC_COMM_TYPE:= groups$REC_COMM_TYPE[ii] ]
				z	 
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_COMM_TYPE, TR_COMM_TYPE, TR_INMIGRATE)
	z[, STAT:='sources']	
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	z		<- NULL
	gc()
	
	#
	#	by migration status (simple) and community status
	#
	
	#	flows
	z		<- copy(mcpi)
	#	reset to simple migration status
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')		
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	set(z, z[, which(TR_INMIGRATE=='resident')], 'TR_INMIGRATE', 'resident/outmigrant')					
	z		<- z[, list(PI=sum(PI)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_INMIGRATE','IT')]
	z		<- z[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by=c('TR_COMM_TYPE','TR_INMIGRATE','REC_COMM_TYPE')]
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE )
	z[, STAT:='joint2']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	z		<- NULL	
	gc()
	
	#	WAIFM matrix
	groups	<- data.table(TR_COMM_TYPE=c('inland','fisherfolk','external'), TR_INMIGRATE=c('resident/outmigrant','resident/outmigrant','inmigrant_from_external'))
	z		<- lapply(1:nrow(groups), function(ii)
			{				
				z		<- copy(mcpi)				
				#	reset to simple migration status
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')	
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')		
				set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
				set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
				set(z, z[, which(TR_INMIGRATE=='resident')], 'TR_INMIGRATE', 'resident/outmigrant')					
				z		<- subset(z, TR_COMM_TYPE==groups$TR_COMM_TYPE[ii] & TR_INMIGRATE==groups$TR_INMIGRATE[ii])
				z		<- z[, list(PI=sum(PI)), by=c('REC_COMM_TYPE','IT')]
				z		<- z[, list(REC_COMM_TYPE=REC_COMM_TYPE, PI=PI/sum(PI)), by=c('IT')]				
				z		<- z[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by=c('REC_COMM_TYPE')]
				z[, TR_COMM_TYPE:= groups$TR_COMM_TYPE[ii] ]
				z[, TR_INMIGRATE:= groups$TR_INMIGRATE[ii] ]
				z
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE)
	z[, STAT:='waifm2']	
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans, z)
	z		<- NULL
	gc()
	
	#	sources
	groups	<- data.table(REC_COMM_TYPE=c('inland','fisherfolk'))
	z		<- lapply(1:nrow(groups), function(ii)
			{	
				z		<- copy(mcpi)				
				#	reset to simple migration status
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')	
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')		
				set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
				set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
				set(z, z[, which(TR_INMIGRATE=='resident')], 'TR_INMIGRATE', 'resident/outmigrant')	
				z		<- subset(z, REC_COMM_TYPE==groups$REC_COMM_TYPE[ii])								
				z		<- z[, list(PI=sum(PI)), by=c('TR_COMM_TYPE','TR_INMIGRATE','IT')]
				z		<- z[, list(TR_COMM_TYPE=TR_COMM_TYPE, TR_INMIGRATE=TR_INMIGRATE, PI=PI/sum(PI)), by=c('IT')]
				z		<- z[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by=c('TR_COMM_TYPE','TR_INMIGRATE')]
				z[, REC_COMM_TYPE:= groups$REC_COMM_TYPE[ii] ]
				z	 
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_COMM_TYPE, TR_COMM_TYPE, TR_INMIGRATE)
	z[, STAT:='sources2']	
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	z		<- NULL
	gc()
	
	#	save
	save(ans, file=paste0(outfile.base,'_flows_fishinlandmigration.rda'))
}

RakaiFull.phylogeography.181006.flows.fishinlandmigrantgender<- function(infile.inference=NULL)
{
	require(data.table)	
	require(Hmisc)	
	
	if(is.null(infile.inference))
	{
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_core_inference_mcmc.rda"
	}			
	outfile.base		<- gsub('.rda$','',gsub('_core_inference','',infile.inference))
	load(infile.inference)
	
	#	thin MCMC output
	tmp			<- seq.int(2, nrow(mc$pars$Z), mc$sweep)
	mc$pars$S	<- mc$pars$S[tmp,,drop=FALSE]
	mc$pars$Z	<- mc$pars$Z[tmp,,drop=FALSE]
	mc$pars$PI	<- mc$pars$PI[tmp,,drop=FALSE]
	mc$pars$N	<- mc$pars$N[tmp,,drop=FALSE]
	gc()
	colnames(mc$pars$S)	<- paste0('S-',1:ncol(mc$pars$S))
	colnames(mc$pars$Z)	<- paste0('Z-',1:ncol(mc$pars$Z))
	colnames(mc$pars$PI)<- paste0('PI-',1:ncol(mc$pars$PI))
	colnames(mc$pars$N)	<- 'N'
	
	#
	#	prepare data.table of proportions
	dc[, REC_COMM_TYPE:= as.character(factor(substr(REC_COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fisherfolk','inland')))]
	dc[, TR_COMM_TYPE:= as.character(factor(substr(TR_COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fisherfolk','inland')))]
	mcpi	<- as.data.table(mc$pars$PI)
	mcpi[, IT:= seq.int(2, by=100, length.out=nrow(mcpi))]
	mcpi	<- melt(mcpi, id.vars='IT', variable.name='COUNT_ID', value.name='PI')
	set(mcpi, NULL, 'COUNT_ID', mcpi[, gsub('([A-Z]+)-([0-9]+)','\\2',COUNT_ID)])
	set(mcpi, NULL, 'COUNT_ID', mcpi[, as.integer(COUNT_ID)])
	tmp		<- subset(dc, select=c(REC_COMM_NUM_A, REC_SEX, REC_COMM_TYPE, REC_INMIGRATE, TR_COMM_NUM_A, TR_SEX, TR_INMIGRATE, TR_COMM_TYPE, TR_OBS, COUNT_ID))
	mcpi	<- merge(tmp, mcpi, by='COUNT_ID')
	
	
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')
	
	#
	#	by migration status (complex), gender, community status
	#
	
	#	flows
	z		<- copy(mcpi)
	#	reset to complex migration status
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')					
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	z		<- z[, list(PI=sum(PI)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_INMIGRATE','TR_SEX','REC_SEX','IT')]
	z		<- z[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by=c('TR_COMM_TYPE','TR_INMIGRATE','TR_SEX','REC_COMM_TYPE','REC_SEX')]
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE+TR_SEX+REC_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE, TR_SEX )
	z[, STAT:='joint3']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- copy(z)	
	z		<- NULL	
	gc()
	
	#	WAIFM matrix
	groups	<- data.table(	TR_COMM_TYPE=c('inland','inland','fisherfolk','fisherfolk','external','inland','inland','fisherfolk','fisherfolk','external'), 
							TR_INMIGRATE=c('resident','outmigrant','resident','outmigrant','inmigrant_from_external','resident','outmigrant','resident','outmigrant','inmigrant_from_external'),
							TR_SEX=c('M','M','M','M','M', 'F','F','F','F','F'))
	z		<- lapply(1:nrow(groups), function(ii)
			{				
				z		<- copy(mcpi)				
				#	reset to complex migration status
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')					
				set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
				set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
				z		<- subset(z, TR_COMM_TYPE==groups$TR_COMM_TYPE[ii] & TR_INMIGRATE==groups$TR_INMIGRATE[ii] & TR_SEX==groups$TR_SEX[ii])
				z		<- z[, list(PI=sum(PI)), by=c('REC_COMM_TYPE','REC_SEX','IT')]
				z		<- z[, list(REC_COMM_TYPE=REC_COMM_TYPE, REC_SEX=REC_SEX, PI=PI/sum(PI)), by=c('IT')]				
				z		<- z[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by=c('REC_COMM_TYPE','REC_SEX')]
				z[, TR_COMM_TYPE:= groups$TR_COMM_TYPE[ii] ]
				z[, TR_INMIGRATE:= groups$TR_INMIGRATE[ii] ]
				z[, TR_SEX:= groups$TR_SEX[ii] ]
				z
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE+TR_SEX+REC_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE, TR_SEX)
	z[, STAT:='waifm3']	
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans, z)
	z		<- NULL
	gc()
	
	#	sources
	groups	<- data.table(	REC_COMM_TYPE=c('inland','fisherfolk','inland','fisherfolk'),
							REC_SEX=c('M','M','F','F'))
	z		<- lapply(1:nrow(groups), function(ii)
			{	
				z		<- copy(mcpi)				
				#	reset to complex migration status
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')					
				set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
				set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')				
				z		<- subset(z, REC_COMM_TYPE==groups$REC_COMM_TYPE[ii] & REC_SEX==groups$REC_SEX[ii])								
				z		<- z[, list(PI=sum(PI)), by=c('TR_COMM_TYPE','TR_SEX','TR_INMIGRATE','IT')]
				z		<- z[, list(TR_COMM_TYPE=TR_COMM_TYPE, TR_SEX=TR_SEX, TR_INMIGRATE=TR_INMIGRATE, PI=PI/sum(PI)), by=c('IT')]
				z		<- z[, list(P=qsn, Q=unname(quantile(PI, p=qs))), by=c('TR_COMM_TYPE','TR_INMIGRATE','TR_SEX')]
				z[, REC_COMM_TYPE:= groups$REC_COMM_TYPE[ii] ]
				z[, REC_SEX:= groups$REC_SEX[ii] ]
				z	 
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE+TR_SEX+REC_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_COMM_TYPE, TR_COMM_TYPE, TR_INMIGRATE, TR_SEX, REC_SEX)
	z[, STAT:='sources3']	
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	z		<- NULL
	gc()
	
	#	save
	save(ans, file=paste0(outfile.base,'_flows_fishinlandmigrationgender.rda'))
}

RakaiFull.phylogeography.180618.popviraemia<- function()
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
	require(rethinking)
	
	outfile.base	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/popviraemia"
	
	#
	#	prepare self reported ART data
	#	
	infile			<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender_age.rda'		
	load(infile)
	#	get those infected before round 17 and who participated in at least one round 15-16
	df		<- subset(de, HIV_1517==1)
	
	#
	#	process self report ART at each visit
	tmp		<- RakaiCirc.epi.get.info.170208()
	ra		<- tmp$ra		
	setnames(ra, 'ARVMED', 'SELFREPORTART')
	
	#	sum by self-reported ART status	
	#	infected at or before round 15
	tmp	<- subset(ra, VISIT==15 & FIRSTPOSDATE<=VISDATE)[, list(N=length(RID)), by=c('VISIT','SEX','COMM_NUM','COMM_NUM_A','SELFREPORTART')]
	#	infected at or before round 15.1
	tmp2<- subset(ra, VISIT==15.1 & FIRSTPOSDATE<=VISDATE)[, list(N=length(RID)), by=c('VISIT','SEX','COMM_NUM','COMM_NUM_A','SELFREPORTART')]
	tmp	<- rbind(tmp, tmp2)
	#	infected at or before round 16
	tmp2<- subset(ra, VISIT==16 & FIRSTPOSDATE<=VISDATE)[, list(N=length(RID)), by=c('VISIT','SEX','COMM_NUM','COMM_NUM_A','SELFREPORTART')]
	dff	<- rbind(tmp, tmp2)	
	set(dff, dff[, which(is.na(SELFREPORTART))], 'SELFREPORTART', 2L)	
	set(dff, NULL, 'SELFREPORTART', dff[, as.character(factor(SELFREPORTART, levels=c(0,1,2), labels=c('SLART_NO','SLART_YES','SLART_UNKNOWN')))])	
	dff		<- dcast.data.table(dff, COMM_NUM+COMM_NUM_A+VISIT+SEX~SELFREPORTART, value.var='N')	
	for(x in c('SLART_NO','SLART_UNKNOWN','SLART_YES'))
		set(dff, which(is.na(dff[[x]])), x, 0L)	
	dff		<- dcast.data.table(dff, COMM_NUM+COMM_NUM_A+VISIT~SEX, value.var=c('SLART_NO','SLART_UNKNOWN','SLART_YES'))
	dff[, HIV_F:= SLART_NO_F+SLART_UNKNOWN_F+SLART_YES_F]
	dff[, HIV_M:= SLART_NO_M+SLART_UNKNOWN_M+SLART_YES_M]
	dff[, SLART_KNOWN_F:= SLART_NO_F+SLART_YES_F]
	dff[, SLART_KNOWN_M:= SLART_NO_M+SLART_YES_M]
	
	
	#
	#	prepare prevalence data
	#
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/community_hivprev_byGenderStudyVisit.rda"
	load(infile)
	df	<- as.data.table(melt(female.negative, varnames=c('COMM_NUM', 'VISIT'), value.name='FEMALE_NEG'))	
	tmp	<- as.data.table(melt(male.negative, varnames=c('COMM_NUM', 'VISIT'), value.name='MALE_NEG'))	
	df	<- merge(df, tmp, by=c('COMM_NUM', 'VISIT'), all=TRUE)	
	tmp	<- as.data.table(melt(female.positive, varnames=c('COMM_NUM', 'VISIT'), value.name='FEMALE_POS'))
	df	<- merge(df, tmp, by=c('COMM_NUM', 'VISIT'), all=TRUE)
	tmp	<- as.data.table(melt(male.positive, varnames=c('COMM_NUM', 'VISIT'), value.name='MALE_POS'))
	df	<- merge(df, tmp, by=c('COMM_NUM', 'VISIT'), all=TRUE)
	#	add geo-locations
	load("~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/community_geography.rda")
	comgps	<- as.data.table(comgps)
	set(comgps, NULL, 'COMM_NUM', comgps[, as.integer(as.character(COMM_NUM))])
	df	<- merge(df, comgps, by='COMM_NUM', all.x=TRUE)	
	set(df, NULL, 'COMM_NUM', df[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',as.character(COMM_NUM)))))])
	for(x in c('FEMALE_NEG','MALE_NEG','FEMALE_POS','MALE_POS'))
		set(df, which(is.na(df[[x]])), x, 0L)
	df	<- subset(df, FEMALE_NEG!=0 | MALE_NEG!=0 | FEMALE_POS!=0 | MALE_POS!=0)
	#	sum by merged communities that are essentially the same
	df	<- df[, list(FEMALE_NEG=sum(FEMALE_NEG), MALE_NEG=sum(MALE_NEG), FEMALE_POS=sum(FEMALE_POS), MALE_POS=sum(MALE_POS), longitude=mean(longitude), latitude=mean(latitude)), by= c('COMM_NUM','VISIT')]
	#	add anonymized ID
	tmp	<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	tmp	<- unique(subset(tmp, select=c(COMM_NUM, COMM_NUM_A, COMM_TYPE)))
	df	<- merge(df, tmp, by='COMM_NUM', all.x=TRUE)
	#	subset(df, is.na(COMM_NUM_A))	#comms with missing anonymized IDs are some historic ones, not relevant
	df	<- subset(df, !is.na(COMM_NUM_A))	
	df[, FEMALE:= FEMALE_NEG+FEMALE_POS]
	df[, MALE:= MALE_NEG+MALE_POS]
	df	<- subset(df, VISIT%in%c(15,15.1,16))
	
	
	#
	#	merge
	#
	df	<- merge(df, dff, by=c('COMM_NUM','COMM_NUM_A','VISIT'))
	df[, POS:= FEMALE_POS+MALE_POS]
	df[, HIV:= HIV_F+HIV_M]
	
	save(df, file=paste0(outfile.base, '_data.rda'))
	#	check HIV makes sense
	#	--> some very small differences, ignore
	subset(df, POS != HIV)
	#	check unknown ART status
	#	--> not known for 2 individuals, ignore
	subset(df, SLART_UNKNOWN_F>0 | SLART_UNKNOWN_M>0)
	
	#	new clean community numbers	
	tmp	<- unique(subset(df, select=c(COMM_NUM)))
	tmp[, COMM_NUM2:= seq_len(nrow(tmp))]
	df	<- merge(df, tmp, by=c('COMM_NUM'))
	
	#	estimate prevalences and prevalence ratio for each community		
	mv.1 	<- map2stan(
			alist(
					# prevalence likelihood
					FEMALE_POS ~ dbinom(FEMALE, hiv_f),
					MALE_POS ~ dbinom(MALE, hiv_m),
					logit(hiv_f) <- logit_hiv_f[COMM_NUM2],
					logit(hiv_m) <- logit_hiv_m[COMM_NUM2],
					# self report ART likelihood
					# TODO make both informative on hiv_f and hiv_m
					SLART_YES_M ~ dbinom(SLART_KNOWN_M, combined_m),
					SLART_YES_F ~ dbinom(SLART_KNOWN_F, combined_f),
					combined_m 	<- art_m * selfreport_sensitivity + (1-art_m)*(1-selfreport_specificity),
					combined_f 	<- art_f * selfreport_sensitivity + (1-art_f)*(1-selfreport_specificity),
					logit(art_f) <- logit_art_f[COMM_NUM2],
					logit(art_m) <- logit_art_m[COMM_NUM2],
					# uninformative priors
					logit_hiv_f[COMM_NUM2] ~ dnorm(0, 10),
					logit_hiv_m[COMM_NUM2] ~ dnorm(0, 10),
					logit_art_f[COMM_NUM2] ~ dnorm(0, 10),
					logit_art_m[COMM_NUM2] ~ dnorm(0, 10),
					# informative priors
					selfreport_sensitivity ~ dunif(0.7, 0.83),
					selfreport_specificity ~ dunif(0.98, 1)	
			),
			data=as.data.frame(df), 
			start=list(	logit_hiv_f=rep(0,length(unique(df$COMM_NUM2))), 
						logit_hiv_m=rep(0,length(unique(df$COMM_NUM2))),
						logit_art_f=rep(0,length(unique(df$COMM_NUM2))), 
						logit_art_m=rep(0,length(unique(df$COMM_NUM2))),
						selfreport_sensitivity=0.77, selfreport_specificity=0.99
						),
			warmup=5e2, iter=5e3, chains=1, cores=4
			)	
	#print(mv.1)
	plot(precis(mv.1, depth=2))	
	#traceplot(mv.1)
	#pairs(mv.1)
	
	#	extract Monte Carlo samples
	post	<- extract.samples(mv.1)
	dff		<- data.table(COMM_NUM2=unique(df$COMM_NUM2))
	dff		<- dff[, list( 	MC_IT= seq_len(nrow(post$logit_hiv_f)),
							PREVHIV_FEMALE= logistic(post$logit_hiv_f[,COMM_NUM2]),
							PREVHIV_MALE= logistic(post$logit_hiv_m[,COMM_NUM2]),
							PROPART_FEMALE= logistic(post$logit_art_f[,COMM_NUM2]),
							PROPART_MALE= logistic(post$logit_art_m[,COMM_NUM2])
							), by='COMM_NUM2']
	tmp		<- unique(subset(df, select=c(COMM_NUM, COMM_NUM_A, COMM_NUM2)))
	dff		<- merge(tmp, dff, by='COMM_NUM2')
	#	define pop viraemia as in NEJM paper
	dff[, POPVIRAEMIA_MALE:= PREVHIV_MALE - 0.9*PROPART_MALE*PREVHIV_MALE]
	dff[, POPVIRAEMIA_FEMALE:= PREVHIV_FEMALE - 0.9*PROPART_FEMALE*PREVHIV_FEMALE]
	
	#	summarize
	dff		<- melt(dff, id.vars=c('COMM_NUM2','COMM_NUM','COMM_NUM_A','MC_IT'), measure.vars=c('PREVHIV_FEMALE','PREVHIV_MALE','POPVIRAEMIA_MALE','POPVIRAEMIA_FEMALE','PROPART_FEMALE','PROPART_MALE'))
	dff		<- dff[	, list(	STAT=c('M','CL','IL','IU','CU'),
							value= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))
							)
					, by=c('COMM_NUM2','COMM_NUM','COMM_NUM_A','variable')]
	set(dff, NULL, 'SEX', dff[, tolower(gsub('([A-Z]+)_([A-Z]+)','\\2',variable))])
	set(dff, NULL, 'variable', dff[, tolower(gsub('([A-Z]+)_([A-Z]+)','\\1',variable))])
	dff		<- dcast.data.table(dff, COMM_NUM2+COMM_NUM+COMM_NUM_A+SEX+variable~STAT, value.var='value')
	dff[, COMM_TYPE:='inland communities']	
	set(dff, dff[, which(substr(COMM_NUM_A,1,1)=='f')], 'COMM_TYPE', 'fishing sites')
	set(dff, NULL, 'SEX', dff[, gsub('male','men',gsub('female','women',SEX))])
	#	reorder
	tmp2	<- subset(dff, SEX=='men' & variable=='popviraemia', c(COMM_NUM, COMM_NUM_A, COMM_TYPE, M))	
	tmp2	<- tmp2[order(COMM_TYPE, M),]
	tmp2[, DUMMY:= seq_len(nrow(tmp2))]
	set(tmp2, NULL, 'COMM_NUM_A', tmp2[, factor(DUMMY, levels=DUMMY, labels=COMM_NUM_A)])
	tmp2	<- subset(tmp2, select=c(COMM_NUM,COMM_NUM_A))
	dff[, COMM_NUM_A:=NULL]
	dff		<- merge(dff, tmp2, by='COMM_NUM')
	set(dff, NULL, 'variable', dff[, factor(variable, levels=c('prevhiv','popviraemia','propart'), labels=c('Population infected','Population with viraemia','Proportion infected on ART'))])
	#	plot
	ggplot(dff, aes(x=COMM_NUM_A)) +
			geom_boxplot(aes(middle=M, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=variable), stat='identity') +
			theme_bw() + 
			facet_grid(SEX~COMM_TYPE, space='free', scale='free') +
			#scale_fill_manual(values=c('women'='hotpink2', 'men'='deepskyblue')) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0)) +
			labs(x='\ncommunity', y='HIV-1 prevalence and population viraemia\n', fill='')
	ggsave(file=paste0(outfile.base,'popviraemia_comm_MF_estimated_180618.pdf'), w=12, h=8)
	ggplot(subset(dff, variable!='Proportion infected on ART'), aes(x=COMM_NUM_A)) +
			geom_boxplot(aes(middle=M, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=variable), stat='identity') +
			theme_bw() + 
			facet_grid(SEX~COMM_TYPE, space='free', scale='free') +			
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,0.6)) +
			labs(x='\ncommunity', y='HIV-1 prevalence and population viraemia\n', fill='')
	ggsave(file=paste0(outfile.base,'popviraemiaB_comm_MF_estimated_180618.pdf'), w=12, h=8)
	
	
	subset(df, COMM_NUM_A=='agf')
	#
	#	STAN model that marginalises out the ART likelihood combinations
	#	of SLART-TP and SLART-FP
	#

	
	mv.2		<- list()
	mv.2$code	<- "
data{
    int<lower=1> N;
    int<lower=1> N_COMM_NUM2;
	int COMM_NUM2[N];
	// for ART likelihood
	int SLART_KNOWN_M[N];
	int SLART_KNOWN_F[N];
    int SLART_YES_M[N];    
    int SLART_YES_F[N];
	// for prevalence likelihood
    int MALE[N];
	int FEMALE[N];
	int MALE_POS[N];
	int FEMALE_POS[N];
}
parameters{
    vector[N_COMM_NUM2] logit_hiv_f;
    vector[N_COMM_NUM2] logit_hiv_m;
    vector[N_COMM_NUM2] logit_art_f;
    vector[N_COMM_NUM2] logit_art_m;
    real<lower=0.7,upper=0.83> selfreport_sensitivity;	
    real<lower=0.98,upper=1> selfreport_specificity;
}
model{
	// transformed parameters that we don t want to store in output
    vector[N] art_m;				// on probability scale
    vector[N] art_f;				// on probability scale
    vector[N] selfreport_art_f;		// on probability scale
    vector[N] selfreport_art_m;		// on probability scale
    vector[N] hiv_m;				// on logit scale
    vector[N] hiv_f;				// on logit scale
    // priors
	// priors for selfreport_specificity and selfreport_sensitivity defined implicitly through lower, upper
    logit_art_m ~ normal( 0 , 10 );
    logit_art_f ~ normal( 0 , 10 );
    logit_hiv_m ~ normal( 0 , 10 );
    logit_hiv_f ~ normal( 0 , 10 );
	// update transformed parameters
    for ( i in 1:N ) {
        art_m[i] = inv_logit( logit_art_m[COMM_NUM2[i]] );
		art_f[i] = inv_logit( logit_art_f[COMM_NUM2[i]] );
		hiv_m[i] = logit_hiv_m[COMM_NUM2[i]];
		hiv_f[i] = logit_hiv_f[COMM_NUM2[i]];
		selfreport_art_m[i] = art_m[i] * selfreport_sensitivity + (1 - art_m[i]) * (1 - selfreport_specificity);
		selfreport_art_f[i] = art_f[i] * selfreport_sensitivity + (1 - art_f[i]) * (1 - selfreport_specificity);
    }
	// likelihood terms
	target += binomial_logit_lpmf( MALE_POS | MALE , hiv_m );
	target += binomial_logit_lpmf( FEMALE_POS | FEMALE , hiv_f );
 	target += binomial_lpmf( SLART_YES_F | SLART_KNOWN_F , selfreport_art_f );
	target += binomial_lpmf( SLART_YES_M | SLART_KNOWN_M , selfreport_art_m );	
}
"
	mv.2$code		<- gsub('\t','',mv.2$code)
	mv.2$codefile	<- paste0(outfile.base,'_stanmv.2.stan')
	#write(m1$code, file=m1$codefile)	
	mv.2$data		<- as.list(df)
	mv.2$data[['N']]<- nrow(df)
	mv.2$data[['N_COMM_NUM2']] 	<- max(df$COMM_NUM2)
	mv.2$fit		<- stan(model_code=mv.2$code, 
			data=mv.2$data, 
			warmup=5e2, iter=5e3, chains=1, cores=4,
			init=list(list(	logit_hiv_f=rep(0,length(unique(df$COMM_NUM2))), 
							logit_hiv_m=rep(0,length(unique(df$COMM_NUM2))),
							logit_art_f=rep(0,length(unique(df$COMM_NUM2))), 
							logit_art_m=rep(0,length(unique(df$COMM_NUM2))),
							selfreport_sensitivity=0.77, selfreport_specificity=0.99)))				
	#
	#	compare inference -- this should be the same for mv.1 and mv.2
	#
	dff	<- as.data.table(precis(mv.1, depth=2, prob=0.95)@output)
	dff[, MODEL:='mv1']
	dff[, VAR:= rownames(precis(mv.1, depth=2, prob=0.95)@output)]
	tmp	<- as.data.table(precis(mv.2$fit, depth=2, prob=0.95)@output)
	tmp[, MODEL:='mv2']
	tmp[, VAR:= rownames(precis(mv.2$fit, depth=2, prob=0.95)@output)]
	dff	<- rbind(dff, tmp)
	setnames(dff, c('lower 0.95','upper 0.95'), c('CL','CU'))
	ggplot(dff, aes(y=Mean, x=VAR, ymin=CL, ymax=CU, colour=MODEL)) +
			geom_point(position=position_dodge(0.9)) +
			geom_errorbar(width=0.4, position=position_dodge(0.9)) +
			theme_bw() +
			coord_flip()
	ggsave(file=paste0(outfile.base,'popviraemia_stanmodelcomparison.pdf'), w=8, h=24)
	
	#dfo	<- copy(df)
	df	<- subset(df, COMM_NUM_A%in%c('asm','abm'))	
	df[, START_IDX_LATENT_M:= 1+c(0,cumsum(SLART_YES_M+1)[-length(SLART_YES_M)])]
	df[, START_IDX_LATENT_F:= 1+c(0,cumsum(SLART_YES_F+1)[-length(SLART_YES_F)])]
	mv.3		<- list()
	mv.3$code	<- '
data{
	int<lower=1> N;
	int<lower=1> N_COMM_NUM2;
	int COMM_NUM2[N];
	// for ART likelihood
	int SLART_KNOWN_M[N];
	int SLART_KNOWN_F[N];
	int SLART_YES_M[N];    
	int SLART_YES_F[N];
	// for prevalence likelihood
	int MALE[N];
	int FEMALE[N];
	int MALE_POS[N];
	int FEMALE_POS[N];
	// for marginalizing out latent SLART_YES_ONART, SLART_YES_NOTONART
 	int<lower=1> N_LATENT_M;
	int<lower=1> N_LATENT_F;
	int<lower=1> N_LATENT;
	int START_IDX_LATENT_M[N];
	int START_IDX_LATENT_F[N];
}
parameters{
	vector[N_COMM_NUM2] comm_hiv_f;
	vector[N_COMM_NUM2] comm_hiv_m;
	vector[N_COMM_NUM2] comm_art_f;
	vector[N_COMM_NUM2] comm_art_m;
	real<lower=0.7,upper=0.83> selfreport_sensitivity;	
	real<lower=0.98,upper=1> selfreport_specificity;
}
transformed parameters{
	// transformed parameters that we don t want to store in output
	vector[N] art_m;				// on probability scale
	vector[N] art_f;				// on probability scale
	vector[N] hiv_m;				// on probability scale
	vector[N] hiv_f;				// on probability scale
	// log probability of model with all latent variables
	vector[N_LATENT] lp;
	real tmp;
	// update transformed parameters
	for ( i in 1:N ) {
		art_m[i] = comm_art_m[COMM_NUM2[i]];
		art_f[i] = comm_art_f[COMM_NUM2[i]];
		hiv_m[i] = comm_hiv_m[COMM_NUM2[i]];
		hiv_f[i] = comm_hiv_f[COMM_NUM2[i]];		
	}
	// likelihood terms and priors for latent variables
	//print("start");
	for(i in 1:N)
	{
		tmp	= SLART_YES_M[i]+1;
		tmp = log(1/tmp);
		// checked: iterating through reported ART values correctly
		//print(SLART_YES_M[i]);
		for(k in 1:(SLART_YES_M[i]+1))	
		{
			// k-1 is latent SLART_YES_ONART_M
			//checked: all these values below make sense
			//print(k-1);
			//print(SLART_YES_M[i]-k+1);
			//print(SLART_KNOWN_M[i]);		
			//print(art_m[i]);	
			//print(selfreport_sensitivity);
			//print(tmp);
			//print(binomial_lpmf( k-1 | SLART_KNOWN_M[i] , art_m[i]*selfreport_sensitivity ));
			//print(binomial_lpmf( SLART_YES_M[i]-k+1 | SLART_KNOWN_M[i] , (1 - art_m[i]) * (1 - selfreport_specificity) ));
			//checked: lp[at idx] is correct sum
			lp[  START_IDX_LATENT_M[i]+k-1  ] = binomial_lpmf( k-1 | SLART_KNOWN_M[i] , art_m[i]*selfreport_sensitivity ) + binomial_lpmf( SLART_YES_M[i]-k+1 | SLART_KNOWN_M[i] , (1 - art_m[i]) * (1 - selfreport_specificity) );
		}
		//print(lp);
		tmp	= SLART_YES_F[i]+1;
		tmp = log(1/tmp);			
		//print(SLART_YES_F[i]);
		for(k in 1:(SLART_YES_F[i]+1))	
		{
			// k-1 is latent SLART_YES_ONART_F, and SLART_KNOWN_F[i]-(k-1) is latent SLART_YES_NOTONART_F
			lp[  N_LATENT_M+START_IDX_LATENT_F[i]+k-1  ] = binomial_lpmf( k-1 | SLART_KNOWN_F[i] , art_f[i]*selfreport_sensitivity ) + binomial_lpmf( SLART_YES_F[i]-k+1 | SLART_KNOWN_F[i] , (1 - art_f[i]) * (1 - selfreport_specificity) );
		}
		//checked: lp correctly filled in
		//print(lp);			
	}
	//print("done");
	// print(lp);
}
model{
	// priors
	// priors for selfreport_specificity and selfreport_sensitivity defined implicitly through lower, upper
	comm_art_m ~ beta( .5, .5);
	comm_art_f ~ beta( .5, .5);
	comm_hiv_m ~ beta( .5, .5);
	comm_hiv_f ~ beta( .5, .5);	
	// add up total log probability
	target += log_sum_exp(lp); 		 
	target += binomial_lpmf( MALE_POS | MALE , hiv_m );
	target += binomial_lpmf( FEMALE_POS | FEMALE , hiv_f );
}		
	
'
	mv.3$code		<- gsub('\t','',mv.3$code)
	mv.3$codefile	<- paste0(outfile.base,'_stanmv.3.stan')
	write(mv.3$code, file=mv.3$codefile)	
	mv.3$data		<- as.list(df)
	mv.3$data[['N']]<- nrow(df)
	mv.3$data[['N_COMM_NUM2']] 		<- max(df$COMM_NUM2)
	mv.3$data[['N_LATENT_M']] 		<- df[, cumsum(SLART_YES_M+1)[length(SLART_YES_M)]]
	mv.3$data[['N_LATENT_F']] 		<- df[, cumsum(SLART_YES_F+1)[length(SLART_YES_F)]]
	mv.3$data[['N_LATENT']]	<- mv.3$data[['N_LATENT_F']]+mv.3$data[['N_LATENT_M']]
	mv.3$fit		<- stan(file=mv.3$codefile, #model_code=mv.3$code, 
			data=mv.3$data, 
			#verbose=TRUE, diagnostic_file= paste0(mv.3$codefile,'.diagnostics.txt'),
			warmup=5e2, iter=5e3, chains=1,
			#warmup=5, iter=10, chains=1,
			init=list(list(	comm_hiv_f=rep(0.5,length(unique(df$COMM_NUM2))), 
							comm_hiv_m=rep(0.5,length(unique(df$COMM_NUM2))),
							comm_art_f=rep(0.5,length(unique(df$COMM_NUM2))), 
							comm_art_m=rep(0.5,length(unique(df$COMM_NUM2))),
							selfreport_sensitivity=0.77, selfreport_specificity=0.99)))				
	pdf(paste0(outfile.base,'popviraemia_stanmodel_try.pdf'), w=8, h=40)
	plot(precis(mv.3$fit, depth=2, prob=0.95, par=c('comm_hiv_f','comm_hiv_m','comm_art_f','comm_art_m','selfreport_sensitivity','selfreport_specificity')))
	dev.off()
	traceplot(mv.3$fit)
}

RakaiFull.phylogeography.180618.gender.mobility.core.inference<- function(infile.inference=NULL)
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
	require(rethinking)
	
	opt.adjust.sequencing.bias			<- 1
	opt.adjust.participation.bias		<- 1
	opt.exclude.onART.from.denominator	<- 0
	opt.set.missing.migloc.to.inland	<- 0
	opt.set.missing.migloc.to.fishing	<- 1-opt.set.missing.migloc.to.inland
	
	indir								<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run"
	
	if(is.null(infile.inference))
	{		
		infile.inference	<- file.path(indir,"todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda")		
	}
	
	outfile.base	<- gsub('180522','180618',gsub('data_with_inmigrants.rda','',infile.inference))
	load(infile.inference)
	
	#	select which inmigration data to work with
	setnames(rtr2, c('TR_INMIGRATE_2YR','REC_INMIGRATE_2YR'), c('TR_INMIGRATE','REC_INMIGRATE'))
	#	double check that TR_INMIGRATE does not have any unknown origin
	tmp	<- rtr2[, which(TR_INMIGRATE=='inmigrant_from_unknown')]
	if(opt.set.missing.migloc.to.inland)
		set(rtr2, tmp, 'TR_INMIGRATE', 'inmigrant_from_inland')
	if(opt.set.missing.migloc.to.fishing)
		set(rtr2, tmp, 'TR_INMIGRATE', 'inmigrant_from_fisherfolk')
	
	
	#	select variables that we need, define missing ones
	rtr2	<- subset(rtr2, select=c(PAIRID, TR_RID, REC_RID, TR_COMM_NUM_A, REC_COMM_NUM_A, TR_SEX, REC_SEX, TR_INMIGRATE, REC_INMIGRATE, TR_BIRTHDATE, REC_BIRTHDATE, TR_COMM_TYPE, REC_COMM_TYPE))
	set(rtr2, NULL, 'TR_COMM_NUM_A', rtr2[, as.character(TR_COMM_NUM_A)])
	set(rtr2, NULL, 'REC_COMM_NUM_A', rtr2[, as.character(REC_COMM_NUM_A)])		
	#	impute missing birth dates	
	set(rtr2, rtr2[, which(is.na(TR_BIRTHDATE))], 'TR_BIRTHDATE', rtr2[, mean(TR_BIRTHDATE, na.rm=TRUE)])
	set(rtr2, rtr2[, which(is.na(REC_BIRTHDATE))], 'REC_BIRTHDATE', rtr2[, mean(REC_BIRTHDATE, na.rm=TRUE)])
	#	set age at midpoint of study period
	rtr2[, TR_AGE_AT_MID:= 2013.25 - TR_BIRTHDATE]
	rtr2[, REC_AGE_AT_MID:= 2013.25 - REC_BIRTHDATE]
	#	stratify age
	rtr2[, TR_AGE_AT_MID_C:= as.character(cut(TR_AGE_AT_MID, breaks=c(10,25,35,65), labels=c('15-24','25-34','35+'), right=FALSE))]
	rtr2[, REC_AGE_AT_MID_C:= as.character(cut(REC_AGE_AT_MID, breaks=c(10,25,35,65), labels=c('15-24','25-34','35+'), right=FALSE))]
	stopifnot( nrow(subset(rtr2, is.na(TR_AGE_AT_MID_C)))==0 )
	stopifnot( nrow(subset(rtr2, is.na(REC_AGE_AT_MID_C)))==0 )
	#	delete what we don t need
	set(rtr2, NULL, c('TR_BIRTHDATE','REC_BIRTHDATE','TR_AGE_AT_MID','REC_AGE_AT_MID'), NULL)
	
	#	sum transmission events by community, gender, age, inmigration status
	dc	<- rtr2[, list(TR_OBS=length(PAIRID)), by=c('TR_COMM_NUM_A','REC_COMM_NUM_A','TR_SEX','REC_SEX','TR_INMIGRATE','REC_INMIGRATE','TR_AGE_AT_MID_C','REC_AGE_AT_MID_C')]
	setkey(dc, TR_COMM_NUM_A, REC_COMM_NUM_A, TR_SEX, TR_INMIGRATE, REC_INMIGRATE, TR_AGE_AT_MID_C, REC_AGE_AT_MID_C )
	# 	use the Berger objective prior with minimal loss compared to marginal Beta reference prior
	#	(https://projecteuclid.org/euclid.ba/1422556416)
	dc[, TR_PRIOR:= 0.8/nrow(dc)]
	
	#
	#	Bayesian model sampling prior: 
	#	define prior for sampling probabilities from previously fitted participation and sequencing models	
	#	sampling probs differ by: age community and sex
	#
	infile.participation	<- file.path(indir,"participation_differences_180322_logisticmodels.rda")	
	infile.sequencing		<- file.path(indir,"sequencing_differences_180322_logisticmodels.rda")
	if(opt.exclude.onART.from.denominator)
		infile.sequencing	<- file.path(indir,"sequencing_differences_180322_exclART_logisticmodels.rda")
		

	load(infile.participation)
	mp1 <- mp2 <- mp4 <- NULL
	#	binarize covariates
	dc[, TR_AGE_YOUNG:= as.integer(TR_AGE_AT_MID_C=='15-24')]
	dc[, TR_AGE_MID:= as.integer(TR_AGE_AT_MID_C=='25-34')]
	dc[, TR_MALE:= as.integer(TR_SEX=='M')]
	dc[, TR_COMM_TYPE_F:= as.integer(substr(TR_COMM_NUM_A,1,1)=='f')]
	dc[, TR_COMM_TYPE_T:= as.integer(substr(TR_COMM_NUM_A,1,1)=='t')]
	dc[, TR_INMIGRANT:= as.integer(grepl('inmigrant',TR_INMIGRATE))]	
	dc[, REC_AGE_YOUNG:= as.integer(REC_AGE_AT_MID_C=='15-24')]
	dc[, REC_AGE_MID:= as.integer(REC_AGE_AT_MID_C=='25-34')]
	dc[, REC_MALE:= as.integer(REC_SEX=='M')]
	dc[, REC_COMM_TYPE_F:= as.integer(substr(REC_COMM_NUM_A,1,1)=='f')]
	dc[, REC_COMM_TYPE_T:= as.integer(substr(REC_COMM_NUM_A,1,1)=='t')]
	dc[, REC_INMIGRANT:= as.integer(grepl('inmigrant',REC_INMIGRATE))]
	#	define community number to match STAN output for participation model
	tmp		<- unique(subset(dg, select=c(COMM_NUM_A,COMM_NUM_B)))
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc		<- merge(dc, tmp, by='TR_COMM_NUM_A')
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc		<- merge(dc, tmp, by='REC_COMM_NUM_A')
	#	define community number to match STAN output for sequencing model
	load(infile.sequencing)
	ms2 <- ms3 <- ms4 <- NULL
	tmp		<- unique(subset(dg, select=c(COMM_NUM_A,COMM_NUM_B)))
	setnames(tmp, 'COMM_NUM_B', 'COMM_NUM_B2')
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc		<- merge(dc, tmp, by='TR_COMM_NUM_A')
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc		<- merge(dc, tmp, by='REC_COMM_NUM_A')
	#	extract Monte Carlo samples from best WAIC participation and best WAIC sequencing models	
	mps			<- extract.samples(mp3)
	mss			<- extract.samples(ms1)
	mc.it		<- 5e4
	gc()
	#mc.it		<- 1e2
	dcb		<- dc[, {
				#	get Monte Carlo samples from posterior distribution of logit(participation probs)
				p.tr.part	<- with(mps, a + comm[TR_COMM_NUM_B] + trading*TR_COMM_TYPE_T + fishing*TR_COMM_TYPE_F + 
											inmigrant*TR_INMIGRANT + inmigrant_young*TR_INMIGRANT*TR_AGE_YOUNG + 
											male*TR_MALE + 
											young_male*TR_AGE_YOUNG*TR_MALE + young_female*TR_AGE_YOUNG*(1-TR_MALE) +
											midage*TR_AGE_MID)
				p.rec.part	<- with(mps, a + comm[REC_COMM_NUM_B] + trading*REC_COMM_TYPE_T + fishing*REC_COMM_TYPE_F + 
											inmigrant*REC_INMIGRANT + inmigrant_young*REC_INMIGRANT*REC_AGE_YOUNG + 
											male*REC_MALE + 
											young_male*REC_AGE_YOUNG*REC_MALE + young_female*REC_AGE_YOUNG*(1-REC_MALE) +
											midage*REC_AGE_MID)
				#	get Monte Carlo samples from posterior distribution of logit(sequencing probs)			
				p.tr.seq	<- with(mss, a + comm[TR_COMM_NUM_B2] + trading*TR_COMM_TYPE_T + fishing*TR_COMM_TYPE_F + 
											inmigrant*TR_INMIGRANT + male*TR_MALE + young*TR_AGE_YOUNG + midage*TR_AGE_MID)
				p.rec.seq	<- with(mss, a + comm[REC_COMM_NUM_B2] + trading*REC_COMM_TYPE_T + fishing*REC_COMM_TYPE_F + 
											inmigrant*REC_INMIGRANT + male*REC_MALE + young*REC_AGE_YOUNG + midage*REC_AGE_MID)
				#	sensitivity analyses			
				if(!opt.adjust.participation.bias)
				{
					p.tr.part <- p.rec.part <- rep(Inf, length(p.tr.part))
				}					
				if(!opt.adjust.sequencing.bias)
				{
					p.tr.seq <- p.rec.seq <- rep(Inf, length(p.tr.seq))
				}				
				#	sample mc.it many, and calculate product of sequencing and participation probs
				mc.idx		<- sample.int(length(p.tr.part), mc.it, replace=TRUE)				
				tmp			<- log1p(exp(-p.tr.part[mc.idx])) + log1p(exp(-p.rec.part[mc.idx]))
				mc.idx		<- sample.int(length(p.tr.seq), mc.it, replace=TRUE)
				tmp			<- tmp + log1p(exp(-p.tr.seq[mc.idx])) + log1p(exp(-p.rec.seq[mc.idx]))
				tmp			<- exp(-tmp)				
				#print(TR_COMM_NUM_B)
				#str(p.tr.seq)
				#str(p.rec.seq)
				#print(tmp)
				#print(summary(tmp))				
				tmp			<- rnbinom(mc.it, TR_OBS, tmp)
				#print(tmp)
				#print(mean(tmp))
				#stop()
				list(	MONTE_CARLO_IT=seq_len(mc.it), 
						TR_PRIOR=TR_PRIOR, 
						TR_OBS=TR_OBS, 
						TR_MISS=tmp)
			}, by=c('TR_COMM_NUM_A','REC_COMM_NUM_A','TR_SEX','REC_SEX','TR_INMIGRATE','REC_INMIGRATE','TR_AGE_AT_MID_C','REC_AGE_AT_MID_C')]				
	#
	#	Bayesian model augmented posterior: 
	#	the augmented likelihood is multinomial, so the augmented posterior is 
	#	analytically tractable, and a Dirichlet posterior
	#
	dcb[, PI_IJ_ALPHA:= TR_OBS+TR_MISS+TR_PRIOR]
	#
	#	this is the end of the source attribution inference on the WAIFM matrix
	#
	tmp		<- paste0(outfile.base,'core_inference.rda')
	if(!opt.adjust.sequencing.bias)
		tmp	<- gsub('core_inference','core_inference_noadjseqbias',tmp)		
	if(!opt.adjust.participation.bias)
		tmp	<- gsub('core_inference','core_inference_noadjpartbias',tmp)
	if(opt.exclude.onART.from.denominator)
		tmp	<- gsub('core_inference','core_inference_seqExclART',tmp)		
	
	save(zm, rtpdm, rtr2, ds, dc, dcb, mps, mss, file=tmp)	
}

RakaiFull.phylogeography.180618.figure.impact.of.nosamplingadjustments<- function()
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
	
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')
	mc.it	<- 1e2
	
	#
	#	geography who infects whom matrix  between fisherfolk and others
	#	process central results
	#
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_inference.rda"
	outfile.base		<- gsub('_inference.rda','',infile.inference)
	load(infile.inference)
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset transmitter location by migration status
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',REC_COMM_TYPE)]	
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\1',variable)]	
	z[, REC_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\2',variable)]		
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE )	
	z[, STAT:='adjusted for variation in participation and sequencing rates']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- copy(z)
	#	estimate sink
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))		
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset transmitter location by migration status
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',REC_COMM_TYPE)]
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z[ , inlanddivfisherfolk:= z[['inland fisherfolk']] / z[['fisherfolk inland']] ]
	z		<- melt(z, id.vars='MONTE_CARLO_IT', measure.vars=c('inlanddivfisherfolk'))
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, SINK_TYPE:= gsub('([^_]+)_([^_]+)','\\1',variable)]
	z[, SINK_SEX:= 'Any']						
	z		<- dcast.data.table(z, SINK_TYPE+SINK_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),' - ',round(CU,d=2),']')]
	z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'-',round(CU,d=2),')')]
	setkey(z, SINK_TYPE, SINK_SEX )
	z[, STAT:='adjusted for variation in participation and sequencing rates']
	z[, DUMMY:= seq_len(nrow(z))]
	ans2	<- copy(z)
	z		<- dcb	<- dc <- NULL
	gc()
	
	
	#
	#	geography who infects whom matrix  between fisherfolk and others
	#	process sensitivity results adjust only for variation in participation
	#
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_inference_noadjseqbias.rda"	
	load(infile.inference)
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset transmitter location by migration status
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',REC_COMM_TYPE)]	
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\1',variable)]	
	z[, REC_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\2',variable)]		
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE )
	z[, STAT:='adjusted for variation in participation rates']
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans,z)	
	#	estimate sink
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))		
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset transmitter location by migration status
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')	
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',REC_COMM_TYPE)]
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z[ , inlanddivfisherfolk:= z[['inland fisherfolk']] / z[['fisherfolk inland']] ]
	z		<- melt(z, id.vars='MONTE_CARLO_IT', measure.vars=c('inlanddivfisherfolk'))
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, SINK_TYPE:= gsub('([^_]+)_([^_]+)','\\1',variable)]
	z[, SINK_SEX:= 'Any']						
	z		<- dcast.data.table(z, SINK_TYPE+SINK_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),' - ',round(CU,d=2),']')]
	z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'-',round(CU,d=2),')')]
	setkey(z, SINK_TYPE, SINK_SEX )
	z[, STAT:='adjusted for variation in participation rates']
	z[, DUMMY:= nrow(ans2)+seq_len(nrow(z))]
	ans2	<- rbind(ans2,z)	
	z		<- dcb	<- dc <- NULL
	gc()
	
	
	#
	#	geography who infects whom matrix  between fisherfolk and others
	#	process sensitivity results adjust only for variation in sequencing
	#
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_inference_noadjpartbias.rda"	
	load(infile.inference)
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset transmitter location by migration status
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',REC_COMM_TYPE)]	
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\1',variable)]	
	z[, REC_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\2',variable)]		
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE )
	z[, STAT:='adjusted for variation in sequencing rates']
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans,z)	
	#	estimate sink
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))		
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset transmitter location by migration status
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')	
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',REC_COMM_TYPE)]
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z[ , inlanddivfisherfolk:= z[['inland fisherfolk']] / z[['fisherfolk inland']] ]
	z		<- melt(z, id.vars='MONTE_CARLO_IT', measure.vars=c('inlanddivfisherfolk'))
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, SINK_TYPE:= gsub('([^_]+)_([^_]+)','\\1',variable)]
	z[, SINK_SEX:= 'Any']						
	z		<- dcast.data.table(z, SINK_TYPE+SINK_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),' - ',round(CU,d=2),']')]
	z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'-',round(CU,d=2),')')]
	setkey(z, SINK_TYPE, SINK_SEX )
	z[, STAT:='adjusted for variation in sequencing rates']
	z[, DUMMY:= nrow(ans2)+seq_len(nrow(z))]
	ans2	<- rbind(ans2,z)	
	z		<- dcb	<- dc <- NULL
	gc()
	
	
	#
	#	geography who infects whom matrix  between fisherfolk and others
	#	process sensitivity results adjust only for variation in sequencing
	#
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_inference_noadjpartbias_noadjseqbias.rda"	
	load(infile.inference)
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset transmitter location by migration status
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',REC_COMM_TYPE)]	
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\1',variable)]	
	z[, REC_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\2',variable)]		
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE )
	z[, STAT:='no adjustments to variation in participation and sequencing rates']
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans,z)	
	#	estimate sink
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))		
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset transmitter location by migration status
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')	
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',REC_COMM_TYPE)]
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z[ , inlanddivfisherfolk:= z[['inland fisherfolk']] / z[['fisherfolk inland']] ]
	z		<- melt(z, id.vars='MONTE_CARLO_IT', measure.vars=c('inlanddivfisherfolk'))
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, SINK_TYPE:= gsub('([^_]+)_([^_]+)','\\1',variable)]
	z[, SINK_SEX:= 'Any']						
	z		<- dcast.data.table(z, SINK_TYPE+SINK_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),' - ',round(CU,d=2),']')]
	z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'-',round(CU,d=2),')')]
	setkey(z, SINK_TYPE, SINK_SEX )
	z[, STAT:='no adjustments to variation in participation and sequencing rates']
	z[, DUMMY:= nrow(ans2)+seq_len(nrow(z))]
	ans2	<- rbind(ans2,z)	
	z		<- dcb	<- dc <- NULL
	gc()
	
	
	#
	#	make plots
	#
	#	FLOWS
	tmp		<- copy(ans)
	setnames(tmp, c('TR_COMM_TYPE','REC_COMM_TYPE'), c('TR_FISHCOMM','REC_FISHCOMM'))
	set(tmp, NULL, 'TR_FISHCOMM', tmp[, gsub('external','External',gsub('inland','Inland',gsub('fisherfolk','Fishing',TR_FISHCOMM)))])
	set(tmp, NULL, 'REC_FISHCOMM', tmp[, gsub('external','External',gsub('inland','Inland',gsub('fisherfolk','Fishing',REC_FISHCOMM)))])		
	tmp[, X:=paste0(TR_FISHCOMM,'->',REC_FISHCOMM)]
	set(tmp, NULL, 'STAT', tmp[, factor(STAT, levels=c('no adjustments to variation in participation and sequencing rates','adjusted for variation in participation rates','adjusted for variation in sequencing rates','adjusted for variation in participation and sequencing rates'))])
	set(tmp, NULL, 'X', tmp[, factor(X, levels=c('Inland->Inland','Inland->Fishing','Fishing->Inland','Fishing->Fishing','External->Inland','External->Fishing'))])	
	ggplot(tmp, aes(x=X, fill=STAT)) +
			geom_bar(aes(y=M), position=position_dodge(0.9), stat='identity') +
			geom_errorbar(aes(ymin=CL, ymax=CU), position=position_dodge(0.9), width=0.4) +			
			theme_bw() + 		
			theme(axis.text.x = element_text(angle=45, hjust=1)) +
			#theme(legend.position='bottom') + guides( fill=guide_legend(ncol=2)) +
			scale_fill_brewer(palette='Spectral') +
			scale_y_continuous(label=scales:::percent, expand=c(0,0), lim=c(0,0.55)) +
			labs(x='', y='Estimated\nTransmission flows\n', fill='') 
	ggsave(file=paste0(outfile.base,'_compare_m5_sampling.pdf'), w=12, h=6)
	#	SINK
	tmp		<- copy(ans2)
	set(tmp, NULL, 'STAT', tmp[, factor(STAT, levels=c('no adjustments to variation in participation and sequencing rates','adjusted for variation in participation rates','adjusted for variation in sequencing rates','adjusted for variation in participation and sequencing rates'))])		
	ggplot(tmp, aes(x=SINK_TYPE, fill=STAT)) +
			geom_hline(yintercept=1, colour='grey50', size=1.5) +
			geom_boxplot(aes(middle=M, ymin=CL, ymax=CU, lower=IL, upper=IU), position=position_dodge(0.9), stat='identity') +						
			theme_bw() + 		
			scale_fill_brewer(palette='Spectral') +
			scale_y_log10(expand=c(0,0), breaks=c(1/10,1/5,1/4,1/3,1/2,2/3,1,2,3,4,5,10), labels=c('1/10','1/5','1/4','1/3','1/2','2/3','1','2','3','4','5','10')) +
			coord_flip(ylim=c(1/13,13)) +						
			labs(x='', y='\nInland->Fishing / Fishing->Inland', fill='') +
			guides(fill='none')
	ggsave(file=paste0(outfile.base,'_compare_m5_sampling_sink.pdf'), w=8, h=2)
}

RakaiFull.phylogeography.180618.flows.fishinlandmigrant<- function(infile.inference=NULL)
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
	
	if(is.null(infile.inference))
	{
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_inference.rda"
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_inference_seqExclART.rda"
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min10_phylogeography_core_inference.rda"
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min20_phylogeography_core_inference.rda"
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min50_phylogeography_core_inference.rda"		
	}
	
	outfile.base		<- gsub('.rda$','',gsub('_inference','',infile.inference))
	load(infile.inference)
	
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')
	
	#
	#	geography who infects whom matrix  between fisherfolk and others account for migration (complex)
	#	adjusted P
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset some TR_INMIGRATE
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_INMIGRATE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',TR_INMIGRATE,' ',REC_COMM_TYPE)]
	#unique(z$FLOW)
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+) ([^ ]+)','\\1',variable)]
	z[, TR_INMIGRATE:= gsub('([^ ]+) ([^ ]+) ([^ ]+)','\\2',variable)]
	z[, REC_COMM_TYPE:= gsub('([^ ]+) ([^ ]+) ([^ ]+)','\\3',variable)]		
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE )
	z[, STAT:='joint']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- copy(z)	
	z		<- NULL
	gc()
	
	#
	#	WAIFM matrix (complex model)
	#
	groups	<- data.table(TR_COMM_TYPE=c('inland','inland','fisherfolk','fisherfolk','external'), TR_INMIGRATE=c('resident','outmigrant','resident','outmigrant','inmigrant_from_external'))
	z		<- lapply(1:nrow(groups), function(ii)
			{				
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
				#set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				#	reset to complex migrant
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')					
				set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
				set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
				#
				z		<- subset(z, TR_COMM_TYPE==groups$TR_COMM_TYPE[ii] & TR_INMIGRATE==groups$TR_INMIGRATE[ii])				
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
				z[, TR_COMM_TYPE:= groups$TR_COMM_TYPE[ii] ]
				z[, TR_INMIGRATE:= groups$TR_INMIGRATE[ii] ] 
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE)
	z[, STAT:='waifm']	
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans, z)
	z		<- NULL
	gc()
	
	#
	#	sources (complex model)
	#
	groups	<- data.table(REC_COMM_TYPE=c('inland','fisherfolk'))
	z		<- lapply(1:nrow(groups), function(ii)
			{				
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
				#set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				#	reset to complex
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')					
				set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
				set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
				#
				z		<- subset(z, REC_COMM_TYPE==groups$REC_COMM_TYPE[ii])				
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('TR_COMM_TYPE','TR_INMIGRATE','MONTE_CARLO_IT')]
				z[, FLOW:= paste0(TR_COMM_TYPE, ' ',TR_INMIGRATE)]
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- FLOW
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	get quantiles
				z		<- melt(z, id.vars=c('MONTE_CARLO_IT'))								
				z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by='variable']
				z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\1',variable)]
				z[, TR_INMIGRATE:= gsub('([^ ]+) ([^ ]+)','\\2',variable)]				
				z[, REC_COMM_TYPE:= groups$REC_COMM_TYPE[ii] ]				 
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_COMM_TYPE, TR_COMM_TYPE, TR_INMIGRATE)
	z[, STAT:='sources']	
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	z		<- NULL
	gc()
	
	
	
	#
	#	geography who infects whom matrix  between fisherfolk and others account for migration (simple)
	#	adjusted P
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
	set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset 
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')		
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	set(z, z[, which(TR_INMIGRATE=='resident')], 'TR_INMIGRATE', 'resident/outmigrant')	
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_INMIGRATE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',TR_INMIGRATE,' ',REC_COMM_TYPE)]
	#unique(z$FLOW)
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+) ([^ ]+)','\\1',variable)]
	z[, TR_INMIGRATE:= gsub('([^ ]+) ([^ ]+) ([^ ]+)','\\2',variable)]
	z[, REC_COMM_TYPE:= gsub('([^ ]+) ([^ ]+) ([^ ]+)','\\3',variable)]		
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE )
	z[, STAT:='joint2']
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	z		<- NULL
	gc()
	
	
	
	#
	#	WAIFM matrix (simple model)
	#
	groups	<- data.table(TR_COMM_TYPE=c('inland','fisherfolk','external'), TR_INMIGRATE=c('resident/outmigrant','resident/outmigrant','inmigrant_from_external'))
	z		<- lapply(1:nrow(groups), function(ii)
			{				
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
				#set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				#	reset to simple
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')	
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')					
				set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
				set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
				set(z, z[, which(TR_INMIGRATE=='resident')], 'TR_INMIGRATE', 'resident/outmigrant')	
				#
				z		<- subset(z, TR_COMM_TYPE==groups$TR_COMM_TYPE[ii] & TR_INMIGRATE==groups$TR_INMIGRATE[ii])				
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
				z[, TR_COMM_TYPE:= groups$TR_COMM_TYPE[ii] ]
				z[, TR_INMIGRATE:= groups$TR_INMIGRATE[ii] ] 
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE)
	z[, STAT:='waifm2']	
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	z		<- NULL
	gc()
	
	#
	#	sources (simple model)
	#
	groups	<- data.table(REC_COMM_TYPE=c('inland','fisherfolk'))
	z		<- lapply(1:nrow(groups), function(ii)
			{				
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
				#set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				#	reset to simple
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
				set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')	
				tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
				set(z, tmp2, 'TR_COMM_TYPE', 'inland')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
				tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
				set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
				set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')					
				set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
				set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
				set(z, z[, which(TR_INMIGRATE=='resident')], 'TR_INMIGRATE', 'resident/outmigrant')	
				#
				z		<- subset(z, REC_COMM_TYPE==groups$REC_COMM_TYPE[ii])				
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('TR_COMM_TYPE','TR_INMIGRATE','MONTE_CARLO_IT')]
				z[, FLOW:= paste0(TR_COMM_TYPE, ' ',TR_INMIGRATE)]
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- FLOW
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	get quantiles
				z		<- melt(z, id.vars=c('MONTE_CARLO_IT'))								
				z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by='variable']
				z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\1',variable)]
				z[, TR_INMIGRATE:= gsub('([^ ]+) ([^ ]+)','\\2',variable)]				
				z[, REC_COMM_TYPE:= groups$REC_COMM_TYPE[ii] ]				 
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_COMM_TYPE, TR_COMM_TYPE, TR_INMIGRATE)
	z[, STAT:='sources2']	
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	z		<- NULL
	gc()
	
	
	
	#	save
	save(ans, file=paste0(outfile.base,'_flows_fishinlandmigration.rda'))
	
	
	#	make table
	df	<- subset(ans, !grepl('sources',STAT))
	df[, STAT2:= gsub('joint2','joint',STAT)]
	set(df, NULL, 'TR_COMM_TYPE', df[, factor(TR_COMM_TYPE, levels=c('inland','fisherfolk','external'))])
	set(df, NULL, 'REC_COMM_TYPE', df[, factor(REC_COMM_TYPE, levels=c('inland','fisherfolk'))])
	set(df, NULL, 'TR_INMIGRATE', df[, factor(TR_INMIGRATE, levels=c('resident/outmigrant','resident','outmigrant','inmigrant_from_external'))])
	df	<- subset(df, !(TR_COMM_TYPE=='external' & STAT=='joint2'))
	df	<- dcast.data.table(df, TR_COMM_TYPE+REC_COMM_TYPE+TR_INMIGRATE~STAT2, value.var='LABEL2')
	setkey(df, TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE)
	write.csv(df, row.names=FALSE, file=paste0(outfile.base, '_table_fishinlandmigration.csv'))
	
	if(0)
	{
		#	make plot
		tmp		<- copy(z)	
		set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','External',gsub('inland','Inland',gsub('fisherfolk','Fishing',TR_COMM_TYPE)))])
		set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','Inland',gsub('fisherfolk','Fishing',REC_COMM_TYPE))])
		set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('_',' ',TR_INMIGRATE)])
		set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('inmigrant from external','',TR_INMIGRATE)])
		tmp[, MODEL:= 'MonteCarloAdjustSampling']	
		tmp[, X:= paste0(TR_COMM_TYPE,'->',REC_COMM_TYPE, ' (', TR_INMIGRATE, ')')]
		ggplot(tmp, aes(x=X, fill=MODEL)) +
				geom_boxplot(aes(ymin=CL, lower=IL, upper=IU, ymax=CU, middle=M), stat='identity') +
				theme_bw() + 
				coord_flip() +
				scale_y_continuous(label=scales:::percent) +
				labs(x='Transmissions from -> to\n', y='\nProportion of transmissions')
		ggsave(file=paste0(outfile.base,'_flows_fishinlandmigration.pdf'), w=8, h=5)
		
		#	plot
		tmp		<- copy(z)	
		set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','External',gsub('inland','Inland',gsub('fisherfolk','Fishing',TR_COMM_TYPE)))])
		set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','Inland',gsub('fisherfolk','Fishing',REC_COMM_TYPE))])
		set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('_',' ',TR_INMIGRATE)])
		set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('inmigrant from external','',TR_INMIGRATE)])
		tmp[, MODEL:= 'MonteCarloAdjustSampling']	
		tmp[, X:= paste0(TR_COMM_TYPE,'->',REC_COMM_TYPE, ' (', TR_INMIGRATE, ')')]
		ggplot(tmp, aes(x=X, fill=MODEL)) +
				geom_boxplot(aes(ymin=CL, lower=IL, upper=IU, ymax=CU, middle=M), stat='identity') +
				theme_bw() + 
				coord_flip() +
				scale_y_continuous(label=scales:::percent) +
				labs(x='Transmissions from -> to\n', y='\nProportion of transmissions by source population')
		ggsave(file=paste0(outfile.base,'_waifm_fishinlandmigration.pdf'), w=8, h=5)
		
		#	plot
		tmp		<- copy(z)	
		set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','External',gsub('inland','Inland',gsub('fisherfolk','Fishing',TR_COMM_TYPE)))])
		set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','Inland',gsub('fisherfolk','Fishing',REC_COMM_TYPE))])
		set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('_',' ',TR_INMIGRATE)])
		set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('^resident$','resident/outmigrant',gsub('inmigrant from external','',TR_INMIGRATE))])
		tmp[, MODEL:= 'MonteCarloAdjustSampling']	
		tmp[, X:= paste0(TR_COMM_TYPE,'->',REC_COMM_TYPE, ' (', TR_INMIGRATE, ')')]
		ggplot(tmp, aes(x=X, fill=MODEL)) +
				geom_boxplot(aes(ymin=CL, lower=IL, upper=IU, ymax=CU, middle=M), stat='identity') +
				theme_bw() + 
				coord_flip() +
				scale_y_continuous(label=scales:::percent) +
				labs(x='Transmissions from -> to\n', y='\nProportion of transmissions')
		ggsave(file=paste0(outfile.base,'_flows_fishinlandmigration2.pdf'), w=8, h=5)
		
		#	plot
		tmp		<- copy(z)	
		set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','External',gsub('inland','Inland',gsub('fisherfolk','Fishing',TR_COMM_TYPE)))])
		set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','Inland',gsub('fisherfolk','Fishing',REC_COMM_TYPE))])
		set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('_',' ',TR_INMIGRATE)])
		set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('inmigrant from external','',TR_INMIGRATE)])
		tmp[, MODEL:= 'MonteCarloAdjustSampling']	
		tmp[, X:= paste0(TR_COMM_TYPE,'->',REC_COMM_TYPE, ' (', TR_INMIGRATE, ')')]
		ggplot(tmp, aes(x=X, fill=MODEL)) +
				geom_boxplot(aes(ymin=CL, lower=IL, upper=IU, ymax=CU, middle=M), stat='identity') +
				theme_bw() + 
				coord_flip() +
				scale_y_continuous(label=scales:::percent) +
				labs(x='Transmissions from -> to\n', y='\nProportion of transmissions by source population')
		ggsave(file=paste0(outfile.base,'_waifm_fishinlandmigration2.pdf'), w=8, h=5)
		
		
		
		#	make alluvial sources plot
		require(ggalluvial)
		df	<- subset(ans, STAT=='joint2', select=c(TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE, STAT, M))
		tmp	<- subset(ans, STAT=='sources2', select=c(TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE, STAT, M, LABEL2))
		ggplot(df, aes(weight= M, axis1= TR_COMM_TYPE, axis2 = REC_COMM_TYPE)) +
				geom_alluvium(aes(fill = M, width = 2/12), alpha=0.8) +
				geom_stratum(width = 2/12, fill='grey', color="black") +
				geom_label(stat = "stratum", label.strata = TRUE) +
				geom_label(stat = "flow", label=rep(tmp$LABEL2,2)) +
				scale_x_continuous(breaks = 1:2, labels = c("Source population", "Recipient population")) +
				scale_y_continuous(labels=scales:::percent, expand=c(0,0)) +
				scale_fill_gradient2(trans='log', breaks=c(0.05, 0.1, 0.2, 0.4), low="deepskyblue", mid="orange", high='red', midpoint=log(0.15)) +
				#scale_fill_manual(values=c('external'='grey50','fisherfolk'='firebrick4','inland'='darkgreen')) +
				theme_minimal()
		ggsave(file=paste0(outfile.base,'_alluvial_fishinlandmigration2.pdf'), w=5, h=8)
		#	https://stackoverflow.com/questions/43053375/weighted-sankey-alluvial-diagram-for-visualizing-discrete-and-continuous-panel
		#dff	<- melt(df, id.vars=c('TR_INMIGRATE'))
		#ggplot(data=d, aes(x = timeperiod, stratum = discretechoice, alluvium = individual, weight = continuouschoice)) +
		#			geom_stratum(aes(fill = discretechoice)) +
		#			geom_flow()
		
		#	make boxplot
		tmp	<- subset(ans, STAT=='joint2' | STAT=='waifm2', select=c(TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE, STAT, M, CL, CU, IL, IU))
		set(tmp, NULL, 'STAT', tmp[, gsub('joint2','Transmission flow\noverall',gsub('waifm2','Transmission flow\nby source',STAT))])
		set(tmp, NULL, 'STAT', tmp[, factor(STAT, levels=c('Transmission flow\noverall','Transmission flow\nby source'))])
		set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','From\nExternal',gsub('inland','From\nInland',gsub('fisherfolk','From\nFishing',TR_COMM_TYPE)))])
		set(tmp, NULL, 'TR_COMM_TYPE', tmp[, factor(TR_COMM_TYPE, levels=c('From\nInland','From\nFishing','From\nExternal'))])
		set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','To Inland',gsub('fisherfolk','To Fishing',REC_COMM_TYPE))])
		set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('_',' ',TR_INMIGRATE)])
		set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('inmigrant from external','',TR_INMIGRATE)])
		#tmp[, X:= paste0(TR_COMM_TYPE,'->',REC_COMM_TYPE, ' (', TR_INMIGRATE, ')')]
		tmp[, X:= REC_COMM_TYPE]
		ggplot(tmp, aes(x=X)) +
				geom_boxplot(aes(ymin=CL, lower=IL, upper=IU, ymax=CU, middle=M, fill=TR_COMM_TYPE), stat='identity') +
				theme_bw() + 			
				scale_y_continuous(label=scales:::percent) +
				scale_fill_manual(values=c('From\nExternal'='grey50', 'From\nInland'='DarkGreen', 'From\nFishing'='firebrick1')) +
				coord_flip() +
				labs(x='', y='\nHIV transmissions between RCCS communities and external locations') +
				facet_grid(TR_COMM_TYPE~STAT, scales='free', switch='y') +
				theme(strip.text.y = element_text(angle=180), strip.placement = "outside") +
				guides(fill='none')
		ggsave(file=paste0(outfile.base,'_alluvial_fishinlandmigration_jointwaifm.pdf'), w=8, h=5)
		
		#	make barplot
		tmp	<- subset(ans, STAT=='joint2' | STAT=='waifm2', select=c(TR_COMM_TYPE, REC_COMM_TYPE, TR_INMIGRATE, STAT, M, CL, CU, IL, IU))
		set(tmp, NULL, 'STAT', tmp[, gsub('joint2','Transmission flow\noverall',gsub('waifm2','Transmission flow\nby source',STAT))])
		set(tmp, NULL, 'STAT', tmp[, factor(STAT, levels=c('Transmission flow\noverall','Transmission flow\nby source'))])
		set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','From\nExternal',gsub('inland','From\nInland',gsub('fisherfolk','From\nFishing',TR_COMM_TYPE)))])
		set(tmp, NULL, 'TR_COMM_TYPE', tmp[, factor(TR_COMM_TYPE, levels=c('From\nInland','From\nFishing','From\nExternal'))])
		set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','To\nInland',gsub('fisherfolk','To\nFishing',REC_COMM_TYPE))])
		set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('_',' ',TR_INMIGRATE)])
		set(tmp, NULL, 'TR_INMIGRATE', tmp[, gsub('inmigrant from external','',TR_INMIGRATE)])
		#tmp[, X:= paste0(TR_COMM_TYPE,'->',REC_COMM_TYPE, ' (', TR_INMIGRATE, ')')]
		tmp[, X:= REC_COMM_TYPE]
		ggplot(tmp, aes(x=X)) +
				geom_bar(aes(y=M, fill=TR_COMM_TYPE), stat='identity') +
				geom_errorbar(aes(ymin=CL, ymax=CU), width=0.5) +
				theme_bw() + 			
				scale_y_continuous(label=scales:::percent) +
				scale_fill_manual(values=c('From\nExternal'='grey50', 'From\nInland'='DarkGreen', 'From\nFishing'='firebrick1')) +			
				labs(x='', y='HIV transmissions\nbetween RCCS communities and external locations\n') +
				facet_grid(STAT~TR_COMM_TYPE, scales='free') +			
				guides(fill='none')
		ggsave(file=paste0(outfile.base,'_alluvial_fishinlandmigration_jointwaifm_bar.pdf'), w=5, h=8)
		
	}
	
}

RakaiFull.phylogeography.180618.flows.fishinlandmigrationgender.netflows<- function(infile.inference=NULL)
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
	if(is.null(infile.inference))
	{
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_inference.rda"
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_inference_seqExclART.rda"
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min10_phylogeography_core_inference.rda"
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min20_phylogeography_core_inference.rda"
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min50_phylogeography_core_inference.rda"		
	}
	
	outfile.base		<- gsub('.rda$','',gsub('_inference','',infile.inference))
	load(infile.inference)
	
	
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')
	
	
	#
	#	geography flows inland->fish / flows fish->inland
	#	start by gender	
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))		
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset to simple
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	set(z, z[, which(TR_INMIGRATE=='resident')], 'TR_INMIGRATE', 'resident/outmigrant')	
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_SEX','REC_SEX','TR_INMIGRATE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',TR_INMIGRATE,' ',TR_SEX,' ',REC_COMM_TYPE,' ',REC_SEX)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	set(z, NULL, c('inland resident/outmigrant F inland M','external inmigrant_from_external F inland M','external inmigrant_from_external F fisherfolk M','fisherfolk resident/outmigrant F fisherfolk M'), NULL)
	set(z, NULL, c('inland resident/outmigrant M inland F','external inmigrant_from_external M inland F','external inmigrant_from_external M fisherfolk F','fisherfolk resident/outmigrant M fisherfolk F'), NULL)
	z[ , inlanddivfisherfolk_M:= z[['inland resident/outmigrant M fisherfolk F']] / z[['fisherfolk resident/outmigrant M inland F']] ]
	z[ , inlanddivfisherfolk_F:= z[['inland resident/outmigrant F fisherfolk M']] / z[['fisherfolk resident/outmigrant F inland M']] ]				
	z		<- melt(z, id.vars='MONTE_CARLO_IT', measure.vars=c('inlanddivfisherfolk_M','inlanddivfisherfolk_F'))
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, SINK_TYPE:= gsub('([^_]+)_([^_]+)','\\1',variable)]
	z[, SINK_SEX:= gsub('([^_]+)_([^_]+)','\\2',variable)]						
	z		<- dcast.data.table(z, SINK_TYPE+SINK_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),' - ',round(CU,d=2),']')]
	z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'-',round(CU,d=2),')')]
	setkey(z, SINK_TYPE, SINK_SEX )
	ans		<- copy(z)
	
	
	
	#
	#	geography flows inland->fish / flows fish->inland
	#	add by migration status	
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))		
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')	
	#	reset to complex migrant
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_INMIGRATE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',TR_INMIGRATE,' ',REC_COMM_TYPE)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	set(z, NULL, c('inland resident inland','inland outmigrant inland', 'external inmigrant_from_external inland','external inmigrant_from_external fisherfolk','fisherfolk resident fisherfolk'), NULL)	
	z[ , inlanddivfisherfolk_resident:= z[['inland resident fisherfolk']] / z[['fisherfolk resident inland']] ]
	z[ , inlanddivfisherfolk_outmigrant:= z[['inland outmigrant fisherfolk']] / z[['fisherfolk outmigrant inland']] ]				
	z		<- melt(z, id.vars='MONTE_CARLO_IT', measure.vars=c('inlanddivfisherfolk_resident','inlanddivfisherfolk_outmigrant'))
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, SINK_TYPE:= gsub('([^_]+)_([^_]+)','\\1',variable)]
	z[, SINK_MIGRATIONSTATUS:= gsub('([^_]+)_([^_]+)','\\2',variable)]						
	z		<- dcast.data.table(z, SINK_TYPE+SINK_MIGRATIONSTATUS~P, value.var='Q')
	z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),' - ',round(CU,d=2),']')]
	z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'-',round(CU,d=2),')')]
	z[, SINK_SEX:= 'Any']	
	setkey(z, SINK_TYPE, SINK_MIGRATIONSTATUS )	
	ans		<- rbind(z, ans, fill=TRUE)
	
	#
	#	add overall
	#
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))		
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	#	reset to simple
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_fisherfolk')]
	set(z, tmp2, 'TR_COMM_TYPE', 'fisherfolk')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='inland' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')], 'TR_INMIGRATE', 'resident')
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='inland' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')	
	tmp2	<- z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_inland')]
	set(z, tmp2, 'TR_COMM_TYPE', 'inland')
	set(z, tmp2, 'TR_INMIGRATE', 'resident/outmigrant')
	set(z, z[, which(TR_COMM_TYPE=='fisherfolk' & REC_COMM_TYPE=='fisherfolk' & TR_INMIGRATE=='inmigrant_from_fisherfolk')], 'TR_INMIGRATE', 'resident')
	set(z, z[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_COMM_TYPE', 'external')
	set(z, z[, which(TR_INMIGRATE=='resident')], 'TR_INMIGRATE', 'resident/outmigrant')	
	#	use stick-breaking property: marginal of Dirichlet is again Dirichlet, self normalising
	z		<- z[, list(	TR_OBS=sum(TR_OBS), 
					TR_MISS=sum(TR_MISS), 
					PI_ST_ALPHA=sum(PI_IJ_ALPHA)), 
			by=c('REC_COMM_TYPE','TR_COMM_TYPE','TR_INMIGRATE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',TR_INMIGRATE,' ',REC_COMM_TYPE)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	set(z, NULL, c('inland resident/outmigrant inland','external inmigrant_from_external inland','external inmigrant_from_external fisherfolk','fisherfolk resident/outmigrant fisherfolk'), NULL)
	z[ , inlanddivfisherfolk:= z[['inland resident/outmigrant fisherfolk']] / z[['fisherfolk resident/outmigrant inland']] ]
	z		<- melt(z, id.vars='MONTE_CARLO_IT', measure.vars=c('inlanddivfisherfolk'))
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, SINK_TYPE:= gsub('([^_]+)_([^_]+)','\\1',variable)]
	z[, SINK_SEX:= 'Any']						
	z		<- dcast.data.table(z, SINK_TYPE+SINK_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M, d=2), '\n[',round(CL,d=2),' - ',round(CU,d=2),']')]
	z[, LABEL2:= paste0(round(M, d=2), ' (',round(CL,d=2),'-',round(CU,d=2),')')]
	setkey(z, SINK_TYPE, SINK_SEX )
	ans		<- rbind(z, ans, fill=TRUE)
	
	#	save
	save(ans, file=paste0(outfile.base,'_flows_sinkfishinlandgender.rda'))
	
	
	#
	if(0)
	{
		ggplot(ans, aes(x=SINK_SEX, middle=M, min=CL, lower=IL, upper=IU, max=CU)) +
				geom_hline(yintercept=1, colour='grey50', size=1.5) +
				geom_boxplot(stat='identity', fill='grey50') +
				theme_bw() + theme(legend.position='bottom') +
				scale_y_log10(expand=c(0,0), breaks=c(1/4,1/3,1/2,2/3,1,3/2,2,3,4,8,16,32), labels=c('1/4','1/3','1/2','2/3','1','3/2','2','3','4','8','16','32')) +
				coord_flip(ylim=c(2/3,40)) +			
				labs(x='Gender\n', y='\nflows inland->fishing / flows fishing->inland') +
				guides(fill='none')
		ggsave(file=paste0(outfile.base,'_flows_sinkfishinlandgender.pdf'), w=6, h=3)	
	}
	
}

RakaiFull.phylogeography.180618.flows.wrapper<- function()
{
	#
	#	make data set
	#
	infiles		<- c("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_withmetadata.rda"
					, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min10_withmetadata.rda"
					, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min20_withmetadata.rda"
					, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min50_withmetadata.rda"
					, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf50_withmetadata.rda"
					, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf55_withmetadata.rda"
					, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_withmetadata.rda"
					, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf65_withmetadata.rda"
					, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf70_withmetadata.rda"
					)
	for(ii in seq_along(infiles))
	{
		infile.inference	<- infiles[[ii]]
		RakaiFull.phylogeography.180521.gender.mobility.data(infile.inference)		
	}
	
		
	#
	#	core inference
	#
	infiles		<- c("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min10_phylogeography_data_with_inmigrants.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min20_phylogeography_data_with_inmigrants.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min50_phylogeography_data_with_inmigrants.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf50_phylogeography_data_with_inmigrants.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf55_phylogeography_data_with_inmigrants.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf65_phylogeography_data_with_inmigrants.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf70_phylogeography_data_with_inmigrants.rda"
				)
	for(ii in seq_along(infiles))
	{
		infile.inference	<- infiles[[ii]]
		RakaiFull.phylogeography.180618.gender.mobility.core.inference(infile.inference)		
	}
			
	#
	#	extract target variables of interest
	#
	infiles		<- c(#"~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_inference.rda"
				#, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_inference_seqExclART.rda"
				#, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min10_phylogeography_core_inference.rda"
				#, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min20_phylogeography_core_inference.rda"
				#, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min50_phylogeography_core_inference.rda"
				"~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf50_phylogeography_core_inference.rda"
				, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf55_phylogeography_core_inference.rda"
				#, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_core_inference.rda"
				#, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf65_phylogeography_core_inference.rda"
				#, "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf70_phylogeography_core_inference.rda"
				)
	for(ii in seq_along(infiles))
	{
		infile.inference	<- infiles[[ii]]
		#RakaiFull.phylogeography.180618.flows.fishinlandgender(infile.inference)
		RakaiFull.phylogeography.180618.flows.fishinlandmigrant(infile.inference)
		RakaiFull.phylogeography.180618.flows.fishinlandmigrationgender.netflows(infile.inference)
	}
		
}

RakaiFull.phylogeography.180618.flows.fishinlandgender<- function(infile.inference)
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
	
	if(is.null(infile.inference))
	{
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_inference.rda"	
	}
	
	outfile.base		<- gsub('.rda$','',gsub('_inference','',infile.inference))
	load(infile.inference)
	
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')
	
	#
	#	geography who infects whom matrix  between fisherfolk and others account for migration (complex)
	#	adjusted P
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
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
	z[, FLOW:=paste0(TR_COMM_TYPE,' ',TR_SEX,' ',REC_COMM_TYPE, ' ', REC_SEX)]
	#unique(z$FLOW)
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+) ([^ ]+) ([^ ]+)','\\1',variable)]
	z[, TR_SEX:= gsub('([^ ]+) ([^ ]+) ([^ ]+) ([^ ]+)','\\2',variable)]
	z[, REC_COMM_TYPE:= gsub('([^ ]+) ([^ ]+) ([^ ]+) ([^ ]+)','\\3',variable)]
	z[, REC_SEX:= gsub('([^ ]+) ([^ ]+) ([^ ]+) ([^ ]+)','\\4',variable)]
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_SEX+REC_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_SEX, REC_SEX )
	z[, STAT:='joint_sex']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- copy(z)		
	z		<- NULL
	gc()
	#
	#	WAIFM matrix (simple model)
	#
	groups	<- data.table(TR_COMM_TYPE=c('inland','inland','fisherfolk','fisherfolk'), TR_SEX=c('M','F','M','F'))
	z		<- lapply(1:nrow(groups), function(ii)
			{				
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
				#set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				#
				z		<- subset(z, TR_COMM_TYPE==groups$TR_COMM_TYPE[ii] & TR_SEX==groups$TR_SEX[ii])				
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('REC_COMM_TYPE','REC_SEX','MONTE_CARLO_IT')]
				z[, FLOW:=paste0(REC_COMM_TYPE, ' ', REC_SEX)]
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- FLOW
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	get quantiles
				z		<- melt(z, id.vars=c('MONTE_CARLO_IT'))								
				z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by='variable']
				z[, REC_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\1',variable)]
				z[, REC_SEX:= gsub('([^ ]+) ([^ ]+)','\\2',variable)]				
				z[, TR_COMM_TYPE:= groups$TR_COMM_TYPE[ii] ]
				z[, TR_SEX:= groups$TR_SEX[ii] ] 
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_SEX+REC_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_SEX, REC_SEX)
	z[, STAT:='waifm_sex']	
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	z		<- NULL
	gc()
	
	
	#
	#	sources
	#	
	groups	<- data.table(REC_COMM_TYPE=c('inland','inland','fisherfolk','fisherfolk'), REC_SEX=c('M','F','M','F'))
	z		<- lapply(1:nrow(groups), function(ii)
			{		
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
				#set(tmp, tmp[, which(COMM_TYPE!='fisherfolk')], 'COMM_TYPE', 'inland')	
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				#
				z		<- subset(z, REC_COMM_TYPE==groups$REC_COMM_TYPE[ii] & REC_SEX==groups$REC_SEX[ii])				
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('TR_COMM_TYPE','TR_SEX','MONTE_CARLO_IT')]
				z[, FLOW:=paste0(TR_COMM_TYPE, ' ', TR_SEX)]
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- FLOW
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	get quantiles
				z		<- melt(z, id.vars=c('MONTE_CARLO_IT'))								
				z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by='variable']
				z[, TR_COMM_TYPE:= gsub('([^ ]+) ([^ ]+)','\\1',variable)]
				z[, TR_SEX:= gsub('([^ ]+) ([^ ]+)','\\2',variable)]				
				z[, REC_COMM_TYPE:= groups$REC_COMM_TYPE[ii] ]
				z[, REC_SEX:= groups$REC_SEX[ii] ] 							 
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE+TR_SEX+REC_SEX~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, TR_COMM_TYPE, REC_COMM_TYPE, TR_SEX, REC_SEX)
	z[, STAT:='sources_sex']	
	z[, DUMMY:= nrow(ans)+seq_len(nrow(z))]
	ans		<- rbind(ans, z)	
	z		<- NULL
	gc()
	
	#	save
	save(ans, file=paste0(outfile.base,'_flows_fishinlandgender.rda'))
	
	
	#	write table
	df	<- copy(ans)
	df	<- dcast.data.table(df, TR_COMM_TYPE+TR_SEX+REC_COMM_TYPE+REC_SEX~STAT, value.var='LABEL2')
	setkey(df, TR_COMM_TYPE, TR_SEX, REC_COMM_TYPE, REC_SEX)
	write.csv(df, row.names=FALSE, file=paste0(outfile.base, '_table_fishinlandgender.csv'))
	
	
	if(0)
	{
		#	make plot
		tmp		<- copy(z)	
		set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','External',gsub('inland','Inland',gsub('from_fisherfolk','Fishing',TR_COMM_TYPE)))])
		set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','Inland',gsub('fisherfolk','Fishing',REC_COMM_TYPE))])
		set(tmp, NULL, 'TR_SEX', tmp[, gsub('M','male',TR_SEX)])
		set(tmp, NULL, 'TR_SEX', tmp[, gsub('F','female',TR_SEX)])
		set(tmp, NULL, 'REC_SEX', tmp[, gsub('M','male',REC_SEX)])
		set(tmp, NULL, 'REC_SEX', tmp[, gsub('F','female',REC_SEX)])
		tmp[, MODEL:= 'MonteCarloAdjustSampling']	
		tmp[, X:= paste0(TR_SEX, ' ', TR_COMM_TYPE,'->',REC_SEX, ' ',REC_COMM_TYPE)]
		ggplot(tmp, aes(x=X, fill=MODEL)) +
				geom_boxplot(aes(ymin=CL, lower=IL, upper=IU, ymax=CU, middle=M), stat='identity') +
				theme_bw() + 
				coord_flip() +
				scale_y_continuous(label=scales:::percent) +
				labs(x='Transmissions from -> to\n', y='\nProportion of transmissions')
		ggsave(file=paste0(outfile.base,'_flows_fishinlandgender.pdf'), w=8, h=5)
		
		#	plot
		tmp		<- copy(z)	
		set(tmp, NULL, 'TR_COMM_TYPE', tmp[, gsub('external','External',gsub('inland','Inland',gsub('fisherfolk','Fishing',TR_COMM_TYPE)))])
		set(tmp, NULL, 'REC_COMM_TYPE', tmp[, gsub('inland','Inland',gsub('fisherfolk','Fishing',REC_COMM_TYPE))])
		set(tmp, NULL, 'TR_SEX', tmp[, gsub('M','male',TR_SEX)])
		set(tmp, NULL, 'TR_SEX', tmp[, gsub('F','female',TR_SEX)])
		set(tmp, NULL, 'REC_SEX', tmp[, gsub('M','male',REC_SEX)])
		set(tmp, NULL, 'REC_SEX', tmp[, gsub('F','female',REC_SEX)])	
		tmp[, MODEL:= 'MonteCarloAdjustSampling']	
		tmp[, X:= paste0(TR_SEX, ' ', TR_COMM_TYPE,'->',REC_SEX, ' ',REC_COMM_TYPE)]
		ggplot(tmp, aes(x=X, fill=MODEL)) +
				geom_boxplot(aes(ymin=CL, lower=IL, upper=IU, ymax=CU, middle=M), stat='identity') +
				theme_bw() + 
				coord_flip() +
				scale_y_continuous(label=scales:::percent) +
				labs(x='Transmissions from -> to\n', y='\nProportion of transmissions by source population')
		ggsave(file=paste0(outfile.base,'_waifm_fishinlandgender.pdf'), w=8, h=5)		
	}
	
}


RakaiFull.phylogeography.180618.figure.overall.source.attribution<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	
	infile			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_flows_fishinlandmigration.rda"	
	outfile.base	<- gsub('.rda$','',gsub('_core_flows_fishinlandmigration','',infile))
	
	#	fish inland by migration
	load(infile)
	ans[, SAMPLING:='regression_180618']
	ans[, ANALYSIS:='by migration status']
	df	<- copy(ans)
	
	#	clean up
	df	<- subset(df, grepl('2$',STAT))
	set(df, NULL, 'STAT', df[, gsub('2$','',STAT)])
	set(df, NULL, 'TR_COMM_TYPE', df[, gsub('from_','',TR_COMM_TYPE)])
	set(df, NULL, 'REC_COMM_TYPE', df[, gsub('to_','',REC_COMM_TYPE)])
	set(df, NULL, 'TR_COMM_TYPE', df[, factor(TR_COMM_TYPE, levels=c('inland','fisherfolk','external'), labels=c('inland communities','fishing sites','external'))])
	set(df, NULL, 'REC_COMM_TYPE', df[, factor(REC_COMM_TYPE, levels=c('inland','fisherfolk','external'), labels=c('inland communities','fishing sites','external'))])
	for(x in c('TR_INMIGRATE'))
		set(df, which(is.na(df[[x]])), x, '')
	set(df, df[, which(TR_INMIGRATE=='resident/outmigrant')], 'TR_INMIGRATE', '')
	set(df, df[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_INMIGRATE', '')
	set(df, NULL, 'STAT', df[, factor(STAT, levels=c("joint","sources", "waifm"), labels=c("joint","sources of", "destinations of"))])
		
	#	plot
	ggplot(subset(df, STAT=='sources of'), aes(x=TR_COMM_TYPE, y=M, ymin=CL, lower=IL, upper=IU, ymax=CU, fill=TR_COMM_TYPE)) +			
			geom_bar(stat='identity') +
			geom_point() +
			geom_errorbar(width=0.5) +
			theme_bw() +
			theme(legend.position='bottom', panel.spacing = unit(2, "lines")) +
			scale_y_continuous(label=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_fill_manual(values=c('inland communities'="#66BD63",'fishing sites'="firebrick1",'external'='deepskyblue1')) +
			facet_grid(.~STAT+REC_COMM_TYPE, scales='free', space='free') +
			coord_flip() +
			labs(x='', y='', fill='source population')
	ggsave(file=paste0(outfile.base,'_central_sources.pdf'), w=7, h=2.3, useDingbats=FALSE)
	ggplot(subset(df, STAT=='destinations of'), aes(x=REC_COMM_TYPE, y=M, ymin=CL, lower=IL, upper=IU, ymax=CU, fill=TR_COMM_TYPE)) +			
			geom_bar(stat='identity') +
			geom_point() +
			geom_errorbar(width=0.5) +
			theme_bw() +
			theme(legend.position='bottom', panel.spacing = unit(2, "lines")) +
			scale_y_continuous(label=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_fill_manual(values=c('inland communities'="#66BD63",'fishing sites'="firebrick1",'external'='deepskyblue1')) +
			facet_grid(.~STAT+TR_COMM_TYPE, scales='free', space='free') +
			coord_flip() +
			labs(x='', y='', fill='source population')
	ggsave(file=paste0(outfile.base,'_central_destinations.pdf'), w=10, h=2.1, useDingbats=FALSE)
	
	curvature	<- -0.2	
	arrow		<- arrow(length=unit(0.05, "npc"), type="closed", angle=15)
	dff	<- subset(df, STAT=='joint')
	set(dff, NULL, c('X','Y','X2','Y2','S'), NA_real_)
	set(dff, dff[, which(TR_COMM_TYPE=='external')], 'X', 0)
	set(dff, dff[, which(TR_COMM_TYPE=='external')], 'Y', 5)
	set(dff, dff[, which(TR_COMM_TYPE=='external')], 'S', 7e3)
	set(dff, dff[, which(TR_COMM_TYPE=='inland communities')], 'X', 1)
	set(dff, dff[, which(TR_COMM_TYPE=='inland communities')], 'Y', 1)
	set(dff, dff[, which(TR_COMM_TYPE=='inland communities')], 'S', 29100)
	set(dff, dff[, which(REC_COMM_TYPE=='inland communities')], 'X2', 1)
	set(dff, dff[, which(REC_COMM_TYPE=='inland communities')], 'Y2', 1)
	set(dff, dff[, which(TR_COMM_TYPE=='fishing sites')], 'X', 3)
	set(dff, dff[, which(TR_COMM_TYPE=='fishing sites')], 'Y', 1)
	set(dff, dff[, which(REC_COMM_TYPE=='fishing sites')], 'X2', 3)
	set(dff, dff[, which(REC_COMM_TYPE=='fishing sites')], 'Y2', 1)	
	set(dff, dff[, which(TR_COMM_TYPE=='fishing sites')], 'S', 8500)
	ggplot(dff) +
			geom_point(	data=unique(subset(dff, select=c(TR_COMM_TYPE, X, Y, S))), 
						aes(x=X, y=Y, size=S, colour=TR_COMM_TYPE)) +
			geom_curve(	data=subset(dff, X!=X2), 
						aes(x=X, xend=X2, y=Y, yend=Y2, colour=TR_COMM_TYPE), 
						curvature=curvature, arrow=arrow, lineend="butt", linejoin='mitre', size=1) +
			geom_curve(	data=subset(dff, X==X2), 
						aes(x=X+0.005, xend=X2+0.03, y=Y-0.005, yend=Y, colour=TR_COMM_TYPE), 
						curvature=1) +			
			geom_text(	data=dff, aes(x=(X+X2)/2, y=(Y+Y2)/2, label=LABEL)) +
			theme_void() +
			scale_size(range=c(5,30)) +
			scale_colour_manual(values=c('inland communities'="#66BD63",'fishing sites'="firebrick1",'external'='deepskyblue1')) +
			scale_x_continuous(lim=c(-.5,3.5)) +
			scale_y_continuous(lim=c(-1,6)) +
			guides(colour='none', size='none')
	ggsave(file=paste0(outfile.base,'_central_flows.pdf'), w=7, h=4, useDingbats=FALSE)
	
	#
	#	load sink estimates 
	#
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_flows_sinkfishinlandgender.rda")
	ans[, SAMPLING:='regression_180618']
	ans[, ANALYSIS:='fishinland']
	df	<- copy(ans)
	
	#	clean up
	df[, X:=NA_character_]
	tmp	<- df[, which(grepl('Any',SINK_SEX) & is.na(SINK_MIGRATIONSTATUS))]
	set(df, tmp, 'ANALYSIS', 'overall')
	set(df, tmp, 'X', 'overall')
	tmp	<- df[, which(grepl('M|F',SINK_SEX) & is.na(SINK_MIGRATIONSTATUS))]
	set(df, tmp, 'ANALYSIS', 'by gender')
	set(df, tmp, 'X', df[tmp, gsub('F','women',gsub('M','men',SINK_SEX))])
	tmp	<- df[, which(grepl('Any',SINK_SEX) & !is.na(SINK_MIGRATIONSTATUS))]
	set(df, tmp, 'ANALYSIS', 'by migration status')
	set(df, tmp, 'X', df[tmp, gsub('resident','resident',gsub('outmigrant','outmigrant',SINK_MIGRATIONSTATUS))])
	set(df, NULL, 'ANALYSIS', df[, factor(ANALYSIS, levels=c("overall","by gender", "by migration status"), labels=c("overall","by\ngender", "by\nmigration status"))])
	
	df[,CU2:= CU]
	set(df, df[, which(CU2>65)], 'CU2', 65)
	
	#	plot
	tmp	<- c(2,4,10,20,50)
	ggplot(df, aes(x=X, middle=M, ymin=CL, lower=IL, upper=IU, ymax=CU2, fill=X)) +
			geom_hline(yintercept=1, colour='grey50') +
			geom_boxplot(stat='identity') +
			theme_bw() +
			theme(legend.position='bottom', strip.text.y = element_text(angle=0)) +
			scale_y_log10(expand=c(0,0), breaks=c(1/tmp, 1, tmp), labels=c(paste0('1/',as.character(tmp)), '1', as.character(tmp)), lim=range(c(1/(tmp*1.31), tmp*1.31))) +
			scale_fill_manual(values=c('overall'='grey60', 'men'='deepskyblue1', 'women'='hotpink1', 'resident'='darkolivegreen3', 'outmigrant'='tan3')) +
			facet_grid(ANALYSIS~., scales='free', space='free') +
			coord_flip() +
			labs(x='', y='\ntransmissions inland->fishing sites  /  transmissions fishing sites->inland', fill='sampling\nmodel') +
			guides(fill='none')
	ggsave(file=paste0(outfile.base,'_central_sink.pdf'), w=7, h=2.2, useDingbats=FALSE)
	
	
}
	

RakaiFull.phylogeography.180618.figure.impact.of.sampling.models<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	outfile.base	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run"
	
	#	overall fish-inland 
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_fishing_inland_results.rda")
	ans[, SAMPLING:='regression_180618']
	ans[, ANALYSIS:='fishinland']
	df	<- copy(ans)
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_seqExclART_fishing_inland_results.rda")
	ans[, SAMPLING:='regression_exclART_180618']
	ans[, ANALYSIS:='fishinland']
	df	<- rbind(df, copy(ans))	
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_core_fishing_inland_results.rda")
	ans[, SAMPLING:='empirical-ratios_180522']
	ans[, ANALYSIS:='fishinland']
	df	<- rbind(df, copy(ans))
	
	
	#	fish inland by gender
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_flows_fishinlandgender.rda")
	ans[, SAMPLING:='regression_180618']
	ans[, ANALYSIS:='by gender']
	df	<- rbind(df, copy(ans), fill=TRUE)
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_seqExclART_flows_fishinlandgender.rda")
	ans[, SAMPLING:='regression_exclART_180618']
	ans[, ANALYSIS:='by gender']
	df	<- rbind(df, copy(ans), fill=TRUE)	
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_core_flows_fishinlandgender.rda")
	ans[, SAMPLING:='empirical-ratios_180522']
	ans[, ANALYSIS:='by gender']
	df	<- rbind(df, copy(ans), fill=TRUE)
	
	#	fish inland by migration
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_flows_fishinlandmigration.rda")
	ans[, SAMPLING:='regression_180618']
	ans[, ANALYSIS:='by migration status']
	df	<- rbind(df, copy(ans), fill=TRUE)
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_seqExclART_flows_fishinlandmigration.rda")
	ans[, SAMPLING:='regression_exclART_180618']
	ans[, ANALYSIS:='by migration status']
	df	<- rbind(df, copy(ans), fill=TRUE)	
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_core_flows_fishinlandmigration.rda")
	ans[, SAMPLING:='empirical-ratios_180522']
	ans[, ANALYSIS:='by migration status']
	df	<- rbind(df, copy(ans), fill=TRUE)
	
	#	clean up
	tmp	<- df[, which(grepl('2$',STAT))]
	set(df, tmp, 'ANALYSIS', df[tmp, paste0(ANALYSIS,'2')])
	set(df, tmp, 'STAT', df[tmp, gsub('2$','',STAT)])
	set(df, NULL, 'STAT', df[, gsub('_sex','',STAT)])
	set(df, NULL, 'TR_COMM_TYPE', df[, gsub('from_','',TR_COMM_TYPE)])
	set(df, NULL, 'REC_COMM_TYPE', df[, gsub('to_','',REC_COMM_TYPE)])
	for(x in c('TR_SEX','REC_SEX','TR_INMIGRATE'))
		set(df, which(is.na(df[[x]])), x, '')
	set(df, df[, which(TR_INMIGRATE=='resident/outmigrant')], 'TR_INMIGRATE', '')
	set(df, df[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_INMIGRATE', '')
	set(df, NULL, 'X', df[, gsub(' +',' ',paste0(TR_COMM_TYPE,' ',TR_SEX,' ',TR_INMIGRATE,' -> ',REC_COMM_TYPE,' ',REC_SEX))])
	df	<- subset(df, ANALYSIS!='fishinland')
	set(df, NULL, 'ANALYSIS', df[, gsub('by migration status2','overall',ANALYSIS)])
	set(df, NULL, 'ANALYSIS', df[, factor(ANALYSIS, levels=c("overall","by gender", "by migration status"))])
	set(df, NULL, 'STAT', df[, factor(STAT, levels=c("joint","sources", "waifm"), labels=c("joint","sources", "destinations"))])
	set(df, NULL, 'SAMPLING', df[, factor(SAMPLING, levels=c("regression_180618","regression_exclART_180618", "empirical-ratios_180522"))])
	#	plot
	ggplot(df, aes(x=X, middle=M, ymin=CL, lower=IL, upper=IU, ymax=CU, fill=SAMPLING)) +			
			geom_boxplot(stat='identity') +
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(label=scales:::percent, expand=c(0,0)) +
			facet_grid(ANALYSIS~STAT, scales='free', space='free') +
			coord_flip() +
			labs(x='', y='', fill='sampling\nmodel')
	ggsave(file=file.path(outfile.base,'samplingmodels_comparison.pdf'), w=15, h=15)
	
	#
	#	load sink estimates 
	#
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_flows_sinkfishinlandgender.rda")
	ans[, SAMPLING:='regression_180618']
	ans[, ANALYSIS:='fishinland']
	df	<- copy(ans)
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_seqExclART_flows_sinkfishinlandgender.rda")
	ans[, SAMPLING:='regression_exclART_180618']
	ans[, ANALYSIS:='fishinland']
	df	<- rbind(df, copy(ans))	
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180522_cl25_d50_prior23_min30_phylogeography_core_flows_sinkfishinlandgender.rda")
	ans[, SAMPLING:='empirical-ratios_180522']
	ans[, ANALYSIS:='fishinland']
	df	<- rbind(df, copy(ans))
	
	#	clean up
	df[, X:=NA_character_]
	set(df, NULL, 'SAMPLING', df[, factor(SAMPLING, levels=c("regression_180618","regression_exclART_180618", "empirical-ratios_180522"))])
	tmp	<- df[, which(grepl('Any',SINK_SEX) & is.na(SINK_MIGRATIONSTATUS))]
	set(df, tmp, 'ANALYSIS', 'overall')
	set(df, tmp, 'X', 'overall')
	tmp	<- df[, which(grepl('M|F',SINK_SEX) & is.na(SINK_MIGRATIONSTATUS))]
	set(df, tmp, 'ANALYSIS', 'by gender')
	set(df, tmp, 'X', df[tmp, gsub('F','women',gsub('M','men',SINK_SEX))])
	tmp	<- df[, which(grepl('Any',SINK_SEX) & !is.na(SINK_MIGRATIONSTATUS))]
	set(df, tmp, 'ANALYSIS', 'by migration status')
	set(df, tmp, 'X', df[tmp, gsub('resident','resident',gsub('outmigrant','outmigrant',SINK_MIGRATIONSTATUS))])
	set(df, NULL, 'ANALYSIS', df[, factor(ANALYSIS, levels=c("overall","by gender", "by migration status"))])
	
	df[,CU2:= CU]
	set(df, df[, which(CU2>65)], 'CU2', 65)
	
	#	plot
	tmp	<- c(2,4,10,20,50)
	ggplot(df, aes(x=X, middle=M, ymin=CL, lower=IL, upper=IU, ymax=CU2, fill=SAMPLING)) +
			geom_hline(yintercept=1, colour='grey50') +
			geom_boxplot(stat='identity') +
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_log10(expand=c(0,0), breaks=c(1/tmp, 1, tmp), labels=c(paste0('1/',as.character(tmp)), '1', as.character(tmp)), lim=range(c(1/(tmp*1.31), tmp*1.31))) +
			facet_grid(ANALYSIS~., scales='free', space='free') +
			coord_flip() +
			labs(x='', y='\ntransmissions inland->fishing sites /\n   transmissions fishing sites->inland', fill='sampling\nmodel')
	ggsave(file=file.path(outfile.base,'samplingmodels_comparison_sink.pdf'), w=15, h=6)
	
}

RakaiFull.phylogeography.180618.figure.impact.of.transmissioneventscutoff<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	outfile.base	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run"
	
	
	#	fish inland by gender
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_flows_fishinlandgender.rda")
	ans[, SAMPLING:='cut60']
	ans[, ANALYSIS:='by gender']
	df	<- copy(ans)	
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf50_phylogeography_core_flows_fishinlandgender.rda")
	ans[, SAMPLING:='cut50']
	ans[, ANALYSIS:='by gender']
	df	<- rbind(df, copy(ans), fill=TRUE)	
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf55_phylogeography_core_flows_fishinlandgender.rda")
	ans[, SAMPLING:='cut55']
	ans[, ANALYSIS:='by gender']
	df	<- rbind(df, copy(ans), fill=TRUE)
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf65_phylogeography_core_flows_fishinlandgender.rda")
	ans[, SAMPLING:='cut65']
	ans[, ANALYSIS:='by gender']
	df	<- rbind(df, copy(ans), fill=TRUE)
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf70_phylogeography_core_flows_fishinlandgender.rda")
	ans[, SAMPLING:='cut70']
	ans[, ANALYSIS:='by gender']
	df	<- rbind(df, copy(ans), fill=TRUE)
	
	#	fish inland by migration
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_flows_fishinlandmigration.rda")
	ans[, SAMPLING:='cut60']
	ans[, ANALYSIS:='by migration status']
	df	<- rbind(df, copy(ans), fill=TRUE)
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf50_phylogeography_core_flows_fishinlandmigration.rda")
	ans[, SAMPLING:='cut50']
	ans[, ANALYSIS:='by migration status']
	df	<- rbind(df, copy(ans), fill=TRUE)	
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf55_phylogeography_core_flows_fishinlandmigration.rda")
	ans[, SAMPLING:='cut55']
	ans[, ANALYSIS:='by migration status']
	df	<- rbind(df, copy(ans), fill=TRUE)
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf65_phylogeography_core_flows_fishinlandmigration.rda")
	ans[, SAMPLING:='cut65']
	ans[, ANALYSIS:='by migration status']
	df	<- rbind(df, copy(ans), fill=TRUE)
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf70_phylogeography_core_flows_fishinlandmigration.rda")
	ans[, SAMPLING:='cut70']
	ans[, ANALYSIS:='by migration status']
	df	<- rbind(df, copy(ans), fill=TRUE)
	
	#	clean up
	tmp	<- df[, which(grepl('2$',STAT))]
	set(df, tmp, 'ANALYSIS', df[tmp, paste0(ANALYSIS,'2')])
	set(df, tmp, 'STAT', df[tmp, gsub('2$','',STAT)])
	set(df, NULL, 'STAT', df[, gsub('_sex','',STAT)])
	set(df, NULL, 'TR_COMM_TYPE', df[, gsub('from_','',TR_COMM_TYPE)])
	set(df, NULL, 'REC_COMM_TYPE', df[, gsub('to_','',REC_COMM_TYPE)])
	for(x in c('TR_SEX','REC_SEX','TR_INMIGRATE'))
		set(df, which(is.na(df[[x]])), x, '')
	set(df, df[, which(TR_INMIGRATE=='resident/outmigrant')], 'TR_INMIGRATE', '')
	set(df, df[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_INMIGRATE', '')
	set(df, NULL, 'X', df[, gsub(' +',' ',paste0(TR_COMM_TYPE,' ',TR_SEX,' ',TR_INMIGRATE,' -> ',REC_COMM_TYPE,' ',REC_SEX))])
	df	<- subset(df, ANALYSIS!='fishinland')
	set(df, NULL, 'ANALYSIS', df[, gsub('by migration status2','overall',ANALYSIS)])
	set(df, NULL, 'ANALYSIS', df[, factor(ANALYSIS, levels=c("overall","by gender", "by migration status"))])
	set(df, NULL, 'STAT', df[, factor(STAT, levels=c("joint","sources", "waifm"), labels=c("joint","sources", "destinations"))])
	set(df, NULL, 'SAMPLING', df[, factor(SAMPLING, levels=c("cut50","cut55", "cut60", "cut65", "cut70"), labels=c("50%","55%", "60% (central)", "65%", "70%"))])
	#	plot
	ggplot(df, aes(x=X, middle=M, ymin=CL, lower=IL, upper=IU, ymax=CU, fill=SAMPLING)) +			
			geom_boxplot(stat='identity') +
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(label=scales:::percent, expand=c(0,0)) +
			facet_grid(ANALYSIS~STAT, scales='free', space='free') +
			coord_flip() +
			labs(x='', y='', fill='minimum proportion of phylogenies\nin support of linkage and directionality')
	ggsave(file=file.path(outfile.base,'trmeventcutoff_comparison.pdf'), w=12, h=14)
	
	#
	#	load sink estimates 
	#
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_flows_sinkfishinlandgender.rda")
	ans[, SAMPLING:='cut60']
	ans[, ANALYSIS:='fishinland']
	df	<- copy(ans)
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf50_phylogeography_core_flows_sinkfishinlandgender.rda")
	ans[, SAMPLING:='cut50']
	ans[, ANALYSIS:='fishinland']
	df	<- rbind(df, copy(ans))	
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf55_phylogeography_core_flows_sinkfishinlandgender.rda")
	ans[, SAMPLING:='cut55']
	ans[, ANALYSIS:='fishinland']
	df	<- rbind(df, copy(ans))
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf65_phylogeography_core_flows_sinkfishinlandgender.rda")
	ans[, SAMPLING:='cut65']
	ans[, ANALYSIS:='fishinland']
	df	<- rbind(df, copy(ans))
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_conf70_phylogeography_core_flows_sinkfishinlandgender.rda")
	ans[, SAMPLING:='cut70']
	ans[, ANALYSIS:='fishinland']
	df	<- rbind(df, copy(ans))
		
	#	clean up
	df[, X:=NA_character_]
	set(df, NULL, 'SAMPLING', df[, factor(SAMPLING, levels=c("cut50","cut55", "cut60", "cut65", "cut70"), labels=c("50%","55%", "60% (central)", "65%", "70%"))])
	tmp	<- df[, which(grepl('Any',SINK_SEX) & is.na(SINK_MIGRATIONSTATUS))]
	set(df, tmp, 'ANALYSIS', 'overall')
	set(df, tmp, 'X', 'overall')
	tmp	<- df[, which(grepl('M|F',SINK_SEX) & is.na(SINK_MIGRATIONSTATUS))]
	set(df, tmp, 'ANALYSIS', 'by gender')
	set(df, tmp, 'X', df[tmp, gsub('F','women',gsub('M','men',SINK_SEX))])
	tmp	<- df[, which(grepl('Any',SINK_SEX) & !is.na(SINK_MIGRATIONSTATUS))]
	set(df, tmp, 'ANALYSIS', 'by migration status')
	set(df, tmp, 'X', df[tmp, gsub('resident','resident',gsub('outmigrant','outmigrant',SINK_MIGRATIONSTATUS))])
	set(df, NULL, 'ANALYSIS', df[, factor(ANALYSIS, levels=c("overall","by gender", "by migration status"))])
	
	df[,CU2:= CU]
	set(df, df[, which(CU2>65)], 'CU2', 65)
	
	#	plot
	tmp	<- c(2,4,10,20,50)
	ggplot(df, aes(x=X, middle=M, ymin=CL, lower=IL, upper=IU, ymax=CU2, fill=SAMPLING)) +
			geom_hline(yintercept=1, colour='grey50') +
			geom_boxplot(stat='identity') +
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_log10(expand=c(0,0), breaks=c(1/tmp, 1, tmp), labels=c(paste0('1/',as.character(tmp)), '1', as.character(tmp)), lim=range(c(1/(tmp*1.31), tmp*1.31))) +
			facet_grid(ANALYSIS~., scales='free', space='free') +
			coord_flip() +
			labs(x='', y='\ntransmissions inland->fishing sites /\n   transmissions fishing sites->inland', fill='minimum proportion of phylogenies\nin support of linkage and directionality')
	ggsave(file=file.path(outfile.base,'trmeventcutoff_comparison_sink.pdf'), w=12, h=5)
	
}

RakaiFull.phylogeography.180618.figure.impact.of.deepseqcutoff<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	outfile.base	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run"
	
	
	#	fish inland by gender
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_flows_fishinlandgender.rda")
	ans[, SAMPLING:='cut30']
	ans[, ANALYSIS:='by gender']
	df	<- copy(ans)	
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min10_phylogeography_core_flows_fishinlandgender.rda")
	ans[, SAMPLING:='cut10']
	ans[, ANALYSIS:='by gender']
	df	<- rbind(df, copy(ans), fill=TRUE)	
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min20_phylogeography_core_flows_fishinlandgender.rda")
	ans[, SAMPLING:='cut20']
	ans[, ANALYSIS:='by gender']
	df	<- rbind(df, copy(ans), fill=TRUE)
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min50_phylogeography_core_flows_fishinlandgender.rda")
	ans[, SAMPLING:='cut50']
	ans[, ANALYSIS:='by gender']
	df	<- rbind(df, copy(ans), fill=TRUE)
	
	#	fish inland by migration
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_flows_fishinlandmigration.rda")
	ans[, SAMPLING:='cut30']
	ans[, ANALYSIS:='by migration status']
	df	<- rbind(df, copy(ans), fill=TRUE)
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min10_phylogeography_core_flows_fishinlandmigration.rda")
	ans[, SAMPLING:='cut10']
	ans[, ANALYSIS:='by migration status']
	df	<- rbind(df, copy(ans), fill=TRUE)	
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min20_phylogeography_core_flows_fishinlandmigration.rda")
	ans[, SAMPLING:='cut20']
	ans[, ANALYSIS:='by migration status']
	df	<- rbind(df, copy(ans), fill=TRUE)
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min50_phylogeography_core_flows_fishinlandmigration.rda")
	ans[, SAMPLING:='cut50']
	ans[, ANALYSIS:='by migration status']
	df	<- rbind(df, copy(ans), fill=TRUE)
	
	#	clean up
	tmp	<- df[, which(grepl('2$',STAT))]
	set(df, tmp, 'ANALYSIS', df[tmp, paste0(ANALYSIS,'2')])
	set(df, tmp, 'STAT', df[tmp, gsub('2$','',STAT)])
	set(df, NULL, 'STAT', df[, gsub('_sex','',STAT)])
	set(df, NULL, 'TR_COMM_TYPE', df[, gsub('from_','',TR_COMM_TYPE)])
	set(df, NULL, 'REC_COMM_TYPE', df[, gsub('to_','',REC_COMM_TYPE)])
	for(x in c('TR_SEX','REC_SEX','TR_INMIGRATE'))
		set(df, which(is.na(df[[x]])), x, '')
	set(df, df[, which(TR_INMIGRATE=='resident/outmigrant')], 'TR_INMIGRATE', '')
	set(df, df[, which(TR_INMIGRATE=='inmigrant_from_external')], 'TR_INMIGRATE', '')
	set(df, NULL, 'X', df[, gsub(' +',' ',paste0(TR_COMM_TYPE,' ',TR_SEX,' ',TR_INMIGRATE,' -> ',REC_COMM_TYPE,' ',REC_SEX))])
	df	<- subset(df, ANALYSIS!='fishinland')
	set(df, NULL, 'ANALYSIS', df[, gsub('by migration status2','overall',ANALYSIS)])
	set(df, NULL, 'ANALYSIS', df[, factor(ANALYSIS, levels=c("overall","by gender", "by migration status"))])
	set(df, NULL, 'STAT', df[, factor(STAT, levels=c("joint","sources", "waifm"), labels=c("joint","sources", "destinations"))])
	set(df, NULL, 'SAMPLING', df[, factor(SAMPLING, levels=c("cut10","cut20", "cut30", "cut50"), labels=c("10X","20X", "30X (central)", "50X"))])
	#	plot
	ggplot(df, aes(x=X, middle=M, ymin=CL, lower=IL, upper=IU, ymax=CU, fill=SAMPLING)) +			
			geom_boxplot(stat='identity') +
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(label=scales:::percent, expand=c(0,0)) +
			facet_grid(ANALYSIS~STAT, scales='free', space='free') +
			coord_flip() +
			labs(x='', y='', fill='minimum\nsequencing depth')
	ggsave(file=file.path(outfile.base,'deepseqcutoff_comparison.pdf'), w=12, h=12)
	
	#
	#	load sink estimates 
	#
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_flows_sinkfishinlandgender.rda")
	ans[, SAMPLING:='cut30']
	ans[, ANALYSIS:='fishinland']
	df	<- copy(ans)
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min10_phylogeography_core_flows_sinkfishinlandgender.rda")
	ans[, SAMPLING:='cut10']
	ans[, ANALYSIS:='fishinland']
	df	<- rbind(df, copy(ans))	
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min20_phylogeography_core_flows_sinkfishinlandgender.rda")
	ans[, SAMPLING:='cut20']
	ans[, ANALYSIS:='fishinland']
	df	<- rbind(df, copy(ans))
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min50_phylogeography_core_flows_sinkfishinlandgender.rda")
	ans[, SAMPLING:='cut50']
	ans[, ANALYSIS:='fishinland']
	df	<- rbind(df, copy(ans))
	
	
	#	clean up
	df[, X:=NA_character_]
	set(df, NULL, 'SAMPLING', df[, factor(SAMPLING, levels=c("cut10","cut20", "cut30", "cut50"), labels=c("10X","20X", "30X (central)", "50X"))])
	tmp	<- df[, which(grepl('Any',SINK_SEX) & is.na(SINK_MIGRATIONSTATUS))]
	set(df, tmp, 'ANALYSIS', 'overall')
	set(df, tmp, 'X', 'overall')
	tmp	<- df[, which(grepl('M|F',SINK_SEX) & is.na(SINK_MIGRATIONSTATUS))]
	set(df, tmp, 'ANALYSIS', 'by gender')
	set(df, tmp, 'X', df[tmp, gsub('F','women',gsub('M','men',SINK_SEX))])
	tmp	<- df[, which(grepl('Any',SINK_SEX) & !is.na(SINK_MIGRATIONSTATUS))]
	set(df, tmp, 'ANALYSIS', 'by migration status')
	set(df, tmp, 'X', df[tmp, gsub('resident','resident',gsub('outmigrant','outmigrant',SINK_MIGRATIONSTATUS))])
	set(df, NULL, 'ANALYSIS', df[, factor(ANALYSIS, levels=c("overall","by gender", "by migration status"))])
	
	df[,CU2:= CU]
	set(df, df[, which(CU2>65)], 'CU2', 65)
	
	#	plot
	tmp	<- c(2,4,10,20,50)
	ggplot(df, aes(x=X, middle=M, ymin=CL, lower=IL, upper=IU, ymax=CU2, fill=SAMPLING)) +
			geom_hline(yintercept=1, colour='grey50') +
			geom_boxplot(stat='identity') +
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_log10(expand=c(0,0), breaks=c(1/tmp, 1, tmp), labels=c(paste0('1/',as.character(tmp)), '1', as.character(tmp)), lim=range(c(1/(tmp*1.31), tmp*1.31))) +
			facet_grid(ANALYSIS~., scales='free', space='free') +
			coord_flip() +
			labs(x='', y='\ntransmissions inland->fishing sites /\n   transmissions fishing sites->inland', fill='minimum\nsequencing depth')
	ggsave(file=file.path(outfile.base,'deepseqcutoff_comparison_sink.pdf'), w=12, h=5)
	
}

RakaiFull.phylogeography.180618.flows.fishinland<- function()
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
	
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_inference.rda"
	infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_180618_cl25_d50_prior23_min30_phylogeography_core_inference_seqExclART.rda"
	outfile.base		<- gsub('.rda$','',gsub('_inference','',infile.inference))
	load(infile.inference)
	
	qs		<- c(0.025,0.25,0.5,0.75,0.975)
	qsn		<- c('CL','IL','M','IU','CU')
	
	#
	#	geography who infects whom matrix  between fisherfolk and others
	#	adjusted P
	tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))
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
	z		<- z[, list(P=qsn, Q=unname(quantile(value, p=qs))), by=c('variable')]
	z[, TR_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\1',variable)]
	z[, REC_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\2',variable)]		
	z		<- dcast.data.table(z, TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	z[, LABEL:= paste0(round(M*100, d=1), '%\n[',round(CL*100,d=1),'% - ',round(CU*100,d=1),'%]')]
	z[, LABEL2:= paste0(round(M*100, d=1), '% (',round(CL*100,d=1),'%-',round(CU*100,d=1),'%)')]
	setkey(z, REC_COMM_TYPE, TR_COMM_TYPE)
	z[, STAT:='joint']
	z[, DUMMY:= seq_len(nrow(z))]
	ans		<- copy(z)	
	
	#
	#	WAIFM
	#
	groups	<- c('inland','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))				
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
	
	
	
	#
	#	sources
	#
	groups	<- c('inland','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- unique(subset(ds, select=c(COMM_NUM_A, COMM_TYPE)))				
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
	
	save(ans, file=paste0(outfile.base,'_fishing_inland_results.rda'))
}

RakaiFull.phylogeography.180322.participitation.bias.figure<- function()
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
	
	infile				<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender_age.rda'
	outfile.base		<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	load(infile)
	
	set(de, de[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
	de[, PART:= de[, factor(PARTICIPATED_ANY_VISIT, levels=c(0,1), labels=c('PART_NEVER','PART_EVER'))]]	
	
	#	full interaction plot
	df	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','SEX','AGE_AT_MID_C','INMIGRANT','PART')]	
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+SEX+AGE_AT_MID_C+INMIGRANT ~ PART, value.var='N')
	set(df, df[, which(is.na(PART_NEVER))], 'PART_NEVER', 0L)
	set(df, df[, which(is.na(PART_EVER))], 'PART_EVER', 0L)
	df	<- df[, list(STAT=c('M','CL','CU'), V=as.numeric(binconf(PART_EVER, PART_EVER+PART_NEVER))), by=c('COMM_NUM','COMM_NUM_A','SEX','AGE_AT_MID_C','INMIGRANT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+SEX+AGE_AT_MID_C+INMIGRANT ~ STAT, value.var='V')
	df[, COMM_TYPE:= factor(substr(COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fishing site','inland community'))]
	set(df, NULL, 'INMIGRANT', df[, factor(INMIGRANT, levels=c(0,1), labels=c('resident','immigrated'))])
	ggplot(df, aes(x=COMM_NUM_A, y=M, colour=AGE_AT_MID_C, shape=INMIGRANT)) +
			geom_point(position=position_dodge(width=0.9)) +
			geom_errorbar(aes(ymin=CL, ymax=CU), position=position_dodge(width=0.9)) +
			facet_grid(SEX~COMM_TYPE, scales='free', space='free') +
			theme_bw()
	
	
	#	gender plot
	de[, STRAT:= SEX]
	df	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','STRAT','PART')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ PART, value.var='N')
	set(df, df[, which(is.na(PART_NEVER))], 'PART_NEVER', 0L)
	set(df, df[, which(is.na(PART_EVER))], 'PART_EVER', 0L)
	df	<- df[, list(STAT=c('M','CL','CU'), V=as.numeric(binconf(PART_EVER, PART_EVER+PART_NEVER))), by=c('COMM_NUM','COMM_NUM_A','STRAT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ STAT, value.var='V')
	df[, COMM_TYPE:= factor(substr(COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fishing site','inland community'))]
	set(df, NULL, 'STRAT', df[, factor(STRAT, levels=c('M','F'), labels=c('male','female'))])
	ggplot(df, aes(x=COMM_NUM_A, fill=STRAT, y=M)) +
			geom_bar(position=position_dodge(0.9), stat='identity') +
			geom_errorbar(aes(ymin=CL, ymax=CU), position=position_dodge(0.9), width=0.4) +
			scale_fill_manual(values=c('female'='hotpink2', 'male'='steelblue2')) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			facet_grid(~COMM_TYPE, space='free', scales='free') +
			labs(x='\ncommunity', y='participation rate\n', fill='gender') +
			theme_bw()
	ggsave(paste0(outfile.base,'participation_differences_gender_barplot_180614.pdf'), w=10, h=5, useDingbats=FALSE)
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+COMM_TYPE~STRAT, value.var=c('CL','CU','M'))
	ggplot(df, aes(x=M_male, y=M_female, colour=COMM_TYPE)) + 
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_errorbar(aes(ymin=CL_female, ymax=CU_female), alpha=0.25) +
			geom_errorbarh(aes(xmin=CL_male, xmax=CU_male), alpha=0.25) +
			geom_point() +			
			scale_x_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_colour_manual(values=c('fishing site'='#F8766D', 'inland community'='#00BFC4')) +
			theme_bw() +
			theme(legend.position = c(0.2, 0.85)) +
			labs(x='\nparticipation rate among men', 
					y='participation rate among women\n',
					colour='')
	ggsave(paste0(outfile.base,'participation_differences_gender_pointplot_180614.pdf'), w=5, h=5, useDingbats=FALSE)
	
	
	#	age plot
	de[, STRAT:= AGE_AT_MID_C]
	df	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','STRAT','PART')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ PART, value.var='N')
	set(df, df[, which(is.na(PART_NEVER))], 'PART_NEVER', 0L)
	set(df, df[, which(is.na(PART_EVER))], 'PART_EVER', 0L)
	df	<- df[, list(STAT=c('M','CL','CU'), V=as.numeric(binconf(PART_EVER, PART_EVER+PART_NEVER))), by=c('COMM_NUM','COMM_NUM_A','STRAT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ STAT, value.var='V')
	df[, COMM_TYPE:= factor(substr(COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fishing site','inland community'))]
	#set(df, NULL, 'STRAT', df[, factor(STRAT, levels=c('M','F'), labels=c('male','female'))])
	ggplot(df, aes(x=COMM_NUM_A, fill=STRAT, y=M)) +
			geom_bar(position=position_dodge(0.9), stat='identity') +
			geom_errorbar(aes(ymin=CL, ymax=CU), position=position_dodge(0.9), width=0.4) +
			scale_fill_manual(values=c('15-24'='#00D2FF', '25-34'='#0096C2','35+'='#007099')) +			
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			facet_grid(~COMM_TYPE, space='free', scales='free') +
			labs(x='\ncommunity', y='participation rate\n', fill='age group') +
			theme_bw()
	ggsave(paste0(outfile.base,'participation_differences_age_barplot_180614.pdf'), w=10, h=5, useDingbats=FALSE)
	set(df, NULL, 'STRAT', df[, gsub('\\+|-','',as.character(STRAT))])
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+COMM_TYPE~STRAT, value.var=c('CL','CU','M'))
	ggplot(df, aes(x=M_1524, y=M_2534, colour=COMM_TYPE)) + 
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_errorbar(aes(ymin=CL_2534, ymax=CU_2534), alpha=0.25) +
			geom_errorbarh(aes(xmin=CL_1524, xmax=CU_1524), alpha=0.25) +
			geom_point() +			
			scale_x_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_colour_manual(values=c('fishing site'='#F8766D', 'inland community'='#00BFC4')) +
			theme_bw() +
			theme(legend.position = c(0.2, 0.85)) +
			labs(x='\nparticipation rate among individuals aged 15-24', 
					y='participation rate among individuals aged 25-34\n',
					colour='')
	ggsave(paste0(outfile.base,'participation_differences_age_pointplot_180614.pdf'), w=5, h=5, useDingbats=FALSE)
	ggplot(df, aes(x=M_1524, y=M_35, colour=COMM_TYPE)) + 
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_errorbar(aes(ymin=CL_35, ymax=CU_35), alpha=0.25) +
			geom_errorbarh(aes(xmin=CL_1524, xmax=CU_1524), alpha=0.25) +
			geom_point() +			
			scale_x_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_colour_manual(values=c('fishing site'='#F8766D', 'inland community'='#00BFC4')) +
			theme_bw() +
			theme(legend.position = c(0.2, 0.85)) +
			labs(x='\nparticipation rate among individuals aged 15-24', 
					y='participation rate among individuals aged 35+\n',
					colour='')
	ggsave(paste0(outfile.base,'participation_differences_age_pointplot2_180614.pdf'), w=5, h=5, useDingbats=FALSE)
	ggplot(df, aes(x=M_2534, y=M_35, colour=COMM_TYPE)) + 
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_errorbar(aes(ymin=CL_35, ymax=CU_35), alpha=0.25) +
			geom_errorbarh(aes(xmin=CL_2534, xmax=CU_2534), alpha=0.25) +
			geom_point() +			
			scale_x_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_colour_manual(values=c('fishing site'='#F8766D', 'inland community'='#00BFC4')) +
			theme_bw() +
			theme(legend.position = c(0.2, 0.85)) +
			labs(x='\nparticipation rate among individuals aged 25-34', 
					y='participation rate among individuals aged 35+\n',
					colour='')
	ggsave(paste0(outfile.base,'participation_differences_age_pointplot3_180614.pdf'), w=5, h=5, useDingbats=FALSE)
	
	
	#	inmigration plot
	de[, STRAT:= INMIGRANT]
	df	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','STRAT','PART')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ PART, value.var='N')
	set(df, df[, which(is.na(PART_NEVER))], 'PART_NEVER', 0L)
	set(df, df[, which(is.na(PART_EVER))], 'PART_EVER', 0L)
	df	<- df[, list(STAT=c('M','CL','CU'), V=as.numeric(binconf(PART_EVER, PART_EVER+PART_NEVER))), by=c('COMM_NUM','COMM_NUM_A','STRAT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ STAT, value.var='V')
	df[, COMM_TYPE:= factor(substr(COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fishing site','inland community'))]
	set(df, NULL, 'STRAT2', df[, factor(STRAT, levels=c(0,1), labels=c('resident','immigrated\nin two years\nbefore visit'))])
	ggplot(df, aes(x=COMM_NUM_A, fill=STRAT2, y=M)) +
			geom_bar(position=position_dodge(0.9), stat='identity') +
			geom_errorbar(aes(ymin=CL, ymax=CU), position=position_dodge(0.9), width=0.4) +
			#scale_fill_manual(values=c('female'='hotpink2', 'male'='steelblue2')) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			facet_grid(~COMM_TYPE, space='free', scales='free') +
			labs(x='\ncommunity', y='participation rate\n', fill='migration status') +
			theme_bw()
	ggsave(paste0(outfile.base,'participation_differences_inmigrants_barplot_180614.pdf'), w=10, h=5, useDingbats=FALSE)
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+COMM_TYPE~STRAT, value.var=c('CL','CU','M'))
	ggplot(df, aes(x=M_0, y=M_1, colour=COMM_TYPE)) + 
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_errorbar(aes(ymin=CL_1, ymax=CU_1), alpha=0.25) +
			geom_errorbarh(aes(xmin=CL_0, xmax=CU_0), alpha=0.25) +
			geom_point() +			
			scale_x_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_colour_manual(values=c('fishing site'='#F8766D', 'inland community'='#00BFC4')) +
			theme_bw() +
			theme(legend.position = c(0.2, 0.85)) +
			labs(x='\nparticipation rate among\nresident individuals', 
					y='participation rate among\nindividuals who inmigrated in two years before visit\n',
					colour='')
	ggsave(paste0(outfile.base,'participation_differences_inmigrants_pointplot_180614.pdf'), w=5, h=5, useDingbats=FALSE)
}

RakaiFull.phylogeography.180322.sequencing.bias.figure<- function()
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
	infile			<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender_age.rda'	
	outfile.base	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	load(infile)
	
	set(de, de[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
	de	<- subset(de, HIV_1517==1)
	set(de, NULL, 'MIN_PNG_OUTPUT', de[, factor(MIN_PNG_OUTPUT, levels=c(0,1), labels=c('SEQ_NO','SEQ_YES'))])
	
	#	gender plot
	de[, STRAT:= SEX]
	df	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','STRAT','MIN_PNG_OUTPUT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ MIN_PNG_OUTPUT, value.var='N')
	set(df, df[, which(is.na(SEQ_NO))], 'SEQ_NO', 0L)
	set(df, df[, which(is.na(SEQ_YES))], 'SEQ_YES', 0L)
	df	<- df[, list(STAT=c('M','CL','CU'), V=as.numeric(binconf(SEQ_YES, SEQ_YES+SEQ_NO))), by=c('COMM_NUM','COMM_NUM_A','STRAT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ STAT, value.var='V')
	df[, COMM_TYPE:= factor(substr(COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fishing site','inland community'))]
	set(df, NULL, 'STRAT', df[, factor(STRAT, levels=c('M','F'), labels=c('male','female'))])
	tmp	<- df[, which(M==0)]
	set(df, tmp, c('CL','CU'), 0)
	ggplot(df, aes(x=COMM_NUM_A, fill=STRAT, y=M)) +
			geom_bar(position=position_dodge(0.9), stat='identity') +
			geom_errorbar(aes(ymin=CL, ymax=CU), position=position_dodge(0.9), width=0.4) +
			scale_fill_manual(values=c('female'='hotpink2', 'male'='steelblue2')) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			facet_grid(~COMM_TYPE, space='free', scales='free') +
			labs(x='\ncommunity', y='sequencing rate\n', fill='gender') +
			theme_bw()
	ggsave(paste0(outfile.base,'sequencing_differences_gender_barplot_180614.pdf'), w=10, h=5, useDingbats=FALSE)
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+COMM_TYPE~STRAT, value.var=c('CL','CU','M'))
	ggplot(df, aes(x=M_male, y=M_female, colour=COMM_TYPE)) + 
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_errorbar(aes(ymin=CL_female, ymax=CU_female), alpha=0.25) +
			geom_errorbarh(aes(xmin=CL_male, xmax=CU_male), alpha=0.25) +
			geom_point() +			
			scale_x_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_colour_manual(values=c('fishing site'='#F8766D', 'inland community'='#00BFC4')) +
			theme_bw() +
			theme(legend.position = c(0.2, 0.85)) +
			labs(x='\nsequencing rate among men', 
					y='sequencing rate among women\n',
					colour='')
	ggsave(paste0(outfile.base,'sequencing_differences_gender_pointplot_180614.pdf'), w=5, h=5, useDingbats=FALSE)
	
	
	#	age plot
	de[, STRAT:= AGE_AT_MID_C]
	df	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','STRAT','MIN_PNG_OUTPUT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ MIN_PNG_OUTPUT, value.var='N')
	set(df, df[, which(is.na(SEQ_NO))], 'SEQ_NO', 0L)
	set(df, df[, which(is.na(SEQ_YES))], 'SEQ_YES', 0L)
	df	<- df[, list(STAT=c('M','CL','CU'), V=as.numeric(binconf(SEQ_YES, SEQ_YES+SEQ_NO))), by=c('COMM_NUM','COMM_NUM_A','STRAT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ STAT, value.var='V')
	df[, COMM_TYPE:= factor(substr(COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fishing site','inland community'))]
	tmp	<- df[, which(M==0)]
	set(df, tmp, c('CL','CU'), 0)	
	#set(df, NULL, 'STRAT', df[, factor(STRAT, levels=c('M','F'), labels=c('male','female'))])
	ggplot(df, aes(x=COMM_NUM_A, fill=STRAT, y=M)) +
			geom_bar(position=position_dodge(0.9), stat='identity') +
			geom_errorbar(aes(ymin=CL, ymax=CU), position=position_dodge(0.9), width=0.4) +
			scale_fill_manual(values=c('15-24'='#00D2FF', '25-34'='#0096C2','35+'='#007099')) +			
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			facet_grid(~COMM_TYPE, space='free', scales='free') +
			labs(x='\ncommunity', y='sequencing rate\n', fill='age group') +
			theme_bw()
	ggsave(paste0(outfile.base,'sequencing_differences_age_barplot_180614.pdf'), w=10, h=5, useDingbats=FALSE)
	set(df, NULL, 'STRAT', df[, gsub('\\+|-','',as.character(STRAT))])
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+COMM_TYPE~STRAT, value.var=c('CL','CU','M'))
	ggplot(df, aes(x=M_1524, y=M_2534, colour=COMM_TYPE)) + 
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_errorbar(aes(ymin=CL_2534, ymax=CU_2534), alpha=0.25) +
			geom_errorbarh(aes(xmin=CL_1524, xmax=CU_1524), alpha=0.25) +
			geom_point() +			
			scale_x_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_colour_manual(values=c('fishing site'='#F8766D', 'inland community'='#00BFC4')) +
			theme_bw() +
			theme(legend.position = c(0.2, 0.85)) +
			labs(x='\nsequencing rate among individuals aged 15-24', 
					y='sequencing rate among individuals aged 25-34\n',
					colour='')
	ggsave(paste0(outfile.base,'sequencing_differences_age_pointplot_180614.pdf'), w=5, h=5, useDingbats=FALSE)
	ggplot(df, aes(x=M_1524, y=M_35, colour=COMM_TYPE)) + 
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_errorbar(aes(ymin=CL_35, ymax=CU_35), alpha=0.25) +
			geom_errorbarh(aes(xmin=CL_1524, xmax=CU_1524), alpha=0.25) +
			geom_point() +			
			scale_x_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_colour_manual(values=c('fishing site'='#F8766D', 'inland community'='#00BFC4')) +
			theme_bw() +
			theme(legend.position = c(0.2, 0.85)) +
			labs(x='\nsequencing rate among individuals aged 15-24', 
					y='sequencing rate among individuals aged 35+\n',
					colour='')
	ggsave(paste0(outfile.base,'sequencing_differences_age_pointplot2_180614.pdf'), w=5, h=5, useDingbats=FALSE)
	ggplot(df, aes(x=M_2534, y=M_35, colour=COMM_TYPE)) + 
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_errorbar(aes(ymin=CL_35, ymax=CU_35), alpha=0.25) +
			geom_errorbarh(aes(xmin=CL_2534, xmax=CU_2534), alpha=0.25) +
			geom_point() +			
			scale_x_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_colour_manual(values=c('fishing site'='#F8766D', 'inland community'='#00BFC4')) +
			theme_bw() +
			theme(legend.position = c(0.2, 0.85)) +
			labs(x='\nsequencing rate among individuals aged 25-34', 
					y='sequencing rate among individuals aged 35+\n',
					colour='')
	ggsave(paste0(outfile.base,'sequencing_differences_age_pointplot3_180614.pdf'), w=5, h=5, useDingbats=FALSE)
	
	
	#	inmigration plot
	de[, STRAT:= INMIGRANT]
	df	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','STRAT','MIN_PNG_OUTPUT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ MIN_PNG_OUTPUT, value.var='N')
	set(df, df[, which(is.na(SEQ_NO))], 'SEQ_NO', 0L)
	set(df, df[, which(is.na(SEQ_YES))], 'SEQ_YES', 0L)
	df	<- df[, list(STAT=c('M','CL','CU'), V=as.numeric(binconf(SEQ_YES, SEQ_YES+SEQ_NO))), by=c('COMM_NUM','COMM_NUM_A','STRAT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ STAT, value.var='V')
	df[, COMM_TYPE:= factor(substr(COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fishing site','inland community'))]
	set(df, NULL, 'STRAT2', df[, factor(STRAT, levels=c(0,1), labels=c('resident','immigrated\nin two years\nbefore visit'))])
	tmp	<- df[, which(M==0)]
	set(df, tmp, c('CL','CU'), 0)	
	ggplot(df, aes(x=COMM_NUM_A, fill=STRAT2, y=M)) +
			geom_bar(position=position_dodge(0.9), stat='identity') +
			geom_errorbar(aes(ymin=CL, ymax=CU), position=position_dodge(0.9), width=0.4) +
			#scale_fill_manual(values=c('female'='hotpink2', 'male'='steelblue2')) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			facet_grid(~COMM_TYPE, space='free', scales='free') +
			labs(x='\ncommunity', y='sequencing rate\n', fill='migration status') +
			theme_bw()
	ggsave(paste0(outfile.base,'sequencing_differences_inmigrants_barplot_180614.pdf'), w=10, h=5, useDingbats=FALSE)
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+COMM_TYPE~STRAT, value.var=c('CL','CU','M'))
	ggplot(df, aes(x=M_0, y=M_1, colour=COMM_TYPE)) + 
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_errorbar(aes(ymin=CL_1, ymax=CU_1), alpha=0.25) +
			geom_errorbarh(aes(xmin=CL_0, xmax=CU_0), alpha=0.25) +
			geom_point() +			
			scale_x_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_colour_manual(values=c('fishing site'='#F8766D', 'inland community'='#00BFC4')) +
			theme_bw() +
			theme(legend.position = c(0.2, 0.85)) +
			labs(x='\nsequencing rate among\nresident individuals', 
					y='sequencing rate among\nindividuals who inmigrated in two years before visit\n',
					colour='')
	ggsave(paste0(outfile.base,'sequencing_differences_inmigrants_pointplot_180614.pdf'), w=5, h=5, useDingbats=FALSE)
	
	
	#	full interaction plot
	df	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','SEX','AGE_AT_MID_C','INMIGRANT','MIN_PNG_OUTPUT')]	
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+SEX+AGE_AT_MID_C+INMIGRANT ~ MIN_PNG_OUTPUT, value.var='N')
	set(df, df[, which(is.na(SEQ_NO))], 'SEQ_NO', 0L)
	set(df, df[, which(is.na(SEQ_YES))], 'SEQ_YES', 0L)
	df	<- df[, list(STAT=c('M','CL','CU'), V=as.numeric(binconf(SEQ_YES, SEQ_YES+SEQ_NO))), by=c('COMM_NUM','COMM_NUM_A','SEX','AGE_AT_MID_C','INMIGRANT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+SEX+AGE_AT_MID_C+INMIGRANT ~ STAT, value.var='V')
	df[, COMM_TYPE:= factor(substr(COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fishing site','inland community'))]
	set(df, NULL, 'INMIGRANT', df[, factor(INMIGRANT, levels=c(0,1), labels=c('resident','immigrated'))])
	tmp	<- df[, which(M==0)]
	set(df, tmp, c('CL','CU'), 0)	
	ggplot(df, aes(x=COMM_NUM_A, y=M, colour=AGE_AT_MID_C, shape=INMIGRANT)) +
			geom_point(position=position_dodge(width=0.9)) +
			geom_errorbar(aes(ymin=CL, ymax=CU), position=position_dodge(width=0.9)) +
			facet_grid(SEX~COMM_TYPE, scales='free', space='free') +
			theme_bw()
	ggsave(paste0(outfile.base,'sequencing_differences_interactions_barplot_180614.pdf'), w=14, h=7, useDingbats=FALSE)
}


RakaiFull.phylogeography.180322.sequencing.bias.exclART.figure<- function()
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
	infile			<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender_age.rda'	
	outfile.base	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	load(infile)
	
	set(de, de[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
	de	<- subset(de, HIV_1517==1)
	de[, DUMMY:=as.integer(SELFREPORTART_AT_FIRST_VISIT!=1 | (SELFREPORTART_AT_FIRST_VISIT==1 & MIN_PNG_OUTPUT==1))]
	
	stopifnot(de[, !any(is.na(SELFREPORTART_AT_FIRST_VISIT))])
	df	<- de[, list( ART_YES= length(which(DUMMY==0)), ART_NO=length(which(DUMMY==1))  ), by='COMM_NUM_A']
	df	<- merge(df, de, by='COMM_NUM_A')
	df	<- df[, list( TYPE='not on ART (self-reported)', STAT=c('M','CL','CU'), V=as.numeric(binconf(sum(MIN_PNG_OUTPUT), round( ART_NO[1]/(ART_NO[1]+ART_YES[1]) * sum(HIV_1517))))  ), by='COMM_NUM_A']
	tmp	<- de[, list( TYPE='all infected', STAT=c('M','CL','CU'), V=as.numeric(binconf(sum(MIN_PNG_OUTPUT), sum(HIV_1517)))  ), by='COMM_NUM_A']
	df	<- rbind(df, tmp)
	df[, COMM_TYPE:= factor(substr(COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fishing site','inland community'))]
	df	<- dcast.data.table(df, COMM_TYPE+COMM_NUM_A+TYPE~STAT, value.var='V')
	tmp	<- df[, which(M==0)]
	set(df, tmp, c('CL','CU'), 0)	
	ggplot(df, aes(x=COMM_NUM_A, colour=TYPE, y=M, ymin=CL, ymax=CU)) +			
			geom_errorbar(position=position_dodge(0.9), width=0.4) +
			geom_point(position=position_dodge(0.9)) +
			scale_colour_brewer(palette='Paired') +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			facet_grid(~COMM_TYPE, space='free', scales='free') +
			labs(x='\ncommunity', y='sequencing rate\n', colour='denominator') +
			theme_bw() +
			theme(legend.position = 'bottom') 
	ggsave(paste0(outfile.base,'sequencingExclART_differences_community_180614.pdf'), w=10, h=6, useDingbats=FALSE)
	
	
	de	<- subset(de, is.na(SELFREPORTART_AT_FIRST_VISIT) | SELFREPORTART_AT_FIRST_VISIT!=1)
	set(de, NULL, 'MIN_PNG_OUTPUT', de[, factor(MIN_PNG_OUTPUT, levels=c(0,1), labels=c('SEQ_NO','SEQ_YES'))])
	
	#	gender plot
	de[, STRAT:= SEX]
	df	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','STRAT','MIN_PNG_OUTPUT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ MIN_PNG_OUTPUT, value.var='N')
	set(df, df[, which(is.na(SEQ_NO))], 'SEQ_NO', 0L)
	set(df, df[, which(is.na(SEQ_YES))], 'SEQ_YES', 0L)
	df	<- df[, list(STAT=c('M','CL','CU'), V=as.numeric(binconf(SEQ_YES, SEQ_YES+SEQ_NO))), by=c('COMM_NUM','COMM_NUM_A','STRAT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ STAT, value.var='V')
	df[, COMM_TYPE:= factor(substr(COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fishing site','inland community'))]
	set(df, NULL, 'STRAT', df[, factor(STRAT, levels=c('M','F'), labels=c('male','female'))])
	tmp	<- df[, which(M==0)]
	set(df, tmp, c('CL','CU'), 0)
	ggplot(df, aes(x=COMM_NUM_A, fill=STRAT, y=M)) +
			geom_bar(position=position_dodge(0.9), stat='identity') +
			geom_errorbar(aes(ymin=CL, ymax=CU), position=position_dodge(0.9), width=0.4) +
			scale_fill_manual(values=c('female'='hotpink2', 'male'='steelblue2')) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			facet_grid(~COMM_TYPE, space='free', scales='free') +
			labs(x='\ncommunity', y='sequencing rate\n', fill='gender') +
			theme_bw()
	ggsave(paste0(outfile.base,'sequencingExclART_differences_gender_barplot_180614.pdf'), w=10, h=5, useDingbats=FALSE)
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+COMM_TYPE~STRAT, value.var=c('CL','CU','M'))
	ggplot(df, aes(x=M_male, y=M_female, colour=COMM_TYPE)) + 
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_errorbar(aes(ymin=CL_female, ymax=CU_female), alpha=0.25) +
			geom_errorbarh(aes(xmin=CL_male, xmax=CU_male), alpha=0.25) +
			geom_point() +			
			scale_x_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_colour_manual(values=c('fishing site'='#F8766D', 'inland community'='#00BFC4')) +
			theme_bw() +
			theme(legend.position = c(0.2, 0.85)) +
			labs(x='\nsequencing rate among men', 
					y='sequencing rate among women\n',
					colour='')
	ggsave(paste0(outfile.base,'sequencingExclART_differences_gender_pointplot_180614.pdf'), w=5, h=5, useDingbats=FALSE)
	
	
	#	age plot
	de[, STRAT:= AGE_AT_MID_C]
	df	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','STRAT','MIN_PNG_OUTPUT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ MIN_PNG_OUTPUT, value.var='N')
	set(df, df[, which(is.na(SEQ_NO))], 'SEQ_NO', 0L)
	set(df, df[, which(is.na(SEQ_YES))], 'SEQ_YES', 0L)
	df	<- df[, list(STAT=c('M','CL','CU'), V=as.numeric(binconf(SEQ_YES, SEQ_YES+SEQ_NO))), by=c('COMM_NUM','COMM_NUM_A','STRAT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ STAT, value.var='V')
	df[, COMM_TYPE:= factor(substr(COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fishing site','inland community'))]
	tmp	<- df[, which(M==0)]
	set(df, tmp, c('CL','CU'), 0)	
	#set(df, NULL, 'STRAT', df[, factor(STRAT, levels=c('M','F'), labels=c('male','female'))])
	ggplot(df, aes(x=COMM_NUM_A, fill=STRAT, y=M)) +
			geom_bar(position=position_dodge(0.9), stat='identity') +
			geom_errorbar(aes(ymin=CL, ymax=CU), position=position_dodge(0.9), width=0.4) +
			scale_fill_manual(values=c('15-24'='#00D2FF', '25-34'='#0096C2','35+'='#007099')) +			
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			facet_grid(~COMM_TYPE, space='free', scales='free') +
			labs(x='\ncommunity', y='sequencing rate\n', fill='age group') +
			theme_bw()
	ggsave(paste0(outfile.base,'sequencingExclART_differences_age_barplot_180614.pdf'), w=10, h=5, useDingbats=FALSE)
	set(df, NULL, 'STRAT', df[, gsub('\\+|-','',as.character(STRAT))])
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+COMM_TYPE~STRAT, value.var=c('CL','CU','M'))
	ggplot(df, aes(x=M_1524, y=M_2534, colour=COMM_TYPE)) + 
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_errorbar(aes(ymin=CL_2534, ymax=CU_2534), alpha=0.25) +
			geom_errorbarh(aes(xmin=CL_1524, xmax=CU_1524), alpha=0.25) +
			geom_point() +			
			scale_x_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_colour_manual(values=c('fishing site'='#F8766D', 'inland community'='#00BFC4')) +
			theme_bw() +
			theme(legend.position = c(0.2, 0.85)) +
			labs(x='\nsequencing rate among individuals aged 15-24', 
					y='sequencing rate among individuals aged 25-34\n',
					colour='')
	ggsave(paste0(outfile.base,'sequencingExclART_differences_age_pointplot_180614.pdf'), w=5, h=5, useDingbats=FALSE)
	ggplot(df, aes(x=M_1524, y=M_35, colour=COMM_TYPE)) + 
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_errorbar(aes(ymin=CL_35, ymax=CU_35), alpha=0.25) +
			geom_errorbarh(aes(xmin=CL_1524, xmax=CU_1524), alpha=0.25) +
			geom_point() +			
			scale_x_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_colour_manual(values=c('fishing site'='#F8766D', 'inland community'='#00BFC4')) +
			theme_bw() +
			theme(legend.position = c(0.2, 0.85)) +
			labs(x='\nsequencing rate among individuals aged 15-24', 
					y='sequencing rate among individuals aged 35+\n',
					colour='')
	ggsave(paste0(outfile.base,'sequencingExclART_differences_age_pointplot2_180614.pdf'), w=5, h=5, useDingbats=FALSE)
	ggplot(df, aes(x=M_2534, y=M_35, colour=COMM_TYPE)) + 
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_errorbar(aes(ymin=CL_35, ymax=CU_35), alpha=0.25) +
			geom_errorbarh(aes(xmin=CL_2534, xmax=CU_2534), alpha=0.25) +
			geom_point() +			
			scale_x_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_colour_manual(values=c('fishing site'='#F8766D', 'inland community'='#00BFC4')) +
			theme_bw() +
			theme(legend.position = c(0.2, 0.85)) +
			labs(x='\nsequencing rate among individuals aged 25-34', 
					y='sequencing rate among individuals aged 35+\n',
					colour='')
	ggsave(paste0(outfile.base,'sequencingExclART_differences_age_pointplot3_180614.pdf'), w=5, h=5, useDingbats=FALSE)
	
	
	#	inmigration plot
	de[, STRAT:= INMIGRANT]
	df	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','STRAT','MIN_PNG_OUTPUT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ MIN_PNG_OUTPUT, value.var='N')
	set(df, df[, which(is.na(SEQ_NO))], 'SEQ_NO', 0L)
	set(df, df[, which(is.na(SEQ_YES))], 'SEQ_YES', 0L)
	df	<- df[, list(STAT=c('M','CL','CU'), V=as.numeric(binconf(SEQ_YES, SEQ_YES+SEQ_NO))), by=c('COMM_NUM','COMM_NUM_A','STRAT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+STRAT ~ STAT, value.var='V')
	df[, COMM_TYPE:= factor(substr(COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fishing site','inland community'))]
	set(df, NULL, 'STRAT2', df[, factor(STRAT, levels=c(0,1), labels=c('resident','immigrated\nin two years\nbefore visit'))])
	tmp	<- df[, which(M==0)]
	set(df, tmp, c('CL','CU'), 0)	
	ggplot(df, aes(x=COMM_NUM_A, fill=STRAT2, y=M)) +
			geom_bar(position=position_dodge(0.9), stat='identity') +
			geom_errorbar(aes(ymin=CL, ymax=CU), position=position_dodge(0.9), width=0.4) +
			#scale_fill_manual(values=c('female'='hotpink2', 'male'='steelblue2')) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			facet_grid(~COMM_TYPE, space='free', scales='free') +
			labs(x='\ncommunity', y='sequencing rate\n', fill='migration status') +
			theme_bw()
	ggsave(paste0(outfile.base,'sequencingExclART_differences_inmigrants_barplot_180614.pdf'), w=10, h=5, useDingbats=FALSE)
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+COMM_TYPE~STRAT, value.var=c('CL','CU','M'))
	ggplot(df, aes(x=M_0, y=M_1, colour=COMM_TYPE)) + 
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_errorbar(aes(ymin=CL_1, ymax=CU_1), alpha=0.25) +
			geom_errorbarh(aes(xmin=CL_0, xmax=CU_0), alpha=0.25) +
			geom_point() +			
			scale_x_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			scale_colour_manual(values=c('fishing site'='#F8766D', 'inland community'='#00BFC4')) +
			theme_bw() +
			theme(legend.position = c(0.2, 0.85)) +
			labs(x='\nsequencing rate among\nresident individuals', 
					y='sequencing rate among\nindividuals who inmigrated in two years before visit\n',
					colour='')
	ggsave(paste0(outfile.base,'sequencingExclART_differences_inmigrants_pointplot_180614.pdf'), w=5, h=5, useDingbats=FALSE)
	
	
	#	full interaction plot
	df	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','SEX','AGE_AT_MID_C','INMIGRANT','MIN_PNG_OUTPUT')]	
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+SEX+AGE_AT_MID_C+INMIGRANT ~ MIN_PNG_OUTPUT, value.var='N')
	set(df, df[, which(is.na(SEQ_NO))], 'SEQ_NO', 0L)
	set(df, df[, which(is.na(SEQ_YES))], 'SEQ_YES', 0L)
	df	<- df[, list(STAT=c('M','CL','CU'), V=as.numeric(binconf(SEQ_YES, SEQ_YES+SEQ_NO))), by=c('COMM_NUM','COMM_NUM_A','SEX','AGE_AT_MID_C','INMIGRANT')]
	df	<- dcast.data.table(df, COMM_NUM+COMM_NUM_A+SEX+AGE_AT_MID_C+INMIGRANT ~ STAT, value.var='V')
	df[, COMM_TYPE:= factor(substr(COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fishing site','inland community'))]
	set(df, NULL, 'INMIGRANT', df[, factor(INMIGRANT, levels=c(0,1), labels=c('resident','immigrated'))])
	tmp	<- df[, which(M==0)]
	set(df, tmp, c('CL','CU'), 0)	
	ggplot(df, aes(x=COMM_NUM_A, y=M, colour=AGE_AT_MID_C, shape=INMIGRANT)) +
			geom_point(position=position_dodge(width=0.9)) +
			geom_errorbar(aes(ymin=CL, ymax=CU), position=position_dodge(width=0.9)) +
			facet_grid(SEX~COMM_TYPE, scales='free', space='free') +
			theme_bw()
	ggsave(paste0(outfile.base,'sequencingExclART_differences_interactions_barplot_180614.pdf'), w=14, h=7, useDingbats=FALSE)
}


RakaiFull.phylogeography.180322.participitation.bias.logisticmodel<- function()
{
	require(data.table)
	require(rethinking)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	#	load des which contains participation and seq counts by comm and gender
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender_age.rda'
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	load(infile)	
	set(de, de[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
	de		<- subset(de, !is.na(AGE_AT_MID_C))
		
	#	quick logistic regression on participation
	dg		<- subset(de, select=c(PARTICIPATED_ANY_VISIT, PERM_ID, COMM_NUM_A, AGE_AT_MID_C, INMIGRANT, SEX))
	#	binarize age, sex
	dg[, AGE_YOUNG:= as.integer(AGE_AT_MID_C=='15-24')]
	dg[, AGE_MID:= as.integer(AGE_AT_MID_C=='25-34')]
	dg[, MALE:= as.integer(SEX=='M')]
	#	binarize community type
	dg[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	dg[, COMM_TYPE_T:= as.integer(substr(COMM_NUM_A,1,1)=='t')]
	#	vanilla community IDs
	tmp		<- data.table(COMM_NUM_A= sort(unique(dg$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(dg$COMM_NUM_A)))
	dg		<- merge(dg, tmp, by='COMM_NUM_A')
	#	aggregate
	dg		<- dg[, list(ELIGIBLE=length(PERM_ID), PART=length(which(PARTICIPATED_ANY_VISIT==1))), by=c('AGE_AT_MID_C','INMIGRANT','SEX','AGE_YOUNG','AGE_MID','MALE','COMM_TYPE_F','COMM_TYPE_T','COMM_NUM_A','COMM_NUM_B')]	
	mp1 	<- map2stan(
			alist(
					PART ~ dbinom(ELIGIBLE, p_part),
					logit(p_part) <- a + comm[COMM_NUM_B] + trading*COMM_TYPE_T + fishing*COMM_TYPE_F + 
										inmigrant*INMIGRANT + male*MALE + young*AGE_YOUNG + midage*AGE_MID,
					a ~ dnorm(0, 100),
					comm ~ dnorm(0, sig_comm),
					sig_comm ~ dcauchy(0,1),
					c(trading, fishing, inmigrant, male, young, midage) ~ dnorm(0,10)
			),
			data=as.data.frame(dg), 
			start=list(	a=0, comm=rep(0,39), sig_comm=1, trading=0, fishing=0, inmigrant=0, male=0, young=0, midage=0),
			warmup=5e2, iter=5e3, chains=1, cores=4
		)		
	plot(mp1)
	plot( precis(mp1, depth=2, prob=0.95) )
	#	posterior check to see if this makes sense
	sims 	<- sim(mp1, n=1e3)
	tmp		<- apply(sims, 2, median)
	dpp 	<- apply(sims, 2, PI, prob=0.95)
	dpp		<- rbind(tmp, dpp)
	rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
	dpp		<- as.data.table(t(dpp))
	dpp		<- cbind(dg, dpp)	
	
	dpp[, mean(PART<predicted_obs_l95 | PART>predicted_obs_u95)]
	#	0.168 -- that s quite high
	dpp[, OFF:= as.integer(PART<predicted_obs_l95 | PART>predicted_obs_u95)]
	dpp[, list(OFF_P= mean(OFF)), by=c('AGE_AT_MID_C','INMIGRANT')]
	#	2:        15-24         0 0.2763158 --> young residents are most frequently off
	dpp[, list(OFF_P= mean(OFF), OFF_N=sum(OFF)), by=c('AGE_AT_MID_C','SEX','INMIGRANT')]
	#	 3:        15-24   F    0 0.44736842 --> especially young female residents are off
	#	           25-34   F    0 0.28947368 --> next one up is midaged female residents

	subset(dpp, PART<predicted_obs_l95 | PART>predicted_obs_u95)
	
	ggplot(dpp, aes(x=COMM_NUM_A)) +
			facet_grid(SEX+AGE_YOUNG+AGE_MID+INMIGRANT~COMM_TYPE_F, scales='free', space='free') +			
			geom_point(aes(y=predicted_obs_median), colour='grey50') +
			geom_errorbar(aes(ymin=predicted_obs_l95, ymax=predicted_obs_u95), colour='grey50') + 
			geom_point(aes(y=PART), colour='red') +
			theme_bw()
	ggsave(file=paste0(outfile.base,'sampling_differences_inmigrants_180322_mp1_pp.pdf'), w=16, h=32)
	#	

	
	mp2 	<- map2stan(
			alist(
					PART ~ dbinom(ELIGIBLE, p_part),
					logit(p_part) <- a + comm[COMM_NUM_B] + trading*COMM_TYPE_T + fishing*COMM_TYPE_F + 
							inmigrant*INMIGRANT + inmigrant_young*INMIGRANT*AGE_YOUNG + 
							male*MALE + 
							young_male*AGE_YOUNG*MALE + young_female*AGE_YOUNG*(1-MALE) +
							midage*AGE_MID,
					a ~ dnorm(0, 100),
					comm ~ dnorm(0, sig_comm),
					sig_comm ~ dcauchy(0,1),
					c(trading, fishing, inmigrant, inmigrant_young, male, young_male, young_female, midage) ~ dnorm(0,10)
			),
			data=as.data.frame(dg), 
			start=list(	a=0, comm=rep(0,39), sig_comm=1, trading=0, fishing=0, inmigrant=0, inmigrant_young=0, male=0, young_male=0, young_female=0, midage=0),
			warmup=5e2, iter=5e3, chains=1, cores=4
		)
	plot( precis(mp2, depth=2, prob=0.95) )
	#	posterior check to see if this makes sense
	sims 	<- sim(mp2, n=1e3)
	tmp		<- apply(sims, 2, median)
	dpp 	<- apply(sims, 2, PI, prob=0.95)
	dpp		<- rbind(tmp, dpp)
	rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
	dpp		<- as.data.table(t(dpp))
	dpp		<- cbind(dg, dpp)		
	dpp[, mean(PART<predicted_obs_l95 | PART>predicted_obs_u95)]
	#	0.1074561 -- better

	dpp[, OFF:= as.integer(PART<predicted_obs_l95 | PART>predicted_obs_u95)]
	dpp[, list(OFF_P= mean(OFF), OFF_N=sum(OFF)), by=c('AGE_AT_MID_C','SEX','INMIGRANT')]
	
	set(dpp, NULL, 'AGE_YOUNG', dpp[,paste0('young_',AGE_YOUNG)])
	set(dpp, NULL, 'AGE_MID', dpp[,paste0('mid_',AGE_MID)])
	set(dpp, NULL, 'INMIGRANT', dpp[,paste0('inm_',INMIGRANT)])
	set(dpp, NULL, 'COMM_TYPE_F', dpp[,paste0('fish_',COMM_TYPE_F)])
	ggplot(dpp, aes(x=COMM_NUM_A)) +
			facet_grid(SEX+AGE_YOUNG+AGE_MID+INMIGRANT~COMM_TYPE_F, scales='free', space='free') +			
			geom_point(aes(y=predicted_obs_median), colour='grey50') +
			geom_errorbar(aes(ymin=predicted_obs_l95, ymax=predicted_obs_u95), colour='grey50') + 
			geom_point(aes(y=PART, colour=factor(OFF))) +
			theme_bw()
	ggsave(file=paste0(outfile.base,'sampling_differences_inmigrants_180322_mp2_pp.pdf'), w=16, h=16)
	
	
	mp3 	<- map2stan(
			alist(
					PART ~ dbetabinom(ELIGIBLE, p_part, dispersion),
					logit(p_part) <- a + comm[COMM_NUM_B] + trading*COMM_TYPE_T + fishing*COMM_TYPE_F + 
							inmigrant*INMIGRANT + inmigrant_young*INMIGRANT*AGE_YOUNG + 
							male*MALE + 
							young_male*AGE_YOUNG*MALE + young_female*AGE_YOUNG*(1-MALE) +
							midage*AGE_MID,
					a ~ dnorm(0, 100),
					comm ~ dnorm(0, sig_comm),
					sig_comm ~ dcauchy(0,1),
					dispersion ~ dexp(1),
					c(trading, fishing, inmigrant, inmigrant_young, male, young_male, young_female, midage) ~ dnorm(0,10)
			),
			data=as.data.frame(dg), 
			start=list(	a=0, dispersion=1, comm=rep(0,39), sig_comm=1, trading=0, fishing=0, inmigrant=0, inmigrant_young=0, male=0, young_male=0, young_female=0, midage=0),
			warmup=5e2, iter=5e3, chains=1, cores=4
	)
	plot( precis(mp3, depth=2, prob=0.95) )
	#	posterior check to see if this makes sense
	sims 	<- sim(mp3, n=1e3)
	tmp		<- apply(sims, 2, median)
	dpp 	<- apply(sims, 2, PI, prob=0.95)
	dpp		<- rbind(tmp, dpp)
	rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
	dpp		<- as.data.table(t(dpp))
	dpp		<- cbind(dg, dpp)		
	dpp[, mean(PART<predicted_obs_l95 | PART>predicted_obs_u95)]
	#	0.01754386
	dpp[, sum(PART<predicted_obs_l95 | PART>predicted_obs_u95)]	
	#	9 / 456

	dpp[, OFF:= as.integer(PART<predicted_obs_l95 | PART>predicted_obs_u95)]
	set(dpp, NULL, 'AGE_YOUNG', dpp[,paste0('young_',AGE_YOUNG)])
	set(dpp, NULL, 'AGE_MID', dpp[,paste0('mid_',AGE_MID)])
	set(dpp, NULL, 'INMIGRANT', dpp[,paste0('inm_',INMIGRANT)])
	set(dpp, NULL, 'COMM_TYPE_F', dpp[,paste0('fish_',COMM_TYPE_F)])
	ggplot(dpp, aes(x=COMM_NUM_A)) +
			facet_grid(SEX+AGE_YOUNG+AGE_MID+INMIGRANT~COMM_TYPE_F, scales='free', space='free') +			
			geom_point(aes(y=predicted_obs_median), colour='grey50') +
			geom_errorbar(aes(ymin=predicted_obs_l95, ymax=predicted_obs_u95), colour='grey50') + 
			geom_point(aes(y=PART, colour=factor(OFF))) +
			theme_bw()
	ggsave(file=paste0(outfile.base,'sampling_differences_inmigrants_180322_mp3_pp.pdf'), w=16, h=16)

	
	
	
	
	mp4 	<- map2stan(
			alist(
					PART ~ dbetabinom(ELIGIBLE, p_part, dispersion),
					logit(p_part) <- a + comm[COMM_NUM_B] + trading*COMM_TYPE_T + fishing*COMM_TYPE_F + 
							inmigrant*INMIGRANT + male*MALE + young*AGE_YOUNG + midage*AGE_MID,
					a ~ dnorm(0, 100),
					comm ~ dnorm(0, sig_comm),
					sig_comm ~ dcauchy(0,1),
					dispersion ~ dexp(1),
					c(trading, fishing, inmigrant, male, young, midage) ~ dnorm(0,10)
			),
			data=as.data.frame(dg), 
			start=list(	a=0, comm=rep(0,39), dispersion=10, sig_comm=1, trading=0, fishing=0, inmigrant=0, male=0, young=0, midage=0),
			warmup=5e2, iter=5e3, chains=1, cores=4
		)		
	plot( precis(mp4, depth=2, prob=0.95) )
	#	posterior check to see if this makes sense
	sims 	<- sim(mp4, n=1e3)
	tmp		<- apply(sims, 2, median)
	dpp 	<- apply(sims, 2, PI, prob=0.95)
	dpp		<- rbind(tmp, dpp)
	rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
	dpp		<- as.data.table(t(dpp))
	dpp		<- cbind(dg, dpp)		
	dpp[, mean(PART<predicted_obs_l95 | PART>predicted_obs_u95)]
	#	0.01096491
	
	dpp[, OFF:= as.integer(PART<predicted_obs_l95 | PART>predicted_obs_u95)]
	set(dpp, NULL, 'AGE_YOUNG', dpp[,paste0('young_',AGE_YOUNG)])
	set(dpp, NULL, 'AGE_MID', dpp[,paste0('mid_',AGE_MID)])
	set(dpp, NULL, 'INMIGRANT', dpp[,paste0('inm_',INMIGRANT)])
	set(dpp, NULL, 'COMM_TYPE_F', dpp[,paste0('fish_',COMM_TYPE_F)])
	ggplot(dpp, aes(x=COMM_NUM_A)) +
			facet_grid(SEX+AGE_YOUNG+AGE_MID+INMIGRANT~COMM_TYPE_F, scales='free', space='free') +			
			geom_point(aes(y=predicted_obs_median), colour='grey50') +
			geom_errorbar(aes(ymin=predicted_obs_l95, ymax=predicted_obs_u95), colour='grey50') + 
			geom_point(aes(y=PART, colour=factor(OFF))) +
			theme_bw()
	ggsave(file=paste0(outfile.base,'sampling_differences_inmigrants_180322_mp4_pp.pdf'), w=16, h=16)
	
	compare(mp1, mp2, mp3, mp4, func=DIC)
	#	mp3 best

	tmp	<-  precis(mp3, depth=2, prob=0.95)@output
	do	<-	as.data.table(tmp)
	do[, VAR:=rownames(tmp)]
	setnames(do, colnames(do), toupper(gsub(' |\\.','_',colnames(do))))
	do[, COMM_NUM_B:= as.integer(gsub('.*\\[([0-9]+)\\].*','\\1',VAR))]
	tmp	<- unique(subset(dg, select=c(COMM_NUM_A, COMM_NUM_B)))
	do	<- merge(do, tmp, by='COMM_NUM_B', all.x=TRUE)
	tmp	<- do[, which(VAR=='dispersion')]
	set(do, tmp, 'MEAN', do[tmp,MEAN/20])
	set(do, tmp, 'LOWER_0_95', do[tmp,LOWER_0_95/20])
	set(do, tmp, 'UPPER_0_95', do[tmp,UPPER_0_95/20])
	tmp	<- do[, which(!is.na(COMM_NUM_A))]
	set(do, tmp, 'VAR', do[tmp,paste0('comm_',COMM_NUM_A)])	
	tmp	<- c("a", "dispersion", "male", "fishing", "trading", "inmigrant", "inmigrant_young", 
			"young_male", "young_female", "midage", "comm_aam", "comm_abj", "comm_abm", 
			"comm_abo", "comm_adi", "comm_aev", "comm_afr", "comm_agf", "comm_agv", "comm_agw", 
			"comm_ait", "comm_akh", "comm_ald", "comm_alr", "comm_amf", "comm_aop", "comm_aqi", 
			"comm_aqj", "comm_aro", "comm_asm", "comm_ati", "comm_aur", "comm_ave", "comm_avl", 
			"comm_avo", "comm_awa", "comm_awr", "comm_aws", "comm_axn", "comm_fno", "comm_fpt", 
			"comm_fus", "comm_fwd", "comm_thd", "comm_tpq", "comm_ttp", "comm_tvm", "comm_twr", "sig_comm")
	set(do, NULL, 'VAR', do[,factor(VAR, levels=tmp)])
	setkey(do, VAR)	
	ggplot(do, aes(x=VAR, y=MEAN)) +
			geom_abline(intercept=0, slope=0, lty=2, colour='grey50') +
			geom_point() +
			geom_errorbar(aes(ymin=LOWER_0_95, ymax=UPPER_0_95), width=NA) +
			theme_bw() +
			theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5)) +
			labs(y='estimate\n(posterior mean, 95% credibility interval)\n', x='')
	ggsave(file=paste0(outfile.base,'sampling_differences_inmigrants_180322_mp3_fit.pdf'), w=12, h=5, useDingbats=FALSE)
	
	save(dg, mp1, mp2, mp3, mp4, file=paste0(outfile.base,'participation_differences_180322_logisticmodels.rda'))
}


RakaiFull.phylogeography.180322.sequencing.bias.logisticmodel<- function()
{
	require(data.table)
	require(rethinking)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	opt.exclude.onART.from.denominator	<- 1
	
	#	load des which contains participation and seq counts by comm and gender
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender_age.rda'
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/sampling_differences_180322_"
	if(opt.exclude.onART.from.denominator)
	{
		outfile.base		<- paste0(outfile.base,'exclART_')
	}
	load(infile)	
	set(de, de[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)	
	de	<- subset(de, HIV_1517==1)
	if(opt.exclude.onART.from.denominator)
	{
		stopifnot(de[, !any(is.na(SELFREPORTART_AT_FIRST_VISIT))])
		de	<- subset(de, SELFREPORTART_AT_FIRST_VISIT!=1 | (SELFREPORTART_AT_FIRST_VISIT==1 & MIN_PNG_OUTPUT==1))		
	}
	
	
	tmp	<- de[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
	de	<- merge(de, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')
	
	#	logistic regression on participation
	#	set up variables for STAN
	dg		<- subset(de, select=c(MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID_C, INMIGRANT, SEX))
	#	binarize age, sex
	dg[, AGE_YOUNG:= as.integer(AGE_AT_MID_C=='15-24')]
	dg[, AGE_MID:= as.integer(AGE_AT_MID_C=='25-34')]
	dg[, MALE:= as.integer(SEX=='M')]
	#	binarize community type
	dg[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	dg[, COMM_TYPE_T:= as.integer(substr(COMM_NUM_A,1,1)=='t')]
	#	vanilla community IDs
	tmp		<- data.table(COMM_NUM_A= sort(unique(dg$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(dg$COMM_NUM_A)))
	dg		<- merge(dg, tmp, by='COMM_NUM_A')
	#	aggregate
	dg		<- dg[, list(HIV=length(PERM_ID), SEQ=length(which(MIN_PNG_OUTPUT==1))), by=c('AGE_AT_MID_C','INMIGRANT','SEX','AGE_YOUNG','AGE_MID','MALE','COMM_TYPE_F','COMM_TYPE_T','COMM_NUM_A','COMM_NUM_B')]
	
	#	run STAN
	ms1 	<- map2stan(
			alist(
					SEQ ~ dbinom(HIV, p_seq),
					logit(p_seq) <- a + comm[COMM_NUM_B] + trading*COMM_TYPE_T + fishing*COMM_TYPE_F + 
							inmigrant*INMIGRANT + male*MALE + young*AGE_YOUNG + midage*AGE_MID,
					a ~ dnorm(0, 100),
					comm ~ dnorm(0, sig_comm),
					sig_comm ~ dcauchy(0,1),
					c(trading, fishing, inmigrant, male, young, midage) ~ dnorm(0,10)
			),
			data=as.data.frame(dg), 
			start=list(	a=0, comm=rep(0,39), sig_comm=1, trading=0, fishing=0, inmigrant=0, male=0, young=0, midage=0),
			warmup=5e2, iter=5e3, chains=1, cores=4
		)		
	plot(ms1)
	plot( precis(ms1, depth=2, prob=0.95) )
	#	posterior check to see if this makes sense
	sims 	<- sim(ms1, n=1e3)
	tmp		<- apply(sims, 2, median)
	dpp 	<- apply(sims, 2, PI, prob=0.95)
	dpp		<- rbind(tmp, dpp)
	rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
	dpp		<- as.data.table(t(dpp))
	dpp		<- cbind(dg, dpp)		
	dpp[, c(mean(SEQ<predicted_obs_l95 | SEQ>predicted_obs_u95), sum(SEQ<predicted_obs_l95 | SEQ>predicted_obs_u95))]
	#	0.01028278, 4 -- that s already super as expected

	tmp		<- extract.samples(ms1)
p.tr.seq	<- with(tmp, a + comm[, TR_COMM_NUM_B2] + trading*TR_COMM_TYPE_T + fishing*TR_COMM_TYPE_F + 
				inmigrant*TR_INMIGRANT + male*TR_MALE + young*TR_AGE_YOUNG + midage*TR_AGE_MID)

	
	ggplot(dpp, aes(x=COMM_NUM_A)) +
			facet_grid(SEX+AGE_YOUNG+AGE_MID+INMIGRANT~COMM_TYPE_F, scales='free', space='free') +			
			geom_point(aes(y=predicted_obs_median), colour='grey50') +
			geom_errorbar(aes(ymin=predicted_obs_l95, ymax=predicted_obs_u95), colour='grey50') + 
			geom_point(aes(y=SEQ), colour='red') +
			theme_bw()
	ggsave(file=paste0(outfile.base,'ms1_pp.pdf'), w=16, h=32)
	#	
	
	
	ms4 	<- map2stan(
			alist(
					SEQ ~ dbetabinom(HIV, p_seq, dispersion),
					logit(p_seq) <- a + comm[COMM_NUM_B] + trading*COMM_TYPE_T + fishing*COMM_TYPE_F + 
							inmigrant*INMIGRANT + male*MALE + young*AGE_YOUNG + midage*AGE_MID,
					a ~ dnorm(0, 100),
					comm ~ dnorm(0, sig_comm),
					sig_comm ~ dcauchy(0,1),
					dispersion ~ dexp(1),
					c(trading, fishing, inmigrant, male, young, midage) ~ dnorm(0,10)
			),
			data=as.data.frame(dg), 
			start=list(	a=0, comm=rep(0,39), dispersion=10, sig_comm=1, trading=0, fishing=0, inmigrant=0, male=0, young=0, midage=0),
			warmup=5e2, iter=5e3, chains=1, cores=4
		)		
	plot( precis(ms4, depth=2, prob=0.95) )
	#	posterior check to see if this makes sense
	sims 	<- sim(ms4, n=1e3)
	tmp		<- apply(sims, 2, median)
	dpp 	<- apply(sims, 2, PI, prob=0.95)
	dpp		<- rbind(tmp, dpp)
	rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
	dpp		<- as.data.table(t(dpp))
	dpp		<- cbind(dg, dpp)		
	dpp[, c(mean(SEQ<predicted_obs_l95 | SEQ>predicted_obs_u95), sum(SEQ<predicted_obs_l95 | SEQ>predicted_obs_u95))]
	#	0 0	
	dpp[, OFF:= as.integer(SEQ<predicted_obs_l95 | SEQ>predicted_obs_u95)]
	set(dpp, NULL, 'AGE_YOUNG', dpp[,paste0('young_',AGE_YOUNG)])
	set(dpp, NULL, 'AGE_MID', dpp[,paste0('mid_',AGE_MID)])
	set(dpp, NULL, 'INMIGRANT', dpp[,paste0('inm_',INMIGRANT)])
	set(dpp, NULL, 'COMM_TYPE_F', dpp[,paste0('fish_',COMM_TYPE_F)])
	ggplot(dpp, aes(x=COMM_NUM_A)) +
			facet_grid(SEX+AGE_YOUNG+AGE_MID+INMIGRANT~COMM_TYPE_F, scales='free', space='free') +			
			geom_point(aes(y=predicted_obs_median), colour='grey50') +
			geom_errorbar(aes(ymin=predicted_obs_l95, ymax=predicted_obs_u95), colour='grey50') + 
			geom_point(aes(y=SEQ, colour=factor(OFF))) +
			theme_bw()
	ggsave(file=paste0(outfile.base,'ms4_pp.pdf'), w=16, h=16)
	



	ms2 	<- map2stan(
			alist(
					SEQ ~ dbinom(HIV, p_seq),
					logit(p_seq) <- a + comm[COMM_NUM_B] + trading*COMM_TYPE_T + fishing*COMM_TYPE_F + 
							inmigrant*INMIGRANT + inmigrant_young*INMIGRANT*AGE_YOUNG + 
							male*MALE + 
							young_male*AGE_YOUNG*MALE + young_female*AGE_YOUNG*(1-MALE) +
							midage*AGE_MID,
					a ~ dnorm(0, 100),
					comm ~ dnorm(0, sig_comm),
					sig_comm ~ dcauchy(0,1),
					c(trading, fishing, inmigrant, inmigrant_young, male, young_male, young_female, midage) ~ dnorm(0,10)
			),
			data=as.data.frame(dg), 
			start=list(	a=0, comm=rep(0,39), sig_comm=1, trading=0, fishing=0, inmigrant=0, inmigrant_young=0, male=0, young_male=0, young_female=0, midage=0),
			warmup=5e2, iter=5e3, chains=1, cores=4
		)
	plot( precis(ms2, depth=2, prob=0.95) )
	#	posterior check to see if this makes sense
	sims 	<- sim(ms2, n=1e3)
	tmp		<- apply(sims, 2, median)
	dpp 	<- apply(sims, 2, PI, prob=0.95)
	dpp		<- rbind(tmp, dpp)
	rownames(dpp)	<-  c('predicted_obs_median','predicted_obs_l95','predicted_obs_u95')
	dpp		<- as.data.table(t(dpp))
	dpp		<- cbind(dg, dpp)		
	dpp[, c(mean(SEQ<predicted_obs_l95 | SEQ>predicted_obs_u95), sum(SEQ<predicted_obs_l95 | SEQ>predicted_obs_u95))]
	#	0.007712082 3.000000000	
	dpp[, OFF:= as.integer(SEQ<predicted_obs_l95 | SEQ>predicted_obs_u95)]
	dpp[, list(OFF_P= mean(OFF), OFF_N=sum(OFF)), by=c('AGE_AT_MID_C','SEX','INMIGRANT')]	
	set(dpp, NULL, 'AGE_YOUNG', dpp[,paste0('young_',AGE_YOUNG)])
	set(dpp, NULL, 'AGE_MID', dpp[,paste0('mid_',AGE_MID)])
	set(dpp, NULL, 'INMIGRANT', dpp[,paste0('inm_',INMIGRANT)])
	set(dpp, NULL, 'COMM_TYPE_F', dpp[,paste0('fish_',COMM_TYPE_F)])
	ggplot(dpp, aes(x=COMM_NUM_A)) +
			facet_grid(SEX+AGE_YOUNG+AGE_MID+INMIGRANT~COMM_TYPE_F, scales='free', space='free') +			
			geom_point(aes(y=predicted_obs_median), colour='grey50') +
			geom_errorbar(aes(ymin=predicted_obs_l95, ymax=predicted_obs_u95), colour='grey50') + 
			geom_point(aes(y=SEQ, colour=factor(OFF))) +
			theme_bw()
	ggsave(file=paste0(outfile.base,'ms2_pp.pdf'), w=16, h=16)
	
	
	compare(ms1, ms2, ms4, func=DIC)
	#	ms1 best
	
	
	tmp	<-  precis(ms1, depth=2, prob=0.95)@output
	do	<-	as.data.table(tmp)
	do[, VAR:=rownames(tmp)]
	setnames(do, colnames(do), toupper(gsub(' |\\.','_',colnames(do))))
	do[, COMM_NUM_B:= as.integer(gsub('.*\\[([0-9]+)\\].*','\\1',VAR))]
	tmp	<- unique(subset(dg, select=c(COMM_NUM_A, COMM_NUM_B)))
	do	<- merge(do, tmp, by='COMM_NUM_B', all.x=TRUE)
	tmp	<- do[, which(!is.na(COMM_NUM_A))]
	set(do, tmp, 'VAR', do[tmp,paste0('comm_',COMM_NUM_A)])	
	tmp	<- c("a", "male", "fishing", "trading", "inmigrant",  
			"young", "midage", "comm_aam", "comm_abj", "comm_abm", 
			"comm_abo", "comm_adi", "comm_aev", "comm_afr", "comm_agf", "comm_agv", "comm_agw", 
			"comm_ait", "comm_akh", "comm_ald", "comm_alr", "comm_amf", "comm_aop", "comm_aqi", 
			"comm_aqj", "comm_aro", "comm_asm", "comm_ati", "comm_aur", "comm_ave", "comm_avl", 
			"comm_avo", "comm_awa", "comm_awr", "comm_aws", "comm_axn", "comm_fno", "comm_fpt", 
			"comm_fus", "comm_fwd", "comm_thd", "comm_tpq", "comm_ttp", "comm_tvm", "comm_twr", "sig_comm")
	set(do, NULL, 'VAR', do[,factor(VAR, levels=tmp)])
	do	<- subset(do, !is.na(VAR))
	setkey(do, VAR)	
	ggplot(do, aes(x=VAR, y=MEAN)) +
			geom_abline(intercept=0, slope=0, lty=2, colour='grey50') +
			geom_point() +
			geom_errorbar(aes(ymin=LOWER_0_95, ymax=UPPER_0_95), width=NA) +
			theme_bw() +
			theme(axis.text.x=element_text(angle = 90, hjust=1, vjust=0.5)) +
			labs(y='estimate\n(posterior mean, 95% credibility interval)\n', x='')
	ggsave(file=paste0(outfile.base,'ms1_fit.pdf'), w=12, h=5, useDingbats=FALSE)
	
	save(dg, ms1, ms2, ms4, file=paste0(outfile.base,'logisticmodels.rda'))
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