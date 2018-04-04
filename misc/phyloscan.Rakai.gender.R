RakaiFull.gender.171122.propfemalepos.communities.merged<- function()
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
	set(df, NULL, 'COMM_NUM', df[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',as.character(COMM_NUM)))))])
	df	<- subset(df, VISIT%in%c(15,15.1,16))
	
	#	select relevant communities
	tmp	<- sort(unique(c(as.character(rtpdm$MALE_COMM_NUM), as.character(rtpdm$MALE_COMM_NUM))))
	tmp	<- data.table(COMM_NUM=tmp)
	df	<- merge(df, tmp, by=c('COMM_NUM'))
	df[, FEMALE:= FEMALE_NEG+FEMALE_POS]
	df[, MALE:= MALE_NEG+MALE_POS]
	df[, POS:= MALE_POS+FEMALE_POS]
	df[, NEG:= MALE_NEG+FEMALE_NEG]
	df[, FEMALE_POS_P:= FEMALE_POS/POS]
	df	<- subset(df, MALE!=0 | FEMALE!=0)
	
	#	divide into groups by proportion of seropos who are female
	tmp	<- df[, list(FEMALE_POS_PM=mean(FEMALE_POS_P)), by='COMM_NUM']
	setkey(tmp, FEMALE_POS_PM)
	tmp[, FEMALE_POS_C:= cut(FEMALE_POS_PM, breaks=c(0.5, 0.55, 0.56, 0.58, 0.68, 0.69,  1), labels=c('a','b','c','d','e','f'))]
	
	ggplot(tmp, aes(x=FEMALE_POS_PM, fill=FEMALE_POS_C)) + geom_histogram()
	df	<- merge(df, subset(tmp, select=c(COMM_NUM, FEMALE_POS_C)), by='COMM_NUM')
	
	#	sum POS and FEMALE_POS by community groups	
	dg	<- df[, list(FEMALE_POS=sum(FEMALE_POS), POS=sum(POS)), by='FEMALE_POS_C']
	
	#	add community groups to rtpdm
	tmp	<- unique(subset(df, select=c(COMM_NUM,FEMALE_POS_C)))
	setnames(tmp, 'COMM_NUM', 'MALE_COMM_NUM')
	tmp	<- merge(subset(rtpdm, MALE_COMM_TYPE==FEMALE_COMM_TYPE), tmp, by='MALE_COMM_NUM')
	tmp	<- tmp[, {
				z<- as.numeric(binconf(length(which(grepl('mf',SELECT))), length(SELECT)))
				list(	MF_TRM	= length(which(grepl('mf',SELECT))),
						MF_TRM_RATIO	= length(which(grepl('mf',SELECT)))/length(which(grepl('fm',SELECT))),
						LMF_TRM_RATIO	= log(length(which(grepl('mf',SELECT)))/length(which(grepl('fm',SELECT)))),
						TRM		= length(SELECT),
						MF_MED	= z[1],
						MF_CL	= z[2],
						MF_CU	= z[3])
			}, by='FEMALE_POS_C']
	dg	<- merge(dg, tmp, by='FEMALE_POS_C')
	tmp	<- dg[, {
				z<- binconf(FEMALE_POS,POS)
				list(FEMALE_POS_M=z[1], FEMALE_POS_CL=z[2], FEMALE_POS_CU=z[3])
			}, by='FEMALE_POS_C']
	dg	<- merge(dg, tmp, by='FEMALE_POS_C')	
	

	mhi.4 	<- map2stan(
			alist(
					MF_TRM ~ dbinom(TRM, ptm),
					logit(ptm) <- a + b*logit(pdiagf[i]),
					FEMALE_POS ~ dbinom(POS, pdiagf),									 
					a ~ dnorm(0,100),
					b ~ dnorm(0,10),
					pdiagf ~ dnorm(0.5,1)
			),
			data=as.data.frame(dg), start=list(a=0, b=0, pdiagf=rep(0.5,6)),
			warmup=5e2, iter=5e3, chains=1, cores=4
	)	
	summary(mhi.4)	# b posterior median nearly 1; credibility intervals wide but this is good
	post	<- extract.samples(mhi.4)
	dummy	<- function(z) quantile(logistic(with(post, a+b*z)), prob=c(0.5, 0.025, 0.975))
	tmp		<- data.table(FEMALE_POS_P=seq(0.5,0.75,0.001))
	tmp		<- tmp[, 	{
				z<- dummy(logit(FEMALE_POS_P))
				list('MF_median'=z[1], 'MF_cl'=z[2], 'MF_cu'=z[3])
			}, by='FEMALE_POS_P']	
	ggplot(dg) +
			geom_ribbon(data=tmp, aes(x=FEMALE_POS_P, ymin=MF_cl, ymax=MF_cu), fill='black', alpha=0.25) +
			geom_line(data=tmp, aes(x=FEMALE_POS_P, y=MF_median), colour='grey50', size=1) +
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_point(aes(x=FEMALE_POS_M, y=MF_MED, colour=FEMALE_POS_C), size=2) +
			geom_errorbar(aes(x=FEMALE_POS_M, ymin=MF_CL, ymax=MF_CU, colour=FEMALE_POS_C), width=0.03, size=0.9) +
			geom_errorbarh(aes(x=FEMALE_POS_M, y=MF_MED, xmin=FEMALE_POS_CL, xmax=FEMALE_POS_CU, colour=FEMALE_POS_C), height=0.05, size=0.9) +
			scale_x_continuous(labels=scales:::percent, breaks=seq(0,1,0.05)) +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1)) +
			scale_colour_brewer(palette='Set2') +
			theme_bw() +
			labs(	x='\nproportion of females among newly diagnosed', 
					y='male to female transmission\namong phylogenetically inferred transmission events\n',
					colour='community group')
	ggsave(file=paste0(outfile.base,'_trmMF_vs_diagFM_by_commgroup.pdf'), w=6, h=4.5)
}

RakaiFull.gender.171122.propfemalepos.communitytype<- function()
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
	set(df, NULL, 'COMM_NUM', df[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',as.character(COMM_NUM)))))])
	df	<- subset(df, VISIT%in%c(15,15.1,16))
	
	tmp	<- sort(unique(c(as.character(rtpdm$MALE_COMM_NUM), as.character(rtpdm$MALE_COMM_NUM))))
	tmp	<- data.table(COMM_NUM=tmp)
	df	<- merge(df, tmp, by=c('COMM_NUM'))
	df	<- df[, list(FEMALE_NEG=sum(FEMALE_NEG), MALE_NEG=sum(MALE_NEG), FEMALE_POS=sum(FEMALE_POS), MALE_POS=sum(MALE_POS)), by='COMM_NUM']
	
	tmp	<- unique(subset(rtpdm, select=c(MALE_COMM_NUM,MALE_COMM_TYPE)))
	setnames(tmp, c('MALE_COMM_NUM','MALE_COMM_TYPE'), c('COMM_NUM','COMM_TYPE'))
	df	<- merge(df, tmp, by='COMM_NUM')
	dg	<- df[, list(FEMALE_NEG=sum(FEMALE_NEG), MALE_NEG=sum(MALE_NEG), FEMALE_POS=sum(FEMALE_POS), MALE_POS=sum(MALE_POS)), by='COMM_TYPE']
	
	dg[, MALE_PROP_DIAG:= MALE_POS/(MALE_NEG+MALE_POS)]
	dg[, FEMALE_PROP_DIAG:= FEMALE_POS/(FEMALE_NEG+FEMALE_POS)]
	dg[, FM_DIAG_RATIO:=FEMALE_PROP_DIAG/MALE_PROP_DIAG]
	dg[, LFM_DIAG_RATIO:=log(FM_DIAG_RATIO)]
	dg[, FEMALE:= FEMALE_NEG+FEMALE_POS]
	dg[, MALE:= MALE_NEG+MALE_POS]
	dg[, POS:= MALE_POS+FEMALE_POS]
	dg[, NEG:= MALE_NEG+FEMALE_NEG]
	tmp	<- dg[, {
				z<- binconf(FEMALE_POS,POS)
				list(FM_DIAG_PROP_M=z[1], FM_DIAG_PROP_CL=z[2], FM_DIAG_PROP_CU=z[3])
			}, by='COMM_TYPE']
	dg	<- merge(dg, tmp, by='COMM_TYPE')	
	dg[, binconf(sum(FEMALE_POS), sum(POS))]				#	0.6008119 0.5901207 0.6114079
	dg[, exp(logit(binconf(sum(FEMALE_POS), sum(POS))))]	#	1.505085  1.439742  1.573393
	
	tmp	<- subset(rtpdm, MALE_COMM_TYPE==FEMALE_COMM_TYPE)
	tmp	<- tmp[, {
				z<- as.numeric(binconf(length(which(grepl('mf',SELECT))), length(SELECT)))
				list(	MF_TRM	= length(which(grepl('mf',SELECT))),
						MF_TRM_RATIO	= length(which(grepl('mf',SELECT)))/length(which(grepl('fm',SELECT))),
						LMF_TRM_RATIO	= log(length(which(grepl('mf',SELECT)))/length(which(grepl('fm',SELECT)))),
						TRM		= length(SELECT),
						MF_MED	= z[1],
						MF_CL	= z[2],
						MF_CU	= z[3])
			}, by='MALE_COMM_TYPE']
	setnames(tmp, 'MALE_COMM_TYPE', 'COMM_TYPE')
	dg	<- merge(dg, tmp, by='COMM_TYPE')
	dg[, COMM_TYPE_2:= as.integer(factor(COMM_TYPE))]
	
	ggplot(dg, aes(x=FM_DIAG_RATIO, y=MF_MED, ymin=MF_CL, ymax=MF_CU)) +
			geom_point() +
			geom_errorbar() +
			theme_bw()
	
	mh.2 	<- map2stan(
			alist(
					MF_TRM ~ dbinom(TRM, ptm),
					logit(ptm) <- a + b*pdiagf[i],
					FEMALE_POS ~ dbinom(POS, pdiagf),									 
					a ~ dnorm(0,100),
					b ~ dnorm(0,10),
					#pdiagf ~ dbeta(1,1),
					pdiagf ~ dnorm(0.5,1)
			),
			data=as.data.frame(dg), start=list(a=0, b=0, pdiagf=rep(0.5,3)),
			warmup=5e2, iter=5e3, chains=1, cores=4
	)	
	mh.3 	<- map2stan(
			alist(
					MF_TRM ~ dbinom(TRM, ptm),
					logit(ptm) <- a + b*logit(pdiagf[i]),
					FEMALE_POS ~ dbinom(POS, pdiagf),									 
					a ~ dnorm(0,100),
					b ~ dnorm(0,10),
					pdiagf ~ dnorm(0.5,1)
			),
			data=as.data.frame(dg), start=list(a=0, b=0, pdiagf=rep(0.5,3)),
			warmup=5e2, iter=5e3, chains=1, cores=4
	)
	
	summary(mh.3)	# b posterior median nearly 1; credibility intervals wide but this is good
	post	<- extract.samples(mh.3)
	dummy	<- function(z) quantile(logistic(with(post, a+b*z)), prob=c(0.5, 0.025, 0.975))
	tmp		<- data.table(FM_DIAG_PROP=seq(0.5,0.75,0.001))
	tmp		<- tmp[, {
				z<- dummy(logit(FM_DIAG_PROP))
				list('MF_median'=z[1], 'MF_cl'=z[2], 'MF_cu'=z[3])
			}, by='FM_DIAG_PROP']
	
	ggplot(dg) +
			geom_ribbon(data=tmp, aes(x=FM_DIAG_PROP, ymin=MF_cl, ymax=MF_cu), fill='black', alpha=0.25) +
			geom_line(data=tmp, aes(x=FM_DIAG_PROP, y=MF_median), colour='grey50', size=1) +
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_point(aes(x=FM_DIAG_PROP_M, y=MF_MED, colour=COMM_TYPE), size=2) +
			geom_errorbar(aes(x=FM_DIAG_PROP_M, ymin=MF_CL, ymax=MF_CU, colour=COMM_TYPE), width=0.03, size=0.9) +
			geom_errorbarh(aes(x=FM_DIAG_PROP_M, y=MF_MED, xmin=FM_DIAG_PROP_CL, xmax=FM_DIAG_PROP_CU, colour=COMM_TYPE), height=0.05, size=0.9) +
			scale_x_continuous(labels=scales:::percent, breaks=seq(0,1,0.05)) +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1)) +
			scale_colour_brewer(palette='Set2') +
			theme_bw() +
			labs(	x='\nproportion of females among newly diagnosed', 
					y='male to female transmission\namong phylogenetically inferred transmission events\n',
					colour='community type')
	ggsave(file=paste0(outfile.base,'_trmMF_vs_diagFM_by_commtype.pdf'), w=6, h=4.5)
}

RakaiFull.gender.171122.hh.are.femalerecipients.transmitters<- function()
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
	
	#zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	#zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
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
	
	dfr	<- subset(rtpdm, grepl('mf',SELECT), select=c(FEMALE_RID, SAMEHH))
	setnames(dfr, 'SAMEHH', 'SAMEHH_OF_FEMALE_RECIPIENT')
	tmp	<- subset(rtpdm, grepl('fm',SELECT))
	dfr	<- merge(dfr, tmp, by='FEMALE_RID', all.x=1)
	tmp	<- dfr[, list(N_TRMS= length(which(!is.na(MALE_RID))) ), by=c('FEMALE_RID','SAMEHH_OF_FEMALE_RECIPIENT')]
	tmp[, table(SAMEHH_OF_FEMALE_RECIPIENT, N_TRMS>=1)]
	
	fisher.test( matrix(c(6,1,83,83),2,2) )
	#	p-value = 0.1185
	#	alternative hypothesis: true odds ratio is not equal to 1
	#95 percent confidence interval:
   	#	0.6984958 278.8704511
	#	sample estimates:
	#odds ratio 
  	#	5.949555 
}

RakaiFull.gender.171122.hh.sources.of.transmitters<- function()
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
	
	#zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	#zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
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
	
	dft	<- subset(rtpdm, grepl('mf',SELECT), select=c(MALE_RID, SAMEHH))
	setnames(dft, 'SAMEHH', 'SAMEHH_OF_MALE_TRANSMITTER')
	tmp	<- subset(rtpdm, grepl('fm',SELECT))
	dft	<- merge(dft, tmp, by='MALE_RID', all.x=1)
	tmp	<- dft[, list(N_TRMS= length(which(!is.na(FEMALE_RID))) ), by=c('MALE_RID','SAMEHH_OF_MALE_TRANSMITTER')]
	tmp[, table(SAMEHH_OF_MALE_TRANSMITTER, N_TRMS>=1)]
	
	fisher.test( matrix(c(6,1,83,83),2,2) )
	#	p-value = 0.1185
	#	alternative hypothesis: true odds ratio is not equal to 1
	#95 percent confidence interval:
	#	0.6984958 278.8704511
	#	sample estimates:
	#odds ratio 
	#	5.949555 
}

RakaiFull.gender.171122.hh.multivariatemodelswithprevalence.stan.with.threshold<- function()
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
	
	#zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	#zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	load(infile)	
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	rtpdm[, PAIR_COMM_TYPE:= FEMALE_COMM_TYPE]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_TYPE!=MALE_COMM_TYPE)], 'PAIR_COMM_TYPE', 'mixed')
	set(rtpdm, rtpdm[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(rtpdm, rtpdm[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')
	
	#rtpdm[, table(PAIR_COMM_TYPE)]
	rtpdm[, COUPLE2:= factor(COUPLE=='no couple', levels=c(TRUE,FALSE), labels=c('no couple','couple'))]
	rtpdm[, SAMEHH:= factor(FEMALE_HH_NUM==MALE_HH_NUM, levels=c(TRUE,FALSE), labels=c('same hh','different hh'))]	
	rtpdm[, PAIR_COMM:= MALE_COMM_NUM_A]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_NUM_A!=MALE_COMM_NUM_A)], 'PAIR_COMM', 'mixed')
	rtpdm[, MALE_SEXP1OUT2:= factor(MALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(MALE_SEXP1OUT=='Unknown')], 'MALE_SEXP1OUT2', 'Unknown')
	rtpdm[, FEMALE_SEXP1OUT2:= factor(FEMALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(FEMALE_SEXP1OUT=='Unknown')], 'FEMALE_SEXP1OUT2', 'Unknown')
	
	
	df		<- subset(rtpdm, FEMALE_HH_NUM==MALE_HH_NUM, select=c(	MALE_RID, FEMALE_RID, PTY_RUN, IDCLU, LINK_MF, POSTERIOR_SCORE_MF, COUPLE2, SAMEHH, PAIR_COMM_TYPE, FEMALE_COMM_NUM,
					MALE_RECENTVL, MALE_SEXP1YR, MALE_SEXP1OUT2, MALE_OCCUP_OLLI, MALE_OCAT, MALE_EDUCAT, MALE_CIRCUM,
					FEMALE_RECENTVL, FEMALE_SEXP1YR, FEMALE_SEXP1OUT2, FEMALE_OCCUP_OLLI, FEMALE_OCAT, FEMALE_EDUCAT 
			))
	#	missing data: fill in
	set(df, df[, which(is.na(MALE_RECENTVL))], 'MALE_RECENTVL', -1)
	set(df, df[, which(is.na(FEMALE_RECENTVL))], 'FEMALE_RECENTVL', -1)
	set(df, df[, which(is.na(MALE_CIRCUM))], 'MALE_CIRCUM', 'Unknown')
	set(df, df[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(df, df[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')
	#	add estimated HIV prevalences (medians)
	infile	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30trmMF_estimated_prevalenceratio.rda"	
	load(infile)
	tmp		<- dcast.data.table(dgg, COMM_NUM~variable, value.var='M')
	setnames(tmp, 'COMM_NUM', 'FEMALE_COMM_NUM')
	df		<- merge(df, tmp, by='FEMALE_COMM_NUM')
	
	#for(x in colnames(df)) print( c(x, any(is.na(df[[x]]))) )
	#	prepare data for STAN
	df[, COMM_FISH:= as.integer(PAIR_COMM_TYPE=='fisherfolk')]
	df[, COMM_AGR:= as.integer(PAIR_COMM_TYPE=='agrarian')]
	df[, COMM_TRAD:= as.integer(PAIR_COMM_TYPE=='trading')]
	df[, COMM_MXD:= as.integer(PAIR_COMM_TYPE=='mixed')]
	df[, FE_NOEDU:= as.integer(FEMALE_EDUCAT=='None')]
	df[, FE_NOEDU_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, MA_NOEDU:= as.integer(MALE_EDUCAT=='None')]
	df[, MA_NOEDU_MISS:= as.integer(MALE_EDUCAT=='Unknown')]	
	df[, MA_CIRCUM:= as.integer(MALE_CIRCUM=='Y')]
	df[, MA_CIRCUM_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, FE_SEXP1YR_G1:= as.integer(FEMALE_SEXP1YR!='1')]
	df[, MA_SEXP1YR_G1:= as.integer(MALE_SEXP1YR!='1')]
	df[, MA_OC_AGRO:= as.integer(MALE_OCAT=='Agro/House')]	
	df[, MA_OC_FISH:= as.integer(MALE_OCAT=='Fishing')]
	df[, MA_OC_TRAD:= as.integer(MALE_OCAT=='Trading/Shop keeper')]
	df[, MA_OC_OTH:= as.integer(MALE_OCAT%in%c('zother','Bar/waitress','Student','Boda/Trucking'))]	
	df[, FE_OC_AGRO:= as.integer(FEMALE_OCAT=='Agro/House')]
	df[, FE_OC_BAR:= as.integer(FEMALE_OCAT=='Bar/waitress')]			
	df[, FE_OC_TRAD:= as.integer(FEMALE_OCAT=='Trading/Shop keeper')]
	df[, FE_OC_OTH:= as.integer(FEMALE_OCAT%in%c('zother','Student'))]		
	df[, PF2:= PF-mean(PF)]
	df[, PM2:= PM-mean(PM)]
	df[, RFM2:= RFM-mean(RFM)]
	#
	#	STAN
	#
	
	#	same household definitely stronger effect than couples	
	ms.1 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + 
							# comm_fish*COMM_FISH is base 
							comm_agr*COMM_AGR + comm_trad*COMM_TRAD + comm_mxd*COMM_MXD +
							noedu_f*FE_NOEDU + noedu_m*MA_NOEDU +
							sexp_f*FE_SEXP1YR_G1 + sexp_m*MA_SEXP1YR_G1 +
							circum_m*MA_CIRCUM +								
							prev_m*PM2 + prev_f*PF2 + prev_ratio_fm*RFM2 +
							# agro occupation is baseline
							mocc_other*MA_OC_OTH + mocc_fish*MA_OC_FISH + mocc_trad*MA_OC_TRAD +										
							focc_other*FE_OC_OTH + focc_bar*FE_OC_BAR +  focc_trad*FE_OC_TRAD,
					base ~ dnorm(0,100),								
					c(comm_agr, comm_trad, comm_mxd) ~ dnorm(0,10),
					c(noedu_f, noedu_m) ~ dnorm(0,10),								
					c(sexp_f, sexp_m, circum_m) ~ dnorm(0,10),
					c(prev_m, prev_f, prev_ratio_fm) ~ dnorm(0,10),
					c(mocc_other, mocc_fish, mocc_trad) ~ dnorm(0,10),
					c(focc_other, focc_bar, focc_trad) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, comm_agr=0, comm_trad=0, comm_mxd=0, noedu_f=0, noedu_m=0, sexp_f=0, sexp_m=0, circum_m=0,
					prev_m=0, prev_f=0, prev_ratio_fm=0,
					focc_other=0, focc_bar=0, focc_trad=0,
					mocc_other=0, mocc_fish=0, mocc_trad=0),			
			warmup=1e3, iter=5e3, chains=1, cores=4
	)		
	#plot(precis(ms.1, prob=0.95))
	#pairs(ms.1)
	post	<- extract.samples(ms.1)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_COMM_FISH=logistic( post$base ),				
			PI_COMM_AGRO=logistic( post$base+post$comm_agr ), 											
			PI_COMM_TRAD=logistic( post$base+post$comm_trad ),							
			PI_COMM_MXD=logistic( post$base+post$comm_mxd ),
			PI_FEDU_NONE=logistic( post$base+post$noedu_f ),
			PI_MEDU_NONE=logistic( post$base+post$noedu_m ),							
			PI_FESEXP1YR_G1=logistic( post$base+post$sexp_f ),
			PI_MASEXP1YR_G1=logistic( post$base+post$sexp_m ),							
			PI_CIRC_YES=logistic( post$base+post$circum_m ),							
			PI_MOCC_OTH=logistic( post$base+post$mocc_other ),
			PI_MOCC_FISH=logistic( post$base+post$mocc_fish ),
			PI_MOCC_TRAD=logistic( post$base+post$mocc_trad ),
			PI_FOCC_OTH=logistic( post$base+post$focc_other ),
			PI_FOCC_BAR=logistic( post$base+post$focc_bar ),
			PI_FOCC_TRAD=logistic( post$base+post$focc_trad ),			
			PI_PREV_MALE=logistic( post$base+post$prev_m ),
			PI_PREV_FEMALE=logistic( post$base+post$prev_f ),
			PI_PREV_RATIOFM=logistic( post$base+post$prev_ratio_fm ),			
			OR_COMM_FISH=exp( post$base ),				
			OR_COMM_AGRO=exp( post$base+post$comm_agr ), 											
			OR_COMM_TRAD=exp( post$base+post$comm_trad ),							
			OR_COMM_MXD=exp( post$base+post$comm_mxd ),
			OR_FEDU_NONE=exp( post$base+post$noedu_f ),
			OR_MEDU_NONE=exp( post$base+post$noedu_m ),							
			OR_FESEXP1YR_G1=exp( post$base+post$sexp_f ),
			OR_MASEXP1YR_G1=exp( post$base+post$sexp_m ),
			OR_CIRC_YES=exp( post$base+post$circum_m ),							
			OR_MOCC_OTH=exp( post$base+post$mocc_other ),
			OR_MOCC_FISH=exp( post$base+post$mocc_fish ),
			OR_MOCC_TRAD=exp( post$base+post$mocc_trad ),
			OR_FOCC_OTH=exp( post$base+post$focc_other ),
			OR_FOCC_BAR=exp( post$base+post$focc_bar ),
			OR_FOCC_TRAD=exp( post$base+post$focc_trad ),			
			OR_PREV_MALE=exp( post$base+post$prev_m ),
			OR_PREV_FEMALE=exp( post$base+post$prev_f ),
			OR_PREV_RATIOFM=exp( post$base+post$prev_ratio_fm ),
			ORX_COMM_AGRO=exp( post$comm_agr ), 											
			ORX_COMM_TRAD=exp( post$comm_trad ),							
			ORX_COMM_MXD=exp( post$comm_mxd ),
			ORX_FEDU_NONE=exp( post$noedu_f ),
			ORX_MEDU_NONE=exp( post$noedu_m ),
			ORX_FESEXP1YR_G1=exp( post$sexp_f ),
			ORX_MASEXP1YR_G1=exp( post$sexp_m ),							
			ORX_CIRC_YES=exp( post$circum_m ),							
			ORX_MOCC_OTH=exp( post$mocc_other ),
			ORX_MOCC_FISH=exp( post$mocc_fish ),
			ORX_MOCC_TRAD=exp( post$mocc_trad ),
			ORX_FOCC_OTH=exp( post$focc_other ),
			ORX_FOCC_BAR=exp( post$focc_bar ),
			ORX_FOCC_TRAD=exp( post$focc_trad ),
			ORX_PREV_MALE=exp( post$prev_m ),
			ORX_PREV_FEMALE=exp( post$prev_f ),
			ORX_PREV_RATIOFM=exp( post$prev_ratio_fm )
	)
	dp		<- melt(dp, id.vars='MC')	
	dss.1	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	dss.1	<- dcast.data.table(dss.1, variable~STAT, value.var='V')
	dss.1[, MODEL:= 'multivariate ms3']
	ds		<- copy(dss.1)
	ds[, LABEL:= paste0(round(MED, d=2), '\n[', round(CL, d=2),'-', round( CU, d=2),']')]
	ds[, STAT:= factor(gsub('^([^_]+)_.*','\\1',variable), levels=c('PI','OR','ORX'), labels=c('proportion MF','odds MF','odds ratio'))]
	ds[, FACTOR:=gsub('^([^_]+)_([^_]+)_([^_]+)','\\3',variable)]
	ds[, MOFA:=gsub('^([^_]+)_([^_]+)_([^_]+)','\\2-\\3',variable)]
	ds[, MODELTYPE:= factor(grepl('univariate',MODEL), levels=c(TRUE,FALSE),labels=c('univariate','multivariate'))]	
	ds		<- subset(ds, STAT=='odds ratio')
	set(ds, NULL, 'MOFA2', ds[, factor(MOFA, levels=rev(c(
											"COMM-AGRO","COMM-TRAD","COMM-MXD",																						
											"FEDU-NONE", "MEDU-NONE",      
											"FESEXP1YR-G1","MASEXP1YR-G1", 
											"CIRC-YES",                   
											"FOCC-BAR","FOCC-TRAD","FOCC-OTH",
											"MOCC-FISH","MOCC-OTH","MOCC-TRAD"   
									)), labels=rev(c("Agrarian community vs. fishing community", "Trading community vs. fishing community", "Mixed vs. fishing community",
											"Female education: None vs. at least primary education", "Male education: None vs. at least primary education",
											"Female sex partners in last year: >1 vs. 1", "Male sex partners in last year: >1 vs. 1",
											"Male circumcised: yes vs. no",
											"Female occupation: bar/waitress vs. agricultural/house", "Female occupation: trading/shopkeeper vs. agricultural/house","Female occupation: other vs. agricultural/house", 
											"Male occupation: fishing vs. agricultural/house", "Male occupation: other vs. agricultural/house", "Male occupation: trading/shopkeeper vs. agricultural/house"
									)))])
	setkey(ds, MOFA)	
	dp		<- subset(ds, !MOFA%in%c('COMM-TRAD','COMM-MXD'))
	dp[, FILL:= '0']
	set(dp, dp[, which(IL>1 | IU<1)], 'FILL', '1')
	set(dp, dp[, which(CL>1 | CU<1)], 'FILL', '2')
	ggplot(dp, aes(x=MOFA2)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/10, 1/4,1/3,1/2,1,2,3,4,10), labels=c('1/10','1/4','1/3','1/2','1','2','3','4','10')) +
			scale_fill_manual(values=c('0'='white', '1'='orange', '2'='red')) +
			coord_flip() +
			guides(fill=FALSE) +
			labs(x='Partner characteristics\n', y='\nOdds ratio in male to female transmission')
	ggsave(file=paste0(outfile.base,'_propmf_factors_univariate_samehh.pdf'), w=12, h=7)		
}

RakaiFull.gender.171122.hh.multivariatemodels.stan.with.threshold<- function()
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
	
	#zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	#zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	load(infile)	
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	rtpdm[, PAIR_COMM_TYPE:= FEMALE_COMM_TYPE]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_TYPE!=MALE_COMM_TYPE)], 'PAIR_COMM_TYPE', 'mixed')
	set(rtpdm, rtpdm[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(rtpdm, rtpdm[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')
	
	#rtpdm[, table(PAIR_COMM_TYPE)]
	rtpdm[, COUPLE2:= factor(COUPLE=='no couple', levels=c(TRUE,FALSE), labels=c('no couple','couple'))]
	rtpdm[, SAMEHH:= factor(FEMALE_HH_NUM==MALE_HH_NUM, levels=c(TRUE,FALSE), labels=c('same hh','different hh'))]	
	rtpdm[, PAIR_COMM:= MALE_COMM_NUM_A]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_NUM_A!=MALE_COMM_NUM_A)], 'PAIR_COMM', 'mixed')
	rtpdm[, MALE_SEXP1OUT2:= factor(MALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(MALE_SEXP1OUT=='Unknown')], 'MALE_SEXP1OUT2', 'Unknown')
	rtpdm[, FEMALE_SEXP1OUT2:= factor(FEMALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(FEMALE_SEXP1OUT=='Unknown')], 'FEMALE_SEXP1OUT2', 'Unknown')
	
	
	df		<- subset(rtpdm, FEMALE_HH_NUM==MALE_HH_NUM, select=c(	MALE_RID, FEMALE_RID, PTY_RUN, IDCLU, LINK_MF, POSTERIOR_SCORE_MF, COUPLE2, SAMEHH, PAIR_COMM_TYPE, FEMALE_COMM_NUM,
					MALE_RECENTVL, MALE_SEXP1YR, MALE_SEXP1OUT2, MALE_OCCUP_OLLI, MALE_OCAT, MALE_EDUCAT, MALE_CIRCUM,
					FEMALE_RECENTVL, FEMALE_SEXP1YR, FEMALE_SEXP1OUT2, FEMALE_OCCUP_OLLI, FEMALE_OCAT, FEMALE_EDUCAT 
			))
	#	missing data: fill in
	set(df, df[, which(is.na(MALE_RECENTVL))], 'MALE_RECENTVL', -1)
	set(df, df[, which(is.na(FEMALE_RECENTVL))], 'FEMALE_RECENTVL', -1)
	set(df, df[, which(is.na(MALE_CIRCUM))], 'MALE_CIRCUM', 'Unknown')
	set(df, df[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(df, df[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')	
	#for(x in colnames(df)) print( c(x, any(is.na(df[[x]]))) )
	#	add estimated HIV prevalences (medians)
	infile	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30trmMF_estimated_prevalenceratio.rda"	
	load(infile)
	tmp		<- dcast.data.table(dgg, COMM_NUM~variable, value.var='M')
	setnames(tmp, 'COMM_NUM', 'FEMALE_COMM_NUM')
	df		<- merge(df, tmp, by='FEMALE_COMM_NUM')

	#	prepare data for STAN
	df[, COMM_FISH:= as.integer(PAIR_COMM_TYPE=='fisherfolk')]
	df[, COMM_AGR:= as.integer(PAIR_COMM_TYPE=='agrarian')]
	df[, COMM_TRAD:= as.integer(PAIR_COMM_TYPE=='trading')]
	df[, COMM_MXD:= as.integer(PAIR_COMM_TYPE=='mixed')]
	df[, FE_NOEDU:= as.integer(FEMALE_EDUCAT=='None')]
	df[, FE_NOEDU_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, MA_NOEDU:= as.integer(MALE_EDUCAT=='None')]
	df[, MA_NOEDU_MISS:= as.integer(MALE_EDUCAT=='Unknown')]	
	df[, MA_CIRCUM:= as.integer(MALE_CIRCUM=='Y')]
	df[, MA_CIRCUM_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, FE_SEXP1YR_G1:= as.integer(FEMALE_SEXP1YR!='1')]
	df[, MA_SEXP1YR_G1:= as.integer(MALE_SEXP1YR!='1')]
	df[, MA_OC_AGRO:= as.integer(MALE_OCAT=='Agro/House')]	
	df[, MA_OC_FISH:= as.integer(MALE_OCAT=='Fishing')]
	df[, MA_OC_TRAD:= as.integer(MALE_OCAT=='Trading/Shop keeper')]
	df[, MA_OC_OTH:= as.integer(MALE_OCAT%in%c('zother','Bar/waitress','Student','Boda/Trucking'))]	
	df[, FE_OC_AGRO:= as.integer(FEMALE_OCAT=='Agro/House')]
	df[, FE_OC_BAR:= as.integer(FEMALE_OCAT=='Bar/waitress')]			
	df[, FE_OC_TRAD:= as.integer(FEMALE_OCAT=='Trading/Shop keeper')]
	df[, FE_OC_OTH:= as.integer(FEMALE_OCAT%in%c('zother','Student'))]		
	df[, RFM_ABOVE_MEAN:=  as.integer( RFM>=subset(dgg, variable=='RFM')[, mean(M)] ) ]
	
	
	#
	#	STAN
	#
	
	#	same household definitely stronger effect than couples	
	ms.1 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + 
							# comm_fish*COMM_FISH is base 
							comm_agr*COMM_AGR + comm_trad*COMM_TRAD + comm_mxd*COMM_MXD +
							noedu_f*FE_NOEDU + noedu_m*MA_NOEDU +
							sexp_f*FE_SEXP1YR_G1 + sexp_m*MA_SEXP1YR_G1 +
							circum_m*MA_CIRCUM +	
							prevratio_fm * RFM_ABOVE_MEAN +
							# agro occupation is baseline
							mocc_other*MA_OC_OTH + mocc_fish*MA_OC_FISH + mocc_trad*MA_OC_TRAD +										
							focc_other*FE_OC_OTH + focc_bar*FE_OC_BAR +  focc_trad*FE_OC_TRAD,
					base ~ dnorm(0,100),								
					c(comm_agr, comm_trad, comm_mxd) ~ dnorm(0,10),
					c(noedu_f, noedu_m) ~ dnorm(0,10),								
					c(sexp_f, sexp_m, circum_m, prevratio_fm) ~ dnorm(0,10),					
					c(mocc_other, mocc_fish, mocc_trad) ~ dnorm(0,10),
					c(focc_other, focc_bar, focc_trad) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, comm_agr=0, comm_trad=0, comm_mxd=0, noedu_f=0, noedu_m=0, sexp_f=0, sexp_m=0, circum_m=0,						
						prevratio_fm=0, focc_other=0, focc_bar=0, focc_trad=0, mocc_other=0, mocc_fish=0, mocc_trad=0),			
			warmup=1e3, iter=5e3, chains=1, cores=4
		)		
	#plot(precis(ms.1, prob=0.95))
	#pairs(ms.1)
	post	<- extract.samples(ms.1)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_COMM_FISH=logistic( post$base ),				
			PI_COMM_AGRO=logistic( post$base+post$comm_agr ), 											
			PI_COMM_TRAD=logistic( post$base+post$comm_trad ),							
			PI_COMM_MXD=logistic( post$base+post$comm_mxd ),
			PI_FEDU_NONE=logistic( post$base+post$noedu_f ),
			PI_MEDU_NONE=logistic( post$base+post$noedu_m ),							
			PI_FESEXP1YR_G1=logistic( post$base+post$sexp_f ),
			PI_MASEXP1YR_G1=logistic( post$base+post$sexp_m ),							
			PI_CIRC_YES=logistic( post$base+post$circum_m ),							
			PI_MOCC_OTH=logistic( post$base+post$mocc_other ),
			PI_MOCC_FISH=logistic( post$base+post$mocc_fish ),
			PI_MOCC_TRAD=logistic( post$base+post$mocc_trad ),
			PI_FOCC_OTH=logistic( post$base+post$focc_other ),
			PI_FOCC_BAR=logistic( post$base+post$focc_bar ),
			PI_FOCC_TRAD=logistic( post$base+post$focc_trad ),
			PI_PREV_FMRATIO=logistic( post$base+post$prevratio_fm ),
			OR_COMM_FISH=exp( post$base ),				
			OR_COMM_AGRO=exp( post$base+post$comm_agr ), 											
			OR_COMM_TRAD=exp( post$base+post$comm_trad ),							
			OR_COMM_MXD=exp( post$base+post$comm_mxd ),
			OR_FEDU_NONE=exp( post$base+post$noedu_f ),
			OR_MEDU_NONE=exp( post$base+post$noedu_m ),							
			OR_FESEXP1YR_G1=exp( post$base+post$sexp_f ),
			OR_MASEXP1YR_G1=exp( post$base+post$sexp_m ),
			OR_CIRC_YES=exp( post$base+post$circum_m ),							
			OR_MOCC_OTH=exp( post$base+post$mocc_other ),
			OR_MOCC_FISH=exp( post$base+post$mocc_fish ),
			OR_MOCC_TRAD=exp( post$base+post$mocc_trad ),
			OR_FOCC_OTH=exp( post$base+post$focc_other ),
			OR_FOCC_BAR=exp( post$base+post$focc_bar ),
			OR_FOCC_TRAD=exp( post$base+post$focc_trad ),	
			OR_PREV_FMRATIO=exp( post$base+post$prevratio_fm ),
			ORX_COMM_AGRO=exp( post$comm_agr ), 											
			ORX_COMM_TRAD=exp( post$comm_trad ),							
			ORX_COMM_MXD=exp( post$comm_mxd ),
			ORX_FEDU_NONE=exp( post$noedu_f ),
			ORX_MEDU_NONE=exp( post$noedu_m ),
			ORX_FESEXP1YR_G1=exp( post$sexp_f ),
			ORX_MASEXP1YR_G1=exp( post$sexp_m ),							
			ORX_CIRC_YES=exp( post$circum_m ),							
			ORX_MOCC_OTH=exp( post$mocc_other ),
			ORX_MOCC_FISH=exp( post$mocc_fish ),
			ORX_MOCC_TRAD=exp( post$mocc_trad ),
			ORX_FOCC_OTH=exp( post$focc_other ),
			ORX_FOCC_BAR=exp( post$focc_bar ),
			ORX_FOCC_TRAD=exp( post$focc_trad ),
			ORX_PREV_FMRATIO=logistic( post$prevratio_fm )
	)
	dp		<- melt(dp, id.vars='MC')	
	dss.1	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	dss.1	<- dcast.data.table(dss.1, variable~STAT, value.var='V')
	dss.1[, MODEL:= 'multivariate ms2']
	ds		<- copy(dss.1)
	ds[, LABEL:= paste0(round(MED, d=2), '\n[', round(CL, d=2),'-', round( CU, d=2),']')]
	ds[, STAT:= factor(gsub('^([^_]+)_.*','\\1',variable), levels=c('PI','OR','ORX'), labels=c('proportion MF','odds MF','odds ratio'))]
	ds[, FACTOR:=gsub('^([^_]+)_([^_]+)_([^_]+)','\\3',variable)]
	ds[, MOFA:=gsub('^([^_]+)_([^_]+)_([^_]+)','\\2-\\3',variable)]
	ds[, MODELTYPE:= factor(grepl('univariate',MODEL), levels=c(TRUE,FALSE),labels=c('univariate','multivariate'))]	
	ds		<- subset(ds, STAT=='odds ratio')
	set(ds, NULL, 'MOFA2', ds[, factor(MOFA, levels=rev(c(
											"PREV-FMRATIO",
											"COMM-AGRO","COMM-TRAD","COMM-MXD",																						
											"FEDU-NONE", "MEDU-NONE",      
											"FESEXP1YR-G1","MASEXP1YR-G1", 
											"CIRC-YES",                   
											"FOCC-BAR","FOCC-TRAD","FOCC-OTH",
											"MOCC-FISH","MOCC-OTH","MOCC-TRAD"   
									)), labels=rev(c("female-to-male prevalence ratio above average vs. below average",
											"Agrarian community vs. fishing community", "Trading community vs. fishing community", "Mixed vs. fishing community",
											"Female education: None vs. at least primary education", "Male education: None vs. at least primary education",
											"Female sex partners in last year: >1 vs. 1", "Male sex partners in last year: >1 vs. 1",
											"Male circumcised: yes vs. no",
											"Female occupation: bar/waitress vs. agricultural/house", "Female occupation: trading/shopkeeper vs. agricultural/house","Female occupation: other vs. agricultural/house", 
											"Male occupation: fishing vs. agricultural/house", "Male occupation: other vs. agricultural/house", "Male occupation: trading/shopkeeper vs. agricultural/house"
									)))])
	setkey(ds, MOFA)	
	dp		<- subset(ds, !MOFA%in%c('COMM-TRAD','COMM-MXD'))
	dp[, FILL:= '0']
	set(dp, dp[, which(IL>1 | IU<1)], 'FILL', '1')
	set(dp, dp[, which(CL>1 | CU<1)], 'FILL', '2')
	ggplot(dp, aes(x=MOFA2)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/20, 1/10, 1/6, 1/3,1/2,1,2,3, 6,10,20), labels=c('1/20','1/10','1/6','1/3','1/2','1','2','3','6','10', '20'), expand=c(0,0)) +
			scale_fill_manual(values=c('0'='white', '1'='orange', '2'='red')) +
			coord_flip(ylim=c(1/20,20)) +
			guides(fill=FALSE) +
			labs(x='Partner characteristics\n', y='\nOdds ratio of male-to-female transmission\nbetween partner characteristics')
	ggsave(file=paste0(outfile.base,'_propmf_factors_univariate_samehh.pdf'), w=8, h=3.5)		
}

RakaiFull.gender.171122.hh.multivariatemodels.stan.with.threshold.v2<- function()
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
	
	#zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	#zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	load(infile)	
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	rtpdm[, PAIR_COMM_TYPE:= FEMALE_COMM_TYPE]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_TYPE!=MALE_COMM_TYPE)], 'PAIR_COMM_TYPE', 'mixed')
	set(rtpdm, rtpdm[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(rtpdm, rtpdm[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')
	
	#rtpdm[, table(PAIR_COMM_TYPE)]
	rtpdm[, COUPLE2:= factor(COUPLE=='no couple', levels=c(TRUE,FALSE), labels=c('no couple','couple'))]
	rtpdm[, SAMEHH:= factor(FEMALE_HH_NUM==MALE_HH_NUM, levels=c(TRUE,FALSE), labels=c('same hh','different hh'))]	
	rtpdm[, PAIR_COMM:= MALE_COMM_NUM_A]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_NUM_A!=MALE_COMM_NUM_A)], 'PAIR_COMM', 'mixed')
	rtpdm[, MALE_SEXP1OUT2:= factor(MALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(MALE_SEXP1OUT=='Unknown')], 'MALE_SEXP1OUT2', 'Unknown')
	rtpdm[, FEMALE_SEXP1OUT2:= factor(FEMALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(FEMALE_SEXP1OUT=='Unknown')], 'FEMALE_SEXP1OUT2', 'Unknown')
	
	#	add immigrant
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender_age.rda")
	tmp		<- subset(de, !is.na(RID), select=c(RID, INMIGRANT))
	infile				<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda"
	load(infile)
	tmp2	<- subset(as.data.table(inmigrant), select=c(RCCS_studyid,inmigrant))
	setnames(tmp2, 'RCCS_studyid', 'RID')
	tmp2	<- tmp2[, list(inmigrant=as.integer(any(inmigrant==1))), by='RID']
	tmp		<- merge(tmp, tmp2, all=TRUE, by='RID')
	#INMIGRANT and inmigrant do not agree.. 
	#subset(tmp, !is.na(inmigrant) & !is.na(INMIGRANT))[, table(inmigrant, INMIGRANT)]
	#...just set those with missing INMIGRANT to inmigrant
	set(tmp, NULL, 'INMIGRANT', tmp[, as.integer(INMIGRANT)])
	tmp2	<- tmp[, which(!is.na(inmigrant) & is.na(INMIGRANT))]
	set(tmp, tmp2, 'INMIGRANT', tmp[tmp2, inmigrant])	
	tmp		<- subset(tmp, select=c(RID, INMIGRANT))	
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	rtpdm	<- merge(rtpdm, tmp, by='MALE_RID', all.x=TRUE)
	setnames(tmp, colnames(tmp), gsub('MALE_','FEMALE_',colnames(tmp)))
	rtpdm	<- merge(rtpdm, tmp, by='FEMALE_RID', all.x=TRUE)
	set(rtpdm, rtpdm[, which(is.na(MALE_INMIGRANT))], 'MALE_INMIGRANT', 0L)
	set(rtpdm, rtpdm[, which(is.na(FEMALE_INMIGRANT))], 'FEMALE_INMIGRANT', 0L)
	
	
	df		<- subset(rtpdm, FEMALE_HH_NUM==MALE_HH_NUM, select=c(	MALE_RID, FEMALE_RID, PTY_RUN, IDCLU, LINK_MF, POSTERIOR_SCORE_MF, COUPLE2, SAMEHH, PAIR_COMM_TYPE, FEMALE_COMM_NUM,
					MALE_RECENTVL, MALE_SEXP1YR, MALE_SEXP1OUT2, MALE_OCCUP_OLLI, MALE_OCAT, MALE_EDUCAT, MALE_CIRCUM, MALE_INMIGRANT,
					FEMALE_RECENTVL, FEMALE_SEXP1YR, FEMALE_SEXP1OUT2, FEMALE_OCCUP_OLLI, FEMALE_OCAT, FEMALE_EDUCAT, FEMALE_INMIGRANT 
			))
	#	missing data: fill in
	set(df, df[, which(is.na(MALE_RECENTVL))], 'MALE_RECENTVL', -1)
	set(df, df[, which(is.na(FEMALE_RECENTVL))], 'FEMALE_RECENTVL', -1)
	set(df, df[, which(is.na(MALE_CIRCUM))], 'MALE_CIRCUM', 'Unknown')
	set(df, df[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(df, df[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')	
	#for(x in colnames(df)) print( c(x, any(is.na(df[[x]]))) )
	#	add estimated HIV prevalences (medians)
	infile	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30trmMF_estimated_prevalenceratio.rda"	
	load(infile)
	tmp		<- dcast.data.table(dgg, COMM_NUM~variable, value.var='M')
	setnames(tmp, 'COMM_NUM', 'FEMALE_COMM_NUM')
	df		<- merge(df, tmp, by='FEMALE_COMM_NUM')
	
	#	prepare data for STAN
	df[, COMM_FISH:= as.integer(PAIR_COMM_TYPE=='fisherfolk')]
	df[, COMM_AGR:= as.integer(PAIR_COMM_TYPE=='agrarian')]
	df[, COMM_TRAD:= as.integer(PAIR_COMM_TYPE=='trading')]
	df[, COMM_MXD:= as.integer(PAIR_COMM_TYPE=='mixed')]
	df[, FE_NOEDU:= as.integer(FEMALE_EDUCAT=='None')]
	df[, FE_NOEDU_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, MA_NOEDU:= as.integer(MALE_EDUCAT=='None')]
	df[, MA_NOEDU_MISS:= as.integer(MALE_EDUCAT=='Unknown')]	
	df[, MA_CIRCUM:= as.integer(MALE_CIRCUM=='Y')]
	df[, MA_CIRCUM_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, FE_SEXP1YR_G1:= as.integer(FEMALE_SEXP1YR!='1')]
	df[, MA_SEXP1YR_G1:= as.integer(MALE_SEXP1YR!='1')]
	df[, MA_OC_AGRO:= as.integer(MALE_OCAT=='Agro/House')]	
	df[, MA_OC_FISH:= as.integer(MALE_OCAT=='Fishing')]
	df[, MA_OC_TRAD:= as.integer(MALE_OCAT=='Trading/Shop keeper')]
	df[, MA_OC_OTH:= as.integer(MALE_OCAT%in%c('zother','Bar/waitress','Student','Boda/Trucking'))]	
	df[, FE_OC_AGRO:= as.integer(FEMALE_OCAT=='Agro/House')]
	df[, FE_OC_BAR:= as.integer(FEMALE_OCAT=='Bar/waitress')]			
	df[, FE_OC_TRAD:= as.integer(FEMALE_OCAT=='Trading/Shop keeper')]
	df[, FE_OC_OTH:= as.integer(FEMALE_OCAT%in%c('zother','Student'))]		
	df[, RFM_ABOVE_MEAN:=  as.integer( RFM>=subset(dgg, variable=='RFM')[, mean(M)] ) ]
	df[, MA_VL_HIGH:= as.integer(MALE_RECENTVL>=1e5)]
	df[, FE_VL_HIGH:= as.integer(FEMALE_RECENTVL>=1e5)]
	#
	#	STAN
	#
	mhh.1 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base +
							male_vlhigh*MA_VL_HIGH + female_vlhigh*FE_VL_HIGH +
							male_noedu*MA_NOEDU + female_noedu*FE_NOEDU +
							male_inmigrant*MALE_INMIGRANT + female_inmigrant*FEMALE_INMIGRANT +							
							female_sexp*FE_SEXP1YR_G1 + male_sexp*MA_SEXP1YR_G1 +															 														
							# occupation: other is baseline
							male_agro*MA_OC_AGRO + male_fish*MA_OC_FISH + male_shopkeeper*MA_OC_TRAD +  
							female_agro*FE_OC_AGRO + female_barworker*FE_OC_BAR + female_shopkeeper*FE_OC_TRAD,
					base ~ dnorm(0,100),													
					c(male_noedu, female_noedu, male_vlhigh, female_vlhigh) ~ dnorm(0, 1),								
					c(male_inmigrant, female_inmigrant, male_sexp, female_sexp) ~ dnorm(0, 1),					
					c(male_fish, male_shopkeeper, male_agro, female_barworker, female_shopkeeper, female_agro) ~ dnorm(0, 1)					
			),
			data=as.data.frame(df), 
			start=list(	base=0, male_noedu=0, female_noedu=0, male_inmigrant=0, female_inmigrant=0, male_sexp=0, female_sexp=0, male_vlhigh=0, female_vlhigh=0,						
						male_fish=0, male_shopkeeper=0, male_agro=0, female_barworker=0, female_shopkeeper=0, female_agro=0),			
			warmup=1e3, iter=1e4, chains=1, cores=4
		)
	precis(mhh.1, prob=0.95)	
	plot(precis(mhh.1, prob=0.95))
	#	get stuff from mhh.1
	post	<- as.data.table(extract.samples(mhh.1))
	post[, MC:= seq_len(nrow(post))]
	post	<- melt(post, id.vars='MC')
	tmp		<- post[, which(!grepl('sig',variable))]
	set(post, tmp, 'value', post[tmp, exp(value)])	
	tmp		<- post[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), 
			by=c('variable')]
	dhh1	<- dcast.data.table(tmp, variable~STAT, value.var='V')
	dhh1[, MODEL:= 'mf in hh']
	dhh1[, STAT:= as.character(factor(grepl('base',variable), levels=c(TRUE,FALSE), labels=c('odds','odds ratio')))]
	set(dhh1, NULL, 'FACTOR', dhh1[, factor(variable, levels=c('base','male_vlhigh','female_vlhigh','male_noedu','female_noedu','male_inmigrant','female_inmigrant',
																'male_sexp','female_sexp','male_fish','male_shopkeeper','male_agro','female_barworker','female_shopkeeper','female_agro'))])
	setkey(dhh1, FACTOR)
	#	save
	save(df, mhh.1, dhh1, file=paste0(outfile.base,'_interaction_runs_mf_followup.rda'))	
	
	#	make plots
	dp		<- subset(dhh1, !grepl('base',variable))	
	dp[, FILL:= '0']
	set(dp, dp[, which(IU<1)], 'FILL', '-1')
	set(dp, dp[, which(CU<1)], 'FILL', '-2')	
	set(dp, dp[, which(IL>1)], 'FILL', '1')
	set(dp, dp[, which(CL>1)], 'FILL', '2')
	ggplot(dp, aes(x=FACTOR)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/20, 1/10, 1/6, 1/3,1/2,1,2,3, 6,10,20), labels=c('1/20','1/10','1/6','1/3','1/2','1','2','3','6','10', '20'), expand=c(0,0)) +			
			scale_fill_manual(values=c('-2'='deepskyblue','-1'='cyan','0'='grey80', '1'='darkorange', '2'='red')) +
			#coord_flip(ylim=c(1/20,20)) +
			coord_flip(ylim=c(1/10,10)) +
			guides(fill=FALSE) +
			labs(x='Partner characteristics\n', y='\nOdds for male-to-female transmission\nin phylogenetically reconstructed transmission pairs') +
			facet_grid(STAT~., scales='free', space='free')	
	ggsave(file=paste0(outfile.base,'_interaction_runs_odds_noage_mf.pdf'), w=5, h=5)	
}

RakaiFull.transmitter.171122.get.data.set<- function()
{
	require(data.table)
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_withmetadata.rda"
	outfile.base			<- gsub('_withmetadata.rda','',infile)
	
	#zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	#zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
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
	rtpdm[, MALE_SEX:='M']
	rtpdm[, FEMALE_SEX:='F']
	#	add immigrant
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender_age.rda")
	tmp		<- subset(de, !is.na(RID), select=c(RID, INMIGRANT))
	infile				<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda"
	load(infile)
	tmp2	<- subset(as.data.table(inmigrant), select=c(RCCS_studyid,inmigrant))
	setnames(tmp2, 'RCCS_studyid', 'RID')
	tmp2	<- tmp2[, list(inmigrant=as.integer(any(inmigrant==1))), by='RID']
	tmp		<- merge(tmp, tmp2, all=TRUE, by='RID')
	#INMIGRANT and inmigrant do not agree.. 
	#subset(tmp, !is.na(inmigrant) & !is.na(INMIGRANT))[, table(inmigrant, INMIGRANT)]
	#...just set those with missing INMIGRANT to inmigrant
	set(tmp, NULL, 'INMIGRANT', tmp[, as.integer(INMIGRANT)])
	tmp2	<- tmp[, which(!is.na(inmigrant) & is.na(INMIGRANT))]
	set(tmp, tmp2, 'INMIGRANT', tmp[tmp2, inmigrant])	
	tmp		<- subset(tmp, select=c(RID, INMIGRANT))	
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	rtpdm	<- merge(rtpdm, tmp, by='MALE_RID', all.x=TRUE)
	setnames(tmp, colnames(tmp), gsub('MALE_','FEMALE_',colnames(tmp)))
	rtpdm	<- merge(rtpdm, tmp, by='FEMALE_RID', all.x=TRUE)
	set(rtpdm, rtpdm[, which(is.na(MALE_INMIGRANT))], 'MALE_INMIGRANT', 0L)
	set(rtpdm, rtpdm[, which(is.na(FEMALE_INMIGRANT))], 'FEMALE_INMIGRANT', 0L)
	setkey(rtpdm, MALE_RID,FEMALE_RID)
	rtpdm[, PAIRID:= seq_len(nrow(rtpdm))]
	
	df		<- subset(rtpdm, select=c(	MALE_RID, FEMALE_RID, PAIRID, PTY_RUN, LINK_MF, COUPLE2, SAMEHH, PAIR_COMM_TYPE, 					
					MALE_SEX, MALE_BIRTHDATE, MALE_COMM_TYPE, MALE_RECENTVL, MALE_SEXP1YR, MALE_SEXP1OUT2, MALE_OCCUP_OLLI, MALE_OCAT, MALE_EDUCAT, FEMALE_INMIGRANT, MALE_CIRCUM,
					FEMALE_SEX, FEMALE_BIRTHDATE, FEMALE_COMM_TYPE, FEMALE_RECENTVL, FEMALE_SEXP1YR, FEMALE_SEXP1OUT2, FEMALE_OCCUP_OLLI, FEMALE_OCAT, FEMALE_EDUCAT, MALE_INMIGRANT 
					))
	#	missing data: fill in
	set(df, df[, which(is.na(MALE_RECENTVL))], 'MALE_RECENTVL', -1)
	set(df, df[, which(is.na(FEMALE_RECENTVL))], 'FEMALE_RECENTVL', -1)
	set(df, df[, which(is.na(MALE_CIRCUM))], 'MALE_CIRCUM', 'Unknown')
	set(df, df[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(df, df[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')
	set(df, df[, which(is.na(MALE_BIRTHDATE))], 'MALE_BIRTHDATE', df[, mean(MALE_BIRTHDATE, na.rm=TRUE)])
	set(df, df[, which(is.na(FEMALE_BIRTHDATE))], 'FEMALE_BIRTHDATE', df[, mean(FEMALE_BIRTHDATE, na.rm=TRUE)])	
	df[, MALE_AGE_AT_MID:= 2013.25-MALE_BIRTHDATE]
	df[, FEMALE_AGE_AT_MID:= 2013.25-FEMALE_BIRTHDATE]
	#for(x in colnames(df)) print( c(x, any(is.na(df[[x]]))) )
	#	add estimated HIV prevalences (medians)
	#infile	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30trmMF_estimated_prevalenceratio.rda"	
	#load(infile)
	#tmp		<- dcast.data.table(dgg, COMM_NUM~variable, value.var='M')
	#setnames(tmp, 'COMM_NUM', 'FEMALE_COMM_NUM')
	#df		<- merge(df, tmp, by='FEMALE_COMM_NUM')
	
	#	cast to TRM==1 and REC==0
	dg		<- subset(df, LINK_MF==1)
	set(dg, NULL, colnames(dg)[grepl('FEMALE',colnames(dg))], NULL)
	setnames(dg, colnames(dg), gsub('^MALE_','',colnames(dg)))
	tmp		<- subset(df, LINK_MF==1)
	set(tmp, NULL, colnames(tmp)[grepl('^MALE',colnames(tmp))], NULL)
	setnames(tmp, colnames(tmp), gsub('^FEMALE_','',colnames(tmp)))
	tmp[, LINK_MF:=0L]
	dg		<- rbind(dg, tmp, fill=TRUE)	
	tmp		<- subset(df, LINK_MF==0)
	set(tmp, NULL, colnames(tmp)[grepl('FEMALE',colnames(tmp))], NULL)
	setnames(tmp, colnames(tmp), gsub('^MALE_','',colnames(tmp)))	
	dg		<- rbind(dg, tmp, fill=TRUE)
	tmp		<- subset(df, LINK_MF==0)
	set(tmp, NULL, colnames(tmp)[grepl('^MALE',colnames(tmp))], NULL)
	setnames(tmp, colnames(tmp), gsub('^FEMALE_','',colnames(tmp)))
	tmp[, LINK_MF:=1L]
	dg		<- rbind(dg, tmp, fill=TRUE)
	setnames(dg, 'LINK_MF','TRANSMITTER')
	dg[, AGE_AT_MID_C:= dg[,as.character(cut(AGE_AT_MID, breaks=c(15,20,25,30,35,40,45,52), right=FALSE, levels=c('15-19','20-24','25-29','30-34','35-39','40-44','45-50')))]]	
	df		<- copy(dg)
	
	save(df, file= paste0(outfile.base, '_transmitterrecipientdata.rda'))
}

RakaiFull.transmitter.171122.gender.explore<- function()
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
	
	infile			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_transmitterrecipientdata.rda"
	outfile.base	<- gsub('data.rda','',infile)
	load(infile)
	#
	#	prepare variables for STAN
	#		
	dg	<- copy(df)
	#	prepare data for STAN
	dg[, COMM_FISH:= as.integer(COMM_TYPE=='fisherfolk')]
	dg[, COMM_AGR:= as.integer(COMM_TYPE=='agrarian')]
	dg[, COMM_TRAD:= as.integer(COMM_TYPE=='trading')]
	#dg[, COMM_MXD:= as.integer(PAIR_COMM_TYPE=='mixed')]
	dg[, HH:= as.integer(SAMEHH=='same hh')]
	dg[, MALE:= as.integer(SEX=='M')]
	dg[, NOEDU:= as.integer(EDUCAT=='None')]
	dg[, NOEDU_MISS:= as.integer(EDUCAT=='Unknown')]
	dg[, CIRCUM2:= as.integer(CIRCUM=='Y')]
	dg[, SEXP1YR_G1:= as.integer(SEXP1YR!='1')]
	dg[, VL_HIGH:= as.integer(RECENTVL>=1e5)]
	#dg[, VL_SUPP:= as.integer(RECENTVL<=200)]
	dg[, OC_AGRO:= as.integer(OCAT=='Agro/House')]
	dg[, OC_BAR:= as.integer(OCAT=='Bar/waitress')]
	dg[, OC_BODA:= as.integer(OCAT=='Boda/Trucking')]
	dg[, OC_FISH:= as.integer(OCAT=='Fishing')]
	dg[, OC_STUD:= as.integer(OCAT=='Student')]
	dg[, OC_TRAD:= as.integer(OCAT=='Trading/Shop keeper')]
	dg[, OC_OTH:= as.integer(OCAT%in%c('zother'))]	
	dg[, AGE_AT_MID_C2:= dg[,as.integer(cut(AGE_AT_MID, breaks=c(15,20,25,30,35,40,45,52), right=FALSE, levels=c('1','2','3','4','5','6','7')))]]
	dg[, AGE_AT_MID_D:= dg[,as.character(cut(AGE_AT_MID, breaks=c(15,20,25,30,35,52), right=FALSE, levels=c('15-19','20-24','25-29','30-34','35-50')))]]
	dg[, AGE_AT_MID_D2:= dg[,as.integer(cut(AGE_AT_MID, breaks=c(15,20,25,30,35,52), right=FALSE, levels=c('1','2','3','4','5')))]]
	#
	
	dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('SAMEHH','SEX')]
	#	         SAMEHH SEX        PT        OT   N
	#1:      same hh   M 0.6666667 2.0000000 126
	#2: different hh   M 0.5329341 1.1410256 167
	#3:      same hh   F 0.3333333 0.5000000 126
	#4: different hh   F 0.4670659 0.8764045 167

	tmp	<- dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('SAMEHH','SEX','AGE_AT_MID_C2','AGE_AT_MID_C')]
	setkey(tmp, SEX, SAMEHH, AGE_AT_MID_C2)
	
	tmp	<- dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('SAMEHH','SEX','AGE_AT_MID_D2','AGE_AT_MID_D')]
	setkey(tmp, SEX, SAMEHH, AGE_AT_MID_D2)
	#	this is great and interesting:
	#	females within HH: bimodal (reaching 1.0); whereas outside HH: increasing
	#	males within HH: pretty constant; whereas outside HH: increasing

	
	#
	#	VIRAL LOAD (yes)
	#
	tmp	<- dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('SAMEHH','SEX','VL_HIGH')]
	setkey(tmp, VL_HIGH, SEX, SAMEHH)
	tmp
	dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('VL_HIGH')]
	#	   VL_HIGH        PT       OT   N
	#1:       0 0.4878049 0.952381 533
	#2:       1 0.6226415 1.650000  53
	#	--> that s great

	#
	#	CIRCUMCISION (no)
	#

	dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('SEX','CIRCUM')]
	#	   SEX  CIRCUM        PT        OT   N
	#1:   M       N 0.6028037 1.5176471 214
	#2:   M       Y 0.5512821 1.2285714  78
	#	--> suggests not to include circumcision

	#
	#	INMIGRANT ( don t include or include with gender)
	#

	dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('INMIGRANT')]
	#	   INMIGRANT        PT       OT   N
	#1:         0 0.5066667 1.027027 450
	#2:         1 0.4779412 0.915493 136	--> not much at all
	tmp	<- dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('SAMEHH','INMIGRANT')]
	setkey(tmp, INMIGRANT, SAMEHH)
	tmp
	#	         SAMEHH INMIGRANT        PT        OT   N
	#1:      same hh         0 0.5130890 1.0537634 191
	#2: different hh         0 0.5019305 1.0077519 259
	#3:      same hh         1 0.4590164 0.8484848  61		--> could be sig
	#4: different hh         1 0.4933333 0.9736842  75
	tmp	<- dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('SEX','INMIGRANT')]
	setkey(tmp, SEX, INMIGRANT)
	tmp		
	#	   SEX INMIGRANT        PT       OT   N
	#1:   F         0 0.4186047 0.720000 215
	#2:   F         1 0.3846154 0.625000  78	--> not much down
	#3:   M         0 0.5872340 1.422680 235
	#4:   M         1 0.6034483 1.521739  58	--> not much up
	#	inmigrant: overall not very strong 
	tmp	<- dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('SEX','SAMEHH','INMIGRANT')]
	setkey(tmp, SEX, SAMEHH, INMIGRANT)
	tmp		
	#	  SEX     SAMEHH INMIGRANT        PT        OT   N
	#1:   F      same hh         0 0.3111111 0.4516129  90
	#2:   F      same hh         1 0.3888889 0.6363636  36
	#3:   F different hh         0 0.4960000 0.9841270 125
	#4:   F different hh         1 0.3809524 0.6153846  42	--> could be sig
	#5:   M      same hh         0 0.6930693 2.2580645 101
	#6:   M      same hh         1 0.5600000 1.2727273  25	--> more down
	#7:   M different hh         0 0.5074627 1.0303030 134
	#8:   M different hh         1 0.6363636 1.7500000  33	--> more up

	#
	#	NO EDU (include with gender)
	#

	tmp	<- dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('SEX','SAMEHH','NOEDU')]
	setkey(tmp, SEX, SAMEHH, NOEDU)
	tmp	
	dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('SEX','NOEDU')]
	#	   SEX NOEDU        PT        OT   N
	#1:   M     0 0.5886792 1.4311927 265
	#2:   M     1 0.6071429 1.5454545  28
	#3:   F     0 0.4220532 0.7302632 263
	#4:   F     1 0.3000000 0.4285714  30			--> interesting
	#	no education shows up among females, sample size but OK ish

	#
	#	SEXP1YR_G1 (include in interaction with gender and household)
	#
	dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('SEXP1YR_G1')]
	#	   SEXP1YR_G1        PT        OT   N
	#1:          0 0.4634146 0.8636364 328
	#2:          1 0.5465116 1.2051282 258	--> this should be up and sig
	dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('SEX','SEXP1YR_G1')]
	#	   SEX SEXP1YR_G1        PT        OT   N
	#1:   M          0 0.5789474 1.3750000 114
	#2:   M          1 0.5977654 1.4861111 179
	#3:   F          0 0.4018692 0.6718750 214
	#4:   F          1 0.4303797 0.7555556  79	--> so the higher overall result is just confounded by more men reporting >1 sex partner
	tmp	<- dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('SEX','SEXP1YR_G1','SAMEHH')]
	setkey(tmp, SEX, SAMEHH, SEXP1YR_G1)
	tmp	
	#1:   F          0      same hh 0.3364486 0.5070423 107
	#2:   F          1      same hh 0.3157895 0.4615385  19
	#3:   F          0 different hh 0.4672897 0.8771930 107
	#4:   F          1 different hh 0.4666667 0.8750000  60
	#5:   M          0      same hh 0.6129032 1.5833333  62
	#6:   M          1      same hh 0.7187500 2.5555556  64	--> this should be sig
	#7:   M          0 different hh 0.5384615 1.1666667  52
	#8:   M          1 different hh 0.5304348 1.1296296 115

	#
	#	OCCUPATION
	#
	tmp	<- dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('SEX','OCAT')]
	setkey(tmp, SEX, OCAT)
	tmp	
	#SEX                OCAT        PT        OT   N
	#1:   F          Agro/House 0.3931034 0.6477273 145
	#2:   F        Bar/waitress 0.4878049 0.9523810  41 --> interesting, higher
	#3:   F             Student 0.2000000 0.2500000   5	--> put into other
	#4:   F Trading/Shop keeper 0.3888889 0.6363636  54
	#5:   F              zother 0.4375000 0.7777778  48
	#6:   M          Agro/House 0.6428571 1.8000000  42
	#7:   M        Bar/waitress 0.0000000 0.0000000   2 --> put into other
	#8:   M       Boda/Trucking 0.4285714 0.7500000   7 --> put into other
	#9:   M             Fishing 0.5606061 1.2758621 132 --> interesting, lower
	#10:   M             Student 0.5714286 1.3333333   7 --> put into other
	#11:   M Trading/Shop keeper 0.6136364 1.5882353  44
	#12:   M              zother 0.6440678 1.8095238  59

	#
	#	explore INMIGRANT COMM_TYPE
	#
	tmp	<- dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('INMIGRANT','COMM_TYPE')]
	setkey(tmp, INMIGRANT, COMM_TYPE)
	tmp	
	#	   INMIGRANT  COMM_TYPE        PT        OT   N
	#1:         0   agrarian 0.4822695 0.9315068 141
	#2:         0 fisherfolk 0.5070922 1.0287770 282
	#3:         0    trading 0.6296296 1.7000000  27
	#4:         1   agrarian 0.5208333 1.0869565  48
	#5:         1 fisherfolk 0.4444444 0.8000000  81
	#6:         1    trading 0.5714286 1.3333333   7
	tmp	<- dg[, list(PT= length(which(TRANSMITTER==1))/length(TRANSMITTER), OT= length(which(TRANSMITTER==1))/length(which(TRANSMITTER==0)), N=length(TRANSMITTER) ), by=c('SEX','INMIGRANT','COMM_TYPE')]
	setkey(tmp, SEX, INMIGRANT, COMM_TYPE)
	tmp	

}

RakaiFull.transmitter.171122.gender.interactionmodel.stan.with.threshold<- function()
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
	
	infile			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_transmitterrecipientdata.rda"
	outfile.base	<- gsub('data.rda','',infile)
	load(infile)
	#
	#	prepare variables for STAN
	#		
	dg	<- copy(df)
	#	prepare data for STAN
	dg[, COMM_FISH:= as.integer(COMM_TYPE=='fisherfolk')]
	dg[, COMM_AGR:= as.integer(COMM_TYPE=='agrarian')]
	dg[, COMM_TRAD:= as.integer(COMM_TYPE=='trading')]
	#dg[, COMM_MXD:= as.integer(PAIR_COMM_TYPE=='mixed')]
	dg[, HH:= as.integer(SAMEHH=='same hh')]
	dg[, MALE:= as.integer(SEX=='M')]
	dg[, NOEDU:= as.integer(EDUCAT=='None')]
	dg[, NOEDU_MISS:= as.integer(EDUCAT=='Unknown')]
	dg[, SEXP1YR_G1:= as.integer(SEXP1YR!='1')]
	dg[, VL_HIGH:= as.integer(RECENTVL>=1e5)]
	dg[, MA_CIRCUM:= as.integer(!is.na(CIRCUM) & CIRCUM=='Y')]
	#dg[, VL_SUPP:= as.integer(RECENTVL<=200)]
	dg[, OC_AGRO:= as.integer(OCAT=='Agro/House')]
	dg[, OC_BAR:= as.integer(OCAT=='Bar/waitress')]
	dg[, OC_BODA:= as.integer(OCAT=='Boda/Trucking')]
	dg[, OC_FISH:= as.integer(OCAT=='Fishing')]
	dg[, OC_STUD:= as.integer(OCAT=='Student')]
	dg[, OC_TRAD:= as.integer(OCAT=='Trading/Shop keeper')]
	dg[, OC_OTH:= as.integer(OCAT%in%c('zother'))]	
	dg[, AGE_AT_MID_C2:= dg[,as.integer(cut(AGE_AT_MID, breaks=c(15,20,25,30,35,40,45,52), right=FALSE, levels=c('1','2','3','4','5','6','7')))]]
	dg[, AGE_AT_MID_D:= dg[,as.character(cut(AGE_AT_MID, breaks=c(15,20,25,30,35,52), right=FALSE, levels=c('15-19','20-24','25-29','30-34','35-50')))]]
	dg[, AGE_AT_MID_D2:= dg[,as.integer(cut(AGE_AT_MID, breaks=c(15,20,25,30,35,52), right=FALSE, levels=c('1','2','3','4','5')))]]
	dg[, MA_OC_AGRO:= as.integer(OCAT=='Agro/House')]	
	dg[, MA_OC_FISH:= as.integer(OCAT=='Fishing')]
	dg[, MA_OC_TRAD:= as.integer(OCAT=='Trading/Shop keeper')]
	dg[, MA_OC_OTH:= as.integer(OCAT%in%c('zother','Bar/waitress','Student','Boda/Trucking'))]	
	dg[, FE_OC_AGRO:= as.integer(OCAT=='Agro/House')]
	dg[, FE_OC_BAR:= as.integer(OCAT=='Bar/waitress')]			
	dg[, FE_OC_TRAD:= as.integer(OCAT=='Trading/Shop keeper')]
	dg[, FE_OC_OTH:= as.integer(OCAT%in%c('zother','Student'))]		
	
	
	
	#	model 1 interaction model 
	#	MALE:HH + VL_HIGH + COMM_TYPE + MALE:NOEDU + MALE:INMIGRANT + MALE:SEXP1YR_G1
	mi.1 <- map2stan(
			alist(
					TRANSMITTER ~ dbinom(1, ptr),
					logit(ptr) <- base + #	is female_in_hh
								  male_in_hh*MALE*HH + 
								  male_extra_hh*MALE*(1-HH) + female_extra_hh*(1-MALE)*(1-HH) +
								  vlhigh*VL_HIGH +
								  comm_fish*COMM_FISH + comm_trad*COMM_TRAD +
								  male_noedu*MALE*NOEDU + female_noedu*(1-MALE)*NOEDU +
								  male_inmigrant*MALE*INMIGRANT + female_inmigrant*(1-MALE)*INMIGRANT +
								  male_sexp*MALE*SEXP1YR_G1 + female_sexp*(1-MALE)*SEXP1YR_G1,
					base ~ dnorm(0,100),										
					c(comm_fish, comm_trad, vlhigh) ~ dnorm(0,10),
					c(male_in_hh, male_extra_hh, female_extra_hh, male_noedu, female_noedu, male_inmigrant, female_inmigrant, male_sexp, female_sexp) ~ dnorm(0,10)													
			),
			data=as.data.frame(dg), 
			start=list(	base=0, comm_fish=0, comm_trad=0, vlhigh=0,
						male_in_hh=0, male_extra_hh=0, female_extra_hh=0, male_noedu=0, female_noedu=0, male_inmigrant=0, female_inmigrant=0, male_sexp=0, female_sexp=0),			
			warmup=1e3, iter=1e4, chains=1, cores=4
		)
	precis(mi.1, prob=0.95)	
	plot(precis(mi.1, prob=0.95))
	
	mi.2 <- map2stan(
			alist(
					TRANSMITTER ~ dbinom(1, ptr),
					logit(ptr) <- male_in_hh*MALE*HH + female_in_hh*(1-MALE)*HH +
							male_extra_hh*MALE*(1-HH) + female_extra_hh*(1-MALE)*(1-HH) +
							vlhigh*VL_HIGH +							
							male_noedu*MALE*NOEDU + female_noedu*(1-MALE)*NOEDU +
							male_inmigrant*MALE*INMIGRANT + female_inmigrant*(1-MALE)*INMIGRANT +
							male_sexp*MALE*SEXP1YR_G1 + female_sexp*(1-MALE)*SEXP1YR_G1,															
					c(vlhigh) ~ dnorm(0,10),
					c(female_in_hh, male_in_hh, male_extra_hh, female_extra_hh, male_noedu, female_noedu, male_inmigrant, female_inmigrant, male_sexp, female_sexp) ~ dnorm(0,10)													
			),
			data=as.data.frame(dg), 
			start=list(	vlhigh=0,
						female_in_hh=0, male_in_hh=0, male_extra_hh=0, female_extra_hh=0, male_noedu=0, female_noedu=0, male_inmigrant=0, female_inmigrant=0, male_sexp=0, female_sexp=0),			
			warmup=1e3, iter=1e4, chains=1, cores=4
		)
	#	mi.2 is good
	mi.4 <- map2stan(
			alist(
					TRANSMITTER ~ dbinom(1, ptr),
					logit(ptr) <- male_in_hh*MALE*HH + female_in_hh*(1-MALE)*HH +
							male_extra_hh*MALE*(1-HH) + female_extra_hh*(1-MALE)*(1-HH) +
							vlhigh*VL_HIGH + 							
							male_noedu*MALE*NOEDU + female_noedu*(1-MALE)*NOEDU +
							male_inmigrant*MALE*INMIGRANT + female_inmigrant*(1-MALE)*INMIGRANT +
							male_sexp*MALE*SEXP1YR_G1 + female_sexp*(1-MALE)*SEXP1YR_G1 +
							# occupation: agro/house is baseline
							male_agro*MALE*MA_OC_AGRO + male_fish*MALE*MA_OC_FISH + male_shopkeeper*MALE*MA_OC_TRAD +  
							female_agro*(1-MALE)*FE_OC_AGRO + female_barworker*(1-MALE)*FE_OC_BAR + female_shopkeeper*(1-MALE)*FE_OC_TRAD,															
					c(vlhigh) ~ dnorm(0,10),
					c(male_fish, male_shopkeeper, male_agro, female_barworker, female_shopkeeper, female_agro) ~ dnorm(0,10),
					c(female_in_hh, male_in_hh, male_extra_hh, female_extra_hh, male_noedu, female_noedu, male_inmigrant, female_inmigrant, male_sexp, female_sexp) ~ dnorm(0,10)													
			),
			data=as.data.frame(dg), 
			start=list(	vlhigh=0,
						male_fish=0, male_shopkeeper=0, male_agro=0, female_barworker=0, female_shopkeeper=0, female_agro=0,
						female_in_hh=0, male_in_hh=0, male_extra_hh=0, female_extra_hh=0, male_noedu=0, female_noedu=0, male_inmigrant=0, female_inmigrant=0, male_sexp=0, female_sexp=0),			
			warmup=1e3, iter=1e4, chains=1, cores=4
		)
	plot(mi.4)
	precis(mi.4, depth=2, prob=0.95)
	plot(precis(mi.4, depth=2, prob=0.95))
	
	
	mi.5 <- map2stan(
			alist(
					TRANSMITTER ~ dbinom(1, ptr),
					logit(ptr) <- male_in_hh*MALE*HH + female_in_hh*(1-MALE)*HH +
							male_extra_hh*MALE*(1-HH) + female_extra_hh*(1-MALE)*(1-HH) +
							vlhigh*VL_HIGH + male_circ*MALE*MA_CIRCUM +							
							male_noedu*MALE*NOEDU + female_noedu*(1-MALE)*NOEDU +
							male_inmigrant*MALE*INMIGRANT + female_inmigrant*(1-MALE)*INMIGRANT +
							male_sexp*MALE*SEXP1YR_G1 + female_sexp*(1-MALE)*SEXP1YR_G1 +
							# occupation: agro/house is baseline
							male_agro*MALE*MA_OC_AGRO + male_fish*MALE*MA_OC_FISH + male_shopkeeper*MALE*MA_OC_TRAD +  
							female_agro*(1-MALE)*FE_OC_AGRO + female_barworker*(1-MALE)*FE_OC_BAR + female_shopkeeper*(1-MALE)*FE_OC_TRAD,															
					c(vlhigh, male_circ) ~ dnorm(0,10),
					c(male_fish, male_shopkeeper, male_agro, female_barworker, female_shopkeeper, female_agro) ~ dnorm(0,10),
					c(female_in_hh, male_in_hh, male_extra_hh, female_extra_hh, male_noedu, female_noedu, male_inmigrant, female_inmigrant, male_sexp, female_sexp) ~ dnorm(0,10)													
			),
			data=as.data.frame(dg), 
			start=list(	vlhigh=0, male_circ=0,
						female_in_hh=0, male_in_hh=0, male_extra_hh=0, female_extra_hh=0, 
						male_fish=0, male_shopkeeper=0, male_agro=0, female_barworker=0, female_shopkeeper=0, female_agro=0,
						male_noedu=0, female_noedu=0, male_inmigrant=0, female_inmigrant=0, male_sexp=0, female_sexp=0),			
			warmup=1e3, iter=1e4, chains=1, cores=4
		)
	precis(mi.5, depth=2, prob=0.95)
	plot(precis(mi.5, depth=2, prob=0.95))
	#	cannot make sense out of male_circ < 0, leave out	
	
	mi.3 <- map2stan(
			alist(
					TRANSMITTER ~ dbinom(1, ptr),
					logit(ptr) <- male_in_hh[AGE_AT_MID_D2]*MALE*HH + female_in_hh[AGE_AT_MID_D2]*(1-MALE)*HH +
							male_extra_hh[AGE_AT_MID_D2]*MALE*(1-HH) + female_extra_hh[AGE_AT_MID_D2]*(1-MALE)*(1-HH) +
							vlhigh*VL_HIGH +							
							male_noedu*MALE*NOEDU + female_noedu*(1-MALE)*NOEDU +
							male_inmigrant*MALE*INMIGRANT + female_inmigrant*(1-MALE)*INMIGRANT +
							male_sexp*MALE*SEXP1YR_G1 + female_sexp*(1-MALE)*SEXP1YR_G1,															
					c(vlhigh) ~ dnorm(0,10),					 
					c(male_noedu, female_noedu, male_inmigrant, female_inmigrant, male_sexp, female_sexp) ~ dnorm(0,10),
					female_in_hh[AGE_AT_MID_D2] ~ dnorm(0, sig_female_in_hh),
					male_in_hh[AGE_AT_MID_D2] ~ dnorm(0, sig_male_in_hh),
					female_extra_hh[AGE_AT_MID_D2] ~ dnorm(0, sig_female_extra_hh),
					male_extra_hh[AGE_AT_MID_D2] ~ dnorm(0, sig_male_extra_hh),
					c(sig_female_in_hh, sig_male_in_hh, sig_female_extra_hh, sig_male_extra_hh) ~ dexp(1)
			),
			data=as.data.frame(dg), 
			start=list(	female_in_hh=rep(0,5), male_in_hh=rep(0,5), male_extra_hh=rep(0,5), female_extra_hh=rep(0,5), 
						sig_female_in_hh=1, sig_male_in_hh=1, sig_female_extra_hh=1, sig_male_extra_hh=1,
						vlhigh=0, male_noedu=0, female_noedu=0, male_inmigrant=0, female_inmigrant=0, male_sexp=0, female_sexp=0),			
			warmup=1e3, iter=1e4, chains=1, cores=4
		)
	plot(mi.3)
	precis(mi.3, depth=2, prob=0.95)
	plot(precis(mi.3, depth=2, prob=0.95))
	
	#	extract stuff from mi.2
	post	<- as.data.table(extract.samples(mi.2))
	post[, MC:= seq_len(nrow(post))]
	post	<- melt(post, id.vars='MC')
	tmp		<- post[, which(!grepl('sig',variable))]
	set(post, tmp, 'value', post[tmp, exp(value)])	
	tmp		<- post[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
								V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), 
					  by=c('variable')]
	dmi2	<- dcast.data.table(tmp, variable~STAT, value.var='V')
	dmi2[, MODEL:= 'no age']
	dmi2[, STAT:= as.character(factor(grepl('extra_hh|in_hh',variable), levels=c(TRUE,FALSE), labels=c('odds','odds ratio')))]
	dmi2	<- dmi2[order(STAT, variable),]
	dmi2[, FACTOR:= seq_len(nrow(dmi2))]
	set(dmi2, NULL, 'FACTOR', dmi2[, factor(FACTOR, levels=FACTOR, labels=variable)])
	#	extract stuff from mi.4
	post	<- as.data.table(extract.samples(mi.4))
	post[, MC:= seq_len(nrow(post))]
	post	<- melt(post, id.vars='MC')
	tmp		<- post[, which(!grepl('sig',variable))]
	set(post, tmp, 'value', post[tmp, exp(value)])	
	tmp		<- post[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), 
			by=c('variable')]
	dmi4	<- dcast.data.table(tmp, variable~STAT, value.var='V')
	dmi4[, MODEL:= 'no age']
	dmi4[, STAT:= as.character(factor(grepl('extra_hh|in_hh',variable), levels=c(TRUE,FALSE), labels=c('odds','odds ratio')))]
	dmi4	<- dmi4[order(STAT, variable),]
	dmi4[, FACTOR:= seq_len(nrow(dmi4))]
	set(dmi4, NULL, 'FACTOR', dmi4[, factor(FACTOR, levels=FACTOR, labels=variable)])
	
	
	#	extract stuff from mi.3
	post	<- as.data.table(extract.samples(mi.3))
	post[, MC:= seq_len(nrow(post))]
	post	<- melt(post, id.vars='MC')
	post[, AGE_AT_MID_D2:= as.integer(gsub('^.*V([0-9])$','\\1',variable))]
	set(post, NULL, 'variable', post[, gsub('\\.V[0-9]','',variable)])
	post	<- merge(post, unique(subset(dg, select=c(AGE_AT_MID_D,AGE_AT_MID_D2))), by='AGE_AT_MID_D2', all.x=TRUE)
	tmp		<- post[, which(!grepl('sig',variable))]
	set(post, tmp, 'value', post[tmp, exp(value)])
	set(post, post[, which(is.na(AGE_AT_MID_D))], 'AGE_AT_MID_D', '')
	tmp		<- post[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
								V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), 
					  by=c('variable','AGE_AT_MID_D')]
	dmi3	<- dcast.data.table(tmp, variable+AGE_AT_MID_D~STAT, value.var='V')
	dmi3[, MODEL:= 'with age']
	dmi3[, STAT:= as.character(factor(grepl('extra_hh|in_hh',variable), levels=c(TRUE,FALSE), labels=c('odds','odds ratio')))]
	set(dmi3, dmi3[, which(grepl('sig_',variable))], 'STAT', 'xhyperparam')
	dmi3	<- dmi3[order(STAT, variable, AGE_AT_MID_D),]
	dmi3[, FACTOR:= seq_len(nrow(dmi3))]
	set(dmi3, NULL, 'FACTOR', dmi3[, factor(FACTOR, levels=FACTOR, labels=paste0(variable,'_',AGE_AT_MID_D))])
	
	#	save
	save(dg, mi.2, mi.3, mi.4, mi.5, dmi2, dmi3, dmi4, file=paste0(outfile.base,'_interaction_runs.rda'))	
	
	#	make plots
	dp		<- copy(dmi2)	
	dp[, FILL:= '0']
	set(dp, dp[, which(IU<1)], 'FILL', '-1')
	set(dp, dp[, which(CU<1)], 'FILL', '-2')	
	set(dp, dp[, which(IL>1)], 'FILL', '1')
	set(dp, dp[, which(CL>1)], 'FILL', '2')
	ggplot(dp, aes(x=FACTOR)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/6, 1/3,1/2,1,2,3, 6), labels=c('1/6','1/3','1/2','1','2','3','6'), expand=c(0,0)) +
			scale_fill_manual(values=c('-2'='deepskyblue','-1'='cyan','0'='grey80', '1'='darkorange', '2'='red')) +
			coord_flip(ylim=c(1/5,5)) +
			guides(fill=FALSE) +
			labs(x='Individual characteristics\n', y='\nOdds for being a transmitter vs. recipient\nin phylogenetically reconstructed transmission pairs') +
			facet_grid(STAT~., scales='free', space='free')	
	ggsave(file=paste0(outfile.base,'_interaction_runs_odds_noage.pdf'), w=5, h=3.5)
	
	dp		<- copy(dmi4)	
	dp[, FILL:= '0']
	set(dp, dp[, which(IU<1)], 'FILL', '-1')
	set(dp, dp[, which(CU<1)], 'FILL', '-2')	
	set(dp, dp[, which(IL>1)], 'FILL', '1')
	set(dp, dp[, which(CL>1)], 'FILL', '2')
	ggplot(dp, aes(x=FACTOR)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/6, 1/3,1/2,1,2,3, 6), labels=c('1/6','1/3','1/2','1','2','3','6'), expand=c(0,0)) +
			scale_fill_manual(values=c('-2'='deepskyblue','-1'='cyan','0'='grey80', '1'='darkorange', '2'='red')) +
			coord_flip(ylim=c(1/5,5)) +
			guides(fill=FALSE) +
			labs(x='Individual characteristics\n', y='\nOdds for being a transmitter vs. recipient\nin phylogenetically reconstructed transmission pairs') +
			facet_grid(STAT~., scales='free', space='free')	
	ggsave(file=paste0(outfile.base,'_interaction_runs_odds_noagewithocc.pdf'), w=5, h=5)
	
	
	dp		<- subset(dmi3, !grepl('sig', FACTOR))	
	dp[, FILL:= '0']
	set(dp, dp[, which(IU<1)], 'FILL', '-1')
	set(dp, dp[, which(CU<1)], 'FILL', '-2')	
	set(dp, dp[, which(IL>1)], 'FILL', '1')
	set(dp, dp[, which(CL>1)], 'FILL', '2')	
	ggplot(dp, aes(x=FACTOR)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/20, 1/10, 1/6, 1/3,1/2,1,2,3, 6,10,20), labels=c('1/20','1/10','1/6','1/3','1/2','1','2','3','6','10', '20'), expand=c(0,0)) +
			scale_fill_manual(values=c('-2'='deepskyblue','-1'='cyan','0'='grey80', '1'='darkorange', '2'='red')) +
			coord_flip(ylim=c(1/20,20)) +
			guides(fill=FALSE) +
			labs(x='Individual characteristics\n', y='\nOdds for being a transmitter vs. recipient\nin phylogenetically reconstructed transmission pairs') +
			facet_grid(STAT~., scales='free', space='free')
	ggsave(file=paste0(outfile.base,'_interaction_runs_odds_withage.pdf'), w=5, h=8)
	
	dp		<- subset(dmi3, !grepl('sig', FACTOR) & grepl('in_hh', FACTOR))	
	dp[, FILL:= '0']
	set(dp, dp[, which(IU<1)], 'FILL', '-1')
	set(dp, dp[, which(CU<1)], 'FILL', '-2')	
	set(dp, dp[, which(IL>1)], 'FILL', '1')
	set(dp, dp[, which(CL>1)], 'FILL', '2')
	dp[, SEX:= gsub('([a-z]+)_([a-z]+_[a-z]+)','\\1',variable)]
	ggplot(dp, aes(x=AGE_AT_MID_D)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/20, 1/10, 1/6, 1/3,1/2,1,2,3, 6,10,20), labels=c('1/20','1/10','1/6','1/3','1/2','1','2','3','6','10', '20'), expand=c(0,0)) +
			scale_fill_manual(values=c('-2'='deepskyblue','-1'='cyan','0'='grey80', '1'='darkorange', '2'='red')) +
			coord_flip(ylim=c(1/10,10)) +
			guides(fill=FALSE) +
			labs(x='Individual characteristics\n', y='\nOdds for being a transmitter vs. recipient\nin phylogenetically reconstructed transmission pairs') +
			facet_grid(STAT+variable~., scales='free', space='free')
	ggsave(file=paste0(outfile.base,'_interaction_runs_odds_onluage.pdf'), w=5, h=4)
}

RakaiFull.transmitter.171122.gender.multivariatemodels.stan.with.threshold<- function()
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
	
	infile			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_transmitterrecipientdata.rda"
	outfile.base	<- gsub('data.rda','',infile)
	load(infile)
	#
	#	prepare variables for STAN
	#		
	dg	<- copy(df)
	#	prepare data for STAN
	dg[, COMM_FISH:= as.integer(COMM_TYPE=='fisherfolk')]
	dg[, COMM_AGR:= as.integer(COMM_TYPE=='agrarian')]
	dg[, COMM_TRAD:= as.integer(COMM_TYPE=='trading')]
	#dg[, COMM_MXD:= as.integer(PAIR_COMM_TYPE=='mixed')]
	dg[, HH:= as.integer(SAMEHH=='same hh')]
	dg[, MALE:= as.integer(SEX=='M')]
	dg[, NOEDU:= as.integer(EDUCAT=='None')]
	dg[, NOEDU_MISS:= as.integer(EDUCAT=='Unknown')]
	dg[, SEXP1YR_G1:= as.integer(SEXP1YR!='1')]
	dg[, VL_HIGH:= as.integer(RECENTVL>=1e5)]
	#dg[, VL_SUPP:= as.integer(RECENTVL<=200)]
	dg[, OC_AGRO:= as.integer(OCAT=='Agro/House')]
	dg[, OC_BAR:= as.integer(OCAT=='Bar/waitress')]
	dg[, OC_BODA:= as.integer(OCAT=='Boda/Trucking')]
	dg[, OC_FISH:= as.integer(OCAT=='Fishing')]
	dg[, OC_STUD:= as.integer(OCAT=='Student')]
	dg[, OC_TRAD:= as.integer(OCAT=='Trading/Shop keeper')]
	dg[, OC_OTH:= as.integer(OCAT%in%c('zother'))]	
	dg[, AGE_AT_MID_C2:= dg[,as.integer(cut(AGE_AT_MID, breaks=c(15,20,25,30,35,40,45,52), right=FALSE, levels=c('1','2','3','4','5','6','7')))]]
	#
	#	look only at transmitters among males
	#
	df	<- subset(dg, MALE==1)		
	#	model 2 with age (exchangeable)
	mma.2 <- map2stan(
			alist(
					TRANSMITTER ~ dbinom(1, ptr),
					logit(ptr) <- base + age[AGE_AT_MID_C2] +   
							# agrarian community is base 
							comm_fish*COMM_FISH + comm_trad*COMM_TRAD + 							
							noedu*NOEDU + sexp*SEXP1YR_G1 +	vlhigh*VL_HIGH + hh*HH + inmigrant*INMIGRANT +							
							# other occupation is baseline
							occ_agro*OC_AGRO + occ_fish*OC_FISH + occ_boda*OC_BODA + occ_stud*OC_STUD + occ_trad*OC_TRAD,
					base ~ dnorm(0,100),
					age[AGE_AT_MID_C2] ~ dnorm(0, age_sig),					
					c(comm_fish, comm_trad) ~ dnorm(0,10),
					c(noedu, sexp, vlhigh, hh, inmigrant) ~ dnorm(0,10),								
					c(occ_agro, occ_boda, occ_stud, occ_trad, occ_fish) ~ dnorm(0,10),
					c(age_sig) ~ dexp(1)
			),
			data=as.data.frame(df), 
			start=list(	base=0, comm_fish=0, comm_trad=0, age=rep(0,7), noedu=0, sexp=0, vlhigh=0, hh=0, inmigrant=0,						
						occ_agro=0, occ_boda=0, occ_stud=0, occ_trad=0, occ_fish=0,
						age_sig=1),			
			warmup=1e3, iter=1e4, chains=1, cores=4
		)
	#plot(precis(mhh.2, prob=0.95, depth=2))
	#pairs(mhh.1)
	post	<- extract.samples(mma.2)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_COMM_AGRO=logistic( post$base ),				
			PI_COMM_FISH=logistic( post$base+post$comm_fish ), 											
			PI_COMM_TRAD=logistic( post$base+post$comm_trad ),										
			PI_EDU_NONE=logistic( post$base+post$noedu ),										
			PI_SEXP1YR_G1=logistic( post$base+post$sexp ),
			PI_VL_HIGH=logistic( post$base+post$vlhigh ),			
			PI_MOBILITY_INMIGRANT=logistic( post$base+post$inmigrant ),
			PI_HH_SAME=logistic( post$base+post$hh ),			
			PI_OCC_AGRO=logistic( post$base+post$occ_agro ),
			#PI_OCC_BAR=logistic( post$base+post$occ_bar ),
			PI_OCC_BODA=logistic( post$base+post$occ_boda ),
			PI_OCC_FISH=logistic( post$base+post$occ_fish ),
			PI_OCC_STUD=logistic( post$base+post$occ_stud ),
			PI_OCC_TRAD=logistic( post$base+post$occ_trad ),
			PI_AGE_1=logistic( post$base+post$age[,1] ),
			PI_AGE_2=logistic( post$base+post$age[,2] ),
			PI_AGE_3=logistic( post$base+post$age[,3] ),
			PI_AGE_4=logistic( post$base+post$age[,4] ),
			PI_AGE_5=logistic( post$base+post$age[,5] ),
			PI_AGE_6=logistic( post$base+post$age[,6] ),
			PI_AGE_7=logistic( post$base+post$age[,7] ),
			OR_COMM_AGRO=exp( post$base ),				
			OR_COMM_FISH=exp( post$base+post$comm_fish ), 											
			OR_COMM_TRAD=exp( post$base+post$comm_trad ),										
			OR_EDU_NONE=exp( post$base+post$noedu ),										
			OR_SEXP1YR_G1=exp( post$base+post$sexp ),
			OR_VL_HIGH=exp( post$base+post$vlhigh ),
			OR_MOBILITY_INMIGRANT=exp( post$base+post$inmigrant ),
			OR_HH_SAME=exp( post$base+post$hh ),			
			OR_OCC_AGRO=exp( post$base+post$occ_agro ),
			#OR_OCC_BAR=exp( post$base+post$occ_bar ),
			OR_OCC_BODA=exp( post$base+post$occ_boda ),
			OR_OCC_FISH=exp( post$base+post$occ_fish ),
			OR_OCC_STUD=exp( post$base+post$occ_stud ),
			OR_OCC_TRAD=exp( post$base+post$occ_trad ),
			OR_AGE_1=exp( post$base+post$age[,1] ),
			OR_AGE_2=exp( post$base+post$age[,2] ),
			OR_AGE_3=exp( post$base+post$age[,3] ),
			OR_AGE_4=exp( post$base+post$age[,4] ),
			OR_AGE_5=exp( post$base+post$age[,5] ),
			OR_AGE_6=exp( post$base+post$age[,6] ),
			OR_AGE_7=exp( post$base+post$age[,7] ),
			ORX_COMM_FISH=exp( post$comm_fish ), 											
			ORX_COMM_TRAD=exp( post$comm_trad ),										
			ORX_EDU_NONE=exp( post$noedu ),										
			ORX_SEXP1YR_G1=exp( post$sexp ),
			ORX_VL_HIGH=exp( post$vlhigh ),
			ORX_MOBILITY_INMIGRANT=exp( post$inmigrant ),
			ORX_HH_SAME=exp( post$hh ),						
			ORX_OCC_AGRO=exp( post$occ_agro ),
			#ORX_OCC_BAR=exp( post$occ_bar ),
			ORX_OCC_BODA=exp( post$occ_boda ),
			ORX_OCC_FISH=exp( post$occ_fish ),
			ORX_OCC_STUD=exp( post$occ_stud ),
			ORX_OCC_TRAD=exp( post$occ_trad ),
			ORX_AGE_1=exp( post$age[,1] ),
			ORX_AGE_2=exp( post$age[,2] ),
			ORX_AGE_3=exp( post$age[,3] ),
			ORX_AGE_4=exp( post$age[,4] ),
			ORX_AGE_5=exp( post$age[,5] ),
			ORX_AGE_6=exp( post$age[,6] ),
			ORX_AGE_7=exp( post$age[,7] )
			)
	dp		<- melt(dp, id.vars='MC')	
	dma.2	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
							V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	dma.2	<- dcast.data.table(dma.2, variable~STAT, value.var='V')
	dma.2[, MODEL:= 'multivariate male with age']
	#
	#	look only at transmitters outside households
	#
	df	<- subset(dg, MALE==0)		
	mfe.2 <- map2stan(
			alist(
					TRANSMITTER ~ dbinom(1, ptr),
					logit(ptr) <- base + age[AGE_AT_MID_C2] +   
							# agrarian community is base 
							comm_fish*COMM_FISH + comm_trad*COMM_TRAD + 							
							noedu*NOEDU + sexp*SEXP1YR_G1 +	vlhigh*VL_HIGH + hh*HH + inmigrant*INMIGRANT +							
							# other occupation is baseline
							occ_agro*OC_AGRO + occ_bar*OC_BAR + occ_stud*OC_STUD + occ_trad*OC_TRAD,
					base ~ dnorm(0,100),
					age[AGE_AT_MID_C2] ~ dnorm(0, age_sig),					
					c(comm_fish, comm_trad) ~ dnorm(0,10),
					c(noedu, sexp, vlhigh, hh, inmigrant) ~ dnorm(0,10),								
					c(occ_agro, occ_bar, occ_stud, occ_trad) ~ dnorm(0,10),
					c(age_sig) ~ dexp(1)
			),
			data=as.data.frame(df), 
			start=list(	base=0, comm_fish=0, comm_trad=0, age=rep(0,7), noedu=0, sexp=0, vlhigh=0, hh=0, inmigrant=0,						
					occ_agro=0, occ_bar=0, occ_stud=0, occ_trad=0,
					age_sig=1),			
			warmup=1e3, iter=1e4, chains=1, cores=4
		)
	#plot(precis(mfe.2, prob=0.95, depth=2))
	#plot(mfe.2)
	post	<- extract.samples(mfe.2)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_COMM_AGRO=logistic( post$base ),				
			PI_COMM_FISH=logistic( post$base+post$comm_fish ), 											
			PI_COMM_TRAD=logistic( post$base+post$comm_trad ),										
			PI_EDU_NONE=logistic( post$base+post$noedu ),										
			PI_SEXP1YR_G1=logistic( post$base+post$sexp ),
			PI_VL_HIGH=logistic( post$base+post$vlhigh ),			
			PI_MOBILITY_INMIGRANT=logistic( post$base+post$inmigrant ),
			PI_HH_SAME=logistic( post$base+post$hh ),			
			PI_OCC_AGRO=logistic( post$base+post$occ_agro ),
			PI_OCC_BAR=logistic( post$base+post$occ_bar ),
			PI_OCC_STUD=logistic( post$base+post$occ_stud ),
			PI_OCC_TRAD=logistic( post$base+post$occ_trad ),
			PI_AGE_1=logistic( post$base+post$age[,1] ),
			PI_AGE_2=logistic( post$base+post$age[,2] ),
			PI_AGE_3=logistic( post$base+post$age[,3] ),
			PI_AGE_4=logistic( post$base+post$age[,4] ),
			PI_AGE_5=logistic( post$base+post$age[,5] ),
			PI_AGE_6=logistic( post$base+post$age[,6] ),
			PI_AGE_7=logistic( post$base+post$age[,7] ),
			OR_COMM_AGRO=exp( post$base ),				
			OR_COMM_FISH=exp( post$base+post$comm_fish ), 											
			OR_COMM_TRAD=exp( post$base+post$comm_trad ),										
			OR_EDU_NONE=exp( post$base+post$noedu ),										
			OR_SEXP1YR_G1=exp( post$base+post$sexp ),
			OR_VL_HIGH=exp( post$base+post$vlhigh ),
			OR_MOBILITY_INMIGRANT=exp( post$base+post$inmigrant ),
			OR_HH_SAME=exp( post$base+post$hh ),			
			OR_OCC_AGRO=exp( post$base+post$occ_agro ),
			OR_OCC_BAR=exp( post$base+post$occ_bar ),
			OR_OCC_STUD=exp( post$base+post$occ_stud ),
			OR_OCC_TRAD=exp( post$base+post$occ_trad ),
			OR_AGE_1=exp( post$base+post$age[,1] ),
			OR_AGE_2=exp( post$base+post$age[,2] ),
			OR_AGE_3=exp( post$base+post$age[,3] ),
			OR_AGE_4=exp( post$base+post$age[,4] ),
			OR_AGE_5=exp( post$base+post$age[,5] ),
			OR_AGE_6=exp( post$base+post$age[,6] ),
			OR_AGE_7=exp( post$base+post$age[,7] ),
			ORX_COMM_FISH=exp( post$comm_fish ), 											
			ORX_COMM_TRAD=exp( post$comm_trad ),										
			ORX_EDU_NONE=exp( post$noedu ),										
			ORX_SEXP1YR_G1=exp( post$sexp ),
			ORX_VL_HIGH=exp( post$vlhigh ),
			ORX_MOBILITY_INMIGRANT=exp( post$inmigrant ),
			ORX_HH_SAME=exp( post$hh ),						
			ORX_OCC_AGRO=exp( post$occ_agro ),
			ORX_OCC_BAR=exp( post$occ_bar ),
			ORX_OCC_STUD=exp( post$occ_stud ),
			ORX_OCC_TRAD=exp( post$occ_trad ),
			ORX_AGE_1=exp( post$age[,1] ),
			ORX_AGE_2=exp( post$age[,2] ),
			ORX_AGE_3=exp( post$age[,3] ),
			ORX_AGE_4=exp( post$age[,4] ),
			ORX_AGE_5=exp( post$age[,5] ),
			ORX_AGE_6=exp( post$age[,6] ),
			ORX_AGE_7=exp( post$age[,7] )
			)
	dp		<- melt(dp, id.vars='MC')	
	dfe.2	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	dfe.2	<- dcast.data.table(dfe.2, variable~STAT, value.var='V')
	dfe.2[, MODEL:= 'multivariate female with age']
	
	
	ds		<- rbind(dfe.2, dma.2)
	ds[, LABEL:= paste0(round(MED, d=2), '\n[', round(CL, d=2),'-', round( CU, d=2),']')]
	ds[, STAT:= factor(gsub('^([^_]+)_.*','\\1',variable), levels=c('PI','OR','ORX'), labels=c('probability transmitter','odds transmitter','odds ratio'))]
	ds[, FACTOR:=gsub('^([^_]+)_([^_]+)_([^_]+)','\\3',variable)]
	ds[, MOFA:=gsub('^([^_]+)_([^_]+)_([^_]+)','\\2-\\3',variable)]		
	#ds		<- subset(ds, STAT=='odds ratio')
	set(ds, NULL, 'MOFA2', ds[, factor(MOFA, levels=rev(c(
											"AGE-1","AGE-2","AGE-3","AGE-4","AGE-5","AGE-6","AGE-7",											
											"COMM-AGRO","COMM-FISH","COMM-TRAD","EDU-NONE","SEXP1YR-G1","VL-HIGH","MOBILITY-INMIGRANT","HH-SAME",
											"OCC-AGRO","OCC-BAR","OCC-BODA","OCC-FISH","OCC-STUD","OCC-TRAD"											
									)), labels=rev(c("age 15-19","age 20-24","age 25-29","age 30-34","age 35-39","age 40-44","age 45-50",											
											"Agrarian community", "Fishing site", "Trading centre",
											"No primary education","More than 1 sex partner in last year", "Viral load above 100,000 cps/ml","Inmigrant","Same household",											
											"Primary occupation: agricultural/house", "Primary occupation: bar worker/waitress","Primary occupation: Boda/trucking",											
											"Primary occupation: fishing", "Primary occupation: student", "Primary occupation: trading/shopkeeper"
									)))])
	setkey(ds, MOFA)	
	#
	#	save
	#
	save(dg, ds, dfe.2, dma.2, mma.2, mfe.2, file=paste0(outfile.base,'_gender_runs.rda'))
	
	dp		<- subset(ds, STAT=='odds transmitter')
	dp[, FILL:= '0']
	set(dp, dp[, which(IL>1 | IU<1)], 'FILL', '1')
	set(dp, dp[, which(CL>1 | CU<1)], 'FILL', '2')
	ggplot(subset(dp, MODEL=='multivariate hh with age'), aes(x=MOFA2)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/20, 1/10, 1/6, 1/3,1/2,1,2,3, 6,10,20), labels=c('1/20','1/10','1/6','1/3','1/2','1','2','3','6','10', '20'), expand=c(0,0)) +
			scale_fill_manual(values=c('0'='white', '1'='orange', '2'='red')) +
			coord_flip(ylim=c(1/20,20)) +
			guides(fill=FALSE) +
			labs(x='Individual characteristics\n', y='\nOdds for being a transmitter vs. recipient\nin phylogenetically reconstructed transmission pairs')
	ggsave(file=paste0(outfile.base,'_odds_samehh_withage.pdf'), w=6, h=8)
	ggplot(subset(dp, MODEL=='multivariate outside hh with age'), aes(x=MOFA2)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/20, 1/10, 1/6, 1/3,1/2,1,2,3, 6,10,20), labels=c('1/20','1/10','1/6','1/3','1/2','1','2','3','6','10', '20'), expand=c(0,0)) +
			scale_fill_manual(values=c('0'='white', '1'='orange', '2'='red')) +
			coord_flip(ylim=c(1/20,20)) +
			guides(fill=FALSE) +
			labs(x='Individual characteristics\n', y='\nOdds for being a transmitter vs. recipient\nin phylogenetically reconstructed transmission pairs')
	ggsave(file=paste0(outfile.base,'_odds_diffhh_withage.pdf'), w=6, h=8)
	
	dp		<- subset(ds, STAT=='odds ratio')
	dp[, FILL:= '0']
	set(dp, dp[, which(IL>1 | IU<1)], 'FILL', '1')
	set(dp, dp[, which(CL>1 | CU<1)], 'FILL', '2')	
	ggplot(subset(dp, MODEL=='multivariate hh no age'), aes(x=MOFA2)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/20, 1/10, 1/6, 1/3,1/2,1,2,3, 6,10,20), labels=c('1/20','1/10','1/6','1/3','1/2','1','2','3','6','10', '20'), expand=c(0,0)) +
			scale_fill_manual(values=c('0'='white', '1'='orange', '2'='red')) +
			coord_flip(ylim=c(1/20,20)) +
			guides(fill=FALSE) +
			labs(x='Individual characteristics\n', y='\nOdds ratio for being a transmitter vs. recipient\nin phylogenetically reconstructed transmission pairs')
	ggsave(file=paste0(outfile.base,'_odds_samehh_noage.pdf'), w=6, h=5)
	ggplot(subset(dp, MODEL=='multivariate outside hh no age'), aes(x=MOFA2)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/20, 1/10, 1/6, 1/3,1/2,1,2,3, 6,10,20), labels=c('1/20','1/10','1/6','1/3','1/2','1','2','3','6','10', '20'), expand=c(0,0)) +
			scale_fill_manual(values=c('0'='white', '1'='orange', '2'='red')) +
			coord_flip(ylim=c(1/20,20)) +
			guides(fill=FALSE) +
			labs(x='Individual characteristics\n', y='\nOdds ratio for being a transmitter vs. recipient\nin phylogenetically reconstructed transmission pairs')
	ggsave(file=paste0(outfile.base,'_odds_diffhh_noage.pdf'), w=6, h=5)
}

RakaiFull.transmitter.171122.hh.multivariatemodels.stan.with.threshold<- function()
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
	
	infile			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_transmitterrecipientdata.rda"
	outfile.base	<- gsub('data.rda','',infile)
	load(infile)
	#
	#	prepare variables for STAN
	#		
	dg	<- copy(df)
	#	prepare data for STAN
	dg[, COMM_FISH:= as.integer(COMM_TYPE=='fisherfolk')]
	dg[, COMM_AGR:= as.integer(COMM_TYPE=='agrarian')]
	dg[, COMM_TRAD:= as.integer(COMM_TYPE=='trading')]
	#dg[, COMM_MXD:= as.integer(PAIR_COMM_TYPE=='mixed')]
	dg[, MALE:= as.integer(SEX=='M')]
	dg[, NOEDU:= as.integer(EDUCAT=='None')]
	dg[, NOEDU_MISS:= as.integer(EDUCAT=='Unknown')]
	dg[, SEXP1YR_G1:= as.integer(SEXP1YR!='1')]
	dg[, VL_HIGH:= as.integer(RECENTVL>=1e5)]
	#dg[, VL_SUPP:= as.integer(RECENTVL<=200)]
	dg[, OC_AGRO:= as.integer(OCAT=='Agro/House')]
	dg[, OC_BAR:= as.integer(OCAT=='Bar/waitress')]
	dg[, OC_BODA:= as.integer(OCAT=='Boda/Trucking')]
	dg[, OC_FISH:= as.integer(OCAT=='Fishing')]
	dg[, OC_STUD:= as.integer(OCAT=='Student')]
	dg[, OC_TRAD:= as.integer(OCAT=='Trading/Shop keeper')]
	dg[, OC_OTH:= as.integer(OCAT%in%c('zother'))]	
	dg[, AGE_AT_MID_C2:= dg[,as.integer(cut(AGE_AT_MID, breaks=c(15,20,25,30,35,40,45,52), right=FALSE, levels=c('1','2','3','4','5','6','7')))]]
	#
	#	look only at transmitters within households
	#
	df	<- subset(dg, SAMEHH=='same hh')		
	#	model 1 no age	
	mhh.1 <- map2stan(
			alist(
					TRANSMITTER ~ dbinom(1, ptr),
					logit(ptr) <- base + 
							# agrarian community is base 
							comm_fish*COMM_FISH + comm_trad*COMM_TRAD + 
							male*MALE + noedu*NOEDU + sexp*SEXP1YR_G1 +	vlhigh*VL_HIGH + inmigrant*INMIGRANT +						
							# other occupation is baseline
							occ_agro*OC_AGRO + occ_fish*OC_FISH + occ_bar*OC_BAR + occ_boda*OC_BODA + occ_stud*OC_STUD + occ_trad*OC_TRAD,
					base ~ dnorm(0,100),								
					c(comm_fish, comm_trad) ~ dnorm(0,10),
					c(male, noedu, sexp, vlhigh, inmigrant) ~ dnorm(0,10),								
					c(occ_agro, occ_bar, occ_boda, occ_stud, occ_trad, occ_fish) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, comm_fish=0, comm_trad=0, male=0, inmigrant=0, noedu=0, sexp=0, vlhigh=0, 						
					occ_agro=0, occ_bar=0, occ_boda=0, occ_stud=0, occ_trad=0, occ_fish=0),			
			warmup=1e3, iter=5e3, chains=1, cores=4
	)
	#plot(precis(mhh.1, prob=0.95))
	#pairs(mhh.1)
	post	<- extract.samples(mhh.1)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_GENDER_MALE=logistic( post$base+post$male ),
			PI_COMM_AGRO=logistic( post$base ),				
			PI_COMM_FISH=logistic( post$base+post$comm_fish ), 											
			PI_COMM_TRAD=logistic( post$base+post$comm_trad ),										
			PI_EDU_NONE=logistic( post$base+post$noedu ),										
			PI_SEXP1YR_G1=logistic( post$base+post$sexp ),
			PI_VL_HIGH=logistic( post$base+post$vlhigh ),
			PI_MOBILITY_INMIGRANT=logistic( post$base+post$inmigrant ),
			PI_OCC_AGRO=logistic( post$base+post$occ_agro ),
			PI_OCC_BAR=logistic( post$base+post$occ_bar ),
			PI_OCC_BODA=logistic( post$base+post$occ_boda ),
			PI_OCC_FISH=logistic( post$base+post$occ_fish ),
			PI_OCC_STUD=logistic( post$base+post$occ_stud ),
			PI_OCC_TRAD=logistic( post$base+post$occ_trad ),
			OR_GENDER_MALE=exp( post$base+post$male ),
			OR_COMM_AGRO=exp( post$base ),				
			OR_COMM_FISH=exp( post$base+post$comm_fish ), 											
			OR_COMM_TRAD=exp( post$base+post$comm_trad ),										
			OR_EDU_NONE=exp( post$base+post$noedu ),										
			OR_SEXP1YR_G1=exp( post$base+post$sexp ),
			OR_VL_HIGH=exp( post$base+post$vlhigh ),	
			OR_MOBILITY_INMIGRANT=exp( post$base+post$inmigrant ),
			OR_OCC_AGRO=exp( post$base+post$occ_agro ),
			OR_OCC_BAR=exp( post$base+post$occ_bar ),
			OR_OCC_BODA=exp( post$base+post$occ_boda ),
			OR_OCC_FISH=exp( post$base+post$occ_fish ),
			OR_OCC_STUD=exp( post$base+post$occ_stud ),
			OR_OCC_TRAD=exp( post$base+post$occ_trad ),
			ORX_GENDER_MALE=exp( post$male ),
			ORX_COMM_FISH=exp( post$comm_fish ), 											
			ORX_COMM_TRAD=exp( post$comm_trad ),										
			ORX_EDU_NONE=exp( post$noedu ),										
			ORX_SEXP1YR_G1=exp( post$sexp ),
			ORX_VL_HIGH=exp( post$vlhigh ),	
			ORX_MOBILITY_INMIGRANT=exp( post$inmigrant ),
			ORX_OCC_AGRO=exp( post$occ_agro ),
			ORX_OCC_BAR=exp( post$occ_bar ),
			ORX_OCC_BODA=exp( post$occ_boda ),
			ORX_OCC_FISH=exp( post$occ_fish ),
			ORX_OCC_STUD=exp( post$occ_stud ),
			ORX_OCC_TRAD=exp( post$occ_trad )			
	)
	dp		<- melt(dp, id.vars='MC')	
	dhh.1	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	dhh.1	<- dcast.data.table(dhh.1, variable~STAT, value.var='V')
	dhh.1[, MODEL:= 'multivariate hh no age']
	
	#	model 2 with age (exchangeable)
	mhh.2 <- map2stan(
			alist(
					TRANSMITTER ~ dbinom(1, ptr),
					logit(ptr) <- base + 
							# agrarian community is base 
							comm_fish*COMM_FISH + comm_trad*COMM_TRAD + 
							male_age[AGE_AT_MID_C2]*MALE + female_age[AGE_AT_MID_C2]*(1-MALE) + 
							noedu*NOEDU + sexp*SEXP1YR_G1 +	vlhigh*VL_HIGH + inmigrant*INMIGRANT +							
							# other occupation is baseline
							occ_agro*OC_AGRO + occ_fish*OC_FISH + occ_bar*OC_BAR + occ_boda*OC_BODA + occ_stud*OC_STUD + occ_trad*OC_TRAD,
					base ~ dnorm(0,100),
					male_age[AGE_AT_MID_C2] ~ dnorm(0, male_age_sig),
					female_age[AGE_AT_MID_C2] ~ dnorm(0, female_age_sig),
					c(comm_fish, comm_trad) ~ dnorm(0,10),
					c(noedu, sexp, vlhigh, inmigrant) ~ dnorm(0,10),								
					c(occ_agro, occ_bar, occ_boda, occ_stud, occ_trad, occ_fish) ~ dnorm(0,10),
					c(male_age_sig, female_age_sig) ~ dcauchy(0,1)
			),
			data=as.data.frame(df), 
			start=list(	base=0, comm_fish=0, comm_trad=0, male_age=rep(0,7), female_age=rep(0,7), noedu=0, inmigrant=0, sexp=0, vlhigh=0, 						
					occ_agro=0, occ_bar=0, occ_boda=0, occ_stud=0, occ_trad=0, occ_fish=0,
					male_age_sig=1, female_age_sig=1),			
			warmup=1e3, iter=5e3, chains=1, cores=4
		)
	#plot(precis(mhh.2, prob=0.95))
	#pairs(mhh.1)
	post	<- extract.samples(mhh.2)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_COMM_AGRO=logistic( post$base ),				
			PI_COMM_FISH=logistic( post$base+post$comm_fish ), 											
			PI_COMM_TRAD=logistic( post$base+post$comm_trad ),										
			PI_EDU_NONE=logistic( post$base+post$noedu ),										
			PI_SEXP1YR_G1=logistic( post$base+post$sexp ),
			PI_VL_HIGH=logistic( post$base+post$vlhigh ),
			PI_MOBILITY_INMIGRANT=logistic( post$base+post$inmigrant ),
			PI_OCC_AGRO=logistic( post$base+post$occ_agro ),
			PI_OCC_BAR=logistic( post$base+post$occ_bar ),
			PI_OCC_BODA=logistic( post$base+post$occ_boda ),
			PI_OCC_FISH=logistic( post$base+post$occ_fish ),
			PI_OCC_STUD=logistic( post$base+post$occ_stud ),
			PI_OCC_TRAD=logistic( post$base+post$occ_trad ),
			PI_MALEAGE_1=logistic( post$base+post$male_age[,1] ),
			PI_MALEAGE_2=logistic( post$base+post$male_age[,2] ),
			PI_MALEAGE_3=logistic( post$base+post$male_age[,3] ),
			PI_MALEAGE_4=logistic( post$base+post$male_age[,4] ),
			PI_MALEAGE_5=logistic( post$base+post$male_age[,5] ),
			PI_MALEAGE_6=logistic( post$base+post$male_age[,6] ),
			PI_MALEAGE_7=logistic( post$base+post$male_age[,7] ),
			PI_FEMALEAGE_1=logistic( post$base+post$female_age[,1] ),
			PI_FEMALEAGE_2=logistic( post$base+post$female_age[,2] ),
			PI_FEMALEAGE_3=logistic( post$base+post$female_age[,3] ),
			PI_FEMALEAGE_4=logistic( post$base+post$female_age[,4] ),
			PI_FEMALEAGE_5=logistic( post$base+post$female_age[,5] ),
			PI_FEMALEAGE_6=logistic( post$base+post$female_age[,6] ),
			PI_FEMALEAGE_7=logistic( post$base+post$female_age[,7] ),			
			OR_COMM_AGRO=exp( post$base ),				
			OR_COMM_FISH=exp( post$base+post$comm_fish ), 											
			OR_COMM_TRAD=exp( post$base+post$comm_trad ),										
			OR_EDU_NONE=exp( post$base+post$noedu ),										
			OR_SEXP1YR_G1=exp( post$base+post$sexp ),
			OR_VL_HIGH=exp( post$base+post$vlhigh ),
			OR_MOBILITY_INMIGRANT=exp( post$base+post$inmigrant ),
			OR_OCC_AGRO=exp( post$base+post$occ_agro ),
			OR_OCC_BAR=exp( post$base+post$occ_bar ),
			OR_OCC_BODA=exp( post$base+post$occ_boda ),
			OR_OCC_FISH=exp( post$base+post$occ_fish ),
			OR_OCC_STUD=exp( post$base+post$occ_stud ),
			OR_OCC_TRAD=exp( post$base+post$occ_trad ),
			OR_MALEAGE_1=exp( post$base+post$male_age[,1] ),
			OR_MALEAGE_2=exp( post$base+post$male_age[,2] ),
			OR_MALEAGE_3=exp( post$base+post$male_age[,3] ),
			OR_MALEAGE_4=exp( post$base+post$male_age[,4] ),
			OR_MALEAGE_5=exp( post$base+post$male_age[,5] ),
			OR_MALEAGE_6=exp( post$base+post$male_age[,6] ),
			OR_MALEAGE_7=exp( post$base+post$male_age[,7] ),
			OR_FEMALEAGE_1=exp( post$base+post$female_age[,1] ),
			OR_FEMALEAGE_2=exp( post$base+post$female_age[,2] ),
			OR_FEMALEAGE_3=exp( post$base+post$female_age[,3] ),
			OR_FEMALEAGE_4=exp( post$base+post$female_age[,4] ),
			OR_FEMALEAGE_5=exp( post$base+post$female_age[,5] ),
			OR_FEMALEAGE_6=exp( post$base+post$female_age[,6] ),
			OR_FEMALEAGE_7=exp( post$base+post$female_age[,7] ),							
			ORX_COMM_FISH=exp( post$comm_fish ), 											
			ORX_COMM_TRAD=exp( post$comm_trad ),										
			ORX_EDU_NONE=exp( post$noedu ),										
			ORX_SEXP1YR_G1=exp( post$sexp ),
			ORX_VL_HIGH=exp( post$vlhigh ),
			ORX_MOBILITY_INMIGRANT=exp( post$inmigrant ),
			ORX_OCC_AGRO=exp( post$occ_agro ),
			ORX_OCC_BAR=exp( post$occ_bar ),
			ORX_OCC_BODA=exp( post$occ_boda ),
			ORX_OCC_FISH=exp( post$occ_fish ),
			ORX_OCC_STUD=exp( post$occ_stud ),
			ORX_OCC_TRAD=exp( post$occ_trad ),
			ORX_MALEAGE_1=exp( post$male_age[,1]-post$female_age[,1] ),
			ORX_MALEAGE_2=exp( post$male_age[,2]-post$female_age[,2] ),
			ORX_MALEAGE_3=exp( post$male_age[,3]-post$female_age[,3] ),
			ORX_MALEAGE_4=exp( post$male_age[,4]-post$female_age[,4] ),
			ORX_MALEAGE_5=exp( post$male_age[,5]-post$female_age[,5] ),
			ORX_MALEAGE_6=exp( post$male_age[,6]-post$female_age[,6] ),
			ORX_MALEAGE_7=exp( post$male_age[,7]-post$female_age[,7] )
	)
	dp		<- melt(dp, id.vars='MC')	
	dhh.2	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	dhh.2	<- dcast.data.table(dhh.2, variable~STAT, value.var='V')
	dhh.2[, MODEL:= 'multivariate hh with age']
	#
	#	look only at transmitters outside households
	#
	df	<- subset(dg, SAMEHH!='same hh')		
	#	model 1 no age	
	moh.1 <- map2stan(
			alist(
					TRANSMITTER ~ dbinom(1, ptr),
					logit(ptr) <- base + 
							# agrarian community is base 
							comm_fish*COMM_FISH + comm_trad*COMM_TRAD + 
							male*MALE + noedu*NOEDU + sexp*SEXP1YR_G1 +	vlhigh*VL_HIGH + inmigrant*INMIGRANT +							
							# other occupation is baseline
							occ_agro*OC_AGRO + occ_fish*OC_FISH + occ_bar*OC_BAR + occ_boda*OC_BODA + occ_stud*OC_STUD + occ_trad*OC_TRAD,
					base ~ dnorm(0,100),								
					c(comm_fish, comm_trad) ~ dnorm(0,10),
					c(male, noedu, sexp, vlhigh, inmigrant) ~ dnorm(0,10),								
					c(occ_agro, occ_bar, occ_boda, occ_stud, occ_trad, occ_fish) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, comm_fish=0, comm_trad=0, male=0, noedu=0, sexp=0, vlhigh=0, inmigrant=0,						
					occ_agro=0, occ_bar=0, occ_boda=0, occ_stud=0, occ_trad=0, occ_fish=0),			
			warmup=1e3, iter=5e3, chains=1, cores=4
		)
	#plot(precis(moh.1, prob=0.95))
	#pairs(moh.1)
	post	<- extract.samples(moh.1)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_GENDER_MALE=logistic( post$base+post$male ),
			PI_COMM_AGRO=logistic( post$base ),				
			PI_COMM_FISH=logistic( post$base+post$comm_fish ), 											
			PI_COMM_TRAD=logistic( post$base+post$comm_trad ),										
			PI_EDU_NONE=logistic( post$base+post$noedu ),										
			PI_SEXP1YR_G1=logistic( post$base+post$sexp ),
			PI_VL_HIGH=logistic( post$base+post$vlhigh ),
			PI_MOBILITY_INMIGRANT=logistic( post$base+post$inmigrant ),
			PI_OCC_AGRO=logistic( post$base+post$occ_agro ),
			PI_OCC_BAR=logistic( post$base+post$occ_bar ),
			PI_OCC_BODA=logistic( post$base+post$occ_boda ),
			PI_OCC_FISH=logistic( post$base+post$occ_fish ),
			PI_OCC_STUD=logistic( post$base+post$occ_stud ),
			PI_OCC_TRAD=logistic( post$base+post$occ_trad ),
			OR_GENDER_MALE=exp( post$base+post$male ),
			OR_COMM_AGRO=exp( post$base ),				
			OR_COMM_FISH=exp( post$base+post$comm_fish ), 											
			OR_COMM_TRAD=exp( post$base+post$comm_trad ),										
			OR_EDU_NONE=exp( post$base+post$noedu ),										
			OR_SEXP1YR_G1=exp( post$base+post$sexp ),
			OR_VL_HIGH=exp( post$base+post$vlhigh ),
			OR_MOBILITY_INMIGRANT=exp( post$base+post$inmigrant ),
			OR_OCC_AGRO=exp( post$base+post$occ_agro ),
			OR_OCC_BAR=exp( post$base+post$occ_bar ),
			OR_OCC_BODA=exp( post$base+post$occ_boda ),
			OR_OCC_FISH=exp( post$base+post$occ_fish ),
			OR_OCC_STUD=exp( post$base+post$occ_stud ),
			OR_OCC_TRAD=exp( post$base+post$occ_trad ),
			ORX_GENDER_MALE=exp( post$male ),
			ORX_COMM_FISH=exp( post$comm_fish ), 											
			ORX_COMM_TRAD=exp( post$comm_trad ),										
			ORX_EDU_NONE=exp( post$noedu ),										
			ORX_SEXP1YR_G1=exp( post$sexp ),
			ORX_VL_HIGH=exp( post$vlhigh ),
			ORX_MOBILITY_INMIGRANT=exp( post$inmigrant ),
			ORX_OCC_AGRO=exp( post$occ_agro ),
			ORX_OCC_BAR=exp( post$occ_bar ),
			ORX_OCC_BODA=exp( post$occ_boda ),
			ORX_OCC_FISH=exp( post$occ_fish ),
			ORX_OCC_STUD=exp( post$occ_stud ),
			ORX_OCC_TRAD=exp( post$occ_trad )			
	)
	dp		<- melt(dp, id.vars='MC')	
	doh.1	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	doh.1	<- dcast.data.table(doh.1, variable~STAT, value.var='V')
	doh.1[, MODEL:= 'multivariate outside hh no age']
	
	#	model 2 with age (exchangeable)
	moh.2 <- map2stan(
			alist(
					TRANSMITTER ~ dbinom(1, ptr),
					logit(ptr) <- base + 
							# agrarian community is base 
							comm_fish*COMM_FISH + comm_trad*COMM_TRAD + 
							male_age[AGE_AT_MID_C2]*MALE + female_age[AGE_AT_MID_C2]*(1-MALE) + 
							noedu*NOEDU + sexp*SEXP1YR_G1 +	vlhigh*VL_HIGH + inmigrant*INMIGRANT +							
							# other occupation is baseline
							occ_agro*OC_AGRO + occ_fish*OC_FISH + occ_bar*OC_BAR + occ_boda*OC_BODA + occ_stud*OC_STUD + occ_trad*OC_TRAD,
					base ~ dnorm(0,100),
					male_age[AGE_AT_MID_C2] ~ dnorm(0, male_age_sig),
					female_age[AGE_AT_MID_C2] ~ dnorm(0, female_age_sig),
					c(comm_fish, comm_trad) ~ dnorm(0,10),
					c(noedu, sexp, vlhigh, inmigrant) ~ dnorm(0,10),								
					c(occ_agro, occ_bar, occ_boda, occ_stud, occ_trad, occ_fish) ~ dnorm(0,10),
					c(male_age_sig, female_age_sig) ~ dcauchy(0,1)
			),
			data=as.data.frame(df), 
			start=list(	base=0, comm_fish=0, comm_trad=0, male_age=rep(0,7), female_age=rep(0,7), noedu=0, sexp=0, vlhigh=0, inmigrant=0, 						
					occ_agro=0, occ_bar=0, occ_boda=0, occ_stud=0, occ_trad=0, occ_fish=0,
					male_age_sig=1, female_age_sig=1),			
			warmup=1e3, iter=5e3, chains=1, cores=4
		)
	#plot(precis(moh.2, prob=0.95, depth=2))
	#pairs(moh.2)
	post	<- extract.samples(moh.2)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_COMM_AGRO=logistic( post$base ),				
			PI_COMM_FISH=logistic( post$base+post$comm_fish ), 											
			PI_COMM_TRAD=logistic( post$base+post$comm_trad ),										
			PI_EDU_NONE=logistic( post$base+post$noedu ),										
			PI_SEXP1YR_G1=logistic( post$base+post$sexp ),
			PI_VL_HIGH=logistic( post$base+post$vlhigh ),
			PI_MOBILITY_INMIGRANT=logistic( post$base+post$inmigrant ),
			PI_OCC_AGRO=logistic( post$base+post$occ_agro ),
			PI_OCC_BAR=logistic( post$base+post$occ_bar ),
			PI_OCC_BODA=logistic( post$base+post$occ_boda ),
			PI_OCC_FISH=logistic( post$base+post$occ_fish ),
			PI_OCC_STUD=logistic( post$base+post$occ_stud ),
			PI_OCC_TRAD=logistic( post$base+post$occ_trad ),
			PI_MALEAGE_1=logistic( post$base+post$male_age[,1] ),
			PI_MALEAGE_2=logistic( post$base+post$male_age[,2] ),
			PI_MALEAGE_3=logistic( post$base+post$male_age[,3] ),
			PI_MALEAGE_4=logistic( post$base+post$male_age[,4] ),
			PI_MALEAGE_5=logistic( post$base+post$male_age[,5] ),
			PI_MALEAGE_6=logistic( post$base+post$male_age[,6] ),
			PI_MALEAGE_7=logistic( post$base+post$male_age[,7] ),
			PI_FEMALEAGE_1=logistic( post$base+post$female_age[,1] ),
			PI_FEMALEAGE_2=logistic( post$base+post$female_age[,2] ),
			PI_FEMALEAGE_3=logistic( post$base+post$female_age[,3] ),
			PI_FEMALEAGE_4=logistic( post$base+post$female_age[,4] ),
			PI_FEMALEAGE_5=logistic( post$base+post$female_age[,5] ),
			PI_FEMALEAGE_6=logistic( post$base+post$female_age[,6] ),
			PI_FEMALEAGE_7=logistic( post$base+post$female_age[,7] ),			
			OR_COMM_AGRO=exp( post$base ),				
			OR_COMM_FISH=exp( post$base+post$comm_fish ), 											
			OR_COMM_TRAD=exp( post$base+post$comm_trad ),										
			OR_EDU_NONE=exp( post$base+post$noedu ),										
			OR_SEXP1YR_G1=exp( post$base+post$sexp ),
			OR_VL_HIGH=exp( post$base+post$vlhigh ),
			OR_MOBILITY_INMIGRANT=exp( post$base+post$inmigrant ),
			OR_OCC_AGRO=exp( post$base+post$occ_agro ),
			OR_OCC_BAR=exp( post$base+post$occ_bar ),
			OR_OCC_BODA=exp( post$base+post$occ_boda ),
			OR_OCC_FISH=exp( post$base+post$occ_fish ),
			OR_OCC_STUD=exp( post$base+post$occ_stud ),
			OR_OCC_TRAD=exp( post$base+post$occ_trad ),
			OR_MALEAGE_1=exp( post$base+post$male_age[,1] ),
			OR_MALEAGE_2=exp( post$base+post$male_age[,2] ),
			OR_MALEAGE_3=exp( post$base+post$male_age[,3] ),
			OR_MALEAGE_4=exp( post$base+post$male_age[,4] ),
			OR_MALEAGE_5=exp( post$base+post$male_age[,5] ),
			OR_MALEAGE_6=exp( post$base+post$male_age[,6] ),
			OR_MALEAGE_7=exp( post$base+post$male_age[,7] ),
			OR_FEMALEAGE_1=exp( post$base+post$female_age[,1] ),
			OR_FEMALEAGE_2=exp( post$base+post$female_age[,2] ),
			OR_FEMALEAGE_3=exp( post$base+post$female_age[,3] ),
			OR_FEMALEAGE_4=exp( post$base+post$female_age[,4] ),
			OR_FEMALEAGE_5=exp( post$base+post$female_age[,5] ),
			OR_FEMALEAGE_6=exp( post$base+post$female_age[,6] ),
			OR_FEMALEAGE_7=exp( post$base+post$female_age[,7] ),							
			ORX_COMM_FISH=exp( post$comm_fish ), 											
			ORX_COMM_TRAD=exp( post$comm_trad ),										
			ORX_EDU_NONE=exp( post$noedu ),										
			ORX_SEXP1YR_G1=exp( post$sexp ),
			ORX_VL_HIGH=exp( post$vlhigh ),
			ORX_MOBILITY_INMIGRANT=exp( post$inmigrant ),
			ORX_OCC_AGRO=exp( post$occ_agro ),
			ORX_OCC_BAR=exp( post$occ_bar ),
			ORX_OCC_BODA=exp( post$occ_boda ),
			ORX_OCC_FISH=exp( post$occ_fish ),
			ORX_OCC_STUD=exp( post$occ_stud ),
			ORX_OCC_TRAD=exp( post$occ_trad ),
			ORX_MALEAGE_1=exp( post$male_age[,1]-post$female_age[,1] ),
			ORX_MALEAGE_2=exp( post$male_age[,2]-post$female_age[,2] ),
			ORX_MALEAGE_3=exp( post$male_age[,3]-post$female_age[,3] ),
			ORX_MALEAGE_4=exp( post$male_age[,4]-post$female_age[,4] ),
			ORX_MALEAGE_5=exp( post$male_age[,5]-post$female_age[,5] ),
			ORX_MALEAGE_6=exp( post$male_age[,6]-post$female_age[,6] ),
			ORX_MALEAGE_7=exp( post$male_age[,7]-post$female_age[,7] )
	)
	dp		<- melt(dp, id.vars='MC')	
	doh.2	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	doh.2	<- dcast.data.table(doh.2, variable~STAT, value.var='V')
	doh.2[, MODEL:= 'multivariate outside hh with age']
	
	
	ds		<- rbind(dhh.1, dhh.2, doh.1, doh.2)
	ds[, LABEL:= paste0(round(MED, d=2), '\n[', round(CL, d=2),'-', round( CU, d=2),']')]
	ds[, STAT:= factor(gsub('^([^_]+)_.*','\\1',variable), levels=c('PI','OR','ORX'), labels=c('probability transmitter','odds transmitter','odds ratio'))]
	ds[, FACTOR:=gsub('^([^_]+)_([^_]+)_([^_]+)','\\3',variable)]
	ds[, MOFA:=gsub('^([^_]+)_([^_]+)_([^_]+)','\\2-\\3',variable)]		
	#ds		<- subset(ds, STAT=='odds ratio')
	set(ds, NULL, 'MOFA2', ds[, factor(MOFA, levels=rev(c("GENDER-MALE",
											"MALEAGE-1","MALEAGE-2","MALEAGE-3","MALEAGE-4","MALEAGE-5","MALEAGE-6","MALEAGE-7",
											"FEMALEAGE-1","FEMALEAGE-2","FEMALEAGE-3","FEMALEAGE-4","FEMALEAGE-5","FEMALEAGE-6","FEMALEAGE-7",
											"COMM-AGRO","COMM-FISH","COMM-TRAD","EDU-NONE","SEXP1YR-G1","VL-HIGH","MOBILITY-INMIGRANT",
											"OCC-AGRO","OCC-BAR","OCC-BODA","OCC-FISH","OCC-STUD","OCC-TRAD"											
									)), labels=rev(c("male","male, age 15-19","male, age 20-24","male, age 25-29","male, age 30-34","male, age 35-39","male, age 40-44","male, age 45-50",
											"female, age 15-19","female, age 20-24","female, age 25-29","female, age 30-34","female, age 35-39","female, age 40-44","female, age 45-50",
											"Agrarian community", "Fishing site", "Trading centre",
											"No primary education","More than 1 sex partner in last year", "Viral load above 100,000 cps/ml", "Inmigrant",											
											"Primary occupation: agricultural/house", "Primary occupation: bar worker/waitress","Primary occupation: Boda/trucking",											
											"Primary occupation: fishing", "Primary occupation: student", "Primary occupation: trading/shopkeeper"
									)))])
	setkey(ds, MOFA)	
	#
	#	save
	#
	save(dg, ds, dhh.1, dhh.2, doh.1, doh.2, moh.1, moh.2, mhh.1, mhh.2, file=paste0(outfile.base,'_hh_runs.rda'))
	
	dp		<- subset(ds, STAT=='odds transmitter')
	dp[, FILL:= '0']
	set(dp, dp[, which(IL>1 | IU<1)], 'FILL', '1')
	set(dp, dp[, which(CL>1 | CU<1)], 'FILL', '2')
	ggplot(subset(dp, MODEL=='multivariate hh with age'), aes(x=MOFA2)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/20, 1/10, 1/6, 1/3,1/2,1,2,3, 6,10,20), labels=c('1/20','1/10','1/6','1/3','1/2','1','2','3','6','10', '20'), expand=c(0,0)) +
			scale_fill_manual(values=c('0'='white', '1'='orange', '2'='red')) +
			coord_flip(ylim=c(1/20,20)) +
			guides(fill=FALSE) +
			labs(x='Individual characteristics\n', y='\nOdds for being a transmitter vs. recipient\nin phylogenetically reconstructed transmission pairs')
	ggsave(file=paste0(outfile.base,'_odds_samehh_withage.pdf'), w=6, h=8)
	ggplot(subset(dp, MODEL=='multivariate outside hh with age'), aes(x=MOFA2)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/20, 1/10, 1/6, 1/3,1/2,1,2,3, 6,10,20), labels=c('1/20','1/10','1/6','1/3','1/2','1','2','3','6','10', '20'), expand=c(0,0)) +
			scale_fill_manual(values=c('0'='white', '1'='orange', '2'='red')) +
			coord_flip(ylim=c(1/20,20)) +
			guides(fill=FALSE) +
			labs(x='Individual characteristics\n', y='\nOdds for being a transmitter vs. recipient\nin phylogenetically reconstructed transmission pairs')
	ggsave(file=paste0(outfile.base,'_odds_diffhh_withage.pdf'), w=6, h=8)
	
	dp		<- subset(ds, STAT=='odds ratio')
	dp[, FILL:= '0']
	set(dp, dp[, which(IL>1 | IU<1)], 'FILL', '1')
	set(dp, dp[, which(CL>1 | CU<1)], 'FILL', '2')	
	ggplot(subset(dp, MODEL=='multivariate hh no age'), aes(x=MOFA2)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/20, 1/10, 1/6, 1/3,1/2,1,2,3, 6,10,20), labels=c('1/20','1/10','1/6','1/3','1/2','1','2','3','6','10', '20'), expand=c(0,0)) +
			scale_fill_manual(values=c('0'='white', '1'='orange', '2'='red')) +
			coord_flip(ylim=c(1/20,20)) +
			guides(fill=FALSE) +
			labs(x='Individual characteristics\n', y='\nOdds ratio for being a transmitter vs. recipient\nin phylogenetically reconstructed transmission pairs')
	ggsave(file=paste0(outfile.base,'_odds_samehh_noage.pdf'), w=6, h=5)
	ggplot(subset(dp, MODEL=='multivariate outside hh no age'), aes(x=MOFA2)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/20, 1/10, 1/6, 1/3,1/2,1,2,3, 6,10,20), labels=c('1/20','1/10','1/6','1/3','1/2','1','2','3','6','10', '20'), expand=c(0,0)) +
			scale_fill_manual(values=c('0'='white', '1'='orange', '2'='red')) +
			coord_flip(ylim=c(1/20,20)) +
			guides(fill=FALSE) +
			labs(x='Individual characteristics\n', y='\nOdds ratio for being a transmitter vs. recipient\nin phylogenetically reconstructed transmission pairs')
	ggsave(file=paste0(outfile.base,'_odds_diffhh_noage.pdf'), w=6, h=5)
}

RakaiFull.gender.171122.nonhh.multivariatemodels.stan.with.threshold<- function()
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
	
	#zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	#zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	load(infile)	
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	rtpdm[, PAIR_COMM_TYPE:= FEMALE_COMM_TYPE]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_TYPE!=MALE_COMM_TYPE)], 'PAIR_COMM_TYPE', 'mixed')
	set(rtpdm, rtpdm[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(rtpdm, rtpdm[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')
	
	#rtpdm[, table(PAIR_COMM_TYPE)]
	rtpdm[, COUPLE2:= factor(COUPLE=='no couple', levels=c(TRUE,FALSE), labels=c('no couple','couple'))]
	rtpdm[, SAMEHH:= factor(FEMALE_HH_NUM==MALE_HH_NUM, levels=c(TRUE,FALSE), labels=c('same hh','different hh'))]	
	rtpdm[, PAIR_COMM:= MALE_COMM_NUM_A]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_NUM_A!=MALE_COMM_NUM_A)], 'PAIR_COMM', 'mixed')
	rtpdm[, MALE_SEXP1OUT2:= factor(MALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(MALE_SEXP1OUT=='Unknown')], 'MALE_SEXP1OUT2', 'Unknown')
	rtpdm[, FEMALE_SEXP1OUT2:= factor(FEMALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(FEMALE_SEXP1OUT=='Unknown')], 'FEMALE_SEXP1OUT2', 'Unknown')
	
	
	df		<- subset(rtpdm, FEMALE_HH_NUM!=MALE_HH_NUM, select=c(	MALE_RID, FEMALE_RID, PTY_RUN, IDCLU, LINK_MF, POSTERIOR_SCORE_MF, COUPLE2, SAMEHH, PAIR_COMM_TYPE, FEMALE_COMM_NUM,
					MALE_RECENTVL, MALE_SEXP1YR, MALE_SEXP1OUT2, MALE_OCCUP_OLLI, MALE_OCAT, MALE_EDUCAT, MALE_CIRCUM,
					FEMALE_RECENTVL, FEMALE_SEXP1YR, FEMALE_SEXP1OUT2, FEMALE_OCCUP_OLLI, FEMALE_OCAT, FEMALE_EDUCAT 
				))
	#	missing data: fill in
	set(df, df[, which(is.na(MALE_RECENTVL))], 'MALE_RECENTVL', -1)
	set(df, df[, which(is.na(FEMALE_RECENTVL))], 'FEMALE_RECENTVL', -1)
	set(df, df[, which(is.na(MALE_CIRCUM))], 'MALE_CIRCUM', 'Unknown')
	set(df, df[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(df, df[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')	
	#for(x in colnames(df)) print( c(x, any(is.na(df[[x]]))) )
	#	prepare data for STAN

	#	add estimated HIV prevalences (medians)
	infile	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30trmMF_estimated_prevalenceratio.rda"	
	load(infile)
	tmp		<- dcast.data.table(dgg, COMM_NUM~variable, value.var='M')
	setnames(tmp, 'COMM_NUM', 'FEMALE_COMM_NUM')
	df		<- merge(df, tmp, by='FEMALE_COMM_NUM')

	df[, COMM_FISH:= as.integer(PAIR_COMM_TYPE=='fisherfolk')]
	df[, COMM_AGR:= as.integer(PAIR_COMM_TYPE=='agrarian')]
	df[, COMM_TRAD:= as.integer(PAIR_COMM_TYPE=='trading')]
	df[, COMM_MXD:= as.integer(PAIR_COMM_TYPE=='mixed')]
	df[, FE_NOEDU:= as.integer(FEMALE_EDUCAT=='None')]
	df[, FE_NOEDU_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, MA_NOEDU:= as.integer(MALE_EDUCAT=='None')]
	df[, MA_NOEDU_MISS:= as.integer(MALE_EDUCAT=='Unknown')]	
	df[, MA_CIRCUM:= as.integer(MALE_CIRCUM=='Y')]
	df[, MA_CIRCUM_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, FE_SEXP1YR_G1:= as.integer(FEMALE_SEXP1YR!='1')]
	df[, MA_SEXP1YR_G1:= as.integer(MALE_SEXP1YR!='1')]
	df[, MA_OC_AGRO:= as.integer(MALE_OCAT=='Agro/House')]	
	df[, MA_OC_FISH:= as.integer(MALE_OCAT=='Fishing')]
	df[, MA_OC_TRAD:= as.integer(MALE_OCAT=='Trading/Shop keeper')]
	df[, MA_OC_OTH:= as.integer(MALE_OCAT%in%c('zother','Bar/waitress','Student','Boda/Trucking'))]	
	df[, FE_OC_AGRO:= as.integer(FEMALE_OCAT=='Agro/House')]
	df[, FE_OC_BAR:= as.integer(FEMALE_OCAT=='Bar/waitress')]			
	df[, FE_OC_TRAD:= as.integer(FEMALE_OCAT=='Trading/Shop keeper')]
	df[, FE_OC_OTH:= as.integer(FEMALE_OCAT%in%c('zother','Student'))]		
	df[, RFM_ABOVE_MEAN:=  as.integer( RFM>=subset(dgg, variable=='RFM')[, mean(M)] ) ]
	#
	#	STAN
	#
	
	#	same household definitely stronger effect than couples	
	ms.2 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + 
							# comm_fish*COMM_FISH is base 
							comm_agr*COMM_AGR + comm_trad*COMM_TRAD + comm_mxd*COMM_MXD +
							noedu_f*FE_NOEDU + noedu_m*MA_NOEDU +
							sexp_f*FE_SEXP1YR_G1 + sexp_m*MA_SEXP1YR_G1 +
							circum_m*MA_CIRCUM +	
							prevratio_fm * RFM_ABOVE_MEAN +
							# agro occupation is baseline
							mocc_other*MA_OC_OTH + mocc_fish*MA_OC_FISH + mocc_trad*MA_OC_TRAD +										
							focc_other*FE_OC_OTH + focc_bar*FE_OC_BAR +  focc_trad*FE_OC_TRAD,
					base ~ dnorm(0,100),								
					c(comm_agr, comm_trad, comm_mxd) ~ dnorm(0,10),
					c(noedu_f, noedu_m) ~ dnorm(0,10),								
					c(sexp_f, sexp_m, circum_m, prevratio_fm) ~ dnorm(0,10),					
					c(mocc_other, mocc_fish, mocc_trad) ~ dnorm(0,10),
					c(focc_other, focc_bar, focc_trad) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, comm_agr=0, comm_trad=0, comm_mxd=0, noedu_f=0, noedu_m=0, sexp_f=0, sexp_m=0, circum_m=0,						
					prevratio_fm=0, focc_other=0, focc_bar=0, focc_trad=0, mocc_other=0, mocc_fish=0, mocc_trad=0),			
			warmup=1e3, iter=5e3, chains=1, cores=4
	)		
	#plot(precis(ms.2, prob=0.95))
	#pairs(ms.2)
	post	<- extract.samples(ms.2)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_COMM_FISH=logistic( post$base ),				
			PI_COMM_AGRO=logistic( post$base+post$comm_agr ), 											
			PI_COMM_TRAD=logistic( post$base+post$comm_trad ),							
			PI_COMM_MXD=logistic( post$base+post$comm_mxd ),
			PI_FEDU_NONE=logistic( post$base+post$noedu_f ),
			PI_MEDU_NONE=logistic( post$base+post$noedu_m ),							
			PI_FESEXP1YR_G1=logistic( post$base+post$sexp_f ),
			PI_MASEXP1YR_G1=logistic( post$base+post$sexp_m ),							
			PI_CIRC_YES=logistic( post$base+post$circum_m ),							
			PI_MOCC_OTH=logistic( post$base+post$mocc_other ),
			PI_MOCC_FISH=logistic( post$base+post$mocc_fish ),
			PI_MOCC_TRAD=logistic( post$base+post$mocc_trad ),
			PI_FOCC_OTH=logistic( post$base+post$focc_other ),
			PI_FOCC_BAR=logistic( post$base+post$focc_bar ),
			PI_FOCC_TRAD=logistic( post$base+post$focc_trad ),	
			PI_PREV_FMRATIO=logistic( post$base+post$prevratio_fm ),
			OR_COMM_FISH=exp( post$base ),				
			OR_COMM_AGRO=exp( post$base+post$comm_agr ), 											
			OR_COMM_TRAD=exp( post$base+post$comm_trad ),							
			OR_COMM_MXD=exp( post$base+post$comm_mxd ),
			OR_FEDU_NONE=exp( post$base+post$noedu_f ),
			OR_MEDU_NONE=exp( post$base+post$noedu_m ),							
			OR_FESEXP1YR_G1=exp( post$base+post$sexp_f ),
			OR_MASEXP1YR_G1=exp( post$base+post$sexp_m ),
			OR_CIRC_YES=exp( post$base+post$circum_m ),							
			OR_MOCC_OTH=exp( post$base+post$mocc_other ),
			OR_MOCC_FISH=exp( post$base+post$mocc_fish ),
			OR_MOCC_TRAD=exp( post$base+post$mocc_trad ),
			OR_FOCC_OTH=exp( post$base+post$focc_other ),
			OR_FOCC_BAR=exp( post$base+post$focc_bar ),
			OR_FOCC_TRAD=exp( post$base+post$focc_trad ),	
			OR_PREV_FMRATIO=exp( post$base+post$prevratio_fm ),
			ORX_COMM_AGRO=exp( post$comm_agr ), 											
			ORX_COMM_TRAD=exp( post$comm_trad ),							
			ORX_COMM_MXD=exp( post$comm_mxd ),
			ORX_FEDU_NONE=exp( post$noedu_f ),
			ORX_MEDU_NONE=exp( post$noedu_m ),
			ORX_FESEXP1YR_G1=exp( post$sexp_f ),
			ORX_MASEXP1YR_G1=exp( post$sexp_m ),							
			ORX_CIRC_YES=exp( post$circum_m ),							
			ORX_MOCC_OTH=exp( post$mocc_other ),
			ORX_MOCC_FISH=exp( post$mocc_fish ),
			ORX_MOCC_TRAD=exp( post$mocc_trad ),
			ORX_FOCC_OTH=exp( post$focc_other ),
			ORX_FOCC_BAR=exp( post$focc_bar ),
			ORX_FOCC_TRAD=exp( post$focc_trad ),
			ORX_PREV_FMRATIO=exp( post$prevratio_fm )
	)
	dp		<- melt(dp, id.vars='MC')	
	dss.2	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	dss.2	<- dcast.data.table(dss.2, variable~STAT, value.var='V')
	dss.2[, MODEL:= 'multivariate ms2']
	ds		<- copy(dss.2)
	ds[, LABEL:= paste0(round(MED, d=2), '\n[', round(CL, d=2),'-', round( CU, d=2),']')]
	ds[, STAT:= factor(gsub('^([^_]+)_.*','\\1',variable), levels=c('PI','OR','ORX'), labels=c('proportion MF','odds MF','odds ratio'))]
	ds[, FACTOR:=gsub('^([^_]+)_([^_]+)_([^_]+)','\\3',variable)]
	ds[, MOFA:=gsub('^([^_]+)_([^_]+)_([^_]+)','\\2-\\3',variable)]
	ds[, MODELTYPE:= factor(grepl('univariate',MODEL), levels=c(TRUE,FALSE),labels=c('univariate','multivariate'))]	
	ds		<- subset(ds, STAT=='odds ratio')
	set(ds, NULL, 'MOFA2', ds[, factor(MOFA, levels=rev(c(
											"PREV-FMRATIO",
											"COMM-AGRO","COMM-TRAD","COMM-MXD",																						
											"FEDU-NONE", "MEDU-NONE",      
											"FESEXP1YR-G1","MASEXP1YR-G1", 
											"CIRC-YES",                   
											"FOCC-BAR","FOCC-TRAD","FOCC-OTH",
											"MOCC-FISH","MOCC-OTH","MOCC-TRAD"   
									)), labels=rev(c("female-to-male prevalence ratio above average vs. below average",
											"Agrarian community vs. fishing community", "Trading community vs. fishing community", "Mixed vs. fishing community",
											"Female education: None vs. at least primary education", "Male education: None vs. at least primary education",
											"Female sex partners in last year: >1 vs. 1", "Male sex partners in last year: >1 vs. 1",
											"Male circumcised: yes vs. no",
											"Female occupation: bar/waitress vs. agricultural/house", "Female occupation: trading/shopkeeper vs. agricultural/house","Female occupation: other vs. agricultural/house", 
											"Male occupation: fishing vs. agricultural/house", "Male occupation: other vs. agricultural/house", "Male occupation: trading/shopkeeper vs. agricultural/house"
									)))])
	setkey(ds, MOFA)	
	dp		<- subset(ds, !MOFA%in%c('COMM-TRAD','COMM-MXD'))
	dp[, FILL:= '0']
	set(dp, dp[, which(IL>1 | IU<1)], 'FILL', '1')
	set(dp, dp[, which(CL>1 | CU<1)], 'FILL', '2')
	ggplot(dp, aes(x=MOFA2)) +
			geom_hline(yintercept=1, colour='grey50', lwd=1) +
			geom_boxplot(aes(middle=MED, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=FILL), stat='identity') +
			theme_bw() +
			scale_y_continuous(trans='log', breaks=c(1/20, 1/10, 1/6, 1/3,1/2,1,2,3, 6,10,20), labels=c('1/20','1/10','1/6','1/3','1/2','1','2','3','6','10', '20'), expand=c(0,0)) +
			scale_fill_manual(values=c('0'='white', '1'='orange', '2'='red')) +
			coord_flip(ylim=c(1/20,20)) +
			guides(fill=FALSE) +
			labs(x='Partner characteristics\n', y='\nOdds ratio of male-to-female transmission\nbetween partner characteristics')
	ggsave(file=paste0(outfile.base,'_propmf_factors_univariate_diffhh.pdf'), w=8, h=3.5)		
}

RakaiFull.gender.171122.couples.multivariateadditivemodels.stan.with.threshold<- function()
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
	
	#zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	#zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	load(infile)	
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	rtpdm[, PAIR_COMM_TYPE:= FEMALE_COMM_TYPE]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_TYPE!=MALE_COMM_TYPE)], 'PAIR_COMM_TYPE', 'mixed')
	set(rtpdm, rtpdm[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(rtpdm, rtpdm[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')
	
	#rtpdm[, table(PAIR_COMM_TYPE)]
	rtpdm[, COUPLE2:= factor(COUPLE=='no couple', levels=c(TRUE,FALSE), labels=c('no couple','couple'))]
	rtpdm[, SAMEHH:= factor(FEMALE_HH_NUM==MALE_HH_NUM, levels=c(TRUE,FALSE), labels=c('same hh','different hh'))]	
	rtpdm[, PAIR_COMM:= MALE_COMM_NUM_A]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_NUM_A!=MALE_COMM_NUM_A)], 'PAIR_COMM', 'mixed')
	rtpdm[, MALE_SEXP1OUT2:= factor(MALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(MALE_SEXP1OUT=='Unknown')], 'MALE_SEXP1OUT2', 'Unknown')
	rtpdm[, FEMALE_SEXP1OUT2:= factor(FEMALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(FEMALE_SEXP1OUT=='Unknown')], 'FEMALE_SEXP1OUT2', 'Unknown')
	
	
	df		<- subset(rtpdm, select=c(	MALE_RID, FEMALE_RID, PTY_RUN, IDCLU, LINK_MF, POSTERIOR_SCORE_MF, COUPLE2, SAMEHH, PAIR_COMM_TYPE,FEMALE_COMM_NUM, 
					MALE_RECENTVL, MALE_SEXP1YR, MALE_SEXP1OUT2, MALE_OCCUP_OLLI, MALE_OCAT, MALE_EDUCAT, MALE_CIRCUM,
					FEMALE_RECENTVL, FEMALE_SEXP1YR, FEMALE_SEXP1OUT2, FEMALE_OCCUP_OLLI, FEMALE_OCAT, FEMALE_EDUCAT 
			))
	#	missing data: fill in
	set(df, df[, which(is.na(MALE_RECENTVL))], 'MALE_RECENTVL', -1)
	set(df, df[, which(is.na(FEMALE_RECENTVL))], 'FEMALE_RECENTVL', -1)
	set(df, df[, which(is.na(MALE_CIRCUM))], 'MALE_CIRCUM', 'Unknown')
	set(df, df[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(df, df[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')
	#for(x in colnames(df)) print( c(x, any(is.na(df[[x]]))) )
	#	prepare data for STAN
	
	#	add estimated HIV prevalences (medians)
	infile	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30trmMF_estimated_prevalenceratio.rda"	
	load(infile)
	tmp		<- dcast.data.table(dgg, COMM_NUM~variable, value.var='M')
	setnames(tmp, 'COMM_NUM', 'FEMALE_COMM_NUM')
	df		<- merge(df, tmp, by='FEMALE_COMM_NUM')
	
	df[, COUPLE3:= as.integer(COUPLE2=='couple')]
	df[, SAMEHH3:= as.integer(SAMEHH=='same hh')]	
	df[, COMM_FISH:= as.integer(PAIR_COMM_TYPE=='fisherfolk')]
	df[, COMM_AGR:= as.integer(PAIR_COMM_TYPE=='agrarian')]
	df[, COMM_TRAD:= as.integer(PAIR_COMM_TYPE=='trading')]
	df[, COMM_MXD:= as.integer(PAIR_COMM_TYPE=='mixed')]
	df[, FE_NOEDU:= as.integer(FEMALE_EDUCAT=='None')]
	df[, FE_NOEDU_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, MA_NOEDU:= as.integer(MALE_EDUCAT=='None')]
	df[, MA_NOEDU_MISS:= as.integer(MALE_EDUCAT=='Unknown')]	
	df[, MA_CIRCUM:= as.integer(MALE_CIRCUM=='Y')]
	df[, MA_CIRCUM_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, FE_SEXP1YR_G1:= as.integer(FEMALE_SEXP1YR!='1')]
	df[, MA_SEXP1YR_G1:= as.integer(MALE_SEXP1YR!='1')]
	df[, MA_OC_AGRO:= as.integer(MALE_OCAT=='Agro/House')]	
	df[, MA_OC_FISH:= as.integer(MALE_OCAT=='Fishing')]
	df[, MA_OC_TRAD:= as.integer(MALE_OCAT=='Trading/Shop keeper')]
	df[, MA_OC_OTH:= as.integer(MALE_OCAT%in%c('zother','Bar/waitress','Student','Boda/Trucking'))]	
	df[, FE_OC_AGRO:= as.integer(FEMALE_OCAT=='Agro/House')]
	df[, FE_OC_BAR:= as.integer(FEMALE_OCAT=='Bar/waitress')]			
	df[, FE_OC_TRAD:= as.integer(FEMALE_OCAT=='Trading/Shop keeper')]
	df[, FE_OC_OTH:= as.integer(FEMALE_OCAT%in%c('zother','Student'))]		
	df[, PF2:= PF-mean(PF)]
	df[, PM2:= PM-mean(PM)]
	df[, RFM2:= RFM-mean(RFM)]
	df[, RFM_ABOVE_MEAN:=  as.integer( RFM>=subset(dgg, variable=='RFM')[, mean(M)] ) ]
	#
	#	STAN
	#
	mm.2 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + couple_b*COUPLE3 + samehh_b*SAMEHH3 +
							# comm_fish*COMM_FISH +	this is base 
							comm_agr*COMM_AGR + comm_trad*COMM_TRAD + comm_mxd*COMM_MXD +
							noedu_f*FE_NOEDU + noedu_m*MA_NOEDU +
							sexp_f*FE_SEXP1YR_G1 + sexp_m*MA_SEXP1YR_G1 +
							circum_m*MA_CIRCUM +
							prevratio_fm * RFM_ABOVE_MEAN +
							# other occupation is baseline
							mocc_agro*MA_OC_AGRO + mocc_fish*MA_OC_FISH + mocc_trad*MA_OC_TRAD +										
							focc_agro*FE_OC_AGRO + focc_bar*FE_OC_BAR +  focc_trad*FE_OC_TRAD,
					base ~ dnorm(0,100),
					c(couple_b, samehh_b) ~ dnorm(0,10),
					c(comm_agr, comm_trad, comm_mxd) ~ dnorm(0,10),
					c(noedu_f, noedu_m) ~ dnorm(0,10),								
					c(sexp_f, sexp_m, circum_m) ~ dnorm(0,10),
					c(prevratio_fm) ~ dnorm(0,10),
					c(mocc_agro, mocc_fish, mocc_trad) ~ dnorm(0,10),
					c(focc_agro, focc_bar, focc_trad) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, couple_b=0, samehh_b=0, comm_agr=0, comm_trad=0, comm_mxd=0, noedu_f=0, noedu_m=0, sexp_f=0, sexp_m=0, circum_m=0,
					prevratio_fm=0,
					focc_agro=0, focc_bar=0, focc_trad=0,
					mocc_agro=0, mocc_fish=0, mocc_trad=0),			
			warmup=1e3, iter=5e3, chains=1, cores=4
		)
	
	#plot(precis(mm.2, prob=0.95))
	#pairs(mm.2)
	post	<- extract.samples(mm.2)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_COMM_FISH=logistic( post$base ),				
			PI_COMM_AGRO=logistic( post$base+post$comm_agr ), 											
			PI_COMM_TRAD=logistic( post$base+post$comm_trad ),							
			PI_COMM_MXD=logistic( post$base+post$comm_mxd ),
			PI_COUPLE_YES=logistic( post$base+post$couple_b ),
			PI_SAMEHH_YES=logistic( post$base+post$samehh_b ),
			PI_FEDU_NONE=logistic( post$base+post$noedu_f ),
			PI_MEDU_NONE=logistic( post$base+post$noedu_m ),							
			PI_FESEXP1YR_G1=logistic( post$base+post$sexp_f ),
			PI_MASEXP1YR_G1=logistic( post$base+post$sexp_m ),							
			PI_CIRC_YES=logistic( post$base+post$circum_m ),							
			PI_MOCC_AGRO=logistic( post$base+post$mocc_agro ),
			PI_MOCC_FISH=logistic( post$base+post$mocc_fish ),
			PI_MOCC_TRAD=logistic( post$base+post$mocc_trad ),
			PI_FOCC_AGRO=logistic( post$base+post$focc_agro ),
			PI_FOCC_BAR=logistic( post$base+post$focc_bar ),
			PI_FOCC_TRAD=logistic( post$base+post$focc_trad ),
			PI_PREVFMRATIO_ABOVE=logistic( post$base+post$prevratio_fm ),
			PI_PREVFMRATIO_BELOW=logistic( post$base ),
			OR_COMM_FISH=exp( post$base ),				
			OR_COMM_AGRO=exp( post$base+post$comm_agr ), 											
			OR_COMM_TRAD=exp( post$base+post$comm_trad ),							
			OR_COMM_MXD=exp( post$base+post$comm_mxd ),
			OR_COUPLE_YES=exp( post$base+post$couple_b ),
			OR_SAMEHH_YES=exp( post$base+post$samehh_b ),
			OR_FEDU_NONE=exp( post$base+post$noedu_f ),
			OR_MEDU_NONE=exp( post$base+post$noedu_m ),							
			OR_FESEXP1YR_G1=exp( post$base+post$sexp_f ),
			OR_MASEXP1YR_G1=exp( post$base+post$sexp_m ),
			OR_CIRC_YES=exp( post$base+post$circum_m ),							
			OR_MOCC_AGRO=exp( post$base+post$mocc_agro ),
			OR_MOCC_FISH=exp( post$base+post$mocc_fish ),
			OR_MOCC_TRAD=exp( post$base+post$mocc_trad ),
			OR_FOCC_AGRO=exp( post$base+post$focc_agro ),
			OR_FOCC_BAR=exp( post$base+post$focc_bar ),
			OR_FOCC_TRAD=exp( post$base+post$focc_trad ),
			OR_PREVFMRATIO_ABOVE=exp( post$base+post$prevratio_fm ),
			OR_PREVFMRATIO_BELOW=exp( post$base ),
			ORX_COMM_AGRO=exp( post$comm_agr ), 											
			ORX_COMM_TRAD=exp( post$comm_trad ),							
			ORX_COMM_MXD=exp( post$comm_mxd ),
			ORX_COUPLE_YES=exp( post$couple_b ),
			ORX_SAMEHH_YES=exp( post$samehh_b ),
			ORX_FEDU_NONE=exp( post$noedu_f ),
			ORX_MEDU_NONE=exp( post$noedu_m ),
			ORX_FESEXP1YR_G1=exp( post$sexp_f ),
			ORX_MASEXP1YR_G1=exp( post$sexp_m ),							
			ORX_CIRC_YES=exp( post$circum_m ),							
			ORX_MOCC_AGRO=exp( post$mocc_agro ),
			ORX_MOCC_FISH=exp( post$mocc_fish ),
			ORX_MOCC_TRAD=exp( post$mocc_trad ),
			ORX_FOCC_AGRO=exp( post$focc_agro ),
			ORX_FOCC_BAR=exp( post$focc_bar ),
			ORX_FOCC_TRAD=exp( post$focc_trad ),
			ORX_PREVFMRATIO_ABOVE=exp( post$prevratio_fm )
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.2	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.2	<- dcast.data.table(ds.2, variable~STAT, value.var='V')
	ds.2[, MODEL:= 'multivariate mm2']	
	#
	#	univariate couples
	#
	mu.1 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + couple_b*COUPLE3,
					base ~ dnorm(0,100),
					couple_b ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, couple_b=0),			
			warmup=500, iter=2e3, chains=1, cores=4
	)				
	post	<- extract.samples(mu.1)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_COUPLE_NO=logistic( post$base ),
			PI_COUPLE_YES=logistic( post$base+post$couple_b ),
			OR_COUPLE_NO=exp( post$base ),
			OR_COUPLE_YES=exp( post$base+post$couple_b ),
			ORX_COUPLE_YES=exp( post$couple_b )							
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.1	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.1	<- dcast.data.table(ds.1, variable~STAT, value.var='V')
	ds.1[, MODEL:= 'univariate couples']	
	#
	#	univariate household
	#
	mu.3 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + samehh_b*SAMEHH3,
					base ~ dnorm(0,100),
					samehh_b ~ dnorm(0,10)					
			),
			data=as.data.frame(df), 
			start=list(	base=0, samehh_b=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.3)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_SAMEHH_NO=logistic( post$base ),
			PI_SAMEHH_YES=logistic( post$base+post$samehh_b ),
			OR_SAMEHH_NO=exp( post$base ),
			OR_SAMEHH_YES=exp( post$base+post$samehh_b ),
			ORX_SAMEHH_YES=exp( post$samehh_b )														
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.3	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.3	<- dcast.data.table(ds.3, variable~STAT, value.var='V')
	ds.3[, MODEL:= 'univariate same household']	
	#
	#	univariate community type
	#
	mu.4 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + comm_agr*COMM_AGR + comm_trad*COMM_TRAD + comm_mxd*COMM_MXD,
					base ~ dnorm(0,100),							
					c(comm_agr, comm_trad, comm_mxd) ~ dnorm(0,10)							
			),
			data=as.data.frame(df), 
			start=list(	base=0, comm_agr=0, comm_trad=0, comm_mxd=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.4)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_COMM_FISH=logistic( post$base ),				
			PI_COMM_AGRO=logistic( post$base+post$comm_agr ), 											
			PI_COMM_TRAD=logistic( post$base+post$comm_trad ),							
			PI_COMM_MXD=logistic( post$base+post$comm_mxd ),
			OR_COMM_FISH=exp( post$base ),				
			OR_COMM_AGRO=exp( post$base+post$comm_agr ), 											
			OR_COMM_TRAD=exp( post$base+post$comm_trad ),							
			OR_COMM_MXD=exp( post$base+post$comm_mxd ),				
			ORX_COMM_AGRO=exp( post$comm_agr ), 											
			ORX_COMM_TRAD=exp( post$comm_trad ),							
			ORX_COMM_MXD=exp( post$comm_mxd )							
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.4	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.4	<- dcast.data.table(ds.4, variable~STAT, value.var='V')
	ds.4[, MODEL:= 'univariate commtype']	
	#
	#	univariate female education
	#
	mu.5 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + noedu_f*FE_NOEDU,
					base ~ dnorm(0,100),
					noedu_f ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, noedu_f=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.5)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_FEDU_YES=logistic( post$base ),
			PI_FEDU_NONE=logistic( post$base+post$noedu_f ),
			OR_FEDU_YES=exp( post$base ),
			OR_FEDU_NONE=exp( post$base+post$noedu_f ),
			ORX_FEDU_NONE=exp( post$noedu_f)	)
	dp		<- melt(dp, id.vars='MC')	
	ds.5	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.5	<- dcast.data.table(ds.5, variable~STAT, value.var='V')
	ds.5[, MODEL:= 'univariate female education']	
	#
	#	univariate male education
	#
	mu.6 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + noedu_m*MA_NOEDU,
					base ~ dnorm(0,100),
					noedu_m ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, noedu_m=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.6)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_MEDU_YES=logistic( post$base ),
			PI_MEDU_NONE=logistic( post$base+post$noedu_m ),
			OR_MEDU_YES=exp( post$base ),
			OR_MEDU_NONE=exp( post$base+post$noedu_m ),
			ORX_MEDU_NONE=exp( post$noedu_m)	)
	dp		<- melt(dp, id.vars='MC')	
	ds.6	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.6	<- dcast.data.table(ds.6, variable~STAT, value.var='V')
	ds.6[, MODEL:= 'univariate male education']	
	#
	#	univariate female sex partners
	#
	mu.7 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + sexp_f*FE_SEXP1YR_G1 ,
					base ~ dnorm(0,100),
					sexp_f ~ dnorm(0,10)					
			),
			data=as.data.frame(df), 
			start=list(	base=0, sexp_f=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.7)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_FESEXP1YR_ONE=logistic( post$base ),
			PI_FESEXP1YR_G1=logistic( post$base+post$sexp_f ),
			OR_FESEXP1YR_ONE=exp( post$base ),
			OR_FESEXP1YR_G1=exp( post$base+post$sexp_f ),
			ORX_FESEXP1YR_G1=exp( post$sexp_f )	)
	dp		<- melt(dp, id.vars='MC')	
	ds.7	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.7	<- dcast.data.table(ds.7, variable~STAT, value.var='V')
	ds.7[, MODEL:= 'univariate female sex partners']	
	#
	#	univariate male sex partners
	#
	mu.8 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + sexp_m*MA_SEXP1YR_G1 ,
					base ~ dnorm(0,100),
					sexp_m ~ dnorm(0,10)					
			),
			data=as.data.frame(df), 
			start=list(	base=0, sexp_m=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.8)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_MASEXP1YR_ONE=logistic( post$base ),
			PI_MASEXP1YR_G1=logistic( post$base+post$sexp_m ),
			OR_MASEXP1YR_ONE=exp( post$base ),
			OR_MASEXP1YR_G1=exp( post$base+post$sexp_m ),
			ORX_MASEXP1YR_G1=exp( post$sexp_m )	)
	dp		<- melt(dp, id.vars='MC')	
	ds.8	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.8	<- dcast.data.table(ds.8, variable~STAT, value.var='V')
	ds.8[, MODEL:= 'univariate male sex partners']	
	#
	#	univariate circumcision
	#
	mu.9 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + circum_m*MA_CIRCUM,
					base ~ dnorm(0,100),
					circum_m ~ dnorm(0,10)					
			),
			data=as.data.frame(df), 
			start=list(	base=0, circum_m=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.9)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_CIRC_NO=logistic( post$base ),
			PI_CIRC_YES=logistic( post$base+post$circum_m ),							
			OR_CIRC_NO=exp( post$base ),
			OR_CIRC_YES=exp( post$base+post$circum_m ),							
			ORX_CIRC_YES=exp( post$circum_m )
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.9	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.9	<- dcast.data.table(ds.9, variable~STAT, value.var='V')
	ds.9[, MODEL:= 'univariate circumcision']	
	#
	#	univariate male occupation
	#
	mu.10 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + mocc_agro*MA_OC_AGRO + mocc_fish*MA_OC_FISH + mocc_trad*MA_OC_TRAD,
					base ~ dnorm(0,100),
					c(mocc_agro, mocc_fish, mocc_trad) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, mocc_agro=0, mocc_fish=0, mocc_trad=0),			
			warmup=1e3, iter=5e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.10)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_MOCC_OTH=logistic( post$base ),
			PI_MOCC_AGRO=logistic( post$base+post$mocc_agro ),
			PI_MOCC_FISH=logistic( post$base+post$mocc_fish ),
			PI_MOCC_TRAD=logistic( post$base+post$mocc_trad ),
			OR_MOCC_OTH=exp( post$base ),
			OR_MOCC_AGRO=exp( post$base+post$mocc_agro ),
			OR_MOCC_FISH=exp( post$base+post$mocc_fish ),
			OR_MOCC_TRAD=exp( post$base+post$mocc_trad ),
			ORX_MOCC_AGRO=exp( post$mocc_agro ),
			ORX_MOCC_FISH=exp( post$mocc_fish ),
			ORX_MOCC_TRAD=exp( post$mocc_trad )	)									
	dp		<- melt(dp, id.vars='MC')	
	ds.10	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.10	<- dcast.data.table(ds.10, variable~STAT, value.var='V')
	ds.10[, MODEL:= 'univariate male occup']	
	#
	#	univariate female occupation
	#
	mu.11 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + focc_agro*FE_OC_AGRO + focc_bar*FE_OC_BAR + focc_trad*FE_OC_TRAD,
					base ~ dnorm(0,100),
					c(focc_agro, focc_bar, focc_trad) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, focc_agro=0, focc_bar=0, focc_trad=0),			
			warmup=1e3, iter=5e3, chains=1, cores=4
		)	
	post	<- extract.samples(mu.11)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_FOCC_OTH=logistic( post$base ),
			PI_FOCC_AGRO=logistic( post$base+post$focc_agro ),
			PI_FOCC_BAR=logistic( post$base+post$focc_bar ),
			PI_FOCC_TRAD=logistic( post$base+post$focc_trad ),
			OR_FOCC_OTH=exp( post$base ),
			OR_FOCC_AGRO=exp( post$base+post$focc_agro ),
			OR_FOCC_BAR=exp( post$base+post$focc_bar ),
			OR_FOCC_TRAD=exp( post$base+post$focc_trad ),
			ORX_FOCC_AGRO=exp( post$focc_agro ),
			ORX_FOCC_BAR=exp( post$focc_bar ),
			ORX_FOCC_TRAD=exp( post$focc_trad )	)									
	dp		<- melt(dp, id.vars='MC')	
	ds.11	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.11	<- dcast.data.table(ds.11, variable~STAT, value.var='V')
	ds.11[, MODEL:= 'univariate female occup']	
	#
	#	univariate female to male ratio
	#
	mu.12 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + prevratio_fm * RFM_ABOVE_MEAN, 
					base ~ dnorm(0,100),
					c(prevratio_fm) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, prevratio_fm=0),			
			warmup=1e3, iter=5e3, chains=1, cores=4
		)
	post	<- extract.samples(mu.12)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_PREVFMRATIO_ABOVE=logistic( post$base+post$prevratio_fm ),
			PI_PREVFMRATIO_BELOW=logistic( post$base ),
			OR_PREVFMRATIO_ABOVE=exp( post$base+post$prevratio_fm ),
			OR_PREVFMRATIO_BELOW=exp( post$base ),
			ORX_PREVFMRATIO_ABOVE=exp( post$prevratio_fm )
			)									
	dp		<- melt(dp, id.vars='MC')	
	ds.12	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.12	<- dcast.data.table(ds.12, variable~STAT, value.var='V')
	ds.12[, MODEL:= 'univariate fm prev ratio']	
	
	#	
	#	collect results
	ds		<- rbind(ds.1,ds.2,ds.3,ds.4,ds.5,ds.6,ds.7,ds.8,ds.9,ds.10,ds.11,ds.12)
	ds[, LABEL:= paste0(round(MED, d=2), ' [', round(CL, d=2),'-', round( CU, d=2),']')]
	ds[, STAT:= factor(gsub('^([^_]+)_.*','\\1',variable), levels=c('PI','OR','ORX'), labels=c('proportion MF','odds MF','odds ratio'))]
	ds[, FACTOR:=gsub('^([^_]+)_([^_]+)_([^_]+)','\\3',variable)]
	ds[, MOFA:=gsub('^([^_]+)_([^_]+)_([^_]+)','\\2-\\3',variable)]
	ds[, MODELTYPE:= factor(grepl('univariate',MODEL), levels=c(TRUE,FALSE),labels=c('univariate','multivariate'))]
	#
	save(ds, file=paste0(outfile.base,'_propmf_factors_univariate.rda'))	
	# 
	ds		<- dcast.data.table(ds, FACTOR+MOFA~MODELTYPE+STAT, value.var='LABEL')	
	set(ds, NULL, 'MOFA', ds[, factor(MOFA, levels=c(
									"PREVFMRATIO-BELOW","PREVFMRATIO-ABOVE",
									"COMM-FISH","COMM-AGRO","COMM-TRAD","COMM-MXD",
									"COUPLE-NO","COUPLE-YES",
									"SAMEHH-NO","SAMEHH-YES",
									"FEDU-YES","FEDU-NONE",
									"MEDU-YES","MEDU-NONE",      
									"FESEXP1YR-ONE","FESEXP1YR-G1",  
									"MASEXP1YR-ONE","MASEXP1YR-G1", 
									"CIRC-NO","CIRC-YES",                   
									"FOCC-OTH","FOCC-AGRO","FOCC-BAR","FOCC-TRAD",
									"MOCC-OTH","MOCC-AGRO","MOCC-FISH","MOCC-TRAD"   
							))])
	setkey(ds, MOFA)	
	for(x in c('univariate_odds ratio','multivariate_proportion MF','multivariate_odds MF','multivariate_odds ratio'))
		set(ds, which(is.na(ds[[x]])), x, '-')
	write.csv(ds, row.names=FALSE, file=paste0(outfile.base,'_propmf_factors_univariate_withprevration.csv'))
}

RakaiFull.gender.171122.couples.multivariateadditivemodels.nocouplesinmodel.stan.with.threshold<- function()
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
	
	#zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	#zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	load(infile)	
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	rtpdm[, PAIR_COMM_TYPE:= FEMALE_COMM_TYPE]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_TYPE!=MALE_COMM_TYPE)], 'PAIR_COMM_TYPE', 'mixed')
	set(rtpdm, rtpdm[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(rtpdm, rtpdm[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')
	
	#rtpdm[, table(PAIR_COMM_TYPE)]
	rtpdm[, COUPLE2:= factor(COUPLE=='no couple', levels=c(TRUE,FALSE), labels=c('no couple','couple'))]
	rtpdm[, SAMEHH:= factor(FEMALE_HH_NUM==MALE_HH_NUM, levels=c(TRUE,FALSE), labels=c('same hh','different hh'))]	
	rtpdm[, PAIR_COMM:= MALE_COMM_NUM_A]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_NUM_A!=MALE_COMM_NUM_A)], 'PAIR_COMM', 'mixed')
	rtpdm[, MALE_SEXP1OUT2:= factor(MALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(MALE_SEXP1OUT=='Unknown')], 'MALE_SEXP1OUT2', 'Unknown')
	rtpdm[, FEMALE_SEXP1OUT2:= factor(FEMALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(FEMALE_SEXP1OUT=='Unknown')], 'FEMALE_SEXP1OUT2', 'Unknown')
	
	
	df		<- subset(rtpdm, select=c(	MALE_RID, FEMALE_RID, PTY_RUN, IDCLU, LINK_MF, POSTERIOR_SCORE_MF, COUPLE2, SAMEHH, PAIR_COMM_TYPE, 
					MALE_RECENTVL, MALE_SEXP1YR, MALE_SEXP1OUT2, MALE_OCCUP_OLLI, MALE_OCAT, MALE_EDUCAT, MALE_CIRCUM,
					FEMALE_RECENTVL, FEMALE_SEXP1YR, FEMALE_SEXP1OUT2, FEMALE_OCCUP_OLLI, FEMALE_OCAT, FEMALE_EDUCAT 
			))
	#	missing data: fill in
	set(df, df[, which(is.na(MALE_RECENTVL))], 'MALE_RECENTVL', -1)
	set(df, df[, which(is.na(FEMALE_RECENTVL))], 'FEMALE_RECENTVL', -1)
	set(df, df[, which(is.na(MALE_CIRCUM))], 'MALE_CIRCUM', 'Unknown')
	set(df, df[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(df, df[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')
	#for(x in colnames(df)) print( c(x, any(is.na(df[[x]]))) )
	#	prepare data for STAN
	df[, COUPLE3:= as.integer(COUPLE2=='couple')]
	df[, SAMEHH3:= as.integer(SAMEHH=='same hh')]
	df[, COMM_FISH:= as.integer(PAIR_COMM_TYPE=='fisherfolk')]
	df[, COMM_AGR:= as.integer(PAIR_COMM_TYPE=='agrarian')]
	df[, COMM_TRAD:= as.integer(PAIR_COMM_TYPE=='trading')]
	df[, COMM_MXD:= as.integer(PAIR_COMM_TYPE=='mixed')]
	df[, FE_NOEDU:= as.integer(FEMALE_EDUCAT=='None')]
	df[, FE_NOEDU_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, MA_NOEDU:= as.integer(MALE_EDUCAT=='None')]
	df[, MA_NOEDU_MISS:= as.integer(MALE_EDUCAT=='Unknown')]	
	df[, MA_CIRCUM:= as.integer(MALE_CIRCUM=='Y')]
	df[, MA_CIRCUM_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, FE_SEXP1YR_G1:= as.integer(FEMALE_SEXP1YR!='1')]
	df[, MA_SEXP1YR_G1:= as.integer(MALE_SEXP1YR!='1')]
	df[, MA_OC_AGRO:= as.integer(MALE_OCAT=='Agro/House')]	
	df[, MA_OC_FISH:= as.integer(MALE_OCAT=='Fishing')]
	df[, MA_OC_TRAD:= as.integer(MALE_OCAT=='Trading/Shop keeper')]
	df[, MA_OC_OTH:= as.integer(MALE_OCAT%in%c('zother','Bar/waitress','Student','Boda/Trucking'))]	
	df[, FE_OC_AGRO:= as.integer(FEMALE_OCAT=='Agro/House')]
	df[, FE_OC_BAR:= as.integer(FEMALE_OCAT=='Bar/waitress')]			
	df[, FE_OC_TRAD:= as.integer(FEMALE_OCAT=='Trading/Shop keeper')]
	df[, FE_OC_OTH:= as.integer(FEMALE_OCAT%in%c('zother','Student'))]		
	
	#
	#	STAN
	#
	
	#	first multivariate model
	mm.1 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + couple_b*COUPLE3 + samehh_b*SAMEHH3 +
							#comm_fish*COMM_FISH +	this is base 
							comm_agr*COMM_AGR + comm_trad*COMM_TRAD + comm_mxd*COMM_MXD +
							noedu_f*FE_NOEDU + noedu_m*MA_NOEDU +
							sexp_f*FE_SEXP1YR_G1 + sexp_m*MA_SEXP1YR_G1 +
							circum_m*MA_CIRCUM,
					base ~ dnorm(0,100),
					c(couple_b, samehh_b) ~ dnorm(0,10),
					c(comm_agr, comm_trad, comm_mxd) ~ dnorm(0,10),
					c(noedu_f, noedu_m) ~ dnorm(0,10),								
					c(sexp_f, sexp_m, circum_m) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(base=0, couple_b=0, samehh_b=0, comm_agr=0, comm_trad=0, comm_mxd=0, noedu_f=0, noedu_m=0, sexp_f=0, sexp_m=0, circum_m=0),			
			warmup=1e3, iter=5e3, chains=1, cores=4
	)
	plot(mm.1)
	precis(mm.1, prob=0.95)
	pairs(mm.1)
	plot(precis(mm.1, prob=0.95))
	#	same household definitely stronger effect than couples
	mm.2 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + couple_b*COUPLE3 + samehh_b*SAMEHH3 +
							# comm_fish*COMM_FISH +	this is base 
							comm_agr*COMM_AGR + comm_trad*COMM_TRAD + comm_mxd*COMM_MXD +
							noedu_f*FE_NOEDU + noedu_m*MA_NOEDU +
							sexp_f*FE_SEXP1YR_G1 + sexp_m*MA_SEXP1YR_G1 +
							circum_m*MA_CIRCUM +										
							# other occupation is baseline
							mocc_agro*MA_OC_AGRO + mocc_fish*MA_OC_FISH + mocc_trad*MA_OC_TRAD +										
							focc_agro*FE_OC_AGRO + focc_bar*FE_OC_BAR +  focc_trad*FE_OC_TRAD,
					base ~ dnorm(0,100),
					c(couple_b, samehh_b) ~ dnorm(0,10),
					c(comm_agr, comm_trad, comm_mxd) ~ dnorm(0,10),
					c(noedu_f, noedu_m) ~ dnorm(0,10),								
					c(sexp_f, sexp_m, circum_m) ~ dnorm(0,10),
					c(mocc_agro, mocc_fish, mocc_trad) ~ dnorm(0,10),
					c(focc_agro, focc_bar, focc_trad) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, couple_b=0, samehh_b=0, comm_agr=0, comm_trad=0, comm_mxd=0, noedu_f=0, noedu_m=0, sexp_f=0, sexp_m=0, circum_m=0, 
					focc_agro=0, focc_bar=0, focc_trad=0,
					mocc_agro=0, mocc_fish=0, mocc_trad=0),			
			warmup=1e3, iter=5e3, chains=1, cores=4
	)
	mm.2b <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + samehh_b*SAMEHH3 +
							# comm_fish*COMM_FISH +	this is base 
							comm_agr*COMM_AGR + comm_trad*COMM_TRAD + comm_mxd*COMM_MXD +
							noedu_f*FE_NOEDU + noedu_m*MA_NOEDU +
							sexp_f*FE_SEXP1YR_G1 + sexp_m*MA_SEXP1YR_G1 +
							circum_m*MA_CIRCUM +										
							# other occupation is baseline
							mocc_agro*MA_OC_AGRO + mocc_fish*MA_OC_FISH + mocc_trad*MA_OC_TRAD +										
							focc_agro*FE_OC_AGRO + focc_bar*FE_OC_BAR +  focc_trad*FE_OC_TRAD,
					base ~ dnorm(0,100),
					c(samehh_b) ~ dnorm(0,10),
					c(comm_agr, comm_trad, comm_mxd) ~ dnorm(0,10),
					c(noedu_f, noedu_m) ~ dnorm(0,10),								
					c(sexp_f, sexp_m, circum_m) ~ dnorm(0,10),
					c(mocc_agro, mocc_fish, mocc_trad) ~ dnorm(0,10),
					c(focc_agro, focc_bar, focc_trad) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, samehh_b=0, comm_agr=0, comm_trad=0, comm_mxd=0, noedu_f=0, noedu_m=0, sexp_f=0, sexp_m=0, circum_m=0, 
					focc_agro=0, focc_bar=0, focc_trad=0,
					mocc_agro=0, mocc_fish=0, mocc_trad=0),			
			warmup=1e3, iter=5e3, chains=1, cores=4
	)		
	#plot(precis(mm.2, prob=0.95))
	#pairs(mm.2)
	post	<- extract.samples(mm.2b)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_COMM_FISH=logistic( post$base ),				
			PI_COMM_AGRO=logistic( post$base+post$comm_agr ), 											
			PI_COMM_TRAD=logistic( post$base+post$comm_trad ),							
			PI_COMM_MXD=logistic( post$base+post$comm_mxd ),
			#PI_COUPLE_YES=logistic( post$base+post$couple_b ),
			PI_SAMEHH_YES=logistic( post$base+post$samehh_b ),
			PI_FEDU_NONE=logistic( post$base+post$noedu_f ),
			PI_MEDU_NONE=logistic( post$base+post$noedu_m ),							
			PI_FESEXP1YR_G1=logistic( post$base+post$sexp_f ),
			PI_MASEXP1YR_G1=logistic( post$base+post$sexp_m ),							
			PI_CIRC_YES=logistic( post$base+post$circum_m ),							
			PI_MOCC_AGRO=logistic( post$base+post$mocc_agro ),
			PI_MOCC_FISH=logistic( post$base+post$mocc_fish ),
			PI_MOCC_TRAD=logistic( post$base+post$mocc_trad ),
			PI_FOCC_AGRO=logistic( post$base+post$focc_agro ),
			PI_FOCC_BAR=logistic( post$base+post$focc_bar ),
			PI_FOCC_TRAD=logistic( post$base+post$focc_trad ),							
			OR_COMM_FISH=exp( post$base ),				
			OR_COMM_AGRO=exp( post$base+post$comm_agr ), 											
			OR_COMM_TRAD=exp( post$base+post$comm_trad ),							
			OR_COMM_MXD=exp( post$base+post$comm_mxd ),
			#OR_COUPLE_YES=exp( post$base+post$couple_b ),
			OR_SAMEHH_YES=exp( post$base+post$samehh_b ),
			OR_FEDU_NONE=exp( post$base+post$noedu_f ),
			OR_MEDU_NONE=exp( post$base+post$noedu_m ),							
			OR_FESEXP1YR_G1=exp( post$base+post$sexp_f ),
			OR_MASEXP1YR_G1=exp( post$base+post$sexp_m ),
			OR_CIRC_YES=exp( post$base+post$circum_m ),							
			OR_MOCC_AGRO=exp( post$base+post$mocc_agro ),
			OR_MOCC_FISH=exp( post$base+post$mocc_fish ),
			OR_MOCC_TRAD=exp( post$base+post$mocc_trad ),
			OR_FOCC_AGRO=exp( post$base+post$focc_agro ),
			OR_FOCC_BAR=exp( post$base+post$focc_bar ),
			OR_FOCC_TRAD=exp( post$base+post$focc_trad ),											
			ORX_COMM_AGRO=exp( post$comm_agr ), 											
			ORX_COMM_TRAD=exp( post$comm_trad ),							
			ORX_COMM_MXD=exp( post$comm_mxd ),
			#ORX_COUPLE_YES=exp( post$couple_b ),
			ORX_SAMEHH_YES=exp( post$samehh_b ),
			ORX_FEDU_NONE=exp( post$noedu_f ),
			ORX_MEDU_NONE=exp( post$noedu_m ),
			ORX_FESEXP1YR_G1=exp( post$sexp_f ),
			ORX_MASEXP1YR_G1=exp( post$sexp_m ),							
			ORX_CIRC_YES=exp( post$circum_m ),							
			ORX_MOCC_AGRO=exp( post$mocc_agro ),
			ORX_MOCC_FISH=exp( post$mocc_fish ),
			ORX_MOCC_TRAD=exp( post$mocc_trad ),
			ORX_FOCC_AGRO=exp( post$focc_agro ),
			ORX_FOCC_BAR=exp( post$focc_bar ),
			ORX_FOCC_TRAD=exp( post$focc_trad )								
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.2	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.2	<- dcast.data.table(ds.2, variable~STAT, value.var='V')
	ds.2[, MODEL:= 'multivariate mm2']	
	#
	#	univariate couples
	#
	mu.1 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + couple_b*COUPLE3,
					base ~ dnorm(0,100),
					couple_b ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, couple_b=0),			
			warmup=500, iter=2e3, chains=1, cores=4
	)				
	post	<- extract.samples(mu.1)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_COUPLE_NO=logistic( post$base ),
			PI_COUPLE_YES=logistic( post$base+post$couple_b ),
			OR_COUPLE_NO=exp( post$base ),
			OR_COUPLE_YES=exp( post$base+post$couple_b ),
			ORX_COUPLE_YES=exp( post$couple_b )							
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.1	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.1	<- dcast.data.table(ds.1, variable~STAT, value.var='V')
	ds.1[, MODEL:= 'univariate couples']	
	#
	#	univariate household
	#
	mu.3 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + samehh_b*SAMEHH3,
					base ~ dnorm(0,100),
					samehh_b ~ dnorm(0,10)					
			),
			data=as.data.frame(df), 
			start=list(	base=0, samehh_b=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.3)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_SAMEHH_NO=logistic( post$base ),
			PI_SAMEHH_YES=logistic( post$base+post$samehh_b ),
			OR_SAMEHH_NO=exp( post$base ),
			OR_SAMEHH_YES=exp( post$base+post$samehh_b ),
			ORX_SAMEHH_YES=exp( post$samehh_b )														
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.3	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.3	<- dcast.data.table(ds.3, variable~STAT, value.var='V')
	ds.3[, MODEL:= 'univariate same household']	
	#
	#	univariate community type
	#
	mu.4 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + comm_agr*COMM_AGR + comm_trad*COMM_TRAD + comm_mxd*COMM_MXD,
					base ~ dnorm(0,100),							
					c(comm_agr, comm_trad, comm_mxd) ~ dnorm(0,10)							
			),
			data=as.data.frame(df), 
			start=list(	base=0, comm_agr=0, comm_trad=0, comm_mxd=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.4)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_COMM_FISH=logistic( post$base ),				
			PI_COMM_AGRO=logistic( post$base+post$comm_agr ), 											
			PI_COMM_TRAD=logistic( post$base+post$comm_trad ),							
			PI_COMM_MXD=logistic( post$base+post$comm_mxd ),
			OR_COMM_FISH=exp( post$base ),				
			OR_COMM_AGRO=exp( post$base+post$comm_agr ), 											
			OR_COMM_TRAD=exp( post$base+post$comm_trad ),							
			OR_COMM_MXD=exp( post$base+post$comm_mxd ),				
			ORX_COMM_AGRO=exp( post$comm_agr ), 											
			ORX_COMM_TRAD=exp( post$comm_trad ),							
			ORX_COMM_MXD=exp( post$comm_mxd )							
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.4	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.4	<- dcast.data.table(ds.4, variable~STAT, value.var='V')
	ds.4[, MODEL:= 'univariate commtype']	
	#
	#	univariate female education
	#
	mu.5 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + noedu_f*FE_NOEDU,
					base ~ dnorm(0,100),
					noedu_f ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, noedu_f=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.5)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_FEDU_YES=logistic( post$base ),
			PI_FEDU_NONE=logistic( post$base+post$noedu_f ),
			OR_FEDU_YES=exp( post$base ),
			OR_FEDU_NONE=exp( post$base+post$noedu_f ),
			ORX_FEDU_NONE=exp( post$noedu_f)	)
	dp		<- melt(dp, id.vars='MC')	
	ds.5	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.5	<- dcast.data.table(ds.5, variable~STAT, value.var='V')
	ds.5[, MODEL:= 'univariate female education']	
	#
	#	univariate male education
	#
	mu.6 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + noedu_m*MA_NOEDU,
					base ~ dnorm(0,100),
					noedu_m ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, noedu_m=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.6)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_MEDU_YES=logistic( post$base ),
			PI_MEDU_NONE=logistic( post$base+post$noedu_m ),
			OR_MEDU_YES=exp( post$base ),
			OR_MEDU_NONE=exp( post$base+post$noedu_m ),
			ORX_MEDU_NONE=exp( post$noedu_m)	)
	dp		<- melt(dp, id.vars='MC')	
	ds.6	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.6	<- dcast.data.table(ds.6, variable~STAT, value.var='V')
	ds.6[, MODEL:= 'univariate male education']	
	#
	#	univariate female sex partners
	#
	mu.7 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + sexp_f*FE_SEXP1YR_G1 ,
					base ~ dnorm(0,100),
					sexp_f ~ dnorm(0,10)					
			),
			data=as.data.frame(df), 
			start=list(	base=0, sexp_f=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.7)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_FESEXP1YR_ONE=logistic( post$base ),
			PI_FESEXP1YR_G1=logistic( post$base+post$sexp_f ),
			OR_FESEXP1YR_ONE=exp( post$base ),
			OR_FESEXP1YR_G1=exp( post$base+post$sexp_f ),
			ORX_FESEXP1YR_G1=exp( post$sexp_f )	)
	dp		<- melt(dp, id.vars='MC')	
	ds.7	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.7	<- dcast.data.table(ds.7, variable~STAT, value.var='V')
	ds.7[, MODEL:= 'univariate female sex partners']	
	#
	#	univariate male sex partners
	#
	mu.8 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + sexp_m*MA_SEXP1YR_G1 ,
					base ~ dnorm(0,100),
					sexp_m ~ dnorm(0,10)					
			),
			data=as.data.frame(df), 
			start=list(	base=0, sexp_m=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.8)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_MASEXP1YR_ONE=logistic( post$base ),
			PI_MASEXP1YR_G1=logistic( post$base+post$sexp_m ),
			OR_MASEXP1YR_ONE=exp( post$base ),
			OR_MASEXP1YR_G1=exp( post$base+post$sexp_m ),
			ORX_MASEXP1YR_G1=exp( post$sexp_m )	)
	dp		<- melt(dp, id.vars='MC')	
	ds.8	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.8	<- dcast.data.table(ds.8, variable~STAT, value.var='V')
	ds.8[, MODEL:= 'univariate male sex partners']	
	#
	#	univariate circumcision
	#
	mu.9 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + circum_m*MA_CIRCUM,
					base ~ dnorm(0,100),
					circum_m ~ dnorm(0,10)					
			),
			data=as.data.frame(df), 
			start=list(	base=0, circum_m=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.9)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_CIRC_NO=logistic( post$base ),
			PI_CIRC_YES=logistic( post$base+post$circum_m ),							
			OR_CIRC_NO=exp( post$base ),
			OR_CIRC_YES=exp( post$base+post$circum_m ),							
			ORX_CIRC_YES=exp( post$circum_m )
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.9	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.9	<- dcast.data.table(ds.9, variable~STAT, value.var='V')
	ds.9[, MODEL:= 'univariate circumcision']	
	#
	#	univariate male occupation
	#
	mu.10 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + mocc_agro*MA_OC_AGRO + mocc_fish*MA_OC_FISH + mocc_trad*MA_OC_TRAD,
					base ~ dnorm(0,100),
					c(mocc_agro, mocc_fish, mocc_trad) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, mocc_agro=0, mocc_fish=0, mocc_trad=0),			
			warmup=1e3, iter=5e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.10)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_MOCC_OTH=logistic( post$base ),
			PI_MOCC_AGRO=logistic( post$base+post$mocc_agro ),
			PI_MOCC_FISH=logistic( post$base+post$mocc_fish ),
			PI_MOCC_TRAD=logistic( post$base+post$mocc_trad ),
			OR_MOCC_OTH=exp( post$base ),
			OR_MOCC_AGRO=exp( post$base+post$mocc_agro ),
			OR_MOCC_FISH=exp( post$base+post$mocc_fish ),
			OR_MOCC_TRAD=exp( post$base+post$mocc_trad ),
			ORX_MOCC_AGRO=exp( post$mocc_agro ),
			ORX_MOCC_FISH=exp( post$mocc_fish ),
			ORX_MOCC_TRAD=exp( post$mocc_trad )	)									
	dp		<- melt(dp, id.vars='MC')	
	ds.10	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.10	<- dcast.data.table(ds.10, variable~STAT, value.var='V')
	ds.10[, MODEL:= 'univariate male occup']	
	#
	#	univariate female occupation
	#
	mu.11 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + focc_agro*FE_OC_AGRO + focc_bar*FE_OC_BAR + focc_trad*FE_OC_TRAD,
					base ~ dnorm(0,100),
					c(focc_agro, focc_bar, focc_trad) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, focc_agro=0, focc_bar=0, focc_trad=0),			
			warmup=1e3, iter=5e3, chains=1, cores=4
	)	
	post	<- extract.samples(mu.11)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_FOCC_OTH=logistic( post$base ),
			PI_FOCC_AGRO=logistic( post$base+post$focc_agro ),
			PI_FOCC_BAR=logistic( post$base+post$focc_bar ),
			PI_FOCC_TRAD=logistic( post$base+post$focc_trad ),
			OR_FOCC_OTH=exp( post$base ),
			OR_FOCC_AGRO=exp( post$base+post$focc_agro ),
			OR_FOCC_BAR=exp( post$base+post$focc_bar ),
			OR_FOCC_TRAD=exp( post$base+post$focc_trad ),
			ORX_FOCC_AGRO=exp( post$focc_agro ),
			ORX_FOCC_BAR=exp( post$focc_bar ),
			ORX_FOCC_TRAD=exp( post$focc_trad )	)									
	dp		<- melt(dp, id.vars='MC')	
	ds.11	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.11	<- dcast.data.table(ds.11, variable~STAT, value.var='V')
	ds.11[, MODEL:= 'univariate female occup']	
	
	#
	#	collect results
	ds		<- rbind(ds.1,ds.2,ds.3,ds.4,ds.5,ds.6,ds.7,ds.8,ds.9,ds.10,ds.11)
	ds[, LABEL:= paste0(round(MED, d=2), '\n[', round(CL, d=2),'-', round( CU, d=2),']')]
	ds[, STAT:= factor(gsub('^([^_]+)_.*','\\1',variable), levels=c('PI','OR','ORX'), labels=c('proportion MF','odds MF','odds ratio'))]
	ds[, FACTOR:=gsub('^([^_]+)_([^_]+)_([^_]+)','\\3',variable)]
	ds[, MOFA:=gsub('^([^_]+)_([^_]+)_([^_]+)','\\2-\\3',variable)]
	ds[, MODELTYPE:= factor(grepl('univariate',MODEL), levels=c(TRUE,FALSE),labels=c('univariate','multivariate'))]
	#
	save(ds, file=paste0(outfile.base,'_propmf_factors_univariate.rda'))	
	# 
	ds		<- dcast.data.table(ds, FACTOR+MOFA~MODELTYPE+STAT, value.var='LABEL')	
	set(ds, NULL, 'MOFA', ds[, factor(MOFA, levels=c(
									"COMM-FISH","COMM-AGRO","COMM-TRAD","COMM-MXD",
									"COUPLE-NO","COUPLE-YES",
									"SAMEHH-NO","SAMEHH-YES",
									"FEDU-YES","FEDU-NONE",
									"MEDU-YES","MEDU-NONE",      
									"FESEXP1YR-ONE","FESEXP1YR-G1",  
									"MASEXP1YR-ONE","MASEXP1YR-G1", 
									"CIRC-NO","CIRC-YES",                   
									"FOCC-OTH","FOCC-AGRO","FOCC-BAR","FOCC-TRAD",
									"MOCC-OTH","MOCC-AGRO","MOCC-FISH","MOCC-TRAD"   
							))])
	setkey(ds, MOFA)	
	for(x in c('univariate_odds ratio','multivariate_proportion MF','multivariate_odds MF','multivariate_odds ratio'))
		set(ds, which(is.na(ds[[x]])), x, '-')
	write.csv(ds, row.names=FALSE, file=paste0(outfile.base,'_propmf_factors_univariate.csv'))
}

RakaiFull.gender.171122.couples.interactionmodels.stan.with.threshold<- function()
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
	
	#zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	#zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	load(infile)	
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm',SELECT))
	rtpdm[, PAIR_COMM_TYPE:= FEMALE_COMM_TYPE]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_TYPE!=MALE_COMM_TYPE)], 'PAIR_COMM_TYPE', 'mixed')
	set(rtpdm, rtpdm[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(rtpdm, rtpdm[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')
	
	#rtpdm[, table(PAIR_COMM_TYPE)]
	rtpdm[, COUPLE2:= factor(COUPLE=='no couple', levels=c(TRUE,FALSE), labels=c('no couple','couple'))]
	rtpdm[, SAMEHH:= factor(FEMALE_HH_NUM==MALE_HH_NUM, levels=c(TRUE,FALSE), labels=c('same hh','different hh'))]	
	rtpdm[, PAIR_COMM:= MALE_COMM_NUM_A]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_NUM_A!=MALE_COMM_NUM_A)], 'PAIR_COMM', 'mixed')
	rtpdm[, MALE_SEXP1OUT2:= factor(MALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(MALE_SEXP1OUT=='Unknown')], 'MALE_SEXP1OUT2', 'Unknown')
	rtpdm[, FEMALE_SEXP1OUT2:= factor(FEMALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(FEMALE_SEXP1OUT=='Unknown')], 'FEMALE_SEXP1OUT2', 'Unknown')
	
	
	#	superspreaders
	tmp		<- subset(rtpdm, grepl('mf',SELECT))[, list(N_REC=length(FEMALE_RID)), by='MALE_RID']
	subset(tmp, N_REC>1)	#15 male "superspreaders"
	tmp		<- subset(rtpdm, grepl('fm',SELECT))[, list(N_REC=length(MALE_RID)), by='FEMALE_RID']
	subset(tmp, N_REC>1)	#5 female "superspreaders"
	
	ans		<- rtpdm[, 	list(	ANA='overall m->f',
					TOTAL=length(MALE_RID),
					STAT=c('CENTRAL','L95','U95'),
					V=as.numeric( binconf( length(MALE_RID[grepl('mf',SELECT)]), length(MALE_RID) ) )
			), ]
	ans		<- dcast.data.table(ans, ANA+TOTAL~STAT, value.var='V')
	#	by community type
	tmp		<- rtpdm[, 	list(	ANA='by commtype m->f',
					TOTAL=length(MALE_RID),
					STAT=c('CENTRAL','L95','U95'),
					V=as.numeric( binconf( length(MALE_RID[grepl('mf',SELECT)]), length(MALE_RID) ) )	
			), 
			by=c('PAIR_COMM_TYPE')]
	tmp		<- dcast.data.table(tmp, ANA+TOTAL+PAIR_COMM_TYPE~STAT, value.var='V')
	ans		<- rbind(ans, tmp, fill=TRUE)
	#	by couple status
	tmp		<- rtpdm[, 	list(	ANA='by couple m->f',
					TOTAL=length(MALE_RID),
					STAT=c('CENTRAL','L95','U95'),
					V=as.numeric( binconf( length(MALE_RID[grepl('mf',SELECT)]), length(MALE_RID) ) )	
			), 
			by=c('COUPLE2')]
	tmp		<- dcast.data.table(tmp, ANA+TOTAL+COUPLE2~STAT, value.var='V')
	ans		<- rbind(ans, tmp, fill=TRUE)
	#	by couple and community status
	tmp		<- rtpdm[, 	list(	ANA='by commtype & couple m->f',
					TOTAL=length(MALE_RID),
					STAT=c('CENTRAL','L95','U95'),
					V=as.numeric( binconf( length(MALE_RID[grepl('mf',SELECT)]), length(MALE_RID) ) )	
			), 
			by=c('PAIR_COMM_TYPE','COUPLE2')]
	tmp		<- dcast.data.table(tmp, ANA+TOTAL+PAIR_COMM_TYPE+COUPLE2~STAT, value.var='V')
	ans		<- rbind(ans, tmp, fill=TRUE)
	#	by "same household" status
	tmp		<- rtpdm[, 	list(	ANA='by household m->f',
					TOTAL=length(MALE_RID),
					STAT=c('CENTRAL','L95','U95'),
					V=as.numeric( binconf( length(MALE_RID[grepl('mf',SELECT)]), length(MALE_RID) ) )	
			), 
			by=c('SAMEHH')]
	tmp		<- dcast.data.table(tmp, ANA+TOTAL+SAMEHH~STAT, value.var='V')
	ans		<- rbind(ans, tmp, fill=TRUE)
	#	by "same household" status and community status
	tmp		<- rtpdm[, 	list(	ANA='by household & couple m->f',
					TOTAL=length(MALE_RID),
					STAT=c('CENTRAL','L95','U95'),
					V=as.numeric( binconf( length(MALE_RID[grepl('mf',SELECT)]), length(MALE_RID) ) )	
			), 
			by=c('PAIR_COMM_TYPE','SAMEHH')]
	tmp		<- dcast.data.table(tmp, ANA+TOTAL+PAIR_COMM_TYPE+SAMEHH~STAT, value.var='V')
	ans		<- rbind(ans, tmp, fill=TRUE)
	#	by community --> sample sizes too small
	#
	# 	by education of female partner
	tmp		<- rtpdm[, 	list(	ANA='by education female m->f',
					TOTAL=length(MALE_RID),
					STAT=c('CENTRAL','L95','U95'),
					V=as.numeric( binconf( length(MALE_RID[grepl('mf',SELECT)]), length(MALE_RID) ) )	
			), 
			by=c('FEMALE_EDUCAT')]
	tmp		<- subset(tmp, !is.na(FEMALE_EDUCAT))
	tmp		<- dcast.data.table(tmp, ANA+TOTAL+FEMALE_EDUCAT~STAT, value.var='V')
	ans		<- rbind(ans, tmp, fill=TRUE)
	#
	# 	by education of male partner
	tmp		<- rtpdm[, 	list(	ANA='by education male m->f',
					TOTAL=length(MALE_RID),
					STAT=c('CENTRAL','L95','U95'),
					V=as.numeric( binconf( length(MALE_RID[grepl('mf',SELECT)]), length(MALE_RID) ) )	
			), 
			by=c('MALE_EDUCAT')]
	tmp		<- subset(tmp, MALE_EDUCAT!='Unknown')
	tmp		<- dcast.data.table(tmp, ANA+TOTAL+MALE_EDUCAT~STAT, value.var='V')
	ans		<- rbind(ans, tmp, fill=TRUE)	
	# 	by sex partners in 1yr male and couple
	tmp		<- rtpdm[, 	list(	ANA='by sex partners in 1yr male m->f',
					TOTAL=length(MALE_RID),
					STAT=c('CENTRAL','L95','U95'),
					V=as.numeric( binconf( length(MALE_RID[grepl('mf',SELECT)]), length(MALE_RID) ) )	
			), 
			by=c('MALE_SEXP1YR','COUPLE2')]	
	tmp		<- dcast.data.table(tmp, ANA+TOTAL+COUPLE2+MALE_SEXP1YR~STAT, value.var='V')
	ans		<- rbind(ans, tmp, fill=TRUE)
	# 	by sex partners out male and couple
	tmp		<- rtpdm[, 	list(	ANA='by sex partners out male m->f',
					TOTAL=length(MALE_RID),
					STAT=c('CENTRAL','L95','U95'),
					V=as.numeric( binconf( length(MALE_RID[grepl('mf',SELECT)]), length(MALE_RID) ) )	
			), 
			by=c('MALE_SEXP1OUT2','COUPLE2')]	
	tmp		<- dcast.data.table(tmp, ANA+TOTAL+COUPLE2+MALE_SEXP1OUT2~STAT, value.var='V')
	ans		<- rbind(ans, tmp, fill=TRUE)
	# 	by sex partners in 1yr female and couple
	tmp		<- rtpdm[, 	list(	ANA='by sex partners in 1yr female m->f',
					TOTAL=length(MALE_RID),
					STAT=c('CENTRAL','L95','U95'),
					V=as.numeric( binconf( length(MALE_RID[grepl('mf',SELECT)]), length(MALE_RID) ) )	
			), 
			by=c('FEMALE_SEXP1YR','COUPLE2')]
	tmp		<- subset(tmp, FEMALE_SEXP1YR!='Unknown')
	tmp		<- dcast.data.table(tmp, ANA+TOTAL+COUPLE2+FEMALE_SEXP1YR~STAT, value.var='V')
	ans		<- rbind(ans, tmp, fill=TRUE)
	#	by male circumcision status
	tmp		<- rtpdm[, 	list(	ANA='by circum couple m->f',
					TOTAL=length(MALE_RID),
					STAT=c('CENTRAL','L95','U95'),
					V=as.numeric( binconf( length(MALE_RID[grepl('mf',SELECT)]), length(MALE_RID) ) )	
			), 
			by=c('MALE_CIRCUM','COUPLE2')]	
	tmp		<- subset(tmp, !is.na(MALE_CIRCUM))
	tmp		<- dcast.data.table(tmp, ANA+TOTAL+COUPLE2+MALE_CIRCUM~STAT, value.var='V')
	ans		<- rbind(ans, tmp, fill=TRUE)
	#	by female primary occupation and couple
	tmp		<- rtpdm[, 	list(	ANA='by female occupation couple m->f',
					TOTAL=length(MALE_RID),
					STAT=c('CENTRAL','L95','U95'),
					V=as.numeric( binconf( length(MALE_RID[grepl('mf',SELECT)]), length(MALE_RID) ) )	
			), 
			by=c('FEMALE_OCCUP_OLLI','COUPLE2')]		
	tmp		<- dcast.data.table(tmp, ANA+TOTAL+COUPLE2+FEMALE_OCCUP_OLLI~STAT, value.var='V')
	ans		<- rbind(ans, tmp, fill=TRUE)
	#	by male primary occupation (olli) and couple
	tmp		<- rtpdm[, 	list(	ANA='by male occupation couple m->f',
					TOTAL=length(MALE_RID),
					STAT=c('CENTRAL','L95','U95'),
					V=as.numeric( binconf( length(MALE_RID[grepl('mf',SELECT)]), length(MALE_RID) ) )	
			), 
			by=c('MALE_OCCUP_OLLI','COUPLE2')]		
	tmp		<- dcast.data.table(tmp, ANA+TOTAL+COUPLE2+MALE_OCCUP_OLLI~STAT, value.var='V')
	ans		<- rbind(ans, tmp, fill=TRUE)
	#	by male primary occupation (Kate) and couple
	tmp		<- rtpdm[, 	list(	ANA='by male occupation (Kate) couple m->f',
					TOTAL=length(MALE_RID),
					STAT=c('CENTRAL','L95','U95'),
					V=as.numeric( binconf( length(MALE_RID[grepl('mf',SELECT)]), length(MALE_RID) ) )	
			), 
			by=c('MALE_OCAT','COUPLE2')]		
	tmp		<- dcast.data.table(tmp, ANA+TOTAL+COUPLE2+MALE_OCAT~STAT, value.var='V')
	ans		<- rbind(ans, tmp, fill=TRUE)
	setkey(ans, ANA, COUPLE2)
	
	#	TODO gradient by MALE_RECENTVL ?
	
	
	require(rethinking)
	
	df		<- subset(rtpdm, select=c(	MALE_RID, FEMALE_RID, PTY_RUN, IDCLU, LINK_MF, POSTERIOR_SCORE_MF, COUPLE2, SAMEHH, PAIR_COMM_TYPE, 
					MALE_RECENTVL, MALE_SEXP1YR, MALE_SEXP1OUT2, MALE_OCCUP_OLLI, MALE_OCAT, MALE_EDUCAT, MALE_CIRCUM,
					FEMALE_RECENTVL, FEMALE_SEXP1YR, FEMALE_SEXP1OUT2, FEMALE_OCCUP_OLLI, FEMALE_OCAT, FEMALE_EDUCAT 
			))
	#	missing data: fill in
	set(df, df[, which(is.na(MALE_RECENTVL))], 'MALE_RECENTVL', -1)
	set(df, df[, which(is.na(FEMALE_RECENTVL))], 'FEMALE_RECENTVL', -1)
	set(df, df[, which(is.na(MALE_CIRCUM))], 'MALE_CIRCUM', 'Unknown')
	set(df, df[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')
	set(df, df[, which(is.na(MALE_EDUCAT))], 'MALE_EDUCAT', 'Unknown')
	#for(x in colnames(df)) print( c(x, any(is.na(df[[x]]))) )
	#	prepare data for STAN
	df[, COUPLE3:= as.integer(COUPLE2=='couple')]
	df[, SAMEHH3:= as.integer(SAMEHH=='same hh')]
	df[, COMM_FISH:= as.integer(PAIR_COMM_TYPE=='fisherfolk')]
	df[, COMM_AGR:= as.integer(PAIR_COMM_TYPE=='agrarian')]
	df[, COMM_TRAD:= as.integer(PAIR_COMM_TYPE=='trading')]
	df[, COMM_MXD:= as.integer(PAIR_COMM_TYPE=='mixed')]
	df[, FE_NOEDU:= as.integer(FEMALE_EDUCAT=='None')]
	df[, FE_NOEDU_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, MA_NOEDU:= as.integer(MALE_EDUCAT=='None')]
	df[, MA_NOEDU_MISS:= as.integer(MALE_EDUCAT=='Unknown')]	
	df[, MA_CIRCUM:= as.integer(MALE_CIRCUM=='Y')]
	df[, MA_CIRCUM_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, FE_SEXP1YR_G1:= as.integer(FEMALE_SEXP1YR!='1')]
	df[, MA_SEXP1YR_G1:= as.integer(MALE_SEXP1YR!='1')]
	df[, MA_OC_AGRO:= as.integer(MALE_OCAT=='Agro/House')]
	df[, MA_OC_BAR:= as.integer(MALE_OCAT=='Bar/waitress')]	
	df[, MA_OC_BODA:= as.integer(MALE_OCAT=='Boda/Trucking')]
	df[, MA_OC_FISH:= as.integer(MALE_OCAT=='Fishing')]
	df[, MA_OC_STUD:= as.integer(MALE_OCAT=='Student')]
	df[, MA_OC_TRAD:= as.integer(MALE_OCAT=='Trading/Shop keeper')]
	df[, MA_OC_OTH:= as.integer(MALE_OCAT=='zother')]	
	df[, FE_OC_AGRO:= as.integer(FEMALE_OCAT=='Agro/House')]
	df[, FE_OC_BAR:= as.integer(FEMALE_OCAT=='Bar/waitress')]			
	df[, FE_OC_STUD:= as.integer(FEMALE_OCAT=='Student')]
	df[, FE_OC_TRAD:= as.integer(FEMALE_OCAT=='Trading/Shop keeper')]
	df[, FE_OC_OTH:= as.integer(FEMALE_OCAT=='zother')]		
	#df[, MA_OC_CAS:= as.integer(MALE_OCCUP_OLLI=='Casual laborer/unemployed')]
	#df[, MA_OC_CONSTR:= as.integer(MALE_OCCUP_OLLI=='Construction/Mechanic')]	
	#df[, MA_OC_GOV:= as.integer(MALE_OCCUP_OLLI=='Government/clerical/related')]
	
	
	#
	#	GLMS
	#
	
	mf.1	<- glm(data=df, LINK_MF~COUPLE2, family=binomial)	
	summary(mf.1)
	logistic(  c(	'no couple'=sum(coef(mf.1)[c('(Intercept)')]),
					'couple'=sum(coef(mf.1)[c('(Intercept)','COUPLE2couple')])					
			))
	mf.3	<- glm(data=df, LINK_MF~COUPLE2:PAIR_COMM_TYPE, family=binomial)
	summary(mf.3)	#no converge
	#
	#	STAN
	#
	
	#	couples effect
	mg.1 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + couple_b*COUPLE3, 
					base ~ dnorm(0,100),
					couple_b ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(base=logit(0.5), couple_b=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)
	post	<- extract.samples(mg.1)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_OVERALL_NOCOUPLE=logistic( post$base ),
			OR_OVERALL_NOCOUPLE=exp( post$base ),
			PI_OVERALL_COUPLE=logistic( post$base+post$couple_b ),
			OR_OVERALL_COUPLE=exp( post$base+post$couple_b ),
			ORC_OVERALL_COUPLE= exp(post$couple_b))
	dp		<- melt(dp, id.vars='MC')	
	ds.1	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.1	<- dcast.data.table(ds.1, variable~STAT, value.var='V')
	ds.1[, MODEL:= 'couple overall']
	#	household effect
	mg.2 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + samehh_b*SAMEHH3, 
					base ~ dnorm(0,100),
					samehh_b ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(base=logit(0.5), samehh_b=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)
	post	<- extract.samples(mg.2)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_OVERALL_DIFFHH=logistic( post$base ),
			OR_OVERALL_DIFFHH=exp( post$base ),
			PI_OVERALL_SAMEHH=logistic( post$base+post$samehh_b ),
			OR_OVERALL_SAMEHH=exp( post$base+post$samehh_b ),
			ORC_OVERALL_SAMEHH= exp(post$samehh_b))
	dp		<- melt(dp, id.vars='MC')	
	ds.2	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.2	<- dcast.data.table(ds.2, variable~STAT, value.var='V')
	ds.2[, MODEL:= 'same household']
	#	couples community interaction model
	mg.3 	<- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- 	fish*COMM_FISH + fish_couple*COMM_FISH*COUPLE3 +
							agr*COMM_AGR + agr_couple*COMM_AGR*COUPLE3 +
							trad*COMM_TRAD + trad_couple*COMM_TRAD*COUPLE3 +
							cmxd*COMM_MXD + cmxd_couple*COMM_MXD*COUPLE3, 
					c(fish,agr,trad,cmxd) ~ dnorm(0,100),
					c(fish_couple,agr_couple,trad_couple,cmxd_couple) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(fish=0,agr=0,trad=0,cmxd=0,fish_couple=0,agr_couple=0,trad_couple=0,cmxd_couple=0),			
			warmup=5e2, iter=4e3, chains=1, cores=4
	)
	post	<- extract.samples(mg.3)
	dp		<- data.table(	MC= seq_along(post$fish),
			PI_FISH_NOCOUPLE= logistic( post$fish ),
			OR_FISH_NOCOUPLE= exp(post$fish),
			PI_FISH_COUPLE=logistic( post$fish+post$fish_couple ),	
			OR_FISH_COUPLE= exp(post$fish+post$fish_couple),							
			PI_AGR_NOCOUPLE=logistic( post$agr ),
			OR_AGR_NOCOUPLE= exp(post$agr),							
			PI_AGR_COUPLE=logistic( post$agr+post$agr_couple ),	
			OR_AGR_COUPLE= exp(post$agr+post$agr_couple),
			PI_TRAD_NOCOUPLE=logistic( post$trad ),
			OR_TRAD_NOCOUPLE=exp( post$trad ),
			PI_TRAD_COUPLE=logistic( post$trad+post$trad_couple ),
			OR_TRAD_COUPLE=exp( post$trad+post$trad_couple ),							
			PI_MXD_NOCOUPLE=logistic( post$cmxd ),
			OR_MXD_NOCOUPLE=exp( post$cmxd ),
			PI_MXD_COUPLE=logistic( post$cmxd+post$cmxd_couple ),	
			OR_MXD_COUPLE=exp( post$cmxd+post$cmxd_couple ),							
			ORX_AGR_NOCOUPLE= exp(post$agr-post$fish),
			ORX_AGR_COUPLE= exp(post$agr+post$agr_couple-post$fish-post$fish_couple),
			ORX_TRAD_NOCOUPLE= exp(post$trad-post$fish),
			ORX_TRAD_COUPLE= exp(post$trad+post$trad_couple-post$fish-post$fish_couple),
			ORX_MXD_NOCOUPLE= exp(post$cmxd-post$fish),
			ORX_MXD_COUPLE= exp(post$cmxd+post$cmxd_couple-post$fish-post$fish_couple),
			ORC_FISH_COUPLE= exp(post$fish_couple),
			ORC_AGR_COUPLE= exp(post$agr_couple),
			ORC_TRAD_COUPLE= exp(post$trad_couple),
			ORC_MXD_COUPLE= exp(post$cmxd_couple)
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.3	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.3	<- dcast.data.table(ds.3, variable~STAT, value.var='V')
	ds.3[, MODEL:= 'couple, community type']
	#	couples female education interaction model
	mg.4 	<- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- 	nocouple*(1-COUPLE3) + nocouple_noedu*(1-COUPLE3)*FE_NOEDU +
							couple*COUPLE3 + couple_noedu*COUPLE3*FE_NOEDU,
					c(couple,nocouple) ~ dnorm(0,100),
					c(couple_noedu,nocouple_noedu) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(couple=0,nocouple=0,couple_noedu=0,nocouple_noedu=0),			
			warmup=5e2, iter=4e3, chains=1, cores=4
	)
	post	<- extract.samples(mg.4)
	dp		<- data.table(	MC= seq_along(post$nocouple),
			PI_YESEDU_NOCOUPLE=logistic( post$nocouple ), 
			OR_YESEDU_NOCOUPLE=exp( post$nocouple ), 
			PI_NOEDU_NOCOUPLE=logistic( post$nocouple+post$nocouple_noedu ),
			OR_NOEDU_NOCOUPLE=exp( post$nocouple+post$nocouple_noedu ),							
			PI_YESEDU_COUPLE=logistic( post$couple ),
			OR_YESEDU_COUPLE=exp( post$couple ),
			PI_NOEDU_COUPLE=logistic( post$couple+post$couple_noedu ),
			OR_NOEDU_COUPLE=exp( post$couple+post$couple_noedu ),							
			ORX_NOEDU_NOCOUPLE=exp( post$nocouple_noedu ), 
			ORX_NOEDU_COUPLE=exp( post$couple_noedu ),
			ORC_YESEDU_COUPLE=exp( post$couple-post$nocouple ),
			ORC_NOEDU_COUPLE=exp( post$couple_noedu-post$nocouple_noedu )
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.4	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.4	<- dcast.data.table(ds.4, variable~STAT, value.var='V')
	ds.4[, MODEL:= 'couple, female education']
	#	couples male education interaction model
	mg.4b 	<- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- 	nocouple*(1-COUPLE3) + nocouple_noedu*(1-COUPLE3)*MA_NOEDU +
							couple*COUPLE3 + couple_noedu*COUPLE3*MA_NOEDU,
					c(couple,nocouple) ~ dnorm(0,100),
					c(couple_noedu,nocouple_noedu) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(couple=0,nocouple=0,couple_noedu=0,nocouple_noedu=0),			
			warmup=5e2, iter=4e3, chains=1, cores=4
	)
	post	<- extract.samples(mg.4b)
	dp		<- data.table(	MC= seq_along(post$nocouple),
			PI_YESMAEDU_NOCOUPLE=logistic( post$nocouple ), 
			OR_YESMAEDU_NOCOUPLE=exp( post$nocouple ), 
			PI_NOMAEDU_NOCOUPLE=logistic( post$nocouple+post$nocouple_noedu ),
			OR_NOMAEDU_NOCOUPLE=exp( post$nocouple+post$nocouple_noedu ),							
			PI_YESMAEDU_COUPLE=logistic( post$couple ),
			OR_YESMAEDU_COUPLE=exp( post$couple ),
			PI_NOMAEDU_COUPLE=logistic( post$couple+post$couple_noedu ),
			OR_NOMAEDU_COUPLE=exp( post$couple+post$couple_noedu ),							
			ORX_NOMAEDU_NOCOUPLE=exp( post$nocouple_noedu ), 
			ORX_NOMAEDU_COUPLE=exp( post$couple_noedu ),
			ORC_YESMAEDU_COUPLE=exp( post$couple-post$nocouple ),
			ORC_NOMAEDU_COUPLE=exp( post$couple_noedu-post$nocouple_noedu )
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.4b	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.4b	<- dcast.data.table(ds.4b, variable~STAT, value.var='V')
	ds.4b[, MODEL:= 'couple, male education']
	# couples circumcision interaction model
	mg.5 	<- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- 	nocouple*(1-COUPLE3) + nocouple_circ*(1-COUPLE3)*MA_CIRCUM +
							couple*COUPLE3 + couple_circ*COUPLE3*MA_CIRCUM,
					c(couple,nocouple) ~ dnorm(0,100),
					c(couple_circ,nocouple_circ) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(couple=0,nocouple=0,nocouple_circ=0,couple_circ=0),			
			warmup=5e2, iter=4e3, chains=1, cores=4
	)
	post	<- extract.samples(mg.5)
	dp		<- data.table(	MC= seq_along(post$nocouple),
			PI_NOCIRC_NOCOUPLE=logistic( post$nocouple ),
			OR_NOCIRC_NOCOUPLE=exp( post$nocouple ),
			PI_CIRC_NOCOUPLE=logistic( post$nocouple+post$nocouple_circ ),	
			OR_CIRC_NOCOUPLE=exp( post$nocouple+post$nocouple_circ ),							
			PI_NOCIRC_COUPLE=logistic( post$couple ),
			OR_NOCIRC_COUPLE=exp( post$couple ),
			PI_CIRC_COUPLE=logistic( post$couple+post$couple_circ ),
			OR_CIRC_COUPLE=exp( post$couple+post$couple_circ ),							
			ORX_CIRC_NOCOUPLE= exp(post$nocouple_circ),
			ORX_CIRC_COUPLE= exp(post$couple_circ),
			ORC_CIRC_COUPLE= exp(post$couple_circ-post$nocouple_circ),
			ORC_NOCIRC_COUPLE= exp(post$couple-post$nocouple)
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.5	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.5	<- dcast.data.table(ds.5, variable~STAT, value.var='V')
	ds.5[, MODEL:= 'couple, male circumcision']
	#	couples "male more sex partners" interaction model
	mg.6 	<- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- 	nocouple*(1-COUPLE3) + nocouple_moresexp*(1-COUPLE3)*MA_SEXP1YR_G1 +
							couple*COUPLE3 + couple_moresexp*COUPLE3*MA_SEXP1YR_G1,
					c(couple,nocouple) ~ dnorm(0,100),
					c(couple_moresexp,nocouple_moresexp) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(couple=0,nocouple=0,nocouple_moresexp=0,couple_moresexp=0),			
			warmup=5e2, iter=4e3, chains=1, cores=4
	)
	post	<- extract.samples(mg.6)
	dp		<- data.table(	MC= seq_along(post$nocouple),
			PI_ONESEXP_NOCOUPLE=logistic( post$nocouple ),
			OR_ONESEXP_NOCOUPLE=exp( post$nocouple ),
			PI_MORESEXP_NOCOUPLE=logistic( post$nocouple+post$nocouple_moresexp ),	
			OR_MORESEXP_NOCOUPLE=exp( post$nocouple+post$nocouple_moresexp ),							
			PI_ONESEXP_COUPLE=logistic( post$couple ),
			OR_ONESEXP_COUPLE=exp( post$couple ),
			PI_MORESEXP_COUPLE=logistic( post$couple+post$couple_moresexp ),
			OR_MORESEXP_COUPLE=exp( post$couple+post$couple_moresexp ),							
			ORX_MORESEXP_NOCOUPLE= exp(post$nocouple_moresexp),
			ORX_MORESEXP_COUPLE= exp(post$couple_moresexp),
			ORC_ONESEXP_COUPLE= exp(post$couple-post$nocouple),
			ORC_MORESEXP_COUPLE= exp(post$couple_moresexp-post$nocouple_moresexp)
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.6	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.6	<- dcast.data.table(ds.6, variable~STAT, value.var='V')
	ds.6[, MODEL:= 'couple, male sex partners']
	#	couples "female more sex partners" interaction model
	mg.6b 	<- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- 	nocouple*(1-COUPLE3) + nocouple_moresexp*(1-COUPLE3)*FE_SEXP1YR_G1 +
							couple*COUPLE3 + couple_moresexp*COUPLE3*FE_SEXP1YR_G1,
					c(couple,nocouple) ~ dnorm(0,100),
					c(couple_moresexp,nocouple_moresexp) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(couple=0,nocouple=0,nocouple_moresexp=0,couple_moresexp=0),			
			warmup=5e2, iter=4e3, chains=1, cores=4
	)
	post	<- extract.samples(mg.6b)
	dp		<- data.table(	MC= seq_along(post$nocouple),
			PI_FEONESEXP_NOCOUPLE=logistic( post$nocouple ),
			OR_FEONESEXP_NOCOUPLE=exp( post$nocouple ),
			PI_FEMORESEXP_NOCOUPLE=logistic( post$nocouple+post$nocouple_moresexp ),	
			OR_FEMORESEXP_NOCOUPLE=exp( post$nocouple+post$nocouple_moresexp ),							
			PI_FEONESEXP_COUPLE=logistic( post$couple ),
			OR_FEONESEXP_COUPLE=exp( post$couple ),
			PI_FEMORESEXP_COUPLE=logistic( post$couple+post$couple_moresexp ),
			OR_FEMORESEXP_COUPLE=exp( post$couple+post$couple_moresexp ),							
			ORX_FEMORESEXP_NOCOUPLE= exp(post$nocouple_moresexp),
			ORX_FEMORESEXP_COUPLE= exp(post$couple_moresexp),
			ORC_FEONESEXP_COUPLE= exp(post$couple-post$nocouple),
			ORC_FEMORESEXP_COUPLE= exp(post$couple_moresexp-post$nocouple_moresexp)
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.6b	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.6b	<- dcast.data.table(ds.6b, variable~STAT, value.var='V')
	ds.6b[, MODEL:= 'couple, female sex partners']
	#	couples "male primary occupation" interaction model
	mg.7 	<- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- 	nocouple_agro*(1-COUPLE3)*MA_OC_AGRO + nocouple_bar*(1-COUPLE3)*MA_OC_BAR +  nocouple_boda*(1-COUPLE3)*MA_OC_BODA +  
							nocouple_fish*(1-COUPLE3)*MA_OC_FISH + nocouple_other*(1-COUPLE3)*MA_OC_OTH +
							nocouple_student*(1-COUPLE3)*MA_OC_STUD + nocouple_trad*(1-COUPLE3)*MA_OC_TRAD +
							couple_agro*COUPLE3*MA_OC_AGRO + couple_bar*COUPLE3*MA_OC_BAR + couple_boda*COUPLE3*MA_OC_BODA + 
							couple_fish*COUPLE3*MA_OC_FISH + couple_other*COUPLE3*MA_OC_OTH +
							couple_student*COUPLE3*MA_OC_STUD + couple_trad*COUPLE3*MA_OC_TRAD,										
					c(nocouple_agro, nocouple_bar, nocouple_boda, nocouple_fish, nocouple_other, nocouple_student, nocouple_trad) ~ dnorm(0,100),
					c(couple_agro, couple_bar, couple_boda, couple_fish, couple_other, couple_student, couple_trad) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	nocouple_agro=0, nocouple_bar=0, nocouple_boda=0, nocouple_fish=0, nocouple_other=0, nocouple_student=0, nocouple_trad=0,
					couple_agro=0, couple_bar=0, couple_boda=0, couple_fish=0, couple_other=0, couple_student=0, couple_trad=0),			
			warmup=5e2, iter=4e3, chains=1, cores=4
	)
	post	<- extract.samples(mg.7)
	dp		<- data.table(	MC= seq_along(post$nocouple_agro),
			PI_AGRO_NOCOUPLE=logistic( post$nocouple_agro ), 
			PI_BAR_NOCOUPLE=logistic( post$nocouple_bar ),	
			PI_BODA_NOCOUPLE=logistic( post$nocouple_boda ),							
			PI_FISH_NOCOUPLE=logistic( post$nocouple_fish ),
			PI_OTH_NOCOUPLE=logistic( post$nocouple_other ),							
			PI_STUD_NOCOUPLE=logistic( post$nocouple_student ),
			PI_TRAD_NOCOUPLE=logistic( post$nocouple_trad ),							
			PI_AGRO_COUPLE=logistic( post$couple_agro ), 
			PI_BAR_COUPLE=logistic( post$couple_bar ),	
			PI_BODA_COUPLE=logistic( post$couple_boda ),							
			PI_FISH_COUPLE=logistic( post$couple_fish ),							
			PI_OTH_COUPLE=logistic( post$couple_other ),							
			PI_STUD_COUPLE=logistic( post$couple_student ),
			PI_TRAD_COUPLE=logistic( post$couple_trad ),							
			OR_AGRO_NOCOUPLE=exp( post$nocouple_agro ), 
			OR_BAR_NOCOUPLE=exp( post$nocouple_bar ),	
			OR_BODA_NOCOUPLE=exp( post$nocouple_boda ),							
			OR_FISH_NOCOUPLE=exp( post$nocouple_fish ),
			OR_OTH_NOCOUPLE=exp( post$nocouple_other ),							
			OR_STUD_NOCOUPLE=exp( post$nocouple_student ),
			OR_TRAD_NOCOUPLE=exp( post$nocouple_trad ),							
			OR_AGRO_COUPLE=exp( post$couple_agro ), 
			OR_BAR_COUPLE=exp( post$couple_bar ),	
			OR_BODA_COUPLE=exp( post$couple_boda ),							
			OR_FISH_COUPLE=exp( post$couple_fish ),
			OR_OTH_COUPLE=exp( post$couple_other ),							
			OR_STUD_COUPLE=exp( post$couple_student ),
			OR_TRAD_COUPLE=exp( post$couple_trad ),					
			ORX_AGRO_NOCOUPLE=exp( post$nocouple_agro - post$nocouple_fish),
			ORX_BAR_NOCOUPLE=exp( post$nocouple_bar - post$nocouple_fish),
			ORX_BODA_NOCOUPLE=exp( post$nocouple_boda - post$nocouple_fish),							
			ORX_OTH_NOCOUPLE=exp( post$nocouple_other - post$nocouple_fish),
			ORX_STUD_NOCOUPLE=exp( post$nocouple_student - post$nocouple_fish),
			ORX_TRAD_NOCOUPLE=exp( post$nocouple_trad - post$nocouple_fish),							
			ORX_AGRO_COUPLE=exp( post$couple_agro - post$couple_fish),
			ORX_BAR_COUPLE=exp( post$couple_bar - post$couple_fish),
			ORX_BODA_COUPLE=exp( post$couple_boda - post$couple_fish),
			ORX_OTH_COUPLE=exp( post$couple_other - post$couple_fish),
			ORX_STUD_COUPLE=exp( post$couple_student - post$couple_fish),
			ORX_TRAD_COUPLE=exp( post$couple_trad - post$couple_fish),							
			ORC_AGRO_COUPLE=exp( post$couple_agro - post$nocouple_agro),
			ORC_BAR_COUPLE=exp( post$couple_bar - post$nocouple_bar),
			ORC_BODA_COUPLE=exp( post$couple_boda - post$nocouple_boda),							
			ORC_FISH_COUPLE=exp( post$couple_fish - post$nocouple_fish),
			ORC_OTH_COUPLE=exp( post$couple_other - post$nocouple_other),
			ORC_STUD_COUPLE=exp( post$couple_student - post$nocouple_student),
			ORC_TRAD_COUPLE=exp( post$couple_trad - post$nocouple_trad)	
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.7	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.7	<- dcast.data.table(ds.7, variable~STAT, value.var='V')
	ds.7[, MODEL:= 'couple, male primary occupation']
	#
	#	couples "male primary occupation" interaction model
	mg.7b 	<- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- 	nocouple_agro*(1-COUPLE3)*FE_OC_AGRO + nocouple_bar*(1-COUPLE3)*FE_OC_BAR +    
							nocouple_other*(1-COUPLE3)*FE_OC_OTH + nocouple_student*(1-COUPLE3)*FE_OC_STUD + nocouple_trad*(1-COUPLE3)*FE_OC_TRAD +
							couple_agro*COUPLE3*MA_OC_AGRO + couple_bar*COUPLE3*MA_OC_BAR +  
							couple_other*COUPLE3*MA_OC_OTH + couple_student*COUPLE3*MA_OC_STUD + couple_trad*COUPLE3*MA_OC_TRAD,										
					c(nocouple_agro, nocouple_bar, nocouple_other, nocouple_student, nocouple_trad) ~ dnorm(0,100),
					c(couple_agro, couple_bar, couple_other, couple_student, couple_trad) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	nocouple_agro=0, nocouple_bar=0, nocouple_other=0, nocouple_student=0, nocouple_trad=0,
					couple_agro=0, couple_bar=0, couple_other=0, couple_student=0, couple_trad=0),			
			warmup=5e2, iter=4e3, chains=1, cores=4
	)
	post	<- extract.samples(mg.7b)
	dp		<- data.table(	MC= seq_along(post$nocouple_agro),
			PI_AGRO_NOCOUPLE=logistic( post$nocouple_agro ), 
			PI_BAR_NOCOUPLE=logistic( post$nocouple_bar ),	
			PI_OTH_NOCOUPLE=logistic( post$nocouple_other ),							
			PI_STUD_NOCOUPLE=logistic( post$nocouple_student ),
			PI_TRAD_NOCOUPLE=logistic( post$nocouple_trad ),							
			PI_AGRO_COUPLE=logistic( post$couple_agro ), 
			PI_BAR_COUPLE=logistic( post$couple_bar ),	
			PI_OTH_COUPLE=logistic( post$couple_other ),							
			PI_STUD_COUPLE=logistic( post$couple_student ),
			PI_TRAD_COUPLE=logistic( post$couple_trad ),							
			OR_AGRO_NOCOUPLE=exp( post$nocouple_agro ), 
			OR_BAR_NOCOUPLE=exp( post$nocouple_bar ),	
			OR_OTH_NOCOUPLE=exp( post$nocouple_other ),							
			OR_STUD_NOCOUPLE=exp( post$nocouple_student ),
			OR_TRAD_NOCOUPLE=exp( post$nocouple_trad ),							
			OR_AGRO_COUPLE=exp( post$couple_agro ), 
			OR_BAR_COUPLE=exp( post$couple_bar ),	
			OR_OTH_COUPLE=exp( post$couple_other ),							
			OR_STUD_COUPLE=exp( post$couple_student ),
			OR_TRAD_COUPLE=exp( post$couple_trad ),					
			ORX_AGRO_NOCOUPLE=exp( post$nocouple_agro - post$nocouple_other),
			ORX_BAR_NOCOUPLE=exp( post$nocouple_bar - post$nocouple_other),
			ORX_OTH_NOCOUPLE=exp( post$nocouple_other - post$nocouple_other),
			ORX_STUD_NOCOUPLE=exp( post$nocouple_student - post$nocouple_other),
			ORX_TRAD_NOCOUPLE=exp( post$nocouple_trad - post$nocouple_other),							
			ORX_AGRO_COUPLE=exp( post$couple_agro - post$couple_other),
			ORX_BAR_COUPLE=exp( post$couple_bar - post$couple_other),							
			ORX_OTH_COUPLE=exp( post$couple_other - post$couple_other),
			ORX_STUD_COUPLE=exp( post$couple_student - post$couple_other),
			ORX_TRAD_COUPLE=exp( post$couple_trad - post$couple_other),							
			ORC_AGRO_COUPLE=exp( post$couple_agro - post$nocouple_agro),
			ORC_BAR_COUPLE=exp( post$couple_bar - post$nocouple_bar),
			ORC_OTH_COUPLE=exp( post$couple_other - post$nocouple_other),
			ORC_STUD_COUPLE=exp( post$couple_student - post$nocouple_student),
			ORC_TRAD_COUPLE=exp( post$couple_trad - post$nocouple_trad)	
	)
	dp		<- melt(dp, id.vars='MC')	
	ds.7b	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.7b	<- dcast.data.table(ds.7b, variable~STAT, value.var='V')
	ds.7b[, MODEL:= 'couple, female primary occupation']
	#
	#	collect results
	ds		<- rbind(ds.1,ds.3,ds.4,ds.4b,ds.5,ds.6,ds.6b,ds.7,ds.7b)
	ds[, LABEL:= paste0(round(MED, d=2), '\n[', round(CL, d=2),'-', round( CU, d=2),']')]
	ds[, COUPLE:= as.character(factor(grepl('_COUPLE',variable), levels=c(TRUE,FALSE), labels=c('couple','no couple')))]
	ds[, STAT:= factor(gsub('^([^_]+)_.*','\\1',variable), levels=c('PI','OR','ORX','ORC'), labels=c('proportion MF','odds ratio MF vs FM','odds ratio factor vs base','odds ratio couple vs casual'))]
	ds[, FACTOR:=gsub('^([^_]+)_([^_]+)_.*','\\2',variable)]
	#
	#	add counts
	tmp		<- df[, list(N=length(LINK_MF)), by=c('FEMALE_OCAT','COUPLE2')]
	set(tmp, NULL, 'FEMALE_OCAT', tmp[, gsub('Student','STUD',gsub('Trading/Shop keeper','TRAD',gsub('zother','OTH',gsub('Bar/waitress','BAR',gsub('Agro/House','AGRO',FEMALE_OCAT)))))])
	setnames(tmp, c('FEMALE_OCAT','COUPLE2'), c('FACTOR','COUPLE'))
	tmp[, MODEL:= 'couple, female primary occupation']
	dss		<- copy(tmp)
	tmp		<- df[, list(N=length(LINK_MF)), by=c('MALE_OCAT','COUPLE2')]
	set(tmp, NULL, 'MALE_OCAT', tmp[, gsub('Boda/Trucking','BODA',gsub('Fishing','FISH',gsub('Student','STUD',gsub('Trading/Shop keeper','TRAD',gsub('zother','OTH',gsub('Bar/waitress','BAR',gsub('Agro/House','AGRO',MALE_OCAT)))))))])
	setnames(tmp, c('MALE_OCAT','COUPLE2'), c('FACTOR','COUPLE'))
	tmp[, MODEL:= 'couple, male primary occupation']
	dss		<- rbind(dss, tmp)
	tmp		<- df[, list(N=length(LINK_MF)), by=c('FE_SEXP1YR_G1','COUPLE2')]
	set(tmp, NULL, 'FE_SEXP1YR_G1', tmp[, gsub('1','FEMORESEXP',gsub('0','FEONESEXP',FE_SEXP1YR_G1))])
	setnames(tmp, c('FE_SEXP1YR_G1','COUPLE2'), c('FACTOR','COUPLE'))
	tmp[, MODEL:= 'couple, female sex partners']
	dss		<- rbind(dss, tmp)
	tmp		<- df[, list(N=length(LINK_MF)), by=c('MA_SEXP1YR_G1','COUPLE2')]
	set(tmp, NULL, 'MA_SEXP1YR_G1', tmp[, gsub('1','MORESEXP',gsub('0','ONESEXP',MA_SEXP1YR_G1))])
	setnames(tmp, c('MA_SEXP1YR_G1','COUPLE2'), c('FACTOR','COUPLE'))
	tmp[, MODEL:= 'couple, male sex partners']
	dss		<- rbind(dss, tmp)	
	tmp		<- df[, list(N=length(LINK_MF)), by=c('MA_CIRCUM','COUPLE2')]
	set(tmp, NULL, 'MA_CIRCUM', tmp[, gsub('1','CIRC',gsub('0','NOCIRC',MA_CIRCUM))])
	setnames(tmp, c('MA_CIRCUM','COUPLE2'), c('FACTOR','COUPLE'))
	tmp[, MODEL:= 'couple, male circumcision']
	dss		<- rbind(dss, tmp)	
	tmp		<- df[, list(N=length(LINK_MF)), by=c('MA_NOEDU','COUPLE2')]
	set(tmp, NULL, 'MA_NOEDU', tmp[, gsub('1','YESMAEDU',gsub('0','NOMAEDU',MA_NOEDU))])
	setnames(tmp, c('MA_NOEDU','COUPLE2'), c('FACTOR','COUPLE'))
	tmp[, MODEL:= 'couple, male education']
	dss		<- rbind(dss, tmp)	
	tmp		<- df[, list(N=length(LINK_MF)), by=c('FE_NOEDU','COUPLE2')]
	set(tmp, NULL, 'FE_NOEDU', tmp[, gsub('1','YESEDU',gsub('0','NOEDU',FE_NOEDU))])
	setnames(tmp, c('FE_NOEDU','COUPLE2'), c('FACTOR','COUPLE'))
	tmp[, MODEL:= 'couple, female education']
	dss		<- rbind(dss, tmp)	
	tmp		<- df[, list(N=length(LINK_MF)), by=c('PAIR_COMM_TYPE','COUPLE2')]
	set(tmp, NULL, 'PAIR_COMM_TYPE', tmp[, gsub('fisherfolk','FISH',gsub('trading','TRAD',gsub('mixed','MXD',gsub('agrarian','AGR',PAIR_COMM_TYPE))))])
	setnames(tmp, c('PAIR_COMM_TYPE','COUPLE2'), c('FACTOR','COUPLE'))
	tmp[, MODEL:= 'couple, community type']
	dss		<- rbind(dss, tmp)	
	tmp		<- df[, list(N=length(LINK_MF)), by=c('COUPLE2')]
	setnames(tmp, c('COUPLE2'), c('COUPLE'))
	tmp[, FACTOR:= 'OVERALL']
	tmp[, MODEL:= 'couple overall']
	dss		<- rbind(dss, tmp)
	dss		<- dcast.data.table(dss, MODEL+FACTOR~COUPLE, value.var='N')
	setnames(dss, c('couple','no couple'), c('N_COUPLE','N_CASUAL'))
	set(dss, dss[, which(is.na(N_CASUAL))], 'N_CASUAL', 0L)
	set(dss, dss[, which(is.na(N_COUPLE))], 'N_COUPLE', 0L)
	#
	save(ds, dss, file=paste0(outfile.base,'_propmf_factors.rda'))	
	# 
	ds		<- dcast.data.table(ds, MODEL+FACTOR~COUPLE+STAT, value.var='LABEL')	
	ds		<- merge(dss, ds, by=c('MODEL','FACTOR'))
	set(ds, NULL, 'MODEL', ds[, factor(MODEL, levels=c("couple overall","couple, community type","couple, male education","couple, female education","couple, male sex partners","couple, female sex partners","couple, male circumcision","couple, male primary occupation","couple, female primary occupation"))])
	set(ds, NULL, 'FACTOR', ds[, factor(FACTOR, levels=c("OVERALL","FISH","OTH","AGR","AGRO","TRAD","MXD","BAR","BODA","STUD","YESMAEDU","NOMAEDU","YESEDU","NOEDU","ONESEXP","MORESEXP","FEONESEXP","FEMORESEXP","NOCIRC","CIRC"))])
	
	setkey(ds, MODEL, FACTOR)
	set(ds, ds[, which(N_CASUAL<10)], c('no couple_proportion MF','no couple_odds ratio MF vs FM','no couple_odds ratio factor vs base'), '?')
	set(ds, ds[, which(N_COUPLE<10)], c('couple_proportion MF','couple_odds ratio MF vs FM','couple_odds ratio factor vs base','couple_odds ratio couple vs casual'), '?')
	for(x in c('couple_odds ratio factor vs base','no couple_odds ratio factor vs base'))
		set(ds, which(is.na(ds[[x]])), x, '-')
	write.csv(ds, row.names=FALSE, file=paste0(outfile.base,'_propmf_factors.csv'))
}


RakaiFull.gender.171122.propfemalepos.communities.merged.into.groups.all<- function()
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
	set(df, NULL, 'COMM_NUM', df[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',as.character(COMM_NUM)))))])
	df	<- subset(df, VISIT%in%c(15,15.1,16))
	
	#	select relevant communities
	tmp	<- sort(unique(c(as.character(rtpdm$MALE_COMM_NUM), as.character(rtpdm$FEMALE_COMM_NUM))))
	tmp	<- data.table(COMM_NUM=tmp)
	df	<- merge(df, tmp, by=c('COMM_NUM'))
	df[, FEMALE:= FEMALE_NEG+FEMALE_POS]
	df[, MALE:= MALE_NEG+MALE_POS]
	df[, POS:= MALE_POS+FEMALE_POS]
	df[, NEG:= MALE_NEG+FEMALE_NEG]
	df[, FEMALE_POS_P:= FEMALE_POS/POS]
	df	<- subset(df, MALE!=0 | FEMALE!=0)
	
	#	divide into groups by proportion of seropos who are female
	tmp	<- df[, list(FEMALE_POS_PM=mean(FEMALE_POS_P)), by='COMM_NUM']
	setkey(tmp, FEMALE_POS_PM)
	tmp[, FEMALE_POS_C:= cut(FEMALE_POS_PM, breaks=c(0.5, 0.55, 0.56, 0.58, 0.68, 0.69,  1), labels=c('a','b','c','d','e','f'))]
	
	ggplot(tmp, aes(x=FEMALE_POS_PM, fill=FEMALE_POS_C)) + geom_histogram()
	df	<- merge(df, subset(tmp, select=c(COMM_NUM, FEMALE_POS_C)), by='COMM_NUM')
	
	#	sum POS and FEMALE_POS by community groups	
	dg	<- df[, list(FEMALE_POS=sum(FEMALE_POS), POS=sum(POS)), by='FEMALE_POS_C']
	
	#	add community groups to rtpdm
	tmp	<- unique(subset(df, select=c(COMM_NUM,FEMALE_POS_C)))
	rtpdm[, COMM_NUM:= FEMALE_COMM_NUM]	
	tmp	<- merge(rtpdm, tmp, by='COMM_NUM')
	tmp	<- tmp[, {
				z<- as.numeric(binconf(length(which(grepl('mf',SELECT))), length(SELECT)))
				list(	MF_TRM	= length(which(grepl('mf',SELECT))),
						MF_TRM_RATIO	= length(which(grepl('mf',SELECT)))/length(which(grepl('fm',SELECT))),
						LMF_TRM_RATIO	= log(length(which(grepl('mf',SELECT)))/length(which(grepl('fm',SELECT)))),
						TRM		= length(SELECT),
						MF_MED	= z[1],
						MF_CL	= z[2],
						MF_CU	= z[3])
			}, by='FEMALE_POS_C']
	dg	<- merge(dg, tmp, by='FEMALE_POS_C')
	tmp	<- dg[, {
				z<- binconf(FEMALE_POS,POS)
				list(FEMALE_POS_M=z[1], FEMALE_POS_CL=z[2], FEMALE_POS_CU=z[3])
			}, by='FEMALE_POS_C']
	dg	<- merge(dg, tmp, by='FEMALE_POS_C')	

	mhi.4 	<- map2stan(
			alist(
					MF_TRM ~ dbinom(TRM, ptm),
					logit(ptm) <- a + b*logit(pdiagf[i]),
					FEMALE_POS ~ dbinom(POS, pdiagf),									 
					a ~ dnorm(0,100),
					b ~ dnorm(0,10),
					pdiagf ~ dnorm(0.5,1)
			),
			data=as.data.frame(dg), start=list(a=0, b=0, pdiagf=rep(0.5,6)),
			warmup=5e2, iter=5e3, chains=1, cores=4
	)	
	summary(mhi.4)	# b posterior median nearly 1; credibility intervals wide but this is good
	post	<- extract.samples(mhi.4)
	dummy	<- function(z) quantile(logistic(with(post, a+b*z)), prob=c(0.5, 0.025, 0.975))
	tmp		<- data.table(FEMALE_POS_P=seq(0.5,0.75,0.001))
	tmp		<- tmp[, 	{
				z<- dummy(logit(FEMALE_POS_P))
				list('MF_median'=z[1], 'MF_cl'=z[2], 'MF_cu'=z[3])
			}, by='FEMALE_POS_P']	
	ggplot(dg) +
			geom_ribbon(data=tmp, aes(x=FEMALE_POS_P, ymin=MF_cl, ymax=MF_cu), fill='black', alpha=0.25) +
			geom_line(data=tmp, aes(x=FEMALE_POS_P, y=MF_median), colour='grey50', size=1) +
			geom_abline(intercept=0, slope=1, lty=2) +
			geom_point(aes(x=FEMALE_POS_M, y=MF_MED, colour=FEMALE_POS_C), size=2) +
			geom_errorbar(aes(x=FEMALE_POS_M, ymin=MF_CL, ymax=MF_CU, colour=FEMALE_POS_C), width=0.03, size=0.9) +
			geom_errorbarh(aes(x=FEMALE_POS_M, y=MF_MED, xmin=FEMALE_POS_CL, xmax=FEMALE_POS_CU, colour=FEMALE_POS_C), height=0.05, size=0.9) +
			scale_x_continuous(labels=scales:::percent, breaks=seq(0,1,0.05)) +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1)) +
			scale_colour_brewer(palette='Set2') +
			theme_bw() +
			labs(	x='\nproportion of females among newly diagnosed', 
					y='male to female transmission\namong phylogenetically inferred transmission events\n',
					colour='community group')
	ggsave(file=paste0(outfile.base,'_trmMF_vs_diagFM_by_commgroup.pdf'), w=6, h=4.5)
}

RakaiFull.gender.171122.highfemaletomaleprevratio.importations<- function()
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
	#	add coordinates
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
	
	#	select relevant communities
	tmp	<- unique(subset(rtpdm, select=c(MALE_COMM_NUM,MALE_COMM_TYPE)))
	setnames(tmp, c('MALE_COMM_NUM','MALE_COMM_TYPE'), c('COMM_NUM','COMM_TYPE'))
	tmp2<- unique(subset(rtpdm, select=c(FEMALE_COMM_NUM,FEMALE_COMM_TYPE)))
	setnames(tmp2, c('FEMALE_COMM_NUM','FEMALE_COMM_TYPE'), c('COMM_NUM','COMM_TYPE'))
	tmp	<- unique(rbind(tmp, tmp2))
	df	<- merge(tmp, df, by='COMM_NUM')
	df[, FEMALE:= FEMALE_NEG+FEMALE_POS]
	df[, MALE:= MALE_NEG+MALE_POS]
	df[, POS:= MALE_POS+FEMALE_POS]
	df[, NEG:= MALE_NEG+FEMALE_NEG]
	#	new clean community numbers
	tmp	<- unique(subset(df, select=c(COMM_NUM,COMM_TYPE)))
	tmp[, COMM_NUM2:= seq_len(nrow(tmp))]
	tmp[, COMM_TYPE2:= as.integer(as.character(factor(COMM_TYPE,levels=c('agrarian','trading','fisherfolk'), labels=c('1','2','3'))))]
	df	<- merge(df, tmp, by=c('COMM_NUM','COMM_TYPE'))
	
	#	estimate prevalence ratio for each community
	dg		<- subset(df, VISIT%in%c(15,15.1,16))	
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
	dgg		<- unique(subset(dg, select=c(COMM_NUM, COMM_TYPE, COMM_NUM2, longitude, latitude)))
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
	#	order communities by ratio, and and make groups
	#
	tmp	<- dcast.data.table(dgg, COMM_NUM+COMM_NUM2+COMM_TYPE~variable, value.var='M')
	#tmp2<- subset(rtpdm, MALE_COMM_NUM==FEMALE_COMM_NUM)[, list(TRM=length(MALE_RID)), by='MALE_COMM_NUM']
	tmp2<- rtpdm[, list(TRM=length(MALE_RID)), by='FEMALE_COMM_NUM']
	setnames(tmp2, 'FEMALE_COMM_NUM', 'COMM_NUM')	
	tmp	<- merge(tmp, tmp2, by='COMM_NUM')		
	tmp	<- tmp[order(RFM),]		
	tmp[, CO_RFM:= seq_len(nrow(tmp))]
	tmp[, CO_RFM_GR:= cut(RFM, breaks=c(1,1.3,1.4,1.62,1.8,2.5,3.4))]
	ggplot(tmp, aes(x=RFM, fill=CO_RFM_GR)) + geom_histogram()
	tmp[, list(sum(TRM)), by='CO_RFM_GR']
	dgg	<- merge(dgg, subset(tmp, select=c(COMM_NUM, CO_RFM_GR)), by='COMM_NUM')
	
	#	calculate extra-comm-transmissions by group
	tmp2<- subset(tmp, select=c(COMM_NUM, CO_RFM_GR))
	rtpdm[, COMM_NUM:=FEMALE_COMM_NUM]
	tmp2<- merge(rtpdm, tmp2, by='COMM_NUM')
	tmp2<- tmp2[, {
				z<- as.numeric(binconf(length(which(FEMALE_COMM_NUM!=MALE_COMM_NUM)), length(SELECT)))
				list(	EXTRACOMM_TRM	= length(which(FEMALE_COMM_NUM!=MALE_COMM_NUM)),
						TRM = length(SELECT),
						EXTRACOMM_TRM_M	= z[1],
						EXTRACOMM_TRM_CL	= z[2],
						EXTRACOMM_TRM_CU	= z[3])
			}, by=c('CO_RFM_GR')]	
	#	calculate avg FM ratio by group and merge
	dh	<- subset(dgg, variable=='RFM')[, list(CL=mean(CL), M=mean(M), CU=mean(CU)), by='CO_RFM_GR']
	dh	<- merge(dh, tmp2, by=c('CO_RFM_GR'))
	
	#	estimate relationship RFM and extra community transmission
	
	dg		<- subset(df, VISIT%in%c(15,15.1,16))	
	tmp		<- rtpdm[, list(EXTRACOMM_TRM= length(which(FEMALE_COMM_NUM!=MALE_COMM_NUM)), TRM= length(SELECT)), by='COMM_NUM']	
	dg		<- merge(dg, tmp, by='COMM_NUM')
	mim.1 	<- map2stan(
			alist(	EXTRACOMM_TRM ~ dbinom(TRM, pimports),
					FEMALE_POS ~ dbinom(FEMALE, seropos_f),
					MALE_POS ~ dbinom(MALE, seropos_m),
					logit(seropos_f) <- afc[COMM_NUM2],
					logit(seropos_m) <- amc[COMM_NUM2],
					logit(pimports) <- aimports + bimports*log_seropos_fm_ratio,
					log_seropos_fm_ratio <- log(1 + exp(-amc[COMM_NUM2])) - log(1 + exp(-afc[COMM_NUM2])),
					afc[COMM_NUM2] ~ dnorm(0, 10),	
					amc[COMM_NUM2] ~ dnorm(0, 10),
					aimports ~ dnorm(0,100),
					bimports ~ dnorm(0,10)
			),
			data=as.data.frame(dg), 
			start=list(aimports=0, bimports=0, afc=rep(0,length(unique(dg$COMM_NUM2))), amc=rep(0,length(unique(dg$COMM_NUM2))) ),
			warmup=5e2, iter=5e3, chains=1, cores=4
		)		
	precis(mim.1, depth=1, prob=0.95)	
	#     Mean StdDev lower 0.95 upper 0.95 n_eff Rhat
	#atm 0.21   0.17      -0.12       0.52  4106    1
	#btm 1.11   0.40       0.35       1.94  3833    1
	#pretty much a straight line ...
	
	post	<- extract.samples(mim.1)
	dummy	<- function(z) quantile(logistic(with(post, aimports+bimports*z)), prob=c(0.5, 0.025, 0.975))
	tmp		<- data.table(FM_PREV_RATIO=seq(1,4,0.01))
	tmp		<- tmp[, {
				z<- dummy(log(FM_PREV_RATIO))
				list('EST_MF_TRM_M'=z[1], 'EST_MF_TRM_CL'=z[2], 'EST_MF_TRM_CU'=z[3])
			}, by='FM_PREV_RATIO']
	tmp2	<- dg[, list(EXTRACOMM_TRM_RAW=mean(EXTRACOMM_TRM/TRM), FM_PREV_RATIO_RAW=mean( (FEMALE_POS/FEMALE)/(MALE_POS/MALE) ) ), by='COMM_NUM']
	ggplot(dh) +
			geom_ribbon(data=tmp, aes(x=FM_PREV_RATIO, ymin=EST_MF_TRM_CL, ymax=EST_MF_TRM_CU), fill='black', alpha=0.25) +
			geom_line(data=tmp, aes(x=FM_PREV_RATIO, y=EST_MF_TRM_M), colour='grey50', size=1) +			
			geom_point(aes(x=M, y=EXTRACOMM_TRM_M, colour=CO_RFM_GR), size=2) +
			geom_errorbar(aes(x=M, ymin=EXTRACOMM_TRM_CL, ymax=EXTRACOMM_TRM_CU, colour=CO_RFM_GR), width=0.1, size=0.9) +
			geom_errorbarh(aes(x=M, y=EXTRACOMM_TRM_M, xmin=CL, xmax=CU, colour=CO_RFM_GR), height=0.03, size=0.9) +
			geom_point(data=tmp2, aes(x=FM_PREV_RATIO_RAW, y=EXTRACOMM_TRM_RAW), colour='black', size=2, pch=17, alpha=1) +
			scale_colour_hue(h = c(-120, 60)) +
			#scale_colour_brewer(palette='Oranges') +
			coord_cartesian( xlim=c(1,3.3), ylim=c(-0.01,1.01) ) +
			theme_bw() +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1), expand=c(0,0)) +
			labs(	x='\nfemale-to-male prevalence ratio', 
					y='proportion of extra-community transmissions\n',
					colour='community group')
	ggsave(file=paste0(outfile.base,'_trmMF_extracommunity_by_commgroup2.pdf'), w=6, h=4.5)
}
	
RakaiFull.gender.171122.prevalanceratio.communities.merged.into.groups<- function()
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
	#	add coordinates
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

	#	select relevant communities
	tmp	<- unique(subset(rtpdm, select=c(MALE_COMM_NUM,MALE_COMM_TYPE)))
	setnames(tmp, c('MALE_COMM_NUM','MALE_COMM_TYPE'), c('COMM_NUM','COMM_TYPE'))
	tmp2<- unique(subset(rtpdm, select=c(FEMALE_COMM_NUM,FEMALE_COMM_TYPE)))
	setnames(tmp2, c('FEMALE_COMM_NUM','FEMALE_COMM_TYPE'), c('COMM_NUM','COMM_TYPE'))
	tmp	<- unique(rbind(tmp, tmp2))
	df	<- merge(tmp, df, by='COMM_NUM')
	df[, FEMALE:= FEMALE_NEG+FEMALE_POS]
	df[, MALE:= MALE_NEG+MALE_POS]
	df[, POS:= MALE_POS+FEMALE_POS]
	df[, NEG:= MALE_NEG+FEMALE_NEG]
	#	new clean community numbers
	tmp	<- unique(subset(df, select=c(COMM_NUM,COMM_TYPE)))
	tmp[, COMM_NUM2:= seq_len(nrow(tmp))]
	tmp[, COMM_TYPE2:= as.integer(as.character(factor(COMM_TYPE,levels=c('agrarian','trading','fisherfolk'), labels=c('1','2','3'))))]
	df	<- merge(df, tmp, by=c('COMM_NUM','COMM_TYPE'))
	
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

	#	estimate prevalence ratio for each community
	dg		<- subset(df, VISIT%in%c(15,15.1,16))	
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
	dgg		<- unique(subset(dg, select=c(COMM_NUM, COMM_TYPE, COMM_NUM2, longitude, latitude)))
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
	ggplot(subset(dgg, variable!='RFM'), aes(x=COMM_NUM)) +
			geom_boxplot(aes(middle=M, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=variable), stat='identity') +
			theme_bw() + 
			facet_grid(~COMM_TYPE, space='free', scale='free') +
			scale_fill_manual(values=c('PF'='hotpink2', 'PM'='deepskyblue')) +
			scale_y_continuous(labels=scales:::percent) +
			labs(x='\ncommunity', y='HIV prevalence\nvisits 15, 15.1, 16\n', fill='gender')
	ggsave(file=paste0(outfile.base,'_trmMF_estimatedprevalence.pdf'), w=12, h=6)
	tmp	<- subset(dgg, variable=='RFM')
	tmp	<- tmp[order(M),]
	tmp[, COMM_NUM3:= seq_len(nrow(tmp))]
	dgg	<- merge(dgg,subset(tmp, select=c(COMM_NUM,COMM_NUM3)), by='COMM_NUM')
	ggplot(subset(dgg, variable=='RFM'), aes(x=COMM_NUM)) +
			geom_boxplot(aes(middle=M, lower=IL, upper=IU, ymin=CL, ymax=CU), stat='identity') +
			theme_bw() + 
			facet_grid(~COMM_TYPE, space='free', scale='free') +
			scale_y_continuous(trans='log', breaks=c(1/2,1,2,3,4,5,6)) +
			labs(x='\ncommunity', y='female to male\nHIV prevalence ratio\nvisits 15, 15.1, 16\n', fill='gender')
	ggsave(file=paste0(outfile.base,'_trmMF_estimatedFMprevalenceratio.pdf'), w=9, h=5)
	
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
	
	#
	#	order communities by ratio, male prev, female prev, and and make groups
	#
	tmp	<- dcast.data.table(dgg, COMM_NUM+COMM_NUM2+COMM_TYPE~variable, value.var='M')
	#tmp2<- subset(rtpdm, MALE_COMM_NUM==FEMALE_COMM_NUM)[, list(TRM=length(MALE_RID)), by='MALE_COMM_NUM']
	tmp2<- rtpdm[, list(TRM=length(MALE_RID)), by='FEMALE_COMM_NUM']
	setnames(tmp2, 'FEMALE_COMM_NUM', 'COMM_NUM')	
	tmp	<- merge(tmp, tmp2, by='COMM_NUM')	
	#	group PF
	tmp	<- tmp[order(PF),]
	ggplot(tmp, aes(x=PF)) + geom_histogram()
	tmp[, CO_PF:= seq_len(nrow(tmp))]	
	tmp[, CO_PF_GR:= cut(PF, breaks=c(0.1,0.15,0.19,0.3,0.7))]
	tmp[, list(sum(TRM)), by='CO_PF_GR']
	#	group PM
	tmp	<- tmp[order(PM),]
	ggplot(tmp, aes(x=PM)) + geom_histogram()
	tmp[, CO_PM:= seq_len(nrow(tmp))]
	tmp[, CO_PM_GR:= cut(PM, breaks=c(0.04,0.092,0.105,0.2,0.4))]
	tmp[, list(sum(TRM)), by='CO_PM_GR']
	#	group RFM
	tmp	<- tmp[order(RFM),]	
	ggplot(tmp, aes(x=RFM)) + geom_histogram()
	tmp[, CO_RFM:= seq_len(nrow(tmp))]
	tmp[, CO_RFM_GR:= cut(RFM, breaks=c(1,1.3,1.4,1.62,1.8,3.4))]
	tmp[, list(sum(TRM)), by='CO_RFM_GR']
	dgg	<- merge(dgg, subset(tmp, select=c(COMM_NUM, CO_PF_GR, CO_PM_GR, CO_RFM_GR)), by='COMM_NUM')
	
	#	save
	save(dgg, file=paste0(outfile.base,'trmMF_estimated_prevalenceratio.rda'))
	
	#	calculate TRM by group
	tmp2<- subset(tmp, select=c(COMM_NUM, CO_PF_GR, CO_PM_GR, CO_RFM_GR))
	rtpdm[, COMM_NUM:=FEMALE_COMM_NUM]
	tmp2<- merge(rtpdm, tmp2, by='COMM_NUM')
	tmp2<- subset(tmp2, select=c(CO_PF_GR, CO_PM_GR, CO_RFM_GR, SELECT))
	tmp2<- melt(tmp2, id.vars='SELECT')	
	tmp2<- tmp2[, {
				z<- as.numeric(binconf(length(which(grepl('mf',SELECT))), length(SELECT)))
				list(	MF_TRM	= length(which(grepl('mf',SELECT))),
						MF_TRM_RATIO	= length(which(grepl('mf',SELECT)))/length(which(grepl('fm',SELECT))),
						LMF_TRM_RATIO	= log(length(which(grepl('mf',SELECT)))/length(which(grepl('fm',SELECT)))),
						TRM		= length(SELECT),
						MF_MED	= z[1],
						MF_CL	= z[2],
						MF_CU	= z[3])
			}, by=c('variable','value')]
	
	#	calculate avg prevalences by group
	dh	<- subset(dgg, variable=='PF')[, list(variable='CO_PF_GR', CL=mean(CL), M=mean(M), CU=mean(CU)), by='CO_PF_GR']
	setnames(dh, 'CO_PF_GR', 'value')
	tmp	<- subset(dgg, variable=='PM')[, list(variable='CO_PM_GR', CL=mean(CL), M=mean(M), CU=mean(CU)), by='CO_PM_GR']
	setnames(tmp, 'CO_PM_GR', 'value')
	dh	<- rbind(dh, tmp)
	tmp	<- subset(dgg, variable=='RFM')[, list(variable='CO_RFM_GR', CL=mean(CL), M=mean(M), CU=mean(CU)), by='CO_RFM_GR']
	setnames(tmp, 'CO_RFM_GR', 'value')
	dh	<- rbind(dh, tmp)
	
	#	merge 
	dh	<- merge(dh, tmp2, by=c('variable','value'))
	
	#	plot with x-axis female to male prev ratio
	tmp	<- subset(dh, variable=='CO_RFM_GR')
	ggplot(tmp) +			
			geom_point(aes(x=M, y=MF_MED, colour=value), size=2) +
			geom_errorbar(aes(x=M, ymin=MF_CL, ymax=MF_CU, colour=value), width=0.03, size=0.9) +
			geom_errorbarh(aes(x=M, y=MF_MED, xmin=CL, xmax=CU, colour=value), height=0.05, size=0.9) +
			scale_colour_brewer(palette='Set2') +			
			theme_bw() +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1)) +
			labs(	x='\nfemale-to-male prevalence ratio', 
					y='probability that male is transmitter\namong phylogenetically inferred transmission events\n',
					colour='community group')
	ggsave(file=paste0(outfile.base,'_trmMF_vs_prevalenceratio_by_commgroup.pdf'), w=6, h=5)
	#	plot with x-axis female prevalence
	tmp	<- subset(dh, variable=='CO_PF_GR')
	ggplot(tmp) +			
			geom_point(aes(x=M, y=MF_MED, colour=value), size=2) +
			geom_errorbar(aes(x=M, ymin=MF_CL, ymax=MF_CU, colour=value), width=0.03, size=0.9) +
			geom_errorbarh(aes(x=M, y=MF_MED, xmin=CL, xmax=CU, colour=value), height=0.05, size=0.9) +
			scale_colour_brewer(palette='Set2') +			
			theme_bw() +
			scale_x_continuous(labels=scales:::percent, breaks=seq(0,1,0.05)) +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1)) +
			labs(	x='\nfemale prevalence', 
					y='probability that male is transmitter\namong phylogenetically inferred transmission events\n',
					colour='community group')
	ggsave(file=paste0(outfile.base,'_trmMF_vs_femaleprevalence_by_commgroup.pdf'), w=6, h=5)
	#	plot with x-axis male prevalence
	tmp	<- subset(dh, variable=='CO_PM_GR')
	ggplot(tmp) +			
			geom_point(aes(x=M, y=MF_MED, colour=value), size=2) +
			geom_errorbar(aes(x=M, ymin=MF_CL, ymax=MF_CU, colour=value), width=0.03, size=0.9) +
			geom_errorbarh(aes(x=M, y=MF_MED, xmin=CL, xmax=CU, colour=value), height=0.05, size=0.9) +
			scale_colour_brewer(palette='Set2') +			
			theme_bw() +
			scale_x_continuous(labels=scales:::percent, breaks=seq(0,1,0.05)) +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1)) +
			labs(	x='\nmale prevalence', 
					y='probability that male is transmitter\namong phylogenetically inferred transmission events\n',
					colour='community group')
	ggsave(file=paste0(outfile.base,'_trmMF_vs_maleprevalence_by_commgroup.pdf'), w=6, h=5)
	#	
	#	estimate relationship between MF transm and FM prev ratio
	#
	dg		<- subset(df, VISIT%in%c(15,15.1,16))
	tmp		<- unique(subset(dgg, select=c(COMM_NUM, CO_PF_GR, CO_PM_GR, CO_RFM_GR)))
	dg2		<- merge(dg, tmp, by='COMM_NUM')
	tmp		<- rtpdm[, list(MF_TRM	= length(which(grepl('mf',SELECT))), TRM= length(SELECT)), by='COMM_NUM']
	dg2		<- merge(dg2, tmp, by='COMM_NUM')
	mpr.2 	<- map2stan(
							alist(
									FEMALE_POS ~ dbinom(FEMALE, seropos_f),
									MALE_POS ~ dbinom(MALE, seropos_m),
									MF_TRM ~ dbinom(TRM, ptm),
									logit(seropos_f) <- afc[COMM_NUM2],
									logit(seropos_m) <- amc[COMM_NUM2],
									logit(ptm) <- atm + btm*log_seropos_fm_ratio,
									log_seropos_fm_ratio <- log(1 + exp(-amc[COMM_NUM2])) - log(1 + exp(-afc[COMM_NUM2])),
									afc[COMM_NUM2] ~ dnorm(0, 10),
									amc[COMM_NUM2] ~ dnorm(0, 10),
									atm ~ dnorm(0,100),
									btm ~ dnorm(0,10)
							),
							data=as.data.frame(dg2), 
							start=list(	atm=0, btm=0,
										afc=rep(0,length(unique(dg$COMM_NUM2))), amc=rep(0,length(unique(dg$COMM_NUM2)))),
							warmup=5e2, iter=5e3, chains=1, cores=4
					)		
	precis(mpr.2, depth=1, prob=0.95)	
	#     Mean StdDev lower 0.95 upper 0.95 n_eff Rhat
	#atm  0.68   0.33       0.05       1.36  1545    1
	#btm -0.85   0.85      -2.54       0.83  1468    1
	#negative,  but not significant
	post	<- extract.samples(mpr.2)
	dummy	<- function(z) quantile(logistic(with(post, atm+btm*z)), prob=c(0.5, 0.025, 0.975))
	tmp		<- data.table(FM_PREV_RATIO=seq(0.5,3.5,0.01))
	tmp		<- tmp[, {
				z<- dummy(log(FM_PREV_RATIO))
				list('EST_MF_TRM_M'=z[1], 'EST_MF_TRM_CL'=z[2], 'EST_MF_TRM_CU'=z[3])
			}, by='FM_PREV_RATIO']	
	tmp2	<- dg2[, list(MF_TRM_RAW=mean(MF_TRM/TRM), FM_PREV_RATIO_RAW=mean( (FEMALE_POS/FEMALE)/(MALE_POS/MALE) ) ), by='COMM_NUM']
	
	ggplot(subset(dh, variable=='CO_RFM_GR')) +
			geom_ribbon(data=tmp, aes(x=FM_PREV_RATIO, ymin=EST_MF_TRM_CL, ymax=EST_MF_TRM_CU), fill='black', alpha=0.25) +
			geom_line(data=tmp, aes(x=FM_PREV_RATIO, y=EST_MF_TRM_M), colour='grey50', size=1) +			
			geom_point(aes(x=M, y=MF_MED, colour=value), size=2) +
			geom_errorbar(aes(x=M, ymin=MF_CL, ymax=MF_CU, colour=value), width=0.1, size=0.9) +
			geom_errorbarh(aes(x=M, y=MF_MED, xmin=CL, xmax=CU, colour=value), height=0.03, size=0.9) +
			geom_point(data=tmp2, aes(x=FM_PREV_RATIO_RAW, y=MF_TRM_RAW), colour='black', size=2, pch=17, alpha=1) +
			scale_colour_hue(h = c(-120, 60)) +
			coord_cartesian( xlim=c(0.95,3.3), ylim=c(-0.01,1.01) ) +
			theme_bw() +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1), expand=c(0,0)) +
			labs(	x='\nfemale-to-male prevalence ratio', 
					y='probability that male is transmitter\namong phylogenetically inferred transmission events\n',
					colour='community group')
	ggsave(file=paste0(outfile.base,'_trmMF_vs_prevalenceratio_by_commgroup_fitted.pdf'), w=6, h=5)
	#	
	#	estimate relationship between MF transm and F prev 
	#
	dg		<- subset(df, VISIT%in%c(15,15.1,16))
	tmp		<- unique(subset(dgg, select=c(COMM_NUM, CO_PF_GR, CO_PM_GR, CO_RFM_GR)))
	dg2		<- merge(dg, tmp, by='COMM_NUM')
	tmp		<- rtpdm[, list(MF_TRM	= length(which(grepl('mf',SELECT))), TRM= length(SELECT)), by='COMM_NUM']
	dg2		<- merge(dg2, tmp, by='COMM_NUM')
	mpr.3 	<- map2stan(
							alist(
									FEMALE_POS ~ dbinom(FEMALE, seropos_f),									
									MF_TRM ~ dbinom(TRM, ptm),
									logit(seropos_f) <- afc[COMM_NUM2],									
									logit(ptm) <- atm + btm*log_seropos_f,
									log_seropos_f <- -log(1 + exp(-afc[COMM_NUM2])),
									afc[COMM_NUM2] ~ dnorm(0, 10),
									atm ~ dnorm(0,100),
									btm ~ dnorm(0,10)
							),
							data=as.data.frame(dg2), 
							start=list(	atm=0, btm=0, afc=rep(0,length(unique(dg$COMM_NUM2)))),
							warmup=5e2, iter=5e3, chains=1, cores=4
						)		
	precis(mpr.3, depth=1, prob=0.95)		
	#   Mean StdDev lower 0.95 upper 0.95 n_eff Rhat
	#	atm  0.15   0.19      -0.22        0.5  4140    1
	#	btm -0.20   0.16      -0.52        0.1  4098    1
	#negative,  but not significant
	post	<- extract.samples(mpr.3)
	dummy	<- function(z) quantile(logistic(with(post, atm+btm*z)), prob=c(0.5, 0.025, 0.975))
	tmp		<- data.table(F_PREV=seq(0.1,0.5,0.005))
	tmp		<- tmp[, {
				z<- dummy(log(F_PREV))
				list('EST_MF_TRM_M'=z[1], 'EST_MF_TRM_CL'=z[2], 'EST_MF_TRM_CU'=z[3])
			}, by='F_PREV']	 
	ggplot(subset(dh, variable=='CO_PF_GR')) +	
			geom_ribbon(data=tmp, aes(x=F_PREV, ymin=EST_MF_TRM_CL, ymax=EST_MF_TRM_CU), fill='black', alpha=0.25) +
			geom_line(data=tmp, aes(x=F_PREV, y=EST_MF_TRM_M), colour='grey50', size=1) +						
			geom_point(aes(x=M, y=MF_MED, colour=value), size=2) +
			geom_errorbar(aes(x=M, ymin=MF_CL, ymax=MF_CU, colour=value), width=0.03, size=0.9) +
			geom_errorbarh(aes(x=M, y=MF_MED, xmin=CL, xmax=CU, colour=value), height=0.05, size=0.9) +
			scale_colour_brewer(palette='Set2') +			
			theme_bw() +
			scale_x_continuous(labels=scales:::percent, breaks=seq(0,1,0.05)) +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1)) +
			labs(	x='\nfemale prevalence', 
					y='probability that male is transmitter\namong phylogenetically inferred transmission events\n',
					colour='community group')
	ggsave(file=paste0(outfile.base,'_trmMF_vs_femaleprevalence_by_commgroup_fitted.pdf'), w=6, h=5)	
}

RakaiFull.gender.171122.propfemalepos.communities.merged.into.groups.noexports<- function()
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
	#	add coordinates
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
	
	#	select relevant communities
	tmp	<- unique(subset(rtpdm, select=c(MALE_COMM_NUM,MALE_COMM_TYPE)))
	setnames(tmp, c('MALE_COMM_NUM','MALE_COMM_TYPE'), c('COMM_NUM','COMM_TYPE'))
	tmp2<- unique(subset(rtpdm, select=c(FEMALE_COMM_NUM,FEMALE_COMM_TYPE)))
	setnames(tmp2, c('FEMALE_COMM_NUM','FEMALE_COMM_TYPE'), c('COMM_NUM','COMM_TYPE'))
	tmp	<- unique(rbind(tmp, tmp2))
	df	<- merge(tmp, df, by='COMM_NUM')
	df[, FEMALE:= FEMALE_NEG+FEMALE_POS]
	df[, MALE:= MALE_NEG+MALE_POS]
	df[, POS:= MALE_POS+FEMALE_POS]
	df[, NEG:= MALE_NEG+FEMALE_NEG]
	#	new clean community numbers
	tmp	<- unique(subset(df, select=c(COMM_NUM,COMM_TYPE)))
	tmp[, COMM_NUM2:= seq_len(nrow(tmp))]
	tmp[, COMM_TYPE2:= as.integer(as.character(factor(COMM_TYPE,levels=c('agrarian','trading','fisherfolk'), labels=c('1','2','3'))))]
	df	<- merge(df, tmp, by=c('COMM_NUM','COMM_TYPE'))
	
	dg	<- subset(df, VISIT%in%c(14,15,15.1,16,17))
	#	the prop of females among HIV+ can vary quite a bit 
	ggplot(dg, aes(x=factor(VISIT), group=COMM_NUM)) +			
			geom_line(aes(y= FEMALE_POS/POS), colour='hotpink2') +			
			geom_point(aes(y= FEMALE_POS/POS), colour='hotpink2') +
			theme_bw() +
			facet_wrap(~COMM_TYPE+COMM_NUM, ncol=4) +
			labs(x='\nvisit', y='gender specific HIV prevalence estimate\n')
	ggsave(file=paste0(outfile.base,'_trmMF_propfemaleprevalence.pdf'), w=9, h=9)
	
	
	#	estimate prop female diagnosed for each community
	dg		<- subset(df, VISIT%in%c(14,15,15.1,16))	
	mfp.1 	<- map2stan(
							alist(
									FEMALE_POS ~ dbinom(POS, seropos_f),					
									logit(seropos_f) <- afc[COMM_NUM2],				
									afc[COMM_NUM2] ~ dnorm(0, 10)								
							),
							data=as.data.frame(dg), 
							start=list(	afc=rep(0,length(unique(dg$COMM_NUM2))) ),
							warmup=5e2, iter=5e3, chains=1, cores=4
						)	
	post	<- extract.samples(mfp.1)			
	dgg		<- unique(subset(dg, select=c(COMM_NUM, COMM_TYPE, COMM_NUM2, longitude, latitude)))
	tmp		<- dgg[, 
			list(	STAT=c('M','CL','IL','IU','CU'),
					PF= as.numeric(quantile(logistic(post$afc[, COMM_NUM2]), prob=c(0.5,0.025,0.25,0.75,0.975)))					
			), 
			by='COMM_NUM2']
	tmp		<- melt(tmp, id.vars=c('COMM_NUM2','STAT'))	
	tmp		<- dcast.data.table(tmp, variable+COMM_NUM2~STAT, value.var='value')
	dgg		<- merge(dgg, tmp, by='COMM_NUM2')
	ggplot(subset(dgg, variable=='PF'), aes(x=COMM_NUM)) +
			geom_boxplot(aes(middle=M, lower=IL, upper=IU, ymin=CL, ymax=CU, fill=variable), stat='identity') +
			theme_bw() + 
			facet_grid(~COMM_TYPE, space='free', scale='free') +
			scale_fill_manual(values=c('PF'='hotpink2')) +
			scale_y_continuous(labels=scales:::percent) +
			labs(x='\ncommunity', y='proportion of female HIV prevalence\nvisits 14, 15, 15.1, 16\n', fill='gender')
	ggsave(file=paste0(outfile.base,'_trmMF_propfemaleprevalence_boxplots.pdf'), w=12, h=6)		
	#
	#	order communities by female prop among prev
	#
	dgg[, FEMALE_POS_C:= cut(M, breaks=c(0.5, 0.55, 0.56, 0.58, 0.69, 0.71,  1))]
	ggplot(dgg, aes(x=M, fill=FEMALE_POS_C)) + geom_histogram()
	
	#	calculate TRM by group
	tmp2<- subset(dgg, select=c(COMM_NUM, FEMALE_POS_C))
	rtpdm[, COMM_NUM:=FEMALE_COMM_NUM]
	tmp2<- merge(rtpdm, tmp2, by='COMM_NUM')
	tmp2<- subset(tmp2, FEMALE_COMM_NUM==MALE_COMM_NUM | (FEMALE_COMM_NUM!=MALE_COMM_NUM & grepl('mf',SELECT)) )	
	tmp2<- tmp2[, {
				z<- as.numeric(binconf(length(which(grepl('mf',SELECT))), length(SELECT)))
				list(	MF_TRM	= length(which(grepl('mf',SELECT))),
						MF_TRM_RATIO	= length(which(grepl('mf',SELECT)))/length(which(grepl('fm',SELECT))),
						LMF_TRM_RATIO	= log(length(which(grepl('mf',SELECT)))/length(which(grepl('fm',SELECT)))),
						TRM		= length(SELECT),
						MF_MED	= z[1],
						MF_CL	= z[2],
						MF_CU	= z[3])
			}, by=c('FEMALE_POS_C')]	
	#	calculate avg female prop by group and merge
	dh	<- dgg[, list(CL=mean(CL), M=mean(M), CU=mean(CU)), by='FEMALE_POS_C']		 
	dh	<- merge(dh, tmp2, by=c('FEMALE_POS_C'))	
	#	
	#	estimate relationship between MF transm and FM prev ratio
	#
	dg		<- subset(df, VISIT%in%c(14,15,15.1,16))
	tmp		<- subset(rtpdm, FEMALE_COMM_NUM==MALE_COMM_NUM | (FEMALE_COMM_NUM!=MALE_COMM_NUM & grepl('mf',SELECT)) )
	tmp		<- tmp[, list(MF_TRM	= length(which(grepl('mf',SELECT))), TRM= length(SELECT)), by='COMM_NUM']
	tmp[, COMM_NUM3:= seq_len(nrow(tmp))]
	dg		<- merge(dg, tmp, by='COMM_NUM')
	mfp.2 	<- map2stan(
			alist(
					FEMALE_POS ~ dbinom(POS, prop_f),					
					MF_TRM ~ dbinom(TRM, ptm),
					logit(prop_f) <- afc[COMM_NUM3],
					logit(ptm) <- atm + btm*logit( 1/(1 + exp(-afc[COMM_NUM3])) ),
					#logit(ptm) <- atm + btm*log( 1/(1 + exp(-afc[COMM_NUM3])) ),
					afc[COMM_NUM3] ~ dnorm(0, 10),					
					atm ~ dnorm(0,100),
					btm ~ dnorm(0,10)
			),
			data=as.data.frame(dg), 
			start=list(	atm=0, btm=0, afc=rep(0,length(unique(dg$COMM_NUM3)))),
			warmup=5e2, iter=5e3, chains=1, cores=4
		)		
	precis(mfp.2, depth=1, prob=0.95)	
	#     Mean StdDev lower 0.95 upper 0.95 n_eff Rhat
	#atm 0.21   0.17      -0.12       0.52  4106    1
	#btm 1.11   0.40       0.35       1.94  3833    1
	#pretty much a straight line ...

	post	<- extract.samples(mfp.2)
	dummy	<- function(z) quantile(logistic(with(post, atm+btm*z)), prob=c(0.5, 0.025, 0.975))
	tmp		<- data.table(FEMALE_PROP=seq(0.4,0.9,0.005))
	tmp		<- tmp[, {
				z<- dummy(logit(FEMALE_PROP))
				list('EST_MF_TRM_M'=z[1], 'EST_MF_TRM_CL'=z[2], 'EST_MF_TRM_CU'=z[3])
			}, by='FEMALE_PROP']
	tmp2	<- dg[, list(MF_TRM_RAW=mean(MF_TRM/TRM), F_PROP_RAW=mean( FEMALE_POS/POS ) ), by='COMM_NUM']
	ggplot(dh) +
			geom_ribbon(data=tmp, aes(x=FEMALE_PROP, ymin=EST_MF_TRM_CL, ymax=EST_MF_TRM_CU), fill='black', alpha=0.25) +
			geom_line(data=tmp, aes(x=FEMALE_PROP, y=EST_MF_TRM_M), colour='grey50', size=1) +			
			geom_point(aes(x=M, y=MF_MED, colour=FEMALE_POS_C), size=2) +
			geom_errorbar(aes(x=M, ymin=MF_CL, ymax=MF_CU, colour=FEMALE_POS_C), width=0.03, size=0.9) +
			geom_errorbarh(aes(x=M, y=MF_MED, xmin=CL, xmax=CU, colour=FEMALE_POS_C), height=0.05, size=0.9) +
			geom_point(data=tmp2, aes(x=F_PROP_RAW, y=MF_TRM_RAW), colour='black', size=2, pch=17, alpha=1) +
			#scale_colour_brewer(palette='Set2') +
			scale_colour_hue(h = c(-120, 60)) +
			coord_cartesian( xlim=c(0.48,0.85), ylim=c(0.3,1.01) ) +
			theme_bw() +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1)) +
			scale_x_continuous(labels=scales:::percent, breaks=seq(0,1,0.1), expand=c(0,0)) +
			labs(	x='\nproportion of females among seropositives', 
					y='probability that male is transmitter\namong phylogenetically inferred transmission events\n',
					colour='community group')
	ggsave(file=paste0(outfile.base,'_trmMF_vs_diagFM_by_commgroup2.pdf'), w=6, h=4.5)	
}


RakaiFull.gender.171122.couples.basicmodel.stan.without.threshold<- function()
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
	
	#zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	#zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	load(infile)	
	setkey(rtp, MALE_RID, FEMALE_RID)
	rtp[, PAIRID:= seq_len(nrow(rtp))]
	rtpdm	<- subset(rtp, grepl('mf|fm|direction not resolved',SELECT))
	rtpdm[, PAIR_COMM_TYPE:= FEMALE_COMM_TYPE]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_TYPE!=MALE_COMM_TYPE)], 'PAIR_COMM_TYPE', 'mixed')
	#rtpdm[, table(PAIR_COMM_TYPE)]
	rtpdm[, COUPLE2:= factor(COUPLE=='no couple', levels=c(TRUE,FALSE), labels=c('no couple','couple'))]
	rtpdm[, SAMEHH:= factor(FEMALE_HH_NUM==MALE_HH_NUM, levels=c(TRUE,FALSE), labels=c('same hh','different hh'))]	
	rtpdm[, PAIR_COMM:= MALE_COMM_NUM_A]
	set(rtpdm, rtpdm[, which(FEMALE_COMM_NUM_A!=MALE_COMM_NUM_A)], 'PAIR_COMM', 'mixed')
	rtpdm[, MALE_SEXP1OUT2:= factor(MALE_SEXP1OUT=='0', levels=c(TRUE,FALSE),labels=c('none','1+'))]
	set(rtpdm, rtpdm[, which(MALE_SEXP1OUT=='Unknown')], 'MALE_SEXP1OUT2', 'Unknown')
	
	
	
	df		<- subset(rtpdm, select=c(MALE_RID, FEMALE_RID, PTY_RUN, IDCLU, LINK_MF, POSTERIOR_SCORE_MF, COUPLE2, SAMEHH, PAIR_COMM_TYPE, MALE_RECENTVL, MALE_CIRCUM, MALE_SEXP1OUT2, MALE_SEXP1YR, MALE_OCCUP_OLLI, FEMALE_EDUCAT, FEMALE_OCCUP_OLLI))
	#	missing data: fill in
	set(df, df[, which(is.na(MALE_RECENTVL))], 'MALE_RECENTVL', -1)	
	set(df, df[, which(is.na(MALE_CIRCUM))], 'MALE_CIRCUM', 'Unknown')
	set(df, df[, which(is.na(FEMALE_EDUCAT))], 'FEMALE_EDUCAT', 'Unknown')	
	#for(x in colnames(df)) print( c(x, any(is.na(df[[x]]))) )
	#	prepare data for STAN
	df[, COUPLE3:= as.integer(COUPLE2=='couple')]
	df[, SAMEHH3:= as.integer(SAMEHH=='same hh')]
	df[, COMM_FISH:= as.integer(PAIR_COMM_TYPE=='fisherfolk')]
	df[, COMM_AGR:= as.integer(PAIR_COMM_TYPE=='agrarian')]
	df[, COMM_TRAD:= as.integer(PAIR_COMM_TYPE=='trading')]
	df[, COMM_MXD:= as.integer(PAIR_COMM_TYPE=='mixed')]
	df[, FE_NOEDU:= as.integer(FEMALE_EDUCAT=='None')]
	df[, FE_NOEDU_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, MA_CIRCUM:= as.integer(MALE_CIRCUM=='Y')]
	df[, MA_CIRCUM_MISS:= as.integer(FEMALE_EDUCAT=='Unknown')]
	df[, MA_SEXP1YR_G1:= as.integer(MALE_SEXP1YR!='1')]
	df[, MA_OC_AGRO:= as.integer(MALE_OCCUP_OLLI=='Agro/House')]
	df[, MA_OC_BAR:= as.integer(MALE_OCCUP_OLLI=='Bar/waitress')]	
	df[, MA_OC_BODA:= as.integer(MALE_OCCUP_OLLI=='Boda/Trucking')]
	df[, MA_OC_CAS:= as.integer(MALE_OCCUP_OLLI=='Casual laborer/unemployed')]
	df[, MA_OC_CONSTR:= as.integer(MALE_OCCUP_OLLI=='Construction/Mechanic')]
	df[, MA_OC_FISH:= as.integer(MALE_OCCUP_OLLI=='Fishing')]
	df[, MA_OC_GOV:= as.integer(MALE_OCCUP_OLLI=='Government/clerical/related')]
	df[, MA_OC_OTH:= as.integer(MALE_OCCUP_OLLI=='Other')]
	df[, MA_OC_STUD:= as.integer(MALE_OCCUP_OLLI=='Student')]
	df[, MA_OC_TRAD:= as.integer(MALE_OCCUP_OLLI=='Trading/Shopkeeper/Hair')]
	
	#
	#	STAN
	#
	
	#	couples effect
	# 	oh shoot, this is far from trivial. I need to work through the STAN manual, but it should be possible.
	mu.1 <- map2stan(
			alist(
					link_mf <- rbinom(length(POSTERIOR_SCORE_MF),1,POSTERIOR_SCORE_MF),
					link_mf ~ dbinom(1, pmf),
					logit(pmf) <- base + couple_b*COUPLE3, 
					base ~ dnorm(0,100),
					couple_b ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(base=logit(0.5), couple_b=0),			
			warmup=5e2, iter=2e3, chains=1, cores=4
	)
	post	<- extract.samples(mg.1)
	dp		<- data.table(	MC= seq_along(post$base),
			PI_OVERALL_NOCOUPLE=logistic( post$base ),
			OR_OVERALL_NOCOUPLE=exp( post$base ),
			PI_OVERALL_COUPLE=logistic( post$base+post$couple_b ),
			OR_OVERALL_COUPLE=exp( post$base+post$couple_b ),
			ORC_OVERALL_COUPLE= exp(post$couple_b))
	dp		<- melt(dp, id.vars='MC')	
	ds.1	<- dp[, list( 	STAT=c('MED','CL','IL','IU','CU'), 
					V= quantile(value, prob=c(0.5,0.025,0.25,0.75,0.975))), by='variable']
	ds.1	<- dcast.data.table(ds.1, variable~STAT, value.var='V')
	ds.1[, MODEL:= 'couple overall']
}

RakaiFull.gender.171122.propfemalepos.community<- function()
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
	set(df, NULL, 'COMM_NUM', df[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',as.character(COMM_NUM)))))])
	df	<- subset(df, VISIT%in%c(15,15.1,16))
	
	tmp	<- sort(unique(c(as.character(rtpdm$MALE_COMM_NUM), as.character(rtpdm$MALE_COMM_NUM))))
	tmp	<- data.table(COMM_NUM=tmp)
	df	<- merge(df, tmp, by=c('COMM_NUM'))
	df	<- df[, list(FEMALE_NEG=sum(FEMALE_NEG), MALE_NEG=sum(MALE_NEG), FEMALE_POS=sum(FEMALE_POS), MALE_POS=sum(MALE_POS)), by='COMM_NUM']
	df[, MALE_PROP_DIAG:= MALE_POS/(MALE_NEG+MALE_POS)]
	df[, FEMALE_PROP_DIAG:= FEMALE_POS/(FEMALE_NEG+FEMALE_POS)]
	df[, FM_DIAG_RATIO:=FEMALE_PROP_DIAG/MALE_PROP_DIAG]
	df[, LFM_DIAG_RATIO:=log(FM_DIAG_RATIO)]
	
	tmp	<- subset(rtpdm, MALE_COMM_NUM==FEMALE_COMM_NUM)
	tmp	<- tmp[, {
				z<- as.numeric(binconf(length(which(grepl('mf',SELECT))), length(SELECT)))
				list(	MF_TRM	= length(which(grepl('mf',SELECT))),
						TRM		= length(SELECT),
						MF_MED	= z[1],
						MF_CL	= z[2],
						MF_CU	= z[3])
			}, by='MALE_COMM_NUM']
	setnames(tmp, 'MALE_COMM_NUM', 'COMM_NUM')
	df	<- merge(df, tmp, by='COMM_NUM')
	ggplot(df, aes(x=MF_DIAG_RATIO, y=MF_MED, ymin=MF_CL, ymax=MF_CU)) +
			geom_point() +
			geom_errorbar() +
			theme_bw()
	
	mh.1 	<- map2stan(
			alist(
					MF_TRM ~ dbinom(TRM, ptmf),
					logit(ptmf) <- a + bdmf*LFM_DIAG_RATIO, 
					a ~ dnorm(0,100),
					bdmf ~ dnorm(0,10)										
			),
			data=as.data.frame(df), start=list(a=0, bdmf=0),
			warmup=5e2, iter=5e3, chains=1, cores=4
	)
	summary(mh.1)	# bdmf not sig above zero
	post	<- extract.samples(mh.1)
	dummy	<- function(dfm) quantile(logistic(with(post, a+bdmf*dfm)), prob=c(0.5, 0.025, 0.975))
	tmp		<- data.table(FM_DIAG_RATIO=seq(1,2.4,0.01))
	tmp		<- tmp[, {
				z<- dummy(log(FM_DIAG_RATIO))
				list('MF_median'=z[1], 'MF_cl'=z[2], 'MF_cu'=z[3])
			}, by='FM_DIAG_RATIO']
	ggplot(df, aes(x=FM_DIAG_RATIO)) +
			geom_point(aes(y=MF_MED), colour='grey50') +
			geom_errorbar(aes(ymin=MF_CL, ymax=MF_CU), colour='grey50') +
			geom_ribbon(data=tmp, aes(ymin=MF_cl, ymax=MF_cu), fill='Darkblue', alpha=0.5) +
			geom_line(data=tmp, aes(y=MF_median), colour='Darkblue') +
			theme_bw()
	ggsave(file=paste0(outfile.base,'_trmMF_vs_diagFM_by_comm.pdf'), w=5, h=5)
	plot(mh.1)	
}