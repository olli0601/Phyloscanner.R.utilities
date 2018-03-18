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

RakaiFull.gender.171122.hhmultivariatemodels.stan.with.threshold<- function()
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
	
	
	df		<- subset(rtpdm, FEMALE_HH_NUM==MALE_HH_NUM, select=c(	MALE_RID, FEMALE_RID, PTY_RUN, IDCLU, LINK_MF, POSTERIOR_SCORE_MF, COUPLE2, SAMEHH, PAIR_COMM_TYPE, 
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
	
	#	same household definitely stronger effect than couples	
	ms.1 <- map2stan(
			alist(
					LINK_MF ~ dbinom(1, pmf),
					logit(pmf) <- base + 
							# comm_fish*COMM_FISH +	this is base 
							comm_agr*COMM_AGR + comm_trad*COMM_TRAD + comm_mxd*COMM_MXD +
							noedu_f*FE_NOEDU + noedu_m*MA_NOEDU +
							sexp_f*FE_SEXP1YR_G1 + sexp_m*MA_SEXP1YR_G1 +
							circum_m*MA_CIRCUM +										
							# other occupation is baseline
							mocc_agro*MA_OC_AGRO + mocc_fish*MA_OC_FISH + mocc_trad*MA_OC_TRAD +										
							focc_agro*FE_OC_AGRO + focc_bar*FE_OC_BAR +  focc_trad*FE_OC_TRAD,
					base ~ dnorm(0,100),								
					c(comm_agr, comm_trad, comm_mxd) ~ dnorm(0,10),
					c(noedu_f, noedu_m) ~ dnorm(0,10),								
					c(sexp_f, sexp_m, circum_m) ~ dnorm(0,10),
					c(mocc_agro, mocc_fish, mocc_trad) ~ dnorm(0,10),
					c(focc_agro, focc_bar, focc_trad) ~ dnorm(0,10)
			),
			data=as.data.frame(df), 
			start=list(	base=0, comm_agr=0, comm_trad=0, comm_mxd=0, noedu_f=0, noedu_m=0, sexp_f=0, sexp_m=0, circum_m=0, 
					focc_agro=0, focc_bar=0, focc_trad=0,
					mocc_agro=0, mocc_fish=0, mocc_trad=0),			
			warmup=1e3, iter=5e3, chains=1, cores=4
	)		
	#plot(precis(ms.1, prob=0.95))
	#pairs(mm.2)
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
											"COMM-AGRO","COMM-TRAD","COMM-MXD",																						
											"FEDU-NONE", "MEDU-NONE",      
											"FESEXP1YR-G1","MASEXP1YR-G1", 
											"CIRC-YES",                   
											"FOCC-AGRO","FOCC-BAR","FOCC-TRAD",
											"MOCC-AGRO","MOCC-FISH","MOCC-TRAD"   
									)), labels=rev(c("Agrarian community vs. fishing community", "Trading community vs. fishing community", "Mixed vs. fishing community",
											"Female education: None vs. at least primary education", "Male education: None vs. at least primary education",
											"Female sex partners in last year: >1 vs. 1", "Male sex partners in last year: >1 vs. 1",
											"Male circumcised: yes vs. no",
											"Female occupation: agricultural/house vs. other", "Female occupation: bar/waitress vs. other", "Female occupation: trading/shopkeeper vs. other",
											"Male occupation: agricultural/house vs. other", "Male occupation: fishing vs. other", "Male occupation: trading/shopkeeper vs. other"
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


RakaiFull.gender.171122.couplesmultivariateadditivemodels.stan.with.threshold<- function()
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

RakaiFull.gender.171122.couplesinteractionmodels.stan.with.threshold<- function()
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
	
	#
	#	plot FM ratio on map
	#	
	require(ggmap)
	zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")	
	#
	#	plot number of observed recipients and observed transmitters
	ggmap(zm) +
			geom_point(data=subset(dgg, variable=='RFM'), aes(x=longitude, y=latitude, pch=COMM_TYPE, colour=M), size=7) +
			geom_text(data=subset(dgg, variable=='RFM'), aes(x=longitude, y=latitude, label=COMM_NUM), nudge_x=0, nudge_y=0, size=3, colour='black') + 
			scale_colour_gradient2(trans='log', breaks=c(1,1.25,1.5,2,2.5,3), low="deepskyblue", mid="orange", high='red', midpoint=log(2.2)) +
			labs(x='\nlongitude',y='latitude\n',colour='female to male\nprevalence ratio', pch='community\ntype')
	ggsave(file=paste0(outfile.base,'trmMF_prevalenceratio_on_map.pdf'), w=7, h=7)
		
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
	tmp		<- data.table(FM_PREV_RATIO=seq(1.01,2.5,0.01))
	tmp		<- tmp[, {
				z<- dummy(log(FM_PREV_RATIO))
				list('EST_MF_TRM_M'=z[1], 'EST_MF_TRM_CL'=z[2], 'EST_MF_TRM_CU'=z[3])
			}, by='FM_PREV_RATIO']
	ggplot(subset(dh, variable=='CO_RFM_GR')) +
			geom_ribbon(data=tmp, aes(x=FM_PREV_RATIO, ymin=EST_MF_TRM_CL, ymax=EST_MF_TRM_CU), fill='black', alpha=0.25) +
			geom_line(data=tmp, aes(x=FM_PREV_RATIO, y=EST_MF_TRM_M), colour='grey50', size=1) +			
			geom_point(aes(x=M, y=MF_MED, colour=value), size=2) +
			geom_errorbar(aes(x=M, ymin=MF_CL, ymax=MF_CU, colour=value), width=0.03, size=0.9) +
			geom_errorbarh(aes(x=M, y=MF_MED, xmin=CL, xmax=CU, colour=value), height=0.05, size=0.9) +
			scale_colour_brewer(palette='Set2') +			
			theme_bw() +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1)) +
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
	tmp		<- data.table(FEMALE_PROP=seq(0.5,0.85,0.005))
	tmp		<- tmp[, {
				z<- dummy(logit(FEMALE_PROP))
				list('EST_MF_TRM_M'=z[1], 'EST_MF_TRM_CL'=z[2], 'EST_MF_TRM_CU'=z[3])
			}, by='FEMALE_PROP']
	ggplot(dh) +
			geom_ribbon(data=tmp, aes(x=FEMALE_PROP, ymin=EST_MF_TRM_CL, ymax=EST_MF_TRM_CU), fill='black', alpha=0.25) +
			geom_line(data=tmp, aes(x=FEMALE_PROP, y=EST_MF_TRM_M), colour='grey50', size=1) +			
			geom_point(aes(x=M, y=MF_MED, colour=FEMALE_POS_C), size=2) +
			geom_errorbar(aes(x=M, ymin=MF_CL, ymax=MF_CU, colour=FEMALE_POS_C), width=0.03, size=0.9) +
			geom_errorbarh(aes(x=M, y=MF_MED, xmin=CL, xmax=CU, colour=FEMALE_POS_C), height=0.05, size=0.9) +
			scale_colour_brewer(palette='Set2') +			
			theme_bw() +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.1)) +
			labs(	x='\nproportion of females among seropositives', 
					y='probability that male is transmitter\namong phylogenetically inferred transmission events\n',
					colour='community group')
	ggsave(file=paste0(outfile.base,'_trmMF_vs_diagFM_by_commgroup2.pdf'), w=6, h=4.5)		
}


RakaiFull.gender.171122.couplesbasicmodel.stan.without.threshold<- function()
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