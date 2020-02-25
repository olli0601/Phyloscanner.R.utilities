date2numeric<- function( x )
{
	if(!class(x)%in%c('Date','character'))	return( x )
	x	<- as.POSIXlt(x)
	tmp	<- x$year + 1900
	x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
	x	
}


vl.get.eligible.round17<- function()
{
	require(data.table)
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/rakai_elibility.rda"
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	load(infile)
	
	#	subset to data of interest
	de	<- as.data.table(eldat)	
	de	<- subset(de, status%in%c('_Participated','Away','Blood refusal','Missing data','Other','Refused','urine sample'))
	de	<- subset(de, visit==17)
	
}


vl.vlprops.by.comm.gender.loc<- function()
{
	require(Hmisc)
	require(data.table)
	require(ggplot2)
	vl.detectable <- 4e2
	vl.suppressed <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	load( infile )
	setalloccol(ds)
	
	# remove HIV+ individuals with missing VLs
	ds <- subset(ds, HIV_STATUS==0 | HIV_AND_VL==1)
	
	# define VL_COPIES for uninfected
	set(ds, NULL, 'VLC', ds$VL_COPIES)
	set(ds, ds[,which(HIV_STATUS==0)], 'VLC', 0)
	
	# define undetectable VL (machine-undetectable)
	# define suppressed VL (according to WHO criteria)	
	set(ds, NULL, 'VLU', ds[, as.integer(VLC<vl.detectable)])
	set(ds, NULL, 'VLS', ds[, as.integer(VLC<vl.suppressed)])
	set(ds, NULL, 'VLD', ds[, as.integer(VLC>=vl.detectable)])
	set(ds, NULL, 'VLNS', ds[, as.integer(VLC>=vl.suppressed)])
	set(ds, NULL, 'HIV_AND_VLD', ds[, as.integer(VLD==1 & HIV_AND_VL==1)])
	
	# reset VLC below machine detectable to 0
	set(ds, ds[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 0)
	setkey(ds, FC, SEX, AGEYRS)
	
	# merge two communities that fully overlap, so we have 40 communities in the end 
	set(ds, ds[, which(COMM_NUM==22)], 'COMM_NUM', 1)
	
	# calculate HIV prevalence and proportion not suppressed of HIV+ by community and gender
	vlc <- ds[, {
				z<- as.vector( unname( binconf( length(which(HIV_STATUS==1)), length(HIV_STATUS) ) ) )
				z2<- as.vector( unname( binconf( length(which(VLNS==1)), length(VLNS) ) ) )
				z3<- as.vector( unname( binconf( length(which(VLNS==1)), length(which(HIV_STATUS==1)) ) ) )
				list(FC=FC[1],
						N= length(HIV_STATUS),
						PHIV_MEAN= z[1],
						PHIV_CL= z[2],
						PHIV_CU= z[3],				 
						PVLNS_MEAN= z2[1],
						PVLNS_CL= z2[2],
						PVLNS_CU= z2[3],
						PVLNSofHIV_MEAN= z3[1],
						PVLNSofHIV_CL= z3[2],
						PVLNSofHIV_CU= z3[3],				 
						VLC_MEAN= mean(VLC))		
			}, by=c('COMM_NUM','SEX')]
	vlc[, PHIV_L:= paste0( round(PHIV_MEAN*100, d=1),' [', round(PHIV_CL*100, d=1),'-', round(PHIV_CU*100, d=1),']' )]
	vlc[, PVLNS_L:= paste0( round(PVLNS_MEAN*100, d=1),' [', round(PVLNS_CL*100, d=1),'-', round(PVLNS_CU*100, d=1),']' )]
	vlc[, PVLNSofHIV_L:= paste0( round(PVLNSofHIV_MEAN*100, d=1),' [', round(PVLNSofHIV_CL*100, d=1),'-', round(PVLNSofHIV_CU*100, d=1),']' )]
	setkey(vlc, SEX, PHIV_MEAN)
	set(vlc, NULL, 'SEX', vlc[, factor(SEX, levels=c('M','F'), labels=c('men','women'))])
	
	ggplot(vlc) +
			scale_x_continuous(labels=scales:::percent) +
			scale_y_continuous(labels=scales:::percent) +
			geom_errorbar(aes(x=PHIV_MEAN, ymin=PVLNSofHIV_CL, ymax=PVLNSofHIV_CU), alpha=0.2) +
			geom_errorbarh(aes(y=PVLNSofHIV_MEAN, xmin=PHIV_CL, xmax=PHIV_CU), alpha=0.2) +
			geom_point(aes(x=PHIV_MEAN, y=PVLNSofHIV_MEAN, colour=FC)) +
			geom_text(aes(x=PHIV_MEAN, y=PVLNSofHIV_MEAN, label=COMM_NUM), size=2) +
			facet_wrap(~SEX, ncol=2) +
			theme_bw() +
			labs(x='\nHIV prevalence', 
					y='proportion unsuppressed HIV among infected\n', 
					colour='location')
	ggsave(file=file.path(prjdir,'results_200220','200220_hivnotsuppofhiv_vs_hivprev_by_gender_fishinland.pdf'), w=9, h=5)
		
	#	write results to file
	write.csv(vlc, file=file.path(prjdir,'results_200220','200220_hivnotsuppofhiv_vs_hivprev_by_gender_fishinland.csv'))
	
}

vl.vlrunningprops.by.gender.loc.age<- function()
{
	require(Hmisc)
	require(data.table)
	require(ggplot2)
	vl.detectable <- 4e2
	vl.suppressed <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	load( infile )
	setalloccol(ds)
	
	# remove HIV+ individuals with missing VLs
	ds <- subset(ds, HIV_STATUS==0 | HIV_AND_VL==1)
	# consider only ARVMED for infected
	set(ds, ds[, which(ARVMED==1 & HIV_STATUS==0)], 'ARVMED', 0) 
	
	# define VL_COPIES for uninfected
	set(ds, NULL, 'VLC', ds$VL_COPIES)
	set(ds, ds[,which(HIV_STATUS==0)], 'VLC', 0)
	
	# define undetectable VL (machine-undetectable)
	# define suppressed VL (according to WHO criteria)	
	set(ds, NULL, 'VLU', ds[, as.integer(VLC<vl.detectable)])
	set(ds, NULL, 'VLS', ds[, as.integer(VLC<vl.suppressed)])
	set(ds, NULL, 'VLD', ds[, as.integer(VLC>=vl.detectable)])
	set(ds, NULL, 'VLNS', ds[, as.integer(VLC>=vl.suppressed)])
	set(ds, NULL, 'HIV_AND_VLD', ds[, as.integer(VLD==1 & HIV_AND_VL==1)])
	
	# reset VLC below machine detectable to 0
	set(ds, ds[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 0)
	setkey(ds, FC, SEX, AGEYRS)
	
	tmp <- seq.int(min(ds$AGEYRS), max(ds$AGEYRS))
	vla <- as.data.table(expand.grid(FC=c('fishing','inland'), SEX=c('M','F'), AGEYRS=tmp))
	ans <- vla[, {		
				z <- which(ds$FC==FC & ds$SEX==SEX & ds$AGEYRS<=(AGEYRS+2) & ds$AGEYRS>=(AGEYRS-2))				
				z2<- as.vector( unname( binconf( length(which(ds$HIV_STATUS[z]==1)), length(z) ) ) )
				z3<- as.vector( unname( binconf( length(which(ds$VLNS[z]==1)), length(z) ) ) )
				z4<- as.vector( unname( binconf( length(which(ds$VLNS[z]==1)), length(which(ds$HIV_STATUS[z]==1)) ) ) )
				z5<- as.vector( unname( binconf( length(which(ds$ARVMED[z]==0 & ds$HIV_STATUS[z]==1 & !is.na(ds$ARVMED[z])), length(which(ds$HIV_STATUS[z]==1 & !is.na(ds$ARVMED[z]))) ) ) )
				list(N= length(z),
						PHIV_MEAN= z2[1],
						PHIV_CL= z2[2],
						PHIV_CU= z2[3],				 
						PVLNS_MEAN= z3[1],
						PVLNS_CL= z3[2],
						PVLNS_CU= z3[3],
						PVLNSofHIV_MEAN= z4[1],
						PVLNSofHIV_CL= z4[2],
						PVLNSofHIV_CU= z4[3],
						PARVofHIV_MEAN= z5[1],
						PARVofHIV_CL= z5[2],
						PARVofHIV_CU= z5[3]
						)				
			}, by=c('FC','SEX','AGEYRS')]
	set(ans, NULL, 'SEX', ans[, factor(SEX, levels=c('M','F'), labels=c('men','women'))])
	
	
	#	HIV prevalence
	ggplot(ans) + 		
			geom_ribbon(aes(x=AGEYRS, ymin=PHIV_CL, ymax=PHIV_CU, group=interaction(SEX,FC)), alpha=0.2) +
			geom_line(aes(x=AGEYRS, y=PHIV_MEAN, colour=SEX)) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('men'='royalblue3','women'='deeppink2')) +
			facet_wrap(~FC, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='HIV prevalence (95% CI)\n', 
					colour='gender', 
					linetype='location')
	ggsave(file=file.path(prjdir,'results_200220','200220_hivprevalence_vs_age_by_gender_fishinland.pdf'), w=6, h=5)
	
	#	HIV unsuppressed viral load
	ggplot(ans) + 		
			geom_ribbon(aes(x=AGEYRS, ymin=PVLNS_CL, ymax=PVLNS_CU, group=interaction(SEX,FC)), alpha=0.2) +
			geom_line(aes(x=AGEYRS, y=PARVofHIV_MEAN, colour=SEX), linetype='dotted') +
			geom_line(aes(x=AGEYRS, y=PVLNS_MEAN, colour=SEX)) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('men'='royalblue3','women'='deeppink2')) +
			facet_wrap(~FC, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='proportion unsuppressed HIV (95% CI)\n', 
					colour='gender', 
					linetype='location')
	ggsave(file=file.path(prjdir,'results_200220','200220_hivnotsupp_vs_age_by_gender_fishinland.pdf'), w=6, h=5)
	
	#	HIV unsuppressed viral load among HIV+
	ggplot(ans) + 		
			geom_ribbon(aes(x=AGEYRS, ymin=PVLNSofHIV_CL, ymax=PVLNSofHIV_CU, group=interaction(SEX,FC)), alpha=0.2) +
			geom_line(aes(x=AGEYRS, y=PVLNSofHIV_MEAN, colour=SEX)) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('men'='royalblue3','women'='deeppink2')) +
			facet_wrap(~FC, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='proportion unsuppressed HIV among infected (95% CI)\n', 
					colour='gender', 
					linetype='location')
	ggsave(file=file.path(prjdir,'results_200220','200220_hivnotsuppofhiv_vs_age_by_gender_fishinland.pdf'), w=6, h=5)
	
	
	#	write results to file
	setkey(ans, FC, SEX, AGEYRS)
	ans[, PHIV_L:= paste0( round(PHIV_MEAN*100, d=1),' [', round(PHIV_CL*100, d=1),'-', round(PHIV_CU*100, d=1),']' )]
	ans[, PVLNS_L:= paste0( round(PVLNS_MEAN*100, d=1),' [', round(PVLNS_CL*100, d=1),'-', round(PVLNS_CU*100, d=1),']' )]
	ans[, PVLNSofHIV_L:= paste0( round(PVLNSofHIV_MEAN*100, d=1),' [', round(PVLNSofHIV_CL*100, d=1),'-', round(PVLNSofHIV_CU*100, d=1),']' )]
	write.csv(ans, file=file.path(prjdir,'results_200220','200220_keystats_by_age_gender_fishinland.csv'))
}

vl.vlrunningmean.by.gender.loc.age<- function()
{
	require(Hmisc)
	require(data.table)
	require(ggplot2)
	vl.detectable <- 4e2
	vl.suppressed <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	load( infile )
	setalloccol(ds)
	
	# remove HIV+ individuals with missing VLs
	ds <- subset(ds, HIV_STATUS==0 | HIV_AND_VL==1)
	
	# define VL_COPIES for uninfected
	set(ds, NULL, 'VLC', ds$VL_COPIES)
	set(ds, ds[,which(HIV_STATUS==0)], 'VLC', 0)
	
	# define undetectable VL (machine-undetectable)
	# define suppressed VL (according to WHO criteria)	
	set(ds, NULL, 'VLU', ds[, as.integer(VLC<vl.detectable)])
	set(ds, NULL, 'VLS', ds[, as.integer(VLC<vl.suppressed)])
	set(ds, NULL, 'VLD', ds[, as.integer(VLC>=vl.detectable)])
	set(ds, NULL, 'VLNS', ds[, as.integer(VLC>=vl.suppressed)])
	set(ds, NULL, 'HIV_AND_VLD', ds[, as.integer(VLD==1 & HIV_AND_VL==1)])
	
	# reset VLC below machine detectable to 0
	set(ds, ds[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 0)
	setkey(ds, FC, SEX, AGEYRS)
	
	tmp <- seq.int(min(ds$AGEYRS), max(ds$AGEYRS))
	vla <- as.data.table(expand.grid(FC=c('fishing','inland'), SEX=c('M','F'), AGEYRS=tmp))
	ans <- vla[, {
				z <- which(ds$FC==FC & ds$SEX==SEX & ds$AGEYRS<=(AGEYRS+2) & ds$AGEYRS>=(AGEYRS-2))
				z2 <- mean( ds$VLC[z] )
				z3 <- sd(ds$VLC[z])/sqrt(length(z))
				list(N= length(z),
					VLCM_M= z2,
					VLCM_CL= z2-1.96*z3,
					VLCM_CU= z2+1.96*z3
					)
			}, by=c('FC','SEX','AGEYRS')]
	set(ans, NULL, 'SEX', ans[, factor(SEX, levels=c('M','F'), labels=c('men','women'))])
	
	ggplot(ans) + 
			#geom_errorbar(aes(x=AGEYRS, ymin=VLCM_CL, ymax=VLCM_CU)) +		
			geom_ribbon(aes(x=AGEYRS, ymin=VLCM_CL, ymax=VLCM_CU, group=interaction(SEX,FC)), alpha=0.2) +
			geom_hline(yintercept=1e3) +
			geom_line(aes(x=AGEYRS, y=VLCM_M, colour=SEX)) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous() +
			scale_colour_manual(values=c('men'='royalblue3','women'='deeppink2')) +
			facet_wrap(~FC, ncol=2, scales='free_y') +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='mean viral load (95% CI)\n', 
					colour='gender', 
					linetype='location')
	ggsave(file=file.path(prjdir,'results_200220','200220_vlmean_vs_age_by_gender_fishinland.pdf'), w=8, h=5)
	
	
	#	loess mean below 0 for some age groups, not a good model
	#	sqrt transformation did not work, gave too low means
	ds[, VLCS:= sqrt(VLC)]
	vlclo <- ds[, loess(VLCS ~ AGEYRS, control=loess.control(trace.hat='approximate'))]	
	ans <- subset(ds, select=c(FC, SEX, AGEYRS, VLC, VLCS))	
	ans[, VLCLO_M:= (vlclo$fitted)^2]
	ggplot(ans) + 
			geom_line(aes(x=AGEYRS, y=VLCLO_M)) +
			scale_x_continuous( expand=c(0,0) )	
	ans <- ds[, {
				vlclo <- loess(VLC ~ AGEYRS, control=loess.control(trace.hat='approximate'))
				list(	VLC= VLC,
						AGEYRS= AGEYRS,
						VLCLO_M= vlclo$fitted 
				)				
			}, by=c('FC','SEX')]	
	predict(vlclo, newdata=NULL, se=TRUE)
}

vl.vldistribution.by.gender.loc<- function()
{
	require(Hmisc)
	require(data.table)
	vl.detectable <- 4e2
	vl.suppressed <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	load( infile )
	setalloccol(ds)
	
	# remove HIV+ individuals with missing VLs
	ds <- subset(ds, HIV_STATUS==0 | HIV_AND_VL==1)
	
	# define VL_COPIES for uninfected
	set(ds, NULL, 'VLC', ds$VL_COPIES)
	set(ds, ds[,which(HIV_STATUS==0)], 'VLC', 0)
	
	# define undetectable VL (machine-undetectable)
	# define suppressed VL (according to WHO criteria)	
	set(ds, NULL, 'VLU', ds[, as.integer(VLC<vl.detectable)])
	set(ds, NULL, 'VLS', ds[, as.integer(VLC<vl.suppressed)])
	set(ds, NULL, 'VLD', ds[, as.integer(VLC>=vl.detectable)])
	set(ds, NULL, 'VLNS', ds[, as.integer(VLC>=vl.suppressed)])
	set(ds, NULL, 'HIV_AND_VLD', ds[, as.integer(VLD==1 & HIV_AND_VL==1)])
	
	# reset VLC below machine detectable to 1e-6 (for plotting)
	set(ds, ds[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 1e-6)
	
	# plot proportion of population with viral load > x
	x <- seq(log(1),log(max(ds$VLC)), length.out=1e3)
	x <- c(0,exp(x))
	vld <- as.data.table(expand.grid(X=x, SEX=c('M','F'), FC=c('fishing','inland'), HIV_AND_VLD=c(0,1)))
	
	ans <- vld[, {
				n <- length(which(ds$SEX==SEX & ds$FC==FC & ds$HIV_AND_VLD>=HIV_AND_VLD))
				k <- length(which(ds$SEX==SEX & ds$FC==FC & ds$HIV_AND_VLD>=HIV_AND_VLD & X<ds$VLC))
				z<- as.vector( unname( binconf(k, n) ) )				
				list(N=n, K=k, P_M= z[1], P_CL=z[2], P_CU=z[3] )
			}, by=c('HIV_AND_VLD','FC','SEX','X')]
	set(ans, NULL, 'HIV_AND_VLD', factor(ans[,HIV_AND_VLD], levels=c(0,1), labels=c('all study participants','infected study participants\nwith detectable viral load')))
	set(ans, NULL, 'SEX', ans[, factor(SEX, levels=c('M','F'), labels=c('men','women'))])
	
	ans <- subset(ans, !(HIV_AND_VLD=='infected study participants\nwith detectable viral load' & X<vl.detectable) )
	ans <- subset(ans, !(HIV_AND_VLD=='all study participants' & X<vl.detectable) )
	
	ggplot(ans) +
			geom_line(aes(x=X, y=P_M, group=interaction(FC,SEX), colour=SEX, linetype=FC)) +
			scale_x_log10() +
			scale_y_continuous(labels=scales:::percent, expand=c(0,0)) +
			scale_colour_manual(values=c('men'='royalblue3','women'='deeppink2')) +
			geom_text(aes(x=1e3, y=P_M * 1.03, label="")) +
			theme_bw() +
			facet_wrap(~HIV_AND_VLD, scales='free', ncol=2) +
			labs(	x='\nviral load\n(copies / ml)', 
					y='proportion of individuals with larger viral load\n',
					colour='gender',
					linetype='location')
	ggsave(file=file.path(prjdir,'results_200220','200220_vldistribution_by_gender_fishinland.pdf'), w=9, h=5)
}

vl.keystats.by.gender.loc<- function()
{
	require(Hmisc)
	require(data.table)
	vl.detectable <- 4e2
	vl.suppressed <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	load( infile )
	setalloccol(ds)
	
	# remove HIV+ individuals with missing VLs
	ds <- subset(ds, HIV_STATUS==0 | HIV_AND_VL==1)
	
	# define VL_COPIES for uninfected
	set(ds, NULL, 'VLC', ds$VL_COPIES)
	set(ds, ds[,which(HIV_STATUS==0)], 'VLC', 0)
	
	# define undetectable VL (machine-undetectable)
	# define suppressed VL (according to WHO criteria)	
	set(ds, NULL, 'VLU', ds[, as.integer(VLC<vl.detectable)])
	set(ds, NULL, 'VLS', ds[, as.integer(VLC<vl.suppressed)])
	set(ds, NULL, 'VLD', ds[, as.integer(VLC>=vl.detectable)])
	set(ds, NULL, 'VLNS', ds[, as.integer(VLC>=vl.suppressed)])
	set(ds, NULL, 'HIV_AND_VLD', ds[, as.integer(VLD==1 & HIV_AND_VL==1)])
	ds[, table(VLD, VLNS)]
	#	   VLNS
	#VLD     0     1
  	#0 17577     0
  	#1    90   976
	#	--> there are only 90 individuals with VL in 4e2-1e3, so setting 4e2 or 1e3 is essentially the same
	
	# reset VLC below machine detectable to 0
	set(ds, ds[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 0)
	
	# entire population: 
	# mean viral load, proportion with DVL
	ds[, mean(VLC)]
	# 2290.494
	binconf( length(which(ds$VLD==1)), nrow(ds) )
	# PointEst   Lower      Upper
	# 0.05717964 0.05393704 0.06060469

	#
	# stratified by men/women inland/fishing
	# 
	ans <- ds[, {
			z<- as.vector( unname( binconf( length(which(HIV_STATUS==1)), length(HIV_STATUS) ) ) )
			z2<- as.vector( unname( binconf( length(which(VLNS==1)), length(VLNS) ) ) )
			z3<- as.vector( unname( binconf( length(which(VLNS==1)), length(which(HIV_STATUS==1)) ) ) )
			list(N= length(HIV_STATUS),
				 PHIV_MEAN= z[1],
				 PHIV_CL= z[2],
				 PHIV_CU= z[3],				 
				 PVLNS_MEAN= z2[1],
				 PVLNS_CL= z2[2],
				 PVLNS_CU= z2[3],
				 PVLNSofHIV_MEAN= z3[1],
				 PVLNSofHIV_CL= z3[2],
				 PVLNSofHIV_CU= z3[3],				 
				 VLC_MEAN= mean(VLC))	
			}, by='SEX']
	ans[, FC:='overall']	
	tmp <- ds[, {
				z<- as.vector( unname( binconf( length(which(HIV_STATUS==1)), length(HIV_STATUS) ) ) )
				z2<- as.vector( unname( binconf( length(which(VLNS==1)), length(VLNS) ) ) )
				z3<- as.vector( unname( binconf( length(which(VLNS==1)), length(which(HIV_STATUS==1)) ) ) )
				list(N= length(HIV_STATUS),
						PHIV_MEAN= z[1],
						PHIV_CL= z[2],
						PHIV_CU= z[3],				 
						PVLNS_MEAN= z2[1],
						PVLNS_CL= z2[2],
						PVLNS_CU= z2[3],
						PVLNSofHIV_MEAN= z3[1],
						PVLNSofHIV_CL= z3[2],
						PVLNSofHIV_CU= z3[3],				 						
						VLC_MEAN= mean(VLC))	
			}, by=c('FC','SEX')]
	ans <- rbind(tmp, ans)
	set(ans, NULL, 'FC', factor(ans$FC, levels=c('overall','fishing','inland')))
	set(ans, NULL, 'SEX', factor(ans$SEX, levels=c('F','M')))
	setkey(ans, FC, SEX)
	
	ans[, PHIV_L:= paste0( round(PHIV_MEAN*100, d=1),' [', round(PHIV_CL*100, d=1),'-', round(PHIV_CU*100, d=1),']' )]
	ans[, PVLNS_L:= paste0( round(PVLNS_MEAN*100, d=1),' [', round(PVLNS_CL*100, d=1),'-', round(PVLNS_CU*100, d=1),']' )]
	ans[, PVLNSofHIV_L:= paste0( round(PVLNSofHIV_MEAN*100, d=1),' [', round(PVLNSofHIV_CL*100, d=1),'-', round(PVLNSofHIV_CU*100, d=1),']' )]
	
	#FC SEX    N  PHIV_MEAN    PHIV_CL   PHIV_CU  PVLD_MEAN    PVLD_CL    PVLD_CU PVLDofHIV_MEAN PVLDofHIV_CL PVLDofHIV_CU  VLC_MEAN           PHIV_L           PVLD_L      PVLDofHIV_L
	#1: overall   F 9984 0.21764824 0.20966345 0.2258502 0.05348558 0.04924138 0.05807325      0.2457432    0.2281006    0.2642832 1376.0877   21.8 [21-22.6]    5.3 [4.9-5.8] 24.6 [22.8-26.4]
	#2: overall   M 8659 0.14943989 0.14208609 0.1571046 0.06143897 0.05657296 0.06669393      0.4111283    0.3846208    0.4381619 3344.8235 14.9 [14.2-15.7]    6.1 [5.7-6.7] 41.1 [38.5-43.8]
	#3: fishing   F 1938 0.44272446 0.42074507 0.4649305 0.11455108 0.10112790 0.12949930      0.2587413    0.2305585    0.2890747 3771.5119 44.3 [42.1-46.5] 11.5 [10.1-12.9] 25.9 [23.1-28.9]
	#4: fishing   M 2108 0.32542694 0.30575906 0.3457299 0.14753321 0.13303557 0.16331313      0.4533528    0.4164628    0.4907623 9052.8373 32.5 [30.6-34.6] 14.8 [13.3-16.3] 45.3 [41.6-49.1]
	#5:  inland   F 8046 0.16343525 0.15551676 0.1716750 0.03877703 0.03477391 0.04322036      0.2372624    0.2150559    0.2609994  799.1138 16.3 [15.6-17.2]    3.9 [3.5-4.3] 23.7 [21.5-26.1]
	#6:  inland   M 6551 0.09281026 0.08602036 0.1000774 0.03373531 0.02962926 0.03838786      0.3634868    0.3262210    0.4024669 1508.0821     9.3 [8.6-10]      3.4 [3-3.8] 36.3 [32.6-40.2]
	
	write.csv(ans, file=file.path(prjdir,'results_200220','200220_keystats_by_gender_fishinland.csv'))
}

vl.vlratio.by.loc<- function()
{
	require(Hmisc)
	require(data.table)
	vl.detectable <- 4e2
	vl.suppressed <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	load( infile )
	setalloccol(ds)
	
	# remove HIV+ individuals with missing VLs
	ds <- subset(ds, HIV_STATUS==0 | HIV_AND_VL==1)
	
	# define VL_COPIES for uninfected
	set(ds, NULL, 'VLC', ds$VL_COPIES)
	set(ds, ds[,which(HIV_STATUS==0)], 'VLC', 0)
	
	# define undetectable VL (machine-undetectable)
	# define suppressed VL (according to WHO criteria)	
	set(ds, NULL, 'VLU', ds[, as.integer(VLC<vl.detectable)])
	set(ds, NULL, 'VLS', ds[, as.integer(VLC<vl.suppressed)])
	set(ds, NULL, 'VLD', ds[, as.integer(VLC>=vl.detectable)])
	set(ds, NULL, 'VLNS', ds[, as.integer(VLC>=vl.suppressed)])
	set(ds, NULL, 'HIV_AND_VLD', ds[, as.integer(VLD==1 & HIV_AND_VL==1)])
	ds[, table(VLD, VLNS)]
	
	# reset VLC below machine detectable to 0
	set(ds, ds[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 0)
	
	
	ans	<- as.data.table(expand.grid(BS= 1:1e3, FC=c('fishing','inland')))
	set.seed(42)
	ans <- ans[, {
				zm <- which(ds$FC==FC & ds$SEX=='M')
				zf <- which(ds$FC==FC & ds$SEX=='F')
				zm <- sample(zm, length(zm), replace=TRUE)
				zf <- sample(zf, length(zf), replace=TRUE)
				list(VLCM_M=mean(ds$VLC[zm]), VLCM_F=mean(ds$VLC[zf])) 
			}, by=c('FC','BS')]
	ans[, VLCR:= VLCM_M/VLCM_F]	
	ans <- ans[, list( V=quantile(VLCR, prob=c(0.5, 0.025, 0.975)),
				P= c('M','CL','CU')
				), by=c('FC')]
	ans <- dcast.data.table(ans, FC~P, value.var='V')
	
	#
	# stratified by men/women inland/fishing
	# 
	ans <- ds[, {
				z<- as.vector( unname( binconf( length(which(HIV_STATUS==1)), length(HIV_STATUS) ) ) )
				z2<- as.vector( unname( binconf( length(which(VLNS==1)), length(VLNS) ) ) )
				z3<- as.vector( unname( binconf( length(which(VLNS==1)), length(which(HIV_STATUS==1)) ) ) )
				list(N= length(HIV_STATUS),
						PHIV_MEAN= z[1],
						PHIV_CL= z[2],
						PHIV_CU= z[3],				 
						PVLNS_MEAN= z2[1],
						PVLNS_CL= z2[2],
						PVLNS_CU= z2[3],
						PVLNSofHIV_MEAN= z3[1],
						PVLNSofHIV_CL= z3[2],
						PVLNSofHIV_CU= z3[3],				 
						VLC_MEAN= mean(VLC))	
			}, by='SEX']
	ans[, FC:='overall']	
	tmp <- ds[, {
				z<- as.vector( unname( binconf( length(which(HIV_STATUS==1)), length(HIV_STATUS) ) ) )
				z2<- as.vector( unname( binconf( length(which(VLNS==1)), length(VLNS) ) ) )
				z3<- as.vector( unname( binconf( length(which(VLNS==1)), length(which(HIV_STATUS==1)) ) ) )
				list(N= length(HIV_STATUS),
						PHIV_MEAN= z[1],
						PHIV_CL= z[2],
						PHIV_CU= z[3],				 
						PVLNS_MEAN= z2[1],
						PVLNS_CL= z2[2],
						PVLNS_CU= z2[3],
						PVLNSofHIV_MEAN= z3[1],
						PVLNSofHIV_CL= z3[2],
						PVLNSofHIV_CU= z3[3],				 						
						VLC_MEAN= mean(VLC))	
			}, by=c('FC','SEX')]
	ans <- rbind(tmp, ans)
	set(ans, NULL, 'FC', factor(ans$FC, levels=c('overall','fishing','inland')))
	set(ans, NULL, 'SEX', factor(ans$SEX, levels=c('F','M')))
	setkey(ans, FC, SEX)
}

vl.age.gender<- function()
{
	require(data.table)
	prjdir	<- '~/Box Sync/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	load( infile )
	
	# drop few infecteds with missing VL
	ds <- subset(ds, HIV_STATUS==0 | (HIV_STATUS==1 & !is.na(VL_COPIES)) )
	# set VL for uninfected to 0, and VL with undetectable VL to 0
	set(ds, ds[, which(HIV_STATUS==0)], 'VL_COPIES', 0)
	set(ds, ds[, which(HIV_STATUS==1 & VL_UNDETECTABLE==1)], 'VL_COPIES', 0)
	
	# 
	# calculate proportion with VL > x among participants
	
	# do general by as characters
	# then determine sort index
	# then calculate empirical quantile
	ds <- ds[order(SEX,VL_COPIES),]
	ds[VL]
	
	ds[, sort(unique(VL_COPIES))]
	#dv <- data.table(VL:= )
}

vl.get.data.round17<- function()
{
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- 'data_raw/ViralLoad_Data_Pangea_Ratmann.rda'
	load( file.path(prjdir, infile) )
	
	#
	# subset to survey round 17
	#
	ds		<- subset(as.data.table(survey_data), visit==17)
	# reset dates from Date format to numeric
	for(x in c('visit_date','lastNegDate','firstPosDate'))
	{
		set(ds, NULL, x, date2numeric(ds[[x]]))
	}
	# make all column names upper case
	setnames(ds, colnames(ds), toupper(colnames(ds)))
	# define FISHING_COMM
	ds[, FC:= as.character(factor(COMM_NUM%in%c(770,771,774,38),levels=c(TRUE,FALSE),labels=c('fishing','inland')))]
	# define ARVMED
	set(ds, ds[, which(ARVMED==8)], 'ARVMED', NA_integer_)
	set(ds, NULL, 'ARVMED', ds[, as.integer(as.character(factor(ARVMED, levels=c(1,2), labels=c('1','0'))))])
	
	#
	# prepare GPS coordinates
	#
	dg	<- as.data.table(gpsdat)	
	# bring dates into common format
	setnames(dg, colnames(dg), gsub('\\.','_',toupper(colnames(dg))))
	tmp	<- which(dg[, grepl('([0-9]+)/([0-9]+)/([0-9]+)',GPS_DATE)])
	set(dg, tmp, 'GPS_DATE', dg[tmp,gsub('([0-9]+)/([0-9]+)/([0-9]+)','\\3-\\1-\\2',GPS_DATE)])
	tmp	<- which(dg[, grepl('([0-9]+)-([A-Za-z]+)-([0-9]+)',GPS_DATE)])
	set(dg, tmp, 'GPS_DATE', dg[tmp,gsub('([0-9]+)-([A-Za-z]+)-([0-9]+)','20\\3-\\2-\\1',GPS_DATE)])	
	set(dg, NULL, 'GPS_DATE', dg[,gsub('Nov','11',gsub('Oct','10',gsub('Sep','09',gsub('Aug','08',gsub('July','07',gsub('Jun','06',gsub('May','05',GPS_DATE)))))))])
	# reset dates from character format to numeric
	set(dg, NULL, 'GPS_DATE', date2numeric(dg[,GPS_DATE]))
	# make households per date unique
	dg	<- unique(dg, by=c('HHID','GPS_DATE'))
	
	#
	# add to surveyed individuals the GPS of their households	
	# 
	tmp	<- unique(subset(ds, select=c(RCCS_STUDYID, VISIT_DATE, HHID)))
	tmp	<- merge(tmp, dg, by='HHID', all.x=TRUE)
	# some households do not have GPS coordinates
	ch	<- subset(tmp, is.na(LATITUDE_JITTER) | is.na(LATITUDE_JITTER))
	if(nrow(ch))
	{
		cat('\nNumber of households without GPS coordinates, n=', nrow(ch))
		write.csv(ch, file=file.path(prjdir,'data/check_missing_coordinates.csv'))
		#	521 households without GPS coordinates
	}
	# for every individual, extract house closest in time
	tmp		<- subset(tmp, !is.na(LATITUDE_JITTER) & !is.na(LATITUDE_JITTER))
	tmp2	<- tmp[, list(GPS_DATE= GPS_DATE[which.min(abs(GPS_DATE-VISIT_DATE))[1]]), by=c('RCCS_STUDYID','VISIT_DATE')]
	tmp		<- merge(tmp, tmp2, by=c('RCCS_STUDYID','VISIT_DATE','GPS_DATE'))
	stopifnot(nrow(tmp)==nrow(tmp2))
	set(tmp, NULL, c('COMM','HOUSE'), NULL)	
	ds		<- merge(ds, tmp, by=c('RCCS_STUDYID','VISIT_DATE','HHID'), all.x=TRUE)
	
	#
	# extract viral loads from round 17
	#
	dvl		<- subset(as.data.table(viralLoads), visit==17)
	setnames(dvl, colnames(dvl), toupper(colnames(dvl)))
	setnames(dvl, c('DATE','COPIES','DONEBY'), c('VL_DATE','VL_COPIES','VL_DONEBY'))
	set(dvl, NULL, 'VL_DATE', date2numeric(dvl[,VL_DATE]))
	stopifnot( !nrow(subset(dvl, is.na(VL_COPIES))) )
	# check if one viral load measurement per person
	tmp <- dvl[, list(N_VL=length(VL_COPIES)), by='RCCS_STUDYID']
	stopifnot( !nrow(subset(tmp, N_VL>1)) )	
	# merge with main data
	set(dvl, NULL, 'VISIT', dvl[, as.numeric(VISIT)])
	set(dvl, NULL, 'VL_DONEBY', NULL)
	ds	<- merge(ds, dvl, by=c('RCCS_STUDYID','VISIT'), all.x=TRUE)
	
	# check if viral load for all infected
	ch	<- subset(ds, HIV_STATUS==1 & is.na(VL_COPIES))
	if(nrow(ch))
	{
		cat('\nFound infected individuals without VL measurement, n=',nrow(ch))
		write.csv(ch, file=file.path(prjdir,'data/check_missing_viralloads.csv'))
		# 13 HIV+ individuals without VL
	}
	ds[, HIV_AND_VL:= as.integer(HIV_STATUS==1 & !is.na(VL_COPIES))]
	
	
	save(ds, file=file.path(prjdir,'data','191101_data_round17_vl_gps.rda'))
}

prop.dectectable.viraemia<- function()
{
	require(data.table)
	require(rgdal)
	require(rgeos)
	library(raster)
	require(RColorBrewer) #Map colours
	
	# load data
	infile	<- '~/Box Sync/OR_Work/2018/2018_RakaiViralLoad/data/merged_round17_vl_gps.rda'
	load(infile)
	
	tmp		<- ds[, list(		HIV_POS= length(which(HIV_STATUS==1)), 
					HIV_NEG= length(which(HIV_STATUS==0)),  
					HIV_PREV= length(which(HIV_STATUS==1))/length(HIV_STATUS)
					), by='COMM_NUM']
			
	thr	<- 1e3		
	tmp2	<- subset(ds, HIV_STATUS==1)[, list(		VL_D= length(which(VL_COPIES>thr)), 
					VL_U= length(which(VL_COPIES<=thr)),  
					VL_DP= length(which(VL_COPIES>thr))/length(VL_COPIES)
			), by='COMM_NUM']
	tmp	<- merge(tmp, tmp2, by='COMM_NUM')
	tmp[, POP_VL_DP:= HIV_PREV*VL_DP]
	ggplot(tmp, aes(y=COMM_NUM, x=POP_VL_DP)) + geom_point()
	
	tmp3	<- subset(ds, HIV_STATUS==1 & COMM_NUM==38)
	ggplot(tmp3, aes(x=VL_COPIES)) + geom_histogram() + facet_grid(~ARVMED)
	
	tmp3	<- subset(ds, HIV_STATUS==1 & COMM_NUM==38 & ARVMED==2)
	tmp3[, VL_COPIES_C:= cut(VL_COPIES, breaks=c(0,1,10,100,1000,1e4,1e5,1e6,1e7,1e10), right=FALSE)]
	tmp3[, table(VL_COPIES_C)]
	
	tmp3	<- subset(ds, HIV_STATUS==1 & COMM_NUM==38 & ARVMED==1)
	tmp3[, VL_COPIES_C:= cut(VL_COPIES, breaks=c(0,1,10,100,1000,1e4,1e5,1e6,1e7,1e10), right=FALSE)]
	tmp3[, table(VL_COPIES_C)]
	
	tmp3	<- subset(ds, HIV_STATUS==1 & COMM_NUM==38)
	tmp3[, table(VL_COPIES>1)]
}

make.map.190129	<- function()
{
	require(data.table)
	require(rgdal)
	require(rgeos)
	library(raster)
	require(RColorBrewer) #Map colours
	
	# load data
	infile	<- '~/Box Sync/OR_Work/2018/2018_RakaiViralLoad/data/merged_round17_vl_gps.rda'
	load(infile)
		
	#convert the data into a data table
	dt<- as.data.table(ds)
	dt<- dt[,.(RCCS_STUDYID, SEX, AGEYRS, HIV_STATUS, LATITUDE_JITTER, LONGITUDE_JITTER, VL_COPIES, VL_UNDETECTABLE)]
	#set the NA VL to 0
	dt[is.na(VL_COPIES), VL_COPIES:=0]
	dt[,VL_DETECTABLE := as.numeric(VL_COPIES>=1000)]
	dt[,RCCS_STUDYID2:= seq_len(nrow(dt)) ]
		
	#################################################### load in maps
	# Load in Uganda Shape files 
	uganda1<-raster::getData('GADM',country="UGA",level=1)# Admin unit 1
	uganda3<- raster::getData('GADM', country='UGA', level=3)
	rakai1<-subset(uganda1, NAME_1=="Rakai")
	rakai3<- subset(uganda3, NAME_1=="Rakai")
	masaka1<-subset(uganda1, NAME_1=="Masaka")
	# Create a smaller Rakai for plotting (not current Rakai region no longer includes kabula subdistrict 3)
	#minirak<-rakai3[which(rakai3$NAME_2!="Kabula" | rakai3$NAME_3=="Lyantonde Tc" | rakai3$NAME_3=="Lyantonde"),]
	minirak<-rakai3[which(rakai3$NAME_2!="Kabula" | rakai3$NAME_3=="Lyantonde Tc"),]
	
	####################################################### Convert the data to meters
	#set the coordinates of the data
	coordinates(dt)<- ~ LONGITUDE_JITTER+LATITUDE_JITTER
	#set coordinate system to match uganda files
	proj4string(dt) <- proj4string(uganda1)
	
	#convert to m in order to build a 30x30m grid
	newcrs <- CRS("+proj=robin +datum=WGS84")
	dtnew<- spTransform(dt, newcrs)
	rakai1trans<- spTransform(rakai1, newcrs)
	miniraktrans<- spTransform(minirak, newcrs)
	masaka1trans<- spTransform(masaka1, newcrs)
	
	###################################################### Build Grid
	#Combine rakai1trans and masaka1trans
	outline<- union(rakai1trans, masaka1trans)
	#find the extent of the data
	exnew<- extent(dtnew)
	#extent of the maps
	exmap<- extent(outline)
	
	#chose extent to cover all the data and rakai district
	
	#With a 30m grid, I think the same individuals are usually entering calculations for a large number of grid points
	#Do we really need a 30m grid? Why not 100m?

	grid<- raster(xmn=min(exnew[1], exmap[1]), xmx= exnew[2], ymn=exmap[3], ymx=exnew[4], res=100 )
	#grid[]<- 1:ncell(grid) #No longer needed
	
	# set the coordinate reference system to match
	proj4string(grid)<- proj4string(dtnew) 
	
	#restrict grid to map
	#gridmask<- mask(grid, outline) #Restrict the map after
	#plot(gridmask)
	
	#consider the grid points in a data frame
	id<- as.data.table(1:ncell(gridmask))
	setnames(id, "V1", "ID")
	griddf<- as.data.table(SpatialPoints(grid))
	griddf<- data.table(id, griddf)
	setnames(griddf, gsub('y','LAT_GRID',gsub('x','LONG_GRID',colnames(griddf))))
	
	bw			<- 3000
	bw2			<- bw*bw
	#require(mvtnorm)
	#dmvnorm( c(3.84,0) )	# ~ 9.996634e-05 
	threshold	<- bw*3.84 	# cut if density is < 1e-4
	threshold	<- threshold*threshold	# square the threshold, to avoid sqrt calculations in loop 		
	norm.const	<- 1/(2*pi*bw2)
	
	tmp			<- griddf[1:1e4,]
	anst<- system.time({
		ans	<- tmp[, {
					z1	<- LONG_GRID - dtnew@coords[,'LONGITUDE_JITTER']
					z2	<- LAT_GRID - dtnew@coords[,'LATITUDE_JITTER']
					z1	<- z1*z1 + z2*z2 		# square distance
					z2	<- which(z1<threshold)	# avoid sqrt on 2e4 entries
					w	<- norm.const*exp(-0.5*z1/bw2)	#now with correct normalising constant
					#	Xiayue
					#z3 <-  z1*z1 + z2*z2
					#z4 <- which(z3<threshold)
					#z <- cbind(matrix(z1[z4],ncol=1),matrix(z2[z4],ncol=1))
					#OR: the source code in Boom seems quite slow, with Cholesky decomposition etc. DIY faster?
					#w <- dmvn(z,mu=c(0,0),bw^2*diag(2))	
					#z2 <- z4
					#	olli
					# z1	<- z1*z1 + z2*z2 		# square distance
					# z2	<- which(z1<threshold)	# avoid sqrt on 2e4 entries
					# # to avoid very large output data, calculate directly all smooths here
					# z1	<- sqrt(z1[z2])			# sqrt on few entries					
					# w	<- dnorm(z1, mean=0, sd=bw) # OR: I agree the normalising constant is not right
					# code assumes @coords and @data has same order. 
					list( 	HIV_STATUS_MEAN=mean( dtnew@data$HIV_STATUS[z2] ),				#no weighting by distance
						HIV_STATUS_KERNEL=sum( dtnew@data$HIV_STATUS[z2]*w )/sum(w),		#Gaussian kernel
						VL_COPIES_KERNEL_GEOMMEAN = exp(sum(w*log(dtnew@data$VL_COPIES[z2]+1))/sum(w))-1, #Geometric Mean Kernel
                                                VL_DETECTABLE_KERNEL = sum( dtnew@data$VL_DETECTABLE[z2]*w)/sum(w) #Detectable Prevelance
      	)
				}, by=c('ID','LONG_GRID','LAT_GRID')]
	})
	 grid[]<- ans[, VL_DETECTABLE_KERNEL]
  	gridmask<- mask(grid, outline)
  	#Breaks chosen by looking at data - need refining
  	plot(gridmask, breaks = c(0, 0.025, 0.05, 0.075, 0.1, 0.5), 
       		col=brewer.pal(11, "RdYlGn")[c(10,9,5,4,3)] ,  axes=FALSE, box=FALSE, ylim= c(exmap[3],-6000), legend=FALSE)
  	plot(outline, add=TRUE)
  	par(xpd=TRUE)
  	legend("right", legend=c("0-2.5","2.5-5","5-7.5","7.5-10", ">10"),fill=brewer.pal(11, "RdYlGn")[c(10,9,5,4,3)],horiz = FALSE, inset=-0.175, title= "Prevelence of \n Detectable \n Viremia (%)",  cex=0.8, box.lty = 0)
  	grid[]<- ans[, VL_COPIES_KERNEL_GEOMMEAN]
  	gridmask<- mask(grid, outline)
  	plot(gridmask, breaks = c(0, 0.8, 1.5, 2.5, 3, 145), col=brewer.pal(11, "RdYlGn")[c(10,9,5,4,3)] ,  axes=FALSE, box=FALSE, ylim= c(exmap[3],-6000), legend=FALSE)
  	plot(outline, add=TRUE)
  	par(xpd=TRUE)
  	legend("right", legend=c("0-0.8","0.8-1.5","1.5-2.5","2.5-3", ">3"),fill=brewer.pal(11, "RdYlGn")[c(10,9,5,4,3)],horiz = FALSE, inset=-0.175, title= "Geometric Mean \n VL (Copies/ml)",  cex=0.8, box.lty = 0)
}


