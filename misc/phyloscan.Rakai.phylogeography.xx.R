Rakai190327.participitation.differences.betabinomialmodel3<- function()
{
	require(data.table)
	require(rstan)
	require(extraDistr)
	
	indir		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis'
	infile.data	<- file.path(indir,"180322_sampling_by_gender_age.rda")
	infile.participation.stan.model <- file.path(indir,"180322_glm_participation.stan")
	outfile.base<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
			
	#	set up variables for STAN
	load(infile.data)
	des		<- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID, SEX, INMIGRANT))	
	
	# remove individuals in the communities where no sequences were obtained successfully
	tmp	<- des[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
	des	<- merge(des, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')
	
	
	# set no INMIGRANT to 0
	set(des, des[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
	
	# age group
	des[, AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,25,35,60), labels = c("15-24", "25-34","35+"))]
	
	#	binarize age, sex
	des[, AGE_YOUNG:= as.integer(AGE_AT_MID_C=='15-24')]
	des[, AGE_MID:= as.integer(AGE_AT_MID_C=='25-34')]
	des[, MALE:= as.integer(SEX=='M')]
	
	#	binarize community type
	des[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	des[, COMM_TYPE_T:= as.integer(substr(COMM_NUM_A,1,1)=='t')]
			
	#	vanilla community IDs
	tmp		<- data.table(COMM_NUM_A= sort(unique(des$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(des$COMM_NUM_A)))
	des		<- merge(des, tmp, by='COMM_NUM_A')
	
	#	aggregate participation
	dp		<- des[, list(TRIAL=length(PERM_ID), SUC=length(which(PARTICIPATED==1))), by=c('AGE_AT_MID_C','SEX','AGE_YOUNG','AGE_MID','INMIGRANT','MALE','COMM_TYPE_F','COMM_TYPE_T','COMM_NUM_A','COMM_NUM_B')]
	dp[, CATEGORY:= paste0(COMM_NUM_A,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
	
	#	run STAN 
	tmp			<- as.list(subset(dp,select=c('COMM_NUM_B','TRIAL','SUC','MALE','AGE_YOUNG','AGE_MID','INMIGRANT','COMM_TYPE_F','COMM_TYPE_T')))
	tmp$N		<- nrow(dp)
	tmp$N_COMM	<- length(unique(dp$COMM_NUM_A))
	fit.par 	<- stan(	file = infile.participation.stan.model, 
						data = tmp, 
						iter = 10e3,
						warmup = 5e2,
						cores = 1,
						chains = 1,
						init = list(list(a=0, comm=rep(0,36), sig_comm=1, dispersion=1, male=0, midage=0, female_young=0, male_young=0, inmigrant_young=0, inmigrant=0, trading=0, fishing=0)))
	
	# assess convergence
	fit.pars	<- c('a','comm','dispersion','sig_comm','female_young','male_young','inmigrant_young','male','midage','inmigrant','trading','fishing')
	any(rhat(fit.par, pars=fit.pars)>1.02)
	any(neff_ratio(fit.par, pars=fit.pars) * 9.5e3 < 500)
	po              <- as.matrix(fit.par) 
	po              <- po[, colnames(po)[!grepl('p_suc|lp__',colnames(po))]]	
	p   <- mcmc_areas(po, pars=colnames(po), prob = 0.95) + 
			geom_vline(xintercept = 0)
	pdf(file=paste0(outfile.base,'190327_participation_model_marginalposteriors.pdf'), w=7, h=20)
	p
	dev.off()
	p   <- mcmc_trace(po, pars=colnames(po), facet_args = list(ncol = 1)) 
	pdf(file=paste0(outfile.base,'190327_participation_model_marginaltraces.pdf'), w=7, h=100)
	p
	dev.off()
	
	# summary csv file
	write.csv(  summary(fit.par, pars= fit.pars, probs = c(0.025, 0.975))$summary,
				file=paste0(outfile.base,'190327_participation_model_summary.csv')
				)
	
	# extract samples for unique strata levels
	nprior		<- 500
	fit.e		<- extract(fit.par)
	set.seed(42)
	tmp			<- sample(length(fit.e$a), nprior)
	dps			<- dp[,	
			{
				z<- with(fit.e, a + comm[,COMM_NUM_B] + male * MALE + 
								trading*COMM_TYPE_T  + fishing*COMM_TYPE_F +
								inmigrant*INMIGRANT + inmigrant_young*INMIGRANT*AGE_YOUNG +
								male_young*AGE_YOUNG*MALE + female_young*AGE_YOUNG*(1-MALE) + midage*AGE_MID)
				list(SAMPLE=1:nprior, ETA=as.numeric(z[tmp]), DISPERSION=as.numeric(fit.e$dispersion[tmp]))
			},	
			by=c('CATEGORY','TRIAL','SUC')]
	dps[, P:= exp(ETA)/(1+exp(ETA))]
	dps[, LP:= dbbinom(SUC, TRIAL, alpha= P*DISPERSION, beta= (1-P)*DISPERSION, log=TRUE)]	
	set(dps, NULL, c('TRIAL','SUC','ETA','DISPERSION'), NULL)
		
	save(dp, dps, fit.par, file=paste0(outfile.base,'190327_participation_model_samples.rda'))
}
Rakai190327.participitation.differences.betabinomialmodel4<- function()
{
	require(data.table)
	require(rstan)
	require(extraDistr)
	
	indir		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis'
	infile.data	<- file.path(indir,"180322_sampling_by_gender_age.rda")
	infile.participation.stan.model <- file.path(indir,"180322_glm_participation_age6.stan")
	outfile.base<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	
	#	set up variables for STAN
	load(infile.data)
	des		<- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID, SEX, INMIGRANT))	
	
	# remove individuals in the communities where no sequences were obtained successfully
	tmp	<- des[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
	des	<- merge(des, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')
	
	
	# set no INMIGRANT to 0
	set(des, des[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
	
	# age group
	des[,AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("20-", "20-24","25-29","30-34","35-39","40-44", "45+"))]
	
	#	binarize age, sex
	des[, AGE1:= as.integer(AGE_AT_MID_C=='20-')]
	des[, AGE2:= as.integer(AGE_AT_MID_C=='20-24')]
	des[, AGE3:= as.integer(AGE_AT_MID_C=='25-29')]
	des[, AGE4:= as.integer(AGE_AT_MID_C=='30-34')]
	des[, AGE5:= as.integer(AGE_AT_MID_C=='35-39')]
	des[, AGE6:= as.integer(AGE_AT_MID_C=='40-44')]
	des[, MALE:= as.integer(SEX=='M')]	
	
	#	binarize community type
	des[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	des[, COMM_TYPE_T:= as.integer(substr(COMM_NUM_A,1,1)=='t')]
	
	#	vanilla community IDs
	tmp		<- data.table(COMM_NUM_A= sort(unique(des$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(des$COMM_NUM_A)))
	des		<- merge(des, tmp, by='COMM_NUM_A')
	
	#	aggregate participation
	dp		<- des[, list(TRIAL=length(PERM_ID), SUC=length(which(PARTICIPATED==1))), by=c('AGE_AT_MID_C','SEX','AGE1','AGE2','AGE3','AGE4','AGE5','AGE6','INMIGRANT','MALE','COMM_TYPE_F','COMM_TYPE_T','COMM_NUM_A','COMM_NUM_B')]
	dp[, CATEGORY:= paste0(COMM_NUM_A,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
	
	#	run STAN 
	tmp			<- as.list(subset(dp,select=c('COMM_NUM_B','TRIAL','SUC','MALE','AGE1','AGE2','AGE3','AGE4','AGE5','AGE6','INMIGRANT','COMM_TYPE_F','COMM_TYPE_T')))
	tmp$N		<- nrow(dp)
	tmp$N_COMM	<- length(unique(dp$COMM_NUM_A))
	fit.par 	<- stan(	file = infile.participation.stan.model, 
			data = tmp, 
			iter = 10e3,
			warmup = 5e2,
			cores = 1,
			chains = 1,
			init = list(list(a=0, comm=rep(0,36), sig_comm=1, dispersion=1, male=0, female_age1=0, male_age1=0, female_age2=0, male_age2=0, age3=0, age4=0, age5=0, age6=0, inmigrant_age1=0, inmigrant_age2=0, inmigrant=0, trading=0, fishing=0)))
	
	# assess convergence
	fit.pars	<- c('a','comm','dispersion','sig_comm','female_age1','male_age1','female_age2','male_age2','inmigrant_age1','inmigrant_age2','age3','age4','age5','age6','male','inmigrant','trading','fishing')
	any(rhat(fit.par, pars=fit.pars)>1.02)
	any(neff_ratio(fit.par, pars=fit.pars) * 9.5e3 < 500)
	po              <- as.matrix(fit.par) 
	po              <- po[, colnames(po)[!grepl('p_suc|lp__',colnames(po))]]	
	p   <- mcmc_areas(po, pars=colnames(po), prob = 0.95) + 
			geom_vline(xintercept = 0)
	pdf(file=paste0(outfile.base,'190327_participation_modelage6_marginalposteriors.pdf'), w=7, h=20)
	p
	dev.off()
	p   <- mcmc_trace(po, pars=colnames(po), facet_args = list(ncol = 1)) 
	pdf(file=paste0(outfile.base,'190327_participation_modelage6_marginaltraces.pdf'), w=7, h=100)
	p
	dev.off()
	
	# summary csv file
	write.csv(  summary(fit.par, pars= fit.pars, probs = c(0.025, 0.975))$summary,
			file=paste0(outfile.base,'190327_participation_modelage6_summary.csv')
	)
	
	# extract samples for unique strata levels
	nprior		<- 1e3
	fit.e		<- extract(fit.par)
	set.seed(42)
	tmp			<- sample(length(fit.e$a), nprior)
	dps			<- dp[,	
			{
				z<- with(fit.e, a + comm[,COMM_NUM_B] + male * MALE + 
								trading*COMM_TYPE_T  + fishing*COMM_TYPE_F +
								inmigrant*INMIGRANT + inmigrant_age1*INMIGRANT*AGE1 + inmigrant_age2*INMIGRANT*AGE2+
								male_age1*AGE1*MALE + female_age1*AGE1*(1-MALE) +
								male_age2*AGE2*MALE + female_age2*AGE2*(1-MALE) +
								age3*AGE3 + age4*AGE4 + age5*AGE5 + age6*AGE6)
				list(SAMPLE=1:nprior, ETA=as.numeric(z[tmp]), DISPERSION=as.numeric(fit.e$dispersion[tmp]))
			},	
			by=c('CATEGORY','TRIAL','SUC')]
	dps[, P:= exp(ETA)/(1+exp(ETA))]
	dps[, LP:= dbbinom(SUC, TRIAL, alpha= P*DISPERSION, beta= (1-P)*DISPERSION, log=TRUE)]	
	set(dps, NULL, c('TRIAL','SUC','ETA','DISPERSION'), NULL)
	
	save(dp, dps, fit.par, file=paste0(outfile.base,'190327_participation_modelage6_samples.rda'))
}


Rakai190327.sequencing.differences.binomialmodel1<- function()
{
	require(data.table)
	require(rstan)
	require(bayesplot)
	
	opt.exclude.onART.from.denominator	<- 1
	indir		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis'
	infile.data	<- file.path(indir,"180322_sampling_by_gender_age.rda")
	infile.sequencing.stan.model <- file.path(indir,"180322_glm_sequencing.stan")
	outfile.base<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	
	#	set up variables for STAN
	load(infile.data)
	des		<- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID, SEX, INMIGRANT))		
	
	# remove individuals in the communities where no sequences were obtained successfully
	tmp	<- des[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
	des	<- merge(des, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')	
	
	# set no INMIGRANT to 0
	set(des, des[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
	
	# age group
	des[, AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,25,35,60), labels = c("15-24", "25-34","35+"))]
	
	#	binarize age, sex
	des[, AGE_YOUNG:= as.integer(AGE_AT_MID_C=='15-24')]
	des[, AGE_MID:= as.integer(AGE_AT_MID_C=='25-34')]
	des[, MALE:= as.integer(SEX=='M')]
	
	#	binarize community type
	des[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	des[, COMM_TYPE_T:= as.integer(substr(COMM_NUM_A,1,1)=='t')]
	
	#	vanilla community IDs
	tmp		<- data.table(COMM_NUM_A= sort(unique(des$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(des$COMM_NUM_A)))
	des		<- merge(des, tmp, by='COMM_NUM_A')
	
	#	aggregate sequencing
	des		<- subset(des, HIV_1517==1)
	if(opt.exclude.onART.from.denominator)
	{
		stopifnot(des[, !any(is.na(SELFREPORTART_AT_FIRST_VISIT))])
		des	<- subset(des, SELFREPORTART_AT_FIRST_VISIT!=1 | (SELFREPORTART_AT_FIRST_VISIT==1 & MIN_PNG_OUTPUT==1))		
	}
	stopifnot( all(unique(sort(des$COMM_NUM_B))==1:36) )
	ds		<- des[, list(TRIAL=length(PERM_ID), SUC=length(which(MIN_PNG_OUTPUT==1))), by=c('AGE_AT_MID_C','SEX','AGE_YOUNG','AGE_MID','INMIGRANT','MALE','COMM_TYPE_F','COMM_TYPE_T','COMM_NUM_A','COMM_NUM_B')]
	ds[, CATEGORY:= paste0(COMM_NUM_A,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
	
	#	run STAN 
	tmp			<- as.list(subset(ds,select=c('COMM_NUM_B','TRIAL','SUC','MALE','AGE_YOUNG','AGE_MID','INMIGRANT','COMM_TYPE_F','COMM_TYPE_T')))
	tmp$N		<- nrow(ds)
	tmp$N_COMM	<- length(unique(ds$COMM_NUM_A))
	fit.seq 	<- stan(	file = infile.sequencing.stan.model, 
			data = tmp, 
			iter = 10e3,
			warmup = 5e2,
			cores = 1,
			chains = 1,
			init = list(list(a=0, comm=rep(0,36), sig_comm=1, male=0, age_mid=0, age_young=0, inmigrant=0, trading=0, fishing=0)))
	
	# assess convergence
	fit.pars	<- c('a','comm','sig_comm','male','age_mid','age_young','inmigrant','trading','fishing')
	any(rhat(fit.seq, pars=fit.pars)>1.02)
	any(neff_ratio(fit.seq, pars=fit.pars) * 9.5e3 < 500)
	po              <- as.matrix(fit.seq) 
	po              <- po[, colnames(po)[!grepl('p_suc_logit|lp__',colnames(po))]]	
	p   <- mcmc_areas(po, pars=colnames(po), prob = 0.95) + 
			geom_vline(xintercept = 0)
	pdf(file=paste0(outfile.base,'190327_sequencing_model_marginalposteriors.pdf'), w=7, h=12)
	p
	dev.off()
	p   <- mcmc_trace(po, pars=colnames(po), facet_args = list(ncol = 1)) 
	pdf(file=paste0(outfile.base,'190327_sequencing_model_marginaltraces.pdf'), w=7, h=100)
	p
	dev.off()
	
	# summary csv file
	write.csv(  summary(fit.seq, pars= fit.pars, probs = c(0.025, 0.975))$summary,
				file=paste0(outfile.base,'190327_sequencing_model_summary.csv')
				)	
	
	# extract samples for unique strata levels
	nprior		<- 500
	fit.e		<- extract(fit.seq)
	set.seed(42)
	tmp			<- sample(length(fit.e$a), nprior)
	dss			<- ds[,	
			{
				z<- with(fit.e, a + comm[,COMM_NUM_B] + male * MALE + 
								trading*COMM_TYPE_T  + fishing*COMM_TYPE_F +
								inmigrant*INMIGRANT + age_young*AGE_YOUNG + age_mid*AGE_MID)
				list(SAMPLE=1:nprior, ETA=as.numeric(z[tmp]))
			},	
			by=c('CATEGORY','TRIAL','SUC')]
	dss[, P:= exp(ETA)/(1+exp(ETA))]
	dss[, LP:= dbinom(SUC, TRIAL, prob=P, log=TRUE)]	
	set(dss, NULL, c('TRIAL','SUC','ETA'), NULL)
	
	save(ds, dss, fit.seq, file=paste0(outfile.base,'190327_sequencing_model_samples.rda'))
}
Rakai190327.sequencing.differences.binomialmodel2<- function()
{
	require(data.table)
	require(rstan)
	require(bayesplot)
	
	opt.exclude.onART.from.denominator	<- 1
	indir		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis'
	infile.data	<- file.path(indir,"180322_sampling_by_gender_age.rda")
	infile.sequencing.stan.model <- file.path(indir,"180322_glm_sequencing_age6.stan")
	outfile.base<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"
	
	#	set up variables for STAN
	load(infile.data)
	des		<- subset(de, select=c(PARTICIPATED, HIV_1517, SELFREPORTART_AT_FIRST_VISIT, MIN_PNG_OUTPUT, PERM_ID, COMM_NUM_A, AGE_AT_MID, SEX, INMIGRANT))		
	
	# remove individuals in the communities where no sequences were obtained successfully
	tmp	<- des[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
	des	<- merge(des, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')	
	
	# set no INMIGRANT to 0
	set(des, des[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
	
	# age group
	des[,AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("20-", "20-24","25-29","30-34","35-39","40-44", "45+"))]
		
	#	binarize age, sex
	des[, AGE1:= as.integer(AGE_AT_MID_C=='20-')]
	des[, AGE2:= as.integer(AGE_AT_MID_C=='20-24')]
	des[, AGE3:= as.integer(AGE_AT_MID_C=='25-29')]
	des[, AGE4:= as.integer(AGE_AT_MID_C=='30-34')]
	des[, AGE5:= as.integer(AGE_AT_MID_C=='35-39')]
	des[, AGE6:= as.integer(AGE_AT_MID_C=='40-44')]
	des[, MALE:= as.integer(SEX=='M')]
	
	#	binarize community type
	des[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	des[, COMM_TYPE_T:= as.integer(substr(COMM_NUM_A,1,1)=='t')]
	
	#	vanilla community IDs
	tmp		<- data.table(COMM_NUM_A= sort(unique(des$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(des$COMM_NUM_A)))
	des		<- merge(des, tmp, by='COMM_NUM_A')
	
	#	aggregate sequencing
	des		<- subset(des, HIV_1517==1)
	if(opt.exclude.onART.from.denominator)
	{
		stopifnot(des[, !any(is.na(SELFREPORTART_AT_FIRST_VISIT))])
		des	<- subset(des, SELFREPORTART_AT_FIRST_VISIT!=1 | (SELFREPORTART_AT_FIRST_VISIT==1 & MIN_PNG_OUTPUT==1))		
	}
	stopifnot( all(unique(sort(des$COMM_NUM_B))==1:36) )
	ds		<- des[, list(TRIAL=length(PERM_ID), SUC=length(which(MIN_PNG_OUTPUT==1))), by=c('AGE_AT_MID_C','SEX','AGE1','AGE2','AGE3','AGE4','AGE5','AGE6','INMIGRANT','MALE','COMM_TYPE_F','COMM_TYPE_T','COMM_NUM_A','COMM_NUM_B')]
	ds[, CATEGORY:= paste0(COMM_NUM_A,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
	
	#	run STAN 
	tmp			<- as.list(subset(ds,select=c('COMM_NUM_B','TRIAL','SUC','MALE','AGE1','AGE2','AGE3','AGE4','AGE5','AGE6','INMIGRANT','COMM_TYPE_F','COMM_TYPE_T')))
	tmp$N		<- nrow(ds)
	tmp$N_COMM	<- length(unique(ds$COMM_NUM_A))
	fit.seq 	<- stan(	file = infile.sequencing.stan.model, 
			data = tmp, 
			iter = 10e3,
			warmup = 5e2,
			cores = 1,
			chains = 1,
			init = list(list(a=0, comm=rep(0,36), sig_comm=1, male=0, age1=0, age2=0, age3=0, age4=0, age5=0, age6=0, inmigrant=0, trading=0, fishing=0)))
	
	# assess convergence
	fit.pars	<- c('a','comm','sig_comm','male','age1','age2','age3','age4','age5','age6','inmigrant','trading','fishing')
	any(rhat(fit.seq, pars=fit.pars)>1.02)
	any(neff_ratio(fit.seq, pars=fit.pars) * 9.5e3 < 500)
	po              <- as.matrix(fit.seq) 
	po              <- po[, colnames(po)[!grepl('p_suc_logit|lp__',colnames(po))]]	
	p   <- mcmc_areas(po, pars=colnames(po), prob = 0.95) + 
			geom_vline(xintercept = 0)
	pdf(file=paste0(outfile.base,'190327_sequencing_modelage6_marginalposteriors.pdf'), w=7, h=12)
	p
	dev.off()
	p   <- mcmc_trace(po, pars=colnames(po), facet_args = list(ncol = 1)) 
	pdf(file=paste0(outfile.base,'190327_sequencing_modelage6_marginaltraces.pdf'), w=7, h=100)
	p
	dev.off()
	
	# summary csv file
	write.csv(  summary(fit.seq, pars= fit.pars, probs = c(0.025, 0.975))$summary,
				file=paste0(outfile.base,'190327_sequencing_modelage6_summary.csv')
				)	
	
	# extract samples for unique strata levels
	nprior		<- 1e3
	fit.e		<- extract(fit.seq)
	set.seed(42)
	tmp			<- sample(length(fit.e$a), nprior)
	dss			<- ds[,	
			{
				z<- with(fit.e, a + comm[,COMM_NUM_B] + male * MALE + 
								trading*COMM_TYPE_T  + fishing*COMM_TYPE_F +
								inmigrant*INMIGRANT + 
								age1*AGE1 + age2*AGE2 + age3*AGE3 + age4*AGE4 + age5*AGE5 + age6*AGE6)
				list(SAMPLE=1:nprior, ETA=as.numeric(z[tmp]))
			},	
			by=c('CATEGORY','TRIAL','SUC')]
	dss[, P:= exp(ETA)/(1+exp(ETA))]
	dss[, LP:= dbinom(SUC, TRIAL, prob=P, log=TRUE)]	
	set(dss, NULL, c('TRIAL','SUC','ETA'), NULL)
	
	save(ds, dss, fit.seq, file=paste0(outfile.base,'190327_sequencing_modelage6_samples.rda'))
}



Rakai190327.RCCStransmissionflows.inference.age3model<- function(infile.inference=NULL, opt=NULL)
{
	require(data.table)	
	require(TransSubpopulation)
	#require(Boom)	
	#require(gtools)	
	
	#logistic<- function(x) 1/(1+exp(-x))
	
	if(is.null(opt))
	{
		opt									<- list()
		opt$adjust.sequencing.bias			<- 1
		opt$adjust.participation.bias		<- 1
		opt$exclude.onART.from.denominator	<- 1
		opt$set.missing.migloc.to.inland	<- 0
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland		
	}
	if(is.null(infile.inference))
	{
		indir				<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run"
		infile.inference	<- file.path(indir,"todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda")
	}	
	indir					<- dirname(infile.inference)	
	outfile.base			<- gsub('data_with_inmigrants.rda','',infile.inference)
	load(infile.inference)
	
	#
	#	prepare data on observed transmission flows
	#
	
	#	subset to variables needed, using RTR3
	rtr	<- subset(rtr3, select=c(	'PAIRID','TR_RID','TR_COMM_NUM','TR_COMM_NUM_A','TR_COMM_NUM_A_MIG',
							'TR_SEX','TR_BIRTHDATE','TR_COMM_TYPE','TR_INMIG_LOC','TR_INMIGRATE_1YR','TR_INMIGRATE_2YR',
							'REC_RID','REC_COMM_NUM','REC_COMM_NUM_A',
							'REC_SEX','REC_BIRTHDATE','REC_COMM_TYPE','REC_INMIGRATE_1YR','REC_INMIGRATE_2YR'))
	
	# add age 
	rtr[,TR_AGE_AT_MID:=2013.25-TR_BIRTHDATE]
	rtr[,REC_AGE_AT_MID:=2013.25-REC_BIRTHDATE]
	# impute age
	tmp	<- which(is.na(rtr$TR_AGE_AT_MID))
	set(rtr, tmp, 'TR_AGE_AT_MID', mean(rtr$TR_AGE_AT_MID[which(!is.na(rtr$TR_AGE_AT_MID))]) )
	tmp	<- which(is.na(rtr$REC_AGE_AT_MID))
	set(rtr, tmp, 'REC_AGE_AT_MID', mean(rtr$REC_AGE_AT_MID[which(!is.na(rtr$REC_AGE_AT_MID))]) )
	# fixup from latest surveillance data
	#de[RID=="C036808",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	#rtr[TR_RID=="C036808",c('TR_COMM_NUM_A','TR_INMIGRANT','TR_AGE_AT_MID','TR_SEX')]
	set(rtr, rtr[,which(TR_RID=="C036808")], 'TR_AGE_AT_MID', 39.946)	
	#de[RID=="D069722",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	#rtr[TR_RID=="D069722",c('TR_COMM_NUM_A','TR_INMIGRANT','TR_AGE_AT_MID','TR_SEX')]
	set(rtr, rtr[,which(TR_RID=="D069722")], 'TR_INMIGRANT', 0L)
	#de[RID=="G036802",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	#rtr[REC_RID=="G036802",c('REC_COMM_NUM_A','REC_INMIGRANT','REC_AGE_AT_MID','REC_SEX')]
	set(rtr, rtr[,which(REC_RID=="G036802")], 'REC_AGE_AT_MID',	44.946)	
	#de[RID=="H103745",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	#rtr[REC_RID=="H103745",c('REC_COMM_NUM_A','REC_INMIGRANT','REC_AGE_AT_MID','REC_SEX')]
	set(rtr, rtr[, which(REC_RID=="H103745")], 'REC_AGE_AT_MID', 20.42)	
	set(rtr, rtr[, which(REC_RID=="C121534")],'REC_AGE_AT_MID', 28.549)
	#	stratify age
	rtr[, TR_AGE_AT_MID_C:= as.character(cut(TR_AGE_AT_MID, breaks=c(10,25,35,65), labels=c('15-24','25-34','35+'), right=FALSE))]
	rtr[, REC_AGE_AT_MID_C:= as.character(cut(REC_AGE_AT_MID, breaks=c(10,25,35,65), labels=c('15-24','25-34','35+'), right=FALSE))]
	stopifnot( nrow(subset(rtr, is.na(TR_AGE_AT_MID_C)))==0 )
	stopifnot( nrow(subset(rtr, is.na(REC_AGE_AT_MID_C)))==0 )
	
		
	# inmigrant status
	setnames(rtr, 'TR_INMIGRATE_2YR', 'TR_INMIGRATE')
	setnames(rtr, 'REC_INMIGRATE_2YR', 'REC_INMIGRATE')
	rtr[, TR_INMIGRANT:= as.integer(TR_INMIGRATE!='resident')]
	rtr[, REC_INMIGRANT:= as.integer(grepl('inmigrant',REC_INMIGRATE))]
	set(rtr, NULL, 'TR_COMM_NUM_A_MIG', rtr[, gsub('[0-9]+','',TR_COMM_NUM_A_MIG)])
	#	set unknown origin to either fishing or inland
	tmp	<- rtr[, which(TR_INMIGRATE=='inmigrant_from_unknown')]
	if(opt$set.missing.migloc.to.inland)
	{
		set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_inland')
		set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'imig')
	}		
	if(opt$set.missing.migloc.to.fishing)
	{
		set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_fisherfolk')
		set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'fmig')
	}
			
	#	build category to match with sampling data tables 
	rtr[, REC_SAMPLING_CATEGORY:= paste0(REC_COMM_NUM_A,':',REC_SEX,':',REC_AGE_AT_MID_C,':',REC_INMIGRANT)]
	rtr[, TR_SAMPLING_CATEGORY:= paste0(TR_COMM_NUM_A,':',TR_SEX,':',TR_AGE_AT_MID_C,':',TR_INMIGRANT)]
	#	build transmission flow category 
	rtr[, REC_TRM_CATEGORY:= paste0(REC_COMM_NUM_A,':',REC_SEX,':',REC_AGE_AT_MID_C,':',REC_INMIGRANT)]
	rtr[, TR_TRM_CATEGORY:= paste0(TR_COMM_NUM_A_MIG,':',TR_SEX,':',TR_AGE_AT_MID_C,':',TR_INMIGRANT)]
	
	#
	#	calculate observed number of transmissions
	#
	dobs	<- rtr[, list( TRM_OBS=length(unique(PAIRID))), by=c('TR_TRM_CATEGORY','REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY')]
	setkey(dobs, TR_TRM_CATEGORY, REC_TRM_CATEGORY )	
	dobs[, TRM_CAT_PAIR_ID:= seq_len(nrow(dobs))]
	
	#	load samples from prior
	infile.participation.prior.samples	<- file.path(indir,"190327_participation_model_samples.rda")
	infile.sequencing.prior.samples		<- file.path(indir,"190327_sequencing_model_samples.rda")
	load(infile.participation.prior.samples)
	load(infile.sequencing.prior.samples)
	dprior	<- merge(dps, dss, by=c('CATEGORY','SAMPLE'))
	dprior[, P:= P.x*P.y]		# multiply participation and sequencing probabilities 
	dprior[, LP:= LP.x+LP.y]	# add log posterior densities
	set(dprior, NULL, c('P.x','LP.x','P.y','LP.y'), NULL)
	setnames(dprior, 'CATEGORY', 'SAMPLING_CATEGORY')
	
	#	run MCMC
	control	<- list(seed=42, mcmc.n=238000*100, verbose=0, outfile=paste0(outfile.base,"SAMCMCv190327_mcmc.rda"))
	source.attribution.mcmc(dobs, dprior, control=control)
	stop()
	
	#	diagnostics
	load(control$outfile)
	
	#	acceptance rate per site
	da	<- subset(mc$it.info, !is.na(PAR_ID) & PAR_ID>0)[, list(ACC_RATE=mean(ACCEPT)), by='PAR_ID']
	setnames(da, 'PAR_ID', 'UPDATE_ID')
	tmp	<- mc$dl[, list(N_TRM_CAT_PAIRS=length(TRM_CAT_PAIR_ID)), by='UPDATE_ID']
	da	<- merge(da, tmp, by='UPDATE_ID')
	ggplot(da, aes(x=N_TRM_CAT_PAIRS, y=ACC_RATE)) + geom_point()
	
	#	traces
	tmp	<- mc$pars$PI
	colnames(tmp)<- paste0('PI-', 1:ncol(tmp))
	p	<- mcmc_trace(tmp, pars=colnames(tmp), facet_args = list(ncol = 1))
	pdf(file=paste0(outfile.base,'samodel_marginaltraces.pdf'), w=7, h=300)
	p
	dev.off()
	
	#	effective sample size on target parameter
	require(coda)
	tmp	<- mcmc(mc$pars$PI)
	de	<- data.table(PI= 1:ncol(tmp), NEFF=as.numeric(effectiveSize(tmp)))
	ggplot(de, aes(x=PI, y=NEFF)) + geom_point()	
	
	p   <- mcmc_trace(po, pars=colnames(po), facet_args = list(ncol = 1)) 
	pdf(file=paste0(outfile.base,'190327_participation_model_marginaltraces.pdf'), w=7, h=100)
	p
	dev.off()
}
Rakai190327.RCCStransmissionflows.inference.age6model<- function(infile.inference=NULL, opt=NULL)
{
	require(data.table)	
	#require(Boom)	
	#require(gtools)	
	
	#logistic<- function(x) 1/(1+exp(-x))
	
	if(is.null(opt))
	{
		opt									<- list()
		opt$adjust.sequencing.bias			<- 1
		opt$adjust.participation.bias		<- 1
		opt$exclude.onART.from.denominator	<- 1
		opt$set.missing.migloc.to.inland	<- 0
		opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland		
	}
	if(is.null(infile.inference))
	{
		infile.inference	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda"
	}	
	indir					<- dirname(infile.inference)	
	outfile.base			<- gsub('data_with_inmigrants.rda','',infile.inference)
	load(infile.inference)
	
	#
	#	prepare data on observed transmission flows
	#
	
	#	subset to variables needed, using RTR3
	rtr	<- subset(rtr3, select=c(	'PAIRID','TR_RID','TR_COMM_NUM','TR_COMM_NUM_A','TR_COMM_NUM_A_MIG',
					'TR_SEX','TR_BIRTHDATE','TR_COMM_TYPE','TR_INMIG_LOC','TR_INMIGRATE_1YR','TR_INMIGRATE_2YR',
					'REC_RID','REC_COMM_NUM','REC_COMM_NUM_A',
					'REC_SEX','REC_BIRTHDATE','REC_COMM_TYPE','REC_INMIGRATE_1YR','REC_INMIGRATE_2YR'))
	
	# add age 
	rtr[,TR_AGE_AT_MID:=2013.25-TR_BIRTHDATE]
	rtr[,REC_AGE_AT_MID:=2013.25-REC_BIRTHDATE]
	# impute age
	tmp	<- which(is.na(rtr$TR_AGE_AT_MID))
	set(rtr, tmp, 'TR_AGE_AT_MID', mean(rtr$TR_AGE_AT_MID[which(!is.na(rtr$TR_AGE_AT_MID))]) )
	tmp	<- which(is.na(rtr$REC_AGE_AT_MID))
	set(rtr, tmp, 'REC_AGE_AT_MID', mean(rtr$REC_AGE_AT_MID[which(!is.na(rtr$REC_AGE_AT_MID))]) )
	# fixup from latest surveillance data
	#de[RID=="C036808",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	#rtr[TR_RID=="C036808",c('TR_COMM_NUM_A','TR_INMIGRANT','TR_AGE_AT_MID','TR_SEX')]
	set(rtr, rtr[,which(TR_RID=="C036808")], 'TR_AGE_AT_MID', 39.946)	
	#de[RID=="D069722",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	#rtr[TR_RID=="D069722",c('TR_COMM_NUM_A','TR_INMIGRANT','TR_AGE_AT_MID','TR_SEX')]
	set(rtr, rtr[,which(TR_RID=="D069722")], 'TR_INMIGRATE_2YR', 'resident')
	#de[RID=="G036802",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	#rtr[REC_RID=="G036802",c('REC_COMM_NUM_A','REC_INMIGRANT','REC_AGE_AT_MID','REC_SEX')]
	set(rtr, rtr[,which(REC_RID=="G036802")], 'REC_AGE_AT_MID',	44.946)	
	#de[RID=="H103745",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	#rtr[REC_RID=="H103745",c('REC_COMM_NUM_A','REC_INMIGRANT','REC_AGE_AT_MID','REC_SEX')]
	set(rtr, rtr[, which(REC_RID=="H103745")], 'REC_AGE_AT_MID', 20.42)	
	set(rtr, rtr[, which(REC_RID=="C121534")],'REC_AGE_AT_MID', 28.549)
	#	stratify age
	rtr[, TR_AGE_AT_MID_C:= cut(TR_AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("20-", "20-24","25-29","30-34","35-39","40-44", "45+"))]
	rtr[, REC_AGE_AT_MID_C:= cut(REC_AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("20-", "20-24","25-29","30-34","35-39","40-44", "45+"))]
	stopifnot( nrow(subset(rtr, is.na(TR_AGE_AT_MID_C)))==0 )
	stopifnot( nrow(subset(rtr, is.na(REC_AGE_AT_MID_C)))==0 )
	
	
	# inmigrant status
	setnames(rtr, 'TR_INMIGRATE_2YR', 'TR_INMIGRATE')
	setnames(rtr, 'REC_INMIGRATE_2YR', 'REC_INMIGRATE')
	rtr[, TR_INMIGRANT:= as.integer(TR_INMIGRATE!='resident')]
	rtr[, REC_INMIGRANT:= as.integer(grepl('inmigrant',REC_INMIGRATE))]
	set(rtr, NULL, 'TR_COMM_NUM_A_MIG', rtr[, gsub('[0-9]+','',TR_COMM_NUM_A_MIG)])
	#	set unknown origin to either fishing or inland
	tmp	<- rtr[, which(TR_INMIGRATE=='inmigrant_from_unknown')]
	if(opt$set.missing.migloc.to.inland)
	{
		set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_inland')
		set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'imig')
	}		
	if(opt$set.missing.migloc.to.fishing)
	{
		set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_fisherfolk')
		set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'fmig')
	}
	
	#	build category to match with sampling data tables 
	rtr[, REC_SAMPLING_CATEGORY:= paste0(REC_COMM_NUM_A,':',REC_SEX,':',REC_AGE_AT_MID_C,':',REC_INMIGRANT)]
	rtr[, TR_SAMPLING_CATEGORY:= paste0(TR_COMM_NUM_A,':',TR_SEX,':',TR_AGE_AT_MID_C,':',TR_INMIGRANT)]
	#	build transmission flow category 
	rtr[, REC_TRM_CATEGORY:= paste0(REC_COMM_NUM_A,':',REC_SEX,':',REC_AGE_AT_MID_C,':',REC_INMIGRANT)]
	rtr[, TR_TRM_CATEGORY:= paste0(TR_COMM_NUM_A_MIG,':',TR_SEX,':',TR_AGE_AT_MID_C,':',TR_INMIGRANT)]
	
	#
	#	calculate observed number of transmissions
	#
	dobs	<- rtr[, list( TRM_OBS=length(unique(PAIRID))), by=c('TR_TRM_CATEGORY','REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY')]
	setkey(dobs, TR_TRM_CATEGORY, REC_TRM_CATEGORY )	
	dobs[, TRM_CAT_PAIR_ID:= seq_len(nrow(dobs))]
	
	#	load samples from prior
	infile.participation.prior.samples	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190327_participation_modelage6_samples.rda"
	infile.sequencing.prior.samples		<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190327_sequencing_modelage6_samples.rda"
	load(infile.participation.prior.samples)
	load(infile.sequencing.prior.samples)
	dprior	<- merge(dps, dss, by=c('CATEGORY','SAMPLE'))
	dprior[, P:= P.x*P.y]		# multiply participation and sequencing probabilities 
	dprior[, LP:= LP.x+LP.y]	# add log posterior densities
	set(dprior, NULL, c('P.x','LP.x','P.y','LP.y'), NULL)
	setnames(dprior, 'CATEGORY', 'SAMPLING_CATEGORY')
	
	#	run MCMC
	control	<- list(seed=42, mcmc.n=238000, verbose=0, outfile="~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/190327_SAMCMCv190327_age6_mcmc.rda")
	mcmc.core.inference(dobs, dprior, control=control)
	
	#	diagnostics
	load(control$outfile)
	
	#	acceptance rate per site
	da	<- subset(mc$it.info, !is.na(PAR_ID) & PAR_ID>0)[, list(ACC_RATE=mean(ACCEPT)), by='PAR_ID']
	setnames(da, 'PAR_ID', 'UPDATE_ID')
	tmp	<- mc$dl[, list(N_TRM_CAT_PAIRS=length(TRM_CAT_PAIR_ID)), by='UPDATE_ID']
	da	<- merge(da, tmp, by='UPDATE_ID')
	ggplot(da, aes(x=N_TRM_CAT_PAIRS, y=ACC_RATE)) + geom_point()
	
	#	traces
	tmp	<- mc$pars$PI
	colnames(tmp)<- paste0('PI-', 1:ncol(tmp))
	p	<- mcmc_trace(tmp, pars=colnames(tmp), facet_args = list(ncol = 1))
	pdf(file=paste0(outfile.base,'samodel_marginaltraces.pdf'), w=7, h=300)
	p
	dev.off()
	
	#	effective sample size on target parameter
	require(coda)
	tmp	<- mcmc(mc$pars$PI)
	de	<- data.table(PI= 1:ncol(tmp), NEFF=as.numeric(effectiveSize(tmp)))
	ggplot(de, aes(x=PI, y=NEFF)) + geom_point()	
	
	p   <- mcmc_trace(po, pars=colnames(po), facet_args = list(ncol = 1)) 
	pdf(file=paste0(outfile.base,'190327_participation_model_marginaltraces.pdf'), w=7, h=100)
	p
	dev.off()
}

xx<- function()
{
	########################################## create dg with stratification variables  ################################################
	##########################################  and induced stan regression covaraites ################################################
	##########################################  the number of success (SEQ) and trials (HIV)################################################
	library(data.table)
	opt.exclude.onART.from.denominator	<- 1
	indir	<- '~/Desktop/data'
	indir	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis'
	load(file.path(indir,"180322_sampling_by_gender_age.rda"))
	# set unknown INMIGRANT with 0
	set(de, de[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L) # MAY LEAD TO BIAS
	
	# # check na
	# de[, lapply(.SD, function(x) sum(is.na(x))), .SDcols = 1:length(colnames(de))]
	
	# take infected individuals
	de	<- subset(de, HIV_1517==1)
	if(opt.exclude.onART.from.denominator)
	{
	  stopifnot(de[, !any(is.na(SELFREPORTART_AT_FIRST_VISIT))])
	  #select not on art or ? NOT CLEAR ABOUT REASON
	  de	<- subset(de, SELFREPORTART_AT_FIRST_VISIT!=1 | (SELFREPORTART_AT_FIRST_VISIT==1 & MIN_PNG_OUTPUT==1))
	}
	
	# remove individuals in the communities where no sequences were obtained successfully
	tmp	<- de[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
	de	<- merge(de, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')
	
	# age group
	de[,AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("20-", "20-24","25-29","30-34","35-39","40-44", "45+"))]
	
	#	set up variables for STAN
	dg		<- subset(de, select=c(MIN_PNG_OUTPUT, COMM_NUM_A, AGE_AT_MID_C, SEX,INMIGRANT))
	
	#	binarize age, sex
	dg[, AGE1:= as.integer(AGE_AT_MID_C=='20-')]
	dg[, AGE2:= as.integer(AGE_AT_MID_C=='20-24')]
	dg[, AGE3:= as.integer(AGE_AT_MID_C=='25-29')]
	dg[, AGE4:= as.integer(AGE_AT_MID_C=='30-34')]
	dg[, AGE5:= as.integer(AGE_AT_MID_C=='35-39')]
	dg[, AGE6:= as.integer(AGE_AT_MID_C=='40-44')]
	dg[, MALE:= as.integer(SEX=='M')]
	#	binarize community type
	dg[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	dg[, COMM_TYPE_T:= as.integer(substr(COMM_NUM_A,1,1)=='t')]
	#	vanilla community IDs
	tmp		<- data.table(COMM_NUM_A= sort(unique(dg$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(dg$COMM_NUM_A)))
	dg		<- merge(dg, tmp, by='COMM_NUM_A')
	
	dg		<- dg[, list(TRIAL=length(MIN_PNG_OUTPUT), SUC=length(which(MIN_PNG_OUTPUT==1))),
	          by=c('AGE_AT_MID_C','SEX','AGE1','AGE2','AGE3','AGE4','AGE5','AGE6','INMIGRANT','MALE','COMM_TYPE_F','COMM_TYPE_T','COMM_NUM_A','COMM_NUM_B')]
	
	dg[,CATEGORY:=paste0(COMM_NUM_A,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
	########################################## create dp with stratification variables  ################################################
	##########################################  and induced stan regression covaraites ################################################
	##########################################  the number of success (PART) and trials (ELIG)################################################
	load(file.path(indir,"180322_sampling_by_gender_age.rda"))
	# set unknown INMIGRANT with 0
	set(de, de[, which(is.na(INMIGRANT))], 'INMIGRANT', 0L)
	
	# remove individuals in the communities where no sequences were obtained successfully
	tmp	<- de[, list(COMM_ANY_MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT, na.rm=TRUE)), by='COMM_NUM_A']
	de	<- merge(de, subset(tmp, COMM_ANY_MIN_PNG_OUTPUT>0, COMM_NUM_A), by='COMM_NUM_A')
	
	# age group
	de[,AGE_AT_MID_C:= cut(AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("20-", "20-24","25-29","30-34","35-39","40-44", "45+"))]
	
	#	set up variables for STAN
	dp		<- subset(de, select=c(PARTICIPATED, PERM_ID, COMM_NUM_A, AGE_AT_MID_C, SEX,INMIGRANT))
	#	binarize age, sex
	dp[, AGE1:= as.integer(AGE_AT_MID_C=='20-')]
	dp[, AGE2:= as.integer(AGE_AT_MID_C=='20-24')]
	dp[, AGE3:= as.integer(AGE_AT_MID_C=='25-29')]
	dp[, AGE4:= as.integer(AGE_AT_MID_C=='30-34')]
	dp[, AGE5:= as.integer(AGE_AT_MID_C=='35-39')]
	dp[, AGE6:= as.integer(AGE_AT_MID_C=='40-44')]
	dp[, MALE:= as.integer(SEX=='M')]
	#	binarize community type
	dp[, COMM_TYPE_F:= as.integer(substr(COMM_NUM_A,1,1)=='f')]
	dp[, COMM_TYPE_T:= as.integer(substr(COMM_NUM_A,1,1)=='t')]
	#	vanilla community IDs
	tmp		<- data.table(COMM_NUM_A= sort(unique(dp$COMM_NUM_A)), COMM_NUM_B= seq_along(unique(dp$COMM_NUM_A)))
	dp		<- merge(dp, tmp, by='COMM_NUM_A')
	
	dp		<- dp[, list(TRIAL=length(PERM_ID), SUC=length(which(PARTICIPATED==1))), by=c('AGE_AT_MID_C','SEX','AGE1','AGE2','AGE3','AGE4','AGE5','AGE6','INMIGRANT','MALE','COMM_TYPE_F','COMM_TYPE_T','COMM_NUM_A','COMM_NUM_B')]
	
	dp[,CATEGORY:=paste0(COMM_NUM_A,':',SEX,':',AGE_AT_MID_C,':',INMIGRANT)]
	
	df.sampling<-list()
	df.sampling[[1]]<-dg
	df.sampling[[2]]<-dp
	
	########################################## create dc with stratification variables ################################################
	########################################## observed transmission counts, count id ################################################
	##########################################      tr_category and rec_category      ################################################
	indir2	<- '~/Desktop/data'
	indir2	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run'
	load(file.path(indir2,"todi_pairs_181006_cl25_d50_prior23_min30_phylogeography_data_with_inmigrants.rda"))
	
	require(data.table)
	# add age into rtr2
	rtr2[,TR_AGE_AT_MID:=2013.25-TR_BIRTHDATE]
	rtr2[,REC_AGE_AT_MID:=2013.25-REC_BIRTHDATE]
	
	# impute age
	rtr2$TR_AGE_AT_MID[which(is.na(rtr2$TR_AGE_AT_MID))]<-mean(rtr2$TR_AGE_AT_MID[which(!is.na(rtr2$TR_AGE_AT_MID))])
	rtr2$REC_AGE_AT_MID[which(is.na(rtr2$REC_AGE_AT_MID))]<-mean(rtr2$REC_AGE_AT_MID[which(!is.na(rtr2$REC_AGE_AT_MID))])
	# MAY LEAD TO BIAS
	
	# inmigrant
	rtr2[,REC_INMIGRANT:=as.numeric(REC_INMIGRATE_2YR!='resident')]
	rtr2[,TR_INMIGRANT:=as.numeric(TR_INMIGRATE_2YR!='resident')]
	
	# some rtr2 doesn't have sampling information
	# correct rtr2 TR_CATEGORY
	#"C036808" "D069722"
	de[RID=="C036808",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	rtr2[TR_RID=="C036808",c('TR_COMM_NUM_A','TR_INMIGRANT','TR_AGE_AT_MID','TR_SEX')]
	rtr2[TR_RID=="C036808",]$TR_AGE_AT_MID<-39.946
	
	de[RID=="D069722",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	rtr2[TR_RID=="D069722",c('TR_COMM_NUM_A','TR_INMIGRANT','TR_AGE_AT_MID','TR_SEX')]
	rtr2[TR_RID=="D069722",]$TR_INMIGRANT<-0
	
	#"G036802" "H103745" "C121534"
	de[RID=="G036802",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	rtr2[REC_RID=="G036802",c('REC_COMM_NUM_A','REC_INMIGRANT','REC_AGE_AT_MID','REC_SEX')]
	rtr2[REC_RID=="G036802",]$REC_AGE_AT_MID<-44.946
	
	de[RID=="H103745",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	rtr2[REC_RID=="H103745",c('REC_COMM_NUM_A','REC_INMIGRANT','REC_AGE_AT_MID','REC_SEX')]
	rtr2[REC_RID=="H103745",]$REC_AGE_AT_MID<-20.42
	
	#####NOT SURE
	de[RID=="C121534",c('COMM_NUM_A','INMIGRANT','AGE_AT_MID','SEX')]
	rtr2[REC_RID=="C121534",c('REC_COMM_NUM_A','REC_INMIGRANT','REC_AGE_AT_MID','REC_SEX')]
	rtr2[REC_RID=="C121534",'REC_AGE_AT_MID']<-28.549
	
	de[COMM_NUM_A=='asm' & AGE_AT_MID=='20-24' & SEX=='M',]
	de[COMM_NUM_A=='asm' &  INMIGRANT==1 & SEX=='M',]
	de[COMM_NUM_A=='asm' &  INMIGRANT==1 & SEX=='M'& AGE_AT_MID_C=='25-29',]
	# age groups
	rtr2[,TR_AGE_AT_MID_C:= cut(TR_AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("20-", "20-24","25-29","30-34","35-39","40-44", "45+"))]
	rtr2[,REC_AGE_AT_MID_C:= cut(REC_AGE_AT_MID, breaks = c(0,20,25,30,35,40,45,60), labels = c("20-", "20-24","25-29","30-34","35-39","40-44", "45+"))]
	
	###OLLI: use rtr3 please
	
	# stratification of transmission flows
	dc	<- rtr2[, list(OBS=length(PAIRID)), by=c('TR_COMM_NUM_A','REC_COMM_NUM_A','TR_SEX','REC_SEX','TR_AGE_AT_MID_C','REC_AGE_AT_MID_C','TR_INMIGRANT','REC_INMIGRANT')]
	setkey(dc,TR_COMM_NUM_A,REC_COMM_NUM_A,TR_SEX,REC_SEX,TR_AGE_AT_MID_C,REC_AGE_AT_MID_C,TR_INMIGRANT)
	
	# add count id, tr_category and rec_category
	dc$COUNT_ID<-1:nrow(dc)
	dc[,TR_CATEGORY:=paste0(TR_COMM_NUM_A,':',TR_SEX,':',TR_AGE_AT_MID_C,':',TR_INMIGRANT)]
	dc[,REC_CATEGORY:=paste0(REC_COMM_NUM_A,':',REC_SEX,':',REC_AGE_AT_MID_C,':',REC_INMIGRANT)]
	
	#
	nrow(dc[substring(TR_COMM_NUM_A, 1, 1)=='f' & substring(REC_COMM_NUM_A, 1, 1)=='f',])
	nrow(dc[substring(TR_COMM_NUM_A, 1, 1)!='f' & substring(REC_COMM_NUM_A, 1, 1)!='f',])
	
	## OLLI: the participation model isn t right, it should be mp3
	
	#setwd("~/TransSubpopulation/data_analysis/complex_model1_complete_code")
	setwd(indir2)
	#	iteration=6e7; sample.method='empirical'; glm.model=c('binomial','beta_binomial'); seed=123
	samples.from.GLM.prior(df.sampling,dc,iteration=6e7,sample.method='empirical',glm.model=c('binomial','beta_binomial'),seed=123)
	load("InputsMCMC_GLM6e+07.rda")
	
	mcmc.core.inference(s.dtl.prior=S.DTL.prior,tr.obs=TR.OBS,seed=123)
	
	load("core_inference_SNPIZ_mcmcEachCount_GLM.rda")
	mcmc.core.inference.diagnostics(mc,burnin=2e3)
	
	samples.from.GLM.prior<-function(df.sampling,dc,iteration,sample.method,glm.model,seed=NULL){
	  library(data.table)
	  library(rstan)
	  library(fitdistrplus)
	  set.seed(seed)
	  
	  # take unique subpopulation that appears in transmitters and recipients in dc
	  tr<-dc[, grep("TR_", names(dc)), with = FALSE]
	  setnames(tr, colnames(tr), gsub('TR_','',colnames(tr)))
	  tr<-unique(tr)
	  rec<-dc[, grep("REC_", names(dc)), with = FALSE]
	  setnames(rec, colnames(rec), gsub('REC_','',colnames(rec)))
	  rec<-unique(rec)
	  unique.group<-unique(rbind(tr,rec))
	  # setkey(unique.group,COMM_NUM_A,SEX,AGE_AT_MID_C,INMIGRANT)
	  setkey(unique.group,CATEGORY)
	  unique.group<-merge(subset(unique.group,select = 'CATEGORY'),dg,by='CATEGORY')
	  unique.group[,CATEGORY_ID:=1:nrow(unique.group)]
	  
	  # take unique pairs of subpopulations for transmission flows
	  unique.pairs<-subset(dc,select=c('COUNT_ID','TR_CATEGORY','REC_CATEGORY'))
	  tmp<-subset(unique.group,select = c('CATEGORY','CATEGORY_ID'))
	  setnames(tmp,colnames(tmp),paste0('TR_',colnames(tmp)))
	  unique.pairs<-merge(unique.pairs,tmp,by='TR_CATEGORY')
	  setnames(tmp,colnames(tmp),gsub('TR_','REC_',colnames(tmp)))
	  unique.pairs<-merge(unique.pairs,tmp,by='REC_CATEGORY')
	  setkey(unique.pairs,COUNT_ID)
	  
	  # number of prior samples
	  #nprior<-ceiling(iteration/(nrow(unique.group)+1))+1
	  #	OLLI: much lower OK
	  nprior<- 500
		
	  S.DTL.prior<-list()
	  fit.list<-list() #model
	  S.DTL.prior$SAM_P<-matrix(1,nrow=nprior,ncol=nrow(unique.group)) # sampling rate
	  S.DTL.prior$SAM_P_LOGD<-matrix(0,nrow=nprior,ncol=nrow(unique.group)) # log density
	  S.DTL.prior$EMP_S<-1
	  
	  for (i in 1L:length(df.sampling)){
	    # record priors of sampling rate for the ith source of bias
	    SAM_P<-matrix(NA_real_,nrow=nprior,ncol=nrow(unique.group))
	    SAM_P_LOGD<-matrix(NA_real_,nrow=nprior,ncol=nrow(unique.group))
	    EMP_S<-(sum(df.sampling[[i]]$SUC)/sum(df.sampling[[i]]$TRIAL))^2
	    
	    # predict sampling rates for stratum a using logistic regression
	    ######################## will depend on covariates??? not general???
	    if (glm.model[i]=='binomial'){
	      data.glm<-as.list(subset(df.sampling[[i]],select=c('COMM_NUM_B','TRIAL','AGE1','AGE2','MALE','SUC','AGE3','AGE4','AGE5','AGE6','INMIGRANT','COMM_TYPE_F','COMM_TYPE_T')))
	      data.glm$N<-nrow(df.sampling[[i]])
	      fit <- stan(file = 'glm_age7sex2comm36inmig2bin.stan', data = data.glm, iter=min(max(1e4,5e2+nprior),120000),warmup = 5e2,
	                  cores = 4,chains = 1,init = list(list(a=0, comm=rep(0,36), sig_comm=1, male=0, age1=0, age2=0,age3=0,age4=0,age5=0,age6=0,inmigrant=0,trading=0,fishing=0)))
	    }else if(glm.model[i]=='beta_binomial'){
	      data.glm<-as.list(subset(df.sampling[[i]],select=c('COMM_NUM_B','TRIAL','AGE1','AGE2','MALE','SUC','AGE3','AGE4','AGE5','AGE6','INMIGRANT','COMM_TYPE_F','COMM_TYPE_T')))
	      data.glm$N<-nrow(df.sampling[[i]])
	      fit <- stan(file = 'glm_age7sex2comm36inmig2betabin.stan', data = data.glm, iter=min(max(1e4,5e2+nprior),120000),warmup = 5e2,
	                  cores = 4,chains = 1,init = list(list(a=0, comm=rep(0,36), sig_comm=1, male=0, age1=0, age2=0,dispersion=1,age3=0,age4=0,age5=0,age6=0,inmigrant=0,trading=0,fishing=0)))
	    }
	    # extract samples for the parameters
	    tmp		<- extract(fit,permute=TRUE)
	    if(sample.method=='empirical'){
	      for (j in 1:nrow(unique.group)){
	        # posterior samples of p_suc
	        # tmp2<-tmp$a + tmp$comm[,unique.group$COMM_NUM_B[j]] + tmp$trading*unique.group$COMM_TYPE_T[j] + tmp$fishing*unique.group$COMM_TYPE_F[j] + tmp$male*unique.group$MALE[j] +
	        #   tmp$age1*unique.group$AGE1[j] + tmp$age2*unique.group$AGE2[j] + tmp$age3*unique.group$AGE3[j] +  tmp$age4*unique.group$AGE4[j] +
	        #   tmp$age5*unique.group$AGE5[j]  + tmp$age6*unique.group$AGE6[j] + tmp$inmigrant*unique.group$INMIGRANT[j]
	        tmp2<-tmp$a + tmp$comm[,unique.group$COMM_NUM_B[j]]  + tmp$male*unique.group$MALE[j] +
	          tmp$age1*unique.group$AGE1[j] + tmp$age2*unique.group$AGE2[j]+tmp$age3*unique.group$AGE3[j] + tmp$age4*unique.group$AGE4[j]+
	          tmp$age5*unique.group$AGE5[j] + tmp$age6*unique.group$AGE6[j]+ tmp$inmigrant*unique.group$INMIGRANT[j]+
	          tmp$trading*unique.group$COMM_TYPE_T[j] + tmp$fishing*unique.group$COMM_TYPE_F[j]
	        tmpp<-exp(tmp2)/(1+exp(tmp2))
	        # draw samples from posterior samples of p_suc
	        SAM_P[,j]<-sample(tmpp,nprior,replace=TRUE)
	        # calculate the log empirical density estimate
	        d.sam<-approxfun(density(tmpp))
	        SAM_P_LOGD[,j]<-log(d.sam(SAM_P[,j]))
	      }
	    }else if(sample.method=='betaapprox'){
	      for (j in 1:nrow(unique.group)){
	        # posterior samples of p_suc
	        # tmp2<-tmp$a + tmp$comm[,unique.group$COMM_NUM_B[j]] + tmp$trading*unique.group$COMM_TYPE_T[j] + tmp$fishing*unique.group$COMM_TYPE_F[j] + tmp$male*unique.group$MALE[j] +
	        #   tmp$age1*unique.group$AGE1[j] + tmp$age2*unique.group$AGE2[j] + tmp$age3*unique.group$AGE3[j] +  tmp$age4*unique.group$AGE4[j] +
	        #   tmp$age5*unique.group$AGE5[j]  + tmp$age6*unique.group$AGE6[j] + tmp$inmigrant*unique.group$INMIGRANT[j]
	        tmp2<-tmp$a + tmp$comm[,unique.group$COMM_NUM_B[j]]  + tmp$male*unique.group$MALE[j] +
	          tmp$age1*unique.group$AGE1[j] + tmp$age2*unique.group$AGE2[j]+tmp$age3*unique.group$AGE3[j] + tmp$age4*unique.group$AGE4[j]+
	          tmp$age5*unique.group$AGE5[j] + tmp$age6*unique.group$AGE6[j]+ tmp$inmigrant*unique.group$INMIGRANT[j]+
	          tmp$trading*unique.group$COMM_TYPE_T[j] + tmp$fishing*unique.group$COMM_TYPE_F[j]
	        tmpp<-exp(tmp2)/(1+exp(tmp2))
	        # fit a beta distribution for p_suc
	        betapar<-fitdist(as.vector(tmpp),"beta")
	        # draw samples from the fitted beta distribution
	        SAM_P[,j]<-rbeta(nprior,betapar$estimate[1],betapar$estimate[2])
	        # calculate the log density from the fitted beta distribution
	        SAM_P_LOGD[,j]<-dbeta(S.DTL.prior$SEQ_P[,j],shape1 = betapar$estimate[1], shape2 = betapar$estimate[2],log=TRUE)
	      }
	    }
	    S.DTL.prior$SAM_P<-S.DTL.prior$SAM_P*SAM_P
	    S.DTL.prior$SAM_P_LOGD<-S.DTL.prior$SAM_P+SAM_P_LOGD
	    S.DTL.prior$EMP_S<-S.DTL.prior$EMP_S*EMP_S
	    fit.list[[i]]<-fit
	  }
	  S.DTL.prior$method<-'GLM'
	  TR.OBS<-data.table(OBS=dc$OBS,TR_CATEGORY_ID=unique.pairs$TR_CATEGORY_ID,
	                     REC_CATEGORY_ID=unique.pairs$REC_CATEGORY_ID)
	  save(S.DTL.prior,fit.list,TR.OBS,file=paste0('InputsMCMC_GLM',iteration,'.rda'))
	}
}
