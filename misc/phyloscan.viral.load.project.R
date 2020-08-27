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

vl.prevalence.by.gender.loc.age.gp<- function()
{
	require(Hmisc)
	require(data.table)
	require(ggplot2)
	require(rstan)
	vl.detectable <- 4e2
	vl.suppressed <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	outdir <- file.path(prjdir,'results_200220')
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
	
	subset(ds, FC=='inland' & SEX=='F' & AGEYRS==25)
	
	tmp <- seq.int(min(ds$AGEYRS), max(ds$AGEYRS))
	vla <- as.data.table(expand.grid(FC=c('fishing','inland'), SEX=c('M','F'), AGEYRS=tmp))
	vla <- vla[, {		
				z <- which(ds$FC==FC & ds$SEX==SEX & ds$AGEYRS==AGEYRS)	
				list(	N= length(z),
						HIV_N= length(which(ds$HIV_STATUS[z]==1)),
						VLNS_N= length(which(ds$VLNS[z]==1)),
						ARV_N= length(which(ds$ARVMED[z]==0 & ds$HIV_STATUS[z]==1 & !is.na(ds$ARVMED[z])))
				)				
			}, by=c('FC','SEX','AGEYRS')]
	setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
	vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
	vla[, SEX:= as.integer(SEX_LABEL=='M')]
	vla[, AGE:= AGE_LABEL-14L]
	vla[, ROW_ID:= seq_len(nrow(vla))]
	
	stan.code <- "
data{	
	int<lower=1> N_predict;
  	real x_predict[N_predict];
  	int<lower=1> N_observed;
  	int<lower=1, upper=N_predict> observed_idx[N_observed];
  	int y_observed[N_observed];
	int total_observed[N_observed];
  	real<lower=0> rho_hyper_par;
  	real<lower=0> alpha_hyper_par;
}

parameters {
	real<lower=0> rho;
	real<lower=0> alpha;
	real sex0_loc0;
  	vector[N_predict] f_tilde;
}

transformed parameters {
	matrix[N_predict, N_predict] cov;
  	matrix[N_predict, N_predict] L_cov;
  	vector[N_predict] logit_p_predict;
	cov = cov_exp_quad(x_predict, alpha, rho) + diag_matrix(rep_vector(1e-10, N_predict));
	L_cov = cholesky_decompose(cov);
	logit_p_predict = sex0_loc0 + L_cov * f_tilde;
}

model {
	rho ~ normal(0, rho_hyper_par);
  	alpha ~ normal(0, alpha_hyper_par);
  	sex0_loc0 ~ normal( 0 , 10 );
  	f_tilde ~ normal(0, 1);
  	y_observed ~ binomial_logit(total_observed, logit_p_predict[observed_idx] );
}

generated quantities {
  vector[N_predict] p_predict = inv_logit(logit_p_predict);  
}
"

	stan.code2 <- "
data{	
	int<lower=1> N_predict;
  	real x_predict[N_predict];
  	int<lower=1> N_observed;
  	int<lower=1, upper=N_predict> observed_idx[N_observed];
  	int y_observed_00[N_observed];
	int y_observed_10[N_observed];
	int y_observed_01[N_observed];	
	int y_observed_11[N_observed];	
	int total_observed_00[N_observed];
	int total_observed_10[N_observed];
	int total_observed_01[N_observed];
	int total_observed_11[N_observed];
  	real<lower=0> rho_hyper_par_00;
	real<lower=0> rho_hyper_par_10;
	real<lower=0> rho_hyper_par_01;
	real<lower=0> rho_hyper_par_11;  	
  	real<lower=0> alpha_hyper_par_00;
	real<lower=0> alpha_hyper_par_10;
	real<lower=0> alpha_hyper_par_01;
	real<lower=0> alpha_hyper_par_11;
}

parameters {
	real<lower=0> rho_00;
	real<lower=0> rho_10;
	real<lower=0> rho_01;
	real<lower=0> rho_11;
	real<lower=0> alpha_00;
	real<lower=0> alpha_10;
	real<lower=0> alpha_01;
	real<lower=0> alpha_11;
	real sex0_loc0;
	real sex1_loc0;
	real sex0_loc1;
	real sex1_loc1;  	
  	vector[N_predict] f_tilde_00;
	vector[N_predict] f_tilde_10;
	vector[N_predict] f_tilde_01;
	vector[N_predict] f_tilde_11;
}

transformed parameters {
  	matrix[N_predict, N_predict] L_cov;
  	vector[N_predict] logit_p_predict_00;
	vector[N_predict] logit_p_predict_10;
	vector[N_predict] logit_p_predict_01;
	vector[N_predict] logit_p_predict_11;
	// GP for 00 and 01 (women)
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_00, rho_00) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_00 = sex0_loc0 + L_cov * f_tilde_00;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_01, rho_01) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_01 = sex0_loc1 + L_cov * f_tilde_01;
	// GP for 10 and 10 (men)
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_10, rho_10) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_10 = sex1_loc0 + L_cov * f_tilde_10;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_11, rho_11) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_11 = sex1_loc1 + L_cov * f_tilde_11;
}

model {
	rho_00 ~ normal(0, rho_hyper_par_00);
	rho_10 ~ normal(0, rho_hyper_par_10);
	rho_01 ~ normal(0, rho_hyper_par_01);
	rho_11 ~ normal(0, rho_hyper_par_11);  	
  	alpha_00 ~ normal(0, alpha_hyper_par_00);
	alpha_10 ~ normal(0, alpha_hyper_par_10);
	alpha_01 ~ normal(0, alpha_hyper_par_01);
	alpha_11 ~ normal(0, alpha_hyper_par_11);
  	sex0_loc0 ~ normal( 0 , 10 );
	sex0_loc1 ~ normal( 0 , 10 );
	sex1_loc0 ~ normal( 0 , 10 );
	sex1_loc1 ~ normal( 0 , 10 );
  	f_tilde_00 ~ normal(0, 1);
	f_tilde_01 ~ normal(0, 1);
	f_tilde_10 ~ normal(0, 1);
	f_tilde_11 ~ normal(0, 1);
  	y_observed_00 ~ binomial_logit(total_observed_00, logit_p_predict_00[observed_idx] );
	y_observed_01 ~ binomial_logit(total_observed_01, logit_p_predict_01[observed_idx] );
	y_observed_10 ~ binomial_logit(total_observed_10, logit_p_predict_10[observed_idx] );
	y_observed_11 ~ binomial_logit(total_observed_11, logit_p_predict_11[observed_idx] );
}

generated quantities {
  	vector[N_predict] p_predict_00;
	vector[N_predict] p_predict_01;
	vector[N_predict] p_predict_10;
	vector[N_predict] p_predict_11;
	p_predict_00 = inv_logit(logit_p_predict_00);
	p_predict_01 = inv_logit(logit_p_predict_01);  
	p_predict_10 = inv_logit(logit_p_predict_10);
	p_predict_11 = inv_logit(logit_p_predict_11);
}
"

	#stan.model <- stan_model(model_name= 'gp_one',model_code = gsub('\t',' ',stan.code))
	stan.model <- stan_model(model_name= 'gp_all',model_code = gsub('\t',' ',stan.code2))
	vla2 <- subset(vla, SEX==0 & LOC==0)
	stan.data <- list()
	stan.data$x_predict <- seq(vla[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
	stan.data$N_predict <- length(stan.data$x_predict)
	stan.data$observed_idx <- which(stan.data$x_predict%%1==0.5)
	stan.data$N_observed <- length(stan.data$observed_idx)
	stan.data$y_observed_00 <- vla[SEX==0 & LOC==0, HIV_N]
	stan.data$y_observed_10 <- vla[SEX==1 & LOC==0, HIV_N]
	stan.data$y_observed_01 <- vla[SEX==0 & LOC==1, HIV_N]
	stan.data$y_observed_11 <- vla[SEX==1 & LOC==1, HIV_N]
	stan.data$total_observed_00 <- vla[SEX==0 & LOC==0, N]
	stan.data$total_observed_10 <- vla[SEX==1 & LOC==0, N]
	stan.data$total_observed_01 <- vla[SEX==0 & LOC==1, N]
	stan.data$total_observed_11 <- vla[SEX==1 & LOC==1, N]
	stan.data$rho_hyper_par_00 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_10 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_01 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_11 <- diff(range(stan.data$x_predict))/3
	stan.data$alpha_hyper_par_00 <- 2
	stan.data$alpha_hyper_par_10 <- 2
	stan.data$alpha_hyper_par_01 <- 2
	stan.data$alpha_hyper_par_11 <- 2
	fit <- sampling(stan.model, data=stan.data, iter=2e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
	save(fit, file=file.path(outdir, "hivprevalence_gp_stanfit_200428f.rda"))
	
	min( summary(fit)$summary[, 'n_eff'] )	
	re <- rstan::extract(fit)
	ps <- c(0.025,0.25,0.5,0.75,0.975)
	
	#	extract hyperparams rho		
	tmp <- cbind( quantile(re$rho_00, probs=ps),
			quantile(re$rho_10, probs=ps),
			quantile(re$rho_01, probs=ps),
			quantile(re$rho_11, probs=ps),
			quantile(re$alpha_00, probs=ps),
			quantile(re$alpha_10, probs=ps),
			quantile(re$alpha_01, probs=ps),
			quantile(re$alpha_11, probs=ps) )			
	colnames(tmp) <- c('rho_00','rho_10','rho_01','rho_11','alpha_00','alpha_10','alpha_01','alpha_11')
	rownames(tmp) <- c('CL','IL','M','IU','CU')
	tmp <- as.data.table(reshape2::melt(tmp))
	setnames(tmp, 'Var2', 'GP_hyper_par')
	tmp[, SEX:= as.integer(gsub('^([a-z]+)_([0-9])([0-9])','\\2',GP_hyper_par))]
	tmp[, LOC:= as.integer(gsub('^([a-z]+)_([0-9])([0-9])','\\3',GP_hyper_par))]
	tmp[, GP_hyper_par:= gsub('^([a-z]+)_([0-9])([0-9])','\\1',GP_hyper_par)]
	tmp <- dcast.data.table(tmp, LOC+SEX+GP_hyper_par~Var1, value.var='value')
	prev.hiv.gp.pars <- merge(unique(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), tmp, by=c('SEX','LOC'))
	ggplot(prev.hiv.gp.pars, aes(x= paste0(GP_hyper_par, ' ', LOC_LABEL, ' ', SEX_LABEL))) +
			geom_point(aes(y=M)) +
			geom_errorbar(aes(ymin=CL, ymax=CU)) +
			coord_flip() +
			theme_bw() +
			labs(x='GP hyperparameter\n', y='')
	ggsave(file=file.path(prjdir,'results_200220','200428f_hivprevalence_gppars.pdf'), w=6, h=3)
	
	
	#	make prevalence plot by age
	tmp <- cbind( apply(re$p_predict_00, 2, quantile, probs=ps),
			apply(re$p_predict_10, 2, quantile, probs=ps),
			apply(re$p_predict_01, 2, quantile, probs=ps),
			apply(re$p_predict_11, 2, quantile, probs=ps)
			)
	rownames(tmp) <- c('CL','IL','M','IU','CU')
	tmp <- as.data.table(reshape2::melt(tmp))
	prev.hiv.by.age <- dcast.data.table(tmp, Var2~Var1, value.var='value')
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	prev.hiv.by.age <- cbind(tmp, prev.hiv.by.age) 
	prev.hiv.by.age <- merge(unique(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), prev.hiv.by.age, by=c('SEX','LOC'))
	ggplot(prev.hiv.by.age) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
			geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			facet_wrap(~LOC_LABEL, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='HIV prevalence (95% credibility interval)\n', 
					colour='gender', 
					linetype='location')
	ggsave(file=file.path(prjdir,'results_200220','200428f_hivprevalence_vs_age_by_gender_fishinland_stan.pdf'), w=6, h=5)
	
	
	#	extract basic prevalence estimates
	ps <- c(0.025, 0.5, 0.975)
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt(tmp))
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict-0.5, SEX=c(0,1), LOC=c(0,1)))
	tmp[, Var2:= seq_len(nrow(tmp))]
	rp <- merge(tmp, rp, by='Var2')
	setnames(rp, c('Var1','value'), c('iterations','P'))
	rp <- merge(rp, vla, by=c('LOC','SEX','AGE_LABEL'))	
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
	rp <- rp[, list(Q= quantile(P, probs=ps), P=c('CL','M','CU')), by=c('LOC','SEX')]
	rp <- dcast.data.table(rp, LOC+SEX~P, value.var='Q')
	rp[, LABEL:= paste0(round(M*100, d=2),'% (',round(CL*100, d=2),'% - ',round(CU*100,d=2),'%)') ]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
	prev.hiv.by.sex.loc <- copy(rp)
	#	   LOC SEX LOC_LABEL SEX_LABEL         CL         CU         M                 LABEL
	#1:   0   0    inland         F 0.15553741 0.17130775 0.1634250 0.16% (0.16% - 0.17%)
	#2:   0   1    inland         M 0.08610792 0.09970375 0.0927843  0.09% (0.09% - 0.1%)
	#3:   1   0   fishing         F 0.42138562 0.46424074 0.4427624 0.44% (0.42% - 0.46%)
	#4:   1   1   fishing         M 0.30626567 0.34474245 0.3254705 0.33% (0.31% - 0.34%)
	
	
	#	extract prevalence ratio female:male and male:female
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt(tmp))
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict-0.5, SEX=c(0,1), LOC=c(0,1)))
	tmp[, Var2:= seq_len(nrow(tmp))]
	rp <- merge(tmp, rp, by='Var2')
	setnames(rp, c('Var1','value'), c('iterations','P'))
	rp <- merge(rp, vla, by=c('LOC','SEX','AGE_LABEL'))	
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
	rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
	rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC','iterations')]
	rp <- melt(rp, id.vars=c('LOC','iterations'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC','variable')]
	rp <- dcast.data.table(rp, LOC+variable~P, value.var='Q')
	rp[, LABEL:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL))), rp, by=c('LOC'))
	prevratio.hiv.by.loc <- copy(rp)
	#	   LOC LOC_LABEL variable        CL        CU         M              LABEL
	#1:   0    inland    PR_FM 1.6122846 1.9240272 1.7615063 1.76 (1.61 - 1.92)
	#2:   0    inland    PR_MF 0.5197432 0.6202379 0.5676960 0.57 (0.52 - 0.62)
	#3:   1   fishing    PR_FM 1.2608839 1.4687378 1.3605713 1.36 (1.26 - 1.47)
	#4:   1   fishing    PR_MF 0.6808567 0.7930944 0.7349854 0.73 (0.68 - 0.79)
	
	
	#	plot prevalence ratio female:male and male:female by age
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt(tmp))
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, Var2:= seq_len(nrow(tmp))]
	rp <- merge(tmp, rp, by='Var2')
	setnames(rp, c('Var1','value'), c('iterations','P'))
	rp <- merge(rp, unique(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), by=c('LOC','SEX'))	
	rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations+AGE_LABEL~SEX_LABEL, value.var='P')
	rp[, PR_FM:= F/M]
	rp[, PR_MF:=M/F]
	rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]
	rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
	prevratio.hiv.by.loc.age <- copy(rp)
	ggplot(subset(prevratio.hiv.by.loc.age, variable=='PR_FM')) + 	
			geom_hline(yintercept=1, linetype=2) +
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=LOC_LABEL), alpha=0.2) +
			geom_line(aes(x=AGE_LABEL, y=M)) +
			scale_x_continuous( expand=c(0,0) ) +
			scale_y_log10() +
			facet_wrap(~LOC_LABEL, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='female to male HIV prevalence ratio\n(95% credibility interval)\n')
	ggsave(file=file.path(prjdir,'results_200220','200428f_hivprevalenceratio_vs_age_by_fishinland_stan.pdf'), w=8, h=5)
	
	#	extract if difference in female:male prevalence risk ratio in fishing vs inland
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt(tmp))
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict-0.5, SEX=c(0,1), LOC=c(0,1)))
	tmp[, Var2:= seq_len(nrow(tmp))]
	rp <- merge(tmp, rp, by='Var2')
	setnames(rp, c('Var1','value'), c('iterations','P'))
	rp <- merge(vla, rp, by=c('LOC','SEX','AGE_LABEL'))
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]	
	rp <- dcast.data.table(rp, LOC_LABEL+iterations~SEX_LABEL, value.var='P')
	rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
	rp <- dcast.data.table(rp, iterations~LOC_LABEL, value.var='PR_FM')
	rp[, PR_FM_D:= inland-fishing]	
	rp[, list(Q= quantile(PR_FM_D, probs=ps), P=c('CL','M','CU'))]
	#1: 0.2182593 CL
	#2: 0.4008081  M
	#3: 0.5894657 CU
	
	save(vla, re, prev.hiv.by.age, prevratio.hiv.by.loc, prev.hiv.by.sex.loc, prevratio.hiv.by.loc.age, file=file.path(outdir, "200428f_hivprevalence.rda"))
	
	
	#	make table version suppressed
	prev.hiv.by.age[, LABEL:= paste0(sprintf('%2.1f',M*100),' (',sprintf('%2.1f',CL*100),'-',sprintf('%2.1f',CU*100),')') ]
	set(prev.hiv.by.age, NULL, 'SEX_LABEL', prev.hiv.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])
	prevratio.hiv.by.loc.age[, LABEL2:= paste0(sprintf('%2.2f',M),' (',sprintf('%2.2f',CL),'-',sprintf('%2.2f',CU),')') ]
	dt <- subset(prev.hiv.by.age, AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5))	
	dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
	tmp <- subset(prevratio.hiv.by.loc.age, variable=="PR_FM" & AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5), c(LOC_LABEL, AGE_LABEL, LABEL2))
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	tmp <- subset(prevratio.hiv.by.loc.age, variable=="PR_MF" & AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5), c(LOC_LABEL, AGE_LABEL, LABEL2))
	setnames(tmp, 'LABEL2', 'PR_MF')
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	write.csv(dt, row.names=FALSE, file=file.path(outdir, "200428f_hivprevalence.csv"))
	
}

vl.prevalence.by.gender.loc.age.icar<- function()
{
	require(Hmisc)
	require(data.table)
	require(ggplot2)
	require(rstan)
	vl.detectable <- 4e2
	vl.suppressed <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	outdir <- file.path(prjdir,'results_200220')
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
	
	subset(ds, FC=='inland' & SEX=='F' & AGEYRS==25)
	
	tmp <- seq.int(min(ds$AGEYRS), max(ds$AGEYRS))
	vla <- as.data.table(expand.grid(FC=c('fishing','inland'), SEX=c('M','F'), AGEYRS=tmp))
	vla <- vla[, {		
				z <- which(ds$FC==FC & ds$SEX==SEX & ds$AGEYRS==AGEYRS)	
				list(	N= length(z),
						HIV_N= length(which(ds$HIV_STATUS[z]==1)),
						VLNS_N= length(which(ds$VLNS[z]==1)),
						ARV_N= length(which(ds$ARVMED[z]==0 & ds$HIV_STATUS[z]==1 & !is.na(ds$ARVMED[z])))
						)				
			}, by=c('FC','SEX','AGEYRS')]
	setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
	vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
	vla[, SEX:= as.integer(SEX_LABEL=='M')]
	vla[, AGE:= AGE_LABEL-14L]
	vla[, ROW_ID:= seq_len(nrow(vla))]
	
	stan.code1 <- "functions{			
	real icar_normal_lpdf( vector phi, int N, int[] node1, int[] node2)
	{			
		return -0.5 * dot_self(phi[node1] - phi[node2]) + normal_lpdf( sum(phi) | 0, 0.001 * N);
	}
}

data{
	int<lower=1> N; 
	int<lower=1> TOTAL[N];    
	int<lower=0> K[N];
	int<lower=0> AGE_N;  
	int<lower=0, upper=AGE_N> AGE[N];  
	int<lower=0,upper=1> SEX[N];
	int<lower=0,upper=1> LOC[N];	 
	int<lower=0> N_edges; 						// number of related age groups
	int<lower=1, upper=AGE_N> node1[N_edges];  	// node1[i] adjacent to node2[i]
	int<lower=1, upper=AGE_N> node2[N_edges];  	// and node1[i] < node2[i]
}
			
parameters{
	real baseline;
	real sex;
	real loc;
	vector[AGE_N] phi; 		
	real<lower=0> sigma;
}
			
transformed parameters{
	vector[N] p_logit;
	for ( i in 1:N ) 
	{
		p_logit[i] = baseline + sex * SEX[i] + loc * LOC[i] + sigma * phi[AGE[i]];
	}
}
			
model{	
	baseline ~ normal( 0 , 10 );
	loc ~ normal( 0 , 3 );
	sex ~ normal( 0 , 3 );	
	phi ~ icar_normal_lpdf(AGE_N , node1 , node2);
	sigma ~ cauchy(0, 1);
	K ~ binomial_logit( TOTAL , p_logit );
}
			
generated quantities{
	vector[N] p;
	p= inv_logit(p_logit);
}"
	stan.code2 <- "functions{			
	real icar_normal_lpdf( vector phi, int N, int[] node1, int[] node2)
	{			
		return -0.5 * dot_self(phi[node1] - phi[node2]) + normal_lpdf( sum(phi) | 0, 0.001 * N);
	}
}
		
data{
	int<lower=1> N; 
	int<lower=1> TOTAL[N];    
	int<lower=0> K[N];
	int<lower=0> AGE_N;  
	int<lower=0, upper=AGE_N> AGE[N];  
	int<lower=0,upper=1> SEX[N];
	int<lower=0,upper=1> LOC[N];	 
	int<lower=0> N_edges; 						// number of related age groups
	int<lower=1, upper=AGE_N> node1[N_edges];  	// node1[i] adjacent to node2[i]
	int<lower=1, upper=AGE_N> node2[N_edges];  	// and node1[i] < node2[i]
}
		
parameters{
	real sex0_loc0;
	real sex1_loc0;
	real sex0_loc1;
	real sex1_loc1;
	vector[AGE_N] phi_sex0_loc0;
	vector[AGE_N] phi_sex1_loc0;
	vector[AGE_N] phi_sex0_loc1;
	vector[AGE_N] phi_sex1_loc1; 			 		
	real<lower=0> sigma_loc0;
	real<lower=0> sigma_loc1;
}
		
transformed parameters{
	vector[N] p_logit;
	for ( i in 1:N ) 
	{
		p_logit[i] = sex0_loc0 * (1-SEX[i]) * (1-LOC[i]) + 
					 sex1_loc0 * SEX[i] * (1-LOC[i]) + 
					 sex0_loc1 * (1-SEX[i]) * LOC[i] + 
					 sex1_loc1 * SEX[i] * LOC[i] +
						sigma_loc0 * phi_sex0_loc0[AGE[i]] * (1-SEX[i]) * (1-LOC[i]) +
						sigma_loc0 * phi_sex1_loc0[AGE[i]] * SEX[i] * (1-LOC[i]) +
						sigma_loc1 * phi_sex0_loc1[AGE[i]] * (1-SEX[i]) * LOC[i] +
						sigma_loc1 * phi_sex1_loc1[AGE[i]] * SEX[i] * LOC[i];
	}
}
		
model{	
	sex0_loc0 ~ normal( 0 , 10 );
	sex1_loc0 ~ normal( 0 , 10 );
	sex0_loc1 ~ normal( 0 , 10 );
	sex1_loc1 ~ normal( 0 , 10 );
	phi_sex0_loc0 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	phi_sex1_loc0 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	phi_sex0_loc1 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	phi_sex1_loc1 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	sigma_loc0 ~ normal(0.0, 1);
	sigma_loc1 ~ normal(0.0, 1);
	//sigma_loc0 ~ cauchy(0, 1);
	//sigma_loc1 ~ cauchy(0, 1);
	K ~ binomial_logit( TOTAL , p_logit );
}
		
generated quantities{
	vector[N] p;
	p= inv_logit(p_logit);
}"
	stan.model <- stan_model(model_name= 'icar_age',model_code = gsub('\t',' ',stan.code))
	stan.model2 <- stan_model(model_name= 'icar_age_interactions',model_code = gsub('\t',' ',stan.code2))
	stan.data <- list()
	stan.data$N <- nrow(vla)
	stan.data$TOTAL <- vla[,N]
	stan.data$K <- vla[,HIV_N]
	stan.data$AGE_N <- vla[, max(AGE)]
	stan.data$AGE <- vla[, AGE]
	stan.data$SEX <- vla[, SEX]
	stan.data$LOC <- vla[, LOC]
	#	second order RW prior
	stan.data$node1 <-  c(vla[, seq.int(1, max(AGE)-1L)], vla[, seq.int(1, max(AGE)-2L)])
	stan.data$node2 <-  c(vla[, seq.int(2, max(AGE))], vla[, seq.int(3, max(AGE))])
	tmp <- sort(stan.data$node1, index.return=TRUE)$ix
	stan.data$node1 <- stan.data$node1[tmp]
	stan.data$node2 <- stan.data$node2[tmp]
	stan.data$N_edges <-  length(stan.data$node1)
	fit <- sampling(stan.model2, data=stan.data, iter=20e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999) )	
	#save(fit, file=file.path(outdir, "hivprevalence_icar_stanfit_200428.rda"))		# trends by age quite rough, using Cauchy prior on sigma
	#save(fit, file=file.path(outdir, "hivprevalence_icar_stanfit_200428c.rda"))	# trends by age still quite rough, using N(0,0.1) prior on sigma
	save(fit, file=file.path(outdir, "hivprevalence_icar_stanfit_200428d.rda"))
	min( summary(fit)$summary[, 'n_eff'] )

	re <- rstan::extract(fit)
	ps <- c(0.025,0.5,0.975)
	
	quantile(re$sigma_loc0, probs=ps)
	#	in 200428b: 0.2230766 0.2943186 0.3861041 
	#	in 200428c: 0.2022420 0.2578030 0.3270586  
	#	in 200428d: 0.4380249 0.5558544 0.7091595
	
	#	make prevalence plot by age
	tmp <- apply(re$p, 2, quantile, probs=ps)
	rownames(tmp) <- c('CL','M','CU')
	tmp <- as.data.table(reshape2::melt(tmp))	
	prev.hiv.by.age <- cbind(vla, dcast.data.table(tmp, Var2~Var1, value.var='value'))	
	ggplot(prev.hiv.by.age) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
			geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			facet_wrap(~LOC_LABEL, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='HIV prevalence (95% credibility interval)\n', 
					colour='gender', 
					linetype='location')
	ggsave(file=file.path(prjdir,'results_200220','200428d_hivprevalence_vs_age_by_gender_fishinland_stan.pdf'), w=6, h=5)
	

	#	extract basic prevalence estimates
	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
	rp <- rp[, list(Q= quantile(P, probs=ps), P=c('CL','M','CU')), by=c('LOC','SEX')]
	rp <- dcast.data.table(rp, LOC+SEX~P, value.var='Q')
	rp[, LABEL:= paste0(round(M, d=2),'% (',round(CL, d=2),'% - ',round(CU,d=2),'%)') ]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
	prev.hiv.by.sex.loc <- copy(rp)
	#	   LOC SEX LOC_LABEL SEX_LABEL         CL         CU         M                 LABEL
	#1:   0   0    inland         F 0.15553741 0.17130775 0.1634250 0.16% (0.16% - 0.17%)
	#2:   0   1    inland         M 0.08610792 0.09970375 0.0927843  0.09% (0.09% - 0.1%)
	#3:   1   0   fishing         F 0.42138562 0.46424074 0.4427624 0.44% (0.42% - 0.46%)
	#4:   1   1   fishing         M 0.30626567 0.34474245 0.3254705 0.33% (0.31% - 0.34%)


	#	extract prevalence ratio female:male and male:female
	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
	rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
	rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC','iterations')]
	rp <- melt(rp, id.vars=c('LOC','iterations'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC','variable')]
	rp <- dcast.data.table(rp, LOC+variable~P, value.var='Q')
	rp[, LABEL:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL))), rp, by=c('LOC'))
	prevratio.hiv.by.loc <- copy(rp)
	#	   LOC LOC_LABEL variable        CL        CU         M              LABEL
	#1:   0    inland    PR_FM 1.6122846 1.9240272 1.7615063 1.76 (1.61 - 1.92)
	#2:   0    inland    PR_MF 0.5197432 0.6202379 0.5676960 0.57 (0.52 - 0.62)
	#3:   1   fishing    PR_FM 1.2608839 1.4687378 1.3605713 1.36 (1.26 - 1.47)
	#4:   1   fishing    PR_MF 0.6808567 0.7930944 0.7349854 0.73 (0.68 - 0.79)
	

	#	plot prevalence ratio female:male and male:female by age
	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations+AGE+AGE_LABEL~SEX_LABEL, value.var='P')
	rp[, PR_FM:= F/M]
	rp[, PR_MF:=M/F]
	rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]
	rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
	prevratio.hiv.by.loc.age <- copy(rp)
	ggplot(subset(prevratio.hiv.by.loc.age, variable=='PR_FM')) + 	
			geom_hline(yintercept=1, linetype=2) +
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=LOC_LABEL), alpha=0.2) +
			geom_line(aes(x=AGE_LABEL, y=M)) +
			scale_x_continuous( expand=c(0,0) ) + 			
			facet_wrap(~LOC_LABEL, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='female to male HIV prevalence ratio\n(95% credibility interval)\n')
	ggsave(file=file.path(prjdir,'results_200220','200428d_hivprevalenceratio_vs_age_by_fishinland_stan.pdf'), w=10, h=5)
	
	#	extract if difference in female:male prevalence risk ratio in fishing vs inland
	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
	rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
	rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
	rp <- dcast.data.table(rp, iterations~LOC_LABEL, value.var='PR_FM')
	rp[, PR_FM_D:= inland-fishing]	
	rp[, list(Q= quantile(PR_FM_D, probs=ps), P=c('CL','M','CU'))]
	#1: 0.2182593 CL
	#2: 0.4008081  M
	#3: 0.5894657 CU
	
	save(vla, re, prev.hiv.by.age, prevratio.hiv.by.loc, prev.hiv.by.sex.loc, prevratio.hiv.by.loc.age, file=file.path(outdir, "hivprevalence_200428.rda"))		
}

vl.meanviralload.by.gender.loc.age.icar<- function()
{
	require(Hmisc)
	require(data.table)
	require(ggplot2)
	require(rstan)
	vl.detectable <- 4e2
	vl.suppressed <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	outdir <- file.path(prjdir,'results_200220')
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
	vla <- vla[, {		
				z <- which(ds$FC==FC & ds$SEX==SEX & ds$AGEYRS==AGEYRS)	
				list(	N= length(z),
						VL_MEAN= mean(ds$VLC[z]),
						VL_SD= sd(ds$VLC[z]),
						VL_MEAN_SD= sd(ds$VLC[z]) / sqrt(length(z)),
						HIV_N= length(which(ds$HIV_STATUS[z]==1)),
						VLNS_N= length(which(ds$VLNS[z]==1)),
						ARV_N= length(which(ds$ARVMED[z]==0 & ds$HIV_STATUS[z]==1 & !is.na(ds$ARVMED[z])))
				)				
			}, by=c('FC','SEX','AGEYRS')]
	setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
	vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
	vla[, SEX:= as.integer(SEX_LABEL=='M')]
	vla[, AGE:= AGE_LABEL-14L]
	vla[, ROW_ID:= seq_len(nrow(vla))]
	
	ggplot(vla) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=VL_MEAN-2*VL_MEAN_SD, ymax=VL_MEAN+2*VL_MEAN_SD, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
			geom_line(aes(x=AGE_LABEL, y=VL_MEAN, colour=SEX_LABEL)) +
			scale_x_continuous( expand=c(0,0) ) + 
			#scale_y_log10() +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			facet_wrap(SEX_LABEL~LOC_LABEL, ncol=2, scales='free') +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='mean viral load\n(95% credibility interval)\n', 
					colour='gender', 
					linetype='location')
	ggsave(file=file.path(prjdir,'results_200220','200428d_mvl_vs_age_by_gender_fishinland_raw.pdf'), w=6, h=5)
	
	stan.code <- "functions{			
	real icar_normal_lpdf( vector phi, int N, int[] node1, int[] node2)
	{			
		return -0.5 * dot_self(phi[node1] - phi[node2]) + normal_lpdf( sum(phi) | 0, 0.001 * N);
	}
}
			
data{
	int<lower=1> N; 
	real MEAN[N];    
	real<lower=0> SD[N];
	int<lower=0> AGE_N;  
	int<lower=0, upper=AGE_N> AGE[N];  
	int<lower=0,upper=1> SEX[N];
	int<lower=0,upper=1> LOC[N];	 
	int<lower=0> N_edges; 						// number of related age groups
	int<lower=1, upper=AGE_N> node1[N_edges];  	// node1[i] adjacent to node2[i]
	int<lower=1, upper=AGE_N> node2[N_edges];  	// and node1[i] < node2[i]
}
			
parameters{
	real sex0_loc0;
	real sex1_loc0;
	real sex0_loc1;
	real sex1_loc1;
	vector[AGE_N] phi_sex0_loc0;
	vector[AGE_N] phi_sex1_loc0;
	vector[AGE_N] phi_sex0_loc1;
	vector[AGE_N] phi_sex1_loc1; 			 		
	real<lower=0> sigma_loc0;
	real<lower=0> sigma_loc1;
}
			
transformed parameters{
	vector[N] mu;
	for ( i in 1:N ) 
	{
		mu[i] = sex0_loc0 * (1-SEX[i]) * (1-LOC[i]) + 
			sex1_loc0 * SEX[i] * (1-LOC[i]) + 
			sex0_loc1 * (1-SEX[i]) * LOC[i] + 
			sex1_loc1 * SEX[i] * LOC[i] +
			sigma_loc0 * phi_sex0_loc0[AGE[i]] * (1-SEX[i]) * (1-LOC[i]) +
			sigma_loc0 * phi_sex1_loc0[AGE[i]] * SEX[i] * (1-LOC[i]) +
			sigma_loc1 * phi_sex0_loc1[AGE[i]] * (1-SEX[i]) * LOC[i] +
			sigma_loc1 * phi_sex1_loc1[AGE[i]] * SEX[i] * LOC[i];
	}
}
			
model{	
	sex0_loc0 ~ normal( 0 , 10 );
	sex1_loc0 ~ normal( 0 , 10 );
	sex0_loc1 ~ normal( 0 , 10 );
	sex1_loc1 ~ normal( 0 , 10 );
	phi_sex0_loc0 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	phi_sex1_loc0 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	phi_sex0_loc1 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	phi_sex1_loc1 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	sigma_loc0 ~ normal(0.0, 1);
	sigma_loc1 ~ normal(0.0, 1);
	MEAN ~ normal( mu , SD );
}
"
	
	stan.model <- stan_model(model_name= 'icar_age_interactions',model_code = gsub('\t',' ',stan.code))
	stan.data <- list()
	stan.data$N <- nrow(vla)
	stan.data$MEAN <- vla[,VL_MEAN]
	stan.data$SD <- pmax(1, vla[,VL_MEAN_SD])
	stan.data$AGE_N <- vla[, max(AGE)]
	stan.data$AGE <- vla[, AGE]
	stan.data$SEX <- vla[, SEX]
	stan.data$LOC <- vla[, LOC]
	#	second order RW prior
	stan.data$node1 <-  c(vla[, seq.int(1, max(AGE)-1L)], vla[, seq.int(1, max(AGE)-2L)])
	stan.data$node2 <-  c(vla[, seq.int(2, max(AGE))], vla[, seq.int(3, max(AGE))])
	tmp <- sort(stan.data$node1, index.return=TRUE)$ix
	stan.data$node1 <- stan.data$node1[tmp]
	stan.data$node2 <- stan.data$node2[tmp]
	stan.data$N_edges <-  length(stan.data$node1)
	fit <- sampling(stan.model, data=stan.data, iter=20e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999) )	
	save(fit, file=file.path(outdir, "mvlinpop_icar_stanfit_200428b.rda"))
	
	min( summary(fit)$summary[, 'n_eff'] )	
	re <- rstan::extract(fit)
	ps <- c(0.025,0.5,0.975)
	
	
	#	make prevalence plot by age
	tmp <- apply(re$mu, 2, quantile, probs=ps)
	rownames(tmp) <- c('CL','M','CU')
	tmp <- as.data.table(reshape2::melt(tmp))	
	mvl.by.age <- cbind(vla, dcast.data.table(tmp, Var2~Var1, value.var='value'))	
	ggplot(mvl.by.age) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
			geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
			scale_x_continuous( expand=c(0,0) ) + 
			#scale_y_log10() +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			facet_wrap(~LOC_LABEL, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='mean viral load\n(95% credibility interval)\n', 
					colour='gender', 
					linetype='location')
	ggsave(file=file.path(prjdir,'results_200220','200428d_mvl_vs_age_by_gender_fishinland_stan.pdf'), w=6, h=5)
	
}

vl.meanviralload.by.gender.loc.age.gp<- function()
{
	require(Hmisc)
	require(data.table)
	require(ggplot2)
	require(rstan)
	vl.detectable <- 4e2
	vl.suppressed <- 1e3
		
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	outdir <- file.path(prjdir,'results_200220')
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
	vla <- vla[, {		
				z <- which(ds$FC==FC & ds$SEX==SEX & ds$AGEYRS==AGEYRS)
				z2 <- which(ds$FC==FC & ds$SEX==SEX & ds$AGEYRS==AGEYRS & ds$VLC>0)
				list(	N= length(z),
						VL_MEAN= mean(ds$VLC[z]),
						VL_MEAN_SD= sd(ds$VLC[z]) / sqrt(length(z)),
						VL_SD= sd(ds$VLC[z]),						
						VLNZ_N= length(z2),
						VLNZ_MEAN= mean(log(ds$VLC[z2])),
						VLNZ_MEAN_SD= sd(log(ds$VLC[z2])) / sqrt(length(z2)),
						VLNZ_S2= var(log(ds$VLC[z2])),
						HIV_N= length(which(ds$HIV_STATUS[z]==1)),
						VLNS_N= length(which(ds$VLNS[z]==1)),
						ARV_N= length(which(ds$ARVMED[z]==0 & ds$HIV_STATUS[z]==1 & !is.na(ds$ARVMED[z])))
				)				
			}, by=c('FC','SEX','AGEYRS')]
	setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
	vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
	vla[, SEX:= as.integer(SEX_LABEL=='M')]
	vla[, AGE:= AGE_LABEL-14L]
	vla[, ROW_ID:= seq_len(nrow(vla))]
	
	
	if(0)
	{			
		stan.code <- "
data{	
	int<lower=1> N_predict;
	real x_predict[N_predict];
	int<lower=1> N_observed;
	int<lower=1, upper=N_predict> observed_idx[N_observed];
	int zero_observed[N_observed];
	int total_observed[N_observed];
	vector<lower=0>[N_observed] mean_observed;
	vector<lower=0>[N_observed] meansd_observed;
	real<lower=0> zero_rho_hyper_par;
	real<lower=0> zero_alpha_hyper_par;
	real<lower=0> mean_rho_hyper_par;
	real<lower=0> mean_alpha_hyper_par;
}
			
parameters {
	real<lower=0> zero_rho;
	real<lower=0> zero_alpha;
	real zero_base;
	vector[N_predict] zero_f_tilde;
	real<lower=0> mean_rho;
	real<lower=0> mean_alpha;			
	real mean_base;
	vector[N_predict] mean_f_tilde;
}
			
transformed parameters {
	matrix[N_predict, N_predict] L_cov;
	vector[N_predict] logit_zero_p_predict;
	vector[N_predict] zero_p_predict;
	vector[N_predict] mean_predict;
	// GP on zeros
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, zero_alpha, zero_rho) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_zero_p_predict = zero_base + L_cov * zero_f_tilde;
	zero_p_predict = inv_logit(logit_zero_p_predict);
	// GP on means
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, mean_alpha, mean_rho) + diag_matrix(rep_vector(1e-10, N_predict)));
	mean_predict = mean_base + L_cov * mean_f_tilde;

}
			
model {
	zero_rho ~ normal(0, zero_rho_hyper_par);
	zero_alpha ~ normal(0, zero_alpha_hyper_par);
	zero_base ~ normal( 0 , 100 );
	zero_f_tilde ~ normal(0, 1);
	
	mean_rho ~ normal(0, mean_rho_hyper_par);
	mean_alpha ~ normal(0, mean_alpha_hyper_par);
	mean_base ~ normal( 0 , 100 );
	mean_f_tilde ~ normal(0, 1);
	
	zero_observed ~ binomial_logit( total_observed, logit_zero_p_predict[observed_idx] );
	mean_observed ~ normal( (1-zero_p_predict[observed_idx]) .* mean_predict[observed_idx], meansd_observed + rep_vector(1e-10, N_observed));
}
			
generated quantities {
	vector[N_predict] zero_inflated_mean_predict = (1-zero_p_predict) .* mean_predict;  
}			
"
		stan.model <- stan_model(model_name= 'gp_one_zero_inflated',model_code = gsub('\t',' ',stan.code))
		vla2 <- subset(vla, SEX==0 & LOC==0)
		stan.data <- list()	
		stan.data$x_predict <- seq(vla2[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
		stan.data$N_predict <- length(stan.data$x_predict)
		stan.data$observed_idx <- which(stan.data$x_predict%%1==0.5)
		stan.data$N_observed <- length(stan.data$observed_idx)
		stan.data$zero_observed <- vla2[,VLC_ZERO] 
		stan.data$total_observed <- vla2[,N]
		stan.data$mean_observed <- vla2[,VL_MEAN]
		stan.data$meansd_observed <- vla2[,VL_MEAN_SD]	
		stan.data$zero_rho_hyper_par <- diff(range(stan.data$x_predict))/3
		stan.data$zero_alpha_hyper_par <- 2
		stan.data$mean_rho_hyper_par <- diff(range(stan.data$x_predict))/3
		stan.data$mean_alpha_hyper_par <- 2
		fit <- sampling(stan.model, data=stan.data, iter=2e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999) )	
	}
	if(0)
	{
		stan.code2 <- "
data{	
	int<lower=1> N_predict;
	real x_predict[N_predict];
	int<lower=1> N_observed;
	int<lower=1, upper=N_predict> observed_idx[N_observed];
	vector<lower=0>[N_observed] y_observed;
	vector<lower=0>[N_observed] sd_observed;
	real<lower=0> rho_hyper_par;
	real<lower=0> alpha_hyper_par;
}
			
parameters {
	real<lower=0> rho;
	real<lower=0> alpha;
	real base;
	vector[N_predict] f_tilde;	
}
			
transformed parameters {
	matrix[N_predict, N_predict] L_cov;
	vector[N_predict] mean_predict;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha, rho) + diag_matrix(rep_vector(1e-10, N_predict)));
	mean_predict = base + L_cov * f_tilde;
}
			
model {
	rho ~ normal(0, rho_hyper_par);
	alpha ~ normal(0, alpha_hyper_par);	
	base ~ normal( 0 , 100 );
	f_tilde ~ normal(0, 1);
		
	y_observed ~ normal( mean_predict[observed_idx], sd_observed );
}				
"
		stan.model <- stan_model(model_name= 'gp_one_logscale',model_code = gsub('\t',' ',stan.code2))
		vla2 <- subset(vla, SEX==0 & LOC==0)
		stan.data <- list()	
		stan.data$x_predict <- seq(vla2[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
		stan.data$N_predict <- length(stan.data$x_predict)
		stan.data$observed_idx <- which(stan.data$x_predict%%1==0.5)
		stan.data$observed_idx <- stan.data$observed_idx[ !is.na(vla2$VLNZ_MEAN_SD) ]
		stan.data$N_observed <- length(stan.data$observed_idx)
		stan.data$y_observed <- vla2[!is.na(VLNZ_MEAN_SD),VLNZ_MEAN] 
		stan.data$sd_observed <- vla2[!is.na(VLNZ_MEAN_SD),VLNZ_MEAN_SD]
		stan.data$rho_hyper_par <- diff(range(stan.data$x_predict))/3
		stan.data$alpha_hyper_par <- 2
		fit <- sampling(stan.model, data=stan.data, iter=2e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999) )
	}
	
	
	
	
	
	
	stan.code3 <- "
data{	
	int<lower=1> N_predict;
	real x_predict[N_predict];
	int<lower=1> N_observed;
	int<lower=1, upper=N_predict> observed_idx[N_observed];
	vector<lower=0>[N_observed] y_observed;
	vector<lower=0>[N_observed] s2_observed; // with denominator (n-1)
	vector<lower=1>[N_observed] n_observed;
	real<lower=0> m_rho_hyper_par_alpha;
	real<lower=0> m_rho_hyper_par_beta;
	real<lower=0> m_alpha_hyper_par;
	real<lower=0> s_rho_hyper_par_alpha;
	real<lower=0> s_rho_hyper_par_beta;
	real<lower=0> s_alpha_hyper_par;
}

transformed data{
	vector[N_observed] n_m1_observed;
	vector[N_observed] inv_sqrt_n_observed;
	vector[N_observed] inv_scaled_s2_observed;

	n_m1_observed = n_observed - rep_vector(1, N_observed);
	inv_sqrt_n_observed = inv(sqrt(n_observed));
	inv_scaled_s2_observed = inv(s2_observed .* n_m1_observed);
}
			
parameters {
	real<lower=0> m_rho;
	real<lower=0> m_alpha;
	real<lower=0> s_rho;
	real<lower=0> s_alpha;
	real m_base;
	real s_base;
	vector[N_predict] m_f_tilde;
	vector[N_predict] s_f_tilde;	
}
			
transformed parameters {
	matrix[N_predict, N_predict] L_cov;
	vector[N_predict] m_predict;
	vector[N_predict] s_predict;
	// GP on mean of lognormal
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, m_alpha, m_rho) + diag_matrix(rep_vector(1e-10, N_predict)));
	m_predict = m_base + L_cov * m_f_tilde;
	// GP on log sigma of lognormal
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, s_alpha, s_rho) + diag_matrix(rep_vector(1e-10, N_predict)));
	s_predict = exp( s_base + L_cov * s_f_tilde);
}
			
model {
	m_rho ~ inv_gamma(m_rho_hyper_par_alpha, m_rho_hyper_par_beta);
	m_alpha ~ normal(0, m_alpha_hyper_par);	
	m_base ~ normal( 0 , 10 );
	s_rho ~ inv_gamma(s_rho_hyper_par_alpha, s_rho_hyper_par_beta);
	s_alpha ~ normal(0, s_alpha_hyper_par);	
	s_base ~ normal( 0 , 10 );	
	m_f_tilde ~ normal(0, 1);
	s_f_tilde ~ normal(0, 1);

	target+= normal_lpdf( y_observed | m_predict[observed_idx], s_predict[observed_idx] .* inv_sqrt_n_observed);
	target+= inv_chi_square_lpdf( s_predict[observed_idx] .* s_predict[observed_idx] .* inv_scaled_s2_observed | n_m1_observed );  	
}

generated quantities {
	vector[N_predict] exp_m_predict = exp( m_predict + s_predict .* s_predict .* rep_vector(0.5, N_predict) );  
}			
"

	
	stan.model <- stan_model(model_name= 'gp_one_logscale',model_code = gsub('\t',' ',stan.code3))	
	vla.select <- vla[, which(SEX==1 & LOC==1)]
	vla.select <- vla[, which(SEX==1 & LOC==0)]
	vla.select <- vla[, which(SEX==0 & LOC==1)]
	stan.data <- list()	
	stan.data$x_predict <- seq(vla[vla.select, min(AGE_LABEL)], vla[vla.select, max(AGE_LABEL)+1], 0.5)
	stan.data$N_predict <- length(stan.data$x_predict)
	stan.data$observed_idx <- which(stan.data$x_predict%%1==0.5)
	stan.data$observed_idx <- stan.data$observed_idx[ vla[vla.select,][,!is.na(VLNZ_S2)] ]
	stan.data$N_observed <- length(stan.data$observed_idx)	
	stan.data$y_observed <- vla[vla.select,][!is.na(VLNZ_S2),VLNZ_MEAN]
	stopifnot( length(stan.data$y_observed)==stan.data$N_observed )
	stan.data$n_observed <- vla[vla.select,][!is.na(VLNZ_S2),VLNZ_N]
	stopifnot( length(stan.data$n_observed)==stan.data$N_observed )
	stan.data$s2_observed <- vla[vla.select,][!is.na(VLNZ_S2),VLNZ_S2]
	stopifnot( length(stan.data$s2_observed)==stan.data$N_observed )	
	cl.target <- diff(range(stan.data$x_predict))/10
	cu.target <- diff(range(stan.data$x_predict))
	inv.gamma.err <- function(pars){ abs(sum( qinvgamma(c(0.025,0.975), exp(pars[1]), exp(pars[2])) - c(cl.target,cu.target) )) }	
	tmp <- optim(log(c(8,30)), inv.gamma.err)
	inv.gamma.pars <- exp(tmp$par)
	stan.data$m_rho_hyper_par_alpha <- inv.gamma.pars[1]
	stan.data$m_rho_hyper_par_beta <- inv.gamma.pars[2]
	stan.data$m_alpha_hyper_par <- 2
	stan.data$s_rho_hyper_par_alpha <- inv.gamma.pars[1]
	stan.data$s_rho_hyper_par_beta <- inv.gamma.pars[2]
	stan.data$s_alpha_hyper_par <- 2
	fit <- sampling(stan.model, data=stan.data, iter=20e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999) )
	#save(fit, file=file.path(outdir, "200428f_hivmeans_gp_stanfit_10.rda"))
	save(fit, file=file.path(outdir, "200428f_hivmeans_gp_stanfit_01.rda"))

	stan.code4 <- "
data{	
	int<lower=1> N_predict;
	real x_predict[N_predict];
	int<lower=1> N_observed_00;
	int<lower=1> N_observed_10;
	int<lower=1> N_observed_01;
	int<lower=1> N_observed_11;
	int<lower=1, upper=N_predict> observed_idx_00[N_observed_00];
	int<lower=1, upper=N_predict> observed_idx_10[N_observed_10];
	int<lower=1, upper=N_predict> observed_idx_01[N_observed_01];
	int<lower=1, upper=N_predict> observed_idx_11[N_observed_11];
	vector<lower=0>[N_observed_00] y_observed_00;
	vector<lower=0>[N_observed_10] y_observed_10;
	vector<lower=0>[N_observed_01] y_observed_01;
	vector<lower=0>[N_observed_11] y_observed_11;
	vector<lower=0>[N_observed_00] s2_observed_00; // with denominator (n-1)
	vector<lower=0>[N_observed_10] s2_observed_10; // with denominator (n-1)
	vector<lower=0>[N_observed_01] s2_observed_01; // with denominator (n-1)
	vector<lower=0>[N_observed_11] s2_observed_11; // with denominator (n-1)
	vector<lower=1>[N_observed_00] n_observed_00;
	vector<lower=1>[N_observed_10] n_observed_10;
	vector<lower=1>[N_observed_01] n_observed_01;
	vector<lower=1>[N_observed_11] n_observed_11;
	real<lower=0> m_rho_hyper_par_00;
	real<lower=0> m_alpha_hyper_par_00;
	real<lower=0> s_rho_hyper_par_00;
	real<lower=0> s_alpha_hyper_par_00;
	real<lower=0> m_rho_hyper_par_10;
	real<lower=0> m_alpha_hyper_par_10;
	real<lower=0> s_rho_hyper_par_10;
	real<lower=0> s_alpha_hyper_par_10;
	real<lower=0> m_rho_hyper_par_01;
	real<lower=0> m_alpha_hyper_par_01;
	real<lower=0> s_rho_hyper_par_01;
	real<lower=0> s_alpha_hyper_par_01;
	real<lower=0> m_rho_hyper_par_11;
	real<lower=0> m_alpha_hyper_par_11;
	real<lower=0> s_rho_hyper_par_11;
	real<lower=0> s_alpha_hyper_par_11;
}

transformed data{
	vector[N_observed_00] n_m1_observed_00;
	vector[N_observed_00] inv_sqrt_n_observed_00;
	vector[N_observed_00] inv_scaled_s2_observed_00;
	vector[N_observed_10] n_m1_observed_10;
	vector[N_observed_10] inv_sqrt_n_observed_10;
	vector[N_observed_10] inv_scaled_s2_observed_10;
	vector[N_observed_01] n_m1_observed_01;
	vector[N_observed_01] inv_sqrt_n_observed_01;
	vector[N_observed_01] inv_scaled_s2_observed_01;
	vector[N_observed_11] n_m1_observed_11;
	vector[N_observed_11] inv_sqrt_n_observed_11;
	vector[N_observed_11] inv_scaled_s2_observed_11;

	n_m1_observed_00 = n_observed_00 - rep_vector(1, N_observed_00);
	inv_sqrt_n_observed_00 = inv(sqrt(n_observed_00));
	inv_scaled_s2_observed_00 = inv(s2_observed_00 .* n_m1_observed_00);
	n_m1_observed_10 = n_observed_10 - rep_vector(1, N_observed_10);
	inv_sqrt_n_observed_10 = inv(sqrt(n_observed_10));
	inv_scaled_s2_observed_10 = inv(s2_observed_10 .* n_m1_observed_10);
	n_m1_observed_01 = n_observed_01 - rep_vector(1, N_observed_01);
	inv_sqrt_n_observed_01 = inv(sqrt(n_observed_01));
	inv_scaled_s2_observed_01 = inv(s2_observed_01 .* n_m1_observed_01);
	n_m1_observed_11 = n_observed_11 - rep_vector(1, N_observed_11);
	inv_sqrt_n_observed_11 = inv(sqrt(n_observed_11));
	inv_scaled_s2_observed_11 = inv(s2_observed_11 .* n_m1_observed_11);
}
			
parameters {
	real<lower=0> m_rho_00;
	real<lower=0> m_alpha_00;
	real<lower=0> s_rho_00;
	real<lower=0> s_alpha_00;
	real m_base_00;
	real s_base_00;
	vector[N_predict] m_f_tilde_00;
	vector[N_predict] s_f_tilde_00;	
	real<lower=0> m_rho_10;
	real<lower=0> m_alpha_10;
	real<lower=0> s_rho_10;
	real<lower=0> s_alpha_10;
	real m_base_10;
	real s_base_10;
	vector[N_predict] m_f_tilde_10;
	vector[N_predict] s_f_tilde_10;	
	real<lower=0> m_rho_01;
	real<lower=0> m_alpha_01;
	real<lower=0> s_rho_01;
	real<lower=0> s_alpha_01;
	real m_base_01;
	real s_base_01;
	vector[N_predict] m_f_tilde_01;
	vector[N_predict] s_f_tilde_01;	
	real<lower=0> m_rho_11;
	real<lower=0> m_alpha_11;
	real<lower=0> s_rho_11;
	real<lower=0> s_alpha_11;
	real m_base_11;
	real s_base_11;
	vector[N_predict] m_f_tilde_11;
	vector[N_predict] s_f_tilde_11;	
}
			
transformed parameters {
	matrix[N_predict, N_predict] L_cov;
	vector[N_predict] m_predict_00;
	vector[N_predict] s_predict_00;
	vector[N_predict] m_predict_10;
	vector[N_predict] s_predict_10;
	vector[N_predict] m_predict_01;
	vector[N_predict] s_predict_01;
	vector[N_predict] m_predict_11;
	vector[N_predict] s_predict_11;
	// GP on mean of lognormal
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, m_alpha_00, m_rho_00) + diag_matrix(rep_vector(1e-10, N_predict)));
	m_predict_00 = m_base_00 + L_cov * m_f_tilde_00;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, m_alpha_10, m_rho_10) + diag_matrix(rep_vector(1e-10, N_predict)));
	m_predict_10 = m_base_10 + L_cov * m_f_tilde_10;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, m_alpha_01, m_rho_01) + diag_matrix(rep_vector(1e-10, N_predict)));
	m_predict_01 = m_base_01 + L_cov * m_f_tilde_01;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, m_alpha_11, m_rho_11) + diag_matrix(rep_vector(1e-10, N_predict)));
	m_predict_11 = m_base_11 + L_cov * m_f_tilde_11;
	// GP on log sigma of lognormal
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, s_alpha_00, s_rho_00) + diag_matrix(rep_vector(1e-10, N_predict)));
	s_predict_00 = exp( s_base_00 + L_cov * s_f_tilde_00);
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, s_alpha_10, s_rho_10) + diag_matrix(rep_vector(1e-10, N_predict)));
	s_predict_10 = exp( s_base_10 + L_cov * s_f_tilde_10);
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, s_alpha_01, s_rho_01) + diag_matrix(rep_vector(1e-10, N_predict)));
	s_predict_01 = exp( s_base_01 + L_cov * s_f_tilde_01);
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, s_alpha_11, s_rho_11) + diag_matrix(rep_vector(1e-10, N_predict)));
	s_predict_11 = exp( s_base_11 + L_cov * s_f_tilde_11);
}
			
model {
	m_rho_00 ~ normal(0, m_rho_hyper_par_00);
	m_alpha_00 ~ normal(0, m_alpha_hyper_par_00);	
	m_base_00 ~ normal( 0 , 10 );
	s_rho_00 ~ normal(0, s_rho_hyper_par_00);
	s_alpha_00 ~ normal(0, s_alpha_hyper_par_00);	
	s_base_00 ~ normal( 0 , 10 );	
	m_f_tilde_00 ~ normal(0, 1);
	s_f_tilde_00 ~ normal(0, 1);

	m_rho_10 ~ normal(0, m_rho_hyper_par_10);
	m_alpha_10 ~ normal(0, m_alpha_hyper_par_10);	
	m_base_10 ~ normal( 0 , 10 );
	s_rho_10 ~ normal(0, s_rho_hyper_par_10);
	s_alpha_10 ~ normal(0, s_alpha_hyper_par_10);	
	s_base_10 ~ normal( 0 , 10 );	
	m_f_tilde_10 ~ normal(0, 1);
	s_f_tilde_10 ~ normal(0, 1);

	m_rho_01 ~ normal(0, m_rho_hyper_par_01);
	m_alpha_01 ~ normal(0, m_alpha_hyper_par_01);	
	m_base_01 ~ normal( 0 , 10 );
	s_rho_01 ~ normal(0, s_rho_hyper_par_01);
	s_alpha_01 ~ normal(0, s_alpha_hyper_par_01);	
	s_base_01 ~ normal( 0 , 10 );	
	m_f_tilde_01 ~ normal(0, 1);
	s_f_tilde_01 ~ normal(0, 1);

	m_rho_11 ~ normal(0, m_rho_hyper_par_11);
	m_alpha_11 ~ normal(0, m_alpha_hyper_par_11);	
	m_base_11 ~ normal( 0 , 10 );
	s_rho_11 ~ normal(0, s_rho_hyper_par_11);
	s_alpha_11 ~ normal(0, s_alpha_hyper_par_11);	
	s_base_11 ~ normal( 0 , 10 );	
	m_f_tilde_11 ~ normal(0, 1);
	s_f_tilde_11 ~ normal(0, 1);

	target+= normal_lpdf( y_observed_00 | m_predict_00[observed_idx_00], s_predict_00[observed_idx_00] .* inv_sqrt_n_observed_00);
	target+= inv_chi_square_lpdf( s_predict_00[observed_idx_00] .* s_predict_00[observed_idx_00] .* inv_scaled_s2_observed_00 | n_m1_observed_00 );
	target+= normal_lpdf( y_observed_10 | m_predict_10[observed_idx_10], s_predict_10[observed_idx_10] .* inv_sqrt_n_observed_10);
	target+= inv_chi_square_lpdf( s_predict_10[observed_idx_10] .* s_predict_10[observed_idx_10] .* inv_scaled_s2_observed_10 | n_m1_observed_10 );
	target+= normal_lpdf( y_observed_01 | m_predict_01[observed_idx_01], s_predict_01[observed_idx_01] .* inv_sqrt_n_observed_01);
	target+= inv_chi_square_lpdf( s_predict_01[observed_idx_01] .* s_predict_01[observed_idx_01] .* inv_scaled_s2_observed_01 | n_m1_observed_01 );
	target+= normal_lpdf( y_observed_11 | m_predict_11[observed_idx_11], s_predict_11[observed_idx_11] .* inv_sqrt_n_observed_11);
	target+= inv_chi_square_lpdf( s_predict_11[observed_idx_11] .* s_predict_11[observed_idx_11] .* inv_scaled_s2_observed_11 | n_m1_observed_11 );  	
}

generated quantities {
	vector[N_predict] exp_m_predict_00;
	vector[N_predict] exp_m_predict_10;
	vector[N_predict] exp_m_predict_01;
	vector[N_predict] exp_m_predict_11; 
	exp_m_predict_00 = exp( m_predict_00 + s_predict_00 .* s_predict_00 .* rep_vector(0.5, N_predict) );
	exp_m_predict_10 = exp( m_predict_10 + s_predict_10 .* s_predict_10 .* rep_vector(0.5, N_predict) );
	exp_m_predict_01 = exp( m_predict_01 + s_predict_01 .* s_predict_01 .* rep_vector(0.5, N_predict) );
	exp_m_predict_11 = exp( m_predict_11 + s_predict_11 .* s_predict_11 .* rep_vector(0.5, N_predict) );  
}	
"
	
	stan.model <- stan_model(model_name= 'gp_all_logscale',model_code = gsub('\t',' ',stan.code4))	
	stan.data <- list()	
	stan.data$x_predict <- seq(vla[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
	stan.data$N_predict <- length(stan.data$x_predict)
	stan.data$observed_idx_00 <- which(stan.data$x_predict%%1==0.5)
	stan.data$observed_idx_00 <- stan.data$observed_idx_00[ subset(vla,SEX==0 & LOC==0)[, !is.na(VLNZ_S2)] ]
	stan.data$N_observed_00 <- length(stan.data$observed_idx_00)	
	stan.data$observed_idx_10 <- which(stan.data$x_predict%%1==0.5)
	stan.data$observed_idx_10 <- stan.data$observed_idx_10[ subset(vla,SEX==1 & LOC==0)[, !is.na(VLNZ_S2)] ]
	stan.data$N_observed_10 <- length(stan.data$observed_idx_10)	
	stan.data$observed_idx_01 <- which(stan.data$x_predict%%1==0.5)
	stan.data$observed_idx_01 <- stan.data$observed_idx_01[ subset(vla,SEX==0 & LOC==1)[, !is.na(VLNZ_S2)] ]
	stan.data$N_observed_01 <- length(stan.data$observed_idx_01)	
	stan.data$observed_idx_11 <- which(stan.data$x_predict%%1==0.5)
	stan.data$observed_idx_11 <- stan.data$observed_idx_11[ subset(vla,SEX==1 & LOC==1)[, !is.na(VLNZ_S2)] ]
	stan.data$N_observed_11 <- length(stan.data$observed_idx_11)
	stan.data$y_observed_00 <- vla[SEX==0 & LOC==0 & !is.na(VLNZ_S2),VLNZ_MEAN]
	stan.data$y_observed_10 <- vla[SEX==1 & LOC==0 & !is.na(VLNZ_S2),VLNZ_MEAN]
	stan.data$y_observed_01 <- vla[SEX==0 & LOC==1 & !is.na(VLNZ_S2),VLNZ_MEAN]
	stan.data$y_observed_11 <- vla[SEX==1 & LOC==1 & !is.na(VLNZ_S2),VLNZ_MEAN]
	stopifnot( length(stan.data$y_observed_00)==stan.data$N_observed_00 )
	stopifnot( length(stan.data$y_observed_10)==stan.data$N_observed_10 )
	stopifnot( length(stan.data$y_observed_01)==stan.data$N_observed_01 )
	stopifnot( length(stan.data$y_observed_11)==stan.data$N_observed_11 )
	stan.data$n_observed_00 <- vla[SEX==0 & LOC==0 & !is.na(VLNZ_S2),VLNZ_N]
	stan.data$n_observed_10 <- vla[SEX==1 & LOC==0 & !is.na(VLNZ_S2),VLNZ_N]
	stan.data$n_observed_01 <- vla[SEX==0 & LOC==1 & !is.na(VLNZ_S2),VLNZ_N]
	stan.data$n_observed_11 <- vla[SEX==1 & LOC==1 & !is.na(VLNZ_S2),VLNZ_N]
	stopifnot( length(stan.data$n_observed_00)==stan.data$N_observed_00 )
	stopifnot( length(stan.data$n_observed_10)==stan.data$N_observed_10 )
	stopifnot( length(stan.data$n_observed_01)==stan.data$N_observed_01 )
	stopifnot( length(stan.data$n_observed_11)==stan.data$N_observed_11 )
	stan.data$s2_observed_00 <- vla[SEX==0 & LOC==0 & !is.na(VLNZ_S2),VLNZ_S2]
	stan.data$s2_observed_10 <- vla[SEX==1 & LOC==0 & !is.na(VLNZ_S2),VLNZ_S2]
	stan.data$s2_observed_01 <- vla[SEX==0 & LOC==1 & !is.na(VLNZ_S2),VLNZ_S2]
	stan.data$s2_observed_11 <- vla[SEX==1 & LOC==1 & !is.na(VLNZ_S2),VLNZ_S2]
	stopifnot( length(stan.data$s2_observed_00)==stan.data$N_observed_00 )
	stopifnot( length(stan.data$s2_observed_10)==stan.data$N_observed_10 )
	stopifnot( length(stan.data$s2_observed_01)==stan.data$N_observed_01 )
	stopifnot( length(stan.data$s2_observed_11)==stan.data$N_observed_11 )
	stan.data$m_rho_hyper_par_00 <- diff(range(stan.data$x_predict))/3
	stan.data$m_alpha_hyper_par_00 <- 2
	stan.data$s_rho_hyper_par_00 <- diff(range(stan.data$x_predict))/3
	stan.data$s_alpha_hyper_par_00 <- 2
	stan.data$m_rho_hyper_par_10 <- diff(range(stan.data$x_predict))/3
	stan.data$m_alpha_hyper_par_10 <- 2
	stan.data$s_rho_hyper_par_10 <- diff(range(stan.data$x_predict))/3
	stan.data$s_alpha_hyper_par_10 <- 2
	stan.data$m_rho_hyper_par_01 <- diff(range(stan.data$x_predict))/3
	stan.data$m_alpha_hyper_par_01 <- 2
	stan.data$s_rho_hyper_par_01 <- diff(range(stan.data$x_predict))/3
	stan.data$s_alpha_hyper_par_01 <- 2
	stan.data$m_rho_hyper_par_11 <- diff(range(stan.data$x_predict))/3
	stan.data$m_alpha_hyper_par_11 <- 2
	stan.data$s_rho_hyper_par_11 <- diff(range(stan.data$x_predict))/3
	stan.data$s_alpha_hyper_par_11 <- 2	
	fit <- sampling(stan.model, data=stan.data, iter=20e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999) )
	save(fit, file=file.path(outdir, "200428f_hivmeans_gp_stanfit.rda"))
	
	stan.code2b <- "
data{	
	int<lower=1> N_predict;
  	real x_predict[N_predict];
  	int<lower=1> N_observed;
  	int<lower=1, upper=N_predict> observed_idx[N_observed];
  	int y_observed_00[N_observed];
	int y_observed_10[N_observed];
	int y_observed_01[N_observed];	
	int y_observed_11[N_observed];	
	int total_observed_00[N_observed];
	int total_observed_10[N_observed];
	int total_observed_01[N_observed];
	int total_observed_11[N_observed];
  	real<lower=0> rho_hyper_par_00;
	real<lower=0> rho_hyper_par_10;
	real<lower=0> rho_hyper_par_01;
	real<lower=0> rho_hyper_par_11;  	
  	real<lower=0> alpha_hyper_par_00;
	real<lower=0> alpha_hyper_par_10;
	real<lower=0> alpha_hyper_par_01;
	real<lower=0> alpha_hyper_par_11;
}

parameters {
	real<lower=0> rho_00;
	real<lower=0> rho_10;
	real<lower=0> rho_01;
	real<lower=0> rho_11;
	real<lower=0> alpha_00;
	real<lower=0> alpha_10;
	real<lower=0> alpha_01;
	real<lower=0> alpha_11;
	real sex0_loc0;
	real sex1_loc0;
	real sex0_loc1;
	real sex1_loc1;  	
  	vector[N_predict] f_tilde_00;
	vector[N_predict] f_tilde_10;
	vector[N_predict] f_tilde_01;
	vector[N_predict] f_tilde_11;
}

transformed parameters {
  	matrix[N_predict, N_predict] L_cov;
  	vector[N_predict] logit_p_predict_00;
	vector[N_predict] logit_p_predict_10;
	vector[N_predict] logit_p_predict_01;
	vector[N_predict] logit_p_predict_11;
	// GP for 00 and 01 (women)
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_00, rho_00) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_00 = sex0_loc0 + L_cov * f_tilde_00;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_01, rho_01) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_01 = sex0_loc1 + L_cov * f_tilde_01;
	// GP for 10 and 10 (men)
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_10, rho_10) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_10 = sex1_loc0 + L_cov * f_tilde_10;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_11, rho_11) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_11 = sex1_loc1 + L_cov * f_tilde_11;
}

model {
	rho_00 ~ normal(0, rho_hyper_par_00);
	rho_10 ~ normal(0, rho_hyper_par_10);
	rho_01 ~ normal(0, rho_hyper_par_01);
	rho_11 ~ normal(0, rho_hyper_par_11);  	
  	alpha_00 ~ normal(0, alpha_hyper_par_00);
	alpha_10 ~ normal(0, alpha_hyper_par_10);
	alpha_01 ~ normal(0, alpha_hyper_par_01);
	alpha_11 ~ normal(0, alpha_hyper_par_11);
  	sex0_loc0 ~ normal( 0 , 10 );
	sex0_loc1 ~ normal( 0 , 10 );
	sex1_loc0 ~ normal( 0 , 10 );
	sex1_loc1 ~ normal( 0 , 10 );
  	f_tilde_00 ~ normal(0, 1);
	f_tilde_01 ~ normal(0, 1);
	f_tilde_10 ~ normal(0, 1);
	f_tilde_11 ~ normal(0, 1);
  	y_observed_00 ~ binomial_logit(total_observed_00, logit_p_predict_00[observed_idx] );
	y_observed_01 ~ binomial_logit(total_observed_01, logit_p_predict_01[observed_idx] );
	y_observed_10 ~ binomial_logit(total_observed_10, logit_p_predict_10[observed_idx] );
	y_observed_11 ~ binomial_logit(total_observed_11, logit_p_predict_11[observed_idx] );
}

generated quantities {
  	vector[N_predict] p_predict_00;
	vector[N_predict] p_predict_01;
	vector[N_predict] p_predict_10;
	vector[N_predict] p_predict_11;
	p_predict_00 = inv_logit(logit_p_predict_00);
	p_predict_01 = inv_logit(logit_p_predict_01);  
	p_predict_10 = inv_logit(logit_p_predict_10);
	p_predict_11 = inv_logit(logit_p_predict_11);
}
"
		
	stan.modelB <- stan_model(model_name= 'gp_all',model_code = gsub('\t',' ',stan.code2b))
	stan.dataB <- list()
	stan.dataB$x_predict <- seq(vla[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
	stan.dataB$N_predict <- length(stan.dataB$x_predict)
	stan.dataB$observed_idx <- which(stan.dataB$x_predict%%1==0.5)
	stan.dataB$N_observed <- length(stan.dataB$observed_idx)
	stan.dataB$y_observed_00 <- vla[SEX==0 & LOC==0, N-VLNZ_N]
	stan.dataB$y_observed_10 <- vla[SEX==1 & LOC==0, N-VLNZ_N]
	stan.dataB$y_observed_01 <- vla[SEX==0 & LOC==1, N-VLNZ_N]
	stan.dataB$y_observed_11 <- vla[SEX==1 & LOC==1, N-VLNZ_N]
	stan.dataB$total_observed_00 <- vla[SEX==0 & LOC==0, N]
	stan.dataB$total_observed_10 <- vla[SEX==1 & LOC==0, N]
	stan.dataB$total_observed_01 <- vla[SEX==0 & LOC==1, N]
	stan.dataB$total_observed_11 <- vla[SEX==1 & LOC==1, N]
	stan.dataB$rho_hyper_par_00 <- diff(range(stan.dataB$x_predict))/3
	stan.dataB$rho_hyper_par_10 <- diff(range(stan.dataB$x_predict))/3
	stan.dataB$rho_hyper_par_01 <- diff(range(stan.dataB$x_predict))/3
	stan.dataB$rho_hyper_par_11 <- diff(range(stan.dataB$x_predict))/3
	stan.dataB$alpha_hyper_par_00 <- 2
	stan.dataB$alpha_hyper_par_10 <- 2
	stan.dataB$alpha_hyper_par_01 <- 2
	stan.dataB$alpha_hyper_par_11 <- 2
	fitB <- sampling(stan.modelB, data=stan.dataB, iter=2e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
	save(fitB, file=file.path(outdir, "200428f_hivzeros_gp_stanfit.rda"))
	
	
	load( file.path(outdir, "200428f_hivzeros_gp_stanfit.rda") )
	load( file.path(outdir, "200428f_hivmeans_gp_stanfit_10.rda") )
	fit10 <- fit
	load( file.path(outdir, "200428f_hivmeans_gp_stanfit_11.rda") )
	fit11 <- fit
	load( file.path(outdir, "200428f_hivmeans_gp_stanfit_01.rda") )
	fit01 <- fit
	load( file.path(outdir, "200428f_hivmeans_gp_stanfit.rda") )
	
	#
	#	what we expect roughly
	ggplot(vla2, aes(x=AGE_LABEL)) +
			geom_point(aes(y= VLNZ_N/N * exp(VLNZ_MEAN+VLNZ_S2/2)))
	
	tmp <- summary(fit)$summary
	tmp[grepl('exp_m_predict',rownames(tmp)),]
	tmp[grepl('m_predict',rownames(tmp)),]
	tmp[grepl('s_predict',rownames(tmp)),]
	
	#
	#	get posterior samples
	reB <- rstan::extract(fitB)
	re01 <- rstan::extract(fit01)
	re10 <- rstan::extract(fit10)
	re11 <- rstan::extract(fit11)
	re <- rstan::extract(fit)
	re$exp_m_predict_10 <- re10$exp_m_predict
	re$exp_m_predict_01 <- re01$exp_m_predict
	re$exp_m_predict_11 <- re11$exp_m_predict
	re$m_rho_10 <- re10$m_rho
	re$m_rho_01 <- re01$m_rho	
	re$m_rho_11 <- re11$m_rho
	re$s_rho_10 <- re10$s_rho
	re$s_rho_01 <- re01$s_rho
	re$s_rho_11 <- re11$s_rho	
	
	
	#
	#	extract hyperparams rho
	ps <- c(0.025,0.25,0.5,0.75,0.975)
	tmp <- cbind( quantile(re$m_rho_00, probs=ps),
			quantile(re$m_rho_10, probs=ps),
			quantile(re$m_rho_01, probs=ps),
			quantile(re$m_rho_11, probs=ps),
			quantile(re$m_alpha_00, probs=ps),
			quantile(re$m_alpha_10, probs=ps),
			quantile(re$m_alpha_01, probs=ps),
			quantile(re$m_alpha_11, probs=ps), 	
			quantile(re$s_rho_00, probs=ps),
			quantile(re$s_rho_10, probs=ps),
			quantile(re$s_rho_01, probs=ps),
			quantile(re$s_rho_11, probs=ps),
			quantile(re$s_alpha_00, probs=ps),
			quantile(re$s_alpha_10, probs=ps),
			quantile(re$s_alpha_01, probs=ps),
			quantile(re$s_alpha_11, probs=ps)			
			)			
	colnames(tmp) <- c(	'm_rho_00','m_rho_10','m_rho_01','m_rho_11','m_alpha_00','m_alpha_10','m_alpha_01','m_alpha_11',
						's_rho_00','s_rho_10','s_rho_01','s_rho_11','s_alpha_00','s_alpha_10','s_alpha_01','s_alpha_11')
	rownames(tmp) <- c('CL','IL','M','IU','CU')
	tmp <- as.data.table(reshape2::melt(tmp))
	setnames(tmp, 'Var2', 'GP_hyper_par')
	tmp[, SEX:= as.integer(gsub('^([a-z]+_[a-z]+)_([0-9])([0-9])','\\2',GP_hyper_par))]
	tmp[, LOC:= as.integer(gsub('^([a-z]+_[a-z]+)_([0-9])([0-9])','\\3',GP_hyper_par))]
	tmp[, GP_hyper_par:= gsub('^([a-z]+_[a-z]+)_([0-9])([0-9])','\\1',GP_hyper_par)]
	tmp <- dcast.data.table(tmp, LOC+SEX+GP_hyper_par~Var1, value.var='value')
	mvl.gp.pars <- merge(unique(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), tmp, by=c('SEX','LOC'))
	ggplot(mvl.gp.pars, aes(x= paste0(GP_hyper_par, ' ', LOC_LABEL, ' ', SEX_LABEL))) +
			geom_point(aes(y=M)) +
			geom_errorbar(aes(ymin=CL, ymax=CU)) +
			coord_flip() +
			theme_bw() +
			labs(x='GP hyperparameter\n', y='')
	ggsave(file=file.path(prjdir,'results_200220','200428f_notsuppAmongInfected_gppars.pdf'), w=6, h=3)
	
	#
	#	extract MVL by gender and location
	ps <- c(0.025,0.5,0.975)
	tmp2 <- sample(seq_len(nrow(reB$p_predict_00)),nrow(re$exp_m_predict_00),replace=TRUE)	
	tmp <- cbind( (1-reB$p_predict_00[tmp2, ]) * re$exp_m_predict_00,
			(1-reB$p_predict_10[tmp2, ]) * re$exp_m_predict_10,
			(1-reB$p_predict_01[tmp2, ]) * re$exp_m_predict_01,
			(1-reB$p_predict_11[tmp2, ]) * re$exp_m_predict_11	)		
	rp <- as.data.table(reshape2::melt(tmp))
	setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, ROW_ID:= seq_len(nrow(tmp))]
	tmp <- merge(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
	rp <- merge(tmp, rp, by=c('ROW_ID'))	
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]
	rp <- rp[, list(Q= quantile(P, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','SEX_LABEL')]
	rp <- dcast.data.table(rp, LOC_LABEL+SEX_LABEL~P, value.var='Q')
	rp[, LABEL:= paste0(sprintf('%2.0f',M),' (',sprintf('%2.0f',CL),'-',sprintf('%2.0f',CU),')') ]
	mvl.by.sex.loc <- copy(rp)
	#	   LOC_LABEL SEX_LABEL        CL         CU         M             LABEL
	#1:   fishing         M 6459.2055 11298.7826 8318.0350 	8318 (6459-11299)
	#2:   fishing         F 2191.5886  4269.4039 2968.2880  2968 (2192-4269)
	#3:    inland         M 1054.5387  1950.9344 1405.3776  1405 (1055-1951)
	#4:    inland         F  559.6118   887.5624  698.1173     698 (560-888)	


	#	extract MVL ratio female:male and male:female
	tmp2 <- sample(seq_len(nrow(reB$p_predict_00)),nrow(re$exp_m_predict_00),replace=TRUE)	
	tmp <- cbind( (1-reB$p_predict_00[tmp2, ]) * re$exp_m_predict_00,
			(1-reB$p_predict_10[tmp2, ]) * re$exp_m_predict_10,
			(1-reB$p_predict_01[tmp2, ]) * re$exp_m_predict_01,
			(1-reB$p_predict_11[tmp2, ]) * re$exp_m_predict_11	)		
	rp <- as.data.table(reshape2::melt( tmp ))
	setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, ROW_ID:= seq_len(nrow(tmp))]
	tmp <- merge(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
	rp <- merge(tmp, rp, by=c('ROW_ID'))	
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]
	rp <- dcast.data.table(rp, LOC_LABEL+iterations~SEX_LABEL, value.var='P')
	rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
	rp <- melt(rp, id.vars=c('LOC_LABEL','iterations'))
	rp <- rp[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','variable')]
	rp <- dcast.data.table(rp, LOC_LABEL+variable~P, value.var='Q')
	rp[, LABEL:= paste0(sprintf('%2.2f',M),' (',sprintf('%2.2f',CL),' - ',sprintf('%2.2f',CU),')') ]	
	mvl.ratio.by.loc <- copy(rp)	
	#	   LOC_LABEL variable        CL        CU         M              LABEL
	#1:   fishing    PR_FM 0.2307099 0.5531012 0.3560961 0.36 (0.23 - 0.55)
	#2:   fishing    PR_MF 1.8079873 4.3344482 2.8082310 2.81 (1.81 - 4.33)
	#3:    inland    PR_FM 0.3349882 0.7226552 0.4958447 0.50 (0.33 - 0.72)
	#4:    inland    PR_MF 1.3837858 2.9851794 2.0167604 2.02 (1.38 - 2.99)


	#
	#	extract MVL by age for gender and location
	ps <- c(0.025,0.5,0.975)
	tmp2 <- sample(seq_len(nrow(reB$p_predict_00)),nrow(re$exp_m_predict_00),replace=TRUE)	
	tmp <- cbind( (1-reB$p_predict_00[tmp2, ]) * re$exp_m_predict_00,
			(1-reB$p_predict_10[tmp2, ]) * re$exp_m_predict_10,
			(1-reB$p_predict_01[tmp2, ]) * re$exp_m_predict_01,
			(1-reB$p_predict_11[tmp2, ]) * re$exp_m_predict_11	)	
	#tmp <- cbind( re$exp_m_predict_00, re$exp_m_predict_10, re$exp_m_predict_01, re$exp_m_predict_11	)	
	tmp <- apply(tmp, 2, quantile, probs=ps)
	rownames(tmp) <- c('CL','M','CU')
	tmp <- as.data.table(reshape2::melt(tmp))	
	tmp <- dcast.data.table(tmp, Var2~Var1, value.var='value')
	mvl.by.age <- cbind( as.data.table(expand.grid(AGE_LABEL= stan.dataB$x_predict, SEX=c(0,1), LOC=c(0,1))), tmp )
	mvl.by.age <- merge(unique(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), mvl.by.age, by=c('LOC','SEX'))
	ggplot(mvl.by.age) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
			geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
			scale_x_continuous( expand=c(0,0) ) + 
			#scale_y_log10() +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			facet_wrap(~LOC_LABEL, ncol=2, scales='free_y') +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='mean viral load\n(95% credibility interval)\n', 
					colour='gender', 
					linetype='location')
	ggsave(file=file.path(prjdir,'results_200220','200428f_mvl_vs_age_by_gender_fishinland_stan.pdf'), w=7, h=5)

	
	#
	#	extract MVL ratio by age for gender and location
	tmp2 <- sample(seq_len(nrow(reB$p_predict_00)),nrow(re$exp_m_predict_00),replace=TRUE)	
	tmp <- cbind( (1-reB$p_predict_00[tmp2, ]) * re$exp_m_predict_00,
			(1-reB$p_predict_10[tmp2, ]) * re$exp_m_predict_10,
			(1-reB$p_predict_01[tmp2, ]) * re$exp_m_predict_01,
			(1-reB$p_predict_11[tmp2, ]) * re$exp_m_predict_11	)	
	rp <- as.data.table(reshape2::melt( tmp ))
	setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, ROW_ID:= seq_len(nrow(tmp))]
	tmp <- merge(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
	rp <- merge(tmp, rp, by=c('ROW_ID'))	
	rp <- dcast.data.table(rp, LOC_LABEL+iterations+AGE_LABEL~SEX_LABEL, value.var='P')
	rp[, PR_FM:= F/M]
	rp[, PR_MF:=M/F]
	rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))
	rp <- rp[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]	
	rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
	rp[, LABEL:= paste0(sprintf('%2.2f',M),' (',sprintf('%2.2f',CL),' - ',sprintf('%2.2f',CU),')') ]	
	mvl.ratio.by.loc.age <- copy(rp)
	
	
		
	save(vla, reB, re, re11, re01, re10, mvl.gp.pars, mvl.by.sex.loc, mvl.ratio.by.loc, mvl.by.age, mvl.ratio.by.loc.age, file=file.path(outdir, "200428f_mvl.rda"))
	
	
	#	make table version suppressed
	mvl.by.age[, LABEL:= paste0(sprintf('%2.0f',M),' (',sprintf('%2.0f',CL),'-',sprintf('%2.0f',CU),')') ]
	set(mvl.by.age, NULL, 'SEX_LABEL', mvl.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])	
	setnames(mvl.ratio.by.loc.age,'LABEL', 'LABEL2')
	dt <- subset(mvl.by.age, AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5))	
	dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
	tmp <- subset(mvl.ratio.by.loc.age, variable=="PR_FM" & AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5), c(LOC_LABEL, AGE_LABEL, LABEL2))
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	tmp <- subset(mvl.ratio.by.loc.age, variable=="PR_MF" & AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5), c(LOC_LABEL, AGE_LABEL, LABEL2))
	setnames(tmp, 'LABEL2', 'PR_MF')
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	write.csv(dt, row.names=FALSE, file=file.path(outdir, "200428f_mvl.csv"))
	
}

vl.suppofinfected.by.gender.loc.age.icar<- function()
{
	require(Hmisc)
	require(data.table)
	require(ggplot2)
	require(rstan)
	vl.detectable <- 4e2
	vl.suppressed <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	outdir <- file.path(prjdir,'results_200220')
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
	
	subset(ds, FC=='inland' & SEX=='F' & AGEYRS==25)
	
	tmp <- seq.int(min(ds$AGEYRS), max(ds$AGEYRS))
	vla <- as.data.table(expand.grid(FC=c('fishing','inland'), SEX=c('M','F'), AGEYRS=tmp))
	vla <- vla[, {		
				z <- which(ds$FC==FC & ds$SEX==SEX & ds$AGEYRS==AGEYRS)	
				list(	N= length(z),
						HIV_N= length(which(ds$HIV_STATUS[z]==1)),
						VLNS_N= length(which(ds$VLNS[z]==1)),
						ARV_N= length(which(ds$ARVMED[z]==0 & ds$HIV_STATUS[z]==1 & !is.na(ds$ARVMED[z])))
				)				
			}, by=c('FC','SEX','AGEYRS')]
	setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
	vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
	vla[, SEX:= as.integer(SEX_LABEL=='M')]
	vla[, AGE:= AGE_LABEL-14L]
	vla[, ROW_ID:= seq_len(nrow(vla))]
	
	
	stan.code2 <- "functions{			
	real icar_normal_lpdf( vector phi, int N, int[] node1, int[] node2)
	{			
		return -0.5 * dot_self(phi[node1] - phi[node2]) + normal_lpdf( sum(phi) | 0, 0.001 * N);
	}
}
			
data{
	int<lower=1> N; 
	int<lower=0> TOTAL[N];    
	int<lower=0> K[N];
	int<lower=0> AGE_N;  
	int<lower=0, upper=AGE_N> AGE[N];  
	int<lower=0,upper=1> SEX[N];
	int<lower=0,upper=1> LOC[N];	 
	int<lower=0> N_edges; 						// number of related age groups
	int<lower=1, upper=AGE_N> node1[N_edges];  	// node1[i] adjacent to node2[i]
	int<lower=1, upper=AGE_N> node2[N_edges];  	// and node1[i] < node2[i]
}
			
parameters{
	real sex0_loc0;
	real sex1_loc0;
	real sex0_loc1;
	real sex1_loc1;
	vector[AGE_N] phi_sex0_loc0;
	vector[AGE_N] phi_sex1_loc0;
	vector[AGE_N] phi_sex0_loc1;
	vector[AGE_N] phi_sex1_loc1; 			 		
	real<lower=0> sigma_loc0;
	real<lower=0> sigma_loc1;
}
			
transformed parameters{
	vector[N] p_logit;
	for ( i in 1:N ) 
	{
		p_logit[i] = sex0_loc0 * (1-SEX[i]) * (1-LOC[i]) + 
		sex1_loc0 * SEX[i] * (1-LOC[i]) + 
		sex0_loc1 * (1-SEX[i]) * LOC[i] + 
		sex1_loc1 * SEX[i] * LOC[i] +
		sigma_loc0 * phi_sex0_loc0[AGE[i]] * (1-SEX[i]) * (1-LOC[i]) +
		sigma_loc0 * phi_sex1_loc0[AGE[i]] * SEX[i] * (1-LOC[i]) +
		sigma_loc1 * phi_sex0_loc1[AGE[i]] * (1-SEX[i]) * LOC[i] +
		sigma_loc1 * phi_sex1_loc1[AGE[i]] * SEX[i] * LOC[i];
	}
}
			
model{	
	sex0_loc0 ~ normal( 0 , 10 );
	sex1_loc0 ~ normal( 0 , 10 );
	sex0_loc1 ~ normal( 0 , 10 );
	sex1_loc1 ~ normal( 0 , 10 );
	phi_sex0_loc0 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	phi_sex1_loc0 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	phi_sex0_loc1 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	phi_sex1_loc1 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	sigma_loc0 ~ normal(0.0, 1);
	sigma_loc1 ~ normal(0.0, 1);
	K ~ binomial_logit( TOTAL , p_logit );
}
			
generated quantities{
	vector[N] p;
	p= inv_logit(p_logit);
}"
	
	stan.model2 <- stan_model(model_name= 'icar_age_interactions',model_code = gsub('\t',' ',stan.code2))
	stan.data <- list()
	stan.data$N <- nrow(vla)
	stan.data$TOTAL <- vla[,HIV_N]
	stan.data$K <- vla[,VLNS_N]
	stan.data$AGE_N <- vla[, max(AGE)]
	stan.data$AGE <- vla[, AGE]
	stan.data$SEX <- vla[, SEX]
	stan.data$LOC <- vla[, LOC]
	#	second order RW prior
	stan.data$node1 <-  c(vla[, seq.int(1, max(AGE)-1L)], vla[, seq.int(1, max(AGE)-2L)])
	stan.data$node2 <-  c(vla[, seq.int(2, max(AGE))], vla[, seq.int(3, max(AGE))])
	tmp <- sort(stan.data$node1, index.return=TRUE)$ix
	stan.data$node1 <- stan.data$node1[tmp]
	stan.data$node2 <- stan.data$node2[tmp]
	stan.data$N_edges <-  length(stan.data$node1)
	fit <- sampling(stan.model2, data=stan.data, iter=20e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999) )	
	#save(fit, file=file.path(outdir, "hivprevalence_icar_stanfit_200428.rda"))		# trends by age quite rough, using Cauchy prior on sigma
	#save(fit, file=file.path(outdir, "hivprevalence_icar_stanfit_200428c.rda"))	# trends by age still quite rough, using N(0,0.1) prior on sigma
	save(fit, file=file.path(outdir, "notsuppAmongInfected_icar_stanfit_200428.rda"))
	
	#
	#	compare to self-report
	#
	stan.data <- list()
	stan.data$N <- nrow(vla)
	stan.data$TOTAL <- vla[,HIV_N]
	stan.data$K <- vla[,ARV_N]
	stan.data$AGE_N <- vla[, max(AGE)]
	stan.data$AGE <- vla[, AGE]
	stan.data$SEX <- vla[, SEX]
	stan.data$LOC <- vla[, LOC]
	#	second order RW prior
	stan.data$node1 <-  c(vla[, seq.int(1, max(AGE)-1L)], vla[, seq.int(1, max(AGE)-2L)])
	stan.data$node2 <-  c(vla[, seq.int(2, max(AGE))], vla[, seq.int(3, max(AGE))])
	tmp <- sort(stan.data$node1, index.return=TRUE)$ix
	stan.data$node1 <- stan.data$node1[tmp]
	stan.data$node2 <- stan.data$node2[tmp]
	stan.data$N_edges <-  length(stan.data$node1)
	fit2 <- sampling(stan.model2, data=stan.data, iter=20e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999) )	
	save(fit2, file=file.path(outdir, "notARVAmongInfected_icar_stanfit_200428.rda"))
	
	
	re <- rstan::extract(fit)
	re2 <- rstan::extract(fit2)
	ps <- c(0.025,0.5,0.975)
	
	
	#	make prevalence plot by age
	tmp <- apply(re$p, 2, quantile, probs=ps)
	rownames(tmp) <- c('CL','M','CU')
	tmp <- as.data.table(reshape2::melt(tmp))	
	nsprev.by.age <- cbind(vla, dcast.data.table(tmp, Var2~Var1, value.var='value'))
	nsprev.by.age[, STAT:='VLNS']
	tmp <- apply(re2$p, 2, quantile, probs=ps)
	rownames(tmp) <- c('CL','M','CU')
	tmp <- as.data.table(reshape2::melt(tmp))	
	naprev.by.age <- cbind(vla, dcast.data.table(tmp, Var2~Var1, value.var='value'))
	naprev.by.age[, STAT:='VLNA']
	tmp <- subset(naprev.by.age, select=c(ROW_ID, M, CL, CU))
	setnames(tmp, c('M','CL','CU'), c('M2','CL2','CU2'))
	tmp <- merge(nsprev.by.age, tmp, by='ROW_ID')	
	ggplot(tmp) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
			geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			facet_wrap(~LOC_LABEL, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='HIV+ individuals with unsuppressed viral load\n(95% credibility interval)\n', 
					colour='gender')
	ggsave(file=file.path(prjdir,'results_200220','200428d_notsuppAmongInfected_vs_age_by_gender_fishinland_stan_v1.pdf'), w=6, h=5)	
	ggplot(tmp) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
			geom_hline(yintercept=c(1-0.9^3, 1-0.95^3)) +
			geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
			geom_line(aes(x=AGE_LABEL, y=M2, colour=SEX_LABEL), linetype=2) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			facet_wrap(~LOC_LABEL, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='HIV+ individuals with unsuppressed viral load\n(95% credibility interval)\n', 
					colour='gender')
	ggsave(file=file.path(prjdir,'results_200220','200428d_notsuppAmongInfected_vs_age_by_gender_fishinland_stan_v2.pdf'), w=6, h=5)	
	ggplot(tmp) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=1-CU, ymax=1-CL, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
			geom_hline(yintercept=c(0.9^3, 0.95^3)) +
			geom_line(aes(x=AGE_LABEL, y=1-M, colour=SEX_LABEL)) +
			geom_line(aes(x=AGE_LABEL, y=1-M2, colour=SEX_LABEL), linetype=2) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			facet_wrap(~LOC_LABEL, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='HIV+ individuals with suppressed viral load\n(95% credibility interval)\n', 
					colour='gender')
	ggsave(file=file.path(prjdir,'results_200220','200428d_suppAmongInfected_vs_age_by_gender_fishinland_stan_v2.pdf'), w=6, h=5)
	tmp <- rbind(nsprev.by.age, naprev.by.age, fill=TRUE)
	ggplot(tmp) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL,STAT)), alpha=0.2) +
			geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL, linetype=STAT)) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			scale_linetype_manual(values=c('VLNS'='solid','VLNA'='dotdash')) +
			facet_grid(SEX_LABEL~LOC_LABEL) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='HIV+ individuals with unsuppressed viral load\n(95% credibility interval)\n', 
					colour='gender')
	ggsave(file=file.path(prjdir,'results_200220','200428d_notsuppAmongInfected_vs_age_by_gender_fishinland_stan_v3.pdf'), w=9, h=8)
	ggplot(tmp) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=1-CU, ymax=1-CL, group=interaction(SEX_LABEL,LOC_LABEL,STAT)), alpha=0.2) +
			geom_line(aes(x=AGE_LABEL, y=1-M, colour=SEX_LABEL, linetype=STAT)) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			scale_linetype_manual(values=c('VLNS'='solid','VLNA'='dotdash')) +
			facet_grid(SEX_LABEL~LOC_LABEL) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='HIV+ individuals with suppressed viral load\n(95% credibility interval)\n', 
					colour='gender')
	ggsave(file=file.path(prjdir,'results_200220','200428d_suppAmongInfected_vs_age_by_gender_fishinland_stan_v3.pdf'), w=9, h=8)
	
	
	#	extract basic not supp estimates
	vla[, list(N=sum(VLNS_N), P=sum(VLNS_N) / sum(HIV_N), N2=sum(HIV_N)-sum(VLNS_N), P2= 1-sum(VLNS_N) / sum(HIV_N)), by=c('LOC_LABEL','SEX_LABEL')]
	#   LOC_LABEL SEX_LABEL   N         P   N2        P2
	#1:   fishing         M 296 0.4314869  390 0.5685131
	#2:    inland         M 203 0.3338816  405 0.6661184
	#3:   fishing         F 198 0.2307692  660 0.7692308
	#4:    inland         F 279 0.2121673 1036 0.7878327

	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
	rp <- rp[, list(Q= quantile(P, probs=ps), P=c('CL','M','CU')), by=c('LOC','SEX')]
	rp <- dcast.data.table(rp, LOC+SEX~P, value.var='Q')
	rp[, LABEL:= paste0(round(M, d=2),'% (',round(CL, d=2),'% - ',round(CU,d=2),'%)') ]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
	nsprev.by.sex.loc <- copy(rp)
	#	   LOC SEX LOC_LABEL SEX_LABEL         CL         CU         M                 LABEL
	#1:   0   0    inland         F 0.2363155 0.3032974 0.2691669  0.27% (0.24% - 0.3%)
	#2:   0   1    inland         M 0.3499670 0.4977837 0.4245407  0.42% (0.35% - 0.5%)
	#3:   1   0   fishing         F 0.2191378 0.2814667 0.2498244 0.25% (0.22% - 0.28%)
	#4:   1   1   fishing         M 0.4360992 0.5266300 0.4826315 0.48% (0.44% - 0.53%)
	
	#	extract risk ratio of unsuppressed VL female:male and male:female
	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
	rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
	rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC','iterations')]
	rp <- melt(rp, id.vars=c('LOC','iterations'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC','variable')]
	rp <- dcast.data.table(rp, LOC+variable~P, value.var='Q')
	rp[, LABEL:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL))), rp, by=c('LOC'))
	nsprev.by.loc <- copy(rp)
	#	   LOC LOC_LABEL variable        CL        CU         M              LABEL
	#1:   0    inland    PR_FM 0.5171747 0.7940446 0.6343892 0.63 (0.52 - 0.79)
	#2:   0    inland    PR_MF 1.2593751 1.9335827 1.5763195 1.58 (1.26 - 1.93)
	#3:   1   fishing    PR_FM 0.4424546 0.6038996 0.5179708  0.52 (0.44 - 0.6)
	#4:   1   fishing    PR_MF 1.6559043 2.2601191 1.9306106 1.93 (1.66 - 2.26)
	

	#	extract risk ratio of unsuppressed VL female:male and male:female by age
	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations+AGE+AGE_LABEL~SEX_LABEL, value.var='P')
	rp[, PR_FM:= F/M]
	rp[, PR_MF:=M/F]
	rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]
	rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
	nsprev.ratio.by.loc.age <- copy(rp)


	#	extract risk ratio of suppressed VL female:male and male:female by age
	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations+AGE+AGE_LABEL~SEX_LABEL, value.var='P')
	rp[, PR_FM:= (1-F)/(1-M)]
	rp[, PR_MF:=(1-M)/(1-F)]
	rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]
	rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
	sprev.ratio.by.loc.age <- copy(rp)
	
	
	#	extract if difference in female:male risk ratio of unsuppressed VL is different in fishing vs inland
	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
	rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
	rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
	rp <- dcast.data.table(rp, iterations~LOC_LABEL, value.var='PR_FM')
	rp[, PR_FM_D:= fishing-inland]	
	rp[, list(Q= quantile(PR_FM_D, probs=ps), P=c('CL','M','CU'))]
	#	            Q  P
	#1: -0.2940328 CL
	#2: -0.1167300  M
	#3:  0.0256620 CU
	

	save(vla, re, sprev.ratio.by.loc.age, nsprev.by.age, naprev.by.age, nsprev.by.loc, nsprev.by.sex.loc, nsprev.ratio.by.loc.age, file=file.path(outdir, "notsuppamonginfected_200428.rda"))
	
	#	make table version unsuppressed
	nsprev.by.age[, LABEL:= paste0(round(M*100, d=1),' (',round(CL*100, d=1),' - ',round(CU*100,d=1),')') ]
	set(nsprev.by.age, NULL, 'SEX_LABEL', nsprev.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])	
	nsprev.ratio.by.loc.age[, LABEL2:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
	dt <- subset(nsprev.by.age, AGE_LABEL%in%c(20,25,30,35,40,45))	
	dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
	tmp <- subset(nsprev.ratio.by.loc.age, variable=="PR_FM" & AGE_LABEL%in%c(20,25,30,35,40,45), c(LOC_LABEL, AGE_LABEL, LABEL2))
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	tmp <- subset(nsprev.ratio.by.loc.age, variable=="PR_MF" & AGE_LABEL%in%c(20,25,30,35,40,45), c(LOC_LABEL, AGE_LABEL, LABEL2))
	setnames(tmp, 'LABEL2', 'PR_MF')
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	write.csv(dt, row.names=FALSE, file=file.path(outdir, "notsuppamonginfected_200428.csv"))
		
	#	make table version suppressed
	nsprev.by.age[, LABEL:= paste0(round((1-M)*100, d=1),' (',round((1-CU)*100, d=1),' - ',round((1-CL)*100,d=1),')') ]
	set(nsprev.by.age, NULL, 'SEX_LABEL', nsprev.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])	
	sprev.ratio.by.loc.age[, LABEL2:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
	dt <- subset(nsprev.by.age, AGE_LABEL%in%c(20,25,30,35,40,45))	
	dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
	tmp <- subset(sprev.ratio.by.loc.age, variable=="PR_FM" & AGE_LABEL%in%c(20,25,30,35,40,45), c(LOC_LABEL, AGE_LABEL, LABEL2))
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	tmp <- subset(sprev.ratio.by.loc.age, variable=="PR_MF" & AGE_LABEL%in%c(20,25,30,35,40,45), c(LOC_LABEL, AGE_LABEL, LABEL2))
	setnames(tmp, 'LABEL2', 'PR_MF')
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	write.csv(dt, row.names=FALSE, file=file.path(outdir, "suppamonginfected_200428.csv"))	
}

vl.suppofinfected.by.gender.loc.age.gp<- function()
{
	require(Hmisc)
	require(data.table)
	require(ggplot2)
	require(rstan)
	vl.detectable <- 4e2
	vl.suppressed <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	outdir <- file.path(prjdir,'results_200220')
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
	
	subset(ds, FC=='inland' & SEX=='F' & AGEYRS==25)
	
	tmp <- seq.int(min(ds$AGEYRS), max(ds$AGEYRS))
	vla <- as.data.table(expand.grid(FC=c('fishing','inland'), SEX=c('M','F'), AGEYRS=tmp))
	vla <- vla[, {		
				z <- which(ds$FC==FC & ds$SEX==SEX & ds$AGEYRS==AGEYRS)	
				list(	N= length(z),
						HIV_N= length(which(ds$HIV_STATUS[z]==1)),
						VLNS_N= length(which(ds$VLNS[z]==1)),
						ARV_N= length(which(ds$ARVMED[z]==0 & ds$HIV_STATUS[z]==1 & !is.na(ds$ARVMED[z])))
				)				
			}, by=c('FC','SEX','AGEYRS')]
	setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
	vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
	vla[, SEX:= as.integer(SEX_LABEL=='M')]
	vla[, AGE:= AGE_LABEL-14L]
	vla[, ROW_ID:= seq_len(nrow(vla))]
	
	stan.code2 <- "
data{	
	int<lower=1> N_predict;
  	real x_predict[N_predict];
  	int<lower=1> N_observed;
  	int<lower=1, upper=N_predict> observed_idx[N_observed];
  	int y_observed_00[N_observed];
	int y_observed_10[N_observed];
	int y_observed_01[N_observed];	
	int y_observed_11[N_observed];	
	int total_observed_00[N_observed];
	int total_observed_10[N_observed];
	int total_observed_01[N_observed];
	int total_observed_11[N_observed];
  	real<lower=0> rho_hyper_par_00;
	real<lower=0> rho_hyper_par_10;
	real<lower=0> rho_hyper_par_01;
	real<lower=0> rho_hyper_par_11;  	
  	real<lower=0> alpha_hyper_par_00;
	real<lower=0> alpha_hyper_par_10;
	real<lower=0> alpha_hyper_par_01;
	real<lower=0> alpha_hyper_par_11;
}

parameters {
	real<lower=0> rho_00;
	real<lower=0> rho_10;
	real<lower=0> rho_01;
	real<lower=0> rho_11;
	real<lower=0> alpha_00;
	real<lower=0> alpha_10;
	real<lower=0> alpha_01;
	real<lower=0> alpha_11;
	real sex0_loc0;
	real sex1_loc0;
	real sex0_loc1;
	real sex1_loc1;  	
  	vector[N_predict] f_tilde_00;
	vector[N_predict] f_tilde_10;
	vector[N_predict] f_tilde_01;
	vector[N_predict] f_tilde_11;
}

transformed parameters {
  	matrix[N_predict, N_predict] L_cov;
  	vector[N_predict] logit_p_predict_00;
	vector[N_predict] logit_p_predict_10;
	vector[N_predict] logit_p_predict_01;
	vector[N_predict] logit_p_predict_11;
	// GP for 00 and 01 (women)
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_00, rho_00) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_00 = sex0_loc0 + L_cov * f_tilde_00;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_01, rho_01) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_01 = sex0_loc1 + L_cov * f_tilde_01;
	// GP for 10 and 10 (men)
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_10, rho_10) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_10 = sex1_loc0 + L_cov * f_tilde_10;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_11, rho_11) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_11 = sex1_loc1 + L_cov * f_tilde_11;
}

model {
	rho_00 ~ normal(0, rho_hyper_par_00);
	rho_10 ~ normal(0, rho_hyper_par_10);
	rho_01 ~ normal(0, rho_hyper_par_01);
	rho_11 ~ normal(0, rho_hyper_par_11);  	
  	alpha_00 ~ normal(0, alpha_hyper_par_00);
	alpha_10 ~ normal(0, alpha_hyper_par_10);
	alpha_01 ~ normal(0, alpha_hyper_par_01);
	alpha_11 ~ normal(0, alpha_hyper_par_11);
  	sex0_loc0 ~ normal( 0 , 10 );
	sex0_loc1 ~ normal( 0 , 10 );
	sex1_loc0 ~ normal( 0 , 10 );
	sex1_loc1 ~ normal( 0 , 10 );
  	f_tilde_00 ~ normal(0, 1);
	f_tilde_01 ~ normal(0, 1);
	f_tilde_10 ~ normal(0, 1);
	f_tilde_11 ~ normal(0, 1);
  	y_observed_00 ~ binomial_logit(total_observed_00, logit_p_predict_00[observed_idx] );
	y_observed_01 ~ binomial_logit(total_observed_01, logit_p_predict_01[observed_idx] );
	y_observed_10 ~ binomial_logit(total_observed_10, logit_p_predict_10[observed_idx] );
	y_observed_11 ~ binomial_logit(total_observed_11, logit_p_predict_11[observed_idx] );
}

generated quantities {
  	vector[N_predict] p_predict_00;
	vector[N_predict] p_predict_01;
	vector[N_predict] p_predict_10;
	vector[N_predict] p_predict_11;
	p_predict_00 = inv_logit(logit_p_predict_00);
	p_predict_01 = inv_logit(logit_p_predict_01);  
	p_predict_10 = inv_logit(logit_p_predict_10);
	p_predict_11 = inv_logit(logit_p_predict_11);
}			
"
		
	stan.model <- stan_model(model_name= 'gp_all',model_code = gsub('\t',' ',stan.code2))	
	stan.data <- list()
	stan.data$x_predict <- seq(vla[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
	stan.data$N_predict <- length(stan.data$x_predict)
	stan.data$observed_idx <- which(stan.data$x_predict%%1==0.5)
	stan.data$N_observed <- length(stan.data$observed_idx)
	stan.data$y_observed_00 <- vla[SEX==0 & LOC==0, HIV_N-VLNS_N]
	stan.data$y_observed_10 <- vla[SEX==1 & LOC==0, HIV_N-VLNS_N]
	stan.data$y_observed_01 <- vla[SEX==0 & LOC==1, HIV_N-VLNS_N]
	stan.data$y_observed_11 <- vla[SEX==1 & LOC==1, HIV_N-VLNS_N]
	stan.data$total_observed_00 <- vla[SEX==0 & LOC==0, HIV_N]
	stan.data$total_observed_10 <- vla[SEX==1 & LOC==0, HIV_N]
	stan.data$total_observed_01 <- vla[SEX==0 & LOC==1, HIV_N]
	stan.data$total_observed_11 <- vla[SEX==1 & LOC==1, HIV_N]
	stan.data$rho_hyper_par_00 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_10 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_01 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_11 <- diff(range(stan.data$x_predict))/3
	stan.data$alpha_hyper_par_00 <- 2
	stan.data$alpha_hyper_par_10 <- 2
	stan.data$alpha_hyper_par_01 <- 2
	stan.data$alpha_hyper_par_11 <- 2
	fit <- sampling(stan.model, data=stan.data, iter=10e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
	save(fit, file=file.path(outdir, "200428f_notsuppAmongInfected_gp_stanfit.rda"))
	
	
	#
	#	compare to self-report
	#
	stan.data <- list()
	stan.data$x_predict <- seq(vla[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
	stan.data$N_predict <- length(stan.data$x_predict)
	stan.data$observed_idx <- which(stan.data$x_predict%%1==0.5)
	stan.data$N_observed <- length(stan.data$observed_idx)
	stan.data$y_observed_00 <- vla[SEX==0 & LOC==0, HIV_N-ARV_N]
	stan.data$y_observed_10 <- vla[SEX==1 & LOC==0, HIV_N-ARV_N]
	stan.data$y_observed_01 <- vla[SEX==0 & LOC==1, HIV_N-ARV_N]
	stan.data$y_observed_11 <- vla[SEX==1 & LOC==1, HIV_N-ARV_N]
	stan.data$total_observed_00 <- vla[SEX==0 & LOC==0, HIV_N]
	stan.data$total_observed_10 <- vla[SEX==1 & LOC==0, HIV_N]
	stan.data$total_observed_01 <- vla[SEX==0 & LOC==1, HIV_N]
	stan.data$total_observed_11 <- vla[SEX==1 & LOC==1, HIV_N]
	stan.data$rho_hyper_par_00 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_10 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_01 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_11 <- diff(range(stan.data$x_predict))/3
	stan.data$alpha_hyper_par_00 <- 2
	stan.data$alpha_hyper_par_10 <- 2
	stan.data$alpha_hyper_par_01 <- 2
	stan.data$alpha_hyper_par_11 <- 2
	fit2 <- sampling(stan.model, data=stan.data, iter=10e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
	save(fit2, file=file.path(outdir, "200428f_notARVAmongInfected_icar_stanfit.rda"))
	
	
	re <- rstan::extract(fit)
	re2 <- rstan::extract(fit2)
	ps <- c(0.025,0.5,0.975)
	tmp <- summary(fit)$summary
	tmp[grepl('^p_predict_',rownames(tmp)),]
	
	#
	#	extract hyperparams rho
	ps <- c(0.025,0.25,0.5,0.75,0.975)
	tmp <- cbind( quantile(re$rho_00, probs=ps),
			quantile(re$rho_10, probs=ps),
			quantile(re$rho_01, probs=ps),
			quantile(re$rho_11, probs=ps),
			quantile(re$alpha_00, probs=ps),
			quantile(re$alpha_10, probs=ps),
			quantile(re$alpha_01, probs=ps),
			quantile(re$alpha_11, probs=ps) )			
	colnames(tmp) <- c('rho_00','rho_10','rho_01','rho_11','alpha_00','alpha_10','alpha_01','alpha_11')
	rownames(tmp) <- c('CL','IL','M','IU','CU')
	tmp <- as.data.table(reshape2::melt(tmp))
	setnames(tmp, 'Var2', 'GP_hyper_par')
	tmp[, SEX:= as.integer(gsub('^([a-z]+)_([0-9])([0-9])','\\2',GP_hyper_par))]
	tmp[, LOC:= as.integer(gsub('^([a-z]+)_([0-9])([0-9])','\\3',GP_hyper_par))]
	tmp[, GP_hyper_par:= gsub('^([a-z]+)_([0-9])([0-9])','\\1',GP_hyper_par)]
	tmp <- dcast.data.table(tmp, LOC+SEX+GP_hyper_par~Var1, value.var='value')
	nsinf.gp.pars <- merge(unique(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), tmp, by=c('SEX','LOC'))
	ggplot(nsinf.gp.pars, aes(x= paste0(GP_hyper_par, ' ', LOC_LABEL, ' ', SEX_LABEL))) +
			geom_point(aes(y=M)) +
			geom_errorbar(aes(ymin=CL, ymax=CU)) +
			coord_flip() +
			theme_bw() +
			labs(x='GP hyperparameter\n', y='')
	ggsave(file=file.path(prjdir,'results_200220','200428f_notsuppAmongInfected_gppars.pdf'), w=6, h=3)
	
	
	#
	#	make prevalence plot by age
	tmp <- cbind( apply(re$p_predict_00, 2, quantile, probs=ps),
			apply(re$p_predict_10, 2, quantile, probs=ps),
			apply(re$p_predict_01, 2, quantile, probs=ps),
			apply(re$p_predict_11, 2, quantile, probs=ps)
			)
	rownames(tmp) <- c('CL','IL','M','IU','CU')
	tmp <- as.data.table(reshape2::melt(tmp))
	nsinf.by.age <- dcast.data.table(tmp, Var2~Var1, value.var='value')
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	nsinf.by.age <- cbind(tmp, nsinf.by.age) 
	nsinf.by.age <- merge(unique(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), nsinf.by.age, by=c('SEX','LOC'))
	nsinf.by.age[, STAT:='VLNS']
	
	tmp <- cbind( apply(re2$p_predict_00, 2, quantile, probs=ps),
			apply(re2$p_predict_10, 2, quantile, probs=ps),
			apply(re2$p_predict_01, 2, quantile, probs=ps),
			apply(re2$p_predict_11, 2, quantile, probs=ps)
			)
	rownames(tmp) <- c('CL','IL','M','IU','CU')
	tmp <- as.data.table(reshape2::melt(tmp))
	nainf.by.age <- dcast.data.table(tmp, Var2~Var1, value.var='value')
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	nainf.by.age <- cbind(tmp, nainf.by.age) 
	nainf.by.age <- merge(unique(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), nainf.by.age, by=c('SEX','LOC'))
	nainf.by.age[, STAT:='VLNA']
	tmp <- subset(nainf.by.age, select=c(SEX, LOC, AGE_LABEL, M, CL, CU))
	setnames(tmp, c('M','CL','CU'), c('M2','CL2','CU2'))
	tmp <- merge(nsinf.by.age, tmp, by=c('SEX','LOC','AGE_LABEL'))			
	ggplot(tmp) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
			geom_hline(yintercept=c(0.9^3, 0.95^3)) +
			geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			facet_wrap(~LOC_LABEL, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='HIV+ individuals with suppressed viral load\n(95% credibility interval)\n', 
					colour='gender', 
					linetype='location')
	ggsave(file=file.path(prjdir,'results_200220','200428f_suppAmongInfected_vs_age_by_gender_fishinland_stan.pdf'), w=6, h=5)		
	ggplot(tmp) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +
			geom_hline(yintercept=c(0.9^3, 0.95^3)) +
			geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
			geom_line(aes(x=AGE_LABEL, y=M2, colour=SEX_LABEL), linetype=2) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			facet_wrap(~LOC_LABEL, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='HIV+ individuals with suppressed viral load\n(95% credibility interval)\n', 
					colour='gender')
	ggsave(file=file.path(prjdir,'results_200220','200428f_suppAmongInfected_vs_age_by_gender_fishinland_stan_v2.pdf'), w=6, h=5)	
	tmp <- rbind(nsinf.by.age, nainf.by.age, fill=TRUE)
	ggplot(tmp) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL,STAT)), alpha=0.2) +
			geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL, linetype=STAT)) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			scale_linetype_manual(values=c('VLNS'='solid','VLNA'='dotdash')) +
			facet_grid(SEX_LABEL~LOC_LABEL) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='HIV+ individuals with unsuppressed viral load\n(95% credibility interval)\n', 
					colour='gender')
	ggsave(file=file.path(prjdir,'results_200220','200428f_suppAmongInfected_vs_age_by_gender_fishinland_stan_v3.pdf'), w=9, h=8)
	
	
	#	extract basic not supp estimates
	vla[, list(N=sum(VLNS_N), P=sum(VLNS_N) / sum(HIV_N), N2=sum(HIV_N)-sum(VLNS_N), P2= 1-sum(VLNS_N) / sum(HIV_N)), by=c('LOC_LABEL','SEX_LABEL')]
	#   LOC_LABEL SEX_LABEL   N         P   N2        P2
	#1:   fishing         M 296 0.4314869  390 0.5685131
	#2:    inland         M 203 0.3338816  405 0.6661184
	#3:   fishing         F 198 0.2307692  660 0.7692308
	#4:    inland         F 279 0.2121673 1036 0.7878327
	

	ps <- c(0.025,0.5,0.975)
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt( tmp ))
	setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, ROW_ID:= seq_len(nrow(tmp))]
	tmp <- merge(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
	rp <- merge(tmp, rp, by=c('ROW_ID'))	
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]
	rp <- rp[, list(Q= quantile(P, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','SEX_LABEL')]
	rp <- dcast.data.table(rp, LOC_LABEL+SEX_LABEL~P, value.var='Q')
	rp[, LABEL:= paste0(round(M*100, d=2),'% (',round(CL*100, d=2),'% - ',round(CU*100,d=2),'%)') ]
	nsinf.by.sex.loc <- copy(rp)
	# LOC_LABEL SEX_LABEL        CL        CU         M                    LABEL
	#	fishing         M 0.4562322 0.5536721 0.5018376 50.18% (45.62% - 55.37%)
	#2: fishing         F 0.7134543 0.7781516 0.7469707  74.7% (71.35% - 77.82%)
	#3: inland         M 0.4910456 0.6555350 0.5703211  57.03% (49.1% - 65.55%)
	#4: inland         F 0.6939593 0.7646014 0.7308165  73.08% (69.4% - 76.46%)



	#	extract risk ratio of suppressed VL female:male and male:female
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt( tmp ))
	setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, ROW_ID:= seq_len(nrow(tmp))]
	tmp <- merge(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
	rp <- merge(tmp, rp, by=c('ROW_ID'))	
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]
	rp <- dcast.data.table(rp, LOC_LABEL+iterations~SEX_LABEL, value.var='P')
	rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
	rp <- melt(rp, id.vars=c('LOC_LABEL','iterations'))
	rp <- rp[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','variable')]
	rp <- dcast.data.table(rp, LOC_LABEL+variable~P, value.var='Q')
	rp[, LABEL:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]	
	nsinf.by.loc <- copy(rp)	
	#   LOC_LABEL variable        CL        CU         M              LABEL
	#1:   fishing    PR_FM 1.3365856 1.6508029 1.4873708 1.49 (1.34 - 1.65)
	#2:   fishing    PR_MF 0.6057659 0.7481750 0.6723273 0.67 (0.61 - 0.75)
	#3:    inland    PR_FM 1.1015949 1.4994209 1.2811266   1.28 (1.1 - 1.5)
	#4:    inland    PR_MF 0.6669241 0.9077748 0.7805630 0.78 (0.67 - 0.91)
	
	
	#	extract risk ratio of unsuppressed VL female:male and male:female by age
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt( tmp ))
	setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, ROW_ID:= seq_len(nrow(tmp))]
	tmp <- merge(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
	rp <- merge(tmp, rp, by=c('ROW_ID'))	
	rp <- dcast.data.table(rp, LOC_LABEL+iterations+AGE_LABEL~SEX_LABEL, value.var='P')
	rp[, PR_FM:= F/M]
	rp[, PR_MF:=M/F]
	rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))
	rp <- rp[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]	
	rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
	rp[, LABEL:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]	
	nsinf.ratio.by.loc.age <- copy(rp)
	
	
	#	extract if difference in female:male risk ratio of unsuppressed VL is different in fishing vs inland
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt( tmp ))
	setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, ROW_ID:= seq_len(nrow(tmp))]
	tmp <- merge(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
	rp <- merge(tmp, rp, by=c('ROW_ID'))	
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]
	rp <- dcast.data.table(rp, LOC_LABEL+iterations~SEX_LABEL, value.var='P')
	rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
	rp <- dcast.data.table(rp, iterations~LOC_LABEL, value.var='PR_FM')
	rp[, PR_FM_D:= fishing-inland]	
	rp[, list(Q= quantile(PR_FM_D, probs=ps), P=c('CL','M','CU'))]
	#	            Q  P
	#1: 1: -0.05655754 CL
	#2:  0.20449222  M
	#3:  0.45032189 CU
	
	
	save(vla, re, re2, nainf.by.age, nsinf.by.age, nsinf.by.sex.loc, nsinf.by.loc, nsinf.ratio.by.loc.age, file=file.path(outdir, "200428f_suppAmongInfected.rda"))
	
	
	#	make table version suppressed
	nsinf.by.age[, LABEL:= paste0(sprintf('%2.1f',M*100),' (',sprintf('%2.1f',CL*100),' - ',sprintf('%2.1f',CU*100),')') ]
	set(nsinf.by.age, NULL, 'SEX_LABEL', nsinf.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])	
	nsinf.ratio.by.loc.age[, LABEL2:= paste0(sprintf('%2.2f',M),' (',sprintf('%2.2f',CL),' - ',sprintf('%2.2f',CU),')') ]
	dt <- subset(nsinf.by.age, AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5))	
	dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
	tmp <- subset(nsinf.ratio.by.loc.age, variable=="PR_FM" & AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5), c(LOC_LABEL, AGE_LABEL, LABEL2))
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	tmp <- subset(nsinf.ratio.by.loc.age, variable=="PR_MF" & AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5), c(LOC_LABEL, AGE_LABEL, LABEL2))
	setnames(tmp, 'LABEL2', 'PR_MF')
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	write.csv(dt, row.names=FALSE, file=file.path(outdir, "200428f_suppamonginfected.csv"))	
}

vl.suppofpop.by.gender.loc.age.gp<- function()
{
	require(Hmisc)
	require(data.table)
	require(ggplot2)
	require(rstan)
	vl.detectable <- 4e2
	vl.suppressed <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	outdir <- file.path(prjdir,'results_200220')
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
	
	subset(ds, FC=='inland' & SEX=='F' & AGEYRS==25)
	
	tmp <- seq.int(min(ds$AGEYRS), max(ds$AGEYRS))
	vla <- as.data.table(expand.grid(FC=c('fishing','inland'), SEX=c('M','F'), AGEYRS=tmp))
	vla <- vla[, {		
				z <- which(ds$FC==FC & ds$SEX==SEX & ds$AGEYRS==AGEYRS)	
				list(	N= length(z),
						HIV_N= length(which(ds$HIV_STATUS[z]==1)),
						VLNS_N= length(which(ds$VLNS[z]==1)),
						ARV_N= length(which(ds$ARVMED[z]==0 & ds$HIV_STATUS[z]==1 & !is.na(ds$ARVMED[z])))
				)				
			}, by=c('FC','SEX','AGEYRS')]
	setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
	vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
	vla[, SEX:= as.integer(SEX_LABEL=='M')]
	vla[, AGE:= AGE_LABEL-14L]
	vla[, ROW_ID:= seq_len(nrow(vla))]
	
	
	stan.code <- "
data{	
	int<lower=1> N_predict;
  	real x_predict[N_predict];
  	int<lower=1> N_observed;
  	int<lower=1, upper=N_predict> observed_idx[N_observed];
  	int y_observed_00[N_observed];
	int y_observed_10[N_observed];
	int y_observed_01[N_observed];	
	int y_observed_11[N_observed];	
	int total_observed_00[N_observed];
	int total_observed_10[N_observed];
	int total_observed_01[N_observed];
	int total_observed_11[N_observed];
  	real<lower=0> rho_hyper_par_00;
	real<lower=0> rho_hyper_par_10;
	real<lower=0> rho_hyper_par_01;
	real<lower=0> rho_hyper_par_11;  	
  	real<lower=0> alpha_hyper_par_00;
	real<lower=0> alpha_hyper_par_10;
	real<lower=0> alpha_hyper_par_01;
	real<lower=0> alpha_hyper_par_11;
}

parameters {
	real<lower=0> rho_00;
	real<lower=0> rho_10;
	real<lower=0> rho_01;
	real<lower=0> rho_11;
	real<lower=0> alpha_00;
	real<lower=0> alpha_10;
	real<lower=0> alpha_01;
	real<lower=0> alpha_11;
	real sex0_loc0;
	real sex1_loc0;
	real sex0_loc1;
	real sex1_loc1;  	
  	vector[N_predict] f_tilde_00;
	vector[N_predict] f_tilde_10;
	vector[N_predict] f_tilde_01;
	vector[N_predict] f_tilde_11;
}

transformed parameters {
  	matrix[N_predict, N_predict] L_cov;
  	vector[N_predict] logit_p_predict_00;
	vector[N_predict] logit_p_predict_10;
	vector[N_predict] logit_p_predict_01;
	vector[N_predict] logit_p_predict_11;
	// GP for 00 and 01 (women)
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_00, rho_00) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_00 = sex0_loc0 + L_cov * f_tilde_00;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_01, rho_01) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_01 = sex0_loc1 + L_cov * f_tilde_01;
	// GP for 10 and 10 (men)
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_10, rho_10) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_10 = sex1_loc0 + L_cov * f_tilde_10;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_11, rho_11) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_11 = sex1_loc1 + L_cov * f_tilde_11;
}

model {
	rho_00 ~ normal(0, rho_hyper_par_00);
	rho_10 ~ normal(0, rho_hyper_par_10);
	rho_01 ~ normal(0, rho_hyper_par_01);
	rho_11 ~ normal(0, rho_hyper_par_11);  	
  	alpha_00 ~ normal(0, alpha_hyper_par_00);
	alpha_10 ~ normal(0, alpha_hyper_par_10);
	alpha_01 ~ normal(0, alpha_hyper_par_01);
	alpha_11 ~ normal(0, alpha_hyper_par_11);
  	sex0_loc0 ~ normal( 0 , 10 );
	sex0_loc1 ~ normal( 0 , 10 );
	sex1_loc0 ~ normal( 0 , 10 );
	sex1_loc1 ~ normal( 0 , 10 );
  	f_tilde_00 ~ normal(0, 1);
	f_tilde_01 ~ normal(0, 1);
	f_tilde_10 ~ normal(0, 1);
	f_tilde_11 ~ normal(0, 1);
  	y_observed_00 ~ binomial_logit(total_observed_00, logit_p_predict_00[observed_idx] );
	y_observed_01 ~ binomial_logit(total_observed_01, logit_p_predict_01[observed_idx] );
	y_observed_10 ~ binomial_logit(total_observed_10, logit_p_predict_10[observed_idx] );
	y_observed_11 ~ binomial_logit(total_observed_11, logit_p_predict_11[observed_idx] );
}

generated quantities {
  	vector[N_predict] p_predict_00;
	vector[N_predict] p_predict_01;
	vector[N_predict] p_predict_10;
	vector[N_predict] p_predict_11;
	p_predict_00 = inv_logit(logit_p_predict_00);
	p_predict_01 = inv_logit(logit_p_predict_01);  
	p_predict_10 = inv_logit(logit_p_predict_10);
	p_predict_11 = inv_logit(logit_p_predict_11);
}
"

	stan.model <- stan_model(model_name= 'gp_all',model_code = gsub('\t',' ',stan.code))	
	stan.data <- list()
	stan.data$x_predict <- seq(vla[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
	stan.data$N_predict <- length(stan.data$x_predict)
	stan.data$observed_idx <- which(stan.data$x_predict%%1==0.5)
	stan.data$N_observed <- length(stan.data$observed_idx)
	stan.data$y_observed_00 <- vla[SEX==0 & LOC==0, VLNS_N]
	stan.data$y_observed_10 <- vla[SEX==1 & LOC==0, VLNS_N]
	stan.data$y_observed_01 <- vla[SEX==0 & LOC==1, VLNS_N]
	stan.data$y_observed_11 <- vla[SEX==1 & LOC==1, VLNS_N]
	stan.data$total_observed_00 <- vla[SEX==0 & LOC==0, N]
	stan.data$total_observed_10 <- vla[SEX==1 & LOC==0, N]
	stan.data$total_observed_01 <- vla[SEX==0 & LOC==1, N]
	stan.data$total_observed_11 <- vla[SEX==1 & LOC==1, N]
	stan.data$rho_hyper_par_00 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_10 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_01 <- diff(range(stan.data$x_predict))/3
	stan.data$rho_hyper_par_11 <- diff(range(stan.data$x_predict))/3
	stan.data$alpha_hyper_par_00 <- 2
	stan.data$alpha_hyper_par_10 <- 2
	stan.data$alpha_hyper_par_01 <- 2
	stan.data$alpha_hyper_par_11 <- 2
	fit <- sampling(stan.model, data=stan.data, iter=10e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
	save(fit, file=file.path(outdir, "200428f_suppAmongPop_gp_stanfit.rda"))
	
	min( summary(fit)$summary[, 'n_eff'] )	
	re <- rstan::extract(fit)		
	
	#
	#	extract hyperparams rho
	ps <- c(0.025,0.25,0.5,0.75,0.975)
	tmp <- cbind( quantile(re$rho_00, probs=ps),
			quantile(re$rho_10, probs=ps),
			quantile(re$rho_01, probs=ps),
			quantile(re$rho_11, probs=ps),
			quantile(re$alpha_00, probs=ps),
			quantile(re$alpha_10, probs=ps),
			quantile(re$alpha_01, probs=ps),
			quantile(re$alpha_11, probs=ps) )			
	colnames(tmp) <- c('rho_00','rho_10','rho_01','rho_11','alpha_00','alpha_10','alpha_01','alpha_11')
	rownames(tmp) <- c('CL','IL','M','IU','CU')
	tmp <- as.data.table(reshape2::melt(tmp))
	setnames(tmp, 'Var2', 'GP_hyper_par')
	tmp[, SEX:= as.integer(gsub('^([a-z]+)_([0-9])([0-9])','\\2',GP_hyper_par))]
	tmp[, LOC:= as.integer(gsub('^([a-z]+)_([0-9])([0-9])','\\3',GP_hyper_par))]
	tmp[, GP_hyper_par:= gsub('^([a-z]+)_([0-9])([0-9])','\\1',GP_hyper_par)]
	tmp <- dcast.data.table(tmp, LOC+SEX+GP_hyper_par~Var1, value.var='value')
	nspop.gp.pars <- merge(unique(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), tmp, by=c('SEX','LOC'))
	ggplot(nspop.gp.pars, aes(x= paste0(GP_hyper_par, ' ', LOC_LABEL, ' ', SEX_LABEL))) +
			geom_point(aes(y=M)) +
			geom_errorbar(aes(ymin=CL, ymax=CU)) +
			coord_flip() +
			theme_bw() +
			labs(x='GP hyperparameter\n', y='')
	ggsave(file=file.path(prjdir,'results_200220','200428f_notsuppAmongPop_gppars.pdf'), w=6, h=3)
	
	
	#
	#	make prevalence plot by age
	ps <- c(0.025,0.5,0.975)
	tmp <- cbind( apply(re$p_predict_00, 2, quantile, probs=ps),
			apply(re$p_predict_10, 2, quantile, probs=ps),
			apply(re$p_predict_01, 2, quantile, probs=ps),
			apply(re$p_predict_11, 2, quantile, probs=ps)
	)
	rownames(tmp) <- c('CL','M','CU')
	tmp <- as.data.table(reshape2::melt(tmp))
	nspop.by.age <- dcast.data.table(tmp, Var2~Var1, value.var='value')
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	nspop.by.age <- cbind(tmp, nspop.by.age) 
	nspop.by.age <- merge(unique(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL))), nspop.by.age, by=c('SEX','LOC'))
	nspop.by.age[, STAT:='VLNS']
	
	ggplot(nspop.by.age) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +			
			geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			facet_wrap(~LOC_LABEL, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='population with unsuppressed viral load\n(95% credibility interval)\n', 
					colour='gender', 
					linetype='location')
	ggsave(file=file.path(prjdir,'results_200220','200428f_nsuppAmongPop_vs_age_by_gender_fishinland_stan.pdf'), w=6, h=5)
	
	
	#	extract basic not supp estimates
	vla[, list(N=sum(VLNS_N), P=sum(VLNS_N) / sum(N), N2=sum(N)-sum(VLNS_N), P2= 1-sum(VLNS_N) / sum(N)), by=c('LOC_LABEL','SEX_LABEL')]		
	#   LOC_LABEL SEX_LABEL   N          P   N2        P2
	#1:   fishing         M 296 0.14041746 1812 0.8595825
	#2:    inland         M 203 0.03098764 6348 0.9690124
	#3:   fishing         F 198 0.10216718 1740 0.8978328
	#4:    inland         F 279 0.03467562 7767 0.9653244
	
	
	ps <- c(0.025,0.5,0.975)
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt( tmp ))
	setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, ROW_ID:= seq_len(nrow(tmp))]
	tmp <- merge(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
	rp <- merge(tmp, rp, by=c('ROW_ID'))	
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]
	rp <- rp[, list(Q= quantile(P, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','SEX_LABEL')]
	rp <- dcast.data.table(rp, LOC_LABEL+SEX_LABEL~P, value.var='Q')
	rp[, LABEL:= paste0(round(M*100, d=2),'% (',round(CL*100, d=2),'% - ',round(CU*100,d=2),'%)') ]
	nspop.by.sex.loc <- copy(rp)
	# 1:   fishing         M 0.12509928 0.15440179 0.13929310 13.93% (12.51% - 15.44%)
	#2:   fishing         F 0.08927186 0.11554352 0.10202937   10.2% (8.93% - 11.55%)
	#3:    inland         M 0.02666003 0.03491415 0.03068550    3.07% (2.67% - 3.49%)
	#4:    inland         F 0.03060328 0.03846557 0.03442258    3.44% (3.06% - 3.85%)


	#	extract risk ratio of suppressed VL female:male and male:female
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt( tmp ))
	setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, ROW_ID:= seq_len(nrow(tmp))]
	tmp <- merge(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
	rp <- merge(tmp, rp, by=c('ROW_ID'))	
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]
	rp <- dcast.data.table(rp, LOC_LABEL+iterations~SEX_LABEL, value.var='P')
	rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
	rp <- melt(rp, id.vars=c('LOC_LABEL','iterations'))
	rp <- rp[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','variable')]
	rp <- dcast.data.table(rp, LOC_LABEL+variable~P, value.var='Q')
	rp[, LABEL:= paste0(sprintf('%2.2f',M),' (',sprintf('%2.2f',CL),' - ',sprintf('%2.2f',CU),')') ]	
	nspop.by.loc <- copy(rp)	
	#   LOC_LABEL variable        CL        CU         M              LABEL
	#1:   fishing    PR_FM 0.6193378 0.8632857 0.7321284 0.73 (0.62 - 0.86)
	#2:   fishing    PR_MF 1.1583654 1.6146279 1.3658807 1.37 (1.16 - 1.61)
	#3:    inland    PR_FM 0.9408229 1.3432295 1.1221810 1.12 (0.94 - 1.34)
	#4:    inland    PR_MF 0.7444744 1.0628992 0.8911218 0.89 (0.74 - 1.06)

	
	#	extract risk ratio of unsuppressed VL female:male and male:female by age
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt( tmp ))
	setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, ROW_ID:= seq_len(nrow(tmp))]
	tmp <- merge(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
	rp <- merge(tmp, rp, by=c('ROW_ID'))	
	rp <- dcast.data.table(rp, LOC_LABEL+iterations+AGE_LABEL~SEX_LABEL, value.var='P')
	rp[, PR_FM:= F/M]
	rp[, PR_MF:=M/F]
	rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))
	rp <- rp[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]	
	rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
	rp[, LABEL:= paste0(sprintf('%2.2f',M),' (',sprintf('%2.2f',CL),' - ',sprintf('%2.2f',CU),')') ]	
	nspop.ratio.by.loc.age <- copy(rp)


	#	extract if difference in female:male risk ratio of unsuppressed VL is different in fishing vs inland
	tmp <- cbind( re$p_predict_00, re$p_predict_10, re$p_predict_01, re$p_predict_11)
	rp <- as.data.table(reshape2::melt( tmp ))
	setnames(rp, 1:3, c('iterations','ROW_ID','P')) 
	tmp <- as.data.table(expand.grid(AGE_LABEL=stan.data$x_predict, SEX=c(0,1), LOC=c(0,1)))
	tmp[, ROW_ID:= seq_len(nrow(tmp))]
	tmp <- merge(subset(vla, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL,N)), tmp, by=c('SEX','LOC','AGE_LABEL'), all=TRUE)
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N]) 
	set(tmp, tmp[, which(is.na(N))], 'N', tmp[which(is.na(N))-1L, N])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(SEX_LABEL))], 'SEX_LABEL', tmp[which(is.na(SEX_LABEL))-1L, SEX_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])
	set(tmp, tmp[, which(is.na(LOC_LABEL))], 'LOC_LABEL', tmp[which(is.na(LOC_LABEL))-1L, LOC_LABEL])	
	rp <- merge(tmp, rp, by=c('ROW_ID'))	
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC_LABEL','SEX_LABEL','iterations')]
	rp <- dcast.data.table(rp, LOC_LABEL+iterations~SEX_LABEL, value.var='P')
	rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
	rp <- dcast.data.table(rp, iterations~LOC_LABEL, value.var='PR_FM')
	rp[, PR_FM_D:= fishing-inland]	
	rp[, list(Q= quantile(PR_FM_D, probs=ps), P=c('CL','M','CU'))]
	#	            Q  P
	#1: -0.6360743 CL
	#2: -0.3903490  M
	#3: -0.1646397 CU


	save(vla, re, nspop.by.age, nspop.by.sex.loc, nspop.ratio.by.loc.age, file=file.path(outdir, "200428f_suppAmongPop.rda"))


	#	make table version suppressed
	nspop.by.age[, LABEL:= paste0(sprintf('%2.1f',M*100),' (',sprintf('%2.1f',CL*100),' - ',sprintf('%2.1f',CU*100),')') ]
	set(nspop.by.age, NULL, 'SEX_LABEL', nspop.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])	
	setnames(nspop.ratio.by.loc.age,'LABEL', 'LABEL2')
	dt <- subset(nspop.by.age, AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5))	
	dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
	tmp <- subset(nspop.ratio.by.loc.age, variable=="PR_FM" & AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5), c(LOC_LABEL, AGE_LABEL, LABEL2))
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	tmp <- subset(nspop.ratio.by.loc.age, variable=="PR_MF" & AGE_LABEL%in%c(20.5,25.5,30.5,35.5,40.5,45.5), c(LOC_LABEL, AGE_LABEL, LABEL2))
	setnames(tmp, 'LABEL2', 'PR_MF')
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	write.csv(dt, row.names=FALSE, file=file.path(outdir, "200428f_suppAmongPop.csv"))

}

vl.suppofpop.by.gender.loc.age.icar<- function()
{
	require(Hmisc)
	require(data.table)
	require(ggplot2)
	require(rstan)
	vl.detectable <- 4e2
	vl.suppressed <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	outdir <- file.path(prjdir,'results_200220')
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
	
	subset(ds, FC=='inland' & SEX=='F' & AGEYRS==25)
	
	tmp <- seq.int(min(ds$AGEYRS), max(ds$AGEYRS))
	vla <- as.data.table(expand.grid(FC=c('fishing','inland'), SEX=c('M','F'), AGEYRS=tmp))
	vla <- vla[, {		
				z <- which(ds$FC==FC & ds$SEX==SEX & ds$AGEYRS==AGEYRS)	
				list(	N= length(z),
						HIV_N= length(which(ds$HIV_STATUS[z]==1)),
						VLNS_N= length(which(ds$VLNS[z]==1)),
						ARV_N= length(which(ds$ARVMED[z]==0 & ds$HIV_STATUS[z]==1 & !is.na(ds$ARVMED[z])))
				)				
			}, by=c('FC','SEX','AGEYRS')]
	setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
	vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
	vla[, SEX:= as.integer(SEX_LABEL=='M')]
	vla[, AGE:= AGE_LABEL-14L]
	vla[, ROW_ID:= seq_len(nrow(vla))]
	
	
	stan.code2 <- "functions{			
			real icar_normal_lpdf( vector phi, int N, int[] node1, int[] node2)
			{			
			return -0.5 * dot_self(phi[node1] - phi[node2]) + normal_lpdf( sum(phi) | 0, 0.001 * N);
			}
			}
			
			data{
			int<lower=1> N; 
			int<lower=0> TOTAL[N];    
			int<lower=0> K[N];
			int<lower=0> AGE_N;  
			int<lower=0, upper=AGE_N> AGE[N];  
			int<lower=0,upper=1> SEX[N];
			int<lower=0,upper=1> LOC[N];	 
			int<lower=0> N_edges; 						// number of related age groups
			int<lower=1, upper=AGE_N> node1[N_edges];  	// node1[i] adjacent to node2[i]
			int<lower=1, upper=AGE_N> node2[N_edges];  	// and node1[i] < node2[i]
			}
			
			parameters{
			real sex0_loc0;
			real sex1_loc0;
			real sex0_loc1;
			real sex1_loc1;
			vector[AGE_N] phi_sex0_loc0;
			vector[AGE_N] phi_sex1_loc0;
			vector[AGE_N] phi_sex0_loc1;
			vector[AGE_N] phi_sex1_loc1; 			 		
			real<lower=0> sigma_loc0;
			real<lower=0> sigma_loc1;
			}
			
			transformed parameters{
			vector[N] p_logit;
			for ( i in 1:N ) 
			{
			p_logit[i] = sex0_loc0 * (1-SEX[i]) * (1-LOC[i]) + 
			sex1_loc0 * SEX[i] * (1-LOC[i]) + 
			sex0_loc1 * (1-SEX[i]) * LOC[i] + 
			sex1_loc1 * SEX[i] * LOC[i] +
			sigma_loc0 * phi_sex0_loc0[AGE[i]] * (1-SEX[i]) * (1-LOC[i]) +
			sigma_loc0 * phi_sex1_loc0[AGE[i]] * SEX[i] * (1-LOC[i]) +
			sigma_loc1 * phi_sex0_loc1[AGE[i]] * (1-SEX[i]) * LOC[i] +
			sigma_loc1 * phi_sex1_loc1[AGE[i]] * SEX[i] * LOC[i];
			}
			}
			
			model{	
			sex0_loc0 ~ normal( 0 , 10 );
			sex1_loc0 ~ normal( 0 , 10 );
			sex0_loc1 ~ normal( 0 , 10 );
			sex1_loc1 ~ normal( 0 , 10 );
			phi_sex0_loc0 ~ icar_normal_lpdf(AGE_N , node1 , node2);
			phi_sex1_loc0 ~ icar_normal_lpdf(AGE_N , node1 , node2);
			phi_sex0_loc1 ~ icar_normal_lpdf(AGE_N , node1 , node2);
			phi_sex1_loc1 ~ icar_normal_lpdf(AGE_N , node1 , node2);
			sigma_loc0 ~ normal(0.0, 1);
			sigma_loc1 ~ normal(0.0, 1);
			K ~ binomial_logit( TOTAL , p_logit );
			}
			
			generated quantities{
			vector[N] p;
			p= inv_logit(p_logit);
			}"
	
	stan.model2 <- stan_model(model_name= 'icar_age_interactions',model_code = gsub('\t',' ',stan.code2))
	stan.data <- list()
	stan.data$N <- nrow(vla)
	stan.data$TOTAL <- vla[,N]
	stan.data$K <- vla[,VLNS_N]
	stan.data$AGE_N <- vla[, max(AGE)]
	stan.data$AGE <- vla[, AGE]
	stan.data$SEX <- vla[, SEX]
	stan.data$LOC <- vla[, LOC]
	#	second order RW prior
	stan.data$node1 <-  c(vla[, seq.int(1, max(AGE)-1L)], vla[, seq.int(1, max(AGE)-2L)])
	stan.data$node2 <-  c(vla[, seq.int(2, max(AGE))], vla[, seq.int(3, max(AGE))])
	tmp <- sort(stan.data$node1, index.return=TRUE)$ix
	stan.data$node1 <- stan.data$node1[tmp]
	stan.data$node2 <- stan.data$node2[tmp]
	stan.data$N_edges <-  length(stan.data$node1)
	fit <- sampling(stan.model2, data=stan.data, iter=20e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999) )	
	save(fit, file=file.path(outdir, "notsuppAmongPop_icar_stanfit_200428.rda"))
	
	
	min( summary(fit)$summary[, 'n_eff'] )	
	re <- rstan::extract(fit)	
	ps <- c(0.025,0.5,0.975)
	
	
	#	make prevalence plot by age
	tmp <- apply(re$p, 2, quantile, probs=ps)
	rownames(tmp) <- c('CL','M','CU')
	tmp <- as.data.table(reshape2::melt(tmp))	
	nspop.by.age <- cbind(vla, dcast.data.table(tmp, Var2~Var1, value.var='value'))	
	ggplot(nspop.by.age) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +			
			geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			facet_wrap(~LOC_LABEL, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='individuals with unsuppressed viral load\n(95% credibility interval)\n', 
					colour='gender')
	ggsave(file=file.path(prjdir,'results_200220','200428_notsuppAmongPop_vs_age_by_gender_fishinland_stan.pdf'), w=6, h=5)	
	
	
	#	extract basic not supp estimates
	vla[, list(N=sum(VLNS_N), P=sum(VLNS_N) / sum(N), N2=sum(N)-sum(VLNS_N), P2= 1-sum(VLNS_N) / sum(N)), by=c('LOC_LABEL','SEX_LABEL')]		
	#   LOC_LABEL SEX_LABEL   N          P   N2        P2
	#1:   fishing         M 296 0.14041746 1812 0.8595825
	#2:    inland         M 203 0.03098764 6348 0.9690124
	#3:   fishing         F 198 0.10216718 1740 0.8978328
	#4:    inland         F 279 0.03467562 7767 0.9653244
	
	
	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
	rp <- rp[, list(Q= quantile(P, probs=ps), P=c('CL','M','CU')), by=c('LOC','SEX')]
	rp <- dcast.data.table(rp, LOC+SEX~P, value.var='Q')
	rp[, LABEL:= paste0(round(M, d=2),'% (',round(CL, d=2),'% - ',round(CU,d=2),'%)') ]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
	nspop.by.sex.loc <- copy(rp)
	#	   LOC SEX LOC_LABEL SEX_LABEL         CL         CU         M                 LABEL
	#1:   0   0    inland         F 0.03080205 0.03877159 0.03465293 0.03% (0.03% - 0.04%)
	#2:   0   1    inland         M 0.02696284 0.03527088 0.03093918 0.03% (0.03% - 0.04%)
	#3:   1   0   fishing         F 0.08900840 0.11607434 0.10193459  0.1% (0.09% - 0.12%)
	#4:   1   1   fishing         M 0.12619340 0.15521611 0.14040214 0.14% (0.13% - 0.16%)
	
	
	#	extract risk ratio of unsuppressed VL female:male and male:female
	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
	rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
	rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC','iterations')]
	rp <- melt(rp, id.vars=c('LOC','iterations'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC','variable')]
	rp <- dcast.data.table(rp, LOC+variable~P, value.var='Q')
	rp[, LABEL:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL))), rp, by=c('LOC'))
	nspop.ratio.by.loc <- copy(rp)
	#	   LOC LOC_LABEL variable        CL        CU         M              LABEL
	#1:   0    inland    PR_FM 0.9375862 1.3330163 1.1200999 1.12 (0.94 - 1.33)
	#2:   0    inland    PR_MF 0.7501784 1.0665686 0.8927775 0.89 (0.75 - 1.07)
	#3:   1   fishing    PR_FM 0.6131129 0.8584486 0.7261857 0.73 (0.61 - 0.86)
	#4:   1   fishing    PR_MF 1.1648921 1.6310210 1.3770582 1.38 (1.16 - 1.63)
	
	
	#	extract if difference in female:male risk ratio of unsuppressed VL is different in fishing vs inland
	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
	rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
	rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
	rp <- dcast.data.table(rp, iterations~LOC_LABEL, value.var='PR_FM')
	rp[, PR_FM_D:= fishing-inland]	
	rp[, list(Q= quantile(PR_FM_D, probs=ps), P=c('CL','M','CU'))]
	#	            Q  P
	#1: -0.6336039 CL
	#2: -0.3931299  M
	#3: -0.1680829 CU
	
	
	#	extract risk ratio of unsuppressed VL female:male and male:female by age
	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations+AGE+AGE_LABEL~SEX_LABEL, value.var='P')
	rp[, PR_FM:= F/M]
	rp[, PR_MF:=M/F]
	rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]
	rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
	nspop.ratio.by.loc.age <- copy(rp)
	
	
	save(vla, re, nspop.by.age, nspop.ratio.by.loc, nspop.ratio.by.loc.age, file=file.path(outdir, "notsuppAmongPop_200428.rda"))
	
	#	make table
	nspop.by.age[, LABEL:= paste0(round(M*100, d=1),' (',round(CL*100, d=1),' - ',round(CU*100,d=1),')') ]
	set(nspop.by.age, NULL, 'SEX_LABEL', nspop.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])	
	nspop.ratio.by.loc.age[, LABEL2:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
	dt <- subset(nspop.by.age, AGE_LABEL%in%c(20,25,30,35,40,45))	
	dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
	tmp <- subset(nspop.ratio.by.loc.age, variable=="PR_FM" & AGE_LABEL%in%c(20,25,30,35,40,45), c(LOC_LABEL, AGE_LABEL, LABEL2))
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	tmp <- subset(nspop.ratio.by.loc.age, variable=="PR_MF" & AGE_LABEL%in%c(20,25,30,35,40,45), c(LOC_LABEL, AGE_LABEL, LABEL2))
	setnames(tmp, 'LABEL2', 'PR_MF')
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	write.csv(dt, row.names=FALSE, file=file.path(outdir, "notsuppAmongPop_200428.csv"))
	
}


vl.suppofpop.by.gender.loc.age.icar<- function()
{
	require(Hmisc)
	require(data.table)
	require(ggplot2)
	require(rstan)
	vl.detectable <- 4e2
	vl.suppressed <- 1e3
	
	
	require(data.table)
	prjdir	<- '~/Box/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- file.path(prjdir,'data','191101_data_round17_vl_gps.rda')
	outdir <- file.path(prjdir,'results_200220')
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
	
	subset(ds, FC=='inland' & SEX=='F' & AGEYRS==25)
	
	tmp <- seq.int(min(ds$AGEYRS), max(ds$AGEYRS))
	vla <- as.data.table(expand.grid(FC=c('fishing','inland'), SEX=c('M','F'), AGEYRS=tmp))
	vla <- vla[, {		
				z <- which(ds$FC==FC & ds$SEX==SEX & ds$AGEYRS==AGEYRS)	
				list(	N= length(z),
						HIV_N= length(which(ds$HIV_STATUS[z]==1)),
						VLNS_N= length(which(ds$VLNS[z]==1)),
						ARV_N= length(which(ds$ARVMED[z]==0 & ds$HIV_STATUS[z]==1 & !is.na(ds$ARVMED[z])))
				)				
			}, by=c('FC','SEX','AGEYRS')]
	setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
	vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
	vla[, SEX:= as.integer(SEX_LABEL=='M')]
	vla[, AGE:= AGE_LABEL-14L]
	vla[, ROW_ID:= seq_len(nrow(vla))]
	
	
	stan.code2 <- "functions{			
	real icar_normal_lpdf( vector phi, int N, int[] node1, int[] node2)
	{			
		return -0.5 * dot_self(phi[node1] - phi[node2]) + normal_lpdf( sum(phi) | 0, 0.001 * N);
	}
}
			
data{
	int<lower=1> N; 
	int<lower=0> TOTAL[N];    
	int<lower=0> K[N];
	int<lower=0> AGE_N;  
	int<lower=0, upper=AGE_N> AGE[N];  
	int<lower=0,upper=1> SEX[N];
	int<lower=0,upper=1> LOC[N];	 
	int<lower=0> N_edges; 						// number of related age groups
	int<lower=1, upper=AGE_N> node1[N_edges];  	// node1[i] adjacent to node2[i]
	int<lower=1, upper=AGE_N> node2[N_edges];  	// and node1[i] < node2[i]
}
			
parameters{
	real sex0_loc0;
	real sex1_loc0;
	real sex0_loc1;
	real sex1_loc1;
	vector[AGE_N] phi_sex0_loc0;
	vector[AGE_N] phi_sex1_loc0;
	vector[AGE_N] phi_sex0_loc1;
	vector[AGE_N] phi_sex1_loc1; 			 		
	real<lower=0> sigma_loc0;
	real<lower=0> sigma_loc1;
}
			
transformed parameters{
	vector[N] p_logit;
	for ( i in 1:N ) 
	{
		p_logit[i] = sex0_loc0 * (1-SEX[i]) * (1-LOC[i]) + 
		sex1_loc0 * SEX[i] * (1-LOC[i]) + 
		sex0_loc1 * (1-SEX[i]) * LOC[i] + 
		sex1_loc1 * SEX[i] * LOC[i] +
		sigma_loc0 * phi_sex0_loc0[AGE[i]] * (1-SEX[i]) * (1-LOC[i]) +
		sigma_loc0 * phi_sex1_loc0[AGE[i]] * SEX[i] * (1-LOC[i]) +
		sigma_loc1 * phi_sex0_loc1[AGE[i]] * (1-SEX[i]) * LOC[i] +
		sigma_loc1 * phi_sex1_loc1[AGE[i]] * SEX[i] * LOC[i];
	}
}
			
model{	
	sex0_loc0 ~ normal( 0 , 10 );
	sex1_loc0 ~ normal( 0 , 10 );
	sex0_loc1 ~ normal( 0 , 10 );
	sex1_loc1 ~ normal( 0 , 10 );
	phi_sex0_loc0 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	phi_sex1_loc0 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	phi_sex0_loc1 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	phi_sex1_loc1 ~ icar_normal_lpdf(AGE_N , node1 , node2);
	sigma_loc0 ~ normal(0.0, 1);
	sigma_loc1 ~ normal(0.0, 1);
	K ~ binomial_logit( TOTAL , p_logit );
}
			
generated quantities{
	vector[N] p;
	p= inv_logit(p_logit);
}"
	
	stan.model2 <- stan_model(model_name= 'icar_age_interactions',model_code = gsub('\t',' ',stan.code2))
	stan.data <- list()
	stan.data$N <- nrow(vla)
	stan.data$TOTAL <- vla[,N]
	stan.data$K <- vla[,VLNS_N]
	stan.data$AGE_N <- vla[, max(AGE)]
	stan.data$AGE <- vla[, AGE]
	stan.data$SEX <- vla[, SEX]
	stan.data$LOC <- vla[, LOC]
	#	second order RW prior
	stan.data$node1 <-  c(vla[, seq.int(1, max(AGE)-1L)], vla[, seq.int(1, max(AGE)-2L)])
	stan.data$node2 <-  c(vla[, seq.int(2, max(AGE))], vla[, seq.int(3, max(AGE))])
	tmp <- sort(stan.data$node1, index.return=TRUE)$ix
	stan.data$node1 <- stan.data$node1[tmp]
	stan.data$node2 <- stan.data$node2[tmp]
	stan.data$N_edges <-  length(stan.data$node1)
	fit <- sampling(stan.model2, data=stan.data, iter=20e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999) )	
	save(fit, file=file.path(outdir, "notsuppAmongPop_icar_stanfit_200428.rda"))
	
		
	min( summary(fit)$summary[, 'n_eff'] )	
	re <- rstan::extract(fit)	
	ps <- c(0.025,0.5,0.975)
	
	
	#	make prevalence plot by age
	tmp <- apply(re$p, 2, quantile, probs=ps)
	rownames(tmp) <- c('CL','M','CU')
	tmp <- as.data.table(reshape2::melt(tmp))	
	nspop.by.age <- cbind(vla, dcast.data.table(tmp, Var2~Var1, value.var='value'))	
	ggplot(nspop.by.age) + 		
			geom_ribbon(aes(x=AGE_LABEL, ymin=CL, ymax=CU, group=interaction(SEX_LABEL,LOC_LABEL)), alpha=0.2) +			
			geom_line(aes(x=AGE_LABEL, y=M, colour=SEX_LABEL)) +
			scale_x_continuous( expand=c(0,0) ) + 
			scale_y_continuous(label=scales:::percent) +
			scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
			facet_wrap(~LOC_LABEL, ncol=2) +
			theme_bw() +
			labs(x='\nage at visit (years)', 
					y='individuals with unsuppressed viral load\n(95% credibility interval)\n', 
					colour='gender')
	ggsave(file=file.path(prjdir,'results_200220','200428_notsuppAmongPop_vs_age_by_gender_fishinland_stan.pdf'), w=6, h=5)	
	
	
	#	extract basic not supp estimates
	vla[, list(N=sum(VLNS_N), P=sum(VLNS_N) / sum(N), N2=sum(N)-sum(VLNS_N), P2= 1-sum(VLNS_N) / sum(N)), by=c('LOC_LABEL','SEX_LABEL')]		
	#   LOC_LABEL SEX_LABEL   N          P   N2        P2
	#1:   fishing         M 296 0.14041746 1812 0.8595825
	#2:    inland         M 203 0.03098764 6348 0.9690124
	#3:   fishing         F 198 0.10216718 1740 0.8978328
	#4:    inland         F 279 0.03467562 7767 0.9653244
	

	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
	rp <- rp[, list(Q= quantile(P, probs=ps), P=c('CL','M','CU')), by=c('LOC','SEX')]
	rp <- dcast.data.table(rp, LOC+SEX~P, value.var='Q')
	rp[, LABEL:= paste0(round(M, d=2),'% (',round(CL, d=2),'% - ',round(CU,d=2),'%)') ]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
	nspop.by.sex.loc <- copy(rp)
	#	   LOC SEX LOC_LABEL SEX_LABEL         CL         CU         M                 LABEL
	#1:   0   0    inland         F 0.03080205 0.03877159 0.03465293 0.03% (0.03% - 0.04%)
	#2:   0   1    inland         M 0.02696284 0.03527088 0.03093918 0.03% (0.03% - 0.04%)
	#3:   1   0   fishing         F 0.08900840 0.11607434 0.10193459  0.1% (0.09% - 0.12%)
	#4:   1   1   fishing         M 0.12619340 0.15521611 0.14040214 0.14% (0.13% - 0.16%)
	

	#	extract risk ratio of unsuppressed VL female:male and male:female
	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
	rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
	rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC','iterations')]
	rp <- melt(rp, id.vars=c('LOC','iterations'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC','variable')]
	rp <- dcast.data.table(rp, LOC+variable~P, value.var='Q')
	rp[, LABEL:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL))), rp, by=c('LOC'))
	nspop.ratio.by.loc <- copy(rp)
	#	   LOC LOC_LABEL variable        CL        CU         M              LABEL
	#1:   0    inland    PR_FM 0.9375862 1.3330163 1.1200999 1.12 (0.94 - 1.33)
	#2:   0    inland    PR_MF 0.7501784 1.0665686 0.8927775 0.89 (0.75 - 1.07)
	#3:   1   fishing    PR_FM 0.6131129 0.8584486 0.7261857 0.73 (0.61 - 0.86)
	#4:   1   fishing    PR_MF 1.1648921 1.6310210 1.3770582 1.38 (1.16 - 1.63)
	
	
	#	extract if difference in female:male risk ratio of unsuppressed VL is different in fishing vs inland
	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- rp[, list(P= sum(P*N)/sum(N) ), by=c('LOC','SEX','iterations')]
	rp <- merge(unique(subset(vla, select=c(LOC,LOC_LABEL,SEX,SEX_LABEL))), rp, by=c('LOC','SEX'))
	rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations~SEX_LABEL, value.var='P')
	rp <- rp[, list(PR_FM= F/M, PR_MF=M/F), by=c('LOC_LABEL','iterations')]
	rp <- dcast.data.table(rp, iterations~LOC_LABEL, value.var='PR_FM')
	rp[, PR_FM_D:= fishing-inland]	
	rp[, list(Q= quantile(PR_FM_D, probs=ps), P=c('CL','M','CU'))]
	#	            Q  P
	#1: -0.6336039 CL
	#2: -0.3931299  M
	#3: -0.1680829 CU
	

	#	extract risk ratio of unsuppressed VL female:male and male:female by age
	rp <- as.data.table(reshape2::melt( re$p ))
	setnames(rp, 2:3, c('ROW_ID','P')) 
	rp <- merge(rp, vla, by='ROW_ID')
	rp <- dcast.data.table(rp, LOC+LOC_LABEL+iterations+AGE+AGE_LABEL~SEX_LABEL, value.var='P')
	rp[, PR_FM:= F/M]
	rp[, PR_MF:=M/F]
	rp <- melt(rp, id.vars=c('LOC_LABEL','AGE_LABEL','iterations'), measure.vars=c('PR_FM','PR_MF'))[, list(Q= quantile(value, probs=ps), P=c('CL','M','CU')), by=c('LOC_LABEL','AGE_LABEL','variable')]
	rp <- dcast.data.table(rp, LOC_LABEL+AGE_LABEL+variable~P, value.var='Q')
	nspop.ratio.by.loc.age <- copy(rp)

	
	save(vla, re, nspop.by.age, nspop.ratio.by.loc, nspop.ratio.by.loc.age, file=file.path(outdir, "notsuppAmongPop_200428.rda"))
	
	#	make table
	nspop.by.age[, LABEL:= paste0(round(M*100, d=1),' (',round(CL*100, d=1),' - ',round(CU*100,d=1),')') ]
	set(nspop.by.age, NULL, 'SEX_LABEL', nspop.by.age[, factor(as.character(SEX_LABEL), levels=c('F','M'))])	
	nspop.ratio.by.loc.age[, LABEL2:= paste0(round(M, d=2),' (',round(CL, d=2),' - ',round(CU,d=2),')') ]
	dt <- subset(nspop.by.age, AGE_LABEL%in%c(20,25,30,35,40,45))	
	dt <- dcast.data.table(dt, LOC_LABEL+AGE_LABEL~SEX_LABEL, value.var='LABEL')
	tmp <- subset(nspop.ratio.by.loc.age, variable=="PR_FM" & AGE_LABEL%in%c(20,25,30,35,40,45), c(LOC_LABEL, AGE_LABEL, LABEL2))
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	tmp <- subset(nspop.ratio.by.loc.age, variable=="PR_MF" & AGE_LABEL%in%c(20,25,30,35,40,45), c(LOC_LABEL, AGE_LABEL, LABEL2))
	setnames(tmp, 'LABEL2', 'PR_MF')
	dt <- merge(dt, tmp, by=c('LOC_LABEL','AGE_LABEL'))
	write.csv(dt, row.names=FALSE, file=file.path(outdir, "notsuppAmongPop_200428.csv"))
	
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
	
	subset(ds, FC=='inland' & SEX=='F' & AGEYRS==25)
	
	tmp <- seq.int(min(ds$AGEYRS), max(ds$AGEYRS))
	vla <- as.data.table(expand.grid(FC=c('fishing','inland'), SEX=c('M','F'), AGEYRS=tmp))
	ans <- vla[, {		
				z <- which(ds$FC==FC & ds$SEX==SEX & ds$AGEYRS<=(AGEYRS+2) & ds$AGEYRS>=(AGEYRS-2))				
				z2<- as.vector( unname( binconf( length(which(ds$HIV_STATUS[z]==1)), length(z) ) ) )
				z3<- as.vector( unname( binconf( length(which(ds$VLNS[z]==1)), length(z) ) ) )
				z4<- as.vector( unname( binconf( length(which(ds$VLNS[z]==1)), length(which(ds$HIV_STATUS[z]==1)) ) ) )
				z5<- as.vector( unname( binconf( length(which(ds$ARVMED[z]==0 & ds$HIV_STATUS[z]==1 & !is.na(ds$ARVMED[z]))), length(which(ds$HIV_STATUS[z]==1 & !is.na(ds$ARVMED[z]))) ) ) )
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
			geom_line(aes(x=AGEYRS, y=PARVofHIV_MEAN, colour=SEX), linetype='dotted') +
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


