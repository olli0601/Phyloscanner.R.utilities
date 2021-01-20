library(data.table)
library(ggplot2)
library(seqinr)
library(ggpubr)
library(tidyverse)
library(igraph)
library(dplyr)
library(rstan)
library(gridExtra)


if(1)
{
  args_dir <- list()
  args_dir[['data_dir']] <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  args_dir[['out_dir']] <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY'
  args_dir[['window_size']] <- 500
  args_dir[['batch_size']] <- 100
  args_dir[['script_dir']] <- '~/phyloscanner/data_analysis_RCCS1519'
  args_dir[['job_tag']] <- 'windowsize500_batchsize100'
  args_dir[['windown']] <- 99
  args_dir[['length_cutoff']] <- 250
}



args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-data_dir')
  stopifnot(args_line[[3]]=='-out_dir')
  stopifnot(args_line[[5]]=='-window_size')
  stopifnot(args_line[[7]]=='-batch_size')
  stopifnot(args_line[[9]]=='-script_dir')
  stopifnot(args_line[[11]]=='-job_tag')
  stopifnot(args_line[[13]]=='-windown')
  stopifnot(args_line[[15]]=='-length_cutoff')
  args_dir <- list()
  args_dir[['data_dir']] <- args_line[[2]]
  args_dir[['out_dir']] <- args_line[[4]]
  args_dir[['window_size']] <- as.integer(args_line[[6]])
  args_dir[['batch_size']] <- as.integer(args_line[[8]])
  args_dir[['script_dir']] <- args_line[[10]]
  args_dir[['job_tag']] <- args_line[[12]]
  args_dir[['windown']] <- as.integer(args_line[[14]])
  args_dir[['length_cutoff']] <- as.integer(args_line[[16]])
} 



# dir
setwd( paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ))
dir.create(file.path(args_dir[['script_dir']], 'plots'))

cat('-------------------make stan data-------------------------')
# load data
load(file.path(args_dir[['script_dir']],'RakaiPangeaMetaData_v2.rda'))
alignment <- read.fasta(file = file.path(args_dir[['data_dir']],"200422_PANGEA2_RCCSMRC_alignment.fasta"))
unique(do.call(c, lapply(alignment, unique)))
nsequence <- length(alignment)
npos <- unique(lengths(alignment))


# windows
windows_last_start <- ceiling(npos/100) * 100 - args_dir[['window_size']] + 1
windows_2last_end <- floor(npos/100) * 100
windows_start <- seq(1,windows_last_start,100)
windows_end <- c(seq(args_dir[['window_size']] ,windows_2last_end,100),npos)
windows <- seq_len(length(windows_start))
dw <- data.table(WINDOW=windows,
                 START=windows_start,
                 END=windows_end)


# map alignments to studyid
dinfo <- data.table(pangea_id=names(alignment))
id.dt <- data.table(read.csv(file.path(args_dir[['data_dir']],'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')))
id.dt <- subset(id.dt,select = c("pt_id","sex","pangea_id"))
id.dt[,pangea_id:=paste0('RCCS_',pangea_id)]
tmp <- data.table(read.csv(file.path(args_dir[['data_dir']],'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')))
tmp <- subset(tmp,select = c("pt_id","sex","pangea_id"))
tmp[,pangea_id:=paste0('MRCUVRI_',pangea_id)]
id.dt <- rbind(id.dt,tmp)
id.dt <- unique(id.dt)
dinfo <- merge(dinfo, id.dt, by="pangea_id", all.x=T)
tmp <- dinfo[is.na(pt_id)]
dinfo <- dinfo[!is.na(pt_id)]
tmp2 <- subset(rccsData,select=c('RCCS_studyid','Pangea.id','SEX'))
colnames(tmp2) <- c("pt_id","pangea_id","sex")
set(tmp, NULL, c("pt_id","sex"), NULL)
tmp <- merge(tmp, tmp2, by='pangea_id',all.x=T)
dinfo <- rbind(dinfo, tmp)
cat('No personal information found for ',nrow(dinfo[is.na(pt_id)]), ' sequences \n')

# process couple, map to studyid
couple <- unique(subset(data.table(coupdat), select=c('male.RCCS_studyid', 'female.RCCS_studyid')))
couple[, COUPLE:=1]
couple[,male.RCCS_studyid:=paste0('RK-',male.RCCS_studyid)]
couple[,female.RCCS_studyid:=paste0('RK-',female.RCCS_studyid)]

tmp_couple <- subset(couple,select=c('male.RCCS_studyid','female.RCCS_studyid', 'COUPLE'))
tmp <- copy(tmp_couple)
setnames(tmp_couple, c('male.RCCS_studyid','female.RCCS_studyid'),c('pt_id1','pt_id2'))
setnames(tmp, c('male.RCCS_studyid','female.RCCS_studyid'),c('pt_id2','pt_id1'))
tmp_couple[,sex1_couple := 'M']
tmp_couple[,sex2_couple := 'F']
tmp[,sex1_couple := 'F']
tmp[,sex2_couple := 'M']
tmp_couple <- rbind(tmp_couple, tmp)
unique(tmp_couple)
tmp_couple  <- tmp_couple[!is.na(pt_id1) & !is.na(pt_id2) & pt_id1!=pt_id2,]

# add similarity per window to couples
for (i in 1:args_dir[['windown']]) {
  # similarity per window
  cat('\n process window ',i , '\n')
  load(paste0('subsequence_window_',i,'_results.rda'))
  #
  distancei <- distancei[LENGTH >= args_dir[['length_cutoff']],]
  # add couple
  setnames(dinfo, colnames(dinfo), paste0(colnames(dinfo),'1'))
  distancei <- merge(distancei,dinfo, by.x=c('TAXA1'), by.y=c('pangea_id1'), all.x=T) 
  setnames(dinfo, colnames(dinfo), gsub('1','2',colnames(dinfo)))
  distancei <- merge(distancei,dinfo, by.x=c('TAXA2'), by.y=c('pangea_id2'), all.x=T) 
  setnames(dinfo, colnames(dinfo), gsub('2','',colnames(dinfo)))
  # remove no pt_id
  distancei <- distancei[!is.na(pt_id1) & !is.na(pt_id2) & pt_id1!=pt_id2,]
  distancei <- merge(distancei,tmp_couple, by=c('pt_id1','pt_id2'), all.x=T)
  distancei <- distancei[COUPLE==1,c('pt_id1','pt_id2','PERC','WINDOW')]
  if(i==1){
    ans <- distancei
  }else{
    ans <- rbind(ans, distancei)
  }
  
}
save(ans, file = file.path(args_dir[['script_dir']],'distance_couple_per_window.rda'))  

# make stan data
load(file.path(args_dir[['script_dir']],'distance_couple_per_window.rda'))
similarity <- list()
for (i in seq_len(max(ans$WINDOW))) {
  similarity[[i]] <- ans[WINDOW==i,]$PERC
}

max_data <- max(lengths(similarity))

for (i in seq_len(max(ans$WINDOW))) {
  if(length(similarity[[i]]) !=max_data){
    similarity[[i]] <- c(similarity[[i]], rep(NA,max_data-length(similarity[[i]])))
  }
}

stan_data <- list()
stan_data$similarity <- do.call(rbind,similarity)
stan_data$Nsimilarity <- apply(stan_data$similarity,1,function(x)sum(!is.na(x)))
stan_data$Nwindow <- length(similarity)
stan_data$max_data <- max_data
stan_data$similarity[is.na(stan_data$similarity)] <- -1

save(stan_data, file = file.path(args_dir[['script_dir']],'stan_data_threshold.rda'))


cat('---------------- fit bimodal model to similarities ---------------------')
load(file.path(args_dir[['script_dir']],'stan_data_threshold.rda'))
stan_data$similarity <- log(stan_data$similarity)
stan_data$similarity[is.na(stan_data$similarity)] <- -1
stan_code <- "
data {
 int<lower = 0> N;
 vector[N] y;
}

parameters {
  ordered[2] mu;
  real<lower=0> sigma[2];
  real<lower=0, upper=1> theta;
}

model {
 sigma ~ exponential(5.0);
 mu ~ normal(0, 1);
 theta ~ beta(5, 5);
 for (n in 1:N)
   target += log_mix(theta,
                     normal_lpdf(y[n] | mu[1], sigma[1]),
                     normal_lpdf(y[n] | mu[2], sigma[2]));
}

"

fit <- list()
for (i in 1:99) {
  cat('fit window ', i, '\n')
  tmp <- list()
  tmp$y <- (stan_data$similarity[i,1:stan_data$Nsimilarity[i]])
  tmp$N <- stan_data$Nsimilarity[i]
  fit[[i]] <- stan(model_code=stan_code,
                         data=tmp,
                         chains=1, seed=42)
  
}

save(fit,file=file.path(args_dir[['script_dir']],'indep_exp5_fit.rda'))


cat('---------------- check fit ---------------------')
load(file.path(args_dir[['script_dir']],'indep_exp5_fit.rda'))

rh= c()
ness= c()
tess= c()
bess= c()
for (i in 1:99) {
  cat('process fit window ', i, '\n')
  pars <- rstan::extract(fit[[i]], pars=names(fit[[i]])[!grepl('lp_',names(fit[[i]]))])
  pars = do.call(cbind,pars)
  
  rh = c(rh,max(summary(fit[[i]])$summary[, "Rhat"]))
  ness = c(ness,min(summary(fit[[i]])$summary[, "n_eff"]))
  tess = c(tess,min(apply(pars, 2, ess_tail)))
  bess = c(bess,min(apply(pars, 2, ess_bulk)))
}

cat('range of rhat is ',range(rh),'\n')
cat('range of ess is ',range(ness),'\n')
cat('range of ess_tail is ',range(tess),'\n')
cat('range of ess_bulk is ',range(bess),'\n')


cat('---------------- extract outputs ---------------------')
df = data.table()
for (i in 1:99) {
  cat('process fit window ', i, '\n')
  tmp = data.table(summary(fit[[i]])$summary[1:5,c("2.5%","50%","97.5%")])
  tmp[,VAR:=rownames(summary(fit[[i]])$summary)[1:5]]
  tmp[,WIN:=i]
  df = rbind(df, tmp)
}


cat('---------------- visualise bimodal fit ---------------------')
load(file.path(args_dir[['script_dir']],'distance_couple_per_window.rda'))
ans[,PERC:=log(PERC)]
plot_list =list()
#  range(ans$PERC)
#  -0.4453817  0.0000000

# plot fit to similarities
for (i in 1:99) {
  cat('process window ',i, '\n')
  x <- seq(-0.45,0,length.out = 90 +1)
  df_plot_fit <- data.table(x=x,
                            c1 = df[VAR=='theta' & WIN ==i,]$`50%` * dnorm(x, df[VAR=='mu[1]' & WIN ==i,]$`50%`, df[VAR=='sigma[1]' & WIN ==i,]$`50%`),
                            c2 = (1-df[VAR=='theta' & WIN ==i,]$`50%`) *dnorm(x, df[VAR=='mu[2]' & WIN ==i,]$`50%`, df[VAR=='sigma[2]' & WIN ==i,]$`50%`))
  max_hist=  max(hist(ans[WINDOW==i,]$PERC,plot=F, breaks = seq(-0.45,0,length.out = 30 + 1), right = FALSE)$counts)
  scale <- max_hist/max(c(df_plot_fit$c1,df_plot_fit$c2))
  g = ggplot(ans[WINDOW==i,],aes(x=PERC))+
    geom_histogram(alpha=0.6)+
    geom_line(data=df_plot_fit, aes(x,c1 * scale), color='red') +
    geom_line(data=df_plot_fit, aes(x,c2 * scale), color='blue') +
    scale_x_continuous(breaks = seq(-0.45, 0,length.out = 30 + 1))+
    theme_bw() +
    labs(x="log similarities", y="counts")+ggtitle(paste0(windows_start[i], ' - ', windows_end[i]))+
    theme(axis.text.x=element_text(angle=45, hjust=1))
  
  plot_list[[i]] = ggplotGrob(g)
}


pdf(file.path(args_dir[['script_dir']],'plots','fitted_histogram_distance_couple_per_window_exp5_indep.pdf'),
    width = 30, height= 45)
do.call("grid.arrange", c(plot_list,ncol= 5))
dev.off()


# plot quantiles of fits
for (i in 1:99) {
  cat('process window ',i, '\n')
  x <- seq(-0.45,0,length.out = 90 +1)
  df_plot_fit <- data.table(x=x,
                            c1 = df[VAR=='theta' & WIN ==i,]$`50%` * dnorm(x, df[VAR=='mu[1]' & WIN ==i,]$`50%`, df[VAR=='sigma[1]' & WIN ==i,]$`50%`),
                            c2 = (1-df[VAR=='theta' & WIN ==i,]$`50%`) *dnorm(x, df[VAR=='mu[2]' & WIN ==i,]$`50%`, df[VAR=='sigma[2]' & WIN ==i,]$`50%`))
  g = ggplot(df_plot_fit)+
    geom_line( aes(x,c1), color='red') +
    geom_line(aes(x,c2), color='blue') +
    geom_vline(aes(xintercept = qnorm(0.05,  df[VAR=='mu[1]' & WIN ==i,]$`50%`, df[VAR=='sigma[1]' & WIN ==i,]$`50%`)),
               linetype=2, color='red') +
    geom_vline(aes(xintercept = qnorm(0.5,  df[VAR=='mu[1]' & WIN ==i,]$`50%`, df[VAR=='sigma[1]' & WIN ==i,]$`50%`)),
               linetype=2, color='red') +
    geom_vline(aes(xintercept = qnorm(0.95,  df[VAR=='mu[1]' & WIN ==i,]$`50%`, df[VAR=='sigma[1]' & WIN ==i,]$`50%`)),
               linetype=2, color='red') +
    geom_vline(aes(xintercept = qnorm(0.05,  df[VAR=='mu[2]' & WIN ==i,]$`50%`, df[VAR=='sigma[2]' & WIN ==i,]$`50%`)),
               linetype=2, color='blue') +
    geom_vline(aes(xintercept = qnorm(0.5,  df[VAR=='mu[2]' & WIN ==i,]$`50%`, df[VAR=='sigma[2]' & WIN ==i,]$`50%`)),
               linetype=2, color='blue') +
    geom_vline(aes(xintercept = qnorm(0.95,  df[VAR=='mu[2]' & WIN ==i,]$`50%`, df[VAR=='sigma[2]' & WIN ==i,]$`50%`)),
               linetype=2, color='blue') +
    theme_bw() +
    labs(x="log similarities", y="counts")+ggtitle(paste0(windows_start[i], ' - ', windows_end[i]))+
    theme(axis.text.x=element_text(angle=45, hjust=1))
  
  plot_list[[i]] = ggplotGrob(g)
  }


pdf(file.path(args_dir[['script_dir']],'plots','density_distance_couple_per_window_quantiles_exp5_indep.pdf'),
    width = 30, height= 45)
do.call("grid.arrange", c(plot_list,ncol= 5))
dev.off()

# mean and sd over windows
ggplot(df[grep('mu',VAR)],aes(x=WIN))+
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
  geom_line(aes(y = `50%`))+
  facet_wrap(VAR~.,ncol = 1) +
  labs(x='window', y='means')
ggsave(file.path(args_dir[['script_dir']],'plots','component_mean_vs_window_exp5_indep.pdf'),
       width = 6, height= 8)

ggplot(df[grep('sigma',VAR)],aes(x=WIN))+
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
  geom_line(aes(y = `50%`))+
  facet_wrap(VAR~.,ncol = 1) +
  labs(x='window', y='standard deviations')
ggsave(file.path(args_dir[['script_dir']],'plots','component_sd_vs_window_exp5_indep.pdf'),
       width = 6, height= 8)


cat('---------------- find threshold ---------------------')
cat('---------------- method 1: intersection ---------------------')

for (i in 1:99) {
  cat('process fit window ', i, '\n')
  if(i==1)threshold = list()
  tmp =c()
  pars = rstan::extract(fit[[i]], pars=c('mu','sigma','theta'))
  for (j in 1:nrow(pars$mu)) {
    f <- function(x) dlnorm(x, pars$mu[j,1], pars$sigma[j,1]) * pars$theta[j] -
      dlnorm(x, pars$mu[j,2], pars$sigma[j,2]) * (1-pars$theta[j])
    tryCatch({
      tmp <- c(tmp, uniroot(f, interval=c(0.7,1))$root)
    }, error=function(e){cat("ERROR:",conditionMessage(e), "\n",
                             'on window ',i, ' iteration ', j , '\n')})
   
  }  
  threshold[[i]] = tmp
}

threshold = lapply(threshold,function(x){quantile(x,probs = c(0.025,0.5,0.975))})
threshold = data.table(do.call(rbind, threshold))
save(threshold, file = file.path(args_dir[['script_dir']],paste0('threshold_dip_indep_exp5.rda')))


cat('---------------- method 2: 5% percentile ---------------------')
for (i in 1:99) {
  cat('process fit window ', i, '\n')
  if(i==1)threshold = list()
  tmp =c()
  pars = rstan::extract(fit[[i]], pars=c('mu','sigma','theta'))
  for (j in 1:nrow(pars$mu)) {
     tmp = c(tmp, qlnorm(0.05, pars$mu[j,2], pars$sigma[j,2]))
     
  }
  threshold[[i]] = tmp
}


threshold = lapply(threshold,function(x){quantile(x,probs = c(0.025,0.5,0.975))})
threshold = data.table(do.call(rbind, threshold))
save(threshold, file = file.path(args_dir[['script_dir']],paste0('threshold_5quantile_indep_exp5.rda')))


cat('---------------- method 3: median ---------------------')
for (i in 1:99) {
  cat('process fit window ', i, '\n')
  if(i==1)threshold = list()
  tmp =c()
  pars = rstan::extract(fit[[i]], pars=c('mu','sigma','theta'))
  for (j in 1:nrow(pars$mu)) {
    tmp = c(tmp, qlnorm(0.5, pars$mu[j,2], pars$sigma[j,2]))
    
  }
  threshold[[i]] = tmp
}


threshold = lapply(threshold,function(x){quantile(x,probs = c(0.025,0.5,0.975))})
threshold = data.table(do.call(rbind, threshold))
save(threshold, file = file.path(args_dir[['script_dir']],paste0('threshold_median_indep_exp5.rda')))
