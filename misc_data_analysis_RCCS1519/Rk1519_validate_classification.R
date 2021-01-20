library(data.table)
library(ggplot2)
library(seqinr)
library(ggpubr)
library(dplyr)
library(pammtools)# plot stepribbon

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
  args_dir[['method']] <- 2
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
  stopifnot(args_line[[17]]=='-method')
  args_dir <- list()
  args_dir[['data_dir']] <- args_line[[2]]
  args_dir[['out_dir']] <- args_line[[4]]
  args_dir[['window_size']] <- as.integer(args_line[[6]])
  args_dir[['batch_size']] <- as.integer(args_line[[8]])
  args_dir[['script_dir']] <- args_line[[10]]
  args_dir[['job_tag']] <- args_line[[12]]
  args_dir[['windown']] <- as.integer(args_line[[14]])
  args_dir[['length_cutoff']] <- as.integer(args_line[[16]])
  args_dir[['method']] <- as.integer(args_line[[18]])
} 

# dir
setwd( paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ))

cat('----------- plot threshold against distance (tree/ average)------------')
# load tree-based distance
tree_distance <- read.csv(file.path(args_dir[['script_dir']] ,'HIV_DistanceNormalisationOverGenome.csv'))
tree_distance <- data.table(tree_distance)

# map nongap positions to the hxb2
all = read.fasta(file.path(args_dir[['script_dir']] ,'HIV1_COM_2015_genome_DNA_ALL.fasta'))
all_hxb2 <- all$B.FR.83.HXB2_LAI_IIIB_BRU.K03455
all_hxb2_nongap <- which(!grepl('-',all$B.FR.83.HXB2_LAI_IIIB_BRU.K03455))
tree_distance[,POSITION_WRT_HXB2_WITH_GAP:=all_hxb2_nongap[tree_distance$POSITION_WRT_HXB2]]

# add missing locations to the hxb2
tree_distance <- rbind(tree_distance,
                  data.table(POSITION_WRT_HXB2_WITH_GAP=setdiff(1:max(tree_distance$POSITION_WRT_HXB2_WITH_GAP),tree_distance$POSITION_WRT_HXB2_WITH_GAP),
                             MEDIAN_PAIRWISE_DISTANCE_BETWEEN_STANDARD_REFS=NA_real_,
                             POSITION_WRT_HXB2=NA_real_))
setkey(tree_distance,POSITION_WRT_HXB2_WITH_GAP)

# load data
alignment <- read.fasta(file = file.path(args_dir[['data_dir']],"200422_PANGEA2_RCCSMRC_alignment.fasta"))
nsequence <- length(alignment)
npos <- unique(lengths(alignment))

# windows
windows_last_start <- ceiling(npos/100) * 100 - args_dir[['window_size']] + 1
windows_2last_end <- floor(npos/100) * 100
windows_start <- seq(1,windows_last_start,100)
windows_end <- c(seq(args_dir[['window_size']] ,windows_2last_end,100),npos)
windows <- seq_len(length(windows_start))

# mean tree_distance per window
tree_distance_per_window <- vector()
for (w in windows){
  tree_distance_per_window <- c( tree_distance_per_window,
                            mean( tree_distance$MEDIAN_PAIRWISE_DISTANCE_BETWEEN_STANDARD_REFS[windows_start[w]:windows_end[w]], na.rm = TRUE ) )
}

# average distance per window
load(paste0('distance_overall_method2_sequence_level_v4.rda'))
average_distance <-  distance%>%
  summarise_if(is.numeric, mean, na.rm = TRUE)
if(args_dir[['method']] == 1){
  load(file.path(args_dir[['script_dir']],paste0('threshold_dip_indep_exp5.rda')))
}else if(args_dir[['method']] == 2){
  load(file.path(args_dir[['script_dir']],paste0('threshold_5quantile_indep_exp5.rda')))
}else if(args_dir[['method']] == 3){
  load(file.path(args_dir[['script_dir']],paste0('threshold_median_indep_exp5.rda')))
}
ans <- data.table(THRES=threshold$`50%`,
                  DISTM=matrix(average_distance[1:99],ncol=1))
saveRDS(ans,file = file.path(args_dir[['script_dir']],'average_distance.rds'))

# compare distance and threshold
ans <- as.data.table(readRDS(file.path(args_dir[['script_dir']],'average_distance.rds')))
ans[, DISTT:=tree_distance_per_window]
colnames(ans) <- c('threshold','average distance','tree distance')
ans[, window:=seq_len(nrow(ans))]

ans <- melt(ans, id.var=c('window','threshold'))

# plot threshold against distances
ggplot(ans, aes(as.numeric(threshold), as.numeric(value), label=factor(window), color=factor(variable)))+
  geom_point()+
  scale_x_continuous(expand = c(0.05,0.05))   +
  scale_y_continuous(expand = c(0.05,0.05)) +
  theme_bw() +
  labs(x='threshold',
       y='distance',
       color='type') +
  scale_color_manual(values=c('tree distance'='royalblue3','average distance'='006400')) +
  theme(legend.position='bottom')
ggsave(file.path(args_dir[['script_dir']],'plots',paste0('distance_vs_threshold_method', args_dir[['method']],'.pdf')),width = 6, height = 6)

y=seq(0.01,0.05,0.005)
x=c(0.022,0.031,0.036,0.041,0.044,0.048,0.051,0.054,0.056)
fit2 <- lm(y~poly(x,4,raw=TRUE))
# plot(x,y)
# xx=seq(0,0.1,0.005)
# lines(xx, predict(fit2, data.frame(x=xx)), col="green")
xx=x
predict(fit2, data.frame(x=xx))-y
coefs = coef(fit2)
f = function(newdist, coefs){
  res <- coefs[1] + (coefs[2] * newdist) + (coefs[3] * newdist^2) + (coefs[4] * newdist^3)+ (coefs[5] * newdist^4)
  return(res)
}
f(0.041,coefs)

if(args_dir[['method']] == 1){
  load(file.path(args_dir[['script_dir']],paste0('threshold_dip_indep_exp5.rda')))
}else if(args_dir[['method']] == 2){
  load(file.path(args_dir[['script_dir']],paste0('threshold_5quantile_indep_exp5.rda')))
}else if(args_dir[['method']] == 3){
  load(file.path(args_dir[['script_dir']],paste0('threshold_median_indep_exp5.rda')))
}
thres = 1-threshold$`50%`
thres_map = sapply(thres,function(x)f(x,coefs))

# identify pol
start_pol <- min(all_hxb2_nongap[2085:5096])
end_pol <- max(all_hxb2_nongap[2085:5096])
start_pol_win <-which( windows_start -start_pol > 0)[1]
end_pol_win <-which(windows_end - end_pol> 0)[1] - 1
quantile(thres[start_pol_win:end_pol_win],probs=c(0.025,0.5,0.975))
# distance
# 0.01615810 0.02545088 0.02954131 
quantile(thres_map[start_pol_win:end_pol_win],probs=c(0.025,0.5,0.975))
# tree
# 0.009986778 0.011284782 0.013984869


# load depth
load(file.path(args_dir[['script_dir']] , 'plots','rccs_mrc_depth.rda'))

ddepth_select <- subset(ddepth,select=c('WINDOW','LEN_LETTER','HIGH_DEPTH','PANGEA_ID'))
ddepth_select <- ddepth_select[HIGH_DEPTH>= args_dir[['length_cutoff']]]
ddepth_select <- unique(ddepth_select)
# tmp = ddepth_select[,length(WINDOW), by='PANGEA_ID']

ddepth_select_w <- ddepth_select[,list(PANGEA_ID=unique(PANGEA_ID)),
                                 by=c('WINDOW')]
ddepth_select_w[,HD:=1]
ddepth_select_w <- data.table::dcast(ddepth_select_w, PANGEA_ID~WINDOW, value.var='HD')


# plot threshold and summarise similarity distribution for close and distant pairs by couple status
load(file.path(args_dir[['script_dir']],'ptid_info.rda'))
df_window_summary <- data.table()
for (i in seq_len(args_dir[['windown']])){
  cat('\n process window ',i , '\n')
  load(paste0('subsequence_window_',i,'_results.rda'))
  #
  distancei <- distancei[LENGTH >= args_dir[['length_cutoff']],]
  
  # select individuals with high depth postions >= 250
  ddepth_selecti = unique(subset(ddepth_select[WINDOW==i,],select='PANGEA_ID'))
  setnames(ddepth_selecti, 'PANGEA_ID', 'pangea_id')
  ddepth_selecti[,HD:=1]
  setnames(ddepth_selecti, colnames(ddepth_selecti), paste0(colnames(ddepth_selecti),'1'))
  distancei <- merge(distancei, ddepth_selecti, by.x='TAXA1', by.y='pangea_id1',all.x=T)
  setnames(ddepth_selecti, colnames(ddepth_selecti), gsub('1','2',colnames(ddepth_selecti)))
  distancei <- merge(distancei, ddepth_selecti, by.x='TAXA1', by.y='pangea_id2',all.x=T)
  setnames(ddepth_selecti, colnames(ddepth_selecti), gsub('2','',colnames(ddepth_selecti)))
  distancei <- distancei[HD1==1 & HD2==1,]
  
  # add couple
  setnames(dinfo, colnames(dinfo), paste0(colnames(dinfo),'1'))
  distancei <- merge(distancei,dinfo, by.x=c('TAXA1'), by.y=c('pangea_id1'), all.x=T) 
  setnames(dinfo, colnames(dinfo), gsub('1','2',colnames(dinfo)))
  distancei <- merge(distancei,dinfo, by.x=c('TAXA2'), by.y=c('pangea_id2'), all.x=T) 
  setnames(dinfo, colnames(dinfo), gsub('2','',colnames(dinfo)))
  
  # remove no pt_id
  distancei <- distancei[!is.na(pt_id1) & !is.na(pt_id2) & pt_id1!=pt_id2,]
  distancei <- merge(distancei,tmp_couple, by=c('pt_id1','pt_id2'), all.x=T)
  distancei[is.na(COUPLE),COUPLE:=0]
  
  # add close
  load(paste0('close_pairs.rda'))
  close_pairs <- subset(close_pairs, select=c('pt_id1','pt_id2'))
  close_pairs[,CLOSE:=1]
  distancei <- merge(distancei, close_pairs, by=c('pt_id1','pt_id2'), all.x=T)
  distancei[is.na(CLOSE),CLOSE:=0]
  # summarise by close and couple per window
  tmp <- distancei[,list(PREC=quantile(PERC, probs=c(0.05, 0.25, 0.5, 0.75, 0.95),na.rm=T),
                         PROB=c(0.05, 0.25, 0.5, 0.75, 0.95)),by=c('COUPLE','WINDOW','CLOSE')]
  tmp <- dcast(tmp,COUPLE+CLOSE+WINDOW~PROB,value.var = 'PREC')
  df_window_summary <- rbind(df_window_summary, tmp)
  
  gc()
}

threshold <- data.table(WINDOW=1:99,
                        THRESHOLD=threshold,
                        POS=windows_start)

df_window_summary <- merge(df_window_summary, threshold, by='WINDOW')

# save 
saveRDS(df_window_summary, file=file.path(args_dir[['script_dir']],paste0('distance_distribution_perwindow_bycoupleclose_method', args_dir[['method']],'.rds')))


df_window_summary <- readRDS(file.path(args_dir[['script_dir']],paste0('distance_distribution_perwindow_bycoupleclose_method', args_dir[['method']],'.rds')))


# plot

ggplot(df_window_summary)+
  geom_step(aes(x=POS,y=`0.5`,col=factor(CLOSE))) +
  geom_step(aes(x=POS,y=`THRESHOLD.50%`),color='black')+
  guides(col=F)+
  geom_stepribbon(aes(x=POS,ymin=`0.25`, ymax=`0.75`,fill=factor(CLOSE)),alpha=0.6)+
  geom_stepribbon(aes(x=POS,ymin=`0.05`, ymax=`0.95`,fill=factor(CLOSE)),alpha=0.2)+
  scale_fill_manual(values = c('1'='#C70039',
                               '0'='#2E86C1'),
                    labels=c('1'='yes',
                             '0'='no')) +
  scale_color_manual(values = c('1'='#C70039',
                                '0'='#2E86C1'),
                     labels=c('1'='yes',
                              '0'='no')) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = 'bottom')+
  labs(x='genome', y='distance',fill='close') +
  facet_grid(factor(COUPLE,levels = c(0,1), labels = c('not couples','couples'))~.)
ggsave(file.path(args_dir[['script_dir']],paste0('similarity_vs_window_method', args_dir[['method']],'.pdf')), width = 6, height = 6)




