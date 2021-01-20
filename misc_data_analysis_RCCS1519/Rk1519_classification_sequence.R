library(data.table)
library(ggplot2)
library(seqinr)
library(ggpubr)
library(tidyverse)
library(igraph)
library(dplyr)
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
dir.create(file.path(args_dir[['script_dir']], 'plots'))


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


# consensus sequence depth
if(exists(file.path(args_dir[['script_dir']] , 'plots','rccs_mrc_depth.rda'))){
  mapping_rccs <- readRDS('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/200422_PANGEA2_RCCS_mapped_samples.rds')
  dconsensus <- subset(mapping_rccs,select=c('PANGEA_ID','CONSENSUS','F'))
  dconsensus <- unique(dconsensus[CONSENSUS!="",])
  dconsensus[,PANGEA_ID:=paste0('RCCS_',PANGEA_ID)]
  mapping_mrc <- readRDS('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC/200422_PANGEA2_MRCUVRI_mapped_samples.rds')
  tmp <- subset(mapping_mrc,select=c('PANGEA_ID','CONSENSUS','F'))
  tmp <- unique(tmp[CONSENSUS!="",])
  tmp[,PANGEA_ID:=paste0('MRCUVRI_',PANGEA_ID)]
  dconsensus = rbind(dconsensus, tmp)
  dconsensus = subset(dconsensus, select=c('PANGEA_ID','F'))
  setnames(dconsensus, c('PANGEA_ID','F'), c('pangea_id','file_name'))
  tmp <-  data.table(pangea_id=names(alignment))
  dconsensus <- merge(dconsensus, tmp, by='pangea_id', all.y=T)
  dconsensus <- unique(dconsensus)
  
  ddepth = data.table()
  for (i in seq_len(nrow(dconsensus))) {
    if(i %% 100 == 1 | i==nrow(dconsensus)){
      cat('processing ', i, 'th out of ', nrow(dconsensus), ' files \n' )
    }
    sequence <- readLines(dconsensus$file[i])
    sequence <- paste(sequence[2:length(sequence)], collapse="")
    sequence <- strsplit(sequence,"")
    
    tmp_df <- dw[,{
      tmp = sequence[[1]][START:END]
      tmp = tmp[!grepl("[^A-Za-z]",tmp)]
      list(
        LEN_LETTER=length(tmp),
        HIGH_DEPTH=sum(tmp == toupper(tmp)),
        PANGEA_ID=dconsensus$pangea_id[i],
        FILE=dconsensus$file[i])
    }, by = c('WINDOW', 'START', 'END')]
    ddepth <- rbind(ddepth, tmp_df)
  }
  
  save(ddepth,file = file.path(args_dir[['script_dir']] , 'plots','rccs_mrc_depth.rda'))
  
}


# load depth
load(file.path(args_dir[['script_dir']] , 'plots','rccs_mrc_depth.rda'))

# high depth on > 250 bps
ddepth_select <- subset(ddepth,select=c('WINDOW','LEN_LETTER','HIGH_DEPTH','PANGEA_ID'))
ddepth_select <- ddepth_select[HIGH_DEPTH>= args_dir[['length_cutoff']]]
ddepth_select <- unique(ddepth_select)
# tmp = ddepth_select[,length(WINDOW), by='PANGEA_ID']

ddepth_select_w <- ddepth_select[,list(PANGEA_ID=unique(PANGEA_ID)),
                                 by=c('WINDOW')]
ddepth_select_w[,HD:=1]
ddepth_select_w <- data.table::dcast(ddepth_select_w, PANGEA_ID~WINDOW, value.var='HD')

tmp = data.table(V=rowSums(!is.na(ddepth_select_w[,2:100])))
ggplot(tmp,aes(V)) +
  geom_histogram() +
  labs(x=paste0('number of windows with high depth postions >= ',args_dir[['length_cutoff']]))+
  theme_bw()
ggsave(filename = file.path(args_dir[['script_dir']] , 'plots','number_of_window_with_high_depth_position_gt250.pdf'),
       width = 6, height = 4)


# load thresholds
if(args_dir[['method']] == 1){
  load(file.path(args_dir[['script_dir']],paste0('threshold_dip_indep_exp5.rda')))
}else if(args_dir[['method']] == 2){
  load(file.path(args_dir[['script_dir']],paste0('threshold_5quantile_indep_exp5.rda')))
}else if(args_dir[['method']] == 3){
  load(file.path(args_dir[['script_dir']],paste0('threshold_median_indep_exp5.rda')))
}
threshold <- threshold$`50%`


# distance & supports per pair per window
for (i in 1:args_dir[['windown']]) {
  # distance per window
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
  
  # add threshold and close status
  distancei[, threshold_bimodal:=threshold[i]]
  distancei[, CLOSE:= as.integer(PERC > threshold_bimodal)]
  
  # average distance per pair
  tmp  <- unique(subset(distancei, select= c('TAXA1','TAXA2','PERC')))
  tmp <- tmp[, list(PERC=mean(PERC, na.rm=T)), by=c('TAXA1','TAXA2')]
  setnames(tmp, 'PERC', paste0('PERC',i))
  if( i ==1 ){
    distance <- copy(tmp)
  }else{
    distance <- merge(distance, tmp, by =c('TAXA1','TAXA2'),all=T)
  }
  gc()
  
  # number of support per pair
  tmp  <- unique(subset(distancei, select= c('TAXA1','TAXA2','CLOSE')))
  tmp <-  tmp[, list(CLOSE=mean(CLOSE,na.rm = T)), by=c('TAXA1','TAXA2')]
  setnames( tmp, 'CLOSE', paste0('CLOSE',i))
  if( i ==1 ){
    dclassif <- copy( tmp)
  }else{
    dclassif <- merge(dclassif,  tmp, by =c('TAXA1','TAXA2'),all=T)
  }
  gc()
}


save(distance,file = paste0('distance_overall_method', args_dir[['method']],'_sequence_level_v4.rda'))
save(dclassif,file = paste0('support_overall_method', args_dir[['method']],'_sequence_level_v4.rda'))
save(dinfo, tmp_couple, file = file.path(args_dir[['script_dir']],'ptid_info.rda'))

