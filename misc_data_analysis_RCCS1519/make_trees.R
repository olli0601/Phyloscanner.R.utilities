library(data.table)
library(dplyr)
library(tidyverse)
library(seqinr)
if(1)
{
  args_dir <- list()
  args_dir[['data_dir']] <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  args_dir[['out_dir']] <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY'
  args_dir[['script_dir']] <- '~/phyloscanner/data_analysis_RCCS1519'
  args_dir[['infile']]  <- 'clusters.rda'
  args_dir[['job_tag']] <- 'windowsize500_batchsize100'
}


args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-data_dir')
  stopifnot(args_line[[3]]=='-out_dir')
  stopifnot(args_line[[5]]=='-script_dir')
  stopifnot(args_line[[7]]=='-infile')
  stopifnot(args_line[[9]]=='-job_tag')
  args_dir <- list()
  args_dir[['data_dir']] <- args_line[[2]]
  args_dir[['out_dir']] <- args_line[[4]]
  args_dir[['script_dir']] <- args_line[[6]]
  args_dir[['infile']] <- args_line[[8]]
  args_dir[['job_tag']] <- args_line[[10]]
} 

setwd( paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ))
# make pty.runs
# find bam files
pty.runs <- data.table(SAMPLE_ID=list.files(args_dir[['data_dir']], pattern='bam$',recursive = TRUE))
pty.runs[,SAMPLE_ID:=paste0(args_dir[['data_dir']],SAMPLE_ID)]
pty.runs = pty.runs[grepl('remap',SAMPLE_ID) & !grepl('PreDedup',SAMPLE_ID),]

# find pangea id per bam file
# mapping_rccs <- readRDS('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/200422_PANGEA2_RCCS_mapped_samples.rds')
mapping_rccs <- readRDS(file.path('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/200422_PANGEA2_RCCS_selected_samples.rds'))
map_bam_to_pangeaid <- subset(mapping_rccs,select=c('PANGEA_ID','REMAP_BAM'))
map_bam_to_pangeaid[,PANGEA_ID:=paste0('RCCS_',PANGEA_ID)]
# mapping_mrc <- readRDS('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC/200422_PANGEA2_MRCUVRI_mapped_samples.rds')
mapping_mrc <- readRDS('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC/200422_PANGEA2_MRCUVRI_selected_samples.rds')
tmp <- subset(mapping_mrc,select=c('PANGEA_ID','REMAP_BAM'))
tmp[,PANGEA_ID:=paste0('MRCUVRI_',PANGEA_ID)]
map_bam_to_pangeaid = rbind(map_bam_to_pangeaid, tmp)
pty.runs = merge(pty.runs, map_bam_to_pangeaid, by.x='SAMPLE_ID',by.y='REMAP_BAM',all.x=T)
# remove if pangea id not found 
cat(nrow(pty.runs[is.na(PANGEA_ID)]), ' PANGEA IDs are NA, remove \n')
# 2774  PANGEA IDs are NA, remove 
pty.runs = pty.runs[!is.na(PANGEA_ID)]

# separate directory in sample id
pty.runs[,SAMPLE_ID:=gsub(args_dir[['data_dir']],'',SAMPLE_ID)]
# set(pty.runs, NULL, 'SAMPLE_ID', pty.runs[, gsub('\\.bam','',SAMPLE_ID)])
# pty.runs[,FOLDER:=dirname(SAMPLE_ID)]
# pty.runs[,SAMPLE_ID:=basename(SAMPLE_ID)]
# pty.runs <- pty.runs[FOLDER!='PANGEA_RCCS_comparisonstudy']
# pty.runs[,PID:=sub('(.*-.*)-.*',"\\1",SAMPLE_ID)]



# add UNIT_ID to pty.runs
id.dt <- data.table(read.csv(file.path(args_dir[['data_dir']],'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')))
id.dt <- subset(id.dt,select = c("pt_id","sex","pangea_id"))
id.dt[,pangea_id:=paste0('RCCS_',pangea_id)]
tmp <- data.table(read.csv(file.path(args_dir[['data_dir']],'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')))
tmp <- subset(tmp,select = c("pt_id","sex","pangea_id"))
tmp[,pangea_id:=paste0('MRCUVRI_',pangea_id)]
id.dt <- rbind(id.dt,tmp)
id.dt <- unique(id.dt)
colnames(id.dt) = c("UNIT_ID","SEX","PANGEA_ID")
pty.runs <- merge(pty.runs, id.dt, by="PANGEA_ID",all.x=T)
cat(nrow(pty.runs[is.na(UNIT_ID)]), ' Person IDs are NA, remove \n')
# 92  Person IDs are NA, remove 
pty.runs <- pty.runs[!is.na(UNIT_ID)]


# anonymisation
aid = read.csv('important_anonymisation_keys_210119.csv')
aid = subset(data.table(aid),select=c('PT_ID','AID'))
setnames(aid,c('PT_ID','AID'),c('UNIT_ID','RENAME_ID'))
pty.runs <- merge(pty.runs, aid, by='UNIT_ID',all.x=T)
pty.runs[,FQ:=seq_len(length(SAMPLE_ID)), by='RENAME_ID']
pty.runs[,FQ:=paste0('fq',FQ)]
pty.runs[,RENAME_ID:= paste0(RENAME_ID,'-',FQ)]
pty.runs[,FQ:=NULL]


# load clusters of close pairs and redefine IDCLU
load(args_dir[['infile']])  
tmp2 = unique(subset(df,select=c('CLU_SIZE','IDCLU')))
setkey(tmp2, CLU_SIZE)
tmp2[, IDCLU2:=seq_len(nrow(tmp2))]  
tmp2[, CLU_SIZE:=NULL]
rtc = copy(df)
rtc	<- subset( merge(rtc, tmp2, by='IDCLU'))
rtc[,IDCLU:=NULL]
setnames(rtc,'IDCLU2','IDCLU')
setkey(rtc, IDCLU)
cat(max(rtc$IDCLU),' clusters of ',length(unique(rtc$ID)), ' individuals, the ten largest cluster sizes are ', 
    unique(subset(rtc,select=c('CLU_SIZE','IDCLU')))$CLU_SIZE[(max(rtc$IDCLU)-10+1):max(rtc$IDCLU)], ', ', nrow(rtc)-length(unique(rtc$ID)),
    'individuals appear in more than one clusters','\n')
# 808  clusters of  2364  individuals, the ten largest cluster sizes are  17 18 21 22 26 37 38 45 50 55 ,  150 individuals appear in more than one clusters 

#	add couples that are not a close pair
load(file.path(args_dir[['script_dir']],'RakaiPangeaMetaData_v2.rda'))
tmp<- sort(unique(as.character(pty.runs[['UNIT_ID']])))	
rp = data.table(unique(subset(coupdat,select=c('male.RCCS_studyid', 'female.RCCS_studyid'))))
setnames(rp, c('male.RCCS_studyid', 'female.RCCS_studyid'), c('MALE_RID','FEMALE_RID'))
rp[,FEMALE_RID:=paste0('RK-',FEMALE_RID)]
rp[,MALE_RID:=paste0('RK-',MALE_RID)]
rp<- unique(subset(rp, FEMALE_RID%in%tmp & MALE_RID%in%tmp, select=c(FEMALE_RID,MALE_RID)))
# rp[FEMALE_RID==MALE_RID]
rp<- rp[FEMALE_RID!=MALE_RID] # not sure why
setnames(rp, 'FEMALE_RID', 'ID')
rp<- merge(rp, subset(rtc, select=c(ID, IDCLU)), all.x=1, by='ID')
setnames(rp, c('ID','IDCLU','MALE_RID'), c('FEMALE_RID', 'FEMALE_IDCLU','ID'))
rp<- merge(rp, subset(rtc, select=c(ID, IDCLU)), all.x=1, by='ID')
setnames(rp, c('ID','IDCLU'), c('MALE_RID', 'MALE_IDCLU'))
tmp<- rp[, which(is.na(FEMALE_IDCLU) & !is.na(MALE_IDCLU))]
set(rp, tmp, 'FEMALE_IDCLU', rp[tmp, MALE_IDCLU])
tmp<- rp[, which(!is.na(FEMALE_IDCLU) & is.na(MALE_IDCLU))]
set(rp, tmp, 'MALE_IDCLU', rp[tmp, FEMALE_IDCLU])
tmp<- rp[, which(is.na(FEMALE_IDCLU) & is.na(MALE_IDCLU))]
set(rp, tmp, c('FEMALE_IDCLU','MALE_IDCLU'), rtc[, max(IDCLU)]+seq_along(tmp))
setnames(rp, c('MALE_IDCLU'), c('IDCLU'))
rp<- subset(melt(rp, id.vars=c('IDCLU'), measure.vars=c('MALE_RID','FEMALE_RID'), value.name='ID'), select=c(ID, IDCLU))
rp <- unique(rp)
set(rtc, NULL, 'CLU_SIZE', NULL)
rtc<- unique(rbind(rtc, rp))
tmp<- rtc[, list(CLU_SIZE=length(ID)), by='IDCLU']
setkey(tmp, CLU_SIZE)
tmp[, IDCLU2:=seq_len(nrow(tmp))]	
rtc<- merge(rtc, tmp, by='IDCLU')
rtc[, IDCLU:=NULL]
setnames(rtc, 'IDCLU2', 'IDCLU')
setkey(rtc,IDCLU)
cat(max(rtc$IDCLU),' clusters of ',length(unique(rtc$ID)), ' individuals, the ten largest cluster sizes are ', 
    unique(subset(rtc,select=c('CLU_SIZE','IDCLU')))$CLU_SIZE[(max(rtc$IDCLU)-10+1):max(rtc$IDCLU)], ', ', nrow(rtc)-length(unique(rtc$ID)),
    'individuals appear in more than one clusters','\n')
# 1059  clusters of  3032  individuals, the ten largest cluster sizes are  21 22 22 30 31 42 49 50 52 66 ,  242 individuals appear in more than one clusters 

#	merge up to 8 individuals into the same phyloscanner run
tn	<- 8
tmp	<- unique(subset(rtc, select=c(IDCLU, CLU_SIZE)))
tmp[, PTY_RUN:= IDCLU]
tmp[, PTY_SIZE:= CLU_SIZE]
# https://stackoverflow.com/questions/49076769/dplyr-r-cumulative-sum-with-reset
sum_reset_at <- function(thresh) {
  function(x) {
    accumulate(x, ~if_else(.x+.y>thresh, .y, .x+.y))
  }  
}
tmp = tmp %>% mutate(PTY_SIZE2 = sum_reset_at(tn)(PTY_SIZE)) 

tmp2 = which(diff(tmp$PTY_SIZE2)<0) + 1
tmp$PTY_RUN2 = 0 
tmp[c(tmp2,max(tmp2):nrow(tmp)),PTY_RUN2:=1]
tmp[,PTY_RUN2:=cumsum(PTY_RUN2)]
# check
stopifnot(sum(tmp[,cumsum(PTY_SIZE),by='PTY_RUN2']$V1 !=tmp$PTY_SIZE2)==0)
tmp[,PTY_SIZE2:=max(PTY_SIZE2),by='PTY_RUN2']
tmp[,PTY_RUN2:=PTY_RUN2+1]
set(tmp,NULL,c('PTY_RUN', 'PTY_SIZE'),NULL)
setnames(tmp,c('PTY_RUN2', 'PTY_SIZE2'), c('PTY_RUN', 'PTY_SIZE'))
rtc = merge(rtc,tmp,by=c('IDCLU', 'CLU_SIZE'))

#	for all individuals in a cluster
#	find 3 closest others

if(file.exists('average_distance_perpair.rds')){
  ddist = readRDS('average_distance_perpair.rds')
}else{
  load(paste0('distance_overall_method2_sequence_level_v4.rda'))
  ddist = subset(distance,select=c('TAXA1','TAXA2'))
  ddist$DIST = rowMeans(distance[,3:ncol(distance)],na.rm = T)
  setnames(id.dt, colnames(id.dt),paste0(colnames(id.dt),'1'))
  ddist = merge(ddist, id.dt, by.x='TAXA1', by.y='PANGEA_ID1',all.x=T)
  setnames(id.dt, colnames(id.dt),gsub('1','2',colnames(id.dt)))
  ddist = merge(ddist, id.dt, by.x='TAXA2', by.y='PANGEA_ID2',all.x=T)
  setnames(id.dt, colnames(id.dt),gsub('2','',colnames(id.dt)))
  ddist <- ddist[,list(DIST=mean(DIST)),by=c('UNIT_ID1','UNIT_ID2')]
  saveRDS(ddist,file='average_distance_perpair.rds') 
}

#
if(file.exists('closest_3_percluster.rds')){
  dcl =readRDS('closest_3_percluster.rds')
}else{
  dcl <- rtc[,{
    tmp = ddist[UNIT_ID1 %in% ID, c('UNIT_ID2','DIST')]
    tmp2 = ddist[UNIT_ID2 %in% ID, c('UNIT_ID1','DIST')]
    tmp = rbind(tmp,tmp2,use.names=F)
    colnames(tmp) = c('UNIT_ID','DIST')
    tmp = tmp[!UNIT_ID%in% ID,]
    tmp = tmp[,list(mean(DIST)),by='UNIT_ID']
    tmp=tmp[order(V1),]
    list(ID_CLOSE1=tmp$UNIT_ID[nrow(tmp)],DISTANCE1=tmp$V1[nrow(tmp)],
         ID_CLOSE2=tmp$UNIT_ID[nrow(tmp)-1],DISTANCE2=tmp$V1[nrow(tmp)-1],
         ID_CLOSE3=tmp$UNIT_ID[nrow(tmp)-2],DISTANCE3=tmp$V1[nrow(tmp)-2])
  },by=c('CLU_SIZE','IDCLU')]
  saveRDS(dcl,file='closest_3_percluster.rds')
}

# add the 3 closest individuals to each cluster
tmp = subset(dcl,select=c('CLU_SIZE', 'IDCLU',  'ID_CLOSE1',  'ID_CLOSE2',  'ID_CLOSE3'))
tmp = melt(tmp, id.vars=c('CLU_SIZE', 'IDCLU'),variable.name= 'CLOSEST', value.name='ID')
tmp[,CLOSEST:=as.integer(gsub('ID_CLOSE','',CLOSEST))]
tmp2 <- unique(subset(rtc,select=c('IDCLU','PTY_RUN')))
set(rtc, NULL ,c('PTY_SIZE', 'PTY_RUN'), NULL)
set(tmp, NULL, 'CLOSEST',NULL)
tmp[, ID_TYPE:='control']
rtc[, ID_TYPE:='target']
rtc <- rbind(rtc,tmp)
rtc <- unique(rtc)
rtc[,CLU_SIZE:=length(ID),by='IDCLU']
rtc <- merge(rtc, tmp2, by='IDCLU')
rtc[,PTY_SIZE:=length(ID),by='PTY_RUN']

cat(max(rtc$IDCLU),' clusters of ',length(unique(rtc$ID)), ' individuals, the ten largest cluster sizes are ', 
    unique(subset(rtc,select=c('CLU_SIZE','IDCLU')))$CLU_SIZE[(max(rtc$IDCLU)-10+1):max(rtc$IDCLU)], ', ', nrow(rtc)-length(unique(rtc$ID)),
    'individuals appear in more than one clusters','\n')
# 1059  clusters of  3253  individuals, the ten largest cluster sizes are  24 25 25 33 34 45 52 53 55 69 ,  3096 individuals appear in more than one clusters 
tmp = unique(subset(rtc,select=c('PTY_RUN','PTY_SIZE')))
table(tmp$PTY_SIZE)
# 6   8   9  10  11  12  13  14  15  17  19  20  22  23  24  25  33  34  45  52 
# 1  38  26  14   7  84   6  38   1  16   1 159   1   1   2   2   1   1   1   1 
# 53  55  69 
# 1   1   1 
save(rtc,pty.runs,file='phyloscanner_inputs.rda')





# ignore the last part
#	for every individual in a cluster
#	find 15 closest others
# # distance from consensus sequences
# 
# mapping_rccs <- readRDS('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/200422_PANGEA2_RCCS_mapped_samples.rds')
# dconsensus <- subset(mapping_rccs,select=c('PANGEA_ID','CONSENSUS','F'))
# dconsensus <- unique(dconsensus[CONSENSUS!="",])
# dconsensus[,PANGEA_ID:=paste0('RCCS_',PANGEA_ID)]
# mapping_mrc <- readRDS('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC/200422_PANGEA2_MRCUVRI_mapped_samples.rds')
# tmp <- subset(mapping_mrc,select=c('PANGEA_ID','CONSENSUS','F'))
# tmp <- unique(tmp[CONSENSUS!="",])
# tmp[,PANGEA_ID:=paste0('MRCUVRI_',PANGEA_ID)]
# dconsensus = rbind(dconsensus, tmp)
# dconsensus = subset(dconsensus, select=c('PANGEA_ID','F'))
# setnames(dconsensus, c('PANGEA_ID','F'), c('pangea_id','file_name'))
# dconsensus <- merge(dconsensus, id.dt, by='pangea_id', all.x=T)
# dconsensus <- subset(dconsensus,select=c('pt_id','file_name'))
# dconsensus <- unique(dconsensus)
# dconsensus[,ID:=seq_len(nrow(dconsensus))]
# # dconsensus[,sum(is.na(file_name))]
# # length(unique(dconsensus$pt_id))
# 
# lconsensus <- list()
# for (i in 1:nrow(dconsensus)) {
#   lconsensus[[i]] = read.fasta(dconsensus$file_name[i])
# }
# 
# dpairs = data.table(combn(dconsensus$ID,2))
# colnames(dpairs) = c('ID1', 'ID2')
# ddist<-dpairs[,{ 
#   seq1 <- lconsensus[[ID1]][[1]]
#   seq2 <- lconsensus[[ID2]][[1]]
#   pos <- which(seq1 !='-' & seq1 !='?' & seq1 !='n' & seq2 !='-' & seq2 !='?' &  seq2 !='n')
#   if (length(pos)==0 ){
#     len <- as.integer(-1.0)
#     match <- -1.0
#   }else{
#     seq1 <- seq1[pos]
#     seq2 <- seq2[pos]
#     len <- length(pos)
#     # cat(TAXA1,'\n',seq1,'\n',TAXA2,'\n',seq2,'\n')
#     # 
#     match <- sum(sapply(1:length(pos),function(pos){is.match(seq1[pos],seq2[pos])}))
#   }
#   list(LENGTH= len,
#        MATCHES= match)
# }, by=c('ID1','ID2')]
# 
# ddist[,SIMILARITY:=MATCHES/LENGTH]
# save(ddist, file = 'consensus_similarity.rda')

