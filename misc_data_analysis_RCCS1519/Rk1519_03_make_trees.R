library(data.table)
library(dplyr)
library(tidyverse)
library(seqinr)

rkuvri.make.anonymised.id <- function(){
  #' input person ids, return ananymise ids 
  
  # file names
  outfile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  infile.ind.rccs <- file.path(data.dir,'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.mrc <- file.path(data.dir,'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
  
  # load person IDs
  id.dt <- data.table(read.csv(infile.ind.rccs))
  id.dt <- subset(id.dt,select = c("pt_id"))
  tmp <- data.table(read.csv(infile.ind.mrc))
  tmp <- subset(tmp,select = c("pt_id"))
  id.dt <- rbind(id.dt,tmp)
  id.dt <- unique(id.dt)
  setkey(id.dt,pt_id)
  setnames(id.dt,colnames(id.dt),toupper(colnames(id.dt)))
  
  # set seed and sample identifier at random
  set.seed(42)
  id.dt[,AID:=sample(1:nrow(id.dt),replace = F)]
  id.dt[,AID:=formatC(AID, width=floor(log10(nrow(tmp))) +1, flag="0")]
  id.dt[,AID:=paste0('AID',AID)]
  write.csv(id.dt,file = outfile.ind.anonymised)
  Sys.chmod(outfile.ind.anonymised,mode = '0444')
}

rkuvri.make.phyloscanner.input.samples <- function()
{
  #' input person ids, selected samples, ananymise ids
  #' returns a data table summarising person id and bam files 
  #' 
  require(data.table)
  
  # file names
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  out.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_'
  infile.rccs <- file.path(data.dir, 'PANGEA2_RCCS/200422_PANGEA2_RCCS_selected_samples.rds')
  infile.mrc <- file.path(data.dir, 'PANGEA2_MRC/200422_PANGEA2_MRCUVRI_selected_samples.rds')
  infile.ind.rccs <- file.path(data.dir,'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.mrc <- file.path(data.dir,'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  
  # read selected files
  mapping_rccs <- readRDS(infile.rccs)
  mapping_rccs <- subset(mapping_rccs,select=c('PANGEA_ID','REMAP_BAM','REMAP_REF_FASTA'))
  mapping_rccs[,PANGEA_ID:=paste0('RCCS_',PANGEA_ID)]
  
  mapping_mrc <- readRDS(infile.mrc)
  mapping_mrc <- subset(mapping_mrc,select=c('PANGEA_ID','REMAP_BAM','REMAP_REF_FASTA'))
  mapping_mrc[,PANGEA_ID:=paste0('MRCUVRI_',PANGEA_ID)]
  
  ds <- rbind(mapping_rccs, mapping_mrc)
  ds <- subset(ds, !is.na(REMAP_BAM))
  
  # simplify REMAP_BAM and REMAP_REF_FASTA
  ds[, REMAP_BAM:=gsub(data.dir,'', REMAP_BAM)]
  ds[, REMAP_REF_FASTA:=gsub(data.dir,'', REMAP_REF_FASTA)]
  
  # double check file names as expected
  set(ds, NULL, 'SAMPLE_ID', ds[, gsub('\\.bam','',REMAP_BAM)])
  set(ds, NULL, 'SAMPLE_ID_CHECK', ds[, gsub('_ref.fasta','',REMAP_REF_FASTA)])
  stopifnot( ds[, all(SAMPLE_ID==SAMPLE_ID_CHECK)] )
  
  # finalise SAMPLE_ID
  set(ds, NULL, c('SAMPLE_ID_CHECK','REMAP_REF_FASTA','REMAP_BAM'), NULL)
  
  # read individual ids
  id.dt <- data.table(read.csv(infile.ind.rccs))
  id.dt <- subset(id.dt,select = c("pt_id","sex","pangea_id"))
  id.dt[,pangea_id:=paste0('RCCS_',pangea_id)]
  tmp <- data.table(read.csv(infile.ind.mrc))
  tmp <- subset(tmp,select = c("pt_id","sex","pangea_id"))
  tmp[,pangea_id:=paste0('MRCUVRI_',pangea_id)]
  id.dt <- rbind(id.dt,tmp)
  id.dt <- unique(id.dt)
  setnames(id.dt, c('pt_id','pangea_id'), c("UNIT_ID","PANGEA_ID"))
  set(id.dt, NULL, 'sex', NULL)
  
  # add anonymised ids
  aid <- as.data.table(read.csv(infile.ind.anonymised, stringsAsFactors = FALSE))
  aid <- subset(data.table(aid),select=c('PT_ID','AID'))
  setnames(aid,c('PT_ID','AID'),c('UNIT_ID','RENAME_ID'))
  id.dt <- merge(id.dt, aid, by='UNIT_ID', all.x=TRUE)
  # note: 92  Person IDs are NA, remove
  
  # merge with sample ids
  ds <- merge(id.dt, ds, by='PANGEA_ID')
  
  
  # finalise RENAME_ID
  ds[,FQ:=seq_len(length(SAMPLE_ID)), by='RENAME_ID']
  ds[,FQ:=paste0('fq',FQ)]
  ds[,RENAME_ID:= paste0(RENAME_ID,'-',FQ)]
  ds[,FQ:=NULL]
  
  # write processed samples 
  cat('\nWriting to file ', paste0(out.base,'phscinput_samples.rds') )
  saveRDS(ds, file=paste0(out.base,'phscinput_samples.rds'))
}

rkuvri.make.phyloscanner.input.runs <- function()
{
  #' input clusters, couples, and three closest persons for each cluster
  #' return a data table allocating each person to a run
  
  require(data.table)
  require(tidyverse)
  
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  out.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_'	
  potential.networks.analysis.dir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100"
  infile.potential.networks <- file.path(potential.networks.analysis.dir, 'clusters.rda')
  infile.clostest.three.per.cluster <- file.path(potential.networks.analysis.dir, 'closest_3_percluster.rds')
  infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
  
  # load data on individuals, sample locations, and rename ids 
  pty.runs <- readRDS(paste0(out.base,'phscinput_samples.rds'))
  
  # load potential transmission networks of close pairs and redefine IDCLU
  load( infile.potential.networks )  
  tmp <- unique(subset(df,select=c('CLU_SIZE','IDCLU')))
  setkey(tmp, CLU_SIZE)
  tmp[, IDCLU2:=seq_len(nrow(tmp))]  
  tmp[, CLU_SIZE:=NULL]
  rtc	<- subset( merge(df, tmp, by='IDCLU'),select=c('ID','IDCLU2','CLU_SIZE'))
  setnames(rtc,'IDCLU2','IDCLU')
  setkey(rtc, IDCLU)
  cat(max(rtc$IDCLU),' potential transmission networks of ',length(unique(rtc$ID)), ' individuals, the ten largest cluster sizes are ', 
      unique(subset(rtc,select=c('CLU_SIZE','IDCLU')))$CLU_SIZE[(max(rtc$IDCLU)-10+1):max(rtc$IDCLU)], ', ', nrow(rtc)-length(unique(rtc$ID)),
      'individuals appear in more than one clusters','\n')
  # 808  potential transmission networks of  2364  individuals, the ten largest cluster sizes are  17 18 21 22 26 37 38 45 50 55 ,  150 individuals appear in more than one clusters 
  
  #
  #	add couples that are not a close pair and redefine IDCLU
  #
  
  #	find which couples are in potential transmission networks
  load(infile.couple)	
  rp <- data.table(unique(subset(coupdat,select=c('male.RCCS_studyid', 'female.RCCS_studyid'))))
  setnames(rp, c('male.RCCS_studyid', 'female.RCCS_studyid'), c('MALE_RID','FEMALE_RID'))
  rp <- subset(rp, MALE_RID!=FEMALE_RID)
  rp[,FEMALE_RID:=paste0('RK-',FEMALE_RID)]
  rp[,MALE_RID:=paste0('RK-',MALE_RID)]
  tmp <- sort(unique(as.character(pty.runs[['UNIT_ID']])))	
  rp <- unique(subset(rp, FEMALE_RID%in%tmp & MALE_RID%in%tmp, select=c(FEMALE_RID,MALE_RID)))
  setnames(rp, 'FEMALE_RID', 'ID')
  rp <- merge(rp, subset(rtc, select=c(ID, IDCLU)), all.x=1, by='ID')
  setnames(rp, c('ID','IDCLU','MALE_RID'), c('FEMALE_RID', 'FEMALE_IDCLU','ID'))
  rp <- merge(rp, subset(rtc, select=c(ID, IDCLU)), all.x=1, by='ID')
  setnames(rp, c('ID','IDCLU'), c('MALE_RID', 'MALE_IDCLU'))
  rp[, COUP_ID := seq_len(nrow(rp))]
  
  #	reassign partners that are not in the same network as their partner
  tmp <- subset(rp, !is.na(FEMALE_IDCLU) & !is.na(MALE_IDCLU))
  tmp2 <- tmp[, list(NOT_IN_SAME = !any(FEMALE_IDCLU == MALE_IDCLU)),by=c('MALE_RID','FEMALE_RID')]
  tmp2 <- subset(tmp2, NOT_IN_SAME)	
  tmp2 <- merge(tmp2, rp, by=c('MALE_RID','FEMALE_RID'))	
  tmp <- rp[, which(COUP_ID %in% tmp2$COUP_ID)]
  set(rp, tmp, 'MALE_IDCLU', rp[tmp, FEMALE_IDCLU])
  set(tmp2, NULL, 'FEMALE_IDCLU', tmp2[,MALE_IDCLU])
  set(tmp2, NULL, 'NOT_IN_SAME', NULL)
  rp <- rbind(rp, tmp2)
  
  #	assign partners that are in no network to the same potential transmission network as their partner
  tmp <- rp[, which(is.na(FEMALE_IDCLU) & !is.na(MALE_IDCLU))]
  set(rp, tmp, 'FEMALE_IDCLU', rp[tmp, MALE_IDCLU])
  tmp <- rp[, which(!is.na(FEMALE_IDCLU) & is.na(MALE_IDCLU))]
  set(rp, tmp, 'MALE_IDCLU', rp[tmp, FEMALE_IDCLU])
  
  #	make a new potential transmission network for couples that are in no network so far
  tmp <- rp[, which(is.na(FEMALE_IDCLU) & is.na(MALE_IDCLU))]
  set(rp, tmp, c('FEMALE_IDCLU','MALE_IDCLU'), rtc[, max(IDCLU)]+seq_along(tmp))
  
  #	recalculate length of clusters
  setnames(rp, c('MALE_IDCLU'), c('IDCLU'))
  rp <- subset(melt(rp, id.vars=c('IDCLU'), measure.vars=c('MALE_RID','FEMALE_RID'), value.name='ID'), select=c(ID, IDCLU))
  rtc <- unique(rbind(subset(rtc,select=c('ID','IDCLU')), rp))
  tmp <- rtc[, list(CLU_SIZE=length(ID)), by='IDCLU']
  setkey(tmp, CLU_SIZE)
  tmp[, IDCLU2:=seq_len(nrow(tmp))]	
  rtc <- subset(merge(rtc, tmp, by='IDCLU'),select=c('ID','IDCLU2','CLU_SIZE'))
  setnames(rtc, 'IDCLU2', 'IDCLU')
  setkey(rtc,IDCLU)
  cat(max(rtc$IDCLU),' potential transmission networks of ',length(unique(rtc$ID)), ' individuals, the ten largest cluster sizes are ', 
      unique(subset(rtc,select=c('CLU_SIZE','IDCLU')))$CLU_SIZE[(max(rtc$IDCLU)-10+1):max(rtc$IDCLU)], ', ', nrow(rtc)-length(unique(rtc$ID)),
      'individuals appear in more than one clusters','\n')
  # 1059  potential transmission networks of  3032  individuals, the ten largest cluster sizes are  22 22 23 31 31 43 49 51 52 67 ,  277 individuals appear in more than one clusters 
  
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
  tmp <- tmp %>% mutate(PTY_SIZE2 = sum_reset_at(tn)(PTY_SIZE))
  tmp <- as.data.table(tmp)	
  tmp2 <- which(diff(tmp$PTY_SIZE2)<0) + 1
  tmp[, PTY_RUN2 := 0]	
  tmp[c(tmp2,max(tmp2):nrow(tmp)), PTY_RUN2:=1]
  tmp[,PTY_RUN2:=cumsum(PTY_RUN2)]
  stopifnot(sum(tmp[,cumsum(PTY_SIZE),by='PTY_RUN2']$V1 !=tmp$PTY_SIZE2)==0)
  tmp[,PTY_SIZE2:=max(PTY_SIZE2),by='PTY_RUN2']
  tmp[,PTY_RUN2:=PTY_RUN2+1]
  set(tmp,NULL,c('PTY_RUN', 'PTY_SIZE'),NULL)
  setnames(tmp,c('PTY_RUN2', 'PTY_SIZE2'), c('PTY_RUN', 'PTY_SIZE'))
  rtc <- merge(rtc,tmp,by=c('IDCLU', 'CLU_SIZE'))
  
  #	for all individuals in a cluster
  #	find 3 closest others
  
  if(file.exists(infile.clostest.three.per.cluster))
  {
    dcl <- readRDS(infile.clostest.three.per.cluster)
  }else{
    infile.average.distance.per.pair <- file.path(potential.networks.analysis.dir, 'average_distance_perpair.rds')
    if(file.exists(infile.average.distance.per.pair))
    {
      ddist <- readRDS(infile.average.distance.per.pair)
    }else{
      infile.average.distance.per.pair.per.window <- file.path(potential.networks.analysis.dir, 'distance_overall_method2_sequence_level_v4.rda')
      load(infile.average.distance.per.pair.per.window)
      ddist <- subset(distance,select=c('TAXA1','TAXA2'))
      ddist$DIST <- rowMeans(distance[,3:ncol(distance)],na.rm = TRUE)
      setnames(id.dt, colnames(id.dt),paste0(colnames(id.dt),'1'))
      ddist <- merge(ddist, id.dt, by.x='TAXA1', by.y='PANGEA_ID1',all.x=TRUE)
      setnames(id.dt, colnames(id.dt),gsub('1','2',colnames(id.dt)))
      ddist <- merge(ddist, id.dt, by.x='TAXA2', by.y='PANGEA_ID2',all.x=TRUE)
      setnames(id.dt, colnames(id.dt),gsub('2','',colnames(id.dt)))
      ddist <- ddist[,list(DIST=mean(DIST)),by=c('UNIT_ID1','UNIT_ID2')]
      saveRDS(ddist,file=infile.average.distance.per.pair) 
    }
    dcl <- rtc[,{
      tmp <- ddist[UNIT_ID1 %in% ID, c('UNIT_ID2','DIST')]
      tmp2 <- ddist[UNIT_ID2 %in% ID, c('UNIT_ID1','DIST')]
      tmp <- rbind(tmp,tmp2,use.names=F)
      colnames(tmp) <- c('UNIT_ID','DIST')
      tmp <- tmp[!UNIT_ID%in% ID,]
      tmp <- tmp[,list(mean(DIST)),by='UNIT_ID']
      tmp <- tmp[order(V1),]
      list(ID_CLOSE1=tmp$UNIT_ID[nrow(tmp)],DISTANCE1=tmp$V1[nrow(tmp)],
           ID_CLOSE2=tmp$UNIT_ID[nrow(tmp)-1],DISTANCE2=tmp$V1[nrow(tmp)-1],
           ID_CLOSE3=tmp$UNIT_ID[nrow(tmp)-2],DISTANCE3=tmp$V1[nrow(tmp)-2])
    },by=c('CLU_SIZE','IDCLU')]
    saveRDS(dcl,file=infile.clostest.three.per.cluster)
  }
  
  # add the 3 closest individuals to each potential transmission network
  tmp <- subset(dcl, select=c('CLU_SIZE', 'IDCLU', 'ID_CLOSE1', 'ID_CLOSE2', 'ID_CLOSE3'))
  tmp <- melt(tmp, id.vars=c('CLU_SIZE', 'IDCLU'), variable.name= 'CLOSEST', value.name='ID')
  tmp[,CLOSEST:=as.integer(gsub('ID_CLOSE','',CLOSEST))]
  tmp2 <- unique(subset(rtc,select=c('IDCLU','PTY_RUN')))
  set(rtc, NULL ,c('PTY_SIZE', 'PTY_RUN'), NULL)
  set(tmp, NULL, 'CLOSEST',NULL)
  tmp[, ID_TYPE:='control']
  rtc[, ID_TYPE:='target']
  rtc <- rbind(rtc,tmp)
  rtc <- unique(rtc)
  rtc[, CLU_SIZE:=length(ID),by='IDCLU']
  rtc <- merge(rtc, tmp2, by='IDCLU')
  rtc[, PTY_SIZE:=length(ID),by='PTY_RUN']
  
  cat(max(rtc$IDCLU),' clusters of ',length(unique(rtc$ID)), ' individuals, the ten largest cluster sizes are ', 
      unique(subset(rtc,select=c('CLU_SIZE','IDCLU')))$CLU_SIZE[(max(rtc$IDCLU)-10+1):max(rtc$IDCLU)], ', ', nrow(rtc)-length(unique(rtc$ID)),
      'individuals appear in more than one clusters','\n')
  # 1059  clusters of  3253  individuals, the ten largest cluster sizes are  25 25 26 34 34 46 52 54 55 70 ,  3131 individuals appear in more than one clusters 
  
  tmp <- unique(subset(rtc,select=c('PTY_RUN','PTY_SIZE')))
  table(tmp$PTY_SIZE)
  #   6   8   9  10  11  12  13  14  16  17  18  19  20  22  23  24  25  26  34  46
  #	1  38  28  15   7  90   5  38   1  15   1   1 156   1   1   1   2   1   2   1
  #	52  54  55  70
  #	1   1   1   1
  
  # combine with sample info
  setnames(rtc, 'ID', 'UNIT_ID')
  pty.runs <- merge(rtc, pty.runs, by='UNIT_ID',all.x=T)
  
  
  # write processed samples 
  cat('\nWriting to file ', paste0(out.base,'phscinput_runs.rds') )
  saveRDS(pty.runs, file=paste0(out.base,'phscinput_runs.rds'))		
}


rkuvri.pre.make.alignments.tree<- function()
{
  #
  #	set up working environment
  require(Phyloscanner.R.utilities)
  library(data.table)
  library(seqinr)
  library(DECIPHER)
  library(msa)
  
  args <- commandArgs(trailingOnly = TRUE)
  print(args)
  job_tag <- args[1]
  # file names
  set.seed(42)
  # downsample <- 50
  downsample <- NULL
  # job_tag <- 'refset1'
  # job_tag <- 'refset2'
  # job_tag <- 'refset3'
  # job_tag <- 'all'
  prog.pty <- '~/phyloscanner/phyloscanner_make_trees.py'
  HOME <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'	
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  potential.networks.analysis.dir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100"
  if(!is.null(downsample)){  
    in.dir <- file.path(HOME,paste0('210211_phsc_input',downsample))		
    work.dir <- file.path(HOME,paste0("210211_phsc_work",downsample))	
    out.dir <- file.path(HOME,paste0("210211_phsc_output",downsample))	
  }else{
    in.dir <- file.path(HOME,paste0('210211_phsc_input_',job_tag))		
    work.dir <- file.path(HOME,paste0("210211_phsc_work_",job_tag))	
    out.dir <- file.path(HOME,paste0("210211_phsc_output_",job_tag))	
  }
  package.dir <- file.path(.libPaths(),'Phyloscanner.R.utilities')
  dir.create(in.dir)
  dir.create(work.dir)
  dir.create(out.dir)
  infile.runs <- paste0(HOME,'210120_RCCSUVRI_phscinput_runs.rds')
  infile.consensus.oneeach <- file.path(in.dir,'2019_New_ConsensusGenomesOneEach_GeneCut.fasta')
  infile.hxb2.package <- file.path(package.dir,'HIV1_compendium_AD_B_CPX_v2.fasta')
  infile.consensus <- file.path(in.dir,'UgandaKenyaTanzaniaGenomes_GeneCut_TreeOrder.FASTA') 
  
  # load runs
  pty.runs <- readRDS(infile.runs)
  max.per.run <- 4900
  cat(system(paste0('cp -r ',gsub(paste0(downsample,'$|_[a-z]+$|_[a-z]+[0-9]$'),'',in.dir),'/2019_New_ConsensusGenomesOneEach_GeneCut.fasta ', 
                    gsub(paste0(downsample,'$|_[a-z]+$|_[a-z]+[0-9]$'),'',in.dir),'/UgandaKenyaTanzaniaGenomes_GeneCut_TreeOrder.FASTA ', in.dir)))
  
  # remove same bam files per run
  setorder(pty.runs,PTY_RUN,-ID_TYPE,UNIT_ID)
  tmp <- pty.runs[,duplicated(SAMPLE_ID),by='PTY_RUN']
  tmp <- tmp[,which(V1)]
  pty.runs <- pty.runs[-tmp,]
  tmp <- pty.runs[,length(SAMPLE_ID)-length(unique(SAMPLE_ID)),by='PTY_RUN']
  stopifnot(all(tmp$V1==0))
  
  # clean taxa names 
  consensus_seq_oneeach<- read.fasta(infile.consensus.oneeach)
  dn <- data.table(RAW_NAME=names(consensus_seq_oneeach)) 
  dn[,NUM:=as.numeric(gsub(".*\\((.*)\\).*", "\\1", RAW_NAME))]
  dn[,NAME:=gsub('\\(.*?\\)', '', RAW_NAME)]
  names(consensus_seq_oneeach) <- dn$NAME
  
  # remove sequences
  consensus_seq<- read.fasta(infile.consensus)
  # compare two sets of consensus, and combine
  tmp = sapply(consensus_seq_oneeach,function(x){any(sapply(consensus_seq, function(y){identical(as.character(x)[as.character(x)!='-'],as.character(y)[as.character(y)!='-'])}))})
  consensus_seq_oneeach <- consensus_seq_oneeach[setdiff(names(consensus_seq_oneeach),c(names(tmp[tmp==T]),'CON_OF_CONS'))]
  
  opn <- grep('^CON_N$|^CON_O$|^CON_P$|CONSENSUS_N$|CONSENSUS_O$|CONSENSUS_P$',names(consensus_seq_oneeach),value=T)
  consensus_seq_oneeach <- consensus_seq_oneeach[setdiff(names(consensus_seq_oneeach),opn)]
  
  # filter sequences
  sample_date <- as.numeric(sapply(strsplit(names( consensus_seq),split='\\.'),function(x)x[3]))
  id <- which(sample_date<=2010|sample_date>=1990)
  consensus_seq<- consensus_seq[names(consensus_seq)[id]]
  
  if(is.null(downsample)){
    infile.consensus.all <- file.path(in.dir,'ConsensusGenomes.fasta')
    if(file.exists(infile.consensus.all)){
      consensus_seq_all <- read.fasta(infile.consensus.all)
    }else{
      if(grepl('refset[0-9]$',job_tag)){
        consensus_seq_all <- read.fasta(file.path(gsub(paste0('_',job_tag,'$'),'',in.dir),'ConsensusGenomes.fasta'))
        tmp <- names(consensus_seq_all)
        if(grepl('1$',job_tag)){
          # consensus (A/D)
          tmp_select <- which(tmp %in%  paste0('REF_',names(consensus_seq_oneeach)))
          tmp_select <- tmp[tmp_select]
          tmp_select2 <- gsub('_|\\-','\\.',tmp_select)
          tmp_select2 <- sapply(strsplit(tmp_select2,split = '\\.'),function(x)x[3])
          tmp_select2 <- grepl('A|D',tmp_select2) | !grepl('[A-Z]+|[a-z]+',tmp_select2)
          tmp_select <- tmp_select[tmp_select2]
          tmp_select <- c(tmp_select, grep('HXB2',tmp,value = T), grep('CON_H',tmp,value = T))
        }else if(grepl('2$',job_tag)){
          # consensus 
          tmp_select <- which(tmp %in%  paste0('REF_',names(consensus_seq_oneeach)))
          tmp_select <- tmp[tmp_select]
          tmp_select <- c(tmp_select, grep('HXB2',tmp,value = T))
        }else if(grepl('3$',job_tag)){
          # consensus + local (A/D)
          tmp_select <- which(!tmp %in%  paste0('REF_',names(consensus_seq_oneeach)))
          tmp_select <- tmp[tmp_select]
          tmp_select2 <- gsub('^REF_','',tmp_select)
          tmp_select2 <- sapply(strsplit(tmp_select2,split = '\\.'),function(x)x[1])
          tmp_select2 <- grepl('A|D',tmp_select2) | !grepl('[A-Z]+|[a-z]+',tmp_select2)
          tmp_select <-  tmp_select[tmp_select2]
          tmp_select <- c(tmp_select, tmp[which(tmp %in%  paste0('REF_',names(consensus_seq_oneeach)))])
          tmp_select <- c(tmp_select, grep('HXB2',tmp,value = T))
        }
        consensus_seq_all <- consensus_seq_all[unique(tmp_select)]
        write.fasta(sequences=consensus_seq_all,
                    names=names(consensus_seq_all),
                    nbchar = 60,
                    file.out=infile.consensus.all)
      }else{
        # tmp <- msa(DNAStringSet(do.call(c,lapply(consensus_seq_all,function(x)paste(x,collapse = '')))),method="Muscle",type='dna')
        tmp <- AlignProfiles(DNAStringSet(do.call(c,lapply(consensus_seq,function(x)paste(x,collapse = '')))),
                             DNAStringSet(do.call(c,lapply(consensus_seq_oneeach,function(x)paste(x,collapse = '')))))
        consensus_seq_all <- as.list(tmp)
        consensus_seq_all <- lapply(consensus_seq_all,function(x){as.character(x)})
        names(consensus_seq_all) <- paste0('REF_',names(consensus_seq_all))
        write.fasta(sequences=consensus_seq_all,
                    names=names(consensus_seq_all),
                    nbchar = 60,
                    file.out=infile.consensus.all)
      }
    }
  }else{
    infile.consensus.all <- file.path(in.dir,paste0('ConsensusGenomes.fasta'))
    if(!file.exists(infile.consensus.all)){
      consensus_seq_all <- read.fasta(file.path(gsub(paste0(downsample,'$'),'',in.dir),'ConsensusGenomes.fasta'))
      hxb2 <- grep('HXB2',names(consensus_seq_all),value = T)
      consensus_seq_all <- consensus_seq_all[c(paste0('REF_',c(names(consensus_seq_oneeach),sample(names(consensus_seq)[!grepl('HXB2',names(consensus_seq))],downsample-1))),hxb2)]
      infile.consensus.all <- file.path(in.dir,paste0('ConsensusGenomes.fasta'))
      write.fasta(sequences=consensus_seq_all,
                  names=names(consensus_seq_all),
                  nbchar = 60,
                  file.out=infile.consensus.all)
    }else{
      consensus_seq_all <- read.fasta(infile.consensus.all)
    }
  }
  
  # take hxb2
  hxb2 <- grep('HXB2',names(consensus_seq_all),value = T)
  hxb2_seq <- consensus_seq_all[[hxb2]]
  
  # compare hxb2 
  package_seq <- read.fasta(infile.hxb2.package)
  package_hxb2_seq <- package_seq[[1]]
  package_hxb2_seq <- as.character(package_hxb2_seq)[as.character(package_hxb2_seq)!='-']
  hxb2_seq <- as.character(hxb2_seq)[as.character(hxb2_seq)!='-']
  stopifnot(all(hxb2_seq==package_hxb2_seq))
  
  # take root
  root.seq <-grep('^REF_CON_M$|REF_CONSENSUS_M$|^REF_CON_H$|REF_CONSENSUS_H$', names(consensus_seq_all),value = T)
  
  # adapt format
  pty.runs[ID_TYPE=='control',UNIT_ID:=paste0('CNTRL-',UNIT_ID)]
  pty.runs[ID_TYPE=='control',RENAME_ID:=paste0('CNTRL-',RENAME_ID)]
  pty.runs[,BAM:=paste0(data.dir,SAMPLE_ID,'.bam')]
  pty.runs[,REF:=paste0(data.dir,SAMPLE_ID,'_ref.fasta')]
  setkey(pty.runs,PTY_RUN,RENAME_ID)
  
  #
  #	define phyloscanner input args to generate read alignments 
  #	for each window and each run
  ptyi <- seq(800,9400,25)		
  pty.c	<- lapply(seq_along(ptyi), function(i)
  {
    cat('---------------------------------------------------------------- \n')
    print(i)
    pty.args <- list(prog.pty=prog.pty, 
                     prog.mafft='mafft', 						 
                     data.dir=data.dir, 
                     work.dir=work.dir, 
                     out.dir=out.dir, 
                     alignments.file=infile.consensus.all,
                     alignments.root=root.seq,
                     alignments.pairwise.to=hxb2,
                     window.automatic= '', 
                     merge.threshold=0, 
                     min.read.count=1, 
                     quality.trim.ends=23, 
                     min.internal.quality=23, 
                     merge.paired.reads=TRUE, 
                     discard.improper.pairs=TRUE,
                     no.trees=TRUE, 
                     dont.check.duplicates=FALSE,
                     dont.check.recombination=TRUE,
                     num.bootstraps=1,
                     all.bootstrap.trees=TRUE,
                     strip.max.len=350, 
                     min.ureads.individual=NA, 
                     win=c(ptyi[i],ptyi[i]+250,25,250),				 				
                     keep.overhangs=FALSE,
                     mem.save=0,
                     verbose=TRUE,					
                     select=NA,
                     default.coord=FALSE
    )											
    pty.c <- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)
    pty.c[, W_FROM:= ptyi[i]]
    # pty.c[, PTY_RUN:= as.integer(sub('.*ptyr([0-9])_.*','\\1',CMD))]
    pty.c
  })
  pty.c	<- do.call('rbind', pty.c)	
  setkey(pty.c,PTY_RUN,W_FROM)
  pty.c[, CASE_ID:= rep(1:max.per.run,times=ceiling(nrow(pty.c)/max.per.run))[1:nrow(pty.c)]]
  pty.c[, JOB_ID:= rep(1:ceiling(nrow(pty.c)/max.per.run),each=max.per.run)[1:nrow(pty.c)]]
  
  #	define PBS variables
  hpc.load			<- "module load intel-suite/2015.1 mpi raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"	# make third party requirements available	 
  hpc.select			<- 1						# number of nodes
  hpc.nproc			<- 1						# number of processors on node
  hpc.walltime		<- 123						# walltime
  hpc.q				<- "pqeelab"				# PBS queue
  hpc.mem				<- "6gb" 					# RAM	
  hpc.array			<- pty.c[, max(CASE_ID)]	# number of runs for job array
  #	define PBS header for job scheduler. this will depend on your job scheduler.
  pbshead		<- "#!/bin/sh"
  tmp			<- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
  pbshead		<- paste(pbshead, tmp, sep = "\n")
  tmp			<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
  pbshead 	<- paste(pbshead, tmp, sep = "\n")
  pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")	
  if(!is.na(hpc.array))
    pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')	
  if(!is.na(hpc.q)) 
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  pbshead 	<- paste(pbshead, hpc.load, sep = "\n")	
  cat(pbshead)
  
  #	create PBS job array
  for (i in 1:pty.c[, max(JOB_ID)]) {
    tmp<-pty.c[JOB_ID==i,]
    cmd<-tmp[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
    cmd<-cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]			
    cmd<-paste(pbshead,cmd,sep='\n')	
    #	submit job
    outfile<-gsub(':','',paste("readali",paste0('job',i),paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
    outfile<-file.path(work.dir, outfile)
    cat(cmd, file=outfile)
    cmd<-paste("qsub", outfile)
    cat(cmd)
    if(i==1){cat(system(cmd, intern= TRUE))}
    # if(i==1|i==pty.c[, max(JOB_ID)]){cat(system(cmd, intern= TRUE))}
  }
  
  
  # tree part
  require(data.table)
  require(Phyloscanner.R.utilities)
  #
  #	produce trees
  #
  # lightweight run
  if(1)	
  {
    hpc.select<- 1; hpc.nproc<- 1; hpc.walltime<- 3; hpc.mem<- "1850mb"; hpc.q<- NA
  }
  # midweight run 
  if(0)	
  {
    hpc.select<- 1; hpc.nproc<- 1; hpc.walltime<- 23; hpc.mem<- "1850mb"; hpc.q<- NA
  }
  # heavyweight run 
  if(0)	
  {
    # hpc.select<- 1; hpc.nproc<- 1; hpc.walltime<- 71; hpc.mem<- "5900mb"; hpc.q<- "pqeelab"
    hpc.select<- 1; hpc.nproc<- 1; hpc.walltime<- 71; hpc.mem<- "63850mb"; hpc.q<- NA
  }
  
  HOME <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'	
  max.per.run <- 4900
  # downsample <- 50
  # downsample <- 100
  # downsample <- 200
  # downsample <- 300
  # downsample <- 400
  # downsample <- 500
  downsample <- NULL
  # job_tag <- 'testnoalign'
  # job_tag <- 'testrealign'
  # job_tag <- 'refset1'
  # job_tag <- 'refset2'
  # job_tag <- 'refset3'
  
  iqtree.pr <- 'iqtree'
  iqtree.args			<- ifelse(hpc.nproc==1, '-m GTR+F+R6 -ntmax 1 -seed 42 -o REF_CON_H', 
                          paste0('-m GTR+F+R6 -ntmax ',hpc.nproc,' -seed 42 -o REF_CON_H'))
  
  if(!is.null(downsample)){  
    in.dir				<- file.path(HOME,paste0('210211_phsc_output',downsample))		
    out.dir				<- in.dir
    work.dir			<- file.path(HOME,paste0("210211_phsc_work",downsample))	
  }else{
    in.dir				<- file.path(HOME,paste0('210211_phsc_output_',job_tag))		
    out.dir				<- in.dir
    work.dir			<- file.path(HOME,paste0("210211_phsc_work_",job_tag))
  }
  
  infiles	<- data.table(FI=list.files(in.dir, pattern='fasta$', full.names=TRUE, recursive=TRUE))
  infiles[, FO:= gsub('.fasta$','',FI)]
  infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FI)))]
  infiles[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FI)))]		
  setkey(infiles, PTY_RUN, W_FROM)	
  
  infiles	<- subset(infiles, PTY_RUN<=14)
  print(infiles)	 		
  
  df<- infiles[, list(CMD=cmd.iqtree(FI, outfile=FO, pr=iqtree.pr, pr.args=iqtree.args)), by=c('PTY_RUN','W_FROM')]
  df[, ID:=ceiling(seq_len(nrow(df))/4)]
  df<- df[, list(CMD=paste(CMD, collapse='\n',sep='')), by='ID']
  
  #	create PBS job array
  if(nrow(df) > max.per.run){
    df[, CASE_ID:= rep(1:max.per.run,times=ceiling(nrow(df)/max.per.run))[1:nrow(df)]]
    df[, JOB_ID:= rep(1:ceiling(nrow(df)/max.per.run),each=max.per.run)[1:nrow(df)]]
    set(df, ID, NULL, NULL)
  }else{
    setnames(df, 'ID', 'CASE_ID')
    df[, JOB_ID:=1]
  }
  
  pbshead	<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=NULL)
  
  
  for (i in 1:df[, max(JOB_ID)]) {
    tmp<-df[JOB_ID==i,]
    hpc.array <- nrow(tmp)
    cmd<-tmp[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
    cmd<-cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]		
    tmp <- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')	
    tmp <- paste(tmp,'module load anaconda3/personal \n source activate phylo' ,sep='\n')
    cmd<-paste(tmp,cmd,sep='\n')	
    #	submit job
    outfile	<- paste("srx",paste0('job',i,'_ds',downsample),paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.')
    outfile<-file.path(work.dir, outfile)
    cat(cmd, file=outfile)
    cmd<-paste("qsub", outfile)
    cat(cmd)
  }
  
  
  library(tidyverse)
  library(ggtree)
  library(data.table)
  library(seqinr)
  library(ggplot2)
  library(viridis)
  library(treedater)
  library(treeio)
  
  # files
  set.seed(42)
  args <- commandArgs(trailingOnly = TRUE)
  print(args)
  # downsample <- as.numeric(args[1])
  # downsample <- 50
  job_tag=args[1]
  downsample <- NULL
  # job_tag <- 'refset1'
  # job_tag <- 'refset2'
  # job_tag <- 'refset3'
  if(is.null(downsample)){
    p <- paste0('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_output_',job_tag,'/ptyr1_trees')
    
  }else{
    p <- paste0('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_output',downsample,'/ptyr1_trees')
  }
  # p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_output/ptyr400_trees'
  dir.create(file.path(p,'plots'))
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  infile.ind.rccs <- file.path('/rds/general/project/ratmann_pangea_deepsequencedata/live/','PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.mrc <- file.path('/rds/general/project/ratmann_pangea_deepsequencedata/live/','PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
  infile.phscinput <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_phscinput_samples.rds'
  infile.consensus <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_input/ConsensusGenomes.fasta'
  infile.consensus.oneeach <- file.path('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_input/2019_New_ConsensusGenomesOneEach_GeneCut.fasta')
  
  # load sampling date
  ds <- readRDS(infile.phscinput)
  dsdt <- data.table(read.csv(infile.ind.rccs))
  dsdt <- subset(dsdt, select=c('pt_id','visit_dt','pangea_id'))
  tmp <- data.table(read.csv(infile.ind.mrc))
  tmp <- subset(tmp, select=c('pt_id','visit_dt','pangea_id'))
  dsdt <- rbind(dsdt, tmp)
  setnames(dsdt, colnames(dsdt), toupper(colnames(dsdt)))
  tmp <- subset(data.table(read.csv(infile.ind.anonymised)),select=c('PT_ID', 'AID'))
  dsdt <- merge(dsdt, tmp, by='PT_ID',all.y=T)
  setnames(dsdt, 'VISIT_DT','SAMPLE_DATE')
  dsdt[grepl('^RK-', PT_ID), PANGEA_ID:=paste0('RCCS_', PANGEA_ID)]
  dsdt[grepl('^MRC-', PT_ID), PANGEA_ID:=paste0('MRCUVRI_', PANGEA_ID)]
  dsdt <- merge(dsdt, subset(ds,select=c('RENAME_ID','PANGEA_ID')),by='PANGEA_ID',all.x=T)
  dsdt <- dsdt[!is.na(RENAME_ID)]
  dsdt[,SAMPLE_DATE:=as.Date(SAMPLE_DATE)]
  
  # consensus sampling dates
  consensus_seq_all <- read.fasta(infile.consensus)
  consensus_seq_oneeach <- read.fasta(infile.consensus.oneeach)
  consensus_seq_oneeach_names <- names(consensus_seq_oneeach)
  consensus_seq_oneeach_names <- gsub('\\(.*?\\)', '', consensus_seq_oneeach_names)
  dsdt_con <- data.table(TAXA=names(consensus_seq_all))
  dsdt_con[,RENAME_ID:=gsub('REF_','',TAXA)]
  dsdt_con[,AID:=RENAME_ID]
  dsdt_con[grepl('\\.',AID), SAMPLE_DATE:=sapply(strsplit(AID,'\\.'),
                                                 function(x){
                                                   n <- x[2:3]
                                                   n <- n[!grepl('\\D',n)]
                                                   return(n)
                                                 })]
  id <- which(lengths(strsplit(dsdt_con$SAMPLE_DATE,''))==2)
  dsdt_con[id,SAMPLE_DATE:=ifelse(as.numeric(SAMPLE_DATE)<20,paste0('20',SAMPLE_DATE),paste0('19',SAMPLE_DATE))]
  dsdt_con[!is.na(SAMPLE_DATE),SAMPLE_DATE:=paste(SAMPLE_DATE,7,1,sep='-')]
  dsdt_con[,SAMPLE_DATE:=as.Date(SAMPLE_DATE)]
  dsdt_con[,LOCAL:=!AID %in% consensus_seq_oneeach_names]
  dsdt_con <- dsdt_con[!is.na(SAMPLE_DATE)]
  dsdt_con[LOCAL==T,ST:=sapply(strsplit( AID,'\\.'),function(x)x[1])]
  dsdt_con[LOCAL==F,ST:=gsub('_|\\-','\\.',AID)]
  dsdt_con[LOCAL==F,ST:=sapply(strsplit( ST,split = '\\.'),function(x)x[2])]
  dsdt_con[,unique(ST)]
  dsdt_con[,INCLU_AD:=grepl('A|D',ST)]
  #
  files <-list.files(p)
  files <- grep('treefile$',files,value = T)
  blen <- c()
  age_root <- matrix(nrow = length(files), ncol = 6)
  palA <- viridis_pal(alpha = 1, begin = 0, end = 1, direction = 1, option = "A")
  
  hivc.db.Date2numeric<- function( x )
  {
    if(!class(x)%in%c('Date','character'))	return( x )
    x	<- as.POSIXlt(x)
    tmp	<- x$year + 1900
    x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
    x	
  }
  
  # 1 2 160 161
  for (i in 1:length(files)) {
    # load tree
    cat('processing ',i,'out of ', length(files),'file \n')
    tmptr <- read.tree(file.path(p,files[i]))
    tmpdf <- data.table(TAXA=tmptr$tip.label)
    tmpdf[,RENAME_ID:=TAXA]
    tmpdf[,RENAME_ID:=gsub('CNTRL-|REF_','',RENAME_ID)]
    tmpdf[!grep('REF_',TAXA),RENAME_ID:=gsub('^([A-Z]+[0-9]+-[a-z]+[0-9]+)_.*','\\1',RENAME_ID)]
    tmpdf[,AID:=RENAME_ID]
    tmpdf[!grep('REF_',TAXA),AID:=gsub('^([A-Z]+[0-9]+)-.*','\\1',RENAME_ID)]    
    
    
    # assign color to nodes
    tmpid <- subset(tmpdf,select=c('AID','TAXA'))
    tmp <- unique(subset(tmpid,select=c('AID')))
    tmp[,REF:=!grepl('^AID',AID)]
    tmp[REF==T,color:= rep('#D3D3D3',sum(tmp[,REF]))]
    tmp[REF==F,color:= palA(nrow(tmp)-sum(tmp[,REF]))]
    tmpid <- merge(tmpid, tmp, by='AID',all.x=T)
    
    # add dates
    tmp <- tmpdf[grep('REF_',TAXA),]
    tmpdf <- tmpdf[!grep('REF_',TAXA),]
    tmp <- merge(tmp, subset(dsdt_con,select=c('RENAME_ID','SAMPLE_DATE')),  by='RENAME_ID')
    tmpdf <- merge(tmpdf, subset(dsdt,select=c('RENAME_ID','SAMPLE_DATE')),  by='RENAME_ID',all.x=T)
    tmpdf <- rbind(tmpdf, tmp)
    cat(nrow(tmpdf[is.na(SAMPLE_DATE)]), ' rows has no sample dates \n')
    
    # set sample time 
    dr <- copy(tmpdf)
    dr[grepl('AID',AID),lower:=SAMPLE_DATE]
    dr[grepl('AID',AID),upper:=SAMPLE_DATE]
    dr[!grepl('AID',AID),lower:=as.Date(paste(substr(SAMPLE_DATE,1,4),1,1,sep = '-'))]
    dr[!grepl('AID',AID),upper:=as.Date(paste(substr(SAMPLE_DATE,1,4),12,31,sep = '-'))]
    dr[,lower:=hivc.db.Date2numeric(lower)]
    dr[,upper:=hivc.db.Date2numeric(upper)]
    sts <- tmpdf$SAMPLE_DATE
    sts <- hivc.db.Date2numeric(sts)
    sts <- setNames(sts, tmpdf$TAXA)
    
    # 
    tmpfa <- read.fasta(file.path(p,gsub('treefile$','fasta',files[i])))
    tmplen <- unique(lengths(tmpfa))
    stopifnot(length(tmplen)==1)
    
    # drop tips
    tmptr2 <- drop.tip(tmptr,setdiff(tmptr$tip.label,c(tmpdf$TAXA,'REF_CON_H')))
    sts.df <- subset(dr,select=c('lower','upper'))
    sts.df <- data.frame(sts.df)
    row.names(sts.df) <- dr$TAXA
    stopifnot(sort(row.names(sts.df)) == sort(names(sts)))
    
    
    # root    
    tmptr2 <- root(tmptr2, 'REF_CON_H')
    tmptr2 <- drop.tip(tmptr2,'REF_CON_H')
    
    tmptrd <- dater(tmptr2, sts=sts,s = tmplen, clock='strict',numStartConditions = 0,estimateSampleTimes=sts.df,omega0 = 0.001)
    
    cat('save root to tip regression plot \n')
    pdf(file.path(p,'plots',paste0(gsub('.treefile','_rttreg.pdf',files[i]))),width = 6,height = 4)
    fit <- rootToTipRegressionPlot(tmptrd)
    dev.off()
    
    
    outliers <- outlierTips( tmptrd)
    cat(nrow(outliers[outliers$q < 0.05,]),'outliers \n')
    age_root[i,] <- c(tmptrd$adjusted.mean.rate,
                      tmptrd$mean.rate,
                      tmptrd$timeOfMRCA,
                      tmptrd$timeToMRCA ,
                      tmptrd$Nnode ,
                      tmptrd$coef_of_variation)
    
    # plot
    if(!is.null(downsample)){
      if(downsample==50){
        tmptr <- rescale_tree(tmptr, branch.length = 'rate')
        blen <- c(blen, sum(tmptr$edge.length))
        g <- ggtree(tmptr) %<+% subset(tmpid,select=c('TAXA','color')) +
          geom_tippoint(aes(fill=I(color), color=I(color)))+
          theme_tree2()+
          theme_bw()+
          geom_tiplab(size=1) +
          theme(legend.position = 'none')+
          scale_y_continuous(expand = c(0,5))
        ggsave(file.path(p,'plots',gsub('.treefile','_tree.pdf',files[i])),g,width = 10,height = 0.05*length(tmptr$tip.label),limitsize = F)
        
      }
    }
  }
  
  if(!is.null(downsample)){
    if(downsample==50){
      save(blen,age_root,file=file.path(p,'plots','treeinfo.rda'))
    }else{
      save(age_root,file=file.path(p,'plots','treeinfo.rda'))
    }
  }else{
    save(age_root,file=file.path(p,'plots','treeinfo.rda'))
  }
  
  library(data.table)
  library(ggplot2)
  library(seqinr)
  ds <- c(50,100,200,300,400,500,'_refset1','_refset2','_refset3','')
  ans <- data.table()
  for (i in 1:length(ds)) {
    p <- paste0('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_output',ds[i],'/ptyr1_trees')
    files <-list.files(p)
    files.fa <- grep('fasta$',files,value = T)
    files.fa.id <- gsub(paste0(gsub('trees$','InWindow',basename(p)),'_'),'',files.fa)
    files.fa.id <- as.numeric(gsub("([0-9]+).*$", "\\1", files.fa.id))
    files <- grep('treefile$',files,value = T)
    files.id <- gsub(paste0(gsub('trees$','InWindow',basename(p)),'_'),'',files)
    files.id <- as.numeric(gsub("([0-9]+).*$", "\\1", files.id))
    
    dmaxgap <- data.table()
    for (k in 1:length(files.fa)) {
      fa <- read.fasta(file.path(p,files.fa[k]))
      fa <- do.call(rbind,fa)
      ngap <- apply(fa,2,function(x)sum(x=='-'))
      allgap <- ngap/nrow(fa) >= 0.8
      allgap_ct <- ave(allgap, cumsum(!allgap), FUN = cumsum)
      tmp <- data.table(mgap = max(allgap_ct), start=files.fa.id[k])
      dmaxgap <- rbind(dmaxgap, tmp)
    }
    if (ds[i]==''){
      load(file.path(p,'plots','s4_treeinfo.rda'))
    }else{
      load(file.path(p,'plots','treeinfo.rda'))}
    tmp <- data.table(start=files.id, 
                      subr=age_root[,1],
                      tmrca=age_root[,3],
                      nnode=age_root[,5], 
                      nbg = ds[i])
    cat(nrow(dmaxgap[mgap<=6,])/nrow(dmaxgap), ' windows satisfy the condition \n')
    dmaxgap <- dmaxgap[mgap <=6]
    tmp <- tmp[start %in% dmaxgap$start]
    ans <- rbind(ans, tmp)
  }
  
  ans[nbg=='_refset1',nbg:='set1']
  ans[nbg=='_refset2',nbg:='set2']
  ans[nbg=='_refset3',nbg:='set3']
  ans[nbg=='',nbg:='all']
  
  g <- ggplot(ans[nbg %in% c(50,100,200,300,400,500,'all')],aes(start,tmrca,color=factor(nbg)))+
    geom_smooth()+
    labs(x='\n genome position',y='time to most recent common ancestor \n ',color='background sequences') +
    theme_bw()+
    theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')+
    scale_x_continuous(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05))+guides(color=guide_legend(nrow=1))
  ggsave(file.path(p,'plots','tmrca_ds.pdf'),g,width = 8,height = 6)
  
  g <- ggplot(ans[!nbg %in% c(50,100,200,300,400,500)],aes(start,tmrca,color=factor(nbg)))+
    geom_smooth()+
    labs(x='\n genome position',y='time to most recent common ancestor \n ',color='background sequences') +
    theme_bw()+
    theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')+
    scale_x_continuous(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05))+guides(color=guide_legend(nrow=1))
  ggsave(file.path(p,'plots','tmrca_refset.pdf'),g,width = 8,height = 6)
  
  tmp <- ans[nbg=='all',]
  ans <- ans[nbg!='all',]
  ans <- merge(ans, tmp, by='start')
  g <- ggplot(ans,aes(start,tmrca.y-tmrca.x,color=factor(nbg.x)))+
    geom_smooth()+
    labs(x='\n genome position',y='difference in time to most recent common ancestor \n ',color='reference sequences') +
    theme_bw()+
    theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')+
    scale_x_continuous(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05))+guides(color=guide_legend(nrow=1))
  ggsave(file.path(p,'plots','tmrca_diff_ds.pdf'),g,width = 8,height = 6)
  
}

consensus <- function(){
  infile1 <- '2019_New_ConsensusGenomesOneEach_GeneCut.fasta'
  infile2 <- 'UgandaKenyaTanzaniaGenomes_GeneCut_TreeOrder.FASTA'
  library(seqinr)
  consensus <- c(read.fasta(infile1),read.fasta(infile2))
  names(consensus) <- paste0('REF_',names(consensus))
  consensus  <- lapply(consensus, function(x){x[x!='-']})
  write.fasta(sequences=consensus,names = names(consensus),file.out = '~/consensus.fasta')
}


rkuvri.make.alignments<-function()
{
  #
  #	set up working environment
  require(Phyloscanner.R.utilities)
  library(data.table)
  library(seqinr)
  # file names
  set.seed(42)
  prog.pty <- '~/phyloscanner/phyloscanner_make_trees.py'
  HOME <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'	
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  potential.networks.analysis.dir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100"
  in.dir <- file.path(HOME,'210325_phsc_input')		
  work.dir <- file.path(HOME,"210325_phsc_work")
  out.dir <- file.path(HOME,"210325_phsc_output")	
  package.dir <- file.path(.libPaths(),'Phyloscanner.R.utilities')
  dir.create(in.dir)
  dir.create(work.dir)
  dir.create(out.dir)
  downsample <- 200
  infile.runs <- paste0(HOME,'210120_RCCSUVRI_phscinput_runs.rds')
  infile.consensus <- file.path(in.dir,'ConsensusGenomes.fasta')
  infile.consensus.oneeach <- file.path(in.dir,'2019_New_ConsensusGenomesOneEach_GeneCut.fasta')
  infile.hxb2.package <- file.path(package.dir,'HIV1_compendium_AD_B_CPX_v2.fasta')
  
  # load runs
  pty.runs <- readRDS(infile.runs)
  max.per.run <- 4900
  
  # remove same bam files per run
  setorder(pty.runs,PTY_RUN,-ID_TYPE,UNIT_ID)
  tmp <- pty.runs[,duplicated(SAMPLE_ID),by='PTY_RUN']
  tmp <- tmp[,which(V1)]
  pty.runs <- pty.runs[-tmp,]
  tmp <- pty.runs[,length(SAMPLE_ID)-length(unique(SAMPLE_ID)),by='PTY_RUN']
  stopifnot(all(tmp$V1==0))
  
  # load consensus 
  consensus_seq<- seqinr::read.fasta(infile.consensus)  
  consensus_seq_names <- names(consensus_seq)

  # take hxb2
  hxb2 <- grep('HXB2',names(consensus_seq),value = T)
  hxb2_seq <- consensus_seq[[hxb2]]
  
  # compare hxb2 
  package_seq <- read.fasta(infile.hxb2.package)
  package_hxb2_seq <- package_seq[[1]]
  package_hxb2_seq <- as.character(package_hxb2_seq)[as.character(package_hxb2_seq)!='-']
  hxb2_seq <- as.character(hxb2_seq)[as.character(hxb2_seq)!='-']
  stopifnot(all((hxb2_seq==package_hxb2_seq)))
  # take root
  root.seq <-grep('^REF_CON_M$|REF_CONSENSUS_M$|^REF_CON_H$|REF_CONSENSUS_H$', names(consensus_seq),value = T)
  
  # adapt format
  pty.runs[ID_TYPE=='control',UNIT_ID:=paste0('CNTRL-',UNIT_ID)]
  pty.runs[ID_TYPE=='control',RENAME_ID:=paste0('CNTRL-',RENAME_ID)]
  pty.runs[,BAM:=paste0(data.dir,SAMPLE_ID,'.bam')]
  pty.runs[,REF:=paste0(data.dir,SAMPLE_ID,'_ref.fasta')]
  setkey(pty.runs,PTY_RUN,RENAME_ID)
  # pty.runs <-  pty.runs[1:2,]
  # ptyi <- c(850,900)
  
  #	define phyloscanner input args to generate read alignments 
  #	for each window and each run
  # ptyi <- seq(800,9400,25)		
  ptyi <- seq(800,9175,25) # exclude ends
  ptyi <- c( ptyi[ptyi <= 6615-250],6825,6850,ptyi[ptyi >= 7636]) # exclude vloop except 2 windows without gaps
  # prog.mafft='mafft', 	
  pty.c	<- lapply(seq_along(ptyi), function(i)
  {
    cat('---------------------------------------------------------------- \n')
    print(i)
    pty.args <- list(prog.pty=prog.pty, 
                     prog.mafft='\" mafft --globalpair --maxiterate 1000 \" ',
                     data.dir=data.dir, 
                     work.dir=work.dir, 
                     out.dir=out.dir, 
                     alignments.file=infile.consensus,
                     alignments.root=root.seq,
                     alignments.pairwise.to=hxb2,
                     window.automatic= '', 
                     merge.threshold=0, 
                     min.read.count=1, 
                     quality.trim.ends=23, 
                     min.internal.quality=23, 
                     merge.paired.reads=TRUE, 
                     discard.improper.pairs=TRUE,
                     no.trees=TRUE, 
                     dont.check.duplicates=FALSE,
                     dont.check.recombination=TRUE,
                     num.bootstraps=1,
                     all.bootstrap.trees=TRUE,
                     strip.max.len=350, 
                     min.ureads.individual=NA, 
                     win=c(ptyi[i],ptyi[i]+250,25,250),				 				
                     keep.overhangs=FALSE,
                     mem.save=0,
                     verbose=TRUE,					
                     select=NA,
                     default.coord=TRUE,
                     realignment=TRUE
    )											
    pty.c <- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)
    pty.c[, W_FROM:= ptyi[i]]
    # pty.c[, PTY_RUN:= as.integer(sub('.*ptyr([0-9])_.*','\\1',CMD))]
    pty.c
  })
  pty.c	<- do.call('rbind', pty.c)	
  setkey(pty.c,PTY_RUN,W_FROM)
  pty.c[, CASE_ID:= rep(1:max.per.run,times=ceiling(nrow(pty.c)/max.per.run))[1:nrow(pty.c)]]
  pty.c[, JOB_ID:= rep(1:ceiling(nrow(pty.c)/max.per.run),each=max.per.run)[1:nrow(pty.c)]]
  save(pty.c,file='~/ptyc_210325.rda')
  
  #	define PBS variables
  hpc.load			<- "module load intel-suite/2015.1 mpi raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"	# make third party requirements available	 
  hpc.select			<- 1						# number of nodes
  hpc.nproc			<- 1						# number of processors on node
  hpc.walltime		<- 123						# walltime
  hpc.q				<- "pqeelab"				# PBS queue
  hpc.mem				<- "6gb" 					# RAM	
  hpc.array			<- pty.c[, max(CASE_ID)]	# number of runs for job array
  #	define PBS header for job scheduler. this will depend on your job scheduler.
  pbshead		<- "#!/bin/sh"
  tmp			<- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
  pbshead		<- paste(pbshead, tmp, sep = "\n")
  tmp			<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
  pbshead 	<- paste(pbshead, tmp, sep = "\n")
  pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")	
  if(!is.na(hpc.array))
    pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')	
  if(!is.na(hpc.q)) 
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  pbshead 	<- paste(pbshead, hpc.load, sep = "\n")	
  cat(pbshead)
  
  #	create PBS job array
  for (i in 1:pty.c[, max(JOB_ID)]) {
    tmp<-pty.c[JOB_ID==i,]
    cmd<-tmp[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
    cmd<-cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]			
    cmd<-paste(pbshead,cmd,sep='\n')	
    #	submit job
    outfile<-gsub(':','',paste("readali",paste0('job',i),paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
    outfile<-file.path(work.dir, outfile)
    cat(cmd, file=outfile)
    cmd<-paste("qsub", outfile)
    cat(cmd)
    cat(system(cmd, intern= TRUE))
  }
  
  library(ape)
  # p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_output/ptyr1_trees'
  dir.create(file.path(p,'plots'))
  files <-list.files(p)
  files <- grep('fasta$',files,value = T)
  for (i in 1:length(files)) {
    tmp <- read.dna(file.path(p,files[i]),format = 'fasta')
    pdf(file.path(p,'plots',gsub('.fasta','_alignment.pdf',files[i])),width = 10, height = 10)
    par(cex=0.5)
    checkAlignment(tmp)
    dev.off()
  }
  
}

gap.rm <- function(){
  library(seqinr)
  HOME ='/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  # in.dir	<- file.path(HOME,'210319_phsc_output/ptyr1_trees')
  in.dir	<- file.path(HOME,'210319_phsc_output/ptyr404_trees')
  fa.file <- list.files(in.dir)
  fa.file <- grep('.fasta',fa.file,value = T)
  # out.dir <- file.path(HOME,'210319_phsc_output_gaprmtr1')
  out.dir <- file.path(HOME,'210319_phsc_output_gaprmtr404')
  dir.create(out.dir)
  for (i in 1:length(fa.file)) {
    cat(i,' out of ', length(fa.file),' files \n')
    fa = seqinr::read.fasta(file.path(in.dir,fa.file[i]))
    fa <- do.call(rbind,fa)
    # dim(fa)
    fa_gap <- apply(fa, 2, function(x)sum(x=='-'))
    fa_gap <- fa_gap/nrow(fa)
    fa <- fa[,which(fa_gap < 0.95)]
    # dim(fa)
    write.fasta(split(fa, row(fa)),names = rownames(fa), file.out = file.path(out.dir, fa.file[i]))
  }
  
}

lgap.rm <- function(){
  library(seqinr)
  HOME ='/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  # in.dir	<- file.path(HOME,'210319_phsc_output/ptyr1_trees')
  in.dir	<- file.path(HOME,'210319_phsc_output/ptyr404_trees')
  fa.file <- list.files(in.dir)
  fa.file <- grep('.fasta',fa.file,value = T)
  # out.dir <- file.path(HOME,'210319_phsc_output_lgaprmtr1')
  out.dir <- file.path(HOME,'210319_phsc_output_lgaprmtr404')
  dir.create(out.dir)
  for (i in 1:length(fa.file)) {
    cat(i,' out of ', length(fa.file),' files \n')
    fa = seqinr::read.fasta(file.path(in.dir,fa.file[i]))
    fa <- do.call(rbind,fa)
    # dim(fa)
    fa_gap <- apply(fa, 2, function(x)sum(x=='-'))
    fa_gap <- fa_gap/nrow(fa)
    allgap <- fa_gap > 0.95
    allgap_num <- ave(allgap, cumsum(!allgap), FUN = cumsum)
    rmid <- c()
    while(max(allgap_num) >= 6){
      tmpid <- which.max(allgap_num)[1]
      tmpid <- seq(tmpid-max(allgap_num)+1,tmpid,by=1)
      rmid <- c(rmid,tmpid)
      allgap_num[tmpid] <- 0 
    }
    if(length(rmid)>0){fa <- fa[,-rmid]}
    # dim(fa)
    write.fasta(split(fa, row(fa)),names = rownames(fa), file.out = file.path(out.dir, fa.file[i]))
  }
  
}

rkuvri.make.trees<- function() 
{
  require(data.table)
  require(Phyloscanner.R.utilities)
  #
  #	produce trees
  #
  # lightweight run
  if(0)	
  {
    hpc.select<- 1; hpc.nproc<- 1; hpc.walltime<- 4; hpc.mem<- "1850mb"; hpc.q<- NA
  }
  # midweight run 
  if(1)	
  {
    hpc.select<- 1; hpc.nproc<- 1; hpc.walltime<- 23; hpc.mem<- "1850mb"; hpc.q<- NA
  }
  # heavyweight run 
  if(0)	
  {
    # hpc.select<- 1; hpc.nproc<- 1; hpc.walltime<- 71; hpc.mem<- "5900mb"; hpc.q<- "pqeelab"
    hpc.select<- 1; hpc.nproc<- 1; hpc.walltime<- 71; hpc.mem<- "63850mb"; hpc.q<- NA
  }
  
  HOME <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'	
  max.per.run <- 4900
  iqtree.pr <- 'iqtree'
  iqtree.args			<- ifelse(hpc.nproc==1, '-m GTR+F+R6 -ntmax 1 -seed 42 -o REF_CON_H', 
                          paste0('-m GTR+F+R6 -ntmax ',hpc.nproc,' -seed 42 -o REF_CON_H'))
  in.dir				<- file.path(HOME,'210325_phsc_output')
  out.dir				<- in.dir
  work.dir			<- file.path(HOME,"210325_phsc_work")
  
  infiles	<- data.table(FI=list.files(in.dir, pattern='_v2.fasta$', full.names=TRUE, recursive=TRUE))
  infiles[, FO:= gsub('.fasta$','',FI)]
  infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FI)))]
  infiles[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FI)))]		
  infiles[is.na(W_FROM),W_FROM:= as.integer(gsub('.*PositionsExcised_([0-9]+)_.*','\\1',basename(FI)))]
  setkey(infiles, PTY_RUN, W_FROM)	
  
  # infiles	<- subset(infiles, PTY_RUN<=102) # job1-6
  # infiles	<- subset(infiles, PTY_RUN>102 & PTY_RUN <=204) # job7-12
  # infiles	<- subset(infiles, PTY_RUN>204 & PTY_RUN <=307) # job13-18
  infiles	<- subset(infiles, PTY_RUN>307) # job18-24
  
  df<- infiles[, list(CMD=cmd.iqtree(FI, outfile=FO, pr=iqtree.pr, pr.args=iqtree.args)), by=c('PTY_RUN','W_FROM')]
  df[, ID:=ceiling(seq_len(nrow(df))/4)]
  df<- df[, list(CMD=paste(CMD, collapse='\n',sep='')), by='ID']
  
  #	create PBS job array
  if(nrow(df) > max.per.run){
    df[, CASE_ID:= rep(1:max.per.run,times=ceiling(nrow(df)/max.per.run))[1:nrow(df)]]
    df[, JOB_ID:= rep(1:ceiling(nrow(df)/max.per.run),each=max.per.run)[1:nrow(df)]]
    # set(df, ID, NULL, NULL)
    df[,ID:=NULL]
  }else{
    setnames(df, 'ID', 'CASE_ID')
    df[, JOB_ID:=1]
  }
  
  pbshead	<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=NULL)
  
  
  for (i in 1:df[, max(JOB_ID)]) {
    tmp<-df[JOB_ID==i,]
    hpc.array <- nrow(tmp)
    cmd<-tmp[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
    cmd<-cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]		
    tmp <- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')	
    tmp <- paste(tmp,'module load anaconda3/personal \n source activate phylo' ,sep='\n')
    cmd<-paste(tmp,cmd,sep='\n')	
    #	submit job
    outfile	<- paste("srx",paste0('job',i),paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.')
    outfile<-file.path(work.dir, outfile)
    cat(cmd, file=outfile)
    cmd<-paste("qsub", outfile)
    cat(cmd)
  }
  
  library(tidyverse)
  library(ggtree)
  library(data.table)
  library(seqinr)
  library(ggplot2)
  library(viridis)
  library(treeio)
  library(treedater)
  library(ape)
  library(ggpubr)
  
  # files
  set.seed(42)
  # p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210308_phsc_output/ptyr1_trees'
  # p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210308_phsc_output_gaprmtr1/'
  # p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210308_phsc_output_gaprmtr404/'
  # p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210308_phsc_output/ptyr404_trees'
  # p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210319_phsc_output_lgaprmtr404/'
  dir.create(file.path(p,'plots'))
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  infile.ind.rccs <- file.path('/rds/general/project/ratmann_pangea_deepsequencedata/live/','PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.mrc <- file.path('/rds/general/project/ratmann_pangea_deepsequencedata/live/','PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
  infile.phscinput <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_phscinput_samples.rds'
  infile.consensus <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_input/ConsensusGenomes.fasta'
  infile.consensus.oneeach <- file.path('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_input/2019_New_ConsensusGenomesOneEach_GeneCut.fasta')
  
  # load sampling date
  ds <- readRDS(infile.phscinput)
  dsdt <- data.table(read.csv(infile.ind.rccs))
  dsdt <- subset(dsdt, select=c('pt_id','visit_dt','pangea_id'))
  tmp <- data.table(read.csv(infile.ind.mrc))
  tmp <- subset(tmp, select=c('pt_id','visit_dt','pangea_id'))
  dsdt <- rbind(dsdt, tmp)
  setnames(dsdt, colnames(dsdt), toupper(colnames(dsdt)))
  tmp <- subset(data.table(read.csv(infile.ind.anonymised)),select=c('PT_ID', 'AID'))
  dsdt <- merge(dsdt, tmp, by='PT_ID',all.y=T)
  setnames(dsdt, 'VISIT_DT','SAMPLE_DATE')
  dsdt[grepl('^RK-', PT_ID), PANGEA_ID:=paste0('RCCS_', PANGEA_ID)]
  dsdt[grepl('^MRC-', PT_ID), PANGEA_ID:=paste0('MRCUVRI_', PANGEA_ID)]
  dsdt <- merge(dsdt, subset(ds,select=c('RENAME_ID','PANGEA_ID')),by='PANGEA_ID',all.x=T)
  dsdt <- dsdt[!is.na(RENAME_ID)]
  dsdt[,SAMPLE_DATE:=as.Date(SAMPLE_DATE)]
  
  # consensus sampling dates
  consensus_seq_all <- read.fasta(infile.consensus)
  consensus_seq_oneeach <- read.fasta(infile.consensus.oneeach)
  consensus_seq_oneeach_names <- names(consensus_seq_oneeach)
  consensus_seq_oneeach_names <- gsub('\\(.*?\\)', '', consensus_seq_oneeach_names)
  dsdt_con <- data.table(TAXA=names(consensus_seq_all))
  dsdt_con[,RENAME_ID:=gsub('REF_','',TAXA)]
  dsdt_con[,AID:=RENAME_ID]
  dsdt_con[grepl('\\.',AID), SAMPLE_DATE:=sapply(strsplit(AID,'\\.'),
                                                 function(x){
                                                   n <- x[2:3]
                                                   n <- n[!grepl('\\D',n)]
                                                   return(n)
                                                 })]
  id <- which(lengths(strsplit(dsdt_con$SAMPLE_DATE,''))==2)
  dsdt_con[id,SAMPLE_DATE:=ifelse(as.numeric(SAMPLE_DATE)<20,paste0('20',SAMPLE_DATE),paste0('19',SAMPLE_DATE))]
  dsdt_con[!is.na(SAMPLE_DATE),SAMPLE_DATE:=paste(SAMPLE_DATE,7,1,sep='-')]
  dsdt_con[,SAMPLE_DATE:=as.Date(SAMPLE_DATE)]
  dsdt_con[,LOCAL:=!AID %in% consensus_seq_oneeach_names]
  dsdt_con <- dsdt_con[!is.na(SAMPLE_DATE)]
  dsdt_con[LOCAL==T,ST:=sapply(strsplit( AID,'\\.'),function(x)x[1])]
  dsdt_con[LOCAL==F,ST:=gsub('_|\\-','\\.',AID)]
  dsdt_con[LOCAL==F,ST:=sapply(strsplit( ST,split = '\\.'),function(x)x[2])]
  dsdt_con[,unique(ST)]
  dsdt_con[,INCLU_AD:=grepl('A|D',ST)]
  #
  files <-list.files(p)
  files <- grep('treefile$',files,value = T)
  blen <- c()
  age_root <- matrix(nrow = length(files), ncol = 6)
  palA <- viridis_pal(alpha = 1, begin = 0, end = 1, direction = 1, option = "A")
  
  hivc.db.Date2numeric<- function( x )
  {
    if(!class(x)%in%c('Date','character'))	return( x )
    x	<- as.POSIXlt(x)
    tmp	<- x$year + 1900
    x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
    x	
  }
  
  # 1 2 160 161
  for (i in 1:length(files)) {
    # load tree
    cat('processing ',i,'out of ', length(files),'file \n')
    tmptr <- read.tree(file.path(p,files[i]))
    tmpdf <- data.table(TAXA=tmptr$tip.label)
    tmpdf[,RENAME_ID:=TAXA]
    tmpdf[,RENAME_ID:=gsub('CNTRL-|REF_','',RENAME_ID)]
    tmpdf[!grep('REF_',TAXA),RENAME_ID:=gsub('^([A-Z]+[0-9]+-[a-z]+[0-9]+)_.*','\\1',RENAME_ID)]
    tmpdf[,AID:=RENAME_ID]
    tmpdf[!grep('REF_',TAXA),AID:=gsub('^([A-Z]+[0-9]+)-.*','\\1',RENAME_ID)]    
    
    
    # assign color to nodes
    tmpid <- subset(tmpdf,select=c('AID','TAXA'))
    tmp <- unique(subset(tmpid,select=c('AID')))
    tmp[,REF:=!grepl('^AID',AID)]
    tmp[REF==T,color:= rep('#D3D3D3',sum(tmp[,REF]))]
    tmp[REF==F,color:= palA(nrow(tmp)-sum(tmp[,REF]))]
    tmpid <- merge(tmpid, tmp, by='AID',all.x=T)
    
    # filter
    tmp <- tmpdf[grep('REF_',TAXA),]
    tmpdf <- tmpdf[!grep('REF_',TAXA),]
    tmp <- merge(tmp, subset(dsdt_con,select=c('RENAME_ID','SAMPLE_DATE')),  by='RENAME_ID')
    tmpdf <- merge(tmpdf, subset(dsdt,select=c('RENAME_ID','SAMPLE_DATE')),  by='RENAME_ID',all.x=T)
    tmpdf <- rbind(tmpdf, tmp)
    cat(nrow(tmpdf[is.na(SAMPLE_DATE)]), ' rows has no sample dates \n')
    
    # set sample time 
    dr <- copy(tmpdf)
    dr[grepl('AID',AID),lower:=SAMPLE_DATE]
    dr[grepl('AID',AID),upper:=SAMPLE_DATE]
    dr[!grepl('AID',AID),lower:=as.Date(paste(substr(SAMPLE_DATE,1,4),1,1,sep = '-'))]
    dr[!grepl('AID',AID),upper:=as.Date(paste(substr(SAMPLE_DATE,1,4),12,31,sep = '-'))]
    dr[,lower:=hivc.db.Date2numeric(lower)]
    dr[,upper:=hivc.db.Date2numeric(upper)]
    sts <- tmpdf$SAMPLE_DATE
    sts <- hivc.db.Date2numeric(sts)
    sts <- setNames(sts, tmpdf$TAXA)
    
    # 
    tmpfa <- read.fasta(file.path(p,gsub('treefile$','fasta',files[i])))
    tmplen <- unique(lengths(tmpfa))
    stopifnot(length(tmplen)==1)
    
    # drop tips
    tmptr2 <- drop.tip(tmptr,setdiff(tmptr$tip.label,c(tmpdf$TAXA,'REF_CON_H')))
    sts.df <- subset(dr,select=c('lower','upper'))
    sts.df <- data.frame(sts.df)
    row.names(sts.df) <- dr$TAXA
    stopifnot(sort(row.names(sts.df)) == sort(names(sts)))
    tmptr2 <- root(tmptr2, 'REF_CON_H')
    tmptr2 <- drop.tip(tmptr2,'REF_CON_H')
    
    # dater
    tmptrd <- dater(tmptr2, sts=sts,s = tmplen, clock='strict',numStartConditions = 0,estimateSampleTimes=sts.df,omega0 = 0.001)
    # 
    cat('save root to tip regression plot \n')
    pdf(file.path(p,'plots',paste0(gsub('.treefile','_rttreg.pdf',files[i]))),width = 6,height = 4)
    fit <- rootToTipRegressionPlot(tmptrd)
    dev.off()
    outliers <- outlierTips( tmptrd)
    cat(nrow(outliers[outliers$q < 0.05,]),'outliers \n')
    age_root[i,] <- c(tmptrd$adjusted.mean.rate,
                      tmptrd$mean.rate,
                      tmptrd$timeOfMRCA,
                      tmptrd$timeToMRCA ,
                      tmptrd$Nnode ,
                      tmptrd$coef_of_variation)
    
    # plot
    tmptr <- rescale_tree(tmptr, branch.length = 'rate')
    blen <- c(blen, sum(tmptr$edge.length))
    g <- ggtree(tmptr) %<+% subset(tmpid,select=c('TAXA','color')) +
      geom_tippoint(aes(fill=I(color), color=I(color)))+
      theme_tree2()+
      theme_bw()+
      geom_tiplab(size=1) +
      theme(legend.position = 'none')+
      scale_y_continuous(expand = c(0,5))
    ggsave(file.path(p,'plots',gsub('.treefile','_tree.pdf',files[i])),g,width = 10,height = 0.05*length(tmptr$tip.label),limitsize = F)
  }
  
  save(blen,age_root,file=file.path(p,'plots',paste0('treeinfo.rda')))
  
  p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_output/ptyr1_trees'

  
  files <-list.files(p)
  files <- grep('treefile$',files,value = T)
  ans <- data.table()
  
  load(file.path(p,'plots',paste0('treeinfo.rda')))
  ans <- data.table(start=files, blen = blen,
                    subr=age_root[,1],tmrca=age_root[,3],nnode=age_root[,5])
  setkey(ans,start)
  
  ans[, PositionsExcised:=F]
  ans[grepl('PositionsExcised',files),PositionsExcised:=T]
  ans[,files:=gsub('_PositionsExcised','',files)]
  ans[,files:=  as.numeric(gsub('.*InWindow_([0-9]+)_.*','\\1', files))]
  setkey(ans, PositionsExcised,files)
  tmp <- ans[tmrca<=1920]
  tmp
  
  ptyi <- seq(800,9175,25) # exclude ends
  ptyi <- c( ptyi[ptyi <= 6615-250],6825,6850,ptyi[ptyi >= 7636]) # exclude vloop except 2 windows without gaps
  
  ggplot(ans,aes(files,tmrca,color=PositionsExcised))+
    # geom_smooth()+
    geom_line()+
    theme_bw()+
    labs(x='genomic windows', y='TMRCA',color='Positions excised')+
    theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')+
    scale_x_continuous(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05))+guides(color=guide_legend(nrow=1))
  # ggsave(file.path(p,'plots','tmrca_excised.pdf'),width = 8,height = 6)
  ggsave(file.path(p,'plots','tmrca_excised_raw.pdf'),width = 8,height = 6)
  
  #combine: if position excision exists, use tmrca from it. 
  # ans2 <- ans[,list(tmrca=tail(tmrca,1)),by='files'] 
  ans2 <- ans %>%
    group_by(files) %>%
    slice(n()) %>%
    ungroup()
  
  ans2 <- data.table(ans2)
  ggplot(ans2,aes(files,tmrca,color=PositionsExcised))+
    geom_point()+
    theme_bw()+
    labs(x='genomic windows', y='TMRCA',color='Positions excised')+
    scale_x_continuous(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05))+
    geom_vline(xintercept=c(6615,7636))+guides(color=guide_legend(nrow=1))+
    theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')+
    geom_hline(yintercept = 1920, linetype=2)
  ggsave(file.path(p,'plots','tmrca_combined.pdf'),width = 8,height = 6)
  
  
  # gap analysis
  ans3 <- copy(ans2)
  ans3[,GAP_REMOVE:='no']
  ans2[,GAP_REMOVE:='yes']
  ans4 <- rbind(ans2,ans3)
  ans2[,GAP_REMOVE:='large chunk removed']
  ans4 <- rbind(ans4,ans2)
  ans <- copy(ans4)
  # tmp <- dcast(ans,files~GAP_REMOVE,value.var = 'tmrca')
  # tmp[,est_b4_1920:=F]
  # tmp[no<1920,est_b4_1920:=T]
  tmp[,diff:=yes-no]
  # ggplot(tmp,aes(files,diff,color=est_b4_1920))+
  #   geom_point()+
  #   theme_bw()+
  #   labs(x='genomic windows', y='differences of TMRCA with and without gaps removed',color='estimates before 1920 with gaps')+
  #   scale_x_continuous(expand = c(0,0.05))+
  #   scale_y_continuous(expand = c(0,0.05))+guides(color=guide_legend(nrow=1))+
  #   theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')
  # ggsave(file.path(p,'plots','tmrca_diff_gap.pdf'),width = 8,height = 6)
  
  g1 <- ggplot(ans,aes(files,tmrca,color=factor(GAP_REMOVE)))+
    geom_point()+
    theme_bw()+
    labs(x='genomic windows', y='TMRCA',color='gaps removed')+
    scale_x_continuous(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05))+
    geom_vline(xintercept=c(6615,7636))+guides(color=guide_legend(nrow=1))+
    theme(legend.position = 'bottom',legend.direction = 'horizontal',legend.box = 'horizontal')+
    geom_hline(yintercept = 1920, linetype=2)+
    facet_grid(PositionsExcised~.)
  
  g2 <- ggplot(ans,aes(factor(GAP_REMOVE),tmrca,color=factor(GAP_REMOVE)))+
    geom_boxplot()+
    theme_bw()+
    labs(x='gaps removed', y='')+
    geom_vline(xintercept=c(6615,7636))+guides(color=guide_legend(nrow=1))+
    theme(legend.position = 'none',legend.direction = 'horizontal',legend.box = 'horizontal')+
    geom_hline(yintercept = 1920, linetype=2)+
    facet_grid(PositionsExcised~.)
  
  p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210319_phsc_output/ptyr404_trees'
  ggarrange(g1,g2, ncol = 2, widths = c(0.6,0.4))%>%
    ggexport(filename = file.path(p,'plots','tmrca_combined_gap_compare.pdf'),width = 10,height = 6)
  
  # 
  # 
  # 
  # ggplot(tmp, aes(x='',y=diff)) +
  #   geom_boxplot(outlier.size = 0,width=0.1) +
  #   theme_bw()+
  #   labs(x='', y='differences of TMRCA with and without gaps removed')+
  #   scale_x_discrete(expand = c(0,0.05))+
  #   scale_y_continuous(expand = c(0,0.05))+
  #   geom_hline(yintercept = 0, linetype=2)
  # 
  # ggsave(file.path(p,'plots','tmrca_diff_gap_boxplot.pdf'),width = 2,height = 4)

  
  # compare run 1 and run 404
  ans3 <- copy(ans2)
  ans3[,RUN:=1]
  ans2[,RUN:=404]
  ans <- rbind(ans2,ans3)
  g1 <- ggplot(ans,aes(files,tmrca,color=factor(RUN)))+
    geom_point()+
    theme_bw()+
    labs(x='genomic windows', y='TMRCA',color='run')+
    scale_x_continuous(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05))+
    geom_vline(xintercept=c(6615,7636))+guides(color=guide_legend(nrow=1))+
    theme(legend.position = 'bottom',legend.direction = 'horizontal',legend.box = 'horizontal')+
    geom_hline(yintercept = 1920, linetype=2)+
    facet_grid(PositionsExcised~.)
  
  g2 <- ggplot(ans,aes(factor(RUN),tmrca,color=factor(RUN)))+
    geom_boxplot()+
    theme_bw()+
    labs(x='run', y='TMRCA')+
    geom_vline(xintercept=c(6615,7636))+guides(color=guide_legend(nrow=1))+
    theme(legend.position = 'bottom',legend.direction = 'horizontal',legend.box = 'horizontal')+
    geom_hline(yintercept = 1920, linetype=2)+
    facet_grid(PositionsExcised~.)
   
  p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210308_phsc_output/ptyr1_trees'
  ggarrange(g1,g2, ncol = 2, widths = c(0.8,0.2))%>%
    ggexport(filename = file.path(p,'plots','tmrca_combined_compare.pdf'),width = 8,height = 6)
  

  # ggplot(tmp)+
  #   geom_point(aes(files,yes,color='yes'))+
  #   geom_point(aes(files,no,color='no'))+
  #   theme_bw()+
  #   labs(x='genomic windows', y='TMRCA',color='gaps removed')+
  #   scale_x_continuous(expand = c(0,0.05))+
  #   scale_y_continuous(expand = c(0,0.05))+guides(color=guide_legend(nrow=1))+
  #   theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')
  # ggsave(file.path(p,'plots','tmrca_gap.pdf'),width = 8,height = 6)
  #

  # # no excision
  # g <- ggplot(ans,aes(start,tmrca))+
  #   geom_point()+
  #   geom_line() +
  #   labs(x='\n genome position',y='time to most recent common ancestor \n ') +
  #   theme_bw()+
  #   theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')+
  #   scale_x_continuous(expand = c(0,0.05))+
  #   scale_y_continuous(expand = c(0,0.05))
  # ggsave(file.path(p,'plots','tmrca.pdf'),g,width = 8,height = 6)
  
  
  # number of nodes per window
  g <- ggplot(ans,aes(start,nnode))+
    geom_bar(position="dodge", stat="identity")+
    labs(x='\n genome position',y='number of nodes\n ') +
    theme_bw()+
    theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')+
    scale_x_discrete(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05))
  ggsave(file.path(p,'plots','nnode.pdf'),g,width = 8,height = 6)
  

  
  p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210308_phsc_output/ptyr1_trees'
  # ggarrange(plotlist=g,ncol = 3)%>%
  #   ggexport(filename = file.path(p,'plots','tmrca_combined_job1.pdf') , width = 6*3, height = 4*2)
  ggarrange(plotlist=g,ncol = 3)%>%
    ggexport(filename = file.path(p,'plots','tmrca_combined_job_last.pdf') , width = 6*3, height = 4*2)
  

  
  p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210319_phsc_output/ptyr404_trees'
  # p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210319_phsc_output/ptyr1_trees'
  tmp <- alignment_gap_check(p)
  tmp[posex==T & order(V1),]
  tmp[posex==F & order(V1),]
  # tmp <- alignment_gap_check('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210308_phsc_output_gaprmtr1')
  ggplot(tmp,aes(pos,V1,color=posex)) +
    geom_point() +
    theme_bw() +
    labs(x='genomic positions', y='lengths of the largest gap chunks',color='Positions excised')+
    scale_x_continuous(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05))+
    guides(color=guide_legend(nrow=1))+
    theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')+
    facet_grid(posex~.)
  ggsave(file.path(p,'plots','gap_vs_win.pdf'),width = 8, height = 6)
  
  
  files <-list.files(p)
  files <- grep('treefile$',files,value = T)
  ans <- data.table()
  
  load(file.path(p,'plots',paste0('treeinfo.rda')))
  ans <- data.table(start=files, blen = blen,
                    subr=age_root[,1],tmrca=age_root[,3],nnode=age_root[,5])
  setkey(ans,start)
  
  ans[, PositionsExcised:=F]
  ans[grepl('PositionsExcised',files),PositionsExcised:=T]
  ans[,files:=gsub('_PositionsExcised','',files)]
  ans[,files:=  as.numeric(gsub('.*InWindow_([0-9]+)_.*','\\1', files))]
  setkey(ans, PositionsExcised,files)

  ans <- merge(ans, tmp, by.x=c('PositionsExcised', 'files'), by.y=c('posex', 'pos'))
  ans[V1<6, range(tmrca)]
  ans[V1<6 & tmrca<1920]
  g <- ggplot(ans,aes(V1,tmrca))+
    geom_point()+
    theme_bw()+
    labs(x='\n lengths of largest consecutive gap per window', y = 'TMRCA')
  ggsave(file.path(p,'plots','lgap_tmrca.pdf'), g, width = 6,height = 4)
  
  # setkey(ans,tmrca)
  # setkey(ans,V1)
  setkey(ans, PositionsExcised,files)
  ans2 <- ans %>%
    group_by(files) %>%
    slice(n()) %>%
    ungroup()
  ans2 <- data.table(ans2)
  setkey(ans2,tmrca)
  setkey(ans2,V1)
  
  tmp2 <- dcast(tmp,pos~posex,value.var = 'V1')
  tmp2 <- tmp2[!is.na(`TRUE`)]
  tmp2[`FALSE`!=`TRUE`]
  

  load(file.path(p,'plots',paste0('treeinfo.rda')))
  files <-list.files(p)
  files <- grep('treefile$',files,value = T)
  files_posex <- grepl('PositionsExcised_',files)
  files_start <- vector(mode = 'numeric', length = length(files))
  files_start[which(files_posex==T)] <- as.numeric(gsub('.*PositionsExcised_([0-9]+)_.*', "\\1", files[which(files_posex==T)]))
  files_start[which(files_posex==F)] <- as.numeric(gsub('.*InWindow_([0-9]+)_.*', "\\1", files[which(files_posex==F)]))
  dblen <- data.table(blen=blen, pos=files_start,posex = files_posex)
  tmp2 <- merge(tmp,dblen,by=c('pos','posex'))
  g <- ggplot(tmp2,aes(V1,blen))+
    geom_point()+
    theme_bw()+
    labs(x='\n lengths of largest consecutive gap per window', y = 'branch lengths')
  ggsave(file.path(p,'plots','lgap_blen.pdf'), g, width = 6,height = 4)
  
  # file <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210308_phsc_input/ConsensusGenomes.fasta'	
  # consensus_seq <- read.fasta(file)
  # tmp <-  consensus_seq[[grep('HXB2',names( consensus_seq))]]
  # tmp <- which(tmp=='-')
  # excoord <- c(823,824,825,892,893,894,907,908,909,1012,1013,1014,1156,1157,1158,1384,1385,1386,1444,1445,1446,1930,1931,1932,1957,1958,1959,2014,2015,2016,2023,2024,2025,2080,2081,2082,2134,2135,2136,2191,2192,2193,2280,2281,2282,2283,2284,2285,2298,2299,2300,2310,2311,2312,2316,2317,2318,2319,2320,2321,2322,2323,2324,2340,2341,2342,2346,2347,2348,2349,2350,2351,2352,2353,2354,2355,2356,2357,2358,2359,2360,2373,2374,2375,2379,2380,2381,2385,2386,2387,2388,2389,2390,2391,2392,2393,2394,2395,2396,2400,2401,2402,2409,2410,2411,2412,2413,2414,2415,2416,2417,2424,2425,2426,2430,2431,2432,2436,2437,2438,2439,2440,2441,2442,2443,2444,2457,2458,2459,2460,2461,2462,2463,2464,2465,2469,2470,2471,2472,2473,2474,2478,2479,2480,2481,2482,2483,2496,2497,2498,2499,2500,2501,2502,2503,2504,2505,2506,2507,2514,2515,2516,2517,2518,2519,2520,2521,2522,2526,2527,2528,2529,2530,2531,2535,2536,2537,2670,2671,2672,2679,2680,2681,2703,2704,2705,2709,2710,2711,2733,2734,2735,2742,2743,2744,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2769,2770,2771,2772,2773,2774,2778,2779,2780,2811,2812,2813,2814,2815,2816,2817,2818,2819,2823,2824,2825,2841,2842,2843,2847,2848,2849,2850,2851,2852,2856,2857,2858,2865,2866,2867,2871,2872,2873,2892,2893,2894,2895,2896,2897,2901,2902,2903,2904,2905,2906,2952,2953,2954,2961,2962,2963,3000,3001,3002,3015,3016,3017,3018,3019,3020,3030,3031,3032,3042,3043,3044,3084,3085,3086,3090,3091,3092,3099,3100,3101,3111,3112,3113,3117,3118,3119,3135,3136,3137,3171,3172,3173,3177,3178,3179,3180,3181,3182,3189,3190,3191,3192,3193,3194,3204,3205,3206,3210,3211,3212,3222,3223,3224,3228,3229,3230,3237,3238,3239,3246,3247,3248,3249,3250,3251,3255,3256,3257,3261,3262,3263,3396,3397,3398,3501,3502,3503,3546,3547,3548,3705,3706,3707,4425,4426,4427,4449,4450,4451,4503,4504,4505,4518,4519,4520,4590,4591,4592,4641,4642,4643,4647,4648,4649,4656,4657,4658,4668,4669,4670,4671,4672,4673,4692,4693,4694,4722,4723,4724,4782,4783,4784,4974,4975,4976,5016,5017,5018,5067,5068,5069, seq(6615,6811,by=1), seq(7110,7636,by=1),7863,7864,7865,7866,7867,7868,7869,7870,7871,7872,7873,7874,7875,7876,7877,7881,7882,7883,7884,7885,7886, seq(9400,9719,by=1)) 
  # ptyi <- seq(800,9175,25) # exclude ends
  # ptyi <- c( ptyi[ptyi <= 6615-250],6825,6850,ptyi[ptyi >= 7636]) # exclude vloop except 2 windows without gaps
  # excoord <- excoord[excoord<= 6615-250 | (excoord >= 7650 & excoord <=9424) | (excoord>=6825 & excoord <= 6850+249)]
  # consensus_seq <- do.call(rbind,consensus_seq)
  # consensus_seq <- consensus_seq[,-tmp]
  # tmp <- apply(consensus_seq, 2, function(x){sum(x=='-')})/nrow(consensus_seq)
  # gapcol <- excoord[which(tmp[excoord]>0.95)]
  # tmp[gapcol]
  # excoord[excoord >=1950 & excoord<=1950+249]
}

check.files <- function(dir,p){
  library(data.table)
  dir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210319_phsc_output/'
  p <- 'ptyr1_trees'
  p <- 'ptyr404_trees'
  # p <- 'ptyr409_trees'
  files <-list.files(file.path(dir,p),include.dirs = F)
  excoord <- c(823,824,825,892,893,894,907,908,909,1012,1013,1014,1156,1157,1158,1384,1385,1386,1444,1445,1446,1930,1931,1932,1957,1958,1959,2014,2015,2016,2023,2024,2025,2080,2081,2082,2134,2135,2136,2191,2192,2193,2280,2281,2282,2283,2284,2285,2298,2299,2300,2310,2311,2312,2316,2317,2318,2319,2320,2321,2322,2323,2324,2340,2341,2342,2346,2347,2348,2349,2350,2351,2352,2353,2354,2355,2356,2357,2358,2359,2360,2373,2374,2375,2379,2380,2381,2385,2386,2387,2388,2389,2390,2391,2392,2393,2394,2395,2396,2400,2401,2402,2409,2410,2411,2412,2413,2414,2415,2416,2417,2424,2425,2426,2430,2431,2432,2436,2437,2438,2439,2440,2441,2442,2443,2444,2457,2458,2459,2460,2461,2462,2463,2464,2465,2469,2470,2471,2472,2473,2474,2478,2479,2480,2481,2482,2483,2496,2497,2498,2499,2500,2501,2502,2503,2504,2505,2506,2507,2514,2515,2516,2517,2518,2519,2520,2521,2522,2526,2527,2528,2529,2530,2531,2535,2536,2537,2670,2671,2672,2679,2680,2681,2703,2704,2705,2709,2710,2711,2733,2734,2735,2742,2743,2744,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2769,2770,2771,2772,2773,2774,2778,2779,2780,2811,2812,2813,2814,2815,2816,2817,2818,2819,2823,2824,2825,2841,2842,2843,2847,2848,2849,2850,2851,2852,2856,2857,2858,2865,2866,2867,2871,2872,2873,2892,2893,2894,2895,2896,2897,2901,2902,2903,2904,2905,2906,2952,2953,2954,2961,2962,2963,3000,3001,3002,3015,3016,3017,3018,3019,3020,3030,3031,3032,3042,3043,3044,3084,3085,3086,3090,3091,3092,3099,3100,3101,3111,3112,3113,3117,3118,3119,3135,3136,3137,3171,3172,3173,3177,3178,3179,3180,3181,3182,3189,3190,3191,3192,3193,3194,3204,3205,3206,3210,3211,3212,3222,3223,3224,3228,3229,3230,3237,3238,3239,3246,3247,3248,3249,3250,3251,3255,3256,3257,3261,3262,3263,3396,3397,3398,3501,3502,3503,3546,3547,3548,3705,3706,3707,4425,4426,4427,4449,4450,4451,4503,4504,4505,4518,4519,4520,4590,4591,4592,4641,4642,4643,4647,4648,4649,4656,4657,4658,4668,4669,4670,4671,4672,4673,4692,4693,4694,4722,4723,4724,4782,4783,4784,4974,4975,4976,5016,5017,5018,5067,5068,5069, seq(6615,6811,by=1), seq(7110,7636,by=1),7863,7864,7865,7866,7867,7868,7869,7870,7871,7872,7873,7874,7875,7876,7877,7881,7882,7883,7884,7885,7886, seq(9400,9719,by=1)) 
  ptyi <- seq(800,9175,25) # exclude ends
  ptyi <- c( ptyi[ptyi <= 6615-250],6825,6850,ptyi[ptyi >= 7636]) # exclude vloop except 2 windows without gaps
  tmp <- data.table(start=ptyi, end=ptyi+249, ID=seq_len(length(ptyi)))
  tmp[,ex:=sum(excoord>=start & excoord <= end),by='ID']
  tmp[,PositionsExcised:=ex!=0]
  df <- rbind(tmp[PositionsExcised==T,c('start','PositionsExcised')],data.table(start= ptyi,PositionsExcised=F))
  
  # fasta
  tmp <- data.table(fa = grep('fasta$',files,value=T))
  tmp[grepl('PositionsExcised',fa),PositionsExcised:=T]
  tmp[is.na(PositionsExcised),PositionsExcised:=F]
  tmp[,fa:=gsub('_PositionsExcised','',fa)]
  tmp[, start := as.numeric(gsub(paste0(gsub('trees$','InWindow_',p),'([0-9]+)_.*'),'\\1',fa))]
  df <- merge(df, tmp, by=c('start','PositionsExcised'),all=T)
  
  # tree
  tmp <- data.table(tr = grep('treefile$',files,value=T))
  tmp[grepl('PositionsExcised',tr),PositionsExcised:=T]
  tmp[is.na(PositionsExcised),PositionsExcised:=F]
  tmp[,tr:=gsub('_PositionsExcised','',tr)]
  tmp[, start := as.numeric(gsub(paste0(gsub('trees$','InWindow_',p),'([0-9]+)_.*'),'\\1',tr))]
  df <- merge(df, tmp, by=c('start','PositionsExcised'),all=T)  
  
  df[is.na(fa)] 
  df[!is.na(fa) & is.na(tr)] 
  df[is.na(start)]
}

alignment_gap <- function(){
  
  
  alignment_gap_check <- function(dir){
    library(seqinr)
    library(data.table)
    library(ggplot2)
    files <-list.files(dir)
    files <- grep('ptyr[0-9]+_InWindow_?[A-Za-z]?{16}?_[0-9]+_to_[0-9]+_v2.fasta',files,value = T)
    files_posex <- grepl('PositionsExcised_',files)
    files_start <- vector(mode='numeric',length = length(files))
    files_start[which(files_posex==F)]<- as.numeric(gsub('.*InWindow_([0-9]+)_.*', "\\1", files[which(files_posex==F)]))
    files_start[which(files_posex==T)]<- as.numeric(gsub('.*PositionsExcised_([0-9]+)_.*', "\\1", files[which(files_posex==T)]))
    df <- data.table()
    setwd(dir)
    for (i in 1:length(files)) {
      skip_to_next <- FALSE
      # tryCatch(
      #   {
          fa <- seqinr::read.fasta(files[i])
          fa <- do.call(rbind,fa)
          tmp <- data.table(ngap = apply(fa,2,function(x)sum(x=='-')),
                            ntotal = nrow(fa),
                            pos=files_start[i],
                            posex=files_posex[i],
                            file=files[i])
          df <- rbind(df,tmp)
      #   }, error = function(e) { skip_to_next <<- TRUE}
      # )
      # if(skip_to_next) { next }
    }
    df[,ID:=seq_along(ntotal),by=c('pos','posex')]
    df[,pgap:=ngap/ntotal]
    df[,allgap:=pgap>0.95]
    df[,count_allgap:=ave(allgap, cumsum(!allgap), FUN = cumsum), by=c('pos','posex')]
    
    
    tmp <- df[,max(count_allgap), by=c('pos','posex')]
    setkey(tmp, V1)
    cat(nrow(tmp[V1>=3]), ' out of ', nrow(tmp), ' windows have large gaps (>=3) \n ')
    cat(nrow(tmp[V1>=6]), ' out of ', nrow(tmp), ' windows have large gaps (>=6) \n ')
    return(tmp)
  }
  alignment_gaps <- list()
  output.dir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_output/'
  dirs <- list.files(output.dir)
  for (i in 1:length(dirs)) {
    p <- file.path(output.dir,dirs[i])
    alignment_gaps[[i]] <- alignment_gap_check(p)
  }
  save(alignment_gaps,file = '~/alignment_lgaps.rda')
  
  
# p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_output/ptyr1_trees'
# tmp <- alignment_gap_check(p)
# setkey(tmp,V1)
# tmp[posex==T & order(V1),]
#   tmp[posex==F & order(V1),]
#   # tmp <- alignment_gap_check('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210308_phsc_output_gaprmtr1')
#   ggplot(tmp,aes(pos,V1,color=posex)) +
#     geom_point() +
#     theme_bw() +
#     labs(x='genomic positions', y='lengths of the largest gap chunks',color='Positions excised')+
#     scale_x_continuous(expand = c(0,0.05))+
#     scale_y_continuous(expand = c(0,0.05))+
#     guides(color=guide_legend(nrow=1))+
#     theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')+
#     facet_grid(posex~.)
#   ggsave(file.path(p,'plots','gap_vs_win.pdf'),width = 8, height = 6)
#   

}

alignment_tmrca <- function(){
  
  tmrca_check <- function(p){
    library(tidyverse)
    library(ggtree)
    library(data.table)
    library(seqinr)
    library(ggplot2)
    library(viridis)
    library(treeio)
    library(treedater)
    library(ape)
    library(ggpubr)
    
    # files
    set.seed(42)
    infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
    infile.ind.rccs <- file.path('/rds/general/project/ratmann_pangea_deepsequencedata/live/','PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
    infile.ind.mrc <- file.path('/rds/general/project/ratmann_pangea_deepsequencedata/live/','PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
    infile.phscinput <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_phscinput_samples.rds'
    infile.consensus <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_input/ConsensusGenomes.fasta'
    infile.consensus.oneeach <- file.path('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_input/2019_New_ConsensusGenomesOneEach_GeneCut.fasta')
    
    # load sampling date
    ds <- readRDS(infile.phscinput)
    dsdt <- data.table(read.csv(infile.ind.rccs))
    dsdt <- subset(dsdt, select=c('pt_id','visit_dt','pangea_id'))
    tmp <- data.table(read.csv(infile.ind.mrc))
    tmp <- subset(tmp, select=c('pt_id','visit_dt','pangea_id'))
    dsdt <- rbind(dsdt, tmp)
    setnames(dsdt, colnames(dsdt), toupper(colnames(dsdt)))
    tmp <- subset(data.table(read.csv(infile.ind.anonymised)),select=c('PT_ID', 'AID'))
    dsdt <- merge(dsdt, tmp, by='PT_ID',all.y=T)
    setnames(dsdt, 'VISIT_DT','SAMPLE_DATE')
    dsdt[grepl('^RK-', PT_ID), PANGEA_ID:=paste0('RCCS_', PANGEA_ID)]
    dsdt[grepl('^MRC-', PT_ID), PANGEA_ID:=paste0('MRCUVRI_', PANGEA_ID)]
    dsdt <- merge(dsdt, subset(ds,select=c('RENAME_ID','PANGEA_ID')),by='PANGEA_ID',all.x=T)
    dsdt <- dsdt[!is.na(RENAME_ID)]
    dsdt[,SAMPLE_DATE:=as.Date(SAMPLE_DATE)]
    
    # consensus sampling dates
    consensus_seq_all <- read.fasta(infile.consensus)
    consensus_seq_oneeach <- read.fasta(infile.consensus.oneeach)
    consensus_seq_oneeach_names <- names(consensus_seq_oneeach)
    consensus_seq_oneeach_names <- gsub('\\(.*?\\)', '', consensus_seq_oneeach_names)
    dsdt_con <- data.table(TAXA=names(consensus_seq_all))
    dsdt_con[,RENAME_ID:=gsub('REF_','',TAXA)]
    dsdt_con[,AID:=RENAME_ID]
    dsdt_con[grepl('\\.',AID), SAMPLE_DATE:=sapply(strsplit(AID,'\\.'),
                                                   function(x){
                                                     n <- x[2:3]
                                                     n <- n[!grepl('\\D',n)]
                                                     return(n)
                                                   })]
    id <- which(lengths(strsplit(dsdt_con$SAMPLE_DATE,''))==2)
    dsdt_con[id,SAMPLE_DATE:=ifelse(as.numeric(SAMPLE_DATE)<20,paste0('20',SAMPLE_DATE),paste0('19',SAMPLE_DATE))]
    dsdt_con[!is.na(SAMPLE_DATE),SAMPLE_DATE:=paste(SAMPLE_DATE,7,1,sep='-')]
    dsdt_con[,SAMPLE_DATE:=as.Date(SAMPLE_DATE)]
    dsdt_con[,LOCAL:=!AID %in% consensus_seq_oneeach_names]
    dsdt_con <- dsdt_con[!is.na(SAMPLE_DATE)]
    dsdt_con[LOCAL==T,ST:=sapply(strsplit( AID,'\\.'),function(x)x[1])]
    dsdt_con[LOCAL==F,ST:=gsub('_|\\-','\\.',AID)]
    dsdt_con[LOCAL==F,ST:=sapply(strsplit( ST,split = '\\.'),function(x)x[2])]
    dsdt_con[,unique(ST)]
    dsdt_con[,INCLU_AD:=grepl('A|D',ST)]
    #
    files <-list.files(p)
    files <- grep('treefile$',files,value = T)
    blen <- c()
    age_root <- matrix(nrow = length(files), ncol = 6)
    palA <- viridis_pal(alpha = 1, begin = 0, end = 1, direction = 1, option = "A")
    
    hivc.db.Date2numeric<- function( x )
    {
      if(!class(x)%in%c('Date','character'))	return( x )
      x	<- as.POSIXlt(x)
      tmp	<- x$year + 1900
      x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
      x	
    }
    

    for (i in 1:length(files)) {
      # load tree
      cat('processing ',i,'out of ', length(files),'file \n')
      tmptr <- read.tree(file.path(p,files[i]))
      tmpdf <- data.table(TAXA=tmptr$tip.label)
      tmpdf[,RENAME_ID:=TAXA]
      tmpdf[,RENAME_ID:=gsub('CNTRL-|REF_','',RENAME_ID)]
      tmpdf[!grep('REF_',TAXA),RENAME_ID:=gsub('^([A-Z]+[0-9]+-[a-z]+[0-9]+)_.*','\\1',RENAME_ID)]
      tmpdf[,AID:=RENAME_ID]
      tmpdf[!grep('REF_',TAXA),AID:=gsub('^([A-Z]+[0-9]+)-.*','\\1',RENAME_ID)]    
      
      
      # assign color to nodes
      tmpid <- subset(tmpdf,select=c('AID','TAXA'))
      tmp <- unique(subset(tmpid,select=c('AID')))
      tmp[,REF:=!grepl('^AID',AID)]
      tmp[REF==T,color:= rep('#D3D3D3',sum(tmp[,REF]))]
      tmp[REF==F,color:= palA(nrow(tmp)-sum(tmp[,REF]))]
      tmpid <- merge(tmpid, tmp, by='AID',all.x=T)
      
      # filter
      tmp <- tmpdf[grep('REF_',TAXA),]
      tmpdf <- tmpdf[!grep('REF_',TAXA),]
      tmp <- merge(tmp, subset(dsdt_con,select=c('RENAME_ID','SAMPLE_DATE')),  by='RENAME_ID')
      tmpdf <- merge(tmpdf, subset(dsdt,select=c('RENAME_ID','SAMPLE_DATE')),  by='RENAME_ID',all.x=T)
      tmpdf <- rbind(tmpdf, tmp)
      cat(nrow(tmpdf[is.na(SAMPLE_DATE)]), ' rows has no sample dates \n')
      
      # set sample time 
      dr <- copy(tmpdf)
      dr[grepl('AID',AID),lower:=SAMPLE_DATE]
      dr[grepl('AID',AID),upper:=SAMPLE_DATE]
      dr[!grepl('AID',AID),lower:=as.Date(paste(substr(SAMPLE_DATE,1,4),1,1,sep = '-'))]
      dr[!grepl('AID',AID),upper:=as.Date(paste(substr(SAMPLE_DATE,1,4),12,31,sep = '-'))]
      dr[,lower:=hivc.db.Date2numeric(lower)]
      dr[,upper:=hivc.db.Date2numeric(upper)]
      sts <- tmpdf$SAMPLE_DATE
      sts <- hivc.db.Date2numeric(sts)
      sts <- setNames(sts, tmpdf$TAXA)
      
      # 
      tmpfa <- read.fasta(file.path(p,gsub('treefile$','fasta',files[i])))
      tmplen <- unique(lengths(tmpfa))
      stopifnot(length(tmplen)==1)
      
      # drop tips
      tmptr2 <- drop.tip(tmptr,setdiff(tmptr$tip.label,c(tmpdf$TAXA,'REF_CON_H')))
      sts.df <- subset(dr,select=c('lower','upper'))
      sts.df <- data.frame(sts.df)
      row.names(sts.df) <- dr$TAXA
      stopifnot(sort(row.names(sts.df)) == sort(names(sts)))
      tmptr2 <- root(tmptr2, 'REF_CON_H')
      tmptr2 <- drop.tip(tmptr2,'REF_CON_H')
      
      # dater
      tmptrd <- treedater:::dater(tmptr2, sts=sts,s = tmplen, clock='strict',numStartConditions = 0,estimateSampleTimes=sts.df,omega0 = 0.001)
      outliers <- outlierTips( tmptrd)
      cat(nrow(outliers[outliers$q < 0.05,]),'outliers \n')
      age_root[i,] <- c(tmptrd$adjusted.mean.rate,
                        tmptrd$mean.rate,
                        tmptrd$timeOfMRCA,
                        tmptrd$timeToMRCA ,
                        tmptrd$Nnode ,
                        tmptrd$coef_of_variation)
      
      tmptr <- rescale_tree(tmptr, branch.length = 'rate')
      blen <- c(blen, sum(tmptr$edge.length))
    }
    return(list(blen,age_root))
  }

  args <- commandArgs(trailingOnly = TRUE)
  print(args)
  num <- as.numeric(args[1])
  output.dir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_output/'
  dirs <- sort(list.files(output.dir,full.names = T))
  blen_list <- list()
  age_root_list <- list()
  tmp <- tmrca_check(dirs[num])
  save(tmp,file = file.path(dirs[num],'alignment_tmrcas.rda'))
    # dirs <- list.files(output.dir,full.names = T)
  # blen_list <- list()
  # age_root_list <- list()
  # for (i in 1:length(dirs)) {
  #   tmp <- tmrca_check(dirs[i])
  #   blen_list[[i]] <- tmp$blen
  #   age_root_list[[i]] <- tmp$age_root
  # }
  # save(blen_list,age_root_list,file = '~/alignment_tmrcas.rda')
  
  # script
  cmd <-" #!/bin/sh \n #PBS -l walltime=71:59:00 \n #PBS -l select=1:ncpus=1:ompthreads=1:mem=20gb \n #PBS -j oe \n #PBS -J 1-409 \n #PBS -q pqeelab\n module load anaconda3/personal \n source activate phylo \n cd /rds/general/user/xx4515/home/ \n case $PBS_ARRAY_INDEX in"
  tmp <- c()
  for (i in 1:500) {
    tmp <- paste0(tmp,i,')\n Rscript 210325_tmrcas_check_batch.R ',i, '\n ;; \n')
  }
  cmd <- paste0(cmd,' \n ', tmp)
  cmd <- paste0(cmd, 'esac \n ')
  #	submit job
  outfile		<- '~/210325_tmrcas_check_batch.qsub'
  cat(cmd, file=outfile)
}
library(data.table)
# load('~/alignment_tmrcas.rda')
load('~/alignment_lgaps.rda')
alignment_gaps <- rbindlist(alignment_gaps,idcol = T)
alignment_gaps <- dcast(alignment_gaps, .id + pos ~ posex, value.var='V1')
tmp <- alignment_gaps[is.na(`TRUE`)]
setnames(tmp,'FALSE','V1')
tmp[,`TRUE`:=NULL]
tmp[,posex:=FALSE]
alignment_gaps <- alignment_gaps[!is.na(`TRUE`)]
setnames(alignment_gaps,'TRUE','V1')
alignment_gaps[,`FALSE`:=NULL]
alignment_gaps[,posex:=TRUE]
alignment_gaps <- rbind(alignment_gaps,tmp)
alignment_gaps <- dcast(alignment_gaps, .id~pos, value.var='V1')
tmp <- as.matrix(alignment_gaps)
# apply(tmp, 2, function(x)paste0(round(median(x,na.rm = T),2), '[', round(quantile(x,0.025,na.rm = T)), '-', round(quantile(x,0.975,na.rm = T)), ']'))
hist(apply(tmp[,2:ncol(tmp)], 1, function(x)max(x,na.rm = T)))
# [1] 12 21 22 23 28 30 21 15 21 21 30 27 24 30 33 27 21 23 21 21 21 26 24 21 42
# [26] 21 21 30 30 24 21 21 21 23 27 21 21 21 17 24 21 19 12 21 21 24 21 12 21 15
# [51] 28 36 24 12 36 28 33 25 33 30 21 36 21 21 30 10 21 33 27 24 24 30 23 21 21
# [76] 45 24 33 42 27 21 21 21 27 21 17 21 21 35 27 21  9 21 25 33 21 30 21 34 33
# [101] 22 24 21 26 21 22 21 21 21 16 21 21 18 21 26 12 25 24 28 22 13 21 30 24 15
# [126] 24 21 33 33 21 36 12 15 21 24 27 63 31 21 36 21 24 21 23 30 39 24 24 27 21
# [151] 21 25 15 18 23 11 21 24 24 22 18 21 21 18 20 21 21 42 21 33 15 12 21 30 30
# [176] 12 21 22 15 21 23 23 15 30  9 12 15 22 35 23 15 21 21 27 25 30 21 20 21 21
# [201] 21 21 24 24 24 24 24 24 19 24 39 24 21 23 30 45 27 10 36 45 19 27 29 33  6
# [226] 12 18 21 33  9 21 24 31 30 42 21 30 21  8 42 20 30 12  9 10 45 15 33  3 18
# [251] 25 21 21 12 24 14 30 24 16 21 15 13 33 21 24 39 24 27 19 24  9 30 12 21 30
# [276] 21 18 11 27 21 20 21 10 21 27 20 29 25 22 18 24 21 12 30 24 27 23 20 30  7
# [301] 21 27 27 33 21 34 30 24 21 30 21 25 21 21 22 16 15 21 36 24 21 19 24 18 18
# [326] 31 25 34 21 27 25 31 36 24 21 34 21 23 21 33 27 27 36 27 30 19 24 27 33 23
# [351] 18 16 21 21 36 22 30 21 30 21 16 42 24 30 25 46 21 75 19 21 36 33 21 21 21
# [376] 22 36 21 32 21 27 30 30 21 15 45 24 21 39 39 22 27 21 21 31 21 24 30 19 25
# [401] 26 15 16 21 27 21 30 16 24


delete_nonexcised <- function(){
  require(data.table)
  
  # files
  HOME <<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  treedir <- file.path(HOME,'210325_phsc_output/')
  df <- data.table(F=list.files(treedir,recursive = T, pattern = '*.treefile'))
  # which to delete
  df[,Ex:=grepl('PositionsExcised',F)]
  df[,base:=basename(F)]
  df[,dir:=dirname(F)]
  df[,wfrom:=sub("^ptyr[0-9]+_InWindow_[A-Za-z]?{16}?_?([0-9]+)_to_([0-9]+)_v2\\.treefile$","\\1",base)]
  df[,wto:=sub("^ptyr[0-9]+_InWindow_[A-Za-z]?{16}?_?([0-9]+)_to_([0-9]+)_v2\\.treefile$","\\2",base)]
  df[,run:=sub("ptyr([0-9]+)_trees","\\1",dir)]
  tmp <- dcast(df, run + wfrom + wto ~ Ex, value.var="F")
  tmp <- data.table(tmp)
  tmp <- tmp[!is.na(`TRUE`) & !is.na(`FALSE`)]
  tmp2 <- unique(subset(df,select='dir'))
  tmp2[,dir:=file.path(treedir,dir)]
  # delete
  tmp[,dir:=dirname(`TRUE`)]
  for (i in 1:nrow(tmp)){
    cat(system(paste0( 'rm ',file.path(treedir,tmp$`FALSE`[i])),intern = T))
  }
  # don't generate trees for those files next time
}
make.phyloscanner<- function()
{
  require(tidyverse)
  require(data.table)
  require(phyloscannerR)
  
  HOME <<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  treedir <- file.path(HOME,'210325_phsc_output/')
  tmpdir	<- file.path(HOME,"210325_phsc_work/")
  # outdir	<- file.path(HOME,'210325_phsc_phscrelationships_02_05/')
  # sa1
  outdir	<- file.path(HOME,'210325_phsc_phscrelationships_02_05_null_min_read/')
  dir.create(outdir)
  prog.phyloscanner_analyse_trees <- '/rds/general/user/xx4515/home/phyloscanner/phyloscanner_analyse_trees.R'

  #	set phyloscanner variables	
  control	<- list()
  control$allow.mt <- TRUE				
  control$alignment.file.directory = NULL 
  control$alignment.file.regex = NULL
  control$blacklist.underrepresented = FALSE	
  control$count.reads.in.parsimony = TRUE
  control$distance.threshold <- '0.02 0.05'
  control$do.dual.blacklisting = FALSE					
  control$duplicate.file.directory = NULL
  control$duplicate.file.regex = NULL
  control$file.name.regex = "^(?:.*\\D)?([0-9]+)_to_([0-9]+).*$"
  control$guess.multifurcation.threshold = FALSE
  control$max.reads.per.host <- NULL
  # control$min.reads.per.host <- 30 
  # sa1
  control$min.reads.per.host <- NULL
  control$min.tips.per.host <- 1	
  control$multifurcation.threshold = 1e-5
  control$multinomial= TRUE
  control$norm.constants = NULL
  control$norm.ref.file.name = "~/normalisation_ByPosition.csv"
  control$norm.standardise.gag.pol = TRUE
  control$no.progress.bars = TRUE
  control$outgroup.name = "REF_CON_H"
  control$output.dir = outdir
  control$parsimony.blacklist.k = 15
  control$prune.blacklist = FALSE
  control$post.hoc.count.blacklisting <- TRUE
  control$ratio.blacklist.threshold = 0.01
  control$raw.blacklist.threshold = 3			
  control$recombination.file.directory = NULL
  control$recombination.file.regex = NULL
  control$relaxed.ancestry = TRUE
  control$sankoff.k = 15
  control$sankoff.unassigned.switch.threshold = 0
  control$seed = 42
  control$splits.rule = 's'
  control$tip.regex = '^(.*)-fq[0-9]+_read_([0-9]+)_count_([0-9]+)'
  control$tree.file.regex = "^(.*)\\.treefile$" # from rscript
  control$treeFileExtension = '.treefile'
  control$use.ff = FALSE
  control$user.blacklist.directory = NULL 
  control$user.blacklist.file.regex = NULL
  control$verbosity = 1	
  

  #	make bash for many files	
  df <- tibble(F=list.files(treedir))
  # ,pattern = '*.treefile$')
  df <- df %>%
    mutate(TYPE:= gsub('ptyr([0-9]+)_(.*)','\\2', F),
           RUN:= as.integer(gsub('ptyr([0-9]+)_(.*)','\\1', F))) %>%
    # mutate(TYPE:= gsub('^([^\\.]+)\\.[a-z]+$','\\1',TYPE)) %>%
    mutate(TYPE:= gsub('^[^\\.]+\\.([a-z]+)$','\\1',TYPE)) %>%
    spread(TYPE, F) %>%
    set_names(~ str_to_upper(.))

  valid.input.args <- cmd.phyloscanner.analyse.trees.valid.args(prog.phyloscanner_analyse_trees)
  cmds <- vector('list',nrow(df))
  
  for(i in seq_len(nrow(df)))
  {
    #	set input args
    control$output.string <- paste0('ptyr',df$RUN[i])
    #	make script
    tree.input <- file.path(treedir, df$TREES[i])
    cmd <- cmd.phyloscanner.analyse.trees(prog.phyloscanner_analyse_trees, 
                                          tree.input, 
                                          control,
                                          valid.input.args=valid.input.args)
    cmds[[i]] <- cmd		
  }	
  cat(cmds[[1]])
  
  #
  # 	submit array job to HPC
  #
  #	make header
  hpc.load			<- "module load anaconda3/personal \n source activate phylo"	# make third party requirements available	 
  hpc.select			<- 1						# number of nodes
  hpc.nproc			<- 1						# number of processors on node
  hpc.walltime		<- 923						# walltime
  if(0)		
  {
    hpc.q			<- NA						# PBS queue
    hpc.mem			<- "36gb" 					# RAM		
  }
  #		or run this block to submit a job array to Oliver's machines
  if(1)
  {
    hpc.q			<- "pqeelab"				# PBS queue
    hpc.mem			<- "12gb" 					# RAM		
  }
  hpc.array			<- length(cmds)	# number of runs for job array	
  pbshead		<- "#!/bin/sh"
  tmp			<- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
  pbshead		<- paste(pbshead, tmp, sep = "\n")
  tmp			<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
  pbshead 	<- paste(pbshead, tmp, sep = "\n")
  pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")	
  if(!is.na(hpc.array))
    pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')	
  if(!is.na(hpc.q)) 
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  pbshead 	<- paste(pbshead, hpc.load, sep = "\n")	
  cat(pbshead)
  #	make array job
  for(i in 1:length(cmds))
    cmds[[i]]<- paste0(i,')\n',cmds[[i]],';;\n')
  cmd		<- paste0('case $PBS_ARRAY_INDEX in\n',paste0(cmds, collapse=''),'esac')	
  cmd		<- paste(pbshead,cmd,sep='\n')	
  #	submit job
  outfile		<- gsub(':','',paste("phsc",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
  outfile		<- file.path(tmpdir, outfile)
  cat(cmd, file=outfile)
  cmd 		<- paste("qsub", outfile)
  cat(cmd)
  cat(system(cmd, intern= TRUE))
  
}
phsc.migrate.transmission.networks<- function()
{
  library(data.table)
  source('~/transmission_network_functions_phyloscanner.R')
  # optional: meta data
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05/'
  # sa1
  # indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/'
  infiles	<- data.table(F=list.files(indir, pattern='*workspace.rda$', full.names=TRUE))
  infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(F)))]
  setkey(infiles, PTY_RUN)
  
  #	prepare all dwin and dc
  dca	<- infiles[, {
    cat(PTY_RUN,'\n')
    load(F)
    dc			
  }, by='PTY_RUN']
  dwina <- infiles[, {
    cat(PTY_RUN,'\n')
    load(F)
    dwin			
  }, by='PTY_RUN']
  control <- list(linked.group='close.and.adjacent.cat',linked.no='not.close.or.nonadjacent',linked.yes='close.and.adjacent', conf.cut=0.6, neff.cut=3,weight.complex.or.no.ancestry=0.5)
  # save(dca,dwina,file = '~/dcdwina_minread30.rda')
  # save(dca,dwina,file = '~/dcdwina_minreadnull.rda')
  load('~/dcdwina_minread30.rda')
  # load('~/dcdwina_minreadnull.rda')
  dca[,CNTRL1:=FALSE]
  dca[,CNTRL2:=FALSE]
  dca[grepl('CNTRL-',host.1),CNTRL1:=TRUE]
  dca[grepl('CNTRL-',host.2),CNTRL2:=TRUE]
  dca[grepl('CNTRL-',host.1),host.1:=gsub('CNTRL-','',host.1)]
  dca[grepl('CNTRL-',host.2),host.2:=gsub('CNTRL-','',host.2)]
  dwina[,CNTRL1:=FALSE]
  dwina[,CNTRL2:=FALSE]
  dwina[grepl('CNTRL-',host.1),CNTRL1:=TRUE]
  dwina[grepl('CNTRL-',host.2),CNTRL2:=TRUE]
  dwina[grepl('CNTRL-',host.1),host.1:=gsub('CNTRL-','',host.1)]
  dwina[grepl('CNTRL-',host.2),host.2:=gsub('CNTRL-','',host.2)]
  
  
  tmp <- dca[!(host.1<host.2)]
  setnames(tmp,c('host.1','host.2','CNTRL1','CNTRL2'),c('host.2','host.1','CNTRL2','CNTRL1'))
  dca <- rbind(dca[(host.1<host.2)],tmp)
  tmp <- dwina[!(host.1<host.2)]
  setnames(tmp,c('host.1','host.2','CNTRL1','CNTRL2','paths12','paths21','nodes1','nodes2'),c('host.2','host.1','CNTRL2','CNTRL1','paths21','paths12','nodes2','nodes1'))
  tmp[,close.and.contiguous.and.directed.cat:=ifelse(close.and.contiguous.and.directed.cat=='12'|close.and.contiguous.and.directed.cat=='21',ifelse(close.and.contiguous.and.directed.cat=='12','21','12'),close.and.contiguous.and.directed.cat)]
  tmp[,close.and.contiguous.and.ancestry.cat:=ifelse(close.and.contiguous.and.ancestry.cat=='12'|close.and.contiguous.and.ancestry.cat=='21',ifelse(close.and.contiguous.and.ancestry.cat=='12','21','12'),close.and.contiguous.and.ancestry.cat)]
  tmp[, close.and.adjacent.and.ancestry.cat:=ifelse( close.and.adjacent.and.ancestry.cat=='12'| close.and.adjacent.and.ancestry.cat=='21',ifelse( close.and.adjacent.and.ancestry.cat=='12','21','12'), close.and.adjacent.and.ancestry.cat)]
  dwina <- rbind(dwina[(host.1<host.2)],tmp)
  
  # tmp <- find.pairs.in.networks(indir, batch.regex='^ptyr([0-9]+)_.*', control=control, verbose=TRUE, dmeta=dmeta)
  library(tidyverse)
  # dwin <- copy(dwina)
  # dc <- copy(dca)
  # verbose=TRUE
  # dmeta=NULL
  tmp <- find.pairs.in.networks(dwina, dca, control=control, verbose=TRUE)
  # dpl <- copy(tmp$linked.pairs)
  dpl <- copy(tmp$network.pairs)
  dc <- copy(tmp$relationship.counts)
  dw <- copy(tmp$windows)
  save(dpl, dc, dw, file=file.path(indir,'Rakai_phscnetworks_allpairs.rda'))
  
  library(glue)
  library(igraph)
  library(RBGL)
  tmp <- find.networks(dc, control=control, verbose=TRUE)
  dnet <- copy(tmp$transmission.networks)
  dchain <- copy(tmp$most.likely.transmission.chains)
  save(dpl, dc, dw, dnet, dchain, file=file.path(indir,'Rakai_phscnetworks.rda'))
  
}


analysis_network <- function(){
  library(data.table)
  # dirs
  # indir			<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05/'
  # # sa1
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/'
  # infile.pairs	<- file.path(indir,'Cluster_strongsupport_allwindows.rda')
  # infile.networks	<- file.path(indir,'Cluster_strongsupport_networksallpairs.rda')
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  infile.ind.rccs <- file.path('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.mrc <- file.path('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
  infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
  load( infile.networks )	# loads rtp, rplkl, rpw

  
  # make table mapping id and aid, aid and genders
  df <- data.table(read.csv(infile.ind.anonymised))
  tmp <- rbind(data.table(read.csv(infile.ind.rccs)),data.table(read.csv(infile.ind.mrc)))
  df2 <- merge(subset(df,select=c('PT_ID', 'AID')),subset(tmp,select=c('pt_id','sex')),by.x='PT_ID',by.y='pt_id',all.x=T)
  df2 <- subset(df2,select=c('AID','sex'))
  colnames(df2) <- c('ID','SEX')
  df2 <- unique(df2)
  
  # classification
  conf.cut		<- 0.6			
  dchain <- as.data.table(dchain)
  dchain[, PHYLOSCANNER_CLASSIFY:= NA_character_]
  set(dchain, dchain[, which(is.na(PTY_RUN))], 'PHYLOSCANNER_CLASSIFY', 'insufficient deep sequence data for at least one individual')
  set(dchain, dchain[, which(!is.na(PTY_RUN) & is.na(LINK_12) & is.na(LINK_21))], 'PHYLOSCANNER_CLASSIFY', 'ph unlinked pair')
  set(dchain, dchain[, which(!is.na(SCORE_LINKED) & SCORE_LINKED<=conf.cut)], 'PHYLOSCANNER_CLASSIFY', 'unclear if pair ph linked or unlinked')
  set(dchain, dchain[, which(!is.na(SCORE_LINKED) & SCORE_LINKED>conf.cut)], 'PHYLOSCANNER_CLASSIFY', 'ph linked pair direction not resolved')
  set(dchain, dchain[, which(!is.na(SCORE_LINKED) & SCORE_LINKED>conf.cut & SCORE_DIR_12>conf.cut)], 'PHYLOSCANNER_CLASSIFY', 'ph linked pair direction 12')
  set(dchain, dchain[, which(!is.na(SCORE_LINKED) & SCORE_LINKED>conf.cut & SCORE_DIR_21>conf.cut)], 'PHYLOSCANNER_CLASSIFY', 'ph linked pair direction 21')
  
  # if several classifications per pair, take the majority
  tmp <- dchain[,length(unique( PHYLOSCANNER_CLASSIFY)),by=c('H1','H2')]
  tmp <- tmp[V1>1]
  tmp[, H:=paste0(H1,'-',H2)]
  tmpm <- tmp$H
  tmp2 <- dchain[paste0(H1,'-',H2)%in% tmpm]
  tmp <- tmp2[,length(SCORE_LINKED),by=c('H1','H2','PHYLOSCANNER_CLASSIFY')]
  setkey(tmp,H1,H2,V1)
  tmp <- tmp[,list(PHYLOSCANNER_CLASSIFY=last(PHYLOSCANNER_CLASSIFY)),by=c('H1','H2')]
  tmp2 <- merge(tmp2, tmp, by=c('H1','H2','PHYLOSCANNER_CLASSIFY'))
  dchain <- dchain[!(paste0(H1,'-',H2)%in% tmpm)]
  dchain <- rbind(dchain,tmp2)
  # dchain[H1=='AID2324'& H2== 'AID2368']
  # 
  dchain <- unique(subset(dchain,select = c('H1','H2','PHYLOSCANNER_CLASSIFY')))
  
  # get genders
  setnames(df2,colnames(df2),paste0(colnames(df2),'1'))
  dchain	<- merge(dchain, df2, by.x='H1',by.y=c('ID1'),all.x=T)
  setnames(df2,colnames(df2),gsub('1','2',colnames(df2)))
  dchain	<- merge(dchain, df2, by.x='H2',by.y=c('ID2'),all.x=T)
  setnames(df2,colnames(df2),gsub('2','',colnames(df2)))
  
  # get couple data
  load(infile.couple)
  couple <- unique(subset(data.table(coupdat), select=c('male.RCCS_studyid', 'female.RCCS_studyid')))
  couple[, COUPLE:=1]
  setnames(couple,c('male.RCCS_studyid','female.RCCS_studyid'),c('PT_ID1','PT_ID2'))
  couple[,PT_ID1:=paste0('RK-',PT_ID1)]
  couple[,PT_ID2:=paste0('RK-',PT_ID2)]
  # couple[,SEX1:='M']
  # couple[,SEX2:='F']
  couple <- merge(couple,subset(df,select=c('PT_ID', 'AID')),by.x='PT_ID1', by.y='PT_ID', all.x=T)
  couple <- merge(couple,subset(df,select=c('PT_ID', 'AID')),by.x='PT_ID2', by.y='PT_ID', all.x=T)
  set(couple,NULL,c('PT_ID1','PT_ID2'),NULL)
  setnames(couple,c('AID.x','AID.y'),c('H1','H2'))
  couple[,H1:=as.character(H1)]
  couple[,H2:=as.character(H2)]
  couple <- couple[!is.na(H1) & !is.na(H2)]
  # tmp <- copy(couple)
  # setnames(tmp,c('H1','H2'),c('H2','H1'))
  # couple <- unique(rbind(couple,tmp))
  
  # order ID
  setnames(dchain,c('H1','H2','SEX1','SEX2'),c('ID1','ID2','ID1_SEX','ID2_SEX'))
  tmp		<- subset(dchain, ID1_SEX=='F' & ID2_SEX=='M')
  setnames(tmp, colnames(tmp), gsub('xx','ID2',gsub('ID2','ID1',gsub('ID1','xx',gsub('xx','12',gsub('12','21',gsub('21','xx',colnames(tmp))))))))
  set(tmp, NULL, 'PHYLOSCANNER_CLASSIFY', tmp[, gsub('xx','12',gsub('12','21',gsub('21','xx',PHYLOSCANNER_CLASSIFY)))])
  dchain	<- rbind(subset(dchain, !(ID1_SEX=='F' & ID2_SEX=='M')), tmp)
  dchain[, PAIR_SEX:= paste0(ID1_SEX,ID2_SEX)]
  unique(dchain$PHYLOSCANNER_CLASSIFY)
  
  # directions
  rtp		<- subset(dchain, !grepl('not resolved|unlinked',PHYLOSCANNER_CLASSIFY))
  # tmp <- rtp[,length(unique( PHYLOSCANNER_CLASSIFY)),by=c('ID1','ID2')]
  # tmp[V1>1]
  rtp[, table(PAIR_SEX)]
  # FF  FU  MF  MM  MU  UF  UM 
  # 95   7 448 113  11   8   5 

  # FF  FU  MF  MM  MU  UF  UM 
  # 146   7 593 137  11  11   7 
  # 
  # check couple
  rtp <- merge(rtp, couple, by.x=c('ID1','ID2'),by.y=c('H1','H2'),all.x=T)
  rtp[!is.na(COUPLE)]
  dchain <- merge(dchain, couple, by.x=c('ID1','ID2'),by.y=c('H1','H2'),all.x=T)  
  cat(nrow(dchain[!is.na(COUPLE) & !grepl('unlinked',PHYLOSCANNER_CLASSIFY)]), ' linked couple \n',
  nrow(dchain[!is.na(COUPLE) & !grepl('not resolved|unlinked',PHYLOSCANNER_CLASSIFY)]), ' linked couple with directions \n',
  nrow(dchain[!is.na(COUPLE)]),' couple \n')
  
  # MIN READ 30 - 
  # 215  linked couple 
  # 168  linked couple with directions 
  # 227  couple 
  # MIN READ NULL - 
  # 291  linked couple 
  # 226  linked couple with directions 
  # 314  couple 

  
  # make dobs
  # dobs <- subset(rtp, select=c('ID1','ID2','LINK_12','LINK_21','ID1_SEX','ID2_SEX','COUPLE'))
  # tmp <- dobs[LINK_12==1]
  # tmp <- subset(tmp,select=c('ID1','ID2','ID1_SEX','ID2_SEX','COUPLE'))
  # setnames(tmp,c('ID1','ID2','ID1_SEX','ID2_SEX'), c('TR_ID','REC_ID','TR_SEX','REC_SEX'))
  # dobs <- dobs[LINK_12==0]
  # dobs <- subset(dobs,select=c('ID1','ID2','ID1_SEX','ID2_SEX','COUPLE'))  
  # setnames(dobs,c('ID1','ID2','ID1_SEX','ID2_SEX'), c('REC_ID','TR_ID','REC_SEX','TR_SEX'))
  # dobs <- rbind(dobs, tmp)
  
  
  # check RCCS
  dobs <- copy(rtp)
  df3 <- subset(df, select=c('PT_ID','AID'))
  df3[,RCCS:=grepl('RK-',PT_ID)]
  df3 <-subset(df3,select=c('AID','RCCS'))
  dobs <- merge(dobs, df3, by.x='ID1',by.y='AID', all.x=T) 
  dobs <- merge(dobs, df3, by.x='ID2',by.y='AID', all.x=T) 
  setnames(dobs,c('RCCS.x','RCCS.y'),c('TR_RCCS','REC_RCCS'))
  
  # get heterosexual 
  tmp <- dobs[ID1_SEX=='M' & ID2_SEX=='F']
  cat(nrow(tmp), ' heterosexual pairs with directions \n ',
  nrow(tmp[TR_RCCS==T & REC_RCCS==T]), ' RCCS')
  # nrow(tmp[(TR_RCCS==T & REC_RCCS==F)|TR_RCCS==F & REC_RCCS==T ])
  
  # 448  heterosexual pairs with directions 
  # 403  RCCS
  # 593  heterosexual pairs with directions 
  # 539  RCCS
}

couple.analysis <- function(){
  library(data.table)
  # load run
  infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
  infile.run='/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_phscinput_runs.rds'
  pty.runs <- data.table(readRDS(infile.run))
  
  # map id and sequences
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  infile.ind.rccs <- file.path(data.dir,'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.mrc <- file.path(data.dir,'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
  id.dt <- data.table(read.csv(infile.ind.rccs))
  id.dt <- subset(id.dt,select = c("pt_id","pangea_id"))
  tmp <- data.table(read.csv(infile.ind.mrc))
  tmp <- subset(tmp,select = c("pt_id","pangea_id"))
  id.dt <- rbind(id.dt,tmp)
  id.dt <- unique(id.dt)
  
  # load couple
  load(infile.couple)
  couple <- data.table(unique(subset(coupdat,select=c('male.RCCS_studyid', 'female.RCCS_studyid'))))
  colnames(couple) <- c('ID1','ID2')
  couple[,COUPLE:=1]
  couple[,ID1:=paste0('RK-',ID1)]
  couple[,ID2:=paste0('RK-',ID2)]
  
  
  # whether couples in run
  cat(nrow(couple), ' couples in total \n',
      nrow(couple[ID1%in% pty.runs$UNIT_ID & ID2%in% pty.runs$UNIT_ID]), ' in runs \n',
      nrow(couple[ID1%in% id.dt$pt_id & ID2%in% id.dt$pt_id]), ' have pangea ids \n')
  tmp <- couple[ID1%in% id.dt$pt_id & ID2%in% id.dt$pt_id & !(ID1%in% pty.runs$UNIT_ID & ID2%in% pty.runs$UNIT_ID)]
  
  # eg
  # WARNING: lose some couples ... check later
  tmp[ID1%in% pty.runs$UNIT_ID & !ID2%in% pty.runs$UNIT_ID]
  merge(tmp[ID1%in% pty.runs$UNIT_ID & !ID2%in% pty.runs$UNIT_ID],id.dt, by.x='ID2',by.y='pt_id',all.x=T)
  tmp[!ID1%in% pty.runs$UNIT_ID & ID2%in% pty.runs$UNIT_ID]
  tmp[!ID1%in% pty.runs$UNIT_ID & !ID2%in% pty.runs$UNIT_ID]
  
  
  # load 19 network
  infile.aid19 <- '~/todi_pairs_171122_cl25_d50_prior23_min30_anonymised_RIDs.csv' 
  infile.pairs19 <- '~/Rakai_phscnetworks_allpairs_190706.rda'
  infile.net19 <- '~/Rakai_phscnetworks_190706.rda'
  aid19 <- data.table(read.csv(infile.aid19))
  aid19[, ID:=paste0('RK-',ID)]
  load(infile.pairs19)
  load(infile.net19)
  dc19 <- copy(dc)
  dw19 <- copy(dw)
  dnet19 <- copy(dnet)
  dchain19 <- copy(dchain)
  
  tmp <- subset(aid19,select=c('ID','AID'))
  dc19 <- merge(dc19, tmp, by.x='H1',by.y='AID',all.x=T)
  dc19 <- merge(dc19, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dc19, c('ID.x','ID.y'),c('ID1','ID2'))
  dw19 <- merge(dw19, tmp, by.x='H1',by.y='AID',all.x=T)
  dw19 <- merge(dw19, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dw19, c('ID.x','ID.y'),c('ID1','ID2'))
  dnet19 <- merge(dnet19, tmp, by.x='H1',by.y='AID',all.x=T)
  dnet19 <- merge(dnet19, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dnet19, c('ID.x','ID.y'),c('ID1','ID2'))
  dchain19 <- merge( dchain19, tmp, by.x='H1',by.y='AID',all.x=T)
  dchain19 <- merge( dchain19, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames( dchain19, c('ID.x','ID.y'),c('ID1','ID2'))
  
  # load network
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/'
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  load( infile.networks )	
  aid <- data.table(read.csv(infile.ind.anonymised))
  
  tmp <- subset(aid,select=c('PT_ID','AID'))
  dc <- merge(dc, tmp, by.x='H1',by.y='AID',all.x=T)
  dc <- merge(dc, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dc, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  dw <- merge(dw, tmp, by.x='H1',by.y='AID',all.x=T)
  dw <- merge(dw, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dw, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  dchain <- merge(dchain, tmp, by.x='H1',by.y='AID',all.x=T)
  dchain <- merge(dchain, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dchain, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  dnet <- merge(dnet, tmp, by.x='H1',by.y='AID',all.x=T)
  dnet <- merge(dnet, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dnet, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  
  # whether couples in dc
  tmp <- unique(subset(dc,select = c('ID1','ID2')))
  tmp[,ID.x:=paste0(ID1,'-',ID2)]
  tmp[,ID.y:=paste0(ID2,'-',ID1)]
  cat(nrow(couple[paste0(ID1,'-',ID2)%in% tmp$ID.x | paste0(ID1,'-',ID2)%in% tmp$ID.y ]), ' couples in networks')
  
  # to dt
  dc19 <- as.data.table(dc19)
  dw19 <- as.data.table(dw19)
  dnet19 <- as.data.table(dnet19)
  dchain19 <- as.data.table(dchain19)
  dchain <- as.data.table(dchain)
  
  # add couple
  tmp <- copy(couple)
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  couple <- data.table(unique(rbind(couple, tmp)))
  
  dc <- merge(dc,couple,by=c('ID1','ID2'),all.x=T)
  dc19 <- merge(dc19,couple,by=c('ID1','ID2'),all.x=T)
  dw <- merge(dw,couple,by=c('ID1','ID2'),all.x=T)
  dw19 <- merge(dw19,couple,by=c('ID1','ID2'),all.x=T)
  dnet <- merge(dnet,couple,by=c('ID1','ID2'),all.x=T)
  dnet19 <- merge(dnet19,couple,by=c('ID1','ID2'),all.x=T)
  dchain <- merge(dchain,couple,by=c('ID1','ID2'),all.x=T)
  dchain19 <- merge(dchain19,couple,by=c('ID1','ID2'),all.x=T)
  
  #
  tmp_dc_coup <- unique(dc[CATEGORISATION=='close.and.adjacent.and.ancestry.cat'&COUPLE==1])
  tmp_dc_coup[,length(N),by=c('ID1','ID2')]
  tmp_dc_coup[ID1=='RK-A035039' & ID2=='RK-D061829']
  tmp <- tmp_dc_coup[,list(N_EFF=max(N_EFF)),by=c('ID1','ID2')]
  tmp_dc_coup <- merge(tmp_dc_coup,tmp, by=c('ID1','ID2','N_EFF'))
  setkey(tmp_dc_coup, ID1, ID2, PTY_RUN, TYPE)
  tmp_dc_coup <- tmp_dc_coup[,head(.SD, 4),by=c('ID1','ID2')]
  tmp_dc_coup <- dcast(tmp_dc_coup,  ID1 + ID2 + COUPLE ~ TYPE, value.var= 'SCORE')
  tmp_dc_coup[,ANALYSIS:=2021]
  
  tmp_dc19_coup <- unique(dc19[CATEGORISATION=='close.and.adjacent.and.ancestry.cat'&COUPLE==1])
  # tmp_dc19_coup[,length(N),by=c('ID1','ID2')]
  # tmp <- tmp_dc19_coup[,list(N_EFF=max(N_EFF)),by=c('ID1','ID2')]
  # tmp_dc19_coup <- merge(tmp_dc19_coup,tmp, by=c('ID1','ID2','N_EFF'))
  tmp_dc19_coup <- dcast(tmp_dc19_coup,  ID1 + ID2 + COUPLE ~ TYPE, value.var= 'SCORE')
  tmp_dc19_coup[,ANALYSIS:=2019]
  
  tmp_dc_coup <- rbind(tmp_dc19_coup, tmp_dc_coup)
  tmp <- tmp_dc_coup[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2','12','21'),c('ID2','ID1','21','12'))
  tmp_dc_coup <- rbind(tmp_dc_coup[(ID1<ID2)],tmp)
  tmp_dc_coup <- melt(tmp_dc_coup,id.vars=c('ID1','ID2', 'COUPLE','ANALYSIS'))
  
  library(dplyr)
  library(ggplot2)
  df <- arrange(tmp_dc_coup, variable, desc(value))
  
  df$ID=paste0(df$ID1,'-',df$ID2)
  df$ID=factor(df$ID, levels = unique(df$ID))
  ggplot(df, aes(ID,value,fill=variable))+
    geom_bar(position="stack", stat="identity")+
    facet_grid(ANALYSIS~ .)+
    labs(x='ID',y='scores',fill='')+ guides(fill = guide_legend(nrow = 1)) + theme(legend.position = "bottom",legend.direction = 'horizontal') +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ 
    scale_x_discrete(limits = rev(levels(df$ID)))
  ggsave(filename = '~/compare_couple.pdf',width = 10, height = 8)
  tmp_dc_coup <- dcast(tmp_dc_coup,ID1+ID2+COUPLE+variable ~ANALYSIS ,value.var='value')
  tmp1 <- tmp_dc_coup[is.na(`2019`) & !is.na(`2021`) ]
  tmp2 <- tmp_dc_coup[!is.na(`2019`) & is.na(`2021`) ]
  tmp3 <- tmp_dc_coup[!is.na(`2019`) & !is.na(`2021`) ]
  nrow(unique(subset(tmp1,select=c('ID1','ID2'))))
  nrow(unique(subset(tmp2,select=c('ID1','ID2'))))
  nrow(unique(subset(tmp3,select=c('ID1','ID2'))))
  # 2021 ONLY
  tmp <- dcast(tmp1,ID1+ID2~variable,value.var='2021')
  tmp[,LINKED:= `12`+`21` +complex.or.no.ancestry>0.6]
  tmp[LINKED==1,DIRECTED:= (`12`+`21`)/(`12`+`21` +complex.or.no.ancestry)>0.6]
  tmp[DIRECTED==1,LINK12:=`12`>`21`]
  nrow(tmp[LINKED==1])/nrow(tmp)
  nrow(tmp[DIRECTED==1])/nrow(tmp[LINKED==1])
  # 2019 ONLY
  tmp <- dcast(tmp2,ID1+ID2~variable,value.var='2019')
  tmp[,LINKED:= `12`+`21` +complex.or.no.ancestry>0.6]
  tmp[LINKED==1,DIRECTED:= (`12`+`21`)/(`12`+`21` +complex.or.no.ancestry)>0.6]
  tmp[DIRECTED==1,LINK12:=`12`>`21`]
  nrow(tmp[LINKED==1])/nrow(tmp)
  nrow(tmp[DIRECTED==1])/nrow(tmp[LINKED==1])
  tmp <- unique(subset(tmp2,select=c('ID1','ID2')))
  tmp[ID1%in% pty.runs$UNIT_ID & ID2%in% pty.runs$UNIT_ID]
  tmp[,ID.x:=paste0(ID1,'-',ID2)]
  tmp[,ID.y:=paste0(ID2,'-',ID1)]
  dw[paste0(ID1,'-',ID2) %in% tmp$ID.x]
  dw[paste0(ID1,'-',ID2) %in% tmp$ID.y]
  load('~/dcdwina_minreadnull.rda')
  dca[,host.1:=gsub('CNTRL-','',host.1)]
  dca[,host.2:=gsub('CNTRL-','',host.2)]
  dcac <- merge(dca, aid, by.x='host.1',by.y='AID',all.x=T)
  dcac <- merge(dcac, aid, by.x='host.2',by.y='AID',all.x=T)
  dcac[paste0(PT_ID.x,'-',PT_ID.y) %in% tmp$ID.x & categorisation == 'close.and.adjacent.and.ancestry.cat']
  dcac[paste0(PT_ID.x,'-',PT_ID.y) %in% tmp$ID.y & categorisation == 'close.and.adjacent.and.ancestry.cat']
  
  # 2019 and 2021
  tmp21 <- dcast(tmp3,ID1+ID2~variable,value.var='2021')
  tmp21[,LINKED_SCORE:=`12`+`21` +complex.or.no.ancestry]
  tmp21[,LINKED:= LINKED_SCORE>0.6]
  tmp21[LINKED==1,DIRECTED_SCORE:=(`12`+`21`)/LINKED_SCORE]
  tmp21[LINKED==1,DIRECTED:= DIRECTED_SCORE>0.6]
  tmp21[DIRECTED==1,LINK12_SCORE:=`12`/(`12`+`21`)]
  tmp21[DIRECTED==1,LINK12:=`12`>`21`]
  nrow(tmp21[LINKED==1])/nrow(tmp21)
  nrow(tmp21[DIRECTED==1])/nrow(tmp21[LINKED==1])
  
  tmp19 <- dcast(tmp3,ID1+ID2~variable,value.var='2019')
  tmp19[,LINKED_SCORE:=`12`+`21` +complex.or.no.ancestry]
  tmp19[,LINKED:= LINKED_SCORE>0.6]
  tmp19[LINKED==1,DIRECTED_SCORE:=(`12`+`21`)/LINKED_SCORE]
  tmp19[LINKED==1,DIRECTED:= DIRECTED_SCORE>0.6]
  tmp19[DIRECTED==1,LINK12_SCORE:=`12`/(`12`+`21`)]
  tmp19[DIRECTED==1,LINK12:=`12`>`21`]
  nrow(tmp19[LINKED==1])/nrow(tmp19)
  nrow(tmp19[DIRECTED==1])/nrow(tmp19[LINKED==1])
  
  tmp <- merge(tmp19,tmp21,by=c('ID1','ID2'))
  tmp[,LINKED_SAME:=LINKED.x==LINKED.y]
  ggplot(tmp,aes(fill=factor(LINKED_SAME,c(TRUE,FALSE),c('yes','no'))))+
    geom_boxplot(aes(x='2019',y=LINKED_SCORE.x)) + 
    geom_boxplot(aes(x='2021',y=LINKED_SCORE.y)) +
    labs(x='analysis',y='linkage scores',fill='same linkage classification')+ 
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",legend.direction = 'horizontal',legend.box = 'horizontal') +
    theme_bw()
  ggsave(filename = '~/compare_linkage_score_couple.pdf',width=6,height=4)
  
  
  tmp[LINKED.x==T & LINKED.y==T,DIRECTED_SAME:=DIRECTED.x==DIRECTED.y]
  tmp[DIRECTED.x==T & DIRECTED.y==T,LINK12_SAME:=LINK12.x==LINK12.y]
  ggplot(tmp[LINKED.x==T & LINKED.y==T],aes(fill=factor(DIRECTED_SAME,c(TRUE,FALSE),c('yes','no'))))+
    geom_boxplot(aes(x='2019',y=DIRECTED_SCORE.x)) + 
    geom_boxplot(aes(x='2021',y=DIRECTED_SCORE.y)) +
    labs(x='analysis',y='direction scores',fill='same directed classification')+ 
    guides(fill = guide_legend(nrow = 1)) + 
    theme(legend.position = "bottom",legend.direction = 'horizontal',legend.box = 'horizontal') +
    theme_bw()
  ggsave(filename = '~/compare_direction_score_couple.pdf',width=6,height=4)
  
  ggplot(tmp[DIRECTED.x==T & DIRECTED.y==T],aes(fill=factor(LINK12_SAME,c(TRUE,FALSE),c('yes','no'))))+
    geom_boxplot(aes(x='2019',y=LINK12_SCORE.x)) + 
    geom_boxplot(aes(x='2021',y=LINK12_SCORE.y)) +
    labs(x='analysis',y='1->2 scores',fill='same direction classification')+ 
    guides(fill = guide_legend(nrow = 1)) + 
    theme(legend.position = "bottom",legend.direction = 'horizontal',legend.box = 'horizontal') +
    theme_bw()
  ggsave(filename = '~/compare_link12_score_couple.pdf',width=6,height=4)
  
  
  # tmp[LINKED.x==LINKED.y,`12.x`+`21.x` +complex.or.no.ancestry.x]
  # tmp[LINKED.x==LINKED.y,`12.y`+`21.y` +complex.or.no.ancestry.y]
  # tmp[LINKED.x!=LINKED.y,`12.x`+`21.x` +complex.or.no.ancestry.x]
  # tmp[LINKED.x!=LINKED.y,`12.y`+`21.y` +complex.or.no.ancestry.y]
  
  # couple plots
  library(ggplot2)
  load('~/dcdwina_minreadnull.rda')
  dwina[,host.1:=gsub('CNTRL-','',host.1)]
  dwina[,host.2:=gsub('CNTRL-','',host.2)]
  dwinac <- merge(dwina , aid, by.x='host.1',by.y='AID',all.x=T)
  dwinac <- merge(dwinac , aid, by.x='host.2',by.y='AID',all.x=T)
  dwinac <- merge(dwinac,couple,by.x=c('PT_ID.x','PT_ID.y'),by.y=c('ID1','ID2'),all.x=T)
  tmp <- unique(subset(dwinac[COUPLE==1],select=c('host.1','host.2','PTY_RUN','tree.id', 'patristic.distance')))
  
  ggplot(tmp, aes(patristic.distance))+
    # geom_density()+theme_bw()+labs(x='patristic distance', y='count')+ coord_cartesian(xlim = c(0,0.05))
    geom_histogram(binwidth = 0.002)+theme_bw()+labs(x='patristic distance', y='count')
  # coord_cartesian(ylim = c(0,1e3))
  ggsave(filename = '~/pdistance_couple.pdf',width = 6, height = 4)
  
  ggplot(tmp, aes(patristic.distance))+
    # geom_density()+theme_bw()+labs(x='patristic distance', y='count')+ coord_cartesian(xlim = c(0,0.05))
    geom_histogram(binwidth = 0.002)+theme_bw()+labs(x='patristic distance', y='count') +coord_cartesian(ylim = c(0,1e3),xlim=c(0,0.1)) 
  # coord_cartesian(ylim = c(0,1e3))
  ggsave(filename = '~/pdistance_couple_zoomin.pdf',width = 6, height = 4)
  
  # pd <- tmp[patristic.distance>0.05 & patristic.distance<0.15,]$patristic.distance
  # library(MASS)
  # fit <- fitdistr(pd, "normal")
  # qnorm(0.05,fit$estimate[1],fit$estimate[2])
  
  # make couple plots 
  dw_couple <- dwinac[COUPLE==1]
  tmp <- dw_couple[,list(M=median(patristic.distance),
                         CL=quantile(patristic.distance,probs = 0.025),
                         CU=quantile(patristic.distance,probs = 0.975)),
                   by=c('host.2','host.1')]
  setkey(tmp,M)
  tmp[,ID:=seq_len(nrow(tmp))]
  ggplot(tmp,aes(ID,M))+
    geom_point()+
    theme_bw()+
    geom_errorbar(aes(ymin=CL,ymax=CU))+
    labs(x='couples',y='patristic distance')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  # +
  #   coord_cartesian(ylim=c(0,0.2))
  ggsave('~/patdist_vs_couple.pdf',width = 10, height = 4)
  
  tmpid <- subset(tmp,select = c('host.1','host.2','ID'))
    
  dcac <- merge(dcac,couple,by.x=c('PT_ID.x','PT_ID.y'),by.y=c('ID1','ID2'),all.x=T)
  dc_couple <- dcac[COUPLE==1 & categorisation == 'close.and.adjacent.and.ancestry.cat']
  tmp <- dc_couple[,list(n.eff=max(n.eff)),by=c('host.1','host.2')]
  tmp <- merge(tmp,dc_couple, by=c('host.1','host.2','n.eff'))
  setkey(tmp, host.1, host.2, PTY_RUN, type)
  tmp <- tmp[,head(.SD, 1),by=c('host.1','host.2','type')]
  tmp <- merge(tmp, tmpid, by=c('host.1','host.2'))
  # tmp2 = tmp[,length(n.eff),by=c('host.1','host.2')]
  # tmp2 = tmp[,as.integer(sum(score)),by=c('host.1','host.2')]
  ggplot(tmp,aes(ID,k.eff,fill=type))+
    theme_bw()+
    geom_bar(position="stack", stat="identity")+
    labs(x='couples',y='number of windows \n supporting subgraph topology')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  ggsave('~/topo_num_vs_couple.pdf',width = 10, height = 4)
  
  ggplot(tmp,aes(ID,score,fill=type))+
    theme_bw()+
    geom_bar(position="stack", stat="identity")+
    labs(x='couples',y='proportion of windows \n supporting subgraph topology')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  ggsave('~/topo_prop_vs_couple.pdf',width = 10, height = 4)
  
  tmp <- dc_couple[,list(n.eff=max(n.eff)),by=c('host.1','host.2')]
  tmp <- merge(tmp,dc_couple, by=c('host.1','host.2','n.eff'))
  setkey(tmp, host.1, host.2, PTY_RUN, type)
  tmp <- tmp[,head(.SD, 1),by=c('host.1','host.2','type')]
  tmp <- merge(tmp, tmpid, by=c('host.1','host.2'))
  tmp <- dcast(tmp, host.1+host.2 + ID~type, value.var = 'score')
  tmp[,score_L:=`12`+ `21`+ complex.or.no.ancestry]
  tmp[,score_D:=(`12`+ `21`)/score_L]
  ggplot(tmp,aes(ID,score_L))+
    geom_point()+
    theme_bw()+
    scale_y_continuous(labels = scales::percent,limits = c(0,1))+
    # geom_errorbar(aes(ymin=CL,ymax=CU))+
    labs(x='couples',y='linkage score')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  ggsave('~/lscore_vs_couple.pdf',width = 10, height = 4)
  
  ggplot(tmp,aes(ID,score_D))+
    geom_point()+
    theme_bw()+
    # geom_errorbar(aes(ymin=CL,ymax=CU))+
    labs(x='couples',y='linkage score')+
    scale_y_continuous(labels = scales::percent,limits = c(0,1))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  ggsave('~/dscore_vs_couple.pdf',width = 10, height = 4)
  
}
