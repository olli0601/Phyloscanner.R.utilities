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
  
  # file names
  set.seed(42)
  # downsample <- 50
  downsample <- NULL
  # job_tag <- 'testnoalign'
  # job_tag <- 'testrealign'
  job_tag <- 'refset1'
  job_tag <- 'refset2'
  job_tag <- 'refset3'
  job_tag <- 'refset4'
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
  cat(system(paste0('cp -r ',gsub(paste0(downsample,'$|_[a-z]+$'),'',in.dir),'/2019_New_ConsensusGenomesOneEach_GeneCut.fasta ', 
                    gsub(paste0(downsample,'$|_[a-z]+$'),'',in.dir),'/UgandaKenyaTanzaniaGenomes_GeneCut_TreeOrder.FASTA ', in.dir)))
  
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
      if(grepl('noalign$',job_tag)){
        consensus_seq_all <- c(consensus_seq_oneeach,consensus_seq)
        names(consensus_seq_all) <- paste0('REF_',names(consensus_seq_all))
        write.fasta(sequences=consensus_seq_all,
                    names=names(consensus_seq_all),
                    nbchar = 60,
                    file.out=infile.consensus.all) 
      }else if (grepl('realign$',job_tag)){
        consensus_seq_all <- c(consensus_seq_oneeach,consensus_seq)
        names(consensus_seq_all) <- paste0('REF_',names(consensus_seq_all))
        unique(unlist(lapply(consensus_seq_all, unique)))
        consensus_seq_all  <- lapply(consensus_seq_all, function(x){x[x!='-']})
        consensus_seq_all  <- DNAStringSet(do.call(c,lapply(consensus_seq_all,function(x)paste(x,collapse = ''))))
        consensus_seq_all  <- msa(consensus_seq_all,method="Muscle",type='dna')
        tmp <- msaConvert(consensus_seq_all,type='bio3d::fasta')
        consensus_seq_all <- tmp$ali
        consensus_seq_all <- lapply(seq_len(nrow(consensus_seq_all)), function(i) consensus_seq_all[i,])
        names(consensus_seq_all) <- tmp$id
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
  
  tmp <- do.call(rbind,consensus_seq_all)
  tmp <- apply(tmp, 2, function(x)sum(x=='-'))
  plot(1:length(tmp), tmp)
  
  library(ggmsa)
  g <- ggmsa(infile.consensus.all, char_width = 0.5, seq_name = T)
  ggsave('~/msa_cons.pdf', g, width = 1000, height = 60, limitsize = F)
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
                     default.coord=TRUE
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
  job_tag <- 'testrealign'
  
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
  # infiles	<- subset(infiles, PTY_RUN>=399)
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
  # args <- commandArgs(trailingOnly = TRUE)
  # print(args)
  # downsample <- as.numeric(args[1])
  # downsample <- 50
  downsample <- NULL
  job_tag <- 'testrealign'
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
  ds <- c(50,100,200,300,400,500,'')
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
  ans[nbg=='',nbg:='all']

  g <- ggplot(ans,aes(start,tmrca,color=factor(nbg)))+
    # geom_point()+
    # geom_line() +
    geom_smooth()+
    labs(x='\n genome position',y='time to most recent common ancester \n ',color='reference sequences') +
    theme_bw()+
    theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')+
    scale_x_continuous(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05))+guides(color=guide_legend(nrow=1))
  ggsave(file.path(p,'plots','tmrca_ds.pdf'),g,width = 8,height = 6)
  
  tmp <- ans[nbg=='all',]
  ans <- ans[nbg!='all',]
  ans <- merge(ans, tmp, by='start')
  g <- ggplot(ans,aes(start,tmrca.y-tmrca.x,color=factor(nbg.x)))+
    # geom_point()+
    # geom_line() +
    geom_smooth()+
    labs(x='\n genome position',y='difference in time to most recent common ancester \n ',color='reference sequences') +
    theme_bw()+
    theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')+
    scale_x_continuous(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05))+guides(color=guide_legend(nrow=1))
  ggsave(file.path(p,'plots','tmrca_diff_ds.pdf'),g,width = 8,height = 6)
  
  
  ans <- data.table()
  ds <- c('_testrealign','')
  for (i in 1:length(ds)) {
    p <- paste0('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_output',ds[i],'/ptyr1_trees')
    files <-list.files(p)
    files <- grep('treefile$',files,value = T)
    files <- gsub(paste0(gsub('trees$','InWindow',basename(p)),'_'),'',files)
    files <- as.numeric(gsub("([0-9]+).*$", "\\1", files))
    if (ds[i]==''){
      load(file.path(p,'plots','s4_treeinfo.rda'))
    }else{
      load(file.path(p,'plots','treeinfo.rda'))}
    tmp <- data.table(start=files, 
                      subr=age_root[,1],
                      tmrca=age_root[,3],
                      nnode=age_root[,5], 
                      nbg = ds[i])
    ans <- rbind(ans, tmp)
  }
  ans[,nbg:=factor(nbg,ds,c('all sequences','two sets'))]
  g <- ggplot(ans,aes(start,tmrca,color=factor(nbg)))+
    geom_point()+
    geom_line() +
    labs(x='\n genome position',y='time to most recent common ancester \n ',color='consensus alignment on') +
    theme_bw()+
    theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')+
    scale_x_continuous(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05))+guides(color=guide_legend(nrow=1))
  ggsave(file.path(p,'plots','tmrca_cons_align.pdf'),g,width = 8,height = 6)
}


rkuvri.make.alignments<- function()
{
  #
  #	set up working environment
  require(Phyloscanner.R.utilities)
  library(data.table)
  library(seqinr)
  library(DECIPHER)
  
  # file names
  set.seed(42)
  prog.pty <- '~/phyloscanner/phyloscanner_make_trees.py'
  HOME <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'	
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  potential.networks.analysis.dir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100"
  in.dir <- file.path(HOME,'210211_phsc_input')		
  work.dir <- file.path(HOME,"210211_phsc_work")
  out.dir <- file.path(HOME,"210211_phsc_output")	
  package.dir <- file.path(.libPaths(),'Phyloscanner.R.utilities')
  dir.create(in.dir)
  dir.create(work.dir)
  dir.create(out.dir)
  use.both.consensus <- TRUE # false options is out of date
  downsample <- 50
  infile.runs <- paste0(HOME,'210120_RCCSUVRI_phscinput_runs.rds')
  infile.consensus.oneeach <- file.path(in.dir,'2019_New_ConsensusGenomesOneEach_GeneCut.fasta')
  infile.hxb2.package <- file.path(package.dir,'HIV1_compendium_AD_B_CPX_v2.fasta')
  if(use.both.consensus){
    # infile.consensus.tree <- file.path(in.dir,'UgandaKenyaTanzaniaGenomes_GeneCut_Tree.tre')
    infile.consensus <- file.path(in.dir,'UgandaKenyaTanzaniaGenomes_GeneCut_TreeOrder.FASTA') 
  }
  
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
  
  # clean taxa names 
  consensus_seq_oneeach<- read.fasta(infile.consensus.oneeach)
  dn <- data.table(RAW_NAME=names(consensus_seq_oneeach)) 
  dn[,NUM:=as.numeric(gsub(".*\\((.*)\\).*", "\\1", RAW_NAME))]
  dn[,NAME:=gsub('\\(.*?\\)', '', RAW_NAME)]
  names(consensus_seq_oneeach) <- dn$NAME
  
  # remove sequences
  if(use.both.consensus){
    consensus_seq<- read.fasta(infile.consensus)
    # compare two sets of consensus, and combine
    tmp = sapply(consensus_seq_oneeach,function(x){any(sapply(consensus_seq, function(y){identical(as.character(x)[as.character(x)!='-'],as.character(y)[as.character(y)!='-'])}))})
    consensus_seq_oneeach <- consensus_seq_oneeach[setdiff(names(consensus_seq_oneeach),c(names(tmp[tmp==T]),'CON_OF_CONS'))]
  }else{
    consensus_seq_oneeach <- consensus_seq_oneeach[setdiff(names(consensus_seq_oneeach),'CON_OF_CONS')]
  }
  opn <- grep('^CON_N$|^CON_O$|^CON_P$|CONSENSUS_N$|CONSENSUS_O$|CONSENSUS_P$',names(consensus_seq_oneeach),value=T)
  consensus_seq_oneeach <- consensus_seq_oneeach[setdiff(names(consensus_seq_oneeach),opn)]

  # filter sequences
  if(use.both.consensus){
     sample_date <- as.numeric(sapply(strsplit(names( consensus_seq),split='\\.'),function(x)x[3]))
     id <- which(sample_date<=2010|sample_date>=1990)
     consensus_seq<- consensus_seq[names(consensus_seq)[id]]
  }
  
  if(0){
    infile.consensus.all <- file.path(in.dir,'ConsensusGenomesNoAlign.fasta')
    consensus_seq_all <- c(consensus_seq_oneeach,consensus_seq)
    names(consensus_seq_all) <- paste0('REF_',names(consensus_seq_all))
    write.fasta(sequences=consensus_seq_all,
                names=names(consensus_seq_all),
                nbchar = 60,
                file.out=infile.consensus.all)
  }
  
  if(0){
    # not tested
    infile.consensus.all <- file.path(in.dir,'ConsensusGenomesReAlign.fasta')
    consensus_seq_all <- c(consensus_seq_oneeach,consensus_seq)
    names(consensus_seq_all) <- paste0('REF_',names(consensus_seq_all))
    write.fasta(sequences=consensus_seq_all,
                names=names(consensus_seq_all),
                nbchar = 60,
                file.out=infile.consensus.all)
  }
  
  # save to file
  if(use.both.consensus){
    if(!exists('downsample')){
      infile.consensus.all <- file.path(in.dir,'ConsensusGenomes.fasta')
      if(!file.exists(infile.consensus.all)){
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
      }else{
        consensus_seq_all <- read.fasta(infile.consensus.all)
      }
    }else{
      infile.consensus.all <- file.path(in.dir,paste0('ConsensusGenomes_ds',downsample,'.fasta'))
      if(!file.exists(infile.consensus.all)){
        for (ds in c(50,100,200,300,400,500)) {
          consensus_seq_ds <- consensus_seq_all[paste0('REF_',c(names(consensus_seq_oneeach),sample(names(consensus_seq),ds)))]
          infile.consensus.ds <- file.path(in.dir,paste0('ConsensusGenomes_ds',ds,'.fasta'))
          write.fasta(sequences=consensus_seq_ds,
                      names=names(consensus_seq_ds),
                      nbchar = 60,
                      file.out=infile.consensus.ds)
        }
      }
    }
    
  }else{
    infile.consensus.all <- file.path(in.dir,'ConsensusGenomes2.fasta')
    if(!file.exists(infile.consensus.all)){
      consensus_seq_all <- copy(consensus_seq_oneeach)
      names(consensus_seq_all) <- paste0('REF_',names(consensus_seq_all))
      write.fasta(sequences=consensus_seq_all,
                  names=names(consensus_seq_all),
                  nbchar = 60,
                  file.out=infile.consensus.all)
    }
    consensus_seq_all <- read.fasta(infile.consensus.all) 
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
  
  if(0){
    out.dir <- file.path(out.dir,'test_consensus_noalign')
    dir.create(out.dir)
  }
  if(exists('downsample')){
    out.dir <- file.path(out.dir,paste0('ds',downsample))
    dir.create(out.dir)
  }
  # out.dir <- gsub('ds','ds',out.dir)
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
                        default.coord=TRUE
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
    # if(i==1|i==pty.c[, max(JOB_ID)]){cat(system(cmd, intern= TRUE))}
  }

  library(ape)
  # p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_output/ptyr1_trees'
  p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_output/ptyr400_trees'
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


rkuvri.make.trees<- function() 
{
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
  downsample <- 50
  iqtree.pr <- 'iqtree'
  iqtree.args			<- ifelse(hpc.nproc==1, '-m GTR+F+R6 -ntmax 1 -seed 42 -o REF_CON_H', 
                         paste0('-m GTR+F+R6 -ntmax ',hpc.nproc,' -seed 42 -o REF_CON_H'))
  in.dir				<- file.path(HOME,'210211_phsc_output')		
  out.dir				<- in.dir
  work.dir			<- file.path(HOME,"210211_phsc_work")
  if(exists('downsample')){
    in.dir <- file.path(in.dir,paste0('ds',downsample))
    out.dir	<- in.dir
  }
  
  infiles	<- data.table(FI=list.files(in.dir, pattern='fasta$', full.names=TRUE, recursive=TRUE))
  infiles[, FO:= gsub('.fasta$','',FI)]
  infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FI)))]
  infiles[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FI)))]		
  setkey(infiles, PTY_RUN, W_FROM)	
  
  infiles	<- subset(infiles, PTY_RUN<=14)
  # infiles	<- subset(infiles, PTY_RUN>=399)
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
  library(treedater)

  # files
  set.seed(42)
  args <- commandArgs(trailingOnly = TRUE)
  print(args)
  method <- as.numeric(args[1])
  method <- 4
  p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_output/ptyr1_trees'
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
    
    # filter
    tmp <- tmpdf[grep('REF_',TAXA),]
    tmpdf <- tmpdf[!grep('REF_',TAXA),]
    if(method==1){
      tmp_dsdt_con <- dsdt_con[LOCAL==F & INCLU_AD==T,]
    }else if(method==2){
      tmp_dsdt_con <- dsdt_con[LOCAL==F,]
    }else if(method==3){
      tmp_dsdt_con <- dsdt_con[LOCAL==F|(LOCAL==T & INCLU_AD==T),]
    }else if(method==4){
      tmp_dsdt_con <- copy(dsdt_con)
    }else if(method==5){
      tmp_dsdt_con <- rbind(dsdt_con[LOCAL==F,],dsdt_con[LOCAL==T,][1:50,])
    }else if(method==6){
      tmp_dsdt_con <- rbind(dsdt_con[LOCAL==F,],dsdt_con[LOCAL==T,][1:100,])
    }else if(method==7){
      tmp_dsdt_con <- rbind(dsdt_con[LOCAL==F,],dsdt_con[LOCAL==T,][1:200,])
    }else if(method==8){
      tmp_dsdt_con <- rbind(dsdt_con[LOCAL==F,],dsdt_con[LOCAL==T,][1:300,])
    }else{
      stop('method should take values from 1 to 8')
    }
    tmp <- merge(tmp, subset(tmp_dsdt_con,select=c('RENAME_ID','SAMPLE_DATE')),  by='RENAME_ID')
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


    # library(phytools)
    # p1 <- drop.tip(tmptr2,'REF_CON_H')
    # p2 <- drop.tip(p1,'REF_CON_H')
    # plotTree(p1,fsize=0.6,lwd=1)
    # plotTree(p2,fsize=0.6,lwd=1)
    # root    
    tmptr2 <- root(tmptr2, 'REF_CON_H')
    tmptr2 <- drop.tip(tmptr2,'REF_CON_H')
    # a <- rtt(tmptr2,sts[names(sts) %in% tmptr2$tip.label])
    # # b <- estimate.dates(tmptr2,sts[names(sts) %in% tmptr2$tip.label])
    # b <- estimate.dates(a, sts[names(sts) %in% tmptr2$tip.label])
    # rootToTipRegressionPlot(a)
    # 
    # tmptr2_t <- tmptr2$tip.label
    # tmptr2_t_ref <-grep('REF_',tmptr2_t,value = T)
    # tmptr2_t_nref <- length(tmptr2_t_ref)
    # tmptr2_t_naid <- length(tmptr2_t) -  tmptr2_t_nref
    # unique(sts[names(sts) %in% tmptr2_t_ref])
    # unique(sts[! names(sts) %in% tmptr2_t_ref])
    # hist(sts[names(sts) %in% tmptr2_t_ref])
      
    # tmptrd <- dater(unroot(tmptr2), sts=sts,s = tmplen, clock='strict',numStartConditions = 0,estimateSampleTimes=sts.df,omega0 = 0.001)
    tmptrd <- dater(tmptr2, sts=sts,s = tmplen, clock='strict',numStartConditions = 0,estimateSampleTimes=sts.df,omega0 = 0.001)
    # tmptrd1 <- dater(tmptr2, sts=sts,s = tmplen, clock='strict',numStartConditions = 0,estimateSampleTimes=sts.df,omega0 = 0.001,temporalConstraints = F)
    # tmptrd2 <- dater(tmptr2, sts=sts,s = tmplen, clock='strict',numStartConditions = 0,estimateSampleTimes=sts.df,omega0 = 0.001,clsSolver = "limSolve")
    # tmptrd3 <- dater(tmptr2, sts=sts,s = tmplen, clock='strict',numStartConditions = 0,estimateSampleTimes=sts.df,omega0 = 0.001,clsSolver = "mgcv")
    # 
    # tmptrd <- dater(tmptr2, sts=sts,s = tmplen, clock='uncorrelated',numStartConditions = 0,estimateSampleTimes=sts.df,omega0 = 0.001)
    # 
    cat('save root to tip regression plot \n')
    pdf(file.path(p,'plots',paste0('s',method,'_',gsub('.treefile','_rttreg.pdf',files[i]))),width = 6,height = 4)
    fit <- rootToTipRegressionPlot(tmptrd)
    
    tmp <- as.data.table(cbind(fit$model,tmptrd$tip.label) )
    # plot(as.numeric(fit$model[,2]),as.numeric(fit$model[,1]),xlim = c(1980,2020), ylim = c(0,1))
    # lm(as.numeric(fit$model[,1])~as.numeric(fit$model[,2]))
    # sort(unique(gsub('^([A-Z]+[0-9]+)-.*','\\1',grep('AID',tmptrd$tip.label,value = T))))
    colnames(tmp) <- c('dist','time','name')
    # tmp[grep('AID4910',name)]
    save(fit,tmp,tmptrd,tmptr2,file = '~/t1w1000.rda')
    save(fit,tmp,tmptrd,tmptr2,file = '~/t1w1025.rda')
    library(data.table)
    library(ape)
    load('~/t1w1000.rda')
    fit1 <- fit
    dt1 <- copy(tmp)
    trd1 <- tmptrd
    tr1 <- tmptr2
    load('~/t1w1025.rda')
    fit2 <- fit   
    dt2 <- copy(tmp)
    trd2 <- tmptrd
    tr2 <- tmptr2
    df <- rbind(data.table(value=node.depth.edgelength(tr1), dated='no', window=1000),
                data.table(value=node.depth.edgelength(trd1), dated='yes', window=1000),
                data.table(value=node.depth.edgelength(tr2), dated='no', window=1025),
                data.table(value=node.depth.edgelength(trd2), dated='yes', window=1025))
    ggplot(df[dated=='yes',],aes(value,fill=factor(window)))+
      geom_density() +
      theme_bw() +
      labs(x='depths',fill='start of window')+
      ggtitle('depth after dating')
    
    ggplot(df[dated=='no',],aes(value,fill=factor(window)))+
      geom_density() +
      theme_bw() +
      labs(x='depths',fill='start of window')+
      ggtitle('depth before dating')
    # dt <- merge(dt1, dt2, by='name', all=T)
    dt1[,W:=1000]
    dt2[,W:=1025]    
    dt <- rbind(dt1,dt2)
    # tmp <- dt[!is.na(time.x)& !is.na(time.y)]
    library(ggplot2)
    ggplot(dt,aes(time,dist,color=factor(W)))+
      geom_point()
    -fit1$coefficients[1]/ fit1$coefficients[2]
    -fit2$coefficients[1]/ fit2$coefficients[2]
    # m1 <- lm(dist~time,dt1)
    # m2 <- lm(dist~time,dt2)
    
    dev.off()
    outliers <- outlierTips( tmptrd)
    cat(nrow(outliers[outliers$q < 0.05,]),'outliers \n')
    age_root[i,] <- c(tmptrd$adjusted.mean.rate,
                      tmptrd$mean.rate,
                      tmptrd$timeOfMRCA,
                      tmptrd$timeToMRCA ,
                      tmptrd$Nnode ,
                      tmptrd$coef_of_variation)
    
    # library(castor)
    # rtd <- get_all_distances_to_root(tmptrd)
    str(tmptrd)
    str(outliers)
    
    # plot
     if(method==1){
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
  
  if(method==1){
    save(blen,age_root,file=file.path(p,'plots',paste0('s',method,'_treeinfo.rda')))
  }else{
    save(age_root,file=file.path(p,'plots',paste0('s',method,'_treeinfo.rda')))
  }
  
  library(data.table)
  library(ggplot2)
  p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_output/ptyr1_trees'
  files <-list.files(p)
  files <- grep('treefile$',files,value = T)
  files <- gsub(paste0(gsub('trees$','InWindow',basename(p)),'_'),'',files)
  files <- as.numeric(gsub("([0-9]+).*$", "\\1", files))
  ans <- data.table()
  for (i in 1:8) {
    load(file.path(p,'plots',paste0('s',i,'_treeinfo.rda')))
    tmp <- data.table(start=files, blen = blen,
                      subr=age_root[,1],tmrca=age_root[,3],nnode=age_root[,5], method=i)
    ans <- rbind(ans, tmp)
  }

  g <- ggplot(ans[method>=1 & method<=4],aes(start,tmrca,color=factor(method,1:4, 
                                               c('A/D subtype representatives that are actual sequences with a sample date',
                                                 'subtype representatives that are actual sequences with a sample date',
                                                 'A/D subtype representatives and Tanzania/Uganda/Kenya background sequences with a sample date',
                                                 'subtype representatives and Tanzania/Uganda/Kenya background sequences with a sample date')))) +
    geom_point()+
    geom_line() +
    labs(x='\n genome position',y='time to most recent common ancester \n ',color='reference sequences') +
    theme_bw()+
    theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')+
    scale_x_continuous(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05)) 
  ggsave(file.path(p,'plots','tmrca.pdf'),g,width = 8,height = 6)
  
  g <- ggplot(ans[method==4],aes(start,tmrca,color=factor(method,1:4,  c('A/D subtype representatives that are actual sequences with a sample date',
                                                                         'subtype representatives that are actual sequences with a sample date',
                                                                         'A/D subtype representatives and Tanzania/Uganda/Kenya background sequences with a sample date',
                                                                         'subtype representatives and Tanzania/Uganda/Kenya background sequences with a sample date')))) +
    geom_point()+
    geom_line() +
    labs(x='\n genome position',y='time to most recent common ancester \n ',color='reference sequences') +
    theme_bw()+
    theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')+
    scale_x_continuous(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05)) 
  ggsave(file.path(p,'plots','tmrca4.pdf'),g,width = 8,height = 6)  
  ans[,list(CL=quantile(tmrca,0.025),
            CU=quantile(tmrca,0.975),
            M=quantile(tmrca,0.5)),by='method']
  # method       CL       CU        M
  # 1:      1 1615.208 1983.304 1912.139
  # 2:      2 1835.966 1977.125 1933.955
  # 3:      3 1892.995 1980.468 1943.616
  # 4:      4 1923.706 1976.923 1958.338
  
  ans[method==4 & tmrca < 1900]
  # start     blen        subr    tmrca nnode method
  # 1:  1000 34.00372 0.001195281 1861.041   752      4
  # 2:  5900 43.26954 0.001449136 1855.946   593      4
  ans[method==4 & start >=800 & start <=1200] 
  ans[method==4 & start >=5700 & start <=6100] 
  g <- ggplot(ans[method>=1 & method<=4],aes(start,nnode,fill=factor(method,1:4,  c('A/D subtype representatives that are actual sequences with a sample date',
                                                             'subtype representatives that are actual sequences with a sample date',
                                                             'A/D subtype representatives and Tanzania/Uganda/Kenya background sequences with a sample date',
                                                             'subtype representatives and Tanzania/Uganda/Kenya background sequences with a sample date')))) +
    geom_bar(position="dodge", stat="identity")+
    labs(x='\n genome position',y='number of nodes\n ',fill='reference sequences') +
    theme_bw()+
    theme(legend.position = 'bottom',legend.direction = 'vertical',legend.box = 'vertical')+
    scale_x_continuous(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05))
  ggsave(file.path(p,'plots','nnode.pdf'),g,width = 8,height = 6)
  
  g <- ggplot(ans[method>=4 & method<=8],aes(start,tmrca,color=factor(method,4:8, 
                                                                      c('all',50,100,200,300)))) +
    geom_point()+
    geom_line() +
    labs(x='\n genome position',y='time to most recent common ancester \n ',color='number of Tanzania/Uganda/Kenya background sequences') +
    theme_bw()+
    theme(legend.position = 'bottom',legend.direction = 'horizontal',legend.box = 'vertical')+
    scale_x_continuous(expand = c(0,0.05))+
    scale_y_continuous(expand = c(0,0.05)) 
  ggsave(file.path(p,'plots','tmrca_ds.pdf'),g,width = 8,height = 4)
  
  # p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_output/ptyr400_trees'
  # files <-list.files(p)
  # files <- grep('treefile$',files,value = T)
  # files <- gsub(paste0(gsub('trees$','InWindow',basename(p)),'_'),'',files)
  # files <- as.numeric(gsub("([0-9]+).*$", "\\1", files))
  # tmp2 <- load(file.path(p,'plots','treeinfo.rda'))
  # ans <- rbind(ans,data.table(start=files, blen = tmp2 , run=400))
  # 

  # g <- ggplot(ans,aes(start,blen,color=factor(run))) +
  #   geom_point()+
  #   geom_line() +
  #   labs(x='\n genome position',y='branch length\n ',color='run') +
  #   theme_bw()+
  #   scale_x_continuous(expand = c(0,0.05))+
  #   scale_y_continuous(expand = c(0,0.05))
  # ggsave(file.path(p,'plots','branchlen.pdf'),g,width = 6,height = 4)
  
  library(data.table)
  out.dir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_output/'
  # p <- 'ptyr1_trees'
  p <- 'ptyr400_trees'
  files <-list.files(file.path(out.dir,p),include.dirs = F)
  df <- data.table(start= seq(800,9400,by=25))
  tmp <- data.table(fa = grep('fasta$',files,value=T))
  tmp[, start := as.numeric(gsub(paste0(gsub('trees$','InWindow_',p),'([0-9]+)_.*'),'\\1',fa))]
  df <- merge(df, tmp, by='start',all.x=T)
  tmp <- data.table(tr = grep('treefile$',files,value=T))
  tmp[, start := as.numeric(gsub(paste0(gsub('trees$','InWindow_',p),'([0-9]+)_.*'),'\\1',tr))]
  df <- merge(df, tmp, by='start',all.x=T)  
  # 'ptyr1_trees'
  df[is.na(fa),start]
  # 7950 9250 9275 9300 9325 9350 9375 9400
  
  df[!is.na(fa)&is.na(tr),start]
  # [1]  975 1050 1075 1175 1275 1375 1475 1650 1675 1750 1775 1825 1850 1875 1925
  # [16] 1950 1975 2025 2050 2075 2125 2150 2175 2250 2275 2350 2375 4375 5075 5175
  # [31] 5275 5375 5450 5475 5575 5675 5850 5875 5975 6075 6175 6275 6375 6450 6475
  # [46] 6550 6575 6650 6675 6775 6875 6975 7075 7150 7175 7250 7275 7375 7450 7475
  # [61] 7550 7575 7775 7875 8000 8100 8200 8300 8400 8475 8500 8600 8700 8800 8900
  # [76] 8975 9000 9100 9200
  
  # 'ptyr400_trees'  
  # numeric(0)
  
  library(seqinr)
  library(data.table)
  library(ape)
  library(ggplot2)
  p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_output/ptyr1_trees'
  files <-list.files(p)
  files <- grep('fasta$',files,value = T)
  files_start <- gsub(paste0(gsub('trees$','InWindow',basename(p)),'_'),'',files)
  files_start <- as.numeric(gsub("([0-9]+).*$", "\\1", files_start))
  df <- data.table()
  setwd(p)
  for (i in 1:length(files)) {
    tmp <- checkAlignment(read.dna(files[i],format = 'fasta'),plot=F)
    fa <- read.fasta(files[i])
    # unique(unlist(lapply(fa, unique)))
    fa <- do.call(rbind,fa)
    tmp <- data.table(ngap = apply(fa,2,function(x)sum(x=='-')),
                      ntotal = nrow(fa),
                      pos=files_start[i],
                      file=files[i])
    df <- rbind(df,tmp)
  }
  df[,ID:=seq_along(ntotal),by='pos']
  save(df,file='~/gap_debug.rda')
  load('~/gap_debug.rda')
  df[,pgap:=ngap/ntotal]
  g <- ggplot(df,aes(ID,pgap))+
    geom_point()+
    facet_wrap(.~factor(pos),ncol=10)+labs(x='genome position',y='proportion of gaps')
  ggsave('~/pgap.pdf', g, width = 20,height = 40)
  
  df[,allgap:=pgap>0.8]
  df[,count_allgap:=ave(allgap, cumsum(!allgap), FUN = cumsum), by='pos']
  tmp <- df[,max(count_allgap), by='pos']
  setkey(tmp, V1)
  tmp[V1>=3]
  tmp[V1>=6]
  # df[pgap==1]
  load(file.path(p,'plots',paste0('s1_treeinfo.rda')))
  files <-list.files(p)
  files <- grep('treefile$',files,value = T)
  files <- gsub(paste0(gsub('trees$','InWindow',basename(p)),'_'),'',files)
  files <- as.numeric(gsub("([0-9]+).*$", "\\1", files))
  dblen <- data.table(blen=blen, pos=files)
  tmp2 <- merge(tmp,dblen,by='pos')
  ggplot(tmp2,aes(V1,blen))+
    geom_point()+
    labs(x='\n longest consecutive gap positions', y = 'branch lengths')
  ggsave('~/lgap_blen.pdf', g, width = 6,height = 4)
}





