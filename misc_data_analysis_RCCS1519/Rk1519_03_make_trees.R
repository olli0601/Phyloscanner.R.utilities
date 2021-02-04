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


rkuvri.make.alignments<- function()
{
  #
  #	set up working environment
  require(Phyloscanner.R.utilities)
  library(data.table)
  library(seqinr)
  library(DECIPHER)
  
  # file names
  prog.pty <- '~/phyloscanner/phyloscanner_make_trees.py'
  HOME <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'	
  in.dir <- file.path(HOME,'210122_phsc_input')		
  work.dir <- file.path(HOME,"210122_phsc_work")
  out.dir <- file.path(HOME,"210122_phsc_output")	
  package.dir <- file.path(.libPaths(),'Phyloscanner.R.utilities')
  dir.create(in.dir)
  dir.create(work.dir)
  dir.create(out.dir)
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  potential.networks.analysis.dir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100"
  infile.runs <- paste0(HOME,'210120_RCCSUVRI_phscinput_runs.rds')
  # infile.consensus.tree <- file.path(in.dir,'UgandaKenyaTanzaniaGenomes_GeneCut_Tree.tre')
  infile.consensus <- file.path(in.dir,'UgandaKenyaTanzaniaGenomes_GeneCut_TreeOrder.FASTA')
  infile.consensus.oneeach <- file.path(in.dir,'2019_New_ConsensusGenomesOneEach_GeneCut.fasta')
  # infile.consensus <- infile.consensus.Tanzania
  infile.hxb2.package <- file.path(package.dir,'HIV1_compendium_AD_B_CPX_v2.fasta')

  # load runs
  pty.runs <- readRDS(infile.runs)
  
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
  consensus_seq<- read.fasta(infile.consensus)
  
  # compare two sets of consensus, and combine
  tmp = sapply(consensus_seq_oneeach,function(x){any(sapply(consensus_seq, function(y){identical(as.character(x)[as.character(x)!='-'],as.character(y)[as.character(y)!='-'])}))})
  consensus_seq_oneeach <- consensus_seq_oneeach[setdiff(names(consensus_seq_oneeach),c(names(tmp[tmp==T]),'CON_OF_CONS'))]
  opn <- grep('^CON_N$|^CON_O$|^CON_P$|CONSENSUS_N$|CONSENSUS_O$|CONSENSUS_P$',names(consensus_seq_oneeach),value=T)
  consensus_seq_oneeach <- consensus_seq_oneeach[setdiff(names(consensus_seq_oneeach),opn)]
  # consensus_seq_all <- c(consensus_seq, consensus_seq_oneeach)
  
  # align
  # tmp <- msa(DNAStringSet(do.call(c,lapply(consensus_seq_all,function(x)paste(x,collapse = '')))),method="Muscle",type='dna')
  tmp <- AlignProfiles(DNAStringSet(do.call(c,lapply(consensus_seq,function(x)paste(x,collapse = '')))),
                       DNAStringSet(do.call(c,lapply(consensus_seq_oneeach,function(x)paste(x,collapse = '')))))
  consensus_seq_all <- as.list(tmp)
  consensus_seq_all <- lapply(consensus_seq_all,function(x){as.character(x)})
  # save
  names(consensus_seq_all) <- paste0('REF_',names(consensus_seq_all))
  infile.consensus.all <- file.path(in.dir,'ConsensusGenomes.fasta')
  write.fasta(sequences=consensus_seq_all,
              names=names(consensus_seq_all),
              nbchar = 60,
              file.out=infile.consensus.all)
  
  consensus_seq_all <- read.fasta(infile.consensus.all) 
  
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
  # tmp <- names(consensus_seq)
  # tmp <- gsub('REF_','',tmp)
  # tmp <- lapply(strsplit(tmp, split = '[.]'), as.numeric)
  # tmp <- lapply(tmp, function(x)x[!is.na(x)])
  # tmp2 <- unlist(lapply(tmp, function(x)x[x>1000 & x<2021]))
  # root.seq <- which(sapply(tmp, function(x){min(tmp2) %in% x}))[1]  
  # root.seq <- names(consensus_seq)[root.seq]
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
                        select=NA	
    )											
    pty.c <- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)
    pty.c[, W_FROM:= ptyi[i]]
    # pty.c[, PTY_RUN:= as.integer(sub('.*ptyr([0-9])_.*','\\1',CMD))]
    pty.c
  })
  pty.c	<- do.call('rbind', pty.c)	
  setkey(pty.c,PTY_RUN,W_FROM)
  # tmp		<- data.table(FO=list.files(out.dir, pattern='ptyr.*fasta$', recursive=TRUE, full.names=TRUE))
  # tmp[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FO)))]
  # tmp[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FO)))]
  # pty.c	<- merge(pty.c, tmp, by=c('PTY_RUN','W_FROM'), all.x=1)		
  # pty.c	<- subset(pty.c, is.na(FO))		
  # pty.c[, CASE_ID:= seq_len(nrow(pty.c))]
  #pty.c	<- subset(pty.c, W_FROM==2350)
  #print(pty.c, n=1e3)
  pty.c[, CASE_ID:= rep(1:4900,times=ceiling(nrow(pty.c)/4900))[1:nrow(pty.c)]]
  pty.c[, JOB_ID:= rep(1:ceiling(nrow(pty.c)/4900),each=4900)[1:nrow(pty.c)]]
  
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

}


rkuvri.make.trees<- function() 
{
  require(data.table)
  require(Phyloscanner.R.utilities)
  #
  #	produce trees
  #
  hpc.load			<- "module load intel-suite/2015.1 mpi"
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
  raxml.pr <- 'iqtree'
  raxml.args			<- ifelse(hpc.nproc==1, '-m GTR+F+R6 -ntmax 1 -seed 42 -o REF_CON_H', 
                         paste0('-m GTR+F+R6 -ntmax ',hpc.nproc,' -seed 42 --o REF_CON_H'))
  in.dir				<- file.path(HOME,'210122_phsc_output')		
  out.dir				<- in.dir
  work.dir			<- file.path(HOME,"210122_phsc_work")
  
  infiles	<- data.table(FI=list.files(in.dir, pattern='fasta$', full.names=TRUE, recursive=TRUE))
  infiles[, FO:= gsub('fasta$','phy',FI)]
  infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FI)))]
  infiles[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FI)))]		
  # tmp		<- data.table(FT=list.files(out.dir, pattern='^ptyr.*tree$', full.names=TRUE, recursive=TRUE))
  # tmp[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FT)))]
  # tmp[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FT)))]
  # infiles	<- merge(infiles, tmp, by=c('PTY_RUN','W_FROM'), all.x=1)
  # infiles	<- subset(infiles, is.na(FT))	
  setkey(infiles, PTY_RUN, W_FROM)	
  
  infiles	<- subset(infiles, PTY_RUN<=14)
  # infiles	<- subset(infiles, PTY_RUN>=399)
  print(infiles)	 		
  
  df<- infiles[, list(CMD=cmd.iqtree(FI, outfile=FO, pr=raxml.pr, pr.args=raxml.args)), by=c('PTY_RUN','W_FROM')]
  df[, ID:=ceiling(seq_len(nrow(df))/5)]
  df<- df[, list(CMD=paste(CMD, collapse='\n',sep='')), by='ID']
  
  #df[1, cat(CMD)]
  #stop()
  invisible(df[,{
    cmd	<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
    cmd <- paste(cmd,'module load anaconda3/personal \n source activate phylo' ,sep='\n')
    cmd	<- paste(cmd,CMD,sep='\n')
    cat(cmd)							
    outfile	<- paste("srx1",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
    cmd.hpccaller(work.dir, outfile, cmd)
    #stop()
  }, by=c('ID')])
}






