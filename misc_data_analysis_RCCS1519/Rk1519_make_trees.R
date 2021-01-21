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
  
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  out.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_'
	potential.networks.analysis.dir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100"
	infile.potential.networks <- file.path(potential.networks.analysis.dir, 'clusters.rda')
	infile.clostest.three.per.cluster <- file.path(potential.networks.analysis.dir, 'closest_3_percluster.rds')
	infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
	
	
  # load clusters of close pairs and redefine IDCLU
	load( infile.potential.networks )  
	tmp <- unique(subset(df,select=c('CLU_SIZE','IDCLU')))
	setkey(tmp, CLU_SIZE)
	tmp[, IDCLU2:=seq_len(nrow(tmp))]  
	tmp[, CLU_SIZE:=NULL]
	rtc	<- subset( merge(df, tmp, by='IDCLU'),select=c('ID','IDCLU2','CLU_SIZE'))
	setnames(rtc,'IDCLU2','IDCLU')
	setkey(rtc, IDCLU)
	cat(max(rtc$IDCLU),' clusters of ',length(unique(rtc$ID)), ' individuals, the ten largest cluster sizes are ', 
			unique(subset(rtc,select=c('CLU_SIZE','IDCLU')))$CLU_SIZE[(max(rtc$IDCLU)-10+1):max(rtc$IDCLU)], ', ', nrow(rtc)-length(unique(rtc$ID)),
			'individuals appear in more than one clusters','\n')
  # 808  clusters of  2364  individuals, the ten largest cluster sizes are  17 18 21 22 26 37 38 45 50 55 ,  150 individuals appear in more than one clusters 
	

	#	add couples that are not a close pair and redefine IDCLU
	load(infile.couple)
	pty.runs <- readRDS(paste0(out.base,'phscinput_samples.rds'))
	rp = data.table(unique(subset(coupdat,select=c('male.RCCS_studyid', 'female.RCCS_studyid'))))
	setnames(rp, c('male.RCCS_studyid', 'female.RCCS_studyid'), c('MALE_RID','FEMALE_RID'))
	rp[,FEMALE_RID:=paste0('RK-',FEMALE_RID)]
	rp[,MALE_RID:=paste0('RK-',MALE_RID)]
	tmp<- sort(unique(as.character(pty.runs[['UNIT_ID']])))	
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
	rtc<- unique(rbind(subset(rtc,select=c('ID','IDCLU')), rp))
	tmp<- rtc[, list(CLU_SIZE=length(ID)), by='IDCLU']
	setkey(tmp, CLU_SIZE)
	tmp[, IDCLU2:=seq_len(nrow(tmp))]	
	rtc<- subset(merge(rtc, tmp, by='IDCLU'),select=c('ID','IDCLU2','CLU_SIZE'))
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
	stopifnot(sum(tmp[,cumsum(PTY_SIZE),by='PTY_RUN2']$V1 !=tmp$PTY_SIZE2)==0)
	tmp[,PTY_SIZE2:=max(PTY_SIZE2),by='PTY_RUN2']
	tmp[,PTY_RUN2:=PTY_RUN2+1]
	set(tmp,NULL,c('PTY_RUN', 'PTY_SIZE'),NULL)
	setnames(tmp,c('PTY_RUN2', 'PTY_SIZE2'), c('PTY_RUN', 'PTY_SIZE'))
	rtc = merge(rtc,tmp,by=c('IDCLU', 'CLU_SIZE'))
	
  #	for all individuals in a cluster
  #	find 3 closest others
	
	if(file.exists(infile.clostest.three.per.cluster)){
	  dcl =readRDS(infile.clostest.three.per.cluster)
	}else{
	  infile.average.distance.per.pair <- file.path(potential.networks.analysis.dir, 'average_distance_perpair.rds')
	  if(file.exists(infile.average.distance.per.pair)){
	    ddist = readRDS(infile.average.distance.per.pair)
	  }else{
	    infile.average.distance.per.pair.per.window <- file.path(potential.networks.analysis.dir, 'distance_overall_method2_sequence_level_v4.rda')
	    load(infile.average.distance.per.pair.per.window)
	    ddist = subset(distance,select=c('TAXA1','TAXA2'))
	    ddist$DIST = rowMeans(distance[,3:ncol(distance)],na.rm = T)
	    setnames(id.dt, colnames(id.dt),paste0(colnames(id.dt),'1'))
	    ddist = merge(ddist, id.dt, by.x='TAXA1', by.y='PANGEA_ID1',all.x=T)
	    setnames(id.dt, colnames(id.dt),gsub('1','2',colnames(id.dt)))
	    ddist = merge(ddist, id.dt, by.x='TAXA2', by.y='PANGEA_ID2',all.x=T)
	    setnames(id.dt, colnames(id.dt),gsub('2','',colnames(id.dt)))
	    ddist <- ddist[,list(DIST=mean(DIST)),by=c('UNIT_ID1','UNIT_ID2')]
	    saveRDS(ddist,file=infile.average.distance.per.pair) 
	  }
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
	  saveRDS(dcl,file=infile.clostest.three.per.cluster)
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
	
	# write processed samples 
	cat('\nWriting to file ', paste0(out.base,'phscinput_runs.rds') )
	saveRDS(rtc, file=paste0(out.base,'phscinput_runs.rds'))
		
}



