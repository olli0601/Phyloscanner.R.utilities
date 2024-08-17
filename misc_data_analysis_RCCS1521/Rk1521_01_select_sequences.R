rk.seq.get.mappeable.samples.Rakai <- function()
{
	#module purge; module add tools/prod; module add R/4.3.2-gfbf-2023a
	require(data.table)
	indir.rk <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS'
	outdir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live'
		
	pangea.sharing.extract.1 <- file.path(indir.rk, '200316_pangea_db_sharing_extract_rakai.csv')
	pangea.sharing.extract.2 <- file.path(indir.rk, '240709_pangea2_extract_longview_olli.csv')
	pangea.files2 <- file.path(indir.rk, '200316_pangea_mappings_rakai.csv')	
	pangea.filesm <- file.path(indir.rk, '200316_unlinkable_files_mappings.csv')
	
	#
	#	read RCCS info from Oxford
	#
	pf <- as.data.table(read.csv(pangea.sharing.extract.1, stringsAsFactors=FALSE))
	pf <- subset(pf, select=c(pt_id, geo_country, cohort_id, main_cohort_id, visit_dt, pangea_id, readnum, readnum_hiv, mapped_num,  duprate, insertfrac_350, lc_cutoff, n_cutoff, length_strict, subtype_bestref))
	pf <- subset(pf, select=c(pt_id, geo_country, cohort_id, main_cohort_id, visit_dt, pangea_id, subtype_bestref))
	setnames(pf, colnames(pf), toupper(colnames(pf)))
	pf <- unique(pf)
	# 7567
	# pf <- subset(pf, GEO_COUNTRY%in%c('Uganda','Masaka'))
		
	#	read sequence IDs and path to data from Oxford
	tmp <- as.data.table(read.csv(pangea.files2, stringsAsFactors=FALSE))
	tmp2 <- as.data.table(read.csv(pangea.filesm, stringsAsFactors=FALSE))
	setnames(tmp, colnames(tmp), toupper(colnames(tmp)))
	set(tmp, NULL, 'DEPRECATED', 'N')
	set(tmp, NULL, 'LINKABLE', 'Y')
	setnames(tmp2, colnames(tmp2), toupper(colnames(tmp2)))
	set(tmp2, NULL, c('CURRENT_SEQUENCE_ID','CURRENT_PATH_TO_DATA'), NULL)
	set(tmp2, NULL, 'LINKABLE', 'N')
	tmp <- rbind(tmp, tmp2)	
	#	update some entries after checking with Laura Thomson, 200318
	tmp2 <- tmp[, which(grepl('PG15-UG504872',PANGEA_ID))]
	stopifnot(length(tmp2)==1)
	set(tmp, tmp2, 'SEQUENCE_ID', 'UG504872_PanRak10-E1A1L1P1S1C1R1')
	set(tmp, tmp2, 'LINKABLE', 'Y')
	set(tmp, tmp2, 'DEPRECATED', 'N')
	tmp <- rbind(tmp, data.table(PANGEA_ID='PG15-UG504872', SEQUENCE_ID='UG504872_PanRak12-E1A1L1P1S1C1R1', LINKABLE='Y', DEPRECATED='Y'), fill=TRUE)	
	tmp2 <- tmp[, which(grepl('PG15-UG504959',PANGEA_ID))]
	stopifnot(length(tmp2)==1)
	set(tmp, tmp2, 'SEQUENCE_ID', 'UG504959_PanRak13-E1A1L1P1S1C1R1')
	set(tmp, tmp2, 'LINKABLE', 'Y')
	set(tmp, tmp2, 'DEPRECATED', 'N')	
	tmp2 <- tmp[, which(grepl('UG504070',SEQUENCE_ID))]
	stopifnot(length(tmp2)==1)
	set(tmp, tmp2, 'PANGEA_ID', 'PG15-UG504070')
	set(tmp, tmp2, 'SEQUENCE_ID', 'UG504070_PanRak7-E1A1L1P1S1C1R1')
	set(tmp, tmp2, 'LINKABLE', 'Y')
	set(tmp, tmp2, 'DEPRECATED', 'N')	
	tmp2 <- tmp[, which(grepl('UG504463',SEQUENCE_ID))]
	stopifnot(length(tmp2)==0)
	tmp <- rbind(tmp, data.table(PANGEA_ID='PG15-UG504463', SEQUENCE_ID='UG504463_PanRak6-E1A1L1P1S0C1R1', LINKABLE='Y', DEPRECATED='N'), fill=TRUE)	
	tmp2 <- tmp[, which(grepl('UG505394',SEQUENCE_ID))]
	stopifnot(length(tmp2)==0)
	tmp <- rbind(tmp, data.table(PANGEA_ID='PG15-UG505394', SEQUENCE_ID='UG505394_PanRak8-E1A1L1P1S1C1R1', LINKABLE='Y', DEPRECATED='N'), fill=TRUE)
	#	
	set(tmp, tmp[, which(SEQUENCE_ID=='NULL')], 'SEQUENCE_ID', NA_character_)
	set(tmp, tmp[, which(PATH_TO_DATA=='NULL')], 'PATH_TO_DATA', NA_character_)	
	tmp[, BATCH_ID:= gsub('.*(PANGEA1_miseq|PANGEA1_hiseq|PanMix1).*','\\1',basename(PATH_TO_DATA))]
	tmp2 <- tmp[, which(grepl('PanRak',SEQUENCE_ID))]
	set(tmp, tmp2, 'BATCH_ID', tmp[tmp2,gsub('.*(PanRak[^-]+).*','\\1',SEQUENCE_ID)])
	# checking 
	unique(subset(tmp, !grepl('PanRak|PANGEA1',BATCH_ID) & BATCH_ID=='TBC'))	
	# 87 are TBC
	unique(subset(tmp, !grepl('PanRak|PANGEA1',BATCH_ID) & grepl('PanMix1',BATCH_ID)))
	# 24 from PanMix1
	unique(subset(tmp, !grepl('PanRak|PANGEA1',BATCH_ID) & is.na(BATCH_ID)))
	# 126 with NA sequence IDs and no data	
	tmp <- unique(tmp)
	# 8359
	# merge data from Oxford
	pf <- merge(tmp, pf, by='PANGEA_ID', all=TRUE)		
	# 8380
	
	#	update a few entries manually, from Laura March 27 2020
	tmp <- pf[, which(grepl('PG15-UG505394',PANGEA_ID))]
	stopifnot(length(tmp)==1)
	set(pf, tmp, 'PT_ID','RK-C068239')
	set(pf, tmp, 'VISIT_DT','2011-02-04')
	set(pf, tmp, 'COHORT_ID','U-JHU-RCCS')
	set(pf, tmp, 'MAIN_COHORT_ID','Rakai')
	set(pf, tmp, 'GEO_COUNTRY','Uganda')
	tmp <- pf[, which(grepl('PG15-UG504872',PANGEA_ID) & DEPRECATED=='N')]
	stopifnot(length(tmp)==1)
	set(pf, tmp, 'PT_ID','RK-D115354')
	set(pf, tmp, 'VISIT_DT','2014-05-14')
	set(pf, tmp, 'COHORT_ID','U-JHU-RCCS')
	set(pf, tmp, 'MAIN_COHORT_ID','Rakai')
	set(pf, tmp, 'GEO_COUNTRY','Uganda')
	tmp <- pf[, which(grepl('PG15-UG504463',PANGEA_ID))]
	stopifnot(length(tmp)==1)
	set(pf, tmp, 'PT_ID','RK-E045395')
	set(pf, tmp, 'VISIT_DT','2015-03-23')
	set(pf, tmp, 'COHORT_ID','U-JHU-RCCS')
	set(pf, tmp, 'MAIN_COHORT_ID','Rakai')
	set(pf, tmp, 'GEO_COUNTRY','Uganda')
	tmp <- pf[, which(grepl('PG15-UG504959',PANGEA_ID))]
	stopifnot(length(tmp)==1)
	set(pf, tmp, 'PT_ID','RK-F102900')
	set(pf, tmp, 'VISIT_DT','2011-06-07')
	set(pf, tmp, 'COHORT_ID','U-JHU-RCCS')
	set(pf, tmp, 'MAIN_COHORT_ID','Rakai')
	set(pf, tmp, 'GEO_COUNTRY','Uganda')
	tmp <- pf[, which(grepl('PG15-UG504070',PANGEA_ID))]
	stopifnot(length(tmp)==1)
	set(pf, tmp, 'PT_ID','RK-J058274')
	set(pf, tmp, 'VISIT_DT','2011-09-07')
	set(pf, tmp, 'COHORT_ID','U-JHU-RCCS')
	set(pf, tmp, 'MAIN_COHORT_ID','Rakai')
	set(pf, tmp, 'GEO_COUNTRY','Uganda')
	
	#
	#	read LONGVIEW sharing extract and add
	#
	tmp <- as.data.table(read.csv(pangea.sharing.extract.2, stringsAsFactors=FALSE))
	tmp <- subset(tmp, select=c(pt_id, geo_country, cohort_id, main_cohort_id, visit_dt, pangea_id, sequence_id))
	setnames(tmp, colnames(tmp), toupper(colnames(tmp)))
	tmp <- unique(tmp)
	tmp[, DEPRECATED := 'N']
	tmp[, LINKABLE := 'Y']
	tmp[, BATCH_ID := 'Longview']
	set(tmp, NULL, 'COHORT_ID', "U-JHU-RCCS")
	pf <- rbindlist(list(pf, tmp), fill = TRUE )
	pf <- subset(pf, !is.na(VISIT_DT))
	
	saveRDS(pf, file=file.path(outdir, '240809_pangea_db_sharing_extract_rakai_combined.rds'))
	
	#
	#	read RCCS data at Imperial
	#
	ds <- data.table(F=list.files(indir.rk, recursive=TRUE, full.names=TRUE))
	#	82961
	

	#save(pf, ds, file= file.path(outdir, 'PANGEA2_RCCS_interim.rda'))
	#load(file.path(outdir, 'PANGEA2_RCCS_interim.rda'))
	
	#
	#	process file ends
	tmp <- '.*(_ref.fasta|.bam|.bam.bai|_ForGlobalAln.fasta)$'
	ds[, FEND:= gsub(tmp,'\\1',F)]
	ds <- subset(ds, grepl(tmp, FEND))
	ds[, PREDEDUP:= as.character(factor(grepl('PreDedup',basename(F)), levels=c(TRUE,FALSE), labels=c('_PreDeDup','')))]
	ds[, REMAP:= as.character(factor(grepl('remap',basename(F)), levels=c(TRUE,FALSE), labels=c('_remap','')))]
	tmp <- ds[, which(grepl('_consensus_',basename(F)))]
	ds[, CONSENSUS:= '']
	set(ds, tmp, 'CONSENSUS', ds[tmp, gsub('.*(_consensus_MinCov_[0-9]+_[0-9]+_ForGlobalAln.fasta).*','\\1',basename(F))])
	ds[, FEND2:= paste0(REMAP,PREDEDUP,FEND)]
	tmp <- ds[, which(nchar(CONSENSUS)>0)]
	set(ds, tmp, 'FEND2', ds[tmp,CONSENSUS])
	
	#
	#	read batch IDs
	ds[, BATCH_ID:= basename(dirname(F))]	
	tmp <- ds[, which(BATCH_ID=='shiver-output')]
	set(ds, tmp, 'BATCH_ID', ds[tmp, basename(dirname(dirname(F)))])
	ds[, table(BATCH_ID, useNA='if')]
	
	#
	#	read sequencing IDs
	ds[, ROW_ID:= seq_len(nrow(ds))]
	set(ds, NULL, c('PREFIX','SEQUENCE_ID','PID'), NULL)
	ds2 <- ds[, list(PREFIX=substr(basename(F), 1, nchar(basename(F))-nchar(FEND2))) , by='ROW_ID']
	#	remove controls
	tmp <- 'Neg|neg|^SC|^ctrl|^100$|^500$|^5000$|^50000$|^500000$|^5000000$|NITP'
	ds2 <- subset(ds2, !grepl(tmp, PREFIX))
	#	check for old Sanger IDs of format 13549_1_84 and set SEQUENCE_ID accordingly
	tmp <- ds2[, which(grepl('.*([0-9]{5}_[0-9]_[0-9]+).*', PREFIX))]
	set(ds2, tmp, 'SEQUENCE_ID', ds2[tmp, gsub('.*([0-9]{5}_[0-9]_[0-9]+).*','\\1',PREFIX)])
	tmp <- ds2[, which(is.na(SEQUENCE_ID))]
	set(ds2, tmp, 'SEQUENCE_ID', ds2[tmp, PREFIX])
	#	build pre-PANGEA ID, this seems to be a mix of PANGEA IDs and Rakai IDs		 
	ds2[, PID:= gsub('.*([A-Z]{2}[0-9]{6}).*', '\\1', PREFIX)]
	tmp <- '.*([A-Z][0-9]{6}R[0-9]{2}).*'
	tmp2 <- ds2[, which(grepl(tmp, PREFIX))]
	set(ds2, tmp2, 'PID', ds2[tmp2,gsub(tmp,'\\1',PREFIX)])
	tmp <- '.*([A-Z][0-9]{6}NEU).*'
	tmp2 <- ds2[, which(grepl(tmp, PREFIX))]
	set(ds2, tmp2, 'PID', ds2[tmp2,gsub(tmp,'\\1',PREFIX)])
	ds <- merge(ds, ds2, by='ROW_ID') # this removes controls
	# 92434
		
	#
	#	prepare linking to PANGEA_IDs from Tanya
	pf2 <- unique(subset(pf, select=c(PANGEA_ID, BATCH_ID, SEQUENCE_ID, PT_ID, VISIT_DT, DEPRECATED, LINKABLE)))
	#	9986	
	ds2 <- unique(subset(ds, select=c(PID, BATCH_ID, SEQUENCE_ID, PREFIX)))
	#	10718

	ds3 <- merge(ds2, pf2, by=c('SEQUENCE_ID','BATCH_ID'), all=TRUE)
	dc <- subset(ds3, !is.na(PT_ID) & !is.na(PID))
	#	matched 9754 entries (this is exactly what we could previously match at this stage plus the new Longview seqs, so all looks good)
	
	ds3 <- subset(ds3, !(!is.na(PT_ID) & !is.na(PID)) )
	ds3 <- subset(ds3, DEPRECATED=='N' | is.na(DEPRECATED) )
	#	could not match 1593 entries

	
	#
	#	manage PANGEA1_miseq / PANGEA1_hiseq with no Oxford link
	#
	dm <- subset(ds3, grepl('PANGEA1_miseq|PANGEA1_hiseq',BATCH_ID) & is.na(PT_ID))
	#	try to isolate full PANGEA ID
	tmp <- '.*(PG[0-9]{2}-[A-Z]{2}[0-9]{6}).*'
	tmp2 <- dm[, which(grepl(tmp, PREFIX))]
	set(dm, tmp2, 'PANGEA_ID', dm[tmp2, gsub(tmp, '\\1', PREFIX)])
	stopifnot( nrow(subset(dm, is.na(PANGEA_ID))) == 0 )
	#	success for all!
	dc <- rbind(dc, dm)
	ds3 <- subset(ds3, !grepl('PANGEA1_miseq|PANGEA1_hiseq',BATCH_ID))
	#	still cannot match 322 entries
	

	#
	#	PIDs with NEU are not from general pop, ignore
	#
	dm <- subset(ds3, grepl('NEU',PID))
	dm[, PT_ID:= paste0('RK-',gsub('.*([A-Z][0-9]{6})NEU.*','\\1',PID))]
	dm <- merge(subset(dm, select=PT_ID), subset(pf, !is.na(PT_ID), select=c(PT_ID, PANGEA_ID, BATCH_ID,SEQUENCE_ID)), by='PT_ID', all.x=TRUE)	
	#write.csv(file=file.path(outdir, '200318_RCCS2_PANGEA_resolvedPANGEAIDs_NEU.csv'), subset(dm, !is.na(PANGEA_ID)), row.names=FALSE)
	#write.csv(file=file.path(outdir, '200318_RCCS2_PANGEA_missing_data_NEU.csv'), subset(dm, is.na(PANGEA_ID)), row.names=FALSE)
	#	success for all except 1!
	#	report to Laura	
	ds3 <- subset(ds3, !grepl('NEU',PID))
	#	still cannot match 138 entries
	
	#
	#	check remaining entries
	#
	dm <- subset(ds3, !is.na(BATCH_ID) & is.na(PT_ID))
	stopifnot(nrow(dm)==0)
	#	success!
	
		
	#
	#	retain samples that can be mapped to PANGEA IDs
	#	
	tmp <- unique(subset(dc, select=c(PREFIX, BATCH_ID, PANGEA_ID)))
	ds <- merge(tmp, ds, by=c('PREFIX','BATCH_ID'), all.x=TRUE)
	#	90725
	saveRDS(ds, file=file.path(outdir, '240809_PANGEA2_RCCS_mapped_samples.rds'))	
}

rk.seq.get.best.mappeable.samples.Rakai <- function()
{
	require(data.table)
	require(ape)
	indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live'
	indir.rk <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS'
	indir.mrc <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC'
	
		
	ds <- readRDS( file.path(indir, '240809_PANGEA2_RCCS_mapped_samples.rds') )
	ds[, length(unique(PANGEA_ID))]
	#	9242 mappeable samples 
	
	#
	#	find samples who also have a consensus sequence
	#
	tmp <- unique(subset(ds, CONSENSUS!='', c(PANGEA_ID, BATCH_ID, PREFIX)))
	tmp[, length(unique(PANGEA_ID))]
	#	9089 mappeable samples with consensus sequence, keep only those
	ds <- merge(tmp, ds, by=c('PANGEA_ID','PREFIX','BATCH_ID'))
	
	
	#
	#	select one consensus sequence for each sample
	#
	tmp <- ds[CONSENSUS!='', list(N=length(CONSENSUS)), by='PANGEA_ID']
	# > tmp[, table(N)]
	#N
	#   1    2    3
	#7967  907  215	
	ds2 <- merge(subset(tmp, N==1,PANGEA_ID), ds, by='PANGEA_ID')	
	ds3 <- merge( subset(tmp, N>1,PANGEA_ID), ds[CONSENSUS!='',], by='PANGEA_ID')
	ds4 <- merge( subset(tmp, N>1,PANGEA_ID), ds, by='PANGEA_ID')
	#	for high seq, ignore anything below 150_300
	ds3 <- rbind( subset(ds3, grepl('hiseq',F) & grepl('MinCov_150_300',F)), subset(ds3, !grepl('hiseq',F)))
	#	2459 consensus files
	#	keep PanRak over hiseq, and hiseq over miseq
	#	for Pan2, keep 5_15 over 15_30
	ds3[, SEQ_TYPE:= NA_integer_]	
	set(ds3, ds3[, which(grepl('PanRak|PanMix',F) & grepl('MinCov_5_15',F))],'SEQ_TYPE',4L)	#order is important, since miseq also has PanRak in name
	set(ds3, ds3[, which(grepl('PanRak|PanMix',F) & grepl('MinCov_15_30',F))],'SEQ_TYPE',3L)
	set(ds3, ds3[, which(grepl('miseq',F) & grepl('MinCov_5_15',F))],'SEQ_TYPE',1L)
	set(ds3, ds3[, which(grepl('miseq',F) & grepl('MinCov_15_30',F))],'SEQ_TYPE',0L)
	set(ds3, ds3[, which(grepl('hiseq',F))],'SEQ_TYPE',2L)	
	tmp <- ds3[, list(SEQ_TYPE=max(SEQ_TYPE)), by=c('PANGEA_ID')]	 
	ds3 <- merge(ds3, tmp, by=c('PANGEA_ID','SEQ_TYPE'))
	ds4 <- rbind( ds4[CONSENSUS=='',], subset(ds3, select=-SEQ_TYPE) )
	#	add chosen samples to ds2
	tmp <- ds3[, list(N=length(CONSENSUS)), by='PANGEA_ID']
	tmp2 <- merge(subset(tmp, N==1,PANGEA_ID), ds3, by='PANGEA_ID')
	tmp2 <- merge( subset(tmp2, select=c(PANGEA_ID, PREFIX, BATCH_ID)), ds4, by=c('PANGEA_ID','PREFIX','BATCH_ID') )
	ds2 <- rbind(ds2, tmp2)
	#	7045 files corresponding to PANGEA_ID with one selected consensus
	ds3 <- merge(subset(tmp, N>1,PANGEA_ID), ds3, by='PANGEA_ID')
	#	1877 files corresponding to PANGEA_IDs with multiple consensus
	#	remaining files are re-sequenced, keep the one with longest consensus
	tmp <- ds3[, list(N_ACTG=nchar(gsub('\\?|\\-','',paste0(as.character(read.dna(F, format='fasta')), collapse='')))), by=c('PANGEA_ID','PREFIX','BATCH_ID')]
	tmp <- tmp[, {
				z<- which.max(N_ACTG)[1]
				list(BATCH_ID= BATCH_ID[z], PREFIX=PREFIX[z])
			}, by='PANGEA_ID']
	tmp2 <- merge( tmp, ds4, by=c('PANGEA_ID','PREFIX','BATCH_ID') )
	ds <- rbind(ds2, tmp2)	
	#	78511 files corresponding to PANGEA_ID with one selected consensus
	#	double-check
	tmp <- subset(ds, CONSENSUS!='', c(PANGEA_ID, PREFIX, BATCH_ID))
	tmp2 <- tmp[, list(N=length(PREFIX)), by='PANGEA_ID']
	stopifnot( nrow(subset(tmp2, N>1))==0 )
	#	success
	

	#
	#	bring into wide format
	#
	ds[, FEND2:= toupper(gsub('^_','',gsub('\\.','_',paste0(REMAP,PREDEDUP,FEND))))]
	tmp <- ds[, list(N=length(F)), by=c('PANGEA_ID','FEND2')]
	stopifnot( nrow(subset(tmp, N>1))==0 )
	ds <- dcast.data.table(ds, PANGEA_ID ~ FEND2, value.var='F')
	#	9089 samples with consensus found
	

	#	add length of consensus
	tmp <- ds[, 
		{
			#print(FORGLOBALALN_FASTA)
			z <- read.dna(FORGLOBALALN_FASTA, format='fasta')
			z <- gsub('a|c|t|g|A|C|T|G','x',paste0(as.character(z)))
			z <- data.table(NUC = as.numeric( z == 'x' ) )
			z[, RUN_ID := rleid(NUC) ] 
			z <- merge(z, z[, list(RUN_LEN = length(NUC)), by = 'RUN_ID'], by = 'RUN_ID')
			z[, NUC_2 := NUC]
			set(z, z[, which( NUC == 0 & RUN_LEN < 2)], 'NUC_2', 1 )
			z[, RUN_ID_2 := rleid(NUC_2) ]
			z2 <- z[ NUC_2 == 1, list(RUN_LEN = length(NUC_2)), by = 'RUN_ID_2']
			
			list( N_NUC = sum(z[, NUC]), 
				  MX_RUN_GRACE_1 = ifelse( nrow(z2) == 0, 0L, max(z2[, RUN_LEN]) ) 
				  )			
		}, 
		by=c('PANGEA_ID')
		]
	ds <- merge(ds, tmp, by = 'PANGEA_ID')
	
		
	# add RCCS study IDs
	pf <- readRDS(file.path(indir, '240809_pangea_db_sharing_extract_rakai_combined.rds'))
	tmp <- unique(subset(pf, select = c(PANGEA_ID, VISIT_DT, PT_ID)))
	ds <- merge(ds, tmp, by = 'PANGEA_ID', all.x = TRUE)	
	write.csv(subset(ds, is.na(PT_ID)), file=file.path(indir, '240809_PANGEA2_RCCS_selected_samples_missing_study_id.csv'))
	saveRDS(subset(ds, is.na(PT_ID)), file=file.path(indir, '240809_PANGEA2_RCCS_selected_samples_missing_study_id.rds'))
	ds <- subset(ds, !is.na(PT_ID))
	set(ds, NULL, 'VISIT_DT', ds[, as.Date(VISIT_DT)])
	setkey(ds, PT_ID, VISIT_DT)
	
	# check multiple samples for each visit
	tmp <- ds[,
		list( N_SEQ = length(FORGLOBALALN_FASTA)), 		 
		by = c('PT_ID', 'VISIT_DT')
		]
	tmp <- merge(ds, subset(tmp, N_SEQ > 1), by = c('PT_ID','VISIT_DT'))
	write.csv(tmp, file = file.path(indir, '240809_PANGEA2_RCCS_selected_samples_multiple_samples_per_visit.csv'))
	
	# for each RCCS individual and each visit, keep sample with longest run
	tmp <- ds[,
		list( N_SEQ = length(FORGLOBALALN_FASTA),  N_NUC = max(N_NUC), MX_RUN_GRACE_1 = max(MX_RUN_GRACE_1)), 		 
		by = c('PT_ID', 'VISIT_DT')
		]
	ds <- merge(ds, subset(tmp, select = -c(N_SEQ, N_NUC)), by = c('PT_ID', 'VISIT_DT', 'MX_RUN_GRACE_1' ))
	# next, for ties, for each RCCS individual and each visit, keep sample with longest nuc
	tmp <- ds[,
		list( N_SEQ = length(FORGLOBALALN_FASTA),  N_NUC = max(N_NUC), MX_RUN_GRACE_1 = max(MX_RUN_GRACE_1)), 		 
		by = c('PT_ID', 'VISIT_DT')
		]
	ds <- merge(ds, subset(tmp, select = -c(N_SEQ, MX_RUN_GRACE_1)), by = c('PT_ID', 'VISIT_DT', 'N_NUC' ))
	# next, for ties, keep one at random
	# TODO is there an error in visit date?
	#subset(ds, PT_ID == 'RK-J053740' & VISIT_DT == '2015-02-26' )
	#Key: <PT_ID, VISIT_DT, N_NUC>
	#        PT_ID   VISIT_DT N_NUC MX_RUN_GRACE_1     PANGEA_ID
	#       <char>     <Date> <num>          <int>        <char>
	#1: RK-J053740 2015-02-26  8968           2989 PG15-UG504451
	#2: RK-J053740 2015-02-26  8968           2989 PG19-UG001719
	ds <- subset(ds, !(PT_ID == 'RK-J053740' & VISIT_DT == '2015-02-26' & PANGEA_ID == 'PG19-UG001719'))
	
	# final check
	tmp <- ds[,
		list( N_SEQ = length(FORGLOBALALN_FASTA)), 		 
		by = c('PT_ID', 'VISIT_DT')
		]
	stopifnot( nrow(subset(tmp, N_SEQ > 1)) == 0 )
	
	# determine number of seqs per person	
	tmp <- ds[, list(VISIT_DT= VISIT_DT, PT_ID_SEQ_NO = seq_along(VISIT_DT)), by = 'PT_ID']
	ds <- merge(ds, tmp, by = c('PT_ID','VISIT_DT'))
	
	ds[, table(cut(N_NUC, c(-1,250,500,1000,20000)))]
	#  (-1,250]     (250,500]   (500,1e+03] (1e+03,2e+04]
    #      996           296           299          7193
         
	ds[, table(cut(MX_RUN_GRACE_1, c(-1,125, 250,500,1000,20000)))]
	#  (-1,125]     (125,250]     (250,500]   (500,1e+03] (1e+03,2e+04]
    #     1037           468          1719           502          5058
	
	# 8784 sequences 
	saveRDS(ds, file=file.path(indir, '240809_PANGEA2_RCCS_final_samples.rds'))
	
	ds2 <- subset(ds, !(MX_RUN_GRACE_1 < 125 & N_NUC < 250) )
	tmp <- ds2[, list(VISIT_DT= VISIT_DT, PT_ID_SEQ_NO = seq_along(VISIT_DT)), by = 'PT_ID']
	ds2 <- merge(subset(ds2, select = -PT_ID_SEQ_NO ), tmp, by = c('PT_ID','VISIT_DT'))
		
	# 7836 sequences 
	saveRDS(ds2, file=file.path(indir, '240809_PANGEA2_RCCS_final_samples_minrun125_minnuc250.rds'))
}

rk.seq.remove.poor.samples.Rakai <- function()
{
	require(data.table)
	require(ape)
	
	indir <- '~/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/PANGEA/2024_PANGEA_Rakai_analyses'
	file <- '240809_PANGEA2_RCCS_final_samples_minrun125_minnuc250.rds'
	ds <- readRDS(file = file.path(indir, file))
	
	# remove poor quality samples if there are multiple sequences per individual
	tmp <- ds[,
		list( N_SEQ = length(FORGLOBALALN_FASTA),
			  CONSIDER_SEQ = N_NUC < 750 & MX_RUN_GRACE_1 < 500,
			  CONSIDER_ALLSEQ = all(N_NUC < 750 & MX_RUN_GRACE_1 < 500),
			  CONSIDER_IND = any(N_NUC < 750 & MX_RUN_GRACE_1 < 500),
			  PANGEA_ID = PANGEA_ID
			  ),			   		 
		by = c('PT_ID')
		]
	tmp <- merge(ds, subset(tmp, N_SEQ > 1), by = c('PT_ID','PANGEA_ID'))
	tmp <- subset(tmp, CONSIDER_IND)
	
	# for manual review
	write.csv(subset(tmp,  CONSIDER_ALLSEQ), file = file.path(indir, '240809_PANGEA2_RCCS_final_samples_minrun125_minnuc250_all_seqs_of_ind_poor.csv'))
	
	# remove all of these
	tmp <- subset(tmp, CONSIDER_IND & !CONSIDER_ALLSEQ & CONSIDER_SEQ)
	tmp[, REMOVE := TRUE]
	ds <- merge(ds, subset(tmp, select = c(PT_ID, PANGEA_ID, REMOVE)), by = c('PT_ID','PANGEA_ID'), all.x = TRUE )
	ds <- subset(ds, is.na(REMOVE), select = -c(REMOVE, PT_ID_SEQ_NO))
	
	# rebuild seq no
	tmp <- ds[, 
				list(VISIT_DT= VISIT_DT, 
					 N_SEQ = length(FORGLOBALALN_FASTA),
					 PT_ID_SEQ_NO = seq_along(VISIT_DT)
					 ), 
				by = 'PT_ID']
	ds <- merge(ds, tmp, by = c('PT_ID','VISIT_DT'))	
	
	# 7691 sequences
	saveRDS(ds, file=file.path(indir, '240809_PANGEA2_RCCS_final_samples_minrun125_minnuc250_minusPoor1.rds'))
}

rkuvri.make.anonymised.id <- function()
{   
  #' input person ids, return ananymise ids

  # file names
  outfile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521_UVRI/important_anonymisation_keys_240817.csv'
  
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  infile.ind.mrc <- file.path(data.dir,'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
  infile.ind.rccs <- file.path(data.dir,'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.rccs2 <- file.path(data.dir, 'PANGEA2_RCCS/240709_pangea2_extract_longview_olli.csv')

  infile.ind.anonymised.2021 <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'	
	

  # load person IDs
  id.dt <- data.table(read.csv(infile.ind.rccs))
  id.dt <- unique(subset(id.dt,select = c("pt_id")))
  
  tmp <- data.table(read.csv(infile.ind.rccs2))
  tmp <- subset(tmp,select = c("pt_id"))
  id.dt <- unique(rbind(id.dt, tmp))
  
  tmp <- data.table(read.csv(infile.ind.mrc))
  tmp <- unique(subset(tmp,select = c("pt_id")))
  id.dt <- rbind(id.dt,tmp)
  
  setkey(id.dt,pt_id)
  setnames(id.dt,colnames(id.dt),toupper(colnames(id.dt)))

  # load previous anonymised IDs
  tmp <- as.data.table(read.csv( file = infile.ind.anonymised.2021 ))
  set(tmp, NULL, "X", NULL)
  set(tmp, NULL, "AID_NO", tmp[, as.integer(gsub('AID','',AID))])
  id.dt <- merge(id.dt, tmp, by = "PT_ID", all.x = TRUE)

  # set seed and sample identifier at random
  set.seed(42)
  tmp <- subset(id.dt, is.na(AID))
  tmp[,AID_NO:=sample(1:nrow(tmp) + id.dt[, max(AID_NO, na.rm = TRUE)],replace = F)]
  tmp[,AID:=formatC(AID_NO, width=floor(log10(nrow(tmp))) +1, flag="0")]
  tmp[,AID:=paste0('AID',AID)]
  id.dt <- rbind(subset(id.dt, !is.na(AID)), tmp)
  id.dt <- subset(id.dt, select = - AID_NO)
  write.csv(id.dt,file = outfile.ind.anonymised)
  Sys.chmod(outfile.ind.anonymised,mode = '0444')
}

rkuvri.seq.make.consensus.alignment <- function()
{
	# TODO update
	require(data.table)
	require(ape)
	indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS'
	infile.rk <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/200422_PANGEA2_RCCS_selected_samples.rds'
	infile.mrc <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC/200422_PANGEA2_MRCUVRI_selected_samples.rds'
	outfile <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/200422_PANGEA2_RCCSMRC_alignment.fasta'
	
	#	add consensus sequences to alignment
	file.remove(outfile)
	dr <- readRDS(infile.rk)
	invisible( dr[, {
				z <- paste0('cat ',FORGLOBALALN_FASTA,' >> ', outfile)
				system(z)
				NULL
			}, by='PANGEA_ID'])
	dm <- readRDS(infile.mrc)
	invisible( dm[, {
						z <- paste0('cat ',FORGLOBALALN_FASTA,' >> ', outfile)
						system(z)
						NULL
					}, by='PANGEA_ID'])	
	
	#	new taxa names
	ds <- read.dna(outfile, format='fasta')
	dr[, LABEL:= paste0('RCCS_',PANGEA_ID)]
	dm[, LABEL:= paste0('MRCUVRI_',PANGEA_ID)]	
	rownames(ds) <- c(dr$LABEL, dm$LABEL)
	write.dna(ds, file=outfile, format='fa', colsep='', nbcol=-1)
}
