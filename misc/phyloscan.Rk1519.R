rk.seq.delete.files.dir.other <- function()
{
	require(data.table)
	indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/other'
	infile <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/200317_RCCS2_PANGEA_delete_from_others.csv'
	
	ds <- data.table(F=list.files(indir, recursive=TRUE, full.names=TRUE))
	tmp <- '(.*)(_ref.fasta|.bam|.bam.bai|_ForGlobalAln.fasta)$'
	ds[, PREFIX:= gsub('_PreDedup','',gsub('_remap','',gsub(tmp,'\\1',basename(F))))]
		
	tmp <- as.data.table(read.csv(infile, stringsAsFactors=FALSE))
	tmp[, DELETE:= TRUE]
	setnames(tmp, 'PID', 'PREFIX')	
	ds <- merge(ds, tmp, by='PREFIX', all.x=TRUE)
	subset(ds, is.na(DELETE))	
	#	only controls left, so we can delete the whole thing
}

rk.seq.delete.files.dir.PanMix1 <- function()
{
	require(data.table)
	indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/PanMix1'
	infile <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/200317_RCCS2_PANGEA_delete_from_PanMix1.csv'
	
	ds <- data.table(F=list.files(indir, recursive=TRUE, full.names=TRUE))
	tmp <- '(.*)(_ref.fasta|.bam|.bam.bai|_ForGlobalAln.fasta)$'
	ds[, PREFIX:= gsub('_PreDedup','',gsub('_remap','',gsub(tmp,'\\1',basename(F))))]
	
	tmp <- as.data.table(read.csv(infile, stringsAsFactors=FALSE))
	tmp[, DELETE:= TRUE]
	ds <- merge(ds, tmp, by='PREFIX', all.x=TRUE)
	tmp <- subset(ds, DELETE)
	tmp[, file.remove(F), by='F']
	
	subset(ds, is.na(DELETE))	
	#	done
}


rk.seq.copy.MRC.samples.in.PANGEA.RCCS.folders <- function()
{
	require(data.table)
	indir.mrc <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC'
	indir.pan1.mi <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/PANGEA1_miseq/shiver-output'
	indir.pan1.hi <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/PANGEA1_hiseq/shiver-output'
	
	
	
	#
	#	read data at Imperial in MRC folder
	#
	ds <- data.table(F=list.files(indir.mrc, recursive=TRUE, full.names=TRUE))
	tmp <- '^(PG[0-9]{2}-UG[0-9]{6}).*'
	stopifnot( nrow(subset(ds, grepl(tmp, basename(F)))) == nrow(ds) )
	ds[, PANGEA_ID:= gsub(tmp,'\\1',basename(F))]
	setnames(ds, 'F', 'F_MRC')
	#	22879
	
	#
	#	read data at Imperial in PANGEA1 folders pulled from Oxford
	#	
	ds2 <- data.table(F=list.files(indir.pan1.mi, recursive=TRUE, full.names=TRUE))
	tmp <- data.table(F=list.files(indir.pan1.hi, recursive=TRUE, full.names=TRUE))
	ds2 <- rbind(ds2, tmp)
	tmp <- '^(PG[0-9]{2}-UG[0-9]{6}).*'
	stopifnot( nrow(subset(ds2, grepl(tmp, basename(F)))) == nrow(ds2) )
	ds2[, PANGEA_ID:= gsub(tmp,'\\1',basename(F))]
	setnames(ds2, 'F', 'F_PAN')
	#	49206 files
	
	#
	#	get files that should be MRC, and copy over from Rakai directory into MRC directory
	#
	ds3 <- merge(unique(subset(ds, select=PANGEA_ID)), ds2, by='PANGEA_ID')
	#	0
	ds3[, F_NEW:= file.path(indir.mrc, basename(F_PAN))]
	ds3[, file.copy(F_PAN, F_NEW, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE), by='F_PAN']
	ds3[, ID:= 1:nrow(ds3)]
	ds3[, file.rename(F_PAN, F_NEW), by='ID']
	#	done
}

rk.seq.make.consensus.alignment <- function()
{
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

rk.seq.get.best.mappeable.samples.Rakai <- function()
{
	require(data.table)
	require(ape)
	indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live'
	indir.rk <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS'
	indir.mrc <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC'
	
		
	ds <- readRDS( file.path(indir, '200422_PANGEA2_RCCS_mapped_samples.rds') )
	ds[, length(unique(PANGEA_ID))]
	#	7615 mappeable samples 
	
	#
	#	find samples who also have a consensus sequence
	#
	tmp <- unique(subset(ds, CONSENSUS!='', c(PANGEA_ID, BATCH_ID, PREFIX)))
	tmp[, length(unique(PANGEA_ID))]
	#	7462 mappeable samples with consensus sequence, keep only those
	ds <- merge(tmp, ds, by=c('PANGEA_ID','PREFIX','BATCH_ID'))
	
	
	#
	#	select one consensus sequence for each sample
	#
	tmp <- ds[CONSENSUS!='', list(N=length(CONSENSUS)), by='PANGEA_ID']	
	ds2 <- merge(subset(tmp, N==1,PANGEA_ID), ds, by='PANGEA_ID')
	#	57032 files corresponding to PANGEA_ID with one available consensus	
	ds3 <- merge( subset(tmp, N>1,PANGEA_ID), ds[CONSENSUS!='',], by='PANGEA_ID')
	ds4 <- merge( subset(tmp, N>1,PANGEA_ID), ds, by='PANGEA_ID')
	#	21531 consensus files corresponding to PANGEA_ID with >1 consensus
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
	#	59067 files corresponding to PANGEA_ID with one selected consensus
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
	#	67122 files corresponding to PANGEA_ID with one selected consensus
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
	#	7462 samples with consensus found
	

	#	add length of consensus
	tmp <- ds[, list(CONS_N_ACTG=nchar(gsub('\\?|\\-','',paste0(as.character(read.dna(FORGLOBALALN_FASTA, format='fasta')), collapse='')))), by=c('PANGEA_ID')]
	ds <- merge(ds, tmp, by='PANGEA_ID')
	ds <- subset(ds, CONS_N_ACTG>0)
	#	6855 samples with non-empty consensus found (although they can be quite small)
	ds[, CONS_N_ACTG_C:= cut(CONS_N_ACTG, breaks=seq(0,10e3,5e2))]
	ds[, table(CONS_N_ACTG_C)]
	#	CONS_N_ACTG_C
    #   	 (0,500]     (500,1e+03] (1e+03,1.5e+03] (1.5e+03,2e+03] (2e+03,2.5e+03]
    #   	     551             255             181             464             291
	#(2.5e+03,3e+03] (3e+03,3.5e+03] (3.5e+03,4e+03] (4e+03,4.5e+03] (4.5e+03,5e+03]
	#            224             175             169             125             254
	#(5e+03,5.5e+03] (5.5e+03,6e+03] (6e+03,6.5e+03] (6.5e+03,7e+03] (7e+03,7.5e+03]
	#            171             162             173             190             264
	#(7.5e+03,8e+03] (8e+03,8.5e+03] (8.5e+03,9e+03] (9e+03,9.5e+03] (9.5e+03,1e+04]
	#            224             246            2237             499               0
	
	saveRDS(ds, file=file.path(indir, '200422_PANGEA2_RCCS_selected_samples.rds'))
}

rk.seq.get.best.mappeable.samples.MRC <- function()
{
	require(data.table)
	require(ape)
	indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live'
	indir.mrc <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC'
	
	
	ds <- readRDS( file.path(indir.mrc, '200422_PANGEA2_MRCUVRI_mapped_samples.rds') )
	ds[, length(unique(PANGEA_ID))]
	#	2504 mappeable samples 
	
	#
	#	find samples who also have a consensus sequence
	#
	tmp <- unique(subset(ds, CONSENSUS!='', c(PANGEA_ID, BATCH_ID, PREFIX)))
	tmp[, length(unique(PANGEA_ID))]
	#	2397 mappeable samples with consensus sequence, keep only those
	ds <- merge(tmp, ds, by=c('PANGEA_ID','PREFIX','BATCH_ID'))
	
	
	#
	#	select one consensus sequence for each sample
	#
	tmp <- ds[CONSENSUS!='', list(N=length(CONSENSUS)), by='PANGEA_ID']	
	ds2 <- merge(subset(tmp, N==1,PANGEA_ID), ds, by='PANGEA_ID')
	#	21546 files corresponding to PANGEA_ID with one available consensus	
	ds3 <- merge( subset(tmp, N>1,PANGEA_ID), ds[CONSENSUS!='',], by='PANGEA_ID')
	ds4 <- merge( subset(tmp, N>1,PANGEA_ID), ds, by='PANGEA_ID')
	#	54 consensus files corresponding to PANGEA_ID with >1 consensus	
	#	files are re-sequenced, keep the one with longest consensus
	tmp <- ds3[, list(N_ACTG=nchar(gsub('\\?|\\-','',paste0(as.character(read.dna(F, format='fasta')), collapse='')))), by=c('PANGEA_ID','PREFIX','BATCH_ID')]
	tmp <- tmp[, {
				z<- which.max(N_ACTG)[1]
				list(BATCH_ID= BATCH_ID[z], PREFIX=PREFIX[z])
			}, by='PANGEA_ID']
	tmp2 <- merge( tmp, ds4, by=c('PANGEA_ID','PREFIX','BATCH_ID') )
	ds <- rbind(ds2, tmp2)	
	#	21573 files corresponding to PANGEA_ID with one selected consensus
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
	#	2397 samples with consensus found
	
	
	#	add length of consensus
	tmp <- ds[, list(CONS_N_ACTG=nchar(gsub('\\?|\\-','',paste0(as.character(read.dna(FORGLOBALALN_FASTA, format='fasta')), collapse='')))), by=c('PANGEA_ID')]
	ds <- merge(ds, tmp, by='PANGEA_ID')
	ds <- subset(ds, CONS_N_ACTG>0)
	#	2042 samples with non-empty consensus found (although they can be quite small)
	ds[, CONS_N_ACTG_C:= cut(CONS_N_ACTG, breaks=seq(0,10e3,5e2))]
	ds[, table(CONS_N_ACTG_C)]
	#	CONS_N_ACTG_C
	#        (0,500]     (500,1e+03] (1e+03,1.5e+03] (1.5e+03,2e+03] (2e+03,2.5e+03]
    #         81              25              22              90              33
	#(2.5e+03,3e+03] (3e+03,3.5e+03] (3.5e+03,4e+03] (4e+03,4.5e+03] (4.5e+03,5e+03]
    #         66              51              53              31              68
	#(5e+03,5.5e+03] (5.5e+03,6e+03] (6e+03,6.5e+03] (6.5e+03,7e+03] (7e+03,7.5e+03]
    #         67              58              69             227             197
	#(7.5e+03,8e+03] (8e+03,8.5e+03] (8.5e+03,9e+03] (9e+03,9.5e+03] (9.5e+03,1e+04]
    #        154             208             426             116               0
	
	saveRDS(ds, file=file.path(indir.mrc, '200422_PANGEA2_MRCUVRI_selected_samples.rds'))
}

rk.seq.get.mappeable.samples.MRC <- function()
{
	require(data.table)	
	indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC'
	outdir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC'
	
	pangea.files <- file.path(indir, '200319_pangea_db_sharing_extract_mrc.csv')
	pangea.files2 <- file.path(indir, '200319_pangea_mappings_mrc.csv')
	pangea.filesm <- file.path(indir, '200323_MRC_data_deprecated.csv')
	
	
	#
	#	read MRC info from Oxford
	#
	pf <- as.data.table(read.csv(pangea.files, stringsAsFactors=FALSE))
	pf <- subset(pf, select=c(pt_id, geo_country, cohort_id, main_cohort_id, visit_dt, pangea_id, readnum, readnum_hiv, mapped_num,  duprate, insertfrac_350, lc_cutoff, n_cutoff, length_strict, subtype_bestref))
	pf <- subset(pf, select=c(pt_id, geo_country, cohort_id, main_cohort_id, visit_dt, pangea_id, subtype_bestref))
	setnames(pf, colnames(pf), toupper(colnames(pf)))	
	pf <- unique(pf)
	# 2931
		
	#	read sequence IDs and path to data from Oxford
	tmp <- as.data.table(read.csv(pangea.files2, stringsAsFactors=FALSE))	
	setnames(tmp, colnames(tmp), toupper(colnames(tmp)))	
	tmp[, DEPRECATED:= 'N']
	tmp2 <- as.data.table(read.csv(pangea.filesm, stringsAsFactors=FALSE))
	setnames(tmp2, colnames(tmp2), toupper(colnames(tmp2)))
	tmp2[, SEQUENCE_ID:= gsub('PG[0-9]{2}-UG[0-9]{6}-([0-9]{5}_[0-9]_[0-9]+)','\\1',PREFIX)]
	tmp <- rbind(tmp, subset(tmp2, select=c(PANGEA_ID, SEQUENCE_ID, DEPRECATED)), fill=TRUE)
	set(tmp, tmp[, which(SEQUENCE_ID=='NULL')], 'SEQUENCE_ID', NA_character_)
	set(tmp, tmp[, which(PATH_TO_DATA=='NULL')], 'PATH_TO_DATA', NA_character_)	
	tmp[, BATCH_ID:= gsub('.*(PANGEA1_miseq|PANGEA1_hiseq).*','\\1',basename(PATH_TO_DATA))]
	tmp <- unique(tmp)
	# 2944
	
	# merge data from Oxford
	pf <- merge(tmp, pf, by='PANGEA_ID', all=TRUE)		
	# 2949	
	
	#
	#	read MRC data at Imperial
	#
	ds <- data.table(F=list.files(indir, recursive=TRUE, full.names=TRUE))
	#	23027
	
	
	#save(pf, ds, file= file.path(outdir, 'PANGEA2_RCCS_interim.rda'))
	#load(file.path(outdir, 'PANGEA2_RCCS_interim.rda'))
	
	#
	#	process file ends
	tmp <- '.*(_ref.fasta|.bam|.bam.bai|_ForGlobalAln.fasta)$'
	ds[, FEND:= gsub(tmp,'\\1',F)]
	#	files to be ignored..
	subset(ds, !grepl(tmp, FEND))
	#	continue
	ds <- subset(ds, grepl(tmp, FEND))
	ds[, PREDEDUP:= as.character(factor(grepl('PreDedup',basename(F)), levels=c(TRUE,FALSE), labels=c('_PreDeDup','')))]
	ds[, REMAP:= as.character(factor(grepl('remap',basename(F)), levels=c(TRUE,FALSE), labels=c('_remap','')))]
	tmp <- ds[, which(grepl('_consensus_',basename(F)))]
	ds[, CONSENSUS:= '']
	set(ds, tmp, 'CONSENSUS', ds[tmp, gsub('.*(_consensus_MinCov_[0-9]+_[0-9]+_ForGlobalAln.fasta).*','\\1',basename(F))])
	ds[, FEND2:= paste0(REMAP,PREDEDUP,FEND)]
	tmp <- ds[, which(nchar(CONSENSUS)>0)]
	set(ds, tmp, 'FEND2', ds[tmp,CONSENSUS])
	#	23025
	
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
	#tmp <- 'Neg|neg|^SC|^ctrl|^100$|^500$|^5000$|^50000$|^500000$|^5000000$|NITP'
	#ds2 <- subset(ds2, !grepl(tmp, PREFIX))
	#	check for old Sanger IDs of format 13549_1_84 and set SEQUENCE_ID accordingly
	#tmp <- ds2[, which(grepl('.*([0-9]{5}_[0-9]_[0-9]+).*', PREFIX))]
	#set(ds2, tmp, 'SEQUENCE_ID', ds2[tmp, gsub('.*([0-9]{5}_[0-9]_[0-9]+).*','\\1',PREFIX)])
	#tmp <- ds2[, which(is.na(SEQUENCE_ID))]
	#set(ds2, tmp, 'SEQUENCE_ID', ds2[tmp, PREFIX])
	#	build PANGEA ID		 
	ds2[, PID:= gsub('^(PG[0-9]{2}-[A-Z]{2}[0-9]{6}).*', '\\1', PREFIX)]
	ds <- merge(ds, ds2, by='ROW_ID') 
	ds[, PANGEA_ID:= PID]
	# 23025
	
	#
	#	prepare linking to PANGEA_IDs from Oxford
	pf2 <- unique(subset(pf, select=c(PANGEA_ID, SEQUENCE_ID, PT_ID, VISIT_DT, DEPRECATED)))
	ds2 <- merge(ds, pf2, by='PANGEA_ID', all=TRUE)
	#	23447
	#	there are entries for which there are no samples at Imperial and Oxford, remove them
	nrow(subset(ds2, is.na(F) & is.na(SEQUENCE_ID)))
	#	419
	ds2 <- subset(ds2, !(is.na(F) & is.na(SEQUENCE_ID)) )	
	#	23028
	

	#
	#	PANGEA IDs at Oxford, which have no bam or reference files on rescomp
	#	likely shiver failed
	ds3 <- subset(ds2, is.na(F), c('PANGEA_ID'))
	write.csv(ds3, file=file.path(outdir, '200326_MRC_likelyshiveraborted.csv'))
	#	3
	ds2 <- subset(ds2, !is.na(F))
	#	23025


	#
	#	remove deprecated files
	ds2 <- subset(ds2, DEPRECATED=='N' | is.na(DEPRECATED))
	#	22863
	
	#	files I have but not in Oxford
	unique(subset(ds2, is.na(PT_ID)), by='PREFIX')
	#	124
	write.csv(unique(subset(ds2, is.na(PT_ID)), by='PREFIX'), file=file.path(outdir,'200323_potentiallyMRC_notmappeable.csv'), row.names=FALSE)	

		
	#	proceed with rest given lockdown at UVRI
	ds2 <- subset(ds2, !is.na(PT_ID))
	saveRDS(ds2, file=file.path(outdir, '200422_PANGEA2_MRCUVRI_mapped_samples.rds'))

	#	PANGEA IDs at Oxford, which may still be in the Rakai folder, copied over
	if(0)
	{
		indir.rk.miseq <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/PANGEA1_miseq'
		indir.rk.hiseq <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/PANGEA1_hiseq'
		indir.mrc.data <- file.path(indir,'PANGEA1')
		ds4 <- rbind( data.table(F=list.files(indir.rk.miseq, recursive=TRUE, full.names=TRUE)), data.table(F=list.files(indir.rk.hiseq, recursive=TRUE, full.names=TRUE)) )
		tmp <- ds3[, list(IN_RCCS_FOLDER= grepl(PANGEA_ID, basename(ds4$F) )), by='PANGEA_ID']
		tmp <- unique(subset(tmp, IN_RCCS_FOLDER, PANGEA_ID))
		#	success, found data for 107 samples
		z<- tmp[, {
					from <- ds4$F[ which(grepl(PANGEA_ID, basename(ds4$F) )) ]
					to <- file.path(indir.mrc.data, basename(from))
					list(OK= file.rename(from, to))				
				}, by='PANGEA_ID']
	}
}

rk.seq.get.mappeable.samples.Rakai <- function()
{
	require(data.table)
	indir.rk <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS'
	outdir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live'
		
	pangea.files <- file.path(indir.rk, '200316_pangea_db_sharing_extract_rakai.csv')
	pangea.files2 <- file.path(indir.rk, '200316_pangea_mappings_rakai.csv')	
	pangea.filesm <- file.path(indir.rk, '200316_unlinkable_files_mappings.csv')
	
	#
	#	read RCCS info from Oxford
	#
	pf <- as.data.table(read.csv(pangea.files, stringsAsFactors=FALSE))
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
	tmp <- '.*([A-Z][0-9]{6}R1[7-9]).*'
	tmp2 <- ds2[, which(grepl(tmp, PREFIX))]
	set(ds2, tmp2, 'PID', ds2[tmp2,gsub(tmp,'\\1',PREFIX)])
	tmp <- '.*([A-Z][0-9]{6}NEU).*'
	tmp2 <- ds2[, which(grepl(tmp, PREFIX))]
	set(ds2, tmp2, 'PID', ds2[tmp2,gsub(tmp,'\\1',PREFIX)])
	ds <- merge(ds, ds2, by='ROW_ID') # this removes controls
	# 81045
		
	#
	#	prepare linking to PANGEA_IDs from Tanya
	pf2 <- unique(subset(pf, select=c(PANGEA_ID, BATCH_ID, SEQUENCE_ID, PT_ID, VISIT_DT, DEPRECATED, LINKABLE)))
	#	8358	
	ds2 <- unique(subset(ds, select=c(PID, BATCH_ID, SEQUENCE_ID, PREFIX)))
	#	9091

	ds3 <- merge(ds2, pf2, by=c('SEQUENCE_ID','BATCH_ID'), all=TRUE)
	dc <- subset(ds3, !is.na(PT_ID) & !is.na(PID))
	#	matched 8127 entries
	
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
	write.csv(file=file.path(outdir, '200318_RCCS2_PANGEA_resolvedPANGEAIDs_NEU.csv'), subset(dm, !is.na(PANGEA_ID)), row.names=FALSE)
	write.csv(file=file.path(outdir, '200318_RCCS2_PANGEA_missing_data_NEU.csv'), subset(dm, is.na(PANGEA_ID)), row.names=FALSE)
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
	#	79336
	saveRDS(ds, file=file.path(outdir, '200422_PANGEA2_RCCS_mapped_samples.rds'))	
}

rk.chase.missing.demographics <- function()
{
	require(data.table)
	infile.missing <- '~/Box/OR_Work/PANGEA/2020_PANGEA_UGanalyses/data_queries/200618_RCCS_missing_demographics.csv'
	dm <- as.data.table(read.csv(infile.missing, stringsAsFactors=FALSE))
	dm[, PID:= gsub('RCCS_','',pangea_id)]
	
	#	the UG9000* was a test plate, remove
	dm <- subset(dm, !grepl('.*UG900.*',PID))
	
	if(0)
	{
		infile.assemblystatus <- '~/Box/OR_Work/PANGEA/PANGEA_data/Rakai_phyloscanner_170704_assemblystatus.rda'
		load(infile.assemblystatus)	
		dm <- merge(subset(dm, select=PID), dc, by='PID')		
	}
	
	infile.metadata <- '~/Box/OR_Work/PANGEA/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda'
	load(infile.metadata)
	dr <- subset(as.data.table(rccsData), select=c(RCCS_studyid,Pangea.id, batch, date, SEX))
	setnames(dr, 'Pangea.id', 'PID')
	merge(dm, dr, by='PID')
}

rk.stage1.close.pairs <- function()
{
	infile <- '~/Box/OR_Work/PANGEA/2020_PANGEA_UGanalyses/BLAST_close_pairs/200518_closepairs_version2.rda'
	load(infile)
}