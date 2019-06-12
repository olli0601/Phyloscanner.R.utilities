couples.find	<- function()
{
	require(data.table)
	
	if(0)
	{
		linked.group	<- 'TYPE_PAIR_TODI2'
		linked.type.yes	<- 'linked'			
		dir.group		<- 'TYPE_DIR_TODI2'		
	}
	if(1)
	{
		linked.group	<- 'TYPE_CHAIN_TODI'
		linked.type.yes	<- 'chain'
		linked.type.no	<- 'distant'
		dir.group		<- 'TYPE_ADJ_DIR_TODI2'		
	}
	
	indir	<- '~/Box Sync/OR_Work/MSCs/2019_Melodie'
	outdir	<- '~/Box Sync/OR_Work/MSCs/2019_Melodie'
	
	# load phyloscanner results for couples
	infile.trmpairs.todi	<- file.path(indir, "data", "todi_pairs_171122_cl25_d50_prior23_min30_withmetadata.rda")							   	
	load(infile.trmpairs.todi)	
	rca	<- subset(rtp, COUPLE!='no couple')
	
	# add anonymised labels (need this to interface with published data)
	infile.anonymised	<- file.path(indir, "data", "todi_pairs_171122_cl25_d50_prior23_min30_anonymised_RIDs.csv")
	dan	<- as.data.table(read.csv(infile.anonymised))	
	setnames(dan, c('ID','AID'), c('MALE_RID','MALE_ARID'))
	rca		<- merge(rca, subset(dan, select=c(MALE_RID, MALE_ARID)), by='MALE_RID')
	setnames(dan, c('MALE_RID','MALE_ARID'), c('FEMALE_RID','FEMALE_ARID'))
	rca		<- merge(rca, subset(dan, select=c(FEMALE_RID, FEMALE_ARID)), by='FEMALE_RID')
	# add sampling dates
	tmp		<- rs[, list(SEQ_DATE=min(SEQ_DATE)), by='RID'] 
	setnames(tmp, c('RID','SEQ_DATE'), c('MALE_RID','MALE_SEQ_DATE'))
	rca		<- merge(rca, tmp, by='MALE_RID')
	setnames(tmp, c('MALE_RID','MALE_SEQ_DATE'), c('FEMALE_RID','FEMALE_SEQ_DATE'))
	rca		<- merge(rca, tmp, by='FEMALE_RID')
		
	
	#
	#	define meta-variables on epidemiologic evidence for direction 
	#	of transmission
	rca[, EXT_TYPE:=NA_character_]
	set(rca, rca[, which(FEMALE_LASTNEGDATE>=MALE_FIRSTPOSDATE)], 'EXT_TYPE', 'serodisc-mf')		
	set(rca, rca[, which(MALE_LASTNEGDATE>=FEMALE_FIRSTPOSDATE)], 'EXT_TYPE', 'serodisc-fm')	
	#	add extra couples in who one has CD4<400 and the other has CD4>800
	#	for male potential recipient - evaluate CD4 around time male first positive
	tmp	<- rca[, which(is.na(EXT_TYPE) &
							MALE_RECENTCD4>800 & abs(MALE_RECENTCD4DATE-MALE_FIRSTPOSDATE)<2 &
							FEMALE_RECENTCD4<400 & abs(FEMALE_RECENTCD4DATE-MALE_FIRSTPOSDATE)<2)]
	set(rca, tmp, 'EXT_TYPE', 'CD4disc-fm')
	#	for female potential recipient - evaluate CD4 around time female first positive
	tmp	<- rca[, which(is.na(EXT_TYPE) &
							FEMALE_RECENTCD4>800 & abs(FEMALE_RECENTCD4DATE-FEMALE_FIRSTPOSDATE)<2 &
							MALE_RECENTCD4<400 & abs(MALE_RECENTCD4DATE-FEMALE_FIRSTPOSDATE)<2)]
	set(rca, tmp, 'EXT_TYPE', 'CD4disc-mf')
	#	separate into EXT_TYPE and EXT_DIR
	set(rca, NULL, 'EXT_DIR', rca[, gsub('^([a-zA-Z0-9]+)-([a-zA-Z]+)$','\\2',EXT_TYPE)])
	set(rca, NULL, 'EXT_TYPE', rca[, gsub('^([a-zA-Z0-9]+)-([a-zA-Z]+)$','\\1',EXT_TYPE)])
	
	
	#
	#	define phylogenetically inferred direction of transmission
	rca[, PHSC_DIR:= NA_character_]
	set(rca, rca[, which(grepl('direction resolved',SELECT) & POSTERIOR_SCORE_MF>POSTERIOR_SCORE_FM)], 'PHSC_DIR', 'mf')
	set(rca, rca[, which(grepl('direction resolved',SELECT) & POSTERIOR_SCORE_MF<POSTERIOR_SCORE_FM)], 'PHSC_DIR', 'fm')
	
	#	save
	save(rca, file=file.path(outdir,'analyses','190611_phsccouples.rda'))
}

couples.csv.for.Tanya<- function()
{
	indir	<- '~/Box Sync/OR_Work/MSCs/2019_Melodie'
	infile.couples <- file.path(indir,'analyses','190611_phsccouples.rda')	
	load(infile.couples)
	
	write.csv( data.table(RID=sort( unique(c( rca[, MALE_RID], rca[, FEMALE_RID] )) )), row.names=FALSE, file=file.path(indir,'analyses','190611_phsccouples_NeedInfTimes.csv')) 
}

couples.distances.from.trees <- function()
{
	require(data.table)
	require(tidyverse)
	require(ape)
	require(phytools)
	
	#	locally
	if(0)
	{
		indir	<- '~/Box Sync/OR_Work/MSCs/2019_Melodie'
		infile.couples <- file.path(indir,'analyses','190611_phsccouples.rda')	
		tree.dir <- "~/sandbox/DeepSeqProjects/RakaiPopSample_deepseqtrees"
	}
	
	#	on HPC
	if(1)
	{
		indir	<- '/rds/general/user/or105/home/WORK/Gates_2014/networks'
		infile.couples <- file.path(indir,'analyses','190611_phsccouples.rda')	
		tree.dir <- "/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS_trees"		
	}
	
	
	#	download the phyloscanner tree of the Rakai population-based analysis
	if(0)
	{
		tmp <- "https://datadryad.org/bitstream/handle/10255/dryad.208473/Dataset_S1.tar?sequence=1"
		#	specify directory to untar public data		
		#	download and untar
		download.file(tmp, destfile="Dataset_S1.tar", method="curl")
		untar("Dataset_S1.tar", exdir=tree.dir, extras='-xvf')	
	}
		
	#	list zipped tree files. One zip file contains the viral trees of individuals in one putative transmission network. 
	df <- tibble(F=list.files(tree.dir))
	df <- df %>% 
			mutate(TYPE:= gsub('ptyr([0-9]+)_(.*)','\\2', F),
					PTY_RUN:= as.integer(gsub('ptyr([0-9]+)_(.*)','\\1', F))) %>%
			mutate(TYPE:= gsub('^([^\\.]+)\\.[a-z]+$','\\1',TYPE)) %>%
			spread(TYPE, F) %>%
			set_names(~ str_to_upper(.))
	
	#	load couples and select those between whom we want to compute tip-to-tip distances
	load(infile.couples)
	set(rca, NULL, 'MALE_ARID', as.character(rca$MALE_ARID))
	set(rca, NULL, 'FEMALE_ARID', as.character(rca$FEMALE_ARID))
	rcl	<- subset(rca, grepl('most likely',SELECT))
	
	#	select phyloscanner runs with linked couples
	tmp	<- unique(subset(rcl, select=PTY_RUN))
	df	<- merge(tmp, as.data.table(df), by='PTY_RUN')
	
	#	for each run, unzip and calculate tip-to-tip distances among linked couples	
	for(run in df$PTY_RUN[-1])
	{
		#run	<- 5	#for devel
		#	unzip
		tmpfile	<- subset(df,PTY_RUN==run)[,file.path(tree.dir, TREES_NEWICK)]
		tmpfile2<- file.path(indir,'analyses',basename(gsub('\\.zip','',tmpfile)))
		unzip(tmpfile, exdir=tmpfile2)
		#	read all trees across genome
		tmp	<- list.files(tmpfile2, pattern='tree$', full.names=TRUE)
		trees <- lapply(tmp, ape:::read.tree)
		names(trees)	<- basename(tmp)
		#	get taxa names between whom we want to compute tip-to-tip calculate distances		
		tmp	<- lapply(seq_along(trees), function(j) {
					data.table(TREE= names(trees)[j], TAXA= trees[[j]][['tip.label']])
				})
		dt	<- do.call('rbind',tmp)
		dt	<- subset(dt, !grepl('^REF', TAXA))
		dt[, ARID:= gsub('^([A-Za-z0-9]+)_.*','\\1',TAXA)]
		rclr<- subset(rcl, PTY_RUN==run)
		dt	<- subset(dt, ARID%in%unique(c(rclr$FEMALE_ARID, rclr$MALE_ARID)))  
		
		ans	<- vector('list', 0)
		for(j in seq_len(nrow(rclr)))
			for(t in unique(dt$TREE))
			{
				#t	<- "ptyr5_InWindow_1000_to_1249.tree" #for devel
				# get data.table of male tips vs female tips
				zm	<- subset(dt, TREE==t & ARID==rclr$MALE_ARID[j])
				zf	<- subset(dt, TREE==t & ARID==rclr$FEMALE_ARID[j])
				setnames(zm, c('TAXA','ARID'), c('MALE_TAXA','MALE_ARID'))
				setnames(zf, c('TAXA','ARID'), c('FEMALE_TAXA','FEMALE_ARID'))
				z	<- merge(zm, zf, by='TREE', allow.cartesian=TRUE)
				# drop tips to calculate distances as fast as possible, calculate distances
				tmp	<- ape:::keep.tip(trees[[t]], subset(dt, TREE==t)[, TAXA])
				tmp	<- adephylo:::distTips( tmp )
				# merge with z to extract those distances that we want
				tmp	<- as.matrix(tmp)
				tmp	<- as.data.table(melt(tmp, as.is=TRUE, varnames=c('MALE_TAXA','FEMALE_TAXA'), value.name = "PD"))
				z	<- merge(z, tmp, by=c('MALE_TAXA','FEMALE_TAXA'))
				# return mean tip-to-tip distance
				ans[[length(ans)+1L]]	<- data.table(TREE=t, MALE_ARID=z$MALE_ARID[1], FEMALE_ARID=z$FEMALE_ARID[1], MPD=z[, mean(PD)])
			}
		ans	<- do.call('rbind', ans)
		ans	<- subset(ans, !is.na(MALE_ARID) & !is.na(FEMALE_ARID))
		write.csv(ans, row.names=FALSE, file=file.path(indir,'analyses',paste0('190611_phsccouples_mpd_',run,'.csv')))
		unlink(tmpfile2, recursive=TRUE)
	}
}