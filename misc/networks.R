#' Make list of RCCS couples with phyloscanner output
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

#' Ask Tanya to provide infection time estimates for RCCS couples
couples.csv.for.Tanya<- function()
{
	indir	<- '~/Box Sync/OR_Work/MSCs/2019_Melodie'
	infile.couples <- file.path(indir,'analyses','190611_phsccouples.rda')	
	load(infile.couples)
	
	write.csv( data.table(RID=sort( unique(c( rca[, MALE_RID], rca[, FEMALE_RID] )) )), row.names=FALSE, file=file.path(indir,'analyses','190611_phsccouples_NeedInfTimes.csv')) 
}

#' Find relationship statistics for RCCS couples
couples.relationship.stats	<- function()
{
	require(data.table)
	require(tidyverse)
	require(ape)
	require(phytools)
	
	#	locally
	if(1)
	{
		indir	<- '~/Box Sync/OR_Work/MSCs/2019_Melodie'
		infile.couples <- file.path(indir,'analyses','190611_phsccouples.rda')	
		rel.dir <- "~/sandbox/DeepSeqProjects/RakaiPopSample_phyloscanner_analysis"
	}
	
	#	list relationship files
	df <- tibble(F=list.files(rel.dir)) %>% 
			filter(grepl('pairwise_relationships.rda',F)) %>%
			mutate(PTY_RUN:= as.integer(gsub('ptyr([0-9]+)_(.*)','\\1', F))) 
			
	#	load couples and select those between whom we want to compute tip-to-tip distances
	load(infile.couples)
	set(rca, NULL, 'MALE_ARID', as.character(rca$MALE_ARID))
	set(rca, NULL, 'FEMALE_ARID', as.character(rca$FEMALE_ARID))
	rcl	<- subset(rca, grepl('most likely',SELECT))
	
	#	select phyloscanner runs with linked couples
	tmp	<- unique(subset(rcl, select=PTY_RUN))
	df	<- merge(tmp, as.data.table(df), by='PTY_RUN')
	
	#	for each run, unzip and calculate tip-to-tip distances among linked couples
	ans	<- vector('list', 0)
	for(run in df$PTY_RUN)
	{
		#run	<- 5	#for devel
		
		#	load relationship file and extract basic statistics for all windows
		tmpfile	<-  subset(df,PTY_RUN==run)[,file.path(rel.dir, F)]
		load(tmpfile)
		dwin 	<- subset(dwin, select=c(SUFFIX, W_FROM, W_TO, ID1, ID2, ADJACENT, CONTIGUOUS, PATHS_12, PATHS_21))
		dwin[, ID1_GENDER:= gsub('.*([MF])$','\\1',ID1)]
		dwin[, ID2_GENDER:= gsub('.*([MF])$','\\1',ID2)]
		dwin 	<- subset(dwin, ID1_GENDER=='F' & ID2_GENDER=='M'  |   ID1_GENDER=='M' & ID2_GENDER=='F')
		tmp		<- subset(dwin, ID1_GENDER=='F' & ID2_GENDER=='M')
		setnames(tmp, c('ID1','ID2','PATHS_12','PATHS_21','ID1_GENDER','ID2_GENDER'), c('ID2','ID1','PATHS_21','PATHS_12','ID2_GENDER','ID1_GENDER'))
		dwin	<- rbind( subset(dwin, ID1_GENDER=='M' & ID2_GENDER=='F'), tmp)
		
		#	select couples
		tmp		<- data.table(ID1= rcl$MALE_ARID, ID2= rcl$FEMALE_ARID)
		dwin	<- merge(tmp, dwin, by=c('ID1','ID2'))
		
		#	return
		dwin[, PTY_RUN:= run]
		ans[[length(ans)+1L]]	<- copy(dwin)
	}
	ans	<- do.call('rbind', ans)
	setkey(ans, PTY_RUN, ID1, ID2)
	write.csv(ans, row.names=FALSE, file=file.path(indir,'analyses','190611_phsccouples_basic_topology_stats.csv'))
}


#' Find relationship statistics for RCCS couples based on the latest version of phyloscanner, version 1.8.0
couples.stats.phyloscanner.1.8.0	<- function()
{
	require(data.table)
	require(tidyverse)
	require(ape)
	require(phytools)
	
	#	locally
	if(1)
	{
		indir	<- '~/Box Sync/OR_Work/MSCs/2019_Melodie'
		infile.couples <- file.path(indir,'analyses','190611_phsccouples.rda')	
		rel.dir <- "~/sandbox/DeepSeqProjects/RakaiPopSample_phsc_out190626"
	}
	
	#	list relationship files
	df <- tibble(F=list.files(rel.dir)) %>% 
			filter(grepl('workspace.rda',F)) %>%
			mutate(PTY_RUN:= as.integer(gsub('ptyr([0-9]+)_(.*)','\\1', F))) 
	
	#	load couples and select those between whom we want to compute tip-to-tip distances and relationships
	load(infile.couples)
	set(rca, NULL, 'MALE_ARID', as.character(rca$MALE_ARID))
	set(rca, NULL, 'FEMALE_ARID', as.character(rca$FEMALE_ARID))
	rcl	<- subset(rca, grepl('most likely',SELECT))
	
	#	select phyloscanner runs with linked couples
	tmp	<- unique(subset(rcl, select=PTY_RUN))
	df	<- merge(tmp, as.data.table(df), by='PTY_RUN')
	
	#	for each run, unzip and calculate tip-to-tip distances among linked couples
	ans	<- vector('list', 0)
	for(run in df$PTY_RUN)
	{
		print(run)
		#run	<- 5	#for devel
		
		#	load phyloscanner output file and extract mean tip to tip distances
		tmpfile	<-  subset(df,PTY_RUN==run)[,file.path(rel.dir, F)]
		load(tmpfile)
		
		#	
		rclr	<- subset(rcl, PTY_RUN==run)
		for(j in seq_len(nrow(rclr)))
			for(k in seq_along(phyloscanner.trees))
			{
				dc	<- phyloscanner.trees[[k]][['classification.results']][['classification']]
				#	extract topology statistics
				dc	<- dc %>% 
						filter( (host.1==rclr[j, FEMALE_ARID] & host.2==rclr[j, MALE_ARID]) | 
								(host.2==rclr[j, FEMALE_ARID] & host.1==rclr[j, MALE_ARID])	)				
				#	extract tip to tip distances
				tmp	<- c( phyloscanner.trees[[k]][['tips.for.hosts']][[rclr[j, FEMALE_ARID]]], phyloscanner.trees[[k]][['tips.for.hosts']][[rclr[j, MALE_ARID]]])
				if(length(tmp)>0 && nrow(dc)>0)
				{
					tmp	<- ape:::keep.tip(phyloscanner.trees[[k]][['tree']], tmp)
					tmp	<- adephylo:::distTips( tmp )
					tmp	<- as.matrix(tmp)
					tmp	<- as_tibble(melt(tmp, as.is=TRUE, varnames=c('MALE_TAXA','FEMALE_TAXA'), value.name = "PD"))
					#	calculate mean tip to tip distance between taxa from male and female
					tmp	<- tmp	%>%
							mutate(	MALE_ARID= gsub('^([^_]+).*','\\1', MALE_TAXA),
									FEMALE_ARID= gsub('^([^_]+).*','\\1', FEMALE_TAXA) ) %>%
							filter( MALE_ARID!=FEMALE_ARID) %>%
							summarise(MPD=mean(PD)) %>%
							as.double()
					dc	<- dc %>%
							mutate(	MPD:= tmp, 
									MPD_STD:= tmp/phyloscanner.trees[[k]][['normalisation.constant']],
									W_FROM:= phyloscanner.trees[[k]][['window.coords']][['start']],
									W_TO:= phyloscanner.trees[[k]][['window.coords']][['end']],
									PTY_RUN:= run)				
					# return 
					ans[[length(ans)+1L]]	<- dc
				}				
			}
	}
		
	
	ans	<- do.call('rbind', ans)
	
	save(ans, file=file.path(indir,'analyses','190627_phsccouples_topology_and_distance.rda'))
	
}


#' Collect initial data set for Melodie
networks.make.data.190726<- function()
{
	require(tidyverse)
	
	infile <- '~/sandbox/DeepSeqProjects/RakaiPopSample_phyloscanner_analysis_190706/Rakai_phscnetworks_190706.rda'
	load(infile)
	
	#	select transmission networks
	idclus <- c(5, 24, 26, 31, 78, 87, 106, 458, 459, 463)
	
	#	get pairs of individuals in transmission networks
	tmp <- dnet %>% 
			filter(IDCLU %in% idclus) %>% 
			select(-c(CATEGORISATION, TYPE, SCORE, H1_SEX, H1_SAMPLING_YEAR, H2_SEX, H2_SAMPLING_YEAR)) %>% 
			distinct() %>%
			arrange(IDCLU, H1, H2)
	#	make all pairwise combinations of individuals
	dnet2 <- tmp %>% select( -c(CLU_SIZE,PTY_RUN) ) %>%
			gather('variable','host',H1,H2) %>%
			select( -variable ) %>%
			arrange(IDCLU,host) %>%
			distinct() %>%
			group_by(IDCLU) 
	dnet2 <- as.data.table(dnet2)[, {
				tmp <- as.data.table(t(combn(host, 2)))
				setnames(tmp, c('V1','V2'), c('H1','H2'))
				tmp
			}, by='IDCLU']
	dnet2 <- as_tibble(dnet2) %>% 
			filter(H1<H2)	
	dnet2 <- dnet2 %>% 
			full_join( tmp %>% select(PTY_RUN,IDCLU) %>% distinct(), by='IDCLU')
	
	
	
			
	#	go back to PTY_RUN data to extract window information
	idruns <- sort( unique( dnet2$PTY_RUN ) )
	indir <- '~/sandbox/DeepSeqProjects/RakaiPopSample_phyloscanner_analysis_190706'
	infiles <- tibble(F=list.files(indir, pattern='workspace.rda$',full.names=TRUE)) %>%
				mutate(PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*$','\\1',basename(F)))) %>%
				filter(PTY_RUN %in% idruns)	
	dwin <- lapply( seq_len(nrow(infiles)), function(i){
				cat('\nprocess', infiles$F[i])
				load( infiles$F[i] )								
				dprs <- dnet2 %>% filter(PTY_RUN==infiles$PTY_RUN[i])
				ans <- vector('list',0)
				for(j in seq_len(nrow(dprs)))
					for(k in seq_along(phyloscanner.trees))
					{
						dc	<- phyloscanner.trees[[k]][['classification.results']][['classification']]
						#	extract topology statistics
						dc	<- dc %>% 
								filter( (host.1==dprs$H1[j] & host.2==dprs$H2[j]) | 
										(host.2==dprs$H1[j] & host.1==dprs$H2[j])	) %>%
								rename( H1:= host.1, H2:=host.2) 
						tmp <- dc %>% 
								filter(H1>H2) %>%
								rename(H1:=H2, H2:=H1, paths12:=paths21, paths21:=paths12, nodes1:=nodes2, nodes2:=nodes1) %>%
								mutate(ancestry:= gsub('tmp$','anc',gsub('Tmp$','Anc',gsub('anc$','desc',gsub('Anc$','Desc',gsub('desc$','tmp',gsub('Desc$','Tmp',ancestry)))))))
						dc <- dc %>% 
								filter(H1<H2) %>%
								bind_rows(tmp) %>%
								inner_join(dprs[j,], by=c('H1','H2'))
						
						#	extract tip to tip distances
						tmp	<- c( phyloscanner.trees[[k]][['tips.for.hosts']][[dprs$H1[j]]], phyloscanner.trees[[k]][['tips.for.hosts']][[dprs$H2[j]]])
						if(length(tmp)>0 && nrow(dc)>0)
						{
							tmp	<- ape:::keep.tip(phyloscanner.trees[[k]][['tree']], tmp)
							tmp	<- adephylo:::distTips( tmp )
							tmp	<- as.matrix(tmp)
							tmp	<- as_tibble(melt(tmp, as.is=TRUE, varnames=c('H1_TAXA','H2_TAXA'), value.name = "PD"))
							#	calculate mean tip to tip distance between taxa from male and female
							tmp	<- tmp	%>%
									mutate(	H1= gsub('^([^_]+).*','\\1', H1_TAXA),
											H2= gsub('^([^_]+).*','\\1', H2_TAXA) ) %>%
									filter( H1<H2) %>%
									summarise(MPD=mean(PD)) %>%
									as.double()
							tmp	<- dc %>%
									mutate(	MPD:= tmp, 
											MPD_STD:= tmp/phyloscanner.trees[[k]][['normalisation.constant']],
											W_FROM:= phyloscanner.trees[[k]][['window.coords']][['start']],
											W_TO:= phyloscanner.trees[[k]][['window.coords']][['end']])				
							# return 
							ans[[length(ans)+1L]]	<- tmp
						}				
					}
				ans <- do.call('rbind',ans)
				ans
			})
	dwin <- do.call('rbind',dwin)
	#	 for each pair only keep results from run with most windows
	tmp <- dwin %>% 
			group_by(IDCLU, PTY_RUN, H1, H2) %>% 
			summarise(N= length(W_FROM)) %>%
			group_by(IDCLU, H1, H2) %>%
			summarise(PTY_RUN:= PTY_RUN[which.max(N)[1]])
	dwin <- tmp %>% inner_join(dwin, by=c('IDCLU','H1','H2','PTY_RUN')) 	
	#	add pairwise combinations without phyloscanner evaluation as
	#	missing
	dwin <- dnet2 %>% 
			select(-PTY_RUN) %>% 
			distinct() %>%
			full_join(dwin, by=c('IDCLU','H1','H2')) 
	
	#	get most likely transmission chain for selected transmission networks
	dchain <- dchain %>% filter( IDCLU %in% idclus) 

	#	get most likely transmission networks
	dnet <- dnet %>% filter( IDCLU %in% idclus)
	
	#	get meta-data for individuals
	hivc.db.Date2numeric<- function( x )
	{
		if(!class(x)%in%c('Date','character'))	return( x )
		x	<- as.POSIXlt(x)
		tmp	<- x$year + 1900
		x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
		x	
	}
	
	dmeta <- dnet2 %>% select(H1,H2) %>% gather() %>% distinct() %>% rename(AID=value) %>% select(-key)	 
	dmeta <- as_tibble(read.csv("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_anonymised_RIDs.csv", stringsAsFactors=FALSE)) %>%
			select(ID,AID) %>% 
			distinct() %>%
			right_join(dmeta, by='AID') %>%
			rename(RID:=ID)
	stopifnot( 0==dmeta %>% filter(is.na(RID)) %>% nrow() )	
	infile <- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender_age.rda'
	load(infile)	
	dmeta <- dmeta %>% 
			left_join(as_tibble(rd), by='RID') %>% 
			select(RID, AID, PID, SEX, BIRTHDATE, LASTNEGDATE, FIRSTPOSDATE, VISIT, DATE, COMM_NUM_A, ARVSTARTDATE, FIRSTSELFREPORTART, EST_DATEDIED) %>%
			rename(VISIT_ROUND:= VISIT, VISIT_DATE:= DATE)
	tmp <- as_tibble(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data shared with consortium/Pangea_Rakai_RCCS_Metadata_11Dec2015.csv', stringsAsFactors=FALSE)) %>%
			select(studyid, Pangea.id, sampleDate) %>%
			rename(RID:= studyid, PID:= Pangea.id, SEQ_DATE:=sampleDate) %>%
			mutate(SEQ_DATE=hivc.db.Date2numeric(SEQ_DATE))
	dmeta <- dmeta %>% left_join(tmp, by=c('RID','PID'))
 	dmeta <- dmeta %>% select(-c(RID,PID))
	
	#	save
	outfile <- '~/Box Sync/OR_Work/MSCs/2019_Melodie/data/Rakai_phscnetworks_190706_selected_for_networkinference.rda'
	save(dnet, dchain, dwin, dmeta, file=outfile)
	
}

#' Extract the phyloscanner summary stats for each individual to estimate infection times
participants.phyloscanner.summarystats<- function()
{
	require(data.table)
	require(tidyverse)
	require(ape)
	require(phytools)
	
	#	locally
	if(1)
	{
		indir	<- '~/Box Sync/OR_Work/MSCs/2019_Melodie'
		#infile.couples <- file.path(indir,'analyses','190611_phsccouples.rda')	
		rel.dir <- "~/Box Sync/OR_Work/MSCs/2019_Melodie/data/RakaiPopSample_phsc_out190626"
	}
	
	#	list relationship files
	df <- tibble(F=list.files(rel.dir)) %>% 
			filter(grepl('workspace.rda',F)) %>%
			mutate(PTY_RUN:= as.integer(gsub('ptyr([0-9]+)_(.*)','\\1', F))) 
	
	ans	<- vector('list', 0)
	for(run in df$PTY_RUN)
	{
		print(run)
		#run	<- 5	#for devel
		tmpfile	<- df %>% filter(PTY_RUN==run) %>% transmute(F:= file.path(rel.dir, F)) %>% as.character()
		#	load phyloscanner output file and extract mean tip to tip distances		
		load(tmpfile)
		ans[[length(ans)+1L]] <- summary.stats
	}
	
	for(run in seq_along(ans))
		ans[[run]]<- dplyr:::select(ans[[run]],setdiff(1:ncol(ans[[run]]),contains('prop.gp')))
	
	ans	<- do.call('rbind', ans)
	
	save(ans, file=file.path(indir,'analyses','190627_phsccouples_summary_stats.rda'))
}

#' Find mean patristic tip to tip distances for RCCS couples
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
	df	<- unique(subset(rcl, select=PTY_RUN)) %>%
			inner_join(df, by='PTY_RUN')
	df	<- tibble(FO=list.files(file.path(indir,'analyses'), pattern='mpd_[0-9]+')) %>%
			mutate(PTY_RUN:= as.integer(gsub('^.*_mpd_([0-9]+).*$','\\1', FO))) %>%
			full_join(df, by='PTY_RUN')			
	df	<- df %>% 
			filter(is.na(FO)) %>%
			arrange(PTY_RUN) 
	df	<- df %>%
			mutate(CASE_ID:= seq_len(nrow(df))) %>%
			mutate(CASE_ID2:= floor(CASE_ID/10))		
	df	<- df	<- df %>% filter(CASE_ID2==0)
	
	df	<- as.data.table(df)
	#	for each run, unzip and calculate tip-to-tip distances among linked couples	
	for(run in df$PTY_RUN)
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

networks.with.high.prop.women<- function()
{
	infile			<- "~/sandbox/DeepSeqProjects/RakaiPopSample_data/Dataset_S2.csv"
	dmeta			<- as_tibble(read.csv(infile, stringsAsFactors=FALSE))
	infile			<- '~/sandbox/DeepSeqProjects/RakaiPopSample_phyloscanner_analysis_190706/Rakai_phscnetworks_190706.rda'
	load(infile)
	
	#	select networks of size > 4 with at most one man 
	#	or size 4 all women
	#	extract females 
	dfem <- dnet %>% 
			filter(IDCLU%in%c(4, 8, 27, 29, 30, 33, 36, 40, 68)) %>% 
			select(H1,H2) %>% 
			gather('key','AID') %>%
			select(AID) %>%
			distinct() %>%
			filter(grepl('F$',AID))
	#	de-anonymise
	infile <- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_anonymised_RIDs.csv'
	da	<- as_tibble(read.csv(infile, stringsAsFactors=FALSE))
	dfem <- dfem %>% left_join(da, by='AID')
	
	#	
	outfile <- '~/Box Sync/OR_Work/2019/2019_RCCS_high_risk_women/191011_femaleonly_primary_selection.csv'
	write.csv(dfem, file=outfile)
	
	#	select networks of size > 3 with at most one man 	
	#	extract females
	dfem <- dnet %>% 
			filter(IDCLU%in%c(50, 54, 60, 61, 64, 65, 69, 72 )) %>% 
			select(H1,H2) %>% 
			gather('key','AID') %>%
			select(AID) %>%
			distinct() %>%
			filter(grepl('F$',AID))
	#	de-anonymise
	infile <- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_anonymised_RIDs.csv'
	da	<- as_tibble(read.csv(infile, stringsAsFactors=FALSE))
	dfem <- dfem %>% left_join(da, by='AID')
	
	#	
	outfile <- '~/Box Sync/OR_Work/2019/2019_RCCS_high_risk_women/191011_femaleonly_backup_selection.csv'
	write.csv(dfem, file=outfile)
	
}