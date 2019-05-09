phsc.migrate.definestage1runs<- function()
{
	require(Phyloscanner.R.utilities)
	require(phyloscannerR)
	indir	<- '~/Box Sync/OR_Work/2018/2018_MunichCluster/BamData'
	dbam	<- tibble(BAM=list.files(indir, pattern='*bam$',full.names=FALSE, recursive=TRUE)) 
	regex.person	<- '^([A-Z0-9]+-[A-Z0-9]+).*$'
	dbam	<- dbam %>% mutate( IND:=  gsub(regex.person,'\\1',BAM) )
	phyloscannerR::define.stage1.analyses(dbam$IND, batch.size=5L)
	
	'~/duke/2018_MunichCluster/Data_180815/15-01402_PRRT.bam'
}

phsc.migrate.findbams<- function()
{
	require(Phyloscanner.R.utilities)
	require(phyloscannerR)
	indir	<- '~/Box Sync/OR_Work/2018/2018_MunichCluster/BamData'
	dbam	<- tibble(BAM=list.files(indir, pattern='*bam$',full.names=FALSE, recursive=TRUE)) 
	regex.person	<- '^([A-Z0-9]+-[A-Z0-9]+).*$'
	dbam	<- dbam %>% mutate( IND:=  gsub(regex.person,'\\1',BAM) )
	phyloscannerR::define.stage1.analyses(dbam$IND, batch.size=5L)
	
	data.dir	<- '~/duke/2018_MunichCluster/Data_180815'	
	regex.person='^([A-Za-z0-9]+-[0-9]+)_.*$'
	regex.bam='^(.*)\\.bam$'
	regex.ref='^(.*)\\.fasta$' 
	verbose=1
	
	tmp	<- phsc.find.bam.and.references(data.dir, regex.person=regex.person, regex.bam='^(.*)\\.bam$', regex.ref='^(.*)_ref\\.fasta$', verbose=1) 
	
}

phsc.migrate.counts.Rakai.branch.mh_devel2<- function()
{
	require(Phyloscanner.R.utilities)
	require(phyloscannerR)
	require(tidyverse)
	
	#	load phyloscanner data	
	file	<- system.file(file.path('extdata','ptyr192_phsc_analyse_trees_output.RData'),package='phyloscannerR')
	load(file)	#loads 'phsc', output from 'phyloscanner.analyse.trees'
	#	use distance thresholds found in analysis of Rakai couples
	close.threshold	<- 0.025
	distant.threshold	<- 0.05	
	#	use relationship types based on adjacency
	#	this also considers linkage etc between individuals who have dual infections, recombinants etc 
	#	..and thus may not have *all* their subgraphs adjacent to each other
	relationship.types	<- c('proximity.3.way',
			'close.and.adjacent',					
			'close.and.adjacent.and.directed',					
			'close.and.adjacent.and.ancestry',			
			'close.and.contiguous',					
			'close.and.contiguous.and.directed',					
			'close.and.contiguous.and.ancestry')		
	dwin	<- classify.pairwise.relationships(phsc, 			 
			close.threshold=close.threshold, 
			distant.threshold=distant.threshold,			
			relationship.types=relationship.types, 
			verbose=TRUE)		
	
	file <- '~/sandbox/DeepSeqProjects/RakaiPopSample_migcode/ptyr192_pairwise_relationships_branch_migrateCounts.rda'
	load(file)
	dwino <- copy(migcnts$dwin)
	
	# check basic classification
	tmp	<- dwin %>% 
			select(host.1, host.2, tree.id, basic.classification, adjacent, contiguous, paths12, paths21, ancestry) %>%
			rename(basic.classification.mhdevel2= basic.classification,
					adjacent.mhdevel2=adjacent, 
					contiguous.mhdevel2=contiguous, 
					paths12.mhdevel2=paths12, 
					paths21.mhdevel2=paths21,
					ancesctry.mhdevel2=ancestry)
	tmp2 <- dwino %>% 
			select(host.1, host.2, tree.id, basic.classification, adjacent, contiguous, paths12, paths21, ancestry) %>%
			rename(basic.classification.migcnts= basic.classification,
					adjacent.migcnts=adjacent, 
					contiguous.migcnts=contiguous, 
					paths12.migcnts=paths12, 
					paths21.migcnts=paths21,
					ancesctry.migcnts=ancestry)
	
	tmp	<- tmp %>% 
			full_join(tmp2, by=c('host.1','host.2','tree.id')) 
	tmp %>%
		filter(basic.classification.mhdevel2!=basic.classification.migcnts)
	tmp %>%
		group_by(basic.classification.mhdevel2, basic.classification.migcnts) %>%
		summarise(k = n()) 

}

phsc.migrate.counts.Rakai.branch.migrateCounts<- function()
{
	require(Phyloscanner.R.utilities)
	require(phyloscannerR)
	require(tidyverse)	
	
	#	load phyloscanner data	
	file	<- system.file(file.path('extdata','ptyr192_phsc_analyse_trees_output.R'),package='phyloscannerR')
	load(file)	#loads 'phsc', output from 'phyloscanner.analyse.trees'
	#	use distance thresholds found in analysis of Rakai couples
	close.threshold	<- 0.025
	distant.threshold	<- 0.05	
	#	use relationship types based on adjacency
	#	this also considers linkage etc between individuals who have dual infections, recombinants etc 
	#	..and thus may not have *all* their subgraphs adjacent to each other
	relationship.types	<- c('proximity.3.way',
			'any.ancestry',
			'close.x.contiguous',			
			'close.and.adjacent',					
			'close.and.adjacent.and.directed',					
			'close.and.adjacent.and.ancestry',			
			'close.and.contiguous',					
			'close.and.contiguous.and.directed',					
			'close.and.contiguous.and.ancestry')	
	use.paths.to.define.relationships <- TRUE
	dwin	<- classify.pairwise.relationships(phsc, 
								allow.mt=TRUE, 
								close.threshold=close.threshold, 
								distant.threshold=distant.threshold,
								use.paths.to.define.relationships=use.paths.to.define.relationships,
								relationship.types=relationship.types, 
								verbose=TRUE)		
	# select windows by read and tip count
	tip.regex <- "^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$"
	min.reads <- 30
	min.tips <- 1
	dwin <- select.windows.by.read.and.tip.count(phsc, dwin, tip.regex, min.reads, min.tips)
	# count phylogenetic relationships across deep-sequence trees
	dc <- count.pairwise.relationships(dwin)
	
	migcnts	<- list(dc=dc, dwin=dwin)
	file <- '~/sandbox/DeepSeqProjects/RakaiPopSample_migcode/ptyr192_pairwise_relationships_branch_migrateCounts.rda'
	save(migcnts, file=file)
	
	
	
	dwinn <- copy(dwin)
	require(data.table)
	file <- '~/sandbox/DeepSeqProjects/RakaiPopSample_phyloscanner_analysis/ptyr192_pairwise_relationships.rda'
	load(file)
	tmp	<- dwinn %>% 
			select(host.1,host.2,tree.id,basic.classification,categorical.distance)
	#%>%
	#		mutate(TYPE_BASIC:= paste0(basic.classification,'_',categorical.distance)) %>%
	#		select(-basic.classification,-categorical.distance)
	
	tmp2 <- subset(dwin, select=c(ID1,ID2,SUFFIX,TYPE_BASIC))
	set(tmp2, NULL, 'ID1', tmp2[, as.character(ID1)])
	set(tmp2, NULL, 'ID2', tmp2[, as.character(ID2)])
	tmp2[, BASIC_CLASSIFICATION:= gsub('\n| ','_',gsub('chain 21','desc',gsub('chain 12','anc',gsub('with intermediate','noncontiguous',gsub('no intermediate','contiguous',gsub('\nclose|\ndistant| close| distant','',TYPE_BASIC))))))]
	tmp2[, CATEGORICAL_DISTANCE:= 'intermediate']
	set(tmp2, tmp2[, which(grepl('close',TYPE_BASIC))], 'CATEGORICAL_DISTANCE', 'close')
	set(tmp2, tmp2[, which(grepl('distant',TYPE_BASIC))], 'CATEGORICAL_DISTANCE', 'distant')
	setnames(tmp2, c('ID1','ID2','SUFFIX'), c('host.1','host.2','tree.id'))	
	tmp	<- tmp %>% full_join(tmp2, by=c('host.1','host.2','tree.id'))
	
	tmp %>% filter( basic.classification!=BASIC_CLASSIFICATION )
 	#
	#	check that results coincide with old code 
	ans	<- multinomial.calculations(phsc, 
								close.threshold=close.threshold,
								tip.regex=tip.regex, 
								allow.mt=TRUE, 
								min.reads=min.reads, 
								min.tips=min.tips, 
								distant.threshold = distant.threshold,
								relationship.types	= c('proximity.3.way',
														'any.ancestry',
														'close.x.contiguous',
														'close.and.contiguous',
														'close.and.contiguous.and.directed',
														'close.and.adjacent.and.directed',
														'close.and.contiguous.and.ancestry.cat',
														'close.and.adjacent.and.ancestry.cat',
														'adjacent.and.proximity.cat'),
								verbose=TRUE)
	dwino<- ans$dwin
	dco<- ans$rplkl
	
	# counts
	tmp		<- dco %>% 
		select(host.1,host.2,categorisation,type,n,n.eff,k,k.eff) %>%
		filter(categorisation=='basic.classification') %>%			
		mutate(basic.classification:= gsub('noancestry','sibling',gsub('complex','intermingled',type))) %>%
		select(-categorisation, -type)	
	tmp2	<- dc %>% 
		select(host.1,host.2,categorisation,categorical.distance,type,n,n.eff,k,k.eff) %>%
		filter(categorisation=='basic.classification') %>%			
		mutate(basic.classification:= paste0(type,'_',categorical.distance)) %>%
		select(-categorisation, -categorical.distance, -type)
	tmp		<- tmp2 %>% inner_join(tmp, by=c('host.1','host.2','basic.classification'))
	stopifnot( tmp %>% filter(n.x!=n.y) %>% nrow() == 0 )
	stopifnot( tmp %>% filter(n.eff.x!=n.eff.y) %>% nrow() == 0 )
	stopifnot( tmp %>% filter(k.x!=k.y) %>% nrow() == 0 )
	stopifnot( tmp %>% filter(k.eff.x!=k.eff.y) %>% nrow() == 0 )
	
	# windows all OK
	tmp	<- dwin %>% select(host.1,host.2,tree.id,basic.classification)
	tmp2<- dwino %>% select(host.1,host.2,tree.id,basic.classification) %>% mutate(basic.classification=gsub('_close|_distant|_intermediate','',basic.classification))
	tmp <- tmp %>% inner_join(tmp2, by=c('host.1','host.2','tree.id'))
	table(tmp$basic.classification.x,tmp$basic.classification.y)
	
	tmp	<- dwin %>% select(host.1,host.2,tree.id,proximity.3.way.cat)
	tmp2<- dwino %>% select(host.1,host.2,tree.id,proximity.3.way)
	tmp <- tmp %>% inner_join(tmp2, by=c('host.1','host.2','tree.id'))
	table(tmp$proximity.3.way.cat,tmp$proximity.3.way)
	
	tmp	<- dwin %>% select(host.1,host.2,tree.id,any.ancestry.cat)
	tmp2<- dwino %>% select(host.1,host.2,tree.id,any.ancestry)
	tmp <- tmp %>% inner_join(tmp2, by=c('host.1','host.2','tree.id'))
	table(tmp$any.ancestry.cat,tmp$any.ancestry)
	
	tmp	<- dwin %>% select(host.1,host.2,tree.id,close.x.contiguous.cat)
	tmp2<- dwino %>% select(host.1,host.2,tree.id,close.x.contiguous)
	tmp <- tmp %>% inner_join(tmp2, by=c('host.1','host.2','tree.id'))
	table(tmp$close.x.contiguous.cat,tmp$close.x.contiguous)
	
	tmp	<- dwin %>% select(host.1,host.2,tree.id,close.and.contiguous.cat)
	tmp2<- dwino %>% select(host.1,host.2,tree.id,close.and.contiguous)
	tmp <- tmp %>% inner_join(tmp2, by=c('host.1','host.2','tree.id'))
	table(tmp$close.and.contiguous.cat,tmp$close.and.contiguous)
	
	tmp	<- dwin %>% select(host.1,host.2,tree.id,close.and.contiguous.and.directed.cat)
	tmp2<- dwino %>% select(host.1,host.2,tree.id,close.and.contiguous.and.directed)
	tmp <- tmp %>% inner_join(tmp2, by=c('host.1','host.2','tree.id'))
	table(tmp$close.and.contiguous.and.directed.cat,tmp$close.and.contiguous.and.directed)
	
	tmp	<- dwin %>% select(host.1,host.2,tree.id,close.and.adjacent.and.directed.cat)
	tmp2<- dwino %>% select(host.1,host.2,tree.id,close.and.adjacent.and.directed)
	tmp <- tmp %>% inner_join(tmp2, by=c('host.1','host.2','tree.id'))
	table(tmp$close.and.adjacent.and.directed.cat,tmp$close.and.adjacent.and.directed)
	
	tmp	<- dwin %>% select(host.1,host.2,tree.id,close.and.contiguous.and.ancestry.cat)
	tmp2<- dwino %>% select(host.1,host.2,tree.id,close.and.contiguous.and.ancestry.cat)
	tmp <- tmp %>% inner_join(tmp2, by=c('host.1','host.2','tree.id'))
	table(tmp$close.and.contiguous.and.ancestry.cat.x,tmp$close.and.contiguous.and.ancestry.cat.y)
	
	tmp	<- dwin %>% select(host.1,host.2,tree.id,close.and.adjacent.and.ancestry.cat)
	tmp2<- dwino %>% select(host.1,host.2,tree.id,close.and.adjacent.and.ancestry.cat)
	tmp <- tmp %>% inner_join(tmp2, by=c('host.1','host.2','tree.id'))
	table(tmp$close.and.adjacent.and.ancestry.cat.x,tmp$close.and.adjacent.and.ancestry.cat.y)
	# windows all OK
	
}

phsc.migrate.analyzetrees.Rakai<- function()
{
	require(Phyloscanner.R.utilities)
	require(phyloscannerR)
	?phyloscanner.analyse.tree	
	
	
	#	extract RCCS example data
	tree.file.zip				<- system.file(file.path('extdata','Rakai_run192_trees.zip'),package='phyloscannerR')
	tree.file.directory			<- '~/Box Sync/OR_Work/2018/2018_Rakai_192'	
	unzip(tree.file.zip, exdir=tree.file.directory, junkpaths=TRUE)
	
	#	arguments as used for RCCS analysis
	file.name.regex				<- "^\\D*([0-9]+)_to_([0-9]+)\\D*$"
	max.reads.per.host			<- 50
	multifurcation.threshold 	<- 1e-5
	norm.ref.file.name			<- system.file('HIV_DistanceNormalisationOverGenome.csv',package='phyloscannerR')	
	outgroup.name				<- "REF_CPX_AF460972"
	raw.blacklist.threshold		<- 20
	sankoff.k					<- 20
	sankoff.unassigned.switch.threshold	<- 0
	seed						<- 42
	splits.rule					<- 's'
	relaxed.ancestry			<- TRUE
	tip.regex					<- "^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$"
	tree.file.regex				<- "^ptyr192_InWindow_([0-9]+_to_[0-9]+)\\.tree$"
	verbosity					<- 1
	
	phsc	<- phyloscanner.analyse.trees(tree.file.directory,
			tree.file.regex = tree.file.regex,
			splits.rule = splits.rule, 
			sankoff.k = sankoff.k,
			sankoff.unassigned.switch.threshold = sankoff.unassigned.switch.threshold,
			outgroup.name = outgroup.name,
			multifurcation.threshold = multifurcation.threshold, 
			guess.multifurcation.threshold = FALSE,
			user.blacklist.directory = NULL, 
			user.blacklist.file.regex = NULL,
			duplicate.file.directory = NULL,
			duplicate.file.regex = NULL,
			recombination.file.directory = NULL,
			recombination.file.regex = NULL,
			alignment.file.directory = NULL, 
			alignment.file.regex = NULL,
			tip.regex = tip.regex,
			file.name.regex = file.name.regex,
			seed = seed, 
			relaxed.ancestry = relaxed.ancestry,
			norm.ref.file.name = NULL,
			norm.standardise.gag.pol = TRUE, 
			norm.constants = NULL,
			parsimony.blacklist.k = sankoff.k, 
			raw.blacklist.threshold = raw.blacklist.threshold,
			ratio.blacklist.threshold = 0, 
			do.dual.blacklisting = FALSE,
			max.reads.per.host = max.reads.per.host, 
			blacklist.underrepresented = FALSE,
			use.ff = FALSE, 
			prune.blacklist = FALSE, 
			count.reads.in.parsimony = TRUE,
			verbosity = verbosity, 
			no.progress.bars = FALSE)
	
	#	save phyloscanner output
	save(phsc, file=file.path(tree.file.directory,'ptyr192_phsc_analyse_trees_output.RData'))	
}


phsc.migrate.analyzetrees.BEEHIVE<- function()
{
	require(Phyloscanner.R.utilities)
	require(phyloscannerR)
	?phyloscanner.analyse.tree	
	
	indir						<- '/Users/Oliver/sandbox/DeepSeqProjects/BEEHIVE_mbeexample_input'
	outdir						<- '/Users/Oliver/sandbox/DeepSeqProjects/BEEHIVE_mbeexample_output'	
	tree.file.directory			<- file.path(indir,'RAxMLfiles')
	alignment.file.directory	<- file.path(indir,'AlignedReads')
	recombination.file.directory<- file.path(indir,'RecombinationData')
			
	raw.blacklist.threshold		<- 10
	splits.rule					<- 's'
	sankoff.k					<- 20
	tree.file.regex				<- "^RAxML_bestTree.InWindow_([0-9]+_to_[0-9]+)\\.tree$"
	alignment.file.regex		<- "^AlignedReads_AsUsed_InWindow_([0-9]+_to_[0-9]+).fasta$"
	recombination.file.regex	<- "^RecombinantReads_InWindow_([0-9]+_to_[0-9]+).csv$"
	tip.regex					<- "^(.*)_read_([0-9]+)_count_([0-9]+)$"
	file.name.regex				<- "^\\D*([0-9]+)_to_([0-9]+)\\D*$"
	seed						<- 42
	outgroup.name				<- "C.BW.00.00BW07621.AF443088"
	
	phsc	<- phyloscanner.analyse.trees(tree.file.directory,
			tree.file.regex = tree.file.regex,
			splits.rule = splits.rule, 
			sankoff.k = 20,
			sankoff.unassigned.switch.threshold = 0,
			outgroup.name = outgroup.name,
			multifurcation.threshold = 1e-5, 
			guess.multifurcation.threshold = FALSE,
			user.blacklist.directory = NULL, 
			user.blacklist.file.regex = NULL,
			duplicate.file.directory = NULL,
			duplicate.file.regex = "^DuplicateReadCountsProcessed_InWindow_([0-9]+_to_[0-9]+).csv$",
			recombination.file.directory = recombination.file.directory,
			recombination.file.regex = recombination.file.regex,
			alignment.file.directory = NULL, 
			alignment.file.regex = alignment.file.regex,
			tip.regex = tip.regex,
			file.name.regex = file.name.regex,
			seed = seed, 
			norm.ref.file.name = NULL,
			norm.standardise.gag.pol = FALSE, 
			norm.constants = NULL,
			parsimony.blacklist.k = sankoff.k, 
			raw.blacklist.threshold = raw.blacklist.threshold,
			ratio.blacklist.threshold = 0, 
			do.dual.blacklisting = FALSE,
			max.reads.per.host = Inf, 
			blacklist.underrepresented = FALSE,
			use.ff = FALSE, 
			prune.blacklist = FALSE, 
			count.reads.in.parsimony = TRUE,
			verbosity = 0, 
			no.progress.bars = FALSE)		
}

