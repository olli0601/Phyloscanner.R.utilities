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

phsc.migrate.phyloscan.plot<- function()
{
	require(Phyloscanner.R.utilities)
	require(phyloscannerR)
	require(tidyverse)
	
	#	load phyloscanner data	
	file	<- system.file(file.path('extdata','ptyr192_phsc_analyse_trees_output.RData'),package='phyloscannerR')
	load(file)	#loads 'phsc', output from 'phyloscanner.analyse.trees'	
	tmp	<- produce.pairwise.graphs(phsc, 
				0.025,	
				hosts=c('RkA05868F','RkA05968M','RkA00369F','RkA01344M'),
				contiguous.pairs = FALSE,
				inclusion = "both")
		
	
	hosts	<- c('RkA05868F','RkA05968M','RkA00369F','RkA01344M')
	inclusion <- "both"# "either"
	tmp	<- produce.pairwise.graphs2(phsc, hosts=hosts, inclusion = "both")
	tmp$graph
	
	
	table( dwin$basic.topology )
}
	
phsc.migrate.counts.Rakai.classify.pairwise.relationships.branch.mh_devel2<- function()
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

	# check close.and.adjacent ...
	tmp	<- dwin %>% 
			select(host.1, host.2, tree.id, close.and.adjacent.cat, close.and.adjacent.and.directed.cat, close.and.adjacent.and.ancestry.cat, adjacent, contiguous, paths12, paths21, ancestry) %>%
			rename(	close.and.adjacent.mhdevel2= close.and.adjacent.cat,
				close.and.adjacent.and.directed.mhdevel2=close.and.adjacent.and.directed.cat,
				close.and.adjacent.and.ancestry.mhdevel2=close.and.adjacent.and.ancestry.cat,
				adjacent.mhdevel2=adjacent, 
				contiguous.mhdevel2=contiguous, 
				paths12.mhdevel2=paths12, 
				paths21.mhdevel2=paths21,
				ancesctry.mhdevel2=ancestry)
	tmp2 <- dwino %>% 
				select(host.1, host.2, tree.id, close.and.adjacent.cat, close.and.adjacent.and.directed.cat, close.and.adjacent.and.ancestry.cat, adjacent, contiguous, paths12, paths21, ancestry) %>%
				rename( close.and.adjacent.migcnts= close.and.adjacent.cat,
						close.and.adjacent.and.directed.migcnts=close.and.adjacent.and.directed.cat,
						close.and.adjacent.and.ancestry.migcnts=close.and.adjacent.and.ancestry.cat,				
						adjacent.migcnts=adjacent, 
						contiguous.migcnts=contiguous, 
						paths12.migcnts=paths12, 
						paths21.migcnts=paths21,
						ancesctry.migcnts=ancestry)
	tmp	<- tmp %>% 
		full_join(tmp2, by=c('host.1','host.2','tree.id'))
	tmp %>%
		filter(close.and.adjacent.mhdevel2!=close.and.adjacent.migcnts)
	tmp %>%
		filter(close.and.adjacent.and.directed.mhdevel2!=close.and.adjacent.and.directed.migcnts)
	tmp %>%
		filter(close.and.adjacent.and.ancestry.mhdevel2!=close.and.adjacent.and.ancestry.migcnts)

	# check close.and.contiguous ...
	tmp	<- dwin %>% 
			select(host.1, host.2, tree.id, close.and.contiguous.cat, close.and.contiguous.and.directed.cat, close.and.contiguous.and.ancestry.cat, adjacent, contiguous, paths12, paths21, ancestry) %>%
			rename(	close.and.contiguous.mhdevel2= close.and.contiguous.cat,
					close.and.contiguous.and.directed.mhdevel2=close.and.contiguous.and.directed.cat,
					close.and.contiguous.and.ancestry.mhdevel2=close.and.contiguous.and.ancestry.cat,
					adjacent.mhdevel2=adjacent, 
					contiguous.mhdevel2=contiguous, 
					paths12.mhdevel2=paths12, 
					paths21.mhdevel2=paths21,
					ancesctry.mhdevel2=ancestry)
	tmp2 <- dwino %>% 
			select(host.1, host.2, tree.id, close.and.contiguous.cat, close.and.contiguous.and.directed.cat, close.and.contiguous.and.ancestry.cat, adjacent, contiguous, paths12, paths21, ancestry) %>%
			rename( close.and.contiguous.migcnts= close.and.contiguous.cat,
					close.and.contiguous.and.directed.migcnts=close.and.contiguous.and.directed.cat,
					close.and.contiguous.and.ancestry.migcnts=close.and.contiguous.and.ancestry.cat,				
					adjacent.migcnts=adjacent, 
					contiguous.migcnts=contiguous, 
					paths12.migcnts=paths12, 
					paths21.migcnts=paths21,
					ancesctry.migcnts=ancestry)
	tmp	<- tmp %>% 
			full_join(tmp2, by=c('host.1','host.2','tree.id'))
	tmp %>%
			filter(close.and.contiguous.mhdevel2!=close.and.contiguous.migcnts)
	tmp %>%
			filter(close.and.contiguous.and.directed.mhdevel2!=close.and.contiguous.and.directed.migcnts)
	tmp %>%
			filter(close.and.contiguous.and.ancestry.mhdevel2!=close.and.contiguous.and.ancestry.migcnts)
	
	#	patristic distance
	tmp	<- dwin %>% 
		select(host.1, host.2, tree.id, patristic.distance, categorical.distance) %>%
		rename(	patristic.distance.mhdevel2= patristic.distance,
				categorical.distance.mhdevel2=categorical.distance)
	tmp2<- dwin %>% 
		select(host.1, host.2, tree.id, patristic.distance, categorical.distance) %>%
		rename(	patristic.distance.migcnts= patristic.distance,
				categorical.distance.migcnts=categorical.distance)
	tmp	<- tmp %>% 
		full_join(tmp2, by=c('host.1','host.2','tree.id'))
	tmp %>%
		filter(categorical.distance.mhdevel2!=categorical.distance.migcnts)
	tmp %>%
		filter(patristic.distance.mhdevel2!=patristic.distance.migcnts)

}

phsc.migrate.counts.Rakai.classify.pairwise.relationships.branch.migrateCounts<- function()
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

phsc.migrate.transmission.networks<- function()
{
	#	get meta-data for individuals in general pop cohort
	infile			<- "~/sandbox/DeepSeqProjects/RakaiPopSample_data/Dataset_S2.csv"
	dmeta			<- as_tibble(read.csv(infile, stringsAsFactors=FALSE))
	
	indir <-  '~/sandbox/DeepSeqProjects/RakaiPopSample_phyloscanner_analysis_190706'
	outfile <- '~/sandbox/DeepSeqProjects/RakaiPopSample_phyloscanner_analysis_190706/Rakai_phscnetworks_allpairs_190706.rda'
	control <- list(linked.group='close.and.adjacent.cat',linked.no='not.close.or.nonadjacent',linked.yes='close.and.adjacent', conf.cut=0.6, neff.cut=3)
	tmp <- find.pairs.in.networks(indir, batch.regex='^ptyr([0-9]+)_.*', control=control, verbose=TRUE, dmeta=dmeta)
	dpl <- copy(tmp$linked.pairs)
	dc <- copy(tmp$relationship.counts)
	dw <- copy(tmp$windows)
	save(dpl, dc, dw, file=outfile)
	
	infile <- '~/sandbox/DeepSeqProjects/RakaiPopSample_phyloscanner_analysis_190706/Rakai_phscnetworks_allpairs_190706.rda'
	outfile <- '~/sandbox/DeepSeqProjects/RakaiPopSample_phyloscanner_analysis_190706/Rakai_phscnetworks_190706.rda'
	load(infile)	
	tmp <- find.networks(dc, control=control, verbose=TRUE)
	dnet <- copy(tmp$transmission.networks)
	dchain <- copy(tmp$most.likely.transmission.chains)
	save(dpl, dc, dw, dnet, dchain, file=outfile)
	
	
	
	infile			<- "~/sandbox/DeepSeqProjects/RakaiPopSample_data/Dataset_S2.csv"
	dmeta			<- as_tibble(read.csv(infile, stringsAsFactors=FALSE))
	infile			<- '~/sandbox/DeepSeqProjects/RakaiPopSample_phyloscanner_analysis_190706/Rakai_phscnetworks_190706.rda'
	load(infile)
	idclus <- sort(unique(dnet$IDCLU))
	di <- copy(dmeta)
	setnames(di, 'ID', 'H')
	df <- dnet %>% 
			filter(IDCLU == 34) %>%
			select(-c(H1_SEX,H2_SEX))		
	control<- list()
	control$point.size = 10
	control$edge.gap = 0.04
	control$edge.size = 2
	control$curvature = -0.2
	control$arrow = arrow(length = unit(0.04, "npc"), type = "open")
	control$curv.shift = 0.06
	control$label.size = 3
	control$node.label = "H" 
	control$node.fill = "SEX"
	control$node.shape = NA_character_
	control$node.shape.values = c(`NA` = 16)
	control$node.fill.values = c(F = "hotpink2", M = "steelblue2")
	control$threshold.linked = 0.6
	p <- plot.network(df, di, control)	
	png(file = paste0(outfile.base, "_trmnetwork_34.png"), width = 6, height = 6, 
			units = "in", res = 400)
	print(p)
	dev.off()
	
	
	control$layout <- p$layout 
	di <- copy(dmeta)			
	setnames(di, 'ID', 'H')
	df <- subset(dchain, IDCLU==34)	
	p2 <- plot.chain(df, di, control)
	
	
	png(file=paste0(outfile.base,'_trmchain_34.png'), width=6, height=6, units='in', res=400)		
	print(p2)
	dev.off()
	
}


phsc.migrate.normalisation<- function()
{
	norm.ref.file.name			<- system.file('HIV_DistanceNormalisationOverGenome.csv',package='phyloscannerR')
	dn	<- as.data.table(read.csv(norm.ref.file.name))
	ggplot(dn, aes(x=POSITION_WRT_HXB2, y=log10(MEDIAN_PAIRWISE_DISTANCE_BETWEEN_STANDARD_REFS))) +
			geom_point()
}

phsc.Rakai.analyzetrees.stage2<- function()
{	
	require(tidyverse)
	require(data.table)
	require(phyloscannerR)
	
	if(0)
	{
		indir	<- '/Users/Oliver/sandbox/DeepSeqProjects/RakaiPopSample_deepseqtrees'
		tmpdir	<- '/Users/Oliver/sandbox/DeepSeqProjects/RakaiPopSample_phsc_work'
		outdir	<- '/Users/Oliver/sandbox/DeepSeqProjects/RakaiPopSample_phsc_out190512'
		prog.phyloscanner_analyse_trees <- '/Users/Oliver/git/phyloscanner/phyloscanner_analyse_trees.R'
		tree.dir <- "RakaiPopSample_deepseqtrees"
		tree.dir <- '/Users/Oliver/sandbox/DeepSeqProjects/RakaiPopSample_deepseqtrees'		
	}
	if(1)
	{
		indir	<- '/Users/Oliver/sandbox/DeepSeqProjects/RakaiPopSample_deepseqtrees'
		tmpdir	<- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS_tmp'
		outdir	<- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS_phsc190706'
		prog.phyloscanner_analyse_trees <- '/rds/general/user/or105/home/libs_sandbox/phyloscanner/phyloscanner_analyse_trees.R'
		tree.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS_trees'		
	}
	
	#	set phyloscanner variables
	#	arguments as used for RCCS analysis
	control	<- list()
	control$allow.mt <- TRUE				
	control$alignment.file.directory = NULL 
	control$alignment.file.regex = NULL
	control$blacklist.underrepresented = FALSE	
	control$count.reads.in.parsimony = TRUE
	control$distance.threshold <- '0.025 0.05'
	control$do.dual.blacklisting = FALSE					
	control$duplicate.file.directory = NULL
	control$duplicate.file.regex = NULL
	control$file.name.regex = "^\\D*([0-9]+)_to_([0-9]+)\\D*$"
	control$guess.multifurcation.threshold = FALSE
	control$max.reads.per.host = 50
	control$min.reads.per.host <- 30
	control$min.tips.per.host <- 1	
	control$multifurcation.threshold = 1e-5
	control$multinomial= TRUE
	control$norm.constants = NULL
	control$norm.ref.file.name = system.file('HIV_DistanceNormalisationOverGenome.csv',package='phyloscannerR')
	control$norm.standardise.gag.pol = TRUE
	control$no.progress.bars = TRUE
	control$outgroup.name = "REF_CPX_AF460972"
	control$output.dir = outdir
	control$parsimony.blacklist.k = 20
	control$prune.blacklist = FALSE
	control$post.hoc.count.blacklisting <- TRUE
	control$ratio.blacklist.threshold = 0 
	control$raw.blacklist.threshold = 20					
	control$recombination.file.directory = NULL
	control$recombination.file.regex = NULL
	control$relaxed.ancestry = TRUE
	control$sankoff.k = 20
	control$sankoff.unassigned.switch.threshold = 0
	control$seed = 42
	control$splits.rule = 's'
	control$tip.regex = "^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$"
	control$tree.file.regex = "^ptyr[0-9]+_InWindow_([0-9]+_to_[0-9]+)\\.tree$"
	control$use.ff = FALSE
	control$user.blacklist.directory = NULL 
	control$user.blacklist.file.regex = NULL
	control$verbosity = 1		
	
	
	#	make bash for many files	
	if(0)
	{
		tmp <- "https://datadryad.org/bitstream/handle/10255/dryad.208473/Dataset_S1.tar?sequence=1"
		download.file(tmp, destfile="Dataset_S1.tar", method="curl")
		untar("Dataset_S1.tar", exdir=tree.dir, extras='-xvf')		
	}
	df <- tibble(F=list.files(tree.dir))
	df <- df %>% 
			mutate(TYPE:= gsub('ptyr([0-9]+)_(.*)','\\2', F),
					RUN:= as.integer(gsub('ptyr([0-9]+)_(.*)','\\1', F))) %>%
			mutate(TYPE:= gsub('^([^\\.]+)\\.[a-z]+$','\\1',TYPE)) %>%
			spread(TYPE, F) %>%
			set_names(~ str_to_upper(.))	
	tmp	<- sort(as.integer(gsub('ptyr([0-9]+)_(.*)','\\1',list.files(outdir, pattern='_workspace.rda$'))))
	df <- df %>% filter(!RUN%in%tmp)
	valid.input.args <- cmd.phyloscanner.analyse.trees.valid.args(prog.phyloscanner_analyse_trees)
	cmds <- vector('list',nrow(df))
	for(i in seq_len(nrow(df)))
	{
		#	set input args
		control$output.string <- paste0('ptyr',df$RUN[i])	
		#	make script
		tree.input <- file.path(tree.dir, df$TREES_NEWICK[i])
		cmd <- cmd.phyloscanner.analyse.trees(prog.phyloscanner_analyse_trees, 
				tree.input, 
				control,
				valid.input.args=valid.input.args)
		cmds[[i]] <- cmd		
	}	
	cat(cmds[[100]])
	
	#
	# 	submit array job to HPC
	#
	#	make header
	hpc.load			<- "module load anaconda3/personal"	# make third party requirements available	 
	hpc.select			<- 1						# number of nodes
	hpc.nproc			<- 1						# number of processors on node
	hpc.walltime		<- 23						# walltime
	if(0)		
	{
		hpc.q			<- NA						# PBS queue
		hpc.mem			<- "48gb" 					# RAM		
	}
	#		or run this block to submit a job array to Oliver's machines
	if(1)
	{
		hpc.q			<- "pqeelab"				# PBS queue
		hpc.mem			<- "36gb" 					# RAM		
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
	
	
	stop()
	
	#	make bash for one file
	prog.phyloscanner_analyse_trees <- '/Users/Oliver/git/phyloscanner/phyloscanner_analyse_trees.R'
	valid.input.args <- cmd.phyloscanner.analyse.trees.valid.args(prog.phyloscanner_analyse_trees)
	tree.input <- system.file(file.path('extdata','Rakai_run192_trees.zip'),package='phyloscannerR')
	control$output.string <- 'Rakai_run192'
	cmd <- cmd.phyloscanner.analyse.trees(prog.phyloscanner_analyse_trees, 
			tree.input, 
			control,
			valid.input.args=valid.input.args)
	cat(cmd)
	
	
	#	locate tree and patient files for each run
	df	<- data.table(F=list.files(indir))	
	df[, TYPE:= gsub('ptyr([0-9]+)_(.*)','\\2', F)]
	set(df, NULL, 'TYPE', df[, gsub('^([^\\.]+)\\.[a-z]+$','\\1',TYPE)])
	df[, RUN:= as.integer(gsub('ptyr([0-9]+)_(.*)','\\1', F))]
	df	<- dcast.data.table(df, RUN~TYPE, value.var='F')
	setnames(df, colnames(df), toupper(colnames(df)))	
	#	call in R
	for(i in seq_len(nrow(df)))
	{
		#	file management: extract tree files, copy patient file		 
		tree.file.directory	<- file.path(tmpdir, paste0('ptyr',df[i,RUN]))
		unzip(file.path(indir, df[i, TREES_NEWICK]), exdir=tree.file.directory, junkpaths=TRUE)		
		file.copy(file.path(indir, df[i, PATIENTS]), tree.file.directory)
		
		control$output.string <- paste0('ptyr',df[i,RUN])
		#	run phyloscanner
		phsc	<- phyloscanner.analyse.trees(tree.file.directory,
						allow.mt= control$allow.mt,
						alignment.file.directory = control$alignment.file.directory, 
						alignment.file.regex = control$alignment.file.regex,
						blacklist.underrepresented = control$blacklist.underrepresented,
						count.reads.in.parsimony = control$count.reads.in.parsimony,
						do.dual.blacklisting = control$do.dual.blacklisting,
						duplicate.file.directory = control$duplicate.file.directory,
						duplicate.file.regex = control$duplicate.file.regex,
						file.name.regex = control$file.name.regex,
						guess.multifurcation.threshold = control$guess.multifurcation.threshold,
						max.reads.per.host = control$max.reads.per.host,
						multifurcation.threshold = control$multifurcation.threshold,
						norm.constants = control$norm.constants,
						norm.ref.file.name = control$norm.ref.file.name,
						norm.standardise.gag.pol = control$norm.standardise.gag.pol,
						no.progress.bars = control$no.progress.bars,
						outgroup.name = control$outgroup.name,
						parsimony.blacklist.k = control$sankoff.k,
						prune.blacklist = control$prune.blacklist,
						ratio.blacklist.threshold = control$ratio.blacklist.threshold, 
						raw.blacklist.threshold = control$raw.blacklist.threshold,			
						recombination.file.directory = control$recombination.file.directory,
						recombination.file.regex = control$recombination.file.regex,
						relaxed.ancestry = control$relaxed.ancestry,
						sankoff.k = control$sankoff.k,
						sankoff.unassigned.switch.threshold = control$sankoff.unassigned.switch.threshold,
						seed = control$seed,
						splits.rule = control$splits.rule,
						tip.regex = control$tip.regex,
						tree.file.regex = control$tree.file.regex,
						use.ff = control$use.ff,
						user.blacklist.directory = control$user.blacklist.directory, 
						user.blacklist.file.regex = control$user.blacklist.file.regex,
						verbosity = control$verbosity 
						)
		#	save output
		save(phsc, file=file.path(outdir,paste0('ptyr',df[i,RUN],'_analyse_trees.RData')))
		#	delete working directory
		sapply( list.files(tree.file.directory, full.names=TRUE), file.remove)
		file.remove(tree.file.directory)
	}		
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
	allow.mt					<- TRUE
	tip.regex					<- "^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$"
	tree.file.regex				<- "^ptyr192_InWindow_([0-9]+_to_[0-9]+)\\.tree$"
	verbosity					<- 1
	
	phsc	<- phyloscanner.analyse.trees(tree.file.directory,
			allow.mt=allow.mt,
			alignment.file.directory = NULL, 
			alignment.file.regex = NULL,
			blacklist.underrepresented = FALSE,
			count.reads.in.parsimony = TRUE,
			do.dual.blacklisting = FALSE,
			duplicate.file.directory = NULL,
			duplicate.file.regex = NULL,
			file.name.regex = file.name.regex,
			guess.multifurcation.threshold = FALSE,
			max.reads.per.host = max.reads.per.host,
			multifurcation.threshold = multifurcation.threshold,
			norm.constants = NULL,
			norm.ref.file.name = norm.ref.file.name,
			norm.standardise.gag.pol = TRUE,
			no.progress.bars = FALSE,
			outgroup.name = outgroup.name,
			parsimony.blacklist.k = sankoff.k,
			prune.blacklist = FALSE,
			ratio.blacklist.threshold = 0, 
			raw.blacklist.threshold = raw.blacklist.threshold,			
			recombination.file.directory = NULL,
			recombination.file.regex = NULL,
			relaxed.ancestry = relaxed.ancestry,
			sankoff.k = sankoff.k,
			sankoff.unassigned.switch.threshold = sankoff.unassigned.switch.threshold,
			seed = seed,
			splits.rule = splits.rule,
			tip.regex = tip.regex,
			tree.file.regex = tree.file.regex,
			use.ff = FALSE,
			user.blacklist.directory = NULL, 
			user.blacklist.file.regex = NULL,
			verbosity = verbosity 
			)
	
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

