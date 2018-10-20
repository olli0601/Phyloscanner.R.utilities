pty.2.infer.phylo.relationships.on.stage2.trees<- function() 
{
	require(Phyloscanner.R.utilities)
	
	#
	#	INPUT ARGS PLATFORM
	if(1)
	{	
		HOME				<<- '~/sandbox/DeepSeqProjects'								
		#	unzip files in Data Set S1 to directory with full path name as in 'in.dir'
		#	Data Set S1 contains batches of deep sequence trees that are to be processed with phyloscanner
		in.dir				<- file.path(HOME,'RakaiPopSample_deepseqtrees')	
		#	define output directory
		out.dir				<- file.path(HOME,"RakaiPopSample_phyloscanner_out")
		#	define directory for temporary files 
		work.dir			<- file.path(HOME,"RakaiPopSample_phyloscanner_work")		
		#	define path to 'phyloscanner_make_trees.py'		
		prog.pty			<- '/Users/Oliver/git/phylotypes/phyloscanner_make_trees.py'
	}	
	
	# create directories if they dont exist
	dir.create(out.dir, showWarnings=FALSE)
	dir.create(work.dir, showWarnings=FALSE)
	
	#
	#	phyloscanner default input arguments to process deep-sequence trees
	pty.args			<- list(	prog.pty=prog.pty, 
			prog.mafft=NA, 
			prog.raxml=NA, 
			data.dir=NA, 
			work.dir=work.dir, 
			out.dir=out.dir, 
			alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX_v2.fasta"),
			alignments.root='REF_CPX_AF460972', 
			alignments.pairwise.to='REF_B_K03455',
			bl.normalising.reference.file=system.file(package="Phyloscanner.R.utilities", "data", "hiv.hxb2.norm.constants.rda"),
			bl.normalising.reference.var='MEDIAN_PWD',														
			window.automatic= '', 
			merge.threshold=0, 
			min.read.count=1, 
			quality.trim.ends=23, 
			min.internal.quality=23, 
			merge.paired.reads=TRUE, 
			no.trees=FALSE, 
			dont.check.duplicates=FALSE,
			dont.check.recombination=TRUE,
			num.bootstraps=1,
			all.bootstrap.trees=TRUE,
			strip.max.len=350, 
			min.ureads.individual=NA, 
			win=c(800,9400,25,250), 				
			keep.overhangs=FALSE,
			use.blacklisters=c('ParsimonyBasedBlacklister','DownsampleReads'),
			tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$',
			roguesubtree.kParam=20,
			roguesubtree.prop.threshold=0,
			roguesubtree.read.threshold=20,
			dwns.maxReadsPerPatient=50,	
			multifurcation.threshold=1e-5,
			split.rule='s',
			split.kParam=20,
			split.proximityThreshold=0,
			split.readCountsMatterOnZeroBranches=TRUE,
			split.pruneBlacklist=FALSE,
			trms.allowMultiTrans=TRUE,
			pw.trmw.min.reads=30,									
			pw.trmw.min.tips=1,
			pw.trmw.close.brl=0.025,
			pw.trmw.distant.brl=0.05,
			pw.prior.keff=2,
			pw.prior.neff=3,
			pw.prior.keff.dir=2,
			pw.prior.neff.dir=3,				
			pw.prior.calibrated.prob=0.66,
			mem.save=0,
			verbose=TRUE,				
			select=NA 
	)	
	save(pty.args, file=file.path(out.dir, 'pty.args.rda'))
		
	
	#
	#	phyloscanner runs on batches of deep-sequence trees	
	if(1)
	{
		#	for each batch of trees to be processed:
		#	find list of patients in input directory
		pty.c	<- data.table(FILE_BAM=list.files(in.dir, pattern='_patients.txt', full.names=TRUE))
		pty.c[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_patients.txt','',basename(FILE_BAM))))]
		#	check which (if any) batches have already been processed, and remove from TODO list
		tmp		<- data.table(FILE_TRMW=list.files(out.dir, pattern='_pairwise_relationships.rda', full.names=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_pairwise_relationships.rda','',basename(FILE_TRMW))))]
		pty.c	<- merge(pty.c, tmp, by='PTY_RUN', all.x=1)
		pty.c	<- subset(pty.c, is.na(FILE_TRMW))
		#	for each batch of deep-sequence trees: create bash script to run phyloscanner 
		setkey(pty.c, PTY_RUN)		
		pty.c	<- pty.c[, { 
					#FILE_BAM<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr1_bam.txt'
					#FILE_BAM<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr1_bam.txt'
					#cat('\n',FILE_BAM)
					prefix.infiles	<- gsub('patients.txt','',FILE_BAM)
					cmd				<- phsc.cmd.phyloscanner.one.resume(prefix.infiles, pty.args)
					list(CMD=cmd) 
				}, by='PTY_RUN']		
		pty.c[1,cat(CMD)]		
	}
	
	#
	# run phyloscanner on each batch
	# option1: process each batch on a standard UNIX desktop
	if(1)
	{
		#	write bash scripts to file
		invisible(pty.c[,	{					
							outfile		<- gsub(':','',paste("phsc",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
							outfile		<- file.path(pty.args[['work.dir']], outfile)
							cat(CMD, file=outfile)
							Sys.chmod(outfile, mode="777")
							Sys.sleep(1)
						}, by='PTY_RUN'])
		#	now run each bash script manually on the command prompt
	}
	
	#
	# run phyloscanner on each batch
	# option2: process each batch on a high performance environment with a job scheduling system	
	# write bash scripts to file and submit to job scheduling system on a high performance environment 
	if(0)
	{
		# define variables for PBS header. this will depend on your job scheduler.
		hpc.load			<- "module load R/3.3.3"	# make R available 
		hpc.select			<- 1						# number of nodes
		hpc.nproc			<- 1						# number of processors on node
		hpc.walltime		<- 15						# walltime
		hpc.q				<- "pqeelab"				# PBS queue
		hpc.mem				<- "6gb" 					# RAM
		
		# define PBS header
		pbshead	<- "#!/bin/sh"
		tmp		<- paste("#PBS -l walltime=", hpc.walltime, ":59:59,pcput=", hpc.walltime, ":45:00", sep = "")
		pbshead	<- paste(pbshead, tmp, sep = "\n")
		tmp		<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
		pbshead 	<- paste(pbshead, tmp, sep = "\n")
		pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")
		if (!is.na(hpc.q)) 
			pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
		pbshead <- paste(pbshead, hpc.load, sep = "\n")
		# for each batch of deep-sequence trees: submit PBS job 
		invisible(pty.c[,	{					
						cmd			<- paste(pbshead,'cd $TMPDIR',sep='\n')
						cmd			<- paste(cmd,CMD,sep='\n')	
						outfile		<- gsub(':','',paste("phsc",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
						outfile		<- file.path(pty.args[['work.dir']], outfile)
						cat(CMD, file=outfile)
						cmd 		<- paste("qsub", outfile)
						cat(cmd)
						cat(system(cmd, intern= TRUE))
						Sys.sleep(1)						
					}, by='PTY_RUN'])
	}
}


pty.3.infer.transmission.networks.from.phylo.relationships<- function()
{
	require(data.table)	
	require(igraph)
	require(sna)
	require(Phyloscanner.R.utilities)
		
	#	setting up workspace
	HOME			<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'
	indir			<- file.path(HOME, 'RakaiPopSample_phyloscanner_analysis')
	outdir			<- indir
	outfile.base	<- file.path(outdir,'phsc_analysis_of_dataset_S1')
	neff.cut		<- 3
	conf.cut		<- 0.6		
	
	#	get meta-data for individuals in general pop cohort
	infile			<- "~/Dropbox (SPH Imperial College)/2017_phyloscanner_validation/Supp_Data/Data_Set_S2.csv"
	dmeta			<- as.data.table(read.csv(infile, stringsAsFactors=FALSE))
	
	
	#	get pairs between whom linkage is not excluded
	tmp				<- phsc.find.linked.pairs(indir, batch.regex='^ptyr([0-9]+)_.*', conf.cut=conf.cut, neff.cut=neff.cut, verbose=TRUE, dmeta=dmeta)
	rtp				<- copy(tmp$linked.pairs)
	rplkl			<- copy(tmp$relationship.counts)
	rpw				<- copy(tmp$windows)
	save(rtp, rplkl, rpw, file=paste0(outfile.base,'_allpairs.rda'))
	
	#	plot phyloscan
	rpw2 <- subset(rpw, ID1=='RkA04565F' & ID2=='RkA05315F')		
	p	<- phsc.plot.phyloscan(rpw2)
	p
	ggsave(file=paste0(outfile.base,'_phyloscan_RkA04565F_RkA05315F.png'), width=6, height=2.8, units='in', dpi=400)
	
	#	reconstruct transmission networks	
	tmp				<- phsc.find.transmission.networks.from.linked.pairs(rtp, rplkl, conf.cut=conf.cut, neff.cut=neff.cut, verbose=TRUE)
	rtn				<- copy(tmp$transmission.networks)
	rtnn			<- copy(tmp$most.likely.transmission.chains)		
	save(rtn, rtnn, file=paste0(outfile.base,'_allnetworks.rda'))
	
	#	plot transmission networks: one example	
	idclus	<- sort(unique(rtn$IDCLU))
	di		<- copy(dmeta)									
	df		<- subset(rtn, IDCLU==idclus[34])
	set(df, NULL, c('ID1_SEX','ID2_SEX'), NULL)
	p		<- phsc.plot.transmission.network(df, di, point.size=10, 
			edge.gap=0.04, 
			edge.size=0.4, 
			curvature= -0.2, 
			arrow=arrow(length=unit(0.04, "npc"), type="open"), 
			curv.shift=0.06, 
			label.size=3, 
			node.label='ID', 			 
			node.fill='SEX', 
			node.shape.values=c('NA'=16), 
			node.fill.values=c('F'='hotpink2', 'M'='steelblue2'),
			threshold.linked=0.6)	
	png(file=paste0(outfile.base,'_trmnetwork_34.png'), width=6, height=6, units='in', res=400)		
	print(p)
	dev.off()
	
	#	plot the corresponding most likely transmission chain
	layout	<- p$layout 
	di		<- copy(dmeta)									
	df		<- subset(rtnn, IDCLU==idclus[34])	
	p2		<- phsc.plot.most.likely.transmission.chain(df, di, point.size=10, 
			edge.gap=0.04, 
			edge.size=0.4, 
			curvature= -0.2, 
			arrow=arrow(length=unit(0.04, "npc"), type="open"), 
			curv.shift=0.06, 
			label.size=3, 
			node.label='ID', 			 
			node.fill='SEX', 
			node.shape.values=c('NA'=16), 
			node.fill.values=c('F'='hotpink2', 'M'='steelblue2'),
			threshold.linked=0.6,
			layout=layout)	
	png(file=paste0(outfile.base,'_trmchain_34.png'), width=6, height=6, units='in', res=400)		
	print(p2)
	dev.off()
}


pty.4.female.female.links<- function()
{
	require(data.table)	
	require(igraph)
	require(sna)
	require(Phyloscanner.R.utilities)
	
	#	setting up workspace
	HOME			<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'
	indir			<- file.path(HOME, 'RakaiPopSample_phyloscanner_analysis')
	outdir			<- indir
	outfile.base	<- file.path(outdir,'phsc_analysis_of_dataset_S1')
	neff.cut		<- 3
	conf.cut		<- 0.6		
	
	#	get meta-data for individuals in general pop cohort
	infile			<- "~/Dropbox (SPH Imperial College)/2017_phyloscanner_validation/Supp_Data/Data_Set_S2.csv"
	dmeta			<- as.data.table(read.csv(infile, stringsAsFactors=FALSE))
	
	
	#	get pairs between whom linkage is not excluded
	tmp				<- phsc.find.linked.pairs(indir, batch.regex='^ptyr([0-9]+)_.*', conf.cut=conf.cut, neff.cut=neff.cut, verbose=TRUE, dmeta=dmeta)
	rtp				<- copy(tmp$linked.pairs)
	rplkl			<- copy(tmp$relationship.counts)
	rpw				<- copy(tmp$windows)
	save(rtp, rplkl, rpw, file=paste0(outfile.base,'_allpairs.rda'))
	
	#	plot phyloscan
	rpw2 <- subset(rpw, ID1=='RkA04565F' & ID2=='RkA05315F')		
	p	<- phsc.plot.phyloscan(rpw2)
	p
	ggsave(file=paste0(outfile.base,'_phyloscan_RkA04565F_RkA05315F.png'), width=6, height=2.8, units='in', dpi=400)
	
	#	reconstruct transmission networks	
	tmp				<- phsc.find.transmission.networks.from.linked.pairs(rtp, rplkl, conf.cut=conf.cut, neff.cut=neff.cut, verbose=TRUE)
	rtn				<- copy(tmp$transmission.networks)
	rtnn			<- copy(tmp$most.likely.transmission.chains)		
	save(rtn, rtnn, file=paste0(outfile.base,'_allnetworks.rda'))
	
	#	plot transmission networks: one example	
	idclus	<- sort(unique(rtn$IDCLU))
	di		<- copy(dmeta)									
	df		<- subset(rtn, IDCLU==idclus[34])
	set(df, NULL, c('ID1_SEX','ID2_SEX'), NULL)
	p		<- phsc.plot.transmission.network(df, di, point.size=10, 
			edge.gap=0.04, 
			edge.size=0.4, 
			curvature= -0.2, 
			arrow=arrow(length=unit(0.04, "npc"), type="open"), 
			curv.shift=0.06, 
			label.size=3, 
			node.label='ID', 			 
			node.fill='SEX', 
			node.shape.values=c('NA'=16), 
			node.fill.values=c('F'='hotpink2', 'M'='steelblue2'),
			threshold.linked=0.6)	
	png(file=paste0(outfile.base,'_trmnetwork_34.png'), width=6, height=6, units='in', res=400)		
	print(p)
	dev.off()
	
	#	plot the corresponding most likely transmission chain
	layout	<- p$layout 
	di		<- copy(dmeta)									
	df		<- subset(rtnn, IDCLU==idclus[34])	
	p2		<- phsc.plot.most.likely.transmission.chain(df, di, point.size=10, 
			edge.gap=0.04, 
			edge.size=0.4, 
			curvature= -0.2, 
			arrow=arrow(length=unit(0.04, "npc"), type="open"), 
			curv.shift=0.06, 
			label.size=3, 
			node.label='ID', 			 
			node.fill='SEX', 
			node.shape.values=c('NA'=16), 
			node.fill.values=c('F'='hotpink2', 'M'='steelblue2'),
			threshold.linked=0.6,
			layout=layout)	
	png(file=paste0(outfile.base,'_trmchain_34.png'), width=6, height=6, units='in', res=400)		
	print(p2)
	dev.off()
}
