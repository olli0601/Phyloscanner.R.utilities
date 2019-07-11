pty.MRC.stage1.generate.read.alignments<- function()
{
	#	set up working environment
	require(Phyloscanner.R.utilities)
	#HOME			<<- '/rds/general/project/ratmann_pangea_analyses_mrc_uvri/live'
	#data.dir		<- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC'
	#prog.pty		<- '/rds/general/user/or105/home/phyloscanner/phyloscanner_make_trees.py'
	HOME			<<- '~/sandbox/DeepSeqProjects'
	data.dir		<- '~/sandbox/DeepSeqProjects/MRCPopSample_data'
	prog.pty		<- '/Users/Oliver/git/phylotypes/phyloscanner_make_trees.py'
	
	in.dir			<- file.path(HOME,'MRCPopSample_phsc_stage1_input')		
	work.dir		<- file.path(HOME,"MRCPopSample_phsc_work")
	out.dir			<- file.path(HOME,"MRCPopSample_phsc_stage1_output")	
	dir.create(in.dir)
	dir.create(work.dir)
	dir.create(out.dir)
	
	#
	#	define phyloscanner runs for all pairs of batches of individuals of the population-based sample
	#
	
	if(0)
	{
		infile	<- 'phsc_runs_MRC_data.csv'	#	file listing all available bam files and corresponding individual IDs
		outfile	<- 'phsc_runs_MRC_stage1_n2531_181026.csv'			
		dbam	<- as.data.table(read.csv(file.path(in.dir,infile)))
		dbam	<- subset(dbam, grepl('remap\\.bam',BAM))
		setnames(dbam,'INDIVIDUAL','IND')
		#	make data.table with just individuals
		dind	<- unique(subset(dbam, select=c(IND)))
		#	create batches and create runs of pairs of batches
		pty.runs<- phsc.define.stage1.analyses(dind$IND, batch.size=50)
		#	add BAM files per individual
		pty.runs<- merge(pty.runs, subset(dbam, select=c(IND, BAM)), by='IND', allow.cartesian=TRUE)
		setkey(pty.runs, PTY_RUN)
		#	save to file
		write.csv(pty.runs, file=file.path(in.dir,outfile), row.names=FALSE)		
	}		
	if(1)
	{
		pty.runs<- as.data.table(read.csv(file.path(in.dir, "phsc_runs_MRC_stage1_n2531_181026.csv")))
	}
	
	#
	#	define phyloscanner input args to generate read alignments
	#
	
	pty.select			<- NA #1:2
	pty.args			<- list(	prog.pty=prog.pty, 
									prog.mafft='mafft',  
									data.dir=data.dir, 
									work.dir=work.dir, 
									out.dir=out.dir, 
									alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX_v2.fasta"),
									alignments.root='REF_CPX_AF460972', 
									alignments.pairwise.to='REF_B_K03455',
									window.automatic='', 
									merge.threshold=2, 
									min.read.count=1, 
									quality.trim.ends=23, 
									min.internal.quality=23, 
									merge.paired.reads=TRUE, 
									no.trees=TRUE, 
									dont.check.duplicates=FALSE,
									dont.check.recombination=TRUE,
									win=c(800,9400,25,250),				 				
									keep.overhangs=FALSE,
									mem.save=0,
									verbose=TRUE,					
									select=pty.select	#of 240
									)	
	save(pty.args, file=file.path(in.dir, 'phsc_args_stage1_create_read_alignments.rda'))
	
	
	
	#
	#	define bash scripts to generate read alignments for all pairs of batches
	#
	
	#	check which (if any) batches have already been processed, and remove from TODO list
	tmp			<- data.table(FILE_FASTA=list.files(out.dir, pattern='^ptyr[0-9]+_trees$', full.names=TRUE))
	tmp[, PTY_RUN:= as.integer(gsub('ptyr([0-9]+)_.*','\\1',basename(FILE_FASTA)))]
	tmp			<- merge(tmp, tmp[, list(FILE_FASTA_N= length(list.files(FILE_FASTA, pattern='fasta$'))), by='PTY_RUN'], by='PTY_RUN')
	pty.runs	<- merge(pty.runs, tmp, by='PTY_RUN', all.x=1)
	pty.runs	<- subset(pty.runs, FILE_FASTA_N==0)
	
	#	search for bam files and references and merge with runs	
	pty.runs	<- subset(pty.runs, select=c(PTY_RUN, IND))	
	tmp			<- phsc.find.bam.and.references(pty.args[['data.dir']], regex.person='^([A-Z0-9]+-[A-Z0-9]+)-.*$')	
	pty.runs	<- merge(pty.runs, tmp, by='IND')
	#	create UNIX bash scripts
	pty.runs	<- unique(pty.runs, by=c('PTY_RUN','UNIT_ID','BAM'))
	pty.runs	<- pty.runs[, list(BAM=BAM, REF=REF, SAMPLE=SAMPLE, RENAME_ID=paste0(IND,'-fq',seq_len(length(BAM)))), by=c('PTY_RUN','IND')]
	setkey(pty.runs, PTY_RUN)		
	setnames(pty.runs, c('IND','SAMPLE'), c('UNIT_ID','SAMPLE_ID'))
	pty.c		<- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)
	pty.c[, CASE_ID:= seq_len(nrow(pty.c))]
	
	#
	#	submit jobs to PBS job scheduler
	#
	
	#	define PBS variables
	hpc.load			<- "module load intel-suite/2015.1 mpi raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"	# make third party requirements available	 
	hpc.select			<- 1						# number of nodes
	hpc.nproc			<- 1						# number of processors on node
	hpc.walltime		<- 171						# walltime
	hpc.q				<- "pqeelab"				# PBS queue
	hpc.mem				<- "12gb" 					# RAM	
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
	cmd		<- pty.c[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
	cmd		<- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]			
	cmd		<- paste(pbshead,cmd,sep='\n')	
	#	submit job
	outfile		<- gsub(':','',paste("readali",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
	outfile		<- file.path(pty.args[['work.dir']], outfile)
	cat(cmd, file=outfile)
	cmd 		<- paste("qsub", outfile)
	cat(cmd)
	cat(system(cmd, intern= TRUE))		
}

pty.MRC.stage1.zip.trees<- function()
{
	require(data.table)
	require(Phyloscanner.R.utilities)
	
	#	set up working environment	
	HOME			<<- '/rds/general/project/ratmann_pangea_analyses_mrc_uvri/live'
	data.dir		<- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC'	
	#in.dir			<- file.path(HOME,'MRCPopSample_phsc_stage1_output')		
	#work.dir		<- file.path(HOME,"MRCPopSample_phsc_work")
	#out.dir			<- file.path(HOME,"MRCPopSample_phsc_stage1_output")	
	
	#
	#	zip trees	
	if(1)
	{
		indirs	<- file.path(HOME,'MRCPopSample_phsc_stage1_output')
		#
		indirs	<- list.files(indirs, pattern='^ptyr[0-9]+_trees$', full.names=TRUE)
		allwin	<- data.table(W_FROM=seq(800,9150,25))
		#allwin	<- data.table(W_FROM=seq(800,9050,125))
		#indirs	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo/ptyr97_trees'
		for(i in seq_along(indirs))
		{
			indir	<- indirs[i]
			pty.run	<- as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(indir)))
			#	check if we have all fasta and tree files
			infiles	<- data.table(F=list.files(indir,pattern='ptyr.*fasta$',full.names=TRUE))
			infiles[, W_FROM:= as.integer(gsub('.*_InWindow_([0-9]+)_.*','\\1',basename(F)))]
			infiles	<- merge(allwin, infiles, by='W_FROM', all.x=1)			 					
			missfs	<- subset(infiles, is.na(F))[, W_FROM]
			if(length(missfs))
				cat('\nIn',indir,'Found missing fasta files for',paste(missfs,collapse=', '))
			infiles	<- data.table(F=list.files(indir,pattern='ptyr.*tree$',full.names=TRUE))
			infiles[, W_FROM:= as.integer(gsub('.*_InWindow_([0-9]+)_.*','\\1',basename(F)))]
			infiles	<- merge(allwin, infiles, by='W_FROM', all.x=1)
			misstrs	<- subset(infiles, is.na(F))[, W_FROM]
			if(length(misstrs))
				cat('\nIn',indir,'Found missing tree files for',paste(misstrs,collapse=', '))
			zipit	<- 0
			if(!length(missfs) & !length(misstrs))
			{
				cat('\nIn',indir,'Found all fasta and tree files')
				zipit	<- 1
			}			
			#if(!length(setdiff(misstrs,missfs)))
			#{ 
			#	cat('\nIn',indir,'Found all tree files for which there is a fasta file')
			#	zipit	<- 1
			#}	
			#zipit	<- 0
			#
			if(zipit)
			{			
				cat('\nProcess',indir)
				#	first combine all zip files into ptyrXXX_otherstuff.zip
				infiles	<- data.table(F=list.files(indir,pattern='zip$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				stopifnot(!nrow(subset(infiles, is.na(PTY_RUN))))		
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_otherstuff.zip'))]
				cat('\nZip to file', tmp,'...\n')
				suppressWarnings(invisible( infiles[, list(RTN= unzip(F, overwrite=FALSE, exdir=file.path(indir,'tmp42'))), by='F'] ))
				invisible( infiles[, list(RTN= file.remove(F)), by='F'] )
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	now zip fasta files
				infiles	<- data.table(F=list.files(indir,pattern='ptyr.*fasta$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				invisible( infiles[, file.rename(F, file.path(indir,'tmp42',basename(F))), by='F'] )
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_trees_fasta.zip'))]
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	now zip tree files
				infiles	<- data.table(F=list.files(indir,pattern='ptyr.*tree$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				invisible( infiles[, file.rename(F, file.path(indir,'tmp42',basename(F))), by='F'] )
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_trees_newick.zip'))]		
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	remove tmp dir
				invisible( file.remove( file.path(indir,'tmp42') ) )
				#	move one level down
				infiles	<- data.table(F=list.files(indir, full.names=TRUE))
				invisible( infiles[, file.rename(F, file.path(dirname(indir),basename(F))), by='F'] )
				cat('\nDone',indir)
			}
			#if(!length(misstrs))
			if(zipit)
				invisible(unlink(indir, recursive=TRUE))
			#	expand again if asked to
			#if(length(misstrs))
			if(0)
			{
				cat('\nExtract',file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_fasta.zip')))
				unzip(file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_fasta.zip')), junkpaths=TRUE, exdir=indir)
				cat('\nExtract',file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_newick.zip')))
				unzip(file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_newick.zip')), junkpaths=TRUE, exdir=indir)
			}
		}					
	}
}

pty.MRC.stage1.generate.trees<- function()
{
	require(data.table)
	require(Phyloscanner.R.utilities)
	
	#	set up working environment	
	HOME			<<- '/rds/general/project/ratmann_pangea_analyses_mrc_uvri/live'
	data.dir		<- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC'
	prog.pty		<- '/rds/general/user/or105/home/phyloscanner/phyloscanner_make_trees.py'
	#HOME			<<- '~/sandbox/DeepSeqProjects'
	#data.dir		<- '~/sandbox/DeepSeqProjects/MRCPopSample_data'
	#prog.pty		<- '/Users/Oliver/git/phylotypes/phyloscanner_make_trees.py'	
	in.dir			<- file.path(HOME,'MRCPopSample_phsc_stage1_output')		
	work.dir		<- file.path(HOME,"MRCPopSample_phsc_work")
	out.dir			<- file.path(HOME,"MRCPopSample_phsc_stage1_output")	
		
	
	#
	#	generate trees	
	infiles	<- data.table(FI=list.files(in.dir, pattern='fasta$', full.names=TRUE, recursive=TRUE))
	infiles[, FO:= gsub('fasta$','tree',FI)]
	infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FI)))]
	infiles[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FI)))]
	#	check which (if any) trees have already been processed, and remove from TODO list
	tmp		<- data.table(FT=list.files(out.dir, pattern='^ptyr.*tree$', full.names=TRUE, recursive=TRUE))
	tmp[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FT)))]
	tmp[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FT)))]
	infiles	<- merge(infiles, tmp, by=c('PTY_RUN','W_FROM'), all.x=1)
	infiles	<- subset(infiles, is.na(FT))	
	setkey(infiles, PTY_RUN, W_FROM)			
	#	determine number of taxa per fasta file, and sort alignments by size
	#	so that trees with smallest size are generated first
	#tmp		<- infiles[, list(N_TAXA=as.integer(system(paste0('grep -c "^>" ', FI), intern=TRUE))), by=c('FI')]
	#infiles	<- merge(infiles, tmp, by='FI')
	#infiles	<- infiles[order(N_TAXA),]	
	
	#
	#	define pbs header, max 10,000 trees (re-start a few times)
	#	create header first because RAXML call depends on single-core / multi-core
	#	TODO: 
	#		30 runs are roughly 10,000 trees
	#		every time you set up a new job array, select the next 30 runs.   
	#		For example, the job array will be 151:180
	#infiles	<- subset(infiles, PTY_RUN%in%c(1101:1400))
	#infiles	<- subset(infiles, PTY_RUN%in%c(1001:1400))
	infiles[, CASE_ID2:= seq_len(nrow(infiles))]
	infiles[, CASE_ID:= ceiling(CASE_ID2/1)]
	print(infiles)
	#infiles[, CASE_ID:= ceiling(CASE_ID2/1)]
	stop()
	hpc.load			<- "module load intel-suite/2015.1 mpi raxml/8.2.9"	# make third party requirements available	 
	hpc.select			<- 1						# number of nodes
	hpc.nproc			<- 1						# number of processors on node
	hpc.walltime		<- 123						# walltime
	#	TODO:
	#		run either this block to submit a job array to college machines 
	#		the choice depends on whether the previous job array on college machines is done, 
	#		or on whether the previous job array on Oliver's machines is done.
	if(0)		
	{
		hpc.q			<- NA						# PBS queue
		hpc.mem			<- "2gb" 					# RAM
		raxml.pr		<- ifelse(hpc.nproc==1, 'raxmlHPC-SSE3', 'raxmlHPC-PTHREADS-SSE3')	#on older machines without AVX instructions
	}
	#		or run this block to submit a job array to Oliver's machines
	if(1)
	{
		hpc.q			<- "pqeelab"				# PBS queue
		hpc.mem			<- "6gb" 					# RAM
		raxml.pr		<- ifelse(hpc.nproc==1, 'raxmlHPC-AVX','raxmlHPC-PTHREADS-AVX')		#on newer machines with AVX instructions
	}
	#
	#	
	hpc.array			<- infiles[, max(CASE_ID)]	# number of runs for job array	
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
	
	#
	#	create UNIX bash script to generate trees with RAxML
	#	
	raxml.args	<- ifelse(hpc.nproc==1, '-m GTRCAT --HKY85 -p 42 -o REF_B_K03455', paste0('-m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42 -o REF_B_K03455'))
	pty.c	<- infiles[, list(CMD=raxml.cmd(FI, outfile=FO, pr=raxml.pr, pr.args=raxml.args)), by=c('CASE_ID','CASE_ID2')]
	pty.c	<- pty.c[, list(CMD=paste(CMD, collapse='\n')), by='CASE_ID']
	pty.c[1,cat(CMD)]
	
	#
	#	make array job
	cmd		<- pty.c[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
	cmd		<- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]			
	cmd		<- paste(pbshead,cmd,sep='\n')	
	#	submit job
	outfile		<- gsub(':','',paste("trs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
	outfile		<- file.path(work.dir, outfile)
	cat(cmd, file=outfile)
	cmd 		<- paste("qsub", outfile)
	cat(cmd)
	cat(system(cmd, intern= TRUE))	
}


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
	require(Phyloscanner.R.utilities)
		
	#	setting up workspace
	HOME			<- '~/sandbox/DeepSeqProjects'
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


pty.4.highlysupported.links<- function()
{
	require(Phyloscanner.R.utilities)
	
	#	setting up workspace
	HOME			<- '~/sandbox/DeepSeqProjects'
	indir			<- file.path(HOME, 'RakaiPopSample_phyloscanner_analysis')
	
	#	load pairs and networks
	infile.pairs	<- file.path(indir,'phsc_analysis_of_dataset_S1_allpairs.rda')
	infile.networks	<- file.path(indir,'phsc_analysis_of_dataset_S1_allnetworks.rda')
	load( infile.pairs )	# loads rtp, rplkl, rpw
	load( infile.networks )	# loads rtn, rtnn
		
	#	classify linkages
	conf.cut		<- 0.6				
	rtnn[, SELECT:= NA_character_]
	set(rtnn, rtnn[, which(is.na(PTY_RUN))], 'SELECT', 'insufficient deep sequence data for at least one individual')
	set(rtnn, rtnn[, which(!is.na(PTY_RUN) & is.na(LINK_12) & is.na(LINK_21))], 'SELECT', 'ph unlinked pair')
	set(rtnn, rtnn[, which(!is.na(POSTERIOR_SCORE_LINKED) & POSTERIOR_SCORE_LINKED<=conf.cut)], 'SELECT', 'unclear if pair ph linked or unlinked')
	set(rtnn, rtnn[, which(!is.na(POSTERIOR_SCORE_LINKED) & POSTERIOR_SCORE_LINKED>conf.cut)], 'SELECT', 'ph linked pair direction not resolved')
	set(rtnn, rtnn[, which(!is.na(POSTERIOR_SCORE_LINKED) & POSTERIOR_SCORE_LINKED>conf.cut & POSTERIOR_SCORE_12>conf.cut)], 'SELECT', 'ph linked pair direction 12')
	set(rtnn, rtnn[, which(!is.na(POSTERIOR_SCORE_LINKED) & POSTERIOR_SCORE_LINKED>conf.cut & POSTERIOR_SCORE_21>conf.cut)], 'SELECT', 'ph linked pair direction 21')	
	
	#	add gender
	tmp		<- subset(rtp, select=c(ID1, ID2, ID1_SEX, ID2_SEX))
	rtnn	<- merge(rtnn, tmp, by=c('ID1','ID2'))
	#	re-order so male is ID1
	tmp		<- subset(rtnn, ID1_SEX=='F' & ID2_SEX=='M')
	setnames(tmp, colnames(tmp), gsub('xx','ID2',gsub('ID2','ID1',gsub('ID1','xx',gsub('xx','12',gsub('12','21',gsub('21','xx',colnames(tmp))))))))
	set(tmp, NULL, 'SELECT', tmp[, gsub('xx','12',gsub('12','21',gsub('21','xx',SELECT)))])
	rtnn	<- rbind(subset(rtnn, !(ID1_SEX=='F' & ID2_SEX=='M')), tmp)
	#	define gender pair
	rtnn[, PAIR_SEX:= paste0(ID1_SEX,ID2_SEX)]
	
	#	select pairs with highly supported linkages
	rtp		<- subset(rtnn, !grepl('unlinked|insufficient',SELECT))
	rtp[, table(PAIR_SEX)]
	#FF  MF  MM 
	#80 376  81
	#	conclusion: high proportion of missed intermediates
	
	
	#	select male-female pairs with highly supported direction
	rtpd	<- subset(rtnn, ID1_SEX=='M' & ID2_SEX=='F' & grepl('direction 12|direction 21',SELECT))
	set(rtpd, NULL, 'SELECT', rtpd[, gsub('ph linked pair direction 21','fm',gsub('ph linked pair direction 12','mf',SELECT))])
	setnames(rtpd, 'SELECT', 'PHYSCANNER_DIR')
	rtpd	<- subset(rtpd, select=c(ID1, ID2, PHYSCANNER_DIR))
	#	293/376

	#	read pairs with clinical data on direction (assuming these are linked)
	infile		<- '~/sandbox/DeepSeqProjects/RakaiPopSample_data/Dataset_S3.csv'
	red			<- as.data.table(read.csv(infile))
	setnames(red, c('MALE_ID','FEMALE_ID'), c('ID1','ID2'))
	rtpd		<- merge(rtpd, red, by=c('ID1','ID2'))
	rtpd[, PHYSCANNER_DIR_CONSISTENT:= as.integer(PHYSCANNER_DIR==EPID_EVIDENCE_DIR)]
	rtpd[, table(PHYSCANNER_DIR_CONSISTENT)]
}
