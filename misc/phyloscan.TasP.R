TasP.calculate.bam.readlength.coverage<- function()
{
	require(ggplot2)
	require(data.table)
	require(Rsamtools)		
	
	pty.data.dir	<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/data'			
	outfile			<- file.path(pty.data.dir,'bam_stats_170928.rda')			
	
	bfiles			<- data.table(FILE=list.files(pty.data.dir, pattern='bam$'))	
	bfiles			<- subset(bfiles, !grepl('Contam',FILE) & !grepl('merged',FILE))
	bfiles[, FILE_ID:=gsub('.bam','',FILE)]	
	#FILE<- "14939_1_80.bam"
	#cat('\nreading',FILE)
	#bam<- scanBam(file.path(pty.data.dir,FILE))[[1]]	#scan entire file
	bam.cov	<- bfiles[,{
				cat('\nCOV reading',FILE)
				z	<- scanBam(file.path(pty.data.dir,FILE), param=ScanBamParam(what=c('qwidth','pos','rname')))[[1]]
				tmp	<- IRanges(start=z$pos, width=z$qwidth)
				tmp <- coverage(tmp)
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths, REF=z$rname[1])
			}, by='FILE_ID']
	#	get lengths of all reads in quality trimmed bam file
	#z	<- scanBam(file.path(pty.data.dir,FILE), param=ScanBamParam(what=c('qwidth','qual')))
	bam.len			<- bfiles[,{
				#FILE<- "14939_1_80.bam"
				cat('\nLEN reading',FILE)
				z	<- scanBam(file.path(pty.data.dir,FILE), param=ScanBamParam(what=c('qwidth')))
				list(QU=seq(0,320,20), CDF=ecdf(z[[1]][['qwidth']])(seq(0,320,20)))
				#list(PR=seq(0.01,0.99,0.01), QU=quantile(z[[1]][['qwidth']], p=seq(0.01,0.99,0.01)))				
			}, by='FILE']
	#
	save(bam.len, bam.cov, file= outfile)
}

TasP.evaluate.window.len<- function()
{	
	tip.regex	<- '(.*)_read_([0-9]+)_count_([0-9]+)'
	min.count	<- 30
	indir		<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch'
	plot.file	<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/171001_evaluate_winlen.pdf'
	infiles		<- data.table(F=list.files(indir, recursive=TRUE, full.names=TRUE, pattern='subgraphs_s_ptyr1_InWindow_.*.rda$'))
	infiles[, WLEN:= as.integer(gsub('.*_winlen([0-9]+).*','\\1',F))]
	infiles[, WFROM:= as.integer(gsub('.*_InWindow_([0-9]+).*','\\1',F))]
	infiles[, PTY_RUN:= as.integer(gsub('.*ptyr([0-9]+)_.*','\\1',F))]
	
	
	dw			<- infiles[, {
				#F			<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/ptyr1_winlen180/subgraphs_s_ptyr1_InWindow_2300_to_2479.rda'
				load(F)
				tmp			<- data.table(TAXA=tree$tip.label)
				tmp			<- subset(tmp, grepl(tip.regex, TAXA))
				tmp[, ID:= gsub(tip.regex,'\\1',TAXA)]
				tmp[, READ:= as.integer(gsub(tip.regex,'\\2',TAXA))]
				tmp[, COUNT:= as.integer(gsub(tip.regex,'\\3',TAXA))]				
				tmp			<- tmp[, list(READ=length(unique(READ)), COUNT=sum(COUNT)),by='ID']
				tmp
			}, by=c('PTY_RUN','WLEN','WFROM')]
	
	tmp			<- subset(dw, COUNT>=min.count)[, list(IDN=length(unique(ID))), by=c('PTY_RUN','WLEN','WFROM')]
	tmp			<- tmp[, list(IDN_L=min(IDN), IDN_M=mean(IDN), IDN_U=max(IDN)), by=c('PTY_RUN','WLEN')]
	
	ggplot(tmp, aes(x=WLEN, y=IDN_M, ymin=IDN_L, ymax=IDN_U)) + 
			geom_point() + 
			geom_errorbar() + 
			scale_y_continuous(lim=c(0,75), breaks=seq(0,75,5), expand=c(0,0)) +
			labs(x='window length', y='number of individuals with min 30 reads')
	ggsave(file=plot.file, w=5, h=6)
	ggplot(tmp, aes(x=WFROM, y=IDN)) + geom_bar(stat='identity') + facet_wrap(~WLEN, scales='free_x', ncol=4)
}

TasP.evaluate.bam.len<- function()
{	
	infile.bam		<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/output/bam_stats_170928.rda'
	infile.runs		<- '/Users/Oliver/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/input/ptyRuns_2batchesTest_OR.csv'
	plot.file		<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/output/bam_len_170928.pdf'
	load(infile.bam)
	
	pty.runs		<- as.data.table(read.csv( infile.runs ))
	
	setnames(bam.len, 'FILE', 'SAMPLE_ID')
	set(bam.len,NULL,'SAMPLE_ID',bam.len[, gsub('\\.bam','',SAMPLE_ID)])
	#bam.len			<- bam.len[, list(X=paste0(QU[-length(QU)],'-',QU[-1]-1), PDF=diff(CDF)), by='SAMPLE_ID']
	
	pty.runs		<- unique(pty.runs, by=c('PTY_RUN','SAMPLE_ID'))
	bam.len			<- merge(bam.len, pty.runs, by='SAMPLE_ID')
	setkey(bam.len, PTY_RUN, SAMPLE_ID, QU)
	
	ggplot(subset(bam.len, QU>=40 & QU<320), aes(y=CDF, x=factor(QU), fill=factor(PTY_RUN))) + 
			geom_boxplot() + 
			scale_fill_brewer(palette='Set1') +
			theme_bw() + 
			labs(y='cumulative frequency\nin one individual\n', x='\nlength of quality-trimmed short reads\n(nt)', fill='sequence run') +
			theme(legend.position='bottom')
	ggsave(file=plot.file, w=8, h=6)
}

TasP.evaluate.bam.coverage<- function()
{	
	infile.bam		<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/output/bam_stats_170928.rda'
	infile.runs		<- '/Users/Oliver/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/input/ptyRuns_2batchesTest_OR.csv'
	plot.file		<- '~/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch/output/bam_coverage_170928.pdf'
	load(infile.bam)
	
	pty.runs		<- as.data.table(read.csv( infile.runs ))
	
	setnames(bam.cov, 'FILE_ID', 'SAMPLE_ID')
	pty.runs		<- unique(pty.runs, by=c('PTY_RUN','SAMPLE_ID'))
	bam.cov			<- merge(bam.cov, pty.runs, by='SAMPLE_ID')
	setkey(bam.cov, PTY_RUN, SAMPLE_ID, POS)
	
	ps				<- lapply(bam.cov[, unique(PTY_RUN)], function(ptyr)
			{
				cat('\nprocess run', ptyr)
				tmp				<- subset(bam.cov, PTY_RUN==ptyr)
				tmp				<- tmp[, {
							z	<- rep(COV,REP)
							list(COV=z, POS=seq_along(z), REF=REF[1])
						}, by=c('SAMPLE_ID','PTY_RUN')]
				p	<- ggplot(tmp, aes(x=POS, y=COV, colour=SAMPLE_ID)) + geom_step() + theme_bw() +
						scale_y_log10() +
						facet_grid(~PTY_RUN) +
						labs(x='genome position\n(relative to reference)', y='read coverage', colour='patient')
				p
			})
	#pdf(file=gsub('\\.rda','_coverage_UG.pdf',infile), w=12, h=5)
	pdf(file=plot.file, w=12, h=5)
	for(i in seq_along(ps))
		print(ps[[i]])	
	dev.off()
}

Tasp.pipeline.testbatch.170925.stage1<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	#
	#	INPUT ARGS PLATFORM
	#	
	if(0)
	{
		#fix Tiago's file
		HOME				<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch'
		in.dir				<- file.path(HOME,"input")
		pty.runs			<- as.data.table(read.table( file.path(in.dir, 'ptyRuns_2batchesTest.txt'), header=TRUE ))
		set(pty.runs, NULL, c('PTY_RUN','RID'), NULL)		
		setnames(pty.runs, c('BATCH','SID','RENAME_SID','RefUsedToMap'), c('PTY_RUN','SAMPLE_ID','RENAME_ID','REFERENCE_ID'))
		set(pty.runs, NULL, 'REFERENCE_ID', pty.runs[, paste0(REFERENCE_ID,'.fasta')])
		set(pty.runs, NULL, 'SAMPLE_ID', pty.runs[, gsub('\\.bam','.sorted.bam',SAMPLE_ID)])
		set(pty.runs, NULL, 'RENAME_ID', pty.runs[, paste0('TasP',gsub('^([0-9]+).*','\\1',SAMPLE_ID))])
		write.csv(pty.runs, row.names=FALSE, file=file.path(in.dir, 'ptyRuns_2batchesTest_OR.csv'))		
	}
	if(1)
	{	
		#olli's paths
		#HOME				<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2018_TasP_phyloscanner/170908_test_batch'
		in.dir				<- file.path(HOME,"input")
		work.dir			<- file.path(HOME,"tmp")		
		out.dir				<- file.path(HOME,"output")				
		pty.runs			<- as.data.table(read.csv( file.path(in.dir, 'ptyRuns_2batchesTest_OR.csv'), header=TRUE, stringsAsFactors=FALSE ))
		hpc.load			<- "module load intel-suite/2015.1 mpi R/3.3.2 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"
		#hpc.nproc			<- 4	
		hpc.nproc			<- 1
		#prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner_make_trees.py'
		#pty.data.dir		<- '/work/or105/PANGEA_mapout/data'
		prog.pty			<- '/Users/Oliver/git/phylotypes/phyloscanner_make_trees.py'
		pty.data.dir		<- file.path(HOME,"data")
		#prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-SSE3 -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-SSE3 -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))
		prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-AVX -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-AVX -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))		
		pty.select			<- 1			
	}					
	#
	#	INPUT ARGS PHYLOSCANNER RUN
	#	
	if(1)
	{				
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft='mafft', 
				prog.raxml=prog.raxml, 
				data.dir=pty.data.dir, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=NA,	#to be specified as part of pty.runs --> BACKGROUND_ID
				alignments.root='B.FR.83.HXB2.K03455', 
				alignments.pairwise.to='B.FR.83.HXB2.K03455',
				bl.normalising.reference.file=system.file(package="Phyloscanner.R.utilities", "data", "hiv.hxb2.norm.constants.rda"),
				bl.normalising.reference.var='MEDIAN_PWD',
				tip.regex='^(.*)_read_([0-9]+)_count_([0-9]+)$',				
				window.automatic= '', 
				merge.threshold=2, 
				min.read.count=1, 				
				merge.paired.reads=TRUE, 
				#no.trees=TRUE,
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				dont.check.recombination=TRUE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				#win=c(2300,2800,125,250),
				#win=c(2300,2800,75,150),
				#win=c(2300,2800,100,200),
				#win=c(2300,2800,100,180),
				#win=c(2300,2800,95,190),
				win=c(2300,2800,85,170),			
				keep.overhangs=FALSE,	
				use.blacklisters=c('ParsimonyBasedBlacklister','DownsampleReads'),				
				roguesubtree.kParam=20,
				roguesubtree.prop.threshold=0,
				roguesubtree.read.threshold=20,
				dwns.maxReadsPerPatient=50,	
				multifurcation.threshold=1e-5,
				split.rule='s',
				split.kParam=20,
				split.proximityThreshold=0.035,	
				split.readCountsMatterOnZeroBranches=TRUE,
				pw.trmw.min.reads=20,									
				pw.trmw.min.tips=1,
				pw.trmw.close.brl=0.035,
				pw.trmw.distant.brl=0.08,
				pw.prior.keff=2,
				pw.prior.neff=3,
				pw.prior.keff.dir=2,
				pw.prior.neff.dir=3,				
				pw.prior.calibrated.prob=0.66,
				mem.save=0,
				verbose=TRUE,
				select=pty.select
		)		 
	}	
	#
	#	RUN PHYLOSCANNER
	#
	if(1)
	{
		pty.c				<- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)		
		#pty.c[1,cat(CMD)]		
		invisible(pty.c[,	{
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=171, hpc.q="pqeelab", hpc.mem="23600mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=171, hpc.q="pqeelab", hpc.mem="5600mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=24, hpc.q=NA, hpc.mem="1890mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=199, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile		<- paste("scRA2",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
							#stop()
						}, by='PTY_RUN'])
		quit('no')
	}	
	if(0)	#check failing runs
	{
		df	<- data.table(F=readLines('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/all_files.txt'))
		df[, PTY_RUN:= as.integer(gsub('ptyr([0-9]+)_.*','\\1',F))]
		tmp	<- df[, list(DONE=any(grepl('pairwise',F))), by='PTY_RUN']
		tmp	<- subset(tmp, !DONE)
		merge(tmp, subset(pty.runs, SID=='15080_1_28', c('PTY_RUN','SID')), by='PTY_RUN', all=TRUE)
		
		subset(tmp, !DONE)[, cat(paste(sort(PTY_RUN), collapse='","'))]
	}
}