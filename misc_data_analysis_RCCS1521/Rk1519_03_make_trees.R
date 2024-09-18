rkuvri.make.alignments<-function()
{
  #
  #	set up working environment
  require(Phyloscanner.R.utilities)
  library(data.table)
  library(seqinr)
  
  #
  set.seed(42)
  prog.pty <- '~/phyloscanner/phyloscanner_make_trees.py'
  
  indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
  HOME <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/'
  netdir <- paste0(HOME, "potential_nets/")
  
  in.dir <- file.path(HOME,'240919_phsc_input')
  work.dir <- file.path(HOME,"240919_phsc_work")
  out.dir <- file.path(HOME,"240919_phsc_output")
  package.dir <- file.path(.libPaths(),'Phyloscanner.R.utilities')
  dir.create(in.dir)
  dir.create(work.dir)
  dir.create(out.dir)
  
  #
  downsample <- 200
  infile.runs <- paste0(HOME,'240809_RCCSUVRI_phscinput_runs.rds')
  infile.consensus <- file.path(in.dir,'ConsensusGenomes.fasta')
  infile.consensus.oneeach <- file.path(in.dir,'2019_New_ConsensusGenomesOneEach_GeneCut.fasta')
  infile.hxb2.package <- file.path(package.dir,'HIV1_compendium_AD_B_CPX_v2.fasta')
  
  # load runs
  pty.runs <- readRDS(infile.runs)
  max.per.run <- 4900
  
  # remove same bam files per run
  setorder(pty.runs,PTY_RUN,-ID_TYPE,UNIT_ID)
  tmp <- pty.runs[,duplicated(SAMPLE_ID),by='PTY_RUN']
  tmp <- tmp[,which(V1)]
  pty.runs <- pty.runs[-tmp,]
  tmp <- pty.runs[,length(SAMPLE_ID)-length(unique(SAMPLE_ID)),by='PTY_RUN']
  stopifnot(all(tmp$V1==0))
  
  # load consensus
  consensus_seq<- seqinr::read.fasta(infile.consensus)
  consensus_seq_names <- names(consensus_seq)
  
  # take hxb2
  hxb2 <- grep('HXB2',names(consensus_seq),value = T)
  hxb2_seq <- consensus_seq[[hxb2]]
  
  # compare hxb2
  package_seq <- read.fasta(infile.hxb2.package)
  package_hxb2_seq <- package_seq[[1]]
  package_hxb2_seq <- as.character(package_hxb2_seq)[as.character(package_hxb2_seq)!='-']
  hxb2_seq <- as.character(hxb2_seq)[as.character(hxb2_seq)!='-']
  stopifnot(all((hxb2_seq==package_hxb2_seq)))
  # take root
  root.seq <-grep('^REF_CON_M$|REF_CONSENSUS_M$|^REF_CON_H$|REF_CONSENSUS_H$', names(consensus_seq),value = T)
  
  # adapt format
  pty.runs[ID_TYPE=='control',UNIT_ID:=paste0('CNTRL-',UNIT_ID)]
  pty.runs[ID_TYPE=='control',RENAME_ID:=paste0('CNTRL-',RENAME_ID)]
  pty.runs[,BAM:=paste0(data.dir,SAMPLE_ID,'.bam')]
  pty.runs[,REF:=paste0(data.dir,SAMPLE_ID,'_ref.fasta')]
  setkey(pty.runs,PTY_RUN,RENAME_ID)
  # pty.runs <-  pty.runs[1:2,]
  # ptyi <- c(850,900)
  
  #	define phyloscanner input args to generate read alignments
  #	for each window and each run
  # ptyi <- seq(800,9400,25)
  ptyi <- seq(800,9175,25) # exclude ends
  ptyi <- c( ptyi[ptyi <= 6615-250],6825,6850,ptyi[ptyi >= 7636]) # exclude vloop except 2 windows without gaps
  # prog.mafft='mafft',
  pty.c	<- lapply(seq_along(ptyi), function(i)
  {
    cat('---------------------------------------------------------------- \n')
    print(i)
    pty.args <- list(prog.pty=prog.pty,
                     prog.mafft='\" mafft --globalpair --maxiterate 1000 \" ',
                     data.dir=data.dir,
                     work.dir=work.dir,
                     out.dir=out.dir,
                     alignments.file=infile.consensus,
                     alignments.root=root.seq,
                     alignments.pairwise.to=hxb2,
                     window.automatic= '',
                     merge.threshold=0,
                     min.read.count=1,
                     quality.trim.ends=23,
                     min.internal.quality=23,
                     merge.paired.reads=TRUE,
                     discard.improper.pairs=TRUE,
                     no.trees=TRUE,
                     dont.check.duplicates=FALSE,
                     dont.check.recombination=TRUE,
                     num.bootstraps=1,
                     all.bootstrap.trees=TRUE,
                     strip.max.len=350,
                     min.ureads.individual=NA,
                     win=c(ptyi[i],ptyi[i]+250,25,250),
                     keep.overhangs=FALSE,
                     mem.save=0,
                     verbose=TRUE,
                     select=NA,
                     default.coord=TRUE,
                     realignment=TRUE
    )
    pty.c <- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)
    pty.c[, W_FROM:= ptyi[i]]
    # pty.c[, PTY_RUN:= as.integer(sub('.*ptyr([0-9])_.*','\\1',CMD))]
    pty.c
  })
  pty.c	<- do.call('rbind', pty.c)
  setkey(pty.c,PTY_RUN,W_FROM)
  pty.c[, CASE_ID:= rep(1:max.per.run,times=ceiling(nrow(pty.c)/max.per.run))[1:nrow(pty.c)]]
  pty.c[, JOB_ID:= rep(1:ceiling(nrow(pty.c)/max.per.run),each=max.per.run)[1:nrow(pty.c)]]
  save(pty.c,file='~/ptyc_240919.rda')
  
  #	define PBS variables
  hpc.load			<- "module load intel-suite/2015.1 mpi raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"	# make third party requirements available
  hpc.select			<- 1						# number of nodes
  hpc.nproc			<- 1						# number of processors on node
  hpc.walltime		<- 123						# walltime
  hpc.q				<- "pqeelab"				# PBS queue
  hpc.mem				<- "6gb" 					# RAM
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
  for (i in 1:pty.c[, max(JOB_ID)]) {
    tmp<-pty.c[JOB_ID==i,]
    cmd<-tmp[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
    cmd<-cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]
    cmd<-paste(pbshead,cmd,sep='\n')
    #	submit job
    outfile<-gsub(':','',paste("readali",paste0('job',i),paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
    outfile<-file.path(work.dir, outfile)
    cat(cmd, file=outfile)
    cmd<-paste("qsub", outfile)
    cat(cmd)
    cat(system(cmd, intern= TRUE))
  }
  
  # library(ape)
  # # p <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210211_phsc_output/ptyr1_trees'
  # dir.create(file.path(p,'plots'))
  # files <-list.files(p)
  # files <- grep('fasta$',files,value = T)
  # for (i in 1:length(files)) {
  #   tmp <- read.dna(file.path(p,files[i]),format = 'fasta')
  #   pdf(file.path(p,'plots',gsub('.fasta','_alignment.pdf',files[i])),width = 10, height = 10)
  #   par(cex=0.5)
  #   checkAlignment(tmp)
  #   dev.off()
  # }
  
}


rkuvri.make.trees<- function()
{
  require(data.table)
  require(Phyloscanner.R.utilities)
  #
  #	produce trees
  #
  # lightweight run
  if(0)
  {
    hpc.select<- 1; hpc.nproc<- 1; hpc.walltime<- 4; hpc.mem<- "1850mb"; hpc.q<- NA
  }
  # midweight run
  if(1)
  {
    hpc.select<- 1; hpc.nproc<- 1; hpc.walltime<- 23; hpc.mem<- "1850mb"; hpc.q<- NA
  }
  # heavyweight run
  if(0)
  {
    # hpc.select<- 1; hpc.nproc<- 1; hpc.walltime<- 71; hpc.mem<- "5900mb"; hpc.q<- "pqeelab"
    hpc.select<- 1; hpc.nproc<- 1; hpc.walltime<- 71; hpc.mem<- "63850mb"; hpc.q<- NA
  }
  
  HOME <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  max.per.run <- 4900
  iqtree.pr <- 'iqtree'
  iqtree.args			<- ifelse(hpc.nproc==1, '-m GTR+F+R6 -ntmax 1 -seed 42 -o REF_CON_H',
                          paste0('-m GTR+F+R6 -ntmax ',hpc.nproc,' -seed 42 -o REF_CON_H'))
  in.dir				<- file.path(HOME,'240919_phsc_output')
  out.dir				<- in.dir
  work.dir			<- file.path(HOME,"240919_phsc_work")
  
  infiles	<- data.table(FI=list.files(in.dir, pattern='_v2.fasta$', full.names=TRUE, recursive=TRUE))
  infiles[, FO:= gsub('.fasta$','',FI)]
  infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FI)))]
  infiles[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FI)))]
  infiles[is.na(W_FROM),W_FROM:= as.integer(gsub('.*PositionsExcised_([0-9]+)_.*','\\1',basename(FI)))]
  setkey(infiles, PTY_RUN, W_FROM)
  
  # infiles	<- subset(infiles, PTY_RUN<=102) # job1-6
  # infiles	<- subset(infiles, PTY_RUN>102 & PTY_RUN <=204) # job7-12
  # infiles	<- subset(infiles, PTY_RUN>204 & PTY_RUN <=307) # job13-18
  infiles	<- subset(infiles, PTY_RUN>307) # job18-24
  
  df<- infiles[, list(CMD=cmd.iqtree(FI, outfile=FO, pr=iqtree.pr, pr.args=iqtree.args)), by=c('PTY_RUN','W_FROM')]
  df[, ID:=ceiling(seq_len(nrow(df))/4)]
  df<- df[, list(CMD=paste(CMD, collapse='\n',sep='')), by='ID']
  
  #	create PBS job array
  if(nrow(df) > max.per.run){
    df[, CASE_ID:= rep(1:max.per.run,times=ceiling(nrow(df)/max.per.run))[1:nrow(df)]]
    df[, JOB_ID:= rep(1:ceiling(nrow(df)/max.per.run),each=max.per.run)[1:nrow(df)]]
    # set(df, ID, NULL, NULL)
    df[,ID:=NULL]
  }else{
    setnames(df, 'ID', 'CASE_ID')
    df[, JOB_ID:=1]
  }
  
  pbshead	<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=NULL)
  
  
  for (i in 1:df[, max(JOB_ID)]) {
    tmp<-df[JOB_ID==i,]
    hpc.array <- nrow(tmp)
    cmd<-tmp[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
    cmd<-cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]
    tmp <- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
    tmp <- paste(tmp,'module load anaconda3/personal \n source activate phylo' ,sep='\n')
    cmd<-paste(tmp,cmd,sep='\n')
    #	submit job
    outfile	<- paste("srx",paste0('job',i),paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.')
    outfile<-file.path(work.dir, outfile)
    cat(cmd, file=outfile)
    cmd<-paste("qsub", outfile)
    cat(cmd)
  }
}