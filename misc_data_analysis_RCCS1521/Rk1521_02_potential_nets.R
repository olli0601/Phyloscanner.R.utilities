rk.seq.make.consensus.alignment <- function(){
  #
  library(ape)
  library(data.table)
  indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
  infile <- file.path(indir, '240809_PANGEA2_RCCS_final_samples_minrun125_minnuc250.rds')
  infile_new <- file.path(indir, "240809_PANGEA2_RCCS_final_samples_minrun125_minnuc250_minusPoor1.rds")
  infile_poor <- file.path(indir,"240809_PANGEA2_RCCS_final_samples_minrun125_minnuc250_all_seqs_of_ind_poor.csv")
  outfile <- file.path(indir, '240809_PANGEA2_RCCS_final_samples_alignment.fasta')
  outfile_new <- file.path(indir, '240809_PANGEA2_RCCS_final_samples_alignment_minusPoor1.fasta')
  
  #
  ds <- readRDS(infile)
  dsn <- readRDS(infile_new)
  dp <- read.csv(infile_poor)
  for(i in 1:nrow(ds)){
    if(i%%10 == 0){cat("processing sequence ", i, '\n')}
    if(i==1){
      fas <- read.dna(ds$FORGLOBALALN_FASTA[i], format='fasta')
    }else{
      fas <- rbind(fas, read.dna(ds$FORGLOBALALN_FASTA[i], format='fasta'))
      if(length(fas[-1,])!=length(fas[1,]))warnings(paste0("lengths do not match for sequence ", i))
    }
  }
  
  write.dna(fas, file=outfile, format='fasta', colsep='', nbcol=-1)
  
  fas <- read.dna(file=outfile, format = "fasta")
  # fix names
  rownames(fas) <- gsub("UG50406_PanRak7","UG504060_PanRak7",rownames(fas))
  fasn <- fas[rownames(fas) %in% gsub("(consensus).*", "\\1",basename(dsn$FORGLOBALALN_FASTA)),]
  write.dna(fasn, file=outfile_new, format='fasta', colsep='', nbcol=-1)
  
  # str(fasn)
  # tmp <- gsub("(consensus).*", "\\1",basename(dsn$FORGLOBALALN_FASTA))
  # tmp[!tmp %in% rownames(fasn)]
  # "UG504060_PanRak7-E1A1L1P1S1C1R1_consensus"
  # which(rownames(fas) == "UG504060_PanRak7-E1A1L1P1S1C1R1_consensus")
  # grep("UG504060",rownames(fas))
  # grep("UG504060",basename(ds$FORGLOBALALN_FASTA))
  # fas[5059,]
  # dsn[grep("UG50406", dsn$FORGLOBALALN_FASTA),]
}

# helper function
is.match <- function(x,y){
  #' returns match score for each position
  #' @param x,y two nucleotides at one positions
  if(x %in% c('a','g','c','t') & y %in% c('a','g','c','t')){
    if(x==y){
      ans = 1
    }else{
      ans = 0
    }
  }else{
    if(((x=='g'|x=='a')&y=='r')|
       ((y=='g'|y=='a')&x=='r')|
       ((x=='t'|x=='c')&y=='y')|
       ((y=='t'|y=='c')&x=='y')|
       ((x=='a'|x=='t')&y=='w')|
       ((y=='a'|y=='t')&x=='w')|
       ((x=='a'|x=='c')&y=='m')|
       ((y=='a'|y=='c')&x=='m')|
       ((x=='g'|x=='c')&y=='s')|
       ((y=='g'|y=='c')&x=='s')|
       ((x=='g'|x=='t')&y=='k')|
       ((y=='g'|y=='t')&x=='k')){
      ans = 1/2
    }else if(((x=='a'|x=='c'|x=='t')&y=='h')|
             ((y=='a'|y=='c'|y=='t')&x=='h')|
             ((x=='a'|x=='g'|x=='t')&y=='d')|
             ((y=='a'|y=='g'|y=='t')&x=='d')|
             ((x=='g'|x=='c'|x=='t')&y=='b')|
             ((y=='g'|y=='c'|y=='t')&x=='b')|
             ((x=='a'|x=='c'|x=='g')&y=='v')|
             ((y=='a'|y=='c'|y=='g')&x=='v')){
      ans = 1/3
    }else{
      ans = 0
    }
  }
  return(ans)
}

# helper function
make.PBS.header <- function(hpc.walltime=47, hpc.select=1, hpc.nproc=1, hpc.mem= "6gb", hpc.load= "module load anaconda3/personal", hpc.q="pqeelab", hpc.array=1 )
{	
  #'  function to make PBS header
  pbshead <- "#!/bin/sh"
  tmp <- paste("#PBS -l walltime=", hpc.walltime, ":59:00", sep = "")
  pbshead <- paste(pbshead, tmp, sep = "\n")
  tmp <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":ompthreads=", hpc.nproc,":mem=", hpc.mem, sep = "")	
  pbshead <- paste(pbshead, tmp, sep = "\n")
  pbshead <- paste(pbshead, "#PBS -j oe", sep = "\n")
  if(hpc.array>1)
  {
    pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
  }				
  if(!is.na(hpc.q))
  {
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  }		
  pbshead	<- paste(pbshead, hpc.load, sep = "\n")
  pbshead
}

rk.seq.make.subsequences <- function()
{
  #' break sequences into subsequences by windows
  #' it takes sequences and break them into overlapping windows of length 500 and defines batches for distance calculation
  
  # break sequences into windows and batches
  window_size <- 500
  batch_size <- 100
  
  #
  indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
  outdir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/'
  infile <- file.path(indir, '240809_PANGEA2_RCCS_final_samples_alignment_minusPoor1.fasta')
  alignment <- read.dna(infile, format='fasta')
  if(!dir.exists(outdir)) dir.create(outdir)
  
  # 
  nsequence <- nrow(alignment)
  npos <- ncol(alignment)
  
  # define windows 
  windows_last_start <- ceiling(npos/100) * 100 - window_size + 1
  windows_2last_end <- floor(npos/100) * 100 
  windows_start <- seq(1,windows_last_start,100)
  windows_end <- c(seq(window_size ,windows_2last_end,100),npos)
  
  # break sequence into windows and save to subsequence
  for (i in 1:length(windows_start)){
    subsequence <- alignment[, windows_start[i]:windows_end[i]]
    write.dna(subsequence, file=paste0(outdir, "subsequence", i, ".fasta"), format='fasta', colsep='', nbcol=-1)
  }
  
  
  # create pairs
  dfo	<- as.data.table(t(combn(rownames(alignment),2)))
  setnames(dfo, c('V1','V2'), c('TAXA1','TAXA2'))
  
  #	create batches
  dfo[, BATCH:= seq_len(nrow(dfo))%%batch_size+1L]
  
  # write batches to rds
  cat('\nWriting to file ', file.path(outdir,'batch_size100.rds') )
  saveRDS(dfo,file = file.path(outdir,'batch_size100.rds'))
  
}

rk.seq.make.distance <- function()
{
  #'  calculate length of valid sequences, total match scores, percentage of matches of pair of sequences in one batch and in one window
  #' it inputs the subsequence and batch info, and returns a data.table with distances
  
  library(seqinr)
  library(data.table)
  if(1)
  {
    args_dir <- list()
    indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
    args_dir[['data_dir']] <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
    args_dir[['out_dir']] <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/'
    args_dir[['batchi']] <- 61
    args_dir[['windowi']] <- 48
  }
  
  args_line <-  as.list(commandArgs(trailingOnly=TRUE))
  if(length(args_line) > 0) 
  {
    stopifnot(args_line[[1]]=='-data_dir')
    stopifnot(args_line[[3]]=='-out_dir')
    stopifnot(args_line[[5]]=='-batchi')
    stopifnot(args_line[[7]]=='-windowi')
    args_dir <- list()
    args_dir[['data_dir']] <- args_line[[2]]
    args_dir[['out_dir']] <- args_line[[4]]
    args_dir[['batchi']] <- as.integer(args_line[[6]])
    args_dir[['windowi']] <- as.integer(args_line[[8]])
  } 
  
  # file names
  infileb <- file.path(args_dir[['out_dir']],paste0('batch_size100.rds'))
  infiles <- file.path(args_dir[['out_dir']], paste0('subsequence', args_dir[['windowi']], '.fasta'))
  
  #	load
  sq<- read.fasta(infiles)
  dfo <- readRDS(infileb)
  
  # take batchi
  dfo	<- subset(dfo, BATCH==args_dir[['batchi']])
  setkey(dfo, TAXA1, TAXA2)
  
  # remove sequence with all n
  sq	<- sq[ dfo[, unique(c(TAXA1, TAXA2))] ]
  sq.n <- unlist(lapply(sq, function(x){sum(x=='?' | x=='-' | x=='n')==500}))
  sq.n.names <- names(sq.n[sq.n==TRUE])
  tmp <- dfo[TAXA1 %in% sq.n.names|TAXA2 %in% sq.n.names]
  tmp[,LENGTH:=-1.0]
  tmp[,MATCHES:=-1.0]
  tmp[,PERC:=-1.0]
  dfo <- dfo[!TAXA1 %in% sq.n.names & !TAXA2 %in% sq.n.names]
  
  
  #	get similarity
  dfo		<- dfo[,{ 
    seq1 <- sq[[TAXA1]]
    seq2 <- sq[[TAXA2]]
    pos <- which(seq1 !='-' & seq1 !='?' & seq1 !='n' & seq2 !='-' & seq2 !='?' &  seq2 !='n')
    if (length(pos)==0 ){
      len <- as.integer(-1.0)
      match <- -1.0
    }else{
      seq1 <- seq1[pos]
      seq2 <- seq2[pos]
      len <- length(pos)
      # cat(TAXA1,'\n',seq1,'\n',TAXA2,'\n',seq2,'\n')
      match <- sum(sapply(1:length(pos),function(pos){is.match(seq1[pos],seq2[pos])}))
    }
    list(LENGTH= len,
         MATCHES= match)
  }, by=c('BATCH','TAXA1','TAXA2')]
  dfo[,PERC:=MATCHES/LENGTH]
  dfo <- rbind(dfo,tmp)
  
  #	save
  saveRDS(dfo, file=paste0(gsub('.fasta$','',infiles),'_batch',args_dir[['batchi']],'.rds')	)
  
}



rkuvri.submit.make.distance <- function(){
  #' write script to calculate distance
  
  window_size <- 500
  batch_size <- 100
  
  # files
  #
  indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
  outdir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/'
  infile <- file.path(indir, '240809_PANGEA2_RCCS_final_samples_alignment_minusPoor1.fasta')
  
  # load data
  alignment <- read.dna(file = infile, format = "fasta")
  # 
  nsequence <- nrow(alignment)
  npos <- ncol(alignment)
  
  # define windows 
  windows_last_start <- ceiling(npos/100) * 100 - window_size + 1
  windows_2last_end <- floor(npos/100) * 100 
  windows_start <- seq(1,windows_last_start,100)
  windows_end <- c(seq(window_size ,windows_2last_end,100),npos)
  
  # make scripts
  windown <- length(windows_start)
  jobn <- windown * batch_size
  windows <- rep(1:windown, times= batch_size)
  batchs <- rep(1:batch_size, each= windown)
  cmd <-  make.PBS.header(hpc.walltime=71, 
                          hpc.select=1, 
                          hpc.nproc=1, 
                          hpc.mem= "5gb", 
                          hpc.load= "module load anaconda3/personal\nsource activate MakePotentialSubpgraphs", 
                          hpc.array= jobn,
                          hpc.q = NA)
  
  cmd <- paste0(cmd ,'\n case $PBS_ARRAY_INDEX in \n')
  for (i in seq_len(jobn)) {
    cmd <- paste0(cmd, i,') \n')
    cmd <- paste0(cmd, 'Rscript ', '~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1521/distance.R -data_dir ',  data.dir, 
                  ' -out_dir ', outdir,
                  ' -batchi ', batchs[i],
                  ' -windowi ', windows[i], '\n')
    cmd <- paste0(cmd, ';; \n')
  }
  cmd <- paste0(cmd, 'esac \n')
  
  # write submission file	
  processing.file <- file.path(potential.networks.analysis.dir, 'processing.sh')
  cat(cmd, file=processing.file)
  
  # set permissions
  Sys.chmod(processing.file, mode='644')	
  
  # run
  cmd <-''
  cmd 	<- paste0(cmd, '\tcd ', dirname(processing.file),'\n')
  cmd 	<- paste0(cmd,'\tqsub ', basename(processing.file),'\n')
  cat(cmd)
  
}