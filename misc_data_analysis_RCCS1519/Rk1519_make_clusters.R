library(seqinr)
require(data.table)

rkuvri.make.subsequences <- function()
{
  # break sequences into windows and batches
  
  # file names
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  out.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_'
  potential.networks.analysis.dir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100"
  window_size <- 500
  batch_size <- 100
  
  # load data
  alignment <- read.fasta(file = file.path(data.dir,"200422_PANGEA2_RCCSMRC_alignment.fasta"))
  nsequence <- length(alignment)
  npos <- unique(lengths(alignment))  
  
  # define windows 
  windows_last_start <- ceiling(npos/100) * 100 - window_size + 1
  windows_2last_end <- floor(npos/100) * 100 
  windows_start <- seq(1,windows_last_start,100)
  windows_end <- c(seq(window_size ,windows_2last_end,100),npos)
  
  # break sequence into windows and save to subsequence
  for (i in 1:length(windows_start)){
    subsequence <- list()
    subsequence.names <- vector()
    for (j in 1:nsequence){
      subsequence[[j]] <-  alignment[[j]][windows_start[i]:windows_end[i]]
      subsequence.names <- c(subsequence.names,names(alignment)[j])
    }
    names(subsequence) <- subsequence.names
    write.fasta(sequences=subsequence,
                names=names(subsequence),
                nbchar = 60,
                file.out=file.path(potential.networks.analysis.dir ,paste0('subsequence_window_',i,'.fasta')))
  }
  

  # create pairs
  dfo		<- as.data.table(t(combn(names(alignment),2)))
  setnames(dfo, c('V1','V2'), c('TAXA1','TAXA2'))
  
  #	create batches
  batchn <- batch_size
  dfo[, BATCH:= seq_len(nrow(dfo))%%batchn+1L]
  
  # write batches to rds
  cat('\nWriting to file ', file.path(potential.networks.analysis.dir,'batch_windowsize500_batchsize100.rds') )
  saveRDS(dfo,file = file.path(potential.networks.analysis.dir,'batch_windowsize500_batchsize100.rds'))
  
}

# define match score for each position
is.match <- function(x,y){
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


rkuvri.make.distance<- function()
{
  # send jobs
  
  if(1)
  {
    args_dir <- list()
    args_dir[['data_dir']] <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
    args_dir[['out_dir']] <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100'
    args_dir[['script_dir']] <- '~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/'
    args_dir[['batchi']] <- 61
    args_dir[['windowi']] <- 48
  }
  
  args_line <-  as.list(commandArgs(trailingOnly=TRUE))
  if(length(args_line) > 0) 
  {
    stopifnot(args_line[[1]]=='-data_dir')
    stopifnot(args_line[[3]]=='-out_dir')
    stopifnot(args_line[[5]]=='-script_dir')
    stopifnot(args_line[[7]]=='-windowi')
    args_dir <- list()
    args_dir[['data_dir']] <- args_line[[2]]
    args_dir[['out_dir']] <- args_line[[4]]
    args_dir[['script_dir']] <- args_line[[6]]
    args_dir[['windowi']] <- as.integer(args_line[[8]])
  } 
  
  #	read sequences	
  inputFile <- file.path(args_dir[['out_dir']], paste0('subsequence_window_', args_dir[['windowi']], '.fasta'))
  sq<- read.fasta(inputFile)
  
  #	prepare and select batch
  dfo <- readRDS(file.path(args_dir[['out_dir']],paste0('batch_windowsize500_batchsize100.rds')))
  
  # take batchi
  dfo		<- subset(dfo, BATCH==args_dir[['batchi']])
  setkey(dfo, TAXA1, TAXA2)
  
  # remove sequence with all n
  sq		<- sq[ dfo[, unique(c(TAXA1, TAXA2))] ]
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
  saveRDS(dfo, file=paste0(gsub('.fasta$','',inputFile),'_batch',args_dir[['batchi']],'.rds')	)
  
}

#	function to make PBS header
make.PBS.header <- function(hpc.walltime=47, hpc.select=1, hpc.nproc=1, hpc.mem= "6gb", hpc.load= "module load anaconda3/personal\nsource activate covid19model", hpc.q="pqcovid19c", hpc.array=1 )
{	
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


rkuvri.submit.make.distance <- function(){
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  out.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_'
  potential.networks.analysis.dir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100"
  window_size <- 500
  
  # load data
  alignment <- read.fasta(file = file.path(data_dir,"200422_PANGEA2_RCCSMRC_alignment.fasta"))
  nsequence <- length(alignment)
  npos <- unique(lengths(alignment))
  
  # windows 
  windows_last_start <- ceiling(npos/100) * 100 - window_size + 1
  windows_2last_end <- floor(npos/100) * 100 
  windows_start <- seq(1,windows_last_start,100)
  windows_end <- c(seq(args_dir[['window_size']] ,windows_2last_end,100),npos)
  
  # make scripts
  windown <- length(windows_start)
  jobn <- windown * batchn
  windows <- rep(1:windown, times= batchn)
  batchs <- rep(1:batchn, each= windown)
  cmd <-  make.PBS.header(hpc.walltime=71, 
                          hpc.select=1, 
                          hpc.nproc=1, 
                          hpc.mem= "5gb", 
                          hpc.load= "module load anaconda3/personal\nsource activate Renv", 
                          hpc.array= jobn,
                          hpc.q = NA)
  
  cmd <- paste0(cmd ,'\n case $PBS_ARRAY_INDEX in \n')
  for (i in seq_len(jobn)) {
    cmd <- paste0(cmd, i,') \n')
    cmd <- paste0(cmd, 'Rscript ', args_dir[['script_dir']], 'distance.R  -data_dir ',  args_dir[['data_dir']], 
                  ' -out_dir ', paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ),
                  ' -script_dir ', args_dir[['script_dir']],
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

rkuvri.process.distance <- function(){
  
}

