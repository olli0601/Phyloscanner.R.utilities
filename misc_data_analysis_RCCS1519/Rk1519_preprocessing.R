library(seqinr)
require(data.table)

if(1)
{
  args_dir <- list()
  args_dir[['data_dir']] <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  args_dir[['out_dir']] <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY'
  args_dir[['window_size']] <- 500
  args_dir[['batch_size']] <- 100
  args_dir[['script_dir']] <- '~/phyloscanner/data_analysis_RCCS1519/'
  args_dir[['job_tag']] <- 'windowsize500_batchsize100'
}


args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-data_dir')
  stopifnot(args_line[[3]]=='-out_dir')
  stopifnot(args_line[[5]]=='-window_size')
  stopifnot(args_line[[7]]=='-batch_size')
  stopifnot(args_line[[9]]=='-script_dir')
  stopifnot(args_line[[11]]=='-job_tag')
  args_dir <- list()
  args_dir[['data_dir']] <- args_line[[2]]
  args_dir[['out_dir']] <- args_line[[4]]
  args_dir[['window_size']] <- as.integer(args_line[[6]])
  args_dir[['batch_size']] <- as.integer(args_line[[8]])
  args_dir[['script_dir']] <- args_line[[10]]
  args_dir[['job_tag']] <- args_line[[12]]
} 

cat("---------------------- prepare windows: -----------------------\n ")
# create dir
dir.create(paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ))
setwd(paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ))

# load data
alignment <- read.fasta(file = file.path(args_dir[['data_dir']],"200422_PANGEA2_RCCSMRC_alignment.fasta"))
nsequence <- length(alignment)
npos <- unique(lengths(alignment))

# windows 
windows_last_start <- ceiling(npos/100) * 100 - args_dir[['window_size']] + 1
windows_2last_end <- floor(npos/100) * 100 
windows_start <- seq(1,windows_last_start,100)
windows_end <- c(seq(args_dir[['window_size']] ,windows_2last_end,100),npos)

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
              file.out=file.path(paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ) ,paste0('subsequence_window_',i,'.fasta')))
}

cat("---------------------- create batches: -----------------------\n ")

# create pairs
dfo		<- as.data.table(t(combn(names(alignment),2)))
setnames(dfo, c('V1','V2'), c('TAXA1','TAXA2'))

#	create batch
batchn <- args_dir[['batch_size']]
dfo[, BATCH:= seq_len(nrow(dfo))%%batchn+1L]
saveRDS(dfo,file = file.path(args_dir[['script_dir']],paste0('batch_', args_dir[['job_tag']] ,'.rds')))


cat("---------------------- create script for similarity calculation: -----------------------\n ")
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
                ' -window_size ', args_dir[['window_size']],
                ' -batch_size ', args_dir[['batch_size']],
                ' -script_dir ', args_dir[['script_dir']],
                ' -job_tag ', args_dir[['job_tag']],
                ' -batchi ', batchs[i],
                ' -windowi ', windows[i], '\n')
  cmd <- paste0(cmd, ';; \n')
}
cmd <- paste0(cmd, 'esac \n')
# write submission file	
processing.file <- file.path(paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ), 'processing.sh')
cat(cmd, file=processing.file)
# set permissions
Sys.chmod(processing.file, mode='644')	
# run
cmd2 <-''
cmd2 	<- paste0(cmd2, '\tcd ', dirname(processing.file),'\n')
cmd2 	<- paste0(cmd2,'\tqsub ', basename(processing.file),'\n')
cat(cmd2)
