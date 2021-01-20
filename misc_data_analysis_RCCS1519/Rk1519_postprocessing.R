require(seqinr)
require(data.table)


if(1)
{
  args_dir <- list()
  args_dir[['data_dir']] <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  args_dir[['out_dir']] <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY'
  args_dir[['window_size']] <- 500
  args_dir[['batch_size']] <- 100
  args_dir[['script_dir']] <-'~/phyloscanner/data_analysis_RCCS1519/'
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

# load data
alignment <- read.fasta(file = file.path(args_dir[['data_dir']],"200422_PANGEA2_RCCSMRC_alignment.fasta"))
nsequence <- length(alignment)
npos <- unique(lengths(alignment))

# windows
windows_last_start <- ceiling(npos/100) * 100 - args_dir[['window_size']] + 1
windows_2last_end <- floor(npos/100) * 100 
windows_start <- seq(1,windows_last_start,100)
windows_end <- c(seq(args_dir[['window_size']] ,windows_2last_end,100),npos)
windows <- seq_len(length(windows_start))


cat("---------------------- create script for postprocess: -----------------------\n ")
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
cmd <-  make.PBS.header(hpc.walltime=71, 
                         hpc.select=1, 
                         hpc.nproc=1, 
                         hpc.mem= "50gb", 
                         hpc.load= "module load anaconda3/personal\nsource activate Renv", 
                         hpc.array= windown,
                         hpc.q=NA)

cmd <- paste0(cmd, '\n case $PBS_ARRAY_INDEX in \n')
for (i in seq_len(windown)) {
  cmd <- paste0(cmd, i,') \n')
  cmd <- paste0(cmd, 'Rscript ', args_dir[['script_dir']], 'postprocessing_util.R  -data_dir ',  args_dir[['data_dir']], 
                ' -out_dir ', args_dir[['out_dir']],
                ' -window_size ', args_dir[['window_size']],
                ' -batch_size ', args_dir[['batch_size']],
                ' -script_dir ', args_dir[['script_dir']],
                ' -job_tag ', args_dir[['job_tag']],
                ' -windowi ', i, '\n')
  cmd <- paste0(cmd, ';; \n')
}



cmd <- paste0(cmd, 'esac \n')
# write submission file	
collect.results.file <- file.path(paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ), 'collect.results.sh')
cat(cmd, file=collect.results.file)
# set permissions
Sys.chmod(collect.results.file, mode='644')	


# 
cmd <-  make.PBS.header(hpc.walltime=71, 
                         hpc.select=1, 
                         hpc.nproc=1, 
                         hpc.mem= "96gb", 
                         hpc.load= "module load anaconda3/personal\nsource activate Renv",
                         hpc.q=NA)
cmd <- paste0(cmd, '\n Rscript ', file.path(args_dir[['script_dir']]), 'threshold.R  -data_dir ',  args_dir[['data_dir']], 
              ' -out_dir ', args_dir[['out_dir']],
              ' -window_size ', args_dir[['window_size']],
              ' -batch_size ', args_dir[['batch_size']],
              ' -script_dir ', args_dir[['script_dir']],
              ' -job_tag ', args_dir[['job_tag']],
              ' -windown ', length(windows_start), '\n',
              ' -length_cutoff ', 250, '\n')

# write submission file	
threshold.file <- file.path(paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ), 'threshold.sh')
cat(cmd, file=threshold.file)
# set permissions
Sys.chmod(threshold.file, mode='644')	

      
#
cmd <-  make.PBS.header(hpc.walltime=71, 
                        hpc.select=1, 
                        hpc.nproc=1, 
                        hpc.mem= "96gb", 
                        hpc.load= "module load anaconda3/personal\nsource activate Renv",
                        hpc.q=NA)
cmd <- paste0(cmd, '\n Rscript ', file.path(args_dir[['script_dir']]), 'classification_sequence.R  -data_dir ',  args_dir[['data_dir']], 
              ' -out_dir ', args_dir[['out_dir']],
              ' -window_size ', args_dir[['window_size']],
              ' -batch_size ', args_dir[['batch_size']],
              ' -script_dir ', args_dir[['script_dir']],
              ' -job_tag ', args_dir[['job_tag']],
              ' -windown ', length(windows_start), '\n',
              ' -length_cutoff ', 250, '\n',
              ' -method ', 2, '\n')

# write submission file	
classification.sequence.file <- file.path(paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ), 'classification.sequence.sh')
cat(cmd, file=classification.sequence.file)
# set permissions
Sys.chmod(classification.sequence.file, mode='644')	


#
cmd <-  make.PBS.header(hpc.walltime=71, 
                        hpc.select=1, 
                        hpc.nproc=1, 
                        hpc.mem= "96gb", 
                        hpc.load= "module load anaconda3/personal\nsource activate Renv",
                        hpc.q=NA)
cmd <- paste0(cmd, '\n Rscript ', file.path(args_dir[['script_dir']]), 'classification_person.R  -data_dir ',  args_dir[['data_dir']], 
              ' -out_dir ', args_dir[['out_dir']],
              ' -window_size ', args_dir[['window_size']],
              ' -batch_size ', args_dir[['batch_size']],
              ' -script_dir ', args_dir[['script_dir']],
              ' -job_tag ', args_dir[['job_tag']],
              ' -windown ', length(windows_start), '\n',
              ' -length_cutoff ', 250, '\n',
              ' -method ', 2, '\n')

# write submission file	
classification.person.file <- file.path(paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ), 'classification.person.sh')
cat(cmd, file=classification.person.file)
# set permissions
Sys.chmod(classification.person.file, mode='644')	



#
cmd <-  make.PBS.header(hpc.walltime=71, 
                        hpc.select=1, 
                        hpc.nproc=1, 
                        hpc.mem= "96gb", 
                        hpc.load= "module load anaconda3/personal\nsource activate Renv",
                        hpc.q=NA)
cmd <- paste0(cmd, '\n Rscript ', file.path(args_dir[['script_dir']]), 'validate_classification.R  -data_dir ',  args_dir[['data_dir']], 
              ' -out_dir ', args_dir[['out_dir']],
              ' -window_size ', args_dir[['window_size']],
              ' -batch_size ', args_dir[['batch_size']],
              ' -script_dir ', args_dir[['script_dir']],
              ' -job_tag ', args_dir[['job_tag']],
              ' -windown ', length(windows_start), '\n',
              ' -length_cutoff ', 250, '\n',
              ' -method ', 2, '\n')

# write submission file	
classification.validate.file <- file.path(paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ), 'classification.validate.sh')
cat(cmd, file=classification.validate.file)
# set permissions
Sys.chmod(classification.validate.file, mode='644')	


# run
cmd <-''
cmd 	<- paste0(cmd, '\tcd ', dirname(collect.results.file),'\n')
cmd 	<- paste0(cmd,'\tqsub ', basename(collect.results.file),'\n')
cmd 	<- paste0(cmd,'\tqsub ', basename(threshold.file),'\n')
cmd 	<- paste0(cmd,'\tqsub ', basename(classification.sequence.file),'\n')
cmd 	<- paste0(cmd,'\tqsub ', basename(classification.person.file),'\n')
cmd 	<- paste0(cmd,'\tqsub ', basename(classification.validate.file),'\n')
cat(cmd)



