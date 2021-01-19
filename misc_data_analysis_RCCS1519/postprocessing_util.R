library(data.table)

if(1)
{
  args_dir <- list()
  args_dir[['data_dir']] <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  args_dir[['out_dir']] <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY'
  args_dir[['window_size']] <- 500
  args_dir[['batch_size']] <- 100
  args_dir[['script_dir']] <- '~/phyloscanner/data_analysis_RCCS1519'
  args_dir[['job_tag']] <- 'windowsize500_batchsize100'
  args_dir[['windowi']] <- 83
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
  stopifnot(args_line[[13]]=='-windowi')
  args_dir <- list()
  args_dir[['data_dir']] <- args_line[[2]]
  args_dir[['out_dir']] <- args_line[[4]]
  args_dir[['window_size']] <- as.integer(args_line[[6]])
  args_dir[['batch_size']] <- as.integer(args_line[[8]])
  args_dir[['script_dir']] <- args_line[[10]]
  args_dir[['job_tag']] <- args_line[[12]]
  args_dir[['windowi']] <- as.integer(args_line[[14]])
} 

setwd(paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ))

# gather outputs for each window
base <- paste0('subsequence_window_',args_dir[['windowi']])
distancei <- data.table()
for (i in 1:args_dir[['batch_size']]) {
  tmp <- readRDS(paste0(base,'_batch',i,'.rds'))
  distancei <- rbind(distancei,tmp)
#  file.remove(paste0(base,'_batch',i,'.rds'))
}

# summarise
metai <- data.table(window=args_dir[['windowi']],
                    all_pairs=nrow(distancei),
                    valid_pairs=nrow(distancei[LENGTH!=-1]),
                    length200=nrow(distancei[LENGTH>=200]),
                    length250=nrow(distancei[LENGTH>=250]),
                    length300=nrow(distancei[LENGTH>=300]),
                    max_length=max(distancei[LENGTH>=0,LENGTH]),
                    min_length=min(distancei[LENGTH>=0,LENGTH]),
                    max_percentage=max(distancei[PERC>=0,PERC]),
                    min_percentage=min(distancei[PERC>=0,PERC]))

# remove invalid pairs
distancei <- distancei[LENGTH>=0]
distancei[,WINDOW:=args_dir[['windowi']]]

save(distancei,metai,file=paste0(base,'_results.rda'))






