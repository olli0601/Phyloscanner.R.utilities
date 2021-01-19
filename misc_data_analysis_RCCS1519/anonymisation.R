library(data.table)
if(1)
{
  args_dir <- list()
  args_dir[['data_dir']] <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  args_dir[['out_dir']] <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY'
  args_dir[['script_dir']] <- '~/phyloscanner/data_analysis_RCCS1519'
  args_dir[['job_tag']] <- 'windowsize500_batchsize100'
}


args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-data_dir')
  stopifnot(args_line[[3]]=='-out_dir')
  stopifnot(args_line[[5]]=='-script_dir')
  stopifnot(args_line[[7]]=='-job_tag')
  args_dir <- list()
  args_dir[['data_dir']] <- args_line[[2]]
  args_dir[['out_dir']] <- args_line[[4]]
  args_dir[['script_dir']] <- args_line[[6]]
  args_dir[['job_tag']] <- args_line[[8]]
} 

setwd( paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ))
# load person IDs
id.dt <- data.table(read.csv(file.path(args_dir[['data_dir']],'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')))
id.dt <- subset(id.dt,select = c("pt_id"))
tmp <- data.table(read.csv(file.path(args_dir[['data_dir']],'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')))
tmp <- subset(tmp,select = c("pt_id"))
id.dt <- rbind(id.dt,tmp)
id.dt <- unique(id.dt)
setkey(id.dt,pt_id)
setnames(id.dt,colnames(id.dt),toupper(colnames(id.dt)))

# set seed and sample identifier at random
set.seed(42)
id.dt[,AID:=sample(1:nrow(id.dt),replace = F)]
id.dt[,AID:=formatC(AID, width=floor(log10(nrow(tmp))) +1, flag="0")]
id.dt[,AID:=paste0('AID',AID)]
write.csv(id.dt,file = 'important_anonymisation_keys_210119.csv')
Sys.chmod('important_anonymisation_keys_210119.csv',mode = '0444')

