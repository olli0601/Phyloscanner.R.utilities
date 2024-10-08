# rk.seq.make.distance
#'  calculate length of valid sequences, total match scores, percentage of matches of pair of sequences in one batch and in one window
#' it inputs the subsequence and batch info, and returns a data.table with distances

library(seqinr)
library(data.table)
# if(1)
# {
#   args_dir <- list()
#   indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
#   args_dir[['data_dir']] <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
#   args_dir[['out_dir']] <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/'
#   args_dir[['batchi']] <- 61
#   args_dir[['windowi']] <- 48
# }

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
