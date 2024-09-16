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
  # fix namesF
  rownames(fas) <- gsub("UG50406_PanRak7","UG504060_PanRak7",rownames(fas))
  fasn <- fas[rownames(fas) %in% gsub("(consensus).*", "\\1",basename(dsn$FORGLOBALALN_FASTA)),]
  # all(gsub("(consensus).*", "\\1",basename(dsn$FORGLOBALALN_FASTA)) == rownames(fasn))
  rownames(fasn) <- dsn$PANGEA_ID
  write.dna(fasn, file=outfile_new, format='fasta', colsep='', nbcol=-1)
  
  # check
  library(seqinr)
  df <- readRDS(infile_new)
  fas <- read.fasta(outfile_new)
  # unique(unlist(lapply(fas, unique)))
  
  # str(fasn)
  # tmp <- gsub("(consensus).*", "\\1",basename(dsn$FORGLOBALALN_FASTA))
  # tmp[!tmp %in% rownames(fasn)]
  # "UG504060_PanRak7-E1A1L1P1S1C1R1_consensus"
  # which(rownames(fas) == "UG504060_PanRak7-E1A1L1P1S1C1R1_consensus")
  # grep("UG504060",rownames(fas))
  # grep("UG504060",basename(ds$FORGLOBALALN_FASTA))
  # fas[5059,]
  # dsn[grep("UG50406", dsn$FORGLOBALALN_FASTA),]
  # 
  
  # # CHECK
  # # eg1 
  # sequence <- readLines(ds$FORGLOBALALN_FASTA[1])
  # sequence <- paste(sequence[2:length(sequence)], collapse="")
  # sequence <- strsplit(sequence,"")
  # unique(sequence[[1]])
  # ds$PANGEA_ID[1]
  # # "PG14-UG503444"
  # ds[1,]
  # 
  # sequence2 <- readLines(dconsensus[dconsensus$pangea_id == "PG14-UG503444",]$file_name[2])
  # sequence2 <- paste(sequence2[2:length(sequence2)], collapse="")
  # sequence2 <- strsplit(sequence2,"")
  # unique(sequence2[[1]])
  
  # # eg2
  sequence2 <- readLines(dconsensus$file_name[1])
  sequence2 <- paste(sequence2[2:length(sequence2)], collapse="")
  sequence2 <- strsplit(sequence2,"")
  unique(sequence2[[1]])
  dconsensus[1,] # PG14-UG500001
  
  sequence <- readLines(ds[ds$PANGEA_ID == "PG14-UG500001"]$FORGLOBALALN_FASTA)
  sequence <- paste(sequence[2:length(sequence)], collapse="")
  sequence <- strsplit(sequence,"")
  unique(sequence[[1]])
  
  # readLines(dconsensus[dconsensus$pangea_id == "PG14-UG503444",])
  
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

rk.seq.make.subsequence <- function()
{
  #' break sequences into subsequences by windows

  library(ape)
  library(data.table)
  
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
    cat(i,";")
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

rk.seq.make.scores <- function()
{
  #' count valid length, match scores and match percentages
  #' 
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

rk.seq.consensus.length <- function(){
  library(seqinr)
  library(data.table)
  windown <- 99
  indir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/"
  infiles <- file.path(indir, paste0('subsequence', 1:windown, '.fasta'))
  dl <- list()
  for (i in seq_len(length(infiles))){
    cat(i, "; ")
    sq <- read.fasta(infiles[i])
    dl[[i]] <-  sapply(sq, function(x)sum(x !='-' & x !='?' & x !='n'))
  }
  dl <- do.call(rbind, lapply(dl, function(x) x[match(names(dl[[1]]), names(x))]))
  
  length_cut <- 250
  vdl <- apply(dl, 2, function(x)sum(x>length_cut))
  length(vdl) # 7691
  sum(vdl!=0) # 7588
  sum(vdl >= 4) # 7364
}

rk.seq.submit.make.scores <- function(){
  
  #' write script to calculate scores
  
  library(ape)
  
  window_size <- 500
  batch_size <- 100
  
  # files
  #
  scriptdir <- '~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1521/'
  indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
  outdir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/'
  infile <- file.path(indir, '240809_PANGEA2_RCCS_final_samples_alignment_minusPoor1.fasta')
  source(paste0(scriptdir,"Rk1521_02_potential_nets.R"))
  
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
                          hpc.load= "module load anaconda3/personal\nsource activate phyloenv", 
                          hpc.array= jobn,
                          hpc.q = NA)
  
  cmd <- paste0(cmd ,'\n case $PBS_ARRAY_INDEX in \n')
  for (i in seq_len(jobn)) {
    cmd <- paste0(cmd, i,') \n')
    cmd <- paste0(cmd, 'Rscript ', '~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1521/distance.R -data_dir ',  outdir, 
                  ' -out_dir ', outdir,
                  ' -batchi ', batchs[i],
                  ' -windowi ', windows[i], '\n')
    cmd <- paste0(cmd, ';; \n')
  }
  cmd <- paste0(cmd, 'esac \n')
  
  # write submission file	
  processing.file <- file.path(outdir, 'processing.sh')
  cat(cmd, file=processing.file)
  
  # set permissions
  Sys.chmod(processing.file, mode='644')	
  
  # run
  cmd <-''
  cmd 	<- paste0(cmd, '\tcd ', dirname(processing.file),'\n')
  cmd 	<- paste0(cmd,'\tqsub ', basename(processing.file),'\n')
  cat(cmd)
  
}



rk.seq.gather.scores <- function(){

  #' gather scores together per chunk of sequences
  #' 
  library(data.table)
  
  indir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/'

  batch_size <- 100
  windown <- 99
  
  infile_base <- paste0(indir, paste0('subsequence',1:windown))
  
  # gather outputs
  for (w in 1:windown) {
    cat(w, "; ")
    base <- infile_base[w]
    df <- data.table()
    for (i in 1:batch_size) {
      filename <- paste0(base,'_batch',i,'.rds')
      # cat(filename)
      if(!file.exists(filename)){warning(paste0(filename, " does not exist!"))}
      tmp <- readRDS(filename)
      df <- rbind(df,tmp)
      # file.remove(paste0(base,'_batch',i,'.rds'))
    }
    
    # summaries
    dmeta <- data.table(WIN = w,
                        N = nrow(df),
                        N_VALID = nrow(df[LENGTH!=-1]),
                        LEN_GR_200 = nrow(df[LENGTH >= 200]),
                        LEN_GR_250 = nrow(df[LENGTH >= 250]),
                        LEN_GR_300 = nrow(df[LENGTH >= 300]),
                        MAX_LEN = max(df[LENGTH>=0,LENGTH]),
                        MIN_LEN = min(df[LENGTH>=0,LENGTH]),
                        MAX_SCORE = max(df[PERC>=0,PERC]),
                        MIN_SCORE = min(df[PERC>=0,PERC]))
    
    # remove invalid pairs
    df <- df[LENGTH >= 0]
    df[, WINDOW := w]
    
    # save combined results per window
    save(df,dmeta,file=paste0(base,'_results.rda'))
    gc()
  }
}

# rk.seq.depths <- function(){
#   
#   library(ape)
#   library(data.table)
#   library(ggplot2)
#   indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
#   infile_new <- file.path(indir, "240809_PANGEA2_RCCS_final_samples_minrun125_minnuc250_minusPoor1.rds")
#   infile_new_seq <- file.path(indir, '240809_PANGEA2_RCCS_final_samples_alignment_minusPoor1.fasta')
#   outdir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/'
#   outfile <- paste0(indir, "240809_PANGEA2_RCCS_final_samples_minrun125_minnuc250_minusPoor1_depths.rds")
#   
#   # load data
#   alignment <- read.dna(file = infile_new_seq, format = "fasta")
#   nsequence <- nrow(alignment)
#   npos <- ncol(alignment)
#   
#   # windows
#   window_size <- 500
#   windows_last_start <- ceiling(npos/100) * 100 - window_size + 1
#   windows_2last_end <- floor(npos/100) * 100
#   windows_start <- seq(1,windows_last_start,100)
#   windows_end <- c(seq(window_size ,windows_2last_end,100),npos)
#   windows <- seq_len(length(windows_start))
#   dw <- data.table(WINDOW=windows,
#                    START=windows_start,
#                    END=windows_end)
#   
#   #
#   dsn <- readRDS(infile_new)
#   ddepth <- data.table()
#   for(i in 1:nrow(dsn)){
#     if(i%%100 == 0){cat("processing sequence ", i, '\n')}
#   
#     sequence <- readLines(dsn$FORGLOBALALN_FASTA[i])
#     sequence <- paste(sequence[2:length(sequence)], collapse="")
#     sequence <- strsplit(sequence,"")
#     
#     tmp_df <- dw[,{
#       tmp = sequence[[1]][START:END]
#       tmp = tmp[!grepl("[^A-Za-z]",tmp)]
#       list(
#         LEN_LETTER = length(tmp),
#         HIGH_DEPTH = sum(tmp == toupper(tmp)),
#         PANGEA_ID = dsn$PANGEA_ID[i],
#         FILE = dsn$FORGLOBALALN_FASTA[i])
#     }, by = c('WINDOW', 'START', 'END')]
#     ddepth <- rbind(ddepth, tmp_df)
#   }
#   
#   #
#   # ddepth[ddepth$PANGEA_ID == "PG15-UG504899" & ddepth$WIN == 1,]
#   saveRDS(ddepth, file = outfile)
#   
#   #
#   ddepth <- readRDS(outfile)
#   # ddepth2 <- readRDS(file.path(indir, "240809_PANGEA2_RCCS_depths.rds"))
#   # ddepth[ddepth$WINDOW == 1 & ddepth$PANGEA_ID == "PG19-UG001819",]
#   # ddepth2[ddepth2$WINDOW == 1 & ddepth2$PANGEA_ID == "PG19-UG001819",]
#   
#   length_cut <- 250
#   ddepth_select <- subset(ddepth,select=c('WINDOW','LEN_LETTER','HIGH_DEPTH','PANGEA_ID'))
#   
#   ddepth_select <- ddepth_select[HIGH_DEPTH >= length_cut, ]
#   
#   
#   ddepth_select_w <- ddepth_select[,list(PANGEA_ID=unique(PANGEA_ID)),
#                                    by=c('WINDOW')]
#   ddepth_select_w[,HD:=1]
#   ddepth_select_w <- data.table::dcast(ddepth_select_w, PANGEA_ID~WINDOW, value.var='HD')
#  
#   length(unique(ddepth$PANGEA_ID)) # 7691
#   nrow(ddepth_select_w) # 7028
#   sum(apply(as.matrix(ddepth_select_w[,-1]), 1, function(x)sum(x,na.rm=T)) >=4)
#   # 6807
#   
#   # plot histogram of number of windows with high depth
#   #
#   tmp = data.table(V=rowSums(!is.na(ddepth_select_w[,2:100])))
#   ggplot(tmp,aes(V)) +
#     geom_histogram() +
#     labs(x=paste0('number of windows with high depth postions >= ',length_cut))+
#     theme_bw()
#   ggsave(filename = file.path(outdir, 'histogram_number_of_window_depth_gt250.pdf'),
#          width = 6, height = 4)
#  #  
#  #  
#  #  ddepth <- readRDS(outfile)
#  #  length_cut <- 250
#  #  ddepth_select <- subset(ddepth,select=c('WINDOW','LEN_LETTER','HIGH_DEPTH','PANGEA_ID'))
#  #  ddepth_select <- ddepth_select[LEN_LETTER >= length_cut, ]
#  #  ddepth_select_w <- ddepth_select[,list(PANGEA_ID=unique(PANGEA_ID)), by=c('WINDOW')]
#  #  ddepth_select_w[,HL:=1]
#  #  ddepth_select_w <- data.table::dcast(ddepth_select_w, PANGEA_ID~WINDOW, value.var='HL')
#  #  length(unique(ddepth$PANGEA_ID)) # 7691
#  #  nrow(ddepth_select_w) # 7599
#  #  sum(apply(as.matrix(ddepth_select_w[,-1]), 1, function(x)sum(x,na.rm=T)) >=4) # 7386
#  #  
#  # #
#  # dcheck <- read.csv("240809_PANGEA2_RCCS_final_samples_minrun125_minnuc250_all_seqs_of_ind_poor.csv")
#  # dcheck <- ddepth[ddepth$PANGEA_ID %in% dcheck$PANGEA_ID,] # 12
#  # dcheck[,sum(HIGH_DEPTH!=0), by ="PANGEA_ID"] # 5
#  # dcheck[,sum(HIGH_DEPTH>length_cut), by ="PANGEA_ID"] # 4
# }



rk.seq.couple.scores <- function(){
  #' find thresholds from couples 

  library(ape)
  library(data.table)
  
  window_size <- 500
  windown <- 99
  
  indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
  outdir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/'
  infiles <- paste0(outdir, 'subsequence',1:windown,'_results.rda')
  infile_sequence <- paste0(indir, "240809_PANGEA2_RCCS_final_samples_alignment_minusPoor1.fasta")
  infile_db <- paste0(indir, "240809_pangea_db_sharing_extract_rakai_combined.rds")
  # couples from last round
  infile_couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'

  # load data
  alignment <- read.dna(file = infile_sequence, format = "fasta")
  dinfo <- data.table(PANGEA_ID=rownames(alignment))
  id.dt <- data.table(readRDS(infile_db))
  id.dt <- subset(id.dt,select = c("PANGEA_ID","PT_ID"))
  id.dt <- unique(id.dt)
  dinfo <- merge(dinfo, id.dt, by="PANGEA_ID", all.x=T)
  # sum(is.na(dinfo$PT_ID)) # 0
  
  # couples
  load(infile_couple)
  couple <- unique(subset(data.table(coupdat), select=c('male.RCCS_studyid', 'female.RCCS_studyid')))
  couple[, COUPLE:=1]
  couple[,male.RCCS_studyid:=paste0('RK-',male.RCCS_studyid)]
  couple[,female.RCCS_studyid:=paste0('RK-',female.RCCS_studyid)]
  
  # 
  tmp_couple <- subset(couple,select=c('male.RCCS_studyid','female.RCCS_studyid', 'COUPLE'))
  tmp <- copy(tmp_couple)
  setnames(tmp_couple, c('male.RCCS_studyid','female.RCCS_studyid'),c('PT_ID1','PT_ID2'))
  setnames(tmp, c('male.RCCS_studyid','female.RCCS_studyid'),c('PT_ID2','PT_ID1'))
  tmp_couple[,sex1_couple := 'M']
  tmp_couple[,sex2_couple := 'F']
  tmp[,sex1_couple := 'F']
  tmp[,sex2_couple := 'M']
  tmp_couple <- rbind(tmp_couple, tmp)
  unique(tmp_couple)
  tmp_couple  <- tmp_couple[!is.na(PT_ID1) & !is.na(PT_ID2) & PT_ID1!=PT_ID2,]
  
  length_cutoff <- 250
  
  # add scores
  #
  for (i in 1:windown) {
    
    cat('\n process window ',i , '\n')
    load(infiles[i])
    
    # 
    df <- df[LENGTH >= length_cutoff,]
    
    # add couple
    setnames(dinfo, colnames(dinfo), paste0(colnames(dinfo),'1'))
    df <- merge(df, dinfo, by.x=c('TAXA1'), by.y=c('PANGEA_ID1'), all.x=T) 
    setnames(dinfo, colnames(dinfo), gsub('1','2',colnames(dinfo)))
    df <- merge(df, dinfo, by.x=c('TAXA2'), by.y=c('PANGEA_ID2'), all.x=T) 
    setnames(dinfo, colnames(dinfo), gsub('2','',colnames(dinfo)))
    
    # keep couples only
    df <- df[!is.na(PT_ID1) & !is.na(PT_ID2) & PT_ID1!=PT_ID2,]
    df <- merge(df, tmp_couple, by=c('PT_ID1','PT_ID2'), all.x=T)
    df <- df[COUPLE==1, c('PT_ID1','PT_ID2','PERC','WINDOW')]
    if(i==1){
      ans <- df
    }else{
      ans <- rbind(ans, df)
    }
  }
  
  # save 
  saveRDS(ans, file = file.path(outdir,'scores_couples.rds'))  

  infile_new <- file.path(indir, "240809_PANGEA2_RCCS_final_samples_minrun125_minnuc250_minusPoor1.rds")
  dpr <- readRDS(infile_new)
  dpr$PT_ID
}



rk.seq.find.thresholds <- function(){

  library(rstan)
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  
  dir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/'
  dir_plot <- paste0(dir, "plots/")
  dir.create(dir_plot)
  infile <- paste0(dir, "scores_couples.rds")
  ans <- readRDS(infile)
  
  similarity <- list()
  for (i in seq_len(max(ans$WINDOW))) {
    similarity[[i]] <- ans[WINDOW==i,]$PERC
  }
  
  max_data <- max(lengths(similarity))
  
  for (i in seq_len(max(ans$WINDOW))) {
    if(length(similarity[[i]]) !=max_data){
      similarity[[i]] <- c(similarity[[i]], rep(NA,max_data-length(similarity[[i]])))
    }
  }
  
  stan_data <- list()
  stan_data$similarity <- do.call(rbind,similarity)
  stan_data$Nsimilarity <- apply(stan_data$similarity,1,function(x)sum(!is.na(x)))
  stan_data$Nwindow <- length(similarity)
  stan_data$max_data <- max_data
  stan_data$similarity[is.na(stan_data$similarity)] <- -1
  
   
  # tmp <- stan_data$similarity
  # tmp <- data.table(t(tmp[seq(1,99,20),]))
  # colnames(tmp) <- paste0(windows_start[seq(1,99,20)],' - ', windows_end[seq(1,99,20)])
  tmp <- data.table(t(stan_data$similarity))
  tmp <- melt(tmp)
  tmp <- tmp[value!=-1]
  
  g <- ggplot(tmp, aes(value))+
    geom_histogram()+
    facet_wrap(~gsub("V","window ",variable), ncol = 10, scales = 'free') +
    theme_bw()+
    labs(x='\n similarity scores', y='count \n') +
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0.05))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          text = element_text(size=20))
  ggsave(paste0(dir_plot, "couples_scores_histogram.pdf"), g, width = 40, height = 30, limitsize = F)
  
  
  # fit bimodal model to scores
  #
  stan_data$similarity <- log(stan_data$similarity)
  stan_data$similarity[is.na(stan_data$similarity)] <- -1
  stan_code <- "
  data {
  int<lower = 0> N;
  vector[N] y;
  }

  parameters {
    ordered[2] mu;
    real<lower=0> sigma[2];
    real<lower=0, upper=1> theta;
  }

  model {
  sigma ~ exponential(5.0);
  mu ~ normal(0, 1);
  theta ~ beta(5, 5);
  for (n in 1:N)
     target += log_mix(theta,
                       normal_lpdf(y[n] | mu[1], sigma[1]),
                       normal_lpdf(y[n] | mu[2], sigma[2]));
  }

  "
  fit <- list()
  
  windown <- 99
  for (i in 1:windown) {
    cat('fit window ', i, '\n')
    tmp <- list()
    tmp$y <- (stan_data$similarity[i,1:stan_data$Nsimilarity[i]])
    tmp$N <- stan_data$Nsimilarity[i]
    fit[[i]] <- stan(model_code=stan_code,
                     data=tmp,
                     chains=1, seed=42, iteration = 2e5)
    
  }
  
  # 
  save(fit,file=file.path(dir,'couples_scores_stan_fit.rda'))
  
  
  load(file.path(dir,'couples_scores_stan_fit.rda'))
  
  # check fit 
  #
  rh= c()
  ness= c()
  tess= c()
  bess= c()
  for (i in 1:windown) {
    cat(i, "; ")
    pars <- rstan::extract(fit[[i]], pars=names(fit[[i]])[!grepl('lp_',names(fit[[i]]))])
    pars = do.call(cbind,pars)
    
    rh = c(rh,max(summary(fit[[i]])$summary[, "Rhat"]))
    ness = c(ness,min(summary(fit[[i]])$summary[, "n_eff"]))
    tess = c(tess,min(apply(pars, 2, ess_tail)))
    bess = c(bess,min(apply(pars, 2, ess_bulk)))
  }
  
  cat('range of rhat is ',range(rh),'\n')
  cat('range of ess is ',range(ness),'\n')
  cat('range of ess_tail is ',range(tess),'\n')
  cat('range of ess_bulk is ',range(bess),'\n')
  
  
  # extract outputs 
  #
  df = data.table()
  for (i in 1:windown) {
    cat(i, "; ")
    tmp = data.table(summary(fit[[i]])$summary[1:5,c("2.5%","50%","97.5%")])
    tmp[,VAR:=rownames(summary(fit[[i]])$summary)[1:5]]
    tmp[,WINDOW:=i]
    df = rbind(df, tmp)
  }
  
  # visualise
  #
  ans[,PERC:=log(PERC)]
  #  range(ans$PERC)
  #  -0.4453817  0.0000000
  
  tmp <- dcast(df, WINDOW~VAR, value.var='50%')
  tmp <- merge(tmp, ans[,list(NOBS=length(PERC)),by=c('WINDOW')],by='WINDOW')
  ans[,WINDOW:=as.integer(WINDOW)]
  bw <- 0.45/30
  tmp <- tmp[,{
    logsim=seq(-0.45,0,length.out=101)
    y1=bw * NOBS * theta * dnorm(logsim,`mu[1]`,`sigma[1]`)
    y2=bw * NOBS * (1-theta) * dnorm(logsim,`mu[2]`,`sigma[2]`)
    list(X=logsim,Y1=y1,Y2=y2)
  },by='WINDOW']
  
  
  g <- ggplot(ans, aes(PERC))+
    geom_histogram(binwidth = bw)+
    geom_line(tmp,mapping=aes(x=X,y=Y1))+
    geom_line(tmp,mapping=aes(x=X,y=Y2))+
    facet_wrap(~WINDOW, ncol = 10, scales = 'free') +
    theme_bw()+
    labs(x='log similarity scores') +
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0.05))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  ggsave(paste0(dir_plot, "couples_scores_fitted_histogram.pdf"), g, width = 40, height = 30, limitsize = F)
  
  # mean and var vs windows
  #
  tmp = df[grep('mu',VAR),]
  tmp = dcast(tmp, WINDOW~VAR,value.var = '50%')
  tmp[,diff:=`mu[2]`-`mu[1]`]

  ggplot(df[grep('mu',VAR)],aes(x=WINDOW))+
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
    geom_line(aes(y = `50%`))+
    facet_wrap(VAR~.,ncol = 1) +
    labs(x='window', y='means')
  ggsave(paste0(dir_plot, "component_means_vs_windows.pdf"), 
         width = 6, height= 8)
  
  ggplot(df[grep('sigma',VAR)],aes(x=WINDOW))+
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
    geom_line(aes(y = `50%`))+
    facet_wrap(VAR~.,ncol = 1) +
    labs(x='window', y='standard deviations')
  ggsave(paste0(dir_plot, "component_sds_vs_windows.pdf"), 
         width = 6, height= 8)
  
  
  # threshold: 5% percentile 
  #
  for (i in 1:windown) {
    cat(i, '; ')
    if(i==1)threshold = list()
    tmp <- c()
    pars <- rstan::extract(fit[[i]], pars=c('mu','sigma','theta'))
    for (j in 1:nrow(pars$mu)) {
      tmp <- c(tmp, qlnorm(0.05, pars$mu[j,2], pars$sigma[j,2]))
      
    }
    threshold[[i]] <- tmp
  }
  
  
  threshold <- lapply(threshold,function(x){quantile(x,probs = c(0.025,0.5,0.975))})
  threshold <- data.table(do.call(rbind, threshold))
  saveRDS(threshold, 
       file = file.path(dir, paste0('couples_scores_threshold_5percentile.rds')))
  
  # plot threshold
  #
  threshold <- readRDS(file.path(dir, paste0('couples_scores_threshold_5percentile.rds')))
  tmp <- data.table(threshold)
  tmp[, WINDOW:=seq_len(nrow(tmp))]
  

  ggplot(tmp,aes(x=WINDOW,y=`50%`,ymin=`2.5%`,ymax=`97.5%`))+
    geom_line()+
    geom_ribbon(alpha=0.6)+
    theme_bw()+
    labs(x='\n genomic windows',y='threshold \n')+
    scale_x_continuous(expand = c(0,0),breaks = seq(0,100,by=10)) +
    scale_y_continuous(expand = c(0,0), labels = scales::percent, limits = c(0.85,1))+
    theme(text = element_text(size=20))
  
  ggsave(paste0(dir_plot, "thresholds_vs_windows.pdf"), 
         width = 6, height= 4)
  
}



rk.seq.classify.pairs <- function(){
  
  library(data.table)
  library(ggplot2)
  
  window_size <- 500
  windown <- 99
  
  indir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/'
  infiles <- paste0(indir, 'subsequence',1:windown,'_results.rda')
  infile_threshold <- paste0(indir, 'couples_scores_threshold_5percentile.rds')
  # infile_depths <- "/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/240809_PANGEA2_RCCS_final_samples_minrun125_minnuc250_minusPoor1_depths.rds"
 
  #
  threshold <- readRDS(infile_threshold)
  threshold <- threshold$`50%`
  length_cut <- 250
  
  # #
  # ddepth <- readRDS(infile_depths)
  # ddepth_select <- subset(ddepth,select=c('WINDOW','LEN_LETTER','HIGH_DEPTH','PANGEA_ID'))
  # ddepth_select <- ddepth_select[HIGH_DEPTH >= length_cut, ]
  
  #
  for (i in 1:windown) {
    # 
    cat(i, "; ")
    load(infiles[i])
    df <- df[LENGTH >= length_cut,]
    
    # # select high depths
    # #
    # ddepth_selecti <- unique(subset(ddepth_select[WINDOW==i, ], select='PANGEA_ID'))
    # ddepth_selecti[,HD:=1]
    # setnames(ddepth_selecti, colnames(ddepth_selecti), paste0(colnames(ddepth_selecti),'1'))
    # df <- merge(df, ddepth_selecti, by.x='TAXA1', by.y='PANGEA_ID1', all.x=T)
    # setnames(ddepth_selecti, colnames(ddepth_selecti), gsub('1','2',colnames(ddepth_selecti)))
    # df <- merge(df, ddepth_selecti, by.x='TAXA1', by.y='PANGEA_ID2', all.x=T)
    # setnames(ddepth_selecti, colnames(ddepth_selecti), gsub('2','',colnames(ddepth_selecti)))
    # df <- df[HD1==1 & HD2==1,]
    
    # add threshold
    #
    df[, threshold_bimodal:=threshold[i]]
    df[, CLOSE:= as.integer(PERC > threshold_bimodal)]
    setnames(df, 'PERC', paste0('PERC',i))
    setnames(df, 'CLOSE', paste0('CLOSE',i))
    
    # length(unique(paste0(df$TAXA1, df$TAXA2)))
    # nrow(df)
  
    if(i ==1){
      dscore <- subset(df, select = c('TAXA1','TAXA2',paste0("PERC",i)))
      dclassif <- subset(df, select = c('TAXA1','TAXA2',paste0("CLOSE",i)))
    }else{
      dscore <- merge(dscore, subset(df, select = c('TAXA1','TAXA2',paste0("PERC",i))), by = c('TAXA1','TAXA2'), all=T)
      dclassif <- merge(dclassif, subset(df, select = c('TAXA1','TAXA2',paste0("CLOSE",i))), by = c('TAXA1','TAXA2'), all=T)
    }
    gc()
  }
  saveRDS(dscore, paste0(indir, "scores.rds"))
  saveRDS(dclassif, paste0(indir, "classification.rds"))
  
  # summarise over windows
  dclassif <- readRDS(paste0(indir, "classification.rds"))
  dclassif_nna <- rowSums(is.na(dclassif[,3:101]))
  dclassif_nnna <- windown - dclassif_nna
  dclassif_nc <- rowSums(dclassif[,3:101], na.rm = T)
  # 
  # at least one common windows 
  #
  ds <- dclassif[,1:2]
  ds$N_ELG <- dclassif_nnna
  ds$N_CLOSE <- dclassif_nc
  # hist(ds$N_ELG)
  window_cut <- 1
  dclassifye <- ds[N_ELG >= window_cut,]
  saveRDS(dclassifye, paste0(indir, "classification_elg_",window_cut,".rds"))
  print(length(unique(c(dclassifye$TAXA1, dclassifye$TAXA2)))) # 7588
  
  #
  dclassifye <- readRDS(paste0(indir, "classification_elg_",window_cut,".rds"))
  dclassifye[, P_CLOSE := N_CLOSE / N_ELG]
  hist(dclassifye$P_CLOSE)
  ggplot(dclassifye,aes(P_CLOSE))+
    geom_histogram()+
    theme_bw()+
    labs(x='\n proportions of close windows','counts of sequence pairs \n ')+
    scale_x_continuous(expand = c(0,0), labels = scales::percent) +
    scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(ylim=c(0,5e4))+
    theme(text = element_text(size=20))
  ggsave(paste0(indir, "plots/histogram_pclose_",window_cut,".pdf"),
         width = 6, height = 4)
  
  dclassifye[,P_CLOSE_C := cut(P_CLOSE, breaks = seq(0,1,by=0.2), include.lowest = T)]
  g <- ggplot(dclassifye, aes(N_ELG, fill = P_CLOSE_C))+
    geom_histogram(bins = 100) +
    scale_x_continuous() +
    labs(x = "number of eligible windows", fill = "proportion of close windows")+
    theme_bw() + theme(legend.position = "bottom")
  ggsave(paste0(indir, "plots/histogram_nelg_by_pclose.pdf"), g,
         width = 6, height = 4)
  
  g <- g + coord_cartesian(ylim=c(0,5e3)) 
  ggsave(paste0(indir, "plots/histogram_nelg_by_pclose_zoomin.pdf"), g,
         width = 6, height = 4)
  #
  dt <- readRDS(paste0(indir, "couples_scores_threshold_5percentile.rds"))
  dt[, WINDOW := seq_len(nrow(dt))]
  dcp <- readRDS(paste0(indir, "scores_couples.rds"))
  dcp <- merge(dcp, subset(dt, select = c("WINDOW","50%")), by ="WINDOW", all.x = T)
  dcp <- dcp[, list(P_CLOSE = sum(PERC >= `50%`)/length(PERC)), by = c("PT_ID1","PT_ID2")]
  ggplot(dcp,aes(P_CLOSE))+
    geom_histogram()+
    theme_bw()+
    labs(x='\n proportions of close windows','counts of sequence pairs \n ')+
    scale_x_continuous(expand = c(0,0), labels = scales::percent) +
    scale_y_continuous(expand = c(0,0)) +
    theme(text = element_text(size=20))
  ggsave(paste0(indir, "plots/histogram_pclose_",window_cut,"_couple.pdf"),
         width = 6, height = 4)
  
  #
  cut_pclose <- 0.8
  dc <- dclassifye[P_CLOSE >= cut_pclose, ]
  length(unique(c(dc$TAXA1, dc$TAXA2))) # 4836; 6974
  saveRDS(dc, paste0(indir, "closepairs_",window_cut,".rds"))
  
  # dh <- dcast(dclassifye, TAXA1 ~ TAXA2, value.var = "P_CLOSE")
  # dh <- as.matrix(dh[,-1])
  # dh[is.na(dh)] <- -1
  # # image(dh,useRaster=TRUE)
  # library(gplots)
  # pdf(paste0(indir, "heatmap_pclose.pdf"), width = 20, height = 20)
  # heatmap.2(dh)
  # dev.off()  
}



rk.seq.make.scores <- function()
{
  #' count valid length, match scores and match percentages
  #' 
  library(seqinr)
  library(data.table)
  
  args_line <-  as.list(commandArgs(trailingOnly=TRUE))
  print(args_line)
  if(length(args_line) > 0) 
  {
    stopifnot(args_line[[1]]=='-jobi')
    i <- as.integer(args_line[[2]])
  } 
  
  #
  indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
  outdir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/'
  infile <- paste0(indir, "240809_PANGEA2_RCCS_final_samples_alignment_minusPoor1.fasta")
  
  #	load
  sq <- read.fasta(infile)
  
  # remove sequence with all n
  # sq	<- sq[ dfo[, unique(c(TAXA1, TAXA2))] ]
  # sq.n <- unlist(lapply(sq, function(x){sum(x=='?' | x=='-' | x=='n')==length(x)}))
  # sq.n.names <- names(sq.n[sq.n==TRUE])
  # tmp <- dfo[TAXA1 %in% sq.n.names|TAXA2 %in% sq.n.names]
  # tmp[,LENGTH:=-1.0]
  # tmp[,MATCHES:=-1.0]
  # tmp[,PERC:=-1.0]
  # dfo <- dfo[!TAXA1 %in% sq.n.names & !TAXA2 %in% sq.n.names]
  
  cat(i,"\n")
  dfo <- data.table(TAXA1 = names(sq)[i],
                    TAXA2 = names(sq)[-i])
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
  }, by=c('TAXA1','TAXA2')]
  dfo[,PERC:=MATCHES/LENGTH]
  saveRDS(dfo, file = paste0(outdir, "avescores_",names(sq)[i], ".rds"))
  
  
  
  
  # make scripts
  jobn <- length(sq)
  cmd <-  make.PBS.header(hpc.walltime=24, 
                          hpc.select=1, 
                          hpc.nproc=1, 
                          hpc.mem= "20gb", 
                          hpc.load= "module load anaconda3/personal\nsource activate phyloenv", 
                          hpc.array= jobn,
                          hpc.q = NA)
  
  cmd <- paste0(cmd ,'\n case $PBS_ARRAY_INDEX in \n')
  for (i in seq_len(jobn)) {
    cmd <- paste0(cmd, i,') \n')
    cmd <- paste0(cmd, 'Rscript ', '~/ave_score.R -jobi ',  i, 
                  '\n')
    cmd <- paste0(cmd, ';; \n')
  }
  cmd <- paste0(cmd, 'esac \n')
  
  # write submission file	
  processing.file <- file.path(outdir, 'processing_avescore.sh')
  cat(cmd, file=processing.file)
  
  # set permissions
  Sys.chmod(processing.file, mode='644')	
  
  # run
  cmd <-''
  cmd 	<- paste0(cmd, '\tcd ', dirname(processing.file),'\n')
  cmd 	<- paste0(cmd,'\tqsub ', basename(processing.file),'\n')
  cat(cmd)
}

rk.seq.cluster <- function(){
  library(data.table)
  library(igraph)
  library(seqinr)
  
  #
  indir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/'
  datadir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
  infile <- file.path(datadir, "240809_PANGEA2_RCCS_final_samples_minrun125_minnuc250_minusPoor1.rds")
  infile_seq <- file.path(datadir, '240809_PANGEA2_RCCS_final_samples_alignment_minusPoor1.fasta')
  window_cut <- 1
  cpfile <- paste0(indir, "closepairs_",window_cut,".rds")
  dc <- readRDS(cpfile)
  dt <- readRDS(infile)

  # compute number of neighbors
  tmp <- copy(dc)
  setnames(tmp, c("TAXA1","TAXA2"), c("TAXA2","TAXA1"))
  dc <- rbind(dc, tmp)
  dcs <- dc[,list(N_IND_CLOSE = length(TAXA2)), by ="TAXA1"]
  setkey(dcs, N_IND_CLOSE)
  cat("Number of neighbours:\n")
  print(hist(dcs$N_IND_CLOSE, plot = F, breaks = seq(0, 100 * ceiling(max(dcs$N_IND_CLOSE)/100) , by = 10)))
  cat("Sequences with >=50 neighbours: \n")
  print(dcs[dcs$N_IND_CLOSE > 50, ])
  # hist(dcs$N_IND_CLOSE, plot=F, breaks = 50)
  # for(taxa in dcs[dcs$N_IND_CLOSE > 50, ]$TAXA){
  #   print(dc[TAXA1==taxa,])
  # }
  cat("Unique sequences in potential networks: ", length(unique(c(dc$TAXA1, dc$TAXA2))), "\n")
  # print(length(unique(c(dc$TAXA1, dc$TAXA2)))) # 6974 (1) # 4836 (4)
 
  # remove pairs from the same individual
  dts <- subset(dt, select = c("PT_ID","PANGEA_ID"))
  dc <- merge(dc, dts, by.x = "TAXA1", by.y = "PANGEA_ID", all.x = T)
  dc <- merge(dc, dts, by.x = "TAXA2", by.y = "PANGEA_ID", all.x = T)
  setnames(dc, c("PT_ID.x","PT_ID.y"), c("PT_ID1", "PT_ID2"))
  # tmp <- dc[!(dc$TAXA1 %in%  dcs[dcs$N_IND_CLOSE > 50, ]$TAXA1) &
  #           !(dc$TAXA2 %in%  dcs[dcs$N_IND_CLOSE > 50, ]$TAXA1)]
  # tmp <- tmp[PT_ID1 != PT_ID2,]
  # length(unique(c(tmp$PT_ID1, tmp$PT_ID2)))# 2405
  # 
  dc <- unique(dc[PT_ID1!=PT_ID2,c("PT_ID1", "PT_ID2")])
  cat("Unique individuals in potential networks after removing pairs from the same individual: ",
  length(unique(c(dc$PT_ID1, dc$PT_ID2)))) # 5094 (1) # 2637 (4)
  cat("Number of close pairs: ", nrow(dc),"\n")
  
  # make graph
  chains <- graph.data.frame(dc[,c("PT_ID1","PT_ID2")], directed=FALSE, vertices=NULL)
  print(length(V(chains)))
  chains2 <- delete_vertices(chains, dts[dts$PANGEA_ID %in% dcs[dcs$N_IND_CLOSE > 50, ]$TAXA,]$PT_ID)
  print(length(V(chains2)))
  print(length(unique(dts[dts$PANGEA_ID %in% dcs[dcs$N_IND_CLOSE > 50, ]$TAXA,]$PT_ID)))
  
  rtc	<- data.table(ID=V(chains2)$name, CLU=clusters(chains2, mode="weak")$membership)
  
  # make clusters 
  tmp	<- rtc[, list(CLU_SIZE=length(ID)), by='CLU']
  setkey(tmp, CLU_SIZE)
  print(tail(tmp))
  tmp[, IDCLU:=rev(seq_len(nrow(tmp)))]
  rtc	<- subset( merge(rtc, tmp, by='CLU'))
  rtc[, CLU:=NULL]
  setkey(rtc, IDCLU)
  print(table(tmp$CLU_SIZE))
  
  
  # take the largest ID = 1, size = 835; 
  # but after removing > 50 neighbors, max size = 303; 
  # slightly larger but comparable than last round
  # remember to add them to not close individuals
  ind <- rtc[IDCLU==1,]$ID
  dcl <- subset(dc, select=c(PT_ID1, PT_ID2))
  dcl <- dcl[PT_ID1 %in% ind & PT_ID2 %in% ind,]
  lcons <- read.fasta(infile_seq)
  names(lcons) <- factor(names(lcons), dts$PANGEA_ID, dts$PT_ID)
  
  # similarity scores over whole genome
  dcl_sc <- dcl[,{
    seq1 <- lcons[[as.character(PT_ID1)]]
    seq2 <- lcons[[as.character(PT_ID2)]]
    # cat("seq1: ", length(seq1), "seq2: ", length(seq2), "\n")
    pos <- which(seq1 !='-' & seq1 !='?' & seq1 !='n' & seq2 !='-' & seq2 !='?' &  seq2 !='n')
    if (length(pos)==0 ){
      len <- as.integer(-1.0)
      match <- -1.0
    }else{
      seq1 <- seq1[pos]
      seq2 <- seq2[pos]
      len <- length(pos)
      match <- sum(sapply(1:length(pos),function(pos){is.match(seq1[pos],seq2[pos])}))
    }
    list(LENGTH = len,
         MATCHES = match)
  }, by=c('PT_ID1','PT_ID2')]

  # save(dcl_sc, rtc, file = paste0("~/cluster_",window_cut,".rda")) # need igraph 2.0.3 - cannot solve envir on hpc
  # load("~/cluster_1.rda")

  # set similar scores to 1
  #
  dcl_sc[, PERC := MATCHES/LENGTH]
  dcl_sc[PERC >= 0.975, PERC := 1]
  dcl_sc <- dcl_sc[,list(PERC = mean(PERC)), by=c('PT_ID1','PT_ID2')]
  chains <- graph.data.frame(dcl_sc, directed = FALSE, vertices = NULL)
  E(chains)$weight <- dcl_sc$PERC
  chains <- simplify(chains)
  E(chains)$weight <- ifelse(E(chains)$weight > 1,E(chains)$weight/2, E(chains)$weight)
  
  # 
  net <- chains
  V(net)$size <- 3
  V(net)$frame.color <- "white"
  V(net)$color <- "orange"
  V(net)$label <- ""
  E(net)$arrow.mode <- 0
  
  pdf(paste0(indir, "plots/largest_cluster"),width = 10, height = 8)
  plot(net)
  dev.off()
  
  # break
  # dts <- readRDS("~/dts.rds")
  set.seed(42)
  comm <- cluster_louvain(chains, resolution = 2.5) # , resolution = 1
  # layout <-layout.fruchterman.reingold(chains)
  # # before max = 50 but now max = 140
  # # but luckily resolution options allows further breaking
  # # resolution = 7 is determined by trying and making sure max is approx 50
  # # repeating this strategy leads to 391 bridging individuals, 
  # # between communities of bridge individuals, 317 individuals - hard to break
  # so try to remove weird inds with lots of neighbors when making clusters
  sort(sapply(1:max(comm$membership),function(x)sum(comm$membership==x)))
  plot(comm, chains, vertex.label=NA, vertex.size = 2)
  
  # df <- data.table(PT_ID=comm$names, MEMBERSHIP = comm$membership)
  # # dts[dts$PANGEA_ID %in% dcs[dcs$N_IND_CLOSE > 30, ]$TAXA,]$PT_ID %in% attr(V(chains), "names")
  # for(m in df[PT_ID %in% dts[dts$PANGEA_ID %in% dcs[dcs$N_IND_CLOSE > 30, ]$TAXA,]$PT_ID, ]$MEMBERSHIP){
  #   cat(m, " : ", sum(comm$membership == m), "\n")
  # }
  # #
  # # option 2
  # dts[dts$PANGEA_ID %in% dcs[dcs$N_IND_CLOSE > 30, ]$TAXA,]$PT_ID %in% names(V(chains))
  # chains2 <- delete_vertices(chains, dts[dts$PANGEA_ID %in% dcs[dcs$N_IND_CLOSE > 30, ]$TAXA,]$PT_ID)
  # tmp	<- data.table(ID=V(chains2)$name, CLU=clusters(chains2, mode="weak")$membership)
  # table(tmp$CLU)
  # 
  # comm <- cluster_louvain(chains2)
  # sort(sapply(1:max(comm$membership),function(x)sum(comm$membership==x)))
  df <- data.table(PT_ID=comm$names, MEMBERSHIP = comm$membership)
  # 
  
  #
  df_graph <- as_long_data_frame(chains)
  df_graph <- data.table(df_graph[,3:5])
  colnames(df_graph) <- c("PERC","PT_ID1","PT_ID2")

  # bridging individuals function
  take_bridging_individual <- function(df_graph, df){
    bridging <- list()
    for (i in 1:max(df$MEMBERSHIP)) {
      tmp1 <- df[MEMBERSHIP==i,]$PT_ID
      tmp2 <- df_graph[PT_ID1 %in% tmp1 & !(PT_ID2 %in% tmp1) & PERC >=0.975,]
      tmp3 <- df_graph[PT_ID2 %in% tmp1 & !(PT_ID1 %in% tmp1) & PERC >=0.975,]
      bridging[[i]] <- unique(c(as.character(tmp2$PT_ID1),
                               as.character(tmp3$PT_ID2)))
    }
    return(bridging)
  }
  
  #
  bridging <- take_bridging_individual(df_graph, df)
  bridging_all <- unique(do.call(c,bridging))
  length(bridging_all)
  # 50
  # might be okay
  
  # # break bridging individuals into clusters
  # id_clu <- max(df$MEMBERSHIP)
  # dcl_sc <- dcl_sc[PT_ID1 %in% bridging_all & PT_ID2 %in% bridging_all]
  # chains <- graph.data.frame(dcl_sc, directed=FALSE, vertices=NULL)
  # E(chains)$weight <- dcl_sc$PERC
  # chains <- simplify(chains)
  # E(chains)$weight <- ifelse(E(chains)$weight > 1,E(chains)$weight/2, E(chains)$weight)
  # tmp <- clusters(chains, mode='weak')
  # # cannot break here
  # 
  # # break bridging individuals into communities
  # PT_ID_bridge <- names(tmp$membership)
  # dcl_sc_bridge <- dcl_sc[PT_ID1 %in% PT_ID_bridge & PT_ID2 %in% PT_ID_bridge]
  # chains <- graph.data.frame(dcl_sc_bridge, directed=FALSE, vertices=NULL)
  # E(chains)$weight <- dcl_sc_bridge$PERC
  # chains <- simplify(chains)
  # E(chains)$weight <- ifelse(E(chains)$weight > 1,E(chains)$weight/2, E(chains)$weight)
  # comm <- cluster_louvain(chains, resolution = 1.5) # again choose resolution to make max ~ 50
  # sapply(1:max(comm$membership),function(x)sum(comm$membership==x))
  # # 41 20 25 18 13 25 25 22 10 25 13 22 11 34 13 47  7

  # replace the largest cluster in rtc with spitted communities
  # and set bridging individuals to a seperate cluster
  df_tmp <- data.table(PT_ID = bridging_all, MEMBERSHIP = max(df$MEMBERSHIP) + 1)
  df <- rbind(df, df_tmp)
  tmp <- df[,list(CLU_SIZE = length(PT_ID)),by = 'MEMBERSHIP']
  df <- merge(df, tmp, by = 'MEMBERSHIP')
  setnames(df, c('MEMBERSHIP','PT_ID'), c('IDCLU','ID'))

  
  # # Take the bridging individuals from communities of bridging individuals and consider them as a separate cluster.
  # #
  # df_graph_tmp <- as_long_data_frame(chains)
  # df_graph_tmp <- data.table(df_graph_tmp[,3:5])
  # colnames(df_graph_tmp) <- c('PERC','PT_ID1','PT_ID2')
  # 
  # bridging <- take_bridging_individual(df_graph_tmp, df_tmp)
  # # total 371 - 316 bridging too many 
  # bridging_all <- unique(do.call(c,bridging))
  # id_clu <- max(df$MEMBERSHIP) + 1
  # df <- rbind(df, data.table(PT_ID = bridging_all, MEMBERSHIP = id_clu))
  # tmp <- df[,list(CLU_SIZE = length(PT_ID)),by = 'MEMBERSHIP']
  # df <- merge(df, tmp, by = 'MEMBERSHIP')
 

  # merge with rtc
  rtc <- rtc[IDCLU != 1, ]
  rtc[, IDCLU := IDCLU + max(df$IDCLU) - 1]
  df <- rbind(rtc, df)
  table(table(df$IDCLU))

  cat('Write clusters to ', paste0(indir,'clusters.rds'),'...\n')
  saveRDS(df, file = paste0(indir,'clusters_potential_nets_1.rds'))
  
  # not in close pairs
  length(unique(df$ID))
  
  #
  df <- readRDS(paste0(indir,'clusters_potential_nets_1.rds'))
  dts <- readRDS("~/dts.rds")
  load("~/cluster_1.rda")
  
  #
  inds_notin <- setdiff(unique(dts$PT_ID), unique(df$ID)) 
  inds_in_iso <- df$ID[df$CLU_SIZE==1]
  
  #
  inds_in_rmneign <- dts[dts$PANGEA_ID %in% dcs[dcs$N_IND_CLOSE > 50, ]$TAXA,]
  
  # for those removed due to massive neighbors, add 10 closest individuals
  nclose <- 10
  minperc <- c()
  for(i in 1:nrow(inds_in_rmneign)){
    id <- inds_in_rmneign$PANGEA_ID[i]
    das <- readRDS(paste0("avescores_",id,".rds"))
    das <- merge(das, dts, by.x = "TAXA1", by.y = "PANGEA_ID", all.x = T)
    das <- merge(das, dts, by.x = "TAXA2", by.y = "PANGEA_ID", all.x = T)
    das$TAXA1 <- das$TAXA2 <- NULL
    setnames(das, c("PT_ID.x","PT_ID.y"), c("TAXA1", "TAXA2"))
    # sort(das$LENGTH)
    # print(sort(das$PERC, decreasing = T))
    # only trust scores if LEN > 500 or median length (some lengths are small)
    # CLOSE if >= 0.975
    # max add 10
    tmp <- das[das$LENGTH >= min(500, median(das$LENGTH)),]
    tmp <- tmp[order(-tmp$PERC),]
    # in case closest few sequences are from the same ind
    tmptx <- unique(tmp$TAXA2)
    tmptx <- tmptx[1:min(nclose, nrow(tmp))]
    # record min scores among these closest ind
    tmp <- tmp[tmp$TAXA2 %in% tmptx,]
    tmp <- tmp[,list(PERC = max(PERC)), by =c("TAXA1","TAXA2")]
    minperc <- c(minperc, min(tmp$PERC))
    # & das$PERC >= 0.9
    if(i == 1){
      df_in_rmneigh <- data.table(ID = c(unique(tmp$TAXA1),tmptx),
                                  IDCLU = i)
    }else{
      df_in_rmneigh <- rbind(df_in_rmneigh, data.table(ID = c(unique(tmp$TAXA1),tmptx),
                                                       IDCLU = i))
    }
  }
  
  range(minperc) # 0.9410569 0.9969040
 
  
  # for those isolated vertex and those not in potential nets, add 3 closest individuals
  nclose <- 3
  minperc <- c()
  # which(inds_in_iso == "RK-K127989")
  for(i in 1:length(inds_in_iso)){
    cat (i,"\n")
    ids <- dts[dts$PT_ID %in% inds_in_iso[i],]$PANGEA_ID
    dass <- list()
    for(id in ids){
      dass[[id]] <- readRDS(paste0("avescores_",id,".rds"))
    }
    das <- rbindlist(dass)
    das <- merge(das, dts, by.x = "TAXA1", by.y = "PANGEA_ID", all.x = T)
    das <- merge(das, dts, by.x = "TAXA2", by.y = "PANGEA_ID", all.x = T)
    das$TAXA1 <- das$TAXA2 <- NULL
    setnames(das, c("PT_ID.x","PT_ID.y"), c("TAXA1", "TAXA2"))
    
    # sort(das$LENGTH)
    # print(sort(das$PERC, decreasing = T))
    # only trust scores if LEN > 500 or median length (some lengths are small)
    # max add 3
    #
    tmp <- das[das$LENGTH >= min(500, median(das$LENGTH)),]
    tmptx <- unique(tmp[order(-tmp$PERC),]$TAXA2)
    tmptx <- tmptx[1:min(nclose, nrow(tmp))]
    #
    tmp <- tmp[tmp$TAXA2 %in% tmptx,]
    tmp <- tmp[,list(PERC = max(PERC)), by =c("TAXA1","TAXA2")]
    minperc <- c(minperc, min(tmp$PERC))
    #
    if(i == 1){
      df_in_iso <- data.table(ID = c(unique(tmp$TAXA1),tmptx),
                              IDCLU = i)
    }else{
      df_in_iso <- rbind(df_in_iso, data.table(ID = c(unique(tmp$TAXA1),tmptx),
                                               IDCLU = i))
    }
  }
  range(minperc) 
  # 0.8681818 0.9833659
  # which(minperc < 0.9) # 193 # check
  # das <- das[das$LENGTH !=-1,]
  # das <- das[das$LENGTH>200,] # there are closer inds but length ~ 250, not sure whether enough
  # das[order(-das$PERC),]
  # hist(das$LENGTH,plot=F)
  # names(which(table(dts$PT_ID) >1))[names(which(table(dts$PT_ID) >1)) %in% inds_in_iso]
  
  # for those not in potential nets, same as isolated vertices
  nclose <- 3
  minperc <- c()
  # which(inds_notin == "RK-K127989")
  for(i in 1:length(inds_notin)){
    cat (i,"\n")
    ids <- dts[dts$PT_ID %in% inds_notin[i],]$PANGEA_ID
    dass <- list()
    for(id in ids){
      dass[[id]] <- readRDS(paste0("avescores_",id,".rds"))
    }
    das <- rbindlist(dass)
    das <- merge(das, dts, by.x = "TAXA1", by.y = "PANGEA_ID", all.x = T)
    das <- merge(das, dts, by.x = "TAXA2", by.y = "PANGEA_ID", all.x = T)
    das$TAXA1 <- das$TAXA2 <- NULL
    setnames(das, c("PT_ID.x","PT_ID.y"), c("TAXA1", "TAXA2"))
    
    # sort(das$LENGTH)
    # print(sort(das$PERC, decreasing = T))
    # only trust scores if LEN > 500 or median length (some lengths are small)
    # max add 3
    #
    tmp <- das[das$LENGTH >= min(500, median(das$LENGTH)),]
    tmptx <- unique(tmp[order(-tmp$PERC),]$TAXA2)
    tmptx <- tmptx[1:min(nclose, nrow(tmp))]
    #
    tmp <- tmp[tmp$TAXA2 %in% tmptx,]
    tmp <- tmp[,list(PERC = max(PERC)), by =c("TAXA1","TAXA2")]
    minperc <- c(minperc, min(tmp$PERC))
    #
    if(i == 1){
      df_notin <- data.table(ID = c(unique(tmp$TAXA1),tmptx),
                              IDCLU = i)
    }else{
      df_notin <- rbind(df_notin, data.table(ID = c(unique(tmp$TAXA1),tmptx),
                                               IDCLU = i))
    }
  }
  range(minperc) # 0.8811189 1.0000000
  
  # combine
  # unique(df$ID)
  df <- df[df$CLU_SIZE !=1, ]
  
  df_in_rmneigh[, IDCLU := IDCLU + max(df$IDCLU)]
  df <- rbind(df[,c("ID", "IDCLU")], df_in_rmneigh)
  
  df_in_iso[, IDCLU := IDCLU + max(df$IDCLU)]
  df <- rbind(df, df_in_iso)
  
  df_notin[, IDCLU := IDCLU + max(df$IDCLU)]
  df <- rbind(df, df_notin)
  
  tmp <- df[,list(CLU_SIZE = length(ID)),by = 'IDCLU']
  df <- merge(df, tmp, by = "IDCLU", all.x=T)
  
  #
  table(table(df$IDCLU))
  
  length(unique(dts$PT_ID))
  length(unique(df$ID)) # 5825
  saveRDS(df, paste0(indir, 'clusters_all_1.rds'))
  

}

rkuvri.make.phyloscanner.input.samples <- function()
{
  #' input person ids, selected samples, ananymise ids
  #' returns a data table summarising person id and bam files
  #'
  require(data.table)
  
  # file names
  indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
  netdir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/'
  infile <- file.path(indir, "240809_PANGEA2_RCCS_final_samples_minrun125_minnuc250_minusPoor1.rds")
  infile_anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521_UVRI/important_anonymisation_keys_240817.csv'
  outbase <- paste0(indir, "240809_RCCSUVRI_")
  #
  ds <- readRDS(infile)
  hivc.db.Date2numeric<- function( x )
  {
    if(!class(x)%in%c('Date','character'))	return( x )
    x	<- as.POSIXlt(x)
    tmp	<- x$year + 1900
    x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
    x	
  }
  ds$VISIT_DT_N <- hivc.db.Date2numeric(ds$VISIT_DT)
  ds <- ds[ , .SD[which.min(VISIT_DT_N)], by = "PT_ID"]
  # ds[order(VISIT_DT_N),]
  
  ds <- subset(ds, select = c("PANGEA_ID", "PT_ID", "REMAP_BAM", "REMAP_REF_FASTA"))
  
  # simplify REMAP_BAM and REMAP_REF_FASTA
  ds[, REMAP_BAM:= gsub(indir,'', REMAP_BAM)]
  ds[, REMAP_REF_FASTA:= gsub(indir,'', REMAP_REF_FASTA)]
  
  # double check file names as expected
  set(ds, NULL, 'SAMPLE_ID', ds[, gsub('\\.bam','',REMAP_BAM)])
  set(ds, NULL, 'SAMPLE_ID_CHECK', ds[, gsub('_ref.fasta','',REMAP_REF_FASTA)])
  # ds[SAMPLE_ID!=SAMPLE_ID_CHECK, ]
  
  # finalise SAMPLE_ID
  set(ds, NULL, c('SAMPLE_ID_CHECK','REMAP_REF_FASTA','REMAP_BAM'), NULL)
  
  # add anonymised ids
  aid <- data.table(read.csv(infile_anonymised, stringsAsFactors = FALSE))
  aid <- subset(aid, select=c('PT_ID','AID'))
  ds <- merge(ds, aid, by.x = "PT_ID", by.y = "PT_ID", all.x = T)
  # ds[is.na(AID),]
  setnames(ds, c('PT_ID','AID'), c('UNIT_ID','RENAME_ID'))


  # finalise RENAME_ID
  ds[,FQ:=seq_len(length(SAMPLE_ID)), by='RENAME_ID']
  ds[,FQ:=paste0('fq',FQ)]
  ds[,RENAME_ID:= paste0(RENAME_ID,'-',FQ)]
  ds[,FQ:=NULL]
  
  # write processed samples
  cat('\nWriting to file ', paste0(outbase,'phscinput_samples.rds') )
  saveRDS(ds, file=paste0(outbase,'phscinput_samples.rds'))
}

rkuvri.make.phyloscanner.input.runs <- function()
{
  #' input clusters, couples, and three closest persons for each cluster
  #' return a data table allocating each person to a run
  
  require(data.table)
  require(tidyverse)
  
  indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
  outbase <- paste0(indir, '240809_RCCSUVRI_')
  netdir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/'
  infile_net <- file.path(netdir, 'clusters_all_1.rds')
  infile_couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
  
  # load data on individuals, sample locations, and rename ids
  pty.runs <- readRDS(paste0(outbase,'phscinput_samples.rds'))
  
  # load potential transmission networks of close pairs and redefine IDCLU
  df <- readRDS(infile_net)
  tmp <- unique(subset(df, select=c('CLU_SIZE','IDCLU')))
  setkey(tmp, CLU_SIZE)
  tmp[, IDCLU2 := seq_len(nrow(tmp))]
  tmp[, CLU_SIZE := NULL]
  
  rtc	<- subset( merge(df, tmp, by='IDCLU'),
                 select = c('ID', 'IDCLU2', 'CLU_SIZE'))
  setnames(rtc, 'IDCLU2', 'IDCLU')
  setkey(rtc, IDCLU)
  cat(max(rtc$IDCLU),' potential transmission networks of ',length(unique(rtc$ID)), ' individuals, \nThe ten largest cluster sizes are ',
      unique(subset(rtc, select=c('CLU_SIZE','IDCLU')))$CLU_SIZE[(max(rtc$IDCLU)-10+1):max(rtc$IDCLU)], "\n")
  # 3754  potential transmission networks of  5825  individuals, 
  # The ten largest cluster sizes are  33 35 39 40 42 46 50 51 51 55 
  
  #	add couples that are not a close pair and redefine IDCLU
  #
  #	find which couples are in potential transmission networks
  load(infile_couple)
  rp <- data.table(unique(subset(coupdat, select=c('male.RCCS_studyid', 'female.RCCS_studyid'))))
  setnames(rp, c('male.RCCS_studyid', 'female.RCCS_studyid'), c('MALE_RID','FEMALE_RID'))
  rp <- subset(rp, MALE_RID!=FEMALE_RID)
  rp[, FEMALE_RID:=paste0('RK-',FEMALE_RID)]
  rp[, MALE_RID:=paste0('RK-',MALE_RID)]
  tmp <- pty.runs$UNIT_ID
  rp <- unique(subset(rp, FEMALE_RID %in% tmp & MALE_RID %in% tmp, select = c(FEMALE_RID,MALE_RID)))

  # 
  rp <- merge(rp, subset(rtc, select = c(ID, IDCLU)), by.x = "FEMALE_RID", by.y = "ID", all.x = T)
  setnames(rp, "IDCLU", "FEMALE_IDCLU")
  rp <- merge(rp, subset(rtc, select = c(ID, IDCLU)), by.x = "MALE_RID", by.y = "ID", all.x = T)
  setnames(rp, "IDCLU", "MALE_IDCLU")
  rp[, COUP_ID := seq_len(nrow(rp))]
  
  #	reassign partners that are not in the same network as their partner
  tmp <- subset(rp, !is.na(FEMALE_IDCLU) & !is.na(MALE_IDCLU)) 
  tmp2 <- tmp[, list(NOT_IN_SAME = !any(FEMALE_IDCLU == MALE_IDCLU)),
              by=c('MALE_RID','FEMALE_RID')]
  tmp2 <- subset(tmp2, NOT_IN_SAME)
  for(i in 1:nrow(tmp2)){
    tmp2_male <- rtc[ID == tmp2$MALE_RID[i],]
    tmp2_male$CAT <- "M"
    tmp2_female <- rtc[ID == tmp2$FEMALE_RID[i],]  
    tmp2_female$CAT <- "F"
    tmp2_all <- rbind(tmp2_male, tmp2_female)
    tmp2_all <-  tmp2_all[order(CLU_SIZE),]
    # add to smallest cluster one of the couple is in
    if(tmp2_all$CAT[1] == "M"){
      rtc <- rbind(rtc, data.table(ID = tmp2$FEMALE_RID[i], IDCLU = tmp2_all$IDCLU[1]), fill = T)
    }else{
      rtc <- rbind(rtc, data.table(ID = tmp2$MALE_RID[i], IDCLU = tmp2_all$IDCLU[1]), fill = T)
    }
  }

  # tmp2 <- merge(tmp2, rp, by = c('MALE_RID','FEMALE_RID'))
  # 
  # tmp <- rp[, which(COUP_ID %in% tmp2$COUP_ID)]
  # set(rp, tmp, 'MALE_IDCLU', rp[tmp, FEMALE_IDCLU])
  # set(tmp2, NULL, 'FEMALE_IDCLU', tmp2[,MALE_IDCLU])
  # set(tmp2, NULL, 'NOT_IN_SAME', NULL)
  # rp <- rbind(rp, tmp2)
  # setnames(rp, c('MALE_IDCLU'), c('IDCLU'))
  # rp <- subset(melt(rp, id.vars=c('IDCLU'), 
  #                   measure.vars=c('MALE_RID','FEMALE_RID'), 
  #                   value.name='ID'), select=c(ID, IDCLU))
  # rtc <- unique(rbind(subset(rtc,select=c('ID','IDCLU')), rp))
  
  #	recalculate length of clusters
  rtc$CLU_SIZE <- NULL
  tmp <- rtc[, list(CLU_SIZE=length(ID)), by='IDCLU']
  setkey(tmp, CLU_SIZE)
  tmp[, IDCLU2 := seq_len(nrow(tmp))]
  rtc <- subset(merge(rtc, tmp, by='IDCLU', all.x = T),
                select=c('ID','IDCLU2','CLU_SIZE'))
  setnames(rtc, 'IDCLU2', 'IDCLU')
  setkey(rtc,IDCLU)
  cat(max(rtc$IDCLU),' potential transmission networks of ',length(unique(rtc$ID)), ' individuals, \nThe ten largest cluster sizes are ',
      unique(subset(rtc, select=c('CLU_SIZE','IDCLU')))$CLU_SIZE[(max(rtc$IDCLU)-10+1):max(rtc$IDCLU)], "\n")
  
  # 3754  potential transmission networks of  5825  individuals, 
  # The ten largest cluster sizes are  35 37 39 40 42 46 50 51 52 55
  # 
  #	merge up to 50 individuals into the same phyloscanner run
  tn	<- 50
  tmp	<- unique(subset(rtc, select=c(IDCLU, CLU_SIZE)))
  tmp[, PTY_RUN:= IDCLU]
  tmp[, PTY_SIZE:= CLU_SIZE]
  # https://stackoverflow.com/questions/49076769/dplyr-r-cumulative-sum-with-reset
  sum_reset_at <- function(thresh) {
    function(x) {
      accumulate(x, ~if_else(.x+.y>thresh, .y, .x+.y))
    }
  }
  tmp <- tmp %>% mutate(PTY_SIZE2 = sum_reset_at(tn)(PTY_SIZE))
  tmp <- as.data.table(tmp)
  tmp2 <- which(diff(tmp$PTY_SIZE2)<0) + 1
  tmp[, PTY_RUN2 := 0]
  tmp[c(tmp2,max(tmp2):nrow(tmp)), PTY_RUN2:=1]
  tmp[,PTY_RUN2:=cumsum(PTY_RUN2)]
  stopifnot(sum(tmp[,cumsum(PTY_SIZE),by='PTY_RUN2']$V1 !=tmp$PTY_SIZE2)==0)
  tmp[,PTY_SIZE2:=max(PTY_SIZE2),by='PTY_RUN2']
  tmp[,PTY_RUN2:=PTY_RUN2+1]
  set(tmp,NULL,c('PTY_RUN', 'PTY_SIZE'),NULL)
  setnames(tmp,c('PTY_RUN2', 'PTY_SIZE2'), c('PTY_RUN', 'PTY_SIZE'))
  rtc <- merge(rtc,tmp,by=c('IDCLU', 'CLU_SIZE'))
  # table(rtc$PTY_SIZE)
  saveRDS(rtc, file=paste0(outbase,'phscinput_runs_without_ctrl.rds'))
  
  #	for all individuals in a cluster
  #	find 3 closest others
  dts <- readRDS("~/dts.rds")
  rtc <- readRDS(paste0(outbase,'phscinput_runs_without_ctrl.rds'))
  nclose <- 3
  for(i in 1:max(rtc$PTY_RUN)){
    cat(i, "\n")
    ptids <- rtc[PTY_RUN == i,]$ID
    pgids <- dts$PANGEA_ID[dts$PT_ID %in% ptids]
    das <- lapply(pgids, function(pgid) {readRDS(paste0(netdir, "avescores_",pgid,".rds"))})
    das <- rbindlist(das)
    das <- merge(das, dts, by.x = "TAXA2", by.y = "PANGEA_ID", all.x = T)
    das <- das[LENGTH >= 500, list(PERC = mean(PERC)), by = c("PT_ID")]
    das <- das[order(-das$PERC),]
    rtc <- rbind(rtc, data.table(ID = das$PT_ID[1:nclose], PTY_RUN = i), fill = T)
    saveRDS(das, paste0(netdir, "avescores_ptyrun",i,".rds"))
  }
  # those with PTY_RUN but without IDCLU are controls
  rtc[, PTY_SIZE:=length(ID), by='PTY_RUN']
  
  tmp <- unique(subset(rtc, select=c('PTY_RUN','PTY_SIZE')))
  table(tmp$PTY_SIZE)
  
  #
  # load("~/phsinput_runs_except_last.rda")
  # combine with sample info
  pty.runs <- merge(rtc, pty.runs, by.x = "ID", by.y ="UNIT_ID", all.x = T)
  pty.runs <- pty.runs[order(pty.runs$PTY_RUN),]
  
  # write processed samples
  cat('\nWriting to file ', paste0(outbase,'phscinput_runs.rds') )
  saveRDS(pty.runs, file = paste0(outbase,'phscinput_runs.rds'))
}

