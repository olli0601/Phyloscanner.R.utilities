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

rk.seq.make.subsequence <- function()
{
  #' break sequences into subsequences by windows

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



rk.seq.submit.make.distance <- function(){
  
  #' write script to calculate distance
  
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



rk.seq.gather.distance <- function(){

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

rk.seq.depths <- function(){
  library(ape)
  library(data.table)
  
  scriptdir <- '~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1521/'
  indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
  infile_sequence <- paste0(indir, "240809_PANGEA2_RCCS_final_samples_alignment_minusPoor1.fasta")
  infile <- paste0(indir, "240809_PANGEA2_RCCS_mapped_samples.rds")
  outfile <- paste0(indir, "240809_PANGEA2_RCCS_depths.rds")
  
  # load data
  alignment <- read.dna(file = infile_sequence, format = "fasta")
  nsequence <- nrow(alignment)
  npos <- ncol(alignment)
  
  #
  mapping_rccs <- readRDS(infile)
  dconsensus <- subset(mapping_rccs,select=c('PANGEA_ID','CONSENSUS','F'))
  dconsensus <- unique(dconsensus[CONSENSUS!="",])
  dconsensus <- subset(dconsensus, select=c('PANGEA_ID','F'))
  setnames(dconsensus, c('PANGEA_ID','F'), c('pangea_id','file_name'))
  tmp <-  data.table(pangea_id=rownames(alignment))
  dconsensus <- merge(dconsensus, tmp, by='pangea_id', all.y=T)
  dconsensus <- unique(dconsensus)
  # sort(table(dconsensus$pangea_id))
  # dconsensus[dconsensus$pangea_id == "PG15-UG504899",]
  
  # windows
  window_size <- 500
  windows_last_start <- ceiling(npos/100) * 100 - window_size + 1
  windows_2last_end <- floor(npos/100) * 100
  windows_start <- seq(1,windows_last_start,100)
  windows_end <- c(seq(window_size ,windows_2last_end,100),npos)
  windows <- seq_len(length(windows_start))
  dw <- data.table(WINDOW=windows,
                   START=windows_start,
                   END=windows_end)
  
  # 
  ddepth <- data.table()
  
  for (i in seq_len(nrow(dconsensus))) {
    if(i %% 100 == 1 | i==nrow(dconsensus)){
      cat('processing ', i, 'th out of ', nrow(dconsensus), ' files \n' )
    }
    sequence <- readLines(dconsensus$file_name[i])
    sequence <- paste(sequence[2:length(sequence)], collapse="")
    sequence <- strsplit(sequence,"")
    
    tmp_df <- dw[,{
      tmp = sequence[[1]][START:END]
      tmp = tmp[!grepl("[^A-Za-z]",tmp)]
      list(
        LEN_LETTER=length(tmp),
        HIGH_DEPTH=sum(tmp == toupper(tmp)),
        PANGEA_ID=dconsensus$pangea_id[i],
        FILE=dconsensus$file[i])
    }, by = c('WINDOW', 'START', 'END')]
    ddepth <- rbind(ddepth, tmp_df)
  }
  #
  ddepth[ddepth$PANGEA_ID == "PG15-UG504899" & ddepth$WIN == 1,]
  saveRDS(ddepth, file = outfile)
  
  
  ddepth <- readRDS(outfile)
  length_cut <- 250
  ddepth_select <- subset(ddepth,select=c('WINDOW','LEN_LETTER','HIGH_DEPTH','PANGEA_ID'))
  ddepth_select <- ddepth_select[HIGH_DEPTH >= length_cut, ]
  
  
  ddepth_select_w <- ddepth_select[,list(PANGEA_ID=unique(PANGEA_ID)),
                                   by=c('WINDOW')]
  ddepth_select_w[,HD:=1]
  ddepth_select_w <- data.table::dcast(ddepth_select_w, PANGEA_ID~WINDOW, value.var='HD')
 
  length(unique(ddepth$PANGEA_ID)) # 7691
  nrow(ddepth_select_w) # 7063
  sum(apply(as.matrix(ddepth_select_w[,-1]), 1, function(x)sum(x,na.rm=T)) >=4)
  # 6844
}



rk.seq.find.threshold <- function(){
  #' find thresholds from couples 

  library(ape)
  library(data.table)
  
  window_size <- 500
  windown <- 99
  
  scriptdir <- '~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1521/'
  indir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/'
  outdir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/potential_nets/'
  infiles <- paste0(outdir, 'subsequence',1:windown,'_results.rda')
  # infiles <- "/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/subsequence1_results_test.rda"
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
  sum(is.na(dinfo$PT_ID))
  
  # couples
  load(infile_couple)
  couple <- unique(subset(data.table(coupdat), select=c('male.RCCS_studyid', 'female.RCCS_studyid')))
  couple[, COUPLE:=1]
  couple[,male.RCCS_studyid:=paste0('RK-',male.RCCS_studyid)]
  couple[,female.RCCS_studyid:=paste0('RK-',female.RCCS_studyid)]
  
  # 
  tmp_couple <- subset(couple,select=c('male.RCCS_studyid','female.RCCS_studyid', 'COUPLE'))
  tmp <- copy(tmp_couple)
  setnames(tmp_couple, c('male.RCCS_studyid','female.RCCS_studyid'),c('pt_id1','pt_id2'))
  setnames(tmp, c('male.RCCS_studyid','female.RCCS_studyid'),c('pt_id2','pt_id1'))
  tmp_couple[,sex1_couple := 'M']
  tmp_couple[,sex2_couple := 'F']
  tmp[,sex1_couple := 'F']
  tmp[,sex2_couple := 'M']
  tmp_couple <- rbind(tmp_couple, tmp)
  unique(tmp_couple)
  tmp_couple  <- tmp_couple[!is.na(pt_id1) & !is.na(pt_id2) & pt_id1!=pt_id2,]
  
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
    df <- merge(df,dinfo, by.x=c('TAXA1'), by.y=c('PANGEA_ID1'), all.x=T) 
    setnames(dinfo, colnames(dinfo), gsub('1','2',colnames(dinfo)))
    df <- merge(df,dinfo, by.x=c('TAXA2'), by.y=c('PANGEA_ID2'), all.x=T) 
    setnames(dinfo, colnames(dinfo), gsub('2','',colnames(dinfo)))
    
    # remove no pt_id
    df <- df[!is.na(PT_ID1) & !is.na(PT_ID2) & PT_ID1!=PT_ID2,]
    df <- merge(df,tmp_couple, by=c('PT_ID1','PT_ID2'), all.x=T)
    df <- df[COUPLE==1,c('PT_ID1','PT_ID2','PERC','WINDOW')]
    if(i==1){
      ans <- df
    }else{
      ans <- rbind(ans, df)
    }
  }
  
  # save 
  saveRDS(ans, file = file.path(outdir,'results_couples.rds'))  
  
  # visualisation
  
  # fit stan
  
  
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
  
  # save(stan_data, file = file.path(potential.networks.analysis.dir,'stan_data_threshold.rda'))
  
  # library(data.table)
  tmp <- stan_data$similarity
  tmp <- data.table(t(tmp[seq(1,99,20),]))
  colnames(tmp) <- paste0(windows_start[seq(1,99,20)],' - ', windows_end[seq(1,99,20)])
  tmp <- melt(tmp)
  tmp <- tmp[value!=-1]
  
  g <- ggplot(tmp, aes(value))+
    geom_histogram()+
    facet_wrap(~variable,ncol = 2,scales = 'free') +
    theme_bw()+
    labs(x='\n similarity scores', y='count \n') +
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0.05))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          text = element_text(size=20))
  ggsave(file.path(potential.networks.analysis.dir,'plots','raw_histogram_distance_couple_per_window_sample5.pdf'),g,width = 8, height= 6)
  
  
  # fit bimodal model to similarities
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
  for (i in 1:windown) {
    cat('fit window ', i, '\n')
    tmp <- list()
    tmp$y <- (stan_data$similarity[i,1:stan_data$Nsimilarity[i]])
    tmp$N <- stan_data$Nsimilarity[i]
    fit[[i]] <- stan(model_code=stan_code,
                     data=tmp,
                     chains=1, seed=42)
    
  }
  
  # save fitted model
  save(fit,file=file.path(potential.networks.analysis.dir,'indep_exp5_fit.rda'))
  
  
  load(file.path(potential.networks.analysis.dir,'indep_exp5_fit.rda'))
  # check fit 
  rh= c()
  ness= c()
  tess= c()
  bess= c()
  for (i in 1:windown) {
    cat('process fit window ', i, '\n')
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
  df = data.table()
  for (i in 1:windown) {
    cat('process fit window ', i, '\n')
    tmp = data.table(summary(fit[[i]])$summary[1:5,c("2.5%","50%","97.5%")])
    tmp[,VAR:=rownames(summary(fit[[i]])$summary)[1:5]]
    tmp[,WIN:=i]
    df = rbind(df, tmp)
  }
  
  
  # visualise bimodal fit 
  ans[,PERC:=log(PERC)]
  plot_list =list()
  #  range(ans$PERC)
  #  -0.4453817  0.0000000
  
  tmp <- dcast(df, WIN~VAR, value.var='50%')
  setnames(tmp,'WIN','WINDOW')
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
    facet_wrap(~factor(WINDOW,1:length(windows_start), paste0(windows_start,' - ',windows_end)),ncol = 8,scales = 'free') +
    theme_bw()+
    labs(x='log similarity scores') +
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0.05))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  ggsave(file.path(potential.networks.analysis.dir,'plots','fitted_histogram_distance_couple_per_window_exp5_indep.pdf'),g,width = 12, height= 16)
  
  g <- ggplot(ans[WINDOW%%20==1], aes(PERC))+
    geom_histogram(binwidth = bw)+
    geom_line(tmp[WINDOW%%20==1],mapping=aes(x=X,y=Y1))+
    geom_line(tmp[WINDOW%%20==1],mapping=aes(x=X,y=Y2))+
    facet_wrap(~factor(WINDOW,1:length(windows_start), paste0(windows_start,' - ',windows_end)),ncol = 2,scales = 'free') +
    theme_bw()+
    labs(x='\n log similarity scores', y='count \n') +
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0.05))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          text = element_text(size=20))
  ggsave(file.path(potential.networks.analysis.dir,'plots','fitted_histogram_distance_couple_per_window_exp5_indep_sample5.pdf'),g,width = 8, height= 6)
  
  
  # plot quantiles of fits
  for (i in 1:windown) {
    cat('process window ',i, '\n')
    x <- seq(-0.45,0,length.out = 90 +1)
    df_plot_fit <- data.table(x=x,
                              c1 = df[VAR=='theta' & WIN ==i,]$`50%` * dnorm(x, df[VAR=='mu[1]' & WIN ==i,]$`50%`, df[VAR=='sigma[1]' & WIN ==i,]$`50%`),
                              c2 = (1-df[VAR=='theta' & WIN ==i,]$`50%`) *dnorm(x, df[VAR=='mu[2]' & WIN ==i,]$`50%`, df[VAR=='sigma[2]' & WIN ==i,]$`50%`))
    g = ggplot(df_plot_fit)+
      geom_line( aes(x,c1), color='red') +
      geom_line(aes(x,c2), color='blue') +
      geom_vline(aes(xintercept = qnorm(0.05,  df[VAR=='mu[1]' & WIN ==i,]$`50%`, df[VAR=='sigma[1]' & WIN ==i,]$`50%`)),
                 linetype=2, color='red') +
      geom_vline(aes(xintercept = qnorm(0.5,  df[VAR=='mu[1]' & WIN ==i,]$`50%`, df[VAR=='sigma[1]' & WIN ==i,]$`50%`)),
                 linetype=2, color='red') +
      geom_vline(aes(xintercept = qnorm(0.95,  df[VAR=='mu[1]' & WIN ==i,]$`50%`, df[VAR=='sigma[1]' & WIN ==i,]$`50%`)),
                 linetype=2, color='red') +
      geom_vline(aes(xintercept = qnorm(0.05,  df[VAR=='mu[2]' & WIN ==i,]$`50%`, df[VAR=='sigma[2]' & WIN ==i,]$`50%`)),
                 linetype=2, color='blue') +
      geom_vline(aes(xintercept = qnorm(0.5,  df[VAR=='mu[2]' & WIN ==i,]$`50%`, df[VAR=='sigma[2]' & WIN ==i,]$`50%`)),
                 linetype=2, color='blue') +
      geom_vline(aes(xintercept = qnorm(0.95,  df[VAR=='mu[2]' & WIN ==i,]$`50%`, df[VAR=='sigma[2]' & WIN ==i,]$`50%`)),
                 linetype=2, color='blue') +
      theme_bw() +
      labs(x="log similarities", y="counts")+ggtitle(paste0(windows_start[i], ' - ', windows_end[i]))+
      theme(axis.text.x=element_text(angle=45, hjust=1))
    
    plot_list[[i]] = ggplotGrob(g)
  }
  
  # save densities of fitted plots and quantiles
  pdf(file.path(potential.networks.analysis.dir,'plots','density_distance_couple_per_window_quantiles_exp5_indep.pdf'),
      width = 30, height= 45)
  do.call("grid.arrange", c(plot_list,ncol= 5))
  dev.off()
  
  tmp = df[grep('mu',VAR),]
  tmp = dcast(tmp, WIN~VAR,value.var = '50%')
  tmp[,diff:=`mu[2]`-`mu[1]`]
  # mean and sd over windows
  ggplot(df[grep('mu',VAR)],aes(x=WIN))+
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
    geom_line(aes(y = `50%`))+
    facet_wrap(VAR~.,ncol = 1) +
    labs(x='window', y='means')
  # save component means across windows 
  ggsave(file.path(potential.networks.analysis.dir,'plots','component_mean_vs_window_exp5_indep.pdf'),
         width = 6, height= 8)
  
  ggplot(df[grep('sigma',VAR)],aes(x=WIN))+
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
    geom_line(aes(y = `50%`))+
    facet_wrap(VAR~.,ncol = 1) +
    labs(x='window', y='standard deviations')
  # save component standard deviations across windows 
  ggsave(file.path(potential.networks.analysis.dir,'plots','component_sd_vs_window_exp5_indep.pdf'),
         width = 6, height= 8)
  
  # method 2: 5% percentile 
  # currently used
  for (i in 1:windown) {
    cat('process fit window ', i, '\n')
    if(i==1)threshold = list()
    tmp =c()
    pars = rstan::extract(fit[[i]], pars=c('mu','sigma','theta'))
    for (j in 1:nrow(pars$mu)) {
      tmp = c(tmp, qlnorm(0.05, pars$mu[j,2], pars$sigma[j,2]))
      
    }
    threshold[[i]] = tmp
  }
  
  
  threshold = lapply(threshold,function(x){quantile(x,probs = c(0.025,0.5,0.975))})
  threshold = data.table(do.call(rbind, threshold))
  save(threshold, file = file.path(potential.networks.analysis.dir,paste0('threshold_5quantile_indep_exp5.rda')))
  
  # plot threshold
  load(file.path(potential.networks.analysis.dir,paste0('threshold_5quantile_indep_exp5.rda')))
  tmp <- data.table(threshold)
  # tmp[,WSTART:=windows_start]
  tmp[,W:=seq_len(nrow(tmp))]
  
  # ggplot(tmp,aes(x=WSTART,y=`50%`,ymin=`2.5%`,ymax=`97.5%`))+
  ggplot(tmp,aes(x=W,y=`50%`,ymin=`2.5%`,ymax=`97.5%`))+
    geom_line()+
    geom_ribbon(alpha=0.6)+
    theme_bw()+
    # labs(x='\n genomic positions',y='threshold \n')+
    labs(x='\n genomic windows',y='threshold \n')+
    scale_x_continuous(expand = c(0,0),breaks = seq(0,100,by=10)) +
    scale_y_continuous(expand = c(0,0), labels = scales::percent, limits = c(0.85,1))+
    theme(text = element_text(size=20))
  
  ggsave(file.path(potential.networks.analysis.dir,'plots',paste0('threshold_vs_window.pdf')),width = 6, height = 4)
  
  # ggsave(file.path(potential.networks.analysis.dir,'plots',paste0('threshold_over_positions.pdf')),width = 6, height = 4)
  
}