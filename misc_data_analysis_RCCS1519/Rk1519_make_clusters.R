library(data.table)
library(ggplot2)
library(seqinr)
library(ggpubr)
library(tidyverse)
library(igraph)
library(dplyr)
library(pammtools)# plot stepribbon


rkuvri.make.subsequences <- function()
{
  #' break sequences into subsequences by windows
  #' it takes sequences and break them into overlapping windows of length 500 and defines batches for distance calculation
  
  # break sequences into windows and batches
  window_size <- 500
  batch_size <- 100
  
  # file names
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  # out.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_'
  potential.networks.analysis.dir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100"
  infile.sequence <- file.path(data.dir,"200422_PANGEA2_RCCSMRC_alignment.fasta")

  
  # load data
  alignment <- read.fasta(file = infile.sequence)
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
  dfo[, BATCH:= seq_len(nrow(dfo))%%batch_size+1L]
  
  # write batches to rds
  cat('\nWriting to file ', file.path(potential.networks.analysis.dir,'batch_windowsize500_batchsize100.rds') )
  saveRDS(dfo,file = file.path(potential.networks.analysis.dir,'batch_windowsize500_batchsize100.rds'))
  
}


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


rkuvri.make.distance<- function()
{
  #'  calculate length of valid sequences, total match scores, percentage of matches of pair of sequences in one batch and in one window
  #' it inputs the subsequence and batch info, and returns a data.table with distances
  
  
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
    stopifnot(args_line[[7]]=='-batchi')
    stopifnot(args_line[[9]]=='-windowi')
    args_dir <- list()
    args_dir[['data_dir']] <- args_line[[2]]
    args_dir[['out_dir']] <- args_line[[4]]
    args_dir[['script_dir']] <- args_line[[6]]
    args_dir[['batchi']] <- as.integer(args_line[[8]])
    args_dir[['windowi']] <- as.integer(args_line[[10]])
  } 
  
  # file names
  infile.batch <- file.path(args_dir[['out_dir']],paste0('batch_windowsize500_batchsize100.rds'))
  infile.subsequence <- file.path(args_dir[['out_dir']], paste0('subsequence_window_', args_dir[['windowi']], '.fasta'))
  
  #	read sequences	
  sq<- read.fasta(infile.subsequence)
  
  #	prepare and select batch
  dfo <- readRDS(infile.batch)
  
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


rkuvri.submit.make.distance <- function(){
  #' write script to calculate distance
  
  window_size <- 500
  batch_size <- 100
  
  # files
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  potential.networks.analysis.dir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100"
  infile.sequence <- file.path(data_dir,"200422_PANGEA2_RCCSMRC_alignment.fasta")

  
  # load data
  alignment <- read.fasta(file = infile.sequence)
  nsequence <- length(alignment)
  npos <- unique(lengths(alignment))
  
  # windows 
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
                          hpc.load= "module load anaconda3/personal\nsource activate Renv", 
                          hpc.array= jobn,
                          hpc.q = NA)
  
  cmd <- paste0(cmd ,'\n case $PBS_ARRAY_INDEX in \n')
  for (i in seq_len(jobn)) {
    cmd <- paste0(cmd, i,') \n')
    cmd <- paste0(cmd, 'Rscript ', '~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/distance.R -data_dir ',  data.dir, 
                  ' -out_dir ', potential.networks.analysis.dir,
                  ' -script_dir ~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/',
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

rkuvri.collect.distance.per.person.pair <- function(){
  #' input similarities per batch and per window
  #' combine similarities to one data.table for each window 
  #' summarise basic infomations per window and clean by removing invalid pairs
  
  batch_size <- 100
  windown <- 99
  
  # file names
  potential.networks.analysis.dir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100"
  infile.dist.base <- file.path(potential.networks.analysis.dir, paste0('subsequence_window_',1:windown))

    # clean results 
  # gather outputs for each window
  for (w in 1:windown) {
    base <- infile.dist.base[w]
    distancei <- data.table()
    for (i in 1:batch_size) {
      tmp <- readRDS(paste0(base,'_batch',i,'.rds'))
      distancei <- rbind(distancei,tmp)
      file.remove(file.path(potential.networks.analysis.dir,paste0(base,'_batch',i,'.rds')))
    }
    
    # summarise
    metai <- data.table(window=k,
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
    distancei[,WINDOW:=k]
    
    # save combined results per window
    save(distancei,metai,file=paste0(infile.dist.base[w],'_results.rda'))
    gc()
  }
}


rkuvri.find.threshold.from.couples<- function(){
  #' input couples, sequences, person ids, similarity results, 
  #' take similarity data for couples, fit bimodal distributions to similarities, define the threshold
  #' the most useful output would be threshold files
  #' Note: eligible sequences should satisfy length >= 250, have person id
  #' eligible pairs of sequences should not be the same person
  
  # file names
  window_size <- 500
  windown <- 99
  
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  out.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_'
  potential.networks.analysis.dir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100"
  infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
  infile.sequence <- file.path(data.dir,"200422_PANGEA2_RCCSMRC_alignment.fasta")
  infile.ind.rccs <- file.path(data.dir,'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.mrc <- file.path(data.dir,'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
  infile.dist <- file.path(potential.networks.analysis.dir, paste0('subsequence_window_',1:windown,'results.rda'))
  
  # load data
  alignment <- read.fasta(file = infile.sequence)
  unique(do.call(c, lapply(alignment, unique)))
  nsequence <- length(alignment)
  npos <- unique(lengths(alignment))
  
  # windows
  windows_last_start <- ceiling(npos/100) * 100 - window_size + 1
  windows_2last_end <- floor(npos/100) * 100
  windows_start <- seq(1,windows_last_start,100)
  windows_end <- c(seq(window_size ,windows_2last_end,100),npos)
  windows <- seq_len(length(windows_start))
  dw <- data.table(WINDOW=windows,
                   START=windows_start,
                   END=windows_end)
  
  
  # map alignments to studyid
  dinfo <- data.table(pangea_id=names(alignment))
  id.dt <- data.table(read.csv(infile.ind.rccs))
  id.dt <- subset(id.dt,select = c("pt_id","sex","pangea_id"))
  id.dt[,pangea_id:=paste0('RCCS_',pangea_id)]
  tmp <- data.table(read.csv(infile.ind.mrc))
  tmp <- subset(tmp,select = c("pt_id","sex","pangea_id"))
  tmp[,pangea_id:=paste0('MRCUVRI_',pangea_id)]
  id.dt <- rbind(id.dt,tmp)
  id.dt <- unique(id.dt)
  dinfo <- merge(dinfo, id.dt, by="pangea_id", all.x=T)
  tmp <- dinfo[is.na(pt_id)]
  dinfo <- dinfo[!is.na(pt_id)]
  tmp2 <- subset(rccsData,select=c('RCCS_studyid','Pangea.id','SEX'))
  colnames(tmp2) <- c("pt_id","pangea_id","sex")
  set(tmp, NULL, c("pt_id","sex"), NULL)
  tmp <- merge(tmp, tmp2, by='pangea_id',all.x=T)
  dinfo <- rbind(dinfo, tmp)
  cat('No personal information found for ',nrow(dinfo[is.na(pt_id)]), ' sequences \n')
  
  # process couple, map to studyid
  load(infile.couple)
  couple <- unique(subset(data.table(coupdat), select=c('male.RCCS_studyid', 'female.RCCS_studyid')))
  couple[, COUPLE:=1]
  couple[,male.RCCS_studyid:=paste0('RK-',male.RCCS_studyid)]
  couple[,female.RCCS_studyid:=paste0('RK-',female.RCCS_studyid)]
  
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
  
  # calculate similarities of pairs among couples on each window
  for (i in 1:windown) {
    # similarity per window
    cat('\n process window ',i , '\n')
    load(infile.dist[i])
    #
    distancei <- distancei[LENGTH >= length_cutoff,]
    # add couple
    setnames(dinfo, colnames(dinfo), paste0(colnames(dinfo),'1'))
    distancei <- merge(distancei,dinfo, by.x=c('TAXA1'), by.y=c('pangea_id1'), all.x=T) 
    setnames(dinfo, colnames(dinfo), gsub('1','2',colnames(dinfo)))
    distancei <- merge(distancei,dinfo, by.x=c('TAXA2'), by.y=c('pangea_id2'), all.x=T) 
    setnames(dinfo, colnames(dinfo), gsub('2','',colnames(dinfo)))
    # remove no pt_id
    distancei <- distancei[!is.na(pt_id1) & !is.na(pt_id2) & pt_id1!=pt_id2,]
    distancei <- merge(distancei,tmp_couple, by=c('pt_id1','pt_id2'), all.x=T)
    distancei <- distancei[COUPLE==1,c('pt_id1','pt_id2','PERC','WINDOW')]
    if(i==1){
      ans <- distancei
    }else{
      ans <- rbind(ans, distancei)
    }
    
  }
  
  # save 
  save(ans, file = file.path(potential.networks.analysis.dir,'distance_couple_per_window.rda'))  
  
  # make stan data
  # load(file.path(potential.networks.analysis.dir,'distance_couple_per_window.rda'))
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
  
  # plot fit to similarities
  for (i in 1:windown) {
    cat('process window ',i, '\n')
    x <- seq(-0.45,0,length.out = 90 +1)
    df_plot_fit <- data.table(x=x,
                              c1 = df[VAR=='theta' & WIN ==i,]$`50%` * dnorm(x, df[VAR=='mu[1]' & WIN ==i,]$`50%`, df[VAR=='sigma[1]' & WIN ==i,]$`50%`),
                              c2 = (1-df[VAR=='theta' & WIN ==i,]$`50%`) *dnorm(x, df[VAR=='mu[2]' & WIN ==i,]$`50%`, df[VAR=='sigma[2]' & WIN ==i,]$`50%`))
    max_hist=  max(hist(ans[WINDOW==i,]$PERC,plot=F, breaks = seq(-0.45,0,length.out = 30 + 1), right = FALSE)$counts)
    scale <- max_hist/max(c(df_plot_fit$c1,df_plot_fit$c2))
    g = ggplot(ans[WINDOW==i,],aes(x=PERC))+
      geom_histogram(alpha=0.6)+
      geom_line(data=df_plot_fit, aes(x,c1 * scale), color='red') +
      geom_line(data=df_plot_fit, aes(x,c2 * scale), color='blue') +
      scale_x_continuous(breaks = seq(-0.45, 0,length.out = 30 + 1))+
      theme_bw() +
      labs(x="log similarities", y="counts")+ggtitle(paste0(windows_start[i], ' - ', windows_end[i]))+
      theme(axis.text.x=element_text(angle=45, hjust=1))
    
    plot_list[[i]] = ggplotGrob(g)
  }
  
  # create dir to save plots related to network analysis
  dir.create(file.path(potential.networks.analysis.dir,'plots'))
  # save histogram and fitted plots
  pdf(file.path(potential.networks.analysis.dir,'plots','fitted_histogram_distance_couple_per_window_exp5_indep.pdf'),
      width = 30, height= 45)
  do.call("grid.arrange", c(plot_list,ncol= 5))
  dev.off()
  
  
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
  
  
  ##### find threshold 
  # method 1: intersection
  # used initially, but not meaningful
  for (i in 1:windown) {
    cat('process fit window ', i, '\n')
    if(i==1)threshold = list()
    tmp =c()
    pars = rstan::extract(fit[[i]], pars=c('mu','sigma','theta'))
    for (j in 1:nrow(pars$mu)) {
      f <- function(x) dlnorm(x, pars$mu[j,1], pars$sigma[j,1]) * pars$theta[j] -
        dlnorm(x, pars$mu[j,2], pars$sigma[j,2]) * (1-pars$theta[j])
      tryCatch({
        tmp <- c(tmp, uniroot(f, interval=c(0.7,1))$root)
      }, error=function(e){cat("ERROR:",conditionMessage(e), "\n",
                               'on window ',i, ' iteration ', j , '\n')})
      
    }  
    threshold[[i]] = tmp
  }
  
  threshold = lapply(threshold,function(x){quantile(x,probs = c(0.025,0.5,0.975))})
  threshold = data.table(do.call(rbind, threshold))
  save(threshold, file = file.path(potential.networks.analysis.dir,paste0('threshold_dip_indep_exp5.rda')))
  
  
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
  
  
  # method 3: median 
  # 
  for (i in 1:windown) {
    cat('process fit window ', i, '\n')
    if(i==1)threshold = list()
    tmp =c()
    pars = rstan::extract(fit[[i]], pars=c('mu','sigma','theta'))
    for (j in 1:nrow(pars$mu)) {
      tmp = c(tmp, qlnorm(0.5, pars$mu[j,2], pars$sigma[j,2]))
      
    }
    threshold[[i]] = tmp
  }
  
  
  threshold = lapply(threshold,function(x){quantile(x,probs = c(0.025,0.5,0.975))})
  threshold = data.table(do.call(rbind, threshold))
  save(threshold, file = file.path(potential.networks.analysis.dir,paste0('threshold_median_indep_exp5.rda')))
  
  
}


rkuvri.summarise.distance.and.support.per.sequence.pair.per.window <- function(){
  #' input couples, sequences, person ids, consensus sequences info (or depth info), similarity results, thresholds
  #' returns average similarities and supports per pair of sequences
  #' Note: eligible sequences should satisfy length >= 250, have person id, have >250 high depth positions
  #' eligible pairs of sequences should not be the same person
  #' 
  window_size <- 500
  method <- 2
  windown <- 99 
  
  # file names
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  out.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_'
  potential.networks.analysis.dir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100"
  infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
  infile.sequence <- file.path(data.dir,"200422_PANGEA2_RCCSMRC_alignment.fasta")
  infile.ind.rccs <- file.path(data.dir,'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.mrc <- file.path(data.dir,'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
  infile.rccs <- file.path(data.dir,'PANGEA2_RCCS/200422_PANGEA2_RCCS_mapped_samples.rds')
  infile.mrc <- file.path(data.dir,'PANGEA2_MRC/200422_PANGEA2_MRCUVRI_mapped_samples.rds')
  infile.depth <- file.path(data.dir,'210121_RCCSMRC_depth.rda')
  infile.dist <- file.path(potential.networks.analysis.dir, paste0('subsequence_window_',1:windown,'results.rda'))
  if(method == 1){
    infile.threshold <- file.path(potential.networks.analysis.dir,paste0('threshold_dip_indep_exp5.rda'))
  }else if(method == 2){
    infile.threshold <- file.path(potential.networks.analysis.dir,paste0('threshold_5quantile_indep_exp5.rda'))
  }else if(method == 3){
    infile.threshold <- file.path(potential.networks.analysis.dir,paste0('threshold_median_indep_exp5.rda'))
  }
  
  # load data
  alignment <- read.fasta(file = infile.sequence)
  unique(do.call(c, lapply(alignment, unique)))
  nsequence <- length(alignment)
  npos <- unique(lengths(alignment))
  
  
  # windows
  windows_last_start <- ceiling(npos/100) * 100 - window_size + 1
  windows_2last_end <- floor(npos/100) * 100
  windows_start <- seq(1,windows_last_start,100)
  windows_end <- c(seq(window_size ,windows_2last_end,100),npos)
  windows <- seq_len(length(windows_start))
  dw <- data.table(WINDOW=windows,
                   START=windows_start,
                   END=windows_end)
  
  # map alignments to studyid
  load(infile.couple)
  dinfo <- data.table(pangea_id=names(alignment))
  id.dt <- data.table(read.csv(infile.ind.rccs))
  id.dt <- subset(id.dt,select = c("pt_id","sex","pangea_id"))
  id.dt[,pangea_id:=paste0('RCCS_',pangea_id)]
  tmp <- data.table(read.csv(infile.ind.mrc))
  tmp <- subset(tmp,select = c("pt_id","sex","pangea_id"))
  tmp[,pangea_id:=paste0('MRCUVRI_',pangea_id)]
  id.dt <- rbind(id.dt,tmp)
  id.dt <- unique(id.dt)
  dinfo <- merge(dinfo, id.dt, by="pangea_id", all.x=T)
  tmp <- dinfo[is.na(pt_id)]
  dinfo <- dinfo[!is.na(pt_id)]
  tmp2 <- subset(rccsData,select=c('RCCS_studyid','Pangea.id','SEX'))
  colnames(tmp2) <- c("pt_id","pangea_id","sex")
  set(tmp, NULL, c("pt_id","sex"), NULL)
  tmp <- merge(tmp, tmp2, by='pangea_id',all.x=T)
  dinfo <- rbind(dinfo, tmp)
  cat('No personal information found for ',nrow(dinfo[is.na(pt_id)]), ' sequences \n')
  
  # process couple, map to studyid
  load(infile.couple)
  couple <- unique(subset(data.table(coupdat), select=c('male.RCCS_studyid', 'female.RCCS_studyid')))
  couple[, COUPLE:=1]
  couple[,male.RCCS_studyid:=paste0('RK-',male.RCCS_studyid)]
  couple[,female.RCCS_studyid:=paste0('RK-',female.RCCS_studyid)]
  
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
  
  
  # consensus sequence depth
  if(file.exists(infile.depth)){
    # load depth
    load(infile.depth)
  }else{
    # read consensus sequences 
    mapping_rccs <- readRDS(infile.rccs)
    dconsensus <- subset(mapping_rccs,select=c('PANGEA_ID','CONSENSUS','F'))
    dconsensus <- unique(dconsensus[CONSENSUS!="",])
    dconsensus[,PANGEA_ID:=paste0('RCCS_',PANGEA_ID)]
    mapping_mrc <- readRDS(infile.mrc)
    tmp <- subset(mapping_mrc,select=c('PANGEA_ID','CONSENSUS','F'))
    tmp <- unique(tmp[CONSENSUS!="",])
    tmp[,PANGEA_ID:=paste0('MRCUVRI_',PANGEA_ID)]
    dconsensus = rbind(dconsensus, tmp)
    dconsensus = subset(dconsensus, select=c('PANGEA_ID','F'))
    setnames(dconsensus, c('PANGEA_ID','F'), c('pangea_id','file_name'))
    tmp <-  data.table(pangea_id=names(alignment))
    dconsensus <- merge(dconsensus, tmp, by='pangea_id', all.y=T)
    dconsensus <- unique(dconsensus)
    
    # valid lengths and high depth position lengths in consensus sequences
    ddepth = data.table()
    for (i in seq_len(nrow(dconsensus))) {
      if(i %% 100 == 1 | i==nrow(dconsensus)){
        cat('processing ', i, 'th out of ', nrow(dconsensus), ' files \n' )
      }
      sequence <- readLines(dconsensus$file[i])
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
    save(ddepth,file = infile.depth)
  }
  

  # high depth on > 250 bps
  ddepth_select <- subset(ddepth,select=c('WINDOW','LEN_LETTER','HIGH_DEPTH','PANGEA_ID'))
  ddepth_select <- ddepth_select[HIGH_DEPTH>= length_cutoff]
  ddepth_select <- unique(ddepth_select)
  # tmp = ddepth_select[,length(WINDOW), by='PANGEA_ID']
  
  ddepth_select_w <- ddepth_select[,list(PANGEA_ID=unique(PANGEA_ID)),
                                   by=c('WINDOW')]
  ddepth_select_w[,HD:=1]
  ddepth_select_w <- data.table::dcast(ddepth_select_w, PANGEA_ID~WINDOW, value.var='HD')
  
  # plot histogram of number of windows with high depth
  tmp = data.table(V=rowSums(!is.na(ddepth_select_w[,2:100])))
  ggplot(tmp,aes(V)) +
    geom_histogram() +
    labs(x=paste0('number of windows with high depth postions >= ',length_cutoff))+
    theme_bw()
  ggsave(filename = file.path(potential.networks.analysis.dir , 'plots',
                              'number_of_window_with_high_depth_position_gt250.pdf'),
         width = 6, height = 4)
  
  
  # load thresholds
  load(infile.threshold)
  threshold <- threshold$`50%`
  
  
  # summarise distance & supports per pair per window into a large matrix type
  # rows pairs
  # columns windows
  for (i in 1:windown) {
    # distance per window
    cat('\n process window ',i , '\n')
    load(infile.dist[i])
    #
    distancei <- distancei[LENGTH >= length_cutoff,]
    
    # select individuals with high depth postions >= 250
    ddepth_selecti = unique(subset(ddepth_select[WINDOW==i,],select='PANGEA_ID'))
    setnames(ddepth_selecti, 'PANGEA_ID', 'pangea_id')
    ddepth_selecti[,HD:=1]
    setnames(ddepth_selecti, colnames(ddepth_selecti), paste0(colnames(ddepth_selecti),'1'))
    distancei <- merge(distancei, ddepth_selecti, by.x='TAXA1', by.y='pangea_id1',all.x=T)
    setnames(ddepth_selecti, colnames(ddepth_selecti), gsub('1','2',colnames(ddepth_selecti)))
    distancei <- merge(distancei, ddepth_selecti, by.x='TAXA1', by.y='pangea_id2',all.x=T)
    setnames(ddepth_selecti, colnames(ddepth_selecti), gsub('2','',colnames(ddepth_selecti)))
    distancei <- distancei[HD1==1 & HD2==1,]
    
    # add threshold and close status
    distancei[, threshold_bimodal:=threshold[i]]
    distancei[, CLOSE:= as.integer(PERC > threshold_bimodal)]
    
    # average distance per pair
    tmp  <- unique(subset(distancei, select= c('TAXA1','TAXA2','PERC')))
    tmp <- tmp[, list(PERC=mean(PERC, na.rm=T)), by=c('TAXA1','TAXA2')]
    setnames(tmp, 'PERC', paste0('PERC',i))
    if( i ==1 ){
      distance <- copy(tmp)
    }else{
      distance <- merge(distance, tmp, by =c('TAXA1','TAXA2'),all=T)
    }
    gc()
    
    # number of support per pair
    tmp  <- unique(subset(distancei, select= c('TAXA1','TAXA2','CLOSE')))
    tmp <-  tmp[, list(CLOSE=mean(CLOSE,na.rm = T)), by=c('TAXA1','TAXA2')]
    setnames( tmp, 'CLOSE', paste0('CLOSE',i))
    if( i ==1 ){
      dclassif <- copy( tmp)
    }else{
      dclassif <- merge(dclassif,  tmp, by =c('TAXA1','TAXA2'),all=T)
    }
    gc()
  }
  
  # save
  save(distance,file = file.path(potential.networks.analysis.paste0('distance_overall_method', method,'_sequence_level_v4.rda')))
  save(dclassif,file = file.path(potential.networks.analysis.paste0('support_overall_method', method,'_sequence_level_v4.rda')))
}

rkuvri.make.close.pairs.and.clusters <- function(){
  #' input couples, sequences, person ids, consensus sequence info, distances and supports per sequence pair
  #' classify sequence pair, and thus individual pairs, to be close or distant, and break individuals into close clusters
  #' main output would be close pairs of individuals and clusters
  #' Note: eligible pairs further requires at least 4 windows have >=250 high-depth positions 
  window_size <- 500
  windown <- 99 
  method <- 2
  
  # file names
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  out.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_'
  potential.networks.analysis.dir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100"
  infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
  infile.sequence <- file.path(data.dir,"200422_PANGEA2_RCCSMRC_alignment.fasta")
  infile.ind.rccs <- file.path(data.dir,'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.mrc <- file.path(data.dir,'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
  infile.rccs <- file.path(data.dir,'PANGEA2_RCCS/200422_PANGEA2_RCCS_mapped_samples.rds')
  infile.mrc <- file.path(data.dir,'PANGEA2_MRC/200422_PANGEA2_MRCUVRI_mapped_samples.rds')
  infile.dist <- file.path(potential.networks.analysis.dir, paste0('subsequence_window_',1:windown,'results.rda'))
  infile.support <- file.path(potential.networks.analysis.dir ,paste0('support_overall_method', method,'_sequence_level_processed_v4.rda'))

  # load
  load(infile.support)
  nrow(df)
  
  # enough overlapping windows with high depths
  df = df[DATA >=4,]
  cat(nrow(df),' pairs of sequences left after requiring at least 4 overlapping windows with data \n')
  
  # take close pairs
  df =df[PROP >= 0.8]
  cat(nrow(df),' pairs of sequences left after requiring 80% windows indicating close \n' )
  
  # add study id
  id.dt <- data.table(read.csv(infile.ind.rccs))
  id.dt <- subset(id.dt,select = c("pt_id","sex","pangea_id"))
  id.dt[,pangea_id:=paste0('RCCS_',pangea_id)]
  tmp <- data.table(read.csv(infile.ind.mrc))
  tmp <- subset(tmp,select = c("pt_id","sex","pangea_id"))
  tmp[,pangea_id:=paste0('MRCUVRI_',pangea_id)]
  id.dt <- rbind(id.dt,tmp)
  id.dt <- unique(id.dt)
  setnames(id.dt, colnames(id.dt), paste0(colnames(id.dt),'1'))
  df <- merge(df,id.dt, by.x=c('TAXA1'), by.y=c('pangea_id1'), all.x=T)
  setnames(id.dt, colnames(id.dt), gsub('1','2',colnames(id.dt)))
  df <- merge(df,id.dt, by.x=c('TAXA2'), by.y=c('pangea_id2'), all.x=T)
  setnames(id.dt, colnames(id.dt), gsub('2','',colnames(id.dt)))
  
  # remove sequence with unknown study ids
  df[is.na(pt_id1) | is.na(pt_id2)]
  df = df[!is.na(pt_id1) & !is.na(pt_id2)]
  
  # take close pairs and remove those composed of one individual
  close_pairs <- unique(subset(df,select=c('pt_id1', 'pt_id2')))
  close_pairs <- close_pairs[pt_id1!=pt_id2]
  close_pairs[,pt_id1:= as.character(pt_id1)]
  close_pairs[,pt_id2:= as.character(pt_id2)]
  
  # save close pairs
  cat('Write close pairs to ',file.path(potential.networks.analysis.dir,'close_pairs.rda'),'...\n')
  save(close_pairs,file=file.path(potential.networks.analysis.dir,'close_pairs.rda'))
  
  # load(file.path(potential.networks.analysis.dir,'close_pairs.rda'))
  # make graph
  chains<- graph.data.frame(close_pairs, directed=FALSE, vertices=NULL)
  rtc	<- data.table(ID=V(chains)$name, CLU=clusters(chains, mode="weak")$membership)
  
  # make clusters 
  tmp	<- rtc[, list(CLU_SIZE=length(ID)), by='CLU']
  setkey(tmp, CLU_SIZE)
  print(tail(tmp))
  tmp[, IDCLU:=rev(seq_len(nrow(tmp)))]
  rtc	<- subset( merge(rtc, tmp, by='CLU'))
  rtc[, CLU:=NULL]
  setkey(rtc, IDCLU)
  rtc
  
  # take individuals from the largest cluster - 291
  id = rtc[IDCLU==1,]$ID
  chains<- subset(close_pairs, select=c(pt_id1, pt_id2))
  chains<- chains[pt_id1 %in% id & pt_id2 %in% id ]
  # id_largest_clu = unique(c(chains$pt_id1,chains$pt_id2))
  
  
  # distance from consensus sequences
  mapping_rccs <- readRDS(infile.rccs)
  dconsensus <- subset(mapping_rccs,select=c('PANGEA_ID','CONSENSUS','F'))
  dconsensus <- unique(dconsensus[CONSENSUS!="",])
  dconsensus[,PANGEA_ID:=paste0('RCCS_',PANGEA_ID)]
  mapping_mrc <- readRDS(infile.mrc)
  tmp <- subset(mapping_mrc,select=c('PANGEA_ID','CONSENSUS','F'))
  tmp <- unique(tmp[CONSENSUS!="",])
  tmp[,PANGEA_ID:=paste0('MRCUVRI_',PANGEA_ID)]
  dconsensus = rbind(dconsensus, tmp)
  dconsensus = subset(dconsensus, select=c('PANGEA_ID','F'))
  setnames(dconsensus, c('PANGEA_ID','F'), c('pangea_id','file_name'))
  dconsensus <- merge(dconsensus, id.dt, by='pangea_id', all.x=T)
  tmp <-  data.table(pt_id=id)
  dconsensus <- merge(dconsensus, tmp, by='pt_id', all.y=T)
  dconsensus <- subset(dconsensus,select=c('pt_id','file_name'))
  dconsensus <- unique(dconsensus)
  dconsensus[,sum(is.na(file_name))]
  lconsensus <- list()
  for (i in 1:nrow(dconsensus)) {
    lconsensus[[as.character(dconsensus$pt_id[i])]] = read.fasta(dconsensus$file_name[i])
  }
  
  
  ddist_copy		<- chains[,{ 
    seq1 <- lconsensus[[as.character(pt_id1)]][[1]]
    seq2 <- lconsensus[[as.character(pt_id2)]][[1]]
    pos <- which(seq1 !='-' & seq1 !='?' & seq1 !='n' & seq2 !='-' & seq2 !='?' &  seq2 !='n')
    if (length(pos)==0 ){
      len <- as.integer(-1.0)
      match <- -1.0
    }else{
      seq1 <- seq1[pos]
      seq2 <- seq2[pos]
      len <- length(pos)
      # cat(TAXA1,'\n',seq1,'\n',TAXA2,'\n',seq2,'\n')
      # 
      match <- sum(sapply(1:length(pos),function(pos){is.match(seq1[pos],seq2[pos])}))
    }
    list(LENGTH= len,
         MATCHES= match)
  }, by=c('pt_id1','pt_id2')]
  
  ddist_copy[,SIMILARITY:=MATCHES/LENGTH]
  
  # pairs from the largest cluster
  chains_pairs <- copy(chains)
  
  # break the largest cluster into communities
  ddist = copy(ddist_copy)
  ddist[SIMILARITY>=0.975, SIMILARITY:=1]
  ddist <- ddist[,list(SIMILARITY = mean(SIMILARITY)),by=c('pt_id1','pt_id2')]
  chains<- graph.data.frame(ddist, directed=FALSE, vertices=NULL)
  E(chains)$weight = ddist$SIMILARITY
  chains <- simplify(chains)
  E(chains)$weight = ifelse(E(chains)$weight > 1,E(chains)$weight/2, E(chains)$weight)
  # clusters(chains, mode='weak')
  comm <-cluster_louvain(chains)
  sapply(1:max(comm$membership),function(x)sum(comm$membership==x))
  # 4 18 10 45 11 15 50 26 38 17 16 22  5 14
  df = data.table(PT_ID=comm$names,MEMBERSHIP = comm$membership)
  
  #
  df_graph = as_long_data_frame(chains)
  df_graph = data.table(df_graph[,3:5])
  colnames(df_graph) = c('SIMILARITY','pt_id1','pt_id2')
  
  
  # take bridging individuals 
  take_bridging_individual <- function(df_graph,df){
    bridging = list()
    for (i in 1:max(df$MEMBERSHIP)) {
      tmp1 = df[MEMBERSHIP==i,]$PT_ID
      tmp2 = df_graph[pt_id1 %in% tmp1 & !(pt_id2 %in% tmp1) & SIMILARITY >=0.975,]
      tmp3 = df_graph[pt_id2 %in% tmp1 & !(pt_id1 %in% tmp1) & SIMILARITY >=0.975,]
      bridging[[i]] = unique(c(as.character(tmp2$pt_id1),  
                               as.character(tmp3$pt_id2)))
    }
    return(bridging)
  }
  #
  bridging <- take_bridging_individual(df_graph,df)  
  bridging_all <- unique(do.call(c,bridging))
  length(bridging_all)
  # 95
  
  # break bridging individuals into clusters
  id_clu <- max(df$MEMBERSHIP) 
  ddist <- ddist[pt_id1 %in% bridging_all & pt_id2 %in% bridging_all]
  chains<- graph.data.frame(ddist, directed=FALSE, vertices=NULL)
  E(chains)$weight = ddist$SIMILARITY
  chains <- simplify(chains)
  E(chains)$weight = ifelse(E(chains)$weight > 1,E(chains)$weight/2, E(chains)$weight)
  tmp = clusters(chains, mode='weak')
  
  # break bridging individuals into communities
  pt_id_tmp = names(tmp$membership)
  ddist_tmp <- ddist[pt_id1 %in% pt_id_tmp & pt_id2 %in% pt_id_tmp]
  chains<- graph.data.frame(ddist_tmp, directed=FALSE, vertices=NULL)
  E(chains)$weight = ddist_tmp$SIMILARITY
  chains <- simplify(chains)
  E(chains)$weight = ifelse(E(chains)$weight > 1,E(chains)$weight/2, E(chains)$weight)
  comm <-cluster_louvain(chains)
  sapply(1:max(comm$membership),function(x)sum(comm$membership==x))
  #  9  5 13 10 21 37
  
  # communities of bridging individuals as clusters
  df_tmp = data.table(PT_ID=comm$names,MEMBERSHIP = comm$membership + id_clu)
  df = rbind(df, df_tmp)
  df_graph_tmp = as_long_data_frame(chains)
  df_graph_tmp = data.table(df_graph_tmp[,3:5])
  colnames(df_graph_tmp) = c('SIMILARITY','pt_id1','pt_id2')
  
  # take bridging individuals of communities of bridging individuals and consider it as a separate cluster
  bridging <- take_bridging_individual(df_graph_tmp,df_tmp)  
  bridging_all <- unique(do.call(c,bridging))
  id_clu <- max(df$MEMBERSHIP) + 1
  df = rbind(df, data.table(PT_ID=bridging_all,MEMBERSHIP = id_clu))
  tmp = df[,list(CLU_SIZE = length(PT_ID)),by='MEMBERSHIP']
  df = merge(df, tmp,by='MEMBERSHIP')
  setnames(df, c('MEMBERSHIP','PT_ID'), c('IDCLU','ID'))
  
  # merge with rtc
  rtc = rtc[IDCLU!=1,]
  rtc[,IDCLU:=IDCLU + max(df$IDCLU)-1]
  df = rbind(rtc,df)
  
  cat('Write clusters to ',file.path(potential.networks.analysis.dir,'clusters.rda'),'...\n')
  save(df,file = file.path(potential.networks.analysis.dir,'clusters.rda'))
  
  length(unique(df$ID))
  tmp = unique(subset(df,select=c('IDCLU','CLU_SIZE')))
  setkey(tmp, IDCLU)
  
  nrow(close_pairs)
  length(unique(c(unique(close_pairs$pt_id1),unique(close_pairs$pt_id2))))
  # 2525 pairs
  # 808 networks
  
  ggplot(tmp,
         aes(x = reorder(IDCLU, -CLU_SIZE), y = CLU_SIZE,label=CLU_SIZE)) +
    geom_bar(stat = "identity",width = 0.1)+
    labs(x='cluster ID', y='cluster size') +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  # save histogram of cluster sizes
  ggsave(file.path(potential.networks.analysis.dir,'plots','cluster_size_plot.pdf'),width = 10, height = 4)
}


rkuvri.validate.classification <- function(){
  #' input tree distances, hxb2 to map, sequences, average distance per pair (or thresholds and distances per window to calculate it)
  #' input depths, couples, person ids, close pairs, distances per pair
  #' not necessary for making clusters
  #' validates the classification and thresholds
  #' 
  window_size <- 500
  windown <- 99
  length_cutoff <- 250
  method <- 2
  
  # file names
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  out.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_'
  potential.networks.analysis.dir <- "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100"
  infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
  infile.sequence <- file.path(data.dir,"200422_PANGEA2_RCCSMRC_alignment.fasta")
  infile.ind.rccs <- file.path(data.dir,'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.mrc <- file.path(data.dir,'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
  infile.tree.distance <- file.path(data.dir ,'HIV_DistanceNormalisationOverGenome.csv')
  infile.average.distance <- file.path(potential.networks.analysis.dir,'average_distance.rds')
  infile.hxb2 <- file.path(data.dir ,'HIV1_COM_2015_genome_DNA_ALL.fasta')
  infile.close.pairs <- file.path(potential.networks.analysis.dir,'close_pairs.rda')
  infile.dist <- file.path(potential.networks.analysis.dir, paste0('subsequence_window_',1:windown,'results.rda'))
  infile.depth <- file.path(data.dir,'210121_RCCSMRC_depth.rda')
  if(method == 1){
    infile.threshold <- file.path(potential.networks.analysis.dir,'threshold_dip_indep_exp5.rda')
  }else if(method == 2){
    infile.threshold <- file.path(potential.networks.analysis.dir,'threshold_5quantile_indep_exp5.rda')
  }else if(method == 3){
    infile.threshold <- file.path(potential.networks.analysis.dir,'threshold_median_indep_exp5.rda')
  }

  
  ### plot threshold vs distances per window
  # load tree-based distance
  tree_distance <- read.csv(infile.tree.distance)
  tree_distance <- data.table(tree_distance)
  
  # map nongap positions to the hxb2
  all = read.fasta(infile.hxb2)
  all_hxb2 <- all$B.FR.83.HXB2_LAI_IIIB_BRU.K03455
  all_hxb2_nongap <- which(!grepl('-',all$B.FR.83.HXB2_LAI_IIIB_BRU.K03455))
  tree_distance[,POSITION_WRT_HXB2_WITH_GAP:=all_hxb2_nongap[tree_distance$POSITION_WRT_HXB2]]
  
  # add missing locations to the hxb2
  tree_distance <- rbind(tree_distance,
                         data.table(POSITION_WRT_HXB2_WITH_GAP=setdiff(1:max(tree_distance$POSITION_WRT_HXB2_WITH_GAP),tree_distance$POSITION_WRT_HXB2_WITH_GAP),
                                    MEDIAN_PAIRWISE_DISTANCE_BETWEEN_STANDARD_REFS=NA_real_,
                                    POSITION_WRT_HXB2=NA_real_))
  setkey(tree_distance,POSITION_WRT_HXB2_WITH_GAP)
  
  # load data
  alignment <- read.fasta(file = infile.sequence)
  nsequence <- length(alignment)
  npos <- unique(lengths(alignment))
  
  # windows
  windows_last_start <- ceiling(npos/100) * 100 - window_size + 1
  windows_2last_end <- floor(npos/100) * 100
  windows_start <- seq(1,windows_last_start,100)
  windows_end <- c(seq(window_size ,windows_2last_end,100),npos)
  windows <- seq_len(length(windows_start))
  
  # mean tree_distance per window
  tree_distance_per_window <- vector()
  for (w in windows){
    tree_distance_per_window <- c( tree_distance_per_window,
                                   mean( tree_distance$MEDIAN_PAIRWISE_DISTANCE_BETWEEN_STANDARD_REFS[windows_start[w]:windows_end[w]], na.rm = TRUE ) )
  }
  
  # average distance per window
  if(file.exists(infile.average.distance)){
    ans <- as.data.table(readRDS(infile.average.distance))
  }else{
    load(infile.distance)
    average_distance <-  distance%>%
      summarise_if(is.numeric, mean, na.rm = TRUE)
    
    load(infile.threshold)
    ans <- data.table(THRES=threshold$`50%`,
                      DISTM=matrix(average_distance[1:windown],ncol=1))
    saveRDS(ans,file = file.path(potential.networks.analysis.dir,'average_distance.rds'))
  }

  
  # compare distance and threshold
  ans[, DISTT:=tree_distance_per_window]
  colnames(ans) <- c('threshold','average distance','tree distance')
  ans[, window:=seq_len(nrow(ans))]
  
  ans <- melt(ans, id.var=c('window','threshold'))
  
  # plot threshold against distances
  ggplot(ans, aes(as.numeric(threshold), as.numeric(value), label=factor(window), color=factor(variable)))+
    geom_point()+
    scale_x_continuous(expand = c(0.05,0.05))   +
    scale_y_continuous(expand = c(0.05,0.05)) +
    theme_bw() +
    labs(x='threshold',
         y='distance',
         color='type') +
    scale_color_manual(values=c('tree distance'='royalblue3','average distance'='006400')) +
    theme(legend.position='bottom')
  ggsave(file.path(potential.networks.analysis.dir,'plots',paste0('distance_vs_threshold_method', method,'.pdf')),width = 6, height = 6)
  
  
  ### check the corresponding tree distance of thresholds in pol 
  # create a function mapping distance to tree distance
  y=seq(0.01,0.05,0.005)
  x=c(0.022,0.031,0.036,0.041,0.044,0.048,0.051,0.054,0.056)
  fit2 <- lm(y~poly(x,4,raw=TRUE))
  # check
  xx=x
  cat('deviations from y: ',predict(fit2, data.frame(x=xx))-y,'\n')
  # function
  coefs = coef(fit2)
  f = function(newdist, coefs){
    res <- coefs[1] + (coefs[2] * newdist) + (coefs[3] * newdist^2) + (coefs[4] * newdist^3)+ (coefs[5] * newdist^4)
    return(res)
  }
  
  load(infile.threshold)
  thres = 1-threshold$`50%`
  thres_map = sapply(thres,function(x)f(x,coefs))
  
  # identify pol
  start_pol <- min(all_hxb2_nongap[2085:5096])
  end_pol <- max(all_hxb2_nongap[2085:5096])
  start_pol_win <-which( windows_start -start_pol > 0)[1]
  end_pol_win <-which(windows_end - end_pol> 0)[1] - 1
  quantile(thres[start_pol_win:end_pol_win],probs=c(0.025,0.5,0.975))
  # distance
  # 0.01615810 0.02545088 0.02954131 
  quantile(thres_map[start_pol_win:end_pol_win],probs=c(0.025,0.5,0.975))
  # tree distance
  # 0.009986778 0.011284782 0.013984869
  
  
  # load depth
  load(infile.depth)
  
  ddepth_select <- subset(ddepth,select=c('WINDOW','LEN_LETTER','HIGH_DEPTH','PANGEA_ID'))
  ddepth_select <- ddepth_select[HIGH_DEPTH>= length_cutoff]
  ddepth_select <- unique(ddepth_select)
  # tmp = ddepth_select[,length(WINDOW), by='PANGEA_ID']
  
  ddepth_select_w <- ddepth_select[,list(PANGEA_ID=unique(PANGEA_ID)),
                                   by=c('WINDOW')]
  ddepth_select_w[,HD:=1]
  ddepth_select_w <- data.table::dcast(ddepth_select_w, PANGEA_ID~WINDOW, value.var='HD')
  
  
  # map alignments to studyid
  dinfo <- data.table(pangea_id=names(alignment))
  id.dt <- data.table(read.csv(infile.ind.rccs))
  id.dt <- subset(id.dt,select = c("pt_id","sex","pangea_id"))
  id.dt[,pangea_id:=paste0('RCCS_',pangea_id)]
  tmp <- data.table(read.csv(infile.ind.mrc))
  tmp <- subset(tmp,select = c("pt_id","sex","pangea_id"))
  tmp[,pangea_id:=paste0('MRCUVRI_',pangea_id)]
  id.dt <- rbind(id.dt,tmp)
  id.dt <- unique(id.dt)
  dinfo <- merge(dinfo, id.dt, by="pangea_id", all.x=T)
  tmp <- dinfo[is.na(pt_id)]
  dinfo <- dinfo[!is.na(pt_id)]
  tmp2 <- subset(rccsData,select=c('RCCS_studyid','Pangea.id','SEX'))
  colnames(tmp2) <- c("pt_id","pangea_id","sex")
  set(tmp, NULL, c("pt_id","sex"), NULL)
  tmp <- merge(tmp, tmp2, by='pangea_id',all.x=T)
  dinfo <- rbind(dinfo, tmp)
  cat('No personal information found for ',nrow(dinfo[is.na(pt_id)]), ' sequences \n')
  
  # process couple, map to studyid
  load(infile.couple)
  couple <- unique(subset(data.table(coupdat), select=c('male.RCCS_studyid', 'female.RCCS_studyid')))
  couple[, COUPLE:=1]
  couple[,male.RCCS_studyid:=paste0('RK-',male.RCCS_studyid)]
  couple[,female.RCCS_studyid:=paste0('RK-',female.RCCS_studyid)]
  
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
  
  
  # save similarities by couples and classifications
  df_window_summary <- data.table()
  for (i in seq_len(windown)){
    cat('\n process window ',i , '\n')
    load(infile.dist[i])
  
    # select individuals with high depth postions >= 250
    distancei <- distancei[LENGTH >= 250,]
    ddepth_selecti = unique(subset(ddepth_select[WINDOW==i,],select='PANGEA_ID'))
    setnames(ddepth_selecti, 'PANGEA_ID', 'pangea_id')
    ddepth_selecti[,HD:=1]
    setnames(ddepth_selecti, colnames(ddepth_selecti), paste0(colnames(ddepth_selecti),'1'))
    distancei <- merge(distancei, ddepth_selecti, by.x='TAXA1', by.y='pangea_id1',all.x=T)
    setnames(ddepth_selecti, colnames(ddepth_selecti), gsub('1','2',colnames(ddepth_selecti)))
    distancei <- merge(distancei, ddepth_selecti, by.x='TAXA1', by.y='pangea_id2',all.x=T)
    setnames(ddepth_selecti, colnames(ddepth_selecti), gsub('2','',colnames(ddepth_selecti)))
    distancei <- distancei[HD1==1 & HD2==1,]
    
    # add pt_id
    setnames(dinfo, colnames(dinfo), paste0(colnames(dinfo),'1'))
    distancei <- merge(distancei,dinfo, by.x=c('TAXA1'), by.y=c('pangea_id1'), all.x=T) 
    setnames(dinfo, colnames(dinfo), gsub('1','2',colnames(dinfo)))
    distancei <- merge(distancei,dinfo, by.x=c('TAXA2'), by.y=c('pangea_id2'), all.x=T) 
    setnames(dinfo, colnames(dinfo), gsub('2','',colnames(dinfo)))
    
    # remove no pt_id
    distancei <- distancei[!is.na(pt_id1) & !is.na(pt_id2) & pt_id1!=pt_id2,]
    
    # add couple status
    distancei <- merge(distancei,tmp_couple, by=c('pt_id1','pt_id2'), all.x=T)
    distancei[is.na(COUPLE),COUPLE:=0]
    
    # add close
    load(infile.close.pairs)
    close_pairs <- subset(close_pairs, select=c('pt_id1','pt_id2'))
    close_pairs[,CLOSE:=1]
    distancei <- merge(distancei, close_pairs, by=c('pt_id1','pt_id2'), all.x=T)
    distancei[is.na(CLOSE),CLOSE:=0]
    # summarise by close and couple per window
    tmp <- distancei[,list(PREC=quantile(PERC, probs=c(0.05, 0.25, 0.5, 0.75, 0.95),na.rm=T),
                           PROB=c(0.05, 0.25, 0.5, 0.75, 0.95)),by=c('COUPLE','WINDOW','CLOSE')]
    tmp <- dcast(tmp,COUPLE+CLOSE+WINDOW~PROB,value.var = 'PREC')
    df_window_summary <- rbind(df_window_summary, tmp)
    
    gc()
  }
  
  threshold <- data.table(WINDOW=1:windown,
                          THRESHOLD=threshold,
                          POS=windows_start)
  
  df_window_summary <- merge(df_window_summary, threshold, by='WINDOW')
  
  # save 
  saveRDS(df_window_summary, file=file.path(potential.networks.analysis.dir,paste0('distance_distribution_perwindow_bycoupleclose_method', method,'.rds')))
  
  
  # df_window_summary <- readRDS(file.path(potential.networks.analysis.dir,paste0('distance_distribution_perwindow_bycoupleclose_method', method,'.rds')))
  # plot similarity distributions over windows by couple and close and thresholds
  
  ggplot(df_window_summary)+
    geom_step(aes(x=POS,y=`0.5`,col=factor(CLOSE))) +
    geom_step(aes(x=POS,y=`THRESHOLD.50%`),color='black')+
    guides(col=F)+
    geom_stepribbon(aes(x=POS,ymin=`0.25`, ymax=`0.75`,fill=factor(CLOSE)),alpha=0.6)+
    geom_stepribbon(aes(x=POS,ymin=`0.05`, ymax=`0.95`,fill=factor(CLOSE)),alpha=0.2)+
    scale_fill_manual(values = c('1'='#C70039',
                                 '0'='#2E86C1'),
                      labels=c('1'='yes',
                               '0'='no')) +
    scale_color_manual(values = c('1'='#C70039',
                                  '0'='#2E86C1'),
                       labels=c('1'='yes',
                                '0'='no')) +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = 'bottom')+
    labs(x='genome', y='distance',fill='close') +
    facet_grid(factor(COUPLE,levels = c(0,1), labels = c('not couples','couples'))~.)
  ggsave(file.path(potential.networks.analysis.dir,'plots',paste0('similarity_vs_window_method', method,'.pdf')), width = 6, height = 6)
  
  
  
  
  
  
}