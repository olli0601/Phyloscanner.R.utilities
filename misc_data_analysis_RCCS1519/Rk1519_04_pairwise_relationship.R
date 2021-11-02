delete_nonexcised <- function(){
  require(data.table)
  # files
  HOME <<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  treedir <- file.path(HOME,'210325_phsc_output/')
  df <- data.table(F=list.files(treedir,recursive = T, pattern = '*.treefile'))
  
  # find the case when the tree files exist with and without position excision
  df[,Ex:=grepl('PositionsExcised',F)]
  df[,base:=basename(F)]
  df[,dir:=dirname(F)]
  df[,wfrom:=sub("^ptyr[0-9]+_InWindow_[A-Za-z]?{16}?_?([0-9]+)_to_([0-9]+)_v2\\.treefile$","\\1",base)]
  df[,wto:=sub("^ptyr[0-9]+_InWindow_[A-Za-z]?{16}?_?([0-9]+)_to_([0-9]+)_v2\\.treefile$","\\2",base)]
  df[,run:=sub("ptyr([0-9]+)_trees","\\1",dir)]
  tmp <- dcast(df, run + wfrom + wto ~ Ex, value.var="F")
  tmp <- data.table(tmp)
  tmp <- tmp[!is.na(`TRUE`) & !is.na(`FALSE`)]
  tmp2 <- unique(subset(df,select='dir'))
  tmp2[,dir:=file.path(treedir,dir)]
  
  # delete the one without position excision
  tmp[,dir:=dirname(`TRUE`)]
  for (i in 1:nrow(tmp)){
    cat(system(paste0( 'rm ',file.path(treedir,tmp$`FALSE`[i])),intern = T))
  }
  # don't generate trees for those files next time
}


make.phyloscanner<- function()
{
  require(tidyverse)
  require(data.table)
  require(phyloscannerR)
  
  # files
  HOME <<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  treedir <- file.path(HOME,'210325_phsc_output/')
  tmpdir	<- file.path(HOME,"210325_phsc_work/")
  # # outdir	<- file.path(HOME,'210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca/') #100 max postcount 30min
  # # outdir	<- file.path(HOME,'210325_phsc_phscrelationships_025_30_min_read_null_max_read_posthoccount_im_mrca/') #null max postcount 30min
  # outdir	<- file.path(HOME,'210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount/')
  # outdir	<- file.path(HOME,'210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount/')
  outdir	<- file.path(HOME,'210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/')
  outdir	<- file.path(HOME,'210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd/')
  # for(k in 1:10){
  outdir	<- paste0(HOME,'/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_seed',k,'/')
  dir.create(outdir)
  prog.phyloscanner_analyse_trees <- '/rds/general/user/xx4515/home/phyloscanner/phyloscanner_analyse_trees.R'
  
  #	set phyloscanner variables	
  control	<- list()
  control$allow.mt <- TRUE
  # control$n.mt <- 0.5
  # control$n.mt <- 3	
  # control$n.mt <- NULL
  control$identify.multifurcation = TRUE
  # control$identify.multifurcation = FALSE
  control$alignment.file.directory = NULL 
  control$alignment.file.regex = NULL
  control$blacklist.underrepresented = FALSE	
  control$count.reads.in.parsimony = TRUE
  # control$distance.threshold <- '0.025'
  control$distance.threshold <- '0.02 0.05'
  control$do.dual.blacklisting = FALSE					
  control$duplicate.file.directory = NULL
  control$duplicate.file.regex = NULL
  control$file.name.regex = "^(?:.*\\D)?([0-9]+)_to_([0-9]+).*$"
  control$guess.multifurcation.threshold = FALSE
  control$max.reads.per.host <- NULL
  # control$min.reads.per.host <- NULL
  # control$max.reads.per.host <- 100
  control$min.reads.per.host <- 30
  control$min.tips.per.host <- 1	
  control$multifurcation.threshold = 1e-5
  control$multinomial= TRUE
  control$norm.constants = NULL
  control$norm.ref.file.name = "~/normalisation_ByPosition.csv"
  control$norm.standardise.gag.pol = TRUE
  control$no.progress.bars = TRUE
  control$outgroup.name = "REF_CON_H"
  control$output.dir = outdir
  control$parsimony.blacklist.k = 15
  control$prune.blacklist = FALSE
  control$post.hoc.count.blacklisting <- TRUE
  # control$post.hoc.count.blacklisting <- NULL
  control$ratio.blacklist.threshold = 0.01
  control$raw.blacklist.threshold = 3			
  control$recombination.file.directory = NULL
  control$recombination.file.regex = NULL
  control$relaxed.ancestry = TRUE
  control$sankoff.k = 15
  control$sankoff.unassigned.switch.threshold = 0
  control$seed = 42
  # control$seed = k
  control$splits.rule = 's'
  control$tip.regex = '^(.*)-fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'
  control$tree.file.regex = "^(.*)\\.treefile$" # from rscript
  control$treeFileExtension = '.treefile'
  control$use.ff = FALSE
  control$user.blacklist.directory = NULL 
  control$user.blacklist.file.regex = NULL
  control$verbosity = 1	
  
  
  #	make bash for many files	
  df <- tibble(F=list.files(treedir,pattern = 'ptyr*'))
  # ,pattern = '*.treefile$')
  df <- df %>%
    mutate(TYPE:= gsub('ptyr([0-9]+)_(.*)','\\2', F),
           RUN:= as.integer(gsub('ptyr([0-9]+)_(.*)','\\1', F))) %>%
    # mutate(TYPE:= gsub('^([^\\.]+)\\.[a-z]+$','\\1',TYPE)) %>%
    mutate(TYPE:= gsub('^[^\\.]+\\.([a-z]+)$','\\1',TYPE)) %>%
    spread(TYPE, F) %>%
    set_names(~ str_to_upper(.))
  
  valid.input.args <- cmd.phyloscanner.analyse.trees.valid.args(prog.phyloscanner_analyse_trees)
  cmds <- vector('list',nrow(df))
  
  for(i in seq_len(nrow(df)))
  {
    #	set input args
    control$output.string <- paste0('ptyr',df$RUN[i])
    #	make script
    tree.input <- file.path(treedir, df$TREES[i])
    cmd <- cmd.phyloscanner.analyse.trees(prog.phyloscanner_analyse_trees, 
                                          tree.input, 
                                          control,
                                          valid.input.args=valid.input.args)
    cmds[[i]] <- cmd		
  }	
  cat(cmds[[1]])
  
  #
  # 	submit array job to HPC
  #
  #	make header
  hpc.load			<- "module load anaconda3/personal \n source activate phylo"	# make third party requirements available	 
  hpc.select			<- 1						# number of nodes
  hpc.nproc			<- 1						# number of processors on node
  hpc.walltime		<- 23						# walltime
  # hpc.walltime		<- 923						# walltime
  if(1)		
  {
    hpc.q			<- NA						# PBS queue
    hpc.mem			<- "36gb" 					# RAM		
  }
  #		or run this block to submit a job array to Oliver's machines
  # if(1)
  # {
  #   hpc.q			<- "pqeelab"				# PBS queue
  #   hpc.mem			<- "12gb" 					# RAM		
  # }
  hpc.array			<- length(cmds)	# number of runs for job array	
  pbshead		<- "#!/bin/sh"
  tmp			<- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
  pbshead		<- paste(pbshead, tmp, sep = "\n")
  tmp			<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
  pbshead 	<- paste(pbshead, tmp, sep = "\n")
  pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")	
  if(!is.na(hpc.array))
    pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')	
  if(!is.na(hpc.q)) 
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  pbshead 	<- paste(pbshead, hpc.load, sep = "\n")	
  cat(pbshead)
  #	make array job
  for(i in 1:length(cmds))
    cmds[[i]]<- paste0(i,')\n',cmds[[i]],';;\n')
  cmd		<- paste0('case $PBS_ARRAY_INDEX in\n',paste0(cmds, collapse=''),'esac')	
  cmd		<- paste(pbshead,cmd,sep='\n')	
  #	submit job
  outfile		<- gsub(':','',paste("phsc",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
  outfile		<- file.path(tmpdir, outfile)
  cat(cmd, file=outfile)
  cmd 		<- paste("qsub", outfile)
  cat(cmd)
  cat(system(cmd, intern= TRUE))
  # }
}


phsc.transmission.networks<- function()
{
  
  make_networks <- function(job_tag)
    {  
      library(data.table)
      library(tidyverse)
      # source('~/transmission_network_functions_phyloscanner.R')
      # optional: meta data
      indir.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships'
      indir	<- paste0(indir.base,job_tag)
      infiles	<- data.table(F=list.files(indir, pattern='*workspace.rda$', full.names=TRUE))
      infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(F)))]
      setkey(infiles, PTY_RUN)
      
      #	prepare all dwin and dc
      dca	<- infiles[, {
        cat(PTY_RUN,'\n')
        load(F)
        dc			
      }, by='PTY_RUN']
      dwina <- infiles[, {
        cat(PTY_RUN,'\n')
        load(F)
        dwin			
      }, by='PTY_RUN']
      save(dca,dwina,file = paste0('~/dcdwina',job_tag,'.rda'))
      
      control <- list(linked.group='close.and.adjacent.cat',linked.no='not.close.or.nonadjacent',linked.yes='close.and.adjacent', conf.cut=0.6, neff.cut=3,weight.complex.or.no.ancestry=0.5)
      indir.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships'
      indir	<- paste0(indir.base,job_tag)
      
      library(tidyverse)
      library(glue)
      library(igraph)
      library(RBGL)
      library(data.table)
      library(phyloscannerR)
      load(paste0('~/dcdwina',job_tag,'.rda'))
      # dwina[host.1=='AID0812' & host.2=='AID1868' & grepl('1100_',tree.id)]
      # control
      dca[,CNTRL1:=FALSE]
      dca[,CNTRL2:=FALSE]
      dca[grepl('CNTRL-',host.1),CNTRL1:=TRUE]
      dca[grepl('CNTRL-',host.2),CNTRL2:=TRUE]
      dca[grepl('CNTRL-',host.1),host.1:=gsub('CNTRL-','',host.1)]
      dca[grepl('CNTRL-',host.2),host.2:=gsub('CNTRL-','',host.2)]
      dwina[,CNTRL1:=FALSE]
      dwina[,CNTRL2:=FALSE]
      dwina[grepl('CNTRL-',host.1),CNTRL1:=TRUE]
      dwina[grepl('CNTRL-',host.2),CNTRL2:=TRUE]
      dwina[grepl('CNTRL-',host.1),host.1:=gsub('CNTRL-','',host.1)]
      dwina[grepl('CNTRL-',host.2),host.2:=gsub('CNTRL-','',host.2)]
      
      # sort dwin
      tmp			<- subset(dwina, host.1>host.2)
      setnames(tmp, c('host.1','host.2','paths12','paths21','nodes1','nodes2','CNTRL1','CNTRL2'),
               c('host.2','host.1','paths21','paths12','nodes2','nodes1','CNTRL2','CNTRL1'))
      set(tmp, NULL, 'close.and.contiguous.and.directed.cat',
          tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',close.and.contiguous.and.directed.cat)))])
      set(tmp, NULL, 'close.and.adjacent.and.directed.cat',
          tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',close.and.adjacent.and.directed.cat)))])
      set(tmp, NULL, 'close.and.contiguous.and.ancestry.cat',
          tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',close.and.contiguous.and.ancestry.cat)))])
      set(tmp, NULL, 'close.and.adjacent.and.ancestry.cat',
          tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',close.and.adjacent.and.ancestry.cat)))])
      dwina		<- rbind(subset(dwina, !(host.1>host.2)), tmp)
      
      # sort dca
      tmp			<- subset(dca, host.1>host.2)
      setnames(tmp, c('host.1','host.2','CNTRL1','CNTRL2'),
               c('host.2','host.1','CNTRL2','CNTRL1'))
      set(tmp, NULL, 'type',
          tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',type)))])
      dca		<- rbind(subset(dca, !(host.1>host.2)), tmp)  
      
      tmp <- unique(subset(dwina,select=c('PTY_RUN','host.1','host.2')))
      # tmp[host.1=='AID0812' & host.2=='AID1868']
      # dwina[host.1=='AID0812' & host.2=='AID1868' & grepl('1100_',tree.id)]
      tmp <- tmp[,list(PTY_RUN=PTY_RUN[1]),by=c('host.1','host.2')]
      dwina <- merge(dwina,tmp, by=c('host.1','host.2'))
      dca <- merge(dca,tmp, by=c('host.1','host.2','PTY_RUN'))
      
      dwina$PTY_RUN.y=NULL
      dca$PTY_RUN.y=NULL
      setnames(dwina,'PTY_RUN.x','PTY_RUN',skip_absent=T)
      setnames(dca,'PTY_RUN.x','PTY_RUN',skip_absent=T)
      
      # find pairs
      tmp <- find.pairs.in.networks(dwina, dca, control=control, verbose=TRUE)
      dpl <- copy(tmp$network.pairs)
      dc <- copy(tmp$relationship.counts)
      dw <- copy(tmp$windows)
      save(dpl, dc, dw, file=file.path(indir,'Rakai_phscnetworks_allpairs.rda'))
      
      # find networks
      tmp <- find.networks(dc, control=control, verbose=TRUE)
      dnet <- copy(tmp$transmission.networks)
      dchain <- copy(tmp$most.likely.transmission.chains)
      save(dpl, dc, dw, dnet, dchain, file=file.path(indir,'Rakai_phscnetworks.rda'))
    }
    job_tag <- '_02_05_30_min_read_100_max_read_posthoccount_im_mrca'
    job_tag <- '_02_05_30_min_read_null_max_read_posthoccount_im_mrca'
    job_tag <- '_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd'
    job_tag <- '_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd'
    job_tag <- '_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd_seed1'#9,10,1
    job_tag <- '_02_05_30_min_read_100_max_read_posthoccount_im_mrca_seed10'
    job_tag <- '_02_05_30_min_read_null_max_read_posthoccount'
    job_tag <- '_02_05_30_min_read_100_max_read_posthoccount'
    job_tag <- '_025_30_min_read_null_max_read_posthoccount_im_mrca'
    job_tag <- '_02_05_30_min_read_100_max_read_posthoccount_seed10'
    make_networks(job_tag)  

    
}



clinical_data <- function(){
  library(data.table)
  infile.sero <- '~/serodata19.rda'
  # infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
  
  # 19
  infile.aid19 <- '~/todi_pairs_171122_cl25_d50_prior23_min30_anonymised_RIDs.csv'
  infile.net19 <- '~/Rakai_phscnetworks_190706.rda'
  
  sero_analysis <- function(infile.net,infile.aid,infile.sero,analysis='2021'){
    # couple
    # load(infile.couple)
    # couple <- data.table(unique(subset(coupdat,select=c('male.RCCS_studyid', 'female.RCCS_studyid'))))
    # colnames(couple) <- c('MALE_RID','FEMALE_RID')
    # couple[,COUPLE:=1]
    # couple[,MALE_RID:=paste0('RK-',MALE_RID)]
    # couple[,FEMALE_RID:=paste0('RK-',FEMALE_RID)]
    
    load(infile.sero)
    
    # 
    aid <- data.table(read.csv(infile.aid))
    # consistent names
    if(analysis=='2021'){
      setnames(aid,'PT_ID','ID')
    }
    if(analysis=='2019'){
      aid[,ID:=paste0('RK-',ID)]   
    }
    
    
    
    load(infile.net)
    dc <- data.table(dc)
    # unique(dc$CATEGORISATION)
    dclkl <- dc[CATEGORISATION=='close.and.adjacent.and.ancestry.cat']
    tmp <- dclkl[,list(MAX_N_EFF=max(N_EFF)),by=c('H1','H2')]
    dclkl <- merge(dclkl,tmp, by=c('H1','H2'))
    dclkl <- dclkl[MAX_N_EFF==N_EFF,]
    tmp <- dclkl[,length(N_EFF),by=c('H1','H2')]
    setkey(tmp,V1)
    dclkl <- dclkl[,head(.SD, 4),by=c('H1','H2')]
    dclkl <- data.table(dcast(dclkl, H1+H2~TYPE,value.var='K_EFF'))
    # revisit
    # colnames(dclkl)
    dclkl[, SCORE_LINK:=(complex.or.no.ancestry+`12`+`21`)/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
    # dclkl[, SCORE_DIR:=(`12`+`21`)/(complex.or.no.ancestry+`12`+`21`)]
    # dclkl[, SCORE_12:=`12`/(`12`+`21`)]
    dclkl[,SCORE_12:=`12`/(complex.or.no.ancestry+`12`+`21`)]
    dclkl[,SCORE_21:=`21`/(complex.or.no.ancestry+`12`+`21`)]
    dclkl[, LINK:=(SCORE_LINK>0.6)]
    # dclkl[LINK==T, DIR:=(SCORE_DIR>0.6)]
    # dclkl[DIR==T, DIR12:=(SCORE_12>0.5)]
    dclkl[LINK==T, DIR12:=(SCORE_12>0.6)]
    dclkl[LINK==T, DIR21:=(SCORE_21>0.6)]
    
    # sort as M F
    if(analysis=='2019'){
      dclkl[,H1S:=substr(H1,9,9)]
      dclkl[,H2S:=substr(H2,9,9)]
      tmp <- dclkl[!(H1S=='M'),]
      setnames(tmp,c('H1','H2','12','21','H1S','H2S'),
               c('H2','H1','21','12','H2S','H1S'))
      tmp[,SCORE_12:=1-SCORE_12]
      tmp[,DIR12:=1-DIR12]
      dclkl <- rbind(tmp,dclkl[(H1S=='M'),])      
      dclkl <- merge(dclkl, subset(aid,select=c('ID','AID')),
                     by.x='H1',by.y='AID',all.x=T)
      dclkl <- merge(dclkl, subset(aid,select=c('ID','AID')),
                     by.x='H2',by.y='AID',all.x=T)
      dclkl <- merge(dclkl, dsero, by.x=c('ID.x','ID.y'),by.y=c('MALE_RID', 'FEMALE_RID'),all.x=T)
    }else{
      dclkl <- merge(dclkl, subset(aid,select=c('ID','AID')),
                     by.x='H1',by.y='AID',all.x=T)
      dclkl <- merge(dclkl, subset(aid,select=c('ID','AID')),
                     by.x='H2',by.y='AID',all.x=T)
      tmp <- merge(dclkl,dsero,by.x=c('ID.x','ID.y'),by.y=c('FEMALE_RID','MALE_RID'))
      dclkl <- merge(dclkl,dsero,by.x=c('ID.x','ID.y'),by.y=c('MALE_RID','FEMALE_RID'))
      setnames(tmp,c('H1','H2','12','21','ID.x','ID.y'),
               c('H2','H1','21','12','ID.y','ID.x'))
      tmp[,SCORE_12:=1-SCORE_12]
      tmp[,DIR12:=1-DIR12]
      dclkl <- rbind(dclkl, tmp)
    }
    
    dclkl <- dclkl[!is.na(EXT_TYPE)]
    dclkl[DIR12==1,PHSC_DIR:='mf']
    dclkl[DIR12==0,PHSC_DIR:='fm']
    
    tmp	<- dclkl[, which(!is.na(PHSC_DIR))]
    set(dclkl, tmp, 'EXT_EVAL', dclkl[tmp, as.character(factor(PHSC_DIR==EXT_DIR, levels=c(TRUE,FALSE),labels=c('correct','incorrect')))])
    dclkl[is.na(DIR12),EXT_EVAL:='couple most likely a pair direction not resolved']
    dclkl[LINK==F,EXT_EVAL:='couple most likely not a pair']
    # dclkl[is.na(EXT_EVAL)]
    dclkl[, table(EXT_EVAL, EXT_TYPE, useNA='if')]
  }
  
  sero_analysis(infile.net19,infile.aid19,infile.sero,analysis='2019')
  # EXT_EVAL                                           serodisc
  # correct                                                28
  # couple most likely a pair direction not resolved        5
  # couple most likely not a pair                          13
  # incorrect                                               4
  
  # 21 minread=1 no downsampling
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/'
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021')
  
  # EXT_EVAL                                           serodisc
  # correct                                                21
  # couple most likely a pair direction not resolved        8
  # couple most likely not a pair                           8
  # incorrect                                               4
  
  # 21 minread=1 downsampling=100
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read_100_max_read/'
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021')
  # EXT_EVAL                                           serodisc
  # correct                                                22
  # couple most likely a pair direction not resolved        6
  # couple most likely not a pair                           6
  # incorrect                                               5
  
  
  # 21 minread=30 downsampling=100
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read/'
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021')
  
  # EXT_EVAL                                           serodisc
  # correct                                                22
  # couple most likely a pair direction not resolved        9
  # couple most likely not a pair                           7
  # incorrect                                               3  
  
  # 21 minread=30 downsampling=100 posthoccount
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount/'
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021')
  # correct                                                21
  # couple most likely a pair direction not resolved       10
  # couple most likely not a pair                           5
  # incorrect                                               3
  
  # 21 minread=30 downsampling=100 0.25
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_025_05_30_min_read_100_max_read/'
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021')
  # EXT_EVAL                                           serodisc
  # correct                                                22
  # couple most likely a pair direction not resolved       12
  # couple most likely not a pair                           5
  # incorrect                                               3
  # 
  # 
  # 21 minread=30 downsampling=100 rawbt=10
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_rawbt10/'
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021')
  # EXT_EVAL                                           serodisc
  # correct                                                19
  # couple most likely a pair direction not resolved       13
  # couple most likely not a pair                           8
  # incorrect                                               3
  
  
  # 2019
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_2019/'
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021')
  
  # correct                                                19
  # couple most likely a pair direction not resolved       13
  # couple most likely not a pair                           8
  # incorrect                                               3
  # 
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca/'
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021')
  # EXT_TYPE
  # EXT_EVAL                                           serodisc
  # correct                                                14
  # couple most likely a pair direction not resolved       20
  # couple most likely not a pair                           5
  
  
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca/'
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021')
  # EXT_TYPE
  # EXT_EVAL                                           serodisc
  # correct                                                16
  # couple most likely a pair direction not resolved       18
  # couple most likely not a pair                           5
  # incorrect                                               1  

  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_maxread_posthoccount/'
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021')
  # correct                                                19
  # couple most likely a pair direction not resolved       13
  # couple most likely not a pair                           6
  # incorrect                                               2  
19+13+2
16+18+1
indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_025_30_min_read_null_max_read_posthoccount_im_mrca/'
infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021')
# EXT_EVAL                                           serodisc
# correct                                                14
# couple most likely a pair direction not resolved       22
# couple most likely not a pair                           5
# incorrect                                               1
  
  # # 21 minread=30 no downsampling posthoccount
  # indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_maxread_posthoccount/'
  # infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  # infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  # sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021')
  # 
  # library(Phyloscanner.R.utilities)
  # 

  # load(infile.sero)
  # dcln <- data.table(read.csv('~/Dataset_S3.csv'))
  # dcln <- dcln[EPID_EVIDENCE_TYPE=='serodisc',] 
  # aid <- data.table(read.csv(infile.aid19))
  # aid[,ID:=paste0('RK-',ID)]   
  # dcln <- merge(dcln,subset(aid,select=c('ID','AID')),by.x=c('MALE_ID'),by.y='AID',all.x=T)
  # dcln <- merge(dcln,subset(aid,select=c('ID','AID')),by.x=c('FEMALE_ID'),by.y='AID',all.x=T)
  # dcp <- merge(dsero,dcln,by.x=c('FEMALE_RID','MALE_RID'),by.y=c('ID.y','ID.x'),all=T)
  # dcp[is.na(EPID_EVIDENCE_DIR)]
  # dcp[is.na(EXT_DIR)]
  # dcp[EXT_DIR!=EPID_EVIDENCE_DIR,]
  
}

downsampling_sero <- function(){
  # 
  sero_analysis <- function(infile.net,infile.aid,infile.sero,analysis='2021',full=T){
    # couple
    # load(infile.couple)
    # couple <- data.table(unique(subset(coupdat,select=c('male.RCCS_studyid', 'female.RCCS_studyid'))))
    # colnames(couple) <- c('MALE_RID','FEMALE_RID')
    # couple[,COUPLE:=1]
    # couple[,MALE_RID:=paste0('RK-',MALE_RID)]
    # couple[,FEMALE_RID:=paste0('RK-',FEMALE_RID)]
    
    load(infile.sero)
    
    # 
    aid <- data.table(read.csv(infile.aid))
    # consistent names
    if(analysis=='2021'){
      setnames(aid,'PT_ID','ID')
    }
    if(analysis=='2019'){
      aid[,ID:=paste0('RK-',ID)]   
      full <- F
    }
    
    
    
    load(infile.net)
    if(full==T){
      dc <- data.table(dca)
      setnames(dc,c('host.1','host.2','categorisation','n.eff','k.eff','type'),
               c('H1','H2','CATEGORISATION','N_EFF','K_EFF','TYPE'))
      dc[,H1 := gsub('CNTRL-','',H1)]
      dc[,H2 := gsub('CNTRL-','',H2)]
    }else{
      dc <- data.table(dc)  
    }
    
    # unique(dc$CATEGORISATION)
    dclkl <- dc[CATEGORISATION=='close.and.adjacent.and.ancestry.cat']
    tmp <- dclkl[,list(MAX_N_EFF=max(N_EFF)),by=c('H1','H2')]
    dclkl <- merge(dclkl,tmp, by=c('H1','H2'))
    dclkl <- dclkl[MAX_N_EFF==N_EFF,]
    tmp <- dclkl[,list(MIN_PTY_RUN=min(PTY_RUN)),by=c('H1','H2')]
    dclkl <- merge(dclkl,tmp, by=c('H1','H2'))
    dclkl <- dclkl[MIN_PTY_RUN==PTY_RUN,]
    
    # tmp <- dclkl[,length(N_EFF),by=c('H1','H2')]
    # setkey(tmp,V1)
    # dclkl <- dclkl[,head(.SD, 4),by=c('H1','H2')]
    dclkl <- data.table(dcast(dclkl, H1+H2~TYPE,value.var='K_EFF'))
    # revisit
    # colnames(dclkl)
    dclkl[, SCORE_LINK:=(complex.or.no.ancestry+`12`+`21`)/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
    # dclkl[, SCORE_DIR:=(`12`+`21`)/(complex.or.no.ancestry+`12`+`21`)]
    # dclkl[, SCORE_12:=`12`/(`12`+`21`)]
    dclkl[,SCORE_12:=`12`/(`12`+`21`)]
    dclkl[,SCORE_21:=`21`/(`12`+`21`)]
    dclkl[, LINK:=(SCORE_LINK>0.6)]
    # dclkl[LINK==T, DIR:=(SCORE_DIR>0.6)]
    # dclkl[DIR==T, DIR12:=(SCORE_12>0.5)]
    dclkl[LINK==T, DIR12:=(SCORE_12>0.6)]
    dclkl[LINK==T, DIR21:=(SCORE_21>0.6)]
    
    # sort as M F
    if(analysis=='2019'){
      dclkl[,H1S:=substr(H1,9,9)]
      dclkl[,H2S:=substr(H2,9,9)]
      tmp <- dclkl[!(H1S=='M'),]
      setnames(tmp,c('H1','H2','12','21','H1S','H2S'),
               c('H2','H1','21','12','H2S','H1S'))
      tmp[,SCORE_12:=1-SCORE_12]
      tmp[,DIR12:=1-DIR12]
      dclkl <- rbind(tmp,dclkl[(H1S=='M'),])      
      dclkl <- merge(dclkl, subset(aid,select=c('ID','AID')),
                     by.x='H1',by.y='AID',all.x=T)
      dclkl <- merge(dclkl, subset(aid,select=c('ID','AID')),
                     by.x='H2',by.y='AID',all.x=T)
      dclkl <- merge(dclkl, dsero, by.x=c('ID.x','ID.y'),by.y=c('MALE_RID', 'FEMALE_RID'),all.x=T)
    }else{
      dclkl <- merge(dclkl, subset(aid,select=c('ID','AID')),
                     by.x='H1',by.y='AID',all.x=T)
      dclkl <- merge(dclkl, subset(aid,select=c('ID','AID')),
                     by.x='H2',by.y='AID',all.x=T)
      tmp <- merge(dclkl,dsero,by.x=c('ID.x','ID.y'),by.y=c('FEMALE_RID','MALE_RID'))
      dclkl <- merge(dclkl,dsero,by.x=c('ID.x','ID.y'),by.y=c('MALE_RID','FEMALE_RID'))
      setnames(tmp,c('H1','H2','12','21','ID.x','ID.y'),
               c('H2','H1','21','12','ID.y','ID.x'))
      tmp[,SCORE_12:=1-SCORE_12]
      tmp[,DIR12:=1-DIR12]
      dclkl <- rbind(dclkl, tmp)
    }
    
    dclkl <- dclkl[!is.na(EXT_TYPE)]
    dclkl[DIR12==1,PHSC_DIR:='mf']
    dclkl[DIR12==0,PHSC_DIR:='fm']
    
    tmp	<- dclkl[, which(!is.na(PHSC_DIR))]
    set(dclkl, tmp, 'EXT_EVAL', dclkl[tmp, as.character(factor(PHSC_DIR==EXT_DIR, levels=c(TRUE,FALSE),labels=c('correct','incorrect')))])
    dclkl[DIR12==F & DIR21==F,EXT_EVAL:='couple most likely a pair direction not resolved']
    dclkl[LINK==F,EXT_EVAL:='couple most likely not a pair']
    # dclkl[is.na(EXT_EVAL)]
    dclkl[, table(EXT_EVAL, EXT_TYPE, useNA='if')]
    return(dclkl)
  }

  sero_analysis_matthew_threshold <- function(infile.net,infile.aid,infile.sero,analysis='2021',full=T){
    # couple
    # load(infile.couple)
    # couple <- data.table(unique(subset(coupdat,select=c('male.RCCS_studyid', 'female.RCCS_studyid'))))
    # colnames(couple) <- c('MALE_RID','FEMALE_RID')
    # couple[,COUPLE:=1]
    # couple[,MALE_RID:=paste0('RK-',MALE_RID)]
    # couple[,FEMALE_RID:=paste0('RK-',FEMALE_RID)]
    
    load(infile.sero)
    aid <- data.table(read.csv(infile.aid))
    # consistent names
    if(analysis=='2021'){
      setnames(aid,'PT_ID','ID')
    }
    if(analysis=='2019'){
      aid[,ID:=paste0('RK-',ID)]   
      full <- F
    }
    load(infile.net)
    if(full==T){
      dc <- data.table(dca)
      setnames(dc,c('host.1','host.2','categorisation','n.eff','k.eff','type'),
               c('H1','H2','CATEGORISATION','N_EFF','K_EFF','TYPE'))
      dc[,H1 := gsub('CNTRL-','',H1)]
      dc[,H2 := gsub('CNTRL-','',H2)]
    }else{
      dc <- data.table(dc)  
    }
    dclkl <- dc[CATEGORISATION=='close.and.adjacent.and.ancestry.cat']
    tmp <- dclkl[,list(MAX_N_EFF=max(N_EFF)),by=c('H1','H2')]
    dclkl <- merge(dclkl,tmp, by=c('H1','H2'))
    dclkl <- dclkl[MAX_N_EFF==N_EFF,]
    tmp <- dclkl[,list(MIN_PTY_RUN=min(PTY_RUN)),by=c('H1','H2')]
    dclkl <- merge(dclkl,tmp, by=c('H1','H2'))
    dclkl <- dclkl[MIN_PTY_RUN==PTY_RUN,]
    dclkl <- data.table(dcast(dclkl, H1+H2~TYPE,value.var='K_EFF'))
    dclkl[, SCORE_LINK:=(complex.or.no.ancestry+`12`+`21`)/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
    dclkl[,SCORE_12:=`12`/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
    dclkl[,SCORE_21:=`21`/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
    dclkl[, LINK:=(SCORE_LINK>0.5)]
    dclkl[LINK==T, DIR12:=(SCORE_12>0.33)]
    dclkl[LINK==T, DIR21:=(SCORE_21>0.33)]
    
    # sort as M F
    if(analysis=='2019'){
      dclkl[,H1S:=substr(H1,9,9)]
      dclkl[,H2S:=substr(H2,9,9)]
      tmp <- dclkl[!(H1S=='M'),]
      setnames(tmp,c('H1','H2','12','21','H1S','H2S'),
               c('H2','H1','21','12','H2S','H1S'))
      tmp[,SCORE_12:=1-SCORE_12]
      tmp[,DIR12:=1-DIR12]
      dclkl <- rbind(tmp,dclkl[(H1S=='M'),])      
      dclkl <- merge(dclkl, subset(aid,select=c('ID','AID')),
                     by.x='H1',by.y='AID',all.x=T)
      dclkl <- merge(dclkl, subset(aid,select=c('ID','AID')),
                     by.x='H2',by.y='AID',all.x=T)
      dclkl <- merge(dclkl, dsero, by.x=c('ID.x','ID.y'),by.y=c('MALE_RID', 'FEMALE_RID'),all.x=T)
    }else{
      dclkl <- merge(dclkl, subset(aid,select=c('ID','AID')),
                     by.x='H1',by.y='AID',all.x=T)
      dclkl <- merge(dclkl, subset(aid,select=c('ID','AID')),
                     by.x='H2',by.y='AID',all.x=T)
      tmp <- merge(dclkl,dsero,by.x=c('ID.x','ID.y'),by.y=c('FEMALE_RID','MALE_RID'))
      dclkl <- merge(dclkl,dsero,by.x=c('ID.x','ID.y'),by.y=c('MALE_RID','FEMALE_RID'))
      setnames(tmp,c('H1','H2','12','21','ID.x','ID.y'),
               c('H2','H1','21','12','ID.y','ID.x'))
      tmp[,SCORE_12:=1-SCORE_12]
      tmp[,DIR12:=1-DIR12]
      dclkl <- rbind(dclkl, tmp)
    }
    
    dclkl <- dclkl[!is.na(EXT_TYPE)]
    dclkl[DIR12==1,PHSC_DIR:='mf']
    dclkl[DIR12==0,PHSC_DIR:='fm']
    
    tmp	<- dclkl[, which(!is.na(PHSC_DIR))]
    set(dclkl, tmp, 'EXT_EVAL', dclkl[tmp, as.character(factor(PHSC_DIR==EXT_DIR, levels=c(TRUE,FALSE),labels=c('correct','incorrect')))])
    dclkl[DIR12==F & DIR21==F,EXT_EVAL:='couple most likely a pair direction not resolved']
    dclkl[LINK==F,EXT_EVAL:='couple most likely not a pair']
    # dclkl[is.na(EXT_EVAL)]
    dclkl[, table(EXT_EVAL, EXT_TYPE, useNA='if')]
    return(dclkl)
  }
  
  library(data.table)
  infile.sero <- '~/serodata19.rda'
  infile.aid19 <- '~/todi_pairs_171122_cl25_d50_prior23_min30_anonymised_RIDs.csv'
  infile.net19 <- '~/Rakai_phscnetworks_190706.rda'
  # load(infile.net19)
  # range(dc$N_EFF)
  
  # dc_19 <- sero_analysis(infile.net19,infile.aid19,infile.sero,analysis='2019',full=F)
  # correct                                                28
  # couple most likely a pair direction not resolved        5
  # couple most likely not a pair                          13
  # incorrect                                               4
  dc_19 <- sero_analysis_matthew_threshold(infile.net19,infile.aid19,infile.sero,analysis='2019',full=F)
  # 
  dc_ds <- list()
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  
  for (k in 1:10) {
    # HOME <<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
    # indir	<- file.path(HOME,paste0('210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_seed',k,'/'))
    # infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
    # job_tag <- paste0('_02_05_30_min_read_100_max_read_posthoccount_seed',k)
    # job_tag <- paste0('_02_05_30_min_read_100_max_read_posthoccount_im_mrca_seed',k)
    job_tag <- paste0('_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd_seed',k)
    infile.networks <-  paste0('~/dcdwina',job_tag,'.rda')
    # dc_ds[[k]] <- sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021',full=T)
    dc_ds[[k]] <- sero_analysis_matthew_threshold(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021',full=T)
    
  }
  job_tag <- paste0('_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd')
  # job_tag <- paste0('_02_05_30_min_read_null_max_read_posthoccount_im_mrca')
  # job_tag <- paste0('_02_05_30_min_read_null_maxread_posthoccount')
  infile.networks <-  paste0('~/dcdwina',job_tag,'.rda')
  # dc_nods <- sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021',full=T)
  dc_nods <- sero_analysis_matthew_threshold(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021',full=T)
  
  tab <- lapply(dc_ds,function(x)as.data.table(table(x$EXT_EVAL)))
  tab <- rbindlist(tab)
  tab[,paste0(round(median(N),2), ' ( ', paste(round(quantile(N,probs=c(0.025,0.975)),2),collapse = ' - '), ' ) '),by=c('V1')]
  # 1:                                          correct 22 ( 21 - 22.78 ) 
  # 2: couple most likely a pair direction not resolved    10 ( 10 - 11 ) 
  # 3:                    couple most likely not a pair    10 ( 10 - 10 ) 
  # 4:                                        incorrect       3 ( 2 - 4 ) 
  table(dc_nods$EXT_EVAL)
  # correct 
  # 21 
  # couple most likely a pair direction not resolved 
  # 13 
  # couple most likely not a pair 
  # 9 
  # incorrect 
  # 2     
  
  table(dc_19$EXT_EVAL) 
  # check why downsampling exactly same???
  job_tag <- paste0('_02_05_30_min_read_100_max_read_posthoccount_seed',1)
  infile.networks <-  paste0('~/dcdwina',job_tag,'.rda')
  load(infile.networks)
  dc1<- data.table(dca)
  job_tag <- paste0('_02_05_30_min_read_100_max_read_posthoccount_seed',2)
  infile.networks <-  paste0('~/dcdwina',job_tag,'.rda')
  load(infile.networks)
  dc2<- data.table(dca)
  dc1 <- dc1[categorisation=='close.and.adjacent.and.ancestry.cat']
  dc2 <- dc2[categorisation=='close.and.adjacent.and.ancestry.cat']
  dc <- merge(dc1,dc2,by=c(  'PTY_RUN',  'host.1',  'host.2','type'))
  dc[k.eff.x!=k.eff.y]
  
  # sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021',full=T)
  # 
  # infile.net= infile.networks
  # infile.aid=infile.ind.anonymised
  # analysis='2021'
  # full=T
  
  job_tag <- paste0('_02_05_30_min_read_100_max_read_posthoccount_im')
  infile.networks <-  paste0('~/dcdwina',job_tag,'.rda')
  dc_im <- sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021',full=T)
  job_tag <- paste0('_02_05_30_min_read_100_max_read_posthoccount')
  infile.networks <-  paste0('~/dcdwina',job_tag,'.rda')
  dc_21 <- sero_analysis(infile.networks,infile.ind.anonymised,infile.sero,analysis='2021',full=T)
  
  dc_21[, table(EXT_EVAL, EXT_TYPE, useNA='if')]
  dc_im[, table(EXT_EVAL, EXT_TYPE, useNA='if')]
  
  tmp1 <- subset(dc_im,select=c('ID.x','ID.y','EXT_EVAL'))
  tmp2 <- subset(dc_21,select=c('ID.x','ID.y','EXT_EVAL'))
  tmp <- tmp1[!(ID.x<ID.y)]
  setnames(tmp,c('ID.x','ID.y'),c('ID.y','ID.x'))
  tmp1 <- rbind(tmp1[(ID.x<ID.y)],tmp)
  tmp <- tmp2[!(ID.x<ID.y)]
  setnames(tmp,c('ID.x','ID.y'),c('ID.y','ID.x'))
  tmp2 <- rbind(tmp2[(ID.x<ID.y)],tmp)
  dc_im21 <- merge(tmp1,tmp2,by=c('ID.x','ID.y'))
  dc_im21[, table(EXT_EVAL.x)]
  dc_im21[, table(EXT_EVAL.y)]
  
  library(data.table)
  job_tag <- paste0('_02_05_30_min_read_100_max_read_posthoccount')
  infile.networks <-  paste0('~/dcdwina',job_tag,'.rda')
  load( infile.networks)
  dca1 <- data.table(dca)
  dwina1 <- data.table(dwina)
  
  job_tag <- paste0('_02_05_30_min_read_100_max_read_posthoccount_im')
  infile.networks <-  paste0('~/dcdwina',job_tag,'.rda')
  load( infile.networks)
  dca2 <- data.table(dca)
  dwina2 <- data.table(dwina)
  
  tmp1 <- dca1[PTY_RUN==249 & categorisation=='close.and.adjacent.and.ancestry.cat',]
  tmp2 <- dca2[PTY_RUN==249 & categorisation=='close.and.adjacent.and.ancestry.cat',]
  tmp <- merge(tmp1,tmp2, by=c('host.1','host.2','categorisation','PTY_RUN','type'),all=T)
  # tmp[,c('host.1','host.2','k.eff.x','n.eff.x','k.eff.y','n.eff.y')]
  tmp[k.eff.x!=k.eff.y]
  tmp1[host.2=='AID3114' & host.1=='AID2722']
  tmp2[host.2=='AID3114' & host.1=='AID2722']
  
  # tmp1 <- dwina1[PTY_RUN==249 ,]
  # tmp2 <- dwina2[PTY_RUN==249 ,]
  tmp1 <- dwina1[PTY_RUN==349 ,]
  tmp2 <- dwina2[PTY_RUN==349 ,]
  tmp <- merge(tmp1, tmp2, by=c('PTY_RUN','host.1','host.2','tree.id'),all=T)
  
  tmp[close.and.adjacent.and.ancestry.cat.x!= close.and.adjacent.and.ancestry.cat.y,
      c('PTY_RUN','host.1','host.2','tree.id','close.and.adjacent.and.ancestry.cat.x','close.and.adjacent.and.ancestry.cat.y')]
  tmp[grepl('multi',ancestry.x),table(close.and.adjacent.and.ancestry.cat.x,close.and.adjacent.and.ancestry.cat.y)]
  tmp[grepl('multi',ancestry.x) & close.and.adjacent.and.ancestry.cat.x!='not.close.or.nonadjacent']
  tmp[grepl('multi',ancestry.x) & close.and.adjacent.and.ancestry.cat.x!='not.close.or.nonadjacent',table(host.1,host.2)] 
  tmp[host.1=='AID1371' & host.2=='AID6853' & grepl('multi',ancestry.x) & close.and.adjacent.and.ancestry.cat.x!='not.close.or.nonadjacent' &
        close.and.adjacent.and.ancestry.cat.x==close.and.adjacent.and.ancestry.cat.y,'tree.id']   
  tmp[host.1=='AID1371' & host.2=='AID6853' & grepl('multi',ancestry.x) & close.and.adjacent.and.ancestry.cat.x!='not.close.or.nonadjacent' &
        close.and.adjacent.and.ancestry.cat.x!=close.and.adjacent.and.ancestry.cat.y,'tree.id']   
  # tmp[host.1=='AID0719' & host.2=='AID1311' & close.and.adjacent.and.ancestry.cat.x!= close.and.adjacent.and.ancestry.cat.y,
  #     c('PTY_RUN','host.1','host.2','tree.id','close.and.adjacent.and.ancestry.cat.x','close.and.adjacent.and.ancestry.cat.y')]
  # tmp[host.1=='AID0719' & host.2=='AID1311' & close.and.adjacent.and.ancestry.cat.x== close.and.adjacent.and.ancestry.cat.y,
  #     c('PTY_RUN','host.1','host.2','tree.id','close.and.adjacent.and.ancestry.cat.x','close.and.adjacent.and.ancestry.cat.y')]
  # tmp[host.1=='AID0719' & host.2=='AID1311',table(ancestry.x,ancestry.y)]
  # tmp[host.1=='AID0719' & host.2=='AID1311' & grepl('multi',ancestry.x)]
  # tmp[host.1=='AID0719' & host.2=='AID1311' & grepl('multi',ancestry.x),table(close.and.adjacent.and.ancestry.cat.x,close.and.adjacent.and.ancestry.cat.y)]
  load('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI//210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im/ptyr249_workspace.rda')
  dwim <- data.table(dwin)
  load('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI//210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount/ptyr249_workspace.rda')
  dw <- data.table(dwin)
  tmp <- merge(dwim, dw, by=c('host.1','host.2','tree.id'),all=T)
  tmp[close.and.adjacent.and.ancestry.cat.x!= close.and.adjacent.and.ancestry.cat.y]
  tmp[ancestry.x=='multiAnc',c(close.and.adjacent.and.ancestry.cat.y)]
  x=   load('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI//210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im/ptyr249_workspace.rda')
}


cost_of_zerolengthadj <- function(){
  library(data.table)
  setwd('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount')
  # setwd('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca')
  setwd('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd')
  setwd('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount')
  # setwd('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca')
  setwd('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd')

  
  load('Rakai_phscnetworks.rda')
  # unique(dc$CATEGORISATION)
  dclkl <- dc[CATEGORISATION=='close.and.adjacent.and.ancestry.cat']
  tmp <- dclkl[,list(MAX_N_EFF=max(N_EFF)),by=c('H1','H2')]
  dclkl <- merge(dclkl,tmp, by=c('H1','H2'))
  dclkl <- dclkl[MAX_N_EFF==N_EFF,]
  tmp <- dclkl[,list(MIN_PTY_RUN=min(PTY_RUN)),by=c('H1','H2')]
  dclkl <- merge(dclkl,tmp, by=c('H1','H2'))
  dclkl <- dclkl[MIN_PTY_RUN==PTY_RUN,]
  
  # tmp <- dclkl[,length(N_EFF),by=c('H1','H2')]
  # setkey(tmp,V1)
  # dclkl <- dclkl[,head(.SD, 4),by=c('H1','H2')]
  dclkl <- data.table(dcast(dclkl, H1+H2~TYPE,value.var='K_EFF'))
  # revisit
  # colnames(dclkl)
  dclkl[, SCORE_LINK:=(complex.or.no.ancestry+`12`+`21`)/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
  # dclkl[, SCORE_DIR:=(`12`+`21`)/(complex.or.no.ancestry+`12`+`21`)]
  # dclkl[, SCORE_12:=`12`/(`12`+`21`)]
  if(0){
    dclkl[,SCORE_12:=`12`/(`12`+`21`+not.close.or.nonadjacent+complex.or.no.ancestry)]
    dclkl[,SCORE_21:=`21`/(`12`+`21`+not.close.or.nonadjacent+complex.or.no.ancestry)]
    dclkl[, LINK:=(SCORE_LINK>0.5)]
    dclkl[LINK==T, DIR12:=(SCORE_12>0.33)]
    dclkl[LINK==T, DIR21:=(SCORE_21>0.33)]
  }
  if(1){
    dclkl[,SCORE_12:=`12`/(`12`+`21`)]
    dclkl[,SCORE_21:=`21`/(`12`+`21`)]
    dclkl[, LINK:=(SCORE_LINK>0.6)]
    dclkl[LINK==T, DIR12:=(SCORE_12>0.6)]
    dclkl[LINK==T, DIR21:=(SCORE_21>0.6)]
  }

  
  # dclkl[LINK==T, DIR:=(SCORE_DIR>0.6)]
  # dclkl[DIR==T, DIR12:=(SCORE_12>0.5)]
  # nrow(dclkl[LINK==F])
  # nrow(dclkl[LINK==T])
  # dclkl[LINK==T, table(DIR)]
  # hist(dclkl[LINK==T, SCORE_DIR])
  # nrow(dclkl[DIR==T])
  
  nrow(dclkl[LINK==T])
  nrow(dclkl[DIR12==T|DIR21==T])
  
  # [1] 930
  # [1] 466
  # [1] 248
  # > 
  
  # [1] 928
  # [1] 639
  # [1] 333
  nrow(dclkl[LINK==T])
  nrow(dclkl[DIR12==T|DIR21==T])
  
  # dclkl[LINK==T, SCORE_12:=`12`/(`12`+`21`+complex.or.no.ancestry)]
  # dclkl[LINK==T, SCORE_21:=`21`/(`12`+`21`+complex.or.no.ancestry)]
  # dclkl[LINK==T, DIR12_M := (SCORE_12>0.33)]
  # dclkl[LINK==T, DIR21_M := (SCORE_21>0.33)]
  # nrow(dclkl[DIR12_M==T])
  # nrow(dclkl[DIR21_M==T])
  # [1] 356
  # [1] 360
  
  # [1] 430
  # [1] 429
  # nrow(dclkl[LINK==T& DIR12_M==F & DIR21_M==F])
}



downsample_ratio <- function(){
  library(purrr)
  library(dplyr)
  library(phyloscannerR)
  library(data.table)

  indir.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships'
  
  # imbalanced pairs no downsampling
  imbalance_read_df <- data.table()
  job_tag <- '_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd'
  indir	<- paste0(indir.base,job_tag)
  for (k in 1:409) {
    try({
      cat('processing ', k, ' th run \n')
      load(file.path(indir,paste0('ptyr',k,'_workspace.rda')))
      host.tips.and.reads <- map(phyloscanner.trees, function(x) phyloscannerR:::get.tip.and.read.counts(x, all.hosts.from.trees(phyloscanner.trees), tip.regex, attr(phyloscanner.trees, 'has.read.counts'), verbose = F))
      host.tips.and.reads <- bind_rows(host.tips.and.reads)
      dwin <- dwin %>% 
        inner_join(host.tips.and.reads, by=c("host.1"="host.id", "tree.id")) %>% 
        rename(tips.1 = tips, reads.1=reads)
      dwin <- dwin %>% 
        inner_join(host.tips.and.reads, by=c("host.2"="host.id", "tree.id")) %>% 
        rename(tips.2 = tips, reads.2=reads)
      dwin <- as.data.table(dwin)
      dwin[,PTY_RUN:=k]
      imbalance_read_df <- rbind(imbalance_read_df,dwin[abs(reads.1-reads.2)>100]) 
    })
  }
  saveRDS(imbalance_read_df,file = file.path(indir,'imbalance_read_df.rds'))
  
  # imbalanced pairs downsampling
  indir.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships'
  
  imbalance_read_df <- data.table()
  job_tag <- '_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd'
  indir	<- paste0(indir.base,job_tag)
  iddf <- readRDS(file.path( paste0(indir.base,'_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd'),'imbalance_read_df.rds'))
  iddf <- subset(iddf,select=c('host.1','host.2','tree.id','tips.1', 'reads.1', 'tips.2', 'reads.2','PTY_RUN'))
  for (k in 1:409) {
    try({  
      cat('processing ', k, ' th run \n')
      load(file.path(indir,paste0('ptyr',k,'_workspace.rda')))
      tmp <- iddf[PTY_RUN==k]
      dwin <- merge(dwin,tmp,by=c('host.1',  'host.2','tree.id'),all.y=T)
      imbalance_read_df <- rbind(imbalance_read_df,dwin)
    })
  }
  saveRDS(imbalance_read_df,file = file.path(indir,'imbalance_read_df.rds'))
  
  
  indir.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships'
  job_tag <- '_02_05_30_min_read_null_max_read_posthoccount_im_mrca'
  indir	<- paste0(indir.base,job_tag)
  imbalance_read_df <- readRDS(file.path(indir,'imbalance_read_df.rds'))
  # comment to consider min=1
  # imbalance_read_df <- imbalance_read_df[reads.1>=30 & reads.2>=30]
  tmp <- imbalance_read_df[(close.and.adjacent.and.ancestry.cat=='12' & reads.1>reads.2)|
                             (close.and.adjacent.and.ancestry.cat=='21' & reads.1<reads.2)]
  nrow(tmp)
  nrow(imbalance_read_df)
  nrow(tmp)/nrow(imbalance_read_df)
  
  job_tag <- '_02_05_30_min_read_100_max_read_posthoccount_im_mrca'
  indir	<- paste0(indir.base,job_tag)
  imbalance_read_df_downsampling <- readRDS(file.path(indir,'imbalance_read_df.rds'))
  # imbalance_read_df_downsampling <- imbalance_read_df_downsampling[reads.1>=30 & reads.2>=30]
  tmp_downsampling <- imbalance_read_df_downsampling[(close.and.adjacent.and.ancestry.cat=='12' & reads.1>reads.2)|
                                                       (close.and.adjacent.and.ancestry.cat=='21' & reads.1<reads.2)]
  nrow(tmp_downsampling)
  nrow(imbalance_read_df_downsampling)
  nrow(tmp_downsampling)/nrow(imbalance_read_df_downsampling)
  # whether to use posthoccount
  nrow(tmp_downsampling)/nrow(tmp)
  
  # pairs with clinical data
  # tmp <- dclnid[EPID_EVIDENCE_TYPE=='serodisc']
  infile.sero <- '~/serodata19.rda'
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  
  load(infile.sero)
  aid <- data.table(read.csv(infile.ind.anonymised))
  aid[,AID:=as.character(AID)]
  tmp <- merge(dsero, aid, by.x='MALE_RID',by.y='PT_ID')
  tmp <- merge(tmp, aid, by.x='FEMALE_RID',by.y='PT_ID')
  setnames(tmp,c('AID.x','AID.y'),c('AID1','AID2'))
  tmp <- subset(tmp,select=c('AID1','AID2','EXT_DIR'))
  tmp[,EXT_DIR := ifelse(EXT_DIR=='mf','12','21')]
  tmp2 <- tmp[!(AID1<AID2)]
  setnames(tmp2,c('AID1','AID2'),c('AID2','AID1'))
  set(tmp2, NULL, 'EXT_DIR', tmp2[, gsub('xx','12',gsub('12','21',gsub('21','xx',EXT_DIR)))])
  tmp <- rbind(tmp[(AID1<AID2)],tmp2)
  
  imbalance_read_df2 <- merge(imbalance_read_df,tmp,by.x=c('host.1','host.2'),by.y=c('AID1','AID2'))
  imbalance_read_df_downsampling2 <- merge(imbalance_read_df_downsampling,tmp,by.x=c('host.1','host.2'),by.y=c('AID1','AID2'))
  
  tmp2 <- imbalance_read_df2[(close.and.adjacent.and.ancestry.cat=='12' & reads.1>reads.2)|
                               (close.and.adjacent.and.ancestry.cat=='21' & reads.1<reads.2)]
  nrow(tmp2)
  nrow(imbalance_read_df2)
  nrow(tmp2)/nrow(imbalance_read_df2)
  
  tmp_downsampling2 <- imbalance_read_df_downsampling2[(close.and.adjacent.and.ancestry.cat=='12' & reads.1>reads.2)|
                                                         (close.and.adjacent.and.ancestry.cat=='21' & reads.1<reads.2)]
  
  nrow(tmp_downsampling2)
  nrow(imbalance_read_df_downsampling2)
  nrow(tmp_downsampling2)/nrow(imbalance_read_df_downsampling2)
  
  nrow(tmp_downsampling2)/nrow(tmp2)
  
  # couple
  infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
  load(infile.couple)
  couple <- data.table(subset(coupdat,select=c('male.RCCS_studyid','female.RCCS_studyid')))
  couple[,male.RCCS_studyid:=paste0('RK-',male.RCCS_studyid)]
  couple[,female.RCCS_studyid:=paste0('RK-',female.RCCS_studyid)]
  aid <- data.table(read.csv(infile.ind.anonymised))
  aid[,AID:=as.character(AID)]
  tmp <- merge(couple, aid, by.x='male.RCCS_studyid',by.y='PT_ID')
  tmp <- merge(tmp, aid, by.x='female.RCCS_studyid',by.y='PT_ID')
  setnames(tmp,c('AID.x','AID.y'),c('AID1','AID2'))
  tmp <- subset(tmp,select=c('AID1','AID2'))
  tmp2 <- tmp[!(AID1<AID2)]
  setnames(tmp2,c('AID1','AID2'),c('AID2','AID1'))
  tmp <- rbind(tmp[(AID1<AID2)],tmp2)
  
  imbalance_read_df2 <- merge(imbalance_read_df,tmp,by.x=c('host.1','host.2'),by.y=c('AID1','AID2'))
  imbalance_read_df_downsampling2 <- merge(imbalance_read_df_downsampling,tmp,by.x=c('host.1','host.2'),by.y=c('AID1','AID2'))
  
  tmp2 <- imbalance_read_df2[(close.and.adjacent.and.ancestry.cat=='12' & reads.1>reads.2)|
                               (close.and.adjacent.and.ancestry.cat=='21' & reads.1<reads.2)]
  nrow(tmp2)
  nrow(imbalance_read_df2)
  nrow(tmp2)/nrow(imbalance_read_df2)
  
  tmp_downsampling2 <- imbalance_read_df_downsampling2[(close.and.adjacent.and.ancestry.cat=='12' & reads.1>reads.2)|
                                                         (close.and.adjacent.and.ancestry.cat=='21' & reads.1<reads.2)]
  
  nrow(tmp_downsampling2)
  nrow(imbalance_read_df_downsampling2)
  nrow(tmp_downsampling2)/nrow(imbalance_read_df_downsampling2)
  
  nrow(tmp_downsampling2)/nrow(tmp2)
  
}

compare_phyloscan_multitran <- function(){
  library(data.table)
  infile.sero <- '~/serodata19.rda'
  # infile.aid19 <- '~/todi_pairs_171122_cl25_d50_prior23_min30_anonymised_RIDs.csv'
  # infile.net19 <- '~/Rakai_phscnetworks_190706.rda'
  indir1	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount/'
  # indir1	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount/'
  
  infile.networks1 <-  file.path(indir1,'Rakai_phscnetworks.rda')
  # indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca/'
  # indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca/'
  indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/'
  # indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd/'
  
  infile.networks2 <-  file.path(indir2,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  
  load(infile.sero)
  aid <- data.table(read.csv(infile.ind.anonymised))
  load(infile.networks1)
  dc1 <- copy(dc)
  dw1 <- copy(dw)
  load(infile.networks2)
  dw2 <- copy(dw)
  dc2 <- copy(dc)
  
  aid[,X:=NULL]
  dsero <- merge(dsero,aid, by.x='FEMALE_RID',by.y='PT_ID',all.x=T)
  dsero <- merge(dsero,aid, by.x='MALE_RID',by.y='PT_ID',all.x=T)
  dsero[,AID.x:=as.character(AID.x)]
  dsero[,AID.y:=as.character(AID.y)]
  dsero[,SERO_DIR:=ifelse(EXT_DIR=='fm','12','21')]
  dsero <- subset(dsero,select=c('AID.x','AID.y','SERO_DIR'))
  
  tmp  <- dsero[!(AID.x<AID.y)]
  setnames(tmp,c('AID.x','AID.y'),c('AID.y','AID.x'))
  tmp[,SERO_DIR:=ifelse(SERO_DIR=='12','21','12')]
  dsero <- rbind(dsero[(AID.x<AID.y)],tmp)
  

  dclkl1 <- dc1[CATEGORISATION=='close.and.adjacent.and.ancestry.cat']
  tmp <- dclkl1[,list(MAX_N_EFF=max(N_EFF)),by=c('H1','H2')]
  dclkl1 <- merge(dclkl1,tmp, by=c('H1','H2'))
  dclkl1 <- dclkl1[MAX_N_EFF==N_EFF,]
  tmp <- dclkl1[,length(N_EFF),by=c('H1','H2')]
  setkey(tmp,V1)
  dclkl1 <- dclkl1[,head(.SD, 4),by=c('H1','H2')]
  dclkl1 <- data.table(dcast(dclkl1, H1+H2~TYPE,value.var='K_EFF'))
  dclkl1[, SCORE_LINK:=(complex.or.no.ancestry+`12`+`21`)/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
  # dclkl1[, SCORE_DIR:=(`12`+`21`)/(complex.or.no.ancestry+`12`+`21`)]
  # dclkl1[, SCORE_12:=`12`/(`12`+`21`)]
  if(0){
    dclkl1[, SCORE_12:=`12`/(`12`+`21`)]
    dclkl1[, SCORE_21:=`21`/(`12`+`21`)]
    dclkl1[, LINK:=(SCORE_LINK>0.6)]
    # dclkl1[LINK==T, DIR:=(SCORE_DIR>0.6)]
    # dclkl1[DIR==T, DIR12:=(SCORE_12>0.5)]
    dclkl1[LINK==T, DIR12:=(SCORE_12>0.6)]
    dclkl1[LINK==T, DIR21:=(SCORE_21>0.6)]
  }
  if(1){
    dclkl1[, SCORE_12:=`12`/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
    dclkl1[, SCORE_21:=`21`/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
    dclkl1[, LINK:=(SCORE_LINK>0.5)]
    dclkl1[LINK==T, DIR12:=(SCORE_12>0.33)]
    dclkl1[LINK==T, DIR21:=(SCORE_21>0.33)]
  }
  
  
  dclkl2 <- dc2[CATEGORISATION=='close.and.adjacent.and.ancestry.cat']
  tmp <- dclkl2[,list(MAX_N_EFF=max(N_EFF)),by=c('H1','H2')]
  dclkl2 <- merge(dclkl2,tmp, by=c('H1','H2'))
  dclkl2 <- dclkl2[MAX_N_EFF==N_EFF,]
  tmp <- dclkl2[,length(N_EFF),by=c('H1','H2')]
  setkey(tmp,V1)
  dclkl2 <- dclkl2[,head(.SD, 4),by=c('H1','H2')]
  dclkl2 <- data.table(dcast(dclkl2, H1+H2~TYPE,value.var='K_EFF'))
  dclkl2[, SCORE_LINK:=(complex.or.no.ancestry+`12`+`21`)/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
  # dclkl2[, SCORE_DIR:=(`12`+`21`)/(complex.or.no.ancestry+`12`+`21`)]
  # dclkl2[, SCORE_12:=`12`/(`12`+`21`)]
  if(0){
    dclkl2[, SCORE_12:=`12`/(`12`+`21`)]
    dclkl2[, SCORE_21:=`21`/(`12`+`21`)]
    dclkl2[, LINK:=(SCORE_LINK>0.6)]
    # dclkl2[LINK==T, DIR:=(SCORE_DIR>0.6)]
    # dclkl2[DIR==T, DIR12:=(SCORE_12>0.5)]
    dclkl2[LINK==T, DIR12:=(SCORE_12>0.6)]
    dclkl2[LINK==T, DIR21:=(SCORE_21>0.6)]
  }
  if(1){
    dclkl2[, SCORE_12:=`12`/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
    dclkl2[, SCORE_21:=`21`/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
    dclkl2[, LINK:=(SCORE_LINK>0.5)]
    dclkl2[LINK==T, DIR12:=(SCORE_12>0.33)]
    dclkl2[LINK==T, DIR21:=(SCORE_21>0.33)]
  }
  dclkl <- merge(dclkl1,dclkl2,by=c('H1','H2'),all=T)
  dclkl <- merge(dclkl, dsero, by.x=c('H1','H2'),by.y=c('AID.x','AID.y'))
  # dclkl[LINK.x!=LINK.y]
  # dclkl[DIR.x!=DIR.y]
  # dclkl[DIR.x!=DIR.y,(DIR12.x==T & SERO_DIR=='12')|(DIR12.x==F & SERO_DIR=='21')]
  # # 2 false 7 true
  # dclkl[DIR12.x!=DIR12.y]
  # # phyloscan for these pairs
  
  # tmp <- dclkl[LINK.x==T & LINK.y==T]
  # tmp[ (DIR12.x==T & SERO_DIR=='12')|(DIR12.x==F & SERO_DIR=='21'),TYPE.x:='consistent']
  # tmp[ !((DIR12.x==T & SERO_DIR=='12')|(DIR12.x==F & SERO_DIR=='21')),TYPE.x:='inconsistent']
  # tmp[DIR.x==F, TYPE.x:='unresolved']
  # tmp[ (DIR12.y==T & SERO_DIR=='12')|(DIR12.y==F & SERO_DIR=='21'),TYPE.y:='consistent']
  # tmp[ !((DIR12.y==T & SERO_DIR=='12')|(DIR12.y==F & SERO_DIR=='21')),TYPE.y:='inconsistent']
  # tmp[DIR.y==F, TYPE.y:='unresolved']
  # tmp[, table(TYPE.x, TYPE.y, useNA='if')]
  # tmp[, table(TYPE.x, useNA='if')]
  # tmp[, table(TYPE.y, useNA='if')]
  tmp <- dclkl[LINK.x==T & LINK.y==T]
  tmp[ (DIR12.x==T & SERO_DIR=='12')|(DIR21.x==T & SERO_DIR=='21'),TYPE.x:='consistent']
  tmp[ (DIR12.x==T & SERO_DIR=='21')|(DIR21.x==T & SERO_DIR=='12'),TYPE.x:='inconsistent']
  tmp[is.na(TYPE.x), TYPE.x:='unresolved']
  tmp[ (DIR12.y==T & SERO_DIR=='12')|(DIR21.y==T & SERO_DIR=='21'),TYPE.y:='consistent']
  tmp[ (DIR12.y==T & SERO_DIR=='21')|(DIR21.y==T & SERO_DIR=='12'),TYPE.y:='inconsistent']
  tmp[is.na(TYPE.y), TYPE.y:='unresolved']
  tmp[, table(TYPE.x, TYPE.y, useNA='if')]
  tmp[, table(TYPE.x, useNA='if')]
  tmp[, table(TYPE.y, useNA='if')]
  tmp <- tmp[TYPE.y=='inconsistent']
  # dclkl1 <- merge(dclkl1, dsero, by.x=c('H1','H2'),by.y=c('AID.x','AID.y'))
  # dclkl2 <- merge(dclkl2, dsero, by.x=c('H1','H2'),by.y=c('AID.x','AID.y'))
  # dclkl1[(DIR12==T & SERO_DIR=='12')|(DIR12==F & SERO_DIR=='21') ]
  # dclkl2[(DIR12==T & SERO_DIR=='12')|(DIR12==F & SERO_DIR=='21') ]
  # 
  # dclkl1[!((DIR12==T & SERO_DIR=='12')|(DIR12==F & SERO_DIR=='21')) ]
  # dclkl2[!((DIR12==T & SERO_DIR=='12')|(DIR12==F & SERO_DIR=='21')) ]
  # tmp <- dclkl[DIR.x!=DIR.y]
   # dir.create('~/compare_multitrans_nodownsampling/')
  
  tmp <- dclkl[DIR12.x!=DIR12.y|DIR21.x!=DIR21.y]
  
  # dw1[H1== "AID1755" & H2=="AID6465" & TREE_ID=="1225_to_1474"]
  # dw2[H1== "AID1755" & H2=="AID6465" & TREE_ID=="1225_to_1474"]
  # load("~/dcdwina_02_05_30_min_read_100_max_read_posthoccount_im_mrca.rda")
  # dwina[host.1== "AID1755" & host.2=="AID6465" & tree.id=="1225_to_1474"]
  
  # dir.create('~/compare_multitrans_downsampling_thresholdA/')
  # dir.create('~/compare_multitrans_nodownsampling_thresholdA/')
  # dir.create('~/compare_multitrans_downsampling_thresholdB/')
  # dir.create('~/compare_multitrans_nodownsampling_thresholdB/')
  dir.create('~/compare_multitrans_fixpd_downsampling_thresholdA/')
  dir.create('~/compare_multitrans_fixpd_nodownsampling_thresholdA/')
  dir.create('~/compare_multitrans_fixpd_downsampling_thresholdB/')
  dir.create('~/compare_multitrans_fixpd_nodownsampling_thresholdB/')
  dir.create('~/compare_multitrans_fixpd_inconsistent_downsampling_thresholdA/')
  dir.create('~/compare_multitrans_fixpd_inconsistent_nodownsampling_thresholdA/')
  control = list(	yintercept_close=0.025,
                  yintercept_dist=1.0,
                  breaks_x=seq(0,1e4,500), 
                  minor_breaks_x=seq(0,1e4,100),
                  breaks_y=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),										
                  limits_y=c(1e-3,0.4),
                  fill.topology=c("12"="deepskyblue1","21"="deepskyblue4","complex.or.no.ancestralestry"= "#FDB863","not.close.or.nonadjacent"="grey80"))
  
  library(phyloscannerR)
  library(tibble)
  library(ggpubr)
  library(ggplot2)
  setnames(dw1, c('H1','H2','TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  setnames(dw2, c('H1','H2','TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  dw2[,PTY_RUN_Y:=NULL]
  setnames(dw2, 'PTY_RUN_X','PTY_RUN',skip_absent = T)
  tmpdw <- dw1[,list(PTY_RUN=min(PTY_RUN)),by=c('host.1', 'host.2')]
  dw1 <- merge(dw1,tmpdw,by=c('host.1', 'host.2','PTY_RUN'))
  tmpdw <- dw2[,list(PTY_RUN=min(PTY_RUN)),by=c('host.1', 'host.2')]
  dw2 <- merge(dw2,tmpdw,by=c('host.1', 'host.2','PTY_RUN'))
  
  setkey(tmp,H1,H2)
  for (i in 1:nrow(tmp)) {
    cat('processing ',i,' out of ', nrow(tmp), ' pairs \n')
    label <- paste0(tmp$DIR12.x[i],'_',tmp$DIR21.x[i],'_',tmp$DIR12.y[i],'_',tmp$DIR21.y[i],'_',tmp$SERO_DIR[i])
    hosts <- c(tmp$H1[i],tmp$H2[i])
    tmpx <- dw1[(host.1==hosts[1] & host.2==hosts[2])|(host.1==hosts[2] & host.2==hosts[1]),'tree.id']
    tmpx <- rbind(tmpx,dw2[(host.1==hosts[1] & host.2==hosts[2])|(host.1==hosts[2] & host.2==hosts[1]),'tree.id'])
    tmpx <- as.integer(gsub('([0-9]+)_to_([0-9]+)','\\1',tmpx$tree.id))
    control$limits_x <- range(tmpx) + c(-25,25) # bin width
    tmpp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw1), inclusion = "both",control=control)
    g1 <- tmpp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    tmpp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw2), inclusion = "both",control=control)
    g2 <- tmpp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    arrange <- ggarrange(g1 + ggtitle('no adjustment for zero branch lengths \n '),
                         g2 + ggtitle('adjustment for zero branch lengths \n '), 
                         ncol = 1, nrow = 2, common.legend = T, legend='bottom')
    # ggsave(paste0('~/compare_multitrans_nodownsampling/',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 8)
    # ggsave(paste0('~/compare_multitrans_downsampling_thresholdA/',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 8)
    # ggsave(paste0('~/compare_multitrans_nodownsampling_thresholdA/',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 8)
    # ggsave(paste0('~/compare_multitrans_downsampling_thresholdB/',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 8)
    # ggsave(paste0('~/compare_multitrans_nodownsampling_thresholdB/',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 8)
    # ggsave(paste0('~/compare_multitrans_fixpd_nodownsampling_thresholdA/',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 8)
    # ggsave(paste0('~/compare_multitrans_fixpd_nodownsampling_thresholdB/',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 8)
    # ggsave(paste0('~/compare_multitrans_fixpd_downsampling_thresholdB/',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 8)
    # ggsave(paste0('~/compare_multitrans_fixpd_downsampling_thresholdA/',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 8)
    # ggsave(paste0('~/compare_multitrans_fixpd_inconsistent_downsampling_thresholdA/',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 8)
    ggsave(paste0('~/compare_multitrans_fixpd_inconsistent_nodownsampling_thresholdA/',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 8)
    
  }
  i=3
  i=4
  i=6
  i=7
  i=9
  hosts <- c(tmp$H1[i],tmp$H2[i])
  tmp1 <- dw1[(host.1==hosts[1] & host.2==hosts[2])|(host.1==hosts[2] & host.2==hosts[1]),]
  tmp2 <- dw2[(host.1==hosts[1] & host.2==hosts[2])|(host.1==hosts[2] & host.2==hosts[1]),]
  tmp3 <- merge(tmp1,tmp2,by='tree.id')
  tmp3[basic.classification.y!=basic.classification.x,c('tree.id','basic.classification.x','basic.classification.y')]
  tmp3[basic.classification.y!=basic.classification.x,table(basic.classification.x)]
  # i=3
  # tmp3[tree.id=='1300_to_1549']
  # setwd(indir1)
  # library(phyloscannerR)
  # library(treeio)
  indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca/'
  indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca/'
  
  indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_maxread_posthoccount/'
  
  setwd(indir2)
  
  load('ptyr249_workspace.rda')
  ptree <- phyloscanner.trees[["1250_to_1499"]]
  
  load('ptyr299_workspace.rda')
  ptree <- phyloscanner.trees[["1300_to_1549"]]
  
  load('ptyr396_workspace.rda')
  ptree <- phyloscanner.trees[["1250_to_1499"]]
  
  load('ptyr209_workspace.rda')
  ptree <- phyloscanner.trees[["1250_to_1499"]]
  
  indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount/'
  
  indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca/'
  indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca/'
  
  indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/'
  indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd/'
  
  setwd(indir2)
  # load('ptyr76_workspace.rda')
  # ptree <- phyloscanner.trees[["1225_to_1474"]]
  load('ptyr251_workspace.rda')
  load('ptyr389_workspace.rda')
  ptree <- phyloscanner.trees[["1100_to_1349"]]
  ptree <- phyloscanner.trees[["1750_to_1999"]]
  a <- ptree$classification.results[[1]]
  a <- data.table(a)
  a[host.1=='AID2168' & host.2=='AID4448']
  pat.1 <- 8
  pat.2 <- 9
  
  pat.1 <- 5
  pat.2 <- 6
  
  
  
  # indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca/'
  indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount/'
  setwd(indir2)
  load('ptyr407_workspace.rda')
  ptree <- phyloscanner.trees[["1225_to_1474"]]
  library(phyloscannerR)
  a <-classify(ptree, allow.mt = T, n.mt=Inf, p.mt= Inf, identify.multifurcation=F, relaxed.ancestry = F,verbose = F, no.progress.bars = F)
  a <-data.table(a[[1]])
  a[host.1=='AID4452' & host.2=='AID6357']
  
  allow.mt = T
  n.mt=Inf
  p.mt= Inf
  # identify.multifurcation=T
  identify.multifurcation=F
  relaxed.ancestry = F
  verbose = F
  no.progress.bars = F

  library(ape)
  library(tibble)
  library(phangorn)
  library(phyloscannerR)
  library(ggtree)
  
  grep('4452',hosts.included)
  grep('6357',hosts.included)
  pat.1 <- 19
  pat.2 <- 40

  library(data.table)
  setwd("/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount")  
  # x=load("ptyr407_workspace.rda")
  # grep('dw',x,value=T)
  load("ptyr407_workspace.rda")
  ptree <- phyloscanner.trees[["1225_to_1474"]]
  a <- ptree$classification.results[[1]]
  a <- data.table(a)
  a[host.1=='AID4452' & host.2=='AID6357']
  classification <-classify(ptree, allow.mt = T, n.mt=Inf, p.mt= Inf, identify.multifurcation=F, relaxed.ancestry = F,verbose = F, no.progress.bars = F)
  a <- classification[[1]]
  a <- data.table(a)
  a[host.1=='AID4452' & host.2=='AID6357']
  dwin <- data.table(dwin)
  dwin[host.1=='AID4452' & host.2=='AID6357' & grepl('1225_',tree.id)]
  
  setwd("/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca")  
  load("ptyr407_workspace.rda")
  ptree <- phyloscanner.trees[["1225_to_1474"]]
  classification <-classify(ptree, allow.mt = T, n.mt=Inf, p.mt= Inf, identify.multifurcation=T, relaxed.ancestry = F,verbose = F, no.progress.bars = F)
  b <- classification[[1]]
  b <- data.table(b)
  b[host.1=='AID4452' & host.2=='AID6357']
  dwin <- data.table(dwin)
  dwin[host.1=='AID4452' & host.2=='AID6357' & grepl('1225_',tree.id)]
  
    # setwd("/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount")  
  # load("ptyr251_workspace.rda")
  # dw[H1=='AID0812' & H2=='AID1868' & grepl('1100_',TREE_ID)]
  # load("Rakai_phscnetworks.rda")
  # library(data.table)
  # dwin <- data.table(dwin)
  # dwin[host.1=="AID0812" & host.2=="AID1868" & grepl('1100_',tree.id)]
  # load("~/dcdwina_02_05_30_min_read_100_max_read_posthoccount.rda")
  # dwina[host.1=="AID0812" & host.2=="AID1868" & grepl('1100_',tree.id)]
  # library(data.table)
  # dwin <- data.table(dwin)
  # dwin[host.1=="AID1755" & host.2=="AID6465" &tree.id=="1225_to_1474"]

  
  library(data.table)
  # setwd("/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount")  
  setwd("/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca")  
  load("ptyr251_workspace.rda")
  # ptree <- phyloscanner.trees[["1100_to_1349"]]
  ptree <- phyloscanner.trees[["3600_to_3849"]]
  a <- ptree$classification.results[[1]]
  a <- data.table(a)
  # a[host.1=='AID0812' & host.2=='AID1868']
  a[host.2=='AID8610' & host.1=='AID8758']
  # classification <-classify(ptree, allow.mt = T, n.mt=Inf, p.mt= Inf, identify.multifurcation=F, relaxed.ancestry = T,verbose = F, no.progress.bars = F)
  # a <- classification[[1]]
  # a <- data.table(a)
  
  library(data.table)
  setwd("/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount")  
  load("ptyr251_workspace.rda")
  ptree <- phyloscanner.trees[["3600_to_3849"]]
  allow.mt = T
  n.mt=Inf
  p.mt= Inf
  identify.multifurcation=T
  # identify.multifurcation=F
  relaxed.ancestry = T
  verbose = F
  no.progress.bars = F
  
  library(ape)
  library(tibble)
  library(phangorn)
  library(phyloscannerR)
  library(ggtree)
  
  grep('8610',hosts.included)
  grep('8758',hosts.included)
  pat.1 <- 2
  pat.2 <- 4
  
  # desc.tips <- ptree$tips.for.hosts[[pat.2.id]]
  # desc.tips[desc.tips%in% desc.nodes[desc.nodes %in% desc.tips]]
  # desc.anc.pd <- abs(depths[desc.nodes] - depths[anc.mrca.node])
  pat.1 <- 1
  pat.2 <- 4
  pat.1 <- 3
  pat.2 <- 5
  ##
  pat.1 <- 2
  pat.2 <- 5
  
  pat.1 <- 2
  pat.2 <- 12
  
  infile.run='/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_phscinput_runs.rds'
  pty.runs <- data.table(readRDS(infile.run))
  for(i in 1:nrow(tmp)){
    tmp1 <- unique(pty.runs[grep(tmp$H1[i],RENAME_ID),]$PTY_RUN)
    tmp2 <- unique(pty.runs[grep(tmp$H2[i],RENAME_ID),]$PTY_RUN)
    tmp3 <- sort(intersect(tmp1,tmp2))[1]
    cat('H1: ',tmp$H1[i],'\n',
        'H2: ',tmp$H2[i],'\n',
        'RUN: ',tmp3,'\n')
  }
  tmp 
}

library(devtools)
install_github("BDI-pathogens/phyloscanner/phyloscannerR", ref="ordev2c" , dependencies=FALSE, INSTALL_opts="--no-staged-install")


desc.tips <- ptree$tips.for.hosts[[pat.2.id]]
anc.tips <- ptree$tips.for.hosts[[pat.1.id]]
# desc.tips.anc <- sapply(desc.tips,
#                         function(desc.tip){
#                           anc.tips <- Ancestors(tree,desc.tip)
#                           anc.tips[which(individual[anc.tips]==pat.1.id)[1]]})
# desc.tips.anc.pd <- depths[desc.tips]-depths[desc.tips.anc]
anc.tips.mrca <- mrca.phylo(tree,anc.tips)
desc.tips.mrca <- mrca.phylo(tree,desc.tips)
desc.tips.anc.pd <- depths[desc.tips]-depths[anc.tips.mrca]

downsample_ratio_seed <- function(){
  library(purrr)
  library(dplyr)
  library(phyloscannerR)
  library(data.table)

  
  indir.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships'
  # job_tag <- '_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd'
  job_tag <- '_02_05_30_min_read_null_max_read_posthoccount'
  indir	<- paste0(indir.base,job_tag)
  imbalance_read_df <- readRDS(file.path(indir,'imbalance_read_df.rds'))
  # comment to consider min=1
  # imbalance_read_df <- imbalance_read_df[reads.1>=30 & reads.2>=30]
  tmp <- imbalance_read_df[(close.and.adjacent.and.ancestry.cat=='12' & reads.1>reads.2)|
                             (close.and.adjacent.and.ancestry.cat=='21' & reads.1<reads.2)]
  nrow(tmp)
  nrow(imbalance_read_df)
  nrow(tmp)/nrow(imbalance_read_df)
  
  ans <- c()
  for (i in 1:10){
  # job_tag <- paste0('_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd_seed',i)
  job_tag <- paste0('_02_05_30_min_read_100_max_read_posthoccount_seed',i)
  indir	<- paste0(indir.base,job_tag)
  imbalance_read_df_downsampling <- readRDS(file.path(indir,'imbalance_read_df.rds'))
  # imbalance_read_df_downsampling <- imbalance_read_df_downsampling[reads.1>=30 & reads.2>=30]
  tmp_downsampling <- imbalance_read_df_downsampling[(close.and.adjacent.and.ancestry.cat=='12' & reads.1>reads.2)|
                                                       (close.and.adjacent.and.ancestry.cat=='21' & reads.1<reads.2)]
  nrow(tmp_downsampling)
  nrow(imbalance_read_df_downsampling)
  nrow(tmp_downsampling)/nrow(imbalance_read_df_downsampling)
  # whether to use posthoccount
  ans <- c(ans,nrow(tmp_downsampling)/nrow(tmp))
  
  }
  round(quantile(ans,probs = c(0.5,0.025,0.975)),4)
  # pairs with clinical data
  # tmp <- dclnid[EPID_EVIDENCE_TYPE=='serodisc']
  
  
  
  indir.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships'
  # job_tag <- '_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd'
  job_tag <- '_02_05_30_min_read_null_max_read_posthoccount'
  indir	<- paste0(indir.base,job_tag)
  imbalance_read_df <- readRDS(file.path(indir,'imbalance_read_df.rds'))
  ans <- c()
  for(i in 1:10){
    # job_tag <- paste0('_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd_seed',i)
    job_tag <- paste0('_02_05_30_min_read_100_max_read_posthoccount_seed',i)
    indir	<- paste0(indir.base,job_tag)
    imbalance_read_df_downsampling <- readRDS(file.path(indir,'imbalance_read_df.rds'))
    
    infile.sero <- '~/serodata19.rda'
    infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
    
    load(infile.sero)
    aid <- data.table(read.csv(infile.ind.anonymised))
    aid[,AID:=as.character(AID)]
    tmp <- merge(dsero, aid, by.x='MALE_RID',by.y='PT_ID')
    tmp <- merge(tmp, aid, by.x='FEMALE_RID',by.y='PT_ID')
    setnames(tmp,c('AID.x','AID.y'),c('AID1','AID2'))
    tmp <- subset(tmp,select=c('AID1','AID2','EXT_DIR'))
    tmp[,EXT_DIR := ifelse(EXT_DIR=='mf','12','21')]
    tmp2 <- tmp[!(AID1<AID2)]
    setnames(tmp2,c('AID1','AID2'),c('AID2','AID1'))
    set(tmp2, NULL, 'EXT_DIR', tmp2[, gsub('xx','12',gsub('12','21',gsub('21','xx',EXT_DIR)))])
    tmp <- rbind(tmp[(AID1<AID2)],tmp2)
    
    imbalance_read_df2 <- merge(imbalance_read_df,tmp,by.x=c('host.1','host.2'),by.y=c('AID1','AID2'))
    imbalance_read_df_downsampling2 <- merge(imbalance_read_df_downsampling,tmp,by.x=c('host.1','host.2'),by.y=c('AID1','AID2'))
    
    tmp2 <- imbalance_read_df2[(close.and.adjacent.and.ancestry.cat=='12' & reads.1>reads.2)|
                                 (close.and.adjacent.and.ancestry.cat=='21' & reads.1<reads.2)]
    nrow(tmp2)
    nrow(imbalance_read_df2)
    nrow(tmp2)/nrow(imbalance_read_df2)
    
    
    tmp_downsampling2 <- imbalance_read_df_downsampling2[(close.and.adjacent.and.ancestry.cat=='12' & reads.1>reads.2)|
                                                           (close.and.adjacent.and.ancestry.cat=='21' & reads.1<reads.2)]
    
    nrow(tmp_downsampling2)
    nrow(imbalance_read_df_downsampling2)
    nrow(tmp_downsampling2)/nrow(imbalance_read_df_downsampling2)
    
    ans <- c(ans, nrow(tmp_downsampling2)/nrow(tmp2))
  }
  round(quantile(ans,probs = c(0.5,0.025,0.975)),4)
  
  
  # couple
  indir.base <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships'
  # job_tag <- '_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd'
  job_tag <- '_02_05_30_min_read_null_max_read_posthoccount'
  indir	<- paste0(indir.base,job_tag)
  imbalance_read_df <- readRDS(file.path(indir,'imbalance_read_df.rds'))
  ans <- c()
  for(i in 1:10){
    # job_tag <- paste0('_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd_seed',i)
    job_tag <- paste0('_02_05_30_min_read_100_max_read_posthoccount_seed',i)
    indir	<- paste0(indir.base,job_tag)
    imbalance_read_df_downsampling <- readRDS(file.path(indir,'imbalance_read_df.rds'))
    
    infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
    load(infile.couple)
    couple <- data.table(subset(coupdat,select=c('male.RCCS_studyid','female.RCCS_studyid')))
    couple[,male.RCCS_studyid:=paste0('RK-',male.RCCS_studyid)]
    couple[,female.RCCS_studyid:=paste0('RK-',female.RCCS_studyid)]
    aid <- data.table(read.csv(infile.ind.anonymised))
    aid[,AID:=as.character(AID)]
    tmp <- merge(couple, aid, by.x='male.RCCS_studyid',by.y='PT_ID')
    tmp <- merge(tmp, aid, by.x='female.RCCS_studyid',by.y='PT_ID')
    setnames(tmp,c('AID.x','AID.y'),c('AID1','AID2'))
    tmp <- subset(tmp,select=c('AID1','AID2'))
    tmp2 <- tmp[!(AID1<AID2)]
    setnames(tmp2,c('AID1','AID2'),c('AID2','AID1'))
    tmp <- rbind(tmp[(AID1<AID2)],tmp2)
    
    imbalance_read_df2 <- merge(imbalance_read_df,tmp,by.x=c('host.1','host.2'),by.y=c('AID1','AID2'))
    imbalance_read_df_downsampling2 <- merge(imbalance_read_df_downsampling,tmp,by.x=c('host.1','host.2'),by.y=c('AID1','AID2'))
    
    tmp2 <- imbalance_read_df2[(close.and.adjacent.and.ancestry.cat=='12' & reads.1>reads.2)|
                                 (close.and.adjacent.and.ancestry.cat=='21' & reads.1<reads.2)]
    nrow(tmp2)
    nrow(imbalance_read_df2)
    nrow(tmp2)/nrow(imbalance_read_df2)
    
    tmp_downsampling2 <- imbalance_read_df_downsampling2[(close.and.adjacent.and.ancestry.cat=='12' & reads.1>reads.2)|
                                                           (close.and.adjacent.and.ancestry.cat=='21' & reads.1<reads.2)]
    
    nrow(tmp_downsampling2)
    nrow(imbalance_read_df_downsampling2)
    nrow(tmp_downsampling2)/nrow(imbalance_read_df_downsampling2)
    
    ans <- c(ans, nrow(tmp_downsampling2)/nrow(tmp2))
  }
  round(quantile(ans,probs = c(0.5,0.025,0.975)),4)
}

networks <- function(){
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/'
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd/'
  setwd(indir)
  load('Rakai_phscnetworks.rda')
  dclkl <- dc[CATEGORISATION=='close.and.adjacent.and.ancestry.cat']
  tmp <- dclkl[,list(MAX_N_EFF=max(N_EFF)),by=c('H1','H2')]
  dclkl <- merge(dclkl,tmp, by=c('H1','H2'))
  dclkl <- dclkl[MAX_N_EFF==N_EFF,]
  tmp <- dclkl[,list(MIN_PTY_RUN=min(PTY_RUN)),by=c('H1','H2')]
  dclkl <- merge(dclkl,tmp, by=c('H1','H2'))
  dclkl <- dclkl[MIN_PTY_RUN==PTY_RUN,]
  
  # tmp <- dclkl[,length(N_EFF),by=c('H1','H2')]
  # setkey(tmp,V1)
  # dclkl <- dclkl[,head(.SD, 4),by=c('H1','H2')]
  dclkl <- data.table(dcast(dclkl, H1+H2~TYPE,value.var='K_EFF'))
  # revisit
  # colnames(dclkl)
  dclkl[, SCORE_LINK:=(complex.or.no.ancestry+`12`+`21`)/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
  # dclkl[, SCORE_DIR:=(`12`+`21`)/(complex.or.no.ancestry+`12`+`21`)]
  # dclkl[, SCORE_12:=`12`/(`12`+`21`)]
  if(1){
    dclkl[,SCORE_12:=`12`/(`12`+`21`)]
    dclkl[,SCORE_21:=`21`/(`12`+`21`)]
    dclkl[, LINK:=(SCORE_LINK>0.6)]
    dclkl[LINK==T, DIR12:=(SCORE_12>0.6)]
    dclkl[LINK==T, DIR21:=(SCORE_21>0.6)]
  }
  
  
  # dclkl[LINK==T, DIR:=(SCORE_DIR>0.6)]
  # dclkl[DIR==T, DIR12:=(SCORE_12>0.5)]
  # nrow(dclkl[LINK==F])
  # nrow(dclkl[LINK==T])
  # dclkl[LINK==T, table(DIR)]
  # hist(dclkl[LINK==T, SCORE_DIR])
  # nrow(dclkl[DIR==T])
  
  nrow(dclkl[LINK==T])
  nrow(dclkl[LINK==T & (DIR12==T|DIR21==T)])
  
  dpairs <- dclkl[LINK==T & (DIR12==T|DIR21==T)]
  tmp <- dpairs[DIR12==F]
  setnames(tmp, c('H1','H2','12','21','SCORE_12','SCORE_21','DIR12','DIR21'),
           c('H2','H1','21','12','SCORE_21','SCORE_12','DIR21','DIR12'))
  dpairs <- rbind(dpairs[DIR12==T],tmp)
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  aid <- data.table(read.csv(infile.ind.anonymised))
  aid <- subset(aid,select=c('PT_ID','AID'))
  # dpairs <- merge(dpairs, aid, by.x='H1', by.y='AID',all.x=T)
  # dpairs <- merge(dpairs, aid, by.x='H2', by.y='AID',all.x=T)
  
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  infile.ind.rccs <- file.path(data.dir,'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.mrc <- file.path(data.dir,'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
  id.dt <- data.table(read.csv(infile.ind.rccs))
  id.dt <- subset(id.dt,select = c("pt_id","sex","age_enrol","geo_country","ever_art"))
  tmp <- data.table(read.csv(infile.ind.mrc))
  tmp <- subset(tmp,select = c("pt_id","sex","age_enrol","geo_country","ever_art"))
  id.dt <- rbind(id.dt,tmp)
  id.dt <- unique(id.dt)
  id.dt <- merge(id.dt,aid,by.y='PT_ID',by.x='pt_id')
  dpairs <- merge(dpairs,id.dt, by.x=c('H1'),by.y='AID',all.x=T)
  dpairs <- merge(dpairs,id.dt, by.x=c('H2'),by.y='AID',all.x=T)
  dpairs[,RK.x:=grepl('RK-',pt_id.x)]
  dpairs[,RK.y:=grepl('RK-',pt_id.y)]
  table(dpairs$sex.x,dpairs$sex.y)
  dpairs[sex.x%in%c('M','F') & sex.y%in%c('M','F')]
  dpairs[sex.x%in%c('M','F') & sex.y%in%c('M','F'),sum(sex.x!=sex.y)]
  dpairs[sex.x%in%c('M','F') & sex.y%in%c('M','F')&sex.x!=sex.y,sum(sex.x=='M')]
  dpairs[sex.x%in%c('M','F') & sex.y%in%c('M','F')&sex.x!=sex.y,sum(sex.x=='F')]
  
  table(dpairs$RK.x,dpairs$RK.y)
  length(unique(c(dpairs$H1,dpairs$H2)))
  
}

compare_with_2019 <- function(){
    library(data.table)
    infile.sero <- '~/serodata19.rda'
    infile.aid <- '~/todi_pairs_171122_cl25_d50_prior23_min30_anonymised_RIDs.csv'
    infile.networks1 <-  '~/Rakai_phscnetworks_190706.rda'
    indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/'
    # indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd/'
    
    infile.networks2 <-  file.path(indir2,'Rakai_phscnetworks.rda')
    infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
    
    load(infile.sero)
    aid <- data.table(read.csv(infile.ind.anonymised))
    aid19 <-  data.table(read.csv(infile.aid))
    aid19 <- subset(aid19,select=c('ID','AID'))
    aid19[,ID:=paste0('RK-',ID)]
    aid19 <- merge(aid19, subset(aid,select=c('PT_ID','AID')), by.x='ID', by.y='PT_ID',all.x=T)
    aid19 <- subset(aid19, select=c('AID.x','AID.y'))
    
    
    load(infile.networks1)
    dc1 <- copy(dc)
    dw1 <- copy(dw)
    dc1 <- merge(dc1, aid19, by.x='H1', by.y='AID.x',all.x=T)
    dc1 <- merge(dc1, aid19, by.x='H2', by.y='AID.x',all.x=T)
    set(dc1,NULL,c('H1','H2'),NULL)
    setnames(dc1, c("AID.y.x","AID.y.y"),c('H1','H2'))
    dc1 <- as.data.table(dc1)
    dw1 <- merge(dw1, aid19, by.x='H1', by.y='AID.x',all.x=T)
    dw1 <- merge(dw1, aid19, by.x='H2', by.y='AID.x',all.x=T)
    set(dw1,NULL,c('H1','H2'),NULL)
    setnames(dw1, c("AID.y.x","AID.y.y"),c('H1','H2'))
    dw1 <- as.data.table(dw1)
    
       
    
    
    
    load(infile.networks2)
    dw2 <- copy(dw)
    dc2 <- copy(dc)
    
    aid[,X:=NULL]
    dsero <- merge(dsero,aid, by.x='FEMALE_RID',by.y='PT_ID',all.x=T)
    dsero <- merge(dsero,aid, by.x='MALE_RID',by.y='PT_ID',all.x=T)
    dsero[,AID.x:=as.character(AID.x)]
    dsero[,AID.y:=as.character(AID.y)]
    dsero[,SERO_DIR:=ifelse(EXT_DIR=='fm','12','21')]
    dsero <- subset(dsero,select=c('AID.x','AID.y','SERO_DIR'))
    
    tmp  <- dsero[!(AID.x<AID.y)]
    setnames(tmp,c('AID.x','AID.y'),c('AID.y','AID.x'))
    tmp[,SERO_DIR:=ifelse(SERO_DIR=='12','21','12')]
    dsero <- rbind(dsero[(AID.x<AID.y)],tmp)
    
    
    dclkl1 <- dc1[CATEGORISATION=='close.and.adjacent.and.ancestry.cat']
    tmp <- dclkl1[,list(MAX_N_EFF=max(N_EFF)),by=c('H1','H2')]
    dclkl1 <- merge(dclkl1,tmp, by=c('H1','H2'))
    dclkl1 <- dclkl1[MAX_N_EFF==N_EFF,]
    tmp <- dclkl1[,length(N_EFF),by=c('H1','H2')]
    setkey(tmp,V1)
    dclkl1 <- dclkl1[,head(.SD, 4),by=c('H1','H2')]
    dclkl1 <- data.table(dcast(dclkl1, H1+H2~TYPE,value.var='K_EFF'))
    dclkl1[,H1:=as.character(H1)]
    dclkl1[,H2:=as.character(H2)]
    tmp <- dclkl1[H1>H2]
    setnames(tmp, c('H1','H2','12','21'),c('H2','H1','21','12'))
    dclkl1 <- rbind(dclkl1[H1<H2],tmp)
      
    dclkl1[, SCORE_LINK:=(complex.or.no.ancestry+`12`+`21`)/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
    # dclkl1[, SCORE_DIR:=(`12`+`21`)/(complex.or.no.ancestry+`12`+`21`)]
    # dclkl1[, SCORE_12:=`12`/(`12`+`21`)]
    if(1){
      dclkl1[, SCORE_12:=`12`/(`12`+`21`)]
      dclkl1[, SCORE_21:=`21`/(`12`+`21`)]
      dclkl1[, LINK:=(SCORE_LINK>0.6)]
      # dclkl1[LINK==T, DIR:=(SCORE_DIR>0.6)]
      # dclkl1[DIR==T, DIR12:=(SCORE_12>0.5)]
      dclkl1[LINK==T, DIR12:=(SCORE_12>0.6)]
      dclkl1[LINK==T, DIR21:=(SCORE_21>0.6)]
    }
    if(0){
      dclkl1[, SCORE_12:=`12`/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
      dclkl1[, SCORE_21:=`21`/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
      dclkl1[, LINK:=(SCORE_LINK>0.5)]
      dclkl1[LINK==T, DIR12:=(SCORE_12>0.33)]
      dclkl1[LINK==T, DIR21:=(SCORE_21>0.33)]
    }
    
    
    dclkl2 <- dc2[CATEGORISATION=='close.and.adjacent.and.ancestry.cat']
    tmp <- dclkl2[,list(MAX_N_EFF=max(N_EFF)),by=c('H1','H2')]
    dclkl2 <- merge(dclkl2,tmp, by=c('H1','H2'))
    dclkl2 <- dclkl2[MAX_N_EFF==N_EFF,]
    tmp <- dclkl2[,length(N_EFF),by=c('H1','H2')]
    setkey(tmp,V1)
    dclkl2 <- dclkl2[,head(.SD, 4),by=c('H1','H2')]
    dclkl2 <- data.table(dcast(dclkl2, H1+H2~TYPE,value.var='K_EFF'))
    dclkl2[, SCORE_LINK:=(complex.or.no.ancestry+`12`+`21`)/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
    # dclkl2[, SCORE_DIR:=(`12`+`21`)/(complex.or.no.ancestry+`12`+`21`)]
    # dclkl2[, SCORE_12:=`12`/(`12`+`21`)]
    if(1){
      dclkl2[, SCORE_12:=`12`/(`12`+`21`)]
      dclkl2[, SCORE_21:=`21`/(`12`+`21`)]
      dclkl2[, LINK:=(SCORE_LINK>0.6)]
      # dclkl2[LINK==T, DIR:=(SCORE_DIR>0.6)]
      # dclkl2[DIR==T, DIR12:=(SCORE_12>0.5)]
      dclkl2[LINK==T, DIR12:=(SCORE_12>0.6)]
      dclkl2[LINK==T, DIR21:=(SCORE_21>0.6)]
    }
    if(0){
      dclkl2[, SCORE_12:=`12`/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
      dclkl2[, SCORE_21:=`21`/(not.close.or.nonadjacent+complex.or.no.ancestry+`12`+`21`)]
      dclkl2[, LINK:=(SCORE_LINK>0.5)]
      dclkl2[LINK==T, DIR12:=(SCORE_12>0.33)]
      dclkl2[LINK==T, DIR21:=(SCORE_21>0.33)]
    }
    dclkl <- merge(dclkl1,dclkl2,by=c('H1','H2'),all=T)
    dclkl <- merge(dclkl, dsero, by.x=c('H1','H2'),by.y=c('AID.x','AID.y'))
    # dclkl[LINK.x!=LINK.y]
    # dclkl[DIR.x!=DIR.y]
    # dclkl[DIR.x!=DIR.y,(DIR12.x==T & SERO_DIR=='12')|(DIR12.x==F & SERO_DIR=='21')]
    # # 2 false 7 true
    # dclkl[DIR12.x!=DIR12.y]
    # # phyloscan for these pairs
    
    # tmp <- dclkl[LINK.x==T & LINK.y==T]
    # tmp[ (DIR12.x==T & SERO_DIR=='12')|(DIR12.x==F & SERO_DIR=='21'),TYPE.x:='consistent']
    # tmp[ !((DIR12.x==T & SERO_DIR=='12')|(DIR12.x==F & SERO_DIR=='21')),TYPE.x:='inconsistent']
    # tmp[DIR.x==F, TYPE.x:='unresolved']
    # tmp[ (DIR12.y==T & SERO_DIR=='12')|(DIR12.y==F & SERO_DIR=='21'),TYPE.y:='consistent']
    # tmp[ !((DIR12.y==T & SERO_DIR=='12')|(DIR12.y==F & SERO_DIR=='21')),TYPE.y:='inconsistent']
    # tmp[DIR.y==F, TYPE.y:='unresolved']
    # tmp[, table(TYPE.x, TYPE.y, useNA='if')]
    # tmp[, table(TYPE.x, useNA='if')]
    # tmp[, table(TYPE.y, useNA='if')]
    tmp <- dclkl[LINK.x==T & LINK.y==T]
    tmp[ (DIR12.x==T & SERO_DIR=='12')|(DIR21.x==T & SERO_DIR=='21'),TYPE.x:='consistent']
    tmp[ (DIR12.x==T & SERO_DIR=='21')|(DIR21.x==T & SERO_DIR=='12'),TYPE.x:='inconsistent']
    tmp[is.na(TYPE.x), TYPE.x:='unresolved']
    tmp[ (DIR12.y==T & SERO_DIR=='12')|(DIR21.y==T & SERO_DIR=='21'),TYPE.y:='consistent']
    tmp[ (DIR12.y==T & SERO_DIR=='21')|(DIR21.y==T & SERO_DIR=='12'),TYPE.y:='inconsistent']
    tmp[is.na(TYPE.y), TYPE.y:='unresolved']
    tmp[, table(TYPE.x, TYPE.y, useNA='if')]
    tmp[, table(TYPE.x, useNA='if')]
    tmp[, table(TYPE.y, useNA='if')]
    # tmp <- tmp[TYPE.y=='inconsistent']
    # dclkl1 <- merge(dclkl1, dsero, by.x=c('H1','H2'),by.y=c('AID.x','AID.y'))
    # dclkl2 <- merge(dclkl2, dsero, by.x=c('H1','H2'),by.y=c('AID.x','AID.y'))
    # dclkl1[(DIR12==T & SERO_DIR=='12')|(DIR12==F & SERO_DIR=='21') ]
    # dclkl2[(DIR12==T & SERO_DIR=='12')|(DIR12==F & SERO_DIR=='21') ]
    # 
    # dclkl1[!((DIR12==T & SERO_DIR=='12')|(DIR12==F & SERO_DIR=='21')) ]
    # dclkl2[!((DIR12==T & SERO_DIR=='12')|(DIR12==F & SERO_DIR=='21')) ]
    # tmp <- dclkl[DIR.x!=DIR.y]
    # dir.create('~/compare_multitrans_nodownsampling/')
    
    tmp <- dclkl[DIR12.x!=DIR12.y|DIR21.x!=DIR21.y]
    
}
