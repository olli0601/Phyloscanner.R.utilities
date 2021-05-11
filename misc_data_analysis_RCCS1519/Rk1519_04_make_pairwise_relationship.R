
delete_nonexcised <- function(){
  require(data.table)
  
  # files
  HOME <<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  treedir <- file.path(HOME,'210325_phsc_output/')
  df <- data.table(F=list.files(treedir,recursive = T, pattern = '*.treefile'))
  # which to delete
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
  # delete
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
  
  HOME <<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  treedir <- file.path(HOME,'210325_phsc_output/')
  tmpdir	<- file.path(HOME,"210325_phsc_work/")
  # outdir	<- file.path(HOME,'210325_phsc_phscrelationships_02_05_null_min_read/') #no min no max
  # outdir	<- file.path(HOME,'210325_phsc_phscrelationships_02_05_30_min_read_100_max_read/') #30 min 100 max
  outdir	<- file.path(HOME,'210325_phsc_phscrelationships_02_05_null_min_read_100_max_read/') #no min 100 max
  dir.create(outdir)
  prog.phyloscanner_analyse_trees <- '/rds/general/user/xx4515/home/phyloscanner/phyloscanner_analyse_trees.R'
  
  #	set phyloscanner variables	
  control	<- list()
  control$allow.mt <- TRUE				
  control$alignment.file.directory = NULL 
  control$alignment.file.regex = NULL
  control$blacklist.underrepresented = FALSE	
  control$count.reads.in.parsimony = TRUE
  control$distance.threshold <- '0.02 0.05'
  control$do.dual.blacklisting = FALSE					
  control$duplicate.file.directory = NULL
  control$duplicate.file.regex = NULL
  control$file.name.regex = "^(?:.*\\D)?([0-9]+)_to_([0-9]+).*$"
  control$guess.multifurcation.threshold = FALSE
  # control$max.reads.per.host <- NULL
  control$min.reads.per.host <- NULL
  control$max.reads.per.host <- 100
  # control$min.reads.per.host <- 30
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
  # control$post.hoc.count.blacklisting <- TRUE
  control$post.hoc.count.blacklisting <- NULL
  control$ratio.blacklist.threshold = 0.01
  control$raw.blacklist.threshold = 3			
  control$recombination.file.directory = NULL
  control$recombination.file.regex = NULL
  control$relaxed.ancestry = TRUE
  control$sankoff.k = 15
  control$sankoff.unassigned.switch.threshold = 0
  control$seed = 42
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
  hpc.walltime		<- 923						# walltime
  if(0)		
  {
    hpc.q			<- NA						# PBS queue
    hpc.mem			<- "36gb" 					# RAM		
  }
  #		or run this block to submit a job array to Oliver's machines
  if(1)
  {
    hpc.q			<- "pqeelab"				# PBS queue
    hpc.mem			<- "12gb" 					# RAM		
  }
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
  
}

phsc.transmission.networks<- function()
{
  library(data.table)
  library(tidyverse)
  source('~/transmission_network_functions_phyloscanner.R')
  # optional: meta data
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/' #no min no max
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read/' #30 min 100 max
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read_100_max_read/' #no min 100 max
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
  control <- list(linked.group='close.and.adjacent.cat',linked.no='not.close.or.nonadjacent',linked.yes='close.and.adjacent', conf.cut=0.6, neff.cut=3,weight.complex.or.no.ancestry=0.5)
  # save(dca,dwina,file = '~/dcdwina_minreadnull.rda')
  # save(dca,dwina,file = '~/dcdwina_minread30_maxread100.rda')
  # save(dca,dwina,file = '~/dcdwina_minreadnull_maxread100.rda')
  # load('~/dcdwina_minreadnull.rda')
  # load('~/dcdwina_minread30_maxread100.rda')
  # load('~/dcdwina_minreadnull_maxread100.rda')
  
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


analysis_network <- function(){
  library(data.table)
  # dirs
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/'
  # indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read/' #30 min 100 max
  # indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read_100_max_read/' #no min 100 max
  # 

  # load 
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  infile.ind.rccs <- file.path('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.mrc <- file.path('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
  infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
  load( infile.networks )	# loads rtp, rplkl, rpw
  
  # meta
  df <- data.table(read.csv(infile.ind.anonymised))
  tmp <- rbind(data.table(read.csv(infile.ind.rccs)),data.table(read.csv(infile.ind.mrc)))
  df2 <- merge(subset(df,select=c('PT_ID', 'AID')),subset(tmp,select=c('pt_id','sex')),by.x='PT_ID',by.y='pt_id',all.x=T)
  df2 <- subset(df2,select=c('AID','sex'))
  colnames(df2) <- c('ID','SEX')
  df2 <- unique(df2)
  
  # classification
  conf.cut		<- 0.6			
  dchain <- as.data.table(dchain)
  dchain[, PHYLOSCANNER_CLASSIFY:= NA_character_]
  set(dchain, dchain[, which(is.na(PTY_RUN))], 'PHYLOSCANNER_CLASSIFY', 'insufficient deep sequence data for at least one individual')
  set(dchain, dchain[, which(!is.na(PTY_RUN) & is.na(LINK_12) & is.na(LINK_21))], 'PHYLOSCANNER_CLASSIFY', 'ph unlinked pair')
  set(dchain, dchain[, which(!is.na(SCORE_LINKED) & SCORE_LINKED<=conf.cut)], 'PHYLOSCANNER_CLASSIFY', 'unclear if pair ph linked or unlinked')
  set(dchain, dchain[, which(!is.na(SCORE_LINKED) & SCORE_LINKED>conf.cut)], 'PHYLOSCANNER_CLASSIFY', 'ph linked pair direction not resolved')
  set(dchain, dchain[, which(!is.na(SCORE_LINKED) & SCORE_LINKED>conf.cut & SCORE_DIR_12>conf.cut)], 'PHYLOSCANNER_CLASSIFY', 'ph linked pair direction 12')
  set(dchain, dchain[, which(!is.na(SCORE_LINKED) & SCORE_LINKED>conf.cut & SCORE_DIR_21>conf.cut)], 'PHYLOSCANNER_CLASSIFY', 'ph linked pair direction 21')
  
  # if several classifications per pair, take the majority
  tmp <- dchain[,length(unique( PHYLOSCANNER_CLASSIFY)),by=c('H1','H2')]
  tmp <- tmp[V1>1]
  tmp[, H:=paste0(H1,'-',H2)]
  tmpm <- tmp$H
  tmp2 <- dchain[paste0(H1,'-',H2)%in% tmpm]
  tmp <- tmp2[,length(SCORE_LINKED),by=c('H1','H2','PHYLOSCANNER_CLASSIFY')]
  setkey(tmp,H1,H2,V1)
  tmp <- tmp[,list(PHYLOSCANNER_CLASSIFY=last(PHYLOSCANNER_CLASSIFY)),by=c('H1','H2')]
  tmp2 <- merge(tmp2, tmp, by=c('H1','H2','PHYLOSCANNER_CLASSIFY'))
  dchain <- dchain[!(paste0(H1,'-',H2)%in% tmpm)]
  dchain <- rbind(dchain,tmp2)
  # 
  dchain <- unique(subset(dchain,select = c('H1','H2','PHYLOSCANNER_CLASSIFY')))
  
  # get genders
  setnames(df2,colnames(df2),paste0(colnames(df2),'1'))
  dchain	<- merge(dchain, df2, by.x='H1',by.y=c('ID1'),all.x=T)
  setnames(df2,colnames(df2),gsub('1','2',colnames(df2)))
  dchain	<- merge(dchain, df2, by.x='H2',by.y=c('ID2'),all.x=T)
  setnames(df2,colnames(df2),gsub('2','',colnames(df2)))
  
  # get couple data
  load(infile.couple)
  couple <- unique(subset(data.table(coupdat), select=c('male.RCCS_studyid', 'female.RCCS_studyid')))
  couple[, COUPLE:=1]
  setnames(couple,c('male.RCCS_studyid','female.RCCS_studyid'),c('PT_ID1','PT_ID2'))
  couple[,PT_ID1:=paste0('RK-',PT_ID1)]
  couple[,PT_ID2:=paste0('RK-',PT_ID2)]
  couple <- merge(couple,subset(df,select=c('PT_ID', 'AID')),by.x='PT_ID1', by.y='PT_ID', all.x=T)
  couple <- merge(couple,subset(df,select=c('PT_ID', 'AID')),by.x='PT_ID2', by.y='PT_ID', all.x=T)
  set(couple,NULL,c('PT_ID1','PT_ID2'),NULL)
  setnames(couple,c('AID.x','AID.y'),c('H1','H2'))
  couple[,H1:=as.character(H1)]
  couple[,H2:=as.character(H2)]
  couple <- couple[!is.na(H1) & !is.na(H2)]
  
  # order ID
  setnames(dchain,c('H1','H2','SEX1','SEX2'),c('ID1','ID2','ID1_SEX','ID2_SEX'))
  tmp		<- subset(dchain, ID1_SEX=='F' & ID2_SEX=='M')
  setnames(tmp, colnames(tmp), gsub('xx','ID2',gsub('ID2','ID1',gsub('ID1','xx',gsub('xx','12',gsub('12','21',gsub('21','xx',colnames(tmp))))))))
  set(tmp, NULL, 'PHYLOSCANNER_CLASSIFY', tmp[, gsub('xx','12',gsub('12','21',gsub('21','xx',PHYLOSCANNER_CLASSIFY)))])
  dchain	<- rbind(subset(dchain, !(ID1_SEX=='F' & ID2_SEX=='M')), tmp)
  dchain[, PAIR_SEX:= paste0(ID1_SEX,ID2_SEX)]
  unique(dchain$PHYLOSCANNER_CLASSIFY)
  
  # directions
  rtp		<- subset(dchain, !grepl('not resolved|unlinked',PHYLOSCANNER_CLASSIFY))
  # tmp <- rtp[,length(unique( PHYLOSCANNER_CLASSIFY)),by=c('ID1','ID2')]
  # tmp[V1>1]
  rtp[, table(PAIR_SEX)]
  # FF  FU  MF  MM  MU  UF  UM 
  # 95   7 448 113  11   8   5 
  
  # FF  FU  MF  MM  MU  UF  UM 
  # 146   7 593 137  11  11   7 
  # 
  # check couple
  rtp <- merge(rtp, couple, by.x=c('ID1','ID2'),by.y=c('H1','H2'),all.x=T)
  rtp[!is.na(COUPLE)]
  dchain <- merge(dchain, couple, by.x=c('ID1','ID2'),by.y=c('H1','H2'),all.x=T)  
  cat(nrow(dchain[!is.na(COUPLE) & !grepl('unlinked',PHYLOSCANNER_CLASSIFY)]), ' linked couple \n',
      nrow(dchain[!is.na(COUPLE) & !grepl('not resolved|unlinked',PHYLOSCANNER_CLASSIFY)]), ' linked couple with directions \n',
      nrow(dchain[!is.na(COUPLE)]),' couple \n')
  
  # MIN READ 30 - 
  # 215  linked couple 
  # 168  linked couple with directions 
  # 227  couple 
  # MIN READ NULL - 
  # 291  linked couple 
  # 226  linked couple with directions 
  # 314  couple 
  
  
  # check RCCS
  dobs <- copy(rtp)
  df3 <- subset(df, select=c('PT_ID','AID'))
  df3[,RCCS:=grepl('RK-',PT_ID)]
  df3 <-subset(df3,select=c('AID','RCCS'))
  dobs <- merge(dobs, df3, by.x='ID1',by.y='AID', all.x=T) 
  dobs <- merge(dobs, df3, by.x='ID2',by.y='AID', all.x=T) 
  setnames(dobs,c('RCCS.x','RCCS.y'),c('TR_RCCS','REC_RCCS'))
  
  # get heterosexual 
  tmp <- dobs[ID1_SEX=='M' & ID2_SEX=='F']
  cat(nrow(tmp), ' heterosexual pairs with directions \n ',
      nrow(tmp[TR_RCCS==T & REC_RCCS==T]), ' RCCS')
  # nrow(tmp[(TR_RCCS==T & REC_RCCS==F)|TR_RCCS==F & REC_RCCS==T ])
  
  # 448  heterosexual pairs with directions 
  # 403  RCCS
  # 593  heterosexual pairs with directions 
  # 539  RCCS
}

couple.analysis <- function(){
  library(data.table)
  # load infiles
  infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
  infile.run='/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_phscinput_runs.rds'
  pty.runs <- data.table(readRDS(infile.run))
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  infile.ind.rccs <- file.path(data.dir,'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.mrc <- file.path(data.dir,'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
  infile.aid19 <- '~/todi_pairs_171122_cl25_d50_prior23_min30_anonymised_RIDs.csv'
  infile.pairs19 <- '~/Rakai_phscnetworks_allpairs_190706.rda'
  infile.net19 <- '~/Rakai_phscnetworks_190706.rda'
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/'
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  
  
  # map id and sequences
  id.dt <- data.table(read.csv(infile.ind.rccs))
  id.dt <- subset(id.dt,select = c("pt_id","pangea_id"))
  tmp <- data.table(read.csv(infile.ind.mrc))
  tmp <- subset(tmp,select = c("pt_id","pangea_id"))
  id.dt <- rbind(id.dt,tmp)
  id.dt <- unique(id.dt)
  
  # load couple
  load(infile.couple)
  couple <- data.table(unique(subset(coupdat,select=c('male.RCCS_studyid', 'female.RCCS_studyid'))))
  colnames(couple) <- c('ID1','ID2')
  couple[,COUPLE:=1]
  couple[,ID1:=paste0('RK-',ID1)]
  couple[,ID2:=paste0('RK-',ID2)]
  
  
  # whether couples in run
  cat(nrow(couple), ' couples in total \n',
      nrow(couple[ID1%in% pty.runs$UNIT_ID & ID2%in% pty.runs$UNIT_ID]), ' in runs \n',
      nrow(couple[ID1%in% id.dt$pt_id & ID2%in% id.dt$pt_id]), ' have pangea ids \n')
  tmp <- couple[ID1%in% id.dt$pt_id & ID2%in% id.dt$pt_id & !(ID1%in% pty.runs$UNIT_ID & ID2%in% pty.runs$UNIT_ID)]
  
  # eg
  # WARNING: lose some couples ... check later
  tmp[ID1%in% pty.runs$UNIT_ID & !ID2%in% pty.runs$UNIT_ID]
  merge(tmp[ID1%in% pty.runs$UNIT_ID & !ID2%in% pty.runs$UNIT_ID],id.dt, by.x='ID2',by.y='pt_id',all.x=T)
  tmp[!ID1%in% pty.runs$UNIT_ID & ID2%in% pty.runs$UNIT_ID]
  tmp[!ID1%in% pty.runs$UNIT_ID & !ID2%in% pty.runs$UNIT_ID]
  
  
  # load 19 network
  aid19 <- data.table(read.csv(infile.aid19))
  aid19[, ID:=paste0('RK-',ID)]
  load(infile.pairs19)
  load(infile.net19)
  dc19 <- copy(dc)
  dw19 <- copy(dw)

  
  tmp <- subset(aid19,select=c('ID','AID'))
  dc19 <- merge(dc19, tmp, by.x='H1',by.y='AID',all.x=T)
  dc19 <- merge(dc19, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dc19, c('ID.x','ID.y'),c('ID1','ID2'))
  dw19 <- merge(dw19, tmp, by.x='H1',by.y='AID',all.x=T)
  dw19 <- merge(dw19, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dw19, c('ID.x','ID.y'),c('ID1','ID2'))
  
  # load network
  load( infile.networks )
  aid <- data.table(read.csv(infile.ind.anonymised))

  tmp <- subset(aid,select=c('PT_ID','AID'))
  dc <- merge(dc, tmp, by.x='H1',by.y='AID',all.x=T)
  dc <- merge(dc, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dc, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  dw <- merge(dw, tmp, by.x='H1',by.y='AID',all.x=T)
  dw <- merge(dw, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dw, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  
  # whether couples in dc
  tmp <- unique(subset(dc,select = c('ID1','ID2')))
  tmp[,ID.x:=paste0(ID1,'-',ID2)]
  tmp[,ID.y:=paste0(ID2,'-',ID1)]
  cat(nrow(couple[paste0(ID1,'-',ID2)%in% tmp$ID.x | paste0(ID1,'-',ID2)%in% tmp$ID.y ]), ' couples in networks')
  
  # to dt
  dc19 <- as.data.table(dc19)
  dw19 <- as.data.table(dw19)
  
  # add couple
  tmp <- copy(couple)
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  couple <- data.table(unique(rbind(couple, tmp)))
  
  dc <- merge(dc,couple,by=c('ID1','ID2'),all.x=T)
  dc19 <- merge(dc19,couple,by=c('ID1','ID2'),all.x=T)
  dw <- merge(dw,couple,by=c('ID1','ID2'),all.x=T)
  dw19 <- merge(dw19,couple,by=c('ID1','ID2'),all.x=T)
  
  #
  tmp_dc_coup <- unique(dc[CATEGORISATION=='close.and.adjacent.and.ancestry.cat'&COUPLE==1])
  tmp_dc_coup[,length(N),by=c('ID1','ID2')]
  tmp <- tmp_dc_coup[,list(N_EFF=max(N_EFF)),by=c('ID1','ID2')]
  tmp_dc_coup <- merge(tmp_dc_coup,tmp, by=c('ID1','ID2','N_EFF'))
  setkey(tmp_dc_coup, ID1, ID2, PTY_RUN, TYPE)
  tmp_dc_coup <- tmp_dc_coup[,head(.SD, 4),by=c('ID1','ID2')]
  tmp_dc_coup <- dcast(tmp_dc_coup,  ID1 + ID2 + COUPLE ~ TYPE, value.var= 'SCORE')
  tmp_dc_coup <- as.data.table(tmp_dc_coup)
  tmp_dc_coup[,ANALYSIS:=2021]
  
  tmp_dc19_coup <- unique(dc19[CATEGORISATION=='close.and.adjacent.and.ancestry.cat'&COUPLE==1])
  # tmp_dc19_coup[,length(N),by=c('ID1','ID2')]
  # tmp <- tmp_dc19_coup[,list(N_EFF=max(N_EFF)),by=c('ID1','ID2')]
  # tmp_dc19_coup <- merge(tmp_dc19_coup,tmp, by=c('ID1','ID2','N_EFF'))
  tmp_dc19_coup <- dcast(tmp_dc19_coup,  ID1 + ID2 + COUPLE ~ TYPE, value.var= 'SCORE')
  tmp_dc19_coup <- as.data.table(tmp_dc19_coup)
  tmp_dc19_coup[,ANALYSIS:=2019]
  
  # order
  tmp <- tmp_dc_coup[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2','12','21'), c('ID2','ID1','21','12'))
  tmp_dc_coup <- rbind(tmp_dc_coup[(ID1<ID2)], tmp)
  tmp <- tmp_dc19_coup[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2','12','21'), c('ID2','ID1','21','12'))
  tmp_dc19_coup <- rbind(tmp_dc19_coup[(ID1<ID2)], tmp)
  tmp_dc_coup <- rbind(tmp_dc19_coup, tmp_dc_coup)
  tmp_dc_coup <- melt(tmp_dc_coup,id.vars=c('ID1','ID2', 'COUPLE','ANALYSIS'))
  
  library(dplyr)
  library(ggplot2)
  df <- arrange(tmp_dc_coup, variable, desc(value))
  
  df$ID=paste0(df$ID1,'-',df$ID2)
  df$ID=factor(df$ID, levels = unique(df$ID))
  ggplot(df, aes(ID,value,fill=variable))+
    geom_bar(position="stack", stat="identity")+
    facet_grid(ANALYSIS~ .)+
    labs(x='ID',y='scores',fill='')+ guides(fill = guide_legend(nrow = 1)) + theme(legend.position = "bottom",legend.direction = 'horizontal') +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    scale_x_discrete(limits = rev(levels(df$ID)))
  ggsave(filename = '~/compare_couple.pdf',width = 10, height = 8)
  
  
  tmp_dc_coup <- dcast(tmp_dc_coup,ID1+ID2+COUPLE+variable ~ANALYSIS ,value.var='value')
  tmp_dc_coup <- as.data.table(tmp_dc_coup)
  tmp1 <- tmp_dc_coup[is.na(`2019`) & !is.na(`2021`) ]
  tmp2 <- tmp_dc_coup[!is.na(`2019`) & is.na(`2021`) ]
  tmp3 <- tmp_dc_coup[!is.na(`2019`) & !is.na(`2021`) ]
  nrow(unique(subset(tmp1,select=c('ID1','ID2'))))
  nrow(unique(subset(tmp2,select=c('ID1','ID2'))))
  nrow(unique(subset(tmp3,select=c('ID1','ID2'))))
  # 2021 ONLY
  tmp <- dcast(tmp1,ID1+ID2~variable,value.var='2021')
  tmp <- as.data.table(tmp)
  tmp[,LINKED:= `12`+`21` +complex.or.no.ancestry>0.6]
  tmp[LINKED==1,DIRECTED:= (`12`+`21`)/(`12`+`21` +complex.or.no.ancestry)>0.6]
  tmp[DIRECTED==1,LINK12:=`12`>`21`]
  nrow(tmp[LINKED==1])/nrow(tmp)
  nrow(tmp[DIRECTED==1])/nrow(tmp[LINKED==1])
  # 2019 ONLY
  tmp <- dcast(tmp2,ID1+ID2~variable,value.var='2019')
  tmp <- as.data.table(tmp)
  tmp[,LINKED:= `12`+`21` +complex.or.no.ancestry>0.6]
  tmp[LINKED==1,DIRECTED:= (`12`+`21`)/(`12`+`21` +complex.or.no.ancestry)>0.6]
  tmp[DIRECTED==1,LINK12:=`12`>`21`]
  nrow(tmp[LINKED==1])/nrow(tmp)
  nrow(tmp[DIRECTED==1])/nrow(tmp[LINKED==1])
  tmp <- unique(subset(tmp2,select=c('ID1','ID2')))
  tmp[ID1%in% pty.runs$UNIT_ID & ID2%in% pty.runs$UNIT_ID]
  tmp[,ID.x:=paste0(ID1,'-',ID2)]
  tmp[,ID.y:=paste0(ID2,'-',ID1)]
  dw[paste0(ID1,'-',ID2) %in% tmp$ID.x]
  dw[paste0(ID1,'-',ID2) %in% tmp$ID.y]
  
  # 
  # 2019 and 2021
  tmp21 <- dcast(tmp3,ID1+ID2~variable,value.var='2021')
  tmp21 <- as.data.table(tmp21)
  tmp21[,LINKED_SCORE:=`12`+`21` +complex.or.no.ancestry]
  tmp21[,LINKED:= LINKED_SCORE>0.6]
  tmp21[LINKED==1,DIRECTED_SCORE:=(`12`+`21`)/LINKED_SCORE]
  tmp21[LINKED==1,DIRECTED:= DIRECTED_SCORE>0.6]
  tmp21[DIRECTED==1,LINK12_SCORE:=`12`/(`12`+`21`)]
  tmp21[DIRECTED==1,LINK12:=`12`>`21`]
  nrow(tmp21[LINKED==1])/nrow(tmp21)
  nrow(tmp21[DIRECTED==1])/nrow(tmp21[LINKED==1])
  
  tmp19 <- dcast(tmp3,ID1+ID2~variable,value.var='2019')
  tmp19 <- as.data.table(tmp19)
  tmp19[,LINKED_SCORE:=`12`+`21` +complex.or.no.ancestry]
  tmp19[,LINKED:= LINKED_SCORE>0.6]
  tmp19[LINKED==1,DIRECTED_SCORE:=(`12`+`21`)/LINKED_SCORE]
  tmp19[LINKED==1,DIRECTED:= DIRECTED_SCORE>0.6]
  tmp19[DIRECTED==1,LINK12_SCORE:=`12`/(`12`+`21`)]
  tmp19[DIRECTED==1,LINK12:=`12`>`21`]
  nrow(tmp19[LINKED==1])/nrow(tmp19)
  nrow(tmp19[DIRECTED==1])/nrow(tmp19[LINKED==1])
  
  tmp <- merge(tmp19,tmp21,by=c('ID1','ID2'))
  tmp[,LINKED_SAME:=LINKED.x==LINKED.y]
  
  ggplot(tmp,aes(fill=factor(LINKED_SAME,c(TRUE,FALSE),c('yes','no'))))+
    geom_boxplot(aes(x='2019',y=LINKED_SCORE.x)) +
    geom_boxplot(aes(x='2021',y=LINKED_SCORE.y)) +
    labs(x='analysis',y='linkage scores',fill='same linkage classification')+
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",legend.direction = 'horizontal',legend.box = 'horizontal') +
    theme_bw()
  ggsave(filename = '~/compare_linkage_score_couple.pdf',width=6,height=4)
  
  
  tmp[LINKED.x==T & LINKED.y==T,DIRECTED_SAME:=DIRECTED.x==DIRECTED.y]
  tmp[DIRECTED.x==T & DIRECTED.y==T,LINK12_SAME:=LINK12.x==LINK12.y]
  ggplot(tmp[LINKED.x==T & LINKED.y==T],aes(fill=factor(DIRECTED_SAME,c(TRUE,FALSE),c('yes','no'))))+
    geom_boxplot(aes(x='2019',y=DIRECTED_SCORE.x)) +
    geom_boxplot(aes(x='2021',y=DIRECTED_SCORE.y)) +
    labs(x='analysis',y='direction scores',fill='same directed classification')+
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",legend.direction = 'horizontal',legend.box = 'horizontal') +
    theme_bw()
  ggsave(filename = '~/compare_direction_score_couple.pdf',width=6,height=4)
  
  ggplot(tmp[DIRECTED.x==T & DIRECTED.y==T],aes(fill=factor(LINK12_SAME,c(TRUE,FALSE),c('yes','no'))))+
    geom_boxplot(aes(x='2019',y=LINK12_SCORE.x)) +
    geom_boxplot(aes(x='2021',y=LINK12_SCORE.y)) +
    labs(x='analysis',y='1->2 scores',fill='same direction classification')+
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",legend.direction = 'horizontal',legend.box = 'horizontal') +
    theme_bw()
  ggsave(filename = '~/compare_link12_score_couple.pdf',width=6,height=4)
  
  
  dc_couple_1921 <- copy(tmp)
  
  require(Phyloscanner.R.utilities)
  require(phyloscannerR)
  require(tidyverse)
  
  dw19_couple <- dw19[COUPLE==1,c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE')]
  aid[,AID:=as.character(AID)]
  dw19_couple <- merge(dw19_couple, aid,by.x='ID1',by.y='PT_ID',all.x=T)
  dw19_couple <- merge(dw19_couple, aid,by.x='ID2',by.y='PT_ID',all.x=T)
  dw19_couple[,X.x:=NULL]
  dw19_couple[,X.y:=NULL]
  dw19_couple[,ID1:=NULL]
  dw19_couple[,ID2:=NULL]
  setnames( dw19_couple, c('AID.x','AID.y'), c('ID1','ID2'))
  tmp <- dw19_couple[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  set(tmp, NULL, 'BASIC_CLASSIFICATION', tmp[,gsub('xx','desc_',gsub('desc_','anc_',gsub('anc_','xx',BASIC_CLASSIFICATION)))])
  dw19_couple <- unique(rbind(dw19_couple[ID1<ID2],tmp))
  
  dw_couple <- dw[COUPLE==1,c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN')]
  setnames( dw_couple, c('H1','H2'), c('ID1','ID2'))
  tmp <- dw_couple[,min(PTY_RUN), by=c('ID1','ID2')]
  
  
  dw_couple <-  merge(dw_couple,tmp, by.x=c('ID1','ID2','PTY_RUN'),by.y=c('ID1','ID2','V1'))
  tmp <- dw_couple[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  set(tmp, NULL, 'BASIC_CLASSIFICATION', tmp[,gsub('xx','desc_',gsub('desc_','anc_',gsub('anc_','xx',BASIC_CLASSIFICATION)))])
  dw_couple <- unique(rbind(dw_couple[ID1<ID2],tmp))


  # pty.runs[grepl('AID5758|AID7973',RENAME_ID)]
  # tmp1 <- dw_couple[ID1=='AID5758' & ID2=='AID7973',]
  # tmp2 <- dw19_couple[ID1=='AID5758' & ID2=='AID7973',]
  # tmp <- merge(tmp1,tmp2,by='TREE_ID',all=T)
  # tmp[BASIC_CLASSIFICATION.x!=BASIC_CLASSIFICATION.y,c('TREE_ID','BASIC_CLASSIFICATION.x','BASIC_CLASSIFICATION.y')]
  # tmp[BASIC_CLASSIFICATION.x!=BASIC_CLASSIFICATION.y,table(paste0(BASIC_CLASSIFICATION.y,'->',BASIC_CLASSIFICATION.x))]
  # 
  # 
  # tmpdw <- dw[(H1=='AID5758' & H2=='AID7973')|(H2=='AID5758' & H1=='AID7973'),c('TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT')]
  # tmpdw19 <- dw19c[(AID.x=='AID5758' & AID.y=='AID7973')|(AID.y=='AID5758' & AID.x=='AID7973'),c('TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT')]
  # tmpdw <- merge(tmpdw19,tmpdw, by='TREE_ID',all=T)
  # tmpdw[is.na( CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT.x )]
  # tmpdw[is.na( CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT.y )]
  # # tmpdw_maxh <- dw_maxhost[(H1=='AID5758' & H2=='AID7973')|(H2=='AID5758' & H1=='AID7973'),c('TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT')]
  # # tmpdw_minmaxh <- dw_minmaxhost[(H1=='AID5758' & H2=='AID7973')|(H2=='AID5758' & H1=='AID7973'),c('TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT')]
  # # tmpdw <- merge(tmpdw,tmpdw_maxh, by='TREE_ID',all=T)
  # # tmpdw <- merge(tmpdw,tmpdw_minmaxh, by='TREE_ID',all=T)
  # # colnames(tmpdw) <- c("TREE_ID","CAT19",'CAT','CAT_MH','CAT_MMH')
  # # tmpdw[is.na(CAT_MH)]
  # # tmpdw[is.na(CAT_MMH)]
  # # tmpdc <- dc[(H1=='AID5758' & H2=='AID7973')|(H2=='AID5758' & H1=='AID7973'),]
  # #

  
  # aid
  dc_couple_1921 <- merge(dc_couple_1921, subset(aid,select=c('PT_ID','AID')),by.x='ID1',by.y='PT_ID',all.x=T)
  dc_couple_1921 <- merge(dc_couple_1921, subset(aid,select=c('PT_ID','AID')),by.x='ID2',by.y='PT_ID',all.x=T)
  dc_couple_1921[,AID.x:=NULL]
  dc_couple_1921[,AID.y:=NULL]
  
  tmp <- merge(dc_couple_1921, aid,by.x='ID1',by.y='PT_ID',all.x=T)
  tmp <- merge(tmp, aid,by.x='ID2',by.y='PT_ID',all.x=T)
  tmp <- tmp[!(LINKED_SAME==T & DIRECTED_SAME==T & LINK12_SAME==T)]
  
  # change column names 
  inclusion <- "both"# "either"
  setnames(dw19_couple, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  setnames(dw_couple, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  
  library(ggpubr) 
  dir.create('~/compare_2019_2021/')
  control = list(	yintercept_close=0.025,
                  yintercept_dist=1.0,
                  breaks_x=seq(0,1e4,500), 
                  minor_breaks_x=seq(0,1e4,100),
                  breaks_y=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),										
                  limits_y=c(1e-3,0.4),
                  fill.topology=c("ancestral"="deepskyblue1","descendant"="deepskyblue4","intermingled"= "#FDB863",'sibling'="#8073AC","other"="grey80"))
  for (i in 1:nrow(tmp)) {
    cat('processing ',i,' out of ', nrow(tmp), ' pairs \n')
    hosts <- c(tmp$AID.x[i],tmp$AID.y[i])
    tmpx <- dw_couple[(host.1==hosts[1] & host.2==hosts[2]) | (host.1==hosts[2] & host.2==hosts[1]),'tree.id']
    tmpx <- rbind(tmpx, dw19_couple[(host.1==hosts[1] & host.2==hosts[2]) | (host.1==hosts[2] & host.2==hosts[1]),'tree.id'])
    tmpx <- as.integer(gsub('([0-9]+)_to_([0-9]+)','\\1',tmpx$tree.id))
    control$limits_x <- range(tmpx) + c(-25,25) # bin width
    tmpp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw19_couple), inclusion = "both",control=control)
    g1 <- tmpp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    tmpp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw_couple), inclusion = "both",control=control)
    g2 <- tmpp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    library(ggpubr)
    arrange <- ggarrange(g1 + ggtitle('2019 analysis \n '),
                         g2 + ggtitle('2021 analysis, \n without minreadhost and maxreadhost'), 
                         ncol = 1, nrow = 2, common.legend = T, legend='bottom')
    ggsave(paste0('~/compare_2019_2021/',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 8)
  }
  
  
  
  # v loop
  unique(subset(dw19_couple,select=c('host.1','host.2')))
  dw19_couple[,start:=as.integer(gsub('([0-9]+)_to_([0-9]+)','\\1',tree.id))]
  tmp <- dw19_couple[start>=6615-250 & start<=7636 & start!=6825 & start!=6850]
  tmp <- unique(subset(tmp,select=c('host.1','host.2')))
  tmp2 <- unique(subset(dw_couple,select=c('host.1','host.2')))
  tmp2[,host.x:=paste0(host.1,'-',host.2)]
  tmp2[,host.y:=paste0(host.2,'-',host.1)]
  tmp <- tmp[paste0(host.1,'-', host.2)%in% tmp2$host.x] 
  control = list(	yintercept_close=0.025,
                  yintercept_dist=1.0,
                  breaks_x=seq(0,1e4,500), 
                  minor_breaks_x=seq(0,1e4,100),
                  breaks_y=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),										
                  limits_y=c(1e-3,0.4),
                  fill.topology=c("ancestral"="deepskyblue1","descendant"="deepskyblue4","intermingled"= "#FDB863",'sibling'="#8073AC","other"="grey80"))
  
  g1 <- list()
  g2 <- list()
  for (i in 1:nrow(tmp)) {
    cat('processing ',i,' out of ', nrow(tmp), ' pairs \n')
    hosts <- c(tmp$host.1[i],tmp$host.2[i])
    tmpx <- dw_couple[(host.1==hosts[1] & host.2==hosts[2]) | (host.1==hosts[2] & host.2==hosts[1]),'tree.id']
    tmpx <- rbind(tmpx, dw19_couple[(host.1==hosts[1] & host.2==hosts[2]) | (host.1==hosts[2] & host.2==hosts[1]),'tree.id'])
    tmpx <- as.integer(gsub('([0-9]+)_to_([0-9]+)','\\1',tmpx$tree.id))
    control$limits_x <- range(tmpx) + c(-25,25)
    tmpp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw19_couple), inclusion = "both",control=control)
    g1[[i]] <- tmpp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+ 
      ggtitle('2019 analysis \n ')
    tmpp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw_couple), inclusion = "both",control=control)
    g2[[i]] <- tmpp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
      ggtitle('2021 analysis, \n without minreadhost and maxreadhost')
  }
  glist <- c(g1,g2)
  library(ggpubr)
  arrange <- ggarrange(plotlist = glist, 
                       ncol = nrow(tmp), nrow = 2, common.legend = T, legend='bottom')
  ggsave(paste0('~/compare_2019_2021/vloop_couples.pdf'),arrange,width = 5*nrow(tmp), height = 8,limitsize = F)
  
  
  
  # compare 2019 and 2021 all
  
  library(ggpubr) 
  dir.create('~/compare_2019_2021_all/')
  # load 2021 maxreadhost
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read/' #30 min 100 max
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  load(infile.networks)
  dw_maxhost <- merge(dw, aid, by.x='H1',by.y='AID',all.x=T)
  dw_maxhost <- merge(dw_maxhost, aid, by.x='H2',by.y='AID',all.x=T)
  setnames(dw_maxhost, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  dw_maxhost <- merge(dw_maxhost,couple,by=c('ID1','ID2'),all.x=T)
  
  dw_maxhost_couple <- dw_maxhost[COUPLE==1,c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN')]
  setnames( dw_maxhost_couple, c('H1','H2'), c('ID1','ID2'))
  tmp <- dw_maxhost_couple[,min(PTY_RUN), by=c('ID1','ID2')]
  dw_maxhost_couple <-  merge(dw_maxhost_couple,tmp, by.x=c('ID1','ID2','PTY_RUN'),by.y=c('ID1','ID2','V1'))
  tmp <- dw_maxhost_couple[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  set(tmp, NULL, 'BASIC_CLASSIFICATION', tmp[,gsub('xx','desc_',gsub('desc_','anc_',gsub('anc_','xx',BASIC_CLASSIFICATION)))])
  dw_maxhost_couple <- unique(rbind(dw_maxhost_couple[ID1<ID2],tmp))
  
  setnames(dw_maxhost_couple, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  
  # load 2021 minmaxreadhost
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read/' #30 min 100 max
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  load(infile.networks)
  # dw[H1=='AID5758' & H2=='AID7973']
  dw_minmaxhost <- merge(dw, aid, by.x='H1',by.y='AID',all.x=T)
  dw_minmaxhost <- merge(dw_minmaxhost, aid, by.x='H2',by.y='AID',all.x=T)
  setnames(dw_minmaxhost, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  dw_minmaxhost <- merge(dw_minmaxhost,couple,by=c('ID1','ID2'),all.x=T)
  
  dw_minmaxhost_couple <- dw_minmaxhost[COUPLE==1,c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN')]
  setnames( dw_minmaxhost_couple, c('H1','H2'), c('ID1','ID2'))
  tmp <- dw_minmaxhost_couple[,min(PTY_RUN), by=c('ID1','ID2')]
  dw_minmaxhost_couple <-  merge(dw_minmaxhost_couple,tmp, by.x=c('ID1','ID2','PTY_RUN'),by.y=c('ID1','ID2','V1'))
  tmp <- dw_minmaxhost_couple[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  set(tmp, NULL, 'BASIC_CLASSIFICATION', tmp[,gsub('xx','desc_',gsub('desc_','anc_',gsub('anc_','xx',BASIC_CLASSIFICATION)))])
  dw_minmaxhost_couple <- unique(rbind(dw_minmaxhost_couple[ID1<ID2],tmp))
  
  setnames(dw_minmaxhost_couple, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  
  # pair selection
  tmp <- merge(dc_couple_1921, aid,by.x='ID1',by.y='PT_ID',all.x=T)
  tmp <- merge(tmp, aid,by.x='ID2',by.y='PT_ID',all.x=T)
  tmp <- tmp[!(LINKED_SAME==T & DIRECTED_SAME==T & LINK12_SAME==T)]
  
  
  
  tmp2 <- unique(subset(dw_maxhost_couple,select=c('host.1','host.2')))
  tmp2[,host.x:=paste0(host.1,'-',host.2)]
  tmp2[,host.y:=paste0(host.2,'-',host.1)]
  tmp <- tmp[paste0(AID.x,'-',AID.y) %in% tmp2$host.x | paste0(AID.x,'-',AID.y) %in% tmp2$host.y]
  tmp2 <- unique(subset(dw_minmaxhost_couple,select=c('host.1','host.2')))
  tmp2[,host.x:=paste0(host.1,'-',host.2)]
  tmp2[,host.y:=paste0(host.2,'-',host.1)]
  tmp <- tmp[paste0(AID.x,'-',AID.y) %in% tmp2$host.x | paste0(AID.x,'-',AID.y) %in% tmp2$host.y]
  
  control = list(	yintercept_close=0.025,
                  yintercept_dist=1.0,
                  breaks_x=seq(0,1e4,500), 
                  minor_breaks_x=seq(0,1e4,100),
                  breaks_y=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),										
                  limits_y=c(1e-3,0.4),
                  fill.topology=c("ancestral"="deepskyblue1","descendant"="deepskyblue4","intermingled"= "#FDB863",'sibling'="#8073AC","other"="grey80"))
  for (i in 1:nrow(tmp)) {
    cat('processing ',i,' out of ', nrow(tmp), ' pairs \n')
    hosts <- c(tmp$AID.x[i],tmp$AID.y[i])
    tmpx <- dw_couple[(host.1==hosts[1] & host.2==hosts[2]) | (host.1==hosts[2] & host.2==hosts[1]),'tree.id']
    tmpx <- rbind(tmpx, dw19_couple[(host.1==hosts[1] & host.2==hosts[2]) | (host.1==hosts[2] & host.2==hosts[1]),'tree.id'])
    tmpx <- as.integer(gsub('([0-9]+)_to_([0-9]+)','\\1',tmpx$tree.id))
    control$limits_x <- range(tmpx) + c(-25,25) # bin width
    tmpp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw19_couple), inclusion = "both",control=control)
    g1 <- tmpp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    tmpp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw_couple), inclusion = "both",control=control)
    g2 <- tmpp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    tmpp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw_maxhost_couple), inclusion = "both",control=control)
    g3 <- tmpp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    tmpp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw_minmaxhost_couple), inclusion = "both",control=control)
    g4 <- tmpp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    # library(ggpubr)
    arrange <- ggarrange(g1 + ggtitle('2019 analysis \n '),
                         g2 + ggtitle('2021 analysis, \n without minreadhost and maxreadhost'),
                         g3 + ggtitle('2021 analysis, \n without minreadhost and with maxreadhost'),
                         g4 + ggtitle('2021 analysis, \n with minreadhost and maxreadhost'),
                         ncol = 1, nrow = 4, common.legend = T, legend='bottom')
    ggsave(paste0('~/compare_2019_2021_all/',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 12)
  }
  
  # couple plots
  library(ggplot2)
  load('~/dcdwina_minreadnull.rda')
  dwina[,host.1:=gsub('CNTRL-','',host.1)]
  dwina[,host.2:=gsub('CNTRL-','',host.2)]
  dwinac <- merge(dwina , aid, by.x='host.1',by.y='AID',all.x=T)
  dwinac <- merge(dwinac , aid, by.x='host.2',by.y='AID',all.x=T)
  dwinac <- merge(dwinac,couple,by.x=c('PT_ID.x','PT_ID.y'),by.y=c('ID1','ID2'),all.x=T)
  tmp <- unique(subset(dwinac[COUPLE==1],select=c('host.1','host.2','PTY_RUN','tree.id', 'patristic.distance')))
  
  ggplot(tmp, aes(patristic.distance))+
    geom_histogram(binwidth = 0.002)+theme_bw()+labs(x='patristic distance', y='count')
  ggsave(filename = '~/pdistance_couple.pdf',width = 6, height = 4)
  
  ggplot(tmp, aes(patristic.distance))+
     geom_histogram(binwidth = 0.002)+theme_bw()+labs(x='patristic distance', y='count') +coord_cartesian(ylim = c(0,1e3),xlim=c(0,0.1))
  ggsave(filename = '~/pdistance_couple_zoomin.pdf',width = 6, height = 4)
  
  # pd <- tmp[patristic.distance>0.05 & patristic.distance<0.15,]$patristic.distance
  # library(MASS)
  # fit <- fitdistr(pd, "normal")
  # qnorm(0.05,fit$estimate[1],fit$estimate[2])
  
  # make couple plots
  dw_couple <- dwinac[COUPLE==1]
  tmp <- dw_couple[,list(M=median(patristic.distance),
                         CL=quantile(patristic.distance,probs = 0.025),
                         CU=quantile(patristic.distance,probs = 0.975)),
                   by=c('host.2','host.1')]
  setkey(tmp,M)
  tmp[,ID:=seq_len(nrow(tmp))]
  ggplot(tmp,aes(ID,M))+
    geom_point()+
    theme_bw()+
    geom_errorbar(aes(ymin=CL,ymax=CU))+
    labs(x='couples',y='patristic distance')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())s
  ggsave('~/patdist_vs_couple.pdf',width = 10, height = 4)
  
  tmpid <- subset(tmp,select = c('host.1','host.2','ID'))
  dcac <- merge(dca,subset(aid,select=c('PT_ID','AID')),by.x='host.1',by.y='AID',all=T)
  dcac <- merge(dcac,subset(aid,select=c('PT_ID','AID')),by.x='host.2',by.y='AID',all=T)
  dcac <- merge(dcac,couple,by.x=c('PT_ID.x','PT_ID.y'),by.y=c('ID1','ID2'),all.x=T)
  dc_couple <- dcac[COUPLE==1 & categorisation == 'close.and.adjacent.and.ancestry.cat']
  tmp <- dc_couple[,list(n.eff=max(n.eff)),by=c('host.1','host.2')]
  tmp <- merge(tmp,dc_couple, by=c('host.1','host.2','n.eff'))
  setkey(tmp, host.1, host.2, PTY_RUN, type)
  tmp <- tmp[,head(.SD, 1),by=c('host.1','host.2','type')]
  tmp <- merge(tmp, tmpid, by=c('host.1','host.2'))
  # tmp2 = tmp[,length(n.eff),by=c('host.1','host.2')]
  # tmp2 = tmp[,as.integer(sum(score)),by=c('host.1','host.2')]
  ggplot(tmp,aes(ID,k.eff,fill=type))+
    theme_bw()+
    geom_bar(position="stack", stat="identity")+
    labs(x='couples',y='number of windows \n supporting subgraph topology')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  ggsave('~/topo_num_vs_couple.pdf',width = 10, height = 4)
  
  ggplot(tmp,aes(ID,score,fill=type))+
    theme_bw()+
    geom_bar(position="stack", stat="identity")+
    labs(x='couples',y='proportion of windows \n supporting subgraph topology')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  ggsave('~/topo_prop_vs_couple.pdf',width = 10, height = 4)
  
  tmp <- dc_couple[,list(n.eff=max(n.eff)),by=c('host.1','host.2')]
  tmp <- merge(tmp,dc_couple, by=c('host.1','host.2','n.eff'))
  setkey(tmp, host.1, host.2, PTY_RUN, type)
  tmp <- tmp[,head(.SD, 1),by=c('host.1','host.2','type')]
  tmp <- merge(tmp, tmpid, by=c('host.1','host.2'))
  tmp <- dcast(tmp, host.1+host.2 + ID~type, value.var = 'score')
  tmp <- data.table(tmp)
  tmp[,score_L:=`12`+ `21`+ complex.or.no.ancestry]
  tmp[,score_D:=(`12`+ `21`)/score_L]
  ggplot(tmp,aes(ID,score_L))+
    geom_point()+
    theme_bw()+
    scale_y_continuous(labels = scales::percent,limits = c(0,1))+
    # geom_errorbar(aes(ymin=CL,ymax=CU))+
    labs(x='couples',y='linkage score')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  ggsave('~/lscore_vs_couple.pdf',width = 10, height = 4)
  
  ggplot(tmp,aes(ID,score_D))+
    geom_point()+
    theme_bw()+
    # geom_errorbar(aes(ymin=CL,ymax=CU))+
    labs(x='couples',y='linkage score')+
    scale_y_continuous(labels = scales::percent,limits = c(0,1))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  ggsave('~/dscore_vs_couple.pdf',width = 10, height = 4)
  
  
}

  
  
test_extra_windows <- function(){
  x <- mget(load("/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_Matthew_lastrun_191011/PANGEA_54wins_Cluster_9_newIDs_workspace.rda", envir=(NE. <- new.env())), envir=NE.)
  y <- mget(load("/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/ptyr1_workspace.rda", envir=(NE. <- new.env())), envir=NE.)
  x$args
  tmp <- lapply(1:length(x$args), function(i)x$args[[i]]==y$args[[i]])
  tmp <- lapply(tmp, function(x)ifelse(length(x)==0,NA,x))
  id <- sapply(tmp, function(x)x==TRUE)
  idd <- which(id!=TRUE|is.na(id))
  for (i in idd) {
    cat(i,' th arguments - ', names(x$args)[i],'\n',
        '2019 values: ',x$args[[i]],'\n',
        '2021 values: ',y$args[[i]],'\n')
  }
  
  require(Phyloscanner.R.utilities)
  require(phyloscannerR)
  require(tidyverse)
  library(data.table)
  
  indir1	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_test_extra_window/'
  infile.networks1 <-  file.path(indir1,'Rakai_phscnetworks.rda')
  indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_test_extra_window2/'
  infile.networks2 <-  file.path(indir2,'Rakai_phscnetworks.rda')
  indir3	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_test_extra_window3/'
  infile.networks3 <-  file.path(indir3,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  
  load( infile.networks1)
  dw1 <- copy(dw)
  dw1 <- dw1[,c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN')]
  setnames( dw1, c('H1','H2'), c('ID1','ID2'))
  tmp <- dw1[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  set(tmp, NULL, 'BASIC_CLASSIFICATION', tmp[,gsub('xx','desc_',gsub('desc_','anc_',gsub('anc_','xx',BASIC_CLASSIFICATION)))])
  dw1 <- unique(rbind(dw1[ID1<ID2],tmp))
  load( infile.networks2)
  dw2 <- copy(dw)
  dw2 <- dw2[,c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN')]
  setnames( dw2, c('H1','H2'), c('ID1','ID2'))
  tmp <- dw2[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  set(tmp, NULL, 'BASIC_CLASSIFICATION', tmp[,gsub('xx','desc_',gsub('desc_','anc_',gsub('anc_','xx',BASIC_CLASSIFICATION)))])
  dw2 <- unique(rbind(dw2[ID1<ID2],tmp))
  load( infile.networks3)
  dw3 <- copy(dw)
  dw3 <- dw3[,c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN')]
  setnames( dw3, c('H1','H2'), c('ID1','ID2'))
  tmp <- dw3[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  set(tmp, NULL, 'BASIC_CLASSIFICATION', tmp[,gsub('xx','desc_',gsub('desc_','anc_',gsub('anc_','xx',BASIC_CLASSIFICATION)))])
  dw3 <- unique(rbind(dw3[ID1<ID2],tmp))
  aid <- data.table(read.csv(infile.ind.anonymised))
  


  # tmp <- subset(aid,select=c('PT_ID','AID'))
  # dc <- merge(dc, tmp, by.x='H1',by.y='AID',all.x=T)
  # dc <- merge(dc, tmp, by.x='H2',by.y='AID',all.x=T)
  # setnames(dc, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  # dw <- merge(dw, tmp, by.x='H1',by.y='AID',all.x=T)
  # dw <- merge(dw, tmp, by.x='H2',by.y='AID',all.x=T)
  # setnames(dw, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))

  inclusion <- "both"# "either"
  setnames(dw1, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  setnames(dw2, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  setnames(dw3, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  
  control = list(	yintercept_close=0.025,
                  yintercept_dist=1.0,
                  breaks_x=seq(0,1e4,500), 
                  minor_breaks_x=seq(0,1e4,100),
                  breaks_y=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),										
                  limits_y=c(1e-3,0.4),
                  fill.topology=c("ancestral"="deepskyblue1","descendant"="deepskyblue4","intermingled"= "#FDB863",'sibling'="#8073AC","other"="grey80"))

    hosts <- c('AID7973','AID5758')
    tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw1), inclusion = "both",control=control)
    g1<- tmp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw2), inclusion = "both",control=control)
    g2<- tmp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw3), inclusion = "both",control=control)
    g3<- tmp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
 
    library(ggpubr)
    arrange <- ggarrange(g1,g2,g3, 
                         ncol = 1, nrow = 3, common.legend = T, legend='bottom')
    ggsave(paste0('~/test_',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 12)
    
    
    #5758_7973
    library(data.table)
    indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/'
    # indir <- '~/hpc210506/'
    infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
    # infile.workspace <- file.path(indir,'ptyr72_workspace.rda')
    infile.workspace <- file.path(indir,'ptyr389_workspace.rda')
    infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
    infile.ind.anonymised <- '~/hpc210506/important_anonymisation_keys_210119.csv'
    
    aid <- data.table(read.csv(infile.ind.anonymised))
    load(infile.workspace)
    tl <- list()
    for (i in 1:length(phyloscanner.trees)) {
      tl[[i]] <- data.table(TIP=phyloscanner.trees[[i]]$tree$tip.label)
    }
    names(tl) <- names(phyloscanner.trees)
    tl <- rbindlist(tl,idcol = T)
    setnames(tl,'.id', 'TREE_ID')
    dt <- tl[grep('AID5758|AID7973',TIP)]
    # dt <- tl[grep('AID2168|AID5977|AID4448',TIP)]
    dt[,AID:=gsub('^(AID[0-9]+)-fq([0-9]+)_read_([0-9]+)_count_([0-9]+)','\\1',TIP)]
    dt[,FQ:=gsub('^(AID[0-9]+)-fq([0-9]+)_read_([0-9]+)_count_([0-9]+)','\\2',TIP)]
    dt[,READ:=gsub('^(AID[0-9]+)-fq([0-9]+)_read_([0-9]+)_count_([0-9]+)','\\3',TIP)]
    dt[,COUNT:=gsub('^(AID[0-9]+)-fq([0-9]+)_read_([0-9]+)_count_([0-9]+)','\\4',TIP)]
    dt[,READ:=as.integer(gsub('_X_CONTAMINANT','',READ))]
    dt[,COUNT:=as.integer(gsub('_X_CONTAMINANT','',COUNT))]
    dt[,AID:=gsub('_X_CONTAMINANT','',AID)]
    # dt[,list(TOTAL_COUNT=sum(COUNT)),by='AID']
    
    dt2 <- unique(dt[COUNT>=30,c('TREE_ID','AID')])
    dt <- dcast(dt2,TREE_ID~AID)
    dt <- dt[!is.na(AID5758) & !is.na(AID7973),'TREE_ID']
       
    load( infile.networks )
    dw <- dw[(H1=='AID5758' & H2=='AID7973')|(H2=='AID5758' & H1=='AID7973')]
    dw <- subset(dw,select=c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN'))
    setnames(dw,c('H1','H2'),c('ID1','ID2'))
    tmp <- dw[!(ID1<ID2)]
    setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
    set(tmp, NULL, 'BASIC_CLASSIFICATION', tmp[,gsub('xx','desc_',gsub('desc_','anc_',gsub('anc_','xx',BASIC_CLASSIFICATION)))])
    dw <- unique(rbind(dw[ID1<ID2],tmp))

    dw <- merge(dw,dt,by='TREE_ID')
    inclusion <- "both"# "either"
    setnames(dw, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
             c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
    dw <- unique(dw)
    # 
    control = list(	yintercept_close=0.025,
                    yintercept_dist=1.0,
                    breaks_x=seq(0,1e4,500), 
                    minor_breaks_x=seq(0,1e4,100),
                    breaks_y=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),										
                    limits_y=c(1e-3,0.4),
                    fill.topology=c("ancestral"="deepskyblue1","descendant"="deepskyblue4","intermingled"= "#FDB863",'sibling'="#8073AC","other"="grey80"))
    
    hosts <- c('AID7973','AID5758')
    tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw), inclusion = "both",control=control)
    g<- tmp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    ggsave(paste0('~/post_minread_30_',hosts[1],'_',hosts[2],'.pdf'),g,width = 8, height = 4)
    
    # 5977-2168
    library(data.table)
    indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/'
    infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
    infile.workspace <- file.path(indir,'ptyr389_workspace.rda')
    infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
    aid <- data.table(read.csv(infile.ind.anonymised))
    load(infile.workspace)
    tl <- list()
    for (i in 1:length(phyloscanner.trees)) {
      tl[[i]] <- data.table(TIP=phyloscanner.trees[[i]]$tree$tip.label)
    }
    names(tl) <- names(phyloscanner.trees)
    tl <- rbindlist(tl,idcol = T)
    setnames(tl,'.id', 'TREE_ID')
    dt <- tl[grep('AID2168|AID5977|AID4448',TIP)]
    dt[,AID:=gsub('^(AID[0-9]+)-fq([0-9]+)_read_([0-9]+)_count_([0-9]+)','\\1',TIP)]
    dt[,FQ:=gsub('^(AID[0-9]+)-fq([0-9]+)_read_([0-9]+)_count_([0-9]+)','\\2',TIP)]
    dt[,READ:=gsub('^(AID[0-9]+)-fq([0-9]+)_read_([0-9]+)_count_([0-9]+)','\\3',TIP)]
    dt[,COUNT:=gsub('^(AID[0-9]+)-fq([0-9]+)_read_([0-9]+)_count_([0-9]+)','\\4',TIP)]
    # dt[AID=='AID2168_X_CONTAMINANT',AID:='AID2168']
    # dt[grep('_X_CONTAMINANT',READ)]
    dt[,READ:=as.integer(gsub('_X_CONTAMINANT','',READ))]
    dt[,COUNT:=as.integer(gsub('_X_CONTAMINANT','',COUNT))]
    dt[,AID:=gsub('_X_CONTAMINANT','',AID)]
    dt_reads <- dt[,list(TORAL_COUNT=sum(as.integer(COUNT))),by=c('AID','TREE_ID')]
    dt_reads[,START:=as.integer(gsub('([0-9]+)_to_([0-9]+)','\\1',TREE_ID))]
    library(ggplot2)
    ggplot(dt_reads,
           aes(START, TORAL_COUNT,fill=AID))  + 
      geom_bar(position="dodge", stat="identity")+
      labs(x='\n genomic windows', y='total counts', fill='ID')+theme_bw() +
      guides(fill = guide_legend(nrow = 1)) + theme(legend.position = "bottom",legend.direction = 'horizontal')
   # ggsave('~/reads_AID5977_AID2168.pdf',width = 10, height = 4) 
    ggsave('~/reads_AID5977_AID2168_AID4448.pdf',width = 10, height = 4) 
}

not.balance.check <- function()
{
  library(data.table)
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/'
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  indir2	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read_100_max_read/'
  infile.networks2 <-  file.path(indir2,'Rakai_phscnetworks.rda')
  infile.workspace <- file.path(indir,paste0('ptyr',1:409,'_workspace.rda'))
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  # infile.ind.anonymised <- '~/hpc210506/important_anonymisation_keys_210119.csv'
  infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
  
  aid <- data.table(read.csv(infile.ind.anonymised))
  tll <- list()
  for(k in 1:409){
    cat('run ',k, '\n')
    try(load(infile.workspace[k]))
    tl <- list()
    for (i in 1:length(phyloscanner.trees)) {
      tl[[i]] <- data.table(TIP=phyloscanner.trees[[i]]$tree$tip.label)
    }
    names(tl) <- names(phyloscanner.trees)
    tl <- rbindlist(tl,idcol = T) 
    setnames(tl,'.id', 'TREE_ID')
    tll[[k]] <- tl
  }
  # error from 
  tll <- rbindlist(tll,idcol = T) 
  setnames(tll,'.id', 'RUN')
  # tll <- tll[! RUN%in% c(99,173,180,182,190)] 
  # 0byte trees
  # save(tll,file = '~/alltips.rda')
  # load('~/alltips.rda')
  load(infile.couple)
  couple <- unique(subset(data.table(coupdat), select=c('male.RCCS_studyid', 'female.RCCS_studyid')))
  couple[, COUPLE:=1]
  setnames(couple,c('male.RCCS_studyid','female.RCCS_studyid'),c('PT_ID1','PT_ID2'))
  couple[,PT_ID1:=paste0('RK-',PT_ID1)]
  couple[,PT_ID2:=paste0('RK-',PT_ID2)]
  couple <- merge(couple,subset(aid,select=c('PT_ID', 'AID')),by.x='PT_ID1', by.y='PT_ID', all.x=T)
  couple <- merge(couple,subset(aid,select=c('PT_ID', 'AID')),by.x='PT_ID2', by.y='PT_ID', all.x=T)
  set(couple,NULL,c('PT_ID1','PT_ID2'),NULL)
  setnames(couple,c('AID.x','AID.y'),c('H1','H2'))
  couple[,H1:=as.character(H1)]
  couple[,H2:=as.character(H2)]
  couple <- couple[!is.na(H1) & !is.na(H2)]
  # tmp <- copy(couple)
  # setnames(tmp,c('H1','H2'),c('H2','H1'))
  # couple <- rbind(couple,tmp)
  
  dt <- tll[grep('AID',TIP)]
  dt[,START:=as.integer(gsub('([0-9]+)_to_([0-9]+)','\\1',TREE_ID))]
  dt <- dt[START<=2292-250]
  dt[grep('_X_CONTAMINANT',TIP),CONTAMINANT:=T]
  dt[,TIP:=gsub('_X_CONTAMINANT','',TIP)]
  dt[grep('CNTRL-',TIP),CNTRL:=T]
  dt[,TIP:=gsub('CNTRL-','',TIP)]
  
  dt[,AID:=gsub('^(AID[0-9]+)-fq([0-9]+)_read_([0-9]+)_count_([0-9]+)','\\1',TIP)]
  dt[,FQ:=gsub('^(AID[0-9]+)-fq([0-9]+)_read_([0-9]+)_count_([0-9]+)','\\2',TIP)]
  dt[,READ:=gsub('^(AID[0-9]+)-fq([0-9]+)_read_([0-9]+)_count_([0-9]+)','\\3',TIP)]
  dt[,COUNT:=gsub('^(AID[0-9]+)-fq([0-9]+)_read_([0-9]+)_count_([0-9]+)','\\4',TIP)]
  # dt[,READ:=as.integer(gsub('_X_CONTAMINANT','',READ))]
  # dt[,COUNT:=as.integer(gsub('_X_CONTAMINANT','',COUNT))]
  # dt[,AID:=gsub('_X_CONTAMINANT','',AID)]
  save(dt,file = '~/readcountgag.rda')
  
  load('~/readcountgag.rda')
  dt <- dt[is.na(CONTAMINANT)]
  dt <- dt[AID%in% unique(c(couple$H1,couple$H2))]
  dta <- dt[,list(TOTAL_COUNT=sum(as.integer(COUNT))),by=c('AID','TREE_ID')]
  dtac <- merge(couple,dta,by.x='H1',by.y='AID',all.x=T)
  dtac <- merge(dtac,dta,by.x=c('H2','TREE_ID'),by.y=c('AID','TREE_ID'),all.x=T)
  dtac <- dtac[!is.na(TOTAL_COUNT.x) & !is.na(TOTAL_COUNT.y)]
  dtac[,TOTAL_COUNT.ad:=abs(TOTAL_COUNT.x-TOTAL_COUNT.y)]
  dtac[,TOTAL_COUNT.xdy:=TOTAL_COUNT.x/TOTAL_COUNT.y]
  dtac[,TOTAL_COUNT.ydx:=TOTAL_COUNT.y/TOTAL_COUNT.x]
  dtac[,TOTAL_COUNT.r:=pmax(TOTAL_COUNT.xdy,TOTAL_COUNT.ydx)]
  tmp <- dtac[!(H1<H2)]
  setnames(tmp,c('H1','H2'),c('H2','H1'))
  dtac <- unique(rbind(dtac[H1<H2],tmp))
  
  # dtac[H1=='AID5758' & H2=='AID7973']
  
  # hist(dtac$TOTAL_COUNT.ad)
  # hist(dtac$TOTAL_COUNT.r)
  tmp <- dtac[,list(NTREE=length(TREE_ID)), by=c('H1','H2')]
  dtac <- merge(dtac,tmp,by=c('H1','H2'))
  
  dtac_ad <- dtac[TOTAL_COUNT.ad>=100]
  tmp <- dtac_ad[,list(NTREE_AD=length(TREE_ID)), by=c('H1','H2')]
  dtac_ad <- merge(dtac_ad,tmp,by=c('H1','H2'))
  dtac_r <- dtac[TOTAL_COUNT.r>=1.5]
  tmp <- dtac_r[,list(NTREE_R=length(TREE_ID)), by=c('H1','H2')]
  dtac_r <- merge(dtac_r,tmp,by=c('H1','H2'))
  # hist(tmp$NTREE)
  
  
  load( infile.networks )
  dw1 <- subset(dw,select=c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN'))
  setnames(dw1,c('H1','H2'),c('ID1','ID2'))
  tmp <- dw1[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  set(tmp, NULL, 'BASIC_CLASSIFICATION', tmp[,gsub('xx','desc_',gsub('desc_','anc_',gsub('anc_','xx',BASIC_CLASSIFICATION)))])
  dw1 <- unique(rbind(dw1[ID1<ID2],tmp))
  tmp <- dw1[,min(PTY_RUN), by=c('ID1','ID2')]
  dw1 <-  merge(dw1,tmp, by.x=c('ID1','ID2','PTY_RUN'),by.y=c('ID1','ID2','V1'))
  
  inclusion <- "both"# "either"
  setnames(dw1, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  
  dw1 <- unique(dw1)
  # 
  load( infile.networks2 )
  dw2 <- subset(dw,select=c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN'))
  setnames(dw2,c('H1','H2'),c('ID1','ID2'))
  tmp <- dw2[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  set(tmp, NULL, 'BASIC_CLASSIFICATION', tmp[,gsub('xx','desc_',gsub('desc_','anc_',gsub('anc_','xx',BASIC_CLASSIFICATION)))])
  dw2 <- unique(rbind(dw2[ID1<ID2],tmp))
  tmp <- dw2[,min(PTY_RUN), by=c('ID1','ID2')]
  dw2 <-  merge(dw2,tmp, by.x=c('ID1','ID2','PTY_RUN'),by.y=c('ID1','ID2','V1'))
  
  
  inclusion <- "both"# "either"
  setnames(dw2, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  dw2 <- unique(dw2)
  # 
  
  tmp_dtac_ad <- unique(subset(dtac_ad,select=c('H1','H2','TREE_ID')))
  tmp_dtac_r <- unique(subset(dtac_r,select=c('H1','H2','TREE_ID')))
  tmp_dw1 <- merge(dw1, tmp_dtac_ad, by.x=c('host.1','host.2','tree.id'),by.y=c('H1','H2','TREE_ID'))
  tmp_dw2 <- merge(dw2, tmp_dtac_ad, by.x=c('host.1','host.2','tree.id'),by.y=c('H1','H2','TREE_ID'))
  tmp_dw <- merge(tmp_dw1,tmp_dw2,by=c('host.1','host.2','tree.id'))
  table(tmp_dw1$basic.classification)
  table(tmp_dw2$basic.classification)
  tmp <- tmp_dw[basic.classification.x!=basic.classification.y,]
  unique(subset(tmp,select=c('host.1','host.2')))
  table(tmp$basic.classification.x)
  table(tmp$basic.classification.y)
  tmp <- tmp[,length(tree.id),by=c('host.1','host.2')]
  setkey(tmp,V1)
  tmp <- tail(tmp,10)
  # tmp <- unique(subset(dw2,select=c('host.1','host.2')))
  # tmp[,host:=paste0(host.1,'-',host.2)]
  # dtac_ad <- dtac_ad[paste0(H1,'-',H2) %in% tmp$host]
  # dtac_r <- dtac_r[paste0(H1,'-',H2) %in% tmp$host]
  # dtac_ads <- unique(subset(dtac_ad,select=c('H1','H2','NTREE','NTREE_AD')))
  # dtac_rs <- unique(subset(dtac_r,select=c('H1','H2','NTREE','NTREE_R')))
  # dtac_ads[,P:=NTREE_AD/NTREE]
  # dtac_rs[,P:=NTREE_R/NTREE]
  # setkey(dtac_ads,P,NTREE)
  # setkey(dtac_rs,P,NTREE)
  # hist(dtac_ads$P)
  # hist(dtac_rs$P)
  # dsac_s <- merge(dtac_ads,dtac_rs,by=c('H1','H2'),all=T)
  # setkey(dsac_s,P.x,P.y,NTREE.x,NTREE.y)
  # tmp <- tail(dsac_s,10)
  # dsac_s[H1=='AID2168' & H2=='AID5977']
  
  library(tidyverse)
  control = list(	yintercept_close=0.025,
                  yintercept_dist=1.0,
                  breaks_x=seq(0,1e4,500),
                  minor_breaks_x=seq(0,1e4,100),
                  breaks_y=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),
                  limits_y=c(1e-3,0.4),
                  fill.topology=c("ancestral"="deepskyblue1","descendant"="deepskyblue4","intermingled"= "#FDB863",'sibling'="#8073AC","other"="grey80"))

  # dir.create('~/compare_2021_imbalance/')
  dir.create('~/compare_2021_imbalance_v2/')
  for (i in 1:nrow(tmp)) {
    cat('processing ',i,' out of ', nrow(tmp), ' pairs \n')
    hosts <- c(tmp$host.1[i],tmp$host.2[i])
    tmpx <- dw1[(host.1==hosts[1] & host.2==hosts[2]) | (host.1==hosts[2] & host.2==hosts[1]),'tree.id']
    tmpx <- rbind(tmpx, dw2[(host.1==hosts[1] & host.2==hosts[2]) | (host.1==hosts[2] & host.2==hosts[1]),'tree.id'])
    tmpx <- as.integer(gsub('([0-9]+)_to_([0-9]+)','\\1',tmpx$tree.id))
    control$limits_x <- range(tmpx) + c(-25,25) # bin width
    tmpp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw1), inclusion = "both",control=control)
    g1<- tmpp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    tmpp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw2), inclusion = "both",control=control)
    g2<- tmpp$graph +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    arrange <- ggarrange(g1 + ggtitle('2021 analysis, \n without minreadhost and maxreadhost'),
                         g2 + ggtitle('2021 analysis, \n without minreadhost and with maxreadhost'),
                         ncol = 1, nrow = 2, common.legend = T, legend='bottom')
    ggsave(paste0('~/compare_2021_imbalance_v2/',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 8)
  }
}

debug_minreadperhost <- function()
{# ifelse(post.hoc.min.counts, 1, min.reads.per.host),
  # ifelse(post.hoc.min.counts, 1, min.tips.per.host),
  
  args <- list()
  args$tree="/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI//210325_phsc_output//ptyr72_trees"
  args$outputString = "ptyr72"
  args$outgroupName <- "REF_CON_H"
  args$multifurcationThreshold <- 1e-5
  args$userBlacklist <- NULL
  args$alignment <- NULL
  args$outputDir <- NULL
  args$verbose <- 1
  args$noProgressBars <- T
  args$tipRegex <-"^(.*)-fq[0-9]+_read_([0-9]+)_count_([0-9]+)"
  # args$fileNameRegex <- "^(?:.*\D)?([0-9]+)_to_([0-9]+).*$"
  args$fileNameRegex <- "^(?:.*\\D)?([0-9]+)_to_([0-9]+).*$"
  
  args$treeFileExtension <- ".treefile"
  args$csvFileExtension = NULL
  args$alignmentFileExtension = NULL
  args$pdfWidth = NULL
  args$pdfRelHeight = NULL
  args$pdfScaleBarWidth =NULL
  args$outputRDA =T
  args$RDAonly =F
  args$seed =42
  args$overwrite =T
  
  args$normRefFileName ="~/normalisation_ByPosition.csv"
  args$normStandardiseGagPol=T
  args$normalisationConstants=NULL
  args$duplicateBlacklist=NULL
  args$parsimonyBlacklistK=15
  args$rawBlacklistThreshold=3
  args$ratioBlacklistThreshold=0.01
  args$dualBlacklist=NULL
  args$minReadsPerHost=30
  args$minTipsPerHost=1
  
  args$postHocCountBlacklisting <- F
  args$maxReadsPerHost <- 100
  args$blacklistUnderrepresented
  args$blacklistReport = F
  args$useff=NULL
  args$splitsRule='s,15'
  args$pruneBlacklist=NULL
  args$readCountsMatterOnZeroLengthBranches=NULL
  args$outputNexusTree=NULL
  args$recombinationFiles=NULL
  
  args$allClassifications = F
  args$collapsedTrees=F
  args$windowThreshold=0.5
  args$distanceThreshold=c(0.02,0.05)
  args$allowMultiTrans=T
  args$relaxedAncestry=T
  args$multinomial=T
  args$directionThreshold=0.33
  args$summaryPlotDimensions=NULL
  args$skipSummaryGraph=F
  
  args$outputNexusTree = F
  
  
  tree.file.directory=tree.directory
  # tree.file.regex = "^RAxML_bestTree.InWindow_([0-9]+_to_[0-9]+)\\.tree$"
  tree.file.regex = tree.file.regex
  splits.rule = reconstruction.mode
  # sankoff.k = 0
  sankoff.unassigned.switch.threshold = 0
  continuation.unassigned.proximity.cost = 1000
  # outgroup.name = NULL
  multifurcation.threshold = m.thresh
  guess.multifurcation.threshold = is.na(m.thresh)
  # user.blacklist.directory = NULL
  # user.blacklist.file.regex = NULL
  # duplicate.file.directory = NULL
  # duplicate.file.regex = "^DuplicateReadCountsProcessed_InWindow_([0-9]+_to_[0-9]+).csv$"
  recombination.file.directory = recomb.file.directory
  recombination.file.regex = recomb.file.regex
  alignment.file.directory = alignment.directory
  # alignment.file.regex = NULL
  # tip.regex = "^(.*)_read_([0-9]+)_count_([0-9]+)$"
  # file.name.regex = "^(?:.*\\D)?([0-9]+)_to_([0-9]+).*$"
  # seed = sample(1:10000000,1)
  # norm.ref.file.name = NULL
  norm.standardise.gag.pol = norm.standardise.gp
  norm.constants = norm.constants.input
  # allow.mt = F
  relaxed.ancestry = F
  parsimony.blacklist.k = par.blacklisting.k
  raw.blacklist.threshold = bl.raw.threshold
  ratio.blacklist.threshold = bl.ratio.threshold
  # do.dual.blacklisting = F
  max.reads.per.host = downsampling.limit
  blacklist.underrepresented = blacklist.ur
  min.reads.per.host = ifelse(post.hoc.min.counts, 1, min.reads.per.host)
  min.tips.per.host = ifelse(post.hoc.min.counts, 1, min.tips.per.host)
  use.ff = useff
  # prune.blacklist = F
  count.reads.in.parsimony = read.counts.matter
  # verbosity = 0
  # no.progress.bars = F
  
  
  do.dual.blacklisting=F
  count.reads.in.parsimony=F
  
  # bl
  parsimony.blacklist.k = if(do.par.blacklisting) parsimony.blacklist.k else 0
  
  
}

additional_individual <- function(){
  # AID5977 - AID2168
  library(ape)
  dir <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI//210325_phsc_output_check_bl3rd/ptyr389_trees'
  files <- list.files(dir,pattern = '*.treefile',full.names = T)
  for (file in files) {
    cat('processing ',file, '...\n')
    tree <- read.tree(file)
    # tree$tip.label[grepl('^AID4448-fq[0-9]+_read_[0-9]+_count_[0-9]+$',tree$tip.label)] <- paste0('REMOVE_',grep('^AID4448-fq[0-9]+_read_[0-9]+_count_[0-9]+$',tree$tip.label,value = T))
    tree$tip.label[grepl('^AID4448-fq[0-9]+_read_[0-9]+_count_[0-9]+$',tree$tip.label)] <- paste0(grep('^AID4448-fq[0-9]+_read_[0-9]+_count_[0-9]+$',tree$tip.label,value = T),'_REMOVE')
     write.tree(tree, file = file)
  }
  
  library(data.table)
  library(tidyverse)
  library(glue)
  library(igraph)
  library(RBGL)
  source('~/transmission_network_functions_phyloscanner.R')
  # optional: meta data
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_output_check_bl3rd'
  infiles	<- list.files(indir, pattern='*workspace*', full.names=TRUE)
  control <- list(linked.group='close.and.adjacent.cat',linked.no='not.close.or.nonadjacent',linked.yes='close.and.adjacent', conf.cut=0.6, neff.cut=3,weight.complex.or.no.ancestry=0.5)
  
  for (infile in infiles) {
    load(infile)
    dca <- copy(data.table(dc))
    dwina <- copy(data.table(dwin))
    if(grepl('ptyr204',infile)){
      run <- 204
    }else if(grepl('ptyr249',infile)){
      run <- 249
    }else if(grepl('ptyr389',infile)){
      run <- 389
    }
    dca[,PTY_RUN:=run]
    dwina[,PTY_RUN:=run]
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
    
    # find pairs
    tmp <- find.pairs.in.networks(dwina, dca, control=control, verbose=TRUE)
    dpl <- copy(tmp$network.pairs)
    dc <- copy(tmp$relationship.counts)
    dw <- copy(tmp$windows)
    save(dpl, dc, dw, file=gsub('workspace','phscnetworks_allpairs',infile))
    
    # find networks
    tmp <- find.networks(dc, control=control, verbose=TRUE)
    dnet <- copy(tmp$transmission.networks)
    dchain <- copy(tmp$most.likely.transmission.chains)
    save(dpl, dc, dw, dnet, dchain, file=gsub('workspace','phscnetworks',infile))
  }

  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_output_check_bl3rd/'
  infile.networks1 <-  gsub('workspace','phscnetworks',infiles[2])
  infile.networks2 <-  gsub('workspace','phscnetworks',infiles[3])
  infile.networks3 <-  gsub('workspace','phscnetworks',infiles[4])
  infile.networks4 <-  gsub('workspace','phscnetworks',infiles[5])
  infile.networks5 <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/Rakai_phscnetworks.rda'
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  
  load( infile.networks1)
  dw1 <- copy(dw)
  dw1 <- dw1[,c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN')]
  setnames( dw1, c('H1','H2'), c('ID1','ID2'))
  tmp <- dw1[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  set(tmp, NULL, 'BASIC_CLASSIFICATION', tmp[,gsub('xx','desc_',gsub('desc_','anc_',gsub('anc_','xx',BASIC_CLASSIFICATION)))])
  dw1 <- unique(rbind(dw1[ID1<ID2],tmp))
  load( infile.networks2)
  dw2 <- copy(dw)
  dw2 <- dw2[,c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN')]
  setnames( dw2, c('H1','H2'), c('ID1','ID2'))
  tmp <- dw2[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  set(tmp, NULL, 'BASIC_CLASSIFICATION', tmp[,gsub('xx','desc_',gsub('desc_','anc_',gsub('anc_','xx',BASIC_CLASSIFICATION)))])
  dw2 <- unique(rbind(dw2[ID1<ID2],tmp))
  load( infile.networks3)
  dw3 <- copy(dw)
  dw3 <- dw3[,c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN')]
  setnames( dw3, c('H1','H2'), c('ID1','ID2'))
  tmp <- dw3[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  set(tmp, NULL, 'BASIC_CLASSIFICATION', tmp[,gsub('xx','desc_',gsub('desc_','anc_',gsub('anc_','xx',BASIC_CLASSIFICATION)))])
  dw3 <- unique(rbind(dw3[ID1<ID2],tmp))
  load( infile.networks4)
  dw4 <- copy(dw)
  dw4 <- dw4[,c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN')]
  setnames( dw4, c('H1','H2'), c('ID1','ID2'))
  tmp <- dw4[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  set(tmp, NULL, 'BASIC_CLASSIFICATION', tmp[,gsub('xx','desc_',gsub('desc_','anc_',gsub('anc_','xx',BASIC_CLASSIFICATION)))])
  dw4 <- unique(rbind(dw4[ID1<ID2],tmp))
  aid <- data.table(read.csv(infile.ind.anonymised))
  load( infile.networks5)
  dw5 <- copy(dw)
  dw5 <- dw5[,c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN')]
  setnames( dw5, c('H1','H2'), c('ID1','ID2'))
  tmp <- dw5[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  set(tmp, NULL, 'BASIC_CLASSIFICATION', tmp[,gsub('xx','desc_',gsub('desc_','anc_',gsub('anc_','xx',BASIC_CLASSIFICATION)))])
  dw5 <- unique(rbind(dw5[ID1<ID2],tmp))
  aid <- data.table(read.csv(infile.ind.anonymised))

  
  inclusion <- "both"# "either"
  setnames(dw1, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  setnames(dw2, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  setnames(dw3, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  setnames(dw4, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  setnames(dw5, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  
  library(Phyloscanner.R.utilities)
  library(phyloscannerR)
  control = list(	yintercept_close=0.025,
                  yintercept_dist=1.0,
                  breaks_x=seq(0,1e4,500),
                  minor_breaks_x=seq(0,1e4,100),
                  breaks_y=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),
                  limits_y=c(1e-3,0.4),
                  fill.topology=c("ancestral"="deepskyblue1","descendant"="deepskyblue4","intermingled"= "#FDB863",'sibling'="#8073AC","other"="grey80"))

  hosts <- c('AID0719','AID1311')
  tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw5), inclusion = "both",control=control)
  g1<- tmp$graph +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw1), inclusion = "both",control=control)
  g2<- tmp$graph +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

  library(ggpubr)
  arrange <- ggarrange(g1,g2,
                       ncol = 1, nrow = 2, common.legend = T, legend='bottom')
  ggsave(paste0('~/test_',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 8)

  
  hosts<- c('AID5977', 'AID2168')
  tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw2), inclusion = "both",control=control)
  g1<- tmp$graph +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw3), inclusion = "both",control=control)
  g2<- tmp$graph +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw4), inclusion = "both",control=control)
  g3<- tmp$graph +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  
  library(ggpubr)
  arrange <- ggarrange(g1,g2,g3, 
                       ncol = 1, nrow = 3, common.legend = T, legend='bottom')
  ggsave(paste0('~/test_',hosts[1],'_',hosts[2],'.pdf'),arrange,width = 8, height = 12)
  
  # pty.runs[grepl('AID0719|AID1311',RENAME_ID)]
  # unique(subset(dw1,select=c('host.1','host.2')))
  # dw1[host.1=='AID0719'&host.2=='AID2168']

}