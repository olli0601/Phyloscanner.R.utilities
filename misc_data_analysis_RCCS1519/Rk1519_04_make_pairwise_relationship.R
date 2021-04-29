
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
  # outdir	<- file.path(HOME,'210325_phsc_phscrelationships_02_05/')
  # sa1
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
  control$max.reads.per.host <- 100
  # control$min.reads.per.host <- 30
  control$min.reads.per.host <- NULL
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
  control$ratio.blacklist.threshold = 0.01
  control$raw.blacklist.threshold = 3			
  control$recombination.file.directory = NULL
  control$recombination.file.regex = NULL
  control$relaxed.ancestry = TRUE
  control$sankoff.k = 15
  control$sankoff.unassigned.switch.threshold = 0
  control$seed = 42
  control$splits.rule = 's'
  control$tip.regex = '^(.*)-fq[0-9]+_read_([0-9]+)_count_([0-9]+)'
  control$tree.file.regex = "^(.*)\\.treefile$" # from rscript
  control$treeFileExtension = '.treefile'
  control$use.ff = FALSE
  control$user.blacklist.directory = NULL 
  control$user.blacklist.file.regex = NULL
  control$verbosity = 1	
  
  
  #	make bash for many files	
  df <- tibble(F=list.files(treedir))
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
phsc.migrate.transmission.networks<- function()
{
  library(data.table)
  source('~/transmission_network_functions_phyloscanner.R')
  # optional: meta data
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05/'
  # sa1
  # indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/' #no min no max
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
  # save(dca,dwina,file = '~/dcdwina_minread30.rda')
  # save(dca,dwina,file = '~/dcdwina_minreadnull.rda')
  # save(dca,dwina,file = '~/dcdwina_minread30_maxread100.rda')
  # save(dca,dwina,file = '~/dcdwina_minreadnull_maxread100.rda')
  # load('~/dcdwina_minread30.rda')
  # load('~/dcdwina_minreadnull.rda')
  load('~/dcdwina_minread30_maxread100.rda')
  load('~/dcdwina_minreadnull_maxread100.rda')
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
  
  tmp			<- subset(dca, host.1>host.2)
  setnames(tmp, c('host.1','host.2','CNTRL1','CNTRL2'),
           c('host.2','host.1','CNTRL2','CNTRL1'))
  set(tmp, NULL, 'type',
      tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',type)))])
  dca		<- rbind(subset(dca, !(host.1>host.2)), tmp)  
  
  # tmp <- dca[!(host.1<host.2)]
  # setnames(tmp,c('host.1','host.2','CNTRL1','CNTRL2'),c('host.2','host.1','CNTRL2','CNTRL1'))
  # dca <- rbind(dca[(host.1<host.2)],tmp)
  # tmp <- dwina[!(host.1<host.2)]
  # setnames(tmp,c('host.1','host.2','CNTRL1','CNTRL2','paths12','paths21','nodes1','nodes2'),c('host.2','host.1','CNTRL2','CNTRL1','paths21','paths12','nodes2','nodes1'))
  # tmp[,close.and.contiguous.and.directed.cat:=ifelse(close.and.contiguous.and.directed.cat=='12'|close.and.contiguous.and.directed.cat=='21',ifelse(close.and.contiguous.and.directed.cat=='12','21','12'),close.and.contiguous.and.directed.cat)]
  # tmp[,close.and.contiguous.and.ancestry.cat:=ifelse(close.and.contiguous.and.ancestry.cat=='12'|close.and.contiguous.and.ancestry.cat=='21',ifelse(close.and.contiguous.and.ancestry.cat=='12','21','12'),close.and.contiguous.and.ancestry.cat)]
  # tmp[, close.and.adjacent.and.ancestry.cat:=ifelse( close.and.adjacent.and.ancestry.cat=='12'| close.and.adjacent.and.ancestry.cat=='21',ifelse( close.and.adjacent.and.ancestry.cat=='12','21','12'), close.and.adjacent.and.ancestry.cat)]
  # dwina <- rbind(dwina[(host.1<host.2)],tmp)
  
  # tmp <- find.pairs.in.networks(indir, batch.regex='^ptyr([0-9]+)_.*', control=control, verbose=TRUE, dmeta=dmeta)
  library(tidyverse)
  # dwin <- copy(dwina)
  # dc <- copy(dca)
  # verbose=TRUE
  # dmeta=NULL
  tmp <- find.pairs.in.networks(dwina, dca, control=control, verbose=TRUE)
  # dpl <- copy(tmp$linked.pairs)
  dpl <- copy(tmp$network.pairs)
  dc <- copy(tmp$relationship.counts)
  dw <- copy(tmp$windows)
  save(dpl, dc, dw, file=file.path(indir,'Rakai_phscnetworks_allpairs.rda'))
  
  library(glue)
  library(igraph)
  library(RBGL)
  tmp <- find.networks(dc, control=control, verbose=TRUE)
  dnet <- copy(tmp$transmission.networks)
  dchain <- copy(tmp$most.likely.transmission.chains)
  save(dpl, dc, dw, dnet, dchain, file=file.path(indir,'Rakai_phscnetworks.rda'))
  
  
  load('~/dcdwina_minreadnull.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  aid <- data.table(read.csv(infile.ind.anonymised))
  
  tmp <- subset(aid,select=c('PT_ID','AID'))
  dwina <- merge(dwina, tmp, by.x='host.1',by.y='AID',all.x=T)
  dwina <- merge(dwina, tmp, by.x='host.2',by.y='AID',all.x=T)
  setnames(dwina, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  
  # X      PT_ID     AID
  # 1: 5231 RK-J104547 AID2168
  # X      PT_ID     AID
  # 1: 1500 RK-C105106 AID5977  
  aid[AID=='AID2168']
  aid[AID=='AID5977']
  tmp <- dwina[ID1=='RK-J104547' & ID2== 'RK-C105106',c('tree.id','close.and.adjacent.and.ancestry.cat','PTY_RUN')]
  tmp[,{
    tb <- table(close.and.adjacent.and.ancestry.cat)
    list(COUNT=as.vector(tb), CAT=names(tb))},by='PTY_RUN']
  
  infile.aid19 <- '~/todi_pairs_171122_cl25_d50_prior23_min30_anonymised_RIDs.csv'
  infile.pairs19 <- '~/Rakai_phscnetworks_allpairs_190706.rda'
  # infile.net19 <- '~/Rakai_phscnetworks_190706.rda'
  aid19 <- data.table(read.csv(infile.aid19))
  aid19[, ID:=paste0('RK-',ID)]
  load(infile.pairs19)
  dw19 <- copy(dw)
  dw19 <- merge(dw19, aid19, by.x='H1',by.y='AID',all.x=T)
  dw19 <- merge(dw19, aid19, by.x='H2',by.y='AID',all.x=T)
  setnames(dw19, c('ID.x','ID.y'),c('ID1','ID2'))
  dw19 <- as.data.table(dw19)
  tmp19 <- dw19[ID2=='RK-J104547' & ID1== 'RK-C105106',c('TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT','PTY_RUN')]
  tmp19[,{
    tb <- table(CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT)
    list(COUNT=as.vector(tb), CAT=names(tb))},by='PTY_RUN']
}


analysis_network <- function(){
  library(data.table)
  # dirs
  # indir			<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05/'
  # # sa1
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/'
  # indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read/' #30 min 100 max
  # indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read_100_max_read/' #no min 100 max
  # 
  # infile.pairs	<- file.path(indir,'Cluster_strongsupport_allwindows.rda')
  # infile.networks	<- file.path(indir,'Cluster_strongsupport_networksallpairs.rda')
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  infile.ind.rccs <- file.path('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.mrc <- file.path('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
  infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
  load( infile.networks )	# loads rtp, rplkl, rpw
  
  
  # make table mapping id and aid, aid and genders
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
  # dchain[H1=='AID2324'& H2== 'AID2368']
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
  # couple[,SEX1:='M']
  # couple[,SEX2:='F']
  couple <- merge(couple,subset(df,select=c('PT_ID', 'AID')),by.x='PT_ID1', by.y='PT_ID', all.x=T)
  couple <- merge(couple,subset(df,select=c('PT_ID', 'AID')),by.x='PT_ID2', by.y='PT_ID', all.x=T)
  set(couple,NULL,c('PT_ID1','PT_ID2'),NULL)
  setnames(couple,c('AID.x','AID.y'),c('H1','H2'))
  couple[,H1:=as.character(H1)]
  couple[,H2:=as.character(H2)]
  couple <- couple[!is.na(H1) & !is.na(H2)]
  # tmp <- copy(couple)
  # setnames(tmp,c('H1','H2'),c('H2','H1'))
  # couple <- unique(rbind(couple,tmp))
  
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
  
  
  # make dobs
  # dobs <- subset(rtp, select=c('ID1','ID2','LINK_12','LINK_21','ID1_SEX','ID2_SEX','COUPLE'))
  # tmp <- dobs[LINK_12==1]
  # tmp <- subset(tmp,select=c('ID1','ID2','ID1_SEX','ID2_SEX','COUPLE'))
  # setnames(tmp,c('ID1','ID2','ID1_SEX','ID2_SEX'), c('TR_ID','REC_ID','TR_SEX','REC_SEX'))
  # dobs <- dobs[LINK_12==0]
  # dobs <- subset(dobs,select=c('ID1','ID2','ID1_SEX','ID2_SEX','COUPLE'))  
  # setnames(dobs,c('ID1','ID2','ID1_SEX','ID2_SEX'), c('REC_ID','TR_ID','REC_SEX','TR_SEX'))
  # dobs <- rbind(dobs, tmp)
  
  
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
  # load run
  infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
  infile.run='/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_phscinput_runs.rds'
  
  pty.runs <- data.table(readRDS(infile.run))
  
  # map id and sequences
  data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  infile.ind.rccs <- file.path(data.dir,'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
  infile.ind.mrc <- file.path(data.dir,'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
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
  infile.aid19 <- '~/todi_pairs_171122_cl25_d50_prior23_min30_anonymised_RIDs.csv'
  infile.pairs19 <- '~/Rakai_phscnetworks_allpairs_190706.rda'
  infile.net19 <- '~/Rakai_phscnetworks_190706.rda'
  aid19 <- data.table(read.csv(infile.aid19))
  aid19[, ID:=paste0('RK-',ID)]
  load(infile.pairs19)
  load(infile.net19)
  dc19 <- copy(dc)
  dw19 <- copy(dw)
  dnet19 <- copy(dnet)
  dchain19 <- copy(dchain)
  
  
  tmp <- subset(aid19,select=c('ID','AID'))
  dc19 <- merge(dc19, tmp, by.x='H1',by.y='AID',all.x=T)
  dc19 <- merge(dc19, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dc19, c('ID.x','ID.y'),c('ID1','ID2'))
  dw19 <- merge(dw19, tmp, by.x='H1',by.y='AID',all.x=T)
  dw19 <- merge(dw19, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dw19, c('ID.x','ID.y'),c('ID1','ID2'))
  dnet19 <- merge(dnet19, tmp, by.x='H1',by.y='AID',all.x=T)
  dnet19 <- merge(dnet19, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dnet19, c('ID.x','ID.y'),c('ID1','ID2'))
  dchain19 <- merge( dchain19, tmp, by.x='H1',by.y='AID',all.x=T)
  dchain19 <- merge( dchain19, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames( dchain19, c('ID.x','ID.y'),c('ID1','ID2'))
  
  # load network
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/'
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
  load( infile.networks )
  aid <- data.table(read.csv(infile.ind.anonymised))
  
  tmp <- subset(aid,select=c('PT_ID','AID'))
  dc <- merge(dc, tmp, by.x='H1',by.y='AID',all.x=T)
  dc <- merge(dc, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dc, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  dw <- merge(dw, tmp, by.x='H1',by.y='AID',all.x=T)
  dw <- merge(dw, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dw, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  dchain <- merge(dchain, tmp, by.x='H1',by.y='AID',all.x=T)
  dchain <- merge(dchain, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dchain, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  dnet <- merge(dnet, tmp, by.x='H1',by.y='AID',all.x=T)
  dnet <- merge(dnet, tmp, by.x='H2',by.y='AID',all.x=T)
  setnames(dnet, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  
  # whether couples in dc
  tmp <- unique(subset(dc,select = c('ID1','ID2')))
  tmp[,ID.x:=paste0(ID1,'-',ID2)]
  tmp[,ID.y:=paste0(ID2,'-',ID1)]
  cat(nrow(couple[paste0(ID1,'-',ID2)%in% tmp$ID.x | paste0(ID1,'-',ID2)%in% tmp$ID.y ]), ' couples in networks')
  
  # to dt
  dc19 <- as.data.table(dc19)
  dw19 <- as.data.table(dw19)
  dnet19 <- as.data.table(dnet19)
  dchain19 <- as.data.table(dchain19)
  dchain <- as.data.table(dchain)
  
  # dw19[ID1=='RK-H180029' & ID2== 'RK-B104396']
  # tmp <- dw19[ID2=='RK-H180029' & ID1== 'RK-B104396',c('TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT')]
  # dw[ID1=='RK-H180029' & ID2== 'RK-B104396']
  # tmp <- merge(tmp, 
  #              dw[ID2=='RK-H180029' & ID1== 'RK-B104396',c('TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT')],
  # by='TREE_ID',all=T)
  # dw19[ID2=='RK-H180029' & ID1== 'RK-B104396',table(CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT)]
  # dw[ID2=='RK-H180029' & ID1== 'RK-B104396',table(CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT)]
  
  # add couple
  tmp <- copy(couple)
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  couple <- data.table(unique(rbind(couple, tmp)))
  
  dc <- merge(dc,couple,by=c('ID1','ID2'),all.x=T)
  dc19 <- merge(dc19,couple,by=c('ID1','ID2'),all.x=T)
  dw <- merge(dw,couple,by=c('ID1','ID2'),all.x=T)
  dw19 <- merge(dw19,couple,by=c('ID1','ID2'),all.x=T)
  dnet <- merge(dnet,couple,by=c('ID1','ID2'),all.x=T)
  dnet19 <- merge(dnet19,couple,by=c('ID1','ID2'),all.x=T)
  dchain <- merge(dchain,couple,by=c('ID1','ID2'),all.x=T)
  dchain19 <- merge(dchain19,couple,by=c('ID1','ID2'),all.x=T)
  
  #
  tmp_dc_coup <- unique(dc[CATEGORISATION=='close.and.adjacent.and.ancestry.cat'&COUPLE==1])
  tmp_dc_coup[,length(N),by=c('ID1','ID2')]
  tmp_dc_coup[ID1=='RK-A035039' & ID2=='RK-D061829']
  tmp <- tmp_dc_coup[,list(N_EFF=max(N_EFF)),by=c('ID1','ID2')]
  tmp_dc_coup <- merge(tmp_dc_coup,tmp, by=c('ID1','ID2','N_EFF'))
  setkey(tmp_dc_coup, ID1, ID2, PTY_RUN, TYPE)
  tmp_dc_coup <- tmp_dc_coup[,head(.SD, 4),by=c('ID1','ID2')]
  tmp_dc_coup <- dcast(tmp_dc_coup,  ID1 + ID2 + COUPLE ~ TYPE, value.var= 'SCORE')
  tmp_dc_coup[,ANALYSIS:=2021]
  
  tmp_dc19_coup <- unique(dc19[CATEGORISATION=='close.and.adjacent.and.ancestry.cat'&COUPLE==1])
  # tmp_dc19_coup[,length(N),by=c('ID1','ID2')]
  # tmp <- tmp_dc19_coup[,list(N_EFF=max(N_EFF)),by=c('ID1','ID2')]
  # tmp_dc19_coup <- merge(tmp_dc19_coup,tmp, by=c('ID1','ID2','N_EFF'))
  tmp_dc19_coup <- dcast(tmp_dc19_coup,  ID1 + ID2 + COUPLE ~ TYPE, value.var= 'SCORE')
  tmp_dc19_coup[,ANALYSIS:=2019]
  
  tmp_dc_coup <- rbind(tmp_dc19_coup, tmp_dc_coup)
  tmp <- tmp_dc_coup[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2','12','21'),c('ID2','ID1','21','12'))
  tmp_dc_coup <- rbind(tmp_dc_coup[(ID1<ID2)],tmp)
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
  tmp1 <- tmp_dc_coup[is.na(`2019`) & !is.na(`2021`) ]
  tmp2 <- tmp_dc_coup[!is.na(`2019`) & is.na(`2021`) ]
  tmp3 <- tmp_dc_coup[!is.na(`2019`) & !is.na(`2021`) ]
  nrow(unique(subset(tmp1,select=c('ID1','ID2'))))
  nrow(unique(subset(tmp2,select=c('ID1','ID2'))))
  nrow(unique(subset(tmp3,select=c('ID1','ID2'))))
  # 2021 ONLY
  tmp <- dcast(tmp1,ID1+ID2~variable,value.var='2021')
  tmp[,LINKED:= `12`+`21` +complex.or.no.ancestry>0.6]
  tmp[LINKED==1,DIRECTED:= (`12`+`21`)/(`12`+`21` +complex.or.no.ancestry)>0.6]
  tmp[DIRECTED==1,LINK12:=`12`>`21`]
  nrow(tmp[LINKED==1])/nrow(tmp)
  nrow(tmp[DIRECTED==1])/nrow(tmp[LINKED==1])
  # 2019 ONLY
  tmp <- dcast(tmp2,ID1+ID2~variable,value.var='2019')
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
  
  load('~/dcdwina_minreadnull.rda')
  dca[,host.1:=gsub('CNTRL-','',host.1)]
  dca[,host.2:=gsub('CNTRL-','',host.2)]
  dcac <- merge(dca, aid, by.x='host.1',by.y='AID',all.x=T)
  dcac <- merge(dcac, aid, by.x='host.2',by.y='AID',all.x=T)
  dcac[paste0(PT_ID.x,'-',PT_ID.y) %in% tmp$ID.x & categorisation == 'close.and.adjacent.and.ancestry.cat']
  dcac[paste0(PT_ID.x,'-',PT_ID.y) %in% tmp$ID.y & categorisation == 'close.and.adjacent.and.ancestry.cat']
  
  # 2019 and 2021
  tmp21 <- dcast(tmp3,ID1+ID2~variable,value.var='2021')
  tmp21[,LINKED_SCORE:=`12`+`21` +complex.or.no.ancestry]
  tmp21[,LINKED:= LINKED_SCORE>0.6]
  tmp21[LINKED==1,DIRECTED_SCORE:=(`12`+`21`)/LINKED_SCORE]
  tmp21[LINKED==1,DIRECTED:= DIRECTED_SCORE>0.6]
  tmp21[DIRECTED==1,LINK12_SCORE:=`12`/(`12`+`21`)]
  tmp21[DIRECTED==1,LINK12:=`12`>`21`]
  nrow(tmp21[LINKED==1])/nrow(tmp21)
  nrow(tmp21[DIRECTED==1])/nrow(tmp21[LINKED==1])
  
  tmp19 <- dcast(tmp3,ID1+ID2~variable,value.var='2019')
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
  # tmp[LINKED.x==LINKED.y,`12.x`+`21.x` +complex.or.no.ancestry.x]
  # tmp[LINKED.x==LINKED.y,`12.y`+`21.y` +complex.or.no.ancestry.y]
  # tmp[LINKED.x!=LINKED.y,`12.x`+`21.x` +complex.or.no.ancestry.x]
  # tmp[LINKED.x!=LINKED.y,`12.y`+`21.y` +complex.or.no.ancestry.y]
  
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
  dw19_couple <- unique(rbind(dw19_couple[ID1<ID2],tmp))
  
  dw_couple <- dw[COUPLE==1,c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN')]
  setnames( dw_couple, c('H1','H2'), c('ID1','ID2'))
  tmp <- dw_couple[,min(PTY_RUN), by=c('ID1','ID2')]
  
  # tmp[(ID1=='AID6347' & ID2=='AID6314')|(ID2=='AID6347' & ID1=='AID6314')]
  # tmp[(ID1=='AID5758' & ID2=='AID7973')|(ID2=='AID5758' & ID1=='AID7973')]
  # tmp[(ID1=='AID5977' & ID2=='AID2168')|(ID2=='AID5977' & ID1=='AID2168')]
  # # ID1     ID2  V1
  # # 1: AID6314 AID6347 119
  # # ID1     ID2 V1
  # # 1: AID5758 AID7973 72
  # # ID1     ID2  V1
  # # 1: AID2168 AID5977 389
  
  dw_couple <-  merge(dw_couple,tmp, by.x=c('ID1','ID2','PTY_RUN'),by.y=c('ID1','ID2','V1'))
  tmp <- dw_couple[!(ID1<ID2)]
  setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
  dw_couple <- unique(rbind(dw_couple[ID1<ID2],tmp))
  dw19c <- merge(dw19, aid,by.x='ID1',by.y='PT_ID',all.x=T)
  dw19c <- merge(dw19c, aid,by.x='ID2',by.y='PT_ID',all.x=T)
  
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read_100_max_read/' #30 min 100 max
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  load(infile.networks)
  dw_maxhost <- merge(dw, aid, by.x='H1',by.y='AID',all.x=T)
  dw_maxhost <- merge(dw_maxhost, aid, by.x='H2',by.y='AID',all.x=T)
  setnames(dw_maxhost, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read/' 
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  load(infile.networks)
  dw_minmaxhost <- merge(dw, aid, by.x='H1',by.y='AID',all.x=T)
  dw_minmaxhost <- merge(dw_minmaxhost, aid, by.x='H2',by.y='AID',all.x=T)
  setnames(dw_minmaxhost, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
  
  dw_couple[(ID1=='AID6347' & ID2=='AID6314')|(ID2=='AID6347' & ID1=='AID6314')]
  dw_couple[(ID1=='AID5758' & ID2=='AID7973')|(ID2=='AID5758' & ID1=='AID7973')] # 72
  dw_couple[(ID1=='AID5977' & ID2=='AID2168')|(ID2=='AID5977' & ID1=='AID2168')] # 389
  
  unique(dc[(H1=='AID5758' & H2=='AID7973')|(H2=='AID5758' & H1=='AID7973'),c('PTY_RUN','N_EFF')])
  unique(dc[(H1=='AID5977' & H2=='AID2168')|(H2=='AID5977' & H1=='AID2168'),c('PTY_RUN','N_EFF')])
  
  pty.runs[grepl('AID6347|AID6314',RENAME_ID)]
  pty.runs[grepl('AID5758|AID7973',RENAME_ID)]
  pty.runs[grepl('AID5977|AID2168',RENAME_ID)]
  
  tmpdw <- dw[(H1=='AID5758' & H2=='AID7973')|(H2=='AID5758' & H1=='AID7973'),c('TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT')]
  tmpdw19 <- dw19c[(AID.x=='AID5758' & AID.y=='AID7973')|(AID.y=='AID5758' & AID.x=='AID7973'),c('TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT')]
  tmpdw <- merge(tmpdw19,tmpdw, by='TREE_ID',all=T)
  tmpdw[is.na( CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT.x )]  
  tmpdw[is.na( CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT.y )] 
  tmpdw_maxh <- dw_maxhost[(H1=='AID5758' & H2=='AID7973')|(H2=='AID5758' & H1=='AID7973'),c('TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT')]
  tmpdw_minmaxh <- dw_minmaxhost[(H1=='AID5758' & H2=='AID7973')|(H2=='AID5758' & H1=='AID7973'),c('TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT')]
  tmpdw <- merge(tmpdw,tmpdw_maxh, by='TREE_ID',all=T)
  tmpdw <- merge(tmpdw,tmpdw_minmaxh, by='TREE_ID',all=T)
  colnames(tmpdw) <- c("TREE_ID","CAT19",'CAT','CAT_MH','CAT_MMH')
  tmpdw[is.na(CAT_MH)]
  tmpdw[is.na(CAT_MMH)]
  tmpdc <- dc[(H1=='AID5758' & H2=='AID7973')|(H2=='AID5758' & H1=='AID7973'),]
  
  tmpdw <- unique(dw[(H1=='AID5977' & H2=='AID2168')|(H2=='AID5977' & H1=='AID2168'),c('TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT','PTY_RUN')])
  tmpdw <- dcast(tmpdw,TREE_ID~PTY_RUN, value.var = 'CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT')
  tmpdw19 <- (dw19c[(AID.x=='AID5977' & AID.y=='AID2168')|(AID.y=='AID5977' & AID.x=='AID2168'),c('TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT','PTY_RUN')])
  tmpdw <- merge(tmpdw19,tmpdw, by=c('TREE_ID'),all=T)
  tmpdw[CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT!=`389`]
  table(tmpdw$`389`)
  table(tmpdw$`404`)
  table(tmpdw$`409`)
  table(tmpdw$CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT)
  unique(subset(tmpdw,select='TREE_ID'))
  
  dw[(H1=='AID5977' & H2=='AID2168')|(H2=='AID5977' & H1=='AID2168'),]
  dw19[(ID1=='RK-J104547' & ID2=='RK-C105106')|(ID2=='RK-J104547' & ID1=='RK-C105106'),]
  dw[(ID1=='RK-J104547' & ID2=='RK-C105106')|(ID2=='RK-J104547' & ID1=='RK-C105106'),]
  
  # add aid
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
  
  
  # AID0011-AID1640
  tmp <- dc[(H1=='AID0011' & H2=='AID1640')|(H2=='AID0011' & H1=='AID1640') ]
  tmp <- tmp[CATEGORISATION == 'close.and.adjacent.and.ancestry.cat']
  setkey(tmp,TYPE)
  # TYPE  K      K_EFF   N N_EFF      SCORE CNTRL1 CNTRL2
  # 1:                       12 33  4.6427137 185  24.8 0.18720620  FALSE  FALSE
  # 2:                       21 92 12.6685821 185  24.8 0.51082992  FALSE  FALSE
  # 3:   complex.or.no.ancestry 53  6.5354749 185  24.8 0.26352721  FALSE  FALSE
  # 4: not.close.or.nonadjacent  7  0.9532293 185  24.8 0.03843667  FALSE  FALSE
  dc19c <- merge(dc19, aid,by.x='ID1',by.y='PT_ID',all.x=T)
  dc19c <- merge(dc19c, aid,by.x='ID2',by.y='PT_ID',all.x=T)
  
  tmp <- dc19c[(AID.x=='AID0011' & AID.y=='AID1640')|(AID.y=='AID0011' & AID.x=='AID1640') ]
  tmp <- tmp[CATEGORISATION == 'close.and.adjacent.and.ancestry.cat']
  setkey(tmp,TYPE)
  
  # 1:                       12 16 2.129091 89  17.9 0.1189436      F
  # 2:                       21 32 3.723636 89  17.9 0.2080244      F
  # 3:   complex.or.no.ancestry 14 2.597273 89  17.9 0.1450990      F
  # 4: not.close.or.nonadjacent 27 9.450000 89  17.9 0.5279330      F
  # 
  pty.runs[grepl('AID0011|AID1640',RENAME_ID)]
  
  tmp1 <- dw[(H1=='AID0011' & H2=='AID1640')|(H2=='AID0011' & H1=='AID1640') ]
  setkey(tmp1, TREE_ID)
  tmp1 = tmp1[,c('TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT')]
  tmp1[,START:=as.numeric(gsub('([0-9]+)_to_([0-9]+)','\\1',TREE_ID))]
  tmp1[START>6000]
  
  
  tmp2 <- dw19c[(AID.x=='AID0011' & AID.y=='AID1640')|(AID.y=='AID0011' & AID.x=='AID1640') ]
  setkey(tmp2, TREE_ID)
  tmp2 = tmp2[,c('TREE_ID','CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT')]
  tmp2[,START:=as.numeric(gsub('([0-9]+)_to_([0-9]+)','\\1',TREE_ID))]
  tmp2[START>6000]
  
  tmp <- merge(tmp2,tmp1,by=c('TREE_ID','START'),all=T)
  setkey(tmp,START)
  write.csv(tmp,file = '~/AID0011_AID1640.csv')
  
  
  
  
  # AID2168-AID5977
  pty.runs[grepl('AID2168|AID5977',RENAME_ID),unique(PTY_RUN)]
  # 389, 404, 409
  # AID5758-AID7973
  pty.runs[grepl('AID5758|AID7973',RENAME_ID),]
  pty.runs[grepl('AID5758|AID7973',RENAME_ID),unique(PTY_RUN)]
  # 72, 204
  
  # v loop
  unique(subset(dw19_couple,select=c('host.1','host.2')))
  # unique(subset(dw_couple,select=c('host.1','host.2')))
  dw19_couple[,start:=as.integer(gsub('([0-9]+)_to_([0-9]+)','\\1',tree.id))]
  tmp <- dw19_couple[start>=6615-250 & start<=7636 & start!=6825 & start!=6850]
  tmp <- unique(subset(tmp,select=c('host.1','host.2')))
  tmp2 <- unique(subset(dw_couple,select=c('host.1','host.2')))
  tmp2[,host.x:=paste0(host.1,'-',host.2)]
  tmp2[,host.y:=paste0(host.2,'-',host.1)]
  tmp <- tmp[paste0(host.1,'-', host.2)%in% tmp2$host.x |  paste0(host.1,'-', host.2)%in% tmp2$host.y]
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
  dw_maxhost_couple <- unique(rbind(dw_maxhost_couple[ID1<ID2],tmp))
  
  setnames(dw_maxhost_couple, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
           c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
  
  # load 2021 minmaxreadhost
  indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read/' #30 min 100 max
  infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
  load(infile.networks)
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
    # geom_density()+theme_bw()+labs(x='patristic distance', y='count')+ coord_cartesian(xlim = c(0,0.05))
    geom_histogram(binwidth = 0.002)+theme_bw()+labs(x='patristic distance', y='count')
  # coord_cartesian(ylim = c(0,1e3))
  ggsave(filename = '~/pdistance_couple.pdf',width = 6, height = 4)
  
  ggplot(tmp, aes(patristic.distance))+
    # geom_density()+theme_bw()+labs(x='patristic distance', y='count')+ coord_cartesian(xlim = c(0,0.05))
    geom_histogram(binwidth = 0.002)+theme_bw()+labs(x='patristic distance', y='count') +coord_cartesian(ylim = c(0,1e3),xlim=c(0,0.1))
  # coord_cartesian(ylim = c(0,1e3))
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
          axis.ticks.x=element_blank())
  # +
  #   coord_cartesian(ylim=c(0,0.2))
  ggsave('~/patdist_vs_couple.pdf',width = 10, height = 4)
  
  tmpid <- subset(tmp,select = c('host.1','host.2','ID'))
  
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


# couple.analysis <- function(){
#   library(data.table)
#   # load run
#   infile.couple <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
#   infile.run='/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_phscinput_runs.rds'
# 
#   pty.runs <- data.table(readRDS(infile.run))
# 
#   # map id and sequences
#   data.dir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
#   infile.ind.rccs <- file.path(data.dir,'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
#   infile.ind.mrc <- file.path(data.dir,'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
#   id.dt <- data.table(read.csv(infile.ind.rccs))
#   id.dt <- subset(id.dt,select = c("pt_id","pangea_id"))
#   tmp <- data.table(read.csv(infile.ind.mrc))
#   tmp <- subset(tmp,select = c("pt_id","pangea_id"))
#   id.dt <- rbind(id.dt,tmp)
#   id.dt <- unique(id.dt)
# 
#   # load couple
#   load(infile.couple)
#   couple <- data.table(unique(subset(coupdat,select=c('male.RCCS_studyid', 'female.RCCS_studyid'))))
#   colnames(couple) <- c('ID1','ID2')
#   couple[,COUPLE:=1]
#   couple[,ID1:=paste0('RK-',ID1)]
#   couple[,ID2:=paste0('RK-',ID2)]
# 
# 
#   # whether couples in run
#   cat(nrow(couple), ' couples in total \n',
#       nrow(couple[ID1%in% pty.runs$UNIT_ID & ID2%in% pty.runs$UNIT_ID]), ' in runs \n',
#       nrow(couple[ID1%in% id.dt$pt_id & ID2%in% id.dt$pt_id]), ' have pangea ids \n')
#   tmp <- couple[ID1%in% id.dt$pt_id & ID2%in% id.dt$pt_id & !(ID1%in% pty.runs$UNIT_ID & ID2%in% pty.runs$UNIT_ID)]
# 
#   # eg
#   # WARNING: lose some couples ... check later
#   tmp[ID1%in% pty.runs$UNIT_ID & !ID2%in% pty.runs$UNIT_ID]
#   merge(tmp[ID1%in% pty.runs$UNIT_ID & !ID2%in% pty.runs$UNIT_ID],id.dt, by.x='ID2',by.y='pt_id',all.x=T)
#   tmp[!ID1%in% pty.runs$UNIT_ID & ID2%in% pty.runs$UNIT_ID]
#   tmp[!ID1%in% pty.runs$UNIT_ID & !ID2%in% pty.runs$UNIT_ID]
# 
# 
#   # load 19 network
#   infile.aid19 <- '~/todi_pairs_171122_cl25_d50_prior23_min30_anonymised_RIDs.csv'
#   infile.pairs19 <- '~/Rakai_phscnetworks_allpairs_190706.rda'
#   infile.net19 <- '~/Rakai_phscnetworks_190706.rda'
#   aid19 <- data.table(read.csv(infile.aid19))
#   aid19[, ID:=paste0('RK-',ID)]
#   load(infile.pairs19)
#   load(infile.net19)
#   dc19 <- copy(dc)
#   dw19 <- copy(dw)
#   dnet19 <- copy(dnet)
#   dchain19 <- copy(dchain)
# 
#   tmp <- subset(aid19,select=c('ID','AID'))
#   dc19 <- merge(dc19, tmp, by.x='H1',by.y='AID',all.x=T)
#   dc19 <- merge(dc19, tmp, by.x='H2',by.y='AID',all.x=T)
#   setnames(dc19, c('ID.x','ID.y'),c('ID1','ID2'))
#   dw19 <- merge(dw19, tmp, by.x='H1',by.y='AID',all.x=T)
#   dw19 <- merge(dw19, tmp, by.x='H2',by.y='AID',all.x=T)
#   setnames(dw19, c('ID.x','ID.y'),c('ID1','ID2'))
#   dnet19 <- merge(dnet19, tmp, by.x='H1',by.y='AID',all.x=T)
#   dnet19 <- merge(dnet19, tmp, by.x='H2',by.y='AID',all.x=T)
#   setnames(dnet19, c('ID.x','ID.y'),c('ID1','ID2'))
#   dchain19 <- merge( dchain19, tmp, by.x='H1',by.y='AID',all.x=T)
#   dchain19 <- merge( dchain19, tmp, by.x='H2',by.y='AID',all.x=T)
#   setnames( dchain19, c('ID.x','ID.y'),c('ID1','ID2'))
# 
#   # load network
#   indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read/'
#   infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
#   infile.ind.anonymised <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv'
#   load( infile.networks )
#   aid <- data.table(read.csv(infile.ind.anonymised))
# 
#   tmp <- subset(aid,select=c('PT_ID','AID'))
#   dc <- merge(dc, tmp, by.x='H1',by.y='AID',all.x=T)
#   dc <- merge(dc, tmp, by.x='H2',by.y='AID',all.x=T)
#   setnames(dc, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
#   dw <- merge(dw, tmp, by.x='H1',by.y='AID',all.x=T)
#   dw <- merge(dw, tmp, by.x='H2',by.y='AID',all.x=T)
#   setnames(dw, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
#   dchain <- merge(dchain, tmp, by.x='H1',by.y='AID',all.x=T)
#   dchain <- merge(dchain, tmp, by.x='H2',by.y='AID',all.x=T)
#   setnames(dchain, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
#   dnet <- merge(dnet, tmp, by.x='H1',by.y='AID',all.x=T)
#   dnet <- merge(dnet, tmp, by.x='H2',by.y='AID',all.x=T)
#   setnames(dnet, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
# 
#   # whether couples in dc
#   tmp <- unique(subset(dc,select = c('ID1','ID2')))
#   tmp[,ID.x:=paste0(ID1,'-',ID2)]
#   tmp[,ID.y:=paste0(ID2,'-',ID1)]
#   cat(nrow(couple[paste0(ID1,'-',ID2)%in% tmp$ID.x | paste0(ID1,'-',ID2)%in% tmp$ID.y ]), ' couples in networks')
# 
#   # to dt
#   dc19 <- as.data.table(dc19)
#   dw19 <- as.data.table(dw19)
#   dnet19 <- as.data.table(dnet19)
#   dchain19 <- as.data.table(dchain19)
#   dchain <- as.data.table(dchain)
# 
#   # add couple
#   tmp <- copy(couple)
#   setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
#   couple <- data.table(unique(rbind(couple, tmp)))
# 
#   dc <- merge(dc,couple,by=c('ID1','ID2'),all.x=T)
#   dc19 <- merge(dc19,couple,by=c('ID1','ID2'),all.x=T)
#   dw <- merge(dw,couple,by=c('ID1','ID2'),all.x=T)
#   dw19 <- merge(dw19,couple,by=c('ID1','ID2'),all.x=T)
#   dnet <- merge(dnet,couple,by=c('ID1','ID2'),all.x=T)
#   dnet19 <- merge(dnet19,couple,by=c('ID1','ID2'),all.x=T)
#   dchain <- merge(dchain,couple,by=c('ID1','ID2'),all.x=T)
#   dchain19 <- merge(dchain19,couple,by=c('ID1','ID2'),all.x=T)
# 
#   #
#   tmp_dc_coup <- unique(dc[CATEGORISATION=='close.and.adjacent.and.ancestry.cat'&COUPLE==1])
#   tmp_dc_coup[,length(N),by=c('ID1','ID2')]
#   tmp_dc_coup[ID1=='RK-A035039' & ID2=='RK-D061829']
#   tmp <- tmp_dc_coup[,list(N_EFF=max(N_EFF)),by=c('ID1','ID2')]
#   tmp_dc_coup <- merge(tmp_dc_coup,tmp, by=c('ID1','ID2','N_EFF'))
#   setkey(tmp_dc_coup, ID1, ID2, PTY_RUN, TYPE)
#   tmp_dc_coup <- tmp_dc_coup[,head(.SD, 4),by=c('ID1','ID2')]
#   tmp_dc_coup <- dcast(tmp_dc_coup,  ID1 + ID2 + COUPLE ~ TYPE, value.var= 'SCORE')
#   tmp_dc_coup[,ANALYSIS:=2021]
# 
#   tmp_dc19_coup <- unique(dc19[CATEGORISATION=='close.and.adjacent.and.ancestry.cat'&COUPLE==1])
#   # tmp_dc19_coup[,length(N),by=c('ID1','ID2')]
#   # tmp <- tmp_dc19_coup[,list(N_EFF=max(N_EFF)),by=c('ID1','ID2')]
#   # tmp_dc19_coup <- merge(tmp_dc19_coup,tmp, by=c('ID1','ID2','N_EFF'))
#   tmp_dc19_coup <- dcast(tmp_dc19_coup,  ID1 + ID2 + COUPLE ~ TYPE, value.var= 'SCORE')
#   tmp_dc19_coup[,ANALYSIS:=2019]
# 
#   tmp_dc_coup <- rbind(tmp_dc19_coup, tmp_dc_coup)
#   tmp <- tmp_dc_coup[!(ID1<ID2)]
#   setnames(tmp,c('ID1','ID2','12','21'),c('ID2','ID1','21','12'))
#   tmp_dc_coup <- rbind(tmp_dc_coup[(ID1<ID2)],tmp)
#   tmp_dc_coup <- melt(tmp_dc_coup,id.vars=c('ID1','ID2', 'COUPLE','ANALYSIS'))
# 
#   library(dplyr)
#   library(ggplot2)
#   df <- arrange(tmp_dc_coup, variable, desc(value))
# 
#   df$ID=paste0(df$ID1,'-',df$ID2)
#   df$ID=factor(df$ID, levels = unique(df$ID))
#   ggplot(df, aes(ID,value,fill=variable))+
#     geom_bar(position="stack", stat="identity")+
#     facet_grid(ANALYSIS~ .)+
#     labs(x='ID',y='scores',fill='')+ guides(fill = guide_legend(nrow = 1)) + theme(legend.position = "bottom",legend.direction = 'horizontal') +
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank())+
#     scale_x_discrete(limits = rev(levels(df$ID)))
#   ggsave(filename = '~/compare_couple.pdf',width = 10, height = 8)
#   tmp_dc_coup <- dcast(tmp_dc_coup,ID1+ID2+COUPLE+variable ~ANALYSIS ,value.var='value')
#   tmp1 <- tmp_dc_coup[is.na(`2019`) & !is.na(`2021`) ]
#   tmp2 <- tmp_dc_coup[!is.na(`2019`) & is.na(`2021`) ]
#   tmp3 <- tmp_dc_coup[!is.na(`2019`) & !is.na(`2021`) ]
#   nrow(unique(subset(tmp1,select=c('ID1','ID2'))))
#   nrow(unique(subset(tmp2,select=c('ID1','ID2'))))
#   nrow(unique(subset(tmp3,select=c('ID1','ID2'))))
#   # 2021 ONLY
#   tmp <- dcast(tmp1,ID1+ID2~variable,value.var='2021')
#   tmp[,LINKED:= `12`+`21` +complex.or.no.ancestry>0.6]
#   tmp[LINKED==1,DIRECTED:= (`12`+`21`)/(`12`+`21` +complex.or.no.ancestry)>0.6]
#   tmp[DIRECTED==1,LINK12:=`12`>`21`]
#   nrow(tmp[LINKED==1])/nrow(tmp)
#   nrow(tmp[DIRECTED==1])/nrow(tmp[LINKED==1])
#   # 2019 ONLY
#   tmp <- dcast(tmp2,ID1+ID2~variable,value.var='2019')
#   tmp[,LINKED:= `12`+`21` +complex.or.no.ancestry>0.6]
#   tmp[LINKED==1,DIRECTED:= (`12`+`21`)/(`12`+`21` +complex.or.no.ancestry)>0.6]
#   tmp[DIRECTED==1,LINK12:=`12`>`21`]
#   nrow(tmp[LINKED==1])/nrow(tmp)
#   nrow(tmp[DIRECTED==1])/nrow(tmp[LINKED==1])
#   tmp <- unique(subset(tmp2,select=c('ID1','ID2')))
#   tmp[ID1%in% pty.runs$UNIT_ID & ID2%in% pty.runs$UNIT_ID]
#   tmp[,ID.x:=paste0(ID1,'-',ID2)]
#   tmp[,ID.y:=paste0(ID2,'-',ID1)]
#   dw[paste0(ID1,'-',ID2) %in% tmp$ID.x]
#   dw[paste0(ID1,'-',ID2) %in% tmp$ID.y]
#   load('~/dcdwina_minreadnull.rda')
#   dca[,host.1:=gsub('CNTRL-','',host.1)]
#   dca[,host.2:=gsub('CNTRL-','',host.2)]
#   dcac <- merge(dca, aid, by.x='host.1',by.y='AID',all.x=T)
#   dcac <- merge(dcac, aid, by.x='host.2',by.y='AID',all.x=T)
#   dcac[paste0(PT_ID.x,'-',PT_ID.y) %in% tmp$ID.x & categorisation == 'close.and.adjacent.and.ancestry.cat']
#   dcac[paste0(PT_ID.x,'-',PT_ID.y) %in% tmp$ID.y & categorisation == 'close.and.adjacent.and.ancestry.cat']
# 
#   # 2019 and 2021
#   tmp21 <- dcast(tmp3,ID1+ID2~variable,value.var='2021')
#   tmp21[,LINKED_SCORE:=`12`+`21` +complex.or.no.ancestry]
#   tmp21[,LINKED:= LINKED_SCORE>0.6]
#   tmp21[LINKED==1,DIRECTED_SCORE:=(`12`+`21`)/LINKED_SCORE]
#   tmp21[LINKED==1,DIRECTED:= DIRECTED_SCORE>0.6]
#   tmp21[DIRECTED==1,LINK12_SCORE:=`12`/(`12`+`21`)]
#   tmp21[DIRECTED==1,LINK12:=`12`>`21`]
#   nrow(tmp21[LINKED==1])/nrow(tmp21)
#   nrow(tmp21[DIRECTED==1])/nrow(tmp21[LINKED==1])
# 
#   tmp19 <- dcast(tmp3,ID1+ID2~variable,value.var='2019')
#   tmp19[,LINKED_SCORE:=`12`+`21` +complex.or.no.ancestry]
#   tmp19[,LINKED:= LINKED_SCORE>0.6]
#   tmp19[LINKED==1,DIRECTED_SCORE:=(`12`+`21`)/LINKED_SCORE]
#   tmp19[LINKED==1,DIRECTED:= DIRECTED_SCORE>0.6]
#   tmp19[DIRECTED==1,LINK12_SCORE:=`12`/(`12`+`21`)]
#   tmp19[DIRECTED==1,LINK12:=`12`>`21`]
#   nrow(tmp19[LINKED==1])/nrow(tmp19)
#   nrow(tmp19[DIRECTED==1])/nrow(tmp19[LINKED==1])
# 
#   tmp <- merge(tmp19,tmp21,by=c('ID1','ID2'))
#   tmp[,LINKED_SAME:=LINKED.x==LINKED.y]
#   ggplot(tmp,aes(fill=factor(LINKED_SAME,c(TRUE,FALSE),c('yes','no'))))+
#     geom_boxplot(aes(x='2019',y=LINKED_SCORE.x)) +
#     geom_boxplot(aes(x='2021',y=LINKED_SCORE.y)) +
#     labs(x='analysis',y='linkage scores',fill='same linkage classification')+
#     guides(fill = guide_legend(nrow = 1)) +
#     theme(legend.position = "bottom",legend.direction = 'horizontal',legend.box = 'horizontal') +
#     theme_bw()
#   ggsave(filename = '~/compare_linkage_score_couple.pdf',width=6,height=4)
# 
# 
#   tmp[LINKED.x==T & LINKED.y==T,DIRECTED_SAME:=DIRECTED.x==DIRECTED.y]
#   tmp[DIRECTED.x==T & DIRECTED.y==T,LINK12_SAME:=LINK12.x==LINK12.y]
#   ggplot(tmp[LINKED.x==T & LINKED.y==T],aes(fill=factor(DIRECTED_SAME,c(TRUE,FALSE),c('yes','no'))))+
#     geom_boxplot(aes(x='2019',y=DIRECTED_SCORE.x)) +
#     geom_boxplot(aes(x='2021',y=DIRECTED_SCORE.y)) +
#     labs(x='analysis',y='direction scores',fill='same directed classification')+
#     guides(fill = guide_legend(nrow = 1)) +
#     theme(legend.position = "bottom",legend.direction = 'horizontal',legend.box = 'horizontal') +
#     theme_bw()
#   ggsave(filename = '~/compare_direction_score_couple.pdf',width=6,height=4)
# 
#   ggplot(tmp[DIRECTED.x==T & DIRECTED.y==T],aes(fill=factor(LINK12_SAME,c(TRUE,FALSE),c('yes','no'))))+
#     geom_boxplot(aes(x='2019',y=LINK12_SCORE.x)) +
#     geom_boxplot(aes(x='2021',y=LINK12_SCORE.y)) +
#     labs(x='analysis',y='1->2 scores',fill='same direction classification')+
#     guides(fill = guide_legend(nrow = 1)) +
#     theme(legend.position = "bottom",legend.direction = 'horizontal',legend.box = 'horizontal') +
#     theme_bw()
#   ggsave(filename = '~/compare_link12_score_couple.pdf',width=6,height=4)
# 
# 
#   dc_couple_1921 <- copy(tmp)
#   # tmp[LINKED.x==LINKED.y,`12.x`+`21.x` +complex.or.no.ancestry.x]
#   # tmp[LINKED.x==LINKED.y,`12.y`+`21.y` +complex.or.no.ancestry.y]
#   # tmp[LINKED.x!=LINKED.y,`12.x`+`21.x` +complex.or.no.ancestry.x]
#   # tmp[LINKED.x!=LINKED.y,`12.y`+`21.y` +complex.or.no.ancestry.y]
# 
#   require(Phyloscanner.R.utilities)
#   require(phyloscannerR)
#   require(tidyverse)
# 
#   dw19_couple <- dw19[COUPLE==1,c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE')]
#   aid[,AID:=as.character(AID)]
#   dw19_couple <- merge(dw19_couple, aid,by.x='ID1',by.y='PT_ID',all.x=T)
#   dw19_couple <- merge(dw19_couple, aid,by.x='ID2',by.y='PT_ID',all.x=T)
#   dw19_couple[,X.x:=NULL]
#   dw19_couple[,X.y:=NULL]
#   dw19_couple[,ID1:=NULL]
#   dw19_couple[,ID2:=NULL]
#   setnames( dw19_couple, c('AID.x','AID.y'), c('ID1','ID2'))
#   tmp <- dw19_couple[!(ID1<ID2)]
#   setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
#   dw19_couple <- unique(rbind(dw19_couple[ID1<ID2],tmp))
# 
#   dw_couple <- dw[COUPLE==1,c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN')]
#   setnames( dw_couple, c('H1','H2'), c('ID1','ID2'))
#   tmp <- dw_couple[,min(PTY_RUN), by=c('ID1','ID2')]
#   dw_couple <-  merge(dw_couple,tmp, by.x=c('ID1','ID2','PTY_RUN'),by.y=c('ID1','ID2','V1'))
#   tmp <- dw_couple[!(ID1<ID2)]
#   setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
#   dw_couple <- unique(rbind(dw_couple[ID1<ID2],tmp))
#   # dw_couple
# 
# 
#   # hosts	<- data.table(PT_ID=unique(c(tmp3$ID1,tmp3$ID2)))
#   # hosts <- merge(hosts, subset(aid,select=c('PT_ID','AID')),by='PT_ID',all.x=T)
#   # hosts <- as.character(unique(hosts$AID))
#   tmp <- dc_couple_1921[LINKED.x==LINKED.y & DIRECTED.x==DIRECTED.y & LINK12.x==LINK12.y,]
#   tmp <- merge(tmp, aid,by.x='ID1',by.y='PT_ID',all.x=T)
#   tmp <- merge(tmp, aid,by.x='ID2',by.y='PT_ID',all.x=T)
#   hosts_same <- unique(c(tmp$AID.x,tmp$AID.y))
#   tmp <-  dc_couple_1921[!(LINKED.x==LINKED.y & DIRECTED.x==DIRECTED.y & LINK12.x==LINK12.y),]
#   tmp <- merge(tmp, aid,by.x='ID1',by.y='PT_ID',all.x=T)
#   tmp <- merge(tmp, aid,by.x='ID2',by.y='PT_ID',all.x=T)
#   hosts_diff <- unique(c(tmp$AID.x,tmp$AID.y))
#   # hosts <- unique(c(tmp3$ID1,tmp3$ID2))
# 
#   # hosts %in% dw19_couple$ID1
#   inclusion <- "both"# "either"
#   setnames(dw19_couple, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
#            c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
#   setnames(dw_couple, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
#            c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
#   tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts_same,dwin=as_tibble(dw19_couple), inclusion = "both")
#   g1s <- tmp$graph
#   tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts_same,dwin=as_tibble(dw_couple), inclusion = "both")
#   g2s <- tmp$graph
#   library(ggpubr)
#   # arranges <- ggarrange(g1s, g2s, ncol = 2, nrow = 1)
#   # ggsave('~/couple_consistent_pairwise_plot.pdf',arranges,width = 16, height = 150, limitsize = F)
#   # # ggsave('~/couple_pairwise_plot.pdf',g2,width = 12, height = 300, limitsize = F)
#   #
#   tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts_diff,dwin=as_tibble(dw19_couple), inclusion = "both")
#   g1d <- tmp$graph
#   tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts_diff,dwin=as_tibble(dw_couple), inclusion = "both")
#   g2d <- tmp$graph
#   # library(ggpubr)
#   # arranged <- ggarrange(g1d, g2d, ncol = 2, nrow = 1)
#   # ggsave('~/couple_inconsistent_pairwise_plot.pdf',arranged,width = 16, height = 150, limitsize = F)
# 
#   # table( tmp$data$basic.topology )
#   indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read/' #30 min 100 max
#   infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
#   load(infile.networks)
#   dw <- merge(dw, aid, by.x='H1',by.y='AID',all.x=T)
#   dw <- merge(dw, aid, by.x='H2',by.y='AID',all.x=T)
#   setnames(dw, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
#   dw <- merge(dw,couple,by=c('ID1','ID2'),all.x=T)
# 
#   dw_couple <- dw[COUPLE==1,c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN')]
#   setnames( dw_couple, c('H1','H2'), c('ID1','ID2'))
#   tmp <- dw_couple[,min(PTY_RUN), by=c('ID1','ID2')]
#   dw_couple <-  merge(dw_couple,tmp, by.x=c('ID1','ID2','PTY_RUN'),by.y=c('ID1','ID2','V1'))
#   tmp <- dw_couple[!(ID1<ID2)]
#   setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
#   dw_couple <- unique(rbind(dw_couple[ID1<ID2],tmp))
# 
#   setnames(dw_couple, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
#            c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
# 
#   tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts_same,dwin=as_tibble(dw_couple), inclusion = "both")
#   g3s <- tmp$graph
#   tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts_diff,dwin=as_tibble(dw_couple), inclusion = "both")
#   g3d <- tmp$graph
# 
#   indir	<- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_null_min_read_100_max_read/' #no min 100 max
#   infile.networks <-  file.path(indir,'Rakai_phscnetworks.rda')
#   load(infile.networks)
#   dw <- merge(dw, aid, by.x='H1',by.y='AID',all.x=T)
#   dw <- merge(dw, aid, by.x='H2',by.y='AID',all.x=T)
#   setnames(dw, c('PT_ID.x','PT_ID.y'),c('ID1','ID2'))
#   dw <- merge(dw,couple,by=c('ID1','ID2'),all.x=T)
# 
#   dw_couple <- dw[COUPLE==1,c('H1','H2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE','PTY_RUN')]
#   setnames( dw_couple, c('H1','H2'), c('ID1','ID2'))
#   tmp <- dw_couple[,min(PTY_RUN), by=c('ID1','ID2')]
#   dw_couple <-  merge(dw_couple,tmp, by.x=c('ID1','ID2','PTY_RUN'),by.y=c('ID1','ID2','V1'))
#   tmp <- dw_couple[!(ID1<ID2)]
#   setnames(tmp,c('ID1','ID2'),c('ID2','ID1'))
#   dw_couple <- unique(rbind(dw_couple[ID1<ID2],tmp))
# 
#   setnames(dw_couple, c('ID1','ID2','TREE_ID','BASIC_CLASSIFICATION','PATRISTIC_DISTANCE'),
#            c('host.1', 'host.2', 'tree.id', 'basic.classification', 'patristic.distance') )
# 
#   tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts_same,dwin=as_tibble(dw_couple), inclusion = "both")
#   g4s <- tmp$graph
#   tmp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts_diff,dwin=as_tibble(dw_couple), inclusion = "both")
#   g4d <- tmp$graph
# 
# 
#   library(ggpubr)
#   arranges <- ggarrange(g1s + ggtitle('2019 analysis \n '),
#                         g3s + ggtitle('2021 analysis, \n with minreadhost and maxreadhost'),
#                         g4s + ggtitle('2021 analysis, \n without minreadhost and with maxreadhost'),
#                         g2s + ggtitle('2021 analysis, \n without minreadhost and maxreadhost'), ncol = 4, nrow = 1)
#   ggsave('~/couple_consistent_pairwise_plot.pdf',arranges,width = 20, height = 150, limitsize = F)
# 
#   arranged <- ggarrange(g1d + ggtitle('2019 analysis \n '),
#                         g3d + ggtitle('2021 analysis, \n with minreadhost and maxreadhost'),
#                         g4d + ggtitle('2021 analysis, \n without minreadhost and with maxreadhost'),
#                         g2d + ggtitle('2021 analysis, \n without minreadhost and maxreadhost'), ncol = 4, nrow = 1)
# 
#   ggsave('~/couple_inconsistent_pairwise_plot.pdf',arranged,width = 20, height = 150, limitsize = F)
# 
#   tmp_same <-  dc_couple_1921[LINKED.x==LINKED.y & DIRECTED.x==DIRECTED.y & LINK12.x==LINK12.y,]
#   tmp_same <-  merge(tmp_same, aid, by.x= 'ID1', by.y='PT_ID',all.x=T)
#   tmp_same <-  merge(tmp_same, aid, by.x= 'ID2', by.y='PT_ID',all.x=T)
#   tmp <- tmp_same[,c('AID.x','AID.y','LINKED.x','LINKED.y','LINKED_SCORE.x','LINKED_SCORE.y',
#                      'DIRECTED.x','DIRECTED.y','DIRECTED_SCORE.x','DIRECTED_SCORE.y',
#                      'LINK12.x','LINK12.y','LINK12_SCORE.x','LINK12_SCORE.y')]
#   write.csv(tmp,'~/couple_consistent_pairwise.csv')
# 
#   tmp_diff <-  dc_couple_1921[!(LINKED.x==LINKED.y & DIRECTED.x==DIRECTED.y & LINK12.x==LINK12.y),]
#   tmp_diff <-  merge(tmp_diff, aid, by.x= 'ID1', by.y='PT_ID',all.x=T)
#   tmp_diff <-  merge(tmp_diff, aid, by.x= 'ID2', by.y='PT_ID',all.x=T)
#   tmp <- tmp_diff[,c('AID.x','AID.y','LINKED.x','LINKED.y','LINKED_SCORE.x','LINKED_SCORE.y',
#                      'DIRECTED.x','DIRECTED.y','DIRECTED_SCORE.x','DIRECTED_SCORE.y',
#                      'LINK12.x','LINK12.y','LINK12_SCORE.x','LINK12_SCORE.y')]
#   write.csv(tmp,'~/couple_inconsistent_pairwise.csv')
# 
#   setkey(tmp, AID.x,AID.y)
#   setkey(tmp_diff, AID.x,AID.y)
#   tmp_diff[c(4,7)]
#   # couple plots
#   library(ggplot2)
#   load('~/dcdwina_minreadnull.rda')
#   dwina[,host.1:=gsub('CNTRL-','',host.1)]
#   dwina[,host.2:=gsub('CNTRL-','',host.2)]
#   dwinac <- merge(dwina , aid, by.x='host.1',by.y='AID',all.x=T)
#   dwinac <- merge(dwinac , aid, by.x='host.2',by.y='AID',all.x=T)
#   dwinac <- merge(dwinac,couple,by.x=c('PT_ID.x','PT_ID.y'),by.y=c('ID1','ID2'),all.x=T)
#   tmp <- unique(subset(dwinac[COUPLE==1],select=c('host.1','host.2','PTY_RUN','tree.id', 'patristic.distance')))
# 
#   ggplot(tmp, aes(patristic.distance))+
#     # geom_density()+theme_bw()+labs(x='patristic distance', y='count')+ coord_cartesian(xlim = c(0,0.05))
#     geom_histogram(binwidth = 0.002)+theme_bw()+labs(x='patristic distance', y='count')
#   # coord_cartesian(ylim = c(0,1e3))
#   ggsave(filename = '~/pdistance_couple.pdf',width = 6, height = 4)
# 
#   ggplot(tmp, aes(patristic.distance))+
#     # geom_density()+theme_bw()+labs(x='patristic distance', y='count')+ coord_cartesian(xlim = c(0,0.05))
#     geom_histogram(binwidth = 0.002)+theme_bw()+labs(x='patristic distance', y='count') +coord_cartesian(ylim = c(0,1e3),xlim=c(0,0.1))
#   # coord_cartesian(ylim = c(0,1e3))
#   ggsave(filename = '~/pdistance_couple_zoomin.pdf',width = 6, height = 4)
# 
#   # pd <- tmp[patristic.distance>0.05 & patristic.distance<0.15,]$patristic.distance
#   # library(MASS)
#   # fit <- fitdistr(pd, "normal")
#   # qnorm(0.05,fit$estimate[1],fit$estimate[2])
# 
#   # make couple plots
#   dw_couple <- dwinac[COUPLE==1]
#   tmp <- dw_couple[,list(M=median(patristic.distance),
#                          CL=quantile(patristic.distance,probs = 0.025),
#                          CU=quantile(patristic.distance,probs = 0.975)),
#                    by=c('host.2','host.1')]
#   setkey(tmp,M)
#   tmp[,ID:=seq_len(nrow(tmp))]
#   ggplot(tmp,aes(ID,M))+
#     geom_point()+
#     theme_bw()+
#     geom_errorbar(aes(ymin=CL,ymax=CU))+
#     labs(x='couples',y='patristic distance')+
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank())
#   # +
#   #   coord_cartesian(ylim=c(0,0.2))
#   ggsave('~/patdist_vs_couple.pdf',width = 10, height = 4)
# 
#   tmpid <- subset(tmp,select = c('host.1','host.2','ID'))
# 
#   dcac <- merge(dcac,couple,by.x=c('PT_ID.x','PT_ID.y'),by.y=c('ID1','ID2'),all.x=T)
#   dc_couple <- dcac[COUPLE==1 & categorisation == 'close.and.adjacent.and.ancestry.cat']
#   tmp <- dc_couple[,list(n.eff=max(n.eff)),by=c('host.1','host.2')]
#   tmp <- merge(tmp,dc_couple, by=c('host.1','host.2','n.eff'))
#   setkey(tmp, host.1, host.2, PTY_RUN, type)
#   tmp <- tmp[,head(.SD, 1),by=c('host.1','host.2','type')]
#   tmp <- merge(tmp, tmpid, by=c('host.1','host.2'))
#   # tmp2 = tmp[,length(n.eff),by=c('host.1','host.2')]
#   # tmp2 = tmp[,as.integer(sum(score)),by=c('host.1','host.2')]
#   ggplot(tmp,aes(ID,k.eff,fill=type))+
#     theme_bw()+
#     geom_bar(position="stack", stat="identity")+
#     labs(x='couples',y='number of windows \n supporting subgraph topology')+
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank())
#   ggsave('~/topo_num_vs_couple.pdf',width = 10, height = 4)
# 
#   ggplot(tmp,aes(ID,score,fill=type))+
#     theme_bw()+
#     geom_bar(position="stack", stat="identity")+
#     labs(x='couples',y='proportion of windows \n supporting subgraph topology')+
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank())
#   ggsave('~/topo_prop_vs_couple.pdf',width = 10, height = 4)
# 
#   tmp <- dc_couple[,list(n.eff=max(n.eff)),by=c('host.1','host.2')]
#   tmp <- merge(tmp,dc_couple, by=c('host.1','host.2','n.eff'))
#   setkey(tmp, host.1, host.2, PTY_RUN, type)
#   tmp <- tmp[,head(.SD, 1),by=c('host.1','host.2','type')]
#   tmp <- merge(tmp, tmpid, by=c('host.1','host.2'))
#   tmp <- dcast(tmp, host.1+host.2 + ID~type, value.var = 'score')
#   tmp[,score_L:=`12`+ `21`+ complex.or.no.ancestry]
#   tmp[,score_D:=(`12`+ `21`)/score_L]
#   ggplot(tmp,aes(ID,score_L))+
#     geom_point()+
#     theme_bw()+
#     scale_y_continuous(labels = scales::percent,limits = c(0,1))+
#     # geom_errorbar(aes(ymin=CL,ymax=CU))+
#     labs(x='couples',y='linkage score')+
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank())
#   ggsave('~/lscore_vs_couple.pdf',width = 10, height = 4)
# 
#   ggplot(tmp,aes(ID,score_D))+
#     geom_point()+
#     theme_bw()+
#     # geom_errorbar(aes(ymin=CL,ymax=CU))+
#     labs(x='couples',y='linkage score')+
#     scale_y_continuous(labels = scales::percent,limits = c(0,1))+
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank())
#   ggsave('~/dscore_vs_couple.pdf',width = 10, height = 4)
# 
# 
# }
# 
