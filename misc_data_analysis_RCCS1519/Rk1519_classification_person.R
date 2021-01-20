library(data.table)
library(ggplot2)
library(seqinr)
library(ggpubr)
library(tidyverse)
library(igraph)
library(dplyr)

is.match <- function(x,y){
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

if(1)
{
  args_dir <- list()
  args_dir[['data_dir']] <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
  args_dir[['out_dir']] <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY'
  args_dir[['window_size']] <- 500
  args_dir[['batch_size']] <- 100
  args_dir[['script_dir']] <- '~/phyloscanner/data_analysis_RCCS1519'
  args_dir[['job_tag']] <- 'windowsize500_batchsize100'
  args_dir[['windown']] <- 99
  args_dir[['length_cutoff']] <- 250
  args_dir[['method']] <- 2
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
  stopifnot(args_line[[13]]=='-windown')
  stopifnot(args_line[[15]]=='-length_cutoff')
  stopifnot(args_line[[17]]=='-method')
  args_dir <- list()
  args_dir[['data_dir']] <- args_line[[2]]
  args_dir[['out_dir']] <- args_line[[4]]
  args_dir[['window_size']] <- as.integer(args_line[[6]])
  args_dir[['batch_size']] <- as.integer(args_line[[8]])
  args_dir[['script_dir']] <- args_line[[10]]
  args_dir[['job_tag']] <- args_line[[12]]
  args_dir[['windown']] <- as.integer(args_line[[14]])
  args_dir[['length_cutoff']] <- as.integer(args_line[[16]])
  args_dir[['method']] <- as.integer(args_line[[18]])
} 



# dir
setwd( paste0(args_dir[['out_dir']],'_',args_dir[['job_tag']] ))

# load
load(paste0('support_overall_method', args_dir[['method']],'_sequence_level_processed_v4.rda'))
nrow(df)

# enough overlapping windows with high depths
df = df[DATA >=4,]
nrow(df)

# take close pairs
df =df[PROP >= 0.8]
nrow(df)

# add study id
id.dt <- data.table(read.csv(file.path(args_dir[['data_dir']],'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')))
id.dt <- subset(id.dt,select = c("pt_id","sex","pangea_id"))
id.dt[,pangea_id:=paste0('RCCS_',pangea_id)]
tmp <- data.table(read.csv(file.path(args_dir[['data_dir']],'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')))
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
save(close_pairs,file='close_pairs.rda')

load('close_pairs.rda')
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
mapping_rccs <- readRDS('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/200422_PANGEA2_RCCS_mapped_samples.rds')
dconsensus <- subset(mapping_rccs,select=c('PANGEA_ID','CONSENSUS','F'))
dconsensus <- unique(dconsensus[CONSENSUS!="",])
dconsensus[,PANGEA_ID:=paste0('RCCS_',PANGEA_ID)]
mapping_mrc <- readRDS('/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC/200422_PANGEA2_MRCUVRI_mapped_samples.rds')
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

save(df,file = 'clusters.rda')
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
  

ggsave('cluster_size_plot.pdf',width = 10, height = 4)