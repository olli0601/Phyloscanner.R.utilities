# The phyloscanner run for TSI estimates - create potential networks

# Preamble
# The set of scripts aims to run phyloscanner for TSI estimates.
# We change the sliding windows from 25 to 10 and
# update the options for buiding and analysing trees.
# This script aims to create potential networks from close pairs.

# Load the required packages
library(data.table)
library(seqinr)
library(tidyverse)
library(igraph)
library(dplyr)

#
# Define input arguments that can be changed by users
#
option_list <- list(
  optparse::make_option(
    c("-v", "--verbose"),
    action = "store_true",
    default = FALSE,
    help = "Print extra output [default]",
    dest = "verbose"
  ),
  optparse::make_option(
    "--seed",
    type = "integer",
    default = 42L,
    help = "Random number seed [default %default]",
    dest = "seed"
  ),
  optparse::make_option(
    "--window_cutoff",
    type = "numeric",
    default = 0.5,
    help = "Cutoff of proportion of windows indicating close relationship [default %default]",
    dest = "window_cutoff"
  ),
  optparse::make_option(
    "--n_control",
    type = "integer",
    default = 5,
    help = "Number of controls added [default %default]",
    dest = "n_control"
  ),
  optparse::make_option(
    "--save_data",
    action = "store_true",
    default = TRUE,
    help = "Save data [default]",
    dest = 'if_save_data'
  ),
  optparse::make_option(
    "--pkg_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to package directory, used as long we don t build an R package [default]",
    dest = 'prj.dir'
  ),
  optparse::make_option(
    "--out_dir_base",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where all output is stored [default]",
    dest = 'out.dir'
  ),
  optparse::make_option(
    "--date",
    type = 'character',
    default = Sys.Date(),
    metavar = '"YYYY-MM-DD"',
    help = 'As of date to extract data from.  Defaults to today.',
    dest = 'date'
  )
  
)

args <-
  optparse::parse_args(optparse::OptionParser(option_list = option_list))

#
# test
#
args <- list(
  verbose = T,
  seed = 42,
  sliding_width = 10L,
  window_cutoff = 0.5,
  n_control = 5,
  if_save_data = T,
  date = Sys.Date(),
  out.dir = NA,
  prj.dir = NA
)

#
# use manually specified directories when args$out.dir is NA
#
tmp <- Sys.info()
if (tmp["user"] == "xx4515" &
    grepl("hpc.ic.ac.uk", tmp["nodename"]))
{
  if (is.na(args$out.dir))
  {
    args$out.dir <-
      "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/"
  }
  if (is.na(args$prj.dir))
  {
    args$prj.dir <-
      "~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/"
  }
}

# if prj.dir and out.dir are not manually set, default to here()
if (is.na(args$prj.dir))
{
  args$prj.dir <- here::here()
  args$out.dir <- here::here()
}

#
# add constants that should not be changed by the user
#
dir.data <-
  '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
dir.net <-
  file.path(args$out.dir, "200929_SIMILARTY_windowsize500_batchsize100")
infile.couple <-
  file.path(args$out.dir,
            'RakaiPangeaMetaData_v2.rda')
infile.ind.rccs <-
  file.path(dir.data,
            'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
infile.ind.mrc <-
  file.path(dir.data,
            'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
infile.support <-
  file.path(dir.net ,
            paste0('support_overall_method2_sequence_level_processed_v4.rda'))
infile.close.paris <-
  file.path(dir.net,
            paste0(
              'close_pairs_cutoffprop',
              gsub('\\.', '', args$window_cutoff),
              '.rda'
            ))
infile.ddist <- file.path(dir.net, 'ddist.rda')
infile.pty.runs <-
  file.path(args$out.dir, '210120_RCCSUVRI_phscinput_samples.rds')


# Define functions just used in this script alone
take_bridging_individual <- function(df_graph, df) {
  bridging = list()
  for (k in 1:max(df$MEMBERSHIP)) {
    tmp1 = df[MEMBERSHIP == k,]$PT_ID
    tmp2 = df_graph[pt_id1 %in% tmp1 &
                      !(pt_id2 %in% tmp1) & SIMILARITY >= 0.975,]
    tmp3 = df_graph[pt_id2 %in% tmp1 &
                      !(pt_id1 %in% tmp1) & SIMILARITY >= 0.975,]
    bridging[[k]] = unique(c(as.character(tmp2$pt_id1),
                             as.character(tmp3$pt_id2)))
  }
  return(bridging)
}
break_large_clusters <- function(idlarge, dflarge, ddist) {
  ans <- data.table()
  for (i in idlarge) {
    idlarge_ind <- dflarge[MEMBERSHIP == i]$PT_ID
    maxsize <- unique(dflarge[MEMBERSHIP == i,]$CLU_SIZE)
    ddist_idlarge_ind <-
      ddist[pt_id1 %in% idlarge_ind & pt_id2 %in% idlarge_ind]
    chain_idlarge_ind <-
      graph.data.frame(ddist_idlarge_ind,
                       directed = FALSE,
                       vertices = NULL)
    E(chain_idlarge_ind)$weight <- ddist_idlarge_ind$SIMILARITY
    chain_idlarge_ind <- simplify(chain_idlarge_ind)
    E(chain_idlarge_ind)$weight = ifelse(
      E(chain_idlarge_ind)$weight > 1,
      E(chain_idlarge_ind)$weight / 2,
      E(chain_idlarge_ind)$weight
    )
    tmp = clusters(chain_idlarge_ind, mode = 'weak')
    stopifnot(all(tmp$membership == 1))
    comm_idlarge_ind <- cluster_louvain(chain_idlarge_ind)
    df_idlarge_ind = data.table(PT_ID = comm_idlarge_ind$names,
                                MEMBERSHIP = comm_idlarge_ind$membership)
    df_graph_idlarge_ind = as_long_data_frame(chain_idlarge_ind)
    df_graph_idlarge_ind = data.table(df_graph_idlarge_ind[, 3:5])
    colnames(df_graph_idlarge_ind) = c('SIMILARITY', 'pt_id1', 'pt_id2')
    bridging <-
      take_bridging_individual(df_graph_idlarge_ind, df_idlarge_ind)
    bridging_all <- unique(do.call(c, bridging))
    if (length(bridging_all) != 0) {
      tmpidclu <- max(comm_idlarge_ind$membership) + 1
      df_idlarge_ind = rbind(df_idlarge_ind,
                             data.table(PT_ID = bridging_all, MEMBERSHIP = tmpidclu))
    }
    tmp = df_idlarge_ind[, list(CLU_SIZE = length(PT_ID)), by = 'MEMBERSHIP']
    df_idlarge_ind = merge(df_idlarge_ind, tmp, by = 'MEMBERSHIP')
    df_idlarge_ind[, IDC := i]
    ans <- rbind(ans, df_idlarge_ind)
  }
  tmp <- unique(subset(ans, select = c('MEMBERSHIP', 'IDC')))
  tmp[, MEMBERSHIP2 := seq_len(nrow(tmp))]
  ans <- merge(ans, tmp, by = c('MEMBERSHIP', 'IDC'))
  ans[, MEMBERSHIP := NULL]
  ans[, IDC := NULL]
  setnames(ans, 'MEMBERSHIP2', 'MEMBERSHIP')
  return(ans)
}

# Set seed
set.seed(args$seed)

cat('----------- Generate clusters ----------- \n')
load(infile.close.paris)
chains <-
  graph.data.frame(close_pairs, directed = FALSE, vertices = NULL)
rtc	<-
  data.table(ID = V(chains)$name,
             CLU = clusters(chains, mode = "weak")$membership)
tmp	<- rtc[, list(CLU_SIZE = length(ID)), by = 'CLU']
setkey(tmp, CLU_SIZE)
if (args$verbose) {
  cat('---------- Largest cluster sizes ----------\n')
  print(tail(tmp))
}
tmp[, IDCLU := rev(seq_len(nrow(tmp)))]
rtc	<- subset(merge(rtc, tmp, by = 'CLU'))
rtc[, CLU := NULL]
setkey(rtc, IDCLU)

cat('----------- Break the largest clusters into sub-clusters of size <=55 ----------- \n')
id = rtc[IDCLU == 1,]$ID
chains <- subset(close_pairs, select = c(pt_id1, pt_id2))
chains <- chains[pt_id1 %in% id & pt_id2 %in% id]
load(infile.ddist)
ddist = copy(ddist_copy)
ddist[SIMILARITY >= 0.975, SIMILARITY := 1]
ddist <-
  ddist[, list(SIMILARITY = mean(SIMILARITY)), by = c('pt_id1', 'pt_id2')]
chains_tmp <-
  graph.data.frame(ddist, directed = FALSE, vertices = NULL)
E(chains_tmp)$weight = ddist$SIMILARITY
chains_tmp <- simplify(chains_tmp)
E(chains_tmp)$weight = ifelse(E(chains_tmp)$weight > 1,
                              E(chains_tmp)$weight / 2,
                              E(chains_tmp)$weight)
comm <- cluster_louvain(chains_tmp)
df = data.table(PT_ID = comm$names, MEMBERSHIP = comm$membership)
df_graph = as_long_data_frame(chains_tmp)
df_graph = data.table(df_graph[, 3:5])
colnames(df_graph) = c('SIMILARITY', 'pt_id1', 'pt_id2')
bridging <- take_bridging_individual(df_graph, df)
bridging_all <- unique(do.call(c, bridging))
id_clu <- max(df$MEMBERSHIP)
tmp <- data.table(PT_ID = bridging_all, MEMBERSHIP = id_clu + 1)
df <- rbind(df, tmp)
df = rbind(df, data.table(PT_ID = bridging_all, MEMBERSHIP = id_clu))
tmp = df[, list(CLU_SIZE = length(PT_ID)), by = 'MEMBERSHIP']
df = merge(df, tmp, by = 'MEMBERSHIP')

# Repeat tricks
dflarge <- df[CLU_SIZE > 55]
df <- df[CLU_SIZE <= 55]
tmp <-
  data.table(MEMBERSHIP2 = seq_len(length(unique(df$MEMBERSHIP))),
             MEMBERSHIP = unique(df$MEMBERSHIP))
df <- merge(df, tmp, by = 'MEMBERSHIP')
df[, MEMBERSHIP := NULL]
setnames(df, 'MEMBERSHIP2', 'MEMBERSHIP')
idlarge <- unique(dflarge$MEMBERSHIP)
breaked_dflarge <- break_large_clusters(idlarge, dflarge, ddist)
while (T) {
  id_clu <- max(df$MEMBERSHIP)
  tmp <- breaked_dflarge[CLU_SIZE <= 55]
  tmp2 <- unique(subset(tmp, select = c('CLU_SIZE', 'MEMBERSHIP')))
  setkey(tmp2, MEMBERSHIP)
  tmp2[, MEMBERSHIP2 := seq_len(nrow(tmp2))]
  tmp <- merge(tmp, tmp2, by = c('CLU_SIZE', 'MEMBERSHIP'))
  tmp[, MEMBERSHIP := NULL]
  setnames(tmp, 'MEMBERSHIP2', 'MEMBERSHIP')
  tmp[, MEMBERSHIP := MEMBERSHIP + id_clu]
  df <- rbind(df, tmp)
  dflarge <- breaked_dflarge[CLU_SIZE > 55]
  idlarge <- unique(dflarge$MEMBERSHIP)
  if (length(idlarge) == 0)
    break
  breaked_dflarge <- break_large_clusters(idlarge, dflarge, ddist)
}

# Merge
setnames(df, c('PT_ID', 'MEMBERSHIP'), c('ID', 'IDCLU'))
rtc = rtc[IDCLU != 1,]
rtc[, IDCLU := IDCLU + max(df$IDCLU) - 1]
df = rbind(rtc, df)

# Save
if (args$if_save_data) {
  cat('Write potential transmission networks to ',
      file.path(dir.net,
                paste0(
                  'clusters_cutoff',
                  gsub('\\.', '', args$window_cutoff),
                  '.rda'
                )),
      '...\n')
  save(df, file = file.path(dir.net,
                            paste0(
                              'clusters_cutoff',
                              gsub('\\.', '', args$window_cutoff),
                              '.rda'
                            )))
}

# Order IDCLU
tmp <- unique(subset(df, select = c('CLU_SIZE', 'IDCLU')))
setkey(tmp, CLU_SIZE)
tmp[, IDCLU2 := seq_len(nrow(tmp))]
tmp[, CLU_SIZE := NULL]
rtc	<-
  subset(merge(df, tmp, by = 'IDCLU'), select = c('ID', 'IDCLU2', 'CLU_SIZE'))
setnames(rtc, 'IDCLU2', 'IDCLU')
setkey(rtc, IDCLU)

cat('----------- Add controls to each of the clusters ----------- \n')
pty.runs <- readRDS(infile.pty.runs)
#	Add couples
load(infile.couple)
rp <-
  data.table(unique(subset(
    coupdat, select = c('male.RCCS_studyid', 'female.RCCS_studyid')
  )))
setnames(rp,
         c('male.RCCS_studyid', 'female.RCCS_studyid'),
         c('MALE_RID', 'FEMALE_RID'))
rp <- subset(rp, MALE_RID != FEMALE_RID)
rp[, FEMALE_RID := paste0('RK-', FEMALE_RID)]
rp[, MALE_RID := paste0('RK-', MALE_RID)]
tmp <- sort(unique(as.character(pty.runs[['UNIT_ID']])))
rp <-
  unique(subset(
    rp,
    FEMALE_RID %in% tmp &
      MALE_RID %in% tmp,
    select = c(FEMALE_RID, MALE_RID)
  ))
setnames(rp, 'FEMALE_RID', 'ID')
rp <-
  merge(rp, subset(rtc, select = c(ID, IDCLU)), all.x = 1, by = 'ID')
setnames(rp,
         c('ID', 'IDCLU', 'MALE_RID'),
         c('FEMALE_RID', 'FEMALE_IDCLU', 'ID'))
rp <-
  merge(rp, subset(rtc, select = c(ID, IDCLU)), all.x = 1, by = 'ID')
setnames(rp, c('ID', 'IDCLU'), c('MALE_RID', 'MALE_IDCLU'))
rp[, COUP_ID := seq_len(nrow(rp))]

#	Reassign partners that are not in the same network as their partner
tmp <- subset(rp, !is.na(FEMALE_IDCLU) & !is.na(MALE_IDCLU))
tmp2 <-
  tmp[, list(NOT_IN_SAME = !any(FEMALE_IDCLU == MALE_IDCLU)), by = c('MALE_RID', 'FEMALE_RID')]
tmp2 <- subset(tmp2, NOT_IN_SAME)
tmp2 <- merge(tmp2, rp, by = c('MALE_RID', 'FEMALE_RID'))
tmp <- rp[, which(COUP_ID %in% tmp2$COUP_ID)]
set(rp, tmp, 'MALE_IDCLU', rp[tmp, FEMALE_IDCLU])
set(tmp2, NULL, 'FEMALE_IDCLU', tmp2[, MALE_IDCLU])
set(tmp2, NULL, 'NOT_IN_SAME', NULL)
rp <- rbind(rp, tmp2)

#	Assign partners that are in no network to the same potential transmission network as their partner
tmp <- rp[, which(is.na(FEMALE_IDCLU) & !is.na(MALE_IDCLU))]
set(rp, tmp, 'FEMALE_IDCLU', rp[tmp, MALE_IDCLU])
tmp <- rp[, which(!is.na(FEMALE_IDCLU) & is.na(MALE_IDCLU))]
set(rp, tmp, 'MALE_IDCLU', rp[tmp, FEMALE_IDCLU])

#	Make a new potential transmission network for couples that are in no network so far
tmp <- rp[, which(is.na(FEMALE_IDCLU) & is.na(MALE_IDCLU))]
set(rp, tmp, c('FEMALE_IDCLU', 'MALE_IDCLU'), rtc[, max(IDCLU)] + seq_along(tmp))

#	Calculate cluster size
setnames(rp, c('MALE_IDCLU'), c('IDCLU'))
rp <-
  subset(melt(
    rp,
    id.vars = c('IDCLU'),
    measure.vars = c('MALE_RID', 'FEMALE_RID'),
    value.name = 'ID'
  ),
  select = c(ID, IDCLU))
rtc <- unique(rbind(subset(rtc, select = c('ID', 'IDCLU')), rp))
tmp <- rtc[, list(CLU_SIZE = length(ID)), by = 'IDCLU']
setkey(tmp, CLU_SIZE)
tmp[, IDCLU2 := seq_len(nrow(tmp))]
rtc <-
  subset(merge(rtc, tmp, by = 'IDCLU'),
         select = c('ID', 'IDCLU2', 'CLU_SIZE'))
setnames(rtc, 'IDCLU2', 'IDCLU')
setkey(rtc, IDCLU)

# Merge small runs
tn	<- 8
tmp	<- unique(subset(rtc, select = c(IDCLU, CLU_SIZE)))
tmp[, PTY_RUN := IDCLU]
tmp[, PTY_SIZE := CLU_SIZE]
sum_reset_at <- function(thresh) {
  function(x) {
    accumulate(x, ~ if_else(.x + .y > thresh, .y, .x + .y))
  }
}
tmp <- tmp %>% mutate(PTY_SIZE2 = sum_reset_at(tn)(PTY_SIZE))
tmp <- as.data.table(tmp)
tmp2 <- which(diff(tmp$PTY_SIZE2) < 0) + 1
tmp[, PTY_RUN2 := 0]
tmp[c(tmp2, max(tmp2):nrow(tmp)), PTY_RUN2 := 1]
tmp[, PTY_RUN2 := cumsum(PTY_RUN2)]
stopifnot(sum(tmp[, cumsum(PTY_SIZE), by = 'PTY_RUN2']$V1 != tmp$PTY_SIZE2) ==
            0)
tmp[, PTY_SIZE2 := max(PTY_SIZE2), by = 'PTY_RUN2']
tmp[, PTY_RUN2 := PTY_RUN2 + 1]
set(tmp, NULL, c('PTY_RUN', 'PTY_SIZE'), NULL)
setnames(tmp, c('PTY_RUN2', 'PTY_SIZE2'), c('PTY_RUN', 'PTY_SIZE'))
rtc <- merge(rtc, tmp, by = c('IDCLU', 'CLU_SIZE'))

# Find n_control closest individuals to all the individuals in each cluster
id.dt <- data.table(read.csv(infile.ind.rccs))
id.dt <- subset(id.dt, select = c("pt_id", "sex", "pangea_id"))
id.dt[, pangea_id := paste0('RCCS_', pangea_id)]
tmp <- data.table(read.csv(infile.ind.mrc))
tmp <- subset(tmp, select = c("pt_id", "sex", "pangea_id"))
tmp[, pangea_id := paste0('MRCUVRI_', pangea_id)]
id.dt <- rbind(id.dt, tmp)
id.dt <- unique(id.dt)
infile.average.distance.per.pair <-
  file.path(dir.net,
            'average_distance_perpair.rds')
if (file.exists(infile.average.distance.per.pair))
{
  ddist <- readRDS(infile.average.distance.per.pair)
} else{
  infile.average.distance.per.pair.per.window <-
    file.path(dir.net,
              'distance_overall_method2_sequence_level_v4.rda')
  load(infile.average.distance.per.pair.per.window)
  ddist <- subset(distance, select = c('TAXA1', 'TAXA2'))
  ddist$DIST <- rowMeans(distance[, 3:ncol(distance)], na.rm = TRUE)
  setnames(id.dt, colnames(id.dt), paste0(colnames(id.dt), '1'))
  ddist <-
    merge(ddist,
          id.dt,
          by.x = 'TAXA1',
          by.y = 'PANGEA_ID1',
          all.x = TRUE)
  setnames(id.dt, colnames(id.dt), gsub('1', '2', colnames(id.dt)))
  ddist <-
    merge(ddist,
          id.dt,
          by.x = 'TAXA2',
          by.y = 'PANGEA_ID2',
          all.x = TRUE)
  setnames(id.dt, colnames(id.dt), gsub('2', '', colnames(id.dt)))
  ddist <-
    ddist[, list(DIST = mean(DIST)), by = c('UNIT_ID1', 'UNIT_ID2')]
  saveRDS(ddist, file = infile.average.distance.per.pair)
}
dcl <- rtc[, {
  tmp <- ddist[UNIT_ID1 %in% ID, c('UNIT_ID2', 'DIST')]
  tmp2 <- ddist[UNIT_ID2 %in% ID, c('UNIT_ID1', 'DIST')]
  tmp <- rbind(tmp, tmp2, use.names = F)
  colnames(tmp) <- c('UNIT_ID', 'DIST')
  tmp <- tmp[!UNIT_ID %in% ID,]
  tmp <- tmp[, list(mean(DIST)), by = 'UNIT_ID']
  tmp <- tmp[order(V1,decreasing = T),]
  list(variable = c(paste0('ID_CLOSE',1:args$n_control), paste0('DISTANCE', 1:args$n_control)),
       value = c(as.character(tmp$UNIT_ID[1:5]), tmp$V1[1:5]))
}, by = c('CLU_SIZE', 'IDCLU')]
dcl <- dcast(dcl, CLU_SIZE + IDCLU ~ variable, value.var='value')

if(args$if_save_data){
  infile.clostest.per.cluster <- file.path(dir.net,paste0(args$n_control, 'closest_individuals_percluster.rds'))
  saveRDS(dcl, file = infile.clostest.per.cluster)  
}

tmp <-
  subset(dcl,
         select = c('CLU_SIZE', 'IDCLU', paste0('ID_CLOSE', 1:args$n_control)))
tmp <-
  melt(
    tmp,
    id.vars = c('CLU_SIZE', 'IDCLU'),
    variable.name = 'CLOSEST',
    value.name = 'ID'
  )
tmp[, CLOSEST := as.integer(gsub('ID_CLOSE', '', CLOSEST))]
tmp2 <- unique(subset(rtc, select = c('IDCLU', 'PTY_RUN')))
set(rtc, NULL , c('PTY_SIZE', 'PTY_RUN'), NULL)
set(tmp, NULL, 'CLOSEST', NULL)
tmp[, ID_TYPE := 'control']
rtc[, ID_TYPE := 'target']
rtc <- rbind(rtc, tmp)
rtc <- unique(rtc)
rtc[, CLU_SIZE := length(ID), by = 'IDCLU']
rtc <- merge(rtc, tmp2, by = 'IDCLU')
rtc[, PTY_SIZE := length(ID), by = 'PTY_RUN']

if (args$verbose) {
  cat(
    max(rtc$IDCLU),
    ' clusters of ',
    length(unique(rtc$ID)),
    ' individuals the ten largest cluster sizes are ',
    unique(subset(rtc, select = c(
      'CLU_SIZE', 'IDCLU'
    )))$CLU_SIZE[(max(rtc$IDCLU) - 10 + 1):max(rtc$IDCLU)]
  )
}

# Combine with sample info
setnames(rtc, 'ID', 'UNIT_ID')
pty.runs <- merge(rtc, pty.runs, by = 'UNIT_ID', all.x = T)

# Write processed samples
if (args$if_save_data) {
  cat('\nWriting to file ',
      file.path(
        args$out.dir,
        paste0(
          'phscinput_runs_cutoff',
          gsub('\\.', '', args$window_cutoff),
          '_ncontrol',args$n_control,
          '.rds'
        )
      ))
  saveRDS(pty.runs,
          file = file.path(
            args$out.dir,
            paste0(
              'phscinput_runs_cutoff',
              gsub('\\.', '', args$window_cutoff),
              '_ncontrol',args$n_control,
              '.rds'
            )
          ))
}
