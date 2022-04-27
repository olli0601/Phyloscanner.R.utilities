cat("\n\n====== TSIcreate_network.R ======\n\n")
# The phyloscanner run - divide individuals into batches for the purpose of TSI estimation.

# Preamble
# The set of scripts aims to run phyloscanner.
# This script aims to break individuals into batches
# based on the pre-calculated similarity scores over the genome

# Load the required packages
library(data.table)
library(seqinr)
library(tidyverse)
library(igraph)
library(dplyr)
library(ggplot2)

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
    "--cluster_size",
    type = "integer",
    default = 100L,
    help = "Maximum cluster size [default %default]",
    dest = "cluster_size"
  ),
  optparse::make_option(
    "--n_control",
    type = "integer",
    default = 0L,
    help = "Number of controls added [default %default]",
    dest = "n_control"
  ),
  optparse::make_option(
    "--n_merge",
    type = "integer",
    default = 100L,
    help = "Maximum cluster sizes to which small clusters were merged [default %default]",
    dest = "n_merge"
  ),
  optparse::make_option(
    "--n_pos",
    type = "integer",
    default = 750L,
    help = "Number of nucleotide locations needed to be considered in phylogenetic analysis [default %default]",
    dest = "n_pos"
  ),
  optparse::make_option(
    "--add_couple_control",
    action = "store_true",
    default = FALSE,
    help = "Add controls by known couples [default]",
    dest = 'if_add_couple_control'
  ),
  optparse::make_option(
    "--save_data",
    action = "store_true",
    default = TRUE,
    help = "Save data [default]",
    dest = 'if_save_data'
  ),
  optparse::make_option(
    "--save_plots",
    action = "store_true",
    default = TRUE,
    help = "Save plots [default]",
    dest = 'if_save_plots'
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
  )
)
#
args <-
  optparse::parse_args(optparse::OptionParser(option_list = option_list))
#
# test
#
if (0) {
  args <- list(
    verbose = T,
    seed = 42,
    cluster_size = 100,
    n_control = 0,
    n_pos = 750,
    n_merge = 100,
    if_add_couple_control = F,
    if_save_data = T,
    if_save_plots = T,
    out.dir = NA,
    prj.dir = NA
  )
}
#
# use manually specified directories when args$out.dir is NA
#
tmp <- Sys.info()
if (tmp["user"] == "xx4515")
{
  if (is.na(args$out.dir))
  {
    # args$out.dir <-
    #   "/rds/general/project/ratmann_deepseq_analyses/live/testTSI/"
    args$out.dir <-
      "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/"
  }
  if (is.na(args$prj.dir))
  {
    args$prj.dir <-
      "~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/"
  }
}

if(tmp["user"] == "andrea")
{
        if(is.na(args$out.dir))
        {
                args$out.dir <- "~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_MRC_TSI/"
        }
        if(is.na(args$prj.dir))
        {
                args$prj.dir <- "/home/andrea/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software"
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
out.dir <- file.path(args$out.dir, 'potential_network')
dir.create(file.path(out.dir, 'plots'))
infile.ind.rccs <-
  file.path(dir.data,
            'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
infile.ind.mrc <-
  file.path(dir.data,
            'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
infile.run <-
  '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_phscinput_samples.rds'
if (args$if_add_couple_control) {
  infile.couple <-
    '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/RakaiPangeaMetaData_v2.rda'
}
#
# Set seed
#
set.seed(args$seed)

if(0)
{
        infile.ind.mrc <- '/home/andrea/Documents/Box/ratmann_pangea_deepsequencedata/PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv'
        infile.ind.rccs <- '/home/andrea/Documents/Box/ratmann_pangea_deepsequencedata/PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv'
}

cat(' ---------- Load files ---------- \n')
files <- list.files(out.dir, full.names = T)
files <- files[grepl('similarity[0-9]+.rds', files)]
df <- data.table()
for (file in files) {
  if (args$verbose) {
    cat('processing ', file, '...\n')
  }
  df <- rbind(df, data.table(readRDS(file)))
}

# check that all similarity have run 
# (IDEA: could use return codes to tell the HPC to rerun previous qsub job!)
stopifnot(df[, .N==choose(length(unique(c(H1, H2))), 2)])

#
# set unknown distance
#
df[LENGTH == -1, PERC := NA]
df[LENGTH <= args$n_pos, PERC := NA]
#

#dfrccs <- data.table(read.csv(infile.ind.rccs))
dfrccs <- fread(infile.ind.rccs)
dfrccs <- unique(subset(dfrccs, select = c('pt_id', 'pangea_id')))
dfrccs[, pangea_id := paste0('RCCS_', pangea_id)]
#dfmrc <- data.table(read.csv(infile.ind.mrc))
dfmrc <- fread(infile.ind.mrc)
dfmrc <- unique(subset(dfmrc, select = c('pt_id', 'pangea_id')))
dfmrc[, pangea_id := paste0('MRCUVRI_', pangea_id)]
dfrccsmrc <- rbind(dfrccs, dfmrc)

df <- merge(df,
            dfrccsmrc,
            by.x = 'H1',
            by.y = 'pangea_id',
            all.x = T)
df <- merge(df,
            dfrccsmrc,
            by.x = 'H2',
            by.y = 'pangea_id',
            all.x = T)

# Swap column names so that 1st col always comes first in lexicographical order
df <- df[, .(pt_id.x, pt_id.y, PERC)]
c("pt_id.x", "pt_id.y")
df[, c("pt_id.x", "pt_id.y") := lapply(.SD, as.character), .SDcols= c("pt_id.x", "pt_id.y")]
tmp <- df[pt_id.x > pt_id.y, list(tmp=pt_id.x, pt_id.x, pt_id.y, PERC)] 
tmp <- tmp[, list(pt_id.x=pt_id.y, pt_id.y=tmp, PERC)]
df <- rbind(df[pt_id.x <= pt_id.y], tmp)
stopifnot(df[pt_id.x > pt_id.y, .N == 0])
df <- df[pt_id.x != pt_id.y]


cat(' ---------- Calculate similarity between individuals ---------- \n')
df <-
  df[, list(SIMILARITY = mean(PERC, na.rm = T)), by = c('pt_id.x', 'pt_id.y')]
pt_id <- unique(c(df$pt_id.x, df$pt_id.y))
tmp <- data.table(expand.grid(pt_id.x = pt_id,
                              pt_id.y = pt_id))
tmp[, pt_id.x := as.character(pt_id.x)]
tmp[, pt_id.y := as.character(pt_id.y)]
tmp <- tmp[pt_id.x <= pt_id.y]
tmp <-
  tmp[!paste0(pt_id.x, '-', pt_id.y) %in% paste0(df$pt_id.x, '-', df$pt_id.y)]
tmp[, SIMILARITY := NA]
df <- rbind(df, tmp)
setkey(df, pt_id.x, pt_id.y)

# Clustering algorithm based on source code here: https://github.com/jmonlong/Hippocamplus/blob/master/content/post/2018-06-09-ClusterEqualSize.Rmd 
cluster.fixed.size <- function(dmat, clsize = 100){

        #         dmat <- copy(df)
        #         clsize=100

        clsize.rle = rle(as.numeric(cut(1:nrow(dmat), ceiling(nrow(dmat)/clsize))))
        clsize = clsize.rle$lengths
        # Performs clustering with clusters of determined sizes.
        # At each step, pick the 
        lab = rep(NA, nrow(dmat))
        cpt = 1
        while(sum(is.na(lab)) > 0){
                lab.ii = which(is.na(lab))
                dmat.m = dmat[lab.ii,lab.ii]
         
                # Pick point which is furthest from remaining unclassified
                ii = which.max(rowSums(dmat.m))

                lab.m = rep(NA, length(lab.ii))
                lab.m[head(order(dmat.m[ii,]), clsize[cpt])] = cpt
                lab[lab.ii] = lab.m
                cpt = cpt + 1
        }
        if(any(is.na(lab))){
                lab[which(is.na(lab))] = cpt
        }
  lab
}


if(1)
{
        cat(' ------ Make Clusters ------ \n')
        df[pt_id.x == pt_id.y, SIMILARITY := 1] 

        df <- dcast(df, pt_id.x ~ pt_id.y, value.var = 'SIMILARITY')
        rownames(df) <- df$pt_id.x
        df[, pt_id.x:=NULL]
        stopifnot(all(rownames(df)==colnames(df)))
        df <- t(df)
        colnames(df) <- rownames(df)

        # Get distance matrix
        mean(is.na(df[lower.tri(df)]))
        all(is.na(df[upper.tri(df)]))
        df[upper.tri(df)] <- df[lower.tri(df)] 
        # isSymmetric.matrix(df)
        df <- sqrt(1-df)
        df[which(is.na(df))] <- min(df, na.rm=T)

        clus <- cluster.fixed.size(df, clsize=100)
        rtc <- data.table(ID=colnames(df), IDCLU=clus)
        # setkey(tmp,cluster, pt_id)
        rtc[, CLU_SIZE:=.N, by='IDCLU']
        
}else{

        cat(' ---------- Make clusters  ---------- \n')
        df <- dcast(df, pt_id.x ~ pt_id.y, value.var = 'SIMILARITY')
        rownames(df) <- df$pt_id.x
        df[, pt_id.x:=NULL]
        stopifnot(all(rownames(df)==colnames(df)))
        df <- t(df)

# AB: I don t think that as dist makes sense here.
# AB: The percentages are already inversely proportional to a pairwise distance
        dist <- as.dist(1 - df)
        dist[is.na(dist)] <- max(dist[!is.na(dist)])
        dist_clu <- hclust(dist, method = 'ward.D')
        if (args$if_save_plots) {
          pdf(file.path(out.dir, 'cluster.pdf'),
              width = 10,
              height = 8)
          plot(dist_clu)
          dev.off()
        }
# cluster
# I did not apply hclust directly in case only one large clusters left and others are of size 1
# the code increases the number of clusters from 2,
# identifies the cluster of size < pre-defined cluster_size,
# deletes the individuals from distance matrix,
# and repeats the procedure to the rest of individuals.
# the procedure stops if there are less than cluster_size individuals left in the distance matrix,
# or all the newly identified clusters of size < cluster_size
#
# number of clusters specified for cutree
        k <- 2
# cluster ID
        count <- 1
        while (T) {
          if (args$verbose) {
            cat('Divided into ', k, ' clusters ... \n')
          }
          dist_cut <- cutree(dist_clu, k)
          if (any(table(dist_cut) < args$cluster_size)) {
            ids <- which(table(dist_cut) < args$cluster_size)
            for (id in ids) {
              id_name <- names(dist_cut)[dist_cut == id]
              if (count == 1) {
                rtc <- data.table(CLU = count,
                                  ID = id_name)
              } else{
                tmp <- data.table(CLU = count,
                                  ID = id_name)
                rtc <- rbind(rtc, tmp)
              }
              tmp <- which(dist_cut == id)
              dist <- as.matrix(dist)
              dist <- dist[-tmp, -tmp]
              if (nrow(dist) <= args$cluster_size) {
                tmp <- data.table(CLU = count + 1,
                                  ID = rownames(dist))
                rtc <- rbind(rtc, tmp)
                break
              }
              dist <- as.dist(dist)
              count <- count + 1
            }
          }
          if (all(table(dist_cut) < args$cluster_size)) {
            break
          }
          k <- k + 1
          dist_clu <- hclust(dist, method = 'ward.D')
        }

        tmp	<- rtc[, list(CLU_SIZE = length(ID)), by = 'CLU']
        setkey(tmp, CLU_SIZE)
        if (args$verbose) {
          cat('---------- Largest cluster sizes ----------\n')
          print(tail(tmp))
          cat('---------- Smallest cluster sizes ----------\n')
          print(head(tmp))
        }
        tmp[, IDCLU := seq_len(nrow(tmp))]
        rtc	<- subset(merge(rtc, tmp, by = 'CLU'))
        rtc[, CLU := NULL]
        setkey(rtc, IDCLU)
}



if (args$if_save_data) {
  cat('Writing files to ',
      file.path(out.dir, 'clusters.rds'),
      '... \n')
  saveRDS(rtc, file = file.path(out.dir, 'clusters.rds'))
}

if (args$if_add_couple_control) {
  cat('----------- Add known couples as controls ----------- \n')
  pty.runs <- readRDS(infile.pty.runs)
  #	Add couples
  load(infile.couple)
  rp <-
    data.table(unique(subset(
      coupdat,
      select = c('male.RCCS_studyid', 'female.RCCS_studyid')
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
  tmp <- subset(rp,!is.na(FEMALE_IDCLU) & !is.na(MALE_IDCLU))
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
  set(rp,
      tmp,
      c('FEMALE_IDCLU', 'MALE_IDCLU'),
      rtc[, max(IDCLU)] + seq_along(tmp))
  
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
}

cat('----------- Merge small runs ----------- \n')
tmp	<- unique(subset(rtc, select = c(IDCLU, CLU_SIZE)))
tmp[, PTY_RUN := IDCLU]
tmp[, PTY_SIZE := CLU_SIZE]
sum_reset_at <- function(thresh) {
  function(x) {
    accumulate(x, ~ if_else(.x + .y > thresh, .y, .x + .y))
  }
}
tmp <-
  tmp %>% mutate(PTY_SIZE2 = sum_reset_at(args$n_merge)(PTY_SIZE))
tmp <- as.data.table(tmp)
tmp2 <- which(diff(tmp$PTY_SIZE2) < 0) + 1
tmp[, PTY_RUN2 := 0]
tmp[c(tmp2, max(tmp2):nrow(tmp)), PTY_RUN2 := 1]
tmp[, PTY_RUN2 := cumsum(PTY_RUN2)]
tmp[, PTY_SIZE2 := max(PTY_SIZE2), by = 'PTY_RUN2']
tmp[, PTY_RUN2 := PTY_RUN2 + 1]
set(tmp, NULL, c('PTY_RUN', 'PTY_SIZE'), NULL)
setnames(tmp, c('PTY_RUN2', 'PTY_SIZE2'), c('PTY_RUN', 'PTY_SIZE'))
rtc <- merge(rtc, tmp, by = c('IDCLU', 'CLU_SIZE'))
tmp	<- rtc[, list(PTY_SIZE = length(ID)), by = 'PTY_RUN']
setkey(tmp, PTY_SIZE)
if (args$verbose) {
  cat('---------- Largest phyloscanner run sizes ----------\n')
  print(tail(tmp))
  cat('---------- Smallest phyloscanner run sizes ----------\n')
  print(head(tmp))
}

# TO ADD N_CONTROL
# AFTER DEVELOPING ANOTHER SCRIPT...

# Combine with sample info
setnames(rtc, 'ID', 'UNIT_ID')
pty.runs <- data.table(readRDS(infile.run))
pty.runs <- merge(rtc, pty.runs, by = 'UNIT_ID', all.x = T)

# Write processed samples
if (args$if_save_data) {
  cat('\nWriting to file ',
      file.path(
        args$out.dir,
        paste0(
          'phscinput_runs_clusize_',
          args$cluster_size,
          '_ncontrol_',
          args$n_control,
          '.rds'
        )
      ))
  saveRDS(pty.runs,
          file =  file.path(
            args$out.dir,
            paste0(
              'phscinput_runs_clusize_',
              args$cluster_size,
              '_ncontrol_',
              args$n_control,
              '.rds'
            )
          ))
}



