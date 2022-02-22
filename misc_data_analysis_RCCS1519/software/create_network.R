# The phyloscanner run  - divide individuals into batches for the purpose of reconstructing transmissions.

# Preamble
# The set of scripts aims to run phyloscanner.
# This script aims to break individuals into batches
# based on the pre-calculated similarity scores over the windows.

# Load the required packages
library(data.table)
library(seqinr)
library(tidyverse)
library(igraph)
library(dplyr)
library(rstan)

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
    c("-o", "--overwrite"),
    action = "store_true",
    default = TRUE,
    help = "Print extra output [default]",
    dest = "overwrite"
  ),
  optparse::make_option(
    "--seed",
    type = "integer",
    default = 42L,
    help = "Random number seed [default %default]",
    dest = "seed"
  ),
  optparse::make_option(
    "--length_cutoff",
    type = "numeric",
    default = 250L,
    help = "Minumum length of genome to be considered [default %default]",
    dest = "length_cutoff"
  ),
  optparse::make_option(
    "--cluster_size",
    type = "integer",
    default = 55L,
    help = "Maximum cluster size [default %default]",
    dest = "cluster_size"
  ),
  optparse::make_option(
    "--n_merge",
    type = "integer",
    default = 8L,
    help = "Maximum cluster sizes to which small clusters were merged [default %default]",
    dest = "n_merge"
  ),
  optparse::make_option(
    "--window_size",
    type = "integer",
    default = 500L,
    help = "Number of pairs per job for the similarity score calculation [default %default]",
    dest = "window_size"
  ),
  optparse::make_option(
    "--sliding_width",
    type = "integer",
    default = 100L,
    help = "Sliding width for the similarity score calculation [default %default]",
    dest = "sliding_width"
  ),
  optparse::make_option(
    "--high_depth_length_cutoff",
    type = "numeric",
    default = 250L,
    help = "Minumum length of genome with high depth [default %default]",
    dest = "depth_cutoff"
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
    default = 3,
    help = "Number of controls added [default %default]",
    dest = "n_control"
  ),
  optparse::make_option(
    "--n_overlapping_windows",
    type = "integer",
    default = 4,
    help = "Minumum number of overlapping windows required [default %default]",
    dest = "n_overlap"
  ),
  optparse::make_option(
    "--n_iteration",
    type = "integer",
    default = 1e4,
    help = "Number of iterations for fitting the mixture model [default %default]",
    dest = "iter"
  ),
  optparse::make_option(
    "--add_couple_control",
    action = "store_true",
    default = FALSE,
    help = "Add controls by known couples [default]",
    dest = 'if_add_couple_control'
  ),
  optparse::make_option(
    "--n_warmup",
    type = "integer",
    default = 2e3,
    help = "Number of warm-up for fitting the mixture model [default %default]",
    dest = "warmup"
  ),
  optparse::make_option(
    "--save_data",
    action = "store_true",
    default = TRUE,
    help = "Save data [default]",
    dest = 'if_save_data'
  ),
  optparse::make_option(
    "--infile",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to data directory where consensus sequences are stored [default]",
    dest = 'infile'
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

args <-
  optparse::parse_args(optparse::OptionParser(option_list = option_list))

#
# test
#
if (0) {
  args <- list(
    verbose = T,
    seed = 42,
    length_cutoff = 250L,
    depth_cutoff = 250L,
    cluster_size = 55L,
    window_cutoff = 0.5,
    sliding_width = 100L,
    window_size = 500L,
    n_control = 3L,
    if_add_couple_control = T,
    n_merge = 8L,
    iter = 1e4,
    warmup = 2e3,
    if_save_data = T,
    date = as.character(Sys.Date()),
    out.dir = NA,
    prj.dir = NA,
    infile = NA
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
    args$out.dir <-
      "/rds/general/project/ratmann_deepseq_analyses/live/test_full/"
  }
  if (is.na(args$prj.dir))
  {
    args$prj.dir <-
      "~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/"
  }
  if (is.na(args$infile))
  {
    args$infile <-
      "/rds/general/project/ratmann_pangea_deepsequencedata/live/200422_PANGEA2_RCCSMRC_alignment.fasta"
  }
}

# if prj.dir and out.dir are not manually set, default to here()
if (is.na(args$prj.dir))
{
  args$prj.dir <- here::here()
}
if (is.na(args$out.dir))
{
  args$out.dir <- here::here()
}
if (is.na(args$infile)) {
  stop("Please input the sequence file. ")
}
#
# add constants that should not be changed by the user
#
dir.data <-
  '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
dir.rccs1519 <-
  '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
infile.couple <-
  file.path(dir.rccs1519,
            'RakaiPangeaMetaData_v2.rda')
infile.ind.rccs <-
  file.path(dir.data,
            'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
infile.ind.mrc <-
  file.path(dir.data,
            'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')
args$infile.rccs <-
  file.path(dir.data,
            'PANGEA2_RCCS/200422_PANGEA2_RCCS_mapped_samples.rds')
args$infile.mrc <-
  file.path(dir.data,
            'PANGEA2_MRC/200422_PANGEA2_MRCUVRI_mapped_samples.rds')

args$infile.pty.runs <-
  file.path(dir.rccs1519, '210120_RCCSUVRI_phscinput_samples.rds')
out.dir <- file.path(args$out.dir, 'potential_network')
infile.threshold <-
  '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/200929_SIMILARTY_windowsize500_batchsize100/threshold_5quantile_indep_exp5.rda'


# Define functions just used in this script alone
take_bridging_individual <- function(df_graph, df) {
  bridging = list()
  for (k in 1:max(df$MEMBERSHIP)) {
    tmp1 = df[MEMBERSHIP == k, ]$PT_ID
    tmp2 = df_graph[pt_id1 %in% tmp1 &
                      !(pt_id2 %in% tmp1) & SIMILARITY >= 0.975, ]
    tmp3 = df_graph[pt_id2 %in% tmp1 &
                      !(pt_id1 %in% tmp1) & SIMILARITY >= 0.975, ]
    bridging[[k]] = unique(c(as.character(tmp2$pt_id1),
                             as.character(tmp3$pt_id2)))
  }
  return(bridging)
}
break_large_clusters <- function(idlarge, dflarge, ddist) {
  ans <- data.table()
  for (i in idlarge) {
    idlarge_ind <- dflarge[MEMBERSHIP == i]$PT_ID
    maxsize <- unique(dflarge[MEMBERSHIP == i, ]$CLU_SIZE)
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

cat(' ---------- Collect similarity scores ---------- \n')
files <- list.files(out.dir, pattern = '*.rds')
if (any(grepl('similarity([0-9]+)_window_([0-9]+).rds', files))) {
  df <- data.table(FILE = files)
  df <- df[grep('similarity([0-9]+)_window_([0-9]+).rds', FILE)]
  df[, BATCH := as.numeric(gsub('similarity([0-9]+)_window_([0-9]+).rds', '\\1', FILE))]
  df[, WINDOW := as.numeric(gsub('similarity([0-9]+)_window_([0-9]+).rds', '\\2', FILE))]
  #
  for (window in unique(df$WINDOW)) {
    if (args$verbose) {
      cat('Process window ', window, '...\n')
    }
    dfbatch <- df[WINDOW == window, {
      readRDS(file.path(
        out.dir,
        paste0('similarity', BATCH, '_window_', window, '.rds')
      ))
    },
    by = 'BATCH']
    dfbatch[, WINDOW := window]
    
    # summary
    sbatch <- dfbatch[, list(
      NPAIR = length(LENGTH),
      NEPAIR = sum(LENGTH != -1L),
      MAXLEN = max(LENGTH[LENGTH > 0]),
      MINLEN = min(LENGTH[LENGTH > 0]),
      MAXPERC = max(PERC[LENGTH > 0]),
      MINPERC = min(PERC[LENGTH > 0])
    )]
    dfbatch <- dfbatch[LENGTH != -1]
    # save
    save(dfbatch, sbatch, file = file.path(out.dir, paste0('similarity', window, '.rda')))
  }
  file.remove(file.path(out.dir, df$FILE))
  infile.similarity <-
    file.path(out.dir, paste0('similarity', unique(df$WINDOW), '.rda'))
}else{
  infile.similarity <-
    list.files(out.dir, pattern = '*.rda', full.names = T)
  infile.similarity <-  grep('similarity[0-9]+', infile.similarity,value = T)
}

#
cat(' ---------- Determine the thresholds ---------- \n')
# couple
load(infile.couple)
couple <-
  unique(subset(
    data.table(coupdat),
    select = c('male.RCCS_studyid', 'female.RCCS_studyid')
  ))
couple[, COUPLE := 1]
couple[, male.RCCS_studyid := paste0('RK-', male.RCCS_studyid)]
couple[, female.RCCS_studyid := paste0('RK-', female.RCCS_studyid)]
couple <-
  subset(couple,
         select = c('male.RCCS_studyid', 'female.RCCS_studyid', 'COUPLE'))
tmp <- copy(couple)
setnames(couple,
         c('male.RCCS_studyid', 'female.RCCS_studyid'),
         c('pt_id1', 'pt_id2'))
setnames(tmp,
         c('male.RCCS_studyid', 'female.RCCS_studyid'),
         c('pt_id2', 'pt_id1'))
couple <- rbind(couple, tmp)
couple  <- couple[pt_id1 != pt_id2, ]

# map ids
ids <- data.table(read.csv(infile.ind.rccs))
ids <- subset(ids, select = c("pt_id", "sex", "pangea_id"))
ids[, pangea_id := paste0('RCCS_', pangea_id)]
tmp <- data.table(read.csv(infile.ind.mrc))
tmp <- subset(tmp, select = c("pt_id", "sex", "pangea_id"))
tmp[, pangea_id := paste0('MRCUVRI_', pangea_id)]
ids <- rbind(ids, tmp)
ids <- unique(ids)

if(!file.exists(file.path(
  out.dir,
  paste0('similarity_couple_minlen_',
         args$length_cutoff,
         '.rds')
))|args$overwrite){
  # similarity among couples
  for (file in infile.similarity) {
    if (args$verbose) {
      cat('Process ', file, '... \n')
    }
    load(file)
    dfbatch <- dfbatch[LENGTH >= args$length_cutoff,]
    # ids
    setnames(ids, colnames(ids), paste0(colnames(ids), '1'))
    dfbatch <-
      merge(
        dfbatch,
        ids,
        by.x = c('H1'),
        by.y = c('pangea_id1'),
        all.x = T
      )
    setnames(ids, colnames(ids), gsub('1', '2', colnames(ids)))
    dfbatch <-
      merge(
        dfbatch,
        ids,
        by.x = c('H2'),
        by.y = c('pangea_id2'),
        all.x = T
      )
    setnames(ids, colnames(ids), gsub('2', '', colnames(ids)))
    dfbatch <-
      dfbatch[!is.na(pt_id1) & !is.na(pt_id2) & pt_id1 != pt_id2,]
    dfbatch <- merge(dfbatch,
                     couple,
                     by = c('pt_id1', 'pt_id2'),
                     all.x = T)
    dfbatch <- dfbatch[COUPLE == 1, c('pt_id1', 'pt_id2', 'PERC')]
    window <-
      as.numeric(gsub('similarity([0-9]+).rda', '\\1', basename(file)))
    dfbatch[, WINDOW := window]
    if (!exists("ans")) {
      ans <- dfbatch
    } else{
      ans <- rbind(ans, dfbatch)
    }
  }
  
  # Save
  if (args$if_save_data) {
    cat('Write similarity among couples to ',
        file.path(
          out.dir,
          paste0('similarity_couple_minlen_',
                 args$length_cutoff,
                 '.rds')
        ),
        '...\n')
    saveRDS(ans, file = file.path(
      out.dir,
      paste0('similarity_couple_minlen_',
             args$length_cutoff,
             '.rds')
    ))
  }
  gc()
}else{
  ans <- data.table(readRDS(file.path(
    out.dir,
    paste0('similarity_couple_minlen_',
           args$length_cutoff,
           '.rds'))))
}



#
if ((nrow(ans) != 0 & !file.exists(file.path(
  out.dir,
  paste0(
    'similarity_couple_fitted_thresholds_minlen_',
    args$length_cutoff,
    '.rda'))))|args$overwrite) {
  setkey(ans, pt_id1,  pt_id2)
  cat('Fit the two-component Gaussian mixture model ... \n')
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
  for (window in 1:max(ans$WINDOW)) {
    if (args$verbose) {
      cat('Process window ', window, ' ...\n')
    }
    tmp <- list()
    tmp$y <- ans[WINDOW == window]$PERC
    tmp$N <- length(tmp$y)
    fit[[window]] <- stan(
      model_code = stan_code,
      data = tmp,
      iter = args$iter,
      warmup = args$warmup,
      chains = 1,
      seed = args$seed
    )
  }
  
  if (args$if_save_data) {
    cat('Write the fitted models to ',
        file.path(
          out.dir,
          paste0(
            'similarity_couple_fitted_model_minlen_',
            args$length_cutoff,
            '.rds'
          )
        ),
        '...\n')
    saveRDS(fit, file = file.path(
      out.dir,
      paste0(
        'similarity_couple_fitted_model_minlen_',
        args$length_cutoff,
        '.rds'
      )
    ))
  }
  
  fit <- readRDS(file.path(
    out.dir,
    paste0(
      'similarity_couple_fitted_model_minlen_',
      args$length_cutoff,
      '.rds'
    )
  ))
  cat(' ---------- Check the fitted models ---------- \n')
  rh = c()
  ness = c()
  tess = c()
  bess = c()
  for (window in 1:max(ans$WINDOW)) {
    cat('Process window ', window, '\n')
    pars <-
      rstan::extract(fit[[window]], pars = names(fit[[window]])[!grepl('lp_', names(fit[[window]]))])
    pars = do.call(cbind, pars)
    rh = c(rh, max(summary(fit[[window]])$summary[, "Rhat"]))
    ness = c(ness, min(summary(fit[[window]])$summary[, "n_eff"]))
  }
  cat('Range of rhat is ', range(rh), '\n')
  cat('Range of ess is ', range(ness), '\n')
  
  
  cat(' ---------- Extract fitted parameters ---------- \n')
  df = data.table()
  for (window in 1:max(ans$WINDOW)) {
    cat('Process window ', window, '\n')
    tmp = data.table(summary(fit[[window]])$summary[1:5, c("2.5%", "50%", "97.5%")])
    tmp[, VARIABLE := rownames(summary(fit[[window]])$summary)[1:5]]
    tmp[, WINDOW := window]
    df = rbind(df, tmp)
  }
  
  if (args$if_save_data) {
    cat('Write the fitted model parameters to ',
        file.path(
          out.dir,
          paste0(
            'similarity_couple_fitted_model_parameters_minlen_',
            args$length_cutoff,
            '.rds'
          )
        ),
        '...\n')
    saveRDS(df, file = file.path(
      out.dir,
      paste0(
        'similarity_couple_fitted_model_parameters_minlen_',
        args$length_cutoff,
        '.rds'
      )
    ))
  }
  
  cat(' ---------- Extract thresholds ---------- \n')
  threshold <- list()
  for (window in 1:max(ans$WINDOW)) {
    cat('Process window ', window, '\n')
    pars <-
      rstan::extract(fit[[window]], pars = c('mu', 'sigma', 'theta'))
    tmp <- c()
    for (j in 1:nrow(pars$mu)) {
      tmp = c(tmp, qlnorm(0.05, pars$mu[j, 2], pars$sigma[j, 2]))
    }
    threshold[[window]] = tmp
  }
  threshold = lapply(threshold, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  })
  threshold = data.table(do.call(rbind, threshold))
  if (args$if_save_data) {
    cat('Write the thresholds to ',
        file.path(
          out.dir,
          paste0(
            'similarity_couple_thresholds_minlen_',
            args$length_cutoff,
            '.rda'
          )
        ),
        '...\n')
    save(threshold, file = file.path(
      out.dir,
      paste0(
        'similarity_couple_fitted_thresholds_minlen_',
        args$length_cutoff,
        '.rda'
      )
    ))
  }
} else if(file.exists(file.path(
  out.dir,
  paste0(
    'similarity_couple_fitted_thresholds_minlen_',
    args$length_cutoff,
    '.rda')))){
  infile.threshold <- file.path(
    out.dir,
    paste0(
      'similarity_couple_fitted_thresholds_minlen_',
      args$length_cutoff,
      '.rda'))
  load(infile.threshold)
}else{
  load(infile.threshold)
}
rm(list=setdiff(ls(), c('couple','ids','threshold','args',
                        'infile.similarity','out.dir')))
gc()

cat(' ----------- Select high depth sequences ---------- ')
alignment <- read.fasta(file = args$infile)
npos <- unique(lengths(alignment))
windows_last_start <-
  ceiling(npos / args$sliding_width) * args$sliding_width - args$window_size + 1
windows_2last_end <-
  floor(npos / args$sliding_width) * args$sliding_width
windows_start <- seq(1, windows_last_start, args$sliding_width)
windows_end <-
  c(seq(args$window_size , windows_2last_end, args$sliding_width),
    npos)
windows <- seq_len(length(windows_start))
dwindows <- data.table(WINDOW = windows,
                       START = windows_start,
                       END = windows_end)
if(!file.exists(gsub('.fasta', '_depths.rds', args$infile)) | args$overwrite){
  mapping_rccs <- readRDS(args$infile.rccs)
  dconsensus <-
    subset(mapping_rccs, select = c('PANGEA_ID', 'CONSENSUS', 'F'))
  dconsensus <- unique(dconsensus[CONSENSUS != "",])
  dconsensus[, PANGEA_ID := paste0('RCCS_', PANGEA_ID)]
  mapping_mrc <- readRDS(args$infile.mrc)
  tmp <-
    subset(mapping_mrc, select = c('PANGEA_ID', 'CONSENSUS', 'F'))
  tmp <- unique(tmp[CONSENSUS != "",])
  tmp[, PANGEA_ID := paste0('MRCUVRI_', PANGEA_ID)]
  dconsensus = rbind(dconsensus, tmp)
  dconsensus = subset(dconsensus, select = c('PANGEA_ID', 'F'))
  setnames(dconsensus, c('PANGEA_ID', 'F'), c('pangea_id', 'file_name'))
  tmp <-  data.table(pangea_id = names(alignment))
  dconsensus <- merge(dconsensus, tmp, by = 'pangea_id', all.y = T)
  dconsensus <- unique(dconsensus)
  
  # valid lengths and high depth position lengths in consensus sequences
  ddepth = data.table()
  for (i in seq_len(nrow(dconsensus))) {
    if (i %% 100 == 1 | i == nrow(dconsensus)) {
      cat('processing ', i, 'th out of ', nrow(dconsensus), ' files \n')
    }
    sequence <- readLines(dconsensus$file_name[i])
    sequence <- paste(sequence[2:length(sequence)], collapse = "")
    sequence <- strsplit(sequence, "")
    
    tmp_df <- dwindows[, {
      tmp = sequence[[1]][START:END]
      tmp = tmp[!grepl("[^A-Za-z]", tmp)]
      list(
        length = length(tmp),
        high_depth = sum(tmp == toupper(tmp)),
        pangea_id = dconsensus$pangea_id[i],
        file_name = dconsensus$file_name[i]
      )
    }, by = c('WINDOW', 'START', 'END')]
    ddepth <- rbind(ddepth, tmp_df)
  }
  #
  if (args$if_save_data) {
    cat('Write the depths to ',
        gsub('.fasta', '_depths.rds', args$infile),
        '...\n')
    saveRDS(ddepth, file =  gsub('.fasta', '_depths.rds', args$infile))
  }
}else{
 ddepth <- data.table(readRDS(gsub('.fasta', '_depths.rds', args$infile)))
}

rm(list=setdiff(ls(), c('couple','ids','threshold','ddepth','args',
                        'windows_start','infile.similarity','out.dir')))
gc()

# select
ddepth <- unique(ddepth[high_depth >= args$depth_cutoff])
ddepth_windows <-
  ddepth[, list(pangea_id = unique(pangea_id)), by = c('WINDOW')]
ddepth_windows[, HD := 1]
# ddepth_windows <-
#   data.table::dcast(ddepth_windows, pangea_id ~ WINDOW, value.var = 'HD')

rm(list=setdiff(ls(), c('couple','ids','threshold','ddepth_windows','args',
                        'windows_start','infile.similarity','out.dir')))
gc()

cat(' ---------- Summarise closeness between sequences by windows ----------')
threshold <- threshold$`50%`
for (window in 1:length(windows_start)) {
  cat('\n Process window ', window , '\n')
  load(infile.similarity[window])
  dfbatch <- dfbatch[LENGTH >= args$length_cutoff,]
  
  # select individuals with high depth
  ddepth_select <-
    unique(subset(ddepth_windows[WINDOW == window,], select = c('pangea_id','HD')))
  setnames(ddepth_select,
           colnames(ddepth_select),
           paste0(colnames(ddepth_select), '1'))
  dfbatch <-
    merge(
      dfbatch,
      ddepth_select,
      by.x = 'H1',
      by.y = 'pangea_id1',
      all.x = T
    )
  setnames(ddepth_select,
           colnames(ddepth_select),
           gsub('1', '2', colnames(ddepth_select)))
  dfbatch <-
    merge(
      dfbatch,
      ddepth_select,
      by.x = 'H2',
      by.y = 'pangea_id2',
      all.x = T
    )
  setnames(ddepth_select,
           colnames(ddepth_select),
           gsub('2', '', colnames(ddepth_select)))
  dfbatch <- dfbatch[HD1 == 1 & HD2 == 1,]
  
  # add thresholds
  dfbatch[, threshold_window := threshold[window]]
  dfbatch[, CLOSE := as.integer(PERC > threshold_window)]
  
  # average similarity per sequence pair
  tmp  <-
    unique(subset(dfbatch, select = c('H1', 'H2', 'PERC')))
  tmp <-
    tmp[, list(PERC = mean(PERC, na.rm = T)), by = c('H1', 'H2')]
  setnames(tmp, 'PERC', paste0('PERC', window))
  if (window == 1) {
    df <- copy(tmp)
  } else{
    df <- merge(df,
                tmp,
                by = c('H1', 'H2'),
                all = T)
  }
  rm(list=setdiff(ls(), c('couple','ids','threshold','ddepth_windows','dfbatch',
                          'df','args','infile.similarity','out.dir','window','dsupport')))
  gc()
  
  # support for close relationship per sequence pair
  tmp  <-
    unique(subset(dfbatch, select = c('H1', 'H2', 'CLOSE')))
  tmp <-
    tmp[, list(CLOSE = mean(CLOSE, na.rm = T)), by = c('H1', 'H2')]
  setnames(tmp, 'CLOSE', paste0('CLOSE', window))
  if (window == 1) {
    dsupport <- copy(tmp)
  } else{
    dsupport <- merge(dsupport,
                      tmp,
                      by = c('H1', 'H2'),
                      all = T)
  }
  rm(list=setdiff(ls(), c('couple','ids','threshold','ddepth_windows',
                          'dsupport','df','args','infile.similarity','out.dir')))
  gc()
}

#
if (args$if_save_data) {
  cat(
    'Write the average similarity scores to ',
    file.path(out.dir, 'average_similarity.rds'),
    '... \n'
  )
  saveRDS(df, file = file.path(out.dir, 'average_similarity.rds'))
  cat(
    'Write the supports for closeness to ',
    file.path(out.dir, 'support_close.rds'),
    '... \n'
  )
  saveRDS(dsupport, file = file.path(out.dir, 'support_close.rds'))
}
rm(list=setdiff(ls(), c('couple','ids','dsupport','args','out.dir')))
gc()

cat('---------- Classify the closeness between sequences ----------')
dsupport <-
  data.table(readRDS(file.path(out.dir, 'support_close.rds')))
tmp <-
  subset(dsupport, select = grep('CLOSE', colnames(dsupport), value = T))
tmp <- as.matrix(tmp)
tmpn <- apply(tmp, 1, function(x)
  sum(!is.na(x)))
tmpp <- apply(tmp, 1, function(x)
  sum(x, na.rm = T))
tmpp <- tmpp / tmpn
dsupport$DATA <- tmpn
dsupport$PROP <- tmpp
dsupport <- dsupport[DATA >= args$n_overlap,]
dsupport <- dsupport[PROP >= args$window_cutoff]

# map
setnames(ids, colnames(ids), paste0(colnames(ids), '1'))
df <-
  merge(
    dsupport,
    ids,
    by.x = c('H1'),
    by.y = c('pangea_id1'),
    all.x = T
  )
setnames(ids, colnames(ids), gsub('1', '2', colnames(ids)))
df <-
  merge(
    df,
    ids,
    by.x = c('H2'),
    by.y = c('pangea_id2'),
    all.x = T
  )
setnames(ids, colnames(ids), gsub('2', '', colnames(ids)))
df = df[!is.na(pt_id1) & !is.na(pt_id2)]
close_pairs <- unique(subset(df, select = c('pt_id1', 'pt_id2')))
close_pairs <- close_pairs[pt_id1 != pt_id2]
close_pairs[, pt_id1 := as.character(pt_id1)]
close_pairs[, pt_id2 := as.character(pt_id2)]

#
if (args$if_save_data) {
  cat('Write the average similarity scores to ',
      file.path(out.dir, paste0(
        'close_pairs_window_cutoff',
        gsub('\\.', '_', args$window_cutoff),
        '.rds'
      )),
      '... \n')
  saveRDS(close_pairs, file = file.path(out.dir, paste0(
    'close_pairs_window_cutoff',
    gsub('\\.', '', args$window_cutoff),
    '.rds'
  )))
}
rm(list=setdiff(ls(), c('couple','ids','close_pairs','args','out.dir')))
gc()

cat(' ----------- Generate clusters ----------- \n')
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
tmp[, IDCLU := seq_len(nrow(tmp))]
rtc	<- subset(merge(rtc, tmp, by = 'CLU'))
rtc[, CLU := NULL]
setkey(rtc, IDCLU)

cat('----------- Break the large clusters ----------- \n')
if (any(tmp$CLU_SIZE >= args$cluster_size)) {
  id = rtc[IDCLU == 1, ]$ID
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
  dflarge <- df[CLU_SIZE > args$cluster_size]
  df <- df[CLU_SIZE <= args$cluster_size]
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
    tmp <- breaked_dflarge[CLU_SIZE <= args$cluster_size]
    tmp2 <-
      unique(subset(tmp, select = c('CLU_SIZE', 'MEMBERSHIP')))
    setkey(tmp2, MEMBERSHIP)
    tmp2[, MEMBERSHIP2 := seq_len(nrow(tmp2))]
    tmp <- merge(tmp, tmp2, by = c('CLU_SIZE', 'MEMBERSHIP'))
    tmp[, MEMBERSHIP := NULL]
    setnames(tmp, 'MEMBERSHIP2', 'MEMBERSHIP')
    tmp[, MEMBERSHIP := MEMBERSHIP + id_clu]
    df <- rbind(df, tmp)
    dflarge <- breaked_dflarge[CLU_SIZE > args$cluster_size]
    idlarge <- unique(dflarge$MEMBERSHIP)
    if (length(idlarge) == 0)
      break
    breaked_dflarge <-
      break_large_clusters(idlarge, dflarge, ddist)
  }
  # Merge
  setnames(df, c('PT_ID', 'MEMBERSHIP'), c('ID', 'IDCLU'))
  rtc = rtc[IDCLU != 1, ]
  rtc[, IDCLU := IDCLU + max(df$IDCLU) - 1]
  df = rbind(rtc, df)
} else{
  df <- copy(rtc)
}

# Save
if (args$if_save_data) {
  cat('Write potential transmission networks to ',
      file.path(out.dir,
                paste0(
                  'clusters_cutoff',
                  gsub('\\.', '', args$window_cutoff),
                  '.rda'
                )),
      '...\n')
  save(df, file = file.path(out.dir,
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
  subset(merge(df, tmp, by = 'IDCLU'),
         select = c('ID', 'IDCLU2', 'CLU_SIZE'))
setnames(rtc, 'IDCLU2', 'IDCLU')
setkey(rtc, IDCLU)

if (args$if_add_couple_control) {
  cat('----------- Add known couples as controls ----------- \n')
  pty.runs <- readRDS(args$infile.pty.runs)
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

# Find n_control closest individuals to all the individuals in each cluster
if (args$n_control != 0) {
  tmp <-
    data.table(readRDS(file.path(out.dir, 'average_similarity.rds')))
  df <- subset(tmp, select = c('H1', 'H2'))
  df$SIMILARITY <- rowMeans(tmp[, 3:ncol(tmp)], na.rm = TRUE)
  setnames(ids, colnames(ids), paste0(colnames(ids), '1'))
  df <-
    merge(df,
          ids,
          by.x = 'H1',
          by.y = 'pangea_id1',
          all.x = TRUE)
  setnames(ids, colnames(ids), gsub('1', '2', colnames(ids)))
  df <-
    merge(df,
          ids,
          by.x = 'H2',
          by.y = 'pangea_id2',
          all.x = TRUE)
  setnames(ids, colnames(ids), gsub('2', '', colnames(ids)))
  df <-
    df[, list(SIMILARITY = mean(SIMILARITY)), by = c('pt_id1', 'pt_id2')]
  
  if (args$if_save_data) {
    cat('Write to ',
        file.path(out.dir, 'whole_genome_average_similarity.rds'),
        '...')
    saveRDS(df ,
            file.path(out.dir, 'whole_genome_average_similarity.rds'))
  }
  dcl <- rtc[, {
    tmp <- df[pt_id1 %in% ID, c('pt_id2', 'SIMILARITY')]
    tmp2 <- df[pt_id2 %in% ID, c('pt_id1', 'SIMILARITY')]
    tmp <- rbind(tmp, tmp2, use.names = F)
    colnames(tmp) <- c('UNIT_ID', 'SIMILARITY')
    tmp <- tmp[!UNIT_ID %in% ID,]
    tmp <- tmp[, list(mean(SIMILARITY)), by = 'UNIT_ID']
    tmp <- tmp[order(V1, decreasing = T),]
    list(
      variable = c(
        paste0('ID_CLOSE', 1:args$n_control),
        paste0('SIMILARITY', 1:args$n_control)
      ),
      value = c(as.character(tmp$UNIT_ID[1:args$n_control]), tmp$V1[1:args$n_control])
    )
  }, by = c('CLU_SIZE', 'IDCLU')]
  dcl <-
    dcast(dcl, CLU_SIZE + IDCLU ~ variable, value.var = 'value')
  tmp <- grep('SIMILARITY', colnames(dcl), value = T)
  dcl[, (tmp) := lapply(.SD, as.numeric), .SDcols = tmp]
  if (args$if_save_data) {
    cat('Write to ', file.path(
      out.dir,
      paste0(args$n_control, 'closest_individuals.rds')
    ), '...\n')
    saveRDS(dcl, file = file.path(
      out.dir,
      paste0(args$n_control, 'closest_individuals.rds')
    ))
  }
  
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
pty.runs <- data.table(readRDS(args$infile.pty.runs))
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
