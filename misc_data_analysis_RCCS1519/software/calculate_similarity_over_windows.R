# The phyloscanner run - calculate similarities between the consensus sequences
# on the overlapping windows

# Preamble
# The set of scripts aims to run phyloscanner.
# This script aims to generate scripts which
# calculate the similarities between the consensus sequences stored in infile
# on the overlapping windows.

# Load the required packages
library(data.table)
library(seqinr)
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
    "--npair_perjob",
    type = "integer",
    default = 10000L,
    help = "Number of pairs per job for the similarity score calculation [default %default]",
    dest = "npair_perjob"
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
    help = "Sliding width [default %default]",
    dest = "sliding_width"
  ),
  optparse::make_option(
    "--save_data",
    action = "store_true",
    default = TRUE,
    help = "Save data [default]",
    dest = 'if_save_data'
  ),
  optparse::make_option(
    "--env_name",
    type = "character",
    default = 'phyloenv',
    help = "Conda environment name [default]",
    dest = 'env_name'
  ),
  optparse::make_option(
    "--walltime",
    type = "integer",
    default = 24,
    help = "Job walltime (hours) [default]",
    dest = 'walltime'
  ),
  optparse::make_option(
    "--memory",
    type = "integer",
    default = 10,
    help = "Job memory (GB) [default]",
    dest = 'memory'
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
    npair_perjob = 10000,
    window_size = 500,
    sliding_width = 100,
    env_name = "phyloenv",
    walltime = 24,
    memory = 10,
    if_save_data = T,
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
      "/rds/general/project/ratmann_deepseq_analyses/live/test/"
  }
  if (is.na(args$prj.dir))
  {
    args$prj.dir <-
      "~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/"
  }
  if (is.na(args$infile))
  {
    args$infile <-
      "/rds/general/project/ratmann_pangea_deepsequencedata/live/200422_PANGEA2_MRC_alignment.fasta"
  }
}

# if prj.dir and out.dir are not manually set, default to here()
if (is.na(args$prj.dir))
{
  args$prj.dir <- here::here()
  args$out.dir <- here::here()
}

if (is.na(args$infile)) {
  stop("Please input the sequence file. ")
}
#
# add constants that should not be changed by the user
#
dir.create(args$out.dir)
out.dir <- file.path(args$out.dir, 'potential_network')
dir.create(out.dir)
max.per.run <- 4900

# Set seed
set.seed(args$seed)

cat(' ----------- Load data ---------- \n')
alignment <- read.fasta(file = args$infile)
nsequence <- length(alignment)
npos <- unique(lengths(alignment))

# define windows
windows_last_start <-
  ceiling(npos / args$sliding_width) * args$sliding_width - args$window_size + 1
windows_2last_end <-
  floor(npos / args$sliding_width) * args$sliding_width
windows_start <- seq(1, windows_last_start, args$sliding_width)
windows_end <-
  c(seq(args$window_size , windows_2last_end, args$sliding_width),
    npos)

cat(' ----------- Break sequences into windows ---------- \n')
for (i in 1:length(windows_start)) {
  subsequence <- list()
  subsequence.names <- vector()
  for (j in 1:nsequence) {
    subsequence[[j]] <-  alignment[[j]][windows_start[i]:windows_end[i]]
    subsequence.names <- c(subsequence.names, names(alignment)[j])
  }
  names(subsequence) <- subsequence.names
  write.fasta(
    sequences = subsequence,
    names = names(subsequence),
    nbchar = 60,
    file.out = file.path(out.dir , paste0('subsequence_window_', i, '.fasta'))
  )
}

cat(' ----------- Create batches of sequence pairs ---------- \n')
npair_perjob <- args$npair_perjob * length(windows_start)
sequences <- read.fasta(args$infile)
names_sequences <- names(sequences)
df <- data.table(t(combn(names_sequences, 2)))
colnames(df) <- c('H1', 'H2')
df[, SCRIPT_ID := rep(seq_len(ceiling(nrow(df) / npair_perjob)), each =
                        npair_perjob)[seq_len(nrow(df))]]
df[, ARRAY_ID := rep(seq_len(npair_perjob), times = ceiling(nrow(df) / npair_perjob))[seq_len(nrow(df))]]
tmp <- unique(subset(df, select = 'SCRIPT_ID'))
tmp[, JOB_ID := rep(seq_len(ceiling(nrow(tmp) / max.per.run)), each = max.per.run)[seq_len(nrow(tmp))]]
df <- merge(df, tmp, by = 'SCRIPT_ID')

cat('----------- Save job info ----------- \n')
if (args$if_save_data) {
  cat('\nWriting to file ', file.path(
    out.dir,
    paste0(
      'similarity_batches_windowsize_',
      args$window_size,
      'batchsize',
      args$npair_perjob,
      '.rds'
    )
  ))
  saveRDS(df, file = file.path(
    out.dir,
    paste0(
      'similarity_batches_windowsize_',
      args$window_size,
      'batchsize',
      args$npair_perjob,
      '.rds'
    )
  ))
  
}

df <- data.table(readRDS(file.path(
  out.dir,
  paste0(
    'similarity_batches_windowsize_',
    args$window_size,
    'batchsize',
    args$npair_perjob,
    '.rds'
  )
)))

cat('----------- Submit scripts for similarity calculation ----------- \n')
header <- paste0(
  "#!/bin/sh \n
#PBS -l walltime=",
  args$walltime,
  ":00:00,pcput=",
  args$walltime - 1,
  ":50:00 \n
#PBS -l select=1:ncpus=1:mem=",
  args$memory,
  "gb \n
#PBS -j oe \n"
)

for (jobid in seq_len(max(df$JOB_ID))) {
  cmd <- 'case $PBS_ARRAY_INDEX in\n'
  script_ids <- unique(df[JOB_ID == jobid,]$SCRIPT_ID)
  count <- 1
  for (scriptid in script_ids) {
    for (windowid in 1:length(windows_start)) {
      cmd <- paste0(
        cmd,
        count,
        ')\n',
        'Rscript calculate_similarity.R --script_id ',
        scriptid,
        ' --window_id ',
        windowid,
        ' --out_dir ',
        out.dir,
        '\n'
      )
      cmd <- paste0(cmd, ';; \n')
      count <- count + 1
    }
  }
  cmd <- paste0(cmd, 'esac')
  cmd <- paste0(
    header,
    '#PBS -J 1-',
    count,
    '\n',
    "module load anaconda3/personal \nsource activate ",
    args$env_name,
    '\n',
    'cd ',
    args$prj.dir,
    '\n',
    cmd
  )
  outfile <- file.path(out.dir,
                       paste0('script_calculate_similarity_job', jobid, '.qsub'))
  cat(cmd, file = outfile)
  cmd <- paste("qsub", outfile)
  cat(cmd)
  cat(system(cmd, intern = TRUE))
}
