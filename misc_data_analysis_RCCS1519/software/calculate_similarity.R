cat("\n\n===== calculate_similarity.R =====\n\n")
# The script for calculating similarity scores along the whole genome

# Preamble
# This script aims to calculate similarity scores
# defined by is.match for a batch of sequence pairs.
# The inputs of the scripts include
# infile - the path to the consensus sequences.
# out_dir - the path to which outputs are stored.
# It should contain a file called similarity_batches.rds
# which divides the pairs into batches.
# script_id - ID of the script.
# window_id - ID of the window. If NA, consider the whole genome.

# Load the required packages
library(data.table)
library(seqinr)

# Define functions just used in this script alone
is.match <- function(x, y) {
  if (x %in% c('a', 'g', 'c', 't') & y %in% c('a', 'g', 'c', 't')) {
    if (x == y) {
      ans = 1
    } else{
      ans = 0
    }
  } else{
    if (((x == 'g' | x == 'a') & y == 'r') |
        ((y == 'g' | y == 'a') & x == 'r') |
        ((x == 't' | x == 'c') & y == 'y') |
        ((y == 't' | y == 'c') & x == 'y') |
        ((x == 'a' | x == 't') & y == 'w') |
        ((y == 'a' | y == 't') & x == 'w') |
        ((x == 'a' | x == 'c') & y == 'm') |
        ((y == 'a' | y == 'c') & x == 'm') |
        ((x == 'g' | x == 'c') & y == 's') |
        ((y == 'g' | y == 'c') & x == 's') |
        ((x == 'g' | x == 't') & y == 'k') |
        ((y == 'g' | y == 't') & x == 'k')) {
      ans = 1 / 2
    } else if (((x == 'a' | x == 'c' | x == 't') & y == 'h') |
               ((y == 'a' | y == 'c' | y == 't') & x == 'h') |
               ((x == 'a' | x == 'g' | x == 't') & y == 'd') |
               ((y == 'a' | y == 'g' | y == 't') & x == 'd') |
               ((x == 'g' | x == 'c' | x == 't') & y == 'b') |
               ((y == 'g' | y == 'c' | y == 't') & x == 'b') |
               ((x == 'a' | x == 'c' | x == 'g') & y == 'v') |
               ((y == 'a' | y == 'c' | y == 'g') & x == 'v')) {
      ans = 1 / 3
    } else{
      ans = 0
    }
  }
  return(ans)
}

#
# Define input arguments that can be changed by users
#
option_list <- list(
  optparse::make_option(
    "--infile",
    type = "character",
    default = NA,
    help = "Absolute file path to sequence fasta file [default]",
    dest = "infile"
  ),
  optparse::make_option(
    "--out_dir",
    type = "character",
    default = NA,
    help = "Absolute file path to base directory where all output is stored [default]",
    dest = "out.dir"
  ),
  optparse::make_option(
    "--script_id",
    type = "integer",
    default = 1L,
    help = "Current script ID [default]",
    dest = "script_id"
  ),
  optparse::make_option(
    "--window_id",
    type = "integer",
    default = NA_integer_,
    help = "Current window ID [default]",
    dest = "window_id"
  ),
  optparse::make_option(
    "--save_data",
    action = "store_true",
    default = TRUE,
    help = "Save data [default]",
    dest = 'if_save_data'
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
  )
)

args <-
  optparse::parse_args(optparse::OptionParser(option_list = option_list))

#
# test
#
if (0) {
  args <-
    list(
      infile = '/rds/general/project/ratmann_pangea_deepsequencedata/live/200422_PANGEA2_MRC_alignment.fasta',
      out.dir = '/rds/general/project/ratmann_deepseq_analyses/live/testTSI/potential_network/',
      script_id = 1L,
      window_id = NA,
      npair_perjob = 1e4,
      window_size = 500
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
      '/rds/general/project/ratmann_deepseq_analyses/live/testTSI/potential_network'
  }
  if (is.na(args$infile))
  {
    args$infile <-
      '/rds/general/project/ratmann_pangea_deepsequencedata/live/200422_PANGEA2_MRC_alignment.fasta'
  }
}

# if infile sequence and out.dir are not manually set, default to here()
if (is.na(args$out.dir))
{
  args$out.dir <- here::here()
}
if (!is.na(args$window_id)) {
  args$infile <-
    file.path(args$out.dir,
              paste0('subsequence_window_', args$window_id, '.fasta'))
}
if (is.na(args$infile)) {
  stop("Please input the sequence file. ")
}

cat(' ---------- Load job information ---------- \n')
if (is.na(args$window_id)) {
  infile_job_info <- file.path(args$out.dir, 'similarity_batches.rds')
} else{
  infile_job_info <- file.path(
    args$out.dir,
    paste0(
      'similarity_batches_windowsize_',
      args$window_size,
      'batchsize',
      args$npair_perjob,
      '.rds'
    )
  )
}

dfo <- data.table(readRDS(infile_job_info))


cat(' ---------- Load sequences ---------- \n')
sq <- read.fasta(args$infile)

cat(' ---------- Select pairs and sequences ---------- \n')
dfo	<- subset(dfo, SCRIPT_ID == args$script_id)
sq	<- sq[dfo[, unique(c(H1, H2))]]

# remove sequence with all gaps
sq		<- sq[dfo[, unique(c(H1, H2))]]
sq.n <-
  unlist(lapply(sq, function(x) {
    sum(x == '?' | x == '-' | x == 'n') == length(x)
  }))
sq.n.names <- names(sq.n[sq.n == TRUE])
tmp <- dfo[H1 %in% sq.n.names | H2 %in% sq.n.names]
tmp[, LENGTH := -1L]
tmp[, MATCHES := -1.0]
tmp[, PERC := -1.0]
dfo <- dfo[!H1 %in% sq.n.names & !H2 %in% sq.n.names]

cat(' ---------- Calculate similarities ---------- \n')
dfo		<- dfo[, {
  seq1 <- sq[[H1]]
  seq2 <- sq[[H2]]
  pos <-
    which(seq1 != '-' &
            seq1 != '?' & seq1 != 'n' & seq2 != '-' & seq2 != '?' &
            seq2 != 'n')
  if (length(pos) == 0) {
    len <- -1L
    match <- -1.0
  } else{
    seq1 <- seq1[pos]
    seq2 <- seq2[pos]
    len <- length(pos)
    match <-
      sum(sapply(1:length(pos), function(pos) {
        is.match(seq1[pos], seq2[pos])
      }))
  }
  list(LENGTH = len,
       MATCHES = match)
}, by = c('ARRAY_ID', 'SCRIPT_ID','JOB_ID','H1', 'H2')]
dfo[, PERC := MATCHES / LENGTH]
dfo[LENGTH == -1L, PERC := -1.0]
dfo <- rbind(dfo, tmp)

cat(' ---------- Save ---------- \n')
if (args$if_save_data) {
  if (is.na(args$window_id)) {
    cat('\nWriting to file ',
        file.path(
          args$out.dir,
          paste0('similarity', args$script_id, '.rds')
        ),
        '...\n')
    saveRDS(dfo, file = file.path(
      args$out.dir,
      paste0('similarity', args$script_id, '.rds')
    ))
  } else{
    cat('\nWriting to file ', file.path(
      args$out.dir,
      paste0(
        'similarity',
        args$script_id,
        '_window_',
        args$window_id,
        '.rds'
      )
    ))
    saveRDS(dfo, file = file.path(
      args$out.dir,
      paste0(
        'similarity',
        args$script_id,
        '_window_',
        args$window_id,
        '.rds'
      )
    ))
  }
}
