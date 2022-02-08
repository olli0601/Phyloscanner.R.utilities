# The script for calculating similarity scores along the whole genome

# Preamble
# This script aims to calculate similarity scores defined in is.match over the whole genome for a batch of pairs.
# The inputs of the scripts include
# infile_sequence - the path to the consensus sequences.
# out_dir - the path to which outputs are stored.
# It should contain a file called similarity_batches.rds which indicates the pairs that belong to one batch.
# script_id - ID of the script that is processed.

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
    "--infile_sequence",
    type = "character",
    default = NA,
    help = "Absolute file path to sequence fasta file [default]",
    dest = "infile_sequence"
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
  )
)

args <-
  optparse::parse_args(optparse::OptionParser(option_list = option_list))

#
# test
#
if (0) {
  args <-
    list(infile_sequence = '/rds/general/project/ratmann_pangea_deepsequencedata/live/200422_PANGEA2_RCCSMRC_alignment.fasta',
         out.dir = '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/potential_network',
         script_id = 1L)
}



#
# use manually specified directories when args$out.dir is NA
#
tmp <- Sys.info()
if (tmp["user"] == "xx4515")
{
  if (is.na(args$out.dir))
  {
    args$out.dir <-'/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/potential_network'
  }
  if (is.na(args$infile_sequence))
  {
    args$infile_sequence <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/200422_PANGEA2_RCCSMRC_alignment.fasta'
  }
}

# if infile sequence and out.dir are not manually set, default to here()
if (is.na(args$out.dir))
{
  args$out.dir <- here::here()
  args$infile_sequence <-
    file.path(here::here(), '200422_PANGEA2_RCCSMRC_alignment.fasta')
}

cat(' ---------- Load job information ---------- ')
infile_job_info <- file.path(args$out.dir, 'similarity_batches.rds')
dfo <- data.table(readRDS(infile_job_info))


cat(' ---------- Load sequences ---------- ')
sq <- read.fasta(args$infile_sequence)

cat(' ---------- Select pairs and sequences ---------- ')
dfo	<- subset(dfo, SCRIPT_ID == args$script_id)
sq	<- sq[dfo[, unique(c(H1, H2))]]

cat(' ---------- Calculate similarities ---------- ')
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
}, by = c('ARRAY_ID', 'H1', 'H2')]
dfo[, PERC := MATCHES / LENGTH]

cat(' ---------- Save ---------- ')
saveRDS(dfo, file = file.path(args$out.dir, paste0('similarity', args$script_id, '.rds')))
