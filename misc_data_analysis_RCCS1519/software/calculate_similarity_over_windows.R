cat("\n\n===== calculate_similarity_over_windows.R =====\n\n")
# The phyloscanner run - calculate similarities between the consensus sequences on the overlapping windows

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
    dest = 'pkg.dir'
  ),
  optparse::make_option(
    "--out_dir_base",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where all output is stored [default]",
    dest = 'out.dir'
  )
)

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

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
    pkg.dir = NA,
    infile = NA
  )
}

#
# use manually specified directories when args$out.dir is NA
#
tmp <- Sys.info()

.f <- function(x, path) 
        fifelse(is.na(x), path, x)

if (tmp["user"] == "xx4515")
{
        args$out.dir <- .f(args$out.dir, "/rds/general/project/ratmann_deepseq_analyses/live/test")
        args$pkg.dir <- .f(args$pkg.dir, "~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software")
        args$infile  <- .f(args$infile, "/rds/general/project/ratmann_pangea_deepsequencedata/live/200422_PANGEA2_MRC_alignment.fasta")
}

if(0)
{
        args$out.dir <- .f(args$out.dir, "~/Dropbox/CONDESAPopSample_phsc_stage1_output_close")
        args$pkg.dir <- .f(args$pkg.dir, "~/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software")
        args$infile <- .f(args$infile, "~/Dropbox/CONDESAPopSample_phsc_stage1_output_close/_condesa_wgs.nextalign.fas")
}

if (is.na(args$infile)) 
        stop("Please input the sequence file. ")

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
with(args,{
        windows_last_start <<- ceiling(npos / sliding_width) * sliding_width - window_size + 1
        windows_2last_end <<- floor(npos / sliding_width) * sliding_width
        windows_start <<- seq(1, windows_last_start, sliding_width)
        windows_end <<- c(seq(window_size , windows_2last_end, sliding_width), npos)
  })

cat(' ----------- Break sequences into windows ---------- \n')

subset.consensus.seqs.into.windows <- function(i)
{
  .extract.ith.window <- function(x, i)
          x[ windows_start[i]:windows_end[i] ]
  subsequence <- lapply(alignment, .extract.ith.window, i=i)

  write.fasta(
    sequences = subsequence, names = names(subsequence),  nbchar = 60,
    file.out = file.path(out.dir , paste0('subsequence_window_', i, '.fasta'))
  )
  return(TRUE)
}
out <- lapply(1:length(windows_start), subset.consensus.seqs.into.windows)
stopifnot( all( unlist(out) ) )

cat(' ----------- Create batches of sequence pairs ---------- \n')
with(args, {
        npair_perjob <<- npair_perjob * length(windows_start)
        sequences <<- read.fasta(infile)
        names_sequences <<- names(sequences)
  })
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

        filename <- file.path(out.dir, paste0( 'similarity_batches_windowsize_', args$window_size, 'batchsize', args$npair_perjob, '.rds'))
        cat('\nWriting to file ', filename, '\n')
        saveRDS(df, file = filename) 
}

cat('----------- Submit scripts for similarity calculation ----------- \n')

# Source functions
source(file.path(args$pkg.dir, "utility.R"))

# DT <- df[JOB_ID == 1]
.write.job <- function(DT, jobid)
{
        scriptids <- DT[ , unique(SCRIPT_ID)]

        # for each script id, for each window, calculate similarities
        path <- file.path(args$pkg.dir, 'calculate_similarity.R')

       .f <- function(w, s)
               paste0('Rscript ', path,' --script_id ', s,' --window_id ', w, ' --out_dir ', out.dir, '\n')

       .g <- function(s)
               lapply(1:length(windows_start), .f, s=s)

       tmp <- lapply(script_ids, .g)
       tmp <- unlist(tmp)
       J <- length(tmp)

       # cmd <- lapply(1:length(windows_start), .f)
       cmd <- lapply(seq_along(tmp), function(i) paste0(i, ')\n', tmp[i]))
       cmd <- paste0(cmd, collapse=';;\n')
       # cat(cmd)
        
       # Add headers
        pbshead <- cmd.hpcwrapper.cx1.ic.ac.uk(
                                               hpc.select = 1,
                                               hpc.nproc = 1,
                                               hpc.walltime = args$walltime,
                                               hpc.mem = args$memory,
                                               hpc.array = J,
                                               hpc.load = paste0("module load anaconda3/personal source activate ", args$env_name))

       cmd <- paste0(pbshead, '\ncase $PBS_ARRAY_INDEX in\n', cmd, '\n;;\nesac')

       outfile <- file.path(out.dir, paste0('script_calculate_similarity_job', unique(jobid), '.sh'))
       cat(cmd, file = outfile)
       outfile
}

djob <- df[, list(PATH=.write.job(.SD, jobid=JOB_ID)), by=JOB_ID]

.qsub <- function(path)
{
        cmd <- paste("qsub", path)
        cat(cmd)
        cat(system(cmd, intern = TRUE))
}
djob[, .qsub(PATH), by=PATH]
