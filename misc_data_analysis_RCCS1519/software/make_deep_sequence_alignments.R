# The phyloscanner run - make alignments

# Preamble
# The set of scripts aims to run phyloscanner.
# This script aims to make alignments in each potential transmission network.

# Load the required packages
library(data.table)
library(seqinr)
library(tidyverse)
library(dplyr)

#
# Define input arguments that can be changed by users
#
option_list <- list(
  optparse::make_option(
    c("-v", "--verbose"),
    action = "store_true",
    default = TRUE,
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
    "--sliding_width",
    type = "integer",
    default = 10L,
    help = "Sliding width [default %default]",
    dest = "sliding_width"
  ),
  optparse::make_option(
    "--window_size",
    type = "integer",
    default = 250L,
    help = "Window size [default %default]",
    dest = "window_size"
  ),
  optparse::make_option(
    "--window_cutoff",
    type = "numeric",
    default = NA,
    help = "Cutoff of proportion of windows indicating close relationship [default %default]",
    dest = "window_cutoff"
  ),
  optparse::make_option(
    "--n_control",
    type = "integer",
    default = 0,
    help = "Number of controls added [default %default]",
    dest = "n_control"
  ),
  optparse::make_option(
    "--cluster_size",
    type = "integer",
    default = 100L,
    help = "Maximum cluster size [default %default]",
    dest = "cluster_size"
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
    "--prog_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where phyloscanner is stored [default]",
    dest = 'prog.dir'
  ),
  optparse::make_option(
    "--date",
    type = 'character',
    default = as.character(Sys.Date()),
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
if(0){
  args <- list(
    verbose = T,
    seed = 42,
    sliding_width = 10L,
    window_size = 250L,
    window_cutoff = NA,
    n_control = 0,
    cluster_size = 100,
    if_save_data = T,
    date = '2022-02-10',
    out.dir = NA,
    prj.dir = NA,
    prog.dir = NA
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
      "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/"
  }
  if (is.na(args$prj.dir))
  {
    args$prj.dir <-
      "~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/"
  }
  if (is.na(args$prog.dir))
  {
    args$prog.dir <-
      "~/phyloscanner"
  }
}

# if prj.dir and out.dir are not manually set, default to here()
if (is.na(args$prj.dir))
{
  args$prj.dir <- here::here()
  args$out.dir <- here::here()
  args$prog.dir <- here::here()
}

#
# Add constants that should not be changed by the user
#
dir.data <-
  '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
dir.net <-
  file.path(args$out.dir, "potential_network")

if(!is.na(args$window_cutoff)){
  infile.runs <- file.path(
    args$out.dir,
    paste0(
      'phscinput_runs_clusize_',
      args$cluster_size,
      '_ncontrol_',
      args$n_control,
      '_windowcutoff_',
      args$window_cutoff,
      '.rds'
    ))
}else{
  infile.runs <- file.path(
    args$out.dir,
    paste0(
      'phscinput_runs_clusize_',
      args$cluster_size,
      '_ncontrol_',
      args$n_control,
      '.rds'
    ))
}

max.per.run <- 4900

args$date <- gsub('-','_',args$date)
# Set default output directories relative to out.dir
args$out.dir.data <-
  file.path(args$out.dir, paste0(args$date, "_phsc_input"))
args$out.dir.work <-
  file.path(args$out.dir, paste0(args$date, "_phsc_work"))
args$out.dir.output <-
  file.path(args$out.dir, paste0(args$date, "_phsc_output"))

# Create directories if needed
ifelse(!dir.exists(args$out.dir.data),
       dir.create(args$out.dir.data),
       FALSE)
ifelse(!dir.exists(args$out.dir.work),
       dir.create(args$out.dir.work),
       FALSE)
ifelse(!dir.exists(args$out.dir.output),
       dir.create(args$out.dir.output),
       FALSE)

# Copy files into input folder
file.copy(
  list.files(file.path('/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/', '210325_phsc_input'), full.names = T),
  args$out.dir.data,
  overwrite = T,
  recursive = FALSE,
  copy.mode = TRUE
)

infile.consensus <-
  file.path(args$out.dir.data, 'ConsensusGenomes.fasta')
infile.consensus.oneeach <-
  file.path(args$out.dir.data,
            '2019_New_ConsensusGenomesOneEach_GeneCut.fasta')


# Source functions
source(file.path(args$prj.dir, "utility.R"))

# Check duplicates
pty.runs <- data.table(readRDS(infile.runs))
if ('ID_TYPE' %in% colnames(pty.runs)) {
  setorder(pty.runs, PTY_RUN, -ID_TYPE, UNIT_ID)
} else{
  setorder(pty.runs, PTY_RUN, UNIT_ID)
}

tmp <- pty.runs[, duplicated(SAMPLE_ID), by = 'PTY_RUN']
tmp <- tmp[, which(V1)]
if(length(tmp)!=0){
  pty.runs <- pty.runs[-tmp, ]
}

tmp <-
  pty.runs[, length(SAMPLE_ID) - length(unique(SAMPLE_ID)), by = 'PTY_RUN']
stopifnot(all(tmp$V1 == 0))

# Load backgrounds
consensus_seq <- seqinr::read.fasta(infile.consensus)
consensus_seq_names <- names(consensus_seq)
hxb2 <- grep('HXB2', names(consensus_seq), value = T)
hxb2_seq <- consensus_seq[[hxb2]]
root_seq <-
  grep(
    '^REF_CON_M$|REF_CONSENSUS_M$|^REF_CON_H$|REF_CONSENSUS_H$',
    names(consensus_seq),
    value = T
  )

# Change the format
if ('ID_TYPE' %in% colnames(pty.runs)) {
  pty.runs[ID_TYPE == 'control', UNIT_ID := paste0('CNTRL-', UNIT_ID)]
  pty.runs[ID_TYPE == 'control', RENAME_ID := paste0('CNTRL-', RENAME_ID)]
}
pty.runs[, BAM := paste0(dir.data, SAMPLE_ID, '.bam')]
pty.runs[, REF := paste0(dir.data, SAMPLE_ID, '_ref.fasta')]
setkey(pty.runs, PTY_RUN, RENAME_ID)

# Remove starts, ends and vloops
ptyi <- seq(800, 9175, args$sliding_width)
ptyi <- c(ptyi[ptyi <= 6615 - args$window_size], 6825, 6850, ptyi[ptyi >= 7636])

pty.c	<- lapply(seq_along(ptyi), function(i)
{
  pty.args <- list(
    prog.pty = file.path(args$prog.dir, "phyloscanner_make_trees.py"),
    prog.mafft = '\" mafft --globalpair --maxiterate 1000 \" ',
    data.dir = args$out.dir.data,
    work.dir = args$out.dir.work,
    out.dir = args$out.dir.output,
    alignments.file = infile.consensus,
    alignments.root = root_seq,
    alignments.pairwise.to = hxb2,
    window.automatic = '',
    merge.threshold = 0,
    min.read.count = 1,
    quality.trim.ends = 23,
    min.internal.quality = 23,
    merge.paired.reads = TRUE,
    discard.improper.pairs = TRUE,
    no.trees = TRUE,
    dont.check.duplicates = FALSE,
    dont.check.recombination = TRUE,
    num.bootstraps = 1,
    all.bootstrap.trees = TRUE,
    strip.max.len = 350,
    min.ureads.individual = NA,
    win = c(ptyi[i], ptyi[i] + args$window_size, args$sliding_width, args$window_size),
    keep.overhangs = FALSE,
    mem.save = 0,
    verbose = TRUE,
    select = NA,
    default.coord = TRUE,
    realignment = TRUE
  )
  pty.c <- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)
  pty.c[, W_FROM := ptyi[i]]
  pty.c
})
pty.c	<- do.call('rbind', pty.c)
setkey(pty.c, PTY_RUN, W_FROM)

print(pty.c)

pty.c[, CASE_ID := rep(1:max.per.run, times = ceiling(nrow(pty.c) / max.per.run))[1:nrow(pty.c)]]
pty.c[, JOB_ID := rep(1:ceiling(nrow(pty.c) / max.per.run), each = max.per.run)[1:nrow(pty.c)]]

#	Define PBS variables
hpc.load			<-
  "module load intel-suite/2015.1 mpi raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"	# make third party requirements available
hpc.select			<- 1
hpc.nproc			<- 1
hpc.walltime		<- 71
hpc.q				<- NA
hpc.mem				<- "6gb"
hpc.array			<- pty.c[, max(CASE_ID)]

print(hpc.array)
#	Define PBS header for job scheduler
pbshead		<- "#!/bin/sh"
tmp			<-
  paste("#PBS -l walltime=",
        hpc.walltime,
        ":59:00,pcput=",
        hpc.walltime,
        ":45:00",
        sep = "")
pbshead		<- paste(pbshead, tmp, sep = "\n")
tmp			<-
  paste("#PBS -l select=",
        hpc.select,
        ":ncpus=",
        hpc.nproc,
        ":mem=",
        hpc.mem,
        sep = "")
pbshead 	<- paste(pbshead, tmp, sep = "\n")
pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")
if (!is.na(hpc.array))
  pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep = '')
if (!is.na(hpc.q))
  pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
pbshead 	<- paste(pbshead, hpc.load, sep = "\n")
cat(pbshead)

print(pty.c[, max(JOB_ID)])
print(max(pty.c$JOB_ID))
#	Create PBS job array
for (i in 1:pty.c[, max(JOB_ID)]) {
  tmp <- pty.c[JOB_ID == i, ]
  cmd <-
    tmp[, list(CASE = paste0(CASE_ID, ')\n', CMD, ';;\n')), by = 'CASE_ID']
  cmd <-
    cmd[, paste0('case $PBS_ARRAY_INDEX in\n',
                 paste0(CASE, collapse = ''),
                 'esac')]
  cmd <- paste(pbshead, cmd, sep = '\n')
  outfile <-
    gsub(':', '', paste(
      "readali",
      paste0('job', i),
      paste(
        strsplit(date(), split = ' ')[[1]],
        collapse = '_',
        sep = ''
      ),
      'sh',
      sep = '.'
    ))
  outfile <- file.path(args$out.dir.work, outfile)
  cat(cmd, file = outfile)
  cmd <- paste("cd ",dirname(outfile),'\n',"qsub ", outfile)
  cat(cmd)
  cat(system(cmd, intern = TRUE))
}
