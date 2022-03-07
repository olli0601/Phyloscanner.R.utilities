cat('\n\n=====  TSIcalculate_similarity_over_genome.R =====\n\n')
# The phyloscanner run - calculate similarities between consensus sequences
# on the whole genome

# Preamble
# The set of scripts aims to run phyloscanner.
# This script aims to generate scripts and
# calculate the similarities between consensus sequences
# which are stored in infile over the whole genome.

# Load the required packages
library(data.table)
library(seqinr)
# library(optparse)
# library(here)
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
    default = 24L,
    help = "Job walltime (hours) [default]",
    dest = 'walltime'
  ),
  optparse::make_option(
    "--memory",
    type = "integer",
    default = 10L,
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
  ),
  optparse::make_option(
    "--controller",
    type = "character",
    default = NA_character_, # Think about adding the controller in the software directory
    help = "Path to sh script directing the full analysis",
    dest = 'controller'
  )
)

args <-
  optparse::parse_args(optparse::OptionParser(option_list = option_list))

args <- list(controller=NA)
if(is.na(args$controller))
{
        args$controller <- file.path(args$prj.dir, args$controller)
}

#
# test
#
if (0) {
  args <- list(
    verbose = T,
    seed = 42L,
    npair_perjob = 10000L,
    env_name = "phyloenv",
    walltime = 24L,
    memory = 10L,
    if_save_data = T,
    out.dir = NA,
    prj.dir = NA,
    infile = NA,
    controller=NA
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
      "/rds/general/project/ratmann_deepseq_analyses/live/testTSI/"
  }
  if (is.na(args$prj.dir))
  {
    args$prj.dir <-
      "~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/"
  }
  if (is.na(args$infile)) {
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

cat('----------- Load sequences ----------- \n')
sequences <- read.fasta(args$infile)
names_sequences <- names(sequences)
df <- data.table(t(combn(names_sequences, 2)))
colnames(df) <- c('H1', 'H2')
df[, SCRIPT_ID := rep(seq_len(ceiling(nrow(df) / args$npair_perjob)), each =
                        args$npair_perjob)[seq_len(nrow(df))]]
df[, ARRAY_ID := rep(seq_len(args$npair_perjob), times = ceiling(nrow(df) /
                                                                   args$npair_perjob))[seq_len(nrow(df))]]
tmp <- unique(subset(df, select = 'SCRIPT_ID'))
tmp[, JOB_ID := rep(seq_len(ceiling(nrow(tmp) / max.per.run)), each = max.per.run)[seq_len(nrow(tmp))]]
df <- merge(df, tmp, by = 'SCRIPT_ID')

cat('----------- Save job info ----------- \n')
saveRDS(df, file = file.path(out.dir, 'similarity_batches.rds'))

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

if(file.exists(args$controller))
{
        footer <- paste0('
               OUT="',out.dir,'"
               N_run=$(ls $OUT/similarity[0-9]*.rds | wc -l)
               N_tot=',max(df$JOB_ID),' 

               # queue next command
               if [ $N_run -ge $N_tot ];then
               qsub -v STEP="net" ',args$controller,'
               fi'
        ) 
}


for (jobid in seq_len(max(df$JOB_ID))) {
  cmd <- 'case $PBS_ARRAY_INDEX in\n'
  script_ids <- unique(df[JOB_ID == jobid,]$SCRIPT_ID)
  count <- 1
  for (scriptid in script_ids) {
    cmd <- paste0(
      cmd,
      count,
      ')\n',
      # While loop to check whether output file is produced: rerun at most 10 times.
      'n=0\n',
      paste0('while [ "$n" -lt 10 ] && [ -f similarity', as.character(scriptid),'.rds ]; do'),
      'n=$(( n + 1 ))\n',
      'Rscript ',args$prj.dir,'/calculate_similarity.R --script_id ',
      scriptid,
      ' --out_dir ',
      out.dir,
      ' --infile ',
      args$infile,
      '\ndone',
      '\n;; \n'


    )
    count <- count + 1
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
    cmd
  )
  cmd <- paste0(cmd, '\n',
                'cp * ', out.dir)

  outfile <- file.path(out.dir,
                       paste0('script_calculate_similarity_job', jobid, '.qsub'))
  cat(cmd, file = outfile)
  cmd <- paste("qsub", outfile)
  cat(cmd)
  cat(system(cmd, intern = TRUE))
}

# in the above command, could add a function which checks whether all the runs have run!
# So by checking whether all of these files exists in the output directory!!!!!! 

