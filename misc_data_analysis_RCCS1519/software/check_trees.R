# The phyloscanner run for TSI estimates - make trees if not exist

# Preamble
# This script aims to make trees using each alignment 
# if the trees were not built in make_trees.R due to the time constrain.

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
    "--save_data",
    action = "store_true",
    default = TRUE,
    help = "Save data [default]",
    dest = 'if_save_data'
  ),
  optparse::make_option(
    "--env_name",
    type = "character",
    default = 'phylor4',
    help = "Conda environment name [default]",
    dest = 'env_name'
  ),
  optparse::make_option(
    "--iqtree_option",
    type = "character",
    default = 'GTR+F+R6',
    help = "Options in IQTREE [default]",
    dest = 'iqtree_option'
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
    default = '2022-02-04',
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
    if_save_data = T,
    date = '2022-02-04',
    env_name = 'phylor4',
    iqtree_option = '-m GTR+F+R6 -o REF_CON_H',
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
      "~/phyloscanner/phyloscanner_make_trees.py"
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
max.per.run <- 4900

# Set default output directories relative to out.dir
args$out.dir.data <-
  file.path(args$out.dir, paste0(args$date, "_phsc_input"))
args$out.dir.work <-
  file.path(args$out.dir, paste0(args$date, "_phsc_work"))
args$out.dir.output <-
  file.path(args$out.dir, paste0(args$date, "_phsc_output"))

# Source functions
source(file.path(args$prj.dir, "utility.R"))
#
#	produce trees
#
if(0)	
{
  hpc.select<- 1; hpc.nproc<- 1; hpc.walltime<- 4; hpc.mem<- "1850mb"; hpc.q<- NA
}
if(0)	
{
  hpc.select<- 1; hpc.nproc<- 1; hpc.walltime<- 23; hpc.mem<- "1850mb"; hpc.q<- NA
}
if(1)	
{
  hpc.select<- 1; hpc.nproc<- 1; hpc.walltime<- 71; hpc.mem<- "63850mb"; hpc.q<- NA
}

iqtree.pr <- 'iqtree'
iqtree.pr <- 'iqtree'
if(!is.na(args$seed)){
  iqtree.args	<- paste0(args$iqtree_option, ' -seed ', args$seed)
}else{
  iqtree.args	<- args$iqtree_option
}

# Load alignments
infiles	<- data.table(FI=list.files(args$out.dir.output, pattern='_v2.fasta$', full.names=TRUE, recursive=TRUE))
infiles[, FO:= gsub('.fasta$','',FI)]
infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FI)))]
infiles[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FI)))]		
infiles[is.na(W_FROM),W_FROM:= as.integer(gsub('.*PositionsExcised_([0-9]+)_.*','\\1',basename(FI)))]
infiles[,PositionsExcised:=grepl('PositionsExcised',FI)]
setkey(infiles, PTY_RUN, W_FROM)
tmp <- infiles[,list(NUM=length(PositionsExcised)),by=c('PTY_RUN', 'W_FROM')]
infiles <- merge(infiles, tmp, by=c('PTY_RUN', 'W_FROM'))
infiles <- infiles[!(PositionsExcised==F & NUM==2),]
infiles[,FO_NAME:=paste0(FO,'.treefile')]
infiles[,FO_EXIST:=file.exists(FO_NAME)]
infiles <- infiles[FO_EXIST==F]
if(args$verbose){
  cat('These trees were not built sucessfully...\n')
  print(head(infiles))
  cat('Build again...\n')
}

# Set up jobs
df<- infiles[, list(CMD=cmd.iqtree(FI, outfile=FO, pr=iqtree.pr, pr.args=iqtree.args)), by=c('PTY_RUN','W_FROM')]
df[, ID:=ceiling(seq_len(nrow(df))/4)]
df<- df[, list(CMD=paste(CMD, collapse='\n',sep='')), by='ID']

#	Create PBS job array
if(nrow(df) > max.per.run){
  df[, CASE_ID:= rep(1:max.per.run,times=ceiling(nrow(df)/max.per.run))[1:nrow(df)]]
  df[, JOB_ID:= rep(1:ceiling(nrow(df)/max.per.run),each=max.per.run)[1:nrow(df)]]
  df[,ID:=NULL]
}else{
  setnames(df, 'ID', 'CASE_ID')
  df[, JOB_ID:=1]
}

pbshead	<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=NULL)

for (i in 1:df[, max(JOB_ID)]) {
  tmp<-df[JOB_ID==i,]
  hpc.array <- nrow(tmp)
  cmd<-tmp[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
  cmd<-cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]		
  tmp <- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')	
  tmp <- paste0(tmp,'\n module load anaconda3/personal \n source activate ',args$env_name)
  cmd<-paste(tmp,cmd,sep='\n')	
  
  #	Submit 
  outfile	<- paste("srx",paste0('job',i),paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.')
  outfile <- gsub(':','_',outfile)
  outfile<-file.path(args$out.dir.work, outfile)
  cat(cmd, file=outfile)
  cmd<-paste("qsub", outfile)
  cat(cmd)
  cat(system(cmd, intern= TRUE))
}