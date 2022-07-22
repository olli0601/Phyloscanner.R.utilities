cat("\n\n===== make_trees.R =====\n\n")
# The phyloscanner run - make trees

# Preamble
# This script aims to make trees using each alignment.

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
    "--iqtree_method",
    type = "character",
    default = "GTR+F+R6",
    help = "Methods in IQTREE [default]",
    dest = 'iqtree_method'
  ),
  optparse::make_option(
    "--iqtree_root",
    type = "character",
    default = NULL,
    help = "Root in IQTREE [default]",
    dest = 'iqtree_root'
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
  ),
  optparse::make_option(
    "--controller",
    type = "character",
    default = NA_character_, 
    help = "Path to sh script directing the full analysis",
    dest = 'controller'
  ),
  optparse::make_option(
    "--date",
    type = 'character',
    default = as.character(Sys.Date()),
    metavar = '"YYYY-MM-DD"',
    help = 'As of date to extract data from.  Defaults to today.',
    dest = 'date'
  ),
  optparse::make_option(
    "--walltime_idx",
    type = "integer",
    default = 2,
    help = "Indicator for amount of resources required by job. Values ranging from 1 (lala) to 3 (lala)",
    dest = "walltime_idx"
  )
)

args <-  optparse::parse_args(optparse::OptionParser(option_list = option_list))

#
# Helpers
#

.check.alignments <- function(infiles)
{
  check <- list.files(args$out.dir.output, 
                      pattern='InWindow_(.*)?[0-9]{2}.fasta$',
                      full.names=TRUE, recursive=TRUE)
  check <- gsub('.fasta', '_v2.fasta', check)
  
  check <- round(mean(check %in% infiles)*100, 2)
  if( check == 100){
    cat('For each v1 align, there is a v2_align')
  }else{
    cat('Only', check, '% of v1 alignments have a corresponding v2 alignment')
  }
}

.check.or.resubmit.incompleted.alignments <- function()
{
  dlogs <- classify.log.errors.and.successes(args$out.dir.work)
  dlogs[, F :=  gsub( paste0(args$out.dir.work,'/'), '', F)]
  tmp <- dlogs[ done == 0, .(conda, bam, algRds, kill, F)]
  tmp[ conda == 0]
  
  daligns <- tmp[algRds == 0] 
  daligns[, SH := paste0(dirname(F), '.sh')]
  daligns[, PBS := as.numeric(   gsub('^.*\\.([0-9]+)$', '\\1', F) )]
  
  cmd <- rewrite_job(daligns, double_walltime = TRUE)
  if(str_count(cmd, ';;') > 20)
  {
    # store in 'readali2'-prefixed .sh files
    time <- paste0(gsub(':', '', strsplit(date(), split = ' ')[[1]]), collapse='_')
    outfile <- paste("readali2", time, 'sh', sep='.')
    outfile <- file.path(args$out.dir.work, outfile)
    
    # change to work directory and submit to queue
    cat(cmd, file = outfile)
    cmd <- paste0("cd ",dirname(outfile),'\n',"qsub ", outfile)
    cat(cmd)
    cat(system(cmd, intern = TRUE))
    stop('Submitted realignment')
  }
}

.make.iqtree.opt <- function(args)
{
  iqtree.pr <<- 'iqtree'
  iqtree.args <<- paste0('-m ',args$iqtree_method)
  
  if(!is.null(args$iqtree_root))
    iqtree.args	<<- aste0(iqtree.args, ' -o ', args$iqtree_root)
  
  if(!is.na(args$seed))
    iqtree.args	<<- paste0(iqtree.args, ' -seed ', args$seed)
}

#
# test
#
if(0){
  args <- list(
    out.dir="/rds/general/project/ratmann_deepseq_analyses/live/seroconverters3_alignXX",
    pkg.dir="/rds/general/user/ab1820/home/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software",
    iqtree_method="GTR+F+R6",
    env_name = 'phylostan',
    date = '2022-07-22',
    seed = 42,
    walltime_idx=1,
    controller="/rds/general/user/ab1820/home/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/runall_TSI_seroconv3.sh"
  )
}

# 
# Chose PBS specifications according to PBS index
# Idea is that this script can be run multiple times with increasing res reqs

list(
  `1`= list(hpc.select = 1, hpc.nproc = 1, hpc.walltime = 4 , hpc.mem = "2gb" ,hpc.q = NA, max.per.run=950),
  `2`= list(hpc.select = 1, hpc.nproc = 1, hpc.walltime = 23, hpc.mem = "2gb" ,hpc.q = NA, max.per.run=475),
  `3`= list(hpc.select = 1, hpc.nproc = 1, hpc.walltime = 71, hpc.mem = "63gb",hpc.q = NA, max.per.run=475)
) -> pbs_headers

tmp <- pbs_headers[[args$walltime_idx]]
cat('selected the following PBS specifications:\n')
print(tmp)
invisible(list2env(tmp,globalenv()))

#
#       Paths + cleaning
#

# Set default output directories relative to out.dir
args$date <- gsub('-','_',args$date)
.f <- function(x)
{
        out <- file.path(args$out.dir, paste0(args$date, x))
        stopifnot(file.exists(out))
        out
}
args$out.dir.data <- .f('_phsc_input')
args$out.dir.work <- .f('_phsc_work')
args$out.dir.output <- .f('_phsc_output')

# Source functions
source(file.path(args$pkg.dir, "utility.R"))

# sort sh.o* files in new directories
move.logs(args$out.dir.work)

#
#	produce trees
#

# iqtree option
.make.iqtree.opt(args)

# Load alignments:
# and check how many did not run
infiles <- list.files(args$out.dir.output, 
                      pattern='InWindow_(.*)?_v2.fasta$',
                      full.names=TRUE, recursive=TRUE)

.check.alignments(infiles)

# If we haven't re-run alignments already, check if some are problematic:
cnd <- length(list.files(args$out.dir.work, pattern='readali2')) == 0
if(cnd){.check.or.resubmit.incompleted.alignments()}

# Extract info from name
.f <- function(reg, rep, x) as.integer(gsub(reg, rep, x))

infiles	<- data.table(FI=infiles)
infiles[, FO:= gsub('.fasta$','',FI)]
infiles[, PTY_RUN:= .f('^ptyr([0-9]+)_.*','\\1',basename(FI))]
infiles[, W_FROM := .f( '^.*_([0-9]+)_to.*$' , '\\1', basename(FI) )]
infiles[,PositionsExcised:=grepl('PositionsExcised',FI)]

# Delete the files without excision if the files with excision exist
setkey(infiles, PTY_RUN, W_FROM)
tmp <- infiles[,list(NUM=length(PositionsExcised)),by=c('PTY_RUN', 'W_FROM')]
infiles <- merge(infiles, tmp, by=c('PTY_RUN', 'W_FROM'))
tmp <- infiles[(PositionsExcised==F & NUM==2),]
file.remove(tmp$FI)

infiles <- infiles[!(PositionsExcised==F & NUM==2),]

# check tree completed: treefiles.
infiles[,FO_NAME:=paste0(FO,'.treefile')]
infiles[,FO_EXIST:=file.exists(FO_NAME)]
infiles <- infiles[FO_EXIST==FALSE]
if(args$verbose){
  cat('These trees were not built sucessfully...\n')
  print(head(infiles))
  cat('Build again...\n')
}

if( nrow(infiles) == 0 )
{
  
  # TODO: check that this is actually working...
  `%+%` <- function(x,y) paste0( x,y )
  cat('All trees were built, move on to next, analyse tree step...\n')
  dir <- dirname(args$controller) 
  
  cmd <- 'cd ' %+% dir %+% '\n'
  cmd <- cmd %+% "qsub -v STEP='atr' " %+% basename(args$controller) %+% '\n'

  cat(cmd)
  cat(system(cmd, intern=TRUE))
  
  stop( 'submitted next step of the analysis: atr)' )
}


# Set up jobs: put 2 iqtrees in a same subjob?
# Maybe randomly permute? before doing it?
djob <- infiles[, list(CMD=cmd.iqtree(FI, outfile=FO, pr=iqtree.pr, pr.args=iqtree.args)), by=c('PTY_RUN','W_FROM')]
djob[, ID:=ceiling(seq_len(.N)/2)]
djob<- djob[, list(CMD=paste(CMD, collapse='\n',sep='')), by='ID']
# cat(djob$CMD[1])

#	Create PBS job array
if(nrow(djob) > max.per.run){
  djob[, CASE_ID:= rep(1:max.per.run,times=ceiling(nrow(djob)/max.per.run))[1:nrow(djob)]]
  djob[, JOB_ID:= rep(1:ceiling(nrow(djob)/max.per.run),each=max.per.run)[1:nrow(djob)]]
  djob[,ID:=NULL]
}else{
  setnames(djob, 'ID', 'CASE_ID')
  djob[, JOB_ID:=1]
}


# now write header + cmd

pbshead	<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=NULL)

indexes <- djob[, unique(JOB_ID)]
for (i in seq_along(indexes)) {

  # Prepare Job command for submission
  tmp<- djob[JOB_ID==i,]

  hpc.array <- nrow(tmp)
  pbsJ_bool <- hpc.array > 1

  cmd <- tmp[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
  
  # Decide whether need -J (arrays) or not
  if(pbsJ_bool)
  {
        cmd<-cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]		
        pbshead1 <- paste0(pbshead, "\n#PBS -J 1-", hpc.array)
  }else{
        cmd <- cmd[, paste0(CASE, collapse='')]
        cmd <- gsub('1)', '', cmd)
        cmd <- gsub(';;', '', cmd)
        pbshead1 <- pbshead
  }
  
  tmp <- paste0(pbshead1,
                '\n module load anaconda3/personal \n source activate ',args$env_name)
  cmd<-paste(tmp,cmd,sep='\n')	
  # cat(substr(cmd, 1, 1000))
  
  #	Prepare sh file  
  date <- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
  date <- gsub(':','_',date)
  outfile <- paste("srx",paste0('job',i), date,'sh',sep='.')
  outfile <-file.path(args$out.dir.work, outfile)
  cat(cmd, file=outfile)

  # Change directory and submit job!
  cmd <- paste('cd', dirname(outfile), '\n')
  cmd <- paste(cmd, "qsub", outfile)
  cat(system(cmd, intern=TRUE))
}

