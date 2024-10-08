cat("\n\n===== make_trees.R =====\n\n")
# The phyloscanner run - make trees

# Preamble
# This script aims to make trees using each alignment.

# Load the required packages
library(data.table)
library(seqinr)
library(tidyverse)
library(dplyr)
library(optparse)

#
# Define input arguments that can be changed by users
#
option_list <- list(
  make_option(
    c("-v", "--verbose"),
    action = "store_true",
    default = TRUE,
    help = "Print extra output [default]",
    dest = "verbose"
  ),
  make_option(
    "--seed",
    type = "integer",
    default = 42L,
    help = "Random number seed [default %default]",
    dest = "seed"
  ),
  make_option(
    "--save_data",
    action = "store_true",
    default = TRUE,
    help = "Save data [default]",
    dest = 'if_save_data'
  ),
  make_option(
    "--env_name",
    type = "character",
    default = 'phyloenv',
    help = "Conda environment name [default]",
    dest = 'env_name'
  ),
  make_option(
    "--iqtree_method",
    type = "character",
    default = "GTR+F+R6",
    help = "Methods in IQTREE [default]",
    dest = 'iqtree_method'
  ),
  make_option(
    "--iqtree_root",
    type = "character",
    default = NULL,
    help = "Root in IQTREE [default]",
    dest = 'iqtree_root'
  ),
  make_option(
    "--pkg_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to package directory, used as long we don t build an R package [default]",
    dest = 'pkg.dir'
  ),
  make_option(
    "--out_dir_base",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where all output is stored [default]",
    dest = 'out.dir'
  ),
  make_option(
    "--controller",
    type = "character",
    default = NA_character_, 
    help = "Path to sh script directing the full analysis",
    dest = 'controller'
  ),
  make_option(
    "--date",
    type = 'character',
    default = as.character(Sys.Date()),
    metavar = '"YYYY-MM-DD"',
    help = 'As of date to extract data from.  Defaults to today.',
    dest = 'date'
  ),
  make_option(
    "--walltime_idx",
    type = "integer",
    default = 2,
    help = "Indicator for amount of resources required by job. Values ranging from 1 (lala) to 3 (lala)",
    dest = "walltime_idx"
  ),
  optparse::make_option(
    "--csv-runners",
    type = "character",
    default = NA_character_,
    help = "Path to csv file containing the specifications for each person running the job. If NA, the jobs will be run by the user running this script. (defaults to NA)",
    dest = "runners"
  ),
  optparse::make_option(
    "--dryrun",
    type = "logical",
    default = FALSE,
    help = "Runs only for the first 2 ptyr",
    dest = "dryrun"
  )
)

args <-  parse_args(OptionParser(option_list = option_list))

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

.write.job <- function(DT)
{
    hpc_array <- max(DT$CASE_ID)
    if(hpc_array <= 1)  hpc_array <- NA
  # Define PBS header for job scheduler
  pbshead <- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, 
                              hpc.nproc=hpc.nproc, hpc.mem=hpc.mem, 
                              hpc.walltime=hpc.walltime,
                              hpc.q=hpc.q,  
                              hpc.array = hpc_array, 
                              hpc.load=paste0('\n module load anaconda3/personal \n source activate ',args$env_name)
                              )
  
  if(is.na(hpc_array)){
    cmd <- DT$CMD
  }else{
    job[, CASE_ID := 1:.N ]
    cmd <- job[, list(CASE = paste0(CASE_ID, ')\n', CMD, ';;\n')), by = 'CASE_ID']
    cmd <- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',
                             paste0(CASE, collapse = ''),
                             'esac')]
  }
  cmd <- paste(pbshead, cmd, sep = '\n')
  cmd
}

#
# test
#
if(0){
  args <- list(
    out.dir = "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS_UVRI_TSI2/",
    pkg.dir="/rds/general/user/ab1820/home/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software",
    iqtree_method="GTR+F+R6",
    env_name = 'phylostan',
    date = '2022-07-26',
    seed = 42,
    walltime_idx=1,
    controller="/rds/general/user/ab1820/home/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/runall_TSI_pairs.sh",
    dryrun=TRUE
  )

  if (Sys.info()['user'] == "andrea") {
    args$out.dir = "/home/andrea/HPC/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521"
    args$pkg.dir = "~/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1521"
    args$prog.dir = "~/git/phyloscanner"
    args$controller= "~/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1521/runall_dryrun.sh"
  }
}

# Source functions
source(file.path(args$pkg.dir, "utility.R"))

# 
# Chose PBS specifications according to PBS index
# Idea is that this script can be run multiple times with increasing res reqs
list(
  `1`= list(hpc.select = 1, hpc.nproc = 1, hpc.walltime = 4 , hpc.mem = "2gb" ,hpc.q = NA, max.per.run=950),
  `2`= list(hpc.select = 1, hpc.nproc = 1, hpc.walltime = 7, hpc.mem = "40gb" ,hpc.q = NA, max.per.run=950),
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

# sort sh.o* files in new directories
move.logs(args$out.dir.work)

#
#	produce trees
#

# iqtree option
.make.iqtree.opt(args)

# Look for RDS file containing subjobs CMDS, if can't find, KEEP GOING
cmds.path <- file.path(args$out.dir.work, 'trees_commands.rds')

if(file.exists(cmds.path))
{
  # Load commands
  move.logs(args$out.dir.work)
  infiles <- readRDS(cmds.path)

}else{

  # Load alignments:
  # and check how many did not run
  infiles <- list.files(args$out.dir.output, 
                        pattern='InWindow_(.*)?_v2.fasta$',
                        full.names=TRUE, recursive=TRUE)

  .check.alignments(infiles)

  # If we haven't re-run alignments already, check if some are problematic:
  # This step should be taken care of automatically
  # if(cnd){.check.or.resubmit.incompleted.alignments()}
  
  # Extract info from name
  .f <- function(reg, rep, x) as.integer(gsub(reg, rep, x))
  
  infiles	<- data.table(FI=infiles)
  infiles[, FO:= gsub('.fasta$','',FI)]
  infiles[, PTY_RUN:= .f('^ptyr([0-9]+)_.*','\\1',basename(FI))]
  infiles[, W_FROM := .f( '^.*_([0-9]+)_to.*$' , '\\1', basename(FI) )]
  infiles[,EXCISED:=grepl('PositionsExcised',FI)]
  
  # Delete the files without excision if the files with excision exist
  setkey(infiles, PTY_RUN, W_FROM)
  tmp <- infiles[,list(NUM=length(EXCISED)),by=c('PTY_RUN', 'W_FROM')]
  infiles <- merge(infiles, tmp, by=c('PTY_RUN', 'W_FROM'))
  tmp <- infiles[(EXCISED==F & NUM==2),]
  file.remove(tmp$FI)
  
  infiles <- infiles[!(EXCISED==F & NUM==2),]
  infiles[,  `:=` (NUM=NULL, EXCISED=NULL) ]
  
  if (args$dryrun){
     cat("Dry run: only running for first 2 ptyr\n")
     infiles <- subset(infiles, PTY_RUN %in% c(1,2))
  }
  
  # check tree completed: treefiles.
  infiles[,FO_NAME:=paste0(FO,'.treefile')]
  
  saveRDS(infiles, cmds.path)
}

infiles[,FO_EXIST:=file.exists(FO_NAME)]
infiles[, mean(FO_EXIST)]
infiles[, sum(!FO_EXIST)]
infiles <- infiles[FO_EXIST==FALSE]

if( nrow(infiles) == 0 )
{
  # ISN'T THIS BEAUTIFUL?
  qsub.next.step(file=args$controller,
                 next_step='atr', 
                 res=1,
                 redo=0
  )
  stop('Building Trees step completed, submitted the following task')
}


# Write command for each alignment
djob <- infiles[, list(CMD=cmd.iqtree(infile.fasta=FI, outfile=FO, pr=iqtree.pr, pr.args=iqtree.args)), by=c('PTY_RUN','W_FROM')]

# put 2 cmds per array job (not sure why XX did this but maybe some are really short)
# at least permute though (so that allocation of long ones is random)
djob <- djob[sample(1:.N),  ]
djob[, ID:=ceiling(seq_len(.N)/2)]
djob<- djob[, list(CMD=paste(CMD, collapse='\n',sep='')), by='ID']
djob[,ID:=NULL]

# group into jobs
n_jobs <- ceiling( nrow(djob) / max.per.run)
idx <- 1:nrow(djob)
djob[, CASE_ID := rep(1:max.per.run, times = n_jobs)[idx] ]
djob[, JOB_ID := rep(1:n_jobs, each = max.per.run)[idx] ]

djob2 <- djob[, .(CMD=.write.job(.SD)), by=JOB_ID]

# Load runners specifications
split_jobs_by_n <- 1
if(!is.na(args$runners))
{
  drunners <- fread(args$runners)
  split_jobs_by_n <- nrow(drunners)
}

if ( split_jobs_by_n == 1 | nrow(djob) <= 1000 )
{
    ids <- submit_jobs_from_djob(djob2, output_type = "id", prefix = "srx")

    # qsub alignment step again, to check whether everything has run...
    qsub.next.step(file=args$controller,
        ids=ids, 
        next_step='btr', 
        res=args$walltime_idx + 1, 
        redo=1
    )
}else{

  # split djob in split_jobs_by_n data.tables
  # then perform the above individually
  # by first taking care of different path specifications among runners

  person_to_run <- rep(1:split_jobs_by_n, length.out = nrow(djob2))
  submit_user_script <- NA_character_

  for (person in 1:split_jobs_by_n){

    # Subset to jobs for specific person
    djob_person <- djob2[person_to_run == person]
    # Adapt the specifications according to the person/runner
    djob_person <- adapt_jobspecs_to_runner(
            djob_person,
            drunners,
            idx = person
        )
    # write the pbs files
    pbs_file_person <- submit_jobs_from_djob(
            djob_person,
            output_type = "outfile",
            prefix = "srx"
        )
    # Append the pbs files to the script that each user can submit to queue them
    submit_user_script <- append_pbs_file_person(
            script = submit_user_script,
            pbs = pbs_file_person,
            usr = drunners[index == person,  user_name]
        )
  }

  # Write the script that each user can submit
  outfile = file.path(args$out.dir.work, 'submit_user_srx.sh')
  write( submit_user_script, file = outfile)

}

cat("End of script\n")

