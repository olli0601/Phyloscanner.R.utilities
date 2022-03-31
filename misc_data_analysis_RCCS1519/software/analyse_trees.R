cat("\n\n===== analyse_trees.R =====\n\n")
# The phyloscanner run for TSI estimates - analyse trees

# Preamble
# This script aims to analyse trees to identify ancestral relationships.

# Load the required packages
library(data.table)
library(seqinr)
library(tidyverse)
library(dplyr)
require(phyloscannerR)

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
    default = 'phyloenv',
    help = "Conda environment name [default]",
    dest = 'env_name'
  ),
  optparse::make_option( 
    "--blacklistReport", 
    action="store_true", 
    help="If present, output a CSV file of blacklisted tips from each tree.",
    dest = 'blacklist.report'
    ),
  optparse::make_option(
    "--distanceThreshold",
    action="store", 
    default = NULL,
    help = "Maximum distance threshold on a window for a relationship to be reconstructed between two hosts on that window.[default]",
    dest = 'distance.threshold'
  ),
  # optparse::make_option(
  #   "--fileNameRegex", 
  #   action="store", 
  #   default="^(?:.*\\D)?([0-9]+)_to_([0-9]+).*$",
  #   help="Regular expression identifying window coordinates in tree file names. Two capture groups: start and end; if the latter is missing then the first group is a single numerical identifier for the window. If absent, input will be assumed to be from the phyloscanner pipeline.",
  #   dest = 'file.name.regex'
  #   ),
  optparse::make_option(
    "--maxReadsPerHost",
    action="store", 
    default = NULL,
    help = "If given, blacklist to downsample read counts (or tip counts if no read counts are identified) from each host to this number. [default]",
    dest = 'max.reads.per.host'
  ),
  optparse::make_option(
    "--minReadsPerHost",
    type = "integer",
    default = 1,
    help = "Blacklist hosts from trees where they have less than this number of reads. [default]",
    dest = 'min.reads.per.host'
  ),
  optparse::make_option(
    "--multinomial", 
    action="store_true", 
    help="Use the adjustment for missing and overlapping windows as described in Ratmann et al., Nature Communications, 2019.",
    dest='multinomial'
  ),
  optparse::make_option( 
    "--normRefFileName",
    action="store", 
    default = NULL,
    help="Name of a file giving a normalisation constant for every genome position. Cannot be used simultaneously with -nc. ",
    dest = 'norm.ref.file.name'
  ),
  optparse::make_option(
    "--normStandardiseGagPol", 
    action="store_true", 
    help="An HIV-specific option: if true, the normalising constants are standardised so that the average on gag+pol equals 1. Otherwise they are standardised so the average on the whole genome equals 1.",
    dest = 'norm.standardise.gag.pol'
  ),
  optparse::make_option(
    "--noProgressBars", 
    action="store_true", 
    help="If --verbose, do not display progress bars",
    dest = 'no.progress.bars'
  ),
  optparse::make_option( 
    "--outgroupName",
    action="store",
    default =NULL,
    help="The name of the tip in the phylogeny/phylogenies to be used as outgroup (if unspecified, trees will be assumed to be already rooted).",
    dest = 'outgroup.name'
  ),
  optparse::make_option( 
    "--outputNexusTree",
    action="store_true",
    help="Standard output of annotated trees are in PDF format. If this option is present, output them as NEXUS instead.",
    dest = 'output.nexus.tree'
  ),
  optparse::make_option(
    "--postHocCountBlacklisting", 
    action="store_true", 
    help="Perform minimum read and tip based blacklisting as a separate step at the end of the analysis. (A legacy option).",
    dest = 'post.hoc.count.blacklisting'
  ),
  optparse::make_option( 
    "--ratioBlacklistThreshold", 
    action="store", 
    default=0, 
    help="Used to specify a read count ratio (between 0 and 1) to be used as a threshold for blacklisting. --parsimonyBlacklistK and/or --duplicateBlacklist must also be used. If --parsimonyBlacklistK is used, a subgraph will be blacklisted if the ratio of its read count to the total read count from the same host (in that tree) is strictly less than this threshold.",
    dest = 'ratio.blacklist.threshold'
  ),
  optparse::make_option(
    "--relaxedAncestry", 
    action="store_true", 
    help="If absent, directionality can be inferred so long as at least one subraph from one host is descended from one from the other, and no pair of subgraphs exist in the opposite direction. Otherwise it is required that every subgraph from one host is descended from one from the other.",
    dest = "relaxed.ancestry"
  ),
  optparse::make_option(
    "--skipSummaryGraph", 
    action="store_true", 
    help="If present, do not output a simplified relationship graph",
    dest = 'skip.summary.graph'
  ),
  optparse::make_option(
    "--zeroLengthAdjustment",
    action="store_true",
    help="If present when allowMultiTrans is switched on, two hosts are classified as complex if their MRCAs are in the same multifurcation.",
    dest = 'zero.length.adjustment'
    ),
  optparse::make_option(
    "--pkg_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to package directory, used as long we don t build an R package [default]",
    dest = 'prj.dir'
  ),
  optparse::make_option(
    "--prog_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to the phyloscanner [default]",
    dest = 'prog.dir'
  ),
  optparse::make_option(
    "--out_dir_base",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where all output is stored [default]",
    dest = 'out.dir'
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
# if(1){
#   args$ort=T
#   args$date = '19037'
#   args$env_name = 'phylostan'
#   args$norm.ref.file.name = "~/phyloscanner/InfoAndInputs/HIV_DistanceNormalisationOverGenome.csv"
#   args$outgroup.name = "B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
#   args$output.nexus.tree = T
#   args$ratio.blacklist.threshold = 0.005
#   args$skip.summary.graph = T
#   args$out.dir = NA
#   args$prj.dir = NA
#   args$prog.dir = NA
# }


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
      '/rds/general/user/xx4515/home/phyloscanner/phyloscanner_analyse_trees.R'
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
if (is.na(args$prog.dir))
{
  args$prog.dir <- here::here()
}
if(is.null(args$distance.threshold)){
  args$distance.threshold <- -1
}
#
# Add constants that should not be changed by the user
#
max.per.run <- 4900

# Set default output directories relative to out.dir
args$out.dir.work <-
  file.path(args$out.dir, paste0(args$date, "_phsc_work"))
args$out.dir.output <-
  file.path(args$out.dir, paste0(args$date, "_phsc_output"))

if(args$verbose){
  print(args)
}

tmp <- setdiff(names(args),c("verbose","if_save_data","env_name","prj.dir","prog.dir",
                             "out.dir","date","help",'out.dir.data','out.dir.work',
                             'out.dir.output','norm.ref.file.name'))
tmpv <- args[tmp]
out.dir.analyse.trees <- file.path(args$out.dir, 
                                   paste0(args$date,'_phsc_phscrelationships_',
                                   paste(gsub('\\.','_',names(tmpv)),
                                         gsub('\\.|-','',tmpv),collapse='_',sep = '_')))
out.dir.analyse.trees <- gsub(' ','_',out.dir.analyse.trees)
# simplify name
out.dir.analyse.trees <- gsub('_TRUE','_T',out.dir.analyse.trees)
out.dir.analyse.trees <- gsub('_seed','_sd',out.dir.analyse.trees)
out.dir.analyse.trees <- gsub('_distance_threshold','_sdt',out.dir.analyse.trees)
out.dir.analyse.trees <- gsub('_max_reads_per_host','_dsl',out.dir.analyse.trees)
out.dir.analyse.trees <- gsub('_min_reads_per_host','_mr',out.dir.analyse.trees)
out.dir.analyse.trees <- gsub('_multinomial','_mlt',out.dir.analyse.trees)
out.dir.analyse.trees <- gsub('_no_progress_bars','_npb',out.dir.analyse.trees)
out.dir.analyse.trees <- gsub('_outgroup_name','_og',out.dir.analyse.trees)
out.dir.analyse.trees <- gsub('_post_hoc_count_blacklisting','_phcb',out.dir.analyse.trees)
out.dir.analyse.trees <- gsub('_ratio_blacklist_threshold','_rtt',out.dir.analyse.trees)
out.dir.analyse.trees <- gsub('_relaxed_ancestry','_rla',out.dir.analyse.trees)
out.dir.analyse.trees <- gsub('_zero_length_adjustment','_zla',out.dir.analyse.trees)
print(out.dir.analyse.trees)
dir.create(out.dir.analyse.trees)

# Source functions
source(file.path(args$prj.dir, "utility.R"))

#	Set phyloscanner variables	
control	<- list()
control$allow.mt <- TRUE
control$alignment.file.directory = NULL 
control$alignment.file.regex = NULL
control$blacklist.report = args$blacklist.report
control$blacklist.underrepresented = FALSE	
control$count.reads.in.parsimony = TRUE
control$distance.threshold <- args$distance.threshold
control$do.dual.blacklisting = FALSE					
control$duplicate.file.directory = NULL
control$duplicate.file.regex = NULL
control$file.name.regex = "^(?:.*\\D)?([0-9]+)_to_([0-9]+).*$"
control$guess.multifurcation.threshold = FALSE
control$max.reads.per.host <- args$max.reads.per.host
control$min.reads.per.host <- args$min.reads.per.host
control$min.tips.per.host <- 1	
control$multifurcation.threshold = 1e-5
control$multinomial= args$multinomial
control$norm.constants = NULL
control$norm.ref.file.name = args$norm.ref.file.name
control$norm.standardise.gag.pol = args$norm.standardise.gag.pol
control$no.progress.bars = args$no.progress.bars
control$overwrite = TRUE
control$outgroup.name = args$outgroup.name
control$output.dir = out.dir.analyse.trees
control$output.nexus.tree = args$output.nexus.tree
control$outputRDA = TRUE
control$parsimony.blacklist.k = 15
control$prune.blacklist = FALSE
control$post.hoc.count.blacklisting <- args$post.hoc.count.blacklisting
control$ratio.blacklist.threshold = args$ratio.blacklist.threshold
control$raw.blacklist.threshold = 3			
control$recombination.file.directory = NULL
control$recombination.file.regex = NULL
control$relaxed.ancestry = args$relaxed.ancestry
control$sankoff.k = 15
control$sankoff.unassigned.switch.threshold = 0
control$seed = args$seed
control$skip.summary.graph = args$skip.summary.graph
control$splits.rule = 's'
control$tip.regex = '^(.*)-fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'
control$tree.file.regex = "^(.*)\\.treefile$" 
control$tree.file.extension = '.treefile'
control$use.ff = FALSE
control$user.blacklist.directory = NULL 
control$user.blacklist.file.regex = NULL
control$verbosity = 1	
control$zero.length.adjustment = args$zero.length.adjustment

#	Make scripts
df <- tibble(F=list.files(args$out.dir.output,pattern = 'ptyr*'))
df <- df %>%
  mutate(TYPE:= gsub('ptyr([0-9]+)_(.*)','\\2', F),
         RUN:= as.integer(gsub('ptyr([0-9]+)_(.*)','\\1', F))) %>%
  mutate(TYPE:= gsub('^[^\\.]+\\.([a-z]+)$','\\1',TYPE)) %>%
  spread(TYPE, F) %>%
  set_names(~ str_to_upper(.))

valid.input.args <- cmd.phyloscanner.analyse.trees.valid.args(args$prog.dir)
cmds <- vector('list',nrow(df))

for(i in seq_len(nrow(df)))
{
  #	Set input args
  control$output.string <- paste0('ptyr',df$RUN[i])
  
  #	Make script
  tree.input <- file.path(args$out.dir.output, df$TREES[i])
  cmd <- cmd.phyloscanner.analyse.trees(args$prog.dir, 
                                        tree.input, 
                                        control,
                                        valid.input.args=valid.input.args)
  cmds[[i]] <- cmd		
}	
if(args$verbose){
  cat(cmds[[1]])
}


#
# Submit array job to HPC
#
#	Make header
hpc.load			<- paste0("module load anaconda3/personal \n source activate ", args$env_name)
hpc.select			<- 1
hpc.nproc			<- 1
hpc.walltime		<- 23
if(1)
{
  hpc.q			<- NA
  hpc.mem			<- "36gb"
}
hpc.array			<- length(cmds)
pbshead		<- "#!/bin/sh"
tmp			<- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
pbshead		<- paste(pbshead, tmp, sep = "\n")
tmp			<- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
pbshead 	<- paste(pbshead, tmp, sep = "\n")
pbshead 	<- paste(pbshead, "#PBS -j oe", sep = "\n")
if(!is.na(hpc.array))
  pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
if(!is.na(hpc.q))
  pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
pbshead 	<- paste(pbshead, hpc.load, sep = "\n")
cat(pbshead)
#	Make array job
for(i in 1:length(cmds))
  cmds[[i]]<- paste0(i,')\n',cmds[[i]],';;\n')
cmd		<- paste0('case $PBS_ARRAY_INDEX in\n',paste0(cmds, collapse=''),'esac')
cmd		<- paste(pbshead,cmd,sep='\n')

#	Submit job
outfile		<- gsub(':','',paste("phsc",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
outfile		<- file.path(args$out.dir.work, outfile)
cat(cmd, file=outfile)
cmd 		<- paste("qsub", outfile)
cat(cmd)
cat(system(cmd, intern= TRUE))


