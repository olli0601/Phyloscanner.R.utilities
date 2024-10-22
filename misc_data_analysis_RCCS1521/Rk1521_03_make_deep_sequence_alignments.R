# The phyloscanner run - make alignments
# Preamble
# The set of scripts aims to run phyloscanner.
# This script aims to make alignments in each potential transmission network.
# TODO: check whether problematic_windows is fixed...
# Also qsub.next...

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
    "--windows_start",
    type = "integer",
    default = 800,
    help = "Left extremity of leftmost window[default]",
    dest = "windows_start"
  ),
  optparse::make_option(
    "--windows_end",
    type = "integer",
    default = 9175,
    help = "Left extremity of rightmost windows[default]",
    dest = "windows_end"
  ),
  optparse::make_option(
    "--sliding_width",
    type = "integer",
    default = NA_integer_,
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
    "--prog_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where phyloscanner is stored [default]",
    dest = 'prog.dir'
  ),
  optparse::make_option(
    "--reference",
    type = "character",
    default = NA_character_,
    help = "path to reference consensus sequences .fasta file required for the analysis. Basename is sufficient if file is 'standard'.",
    dest = 'reference'
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
    "--tsi_analysis",
    type = 'logical',
    default = FALSE,
    help = 'Indicator on whether we want to perform a Time Since Infection analysis[default]. Close to being deprecated',
    dest = 'tsi_analysis'
  ),
  optparse::make_option(
    "--mafft",
    type = 'character',
    default = '--globalpair --maxiterate 1000',
    help = 'options for alignment program mafft',
    dest = 'mafft.opt'
  ),
  optparse::make_option(
    "--rm_vloops",
    action = "store_true",
    default = TRUE,
    help = "Indicator on whether to avoid alignments on vloop region.[default]",
    dest = 'rm_vloops'
  ),
  optparse::make_option(
    "--controller",
    type = "character",
    default = NA_character_,
    help = "Path to sh script irecting the full anysis",
    dest = 'controller'
  ),
  optparse::make_option(
    "--walltime_idx",
    type = "integer",
    default = 2,
    help = "Indicator for amount of resources required by job. Values ranging from 1 (lala) to 3 (lala)",
    dest = "walltime_idx"
  ),
  optparse::make_option(
    "--pqeelab",
    type = "logical",
    default = FALSE,
    help = "Indicator for whether to submit (some) jobs to pqeelab queue. [default FALSE]",
    dest = "pqeelab"
  ),
  optparse::make_option(
    "--phsc-runs",
    type = "character",
    default = NA_character_,
    help = "Path to RDS file containing phyloscanner runs",
    dest = 'infile.runs'
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
args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))


# stop if arguments are not provided:
stopifnot("Path to infile.runs is misspecified" = file.exists(args$infile.runs))
if( is.na(args$sliding_width) ) stop('No sliding_width provided')

#
# Helpers
#

.percent <- function(x){paste0(format(round(x*100, 2), nsmall=2), '%') }

.check.existing.outputs <- function(regex, outdir, pattern='')
{
  # returns a vector of TRUE or FALSE based on whether regex is found in outdir
  if( ! file.exists(outdir) ){return(rep(NA_character_, length(regex)))}

  files <- list.files(outdir, pattern=pattern)
  # .f <- function(rgx) any(grepl(x=files,rgx))
  .f <- function(rgx) grep(x=files,rgx, value=T)[1]
  sapply(regex, .f) -> tmp
  unlist(tmp)
}

.remove_problematic_windows <- function(DT)
{
  # check if file problematic_windows.csv exists
  tmp <- list.files(args$out.dir.output, pattern='problematic_windows.csv',
                    recursive = T, full.names = T)
  if( length(tmp) == 0 )
  {
    cat('No problematic windows files found\n')
    return(DT)
  }
  names(tmp) <- basename(dirname(tmp))
  dprob <- lapply(tmp, fread)
  dprob <- rbindlist(dprob, idcol = 'PTY_RUN')
  names(dprob) <- c('PTY_RUN', 'FROM', 'TO', 'PBS_ARRAY_INDEX', 'SH', 'QUEUE')
  dprob[,  PTY_RUN:=as.integer(gsub('[A-z]|_', '', PTY_RUN))]

  # dprob <- dprob[!SH %in% c('readali.job1.Sat_Jul_23_205803_2022.sh', 'readali.job1.Sun_Jul_24_133928_2022.sh')]

  # This below is more for debugging purposes
  # names of logfiles
  .f <- function(x)
  {
    # finds prefix for LOG, needs only PBS_ARRAY_INDEX
    dir <- file.path(args$out.dir.work, gsub('.sh$','',x))
    fil <- list.files(dir,  full.names = T)
    pre <- unique(gsub('\\.[0-9]+$','.',fil))
    if(length(pre) == 1) return(pre)
    NA_character_
  }
  dprob[, LOG:=paste0(.f(SH), PBS_ARRAY_INDEX), by='SH']

  # compare with outputs
  merge(
    DT[, .(W_FROM, PTY_RUN, OUT1, OUT2)],
    dprob[, .(W_FROM = FROM, PTY_RUN, LOG)],
    all.x=TRUE,
    by=c('W_FROM', 'PTY_RUN')
  ) -> dlogs

  tmp <- dlogs[ is.na(LOG), .(W_FROM, PTY_RUN, OUT2) ]
  # stopifnot(tmp[, all(is.na(OUT2))])
  tmp$OUT2 <- NULL

  merge(DT, tmp, by=c('W_FROM', 'PTY_RUN'))
}

.modify_cmd_for_existing_v1 <- function(cmd, outdir, out1)
{
  if(length(out1) == 0){ stop('error: empty input to .modify_cmd_for_existing_v1')}
  # cmd <- pty.c[, CMD[1]]
  out1 <- basename(out1)
  out2 <- gsub('\\.fasta', '\\_v2.fasta', out1)
  lines <- strsplit(cmd, '\n')[[1]]

  # remove useless parts
  rm1_start <- grep('phyloscanner_make_trees', lines)
  rm1_end <- grep('Performing realignment', lines) - 1
  rm1 <- rm1_start:rm1_end
  rm2 <- grep("(?=.*mv.*)(?=.*AlignedReads.*)", lines, perl = TRUE)
  rm2 <- (rm2-1):(rm2+1)
  rm3 <- grep('problematic_windows.csv', lines)
  lines <- lines[ -c(rm1, rm2, rm3)]

  # copy out1 to work dir
  cp_pos <- grep('^cp', lines)[1]
  cp <- lines[cp_pos]
  cp <- gsub('^cp "(.*?)" (.*?)$', paste0('cp "', file.path(outdir, out1)  ,'" \\2'), cp)
  lines[cp_pos] <- cp

  # substitute AlignedReads* with the name of our file
  mafft_pos <- grep('mafft', lines)[1]
  mafft <- lines[mafft_pos]
  mafft2 <- gsub('\\\t', '', mafft)
  mafft2 <- gsub('\\$file', out1, mafft2)
  mafft2 <- gsub('\\$\\{file//.fasta/_v2.fasta\\}', out2, mafft2)
  lines[mafft_pos] <- mafft2
  rm4 <- c(mafft_pos + 1, mafft_pos -1)
  lines <- lines[ -rm4]

  cmd <- paste0(lines, collapse='\n')
  return(cmd)
}

.write.job <- function(DT,path)
{
  # Define PBS header for job scheduler
  pbshead <- cmd.hpcwrapper.cx1.ic.ac.uk(
    hpc.select = hpc.select,
    hpc.nproc = hpc.nproc,
    hpc.walltime = hpc.walltime,
    hpc.q = hpc.q,
    hpc.mem = hpc.mem,
    hpc.array = DT[, max(CASE_ID)],
    hpc.load = "module load intel-suite/2015.1 mpi raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"
  )

  #cmd <- DT[, list(CASE = paste0(CASE_ID, ')\n', CMD, ';;\n')), by = 'CASE_ID']
  cmd <- paste(
    #paste0('\nLIST_NAME=$EPHEMERAL/input_list.txt'),
    paste0('\nLIST_NAME=',args$out.dir.work,'/input_list.txt'),
    'INPUTS=$(head -n $PBS_ARRAY_INDEX $LIST_NAME | tail -1)',
    'PREFIX=$(awk \'{print $1}\' <<< "$INPUTS")',
    'WINDOW=$(awk \'{print $2}\' <<< "$INPUTS")',
    'WINDOW_START=$(echo $WINDOWS | awk -F\',\' \'{print $1}\')',
    'WINDOW_END=$(echo $WINDOWS | awk -F\',\' \'{print $2}\')',
    '\n',
    'echo "PBS Job Id PBS_JOBID is ${PBS_JOBID}"',
    'echo "PBS job array index PBS_ARRAY_INDEX value is ${PBS_ARRAY_INDEX}"',
    'echo "PTYR is ${PREFIX}"',
    'echo "WINDOW ${WINDOW}"',
    '\n',
    paste0('EXCISION_COORDS=$(cat ',path,'/excision_coords.txt)'),
    sep = "\n"
  )
  cmd <- paste0(cmd,'\n')
  cmd <- paste0(cmd, DT$CMD[1]) # just print the first generic command
  #cmd <-    cmd[, paste0('case $PBS_ARRAY_INDEX in\n',
  #                       paste0(CASE, collapse = ''),
  #                       'esac')]
  cmd <- paste(pbshead, cmd, sep = '\n')
  cmd
}

#
# test
#
if(0)
{
  args <- list(
    verbose = T,
    seed = 42,
    windows_start=550L,
    windows_end=9500L,
    sliding_width = 10L,
    window_size = 250L,
    window_cutoff = NA,
    n_control = 0,
    cluster_size = 50,
    if_save_data = T,
    date = '2024-09-23',
    out.dir = "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521",
    pkg.dir = "/rds/general/user/ab1820/home/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1521",
    prog.dir = "/rds/general/user/ab1820/home/git/phyloscanner",
    reference = 'ConsensusGenomes.fasta',
    tsi_analysis=FALSE,
    rm_vloops=TRUE,
    mafft.opt='--globalpair --maxiterate 1000',
    controller='/rds/general/user/ab1820/home/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/runall_TSI_pairs2.sh',
    walltime_idx = 1,
    infile.runs = "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/240923_RCCSUVRI_phscinput_runs.rds",
    dryrun = TRUE
  )

  if (Sys.info()['user'] == "andrea") {
    args$out.dir = "/home/andrea/HPC/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521"
    args$pkg.dir = "~/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1521"
    args$prog.dir = "~/git/phyloscanner"
    args$controller= "~/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1521/runall_dryrun.sh"
    args$infile.runs = "/home/andrea/HPC/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1521/240923_RCCSUVRI_phscinput_runs.rds"
  }
}

# Source functions
source(file.path(args$pkg.dir, "utility.R"))

# I can run at most 1000 simultaneous jobs on the short q.
# however, 10'000 can be put in the same job.
list(
  `1`= list(hpc.select = 1, hpc.nproc = 1, hpc.walltime = 7 , hpc.mem = "4gb",hpc.q = NA, max.per.run=10000),
  `2`= list(hpc.select = 1, hpc.nproc = 1, hpc.walltime = 24, hpc.mem = "64gb",hpc.q = NA, max.per.run=10000),
  `3`= list(hpc.select = 1, hpc.nproc = 1, hpc.walltime = 71, hpc.mem = "64gb",hpc.q = NA, max.per.run=10000)
) -> pbs_headers
tmp <- pbs_headers[[args$walltime_idx]]
cat('selected the following PBS specifications:\n')
print(tmp)
invisible(list2env(tmp,globalenv()))

# Add constants that should not be changed by the user
#
dir.data <-  '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
dir.analyses <- '/rds/general/project/ratmann_deepseq_analyses/live'
dir.net <- file.path(args$out.dir, "potential_nets")

if (Sys.info()['user'] == "andrea") {
  dir.data <-  "/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live"
  dir.analyses <- "/home/andrea/HPC/project/ratmann_deepseq_analyses/live"
}

ifelse(
  !is.na(args$window_cutoff),
  paste0('_n_control_', args$n_control), ''
) -> tmp


# Set default output directories relative to out.dir
args$date <- gsub('-','_',args$date)

.f <- function(x)
{
  dir <- file.path(args$out.dir, paste0(args$date, x))

  # If date has been passed as par, (ie. args$date != today's date)
  # check that everything exists...
  cnd <-  as.Date(args$date, format='%Y_%m_%d') != as.character(Sys.Date())

  if(!dir.exists(dir))
  {
    if ( cnd )
      stop('No existing directories prefixed by the input --date')
    dir.create(dir)
  }
  dir
}
args$out.dir.data <- .f('_phsc_input')
args$out.dir.work <- .f('_phsc_work')
args$out.dir.output <- .f('_phsc_output')


# Look for RDS file containing subjobs CMDS, if can't find, KEEP GOING
cmds.path <- file.path(args$out.dir.work, 'align_commands.rds')
#cmds.path <- file.path(args$out.dir.work, 'input_list.txt')

# write excision coords file
default.coord = TRUE
if(!is.na(default.coord)){
  if(default.coord){
    coords <- paste0('823,824,825,892,893,894,907,908,909,1012,1013,1014,1156,1157,1158,1384,1385,1386,1444,1445,1446,1930,1931,1932,1957,1958,1959,2014,2015,2016,2023,2024,2025,2080,2081,2082,2134,2135,2136,2191,2192,2193,2280,2281,2282,2283,2284,2285,2298,2299,2300,2310,2311,2312,2316,2317,2318,2319,2320,2321,2322,2323,2324,2340,2341,2342,2346,2347,2348,2349,2350,2351,2352,2353,2354,2355,2356,2357,2358,2359,2360,2373,2374,2375,2379,2380,2381,2385,2386,2387,2388,2389,2390,2391,2392,2393,2394,2395,2396,2400,2401,2402,2409,2410,2411,2412,2413,2414,2415,2416,2417,2424,2425,2426,2430,2431,2432,2436,2437,2438,2439,2440,2441,2442,2443,2444,2457,2458,2459,2460,2461,2462,2463,2464,2465,2469,2470,2471,2472,2473,2474,2478,2479,2480,2481,2482,2483,2496,2497,2498,2499,2500,2501,2502,2503,2504,2505,2506,2507,2514,2515,2516,2517,2518,2519,2520,2521,2522,2526,2527,2528,2529,2530,2531,2535,2536,2537,2670,2671,2672,2679,2680,2681,2703,2704,2705,2709,2710,2711,2733,2734,2735,2742,2743,2744,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2769,2770,2771,2772,2773,2774,2778,2779,2780,2811,2812,2813,2814,2815,2816,2817,2818,2819,2823,2824,2825,2841,2842,2843,2847,2848,2849,2850,2851,2852,2856,2857,2858,2865,2866,2867,2871,2872,2873,2892,2893,2894,2895,2896,2897,2901,2902,2903,2904,2905,2906,2952,2953,2954,2961,2962,2963,3000,3001,3002,3015,3016,3017,3018,3019,3020,3030,3031,3032,3042,3043,3044,3084,3085,3086,3090,3091,3092,3099,3100,3101,3111,3112,3113,3117,3118,3119,3135,3136,3137,3171,3172,3173,3177,3178,3179,3180,3181,3182,3189,3190,3191,3192,3193,3194,3204,3205,3206,3210,3211,3212,3222,3223,3224,3228,3229,3230,3237,3238,3239,3246,3247,3248,3249,3250,3251,3255,3256,3257,3261,3262,3263,3396,3397,3398,3501,3502,3503,3546,3547,3548,3705,3706,3707,4425,4426,4427,4449,4450,4451,4503,4504,4505,4518,4519,4520,4590,4591,4592,4641,4642,4643,4647,4648,4649,4656,4657,4658,4668,4669,4670,4671,4672,4673,4692,4693,4694,4722,4723,4724,4782,4783,4784,4974,4975,4976,5016,5017,5018,5067,5068,5069,', paste(seq(6615,6811,by=1),collapse = ','),  ",", paste(seq(7110,7636,by=1),collapse = ',') ,',7863,7864,7865,7866,7867,7868,7869,7870,7871,7872,7873,7874,7875,7876,7877,7881,7882,7883,7884,7885,7886,', paste(seq(9400,9719,by=1),collapse=','))
  }else{
    coords <- paste0('823,824,825,892,893,894,907,908,909,1012,1013,1014,1156,1157,1158,1384,1385,1386,1444,1445,1446,1930,1931,1932,1957,1958,1959,2014,2015,2016,2023,2024,2025,2080,2081,2082,2134,2135,2136,2191,2192,2193,2280,2281,2282,2283,2284,2285,2298,2299,2300,2310,2311,2312,2316,2317,2318,2319,2320,2321,2322,2323,2324,2340,2341,2342,2346,2347,2348,2349,2350,2351,2352,2353,2354,2355,2356,2357,2358,2359,2360,2373,2374,2375,2379,2380,2381,2385,2386,2387,2388,2389,2390,2391,2392,2393,2394,2395,2396,2400,2401,2402,2409,2410,2411,2412,2413,2414,2415,2416,2417,2424,2425,2426,2430,2431,2432,2436,2437,2438,2439,2440,2441,2442,2443,2444,2457,2458,2459,2460,2461,2462,2463,2464,2465,2469,2470,2471,2472,2473,2474,2478,2479,2480,2481,2482,2483,2496,2497,2498,2499,2500,2501,2502,2503,2504,2505,2506,2507,2514,2515,2516,2517,2518,2519,2520,2521,2522,2526,2527,2528,2529,2530,2531,2535,2536,2537,2670,2671,2672,2679,2680,2681,2703,2704,2705,2709,2710,2711,2733,2734,2735,2742,2743,2744,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2769,2770,2771,2772,2773,2774,2778,2779,2780,2811,2812,2813,2814,2815,2816,2817,2818,2819,2823,2824,2825,2841,2842,2843,2847,2848,2849,2850,2851,2852,2856,2857,2858,2865,2866,2867,2871,2872,2873,2892,2893,2894,2895,2896,2897,2901,2902,2903,2904,2905,2906,2952,2953,2954,2961,2962,2963,3000,3001,3002,3015,3016,3017,3018,3019,3020,3030,3031,3032,3042,3043,3044,3084,3085,3086,3090,3091,3092,3099,3100,3101,3111,3112,3113,3117,3118,3119,3135,3136,3137,3171,3172,3173,3177,3178,3179,3180,3181,3182,3189,3190,3191,3192,3193,3194,3204,3205,3206,3210,3211,3212,3222,3223,3224,3228,3229,3230,3237,3238,3239,3246,3247,3248,3249,3250,3251,3255,3256,3257,3261,3262,3263,3396,3397,3398,3501,3502,3503,3546,3547,3548,3705,3706,3707,4425,4426,4427,4449,4450,4451,4503,4504,4505,4518,4519,4520,4590,4591,4592,4641,4642,4643,4647,4648,4649,4656,4657,4658,4668,4669,4670,4671,4672,4673,4692,4693,4694,4722,4723,4724,4782,4783,4784,4974,4975,4976,5016,5017,5018,5067,5068,5069,7863,7864,7865,7866,7867,7868,7869,7870,7871,7872,7873,7874,7875,7876,7877,7881,7882,7883,7884,7885,7886',sep='')
  }
  writeLines(paste0(coords), file.path(args$out.dir.work, 'excision_coords.txt'))
}

if(file.exists(cmds.path))
{
  # Load commands
  cat('Cleaning work directory...\n')
  move.logs(args$out.dir.work)
  pty.c <- readRDS(cmds.path)

}else{

  # Copy files into input folder
  tmp  <- file.path(dir.data, 'PANGEA2_RCCS/phyloscanner_input_data')
  tmp1 <- file.path(dir.data, 'PANGEA2_RCCS/220419_reference_set_for_PARTNERS_mafft.fasta')
  tmp <- c(list.files( tmp, full.names=T), tmp1)
  file.copy(
    tmp,  args$out.dir.data,
    overwrite = T,
    recursive = FALSE,
    copy.mode = TRUE
  )

  # Set consensus sequences.
  # (default consensus/reference for tsi and pair analyses if no arg is passed)
  # not really sure whether oneeach is needed anywhere

  cat("Load consensus sequences\n")
  infile.consensus <- args$reference

  if(is.na(infile.consensus))
  {
    # if NA, look in args$out.dir.data: TODO:check below
    tmp <- ifelse(args$tsi_analysis,
                  yes='220419_reference_set_for_PARTNERS_mafft.fasta',
                  no='ConsensusGenomes.fasta')
    infile.consensus <- file.path(args$out.dir.data, tmp)
  }else{
    # If does not exists, look within args$out.dir.data
    if(! file.exists(infile.consensus) )
    {
      infile.consensus <- file.path(args$out.dir.data, basename(infile.consensus))
    }
  }
  stopifnot(file.exists(infile.consensus))

  cat("Load sequences and remove duplicates if existing...\n")
  pty.runs <- data.table(readRDS(args$infile.runs))

  if (args$dryrun){
    cat("Dry run: only running for first 2 ptyr\n")
    pty.runs <- subset(pty.runs, PTY_RUN %in% c(1,2))
  }

  if ('ID_TYPE' %in% colnames(pty.runs)) {
    setorder(pty.runs, PTY_RUN, -ID_TYPE, ID)
  } else{
    setorder(pty.runs, PTY_RUN, ID)
  }

  tmp <- pty.runs[, list(idx=duplicated(SAMPLE_ID)), by = 'PTY_RUN']
  tmp <- tmp[, which(idx)]
  if(length(tmp)!=0){
    pty.runs <- pty.runs[-tmp, ]
  }

  tmp <- pty.runs[, uniqueN(SAMPLE_ID) == .N, by='PTY_RUN' ]
  stopifnot( all(tmp$V1) )

  cat("Load backgrounds and extract HXB2 for pairwise MSA... \n")
  consensus_seq <- seqinr::read.fasta(infile.consensus)
  consensus_seq_names <- names(consensus_seq)
  hxb2 <- grep('HXB2', names(consensus_seq), value = T)
  hxb2_seq <- consensus_seq[[hxb2]]
  # removed root_seq as it seemed useless

  # Change the format
  if ('ID_TYPE' %in% colnames(pty.runs)) {
    pty.runs[ID_TYPE == 'control', UNIT_ID := paste0('CNTRL-', UNIT_ID)]
    pty.runs[ID_TYPE == 'control', RENAME_ID := paste0('CNTRL-', RENAME_ID)]
  }
  pty.runs[, BAM := paste0(dir.data, SAMPLE_ID, '.bam')]
  pty.runs[, REF := paste0(dir.data, SAMPLE_ID, '_ref.fasta')]
  setkey(pty.runs, PTY_RUN, RENAME_ID)

  cat("Set the alignment options...\n")
  # MAFFT: reformat options so they are readily pasted in sh command.
  args$mafft.opt <- gsub('mafft|"', '', args$mafft.opt)
  args$mafft.opt <- paste0('"mafft ', args$mafft.opt, '"')

  # GENOMIC WINDOWS
  # excision.default will excise more positions, atm I group together with remove vloops
  stopifnot(args$windows_start <= args$window_end)
  ptyi <- seq(args$windows_start, args$windows_end, by=args$sliding_width)

  excision.default.bool <- args$rm_vloops
  if(excision.default.bool)
  {
    cat("Remove vloops...\n")
    ptyi <- c(ptyi[ptyi <= 6615 - args$window_size], 6825, 6850, ptyi[ptyi >= 7636])
  }

  # Now write command
  write.pty.command <- function(w_from) {
    pty.args <- list(
      prog.pty = file.path(args$prog.dir, "phyloscanner_make_trees.py"),
      prog.mafft = args$mafft.opt,
      data.dir = args$out.dir.data,
      work.dir = args$out.dir.work,
      out.dir = args$out.dir.output,
      alignments.file = infile.consensus,
      # alignments.root = root_seq, # is this even doing anything?
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
      win = c(w_from, w_from + args$window_size, args$sliding_width, args$window_size),
      keep.overhangs = FALSE,
      mem.save = 0,
      verbose = TRUE,
      select = NA,
      default.coord = TRUE,
      realignment = TRUE,
      excision.coords = "$EXCISION_COORDS"  # Assuming $EXCISION_COORDS is defined in the shell environment
    )
    pty.c <- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)
    pty.c[, W_FROM := w_from ]
  }

  pty.c	<- lapply(ptyi, write.pty.command)
  pty.c	<- do.call('rbind', pty.c)
  setkey(pty.c, PTY_RUN, W_FROM)
  stopifnot(pty.c[, .N, by=c('PTY_RUN', 'W_FROM')][, all(N == 1)])

  pty.c[ , OUTDIR := file.path(args$out.dir.output, paste0('ptyr', PTY_RUN, '_trees')) ]
  pty.c[, OUT_REGEX_V1 := paste0( W_FROM, '_to_', W_FROM + args$window_size - 1,  '.fasta') ]
  pty.c[, OUT_REGEX_V2 := paste0( W_FROM, '_to_', W_FROM + args$window_size - 1,  '_v2.fasta') ]

  saveRDS(pty.c, cmds.path)
}

# Now check which outputs are not ready yet and need re-running
cols <- c('OUT1', 'OUT2')
pty.c[, (cols) := lapply(.SD, .check.existing.outputs, outdir=OUTDIR, pattern='fasta' ) , by=OUTDIR, .SDcols= names(pty.c) %like% 'REGEX']

# print statments documenting completed things
pty.c[, cat(sum(!is.na(OUT2)),
            '(', .percent(mean(!is.na(OUT2))),')',
            'of desired outputs were generated in the previous runs\n')]
pty.c[, cat(sum(!is.na(OUT1)),
            '(', .percent(mean(!is.na(OUT1))),')',
            'of intermediary outputs were generated in the previous runs\n')]

# Select jobs with missing final outcome.
# remove problematic windows as informed by the shell scripts output
pty.c <- pty.c[ is.na(OUT2)]
pty.c <- .remove_problematic_windows(pty.c)

# if intermediary output is there, change CMD.
any_existing_v1 <- pty.c[ !is.na(OUT1), .N > 0]
if( any_existing_v1)
{
  pty.c[ ! is.na(OUT1) , .modify_cmd_for_existing_v1(cmd=CMD, outdir=OUTDIR, out1=OUT1)  , by='OUT1' ]
}

#
# IF ALL OUTPUTS WERE CREATED, GO TO NEXT STEP IN THE ANALYSIS
#

if( nrow(pty.c) == 0 )
{
  # ISN'T THIS BEAUTIFUL?
  qsub.next.step(file=args$controller,
                 next_step='btr',
                 res=1,
                 redo=0
  )
  stop('Alignment step completed, submitted the following task')
}

#
# NOW THAT WE HAVE THE CMDs to RUN, WE NEED TO AGGREGATE THEM
#

# print  the first cmd to console for checking purposes
cat(pty.c$CMD[1])

# aggregate jobs and run
n_jobs <- ceiling( nrow(pty.c) / max.per.run)
idx <- 1:nrow(pty.c)

pty.c[, CASE_ID := rep(1:max.per.run, times = n_jobs)[idx] ]
pty.c[, JOB_ID := rep(1:n_jobs, each = max.per.run)[idx] ]

# Write and submit:
djob <- pty.c[, .(CMD=.write.job(.SD,path=args$out.dir.work)), by=JOB_ID]
if(args$walltime_idx == 3 & args$pqeelab){
  # Assign one every 5 jobs to pqeelab
  djob[, Q := {
    idx <- seq(5, .N, 5)
    z <- rep('', .N)
    z[idx] <- 'pqeelab'
    z
  }]
}

# save input_list for incomplete jobs
ephemeral_dir <- Sys.getenv("EPHEMERAL")
file_name <- "input_list.txt"
#file_path <- file.path(ephemeral_dir, file_name)
file_path <- file.path(args$out.dir.work, file_name)
input_list <- pty.c[, c('PTY_RUN','W_FROM')]
input_list[, PTY_RUN:= paste0('ptyr', PTY_RUN)]
input_list[, WINDOW:= paste(W_FROM,W_FROM + args$window_size - 1, sep=',')]
input_list <- input_list[, c('PTY_RUN','WINDOW')]
fwrite(input_list, file_path, sep = " ", col.names=FALSE)
cat("File saved to:", file_path, "\n")

# Load runners specifications
split_jobs_by_n <- 1
if(!is.na(args$runners))
{
  drunners <- fread(args$runners)
  split_jobs_by_n <- nrow(drunners)
}

if (split_jobs_by_n == 1 | nrow(pty.c) <= 1000 ){
  ids <- submit_jobs_from_djob(djob, output_type = "id")
  # qsub alignment step again,
  # to check whether everything has run...
  qsub.next.step(file=args$controller,
                 ids=ids,
                 next_step='ali',
                 res=args$walltime_idx + 1,
                 redo=0
  )
}else{
  # split djob in split_jobs_by_n data.tables
  # then perform the above individually
  person_to_run <- rep(1:split_jobs_by_n, length.out = nrow(pty.c))
  submit_user_script <- NA_character_

  for (person in 1:split_jobs_by_n){

    # Subset to jobs for specific person
    djob_person <- pty.c[person_to_run == person]
    djob_person[,CMD:= djob[1,CMD]]
    # Adapt the job specifications according to the person/runner
    djob_person <- adapt_jobspecs_to_runner(
      djob_person,
      drunners,
      djob,
      idx = person
    )
    # write the pbs files
    djob_person <- djob_person[, .(CMD = CMD[1]), by = 'JOB_ID']

    pbs_file_person <- submit_jobs_from_djob(djob_person, output_type = "outfile")
    # Append the pbs files to the script that each user can submit to queue them
    submit_user_script <- append_pbs_file_person(
      script = submit_user_script,
      pbs = pbs_file_person,
      usr = drunners[index == person,  user_name]
    )
  }

  # Write the script that each user can submit
  outfile = file.path(args$out.dir.work, 'submit_user_readali.sh')
  write( submit_user_script, file = outfile)

}

cat("End of script\n")
