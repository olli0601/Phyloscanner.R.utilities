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
  ),
  optparse::make_option(
    "--tsi_analysis",
    type = 'logical',
    default = FALSE,  
    help = 'Indicator on whether we want to perform a Time Since Infection analysis',
    dest = 'tsi_analysis'
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

file.copy(
          file.path(dir.data,'PANGEA2_RCCS/220419_reference_set_for_PARTNERS_mafft.fasta'),
          args$out.dir.data,
          overwrite = T,
          recursive = FALSE,
          copy.mode = TRUE
)

# Set consensus sequences.

if(args$tsi_analysis)
{
        infile.consensus <- file.path(args$out.dir.data,'220419_reference_set_for_PARTNERS_mafft.fasta')
}else{
        infile.consensus <-
          file.path(args$out.dir.data, 'ConsensusGenomes.fasta')
}
infile.consensus.oneeach <-  file.path(args$out.dir.data, '2019_New_ConsensusGenomesOneEach_GeneCut.fasta')


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

args <- list(
        window_size=250,
        sliding_width=10
)

# TODO: if standard:
if(args$tsi_analysis)
{
        ptyi <- (52:949)*10
        mafft.opt <- '\" mafft \"'
        excision.default.bool <- FALSE
}else{
        ptyi <- seq(800, 9175, args$sliding_width)
        ptyi <- c(ptyi[ptyi <= 6615 - args$window_size], 6825, 6850, ptyi[ptyi >= 7636])
        mafft.opt <- '\" mafft --globalpair --maxiterate 1000 \" ',
        excision.default.bool <- TRUE
}

if(0)
{
        infile.consensus <- root_seq <- hxb2 <- 'LALALA'
        phsc.cmd.phyloscanner.multi
}
# function(pty.runs, pty.args) 		
# {
#   stopifnot( any(colnames(pty.runs)=='PTY_RUN') )
#   stopifnot( any(colnames(pty.runs)=='SAMPLE_ID') )
#   stopifnot( any(colnames(pty.runs)=='BAM') )
#   stopifnot( any(colnames(pty.runs)=='REF') )	
#   ptyd	<- copy(pty.runs)		
#   if(!any(is.na(pty.args[['select']])))
#     ptyd<- subset(ptyd, PTY_RUN%in%pty.args[['select']])			
#   #	if background alignment is specified in pty.runs, use it
#   if(any(colnames(ptyd)=='BACKGROUND_ID'))
#   {
#     tmp	<- unique(data.table(BACKGROUND_ID=ptyd[, BACKGROUND_ID], FILE=file.path(pty.args[['data.dir']], ptyd[, BACKGROUND_ID])))
#     #	check if files exist
#     tmp	<- tmp[, list(EXISTS=file.exists(FILE)), by=c('BACKGROUND_ID','FILE')]
#     if(any(tmp[, !EXISTS]))
#       warning('\nCould not find location of BACKGROUND files for all runs in pty.runs, n=', tmp[, length(which(!EXISTS))],'\nRuns with missing background alignments are ignored. Please check.')
#     ptyd <- merge(ptyd, subset(tmp, EXISTS, c(BACKGROUND_ID, FILE)), by='BACKGROUND_ID')
#     set(ptyd, NULL, 'BACKGROUND_ID', NULL)
#     setnames(ptyd, 'FILE', 'BACKGROUND_ID')
#   }
#   setkey(ptyd, PTY_RUN)
#   #	add legacy arguments
#   pty.args['all.bootstrap.trees']	<- FALSE
#   pty.args['num.bootstraps']	<- 0		
#   #	write pty.run files and get pty command lines
#   pty.c		<- ptyd[, {
#     #	PTY_RUN<- z <- 1; BAM<- subset(ptyd, PTY_RUN==z)[, BAM]; REF<- subset(ptyd, PTY_RUN==z)[, REF]
#     #	SAMPLE_ID<- subset(ptyd, PTY_RUN==z)[, SAMPLE_ID]; RENAME_ID<- subset(ptyd, PTY_RUN==z)[, RENAME_ID]; BACKGROUND_ID<- subset(ptyd, PTY_RUN==z)[, BACKGROUND_ID]
#     file.input		<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_input.csv',sep=''))
#     tmp				<- cbind(BAM[!is.na(BAM)&!is.na(REF)], REF[!is.na(BAM)&!is.na(REF)])
#     if(exists('RENAME_ID'))
#       tmp			<- cbind(tmp, RENAME_ID[!is.na(BAM)&!is.na(REF)])
#     write.table(tmp, file=file.input, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')
#     file.patient	<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_patients.txt',sep=''))
#     tmp				<- paste0(SAMPLE_ID[!is.na(BAM)&!is.na(REF)],'.bam')
#     if(exists('RENAME_ID'))
#       tmp			<- unique(RENAME_ID[!is.na(BAM)&!is.na(REF)])
#     if(exists('UNIT_ID'))
#       tmp			<- unique(UNIT_ID[!is.na(BAM)&!is.na(REF)])				
#     cat( paste(tmp,collapse='\n'), file= file.patient	)
#     if(exists('BACKGROUND_ID'))
#     {
#       if(length(unique(BACKGROUND_ID))!=1)
#         stop('\nrun',PTY_RUN,' Expected one background alignment file, found: ',paste(unique(BACKGROUND_ID), collapse=''),'. Please check.')
#       pty.args$alignments.file<- unique(BACKGROUND_ID)
#     }				
#     cmd			<- phsc.cmd.phyloscanner.one(	pty.args, 
#                                         file.input, 
#                                         file.patient)
#     #cmd			<- paste(cmd, pty.cmd.evaluate.fasta(pty.args[['out.dir']], strip.max.len=pty.args[['strip.max.len']], select=paste('^ptyr',PTY_RUN,'_In',sep=''), min.ureads.individual=pty.args[['min.ureads.individual']]), sep='')
#     #cat(cmd)
#     list(CMD= cmd)				
#   },by='PTY_RUN']
#   pty.c
# }
# <bytecode: 0x55f1cb642688>


###
#phsc.cmd.phyloscanner.one
###

# function(pty.args, file.input, file.patient)
# {	
#   stopifnot(is.character(file.input),is.character(file.patient))
#   #	copy input variables into namespace	 
#   attach(pty.args)	
#   #	sense checks
#   stopifnot(is.logical(all.bootstrap.trees))	
#   #	define window coordinates
#   window.coord			<- integer(0)
#   if(!nchar(pty.args[['window.automatic']]))
#   {
#     stopifnot(length(pty.args[['win']])==4)
#     window.coord		<- seq(pty.args[['win']][1], pty.args[['win']][2]-max(pty.args[['win']][3:4]),pty.args[['win']][3])
#     window.coord		<- as.vector(rbind( window.coord,window.coord-1+pty.args[['win']][4] ))											
#   }	
#   #	
#   if(!nchar(window.automatic))	stopifnot( is.numeric(window.coord), !length(window.coord)%%2)
#   if(nchar(window.automatic))		stopifnot( !length(window.coord) )
#   merge.paired.reads			<- ifelse(!is.na(merge.paired.reads) & merge.paired.reads, '--merge-paired-reads', NA_character_)
#   keep.overhangs				<- ifelse(!is.na(keep.overhangs) & keep.overhangs, '--keep-overhangs', NA_character_)
#   no.trees					<- ifelse(!is.na(no.trees) & no.trees, '--no-trees', NA_character_)
#   num.bootstraps				<- ifelse(is.na(no.trees) & !is.na(num.bootstraps) & is.numeric(num.bootstraps) & num.bootstraps>1, as.integer(num.bootstraps), NA_integer_)
#   dont.check.duplicates		<- ifelse(!is.na(dont.check.duplicates) & dont.check.duplicates, '--dont-check-duplicates', NA_character_)
#   dont.check.recombination	<- ifelse(!is.na(dont.check.recombination) & dont.check.recombination==FALSE, '--check-recombination', NA_character_)	
#   #	create local tmp dir
#   cmd		<- paste("CWD=$(pwd)\n",sep='\n')
#   cmd		<- paste(cmd,"echo $CWD\n",sep='')
#   tmpdir	<- paste('pty','_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')	
#   tmpdir	<- paste("$CWD/",tmpdir,sep='')
#   cmd		<- paste(cmd,'mkdir -p "',tmpdir,'"\n',sep='')
#   #	copy files to local tmp dir
#   cmd		<- paste(cmd,'cp "',file.input,'" "',tmpdir,'"\n',sep='')	
#   cmd		<- paste(cmd,'cp "',file.patient,'" "',tmpdir,'"\n',sep='')
#   #	cd to tmp dir
#   cmd		<- paste(cmd, 'cd "',tmpdir,'"\n', sep='')	
#   cmd		<- paste(cmd, prog.pty,' "',basename(file.input),'" ',sep='')	
#   cmd		<- paste(cmd, '--merging-threshold-a', merge.threshold,'--min-read-count',min.read.count,'--quality-trim-ends', quality.trim.ends, '--min-internal-quality',min.internal.quality,'--keep-output-together')
#   if(!is.na(merge.paired.reads))
#     cmd	<- paste(cmd,' ',merge.paired.reads,sep='')	
#   if(!is.na(dont.check.duplicates))
#     cmd	<- paste(cmd,' ',dont.check.duplicates,sep='')
#   if(!is.na(dont.check.recombination))
#     cmd	<- paste(cmd,' ',dont.check.recombination,sep='')	
#   if(!is.na(alignments.pairwise.to))		
#     cmd	<- paste(cmd,' --pairwise-align-to ',alignments.pairwise.to,sep='')
#   if(!is.na(default.coord)){cmd <- paste(cmd, " --excision-ref ", alignments.pairwise.to,sep='')}
#   if(nchar(window.automatic))
#     cmd	<- paste(cmd,' --auto-window-params ', window.automatic,sep='')
#   if(!nchar(window.automatic))
#     cmd	<- paste(cmd,' --windows ', paste(as.character(window.coord), collapse=','),sep='')
#   if(!is.na(alignments.file))
#     cmd	<- paste(cmd,' --alignment-of-other-refs "',alignments.file,'"',sep='')	
#   if(!is.na(no.trees))		
#     cmd	<- paste(cmd, no.trees)
#   if(!is.na(num.bootstraps))
#     cmd	<- paste(cmd, ' --num-bootstraps ', num.bootstraps, sep='')
#   if(!is.na(discard.improper.pairs) & discard.improper.pairs==TRUE)
#     cmd	<- paste(cmd, ' --discard-improper-pairs ', sep = '')
#   if(!is.na(keep.overhangs))
#     cmd	<- paste(cmd, keep.overhangs)
#   cmd		<- paste(cmd, '--x-mafft',prog.mafft)
#   if(is.na(no.trees))
#     cmd	<- paste(cmd, '--x-raxml',prog.raxml)	
#   if(!is.na(default.coord)){
#     if(default.coord){	  
#       cmd	<- paste0(cmd, " --excision-coords \'823,824,825,892,893,894,907,908,909,1012,1013,1014,1156,1157,1158,1384,1385,1386,1444,1445,1446,1930,1931,1932,1957,1958,1959,2014,2015,2016,2023,2024,2025,2080,2081,2082,2134,2135,2136,2191,2192,2193,2280,2281,2282,2283,2284,2285,2298,2299,2300,2310,2311,2312,2316,2317,2318,2319,2320,2321,2322,2323,2324,2340,2341,2342,2346,2347,2348,2349,2350,2351,2352,2353,2354,2355,2356,2357,2358,2359,2360,2373,2374,2375,2379,2380,2381,2385,2386,2387,2388,2389,2390,2391,2392,2393,2394,2395,2396,2400,2401,2402,2409,2410,2411,2412,2413,2414,2415,2416,2417,2424,2425,2426,2430,2431,2432,2436,2437,2438,2439,2440,2441,2442,2443,2444,2457,2458,2459,2460,2461,2462,2463,2464,2465,2469,2470,2471,2472,2473,2474,2478,2479,2480,2481,2482,2483,2496,2497,2498,2499,2500,2501,2502,2503,2504,2505,2506,2507,2514,2515,2516,2517,2518,2519,2520,2521,2522,2526,2527,2528,2529,2530,2531,2535,2536,2537,2670,2671,2672,2679,2680,2681,2703,2704,2705,2709,2710,2711,2733,2734,2735,2742,2743,2744,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2769,2770,2771,2772,2773,2774,2778,2779,2780,2811,2812,2813,2814,2815,2816,2817,2818,2819,2823,2824,2825,2841,2842,2843,2847,2848,2849,2850,2851,2852,2856,2857,2858,2865,2866,2867,2871,2872,2873,2892,2893,2894,2895,2896,2897,2901,2902,2903,2904,2905,2906,2952,2953,2954,2961,2962,2963,3000,3001,3002,3015,3016,3017,3018,3019,3020,3030,3031,3032,3042,3043,3044,3084,3085,3086,3090,3091,3092,3099,3100,3101,3111,3112,3113,3117,3118,3119,3135,3136,3137,3171,3172,3173,3177,3178,3179,3180,3181,3182,3189,3190,3191,3192,3193,3194,3204,3205,3206,3210,3211,3212,3222,3223,3224,3228,3229,3230,3237,3238,3239,3246,3247,3248,3249,3250,3251,3255,3256,3257,3261,3262,3263,3396,3397,3398,3501,3502,3503,3546,3547,3548,3705,3706,3707,4425,4426,4427,4449,4450,4451,4503,4504,4505,4518,4519,4520,4590,4591,4592,4641,4642,4643,4647,4648,4649,4656,4657,4658,4668,4669,4670,4671,4672,4673,4692,4693,4694,4722,4723,4724,4782,4783,4784,4974,4975,4976,5016,5017,5018,5067,5068,5069,", paste(seq(6615,6811,by=1),collapse = ','),  ",", paste(seq(7110,7636,by=1),collapse = ',') ,",7863,7864,7865,7866,7867,7868,7869,7870,7871,7872,7873,7874,7875,7876,7877,7881,7882,7883,7884,7885,7886,", paste(seq(9400,9719,by=1),collapse=','),"\' ")
#     }else{
#       cmd	<- paste(cmd, " --excision-coords \'823,824,825,892,893,894,907,908,909,1012,1013,1014,1156,1157,1158,1384,1385,1386,1444,1445,1446,1930,1931,1932,1957,1958,1959,2014,2015,2016,2023,2024,2025,2080,2081,2082,2134,2135,2136,2191,2192,2193,2280,2281,2282,2283,2284,2285,2298,2299,2300,2310,2311,2312,2316,2317,2318,2319,2320,2321,2322,2323,2324,2340,2341,2342,2346,2347,2348,2349,2350,2351,2352,2353,2354,2355,2356,2357,2358,2359,2360,2373,2374,2375,2379,2380,2381,2385,2386,2387,2388,2389,2390,2391,2392,2393,2394,2395,2396,2400,2401,2402,2409,2410,2411,2412,2413,2414,2415,2416,2417,2424,2425,2426,2430,2431,2432,2436,2437,2438,2439,2440,2441,2442,2443,2444,2457,2458,2459,2460,2461,2462,2463,2464,2465,2469,2470,2471,2472,2473,2474,2478,2479,2480,2481,2482,2483,2496,2497,2498,2499,2500,2501,2502,2503,2504,2505,2506,2507,2514,2515,2516,2517,2518,2519,2520,2521,2522,2526,2527,2528,2529,2530,2531,2535,2536,2537,2670,2671,2672,2679,2680,2681,2703,2704,2705,2709,2710,2711,2733,2734,2735,2742,2743,2744,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2769,2770,2771,2772,2773,2774,2778,2779,2780,2811,2812,2813,2814,2815,2816,2817,2818,2819,2823,2824,2825,2841,2842,2843,2847,2848,2849,2850,2851,2852,2856,2857,2858,2865,2866,2867,2871,2872,2873,2892,2893,2894,2895,2896,2897,2901,2902,2903,2904,2905,2906,2952,2953,2954,2961,2962,2963,3000,3001,3002,3015,3016,3017,3018,3019,3020,3030,3031,3032,3042,3043,3044,3084,3085,3086,3090,3091,3092,3099,3100,3101,3111,3112,3113,3117,3118,3119,3135,3136,3137,3171,3172,3173,3177,3178,3179,3180,3181,3182,3189,3190,3191,3192,3193,3194,3204,3205,3206,3210,3211,3212,3222,3223,3224,3228,3229,3230,3237,3238,3239,3246,3247,3248,3249,3250,3251,3255,3256,3257,3261,3262,3263,3396,3397,3398,3501,3502,3503,3546,3547,3548,3705,3706,3707,4425,4426,4427,4449,4450,4451,4503,4504,4505,4518,4519,4520,4590,4591,4592,4641,4642,4643,4647,4648,4649,4656,4657,4658,4668,4669,4670,4671,4672,4673,4692,4693,4694,4722,4723,4724,4782,4783,4784,4974,4975,4976,5016,5017,5018,5067,5068,5069,7863,7864,7865,7866,7867,7868,7869,7870,7871,7872,7873,7874,7875,7876,7877,7881,7882,7883,7884,7885,7886\' \\\" ",sep='')
#     }
#   }
#   cmd		<- paste(cmd, '\n')
#   run.id	<- gsub('_input.csv','',basename(file.input))
#   #	process RAxML files
#   if(is.na(no.trees) & (is.na(num.bootstraps) | (!is.na(num.bootstraps) & all.bootstrap.trees)))
#     cmd	<- paste(cmd, 'for file in RAxML_bestTree*.tree; do\n\tmv "$file" "${file//RAxML_bestTree\\./',run.id,'_}"\ndone\n',sep='')
#   if(is.na(no.trees) & !is.na(num.bootstraps) & !all.bootstrap.trees)
#     cmd	<- paste(cmd, 'for file in RAxML_bipartitions.MLtreeWbootstraps*.tree; do\n\tmv "$file" "${file//RAxML_bipartitions.MLtreeWbootstraps/',run.id,'_}"\ndone\n',sep='')	
#   #cmd	<- paste(cmd, "for file in AlignedReads*.fasta; do\n\tsed 's/<unknown description>//' \"$file\" > \"$file\".sed\n\tmv \"$file\".sed \"$file\"\ndone\n",sep='')		
#   if(!is.na(alignments.file) & !is.na(keep.overhangs))
#   {
#     cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\tcat "$file" | awk \'{if (substr($0,1,4) == ">REF") censor=1; else if (substr($0,1,1) == ">") censor=0; if (censor==0) print $0}\' > NoRef$file\ndone\n', sep='')		
#     cmd	<- paste(cmd, 'for file in NoRefAlignedReads*.fasta; do\n\t',phsc.cmd.mafft.add(alignments.file,'"$file"','Ref"$file"', options='--keeplength --memsave --parttree --retree 1'),'\ndone\n',sep='')		
#     cmd	<- paste(cmd, 'for file in RefNoRefAlignedReads*.fasta; do\n\t','mv "$file" "${file//RefNoRefAlignedReads/',run.id,'_}"\ndone\n',sep='')		
#   }
#   if(realignment==TRUE){
#     cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\t mafft --globalpair --maxiterate 1000 "$file" > "${file//.fasta/_v2.fasta}" \n done \n', sep='')
#   }
#   if(is.na(alignments.file) || is.na(keep.overhangs))
#   {
#     cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\tmv "$file" "${file//AlignedReads/',run.id,'_}"\ndone\n',sep='')	
#   }	
#   #	move Duplicate Read Counts - only for backward compatibility
#   #cmd	<- paste(cmd, 'for file in DuplicateReadCountsProcessed_*.csv; do\n\tmv "$file" "${file//DuplicateReadCountsProcessed_/',run.id,'_DuplicateReadCounts_}"\ndone',sep='')
#   #	run phyloscanner tools and compress output
#   if(is.na(no.trees))
#     cmd	<- paste(cmd, phsc.cmd.process.phyloscanner.output.in.directory(tmpdir, file.patient, pty.args), sep='\n')
#   #
#   out.dir2<- out.dir	
#   if(!is.na(no.trees))
#   {
#     out.dir2<- file.path(out.dir,paste0(run.id,'_trees'))
#     cmd		<- paste(cmd,'\nmkdir -p ',out.dir2)		
#   }
#   cmd		<- paste(cmd, '\ncp ',run.id,'* "',out.dir2,'"\n',sep='')	
#   #	zip up everything else
#   tmp		<- ''
#   if(length(window.coord)==2)
#     tmp	<- window.coord[1]
#   if(is.null(mem.save) || is.na(mem.save) || mem.save==0)
#   {
#     cmd		<- paste(cmd, 'for file in *; do\n\tzip -ur9XTjq ',paste(run.id,'_otherstuff',tmp,'.zip',sep=''),' "$file"\ndone\n',sep='')
#     cmd		<- paste(cmd, 'cp ',paste(run.id,'_otherstuff',tmp,'.zip',sep=''),' "',out.dir2,'"\n',sep='')		
#   }
#   #	clean up
#   cmd		<- paste(cmd,'cd $CWD\nrm -r "',tmpdir,'"\n',sep='')
#   detach(pty.args)
#   cmd
# }



pty.c	<- lapply(seq_along(ptyi), function(i)
{
  pty.args <- list(
    prog.pty = file.path(args$prog.dir, "phyloscanner_make_trees.py"),
    prog.mafft = mafft.opt,
    data.dir = args$out.dir.data,
    work.dir = args$out.dir.work,
    out.dir = args$out.dir.output,
    alignments.file = infile.consensus,
    alignments.root = root_seq, # is this even doing anything?
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
    default.coord = excision.default.bool,
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
