cat("\n\n===== TSI_initialise.R =====\n\n")

require(data.table)

option_list <- list(
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
    "--cluster_size",
    type = "integer",
    default = 50L,
    help = "Minimum cluster size for host ids groupings[default %default]",
    dest = "cluster_size"
  ),
  optparse::make_option(
    "--controller",
    type = "character",
    default = NA_character_, 
    help = "Path to sh script directing the full analysis",
    dest = 'controller'
  ),
  optparse::make_option(
    "--transmission_chains",
    type = "character",
    default = NA_character_,
    help = "Optional: absolute file path to `phscnetwork.rda` containing individuals in potential transmission pairs",
    dest = 'file.path.chains'
  ),
  optparse::make_option(
    "--include_input",
    type = "character",
    default = NA_character_,
    help = "Optional: path to phscinput*rds file of individuals to include in the analysis (eg seroconverters)",
    dest = 'include.input'
  )
)

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

usr <- Sys.info()[['user']]
if (usr == 'andrea')
{
        args$pkg.dir <- '~/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/'
        indir.deepsequence.analyses <- '~/Documents/Box/ratmann_deepseq_analyses/live'
        indir.deepsequence.xiaoyue <- '~/Documents/Box/ratmann_xiaoyue_jrssc2022_analyses/live'
        indir.deepsequencedata <- '~/Documents/Box/ratmann_pangea_deepsequencedata'
}else{
        args$pkg.dir <- '~/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/'
        indir.deepsequence.analyses <- '/rds/general/project/ratmann_deepseq_analyses/live'
        indir.deepsequence.xiaoyue <- "/rds/general/project/ratmann_xiaoyue_jrssc2022_analyses/live"
        indir.deepsequencedata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live'
}
indir.deepsequence.analyses.old <- file.path(indir.deepsequence.xiaoyue, 'PANGEA2_RCCS1519_UVRI')
file.phsc.input.samples.bf<- file.path(indir.deepsequence.analyses.old, '220331_RCCSUVRI_phscinput_samples_with_bf.rds' )

tmp <- c(indir.deepsequence.analyses.old,
         file.phsc.input.samples.bf,
         file.path.chains)
if( ! is.na(args$file.path.chains) ) tmp <- c(tmp, args$file.path.chaings)
if( ! is.na(args$include.input) ) tmp <- c(tmp, args$include.input)
stopifnot(all(file.exists(tmp)))


###
# Main 
###

# Create ouput directory
if(! dir.exists(args$out.dir))
{
        cat('Creating output directory\n')
        dir.create(args$out.dir)
        dir.create(file.path(args$out.dir, 'potential_network'))
}else{
        cat('Warning: output directory already exists\n')
}

# load individuals in potential transmission pairs if given
#__________________________________________________________

if( ! is.na(args$file.path.chains) )
{
        tmp <- new.env()
        load(args$file.path.chains, envir=tmp )
        dchain <- as.data.table(tmp$dchain)
        rm(tmp)
        include_pairs_aid <- dchain[, unique(c(H1, H2))]
}

if( ! is.na(args$include.input))
{
        tmp <- readRDS(args$include.input)
        include_rename_id <- unique(tmp$RENAME_ID)
        rm(tmp)
}

# Load old phsc input samples, and subset as required
#____________________________________________________

phsc_samples <- readRDS(file.phsc.input.samples.bf)
phsc_samples[, AID := gsub('-fq[0-9]+$', '', RENAME_ID)]
phsc_samples[, INCLUDE := TRUE]

if( ! is.na(args$file.path.chains))
{
        stopifnot( all(include_pairs_aid %in% phsc_samples$AID ))
        phsc_samples[! AID %in% include_pairs_aid, INCLUDE := FALSE]
}
if( ! is.na(args$include.input))
{
        stopifnot( all(include_rename_id %in% phsc_samples$RENAME_ID ))
        phsc_samples[ RENAME_ID %in% include_rename_id, INCLUDE := TRUE]
}

filename=file.path(args$out.dir, basename(file.phsc.input.samples.bf))
if( ! is.na(args$file.path.chains) | ! is.na(args$include.input) )
{
        tmp <- basename(filename)        
        date <- format(Sys.Date(), '%y%m%d')
        tmp <- gsub('^[0-9]+', date, tmp)
        tmp <- gsub('\\.rds$', '_subset.rds', tmp)
        filename <- file.path(dirname(filename), tmp)
}

phsc_samples <- phsc_samples[INCLUDE == TRUE]
phsc_samples[, `:=` (RENAME_ID=NULL, INCLUDE=NULL)]
saveRDS(phsc_samples, filename)


# Make clusters.rds
# ______________________________
set.seed(42)
NCLU = phsc_samples[, floor(uniqueN(UNIT_ID) / args$cluster_size) ]
dclus <- phsc_samples[, list(ID=sample(UNIT_ID), IDCLU= 1:NCLU)]
dclus[, CLU_SIZE:=.N, by='IDCLU']
setkey(dclus, IDCLU)
filename=file.path(args$out.dir, 'potential_network', 'clusters.rds')
saveRDS(dclus, filename)


# make phscinput_runs_clusize...
# ______________________________
suffix <- phsc_samples[, .(UNIT_ID,PANGEA_ID, RENAME_ID, SAMPLE_ID)]
dclus[, `:=` (PTY_RUN=IDCLU,  PTY_SIZE=CLU_SIZE) ]
dclus <- merge(dclus, suffix, by.x='ID', by.y='UNIT_ID')
setnames(dclus, 'ID', 'UNIT_ID')
filename=file.path(args$out.dir,
                   paste0('phscinput_runs_clusize_', args$cluster_size,'_ncontrol_0.rds'))
setkey(dclus, IDCLU)
saveRDS(dclus, filename)


# Source functions
source(file.path(args$pkg.dir, "utility.R"))
qsub.next.step(file=args$controller,
               next_step='ali', 
               res=1, 
               redo=0
)
