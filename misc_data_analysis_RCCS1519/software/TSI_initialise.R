cat("\n\n===== TSI_initialise.R =====\n\n")

require(data.table)
require(ggplot2)

option_list <- list(
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
  )
)

args <-
  optparse::parse_args(optparse::OptionParser(option_list = option_list))

usr <- Sys.info()[['user']]
if (usr == 'andrea')
{
        indir.deepsequence_analyses <- '~/Documents/Box/ratmann_deepseq_analyses/live'
        indir.deepsequencedata <- '~/Documents/Box/ratmann_pangea_deepsequencedata'
        tanya.rakai.dir <- '~/git/HIV-phyloTSI-main/RakExample_Tanya'
}else{
        indir.deepsequence_analyses <- '/rds/general/project/ratmann_deepseq_analyses/live'
        indir.deepsequencedata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live'
        tanya.rakai.dir <- '~/git/HIV-phyloTSI/RakExample_Tanya'
}


# other paths 
indir.deepsequence_analyses_old <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS1519_UVRI')
file.db.sharing <- file.path(indir.deepsequencedata,"/PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv")
file.anonymisation.keys <- file.path(indir.deepsequence_analyses_old,'important_anonymisation_keys_210119.csv')
file.phsc.input.samples.bf<- file.path(indir.deepsequence_analyses_old, '220331_RCCSUVRI_phscinput_samples_with_bf.rds' )

###
# Main 
###

# 
if(! dir.exists(args$out.dir))
{
        cat('Creating output directory\n')
        dir.create(args$out.dir)
}else{
        cat('Warning: output directory already exists\n')
}

# copy phsc_input_samples from old directory to new:
phsc_samples <- readRDS(file.phsc.input.samples.bf)
filename=file.path(args$out.dir, basename(file.phsc.input.samples.bf))
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
                   paste0('phscinput_runs_clusize_',max(dclus[,CLU_SIZE]),'_ncontrol_0.rds'))
setkey(dclus, IDCLU)
saveRDS(dclus, filename)

if(file.exists(args$controller))
{
        cat('\nTry to set up next step step:\n')

        cmd <- paste0('cd ', dirname(args$controller),'\n')
        cmd <- paste0(cmd, 'qsub -v STEP="ali" ',args$controller, '\n')
        cat(system(cmd), TRUE)
}
