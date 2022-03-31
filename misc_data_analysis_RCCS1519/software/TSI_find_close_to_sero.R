cat("\n\n===== TSI_find_close_to_sero.R =====\n\n")

# The objective of this script is to find a small number of 
# participants that were close to some of the seroconverters
# for which we had "bad" Time Since Infection Estimates. We 
# can then group these and run phyloscanner on them. This is
# achieved via the successive scripts, but the make
# alignment script is called at the end of this one.

library(data.table)

option_list <- list(
  optparse::make_option(
    "--cluster_size",
    type = "integer",
    default = 6L,
    help = "Maximum cluster size [default %default]",
    dest = "cluster_size"
  ),
  optparse::make_option(
    "--out_dir_base",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where all output is stored [default]",
    dest = 'out.dir'
  ),
  optparse::make_option(
    "--sero_tsi",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to TSI results of seroconverters",
    dest = 'file.path.sero'
  ),
  optparse::make_option(
    "--path_similarities",
    type = "character",
    default = NA_character_,  
    help = "Absolute file path to average similarities across genome",
    dest = 'file.path.similarities' 
  ),
  optparse::make_option(
    "--controller",
    type = "character",
    default = NA_character_, # Think about adding the controller in the software directory
    help = "Path to sh script directing the full analysis",
    dest = 'controller'
  )
)

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

#
# add constants that should not be changed by the user
#
infile.run <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_phscinput_samples.rds'

infile.run2 <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/220331_RCCSUVRI_phscinput_samples_with_bf.rds'

# Testing 
user <- Sys.info()[['user']]
if(user == 'andrea')
{
        args$file.path.sero <- '/home/andrea/Documents/Box/ratmann_pangea_deepsequencedata/PANGEA2_RCCS/220329_TSI_seroconverters.csv'
        args$file.path.similarities <- '/home/andrea/Documents/Box/ratmann_deepseq_analyses/live/test_full/potential_network'
}

if(is.na(args$file.path.sero))
{
        args$file.path.sero <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/220329_TSI_seroconverters.csv'
        args$file.path.similarities <- '/rds/general/project/ratmann_deepseq_analyses/live/test_full/potential_network/whole_genome_average_similarity.rds'
}

# For simplicity, we want to pick individuals who
# only have one base frequency file
phsc.input <- readRDS(infile.run2)
dbf <- phsc.input[, .(study_id=UNIT_ID, 
                      PANGEA_ID=gsub('^.*?_','', PANGEA_ID),
                      BF_exists=!is.na(BF),
                      N=length(RENAME_ID)),
                    by='UNIT_ID']
dbf <- dbf[N == 1 & BF_exists==TRUE, ][, `:=` (N=NULL, BF_exists=NULL, UNIT_ID=NULL)]
bf_inds <- dbf[grepl('RK',study_id), study_id]

# load all genome similarities and only consider Rakai inds
dgensim <- setDT(readRDS(args$file.path.similarities))
dgensim <- dgensim[ pt_id1 %in% bf_inds & pt_id2 %in% bf_inds  ]
dgensim_ids <- dgensim[, unique(c(pt_id1, pt_id2))]

# load 10 'interesting' seroconverters for which we have bad predictions,
# have known TSI < 4 (E 8 years until symptoms), have BF files and
# genome similarities exist in above file (maybe look at average_similarity.rds instead).
dsero <- fread(args$file.path.sero)
dsero <- dsero[knownTSI < 4 & BF_exists==TRUE & study_id %in% dgensim_ids]
setorder(dsero, -'bias')
sero_ids <- dsero[1:10, study_id]
stopifnot(all(sero_ids %in% dbf$study_id))

# can slim this a bit before selecting close inds
dgensim <- dgensim[ pt_id1 %in% sero_ids | pt_id2 %in% sero_ids ]

dclusters <- data.table(ID=as.character(), IDCLU=as.integer(), CLU_SIZE=as.integer())
for (id in sero_ids)
{
        exclude_ids <-  sero_ids[ sero_ids != id]
        tmp <- dgensim[ pt_id1 == id | pt_id2 == id ] 
        tmp <- tmp[! pt_id1 %in% exclude_ids & ! pt_id2 %in% exclude_ids]
        setorder(tmp, -SIMILARITY)
        close_ids <- tmp[1:(args$cluster_size - 1), c(pt_id1, pt_id2)]
        close_ids <- close_ids[close_ids != id]

        tmp <- data.table(ID = c(id, close_ids), 
                          IDCLU = which(sero_ids == id), 
                          CLU_SIZE= length(c(id, close_ids))
        )
        dclusters <- rbind(dclusters, tmp)
}
# Save clusters.csv in the potential_network directory
filename <- file.path(args$out.dir, 'potential_network', 'clusters.rds')
saveRDS(dclusters, filename)

# To build: phscinput_runs_clusize_6_ncontrol0.rds
pty.runs <- copy(dclusters)
setnames(pty.runs, 'ID', 'UNIT_ID')
tmp <- phsc.input[, .(UNIT_ID, PANGEA_ID, RENAME_ID, SAMPLE_ID)]
pty.runs <- merge(pty.runs, tmp, by='UNIT_ID')
pty.runs[, `:=` (PTY_RUN=IDCLU, PTY_SIZE=CLU_SIZE)]
setcolorder(pty.runs, c( "UNIT_ID","IDCLU", "CLU_SIZE", "PTY_SIZE", "PTY_RUN", "PANGEA_ID", "RENAME_ID", "SAMPLE_ID"))

# Write processed samples
outfile <- paste0('phscinput_runs_clusize_', args$cluster_size,
                  '_ncontrol_',args$n_control, '.rds')
outfile <- file.path(args$out.dir, outfile)
cat('\nWriting to file ', outfile)
saveRDS(pty.runs, file=outfile)

# Set up next step in the analysis as well!
if(file.exists(args$controller))
{
        cmd <- paste0('cd ', dirname(args$controller), '\n')
        cmd <- paste0(cmd,
                      'qsub -v STEP="ali" ', args$controller
        )
        cat(system(cmd, intern= TRUE))
}
