cat('\n\n=====  TSI_estimate_dates.R =====\n\n')

# Objective of this script is to load every TSI prediction produced in the previous step
# And then to attribute estimates of dates of infectios.
# To this end, we require a dataframe containing the dates of collection of each sample.
library(data.table)
library(lubridate)

option_list <- list(
  optparse::make_option(
    "--out_dir_base",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where all output is stored [default]",
    dest = 'out.dir'
  ),
  optparse::make_option(
    "--relationship_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to directory containing analyse_trees.R results", 
    dest = 'rel.dir'
  ),
  optparse::make_option(
    "--pkg_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to package directory, used as long we don t build an R package [default]",
    dest = 'pkg.dir'
  ),
  optparse::make_option(
    "--input_samples",
    type = "character",
    default = NA_character_, 
# '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_phscinput_samples.rds',
    help = "Absolute file path to rds containing PANGEA_IDs and RENAME_IDs", 
    dest = 'phsc.samples'
  ),
#  optparse::make_option(
#    "--coll_dates",
#    type = "character",
#    default = NA_character_,
#    help = "Absolute file path to dataset containing times of collection", 
#    dest = 'coll.dates'
#  ),
  optparse::make_option(
    "--date",
    type = 'character',
    default = NA,
    metavar = '"YYYY-MM-DD"',
    help = 'As of date to extract data from.  Defaults to today.',
    dest = 'date'
  ),
  optparse::make_option(
    "--controller",
    type = "character",
    default = NA_character_, # Think about adding the controller in the software directory
    help = "Path to sh script directing the full analysis",
    dest = 'controller'
  )

)


args <-  optparse::parse_args(optparse::OptionParser(option_list = option_list))


################
# Helpers
################

.get.dates <- function(files)
{
        .f <- function(file)
        {
                ddates <- fread(file)
                ddates <- unique(ddates[, .(pangea_id, visit_dt)])
                ddates[, as.Date(visit_dt, format="%Y-%m-%d")]
                ddates
        }
        ddates <- lapply(files, .f)
        ddates <- rbindlist(ddates)
        stopifnot(ddates[, anyDuplicated(pangea_id) == 0,])
        ddates
}


################
# Testing
################

user <- Sys.info()[['user']]
if(user == 'andrea')
{
        args <- list(
                     out.dir='~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS_MRC_UVRI_TSI/',
                     rel.dir= "~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS_MRC_UVRI_TSI/2022_08_22_phsc_phscTSI_sd_42_sdt_002_005_dsl_100_mr_30_mlt_T_npb_T_og_REF_BFR83HXB2_LAI_IIIB_BRU_K03455_phcb_T_rtt_00/",
                     phsc.samples="~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS_MRC_UVRI_TSI/220331_RCCSUVRI_phscinput_samples_with_bf.rds",
                     date='2022-08-22'
        )
}

################
# MAIN
################

args$date <- gsub('-','_',args$date)
if( ! grepl('output$', args$out.dir))
{
        args$out.dir <- file.path(args$out.dir, paste0(args$date, "_phsc_output"))
}

stopifnot(dir.exists(args$rel.dir))
stopifnot(dir.exists(args$out.dir))

# Load PANGEA_ID for each sample and collection dates
dsamples <- setDT(readRDS(args$phsc.samples))
dsamples <- unique(dsamples[, .(PANGEA_ID, RENAME_ID)])


file.basefreqs <- list.files(dirname(args$out.dir), 'samples_with_bf.rds', full.names=T)
file.basefreqs_old <- list.files(args$rel.dir, pattern='basefreqs_used.csv$', full.names=TRUE)

if(length(file.basefreqs))
{
        dfiles <- list.files(args$rel.dir, pattern='_tsi.csv$', full.names = TRUE)
        dpreds <- lapply(dfiles, fread)
        names(dpreds) <- .gs(dfiles)
        cols <- grep('host.id|^RF', colnames(dpreds[[1]]), value=T)
        dpreds <- lapply(dpreds, function(DT) subset(DT, select=cols) )
        dpreds <- rbindlist(dpreds, idcol = 'PTY')
        dpreds[, AID := gsub('-fq[0-9]$','',host.id)]
        
        setnames(dpreds, 'host.id', 'RENAME_ID')
        if(dpreds[1, RENAME_ID == AID])
                dpreds[, RENAME_ID := NULL]

}else if( length(file.basefreqs_old) ) 
{
        # Load results from prev
        tmp <- list.files(args$rel.dir, pattern='_tsi.csv$', full.names = TRUE)
        tmp <- data.table(tsi.path=tmp, 
                          pty=gsub('^.*?ptyr|_tsi.csv', '', tmp))
        tmp1 <- data.table(bfs.paths=tmp1, 
                          pty=gsub('^.*?ptyr|_basefreqs_used.csv', '', tmp1))
        dfiles <- merge(tmp, tmp1, by='pty')
        dfiles[, pty:=as.integer(pty)]
        setkey(dfiles, pty)

        dbfs <- dfiles[,{
                bfs <- read.csv(bfs.paths)[, 2]
                list(RENAME_ID=bfs[which(bfs != "pos")])
        }, by=pty]
        dbfs <- merge(dbfs, dsamples, by='RENAME_ID')
        dbfs[, AID:=gsub('-fq[0-9]$', '', RENAME_ID)]
        stopifnot(dbfs[, .N == 1 , by=c('pty', 'AID')][, all(V1)])

        dpreds <- dfiles[,{
                tsi <- fread(tsi.path)
                cols <- grep('^RF', colnames(tsi), value=T)
                tsi[, c('host.id', cols)]
        }, by=pty]
        dpreds <- merge(dbfs[, .(pty, AID, RENAME_ID)],
              dpreds,
              by.y=c('pty', 'host.id'),
              by.x=c('pty', 'AID'))
}


# Check if results for same sequence are consistent among ptys.
# Are median prediction always in the interestion of all CrInt?
if(0)
{
        tmp <- dpreds[, .N , by='RENAME_ID'][N>1, RENAME_ID]
        tmp <- dpreds[ RENAME_ID %in% tmp, ]
        setkey(tmp, RENAME_ID)
        tmp[, `:=` (MAX=min(RF_pred_max_linear),
                    MIN=max(RF_pred_min_linear)), by='RENAME_ID' ]
        tmp1 <- tmp[,RF_pred_linear <= MAX && RF_pred_linear >= MIN,by='RENAME_ID']
        # in my first analysis, 7% fell outside the intersections of the CIs
        # can check if correlates with RF_pred_MAE
}

# Take median of predictions and cc's in sqrt space, then transform to linear space
# Alternatively could pick prediciton with minimum MAE
if( dpreds[, .N,by='RENAME_ID'][, any(N>1)] ) 
{
        cols <- grep('RF',colnames(dpreds), value=T)
        cols <- grep('linear', cols, value=T, invert = TRUE)
        dpreds <- dpreds[, lapply(.SD, median) ,by='RENAME_ID', .SDcols=cols]
        cols = grep('pred_sqrt|cc',cols, value=T)
        cols1 <- gsub('sqrt', 'linear',cols)
        cols1 <- gsub('cc025', 'pred_min_linear',cols1)
        cols1 <- gsub('cc975', 'pred_max_linear',cols1)
        dpreds [, (cols1):=lapply(.SD, function(x) x^2), .SDcols=cols, by='RENAME_ID']
}


# if dataset with dates of collection exists, also compute date of infection estimates!
if(user != 'andrea')
{
        db.sharing.path.rccs <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv'
        db.sharing.path.mrc <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv'
}else{
        db.sharing.path.rccs <- '/home/andrea/Documents/Box/ratmann_pangea_deepsequencedata/PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv'
        db.sharing.path.mrc <-  '/home/andrea/Documents/Box/ratmann_pangea_deepsequencedata/PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv' 
}

sampling.dates.available <- all(file.exists(c(db.sharing.path.mrc, db.sharing.path.rccs)))
if(sampling.dates.available)
{

        files <- c(db.sharing.path.rccs, db.sharing.path.mrc)
        ddates <- .get.dates(files)

        tmp <- dsamples[, .(RENAME_ID, pangea_id)]
        ddates <- merge(tmp, ddates, by='pangea_id', all.x=T)

        dpreds <- merge(dpreds, ddates[, .(RENAME_ID, visit_dt)], by='RENAME_ID', all.x=TRUE)
        cat(dpreds[, sum(is.na(visit_dt))], 'base frequency files have no associated sample collection date\n')

        # Transform TSI estimates in day and then estimate DateofInfection
        cols <- grep('linear', colnames(dpreds), value=T)
        cols1 <- gsub('linear', 'days', cols)
        tmp <- dpreds[, lapply(.SD, function(x) as.integer(x*365)) , .SDcols=cols, by=c('RENAME_ID', 'visit_dt')]
        setnames(tmp, cols, cols1)
        tmp <- tmp[, lapply(.SD, function(x){visit_dt - x}) , .SDcols=cols1, by='RENAME_ID']
        setnames(tmp, grep('pred_min', names(tmp),value=T), 'pred_doi_max')
        setnames(tmp, grep('pred_max', names(tmp),value=T), 'pred_doi_min')
        setnames(tmp, grep('RF', names(tmp),value=T), 'pred_doi_mid')

        dpreds <- merge(dpreds, tmp, by='RENAME_ID')
}

filename <- ifelse(sampling.dates.available,
               file.path(args$rel.dir,'aggregated_TSI_with_estimated_dates.csv'),
               file.path(args$rel.dir,'aggregated_TSI.csv')
)
cat('Saving aggregated estimates to:\n', filename, '\n')
fwrite(dpreds, filename, quote=FALSE)


# Queue next step
#________________

source(file.path(args$pkg.dir, "utility.R"))
qsub.next.step(file=args$controller,
               next_step='pst')

cat('\n\n=====  end of script =====\n\n')
