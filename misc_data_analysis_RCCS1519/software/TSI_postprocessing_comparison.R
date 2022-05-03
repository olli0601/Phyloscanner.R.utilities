# We need to compare the predictions obtained with Tanya's run and with ours
# As such, I want to be able to produce two plots systematically:
# - the first one comparing the predictors and predictions between methods
# - the second one testing the predictions on the known seroconversion intervals

# I need to load the aggregated predictions for our analysis
# And Tanyas


# Packages
#_________

require(ggplot2)
require(data.table)
require(gridExtra)


# Args
#_____
cat("arguments...\n")


option_list <- list(
  optparse::make_option(
    "--relationship_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to directory containing analyse_trees.R results", 
    dest = 'rel.dir'
  ),
  optparse::make_option(
    "--TSI_dir",
    type = "character",
    default = '~/git/HIV-phyloTSI-main/',
    help = "Absolute file path to HIV-phylo-TSI-main repository", 
    dest = 'TSI.dir'
  ),
  optparse::make_option(
    "--input_samples",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to input samples rds containing PANGEA_IDs and RENAME_IDs", 
    dest = 'phsc.samples'
  )
)

args <-
  optparse::parse_args(optparse::OptionParser(option_list = option_list))

args[['help']] <- NULL

print(args)

tmp <- file.exists(unlist(args))
stopifnot(all(tmp))

# Paths
#________ 
cat("...loading paths\n")

usr <- Sys.info()[['user']]
if (usr == 'andrea')
{
        indir.deepsequence_analyses <- '~/Documents/Box/ratmann_deepseq_analyses/live'
        indir.deepsequencedata <- '~/Documents/Box/ratmann_pangea_deepsequencedata'
        tanya.rakai.dir <- '~/git/HIV-phyloTSI-main/RakExample_Tanya'

        args$rel.dir <-  file.path(indir.deepsequence_analyses, 'seroconverters2/phscrel')
        tmp <- dirname(args$rel.dir)
        args$phsc.samples <- list.files(tmp, pattern='phscinput_samples', full.names=T)
}else{
        indir.deepsequence_analyses <- '/rds/general/project/ratmann_deepseq_analyses/live'
        indir.deepsequencedata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live'
        tanya.rakai.dir <- file.path(args$TSI.dir, 'RakExample_Tanya')
}

# Analysis results:
dtsi.all.path <- list.files(args$rel.dir, pattern='tsi', full.names=TRUE)
dtsi.aggregated.path <- file.path(args$rel.dir, "aggregated_TSI_with_estimated_dates.csv")

# Meta data
indir.deepsequence_analyses_old <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS1519_UVRI')
file.db.sharing <- file.path(indir.deepsequencedata,"/PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv")
file.anonymisation.keys <- file.path(indir.deepsequence_analyses_old,'important_anonymisation_keys_210119.csv')
file.seroconverters <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', '220329_TSI_seroconverters.csv')

tmp <-file.exists(indir.deepsequence_analyses_old, file.db.sharing,
                  file.anonymisation.keys, file.seroconverters)
)
print(tmp)
stopifnot(all(tmp))






# Microproccesing + helper functions
#___________________________________
cat("helper functions...\n")

predictors.hiv.phylo.tsi <-  c("gag_lrtt", "pol_lrtt", "gp120_lrtt", "gag_maf3c", "gp41_maf3c", "gp41_maf12c", "gag_tips", "gp41_tips", "gp120_tips")

aik <- fread(file.anonymisation.keys, header=TRUE)
aik[, V1:=NULL]

.aid2studyid <- function(x, rm_prefix=FALSE, with_fq=FALSE)
{
        x <- data.table(RENAME_ID = x)
        
        cols <- c('RENAME_ID', 'UNIT_ID')
        tmp <- unique(dsamples[, ..cols])
        if(! with_fq)
        {
                tmp[, RENAME_ID:=gsub('-fq[0-9]$', '', RENAME_ID)]
        }

        x <- merge(x, tmp, by='RENAME_ID', all.x=TRUE)$UNIT_ID
        if(rm_prefix)
        {
                x <- gsub('^.*?-','',x)
        }
        x
}

.sid2pid.and.date <- function(x)
{
        y1 <- gsub('^(.*?)-(.*?)$', '\\1', x)
        y2 <- gsub('^(.*?)-(.*?)$', '\\2', x)

        y2 <- as.numeric(y2)
        year <- y2 %/% 1 
        date <- as.Date(paste0('01-01-', year), '%d-%m-%Y')
        tmp1 <- round(y2 %% 1 *365)
        date <- date + tmp1

        list(PID=y1,date=date)
}

preprocess.tanya <- function()
{
        # Tanya's analysis
        maf_tanya <- fread(file.path(tanya.rakai.dir, 'rak_maf.csv'), header=TRUE)
        ps_tanya <- fread(file.path(tanya.rakai.dir, 'rak_ps.csv'), header=TRUE)
        out_tanya <- fread(file.path(tanya.rakai.dir, 'rak_out.csv'), header=TRUE)

        # extract PID and sample date from host_id
        newcols <- c('PT_ID', 'sample_date')
        maf_tanya[, (newcols):=.sid2pid.and.date(sid)]
        ps_tanya[, (newcols):=.sid2pid.and.date(sid)]
        out_tanya[, (newcols):=.sid2pid.and.date(host.id)]
        out_tanya
        
        # extract predictors and predictions
        tmp <- grep('sqrt|cc', colnames(out_tanya), value=TRUE)
        cols <- c("PT_ID", "sample_date", predictors.hiv.phylo.tsi, tmp)
        predictors.tanya <- out_tanya[, .SD, .SDcols=cols]
        predictors.tanya
}

preprocess.ours <- function()
{
        # Load all outputs from HIVphyloTSI
        dtsi <- fread(dtsi.aggregated.path)
        dtsi[, V1 := NULL]
        tmp <- lapply(dtsi.all.path, fread)
        dtsi.all <- rbindlist(tmp, use.names=TRUE)
        dtsi.all <- merge(aik, dtsi.all, by.y='host.id',by.x='AID')

        if(0)
        {
                # It seems like lrtt vary for different runs of the same ID
                # On the other hand the mafs dont. This makes sense
                dtsi.all.counts <- dtsi.all[, .(N=.N, 
                                                Nmaf=uniqueN(genome_maf3c),
                                                Nlrtt=uniqueN(genome_lrtt)), by='host.id']
        }

        # extract predictors and predictions
        tmp <- grep('sqrt|cc', colnames(dtsi.all), value=TRUE)
        cols <- c("AID", "PT_ID", predictors.hiv.phylo.tsi, tmp)
        predictors.ours <- dtsi.all[, .SD, .SDcols=cols]
        predictors.ours <- predictors.ours[grepl('^RK',PT_ID)]
        predictors.ours[, PT_ID := gsub('RK-', '',PT_ID)] 
        predictors.ours
}

compare.predictors <- function(predictors.ours, predictors.tanya)
{

        # merge data and sort by differences in preds
        cols <- colnames(predictors.ours)
        cols <- cols[cols %in% colnames(predictors.tanya)]
        predictors.merged <- merge(predictors.ours[,..cols], predictors.tanya[,..cols], by='PT_ID')
        predictors.merged[, error := (RF_pred_sqrt.x - RF_pred_sqrt.y)]

        # plot comparisons
        .plot.pred.comparison <- function(tx)
        {
                ty <- gsub('x$', 'y', tx)
                t <- c(tx,ty, 'error')

                dplot <- predictors.merged[, .SD, .SDcols=t]
                colnames(dplot) <- c('x','y','c')
                gg <- ggplot(dplot, aes(x=x, y=y, color=c)) +
                        geom_abline(slope=1, color='black') + 
                        geom_smooth(method='lm', formula= y~x, color="grey50", linetype='dotted', se=FALSE)+ 
                        geom_point() +
                        guides(colour='none') +
                        # scale_colour_gradient(low="blue", high="red", midpoint=0)
                        scale_colour_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 0) +
                        labs(x="our", y="tanya's",                             title=gsub('.x$','',tx)) +
                        theme_bw() +
                        expand_limits(x = 0, y = 0) +
                        theme(plot.title = element_text(hjust = 0.5))
                gg
        }
        cols <- colnames(predictors.merged)
        tmp <- grep('\\.x$',cols, value=T)
        gg_list <- lapply(tmp, .plot.pred.comparison)

        if(1)
        {
                predictors.merged[, .N]
                do.call("grid.arrange", c(gg_list, ncol=3)) -> gg1
                filename=file.path(args$rel.dir , "inputs_outputs_comparison.png")
                ggsave(filename , gg, w=10, h=14)
        }

        
        .plot.compare.95ci <- function(predictors.merged)
        {
                grep('RF',colnames(predictors.merged), value=T) -> cols
                tmp <- predictors.merged[, .SD, .SDcols=cols]
                gg <- ggplot(tmp, aes(x=RF_pred_sqrt.x, y=RF_pred_sqrt.y)) + 
                        geom_abline(slope=1, linetype='dotted', color='red') + 
                        geom_point() + 
                        geom_crossbar(aes(ymin=RF_cc025.y, ymax=RF_cc975.y)) + 
                        geom_linerange(aes(xmin=RF_cc025.x, xmax=RF_cc975.x))
                gg
        }
        gg2 <- .plot.compare.95ci(predictors.merged)
        
        return(list(gg1, gg2))
}


plot.sero.predictions <- function(data, title='')
{       
        data[, in_sero_interval := (RF_pred_linear < knownTSI_max) & (RF_pred_linear > knownTSI_min)]
        data <- data[!is.na(knownTSI) & !is.na(RF_pred_linear)]
        lab <- data[knownTSI < 5, cor(knownTSI, RF_pred_linear)]
        lab <- round(lab, 2)
        lab <- paste0('corr = ', lab, '\n')
        tmp <- data[, round(mean(in_sero_interval), 2)]
        lab <- paste0(lab, 'coherent predictions = ', tmp)
        y <- data[, 0.9* max(knownTSI-RF_pred_linear)]

        cat(data[,.N], ' datapoints plotted \n' )
        ggplot(data, aes(x=knownTSI, y=knownTSI - RF_pred_linear)) + 
                geom_point(aes(color=in_sero_interval)) +
                geom_hline(aes(yintercept=0), color='red', linetype='dotted') + 
                theme_bw() +
                theme(legend.position='bottom') +
                geom_label(aes(x=3, y=y, label=lab)) + 
                labs(x="Known TSI (sample collection date - midpoint)", y="Known TSI - estimated TSI", color='Prediction in seroconversion interval') +
                ggtitle(title)

}

show.corrs <- function(data, title='')
{      
        data <- data[!is.na(knownTSI) & !is.na(RF_pred_linear)]
        lab <- round( cor(data$knownTSI, data$RF_pred_linear), 2 )
        lab <- paste0('corr =', lab)
        tmp <- data[, 0.9* max(RF_pred_linear)]

        cat(data[,.N], ' datapoints plotted \n' )
        ggplot(data, aes(x=knownTSI, y= RF_pred_linear)) + 
                geom_point() + 
                geom_hline(aes(yintercept=0), color='red', linetype='dotted') + 
                geom_label(aes(x=2.5, y=tmp, label=lab)) + 
                theme_bw() + 
                theme(legend.position='bottom')+
                labs(x="Known TSI" , y="Estimated TSI") +
                ggtitle(title)
}

# Main
#_____
cat(" main...\n")

predictors.tanya <- preprocess.tanya() 
predictors.ours <- preprocess.ours()

# Save 1st plot of interest
x <- compare.predictors(predictors.ours, predictors.tanya)

# 2nd plot of interest...

dsamples <- readRDS(args$phsc.samples)
dsamples[, pangea_id:= gsub("^.*?_", "",PANGEA_ID)]

dsero <- fread(file.seroconverters)
dsero <- dsero[pangea_id %in% dsamples$pangea_id]
cols <- c('date_last_negative', 'date_first_positive')
new_cols <- c('knownTSI_max', 'knownTSI_min')
dsero[, (new_cols) := lapply(.SD , function(x) (visit_dt - x)/365), .SDcols=cols ]
dsero <- dsero[, .(study_id, visit_dt, midpoint, knownTSI, knownTSI_min, knownTSI_max)]

dtsi <- fread(dtsi.aggregated.path, drop='V1')
dtsi[, study_id:=.aid2studyid(RENAME_ID, with_fq=TRUE),]

# TODO: some of the inds which Tanya treated as seroconverters 
# are not seroconverters in our eyes:
# if(0)
# {
#         ?write.csv
#         dtsi[, study_id[! study_id %in% dsero$study_id] ] |>
#                 write.table('220428_tanyaseroc_not_ourserconv.csv', 
#                             row.names=F, sep=',', col.names=F)
# }


dtsi <- merge(dtsi, dsero, by='study_id')
dtsi <- unique(dtsi)

ttl=gsub('^.*?live/','',args$rel.dir)
ttl=gsub('/',' ',ttl)
ttl=gsub('_phsc_.*?$',' ',ttl)

g <- plot.sero.predictions(dtsi, title=ttl)
filename <- file.path(args$rel.dir, 'sero_predictions_vs_known.png')
ggsave(filename,g, width=10, height=8)
