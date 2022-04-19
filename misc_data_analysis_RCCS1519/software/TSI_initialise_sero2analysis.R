cat("\n\n===== TSI_initialise_sero2analysis.R =====\n\n")

require(data.table)
require(ggplot2)

usr <- Sys.info()[['user']]
if (usr == 'andrea')
{
        indir.deepsequence_analyses <- '~/Documents/Box/ratmann_deepseq_analyses/live'
        indir.deepsequencedata <- '~/Documents/Box/ratmann_pangea_deepsequencedata'
        tanya.rakai.dir <- '~/git/HIV-phyloTSI-main/RakExample_Tanya'
}else{
        indir.deepsequence_analyses <- '/rds/general/project/ratmann_deepseq_analyses/live'
        indir.deepsequencedata <- '/rds/general/project/ratmann_pangea_deepsequencedata'
        tanya.rakai.dir <- '~/git/HIV-phyloTSI/RakExample_Tanya'
}

args <- list( out.dir=file.path(indir.deepsequence_analyses, 'seroconverters2'))

###
# DATA
###

# Tanya's analysis
maf_tanya <- fread(file.path(tanya.rakai.dir, 'rak_maf.csv'), header=TRUE)
ps_tanya <- fread(file.path(tanya.rakai.dir, 'rak_ps.csv'), header=TRUE)
out_tanya <- fread(file.path(tanya.rakai.dir, 'rak_out.csv'), header=TRUE)

# Our previous results 
indir.deepsequence_analyses_old <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS1519_UVRI')
tmp <- file.path(indir.deepsequence_analyses_old, "19037_phsc_phscrelationships_seed_42_blacklist_report_TRUE_distance_threshold_1_min_reads_per_host_1_multinomial_TRUE_outgroup_name_BFR83HXB2_LAI_IIIB_BRUK03455_output_nexus_tree_TRUE_ratio_blacklist_threshold_0005_skip_summary_graph_TRUE")
dtsi.all.path <- list.files(tmp, pattern='tsi')
dtsi.all.path <- file.path(tmp, dtsi.all.path)
dtsi.1.path <- file.path(tmp, "aggregated_TSI_with_estimated_dates.csv")

# Meta data
file.db.sharing <- file.path(indir.deepsequencedata,"/PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv")
file.anonymisation.keys <- file.path(indir.deepsequence_analyses_old,'important_anonymisation_keys_210119.csv')
file.phsc.input.samples <- file.path(indir.deepsequence_analyses_old, '210120_RCCSUVRI_phscinput_samples.rds' )



###
# Microprocessing + helper function
###

aik <- fread(file.anonymisation.keys, header=TRUE)
aik[, V1:=NULL]

# note we are not considering the is_mrc variable
predictors.hiv.phylo.tsi <-  c("gag_lrtt", "pol_lrtt", "gp120_lrtt", "gag_maf3c", "gp41_maf3c", "gp41_maf12c", "gag_tips", "gp41_tips", "gp120_tips")

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


###
# Main functions
###

preprocess.tanya <- function()
{
        if(0)
        {
                t1 <- maf_tanya[,unique(sid)]
                t2 <- ps_tanya[,unique(sid)]
                t3 <- out_tanya[, unique(host.id)]
                stopifnot( all(sort(t1) == sort(t2)) & all(sort(t3) == sort(t2)) )
                rm(t1,t2,t3)
        }
        
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
        dtsi <- fread(dtsi.1.path)
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

        if(user=='andrea')
        {
                predictors.merged[, .N]
                require(gridExtra)
                do.call("grid.arrange", c(gg_list, ncol=3)) -> gg
                ggsave(file="220408_inputs_outputs_comparison.png", gg, w=10, h=14)
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

        predictors.ours[, range(gp120_tips)]
        predictors.tanya[, range(gp120_tips)]
        colnames(predictors.ours)

}


find.common.seroconverters.pangeaids <- function(predictors.ours, predictors.tanya)
{

        predictors.ours

        # get Tanya dates
        ddates.tanya <- predictors.tanya[, .(pt_id=paste0('RK-', PT_ID), sample_date)]

        # get our sample collection dates from Kate's db sharing
        phsc.input.samples <- readRDS(file.phsc.input.samples)
        phsc.input.samples[, PANGEA_ID:=gsub('^.*?_', '', PANGEA_ID)]
        phsc.input.samples[, SAMPLE_ID:=NULL]
        phsc.input.samples

        ddates.ours <- fread(file.db.sharing, select=c('pt_id', 'pangea_id', 'visit_dt'))
        ddates.ours <- unique(ddates.ours)

        # Merge
        ddates <- merge(ddates.ours, ddates.tanya, by='pt_id')

        # check that we can always find matching Pangea id
        cols <- c('visit_dt','sample_date')
        ddates[, (cols):=lapply(.SD, as.Date, format='%Y-%m-%d'), .SDcols=cols]
        stopifnot(ddates[, any(abs(visit_dt - sample_date)<=10) ,by='pt_id'][, all(V1)])
        ddates <- ddates[abs(visit_dt - sample_date)<=10, ,] 
        
        # There are 12 pangea_ids sharing same pt_id and visit_dt...
        if(0) 
        {
                ddates[, any(.N > 1),by='pt_id'][V1==TRUE, pt_id] -> tmp
                tmp <- ddates[pt_id %in%tmp]
                tmp[, sample_date:=NULL]
                fwrite(tmp, file='220411_pangea_ids_with_same_dates.csv')

                ddates[, table(gsub('-.*?$', '', pangea_id))]
                # From the original source, I would assume not.
                # I still think it could make sense to keep them.
                # In particular, I would only keep the PG15 as it seems that those were used in the analysis.
                ddb <- fread(file.db.sharing, drop='shiver_consensus')
                ddb[pangea_id %in% tmp$pangea_id]
                
        }
        
        # Only pt_ids with 2 such pangea_ids have a pgid starting with PG19, I ll remove those
        cat('Removing PANGEA_IDs starting with PG19 \n')
        ddates <- ddates[!grepl('^PG19', pangea_id)] 

        tmp_ids <- ddates[sample_date != visit_dt, pt_id ]
        tmp1 <- phsc.input.samples[UNIT_ID %in% tmp_ids]
        tmp1[, .N, by='UNIT_ID']
        cat('There are', length(tmp_ids),'out of', ddates[, .N], 'individuals with disagreeing dates\n')
        cat('Those inds do not appear in HIV or ALLHIV...\n')

        # ONCE WE FIND INDIVIDUALS WITH MATCHING DATES ETC... 
        # WHAT DO I NEED TO DO TO RUN THE NEW ANALYSES???
        ddates[! pt_id %in% tmp_ids]
}


make_new_input_samples <- function(ddates)
{
        phsc_input_samples_path <- file.path(indir.deepsequence_analyses_old,"210120_RCCSUVRI_phscinput_samples.rds")
        phsc_input_samples <- readRDS(phsc_input_samples_path)
        phsc_input_samples[, pangea_id := gsub('^.*?_', '', PANGEA_ID)]

        phsc_input_new <- phsc_input_samples[pangea_id %in% ddates$pangea_id]
        phsc_input_new[, pangea_id :=NULL]
        
        msg <- 'all pangeaIDs end by fq1\n'
        if(phsc_input_new[, all(gsub('^.*?-','',RENAME_ID)=='fq1')])    cat(msg)
        phsc_input_new
}

write.hpc.input <- function(ddates, out.dir)
{
        
        # Stores on the HPC everything is needed to initialise the analysis.

        # make phscinput_samples.rds
        # ______________________________
        phsc_samples <- make_new_input_samples(ddates)
        filename=file.path(out.dir, paste0('220419_phscinput_samples.rds'))
        saveRDS(phsc_samples, filename)

        # Make clusters.rds
        # ______________________________
        tmp <- ddates[, pt_id]
        dclus <- data.table(ID=tmp, IDCLU=1:3)
        dclus[,CLUSIZE:=.N,by='IDCLU']
        setkey(dclus, IDCLU)
        filename=file.path(out.dir, 'potential_network', 'clusters.rds')
        saveRDS(dclus, filename)


        # make phscinput_runs_clusize...
        # ______________________________
        suffix <- phsc_samples[, .(UNIT_ID,PANGEA_ID, RENAME_ID, SAMPLE_ID)]
        dclus[, `:=` (PTY_RUN=IDCLU,  PTY_SIZE=CLUSIZE) ]
        merge(dclus, suffix, by.x='ID', by.y='UNIT_ID')
        filename=file.path(out.dir, paste0('phscinput_runs_clusize_',max(dclus(CLUSIZE)),'_ncontrol_0.rds'))
        saveRDS(dclus, filename)
}


########################
# MAIN
########################
dir.create(args$out.dir)
dir.create(file.path(args$out.dir, 'potential_network'))

predictors.ours <- preprocess.ours()
predictors.tanya <- preprocess.tanya() 
ddates <- find.common.seroconverters.pangeaids(predictors.ours, predictors.tanya)
write.hpc.input(ddates, args$out.dir)
