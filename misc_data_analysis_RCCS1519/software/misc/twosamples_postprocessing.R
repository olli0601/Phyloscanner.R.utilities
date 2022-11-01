# AIMS:
# Once having obtained the output from the two sample analysis,
# we want to compare wheter the estimates obtained from those samples are comparable.
# Also, we would like to know which sample to select to obtain maximal prediction accuracy.
# A priori, the best one should be the first positive sample ever, because this suffers of a lower degree of "left" censoring.

################
# DEPENDENCIES #
################

library(data.table)
library(ggplot2)

################
#    PATHS     #
################

usr <- Sys.info()[['user']]
if(usr == 'andrea')
{
        indir.deepdata <- '/home/andrea/Documents/Box/ratmann_pangea_deepsequencedata/'
        indir.deepanalyses <- '/home/andrea/Documents/Box/ratmann_deepseq_analyses/live/'

}else{
        indir.deepdata <- '/home/andrea/Documents/Box/ratmann_pangea_deepsequencedata/'
        indir.deepanalyses <- '/home/andrea/Documents/Box/ratmann_deepseq_analyses/live/'
}


indir.twosamples <- file.path(indir.deepanalyses, 'twosamples', "2022_09_02_phsc_phscTSI_sd_42_sdt_002_005_dsl_100_mr_30_mlt_T_npb_T_og_REF_BFR83HXB2_LAI_IIIB_BRU_K03455_phcb_T_rtt_001_rla_T_zla_T" )

path.tsi <- list.files(indir.twosamples, pattern='tsi.csv$', full.names = TRUE)
path.input.samples <- file.path(indir.deepanalyses, 'twosamples', "220331_RCCSUVRI_phscinput_samples_with_bf.rds")


################
#    HELPERS   #
################


################
#     MAIN     #
################

dsamples <- readRDS(path.input.samples)
dsamples <- as.data.table(dsamples)
ddates <- dsamples[, .(host.id=RENAME_ID, visit_dt)]

dtsi <- lapply(path.tsi, fread)

dtsi <- rbindlist(dtsi, idcol = TRUE)

cols <- grep('id$|^RF',names(dtsi), value=TRUE)
dcomp <- dtsi[, ..cols]

# define person id
dcomp[, study_id := gsub('-fq[0-9]$', '', host.id)]

# merge ddates
dcomp <- merge(dcomp, ddates)
dcomp[, `:=` (
            pred_doi_mid = visit_dt - as.integer(RF_pred_linear*365),
            pred_doi_max = visit_dt - as.integer(RF_pred_min_linear*365),
            pred_doi_min = visit_dt - as.integer(RF_pred_max_linear*365)
)]

setkey(dcomp, study_id, visit_dt)
dcomp[, visit := fifelse(.id <= 2, 'first', 'second')]

dcomp[, {
        z1 <- which(visit=='first')
        z2 <- which(visit=='second')
        list(visit_dt[z1] <= visit_dt[z2])
        },by='study_id'][, all(V1)] |> stopifnot()


cols <- grep('RF|study_id|pred_doi|visit_dt', names(dcomp), value=TRUE)

tmp <- merge(
        dcomp[visit=='first', ..cols],
        dcomp[visit=='second', ..cols],
        by='study_id')

names(tmp) <- gsub('x$', '1', names(tmp))
names(tmp) <- gsub('y$', '2', names(tmp))
tmp[, deltaT := as.integer(visit_dt.2 - visit_dt.1)]

# We expect to underestimate TSI when it increases.
# This means that the TSI for later samples will tend to be smaller than they should be on average
# As a consequence, we expect the estimated date of infection to be later for later samples
# This means most points should be above the diagonal here

.d <- function(x,y) mean(x-y)^2
tmp[, cor(RF_pred_linear.1, RF_pred_linear.2)]
tmp[, cor(RF_pred_linear.1, RF_pred_linear.2 - deltaT)]
tmp[, .d(RF_pred_linear.1, RF_pred_linear.2)]
# tmp[, .d(RF_pred_linear.1, RF_pred_linear.2 - deltaT)]
# tmp[, .d(RF_pred_linear.1, RF_pred_linear.2 + deltaT)]


ggplot(tmp, aes(x=RF_pred_linear.1, y=RF_pred_linear.2)) +
        geom_point() + 
        geom_abline(slope=1, linetype='dashed', color='red') +
        stat_cor() +
        theme_bw() 

ggplot(tmp, aes(x=RF_pred_linear.1, y=RF_pred_linear.2 - deltaT)) +
        geom_point() + 
        geom_abline(slope=1, linetype='dashed', color='red') +
        stat_cor() +
        theme_bw() 

p1 <- ggplot(tmp, aes(x=pred_doi_mid.1, y=pred_doi_mid.2)) +
        geom_point() + 
        geom_abline(slope=1, linetype='dashed', color='red') +
        stat_cor() +
        theme_bw() +
        labs(x='Estimated date of infection with first sample', 
             y='Estimated date of infection with second sample',
             title='No "adjustment" (ie: as it should be)')

.f <- function(x) as.integer(x)
tmp[, cor(.f(pred_doi_mid.1), .f(pred_doi_mid.2))]
tmp[, cor(.f(pred_doi_mid.1), .f(pred_doi_mid.2))]

p2 <- ggplot(tmp, aes(x=pred_doi_mid.1, y=pred_doi_mid.2 - deltaT)) +
        geom_point() + 
        geom_abline(slope=1, linetype='dashed', color='red') + 
        stat_cor() +
        theme_bw() +
        labs(x='Estimated date of infection with first sample', 
             y='Estimated date of infection with second sample',
             title='With DeltaT "adjustment"')

p <- ggarrange(p1, p2, nrow=1, ncol=2)

filename <- file.path(indir.twosamples, 'twosamples_scatter_doi_first_vs_second_with_without_adj.png')
ggsave(p, file=filename, w=12, h=7)




