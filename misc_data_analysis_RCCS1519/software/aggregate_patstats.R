cat('\n\n=====  aggregate_patstats.R =====\n\n')

# This script simply aggregates the patstats files into a unique file
# TODO: could add some checks, or at least print statements.
# TODO: do I need to write a sh script for this?

library(data.table)

option_list <- list(
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
    "--controller",
    type = "character",
    default = NA_character_, # Think about adding the controller in the software directory
    help = "Path to sh script directing the full analysis",
    dest = 'controller'
  )
)

args <-
  optparse::parse_args(optparse::OptionParser(option_list = option_list))

#
# test
#
if (0) {
  args <- list(
        prj.dir=NA,
        out.dir="~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/19037_phsc_phscrelationships_seed_42_blacklist_report_TRUE_distance_threshold_1_min_reads_per_host_1_multinomial_TRUE_outgroup_name_BFR83HXB2_LAI_IIIB_BRUK03455_output_nexus_tree_TRUE_ratio_blacklist_threshold_0005_skip_summary_graph_TRUE",
        controller=NA
  )
}

# args$out.dir <- '~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/19037_phsc_phscrelationships_seed_42_blacklist_report_TRUE_distance_threshold_1_min_reads_per_host_1_multinomial_TRUE_outgroup_name_BFR83HXB2_LAI_IIIB_BRUK03455_output_nexus_tree_TRUE_ratio_blacklist_threshold_0005_skip_summary_graph_TRUE'

zip.files <- list.files(args$out.dir, pattern='.zip$', full.names=TRUE)
patstats <- lapply(zip.files, function(x){
        csv.name <- unzip(x, list = TRUE)$Name
        csv.name <- grep('_patStats.csv$',csv.name,value = T)
        patstat <- data.table(read.csv(unz(x, csv.name), header = TRUE, sep = ","))
})
zip.files <- gsub('^.*?ptyr','',zip.files)
zip.files <- gsub('_otherstuff.zip','',zip.files)
names(patstats) <- zip.files
ans <- rbindlist(patstats,use.names = T, idcol = 'PTY_RUN',fill=T)
gc()

write.csv(ans, file = file.path(args$out.dir,'patstats.csv'))
#


