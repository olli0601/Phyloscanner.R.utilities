# The phyloscanner run - find transmission networks 

# Preamble This script aims to find transmission networks and the most likely transmission chains using the phyloscanner outputs.

# Load the required packages
library(data.table)
library(tidyverse)
library(dplyr)
library(glue)
library(igraph)
library(RBGL)
library(phyloscannerR)

#
# Define input arguments that can be changed by users
#
option_list <- list(
  optparse::make_option(
    c("-v", "--verbose"),
    action = "store_true",
    default = FALSE,
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
    "--save_data",
    action = "store_true",
    default = TRUE,
    help = "Save data [default]",
    dest = 'if_save_data'
  ),
  optparse::make_option(
    "--env_name",
    type = "character",
    default = 'phyloenv',
    help = "Conda environment name [default]",
    dest = 'env_name'
  ),
  optparse::make_option(
    "--classification_rule",
    type = "character",
    default = 'o',
    help = "Rules for classifying linked and directed pairs. It takes values o or m. 
    o: a pair is linked if the linkage score > 0.6, and directed if the ancestral / (ancestral + descedant) > 0.6.
    m: a pair is linked if the linkage score > 0.5, and directed if the ancestral / all > 0.33. 
    b: both [default]",
    dest = 'classif_rule'
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
    help = "Absolute file path to base directory where all the tree outputs are stored [default]",
    dest = 'out.dir'
  ),
  optparse::make_option(
    "--phylo_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where the phyloscanner outputs are stored [default]",
    dest = 'phylo.dir'
  ),
  optparse::make_option(
    "--meta_data",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to meta.data containing serohistory + demographic information [default]",
    dest = 'meta.data'
  ),
  optparse::make_option(
    "--date",
    type = 'character',
    default = '2022-02-04',
    metavar = '"YYYY-MM-DD"',
    help = 'As of date to extract data from.  Defaults to today.',
    dest = 'date'
  )
)

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

#
# Helpers
#

indir.deepsequencedata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live'
indir.deepanalyses.xiaoyue <- '/home/andrea/HPC/project/ratmann_xiaoyue_jrssc2022_analyses/live'

get.sample.collection.dates <- function(select_aid=NULL, get_first_visit=FALSE)
{
    # get collection dates 
    path.sdates.rccs <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS','200316_pangea_db_sharing_extract_rakai.csv')
    path.sdates.mrc <- file.path(indir.deepsequencedata, 'PANGEA2_MRC','200319_pangea_db_sharing_extract_mrc.csv')

    files <- c(path.sdates.rccs, path.sdates.mrc)
    cols <- c('pt_id', 'pangea_id', 'visit_dt')
    ddates <- rbindlist(lapply(files, fread, select=cols))
    ddates <- unique(ddates)
    ddates <- merge(ddates, aik, by.x='pt_id', by.y='PT_ID')
    ddates[, pt_id := NULL]
    stopifnot(ddates[, uniqueN(pangea_id)==.N])
    setnames(ddates, 'AID', 'aid')

    if(!is.null(select_aid))
        ddates <- ddates[aid %in% select_aid]

    if(get_first_visit)
        ddates <- ddates[, .(date_collection=min(visit_dt)),by='aid']
    ddates
}

get.infection.range.from.testing <- function()
{
    # get maximum and minimum dates

    # for missing 1st pos, use date collection instead (also store contradicting first pos)
    tmp <- meta[, get.sample.collection.dates(aid, get_first_visit = TRUE)]
    meta <- merge(meta, tmp, by='aid', all.x=TRUE)

    # get 15th birthdate and check whether anyone was 100% infected previously
    meta[, date15th := date_birth + as.integer(365.25*15)]
    stopifnot(meta[date15th > date_first_positive, .N == 0])
    # don't min date more than 15 years prior first positive
    # meta[, mean(lowb15lastneg > date_last_negative), ]

    contradict_firstpos_datecoll <<- meta[date_collection < date_first_positive]
    drange <- meta[, .(AID=aid, 
                       MIN=pmax(date_last_negative,  date15th, na.rm=TRUE),
                       MAX=pmin(date_first_positive - 30, date_collection - 30, na.rm=TRUE))]
    drange[, lowb15lastneg := MAX - as.integer(365.25*15) ]
    drange[  lowb15lastneg > MIN, MIN:=lowb15lastneg ]
    drange[, lowb15lastneg := NULL]
    drange

}


#
# test
#

if(0)
{
    args$out.dir <- '~/Dropbox/CONDESAPopSample_phsc_stage1_output_close/'
    args$pkg.dir <- '~/git/phyloscanner/misc_data_analysis_RCCS1519/software'
    args$phylo.dir <- '~/Dropbox/CONDESAPopSample_phsc_stage1_output_close/2022-09-22_phsc_phscrelationships_posthoccount_im_mrca_fixpd/'
    args$date <- '2022-09-22'
}

if(1)
{
    # NEW THINGS HERE
    args$phylo.dir <- '/home/andrea/HPC/project/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI/211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd'
    args$pkg.dir <- '~/git/phyloscanner/misc_data_analysis_RCCS1519/software'
    args$out.dir <- '/home/andrea/HPC/project/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI'
    args$date <- '2021-12-20'
    args$meta.data <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live/RCCS_R15_R18/Rakai_Pangea2_RCCS_Metadata_20220329.RData'
}

#
# Add constants that should not be changed by the user
#
max.per.run <- 4900

# Set default output directories relative to out.dir
args$date <- gsub('-','_',args$date)
.f <- function(x) file.path(args$out.dir, paste0(args$date, x))
args$out.dir.data <- .f('_phsc_input')
args$out.dir.work <- .f('_phsc_work')
args$out.dir.output <- .f('_phsc_output')


cat('Load phyloscanner outputs... \n')
stopifnot(dir.exists(args$phylo.dir))
infiles	<- data.table(F=list.files(args$phylo.dir, pattern='*workspace.rda$', full.names=TRUE))
infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(F)))]
setkey(infiles, PTY_RUN)

# 'dwin' : pairwise relationships between the host in each tree
# 'dc' :  summarises pairwise relationships between hosts across ALL trees. 
# 'a' suffix stands for aggregated
dca <-   infiles[, { cat(PTY_RUN,'\n'); load(F); dc }, by='PTY_RUN']
dwina <- infiles[, { cat(PTY_RUN,'\n'); load(F); dwin }, by='PTY_RUN']
DCA1 <- copy(dca)
DWINA <- copy(dwina)

# Change the format
.format.controls <- function(DT)
{
        DT[, `:=` (CNTRL1=FALSE, CNTRL2=FALSE)]
        DT[ host.1 %like% 'CNTRL-', `:=` (CNTRL1=TRUE, host.1=gsub('CNTRL-', '', host.1))]
        DT[ host.2 %like% 'CNTRL-', `:=` (CNTRL2=TRUE, host.2=gsub('CNTRL-', '', host.2))]
}
.format.controls(dca)
.format.controls(dwina)

# Sort, so host.1 < host.2 
.reorder.host.labels <- function(DT)
{
        tmp <- subset(DT, host.1 > host.2)

        cols1 <- c('host.1','host.2','paths12','paths21','nodes1','nodes2','CNTRL1','CNTRL2')
        cols2 <- c('host.2','host.1','paths21','paths12','nodes2','nodes1','CNTRL2','CNTRL1')
        cols1 <- cols1[cols1 %in% names(DT)]
        cols2 <- cols2[cols2 %in% names(DT)]
        setnames(tmp, cols1, cols2)

        .gs <- function(x)
                gsub('xx','21',gsub('21','12',gsub('12','xx',x)))

        cols <- c('close.and.contiguous.and.directed.cat',
                  'close.and.adjacent.and.directed.cat',
                  'close.and.contiguous.and.ancestry.cat',
                  'close.and.adjacent.and.ancestry.cat', 
                  'type')
        cols <- cols[ cols %in% names(DT)]

        tmp[, (cols) := lapply(.SD, .gs) , .SDcols=cols]

        DT <- rbind(DT[!(host.1>host.2)], tmp)
        return(DT)
}
dca   <- .reorder.host.labels(dca)
dwina <- .reorder.host.labels(dwina)
setkey(dwina, PTY_RUN, host.1, host.2)
setkey(dca, PTY_RUN, host.1, host.2)

# subset to phylogenies with strongest data on linkage 
idx <- dca[ categorisation %in% 'close.and.adjacent.cat',
           .(n.eff=unique(n.eff)),
           by=c('host.1', 'host.2', 'PTY_RUN')]

idx <- idx[ ,
           .(PTY_RUN=PTY_RUN[which.max(n.eff)]),
           by=c('host.1', 'host.2')]
setcolorder(idx, c('PTY_RUN','host.1', 'host.2'))
dca <- merge(idx, dca)
dwina <- merge(idx, dwina)

# 
if(file.exists(args$meta.data))
{   # Now, whenever a pairwise transmission is not supported by the serohistory, we need to "switch" the counts supporting that direction.
    cat('Using serohistory + demographic data to reweight evidence of direction\n')
    
    # Load anon. keys
    file.anonymisation.keys <- file.path(indir.deepanalyses.xiaoyue, 'PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv')
    aik <- fread(file.anonymisation.keys, header = TRUE, select=c('PT_ID', 'AID'))

    # load meta data on infection times
    meta_env <- new.env()
    load(args$meta.data, envir=meta_env)
    meta <- subset(meta_env$meta_data,
                   select=c('aid', 'study_id', 'sex', 'date_birth', 'date_first_positive', 'date_last_negative'))
    meta <- unique(meta[!is.na(aid)])
    # idx <- meta[, uniqueN(comm) == 2, by=aid][V1==TRUE, aid]
    stopifnot(meta[, uniqueN(aid) == .N])

    # IFNEEDED save ids herefor which we do not have a first positive date, so we can ask Joseph.

    # get pairs for which one direction of transmission is unsupported according to demographic and serohistory data.
    drange <- get.infection.range.from.testing()
    dexclude <- dca[, .(host.1, host.2)] |> unique()
    dexclude <- merge(dexclude, drange[, .(host.2=AID, MIN.H2=MIN, MAX.H2=MAX)] , by='host.2', all.x=TRUE)
    dexclude <- merge(dexclude, drange[, .(host.1=AID, MIN.H1=MIN, MAX.H1=MAX)] , by='host.1', all.x=TRUE)
    dexclude[ MAX.H1 <= MIN.H2, EXCLUDE := '21']
    dexclude[ MAX.H2 <= MIN.H1, EXCLUDE := '12']
    dexclude[, table(EXCLUDE, useNA = 'always' )/.N * 100 ]

    # check that for each '12' entry, there is a '21' entry in the count data.table
    fill.in.missing.reverse.direction <- function(DCA)
    {

        # the all.x and all.y options allow us to know which entries weren't previously in DCA
        cols <- c('PTY_RUN', 'host.1', 'host.2', 'categorisation', 'categorical.distance', 'CNTRL1', 'CNTRL2')
        merge(
            DCA[type == '12', ..cols],
            DCA[type == '21', ..cols],
            by=cols, all.x=TRUE, all.y=TRUE
        ) -> tmp
        tmp.12 <- merge(tmp, DCA[type == '12'], by=cols, all.x=TRUE)
        tmp.12[is.na(type), type := '12']
        tmp.21 <- merge(tmp, DCA[type == '21'], by=cols, all.x=TRUE)
        tmp.21[is.na(type), type := '21']
        tmp <- rbind(tmp.12, tmp.21)

        # check
        idx <- tmp[!is.na(n),uniqueN(n), by=cols][ V1 > 1, ..cols]
        stopifnot(nrow(idx)==0)

        # and update: for each "by" group, n & n_eff should be constant.
        # some scores will still be NA
        tmp[, n:=na.omit(n)[1], by=cols] 
        tmp[, n.eff:=na.omit(n.eff)[1], by=cols] 
        tmp[, score:= k.eff/n.eff]
        
        out <- rbind(DCA[! type %in% c('12', '21')], tmp)
        setkey(out, PTY_RUN, host.1, host.2, categorisation, categorical.distance, type)
        out
    }
    dca <- fill.in.missing.reverse.direction(dca)

    # Now, modify dca and dwina accordingly, by removing the counts supporting that direction
    categories.to.update <- c("close.and.adjacent.and.ancestry.cat", "close.and.adjacent.and.directed.cat",
                              "close.and.contiguous.and.ancestry.cat", "close.and.contiguous.and.directed.cat")
    dca_sero_only <- dca[categorisation %in% categories.to.update, ]

    # label.to.update <- function(ctg)
    #   fcase(ctg %like% 'directed', NA_character_, ctg %like% 'ancestry.cat', 'complex.or.no.ancestry')
    update.category.counts <- function(DEXCLUDE, DCA)
    {
        # DEXCLUDE <- copy(dexclude); DCA <- copy(DCA1)
        idx.12 <- DEXCLUDE[EXCLUDE == '12', .(host.1, host.2)]
        idx.21 <- DEXCLUDE[EXCLUDE == '21', .(host.1, host.2)]
        setkey(DCA, host.1, host.2); setkey(dwina, host.1, host.2);
        setkey(idx.12, host.1, host.2);setkey(idx.21, host.1, host.2);
        tmp.12 <- merge(DCA, idx.12, by=c('host.1', 'host.2'))
        tmp.21 <- merge(DCA, idx.21, by=c('host.1', 'host.2'))

        update.category.counts.by.unsupported.direction <- function(ctg, DT, dir)
        {
            .update <- function(TMP)
            {
                with(TMP,{
                    # Find values associated with 'impossible' label
                    # cat(unique(host.1), unique(host.2),'\n')
                    new.dir <- setdiff(c('12', '21'), dir)
                    which.dir <- type == dir
                    which.new.dir <- type == new.dir

                    TMP$k[which.new.dir] <- k1 <- sum(TMP$k)
                    TMP$k.eff[which.new.dir] <- k1.eff <- sum(TMP$k.eff)
                    TMP$k[!which.new.dir] <- 0
                    TMP$k.eff[!which.new.dir] <- 0
                    TMP$score <- TMP$k.eff / TMP$n.eff

                    if(! any(which.new.dir))
                    {
                        tmp_row <- TMP[1, ]
                        tmp_row$type <- new.dir
                        tmp_row$k <- k1
                        tmp_row$k.eff <- k1.eff
                        tmp_row$score <- tmp_row$k.eff / tmp_row$n.eff
                        TMP <- rbind(tmp_row, TMP)
                        return(TMP)
                    }

                    return(TMP)
                })
            }
            cols <- setdiff(names(DT),'categorisation')
            out <- DT[categorisation == ctg, .update(.SD) , by=c('host.1', 'host.2', 'PTY_RUN', 'categorisation'), .SDcols=cols]
            return(out)
        }
        
        
        dca.update.12 <- lapply(categories.to.update, update.category.counts.by.unsupported.direction, tmp.12, '12')
        dca.update.21 <- lapply(categories.to.update, update.category.counts.by.unsupported.direction, tmp.21, '21')
        dca.update <- rbind(rbindlist(dca.update.12), rbindlist(dca.update.21))
        # for some reason some cols are duplicated
        dca.update <- subset(dca.update, select=!duplicated(names(dca.update)))
        setcolorder(dca.update,names(dca))
        dca.update[! k %in% c('12', '21'), sum(k)]

        # lapply(dca.update.12, function(DT) DT[type=='21', sum(k)])
        # lapply(dca.update.21, function(DT) DT[type=='21', sum(k)])

        setkey(dca.update, host.1, host.2, PTY_RUN, categorisation)
        setkey(DCA, host.1, host.2, PTY_RUN, categorisation)
        DCA <- rbind(DCA[!dca.update], dca.update)
    }

    dca_sero_only <- update.category.counts(dexclude, dca_sero_only)
    dca_sero_only[, categorisation:=paste0(categorisation, '.sero')]
    dca <- rbind(dca, dca_sero_only)

    if(0)
    {
        classes <- dca[categorisation %like% '.sero', unique(categorisation)]

        for ( class in classes )
        {
            dca_class <- dca[ categorisation %in% class | categorisation %in% gsub( '.sero', '', class) ]
            dca_class <- dca_class[ type %in% c('12', '21')]
            dca_class[, .N , by=c('host.1', 'host.2')]
            dca_class[, .(N=uniqueN(score), anyNA=any(is.na(score))) , by=c('host.1', 'host.2', 'type')][, table(N, anyNA)]
            any(is.na(c(NA, 1, 2)))
        }
    }

    setkey(dca, host.1, host.2)
    # args$phylo.dir <- '/home/andrea/Downloads'
}

if(0)
{
    # checking that redistribution of count has beem done in the correct direction
    tmp <- dca_sero_only[type %in% c('12', '21') & categorisation %like% 'sero', .(host.1, host.2, type, score)]
    tmp <- merge(tmp, drange[, .(host.1=AID, MIN.1=MIN, MAX.1=MAX)], by='host.1')
    tmp <- merge(tmp, drange[, .(host.2=AID, MIN.2=MIN, MAX.2=MAX)], by='host.2')
    tmp[MAX.1 < MIN.2]
    stopifnot(tmp[MAX.1 < MIN.2, mean(score[type=='12'], na.rm=TRUE) == 1 & mean(score[type=='21'], na.rm=TRUE) == 0])
    stopifnot(tmp[MAX.2 < MIN.1, mean(score[type=='12'], na.rm=TRUE) == 0 & mean(score[type=='21'], na.rm=TRUE) == 1])
}

if(0)
{
    # get pairs whose direction was changed according to serohistory
    class <- 'close.and.adjacent.and.directed.cat.sero'
    dca_class <- dca[categorisation == class | categorisation == gsub( '.sero', '', class)]
    # in the majority of cases, scores sum to 0
    summed_scores <- dca_class[, sum( score, na.rm=T), by=c('host.1', 'host.2', 'categorisation')]
    summed_scores[, table(V1)/.N*100 ]
    # in 8 cases, the sum of the scores went from 1 to 0... (why?)
    summed_scores |>
        dcast(host.1 + host.2 ~ categorisation) |> 
        with( table(close.and.adjacent.and.directed.cat, close.and.adjacent.and.directed.cat.sero)) 
    summed_scores[, table(), by=c('host.1', 'host.2')][, table(V1)]
    
    # subset to cases where summed scores were always one ? 
    idx <- summed_scores[ , all(V1==1) , by=c('host.1', 'host.2') ][V1 == TRUE, .(host.1, host.2)]

    idx <- dca_class[idx][, uniqueN(score), by=c('host.1', 'host.2', 'type')][V1>1, .(host.1, host.2)]
    stopifnot( uniqueN(idx) == nrow(idx)/2 )
    idx <- unique(idx)

    # when changes, mass is moving towards 0
    dca_class[idx][categorisation %like% 'sero', table(score)]

}

# find pairs according to classification rule and thresholds.
# classification rule o: Oliver Ratmann's
# classification rule m: Matthew Hall's
# classification rule b: both

if(args$classif_rule=='o'|args$classif_rule=='b')
{
    dir_group <- dca[categorisation %like% 'close.and.adjacent.and.directed.cat', unique(categorisation)]
    idx <- dir_group %like% 'sero'
    dir_group <- ifelse(any(idx), dir_group[idx], dir_group[1])

    control <- list(linked.group='close.and.adjacent.cat',
                    linked.no='not.close.or.nonadjacent',
                    linked.yes='close.and.adjacent', 
                    dir.group = dir_group,
                    conf.cut=0.6, 
                    neff.cut=3,
                    weight.complex.or.no.ancestry=0.5)
    # Find pairs
    tmp <- find.pairs.in.networks(dwina, dca, control=control, verbose=TRUE)
    dpl <- copy(tmp$network.pairs)
    dc <- copy(tmp$relationship.counts)
    dw <- copy(tmp$windows)

    summarise_serohistory_impact_on_pairs <- function(DC, categ = 'close.and.adjacent.and.directed')
    {
        # DC <- copy(dc); categ = 'close.and.adjacent.and.directed'

        # check that both cat and cat.sero categories are available.
        categ <- gsub( '.sero', '', categ)
        stopifnot( DC[, unique(CATEGORISATION) %like% categ |> sum() == 2 ] ) 

        # subset of interest
        dc_cat <- DC[ CATEGORISATION %like% categ ]
        cat( "in how many cases were the scores changed?\n")
        dc_changed <- dc_cat[, uniqueN(SCORE) > 1, by=c('H1', 'H2', 'TYPE') ]
        dc_changed <- dc_changed[, .(CHANGED_SCORE = any(V1)) , by=c('H1', 'H2')]
        dc_changed[, table(CHANGED_SCORE)] |> knitr::kable() |> print()

        cat( "in how many cases were the directions changed?\n")
        idx <- dc_changed[CHANGED_SCORE == TRUE, .(H1, H2)]
        dc_changed_dir <- dc_cat[idx]
        # dc_changed_dir[is.na(SCORE)]

        # for each pair and class type, evaluate whether score in dir 1->2 is larger
        dc_changed_dir <- dc_changed_dir[ , {
            dir12 <- TYPE %like% '12'
            dir21 <- TYPE %like% '21'
            list(DIR12 = SCORE[dir12] >= SCORE[dir21])
        } , by=c('H1', 'H2', 'CATEGORISATION') ] 

        # for each pair, look whether categorisations disagree.
        dc_changed_dir <- dc_changed_dir[, .(CHANGED_DIR = uniqueN(DIR12) == 2 ) , by=c('H1', 'H2')]
        dc_changed_dir[, table(CHANGED_DIR)] 

        # now let's output a DT with the pairs, indicating whether 
        # - serohistory changed scores
        # - serohistory changed directions! 
        .f <- function(DT1, DT2) merge(DT1, DT2, all.x = TRUE, all.y=TRUE)
        out <- list( unique(dc_cat[, .(H1, H2)]), dc_changed, dc_changed_dir) |> Reduce(f=.f)
        out[is.na(CHANGED_DIR), CHANGED_DIR := FALSE]
        out[, table(CHANGED_SCORE, CHANGED_DIR, useNA='ifany')]

        return(out)
    }

    filename <- file.path(args$phylo.dir, paste0('Rakai_phscnetworks_allpairs_ruleo_sero.rda'))
    if(! file.exists(filename))
    {
        save(dpl, dc, dw, file=filename)
    }

    # Find chains
    tmp <- find.networks(dc, control=control, verbose=TRUE)
    dnet <- copy(tmp$transmission.networks)
    dchain <- copy(tmp$most.likely.transmission.chains)
    filename <- file.path('~/Downloads', paste0('Rakai_phscnetworks_ruleo_sero.rda'))
    if(! file.exists(filename ))
    {
        save(dpl, dc, dw, dnet, dchain, file=filename)
    }
}

if(args$classif_rule=='m'|args$classif_rule=='b')
{
    dir_group <- dca[categorisation %like% 'close.and.adjacent.and.ancestry.cat', unique(categorisation)]
    idx <- dir_group %like% 'sero'
    dir_group <- ifelse(any(idx), dir_group[idx], dir_group[1])

    # Find pairs
    control <- list(linked.group='close.and.adjacent.cat',
                    linked.no='not.close.or.nonadjacent',
                    linked.yes='close.and.adjacent',
                    dir.group="close.and.adjacent.and.ancestry.cat",
                    conf.cut=0.5, 
                    neff.cut=3,
                    weight.complex.or.no.ancestry=0.5)
    tmp <- find.pairs.in.networks(dwina, dca, control=control)
    dpl <- copy(tmp$network.pairs)
    dc <- copy(tmp$relationship.counts)
    dw <- copy(tmp$windows)
    save(dpl, dc, dw, file=file.path('~/Downloads','Rakai_phscnetworks_allpairs_rulem.rda'))
    # Find chains
    control<- list(linked.group='close.and.adjacent.cat',
                   linked.no='not.close.or.nonadjacent',
                   linked.yes='close.and.adjacent',
                   dir.group="close.and.adjacent.and.ancestry.cat", 
                   neff.cut=3, 
                   weight.complex.or.no.ancestry=0.5)
    tmp <- find.networks(dc, control=control, verbose=TRUE)
    dnet <- copy(tmp$transmission.networks)
    dchain <- copy(tmp$most.likely.transmission.chains)
    save(dpl, dc, dw, dnet, dchain, file=file.path(args$phylo.dir,'Rakai_phscnetworks_rulem_sero.rda'))
}

if(!args$classif_rule %in% c('o','m','b'))
{
    stop('Please input --classification_rule as o, m or b')
}
