# kate cannot find round 15 samples on the selected sequences PANGEA IDs.
require(data.table)
library(readxl)

indir.deepdata.rccs <- '/home/andrea/Documents/Box/ratmann_pangea_deepsequencedata/PANGEA2_RCCS'
indir.deepdata.mrc <- '/home/andrea/Documents/Box/ratmann_pangea_deepsequencedata/PANGEA2_MRC'
indir.deepdata <- '/home/andrea/Documents/Box/ratmann_pangea_deepsequencedata'

path.sel.sample.rccs <- file.path(indir.deepdata.rccs, '200422_PANGEA2_RCCS_selected_samples.rds')
path.map.sample.rccs <- file.path(indir.deepdata.rccs, '200422_PANGEA2_RCCS_mapped_samples.rds')
path.mappings.rccs <- file.path(indir.deepdata.rccs, '200316_pangea_mappings_rakai.csv')
path.db.sharing.rccs <- file.path(indir.deepdata.rccs, '200316_pangea_db_sharing_extract_rakai.csv')
path.metadata.rccs <- file.path(indir.deepdata, 'RCCS_R15_R18', 'Rakai_Pangea2_RCCS_Metadata__12Nov2019.csv')

path.sel.sample.mrc <- file.path(indir.deepdata.mrc, '200422_PANGEA2_MRCUVRI_selected_samples.rds')
path.map.sample.mrc <- file.path(indir.deepdata.mrc, '200422_PANGEA2_MRCUVRI_mapped_samples.rds')
path.db.sharing.mrc <- file.path(indir.deepdata.mrc, '200319_pangea_db_sharing_extract_mrc.csv')

path.drugres <- file.path(indir.deepdata, 'RCCS_R15_R18', 'Rakai_Drug_Resistance_20220921.xlsx')

###########
# HELPERS #
###########

.get.mapped.samples <- function()
{
        # Note MRC also has VISIT_DT!
        cols <- c("PANGEA_ID","PREFIX","PID","SEQUENCE_ID")
        tmp.rccs <- readRDS(path.map.sample.rccs)
        tmp.mrc  <- readRDS(path.map.sample.mrc)
        
        names(tmp.mrc)
        names(tmp.rccs)

        tmp.rccs <- subset(tmp.rccs, select=cols)
        tmp.mrc  <- subset(tmp.mrc, select=cols)

        tmp <- unique(rbind(tmp.rccs, tmp.mrc))
        tmp
}

.get.approx.round <- function(sID)
{
        tmp <- gsub('^.*R([0-9][0-9])_.*$','\\1', sID)
        tmp <- gsub('^.*NEU.*$','NEU', tmp)
        tmp <- gsub('^.*PG14.*$','PG14', tmp)
        tmp <- gsub('^.*PG15.*$','PG15', tmp)
        tmp <- gsub('^*UG.*$','UG', tmp)
        tmp
}


.get.visit.dates.and.round <- function(nopangeaidNAs=T)
{
        tmp <- lapply(c(path.db.sharing.rccs, path.db.sharing.mrc), fread)
        tmp <- lapply(tmp, subset, select=c('pangea_id', 'pt_id', 'visit_dt'))
        ddates <- unique(rbindlist(tmp))

        tmp <- fread(path.metadata.rccs)
        tmp <- subset(tmp, select=c('study_id', 'round', 'sample_date'))
        tmp[, `:=` (pt_id=paste0('RK-', study_id), 
                    visit_dt=sample_date, 
                    study_id=NULL,
                    sample_date=NULL)]

        setkey(tmp, pt_id, visit_dt)
        setkey(ddates, pt_id, visit_dt)
        ddates <- merge(ddates, tmp, all.x=T, all.y=T)
        if(nopangeaidNAs)
                ddates <- ddates[!is.na(pangea_id)]
        ddates
}


# Study blood samples
# ___________________

dmap <- .get.mapped.samples()
ddates <- .get.visit.dates.and.round()

dsel <- readRDS(path.sel.sample.rccs) |> as.data.table()
id_sel <- dsel[, unique(PANGEA_ID)]
drounds[ pangea_id %in% id_sel][, table(round)]

# Now do this for codebook 
ddrug <- read_xlsx(path.drugres)
setDT(ddrug)

if(0)
{
        idx <- unique(ddrug$sampleID)
        idx <- gsub('_remap.*?$','', idx)
        .f <- function(x) mean(idx %in% x)
        lapply(dmap, .f)
}

dkate <- ddrug[,.(
                sampleID,
                PREFIX=gsub('_remap.*?$','', sampleID),
                round_approx=.get.approx.round(sampleID)
                )]

.check.by.col <- function(DT, col='PREFIX')
        DT[, .N, by=col][, all(N == 1)]


setkey(dmap, PREFIX); setkey(dkate, PREFIX)
tmp <- merge(dkate, dmap, all.x=TRUE)

tmp[is.na(PANGEA_ID), table(round_approx)]
tmp <- merge(tmp, ddates, by.y='pangea_id', by.x='PANGEA_ID', all.x=T)
tmp[, table(.(round_approx, round), useNA = 'ifany')]

nrow(ddrug); nrow(dkate); nrow(tmp)

tmp[, sum(is.na(round))]
tmp[, table(round)]

if(0) # study participants without PANGEA_ID
{
        filename <- ''
        tmp[is.na(PANGEA_ID), .(sampleID, PREFIX)]
        fwrite()


}
cols <- c('SEQUENCE_ID', 'PREFIX')

tmp[! SEQUENCE_ID %in% PREFIX, ..cols]
tmp[, sum(!is.na(visit_dt))]

# save output
cols <- c('PREFIX', 'PANGEA_ID', 'round', 'visit_dt')
tmp0 <- subset(tmp, select=cols)
names(tmp0) <- toupper(cols)
filename=file.path(indir.deepdata, '221007_prefixes2pangeaid.csv')
if('visit_dt' %in% cols)
        filename <- gsub('.csv', '_withdates.csv', filename)
fwrite(tmp0, file.path(indir.deepdata, '221007_prefixes2pangeaid.csv'))


