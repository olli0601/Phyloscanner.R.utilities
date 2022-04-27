cat('\n\n=====  TSI_run_predictions.R =====\n\n')

library(lubridate)
library(data.table)

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
    "--TSI_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to HIV-phylo-TSI-main repository", 
    dest = 'TSI.dir'
  ),
  optparse::make_option(
    "--env_name",
    type = "character",
    default = 'hivphylotsi',
    help = "Conda environment name to run HIV-phylo-TSI analyses [default] ",
    dest = 'env_name'
  ),
  optparse::make_option(
    "--input_samples",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to input samples rds containing PANGEA_IDs and RENAME_IDs", 
    dest = 'phsc.samples'
  ),
  optparse::make_option(
    "--walltime",
    type = "integer",
    default = 4L,
    help = "Job walltime (hours) [default]",
    dest = 'walltime'
  ),
  optparse::make_option(
    "--memory",
    type = "integer",
    default = 2L,
    help = "Job memory (GB) [default]",
    dest = 'memory'
  ),
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

args <-
  optparse::parse_args(optparse::OptionParser(option_list = option_list))


###############
# helper f's
###############

.unzip.patstats <- function(x){
        csv.name <- unzip(x, list = TRUE)$Name
        csv.name <- grep('_patStats.csv$',csv.name,value = T)
        patstat <- data.table(read.csv(unz(x, csv.name), header = TRUE, sep = ","))
        csv.name <- file.path(dirname(x), csv.name)
        write.csv(patstat, file=csv.name)
        return(csv.name)
}

generate.sample <- function(maf) {

        cat('\n', 'Running generate.sample()...', '\n')
        # Lele's function to compute Minor Allele Frequencies. 
        MAF_matrix<-matrix(0,ncol=10000,nrow=(nrow(maf)+1))
        MAF_matrix[1,]<-seq(1,10000)
        rownames(MAF_matrix)<-c('pos', sort(maf$SAMPLE_ID))
        names_nobf <- maf[HXB2_EXISTS==FALSE, SAMPLE_ID] 
        idx_nobf <- which(rownames(MAF_matrix) %in% names_nobf)
        MAF_matrix[idx_nobf, ] <- NA
        spls <- names(which(!is.na(MAF_matrix[,1])))
        spls <- spls[spls != 'pos']

        for(sp in spls){
                sp_file <- maf[SAMPLE_ID == sp, HXB2_PATH]
                # cat(paste0('Computing MAF for: ', sp_file,'\n'))
                basefile<-read.csv(sp_file, header=T, stringsAsFactors=FALSE)
                indexes_HXB2_pos<-which(basefile$Position.in.B.FR.83.HXB2_LAI_IIIB_BRU.K03455!='-')
                sample_HXB2_pos<-as.numeric(basefile[indexes_HXB2_pos,1])
                sample_MAFs<-apply(basefile[indexes_HXB2_pos,4:7], 1, function(x) 1-(max(x,na.rm = T)/sum(x,na.rm = T)))
                sample_MAFs<-as.numeric(gsub(NaN, 0,sample_MAFs))
                MAF_matrix[sp,sample_HXB2_pos]<-sample_MAFs 
        }

        return(MAF_matrix)
}

get.sampling.dates <- function(phsc.samples = args$phsc.samples)
{
        # Find files containing all sample collection dates 
        if(user != 'andrea')
        {
                db.sharing.path.rccs <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv'
                db.sharing.path.mrc <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv'
        }else{
                db.sharing.path.rccs <- '/home/andrea/Documents/Box/ratmann_pangea_deepsequencedata/PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv'
                db.sharing.path.mrc <-  '/home/andrea/Documents/Box/ratmann_pangea_deepsequencedata/PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv' 
        }

        tmp <- c(db.sharing.path.rccs,db.sharing.path.mrc, phsc.samples)
        stopifnot(all(file.exists(tmp)))

        dsamples <- setDT(readRDS(phsc.samples))
        dsamples <- unique(dsamples[, .(PANGEA_ID, RENAME_ID)])
        dsamples[, PANGEA_ID:=gsub('^.*?_','',PANGEA_ID)]

        ddates <- setDT(read.csv(db.sharing.path.mrc))
        ddates <- unique(ddates[, .(pangea_id, visit_dt)])
        tmp <- as.data.table(read.csv(db.sharing.path.rccs))
        tmp <- unique(tmp[, .(pangea_id, visit_dt)])
        ddates <- rbind(tmp, ddates)
        ddates[, visit_dt:=as.Date(visit_dt, format="%Y-%m-%d")]
        stopifnot(ddates[, anyDuplicated(pangea_id) == 0,])

        ddates <- merge(dsamples, ddates, all.x=TRUE,
                        by.x='PANGEA_ID', by.y='pangea_id')
        ddates[, PANGEA_ID := NULL]
        
        # Order based on sampling dates 
        setnames(ddates, 'RENAME_ID', 'SAMPLE_ID')
        ddates[, AID := gsub('-fq.*?$','', SAMPLE_ID)]
        setorder(ddates, AID, -visit_dt)

        return(ddates)
}

write.mafs.and.cmds <-function(pty_idx)
{
        # For each PTY index, computes the MAF matrix and stores it in the phscrel directory
        # Also, writes bash command to set up HIV-phylo-TSI

        cat('Processing PTY index: ', pty_idx, '...\n')

        # Load patstats
        files_pty <- as.vector(dfiles[pty==pty_idx]) 
        patstats <- as.data.table(read.csv(files_pty$pat.path, header=TRUE, sep=",", stringsAsFactors=F))
        if(any((patstats$xcoord %% 1) != 0)){
                patstats[, xcoord:=ceiling(xcoord)]
                write.csv(patstats, files_pty$pat.path)
        }
        ph.input <- as.data.table(read.csv(files_pty$phi.path, header=FALSE, sep=",",stringsAsFactors=F))
        colnames(ph.input) <- c('BAM_PATH','FASTA_PATH', 'SAMPLE_ID')
        ph.input[, AID:=gsub('-fq.*?$','', SAMPLE_ID)]
        stopifnot( all(patstats[, unique(host.id)] %in% ph.input[, unique(AID)]) )

        # (Checked that BAM and FASTA files are consistent in terms of namings and locs)

        # Find BAM_PATH then get the MAF
        ph.input[, HXB2_PATH := gsub('.bam$','_BaseFreqs_WithHXB2.csv', basename(BAM_PATH))]
        if(user == 'andrea')
        {
                fraser.dir <- '~/Documents/Box/ratmann_pangea_deepsequencedata/fraserupload-KoPKvNP7dgmnfnxE/well/fraser/DATA/processing'
                hxb2files <- as.character(list.files(fraser.dir, recursive=TRUE, pattern='^.*?HXB2.csv', full.names = TRUE))
                hxb2files <- data.table(FULL=hxb2files, BASE=basename(hxb2files))
                ph.input <- merge(ph.input, hxb2files, by.x='HXB2_PATH', by.y='BASE', all.x=TRUE)
                ph.input[, `:=`(HXB2_PATH=FULL, FULL=NULL)]
                ph.input[, HXB2_EXISTS:=file.exists(HXB2_PATH) ]
        }else{
                ph.input[, HXB2_PATH := file.path(dirname(BAM_PATH), HXB2_PATH)]
                ph.input[, HXB2_EXISTS := file.exists(HXB2_PATH)]
        }
        maf <- ph.input[, .(SAMPLE_ID, HXB2_PATH, HXB2_EXISTS)]
        maf_mat <- generate.sample(maf)
        cat(maf_mat[2, 1:10], '\n')

        # If there are multiple sequences associated to one AID:
        # take sequence with associated BaseFreq file ("HXB2")
        # with latest collection date  
        tmp1 <- gsub('-fq.*?$','',rownames(maf_mat))
        tmp1 <- unique(tmp1[duplicated(tmp1)])
        tmp1 <- data.table(AID = tmp1)
        if(tmp1[, .N>1])
        {
                tmp1 <- tmp1[, list(FQ=grep(AID, rownames(maf_mat), value=T)),by=AID]
                # which do not have HXB2?
                tmp1 <- tmp1[, list(HXB2 = !is.na(maf_mat[FQ, 1])), by=c("AID","FQ") ]
                tmp1 <- merge(tmp1, ddates, by.x=c('AID', 'FQ'), by.y=c("AID", "SAMPLE_ID"))
                setorder(ddates, AID, -visit_dt)
                tmp1 <- tmp1[, {
                        z <- which(HXB2 == TRUE)[1];
                        z <- ifelse(is.na(z), 1, z)
                        list(FQ=FQ[z], visit_dt=visit_dt[z])
                }, by='AID']


                rows_to_del <- rownames(maf_mat)[which(grepl(paste0(tmp1$AID, collapse='|'), rownames(maf_mat) ))]
                rows_to_del <- rows_to_del[! rows_to_del %in% tmp1$FQ]
                
                # Store -fq used so we can check exact sampling date
                maf_mat <- maf_mat[! rownames(maf_mat) %in% rows_to_del, ]
        }else{
                rows_to_del <- c()
        }
        filename=paste0('ptyr', pty_idx, '_basefreqs_used.csv')

        write.csv(rownames(maf_mat), 
                  file=file.path(dirname(files_pty$pat.path), filename))

        rownames(maf_mat) <- gsub('-fq[0-9]$', '', rownames(maf_mat))
        rm(tmp1, rows_to_del)
        maf_file=paste0('ptyr', pty_idx, '_maf.csv')
        maf_file=file.path(dirname(files_pty$pat.path), maf_file)
        cat(maf_mat[2, 1:10], '\n')
        write.table(maf_mat, file=maf_file, quote=F, col.names=FALSE, sep=',')


        # writing command
        #________________
        cmd <- paste0('python ',args$TSI.dir,'/HIVPhyloTSI.py \\\n',
                      ' -d ', args$TSI.dir, '/Model \\\n',
                      ' -p ', files_pty$pat.path,' \\\n',
                      ' -m ', files_pty$maf.path,' \\\n',
                      ' -o ', files_pty$tsi.path, '\n'
                     )
        dfiles[pty==pty_idx, CMD:=cmd]

        # guess I may need to free up some memory here
        rm(maf_mat, ph.input, patstats)
        gc()
}

###############
# testing
###############
user <- Sys.info()[['user']]
if (user=='andrea') {
        args <- list(
                out.dir='~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/19037_phsc_output',
                rel.dir= "~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/19037_phsc_phscrelationships_seed_42_blacklist_report_TRUE_distance_threshold_1_min_reads_per_host_1_multinomial_TRUE_outgroup_name_BFR83HXB2_LAI_IIIB_BRUK03455_output_nexus_tree_TRUE_ratio_blacklist_threshold_0005_skip_summary_graph_TRUE/",
                TSI.dir='~/git/HIV-phyloTSI-main',
                env_name='hivphylotsi',
                walltime = 3L,
                memory = 2L,
                controller=NA,
                phsc.samples="~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210120_RCCSUVRI_phscinput_samples.rds"
        )
}


if(0)
{
        args <- list(
                out.dir='/rds/general/project/ratmann_deepseq_analyses/live/seroconverters2/2022_04_25_phsc_output',
                rel.dir= "/rds/general/project/ratmann_deepseq_analyses/live/seroconverters2/2022_04_25_phsc_phscrelationships_sd_42_blacklist_report_T_mr_1_og_A1UGANDA2007p191845JX236671_output_nexus_tree_T_rtt_0005_skip_summary_graph_T_sdt_1",
                TSI.dir='~/git/HIV-phyloTSI-main',
                env_name='hivphylotsi',
                walltime = 3L,
                memory = 2L,
                controller=NA,
                phsc.samples="/rds/general/project/ratmann_deepseq_analyses/live/seroconverters2/210419_phscinput_samples.rds"
        )
}

################
# main
################

args$date <- gsub('-','_',args$date)
if( ! grepl('output$', args$out.dir))
{
        args$out.dir <- 
                file.path(args$out.dir, paste0(args$date, "_phsc_output"))
}

work.dir <- gsub('_output$','_work',args$out.dir)

stopifnot(dir.exists(work.dir))
stopifnot(dir.exists(args$rel.dir))
stopifnot(dir.exists(args$out.dir))
stopifnot(file.exists(args$phsc.samples))

# dates of collection for samples.
ddates <- get.sampling.dates(phsc.samples=args$phsc.samples)

# Collect pty files allowing to run Tanya's algorithm
patstats_zipped <- list.files(args$rel.dir, pattern='zip$', full.name=TRUE)
tmp <- gsub('^.*?ptyr|_otherstuff.zip','',patstats_zipped)
phsc_inputs <- file.path(args$out.dir, 
                         paste0('ptyr', tmp, '_trees'), 
                         paste0('ptyr', tmp, '_input.csv' ))

dfiles <- data.table(pty=as.integer(tmp), 
                     zip.path=patstats_zipped, 
                     phi.path=phsc_inputs)
setkey(dfiles, pty)
dfiles <- dfiles[ file.exists(zip.path) & file.exists(phi.path)]
dfiles[, pat.path:=.unzip.patstats(zip.path), by='pty']
dfiles[, maf.path:=gsub('patStats.csv$','maf.csv',pat.path)]
dfiles[, tsi.path:=gsub('patStats.csv$','tsi.csv',pat.path)]
dfiles[, CMD:=NA_character_]

lapply(dfiles$pty, write.mafs.and.cmds)

##################################
# submit jobs 
##################################

dfiles[, IDX := 1:.N, ]
dfiles[, CMD := paste0(IDX, ')\n',CMD, ';;\n')]
cmd <- paste0(dfiles$CMD, collapse='\n')
cmd <- paste0('case $PBS_ARRAY_INDEX in\n', cmd, 'esac\n')

cmd <- paste0('\n module load anaconda3/personal \n',
              'source activate ',args$env_name, '\n',
              cmd)

#make header
header <- paste0(
                 "#!/bin/sh \n",
                 "#PBS -l walltime=",  args$walltime,":00:00,pcput=",args$walltime - 1,":50:00 \n",
                 "#PBS -l select=1:ncpus=1:mem=",  args$memory,  "gb \n",
                 "#PBS -j oe \n",
                 "#PBS -J 1-", dfiles[, max(IDX)] ,"\n"
)

# Patch together and write 
cmd <- paste0(header, cmd)
datetime <- paste(strsplit(date(), split = ' ')[[1]], collapse = '_', sep = '')
datetime <- gsub(':','',datetime)
outfile <- file.path(work.dir, paste0('phylo_tsi_',datetime,'.sh'))
cat(cmd, file=outfile)

# Run the command
cmd <- paste("cd ",dirname(outfile),'\n',"qsub ", outfile)
cat(cmd)
cat(system(cmd, intern = TRUE))
