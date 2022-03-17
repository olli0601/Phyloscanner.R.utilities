cat('\n\n=====  TSI_generate_maf.R =====\n\n')
# TODO: note we should redo some HXB2 pushing

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
    "--walltime",
    type = "integer",
    default = 3L,
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

.read <- function(x){
  if(grepl('.csv$', x)){return(as.data.table(read.csv(x)))}
  if(grepl('.rds$|.RDS$',x)){return(as.data.table(readRDS(x)))}
}

.unzip.patstats <- function(x){
        csv.name <- unzip(x, list = TRUE)$Name
        csv.name <- grep('_patStats.csv$',csv.name,value = T)
        patstat <- data.table(read.csv(unz(x, csv.name), header = TRUE, sep = ","))
        csv.name <- file.path(dirname(x), csv.name)
        write.csv(patstat, file=csv.name)
        return(csv.name)
}

generate_sample <- function(maf) {

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
                cat(paste0('Computing MAF for: ', sp_file,'\n'))
                basefile<-read.csv(sp_file, header=T)
                indexes_HXB2_pos<-which(basefile$Position.in.B.FR.83.HXB2_LAI_IIIB_BRU.K03455!='-')
                sample_HXB2_pos<-as.numeric(basefile[indexes_HXB2_pos,1])
                sample_MAFs<-apply(basefile[indexes_HXB2_pos,4:7], 1, function(x) 1-(max(x,na.rm = T)/sum(x,na.rm = T)))
                sample_MAFs<-as.numeric(gsub(NaN, 0,sample_MAFs))
                MAF_matrix[sp,sample_HXB2_pos]<-sample_MAFs 
        }

        return(MAF_matrix)
}

###############
# testing
###############
if (Sys.info()[['user']]=='andrea') {
        args <- list(
                out.dir='~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/19037_phsc_output',
                rel.dir= "~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/19037_phsc_phscrelationships_seed_42_blacklist_report_TRUE_distance_threshold_1_min_reads_per_host_1_multinomial_TRUE_outgroup_name_BFR83HXB2_LAI_IIIB_BRUK03455_output_nexus_tree_TRUE_ratio_blacklist_threshold_0005_skip_summary_graph_TRUE/",
                TSI.dir='~/git/HIV-phyloTSI-main',
                env_name='hivphylotsi',
                walltime = 3L,
                memory = 2L,
                controller=NA,
                file.bf.locs="~/Documents/Box/2021/phyloTSI/bfloc2hpc_20220103.rds"
        )
}

################
# main
################

work.dir <- gsub('_output$','_work',args$out.dir)
stopifnot(dir.exists(work.dir))
stopifnot(dir.exists(args$rel.dir))
stopifnot(dir.exists(args$out.dir))

# These determine the pty's for which we can run Tanya's algorithm
patstats_zipped <- list.files(args$rel.dir, pattern='zip$', full.name=TRUE)
phsc_inputs <- list.files(args$out.dir, pattern='input.csv$', full.name=TRUE, recursive=TRUE)

# dfiles containing paths of interest + unzipping PatStats 
tmp <- gsub('^.*?ptyr|_otherstuff.zip','',patstats_zipped)
dfiles <- data.table(pty=as.integer(tmp))
setkey(dfiles, pty)
dfiles[, zip.path:=grep( paste0('ptyr', pty, '_'), patstats_zipped, value=T) ,by=pty]
dfiles[, phi.path:=grep( paste0('ptyr', pty, '_'), phsc_inputs, value=T) ,by=pty]
dfiles <- dfiles[ file.exists(zip.path) & file.exists(phi.path)]
dfiles[, pat.path:=.unzip.patstats(zip.path), by='pty']
dfiles[, maf.path:=gsub('patStats.csv$','maf.csv',pat.path)]
dfiles[, tsi.path:=gsub('patStats.csv$','tsi.csv',pat.path)]
dfiles[, CMD:=NA_character_]

for (pty_idx in dfiles$pty)
{
        # Load patstats
        files_pty <- as.vector(dfiles[pty==pty_idx]) 
        patstats <- as.data.table(read.csv(files_pty$pat.path, header=TRUE, sep=","))
        if(any((patstats$xcoord %% 1) != 0)){
                patstats[, xcoord:=ceiling(xcoord)]
                write.csv(patstats, files_pty$pat.path)
        }
        ph.input <- as.data.table(read.csv(files_pty$phi.path, header=FALSE, sep=","))
        colnames(ph.input) <- c('BAM_PATH','FASTA_PATH', 'SAMPLE_ID')
        ph.input[, AID:=gsub('-fq.*?$','', SAMPLE_ID)]
        stopifnot( all(patstats[, unique(host.id)] %in% ph.input[, unique(AID)]) )

        if(0) # Check that bam and fasta path make sense
        {
                tmp1 <- ph.input[, .(BAM_PATH, FASTA_PATH) ]
                tmp1[, lapply(.SD, dirname)][BAM_PATH != FASTA_PATH]
                tmp1[, colnames(tmp1):=lapply(.SD, basename)]
                .f <- function(x){gsub('_ref.fasta$|.bam$|.fasta$', '', x)}
                tmp1[, colnames(tmp1):=lapply(.SD, .f)]
                tmp1[BAM_PATH != FASTA_PATH]
        }

        # Find BAM_PATH then get the MAF
        ph.input[, HXB2_PATH := gsub('.bam$','_BaseFreqs_WithHXB2.csv', basename(BAM_PATH))]
        if(Sys.info()[['user']] == 'andrea')
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
        maf <- maf[, .(SAMPLE_ID, HXB2_PATH, HXB2_EXISTS)]
        maf_mat <- generate_sample(maf)

        # If there are multiple sequences associated to one AID
        # take first non-NA -fq[0-9]
        # then save the maf 
        tmp1 <- gsub('-fq.*?$','',rownames(maf_mat))
        tmp1 <- unique(tmp1[duplicated(tmp1)])
        tmp1 <- data.table(AID = tmp1)
        tmp1 <- tmp1[, list(FQ=grep(AID, rownames(maf_mat), value=T)),by=AID]
        # remove the NAs
        tmp1 <- tmp1[, list(IS_NA = is.na(maf_mat[FQ, 1])), by=c("AID","FQ") ]
        tmp1[, IDX:=which(IS_NA == TRUE)[1] , by="AID"]
        tmp1[is.na(IDX), IDX:=1]
        tmp1 <- unique(tmp1[, list(FQ=FQ[IDX]), by='AID'])

        rows_to_del <- rownames(maf_mat)[which(grepl(paste0(tmp1$AID, collapse='|'), rownames(maf_mat) ))]
        rows_to_del <- rows_to_del[! rows_to_del %in% tmp1$FQ]
        
        # Store -fq used so we can check exact sampling date
        maf_mat <- maf_mat[! rownames(maf_mat) %in% rows_to_del, ]
        filename=paste0('ptyr', pty_idx, '_basefreqs_used.csv')
        write.csv(rownames(maf_mat), 
                  file=file.path(dirname(files_pty$pat.path), filename))
        rownames(maf_mat) <- gsub('-fq[0-9]$', '', rownames(maf_mat))
        rm(tmp1, rows_to_del)

        maf_file=paste0('ptyr', pty_idx, '_maf.csv')
        maf_file=file.path(dirname(files_pty$pat.path), maf_file)
        write.table(maf_mat, file=maf_file, quote=F, col.names=FALSE, sep=',')

        cmd <- paste0('python ',args$TSI.dir,'/HIVPhyloTSI.py \\\n',
                      ' -d ', args$TSI.dir, '/Model \\\n',
                      ' -p ', files_pty$pat.path,' \\\n',
                      ' -m ', files_pty$maf.path,' \\\n',
                      ' -o ', files_pty$tsi.path, '\n'
                     )
        # guess I may need to free up some memory here
        rm(maf_mat, ph.input, patstats)
        gc()
}


##################################
# submit jobs 
##################################

dfiles[, IDX := 1:.N, ]
dfiles[, CMD := paste0(IDX, ')\n',CMD, ';;\n')]
cmd <- paste0(dfiles$CMD, collapse='\n')
cmd <- paste0('case $PBS_array_index in\n', cmd, 'esac\n')

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
