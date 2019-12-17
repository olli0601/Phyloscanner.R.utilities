get.files.to.share<- function()
{
	home <- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis'
	infile <- file.path(home,'Rakai_phyloscanner_170704_stagethree.csv')
	df <- as.data.table(read.csv(infile, stringsAsFactors=FALSE))
	
	random.sample <- c(194, 118, 220, 197,  35, 336,  41, 190,  96,  40)
	df <- subset(df, PTY_RUN%in%random.sample)
	
	bamdir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA_mapout/data'
	outdir <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA_RCCS_comparisonstudy'
	
	dd <- data.table(F=list.files(bamdir))
	
}