PR.PACKAGE			<- "phyloscan" 
PR.read.processed.phyloscanner.output.in.directory		<- paste('Rscript', system.file(package=PR.PACKAGE, "phsc.read.processed.phyloscanner.output.in.directory.Rscript"))


phsc.cmd.mafft.add<- function(infile, reffile, outfile, options='')
{
	#mafft --reorder --anysymbol --add new_sequences --auto input
	tmp		<- c( 	gsub(' ','\\ ',gsub('(','\\(',gsub(')','\\)',infile,fixed=T),fixed=T),fixed=T),
			gsub(' ','\\ ',gsub('(','\\(',gsub(')','\\)',reffile,fixed=T),fixed=T),fixed=T),
			gsub(' ','\\ ',gsub('(','\\(',gsub(')','\\)',outfile,fixed=T),fixed=T),fixed=T)
	)
	cmd		<- paste('mafft --anysymbol ',options,' --add ',tmp[1],' --auto ',tmp[2],' > ',tmp[3], sep='')
	cmd
}

phsc.cmd.blacklist.reads<- function(pr, inputFileName, outputFileName, rawThreshold=1, ratioThreshold=1/200, tipRegex=NA)
{
	cmd	<- paste(pr, rawThreshold, ratioThreshold, inputFileName, outputFileName)
	if(!is.na(tipRegex))
		cmd	<- paste(cmd, '--tipRegex', tipRegex)
	cmd
}

phsc.cmd.blacklist.dualinfections<- function(pr, inputFileNameDuals, outputFileName, blacklistFileName=NA, summaryFileName=NA, treeFileName=NA, dual.prop.threshold=0.01, windowCount=NA)
{
	cmd	<- paste(pr, dual.prop.threshold, inputFileNameDuals, outputFileName)
	if(!is.na(treeFileName))
		cmd	<- paste(cmd, '--treePrefix', treeFileName)		
	if(!is.na(windowCount))
		cmd	<- paste(cmd, '--windowCount', windowCount)
	if(!is.na(blacklistFileName))
		cmd	<- paste(cmd, '--existingBlacklistsPrefix', blacklistFileName)
	if(!is.na(summaryFileName))
		cmd	<- paste(cmd, '--summaryFile', summaryFileName)
	cmd
}

phsc.cmd.blacklist.parsimonybased<- function(pr, scriptdir, inputFileName, outputFileName, dualCandidatesOutputFileName=NA, rawThreshold=1, ratioThreshold=1/200, sankhoffK=20, outgroupName=NA, tipRegex=NA)
{
	cmd	<- paste(pr, '--scriptdir',scriptdir, rawThreshold, ratioThreshold, sankhoffK, inputFileName, outputFileName)
	if(!is.na(dualCandidatesOutputFileName))
		cmd	<- paste(cmd, '--dualsOutputFile', dualCandidatesOutputFileName)			
	if(!is.na(outgroupName))
		cmd	<- paste(cmd, '--outgroupName', outgroupName)		
	if(!is.na(tipRegex))
		cmd	<- paste(cmd, '--tipRegex', tipRegex)
	cmd
}

phsc.cmd.blacklist.rogue.geneticdistance<- function(pr, scriptdir, inputTreeFileName, outputFileName, longestBranchLength=0.04, dropProportion=1/100, blackListFileName=NA, outgroupName=NA, tipRegex=NA)
{
	cmd	<- paste(pr, '--scriptdir',scriptdir, dropProportion, longestBranchLength, inputTreeFileName, outputFileName)
	if(!is.na(tipRegex))
		cmd	<- paste(cmd, '--tipRegex', tipRegex)
	if(!is.na(outgroupName))
		cmd	<- paste(cmd, '--outgroupName', outgroupName)	
	if(!is.na(blackListFileName))
		cmd	<- paste(cmd, '--blacklist', blackListFileName)
	cmd
}

phsc.cmd.blacklist.rogue.extremeprob<- function(pr, scriptdir, inputTreeFileName, outputFileName, probThreshold=0.001, longestBranchLength=0.04, dropProportion=1/100, blackListFileName=NA, outgroupName=NA, tipRegex=NA)
{
	cmd	<- paste(pr, '--scriptdir',scriptdir, dropProportion, probThreshold, longestBranchLength, inputTreeFileName, outputFileName)
	if(!is.na(tipRegex))
		cmd	<- paste(cmd, '--tipRegex', tipRegex)
	if(!is.na(outgroupName))
		cmd	<- paste(cmd, '--outgroupName', outgroupName)	
	if(!is.na(blackListFileName))
		cmd	<- paste(cmd, '--blacklist', blackListFileName)
	cmd
}

phsc.cmd.blacklist.downsample<- function(pr, scriptdir, inputTreeFileName, outputFileName, maxReadsPerPatient=200, blackListFileName=NA, outgroupName=NA, tipRegex=NA)
{
	cmd	<- paste(pr, ' --scriptdir ',scriptdir, maxReadsPerPatient, inputTreeFileName, outputFileName)
	if(!is.na(outgroupName))
		cmd	<- paste(cmd, '--outgroupName', outgroupName)	
	if(!is.na(blackListFileName))
		cmd	<- paste(cmd, '--blacklist', blackListFileName)
	cmd
}


phsc.cmd.SummaryStatistics<- function(pr, scriptdir, outgroupName, file.patients, treeFiles, splitsFiles, outputBaseName, blacklistFiles=NA)
{
	cmd	<- paste(pr,' --scriptdir ',scriptdir,' --outgroupName ', outgroupName, ' "', file.patients, '" "', treeFiles, '" "',splitsFiles, '" "',outputBaseName, '"', sep='')
	if(!is.na(blacklistFiles))
		cmd	<- paste(cmd, ' --blacklists "',blacklistFiles,'"', sep='')
	cmd
}


phsc.cmd.read.processed.phyloscanner.output.in.directory<- function(prefix.infiles, save.file.base, read.likelytransmissions=TRUE, read.trees=TRUE, read.subtrees=TRUE, resume=FALSE, zip=FALSE, pr=PR.read.processed.phyloscanner.output.in.directory)	
{
	cmd	<- paste(pr,' --prefix.infiles "',prefix.infiles,'" --save.file.base "',save.file.base,'"',sep='')
	if(read.likelytransmissions)
		cmd	<- paste(cmd, '--read.likelytransmissions')
	if(read.trees)
		cmd	<- paste(cmd, '--read.trees')
	if(read.subtrees)
		cmd	<- paste(cmd, '--read.subtrees')	
	if(resume)
		cmd	<- paste(cmd,'--resume')
	if(zip)
		cmd	<- paste(cmd,'--zip')
	cmd
}

phsc.cmd.SplitPatientsToSubtrees<- function(pr, scriptdir, infile, outputdir=NA, blacklistFiles=NA, outputFileIdentifier=NA, outgroupName=NA, splitsRule=NA, sankhoff.k=NA, pdfwidth=30, pdfrelheight=0.15)	
{
	cmd	<- paste(pr, ' --inputFile "', infile, '"', ' --scriptdir "', scriptdir,'"', sep='')
	if(!is.na(blacklistFiles))
		cmd	<- paste(cmd, ' --blacklist "', blacklistFiles,'"', sep='')
	if(!is.na(outputdir))
		cmd	<- paste(cmd, ' --outputdir "', outputdir,'"', sep='')	
	if(!is.na(outputFileIdentifier))
		cmd	<- paste(cmd, ' --outputfileid "', outputFileIdentifier, '"', sep='')
	if(!is.na(outgroupName))
		cmd	<- paste(cmd, ' --outgroupName ', outgroupName, sep='')
	if(!is.na(splitsRule))
		cmd	<- paste(cmd, ' --splitsRule ', splitsRule, sep='')
	if(!is.na(sankhoff.k))
		cmd	<- paste(cmd, ' --kParam ', sankhoff.k, sep='')	
	if(!is.na(pdfwidth))
		cmd	<- paste(cmd, ' --pdfwidth ', pdfwidth, sep='')
	if(!is.na(pdfrelheight))
		cmd	<- paste(cmd, ' --pdfrelheight ', pdfrelheight, sep='')
	cmd	
}
	
phsc.cmd.LikelyTransmissions<- function(pr, scriptdir, file.tree, file.splits, file.out)
{
	cmd<- paste(pr, ' "', file.tree, '" "', file.splits, '" "', file.out,'" ',' --scriptdir ', scriptdir, sep='')
	cmd
}

phsc.cmd.LikelyTransmissionsSummary<- function(pr, scriptdir, file.patients, file.summary, file.lkl, file.out, file.detailed.out=NA, min.threshold=1, allow.splits=FALSE)
{
	cmd<- paste(pr, ' "',file.patients, '" "', file.lkl, '" "', file.out,'" --scriptdir ', scriptdir, ' --summaryFile "', file.summary, '" --minThreshold ', min.threshold, sep='')
	if(!is.na(file.detailed.out))
		cmd	<- paste(cmd, ' --detailedOutput "', file.detailed.out, '"', sep='')
	if(allow.splits)
		cmd	<- paste(cmd, ' --allowSplits', sep='')	
	cmd
}

#' @export
#' @title Combine output from multiple phyloscanner runs
#' @param in.dir Full directory path to output of phyloscanner runs
#' @param save.file If not NA, the combined output is saved to file.
#' @param postfix.trees Pattern at end of file names to identify short read tree output generated by the phyloscanner tools. By default 'trees.rda'.
#' @param postfix.trmwindowstats Pattern at end of file names to identify transmission statistics output per window, generated by the phyloscanner tools. By default 'trmStatsPerWindow.rda'.
#' @param regex.ind Regular expression that identifies the ID of query individuals
#' @param trmw.min.reads Minimum number of reads per individuals in order to make a transmission assignment that involves this individual
#' @return list of three R objects 'phs', 'dtrees', 'dtrms'. See description.
#' @description This function generates three R objects. 'phs' is a list of all short read trees in ape format.
#' 	  'dtrees' is a data.table that provides info on each short read tree. Columns are 'PTY_RUN' (phyloscanner run id), 'W_FROM' (start of window), 'W_TO' (end of window), 'IDX' (index of tree in phs).
#' 	  'dtrms' is a data.table that provides info on transmission assignments. Columns are 'PTY_RUN' (phyloscanner run id), 'ID1' (identifier of first individual), 'ID2' (identifier of second individual), 'TYPE' (window assignment), 
#'	  'WIN_OF_TYPE' (number of windows with that assignment), 'WIN_TOTAL' (Number of windows where reads from both individuals are present), 'PAIR_ID' (unique ID of each combination of PTY_RUN,ID1, ID2)
phsc.combine.phyloscanner.output<- function(in.dir, save.file=NA, postfix.trees='trees.rda', postfix.trmwindowstats='trmStatsPerWindow.rda', regex.ind="^[0-9]+_[0-9]_[0-9]+", trmw.min.reads=1, trmw.min.tips=1, trmw.close.brl=Inf, trmw.distant.brl=Inf)
{
	#	read trees
	cat('\nread trees')
	ptyr	<- data.table(FILE=list.files(in.dir, pattern=postfix.trees, full.names=TRUE))	
	phs		<- lapply(seq_len(nrow(ptyr)), function(i)
			{
				load(ptyr[i, FILE])
				phs
			})
	phs		<- do.call(c, phs)
	#	read tree info
	dtrees	<- lapply(seq_len(nrow(ptyr)), function(i)
			{
				load(ptyr[i, FILE])
				dfr
			})
	dtrees	<- do.call('rbind', dtrees)
	dtrees[, IDX:=seq_len(nrow(dtrees))]
	#
	#	read transmission window stats
	#	
	cat('\nread transmission window stats')
	ptyr	<- data.table(FILE=list.files(in.dir, pattern=postfix.trmwindowstats, full.names=TRUE))
	ptyr[, PTY_RUN:=ptyr[, as.integer(gsub('ptyr','',regmatches(FILE, regexpr('ptyr[0-9]+',FILE))))]]
	dwin		<- lapply(seq_len(nrow(ptyr)), function(i)
			{
				load(ptyr[i, FILE])
				tt[, PTY_RUN:= ptyr[i,PTY_RUN]]
				tt
			})
	dwin	<- do.call('rbind', dwin)
	setnames(dwin, 	c('PAT.1','PAT.2','PAT.1_LEAVES','PAT.2_LEAVES','PAT.1_READS','PAT.2_READS','PATHS.12','PATHS.21'),
					c('ID1','ID2','ID1_L','ID2_L','ID1_R','ID2_R','PATHS_12','PATHS_21'))
	set(dwin, NULL, 'ID1', dwin[, regmatches(ID1, regexpr(regex.ind, ID1))])
	set(dwin, NULL, 'ID2', dwin[, regmatches(ID2, regexpr(regex.ind, ID2))])
	set(dwin, NULL, 'PATRISTIC_DISTANCE', dwin[, as.numeric(PATRISTIC_DISTANCE)])	
	#
	#	build transmission summary stats based on selection criteria
	#	
	cat('\nreduce transmission window stats to windows with at least',trmw.min.reads,'reads and at least',trmw.min.tips,'tips')
	dwin	<- subset(dwin, ID1_R>=trmw.min.reads & ID2_R>=trmw.min.reads & ID1_L>=trmw.min.tips & ID2_L>=trmw.min.tips)
	cat('\ntotal number of windows with trm assignments is',nrow(dwin))		
	#
	#	merge assignments (keeping as much detail as needed for full interpretation)
	#
	dwin[, TYPE_RAW:= TYPE]
	#	chains with no intermediate
	tmp		<- dwin[, which(TYPE=="anc_12" & CONTIGUOUS)]
	cat('\nFound contiguous anc_12, n=', length(tmp),'--> chain with no intermediate')
	set(dwin, tmp, 'TYPE', 'chain_12_nointermediate')
	tmp		<- dwin[, which(TYPE=="multi_anc_12" & CONTIGUOUS)]
	cat('\nFound contiguous multi_anc_12, n=', length(tmp),'--> chain with no intermediate')
	set(dwin, tmp, 'TYPE', 'chain_12_nointermediate')
	#	chains with no intermediate
	tmp		<- dwin[, which(TYPE=="anc_21" & CONTIGUOUS)]
	cat('\nFound contiguous anc_21, n=', length(tmp),'--> chain with no intermediate')
	set(dwin, tmp, 'TYPE', 'chain_21_nointermediate')
	tmp		<- dwin[, which(TYPE=="multi_anc_21" & CONTIGUOUS)]
	cat('\nFound contiguous multi_anc_21, n=', length(tmp),'--> chain with no intermediate')
	set(dwin, tmp, 'TYPE', 'chain_21_nointermediate')	
	#
	#	chains with intermediate
	tmp		<- dwin[, which(TYPE=="anc_12" & !CONTIGUOUS)]
	cat('\nFound discontiguous anc_12, n=', length(tmp),'--> chain with intermediate')
	set(dwin, tmp, 'TYPE', 'chain_12_withintermediate')
	tmp		<- dwin[, which(TYPE=="multi_anc_12" & !CONTIGUOUS)]
	cat('\nFound discontiguous multi_anc_12, n=', length(tmp),'--> chain with intermediate')
	set(dwin, tmp, 'TYPE', 'chain_12_withintermediate')
	#	chains with intermediate
	tmp		<- dwin[, which(TYPE=="anc_21" & !CONTIGUOUS)]
	cat('\nFound discontiguous anc_21, n=', length(tmp),'--> chain with intermediate')
	set(dwin, tmp, 'TYPE', 'chain_21_withintermediate')
	tmp		<- dwin[, which(TYPE=="multi_anc_21" & !CONTIGUOUS)]
	cat('\nFound discontiguous multi_anc_21, n=', length(tmp),'--> chain with intermediate')
	set(dwin, tmp, 'TYPE', 'chain_21_withintermediate')
	#
	#	intermingled with no intermediate
	tmp		<- dwin[, which(TYPE=="conflict" & CONTIGUOUS)]
	cat('\nFound contiguous conflict, n=', length(tmp),'--> intermingled with no intermediate')
	set(dwin, tmp, 'TYPE', 'intermingled_nointermediate')
	#	intermingled with intermediate
	tmp		<- dwin[, which(TYPE=="conflict" & !CONTIGUOUS)]
	cat('\nFound not contiguous conflict, n=', length(tmp),'--> intermingled with intermediate')
	set(dwin, tmp, 'TYPE', 'intermingled_withintermediate')
	#
	#	other	
	tmp		<- dwin[, which(CONTIGUOUS & PATHS_12==0 & PATHS_21==0)]
	cat('\nFound contiguous with no paths, n=', length(tmp),'--> other')
	set(dwin, tmp, 'TYPE', 'other')
	tmp		<- dwin[, which(!CONTIGUOUS & PATHS_12==0 & PATHS_21==0)]
	cat('\nFound discontiguous with no assignment, n=', length(tmp),'--> other')
	set(dwin, tmp, 'TYPE', 'other')
	#	check
	stopifnot( !nrow(subset(dwin, TYPE=='none'))	)
	#
	#	merge assignments (short for pairs)
	#
	dwin[, TYPE_PAIR:= TYPE]
	cat('\nBuilding TYPE_PAIR')
	tmp		<- dwin[, which(grepl('nointermediate', TYPE_PAIR))]
	cat('\nFound pairs, n=', length(tmp),'--> pair')
	set(dwin, tmp, 'TYPE_PAIR', 'pair')
	tmp		<- dwin[, which(grepl('withintermediate', TYPE_PAIR))]
	cat('\nFound individuals with intermediates, n=', length(tmp),'--> withintermediate')
	set(dwin, tmp, 'TYPE_PAIR', 'withintermediate')
	tmp		<- dwin[, which(grepl('other', TYPE_PAIR))]
	cat('\nFound other relationships, n=', length(tmp),'--> other')
	set(dwin, tmp, 'TYPE_PAIR', 'other')
	#
	#	add distance as second dimension
	#
	if(!is.na(trmw.close.brl) & is.finite(trmw.close.brl))
	{
		cat('\nidentifying close pairwise assignments using distance=',trmw.close.brl)
		tmp		<- dwin[, which(PATRISTIC_DISTANCE<trmw.close.brl)]
		cat('\nFound close, n=', length(tmp))
		set(dwin, tmp, 'TYPE', dwin[tmp, paste0(TYPE,'_close')])
		set(dwin, tmp, 'TYPE_PAIR', dwin[tmp, paste0(TYPE_PAIR,'_close')])
	}
	if(!is.na(trmw.distant.brl) & is.finite(trmw.distant.brl))
	{
		cat('\nidentifying distant pairwise assignments using distance=',trmw.distant.brl)
		tmp		<- dwin[, which(PATRISTIC_DISTANCE>=trmw.distant.brl)]
		cat('\nFound distant, n=', length(tmp))
		set(dwin, tmp, 'TYPE', dwin[tmp, paste0(TYPE,'_distant')])
		set(dwin, tmp, 'TYPE_PAIR', dwin[tmp, paste0(TYPE_PAIR,'_distant')])
	}
	#	
	#	summarise transmission stats
	#	
	dtrms	<- dwin[, list(WIN_OF_TYPE=length(W_FROM), ID1_R=ID1_R[1], ID1_L=ID1_L[1], ID2_R=ID2_R[1], ID2_L=ID2_L[1]), by=c('PTY_RUN','ID1','ID2','TYPE','TYPE_PAIR')]
	tmp		<- dtrms[, list(WIN_TOTAL=sum(WIN_OF_TYPE)), by=c('PTY_RUN','ID1','ID2')]
	dtrms	<- merge(dtrms, tmp, by=c('PTY_RUN','ID1','ID2'))
	#	set pair id	
	tmp		<- dtrms[, list(SCORE=sum(WIN_OF_TYPE[grepl('close|chain_12_nointermediate|chain_21_nointermediate|intermingled_nointermediate',TYPE)])), by=c('ID1','ID2','PTY_RUN')]
	dtrms	<- merge(dtrms, tmp, by=c('ID1','ID2','PTY_RUN'))
	#	give every pair an ID
	tmp		<- unique(dtrms, by=c('ID1','ID2'))
	tmp		<- tmp[order(-SCORE),]
	tmp[, PAIR_ID:= seq_len(nrow(tmp))]	
	tmp2	<- unique(dtrms, by=c('ID1','ID2','PTY_RUN'))
	tmp		<- merge(tmp2, subset(tmp, select=c(ID1,ID2,PAIR_ID)), by=c('ID1','ID2'))
	tmp		<- tmp[, list(PAIR_ID= paste(PAIR_ID,'-',seq_along(PAIR_ID),sep=''), PTY_RUN=PTY_RUN), by=c('ID1','ID2')]	
	dtrms	<- merge(dtrms, tmp, by=c('ID1','ID2', 'PTY_RUN'))	
	setkey(dtrms, PAIR_ID)
	#	save to file
	if(!is.na(save.file))
	{
		cat('\nwrite to file', save.file)
		save(phs, dtrees, dtrms, dwin, file=save.file)		
	}
	list(phs=phs, dtrees=dtrees, dtrms=dtrms, dwin=dwin)
}	

#' @export
#' @title Generate bash command to resume a single phyloscanner run
#' @param prefix.infiles File name that points phyloscanner output.
#' @param pty.args List of phyloscanner input variables. See examples.
#' @return Character string of phyloscanner commands.
#' @description This function generates bash commands to resume a single phyloscanner run, from the point where all read alignments and read phylogenies were created. The bash script can be called via 'system' in R, or written to file to run on a UNIX system.
#' @example example/ex.cmd.phyloscanner.one.resume.R      
phsc.cmd.phyloscanner.one.resume<- function(prefix.infiles, pty.args)
{	
	#	create local tmp dir
	cmd		<- paste("CWD=$(pwd)\n",sep='\n')
	cmd		<- paste(cmd,"echo $CWD\n",sep='')
	tmpdir	<- paste('pty','_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')	
	tmpdir	<- paste("$CWD/",tmpdir,sep='')
	cmd		<- paste(cmd,'mkdir -p "',tmpdir,'"\n',sep='')
	#	copy required files to local tmp dir	
	file.bam<- list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*bam.txt',sep=''), full.names=TRUE)
	stopifnot(length(file.bam)==1)	
	cmd		<- paste(cmd,'cp "',file.bam,'" "',tmpdir,'"\n',sep='')
	if(1)	#OLD CODE (as long as we work with prev generated zip files)
	{
		tmp		<- list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*otherstuff.zip',sep=''), full.names=TRUE)	
	}	
	if(0)	#NEW CODE TODO 
	{
		tmp		<- list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*DuplicateReadCounts.zip',sep=''), full.names=TRUE)	
	}
	stopifnot(length(tmp)==1)
	cmd		<- paste(cmd,'unzip "',tmp,'" -d "',tmpdir,'"\n',sep='')
	tmp		<- list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*fasta.zip',sep=''), full.names=TRUE)
	stopifnot(length(tmp)==1)
	cmd		<- paste(cmd,'unzip "',tmp,'" -d "',tmpdir,'"\n',sep='')
	tmp		<- list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*newick.zip',sep=''), full.names=TRUE)
	stopifnot(length(tmp)==1)
	cmd		<- paste(cmd,'unzip "',tmp,'" -d "',tmpdir,'"\n',sep='')
	#	cd to tmp dir
	cmd		<- paste(cmd, 'cd "',tmpdir,'"\n', sep='')	
	#	add all toolkit commands according to pty.args
	cmd		<- paste(cmd, phsc.cmd.process.phyloscanner.output.in.directory(tmpdir, file.bam, pty.args), collapse='\n',sep='')
	#	move all files starting with current run ID
	run.id	<- gsub('_bam.txt','',basename(file.bam))
	cmd		<- paste(cmd, '\nmv ',run.id,'* "',pty.args$out.dir,'"\n',sep='')	
	#	zip up everything else
	cmd		<- paste(cmd, 'for file in *; do\n\tzip -ur9XTj ',paste(run.id,'_otherstuff.zip',sep=''),' "$file"\ndone\n',sep='')
	cmd		<- paste(cmd, 'mv ',paste(run.id,'_otherstuff.zip',sep=''),' "',pty.args$out.dir,'"\n',sep='')
	#	clean up
	cmd		<- paste(cmd,'cd $CWD\nrm -r "',tmpdir,'"\n',sep='')		
	cmd
}


#' @export
#' @title Generate bash command for a single phyloscanner run
#' @param file.bam File name of the file that contains the list of bam files.
#' @param file.ref File name of the file that contains the list of reference files.
#' @param pty.args List of phyloscanner input variables. See examples.
#' @return Character string of phyloscanner commands.
#' @description This function generates bash commands for a single phyloscanner run, that can be called via 'system' in R, or written to file to run on a UNIX system.
#' @example example/ex.cmd.phyloscanner.one.R    
phsc.cmd.phyloscanner.one<- function(file.bam, file.ref, pty.args)
{	
	stopifnot(is.character(file.bam),is.character(file.ref))
	#	copy input variables into namespace	 
	attach(pty.args)	
	#	sense checks
	stopifnot(is.logical(all.bootstrap.trees))	
	#	define window coordinates
	window.coord			<- integer(0)
	if(!nchar(pty.args[['window.automatic']]))
	{
		stopifnot(length(pty.args[['win']])==4)
		window.coord		<- seq(pty.args[['win']][1], pty.args[['win']][2]-max(pty.args[['win']][3:4]),pty.args[['win']][3])
		window.coord		<- as.vector(rbind( window.coord,window.coord-1+pty.args[['win']][4] ))											
	}	
	#	
	if(!nchar(window.automatic))	stopifnot( is.numeric(window.coord), !length(window.coord)%%2)
	if(nchar(window.automatic))		stopifnot( !length(window.coord) )
	merge.paired.reads			<- ifelse(!is.na(merge.paired.reads) & merge.paired.reads, '--merge-paired-reads', NA_character_)
	keep.overhangs				<- ifelse(!is.na(keep.overhangs) & keep.overhangs, '--keep-overhangs', NA_character_)
	no.trees					<- ifelse(!is.na(no.trees) & no.trees, '--no-trees', NA_character_)
	num.bootstraps				<- ifelse(is.na(no.trees) & !is.na(num.bootstraps) & is.numeric(num.bootstraps) & num.bootstraps>1, as.integer(num.bootstraps), NA_integer_)
	dont.check.duplicates		<- ifelse(!is.na(dont.check.duplicates) & dont.check.duplicates, '--dont-check-duplicates', NA_character_)
	#	create local tmp dir
	cmd		<- paste("CWD=$(pwd)\n",sep='\n')
	cmd		<- paste(cmd,"echo $CWD\n",sep='')
	tmpdir	<- paste('pty','_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')	
	tmpdir	<- paste("$CWD/",tmpdir,sep='')
	cmd		<- paste(cmd,'mkdir -p "',tmpdir,'"\n',sep='')
	#	copy files to local tmp dir
	cmd		<- paste(cmd,'cp "',file.bam,'" "',tmpdir,'"\n',sep='')
	cmd		<- paste(cmd,'cp "',file.ref,'" "',tmpdir,'"\n',sep='')
	#	cd to tmp dir
	cmd		<- paste(cmd, 'cd "',tmpdir,'"\n', sep='')	
	cmd		<- paste(cmd, prog.pty,' "',basename(file.bam),'" "',basename(file.ref),'" ',sep='')	
	cmd		<- paste(cmd, '--merging-threshold', merge.threshold,'--min-read-count',min.read.count,'--quality-trim-ends', quality.trim.ends, '--min-internal-quality',min.internal.quality)
	if(!is.na(merge.paired.reads))
		cmd	<- paste(cmd,' ',merge.paired.reads,sep='')	
	if(!is.na(dont.check.duplicates))
		cmd	<- paste(cmd,' ',dont.check.duplicates,sep='')
	if(!is.na(alignments.pairwise.to))		
		cmd	<- paste(cmd,' --pairwise-align-to ',alignments.pairwise.to,sep='')
	if(nchar(window.automatic))
		cmd	<- paste(cmd,' --auto-window-params ', window.automatic,sep='')
	if(!nchar(window.automatic))
		cmd	<- paste(cmd,' --windows ', paste(as.character(window.coord), collapse=','),sep='')
	if(!is.na(alignments.file))
		cmd	<- paste(cmd,' --alignment-of-other-refs ',alignments.file,sep='')	
	if(!is.na(alignments.root))
		cmd	<- paste(cmd,' --ref-for-rooting ',alignments.root,sep='')	
	if(!is.na(no.trees))		
		cmd	<- paste(cmd, no.trees)
	if(!is.na(num.bootstraps))
		cmd	<- paste(cmd, ' --num-bootstraps ', num.bootstraps, sep='')
	if(!is.na(keep.overhangs))
		cmd	<- paste(cmd, keep.overhangs)
	cmd		<- paste(cmd, '--x-raxml',prog.raxml,'--x-mafft',prog.mafft,'\n')
	run.id	<- gsub('_bam.txt','',basename(file.bam))
	#	process RAxML files
	if(is.na(no.trees) & (is.na(num.bootstraps) | (!is.na(num.bootstraps) & all.bootstrap.trees)))
		cmd	<- paste(cmd, 'for file in RAxML_bestTree*.tree; do\n\tmv "$file" "${file//RAxML_bestTree\\./',run.id,'_}"\ndone\n',sep='')
	if(is.na(no.trees) & !is.na(num.bootstraps) & !all.bootstrap.trees)
		cmd	<- paste(cmd, 'for file in RAxML_bipartitions.MLtreeWbootstraps*.tree; do\n\tmv "$file" "${file//RAxML_bipartitions.MLtreeWbootstraps/',run.id,'_}"\ndone\n',sep='')	
	#cmd	<- paste(cmd, "for file in AlignedReads*.fasta; do\n\tsed 's/<unknown description>//' \"$file\" > \"$file\".sed\n\tmv \"$file\".sed \"$file\"\ndone\n",sep='')		
	if(!is.na(alignments.file) & !is.na(keep.overhangs))
	{
		cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\tcat "$file" | awk \'{if (substr($0,1,4) == ">REF") censor=1; else if (substr($0,1,1) == ">") censor=0; if (censor==0) print $0}\' > NoRef$file\ndone\n', sep='')
		#cmd	<- paste(cmd, 'for file in NoRefAlignedReads*.fasta; do\n\tcp "$file" "${file//NoRefAlignedReads/NoRef',run.id,'_}"\n\tmv NoRef', run.id,'*fasta "',out.dir,'"\ndone\n',sep='')
		cmd	<- paste(cmd, 'for file in NoRefAlignedReads*.fasta; do\n\t',phsc.cmd.mafft.add(alignments.file,'"$file"','Ref"$file"', options='--keeplength --memsave --parttree --retree 1'),'\ndone\n',sep='')		
		cmd	<- paste(cmd, 'for file in RefNoRefAlignedReads*.fasta; do\n\t','mv "$file" "${file//RefNoRefAlignedReads/',run.id,'_}"\ndone\n',sep='')		
	}
	if(is.na(alignments.file) || is.na(keep.overhangs))
	{
		cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\tmv "$file" "${file//AlignedReads/',run.id,'_}"\ndone\n',sep='')	
	}
	#	move Duplicate Read Counts
	cmd		<- paste(cmd, 'for file in DuplicateReadCountsProcessed_*.csv; do\n\tmv "$file" "${file//DuplicateReadCountsProcessed_/',run.id,'_DuplicateReadCounts_}"\ndone',sep='')
	#	run phyloscanner tools and compress output
	cmd		<- paste(cmd, phsc.cmd.process.phyloscanner.output.in.directory(tmpdir, file.bam, pty.args), sep='\n')
	#
	cmd		<- paste(cmd, '\nmv ',run.id,'* "',out.dir,'"\n',sep='')	
	#	zip up everything else
	cmd		<- paste(cmd, 'for file in *; do\n\tzip -ur9XTjq ',paste(run.id,'_otherstuff.zip',sep=''),' "$file"\ndone\n',sep='')
	cmd		<- paste(cmd, 'mv ',paste(run.id,'_otherstuff.zip',sep=''),' "',out.dir,'"\n',sep='')
	#	clean up
	cmd		<- paste(cmd,'cd $CWD\nrm -r "',tmpdir,'"\n',sep='')
	detach(pty.args)
	cmd
}

#' @export
#' @title Generate bash commands for a multiple phyloscanner runs
#' @param pty.runs Data.table of individual assignments to phyloscanner runs, with columns 'PTY_RUN' (run id) and 'IND_ID' (ID of individuals that are assigned to that run).
#' @param pty.args List of phyloscanner input variables. See examples.
#' @return Data.table with columns 'PTY_RUN' (run id) and 'CMD' (bash commands for that run). 
#' @description This function generates bash commands for multiple phyloscanner runs, that can be called via 'system' in R, or written to file to run on a UNIX system.
#' @example example/ex.cmd.phyloscanner.multi.R  
phsc.cmd.phyloscanner.multi <- function(pty.runs, pty.args) 		
{
	#
	#	associate BAM and REF files with each scheduled phylotype run
	#	
	#	get available Bam files
	ptyd		<- data.table(FILE=list.files(pty.args[['data.dir']], full.names=TRUE))
	ptyd[, TYPE:=NA_character_]
	set(ptyd, ptyd[, which(grepl('_ref.fasta$',FILE))], 'TYPE', 'REF')
	set(ptyd, ptyd[, which(grepl('.bam$',FILE))], 'TYPE', 'BAM')
	ptyd		<- subset(ptyd, !is.na(TYPE))
	ptyd[, IND_ID:= gsub('\\.bam|_ref\\.fasta','',basename(FILE))]
	ptyd		<- dcast.data.table(ptyd, IND_ID~TYPE, value.var='FILE')
	#	merge
	ptyd	<- merge(pty.runs, ptyd, by='IND_ID', all.x=1)
	if(!any(is.na(pty.args[['select']])))
		ptyd<- subset(ptyd, PTY_RUN%in%pty.args[['select']])		
	if(ptyd[,any(is.na(BAM))])
		warning('\nCould not find location of BAM files for all individuals in pty.runs, n=', ptyd[, length(which(is.na(BAM)))],'\nMissing individuals are ignored. Please check.')	
	if(ptyd[,any(is.na(REF))])
		warning('\nCould not find location of reference files for all individuals in pty.runs, n=', ptyd[, length(which(is.na(REF)))],'\nMissing individuals are ignored. Please check.')	
	setkey(ptyd, PTY_RUN)
	#
	#	write pty.run files and get pty command lines
	#
	pty.c		<- ptyd[, {
				#PTY_RUN<- z <- 5; BAM<- subset(ptyd, PTY_RUN==z)[, BAM]; REF<- subset(ptyd, PTY_RUN==z)[, REF]									
				file.bam	<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_bam.txt',sep=''))
				file.ref	<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_ref.txt',sep=''))
				cat( paste(BAM[!is.na(BAM)],collapse='\n'), file= file.bam	)
				cat( paste(REF[!is.na(REF)],collapse='\n'), file= file.ref	)
				cmd			<- phsc.cmd.phyloscanner.one(file.bam, file.ref, pty.args)
				#cmd			<- paste(cmd, pty.cmd.evaluate.fasta(pty.args[['out.dir']], strip.max.len=pty.args[['strip.max.len']], select=paste('^ptyr',PTY_RUN,'_In',sep=''), min.ureads.individual=pty.args[['min.ureads.individual']]), sep='')
				#cat(cmd)
				list(CMD= cmd)				
			},by='PTY_RUN']
	pty.c
}	

#' @import data.table
#' @title Generate bash commands to process phyloscanner output
#' @description This function generates bash commands that combine the various Rscripts in the phyloscanner toolkit   
#' @param tmp.dir Directory with phyloscanner output.
#' @param file.bam File name of the file that contains the list of bam files.
#' @param pty.args List of phyloscanner input variables. 
#' @return character string of bash commands.
phsc.cmd.process.phyloscanner.output.in.directory<- function(tmp.dir, file.bam, pty.args)
{
	#run.id		<- 'ptyr5'; tmp.dir		<- '$CWD/pty_16-09-08-07-32-26'; file.bam	<- '/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptinput/ptyr5_bam.txt'
			
	#	define variables
	prog.pty						<- pty.args[['prog.pty']]
	root.name						<- pty.args[['alignments.root']]	
	sankhoff.k						<- pty.args[['sankhoff.k']]
	use.blacklisters				<- pty.args[['use.blacklisters']]
	contaminant.read.threshold		<- pty.args[['contaminant.read.threshold']]
	contaminant.prop.threshold		<- pty.args[['contaminant.prop.threshold']]
	dual.prop.threshold				<- pty.args[['dual.prop.threshold']]
	roguesubtree.prop.threshold		<- pty.args[['roguesubtree.prop.threshold']]
	roguesubtree.read.threshold		<- pty.args[['roguesubtree.read.threshold']]
	rogue.prop.threshold			<- pty.args[['rogue.prop.threshold']]
	rogue.longestBranchLength		<- pty.args[['rogue.longestBranchLength']]
	rogue.probThreshold				<- pty.args[['rogue.probThreshold']]	
	dwns.maxReadsPerPatient			<- pty.args[['dwns.maxReadsPerPatient']]	
	pty.tools.dir					<- file.path(dirname(prog.pty),'tools')
	prog.pty.readblacklist			<- paste('Rscript ',file.path(pty.tools.dir,'MakeReadBlacklist.R'),sep='')
	prog.pty.readblacklistsankoff	<- paste('Rscript ',file.path(pty.tools.dir,'ParsimonyBasedBlacklister.R'),sep='')
	prog.pty.rogueblacklist			<- paste('Rscript ',file.path(pty.tools.dir,'MakeRogueBlacklist.R'),sep='')
	prog.pty.roguewblacklist		<- paste('Rscript ',file.path(pty.tools.dir,'MakeRogueBlacklistWeibull.R'),sep='')
	prog.pty.dualblacklist			<- paste('Rscript ',file.path(pty.tools.dir,'DualPatientBlacklister.R'),sep='')
	prog.pty.downsample				<- paste('Rscript ',file.path(pty.tools.dir,'DownsampleReads.R'),sep='')
	prog.pty.split					<- paste('Rscript ',file.path(pty.tools.dir,'SplitPatientsToSubtrees.R'),sep='')
	prog.pty.smry					<- paste('Rscript ',file.path(pty.tools.dir,'SummaryStatistics.R'),sep='')
	prog.pty.lkltrm					<- paste('Rscript ',file.path(pty.tools.dir,'ClassifyRelationshipsExperimental.R'),sep='')	
	prog.pty.lkl.smry				<- paste('Rscript ',file.path(pty.tools.dir,'TransmissionSummaryExperimental.R'),sep='')	
	run.id							<- gsub('_bam.txt','',basename(file.bam))
	run.id_							<- ifelse(grepl('[a-z0-9]',substring(run.id, nchar(run.id))), paste(run.id,'_',sep=''), run.id)
	blacklistFiles					<- NA_character_
	#
	stopifnot(	!is.null(use.blacklisters) & !is.na(use.blacklisters)	)
	if(any(grepl('MakeReadBlacklist', use.blacklisters)))
		stopifnot(	!is.null(contaminant.prop.threshold) & !is.na(contaminant.prop.threshold), 
					!is.null(contaminant.read.threshold) & !is.na(contaminant.read.threshold))
	if(any(grepl('ParsimonyBasedBlacklister', use.blacklisters)))
		stopifnot(	!any(grepl('MakeRogueBlacklist', use.blacklisters)),	
					!is.null(roguesubtree.prop.threshold) & !is.na(roguesubtree.prop.threshold),
					!is.null(roguesubtree.read.threshold) & !is.na(roguesubtree.read.threshold),
					!is.null(sankhoff.k) & !is.na(sankhoff.k))		
	if(any(grepl('MakeRogueBlacklist', use.blacklisters))) 
		stopifnot(	!any(grepl('MakeRogueBlacklistWeibull', use.blacklisters)),
					!is.null(rogue.prop.threshold) & !is.na(rogue.prop.threshold),
					!is.null(rogue.longestBranchLength) & !is.na(rogue.longestBranchLength))
	if(any(grepl('MakeRogueBlacklistWeibull', use.blacklisters))) 
		stopifnot(	!any(grepl('ParsimonyBasedBlacklister', use.blacklisters)),
					!is.null(rogue.prop.threshold) & !is.na(rogue.prop.threshold),
					!is.null(rogue.longestBranchLength) & !is.na(rogue.longestBranchLength),
					!is.null(rogue.probThreshold) & !is.na(rogue.probThreshold))
	if(any(grepl('DualPatientBlacklister', use.blacklisters))) 
		stopifnot(!is.null(dual.prop.threshold) & !is.na(dual.prop.threshold))
	if(any(grepl('DownsampleReads', use.blacklisters))) 
		stopifnot(!is.null(dwns.maxReadsPerPatient) & !is.na(dwns.maxReadsPerPatient))
	use.sankhoff.blacklister		<- any(grepl('ParsimonyBasedBlacklister', use.blacklisters))
	#
	cmd					<- ''
	#
	if(1)	# OLD CODE (as long as we work with prev generated zip files) TODO: swith off
	{
		cmd				<- paste(cmd, 'for file in DuplicateReadCountsProcessed_*.csv; do\n\tmv "$file" "${file//DuplicateReadCountsProcessed_/',run.id_,'DuplicateReadCounts_}"\ndone',sep='')	
	}	
	#
	#	bash command to blacklist taxa with duplicate counts that suggest contaminants
	#	
	if(any(grepl('MakeReadBlacklist', use.blacklisters)))
	{
		cmd				<- paste(cmd, '\nfor file in ', run.id_,'DuplicateReadCounts_*.csv; do\n\t',sep='')
		tmp				<- phsc.cmd.blacklist.reads(prog.pty.readblacklist, '"$file"', paste('"${file//DuplicateReadCounts/blacklist}"',sep=''), contaminant.read.threshold, contaminant.prop.threshold, tipRegex=NA)
		cmd				<- paste(cmd, tmp, '\ndone', sep='')	
		blacklistFiles	<- file.path(tmp.dir, paste(run.id_,'blacklist_',sep=''))
	}
	#
	#	bash command to blacklist taxa with duplicate counts that suggest contaminants
	#		and identifying contaminants through a Sankhoff parsimony reconstruction
	#
	if(any(grepl('ParsimonyBasedBlacklister', use.blacklisters)))	
	{
		cmd				<- paste(cmd, '\nfor file in ', run.id_,'InWindow_*.tree; do\n\t',sep='')
		cmd				<- paste(cmd,'TMP=${file//tree/csv}\n\t',sep='')
		tmp				<- phsc.cmd.blacklist.parsimonybased( 	prog.pty.readblacklistsankoff, 
																	pty.tools.dir,
																	'"$file"', 
																	paste('"${TMP//InWindow/blacklistsank_InWindow}"',sep=''), 
																	dualCandidatesOutputFileName=paste('"${TMP//InWindow/duallistsank_InWindow}"',sep=''),
																	rawThreshold=roguesubtree.read.threshold, 
																	ratioThreshold=roguesubtree.prop.threshold, 
																	sankhoffK=sankhoff.k,	
																	outgroupName=root.name,
																	tipRegex=NA)
		cmd				<- paste(cmd, tmp, '\ndone', sep='')	
		blacklistFiles	<- file.path(tmp.dir, paste(run.id_,'blacklistsank_',sep=''))
	}		
	#
	#	bash command to make blacklists of rogue taxa for each window based on branch lengths
	#
	if(any(grepl('MakeRogueBlacklist', use.blacklisters)))	
	{
		cmd				<- paste(cmd, '\nfor file in ', file.path(tmp.dir,run.id_),'*tree; do\n\t',sep='')
		cmd				<- paste(cmd,'TMP=${file//tree/csv}\n\t',sep='')
		tmp				<- ifelse(is.na(blacklistFiles), NA_character_, '"${TMP//InWindow/blacklist_InWindow}"')
		tmp				<- phsc.cmd.blacklist.rogue.geneticdistance(	prog.pty.rogueblacklist, 
															pty.tools.dir, 
															'"$file"', 
															'"${TMP//InWindow/blacklistrogue_InWindow}"', 
															longestBranchLength=rogue.longestBranchLength, 
															dropProportion=rogue.prop.threshold, 
															blackListFileName=tmp, 
															outgroupName=root.name, 
															tipRegex=NA)				
		cmd				<- paste(cmd, tmp, '\n\t','mv "${TMP//InWindow/blacklistrogue_InWindow}" "${TMP//InWindow/blacklist_InWindow}"','\n','done', sep='')
		blacklistFiles	<- file.path(tmp.dir, paste(run.id_,'blacklist_',sep=''))
	}
	#
	#	bash command to make blacklists of rogue taxa for each window based on Weibull extreme value probability of branch lengths
	#
	if(any(grepl('MakeRogueBlacklistWeibull', use.blacklisters)))	
	{
		cmd				<- paste(cmd, '\nfor file in ', file.path(tmp.dir,run.id_),'*tree; do\n\t',sep='')
		cmd				<- paste(cmd,'TMP=${file//tree/csv}\n\t',sep='')
		tmp				<- ifelse(is.na(blacklistFiles), NA_character_, '"${TMP//InWindow/blacklist_InWindow}"')
		tmp				<- phsc.cmd.blacklist.rogue.extremeprob(	prog.pty.roguewblacklist, 
															pty.tools.dir, 
															'"$file"', 
															'"${TMP//InWindow/blacklistrogue_InWindow}"', 
															probThreshold=rogue.probThreshold,
															longestBranchLength=rogue.longestBranchLength, 
															dropProportion=rogue.prop.threshold, 
															blackListFileName=tmp, 
															outgroupName=root.name, 
															tipRegex=NA)				
		cmd				<- paste(cmd, tmp, '\n\t','mv "${TMP//InWindow/blacklistrogue_InWindow}" "${TMP//InWindow/blacklist_InWindow}"','\n','done', sep='')
		blacklistFiles	<- file.path(tmp.dir, paste(run.id_,'blacklist_',sep=''))
	}
	#
	#	bash command to make blacklists of duplicate taxa based on candidate duplicates output from ParsimonyBlacklist.R
	#
	if(any(grepl('DualPatientBlacklister', use.blacklisters)))	
	{			
		tmp				<- phsc.cmd.blacklist.dualinfections(	prog.pty.dualblacklist, 																 
																file.path(tmp.dir,paste0(run.id_,'duallistsank_')),
																file.path(tmp.dir,paste0(run.id_,'blacklistdual_')),
																treeFileName=file.path(tmp.dir,paste0(run.id_,'.*tree')),
																summaryFileName=file.path(tmp.dir,paste0(run.id_,'dualsummary.csv')),
																blacklistFileName=blacklistFiles, 
																dual.prop.threshold=dual.prop.threshold)
		cmd				<- paste(cmd, tmp, sep='\n')
		cmd				<- paste0(cmd, '\n','for file in ', basename(blacklistFiles),'*csv; do\n\t','cp "$file" "${file//',basename(blacklistFiles),'/',paste0(run.id_,'blacklistfinal_'),'}"','\n','done')
		cmd				<- paste0(cmd, '\n','for file in ', paste0(run.id_,'blacklistdual_'),'*csv; do\n\t','cp "$file" "${file//blacklistdual/blacklistfinal}"','\n','done')
		blacklistFiles	<- file.path(tmp.dir, paste(run.id_,'blacklistfinal_',sep=''))
	}
	#
	#	bash command to downsample tips (add to blacklist)
	#
	if(any(grepl('DownsampleReads', use.blacklisters)))
	{
		tmp				<- phsc.cmd.blacklist.downsample(		prog.pty.downsample, 
																pty.tools.dir, 
																file.path(tmp.dir,run.id_),																
																file.path(tmp.dir,paste(run.id_,'blacklistdwns_',sep='')), 
																maxReadsPerPatient=dwns.maxReadsPerPatient, 
																blackListFileName=blacklistFiles, 
																outgroupName=root.name, 
																tipRegex=NA)
		cmd				<- paste(cmd, tmp, sep='\n')
		blacklistFiles	<- file.path(tmp.dir, paste(run.id_,'blacklistdwns_',sep=''))
		#cmd				<- paste(cmd, '\n','for file in ', file.path(tmp.dir,paste(run.id_,'blacklistdwns_*',sep='')),'; do\n\t',sep='')
		#cmd				<- paste(cmd, 'mv "$file" "${file//blacklistdwns/blacklist}"\ndone',sep='')		
	}
	#
	#	bash command to plot trees and calculate splits
	#							
	tmp				<- phsc.cmd.SplitPatientsToSubtrees(	prog.pty.split,
															pty.tools.dir,
															file.path(tmp.dir,run.id_),
															blacklistFiles=blacklistFiles,
															outputdir=tmp.dir,																													
															outgroupName=root.name,
															splitsRule= ifelse(use.sankhoff.blacklister, 's', 'r'),
															sankhoff.k=sankhoff.k,
															pdfwidth=30, pdfrelheight=0.15)
	cmd				<- paste(cmd, tmp, sep='\n')
	file.patients	<- paste(run.id_,'patients.txt',sep='')	
	cmd				<- paste(cmd,"\nsed 's/.*\\///' \"", file.path(tmp.dir,basename(file.bam)), '" > "',file.path(tmp.dir,file.patients),'"', sep='')	
	#
	#	bash command to calculate patient stats
	#							
	tmp				<- phsc.cmd.SummaryStatistics( 	prog.pty.smry, 
													pty.tools.dir, 
													root.name, 
													file.path(tmp.dir, file.patients), 
													file.path(tmp.dir, paste(run.id_,'InWindow_',sep='')), 
													file.path(tmp.dir, paste('Subtrees_',ifelse(use.sankhoff.blacklister, 's', 'r'),'_',run.id_,'InWindow_',sep='')), 
													file.path(tmp.dir, substr(run.id_,1,nchar(run.id_)-1)),
													blacklistFiles=blacklistFiles
													)
	cmd				<- paste(cmd, tmp, sep='\n')
	#
	#	bash command to get likely transmissions 
	#
	tmp				<- phsc.cmd.LikelyTransmissions(	prog.pty.lkltrm, 
														pty.tools.dir, 
														file.path(tmp.dir,paste('ProcessedTree_',ifelse(use.sankhoff.blacklister, 's', 'r'),'_',run.id_,sep='')), 
														file.path(tmp.dir,paste('Subtrees_',ifelse(use.sankhoff.blacklister, 's', 'r'),'_',run.id_,sep='')), 
														file.path(tmp.dir,run.id_))
	cmd				<- paste(cmd, tmp, sep='\n')
	#
	#	add bash command to get likely transmissions summary
	#						
	tmp				<- phsc.cmd.LikelyTransmissionsSummary(	prog.pty.lkl.smry, 
															pty.tools.dir,
															file.path(tmp.dir, file.patients),												
															file.path(tmp.dir, paste(run.id_,'patStatsFull.csv',sep='')),
															file.path(tmp.dir, paste('ProcessedTree_',ifelse(use.sankhoff.blacklister, 's', 'r'),'_',run.id_,sep='')), 
															file.path(tmp.dir, paste(run.id_,'trmStats.csv',sep='')),
															file.path(tmp.dir, paste(run.id_,'trmStatsPerWindow.rda',sep='')),
															min.threshold=1, 
															allow.splits=TRUE)
	cmd				<- paste(cmd, tmp, sep='\n')
	#
	#	add bash command to compress phyloscanner output
	#							
	tmp				<- phsc.cmd.read.processed.phyloscanner.output.in.directory(file.path(tmp.dir, run.id_), 
																				file.path(tmp.dir, run.id_), 
																				read.likelytransmissions=TRUE, 
																				read.trees=TRUE, 
																				read.subtrees=TRUE, 
																				resume=FALSE, 
																				zip=TRUE)
	cmd				<- paste(cmd, tmp, sep='\n')
	cmd
}

#' @title Read processed phyloscanner output
#' @description This function creates R data.tables from processed phyloscanner output, for further data analysis in R.
#' @param prefix.infiles Full path name to processed phyloscanner output
#' @param save.file.base Output will be stored to files that start with 'save.file.base'.
#' @param read.likelytransmissions If TRUE, read and process likely transmissions 
#' @param read.trees If TRUE, read and process trees
#' @param read.subtrees If TRUE, read and process subtree files
#' @param resume If TRUE, the function does not process existing rda files.
#' @param zip If TRUE, the function zips processed phyloscanner output, and then deletes the zipped, processed phyloscanner output files.     
#' @return Nothing, rda objects written to file.
phsc.read.processed.phyloscanner.output.in.directory<- function(prefix.infiles, save.file.base, read.likelytransmissions=TRUE, read.trees=TRUE, read.subtrees=TRUE, resume=FALSE, zip=FALSE)
{
	#
	#	read all likely transmissions in indir (these are files ending in _trmStats.csv)
	#
	if(read.likelytransmissions)
	{
		save.file		<- paste(save.file.base,'trmStats.rda',sep='')
		plot.file		<- gsub('\\.rda','\\.pdf',save.file)
		phsc.read.likelytransmissions(prefix.infiles, prefix.run='ptyr', regexpr.lklsu='trmStats.csv$', regexpr.patient='^[0-9]+_[0-9]+_[0-9]+', save.file=save.file, plot.file=plot.file, resume=resume, zip=zip)				
	}
	#
	#	read trees
	#
	if(read.trees)
	{
		save.file		<- paste(save.file.base,'trees.rda',sep='')
		tmp				<- phsc.read.trees(prefix.infiles, prefix.run='ptyr', regexpr.trees='Subtrees_[crs]_.*\\.rda$', prefix.wfrom='Window_', prefix.wto='Window_[0-9]+_to_', save.file=save.file, resume=resume, zip=zip)
		tmp				<- NULL		
	}
	#
	#	read subtrees files and save (must be run after phsc.read.trees!)
	#
	if(read.subtrees)
	{
		save.file		<- paste(save.file.base,'subtrees_r.rda',sep='')
		stat.subtrees	<- phsc.read.subtrees(prefix.infiles, prefix.run='ptyr', regexpr.subtrees='Subtrees_[crs]_.*\\.rda$', prefix.wfrom='Window_', prefix.wto='Window_[0-9]+_to_', save.file=save.file, resume=resume, zip=zip)
		stat.subtrees	<- NULL				
	}
	NULL
}

#' @export
#' @import data.table grid ggtree
#' @title Plot all short read phylogenies including two individuals
#' @description This function plots short read phylogenies and highlights the clades of two individuals in red and blue.  
#' @param phs List of trees in ape format
#' @param dfs data.table with mandatory column 'IDX' and optional column 'TITLE'. IDX is the index of all phylogenies in 'phs' that are to be plotted. TITLE is a title for each sub-plot, for example specifying a window.
#' @param id1	Regular expression that identifies the first individual, to be plotted in red.  
#' @param id2	Regular expression that identifies the first individual, to be plotted in blue.
#' @param pdf.h	Height of the pdf file in inches.
#' @param pdf.rw Relative width of the pdf file, internally multiplied by the number of phylogenies to give the total width in inches.
#' @param pdf.ntrees Number of trees per pdf.
#' @param pdf.title.size Size of pdf title in inches.
#' @param plot.file If not missing, the phylogenies will be printed to file.	
#' @return List of ggtree objects, ready for printing.
phsc.plot.selected.pairs<- function(phs, dfs, id1, id2, plot.file=NA, pdf.h=50, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40)
{
	suppressMessages(require(grid))
	#	determine which phylogenies contain both individuals
	tmp		<- copy(dfs)
	tmp		<- merge(tmp, tmp[, {
				ph	<- phs[[ IDX ]]
				list(HAS_BOTH_IND= any(grepl(id1, attr(ph, "INDIVIDUAL"))) & any(grepl(id2, attr(ph, "INDIVIDUAL"))))  				
			}, by='IDX'], by='IDX')
	tmp		<- subset(tmp, HAS_BOTH_IND)
	phps	<- lapply(seq_len(nrow(tmp)), function(i){
				ph.title	<- NULL
				if('TITLE'%in%colnames(tmp))
					ph.title	<- tmp[i, TITLE]										
				ph			<- phs[[ tmp[i, IDX] ]]
				col			<- rep('grey50', length(attr(ph, "INDIVIDUAL"))) 
				col[ grepl(id1, attr(ph, "INDIVIDUAL")) ]	<- 'red'
				col[ grepl(id2, attr(ph, "INDIVIDUAL")) ]	<- 'blue'
				attr(ph, 'COLOUR')	<- col								
				#						
				p 			<- ggtree(ph, aes(color=I(COLOUR))) +
						geom_point2(shape = 16, size=3, aes(subset=NODE_SHAPES)) +
						scale_fill_hue(na.value="black") +								
						theme(legend.position="none") +
						geom_tiplab(aes(col=I(COLOUR))) +
						theme_tree2() +
						theme(legend.position="bottom", plot.title = element_text(size=pdf.title.size)) + 
						ggplot2::xlim(0, max(node.depth.edgelength(ph)[1:Ntip(ph)])*1.3) +
						labs(x='subst/site', title=ph.title)						
				p
			})
	#
	#	single page plot
	#		
	if(!is.na(plot.file))					
	{
		if(length(phps)<=pdf.ntrees)
		{
			cat('Plotting to file', plot.file,'...\n')
			pdf(file=plot.file, w=pdf.rw*length(phps), h=pdf.h)
			grid.newpage()
			pushViewport(viewport(layout=grid.layout(1, length(phps))))
			for(i in seq_along(phps))
				print(phps[[i]], vp = viewport(layout.pos.row=1, layout.pos.col=i))
			dev.off()
		}
		if(length(phps)>pdf.ntrees)
		{
			pi	<- data.table(IDX=seq_along(phps))
			pi[, PLOT:= ceiling(IDX/pdf.ntrees)]
			pi[, PLOT_IDX:= (IDX-1)%%pdf.ntrees+1]
			pi[,{
						cat('Plotting to file', gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file),'...\n')
						pdf(file=gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file), w=pdf.rw*pdf.ntrees, h=pdf.h)
						grid.newpage()
						pushViewport(viewport(layout=grid.layout(1, pdf.ntrees)))
						for(i in seq_along(IDX))
							print(phps[[IDX[i]]], vp = viewport(layout.pos.row=1, layout.pos.col=PLOT_IDX[i]))
						dev.off()
					}, by='PLOT']
		}
	}
	phps	
}

#' @export
#' @import data.table grid ggtree
#' @title Plot short read phylogenies and highlight individuals
#' @description This function plots short read phylogenies and highlights the clades of two individuals in red and blue.  
#' @param phs List of trees in ape format
#' @param dfs data.table with mandatory column 'IDX' and optional column 'TITLE'. IDX is the index of all phylogenies in 'phs' that are to be plotted. TITLE is a title for each sub-plot, for example specifying a window.
#' @param ids Vector of regular expressions that identify individuals to be highlighted in colour.
#' @param plot.cols Vector of colours for each individual  
#' @param pdf.h	Height of the pdf file in inches.
#' @param pdf.rw Relative width of the pdf file, internally multiplied by the number of phylogenies to give the total width in inches.
#' @param pdf.ntrees Number of trees per pdf.
#' @param pdf.title.size Size of pdf title in inches.
#' @param plot.file If not missing, the phylogenies will be printed to file.	
#' @return List of ggtree objects, ready for printing.
phsc.plot.selected.individuals<- function(phs, dfs, ids, plot.cols=rainbow(length(ids)), plot.file=NA, group.redo=FALSE, pdf.h=50, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40)
{
	suppressMessages(require(grid))
	#	determine which phylogenies contain at least one of the requested individuals
	tmp		<- copy(dfs)
	tmp		<- merge(tmp, tmp[, {
						ph	<- phs[[ IDX ]]
						list(HAS_ANY_IND=any(sapply(ids, function(id) any(grepl(id, ph$tip.label)) )))						  				
					}, by='IDX'], by='IDX')
	tmp		<- subset(tmp, HAS_ANY_IND)
	phps	<- lapply(seq_len(nrow(tmp)), function(i){
				ph.title	<- NULL
				if('TITLE'%in%colnames(tmp))
					ph.title	<- tmp[i, TITLE]										
				ph			<- phs[[ tmp[i, IDX] ]]
				#
				#
				if(group.redo)
				{
					phb			<- data.table(IDX=seq_along(ph$tip.label), TAXA=ph$tip.label, ID='none')
					for(id in ids)
						set(phb, phb[, which(grepl(id,TAXA))], 'ID', id)
					tmp				<- lapply( phb[, unique(ID)], function(x)	subset(phb, ID==x)[, TAXA]	)
					names(tmp)		<- phb[, unique(ID)]
					ph				<- groupOTU(ph, tmp, group='INDIVIDUAL')
					#z	<- merge(data.table(FROM=ph$edge[,1],IDX=ph$edge[,2]), phb, by='IDX', all=1)
					#z[, GROUP:= attr(ph,'INDIVIDUAL')[1:nrow(ph$edge)]]
					#z	<- unique(subset(z, !is.na(ID), select=c(ID, GROUP)))
					#attr(ph,'INDIVIDUAL')	<- factor(attr(ph,'INDIVIDUAL'), levels=c(0,z[,as.character(GROUP)]), labels=c('not characterized',z[,FILE_ID]))
					#ph									
				}				
				#
				#
				cols		<- rep('grey50', length(attr(ph, "INDIVIDUAL"))) 
				for(i in seq_along(ids))
					cols[ grepl(ids[i], attr(ph, "INDIVIDUAL")) ]<- plot.cols[i]
				attr(ph, 'COLOUR')	<- cols								
				#						
				p 			<- ggtree(ph, aes(color=I(COLOUR))) +
						geom_point2(shape = 16, size=3, aes(subset=NODE_SHAPES)) +
						scale_fill_hue(na.value="black") +								
						theme(legend.position="none") +
						geom_tiplab(aes(col=I(COLOUR))) +
						theme_tree2() +
						theme(legend.position="bottom", plot.title = element_text(size=pdf.title.size)) + 
						ggplot2::xlim(0, max(node.depth.edgelength(ph)[1:Ntip(ph)])*1.3) +
						labs(x='subst/site', title=ph.title)
				#ggsave(plot.file, w=10, h=200, limitsize = FALSE)
				p
			})
	#
	#	single page plot
	#		
	if(!is.na(plot.file))					
	{
		if(length(phps)<=pdf.ntrees)
		{
			cat('Plotting to file', plot.file,'...\n')
			pdf(file=plot.file, w=pdf.rw*length(phps), h=pdf.h)
			grid.newpage()
			pushViewport(viewport(layout=grid.layout(1, length(phps))))
			for(i in seq_along(phps))
				print(phps[[i]], vp = viewport(layout.pos.row=1, layout.pos.col=i))
			dev.off()
		}
		if(length(phps)>pdf.ntrees)
		{
			pi	<- data.table(IDX=seq_along(phps))
			pi[, PLOT:= ceiling(IDX/pdf.ntrees)]
			pi[, PLOT_IDX:= (IDX-1)%%pdf.ntrees+1]
			pi[,{
						cat('Plotting to file', gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file),'...\n')
						pdf(file=gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file), w=pdf.rw*pdf.ntrees, h=pdf.h)
						grid.newpage()
						pushViewport(viewport(layout=grid.layout(1, pdf.ntrees)))
						for(i in seq_along(IDX))
							print(phps[[IDX[i]]], vp = viewport(layout.pos.row=1, layout.pos.col=PLOT_IDX[i]))
						dev.off()
					}, by='PLOT']
		}
	}
	phps	
}

#' @import data.table ggplot2
#' @title Read likely transmissions summary files into a data.table
#' @description This function reads likely transmissions summary files from the phyloscanner toolkit. 
#' @param prefix.infiles Full path name identifying likely transmissions summary files
#' @param prefix.run Character string to identify separate phyloscanner runs. After the prefix, an integer number is expected. For example, if prefix.run='ptyr', then files are expected to start like 'ptyr123_'.
#' @param regexpr.lklsu Regular expression that identifies likely transmissions summary files in the directory.  
#' @param regexpr.patient Regular expression that identifies individuals of interest that are reported in the likely transmissions summary files.	
#' @param save.file If not missing, function output (a data.table) will be stored to 'save.file', input files will be zipped, and input files will be deleted.   
#' @param plot.file If not missing, types of evidence of transmission are plotted to this file	
#' @param resume If TRUE and save.file is not missing, the function loads and returns trees stored in save.file.
#' @return Data table with columns PAIR_ID, ID1, ID2, TYPE, WIN_OF_TYPE, PTY_RUN, WIN_TOTAL, SCORE 
phsc.read.likelytransmissions<- function(prefix.infiles, prefix.run='ptyr', regexpr.lklsu='trmStats.csv$', regexpr.patient='^[0-9]+_[0-9]+_[0-9]+', save.file=NA, plot.file=NA, resume=FALSE, zip=FALSE)
{
	if(!is.na(save.file) & resume)
	{
		options(show.error.messages = FALSE)		
		tmp		<- try(suppressWarnings(load(save.file)))
		options(show.error.messages = TRUE)
		if(!inherits(tmp, "try-error"))
			return(df)
	}		
	cat('\nReading files starting', prefix.infiles,'...\n')
	dfr	<- data.table(FILE_TRSU= list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles), '.*', regexpr.lklsu, sep=''), full.names=TRUE))
	dfr[, PTY_RUN:= as.integer(gsub(prefix.run,'',regmatches(FILE_TRSU,regexpr(paste(prefix.run,'[0-9]+',sep=''), FILE_TRSU))))]
	setkey(dfr, PTY_RUN)
	cat('\nFound likely transmission summary files, n=', nrow(dfr),'...\n')	
	if(zip & !is.na(save.file))
	{		
		tmp	<- copy(dfr)
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE_TRSU, flags = "-ur9XTjq")), by='FILE_TRSU'] )
			invisible( file.remove( dfr[, FILE_TRSU] ) )
		}
		tmp	<- data.table(FILE= list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*_LikelyTransmissions.csv$',sep=''), full.names=TRUE))
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'_LikelyTransmissions.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE, flags = "-ur9XTjq")), by='FILE'] )
			invisible( file.remove( tmp[, FILE] ) )			
		}
	}	
	#	could still do the plotting using the per window assignments
	if(0 & !is.na(plot.file))
	{		
		cat('\nPlot to file', plot.file,'...\n')
		setkey(df, ID1, ID2)
		tmp	<- unique(df)
		tmp	<- tmp[order(-PAIR_ID),]
		tmp[, LABEL:= factor(PAIR_ID, levels=PAIR_ID, labels=paste('Pair',PAIR_ID, ' (', ID1,'<->', ID2,')',sep=''))]
		tmp	<- merge(subset(tmp, select=c(PAIR_ID, LABEL)), df, by='PAIR_ID')
		set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, 	levels=c('anc_12','anc_21','sibling','int','unint','disconnected'), 
													labels=c('from 1 to 2','from 2 to 1','1, 2 are siblings','1, 2 are intermingled','1, 2 are unint','1, 2 are disconnected'))])
		ggplot(tmp, aes(x=LABEL, y=WIN_OF_TYPE, fill=TYPE)) +
				geom_bar(stat='identity', position='stack') +
				coord_flip() +
				labs(x='', y='number of read windows', fill='topology of clades from reads\nbetween patient pairs') +
				scale_fill_manual(values=c('from 1 to 2'="#9E0142",'from 2 to 1'="#F46D43",'1, 2 are siblings'="#ABDDA4",'1, 2 are intermingled'="#3288BD", '1, 2 are unint'='grey70', '1, 2 are disconnected'='grey50')) +
				theme_bw() + theme(legend.position='bottom') +
				guides(fill=guide_legend(ncol=2))
		ggsave(plot.file, w=10, h=0.15*nrow(unique(df)), limitsize = FALSE)
	}	
}

#' @import data.table ape
#' @title Read subtree information
#' @description This function reads the subtree information that is generated with the phyloscanner toolkit. 
#' @param prefix.infiles Full path name that identifies subtree files
#' @param prefix.run Character string to identify separate phyloscanner runs. After the prefix, an integer number is expected. For example, if prefix.run='ptyr', then files are expected to start like 'ptyr123_'.
#' @param regexpr.subtrees Regular expression that identifies subtree files in the directory. By default, this is 'Subtrees_r_.*\\.rda$' or 'Subtrees_c_.*\\.rda$' from the phyloscanner toolkit.
#' @param prefix.wfrom Character string to identify the start of a short read window. After the prefix, an integer number is expected. For example, if prefix.wfrom='Window_', then 'Window_800_to_1100' has start coordinate 800.
#' @param prefix.wto Character string to identify the end of a short read window. After the prefix, an integer number is expected. For example, if prefix.wto='Window_[0-9]+_to_', then 'Window_800_to_1100' has end coordinate 1100.
#' @param save.file If not missing, function output (a data.table) will be stored to 'save.file', input files will be zipped, and input files will be deleted.
#' @param resume If TRUE and save.file is not missing, the function loads and returns subtree info stored in save.file.
#' @param zip If TRUE and save.file is not missing, the function zips and removes subtree files that match the regular expression for subtrees.     
#' @return data.table with columns PTY_RUN W_FROM W_TO orig.patients patient.splits tip.names. 
phsc.read.subtrees<- function(prefix.infiles, prefix.run='ptyr', regexpr.subtrees='Subtrees_r_.*\\.rda$', prefix.wfrom='Window_', prefix.wto='Window_[0-9]+_to_', save.file=NA, resume=FALSE, zip=FALSE)
{
	if(!is.na(save.file) & resume)
	{
		options(show.error.messages = FALSE)		
		tmp		<- try(suppressWarnings(load(save.file)))
		options(show.error.messages = TRUE)
		if(!inherits(tmp, "try-error"))
		{
			cat('\nLoaded from file', save.file,'...\n')			
			return(rs.subtrees)
		}			
	}
	
	cat('\nReading tree summary files starting', prefix.infiles,'...\n')	
	dfr	<- data.table(FILE_TR= list.files(dirname(prefix.infiles), pattern=regexpr.subtrees, full.names=TRUE))
	dfr	<- subset(dfr, grepl(basename(prefix.infiles), basename(FILE_TR), fixed=1))
	dfr[, PTY_RUN:= as.integer(gsub(prefix.run,'',regmatches(FILE_TR,regexpr(paste(prefix.run,'[0-9]+',sep=''), FILE_TR))))]	
	dfr[, W_FROM:= as.integer(gsub(prefix.wfrom,'',regmatches(FILE_TR,regexpr(paste(prefix.wfrom,'[0-9]+',sep=''), FILE_TR))))]
	dfr[, W_TO:= as.integer(gsub(prefix.wto,'',regmatches(FILE_TR,regexpr(paste(prefix.wto,'[0-9]+',sep=''), FILE_TR))))]	
	setkey(dfr, PTY_RUN, W_FROM, W_TO)	
	cat('\nFound tree summary files, n=', nrow(dfr),'...\n')	
	rs.subtrees	<- lapply(seq_len(nrow(dfr)), function(i){				 
				#z	<- '/Users/Oliver/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/pty_Rakai_160825/ptyr1_InWindow_1050_to_1299_subtrees_r.rda'				
				load(dfr[i,FILE_TR])	
				rs.subtrees	<- as.data.table(rs.subtrees)
				rs.subtrees[, PTY_RUN:= dfr[i, PTY_RUN]]
				rs.subtrees[, W_FROM:= dfr[i, W_FROM]]
				rs.subtrees[, W_TO:= dfr[i, W_TO]]
				rs.subtrees
			})
	rs.subtrees	<- do.call('rbind', rs.subtrees)	
	if(!is.na(save.file))
	{
		cat('\nSave to file', save.file,'...\n')
		save(rs.subtrees, file=save.file)
	}
	if(zip & !is.na(save.file))
	{
		tmp	<- copy(dfr)		
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'_rda.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE_TR, flags = "-ur9XTjq")), by='FILE_TR'] )
			invisible( file.remove( tmp[, FILE_TR] ) )			
		}
		tmp	<- data.table(FILE= list.files(dirname(prefix.infiles), pattern=gsub('\\.rda','\\.csv',regexpr.subtrees), full.names=TRUE))
		tmp	<- subset(tmp, grepl(basename(prefix.infiles), basename(FILE), fixed=1))
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'_csv.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE, flags = "-ur9XTjq")), by='FILE'] )			
			invisible( file.remove( tmp[, FILE] ) )			
		}						
	}
	rs.subtrees
}

#' @import data.table ape
#' @title Read short read trees
#' @description This function reads short read trees that are generated with the phyloscanner toolkit. 
#' @param prefix.infiles Full path name identifying likely transmissions summary files
#' @param prefix.run Character string to identify separate phyloscanner runs. After the prefix, an integer number is expected. For example, if prefix.run='ptyr', then files are expected to start like 'ptyr123_'.
#' @param regexpr.trees Regular expression that identifies tree summary files in the directory. By default, this is 'Subtrees_r_.*\\.rda$' or 'Subtrees_c_.*\\.rda$' from the phyloscanner toolkit.
#' @param prefix.wfrom Character string to identify the start of a short read window. After the prefix, an integer number is expected. For example, if prefix.wfrom='Window_', then 'Window_800_to_1100' has start coordinate 800.
#' @param prefix.wto Character string to identify the end of a short read window. After the prefix, an integer number is expected. For example, if prefix.wto='Window_[0-9]+_to_', then 'Window_800_to_1100' has end coordinate 1100.
#' @param save.file If not missing, function output (a data.table) will be stored to 'save.file', input files will be zipped, and input files will be deleted.
#' @param resume If TRUE and save.file is not missing, the function loads and returns trees stored in save.file.
#' @param zip If TRUE and save.file is not missing, the function zips and removes trees and tree pdfs that match the regular expression for trees.     
#' @return list of named trees in ape format. 
phsc.read.trees<- function(prefix.infiles, prefix.run='ptyr', regexpr.trees='Subtrees_r_.*\\.rda$', prefix.wfrom='Window_', prefix.wto='Window_[0-9]+_to_', save.file=NA, resume=FALSE, zip=FALSE)
{
	if(!is.na(save.file) & resume)
	{
		options(show.error.messages = FALSE)		
		tmp		<- try(suppressWarnings(load(save.file)))
		options(show.error.messages = TRUE)
		if(!inherits(tmp, "try-error"))
		{
			cat('\nLoaded from file', save.file,'...\n')
			ans	<- list(phs=phs, dfr=dfr)
			return(ans)
		}			
	}
		
	cat('\nReading tree summary files starting', prefix.infiles,'...\n')	
	dfr	<- data.table(FILE_TR= list.files(dirname(prefix.infiles), pattern=regexpr.trees, full.names=TRUE))
	dfr	<- subset(dfr, grepl(basename(prefix.infiles), basename(FILE_TR),fixed=1))
	dfr[, PTY_RUN:= as.integer(gsub(prefix.run,'',regmatches(FILE_TR,regexpr(paste(prefix.run,'[0-9]+',sep=''), FILE_TR))))]	
	dfr[, W_FROM:= as.integer(gsub(prefix.wfrom,'',regmatches(FILE_TR,regexpr(paste(prefix.wfrom,'[0-9]+',sep=''), FILE_TR))))]
	dfr[, W_TO:= as.integer(gsub(prefix.wto,'',regmatches(FILE_TR,regexpr(paste(prefix.wto,'[0-9]+',sep=''), FILE_TR))))]	
	setkey(dfr, PTY_RUN, W_FROM, W_TO)
	dfr[, IDX:= seq_len(nrow(dfr))]
	cat('\nFound tree summary files, n=', nrow(dfr),'...\n')	
	phs	<- lapply(seq_len(nrow(dfr)), function(i){				 
				#z	<- '/Users/Oliver/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/pty_Rakai_160825/ptyr1_InWindow_1050_to_1299_subtrees_r.rda'				
				load(dfr[i,FILE_TR])
				tree				
			})
	names(phs)	<- dfr[, paste(PTY_RUN,'_',W_FROM,'_',W_TO,sep='')]
	if(!is.na(save.file))
	{
		cat('\nSave to file', save.file,'...\n')
		save(phs, dfr, file=save.file)
	}
	if(zip & !is.na(save.file))
	{
		tmp	<- data.table(FILE= list.files(dirname(prefix.infiles), pattern=gsub('\\.rda','\\.pdf',gsub("subtrees","tree",gsub("Subtrees","Tree",regexpr.trees))), full.names=TRUE))
		tmp	<- subset(tmp, grepl(basename(prefix.infiles), basename(FILE),fixed=1))
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'_pdf.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE, flags = "-ur9XTjq")), by='FILE'] )
			invisible( file.remove( tmp[, FILE] ) )			
		}
		tmp	<- data.table(FILE= list.files(dirname(prefix.infiles), pattern=paste('^',basename(prefix.infiles),'.*tree$',sep=''), full.names=TRUE))
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'_newick.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE, flags = "-ur9XTjq")), by='FILE'] )			
			invisible( file.remove( tmp[, FILE] ) )			
		}
		tmp	<- data.table(FILE= list.files(dirname(prefix.infiles), pattern=paste('^ProcessedTree_.*',basename(prefix.infiles),'.*tree$',sep=''), full.names=TRUE))
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'_processed_nexus.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE, flags = "-ur9XTjq")), by='FILE'] )			
			invisible( file.remove( tmp[, FILE] ) )			
		}		
		tmp	<- data.table(FILE= list.files(dirname(prefix.infiles), pattern=paste('^',basename(prefix.infiles),'.*fasta$',sep=''), full.names=TRUE))
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'_fasta.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE, flags = "-ur9XTjq")), by='FILE'] )
			invisible( file.remove( tmp[, FILE] ) )	
		}	
		tmp	<- data.table(FILE= list.files(dirname(prefix.infiles), pattern=paste('^',basename(prefix.infiles),'.*DuplicateReadCounts.*csv$',sep=''), full.names=TRUE))
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'_DuplicateReadCounts.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE, flags = "-ur9XTjq")), by='FILE'] )			
			invisible( file.remove( tmp[, FILE] ) )			
		}		
		tmp	<- data.table(FILE= list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*blacklist.*csv$',sep=''), full.names=TRUE))
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'_blacklist.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE, flags = "-ur9XTjq")), by='FILE'] )			
			invisible( file.remove( tmp[, FILE] ) )			
		}
		tmp	<- data.table(FILE= list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*duallist.*csv$',sep=''), full.names=TRUE))
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'_duallist.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE, flags = "-ur9XTjq")), by='FILE'] )			
			invisible( file.remove( tmp[, FILE] ) )			
		}
		tmp	<- data.table(FILE= list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*collapsed\\.csv$',sep=''), full.names=TRUE))
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'_collapsed.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE, flags = "-ur9XTjq")), by='FILE'] )			
			invisible( file.remove( tmp[, FILE] ) )			
		}		
	}
	list(phs=phs, dfr=dfr)
}

#' @import data.table 
#' @title Get transmission assignments per window
#' @description This function extracts all transmission assignments per window for two individuals. 
#' @param id1 Identifier of the first individual. This is the entry for that individual in the phyloscanner .bam file.
#' @param id2 Identifier of the second individual. 
#' @param infiles One or more file names of transmission stats rda files. These end typically with '_trmStatsPerWindow.rda'.  
#' @return data.table with columns ID1 ID2 W_FROM (window start) W_TO (window end), and further columns for each infile that give the transmission assignments per window. 
phsc.get.assignments.by.window.for.couple<- function(id1, id2, infiles)
{
	tmp	<- data.table(F= infiles)
	tmp[, RUN:=tmp[, basename(dirname(F))]]
	df		<- tmp[, {
				#load( '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d20_rerun/ptyr66_trmStatsPerWindow.rda' )
				load(F)
				df		<- subset(window.table, (pat.1==id1 & pat.2==id2) | (pat.1==id2 & pat.2==id1))
				setnames(df, c('pat.1','pat.2','type'), c('ID1','ID2','TYPE'))
				set(df, NULL, 'SOURCE_FILE', NULL)
				set(df, df[, which(is.na(TYPE))], 'TYPE', 'disconnected')
				set(df, df[, which(ID1==id1 & TYPE=='anc')], 'TYPE', 'anc_12')
				set(df, df[, which(ID1==id2 & TYPE=='anc')], 'TYPE', 'anc_21')
				set(df, NULL, 'ID1', id1)
				set(df, NULL, 'ID2', id2)
				df
			}, by='RUN']
	setkey(df, ID1, ID2, RUN, W_FROM)
	df		<- dcast.data.table(df, ID1+ID2+W_FROM+W_TO~RUN, value.var='TYPE')
	df
}