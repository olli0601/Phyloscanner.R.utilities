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
		cmd	<- paste0(cmd, ' --tipRegex "', tipRegex,'"')
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

phsc.cmd.blacklist.parsimonybased<- function(pr, scriptdir, inputFileName, outputFileName, dualCandidatesOutputFileName=NA, blackListFileName=NA, rawThreshold=1, ratioThreshold=1/200, sankhoffK=20, outgroupName=NA, tipRegex=NA)
{
	cmd	<- paste(pr, '--scriptdir',scriptdir, rawThreshold, ratioThreshold, sankhoffK, inputFileName, outputFileName)
	if(!is.na(dualCandidatesOutputFileName))
		cmd	<- paste(cmd, '--dualsOutputFile', dualCandidatesOutputFileName)			
	if(!is.na(outgroupName))
		cmd	<- paste(cmd, '--outgroupName', outgroupName)
	if(!is.na(blackListFileName))
		cmd	<- paste(cmd, '--blacklist', blackListFileName)	
	if(!is.na(tipRegex))
		cmd	<- paste0(cmd, ' --tipRegex "', tipRegex, '"')
	cmd
}

phsc.cmd.blacklist.rogue.geneticdistance<- function(pr, scriptdir, inputTreeFileName, outputFileName, longestBranchLength=0.04, dropProportion=1/100, blackListFileName=NA, outgroupName=NA, tipRegex=NA)
{
	cmd	<- paste(pr, '--scriptdir',scriptdir, dropProportion, longestBranchLength, inputTreeFileName, outputFileName)
	if(!is.na(tipRegex))
		cmd	<- paste0(cmd, ' --tipRegex "', tipRegex,'"')
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
		cmd	<- paste0(cmd, ' --tipRegex "', tipRegex,'"')
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


phsc.cmd.SummaryStatistics<- function(pr, scriptdir, outgroupName, file.patients, treeFiles, splitsFiles, outputBaseName, tipRegex=NA, blacklistFiles=NA)
{
	cmd	<- paste(pr,' --scriptdir ',scriptdir,' --outgroupName ', outgroupName, ' "', file.patients, '" "', treeFiles, '" "',splitsFiles, '" "',outputBaseName, '"', sep='')	
	if(!is.na(tipRegex))
		cmd	<- paste(cmd, ' --tipRegex "',tipRegex,'"', sep='')	
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

phsc.cmd.SplitPatientsToSubtrees<- function(pr, scriptdir, infile, outputdir=NA, blacklistFiles=NA, outputFileIdentifier=NA, outgroupName=NA, splitsRule=NA, sankhoff.k=NA, tipRegex=NA, pdfwidth=30, pdfrelheight=0.15)	
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
	if(!is.na(tipRegex))
		cmd	<- paste(cmd, ' --tipRegex "', tipRegex, '"', sep='')		
	if(!is.na(pdfwidth))
		cmd	<- paste(cmd, ' --pdfwidth ', pdfwidth, sep='')
	if(!is.na(pdfrelheight))
		cmd	<- paste(cmd, ' --pdfrelheight ', pdfrelheight, sep='')
	cmd	
}
	
phsc.cmd.LikelyTransmissions<- function(pr, scriptdir, file.tree, file.splits, file.out, tipRegex=NA)
{
	cmd<- paste(pr, ' "', file.tree, '" "', file.splits, '" "', file.out,'" ',' --scriptdir ', scriptdir, sep='')
	if(!is.na(tipRegex))
		cmd	<- paste(cmd, ' --tipRegex "', tipRegex,'"', sep='')			
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
	if(!any(colnames(dwin)=='UNINTERRUPTED'))	#for backward compatibility
		dwin[, UNINTERRUPTED:=FALSE]
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
	tmp		<- dwin[, which(TYPE=="anc_12" & (CONTIGUOUS|UNINTERRUPTED))]
	cat('\nFound contiguous anc_12, n=', length(tmp),'--> chain with no intermediate')
	set(dwin, tmp, 'TYPE', 'chain_12_nointermediate')
	tmp		<- dwin[, which(TYPE=="multi_anc_12" & (CONTIGUOUS|UNINTERRUPTED))]
	cat('\nFound contiguous multi_anc_12, n=', length(tmp),'--> chain with no intermediate')
	set(dwin, tmp, 'TYPE', 'chain_12_nointermediate')
	#	chains with no intermediate
	tmp		<- dwin[, which(TYPE=="anc_21" & (CONTIGUOUS|UNINTERRUPTED))]
	cat('\nFound contiguous anc_21, n=', length(tmp),'--> chain with no intermediate')
	set(dwin, tmp, 'TYPE', 'chain_21_nointermediate')
	tmp		<- dwin[, which(TYPE=="multi_anc_21" & (CONTIGUOUS|UNINTERRUPTED))]
	cat('\nFound contiguous multi_anc_21, n=', length(tmp),'--> chain with no intermediate')
	set(dwin, tmp, 'TYPE', 'chain_21_nointermediate')	
	#
	#	chains with intermediate
	tmp		<- dwin[, which(TYPE=="anc_12" & !(CONTIGUOUS|UNINTERRUPTED))]
	cat('\nFound discontiguous anc_12, n=', length(tmp),'--> chain with intermediate')
	set(dwin, tmp, 'TYPE', 'chain_12_withintermediate')
	tmp		<- dwin[, which(TYPE=="multi_anc_12" & !(CONTIGUOUS|UNINTERRUPTED))]
	cat('\nFound discontiguous multi_anc_12, n=', length(tmp),'--> chain with intermediate')
	set(dwin, tmp, 'TYPE', 'chain_12_withintermediate')
	#	chains with intermediate
	tmp		<- dwin[, which(TYPE=="anc_21" & !(CONTIGUOUS|UNINTERRUPTED))]
	cat('\nFound discontiguous anc_21, n=', length(tmp),'--> chain with intermediate')
	set(dwin, tmp, 'TYPE', 'chain_21_withintermediate')
	tmp		<- dwin[, which(TYPE=="multi_anc_21" & !(CONTIGUOUS|UNINTERRUPTED))]
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
	tmp		<- dwin[, which((CONTIGUOUS|UNINTERRUPTED) & PATHS_12==0 & PATHS_21==0)]
	cat('\nFound contiguous with no paths, n=', length(tmp),'--> other')
	set(dwin, tmp, 'TYPE', 'other')
	tmp		<- dwin[, which(!(CONTIGUOUS|UNINTERRUPTED) & PATHS_12==0 & PATHS_21==0)]
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
	cmd			<- paste("CWD=$(pwd)\n",sep='\n')
	cmd			<- paste(cmd,"echo $CWD\n",sep='')
	tmpdir		<- paste('pty','_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')	
	tmpdir		<- paste("$CWD/",tmpdir,sep='')
	cmd			<- paste(cmd,'mkdir -p "',tmpdir,'"\n',sep='')
	#	copy required files to local tmp dir	
	file.patient<- list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*patients.txt',sep=''), full.names=TRUE)
	stopifnot(length(file.patient)==1)	
	cmd		<- paste(cmd,'cp "',file.patient,'" "',tmpdir,'"\n',sep='')	
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
	cmd		<- paste(cmd, phsc.cmd.process.phyloscanner.output.in.directory(tmpdir, file.patient, pty.args), collapse='\n',sep='')
	#	move all files starting with current run ID
	run.id	<- gsub('_patients.txt','',basename(file.patient))
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
#' @param pty.args List of phyloscanner input variables. See examples.
#' @param file.bam File name of the file that contains the list of bam files.
#' @param file.ref File name of the file that contains the list of reference files.
#' @param file.patient File name of the file that contains the list of unique individuals/units to which the bam files correspond. Multiple bam files are allowed per individual/unit.
#' @param file.rename File name of the file that contains new ID names for each bam/ref file. Can be set to NA. 
#' @return Character string of phyloscanner commands.
#' @description This function generates bash commands for a single phyloscanner run, that can be called via 'system' in R, or written to file to run on a UNIX system.
#' @example example/ex.cmd.phyloscanner.one.R    
phsc.cmd.phyloscanner.one<- function(pty.args, file.bam, file.ref, file.patient, file.rename=NA)
{	
	stopifnot(is.character(file.bam),is.character(file.ref),is.character(file.patient))
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
	cmd		<- paste(cmd,'cp "',file.patient,'" "',tmpdir,'"\n',sep='')
	if(!is.na(file.rename))
		cmd	<- paste(cmd,'cp "',file.rename,'" "',tmpdir,'"\n',sep='')
	#	cd to tmp dir
	cmd		<- paste(cmd, 'cd "',tmpdir,'"\n', sep='')	
	cmd		<- paste(cmd, prog.pty,' "',basename(file.bam),'" "',basename(file.ref),'" ',sep='')	
	cmd		<- paste(cmd, '--merging-threshold', merge.threshold,'--min-read-count',min.read.count,'--quality-trim-ends', quality.trim.ends, '--min-internal-quality',min.internal.quality)
	if(!is.na(merge.paired.reads))
		cmd	<- paste(cmd,' ',merge.paired.reads,sep='')	
	if(!is.na(dont.check.duplicates))
		cmd	<- paste(cmd,' ',dont.check.duplicates,sep='')
	if(!is.na(file.rename))	
		cmd	<- paste(cmd,' --renaming-file "',basename(file.rename),'"',sep='')
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
	cmd		<- paste(cmd, phsc.cmd.process.phyloscanner.output.in.directory(tmpdir, file.patient, pty.args), sep='\n')
	#
	cmd		<- paste(cmd, '\nmv ',run.id,'* "',out.dir,'"\n',sep='')	
	#	zip up everything else
	if(is.null(mem.save) || is.na(mem.save) || mem.save==0)
	{
		cmd		<- paste(cmd, 'for file in *; do\n\tzip -ur9XTjq ',paste(run.id,'_otherstuff.zip',sep=''),' "$file"\ndone\n',sep='')
		cmd		<- paste(cmd, 'mv ',paste(run.id,'_otherstuff.zip',sep=''),' "',out.dir,'"\n',sep='')		
	}
	#	clean up
	cmd		<- paste(cmd,'cd $CWD\nrm -r "',tmpdir,'"\n',sep='')
	detach(pty.args)
	cmd
}

#' @export
#' @title Generate bash commands for a multiple phyloscanner runs
#' @param pty.runs Data.table of individual assignments to phyloscanner runs, with columns 'PTY_RUN' (run id), 'SAMPLE_ID' (ID of individuals that are assigned to that run). Optional columns: 'RENAME_ID' (new ID for each bam file in phyloscanner output).
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
	ptyd[, SAMPLE_ID:= gsub('\\.bam|_ref\\.fasta','',basename(FILE))]
	ptyd		<- dcast.data.table(ptyd, SAMPLE_ID~TYPE, value.var='FILE')
	#	merge
	ptyd	<- merge(pty.runs, ptyd, by='SAMPLE_ID', all.x=1)
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
				#	PTY_RUN<- z <- 1; BAM<- subset(ptyd, PTY_RUN==z)[, BAM]; REF<- subset(ptyd, PTY_RUN==z)[, REF]									
				file.bam		<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_bam.txt',sep=''))
				file.ref		<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_ref.txt',sep=''))				
				cat( paste(BAM[!is.na(BAM)&!is.na(REF)],collapse='\n'), file= file.bam	)
				cat( paste(REF[!is.na(BAM)&!is.na(REF)],collapse='\n'), file= file.ref	)				
				file.rename		<- NA_character_
				if(!is.null(RENAME_ID))
				{
					file.rename	<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_rename.txt',sep=''))
					cat( paste(RENAME_ID[!is.na(BAM)&!is.na(REF)],collapse='\n'), file= file.rename	)					
				}	
				file.patient	<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_patient.txt',sep=''))
				tmp				<- SAMPLE_ID[!is.na(BAM)&!is.na(REF)]
				if(!is.null(UNIT_ID))
					tmp			<- unique(UNIT_ID[!is.na(BAM)&!is.na(REF)])
				cat( paste(tmp,collapse='\n'), file= file.patient	)
				cmd			<- phsc.cmd.phyloscanner.one(	pty.args, 
															file.bam, 
															file.ref,
															file.patient,
															file.rename)
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
#' @param file.patients File name of the file that contains the list of unique individuals/units that the bam files correspond to, possibly after re-naming.
#' @param pty.args List of phyloscanner input variables. 
#' @return character string of bash commands.
phsc.cmd.process.phyloscanner.output.in.directory<- function(tmp.dir, file.patients, pty.args)
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
	tip.regex						<- pty.args[['tip.regex']]
	mem.save						<- pty.args[['mem.save']]
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
	run.id							<- gsub('_rename.txt|_bam.txt|_patient.txt|_patients.txt','',basename(file.patients))
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
	if(is.null(tip.regex))
		tip.regex					<- NA_character_
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
		tmp				<- phsc.cmd.blacklist.reads(	prog.pty.readblacklist, 
														'"$file"', 
														paste('"${file//DuplicateReadCounts/blacklist}"',sep=''), 
														contaminant.read.threshold, 
														contaminant.prop.threshold, 
														tipRegex=tip.regex)
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
		tmp				<- ifelse(any(is.na(blacklistFiles)), NA_character_, '"${TMP//InWindow/blacklist_InWindow}"')
		tmp				<- phsc.cmd.blacklist.parsimonybased( 	prog.pty.readblacklistsankoff, 
																	pty.tools.dir,
																	'"$file"', 
																	paste('"${TMP//InWindow/blacklistsank_InWindow}"',sep=''), 
																	dualCandidatesOutputFileName=paste('"${TMP//InWindow/duallistsank_InWindow}"',sep=''),
																	blackListFileName=tmp,
																	rawThreshold=roguesubtree.read.threshold, 
																	ratioThreshold=roguesubtree.prop.threshold, 
																	sankhoffK=sankhoff.k,	
																	outgroupName=root.name,
																	tipRegex=tip.regex)
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
		tmp				<- ifelse(any(is.na(blacklistFiles)), NA_character_, '"${TMP//InWindow/blacklist_InWindow}"')
		tmp				<- phsc.cmd.blacklist.rogue.geneticdistance(	prog.pty.rogueblacklist, 
															pty.tools.dir, 
															'"$file"', 
															'"${TMP//InWindow/blacklistrogue_InWindow}"', 
															longestBranchLength=rogue.longestBranchLength, 
															dropProportion=rogue.prop.threshold, 
															blackListFileName=tmp, 
															outgroupName=root.name, 
															tipRegex=tip.regex)				
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
															tipRegex=tip.regex)				
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
																tipRegex=tip.regex)
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
															tipRegex=tip.regex,
															pdfwidth=30, 
															pdfrelheight=0.15)
	cmd				<- paste(cmd, tmp, sep='\n')
	#
	#	bash command to calculate patient stats
	#							
	tmp				<- phsc.cmd.SummaryStatistics( 	prog.pty.smry, 
													pty.tools.dir, 
													root.name, 
													file.path(tmp.dir, basename(file.patients)), 
													file.path(tmp.dir, paste(run.id_,'InWindow_',sep='')), 
													file.path(tmp.dir, paste('Subtrees_',ifelse(use.sankhoff.blacklister, 's', 'r'),'_',run.id_,'InWindow_',sep='')), 
													file.path(tmp.dir, substr(run.id_,1,nchar(run.id_)-1)),
													tipRegex=tip.regex,
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
															file.path(tmp.dir, basename(file.patients)),												
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
	#
	#	if mem-save output delete: .*subtrees_[scr]_csv.zip  .*subtrees_[scr]_rda.zip .*DuplicateReadCounts.zip .*_blacklist.zip .*_duallist.zip .*_collapsed.zip .*_LikelyTransmissions.zip
	#	
	if(!is.null(mem.save) & mem.save==1)
	{
		tmp				<- c('*_subtrees_r.rda','*subtrees_[scr]_rda.zip','*DuplicateReadCounts.zip','*_blacklist.zip','*_duallist.zip','*_collapsed.zip','*_LikelyTransmissions.zip')
		tmp				<- paste('rm ',tmp, sep='',collapse='\n')
		cmd				<- paste(cmd, tmp, sep='\n')	
	}	
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
		phsc.read.likelytransmissions(prefix.infiles, prefix.run='ptyr', regexpr.lklsu='trmStats.csv$', save.file=save.file, resume=resume, zip=zip)				
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
phsc.plot.phy.selected.pairs<- function(phs, dfs, id1, id2, plot.file=NA, pdf.h=50, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40)
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
#' @import data.table
#' @title Default colours for relationship types
#' @description This function returns default colours for each relationship type in one of the defined relationship groups.  
#' @return Named list with names corresponding to RELATIONSHIP GROUPS. For each relationship group, the list contains a named vector. Entries in this vector are colours, and names are relationship types.
phsc.plot.default.colours.for.relationtypes<- function()
{
	cols.type	<- list()
	tmp2		<- do.call('rbind',list(
					data.table(	TYPE= c("chain fm\nno intermediate\nclose","chain fm\nno intermediate","chain fm\nno intermediate\ndistant"),
							COLS= brewer.pal(11, 'PiYG')[c(1,2,4)]),
					data.table(	TYPE= c("chain mf\nno intermediate\nclose","chain mf\nno intermediate","chain mf\nno intermediate\ndistant"),
							COLS= brewer.pal(11, 'PuOr')[c(1,2,4)]),
					data.table(	TYPE= c("intermingled\nno intermediate\nclose","intermingled\nno intermediate","intermingled\nno intermediate\ndistant"),
							COLS= brewer.pal(11, 'PRGn')[c(1,2,4)]),
					data.table(	TYPE= c("chain fm\nwith intermediate\nclose","chain fm\nwith intermediate","chain fm\nwith intermediate\ndistant"),
							COLS= rev(brewer.pal(11, 'BrBG'))[c(3,4,5)]),
					data.table(	TYPE= c("chain mf\nwith intermediate\nclose","chain mf\nwith intermediate","chain mf\nwith intermediate\ndistant"),
							COLS= rev(brewer.pal(11, 'PRGn'))[c(3,4,5)]),
					data.table(	TYPE= c("intermingled\nwith intermediate\nclose","intermingled\nwith intermediate","intermingled\nwith intermediate\ndistant"),
							COLS= rev(brewer.pal(11, 'RdBu'))[c(3,4,5)]),
					data.table(	TYPE= c("other close","other","other distant"),
							COLS= rev(brewer.pal(11, 'RdGy'))[c(3,4,5)])))
	cols.type[['TYPE_DIR_TODI7x3']]	<- { tmp<- tmp2[, COLS]; names(tmp) <- tmp2[, TYPE]; tmp }
	tmp2		<- do.call('rbind',list(
					data.table(	TYPE= c("pair close","pair","pair distant"),
							COLS= brewer.pal(11, 'PuOr')[c(1,2,4)]),
					data.table(	TYPE= c("withintermediate close","withintermediate","withintermediate distant"),
							COLS= rev(brewer.pal(11, 'RdBu'))[c(3,4,5)]),
					data.table(	TYPE= c("other close","other","other distant"),
							COLS= rev(brewer.pal(11, 'RdGy'))[c(3,4,5)])))
	tmp2		<- { tmp<- tmp2[, COLS]; names(tmp) <- tmp2[, TYPE]; tmp }
	cols.type[['TYPE_PAIR_TODI3x3']]	<- tmp2	
	tmp2		<- do.call('rbind',list(
					data.table(	TYPE= c('close ancestral/\nintermingled', 'not close ancestral/\nintermingled'),
							COLS= brewer.pal(11, 'PRGn')[c(2,4)]),
					data.table(	TYPE= c('close other','not close other'),
							COLS= rev(brewer.pal(11, 'RdGy'))[c(3,5)])))
	tmp2		<- { tmp<- tmp2[, COLS]; names(tmp) <- tmp2[, TYPE]; tmp }
	cols.type[['TYPE_PAIR_TODI2x2']]	<- tmp2	
	tmp2		<- do.call('rbind',list(
					data.table(	TYPE= c('likely pair', 'chain'),
							COLS= brewer.pal(11, 'PRGn')[c(2,4)]),
					data.table(	TYPE= c('intermdiate distance','distant'),
							COLS= rev(brewer.pal(11, 'RdGy'))[c(3,5)])))
	tmp2		<- { tmp<- tmp2[, COLS]; names(tmp) <- tmp2[, TYPE]; tmp }
	cols.type[['TYPE_PAIRCHAIN_TODI']]	<- tmp2	
	tmp2		<- data.table(	TYPE= c("no intermediate\n and close", "no intermediate\n but not close", "with intermediate\nor distant"),
			COLS= c(brewer.pal(11, 'RdBu')[c(2,4)], rev(brewer.pal(11, 'RdGy'))[4]))					
	tmp2		<- { tmp<- tmp2[, COLS]; names(tmp) <- tmp2[, TYPE]; tmp }
	cols.type[['TYPE_PAIR_TODI3']]	<- tmp2
	tmp2		<- data.table(	TYPE= c("close", "intermediate\ndistance", "distant"),
			COLS= c(rev(brewer.pal(11, 'RdBu'))[c(2,4)], rev(brewer.pal(11, 'RdGy'))[4]))					
	tmp2		<- { tmp<- tmp2[, COLS]; names(tmp) <- tmp2[, TYPE]; tmp }
	cols.type[['TYPE_PAIR_DI']]	<- tmp2
	tmp2		<- data.table(	TYPE= c("ancestral/\nintermingled", "other"),
			COLS= c(rev(brewer.pal(9, 'Greens'))[2], rev(brewer.pal(11, 'RdGy'))[4]))					
	tmp2		<- { tmp<- tmp2[, COLS]; names(tmp) <- tmp2[, TYPE]; tmp }
	cols.type[['TYPE_PAIR_TO']]	<- tmp2
	
	cols.type
}

#' @export
#' @import data.table ggplot2 scales
#' @title Compare window assignments by different phyloscanner runs
#' @description This function plots the window assignments for each pair of individuals and for several phyloscanner runs.  
#' @param rplkl2 Posterior probability assignments for pairs of individuals. This is specified as a data.table with columns RUN, PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, TYPE, KEFF, POSTERIOR_ALPHA, POSTERIOR_BETA 
#' @param plot.file Name of output file
#' @param cols Colours for phyloscanner runs. This is specified as a named vector which must not be missing.
#' @param cols Colours for relationship types. This is specified as a named vector, which can be null. 
#' @param group Relationship group name. This is used to define default colours if cols==NULL.   
#' @param height Total height of pdf in inches. This is overwritten if group is specified.
#' @return Plots to file.
phsc.plot.windowassignments.by.runs <- function(rplkl2, plot.file, plot.prob.select, cols.run, cols=NULL, group=NA, height=40)
{
	#	height<- 40
	g_legend<-function(a.gplot)
	{
		tmp <- ggplot_gtable(ggplot_build(a.gplot))
		leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
		legend <- tmp$grobs[[leg]]
		legend
	}
	#	define colours to be used
	if(is.null(cols))
	{
		stopifnot(!is.na(group))		
		cols.type	<- phsc.plot.default.colours.for.relationtypes()
		stopifnot(group%in%names(cols.type))
		cols		<- cols.type[[group]]
	}
	#	re-order pairs
	setkey(rplkl2, LABEL_SH)	
	tmp		<- subset(rplkl2, TYPE%in%plot.prob.select)
	tmp2	<- tmp[, list(D=mean(POSTERIOR_ALPHA/(POSTERIOR_ALPHA+POSTERIOR_BETA))), by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','LABEL_SH')]
	tmp2	<- tmp2[order(-D),]
	tmp2[, PLOT_ID:= 1:nrow(tmp2)]
	set(tmp2, NULL, 'PLOT_ID', tmp2[, factor(PLOT_ID, levels=PLOT_ID, labels=LABEL_SH)])
	setnames(tmp2, c('PLOT_ID','LABEL_SH'),c('LABEL_SH','OLD_LABEL_SH'))
	set(tmp2, NULL, c('OLD_LABEL_SH','D'), NULL)
	set(rplkl2, NULL, c('LABEL_SH'), NULL)
	rplkl2	<- merge(rplkl2, tmp2, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID'))
	setkey(rplkl2, LABEL)
	p1		<- ggplot(rplkl2, aes(x=RUN, y=KEFF, fill=TYPE, colour=RUN)) + 
			geom_bar(stat='identity',position='stack') +
			scale_y_continuous(expand=c(0,0)) +
			scale_fill_manual(values=cols) +
			scale_colour_manual(values=cols.run) +
			theme_bw() + 
			theme(	legend.position='bottom', axis.text.y=element_blank(),
					axis.ticks.y=element_blank(),
					panel.spacing=unit(0.4, "lines"), strip.text.y = element_text(angle=180),
					strip.background=element_rect(fill="transparent", colour="transparent"),
					panel.border=element_rect(color="transparent")) +
			facet_grid(LABEL_SH~., switch='y') +
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE), colour=guide_legend(ncol=1, byrow = TRUE)) +
			labs(	x='', 
					y='non-overlapping windows\n(number)\n',
					fill='phylogenetic\nrelationship')
	tmp		<- subset(rplkl2, TYPE%in%plot.prob.select)
	p2		<- ggplot(tmp, aes(x=RUN, fill=TYPE, colour=RUN, 
							middle= qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA), 
							lower=qbeta(0.25, POSTERIOR_ALPHA, POSTERIOR_BETA), 
							upper=qbeta(0.75, POSTERIOR_ALPHA, POSTERIOR_BETA), 
							ymin=qbeta(0.025, POSTERIOR_ALPHA, POSTERIOR_BETA), 
							ymax=qbeta(0.975, POSTERIOR_ALPHA, POSTERIOR_BETA))) + 
			geom_boxplot(stat='identity') +
			scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.2), labels=percent) +
			scale_fill_manual(values=cols) +
			scale_colour_manual(values=cols.run) +
			theme_bw() + 
			theme(	legend.position='bottom', axis.text.y=element_blank(),
					axis.ticks.y=element_blank(),
					panel.spacing=unit(0.4, "lines"), strip.text.y=element_blank(),
					strip.background=element_blank(),
					panel.border=element_rect(color="transparent")) +
			facet_grid(LABEL_SH~., switch='y') +
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE), colour=guide_legend(ncol=1, byrow = TRUE)) +
			labs(	x='', 
					y='posterior probability\n\n',
					fill='phylogenetic\nrelationship')
	p3		<- g_legend(p1)	
	p3$vp	<- viewport(layout.pos.row=2, layout.pos.col=1:2)	
	pdf(file=plot.file, w=12, h=height)
	grid.newpage()	
	pushViewport(viewport(layout = grid.layout(2, 2, heights=unit(c(40,1), "null"), widths=unit(c(7, 3), "null"))))   	
	print(p1+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
	print(p2+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
	grid.draw(p3)
	#pushViewport(viewport(layout = grid.layout(1, 2, heights=unit(c(10,1), "null"), widths=unit(c(7, 3), "null"))))   	
	#print(p1+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
	#print(p2+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))	
	dev.off()	
}

#' @export
#' @import data.table ggplot2 scales
#' @title Plot window summaries for pairs of individuals
#' @description This function plots the window summaries for each pair of individuals.  
#' @param plot.select Select pairs for plotting. This is specified as a data.table with columns PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL. Each pair has assigned a label that is used as title.
#' @param rpw2 Window assignments for pairs of individuals. This is specified as a data.table with columns PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, W_FROM, ID_R_MAX, ID_R_MIN, TYPE.  
#' @param rplkl2 Posterior probability assignments for pairs of individuals. This is specified as a data.table with columns PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, TYPE, KEFF, POSTERIOR_ALPHA, POSTERIOR_BETA 
#' @param plot.file Name of output file
#' @param cols Colours for relationship types. This is specified as a named vector, which can be null. 
#' @param group Relationship group name. This is used to define default colours if cols==NULL.   
#' @param widhts Widths of subfigures in grid.layout. This is overwritten if group is specified. 
#' @param heights Heights of subfigures in grid.layout. This is overwritten if group is specified. 
#' @param height Total height of pdf in inches. This is overwritten if group is specified.
#' @return Plots to file.
phsc.plot.windowsummaries.for.pairs<- function(plot.select, rpw2, rplkl2, plot.file, cols=NULL, group=NA, widths=unit(c(4, 6), "null"), heights=unit(c(2, 3.5, 4, 5), "null"), height=9)
{
	#	widths	<- unit(c(4, 6), "null"); heights	<- unit(c(2, 3.5, 4, 5), "null"); height	<- 9
	
	#	define colours to be used
	if(is.null(cols))
	{
		stopifnot(!is.na(group))		
		cols.type	<- phsc.plot.default.colours.for.relationtypes()
		stopifnot(group%in%names(cols.type))
		cols		<- cols.type[[group]]
	}
	#	re-define dimensions if group specified
	if(!is.na(group))
	{
		if(group%in%c('TYPE_DIR_TODI7x3'))
		{
			widths	<- unit(c(4, 6), "null")
			heights	<- unit(c(2, 3.5, 4, 15), "null")
			height	<- 17
		}		
		if(group%in%c('TYPE_PAIR_TODI2x2','TYPE_PAIRCHAIN_TODI'))
		{
			heights	<- unit(c(2, 3.5, 4, 3.75), "null")
			height	<- 8
		}
		if(group%in%c('TYPE_PAIR_TODI3','TYPE_PAIR_DI'))
		{
			heights	<- unit(c(2, 3.5, 4, 3.5), "null")
			height	<- 7
		}	
	}
	pdf(file=plot.file, w=10, h=height)		
	setkey(plot.select, LABEL)		
	for(i in seq_len(nrow(plot.select)))
	{
		#i<- 85
		#pty_run	<- 38; id1		<- '16016_1_4'; id2		<- '15105_1_35'
		pty_run	<- plot.select[i, PTY_RUN]; id1		<- plot.select[i, FEMALE_SANGER_ID]; id2		<- plot.select[i, MALE_SANGER_ID]		
		tmp		<- subset(rpw2, PTY_RUN==pty_run & FEMALE_SANGER_ID==id1 & MALE_SANGER_ID==id2)
		p1		<- ggplot(tmp, aes(x=W_FROM)) +			
				geom_bar(aes(y=ID_R_MAX, colour=TYPE), stat='identity', fill='transparent') +
				geom_bar(aes(y=ID_R_MIN, fill=TYPE), stat='identity', colour='transparent') +
				labs(x='', y='number of reads', fill='phylogenetic\nrelationship\n', colour='phylogenetic\nrelationship\n') +
				scale_fill_manual(values=cols) +
				scale_colour_manual(values=cols) +
				scale_x_continuous(breaks=seq(0,1e4,500), minor_breaks=seq(0,1e4,100), limits=c(rpw2[, min(W_FROM)-100], rpw2[, max(W_FROM)+100])) +
				scale_y_log10(breaks=c(10,100,1000,1e4,1e5)) +
				theme_bw() + theme(legend.position='left') +			
				guides(fill=FALSE, colour=FALSE)
		p2		<- ggplot(tmp, aes(x=W_FROM, y=PATRISTIC_DISTANCE)) +
				geom_point(size=1) +					
				labs(x='window start\n\n', y='patristic distance') +
				scale_x_continuous(breaks=seq(0,1e4,500), minor_breaks=seq(0,1e4,100), limits=c(rpw2[, min(W_FROM)], rpw2[, max(W_FROM)])) +
				scale_y_log10(labels=percent, limits=c(0.001, 0.7), expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
				theme_bw() + theme(legend.position='left')
		tmp		<- subset(rplkl2, PTY_RUN==pty_run & FEMALE_SANGER_ID==id1 & MALE_SANGER_ID==id2)
		p3		<- ggplot(tmp, aes(x=TYPE, y=KEFF, fill=TYPE)) + geom_bar(stat='identity') +
				scale_fill_manual(values=cols) +
				theme_bw() + theme(legend.position='bottom') +
				coord_flip() + guides(fill=FALSE) +			
				labs(x='', y='\nnon-overlapping windows\n(number)', fill='phylogenetic\nrelationship\n')
		p4		<- ggplot(tmp, aes(x=TYPE, 	middle=qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA), 
								ymin=qbeta(0.025, POSTERIOR_ALPHA, POSTERIOR_BETA), 
								ymax=qbeta(0.975, POSTERIOR_ALPHA, POSTERIOR_BETA), 
								lower=qbeta(0.25, POSTERIOR_ALPHA, POSTERIOR_BETA), 
								upper=qbeta(0.75, POSTERIOR_ALPHA, POSTERIOR_BETA), 
								fill=TYPE)) + 
				geom_boxplot(stat='identity') +
				scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), limits=c(0,1), expand=c(0,0)) +
				scale_fill_manual(values=cols) +
				theme_bw() + theme(legend.position='right', legend.margin=margin(0, .1, 0, 1, "cm")) +
				coord_flip() + guides(fill=guide_legend(ncol=1)) +
				labs(x='', y='\nposterior probability\n', fill='phylogenetic\nrelationship\n')				
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(4, 2, heights=heights, widths=widths)))   
		grid.text(tmp[1,LABEL], gp=gpar(fontsize=10), vp=viewport(layout.pos.row = 1, layout.pos.col = 1:2))
		print(p1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
		print(p2, vp = viewport(layout.pos.row = 3, layout.pos.col = 1:2))         
		print(p3, vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
		print(p4, vp = viewport(layout.pos.row = 4, layout.pos.col = 2))
	}
	dev.off()
}

phsc.extract.clade <- function(phy, node, root.edge = 0, interactive = FALSE)
{
	n <- length(phy$tip.label)
	if (interactive) {
		cat("Click close to the node...\n")
		node <- identify(phy)$nodes
	} else {
		if (length(node) > 1) {
			node <- node[1]
			warning("only the first value of 'node' has been considered")
		}
		if (is.character(node)) {
			if (is.null(phy$node.label))
				stop("the tree has no node labels")
			node <- match(node, phy$node.label) + n
			if (is.na(node)) stop("'node' not among the node labels.")
		}
		if (node <= n)
			stop("node number must be greater than the number of tips")
	}
	if (node == n + 1L) return(phy)
	keep <- prop.part(phy)[[node - n]]
	drop.tip(phy, (1:n)[-keep], root.edge = root.edge, rooted = TRUE)
}

#' @import data.table phangorn ape
phsc.phy.collapse.monophyletic.clades<- function(ph, drop.blacklisted=TRUE, tip.regex='(.*)_read_([0-9]+)_count_([0-9]+)')
{
	#ph<- phs[[1]]
	#references.pattern<- 'REF'
	phb			<- data.table(	IDX=seq_along(ph$tip.label), 
								BAM=ph$tip.label, 
								IND= gsub(tip.regex,'\\1',ph$tip.label),
								COUNT=sub(tip.regex,'\\3',ph$tip.label),
								BLACKLISTED= is.na(attr(ph, 'INDIVIDUAL')[seq_len(Ntip(ph))])
								)
	set(phb, phb[, which(IND==COUNT)], 'COUNT', NA_character_)
	set(phb, NULL, 'COUNT', phb[, as.integer(COUNT)])
	set(phb, NULL, 'REF', phb[, is.na(COUNT)])	
	set(phb, phb[, which(REF)],'IND','REFERENCE')
	tmp			<- subset(phb, BLACKLISTED & !REF)[,IDX]
	if(drop.blacklisted & length(tmp))
	{
		ph		<- phsc.drop.tip(ph, tmp)
		phb		<- data.table(	IDX=seq_along(ph$tip.label), 
								BAM=ph$tip.label, 
								IND= gsub(tip.regex,'\\1',ph$tip.label),
								COUNT=sub(tip.regex,'\\3',ph$tip.label)		)
		set(phb, phb[, which(IND==COUNT)], 'COUNT', NA_character_)
		set(phb, NULL, 'COUNT', phb[, as.integer(COUNT)])
		set(phb, NULL, 'REF', phb[, is.na(COUNT)])	
		set(phb, phb[, which(REF)],'IND','REFERENCE')		
	}
	#	for each patient define mrca (MRCA)
	#	for each patient, 
	#		determine MRCA (mrca)
	#		calculate number of other individuals in clade below MRCA (diff)	
	z			<- phb[, {
				#IND<- "15861_1_37.bam"
				#IDX<- subset(phb, IND=="15861_1_37.bam")[, IDX]
				#print(IND)
				mrca	<- IDX		
				diff	<- 0L
				if(length(IDX)>1)
				{
					mrca	<- as.integer(getMRCA(ph, IDX))
					tmp		<- phsc.extract.clade(ph, mrca, root.edge=1)
					#print(tmp)
					diff	<- length(setdiff(gsub(tip.regex,'\\1',tmp$tip.label), IND))								
				}												
				list(MRCA=mrca, DIFF_IND=diff)							
			}, by=c('IND')]	
	z			<- subset(z, IND!='REFERENCE')
	#	for all patients that are not monophyletic, there could be several clades
	#		trace back all ancestors between tips of other individuals and MRCA
	#		for each such tip, construct the unique path to MRCA that is not on a previous path
	#		find change points on these unique paths below which the subtree contains reads from the current patient
	#		the idea is that all monophyletic clades of this patient must break off from one of these paths	
	tmp			<- subset(z, DIFF_IND>0, c(IND, MRCA))
	if(!nrow(tmp))
		z[, CLADE:=NA_integer_]
	if(nrow(tmp))
	{
		#	find all tips that belong to another patient than the current individual
		tmp			<- tmp[, {
					if(MRCA<=Ntip(ph))
						tmp	<- ph$tip.label[MRCA]
					if(MRCA>Ntip(ph))
						tmp	<- phsc.extract.clade(ph,MRCA)$tip.label
					list(MRCA=MRCA, MRCA_IND=IND, BAM=tmp)	
				}, by='IND']
		tmp[, IND:=NULL]
		zz			<- merge(tmp, subset(phb, select=c(BAM, IND, IDX)), by='BAM')
		zz			<- subset(zz, MRCA_IND!=IND)
		#	determine change points
		zz			<- zz[, {
					#IDX<- c(980,910,912,950); MRCA<- 1580; IND<- 'R1_RES669_S20_L001'; MRCA_IND<- 'R1_RES827_S23_L001'					
					#	determine paths to MRCA from each tip
					anc.group		<- Ancestors(ph, IDX)	
					if(!is.list(anc.group))
						anc.group	<- list(anc.group)
					anc.group	<- lapply(seq_along(anc.group), function(i)  anc.group[[i]][ seq_len(which(anc.group[[i]]==MRCA[1])-1)] )							
					#	determine unique paths until we hit a path that is already visited
					anc.join	<- lapply(seq_along(anc.group), function(i){	unique(unlist(lapply(seq_len(i), function(j) anc.group[[j]])))	})
					anc.join	<- c(NA,anc.join)
					anc.group	<- lapply(seq_along(anc.group), function(i)	setdiff(anc.group[[i]],anc.join[[i]])	)
					#	check which clades defined by the mrcas on the ancestor path contain at least one read from MRCA_IND
					tmp			<- lapply(seq_along(anc.group), function(i) sapply(anc.group[[i]], function(j)	any(grepl(MRCA_IND, phsc.extract.clade(ph,j)$tip.label))		)	)				
					#	determine lowest node (counting from tips) which contains at least one read from MRCA_IND
					tmp			<- lapply(seq_along(tmp), function(i){
								ans		<- NA_integer_	
								tmp2	<- integer(0)
								if(length(tmp[[i]]))
									tmp2<- which(tmp[[i]]) 
								if(length(tmp2))
									ans	<- as.integer(anc.group[[i]][tmp2[1]])										
								ans
							})				
					list(IDX=IDX, CHANGE_NODE=unlist(tmp))
				},by=c('IND','MRCA','MRCA_IND')]
		zz			<- subset(zz, !is.na(CHANGE_NODE))
		#	each ancestor before a change node could have as one of its children a monophyletic clade from this patient		
		zz	<- unique(zz, by=c('MRCA_IND','CHANGE_NODE'))
		if(!nrow(zz))
		{
			set(zz, NULL, c('IDX','MRCA','IND'), NULL)
			zz[, CLADE:= rep(NA_integer_,0)]
		}			
		if(nrow(zz))
			zz	<- zz[, {
						#CHANGE_NODE<- 2053; MRCA<- 1580; MRCA_IND<- 'R1_RES827_S23_L001'
						#MRCA<- 1212; MRCA_IND<- 'R1_RES669_S20_L001'; CHANGE_NODE<- 1218
						#	the first two potentially monophyletic clades are the children of the CHANGE_NODE
						tmp			<- Children(ph, CHANGE_NODE)					
						#	define path from change node to just before mrca 					
						path		<- c(CHANGE_NODE, setdiff( Ancestors(ph,CHANGE_NODE), c(MRCA,Ancestors(ph,MRCA)) ))
						#	monophyletic clades could break off the path ending at MRCA
						#	the mrca's of these clades are the siblings of the path constructed above 
						pot.clade	<- c(tmp, unlist(Siblings(ph, path)))
						#	check if potential clades are monophyletic
						tmp			<- sapply(pot.clade, function(x){
									if(x<=Ntip(ph))
										ans	<- grepl(MRCA_IND, ph$tip.label[x])
									if(x>Ntip(ph))
										ans	<- all(grepl(MRCA_IND, phsc.extract.clade(ph,x)$tip.label))
									ans
								})	
						#	return monophyletic clades
						list(CLADE=as.integer(pot.clade[tmp]))
					}, by=c('MRCA_IND','CHANGE_NODE')]
		zz[, CHANGE_NODE:=NULL]
		setnames(zz,'MRCA_IND','IND')
		#	merge all monophyletic clades
		z	<- merge(z, unique(zz), all=1, by='IND', allow.cartesian=TRUE)
	}	
	tmp	<- z[, which(DIFF_IND==0)]
	set(z, tmp, 'CLADE', z[tmp,MRCA])
	z	<- subset(z, !is.na(CLADE))
	#	double check CLADEs (the MONOPH_PA condition is key)
	if(nrow(z))
	{
		tmp	<- subset(z, !is.na(CLADE))[,{
					if(CLADE<=Ntip(ph))
						tmp	<- grepl(IND, ph$tip.label[CLADE])
					if(CLADE>Ntip(ph))
						tmp	<- all(grepl(IND, phsc.extract.clade(ph,CLADE)$tip.label))
					list(MONOPH=tmp, MONOPH_PA=all(grepl(IND, phsc.extract.clade(ph,Ancestors(ph, CLADE, type="parent"))$tip.label)))
				}, by=c('IND','CLADE')]
		stopifnot(tmp[, all(MONOPH)], tmp[, !any(MONOPH_PA)])
		#
		#	collapse clades
		#
		zz	<- z[, {	
					#CLADE<- 292
					list(IDX=Descendants(ph, CLADE, type='tips')[[1]])			
				}, by='CLADE']
		zz	<- merge(phb, zz, by='IDX')
		zz	<- merge(zz, zz[, list(CLADE_N_TAXA=length(COUNT), CLADE_LABEL=paste0(IND[1],'_read_',CLADE[1],'_count_',sum(COUNT))), by='CLADE'],by='CLADE')
		zz	<- subset(zz, CLADE_N_TAXA>1)
		for(clu.id in zz[, unique(CLADE)])
		{
			tmp			<- subset(zz, CLADE==clu.id)[, BAM]			
			z			<- phsc.drop.tip(ph, tip=tmp, subtree=TRUE)
			if(!is.null(z))
			{
				ph		<- z
				ph$tip.label[ which(grepl("[",ph$tip.label,fixed=TRUE)) ]	<- subset(zz, CLADE==clu.id)[1,CLADE_LABEL]	
			}			 				
		}
	}		
	ph
}

phsc.drop.tip<- function (phy, tip, trim.internal = TRUE, subtree = FALSE, root.edge = 0, rooted = is.rooted(phy)) 
{
	#	code from ape
	if (!inherits(phy, "phylo")) 
		stop("object \"phy\" is not of class \"phylo\"")
	Ntip <- length(phy$tip.label)
	if (is.character(tip)) 
			tip <- which(phy$tip.label %in% tip)	
	out.of.range <- tip > Ntip
	if (any(out.of.range)) 
	{
		warning("some tip numbers were larger than the number of tips: they were ignored")
		tip <- tip[!out.of.range]
	}
	if (!length(tip)) 
		return(phy)
	if (length(tip) == Ntip) 
	{
		warning("drop all tips of the tree: returning NULL")
		return(NULL)
	}
	wbl <- !is.null(phy$edge.length)
	if (!rooted && subtree) 
	{
		phy <- root(phy, (1:Ntip)[-tip][1])
		root.edge <- 0
	}
	phy <- reorder(phy)
	NEWROOT <- ROOT <- Ntip + 1
	Nnode <- phy$Nnode
	Nedge <- dim(phy$edge)[1]
	if (subtree) 
	{
		trim.internal <- TRUE
		tr <- reorder(phy, "postorder")
		N <- .C(node_depth, as.integer(Ntip), as.integer(Nnode), 
				as.integer(tr$edge[, 1]), as.integer(tr$edge[, 2]), 
				as.integer(Nedge), double(Ntip + Nnode), 1L)[[6]]
	}
	edge1 <- phy$edge[, 1]
	edge2 <- phy$edge[, 2]
	keep <- !logical(Nedge)
	keep[match(tip, edge2)] <- FALSE
	if (trim.internal) 
	{
		ints <- edge2 > Ntip
		repeat 
		{
			sel <- !(edge2 %in% edge1[keep]) & ints & keep
			if (!sum(sel)) 
				break
			keep[sel] <- FALSE
		}
		if (subtree) 
		{
			subt <- edge1 %in% edge1[keep] & edge1 %in% edge1[!keep]
			keep[subt] <- TRUE
		}
		if (root.edge && wbl) 
		{
			degree <- tabulate(edge1[keep])
			if (degree[ROOT] == 1) 
			{
				j <- integer(0)
				repeat 
				{
					i <- which(edge1 == NEWROOT & keep)
					j <- c(i, j)
					NEWROOT <- edge2[i]
					degree <- tabulate(edge1[keep])
					if (degree[NEWROOT] > 1) 
						break
				}
				keep[j] <- FALSE
				if (length(j) > root.edge) 
					j <- 1:root.edge
				NewRootEdge <- sum(phy$edge.length[j])
				if (length(j) < root.edge && !is.null(phy$root.edge)) 
					NewRootEdge <- NewRootEdge + phy$root.edge
				phy$root.edge <- NewRootEdge
			}
		}
	}
	if (!root.edge) 
		phy$root.edge <- NULL
	phy$edge <- phy$edge[keep, ]
	if (wbl) 
		phy$edge.length <- phy$edge.length[keep]
	TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])	#all nodes that do not appear in first column are new tips
	oldNo.ofNewTips <- phy$edge[TERMS, 2]
	if (subtree) 
	{
		i <- which(tip %in% oldNo.ofNewTips)
		if (length(i)) 
		{
			phy$tip.label[tip[i]] <- "[1_tip]"
			tip <- tip[-i]
		}
	}
	n <- length(oldNo.ofNewTips)
	phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])	
	phy$tip.label <- phy$tip.label[-tip]
	if (subtree || !trim.internal) 
	{
		node2tip <- oldNo.ofNewTips[oldNo.ofNewTips > Ntip]
		new.tip.label <- if(subtree) 
				{
					paste("[", N[node2tip], "_tips]", sep = "")
				}
				else 
				{
					if (is.null(phy$node.label)) 
						rep("NA", length(node2tip))
					else phy$node.label[node2tip - Ntip]
				}
		phy$tip.label <- c(phy$tip.label, new.tip.label)
		#	the order of oldNo.ofNewTips does not mean anything. Need to sort.
		#	also: node2tip may not have a colour, whereas it should have the colour of one of the tips, so use 1 tip
		oldNo.ofNewTips	<- c( sort(oldNo.ofNewTips[-which(oldNo.ofNewTips > Ntip)]), tip[1] )
	}
	phy$Nnode <- dim(phy$edge)[1] - n + 1L
	newNb <- integer(Ntip + Nnode)
	newNb[NEWROOT] <- n + 1L
	sndcol <- phy$edge[, 2] > n
	newNb[sort(phy$edge[sndcol, 2])] <- (n + 2):(n + phy$Nnode)
	phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]]
	phy$edge[, 1] <- newNb[phy$edge[, 1]]
	storage.mode(phy$edge) <- "integer"
	if (!is.null(phy$node.label)) 
		phy$node.label <- phy$node.label[which(newNb > 0) - Ntip]	
	#	add new tip numbers to newNb
	#	want to reorder oldNo.ofNewTips so this corresponds to new tips
	newNb[oldNo.ofNewTips] 	<- seq_along(oldNo.ofNewTips)
	#	select remaining nodes and put into right order
	newNb					<- subset(data.table(IDX_O=seq_along(newNb), IDX_N=newNb), IDX_N>0)
	setkey(newNb, IDX_N)
	attr(phy, "INDIVIDUAL")	<- attr(phy, "INDIVIDUAL")[ newNb[, IDX_O] ]
	attr(phy, "NODE_SHAPES")<- attr(phy, "NODE_SHAPES")[ newNb[, IDX_O] ]
	attr(phy, "SPLIT")		<- attr(phy, "SPLIT")[ newNb[, IDX_O] ]
	#	collapse singles
	#	code from ape
	n 	<- length(phy$tip.label)
	e1 	<- phy$edge[, 1]
	e2 	<- phy$edge[, 2]
	tab <- tabulate(e1)
	if(all(tab > 1)) 
		return(phy)
	if(is.null(phy$edge.length)) 
	{
		root.edge <- FALSE
		wbl <- FALSE
	}
	if(!is.null(phy$edge.length)) 
	{
		wbl <- TRUE
		el <- phy$edge.length
	}
	if (root.edge) 
		ROOTEDGE <- 0
	ROOT <- n + 1L
	while (tab[ROOT] == 1) 
	{
		i <- which(e1 == ROOT)
		ROOT <- e2[i]
		if (wbl) {
			if (root.edge) 
				ROOTEDGE <- ROOTEDGE + el[i]
			el <- el[-i]
		}
		e1 <- e1[-i]
		e2 <- e2[-i]
	}
	singles <- which(tabulate(e1) == 1)
	while (length(singles)) 
	{
		i <- which(e1 == singles[1])
		j <- which(e2 == e1[i])
		e2[j] <- e2[i]
		if (wbl) {
			el[j] <- el[j] + el[i]
			el <- el[-i]
		}
		e1 <- e1[-i]
		e2 <- e2[-i]
		singles <- which(tabulate(e1) == 1)
	}
	Nnode <- length(e1) - n + 1L
	oldnodes <- unique(e1)
	if (!is.null(phy$node.label)) 
		phy$node.label <- phy$node.label[oldnodes - n]
	newNb <- integer(max(oldnodes))
	newNb[ROOT] <- n + 1L
	sndcol <- e2 > n
	e2[sndcol] <- newNb[e2[sndcol]] <- n + 2:Nnode
	e1 <- newNb[e1]
	phy$edge <- cbind(e1, e2, deparse.level = 0)
	phy$Nnode <- Nnode
	if (wbl) 
	{
		if (root.edge) 
			phy$root.edge <- ROOTEDGE
		phy$edge.length <- el
	}
	#	select remaining nodes and put into right order
	attr(phy, "INDIVIDUAL")		<- attr(phy, "INDIVIDUAL")[c(seq_len(Ntip(phy)), oldnodes)]
	attr(phy, "NODE_SHAPES")	<- attr(phy, "NODE_SHAPES")[c(seq_len(Ntip(phy)), oldnodes)]
	attr(phy, "SPLIT")			<- attr(phy, "SPLIT")[c(seq_len(Ntip(phy)), oldnodes)]		
	phy	
}

drop.tip2 <- function (phy, tip, trim.internal=TRUE, subtree=FALSE, root.edge=0, rooted=is.rooted(phy)) 
{
	if (!inherits(phy, "phylo")) 
		stop("object \"phy\" is not of class \"phylo\"")
	Ntip <- length(phy$tip.label)
	if (is.character(tip)) 
		tip <- which(phy$tip.label %in% tip)
	if (!rooted && subtree) {
		phy <- root(phy, (1:Ntip)[-tip][1])
		root.edge <- 0
	}
	phy <- reorder(phy)
	NEWROOT <- ROOT <- Ntip + 1
	Nnode <- phy$Nnode
	Nedge <- dim(phy$edge)[1]
	if (subtree) {
		trim.internal <- TRUE
		tr <- reorder(phy, "pruningwise")
		N <- .C("node_depth", as.integer(Ntip), as.integer(Nnode), 
				as.integer(tr$edge[, 1]), as.integer(tr$edge[, 2]), 
				as.integer(Nedge), double(Ntip + Nnode), DUP = FALSE, 
				PACKAGE = "ape")[[6]]
	}
	wbl <- !is.null(phy$edge.length)
	edge1 <- phy$edge[, 1]
	edge2 <- phy$edge[, 2]
	keep <- !logical(Nedge)
	if (is.character(tip)) 
		tip <- which(phy$tip.label %in% tip)
	if (!rooted && subtree) {
		phy <- root(phy, (1:Ntip)[-tip][1])
		root.edge <- 0
	}
	keep[match(tip, edge2)] <- FALSE
	if (trim.internal) {
		ints <- edge2 > Ntip
		repeat {
			sel <- !(edge2 %in% edge1[keep]) & ints & keep
			if (!sum(sel)) 
				break
			keep[sel] <- FALSE
		}
		if (subtree) {
			subt <- edge1 %in% edge1[keep] & edge1 %in% edge1[!keep]
			keep[subt] <- TRUE
		}
		if (root.edge && wbl) {
			degree <- tabulate(edge1[keep])
			if (degree[ROOT] == 1) {
				j <- integer(0)
				repeat {
					i <- which(edge1 == NEWROOT & keep)
					j <- c(i, j)
					NEWROOT <- edge2[i]
					degree <- tabulate(edge1[keep])
					if (degree[NEWROOT] > 1) 
						break
				}
				keep[j] <- FALSE
				if (length(j) > root.edge) 
					j <- 1:root.edge
				NewRootEdge <- sum(phy$edge.length[j])
				if (length(j) < root.edge && !is.null(phy$root.edge)) 
					NewRootEdge <- NewRootEdge + phy$root.edge
				phy$root.edge <- NewRootEdge
			}
		}
	}
	if (!root.edge) 
		phy$root.edge <- NULL
	phy$edge <- phy$edge[keep, ]
	if (wbl) 
		phy$edge.length <- phy$edge.length[keep]
	TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])
	oldNo.ofNewTips <- phy$edge[TERMS, 2]
	n <- length(oldNo.ofNewTips)
	phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
	if (subtree || !trim.internal) {
		tips.kept <- oldNo.ofNewTips <= Ntip & !(oldNo.ofNewTips %in% 
					tip)
		new.tip.label <- character(n)
		new.tip.label[tips.kept] <- phy$tip.label[-tip]
		node2tip <- oldNo.ofNewTips[!tips.kept]
		new.tip.label[!tips.kept] <- if (subtree) {
					paste("[", N[node2tip], "_tips]", sep = "")
				}
				else {
					if (is.null(phy$node.label)) 
						rep("NA", length(node2tip))
					else phy$node.label[node2tip - Ntip]
				}
		if (!is.null(phy$node.label)) 
			phy$node.label <- phy$node.label[-(node2tip - Ntip)]
		phy$tip.label <- new.tip.label
	}
	if (!(subtree || !trim.internal))
		phy$tip.label <- phy$tip.label[-tip]
	if (!is.null(phy$node.label)) 
		phy$node.label <- phy$node.label[sort(unique(phy$edge[, 
										1])) - Ntip]
	phy$Nnode <- dim(phy$edge)[1] - n + 1L
	newNb <- integer(n + phy$Nnode)
	newNb[NEWROOT] <- n + 1L
	sndcol <- phy$edge[, 2] > n
	phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]] <- (n + 
				2):(n + phy$Nnode)
	phy$edge[, 1] <- newNb[phy$edge[, 1]]
	storage.mode(phy$edge) <- "integer"
	phy <- collapse.singles(phy)	
	# handle non-standard node attributes
	drop <- edge1[which(!keep)]
	z	<- attr(phy, "INDIVIDUAL")
	z[-drop]
	attr(phy, "INDIVIDUAL")[-drop]

	sdtnames <- c("edge", "edge.length", "Nnode", "tip.label", 		"node.label")
	if (!all(names(phy) %in% sdtnames)){
		id <- which(!names(phy) %in% sdtnames)
		drop <- edge1[which(!keep)] - Ntip
		for (i in id){
			phy[[i]] <- phy[[i]][-drop]
		}
	}	
	phy    
}

#' @export
#' @import data.table grid ggtree
#' @title Plot short read phylogenies and highlight individuals
#' @description This function plots short read phylogenies and highlights the clades of two individuals in red and blue.  
#' @param phs List of trees in ape format
#' @param dfs data.table with mandatory column 'IDX' and optional column 'TITLE'. IDX is the index of all phylogenies in 'phs' that are to be plotted. TITLE is a title for each sub-plot, for example specifying a window.
#' @param ids Vector of regular expressions that identify individuals to be highlighted in colour.
#' @param plot.cols Vector of colours for each individual
#' @param group.redo Logical, indicating if the colour groups should be recalculated from ids.
#' @param drop.blacklisted Logical, indicating if all blacklisted taxa should be dropped prior to plotting.  
#' @param pdf.h	Height of the pdf file in inches.
#' @param pdf.rw Relative width of the pdf file, internally multiplied by the number of phylogenies to give the total width in inches.
#' @param pdf.ntrees Number of trees per pdf.
#' @param pdf.title.size Size of pdf title in inches.
#' @param plot.file If not missing, the phylogenies will be printed to file.	
#' @return List of ggtree objects, ready for printing.
phsc.plot.phycollapsed.selected.individuals<- function(phs, dfs, ids, plot.cols=rainbow(length(ids)), plot.cols.background=function(n){ rainbow_hcl(n, start = 60, end = 240, c=30, l=80) }, drop.blacklisted=TRUE, tip.regex='(.*)_read_([0-9]+)_count_([0-9]+)', plot.file=NA, drop.less.than.n.ids=2, pdf.h=50, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40)
{
	require(colorspace)
	suppressMessages(require(grid))	
	#	determine which phylogenies contain at least one of the requested individuals
	tmp		<- copy(dfs)
	tmp		<- merge(tmp, tmp[, {
						ph	<- phs[[ IDX ]]
						z	<- ph$tip.label
						if(drop.blacklisted)
							z	<- as.character(attr(ph, "INDIVIDUAL"))
						list(HAS_N_IND=length(which(sapply(ids, function(id) any(grepl(id, z)) ))))						  				
					}, by='IDX'], by='IDX')
	tmp		<- subset(tmp, HAS_N_IND>=drop.less.than.n.ids)
	#	define colours so that they don t change across windows
	cols	<- tmp[, {
				ph	<- phs[[ IDX ]]
				phb			<- data.table(	TAXA=ph$tip.label,
											IDX=seq_along(ph$tip.label), 											 
											ID= gsub(tip.regex,'\\1',ph$tip.label),
											COUNT=sub(tip.regex,'\\3',ph$tip.label)		)
				set(phb, phb[, which(ID==COUNT)], 'COUNT', NA_character_)
				set(phb, NULL, 'COUNT', phb[, as.integer(COUNT)])
				set(phb, phb[, which(is.na(COUNT))], 'ID','REFERENCE')	
				list(ID= phb[, unique(ID)])				
			}, by='IDX']
	cols	<- unique(subset(cols, select=ID))
	cols[,COL:= plot.cols.background(nrow(cols))]	
	set(cols, cols[,which(grepl('REFERENCE',ID))], 'COL', 'grey50')
	for(i in seq_along(ids))
		set(cols, cols[, which(grepl(ids[i], ID))], 'COL', plot.cols[i])
	#	
	phps	<- lapply(seq_len(nrow(tmp)), function(i){
				cat(i,'\n')
				ph.title	<- NULL
				if('TITLE'%in%colnames(tmp))
					ph.title	<- tmp[i, TITLE]										
				ph			<- phs[[ tmp[i, IDX] ]]
				ph			<- phsc.phy.collapse.monophyletic.clades(ph, drop.blacklisted=TRUE, tip.regex=tip.regex)
				#
				#	define IDs, find COUNTS, find references
				phb			<- data.table(	TAXA=ph$tip.label,
											IDX=seq_along(ph$tip.label), 											 
											ID= gsub(tip.regex,'\\1',ph$tip.label),
											COUNT=sub(tip.regex,'\\3',ph$tip.label)		)
				set(phb, phb[, which(ID==COUNT)], 'COUNT', NA_character_)
				set(phb, NULL, 'COUNT', phb[, as.integer(COUNT)])
				set(phb, phb[, which(is.na(COUNT))], 'ID','REFERENCE')	
				set(phb, phb[, which(is.na(COUNT))], 'COUNT', 1L)
				#	define cols for this tree
				tmp					<- data.table(IDX= seq_along(attr(ph, "INDIVIDUAL")), ID= attr(ph, "INDIVIDUAL"))
				set(tmp, tmp[, which(is.na(ID))], 'ID', 'REFERENCE')
				tmp					<- merge(tmp, cols, by='ID')
				attr(ph, 'COLOUR')	<- tmp[order(IDX),][,COL]	
				#						
				p 			<- ggtree(ph, aes(color=I(COLOUR))) %<+% phb +
						geom_tippoint(aes(size=COUNT), shape=18) +
						#geom_point2(shape=16, aes(size=NODE_S, subset=NODE_S)) +
						scale_fill_hue(na.value="black") +								
						theme(legend.position="none") +
						geom_tiplab(align=T, aes(col=I(COLOUR)), size=1, linetype=NA, linesize=NA) +
						guides(size='none') +
						theme_tree2() +
						theme(	axis.line.x=element_line(),
								panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3),
								legend.position="bottom",								
								plot.title = element_text(hjust = 0.5, size=pdf.title.size)) + 
						ggplot2::xlim(0, max(node.depth.edgelength(ph)[1:Ntip(ph)])*1.05) +
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

#' @export
#' @import data.table grid ggtree
#' @title Plot short read phylogenies and highlight individuals
#' @description This function plots short read phylogenies and highlights the clades of two individuals in red and blue.  
#' @param phs List of trees in ape format
#' @param dfs data.table with mandatory column 'IDX' and optional column 'TITLE'. IDX is the index of all phylogenies in 'phs' that are to be plotted. TITLE is a title for each sub-plot, for example specifying a window.
#' @param ids Vector of regular expressions that identify individuals to be highlighted in colour.
#' @param plot.cols Vector of colours for each individual
#' @param group.redo Logical, indicating if the colour groups should be recalculated from ids.
#' @param drop.blacklisted Logical, indicating if all blacklisted taxa should be dropped prior to plotting.  
#' @param pdf.h	Height of the pdf file in inches.
#' @param pdf.rw Relative width of the pdf file, internally multiplied by the number of phylogenies to give the total width in inches.
#' @param pdf.ntrees Number of trees per pdf.
#' @param pdf.title.size Size of pdf title in inches.
#' @param plot.file If not missing, the phylogenies will be printed to file.	
#' @return List of ggtree objects, ready for printing.
phsc.plot.phy.selected.individuals<- function(phs, dfs, ids, plot.cols=rainbow(length(ids)), plot.file=NA, group.redo=FALSE, drop.less.than.n.ids=2, drop.blacklisted=FALSE, pdf.h=50, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40)
{
	suppressMessages(require(grid))
	#	determine which phylogenies contain at least one of the requested individuals
	tmp		<- copy(dfs)
	tmp		<- merge(tmp, tmp[, {
						ph	<- phs[[ IDX ]]
						z	<- ph$tip.label
						if(drop.blacklisted)
							z	<- as.character(attr(ph, "INDIVIDUAL"))
						list(HAS_N_IND=length(which(sapply(ids, function(id) any(grepl(id, z)) ))))						  				
					}, by='IDX'], by='IDX')
	tmp		<- subset(tmp, HAS_N_IND>=drop.less.than.n.ids)
	phps	<- lapply(seq_len(nrow(tmp)), function(i){
				ph.title	<- NULL
				if('TITLE'%in%colnames(tmp))
					ph.title	<- tmp[i, TITLE]										
				ph			<- phs[[ tmp[i, IDX] ]]
				#
				#
				if(drop.blacklisted)
				{
					z						<- data.table(FROM=ph$edge[,1],TO=ph$edge[,2])
					z						<- merge(z, data.table(TO= seq_len(length(attr(ph,'INDIVIDUAL'))), GROUP= attr(ph,'INDIVIDUAL')), by='TO')
					z						<- subset(z, is.na(GROUP) & TO<=Ntip(ph))[, TO]
					ph						<- drop.tip(ph, z)
					attr(ph,'NODE_SHAPES')	<- rep(FALSE, Nnode(ph, internal.only=FALSE))
					attr(ph,'INDIVIDUAL')	<- NULL
					attr(ph,'SPLIT')		<- NULL
				}
				#
				#
				if(group.redo||drop.blacklisted)
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
						#geom_point2(shape = 16, size=3, aes(subset=NODE_SHAPES)) +
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
#' @param save.file If not missing, function output (a data.table) will be stored to 'save.file', input files will be zipped, and input files will be deleted.   	
#' @param resume If TRUE and save.file is not missing, the function loads and returns trees stored in save.file.
#' @return Data table with columns PAIR_ID, ID1, ID2, TYPE, WIN_OF_TYPE, PTY_RUN, WIN_TOTAL, SCORE 
phsc.read.likelytransmissions<- function(prefix.infiles, prefix.run='ptyr', regexpr.lklsu='trmStats.csv$', save.file=NA, resume=FALSE, zip=FALSE)
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