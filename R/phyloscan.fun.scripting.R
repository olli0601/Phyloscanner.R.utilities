#' @export
#' @import data.table
#' @title Generate bash commands to calculate read distribution in bam file
#' @description This function generates bash commands to calculate read distribution in bam file   
#' @param pty.runs data.table with columns SAMPLE_ID and PTY_RUN
#' @param pty.args List of input variables containing the fields "prog.bam.distr.calculator", "data.dir", "out.dir", "work.dir"
#' @return character string of bash commands.
phsc.cmd.bam.calculate.read.distribution <- function(pty.runs, pty.args) 		
{
	#
	#	associate BAM and REF files with each scheduled phylotype run
	#	
	set(pty.runs, NULL, 'SAMPLE_ID',  pty.runs[, gsub('\\.bam$','',SAMPLE_ID)])
	#	get available Bam files
	ptyd		<- data.table(FILE=list.files(pty.args[['data.dir']], full.names=TRUE))
	ptyd[, TYPE:=NA_character_]
	set(ptyd, ptyd[, which(grepl('.bam$',FILE))], 'TYPE', 'BAM')
	#	get Reference files
	#	if reference files that were used in assembly are not specified in pty.runs, search for SAMPLE_ID+'_ref.fasta'
	if(!any(colnames(pty.runs)=='REFERENCE_ID'))
	{		
		set(ptyd, ptyd[, which(grepl('_ref.fasta$',FILE))], 'TYPE', 'REF')		
		ptyd		<- subset(ptyd, !is.na(TYPE))
		ptyd[, SAMPLE_ID:= gsub('\\.bam|_ref\\.fasta','',basename(FILE))]
		ptyd		<- dcast.data.table(ptyd, SAMPLE_ID~TYPE, value.var='FILE')		
	}
	#	if reference files that were used in assembly are specified in pty.runs, use these
	if(any(colnames(pty.runs)=='REFERENCE_ID'))
	{
		tmp			<- subset(ptyd, is.na(TYPE))
		tmp[, REFERENCE_ID:= gsub('\\.bam','',basename(FILE))]
		tmp			<- merge(subset(pty.runs, select=c(SAMPLE_ID, REFERENCE_ID)), subset(tmp, select=c(REFERENCE_ID, FILE)), by='REFERENCE_ID')
		tmp[, TYPE:='REF']
		tmp[, REFERENCE_ID:=NULL]
		ptyd		<- subset(ptyd, !is.na(TYPE))
		ptyd[, SAMPLE_ID:= gsub('\\.bam$','',basename(FILE))]
		ptyd		<- rbind(ptyd, tmp)
		ptyd		<- dcast.data.table(ptyd, SAMPLE_ID~TYPE, value.var='FILE')
	}
	#	merge
	ptyd	<- merge(pty.runs, ptyd, by='SAMPLE_ID', all.x=1)	
	if(ptyd[,any(is.na(BAM))])
		warning('\nCould not find location of BAM files for all individuals in pty.runs, n=', ptyd[, length(which(is.na(BAM)))],'\nMissing individuals are ignored. Please check.')	
	if(ptyd[,any(is.na(REF))])
		warning('\nCould not find location of reference files for all individuals in pty.runs, n=', ptyd[, length(which(is.na(REF)))],'\nMissing individuals are ignored. Please check.')
	ptyd	<- subset(ptyd, !is.na(BAM) & !is.na(REF))
	setkey(ptyd, PTY_RUN)
	#
	#	write pty.run files and get pty command lines
	#
	pty.c		<- ptyd[, {
				#	PTY_RUN<- z <- 1; BAM<- subset(ptyd, PTY_RUN==z)[, BAM]; REF<- subset(ptyd, PTY_RUN==z)[, REF]				
				file.input		<- file.path(pty.args[['work.dir']], paste('bamr',PTY_RUN,'_input.csv',sep=''))
				file.output		<- file.path(pty.args[['out.dir']], paste('bamr',PTY_RUN,'_read_distributions.csv',sep=''))
				tmp				<- cbind(BAM[!is.na(BAM)&!is.na(REF)], REF[!is.na(BAM)&!is.na(REF)])
				write.table(tmp, file=file.input, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')
				cmd				<- paste0( 	pty.args[['prog.bam.distr.calculator']],
						' "',file.input,'"',
						' --out-filename "',file.output,'"',
						' --overlapping-insert-sizes',
						' --dont-plot')				
				list(CMD= cmd)				
			},by='PTY_RUN']
	pty.c
}	

#' @export
#' @title Generate bash commands to generate phylogenies with RAxML  
#' @description This function generates bash commands to generate one phylogeny with RAxML, that can be called via 'system' in R, or written to file to run on a UNIX system.
#' @param infile.fasta Full path name to input fasta file.
#' @param outfile Full path name to output file. All RAxML output other than the best tree will be zipped and returned as 'outfile.zip'.
#' @param pr Full path name to RAXML program.
#' @param pr.args RAXML arguments. 
#' @return	Character string
raxml.cmd<- function(infile.fasta, outfile=paste(infile.fasta,'.newick',sep=''), pr='raxmlHPC-SSE3', pr.args='-m GTRCAT --HKY85 -p 42')
{		
	cmd				<- paste("#######################################################
# start: RAXML
#######################################################\n",sep='')												
	cmd				<- paste(cmd,"CWD=$(pwd)\n",sep='')
	cmd				<- paste(cmd,"echo $CWD\n",sep='')	
	tmpdir.prefix	<- paste('rx_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
	tmpdir			<- paste("$CWD/",tmpdir.prefix,sep='')
	tmp.in			<- basename(infile.fasta)
	tmp.out			<- basename(outfile)	
	cmd				<- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
	cmd				<- paste(cmd,'cp "',infile.fasta,'" ',file.path(tmpdir,tmp.in),'\n', sep='')	
	cmd				<- paste(cmd,'cd "',tmpdir,'"\n', sep='')	
	cmd				<- paste(cmd, pr,' ',pr.args,' -s ', tmp.in,' -n ', tmp.out,'\n', sep='')
	cmd				<- paste(cmd, "rm ", tmp.in,'\n',sep='')	
	cmd				<- paste(cmd, 'cp RAxML_bestTree.',basename(outfile),' "',outfile,'"\n',sep='')
	cmd				<- paste(cmd, 'for file in *; do\n\tzip -ur9XTjq ',basename(outfile),'.zip "$file"\ndone\n',sep='')	
	cmd				<- paste(cmd, 'cp ',basename(outfile),'.zip "',dirname(outfile),'"\n',sep='')
	cmd				<- paste(cmd,'cd $CWD\n', sep='')
	cmd				<- paste(cmd, "rm -r ", tmpdir,'\n',sep='')
	cmd				<- paste(cmd, "#######################################################
# end: RAXML
#######################################################\n",sep='')
	cmd
}

#' @export
#' @title Generate bash commands for multiple phyloscanner runs
#' @param pty.runs Data.table of individual assignments to phyloscanner runs, with columns 'PTY_RUN' (run id), 'SAMPLE_ID' (ID of individuals that are assigned to that run). Optional columns: 'RENAME_ID' (new ID for each bam file in phyloscanner output).
#' @param pty.args List of phyloscanner input variables. See examples.
#' @return Data.table with columns 'PTY_RUN' (run id) and 'CMD' (bash commands for that run). 
#' @description This function generates bash commands for multiple phyloscanner runs, that can be called via 'system' in R, or written to file to run on a UNIX system.
#' @example example/ex.cmd.phyloscanner.multi.R  
phsc.cmd.phyloscanner.multi <- function(pty.runs, pty.args) 		
{
	stopifnot( any(colnames(pty.runs)=='PTY_RUN') )
	stopifnot( any(colnames(pty.runs)=='SAMPLE_ID') )
	stopifnot( any(colnames(pty.runs)=='BAM') )
	stopifnot( any(colnames(pty.runs)=='REF') )	
	ptyd	<- copy(pty.runs)		
	if(!any(is.na(pty.args[['select']])))
		ptyd<- subset(ptyd, PTY_RUN%in%pty.args[['select']])			
	#	if background alignment is specified in pty.runs, use it
	if(any(colnames(ptyd)=='BACKGROUND_ID'))
	{
		tmp	<- unique(data.table(BACKGROUND_ID=ptyd[, BACKGROUND_ID], FILE=file.path(pty.args[['data.dir']], ptyd[, BACKGROUND_ID])))
		#	check if files exist
		tmp	<- tmp[, list(EXISTS=file.exists(FILE)), by=c('BACKGROUND_ID','FILE')]
		if(any(tmp[, !EXISTS]))
			warning('\nCould not find location of BACKGROUND files for all runs in pty.runs, n=', tmp[, length(which(!EXISTS))],'\nRuns with missing background alignments are ignored. Please check.')
		ptyd <- merge(ptyd, subset(tmp, EXISTS, c(BACKGROUND_ID, FILE)), by='BACKGROUND_ID')
		set(ptyd, NULL, 'BACKGROUND_ID', NULL)
		setnames(ptyd, 'FILE', 'BACKGROUND_ID')
	}
	setkey(ptyd, PTY_RUN)
	#	add legacy arguments
	pty.args['all.bootstrap.trees']	<- FALSE
	pty.args['num.bootstraps']	<- 0		
	#	write pty.run files and get pty command lines
	pty.c		<- ptyd[, {
				#	PTY_RUN<- z <- 1; BAM<- subset(ptyd, PTY_RUN==z)[, BAM]; REF<- subset(ptyd, PTY_RUN==z)[, REF]
				#	SAMPLE_ID<- subset(ptyd, PTY_RUN==z)[, SAMPLE_ID]; RENAME_ID<- subset(ptyd, PTY_RUN==z)[, RENAME_ID]; BACKGROUND_ID<- subset(ptyd, PTY_RUN==z)[, BACKGROUND_ID]
				file.input		<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_input.csv',sep=''))
				tmp				<- cbind(BAM[!is.na(BAM)&!is.na(REF)], REF[!is.na(BAM)&!is.na(REF)])
				if(exists('RENAME_ID'))
					tmp			<- cbind(tmp, RENAME_ID[!is.na(BAM)&!is.na(REF)])
				write.table(tmp, file=file.input, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')
				file.patient	<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_patients.txt',sep=''))
				tmp				<- paste0(SAMPLE_ID[!is.na(BAM)&!is.na(REF)],'.bam')
				if(exists('RENAME_ID'))
					tmp			<- unique(RENAME_ID[!is.na(BAM)&!is.na(REF)])
				if(exists('UNIT_ID'))
					tmp			<- unique(UNIT_ID[!is.na(BAM)&!is.na(REF)])				
				cat( paste(tmp,collapse='\n'), file= file.patient	)
				if(exists('BACKGROUND_ID'))
				{
					if(length(unique(BACKGROUND_ID))!=1)
						stop('\nrun',PTY_RUN,' Expected one background alignment file, found: ',paste(unique(BACKGROUND_ID), collapse=''),'. Please check.')
					pty.args$alignments.file<- unique(BACKGROUND_ID)
				}				
				cmd			<- phsc.cmd.phyloscanner.one(	pty.args, 
						file.input, 
						file.patient)
				#cmd			<- paste(cmd, pty.cmd.evaluate.fasta(pty.args[['out.dir']], strip.max.len=pty.args[['strip.max.len']], select=paste('^ptyr',PTY_RUN,'_In',sep=''), min.ureads.individual=pty.args[['min.ureads.individual']]), sep='')
				#cat(cmd)
				list(CMD= cmd)				
			},by='PTY_RUN']
	pty.c
}	

#' @export
#' @title Generate bash command for a single phyloscanner run
#' @param pty.args List of phyloscanner input variables. See examples.
#' @param file.input File name of the file that contains the list of bam files, reference files, and potentially aliases
#' @param file.patient File name of the file that contains the list of unique individuals/units to which the bam files correspond. Multiple bam files are allowed per individual/unit. 
#' @return Character string of phyloscanner commands.
#' @description This function generates bash commands for a single phyloscanner run, that can be called via 'system' in R, or written to file to run on a UNIX system.
#' @example example/ex.cmd.phyloscanner.one.R    
phsc.cmd.phyloscanner.one<- function(pty.args, file.input, file.patient)
{	
	stopifnot(is.character(file.input),is.character(file.patient))
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
	dont.check.recombination	<- ifelse(!is.na(dont.check.recombination) & dont.check.recombination==FALSE, '--check-recombination', NA_character_)	
	#	create local tmp dir
	cmd		<- paste("CWD=$(pwd)\n",sep='\n')
	cmd		<- paste(cmd,"echo $CWD\n",sep='')
	tmpdir	<- paste('pty','_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')	
	tmpdir	<- paste("$CWD/",tmpdir,sep='')
	cmd		<- paste(cmd,'mkdir -p "',tmpdir,'"\n',sep='')
	#	copy files to local tmp dir
	cmd		<- paste(cmd,'cp "',file.input,'" "',tmpdir,'"\n',sep='')	
	cmd		<- paste(cmd,'cp "',file.patient,'" "',tmpdir,'"\n',sep='')
	#	cd to tmp dir
	cmd		<- paste(cmd, 'cd "',tmpdir,'"\n', sep='')	
	cmd		<- paste(cmd, prog.pty,' "',basename(file.input),'" ',sep='')	
	cmd		<- paste(cmd, '--merging-threshold-a', merge.threshold,'--min-read-count',min.read.count,'--quality-trim-ends', quality.trim.ends, '--min-internal-quality',min.internal.quality,'--keep-output-together')
	if(!is.na(merge.paired.reads))
		cmd	<- paste(cmd,' ',merge.paired.reads,sep='')	
	if(!is.na(dont.check.duplicates))
		cmd	<- paste(cmd,' ',dont.check.duplicates,sep='')
	if(!is.na(dont.check.recombination))
		cmd	<- paste(cmd,' ',dont.check.recombination,sep='')	
	if(!is.na(alignments.pairwise.to))		
		cmd	<- paste(cmd,' --pairwise-align-to ',alignments.pairwise.to,sep='')
	if(nchar(window.automatic))
		cmd	<- paste(cmd,' --auto-window-params ', window.automatic,sep='')
	if(!nchar(window.automatic))
		cmd	<- paste(cmd,' --windows ', paste(as.character(window.coord), collapse=','),sep='')
	if(!is.na(alignments.file))
		cmd	<- paste(cmd,' --alignment-of-other-refs "',alignments.file,'"',sep='')	
	if(!is.na(no.trees))		
		cmd	<- paste(cmd, no.trees)
	if(!is.na(num.bootstraps))
		cmd	<- paste(cmd, ' --num-bootstraps ', num.bootstraps, sep='')
	if(!is.na(keep.overhangs))
		cmd	<- paste(cmd, keep.overhangs)
	cmd		<- paste(cmd, '--x-mafft',prog.mafft)
	if(is.na(no.trees))
		cmd	<- paste(cmd, '--x-raxml',prog.raxml)	
	cmd		<- paste(cmd, '\n')
	run.id	<- gsub('_input.csv','',basename(file.input))
	#	process RAxML files
	if(is.na(no.trees) & (is.na(num.bootstraps) | (!is.na(num.bootstraps) & all.bootstrap.trees)))
		cmd	<- paste(cmd, 'for file in RAxML_bestTree*.tree; do\n\tmv "$file" "${file//RAxML_bestTree\\./',run.id,'_}"\ndone\n',sep='')
	if(is.na(no.trees) & !is.na(num.bootstraps) & !all.bootstrap.trees)
		cmd	<- paste(cmd, 'for file in RAxML_bipartitions.MLtreeWbootstraps*.tree; do\n\tmv "$file" "${file//RAxML_bipartitions.MLtreeWbootstraps/',run.id,'_}"\ndone\n',sep='')	
	#cmd	<- paste(cmd, "for file in AlignedReads*.fasta; do\n\tsed 's/<unknown description>//' \"$file\" > \"$file\".sed\n\tmv \"$file\".sed \"$file\"\ndone\n",sep='')		
	if(!is.na(alignments.file) & !is.na(keep.overhangs))
	{
		cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\tcat "$file" | awk \'{if (substr($0,1,4) == ">REF") censor=1; else if (substr($0,1,1) == ">") censor=0; if (censor==0) print $0}\' > NoRef$file\ndone\n', sep='')		
		cmd	<- paste(cmd, 'for file in NoRefAlignedReads*.fasta; do\n\t',phsc.cmd.mafft.add(alignments.file,'"$file"','Ref"$file"', options='--keeplength --memsave --parttree --retree 1'),'\ndone\n',sep='')		
		cmd	<- paste(cmd, 'for file in RefNoRefAlignedReads*.fasta; do\n\t','mv "$file" "${file//RefNoRefAlignedReads/',run.id,'_}"\ndone\n',sep='')		
	}
	if(is.na(alignments.file) || is.na(keep.overhangs))
	{
		cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\tmv "$file" "${file//AlignedReads/',run.id,'_}"\ndone\n',sep='')	
	}	
	#	move Duplicate Read Counts - only for backward compatibility
	#cmd	<- paste(cmd, 'for file in DuplicateReadCountsProcessed_*.csv; do\n\tmv "$file" "${file//DuplicateReadCountsProcessed_/',run.id,'_DuplicateReadCounts_}"\ndone',sep='')
	#	run phyloscanner tools and compress output
	if(is.na(no.trees))
		cmd	<- paste(cmd, phsc.cmd.process.phyloscanner.output.in.directory(tmpdir, file.patient, pty.args), sep='\n')
	#
	out.dir2<- out.dir	
	if(!is.na(no.trees))
	{
		out.dir2<- file.path(out.dir,paste0(run.id,'_trees'))
		cmd		<- paste(cmd,'\nmkdir -p ',out.dir2)		
	}
	cmd		<- paste(cmd, '\ncp ',run.id,'* "',out.dir2,'"\n',sep='')	
	#	zip up everything else
	tmp		<- ''
	if(length(window.coord)==2)
		tmp	<- window.coord[1]
	if(is.null(mem.save) || is.na(mem.save) || mem.save==0)
	{
		cmd		<- paste(cmd, 'for file in *; do\n\tzip -ur9XTjq ',paste(run.id,'_otherstuff',tmp,'.zip',sep=''),' "$file"\ndone\n',sep='')
		cmd		<- paste(cmd, 'cp ',paste(run.id,'_otherstuff',tmp,'.zip',sep=''),' "',out.dir2,'"\n',sep='')		
	}
	#	clean up
	cmd		<- paste(cmd,'cd $CWD\nrm -r "',tmpdir,'"\n',sep='')
	detach(pty.args)
	cmd
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
	#tmp		<- list.files(dirname(prefix.infiles), pattern=paste('^',basename(prefix.infiles),'.*fasta.zip',sep=''), full.names=TRUE)
	#stopifnot(length(tmp)==1)
	#cmd		<- paste(cmd,'unzip "',tmp,'" -d "',tmpdir,'"\n',sep='')
	tmp		<- list.files(dirname(prefix.infiles), pattern=paste('^',basename(prefix.infiles),'.*newick.zip',sep=''), full.names=TRUE)
	stopifnot(length(tmp)==1)
	cmd		<- paste(cmd,'unzip "',tmp,'" -d "',tmpdir,'"\n',sep='')
	#	cd to tmp dir
	cmd		<- paste(cmd, 'cd "',tmpdir,'"\n', sep='')	
	#	add all toolkit commands according to pty.args
	cmd		<- paste(cmd, phsc.cmd.process.phyloscanner.output.in.directory(tmpdir, file.patient, pty.args), collapse='\n',sep='')
	#	move all files starting with current run ID
	run.id	<- gsub('_patients.txt','',basename(file.patient))
	cmd		<- paste(cmd, '\ncp ',run.id,'* "',pty.args$out.dir,'"\n',sep='')	
	#	zip up everything else	
	tmp	<- paste(run.id,'_otherstuff.zip',sep='')
	cmd		<- paste(cmd, 'for file in *; do\n\tzip -ur9XTj ',tmp,' "$file"\ndone\n',sep='')
	cmd		<- paste(cmd, 'cp ',tmp,' "',pty.args$out.dir,'"\n',sep='')
	#	clean up
	cmd		<- paste(cmd,'cd $CWD\nrm -r "',tmpdir,'"\n',sep='')		
	cmd
}

