PR.PACKAGE			<- "phyloscan" 
PR.read.processed.phyloscanner.output.in.directory		<- paste('Rscript', system.file(package=PR.PACKAGE, "phsc.read.processed.phyloscanner.output.in.directory.Rscript"))


#' @export
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

#' @export
phsc.cmd.SummaryStatistics<- function(pr, scriptdir, outgroupName, file.patients, treeFiles, fastaFiles, splitsFiles, outputBaseName)
{
	paste(pr,' --scriptdir ',scriptdir,' --outgroupName ', outgroupName, ' "', file.patients, '" "', treeFiles, '" "',fastaFiles, '" "',splitsFiles, '" "',outputBaseName, '"', sep='')
}


#' @export
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

#' @export
phsc.cmd.SplitPatientsToSubtrees<- function(pr, scriptdir, infile, outputdir=NA, outputFileIdentifier=NA, outgroupName=NA, pdfwidth=30, pdfrelheight=0.15)	
{
	cmd	<- paste(pr, ' "', infile, '"', ' --scriptdir "', scriptdir,'"', sep='')
	if(!is.na(outputdir))
		cmd	<- paste(cmd, ' --outputdir "', outputdir,'"', sep='')	
	if(!is.na(outputFileIdentifier))
		cmd	<- paste(cmd, ' --outputfileid "', outputFileIdentifier, '"', sep='')
	if(!is.na(outgroupName))
		cmd	<- paste(cmd, ' --outgroupName ', outgroupName, sep='')
	if(!is.na(pdfwidth))
		cmd	<- paste(cmd, ' --pdfwidth ', pdfwidth, sep='')
	if(!is.na(pdfrelheight))
		cmd	<- paste(cmd, ' --pdfrelheight ', pdfrelheight, sep='')
	cmd	
}
	
#' @export
phsc.cmd.LikelyTransmissions<- function(pr, scriptdir, file.tree, file.splits, file.out, root.name=NA, zeroLengthTipsCount=FALSE, romeroSeverson=FALSE, dual.inf.thr=NA)
{
	cmd<- paste(pr, ' "', file.tree, '" "', file.splits, '" "', file.out,'" ',' --scriptdir ', scriptdir, ' --outgroupName ',root.name, sep='')
	if(!is.na(dual.inf.thr))
		cmd	<- paste(cmd, '--dualInfectionThreshold', dual.inf.thr)
	if(zeroLengthTipsCount)
		cmd	<- paste(cmd, '--zeroLengthTipsCount')
	if(romeroSeverson)
		cmd	<- paste(cmd, '--romeroSeverson')
	cmd
}

#' @export
phsc.cmd.LikelyTransmissionsSummary<- function(pr, scriptdir, file.patients, file.summary, file.lkl, file.out, min.threshold=1, allow.splits=FALSE)
{
	cmd<- paste(pr, ' "',file.patients, '" "', file.lkl, '" "', file.out,'" --scriptdir ', scriptdir, ' --summaryFile "', file.summary, '" --minThreshold ', min.threshold, sep='')
	if(allow.splits)
		cmd	<- paste(cmd, ' --allowSplits')	
	cmd
}

#' @export
#' @title Generate bash command for a single phyloscanner run
#' @description This function generates bash commands for a single phyloscanner run, that can be called via 'system' in R, or written to file to run on a UNIX system.  
phsc.cmd.phyloscanner.one<- function(file.bam, file.ref, window.coord=integer(0), window.automatic='', prog=PROG.PTY, prog.raxml='raxmlHPC-AVX', prog.mafft='mafft', merge.threshold=1, min.read.count=2, quality.trim.ends=30, min.internal.quality=2, merge.paired.reads=TRUE, num.bootstraps=1, all.bootstrap.trees=TRUE, no.trees=FALSE,keep.overhangs=FALSE, dont.check.duplicates=FALSE, file.alignments=NA, root=NA, align.pairwise.to=NA, out.dir='.')
{	
	stopifnot(is.character(file.bam),is.character(file.ref),is.logical(all.bootstrap.trees))
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
	cmd		<- paste(cmd, prog,' "',basename(file.bam),'" "',basename(file.ref),'" ',sep='')	
	cmd		<- paste(cmd, '--merging-threshold', merge.threshold,'--min-read-count',min.read.count,'--quality-trim-ends', quality.trim.ends, '--min-internal-quality',min.internal.quality)
	if(!is.na(merge.paired.reads))
		cmd	<- paste(cmd,' ',merge.paired.reads,sep='')	
	if(!is.na(dont.check.duplicates))
		cmd	<- paste(cmd,' ',dont.check.duplicates,sep='')
	if(!is.na(align.pairwise.to))		
		cmd	<- paste(cmd,' --pairwise-align-to ',align.pairwise.to,sep='')
	if(nchar(window.automatic))
		cmd	<- paste(cmd,' --auto-window-params ', window.automatic,sep='')
	if(!nchar(window.automatic))
		cmd	<- paste(cmd,' --windows ', paste(as.character(window.coord), collapse=','),sep='')
	if(!is.na(file.alignments))
		cmd	<- paste(cmd,' --alignment-of-other-refs ',file.alignments,sep='')	
	if(!is.na(root))
		cmd	<- paste(cmd,' --ref-for-rooting ',root,sep='')	
	if(!is.na(no.trees))		
		cmd	<- paste(cmd, no.trees)
	if(!is.na(num.bootstraps))
		cmd	<- paste(cmd, ' --num-bootstraps ', num.bootstraps, sep='')
	if(!is.na(keep.overhangs))
		cmd	<- paste(cmd, keep.overhangs)
	cmd		<- paste(cmd, '--x-raxml',prog.raxml,'--x-mafft',prog.mafft,'\n')
	run.id	<- gsub('_bam.txt','',basename(file.bam))	
	if(is.na(no.trees) & (is.na(num.bootstraps) | (!is.na(num.bootstraps) & all.bootstrap.trees)))
		cmd	<- paste(cmd, 'for file in RAxML_bestTree*.tree; do\n\tmv "$file" "${file//RAxML_bestTree\\./',run.id,'_}"\ndone\n',sep='')
	if(is.na(no.trees) & !is.na(num.bootstraps) & !all.bootstrap.trees)
		cmd	<- paste(cmd, 'for file in RAxML_bipartitions.MLtreeWbootstraps*.tree; do\n\tmv "$file" "${file//RAxML_bipartitions.MLtreeWbootstraps/',run.id,'_}"\ndone\n',sep='')	
	#cmd	<- paste(cmd, "for file in AlignedReads*.fasta; do\n\tsed 's/<unknown description>//' \"$file\" > \"$file\".sed\n\tmv \"$file\".sed \"$file\"\ndone\n",sep='')		
	if(!is.na(file.alignments) & !is.na(keep.overhangs))
	{
		cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\tcat "$file" | awk \'{if (substr($0,1,4) == ">REF") censor=1; else if (substr($0,1,1) == ">") censor=0; if (censor==0) print $0}\' > NoRef$file\ndone\n', sep='')
		#cmd	<- paste(cmd, 'for file in NoRefAlignedReads*.fasta; do\n\tcp "$file" "${file//NoRefAlignedReads/NoRef',run.id,'_}"\n\tmv NoRef', run.id,'*fasta "',out.dir,'"\ndone\n',sep='')
		cmd	<- paste(cmd, 'for file in NoRefAlignedReads*.fasta; do\n\t',phsc.cmd.mafft.add(file.alignments,'"$file"','Ref"$file"', options='--keeplength --memsave --parttree --retree 1'),'\ndone\n',sep='')		
		cmd	<- paste(cmd, 'for file in RefNoRefAlignedReads*.fasta; do\n\t','mv "$file" "${file//RefNoRefAlignedReads/',run.id,'_}"\ndone\n',sep='')		
	}
	if(is.na(file.alignments) || is.na(keep.overhangs))
	{
		cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\tmv "$file" "${file//AlignedReads/',run.id,'_}"\ndone\n',sep='')	
	}
	#	run phyloscanner tools and compress output
	cmd		<- paste(cmd, phsc.cmd.process.phyloscanner.output.in.directory(prog, tmpdir, file.bam, root), sep='\n')
	#
	cmd		<- paste(cmd, '\nmv ',run.id,'* "',out.dir,'"\n',sep='')	
	#	zip up everything else
	cmd		<- paste(cmd, 'for file in *; do\n\tzip -ur9XTj ',paste(run.id,'_otherstuff.zip',sep=''),' "$file"\ndone\n',sep='')
	cmd		<- paste(cmd, 'mv ',paste(run.id,'_otherstuff.zip',sep=''),' "',out.dir,'"\n',sep='')
	#	clean up
	cmd		<- paste(cmd,'cd $CWD\nrm -r "',tmpdir,'"\n',sep='')
	cmd
}

#' @export
#' @title Generate bash commands for a multiple phyloscanner runs
#' @description This function generates bash commands for multiple phyloscanner runs, that can be called via 'system' in R, or written to file to run on a UNIX system.  
phsc.cmd.phyloscanner.multi <- function(pty.runs, pty.args) 		
{
	stopifnot(length(pty.args[['win']])==4)
	#
	#	associate BAM and REF files with each scheduled phylotype run
	#	
	#	get available Bam files
	ptyd		<- data.table(FILE=list.files(pty.args[['data.dir']]))
	ptyd[, TYPE:=NA_character_]
	set(ptyd, ptyd[, which(grepl('_ref.fasta$',FILE))], 'TYPE', 'REF')
	set(ptyd, ptyd[, which(grepl('.bam$',FILE))], 'TYPE', 'BAM')
	ptyd		<- subset(ptyd, !is.na(TYPE))
	ptyd[, FILE_ID:= gsub('\\.bam|_ref\\.fasta','',FILE)]
	ptyd		<- dcast.data.table(ptyd, FILE_ID~TYPE, value.var='FILE')
	#	merge
	pty.runs	<- merge(pty.runs, ptyd, by='FILE_ID', all.x=1)
	if(!any(is.na(pty.args[['select']])))
		pty.runs<- subset(pty.runs, PTY_RUN%in%pty.args[['select']])
	tmp			<- subset(pty.runs, is.na(BAM) | is.na(REF))
	if(nrow(tmp))
	{
		print(tmp)
		stop()	#check we have all BAM files		
	}
	pty.runs	<- subset(pty.runs, !is.na(BAM) & !is.na(REF)) 	
	setkey(pty.runs, PTY_RUN)
	#	get alignments file, and root, and taxon against which bams are pairwise aligned to
	alignments.file 			<- system.file(package="phyloscan", "HIV1_compendium_C_B_CPX.fasta")
	if(!is.na(pty.args[['alignments.file']]))
		alignments.file			<- pty.args[['alignments.file']]
	alignments.root				<- PR.ALIGNMENT.ROOT
	if(!is.na(pty.args[['alignments.root']]))		
		alignments.root			<- pty.args[['alignments.root']]	
	tmp							<- rownames(read.dna(alignments.file,format='fa'))
	alignments.root				<- tmp[grepl(alignments.root,tmp)]
	alignments.pairwise.to		<- PR.ALIGNMENT.TO
	if(!is.na(pty.args[['alignments.pairwise.to']]))
		alignments.pairwise.to	<- pty.args[['alignments.pairwise.to']]
	alignments.pairwise.to		<- tmp[grepl(alignments.pairwise.to,tmp)]
	#
	#	write pty.run files and get pty command lines
	#
	pty.c		<- pty.runs[, {
				if(0)
				{
					PTY_RUN		<- z <- 5
					BAM			<- subset(pty.runs, PTY_RUN==z)[, BAM]
					REF			<- subset(pty.runs, PTY_RUN==z)[, REF]					
				}
				file.bam	<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_bam.txt',sep=''))
				file.ref	<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_ref.txt',sep=''))
				cat( paste(file.path(pty.args[['data.dir']],BAM),collapse='\n'), file= file.bam	)
				cat( paste(file.path(pty.args[['data.dir']],REF),collapse='\n'), file= file.ref	)
				windows		<- integer(0)
				if(!nchar(pty.args[['window.automatic']]))
				{
					windows		<- seq(pty.args[['win']][1],pty.args[['win']][2]-pty.args[['win']][3],pty.args[['win']][3])
					windows		<- windows[seq.int(1, length(windows), by=pty.args[['win']][4])]
					windows		<- as.vector(rbind( windows,windows-1+pty.args[['win']][3] ))									
				}
				cmd			<- phsc.cmd.phyloscanner.one(	file.bam, 
										file.ref, 										
										window.coord=windows, 
										window.automatic=pty.args[['window.automatic']], 
										file.alignments=alignments.file, 
										root=alignments.root, 
										align.pairwise.to=alignments.pairwise.to,
										prog=pty.args[['prog']], 
										prog.raxml=pty.args[['raxml']], 
										prog.mafft=pty.args[['mafft']], 
										merge.threshold=pty.args[['merge.threshold']], 
										min.read.count=pty.args[['min.read.count']],
										dont.check.duplicates=pty.args[['dont.check.duplicates']],
										quality.trim.ends=pty.args[['quality.trim.ends']], 
										min.internal.quality=pty.args[['min.internal.quality']], 
										merge.paired.reads=pty.args[['merge.paired.reads']], 
										no.trees=pty.args[['no.trees']], 
										num.bootstraps=pty.args[['num.bootstraps']],
										all.bootstrap.trees=pty.args[['all.bootstrap.trees']],
										keep.overhangs=pty.args[['keep.overhangs']],
										out.dir=pty.args[['out.dir']])
				#cmd			<- paste(cmd, pty.cmd.evaluate.fasta(pty.args[['out.dir']], strip.max.len=pty.args[['strip.max.len']], select=paste('^ptyr',PTY_RUN,'_In',sep=''), min.ureads.individual=pty.args[['min.ureads.individual']]), sep='')
				#cat(cmd)
				list(CMD= cmd)				
			},by='PTY_RUN']
	pty.c
}	

#' @export
#' @import data.table
#' @title Generate bash commands to process phyloscanner output
#' @description This function generates bash commands that combine the various Rscripts in the phyloscanner toolkit   
#' @param tmp.dir Directory with phyloscanner output.
#' @param file.bam File name of the file that contains the list of bam files.
#' @param root.name Name of root taxon. 
#' @return character string of bash commands.
phsc.cmd.process.phyloscanner.output.in.directory<- function(pty.prog, tmp.dir, file.bam, root.name)
{
	#run.id		<- 'ptyr5'; tmp.dir		<- '$CWD/pty_16-09-08-07-32-26'; file.bam	<- '/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptinput/ptyr5_bam.txt'
			
	#	define variables	
	pty.tools.dir		<- file.path(dirname(pty.prog),'tools')
	pty.prog.split		<- paste('Rscript ',file.path(pty.tools.dir,'SplitPatientsToSubtrees.R'),sep='')
	pty.prog.smry		<- paste('Rscript ',file.path(pty.tools.dir,'SummaryStatistics.R'),sep='')
	pty.prog.lkltrm		<- paste('Rscript ',file.path(pty.tools.dir,'LikelyTransmissions.R'),sep='')	
	pty.prog.lkl.smry	<- paste('Rscript ',file.path(pty.tools.dir,'TransmissionSummary.R'),sep='')	
	run.id				<- gsub('_bam.txt','',basename(file.bam))
	run.id_				<- ifelse(grepl('[a-z0-9]',substring(run.id, nchar(run.id))), paste(run.id,'_',sep=''), run.id)	
	#
	#	bash command to plot trees and calculate splits
	#							
	cmd				<- phsc.cmd.SplitPatientsToSubtrees(	pty.prog.split,
															pty.tools.dir,
															file.path(tmp.dir,run.id_),
															outputdir=tmp.dir,																													
															outgroupName=root.name, 
															pdfwidth=30, pdfrelheight=0.15)
	file.patients	<- paste(run.id_,'patients.txt',sep='')	
	cmd				<- paste(cmd,"\nsed 's/.*\\///' \"", file.path(tmp.dir,basename(file.bam)), '" > "',file.path(tmp.dir,file.patients),'"', sep='')
	#
	#	bash command to calculate patient stats
	#							
	tmp				<- phsc.cmd.SummaryStatistics( 	pty.prog.smry, 
													pty.tools.dir, 
													root.name, 
													file.path(tmp.dir, file.patients), 
													file.path(tmp.dir, paste(run.id_,'InWindow_',sep='')), 
													file.path(tmp.dir, paste(run.id_,'InWindow_',sep='')),
													file.path(tmp.dir, paste('Subtrees_r_',run.id_,'InWindow_',sep='')), 
													file.path(tmp.dir, substr(run.id_,1,nchar(run.id_)-1)))
	cmd				<- paste(cmd, tmp, sep='\n')
	#
	#	bash command to get likely transmissions 
	#
	tmp		<- phsc.cmd.LikelyTransmissions(	pty.prog.lkltrm, 
												pty.tools.dir, 
												file.path(tmp.dir,run.id_), 
												file.path(tmp.dir,paste('Subtrees_r_',run.id_,sep='')), 
												file.path(tmp.dir,run.id_),
												root.name=root.name, 
												zeroLengthTipsCount=FALSE, 
												dual.inf.thr=NA,
												romeroSeverson=TRUE)
	cmd		<- paste(cmd, tmp, sep='\n')
	#
	#	add bash command to get likely transmissions summary
	#						
	tmp	<- phsc.cmd.LikelyTransmissionsSummary(	pty.prog.lkl.smry, 
												pty.tools.dir,
												file.path(tmp.dir, file.patients),												
												file.path(tmp.dir, paste(run.id_,'patStatsFull.csv',sep='')),
												file.path(tmp.dir, run.id_), 
												file.path(tmp.dir, paste(run.id_,'trmStats.csv',sep='')),
												min.threshold=1, 
												allow.splits=TRUE)
	cmd	<- paste(cmd, tmp, sep='\n')
	#
	#	add bash command to compress phyloscanner output
	#							
	tmp	<- phsc.cmd.read.processed.phyloscanner.output.in.directory(file.path(tmp.dir, run.id_), 
																	file.path(tmp.dir, run.id_), 
																	read.likelytransmissions=TRUE, 
																	read.trees=TRUE, 
																	read.subtrees=TRUE, 
																	resume=FALSE, 
																	zip=TRUE)
	cmd	<- paste(cmd, tmp, sep='\n')
	cmd
}

#' @export
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
		stat.lkltrm		<- phsc.read.likelytransmissions(prefix.infiles, prefix.run='ptyr', regexpr.lklsu='trmStats.csv$', regexpr.patient='^[0-9]+_[0-9]+_[0-9]+', save.file=save.file, plot.file=plot.file, resume=resume, zip=zip)
		stat.lkltrm		<- NULL		
	}
	#
	#	read trees
	#
	if(read.trees)
	{
		save.file		<- paste(save.file.base,'trees.rda',sep='')
		tmp				<- phsc.read.trees(prefix.infiles, prefix.run='ptyr', regexpr.trees='Subtrees_r_.*\\.rda$', prefix.wfrom='Window_', prefix.wto='Window_[0-9]+_to_', save.file=save.file, resume=resume, zip=zip)
		tmp				<- NULL		
	}
	#
	#	read subtrees files and save (must be run after phsc.read.trees!)
	#
	if(read.subtrees)
	{
		save.file		<- paste(save.file.base,'subtrees_r.rda',sep='')
		stat.subtrees	<- phsc.read.subtrees(prefix.infiles, prefix.run='ptyr', regexpr.subtrees='Subtrees_r_.*\\.rda$', prefix.wfrom='Window_', prefix.wto='Window_[0-9]+_to_', save.file=save.file, resume=resume, zip=zip)
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
#' @param plot.file If not missing, the phylogenies will be printed to file.	
#' @return List of ggtree objects, ready for printing.
phsc.plot.selected.pairs<- function(phs, dfs, id1, id2, plot.file=NA, pdf.h=50, pdf.rw=10)
{
	require(grid)
	phps	<- lapply(seq_len(nrow(dfs)), function(i){
				ph.title	<- NULL
				if('TITLE'%in%colnames(dfs))
					ph.title	<- dfs[i, TITLE]										
				ph			<- phs[[ dfs[i, IDX] ]]
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
						theme(legend.position="bottom") + 
						labs(x='subst/site', title=ph.title)						
				p
			})
	#
	#	single page plot
	#		
	if(!is.na(plot.file))					
	{
		cat('Plotting to file',plot.file,'...\n')
		pdf(file=plot.file, w=pdf.rw*length(phps), h=pdf.h)
		grid.newpage()
		pushViewport(viewport(layout=grid.layout(1, length(phps))))
		for(i in seq_along(phps))
			print(phps[[i]], vp = viewport(layout.pos.row=1, layout.pos.col=i))
		dev.off()
	}
	phps	
}

#' @export
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
	df	<- lapply(seq_len(nrow(dfr)), function(i){				 
				#z	<- '/Users/Oliver/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/pty_Rakai_160825/ptyr1_trmStats.csv'				
				z	<- as.data.table(read.csv(dfr[i,FILE_TRSU], stringsAsFactors=FALSE))
				z[, PTY_RUN:= dfr[i, PTY_RUN]]
				z
			})
	df	<- do.call('rbind',df)
	cat('Found likely transmission entries, n=', nrow(df),'...\n')
	setnames(df, colnames(df), gsub('\\.','_',toupper(colnames(df))))
	setnames(df, 'WINDOWS', 'WIN_OF_TYPE')
	df[, WIN_TOTAL:= as.integer(gsub('/','',regmatches(FRACTION,regexpr('/[0-9]+', FRACTION))))]
	#	reduce to patients that match regexpr.patient
	df	<- subset(df, grepl(regexpr.patient, PAT_1) & grepl(regexpr.patient, PAT_2)) 
	cat('Found likely transmission entries for individuals that meet regexpr, n=', nrow(df),'...\n')	
	df[, ID1:= regmatches(PAT_1,regexpr(regexpr.patient, PAT_1))]
	df[, ID2:= regmatches(PAT_2,regexpr(regexpr.patient, PAT_2))]
	#	reduce from ordered pairs (p1,p2) where order indicates transmission to 
	#		pairs (p1,p2) where the TYPE variable indicates if trm is from p1->p2 or vice versa	
	set(df, df[, which(TYPE=='anc')], 'TYPE', 'anc_12')
	tmp	<- copy(df)
	setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
	set(tmp, tmp[, which(TYPE=='anc_12')], 'TYPE', 'anc_21')
	df	<- rbind(df, tmp)
	#	we have now doubled the data set to (p1,p2) and (p2,p1). 
	#	Can now simply get unique non-ordered pairs by p1<p2, and the new TYPE variable indicates directionality:
	df	<- subset(df, ID1<ID2)
	setkey(df, ID1, ID2)
	cat('Found pairs, n=', nrow(unique(df)),'...\n')
	cat('Found pairs with ancestral relationships, n=', nrow(unique(subset(df, TYPE=='anc_12' | TYPE=='anc_21'))),'...\n')
	#	check total transmissions
	tmp	<- df[, list(OK= sum(WIN_OF_TYPE[TYPE=='anc_12'|TYPE=='anc_21'])==TOTAL_TRANS[1]), by=c('ID1','ID2')]
	stopifnot( nrow(subset(tmp, !OK))==0 )
	set(df, NULL, c('PAT_1','PAT_2','TOTAL_TRANS','FRACTION'), NULL)
	#	determine number of unresolved windows and add as type
	tmp	<- df[, list(WIN_OF_TYPE=WIN_TOTAL[1]-sum(WIN_OF_TYPE), TYPE='disconnected', PTY_RUN=PTY_RUN[1], WIN_TOTAL=WIN_TOTAL[1]), by=c('ID1','ID2')]
	df	<- rbind(tmp, df, use.names=TRUE)
	#	to plot, set overall 'score' to number of windows with anc_12 or anc_21
	tmp	<- df[, list(SCORE=sum(WIN_OF_TYPE[TYPE=='anc_12'|TYPE=='anc_21'])), by=c('ID1','ID2')]
	df	<- merge(df, tmp, by=c('ID1','ID2'))
	#	give every pair an ID
	setkey(df, ID1, ID2)
	tmp	<- unique(df)
	tmp	<- tmp[order(-SCORE),]
	tmp[, PAIR_ID:= seq_len(nrow(tmp))]
	df	<- merge(df, subset(tmp, select=c(ID1,ID2,PAIR_ID)), by=c('ID1','ID2'))
	setkey(df, PAIR_ID)
	#	if save.file, zip summaries, delete individual files and save rda
	if(!is.na(save.file))
	{		
		save(df, file=save.file)				
	}	
	if(zip & !is.na(save.file))
	{		
		tmp	<- copy(dfr)
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE_TRSU, flags = "-ur9XTj")), by='FILE_TRSU'] )
			invisible( file.remove( dfr[, FILE_TRSU] ) )
		}
		tmp	<- data.table(FILE= list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*_LikelyTransmissions.csv$',sep=''), full.names=TRUE))
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'_LikelyTransmissions.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE, flags = "-ur9XTj")), by='FILE'] )
			invisible( file.remove( tmp[, FILE] ) )			
		}
	}	
	if(!is.na(plot.file))
	{		
		cat('\nPlot to file', plot.file,'...\n')
		setkey(df, ID1, ID2)
		tmp	<- unique(df)
		tmp	<- tmp[order(-PAIR_ID),]
		tmp[, LABEL:= factor(PAIR_ID, levels=PAIR_ID, labels=paste('Pair',PAIR_ID, ' (', ID1,'<->', ID2,')',sep=''))]
		tmp	<- merge(subset(tmp, select=c(PAIR_ID, LABEL)), df, by='PAIR_ID')
		set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('anc_12','anc_21','sib','int','disconnected'), labels=c('from 1 to 2','from 2 to 1','1, 2 are siblings','1, 2 are intermingled','1, 2 are disconnected'))])
		ggplot(tmp, aes(x=LABEL, y=WIN_OF_TYPE, fill=TYPE)) +
				geom_bar(stat='identity', position='stack') +
				coord_flip() +
				labs(x='', y='number of read windows', fill='topology of clades from reads\nbetween patient pairs') +
				scale_fill_manual(values=c('from 1 to 2'="#9E0142",'from 2 to 1'="#F46D43",'1, 2 are siblings'="#ABDDA4",'1, 2 are intermingled'="#3288BD",'1, 2 are disconnected'='grey50')) +
				theme_bw() + theme(legend.position='bottom') +
				guides(fill=guide_legend(ncol=2))
		ggsave(plot.file, w=10, h=0.15*nrow(unique(df)), limitsize = FALSE)
	}
	subset(df, select=c(PAIR_ID, ID1, ID2, TYPE, WIN_OF_TYPE, WIN_TOTAL, PTY_RUN, SCORE))
}

#' @export
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
			invisible( tmp[, list(RTN= zip( tmp2, FILE_TR, flags = "-ur9XTj")), by='FILE_TR'] )
			invisible( file.remove( tmp[, FILE_TR] ) )			
		}
		tmp	<- data.table(FILE= list.files(dirname(prefix.infiles), pattern=gsub('\\.rda','\\.csv',regexpr.subtrees), full.names=TRUE))
		tmp	<- subset(tmp, grepl(basename(prefix.infiles), basename(FILE), fixed=1))
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'_csv.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE, flags = "-ur9XTj")), by='FILE'] )			
			invisible( file.remove( tmp[, FILE] ) )			
		}						
	}
	rs.subtrees
}

#' @export
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
			invisible( tmp[, list(RTN= zip( tmp2, FILE, flags = "-ur9XTj")), by='FILE'] )
			invisible( file.remove( tmp[, FILE] ) )			
		}
		tmp	<- data.table(FILE= list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*tree$',sep=''), full.names=TRUE))
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'_newick.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE, flags = "-ur9XTj")), by='FILE'] )			
			invisible( file.remove( tmp[, FILE] ) )			
		}
		tmp	<- data.table(FILE= list.files(dirname(prefix.infiles), pattern=paste(basename(prefix.infiles),'.*fasta$',sep=''), full.names=TRUE))
		if(nrow(tmp))
		{
			tmp2	<- paste(gsub('\\.rda','',save.file),'_fasta.zip',sep='')
			cat('\nZip to file', tmp2,'...\n')
			invisible( tmp[, list(RTN= zip( tmp2, FILE, flags = "-ur9XTj")), by='FILE'] )
			invisible( file.remove( tmp[, FILE] ) )	
		}				
	}
	list(phs=phs, dfr=dfr)
}

