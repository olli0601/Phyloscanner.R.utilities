PR.PACKAGE			<- "phyloscan" 
PR.EVAL.FASTA		<- paste('Rscript', system.file(package=PR.PACKAGE, "phyloscan.evaluate.fasta.Rscript"))
PR.EVAL.EXAML		<- paste('Rscript', system.file(package=PR.PACKAGE, "phyloscan.read.trees.Rscript"))
PR.SCAN.STATISTICS	<- paste('Rscript', system.file(package=PR.PACKAGE, "phyloscan.scan.statistics.Rscript"))
PR.SCAN.SUPERINF	<- paste('Rscript', system.file(package=PR.PACKAGE, "phyloscan.scan.superinfections.Rscript"))
PR.ALIGNMENT.FILE	<- system.file(package=PR.PACKAGE, "HIV1_compendium_C_B_CPX.fasta")
PR.ALIGNMENT.ROOT	<- "R0_CPX_AF460972_read_1_count_0"
EPS					<- 1e-12

#' @export
pty.cmd.evaluate.fasta<- function(indir, outdir=indir, strip.max.len=Inf, select='', min.ureads.individual=NA)
{
	cmd		<- paste('\n',PR.EVAL.FASTA, ' -indir=',indir, ' -strip.max.len=',strip.max.len, sep='')
	if(select!='')
		cmd	<- paste(cmd,' -select=',select,sep='')
	if(!is.na(min.ureads.individual))
		cmd	<- paste(cmd,' -min.ureads.individual=',min.ureads.individual,sep='')	
	cmd	<- paste(cmd,'\n',sep='')
	cmd
}

#' @export
pty.cmd.scan.statistics<- function(indir, outdir=indir, select='')
{
	cmd		<- paste('\n',PR.SCAN.STATISTICS, ' -indir=',indir, sep='')
	if(select!='')
		cmd	<- paste(cmd,' -select=',select,sep='')
	if(!is.na(outdir))
		cmd	<- paste(cmd,' -outdir=',outdir,sep='')	
	cmd	<- paste(cmd,'\n',sep='')
	cmd
}


#' @export
pty.cmd.scan.superinfections<- function(indir, outdir=indir, select='', references.pattern='REF', run.pattern='ptyr', plot.max.clade=0)
{
	cmd		<- paste('\n',PR.SCAN.SUPERINF, ' -indir=',indir, ' -plot.max.clade=', plot.max.clade, sep='')
	if(select!='')
		cmd	<- paste(cmd,' -select=',select,sep='')
	if(references.pattern!='')
		cmd	<- paste(cmd,' -references.pattern=',references.pattern,sep='')
	if(run.pattern!='')
		cmd	<- paste(cmd,' -run.pattern=',run.pattern,sep='')
	cmd	<- paste(cmd,'\n',sep='')
	cmd
}

#' @export
pty.cmd.evaluate.examl<- function(indir, outdir=indir, select='', outgroup=NA, references.pattern='REF', run.pattern='ptyr', rm.newick=0, rm.fasta=0, infile=NA)
{
	cmd		<- paste('\n',PR.EVAL.EXAML, ' -indir=',indir, ' -rm.newick=', rm.newick, ' -rm.fasta=', rm.fasta, sep='')
	if(!is.na(infile))
		cmd	<- paste(cmd,' -infile=',infile,sep='')		
	if(select!='')
		cmd	<- paste(cmd,' -select=',select,sep='')
	if(references.pattern!='')
		cmd	<- paste(cmd,' -references.pattern=',references.pattern,sep='')
	if(run.pattern!='')
		cmd	<- paste(cmd,' -run.pattern=',run.pattern,sep='')
	if(!is.na(outgroup))
		cmd	<- paste(cmd,' -outgroup=',outgroup,sep='')	
	cmd	<- paste(cmd,'\n',sep='')
	cmd
}

#' @export
pty.cmd.mafft.add<- function(infile, reffile, outfile, options='')
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
pty.cmd<- function(file.bam, file.ref, window.coord=integer(0), window.automatic='', prog=PROG.PTY, prog.raxml='raxmlHPC-AVX', prog.mafft='mafft', merge.threshold=1, min.read.count=2, quality.trim.ends=30, min.internal.quality=2, merge.paired.reads='-P',no.trees='',keep.overhangs='', ref.file=PR.ALIGNMENT.FILE, ref.root=PR.ALIGNMENT.ROOT, out.dir='.')
{	
	stopifnot(is.character(file.bam),is.character(file.ref))
	if(!nchar(window.automatic))	stopifnot( is.numeric(window.coord), !length(window.coord)%%2)
	if(nchar(window.automatic))		stopifnot( !length(window.coord) )
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
	cmd		<- paste(cmd, prog,' ',merge.threshold,' ',min.read.count,' "',basename(file.bam),'" "',basename(file.ref),'" ',paste(as.character(window.coord), collapse=' '),sep='')	
	cmd		<- paste(cmd, '-Q1', quality.trim.ends, '-Q2',min.internal.quality, merge.paired.reads)	
	if(nchar(window.automatic))
		cmd	<- paste(cmd,' --auto-window-params ', window.automatic,sep='')
	if(nchar(ref.file) & keep.overhangs!='--keep-overhangs')
		cmd	<- paste(cmd,' --alignment-of-other-refs ',ref.file,sep='')	
	if(nchar(ref.root) & keep.overhangs!='--keep-overhangs')
		cmd	<- paste(cmd,' --ref-for-rooting ',ref.root,sep='')
	if(nchar(no.trees))		
		cmd	<- paste(cmd, no.trees)	
	if(nchar(keep.overhangs))
		cmd	<- paste(cmd, keep.overhangs)
	cmd		<- paste(cmd, '--x-raxml',prog.raxml,'--x-mafft',prog.mafft,'\n')
	tmp		<- gsub('_bam.txt','',basename(file.bam))	
	cmd		<- paste(cmd, 'for file in RAxML_bestTree\\.*.tree; do\n\tmv "$file" "${file//RAxML_bestTree\\./',tmp,'_}"\ndone\n',sep='')
	cmd		<- paste(cmd, "for file in AlignedReads*.fasta; do\n\tsed 's/<unknown description>//' \"$file\" > \"$file\".sed\n\tmv \"$file\".sed \"$file\"\ndone\n",sep='')		
	if(nchar(ref.file) & keep.overhangs=='--keep-overhangs')
	{
		cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\t',pty.cmd.mafft.add(ref.file,'"$file"','Ref"$file"', options='--keeplength'),'\ndone\n',sep='')		
		cmd	<- paste(cmd, 'for file in RefAlignedReads*.fasta; do\n\t','mv "$file" "${file//RefAlignedReads/',tmp,'_}"\ndone\n',sep='')		
	}
	if(!nchar(ref.file) || keep.overhangs!='--keep-overhangs')
	{
		cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\tmv "$file" "${file//AlignedReads/',tmp,'_}"\ndone\n',sep='')	
	}	
	cmd		<- paste(cmd, 'mv ',tmp,'*tree "',out.dir,'"\n',sep='')	
	cmd		<- paste(cmd, 'mv ',tmp,'*fasta "',out.dir,'"\n',sep='')
	#	clean up
	cmd		<- paste(cmd,'cd $CWD\nrm -r "',tmpdir,'"\n',sep='')
	cmd
}

#' @export
pty.cmdwrap.fasta <- function(pty.runs, pty.args) 		
{
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
	tmp			<- subset(pty.runs, is.na(BAM) | is.na(REF))
	if(nrow(tmp))
	{
		print(tmp)
		stop()	#check we have all BAM files		
	}
	pty.runs	<- subset(pty.runs, !is.na(BAM) & !is.na(REF)) 	
	#	determine length of each ref file
	tmp			<- unique(subset(pty.runs, select=REF))
	tmp			<- tmp[, list(REF_LEN=ncol(read.dna(file.path(pty.args[['data.dir']],REF), format='fasta'))), by='REF']
	pty.runs	<- merge(pty.runs, tmp, by='REF')
	setkey(pty.runs, PTY_RUN)
	#
	#	write pty.run files and get pty command lines
	#
	pty.win		<- pty.args[['win']]
	pty.cmd		<- pty.runs[, {
				if(0)
				{
					PTY_RUN		<- 2
					BAM			<- subset(pty.runs, PTY_RUN==2)[, BAM]
					REF			<- subset(pty.runs, PTY_RUN==2)[, REF]
					REF_LEN		<- subset(pty.runs, PTY_RUN==2)[, REF_LEN]
				}
				file.bam	<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_bam.txt',sep=''))
				file.ref	<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_ref.txt',sep=''))
				cat( paste(file.path(pty.args[['data.dir']],BAM),collapse='\n'), file= file.bam	)
				cat( paste(file.path(pty.args[['data.dir']],REF),collapse='\n'), file= file.ref	)
				windows		<- integer(0)
				if(!nchar(pty.args[['window.automatic']]))
				{
					windows		<- seq(1,by=pty.win,len=ceiling(max(REF_LEN)/pty.win))
					windows		<- as.vector(rbind( windows,windows-1+pty.win ))									
				}
				cmd			<- pty.cmd(	file.bam, file.ref, 
										window.coord=windows, window.automatic=pty.args[['window.automatic']],
										prog=pty.args[['prog']], prog.raxml=pty.args[['raxml']], prog.mafft=pty.args[['mafft']], 
										merge.threshold=pty.args[['merge.threshold']], min.read.count=pty.args[['min.read.count']], quality.trim.ends=pty.args[['quality.trim.ends']], min.internal.quality=pty.args[['min.internal.quality']], 
										merge.paired.reads=pty.args[['merge.paired.reads']], no.trees=pty.args[['no.trees']], keep.overhangs=pty.args[['keep.overhangs']],
										out.dir=pty.args[['out.dir']])
				cmd			<- paste(cmd, pty.cmd.evaluate.fasta(pty.args[['out.dir']], strip.max.len=pty.args[['strip.max.len']], select=paste('^ptyr',PTY_RUN,'_In',sep=''), min.ureads.individual=pty.args[['min.ureads.individual']]), sep='')
				#cat(cmd)
				list(CMD= cmd)				
			},by='PTY_RUN']
	pty.cmd
}	

#' @export
project.dualinfecions.phylotypes.mltrees.160115<- function() 
{
	require(ggtree)
	require(phytools)
	
	pty.infile		<- file.path(HOME,"data", "PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_ptyrunsinput.rda")
	indir.tr		<- file.path(HOME,"phylotypes_160120")
	pty.evaluate.tree(pty.infile, indir.tr)
	#	can now delete all newick files	
	
	pty.infile		<- file.path(HOME,"data", "PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_ptyrunstrees.rda")
	#
	#	get statistics
	#
	#	node heights
	tmp					<- ptyfiles[,{										
										ph		<- pty.ph[[FILE]]
										tmp		<- node.depth.edgelength(ph)[1:Ntip(ph)]
										list(BAM=ph$tip.label, HEIGHT=tmp)
									}, by='FILE']
	pty.stat			<- merge(ptyfiles, tmp, by='FILE')
	#	extract individual from read name
	pty.stat[, IND:= gsub('_read.*','',BAM)]
	#	extract number of identical reads from read name 
	pty.stat[, READ_N:= as.integer(gsub('count_','',regmatches(BAM, regexpr('count_[0-9]+',BAM))))]
	#	get number of reads per individual
	pty.stat			<- merge(pty.stat,pty.stat[, list(IND_N=sum(READ_N)),by=c('FILE','IND')],by=c('FILE','IND'))
	#	check if reads from same individual are monophyletic
	pty.stat			<- merge(pty.stat, pty.stat[,	list(MONOPH=is.monophyletic(pty.ph[[FILE]], BAM)), by=c('FILE','IND')], by=c('FILE','IND'))
	set(pty.stat, NULL, 'MONOPH', pty.stat[, factor(MONOPH, levels=c(FALSE,TRUE),labels=c('N','Y'))])
	
	
	
	ptyfiles				<- merge(ptyfiles, pty.stat[, list(HMX= max(HEIGHT)), by='FILE'], by='FILE')
	
	
	# For each patient: record whether his/her tips are monophyletic, find the
	# pairwise patristic distances between the tips - the 'cophenetic distances' -
	# and characterise those distances. 
	dummy.p.value<-0
	
	for (i in 1:num.ids) {
		id <- ids[i]
		num.leaves <- length(patient.tips[[id]])
		num.reads<-0
		for (tip in patient.tips[[id]]) num.reads <- num.reads + as.numeric(unlist(strsplit(tip,"count_"))[2])
		if (num.leaves>0) {
			monophyletic <- as.numeric(is.monophyletic(tree, patient.tips[[id]]))
			if (num.leaves > 1) {
				subtree <- drop.tip(phy=tree,
						tip=tree$tip.label[!(tree$tip.label %in% patient.tips[[id]])])
				subtree.dist.matrix <- cophenetic(subtree)
				subtree.dist <- subtree.dist.matrix[lower.tri(subtree.dist.matrix)]
				mean.size <- mean(subtree.dist)
				variance <- var(subtree.dist)
				coeff.of.var.size <- ifelse(num.leaves > 2, sqrt(variance)/mean.size, 0)
				## Distances are not weighted for the number of reads for each tip. 
				## Corrections are included because mean and variance of distance matrices
				## include zeros on diagonal, but shouldn't.
				#subtree.dist.cf <- as.vector(subtree.dist.matrix)
				#mean.size.cf <- mean(subtree.dist.cf)/(1-1/num.leaves)
				#variance.cf <- var(subtree.dist.cf)*(1+1/num.leaves)-mean.size.cf^2/num.leaves
				#cat("CF mean & var: ", mean.size.cf, variance.cf, '\n')
				#cat("CW mean & var: ", mean.size.cw, variance.cw, '\n')
				#cat('\n')
			} else {
				mean.size <- NA
				coeff.of.var.size <- NA
			}
			root.to.tip<-0
			for (i in 1:length(subtree$tip.label)) root.to.tip<-root.to.tip+nodeheight(subtree,i)
			root.to.tip<-root.to.tip/length(subtree$tip.label)
		} else {
			monophyletic<-NA
			mean.size<-NA
			coeff.of.var.size<-NA
			root.to.tip<-NA
		}
		pat.stats <- rbind(pat.stats, c(id, window, num.leaves, num.reads, monophyletic,
						mean.size, coeff.of.var.size,root.to.tip))
	}	
}

#' @export
pty.evaluate.tree.root.nofill<- function(ptyfiles, pty.ph, pty.runs)
{
	pty.root	<- lapply(ptyfiles[, unique(PTY_RUN)], function(ptyr){	
				#print(ptyr)	
				#ptyr<- 15
				tmp					<- subset(ptyfiles, PTY_RUN==ptyr)[,FILE]
				phs					<- lapply(tmp, function(x) pty.ph[[x]])	
				names(phs)			<- tmp
				#	number of read counts per individual				
				phpd				<- do.call('rbind',lapply(seq_along(phs),function(i){
									#i			<- 1
									#print(i)		
									ph			<- phs[[i]]
									phb			<- data.table(IDX=seq_along(ph$tip.label), BAM=ph$tip.label, FILE_ID= gsub('_read.*','',ph$tip.label))
									phb[, COUNT:= as.numeric(gsub('count_','',regmatches(BAM, regexpr('count_[0-9]+',BAM))))]
									phb			<- merge(phb, subset(pty.runs[pty.runs$PTY_RUN==ptyr,],select=FILE_ID), by='FILE_ID', all=1)
									set(phb, phb[, which(is.na(COUNT))], 'COUNT', 0)
									phb			<- phb[, list(COUNT=sum(COUNT)), by='FILE_ID']
									phb[, FILE:= names(phs)[i]]
									phb					
								}))
				phpd				<- merge(phpd, subset(ptyfiles, PTY_RUN==ptyr, select=c(FILE, W_FROM)), by='FILE')
				setkey(phpd, W_FROM)
				outgroup.order		<- phpd[, list(COUNT=mean(COUNT)), by='FILE_ID']				
				outgroup.order		<- merge(outgroup.order, phpd[, list(PRESENT=length(which(COUNT>0))), by='FILE_ID'], by='FILE_ID')
				setkey(outgroup.order, COUNT)
				outgroup.order		<- subset(outgroup.order, PRESENT>0)
				#	ususally not all present throughout; pick in order
				ans					<- do.call('rbind',lapply(seq_along(phs),function(i){
							#i			<- 147
							#print(i)
							ph				<- phs[[i]]							
							phb				<- data.table( IDX=seq_along(ph$tip.label), BAM=ph$tip.label, FILE_ID= gsub('_read.*','',ph$tip.label))				
							phb				<- merge(phb, outgroup.order, by='FILE_ID')
							outgroup.ind	<- phb[which.min(COUNT), FILE_ID]							
							phgd			<- cophenetic.phylo(ph)
							tmp				<- subset(phb, FILE_ID==outgroup.ind)[, BAM]
							if(length(tmp)<nrow(phgd))
								phgd		<- phgd[ tmp, setdiff(rownames(phgd), tmp), drop=FALSE]
							outgroup.seq	<- rownames(phgd)[ which(min(phgd)==phgd, arr.ind=TRUE)[1,1] ]
							data.table(FILE=names(phs)[i], ROOT=outgroup.seq)					
						}))
				ans[, PTY_RUN:= ptyr]
				ans				
			})
	pty.root	<- do.call('rbind',pty.root)
	pty.root
}

#' @export
pty.evaluate.tree.root.withfill<- function(ptyfiles, pty.ph, pty.runs)
{
	pty.root	<- lapply(ptyfiles[, unique(PTY_RUN)], function(ptyr){	
				#print(ptyr)	
				#ptyr<- 15
				tmp			<- subset(ptyfiles, PTY_RUN==ptyr)[,FILE]
				phs			<- lapply(tmp, function(x) pty.ph[[x]])	
				names(phs)	<- tmp
				#get patristic distances between target and fill
				#if no target, select target_ind with lowest taxon index (just makes an unambiguous selection when windows are considered one by one)
				phpd		<- do.call('rbind',lapply(seq_along(phs),function(i){
									#i			<- 1
									print(i)
									ph			<- phs[[i]]							
									phb			<- data.table(IDX=seq_along(ph$tip.label), BAM=ph$tip.label, FILE_ID= gsub('_read.*','',ph$tip.label))				
									phb			<- merge(phb, pty.runs[pty.runs$PTY_RUN==ptyr,], by='FILE_ID', all=1)
									phgd		<- cophenetic.phylo(ph)
									tmp			<- subset(phb, !FILL & !is.na(BAM))[, BAM]
									if(length(tmp)==0)
									{
										tmp		<- subset(phb, !is.na(BAM))[, FILE_ID[which.min(TX_IDX)]]
										tmp		<- subset(phb, !is.na(BAM) & FILE_ID==tmp)[, BAM]
									}								
									tmp			<- phgd[ tmp, setdiff(rownames(phgd), tmp), drop=FALSE]
									stopifnot(ncol(tmp)>0)								
									ans			<- as.data.table(melt(tmp, value.name='PATR'))
									setnames(ans, c('Var1','Var2'), c('TARGET','FILL'))
									ans[, FILE:= names(phs)[i]]	
									ans[, IDX:=i]							
									ans							
								}))
				#print(phpd[, unique(IDX)])				
				phpd[, FILL_IND:= gsub('_read.*','',FILL)]				
				#print(phpd)		
				#	try consider as root only individual present across all BAM files
				tmp		<- phpd[, list(FILL_IND_N=length(FILL)), by=c('FILL_IND','FILE')]
				tmp		<- dcast.data.table(tmp, FILE~FILL_IND, value.var='FILL_IND_N')
				tmp		<- apply(tmp[,-1, with=FALSE], 2, function(x) !any(is.na(x))	)
				tmp		<- data.table(FILL_IND=names(tmp)[tmp])
				if(nrow(tmp)==0)	# 	may be empty
				{
					print(PTY_RUN)
					stop()
				}				
				tmp		<- merge(phpd, tmp, by='FILL_IND')
				#	select individual with average largest distance
				tmp		<- tmp[, list(PATR= median(PATR)), by='FILL_IND']	
				root	<- tmp[which.max(PATR),][,FILL_IND]
				ans		<- subset(phpd, FILL_IND==root)[, list(ROOT=FILL[which.max(PATR)]), by=c('IDX','FILE')]
				tmp		<- ans[, list(CHECK= ROOT%in%phs[[IDX]]$tip.label) , by='IDX']
				stopifnot( tmp[, all(CHECK)] )		
				ans[, IDX:=NULL]
				ans[, PTY_RUN:= ptyr]
				ans
			})
	pty.root	<- do.call('rbind',pty.root)
	pty.root
}

#' @export
pty.stat.different<- function(ph)
{
	phb			<- data.table(IDX=seq_along(ph$tip.label), BAM=ph$tip.label, FILE_ID= gsub('_read.*','',ph$tip.label))				
	phb[, COUNT:= as.numeric(gsub('count_','',regmatches(BAM, regexpr('count_[0-9]+',BAM))))]				
	phm			<- phb[, list(MRCA=as.numeric(getMRCA(ph, IDX))), by='FILE_ID']
	#	find longest branch and diversity in clade below
	phm			<- phm[, {
				#FILE_ID<- '13554_1_18'
				#MRCA<- 2331
				z	<- extract.clade(ph, MRCA, root.edge=1)
				#plot(z, show.tip.label=0)
				#add.scale.bar()
				# tips in same clade have consecutive tip indices in ape format
				df	<- data.table(IDX=seq_along(z$tip.label), BAM=z$tip.label, FILE_ID= gsub('_read.*','',z$tip.label))
				df	<- merge(df, data.table(IDX= df$IDX[which(df$FILE_ID!=FILE_ID)]), by='IDX')
				tmp	<- c(nrow(df), df[, length(unique(FILE_ID))])
				if(tmp[1])
				{
					
					df[, GROUP:= cumsum(c(1,as.numeric(diff(IDX)>1)))]
					df[, GROUP:= df[, as.numeric(factor(paste(FILE_ID,'-',GROUP,sep='')))]]
					df[, GROUP:= df[, cumsum(c(1,as.numeric(diff(GROUP)>1)))]]
					tmp[1]	<- df[, length(unique(GROUP))]					
				}
				list(DIFF=tmp[1], DIFF_IND=tmp[2] )								
			}, by='FILE_ID']
	phm
}

#' @export
pty.stat.maxlocalsep<- function(ph, threshold.min.stem=0.01)
{
	phb			<- data.table(IDX=seq_along(ph$tip.label), BAM=ph$tip.label, FILE_ID= gsub('_read.*','',ph$tip.label))				
	phb[, COUNT:= as.numeric(gsub('count_','',regmatches(BAM, regexpr('count_[0-9]+',BAM))))]				
	phm			<- phb[, list(MRCA=as.numeric(getMRCA(ph, IDX))), by='FILE_ID']
	#	find longest branch and diversity in clade below
	phm			<- phm[, {
				#FILE_ID<- 'R3_res567_S22_L001'; MRCA<- 1213
				#print(MRCA)
				z		<- extract.clade(ph, MRCA, root.edge=1)
				z		<- drop.tip(z, z$tip.label[ !grepl(FILE_ID, z$tip.label) ], root.edge=1)
				#plot(z, cex=0.2)
				#add.scale.bar()
				# find largest internal branch							
				df		<- data.table(E_IDX=which(z$edge[,2]>Ntip(z)))
				df[, LEN:= z$edge.length[E_IDX]]
				df[, TO:= z$edge[E_IDX,2]]
				df		<- subset( df[order(-LEN),], LEN>threshold.min.stem )
				ans		<- data.table(E_IDX=NA_integer_, CL_MX_LOCAL_SEP=NA_real_, CL_AVG_HEIGHT_UNIQUE=NA_real_, CL_AVG_HEIGHT_ALL=NA_real_, CL_TAXA_N=NA_integer_, CL_COUNT_N=NA_real_)
				if(nrow(df))
					ans	<- df[, {
							tmp	<- extract.clade(z, TO, root.edge=0)
							tmp2<- node.depth.edgelength(tmp)[seq_len(Ntip(tmp))]	
							tmp3<- merge(phb,data.table(BAM=tmp$tip.label),by='BAM')[, COUNT]
							list(	CL_MX_LOCAL_SEP=LEN, 
									CL_AVG_HEIGHT_UNIQUE=mean(tmp2),
									CL_AVG_HEIGHT_ALL=mean(rep(tmp2,tmp3)),
									CL_TAXA_N=Ntip(tmp),
									CL_COUNT_N=sum(tmp3))
						}, by='E_IDX']
				ans[, MRCA:=MRCA]				 
				ans[, TAXA_N:=Ntip(z)]
				ans[, COUNT_N:=merge(phb,data.table(BAM=z$tip.label),by='BAM')[, sum(COUNT)]]
				ans
			}, by='FILE_ID']
	subset(phm, !is.na(E_IDX))
}

#' @export
pty.stat.withinhost.diversity<- function(ph, coi.div.probs=c(0.02,0.25,0.5,0.75,0.98))
{
	phb			<- data.table( BAM=ph$tip.label, FILE_ID= gsub('_read.*','',ph$tip.label))				
	phb[, COUNT:= as.numeric(gsub('count_','',regmatches(BAM, regexpr('count_[0-9]+',BAM))))]	
	phgd		<- cophenetic.phylo(ph)	
	phgd		<- as.data.table(melt(phgd, value.name='PATR'))
	setnames(phgd, c('Var1','Var2'), c('BAM','BAM2'))
	phgd		<- merge(phb, phgd, by='BAM')
	setnames(phgd, c('BAM','FILE_ID','COUNT','BAM2'), c('BAM1','FILE_ID1','COUNT1','BAM'))
	phgd		<- merge(phb, phgd, by='BAM')
	setnames(phgd, c('BAM','FILE_ID','COUNT'), c('BAM2','FILE_ID','COUNT2'))
	phgd		<- subset(phgd, FILE_ID1==FILE_ID)
	phgd[, FILE_ID1:= NULL]
	tmp			<- phgd[, list(READS='UNIQUE', P=coi.div.probs, QU=quantile(PATR,p=coi.div.probs)), by='FILE_ID']
	ans			<- phgd[, list(READS='ALL', P=coi.div.probs, QU=quantile(rep(PATR,COUNT1*COUNT2),p=coi.div.probs)), by='FILE_ID']
	ans			<- rbind(ans,tmp)
	set(ans, NULL, 'P', ans[, paste('WHD_p',P*100,sep='')])
	ans			<- dcast.data.table( ans, READS+FILE_ID~P, value.var='QU')
	tmp			<- phb[, list(READS='ALL', N=sum(COUNT)), by='FILE_ID']
	tmp			<- rbind(tmp, phb[, list(READS='UNIQUE', N=length(BAM)), by='FILE_ID'])
	ans			<- merge(ans,tmp,by=c('READS','FILE_ID'))	
}

#' @export
pty.stat.collect	<- function(indir, outdir=indir, outfile=file.path(outdir,'ptyr_examl_stat.rda'))
{
	infiles		<- data.table(FILE=list.files(indir, pattern='examl.rda$'))
	infiles[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]
	tmp			<- data.table(FILE_STAT=list.files(indir, pattern='examl_stat.rda$'))
	tmp[, FILE:= gsub('_stat','',FILE_STAT)]
	infiles		<- merge(infiles, tmp, by='FILE',all.x=1)	
	setkey(infiles, PTY_RUN)
	if( infiles[, any(is.na(FILE_STAT))] )
		cat('\nwarning stat files not available for=',paste( subset(infiles, is.na(FILE_STAT))[, FILE], collapse=', '))	
	infiles		<- subset(infiles, !is.na(FILE_STAT))	
	pty.stat	<- do.call('rbind',lapply(seq_len(nrow(infiles)), function(i){
						load( gsub('\\.rda','_stat\\.rda', file.path(indir, infiles[i,FILE])) )				
						pty.stat
					}))	
	save(pty.stat, file=outfile)
}


#' @export
pty.stat.monophyletic.clades<- function(ph)
{
	phb			<- data.table(IDX=seq_along(ph$tip.label), BAM=ph$tip.label, IND= gsub('_read.*','',ph$tip.label), REF=grepl(references.pattern,ph$tip.label))
	set(phb, phb[, which(REF)],'IND','REFERENCE')
	#	for each patient define mrca (MRCA)
	#	for each patient, 
	#		determine MRCA (mrca)
	#		calculate number of other individuals in clade below MRCA (diff)	
	z			<- phb[, {
				#print(GROUP)
				mrca	<- IDX		
				diff	<- 0L
				if(length(IDX)>1)
				{
					mrca	<- as.integer(getMRCA(ph, IDX))
					tmp		<- extract.clade(ph, mrca, root.edge=1)
					#print(tmp)
					diff	<- length(setdiff(gsub('_read.*','',tmp$tip.label), IND))								
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
	if(nrow(tmp))
	{
		#	find all tips that belong to another patient than the current individual
		tmp			<- tmp[, {
					if(MRCA<=Ntip(ph))
						tmp	<- ph$tip.label[MRCA]
					if(MRCA>Ntip(ph))
						tmp	<- extract.clade(ph,MRCA)$tip.label
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
					tmp			<- lapply(seq_along(anc.group), function(i) sapply(anc.group[[i]], function(j)	any(grepl(MRCA_IND, extract.clade(ph,j)$tip.label))		)	)				
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
		setkey(zz, MRCA_IND, CHANGE_NODE)
		zz	<- unique(zz)
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
									ans	<- all(grepl(MRCA_IND, extract.clade(ph,x)$tip.label))
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
	#	double check CLADEs (the MONOPH_PA condition is key)
	tmp	<- z[,{
				if(CLADE<=Ntip(ph))
					tmp	<- grepl(IND, ph$tip.label[CLADE])
				if(CLADE>Ntip(ph))
					tmp	<- all(grepl(IND, extract.clade(ph,CLADE)$tip.label))
				list(MONOPH=tmp, MONOPH_PA=all(grepl(IND, extract.clade(ph,Ancestors(ph, CLADE, type="parent"))$tip.label)))
			}, by=c('IND','CLADE')]
	stopifnot(tmp[, all(MONOPH)], tmp[, !any(MONOPH_PA)])
	z	
}


pty.stat.monophyletic.clades.using.mxpars<- function(ph)
{
	phb			<- data.table(IDX=seq_along(ph$tip.label), BAM=ph$tip.label, IND= gsub('_read.*','',ph$tip.label), REF=grepl(references.pattern,ph$tip.label))
	set(phb, phb[, which(REF)],'IND','REFERENCE')
	#phb[, COUNT:= as.numeric(gsub('count_','',regmatches(BAM, regexpr('count_[0-9]+',BAM))))]
	#	consecutive tips of the same individual could define separate clades
	phb[, GROUP:=cumsum(c(1,as.numeric(diff(as.numeric(factor(IND)))!=0)))]
	#	for each potential clade, 
	#		determine MRCA (mrca)
	#		calculate number of other individuals in clade below MRCA (diff)
	#		get node index of ancestor of MRCA (mrcaa)
	z			<- phb[, {
				#print(GROUP)
				mrca	<- IDX		
				diff	<- 0L
				mrcaa	<- NA_real_
				if(length(IDX)>1)
				{
					mrca	<- as.integer(getMRCA(ph, IDX))
					tmp		<- extract.clade(ph, mrca, root.edge=1)
					#print(tmp)
					diff	<- length(setdiff(gsub('_read.*','',tmp$tip.label), IND))								
				}												
				#ancestor of MRCA
				mrcaa	<- ph$edge[which(ph$edge[,2]==mrca),1]																
				list(MRCA=mrca, DIFF_IND=diff, MRCA_ANC=mrcaa)							
			}, by=c('GROUP','IND')]
	z			<- subset(z, IND!='REFERENCE')
	#	get max parsimony assignment for edge ending in ancestor of MRCA
	set(z, NULL, 'MRCA_ANC_EDGE', z[, attr(ph, "INDIVIDUAL")[MRCA_ANC]])
	set(z, NULL, 'MONOPH', 0L)
	#	if DIFF_IND==0 and MRCA_ANC_EDGE!=IND, then we have a unique clade
	set(z, z[, which(DIFF_IND==0L & MRCA_ANC_EDGE!=IND)], 'MONOPH', 1L)
	#	no double counting for every individual: ignore nested clades
	tmp			<- z[, {
				tmp<- Descendants(ph, min(MRCA), type='all')
				list(MRCA=MRCA, SUBCLADE= MRCA%in%tmp)
			}, by='IND']
	z			<- merge(z, subset(tmp, !SUBCLADE), by=c('IND','MRCA'))
	z[, SUBCLADE:=NULL]	
	#	for each potential clade with >=1 other individuals within, there could be several clades
	#		trace back all ancestors between tips of other individuals and MRCA
	#		for each such tip, construct the unique path to MRCA that is not on a previous path
	#		find change points in max pars edge assignments along the unique paths
	#		the idea is that this finds the breaking edges that result in monophylies
	tmp			<- subset(z, DIFF_IND>0, c(GROUP, IND, MRCA))
	if(nrow(tmp))
	{
		tmp			<- tmp[, {
					if(MRCA<=Ntip(ph))
						tmp	<- ph$tip.label[MRCA]
					if(MRCA>Ntip(ph))
						tmp	<- extract.clade(ph,MRCA)$tip.label
					list(MRCA=MRCA, MRCA_IND=IND, BAM=tmp)	
				}, by='GROUP']
		zz			<- merge(tmp, subset(phb, select=c(BAM, IND, IDX)), by='BAM')
		zz			<- subset(zz, MRCA_IND!=IND)
		zz			<- zz[, {
					#IDX<- c(980,910); MRCA<- 1580; IND<- 'R1_RES669_S20_L001'; MRCA_IND<- 'R1_RES827_S23_L001'
					#IDX<- c(917); MRCA<- 2690; IND<- 'R10_105205_S83_L001'; MRCA_IND<- 'R1_RES669_S20_L001'
					#	determine paths to MRCA from each tip
					anc.group		<- Ancestors(ph, IDX)	
					if(!is.list(anc.group))
						anc.group	<- list(anc.group)
					anc.group	<- lapply(seq_along(anc.group), function(i)  anc.group[[i]][ seq_len(which(anc.group[[i]]==MRCA[1])-1)] )							
					#	determine unique paths until we hit a path that is already visited
					anc.join	<- lapply(seq_along(anc.group), function(i){	unique(unlist(lapply(seq_len(i), function(j) anc.group[[j]])))	})
					anc.join	<- c(NA,anc.join)
					anc.group	<- lapply(seq_along(anc.group), function(i)	setdiff(anc.group[[i]],anc.join[[i]])	)
					#	determine max parsimony assignment along unique branches
					#	to find patient-clades within the subtree containing all patient taxa, it does not matter if branches change *between* other individuals
					mxp.assign	<- as.character(attr(ph, "INDIVIDUAL"))
					mxp.assign[ mxp.assign!=MRCA_IND[1] ]	<- 'OTHER'
					tmp			<- ifelse(IND==MRCA_IND[1], IND, 'OTHER')
					tmp			<- lapply(seq_along(anc.group), function(i)	c(tmp, mxp.assign[ anc.group[[i]] ]) )
					#	count changes in max parsimony assignment along each branch
					mxp.changes	<- lapply(seq_along(tmp), function(i)	sum( as.numeric(diff(as.numeric(factor(tmp[[i]])))!=0) )	)							
					#	determine lowest node at which change occurs
					tmp			<- lapply(seq_along(tmp), function(i){
								tmp2<- tail(which(diff(as.numeric(factor(tmp[[i]])))!=0),1) #index of ancestor *below* which change occurs
								if(length(tmp2))
									tmp2<- as.integer(anc.group[[i]][tmp2])										
								if(!length(tmp2))
									tmp2<- NA_integer_
								tmp2
							})				
					list(IDX=IDX, MRCA_IND=MRCA_IND, CHANGES_ON_UNIQUE_PATH=unlist(mxp.changes), CHANGE_NODE=unlist(tmp))
				},by=c('GROUP','IND','MRCA')]
		zz			<- subset(zz, CHANGES_ON_UNIQUE_PATH>0)
		#
		#	each ancestor before a change node could have as one of its children a monophyletic clade
		#	
		setkey(zz, GROUP, CHANGE_NODE)
		zz	<- unique(zz)
		zz	<- zz[, {
					#CHANGE_NODE<- 2053; MRCA<- 1580; MRCA_IND<- 'R1_RES827_S23_L001'
					#cat(CHANGE_NODE)
					#	path from change node to just before mrca; we will be checking the sibling of these so this is why the mrca is excluded
					path		<- c(CHANGE_NODE, setdiff( Ancestors(ph,CHANGE_NODE), c(MRCA,Ancestors(ph,MRCA)) ))
					#	the first potentially monophyletic clade must be calculated in a different way							
					tmp			<- ph$edge[ph$edge[,1]==CHANGE_NODE,2]
					tmp2		<- which(as.character(attr(ph, "INDIVIDUAL"))[tmp]==MRCA_IND)
					#	collect the potentially monophyletic clades along the path
					pot.clade	<- c(tmp[tmp2], unlist(Siblings(ph, path)))
					#	check if potential clades are monophyletic
					tmp			<- sapply(pot.clade, function(x){
								if(x<=Ntip(ph))
									ans	<- grepl(MRCA_IND, ph$tip.label[x])
								if(x>Ntip(ph))
									ans	<- all(grepl(MRCA_IND, extract.clade(ph,x)$tip.label))
								ans
							})	
					#	return monophyletic clades
					list(CLADE=as.integer(pot.clade[tmp]))
				}, by=c('GROUP','CHANGE_NODE')]
		zz[, CHANGE_NODE:=NULL]
		#	merge all monophyletic clades
		z	<- merge(z, unique(zz), all=1, by='GROUP', allow.cartesian=TRUE)
	}	
	tmp	<- z[, which(MONOPH==1)]
	set(z, tmp, 'CLADE', z[tmp,MRCA])
	subset(z, select=c(IND, MRCA_ANC_EDGE, DIFF_IND, CLADE))
}

#' @export
pty.stat.avgnodedepth.and.counts<- function(ph, mrca, select=NA)
{
	if(mrca<=Ntip(ph))
	{
		phb		<- data.table(	BAM=ph$tip.label[mrca],
								IND=gsub('_read.*','',ph$tip.label[mrca]),
								COUNT=as.numeric(gsub('count_','',regmatches(ph$tip.label[mrca], regexpr('count_[0-9]+',ph$tip.label[mrca])))),											
								DEPTH=0	)
	}
	if(mrca>Ntip(ph))
	{
		ph		<- extract.clade(ph,mrca)
		phb		<- data.table(	BAM=ph$tip.label,
								IND=gsub('_read.*','',ph$tip.label),
								COUNT= as.numeric(gsub('count_','',regmatches(ph$tip.label, regexpr('count_[0-9]+',ph$tip.label))))	)
		phb[, DEPTH:=node.depth.edgelength(ph)[seq_len(Ntip(ph))]]
	}		
	if(!is.na(select))
		phb		<- subset(phb, grepl(select,IND))
	list(AVG_ROOT2TIP=phb[, sum(DEPTH*COUNT) / sum(COUNT)], READ_N=phb[, sum(COUNT)], UREAD_N=phb[, length(COUNT)])	
}

#' @export
pty.stat.superinfections.160208<- function(pty.ph, ptyfiles, references.pattern='REF')
{	
	#	@CF: with data.table, the df[, {}, by=BY] syntax specifies the columns BY to loop over
	#	@CF: the return variables of a data.table block {} are either a named list, or another data.table
	#	@CF: one difference to data.frame is that the column names are available as variable names inside the block {}
	
	#	For each phylotype run (PTY_RUN), each patient (IND), each window (W_FROM, W_TO):
	#	- find the root node (CLADE) defining all monophyletic clades		
	stat.clades	<- ptyfiles[, {
				#FILE		<- subset(ptyfiles, W_FROM==661)[,FILE]
				#FILE<- 'ptyr22_InWindow_481_to_540_dophy_examl.newick'
				cat('\n',FILE)
				ph			<- pty.ph[[FILE]]
				pty.stat.monophyletic.clades(ph)							
			}, by=c('PTY_RUN','W_FROM','W_TO','FILE')]
	gc()	
	#	For each monophyletic clade:
	#	- calculate the mean node root-to-tip within clade, weighted by read count
	#	- calculate the number of reads
	#	- calculate the number of unique reads
	tmp			<- stat.clades[, {
				#tmp<- subset(stat.clades, W_FROM==661); FILE<- tmp[,FILE]
				#FILE<- 'ptyr22_InWindow_661_to_720_dophy_examl.newick'; CLADE<- 2066
				#cat('\n',FILE)
				ph			<- pty.ph[[FILE]]
				tmp			<- pty.stat.avgnodedepth.and.counts(ph, CLADE, select=NA)	
				names(tmp)	<- paste(names(tmp),'_CLADE',sep='')				
				tmp
			}, by=c('PTY_RUN','W_FROM','W_TO','FILE','CLADE')]
	stat.clades	<- merge(stat.clades, tmp, by=c('PTY_RUN','W_FROM','W_TO','FILE','CLADE'))
	gc()	
	#	For each phylotype run (PTY_RUN), each patient (IND), each window (W_FROM, W_TO):
	#	- calculate the average root to tip distance, weighted by read count
	#	  this requires the MRCA which can be expensive but has already been computed
	#	- calculate the total number of reads 
	#	- calculate the number of unique reads
	setkey(stat.clades, PTY_RUN, W_FROM, W_TO, FILE, IND)
	tmp			<- subset(unique(stat.clades), select=c(PTY_RUN, W_FROM, W_TO, FILE, IND, MRCA))
	stat.ind	<- tmp[,{
				ph			<- pty.ph[[FILE]]
				pty.stat.avgnodedepth.and.counts(ph, MRCA, select=IND)				
			}, by=c('PTY_RUN','W_FROM','W_TO','FILE','IND')]
	set(stat.ind, NULL, 'UREAD_N', stat.ind[, as.numeric(UREAD_N)])
	#	For each phylotype run (PTY_RUN), each patient (IND), each window (W_FROM, W_TO):
	#	- calculate number of clades
	tmp			<- stat.clades[, list(CLADE_N=as.numeric(length(CLADE))), by=c('PTY_RUN','W_FROM','W_TO','FILE','IND')]	
	stat.ind	<- merge(stat.ind, tmp, by=c('PTY_RUN','W_FROM','W_TO','FILE','IND'))
	#	For each monophyletic clade:
	#	- calculate the proportion of reads in the clade
	#	- calculate clade order by proportion
	stat.clades	<- merge(stat.clades, subset(stat.ind, select=c(PTY_RUN, W_FROM, W_TO, FILE, IND, READ_N)), by=c('PTY_RUN','W_FROM','W_TO','FILE','IND'))	
	stat.clades[, READ_P_CLADE:= READ_N_CLADE/READ_N]
	stat.clades[, READ_N:= NULL]
	stat.clades	<- stat.clades[ order(PTY_RUN, W_FROM, W_TO, FILE, IND, -READ_P_CLADE), ]	
	tmp			<- stat.clades[, list(ORDER_CLADE= seq_along(CLADE), CLADE=CLADE), by=c('PTY_RUN','W_FROM','W_TO','FILE','IND')]
	stat.clades	<- merge(stat.clades, tmp, by=c('PTY_RUN','W_FROM','W_TO','FILE','IND','CLADE'))
	#
	#	check that all reads are in one clade
	#
	tmp			<- stat.clades[, list(READ_N_CHECK=sum(READ_N_CLADE)), by=c('PTY_RUN','W_FROM','W_TO','FILE','IND')]	
	tmp			<- merge(tmp, subset(stat.ind, select=c(PTY_RUN, W_FROM, W_TO, FILE, IND, READ_N)), by=c('PTY_RUN','W_FROM','W_TO','FILE','IND'))
	stopifnot( !nrow(subset(tmp, READ_N_CHECK<READ_N))	)
	list(stat.clades=stat.clades, stat.ind=stat.ind)
}

#' @export
pty.stat.superinfections.160208.plot<- function(stat.clades, stat.ind, outfile, plot.max.clade=5)
{			
	#
	#	print by patient
	#
	
	#	@CF: To get a sub-plot for each statistic using 'facet_grid', I melt the data.table
	#	@CF: the 'set' command looks a bit weird, but is great because it avoids creating new R objects via '<-'
	
	#	melt the statistics on individuals
	stat.plot	<- melt(stat.ind, measure.vars=c('READ_N','UREAD_N','CLADE_N','AVG_ROOT2TIP'), id.vars=c('PTY_RUN','W_FROM','W_TO','FILE','IND'), value.name='V', variable.name='STAT')
	#	melt the statistics on clades
	tmp			<- melt(stat.clades, measure.vars=c('READ_P_CLADE','AVG_ROOT2TIP_CLADE'), id.vars=c('PTY_RUN','W_FROM','W_TO','FILE','IND','ORDER_CLADE'), value.name='V', variable.name='STAT')
	set(tmp, NULL, 'STAT', tmp[,gsub('_CLADE','',STAT)])
	#	select most prevalent clades and put all together
	stat.plot	<- rbind( stat.plot, subset(tmp, ORDER_CLADE<=plot.max.clade), use.names=TRUE, fill=TRUE )
	#	give things nice names
	set(stat.plot, NULL, 'ORDER_CLADE', stat.plot[, as.character(ORDER_CLADE)])
	tmp			<- stat.plot[, which(!is.na(ORDER_CLADE))]
	set(stat.plot, tmp, 'ORDER_CLADE', stat.plot[tmp, paste('Clade',ORDER_CLADE)])
	tmp			<- stat.plot[, sort(na.omit(unique(ORDER_CLADE)))]
	set(stat.plot, stat.plot[, which(is.na(ORDER_CLADE))], 'ORDER_CLADE', 'Individual')
	set(stat.plot, NULL, 'ORDER_CLADE', stat.plot[, factor(ORDER_CLADE, levels=c('Individual',tmp), labels=c('Individual',tmp))])	
	set(stat.plot, NULL, 'STAT', stat.plot[,factor(STAT, levels=c('READ_N','UREAD_N','CLADE_N','READ_P','AVG_ROOT2TIP'), labels=c('Reads\n(#)', 'Tips\n(#)','Monophyletic\nclades\n(#)','Reads\n(%)','Root-to-tip\npatristic distance\n(mean subst/site)'))])
	#	
	#	setup plots per patient	
	#
	ps	<- lapply(stat.plot[, unique(IND)], function(x){
				#x<- 'R1_RES669_S20_L001'
				tmp			<- subset(stat.plot, IND==x)
				set(tmp, NULL, 'IND', tmp[,paste('run:',PTY_RUN,', individual:',IND)])
				col			<- c('black',rainbow_hcl(tmp[, length(unique(ORDER_CLADE))-1], start = 270, end = -30, c=100, l=50))
				names(col)	<- tmp[, unique(ORDER_CLADE)]	
				setkey(tmp, PTY_RUN, W_FROM, IND, ORDER_CLADE)
				p			<- ggplot(tmp , aes(x=W_FROM, y=V)) +
						geom_point(data=subset(tmp, !grepl('Reads\n(%)',STAT)), size=0.5, aes(colour=ORDER_CLADE)) +
						geom_bar(data=subset(tmp, STAT=='Reads\n(%)'), stat='identity', aes(fill=ORDER_CLADE)) +
						geom_line(data=subset(tmp, grepl('Root-to-tip',STAT)), aes(colour=ORDER_CLADE, group=ORDER_CLADE)) +
						scale_colour_manual(values=col) +
						scale_fill_manual(values=col, guide=FALSE) +
						scale_x_continuous(breaks=seq(0,12e3,1e3), minor_breaks=seq(0,12e3,2e2), expand=c(0,0), lim=c(1,stat.plot[, max(W_TO)+1])) +
						theme_bw() + labs(x='Genome location', y='', colour='', fill='') +
						facet_grid(STAT~IND, scale='free_y')				
			})
	cat('\nPlot to',outfile)
	pdf(file=outfile, w=9, h=9)			
	for(p in ps)
		print(p)					
	dev.off()	
}

#' @export
pty.stat.all.160128	<- function(pty.ph, ptyfiles)
{
	#
	#	quantiles of within individual diversity (patristic distance)
	#	calculated among unique reads, and all reads (expand by count of each read)		
	coi.div		<- ptyfiles[, {
				#FILE		<- subset(ptyfiles, W_FROM==6961)[,FILE]
				ph			<- pty.ph[[FILE]]
				pty.stat.withinhost.diversity(ph)
			}, by=c('PTY_RUN','W_FROM','W_TO','FILE')]
	coi.div[, N:=NULL]
	gc()
	tmp			<- subset(coi.div, READS=='ALL')
	setnames(tmp, c("WHD_p2","WHD_p25","WHD_p50","WHD_p75","WHD_p98"), c("WHDA_p2","WHDA_p25","WHDA_p50","WHDA_p75","WHDA_p98"))
	tmp[, READS:=NULL]
	coi.div		<- subset(coi.div, READS=='UNIQUE')
	setnames(coi.div, c("WHD_p2","WHD_p25","WHD_p50","WHD_p75","WHD_p98"), c("WHDU_p2","WHDU_p25","WHDU_p50","WHDU_p75","WHDU_p98"))
	coi.div[, READS:=NULL]
	coi.div		<- merge(coi.div, tmp, by=c('PTY_RUN','W_FROM','W_TO','FILE','FILE_ID'))
	#
	#	maximum branch length within individual, and mean path length of subtending clade
	#	calculated among unique reads, and all reads (expand by count of each read)
	coi.lsep	<- ptyfiles[, {
				#FILE		<- subset(ptyfiles, W_FROM==1021)[,FILE]						
				ph			<- pty.ph[[FILE]]
				pty.stat.maxlocalsep(ph)														
			}, by=c('PTY_RUN','W_FROM','W_TO','FILE')]
	coi.div		<- merge(coi.div, coi.lsep, by=c('PTY_RUN','W_FROM','W_TO','FILE','FILE_ID'))
	coi.lsep	<- NULL
	gc()
	#
	#	Different: number of clades from different individual
	coi.diff	<- ptyfiles[, {
				#FILE		<- subset(ptyfiles, W_FROM==1621)[,FILE]
				ph			<- pty.ph[[FILE]]
				pty.stat.different(ph)								
			}, by=c('PTY_RUN','W_FROM','W_TO','FILE')]
	pty.stat	<- merge(coi.div, coi.diff, by=c('PTY_RUN','W_FROM','W_TO','FILE','FILE_ID'))
	coi.diff	<- coi.div	<- NULL
	gc()
	pty.stat
}

#' @export
pty.evaluate.tree.collapse.clusters<- function(ph, thresh.brl=8e-6)
{
	#	cluster by genetic distance
	dist.brl	<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)							
	clustering	<- hivc.clu.clusterbythresh(ph, thresh.brl=thresh.brl, dist.brl=dist.brl, retval="all")			
	#print(clustering)	
	df			<- data.table(BAM=ph$tip.label, FILE_ID= gsub('_read.*','',ph$tip.label), COUNT= as.numeric(gsub('count_','',regmatches(ph$tip.label, regexpr('count_[0-9]+',ph$tip.label)))),  CLU_ID=clustering$clu.mem[seq_len(Ntip(ph))])
	#df			<- merge(df, data.table(CLU_ID=seq_along(clustering$clu.idx), MRCA=clustering$clu.idx), by='CLU_ID', all.x=1)
	df			<- subset(df,!is.na(CLU_ID))
	phc			<- copy(ph)
	if(nrow(df))
	{
		#	only collapse sequences from same individual
		tmp		<- subset(df[, list(CLU_N=length(unique(FILE_ID))), by='CLU_ID'], CLU_N==1)
		df		<- merge(df, tmp, by='CLU_ID')			
		df		<- merge(df, df[, list(TAXA_N=length(FILE_ID), BAMCLU= paste( FILE_ID[1],'_clu_',CLU_ID,'_count_',sum(COUNT), sep='' ) ), by='CLU_ID'], by='CLU_ID')
		#	do not collapse if this results in singleton
		if(Ntip(phc)==nrow(df) && df[, length(unique(CLU_ID))]==1)
			df	<- subset(df, CLU_ID<0)
		for(clu.id in df[, unique(CLU_ID)])
		{
			tmp			<- subset(df, CLU_ID==clu.id)			
			z			<- drop.tip(phc, tip=tmp[, BAM], subtree=TRUE)
			if(!is.null(z))
			{
				phc		<- z
				phc$tip.label[ which(grepl("[",phc$tip.label,fixed=TRUE)) ]	<- tmp[1,BAMCLU]	
			}			 				
		}
	}	
	phc
}

pty.evaluate.tree.groupindividuals<- function(ph, phb)
{
	tmp				<- lapply( phb[, unique(FILE_ID)], function(x)	subset(phb, FILE_ID==x)[, BAM]	)
	names(tmp)		<- phb[, unique(FILE_ID)]
	ph				<- groupOTU(ph, tmp, group='INDIVIDUAL')		
	z	<- merge(data.table(FROM=ph$edge[,1],IDX=ph$edge[,2]), phb, by='IDX', all=1)
	z[, GROUP:= attr(ph,'INDIVIDUAL')[1:nrow(ph$edge)]]
	z	<- unique(subset(z, !is.na(FILE_ID), select=c(FILE_ID, GROUP)))
	attr(ph,'INDIVIDUAL')	<- factor(attr(ph,'INDIVIDUAL'), levels=c(0,z[,as.character(GROUP)]), labels=c('not characterized',z[,FILE_ID]))
	ph
}

pty.evaluate.tree.groupcandidates<- function(ph, phb)
{
	set(phb, NULL, 'FILL', phb[, factor(FILL, levels=c(0,1), labels=c('target','filler'))])
	tmp				<- lapply( phb[, unique(FILL)], function(x)	subset(phb, FILL==x)[, BAM]	)
	names(tmp)		<- phb[, unique(FILL)]
	ph				<- groupOTU(ph, tmp, group='TYPE')				
	z	<- merge(data.table(FROM=ph$edge[,1],IDX=ph$edge[,2]), phb, by='IDX', all=1)
	z[, GROUP:= attr(ph,'TYPE')[1:nrow(ph$edge)]]
	z	<- unique(subset(z, !is.na(FILE_ID), select=c(FILL, GROUP)))
	attr(ph,'TYPE')	<- factor(attr(ph,'TYPE'), levels=c(0,z[,as.character(GROUP)]), labels=c('not characterized',z[,as.character(FILL)]))
	ph
}

#' @export
pty.evaluate.tree.collapse<- function(pty.runs, ptyfiles, pty.phc, outdir, thresh.brl=8e-6)
{	
	#pty.phc<- copy(pty.ph)
	#length(pty.phc)<- 10
	#	collapse nodes to make browsing of plots easier
	cat('\ncalculating cladograms')
	for(i in seq_along(pty.phc))
	{
		print(i) 
		ph				<- pty.phc[[i]]		
		#FILE<- "ptyr1_InWindow_481_to_540_dophy_examl.newick"
		#ph		<- pty.phc[[FILE]]		
		ph				<- pty.evaluate.tree.collapse.clusters(ph, thresh.brl=thresh.brl)
		phb				<- data.table(IDX=seq_along(ph$tip.label), BAM=ph$tip.label, FILE_ID= gsub('_read.*|_clu.*','',ph$tip.label))
		#	group edges by individual
		ph				<- pty.evaluate.tree.groupindividuals(ph, phb)
		#	group edges by FILL
		tmp				<- as.numeric(gsub('ptyr','',regmatches(names(pty.ph)[i], regexpr('ptyr[0-9]+',names(pty.ph)[i]))))
		phb				<- merge(phb, subset(pty.runs, PTY_RUN==tmp), by='FILE_ID')
		ph				<- pty.evaluate.tree.groupcandidates(ph, phb)	
		pty.phc[[i]]	<- ph		
	}
	#	save
	tmp		<- file.path( outdir, ptyfiles[1, gsub('\\.newick','\\collapsed.rda',gsub('_dophy','',gsub('_InWindow_[0-9]+_to_[0-9]+','',FILE)))] )
	cat('\nsave to file', tmp)
	save(pty.phc, file=tmp)
	#	need node heights for plotting
	tmp				<- ptyfiles[,	{
								#print(FILE)
								#FILE<- "ptyr1_InWindow_481_to_540_dophy_examl.newick"
								ph		<- pty.phc[[FILE]]
								tmp		<- node.depth.edgelength(ph)[1:Ntip(ph)]
								list(BAM=ph$tip.label, HEIGHT=tmp)
							}, by='FILE']	
	ptyfiles		<- merge(ptyfiles, tmp[, list(HMX= max(HEIGHT)), by='FILE'], by='FILE')
	setkey(ptyfiles, PTY_RUN, W_FROM)	
	#ptyfiles			<- subset(ptyfiles, PTY_RUN==1)
	for(ptyr in ptyfiles[, unique(PTY_RUN)])
	{
		#ptyr<- 15
		cat('\nplotting cladograms to R for',ptyr)
		tmp			<- subset(ptyfiles, PTY_RUN==ptyr)
		#	title
		tmp[, TITLE:=paste('run',PTY_RUN,', window [',W_FROM,',',W_TO,']',sep='')]
		setkey(tmp, W_FROM)
		phs			<- lapply(tmp[, FILE], function(x) pty.phc[[x]])
		names(phs)	<- tmp[, TITLE] 
		#	colours	
		tmp2		<- setdiff(sort(unique(unlist(lapply(seq_along(phs), function(i)	levels(attr(phs[[i]],'INDIVIDUAL'))	)))),'not characterized')
		col			<- c('black',rainbow_hcl(length(tmp2), start = 270, end = -30, c=100, l=50))
		names(col)	<- c('not characterized',tmp2)	
		phps		<- lapply(seq_along(phs), function(i){					
					max.node.height	<- tmp[i,][, HMX]
					df				<- data.table(	BAM=phs[[i]]$tip.label, IDX=seq_along(phs[[i]]$tip.label), 
													COUNT= as.numeric(gsub('count_','',regmatches(phs[[i]]$tip.label, regexpr('count_[0-9]+',phs[[i]]$tip.label)))), 
													CLU=grepl('_clu',phs[[i]]$tip.label), 
													FILE_ID= gsub('_read.*|_clu.*','',phs[[i]]$tip.label))
					p				<- ggtree(phs[[i]], aes(color=INDIVIDUAL, linetype=TYPE)) %<+% df + 
													geom_nodepoint(size=phs[[i]]$node.label/100*3) +
													geom_tippoint(aes(size=COUNT, shape=CLU)) +
													geom_tiplab(size=1.2,  hjust=-.1) +							 
													scale_color_manual(values=col, guide = FALSE) +
													scale_shape_manual(values=c(20,18), guide=FALSE) +
													scale_size_area(guide=FALSE) +							
													scale_linetype_manual(values=c('target'='solid','filler'='dotted'),guide = FALSE) +
													theme_tree2() +
													theme(legend.position="bottom") + ggplot2::xlim(0, max.node.height*1.3) +
													labs(x='subst/site', title=names(phs)[i])
					p					
				})	
		names(phps)	<- names(phs)
		file	<- file.path( outdir, tmp[1,gsub('.newick','collapsed.pdf',gsub('_dophy','',gsub('_InWindow_[0-9]+_to_[0-9]+','',FILE)))] )
		cat('\nplotting cladograms to file',file)
		if(1)		#for window length 60 (multiple pages)
		{				
			tmp			<- seq_len(ceiling(length(phps)/10))
			pdf(file=file, w=20, h=40)		#for win=60
			for(i in tmp)
			{		
				grid.newpage()
				pushViewport(viewport(layout=grid.layout(2, 5)))
				z	<- intersect(seq.int((i-1)*10+1, i*10), seq_len(length(phps)))
				for(j in z)
					print(phps[[j]], vp = viewport(layout.pos.row=(ceiling(j/5)-1)%%2+1, layout.pos.col=(j-1)%%5+1))				
			}
			dev.off()	
		}	
	}
}

#' @import ape data.table gridExtra colorspace ggtree
#' @export
pty.evaluate.tree<- function(indir, pty.runs=NULL, outdir=indir, select='', outgroup=NA, references.pattern='REF', run.pattern='ptyr', tree.pattern='newick$', rm.newick=FALSE, rm.fasta=FALSE, plot.trees.per.page=10, plot.w=20, plot.h=40)
{	
	if(0)
	{
		indir		<- "~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/coinf_ptoutput_150121"
		indir		<- "~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/coinf_ptoutput_150201"
		outdir		<- indir
		pty.infile	<- "~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_coinfrunsinput.rda"
		select		<- '^ptyr22_In'
		outgroup	<- NA
		outgroup	<- "CPX_AF460972"
	}
	stopifnot(!is.null(pty.runs) || (is.null(pty.runs) & !is.na(outgroup)))
	#
	#	collect ML tree files
	#
	ptyfiles		<- data.table(FILE=list.files(indir, tree.pattern))
	if(nchar(run.pattern))
		ptyfiles[, PTY_RUN:= as.numeric(gsub(run.pattern,'',sapply(strsplit(FILE,'_'),'[[',1)))]
	if(!nchar(run.pattern))
		ptyfiles[, PTY_RUN:=1L]
	ptyfiles[, W_FROM:= as.numeric(gsub('InWindow_','',regmatches(FILE,regexpr('InWindow_[0-9]+',FILE))))] 
	ptyfiles[, W_TO:= as.numeric(gsub('to_','',regmatches(FILE,regexpr('to_[0-9]+',FILE))))]
	ptyfiles		<- subset(ptyfiles, grepl(select,FILE))
	setkey(ptyfiles, PTY_RUN, W_FROM)
	cat('\nfound tree files, n=',nrow(ptyfiles))
	#ptyfiles		<- subset(ptyfiles, W_FROM<9000)
	#
	#ptyfiles[, table(PTY_RUN)]
	#	raw trees w/o any attributes
	pty.ph		<- lapply( seq_len(nrow(ptyfiles)), function(i)
			{				
				ph			<- read.tree(file.path(indir,ptyfiles[i, FILE]))
				#	node labels
				if(is.null(ph$node.label))
					ph$node.label	<- rep('0',Nnode(ph))
				tmp				<- ph$node.label
				tmp[which(tmp=='Root'|tmp=='')]	<- '0'
				ph$node.label	<- as.numeric(tmp)
				ph
			})
	names(pty.ph)<- ptyfiles[, FILE]
	#
	ptyfiles[, IND_N:= sapply(seq_along(pty.ph), function(i)  length(unique(gsub('_read.*','',pty.ph[[i]]$tip.label)))		)]
	#	rm trees with just one individual if no outgroup specified	
	if(is.na(outgroup) & pty.runs[, any(FILL==1)])
	{
		cat('\nignoring trees with only one individual (only when no root specified)', subset(ptyfiles, IND_N==1)[, paste(unique(FILE),collapse=', ')])
		ptyfiles		<- subset(ptyfiles, IND_N>1)
		pty.ph			<- lapply( ptyfiles[, FILE], function(x)	pty.ph[[x]]		)
		names(pty.ph)	<- ptyfiles[, FILE]			
	}	
	ptyfiles[, IDX:=seq_len(nrow(ptyfiles))]
	#
	#	GET ROOTS
	#	
	#	root specified; check if in tree
	if(!is.na(outgroup))
	{
		cat('\nroot and evaluate trees')
		pty.root	<- ptyfiles[, {
					#FILE<- 'ptyr22_InWindow_241_to_300_dophy_examl.newick'
					ans		<- NA_character_
					tmp		<- which(grepl(outgroup,pty.ph[[FILE]]$tip.label))
					if(length(tmp)==0)
						cat('\ncannot find root', outgroup,'in tree',FILE,', ignoring file')
					if(length(tmp)==1)
						ans	<- pty.ph[[FILE]]$tip.label[tmp]
					if(length(tmp)>1)
						cat('\nfound >1 roots', outgroup,'in tree',FILE,', ignoring file')
					list(ROOT=ans)
				}, by= c('FILE','PTY_RUN')]
		pty.root	<- subset(pty.root, !is.na(ROOT))	
		ptyfiles	<- merge(ptyfiles,pty.root, by=c('FILE','PTY_RUN'))
	}
	#	root not specified; determine root for each run: find taxon with largest distance from BAM of selected individuals (these have FILL==0)
	if(is.na(outgroup) & pty.runs[, any(FILL==1)])
	{
		cat('\ndetermine root among fillers')
		pty.root	<- pty.evaluate.tree.root.withfill(ptyfiles, pty.ph, pty.runs)
		ptyfiles	<- merge(ptyfiles,pty.root, by=c('FILE','PTY_RUN'))
	}
	#	root not specified; determine root for each run: 
	if(is.na(outgroup) & pty.runs[, all(FILL==0)])
	{
		cat('\ndetermine root among all individuals')
		pty.root	<- pty.evaluate.tree.root.nofill(ptyfiles, pty.ph, pty.runs)
		ptyfiles	<- merge(ptyfiles,pty.root, by=c('FILE','PTY_RUN'))
	}	
	#
	#	trees: root, ladderize, group
	#
	#pty.ph.cp			<- copy(pty.ph)	
	#	pty.ph	<- copy(pty.ph.cp)
	cat('\nroot, ladderize, group trees')
	for(i in seq_along(pty.ph))
	{
		cat('\nprocess tree ',names(pty.ph)[i])
		#print(i)
		#i<- 44
		#i<- 5
		root			<- subset(ptyfiles, FILE==names(pty.ph)[i])[, ROOT]
		ph				<- pty.ph[[i]]
		tmp										<- which(ph$tip.label==root)
		ph										<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
		ph$node.label[ph$node.label=='Root']	<- 0
		ph$node.label							<- as.numeric(ph$node.label)			
		ph				<- ladderize(ph)
		phb				<- data.table(IDX=seq_along(ph$tip.label), BAM=ph$tip.label, FILE_ID= gsub('_read.*','',ph$tip.label), REF=grepl(references.pattern,ph$tip.label))
		set(phb, phb[, which(REF)],'FILE_ID','REFERENCE')
		tmp				<- phb[, which(!grepl('_read_[0-9]+_count_[0-9]+',BAM))]
		#	homogenize taxon labels
		if(length(tmp))
		{
			stopifnot( !nrow(subset(phb[tmp, ],!REF))	)	#incorrect taxon label
			set(phb, tmp, 'BAM', phb[tmp, paste(BAM,'_read_1_count_0',sep='')])
			setkey(phb, IDX)
			ph$tip.label	<- phb[,BAM]			
		}
		#	group edges by individual
		ph				<- pty.evaluate.tree.groupindividuals(ph, phb)
		#	group edges by FILL
		#tmp			<- as.numeric(gsub('ptyr','',regmatches(names(pty.ph)[i], regexpr('ptyr[0-9]+',names(pty.ph)[i]))))
		#phb			<- merge(phb, subset(pty.runs, PTY_RUN==tmp), by='FILE_ID')
		#ph				<- pty.evaluate.tree.groupindividuals(ph, phb)	
		pty.ph[[i]]		<- ph
	}		
	#
	#	save trees
	#	
	tmp		<- ptyfiles[1, gsub(paste('.',tree.pattern,sep=''),'_preprtr.rda',gsub('_dophy','',gsub('\\.InWindow_[0-9]+_to_[0-9]+|_InWindow_[0-9]+_to_[0-9]+','',FILE)))]
	cat('\nsave trees to file ',tmp)
	save(pty.ph, ptyfiles, file=file.path(outdir,tmp))
	#
	#	remove newick / fasta files
	#
	if(rm.newick)
	{
		invisible( ptyfiles[, list(SUCCESS=file.remove(file.path(indir,FILE))), by='FILE'] )
		invisible( ptyfiles[, list(SUCCESS=file.remove(file.path(indir,gsub('newick','txt',FILE)))), by='FILE'] )
	}
	if(rm.fasta)
	{
		invisible( ptyfiles[, list(SUCCESS=file.remove(file.path(indir,gsub('_examl\\.newick','\\.fasta',FILE)))), by='FILE'] )
	}
	#
	#	plot full trees
	#
	#	need node heights to determine x-axis for plotting
	tmp					<- ptyfiles[,{										
				ph		<- pty.ph[[FILE]]
				tmp		<- node.depth.edgelength(ph)[1:Ntip(ph)]
				list(BAM=ph$tip.label, HEIGHT=tmp)
			}, by='FILE']	
	ptyfiles				<- merge(ptyfiles, tmp[, list(HMX= max(HEIGHT)), by='FILE'], by='FILE')
	setkey(ptyfiles, PTY_RUN, W_FROM)	
	#	plot by run.pattern
	#ptyfiles			<- subset(ptyfiles, PTY_RUN==1)
	for(ptyr in ptyfiles[, unique(PTY_RUN)])
	{
		#ptyr<- 22
		cat('\nplot trees with run.pattern',ptyr,'\n\n')
		tmp			<- subset(ptyfiles, PTY_RUN==ptyr)
		#	setup title
		tmp[, TITLE:=paste('run',PTY_RUN,', window [',W_FROM,',',W_TO,']',sep='')]
		setkey(tmp, W_FROM)
		phs			<- lapply(tmp[, FILE], function(x) pty.ph[[x]])
		names(phs)	<- tmp[, TITLE] 
		#	setup colours	
		tmp2		<- setdiff(sort(unique(unlist(lapply(seq_along(phs), function(i)	levels(attr(phs[[i]],'INDIVIDUAL'))	)))),c('not characterized'))
		if(!'REFERENCE'%in%tmp2)
		{
			col			<- c('black',rainbow_hcl(length(tmp2), start = 270, end = -30, c=100, l=50))
			names(col)	<- c('not characterized',tmp2)						
		}
		if('REFERENCE'%in%tmp2)
		{
			tmp2		<- setdiff(tmp2,'REFERENCE')
			col			<- c('grey50','black',rainbow_hcl(length(tmp2), start = 270, end = -30, c=100, l=50))
			names(col)	<- c('not characterized','REFERENCE',tmp2)						
		}		
		#	setup ggtrees
		phps		<- lapply(seq_along(phs), function(i){
					cat('\nsetup tree ',names(phs)[i])
					max.node.height	<- tmp[i,][, HMX]
					ph				<- phs[[i]]					
					df			<- data.table(	BAM=ph$tip.label, IDX=seq_along(ph$tip.label), 
							COUNT= as.numeric(gsub('count_','',regmatches(ph$tip.label, regexpr('count_[0-9]+',ph$tip.label)))), 
							CLU=grepl('_clu',ph$tip.label), 
							FILE_ID= gsub('_read.*|_clu.*','',ph$tip.label))
					set(df, df[, which(grepl(references.pattern,FILE_ID))],'FILE_ID','REFERENCE')
					p			<- ggtree(ph, aes(color=INDIVIDUAL)) %<+% df +
							geom_nodepoint(size=ph$node.label/100*3) +
							geom_tippoint(aes(size=COUNT, shape=CLU)) +
							#geom_point2(aes(subset=(node==42)), size=5, shape=23, fill='#068C00') +
							#geom_point2(data=df, aes(subset=(node==IDX)), size=5, shape=23, fill=df[,COL]) +
							#geom_text(aes(label=label), size=1.5,  hjust=-.1) +
							geom_tiplab(size=1.2,  hjust=-.1) +
							scale_color_manual(values=col, guide=FALSE) +
							scale_shape_manual(values=c(20,18), guide=FALSE) +
							scale_size_area(guide=FALSE) +							
							theme_tree2() +
							theme(legend.position="bottom") + ggplot2::xlim(0, max.node.height*1.3) +
							labs(x='subst/site', title=names(phs)[i])
					p
				})	
		names(phps)	<- names(phs)
		file	<- file.path( indir, tmp[1,gsub(paste('.',tree.pattern,sep=''),'.pdf',gsub('_dophy','',gsub('.InWindow_[0-9]+_to_[0-9]+|_InWindow_[0-9]+_to_[0-9]+','',FILE)))] )
		cat('\nplot trees to file ',file)		
		tmp			<- seq_len(ceiling(length(phps)/plot.trees.per.page))
		#	single-page plot
		if(plot.trees.per.page==1)
		{
			pdf(file=file, w=plot.w, h=plot.h)		#for win=60
			for(i in tmp)
				print(phps[[i]])
			dev.off()
		}
		#	multi-page plot
		if(plot.trees.per.page>1)
		{
			stopifnot( plot.trees.per.page%%2==0 )			
			pdf(file=file, w=plot.w, h=plot.h)		#for win=60
			for(i in tmp)
			{		
				grid.newpage()
				pushViewport(viewport(layout=grid.layout(2, plot.trees.per.page/2)))
				z	<- intersect(seq.int((i-1)*plot.trees.per.page+1, i*plot.trees.per.page), seq_len(length(phps)))
				for(j in z)
					print(phps[[j]], vp = viewport(layout.pos.row=(ceiling(j/(plot.trees.per.page/2))-1)%%2+1, layout.pos.col=(j-1)%%(plot.trees.per.page/2)+1))				
			}
			dev.off()
		}					
	}
}	

pty.evaluate.tree.plot.develstuff<- function()
{
	if(0)	#devel
	{
		tmp	 	<- matrix(seq(1, 2*ceiling(length(phps)/2)), ncol=ceiling(length(phps)/2), nrow=2)
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(tmp), ncol(tmp))))
		for(i in seq_along(phps))
			print(phps[[i]], vp = viewport(layout.pos.row=which(tmp==i, arr.ind=TRUE)[1,'row'], layout.pos.col=which(tmp==i, arr.ind=TRUE)[1,'col']))	
		dev.off()
	}
	if(0)	#devel
	{	
		ph				<- pty.ph[[30]]
		ph				<- pty.evaluate.tree.collapse.clusters(ph, thresh.brl=8e-6)
		phb				<- data.table(IDX=seq_along(ph$tip.label), BAM=ph$tip.label, FILE_ID= gsub('_read.*|_clu.*','',ph$tip.label))
		#	group edges by individual
		ph				<- pty.evaluate.tree.groupindividuals(ph, phb)
		#	group edges by FILL
		tmp				<- as.numeric(gsub('ptyr','',regmatches(names(pty.ph)[i], regexpr('ptyr[0-9]+',names(pty.ph)[i]))))
		phb				<- merge(phb, subset(pty.runs, PTY_RUN==tmp), by='FILE_ID')
		ph				<- pty.evaluate.tree.groupcandidates(ph, phb)	
		pty.phc[[i]]	<- ph
		
		
		tmp2		<- levels(attr(ph,'INDIVIDUAL'))
		col			<- c('black',rainbow_hcl(length(tmp2)-1, start = 270, end = -30, c=100, l=50))
		names(col)	<- tmp2	
		df			<- data.table(	BAM=ph$tip.label, IDX=seq_along(ph$tip.label), 
				COUNT= as.numeric(gsub('count_','',regmatches(ph$tip.label, regexpr('count_[0-9]+',ph$tip.label)))), 
				CLU=grepl('_clu',ph$tip.label), 
				FILE_ID= gsub('_read.*|_clu.*','',ph$tip.label))
		set(df, df[, which(grepl(references.pattern,FILE_ID))],'FILE_ID','REFERENCE')
		ggtree(ph, aes(color=INDIVIDUAL)) %<+% df +
				geom_nodepoint(size=ph$node.label/100*3) +
				geom_tippoint(aes(size=COUNT, shape=CLU)) +
				#geom_point2(aes(subset=(node==42)), size=5, shape=23, fill='#068C00') +
				#geom_point2(data=df, aes(subset=(node==IDX)), size=5, shape=23, fill=df[,COL]) +
				#geom_text(aes(label=label), size=1.5,  hjust=-.1) +
				geom_tiplab(size=1.2,  hjust=-.1) +
				scale_color_manual(values=col, guide=FALSE) +
				scale_shape_manual(values=c(20,18), guide=FALSE) +
				scale_size_area(guide=FALSE) +
				#scale_linetype_manual(values=c('target'='solid','filler'='dotted'),guide = FALSE) +
				theme_tree2() +
				theme(legend.position="bottom") + ggplot2::xlim(0, 0.3) +
				labs(x='subst/site')
		
	}
}

pty.drop.tip<- function (phy, tip, trim.internal = TRUE, subtree = FALSE, root.edge = 0, 
		rooted = is.rooted(phy), interactive = FALSE, subtree.label=NA) 
{
	if (!inherits(phy, "phylo")) 
		stop("object \"phy\" is not of class \"phylo\"")
	Ntip <- length(phy$tip.label)
	if (interactive) {
		cat("Left-click close to the tips you want to drop; right-click when finished...\n")
		xy <- locator()
		nToDrop <- length(xy$x)
		tip <- integer(nToDrop)
		lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
		for (i in 1:nToDrop) {
			d <- sqrt((xy$x[i] - lastPP$xx)^2 + (xy$y[i] - lastPP$yy)^2)
			tip[i] <- which.min(d)
		}
	}
	else {
		if (is.character(tip)) 
			tip <- which(phy$tip.label %in% tip)
	}
	out.of.range <- tip > Ntip
	if (any(out.of.range)) {
		warning("some tip numbers were larger than the number of tips: they were ignored")
		tip <- tip[!out.of.range]
	}
	if (!length(tip)) 
		return(phy)
	if (length(tip) == Ntip) {
		warning("drop all tips of the tree: returning NULL")
		return(NULL)
	}
	wbl <- !is.null(phy$edge.length)
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
		tr <- reorder(phy, "postorder")
		N <- .C(node_depth, as.integer(Ntip), as.integer(Nnode), 
				as.integer(tr$edge[, 1]), as.integer(tr$edge[, 2]), 
				as.integer(Nedge), double(Ntip + Nnode), 1L)[[6]]
	}
	edge1 <- phy$edge[, 1]
	edge2 <- phy$edge[, 2]
	keep <- !logical(Nedge)
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
	if (subtree) {
		i <- which(tip %in% oldNo.ofNewTips)
		if (length(i)) {
			phy$tip.label[tip[i]] <- "[1_tip]"
			tip <- tip[-i]
		}
	}
	n <- length(oldNo.ofNewTips)
	phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
	phy$tip.label <- phy$tip.label[-tip]
	if (subtree || !trim.internal) {
		node2tip <- oldNo.ofNewTips[oldNo.ofNewTips > Ntip]
		new.tip.label <- if (subtree) {
					paste("[", subtree.label[node2tip], "_tips]", sep = "")
				}
				else {
					if (is.null(phy$node.label)) 
						rep("NA", length(node2tip))
					else phy$node.label[node2tip - Ntip]
				}
		phy$tip.label <- c(phy$tip.label, new.tip.label)
	}
	phy$Nnode <- dim(phy$edge)[1] - n + 1L
	newNb <- integer(Ntip + Nnode)
	newNb[NEWROOT] <- n + 1L
	sndcol <- phy$edge[, 2] > n
	phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]] <- (n + 
				2):(n + phy$Nnode)
	phy$edge[, 1] <- newNb[phy$edge[, 1]]
	storage.mode(phy$edge) <- "integer"
	if (!is.null(phy$node.label)) 
		phy$node.label <- phy$node.label[which(newNb > 0) - Ntip]
	collapse.singles(phy)
}

#' @import ape zoo data.table
#' @export
pty.evaluate.fasta<- function(indir, outdir=indir, strip.max.len=Inf, select='', min.ureads.individual=NA, verbose=1)
{
	#select			<- '^ptyr1_InWindow'
	#strip.max.len	<- 350		
	infiles			<- data.table(FILE=list.files(indir, pattern='fasta$'))
	infiles			<- subset(infiles, !grepl('*',FILE,fixed=1))
	infiles			<- subset(infiles, grepl(select,FILE))
	set(infiles, NULL, 'PTY_RUN', as.numeric(gsub('ptyr','',sapply(strsplit(infiles$FILE,'_'),'[[',1))))
	set(infiles, NULL, 'W_FROM', as.numeric(gsub('InWindow_','',regmatches(infiles$FILE,regexpr('InWindow_[0-9]+',infiles$FILE))))) 
	set(infiles, NULL, 'W_TO', as.numeric(gsub('to_','',regmatches(infiles$FILE,regexpr('to_[0-9]+',infiles$FILE)))))	
	setkey(infiles, W_FROM)	
	#	read all fasta into R
	pty.seq.rw			<- lapply( infiles[, FILE], function(x) read.dna(file.path(indir,x), format='fasta')	)
	names(pty.seq.rw)	<- infiles[, FILE]
	pty.seq				<- lapply( infiles[, FILE], function(x) seq.strip.gap(pty.seq.rw[[x]], strip.max.len=strip.max.len)	)
	names(pty.seq)		<- infiles[, FILE]
	#	get statistics
	if(verbose==1)
		cat('\ncollecting alignment statistics\t',select,'\n')
	seqd	<- infiles[,{
				if(verbose==1)	cat('\n',FILE)
				#FILE	<- 'ptyr1_InWindow_1_to_60.fasta'
				seq		<- pty.seq[[FILE]]					
				tmp		<- rownames(seq)
				seqd	<- data.table(READ=tmp, FILE_ID= gsub('_read.*','',tmp)) 
				seqd[, READ_N:=as.integer(gsub('_count_','',regmatches(tmp,regexpr('_count_[0-9]+',tmp))))]			
				tmp		<- as.character(seq)
				seqd[, FIRST:= apply( tmp, 1, function(x) which(x!='-')[1] )]
				seqd[, LAST:= ncol(tmp)-apply( tmp, 1, function(x) which(rev(x)!='-')[1] ) + 1L]
				seqd[, LEN:= LAST-FIRST+1L]
				seqd	<- merge(seqd, seqd[, list(UNIQUE_N=length(READ)),by='FILE_ID'], by='FILE_ID')
				seqd[, TOTAL_N:= nrow(seqd)]
				seqd
			},by='FILE']
	seqd	<- merge(infiles, seqd, by='FILE')
	#tmp		<- subset(si, select=c(SANGER_ID, PANGEA_ID))
	#set(tmp, NULL, 'PANGEA_ID', tmp[, gsub('-','_',PANGEA_ID)])
	#setnames(tmp, c('PANGEA_ID','SANGER_ID'), c('TAXA','FILE_ID'))
	#seqd	<- merge(seqd, tmp, by='FILE_ID',all.x=1)
	#tmp		<- seqd[, which(is.na(TAXA))]
	#set(seqd, tmp, 'TAXA', seqd[tmp, FILE_ID])
	#seqd	<- merge(seqd, subset(pty.runs, select=c(FILE_ID,PTY_RUN,FILL)), by=c('PTY_RUN','FILE_ID'), all.x=1)
	#set(seqd, NULL, 'FILL', seqd[, factor(FILL, levels=c(0,1), labels=c('candidate','filler'))])	
	tmp		<- file.path(outdir, gsub('.fasta','_alignments.rda',gsub('_InWindow_[0-9]+_to_[0-9]+','',infiles[1,][,FILE])))
	if(verbose)
		cat('\nsave to',tmp)
	#	save and delete old fasta files
	save(pty.seq.rw, pty.seq, seqd, file=tmp)	
	invisible( seqd[, file.remove(file.path(indir,FILE))] )
	#
	#	prepare examl files
	#
	setkey(seqd, W_FROM)
	#	select
	if(!is.na(min.ureads.individual))	#	select individuals with x unique reads in each window
	{
		seqd	<- subset(seqd, grepl('REF',FILE_ID) || UNIQUE_N>=min.ureads.individual)
		seqd[, TOTAL_N:=NULL]
		seqd	<- merge(seqd, seqd[, list(TOTAL_N=length(READ)), by='FILE'], by='FILE')
	}
	#	write fasta
	pty.fa	<- seqd[, {	
				#FILE<- 'ptyr1_InWindow_1_to_60.fasta'
				#READ	<- tmp$READ[which(tmp$FILE==FILE)]
				z		<- file.path(outdir, gsub('\\.fasta','_dophy\\.fasta',FILE))
				if(!file.exists(z))
				{
					cat('\nwrite to',z)
					write.dna(pty.seq[[FILE]][READ,], file=z, format='fasta', colsep='', nbcol=-1)
				}										
				list(FILE=z, TAXA_N=length(READ), PTY_RUN=PTY_RUN[1], W_FROM=W_FROM[1], W_TO=W_TO[1])
			}, by='FILE']
	cat('\nwrote fasta files for trees, min taxa=',pty.fa[,min(TAXA_N)],'median taxa=',pty.fa[,median(TAXA_N)],'max taxa=',pty.fa[,max(TAXA_N)])
}

#' @export
project.dualinfecions.phylotypes.countbam.150120<- function()
{
	require(ggplot2)
	require(data.table)
	require(Rsamtools)
	pty.infile	<- file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_ptyrunsinput.rda")
	pty.data.dir	<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'
	#pty.data.dir	<- '/work/or105/PANGEA_mapout/data'
		
	load( pty.infile )	
	tmp			<- subset(si, select=c(SANGER_ID, PANGEA_ID))
	set(tmp, NULL, 'PANGEA_ID', tmp[, gsub('-','_',PANGEA_ID)])
	setnames(tmp, c('PANGEA_ID','SANGER_ID'), c('TAXA','FILE_ID'))
	pty.runs	<- merge(pty.runs, tmp, by='TAXA', all.x=1)
	tmp			<- pty.runs[, which(is.na(FILE_ID))]
	set(pty.runs, tmp,'FILE_ID', pty.runs[tmp, TAXA])	
	tmp			<- subset(pty.runs, FILL==0)	
	pty.runs	<- subset(pty.runs, FILL==1)
	setkey(pty.runs, FILE_ID)
	pty.runs	<- unique(pty.runs)
	pty.runs	<- pty.runs[ setdiff(pty.runs$FILE_ID, tmp$FILE_ID), ] 
	pty.runs	<- rbind(tmp, pty.runs)
	
	bfiles			<- data.table(FILE=list.files(pty.data.dir, pattern='bam$'))	
	bfiles			<- subset(bfiles, !grepl('Contam',FILE))
	#	get lengths of all reads in quality trimmed bam file
	#z	<- scanBam(file.path(pty.data.dir,FILE), param=ScanBamParam(what=c('qwidth','qual')))
	bam.len			<- bfiles[,{
				z	<- scanBam(file.path(pty.data.dir,FILE), param=ScanBamParam(what=c('qwidth')))
				list(QU=seq(0,320,20), CDF=ecdf(z[[1]][['qwidth']])(seq(0,320,20)))
				#list(PR=seq(0.01,0.99,0.01), QU=quantile(z[[1]][['qwidth']], p=seq(0.01,0.99,0.01)))				
			}, by='FILE']
	bam.len[, FILE_ID:= gsub('.bam','',FILE)]
	bam.len		<- merge(bam.len, subset(pty.runs, select=c(TAXA, FILL, FILE_ID)), by='FILE_ID', all.x=1 )
	set(bam.len, bam.len[, which(is.na(FILL))], 'FILL', 2)
	set(bam.len, NULL, 'FILL', bam.len[, factor(FILL, levels=c(0,1, 2), labels=c('candidate','filler','no consensus'))])
	setnames(bam.len, 'FILL', 'TYPE')
	save(bam.len, file= gsub('ptyrunsinput','bamlen',pty.infile))	
	#
	ggplot(bam.len, aes(y=CDF, x=factor(QU), fill=TYPE)) + geom_boxplot() + 
			theme_bw() + labs(y='cumulative frequency\nin one individual\n', x='\nlength of quality-trimmed short reads\n(nt)', fill='individual') +
			theme(legend.position='bottom')
	ggsave(file= gsub('ptyrunsinput.rda','bamlen.pdf',pty.infile), w=10, h=7)	
}

#' @export
pty.cmdwrap.examl<- function(pty.args)
{
	indir					<- pty.args[['out.dir']]
	outdir					<- indir	
	stopifnot( pty.args[['exa.n.per.run']]>=0, pty.args[['bs.n']]>=0, !is.na(pty.args[['min.ureads.individual']]) | !is.na(pty.args[['min.ureads.candidate']])	)	
	infiles		<- data.table(FILE=list.files(indir, pattern='_alignments.rda$'))
	infiles[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]	
	tmp				<- data.table(OUTFILE=list.files(indir, pattern='.newick$'))						
	tmp[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(OUTFILE,'_'),'[[',1)))]
	infiles			<- merge(infiles, tmp, by='PTY_RUN', all.x=1, allow.cartesian=TRUE)		
	setkey(infiles, PTY_RUN)
	infiles			<- subset(unique(infiles), is.na(OUTFILE))	
	infiles		<- subset(infiles, PTY_RUN%in%c(1,5,7,23,33,34,44,48))
	
	pty.fa		<- do.call('rbind',lapply( seq_len(nrow(infiles)), function(i)
					{
						#i<- 1							
						load(file.path(indir,infiles[i,FILE]))	#loads "pty.seq.rw" "pty.seq"    "seqd"
						setkey(seqd, W_FROM)
						#	select
						if(!is.na(pty.args[['min.ureads.individual']]))	#	select individuals with x unique reads in each window
						{
							seqd	<- subset(seqd, grepl('REF',FILE_ID) || UNIQUE_N>=pty.args[['min.ureads.individual']])
							seqd[, TOTAL_N:=NULL]
							seqd	<- merge(seqd, seqd[, list(TOTAL_N=length(READ)), by='FILE'], by='FILE')
						}
						#	write fasta
						pty.fa		<- seqd[, {	
									#FILE<- 'ptyr1_InWindow_1_to_60.fasta'
									#READ	<- tmp$READ[which(tmp$FILE==FILE)]
									z		<- file.path(outdir,gsub('\\.fasta','_dophy\\.fasta',FILE))
									if(!file.exists(z))
									{
										write.dna(pty.seq[[FILE]][READ,], file=z, format='fasta', colsep='', nbcol=-1)
									}										
									list(TAXA_N=length(READ), PTY_RUN=PTY_RUN[1], W_FROM=W_FROM[1], W_TO=W_TO[1])
								}, by='FILE']
						cat('\nusing fasta files for trees',infiles[i,FILE],', min taxa=',pty.fa[,min(TAXA_N)],'median taxa=',pty.fa[,median(TAXA_N)],'max taxa=',pty.fa[,max(TAXA_N)])
						pty.fa
					}))	
	#
	#	get commands to reconstruct tree
	#
	if(pty.args[['bs.n']]>0)	#	bootstrap on one machine version 
	{
		exa.cmd			<- pty.fa[,	list(CMD=cmd.examl.bootstrap.on.one.machine(outdir, sub("\\.fasta$", "_dophy",FILE), bs.from=0, bs.to=pty.args[['bs.n']]-1, bs.n=pty.args[['bs.n']], opt.bootstrap.by="nucleotide", args.examl=pty.args[['args.examl']])), by='FILE']
		exa.cmd[, RUN_ID:=seq_len(nrow(exa.cmd))]
	}
	if(pty.args[['bs.n']]==0)	#	no bootstrap version
	{		
		exa.cmd			<- pty.fa[,{				 					
					cmd			<- cmd.examl.single(outdir, sub("\\.fasta$", "_dophy",FILE), args.examl=pty.args[['args.examl']])					
					list(CMD=cmd)					
				}, by=c('PTY_RUN','FILE')]
		if(!is.na(pty.args[['exa.n.per.run']]))
			exa.cmd[, RUN_ID:= ceiling(seq_len(nrow(exa.cmd))/pty.args[['exa.n.per.run']])]		
		if(is.na(pty.args[['exa.n.per.run']]))
			exa.cmd[, RUN_ID:=PTY_RUN]
		exa.cmd			<- exa.cmd[,	list(CMD=paste(CMD,collapse='\n',sep='\n')), 	by='RUN_ID']
	}
	exa.cmd
}	

#' @export
pty.get.taxa.combinations<- function(pty.taxan, pty.sel.n)
{
	pty.maxn	<- ceiling( (pty.taxan-pty.sel.n) / (pty.sel.n-1) )
	pty.maxn	<- pty.maxn*(pty.sel.n-1)
	pty.idx		<- matrix( rev(seq_len(pty.maxn))+pty.sel.n, nrow=pty.sel.n-1 )
	#	for each col with lowest number x, repeat x-1 times
	tmp			<- lapply( seq_len(ncol(pty.idx)), function(j)
			{
				z	<- seq_len(min(pty.idx[,j])-1)
				z	<- t(matrix(data=z, nrow=length(z), ncol=nrow(pty.idx)+1))
				z[seq_len(nrow(pty.idx))+1,]	<- NA				
				z[seq_len(nrow(pty.idx))+1,]	<- pty.idx[,j]
				t(z)
			})
	pty.idx		<- do.call('rbind',tmp)
	pty.idx		<- rbind(matrix(seq_len(pty.sel.n),nrow=1), pty.idx)
	pty.idx		<- do.call('rbind',lapply( seq_len(nrow(pty.idx)), function(i) data.table(RUN=i, TX_IDX=pty.idx[i,])))
	pty.idx
}



#' @export
pty.various<- function()
{		
	if(1)
	{
		project.dualinfecions.phylotypes.evaluatereads.150119()
	}	
}