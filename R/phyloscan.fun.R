PR.PACKAGE		<- "phyloscan" 
PR.EVAL.FASTA	<- paste('Rscript', system.file(package=PR.PACKAGE, "phyloscan.evaluate.fasta.Rscript"))

#' @export
pty.cmd.evaluate.fasta<- function(indir, outdir=indir, strip.max.len=Inf, select='', min.ureads.individual=NA, verbose=1)
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
pty.cmd<- function(file.bam, file.ref, window.coord, prog=PROG.PTY, prog.raxml='raxmlHPC-AVX', prog.mafft='mafft', merge.threshold=1, min.read.count=2, quality.trim.ends=30, min.internal.quality=2, merge.paired.reads='-P',no.trees='',keep.overhangs='',out.dir='.')
{	
	stopifnot(is.character(file.bam),is.character(file.ref),is.numeric(window.coord), !length(window.coord)%%2)
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
	if(nchar(no.trees))		
		cmd	<- paste(cmd, no.trees)	
	if(nchar(keep.overhangs))
		cmd	<- paste(cmd, keep.overhangs)
	cmd		<- paste(cmd, '--x-raxml',prog.raxml,'--x-mafft',prog.mafft,'\n')
	tmp		<- gsub('_bam.txt','',basename(file.bam))	
	cmd		<- paste(cmd, 'for file in RAxML_bestTree\\.*.tree; do\n\tmv "$file" "${file//RAxML_bestTree\\./',tmp,'_}"\ndone\n',sep='')
	cmd		<- paste(cmd, "for file in AlignedReads*.fasta; do\n\tsed 's/<unknown description>//' \"$file\" > \"$file\".sed\n\tmv \"$file\".sed \"$file\"\ndone\n",sep='')	
	cmd		<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\tmv "$file" "${file//AlignedReads/',tmp,'_}"\ndone\n',sep='')
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
				windows		<- seq(1,by=pty.win,len=ceiling(max(REF_LEN)/pty.win))
				windows		<- as.vector(rbind( windows,windows-1+pty.win ))
				cmd			<- pty.cmd(	file.bam, file.ref, window.coord=windows, 
										prog=pty.args[['prog']], prog.raxml=pty.args[['raxml']], prog.mafft=pty.args[['mafft']], 
										merge.threshold=pty.args[['merge.threshold']], min.read.count=pty.args[['min.read.count']], quality.trim.ends=pty.args[['quality.trim.ends']], min.internal.quality=pty.args[['min.internal.quality']], 
										merge.paired.reads=pty.args[['merge.paired.reads']], no.trees=pty.args[['no.trees']], keep.overhangs=pty.args[['keep.overhangs']],
										out.dir=pty.args[['out.dir']])
				cmd			<- paste(cmd, pty.cmd.evaluate.fasta(pty.args[['out.dir']], strip.max.len=pty.args[['strip.max.len']], select=paste('^ptyr',PTY_RUN,'_In',sep=''), min.ureads.individual=pty.args[['min.ureads.individual']], verbose=1), sep='')
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
pty.evaluate.tree<- function(pty.infile, indir.tr)
{
	load( pty.infile )
	tmp				<- subset(si, select=c(SANGER_ID,PANGEA_ID))
	set(tmp, NULL, 'SANGER_ID', tmp[, gsub('-','_',SANGER_ID)])
	set(tmp, NULL, 'PANGEA_ID', tmp[, gsub('-','_',PANGEA_ID)])
	setnames(tmp, c('SANGER_ID','PANGEA_ID'),c('FILE_ID','TAXA'))
	pty.runs		<- merge(pty.runs, tmp, by='TAXA',all.x=1)
	tmp				<- pty.runs[, which(is.na(FILE_ID))]
	set(pty.runs, tmp, 'FILE_ID', pty.runs[tmp, TAXA])
	#
	#	collect ML tree files
	#
	ptyfiles		<- data.table(FILE=list.files(indir.tr, 'newick$'))
	ptyfiles[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]
	ptyfiles[, W_FROM:= as.numeric(gsub('InWindow_','',regmatches(FILE,regexpr('InWindow_[0-9]+',FILE))))] 
	ptyfiles[, W_TO:= as.numeric(gsub('to_','',regmatches(FILE,regexpr('to_[0-9]+',FILE))))]		
	ptyfiles		<- subset(ptyfiles, W_FROM<9000)
	#
	ptyfiles[, table(PTY_RUN)]
	#	raw trees w/o any attributes
	pty.ph		<- lapply( seq_len(nrow(ptyfiles)), function(i)
			{				
				ph			<- read.tree(file.path(indir.tr,ptyfiles[i, FILE]))
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
	#	rm trees with just one individual
	#
	ptyfiles[, IND_N:= sapply(seq_along(pty.ph), function(i)  length(unique(gsub('_read.*','',pty.ph[[i]]$tip.label)))		)]
	ptyfiles			<- subset(ptyfiles, IND_N>1)
	ptyfiles[, IDX:=seq_len(nrow(ptyfiles))]
	pty.ph			<- lapply( ptyfiles[, FILE], function(x)	pty.ph[[x]]		)
	names(pty.ph)	<- ptyfiles[, FILE]	
	#
	#	determine root for each run: find taxon with largest distance from BAM of select individuals
	#
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
	ptyfiles		<- merge(ptyfiles,pty.root, by=c('FILE','PTY_RUN'))
	#
	#	root and ladderize trees
	#
	pty.ph.cp						<- copy(pty.ph)	
	#	pty.ph	<- copy(pty.ph.cp)
	for(i in seq_along(pty.ph))
	{
		print(i)
		#i<- 44
		#i<- 5
		root			<- subset(ptyfiles, FILE==names(pty.ph)[i])[, ROOT]
		ph				<- pty.ph[[i]]
		tmp										<- which(ph$tip.label==root)
		ph										<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])
		ph$node.label[ph$node.label=='Root']	<- 0
		ph$node.label							<- as.numeric(ph$node.label)			
		ph				<- ladderize(ph)
		phb				<- data.table(IDX=seq_along(ph$tip.label), BAM=ph$tip.label, FILE_ID= gsub('_read.*','',ph$tip.label))
		#	group edges by individual				
		tmp				<- lapply( phb[, unique(FILE_ID)], function(x)	subset(phb, FILE_ID==x)[, BAM]	)
		names(tmp)		<- phb[, unique(FILE_ID)]
		ph				<- groupOTU(ph, tmp, group='INDIVIDUAL')		
		z	<- merge(data.table(FROM=ph$edge[,1],IDX=ph$edge[,2]), phb, by='IDX', all=1)
		z[, GROUP:= attr(ph,'INDIVIDUAL')[1:nrow(ph$edge)]]
		z	<- unique(subset(z, !is.na(FILE_ID), select=c(FILE_ID, GROUP)))
		attr(ph,'INDIVIDUAL')	<- factor(attr(ph,'INDIVIDUAL'), levels=c(0,z[,as.character(GROUP)]), labels=c('not characterized',z[,FILE_ID]))
		#	group edges by FILL
		tmp				<- as.numeric(gsub('ptyr','',regmatches(names(pty.ph)[i], regexpr('ptyr[0-9]+',names(pty.ph)[i]))))
		phb				<- merge(phb, subset(pty.runs, PTY_RUN==tmp), by='FILE_ID')
		set(phb, NULL, 'FILL', phb[, factor(FILL, levels=c(0,1), labels=c('target','filler'))])
		tmp				<- lapply( phb[, unique(FILL)], function(x)	subset(phb, FILL==x)[, BAM]	)
		names(tmp)		<- phb[, unique(FILL)]
		ph				<- groupOTU(ph, tmp, group='TYPE')				
		z	<- merge(data.table(FROM=ph$edge[,1],IDX=ph$edge[,2]), phb, by='IDX', all=1)
		z[, GROUP:= attr(ph,'TYPE')[1:nrow(ph$edge)]]
		z	<- unique(subset(z, !is.na(FILE_ID), select=c(FILL, GROUP)))
		attr(ph,'TYPE')	<- factor(attr(ph,'TYPE'), levels=c(0,z[,as.character(GROUP)]), labels=c('not characterized',z[,as.character(FILL)]))
		pty.ph[[i]]		<- ph
	}		
	#
	#	save trees
	#
	save(pty.ph, ptyfiles, file=gsub('ptyrunsinput','ptyrunstrees',pty.infile) )	
	#	need node heights for plotting
	tmp					<- ptyfiles[,{										
				ph		<- pty.ph[[FILE]]
				tmp		<- node.depth.edgelength(ph)[1:Ntip(ph)]
				list(BAM=ph$tip.label, HEIGHT=tmp)
			}, by='FILE']	
	ptyfiles				<- merge(ptyfiles, tmp[, list(HMX= max(HEIGHT)), by='FILE'], by='FILE')
	#
	#	plot trees
	#
	require(gridExtra)
	require(colorspace)
	setkey(ptyfiles, PTY_RUN, W_FROM)	
	#ptyfiles			<- subset(ptyfiles, PTY_RUN==1)
	for(ptyr in ptyfiles[, unique(PTY_RUN)])
	{
		#ptyr<- 15
		tmp			<- subset(ptyfiles, PTY_RUN==ptyr)
		#	title
		tmp[, TITLE:=paste('run',PTY_RUN,', window [',W_FROM,',',W_TO,']',sep='')]
		setkey(tmp, W_FROM)
		phs			<- lapply(tmp[, FILE], function(x) pty.ph[[x]])
		names(phs)	<- tmp[, TITLE] 
		#	colours	
		tmp2		<- unique(unlist(lapply(seq_along(phs), function(i)	levels(attr(phs[[i]],'INDIVIDUAL'))	)))
		col			<- c('black',rainbow_hcl(length(tmp2)-1, start = 270, end = -30, c=100, l=50))
		names(col)	<- tmp2	
		phps		<- lapply(seq_along(phs), function(i){
					max.node.height	<- tmp[i,][, HMX]
					p				<- ggtree(phs[[i]], aes(color=INDIVIDUAL, linetype=TYPE)) + 
							geom_nodepoint(size=phs[[i]]$node.label/100*3) +
							geom_tiplab(size=1.2,  hjust=-.1) +							 
							scale_color_manual(values=col, guide = FALSE) +											 
							scale_linetype_manual(values=c('target'='solid','filler'='dotted'),guide = FALSE) +
							theme_tree2() +
							theme(legend.position="bottom") + ggplot2::xlim(0, max.node.height*1.3) +
							labs(x='subst/site', title=names(phs)[i])
					p
				})	
		file	<- file.path( indir.tr, tmp[1,gsub('newick','pdf',gsub('_InWindow_[0-9]+_to_[0-9]+','',FILE))] )
		if(0)
		{
			pdf(file=file, w=120, h=40)	#for win=300
			tmp	 	<- matrix(seq(1, 2*ceiling(length(phps)/2)), ncol=ceiling(length(phps)/2), nrow=2)
			grid.newpage()
			pushViewport(viewport(layout = grid.layout(nrow(tmp), ncol(tmp))))
			for(i in seq_along(phps))
				print(phps[[i]], vp = viewport(layout.pos.row=which(tmp==i, arr.ind=TRUE)[1,'row'], layout.pos.col=which(tmp==i, arr.ind=TRUE)[1,'col']))	
			dev.off()			
		}
		if(1)
		{
			tmp		<- seq_len(ceiling(length(phps)/10))		
			for(i in tmp)
			{			
				pdf(file=gsub('pdf',paste(i,'.pdf',sep=''),file), w=20, h=40)		#for win=60
				grid.newpage()
				pushViewport(viewport(layout=grid.layout(2, 5)))
				z	<- intersect(seq.int((i-1)*10+1, i*10), seq_len(length(phps)))
				for(j in z)
					print(phps[[j]], vp = viewport(layout.pos.row=(j+1)%%2+1, layout.pos.col=(j-1)%%5+1))
				dev.off()	
			}
		}
		if(0)	#devel
		{
			tmp	 	<- matrix(seq(1, 2*ceiling(length(phps)/2)), ncol=ceiling(length(phps)/2), nrow=2)
			grid.newpage()
			pushViewport(viewport(layout = grid.layout(nrow(tmp), ncol(tmp))))
			for(i in seq_along(phps))
				print(phps[[i]], vp = viewport(layout.pos.row=which(tmp==i, arr.ind=TRUE)[1,'row'], layout.pos.col=which(tmp==i, arr.ind=TRUE)[1,'col']))	
			dev.off()
			
			ggtree(ph, aes(color=INDIVIDUAL, linetype=TYPE)) + 
					geom_nodepoint(size=ph$node.label/100*3) +
					#geom_text(aes(label=label), size=1.5,  hjust=-.1) +
					geom_tiplab(size=1.2,  hjust=-.1) +
					scale_color_manual(values=col, guide = FALSE) +											 
					scale_linetype_manual(values=c('target'='solid','filler'='dotted'),guide = FALSE) +
					theme_tree2() +
					theme(legend.position="bottom") + ggplot2::xlim(0, 0.3) +
					labs(x='subst/site')
		}
	}
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
		seqd	<- subset(seqd, UNIQUE_N>=min.ureads.individual)
	}
	#	write fasta
	pty.fa	<- seqd[, {	
				#FILE<- 'ptyr1_InWindow_1_to_60.fasta'
				#READ	<- tmp$READ[which(tmp$FILE==FILE)]
				z		<- file.path(outdir, gsub('\\.fasta','_dophy\\.fasta',FILE))
				if(!file.exists(z))
				{
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
	#infiles		<- subset(infiles, PTY_RUN%in%c(15,17,22))
	
	pty.fa		<- do.call('rbind',lapply( seq_len(nrow(infiles)), function(i)
					{
						#i<- 1							
						load(file.path(indir,infiles[i,FILE]))	#loads "pty.seq.rw" "pty.seq"    "seqd"
						setkey(seqd, W_FROM)
						#
						#	select
						#
						if(!is.na(pty.args[['min.ureads.individual']]))	#	select individuals with x unique reads in each window
							seqd	<- subset(seqd, UNIQUE_N>=pty.args[['min.ureads.individual']])
						#	write fasta
						pty.fa		<- seqd[, {	
									#FILE<- 'ptyr1_InWindow_1_to_60.fasta'
									#READ	<- tmp$READ[which(tmp$FILE==FILE)]
									z		<- file.path(outdir,gsub('\\.fasta','_dophy\\.fasta',FILE))
									if(!file.exists(z))
									{
										write.dna(pty.seq[[FILE]][READ,], file=z, format='fasta', colsep='', nbcol=-1)
									}										
									list(FILE=z, TAXA_N=length(READ), PTY_RUN=PTY_RUN[1], W_FROM=W_FROM[1], W_TO=W_TO[1])
								}, by='FILE']
						cat('\nusing fasta files for trees',i,', min taxa=',pty.fa[,min(TAXA_N)],'median taxa=',pty.fa[,median(TAXA_N)],'max taxa=',pty.fa[,max(TAXA_N)])
						pty.fa
					}))	
	#
	#	get commands to reconstruct tree
	#
	if(pty.args[['bs.n']]>0)	#	bootstrap on one machine version 
	{
		exa.cmd			<- pty.fa[,	list(CMD=cmd.examl.bootstrap.on.one.machine(outdir, sub("\\.[^.]*$", "",FILE), bs.from=0, bs.to=pty.args[['bs.n']]-1, bs.n=pty.args[['bs.n']], opt.bootstrap.by="nucleotide", args.examl=pty.args[['args.examl']])), by='FILE']
		exa.cmd[, RUN_ID:=seq_len(nrow(exa.cmd))]
	}
	if(pty.args[['bs.n']]==0)	#	no bootstrap version
	{		
		exa.cmd			<- pty.fa[,{				 					
					cmd			<- cmd.examl.single(outdir, sub("\\.[^.]*$", "",FILE), args.examl=pty.args[['args.examl']])					
					list(CMD=cmd)					
				}, by=c('PTY_RUN','FILE')]
		exa.cmd[, RUN_ID:= ceiling(seq_len(nrow(exa.cmd))/pty.args[['exa.n.per.run']])]
		exa.cmd			<- exa.cmd[,	list(CMD=paste(CMD,collapse='\n',sep='\n')), 	by='RUN_ID']
	}
	exa.cmd
}	

#' @export
pty.pipeline.examl<- function() 
{
	require(big.phylo)
	#
	#	input args
	#	(used function project.dualinfecions.phylotypes.pipeline.fasta.160110 to create all fasta files)
	#
	#indir			<- "/Users/Oliver/duke/2016_PANGEAphylotypes/phylotypes"
	if(0)	#test on Mac
	{
		work.dir		<- file.path(HOME,"ptyruns")
		out.dir			<- file.path(HOME,"phylotypes_160119")
		hpc.load		<- ''
	}
	if(1)	#coinfections on HPC
	{				
		work.dir		<- file.path(HOME,"coinf_ptinput")
		out.dir			<- file.path(HOME,"coinf_ptoutput_150121")		
		hpc.load		<- "module load intel-suite/2015.1 mpi R/3.2.0"		
	}
	#	get alignment rda files
	if(0)
	{
		infiles			<- data.table(FILE=list.files(out.dir, pattern='fasta$'))
		infiles			<- subset(infiles, !grepl('*',FILE,fixed=1) & !grepl('dophy\\.fasta',FILE))				
		infiles[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]
		invisible(infiles[, {
					cmd			<- pty.cmd.evaluate.fasta(out.dir, strip.max.len=350, select=paste('ptyr',PTY_RUN,'_In',sep=''))
					cat(cmd)
					cmd			<- cmd.hpcwrapper(cmd, hpc.walltime=1, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=1, hpc.load=hpc.load)										
					outfile		<- paste("ptye",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
					cmd.hpccaller(work.dir, outfile, cmd)
					NULL
				}, by='PTY_RUN'])
		stop()
	}
	#	get CMD
	pty.args		<- list(	out.dir=out.dir, work.dir=work.dir, 
								min.ureads.individual=20, min.ureads.candidate=NA, 
								args.examl="-f d -D -m GAMMA", bs.n=0, exa.n.per.run=10)	
	exa.cmd			<- pty.cmdwrap.examl(pty.args)
	cat( exa.cmd[1, cat(CMD)] )		
	stop()
	invisible(exa.cmd[,	{		
					cmd			<- cmd.hpcwrapper(CMD, hpc.walltime=20, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=1, hpc.load=hpc.load)					
					#cmd		<- cmd.hpcwrapper(cmd, hpc.walltime=10, hpc.q="pqeph", hpc.mem="1800mb",  hpc.nproc=1, hpc.load=hpc.load)
					#cat(cmd)
					outfile		<- paste("pexa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
					cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
				}, by='RUN_ID'])		
}

#' @export
pty.pipeline.fasta<- function() 
{
	require(big.phylo)
	#
	#	input args
	#	(used function project.dualinfecions.phylotypes.setup.ZA.160110 to select bam files for one run)
	#
	if(0)	#trm pairs on Mac
	{
		load( file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_ptyrunsinput.rda") )	
		pty.data.dir		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'
		work.dir			<- '/Users/Oliver/duke/2016_PANGEAphylotypes/ptyruns'
		out.dir				<- file.path(HOME,"phylotypes")
		pty.prog			<- '/Users/Oliver/git/phylotypes/phylotypes.py'
		raxml				<- 'raxmlHPC-AVX'
		no.trees			<- '-T'
		
	}	
	if(0)	#trm pairs on HPC
	{
		load( file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_ptyrunsinput.rda") )
		pty.data.dir	<- '/work/or105/PANGEA_mapout/data'
		work.dir		<- file.path(HOME,"ptyruns")
		out.dir			<- file.path(HOME,"phylotypes")
		pty.prog		<- '/work/or105/libs/phylotypes/phylotypes.py'
		raxml			<- 'raxml'
		no.trees		<- '-T'
		hpc.load		<- "module load intel-suite/2015.1 mpi R/3.2.0 raxml/8.2.4 mafft/7 anaconda/2.3.0 samtools"
	}
	if(1)	#coinfections on HPC
	{
		load( file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_coinfrunsinput.rda") )
		pty.data.dir	<- '/work/or105/PANGEA_mapout/data'
		work.dir		<- file.path(HOME,"coinf_ptinput")
		out.dir			<- file.path(HOME,"coinf_ptoutput_150121")
		pty.prog		<- '/work/or105/libs/phylotypes/phylotypes.py'
		raxml			<- 'raxml'
		no.trees		<- '-T'
		hpc.load		<- "module load intel-suite/2015.1 mpi R/3.2.0 raxml/8.2.4 mafft/7 anaconda/2.3.0 samtools"
	}
	#
	#	set up all temporary files and create bash commands
	#
	#	run 160115	window length 300
	if(0)
	{
		pty.args			<- list(	prog=pty.prog, mafft='mafft', raxml=raxml, data.dir=pty.data.dir, work.dir=work.dir, out.dir=out.dir,
										merge.threshold=1, min.read.count=2, quality.trim.ends=30, min.internal.quality=2, merge.paired.reads='-P',no.trees=no.trees, win=300, keep.overhangs='',
										strip.max.len=350, min.ureads.individual=NA)		
		pty.c				<- pty.cmdwrap.fasta(pty.runs, pty.args)		
	}
	#	run 160118	window length 60 & Q1 18 & keep overhangs
	if(1)
	{		
		pty.args			<- list(	prog=pty.prog, mafft='mafft', raxml=raxml, data.dir=pty.data.dir, work.dir=work.dir, out.dir=out.dir,
										merge.threshold=1, min.read.count=2, quality.trim.ends=18, min.internal.quality=2, merge.paired.reads='-P',no.trees=no.trees, win=60, keep.overhangs='--keep-overhangs',
										strip.max.len=350, min.ureads.individual=20)
		pty.c				<- pty.cmdwrap.fasta(pty.runs, pty.args)			
	}
	#	run 160119	window length 60 & Q1 18 & keep overhangs & merge.threshold=3
	if(0)
	{		
		pty.args			<- list(	prog=pty.prog, mafft='mafft', raxml=raxml, data.dir=pty.data.dir, work.dir=work.dir, out.dir=out.dir,
										merge.threshold=3, min.read.count=2, quality.trim.ends=18, min.internal.quality=2, merge.paired.reads='-P',no.trees=no.trees, win=60, keep.overhangs='--keep-overhangs',
										strip.max.len=350, min.ureads.individual=NA)
		pty.c				<- pty.cmdwrap.fasta(pty.runs, pty.args)
		pty.c				<- subset(pty.c, PTY_RUN%in%c(24,26,2,31,34,36,38,44,46,60,69,77,70,78,8))
	}
	if(no.trees=='-T')
	{
		invisible(pty.c[,	{					
					#cmd			<- cmd.hpcwrapper(CMD, hpc.walltime=1, hpc.q="pqeelab", hpc.mem="4800mb",  hpc.nproc=1, hpc.load=hpc.load)
					cmd			<- cmd.hpcwrapper(CMD, hpc.walltime=4, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=1, hpc.load=hpc.load)
					#cat(cmd)					
					outfile		<- paste("pty",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
					cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
				}, by='PTY_RUN'])
	}
	if(0)
	{
		#
		#	add HPC header and submit
		#
		invisible(pty.c[,	{
							#cmd		<- cmd.hpcwrapper(CMD, hpc.walltime=5, hpc.q="pqeelab", hpc.mem="5000mb",  hpc.nproc=1, hpc.load=hpc.load)
							cmd			<- cmd.hpcwrapper(CMD, hpc.walltime=1, hpc.q="pqeph", hpc.mem="1800mb",  hpc.nproc=2, hpc.load=hpc.load)
							outfile		<- paste("pty",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
							stop()
						}, by='PTY_RUN'])
	}	
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