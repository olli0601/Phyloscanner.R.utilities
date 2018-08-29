project.dual<- function()
{	
	CODE.HOME	<<- "/work/or105/libs/Phyloscanner.R.utilities"
	HOME		<<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA'
	HOME		<<- '/work/or105/MUNICH'
	#HOME		<<- "~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA"
	#project.readlength.count.bam.150218()
	#project.readlength.count.bam.171208()
	#project.dual.distances.231015()
	#project.dual.examl.231015()
	#pty.pipeline.fasta()
	#pty.pipeline.phyloscanner.160825()	
	#pty.pipeline.phyloscanner.160915.couples()
	#pty.pipeline.phyloscanner.160915.couples.rerun()
	#pty.pipeline.phyloscanner.170301.firstbatchofall()	
	#pty.pipeline.phyloscanner.170301.firstbatchofall.rerun()	
	#pty.pipeline.phyloscanner.170301.secondstage() 
	#pty.pipeline.phyloscanner.170301.secondstage.ptyr1()		
	#pty.pipeline.phyloscanner.170301.firstbatchsecondbatchofall.fix()	
	#pty.pipeline.phyloscanner.170301.secondstage.ptyrtrees() 	
	#pty.pipeline.phyloscanner.170301.secondstage.rerun()
	#project.bam.read.distribution.calculate()
	#pty.pipeline.phyloscanner.170301.thirdstage() 
	#project.Rakai.ExaMLTree.170601()		
	#pty.pipeline.phyloscanner.170301.secondbatchofall()	
	#project.Rakai.FastTree.170704()
	#project.RakaiAll.setup.RAxMLmodel.170301()
	#pty.pipeline.compress.phyloscanner.output()
	#pty.pipeline.examl()	
	#pty.pipeline.coinfection.statistics()
	#project.dualinfecions.phylotypes.evaluatereads.150119()
	#project.RakaiAll.setup.phyloscanner.170301.stagetwo.preprocessing2()
	#pty.pipeline.phyloscanner.180302.beehive67.process()
	pty.pipeline.phyloscanner.180605.MunichCluster.process()
	#	various   
	if(0) 
	{
		require(big.phylo)
		cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=20, hpc.q="pqeelab", hpc.mem="11gb",  hpc.nproc=1, hpc.load="module load intel-suite/2015.1 mpi R/3.3.3 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools")
		#cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=10, hpc.q=NA, hpc.mem="96gb",  hpc.nproc=1, hpc.load="module load intel-suite/2015.1 mpi R/3.3.3 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools")
		tmp		<- paste('Rscript ',file.path(CODE.HOME, "misc/phyloscan.startme.Rscript"), ' -exe=VARIOUS', '\n', sep='')
		cmd		<- paste(cmd,tmp,sep='\n')
		cat(cmd)	 								
		outfile		<- paste("pv",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(file.path(HOME,"ptyruns"), outfile, cmd)
		quit("no")	
	}			 
} 

pty.various	<- function()
{
	#project.scan.superinfections()
	#project.scan.contaminants()
	#project.readlength.count.all()
	#project.readlength.count.bam.150218()
	#project.readlength.count.bam.171018()
	#project.check.bam.integrity()
	#RakaiFull.divergent.clades.calculate.170811()
	#project.readlength.count.bam.171208()
	#project.readlength.calculate.coverage.171208()
	#project.RakaiAll.setup.phyloscanner.170301.stagetwo.preprocessing()
	project.RakaiAll.setup.phyloscanner.170301.stagetwo.preprocessing2()
}

project.dual.alignments.missing<- function()
{
	infiles	<- data.table(FILE=list.files(out.dir, pattern='alignments.rda'))
	infiles[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]
	setdiff(1:52, infiles[, PTY_RUN])
}

project.dual.alignments.reference<- function()
{
	require(ape)
	
	file			<- '~/Dropbox (SPH Imperial College)/PANGEA_data/HIV1_COM_2012_genome_DNA.fasta'
	outfile			<- '~/git/Phyloscanner.R.utilities/inst/HIV1_compendium_C_B_CPX.fasta'
	ref				<- read.dna(file, format='fasta')
	df				<- data.table(TAXA=rownames(ref))
	df[,SUBTYPE:= toupper(gsub('^[0-9]+_','',regmatches(TAXA,regexpr('^[^\\.]+',TAXA))))]
	df				<- subset(df, grepl('AF460972',TAXA) | grepl('HXB2',TAXA) | SUBTYPE%in%c('C','BC','CD'))
	df[,GENBANK:= regmatches(TAXA,regexpr('[^\\.]+$',TAXA))]
	df[,TAXA_NEW:= df[,paste('R0_REF_',SUBTYPE,'_',GENBANK,'_read_1_count_0',sep='')]]
	ref				<- ref[ df[, TAXA], ]
	rownames(ref)	<- df[,TAXA_NEW]
	write.dna(ref, file=outfile, format='fasta', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/PANGEA_data/HIV1_COM_2012_genome_DNA.fasta'
	outfile			<- '~/git/Phyloscanner.R.utilities/inst/HIV1_compendium_AD_B_CPX.fasta'
	ref				<- read.dna(file, format='fasta')
	df				<- data.table(TAXA=rownames(ref))
	df[,SUBTYPE:= toupper(gsub('^[0-9]+_','',regmatches(TAXA,regexpr('^[^\\.]+',TAXA))))]
	df				<- subset(df, grepl('AF460972|GQ477441|JX140646|AF193276|AF193253',TAXA) | grepl('HXB2',TAXA) | grepl('^A1|^A2|^AD',SUBTYPE) | SUBTYPE=='D' )
	df[,GENBANK:= regmatches(TAXA,regexpr('[^\\.]+$',TAXA))]
	df[,TAXA_NEW:= df[,paste('REF_',SUBTYPE,'_',GENBANK,'_read_1_count_0',sep='')]]
	ref				<- ref[ df[, TAXA], ]
	rownames(ref)	<- df[,TAXA_NEW]
	write.dna(ref, file=outfile, format='fasta', colsep='', nbcol=-1)
	
	file			<- '~/Dropbox (SPH Imperial College)/PANGEA_data/HIV1_COM_2012_genome_DNA.fasta'
	outfile			<- '~/git/Phyloscanner.R.utilities/inst/HIV1_compendium_B.fasta'
	ref				<- read.dna(file, format='fasta')
	df				<- data.table(TAXA=rownames(ref))
	df[,SUBTYPE:= toupper(gsub('^[0-9]+_','',regmatches(TAXA,regexpr('^[^\\.]+',TAXA))))]
	df				<- subset(df, SUBTYPE=='B')
	df[,GENBANK:= regmatches(TAXA,regexpr('[^\\.]+$',TAXA))]
	df[,TAXA_NEW:= df[,paste('REF_',SUBTYPE,'_',GENBANK,sep='')]]
	ref				<- ref[ df[, TAXA], ]
	rownames(ref)	<- df[,TAXA_NEW]
	write.dna(ref, file=outfile, format='fasta', colsep='', nbcol=-1)
}

project.dual.alignments.160110<- function()
{
	outdir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110'
	#	read info
	#file	<- "~/Dropbox (SPH Imperial College)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	#si		<- as.data.table(read.csv(file, stringsAsFactors=FALSE))
	#setnames(si, colnames(si), toupper(gsub('.','_',colnames(si),fixed=1))) 
	#set(si, NULL, 'PANGEA_ID', si[, gsub(' ','',PANGEA_ID)])
	#setnames(si, 'CLINICAL_GENOME_COVERAGE', 'COV')
	
	#	read global PANGEA alignment w SA seqs and split by site	
	file			<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/AfricaCentreSeqs/GlobalAln_PlusAllPANGEA.fasta'
	sq				<- read.dna(file, format='fasta')
	sqi				<- data.table(TAXA=rownames(sq), DUMMY=seq_len(nrow(sq)))
	tmp				<- sqi[, which(duplicated(TAXA))]
	set(sqi, tmp, 'TAXA', sqi[tmp, paste(TAXA,'-R2',sep='')])
	set(sqi, NULL, 'TAXA', sqi[, gsub('_BaseFreqs','',TAXA)])
	set(sqi, NULL, 'TAXA', sqi[, gsub('-','_',TAXA)])
	setkey(sqi, DUMMY)
	rownames(sq)	<- sqi[,TAXA]	
	tmp				<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
	sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]	
	sqi[, PNG:= sqi[, factor(grepl('^PG|^R[0-9]',TAXA),levels=c(TRUE,FALSE),labels=c('Y','N'))]]		
	sqi[, SITE:= NA_character_]
	tmp				<- sqi[, which(PNG=='Y' & grepl('^PG',TAXA))]
	set(sqi, tmp, 'SITE', sqi[tmp, substring(sapply(strsplit(TAXA,'_'),'[[',2),1,2)])
	tmp				<- sqi[, which(PNG=='Y' & grepl('^R[0-9]',TAXA))]
	set(sqi, tmp, 'SITE','ZA')
	sqi[, SEQLOC:= 'LosAlamos']
	set(sqi, sqi[, which(grepl('^PG',TAXA))], 'SEQLOC','Sanger')
	set(sqi, sqi[, which(grepl('^R[0-9]',TAXA))], 'SEQLOC','AfricaCentre')	
	sqi[, PANGEA_ID:= gsub('_R[0-9]+$','',TAXA)]
	save(sqi, sq, file=paste(outdir, '/PANGEA_HIV_n5003_Imperial_v160110_GlobalAlignment.rda',sep=''))
	
	sqi				<- subset(sqi, COV>0)
	seq				<- sq[ subset(sqi, SITE=='UG' | PNG=='N')[, TAXA], ]	
	write.dna( seq, file=paste(outdir,'/PANGEA_HIV_n5003_Imperial_v160110_UG.fasta',sep=''), format='fasta', colsep='', nbcol=-1)	
	save( seq, file=paste(outdir,'/PANGEA_HIV_n5003_Imperial_v160110_UG.R',sep=''))
	seq		<- sq[ subset(sqi, SITE=='BW' | PNG=='N')[, TAXA], ]
	write.dna( seq, file=paste(outdir,'/PANGEA_HIV_n5003_Imperial_v160110_BW.fasta',sep=''), format='fasta', colsep='', nbcol=-1)	
	save( seq, file=paste(outdir,'/PANGEA_HIV_n5003_Imperial_v160110_BW.R',sep=''))	
	seq		<- sq[ subset(sqi, SITE=='ZA' | PNG=='N')[, TAXA], ]
	write.dna( seq, file=paste(outdir,'/PANGEA_HIV_n5003_Imperial_v160110_ZA.fasta',sep=''), format='fasta', colsep='', nbcol=-1)	
	save( seq, file=paste(outdir,'/PANGEA_HIV_n5003_Imperial_v160110_ZA.R',sep=''))
}

project.reference.trees<- function()
{
	require(phytools)
	require(phangorn)
	require(data.table)
	require(ape)
	require(ggtree)
	
	indir	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/reference_trees'
	phd		<- data.table(F=list.files(indir, pattern='RAxML_bestTree', full.names=TRUE))
	phd[, W_FROM:= as.integer(gsub('.*bestTree\\.([0-9]+)_.*','\\1',F))]
	phd[, W_TO:= as.integer(gsub('.*bestTree\\.[0-9]+_([0-9]+)\\.tree','\\1',F))]
	phd[, RUN:=1]
	set(phd, phd[, which(grepl('TAKE2',F))],'RUN',2)
	setkey(phd, W_FROM)
	phd[, IDX:= seq_len(nrow(phd))]
	
	#	read trees
	phs		<- lapply(seq_len(nrow(phd)), function(i)
			{
				ph	<- read.tree(phd[i,F])
			})
	#	re-root half way at MRCA of subtype A taxa	
	ph		<- phs[[1]]
	phi		<- data.table(TAXA=ph$tip.label)
	phi[, ST:= gsub('^([^\\.]+)\\..*','\\1',TAXA)]
	phi[, ST2:= gsub('^[0-9]+_','',ST)]
	root.taxa	<- subset(phi, ST2%in%c('A1','A2'))[, TAXA]
	phs		<- lapply(seq_along(phs), function(i)
			{
				#cat('\n',i)
				ph			<- phs[[i]]
				root.parent	<- getMRCA(ph, root.taxa)				
				tryCatch({
							ans			<- phytools::reroot(ph, root.parent, position=ph$edge.length[which(ph$edge[, 2]==root.parent)] / 2)			
						}, error = function(e)
						{
							cat('\ntree',i,'use child of root.parent to root')
							root.parent	<- max(Children(ph, root.parent))
							ans			<<- phytools::reroot(ph, root.parent, position=ph$edge.length[which(ph$edge[, 2]==root.parent)] / 2)
						})				
				ans				
			})
	names(phs)	<- phd[, paste0(W_FROM,'-',W_TO,'-',RUN)] 
	# 	plot trees
	class(phs) 	<- "multiPhylo"			
	p			<- ggtree(phs, size=0.4) +
		theme_classic() + 
		scale_x_continuous(breaks=seq(0,5,0.1)) +
		labs(x='\nsubst/site') +
		theme(axis.line.x=element_line(), axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
		theme(panel.grid.major.x=element_line(color="grey80", size=.3), panel.grid.major.y=element_blank()) +
		facet_wrap(~.id, ncol=5)
	pdf(file=file.path(indir, 'groupM_reference_trees.pdf'), w=40, h=100*20)
	print(p)
	dev.off()	
	#
	#	there are crazy outliers
	#	remove any tip branches longer than 0.2 or grubbs?
	#	
	max.tip.divergence	<- 0.3
	phs					<- unclass(phs)
	phs.mtd				<- lapply(phs, function(ph)
			{	
				#ph			<- phs[[1]]
				repeat{
					tmp			<- which( ph$edge[,2] <= Ntip(ph) )
					outliers	<- which( ph$edge.length[tmp]>max.tip.divergence )
					if(length(outliers))
						ph			<- drop.tip(ph, ph$edge[ tmp[outliers[1]], 2 ])
					if(length(outliers)<2)
						break
				}								
				ph				
			})
	# 	plot trees
	class(phs.mtd) 	<- "multiPhylo"		 	
	p			<- ggtree(phs.mtd, size=0.4) +
			theme_classic() + 
			scale_x_continuous(breaks=seq(0,5,0.1)) +
			labs(x='\nsubst/site') +
			theme(axis.line.x=element_line(), axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
			theme(panel.grid.major.x=element_line(color="grey80", size=.3), panel.grid.major.y=element_blank()) +
			facet_wrap(~.id, ncol=5)
	pdf(file=file.path(indir, 'groupM_reference_trees_maxtipdivergence30.pdf'), w=40, h=250)
	print(p)
	dev.off()	
	#
	#	still not satisfactory
	#	try Grubbs test, sequentially applied
	#		
	require(outliers)
	grubss.pval			<- 0.01
	phs					<- unclass(phs)
	phs.grb				<- lapply(phs, function(ph)
			{	
				#ph			<- phs[[1]]
				repeat{
					tmp			<- which( ph$edge[,2] <= Ntip(ph) )
					gr 			<- grubbs.test(ph$edge.length[tmp])
					outlier		<- NA
					if(gr$p.value<grubss.pval & grepl('highest',gr$alternative))
						outlier	<- which.max(ph$edge.length[tmp])
					if(!is.na(outlier))					
						ph			<- drop.tip(ph, ph$edge[ tmp[outlier], 2 ])
					if(is.na(outlier))
						break
				}								
				ph				
			})
	# 	plot trees
	class(phs.grb) 	<- "multiPhylo"		 	
	p			<- ggtree(phs.grb, size=0.4) +
			theme_classic() + 
			scale_x_continuous(breaks=seq(0,5,0.1)) +
			labs(x='\nsubst/site') +
			theme(axis.line.x=element_line(), axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
			theme(panel.grid.major.x=element_line(color="grey80", size=.3), panel.grid.major.y=element_blank()) +
			facet_wrap(~.id, ncol=5)
	pdf(file=file.path(indir, paste0('groupM_reference_trees_grubbs',grubss.pval*100,'.pdf')), w=40, h=250*15)
	print(p)
	dev.off()		
	#
	#	keep 5860-2, 8460-2, 9060-1, 9140-2
	#	
	phd[, USE:=1]
	set(phd, phd[, which(RUN==1 & W_FROM%in%c(5860, 8460, 9140))], 'USE', 0)
	set(phd, phd[, which(RUN==2 & W_FROM%in%c(9060))], 'USE', 0)
	
	#	calculate stats on trees
	tmp		<- lapply(seq_along(phs.grb), function(i){
				#i	<- 1
				ph	<- phs.grb[[i]]
				tmp	<- cophenetic.phylo(ph)
				tmp[upper.tri(tmp, diag=TRUE)]	<- NA
				tmp	<- as.numeric(na.omit(as.vector(tmp)))
				data.table(IDX=i, MEAN_PWD=mean(tmp), MEDIAN_PWD=median(tmp), MAX_PWD=max(tmp), SUM_BRL=sum(ph$edge.length))				
			})
	tmp		<- do.call('rbind',tmp)
	phd		<- merge(phd, tmp, by='IDX')
	#	plot stats
	phdp	<- melt(phd, measure.vars=c('MEAN_PWD','MEDIAN_PWD','MAX_PWD','SUM_BRL'))	
	ggplot(phdp, aes(x=W_FROM, y=value)) + 
			geom_point(size=1, colour='grey50') +
			facet_grid(variable~., scales='free_y') +
			theme_bw() +
			scale_x_continuous(breaks=seq(0,10e3,250)) +
			labs(x='\nwindow start', y='') #+ 			
			#geom_smooth(method='loess', span=0.05, se=FALSE, colour='black', size=0.5)
	ggsave(file=file.path(indir, paste0('groupM_reference_trees_grubbs',grubss.pval*100,'_stats.pdf')), w=15, h=10)
	#	save to file
	write.csv(subset(phd, select=c(W_FROM, W_TO, MEAN_PWD, MEDIAN_PWD, MAX_PWD, SUM_BRL)), file=file.path(indir, paste0('groupM_reference_trees_grubbs',grubss.pval*100,'_stats.csv')), row.names=FALSE)	
}

project.dual.alignments.151023<- function()
{
	outdir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/alignments_151113'
	#	read info
	file	<- "~/Dropbox (SPH Imperial College)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	si		<- as.data.table(read.csv(file, stringsAsFactors=FALSE))
	setnames(si, colnames(si), toupper(gsub('.','_',colnames(si),fixed=1))) 
	set(si, NULL, 'PANGEA_ID', si[, gsub(' ','',PANGEA_ID)])
	setnames(si, 'CLINICAL_GENOME_COVERAGE', 'COV')
	
	#	read global PANGEA alignment and split by site	
	file			<- "~/Dropbox (SPH Imperial College)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_GlobalAlignment.fasta"
	file			<- '~/Dropbox (SPH Imperial College)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.fasta'
	sq				<- read.dna(file, format='fasta')
	sqi				<- data.table(TAXA=rownames(sq), DUMMY=seq_len(nrow(sq)))
	tmp				<- sqi[, which(duplicated(TAXA))]
	set(sqi, tmp, 'TAXA', sqi[tmp, paste(TAXA,'-R2',sep='')])
	setkey(sqi, DUMMY)
	rownames(sq)	<- sqi[,TAXA]	
	tmp				<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
	sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]	
	sqi[, PNG:= sqi[, factor(grepl('PG',TAXA),levels=c(TRUE,FALSE),labels=c('Y','N'))]]		
	sqi[, SITE:= NA_character_]
	tmp				<- sqi[, which(PNG=='Y')]
	set(sqi, tmp, 'SITE', sqi[tmp, substring(sapply(strsplit(TAXA,'-'),'[[',2),1,2)])
	sqi[, PANGEA_ID:= gsub('-R[0-9]+','',TAXA)]
	save(sqi, sq, si, file=paste(outdir, '/', gsub('fasta','rda',basename(file)),sep=''))
	
	sqi				<- subset(sqi, COV>0)
	sqi				<- merge(sqi, si, by=c('PANGEA_ID','COV'), all.x=1)
	seq				<- sq[ subset(sqi, SITE=='UG' | PNG=='N')[, TAXA], ]	
	write.dna( seq, file=paste(outdir,'/PANGEAconsensuses_2015-09_Imperial_UG_151113.fasta',sep=''), format='fasta', colsep='', nbcol=-1)	
	save( seq, file=paste(outdir,'/PANGEAconsensuses_2015-09_Imperial_UG_151113.R',sep=''))
	seq		<- sq[ subset(sqi, SITE=='BW' | PNG=='N')[, TAXA], ]
	write.dna( seq, file=paste(outdir,'/PANGEAconsensuses_2015-09_Imperial_BW_151113.fasta',sep=''), format='fasta', colsep='', nbcol=-1)	
	save( seq, file=paste(outdir,'/PANGEAconsensuses_2015-09_Imperial_BW_151113.R',sep=''))
	
	
	#	read contig alignment and split by site
	file	<- "~/Dropbox (SPH Imperial College)/pangea_data/PANGEAcontigs_2015-09_Imperial/contigs_cnsalign_PNGIDn3366_CNTGSn6120_stripped99.fasta"
	cr		<- read.dna(file, format='fasta')
	cri		<- data.table(TAXA=rownames(cr))
	cri[, PNG:= cri[, factor(grepl('^[0-9]+_[0-9]+_[0-9]+.*',gsub('^\\.', '', TAXA)),levels=c(TRUE,FALSE),labels=c('Y','N'))]]
	cri[, SANGER_ID:=NA_character_]
	tmp		<- cri[, which(PNG=='Y')]
	set(cri, tmp, 'SANGER_ID', cri[tmp, sapply(strsplit( gsub('^\\.', '', TAXA), '.', fixed=1),'[[',1)] )	
	tmp		<- subset(cri, PNG=='Y')[, list(TAXA=TAXA, CNTG_ID_NEW=seq_along(TAXA), CONTG_ID= gsub(SANGER_ID,'',gsub('^\\.', '', TAXA))), by='SANGER_ID']
	cri		<- merge(cri, tmp, all.x=1, by=c('SANGER_ID','TAXA'))
	setnames(sqi, 'TAXA', 'PANGEA_ID_WDUP')
	cri[, PNG:=NULL]
	cri		<- merge(cri, subset(sqi, !is.na(SANGER_ID)), by='SANGER_ID', all.x=1)
	stopifnot( nrow(subset(cri, is.na(SANGER_ID) & PNG=='Y'))==0 )
	cri[, SITE:= NA_character_]
	tmp		<- cri[, which(PNG=='Y')]
	set(cri, tmp, 'SITE', cri[tmp, substring(sapply(strsplit(PANGEA_ID,'-'),'[[',2),1,2)])		
	tmp		<- cri[, list(TAXA_NEW= ifelse( is.na(CNTG_ID_NEW), paste('Ref.',TAXA,sep=''), paste(PANGEA_ID_WDUP,'-C',CNTG_ID_NEW,sep='') )), by='TAXA']
	cri		<- merge(cri, tmp, by='TAXA')
	setkey(cri, TAXA)
	rownames(cr)	<- cri[rownames(cr),][, TAXA_NEW]
	seq		<- cr[ subset(cri, SITE=='UG' | PNG=='N')[, TAXA_NEW], ]
	write.dna( seq, file=paste(outdir,'/PANGEAcontigs_2015-09_Imperial_UG_151023.fasta',sep=''), format='fasta', colsep='', nbcol=-1)
	save( seq, file=paste(outdir,'/PANGEAcontigs_2015-09_Imperial_UG_151023.R',sep=''))
	seq		<- cr[ subset(cri, SITE=='BW' | PNG=='N')[, TAXA_NEW], ]
	write.dna( seq, file=paste(outdir,'/PANGEAcontigs_2015-09_Imperial_BW_151023.fasta',sep=''), format='fasta', colsep='', nbcol=-1)	
	save( seq, file=paste(outdir,'/PANGEAcontigs_2015-09_Imperial_BW_151023.R',sep=''))
	
	#	save info on consensus and contigs
	save( sqi, cri, file=paste(outdir,'/PANGEAinfo_2015-09_Imperial.R',sep=''))
	
	#	next: distances
}


project.dualinfecions.dev<- function()
{
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/coinf_ptoutput_UG60'
	infiles	<- data.table(FILE=list.files(indir, pattern='fasta$'))
	infiles[, FROM:= as.numeric(gsub('NoRefAlignedReadsInWindow_','',regmatches(FILE, regexpr('NoRefAlignedReadsInWindow_[0-9]+', FILE))))]
	
	tmp		<- infiles[, {
				s	<- read.dna(file.path(indir,FILE),format='fa')
				list(TAXA_N=nrow(s), COV=ncol(s))
			}, by='FILE']
	infiles	<- merge(tmp, infiles, by='FILE')
	ggplot(infiles, aes(x=FROM, y=TAXA_N)) + geom_step()
	ggplot(infiles, aes(x=FROM, y=COV)) + geom_step()
}


project.dualinfecions.copy.bam<- function()
{
	tmp		<- subset(si, grepl('^PG14-ZA', PANGEA_ID), c(PANGEA_ID, SANGER_ID))[, SANGER_ID]
	tmp		<- subset(pty.clu, !is.na(SANGER_ID))[, SANGER_ID]
	paste('zip SangerAC.zip ',paste(tmp,'_ref* ',collapse='',sep=''),paste(tmp,'.bam ',collapse='',sep=''),sep='')
}


project.dualinfecions.phylotypes.setup.trmpairs.ZA.160110<- function()
{
	#
	#	input args
	#
	pty.gd		<- 0.2
	pty.sel.n	<- 20
	
	infile		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/PANGEA_HIV_n5003_Imperial_v160110_GlobalAlignment.rda'
	load(infile)	#loads sqi, sq
	indir		<- "~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/data"
	outdir		<- indir
	infile		<-  "PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500.rda"
	load(paste(indir,infile,sep='/'))	#loads "ph" "dist.brl" "ph.gdtr"  "ph.mrca"
	#	add SANGER_ID
	infile.s	<- "~/Dropbox (SPH Imperial College)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	si			<- as.data.table(read.csv(infile.s, stringsAsFactors=FALSE))
	setnames(si, colnames(si), toupper(gsub('.','_',colnames(si),fixed=1))) 
	set(si, NULL, 'PANGEA_ID', si[, gsub(' ','',PANGEA_ID)])
	setnames(si, 'CLINICAL_GENOME_COVERAGE', 'COV')
	tmp			<- subset(si, grepl('^PG14-ZA', PANGEA_ID), c(PANGEA_ID, SANGER_ID))
	set(tmp, NULL, 'PANGEA_ID', tmp[,gsub('-','_',PANGEA_ID)])
	sqi			<- merge(sqi, tmp, by='PANGEA_ID', all.x=1)	
	#	delete duplicates identified by Tulio from tree
	infile.dup	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/PANGEA_HIV_n5003_Imperial_v160110_ZA_Duplicates.csv'
	dupl		<- as.data.table(read.csv(infile.dup, stringsAsFactors=FALSE))
	dupl		<- subset(dupl, Duplicated==1, select=c(strains, Duplicated))
	dupl[, DUP_ID:= seq_len(nrow(dupl))]
	set(dupl, NULL, 'strains', dupl[, gsub(' +$','',gsub('^ +','',strains))])	
	dupl		<- dupl[, list(TAXA= strsplit(strains,' ')[[1]]), by='DUP_ID']
	dupl		<- merge(dupl, sqi, by='TAXA', all.x=1)
	setkey(dupl, DUP_ID)
	#write.dna( sq[dupl[, TAXA],], format='fasta', file=gsub('.csv','.fasta',infile.dup))
	dupl		<- dupl[, list(TAXA= TAXA[which.min(COV)]),by='DUP_ID']
	ph$tip.label<- gsub('-','_',ph$tip.label) 
	ph			<- drop.tip(ph, dupl[,TAXA])
	#	need to recalculate stats..
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)	
	ph.gdtr							<- cophenetic.phylo(ph)
	ph.mrca							<- mrca(ph)	
	#
	#	determine large clusters
	#
	clustering	<- hivc.clu.clusterbythresh(ph, thresh.brl=pty.gd, dist.brl=dist.brl, retval="all")	
	pty.clu		<- subset(data.table(TX_IDX=seq_len(Ntip(ph)), TAXA= ph$tip.label, CLU_ID=clustering$clu.mem[seq_len(Ntip(ph))]), !is.na(CLU_ID))	
	#	reduce to clusters containing at least one ZA sequence
	pty.clu		<- merge(pty.clu, sqi,by='TAXA')
	tmp			<- pty.clu[, list(ANY_NOT_INCOUNTRY= any(is.na(SITE) | SITE!='ZA'), ALL_NOT_INCOUNTRY= all(is.na(SITE) | SITE!='ZA')), by='CLU_ID']
	cat('\nInspecting clusters if not in-country')
	print(tmp[, table(ANY_NOT_INCOUNTRY, ALL_NOT_INCOUNTRY)])
	tmp			<- subset(tmp, !ALL_NOT_INCOUNTRY)
	pty.clu		<- merge(pty.clu, subset(tmp, select=c(CLU_ID)), by='CLU_ID')
	#	reduce to clusters of ZA sequences
	pty.clu		<- subset(pty.clu, PNG=='Y')
	pty.clu		<- merge(pty.clu, pty.clu[, list(CLU_N= length(TAXA)), by='CLU_ID'], by='CLU_ID')
	pty.clu		<- subset(pty.clu, CLU_N>1)
	cat('\nFound in-country clusters of size:')	
	print( unique(subset(pty.clu, select=c(CLU_ID, CLU_N)))[, table(CLU_N)] )
	#
	#	get all combinations within each cluster
	#
	pty.fill	<- which(grepl('^R[0-9]+|^PG',ph$tip.label))	
	pty.runs	<- pty.clu[,{
				#CLU_ID	<- 1; TX_IDX	<- subset(pty.clu, CLU_ID==1)[, TX_IDX]
				ans	<- pty.get.taxa.combinations(length(TX_IDX), pty.sel.n)	# get all combinations within cluster
				setnames(ans, 'TX_IDX', 'IDX')
				tx	<- TX_IDX
				tmp	<- pty.sel.n-length(TX_IDX)				
				if(tmp>0)
					tx	<- c(tx, sample(setdiff(pty.fill,TX_IDX), tmp))
				ans	<- merge(ans, data.table(TX_IDX=tx, IDX=seq_along(tx), FILL= c(rep(0,length(TX_IDX)),rep(1,tmp))), by='IDX')
				subset(ans, select=c(RUN, TX_IDX, FILL))
			}, by='CLU_ID']
	setnames(pty.runs, 'RUN', 'RUN_OF_CLU')
	tmp			<- unique(subset(pty.runs, select=c(CLU_ID,RUN_OF_CLU)))
	setkey(tmp, CLU_ID, RUN_OF_CLU)
	tmp[, PTY_RUN:= seq_len(nrow(tmp))]
	pty.runs	<- merge(tmp, pty.runs,by=c('CLU_ID','RUN_OF_CLU'))
	pty.runs[, RUN_OF_CLU:=NULL]
	pty.runs[, TAXA:= ph$tip.label[ TX_IDX] ]
	#
	tmp			<- subset(si, select=c(SANGER_ID, PANGEA_ID))
	set(tmp, NULL, 'PANGEA_ID', tmp[, gsub('-','_',PANGEA_ID)])
	setnames(tmp, c('PANGEA_ID','SANGER_ID'), c('TAXA','FILE_ID'))
	pty.runs	<- merge(pty.runs, tmp, by='TAXA', all.x=1)
	tmp			<- pty.runs[, which(is.na(FILE_ID))]
	set(pty.runs, tmp,'FILE_ID', pty.runs[tmp, TAXA])
	setkey(pty.runs, PTY_RUN)	
	#	
	cat('\nNumber of clusters=', pty.runs[, length(unique(CLU_ID))])
	cat('\nNumber of scheduled phylotype runs=', pty.runs[, max(PTY_RUN)])
	cat('\nNumber of selected taxa=', subset(pty.runs, !FILL)[, length(unique(TAXA))])
	cat('\nLargest cluster size=', unique(subset(pty.clu, select=c(CLU_ID, CLU_N)))[,  max(as.numeric(names(table(CLU_N))))] )
	tmp			<- paste(indir, '/', gsub('\\.rda','_ptyrunsinput\\.rda',infile), sep='')
	save(pty.runs, pty.clu, ph, dist.brl, ph.gdtr, ph.mrca, clustering, sqi, sq, file= tmp)	
}

project.dualinfecions.phylotypes.setup.coinfections.ZA.160110<- function()
{
	#
	#	input args
	#	
	pty.gd		<- 0.2
	pty.sel.n	<- 15
		
	infile		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/PANGEA_HIV_n5003_Imperial_v160110_GlobalAlignment.rda'
	load(infile)	#loads sqi, sq
	indir		<- "~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/data"
	outdir		<- indir
	infile		<-  "PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500.rda"
	load(paste(indir,infile,sep='/'))	#loads "ph" "dist.brl" "ph.gdtr"  "ph.mrca"
	#	add SANGER_ID
	infile.s	<- "~/Dropbox (SPH Imperial College)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	si			<- as.data.table(read.csv(infile.s, stringsAsFactors=FALSE))
	setnames(si, colnames(si), toupper(gsub('.','_',colnames(si),fixed=1))) 
	set(si, NULL, 'PANGEA_ID', si[, gsub(' ','',PANGEA_ID)])
	setnames(si, 'CLINICAL_GENOME_COVERAGE', 'COV')
	tmp			<- subset(si, grepl('^PG14-ZA', PANGEA_ID), c(PANGEA_ID, SANGER_ID))
	set(tmp, NULL, 'PANGEA_ID', tmp[,gsub('-','_',PANGEA_ID)])
	sqi			<- merge(sqi, tmp, by='PANGEA_ID', all.x=1)	
	#	delete duplicates identified by Tulio from tree
	infile.dup	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/PANGEA_HIV_n5003_Imperial_v160110_ZA_Duplicates.csv'
	dupl		<- as.data.table(read.csv(infile.dup, stringsAsFactors=FALSE))
	dupl		<- subset(dupl, Duplicated==1, select=c(strains, Duplicated))
	dupl[, DUP_ID:= seq_len(nrow(dupl))]
	set(dupl, NULL, 'strains', dupl[, gsub(' +$','',gsub('^ +','',strains))])	
	dupl		<- dupl[, list(TAXA= strsplit(strains,' ')[[1]]), by='DUP_ID']
	dupl		<- merge(dupl, sqi, by='TAXA', all.x=1)
	setkey(dupl, DUP_ID)
	#write.dna( sq[dupl[, TAXA],], format='fasta', file=gsub('.csv','.fasta',infile.dup))
	dupl		<- dupl[, list(TAXA= TAXA[which.min(COV)]),by='DUP_ID']
	ph$tip.label<- gsub('-','_',ph$tip.label) 
	ph			<- drop.tip(ph, dupl[,TAXA])
	#	need to recalculate stats..
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)	
	ph.gdtr							<- cophenetic.phylo(ph)
	ph.mrca							<- mrca(ph)	
	#
	#	select closest 15 individuals
	#
	clustering	<- hivc.clu.clusterbythresh(ph, thresh.brl=pty.gd, dist.brl=dist.brl, retval="all")	
	pty.clu		<- subset(data.table(TX_IDX=seq_len(Ntip(ph)), TAXA= ph$tip.label, CLU_ID=clustering$clu.mem[seq_len(Ntip(ph))]), !is.na(CLU_ID))	
	#	reduce to clusters containing at least one ZA sequence
	pty.clu		<- merge(pty.clu, sqi,by='TAXA')
	tmp			<- pty.clu[, list(ANY_NOT_INCOUNTRY= any(is.na(SITE) | SITE!='ZA'), ALL_NOT_INCOUNTRY= all(is.na(SITE) | SITE!='ZA')), by='CLU_ID']
	cat('\nInspecting clusters if not in-country')
	print(tmp[, table(ANY_NOT_INCOUNTRY, ALL_NOT_INCOUNTRY)])
	tmp			<- subset(tmp, !ALL_NOT_INCOUNTRY)
	pty.clu		<- merge(pty.clu, subset(tmp, select=c(CLU_ID)), by='CLU_ID')
	#	reduce to clusters of ZA sequences
	pty.clu		<- subset(pty.clu, PNG=='Y')
	pty.clu		<- merge(pty.clu, pty.clu[, list(CLU_N= length(TAXA)), by='CLU_ID'], by='CLU_ID')
	pty.clu		<- subset(pty.clu, CLU_N>1)
	cat('\nFound in-country clusters of size:')	
	print( unique(subset(pty.clu, select=c(CLU_ID, CLU_N)))[, table(CLU_N)] )
	#	get distance data.table
	ph.gdf		<- as.data.table(melt(ph.gdtr,value.name="GD"))
	setnames(ph.gdf, c('Var1','Var2'),c('TAXA','TAXA2'))
	tmp			<- data.table(TX_IDX=seq_len(Ntip(ph)), TAXA= ph$tip.label)
	tmp			<- merge(tmp, sqi,by='TAXA')
	tmp			<- subset(tmp, SITE=='ZA')
	tmp[, DUMMY:=NULL]
	ph.gdf		<- merge(ph.gdf, subset(tmp, select=TAXA), by='TAXA')
	ph.gdf		<- merge(ph.gdf, data.table(TAXA2=tmp[,TAXA]), by='TAXA2')
	ph.gdf		<- subset(ph.gdf, TAXA!=TAXA2)
	setkey(ph.gdf, TAXA, GD)	
	#	initialize pty.runs
	setkey(pty.clu, CLU_ID)	
	tmp			<- unique(pty.clu)
	tmp			<- tmp[order(-CLU_N),]	
	tx.seeds	<- c(tmp[, TAXA], setdiff(subset(sqi, SITE=='ZA')[,TAXA], tmp[, TAXA]))	
	#	not all seeds in tree (ie no consensus)
	tx.seeds	<- intersect( tx.seeds, ph$tip.label )
	cat('\ninitial seeds n=',length(tx.seeds))
	ptyr		<- 1L	
	pty.runs	<- data.table(PTY_RUN=NA_integer_, TAXA=NA_character_, GD=NA_real_)
	#	fill pty.runs
	while(length(tx.seeds))
	{
		seed		<- tx.seeds[1]
		tmp			<- data.table(TAXA=c(seed, subset(ph.gdf, TAXA==seed )[seq_len(pty.sel.n-1), TAXA2]), GD=c(0,subset(ph.gdf, TAXA==seed )[seq_len(pty.sel.n-1), GD]))
		tmp[, PTY_RUN:=ptyr]
		tx.seeds	<- setdiff( tx.seeds, tmp[,TAXA] )
		ph.gdf		<- subset(ph.gdf, !TAXA%in%tmp[,TAXA] & !TAXA2%in%tmp[,TAXA]) 
		cat('\nrun',ptyr,'select n=',nrow(tmp),'remaining seeds n=',length(tx.seeds), 'remaining distances n=', nrow(ph.gdf))
		pty.runs	<- rbind(tmp, pty.runs)
		ptyr		<- ptyr+1L
		
	}
	pty.runs	<- subset(pty.runs, !is.na(TAXA))
	#pty.runs[, table(PTY_RUN)]
	#	NOTE!!! merge last two since last v small
	set(pty.runs, pty.runs[, which(PTY_RUN==53L)], 'PTY_RUN', 52L)
	pty.runs[, FILL:=0L]
	pty.runs[, TX_IDX:=match(pty.runs[,TAXA],ph$tip.label)]
	pty.runs[, GD:=NULL]
	#
	tmp			<- subset(si, select=c(SANGER_ID, PANGEA_ID))
	set(tmp, NULL, 'PANGEA_ID', tmp[, gsub('-','_',PANGEA_ID)])
	setnames(tmp, c('PANGEA_ID','SANGER_ID'), c('TAXA','FILE_ID'))
	pty.runs	<- merge(pty.runs, tmp, by='TAXA', all.x=1)
	tmp			<- pty.runs[, which(is.na(FILE_ID))]
	set(pty.runs, tmp,'FILE_ID', pty.runs[tmp, TAXA])
	setkey(pty.runs, PTY_RUN)
	pty.runs[, FILL:=0]
	#	
	cat('\nNumber of clusters=', pty.runs[, length(unique(CLU_ID))])
	cat('\nNumber of scheduled phylotype runs=', pty.runs[, max(PTY_RUN)])
	cat('\nNumber of selected taxa=', subset(pty.runs, !FILL)[, length(unique(TAXA))])	
	tmp			<- paste(indir, '/', gsub('\\.rda','_coinfrunsinput\\.rda',infile), sep='')
	save(pty.runs, pty.clu, ph, dist.brl, ph.gdtr, ph.mrca, clustering, sqi, sq, file= tmp)	
}

project.dualinfecions.UG.selectsamplesbycoverage.160219<- function()
{	
	min.coverage	<- 600
	min.depth		<- 500
	infile			<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/readlengths/bam_stats_150218.rda'
	load(infile)
	infile			<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/PANGEA_HIV_n5003_Imperial_v160110_GlobalAlignment.rda'
	load(infile)	#loads sqi, sq
	outdir			<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/data'
	#	add SANGER_ID
	infile.s	<- "~/Dropbox (SPH Imperial College)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	si			<- as.data.table(read.csv(infile.s, stringsAsFactors=FALSE))
	setnames(si, colnames(si), toupper(gsub('.','_',colnames(si),fixed=1))) 
	set(si, NULL, 'PANGEA_ID', si[, gsub(' ','',PANGEA_ID)])
	setnames(si, 'CLINICAL_GENOME_COVERAGE', 'COV')
	tmp			<- subset(si, grepl('^PG14-UG', PANGEA_ID), c(PANGEA_ID, SANGER_ID))
	set(tmp, NULL, 'PANGEA_ID', tmp[,gsub('-','_',PANGEA_ID)])
	#	of the duplicate PANGEA_IDs, consider only those with larger coverage	
	sqi			<- merge(sqi, tmp, by='PANGEA_ID', all.x=1, allow.cartesian=TRUE)
	tmp			<- sqi[, list(SANGER_ID=SANGER_ID[which.max(COV)]), by='PANGEA_ID']
	sqi			<- merge(sqi, tmp, by=c('PANGEA_ID','SANGER_ID'))
	#
	#	consider reads only from Uganda that have coverage at least 600 bp
	#	
	tmp			<- subset(sqi, SITE=='UG' & COV>=min.coverage, c(PANGEA_ID, SANGER_ID, COV))
	setnames(bam.cov, c('FILE_ID','COV'), c('SANGER_ID','DEPTH'))
	bam.cov		<- merge(bam.cov, tmp, by='SANGER_ID')
	#	total individuals from UG: 3628
	#	total individuals from UG with min coverage: 3074
	#
	#	define chunks of the genomes with more than 100, 200, 500, 1000 read depth 
	#
	bam.ch		<- do.call('rbind',lapply(c(100,200,500,1000), function(x)
			{
				bam.ch		<- subset(bam.cov, DEPTH>=x)
				bam.ch[, POS_NEXT:= POS+REP]	
				bam.ch		<- bam.ch[, list(POS=POS, DEPTH=DEPTH, REP=REP, CHUNK=cumsum(as.numeric(c(TRUE, POS[-1]!=POS_NEXT[-length(POS_NEXT)])))), by='SANGER_ID']
				bam.ch		<- bam.ch[, list(POS_CH=min(POS), REP_CH=sum(REP), DEPTH_CH= sum(DEPTH*REP)/sum(REP) ), by=c('SANGER_ID','CHUNK')]
				bam.ch[, DEPTH_MIN:=x]
				bam.ch
			}))
	#
	# select chunks with at least 300 sites and individuals with at least 600 selected sites (can be in same chunk, we just need to be able to define at least two windows)
	#
	bam.chs		<- subset(bam.ch, REP_CH>=300)
	bam.chs		<- merge(bam.chs, bam.chs[, list(REP_IND=sum(REP_CH)), by=c('DEPTH_MIN','SANGER_ID')], by=c('DEPTH_MIN','SANGER_ID'))
	bam.chs		<- subset(bam.chs, REP_IND>min.coverage)
	#
	# count how many individuals by the different read depth thresholds
	#
	bam.chs[, list(IND_SELECTED=length(unique(SANGER_ID)), COV_AVG=mean(unique(REP_IND))), by='DEPTH_MIN']
	#   DEPTH_MIN IND_SELECTED COV_AVG
	#1:       100         2858 4799.172
	#2:       200         2819 4655.781
	#3:       500         2706 4297.860
	#4:      1000         2558 4000.289
	# 	OK, try a minimum depth of 500. We keep 88% of individuals with at least 600bp coverage and 74% of individuals with some coverage 
	bam.chs		<- subset(bam.chs, DEPTH_MIN==min.depth)	
	setkey(bam.chs, POS_CH)
	set(bam.chs, NULL, 'SANGER_ID', bam.chs[, factor(SANGER_ID, levels=unique(SANGER_ID), labels=unique(SANGER_ID))])
	ggplot(bam.chs, aes(x=SANGER_ID, xend=SANGER_ID, y=POS_CH, yend=POS_CH+REP_CH-1L)) +
			scale_y_continuous(expand=c(0,0),limits=c(0,10e3), breaks=seq(0,10e3,500), minor_breaks=seq(0,10e3,100)) +
			geom_segment() + theme_bw() + labs(y='genome position', x='Patient') + coord_flip() +
			theme(axis.text.y=element_text(size=2), axis.ticks.y=element_line(size=0.1))
	ggsave(file=file.path(outdir,'PANGEA_HIV_n5003_Imperial_v160110_UG_selectedchunks.pdf'), w=10, h=60, limitsize = FALSE)
	# 	OK, build two alignments of consensus sequences in two deep reed regions: 1-1700 and 1700-9000
	bam.chs			<- merge(bam.chs, bam.chs[, {
						z	<- which.max(POS_CH)
						list(POS_INDMIN=min(POS_CH), POS_INDMAX=POS_CH[z]+REP_CH[z]-1L)	
					}, by='SANGER_ID'], by='SANGER_ID')
	tmp				<- merge(sqi, unique(subset(bam.chs, POS_INDMIN>1700, SANGER_ID)), by='SANGER_ID')[, TAXA]
	tmp				<- c(tmp, subset(sqi, is.na(SITE))[,TAXA])
	tmp				<- as.character(sq[tmp,2001:ncol(sq)])
	tmp[tmp=='?']	<- '-'
	tmp				<- as.DNAbin(tmp)
	write.dna(tmp , format='fa', file=file.path(outdir, 'PANGEA_HIV_n5003_Imperial_v160110_UG_p15toend.fasta'))
	tmp				<- merge(sqi, unique(subset(bam.chs, POS_INDMIN<=1700, SANGER_ID)), by='SANGER_ID')[, TAXA]
	tmp				<- c(tmp, subset(sqi, is.na(SITE))[,TAXA])
	tmp				<- as.character(sq[tmp,1:2000])
	tmp[tmp=='?']	<- '-'
	tmp				<- as.DNAbin(tmp)	
	write.dna( tmp, format='fa', file=file.path(outdir, 'PANGEA_HIV_n5003_Imperial_v160110_UG_gag.fasta'))
	#	run FastTree2 /Users/Oliver/git/big.phylo/inst/FastTree -nt -gtr -gamma
	#	OK, check which individuals have coverage in which segments
	bam.chs[, POS_INDMIN_C:= cut(POS_INDMIN, breaks=c(-Inf,50,650,2000,3900,Inf), labels=c('1','50','650','2000','3900'))]
	bam.chs[, POS_INDMAX_C:= cut(POS_INDMAX, breaks=c(-Inf,50,650,2000,3900,Inf), labels=c('50','650','2000','3900','9000'))]
	bam.chs[, GROUP:=paste(POS_INDMIN_C,'-',POS_INDMAX_C,sep='')]
	setkey(bam.chs, SANGER_ID)
	bam.inds		<- subset(unique(bam.chs), select=c(SANGER_ID, POS_INDMIN_C, POS_INDMAX_C, GROUP))
	bam.inds[, table(POS_INDMIN_C, POS_INDMAX_C)]
	#				POS_INDMAX_C
	#POS_INDMIN_C   	50  650 2000 3900  9000
	#1       			0    0  614  163   1646
	#50     			0    0   53    5   21
	#650     			0    0   80   32   57
	#2000    			0    0    0    0   1
	#3900    			0    0    0    0   34	
	save(bam.chs, bam.inds, file=file.path(outdir,'PANGEA_HIV_n5003_Imperial_v160110_UG_selectedbycoverage_160219.rda'))
}

project.dualinfecions.UG.mltree.160219<- function()
{	
	#
	#	UG consensus BEST tree on gag
	#
	indir		<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/data'	
	infile		<- "PANGEA_HIV_n5003_Imperial_v160110_UG_gag_fasttree.newick"
	ph			<- read.tree( paste(indir,'/',infile,sep='') )	
	#	reroot at SIV		
	tmp			<- which(grepl('AF103818',ph$tip.label))
	ph			<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])	
	ph 			<- ladderize( ph )
	ph.node.bs						<- ph$node.label
	ph.node.bs[ph.node.bs=='Root']	<- NA
	ph.node.bs						<- as.numeric(ph.node.bs)
	ph.node.bs[is.na(ph.node.bs)]	<- 0
	ph$node.label					<- ph.node.bs
	ph$node.label[ph$node.label>1]	<- 1
	#	plot just the phylogeny
	pdf(file=file.path(indir,gsub('newick','pdf',infile)), w=20, h=120)
	plot(ph, cex=0.2)
	dev.off()
	
	stat.fun						<- hivc.clu.min.transmission.cascade
	#stat.fun						<- max
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)
	thresh.brl						<- 0.07
	thresh.bs						<- 0.75	
	clustering						<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
	table(clustering$size.tips)
	tip.color						<- rep('black',Ntip(ph))
	tip.color[ grepl('PG14',ph$tip.label) ]	<- 'DarkRed'
	file							<- paste(outdir,"/", gsub('\\.newick',paste('_gd',100*thresh.brl,'bs',thresh.bs*100,'\\.pdf',sep=''), infile), sep='')
	invisible(hivc.clu.plot(ph, clustering[["clu.mem"]], cex.edge.incluster=3, tip.color=tip.color, file=file, pdf.scaley=100, show.tip.label=TRUE, pdf.width=30))	
	ph.gdtr							<- cophenetic.phylo(ph)
	ph.mrca							<- mrca(ph)
	save(ph, dist.brl, ph.gdtr, ph.mrca, clustering, file=paste(indir,'/',gsub('\\.newick','\\.rda',infile),sep=''))
	#
	#
	#
	infile		<- "PANGEA_HIV_n5003_Imperial_v160110_UG_p15toend_fasttree.newick"
	ph			<- read.tree( paste(indir,'/',infile,sep='') )	
	#	reroot at SIV		
	tmp			<- which(grepl('AF103818',ph$tip.label))
	ph			<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])	
	ph 			<- ladderize( ph )
	ph.node.bs						<- ph$node.label
	ph.node.bs[ph.node.bs=='Root']	<- NA
	ph.node.bs						<- as.numeric(ph.node.bs)
	ph.node.bs[is.na(ph.node.bs)]	<- 0
	ph$node.label					<- ph.node.bs
	ph$node.label[ph$node.label>1]	<- 1
	#	plot just the phylogeny
	pdf(file=file.path(indir,gsub('newick','pdf',infile)), w=20, h=60)
	plot(ph, cex=1)
	dev.off()
	
	stat.fun						<- hivc.clu.min.transmission.cascade
	#stat.fun						<- max
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)
	thresh.brl						<- 0.14
	thresh.bs						<- 0.75	
	clustering						<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")
	table(clustering$size.tips)
	tip.color						<- rep('black',Ntip(ph))
	tip.color[ grepl('PG14',ph$tip.label) ]	<- 'DarkRed'
	file							<- paste(outdir,"/", gsub('\\.newick',paste('_gd',100*thresh.brl,'bs',thresh.bs*100,'\\.pdf',sep=''), infile), sep='')
	invisible(hivc.clu.plot(ph, clustering[["clu.mem"]], cex.edge.incluster=3, tip.color=tip.color, file=file, pdf.scaley=3, show.tip.label=TRUE, pdf.width=20))	
	ph.gdtr							<- cophenetic.phylo(ph)
	ph.mrca							<- mrca(ph)
	save(ph, dist.brl, ph.gdtr, ph.mrca, clustering, file=paste(indir,'/',gsub('\\.newick','\\.rda',infile),sep=''))
}

project.dualinfecions.UG.mltree.160217<- function()
{	
	#
	#	UG consensus BEST tree
	#
	indir		<- "/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110"
	infile		<- "PANGEA_HIV_n5003_Imperial_v160110_UG_NoQ_fast2.newick"
	ph			<- read.tree( paste(indir,'/',infile,sep='') )	
	#	reroot at SIV		
	tmp			<- which(grepl('AF103818',ph$tip.label))
	ph			<- reroot(ph, tmp, ph$edge.length[which(ph$edge[,2]==tmp)])	
	ph 			<- ladderize( ph )
	ph.node.bs						<- ph$node.label
	ph.node.bs[ph.node.bs=='Root']	<- NA
	ph.node.bs						<- as.numeric(ph.node.bs)
	ph.node.bs[is.na(ph.node.bs)]	<- 0
	ph$node.label					<- ph.node.bs
	ph$node.label[ph$node.label>1]	<- 1
	#	plot just the phylogeny
	pdf(file=file.path(indir,gsub('newick','pdf',infile)), w=20, h=120)
	plot(ph, cex=0.2)
	dev.off()
	
	stat.fun						<- hivc.clu.min.transmission.cascade
	#stat.fun						<- max
	dist.brl						<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=stat.fun)
	thresh.brl						<- 0.1
	thresh.bs						<- 0.7	
	clustering						<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs,retval="all")	
	tip.color						<- rep('black',Ntip(ph))
	tip.color[ grepl('PG14',ph$tip.label) ]	<- 'DarkRed'
	file							<- paste(outdir,"/", gsub('\\.newick',paste('_gd',100*thresh.brl,'bs',thresh.bs*100,'\\.pdf',sep=''), infile), sep='')
	invisible(hivc.clu.plot(ph, clustering[["clu.mem"]], cex.edge.incluster=3, tip.color=tip.color, file=file, pdf.scaley=100, show.tip.label=TRUE, pdf.width=30))
	
	ph.gdtr							<- cophenetic.phylo(ph)
	ph.mrca							<- mrca(ph)
	save(ph, dist.brl, ph.gdtr, ph.mrca, clustering, file=paste(indir,'/',gsub('\\.newick','\\.rda',infile),sep=''))
}

project.dualinfecions.UG.setup.windowlength<- function()
{
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/window_len'
	infiles	<- data.table(FILE= list.files(indir, pattern='^AlignedReadsInWindow', full.names=T, recursive=T))
	
	sd		<- infiles[,{
				#FILE<- 'AlignedReadsInWindow_800_to_1099.fasta'
				s	<- read.dna(FILE, format='fa')
				sd	<- data.table(TAXA=rownames(s))
				sd[, IND:= gsub('.bam','',gsub('_read.*','',TAXA))]
				sd[, COUNT:= as.numeric(regmatches(TAXA,regexpr('[0-9]+$',TAXA)))]	
				sd
			}, by=c('FILE')]
	sd		<- sd[, list(COUNT=sum(COUNT), UNIQUE=length(COUNT)), by=c('FILE','IND')]
	sd		<- subset(sd, !grepl('REF',IND))
	sd[, TYPE:= gsub('window_len/','',regmatches(FILE, regexpr('window_len/[^/]+', FILE))) ]
	set(sd, NULL, 'FILE', sd[,basename(FILE)])		
	sd[, FROM:= as.numeric(gsub('AlignedReadsInWindow_','',regmatches(FILE, regexpr('AlignedReadsInWindow_[0-9]+', FILE))))]
	sd		<-melt(sd, measure.vars=c('COUNT','UNIQUE'))
	
	ggplot(sd, aes(x=FROM, y=value, colour=TYPE, group=TYPE)) + geom_step() + facet_grid(variable~IND, scales='free_y') + scale_y_log10()
	ggsave(file='~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_UG_gag_selecthelp_windows.pdf', w=50, h=12, limitsize = FALSE)
}


project.RakaiAll.setup.prior.170301<- function()
{
	df	<- as.data.table(expand.grid(T=c(2,5,8), N0=c(0.1,1,10,100)))
	df	<- df[, list(X=seq(1e-3,1-1e-3,by=1e-3), Y=dbeta(seq(1e-3,1-1e-3,by=1e-3), N0/T, N0-N0/T)), by=c('T','N0')]
	set(df, NULL, 'T', df[,paste0('T=',T)])
	set(df, NULL, 'N0', df[,paste0('n0=',N0)])
	set(df, NULL, 'N0', df[,factor(N0, levels=c('n0=0.1','n0=1','n0=10','n0=100'))])
	ggplot(df,aes(x=X,y=Y)) + 
			geom_line() +
			theme_bw() +
			facet_grid(T~N0) +
			labs(x='pi_t', y='prior density')
	ggsave(file='~/Dropbox (SPH Imperial College)/OR_Work/2016/2016_Rakai_Couples/prior_beta.pdf', w=8, h=8)
	
	n0<- 3; T<- 3; n<- c(3,4,5)	
	a<- (n0+T*n)/T; b<- ((T-1)*n0+T*(sum(n)-n))/T
	pbeta(1/T,a,b,lower.tail=FALSE)
	#0.2611934 0.4755005 0.6898075
	
	n0<- 3; T<- 3; n<- c(1e3,1e3,1e3)	
	a<- (n0+T*n)/T; b<- ((T-1)*n0+T*(sum(n)-n))/T
	pbeta(1/T,a,b,lower.tail=FALSE)
	#0.498284 0.498284 0.498284
	
	n0<- 3; T<- 3; n<- c(1,0,0)	
	a<- (n0+T*n)/T; b<- ((T-1)*n0+T*(sum(n)-n))/T
	pbeta(1/T,a,b,lower.tail=FALSE)
	#0.7407407 0.2962963 0.2962963
	
	n0<- c(0.1,1,3,5,100); T<- 3; n<- c(1)	
	pbeta(1/T, (n0+T*n)/T, (T-1)*n0/T, lower.tail=FALSE)
	#0.9749602 0.8457642 0.7407407 0.6952502 0.5468373
	
	n0<- 3; T<- c(2,3,4,5); n<- c(1)	
	pbeta(1/T, (n0+T*n)/T, (T-1)*n0/T, lower.tail=FALSE)
	#0.7122066 0.7407407 0.7654499 0.7857324
	n0<- 3; T<- c(2,3,4,5); n<- c(1)	
	pbeta(.5, (n0+T*n)/T, (T-1)*n0/T, lower.tail=FALSE)
	#0.7122066 0.5000000 0.3904129 0.3275931
	
	
	phsc.find.n0.aux<- function(n0, phsc.T=3, phsc.n=1)
	{
		abs( pbeta(1/phsc.T+(1-1/phsc.T)/(phsc.T+1), (n0+phsc.T*phsc.n)/phsc.T, (phsc.T-1)*n0/phsc.T, lower.tail=FALSE)-0.5 ) 	
	}
	phsc.find.n0.aux(c(0.1,1,3,5,10,15,20), 3, 1)
	optimize(phsc.find.n0.aux, c(.1,10), phsc.T=3, phsc.n=1)
	
	N_TYPE<-3; NEFF<- 3; KEFF<- NEFF*2/3; n0<- c(0.1,0.5,1,2,3,4)
	pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), n0/N_TYPE+KEFF, n0*(1-1/N_TYPE)+NEFF-KEFF, lower.tail=FALSE)
	
	
	n0	<- 1; T<- 3; n<- seq(3,6,.5); k<- n-1	
	pbeta(1/T, (n0+T*k)/T, ((T-1)*n0+T*(n-k))/T, lower.tail=FALSE)
	
	n0	<- 2; T<- 2; n<- seq(3,6,1); k<- n-1
	a<- (n0+T*k)/T
	b<- ((T-1)*n0+T*(n-k))/T
	(a-1)/(a+b-2)
	(n0+T*(k-1))/ (T*(n+n0-2))
	
	n.states<- 2; keff<- 3; neff<- 4; confidence.cut<- 0.66
	phsc.find.n0.aux<- function(n0, n.states, keff, neff, confidence.cut)
	{
		abs( (n0+n.states*(keff-1))/ (n.states*(neff+n0-2)) - confidence.cut )
	}	
	ans	<- optimize(phsc.find.n0.aux, c(.001,100), n.states=n.states, keff=keff, neff=neff, confidence.cut=confidence.cut)
	ans	<- round(ans$minimum, d=4)
	plot(seq(0,1,1e-3), dbeta(seq(0,1,1e-3),1/n.states*ans,1/n.states*ans), type='l')
}

project.RakaiAll.setup.consensuses<- function()
{
	#get data to Tanya to calculate missing consensuses
	indir	<- '/Volumes/exfat/PANGEA/PANGEA_mapout/RCCS_597_FastQZ'
	infiles	<- data.table(F=list.files(indir))
	infiles[, SID:= gsub('^([0-9]+_[0-9]+_[0-9]+)_.*','\\1',F)]
	tmp		<- unique(subset(infiles, select=SID))
		
	
	indir	<- '/Volumes/exfat/PANGEA/PANGEA_mapout/data_csvetc'
	infiles	<- data.table(F=list.files(indir, full.names=TRUE, pattern='BaseFreqs|coords'))
	infiles	<- subset(infiles, grepl('^[0-9]+_[0-9]+_[0-9]+', basename(F)) & !grepl('remap',basename(F)) )
	infiles[, SID:= gsub('^([0-9]+_[0-9]+_[0-9]+)_.*','\\1',basename(F))]
	infiles[, BASEFREQ:= factor(grepl('BaseFreqs', basename(F)), levels=c(TRUE,FALSE), labels=c('BaseFreq','Coord'))]
	#	missing BaseFreq files
	tmp2	<- data.table(SID=sort(unique(setdiff(tmp[, SID], infiles[, SID]))))
	#	all Rakai 597 BaseFreq files that I have locally
	infiles	<- merge(tmp, infiles, by='SID')
	
	#	check missing have no bam
	indir	<- '/Volumes/exfat/PANGEA/PANGEA_mapout/data'
	tmp		<- data.table(F=list.files(indir, full.names=TRUE, pattern='bam$'))
	tmp		<- subset(tmp, grepl('^[0-9]+_[0-9]+_[0-9]+', basename(F)) & !grepl('remap',basename(F)) )
	tmp[, SID:= gsub('^([0-9]+_[0-9]+_[0-9]+).*','\\1',basename(F))]
	stopifnot(!nrow(merge(tmp2, tmp, by='SID')))
	
	#	zip
	outfile	<- '/Volumes/exfat/PANGEA/PANGEA_mapout/RCCS_597_BaseFreqs.zip'
	invisible( subset(infiles, BASEFREQ=='BaseFreq')[, file.copy(F, file.path(dirname(outfile),'tmp42',basename(F))), by='F'] )
	invisible( zip( outfile, file.path(dirname(outfile),'tmp42'), flags = "-umr9XTjq") )
	outfile	<- '/Volumes/exfat/PANGEA/PANGEA_mapout/RCCS_597_Coord.zip'
	invisible( subset(infiles, BASEFREQ=='Coord')[, file.copy(F, file.path(dirname(outfile),'tmp42',basename(F))), by='F'] )
	invisible( zip( outfile, file.path(dirname(outfile),'tmp42'), flags = "-umr9XTjq") )
	
}

project.RakaiAll.setup.batchno.170301<- function()
{
	require(big.phylo)
	require(data.table)
	
	#	batch=15 merge=2
	indir	<- '/Users/Oliver/duke/tmp/pty_17-03-08-11-12-50'
	infiles	<- list.files(indir, pattern='^ptyr1_.*tree', full.names=TRUE)	
	ntips	<- sapply(infiles, function(x) Ntip(read.tree(x)))
	#	batch=15 merge=0
	indir	<- '/Users/Oliver/duke/tmp/pty_17-03-08-11-30-25'
	infiles	<- list.files(indir, pattern='^ptyr1_.*tree', full.names=TRUE)	
	ntips	<- sapply(infiles, function(x) Ntip(read.tree(x)))
	#	batch=75 merge=2
	indir	<- '/Users/Oliver/duke/tmp/pty_17-03-08-11-52-06'
	infiles	<- list.files(indir, pattern='^ptyr1_.*tree', full.names=TRUE)	
	ntips	<- sapply(infiles, function(x) Ntip(read.tree(x)))
	
}

project.RakaiAll.setup.RAxMLmodel.submit.170301<- function()
{
	require(big.phylo)
	require(data.table)
	hpc.load	<- "module load intel-suite/2015.1 mpi R/3.3.2"
	hpc.nproc	<- 8
	hpc.mem		<- "5900mb"	
	hpc.walltime<- 300 
	
	#indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/RAxML_model_test'
	indir	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RAxML_model_test'
	
	infiles	<- data.table(F=list.files(indir, pattern='fasta$', full.names=TRUE))
	infiles	<- subset(infiles, !grepl('ptyr1_',F))
	infiles[, {				
				cmd			<- cmd.jmodeltest(F, pr.args='-f -i -g 4 -s 3 -DT -S NNI -t ML', nproc=hpc.nproc)
				tmp			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=hpc.walltime, hpc.q="pqeelab", hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=hpc.load)							
				cmd			<- paste(tmp,cmd,sep='\n')
				cat(cmd)
				stop()
				outfile		<- paste("jm",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
				cmd.hpccaller(indir, outfile, cmd)				
			}, by='F']
		
}

project.RakaiAll.setup.RAxMLmodel.evaluate.170301<- function()
{
	require(big.phylo)
	require(data.table)
	require(ggplot2)
	require(scales)
	require(gamlss)
	
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/RAxML_model_test'
	infiles	<- data.table(F=list.files(indir, pattern='jmodeltest$', full.names=TRUE))
	infiles[, PTY_RUN:= gsub('.*ptyr([0-9]+)_InWindow_([0-9]+)_to_([0-9]+).*','\\1',F)]
	infiles[, W_FROM:= gsub('.*ptyr([0-9]+)_InWindow_([0-9]+)_to_([0-9]+).*','\\2',F)]
	infiles[, W_TO:= gsub('.*ptyr([0-9]+)_InWindow_([0-9]+)_to_([0-9]+).*','\\3',F)]
	
	dt	<- infiles[, {
				#F	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/RAxML_model_test/ptyr1_InWindow_2250_to_2499.fasta.jmodeltest'
				cat('process',basename(F),'\n')
				dt	<- readLines(F)
				taxa<- as.numeric(gsub('.*number of sequences: ([0-9]+).*','\\1',dt[which(grepl('number of sequences',dt))]))
				
				
				dt	<- dt[which(grepl('* DT MODEL SELECTION : Selection uncertainty',dt)):which(grepl('* DT MODEL SELECTION : Confidence interval',dt))]
				dt	<- dt[(which(grepl('------------------------------------',dt))[1]+1): (which(grepl('------------------------------------',dt))[2]-1)]
				dt	<- as.data.table(t(sapply(strsplit(dt,' +'), function(x) x)))
				setnames(dt, paste0('V',1:7), toupper(c('Model','lnL','K','DT','delta','weight','cumWeight')))
				for(x in setdiff(colnames(dt),'MODEL'))
					set(dt, NULL, x, as.numeric(dt[[x]]))
				dt[, N_TAXA:=taxa]
				dt[, BEST_MODEL:= MODEL[which.max(WEIGHT)]]
				dt
			}, by=c('PTY_RUN','W_FROM','W_TO')]
		
	ggplot(dt, aes(x=W_FROM, y=N_TAXA, colour=BEST_MODEL)) + 
			geom_point() + 			
			theme_bw() 
	ggsave(file=file.path(indir,'Best_Subst_Model_acrossgenome_by_ntaxa.pdf'), w=5, h=4)
	tmp	<- unique(dt,by=c('PTY_RUN','W_FROM','W_TO'))
	ggplot(tmp, aes(x=BEST_MODEL)) +
			geom_bar() + 			
			theme_bw() 
	ggsave(file=file.path(indir,'Best_Subst_Model_number.pdf'), w=5, h=4)
	ggplot(tmp, aes(x=W_FROM, y=WEIGHT, colour=BEST_MODEL)) + 
			geom_point() + 	
			scale_y_continuous(labels=percent) +
			labs(y='model prob in model averaging') +
			theme_bw() 
	ggsave(file=file.path(indir,'Best_Subst_Model_acrossgenome_by_modelprob.pdf'), w=5, h=4)
	ggplot(tmp, aes(x=W_FROM, y=PTY_RUN, colour=BEST_MODEL)) + 
			geom_point() + 				
			labs(y='phyloscanner run (different patients)') +
			theme_bw() 
	ggsave(file=file.path(indir,'Best_Subst_Model_acrossgenome_by_selectedpatients.pdf'), w=5, h=4)	
	tmp	<- dt[, list(SW=sum(WEIGHT)), by='MODEL']
	ggplot(tmp, aes(x=MODEL, y=SW)) +
			geom_bar(stat='identity') + 			
			theme_bw() +
			coord_flip() +
			labs(y='sum of model probabilities in model averaging over all runs')
	ggsave(file=file.path(indir,'Best_Subst_Model_summodelprob.pdf'), w=5, h=4)	
	tmp	<- subset(dt, MODEL%in%c('HKY+G','GTR+G'))
	ggplot(tmp, aes(x=N_TAXA, y=WEIGHT, colour=MODEL)) + 
			geom_point() +
			theme_bw() +
			scale_y_continuous(labels=percent) +
			labs(y='model prob in model averaging') +
			geom_smooth(se=TRUE, fill='grey80', span=2)
	ggsave(file=file.path(indir,'Best_Subst_Model_trends_by_ntaxa.pdf'), w=5, h=4)	 
}

project.RakaiAll.setup.phyloscanner.170301.stagetwo.preprocessing	<- function()
{
	require(data.table)
	#
	#	from every phyloscanner run, select pairs that are closely related 
	#		
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_tb_blnormed'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170421.rda'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl5_blnormed'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170428.rda'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170428_cl3.rda'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s25_allbatch_sk20_tb_blnormed'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170704_cl35.rda'
	indir	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s25_allbatch_sk20_tb_blnormed_save'
	outfile	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/close_pairs_170704_cl35.rda'

	
	infiles	<- data.table(F=list.files(indir, pattern='pairwise_relationships.rda', full.names=TRUE))
	infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(F)))]
	setkey(infiles, PTY_RUN)
	rtp	<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s25_allbatch_sk20_tb_blnormed/ptyr667_pairwise_relationships.rda'
				#cat(PTY_RUN,'\n')
				load(F)
				group	<- 'TYPE_PAIR_DI2'
				rtp		<- subset(rplkl, GROUP==group)[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, TYPE_MLE=='close')
				rtp		<- merge(rtp, subset(rplkl, GROUP==group & TYPE=='close'), by=c('ID1','ID2'), all.x=1)
				#tmp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIRSCORE_DI' & TYPE=='close'), by=c('ID1','ID2'), all.x=1)
				#rtp		<- rbind(rtp, tmp)				
				#rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				rtp[, POSTERIOR_SCORE:=(POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-2)]
				rtp				
			}, by=c('PTY_RUN')]		
	rn	<- infiles[, {
				load(F)
				group	<- 'TYPE_PAIR_DI2'
				rtp		<- subset(rplkl, GROUP==group & TYPE=='close', select=c(ID1, ID2, NEFF))
				rtp				
			}, by=c('PTY_RUN')]			
	dmin	<- infiles[,{
				#cat(PTY_RUN,'\n')
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl5_blnormed/ptyr1_pairwise_relationships.rda'
				load(F)
				tmp		<- dwin[, list(		PATRISTIC_DISTANCE_MIN=min(PATRISTIC_DISTANCE), 
								PATRISTIC_DISTANCE_MEAN=mean(PATRISTIC_DISTANCE),
								PATRISTIC_DISTANCE_MEDIAN=median(PATRISTIC_DISTANCE)), by=c('ID1','ID2')]
				tmp
			}, by='PTY_RUN']	
	save(rtp, rn, dmin, file=outfile)	
}

project.RakaiAll.setup.phyloscanner.170301.stagetwo.preprocessing2	<- function()
{
	require(data.table)
	#
	#	from every phyloscanner run, select pairs that are closely related 
	#		
	indir	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s25_allbatch_sk20_tb_blnormed_save'
	outfile	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/close_pairs_170704_cl35_withtopology.rda'
		
	infiles	<- data.table(F=list.files(indir, pattern='pairwise_relationships.rda', full.names=TRUE))
	infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(F)))]
	setkey(infiles, PTY_RUN)
	
	#	find runs in which ID1 and ID2 have highest N
	min.read	<- 30
	neff.cut	<- 3
	rtp		<- infiles[, 	{
								load(F)
								ans	<- dwin[, list(N=length(which(ID1_R>=min.read & ID2_R>=min.read))), by=c('ID1','ID2')]
								ans				
							}, by=c('PTY_RUN')]	
	rtp		<- rtp[, list(PTY_RUN= PTY_RUN[which.max(N)], N=max(N)), by=c('ID1','ID2')]
	rtp		<- subset(rtp, N>=neff.cut)
	#	for each pair, now calculate the number of ancestral, intermingled etc topologies,
	#	and the median distance
	tmp		<- unique(rtp$PTY_RUN)
	infiles	<- subset(infiles, PTY_RUN%in%tmp)
	ans		<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s25_allbatch_sk20_tb_blnormed/ptyr667_pairwise_relationships.rda'
				#cat(PTY_RUN,'\n')
				load(F)
				pty.run <- PTY_RUN
				tmp		<- subset(rtp, PTY_RUN==pty.run)
				df		<- merge(dwin, tmp, by=c('ID1','ID2'))
				#	define topology types
				df[, TYPE_BASIC:=NA_character_]
				tmp		<- df[, which(PATHS_12>0 & PATHS_21==0 & ADJACENT)]				
				set(df, tmp, 'TYPE_BASIC', 'chain_12')
				tmp		<- df[, which(PATHS_12==0 & PATHS_21>0 & ADJACENT)]				
				set(df, tmp, 'TYPE_BASIC', 'chain_21')	
				tmp		<- df[, which(PATHS_12>0 & PATHS_21>0 & ADJACENT)]				
				set(df, tmp, 'TYPE_BASIC', 'intermingled')
				tmp		<- df[, which(PATHS_12==0 & PATHS_21==0 & ADJACENT)]
				set(df, tmp, 'TYPE_BASIC', 'sibling')
				tmp		<- df[, which(PATHS_12==0 & PATHS_21==0 & !ADJACENT )]				
				set(df, tmp, 'TYPE_BASIC', 'other')
				tmp		<- df[, which(PATHS_12>0 & PATHS_21==0 & !ADJACENT )]				
				set(df, tmp, 'TYPE_BASIC', 'other')	
				tmp		<- df[, which(PATHS_12==0 & PATHS_21>0 & !ADJACENT )]				
				set(df, tmp, 'TYPE_BASIC', 'other')	
				tmp		<- df[, which(PATHS_12>0 & PATHS_21>0 & !ADJACENT )]				
				set(df, tmp, 'TYPE_BASIC', 'other')	
				stopifnot( !any(is.na(df$TYPE_BASIC)) )
				#	compute median distance
				tmp		<- df[, list(N=length(W_FROM), PATRISTIC_DISTANCE_MEDIAN=median(PATRISTIC_DISTANCE)), by=c('ID1','ID2')]
				#	count each type
				df		<- df[, list(K=length(W_FROM)), by=c('ID1','ID2','TYPE_BASIC')]
				df		<- merge(df, tmp, by=c('ID1','ID2'))
				list(ID1= df$ID1, ID2=df$ID2, TYPE_BASIC=df$TYPE_BASIC, K=df$K, N=df$N, PATRISTIC_DISTANCE_MEDIAN=df$PATRISTIC_DISTANCE_MEDIAN)				
			}, by=c('PTY_RUN')]				
	save(ans, file=outfile)	
}

project.RakaiAll.setup.phyloscanner.170301.stagetwo	<- function()
{
	require(igraph)
	
	infile		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170421.rda'
	infile		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170428_cl3.rda'
	load(infile)	
	rtp			<- subset(rtp, GROUP=='TYPE_PAIR_DI' & POSTERIOR_SCORE>0.5)
	
	#	get transmission chains with igraph
	tmp			<- subset(rtp, select=c(ID1, ID2))			
	tmp			<- graph.data.frame(tmp, directed=FALSE, vertices=NULL)
	rtc			<- data.table(ID=V(tmp)$name, CLU=clusters(tmp, mode="weak")$membership)	
	tmp2		<- rtc[, list(CLU_SIZE=length(ID)), by='CLU']
	setkey(tmp2, CLU_SIZE)
	tmp2[, IDCLU:=rev(seq_len(nrow(tmp2)))]
	rtc			<- subset( merge(rtc, tmp2, by='CLU'))
	rtc[, CLU:=NULL]
	setkey(rtc, IDCLU)
	
	#tmp			<- subset(rtc, IDCLU==1)[, ID]
	#tmp			<- subset(rtp, ID1%in%tmp | ID2%in%tmp, select=c(ID1, ID2))
	#tmp			<- graph.data.frame(tmp, directed=FALSE, vertices=NULL)
	#par(mai=c(0,0,1,0))
	#plot(tmp, vertex.size=2, vertex.label.cex=0.25, edge.arrow.size=0.5, layout=layout.fruchterman.reingold(tmp, niter=1e3) )	
	
	#	add in couples that are not a close pair
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75.rda')
	tmp		<- copy(pty.runs)	
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75_part2.rda')
	tmp		<- rbind(tmp, pty.runs)
	tmp		<- sort(unique(as.character(tmp[['RID']])))	
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	rp		<- unique(subset(rp, FEMALE_RID%in%tmp & MALE_RID%in%tmp, select=c(FEMALE_RID,MALE_RID)))
	setnames(rp, 'FEMALE_RID', 'ID')
	rp		<- merge(rp, subset(rtc, select=c(ID, IDCLU)), all.x=1, by='ID')
	setnames(rp, c('ID','IDCLU','MALE_RID'), c('FEMALE_RID', 'FEMALE_IDCLU','ID'))
	rp		<- merge(rp, subset(rtc, select=c(ID, IDCLU)), all.x=1, by='ID')
	setnames(rp, c('ID','IDCLU'), c('MALE_RID', 'MALE_IDCLU'))
	tmp		<- rp[, which(is.na(FEMALE_IDCLU) & !is.na(MALE_IDCLU))]
	set(rp, tmp, 'FEMALE_IDCLU', rp[tmp, MALE_IDCLU])
	tmp		<- rp[, which(!is.na(FEMALE_IDCLU) & is.na(MALE_IDCLU))]
	set(rp, tmp, 'MALE_IDCLU', rp[tmp, FEMALE_IDCLU])
	tmp		<- rp[, which(is.na(FEMALE_IDCLU) & is.na(MALE_IDCLU))]
	set(rp, tmp, c('FEMALE_IDCLU','MALE_IDCLU'), rtc[, max(IDCLU)]+seq_along(tmp))
	setnames(rp, c('MALE_IDCLU'), c('IDCLU'))
	rp		<- subset(melt(rp, id.vars=c('IDCLU'), measure.vars=c('MALE_RID','FEMALE_RID'), value.name='ID'), select=c(ID, IDCLU))
	set(rtc, NULL, 'CLU_SIZE', NULL)
	rtc		<- unique(rbind(rtc, rp))
	tmp		<- rtc[, list(CLU_SIZE=length(ID)), by='IDCLU']
	setkey(tmp, CLU_SIZE)
	tmp[, IDCLU2:=rev(seq_len(nrow(tmp)))]	
	rtc		<- merge(rtc, tmp, by='IDCLU')
	rtc[, IDCLU:=NULL]
	setnames(rtc, 'IDCLU2', 'IDCLU')
	
	#	1547 close individuals 
	#	added 493 unrelated individuals that are in couple
	
	#	merge up to 8 individuals into the same phyloscanner run
	tn			<- 8
	tmp			<- unique(subset(rtc, select=c(IDCLU, CLU_SIZE)))		
	tmp[, PTY_RUN:= IDCLU]
	tmp[, PTY_SIZE:= CLU_SIZE]
	for(i in rev(seq_len(nrow(tmp))))
	{
		if(tmp[i,PTY_SIZE<tn])
		{
			for(j in 1:i)
			{
				if( (tmp[j,PTY_SIZE]+tmp[i,PTY_SIZE])<=tn)
				{
					set(tmp, i, 'PTY_RUN', tmp[j, PTY_RUN])					
					set(tmp, j, 'PTY_SIZE',  tmp[i, PTY_SIZE]+ tmp[j, PTY_SIZE])
					set(tmp, i, 'PTY_SIZE', 0L)
					break
				}
			}
		}
	}
	# manually add IDCLUD 591 to PTY_RUN 592	
	set(tmp, tmp[, which(PTY_RUN==591)], 'PTY_RUN', 592L)	
	tmp			<- tmp[, list(IDCLU=IDCLU, CLU_SIZE=CLU_SIZE, PTY_SIZE=sum(CLU_SIZE)), by=c('PTY_RUN')]
	rtc			<- merge(rtc, tmp, by=c('IDCLU','CLU_SIZE'))
	#	for every individual in a cluster
	#	find 15 closest others
	sn		<- 15
	setkey(dmin, ID1, ID2)	
	dcl		<- rtc[, {
				tmp		<- subset(dmin, ID1==ID)
				tmp2	<- subset(dmin, ID2==ID)
				setnames(tmp2, c('ID1','ID2'), c('ID2','ID1'))
				tmp		<- rbind(tmp, tmp2)
				tmp		<- tmp[, list(PATRISTIC_DISTANCE=mean(PATRISTIC_DISTANCE_MIN)), by=c('ID1','ID2')]
				setkey(tmp, PATRISTIC_DISTANCE)
				tmp		<- tmp[seq_len(sn),]
				list(ID_CLOSE= tmp[,ID2], PATRISTIC_DISTANCE=tmp[, PATRISTIC_DISTANCE])				
			}, by=c('ID','CLU_SIZE','IDCLU')]
	#	find as many close other individuals as needed to fill the phyloscanner run
	#	up to 15 individuals
	dcl		<- subset(dcl, !is.na(ID_CLOSE))
	tmp		<- dcl[, list(ID_CLOSE=setdiff(ID_CLOSE, ID)), by='IDCLU']
	tmp		<- merge(unique(subset(dcl, select=c('IDCLU','ID_CLOSE','CLU_SIZE','PATRISTIC_DISTANCE'))), tmp, by=c('IDCLU','ID_CLOSE'))
	tmp		<- tmp[, list(PATRISTIC_DISTANCE=min(PATRISTIC_DISTANCE)), by=c('IDCLU','ID_CLOSE','CLU_SIZE')]
	tmp		<- merge( unique(subset(rtc, select=c('IDCLU','PTY_RUN','PTY_SIZE'))), tmp, by='IDCLU' )
	setkey(tmp, PTY_RUN, IDCLU, PATRISTIC_DISTANCE)	
	tmp		<- merge(tmp, tmp[, list(NCLU=length(unique(IDCLU))), by='PTY_RUN'], by='PTY_RUN')	
	tmp		<- tmp[, {
				z<- seq_len( max(0,ceiling((sn-PTY_SIZE[1])/NCLU[1])) )
				list(ID_CLOSE=ID_CLOSE[z], PATRISTIC_DISTANCE=PATRISTIC_DISTANCE[z])
			}, by=c('PTY_RUN','IDCLU','CLU_SIZE')]
	setnames(tmp, 'ID_CLOSE', 'ID')
	tmp[, ID_TYPE:='background']
	rtc[, ID_TYPE:='target']
	
	rtc[, PTY_SIZE:=NULL]
	rtc		<- rbind(rtc, tmp, fill=TRUE, use.names=TRUE)
	rtc		<- merge(rtc, rtc[, list(PTY_SIZE=length(ID)), by='PTY_RUN'], by='PTY_RUN')
	tmp		<- unique(subset(rtc, select=c(PTY_RUN, PTY_SIZE)))
	tmp[, PTY_RUN2:= seq_len(nrow(tmp))]
	rtc		<- merge(rtc, subset(tmp, select=c('PTY_RUN','PTY_RUN2')), by='PTY_RUN' )
	rtc[, PTY_RUN:=NULL]
	setnames(rtc, 'PTY_RUN2', 'PTY_RUN')
	
	#	in total 2320 individuals
	
	#
	#	set up new phyloscanner runs object
	#	find Sanger IDs in previous pty.runs
	#
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75.rda')
	tmp		<- copy(pty.runs)	
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75_part2.rda')
	tmp		<- rbind(tmp, pty.runs)
	set(tmp, NULL, c('BATCH','PTY_RUN'), NULL)
	tmp		<- unique(tmp)
	pty.runs<- NULL
	setnames(rtc, c('IDCLU','CLU_SIZE','ID','ID_TYPE','PATRISTIC_DISTANCE','PTY_SIZE'), c('TRCLU_ID','TRCLU_SIZE','RID','RID_TYPE','RID_PATRISTIC_DISTANCE','PTY_SIZE_RID'))
	pty.runs<- merge(rtc, tmp, by='RID')
	setkey(pty.runs, PTY_RUN, RID_TYPE, RID)
	pty.runs<- unique(pty.runs, by=c('RID','SID','PTY_RUN'))	
	pty.runs<- merge(pty.runs, pty.runs[, list(PTY_SIZE_SID=length(SID)), by='PTY_RUN'], by='PTY_RUN')	
	pty.runs<- merge(pty.runs, pty.runs[, list(SID=SID, RENAME_SID=paste0(RID,'_fq',seq_along(SID))), by=c('PTY_RUN','RID')], by=c('PTY_RUN','RID','SID'))
	
	save(pty.runs, file= '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_stagetwo.rda')
}

project.RakaiAll.setup.phyloscanner.170704.stagethree.andtwo	<- function()
{
	require(igraph)
	#
	#	revisit the huge cluster 1 because we cannot process it
	#	the main point is that all pairs with NEFF<3 are deemed as having insufficient data
	#	this introduces additional sparseness
	#
	infile	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170704_stagetwo.rda'
	load(infile)	
	run1	<- subset(pty.runs, PTY_RUN==1)
	
	infile		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170421.rda'
	infile		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170428_cl3.rda'
	infile		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170704_cl35.rda'
	load(infile)	
	
	rtp			<- subset(rtp, GROUP=='TYPE_PAIR_DI2' & POSTERIOR_SCORE>=0.6)
	tmp			<- unique(run1$RID)
	tmp2		<- subset(rtp, !(ID1%in%tmp | ID2%in%tmp))
	rtp			<- subset(rtp, (ID1%in%tmp | ID2%in%tmp) & NEFF>=3)
	rtp			<- rbind(rtp, tmp2)
	
	#	how many pairs remaining
	tmp			<- unique(subset(rtp, select=c(ID1, ID2)))
	stopifnot(nrow(subset(tmp, ID1>ID2))==0)
	
	#	get gender of these
	infile.seqinfo	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75_part2.rda'
	load(infile.seqinfo)	
	tmp				<- merge(tmp, unique(data.table(ID1=dc$RID, ID1_SEX=dc$SEX)), by='ID1', all.x=TRUE)
	tmp				<- merge(tmp, unique(data.table(ID2=dc$RID, ID2_SEX=dc$SEX)), by='ID2', all.x=TRUE)
	tmp[, SXO:= paste0(ID1_SEX,ID2_SEX)]
	
	#	get transmission chains with igraph
	tmp			<- subset(rtp, select=c(ID1, ID2))			
	tmp			<- graph.data.frame(tmp, directed=FALSE, vertices=NULL)
	rtc			<- data.table(ID=V(tmp)$name, CLU=clusters(tmp, mode="weak")$membership)	
	tmp2		<- rtc[, list(CLU_SIZE=length(ID)), by='CLU']
	setkey(tmp2, CLU_SIZE)
	tmp2[, IDCLU:=rev(seq_len(nrow(tmp2)))]
	rtc			<- subset( merge(rtc, tmp2, by='CLU'))
	rtc[, CLU:=NULL]
	setkey(rtc, IDCLU)	
	rtc[, max(IDCLU)]	#482
	
}

project.RakaiAll.setup.phyloscanner.170704.stagethree	<- function()
{
	require(igraph)
	#
	#	revisit the huge cluster 1 because we cannot process it
	#	the main point is that all pairs with NEFF<3 are deemed as having insufficient data
	#	this introduces additional sparseness
	#
	infile	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170704_stagetwo.rda'
	load(infile)	
	run1	<- subset(pty.runs, PTY_RUN==1)
	
	infile		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170421.rda'
	infile		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170428_cl3.rda'
	infile		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170704_cl35.rda'
	load(infile)	
	
	rtp			<- subset(rtp, GROUP=='TYPE_PAIR_DI2' & POSTERIOR_SCORE>=0.6)
	tmp			<- unique(run1$RID)
	rtp			<- subset(rtp, ID1%in%tmp | ID2%in%tmp)
	rtp			<- subset(rtp, NEFF>=3)
		
	#	get transmission chains with igraph
	tmp			<- subset(rtp, select=c(ID1, ID2))			
	tmp			<- graph.data.frame(tmp, directed=FALSE, vertices=NULL)
	rtc			<- data.table(ID=V(tmp)$name, CLU=clusters(tmp, mode="weak")$membership)	
	tmp2		<- rtc[, list(CLU_SIZE=length(ID)), by='CLU']
	setkey(tmp2, CLU_SIZE)
	tmp2[, IDCLU:=rev(seq_len(nrow(tmp2)))]
	rtc			<- subset( merge(rtc, tmp2, by='CLU'))
	rtc[, CLU:=NULL]
	setkey(rtc, IDCLU)
	
	#	add in couples that are not a close pair
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75.rda')
	tmp		<- copy(pty.runs)	
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75_part2.rda')
	tmp		<- rbind(tmp, pty.runs)
	tmp		<- sort(unique(as.character(tmp[['RID']])))	
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	rp		<- unique(subset(rp, FEMALE_RID%in%tmp & MALE_RID%in%tmp, select=c(FEMALE_RID,MALE_RID)))
	setnames(rp, 'FEMALE_RID', 'ID')
	rp		<- merge(rp, subset(rtc, select=c(ID, IDCLU)), all.x=1, by='ID')
	setnames(rp, c('ID','IDCLU','MALE_RID'), c('FEMALE_RID', 'FEMALE_IDCLU','ID'))
	rp		<- merge(rp, subset(rtc, select=c(ID, IDCLU)), all.x=1, by='ID')
	setnames(rp, c('ID','IDCLU'), c('MALE_RID', 'MALE_IDCLU'))
	tmp		<- rp[, which(is.na(FEMALE_IDCLU) & !is.na(MALE_IDCLU))]
	set(rp, tmp, 'FEMALE_IDCLU', rp[tmp, MALE_IDCLU])
	tmp		<- rp[, which(!is.na(FEMALE_IDCLU) & is.na(MALE_IDCLU))]
	set(rp, tmp, 'MALE_IDCLU', rp[tmp, FEMALE_IDCLU])
	tmp		<- rp[, which(is.na(FEMALE_IDCLU) & is.na(MALE_IDCLU))]
	set(rp, tmp, c('FEMALE_IDCLU','MALE_IDCLU'), rtc[, max(IDCLU)]+seq_along(tmp))
	setnames(rp, c('MALE_IDCLU'), c('IDCLU'))
	rp		<- subset(melt(rp, id.vars=c('IDCLU'), measure.vars=c('MALE_RID','FEMALE_RID'), value.name='ID'), select=c(ID, IDCLU))
	set(rtc, NULL, 'CLU_SIZE', NULL)
	rtc		<- unique(rbind(rtc, rp))
	tmp		<- rtc[, list(CLU_SIZE=length(ID)), by='IDCLU']
	setkey(tmp, CLU_SIZE)
	tmp[, IDCLU2:=rev(seq_len(nrow(tmp)))]	
	rtc		<- merge(rtc, tmp, by='IDCLU')
	rtc[, IDCLU:=NULL]
	setnames(rtc, 'IDCLU2', 'IDCLU')
		
	#	merge up to 8 individuals into the same phyloscanner run
	tn			<- 8
	tmp			<- unique(subset(rtc, select=c(IDCLU, CLU_SIZE)))		
	tmp[, PTY_RUN:= IDCLU]
	tmp[, PTY_SIZE:= CLU_SIZE]
	for(i in rev(seq_len(nrow(tmp))))
	{
		if(tmp[i,PTY_SIZE<tn])
		{
			for(j in 1:i)
			{
				if( (tmp[j,PTY_SIZE]+tmp[i,PTY_SIZE])<=tn)
				{
					set(tmp, i, 'PTY_RUN', tmp[j, PTY_RUN])					
					set(tmp, j, 'PTY_SIZE',  tmp[i, PTY_SIZE]+ tmp[j, PTY_SIZE])
					set(tmp, i, 'PTY_SIZE', 0L)
					break
				}
			}
		}
	}
	tmp			<- tmp[, list(IDCLU=IDCLU, CLU_SIZE=CLU_SIZE, PTY_SIZE=sum(CLU_SIZE)), by=c('PTY_RUN')]
	rtc			<- merge(rtc, tmp, by=c('IDCLU','CLU_SIZE'))
	#	for every individual in a cluster
	#	find 15 closest others
	sn		<- 15
	setkey(dmin, ID1, ID2)	
	dcl		<- rtc[, {
				tmp		<- subset(dmin, ID1==ID)
				tmp2	<- subset(dmin, ID2==ID)
				setnames(tmp2, c('ID1','ID2'), c('ID2','ID1'))
				tmp		<- rbind(tmp, tmp2)
				tmp		<- tmp[, list(PATRISTIC_DISTANCE=mean(PATRISTIC_DISTANCE_MIN)), by=c('ID1','ID2')]
				setkey(tmp, PATRISTIC_DISTANCE)
				tmp		<- tmp[seq_len(sn),]
				list(ID_CLOSE= tmp[,ID2], PATRISTIC_DISTANCE=tmp[, PATRISTIC_DISTANCE])				
			}, by=c('ID','CLU_SIZE','IDCLU')]
	#	find as many close other individuals as needed to fill the phyloscanner run
	#	up to 15 individuals
	dcl		<- subset(dcl, !is.na(ID_CLOSE))
	tmp		<- dcl[, list(ID_CLOSE=setdiff(ID_CLOSE, ID)), by='IDCLU']
	tmp		<- merge(unique(subset(dcl, select=c('IDCLU','ID_CLOSE','CLU_SIZE','PATRISTIC_DISTANCE'))), tmp, by=c('IDCLU','ID_CLOSE'))
	tmp		<- tmp[, list(PATRISTIC_DISTANCE=min(PATRISTIC_DISTANCE)), by=c('IDCLU','ID_CLOSE','CLU_SIZE')]
	tmp		<- merge( unique(subset(rtc, select=c('IDCLU','PTY_RUN','PTY_SIZE'))), tmp, by='IDCLU' )
	setkey(tmp, PTY_RUN, IDCLU, PATRISTIC_DISTANCE)	
	tmp		<- merge(tmp, tmp[, list(NCLU=length(unique(IDCLU))), by='PTY_RUN'], by='PTY_RUN')	
	tmp		<- tmp[, {
				z<- seq_len( max(0,ceiling((sn-PTY_SIZE[1])/NCLU[1])) )
				list(ID_CLOSE=ID_CLOSE[z], PATRISTIC_DISTANCE=PATRISTIC_DISTANCE[z])
			}, by=c('PTY_RUN','IDCLU','CLU_SIZE')]
	setnames(tmp, 'ID_CLOSE', 'ID')
	tmp[, ID_TYPE:='background']
	rtc[, ID_TYPE:='target']
	
	rtc[, PTY_SIZE:=NULL]
	rtc		<- rbind(rtc, tmp, fill=TRUE, use.names=TRUE)
	rtc		<- merge(rtc, rtc[, list(PTY_SIZE=length(ID)), by='PTY_RUN'], by='PTY_RUN')
	tmp		<- unique(subset(rtc, select=c(PTY_RUN, PTY_SIZE)))
	tmp[, PTY_RUN2:= seq_len(nrow(tmp))]
	rtc		<- merge(rtc, subset(tmp, select=c('PTY_RUN','PTY_RUN2')), by='PTY_RUN' )
	rtc[, PTY_RUN:=NULL]
	setnames(rtc, 'PTY_RUN2', 'PTY_RUN')
	
	#	in total 2320 individuals
	
	#
	#	set up new phyloscanner runs object
	#	find Sanger IDs in previous pty.runs
	#
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75.rda')
	tmp		<- copy(pty.runs)	
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75_part2.rda')
	tmp		<- rbind(tmp, pty.runs)
	set(tmp, NULL, c('BATCH','PTY_RUN'), NULL)
	tmp		<- unique(tmp)
	pty.runs<- NULL
	setnames(rtc, c('IDCLU','CLU_SIZE','ID','ID_TYPE','PATRISTIC_DISTANCE','PTY_SIZE'), c('TRCLU_ID','TRCLU_SIZE','RID','RID_TYPE','RID_PATRISTIC_DISTANCE','PTY_SIZE_RID'))
	pty.runs<- merge(rtc, tmp, by='RID')
	setkey(pty.runs, PTY_RUN, RID_TYPE, RID)
	pty.runs<- unique(pty.runs, by=c('RID','SID','PTY_RUN'))	
	pty.runs<- merge(pty.runs, pty.runs[, list(PTY_SIZE_SID=length(SID)), by='PTY_RUN'], by='PTY_RUN')
	pty.runs[, RENAME_SID:=NULL]
	pty.runs<- merge(pty.runs, pty.runs[, list(SID=SID, RENAME_SID=paste0(RID,'_fq',seq_along(SID))), by=c('PTY_RUN','RID')], by=c('PTY_RUN','RID','SID'))
	set(pty.runs, NULL, 'PTY_RUN', pty.runs[, 185L+PTY_RUN])
	
	save(pty.runs, file= '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170704_stagethree.rda')
	
}


project.RakaiAll.setup.phyloscanner.170704.stagetwo	<- function()
{
	require(igraph)
	
	infile		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170421.rda'
	infile		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170428_cl3.rda'
	infile		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170704_cl35.rda'
	load(infile)	
	#rtp			<- subset(rtp, GROUP=='TYPE_PAIR_DI' & POSTERIOR_SCORE>0.5)
	rtp			<- subset(rtp, GROUP=='TYPE_PAIR_DI2' & POSTERIOR_SCORE>=0.6)
	
	#	get transmission chains with igraph
	tmp			<- subset(rtp, select=c(ID1, ID2))			
	tmp			<- graph.data.frame(tmp, directed=FALSE, vertices=NULL)
	rtc			<- data.table(ID=V(tmp)$name, CLU=clusters(tmp, mode="weak")$membership)	
	tmp2		<- rtc[, list(CLU_SIZE=length(ID)), by='CLU']
	setkey(tmp2, CLU_SIZE)
	tmp2[, IDCLU:=rev(seq_len(nrow(tmp2)))]
	rtc			<- subset( merge(rtc, tmp2, by='CLU'))
	rtc[, CLU:=NULL]
	setkey(rtc, IDCLU)
	
	#tmp			<- subset(rtc, IDCLU==1)[, ID]
	#tmp			<- subset(rtp, ID1%in%tmp | ID2%in%tmp, select=c(ID1, ID2))
	#tmp			<- graph.data.frame(tmp, directed=FALSE, vertices=NULL)
	#par(mai=c(0,0,1,0))
	#plot(tmp, vertex.size=2, vertex.label.cex=0.25, edge.arrow.size=0.5, layout=layout.fruchterman.reingold(tmp, niter=1e3) )	
	
	#	add in couples that are not a close pair
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75.rda')
	tmp		<- copy(pty.runs)	
	#load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75_part2.rda')
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170704_b75_part2.rda')
	tmp		<- rbind(tmp, pty.runs)
	tmp		<- sort(unique(as.character(tmp[['RID']])))	
	#load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_info.rda")	
	rp		<- unique(subset(rp, FEMALE_RID%in%tmp & MALE_RID%in%tmp, select=c(FEMALE_RID,MALE_RID)))
	setnames(rp, 'FEMALE_RID', 'ID')
	rp		<- merge(rp, subset(rtc, select=c(ID, IDCLU)), all.x=1, by='ID')
	setnames(rp, c('ID','IDCLU','MALE_RID'), c('FEMALE_RID', 'FEMALE_IDCLU','ID'))
	rp		<- merge(rp, subset(rtc, select=c(ID, IDCLU)), all.x=1, by='ID')
	setnames(rp, c('ID','IDCLU'), c('MALE_RID', 'MALE_IDCLU'))
	tmp		<- rp[, which(is.na(FEMALE_IDCLU) & !is.na(MALE_IDCLU))]
	set(rp, tmp, 'FEMALE_IDCLU', rp[tmp, MALE_IDCLU])
	tmp		<- rp[, which(!is.na(FEMALE_IDCLU) & is.na(MALE_IDCLU))]
	set(rp, tmp, 'MALE_IDCLU', rp[tmp, FEMALE_IDCLU])
	tmp		<- rp[, which(is.na(FEMALE_IDCLU) & is.na(MALE_IDCLU))]
	set(rp, tmp, c('FEMALE_IDCLU','MALE_IDCLU'), rtc[, max(IDCLU)]+seq_along(tmp))
	setnames(rp, c('MALE_IDCLU'), c('IDCLU'))
	rp		<- subset(melt(rp, id.vars=c('IDCLU'), measure.vars=c('MALE_RID','FEMALE_RID'), value.name='ID'), select=c(ID, IDCLU))
	set(rtc, NULL, 'CLU_SIZE', NULL)
	rtc		<- unique(rbind(rtc, rp))
	tmp		<- rtc[, list(CLU_SIZE=length(ID)), by='IDCLU']
	setkey(tmp, CLU_SIZE)
	tmp[, IDCLU2:=rev(seq_len(nrow(tmp)))]	
	rtc		<- merge(rtc, tmp, by='IDCLU')
	rtc[, IDCLU:=NULL]
	setnames(rtc, 'IDCLU2', 'IDCLU')
	
	#	1547 close individuals 
	#	added 493 unrelated individuals that are in couple
	
	#	merge up to 8 individuals into the same phyloscanner run
	tn			<- 8
	tmp			<- unique(subset(rtc, select=c(IDCLU, CLU_SIZE)))		
	tmp[, PTY_RUN:= IDCLU]
	tmp[, PTY_SIZE:= CLU_SIZE]
	for(i in rev(seq_len(nrow(tmp))))
	{
		if(tmp[i,PTY_SIZE<tn])
		{
			for(j in 1:i)
			{
				if( (tmp[j,PTY_SIZE]+tmp[i,PTY_SIZE])<=tn)
				{
					set(tmp, i, 'PTY_RUN', tmp[j, PTY_RUN])					
					set(tmp, j, 'PTY_SIZE',  tmp[i, PTY_SIZE]+ tmp[j, PTY_SIZE])
					set(tmp, i, 'PTY_SIZE', 0L)
					break
				}
			}
		}
	}
	# manually add IDCLUD 591 to PTY_RUN 592	
	#set(tmp, tmp[, which(PTY_RUN==591)], 'PTY_RUN', 592L)	
	tmp			<- tmp[, list(IDCLU=IDCLU, CLU_SIZE=CLU_SIZE, PTY_SIZE=sum(CLU_SIZE)), by=c('PTY_RUN')]
	rtc			<- merge(rtc, tmp, by=c('IDCLU','CLU_SIZE'))
	#	for every individual in a cluster
	#	find 15 closest others
	sn		<- 15
	setkey(dmin, ID1, ID2)	
	dcl		<- rtc[, {
				tmp		<- subset(dmin, ID1==ID)
				tmp2	<- subset(dmin, ID2==ID)
				setnames(tmp2, c('ID1','ID2'), c('ID2','ID1'))
				tmp		<- rbind(tmp, tmp2)
				tmp		<- tmp[, list(PATRISTIC_DISTANCE=mean(PATRISTIC_DISTANCE_MIN)), by=c('ID1','ID2')]
				setkey(tmp, PATRISTIC_DISTANCE)
				tmp		<- tmp[seq_len(sn),]
				list(ID_CLOSE= tmp[,ID2], PATRISTIC_DISTANCE=tmp[, PATRISTIC_DISTANCE])				
			}, by=c('ID','CLU_SIZE','IDCLU')]
	#	find as many close other individuals as needed to fill the phyloscanner run
	#	up to 15 individuals
	dcl		<- subset(dcl, !is.na(ID_CLOSE))
	tmp		<- dcl[, list(ID_CLOSE=setdiff(ID_CLOSE, ID)), by='IDCLU']
	tmp		<- merge(unique(subset(dcl, select=c('IDCLU','ID_CLOSE','CLU_SIZE','PATRISTIC_DISTANCE'))), tmp, by=c('IDCLU','ID_CLOSE'))
	tmp		<- tmp[, list(PATRISTIC_DISTANCE=min(PATRISTIC_DISTANCE)), by=c('IDCLU','ID_CLOSE','CLU_SIZE')]
	tmp		<- merge( unique(subset(rtc, select=c('IDCLU','PTY_RUN','PTY_SIZE'))), tmp, by='IDCLU' )
	setkey(tmp, PTY_RUN, IDCLU, PATRISTIC_DISTANCE)	
	tmp		<- merge(tmp, tmp[, list(NCLU=length(unique(IDCLU))), by='PTY_RUN'], by='PTY_RUN')	
	tmp		<- tmp[, {
				z<- seq_len( max(0,ceiling((sn-PTY_SIZE[1])/NCLU[1])) )
				list(ID_CLOSE=ID_CLOSE[z], PATRISTIC_DISTANCE=PATRISTIC_DISTANCE[z])
			}, by=c('PTY_RUN','IDCLU','CLU_SIZE')]
	setnames(tmp, 'ID_CLOSE', 'ID')
	tmp[, ID_TYPE:='background']
	rtc[, ID_TYPE:='target']
	
	rtc[, PTY_SIZE:=NULL]
	rtc		<- rbind(rtc, tmp, fill=TRUE, use.names=TRUE)
	rtc		<- merge(rtc, rtc[, list(PTY_SIZE=length(ID)), by='PTY_RUN'], by='PTY_RUN')
	tmp		<- unique(subset(rtc, select=c(PTY_RUN, PTY_SIZE)))
	tmp[, PTY_RUN2:= seq_len(nrow(tmp))]
	rtc		<- merge(rtc, subset(tmp, select=c('PTY_RUN','PTY_RUN2')), by='PTY_RUN' )
	rtc[, PTY_RUN:=NULL]
	setnames(rtc, 'PTY_RUN2', 'PTY_RUN')
	setnames(rtc, c('IDCLU','CLU_SIZE','ID','ID_TYPE','PATRISTIC_DISTANCE','PTY_SIZE'), c('TRCLU_ID','TRCLU_SIZE','RID','RID_TYPE','RID_PATRISTIC_DISTANCE','PTY_SIZE_RID'))
	#
	#	set up new phyloscanner runs object
	#	find Sanger IDs
	#
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170704_assemblystatus.rda')
	dc		<- unique(subset(dc, !is.na(SID), select=c(RID,SID,PID)))
	pty.runs<- merge(rtc, dc, by='RID')
	setkey(pty.runs, PTY_RUN, RID_TYPE, RID)
	pty.runs<- unique(pty.runs, by=c('RID','SID','PTY_RUN'))	
	pty.runs<- merge(pty.runs, pty.runs[, list(PTY_SIZE_SID=length(SID)), by='PTY_RUN'], by='PTY_RUN')	
	#unique(subset(pty.runs, select=c(PTY_RUN,PTY_SIZE_SID)))[, table(PTY_SIZE_SID)]
	save(pty.runs, file= '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170704_stagetwo.rda')
}

project.RakaiAll.setup.phyloscanner.170301.stagethree	<- function()
{
	require(igraph)
	
	#	load dmin for filling close pairs
	infile		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170428_cl3.rda'
	load(infile)	
	#	look at run 1
	infile		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun/ptyr1_pairwise_relationships.rda'
	load(infile)	
	group		<- 'TYPE_PAIR_TO'
	state.ml	<- 'ancestral/\nintermingled'
	
	rtp			<- subset(rplkl, GROUP==group)
	#tmp			<- rtp[, which(TYPE==state.ml)]
	#set(rtp, tmp, 'KEFF', rtp[tmp, KEFF*1.25])
	#	determine most likely relationship type
	rtp			<- rtp[, 	{
								z<- which.max(KEFF)
								list(N=N[1], NEFF=NEFF[1], TYPE=TYPE[z], K=K[z], KEFF=KEFF[z])
							}, by=c('ID2','ID1')]
	rtp			<- subset(rtp, TYPE==state.ml)
	
	#	get transmission chains with igraph
	tmp			<- subset(rtp, select=c(ID1, ID2))			
	tmp			<- graph.data.frame(tmp, directed=FALSE, vertices=NULL)
	rtc			<- data.table(ID=V(tmp)$name, CLU=clusters(tmp, mode="weak")$membership)	
	tmp2		<- rtc[, list(CLU_SIZE=length(ID)), by='CLU']
	setkey(tmp2, CLU_SIZE)
	tmp2[, IDCLU:=rev(seq_len(nrow(tmp2)))]
	rtc			<- subset( merge(rtc, tmp2, by='CLU'))
	rtc[, CLU:=NULL]
	setkey(rtc, IDCLU)
		
	#	add in couples that were in run1 and are now disconnected
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")	
	rp		<- unique(subset(rp, select=c(FEMALE_RID,MALE_RID)))
	setnames(rp, c('FEMALE_RID', 'MALE_RID'), c('ID2', 'ID1'))
	tmp		<- copy(rp)
	setnames(tmp, c('ID2', 'ID1'), c('ID1', 'ID2'))
	rp		<- rbind(rp, tmp)
	rp		<- merge(rp, unique(subset(rplkl, select=c(ID1,ID2))), by=c('ID1', 'ID2'))
	setnames(rp, 'ID1', 'ID')
	rp		<- merge(rp, subset(rtc, select=c(ID, IDCLU)), all.x=1, by='ID')
	setnames(rp, c('ID','IDCLU','ID2'), c('ID1', 'ID1_IDCLU','ID'))
	rp		<- merge(rp, subset(rtc, select=c(ID, IDCLU)), all.x=1, by='ID')
	setnames(rp, c('ID','IDCLU'), c('ID2', 'ID2_IDCLU'))
	tmp		<- rp[, which(is.na(ID2_IDCLU) & !is.na(ID1_IDCLU))]
	set(rp, tmp, 'ID2_IDCLU', rp[tmp, ID1_IDCLU])
	tmp		<- rp[, which(is.na(ID1_IDCLU) & !is.na(ID2_IDCLU))]
	set(rp, tmp, 'ID1_IDCLU', rp[tmp, ID2_IDCLU])	
	tmp		<- rp[, which(is.na(ID1_IDCLU) & is.na(ID2_IDCLU))]
	set(rp, tmp, c('ID1_IDCLU','ID2_IDCLU'), rtc[, max(IDCLU)]+seq_along(tmp))	
	setnames(rp, c('ID1_IDCLU'), c('IDCLU'))
	rp		<- subset(melt(rp, id.vars=c('IDCLU'), measure.vars=c('ID1','ID2'), value.name='ID'), select=c(ID, IDCLU))
	set(rtc, NULL, 'CLU_SIZE', NULL)
	rtc		<- unique(rbind(rtc, rp))
	tmp		<- rtc[, list(CLU_SIZE=length(ID)), by='IDCLU']
	setkey(tmp, CLU_SIZE)
	tmp[, IDCLU2:=rev(seq_len(nrow(tmp)))]	
	rtc		<- merge(rtc, tmp, by='IDCLU')
	rtc[, IDCLU:=NULL]
	setnames(rtc, 'IDCLU2', 'IDCLU')
	
	#	85 individuals that are most likely ancestral/intermingled 
	#	added 6 couples back in that are not most likely ancestral/intermingled 
	
	#	merge up to 8 individuals into the same phyloscanner run
	tn			<- 9
	tmp			<- unique(subset(rtc, select=c(IDCLU, CLU_SIZE)))		
	tmp[, PTY_RUN:= IDCLU]
	tmp[, PTY_SIZE:= CLU_SIZE]
	setkey(tmp, PTY_RUN)
	for(i in rev(seq_len(nrow(tmp))))
	{
		if(tmp[i,PTY_SIZE<tn])
		{
			for(j in 1:i)
			{
				if( (tmp[j,PTY_SIZE]+tmp[i,PTY_SIZE])<=tn)
				{
					set(tmp, i, 'PTY_RUN', tmp[j, PTY_RUN])					
					set(tmp, j, 'PTY_SIZE',  tmp[i, PTY_SIZE]+ tmp[j, PTY_SIZE])
					set(tmp, i, 'PTY_SIZE', 0L)
					break
				}
			}
		}
	}
	# manually add IDCLUD 9 to PTY_RUN 8	
	set(tmp, tmp[, which(PTY_RUN==9)], 'PTY_RUN', 8L)	
	tmp			<- tmp[, list(IDCLU=IDCLU, CLU_SIZE=CLU_SIZE, PTY_SIZE=sum(CLU_SIZE)), by=c('PTY_RUN')]
	rtc			<- merge(rtc, tmp, by=c('IDCLU','CLU_SIZE'))
	#	for every individual in a cluster
	#	find 15 closest others
	sn		<- 15
	setkey(dmin, ID1, ID2)	
	dcl		<- rtc[, {
				tmp		<- subset(dmin, ID1==ID)
				tmp2	<- subset(dmin, ID2==ID)
				setnames(tmp2, c('ID1','ID2'), c('ID2','ID1'))
				tmp		<- rbind(tmp, tmp2)
				tmp		<- tmp[, list(PATRISTIC_DISTANCE=mean(PATRISTIC_DISTANCE_MIN)), by=c('ID1','ID2')]
				setkey(tmp, PATRISTIC_DISTANCE)
				tmp		<- tmp[seq_len(sn),]
				list(ID_CLOSE= tmp[,ID2], PATRISTIC_DISTANCE=tmp[, PATRISTIC_DISTANCE])				
			}, by=c('ID','CLU_SIZE','IDCLU')]
	#	find as many close other individuals as needed to fill the phyloscanner run
	#	up to 15 individuals
	dcl		<- subset(dcl, !is.na(ID_CLOSE))
	tmp		<- dcl[, list(ID_CLOSE=setdiff(ID_CLOSE, ID)), by='IDCLU']
	tmp		<- merge(unique(subset(dcl, select=c('IDCLU','ID_CLOSE','CLU_SIZE','PATRISTIC_DISTANCE'))), tmp, by=c('IDCLU','ID_CLOSE'))
	tmp		<- tmp[, list(PATRISTIC_DISTANCE=min(PATRISTIC_DISTANCE)), by=c('IDCLU','ID_CLOSE','CLU_SIZE')]
	tmp		<- merge( unique(subset(rtc, select=c('IDCLU','PTY_RUN','PTY_SIZE'))), tmp, by='IDCLU' )
	setkey(tmp, PTY_RUN, IDCLU, PATRISTIC_DISTANCE)	
	tmp		<- merge(tmp, tmp[, list(NCLU=length(unique(IDCLU))), by='PTY_RUN'], by='PTY_RUN')	
	tmp		<- tmp[, {
				z<- seq_len( max(0,ceiling((sn-PTY_SIZE[1])/NCLU[1])) )
				list(ID_CLOSE=ID_CLOSE[z], PATRISTIC_DISTANCE=PATRISTIC_DISTANCE[z])
			}, by=c('PTY_RUN','IDCLU','CLU_SIZE')]
	setnames(tmp, 'ID_CLOSE', 'ID')
	tmp[, ID_TYPE:='background']
	rtc[, ID_TYPE:='target']
	
	rtc[, PTY_SIZE:=NULL]
	rtc		<- rbind(rtc, tmp, fill=TRUE, use.names=TRUE)
	rtc		<- merge(rtc, rtc[, list(PTY_SIZE=length(ID)), by='PTY_RUN'], by='PTY_RUN')
	tmp		<- unique(subset(rtc, select=c(PTY_RUN, PTY_SIZE)))
	tmp[, PTY_RUN2:= seq_len(nrow(tmp))]
	rtc		<- merge(rtc, subset(tmp, select=c('PTY_RUN','PTY_RUN2')), by='PTY_RUN' )
	rtc[, PTY_RUN:=NULL]
	setnames(rtc, 'PTY_RUN2', 'PTY_RUN')
	
	#	in total 143 individuals
	
	#
	#	set up new phyloscanner runs object
	#	find Sanger IDs in previous pty.runs
	#
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75.rda')
	tmp		<- copy(pty.runs)	
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75_part2.rda')
	tmp		<- rbind(tmp, pty.runs)
	set(tmp, NULL, c('BATCH','PTY_RUN'), NULL)
	tmp		<- unique(tmp)
	pty.runs<- NULL
	setnames(rtc, c('IDCLU','CLU_SIZE','ID','ID_TYPE','PATRISTIC_DISTANCE','PTY_SIZE'), c('TRCLU_ID','TRCLU_SIZE','RID','RID_TYPE','RID_PATRISTIC_DISTANCE','PTY_SIZE_RID'))
	pty.runs<- merge(rtc, tmp, by='RID')
	setkey(pty.runs, PTY_RUN, RID_TYPE, RID)
	pty.runs<- unique(pty.runs, by=c('RID','SID','PTY_RUN'))	
	pty.runs<- merge(pty.runs, pty.runs[, list(PTY_SIZE_SID=length(SID)), by='PTY_RUN'], by='PTY_RUN')
	
	#
	#	add 240 so we continue previous enumeration
	set(pty.runs, NULL, 'PTY_RUN', pty.runs[,PTY_RUN+240L])
	save(pty.runs, file= '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_stagethree.rda')
}

######################################################################################
project.Rakai.ExaMLTree.170601<- function()
{	
	require(big.phylo)
	if(0)	#	MAC
	{
		indir		<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_sequences_latest'
		infile		<- "ATHENA_1610_Sequences_LANL_codonaligned_noDRM_subtype_B_wOutgroup.fasta"
		outdir		<- '~/Dropbox (SPH Imperial College)/2016_ATHENA_Oct_Update/processed_trees'
	}
	if(1)	#	HPC
	{
		indir		<- file.path(HOME,'alignments_170603')
		infile		<- 'PANGEA_HIV_n4232_Imperial_v170601_UgandaAlignment_minCov700.fasta'
		outdir		<- indir
	}
	
	args.parser	<- paste("-m DNA")
	#	ExaML bootstrap args
	bs.from		<- 0
	bs.to		<- 0
	bs.n		<- 100
	
	args.examl	<- "-f d -D -m GAMMA"	#	 -- this is the default that worked in 24 hours	
	cmd			<- cmd.examl.bootstrap(indir, infile, bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, opt.bootstrap.by="codon", args.parser=args.parser, args.examl=args.examl, tmpdir.prefix="examl")					
	dummy		<- lapply(cmd, function(x)
			{				
				#x		<- cmd.hpcwrapper(x, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
				#x		<- cmd.hpcwrapper(x, hpc.walltime=71, hpc.q="pqeph", hpc.mem="3600mb", hpc.nproc=1)
				#x		<- cmd.hpcwrapper(x, hpc.walltime=79, hpc.q="pqeelab", hpc.mem="5800mb", hpc.nproc=8)
				cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=79, hpc.q="pqeelab", hpc.mem="5800mb",  hpc.nproc=1)
				cmd		<- paste(cmd,x,sep='\n')					
				signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
				outfile	<- paste("exa",signat,sep='.')
				cat(cmd)
				cmd.hpccaller(outdir, outfile, cmd)
				Sys.sleep(1)
			})	
}
######################################################################################
project.Rakai.FastTree.170628<- function()
{
	require(ape)
	
	sq				<- read.dna("~/Dropbox (SPH Imperial College)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.fasta", format='fasta')
	#	handle repeats
	sqi				<- data.table(TAXA=rownames(sq), DUMMY=seq_len(nrow(sq)))
	tmp				<- sqi[, which(duplicated(TAXA))]
	cat('Repeats', paste(sqi[tmp,TAXA], collapse=', '))
	set(sqi, tmp, 'TAXA', sqi[tmp, paste(TAXA,'-R2',sep='')])
	setkey(sqi, DUMMY)
	rownames(sq)	<- sqi[,TAXA]
	#	get coverage
	tmp				<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
	sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]	
	sqi[, PNG:= sqi[, factor(grepl('PG',TAXA),levels=c(TRUE,FALSE),labels=c('Y','N'))]]
	sqi[, PIDF:= gsub('-R[0-9]+','',TAXA)]
	#	for each PIDF keep taxon with highest coverage
	sqi		<- merge(sqi, sqi[, list(DUMMY=DUMMY[which.max(COV)]), by='PIDF'], by=c('PIDF','DUMMY'))
	#	select
	sq				<- sq[ sqi$TAXA, ]
	rownames(sq)	<- gsub('-R[0-9]+','',rownames(sq))
	#	add new consensuses
	infile	<- '/Users/Oliver/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017/Rakai_Mar2017_GlobalAlignment_10040_BestSeqPerSample.fasta'
	sn		<- read.dna(infile, format='fa')
	sq		<- rbind(sq,sn)
	#	handle repeats
	sqi				<- data.table(TAXA=rownames(sq), DUMMY=seq_len(nrow(sq)))
	tmp				<- sqi[, which(duplicated(TAXA))]
	cat('Repeats', paste(sqi[tmp,TAXA], collapse=', '))
	set(sqi, tmp, 'TAXA', sqi[tmp, paste(TAXA,'-R2',sep='')])
	setkey(sqi, DUMMY)
	rownames(sq)	<- sqi[,TAXA]
	#	get coverage
	tmp				<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
	sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]	
	sqi[, PNG:= sqi[, factor(grepl('PG',TAXA),levels=c(TRUE,FALSE),labels=c('Y','N'))]]
	sqi[, PIDF:= gsub('-R[0-9]+','',TAXA)]
	#	for each PIDF keep taxon with highest coverage
	sqi				<- merge(sqi, sqi[, list(DUMMY=DUMMY[which.max(COV)]), by='PIDF'], by=c('PIDF','DUMMY'))
	#	select
	sq				<- sq[ sqi$TAXA, ]
	rownames(sq)	<- gsub('-R[0-9]+','',rownames(sq))
	#
	#	update stragglers
	infile			<- '~/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017/stragglers_GlobalAlignment_10040.fasta'
	sn				<- read.dna(infile, format='fa')
	rownames(sn)	<- gsub('_consensus_MinCov_10_50','',rownames(sn))	
	#	keep only longest straggler for each PANGEA ID
	stg				<- data.table(SID=rownames(sn), DUMMY=seq_len(nrow(sn)))
	tmp				<- sapply(seq_len(nrow(sn)), function(i) base.freq(sn[i,], all=TRUE, freq=TRUE))
	stg[, COV:=ncol(sn)-apply( tmp[c('-','?'),], 2, sum	)]		
	load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_assemblystatus.rda')
	tmp				<- merge(sqi, subset(dc, !is.na(SID), select=c(PIDF, SID)), by='PIDF',all.x=TRUE)
	stg				<- merge(stg, subset(tmp, select=c(SID, PIDF)), by='SID',all.x=1)
	stg				<- merge(stg, stg[, list(DUMMY=DUMMY[which.max(COV)]), by='PIDF'], by=c('PIDF','DUMMY'))
	#	add to re-processed stragglers those taxa that are now junk
	tmp				<- data.table(PIDF=c('PG14-UG501871-S03875','PG14-UG500718-S02722','PG14-UG503751-S05755','PG14-UG501956-S03960','PG14-UG503966-S05970','PG14-UG503944-S05948','PG14-UG503322-S05326','PG14-UG501838-S03842','PG14-UG503896-S05900','PG14-UG500743-S02747','PG14-UG503259-S05263','PG14-UG502371-S04375','PG14-UG503509-S05513','PG14-UG503648-S05652','PG14-UG503750-S05754','PG14-UG503518-S05522','PG14-UG501117-S03121','PG14-UG501506-S03510','PG14-UG503995-S05999','PG14-UG502131-S04135','PG14-UG501758-S03762','PG14-UG503968-S05972','PG14-UG503565-S05569','PG14-UG502024-S04028','PG14-UG503929-S05933','PG14-UG503673-S05677','PG14-UG502418-S04422','PG14-UG503643-S05647','PG14-UG501140-S03144','PG14-UG503190-S05194','PG14-UG500711-S02715','PG14-UG503891-S05895','PG14-UG503697-S05701','PG14-UG501743-S03747','PG14-UG503596-S05600','PG14-UG502370-S04374','PG14-UG501271-S03275','PG14-UG503555-S05559','PG14-UG502472-S04476','PG14-UG503500-S05504'))
	stg				<- merge(tmp, stg, by='PIDF', all.x=1)
	#write.csv(stg, file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/stragglers_no1.csv')
	#	reduce central alignment to taxa not in stg	
	sq				<- sq[setdiff(rownames(sq),stg[, PIDF]),]
	tmp				<- subset(stg, !is.na(SID))
	sn				<- sn[tmp[,SID],]
	rownames(sn)	<- tmp[,PIDF]
	sq				<- rbind(sq,sn)	
	sqi				<- merge(sqi,data.table(PIDF=rownames(sq)), by='PIDF')	
	sqi[, COV:=NULL]
	tmp				<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
	sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]	
	#
	#	add in new stragglers from 170628
	#
	infile			<- '/Users/Oliver/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017/BlastfiltTest_10_50_ForGlobalAln_10040.fasta'
	sn				<- read.dna(infile, format='fa')
	rownames(sn)	<- gsub('_consensus_MinCov_10_50','',rownames(sn))	
	#	keep only longest straggler for each PANGEA ID
	stg				<- data.table(SID=rownames(sn), DUMMY=seq_len(nrow(sn)))
	tmp				<- sapply(seq_len(nrow(sn)), function(i) base.freq(sn[i,], all=TRUE, freq=TRUE))
	stg[, COV:=ncol(sn)-apply( tmp[c('-','?'),], 2, sum	)]		
	load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_assemblystatus.rda')	
	stg				<- merge(stg, subset(dc, !is.na(SID), select=c(PIDF, SID)), by='SID',all.x=TRUE)	
	stg				<- merge(stg, stg[, list(SID=SID, PIDF2=paste0('TANYA0628-',PIDF,'-R',seq_along(SID))), by='PIDF'], by=c('SID','PIDF'))
	sn				<- sn[stg[,SID],]
	rownames(sn)	<- stg[,PIDF2]	
	sq				<- rbind(sq,sn)		
	#
	# update sqi
	sqi				<- data.table(TAXA=rownames(sq), DUMMY=seq_len(nrow(sq)))	
	tmp				<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
	sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]	
	sqi[, PNG:= sqi[, factor(grepl('PG',TAXA),levels=c(TRUE,FALSE),labels=c('Y','N'))]]
	sqi[, PIDF:= gsub('TANYA0628-','',gsub('-R[0-9]+','',TAXA))]
	sqi[, GENBANK:= NA_character_]
	tmp				<- sqi[, which(PNG=='N')]
	set(sqi,tmp, 'GENBANK', sqi[tmp,gsub('.*\\.([^\\.]+)$','\\1',PIDF)])
	set(sqi,sqi[, which(PNG=='Y')],'GENBANK','PANGEA')	
	#
	#	save UG
	tmp				<- read.dna('/Users/Oliver/git/Phyloscanner.R.utilities/inst/HIV1_compendium_AD_B_CPX_v2.fasta', format='fa')
	tmp				<- rownames(tmp)
	tmp				<- gsub('.*_([^_]+)$','\\1',tmp[!grepl('CPX',tmp)])
	sqi				<- merge(sqi,data.table(GENBANK=c(tmp,'PANGEA')),by='GENBANK')		
	sq				<- sq[sqi$TAXA,]
	sqi				<- subset(sqi, PNG=='N' | (PNG=='Y' & grepl('PG[0-9][0-9]-UG.*',TAXA)))
	sq				<- sq[sqi$TAXA,]	
	tmp				<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017','PANGEA_HIV_n5046_Imperial_v170628_UgandaAlignment.fasta')
	write.dna(sq, file=tmp, format='fa')	
	
	#	tree as requested by Tanya
	tmp2			<- sq[subset(sqi, COV>700)[,TAXA],]
	tmp				<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017','PANGEA_HIV_n4216_Imperial_v170628_UgandaAlignment_minCov700.fasta')	
	write.dna(tmp2, file=tmp, format='fa')	
	infile.fasta	<- '/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_n4216_Imperial_v170628_UgandaAlignment_minCov700.fasta'
	outfile.ft		<- gsub('\\.fasta','_ft.newick',infile.fasta)
	tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma')
	cat(tmp)
	
	#	take out samples that are covered in updated stragglers
	tmp				<- stg[, unique(PIDF)]
	tmp3			<- tmp2[setdiff(rownames(tmp2), tmp),]
	tmp				<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017','PANGEA_HIV_n4201_Imperial_v170628_UgandaAlignment_minCov700.fasta')	
	write.dna(tmp3, file=tmp, format='fa')	
	infile.fasta	<- '/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_n4201_Imperial_v170628_UgandaAlignment_minCov700.fasta'
	outfile.ft		<- gsub('\\.fasta','_ft.newick',infile.fasta)
	tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma')
	cat(tmp)
	
	tmp				<- merge( subset(sqi, COV>700), unique(subset(stg, select=PIDF)), by='PIDF' )
	tmp				<- merge(tmp,tmp[, list(N=length(TAXA)),by='PIDF'],by='PIDF')	
	write.csv(subset(tmp, N>1), file='/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_n4216_Imperial_v170628_UgandaAlignment_minCov700_compareTANYAtaxa.csv')
	
	
	infile			<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_n4199_Imperial_v170601_UgandaAlignment_minCov700_ft.newick'
	infile			<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_n4216_Imperial_v170628_UgandaAlignment_minCov700_ft.newick'
	infile			<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_n4201_Imperial_v170628_UgandaAlignment_minCov700_ft.newick'
	ph				<- read.tree(infile)
	tmp				<- as.numeric(ph$node.label)
	tmp[is.na(tmp)]	<- 0
	ph$node.label	<- tmp
	tmp				<- which(ph$tip.label=='B.FR.83.HXB2_LAI_IIIB_BRU.K03455')
	ph				<- phytools::reroot(ph, tmp, position=ph$edge.length[which(ph$edge[, 2]==tmp)] / 2)
	ph				<- ladderize(ph)
	pdf(gsub('newick','pdf',infile), w=15, h=500)
	plot(ph, cex=0.5, show.node.label=TRUE)
	dev.off()
	
	
	tmp				<- extract.clade(ph,getMRCA(ph,c("PG14-UG502292-S04296","PG14-UG502476-S04480")))
	tmp				<- merge(sqi, data.table(TAXA=tmp$tip.label), by='TAXA')
	tmp				<- merge(tmp, dc, by='PIDF')
	write.table(tmp, row.names=FALSE, file='/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/stragglers_no3.csv')
	
}

######################################################################################
project.Rakai.FastTree.170630<- function()
{
	require(ape)	
	sq				<- read.dna("~/Dropbox (SPH Imperial College)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.fasta", format='fasta')
	#	handle repeats
	sqi				<- data.table(TAXA=rownames(sq), DUMMY=seq_len(nrow(sq)))
	tmp				<- sqi[, which(duplicated(TAXA))]
	cat('Repeats', paste(sqi[tmp,TAXA], collapse=', '))
	set(sqi, tmp, 'TAXA', sqi[tmp, paste(TAXA,'-R2',sep='')])
	setkey(sqi, DUMMY)
	rownames(sq)	<- sqi[,TAXA]
	#	get coverage
	tmp				<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
	sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]	
	sqi[, PNG:= sqi[, factor(grepl('PG',TAXA),levels=c(TRUE,FALSE),labels=c('Y','N'))]]
	sqi[, PIDF:= gsub('-R[0-9]+','',TAXA)]
	#	for each PIDF keep taxon with highest coverage
	sqi		<- merge(sqi, sqi[, list(DUMMY=DUMMY[which.max(COV)]), by='PIDF'], by=c('PIDF','DUMMY'))
	#	select
	sq				<- sq[ sqi$TAXA, ]
	rownames(sq)	<- gsub('-R[0-9]+','',rownames(sq))
	#	add new consensuses
	infile	<- '/Users/Oliver/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017/Rakai_Mar2017_GlobalAlignment_10040_BestSeqPerSample.fasta'
	sn		<- read.dna(infile, format='fa')
	sq		<- rbind(sq,sn)
	#	handle repeats
	sqi				<- data.table(TAXA=rownames(sq), DUMMY=seq_len(nrow(sq)))
	tmp				<- sqi[, which(duplicated(TAXA))]
	cat('Repeats', paste(sqi[tmp,TAXA], collapse=', '))
	set(sqi, tmp, 'TAXA', sqi[tmp, paste(TAXA,'-R2',sep='')])
	setkey(sqi, DUMMY)
	rownames(sq)	<- sqi[,TAXA]
	#	get coverage
	tmp				<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
	sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]	
	sqi[, PNG:= sqi[, factor(grepl('PG',TAXA),levels=c(TRUE,FALSE),labels=c('Y','N'))]]
	sqi[, PIDF:= gsub('-R[0-9]+','',TAXA)]
	#	for each PIDF keep taxon with highest coverage
	sqi				<- merge(sqi, sqi[, list(DUMMY=DUMMY[which.max(COV)]), by='PIDF'], by=c('PIDF','DUMMY'))
	#	select
	sq				<- sq[ sqi$TAXA, ]
	rownames(sq)	<- gsub('-R[0-9]+','',rownames(sq))
	
	#	add to sqi list of taxa that were re-processed
	load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_assemblystatus.rda')	
	stg	<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017/OddCladeCheckMiseq.txt', header=FALSE, col.names='SID',stringsAsFactors=FALSE))
	stg	<- merge(dc, stg, by='SID')
	#sn				<- read.dna(file='~/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017/OddCladeCheckMiseq_10_50_ForGlobalAln_10040.fasta',format='fasta')
	sn				<- read.dna(file='~/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017/OddCladeCheck_WS17_Miseq_10_50_ForGlobalAln_10040.fasta',format='fasta')
	rownames(sn)	<- gsub('_consensus_MinCov_10_50','',rownames(sn))
	tmp2	<- data.table(SID=rownames(sn), DUMMY=seq_len(nrow(sn)))
	tmp		<- sapply(seq_len(nrow(sn)), function(i) base.freq(sn[i,], all=TRUE, freq=TRUE))
	tmp2[, COV:=ncol(sn)-apply( tmp[c('-','?'),], 2, sum	)]		
	stg		<- merge(stg, tmp2, by='SID',all.x=1)
	set(stg, stg[, which(is.na(COV))],'COV',0)
	#	remove from sqi old taxa that do not appear in the new list
	tmp		<- stg[, unique(PIDF)]
	sqi[, STG0630:='n']
	set(sqi, sqi[, which(PIDF%in%tmp)],'STG0630','y')
	tmp		<- subset(stg, COV>0)[, unique(PIDF)]
	sqi[, STG0630_COV:='n']
	set(sqi, sqi[, which(PIDF%in%tmp)],'STG0630_COV','y')
	sqi		<- subset(sqi, !(STG0630=='y' & STG0630_COV=='n'))
	#	change names
	tmp		<- sqi[, which(PNG=='Y' & STG0630=='n')]
	set(sqi, tmp, 'TAXA', sqi[tmp, paste0('NOCHANGE-',TAXA)])
	tmp		<- sqi[, which(PNG=='Y' & STG0630=='y')]
	set(sqi, tmp, 'TAXA', sqi[tmp, paste0('PREVIOUS-',TAXA)])
	set(sqi, NULL, c('STG0630','STG0630_COV','DUMMY'), NULL)
	sq				<- sq[ sqi$PIDF, ]
	rownames(sq)	<- sqi$TAXA
	stg		<- subset(stg, COV>0, c(SID,PIDF,COV))
	stg[, TAXA:= paste0('TANYA0630-',PIDF,'-',SID)]
	stg[, PNG:='Y']
	sn				<- sn[ stg$SID, ]
	rownames(sn)	<- stg$TAXA
	#	put all sequences together
	sqi				<- rbind(sqi, stg, fill=TRUE)
	sq				<- rbind(sq, sn)
	#	add GENBANK
	sqi[, GENBANK:= NA_character_]
	tmp				<- sqi[, which(PNG=='N')]
	set(sqi,tmp, 'GENBANK', sqi[tmp,gsub('.*\\.([^\\.]+)$','\\1',PIDF)])
	set(sqi,sqi[, which(PNG=='Y')],'GENBANK','PANGEA')	
	
	#	save UG
	tmp				<- read.dna('/Users/Oliver/git/Phyloscanner.R.utilities/inst/HIV1_compendium_AD_B_CPX_v2.fasta', format='fa')
	tmp				<- rownames(tmp)
	tmp				<- gsub('.*_([^_]+)$','\\1',tmp[!grepl('CPX',tmp)])
	sqi				<- merge(sqi,data.table(GENBANK=c(tmp,'PANGEA')),by='GENBANK')		
	sqi				<- subset(sqi, PNG=='N' | (PNG=='Y' & grepl('PG[0-9][0-9]-UG.*',TAXA)))
	sq				<- sq[sqi$TAXA,]	
	tmp				<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017','PANGEA_HIV_n5194_Imperial_v170630b_UgandaAlignment.fasta')
	write.dna(sq, file=tmp, format='fa')		
	#	tree as requested by Tanya
	tmp2			<- sq[subset(sqi, COV>700)[,TAXA],]
	tmp				<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017','PANGEA_HIV_n4288_Imperial_v170630b_UgandaAlignment_minCov700.fasta')	
	write.dna(tmp2, file=tmp, format='fa')	
	#	take out samples that are covered in updated stragglers
	tmp3			<- sq[subset(sqi, COV>700 & !grepl('PREVIOUS',TAXA))[,TAXA],]
	tmp				<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017','PANGEA_HIV_n4171_Imperial_v170630b_UgandaAlignment_minCov700.fasta')	
	write.dna(tmp3, file=tmp, format='fa')	
	
	#
	infile.fasta	<- '/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_n4171_Imperial_v170630b_UgandaAlignment_minCov700.fasta'
	#infile.fasta	<- '/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_n4288_Imperial_v170630_UgandaAlignment_minCov700.fasta'
	outfile.ft		<- gsub('\\.fasta','_ft.newick',infile.fasta)
	tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma')
	cat(tmp)
	
	
	
	infile			<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_n4199_Imperial_v170601_UgandaAlignment_minCov700_ft.newick'
	infile			<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_n4216_Imperial_v170628_UgandaAlignment_minCov700_ft.newick'
	infile			<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_n4201_Imperial_v170628_UgandaAlignment_minCov700_ft.newick'
	infile			<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_n4171_Imperial_v170630b_UgandaAlignment_minCov700_ft.newick'
	ph				<- read.tree(infile)
	tmp				<- as.numeric(ph$node.label)
	tmp[is.na(tmp)]	<- 0
	ph$node.label	<- tmp
	tmp				<- which(ph$tip.label=='B.FR.83.HXB2_LAI_IIIB_BRU.K03455')
	ph				<- phytools::reroot(ph, tmp, position=ph$edge.length[which(ph$edge[, 2]==tmp)] / 2)
	ph				<- ladderize(ph)
	pdf(gsub('newick','pdf',infile), w=15, h=500)
	plot(ph, cex=0.5, show.node.label=TRUE)
	dev.off()
	
	
	tmp				<- extract.clade(ph,getMRCA(ph,c("PG14-UG502292-S04296","PG14-UG502476-S04480")))
	tmp				<- merge(sqi, data.table(TAXA=tmp$tip.label), by='TAXA')
	tmp				<- merge(tmp, dc, by='PIDF')
	write.table(tmp, row.names=FALSE, file='/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/stragglers_no3.csv')
	
}

######################################################################################
project.Rakai.FastTree.170704<- function()
{
	require(ape)	
	require(data.table)
	require(big.phylo)
	
	#
	#	produce alignment
	if(0)
	{
		sq				<- read.dna("~/Dropbox (SPH Imperial College)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.fasta", format='fasta')
		#	handle repeats
		sqi				<- data.table(TAXA=rownames(sq), DUMMY=seq_len(nrow(sq)))
		tmp				<- sqi[, which(duplicated(TAXA))]
		cat('Repeats', paste(sqi[tmp,TAXA], collapse=', '))
		set(sqi, tmp, 'TAXA', sqi[tmp, paste(TAXA,'-R2',sep='')])
		setkey(sqi, DUMMY)
		rownames(sq)	<- sqi[,TAXA]
		#	get coverage
		tmp				<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
		sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]	
		sqi[, PNG:= sqi[, factor(grepl('PG',TAXA),levels=c(TRUE,FALSE),labels=c('Y','N'))]]
		sqi[, PIDF:= gsub('-R[0-9]+','',TAXA)]
		#	for each PIDF keep taxon with highest coverage
		sqi		<- merge(sqi, sqi[, list(DUMMY=DUMMY[which.max(COV)]), by='PIDF'], by=c('PIDF','DUMMY'))
		#	select
		sq				<- sq[ sqi$TAXA, ]
		rownames(sq)	<- gsub('-R[0-9]+','',rownames(sq))
		#
		#	add new consensuses
		infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/RakaiMiseqRerun_GlobalAln_10040.fasta'
		sn		<- read.dna(infile, format='fa')	
		infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/RakaiRerunKate_GlobalAln_10040.fasta'
		tmp		<- read.dna(infile, format='fa')
		sn		<- rbind(tmp,sn)
		infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/RakaiHiseqRerun_GlobalAln_10040.fasta'
		tmp		<- read.dna(infile, format='fa')
		sn		<- rbind(tmp,sn)
		rownames(sn)	<- gsub('_consensus_MinCov_10_50','',rownames(sn))
		rownames(sn)	<- gsub('_consensus','',rownames(sn))
		#	add to sqi list of taxa that were re-processed
		load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_assemblystatus.rda')	
		sqi		<- data.table(SID=rownames(sn))
		sqi		<- merge(dc, sqi, by='SID')
		sqi		<- subset(sqi, select=c(SID,PIDF))
		sqi[, TAXA:= paste0(PIDF,'_WTSI_',SID)]
		sn		<- sn[ sqi$SID, ]
		rownames(sn)	<- sqi$TAXA
		#
		sq		<- rbind(sq,sn)
		#
		#	get coverage
		sqi		<- data.table(TAXA=rownames(sq), DUMMY=seq_len(nrow(sq)))
		tmp		<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
		sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]
		sqi[, PNG:= sqi[, factor(grepl('PG',TAXA),levels=c(TRUE,FALSE),labels=c('Y','N'))]]	
		sqi[, PIDF:= gsub('_WTSI_.*','',TAXA)]
		#	for each PIDF identify taxon with highest coverage
		sqi		<- merge(sqi, sqi[, list(DUMMY=DUMMY[which.max(COV)], HIGHEST_COV='Y'), by='PIDF'], by=c('PIDF','DUMMY'), all.x=TRUE)
		set(sqi, sqi[, which(is.na(HIGHEST_COV))], 'HIGHEST_COV', 'N')
		#
		#	add GENBANK
		sqi[, GENBANK:= NA_character_]
		tmp				<- sqi[, which(PNG=='N')]
		set(sqi,tmp, 'GENBANK', sqi[tmp,gsub('.*\\.([^\\.]+)$','\\1',PIDF)])
		set(sqi,sqi[, which(PNG=='Y')],'GENBANK','PANGEA')	
		#	determine which references for UG alignment
		tmp				<- read.dna('/Users/Oliver/git/Phyloscanner.R.utilities/inst/HIV1_compendium_AD_B_CPX_v2.fasta', format='fa')
		tmp				<- rownames(tmp)
		tmp				<- gsub('.*_([^_]+)$','\\1',tmp[!grepl('CPX',tmp)])
		sqi[, IN_UG_ALI:='N']
		set(sqi, sqi[, which(GENBANK%in%tmp)], 'IN_UG_ALI', 'Y')
		set(sqi, sqi[, which(grepl('PG[0-9][0-9]-UG.*', TAXA))], 'IN_UG_ALI', 'Y')
		#
		#	write to file: Uganda, all SIDs, COV>700
		tmp			<- sq[ subset(sqi, COV>700 & IN_UG_ALI=='Y')[, TAXA], ]
		outfile		<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus','PANGEA_HIV_Imperial_v170704_UG_allSID_cov700.fasta')
		write.dna(tmp, file=outfile, format='fa')		
		#	write to file: Uganda, all SIDs, COV>300
		tmp			<- sq[ subset(sqi, COV>300 & IN_UG_ALI=='Y')[, TAXA], ]
		outfile		<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus','PANGEA_HIV_Imperial_v170704_UG_allSID_cov300.fasta')
		write.dna(tmp, file=outfile, format='fa')		
		#	write to file: Uganda, highest cov SIDs, COV>700
		tmp			<- sq[ subset(sqi, COV>700 & IN_UG_ALI=='Y' & HIGHEST_COV=='Y')[, TAXA], ]
		outfile		<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus','PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700.fasta')
		write.dna(tmp, file=outfile, format='fa')		
		#	write to file: All, highest cov SIDs, COV>0
		tmp			<- sq[ subset(sqi, COV>0 & HIGHEST_COV=='Y')[, TAXA], ]
		outfile		<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus','PANGEA_HIV_Imperial_v170704_ALL_bestcov_cov1.fasta')
		write.dna(tmp, file=outfile, format='fa')		
	}
	#
	#	run one FastTree to check
	if(0)
	{
		require(big.phylo)
		infile.fasta	<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus','PANGEA_HIV_Imperial_v170704_UG_allSID_cov700.fasta')
		infile.fasta	<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus','PANGEA_HIV_Imperial_v170704_UG_allSID_cov300.fasta')
		infile.fasta	<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus','PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700.fasta')	
		outfile.ft		<- gsub('\\.fasta','_ft.newick',infile.fasta)
		tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma')
		cat(tmp)
		
		
		#
		#	check FastTree
		infile			<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus','PANGEA_HIV_Imperial_v170704_UG_allSID_cov700_ft.newick')
		infile			<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus','PANGEA_HIV_Imperial_v170704_UG_allSID_cov300_ft.newick')
		infile			<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus','PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_ft.newick')	
		ph				<- read.tree(infile)
		tmp				<- as.numeric(ph$node.label)
		tmp[is.na(tmp)]	<- 0
		ph$node.label	<- tmp
		tmp				<- which(ph$tip.label=='B.FR.83.HXB2_LAI_IIIB_BRU.K03455')
		ph				<- phytools::reroot(ph, tmp, position=ph$edge.length[which(ph$edge[, 2]==tmp)] / 2)
		ph				<- ladderize(ph)
		pdf(gsub('newick','pdf',infile), w=15, h=500)
		plot(ph, cex=0.5, show.node.label=TRUE)
		dev.off()
		
		
		tmp				<- extract.clade(ph,getMRCA(ph,c("PG14-UG502292-S04296","PG14-UG502476-S04480")))
		tmp				<- merge(sqi, data.table(TAXA=tmp$tip.label), by='TAXA')
		tmp				<- merge(tmp, dc, by='PIDF')
		write.table(tmp, row.names=FALSE, file='/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/stragglers_no3.csv')		
	}
	#
	#	get 100 trees on HPC
	if(0)
	{			
		#indir			<- '/Users/Oliver/Dropbox (SPH Imperial College)/OR_Work/2017/2017_Seattle'
		indir			<- '/work/or105/Gates_2014/Rakai'
		#infile.fasta	<- file.path(indir,'PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700.fasta')
		infile.fasta	<- file.path(indir,'PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_couples.fasta')
		
		bs.n			<- 100	
		#bs.dir			<- file.path(indir,'PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_bootstrap_trees')
		bs.dir			<- file.path(indir,'PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_couples_bootstrap_trees')
		outfile.ft		<- gsub('\\.fasta',paste0('_ft_bs',bs.n,'.newick'),infile.fasta)
		tmp				<- cmd.fasttree.many.bootstraps(infile.fasta, bs.dir, bs.n, outfile.ft, pr.args='-nt -gtr -gamma', opt.bootstrap.by='nucleotide')
		#tmp				<- cmd.fasttree.many.bootstraps(infile.fasta, bs.dir, bs.n, bs.from=42, outfile.ft, pr.args='-nt -gtr -gamma', opt.bootstrap.by='nucleotide')
		
		
		#	run on HPC
		cmd				<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=998, hpc.q="pqeelab", hpc.mem="5900mb",  hpc.nproc=1, hpc.load='module load R/3.3.3')
		cmd				<- paste(cmd,tmp,sep='\n')
		cat(cmd)					
		outfile.cmd		<- paste("rk",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(indir, outfile.cmd, cmd)	
	}
	#
	#	get 100 ExaML trees on HPC
	if(1)
	{					
		indir			<- '/work/or105/Gates_2014/Rakai'
		infile			<- 'PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700'
		infile			<- 'PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_couples'
		#	ExaML bootstrap args
		bs.from		<- 0
		bs.to		<- 2
		bs.n		<- 100
		outdir		<- indir
		args.parser	<- paste("-m DNA")
		args.examl	<- "-f d -D -m GAMMA"	#	 -- this is the default that worked in 24 hours
		#args.examl	<- "-f d -D -m GTRCAT"			
		cmd			<- cmd.examl.bootstrap(indir, infile, bs.from=bs.from, bs.to=bs.to, bs.n=bs.n, opt.bootstrap.by="nucleotide", args.parser=args.parser, args.examl=args.examl, tmpdir.prefix="examl")					
		dummy		<- lapply(cmd, function(x)
				{				
					#x		<- cmd.hpcwrapper(x, hpc.walltime=21, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
					#x		<- cmd.hpcwrapper(x, hpc.walltime=79, hpc.q="pqeph", hpc.mem="1800mb", hpc.nproc=1)
					cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=998, hpc.q="pqeelab", hpc.mem="5900mb",  hpc.nproc=1, hpc.load='module load R/3.3.3 intel-suite mpi libgmp/4.2.1') 
					cmd		<- paste(cmd,x,sep='\n')
					cat(cmd)					
					outfile.cmd		<- paste("rke",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
					cmd.hpccaller(indir, outfile.cmd, cmd)
					Sys.sleep(1)
				})			
	}
	#
	#	get couples-only alignment 
	if(0)
	{	
		infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700.fasta'
		sq		<- read.dna(infile, format='fasta')
		dfi		<- data.table(F=rownames(sq))
		dfi[, PIDF:= gsub('_WTSI.*','',F)]
		dfi[, PNG:= grepl('^PG[0-9]+-', PIDF)]
		
		#	load couples to search for in phyloscanner output
		load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_info.rda")
		rps		<- unique(subset(rp, select=c(FEMALE_RID, MALE_RID, FEMALE_TAXA, MALE_TAXA, COUPID)))
		rps		<- data.table(PIDF=rps[, unique(na.omit(c(FEMALE_TAXA, MALE_TAXA)))])
		rps		<- merge(rps, dfi, by='PIDF')
		dfi		<- rbind(rps, subset(dfi, !PNG), fill=TRUE)		
		sq		<- sq[ dfi$F, ]
		
		write.dna(sq, file=gsub('\\.fasta','_couples.fasta',infile), format='fasta')		
	}
	
}

######################################################################################
project.Rakai.FastTree.170601<- function()
{
	require(ape)
	
	sq				<- read.dna("~/Dropbox (SPH Imperial College)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.fasta", format='fasta')
	#	handle repeats
	sqi				<- data.table(TAXA=rownames(sq), DUMMY=seq_len(nrow(sq)))
	tmp				<- sqi[, which(duplicated(TAXA))]
	cat('Repeats', paste(sqi[tmp,TAXA], collapse=', '))
	set(sqi, tmp, 'TAXA', sqi[tmp, paste(TAXA,'-R2',sep='')])
	setkey(sqi, DUMMY)
	rownames(sq)	<- sqi[,TAXA]
	#	get coverage
	tmp				<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
	sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]	
	sqi[, PNG:= sqi[, factor(grepl('PG',TAXA),levels=c(TRUE,FALSE),labels=c('Y','N'))]]
	sqi[, PIDF:= gsub('-R[0-9]+','',TAXA)]
	#	for each PIDF keep taxon with highest coverage
	sqi		<- merge(sqi, sqi[, list(DUMMY=DUMMY[which.max(COV)]), by='PIDF'], by=c('PIDF','DUMMY'))
	#	select
	sq				<- sq[ sqi$TAXA, ]
	rownames(sq)	<- gsub('-R[0-9]+','',rownames(sq))
	#	add new consensuses
	infile	<- '/Users/Oliver/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017/Rakai_Mar2017_GlobalAlignment_10040_BestSeqPerSample.fasta'
	sn		<- read.dna(infile, format='fa')
	sq		<- rbind(sq,sn)
	#	handle repeats
	sqi				<- data.table(TAXA=rownames(sq), DUMMY=seq_len(nrow(sq)))
	tmp				<- sqi[, which(duplicated(TAXA))]
	cat('Repeats', paste(sqi[tmp,TAXA], collapse=', '))
	set(sqi, tmp, 'TAXA', sqi[tmp, paste(TAXA,'-R2',sep='')])
	setkey(sqi, DUMMY)
	rownames(sq)	<- sqi[,TAXA]
	#	get coverage
	tmp				<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
	sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]	
	sqi[, PNG:= sqi[, factor(grepl('PG',TAXA),levels=c(TRUE,FALSE),labels=c('Y','N'))]]
	sqi[, PIDF:= gsub('-R[0-9]+','',TAXA)]
	#	for each PIDF keep taxon with highest coverage
	sqi				<- merge(sqi, sqi[, list(DUMMY=DUMMY[which.max(COV)]), by='PIDF'], by=c('PIDF','DUMMY'))
	#	select
	sq				<- sq[ sqi$TAXA, ]
	rownames(sq)	<- gsub('-R[0-9]+','',rownames(sq))
	#
	#	update stragglers
	infile			<- '~/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017/stragglers_GlobalAlignment_10040.fasta'
	sn				<- read.dna(infile, format='fa')
	rownames(sn)	<- gsub('_consensus_MinCov_10_50','',rownames(sn))	
	#	keep only longest straggler for each PANGEA ID
	stg				<- data.table(SID=rownames(sn), DUMMY=seq_len(nrow(sn)))
	tmp				<- sapply(seq_len(nrow(sn)), function(i) base.freq(sn[i,], all=TRUE, freq=TRUE))
	stg[, COV:=ncol(sn)-apply( tmp[c('-','?'),], 2, sum	)]		
	load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_assemblystatus.rda')
	tmp				<- merge(sqi, subset(dc, !is.na(SID), select=c(PIDF, SID)), by='PIDF',all.x=TRUE)
	stg				<- merge(stg, subset(tmp, select=c(SID, PIDF)), by='SID',all.x=1)
	stg				<- merge(stg, stg[, list(DUMMY=DUMMY[which.max(COV)]), by='PIDF'], by=c('PIDF','DUMMY'))
	#	add to re-processed stragglers those taxa that are now junk
	tmp				<- data.table(PIDF=c('PG14-UG501871-S03875','PG14-UG500718-S02722','PG14-UG503751-S05755','PG14-UG501956-S03960','PG14-UG503966-S05970','PG14-UG503944-S05948','PG14-UG503322-S05326','PG14-UG501838-S03842','PG14-UG503896-S05900','PG14-UG500743-S02747','PG14-UG503259-S05263','PG14-UG502371-S04375','PG14-UG503509-S05513','PG14-UG503648-S05652','PG14-UG503750-S05754','PG14-UG503518-S05522','PG14-UG501117-S03121','PG14-UG501506-S03510','PG14-UG503995-S05999','PG14-UG502131-S04135','PG14-UG501758-S03762','PG14-UG503968-S05972','PG14-UG503565-S05569','PG14-UG502024-S04028','PG14-UG503929-S05933','PG14-UG503673-S05677','PG14-UG502418-S04422','PG14-UG503643-S05647','PG14-UG501140-S03144','PG14-UG503190-S05194','PG14-UG500711-S02715','PG14-UG503891-S05895','PG14-UG503697-S05701','PG14-UG501743-S03747','PG14-UG503596-S05600','PG14-UG502370-S04374','PG14-UG501271-S03275','PG14-UG503555-S05559','PG14-UG502472-S04476','PG14-UG503500-S05504'))
	stg				<- merge(tmp, stg, by='PIDF', all.x=1)
	write.csv(stg, file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/stragglers_no1.csv')
	#	reduce central alignment to taxa not in stg	
	sq				<- sq[setdiff(rownames(sq),stg[, PIDF]),]
	tmp				<- subset(stg, !is.na(SID))
	sn				<- sn[tmp[,SID],]
	rownames(sn)	<- tmp[,PIDF]
	sq				<- rbind(sq,sn)	
	sqi				<- merge(sqi,data.table(PIDF=rownames(sq)), by='PIDF')	
	sqi[, COV:=NULL]
	tmp				<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
	sqi[, COV:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum	)]	
	
	#	
	#	
	#	save all
	tmp				<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017','PANGEA_HIV_n5602_Imperial_v170601_GlobalAlignment.fasta')
	write.dna(sq, file=tmp, format='fa')	
	#
	#	save UG
	sqi[, GENBANK:= NA_character_]
	tmp				<- sqi[, which(PNG=='N')]
	set(sqi,tmp, 'GENBANK', sqi[tmp,gsub('.*\\.([^\\.]+)$','\\1',PIDF)])
	set(sqi,sqi[, which(PNG=='Y')],'GENBANK','PANGEA')
	tmp				<- read.dna('/Users/Oliver/git/Phyloscanner.R.utilities/inst/HIV1_compendium_AD_B_CPX_v2.fasta', format='fa')
	tmp				<- rownames(tmp)
	tmp				<- gsub('.*_([^_]+)$','\\1',tmp[!grepl('CPX',tmp)])
	sqi				<- merge(sqi,data.table(GENBANK=c(tmp,'PANGEA')),by='GENBANK')		
	sq				<- sq[sqi$PIDF,]
	sqi				<- subset(sqi, PNG=='N' | (PNG=='Y' & grepl('^PG[0-9][0-9]-UG.*',TAXA)))
	sq				<- sq[sqi$PIDF,]
	tmp				<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017','PANGEA_HIV_n4990_Imperial_v170601_UgandaAlignment.fasta')
	write.dna(sq, file=tmp, format='fa')	
	
	#
	#	first check FastTree
	if(0)
	{
		infile.fasta	<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017','PANGEA_HIV_n5014_Imperial_v170601_UgandaAlignment.fasta')
		outfile.ft		<- gsub('\\.fasta','_ft.newick',infile.fasta)
		tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma')
		outfile.ft		<- gsub('\\.fasta','_ft-pseudo.newick',infile.fasta)
		tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma -pseudo')
		outfile.ft		<- gsub('\\.fasta','_ft-pseudo-spr4-mlacc2-slownni.newick',infile.fasta)
		tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma -pseudo -spr 4 -mlacc 2 -slownni')
		#	not much difference between the options, with standard tree having reasonable LKL
		#	all trees show very long branches sticking out, which could be due to missing nucleotides		
	}
	
	tmp2			<- sq[subset(sqi, COV>700)[,PIDF],]
	tmp				<- file.path('/Users/Oliver/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017','PANGEA_HIV_n4199_Imperial_v170601_UgandaAlignment_minCov700.fasta')	
	write.dna(tmp2, file=tmp, format='fa')	
	infile.fasta	<- '/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_n4199_Imperial_v170601_UgandaAlignment_minCov700.fasta'
	outfile.ft		<- gsub('\\.fasta','_ft.newick',infile.fasta)
	tmp				<- cmd.fasttree(infile.fasta, outfile=outfile.ft, pr.args='-nt -gtr -gamma')
	
	infile			<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_n4199_Imperial_v170601_UgandaAlignment_minCov700_ft.newick'
	ph				<- read.tree(infile)
	tmp				<- as.numeric(ph$node.label)
	tmp[is.na(tmp)]	<- 0
	ph$node.label	<- tmp
	ph				<- ladderize(ph)
	pdf(gsub('newick','pdf',infile), w=15, h=500)
	plot(ph, cex=0.5, show.node.label=TRUE)
	dev.off()
}

project.RakaiAll.forTanya<- function()
{
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75.rda')
	tmp		<- copy(pty.runs)	
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75_part2.rda')
	pty.runs<- rbind(tmp, pty.runs)
	
	dc
	
	#	load demographic info on all individuals
	tmp		<- RakaiCirc.epi.get.info.170208()	
	rd		<- tmp$rd
	rn		<- tmp$rn	
	rd		<- rbind(rd, rn, use.names=TRUE, fill=TRUE)	
	#	remove controls
	rd			<- subset(rd, !(HIVPREV==0 & PANGEA==1))
	#	fixup rd: 
	#	missing first pos date
	tmp			<- rd[, which(is.na(FIRSTPOSDATE) & !is.na(RECENTVLDATE) & TIMESINCEVL==0)]
	set(rd, tmp, 'FIRSTPOSDATE', rd[tmp, RECENTVLDATE])	
	stopifnot(!nrow(subset(rd, is.na(FIRSTPOSDATE))))	
	#	fixup rd: 
	#	there are duplicate RID entries with missing FIRSTPOSDATE, and ambiguous ARVSTARTDATE; or inconsistent across VISIT entries
	#	missing FIRSTPOSDATE -> delete
	#	ambiguous ARVSTARTDATE -> keep earliest	
	tmp		<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
	tmp[, DUMMY:=seq_len(nrow(tmp))]
	tmp		<- merge(tmp, tmp[, {
						ans	<- is.na(FIRSTPOSDATE)	
						if(any(!is.na(ARVSTARTDATE)))
							ans[!is.na(ARVSTARTDATE) & ARVSTARTDATE!=min(ARVSTARTDATE, na.rm=TRUE)]	<- TRUE	
						list(DUMMY=DUMMY, DELETE=ans)		
					}, by=c('RID')], by=c('RID','DUMMY'))
	tmp		<- subset(tmp, !DELETE)
	set(tmp, NULL,c('DUMMY','DELETE'), NULL)
	set(rd, NULL, c('BIRTHDATE','LASTNEGDATE','FIRSTPOSVIS','FIRSTPOSDATE','ARVSTARTDATE','EST_DATEDIED'), NULL)
	rd		<- merge(rd, tmp, by='RID')	
	rd		<- unique(subset(rd, select=c("RID","VISIT","DATE","SEX","LASTNEGVIS","RECENTCD4","RECENTCD4DATE","RECENTVL","RECENTVLDATE","SELFREPORTART","EVERSELFREPORTART","FIRSTSELFREPORTART","RELIGION","COHORT","BIRTHDATE","LASTNEGDATE","FIRSTPOSVIS","FIRSTPOSDATE","ARVSTARTDATE","EST_DATEDIED")), by=c('RID','VISIT'))
	 
	
	save(rd, dc, pty.runs, file='/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_forTanya.rda')
	
	
	
}

project.RakaiAll.setup.phyloscanner.170704.cleanup.old.RakaiBamFiles<- function()
{
	# get all taxa processed
	# these appear in the global alignments, some may just have ??? throughout
	
	# load global alignments
	infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/RakaiMiseqRerun_GlobalAln_10040.fasta'
	sn		<- read.dna(infile, format='fa')	
	infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/RakaiRerunKate_GlobalAln_10040.fasta'
	tmp		<- read.dna(infile, format='fa')
	sn		<- rbind(tmp,sn)
	infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/RakaiHiseqRerun_GlobalAln_10040.fasta'
	tmp		<- read.dna(infile, format='fa')
	sn		<- rbind(tmp,sn)
	rownames(sn)	<- gsub('_consensus_MinCov_10_50','',rownames(sn))
	rownames(sn)	<- gsub('_consensus','',rownames(sn))
	#	get coverage
	sqi				<- data.table(SID=rownames(sn))
	tmp				<- sapply(seq_len(nrow(sn)), function(i) base.freq(sn[i,], all=TRUE, freq=TRUE))
	sqi[, COV:=ncol(sn)-apply( tmp[c('-','?'),], 2, sum	)]	
	
	#	make table: all old files, which newly processed, which newly processed no coverage w/o LTR
	di		<- as.data.table(read.table(file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/listOfOldBams.txt', col.names='BAM', stringsAsFactors=FALSE))
	di[, SID:= gsub('.bam','',BAM)]	
	tmp		<- as.data.table(read.table(file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/listOfOldRefs.txt', col.names='REF', stringsAsFactors=FALSE))
	tmp[, SID:= gsub('_ref.fasta','',REF)]
	di		<- merge(di,tmp,by='SID')
	tmp		<- as.data.table(read.table(file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/listOfBams.txt', col.names='BAM_NEW', stringsAsFactors=FALSE))
	tmp[, SID:= gsub('.bam|_remap.bam','',BAM_NEW)]	
	stopifnot( !nrow(subset(tmp[,list(N_BAM=length(BAM_NEW)),by='SID'], N_BAM!=2)) )
	tmp		<- subset(tmp, !grepl('remap',BAM_NEW))	
	tmp2	<- as.data.table(read.table(file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/listOfRefs.txt', col.names='REF_NEW', stringsAsFactors=FALSE))
	tmp2[, SID:= gsub('_ref.fasta|_remap_ref.fasta','',REF_NEW)]	
	stopifnot( !nrow(subset(tmp2[,list(N_REF=length(REF_NEW)),by='SID'], N_REF!=2)) )
	tmp2	<- subset(tmp2, !grepl('remap',REF_NEW))
	tmp		<- merge(tmp,tmp2,by='SID')
	tmp		<- merge(tmp, sqi, by='SID')
	di		<- merge(di,tmp,by='SID',all=TRUE)
	di[, WTS:= grepl('^[0-9]+_[0-9]_[0-9]+$',SID)]
	
	
	#	only valid Sanger files
	di		<- subset(di, WTS)
	#	old files with an update
	subset(di, !is.na(BAM) & !is.na(BAM_NEW) & COV>0)
	#	old files to be deleted
	cat(subset(di, !is.na(BAM) & !is.na(BAM_NEW) & COV==0, SID)[, paste(SID, collapse=' ')], file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/listOfOldSIDsToBeDeleted.txt')
	#write.table(subset(di, !is.na(BAM) & !is.na(BAM_NEW) & COV==0, SID), row.names=FALSE, col.names=FALSE, file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/listOfOldSIDsToBeDeleted.txt')
}

project.RakaiAll.setup.phyloscanner.170301<- function()
{
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301'
	#	start with latest Sanger IDs
	dc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/WTSI_PANGEA_InfoFind_2017-02-14.csv', header=TRUE, stringsAsFactors=FALSE))
	setnames(dc, c('Lane','Public'), c('SID','PIDF'))
	dc		<- subset(dc, select=c(SID,PIDF))
	set(dc, NULL, 'SID', dc[, gsub('#','_',SID)])
	set(dc, NULL, 'PID', dc[, gsub('-S[0-9]+$','',PIDF)])
	#
	#	list of files that Chris has
	#
	infiles	<- data.table(F=list.files(indir, pattern='^CW_PANGEA',full.names=TRUE))
	tmp		<- do.call('rbind',lapply(infiles[, F], function(file){
							tmp<- as.data.table(read.csv(file, header=FALSE, col.names='PID', stringsAsFactors=FALSE))
							tmp[, F:=file]
							tmp
					}))
	set(tmp, NULL, 'PROC_STATUS', tmp[, gsub('CW_PANGEA_Rakai_','',gsub('\\.txt','',basename(F)))])	
	set(tmp, NULL, 'F', NULL)
	#	all files for which we have some PID
	dc		<- merge(dc, tmp, by='PID',all=1)	
	stopifnot(!nrow(subset(dc, is.na(SID) & PROC_STATUS!='ThoseWithoutFastqs')))
	dc		<- subset(dc, is.na(SID) | SID!='15351_1_1')		#	remove 15351_1_1  with no PANGEA info in WTSI file
	dc		<- subset(dc, is.na(SID) | SID!='15430_1_75')	#	remove 15430_1_75  with no PANGEA info in WTSI file
	stopifnot(!nrow(subset(dc, is.na(PIDF) & !is.na(SID))))
	dc[, PART:= as.numeric(gsub('^[0-9]+_([0-9])_[0-9]+','\\1',SID))]
	dc[, DUMMY:= gsub('^([0-9]+)_[0-9]_([0-9]+)','\\1_x_\\2',SID)]	
	dc		<- merge(dc, dc[, list(N_PART=length(PART), S_PART=sum(PART)), by='DUMMY'], by='DUMMY')	
	#	check if we always have _1_ and _2_
	stopifnot(!nrow(subset(dc, N_PART==2 & S_PART!=3)))
	#	assume Tanya merges to _3_
	tmp		<- dc[, which(N_PART==2)]
	set(dc, tmp, 'SID', dc[tmp, gsub('^([0-9]+)_[0-9]_([0-9]+)','\\1_3_\\2',SID)])
	set(dc, NULL, c('PART','N_PART','S_PART','DUMMY'),NULL)
	dc		<- unique(dc)
	tmp		<- subset(dc, !is.na(SID))[, list(N_SID=length(SID)), by='PIDF']
	dc		<- merge(dc, tmp, by='PIDF',all.x=1)
	#	define controls
	set(dc, dc[, which(grepl('neg',PID))], 'PROC_STATUS', 'NegControl')
	#	define not in Chris census
	set(dc, dc[, which(is.na(PROC_STATUS))], 'PROC_STATUS', 'NotTrackedByChris')	
	#	add extra category to PROC_STATUS: not processed by Kate
	tmp		<- as.data.table(read.table("~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/KG_PANGEA_Processed_597.txt", header=TRUE,stringsAsFactors=FALSE))
	setnames(tmp, c('SampleID','LaneID'), c('PID','SID'))
	tmp[, KATE_PROC:='Y']
	dc		<- merge(dc, tmp, by=c('PID','SID'), all.x=1)	
	set(dc, dc[, which(PROC_STATUS=='ThoseWithFastqs_WithKateShiverOutput' & is.na(KATE_PROC))], 'PROC_STATUS','ThoseWithFastqs_KateNotProcessed') 
	set(dc, NULL, 'KATE_PROC', NULL)	
	#	add RIDs
	load('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda')
	tmp		<- subset(as.data.table(rccsData), select=c(RCCS_studyid,Pangea.id, batch, date, SEX))
	setnames(tmp, c('RCCS_studyid','Pangea.id','batch','date'), c('RID','PID','RCCS_SHIP_BATCH','SAMPLE_DATE'))
	tmp		<- subset(tmp, !is.na(PID))
	tmp		<- unique(tmp, by=c('RID','PID'))
	tmp2	<- subset(as.data.table(neuroData), select=c(studyid, Pangea.id, sampleDate, gender))
	setnames(tmp2, c('studyid','Pangea.id','sampleDate','gender'), c('RID','PID','SAMPLE_DATE','SEX'))
	tmp2[, RCCS_SHIP_BATCH:='neuro']
	tmp		<- rbind(tmp, tmp2, use.names=TRUE)
	set(tmp, NULL, 'RID', tmp[, as.character(RID)])
	set(tmp, NULL, 'PID', tmp[, as.character(PID)])
	set(tmp, NULL, 'SEX', tmp[, as.character(SEX)])
	dc		<- merge(dc, tmp, by='PID',all.x=1)
	#	flag test plate
	set(dc, dc[, which(grepl('PG14-UG9000[0-9][0-9]',PID))], 'RCCS_SHIP_BATCH', 'test')
	#	see if on HPC	
	tmp		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/HPC_census_bams.txt', header=FALSE, col.names='SID', stringsAsFactors=FALSE))
	tmp[, HPC_BAM:='Y']
	set(tmp, NULL, 'SID', tmp[, gsub('\\.bam','',SID)])
	tmp		<- subset(tmp, grepl('^[0-9]+_[0-9]_[0-9]+',SID))	#	reduce to bams from SANGER
	stopifnot( !length(setdiff( tmp[, SID], dc[, SID] )) )
	dc		<- merge(dc, tmp, by='SID',all.x=TRUE)
	set(dc, dc[, which(is.na(HPC_BAM))], 'HPC_BAM', 'N')	
	tmp		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/HPC_census_refs.txt', header=FALSE, col.names='SID', stringsAsFactors=FALSE))
	tmp[, HPC_REF:='Y']
	set(tmp, NULL, 'SID', tmp[, gsub('_ref.fasta','',SID)])
	tmp		<- subset(tmp, grepl('^[0-9]+_[0-9]_[0-9]+',SID))	#	reduce to refs from SANGER	
	stopifnot( !length(setdiff( tmp[, SID], dc[, SID] )) )	
	dc		<- merge(dc, tmp, by='SID',all.x=TRUE)
	set(dc, dc[, which(is.na(HPC_REF))], 'HPC_REF', 'N')
	tmp		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/HPC_census_fastq.txt', header=FALSE, col.names='SID', stringsAsFactors=FALSE))
	tmp[, HPC_FASTQ:='Y']
	set(tmp, NULL, 'SID', tmp[, gsub('_[0-9].fastq.gz','',SID)])	
	tmp		<- subset(tmp, grepl('^[0-9]+_[0-9]_[0-9]+',SID))	#	reduce to fastqz's from SANGER	
	dc		<- merge(dc, unique(tmp), by='SID',all.x=TRUE)
	set(dc, dc[, which(is.na(HPC_FASTQ))], 'HPC_FASTQ', 'N')	
	#	add latest PANGEA stats
	tmp		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/PANGEA_data/2016-07-07_PANGEA_stats_by_sample.csv', stringsAsFactors=FALSE))
	setnames(tmp, 	c('Status','Submitted','DayDiff','ProjectID','Cohort','Submitted.1','Sequenced','Assembled','HIVcontig'), 
					c('WTSI_STATUS','WTSI_SUBMITTED_DATE','DayDiff','PID','Cohort','WTSI_SUBMITTED','WTSI_SEQUENCED','WTSI_ASSEMBLED','WTSI_HIVCONTIG'))
	tmp		<- subset(tmp, PID!='')
	set(tmp, NULL, 'WTSI_SUBMITTED', tmp[, as.character(factor(WTSI_SUBMITTED=='',levels=c(TRUE,FALSE),labels=c('N','Y')))])
	set(tmp, NULL, 'WTSI_SEQUENCED', tmp[, as.character(factor(WTSI_SEQUENCED=='',levels=c(TRUE,FALSE),labels=c('N','Y')))])
	set(tmp, NULL, 'WTSI_ASSEMBLED', tmp[, as.character(factor(WTSI_ASSEMBLED=='',levels=c(TRUE,FALSE),labels=c('N','Y')))])
	set(tmp, NULL, 'WTSI_HIVCONTIG', tmp[, as.character(factor(WTSI_HIVCONTIG=='',levels=c(TRUE,FALSE),labels=c('N','Y')))])
	set(tmp, NULL, c('DayDiff','Cohort'), NULL)
	dc		<- merge(dc, tmp, by='PID', all.x=1)
	#	check what Dan assembled HISEQ
	tmp		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/PANGEA_data/PANGEA_UCL_Feb2017_collated_stats_all_genomes_UCL_release_Feb2017_allhitodate.csv', stringsAsFactors=FALSE))
	setnames(tmp, c('PANGEA_ID','PGID_full','WTSI_ID','Length', 'Cohort'), c('PID','PIDF','SID','UCL_LEN', 'COHORT'))
	tmp		<- subset(tmp, select=c('PID','PIDF','SID','UCL_LEN', 'COHORT'))
	#	convert SID to _3_ to match our convention
	set(tmp, NULL, 'SID', tmp[, gsub('^([0-9]+)_[0-9]_([0-9]+)','\\1_3_\\2',SID)])
	#	check what Dan assembled MISEQ
	tmp2	<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/PANGEA_data/PANGEA_UCL_Feb2017_collated_stats_all_genomes_UCL_release_Feb2017_allmitodate.csv', stringsAsFactors=FALSE))
	setnames(tmp2, c('PANGEA_ID','PGID_full','WTSI_lane_ID','Length', 'Cohort'), c('PID','PIDF','SID','UCL_LEN', 'COHORT'))
	tmp2	<- subset(tmp2, select=c('PID','PIDF','SID','UCL_LEN', 'COHORT'))	
	tmp		<- rbind(tmp, tmp2)
	stopifnot( !length(setdiff(tmp[, SID], dc[, SID])) )	
	#	n=6 plus test plate
	if(0)
	{
		tmp		<- subset(tmp, grepl('Rakai',COHORT))
		tmp2	<- setdiff( tmp[, unique(sort(PID))], dc[, unique(sort(PID))] )
		write.csv(data.table(ID=tmp2), file=file.path(indir,'IDs_that_Dan_flags_from_Rakai_but_not_in_Chris_census.csv'), row.names=FALSE)			
	}
	set(tmp, NULL, 'COHORT', NULL)		
	dc		<- merge(dc, tmp, by=c('PID','PIDF','SID'),all.x=1)	
	#	resolve 60 from Chris that he confirmed he has not processed -- at least one SID from said individual has been processed
	set(dc, dc[, which(PROC_STATUS=='ThoseWithFastqs_WithChrisShiverOutput' & HPC_BAM=='N')], 'PROC_STATUS','ThoseWithFastqs_ChrisNotProcessed')
	#
	#	subset to RAKAI
	#
	dc		<- subset(dc, !is.na(RID))
	#	define RIDs from whom all SIDs are complete and on HPC
	tmp		<- dc[, list(	HPC_ALL_SID_FOR_RID= all(HPC_BAM=='Y') & all(HPC_REF=='Y')), by='RID']	
	set(tmp, NULL, 'HPC_ALL_SID_FOR_RID', tmp[,as.character(factor(HPC_ALL_SID_FOR_RID, levels=c(TRUE,FALSE), labels=c('Y','N')))])
	dc		<- merge(dc, tmp, by='RID',all.x=1)
	#	define Sampling Time
	set(dc, NULL, 'WTSI_SUBMITTED_DATE', dc[, as.Date(WTSI_SUBMITTED_DATE)])	
	stopifnot(!nrow(subset(dc, HPC_BAM!=HPC_REF)))
	#	add column consensus available + length consensus	
	tmp		<- '/Users/Oliver/Dropbox (SPH Imperial College)/PANGEA_data/Rakai_Shiver_Mar2017/PANGEA_HIV_n5014_Imperial_v170601_UgandaAlignment.fasta'
	sq		<- read.dna(tmp, format='fasta')	
	sqi		<- data.table(PIDF=rownames(sq))	
	tmp		<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
	sqi[, IMPERIAL_LEN:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum)]
	dc		<- merge(dc, sqi, by='PIDF', all.x=1)	
	save(dc, file=file.path(indir,'Rakai_phyloscanner_170301_assemblystatus.rda'))
	
	if(0)	#consensus stats
	{
		subset(dc, !is.na(SID) & is.na(IMPERIAL_LEN))	#213
		subset(dc, !is.na(SID) & !is.na(IMPERIAL_LEN) & IMPERIAL_LEN<700) #766
		subset(dc, !is.na(SID) & IMPERIAL_LEN*1.2<UCL_LEN)	#1493
		subset(dc, !is.na(SID) & IMPERIAL_LEN>UCL_LEN*1.2)	#19
	}
	if(0)
	{
		#	add PANGEA2 IDs
		#tmp		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/OR_Work/2016/2016_Rakai_Couples/Old_New_PANGEA_ID_Linkage_Table.csv', stringsAsFactors=FALSE))
		#setnames(tmp, c('Pangea_ID','New_ID'), c('PID','PIDF2'))
		#set(tmp, NULL, 'PID2', tmp[,gsub('-[0-9]+$','',PIDF2)])	
		#dc		<- merge(dc, unique(subset(tmp, select=c(PID,PID2))), by='PID', all.x=1)		
	}
	if(0)
	{
		#
		#	Kate to check if we have RIDs= 240
		#	resolved --> these are neuro
		write.csv(	subset(dc, is.na(RID)), 
				file=file.path(indir,'Kate_RIDs_needed.csv'), row.names=FALSE)
		#	Kate to check samples that Dan could assemble
		#	n=134
		#	NOT YET RESOLVED
		write.csv(	subset(dc, !is.na(UCL_LEN) & PROC_STATUS=='ThoseWithFastqs_WithKateShiverOutput' & HPC_BAM=='N'), 
				file=file.path(indir,'Kate_check_unassembled_among_597_that_Dan_could_assemble.csv'), row.names=FALSE)
		#	WTSI gone missing=16 -- check with Swee Hoe
		#	Anne suggests they 'failed' sequencing. I am not sure what this means
		write.csv(subset(dc, !is.na(RID) & is.na(SID) & WTSI_STATUS=='Assume sequencing failed'), file=file.path(indir,'WTSI_gone_missing.csv'), row.names=FALSE)		
		#	New list for Tanya of unprocessed SIDs that are not in Kate s batch
		#	n=1617
		write.csv(subset(dc, PROC_STATUS%in%c('ThoseWithFastqs_KateNotProcessed','ThoseWithFastqs_ChrisNotProcessed','ThoseWithFastqs_WithoutShiverOutput')), file=file.path(indir,'SHIVER_to_run_on_existing_FASTQ.csv'), row.names=FALSE)
	}
	#
	#	set up first batches
	#
	if(0)
	{
		batch.n	<- 15
		dc.it1	<- subset(dc, HPC_ALL_SID_FOR_RID=='Y')
		tmp		<- unique(subset(dc.it1, select=RID))	
		set.seed(42)
		tmp[, RID_B:= sample(nrow(tmp),nrow(tmp))]
		dc.it1	<- merge(dc.it1, tmp, by='RID')
		setkey(dc.it1, RID_B)
		dc.it1[, BATCH:= ceiling( RID_B/batch.n )]
		#
		#	set up first pty.runs
		#		
		pty.runs	<- as.data.table(t(combn(dc.it1[, unique(BATCH)],2)))
		setnames(pty.runs, c('V1','V2'), c('BATCH','BATCH2'))	
		pty.runs[, PTY_RUN:= seq_len(nrow(pty.runs))]
		pty.runs	<- melt(pty.runs, id.vars='PTY_RUN', variable.name='DUMMY', value.name='BATCH')
		set(pty.runs, NULL, 'DUMMY', NULL)
		setkey(pty.runs, PTY_RUN)
		tmp			<- subset(dc.it1, select=c(BATCH, RID, SID))
		tmp			<- tmp[, list(SID=SID, RENAME_SID=paste0(RID,'_fq',seq_along(SID))), by=c('BATCH','RID')]
		pty.runs	<- merge(pty.runs, tmp, by='BATCH',allow.cartesian=TRUE)
		
		#tmp	<- pty.runs[, list(N_SID=length(SID)), by=c('PTY_RUN','RID')]	
		save(dc, dc.it1, pty.runs, file=file.path(indir,'Rakai_phyloscanner_170301.rda'))		
	}
	
	#
	#	set up first batches
	#
	if(0)
	{
		#	input: dc is a large data.table with column RID (ID of individuals) and other columns 
		batch.n	<- 75			#number of individuals per phyloscanner run
		dc.it1	<- subset(dc, HPC_ALL_SID_FOR_RID=='Y')		#select indivivuals for whom we have already bam files 
		tmp		<- unique(subset(dc.it1, select=RID))	#tmp has one row with column RID
		set.seed(42)
		tmp[, RID_B:= sample(nrow(tmp),nrow(tmp))]	# randomly associate a batch ID to each RID 	
		dc.it1	<- merge(dc.it1, tmp, by='RID')
		setkey(dc.it1, RID_B)
		dc.it1[, BATCH:= ceiling( RID_B/batch.n )]	# convert RID_B into a batch ID, grouping batch.n individuals into one BATCH	
		#
		#	set up first pty.runs
		#		
		pty.runs	<- as.data.table(t(combn(dc.it1[, unique(BATCH)],2)))	# take all combinations of BATCHes
		setnames(pty.runs, c('V1','V2'), c('BATCH','BATCH2'))	
		pty.runs[, PTY_RUN:= seq_len(nrow(pty.runs))]	# assign phyloscanner IDs to all combinations of BATCHes
		pty.runs	<- melt(pty.runs, id.vars='PTY_RUN', variable.name='DUMMY', value.name='BATCH')	# melt into long format
		set(pty.runs, NULL, 'DUMMY', NULL)
		setkey(pty.runs, PTY_RUN)
		tmp			<- subset(dc.it1, select=c(BATCH, RID, SID))
		tmp			<- tmp[, list(SID=SID, RENAME_SID=paste0(RID,'_fq',seq_along(SID))), by=c('BATCH','RID')]	# at my end, individual IDs and bam IDs are different, and I want to rename the short reads in each bam ID too. Here I define the RENAME_SID field
		pty.runs	<- merge(pty.runs, tmp, by='BATCH',allow.cartesian=TRUE)
		
		#tmp	<- pty.runs[, list(N_SID=length(SID)), by=c('PTY_RUN','RID')]	
		save(dc, dc.it1, pty.runs, file=file.path(indir,'Rakai_phyloscanner_170301_b75.rda'))
		
		subset(pty.runs, PTY_RUN==1)[, paste0(SID,'.bam', collapse=' ')]
		subset(pty.runs, PTY_RUN==1)[, paste0(SID,'_ref.fasta', collapse=' ')]
	}
	#
	#	set up second batch
	#
	if(1)
	{
		dc.new	<- copy(dc)
		#	update batches
		load( file.path(indir,'Rakai_phyloscanner_170301_b75.rda') )		
		stopifnot(!length(setdiff(dc[, RID], dc.new[, RID])))
		batch.n	<- 75
		batch.n2<- 55	# use smaller number because we typically have multiple seqs from these individuals
		dc.it2	<- subset(dc.new, HPC_BAM=='Y' & HPC_REF=='Y')
		tmp		<- unique(subset(dc.it2, select=RID))	
		tmp		<- merge(tmp, unique(subset(dc.it1, select=c(RID, RID_B, BATCH))), by='RID', all.x=1)
		tmp2	<- tmp[, which(is.na(RID_B))]
		set.seed(42)
		tmp3	<- tmp[, ceiling(max(RID_B, na.rm=TRUE) / batch.n)*batch.n]
		set(tmp, tmp2, 'RID_B', as.integer(sample(length(tmp2),length(tmp2)) + tmp3) )
		dc.it2	<- merge(dc.it2, tmp, by='RID')
		setkey(dc.it2, RID_B)
		tmp2	<- dc.it2[, which(is.na(BATCH))]
		tmp3	<- min( dc.it2[tmp2, ceiling(RID_B/batch.n2)] ) - dc.it1[, max(BATCH)] - 1L		
		set(dc.it2, tmp2, 'BATCH', dc.it2[tmp2, ceiling(RID_B/batch.n2)] - tmp3)
		set(dc.it2, dc.it2[, which(BATCH==max(BATCH))], 'BATCH', dc.it2[, max(BATCH)-1L])		
		#
		#	set up remaining pty.runs
		#		
		tmp		<- as.data.table(t(combn(dc.it2[, unique(BATCH)],2)))
		setnames(tmp, c('V1','V2'), c('BATCH','BATCH2'))
		pty.runs<- unique(subset(pty.runs, select=c('BATCH','PTY_RUN')))
		pty.runs<- pty.runs[, list(BATCH=BATCH[1], BATCH2=BATCH[2]), by='PTY_RUN']		
		pty.runs<- merge(tmp, pty.runs, by=c('BATCH','BATCH2'), all.x=1)
		tmp		<- pty.runs[, which(is.na(PTY_RUN))]
		set(pty.runs, tmp, 'PTY_RUN', seq_along(tmp)+pty.runs[, max(PTY_RUN, na.rm=1)] )
		#	select only those that are not yet processed
		tmp			<- dc.it1[, max(BATCH)]
		pty.runs	<- subset(pty.runs, BATCH>tmp | BATCH2>tmp)		
		
		pty.runs	<- melt(pty.runs, id.vars='PTY_RUN', variable.name='DUMMY', value.name='BATCH')
		set(pty.runs, NULL, 'DUMMY', NULL)
		setkey(pty.runs, PTY_RUN)
		tmp			<- subset(dc.it2, select=c(BATCH, RID, SID))
		tmp			<- tmp[, list(SID=SID, RENAME_SID=paste0(RID,'_fq',seq_along(SID))), by=c('BATCH','RID')]
		pty.runs	<- merge(pty.runs, tmp, by='BATCH',allow.cartesian=TRUE)
		#	save
		dc			<- copy(dc.new)
		save(dc, dc.it2, pty.runs, file=file.path(indir,'Rakai_phyloscanner_170301_b75_part2.rda'))		
	}

}

project.RakaiAll.setup.phyloscanner.171219<- function()
{
	#	load all UVRI files on HPC
	
	du		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/UVRI_bams.txt', header=FALSE, col.names='BAM', stringsAsFactors=FALSE))
	set(du, NULL, 'SID', du[, gsub('\\.bam','',basename(BAM))])
	
	tmp[, HPC_BAM:='Y']
	
	tmp		<- subset(tmp, grepl('^[0-9]+_[0-9]_[0-9]+',SID))	#	reduce to bams from SANGER
	stopifnot( !length(setdiff( tmp[, SID], dc[, SID] )) )
	dc		<- merge(dc, tmp, by='SID',all.x=TRUE)
	set(dc, dc[, which(is.na(HPC_BAM))], 'HPC_BAM', 'N')	
	tmp		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/listOfRefs170704.txt', header=FALSE, col.names='SID', stringsAsFactors=FALSE))
	tmp[, HPC_REF:='Y']
	set(tmp, NULL, 'SID', tmp[, gsub('_ref.fasta','',SID)])
	tmp		<- subset(tmp, grepl('^[0-9]+_[0-9]_[0-9]+',SID))	#	reduce to refs from SANGER	
	stopifnot( !length(setdiff( tmp[, SID], dc[, SID] )) )	
	dc		<- merge(dc, tmp, by='SID',all.x=TRUE)
	set(dc, dc[, which(is.na(HPC_REF))], 'HPC_REF', 'N')
	
}

project.RakaiAll.setup.phyloscanner.170704<- function()
{
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301'
	#	start with latest Sanger IDs
	dc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/WTSI_PANGEA_InfoFind_2017-02-14.csv', header=TRUE, stringsAsFactors=FALSE))
	setnames(dc, c('Lane','Public'), c('SID','PIDF'))
	dc		<- subset(dc, select=c(SID,PIDF))
	set(dc, NULL, 'SID', dc[, gsub('#','_',SID)])
	set(dc, NULL, 'PID', dc[, gsub('-S[0-9]+$','',PIDF)])
	#
	#	list of files that Chris has
	#
	infiles	<- data.table(F=list.files(indir, pattern='^CW_PANGEA',full.names=TRUE))
	tmp		<- do.call('rbind',lapply(infiles[, F], function(file){
						tmp<- as.data.table(read.csv(file, header=FALSE, col.names='PID', stringsAsFactors=FALSE))
						tmp[, F:=file]
						tmp
					}))
	set(tmp, NULL, 'PROC_STATUS', tmp[, gsub('CW_PANGEA_Rakai_','',gsub('\\.txt','',basename(F)))])	
	set(tmp, NULL, 'F', NULL)
	#	all files for which we have some PID
	dc		<- merge(dc, tmp, by='PID',all=1)	
	stopifnot(!nrow(subset(dc, is.na(SID) & PROC_STATUS!='ThoseWithoutFastqs')))
	dc		<- subset(dc, is.na(SID) | SID!='15351_1_1')		#	remove 15351_1_1  with no PANGEA info in WTSI file
	dc		<- subset(dc, is.na(SID) | SID!='15430_1_75')	#	remove 15430_1_75  with no PANGEA info in WTSI file
	stopifnot(!nrow(subset(dc, is.na(PIDF) & !is.na(SID))))
	dc[, PART:= as.numeric(gsub('^[0-9]+_([0-9])_[0-9]+','\\1',SID))]
	dc[, DUMMY:= gsub('^([0-9]+)_[0-9]_([0-9]+)','\\1_x_\\2',SID)]	
	dc		<- merge(dc, dc[, list(N_PART=length(PART), S_PART=sum(PART)), by='DUMMY'], by='DUMMY')	
	#	check if we always have _1_ and _2_
	stopifnot(!nrow(subset(dc, N_PART==2 & S_PART!=3)))
	#	assume Tanya merges to _3_
	tmp		<- dc[, which(N_PART==2)]
	set(dc, tmp, 'SID', dc[tmp, gsub('^([0-9]+)_[0-9]_([0-9]+)','\\1_3_\\2',SID)])
	set(dc, NULL, c('PART','N_PART','S_PART','DUMMY'),NULL)
	dc		<- unique(dc)
	tmp		<- subset(dc, !is.na(SID))[, list(N_SID=length(SID)), by='PIDF']
	dc		<- merge(dc, tmp, by='PIDF',all.x=1)
	#	define controls
	set(dc, dc[, which(grepl('neg',PID))], 'PROC_STATUS', 'NegControl')
	#	define not in Chris census
	set(dc, dc[, which(is.na(PROC_STATUS))], 'PROC_STATUS', 'NotTrackedByChris')	
	#	add extra category to PROC_STATUS: not processed by Kate
	tmp		<- as.data.table(read.table("~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/KG_PANGEA_Processed_597.txt", header=TRUE,stringsAsFactors=FALSE))
	setnames(tmp, c('SampleID','LaneID'), c('PID','SID'))
	tmp[, KATE_PROC:='Y']
	dc		<- merge(dc, tmp, by=c('PID','SID'), all.x=1)	
	set(dc, dc[, which(PROC_STATUS=='ThoseWithFastqs_WithKateShiverOutput' & is.na(KATE_PROC))], 'PROC_STATUS','ThoseWithFastqs_KateNotProcessed') 
	set(dc, NULL, 'KATE_PROC', NULL)	
	#	add RIDs
	load('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda')
	tmp		<- subset(as.data.table(rccsData), select=c(RCCS_studyid,Pangea.id, batch, date, SEX))
	setnames(tmp, c('RCCS_studyid','Pangea.id','batch','date'), c('RID','PID','RCCS_SHIP_BATCH','SAMPLE_DATE'))
	tmp		<- subset(tmp, !is.na(PID))
	tmp		<- unique(tmp, by=c('RID','PID'))
	tmp2	<- subset(as.data.table(neuroData), select=c(studyid, Pangea.id, sampleDate, gender))
	setnames(tmp2, c('studyid','Pangea.id','sampleDate','gender'), c('RID','PID','SAMPLE_DATE','SEX'))
	tmp2[, RCCS_SHIP_BATCH:='neuro']
	tmp		<- rbind(tmp, tmp2, use.names=TRUE)
	set(tmp, NULL, 'RID', tmp[, as.character(RID)])
	set(tmp, NULL, 'PID', tmp[, as.character(PID)])
	set(tmp, NULL, 'SEX', tmp[, as.character(SEX)])
	dc		<- merge(dc, tmp, by='PID',all.x=1)
	#	flag test plate
	set(dc, dc[, which(grepl('PG14-UG9000[0-9][0-9]',PID))], 'RCCS_SHIP_BATCH', 'test')
	#	see if on HPC	
	tmp		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/listOfBams170704.txt', header=FALSE, col.names='SID', stringsAsFactors=FALSE))
	tmp[, HPC_BAM:='Y']
	set(tmp, NULL, 'SID', tmp[, gsub('\\.bam','',SID)])
	tmp		<- subset(tmp, grepl('^[0-9]+_[0-9]_[0-9]+',SID))	#	reduce to bams from SANGER
	stopifnot( !length(setdiff( tmp[, SID], dc[, SID] )) )
	dc		<- merge(dc, tmp, by='SID',all.x=TRUE)
	set(dc, dc[, which(is.na(HPC_BAM))], 'HPC_BAM', 'N')	
	tmp		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/listOfRefs170704.txt', header=FALSE, col.names='SID', stringsAsFactors=FALSE))
	tmp[, HPC_REF:='Y']
	set(tmp, NULL, 'SID', tmp[, gsub('_ref.fasta','',SID)])
	tmp		<- subset(tmp, grepl('^[0-9]+_[0-9]_[0-9]+',SID))	#	reduce to refs from SANGER	
	stopifnot( !length(setdiff( tmp[, SID], dc[, SID] )) )	
	dc		<- merge(dc, tmp, by='SID',all.x=TRUE)
	set(dc, dc[, which(is.na(HPC_REF))], 'HPC_REF', 'N')
	tmp		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/HPC_census_fastq.txt', header=FALSE, col.names='SID', stringsAsFactors=FALSE))
	tmp[, HPC_FASTQ:='Y']
	set(tmp, NULL, 'SID', tmp[, gsub('_[0-9].fastq.gz','',SID)])	
	tmp		<- subset(tmp, grepl('^[0-9]+_[0-9]_[0-9]+',SID))	#	reduce to fastqz's from SANGER	
	dc		<- merge(dc, unique(tmp), by='SID',all.x=TRUE)
	set(dc, dc[, which(is.na(HPC_FASTQ))], 'HPC_FASTQ', 'N')	
	#	add latest PANGEA stats
	tmp		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/PANGEA_data/2016-07-07_PANGEA_stats_by_sample.csv', stringsAsFactors=FALSE))
	setnames(tmp, 	c('Status','Submitted','DayDiff','ProjectID','Cohort','Submitted.1','Sequenced','Assembled','HIVcontig'), 
			c('WTSI_STATUS','WTSI_SUBMITTED_DATE','DayDiff','PID','Cohort','WTSI_SUBMITTED','WTSI_SEQUENCED','WTSI_ASSEMBLED','WTSI_HIVCONTIG'))
	tmp		<- subset(tmp, PID!='')
	set(tmp, NULL, 'WTSI_SUBMITTED', tmp[, as.character(factor(WTSI_SUBMITTED=='',levels=c(TRUE,FALSE),labels=c('N','Y')))])
	set(tmp, NULL, 'WTSI_SEQUENCED', tmp[, as.character(factor(WTSI_SEQUENCED=='',levels=c(TRUE,FALSE),labels=c('N','Y')))])
	set(tmp, NULL, 'WTSI_ASSEMBLED', tmp[, as.character(factor(WTSI_ASSEMBLED=='',levels=c(TRUE,FALSE),labels=c('N','Y')))])
	set(tmp, NULL, 'WTSI_HIVCONTIG', tmp[, as.character(factor(WTSI_HIVCONTIG=='',levels=c(TRUE,FALSE),labels=c('N','Y')))])
	set(tmp, NULL, c('DayDiff','Cohort'), NULL)
	dc		<- merge(dc, tmp, by='PID', all.x=1)
	#	check what Dan assembled HISEQ
	tmp		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/PANGEA_data/PANGEA_UCL_Feb2017_collated_stats_all_genomes_UCL_release_Feb2017_allhitodate.csv', stringsAsFactors=FALSE))
	setnames(tmp, c('PANGEA_ID','PGID_full','WTSI_ID','Length', 'Cohort'), c('PID','PIDF','SID','UCL_LEN', 'COHORT'))
	tmp		<- subset(tmp, select=c('PID','PIDF','SID','UCL_LEN', 'COHORT'))
	#	convert SID to _3_ to match our convention
	set(tmp, NULL, 'SID', tmp[, gsub('^([0-9]+)_[0-9]_([0-9]+)','\\1_3_\\2',SID)])
	#	check what Dan assembled MISEQ
	tmp2	<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/PANGEA_data/PANGEA_UCL_Feb2017_collated_stats_all_genomes_UCL_release_Feb2017_allmitodate.csv', stringsAsFactors=FALSE))
	setnames(tmp2, c('PANGEA_ID','PGID_full','WTSI_lane_ID','Length', 'Cohort'), c('PID','PIDF','SID','UCL_LEN', 'COHORT'))
	tmp2	<- subset(tmp2, select=c('PID','PIDF','SID','UCL_LEN', 'COHORT'))	
	tmp		<- rbind(tmp, tmp2)
	stopifnot( !length(setdiff(tmp[, SID], dc[, SID])) )	
	#	n=6 plus test plate
	if(0)
	{
		tmp		<- subset(tmp, grepl('Rakai',COHORT))
		tmp2	<- setdiff( tmp[, unique(sort(PID))], dc[, unique(sort(PID))] )
		write.csv(data.table(ID=tmp2), file=file.path(indir,'IDs_that_Dan_flags_from_Rakai_but_not_in_Chris_census.csv'), row.names=FALSE)			
	}
	set(tmp, NULL, 'COHORT', NULL)		
	dc		<- merge(dc, tmp, by=c('PID','PIDF','SID'),all.x=1)	
	#	resolve 60 from Chris that he confirmed he has not processed -- at least one SID from said individual has been processed
	set(dc, dc[, which(PROC_STATUS=='ThoseWithFastqs_WithChrisShiverOutput' & HPC_BAM=='N')], 'PROC_STATUS','ThoseWithFastqs_ChrisNotProcessed')
	#
	#	subset to RAKAI
	#
	dc		<- subset(dc, !is.na(RID))
	#	define RIDs from whom all SIDs are complete and on HPC
	tmp		<- dc[, list(	HPC_ALL_SID_FOR_RID= all(HPC_BAM=='Y') & all(HPC_REF=='Y')), by='RID']	
	set(tmp, NULL, 'HPC_ALL_SID_FOR_RID', tmp[,as.character(factor(HPC_ALL_SID_FOR_RID, levels=c(TRUE,FALSE), labels=c('Y','N')))])
	dc		<- merge(dc, tmp, by='RID',all.x=1)
	stopifnot(!nrow(subset(dc, HPC_BAM!=HPC_REF)))
	#	define Sampling Time
	set(dc, NULL, 'WTSI_SUBMITTED_DATE', dc[, as.Date(WTSI_SUBMITTED_DATE)])		
	#	add column consensus available + length consensus
	tmp		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_Imperial_v170704_ALL_bestcov_cov1.fasta'
	sq		<- read.dna(tmp, format='fasta')	
	sqi		<- data.table(PIDF=rownames(sq))	
	tmp		<- sapply(seq_len(nrow(sq)), function(i) base.freq(sq[i,], all=TRUE, freq=TRUE))
	sqi[, IMPERIAL_LEN:=ncol(sq)-apply( tmp[c('-','?'),], 2, sum)]
	set(sqi, NULL, 'PIDF', sqi[,gsub('_WTSI_*','',PIDF)])
	dc		<- merge(dc, sqi, by='PIDF', all.x=1)	
	save(dc, file=file.path(indir,'Rakai_phyloscanner_170704_assemblystatus.rda'))
	
	#
	#	set up second batch
	#
	if(1)
	{
		dc.new	<- copy(dc)
		#	update batches
		load( file.path(indir,'Rakai_phyloscanner_170301_b75.rda') )		
		stopifnot(!length(setdiff(dc[, RID], dc.new[, RID])))
		batch.n	<- 75
		batch.n2<- 55	# use smaller number because we typically have multiple seqs from these individuals
		dc.it2	<- subset(dc.new, HPC_BAM=='Y' & HPC_REF=='Y')
		tmp		<- unique(subset(dc.it2, select=RID))	
		tmp		<- merge(tmp, unique(subset(dc.it1, select=c(RID, RID_B, BATCH))), by='RID', all.x=1)
		tmp2	<- tmp[, which(is.na(RID_B))]
		set.seed(42)
		tmp3	<- tmp[, ceiling(max(RID_B, na.rm=TRUE) / batch.n)*batch.n]
		set(tmp, tmp2, 'RID_B', as.integer(sample(length(tmp2),length(tmp2)) + tmp3) )
		dc.it2	<- merge(dc.it2, tmp, by='RID')
		setkey(dc.it2, RID_B)
		tmp2	<- dc.it2[, which(is.na(BATCH))]
		tmp3	<- min( dc.it2[tmp2, ceiling(RID_B/batch.n2)] ) - dc.it1[, max(BATCH)] - 1L		
		set(dc.it2, tmp2, 'BATCH', dc.it2[tmp2, ceiling(RID_B/batch.n2)] - tmp3)
		set(dc.it2, dc.it2[, which(BATCH==max(BATCH))], 'BATCH', dc.it2[, max(BATCH)-1L])		
		#
		#	set up remaining pty.runs
		#		
		tmp		<- as.data.table(t(combn(dc.it2[, unique(BATCH)],2)))
		setnames(tmp, c('V1','V2'), c('BATCH','BATCH2'))
		pty.runs<- unique(subset(pty.runs, select=c('BATCH','PTY_RUN')))
		pty.runs<- pty.runs[, list(BATCH=BATCH[1], BATCH2=BATCH[2]), by='PTY_RUN']		
		pty.runs<- merge(tmp, pty.runs, by=c('BATCH','BATCH2'), all.x=1)
		tmp		<- pty.runs[, which(is.na(PTY_RUN))]
		set(pty.runs, tmp, 'PTY_RUN', seq_along(tmp)+pty.runs[, max(PTY_RUN, na.rm=1)] )
		#	select only those that are not yet processed
		tmp			<- dc.it1[, max(BATCH)]
		pty.runs	<- subset(pty.runs, BATCH>tmp | BATCH2>tmp)		
		
		pty.runs	<- melt(pty.runs, id.vars='PTY_RUN', variable.name='DUMMY', value.name='BATCH')
		set(pty.runs, NULL, 'DUMMY', NULL)
		setkey(pty.runs, PTY_RUN)
		tmp			<- subset(dc.it2, select=c(BATCH, RID, SID))
		tmp			<- tmp[, list(SID=SID, RENAME_SID=paste0(RID,'_fq',seq_along(SID))), by=c('BATCH','RID')]
		pty.runs	<- merge(pty.runs, tmp, by='BATCH',allow.cartesian=TRUE)
		#	save
		dc			<- copy(dc.new)
		save(dc, dc.it2, pty.runs, file=file.path(indir,'Rakai_phyloscanner_170704_b75_part2.rda'))		
	}
	
}

project.dualinfecions.UG.setup.coinfections.160219<- function()
{
	#
	#	input args
	#	
	pty.sel.n		<- 25
	indir		 	<- "~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/data"
	outdir			<- indir
	infile.bam		<- file.path(outdir,'PANGEA_HIV_n5003_Imperial_v160110_UG_selectedbycoverage_160219.rda')
	infile.gagtree	<- file.path(outdir,'PANGEA_HIV_n5003_Imperial_v160110_UG_gag_fasttree.rda')
	infile.endtree	<- file.path(outdir,'PANGEA_HIV_n5003_Imperial_v160110_UG_p15toend_fasttree.rda')
	sqi
	
	load(infile.bam)
	load(infile.gagtree)
	#
	#	set groups 	
	#
	set(bam.inds, bam.inds[, which(GROUP=='2000-9000')],'GROUP','3900-9000')	
	set(bam.inds, bam.inds[, which(GROUP=='50-9000')],'GROUP','1-9000')
	set(bam.inds, bam.inds[, which(GROUP=='50-3900')],'GROUP','1-3900')
	set(bam.inds, bam.inds[, which(GROUP=='50-2000')],'GROUP','1-2000')
	set(bam.inds, bam.inds[, which(GROUP=='650-3900')],'GROUP','650-2000')
	bam.inds	<- merge(bam.inds, subset(sqi, select=c(TAXA, SANGER_ID)), by='SANGER_ID')	
	#	get distance data.table
	ph.gdf		<- as.data.table(melt(ph.gdtr))
	setnames(ph.gdf, c('X1','X2','value'),c('TAXA','TAXA2','GD'))
	tmp			<- data.table(TX_IDX=seq_len(Ntip(ph)), TAXA= ph$tip.label)
	tmp			<- merge(tmp, sqi,by='TAXA')
	tmp			<- subset(tmp, SITE=='UG')
	tmp[, DUMMY:=NULL]
	ph.gdf		<- merge(ph.gdf, subset(tmp, select=TAXA), by='TAXA')
	ph.gdf		<- merge(ph.gdf, data.table(TAXA2=tmp[,TAXA]), by='TAXA2')
	ph.gdf		<- subset(ph.gdf, TAXA!=TAXA2)
	setkey(ph.gdf, TAXA, GD)	
	#
	#	select very close individuals
	#
	ph			<- drop.tip(ph,setdiff(which(!grepl('^PG',ph$tip.label)),which(grepl('AF103818',ph$tip.label))))
	#	TODO the cascade clustering does NOT work
	#dist.brl	<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)
	dist.brl	<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=max)	
	thresh.brl	<- 0.05
	thresh.bs	<- 0.01				
	clustering	<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph$node.label,retval="all")	
	table(clustering$size.tips)
	file		<- paste(gsub('\\.rda',paste('_gd',100*thresh.brl,'bs',thresh.bs*100,'\\.pdf',sep=''), infile.gagtree), sep='')
	invisible(hivc.clu.plot(ph, clustering[["clu.mem"]], cex.edge.incluster=3, tip.color=tip.color, file=file, pdf.scaley=100, show.tip.label=TRUE, pdf.width=30))
	pty.clu		<- subset(data.table(TX_IDX=seq_len(Ntip(ph)), TAXA= ph$tip.label, CLU_ID=clustering$clu.mem[seq_len(Ntip(ph))]), !is.na(CLU_ID))
	pty.clu[, PANGEA:= factor(grepl('^PG14',TAXA), levels=c(TRUE,FALSE),labels=c('Y','N'))]
	pty.clu		<- merge(pty.clu, bam.inds, by='TAXA')	
	pty.clu		<- merge(pty.clu, pty.clu[, list(CLU_N= length(TAXA)), by='CLU_ID'], by='CLU_ID')	
	#	initialize pty.runs: for every cluster, go down ancestor path until we get a clade of at most x individuals
	setkey(pty.clu, CLU_ID)	
	tmp			<- unique(pty.clu)
	tmp			<- tmp[order(-CLU_N),]		
	pty.clu		<- tmp[, {
				clu.a	<- Ancestors(ph, TX_IDX)
				clu.n	<- sapply( Descendants(ph, clu.a, type="tips"), length)
				z		<- length(clu.n[clu.n<25])
				ans		<- NA_integer_
				if(z>0)
					ans	<- Descendants(ph, clu.a[z])[[1]]
				list(IDX=ans)				
			}, by=c('CLU_ID','TX_IDX')]
	pty.clu[, PTY_RUN:=CLU_ID]
	#	reduce pty.runs: keep only clades that are not contained in other clades
	for(x in pty.clu[, unique(TX_IDX)])
	{		
		tmp			<- merge(pty.clu, subset(pty.clu, IDX==x, PTY_RUN), by='PTY_RUN')
		if( all(subset(tmp, TX_IDX==x)[,IDX]%in%subset(tmp, TX_IDX!=x)[,IDX]) )
			pty.clu		<- subset(pty.clu, TX_IDX!=x)		
	}
	pty.clu[, table(PTY_RUN)]
	#	add taxa not yet included: for each taxa, determine closest runs
	ph.gdf		<- merge(ph.gdf, data.table(TAXA=ph$tip.label, IDX=seq_along(ph$tip.label)), by='TAXA')
	ph.gdf		<- merge(ph.gdf, subset(pty.clu, select=c(IDX, PTY_RUN)), by='IDX', all.x=1)
	tmp			<- data.table(IDX= setdiff(seq_len(Ntip(ph)), pty.clu[, IDX]))	
	tmp			<- tmp[, {				
				z	<- ph$tip.label[IDX]
				z	<- subset(ph.gdf, !is.na(PTY_RUN) & TAXA2==z)				
				setkey(z, GD)
				list(PTY_RUN=z[1:100,][, unique(PTY_RUN)][1:10], ORDER=1:10)				
			}, by='IDX']	
	tmp			<- subset(tmp, !is.na(PTY_RUN))	#this excludes the root

	#	add taxa not yet included: determine closest run that is not yet full and add
	for(x in tmp[, unique(IDX)])
	{
		z		<- merge(pty.clu, subset(tmp,IDX==x,c(PTY_RUN,ORDER)), by='PTY_RUN')
		z		<- subset(z[, list(N=length(IDX),ORDER=ORDER[1]), by='PTY_RUN'], N<20)
		setkey(z, ORDER)
		#cat('\n',x,z[1,PTY_RUN])
		pty.clu	<- rbind(pty.clu, data.table(IDX=x, PTY_RUN=z[1,PTY_RUN]), use.names=TRUE, fill=TRUE)
	}
	pty.clu[, table(PTY_RUN, useNA='if')]
	#
	#	prepare to plot coverage by current run
	#
	pty.runs	<- copy(pty.clu)
	pty.runs[, TAXA:= ph$tip.label[IDX]]
	pty.runs	<- merge(pty.runs, subset(sqi, select=c(TAXA, SANGER_ID)), by='TAXA')
	setnames(pty.runs, 'SANGER_ID','FILE_ID')
	save(pty.runs, file=gsub('_fasttree\\.rda','_selecthelp.rda',infile.gagtree))
	#
	#	OK, just merge closest groups
	#
	ph.gdtrn	<- dist.nodes(ph)
	pty.clu		<- merge(pty.clu, pty.clu[, list(N=length(IDX), MRCA=getMRCA(ph, IDX)), by='PTY_RUN'], by='PTY_RUN')	
	setkey(pty.clu, PTY_RUN)
	pty.small	<- subset(unique(pty.clu), N<20, c(PTY_RUN, MRCA, N))
	pty.small[, GD:=ph.gdtrn[pty.small[1, as.character(MRCA)], pty.small[, as.character(MRCA)]]]
	pty.small[, N_NEW:= c(pty.small[1,N],pty.small[1,N]+pty.small[-1,N])]
	pty.small	<- subset(pty.small, N_NEW<pty.sel.n)	
	setkey(pty.small, GD)
	while(nrow(pty.small)>1)
	{
		cat('\nmerge',pty.small[1,PTY_RUN],pty.small[2,PTY_RUN])
		set(pty.clu, pty.clu[, which(PTY_RUN==pty.small$PTY_RUN[2])], 'PTY_RUN', pty.small$PTY_RUN[1])
		set(pty.clu, pty.clu[, which(PTY_RUN%in%pty.small$PTY_RUN[1:2])], 'N', pty.small$N_NEW[2])
		setkey(pty.clu, PTY_RUN)
		pty.small	<- subset(unique(pty.clu), N<20, c(PTY_RUN, MRCA, N))
		pty.small[, GD:=ph.gdtrn[pty.small[1, as.character(MRCA)], pty.small[, as.character(MRCA)]]]
		pty.small[, N_NEW:= c(pty.small[1,N],pty.small[1,N]+pty.small[-1,N])]
		pty.small	<- subset(pty.small, N_NEW<pty.sel.n)	
		setkey(pty.small, GD)		
	}
	#	clean up
	set(pty.clu, NULL, c('TX_IDX','N','MRCA'), NULL)
	set(pty.clu, NULL, 'PTY_RUN', pty.clu[, as.integer(factor(as.character(PTY_RUN)))])
	setkey(pty.clu, PTY_RUN)
	pty.clu[, TAXA:=ph$tip.label[IDX]]
	pty.runs	<- merge(pty.clu, subset(sqi, select=c(TAXA, SANGER_ID)), by='TAXA')
	setnames(pty.runs, 'SANGER_ID','FILE_ID')
	#	save
	outfile		<- gsub('fasttree','coinfinput_160825',infile.gagtree)
	save(pty.runs, ph, dist.brl, ph.gdtr, ph.mrca, clustering, sqi, file=outfile)	
	#	plot runs on tree
	phd			<- data.table(IDX=seq_along(ph$tip.label))
	phd			<- merge(phd, pty.runs, by='IDX', all.x=1)
	phd			<- merge(phd, subset(bam.inds, select=c(TAXA, GROUP)), by='TAXA', all.x=1)
	phd[, TN:= paste(PTY_RUN,'-CLU-',CLU_ID,'-GR-',GROUP,'-IDX-',IDX,sep='')]	
	tmp			<- data.table(PTY_RUN=phd[, unique(PTY_RUN)])
	tmp[, PTY_COL:= rainbow(nrow(tmp))]
	phd			<- merge(phd, tmp, by='PTY_RUN')
	setkey(phd, IDX)
	tmp				<- copy(ph)
	tmp$tip.label	<- phd[, TN]	
	tip.color		<- phd[, PTY_COL]
	write.tree(tmp, file=gsub('.rda','_tree.newick',outfile))
	invisible(hivc.clu.plot(tmp, clustering[["clu.mem"]], cex.edge.incluster=3, tip.color=tip.color, file=gsub('.rda','_tree.pdf',outfile), pdf.scaley=100, show.tip.label=TRUE, pdf.width=30))
	#	..looking good!	
}

project.dualinfecions.UG.setup.coinfections.160217<- function()
{
	#
	#	input args
	#	
	thresh.brl	<- 0.1
	thresh.bs	<- 0.7			
	pty.sel.n	<- 25
	
	infile		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/PANGEA_HIV_n5003_Imperial_v160110_GlobalAlignment.rda'
	load(infile)	#loads sqi, sq	
	indir		<- "~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/data"
	outdir		<- indir
	infile		<-  "PANGEA_HIV_n5003_Imperial_v160110_UG_NoQ_fast2.rda"
	load(file.path(indir,infile))	#loads "ph" "dist.brl" "ph.gdtr"  "ph.mrca"	"clustering"
	#	select Uganda samples
	sqi			<- subset(sqi, SITE=='UG' | SEQLOC=='LosAlamos')
	#	add SANGER_ID
	infile.s	<- "~/Dropbox (SPH Imperial College)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	si			<- as.data.table(read.csv(infile.s, stringsAsFactors=FALSE))
	setnames(si, colnames(si), toupper(gsub('.','_',colnames(si),fixed=1))) 
	set(si, NULL, 'PANGEA_ID', si[, gsub(' ','',PANGEA_ID)])
	setnames(si, 'CLINICAL_GENOME_COVERAGE', 'COV')
	tmp			<- subset(si, grepl('^PG14-UG', PANGEA_ID), c(PANGEA_ID, SANGER_ID))
	set(tmp, NULL, 'PANGEA_ID', tmp[,gsub('-','_',PANGEA_ID)])
	#	of the duplicate PANGEA_IDs, consider only those with larger coverage	
	sqi			<- merge(sqi, tmp, by='PANGEA_ID', all.x=1, allow.cartesian=TRUE)
	tmp			<- sqi[, list(SANGER_ID=SANGER_ID[which.max(COV)]), by='PANGEA_ID']
	sqi			<- merge(sqi, tmp, by=c('PANGEA_ID','SANGER_ID'))
	#
	#	select closest 15 individuals
	#
	clustering	<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=thresh.bs, thresh.brl=thresh.brl, dist.brl=dist.brl, nodesupport=ph$node.label,retval="all")	
	table(clustering$size.tips)
	
	pty.clu		<- subset(data.table(TX_IDX=seq_len(Ntip(ph)), TAXA= ph$tip.label, CLU_ID=clustering$clu.mem[seq_len(Ntip(ph))]), !is.na(CLU_ID))	
	#	reduce to clusters containing at least one ZA sequence
	pty.clu		<- merge(pty.clu, sqi,by='TAXA')
	tmp			<- pty.clu[, list(ANY_NOT_INCOUNTRY= any(is.na(SITE) | SITE!='UG'), ALL_NOT_INCOUNTRY= all(is.na(SITE) | SITE!='UG')), by='CLU_ID']
	cat('\nInspecting clusters if not in-country')
	print(tmp[, table(ANY_NOT_INCOUNTRY, ALL_NOT_INCOUNTRY)])
	tmp			<- subset(tmp, !ALL_NOT_INCOUNTRY)
	pty.clu		<- merge(pty.clu, subset(tmp, select=c(CLU_ID)), by='CLU_ID')
	#	reduce to clusters of ZA sequences
	pty.clu		<- subset(pty.clu, PNG=='Y')
	pty.clu		<- merge(pty.clu, pty.clu[, list(CLU_N= length(TAXA)), by='CLU_ID'], by='CLU_ID')
	pty.clu		<- subset(pty.clu, CLU_N>1)
	cat('\nFound in-country clusters of size:')	
	print( unique(subset(pty.clu, select=c(CLU_ID, CLU_N)))[, table(CLU_N)] )
	#	get distance data.table
	ph.gdf		<- as.data.table(melt(ph.gdtr,value.name="GD"))
	setnames(ph.gdf, c('Var1','Var2'),c('TAXA','TAXA2'))
	tmp			<- data.table(TX_IDX=seq_len(Ntip(ph)), TAXA= ph$tip.label)
	tmp			<- merge(tmp, sqi,by='TAXA')
	tmp			<- subset(tmp, SITE=='UG')
	tmp[, DUMMY:=NULL]
	ph.gdf		<- merge(ph.gdf, subset(tmp, select=TAXA), by='TAXA')
	ph.gdf		<- merge(ph.gdf, data.table(TAXA2=tmp[,TAXA]), by='TAXA2')
	ph.gdf		<- subset(ph.gdf, TAXA!=TAXA2)
	setkey(ph.gdf, TAXA, GD)	
	#	initialize pty.runs
	setkey(pty.clu, CLU_ID)	
	tmp			<- unique(pty.clu)
	tmp			<- tmp[order(-CLU_N),]	
	tx.seeds	<- c(tmp[, TAXA], setdiff(subset(sqi, SITE=='UG')[,TAXA], tmp[, TAXA]))	
	#	not all seeds in tree (ie no consensus)
	tx.seeds	<- intersect( tx.seeds, ph$tip.label )
	cat('\ninitial seeds n=',length(tx.seeds))
	ptyr		<- 1L	
	pty.runs	<- data.table(PTY_RUN=NA_integer_, TAXA=NA_character_, GD=NA_real_)
	#	fill pty.runs
	while(length(tx.seeds))
	{
		seed		<- tx.seeds[1]
		tmp			<- data.table(TAXA=c(seed, subset(ph.gdf, TAXA==seed )[seq_len(pty.sel.n-1), TAXA2]), GD=c(0,subset(ph.gdf, TAXA==seed )[seq_len(pty.sel.n-1), GD]))
		tmp[, PTY_RUN:=ptyr]
		tx.seeds	<- setdiff( tx.seeds, tmp[,TAXA] )
		ph.gdf		<- subset(ph.gdf, !TAXA%in%tmp[,TAXA] & !TAXA2%in%tmp[,TAXA]) 
		cat('\nrun',ptyr,'select n=',nrow(tmp),'remaining seeds n=',length(tx.seeds), 'remaining distances n=', nrow(ph.gdf))
		pty.runs	<- rbind(tmp, pty.runs)
		ptyr		<- ptyr+1L	
	}
	pty.runs	<- subset(pty.runs, !is.na(TAXA))
	#pty.runs[, table(PTY_RUN)]
	#	NOTE!!! merge last two since last v small
	set(pty.runs, pty.runs[, which(PTY_RUN==146L)], 'PTY_RUN', 145L)
	pty.runs[, FILL:=0L]
	pty.runs[, TX_IDX:=match(pty.runs[,TAXA],ph$tip.label)]
	pty.runs[, GD:=NULL]
	#
	tmp			<- subset(sqi, select=c(PANGEA_ID, SANGER_ID))
	setnames(tmp, c('PANGEA_ID','SANGER_ID'), c('TAXA','FILE_ID') )
	pty.runs	<- merge(pty.runs, tmp, by='TAXA')	
	setkey(pty.runs, PTY_RUN)	
	#	
	cat('\nNumber of clusters=', pty.runs[, length(unique(CLU_ID))])
	cat('\nNumber of scheduled phylotype runs=', pty.runs[, max(PTY_RUN)])
	cat('\nNumber of selected taxa=', subset(pty.runs, !FILL)[, length(unique(TAXA))])	
	tmp			<- paste(indir, '/', gsub('\\.rda','_coinfrunsinput\\.rda',infile), sep='')
	save(pty.runs, pty.clu, ph, dist.brl, ph.gdtr, ph.mrca, clustering, sqi, sq, file= tmp)	
}

project.dualinfecions.phylotypes.test<- function()
{
	tfd	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/phylotypes'
	tfd	<- '~/duke/2016_PANGEAphylotypes/phylotypes'
	tf	<- data.table(FILE=list.files(tfd, pattern='^ptyr1_*'))
	tf	<- subset(tf, !grepl('fasta2',FILE))
	tf[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]
	tf[, W_FROM:=gsub('InWindow_','',regmatches(FILE,regexpr('InWindow_[0-9]+',FILE)))] 
	tf[, W_TO:=gsub('to_','',regmatches(FILE,regexpr('to_[0-9]+',FILE)))]
	tf[, Q1:=2]
	set(tf, tf[, which(grepl('Q1-30',FILE))], 'Q1', 30)
	invisible(	tf[,{
						system(paste("sed 's/<unknown description>//' ",file.path(tfd,FILE)," > ",file.path(tfd,paste(FILE,'2',sep='')),sep=''))
					}	, by='FILE']	)
	
	
	tmp	<- tf[, list(BAM= rownames(read.dna(file.path(tfd,paste(FILE,'2',sep='')),format='fasta'))), by='FILE']
	tf	<- merge(tf, tmp, by='FILE')	
	tf[, FILE_ID:= gsub('_read.*','',BAM)]	
	tf	<- merge(unique(subset(pty.runs, select=c( TAXA ,FILE_ID))), tf, by='FILE_ID')
	tmp	<- tf[, list(BAM_N=length(BAM)), by=c('W_FROM','Q1','FILE_ID')]
	tmp	<- dcast.data.table(tmp, W_FROM+FILE_ID~Q1, value.var='BAM_N')
	setnames(tmp, c('2','30'), c('Q1_2','Q1_30'))
	subset(tmp, Q1_2<Q1_30)
	
	ggplot(tmp, aes(x=as.numeric(W_FROM),y=BAM_N,colour=FILE_ID,group=FILE_ID)) + geom_line() + facet_grid(~Q1)
	#no R8_RES486_S8_L001 at all!

	load("/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_ptyrunsinput.rda")
}

project.test.BEEHIVEtree<- function()
{
	indir				<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/beehive_test'
	indir				<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/beehive_test/JuicyClade1_NoMerge'
	indir				<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/beehive_test/test2'
	pty.runs			<- NULL 
	outdir				<- indir
	select				<- ''
	outgroup			<- 'HXB2'
	references.pattern	<- 'B\\.|C\\.'
	run.pattern			<- ''
	rm.newick			<- FALSE
	rm.fasta			<- FALSE	
	tree.pattern		<- 'tree$'
	plot.trees.per.page	<- 1
	plot.w				<- 20
	plot.h				<- 40
	pty.evaluate.tree(indir, pty.runs=pty.runs, outdir=indir, select=select, outgroup=outgroup, references.pattern=references.pattern, run.pattern=run.pattern, rm.newick=rm.newick, rm.fasta=rm.fasta)
}
	
project.scan.superinfections.160211	<- function()
{
	indir		<- file.path(HOME,'coinf_ptoutput_150201')
	infiles		<- data.table(FILE=list.files(indir, pattern='_stat.rda$'))
	infiles		<- subset(infiles, grepl('^ptyr[0-9]+_',FILE))
	infiles[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]
	setkey(infiles, PTY_RUN)
	stat.clades	<- infiles[, {
				load(file.path(indir,FILE))
				stat.clades
			}, by='PTY_RUN']
	stat.clades[, PTY_RUN:=NULL]
	stat.ind	<- infiles[, {
				load(file.path(indir,FILE))
				stat.ind
			}, by='PTY_RUN']
	stat.ind[, PTY_RUN:=NULL]
	#	plot all patients
	outfile		<- file.path(indir,'scan_superinfections_160211_all.pdf')
	pty.stat.superinfections.160208.plot(stat.clades, stat.ind, outfile, plot.max.clade=5)
	#	previously identified cases
	#	
	prev			<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/coinf_ptoutput_150121/superinfections_clcountn_100_taxan_40_diffind_1_winn_3_curated.csv'))
	setnames(prev, 'FILE_ID', 'IND')
	prev.super		<- subset(prev, grepl('^Super|^super',THOUGHTS_OR), c(IND, PTY_RUN, WINDOW_N))	
	stat.clades.spr	<- merge(stat.clades, subset(prev.super, select=IND), by='IND')
	stat.ind.spr	<- merge(stat.ind, subset(prev.super, select=IND), by='IND')
	outfile			<- file.path(indir,'scan_superinfections_160211_super.pdf')
	pty.stat.superinfections.160208.plot(stat.clades.spr, stat.ind.spr, outfile, plot.max.clade=5)
	#	plot unclear prev candidates
	tmp				<- unique(as.character(prev.super$IND))
	set(prev, NULL, 'IND', prev[, as.character(IND)])
	prev.unclear	<- subset(prev, !IND%in%tmp, c(IND, PTY_RUN, WINDOW_N))
	
	stat.clades.u	<- merge(stat.clades, subset(prev.unclear, select=IND), by='IND')
	stat.ind.u		<- merge(stat.ind, subset(prev.unclear, select=IND), by='IND')
	outfile			<- file.path(indir,'scan_superinfections_160211_unclear.pdf')
	pty.stat.superinfections.160208.plot(stat.clades.u, stat.ind.u, outfile, plot.max.clade=5)
	
}

project.scan.superinfections.160203	<- function()
{
	#
	# input files
	#stat.infile		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/coinf_ptoutput_150121/ptyr_examl_stat.rda'
	stat.infile	<- file.path(HOME,'coinf_ptoutput_150121/ptyr_examl_stat.rda')
	tree.indir	<- file.path(HOME,'coinf_ptoutput_150121')
	#	should have strong evidence for superinfections: separating individuals
	if(0)	
	{
		coi.args<- list(clcountn=100, taxan=40, diffind=1, winn=3)
	}
	#	weaker evidence for superinfections: no separating individuals
	if(1)	
	{
		coi.args<- list(clcountn=100, taxan=40, diffind=0, winn=3)
	}
	
	tmp			<- paste(sapply(seq_along(coi.args),function(i) paste(names(coi.args)[i],coi.args[i],sep='_')	),collapse='_')
	outfile		<- file.path(dirname(stat.infile),paste('superinfections_',tmp,'.csv',sep=''))
	load(stat.infile)	
	#
	#	long branches. focus on subtending clades with more than 100 reads
	pty.lsep	<- subset(pty.stat, CL_COUNT_N>=coi.args[['clcountn']] & TAXA_N>=coi.args[['taxan']])
	setkey(pty.lsep, PTY_RUN, W_FROM, W_TO, FILE_ID)
	pty.lsep	<- unique(pty.lsep)	
	pty.lsep[, THR_WHER:= 0.1]
	set(pty.lsep, pty.lsep[,which(W_FROM>5200)],'THR_WHER',0.3)
	tmp	<- pty.lsep[, list(CL_MX_LOCAL_SEP_p50=quantile(CL_MX_LOCAL_SEP,p=0.5), CL_MX_LOCAL_SEP_p025=quantile(CL_MX_LOCAL_SEP,p=0.025), CL_MX_LOCAL_SEP_p10=quantile(CL_MX_LOCAL_SEP,p=0.1), CL_MX_LOCAL_SEP_p90=mean(CL_MX_LOCAL_SEP,p=0.9), CL_MX_LOCAL_SEP_p95=quantile(CL_MX_LOCAL_SEP,p=0.95) ), by='W_FROM']	
	ggplot(pty.lsep, aes(x=W_FROM)) +  geom_point(aes(y=CL_MX_LOCAL_SEP, colour=CL_MX_LOCAL_SEP, size=COUNT_N), alpha=0.2 ) + 
			geom_ribbon(data=tmp, aes(ymin=CL_MX_LOCAL_SEP_p10, ymax=CL_MX_LOCAL_SEP_p90), fill='grey50', alpha=0.7) +
			geom_ribbon(data=tmp, aes(ymin=CL_MX_LOCAL_SEP_p90, ymax=CL_MX_LOCAL_SEP_p95), fill='grey80', alpha=0.7) +
			#geom_ribbon(data=tmp, aes(ymin=WHDA_p025, ymax=WHDA_p10), fill='grey80', alpha=0.5) +
			scale_colour_continuous(guide=FALSE) + 
			scale_x_continuous(breaks=seq(0,10e3,500)) +
			geom_step(aes(y=THR_WHER), colour='red') +
			geom_step(aes(y=THR_WHER/2), colour='DarkRed') +
			#geom_line(data=tmp, aes(y=CL_MX_LOCAL_SEP_p95), colour='red') +
			#geom_line(data=tmp, aes(y=CL_MX_LOCAL_SEP_p90), colour='DarkRed') +
			geom_line(data=tmp, aes(y=CL_MX_LOCAL_SEP_p50)) + theme_bw() + theme(legend.position='bottom') +
			labs(x='\ngenome position of window start in each run\n(bp)',y='stem length of clades with >100 reads\n(subst/site)\n', size='quality trimmed short reads per individual\n(#)')
	ggsave(file=gsub('\\.csv','_scan\\.pdf',outfile), w=12, h=6)
	#
	#	get candidates
	#tmp		<- subset(pty.lsep, CL_MX_LOCAL_SEP>THR_WHER/2)[, list(PTY_RUN=PTY_RUN[1], FILE_ID_N=length(W_FROM), CL_MX_LOCAL_SEP_avg=mean(CL_MX_LOCAL_SEP)), by='FILE_ID']
	#tmp		<- tmp[order(-FILE_ID_N),]	
	tmp		<- subset(pty.lsep, CL_MX_LOCAL_SEP>THR_WHER/2)[, list(PTY_RUN=PTY_RUN[1], FILE_ID_N=length(W_FROM), DIFF_IND= max(DIFF_IND)), by='FILE_ID']
	if(coi.args[['diffind']]>0)
		tmp	<- subset(tmp, FILE_ID_N>=coi.args[['winn']] & DIFF_IND>=coi.args[['diffind']])
	if(coi.args[['diffind']]==0)
		tmp	<- subset(tmp, FILE_ID_N>=coi.args[['winn']] & DIFF_IND==0)
	tmp		<- tmp[order(-FILE_ID_N),]	
	write.csv(tmp, file=outfile, row.names=FALSE)
	#
	#	plot trees for each window of candidates
	tmp2	<- merge(pty.lsep, subset(tmp, select=FILE_ID), by='FILE_ID')
	tmp2	<- subset(tmp2, CL_MX_LOCAL_SEP>THR_WHER/2 )
	infiles		<- data.table(FILE=list.files(tree.indir, pattern='examl.rda$'))
	infiles[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]		
	phps	<- lapply(seq_len(nrow(tmp2)), function(i)
			{
				#i				<-1
				cat('\nprepare plot',i,'/',nrow(tmp2))
				file			<- file.path(tree.indir, infiles[ PTY_RUN==tmp2[i,PTY_RUN],	FILE])
				load(file)
				ph				<- pty.ph[[ tmp2[i,FILE]  ]]
				max.node.height	<- max(node.depth.edgelength(ph)[1:Ntip(ph)])
				tmp				<- c(c(tmp2[i,FILE_ID],'not characterized'), sort(setdiff(sort(levels(attr(ph,'INDIVIDUAL'))),c(tmp2[i,FILE_ID],'not characterized'))))
				col				<- c('black','grey50',rainbow_hcl(length(tmp)-2, start = 270, end = -30, c=100, l=50))
				names(col)		<- tmp			
				ph.title		<- tmp2[i, paste( FILE_ID,'\nrun=',PTY_RUN,', win=',W_FROM,'-',W_TO,'\ndiffind', DIFF_IND,' diff',DIFF, ' lsep', round(CL_MX_LOCAL_SEP,d=3), sep='')]
				p				<- ggtree(ph, aes(color=INDIVIDUAL, linetype=TYPE)) + 
						geom_nodepoint(size=ph$node.label/100*3) +
						geom_tiplab(size=1.2,  hjust=-.1) +							 
						scale_color_manual(values=col, guide = FALSE) +											 
						scale_linetype_manual(values=c('target'='solid','filler'='dotted'),guide = FALSE) +
						theme_tree2() +
						theme(legend.position="bottom") + ggplot2::xlim(0, max.node.height*1.3) +
						labs(x='subst/site', title=ph.title)
				p
			})	
	tmp			<- seq_len(ceiling(length(phps)/10))
	pdf(file=gsub('\\.csv','_trees\\.pdf',outfile), w=20, h=40)		#for win=60
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

project.scan.contaminants	<- function()
{
	#
	# input files
	#stat.infile		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/coinf_ptoutput_150121/ptyr_examl_stat.rda'
	stat.infile	<- file.path(HOME,'coinf_ptoutput_150121/ptyr_examl_stat.rda')
	tree.indir	<- file.path(HOME,'coinf_ptoutput_150121')
	#	should have strong evidence for superinfections: separating individuals
	if(1)	
	{
		cnt.args<- list(countn=100, taxan=40, diffind=1, diff=10, winn=10, winx=10)
	}
	
	
	tmp			<- paste(sapply(seq_along(cnt.args),function(i) paste(names(cnt.args)[i],cnt.args[i],sep='_')	),collapse='_')
	outfile		<- file.path(dirname(stat.infile),paste('contaminants_',tmp,'.csv',sep=''))
	load(stat.infile)	
	#
	#	select individuals with at least one other individual inbetween and many separate clades of the other individual
	pty.cm		<- subset(pty.stat, DIFF_IND>=cnt.args[['diffind']] & DIFF>cnt.args[['diff']])
	setkey(pty.cm, PTY_RUN, W_FROM, W_TO, FILE_ID)
	pty.cm		<- unique(pty.cm)
	tmp			<- pty.cm[, list(W_N=length(W_FROM)), by='FILE_ID']
	tmp			<- tmp[order(-W_N),]
	tmp			<- subset(tmp, W_N>=cnt.args[['winn']])
	write.csv(tmp, file=outfile, row.names=FALSE)
	#
	#	plot trees for each window of candidates
	tmp2	<- merge(pty.cm, subset(tmp, select=FILE_ID), by='FILE_ID')
	tmp2	<- tmp2[order(FILE_ID, -COUNT_N),]
	#	select at max winx trees per individual
	tmp2	<- merge(tmp2, tmp2[, list(W_FROM=W_FROM[1:min(cnt.args[['winx']],length(W_FROM))]), by='FILE_ID'], by=c('FILE_ID','W_FROM'))	
	infiles		<- data.table(FILE=list.files(tree.indir, pattern='examl.rda$'))
	infiles[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]		
	phps	<- lapply(seq_len(nrow(tmp2)), function(i)
			{
				#i				<-1
				cat('\nprepare plot',i,'/',nrow(tmp2))
				file			<- file.path(tree.indir, infiles[ PTY_RUN==tmp2[i,PTY_RUN],	FILE])
				load(file)
				ph				<- pty.ph[[ tmp2[i,FILE]  ]]
				max.node.height	<- max(node.depth.edgelength(ph)[1:Ntip(ph)])
				tmp				<- c(c(tmp2[i,FILE_ID],'not characterized'), sort(setdiff(sort(levels(attr(ph,'INDIVIDUAL'))),c(tmp2[i,FILE_ID],'not characterized'))))
				col				<- c('black','grey50',rainbow_hcl(length(tmp)-2, start = 270, end = -30, c=100, l=50))
				names(col)		<- tmp			
				ph.title		<- tmp2[i, paste( FILE_ID,'\nrun=',PTY_RUN,', win=',W_FROM,'-',W_TO,'\ndiffind', DIFF_IND,' diff',DIFF, sep='')]
				p				<- ggtree(ph, aes(color=INDIVIDUAL, linetype=TYPE)) + 
						geom_nodepoint(size=ph$node.label/100*3) +
						geom_tiplab(size=1.2,  hjust=-.1) +							 
						scale_color_manual(values=col, guide = FALSE) +											 
						scale_linetype_manual(values=c('target'='solid','filler'='dotted'),guide = FALSE) +
						theme_tree2() +
						theme(legend.position="bottom") + ggplot2::xlim(0, max.node.height*1.3) +
						labs(x='subst/site', title=ph.title)
				p
			})	
	tmp			<- seq_len(ceiling(length(phps)/10))
	pdf(file=gsub('\\.csv','_trees\\.pdf',outfile), w=20, h=40)		#for win=60
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
		pty.infile		<- file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_coinfrunsinput.rda")
		work.dir		<- file.path(HOME,"ptyruns")
		out.dir			<- file.path(HOME,"phylotypes_160119")
		hpc.load		<- ''
	}
	if(0)	#coinfections ZA on HPC
	{			
		pty.infile		<- file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_coinfrunsinput.rda")
		work.dir		<- file.path(HOME,"coinf_ptinput")
		#out.dir			<- file.path(HOME,"coinf_ptoutput_150121")
		out.dir			<- file.path(HOME,"coinf_ptoutput_150201")
		hpc.load		<- "module load intel-suite/2015.1 mpi R/3.3.2"		
	}
	if(1)	#coinfections UG on HPC
	{			
		pty.infile		<- file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_UG_gag_coinfinput_160219.rda")		
		work.dir		<- file.path(HOME,"coinf_ptinput_UG60")
		out.dir			<- file.path(HOME,"coinf_ptoutput_UG60")
		hpc.load		<- "module load intel-suite/2015.1 mpi R/3.3.2"		
	}	
	#	get alignment rda files
	if(0)
	{
		infiles			<- data.table(FILE=list.files(out.dir, pattern='fasta$'))
		infiles			<- subset(infiles, !grepl('*',FILE,fixed=1) & !grepl('dophy\\.fasta',FILE))				
		infiles[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]
		invisible(infiles[, {
							cmd			<- pty.cmd.evaluate.fasta(out.dir, strip.max.len=350, select=paste('ptyr',PTY_RUN,'_',sep=''))
							cat(cmd)
							cmd			<- cmd.hpcwrapper(cmd, hpc.walltime=1, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=1, hpc.load=hpc.load)										
							outfile		<- paste("ptye",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(work.dir, outfile, cmd)
							NULL
						}, by='PTY_RUN'])
		stop()
	}
	#	run ExaML without bootstrap
	if(0)
	{
		pty.args		<- list(	out.dir=out.dir, work.dir=work.dir, 
									outgroup="CPX_AF460972",
									min.ureads.individual=20, min.ureads.candidate=NA, 
									args.examl="-f d -D -m GAMMA", bs.n=0, exa.n.per.run=2)								
		exa.cmd			<- pty.cmdwrap.examl(pty.args)
		#cat( exa.cmd[1, cat(CMD)] )		
		#stop()
		invisible(exa.cmd[,	{		
							cmd			<- cmd.hpcwrapper(CMD, hpc.walltime=40, hpc.q="pqeelab", hpc.mem="5600mb",  hpc.nproc=1, hpc.load=hpc.load)					
							#cmd		<- cmd.hpcwrapper(cmd, hpc.walltime=40, hpc.q="pqeph", hpc.mem="1800mb",  hpc.nproc=1, hpc.load=hpc.load)
							#cat(cmd)
							outfile		<- paste("pexa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
						}, by='RUN_ID'])	
	}	
	#	run ExaML with bootstrap
	if(0)
	{
		pty.args		<- list(	out.dir=out.dir, work.dir=work.dir, 
									outgroup="CPX_AF460972",
									min.ureads.individual=NA, min.ureads.candidate=NA, 
									args.examl="-f d -D -m GAMMA", bs.n=10, exa.n.per.run=1, select=c(7))								
		exa.cmd			<- pty.cmdwrap.examl(pty.args)
		#cat( exa.cmd[1, cat(CMD)] )		
		#stop()
		invisible(exa.cmd[,	{		
							cmd			<- cmd.hpcwrapper(CMD, hpc.walltime=400, hpc.q="pqeelab", hpc.mem="5600mb",  hpc.nproc=1, hpc.load=hpc.load)					
							#cmd		<- cmd.hpcwrapper(cmd, hpc.walltime=10, hpc.q="pqeph", hpc.mem="1800mb",  hpc.nproc=1, hpc.load=hpc.load)
							#cat(cmd)
							outfile		<- paste("pexa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
						}, by='RUN_ID'])	
	}
	#	process newick output into examl.rda files
	if(1)
	{
		pty.args		<- list(	out.dir=out.dir, references.pattern='REF', run.pattern='ptyr',
									outgroup="CPX_AF460972", tree.pattern='_dtbs10.newick$',
									plot.trees.per.page=11, plot.w=150, plot.h=100,
									rm.newick=0, rm.fasta=0)										
		infiles			<- data.table(FILE=list.files(out.dir, pattern= pty.args[['tree.pattern']]))						
		infiles[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]		
		tmp				<- data.table(OUTFILE=list.files(out.dir, pattern='_preprtr.rda$'))						
		tmp[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(OUTFILE,'_'),'[[',1)))]
		infiles			<- merge(infiles, tmp, by='PTY_RUN', all.x=1)		
		setkey(infiles, PTY_RUN)
		infiles			<- subset(unique(infiles), is.na(OUTFILE))
		
		invisible(infiles[, {
							cmd			<- pty.cmd.evaluate.examl(out.dir, 	select=paste('ptyr',PTY_RUN,'_',sep=''), 
																			outgroup=pty.args[['outgroup']],
																			references.pattern=pty.args[['references.pattern']],
																			run.pattern=pty.args[['run.pattern']],
																			tree.pattern=pty.args[['tree.pattern']],
																			rm.newick=pty.args[['rm.newick']],
																			rm.fasta=pty.args[['rm.fasta']],
																			plot.trees.per.page= pty.args[['plot.trees.per.page']],
																			plot.w= pty.args[['plot.w']],
																			plot.h= pty.args[['plot.h']]	)							
							cat(cmd)							
							cmd			<- cmd.hpcwrapper(cmd, hpc.walltime=5, hpc.q="pqeelab", hpc.mem="5600mb",  hpc.nproc=1, hpc.load=hpc.load)										
							outfile		<- paste("ptye",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(work.dir, outfile, cmd)
							#stop()
							NULL
						}, by='PTY_RUN'])
		stop()
	}	
	#	collapse trees
	if(0)
	{
		infiles			<- data.table(FILE=list.files(out.dir, pattern='examl\\.rda$'))						
		infiles[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]
		tmp				<- data.table(OUT=list.files(out.dir, pattern='examlcollapsed\\.rda$'))
		tmp[,FILE:=gsub('collapsed','',OUT)]
		infiles			<- merge(infiles, tmp, by='FILE',all.x=1)
		infiles			<- subset(infiles, is.na(OUT))
		
		invisible(infiles[, {
							cmd			<- pty.cmd.evaluate.examl(pty.infile, out.dir, select=paste('ptyr',PTY_RUN,'_',sep=''))							
							cat(cmd)							
							cmd			<- cmd.hpcwrapper(cmd, hpc.walltime=5, hpc.q="pqeelab", hpc.mem="3600mb",  hpc.nproc=1, hpc.load=hpc.load)										
							outfile		<- paste("ptye",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(work.dir, outfile, cmd)
							#stop()
							NULL
						}, by='PTY_RUN'])
		stop()
	}	
}

pty.pipeline.fasta.160217<- function() 
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
		prog.pty			<- '/Users/Oliver/git/phylotypes/phylotypes.py'
		prog.raxml			<- 'raxmlHPC-AVX'
		no.trees			<- '-T'
		
	}	
	if(0)	#trm pairs on HPC
	{
		load( file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_ptyrunsinput.rda") )
		pty.data.dir	<- '/work/or105/PANGEA_mapout/data'
		work.dir		<- file.path(HOME,"ptyruns")
		out.dir			<- file.path(HOME,"phylotypes")
		prog.pty		<- '/work/or105/libs/phylotypes/phylotypes.py'
		prog.raxml		<- 'raxml'
		no.trees		<- '-T'
		hpc.load		<- "module load intel-suite/2015.1 mpi R/3.3.2 raxml/8.2.4 mafft/7.271 anaconda/2.3.0 samtools"
	}
	if(0)	#coinfections ZA on Mac
	{
		load( file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_coinfrunsinput.rda") )
		pty.data.dir	<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'
		work.dir		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/ptyruns'
		out.dir			<- file.path(HOME,"phylotypes")
		prog.pty		<- '/Users/Oliver/git/phylotypes/phylotypes.py'
		prog.raxml		<- 'raxml'
		no.trees		<- '-T'
		hpc.load		<- ""
	}
	if(0)	#coinfections UG on Mac
	{		
		load( file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_UG_gag_coinfinput_160219.rda") )
		pty.data.dir	<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'
		work.dir		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/coinf_ptinput_UG'
		out.dir			<- file.path(HOME,"coinf_ptoutput_UG60")
		prog.pty		<- '/Users/Oliver/git/phylotypes/phylotypes.py'
		prog.raxml		<- 'raxml'
		no.trees		<- '-T'
		hpc.load		<- ""
	}
	if(1)	#coinfections UG on HPC
	{		
		load( file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_UG_gag_coinfinput_160219.rda") )
		pty.data.dir	<- '/work/or105/PANGEA_mapout/data'
		work.dir		<- file.path(HOME,"coinf_ptinput_UG60")
		out.dir			<- file.path(HOME,"coinf_ptoutput_UG60")
		prog.pty		<- '/work/or105/libs/phylotypes/phylotypes.py'
		prog.raxml		<- 'raxml'
		no.trees		<- '-T'
		hpc.load		<- "module load intel-suite/2015.1 mpi R/3.3.2 raxml/8.2.4 mafft/7.271 anaconda/2.3.0 samtools"
	}
	if(0)	#coinfections ZA on HPC
	{
		load( file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_coinfrunsinput.rda") )
		pty.data.dir	<- '/work/or105/PANGEA_mapout/data'
		work.dir		<- file.path(HOME,"coinf_ptinput")
		out.dir			<- file.path(HOME,"coinf_ptoutput_150201")
		prog.pty		<- '/work/or105/libs/phylotypes/phylotypes.py'
		prog.raxml		<- 'raxml'
		no.trees		<- '-T'
		hpc.load		<- "module load intel-suite/2015.1 mpi R/3.3.2 raxml/8.2.4 mafft/7.271 anaconda/2.3.0 samtools"
	}
	#
	#	set up all temporary files and create bash commands
	#
	#	run 160115	window length 300
	if(0)
	{
		pty.args			<- list(	prog.pty=prog.pty, prog.mafft='mafft', prog.raxml=prog.raxml, data.dir=pty.data.dir, work.dir=work.dir, out.dir=out.dir,
				merge.threshold=1, min.read.count=2, quality.trim.ends=30, min.internal.quality=2, merge.paired.reads='-P',no.trees=no.trees, win=300, keep.overhangs='',
				strip.max.len=350, min.ureads.individual=NA)		
		pty.c				<- pty.cmdwrap.fasta(pty.runs, pty.args)		
	}
	#	run 160118	window length 60 & Q1 18 & keep overhangs
	if(0)
	{		
		pty.args			<- list(	prog.pty=prog.pty, prog.mafft='mafft', prog.raxml=prog.raxml, data.dir=pty.data.dir, work.dir=work.dir, out.dir=out.dir,
				merge.threshold=1, min.read.count=2, quality.trim.ends=18, min.internal.quality=2, merge.paired.reads='-P',no.trees=no.trees, win=60, keep.overhangs='--keep-overhangs',
				strip.max.len=350, min.ureads.individual=20)
		pty.c				<- pty.cmdwrap.fasta(pty.runs, pty.args)	
		pty.c				<- subset(pty.c, PTY_RUN%in%c(3, 9, 12, 15))
	}
	#	run 160119	window length 60 & Q1 18 & keep overhangs & merge.threshold=3
	if(0)
	{		
		pty.args			<- list(	prog.pty=prog.pty, prog.mafft='mafft', prog.raxml=prog.raxml, data.dir=pty.data.dir, work.dir=work.dir, out.dir=out.dir,
				merge.threshold=3, min.read.count=2, quality.trim.ends=18, min.internal.quality=2, merge.paired.reads='-P',no.trees=no.trees, win=60, keep.overhangs='--keep-overhangs',
				strip.max.len=350, min.ureads.individual=NA)
		pty.c				<- pty.cmdwrap.fasta(pty.runs, pty.args)
		pty.c				<- subset(pty.c, PTY_RUN%in%c(24,26,2,31,34,36,38,44,46,60,69,77,70,78,8))
	}
	#	run 160201	window length 60 & Q1 18 & keep overhangs & alignment specified
	if(0)
	{		
		pty.args			<- list(	prog.pty=prog.pty, prog.mafft='mafft', prog.raxml=prog.raxml, data.dir=pty.data.dir, work.dir=work.dir, out.dir=out.dir,
				window.automatic= '', merge.threshold=1, min.read.count=2, quality.trim.ends=18, min.internal.quality=2, merge.paired.reads='-P',no.trees=no.trees, win=60, keep.overhangs='--keep-overhangs',
				strip.max.len=350, min.ureads.individual=20)
		pty.c				<- pty.cmdwrap.fasta(pty.runs, pty.args)
		pty.c[1,cat(CMD)]
		#stop()
		#pty.c				<- subset(pty.c, PTY_RUN%in%c(3, 9, 12, 15))
	}
	#	run 160219	window length 250 & Q1 18 & keep overhangs & alignment specified
	if(1)
	{		
		pty.args			<- list(	prog.pty=prog.pty, prog.mafft='mafft', prog.raxml=prog.raxml, data.dir=pty.data.dir, work.dir=work.dir, out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX.fasta"),
				alignments.root='REF_CPX_AF460972_read_1_count_0', alignments.pairwise.to='REF_B_K03455_read_1_count_0',
				window.automatic= '', merge.threshold=1, min.read.count=2, quality.trim.ends=20, min.internal.quality=2, merge.paired.reads='-P', no.trees=no.trees, 
				strip.max.len=350, min.ureads.individual=NA, win=c(800,9400,60,3), keep.overhangs='--keep-overhangs',
				select=NA)
		pty.c				<- pty.cmdwrap.fasta(pty.runs, pty.args)
		pty.c[1,cat(CMD)]
		#stop()
		#pty.c				<- subset(pty.c, PTY_RUN%in%c(3, 9, 12, 15))
	}
	if(no.trees=='-T')
	{
		invisible(pty.c[,	{					
							cmd			<- cmd.hpcwrapper(CMD, hpc.walltime=30, hpc.q="pqeelab", hpc.mem="11600mb",  hpc.nproc=1, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper(CMD, hpc.walltime=4, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=1, hpc.load=hpc.load)
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

pty.check.assignments.161005<- function()
{
	#	load all read trees
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/160930/RCCS_160930_w270_phscout.rda")
	#	load per window assignments
	load("~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160930_couples_w270_rerun/ptyr1_trmStatsPerWindow.rda")	
	dt	<- subset(window.table, select=c(pat.1, pat.2, W_FROM, W_TO, type))
	dt[, PTY_RUN:=1]
	setnames(dt, c('pat.1','pat.2','type'), c('IND1','IND2','TYPE'))
	#	check on read trees which windows have both individuals
	dt	<- merge(dt, dtrees, by=c('PTY_RUN','W_FROM','W_TO'))
	dt[, FILE_TR:=NULL]	
	tmp	<- dt[, list(HAS_BOTH_IND=	any(grepl(IND1, phs[[IDX]][['tip.label']], fixed=1)) & 
									any(grepl(IND2, phs[[IDX]][['tip.label']], fixed=1))), by=c('PTY_RUN','W_FROM','W_TO','IND1','IND2')]
	dt	<- merge(dt, subset(tmp, HAS_BOTH_IND, select=c('PTY_RUN','W_FROM','W_TO','IND1','IND2')), by=c('PTY_RUN','W_FROM','W_TO','IND1','IND2'))
	set(dt, dt[, which(is.na(TYPE))], 'TYPE', 'disconnected')
	#	now check pair
	dtp	<- subset(dt, IND1=='15861_1_13.bam' & IND2=='15861_1_18.bam')
	write.csv(dtp, row.names=FALSE, file="~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/160930/RCCS_160919_w270-phsc-serodiscpairs-M->F-M-A035225-F-A058111-1845_trmassignments.csv")
}

pty.process.160901<- function()
{
	#
	#	INPUT ARGS PLATFORM
	#	
	if(0)	#UG on Mac
	{				
		in.dir			<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w270'
		save.file.base	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w270/RCCS_160902_w270_phscout.rda'
		in.dir			<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w280'
		save.file.base	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w280/RCCS_160902_w280_phscout.rda'
		in.dir			<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/RCCS_160919_w270_phscout.rda'
		in.dir			<- '~/duke/tmp/Rakai_ptoutput_160930_couples_w270_rerun'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/RCCS_160930_w270_phscout.rda'
		in.dir			<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161007_couples_w270_rerun'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/RCCS_161007_w270_phscout.rda'		
		in.dir			<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d20_rerun'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/RCCS_161027_w270_d20_phscout.rda'
		in.dir			<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d50_rerun'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/RCCS_161027_w270_d50_phscout.rda'
		in.dir			<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d200_rerun'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/RCCS_161027_w270_d200_phscout.rda'
		in.dir			<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d50_r004_rerun'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/RCCS_161027_w270_d50_r004_phscout.rda'
		tmp				<- phsc.combine.phyloscanner.output(in.dir, save.file=save.file)
	}	
	
	#	161107
	if(0)
	{
		in.dir			<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/RCCS_161110_w270'		
		#for(opt in c('d20_r004_rerun'))
		#for(opt in c('d20_r004_rerun','d50_r004_rerun','d200_r004_rerun'))
		for(opt in c('d50_p001_rerun'))		
		{
			#for(trmw.min.reads in c(1, 20, 30, 40, 50, 100))
			for(trmw.min.reads in c(1, 20, 30, 40, 50, 100))
			{
				tmp.min.tips<- ifelse(trmw.min.reads==1, 1, 2)
				tmp.in		<- paste(in.dir, opt, sep='_')
				tmp.out		<- paste(save.file, '_', gsub('_rerun','', opt),'_mr', trmw.min.reads, '_mt', tmp.min.tips, '_phscout.rda', sep='')
				cat('\n',tmp.in,' -> ',tmp.out)
				invisible(phsc.combine.phyloscanner.output(tmp.in, save.file=tmp.out, trmw.min.reads=trmw.min.reads, trmw.min.tips=tmp.min.tips))		
			}
		}
	}
	#	161213
	if(1)
	{
		infile.base		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161213_couples_w270'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/161213/RCCS_161213_w270'
		for(opt in c('d50_p001_rerun'))		
			for(trmw.min.reads in c(20))
				for(trmw.close.brl in c(0.01, 0.02, 0.04, 100, 101))				
				{
					trmw.min.tips		<- ifelse(trmw.min.reads==1, 1, 2)
					trmw.distant.brl	<- ifelse(trmw.close.brl>=101, Inf, 0.07)
					trmw.close.brl		<- ifelse(trmw.close.brl>=100, Inf, trmw.close.brl)					
					tmp.in		<- paste(infile.base, opt, sep='_')
					tmp.out		<- paste(save.file, '_', gsub('_rerun','', opt),'_mr', trmw.min.reads, '_mt', trmw.min.tips, '_cl', 100*trmw.close.brl, '_d', 100*trmw.distant.brl, '_phscout.rda', sep='')
					cat('\n',tmp.in,' -> ',tmp.out)
					invisible(phsc.combine.phyloscanner.output(tmp.in, save.file=tmp.out, trmw.min.reads=trmw.min.reads, trmw.min.tips=trmw.min.tips, trmw.close.brl=trmw.close.brl, trmw.distant.brl=trmw.distant.brl))		
				}		
	}
	#	161219
	if(1)
	{
		infile.base		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161213_couples_w270'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270'
		for(opt in c('d50_p001_rerun'))		
			for(trmw.min.reads in c(20))
				for(trmw.close.brl in c(0.01, 0.02, Inf))
					for(trmw.distant.brl in c(0.05, Inf))
					{
						trmw.min.tips		<- 1
						tmp.in		<- paste(infile.base, opt, sep='_')
						tmp.out		<- paste(save.file, '_', gsub('_rerun','', opt),'_mr', trmw.min.reads, '_mt', trmw.min.tips, '_cl', 100*trmw.close.brl, '_d', 100*trmw.distant.brl, '_phscout.rda', sep='')
						cat('\n',tmp.in,' -> ',tmp.out)
						invisible(phsc.combine.phyloscanner.output(tmp.in, save.file=tmp.out, trmw.min.reads=trmw.min.reads, trmw.min.tips=trmw.min.tips, trmw.close.brl=trmw.close.brl, trmw.distant.brl=trmw.distant.brl))		
					}		
	}
	#	170227
	if(1)
	{
		infile.base		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_170227_couples_w270'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170227/RCCS_170227_w270'		
		for(opt in c('d50_st20_rerun'))		
			for(trmw.min.reads in c(20))
				for(trmw.close.brl in c(0.02))
					for(trmw.distant.brl in c(0.05))
					{
						trmw.min.tips		<- 1
						tmp.in		<- paste(infile.base, opt, sep='_')
						tmp.out		<- paste(save.file, '_', gsub('_rerun','', opt),'_mr', trmw.min.reads, '_mt', trmw.min.tips, '_cl', 100*trmw.close.brl, '_d', 100*trmw.distant.brl, '_phscout.rda', sep='')
						cat('\n',tmp.in,' -> ',tmp.out)
						invisible(phsc.combine.phyloscanner.output(tmp.in, save.file=tmp.out, trmw.min.reads=trmw.min.reads, trmw.min.tips=trmw.min.tips, trmw.close.brl=trmw.close.brl, trmw.distant.brl=trmw.distant.brl))		
					}		
	}
	#	170309
	if(1)
	{
		infile.base		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_170309_couples_w270'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170309/RCCS_170309_w250'		
		for(opt in c('d50_st20_rerun'))		
			for(trmw.min.reads in c(20))
				for(trmw.close.brl in c(0.02))
					for(trmw.distant.brl in c(0.05))
					{
						trmw.min.tips		<- 1
						tmp.in		<- paste(infile.base, opt, sep='_')
						tmp.out		<- paste(save.file, '_', gsub('_rerun','', opt),'_mr', trmw.min.reads, '_mt', trmw.min.tips, '_cl', 100*trmw.close.brl, '_d', 100*trmw.distant.brl, '_phscout.rda', sep='')
						cat('\n',tmp.in,' -> ',tmp.out)
						invisible(phsc.combine.phyloscanner.output(tmp.in, save.file=tmp.out, trmw.min.reads=trmw.min.reads, trmw.min.tips=trmw.min.tips, trmw.close.brl=trmw.close.brl, trmw.distant.brl=trmw.distant.brl))		
					}		
	}
	#	170320
	if(1)
	{
		infile.base		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_170309_couples_w270'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170309/RCCS_170309_w250'		
		for(opt in c('d50_st20_rerun'))		
			for(trmw.min.reads in c(20))
				for(trmw.close.brl in c(0.02))
					for(trmw.distant.brl in c(0.05))
					{
						trmw.min.tips		<- 1
						tmp.in		<- paste(infile.base, opt, sep='_')
						tmp.out		<- paste(gsub('170309','170320',save.file), '_', gsub('_rerun','', opt),'_mr', trmw.min.reads, '_mt', trmw.min.tips, '_cl', 100*trmw.close.brl, '_d', 100*trmw.distant.brl, '_phscout.rda', sep='')
						cat('\n',tmp.in,' -> ',tmp.out)
						invisible(phsc.combine.phyloscanner.output(tmp.in, save.file=tmp.out, trmw.min.reads=trmw.min.reads, trmw.min.tips=trmw.min.tips, trmw.close.brl=trmw.close.brl, trmw.distant.brl=trmw.distant.brl))		
					}		
	}
	#	170324
	if(1)
	{
		infile.base		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_170322_couples_w250'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170322/RCCS_170322_w250'		
		#for(opt in c('d50_st20_trA_rerun','d50_st20_trU_rerun','d50_st20_trB_rerun'))
		for(opt in c('d50_st20_trA_rerun'))
			for(trmw.min.reads in c(20))
				for(trmw.close.brl in c(0.02))
					for(trmw.distant.brl in c(0.05))
					{
						trmw.min.tips		<- 1
						tmp.in		<- paste(infile.base, opt, sep='_')
						tmp.out		<- paste(save.file, '_', gsub('_rerun','', opt),'_mr', trmw.min.reads, '_mt', trmw.min.tips, '_cl', 100*trmw.close.brl, '_d', 100*trmw.distant.brl, '_phscout.rda', sep='')
						cat('\n',tmp.in,' -> ',tmp.out)
						invisible(phsc.combine.phyloscanner.output(tmp.in, save.file=tmp.out, trmw.min.reads=trmw.min.reads, trmw.min.tips=trmw.min.tips, trmw.close.brl=trmw.close.brl, trmw.distant.brl=trmw.distant.brl))		
					}		
	}
	#	170410
	if(1)
	{
		infile.base		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_170405_couples_w250'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170405/RCCS_170405_w250'
		opt 			<- 'd50_st20_trB_rerun'
		trmw.min.reads 	<- 20
		trmw.close.brl 	<- 0.02
		trmw.distant.brl<- 0.05
		trmw.min.tips	<- 1
		in.dir			<- paste(infile.base, opt, sep='_')
		save.file		<- paste(save.file, '_', gsub('_rerun','_blNormedOnFly', opt),'_mr', trmw.min.reads, '_mt', trmw.min.tips, '_cl', 100*trmw.close.brl, '_d', 100*trmw.distant.brl, '_phscout.rda', sep='')		
		invisible(phsc.combine.phyloscanner.output(	in.dir, save.file=save.file, 
													trmw.min.reads=trmw.min.reads, trmw.min.tips=trmw.min.tips, trmw.close.brl=trmw.close.brl, trmw.distant.brl=trmw.distant.brl,
													norm.file.name="/Users/Oliver/Library/R/3.3/library/Phyloscanner.R.utilities/data/hiv.hxb2.norm.constants.rda"))
				
	}
	#	170410
	if(1)
	{
		infile.base		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_170405_couples_w250'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/RCCS_170410_w250'
		opt 			<- 'd50_st20_trB_rerun'
		trmw.min.reads 	<- 20
		trmw.close.brl 	<- 0.035
		trmw.distant.brl<- 0.08
		trmw.min.tips	<- 1
		in.dir			<- paste(infile.base, opt, sep='_')
		save.file		<- paste(save.file, '_', gsub('_rerun','_blNormedOnFly', opt),'_mr', trmw.min.reads, '_mt', trmw.min.tips, '_cl', 100*trmw.close.brl, '_d', 100*trmw.distant.brl, '_phscout.rda', sep='')		
		invisible(phsc.combine.phyloscanner.output(	in.dir, save.file=save.file, 
						trmw.min.reads=trmw.min.reads, trmw.min.tips=trmw.min.tips, trmw.close.brl=trmw.close.brl, trmw.distant.brl=trmw.distant.brl,
						norm.file.name="/Users/Oliver/Library/R/3.3/library/Phyloscanner.R.utilities/data/hiv.hxb2.norm.constants.rda"))
		norm.file.name	<- NA
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/RCCS_170410_w250'
		for(opt in c('d50_st20_trC_rerun','d50_st20_trU_rerun','d50_st20_trB_rerun'))
		{					
				tmp.in		<- paste(infile.base, opt, sep='_')
				tmp.out		<- paste(save.file, '_', gsub('_rerun','', opt),'_mr', trmw.min.reads, '_mt', trmw.min.tips, '_cl', 100*trmw.close.brl, '_d', 100*trmw.distant.brl, '_phscout.rda', sep='')
				cat('\n',tmp.in,' -> ',tmp.out)
				invisible(phsc.combine.phyloscanner.output(tmp.in, save.file=tmp.out, trmw.min.reads=trmw.min.reads, trmw.min.tips=trmw.min.tips, trmw.close.brl=trmw.close.brl, trmw.distant.brl=trmw.distant.brl))		
		}	
		for(opt in c('d50_st20_trB_blInScriptNormed_rerun'))
		{					
			tmp.in		<- paste(infile.base, opt, sep='_')
			tmp.out		<- paste(save.file, '_', gsub('_rerun','', opt),'_mr', trmw.min.reads, '_mt', trmw.min.tips, '_cl', 100*trmw.close.brl, '_d', 100*trmw.distant.brl, '_phscout.rda', sep='')
			cat('\n',tmp.in,' -> ',tmp.out)
			invisible(phsc.combine.phyloscanner.output(tmp.in, save.file=tmp.out, trmw.min.reads=trmw.min.reads, trmw.min.tips=trmw.min.tips, trmw.close.brl=trmw.close.brl, trmw.distant.brl=trmw.distant.brl))		
		}	
	}
	#	170426
	if(1)
	{
		infile.base		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_170426_couples_w250'
		save.file		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170426/RCCS_170426_w250'				
		for(opt in c('d50_p10_blNormed_rerun','d50_p35_blNormed_rerun'))
			for(trmw.min.reads in c(20))
				for(trmw.close.brl in c(0.035))
					for(trmw.distant.brl in c(0.08))
					{
						trmw.min.tips		<- 1
						tmp.in		<- paste(infile.base, opt, sep='_')
						tmp.out		<- paste(save.file, '_', gsub('_rerun','', opt),'_mr', trmw.min.reads, '_mt', trmw.min.tips, '_cl', 100*trmw.close.brl, '_d', 100*trmw.distant.brl, '_phscout.rda', sep='')
						cat('\n',tmp.in,' -> ',tmp.out)
						invisible(phsc.combine.phyloscanner.output(tmp.in, save.file=tmp.out, trmw.min.reads=trmw.min.reads, trmw.min.tips=trmw.min.tips, trmw.close.brl=trmw.close.brl, trmw.distant.brl=trmw.distant.brl))		
					}		
	}
}

pty.pipeline.compress.phyloscanner.output<- function()
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	#
	#	INPUT ARGS PLATFORM
	#	
	if(0)	#dev on Mac
	{
		hpc.load		<- ""
		work.dir		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptinput'
		ptyf			<- data.table(DIR= c('~/duke/2016_PANGEAphylotypes/Rakai_ptoutput_160902_w250'))
		ptyf[, SAVE_FILE_BASE:= file.path(DIR, ptyf[, paste(gsub('_done','',gsub('Rakai_ptoutput','RCCS', basename(DIR))),'_',sep='')])]				
	}
	if(1)	#HPC
	{
		hpc.load		<- "module load R/3.3.2"
		work.dir		<- file.path(HOME,"Rakai_ptinput_160825")
		ptyf			<- data.table(DIR= c(	# '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w250_done',
												# '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w220_done',
												'/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w200_done'))
		ptyf[, SAVE_FILE_BASE:= file.path(DIR, ptyf[, paste(gsub('_done','',gsub('Rakai_ptoutput','RCCS', basename(DIR))),'_',sep='')])]		
	}
	
	cmds	<- ptyf[, {
					cmd	<- phsc.cmd.read.processed.phyloscanner.output.in.directory(DIR, SAVE_FILE_BASE, resume=FALSE, zip=TRUE)
					list(CMD=cmd)
			}, by='DIR']
	#
	#	submit
	#
	invisible(cmds[,	{					
						cmd			<- cmd.hpcwrapper(CMD, hpc.walltime=5, hpc.q="pqeelab", hpc.mem="5900mb",  hpc.nproc=1, hpc.load=hpc.load)
						#cmd		<- cmd.hpcwrapper(CMD, hpc.walltime=4, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=1, hpc.load=hpc.load)
						cat(cmd)					
						outfile		<- paste("pty", paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
						cmd.hpccaller(work.dir, outfile, cmd)
						#stop()
					}, by='DIR'])
	quit('no')	
}

pty.pipeline.dualparameter<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	load( "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	
	#
	#	get dual infection candidates
	#
	indir		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_170208_couples_w270_d50_p50_rerun'
	outfile		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170222_duals/170222-phsc-dualcandidate-'
	infiles		<- data.table(F=list.files(indir, pattern='dualsummary.csv$',full.names=TRUE))
	infiles[, PTY_RUN:= as.integer(gsub('.*ptyr([0-9]+)_.*','\\1',F))]
	setkey(infiles, PTY_RUN)	
	#setdiff(pty.runs[, unique(PTY_RUN)], infiles[, unique(PTY_RUN)])
	#	missing runs: 22  43  52  73  82 115 (phylos never generated in first place)	
	dd			<- infiles[, as.data.table(read.csv(F, stringsAsFactors=FALSE)), by=c('PTY_RUN')]
	setnames(dd, colnames(dd), toupper(colnames(dd)))
	setnames(dd, 'PATIENT', 'SANGER_ID')
	set(dd, NULL, 'SANGER_ID', dd[, gsub('\\.bam','',SANGER_ID)])
	dd.cand		<- unique(subset(dd, PROPORTION>0.2),by='SANGER_ID')	
	#
	#	read all trees and read tree info
	#
	ptyr	<- data.table(FILE=list.files(indir, pattern='trees.rda', full.names=TRUE))	
	phs		<- lapply(seq_len(nrow(ptyr)), function(i)
				{
					load(ptyr[i, FILE])
					phs
				})
	phs		<- do.call(c, phs)	
	dtrees	<- lapply(seq_len(nrow(ptyr)), function(i)
				{
					load(ptyr[i, FILE])
					dfr
				})
	dtrees	<- do.call('rbind', dtrees)
	dtrees[, IDX:=seq_len(nrow(dtrees))]
	#
	#	print putative duals
	#
	invisible(sapply(seq_len(nrow(dd.cand)), function(ii)
					{	
						#ii<- 14
						id			<- dd.cand[ii, SANGER_ID]
						pty.run		<- dd.cand[ii, PTY_RUN]
						dfs			<- subset(dtrees, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
						dfs[, TITLE:= dfs[, paste('potential dual', id, '\ncount', dd.cand[ii, COUNT], ' proportion ', dd.cand[ii, round(PROPORTION, d=2)], '\nrun ', pty.run, '\nwindow ', W_FROM,'-', W_TO,sep='')]]			
						plot.file	<- paste0(outfile,dd.cand[ii, round(PROPORTION*100, d=0)],'-', id,'.pdf')
						invisible(phsc.plot.selected.individuals(phs, dfs, id, plot.cols='red', group.redo=TRUE, plot.file=plot.file, pdf.h=150, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40))					
					}))	
	#
	#	make sliding dual plots for dual infection candidates
	#
	#	read dual info per window
	dds		<- data.table(FD=list.files(indir, pattern='_dualsummary.rda', full.names=TRUE))
	dds[, PTY_RUN:= as.integer(gsub('.*ptyr([0-9]+)_.*','\\1',FD))]
	dds		<- dds[, {
				#FD		<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_170208_couples_w270_d50_p50_rerun/ptyr1_dualsummary.rda'
				load(FD)
				dd[, W_FROM:= as.integer(gsub('InWindow_([0-9]+)_to_([0-9]+)','\\1',W_INFO))]
				dd[, W_TO:= as.integer(gsub('InWindow_([0-9]+)_to_([0-9]+)','\\2',W_INFO))]
				set(dd, NULL, c('W_INFO','W_POTENTIAL_DUAL'), NULL)
				setkey(dd, PATIENT, W_FROM, READS_IN_SUBTREE)
				dd		<- dd[, list(SUBTREE= rev(seq_along(READS_IN_SUBTREE)), READS_IN_SUBTREE=READS_IN_SUBTREE, TIPS_IN_SUBTREE=TIPS_IN_SUBTREE), by=c('PATIENT','W_FROM','W_TO')]				
			}, by='PTY_RUN']
	setnames(dds, 'PATIENT', 'SANGER_ID')
	set(dds, NULL, 'SANGER_ID', dds[,gsub('\\.bam','',SANGER_ID)])
	dds		<- merge(dds, dd.cand, by=c('PTY_RUN','SANGER_ID'))
	#	read #reads and #tips per patient in all windows
	tmp		<- dtrees[, {
				ph	<- phs[[IDX]]
				phb	<- subset(data.table(TAXA=ph$tip.label), !grepl('REF',TAXA))
				phb[, SANGER_ID:= gsub('^([0-9]+_[0-9]+_[0-9]+).*','\\1',TAXA)]
				phb[, READS_IN_TREE:= as.integer(gsub('.*_count_([0-9]+)','\\1',TAXA))]
				phb	<- phb[, list(READS_IN_TREE=sum(READS_IN_TREE), TIPS_IN_TREE=length(READS_IN_TREE)), by='SANGER_ID']
			}, by=c('PTY_RUN','W_FROM','W_TO')]
	tmp		<- merge(tmp, unique(subset(dds, select=c('PTY_RUN','SANGER_ID'))), by=c('PTY_RUN','SANGER_ID')) 	
	dds		<- merge(tmp, dds, by=c('PTY_RUN','SANGER_ID','W_FROM','W_TO'), all.x=1)
	set(dds, dds[, which(is.na(SUBTREE))], 'SUBTREE', 1L)
	tmp		<- dds[, which(is.na(READS_IN_SUBTREE))]
	set(dds, tmp, 'READS_IN_SUBTREE', dds[tmp, READS_IN_TREE])
	tmp		<- dds[, which(is.na(TIPS_IN_SUBTREE))]
	set(dds, tmp, 'TIPS_IN_SUBTREE', dds[tmp, TIPS_IN_TREE])	
	tmp		<- unique(subset(dds, !is.na(COUNT), c(PTY_RUN, SANGER_ID, COUNT, PROPORTION)))
	set(dds, NULL, c('COUNT','PROPORTION'), NULL)
	dds		<- merge(dds, tmp, by=c('PTY_RUN','SANGER_ID'))
	set(dds, NULL, 'SUBTREE', dds[, as.character(factor(SUBTREE))])
	tmp		<- unique(subset(dds, select=c(PTY_RUN, SANGER_ID, COUNT, PROPORTION)))
	tmp		<- tmp[order(-PROPORTION)]
	tmp[, PLOTID:= seq_len(nrow(tmp))]
	tmp[, TITLE:= factor(PLOTID, levels=PLOTID, labels=paste0(SANGER_ID,'\nrun: ',PTY_RUN,' dual count: ',COUNT,' dual prop:  ',round(PROPORTION,d=2),'\n'))]
	dds		<- merge(dds, subset(tmp, select=c(PTY_RUN,SANGER_ID,TITLE)), by=c('PTY_RUN','SANGER_ID'))
	#	plot
	ggplot(dds, aes(x=W_FROM, y=READS_IN_SUBTREE, fill=SUBTREE)) +
			geom_bar(position='stack', stat='identity') +
			scale_x_continuous(breaks=seq(0,1e4,500), expand=c(0,0)) +
			scale_y_log10(breaks=c(10,1e2,1e3,1e4,1e5, 1e6), expand=c(0,0)) +
			geom_hline(yintercept=20) +
			labs(x='window start', y='reads') + facet_wrap(~TITLE, ncol=1) +
			theme_bw()
	ggsave(file=paste0(outfile, 'subtree_scan.pdf'), w=15, h=40)
	#	write csv
	write.csv(subset(tmp, select=c(PTY_RUN, SANGER_ID, COUNT, PROPORTION)), file=paste0(outfile, 'assessments.csv'), row.names=FALSE)
	
}

pty.pipeline.phyloscanner.test<- function() 
{
	require(Phyloscanner.R.utilities)
	#
	#	INPUT ARGS PLATFORM
	#	
	#HOME		<<- "~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA"
	load( file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_UG_gag_coinfinput_160219.rda") )
	setnames(pty.runs, c('FILE_ID'), c('SAMPLE_ID'))
	pty.data.dir		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'
	work.dir			<- '/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptinput'
	out.dir				<- '/Users/Oliver/duke/2016_PANGEAphylotypes/test_ptoutput'	
	prog.pty			<- '/Users/Oliver/git/phylotypes/phyloscanner_make_trees.py'
	prog.raxml			<- '"raxmlHPC-AVX -m GTRCAT -T 1 -p 42"'
	pty.select			<- c(22)
	hpc.load			<- ""
	hpc.nproc			<- 1
	#
	#	INPUT ARGS TEST RUN
	#	
	if(1)
	{
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft='mafft', 
				prog.raxml=prog.raxml, 
				data.dir=pty.data.dir, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file="/Users/Oliver/git/Phyloscanner.R.utilities/inst/HIV1_compendium_AD_B_CPX_v2.fasta",
				alignments.root='REF_CPX_AF460972', 
				bl.normalising.reference.file='/Users/Oliver/git/Phyloscanner.R.utilities/data/hiv.hxb2.norm.constants.rda',
				bl.normalising.reference.var='MEDIAN_PWD',
				alignments.pairwise.to='REF_B_K03455',
				window.automatic= '', 
				merge.threshold=1, 
				min.read.count=1, 
				quality.trim.ends=20, 
				min.internal.quality=2, 
				merge.paired.reads=TRUE, 
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				dont.check.recombination=TRUE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				win=c(2500,3000,250,250), 
				keep.overhangs=FALSE,	
				use.blacklisters=c('ParsimonyBasedBlacklister','DownsampleReads'),				
				roguesubtree.kParam=20,
				roguesubtree.prop.threshold=0,
				roguesubtree.read.threshold=20,
				dwns.maxReadsPerPatient=50,	
				multifurcation.threshold=1e-5,
				split.rule='s',
				split.kParam=20,
				split.proximityThreshold=0.035,
				split.readCountsMatterOnZeroBranches=TRUE,
				pw.trmw.min.reads=20,									
				pw.trmw.min.tips=1,
				pw.trmw.close.brl=0.035,
				pw.trmw.distant.brl=0.08,
				pw.prior.keff=2,
				pw.prior.neff=3,
				pw.prior.keff.dir=2,
				pw.prior.neff.dir=3,				
				pw.prior.calibrated.prob=0.66,
				mem.save=0,
				verbose=TRUE,
				select=pty.select)														
		pty.c				<- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)		
		pty.c[1,cat(CMD)]	
	}
	if(0)
	{
		ptyi	<- seq(2500,2525,25) 
		pty.c	<- lapply(seq_along(ptyi), function(i){
					pty.args			<- list(	prog.pty=prog.pty, 
							prog.mafft='mafft', 
							prog.raxml=prog.raxml, 
							data.dir=pty.data.dir, 
							work.dir=work.dir, 
							out.dir=out.dir, 
							alignments.file="/Users/Oliver/git/Phyloscanner.R.utilities/inst/HIV1_compendium_AD_B_CPX_v2.fasta",
							alignments.root='REF_CPX_AF460972', 
							bl.normalising.reference.file='/Users/Oliver/git/Phyloscanner.R.utilities/data/hiv.hxb2.norm.constants.rda',
							bl.normalising.reference.var='MEDIAN_PWD',
							alignments.pairwise.to='REF_B_K03455',
							window.automatic= '', 
							merge.threshold=1, 
							min.read.count=1, 
							quality.trim.ends=20, 
							min.internal.quality=2, 
							merge.paired.reads=TRUE, 
							no.trees=TRUE, 
							dont.check.duplicates=FALSE,
							num.bootstraps=1,
							all.bootstrap.trees=TRUE,
							strip.max.len=350, 
							min.ureads.individual=NA, 
							#win=c(2500,3000,250,250),
							win=c(ptyi[i],ptyi[i]+250,25,250),
							keep.overhangs=FALSE,					
							verbose=TRUE,
							mem.save=0,
							select=pty.select)														
					pty.c				<- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)
					pty.c[, W_FROM:= ptyi[i]]
					pty.c
				})
		pty.c	<- do.call('rbind', pty.c)		
		pty.c[1,cat(CMD)]	
	}
		
		
}

pty.pipeline.phyloscanner.160825<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	#
	#	INPUT ARGS PLATFORM
	#	
	if(0)	#coinfections UG on Mac
	{		
		load( file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_UG_gag_coinfinput_160219.rda") )
		pty.data.dir		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'
		work.dir			<- '/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptinput'
		out.dir				<- '/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptoutput'
		out.dir				<- '/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptoutput_160902_w250'
		#out.dir				<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/pty_Rakai_160825'
		prog.pty			<- '/Users/Oliver/git/phylotypes/phyloscanner.py'
		prog.raxml			<- '"raxmlHPC-AVX -m GTRCAT -T 1"'
		pty.select			<- c(5,22,99,115)
		hpc.load			<- ""
		hpc.nproc			<- 1
	}
	if(1)	#coinfections UG on HPC
	{		
		load( file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_UG_gag_coinfinput_160219.rda") )
		hpc.load			<- "module load intel-suite/2015.1 mpi R/3.3.2 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"
		hpc.nproc			<- 1
		pty.data.dir		<- '/work/or105/PANGEA_mapout/data'
		work.dir			<- file.path(HOME,"Rakai_ptinput_160825")		
		prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner.py'
		prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-AVX -m GTRCAT"', paste('"raxmlHPC-PTHREADS-AVX -m GTRCAT -T ',hpc.nproc,'"',sep='')) 
		pty.select			<- NA
		pty.select			<- c(22,62,49,85,72)
		#pty.select			<- c(3,84,96)
	}	
	#
	#	INPUT ARGS THIS RUN
	#	
	if(0)
	{
		#	run 160825:	window length 250, min internal 2, trim ends 20, merge.threshold 1
		#				no overhangs, no bootstrap
		out.dir				<- file.path(HOME,"Rakai_ptoutput_160825")
		pty.args			<- list(	prog.pty=prog.pty, 
										prog.mafft='mafft', 
										prog.raxml=prog.raxml, 
										data.dir=pty.data.dir, 
										work.dir=work.dir, 
										out.dir=out.dir, 
										alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX.fasta"),
										alignments.root='REF_CPX_AF460972_read_1_count_0', 
										alignments.pairwise.to='REF_B_K03455_read_1_count_0',
										window.automatic= '', 
										merge.threshold=1, 
										min.read.count=1, 
										quality.trim.ends=20, 
										min.internal.quality=2, 
										merge.paired.reads=TRUE, 
										no.trees=FALSE, 
										dont.check.duplicates=FALSE,
										num.bootstraps=1,
										all.bootstrap.trees=TRUE,
										strip.max.len=350, 
										min.ureads.individual=NA, 
										win=c(800,9400,250,1), 
										keep.overhangs=FALSE, 
										select=pty.select)
	}
	if(0)
	{
		#	run 160902_w280:	window length 280, min internal 23, trim ends 23, merge.threshold 0
		#						no overhangs, bootstrap=1
		out.dir				<- file.path(HOME,"Rakai_ptoutput_160902_w280")		
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft='mafft', 
				prog.raxml=prog.raxml, 
				data.dir=pty.data.dir, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX.fasta"),
				alignments.root='REF_CPX_AF460972_read_1_count_0', 
				alignments.pairwise.to='REF_B_K03455_read_1_count_0',
				window.automatic= '', 
				merge.threshold=0, 
				min.read.count=1, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				merge.paired.reads=TRUE, 
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				win=c(800,9400,280,1), 
				keep.overhangs=FALSE, 
				select=pty.select)
	}
	if(1)
	{
		#	run 160902_w270:	window length 270, min internal 23, trim ends 23, merge.threshold 0
		#						no overhangs, bootstrap=1
		out.dir				<- file.path(HOME,"Rakai_ptoutput_160902_w270")		
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft='mafft', 
				prog.raxml=prog.raxml, 
				data.dir=pty.data.dir, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX.fasta"),
				alignments.root='REF_CPX_AF460972_read_1_count_0', 
				alignments.pairwise.to='REF_B_K03455_read_1_count_0',
				window.automatic= '', 
				merge.threshold=0, 
				min.read.count=1, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				merge.paired.reads=TRUE, 
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				win=c(800,9400,270,1), 
				keep.overhangs=FALSE, 
				select=pty.select)
	}
	if(0)
	{
		#	run 160902_w250:	window length 250, min internal 23, trim ends 23, merge.threshold 0
		#					no overhangs, bootstrap=10
		out.dir				<- file.path(HOME,"Rakai_ptoutput_160902_w250")		
		pty.args			<- list(	prog.pty=prog.pty, 
										prog.mafft='mafft', 
										prog.raxml=prog.raxml, 
										data.dir=pty.data.dir, 
										work.dir=work.dir, 
										out.dir=out.dir, 
										alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX.fasta"),
										alignments.root='REF_CPX_AF460972_read_1_count_0', 
										alignments.pairwise.to='REF_B_K03455_read_1_count_0',
										window.automatic= '', 
										merge.threshold=0, 
										min.read.count=1, 
										quality.trim.ends=23, 
										min.internal.quality=23, 
										merge.paired.reads=TRUE, 
										no.trees=FALSE, 
										dont.check.duplicates=FALSE,
										num.bootstraps=1,#NOTE
										all.bootstrap.trees=TRUE,
										strip.max.len=350, 
										min.ureads.individual=NA, 
										win=c(800,9400,250,1), 
										keep.overhangs=FALSE, 
										select=pty.select)
	}
	if(0)
	{
		#	run 160902_w220:	window length 220, min internal 23, trim ends 23, merge.threshold 0
		#						no overhangs, bootstrap=10
		out.dir				<- file.path(HOME,"Rakai_ptoutput_160902_w220")		
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft='mafft', 
				prog.raxml=prog.raxml, 
				data.dir=pty.data.dir, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX.fasta"),
				alignments.root='REF_CPX_AF460972_read_1_count_0', 
				alignments.pairwise.to='REF_B_K03455_read_1_count_0',
				window.automatic= '', 
				merge.threshold=0, 
				min.read.count=1, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				merge.paired.reads=TRUE, 
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				num.bootstraps=1,#NOTE
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				win=c(800,9400,220,1), 
				keep.overhangs=FALSE, 
				select=pty.select)
	}
	if(0)
	{
		#	run 160902_w200:	window length 200, min internal 23, trim ends 23, merge.threshold 0
		#						no overhangs, bootstrap=10
		out.dir				<- file.path(HOME,"Rakai_ptoutput_160902_w200")		
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft='mafft', 
				prog.raxml=prog.raxml, 
				data.dir=pty.data.dir, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX.fasta"),
				alignments.root='REF_CPX_AF460972_read_1_count_0', 
				alignments.pairwise.to='REF_B_K03455_read_1_count_0',
				window.automatic= '', 
				merge.threshold=0, 
				min.read.count=1, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				merge.paired.reads=TRUE, 
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				num.bootstraps=1,#NOTE
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				win=c(800,9400,200,1), 
				keep.overhangs=FALSE, 
				select=pty.select)
	}
	#
	#	RUN PHYLOSCANNER
	#
	if(1)
	{
		pty.c				<- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)		
		#pty.c[1,cat(CMD)]		
		invisible(pty.c[,	{					
							cmd			<- cmd.hpcwrapper(CMD, hpc.walltime=400, hpc.q="pqeelab", hpc.mem="5900mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper(CMD, hpc.walltime=400, hpc.q="pqeelab", hpc.mem="13900mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd		<- cmd.hpcwrapper(CMD, hpc.walltime=4, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=1, hpc.load=hpc.load)
							cat(cmd)					
							outfile		<- paste("pty",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
							#stop()
						}, by='PTY_RUN'])
		quit('no')
	}	
}

pty.pipeline.phyloscanner.160915.example<- function() 
{
	#set up example
	load( "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	setnames(pty.runs, 'FILE_ID', 'IND_ID')
	runs_rakai_couples	<- subset(pty.runs, select=c(PTY_RUN, IND_ID, TAXA, PARTNER, COUPID))
	save(runs_rakai_couples, file='~/git/Phyloscanner.R.utilities/data/runs_rakai_couples.rda')
	
	# 	load data.table of predefined phyloscanner runs
	require(Phyloscanner.R.utilities)
	data(runs_rakai_couples)
	#	setup input, output and working directories
	pty.data.dir		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'
	work.dir			<- getwd()
	out.dir				<- getwd()
	#	link to phyloscanner.py
	prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner.py'
	#	produce bash scripts for the first two phyloscanner runs only
	pty.select			<- 1:2	
	#	define phyloscanner arguments
	pty.args			<- list(	prog.pty=prog.pty, 
			prog.mafft='mafft', 
			prog.raxml='"raxmlHPC-AVX -m GTRCAT -p 42"', 
			data.dir=pty.data.dir, 
			work.dir=work.dir, 
			out.dir=out.dir, 
			alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX_v2.fasta"),
			alignments.root='REF_CPX_AF460972', 
			alignments.pairwise.to='REF_B_K03455',
			window.automatic= '', 
			merge.threshold=0, 
			min.read.count=1, 
			quality.trim.ends=23, 
			min.internal.quality=23, 
			merge.paired.reads=TRUE, 
			no.trees=FALSE, 
			dont.check.duplicates=FALSE,
			num.bootstraps=1,
			all.bootstrap.trees=TRUE,										 
			min.ureads.individual=NA, 
			win=c(800,9400,25,250), 
			keep.overhangs=FALSE, 
			duplicated.raw.threshold=3,
			duplicated.ratio.threshold=1/200,				
			select=pty.select)
	#	generate bash scripts						
	pty.c				<- phsc.cmd.phyloscanner.multi(runs_rakai_couples, pty.args)
	#	print first bash script to screen
	pty.c[1,cat(CMD)]
}

pty.pipeline.phyloscanner.160915.couples<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	#
	#	INPUT ARGS PLATFORM
	#	
	if(1)
	{		
		#load( "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
		load( file.path(HOME,"data","Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda") )
		setnames(pty.runs, 'FILE_ID', 'IND_ID')
		hpc.load			<- "module load intel-suite/2015.1 mpi R/3.3.2 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"
		hpc.nproc			<- 4
		hpc.mem				<- "5900mb"
		pty.data.dir		<- '/work/or105/PANGEA_mapout/data'
		work.dir			<- file.path(HOME,"Rakai_ptinput_160915_couples")
		out.dir				<- file.path(HOME,"Rakai_ptoutput_160915_couples_w270")
		prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner.py'
		prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-AVX -m GTRCAT -p 42"', paste('"raxmlHPC-PTHREADS-AVX -m GTRCAT -T ',hpc.nproc,' -p 42"',sep='')) 
		#pty.select			<- 11:122
		pty.select			<- c(11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 37, 38, 39, 40, 42, 43, 44, 45, 48, 49, 51, 52, 53, 56, 58, 59, 60, 61, 62, 63, 64, 66, 70, 71, 72, 73, 76, 77, 79, 80, 81, 82, 83, 85, 86, 87, 88, 89, 90, 91, 93, 94, 95, 96, 97, 98, 100, 101, 102, 103, 104, 105, 106, 107, 109, 110, 111, 113, 116, 117, 120, 121)
		#pty.select			<- c(22,62,49,85,72)
		pty.select			<- c(3)
	}	
	#
	#	INPUT ARGS PHYLOSCANNER RUN
	#	
	if(1)
	{				
		pty.args			<- list(	prog.pty=prog.pty, 
										prog.mafft='mafft', 
										prog.raxml=prog.raxml, 
										data.dir=pty.data.dir, 
										work.dir=work.dir, 
										out.dir=out.dir, 
										alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX_v2.fasta"),
										alignments.root='REF_CPX_AF460972', 
										alignments.pairwise.to='REF_B_K03455',
										window.automatic= '', 
										merge.threshold=0, 
										min.read.count=1, 
										quality.trim.ends=23, 
										min.internal.quality=23, 
										merge.paired.reads=TRUE, 
										no.trees=FALSE, 
										dont.check.duplicates=FALSE,
										num.bootstraps=1,
										all.bootstrap.trees=TRUE,
										strip.max.len=350, 
										min.ureads.individual=NA, 
										win=c(800,9400,25,250), 
										keep.overhangs=FALSE, 
										duplicated.raw.threshold=3,
										duplicated.ratio.threshold=1/200,				
										select=pty.select)
	}	
	#
	#	RUN PHYLOSCANNER
	#
	if(1)
	{
		pty.c				<- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)		
		#pty.c[1,cat(CMD)]		
		invisible(pty.c[,	{					
							cmd			<- cmd.hpcwrapper(CMD, hpc.walltime=71, hpc.q="pqeelab", hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper(CMD, hpc.walltime=400, hpc.q="pqeelab", hpc.mem="13900mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd		<- cmd.hpcwrapper(CMD, hpc.walltime=4, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=1, hpc.load=hpc.load)
							cat(cmd)					
							outfile		<- paste("pty",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
							#stop()
						}, by='PTY_RUN'])
		quit('no')
	}	
}

pty.pipeline.phyloscanner.170301.firstbatchofall<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	#
	#	INPUT ARGS PLATFORM
	#	
	if(1)
	{	
		#HOME				<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'
		in.dir				<- file.path(HOME,"RakaiAll_input_170301")
		work.dir			<- file.path(HOME,"RakaiAll_work_170301")
		out.dir				<- file.path(HOME,"RakaiAll_output_170301_w250_s25")
		#load( file.path(in.dir, 'Rakai_phyloscanner_170301.rda') )
		load( file.path(in.dir, 'Rakai_phyloscanner_170301_b75.rda') )
		setnames(pty.runs, c('SID','RENAME_SID','RID'), c('SAMPLE_ID','RENAME_ID','UNIT_ID'))
		hpc.load			<- "module load intel-suite/2015.1 mpi R/3.3.2 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"
		hpc.nproc			<- 4
		hpc.mem				<- "5900mb"						
		prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner_make_trees.py'
		pty.data.dir		<- '/work/or105/PANGEA_mapout/data'
		pty.select			<- 1:2
		prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-AVX -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-AVX -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))
	}	
	if(0)
	{
		prog.pty			<- '/Users/Oliver/git/phylotypes/phyloscanner.py'
		pty.data.dir		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'
		pty.select			<- c(1)
	}			
	#
	#	INPUT ARGS PHYLOSCANNER RUN
	#	
	if(1)
	{				
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft='mafft', 
				prog.raxml=prog.raxml, 
				data.dir=pty.data.dir, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX_v2.fasta"),
				alignments.root='REF_CPX_AF460972', 
				alignments.pairwise.to='REF_B_K03455',
				window.automatic= '', 
				merge.threshold=2, 
				min.read.count=1, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				merge.paired.reads=TRUE, 
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				#win=c(800,9400,125,250), 
				win=c(1000,2000,250,250), #TEST RUN
				keep.overhangs=FALSE,
				use.blacklisters=c('ParsimonyBasedBlacklister','DownsampleReads'),
				sankhoff.k=20,
				roguesubtree.prop.threshold=0,
				roguesubtree.read.threshold=20,
				dwns.maxReadsPerPatient=50,		
				tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$',
				mem.save=0,
				select=pty.select)		
	}	
	#
	#	RUN PHYLOSCANNER
	#
	if(1)
	{
		pty.c				<- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)		
		#pty.c[1,cat(CMD)]		
		invisible(pty.c[,	{
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=71, hpc.q="pqeelab", hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=hpc.load)							
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile		<- paste("scRA",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
							#stop()
						}, by='PTY_RUN'])
		quit('no')
	}	
	#
	#	combine all the data of failed runs
	#	
	if(1)
	{
		indir	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s25_allbatch_sk20_tb_blnormed'
		indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s25_allbatch_sk20_tb_blnormed'
		infiles	<- data.table(TF=list.files(indir, pattern='tree$', full.names=TRUE))
		infiles[, PTY_RUN:=as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(TF)))]
		#	move incomplete files into directory
		#	or delete them!
		invisible(lapply(infiles[, unique(PTY_RUN)], function(pty.run)
			{
				cat(pty.run,'\n')
				tmp	<- data.table(F=list.files(indir, pattern=paste0('^ptyr',pty.run,'_*'), full.names=TRUE))
				file.remove(tmp[,F])
				if(0)
				{
					tmp[, OF:= file.path(dirname(F), file.path(paste0('ptyr',pty.run,'_trees'),basename(F)))]
					dir.create(file.path(indir,paste0('ptyr',pty.run,'_trees')))
					file.rename(tmp[,F], tmp[,OF])	
				}				
			}))
		
		indirs	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s25_allbatch_sk20_tb_blnormed'
		#
		indirs	<- list.files(indirs, pattern='^ptyr[0-9]+_trees$', full.names=TRUE)
		allwin	<- data.table(W_FROM=seq(800,9050,250))		
		#indirs	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo/ptyr97_trees'
		for(i in seq_along(indirs))
		{
			indir	<- indirs[i]
			pty.run	<- as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(indir)))
			#	check if we have all fasta and tree files
			infiles	<- data.table(F=list.files(indir,pattern='ptyr.*fasta$',full.names=TRUE))
			infiles[, W_FROM:= as.integer(gsub('.*_InWindow_([0-9]+)_.*','\\1',basename(F)))]
			infiles	<- merge(allwin, infiles, by='W_FROM', all.x=1)			 					
			tmp		<- subset(infiles, is.na(F))[, W_FROM]
			if(length(tmp))
				cat('\nIn',indir,'Found missing fasta files for',paste(tmp,collapse=', '))
			infiles	<- data.table(F=list.files(indir,pattern='ptyr.*tree$',full.names=TRUE))
			infiles[, W_FROM:= as.integer(gsub('.*_InWindow_([0-9]+)_.*','\\1',basename(F)))]
			infiles	<- merge(allwin, infiles, by='W_FROM', all.x=1)
			misstrs	<- subset(infiles, is.na(F))[, W_FROM]
			if(length(misstrs))
				cat('\nIn',indir,'Found missing tree files for',paste(misstrs,collapse=', '))
			if(!length(tmp) & !length(misstrs))
				cat('\nIn',indir,'Found all fasta and tree files')
			zipit	<- !length(tmp) & !length(misstrs)
			#zipit	<- 1
			if(zipit)
			{			
				cat('\nProcess',indir)
				#	first combine all zip files into ptyrXXX_otherstuff.zip
				infiles	<- data.table(F=list.files(indir,pattern='zip$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				stopifnot(!nrow(subset(infiles, is.na(PTY_RUN))))		
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_otherstuff.zip'))]
				cat('\nZip to file', tmp,'...\n')
				suppressWarnings(invisible( infiles[, list(RTN= unzip(F, overwrite=FALSE, exdir=file.path(indir,'tmp42'))), by='F'] ))
				invisible( infiles[, list(RTN= file.remove(F)), by='F'] )
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	now zip fasta files
				infiles	<- data.table(F=list.files(indir,pattern='ptyr.*fasta$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				invisible( infiles[, file.rename(F, file.path(indir,'tmp42',basename(F))), by='F'] )
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_trees_fasta.zip'))]
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	now zip tree files
				infiles	<- data.table(F=list.files(indir,pattern='ptyr.*tree$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				invisible( infiles[, file.rename(F, file.path(indir,'tmp42',basename(F))), by='F'] )
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_trees_newick.zip'))]		
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	remove tmp dir
				invisible( file.remove( file.path(indir,'tmp42') ) )
				#	move one level down
				infiles	<- data.table(F=list.files(indir, full.names=TRUE))
				invisible( infiles[, file.rename(F, file.path(dirname(indir),basename(F))), by='F'] )
				cat('\nDone',indir)
			}
			if(!length(misstrs))
				invisible(unlink(indir, recursive=TRUE))
			#	expand again if asked to
			if(length(misstrs))
			{
				cat('\nExtract',file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_fasta.zip')))
				unzip(file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_fasta.zip')), junkpaths=TRUE, exdir=indir)
				cat('\nExtract',file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_newick.zip')))
				unzip(file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_newick.zip')), junkpaths=TRUE, exdir=indir)
			}
		}					
	}
}

pty.pipeline.phyloscanner.170301.secondbatchofall<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	#
	#	INPUT ARGS PLATFORM
	#	
	if(1)
	{	
		#HOME				<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'
		in.dir				<- file.path(HOME,"RakaiAll_input_170301")
		#work.dir			<- file.path(HOME,"RakaiAll_work_170301")
		work.dir			<- file.path(HOME,"RakaiAll_work_170704")
		#out.dir				<- file.path(HOME,"RakaiAll_output_170301_w250_s25_secondbatch_sk20_tb_blnormed")
		out.dir				<- file.path(HOME,"RakaiAll_output_170704_w250_s25_secondbatch_sk20_tb_blnormed")		
		#load( file.path(in.dir, 'Rakai_phyloscanner_170301_b75_part2.rda') )
		load( file.path(in.dir, 'Rakai_phyloscanner_170704_b75_part2.rda') )		
		setnames(pty.runs, c('SID','RENAME_SID','RID'), c('SAMPLE_ID','RENAME_ID','UNIT_ID'))
		hpc.load			<- "module load intel-suite/2015.1 mpi R/3.3.2 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"
		#hpc.nproc			<- 4	
		hpc.nproc			<- 1
		prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner_make_trees.py'
		pty.data.dir		<- '/work/or105/PANGEA_mapout/data'
		#prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-SSE3 -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-SSE3 -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))
		prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-AVX -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-AVX -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))
		#pty.select			<- 667:1100
		pty.select			<- 1667:1781
		pty.select			<- 1822:1891
		pty.select			<- 1782:1821
		pty.select			<- c(1811, 1825)
		#pty.select			<- c(1107, 1117, 1123, 1132, 1142, 667, 682, 692, 717)
		#pty.select			<- 910:1100	
	}	
	if(0)
	{
		prog.pty			<- '/Users/Oliver/git/phylotypes/phyloscanner.py'
		pty.data.dir		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'		
	}			
	#
	#	INPUT ARGS PHYLOSCANNER RUN
	#	
	if(1)
	{				
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft='mafft', 
				prog.raxml=prog.raxml, 
				data.dir=pty.data.dir, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX_v2.fasta"),
				alignments.root='REF_CPX_AF460972', 
				alignments.pairwise.to='REF_B_K03455',
				bl.normalising.reference.file=system.file(package="Phyloscanner.R.utilities", "data", "hiv.hxb2.norm.constants.rda"),
				bl.normalising.reference.var='MEDIAN_PWD',
				tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$',				
				window.automatic= '', 
				merge.threshold=2, 
				min.read.count=1, 				
				merge.paired.reads=TRUE, 
				no.trees=TRUE,
				#no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				dont.check.recombination=TRUE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				win=c(800,9400,125,250), 				
				keep.overhangs=FALSE,	
				use.blacklisters=c('ParsimonyBasedBlacklister','DownsampleReads'),				
				roguesubtree.kParam=20,
				roguesubtree.prop.threshold=0,
				roguesubtree.read.threshold=20,
				dwns.maxReadsPerPatient=50,	
				multifurcation.threshold=1e-5,
				split.rule='s',
				split.kParam=20,
				split.proximityThreshold=0.035,	#!!
				split.readCountsMatterOnZeroBranches=TRUE,
				pw.trmw.min.reads=20,									
				pw.trmw.min.tips=1,
				pw.trmw.close.brl=0.035,
				pw.trmw.distant.brl=0.08,
				pw.prior.keff=2,
				pw.prior.neff=3,
				pw.prior.keff.dir=2,
				pw.prior.neff.dir=3,				
				pw.prior.calibrated.prob=0.66,
				mem.save=0,
				verbose=TRUE,
				select=pty.select
				)		 
	}	
	#
	#	RUN PHYLOSCANNER
	#
	if(1)
	{
		pty.c				<- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)		
		#pty.c[1,cat(CMD)]		
		invisible(pty.c[,	{
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=171, hpc.q="pqeelab", hpc.mem="23600mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=171, hpc.q="pqeelab", hpc.mem="5600mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=24, hpc.q=NA, hpc.mem="1890mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=199, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile		<- paste("scRA2",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
							#stop()
						}, by='PTY_RUN'])
		quit('no')
	}	
	if(0)	#check failing runs
	{
		df	<- data.table(F=readLines('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/all_files.txt'))
		df[, PTY_RUN:= as.integer(gsub('ptyr([0-9]+)_.*','\\1',F))]
		tmp	<- df[, list(DONE=any(grepl('pairwise',F))), by='PTY_RUN']
		tmp	<- subset(tmp, !DONE)
		merge(tmp, subset(pty.runs, SID=='15080_1_28', c('PTY_RUN','SID')), by='PTY_RUN', all=TRUE)
		
		subset(tmp, !DONE)[, cat(paste(sort(PTY_RUN), collapse='","'))]
	}
}

pty.pipeline.phyloscanner.170301.firstbatchsecondbatchofall.fix	<- function()
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	
	trmw.min.reads					<- 20
	trmw.min.tips					<- 1
	trmw.close.brl 					<- 0.035
	trmw.distant.brl				<- 0.08		
	prior.keff						<- 2
	prior.neff						<- 3
	prior.calibrated.prob			<- 0.5
	relationship.types				<- c('TYPE_PAIR_DI','TYPE_PAIR_TO','TYPE_PAIR_TODI2x2','TYPE_DIR_TODI3','TYPE_DIRSCORE_TODI3','TYPE_PAIR_TODI2','TYPE_PAIR_TODI','TYPE_PAIRSCORE_TODI','TYPE_CHAIN_TODI')
	
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed'
	#indir	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed'
	infiles	<- data.table(FI=list.files(indir, pattern='_trmStatsPerWindow.rda$', full.names=TRUE))
	infiles[, FO:=gsub('_trmStatsPerWindow.rda','_pairwise_relationships.rda',FI)]
	infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FI)))]
	setkey(infiles, PTY_RUN)
	work.dir<- file.path(HOME,"RakaiAll_work_170301")
	
	infiles	<- infiles[, list(CMD= phsc.cmd.pairwise.relationships(	FI, 
																	FO, 
																	trmw.min.reads=trmw.min.reads, 
																	trmw.min.tips=trmw.min.tips, 
																	trmw.close.brl=trmw.close.brl, 
																	trmw.distant.brl=trmw.distant.brl, 
																	prior.keff=prior.keff, 
																	prior.neff=prior.neff, 
																	prior.calibrated.prob=prior.calibrated.prob, 
																	rel.likely.pair=TRUE, rel.likely.pair.by.distance.only=TRUE,rel.likely.pair.by.topology.only=TRUE,rel.likely.pair.by.cross.table=TRUE,rel.direction=TRUE,rel.chain=TRUE)), by='PTY_RUN']
	infiles[, BATCH:=ceiling(PTY_RUN/20)]										
	invisible(infiles[,	{
						#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=171, hpc.q="pqeelab", hpc.mem="5900mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
						cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=0, hpc.q=NA, hpc.mem="1899mb",  hpc.nproc=1, hpc.load="module load R/3.3.2")							
						cmd			<- paste(cmd, paste(CMD, collapse='\n'),sep='\n')
						cat(cmd)	
						stop()
						outfile		<- paste("scRA5",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
						cmd.hpccaller(work.dir, outfile, cmd)
						#stop()
					}, by='BATCH'])
	quit('no')	
}

pty.pipeline.phyloscanner.170301.secondstage.ptyrtrees<- function() 
{
	require(data.table)
	require(big.phylo)
	#
	#	produce trees
	#
	if(0)
	{
		#HOME				<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'	
		hpc.load			<- "module load intel-suite/2015.1 mpi raxml/8.2.9"
		if(0)	#first lightweight run to handle most read alignments
		{
			#hpc.select<- 1; hpc.nproc<- 1; 	hpc.walltime<- 3; hpc.mem<- "1850mb"; hpc.q<- NA
			hpc.select<- 1; hpc.nproc<- 1; 	hpc.walltime<- 3; hpc.mem<- "5900mb"; hpc.q<- "pqeelab"
		}
		if(1)	#second midweight run to handle the remaining read alignments
		{
			#hpc.select<- 1; hpc.nproc<- 1; 	hpc.walltime<- 998; hpc.mem<- "3600mb"; hpc.q<- "pqeph"
			#hpc.select<- 1; hpc.nproc<- 1; 	hpc.walltime<- 998; hpc.mem<- "5900mb"; hpc.q<- "pqeelab"			
			hpc.select<- 1; hpc.nproc<- 1; 	hpc.walltime<- 71; hpc.mem<- "1800mb"; hpc.q<- NA
		}
		if(0)	#third heavyweight run to handle the remaining read alignments
		{
			#hpc.select<- 1; hpc.nproc<- 8; 	hpc.walltime<- 71; hpc.mem<- "7850mb"; hpc.q<- NA
			hpc.select<- 1; hpc.nproc<- 1; 	hpc.walltime<- 71; hpc.mem<- "63850mb"; hpc.q<- NA
		}
		
		raxml.pr			<- ifelse(hpc.nproc==1, 'raxmlHPC-SSE3', 'raxmlHPC-PTHREADS-SSE3')	
		#raxml.pr			<- ifelse(hpc.nproc==1, 'raxmlHPC-AVX','raxmlHPC-PTHREADS-AVX')
		raxml.args			<- ifelse(hpc.nproc==1, '-m GTRCAT --HKY85 -p 42 -o REF_CPX_AF460972', paste0('-m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42 -o REF_CPX_AF460972'))
		#in.dir				<- file.path(HOME,'RakaiAll_output_170301_w250_s20_p35_stagetwo')
		in.dir				<- file.path(HOME,'RakaiAll_output_170704_w250_s20_p35_stagetwo')		
		out.dir				<- in.dir
		work.dir			<- file.path(HOME,"RakaiAll_work_170704")
		
		infiles	<- data.table(FI=list.files(in.dir, pattern='fasta$', full.names=TRUE, recursive=TRUE))
		infiles[, FO:= gsub('fasta$','tree',FI)]
		infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FI)))]
		infiles[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FI)))]		
		tmp		<- data.table(FT=list.files(out.dir, pattern='^ptyr.*tree$', full.names=TRUE, recursive=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FT)))]
		tmp[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FT)))]
		infiles	<- merge(infiles, tmp, by=c('PTY_RUN','W_FROM'), all.x=1)
		infiles	<- subset(infiles, is.na(FT))	
		setkey(infiles, PTY_RUN, W_FROM)	
		 
		#infiles	<- subset(infiles, PTY_RUN>=240 & PTY_RUN<300)
		infiles	<- subset(infiles, PTY_RUN>=186 & PTY_RUN<239)
		print(infiles)	 		
		
		df		<- infiles[, list(CMD=cmd.raxml(FI, outfile=FO, pr=raxml.pr, pr.args=raxml.args)), by=c('PTY_RUN','W_FROM')]
		df[, ID:=ceiling(seq_len(nrow(df))/5)]
		df		<- df[, list(CMD=paste(CMD, collapse='\n',sep='')), by='ID']
		
		#df[1, cat(CMD)]
		#stop()
		invisible(df[,	{
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)							
							outfile		<- paste("srx1",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(work.dir, outfile, cmd)
							#stop()
						}, by=c('ID')])
	}
	#
	#	combine all the data
	#	
	if(1)
	{
		indirs 	<- '/Users/Oliver/duke/tmp/ptyr143_trees'
		indirs	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo'
		indirs	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s25_secondbatch_sk20_tb_blnormed'
		indirs	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo'
		#
		indirs	<- list.files(indirs, pattern='^ptyr[0-9]+_trees$', full.names=TRUE)
		allwin	<- data.table(W_FROM=seq(800,9150,25))
		#allwin	<- data.table(W_FROM=seq(800,9050,125))
		#indirs	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo/ptyr97_trees'
		for(i in seq_along(indirs))
		{
			indir	<- indirs[i]
			pty.run	<- as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(indir)))
			#	check if we have all fasta and tree files
			infiles	<- data.table(F=list.files(indir,pattern='ptyr.*fasta$',full.names=TRUE))
			infiles[, W_FROM:= as.integer(gsub('.*_InWindow_([0-9]+)_.*','\\1',basename(F)))]
			infiles	<- merge(allwin, infiles, by='W_FROM', all.x=1)			 					
			missfs	<- subset(infiles, is.na(F))[, W_FROM]
			if(length(missfs))
				cat('\nIn',indir,'Found missing fasta files for',paste(missfs,collapse=', '))
			infiles	<- data.table(F=list.files(indir,pattern='ptyr.*tree$',full.names=TRUE))
			infiles[, W_FROM:= as.integer(gsub('.*_InWindow_([0-9]+)_.*','\\1',basename(F)))]
			infiles	<- merge(allwin, infiles, by='W_FROM', all.x=1)
			misstrs	<- subset(infiles, is.na(F))[, W_FROM]
			if(length(misstrs))
				cat('\nIn',indir,'Found missing tree files for',paste(misstrs,collapse=', '))
			zipit	<- 0
			if(!length(missfs) & !length(misstrs))
			{
				cat('\nIn',indir,'Found all fasta and tree files')
				zipit	<- 1
			}				
			if(!length(setdiff(misstrs,missfs)))
			{ 
				cat('\nIn',indir,'Found all tree files for which there is a fasta file')
				zipit	<- 1
			}	
			#
			if(zipit)
			{			
				cat('\nProcess',indir)
				#	first combine all zip files into ptyrXXX_otherstuff.zip
				infiles	<- data.table(F=list.files(indir,pattern='zip$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				stopifnot(!nrow(subset(infiles, is.na(PTY_RUN))))		
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_otherstuff.zip'))]
				cat('\nZip to file', tmp,'...\n')
				suppressWarnings(invisible( infiles[, list(RTN= unzip(F, overwrite=FALSE, exdir=file.path(indir,'tmp42'))), by='F'] ))
				invisible( infiles[, list(RTN= file.remove(F)), by='F'] )
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	now zip fasta files
				infiles	<- data.table(F=list.files(indir,pattern='ptyr.*fasta$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				invisible( infiles[, file.rename(F, file.path(indir,'tmp42',basename(F))), by='F'] )
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_trees_fasta.zip'))]
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	now zip tree files
				infiles	<- data.table(F=list.files(indir,pattern='ptyr.*tree$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				invisible( infiles[, file.rename(F, file.path(indir,'tmp42',basename(F))), by='F'] )
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_trees_newick.zip'))]		
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	remove tmp dir
				invisible( file.remove( file.path(indir,'tmp42') ) )
				#	move one level down
				infiles	<- data.table(F=list.files(indir, full.names=TRUE))
				invisible( infiles[, file.rename(F, file.path(dirname(indir),basename(F))), by='F'] )
				cat('\nDone',indir)
			}
			#if(!length(misstrs))
			if(zipit)
				invisible(unlink(indir, recursive=TRUE))
			#	expand again if asked to
			#if(length(misstrs))
			if(0)
			{
				cat('\nExtract',file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_fasta.zip')))
				unzip(file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_fasta.zip')), junkpaths=TRUE, exdir=indir)
				cat('\nExtract',file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_newick.zip')))
				unzip(file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_newick.zip')), junkpaths=TRUE, exdir=indir)
			}
		}					
	}
	quit('no')	
}


pty.pipeline.phyloscanner.180302.beehive67.process<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	
	# set up pty.runs file
	# with cols	'SAMPLE_ID','RENAME_ID','UNIT_ID'
	if(0)	
	{
		infile	<- '~/Dropbox (SPH Imperial College)/2017_phyloscanner_beehive_cluster/Cluster_67_wClose_in.csv'
		outfile	<- '~/Dropbox (SPH Imperial College)/2017_phyloscanner_beehive_cluster/Cluster_67_wClose_in.rda'
		df		<- as.data.table(read.csv(infile, header=FALSE, stringsAsFactors=FALSE))
		setnames(df, c('V1','V2','V3'), c('BAM','REF','RENAME_ID'))
		df[, SAMPLE_ID:= gsub('_remap_dedup.bam','',basename(BAM))]
		df[, UNIT_ID:= gsub('^([^-]+)-.*','\\1',RENAME_ID)]
		df		<- df[, list(SAMPLE_ID=SAMPLE_ID, RENAME_ID=paste0(UNIT_ID,'-',1:length(SAMPLE_ID))), by= 'UNIT_ID']		
		df[, PTY_RUN:=67]
		pty.runs<- subset(df, select=c(PTY_RUN, SAMPLE_ID, RENAME_ID, UNIT_ID  ))		
		save(pty.runs, file=outfile)
	}
	#
	#	INPUT ARGS
	if(1)
	{	
		#HOME				<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'
		in.dir				<- file.path(HOME,"BEEHIVE_67_180302_in")
		work.dir			<- file.path(HOME,"BEEHIVE_67_180302_work")
		out.dir				<- file.path(HOME,"BEEHIVE_67_180302_out")		
		load( file.path(in.dir, 'Cluster_67_wClose_in.rda') )
		print(pty.runs)
		hpc.load			<- "module load intel-suite/2015.1 mpi R/3.3.3 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"
		hpc.nproc			<- 1							
		prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner_make_trees.py'
		pty.data.dir		<- '/work/cw109/BEEHIVE/BamsDeduped'
		#prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-SSE3 -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-SSE3 -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))
		prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-AVX -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-AVX -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))
		pty.select			<- 67		
	}	
	
	#
	# generate read alignments
	if(0)
	{		
		ptyi		<- seq(800,9150,25) 
		pty.c		<- lapply(seq_along(ptyi), function(i)
				{
					pty.args			<- list(	prog.pty=prog.pty, 
							prog.mafft='mafft', 
							prog.raxml=prog.raxml, 
							data.dir=pty.data.dir, 
							work.dir=work.dir, 
							out.dir=out.dir, 
							alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX_v2.fasta"),
							alignments.root='REF_CPX_AF460972', 
							alignments.pairwise.to='REF_B_K03455',
							window.automatic= '', 
							merge.threshold=0, 
							min.read.count=1, 
							quality.trim.ends=23, 
							min.internal.quality=23, 
							merge.paired.reads=TRUE, 
							no.trees=TRUE, 
							dont.check.duplicates=FALSE,
							dont.check.recombination=TRUE,
							num.bootstraps=1,
							all.bootstrap.trees=TRUE,
							strip.max.len=350, 
							min.ureads.individual=NA, 
							win=c(ptyi[i],ptyi[i]+250,25,250),				 				
							keep.overhangs=FALSE,
							mem.save=0,
							verbose=TRUE,					
							select=pty.select	#of 240
					)											
					pty.c				<- phsc.cmd.phyloscanner.multi(pty.runs, pty.args,regex.ref='_ref.fasta$', postfix.sample.id='_remap_dedup\\.bam|_remap_ref\\.fasta')
					#cat(pty.c$CMD)
					#stop()
					pty.c[, W_FROM:= ptyi[i]]
					pty.c
				})
		pty.c	<- do.call('rbind', pty.c)	
		tmp		<- data.table(FO=list.files(out.dir, pattern='ptyr.*fasta$', recursive=TRUE, full.names=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FO)))]
		tmp[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FO)))]
		pty.c	<- merge(pty.c, tmp, by=c('PTY_RUN','W_FROM'), all.x=1)
		pty.c	<- subset(pty.c, is.na(FO))		
		#print(pty.c)
		#stop()		 
		invisible(pty.c[,	{
							cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=18, hpc.q="pqeelab", hpc.mem="6gb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd	<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=71, hpc.q=NA, hpc.mem="1850mb",  hpc.nproc=1, hpc.load=hpc.load)
							#cmd	<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=71, hpc.q=NA, hpc.mem="63800mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							cmd		<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile	<- paste("bee",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(work.dir, outfile, cmd)
							#stop()
						}, by=c('PTY_RUN','W_FROM')])
		quit('no')
	}
	#
	# generate trees
	if(0)
	{
		#HOME		<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'	
		hpc.load	<- "module load intel-suite/2015.1 mpi raxml/8.2.9"
		#hpc.select	<- 1; hpc.nproc<- 1; 	hpc.walltime<- 167; hpc.mem<-"6gb"; hpc.q<- "pqeelab"
		hpc.select	<- 1; hpc.nproc<- 16; 	hpc.walltime<- 23; hpc.mem<-"124gb"; hpc.q<- NA
		#raxml.pr	<- ifelse(hpc.nproc==1, 'raxmlHPC-SSE3', 'raxmlHPC-PTHREADS-SSE3')	
		raxml.pr	<- ifelse(hpc.nproc==1, 'raxmlHPC-AVX','raxmlHPC-PTHREADS-AVX')
		raxml.args	<- ifelse(hpc.nproc==1, '-m GTRCAT --HKY85 -p 42 -o REF_CPX_AF460972', paste0('-m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42 -o REF_CPX_AF460972'))
		in.dir		<- file.path(HOME,"BEEHIVE_67_180302_out")
		out.dir		<- in.dir
		work.dir	<- file.path(HOME,"BEEHIVE_67_180302_work")
		
		infiles	<- data.table(FI=list.files(in.dir, pattern='fasta$', full.names=TRUE, recursive=TRUE))
		infiles[, FO:= gsub('fasta$','tree',FI)]
		infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FI)))]
		infiles[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FI)))]		
		tmp		<- data.table(FT=list.files(out.dir, pattern='^ptyr.*tree$', full.names=TRUE, recursive=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FT)))]
		tmp[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FT)))]
		infiles	<- merge(infiles, tmp, by=c('PTY_RUN','W_FROM'), all.x=1)
		infiles	<- subset(infiles, is.na(FT))	
		setkey(infiles, PTY_RUN, W_FROM)			
		#infiles	<- subset(infiles, PTY_RUN>=240 & PTY_RUN<300)
		#infiles	<- subset(infiles, PTY_RUN>=186 & PTY_RUN<239)
		#print(infiles)	 		
		
		df		<- infiles[, list(CMD=cmd.raxml(FI, outfile=FO, pr=raxml.pr, pr.args=raxml.args)), by=c('PTY_RUN','W_FROM')]
		df[, ID:=ceiling(seq_len(nrow(df))/1)]
		#df		<- df[, list(CMD=paste(CMD, collapse='\n',sep='')), by='ID']
		
		#df[1, cat(CMD)]
		#stop()
		invisible(df[,	{
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)							
							outfile		<- paste("srx1",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(work.dir, outfile, cmd)
							#stop()
						}, by=c('ID')])
	}
	#
	#	combine all the data	
	if(0)
	{
		indirs 	<- '/Users/Oliver/duke/tmp/ptyr143_trees'
		indirs	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/BEEHIVE_67_180302_out'
		#
		indirs	<- list.files(indirs, pattern='^ptyr[0-9]+_trees$', full.names=TRUE)
		allwin	<- data.table(W_FROM=seq(800,9150,25))
		#allwin	<- data.table(W_FROM=seq(800,9050,125))
		#indirs	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo/ptyr97_trees'
		for(i in seq_along(indirs))
		{
			indir	<- indirs[i]
			pty.run	<- as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(indir)))
			#	check if we have all fasta and tree files
			infiles	<- data.table(F=list.files(indir,pattern='ptyr.*fasta$',full.names=TRUE))
			infiles[, W_FROM:= as.integer(gsub('.*_InWindow_([0-9]+)_.*','\\1',basename(F)))]
			infiles	<- merge(allwin, infiles, by='W_FROM', all.x=1)			 					
			missfs	<- subset(infiles, is.na(F))[, W_FROM]
			if(length(missfs))
				cat('\nIn',indir,'Found missing fasta files for',paste(missfs,collapse=', '))
			infiles	<- data.table(F=list.files(indir,pattern='ptyr.*tree$',full.names=TRUE))
			infiles[, W_FROM:= as.integer(gsub('.*_InWindow_([0-9]+)_.*','\\1',basename(F)))]
			infiles	<- merge(allwin, infiles, by='W_FROM', all.x=1)
			misstrs	<- subset(infiles, is.na(F))[, W_FROM]
			if(length(misstrs))
				cat('\nIn',indir,'Found missing tree files for',paste(misstrs,collapse=', '))
			zipit	<- 0
			if(!length(missfs) & !length(misstrs))
			{
				cat('\nIn',indir,'Found all fasta and tree files')
				zipit	<- 1
			}				
			if(!length(setdiff(misstrs,missfs)))
			{ 
				cat('\nIn',indir,'Found all tree files for which there is a fasta file')
				zipit	<- 1
			}	
			#
			if(zipit)
			{			
				cat('\nProcess',indir)
				#	first combine all zip files into ptyrXXX_otherstuff.zip
				infiles	<- data.table(F=list.files(indir,pattern='zip$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				stopifnot(!nrow(subset(infiles, is.na(PTY_RUN))))		
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_otherstuff.zip'))]
				cat('\nZip to file', tmp,'...\n')
				suppressWarnings(invisible( infiles[, list(RTN= unzip(F, overwrite=FALSE, exdir=file.path(indir,'tmp42'))), by='F'] ))
				invisible( infiles[, list(RTN= file.remove(F)), by='F'] )
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	now zip fasta files
				infiles	<- data.table(F=list.files(indir,pattern='ptyr.*fasta$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				invisible( infiles[, file.rename(F, file.path(indir,'tmp42',basename(F))), by='F'] )
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_trees_fasta.zip'))]
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	now zip tree files
				infiles	<- data.table(F=list.files(indir,pattern='ptyr.*tree$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				invisible( infiles[, file.rename(F, file.path(indir,'tmp42',basename(F))), by='F'] )
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_trees_newick.zip'))]		
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	remove tmp dir
				invisible( file.remove( file.path(indir,'tmp42') ) )
				#	move one level down
				infiles	<- data.table(F=list.files(indir, full.names=TRUE))
				invisible( infiles[, file.rename(F, file.path(dirname(indir),basename(F))), by='F'] )
				cat('\nDone',indir)
			}
			#if(!length(misstrs))
			if(zipit)
				invisible(unlink(indir, recursive=TRUE))
			#	expand again if asked to
			#if(length(misstrs))
			if(0)
			{
				cat('\nExtract',file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_fasta.zip')))
				unzip(file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_fasta.zip')), junkpaths=TRUE, exdir=indir)
				cat('\nExtract',file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_newick.zip')))
				unzip(file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_newick.zip')), junkpaths=TRUE, exdir=indir)
			}
		}					
	}
	#
	#	Input args 2
	if(1)	
	{	
		in.dir				<- file.path(HOME,"BEEHIVE_67_180302_out")
		out.dir				<- in.dir
		work.dir			<- file.path(HOME,"BEEHIVE_67_180302_work")						
		dir.create(out.dir, showWarnings=FALSE)
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft=NA, 
				prog.raxml=NA, 
				data.dir=NA, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX_v2.fasta"),
				alignments.root='REF_CPX_AF460972', 
				alignments.pairwise.to='REF_B_K03455',
				bl.normalising.reference.file=system.file(package="Phyloscanner.R.utilities", "data", "hiv.hxb2.norm.constants.rda"),
				bl.normalising.reference.var='MEDIAN_PWD',														
				window.automatic= '', 
				merge.threshold=0, 
				min.read.count=1, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				merge.paired.reads=TRUE, 
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				dont.check.recombination=TRUE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				win=c(800,9400,25,250), 				
				keep.overhangs=FALSE,
				use.blacklisters=c('ParsimonyBasedBlacklister','DownsampleReads'),
				tip.regex='^(.*)-[0-9]+_read_([0-9]+)_count_([0-9]+)$',
				roguesubtree.kParam=20,
				roguesubtree.prop.threshold=0,
				roguesubtree.read.threshold=20,
				dwns.maxReadsPerPatient=50,	
				multifurcation.threshold=1e-5,
				split.rule='s',
				split.kParam=20,
				split.proximityThreshold=0,
				split.readCountsMatterOnZeroBranches=TRUE,
				split.pruneBlacklist=FALSE,
				trms.allowMultiTrans=TRUE,
				pw.trmw.min.reads=30,									
				pw.trmw.min.tips=1,
				pw.trmw.close.brl=0.025,
				pw.trmw.distant.brl=0.05,
				pw.prior.keff=2,
				pw.prior.neff=3,
				pw.prior.keff.dir=2,
				pw.prior.neff.dir=3,				
				pw.prior.calibrated.prob=0.6,
				mem.save=0,
				verbose=TRUE,				
				select=NA 
		)	
		save(pty.args, file=file.path(out.dir, 'pty.args.rda'))
	
		pty.c	<- data.table(FILE_BAM=list.files(in.dir, pattern='_patients.txt', full.names=TRUE))
		pty.c[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_patients.txt','',basename(FILE_BAM))))]				
		tmp		<- data.table(FILE_TRMW=list.files(out.dir, pattern='_pairwise_relationships.rda', full.names=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_pairwise_relationships.rda','',basename(FILE_TRMW))))]
		pty.c	<- merge(pty.c, tmp, by='PTY_RUN', all.x=1)
		pty.c	<- subset(pty.c, is.na(FILE_TRMW))
		setkey(pty.c, PTY_RUN)		
		pty.c	<- pty.c[, { 
					#FILE_BAM<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr1_bam.txt'
					#FILE_BAM<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr1_bam.txt'
					#cat('\n',FILE_BAM)
					prefix.infiles	<- gsub('patients.txt','',FILE_BAM)
					print(prefix.infiles)
					cmd				<- phsc.cmd.phyloscanner.one.resume(prefix.infiles, pty.args)
					list(CMD=cmd) 
				}, by='PTY_RUN']		
		pty.c[1,cat(CMD)]		
		#stop()		
		invisible(pty.c[,	{					
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=15, hpc.q="pqeelab", hpc.mem="6gb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=168, hpc.q="pqeelab", hpc.mem="48gb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=15, hpc.q="pqeph", hpc.mem="3.6gb",  hpc.nproc=1, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=23, hpc.q=NA, hpc.mem="96gb",  hpc.nproc=1, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=23, hpc.q=NA, hpc.mem="190gb",  hpc.nproc=1, hpc.load=hpc.load)
							cmd			<- paste(cmd,'cd $TMPDIR',sep='\n')
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile		<- paste("scRAr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
						}, by='PTY_RUN'])		
		quit('no')		
	}
}


pty.pipeline.phyloscanner.180605.MunichCluster.readlengths<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	infile.base	<- '~/Box Sync/OR_Work/2018/2018_MunichCluster/mergedfragmentlen__'
	infile.bams	<- paste0(infile.base,c('cov175.rda','cov200.rda','cov225.rda','cov250.rda','cov275.rda','cov300.rda','cov325.rda','cov350.rda','cov375.rda','lendist.rda'))		
	outfile.base<- '~/Box Sync/OR_Work/2018/2018_MunichCluster/MunichCluster_bamstats_'						
		
	for(i in seq_along(infile.bams))		
		load(infile.bams[i])	
	setnames(bam.cov175, 'FILE_ID', 'SID')
	bam.cov175	<- subset(bam.cov175, !is.na(COV))
	setnames(bam.cov200, 'FILE_ID', 'SID')
	bam.cov200	<- subset(bam.cov200, !is.na(COV))
	setnames(bam.cov225, 'FILE_ID', 'SID')
	bam.cov225	<- subset(bam.cov225, !is.na(COV))	
	setnames(bam.cov250, 'FILE_ID', 'SID')
	bam.cov250	<- subset(bam.cov250, !is.na(COV))
	setnames(bam.cov275, 'FILE_ID', 'SID')
	bam.cov275	<- subset(bam.cov275, !is.na(COV))
	setnames(bam.cov300, 'FILE_ID', 'SID')
	bam.cov300	<- subset(bam.cov300, !is.na(COV))
	setnames(bam.cov325, 'FILE_ID', 'SID')
	bam.cov325	<- subset(bam.cov325, !is.na(COV))	
	setnames(bam.cov350, 'FILE_ID', 'SID')
	bam.cov350	<- subset(bam.cov350, !is.na(COV))
	setnames(bam.cov375, 'FILE_ID', 'SID')
	bam.cov375	<- subset(bam.cov375, !is.na(COV))	
	setnames(bam.len, 'FILE_ID', 'SID')
	bam.len[, LONGER:= factor(CDF>1-1e-4, levels=c(TRUE,FALSE), labels=c('no','yes'))]
	
	
	ggplot(subset(bam.len,QU>=40 & QU<320), aes(y=1-CDF, x=factor(QU))) + 
			geom_boxplot() +
			scale_y_continuous(lab=scales:::percent, expand=c(0,0)) +
			theme_bw() + labs(y='proportion of reads longer than x in one individual\n(boxplot over all individuals)', x='\nx=length of quality-trimmed short reads\n(nt)', fill='sequence run') +
			theme(legend.position='bottom')
	ggsave(file=paste0(outfile.base,'bamlen_longer_than_xnt.pdf'), w=6, h=6)
	#	comparing to beehive, perhaps a good threshold is 250	
	
		
	#
	df			<- copy(bam.cov175)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=175L]
	ans			<- copy(tmp)		
	#
	df			<- copy(bam.cov200)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=200L]
	ans			<- rbind(ans, tmp)	
	#
	df			<- copy(bam.cov225)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=225L]
	ans			<- rbind(ans, tmp)	
	#
	df			<- copy(bam.cov250)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=250L]
	ans			<- rbind(ans, tmp)		
	#
	df			<- copy(bam.cov275)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=275L]
	ans			<- rbind(ans, tmp)	
	#
	df			<- copy(bam.cov300)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=300L]
	ans			<- rbind(ans, tmp)	
	#
	df			<- copy(bam.cov325)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=325L]
	ans			<- rbind(ans, tmp)	
	#
	df			<- copy(bam.cov350)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=350L]
	ans			<- rbind(ans, tmp)	
	#
	df			<- copy(bam.cov375)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{						
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- subset(bam.covm, HCOV>700/9719)
	tmp			<- bam.covm[, list(	N=length(SID), 
					XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
					HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
			), by='COV_MIN']
	tmp[, READL:=375L]
	ans			<- rbind(ans, tmp)	
	#	make table
	ans[, XCOV_LABEL:= paste0(round(XCOV_MEAN,d=0),'x ( ', round(XCOV_QL, d=1),'x - ',round(XCOV_QU, d=0),'x )')]
	ans[, HCOV_LABEL:= paste0(round(100*HCOV_MEAN,d=1),'% (', round(100*HCOV_QL, d=1),'% - ',round(100*HCOV_QU, d=1),'%)')]
	ans[, P:= paste0( round(100*N / ans[COV_MIN==1 & READL==1, N], d=1), '%')]
	set(ans, NULL, 'COV_MIN', ans[, factor(COV_MIN, levels=c(1,10,20,30,50,100), labels=c('1X','10X','20X','30X','50X','100X'))])
	ans			<- subset(ans, select=c(COV_MIN, READL, N, P, HCOV_LABEL, HCOV_MEAN))
	setkey(ans, COV_MIN, READL)
	write.csv(ans, row.names=FALSE, file=paste0(outfile.base,'NGSoutput_info.csv'))
	
	ggplot(subset(ans, READL>1), aes(x=READL, y=N/max(N), group=COV_MIN, pch=COV_MIN)) + 
			geom_line(colour='grey50') +
			geom_point() + 			 
			scale_y_continuous(label=scales:::percent) +
			theme_bw() + theme(legend.position='bottom') +
			labs(x='\nmin read length', y='subjects retained\n',pch='min sequencing\ndepth')
	ggsave(file=paste0(outfile.base,'subjectsretained_by_minreads.pdf'), w=3.5, h=6)
	ggplot(subset(ans, READL>1), aes(x=READL, y=HCOV_MEAN*9719, group=COV_MIN, pch=COV_MIN)) + 
			geom_line(colour='grey50') +
			geom_point() + 			 
			scale_y_continuous() +
			theme_bw() + theme(legend.position='bottom') +
			labs(x='\nmin read length', y='average coverage of HIV-1 genome\n(nt)',pch='min sequencing\ndepth')
	ggsave(file=paste0(outfile.base,'coverage_by_minreads.pdf'), w=3.5, h=6)
	
	
	#	roughly where are these 250 bp reads?
	bam.cov250.30	<- subset(bam.cov250, COV>=30)
	bam.cov250.30[, RID:= SID]	
	tmp				<- bam.cov250.30[, list(SUM_REP=sum(REP)), by=c('SID','RID')]
	tmp				<- subset(tmp, SUM_REP>700)
	tmp				<- tmp[, list(SID=SID[which.max(SUM_REP)]), by='RID']
	bam.cov250.30	<- merge(tmp, bam.cov250.30, by=c('RID','SID'))
	setkey(bam.cov250.30, RID, POS)
	tmp				<- bam.cov250.30[, 	{
				z	<- rep(COV,REP)
				list(COV=z, POS2=POS+seq_along(z)-1L)
			}, by=c('SID','RID','POS')]
	tmp			<- tmp[, list(N=length(RID)), by='POS2']
	ggplot(tmp, aes(x=POS2, y=N)) + geom_area() +
			theme_bw() +
			scale_x_continuous(expand=c(0,0)) +
			scale_y_continuous(expand=c(0,0)) +
			#coord_cartesian(ylim=c(0,NA)) +
			labs(x='\nposition on HIV-1 genome', y='individuals\nwith NGS reads at position\n') +
			theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
	ggsave(file=paste0(outfile.base,'horizontal_coverage_histogram_minNGSoutput_reads250.pdf'), w=8, h=2)
	
	#	who are we losing due to some drop off in quality?
	#
	df			<- copy(bam.cov350)
	bam.covm	<- do.call('rbind',lapply(c(1,10,20,30,50,100), function(x)
					{
						bam.covm	<- subset(df, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	dcast.data.table(bam.covm, SID~COV_MIN, value.var='HCOV')
}

pty.pipeline.phyloscanner.180605.MunichCluster.process<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	
	# set up pty.runs file
	# with cols	'SAMPLE_ID','RENAME_ID','UNIT_ID'
	if(0)	
	{
		pty.runs		<- data.table(SAMPLE_ID=list.files('~/duke/2018_MunichCluster/Data', pattern='bam$'))
		set(pty.runs, NULL, 'SAMPLE_ID', pty.runs[, gsub('\\.bam','',SAMPLE_ID)])
		pty.runs[, RENAME_ID:= paste0('MC-',SAMPLE_ID)]
		pty.runs[, UNIT_ID:= paste0('MC-',SAMPLE_ID)]
		pty.runs[, PTY_RUN:=1]		
		outfile	<- '~/duke/2018_MunichCluster/Data/MunichCluster_180605.rda'
		save(pty.runs, file=outfile)
	}
	# set up pty.runs file -- version 2
	# with cols	'SAMPLE_ID','RENAME_ID','UNIT_ID'
	if(0)	
	{
		pty.runs		<- data.table(SAMPLE_ID=list.files('~/duke/2018_MunichCluster/Data_180607', pattern='bam$'))
		set(pty.runs, NULL, 'SAMPLE_ID', pty.runs[, gsub('\\.bam','',SAMPLE_ID)])
		pty.runs[, RENAME_ID:= paste0('MC-',SAMPLE_ID)]
		pty.runs[, UNIT_ID:= paste0('MC-',gsub('_INT|_PRRT','',SAMPLE_ID))]		
		pty.runs[, PTY_RUN:=1]		
		outfile	<- '~/duke/2018_MunichCluster/Data/MunichCluster_180607.rda'
		save(pty.runs, file=outfile)
	}
	# set up pty.runs file -- version 3
	# with cols	'SAMPLE_ID','RENAME_ID','UNIT_ID'
	if(0)	
	{
		pty.runs		<- data.table(SAMPLE_ID=list.files('~/duke/2018_MunichCluster/Data_180618', pattern='bam$'))
		set(pty.runs, NULL, 'SAMPLE_ID', pty.runs[, gsub('\\.bam','',SAMPLE_ID)])
		pty.runs[, RENAME_ID:= paste0('MC-',SAMPLE_ID)]
		pty.runs[, UNIT_ID:= paste0('MC-',gsub('_INT|_PRRT','',SAMPLE_ID))]		
		pty.runs[, PTY_RUN:=1]		
		outfile	<- '~/duke/2018_MunichCluster/Data_180618/MunichCluster_180618.rda'
		save(pty.runs, file=outfile)
	}
	# set up pty.runs file -- version 4
	# with cols	'SAMPLE_ID','RENAME_ID','UNIT_ID'
	# 39 individuals is too many
	if(0)	
	{
		pty.runs		<- data.table(SAMPLE_ID=list.files('~/duke/2018_MunichCluster/Data_180815', pattern='bam$'))
		set(pty.runs, NULL, 'SAMPLE_ID', pty.runs[, gsub('\\.bam','',SAMPLE_ID)])
		tmp	<- pty.runs[, which(!grepl('Control',SAMPLE_ID))]
		set(pty.runs, tmp, 'RENAME_ID', pty.runs[tmp,paste0('MC-',SAMPLE_ID)])
		tmp	<- pty.runs[, which(grepl('Control',SAMPLE_ID))]
		set(pty.runs, tmp, 'RENAME_ID', pty.runs[tmp,gsub('Control','CNTRL-',SAMPLE_ID)])
		pty.runs[, UNIT_ID:= gsub('_INT|_PRRT','',RENAME_ID)]	
		pty.runs[, PTY_RUN:=1]		
		outfile	<- '~/duke/2018_MunichCluster/Data_180815/MunichCluster_180815.rda'
		save(pty.runs, file=outfile)
	}
	#
	#	INPUT ARGS
	if(1)
	{	
		#HOME				<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'
		HOME				<<- '/work/or105/MUNICH'
		in.dir				<- file.path(HOME,"MunichCluster_180815_in")
		work.dir			<- file.path(HOME,"MunichCluster_180815_work")
		out.dir				<- file.path(HOME,"MunichCluster_180815_out")		
		load( file.path(in.dir, 'MunichCluster_180815.rda') )
		print(pty.runs)
		hpc.load			<- "module load intel-suite/2015.1 mpi R/3.3.3 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"
		hpc.nproc			<- 1
		#prog.pty			<- '/Users/Oliver/git/phylotypes/phyloscanner_make_trees.py'
		#pty.data.dir		<- '~/duke/2018_MunichCluster/Data'
		#prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-SSE3 -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-SSE3 -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))
		#prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner_make_trees.py'
		prog.pty			<- '/work/or105/libs/phyloscanner/phyloscanner_make_trees.py'
		#pty.data.dir		<- '/work/or105/MUNICH/data'
		pty.data.dir		<- '/work/or105/MUNICH/data_180815'
		prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-AVX -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-AVX -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))
		pty.select			<- 1		
	}	
	
	#
	# generate read alignments
	if(0)
	{		
		#ptyi		<- seq(800,9150,25)
		ptyi		<- seq(2000,5500,25)		
		#ptyi		<- seq(2475,2525,25)
		pty.c		<- lapply(seq_along(ptyi), function(i)
				{
					pty.args			<- list(	prog.pty=prog.pty, 
							prog.mafft='mafft', 
							prog.raxml=prog.raxml, 
							data.dir=pty.data.dir, 
							work.dir=work.dir, 
							out.dir=out.dir, 
							alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_B.fasta"),
							alignments.root='REF_B_K03455', 
							alignments.pairwise.to='REF_B_K03455',
							window.automatic= '', 
							merge.threshold=1, 
							min.read.count=1, 
							quality.trim.ends=23, 
							min.internal.quality=23, 
							merge.paired.reads=TRUE, 
							no.trees=TRUE, 
							dont.check.duplicates=FALSE,
							dont.check.recombination=TRUE,
							num.bootstraps=1,
							all.bootstrap.trees=TRUE,
							strip.max.len=350, 
							min.ureads.individual=NA, 
							win=c(ptyi[i],ptyi[i]+250,25,250),				 				
							keep.overhangs=FALSE,
							mem.save=0,
							verbose=TRUE,					
							select=pty.select	#of 240
					)											
					pty.c				<- phsc.cmd.phyloscanner.multi(pty.runs, pty.args,regex.ref='\\.fasta$', postfix.sample.id='\\.bam|\\.fasta')
					#cat(pty.c$CMD)
					#stop()
					pty.c[, W_FROM:= ptyi[i]]
					pty.c
				})
		pty.c	<- do.call('rbind', pty.c)	
		tmp		<- data.table(FO=list.files(out.dir, pattern='ptyr.*fasta$', recursive=TRUE, full.names=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FO)))]
		tmp[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FO)))]
		pty.c	<- merge(pty.c, tmp, by=c('PTY_RUN','W_FROM'), all.x=1)		
		pty.c	<- subset(pty.c, is.na(FO))
		#pty.c	<- subset(pty.c, W_FROM==2350)
		print(pty.c, n=1e3)
		#	make array job
		pty.c[, CASE_ID:= seq_len(nrow(pty.c))]
		tmp		<- pty.c[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
		tmp		<- tmp[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]		
		cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=18, hpc.q="pqeelab", hpc.mem="6gb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load, hpc.array=pty.c[, max(CASE_ID)])
		cmd		<- paste(cmd,tmp,sep='\n')
		outfile	<- paste("mc",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(work.dir, outfile, cmd)
		
		#print(pty.c)
		#stop()		 
		#invisible(pty.c[,	{
		#					cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=18, hpc.q="pqeelab", hpc.mem="6gb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd	<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=71, hpc.q=NA, hpc.mem="1850mb",  hpc.nproc=1, hpc.load=hpc.load)
							#cmd	<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=71, hpc.q=NA, hpc.mem="63800mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
		#					cmd		<- paste(cmd,CMD,sep='\n')
		#					cat(cmd)					
		#					outfile	<- paste("mc",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		#					cmd.hpccaller(work.dir, outfile, cmd)
							#stop()
		#				}, by=c('PTY_RUN','W_FROM')])
		quit('no')
	}
	#
	# generate trees
	if(0)
	{
		#HOME		<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'	
		hpc.load	<- "module load intel-suite/2015.1 mpi raxml/8.2.9"
		hpc.select	<- 1; hpc.nproc<- 1; 	hpc.walltime<- 71; hpc.mem<-"6gb"; hpc.q<- "pqeelab"
		#hpc.select	<- 1; hpc.nproc<- 16; 	hpc.walltime<- 23; hpc.mem<-"124gb"; hpc.q<- NA
		#raxml.pr	<- ifelse(hpc.nproc==1, 'raxmlHPC-SSE3', 'raxmlHPC-PTHREADS-SSE3')	
		raxml.pr	<- ifelse(hpc.nproc==1, 'raxmlHPC-AVX','raxmlHPC-PTHREADS-AVX')
		raxml.args	<- ifelse(hpc.nproc==1, '-m GTRCAT --HKY85 -p 42 -o REF_B_K03455', paste0('-m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42 -o REF_B_K03455'))
		in.dir		<- file.path(HOME,"MunichCluster_180815_out")
		out.dir		<- in.dir
		work.dir	<- file.path(HOME,"MunichCluster_180815_work")
		
		infiles	<- data.table(FI=list.files(in.dir, pattern='fasta$', full.names=TRUE, recursive=TRUE))
		infiles[, FO:= gsub('fasta$','tree',FI)]
		infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FI)))]
		infiles[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FI)))]		
		tmp		<- data.table(FT=list.files(out.dir, pattern='^ptyr.*tree$', full.names=TRUE, recursive=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FT)))]
		tmp[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FT)))]
		infiles	<- merge(infiles, tmp, by=c('PTY_RUN','W_FROM'), all.x=1)
		infiles	<- subset(infiles, is.na(FT))	
		setkey(infiles, PTY_RUN, W_FROM)			
		#infiles	<- subset(infiles, PTY_RUN>=240 & PTY_RUN<300)
		#infiles	<- subset(infiles, PTY_RUN>=186 & PTY_RUN<239)
		#print(infiles)	 		
		
		#	make raxml run
		pty.c	<- infiles[, list(CMD=cmd.raxml(FI, outfile=FO, pr=raxml.pr, pr.args=raxml.args)), by=c('PTY_RUN','W_FROM')]		
		print(pty.c, n=1e3)
		
		#	make array job
		pty.c[, CASE_ID:= seq_len(nrow(pty.c))]
		tmp		<- pty.c[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
		tmp		<- tmp[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]		
		cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=hpc.load, hpc.array=pty.c[, max(CASE_ID)])
		cmd		<- paste(cmd,tmp,sep='\n')
		outfile	<- paste("mct",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(work.dir, outfile, cmd)				
	}
	#
	#	combine all the data	
	if(0)
	{
		indirs 	<- '/Users/Oliver/duke/tmp/ptyr143_trees'
		indirs	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/BEEHIVE_67_180302_out'
		indirs	<- '/work/or105/MUNICH/MunichCluster_180815_out'
		#
		indirs	<- list.files(indirs, pattern='^ptyr[0-9]+_trees$', full.names=TRUE)
		allwin	<- data.table(W_FROM=seq(2050,5075,25))
		#allwin	<- data.table(W_FROM=seq(800,9150,25))		
		#allwin	<- data.table(W_FROM=seq(800,9050,125))
		#indirs	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo/ptyr97_trees'
		for(i in seq_along(indirs))
		{
			indir	<- indirs[i]
			pty.run	<- as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(indir)))
			#	check if we have all fasta and tree files
			infiles	<- data.table(F=list.files(indir,pattern='ptyr.*fasta$',full.names=TRUE))
			infiles[, W_FROM:= as.integer(gsub('.*_InWindow_([0-9]+)_.*','\\1',basename(F)))]
			infiles	<- merge(allwin, infiles, by='W_FROM', all.x=1)			 					
			missfs	<- subset(infiles, is.na(F))[, W_FROM]
			if(length(missfs))
				cat('\nIn',indir,'Found missing fasta files for',paste(missfs,collapse=', '))
			infiles	<- data.table(F=list.files(indir,pattern='ptyr.*tree$',full.names=TRUE))
			infiles[, W_FROM:= as.integer(gsub('.*_InWindow_([0-9]+)_.*','\\1',basename(F)))]
			infiles	<- merge(allwin, infiles, by='W_FROM', all.x=1)
			misstrs	<- subset(infiles, is.na(F))[, W_FROM]
			if(length(misstrs))
				cat('\nIn',indir,'Found missing tree files for',paste(misstrs,collapse=', '))
			zipit	<- 0
			if(!length(missfs) & !length(misstrs))
			{
				cat('\nIn',indir,'Found all fasta and tree files')
				zipit	<- 1
			}				
			if(!length(setdiff(misstrs,missfs)))
			{ 
				cat('\nIn',indir,'Found all tree files for which there is a fasta file')
				zipit	<- 1
			}	
			#zipit	<- 1
			#
			if(zipit)
			{			
				cat('\nProcess',indir)
				#	first combine all zip files into ptyrXXX_otherstuff.zip
				infiles	<- data.table(F=list.files(indir,pattern='zip$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				stopifnot(!nrow(subset(infiles, is.na(PTY_RUN))))		
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_otherstuff.zip'))]
				cat('\nZip to file', tmp,'...\n')
				suppressWarnings(invisible( infiles[, list(RTN= unzip(F, overwrite=FALSE, exdir=file.path(indir,'tmp42'))), by='F'] ))
				invisible( infiles[, list(RTN= file.remove(F)), by='F'] )
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	now zip fasta files
				infiles	<- data.table(F=list.files(indir,pattern='ptyr.*fasta$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				invisible( infiles[, file.rename(F, file.path(indir,'tmp42',basename(F))), by='F'] )
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_trees_fasta.zip'))]
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	now zip tree files
				infiles	<- data.table(F=list.files(indir,pattern='ptyr.*tree$',full.names=TRUE))
				infiles[, PTY_RUN:= gsub('^ptyr([0-9]+)_.*','\\1',basename(F))]
				invisible( infiles[, file.rename(F, file.path(indir,'tmp42',basename(F))), by='F'] )
				tmp		<- infiles[1, file.path(dirname(F),paste0('ptyr',PTY_RUN,'_trees_newick.zip'))]		
				invisible( zip( tmp, file.path(indir,'tmp42'), flags = "-umr9XTjq") )
				#	remove tmp dir
				invisible( file.remove( file.path(indir,'tmp42') ) )
				#	move one level down
				infiles	<- data.table(F=list.files(indir, full.names=TRUE))
				invisible( infiles[, file.rename(F, file.path(dirname(indir),basename(F))), by='F'] )
				cat('\nDone',indir)
			}
			#if(!length(misstrs))
			if(zipit)
				invisible(unlink(indir, recursive=TRUE))
			#	expand again if asked to
			#if(length(misstrs))
			if(0)
			{
				cat('\nExtract',file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_fasta.zip')))
				unzip(file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_fasta.zip')), junkpaths=TRUE, exdir=indir)
				cat('\nExtract',file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_newick.zip')))
				unzip(file.path(dirname(indir),paste0('ptyr',pty.run,'_trees_newick.zip')), junkpaths=TRUE, exdir=indir)
			}
		}					
	}
	#
	#	Input args 2
	if(1)	
	{	
		hpc.select	<- 1; hpc.nproc<- 1; 	hpc.walltime<- 71; hpc.mem<-"6gb"; hpc.q<- "pqeelab"
		in.dir				<- '/work/or105/MUNICH/MunichCluster_180815_out'
		out.dir				<- in.dir
		work.dir			<- '/work/or105/MUNICH/MunichCluster_180815_work' 			
		dir.create(out.dir, showWarnings=FALSE)
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft=NA, 
				prog.raxml=NA, 
				data.dir=NA, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_B.fasta"),
				alignments.root='REF_B_K03455', 
				alignments.pairwise.to='REF_B_K03455',				
				bl.normalising.reference.file=system.file(package="Phyloscanner.R.utilities", "data", "hiv.hxb2.norm.constants.rda"),
				bl.normalising.reference.var='MEDIAN_PWD',														
				window.automatic= '', 
				merge.threshold=1, 
				min.read.count=1, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				merge.paired.reads=TRUE, 
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				dont.check.recombination=TRUE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				win=c(2000,5750,25,250), 				
				keep.overhangs=FALSE,
				use.blacklisters=c('ParsimonyBasedBlacklister'),	#,'DownsampleReads'
				tip.regex='^(.*)-[0-9]-[0-9]+_read_([0-9]+)_count_([0-9]+)$',
				roguesubtree.kParam=20,
				roguesubtree.prop.threshold=0,
				roguesubtree.read.threshold=20,
				#dwns.maxReadsPerPatient=50,	
				multifurcation.threshold=1e-5,
				split.rule='s',
				split.kParam=20,
				split.proximityThreshold=0,
				split.readCountsMatterOnZeroBranches=TRUE,
				split.pruneBlacklist=FALSE,
				trms.allowMultiTrans=TRUE,
				pw.trmw.min.reads=30,									
				pw.trmw.min.tips=1,
				pw.trmw.close.brl=0.025,
				pw.trmw.distant.brl=0.05,
				pw.prior.keff=2,
				pw.prior.neff=3,
				pw.prior.keff.dir=2,
				pw.prior.neff.dir=3,				
				pw.prior.calibrated.prob=0.6,
				mem.save=0,
				verbose=TRUE,				
				select=NA 
		)	
		save(pty.args, file=file.path(out.dir, 'pty.args.rda'))
		
		pty.c	<- data.table(FILE_BAM=list.files(in.dir, pattern='_patients.txt', full.names=TRUE))
		pty.c[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_patients.txt','',basename(FILE_BAM))))]				
		tmp		<- data.table(FILE_TRMW=list.files(out.dir, pattern='_pairwise_relationships.rda', full.names=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_pairwise_relationships.rda','',basename(FILE_TRMW))))]
		pty.c	<- merge(pty.c, tmp, by='PTY_RUN', all.x=1)
		pty.c	<- subset(pty.c, is.na(FILE_TRMW))
		setkey(pty.c, PTY_RUN)		
		pty.c	<- pty.c[, { 
					#FILE_BAM<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr1_bam.txt'
					#FILE_BAM<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr1_bam.txt'
					#cat('\n',FILE_BAM)
					prefix.infiles	<- gsub('patients.txt','',FILE_BAM)
					print(prefix.infiles)
					cmd				<- phsc.cmd.phyloscanner.one.resume(prefix.infiles, pty.args)
					list(CMD=cmd) 
				}, by='PTY_RUN']		
		pty.c[1,cat(CMD)]		
		#stop()		
		invisible(pty.c[,	{
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=hpc.load)							
							cmd			<- paste(cmd,'cd $TMPDIR',sep='\n')
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile		<- paste("scRAr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')														
							cmd.hpccaller(work.dir, outfile, cmd)														
						}, by='PTY_RUN'])		
		quit('no')		
	}
}

pty.pipeline.phyloscanner.170301.secondstage<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	#
	#	INPUT ARGS PLATFORM
	#	
	if(1)
	{	
		#HOME				<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'
		in.dir				<- file.path(HOME,"RakaiAll_input_170301")
		work.dir			<- file.path(HOME,"RakaiAll_work_170704")
		out.dir				<- file.path(HOME,"RakaiAll_output_170704_w250_s20_p35_stagetwo")		
		#load( file.path(in.dir, 'Rakai_phyloscanner_170301.rda') )
		#load( file.path(in.dir, 'Rakai_phyloscanner_170301_stagetwo.rda') )
		load( file.path(in.dir, 'Rakai_phyloscanner_170704_stagethree.rda') )
		#print(pty.runs)
		setnames(pty.runs, c('SID','RENAME_SID','RID'), c('SAMPLE_ID','RENAME_ID','UNIT_ID'))
		hpc.load			<- "module load intel-suite/2015.1 mpi R/3.3.3 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"
		hpc.nproc			<- 1							
		prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner_make_trees.py'
		pty.data.dir		<- '/work/or105/PANGEA_mapout/data'
		#prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-SSE3 -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-SSE3 -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))
		prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-AVX -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-AVX -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))
		pty.select			<- NA		
	}	
	if(0)
	{
		prog.pty			<- '/Users/Oliver/git/phylotypes/phyloscanner.py'
		pty.data.dir		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'		
	}			
	#
	#	INPUT ARGS PHYLOSCANNER RUN
	#	
	if(0)	#run read alignments+trees+Rscripts
	{				
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft='mafft', 
				prog.raxml=prog.raxml, 
				data.dir=pty.data.dir, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX_v2.fasta"),
				alignments.root='REF_CPX_AF460972', 
				alignments.pairwise.to='REF_B_K03455',
				bl.normalising.reference.file=system.file(package="Phyloscanner.R.utilities", "data", "hiv.hxb2.norm.constants.rda"),
				bl.normalising.reference.var='MEDIAN_PWD',														
				window.automatic= '', 
				merge.threshold=0, 
				min.read.count=1, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				merge.paired.reads=TRUE, 
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				win=c(800,9400,25,250),
				#win=c(800,9400,125,250), 				
				keep.overhangs=FALSE,
				use.blacklisters=c('ParsimonyBasedBlacklister','DownsampleReads'),
				tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$',				
				roguesubtree.kParam=20,
				roguesubtree.prop.threshold=0,
				roguesubtree.read.threshold=20,
				dwns.maxReadsPerPatient=50,	
				multifurcation.threshold=1e-5,
				split.rule='s',
				split.kParam=20,
				split.proximityThreshold=0.035,
				split.readCountsMatterOnZeroBranches=TRUE,
				pw.trmw.min.reads=20,									
				pw.trmw.min.tips=1,
				pw.trmw.close.brl=0.035,
				pw.trmw.distant.brl=0.08,
				pw.prior.keff=2,
				pw.prior.neff=3,
				pw.prior.calibrated.prob=0.5,
				mem.save=0,
				verbose=TRUE,					
				select=pty.select	#of 240
		)		 
	}	
	if(1)	#run read alignments
	{				
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft='mafft', 
				prog.raxml=prog.raxml, 
				data.dir=pty.data.dir, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX_v2.fasta"),
				alignments.root='REF_CPX_AF460972', 
				alignments.pairwise.to='REF_B_K03455',
				window.automatic= '', 
				merge.threshold=0, 
				min.read.count=1, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				merge.paired.reads=TRUE, 
				no.trees=TRUE, 
				dont.check.duplicates=FALSE,
				dont.check.recombination=TRUE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				win=c(800,9400,25,250),				 				
				keep.overhangs=FALSE,
				mem.save=0,
				verbose=TRUE,					
				select=pty.select	#of 240
		)		 
	}	
	#
	#	RUN PHYLOSCANNER
	#
	if(1)
	{
		pty.c				<- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)
		#pty.c				<- subset(pty.c, PTY_RUN!=1)
		#pty.c[1,cat(CMD)]		
		invisible(pty.c[,	{
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=998, hpc.q="pqeelab", hpc.mem="5900mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=71, hpc.q=NA, hpc.mem="63850mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=998, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile		<- paste("scRAa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
							#stop()
						}, by='PTY_RUN'])
		quit('no')
	}	
	#
	#	run read alignments, one run per window
	#
	if(0)
	{
		pty.select	<- 1		
		print(pty.select)		
		ptyi		<- seq(800,9150,25) 
		pty.c		<- lapply(seq_along(ptyi), function(i)
				{
					pty.args			<- list(	prog.pty=prog.pty, 
													prog.mafft='mafft', 
													prog.raxml=prog.raxml, 
													data.dir=pty.data.dir, 
													work.dir=work.dir, 
													out.dir=out.dir, 
													alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX_v2.fasta"),
													alignments.root='REF_CPX_AF460972', 
													alignments.pairwise.to='REF_B_K03455',
													window.automatic= '', 
													merge.threshold=0, 
													min.read.count=1, 
													quality.trim.ends=23, 
													min.internal.quality=23, 
													merge.paired.reads=TRUE, 
													no.trees=TRUE, 
													dont.check.duplicates=FALSE,
													dont.check.recombination=TRUE,
													num.bootstraps=1,
													all.bootstrap.trees=TRUE,
													strip.max.len=350, 
													min.ureads.individual=NA, 
													win=c(ptyi[i],ptyi[i]+250,25,250),				 				
													keep.overhangs=FALSE,
													mem.save=0,
													verbose=TRUE,					
													select=pty.select	#of 240
													)											
					pty.c				<- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)
					pty.c[, W_FROM:= ptyi[i]]
					pty.c
				})
		pty.c	<- do.call('rbind', pty.c)			
		tmp		<- data.table(FO=list.files(out.dir, pattern='ptyr.*fasta$', recursive=TRUE, full.names=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FO)))]
		tmp[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FO)))]
		pty.c	<- merge(pty.c, tmp, by=c('PTY_RUN','W_FROM'), all.x=1)
		pty.c	<- subset(pty.c, is.na(FO))		
		print(pty.c)
		#stop()		 
		invisible(pty.c[,	{
							#cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=998, hpc.q="pqeelab", hpc.mem="5900mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=71, hpc.q=NA, hpc.mem="1850mb",  hpc.nproc=1, hpc.load=hpc.load)
							cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=71, hpc.q=NA, hpc.mem="63800mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile		<- paste("snorc",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(work.dir, outfile, cmd)
							#stop()
						}, by=c('PTY_RUN','W_FROM')])
		quit('no')										
	}
	#
	#	resubmit extra long jobs to create read fasta files:
	#   re-run with "--dont-check-recombination"
	#
	if(0) 
	{ 
		infile.xlong	<- '/Users/Oliver/duke/tmp/scRaa_extralong.txt'
		indir.qsub		<- '/Users/Oliver/duke/tmp'
		infile.xlong	<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/scRaa_extralong.txt'
		indir.qsub		<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_work_170301'
		df				<- data.table(FQ=gsub(' ','',readLines(infile.xlong)))
		set(df, NULL, 'FQ', df[, file.path(indir.qsub,FQ)])
		#df[, FO:= gsub('scRaa','smt1',FQ)]
		df[, FO:= gsub('scRaa','snorc',FQ)]
		df[, {
					#FQ<- '/Users/Oliver/duke/tmp/scRaa.Fri_May__5_23:38:32_2017.qsub'
					z	<- paste(readLines(FQ), collapse='\n')
					#z	<- gsub('--merging-threshold 0','--merging-threshold 1',z)
					z	<- gsub('--merging-threshold 0','--merging-threshold 0 --dont-check-recombination',z)
					cat(z, file=FO)
					cat(z)
					cmd	<- paste("qsub",FO,'\n')					
					cat( system(cmd, intern=TRUE) )
					#stop()					
					Sys.sleep(1)
					NULL					
				}, by='FQ']
	}
	if(0)	#check failing runs
	{		
		df		<- data.table(F=readLines('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/all_files.txt'))
		df		<- data.table(F=list.files('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun', pattern='^ptyr',full.names=TRUE))
		df[, PTY_RUN:= as.integer(gsub('ptyr([0-9]+)_.*','\\1',basename(F)))]
		tmp	<- df[, list(DONE=any(grepl('pairwise',F))), by='PTY_RUN']
		setdiff(1:240, tmp[, PTY_RUN])
		tmp	<- subset(tmp, !DONE)
		merge(tmp, subset(pty.runs, SID=='15080_1_28', c('PTY_RUN','SID')), by='PTY_RUN', all=TRUE)
		
		subset(tmp, !DONE)[, cat(paste(sort(PTY_RUN), collapse='","'))]
	}
}

pty.pipeline.phyloscanner.170301.thirdstage<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	#
	#	INPUT ARGS PLATFORM
	#	
	if(1)
	{	
		#HOME				<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'
		in.dir				<- file.path(HOME,"RakaiAll_input_170301")
		work.dir			<- file.path(HOME,"RakaiAll_work_170301")
		out.dir				<- file.path(HOME,"RakaiAll_output_170301_w250_s20_p35_stagethree")		
		load( file.path(in.dir, 'Rakai_phyloscanner_170301_stagethree.rda') )
		setnames(pty.runs, c('SID','RENAME_SID','RID'), c('SAMPLE_ID','RENAME_ID','UNIT_ID'))
		hpc.load			<- "module load intel-suite/2015.1 mpi R/3.3.2 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"
		hpc.nproc			<- 1							
		prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner.py'
		pty.data.dir		<- '/work/or105/PANGEA_mapout/data'
		#prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-SSE3 -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-SSE3 -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))
		prog.raxml			<- ifelse(hpc.nproc==1, '"raxmlHPC-AVX -m GTRCAT --HKY85 -p 42"', paste('"raxmlHPC-PTHREADS-AVX -m GTRCAT --HKY85 -T ',hpc.nproc,' -p 42"',sep=''))
		pty.select			<- NA		
	}	
	if(0)
	{
		prog.pty			<- '/Users/Oliver/git/phylotypes/phyloscanner.py'
		pty.data.dir		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'		
	}			
	#
	#	INPUT ARGS PHYLOSCANNER RUN
	#	
	if(1)	#run read alignments+trees+Rscripts
	{				
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft='mafft', 
				prog.raxml=prog.raxml, 
				data.dir=pty.data.dir, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX_v2.fasta"),
				alignments.root='REF_CPX_AF460972', 
				alignments.pairwise.to='REF_B_K03455',
				bl.normalising.reference.file=system.file(package="Phyloscanner.R.utilities", "data", "hiv.hxb2.norm.constants.rda"),
				bl.normalising.reference.var='MEDIAN_PWD',														
				window.automatic= '', 
				merge.threshold=0, 
				min.read.count=1, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				merge.paired.reads=TRUE, 
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				dont.check.recombination=TRUE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				win=c(800,9400,25,250),
				#win=c(800,9400,125,250), 				
				keep.overhangs=FALSE,
				use.blacklisters=c('ParsimonyBasedBlacklister','DownsampleReads'),
				tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$',				
				roguesubtree.kParam=20,
				roguesubtree.prop.threshold=0,
				roguesubtree.read.threshold=20,
				dwns.maxReadsPerPatient=50,	
				multifurcation.threshold=1e-5,
				split.rule='s',
				split.kParam=20,
				split.proximityThreshold=0.035,
				split.readCountsMatterOnZeroBranches=TRUE,
				pw.trmw.min.reads=20,									
				pw.trmw.min.tips=1,
				pw.trmw.close.brl=0.035,
				pw.trmw.distant.brl=0.08,
				pw.prior.keff=2,
				pw.prior.neff=3,
				pw.prior.calibrated.prob=0.5,
				mem.save=0,
				verbose=TRUE,					
				select=pty.select	#of 8
		)		 
	}		
	#
	#	RUN PHYLOSCANNER
	#
	if(1)
	{
		pty.c				<- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)
		#pty.c				<- subset(pty.c, PTY_RUN!=1)
		#pty.c[1,cat(CMD)]		
		invisible(pty.c[,	{
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=998, hpc.q="pqeelab", hpc.mem="5900mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=71, hpc.q=NA, hpc.mem="7850mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=998, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile		<- paste("scRAa",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
							#stop()
						}, by='PTY_RUN'])
		quit('no')
	}	
	if(0)	#check failing runs
	{		
		df		<- data.table(F=readLines('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/all_files.txt'))
		df		<- data.table(F=list.files('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun', pattern='^ptyr',full.names=TRUE))
		df[, PTY_RUN:= as.integer(gsub('ptyr([0-9]+)_.*','\\1',basename(F)))]
		tmp	<- df[, list(DONE=any(grepl('pairwise',F))), by='PTY_RUN']
		setdiff(1:240, tmp[, PTY_RUN])
		tmp	<- subset(tmp, !DONE)
		merge(tmp, subset(pty.runs, SID=='15080_1_28', c('PTY_RUN','SID')), by='PTY_RUN', all=TRUE)
		
		subset(tmp, !DONE)[, cat(paste(sort(PTY_RUN), collapse='","'))]
	}
}

pty.pipeline.phyloscanner.170301.secondstage.rerun<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	#
	#	INPUT ARGS PLATFORM
	#	
	if(1)
	{	
		#HOME				<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'								
		hpc.load			<- "module load intel-suite/2015.1 mpi R/3.3.3 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"
		hpc.nproc			<- 1		
		#in.dir				<- file.path(HOME,"RakaiAll_output_170301_w250_s20_p35_stagetwo")		
		#out.dir				<- file.path(HOME,"RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun")
		#out.dir				<- file.path(HOME,"RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun34")		
		#out.dir				<- file.path(HOME,"RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun23")
		#work.dir			<- file.path(HOME,"RakaiAll_work_170301")
		#prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner.py'
		in.dir				<- file.path(HOME,'RakaiAll_output_170704_w250_s20_p35_stagetwo')
		#in.dir				<- file.path(HOME,'RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun23')		
		out.dir				<- file.path(HOME,"RakaiAll_output_170704_w250_s20_p25_d50_stagetwo_rerun23_min30_adj_chain_mean_rogthr20")
		work.dir			<- file.path(HOME,"RakaiAll_work_170704")
		prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner_make_trees.py'		
		#prog.pty			<- '/Users/Oliver/git/phylotypes/phyloscanner_make_trees.py.py'
	}	
	if(0)
	{
		prog.pty			<- '/Users/Oliver/git/phylotypes/phyloscanner.py'
		in.dir				<- '/Users/Oliver/duke/tmp'
		out.dir				<- '/Users/Oliver/duke/tmp'
		work.dir			<- '/Users/Oliver/duke/tmp' 
	}			
	#
	#	INPUT ARGS PHYLOSCANNER RERUN
	#	
	if(1)	
	{	
		dir.create(out.dir, showWarnings=FALSE)
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft=NA, 
				prog.raxml=NA, 
				data.dir=NA, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX_v2.fasta"),
				alignments.root='REF_CPX_AF460972', 
				alignments.pairwise.to='REF_B_K03455',
				bl.normalising.reference.file=system.file(package="Phyloscanner.R.utilities", "data", "hiv.hxb2.norm.constants.rda"),
				bl.normalising.reference.var='MEDIAN_PWD',														
				window.automatic= '', 
				merge.threshold=0, 
				min.read.count=1, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				merge.paired.reads=TRUE, 
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				dont.check.recombination=TRUE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				win=c(800,9400,25,250), 				
				keep.overhangs=FALSE,
				use.blacklisters=c('ParsimonyBasedBlacklister','DownsampleReads'),
				tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$',
				roguesubtree.kParam=20,
				roguesubtree.prop.threshold=0,
				roguesubtree.read.threshold=20,
				dwns.maxReadsPerPatient=50,	
				multifurcation.threshold=1e-5,
				split.rule='s',
				split.kParam=20,
				split.proximityThreshold=0,
				split.readCountsMatterOnZeroBranches=TRUE,
				split.pruneBlacklist=FALSE,
				trms.allowMultiTrans=TRUE,
				pw.trmw.min.reads=30,									
				pw.trmw.min.tips=1,
				pw.trmw.close.brl=0.025,
				pw.trmw.distant.brl=0.05,
				pw.prior.keff=2,
				pw.prior.neff=3,
				pw.prior.keff.dir=2,
				pw.prior.neff.dir=3,				
				pw.prior.calibrated.prob=0.66,
				mem.save=0,
				verbose=TRUE,				
				select=NA 
				)	
		save(pty.args, file=file.path(out.dir, 'pty.args.rda'))
	}
	#
	#	RERUN PHYLOSCANNER
	#	
	if(1)
	{
		pty.c	<- data.table(FILE_BAM=list.files(in.dir, pattern='_patients.txt', full.names=TRUE))
		pty.c[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_patients.txt','',basename(FILE_BAM))))]				
		tmp		<- data.table(FILE_TRMW=list.files(out.dir, pattern='_pairwise_relationships.rda', full.names=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_pairwise_relationships.rda','',basename(FILE_TRMW))))]
		pty.c	<- merge(pty.c, tmp, by='PTY_RUN', all.x=1)
		pty.c	<- subset(pty.c, is.na(FILE_TRMW))
		setkey(pty.c, PTY_RUN)		
		pty.c	<- pty.c[, { 
					#FILE_BAM<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr1_bam.txt'
					#FILE_BAM<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr1_bam.txt'
					#cat('\n',FILE_BAM)
					prefix.infiles	<- gsub('patients.txt','',FILE_BAM)
					print(prefix.infiles)
					cmd				<- phsc.cmd.phyloscanner.one.resume(prefix.infiles, pty.args)
					list(CMD=cmd) 
				}, by='PTY_RUN']		
		#pty.c[1,cat(CMD)]		
		#stop()
		tmp		<- subset(pty.c, !PTY_RUN%in%c(11,94,208,139, 164, 177, 315, 60, 104, 116, 120, 136, 141, 144, 148, 151, 156, 165, 170, 178, 185, 218, 273, 319))
		invisible(tmp[,	{					
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=15, hpc.q="pqeelab", hpc.mem="6gb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)							
							#cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=15, hpc.q="pqeph", hpc.mem="3.6gb",  hpc.nproc=1, hpc.load=hpc.load)
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=15, hpc.q=NA, hpc.mem="96gb",  hpc.nproc=1, hpc.load=hpc.load)
							cmd			<- paste(cmd,'cd $TMPDIR',sep='\n')
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile		<- paste("scRAr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
							#stop()
						}, by='PTY_RUN'])
		tmp		<- subset(pty.c, PTY_RUN%in%c(11,94,208,139, 164, 177, 315, 60, 104, 116, 120, 136, 141, 144, 148, 151, 156, 165, 170, 178, 185, 218, 273, 319))
		invisible(tmp[,	{					
							#cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=15, hpc.q="pqeelab", hpc.mem="5.8gb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)							
							#cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=15, hpc.q="pqeph", hpc.mem="3.6gb",  hpc.nproc=1, hpc.load=hpc.load)
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=15, hpc.q=NA, hpc.mem="96gb",  hpc.nproc=1, hpc.load=hpc.load)
							cmd			<- paste(cmd,'cd $TMPDIR',sep='\n')
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile		<- paste("scRAr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
							#stop()
						}, by='PTY_RUN'])
		
		quit('no')		
	}
	if(0)	#check runs
	{
		df	<- data.table(F=list.files('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_adj_chain_mean', pattern='_pairwise_relationships.rda$', full.names=TRUE))
		df[, PTY_RUN:= as.integer(gsub('ptyr([0-9]+)_.*','\\1',basename(F)))]
		paste(setdiff( 1:346, df[, PTY_RUN]), collapse=', ')		
	}
}

pty.pipeline.phyloscanner.170301.firstbatchofall.rerun<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	#
	#	INPUT ARGS PLATFORM
	#	
	if(1)
	{	
		#HOME				<<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'
		hpc.load			<- "module load intel-suite/2015.1 mpi R/3.3.3 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"		
		work.dir			<- file.path(HOME,"RakaiAll_work_170704")
		in.dir				<- file.path(HOME,"RakaiAll_output_170704_w250_s25_allbatch_sk20_tb_blnormed")
		out.dir				<- file.path(HOME,"RakaiAll_output_170704_w250_s25_allbatch_sk20_tb_blnormed_rerun")		
		prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner_make_trees.py'
		pty.select			<- NA
		pty.select			<- c(1811, 1825)
		#prog.pty			<- '/Users/Oliver/git/phylotypes/phyloscanner.py'				
	}	
	#
	#	INPUT ARGS PHYLOSCANNER RUN
	#	
	if(1)
	{			
		pty.args			<- list(	prog.pty=prog.pty, 
				prog.mafft=NA, 
				prog.raxml=NA, 
				data.dir=NA, 
				work.dir=work.dir, 
				out.dir=out.dir, 
				alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX_v2.fasta"),
				alignments.root='REF_CPX_AF460972', 
				alignments.pairwise.to='REF_B_K03455',
				bl.normalising.reference.file=system.file(package="Phyloscanner.R.utilities", "data", "hiv.hxb2.norm.constants.rda"),
				bl.normalising.reference.var='MEDIAN_PWD',
				tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$',				
				window.automatic= '', 
				merge.threshold=2, 
				min.read.count=1, 				
				merge.paired.reads=TRUE, 
				no.trees=FALSE, 
				dont.check.duplicates=FALSE,
				dont.check.recombination=TRUE,
				num.bootstraps=1,
				all.bootstrap.trees=TRUE,
				strip.max.len=350, 
				min.ureads.individual=NA, 
				quality.trim.ends=23, 
				min.internal.quality=23, 
				win=c(800,9400,125,250), 				
				keep.overhangs=FALSE,	
				use.blacklisters=c('ParsimonyBasedBlacklister','DownsampleReads'),				
				roguesubtree.kParam=20,
				roguesubtree.prop.threshold=0,
				roguesubtree.read.threshold=20,
				dwns.maxReadsPerPatient=50,	
				multifurcation.threshold=1e-5,
				split.rule='s',
				split.kParam=20,
				split.proximityThreshold=0.035,	
				split.readCountsMatterOnZeroBranches=TRUE,
				pw.trmw.min.reads=20,									
				pw.trmw.min.tips=1,
				pw.trmw.close.brl=0.035,
				pw.trmw.distant.brl=0.08,
				pw.prior.keff=2,
				pw.prior.neff=3,
				pw.prior.keff.dir=2,
				pw.prior.neff.dir=3,				
				pw.prior.calibrated.prob=0.66,
				mem.save=0,
				verbose=TRUE,
				select=pty.select				
				)		
	}	
	#
	#	RE-RUN PHYLOSCANNER
	#
	if(1)
	{
		pty.c	<- data.table(FILE_BAM=list.files(in.dir, pattern='_patients.txt', full.names=TRUE))
		pty.c[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_patients.txt','',basename(FILE_BAM))))]				
		tmp		<- data.table(FILE_TRMW=list.files(in.dir, pattern='_pairwise_relationships.rda', full.names=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_pairwise_relationships.rda','',basename(FILE_TRMW))))]
		pty.c	<- merge(pty.c, tmp, by='PTY_RUN', all.x=1)
		pty.c	<- subset(pty.c, is.na(FILE_TRMW))
		setkey(pty.c, PTY_RUN)	
		pty.c	<- subset(pty.c, PTY_RUN%in%pty.select)
		pty.c	<- pty.c[, {
					#FILE_BAM<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr1_bam.txt'
					#FILE_BAM<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr1_bam.txt'
					#cat('\n',FILE_BAM)
					prefix.infiles	<- gsub('patients.txt','',FILE_BAM)
					print(prefix.infiles)
					cmd				<- phsc.cmd.phyloscanner.one.resume(prefix.infiles, pty.args)
					list(CMD=cmd)
				}, by='PTY_RUN']		
		#pty.c[1,cat(CMD)]
		invisible(pty.c[,	{					
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=21, hpc.q="pqeelab", hpc.mem="10900mb",  hpc.nproc=1, hpc.load=hpc.load)							
							#cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=21, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=1, hpc.load=hpc.load)
							#cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=1, hpc.q=NA, hpc.mem="1890mb",  hpc.nproc=1, hpc.load=hpc.load)
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile		<- paste("scRAr",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)							
						}, by='PTY_RUN'])
		quit('no')		
	}
	if(0)	#check runs
	{
		df	<- data.table(F=list.files('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s25_allbatch_sk20_tb_blnormed', pattern='_pairwise_relationships.rda$', full.names=TRUE))
		df[, PTY_RUN:= as.integer(gsub('ptyr([0-9]+)_.*','\\1',basename(F)))]
		paste(setdiff( 1:1891, df[, PTY_RUN]), collapse=', ')
		
		tmp	<- df[, list(DONE=any(grepl('pairwise',F))), by='PTY_RUN']
		tmp	<- subset(tmp, !DONE)
		merge(tmp, subset(pty.runs, SID=='15080_1_28', c('PTY_RUN','SID')), by='PTY_RUN', all=TRUE)
		
		subset(tmp, !DONE)[, cat(paste(sort(PTY_RUN), collapse='","'))]
	}
}

pty.pipeline.phyloscanner.160915.couples.rerun<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)  
	#
	#	INPUT ARGS PLATFORM
	#	
	if(1)
	{		
		#HOME<<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'
		hpc.load			<- "module load intel-suite/2015.1 mpi R/3.3.2 raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"
		hpc.nproc			<- 1
		hpc.mem				<- "5900mb"
		work.dir			<- file.path(HOME,"Rakai_ptinput_160915_couples")
		in.dir				<- file.path(HOME,"Rakai_ptoutput_160915_couples_w270")						
		#out.dir				<- file.path(HOME,"Rakai_ptoutput_170405_couples_w270_d50_st20_trU_rerun")
		#out.dir				<- file.path(HOME,"Rakai_ptoutput_170405_couples_w270_d50_st20_trC_rerun")
		#out.dir				<- file.path(HOME,"Rakai_ptoutput_170405_couples_w270_d50_st20_trB_rerun")
		#out.dir				<- file.path(HOME,"Rakai_ptoutput_170405_couples_w270_d50_st20_trB_blNormed_rerun")
		out.dir				<- file.path(HOME,"Rakai_ptoutput_170426_couples_w250_d50_p35_blNormed_rerun")
		#prog.pty			<- '/Users/Oliver/git/phylotypes/phyloscanner.py'		
		prog.pty			<- '/work/or105/libs/phylotypes_new/phyloscanner.py'						
	}	
	#
	#	INPUT ARGS PHYLOSCANNER RUN
	#	
	if(1) 
	{				
		pty.args			<- list(	prog.pty=prog.pty, 
										prog.mafft=NA, 
										prog.raxml=NA, 
										data.dir=NA, 
										work.dir=work.dir, 
										out.dir=out.dir, 
										alignments.file=system.file(package="Phyloscanner.R.utilities", "HIV1_compendium_AD_B_CPX_v2.fasta"),
										alignments.root='REF_CPX_AF460972',
										bl.normalising.reference.file=system.file(package="Phyloscanner.R.utilities", "data", "hiv.hxb2.norm.constants.rda"),
										bl.normalising.reference.var='MEDIAN_PWD',										
										alignments.pairwise.to='REF_B_K03455',
										window.automatic= '', 
										merge.threshold=0, 
										min.read.count=1, 
										quality.trim.ends=23, 
										min.internal.quality=23, 
										merge.paired.reads=TRUE, 
										no.trees=FALSE, 
										dont.check.duplicates=FALSE,
										num.bootstraps=1,
										all.bootstrap.trees=TRUE,
										strip.max.len=350, 
										min.ureads.individual=NA, 
										win=c(800,9400,25,250), 
										keep.overhangs=FALSE,	
										multifurcation.threshold=1e-5,
										use.blacklisters=c('ParsimonyBasedBlacklister','DownsampleReads'),
										roguesubtree.kParam=20,
										roguesubtree.prop.threshold=0,
										roguesubtree.read.threshold=20,
										dwns.maxReadsPerPatient=50,	
										multifurcation.threshold=1e-5,
										split.rule='s',
										split.kParam=20,
										split.proximityThreshold=0.035,
										split.readCountsMatterOnZeroBranches=TRUE,
										pw.trmw.min.reads=20,									
										pw.trmw.min.tips=1,
										pw.trmw.close.brl=0.035,
										pw.trmw.distant.brl=0.08,
										pw.prior.keff=2,
										pw.prior.neff=3,
										pw.prior.calibrated.prob=0.5,
										mem.save=0,
										verbose=TRUE,
										select=NA)
	}	
	#
	#	RUN PHYLOSCANNER RESUME
	#
	if(1)
	{
		pty.c	<- data.table(FILE_BAM=list.files(in.dir, pattern='_bam.txt', full.names=TRUE))
		pty.c[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_bam.txt','',basename(FILE_BAM))))]
		pty.c	<- subset(pty.c, PTY_RUN!=115)	#what happened to run 115??
		#pty.c	<- subset(pty.c, PTY_RUN%in%c(48))
		tmp		<- data.table(FILE_TRMW=list.files(out.dir, pattern='_trmStatsPerWindow.rda', full.names=TRUE))
		tmp[, PTY_RUN:= as.integer(gsub('ptyr','',gsub('_trmStatsPerWindow.rda','',basename(FILE_TRMW))))]
		pty.c	<- merge(pty.c, tmp, by='PTY_RUN', all.x=1)
		pty.c	<- subset(pty.c, is.na(FILE_TRMW))
		setkey(pty.c, PTY_RUN)		
		pty.c	<- pty.c[, {
					#FILE_BAM<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr1_bam.txt'
					#FILE_BAM<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr1_bam.txt'
					#cat('\n',FILE_BAM)
					prefix.infiles	<- gsub('bam.txt','',FILE_BAM)
					cmd				<- phsc.cmd.phyloscanner.one.resume(prefix.infiles, pty.args)
					list(CMD=cmd)
				}, by='PTY_RUN']		
		#pty.c[1,cat(CMD)]		
		invisible(pty.c[,	{					
							cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=21, hpc.q="pqeelab", hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=hpc.load)							
							#cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=21, hpc.q="pqeph", hpc.mem="3600mb",  hpc.nproc=1, hpc.load=hpc.load)
							cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=3, hpc.q=NA, hpc.mem="1890mb",  hpc.nproc=1, hpc.load=hpc.load)
							cmd			<- paste(cmd,CMD,sep='\n')
							cat(cmd)					
							outfile		<- paste("pty",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
							cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
							stop()
						}, by='PTY_RUN'])
		quit('no')
	}	
	#
	#
	#
	if(1)
	{
		trmw.min.reads					<- 20
		trmw.min.tips					<- 1
		trmw.close.brl 					<- 0.035
		trmw.distant.brl				<- 0.08		
		prior.keff						<- 2
		prior.neff						<- 3
		prior.calibrated.prob			<- 0.5
		relationship.types				<- c('TYPE_PAIR_DI','TYPE_PAIR_TO','TYPE_PAIR_TODI2x2','TYPE_DIR_TODI3','TYPE_DIRSCORE_TODI3','TYPE_PAIR_TODI2','TYPE_PAIR_TODI','TYPE_PAIRSCORE_TODI','TYPE_CHAIN_TODI')
		
		indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_170426_couples_w250_d50_p35_blNormed_rerun'
		indir	<- '/Volumes/exfat/PANGEA/PANGEA_phyloscanner/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed'		
		infiles	<- data.table(FI=list.files(indir, pattern='_trmStatsPerWindow.rda$', full.names=TRUE))
		infiles[, FO:=gsub('_trmStatsPerWindow.rda','_pairwise_relationships.rda',FI)]
		infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FI)))]
		setkey(infiles, PTY_RUN)
		
		invisible(infiles[, {	
								#FI	<- '/Volumes/exfat/PANGEA/PANGEA_phyloscanner/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr1_trmStatsPerWindow.rda'
								cat('\nprocess file', FI,'...')
								tmp		<- load(FI)
								if(length(tmp)!=1)	stop("Unexpected length of loaded objects in infile")
								eval(parse(text=paste("dwin<- copy(",tmp,")",sep='')))		
								eval(parse(text=paste(tmp,"<- NULL",sep='')))
								tmp		<- phsc.pairwise.relationships(dwin, trmw.min.reads, trmw.min.tips, trmw.close.brl, trmw.distant.brl, prior.keff, prior.neff, prior.calibrated.prob, relationship.types)
								dwin	<- copy(tmp$dwin)
								rplkl	<- copy(tmp$rplkl)
								cat('\nwrite to file', FO,'...')								
								save(dwin, rplkl, file=FO)
								NULL
							}, by='PTY_RUN'])
		
	}

}

pty.scan.explore	<- function()
{
	pty.infile	<- "~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_coinfrunsinput.rda"
	indir		<- "~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/coinf_ptoutput_150121"
	outdir		<- indir
	if(0)	#	get phylo statistics
	{
		infiles		<- data.table(FILE=list.files(indir, pattern='examl.rda$'))
		infiles[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]
		setkey(infiles, PTY_RUN)
		cat('\nno examl for runs=',paste( setdiff(seq.int(1,infiles[,max(PTY_RUN)]), infiles[,PTY_RUN]), collapse=',' ))		
		#pty.stat	<- do.call('rbind',lapply(seq_len(nrow(infiles)), function(i)
		for(i in seq_len(nrow(infiles)))	
			for(i in 24:52)
			{
				cat('\nprocess',infiles[i,FILE])
				file		<- file.path(indir,infiles[i,FILE])
				load( file )	#loads "pty.ph"   "ptyfiles"
				pty.stat	<- pty.stat.all.160128(pty.ph, ptyfiles)		
				save(pty.stat, file=gsub('\\.rda','_stat\\.rda', file))
				pty.stat<- coi.div<- coi.lsep<- coi.diff<- NULL
				gc()
				#coi.div
			}
	}
	if(0)	#collect all results
	{
		pty.stat.collect(indir)	
	}
	if(0)	#devel
	{
		pty.stat.file	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/coinf_ptoutput_150121/ptyr_examl_stat.rda'
		load(pty.stat.file)
		
		#	how many windows per individual ?
		tmp	<- pty.stat[, list(FILE_ID_N=length(unique(W_FROM))), by='FILE_ID']
		ggplot(tmp, aes(x=FILE_ID_N)) + geom_histogram()
		#	within host diversity
		tmp	<- pty.stat[, list(WHDA_p50=quantile(WHDA_p50,p=0.5), WHDA_p025=quantile(WHDA_p50,p=0.025), WHDA_p10=quantile(WHDA_p50,p=0.1), WHDA_p90=mean(WHDA_p50,p=0.9), WHDA_p95=quantile(WHDA_p50,p=0.95) ), by='W_FROM']	
		ggplot(pty.stat, aes(x=W_FROM)) +  geom_point(aes(y=WHDA_p50, colour=WHDA_p50, size=COUNT_N), alpha=0.2 ) + 
				geom_ribbon(data=tmp, aes(ymin=WHDA_p10, ymax=WHDA_p90), fill='grey50', alpha=0.7) +
				geom_ribbon(data=tmp, aes(ymin=WHDA_p90, ymax=WHDA_p95), fill='grey80', alpha=0.7) +
				#geom_ribbon(data=tmp, aes(ymin=WHDA_p025, ymax=WHDA_p10), fill='grey80', alpha=0.5) +
				scale_colour_continuous(guide=FALSE) + 
				scale_x_continuous(breaks=seq(0,10e3,500)) +
				geom_line(data=tmp, aes(y=WHDA_p95), colour='red') +
				geom_line(data=tmp, aes(y=WHDA_p90), colour='DarkRed') +
				geom_line(data=tmp, aes(y=WHDA_p50)) + theme_bw() + theme(legend.position='bottom') +
				labs(x='\ngenome position of window start in each run\n(bp)',y='mean within-host diversity\n(patristic distance, subst/site)\n', size='quality trimmed short reads per individual\n(#)')
		ggsave(file=gsub('\\.rda','_whdivallreads\\.pdf',pty.stat.file), w=12, h=6)
		#	get candidates by wh diversity	
		#	not clear where to make cut.off. Just keep them all with decreasing priority, and plot for inspection.
		pty.div	<- merge(subset(pty.stat, select=c(PTY_RUN, W_FROM, W_TO, FILE, FILE_ID, WHDA_p50, TAXA_N, COUNT_N)), subset(tmp, select=c(W_FROM, WHDA_p90, WHDA_p95)), by='W_FROM')
		tmp		<- subset(pty.div, WHDA_p50>WHDA_p95 & TAXA_N>40)[, list(PTY_RUN=PTY_RUN[1], FILE_ID_N=length(W_FROM), WHDA_p50_avg=mean(WHDA_p50)), by='FILE_ID']
		tmp		<- tmp[order(-FILE_ID_N),]	
		write.csv(tmp, file=gsub('\\.rda','_whdivallreads_g90qu\\.csv',pty.stat.file), row.names=FALSE)
		#	160126
		#
		#	R11_GA0053_S71_L001		23		May be a root issue. There s a lot of diversity but no clearly defined separate clade
		#	R1_RES301_S13_L001		14		** Could be! There s a deep split that reproduces
		#	R2_RES142_S40_L001		31		Not clear. Consecutive phylogenies should give similar results -- is this the case??
		#	R1_RES452_S17_L001		27		** There are deep splits, occasionally intermingled with R5_RES414. Again, not super clear.
		#	R1_RES124_S6_L001		23		Not clear. What complicates: x-axes across windows not on same scale.
		#	13554_1_43				37		** This one is super diverse throughout. Not sure what this means though.. at 700 also mixup with another individual
		#	13557_1_93				18
		#	R1_RES282_S12_L001		33
		#	R6_RES025_S58_L001		16
		#	R11_PA0012_S28_L001		14		Root issue?
		
		#
		#	long branches. focus on subtending clades with more than 100 reads
		pty.lsep	<- subset(pty.stat, CL_COUNT_N>100)
		pty.lsep[, THR_WHER:= 0.1]
		set(pty.lsep, pty.lsep[,which(W_FROM>5200)],'THR_WHER',0.3)
		tmp	<- pty.lsep[, list(CL_MX_LOCAL_SEP_p50=quantile(CL_MX_LOCAL_SEP,p=0.5), CL_MX_LOCAL_SEP_p025=quantile(CL_MX_LOCAL_SEP,p=0.025), CL_MX_LOCAL_SEP_p10=quantile(CL_MX_LOCAL_SEP,p=0.1), CL_MX_LOCAL_SEP_p90=mean(CL_MX_LOCAL_SEP,p=0.9), CL_MX_LOCAL_SEP_p95=quantile(CL_MX_LOCAL_SEP,p=0.95) ), by='W_FROM']	
		ggplot(pty.lsep, aes(x=W_FROM)) +  geom_point(aes(y=CL_MX_LOCAL_SEP, colour=CL_MX_LOCAL_SEP, size=COUNT_N), alpha=0.2 ) + 
				geom_ribbon(data=tmp, aes(ymin=CL_MX_LOCAL_SEP_p10, ymax=CL_MX_LOCAL_SEP_p90), fill='grey50', alpha=0.7) +
				geom_ribbon(data=tmp, aes(ymin=CL_MX_LOCAL_SEP_p90, ymax=CL_MX_LOCAL_SEP_p95), fill='grey80', alpha=0.7) +
				#geom_ribbon(data=tmp, aes(ymin=WHDA_p025, ymax=WHDA_p10), fill='grey80', alpha=0.5) +
				scale_colour_continuous(guide=FALSE) + 
				scale_x_continuous(breaks=seq(0,10e3,500)) +
				geom_step(aes(y=THR_WHER), colour='red') +
				geom_step(aes(y=THR_WHER/2), colour='DarkRed') +
				#geom_line(data=tmp, aes(y=CL_MX_LOCAL_SEP_p95), colour='red') +
				#geom_line(data=tmp, aes(y=CL_MX_LOCAL_SEP_p90), colour='DarkRed') +
				geom_line(data=tmp, aes(y=CL_MX_LOCAL_SEP_p50)) + theme_bw() + theme(legend.position='bottom') +
				labs(x='\ngenome position of window start in each run\n(bp)',y='stem length of clades with >100 reads\n(subst/site)\n', size='quality trimmed short reads per individual\n(#)')
		ggsave(file=gsub('\\.rda','_whlsepallreads100\\.pdf',pty.stat.file), w=12, h=6)
		#	get candidates
		tmp		<- subset(pty.lsep, CL_MX_LOCAL_SEP>THR_WHER/2 & TAXA_N>40)[, list(PTY_RUN=PTY_RUN[1], FILE_ID_N=length(W_FROM), CL_MX_LOCAL_SEP_avg=mean(CL_MX_LOCAL_SEP)), by='FILE_ID']
		tmp		<- tmp[order(-FILE_ID_N),]	
		
		
		tmp		<- subset(pty.lsep, CL_MX_LOCAL_SEP>THR_WHER/2 & TAXA_N>40)[, list(FILE_ID_N=length(W_FROM), DIFF_IND= max(DIFF_IND)), by='FILE_ID']
		tmp		<- subset(tmp, FILE_ID_N>2 & DIFF_IND>0)		
		tmp		<- tmp[order(-FILE_ID_N),]	
		write.csv(tmp, file=gsub('\\.rda','_whlsepallreads100\\.csv',pty.stat.file), row.names=FALSE)
		
		tmp2	<- merge(pty.lsep, subset(tmp, select=FILE_ID), by='FILE_ID')
		tmp2	<- subset(tmp2, CL_MX_LOCAL_SEP>THR_WHER/2 & TAXA_N>40)
		
		
		infiles		<- data.table(FILE=list.files(indir, pattern='examl.rda$'))
		infiles[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]		
		phps	<- lapply(seq_len(nrow(tmp2)), function(i)
				{
					#i				<-1
					file			<- file.path(indir, infiles[ PTY_RUN==tmp2[i,PTY_RUN],	FILE])
					load(file)
					ph				<- pty.ph[[ tmp2[i,FILE]  ]]
					max.node.height	<- max(node.depth.edgelength(ph)[1:Ntip(ph)])
					tmp				<- c(c(tmp2[i,FILE_ID],'not characterized'), sort(setdiff(sort(levels(attr(ph,'INDIVIDUAL'))),c(tmp2[i,FILE_ID],'not characterized'))))
					col				<- c('black','grey50',rainbow_hcl(length(tmp)-2, start = 270, end = -30, c=100, l=50))
					names(col)		<- tmp			
					ph.title		<- tmp2[i, paste( FILE_ID,'\nrun=',PTY_RUN,', win=',W_FROM,'-',W_TO, sep='')]
					p				<- ggtree(ph, aes(color=INDIVIDUAL, linetype=TYPE)) + 
							geom_nodepoint(size=ph$node.label/100*3) +
							geom_tiplab(size=1.2,  hjust=-.1) +							 
							scale_color_manual(values=col, guide = FALSE) +											 
							scale_linetype_manual(values=c('target'='solid','filler'='dotted'),guide = FALSE) +
							theme_tree2() +
							theme(legend.position="bottom") + ggplot2::xlim(0, max.node.height*1.3) +
							labs(x='subst/site', title=ph.title)
					p
				})
		file	<- file.path( indir, '160128_coinf_lsep05_txn40_fileidn2_diffind1')
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
		#subset(tmp2, FILE_ID=='13557_1_86', c(PTY_RUN, W_FROM, W_TO, FILE, FILE_ID, DIFF_IND, CL_MX_LOCAL_SEP, CL_AVG_HEIGHT_UNIQUE, CL_AVG_HEIGHT_ALL, CL_TAXA_N, CL_COUNT_N,TAXA_N, COUNT_N))
		#
		#	160127	
		#
		#	13557_1_82			6		Two separate clusters split by others. Promising.
		#	13557_1_86			36		Seems highly divergent. Some branches are v long. Perhaps minority variant??
		#	13557_1_77			34		Two separate clusters split by others in pol; else long branches. Promising. 
		#	13557_1_60			34		Has long branches but not split by others
		#	13557_1_56			11		Little support. Some phylogenies look weird. There are duplicates, occasionally breaking clades.
		#	R1_RES406_S14_L001	19		Several separate clades, broken by others. Promising.
		#	13554_1_24			44		Separate clusters. Occasionally split by others. Promising. 
		
		#	looking for high density trees that have two clades broken by others		--> super-infection?
		#	unclear how to intepret separate clades that are rooted in same individual 	--> co-infection?
		
		subset(pty.lsep, CL_MX_LOCAL_SEP>.5)	#	entry 12 looks like an error!
		ggplot(pty.lsep, aes(x=COUNT_N)) + geom_histogram()
		
		ggsave(file=gsub('\\.rda','_whdivallreads\\.pdf',pty.stat.file), w=12, h=6)
		ggplot(tmp, aes(x=FILE_ID_N)) + geom_histogram(binwidth=1)
		
		ggplot(subset(coi.lsep, CL_TAXA_N>5), aes(x=W_FROM, fill=FILE_ID, colour=FILE_ID, group=FILE_ID)) + 			 
				geom_point(aes(y=CL_MX_LOCAL_SEP/CL_AVG_HEIGHT_ALL, size=CL_COUNT_N)) + theme_bw()
	}
	
}

pty.pipeline.coinfection.statistics<- function()
{	
	require(big.phylo)
	#indir			<- file.path(HOME,"coinf_ptoutput_150121")
	indir			<- file.path(HOME,"coinf_ptoutput_150201")
	work.dir		<- file.path(HOME,"coinf_ptinput")
	hpc.load		<- "module load R/3.3.2"
	resume			<- 1
	
	infiles			<- data.table(FILE=list.files(indir, pattern='preprtr.rda$'))
	infiles[, PTY_RUN:= as.numeric(gsub('ptyr','',sapply(strsplit(FILE,'_'),'[[',1)))]
	cat('\nno examl for runs=',paste( setdiff(seq.int(1,infiles[,max(PTY_RUN)]), infiles[,PTY_RUN]), collapse=',' ))
	tmp			<- data.table(FILE_STAT=list.files(indir, pattern='stat.rda$'))
	tmp[, FILE:= gsub('_stat','',FILE_STAT)]
	infiles		<- merge(infiles, tmp, by='FILE',all.x=1)
	if(resume)
		infiles	<- subset(infiles, is.na(FILE_STAT))
	setkey(infiles, PTY_RUN)
	
	invisible(infiles[, {
						cmd			<- pty.cmd.scan.superinfections(indir, select=paste('ptyr',PTY_RUN,'_',sep=''), references.pattern='REF', run.pattern='ptyr', plot.max.clade=0)													
						cat(cmd)							
						#cmd			<- cmd.hpcwrapper(cmd, hpc.walltime=2, hpc.q="pqeelab", hpc.mem="15600mb",  hpc.nproc=1, hpc.load=hpc.load)
						cmd			<- cmd.hpcwrapper(cmd, hpc.walltime=2, hpc.q=NA, hpc.mem="63000mb",  hpc.nproc=1, hpc.load=hpc.load)										
						outfile		<- paste("pts",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
						cmd.hpccaller(work.dir, outfile, cmd)
						NULL
					}, by='PTY_RUN'])
}

project.readlength.count.all<- function()
{	
	#	for each individual, we first create a file that counts how often a read is in the fasq/bam file (COUNT) and how long the read is (LENGTH)
	#	the file format is: rows->unique read. cols separated by whitespace-> col1 is COUNT, col2 is LENGTH 
	#	if you produce similar files, then the scripts below can be used
	
	#	read the read length files and store in data table 'rl'	
	#
	#	Please note:
	#	1) currently, all files start with 'R[0-9]+_' eg 'R13_' to indicate run 13
	#	2) the above files are also stored in three directories 'RawFastqs','Trimmed_AdaptersPrimersOnly','Trimmed_QualToo'
	#	   these directory names are used to determine the processing stage, which is stored in the column 'TYPE' below
	#
	#	If 1 and 2 are not the case, the following will need a bit of modification :-)
	require(data.table)
	if(0)
	{
		indir		<- '/work/cw109/SAseqs/ReadLengths'
		outfile		<- file.path(HOME,'readlen.rda')
		infiles		<- data.table(FILE=list.files(indir, pattern='dat$', recursive=TRUE))
		cat('\nFound files, n=', nrow(infiles))
		cat('\nWill save to', outfile)
		
		infiles[, TYPE:=dirname(FILE)]
		infiles[, ID:=gsub('\\.dat','',basename(FILE))]
		rl			<- infiles[, {
					stopifnot(length(FILE)==1)				
					#z	<- as.data.table(read.table('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/readlengths/RawFastqs/R13_X84265_S17_L001_2.dat',sep='',stringsAsFactors=FALSE))
					z	<- as.data.table(read.table(file.path(indir,FILE),sep='',stringsAsFactors=FALSE))
					setnames(z, c('V1','V2'),c('COUNT','LEN'))
					z
				}, by=c('TYPE','ID')]
		cat('\nRead reads, n=', nrow(rl))
		save(rl, file=outfile)
		cat('\nSaved to',outfile)
	}
	#
	#	produce various plots on read length
	#
	require(ggplot2)
	if(1)
	{
		indir		<- file.path(HOME,'readlengths')
		file		<- file.path(indir,'readlen.rda')
		load(file)
		rl[, RUN:= regmatches(ID,regexpr('^R[0-9]+', ID))]
		#
		rle			<- rl[, {
					zz	<- rep(LEN,COUNT)
					z	<- ecdf(zz)(seq(0,320,20))
					list(QU=seq(0,320,20), CDF=z, CUM_COUNT=z*length(zz), PDF=c(0,diff(z)))	
				}, by=c('RUN','ID','TYPE')]
		setkey(rle, RUN, ID, TYPE)
		#	files without reads after quality trimming
		tmp			<- unique(rle)[, list(TN=length(TYPE), TYPE=TYPE), by=c('RUN','ID')]
		tmp			<- subset(tmp, TN!=3)
		setkey(tmp, RUN, ID)
		write.csv(tmp, file=gsub('\\.rda','_NoReadsAfterQualTrim\\.csv',file), row.names=FALSE)
		#	total number of reads with len<240 increases for some individuals?
		tmp			<- subset(rle, QU==320 & TYPE!='Trimmed_AdaptersPrimersOnly')
		tmp			<- dcast.data.table(tmp, RUN+ID~TYPE, value.var='CUM_COUNT')
		subset(tmp, RawFastqs<Trimmed_QualToo)	#none, phew
		tmp			<- subset(rle, QU==320)
		ggplot(tmp, aes(y=CUM_COUNT, x=RUN, fill=TYPE)) + geom_boxplot() + 
				theme_bw() + labs(y='total reads per individual\n', x='\nsequence run', fill='processing stage') +
				theme(legend.position='bottom')
		ggsave(file=gsub('\\.rda','_TotalCounts\\.pdf',file), w=7,h=7)
		#
		#	plot by individual
		#		
		ggplot(rle, aes(y=CDF, x=factor(QU), colour=TYPE, group=TYPE)) + geom_line() + 
				theme_bw() + labs(y='cumulative frequency\nin one individual\n', x='\nlength of reads\n(nt)', colour='processing stage') +
				theme(legend.position='bottom') + facet_wrap(~ID,ncol=5)
		ggsave(file=gsub('\\.rda','_byInd\\.pdf',file), w=15,h=300,limitsize = FALSE)
		#	after adaptors or primers are trimmed, the reads should be shorter 
		#	for some individuals, the reads are quite a bit longer!? 
		#	suggesting that many reads are just primer+adaptor?
		tmp			<- subset(rle, QU==100 & TYPE!='Trimmed_QualToo')
		tmp			<- subset(dcast.data.table(tmp, RUN+ID~TYPE, value.var='CDF'), RawFastqs-0.07>Trimmed_AdaptersPrimersOnly)
		setnames(tmp, c('RawFastqs','Trimmed_AdaptersPrimersOnly'), c('RawFastqs_PropOfReads<100bp','Trimmed_AdaptersPrimersOnly_PropOfReads<100bp'))
		write.csv(tmp, file=gsub('\\.rda','_LongerReadsAfterAdaptorsPrimersTrimmed\\.csv',file), row.names=FALSE)		
		#
		#	plot by run
		#		
		ggplot(rle, aes(y=CDF, x=factor(QU), fill=TYPE)) + geom_boxplot() + 
				theme_bw() + labs(y='cumulative frequency\nin one individual\n', x='\nlength of reads\n(nt)', fill='processing stage') +
				theme(legend.position='bottom') + facet_wrap(~RUN,ncol=3)
		ggsave(file=gsub('\\.rda','_byRun\\.pdf',file), w=15,h=15)
		#
		#	plot of cumulative counts by run
		#		
		ggplot(rle, aes(y=CUM_COUNT, x=factor(QU), fill=TYPE)) + geom_boxplot() + 
				theme_bw() + labs(y='total reads with length <x\nin one individual\n', x='\nlength of reads\n(nt)', fill='processing stage') +
				theme(legend.position='bottom') + facet_wrap(~RUN,ncol=3, scales='free')
		ggsave(file=gsub('\\.rda','_byRun_TotalCounts\\.pdf',file), w=15,h=15)	
	}
}

project.TillHIV2.power<- function()
{
	require(TrialSize)
	
	alpha	<- 0.05
	beta	<- 1-0.8
	p1		<- 0.9		#	actual success DBS
	p2		<- 0.95		#	success blood
	k		<- 1		#	equal group size
	diff	<- p1-p2
	margin	<- 0.2		#	delta, noninferiority threshold is -delta
	
	# ( qnorm(1-alpha) + qnorm(1-beta) )^2 / (p2-p1-margin)^2 * (p1*(1-p1)/k + p2*(1-p2))	
	TwoSampleProportion.NIS(0.05,0.2,0.9,0.95,1,-0.05,0.2)
	TwoSampleProportion.NIS(0.05,0.2,0.9,0.95,1,-0.05,0.15)
	TwoSampleProportion.NIS(0.05,0.2,0.9,0.95,1,-0.05,0.1)
	
	TwoSampleProportion.NIS(0.05,0.2,0.85,0.9,1,-0.05,0.2)
	TwoSampleProportion.NIS(0.05,0.2,0.85,0.9,1,-0.05,0.15)
	TwoSampleProportion.NIS(0.05,0.2,0.85,0.9,1,-0.05,0.1)
	
	TwoSampleProportion.NIS(0.05,0.2,0.75,0.8,1,-0.05,0.2)
	TwoSampleProportion.NIS(0.05,0.2,0.75,0.8,1,-0.05,0.15)
	TwoSampleProportion.NIS(0.05,0.2,0.75,0.8,1,-0.05,0.1)
}

project.bam.check.integrity<- function()
{
	require(Rsamtools)
	require(data.table)
	
	#indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/check1855/check1855'
	indir	<- '/work/or105/PANGEA_mapout/data'
	
	infiles	<- data.table(F=list.files(indir, pattern='bam$', full.names=TRUE))
	infiles	<- infiles[, {
				#F<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/check1855/check1855/15080_1_28.bam'
				tmp	<- suppressWarnings( try( scanBam(F, param=ScanBamParam(what=c('qwidth','pos','rname')))[[1]], silent=TRUE ) )
				list(STATUS= ifelse(class(tmp)=='try-error', 'corrupt', 'ok'))			
			}, by='F']
	
	save(infiles, file=file.path(indir, 'bam_status.rda'))
	
	#load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/bam_status.rda')
	#subset(infiles, STATUS!='ok')
	#subset(pty.runs, SID=='15080_1_28')[, cat(paste(PTY_RUN, collapse="','"))]
	#these are 61 runs:
	#c('682','707','732','757','782','807','832','857','882','907','932','957','982','1007','1032','1057','1082','1107','1132','1157','1182','1207','1232','1257','1282','1307','1332','1357','1382','1407','1432','1457','1482','1507','1532','1557','1582','1606','1629','1651','1672','1692','1711','1729','1746','1762','1777','1791','1804','1816','1827','1837','1847','1848','1849','1850','1851','1852','1853','1854','1855')
}

project.bam.read.distribution.readcsvs<- function() 
{
	outfile	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/bam_readdistr.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_bamdistr'
	df		<- data.table(F=list.files(indir, pattern='csv$', full.names=TRUE))
	df[, PTY_RUN:= as.integer(gsub('^bamr([0-9]+)_.*','\\1',basename(F)))]
	paste( setdiff(1:1891, df[,PTY_RUN]), collapse=',' )
		
	db		<- df[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_bamdistr/bamr10_read_distributions.csv.csv'
				ans	<- as.data.table(read.csv(F))
				setnames(ans, c('Bam.file','Size.of.overlapping.read.pair.or.length.of.read.in.non.overlapping.pair','Count'), c('BAM','LEN','COUNT'))
				set(ans, NULL, 'BAM', ans[,gsub('\\.bam','',as.character(BAM))])
				ans
			}, by='PTY_RUN']
	db[, DATA:='Rakai']
	tmp		<- unique(subset(db, select=c(BAM,PTY_RUN)))[, list(PTY_RUN=min(PTY_RUN)), by='BAM']
	db		<- merge(tmp, db, by=c('BAM','PTY_RUN'))
	set(db, NULL, 'PTY_RUN',NULL)
	
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/BEEHIVE_OverlappingInsertSizes'
	df		<- data.table(F=list.files(indir, full.names=TRUE))	
	tmp		<- df[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/BEEHIVE_OverlappingInsertSizes/BEE0001-1.csv'
				ans	<- as.data.table(read.csv(F))
				setnames(ans, c('Bam.file','Size.of.overlapping.read.pair.or.length.of.read.in.non.overlapping.pair','Count'), c('BAM','LEN','COUNT'))
				set(ans, NULL, 'BAM', ans[,gsub('\\.bam','',as.character(BAM))])
				ans
			}, by='F']
	set(tmp, NULL, 'F',NULL)
	tmp[, DATA:='BEEHIVE']
	db		<- rbind(db, tmp)
	save(db, file=outfile)
}

project.bam.read.distribution.evaluate<- function() 
{
	infile	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/bam_readdistr.rda'
	load(infile)
	#	quantiles
	qs		<- c(0.01,seq(0.1,0.9,0.1),0.99)
	ans		<- db[, list(	Q= as.numeric(quantile(rep(LEN, COUNT),prob=qs)), P= paste0(100*qs,'%')), by=c('DATA','BAM')]
	#	sort individuals by 50% quantile
	tmp		<- subset(ans, P=='50%')
	setkey(tmp, DATA, Q)
	tmp[, SORTID:= seq_len(nrow(tmp))]
	ans		<- merge(ans, subset(tmp, select=c(DATA,BAM,SORTID)), by=c('DATA','BAM'))	
	#	sort individuals by each quantile
	setkey(ans, DATA, P, Q)
	tmp		<- ans[, list(BAM=BAM, SORTID2=seq_along(BAM)), by=c('DATA','P')]
	ans		<- merge(ans, tmp, by=c('DATA','BAM','P'))
	setkey(ans, DATA, P, Q)
	#	make plot
	tmp		<- subset(ans, P%in%c('1%','10%','30%','50%','70%','90%','99%'))
	tmp		<- subset(ans, P%in%c('50%','90%','99%'))
	set(tmp, NULL, 'P', tmp[, factor(P, levels=c('50%','90%','99%'), labels=c('50%','10%','1%'))])
	ggplot(tmp, aes(y=SORTID2, x=Q, colour=P, group=interaction(DATA,P))) +
			geom_line() +
			theme_bw() + 
			labs(y='number of samples\n',x='\nlength of merged read fragments',colour='proportion\nof reads\nper sample') +
			facet_wrap(~DATA, scales='free_x', ncol=2)
	ggsave(file=gsub('.rda','_propreadsatleastx.pdf',outfile), w=8, h=6)
	
	
	#	read len cutoffs
	ps		<- seq(200, 350, 25)
	ans2	<- db[, list( P=ps,	LEN=ecdf(rep(LEN, COUNT))(ps)), by=c('DATA','BAM')]
	#set(ans2, NULL, 'P', ans2[, as.integer(gsub('%','',P))/100])
	tmp		<- ans2[, list(	STAT=c('ymin','lower','middle','upper','ymax'), 
							Q=quantile(LEN, p=c(0.025,0.25,0.5,0.75,0.975))
							), by=c('DATA','P')]
	tmp		<- dcast.data.table(tmp, DATA+P~STAT, value.var='Q')
	ggplot(tmp, aes(x=P, ymin=ymin, lower=lower, middle=middle, upper=upper, ymax=ymax, fill=DATA)) + 
			geom_boxplot(stat='identity') +
			theme_bw() + 
			scale_y_continuous(labels=scales:::percent, expand=c(0,0), lim=c(0,1)) +
			labs(x='min length of\nmerged read fragments\n',y='reads per sample discarded\n(%)\n', fill='study')
	ggsave(file=gsub('.rda','_propreadsdiscarded.pdf',outfile), w=5, h=5)
	tmp		<- ans2[, list( Q= length(which(LEN>0.9)) ), by=c('DATA','P')]	
	ggplot(tmp, aes(x=P, y=Q, fill=DATA)) + 
			geom_bar(stat='identity', position='dodge') +
			theme_bw() + 
			scale_y_continuous(expand=c(0,0), lim=c(0,2700)) +
			labs(x='min length of\nmerged read fragments\n',y='samples with >90% of reads discarded\n(#)\n', fill='study')
	ggsave(file=gsub('.rda','_samplesreadsdiscarded.pdf',outfile), w=5, h=5)			
}

project.bam.read.distribution.calculate<- function() 
{
	require(big.phylo)
	require(Phyloscanner.R.utilities)
	#
	#	INPUT ARGS PLATFORM
	#
	
	#
	#	TEST ARGS
	#		
	if(0)
	{
		#HOME		<<- "~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA"
		load( file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_UG_gag_coinfinput_160219.rda") )
		setnames(pty.runs, c('FILE_ID'), c('SAMPLE_ID'))
		pty.runs					<- subset(pty.runs, select=c(SAMPLE_ID, PTY_RUN))
		pty.data.dir				<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'
		work.dir					<- '/Users/Oliver/duke/2016_PANGEAphylotypes/Rakai_ptinput'
		out.dir						<- '/Users/Oliver/duke/2016_PANGEAphylotypes/test_ptoutput'	
		prog.bam.distr.calculator	<- '/Users/Oliver/git/phylotypes/tools/EstimateReadCountPerWindow.py'				
		hpc.load					<- ""		
	}
	#
	#	INPUT ARGS
	#
	if(1)
	{
		in.dir						<- file.path(HOME,"RakaiAll_input_170301")
		load( file.path(in.dir, 'Rakai_phyloscanner_170301_b75.rda') )
		tmp							<- subset(pty.runs, select=c(SID, PTY_RUN))
		load( file.path(in.dir, 'Rakai_phyloscanner_170704_b75_part2.rda') )
		pty.runs					<- subset(pty.runs, select=c(SID, PTY_RUN))
		pty.runs					<- rbind(tmp, pty.runs)
		setnames(pty.runs, c('SID'), c('SAMPLE_ID'))
		work.dir					<- file.path(HOME,"RakaiAll_work_170704")
		out.dir						<- file.path(HOME,"RakaiAll_output_170704_bamdistr")		
		prog.bam.distr.calculator	<- '/work/or105/libs/phylotypes/tools/EstimateReadCountPerWindow.py'
		pty.data.dir				<- '/work/or105/PANGEA_mapout/data'
		hpc.load					<- "module load intel-suite/2015.1 anaconda/2.3.0 samtools"										
	}
	#
	#	make scripts
	#
	pty.args			<- list( 	prog.bam.distr.calculator=prog.bam.distr.calculator,
									data.dir=pty.data.dir,
									out.dir=out.dir,
									work.dir=work.dir#,
									#quality.trim.ends=23, 
									#min.internal.quality=23 
									)					
	pty.c				<- phsc.cmd.bam.calculate.read.distribution(pty.runs, pty.args)
	pty.c				<- subset(pty.c, PTY_RUN%in%c(774,775,776,778,779,780,781,782,783,784,785,786,787,788,789,811,814,838,844,855,856,859,860,861,865,867,870,871,872,877,887,888,896,897,906,907,910,945,946,947,949,951,953,955,963,964,967,968,969,972,974,975,981,982,983,987,988,989,995,996,999,1004,1007,1009,1010,1013,1014,1015,1016,1018,1019,1021,1022,1023,1024,1025,1026,1027,1031,1032,1033,1035,1036,1037,1038,1039,1040,1042,1046,1047,1048,1049,1050,1051,1052,1053,1054,1055,1056,1570,1574,1575,1580,1582,1584,1585,1587,1589,1591,1592,1593,1595,1599,1600,1601,1602,1603,1604,1606,1609,1610,1611,1613,1614,1615,1616,1621,1624,1626,1627,1628,1629,1635,1636,1654,1655,1675,1686,1692,1695,1701,1702,1708,1709,1710,1721,1724,1726,1732,1733,1735,1736,1737,1738,1739,1740,1746,1747,1753,1755,1756,1757,1758,1759,1760,1761,1764,1765,1767,1769,1770,1773,1774,1775,1778,1779,1781,1782,1783,1785,1786,1787,1788,1789,1790,1791,1793,1794,1795,1796,1797,1798,1799,1800,1801,1802,1804,1805,1806,1807,1808,1809,1811,1812,1813,1814,1815,1816,1817,1818,1819,1820,1821,1822,1823,1824,1825,1826,1828,1829,1830,1831,1832,1833,1834,1835,1836,1837,1838,1839,1840,1841,1843,1844,1845,1846,1847,1848,1849,1850,1852,1853,1854,1855,1856,1857,1858,1859,1860,1861,1862,1863,1864,1865,1866,1867,1868,1869,1870,1871,1872,1873,1874,1875,1876,1877,1878,1879,1880))
	#pty.c[1,cat(CMD)]	
	#
	#	submit
	#
	invisible(pty.c[,	{
						cmd			<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=5, hpc.q="pqeelab", hpc.mem="5900mb",  hpc.nproc=1, hpc.load=hpc.load)
						#cmd		<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.walltime=71, hpc.q=NA, hpc.mem="63850mb",  hpc.nproc=hpc.nproc, hpc.load=hpc.load)						
						cmd			<- paste(cmd,CMD,sep='\n')
						cat(cmd)					
						outfile		<- paste("bd",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
						cmd.hpccaller(pty.args[['work.dir']], outfile, cmd)
						#stop()
					}, by='PTY_RUN'])
	quit('no')	
}

project.readlength.calculate.coverage.171208<- function()
{
	require(ggplot2)
	require(data.table)
	require(Rsamtools)		
	#indir			<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_bamdistr'		
	#outfile.base	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/bam_stats_171208'	
	#indir			<- '/work/or105/PANGEA_mapout/mergedfragmentlen_BEE'
	#outfile.base	<- '/work/or105/PANGEA_mapout/mergedfragmentlen_BEE_'
	indir			<- '~/Box Sync/OR_Work/2018/2018_MunichCluster/mergedfragmentlen'
	outfile.base	<- '~/Box Sync/OR_Work/2018/2018_MunichCluster/mergedfragmentlen_'
	
	
	bfiles			<- data.table(FILE=list.files(indir, pattern='_mergedfragmentlen.rda$', full.names=TRUE))		
	bfiles[, FILE_ID:=gsub('_mergedfragmentlen.rda','',basename(FILE))]	
	print(bfiles)
	
	bam.cov		<- bfiles[,	{
				#FILE<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_bamdistr/12559_1_1_mergedfragmentlen.rda'
				load(FILE)
				tmp	<- coverage(IRanges(start=dlen$POS, width=dlen$LEN))				
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths)	
			}, by='FILE_ID']
	save(bam.cov, file=paste0(outfile.base,'_covall.rda'))
	bam.cov		<- NULL
	gc()
	
	bam.cov175	<- bfiles[,{
				load(FILE)
				dlen<- subset(dlen, LEN>=175)
				tmp	<- coverage(IRanges(start=dlen$POS, width=dlen$LEN))				
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths)				
			}, by='FILE_ID']
	save(bam.cov175, file=paste0(outfile.base,'_cov175.rda'))
	bam.cov175	<- NULL
	gc()
	
	bam.cov200	<- bfiles[,{
				load(FILE)
				dlen<- subset(dlen, LEN>=200)
				tmp	<- coverage(IRanges(start=dlen$POS, width=dlen$LEN))				
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths)
			}, by='FILE_ID']
	save(bam.cov200, file=paste0(outfile.base,'_cov200.rda'))
	bam.cov200	<- NULL
	gc()
	
	bam.cov225	<- bfiles[,{
				load(FILE)
				dlen<- subset(dlen, LEN>=225)
				tmp	<- coverage(IRanges(start=dlen$POS, width=dlen$LEN))				
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths)
			}, by='FILE_ID']
	save(bam.cov225, file=paste0(outfile.base,'_cov225.rda'))
	bam.cov225	<- NULL
	gc()
	
	bam.cov250	<- bfiles[,{
				load(FILE)
				dlen<- subset(dlen, LEN>=250)
				tmp	<- coverage(IRanges(start=dlen$POS, width=dlen$LEN))				
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths)
			}, by='FILE_ID']
	save(bam.cov250, file=paste0(outfile.base,'_cov250.rda'))
	bam.cov250	<- NULL
	gc()
	
	bam.cov275	<- bfiles[,{
				load(FILE)
				dlen<- subset(dlen, LEN>=275)
				tmp	<- coverage(IRanges(start=dlen$POS, width=dlen$LEN))				
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths)
			}, by='FILE_ID']
	save(bam.cov275, file=paste0(outfile.base,'_cov275.rda'))
	bam.cov275	<- NULL
	gc()
	
	bam.cov300	<- bfiles[,{
				load(FILE)
				dlen<- subset(dlen, LEN>=300)
				tmp	<- coverage(IRanges(start=dlen$POS, width=dlen$LEN))				
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths)
			}, by='FILE_ID']
	save(bam.cov300, file=paste0(outfile.base,'_cov300.rda'))
	bam.cov300	<- NULL
	gc()
	
	#
	bam.cov325	<- bfiles[,{
				load(FILE)
				dlen<- subset(dlen, LEN>=325)
				tmp	<- coverage(IRanges(start=dlen$POS, width=dlen$LEN))				
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths)
			}, by='FILE_ID']
	save(bam.cov325, file=paste0(outfile.base,'_cov325.rda'))
	bam.cov325	<- NULL
	gc()
	
	bam.cov350	<- bfiles[,{
				load(FILE)
				dlen<- subset(dlen, LEN>=350)
				tmp	<- coverage(IRanges(start=dlen$POS, width=dlen$LEN))				
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths)
			}, by='FILE_ID']
	save(bam.cov350, file=paste0(outfile.base,'_cov350.rda'))
	bam.cov350	<- NULL
	gc()
	
	bam.cov375	<- bfiles[,{
				load(FILE)
				dlen<- subset(dlen, LEN>=375)
				tmp	<- coverage(IRanges(start=dlen$POS, width=dlen$LEN))				
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths)
			}, by='FILE_ID']
	save(bam.cov375, file=paste0(outfile.base,'_cov375.rda'))
	bam.cov375	<- NULL
	gc()
	
	bam.len		<- bfiles[,{
				#FILE<- "14939_1_80.bam"				
				load(FILE)				
				list(QU=seq(0,320,20), CDF=ecdf(dlen$LEN)(seq(0,320,20)))								
			}, by='FILE_ID']	
	save(bam.len, file=paste0(outfile.base,'_lendist.rda'))
	bam.len		<- NULL
	gc()
}
	
project.readlength.count.bam.171208<- function()
{
	require(ggplot2)
	require(data.table)
	require(Rsamtools)		
	#pty.data.dir	<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'	
	pty.data.dir	<- '/work/or105/PANGEA_mapout/data'
	out.dir			<- '/work/or105/PANGEA_mapout/mergedfragmentlen'
	pty.data.dir	<- '~/duke/2018_MunichCluster/Data'
	out.dir			<- '~/Box Sync/OR_Work/2018/2018_MunichCluster/mergedfragmentlen'
				
		
	bfiles			<- data.table(FILE=list.files(pty.data.dir, pattern='bam$'))	
	bfiles			<- subset(bfiles, !grepl('Contam',FILE))
	bfiles[, FILE_ID:=gsub('.bam','',FILE)]	
	tmp				<- data.table(FILE_OUT=list.files(out.dir, pattern='_mergedfragmentlen.rda', full.names=TRUE))
	tmp[, FILE_ID:= gsub('_mergedfragmentlen.rda','',basename(FILE_OUT))]
	bfiles	<- merge(bfiles, tmp, by='FILE_ID', all.x=1)
	print(bfiles)
	bfiles	<- subset(bfiles, is.na(FILE_OUT))
	setkey(bfiles, FILE_ID)		
	print(bfiles)
	
	invisible(bfiles[,{
				cat('\nreading',FILE)
				#FILE<- '15-01402.bam'
				dlen	<- phsc.bam.get.length.and.pos.of.mergedreads( file.path(pty.data.dir,FILE), error.strict=FALSE )
				dlen[, QNAME:=NULL]				
				save(dlen, file=file.path(out.dir,paste0(gsub('\\.bam','_mergedfragmentlen.rda',FILE))))
				NULL
			}, by='FILE_ID'])
}

project.readlength.count.bam.171018<- function()
{
	require(ggplot2)
	require(data.table)
	require(Rsamtools)		
	#pty.data.dir	<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'	
	pty.data.dir	<- '/work/or105/PANGEA_mapout/data'
	outfile			<- file.path(pty.data.dir,'bam_stats_171018.rda')			
	
	
	bfiles			<- data.table(FILE=list.files(pty.data.dir, pattern='bam$'))	
	bfiles			<- subset(bfiles, !grepl('Contam',FILE))
	bfiles[, FILE_ID:=gsub('.bam','',FILE)]	
	
	#FILE<- "14939_1_80.bam"
	#cat('\nreading',FILE)
	#bam<- scanBam(file.path(pty.data.dir,FILE))[[1]]	#scan entire file
	bam.cov		<- bfiles[,{
				cat('\nCOV reading',FILE)
				z	<- scanBam(file.path(pty.data.dir,FILE), param=ScanBamParam(what=c('qwidth','pos','rname')))[[1]]
				tmp	<- IRanges(start=z$pos, width=z$qwidth)
				tmp <- coverage(tmp)
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths, REF=z$rname[1])
			}, by='FILE_ID']	
	bam.cov175	<- bfiles[,{
				cat('\nCOV reading',FILE)
				z	<- scanBam(file.path(pty.data.dir,FILE), param=ScanBamParam(what=c('qwidth','pos','rname')))[[1]]
				tmp	<- z[['qwidth']]>=175
				z	<- lapply(z, function(x) x[tmp])					
				tmp	<- IRanges(start=z$pos, width=z$qwidth)
				tmp <- coverage(tmp)
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths, REF=z$rname[1])
			}, by='FILE_ID']
	bam.cov200	<- bfiles[,{
				cat('\nCOV reading',FILE)
				z	<- scanBam(file.path(pty.data.dir,FILE), param=ScanBamParam(what=c('qwidth','pos','rname')))[[1]]
				tmp	<- z[['qwidth']]>=200
				z	<- lapply(z, function(x) x[tmp])					
				tmp	<- IRanges(start=z$pos, width=z$qwidth)
				tmp <- coverage(tmp)
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths, REF=z$rname[1])
			}, by='FILE_ID']	
	bam.cov250	<- bfiles[,{
				cat('\nCOV reading',FILE)
				z	<- scanBam(file.path(pty.data.dir,FILE), param=ScanBamParam(what=c('qwidth','pos','rname')))[[1]]
				tmp	<- z[['qwidth']]>=250
				z	<- lapply(z, function(x) x[tmp])					
				tmp	<- IRanges(start=z$pos, width=z$qwidth)
				tmp <- coverage(tmp)
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths, REF=z$rname[1])
			}, by='FILE_ID']
	bam.cov280	<- bfiles[,{
				cat('\nCOV reading',FILE)
				z	<- scanBam(file.path(pty.data.dir,FILE), param=ScanBamParam(what=c('qwidth','pos','rname')))[[1]]
				tmp	<- z[['qwidth']]>=280
				z	<- lapply(z, function(x) x[tmp])					
				tmp	<- IRanges(start=z$pos, width=z$qwidth)
				tmp <- coverage(tmp)
				list(POS= cumsum(c(1L,tmp@lengths[-length(tmp@lengths)])), COV=tmp@values, REP=tmp@lengths, REF=z$rname[1])
			}, by='FILE_ID']
	#	get lengths of all reads in quality trimmed bam file
	#z	<- scanBam(file.path(pty.data.dir,FILE), param=ScanBamParam(what=c('qwidth','qual')))
	bam.len			<- bfiles[,{
				#FILE<- "14939_1_80.bam"
				cat('\nLEN reading',FILE)
				z	<- scanBam(file.path(pty.data.dir,FILE), param=ScanBamParam(what=c('qwidth')))
				list(QU=seq(0,320,20), CDF=ecdf(z[[1]][['qwidth']])(seq(0,320,20)))
				#list(PR=seq(0.01,0.99,0.01), QU=quantile(z[[1]][['qwidth']], p=seq(0.01,0.99,0.01)))				
			}, by='FILE']	
	#
	save(bam.len, bam.cov, bam.cov175, bam.cov200, bam.cov250, bam.cov280, file= outfile)
}

pty.pipeline.runs.coverage	<- function()
{
	infile.bam		<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/readlengths/bam_stats_150218.rda'
	load(infile.bam)
	infile.runs		<- file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_UG_NoQ_fast2_coinfrunsinput.rda")
	infile.runs		<- "/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_UG_gag_selecthelp.rda"
	infile.runs		<- "/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_UG_gag_coinfinput_160219.rda"
	infile.runs		<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda"
	load( infile.runs )
	
	setkey(pty.runs, PTY_RUN, FILE_ID)
	pty.runs		<- unique(pty.runs)
	bam.cov			<- merge(bam.cov, pty.runs, by='FILE_ID')
	setkey(bam.cov, PTY_RUN, FILE_ID, POS)
	
	ps				<- lapply(bam.cov[, unique(PTY_RUN)], function(ptyr)
			{
				cat('\nprocess run', ptyr)
				tmp				<- subset(bam.cov, PTY_RUN==ptyr)
				tmp				<- tmp[, {
							z	<- rep(COV,REP)
							list(COV=z, POS=seq_along(z), REF=REF[1])
						}, by=c('FILE_ID','PTY_RUN')]
				p	<- ggplot(tmp, aes(x=POS, y=COV, colour=FILE_ID)) + geom_step() + theme_bw() +
						scale_y_log10() +
						facet_grid(~PTY_RUN) +
						labs(x='genome position\n(relative to reference)', y='read coverage', colour='patient')
				p
			})
	#pdf(file=gsub('\\.rda','_coverage_UG.pdf',infile), w=12, h=5)
	pdf(file=gsub('\\.rda','_vertical_read_coverage.pdf',infile.runs), w=12, h=5)
	for(i in seq_along(ps))
		print(ps[[i]])	
	dev.off()		
	
	#bam.min	<- subset(bam.cov, POS>10 & POS<8000)[, list(MC= min(COV)), by='FILE_ID']
	#setkey(bam.min, MC)
}

project.readlength.ZA.count.bam.150120<- function()
{
	require(ggplot2)
	require(data.table)
	require(Rsamtools)
	pty.infile		<- file.path(HOME,"data","PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_ptyrunsinput.rda")
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
	#
	bam.len[, SEQ_LOC:='AfricaCentre']
	set(bam.len, bam.len[, which(!grepl('^R[0-9]+_',FILE))],'SEQ_LOC','Sanger')
	bam.len[, SEQ_RUN:=0L]
	tmp			<- bam.len[, which(SEQ_LOC=='AfricaCentre')]
	set( bam.len, tmp, 'SEQ_RUN', bam.len[tmp, as.integer(gsub('R','',regmatches(FILE,regexpr('^R[0-9]+',FILE))))] )
	set(bam.len,NULL,'SEQ_RUN', bam.len[, paste(SEQ_LOC,SEQ_RUN,sep='-')])
	#
	save(bam.len, file= gsub('ptyrunsinput','bamlen',pty.infile))
	#
	ggplot(subset(bam.len,QU>=40 & QU<320), aes(y=CDF, x=factor(QU), fill=SEQ_RUN)) + geom_boxplot() + 
			theme_bw() + labs(y='cumulative frequency\nin one individual\n', x='\nlength of quality-trimmed short reads\n(nt)', fill='sequence run') +
			theme(legend.position='bottom')
	ggsave(file= gsub('ptyrunsinput.rda','bamlen_byrun.pdf',pty.infile), w=14, h=7)
	#
	ggplot(bam.len, aes(y=CDF, x=factor(QU), fill=TYPE)) + geom_boxplot() + 
			theme_bw() + labs(y='cumulative frequency\nin one individual\n', x='\nlength of quality-trimmed short reads\n(nt)', fill='individual') +
			theme(legend.position='bottom')
	ggsave(file= gsub('ptyrunsinput.rda','bamlen.pdf',pty.infile), w=10, h=7)
	#
	#	compare against horizontal coverage
	#
	sqi[, FILE_ID:= TAXA]
	tmp		<- sqi[, which(!is.na(SANGER_ID))]
	set(sqi, tmp, 'FILE_ID', sqi[tmp, SANGER_ID])
	bam.len	<- merge(bam.len, subset(sqi, select=c(FILE_ID, COV)), by='FILE_ID', all.x=1)
	tmp		<- subset(bam.len, QU==140)
	ggplot(tmp, aes(y=100*CDF, x=COV, colour=SEQ_RUN)) + 
			geom_point() + 
			labs(x='horizontal coverage of consensus\n(bp)', y='proportion of reads with len<=140bp\n(%)', colour='sequence run') +
			theme_bw()
	ggsave(file= gsub('ptyrunsinput.rda','bamlen_vs_horizontalcoverage.pdf',pty.infile), w=7, h=7)
	write.csv(tmp, file=gsub('ptyrunsinput.rda','bamlen_vs_horizontalcoverage.csv',pty.infile),  row.names=FALSE)
}

project.dualinfecions.phylotypes.evaluatereads.150119<- function()
{	
	HOME			<<- "~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA"			
	indir			<- file.path(HOME,"phylotypes_160119")		
	pty.evaluate.fasta(indir, strip.max.len=350, select='^ptyr1')
	
	#pty.infile		<- file.path(HOME,"data", "PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_ptyrunsinput.rda")
	if(0)
	{		
		indir			<- file.path(HOME,"phylotypes_160120")
		infiles	<- list.files(indir, pattern='rda$')
		seqd	<- do.call('rbind',lapply(seq_along(infiles),function(i){
					load(file.path(indir,infiles[i]))
					tmp		<- seqd[, list(LEN_MED=as.numeric(median(LEN)), LEN_QL=as.numeric(quantile(LEN,p=0.025)), LEN_QU=as.numeric(quantile(LEN,p=0.975))), by=c('PTY_RUN','IND','W_FROM')]
					tmp[, LEN_MED_RM:=tmp[, rollapply(LEN_MED, width=10, FUN=mean, align="center", partial=TRUE)]]
					setkey(seqd, PTY_RUN, W_FROM, IND)
					seqd	<- unique(seqd)
					seqd	<- merge(seqd, tmp, by=c('PTY_RUN','IND','W_FROM'))
					seqd
				}))
		seqd[, IND_L:= paste('(',substring(FILL,1,1),') ',IND,sep='')]
		#
		#	number of unique reads per window and individual		
		tmp		<- subset(seqd,FILL=='candidate')[, list(CAND_N=length(IND)), by=c('PTY_RUN','W_FROM')]
		ggplot(tmp, aes(x=W_FROM, y=CAND_N)) + geom_step() +
				scale_x_continuous(breaks=seq(0,10000,500)) +			
				labs(x='window start', y='candidate individuals\n(#)') +
				facet_grid(PTY_RUN~., scales='free_y') + theme_bw()
		ggsave(file=file.path(indir,paste('pty_candidate_individuals.pdf')), w=10, h=80,limitsize = FALSE)
		#
		ggplot(seqd,aes(x=W_FROM, y=UNIQUE_N, fill=IND_L, colour=FILL, alpha=FILL)) + geom_bar(stat='identity') + theme_bw() +
				scale_x_continuous(breaks=seq(0,10000,500)) +
				scale_fill_discrete(guide=FALSE) +
				scale_colour_manual(values=c('candidate'='black', 'filler'='transparent'), guide=FALSE) +
				scale_alpha_manual(values=c('candidate'=1, 'filler'=0.2))+
				labs(x='window start', y='unique reads\n(#)', fill='individuals') +
				facet_grid(PTY_RUN~., scales='free_y') + theme_bw()
		ggsave(file=file.path(indir,'pty_unique_reads.pdf'), w=10, h=80,limitsize = FALSE)
		#
		tmp		<- subset(seqd, FILL=='candidate')[, list(TOTAL_N=sum(UNIQUE_N)), by=c('PTY_RUN','IND')]
		ggplot(tmp, aes(y=IND, yend=IND, x=0, xend=TOTAL_N, colour=factor(TOTAL_N<1e2, levels=c(TRUE,FALSE),labels=c('<100','>=100')))) + geom_segment(size=3) +
				labs(y='candidate individuals', x='sum of unique reads\n(#)', colour='') + theme_bw() +
				scale_colour_brewer(palette='Set1') + facet_grid(PTY_RUN~., scales='free_y', space='free_y') +
				theme(legend.position='bottom')
		ggsave(file=file.path(indir,'pty_reads_from_candidate.pdf'), w=5, h=30, limitsize = FALSE)
		#		
		ggplot(seqd, aes(x=W_FROM, y=LEN_MED_RM, alpha=FILL, colour=IND_L, group=IND_L)) + 
				geom_line() + theme_bw() + scale_x_continuous(breaks=seq(0,10000,500)) +
				scale_alpha_manual(values=c('candidate'=1, 'filler'=0.2))+
				labs(x='window start', y='median read length\n(rolling mean over 10 windows)') + 
				facet_grid(PTY_RUN~., scales='free_y') + guides(colour=FALSE) 
		ggsave(file=file.path(indir,'pty_read_lengths.pdf'), w=10, h=60, limitsize = FALSE)
		#	how many windows have more than 10 unique reads for all candidate individuals?
		tmp		<- seqd[, list(SELECT=UNIQUE_N[FILL=='candidate']), by=c('PTY_RUN','W_FROM')]
		tmp		<- tmp[, list(SELECT=min(SELECT)), by=c('PTY_RUN','W_FROM')]		
		ggplot(tmp, aes(x=W_FROM, y=SELECT, group=PTY_RUN)) + geom_line() + scale_y_log10() + 
			facet_grid(PTY_RUN~.) + theme_bw() + labs(x='window start', y='min unique reads among all candidates')
		ggsave(file=file.path(indir,'pty_minreads_from_candidate.pdf'), w=5, h=30, limitsize = FALSE)	
		tmp		<- subset(tmp, SELECT>10)[, list(SELECT_N=length(SELECT)), by='PTY_RUN']
		ggplot(tmp, aes(y=PTY_RUN, yend=PTY_RUN, x=0, xend=SELECT_N)) + geom_segment(size=2) +
				scale_y_continuous(breaks=seq.int(1,300)) + theme_bw() +
				labs(y='run', x='number of selected windows')
		ggsave(file=file.path(indir,'pty_selected_windows.pdf'), w=5, h=20, limitsize = FALSE)
		
		
	}
	
}

