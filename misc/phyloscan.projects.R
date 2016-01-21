project.dual<- function()
{
	HOME		<<- '/work/or105/Gates_2014/2015_PANGEA_DualPairsFromFastQIVA'
	#HOME		<<- "~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA"	
	#project.dual.distances.231015()
	#project.dual.examl.231015()
	pty.pipeline.fasta()
	#project.dualinfecions.phylotypes.pipeline.examl.160110()
	#project.dualinfecions.phylotypes.evaluatereads.150119()
	
	#	various
	if(0)
	{
		cmd			<- hivc.cmd.various()
		cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeph', hpc.walltime=3, hpc.mem="3600mb")
		cat(cmd)		
		outdir		<- file.path(HOME,"ptyruns")
		outfile		<- paste("pv",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		hivc.cmd.hpccaller(outdir, outfile, cmd)
		quit("no")	
	}			
}


project.dual.alignments.160110<- function()
{
	outdir	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110'
	#	read info
	#file	<- "~/Dropbox (Infectious Disease)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	#si		<- as.data.table(read.csv(file, stringsAsFactors=FALSE))
	#setnames(si, colnames(si), toupper(gsub('.','_',colnames(si),fixed=1))) 
	#set(si, NULL, 'PANGEA_ID', si[, gsub(' ','',PANGEA_ID)])
	#setnames(si, 'CLINICAL_GENOME_COVERAGE', 'COV')
	
	#	read global PANGEA alignment w SA seqs and split by site	
	file			<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/AfricaCentreSeqs/GlobalAln_PlusAllPANGEA.fasta'
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

project.dual.alignments.151023<- function()
{
	outdir	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_151113'
	#	read info
	file	<- "~/Dropbox (Infectious Disease)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	si		<- as.data.table(read.csv(file, stringsAsFactors=FALSE))
	setnames(si, colnames(si), toupper(gsub('.','_',colnames(si),fixed=1))) 
	set(si, NULL, 'PANGEA_ID', si[, gsub(' ','',PANGEA_ID)])
	setnames(si, 'CLINICAL_GENOME_COVERAGE', 'COV')
	
	#	read global PANGEA alignment and split by site	
	file			<- "~/Dropbox (Infectious Disease)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_GlobalAlignment.fasta"
	file			<- '~/Dropbox (Infectious Disease)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.fasta'
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
	file	<- "~/Dropbox (Infectious Disease)/pangea_data/PANGEAcontigs_2015-09_Imperial/contigs_cnsalign_PNGIDn3366_CNTGSn6120_stripped99.fasta"
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
	
	infile		<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/PANGEA_HIV_n5003_Imperial_v160110_GlobalAlignment.rda'
	load(infile)	#loads sqi, sq
	indir		<- "~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data"
	outdir		<- indir
	infile		<-  "PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500.rda"
	load(paste(indir,infile,sep='/'))	#loads "ph" "dist.brl" "ph.gdtr"  "ph.mrca"
	#	add SANGER_ID
	infile.s	<- "~/Dropbox (Infectious Disease)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	si			<- as.data.table(read.csv(infile.s, stringsAsFactors=FALSE))
	setnames(si, colnames(si), toupper(gsub('.','_',colnames(si),fixed=1))) 
	set(si, NULL, 'PANGEA_ID', si[, gsub(' ','',PANGEA_ID)])
	setnames(si, 'CLINICAL_GENOME_COVERAGE', 'COV')
	tmp			<- subset(si, grepl('^PG14-ZA', PANGEA_ID), c(PANGEA_ID, SANGER_ID))
	set(tmp, NULL, 'PANGEA_ID', tmp[,gsub('-','_',PANGEA_ID)])
	sqi			<- merge(sqi, tmp, by='PANGEA_ID', all.x=1)	
	#	delete duplicates identified by Tulio from tree
	infile.dup	<- '/Users/Oliver/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/PANGEA_HIV_n5003_Imperial_v160110_ZA_Duplicates.csv'
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
		
	infile		<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/PANGEA_HIV_n5003_Imperial_v160110_GlobalAlignment.rda'
	load(infile)	#loads sqi, sq
	indir		<- "~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data"
	outdir		<- indir
	infile		<-  "PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500.rda"
	load(paste(indir,infile,sep='/'))	#loads "ph" "dist.brl" "ph.gdtr"  "ph.mrca"
	#	add SANGER_ID
	infile.s	<- "~/Dropbox (Infectious Disease)/pangea_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	si			<- as.data.table(read.csv(infile.s, stringsAsFactors=FALSE))
	setnames(si, colnames(si), toupper(gsub('.','_',colnames(si),fixed=1))) 
	set(si, NULL, 'PANGEA_ID', si[, gsub(' ','',PANGEA_ID)])
	setnames(si, 'CLINICAL_GENOME_COVERAGE', 'COV')
	tmp			<- subset(si, grepl('^PG14-ZA', PANGEA_ID), c(PANGEA_ID, SANGER_ID))
	set(tmp, NULL, 'PANGEA_ID', tmp[,gsub('-','_',PANGEA_ID)])
	sqi			<- merge(sqi, tmp, by='PANGEA_ID', all.x=1)	
	#	delete duplicates identified by Tulio from tree
	infile.dup	<- '/Users/Oliver/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/alignments_160110/PANGEA_HIV_n5003_Imperial_v160110_ZA_Duplicates.csv'
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
	#	
	cat('\nNumber of clusters=', pty.runs[, length(unique(CLU_ID))])
	cat('\nNumber of scheduled phylotype runs=', pty.runs[, max(PTY_RUN)])
	cat('\nNumber of selected taxa=', subset(pty.runs, !FILL)[, length(unique(TAXA))])	
	tmp			<- paste(indir, '/', gsub('\\.rda','_coinfrunsinput\\.rda',infile), sep='')
	save(pty.runs, pty.clu, ph, dist.brl, ph.gdtr, ph.mrca, clustering, sqi, sq, file= tmp)	
}


project.dualinfecions.phylotypes.test<- function()
{
	tfd	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/phylotypes'
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

	load("/Users/Oliver/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_ptyrunsinput.rda")
}

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

project.dualinfecions.phylotypes.evaluatereads.150119<- function()
{
	
	#HOME		<<- "~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA"
	pty.infile		<- file.path(HOME,"data", "PANGEA_HIV_n5003_Imperial_v160110_ZA_examlbs500_ptyrunsinput.rda")		
	indir			<- file.path(HOME,"phylotypes_160119")	
	pty.evaluate.fasta(pty.infile, indir)
	
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
		
		tmp		<- seqd[, list(SELECT=all(UNIQUE_N[FILL=='candidate']>10)), by=c('PTY_RUN','W_FROM')]
		subset(tmp, SELECT)[,]
		
		seqd[, table(UNIQUE_N>10, PTY_RUN)]
	}
	
}

