RakaiCouples.setup.phyloscan.runs<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ape)
	
	wdir	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples"
	#
	#	load table of sanger ids and pangea ids
	load("~/Dropbox (SPH Imperial College)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda")
	#	we need from load: dm,  sqi
	#
	#	load couples data
	infile	<- '~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/Pangea_Couples.csv'
	rc		<- as.data.table(read.csv(infile, stringsAsFactors=FALSE))
	setnames(rc, c('female.PANGEA.ID','male.PANGEA.ID'), c('female.TAXA','male.TAXA'))
	setnames(rc, colnames(rc), gsub('\\.','_',toupper(colnames(rc))))
	setnames(rc, c('MALE_RCCS_STUDYID','FEMALE_RCCS_STUDYID'), c('MALE_RID','FEMALE_RID'))
	set(rc, NULL, 'MALE_DATE', rc[, hivc.db.Date2numeric(as.Date(MALE_DATE))])
	set(rc, NULL, 'FEMALE_DATE', rc[, hivc.db.Date2numeric(as.Date(FEMALE_DATE))])
	set(rc, NULL, 'MALE_LASTNEGDATE', rc[, hivc.db.Date2numeric(as.Date(MALE_LASTNEGDATE))])
	set(rc, NULL, 'FEMALE_LASTNEGDATE', rc[, hivc.db.Date2numeric(as.Date(FEMALE_LASTNEGDATE))])	
	set(rc, NULL, 'MALE_FIRSTPOSDATE', rc[, hivc.db.Date2numeric(as.Date(MALE_FIRSTPOSDATE))])
	set(rc, NULL, 'FEMALE_FIRSTPOSDATE', rc[, hivc.db.Date2numeric(as.Date(FEMALE_FIRSTPOSDATE))])
	#
	#	see if couples in same phylotype runs that were already set up
	infile	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_UG_gag_coinfinput_160219.rda'
	load(infile)
	#	prepare patristic distance matrix
	ph.gdtr	<- as.data.table(melt(ph.gdtr, varnames=c('TAXA1','TAXA2')))
	setnames(ph.gdtr, 'value', 'PD')
	ph.gdtr	<- subset(ph.gdtr, TAXA1!=TAXA2)
	set(ph.gdtr, NULL, 'TAXA1', ph.gdtr[, gsub('_','-',as.character(TAXA1))])
	set(ph.gdtr, NULL, 'TAXA2', ph.gdtr[, gsub('_','-',as.character(TAXA2))])
	tmp		<- subset(pty.runs, select=c(PTY_RUN, TAXA))
	setnames(tmp, c('TAXA','PTY_RUN'), c('TAXA2','PTY_RUN'))
	ph.gdtr	<- merge(ph.gdtr, tmp, all.x=1, by='TAXA2')
	#	load genetic distance matrix with overlap
	infile		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/circumcision/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_gd.rda'
	load(infile)	#loads sq.gd
	#	reduce genetic distances to Rakai seqs
	tmp			<- subset(dm, COHORT=='RCCS', TAXA) 
	setnames(tmp, 'TAXA', 'TAXA1')
	sq.gd		<- merge(sq.gd, tmp, by='TAXA1')
	setnames(tmp, 'TAXA1', 'TAXA2')
	sq.gd		<- merge(sq.gd, tmp, by='TAXA2')
	#	add PTY_RUN
	tmp			<- subset(pty.runs, select=c(PTY_RUN, TAXA))
	setnames(tmp, c('TAXA','PTY_RUN'), c('TAXA2','PTY_RUN'))
	set(tmp, NULL, 'TAXA2', tmp[, gsub('_','-',as.character(TAXA2))])
	sq.gd		<- merge(sq.gd, tmp, all.x=1, by='TAXA2')
	#
	#
	#
	tmp		<- copy(pty.runs)
	setnames(tmp, c('FILE_ID','PTY_RUN','IDX'), c('FEMALE_SANGER_ID','FEMALE_PTY_RUN','FEMALE_PHIDX'))
	tmp		<- subset(tmp, select=c('FEMALE_SANGER_ID','FEMALE_PTY_RUN','FEMALE_PHIDX'))
	setkey(tmp, FEMALE_SANGER_ID)
	rc		<- merge(rc, unique(tmp), by='FEMALE_SANGER_ID', all.x=1)
	setnames(tmp, colnames(tmp), gsub('^FE','',colnames(tmp)))
	rc		<- merge(rc, unique(tmp), by='MALE_SANGER_ID', all.x=1)
	#	get all possible unique pairings of SANGER_IDs from couples
	rp		<- subset(rc, !is.na(MALE_SANGER_ID), c(COUPID, MALE_RID, MALE_SANGER_ID, MALE_TAXA, MALE_PTY_RUN ))
	tmp		<- subset(rc, !is.na(FEMALE_SANGER_ID), c(COUPID, FEMALE_RID, FEMALE_SANGER_ID, FEMALE_TAXA, FEMALE_PTY_RUN ))	
	rp		<- merge(rp, tmp, by='COUPID')
	setkey(rp, FEMALE_SANGER_ID, MALE_SANGER_ID)
	rp		<- unique(rp)
	#	exclude couple F108560:F108560
	rp		<- subset(rp, FEMALE_SANGER_ID!=MALE_SANGER_ID)
	#
	cat('\npairings that were shipped', nrow(rp))
	#
	#	drop pairings including SANGER_IDs for which I don t have bam files 
	rp		<- merge(rp, data.table(MALE_TAXA=sqi[, TAXA]), by='MALE_TAXA')	
	rp		<- merge(rp, data.table(FEMALE_TAXA=sqi[, TAXA]), by='FEMALE_TAXA')
	#
	cat('\npairings for which I have bam files right now', nrow(rp))
	#
	#	assing run id to partners with no run id
	tmp		<- rp[, which(!is.na(MALE_PTY_RUN) & is.na(FEMALE_PTY_RUN))]
	set(rp, tmp, 'FEMALE_PTY_RUN', rp[tmp, MALE_PTY_RUN])
	tmp		<- rp[, which(is.na(MALE_PTY_RUN) & !is.na(FEMALE_PTY_RUN))]
	set(rp, tmp, 'MALE_PTY_RUN', rp[tmp, FEMALE_PTY_RUN])
	#	fixup inconsistent phylotype run COMPLEX_TOPOLOGYs
	tmp		<- subset(rp, FEMALE_PTY_RUN!=MALE_PTY_RUN)
	set(tmp, NULL, 'FEMALE_PTY_RUN', tmp[, MALE_PTY_RUN] )
	tmp2	<- rp[, which(FEMALE_PTY_RUN!=MALE_PTY_RUN)]
	set(rp, tmp2, 'MALE_PTY_RUN', rp[tmp2, FEMALE_PTY_RUN])
	rp		<- rbind(rp, tmp)
	#	assign run id to pairs with no run id; discard pairs with no reasonable overlap to anything else
	tmp		<- subset(rp, is.na(MALE_PTY_RUN))
	tmp2	<- subset(sq.gd, is.finite(PD))
	tmp		<- tmp[, {
				z	<- subset(tmp2, (TAXA1==MALE_TAXA|TAXA1==FEMALE_TAXA) & OVERLAP>400 & !is.na(PTY_RUN))
				if(nrow(z)>0)
					ans	<- z[which.min(PD), PTY_RUN]				
				if(nrow(z)==0 & nrow(subset(tmp2, (TAXA1==MALE_TAXA|TAXA1==FEMALE_TAXA) & OVERLAP>400)))
					ans	<- 0L
				if(nrow(z)==0)
					ans	<- -1L
				list(PTY_RUN= ans)
			}, by=c('COUPID','MALE_TAXA','FEMALE_TAXA')]
	rp		<- merge(rp, subset(tmp, PTY_RUN>0), by=c('COUPID','MALE_TAXA','FEMALE_TAXA'), all.x=1)
	tmp		<- rp[, which(is.na(MALE_PTY_RUN) & !is.na(PTY_RUN))]
	set(rp, tmp, c('FEMALE_PTY_RUN','MALE_PTY_RUN'), rp[tmp,PTY_RUN])
	set(rp, NULL, 'PTY_RUN', NULL)
	cat('\nnot enough sequence overlap to find good taxa (ie couple cons sequence too short)',subset(rp, is.na(MALE_PTY_RUN))[, paste(COUPID, collapse=', ')])
	rp		<- subset(rp, !is.na(MALE_PTY_RUN))
	setkey(rp, COUPID, FEMALE_TAXA, MALE_TAXA, MALE_PTY_RUN)
	rp		<- unique(rp)	
	#	make sure that multiple sequences from same male are never in same phylotype run
	tmp		<- rp[, list(N=length(MALE_SANGER_ID)), by=c('MALE_RID','MALE_PTY_RUN')]
	tmp		<- subset(tmp, N>1)
	tmp		<- merge(tmp, rp, by=c('MALE_RID','MALE_PTY_RUN'))
	tmp		<- tmp[, {
				#MALE_TAXA<- 'PG14-UG501971-S03975'; MALE_PTY_RUN<- c(64)
				cat('\n',MALE_RID)
				z	<- subset(tmp2, TAXA1%in%MALE_TAXA & OVERLAP>400 & !is.na(PTY_RUN) & !PTY_RUN%in%unique(MALE_PTY_RUN))
				if(nrow(z)==0)
					z	<- subset(tmp2, TAXA1%in%MALE_TAXA & OVERLAP>200 & !is.na(PTY_RUN) & !PTY_RUN%in%unique(MALE_PTY_RUN))
				if(nrow(z)==0)
					z	<- subset(tmp2, TAXA1%in%MALE_TAXA & OVERLAP>150 & !is.na(PTY_RUN) & !PTY_RUN%in%unique(MALE_PTY_RUN))				
				stopifnot(nrow(z)>0)
				z	<- z[, list(PD=min(PD)), by='PTY_RUN']
				setkey(z, PD)
				ans						<- MALE_PTY_RUN
				ans[duplicated(ans)]	<- z[seq_len( length(MALE_PTY_RUN)-length(unique(ans)) ), PTY_RUN]
				list(MALE_TAXA=MALE_TAXA, FEMALE_TAXA=FEMALE_TAXA, MALE_PTY_RUN=MALE_PTY_RUN, PTY_RUN= ans)
			}, by=c('MALE_RID')]
	rp		<- merge(rp, tmp, by=c('MALE_RID','MALE_TAXA','FEMALE_TAXA','MALE_PTY_RUN'), all.x=1)
	tmp		<- rp[, which(!is.na(PTY_RUN))]
	set(rp, tmp, c('FEMALE_PTY_RUN','MALE_PTY_RUN'), rp[tmp,PTY_RUN])
	set(rp, NULL, 'PTY_RUN', NULL)
	#	make sure that multiple sequences from same female are never in same phylotype run
	tmp		<- rp[, list(N=length(FEMALE_SANGER_ID)), by=c('FEMALE_RID','FEMALE_PTY_RUN')]
	tmp		<- subset(tmp, N>1)
	tmp		<- merge(tmp, rp, by=c('FEMALE_RID','FEMALE_PTY_RUN'))
	tmp		<- tmp[, {
				#FEMALE_TAXA<- 'PG14-UG501488-S03492'; FEMALE_PTY_RUN<- c(21)
				cat('\n',FEMALE_RID)
				z	<- subset(tmp2, TAXA1%in%FEMALE_TAXA & OVERLAP>400 & !is.na(PTY_RUN) & !PTY_RUN%in%unique(FEMALE_PTY_RUN))
				if(nrow(z)==0)
					z	<- subset(tmp2, TAXA1%in%FEMALE_TAXA & OVERLAP>200 & !is.na(PTY_RUN) & !PTY_RUN%in%unique(FEMALE_PTY_RUN))
				if(nrow(z)==0)
					z	<- subset(tmp2, TAXA1%in%FEMALE_TAXA & OVERLAP>149 & !is.na(PTY_RUN) & !PTY_RUN%in%unique(FEMALE_PTY_RUN))				
				if(nrow(z)==0)
					z	<- subset(tmp2, TAXA1%in%FEMALE_TAXA & !is.na(PTY_RUN) & !PTY_RUN%in%unique(FEMALE_PTY_RUN))
				#	exclude also run ids of the male partner
				ans		<- unique(rp[['MALE_PTY_RUN']][ which(rp[['MALE_RID']]%in%MALE_RID) ])
				z		<- subset(z, !PTY_RUN%in%ans)
				#	find min distance for each run
				z	<- z[, list(PD=min(PD)), by='PTY_RUN']
				setkey(z, PD)
				ans						<- FEMALE_PTY_RUN
				ans[duplicated(ans)]	<- z[seq_len( length(FEMALE_PTY_RUN)-length(unique(ans)) ), PTY_RUN]
				list(MALE_TAXA=MALE_TAXA, FEMALE_TAXA=FEMALE_TAXA, FEMALE_PTY_RUN=FEMALE_PTY_RUN, PTY_RUN= ans)
			}, by=c('FEMALE_RID')]
	rp		<- merge(rp, tmp, by=c('FEMALE_RID','MALE_TAXA','FEMALE_TAXA','FEMALE_PTY_RUN'), all.x=1)
	tmp		<- rp[, which(!is.na(PTY_RUN))]
	set(rp, tmp, c('FEMALE_PTY_RUN','MALE_PTY_RUN'), rp[tmp,PTY_RUN])
	set(rp, NULL, 'PTY_RUN', NULL)
	#	check
	stopifnot( nrow(subset(rp[, list(N=length(FEMALE_SANGER_ID)), by=c('FEMALE_RID','FEMALE_PTY_RUN')], N>1))==0 )
	stopifnot( nrow(subset(rp[, list(N=length(MALE_SANGER_ID)), by=c('MALE_RID','MALE_PTY_RUN')], N>1 ))==0 )
	#
	#	done with pair assignments
	#	
	#	condense pairs to list
	tmp		<- subset(rp, select=c(COUPID, MALE_TAXA, MALE_SANGER_ID, MALE_PTY_RUN))
	tmp[, PARTNER:='Male']
	setnames(tmp, c('MALE_TAXA','MALE_SANGER_ID','MALE_PTY_RUN'),c('TAXA','FILE_ID','PTY_RUN'))
	rp		<- subset(rp, select=c(COUPID, FEMALE_TAXA, FEMALE_SANGER_ID, FEMALE_PTY_RUN))
	rp[, PARTNER:='Female']
	setnames(rp, c('FEMALE_TAXA','FEMALE_SANGER_ID','FEMALE_PTY_RUN'),c('TAXA','FILE_ID','PTY_RUN'))
	rp		<- rbind(rp, tmp)
	setkey(rp, COUPID, FILE_ID, PTY_RUN)	
	# 	complete each run with individuals from outside couple so we have 22 individuals in each run
	tmp2	<- subset(sq.gd, is.finite(PD))
	ans		<- rp[, {
				cat('\n',PTY_RUN)
				#	select study participants not in couple 
				tmp	<- unlist(strsplit(COUPID,':',fixed=1))
				tmp	<- subset(dm, COHORT=='RCCS' & !STUDY_ID%in%tmp, TAXA)
				setnames(tmp, 'TAXA', 'TAXA2')
				tmp	<- merge(tmp2, tmp, by='TAXA2')				
				#	get those that are close to individuals in the current run
				tmp	<- subset(tmp, TAXA1%in%TAXA)
				setkey(tmp, PD)
				z	<- subset(tmp, OVERLAP>400)
				if(nrow(z)==0)
					z	<- subset(tmp, OVERLAP>200)
				if(nrow(z)==0)
					z	<- subset(tmp, OVERLAP>149)
				stopifnot(nrow(z)>0)
				z	<- z[, list(PD=min(PD)), by='TAXA2']
				setkey(z, PD)				
				tmp	<- integer(0)
				if(length(TAXA)<22L)
					tmp	<- seq_len(22L-length(TAXA))
				list(	TAXA= c(TAXA, z[tmp, TAXA2]), 
						PARTNER= c(PARTNER, rep('Other',length(tmp))), 
						COUPID=c(COUPID,rep('Other',length(tmp))) )				
			}, by='PTY_RUN']
	ans			<- merge(ans, subset(dm, select=c(TAXA, SANGER_ID)), by='TAXA')
	setnames(ans, 'SANGER_ID', 'FILE_ID')	
	#	to each run, add one individual with super coverage: 15172_1_43, PG14-UG500280-S02284
	ans			<- rbind(ans, as.data.table(expand.grid(TAXA='PG14-UG500280-S02284', FILE_ID='15172_1_43', PARTNER='Other', COUPID='Other', PTY_RUN=ans[, sort(unique(PTY_RUN))])))
	setkey(ans, PTY_RUN, TAXA)	
	pty.runs	<- copy(ans)
	save(pty.runs, file=file.path('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples','Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda'))
	#	add second dummy individual
	ans			<- rbind(ans, as.data.table(expand.grid(TAXA='PG14-UG501310-S03314', FILE_ID='15777_1_5', PARTNER='Other', COUPID='Other', PTY_RUN=ans[, sort(unique(PTY_RUN))])))
	setkey(ans, PTY_RUN, TAXA)
	pty.runs	<- copy(ans)
	save(pty.runs, file=file.path('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples','Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns_2dummy.rda'))
}
######################################################################################
project.Rakai.checkMissingRakai.150307<- function()
{
	png.f	<- '~/Dropbox (SPH Imperial College)/PANGEA_data/2016-01-20_PANGEA_stats_by_sample.csv'	
	#infile	<- '~/Dropbox (SPH Imperial College)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/CheckPangeaId.csv'
	#infile	<- '~/Dropbox (SPH Imperial College)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/CheckPangeaId2.csv'
	infile	<- '~/Dropbox (SPH Imperial College)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/SequenceDataNOTAvailable_Documented.csv'
	infile	<- '~/Dropbox (SPH Imperial College)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/SequenceDataAvailable_NotDocumented.csv'
	
	#
	#
	#
	dpng	<- as.data.table(read.csv(png.f))
	setnames(dpng, 'ProjectID', 'PANGEA_ID')
	df		<- as.data.table(read.csv(infile))
	setnames(df, 'x', 'PANGEA_ID')
	df		<- merge(dpng, df, by='PANGEA_ID')	
	write.csv(df, row.names=FALSE, file=gsub('.csv','_info.csv',infile))
	
	png.f	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/PANGEA_orig/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
	load(png.f)	#sq, sqi, si
	set(si, NULL, 'PANGEA_ID', si[, regmatches(PANGEA_ID, regexpr('PG[0-9]+-[^-]+',PANGEA_ID))])
	merge(si, df, by='PANGEA_ID')		
}
#
#	plots of individual histories
#	
RakaiCirc.circ.timelines.plots<- function(rt, wdir)
{
#	plot total circumcision by year among individuals in follow up
	setkey(rt, RID, ROUND)
	ggplot(unique(rt), aes(x=factor(ROUND_ST), fill=CIRC)) + geom_bar() + 
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(breaks=seq(0,1e4, 5e2)) +
			facet_grid(SEX~.) +
			labs(x='\nsurvey rounds\n(start date)', y='HIV infected RCCS participants\n(reported alive, incl temporary loss to follow-up)\n', fill='Circumcision\nstatus')
	ggsave(file.path(wdir, 'RCCSinf_by_circ.pdf'), w=10, h=8)
	#	plot proportion circumcision by year among individuals in follow up
	ggplot(subset(unique(rt), SEX=='male'), aes(x=factor(ROUND_ST), fill=CIRC)) + geom_bar(position='fill') + 
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(breaks=seq(0,1, 0.1), labels = scales::percent) +
			labs(x='\nsurvey rounds\n(start date)', y='HIV infected male RCCS participants\n(reported alive, incl temporary loss to follow-up)\n', fill='Circumcision\nstatus')
	ggsave(file.path(wdir, 'RCCSinf_by_circ_proportions_men.pdf'), w=10, h=5)	
	#	plot sequence coverage among HIV infected participants
	tmp		<- unique(rt)
	set(tmp, NULL, 'SEQ_TYPE', tmp[, factor(SEQ_TYPE, levels=c('gag_or_partial_gag_PANGEA','gag_or_partial_gag_historic_only','Sanger_completed_with_IVA','Sanger_completed_without_IVA','Sanger_not_started_by_Jul2016','other_gene_historic_only','Sanger_failed','no_sequence'))])
	ggplot(tmp, aes(x=factor(ROUND_ST), fill=SEQ_TYPE)) + geom_bar() + 
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(breaks=seq(0,1e4, 5e2)) +
			scale_fill_brewer(palette='Set2') +
			labs(x='\nsurvey rounds\n(start date)', y='HIV infected RCCS participants\n(reported alive, incl temporary loss to follow-up)\n', fill='Sequence\nstatus')
	ggsave(file.path(wdir, 'RCCSinf_by_sequence_coverage.pdf'), w=10, h=5)
	#	plot proportion sequence coverage among HIV infected participants
	ggplot(tmp, aes(x=factor(ROUND_ST), fill=SEQ_TYPE)) + geom_bar(position='fill') + 
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(breaks=seq(0,1, 0.1), labels = scales::percent) +
			scale_fill_brewer(palette='Set3') +
			labs(x='\nsurvey rounds\n(start date)', y='HIV infected RCCS participants\n(reported alive, incl temporary loss to follow-up)\n', fill='Sequence\nstatus')
	ggsave(file.path(wdir, 'RCCSinf_by_sequence_coverage_proportions.pdf'), w=10, h=5)	
	#	plot circumcised by religion
	tmp		<- subset(unique(rt), EVERCIRC=='Y')
	ggplot(tmp, aes(x=factor(ROUND_ST), fill=RELIGION)) + geom_bar(position='fill') + 
			theme_bw() +
			theme(legend.position='bottom') +
			scale_y_continuous(breaks=seq(0,1, 0.1), labels = scales::percent) +
			scale_fill_brewer(palette='Set2') +
			labs(x='\nsurvey rounds\n(start date)', y='circumcised male \nHIV infected RCCS participants\n(reported alive, incl temporary loss to follow-up)\n', fill='Faith')
	ggsave(file.path(wdir, 'RCCScircmen_by_religion_proportions.pdf'), w=10, h=5)		
}

RakaiCirc.various<- function()
{
	require(ape)
	require(data.table)
	if(0)
	{
		wdir	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/circumcision"
		wdir	<- '/work/or105/Gates_2014/Rakai'
		load(file.path(wdir,'RCCS_PhInfo_160825.rda'))
		
		#	add index for speed
		tmp			<- data.table(SEQIDc=rownames(sq), IDX=seq_len(nrow(sq)))
		gdgag		<- merge(tmp, gdgag, by='SEQIDc')
		setnames(tmp, c('SEQIDc','IDX'), c('SEQIDc2','IDX2'))
		gdgag		<- merge(tmp, gdgag, by='SEQIDc2')
		#	
		tmp			<- (as.character(sq)!='-')
		setkey(gdgag, IDX, IDX2)
		tmp2		<- gdgag[, list(OVERLAP= sum(apply(tmp[c(IDX,IDX2),],2,all)) ), by=c('IDX','IDX2')]
		gdgag		<- merge(gdgag, tmp2, by=c('SEQIDc','SEQIDc2'))
		
		save(phf, php, phfi, phpi, sq, gdgag, file=file.path(wdir,'RCCS_PhInfo_2_160825.rda'))
	}
	if(1)
	{
		require(ape)
		require(data.table)
		
		wdir	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/circumcision"
		wdir	<- '/work/or105/Gates_2014/Rakai'
		infile	<- file.path(wdir, "PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda")
		load(infile)
		#	prepare genetic distance matrix
		sq			<- as.character(sq)
		sq[sq=='?']	<- '-'
		sq			<- as.DNAbin(sq)
		sq.gd		<- dist.dna(sq, pairwise.deletion=TRUE)
		sq.gd		<- as.data.table(melt(as.matrix(sq.gd), varnames=c('TAXA1','TAXA2')))
		setnames(sq.gd, 'value', 'PD')
		sq.gd		<- subset(sq.gd, TAXA1!=TAXA2)
		set(sq.gd, NULL, 'TAXA1', sq.gd[, as.character(TAXA1)])
		set(sq.gd, NULL, 'TAXA2', sq.gd[, as.character(TAXA2)])		
		#	add overlap	
		setkey(sq.gd, TAXA1, TAXA2)
		tmp			<- as.character(sq)
		tmp			<- !( tmp=='-' | tmp=='?' | tmp=='n' )
		tmp[]		<- as.integer(tmp) 
		tmp2		<- sq.gd[, list(OVERLAP= sum(bitwAnd(tmp[TAXA1,], tmp[TAXA2,])) ), by=c('TAXA1','TAXA2')]	
		sq.gd		<- merge(sq.gd, tmp2, by=c('TAXA1','TAXA2'))
		
		save(sq.gd, file=gsub('\\.rda','_gd.rda',infile))		
	}	
}

RakaiCirc.seq.get.phylogenies<- function()
{
	require(ape)
	require(data.table)
	wdir	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/circumcision"
	#
	#	start with Kate s GP24 FastTree
	#
	infile	<- '~/Dropbox (SPH Imperial College)/PANGEA_alignments/Regional phylogenetic analyses/Region 1 gag analysis/gagAllSmall.nwk'
	phf		<- read.tree(infile)
	phfi	<- data.table(SEQIDb= phf$tip.label)
	#	get patristic distances
	tmp								<- cophenetic.phylo(phf)
	tmp								<- as.matrix(tmp)
	tmp[upper.tri(tmp, diag=TRUE)]	<- NA_real_
	tmp								<- as.data.table(melt(tmp))
	setnames(tmp, c('Var1','Var2','value'), c('SEQIDb','SEQIDb2','PD'))
	phfi							<- subset(tmp, !is.na(PD))
	#
	#	also use the PANGEA ExaML tree
	#
	infile	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_UG_gag_fasttree.rda'
	load(infile)	#loads "ph", "dist.brl", "ph.gdtr", "ph.mrca", "clustering"
	php		<- ph
	tmp		<- copy(ph.gdtr)
	tmp[upper.tri(tmp, diag=TRUE)]	<- NA_real_
	tmp								<- as.data.table(melt(tmp))
	setnames(tmp, c('Var1','Var2','value'), c('PID','PID2','PD'))
	phpi							<- subset(tmp, !is.na(PD))
	#
	#	get raw genetic distances on latest Region1 alignment
	#
	infile	<- '~/Dropbox (SPH Imperial College)/PANGEA_alignments/Regional Alignments/150825_Region1UG_codonaligned_p_PANGEA151113_p_HXB2.rda'
	load(infile)	
	sq			<- as.character(sq)
	sq[sq=='?']	<- '-'
	sq			<- as.DNAbin(sq)
	tmp			<- dist.dna(sq, model='raw', pairwise.deletion=TRUE)
	tmp								<- as.matrix(tmp)
	tmp[upper.tri(tmp, diag=TRUE)]	<- NA_real_
	tmp								<- as.data.table(melt(tmp))
	setnames(tmp, c('Var1','Var2','value'), c('SEQIDc','SEQIDc2','GDRW'))
	gdgag		<- subset(tmp, !is.na(GDRW))
	set(gdgag, NULL, 'SEQIDc', gdgag[, as.character(SEQIDc)])
	set(gdgag, NULL, 'SEQIDc2', gdgag[, as.character(SEQIDc2)])
	#	add index for speed
	tmp			<- data.table(SEQIDc=rownames(sq), IDX=seq_len(nrow(sq)))
	gdgag		<- merge(tmp, gdgag, by='SEQIDc')
	setnames(tmp, c('SEQIDc','IDX'), c('SEQIDc2','IDX2'))
	gdgag		<- merge(tmp, gdgag, by='SEQIDc2')
	#	add overlap
	#	this takes about 2 days.. ..ran on cluster.
	tmp			<- (as.character(sq)!='-')
	setkey(gdgag, IDX, IDX2)
	tmp2		<- gdgag[, list(OVERLAP= sum(apply(tmp[c(IDX,IDX2),],2,all)) ), by=c('IDX','IDX2')]
	gdgag		<- merge(gdgag, tmp2, by=c('SEQIDc','SEQIDc2'))
	save(phf, php, phfi, phpi, sq, gdgag, file=file.path(wdir,'RCCS_PhInfo_160825.rda'))	
	#
	#	plot nucleotide frequencies (unfinished)
	#
	if(0)
	{
		sqf			<- data.table(SITE=seq_len(ncol(sq)))
		sqf			<- sqf[, {
					tmp	<- base.freq(sq[,SITE], freq=TRUE, all=TRUE)
					list(NT=paste('NT',names(tmp),sep=''), N=tmp)
				}, by='SITE']
		sqf			<- subset(sqf, N>0)
		tmp			<- subset(sqf,!NT%in%c('NTn','NT-','NT?'))[, list(TOTAL_SEQ=sum(N)), by='SITE']
		sqf			<- merge(sqf, tmp, all.x=1, by='SITE')
		sqf[, NTc:= 'NTother']
		tmp			<- sqf[, which(NT%in%c('NTa','NTt','NTg','NTc','NT-','NT?'))]
		set(sqf, tmp, 'NTc', sqf[tmp, NT])
		set(sqf, sqf[, which(NT=='NTn')], 'NTc', 'NT?')
		ggplot(subset(sqf, SITE<1000 & NTc!='NT?' & NTc!='NT-'), aes(x=SITE, fill=NTc)) + 
				geom_bar(position='fill') +
				scale_fill_brewer(palette='Set1')
	}
}

RakaiFull.preprocess.closepairs.calculate.170421	<- function()
{
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
	dmin	<- infiles[,{
				#cat(PTY_RUN,'\n')
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl5_blnormed/ptyr1_pairwise_relationships.rda'
				load(F)
				tmp		<- dwin[, list(PATRISTIC_DISTANCE_MIN=min(PATRISTIC_DISTANCE), PATRISTIC_DISTANCE_MEAN=mean(PATRISTIC_DISTANCE)), by=c('ID1','ID2')]
				tmp
			}, by='PTY_RUN']	
	save(rtp, dmin, file=outfile)	
}

RakaiFull.preprocess.couples.todi.phyloscanneroutput.170421<- function()
{
	#	load sequence data
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/circumcision/RCCS_SeqInfo_170505.rda")
	setnames(rs, 'SAMPLE_DATE', 'SEQ_DATE')
	use.posterior.mode		<- 1
	#use.direction.prior.23	<- 1
	#
	#	load demographic info on all individuals
	tmp		<- RakaiCirc.epi.get.info.170208()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#rn		<- tmp$rn
	ra		<- tmp$ra
	#rd		<- rbind(rd, rn, use.names=TRUE, fill=TRUE)
	set(rd, NULL, c('PID','SID'), NULL)	
	set(rd, NULL, 'SEX', rd[, as.character(SEX)])
	set(rd, NULL, 'RECENTVL', rd[, as.numeric(gsub('< 150','1',gsub('> ','',gsub('BD','',gsub(',','',as.character(RECENTVL))))))])
	set(rd, NULL, 'CAUSE_OF_DEATH', rd[, as.character(CAUSE_OF_DEATH)])
	#	fixup rd: 
	#	remove HIV reverters without sequence
	rd		<- subset(rd, !RID%in%c("C117824","C119303","E118889","K067249"))
	#	fixup complex serology
	set(rd, rd[, which(RID=='B106184')], 'FIRSTPOSDATE', rd[which(RID=='B106184'),DATE])
	set(rd, rd[, which(RID=='B106184')], c('LASTNEGVIS','LASTNEGDATE'), NA_real_)
	set(rd, rd[, which(RID=='B106184')], c('HIVPREV'), 1)
	set(rd, rd[, which(RID=='A008742')], 'FIRSTPOSDATE', rd[which(RID=='A008742'),DATE])
	set(rd, rd[, which(RID=='A008742')], c('HIVPREV'), 1)
	#	fixup rd: 
	#	missing first pos date
	rd		<- subset(rd, RID!='A038432')	#has missing firstposdate and not in PANGEA anyway
	rd		<- subset(rd, RID!='H013226')	#has missing firstposdate and not in PANGEA anyway
	rd		<- subset(rd, RID!='K008173')	#has missing firstposdate and not in PANGEA anyway
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
						if(any(!is.na(FIRSTPOSVIS)))
							ans[is.na(FIRSTPOSVIS) | (!is.na(FIRSTPOSVIS) & FIRSTPOSVIS!=min(FIRSTPOSVIS, na.rm=TRUE))]	<- TRUE							
						list(DUMMY=DUMMY, DELETE=ans)		
					}, by=c('RID')], by=c('RID','DUMMY'))
	tmp		<- subset(tmp, !DELETE)
	set(tmp, NULL,c('DUMMY','DELETE'), NULL)
	set(rd, NULL, c('BIRTHDATE','LASTNEGDATE','FIRSTPOSVIS','FIRSTPOSDATE','ARVSTARTDATE','EST_DATEDIED'), NULL)
	rd		<- merge(rd, tmp, by='RID')	
	tmp		<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
	stopifnot(!nrow(merge(subset(tmp[, length(BIRTHDATE), by='RID'], V1>1), tmp, by='RID')))	
	
		
	#
	#	load couples to search for in phyloscanner output
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_info.rda")
	rps		<- unique(subset(rp, select=c(FEMALE_RID, MALE_RID, COUPID)))
	rc		<- melt(rps, id.vars='COUPID', measure.vars=c('FEMALE_RID','MALE_RID'), value.name='RID', variable.name='SEX')
	set(rc, NULL, 'SEX', rc[, substr(SEX,1,1)])
	tmp		<- subset(rps, select=c('FEMALE_RID','MALE_RID'))
	setnames(tmp, c('FEMALE_RID','MALE_RID'), c('MALE_RID','FEMALE_RID'))
	rps		<- rbind(subset(rps, select=c('FEMALE_RID','MALE_RID')), tmp)
	setnames(rps, c('FEMALE_RID','MALE_RID'), c('ID1','ID2'))
		
	#
	#	from every phyloscanner run, select pairs that are closely related 
	#
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_couples_170428_cl3.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun34'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_couples_170516_cl3.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun23'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170610/todi_couples_170610_cl3_prior23.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun34d23'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170610/todi_couples_170610_cl3_prior34d23.rda'
	
	
	infiles	<- data.table(F=list.files(indir, pattern='pairwise_relationships.rda', full.names=TRUE))
	infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(F)))]
	setkey(infiles, PTY_RUN)
	rtp.todi2.basic	<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23/ptyr2_pairwise_relationships.rda'
				cat(PTY_RUN,'\n')
				load(F)
				#	search for couples in this run, and keep if these are most likely not a couple
				rtp		<- merge(rps, subset(rplkl, GROUP=='TYPE_PAIR_TODI2'), by=c('ID1','ID2'))
				rtp		<- rtp[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, TYPE_MLE=='other')
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIR_TODI2' & TYPE=='other'), by=c('ID1','ID2'), all.x=1)
				if(!use.posterior.mode)
					rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				if(use.posterior.mode)
					rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]																
				rtp[, SELECT:= 'couple most likely not a pair']
				ans		<- copy(rtp)
				#	ML likely transmission pairs by distance
				rtp		<- merge(rps, subset(rplkl, GROUP=='TYPE_PAIR_DI'), by=c('ID1','ID2'))
				rtp		<- rtp[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, TYPE_MLE=='close')
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIR_DI' & TYPE=='close'), by=c('ID1','ID2'), all.x=1)
				if(!use.posterior.mode)
					rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				if(use.posterior.mode)
					rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]																
				rtp[, SELECT:= 'couple most likely a pair']
				ans		<- rbind(ans, rtp, use.names=TRUE)
				#	ML likely transmission pairs by topology
				rtp		<- merge(rps, subset(rplkl, GROUP=='TYPE_PAIR_TO'), by=c('ID1','ID2'))
				rtp		<- rtp[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, TYPE_MLE=='ancestral/\nintermingled')
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIR_TO' & TYPE=='ancestral/\nintermingled'), by=c('ID1','ID2'), all.x=1)
				if(!use.posterior.mode)
					rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				if(use.posterior.mode)
					rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]																
				rtp[, SELECT:= 'couple most likely a pair']
				ans		<- rbind(ans, rtp, use.names=TRUE)
				#	find most likely trm pairs based on TODI2
				rtp		<- merge(rps, subset(rplkl, GROUP=='TYPE_PAIR_TODI2'), by=c('ID1','ID2'))
				rtp		<- rtp[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, TYPE_MLE=='likely pair')
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIR_TODI2' & TYPE=='likely pair'), by=c('ID1','ID2'), all.x=1)
				if(!use.posterior.mode)
					rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				if(use.posterior.mode)
					rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]																
				rtp[, SELECT:= 'couple most likely a pair']
				ans		<- rbind(ans, rtp, use.names=TRUE)
				#	ML directed likely transmission pairs by distance + topology
				#	among those pairs that are likely transmission pairs
				rtp		<- merge(subset(rtp, select=c('ID1','ID2')), subset(rplkl, GROUP=='TYPE_DIR_TODI3'), by=c('ID1','ID2'))
				rtp		<- rtp[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, TYPE_MLE!='ambiguous')
				rtp[, TYPE:=TYPE_MLE]				
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_DIRSCORE_TODI3'), by=c('ID1','ID2','TYPE'), all.x=1)
				#if(use.direction.prior.23)
				#{
				#	set(rtp, NULL, c('PAR_PRIOR','POSTERIOR_ALPHA','POSTERIOR_BETA'), NULL)
				#	rtp[, PAR_PRIOR:= phsc.get.prior.parameter.n0(2, keff=2, neff=3, confidence.cut=0.66)]					
				#	rtp[, POSTERIOR_ALPHA:= PAR_PRIOR/N_TYPE+KEFF]
				#	rtp[, POSTERIOR_BETA:= PAR_PRIOR*(1-1/N_TYPE)+NEFF-KEFF]						
				#}									
				if(!use.posterior.mode)
					rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				if(use.posterior.mode)
					rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]																
				rtp[, SELECT:= 'couple most likely a pair with resolved direction']				
				ans		<- rbind(ans, rtp, use.names=TRUE)
				ans				
			}, by=c('PTY_RUN')]		
	#rtp.todi2.basic	<- copy(rtp.todi2)
	#
	#	re-arrange to male-female
	#
	tmp			<- unique(subset(rd, select=c('RID','SEX')))
	setnames(tmp, colnames(tmp), paste0('ID1_',colnames(tmp)))
	setnames(tmp, c('ID1_RID'), c('ID1'))
	stopifnot( !length(setdiff(rtp.todi2[, ID1], tmp[, ID1])) )
	rtp.todi2	<- merge(rtp.todi2.basic, tmp, by=c('ID1'))	
	setnames(tmp, colnames(tmp), gsub('ID1','ID2',colnames(tmp)))
	stopifnot( !length(setdiff(rtp.todi2[, ID2], tmp[, ID2])) )	
	rtp.todi2	<- merge(rtp.todi2, tmp, by=c('ID2'))
	tmp			<- rtp.todi2[, which(TYPE=='12')]
	set(rtp.todi2, tmp, 'TYPE', rtp.todi2[tmp, tolower(paste0(ID1_SEX,ID2_SEX))])
	tmp			<- rtp.todi2[, which(TYPE=='21')]
	set(rtp.todi2, tmp, 'TYPE', rtp.todi2[tmp, tolower(paste0(ID2_SEX,ID1_SEX))])
	#
	#	prepare just the dwin and rplkl for couples
	#
	rplkl		<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr1_pairwise_relationships.rda'
				cat(PTY_RUN,'\n')
				load(F)
				ans	<- merge(unique(subset(rps, select=c('ID1','ID2'))), rplkl, by=c('ID1','ID2'))
				ans			
			}, by='PTY_RUN']
	rpw		<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr1_pairwise_relationships.rda'
				cat(PTY_RUN,'\n')
				load(F)
				ans	<- merge(unique(subset(rps, select=c('ID1','ID2'))), dwin, by=c('ID1','ID2'))
				ans			
			}, by='PTY_RUN']
	rpw			<- melt(rpw, variable.name='GROUP', value.name='TYPE', measure.vars=c("TYPE_RAW","TYPE_BASIC","TYPE_DIR_TODI3","TYPE_DIRSCORE_TODI3","TYPE_DIR_TODI4","TYPE_PAIR_TODI2","TYPE_CHAIN_TODI","TYPE_PAIR_DI2","TYPE_PAIR_TO","TYPE_PAIR_TODI2x2"))
	set(rpw, NULL, 'ID_R_MAX', rpw[, pmax(ID1_R,ID2_R)])
	set(rpw, NULL, 'ID_R_MIN', rpw[, pmin(ID1_R,ID2_R)])			
	#	save
	save(rp, rd, rh, ra, rs, rtp.todi2.basic, rtp.todi2, rplkl, rpw, file=outfile)	
}

RakaiFull.preprocess.couples.todi.phyloscanneroutput.170811<- function()
{
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_couples_170428_cl3.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun34'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_couples_170516_cl3.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun23'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170610/todi_couples_170610_cl3_prior23.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun34d23'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170610/todi_couples_170610_cl3_prior34d23.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min10'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min10_zbl'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10_zbl.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min50'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min50.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p25_stagetwo_rerun23_min30'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl25_prior23_min30.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p45_stagetwo_rerun23_min30'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl45_prior23_min30.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_d30'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_d30.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_d100'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_d100.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_d1000'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_d1000.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s10_p35_stagetwo_rerun23_min30'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_s10.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s40_p35_stagetwo_rerun23_min30'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_s40.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_adj'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_adj.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_mf3'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_mf3.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_mf4'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_mf4.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_mf6'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_mf6.rda'	
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_mt'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_mt.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_prt'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_prt.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_rg1'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_rg1.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_rg20'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_rg20.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_zbl'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_zbl.rda'
	
	
	#	load sequence data
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/circumcision/RCCS_SeqInfo_170505.rda")
	setnames(rs, 'SAMPLE_DATE', 'SEQ_DATE')
	use.posterior.mode		<- 1
	#use.direction.prior.23	<- 1
	#
	#	load demographic info on all individuals
	tmp		<- RakaiCirc.epi.get.info.170208()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#rn		<- tmp$rn
	ra		<- tmp$ra
	#rd		<- rbind(rd, rn, use.names=TRUE, fill=TRUE)
	set(rd, NULL, c('PID','SID'), NULL)	
	set(rd, NULL, 'SEX', rd[, as.character(SEX)])
	set(rd, NULL, 'RECENTVL', rd[, as.numeric(gsub('< 150','1',gsub('> ','',gsub('BD','',gsub(',','',as.character(RECENTVL))))))])
	set(rd, NULL, 'CAUSE_OF_DEATH', rd[, as.character(CAUSE_OF_DEATH)])
	#	fixup rd: 
	#	remove HIV reverters without sequence
	rd		<- subset(rd, !RID%in%c("C117824","C119303","E118889","K067249"))
	#	fixup complex serology
	set(rd, rd[, which(RID=='B106184')], 'FIRSTPOSDATE', rd[which(RID=='B106184'),DATE])
	set(rd, rd[, which(RID=='B106184')], c('LASTNEGVIS','LASTNEGDATE'), NA_real_)
	set(rd, rd[, which(RID=='B106184')], c('HIVPREV'), 1)
	set(rd, rd[, which(RID=='A008742')], 'FIRSTPOSDATE', rd[which(RID=='A008742'),DATE])
	set(rd, rd[, which(RID=='A008742')], c('HIVPREV'), 1)
	#	fixup rd: 
	#	missing first pos date
	rd		<- subset(rd, RID!='A038432')	#has missing firstposdate and not in PANGEA anyway
	rd		<- subset(rd, RID!='H013226')	#has missing firstposdate and not in PANGEA anyway
	rd		<- subset(rd, RID!='K008173')	#has missing firstposdate and not in PANGEA anyway
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
						if(any(!is.na(FIRSTPOSVIS)))
							ans[is.na(FIRSTPOSVIS) | (!is.na(FIRSTPOSVIS) & FIRSTPOSVIS!=min(FIRSTPOSVIS, na.rm=TRUE))]	<- TRUE							
						list(DUMMY=DUMMY, DELETE=ans)		
					}, by=c('RID')], by=c('RID','DUMMY'))
	tmp		<- subset(tmp, !DELETE)
	set(tmp, NULL,c('DUMMY','DELETE'), NULL)
	set(rd, NULL, c('BIRTHDATE','LASTNEGDATE','FIRSTPOSVIS','FIRSTPOSDATE','ARVSTARTDATE','EST_DATEDIED'), NULL)
	rd		<- merge(rd, tmp, by='RID')	
	tmp		<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
	stopifnot(!nrow(merge(subset(tmp[, length(BIRTHDATE), by='RID'], V1>1), tmp, by='RID')))	
	
	
	#
	#	load couples to search for in phyloscanner output
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_info.rda")
	rps		<- unique(subset(rp, select=c(FEMALE_RID, MALE_RID, COUPID)))
	rc		<- melt(rps, id.vars='COUPID', measure.vars=c('FEMALE_RID','MALE_RID'), value.name='RID', variable.name='SEX')
	set(rc, NULL, 'SEX', rc[, substr(SEX,1,1)])
	tmp		<- subset(rps, select=c('FEMALE_RID','MALE_RID'))
	setnames(tmp, c('FEMALE_RID','MALE_RID'), c('MALE_RID','FEMALE_RID'))
	rps		<- rbind(subset(rps, select=c('FEMALE_RID','MALE_RID')), tmp)
	setnames(rps, c('FEMALE_RID','MALE_RID'), c('ID1','ID2'))
	
	#
	#	from every phyloscanner run, select pairs that are closely related 
	#	
	tmp		<- data.table(F=list.files(indir, pattern='_pairwise_relationships.rda$', full.names=TRUE))
	tmp[, PTY_RUN:= as.integer(gsub('ptyr([0-9]+)_.*','\\1',basename(F)))]
	paste(setdiff( 2:346, tmp[, PTY_RUN]), collapse=', ')
	
	infiles	<- data.table(F=list.files(indir, pattern='pairwise_relationships.rda', full.names=TRUE))
	infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(F)))]
	setkey(infiles, PTY_RUN)
	rtp.todi2.basic	<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min10/ptyr346_pairwise_relationships.rda'
				cat(PTY_RUN,'\n')
				load(F)
				#	search for couples in this run, and keep if these are most likely not a couple
				rtp		<- merge(rps, subset(rplkl, GROUP=='TYPE_PAIR_TODI2'), by=c('ID1','ID2'))
				rtp		<- rtp[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, TYPE_MLE=='unlinked')
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIR_TODI2' & TYPE=='unlinked'), by=c('ID1','ID2'), all.x=1)
				if(!use.posterior.mode)
					rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				if(use.posterior.mode)
					rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]																
				rtp[, SELECT:= 'couple most likely not a pair']
				ans		<- copy(rtp)
				#	ML likely transmission pairs by distance
				rtp		<- merge(rps, subset(rplkl, GROUP=='TYPE_PAIR_DI2'), by=c('ID1','ID2'))
				rtp		<- rtp[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, TYPE_MLE=='close')
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIR_DI2' & TYPE=='close'), by=c('ID1','ID2'), all.x=1)
				if(!use.posterior.mode)
					rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				if(use.posterior.mode)
					rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]																
				rtp[, SELECT:= 'couple most likely a pair']
				ans		<- rbind(ans, rtp, use.names=TRUE)
				#	ML likely transmission pairs by topology
				rtp		<- merge(rps, subset(rplkl, GROUP=='TYPE_PAIR_TO'), by=c('ID1','ID2'))
				rtp		<- rtp[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, TYPE_MLE=='ancestral/\nintermingled')
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIR_TO' & TYPE=='ancestral/\nintermingled'), by=c('ID1','ID2'), all.x=1)
				if(!use.posterior.mode)
					rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				if(use.posterior.mode)
					rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]																
				rtp[, SELECT:= 'couple most likely a pair']
				ans		<- rbind(ans, rtp, use.names=TRUE)
				#	find most likely trm pairs based on TODI2
				rtp		<- merge(rps, subset(rplkl, GROUP=='TYPE_PAIR_TODI2'), by=c('ID1','ID2'))
				rtp		<- rtp[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, TYPE_MLE=='linked')
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIR_TODI2' & TYPE=='linked'), by=c('ID1','ID2'), all.x=1)
				if(!use.posterior.mode)
					rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				if(use.posterior.mode)
					rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]																
				rtp[, SELECT:= 'couple most likely a pair']
				ans		<- rbind(ans, rtp, use.names=TRUE)
				#	ML directed likely transmission pairs by distance + topology
				#	among those pairs that are likely transmission pairs
				rtp		<- merge(subset(rtp, select=c('ID1','ID2')), subset(rplkl, GROUP=='TYPE_NETWORK_SCORES'), by=c('ID1','ID2'))
				rtp		<- rtp[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, !TYPE_MLE%in%c('ambiguous','not close/disconnected'))
				rtp[, TYPE:=TYPE_MLE]				
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_DIR_TODI2'), by=c('ID1','ID2','TYPE'), all.x=1)
				#if(use.direction.prior.23)
				#{
				#	set(rtp, NULL, c('PAR_PRIOR','POSTERIOR_ALPHA','POSTERIOR_BETA'), NULL)
				#	rtp[, PAR_PRIOR:= phsc.get.prior.parameter.n0(2, keff=2, neff=3, confidence.cut=0.66)]					
				#	rtp[, POSTERIOR_ALPHA:= PAR_PRIOR/N_TYPE+KEFF]
				#	rtp[, POSTERIOR_BETA:= PAR_PRIOR*(1-1/N_TYPE)+NEFF-KEFF]						
				#}									
				if(!use.posterior.mode)
					rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				if(use.posterior.mode)
					rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]																
				rtp[, SELECT:= 'couple most likely a pair with resolved direction']				
				ans		<- rbind(ans, rtp, use.names=TRUE)
				ans				
			}, by=c('PTY_RUN')]		
	#rtp.todi2.basic	<- copy(rtp.todi2) ZZZ
	#
	#	re-arrange to male-female
	#
	tmp			<- unique(subset(rd, select=c('RID','SEX')))
	setnames(tmp, colnames(tmp), paste0('ID1_',colnames(tmp)))
	setnames(tmp, c('ID1_RID'), c('ID1'))
	stopifnot( !length(setdiff(rtp.todi2.basic[, ID1], tmp[, ID1])) )
	rtp.todi2	<- merge(rtp.todi2.basic, tmp, by=c('ID1'))	
	setnames(tmp, colnames(tmp), gsub('ID1','ID2',colnames(tmp)))
	stopifnot( !length(setdiff(rtp.todi2[, ID2], tmp[, ID2])) )	
	rtp.todi2	<- merge(rtp.todi2, tmp, by=c('ID2'))
	tmp			<- rtp.todi2[, which(TYPE=='12')]
	set(rtp.todi2, tmp, 'TYPE', rtp.todi2[tmp, tolower(paste0(ID1_SEX,ID2_SEX))])
	tmp			<- rtp.todi2[, which(TYPE=='21')]
	set(rtp.todi2, tmp, 'TYPE', rtp.todi2[tmp, tolower(paste0(ID2_SEX,ID1_SEX))])
	
	#
	#	prepare just the dwin and rplkl for couples
	#
	rplkl		<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr1_pairwise_relationships.rda'
				cat(PTY_RUN,'\n')
				load(F)
				ans	<- merge(unique(subset(rps, select=c('ID1','ID2'))), rplkl, by=c('ID1','ID2'))
				ans			
			}, by='PTY_RUN']
	rpw		<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr1_pairwise_relationships.rda'
				cat(PTY_RUN,'\n')
				load(F)
				ans	<- merge(unique(subset(rps, select=c('ID1','ID2'))), dwin, by=c('ID1','ID2'))
				ans			
			}, by='PTY_RUN']
	rpw			<- melt(rpw, variable.name='GROUP', value.name='TYPE', measure.vars=c("TYPE_RAW","TYPE_BASIC","TYPE_PAIR_DI2","TYPE_PAIR_TO","TYPE_PAIR_TODI2x2","TYPE_PAIR_TODI2","TYPE_DIR_TODI2","TYPE_NETWORK_SCORES","TYPE_CHAIN_TODI"))
	set(rpw, NULL, 'ID_R_MAX', rpw[, pmax(ID1_R,ID2_R)])
	set(rpw, NULL, 'ID_R_MIN', rpw[, pmin(ID1_R,ID2_R)])			
	#	save
	save(rp, rd, rh, ra, rs, rtp.todi2.basic, rtp.todi2, rplkl, rpw, file=outfile)	
}

RakaiFull.preprocess.trmpairs.community.sampling.170421<- function()
{
	require(data.table)
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/community_sampling_170522.rda'
	
	#
	#	load couples to search for in phyloscanner output
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_info.rda")
	#
	#	load sequence data
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/circumcision/RCCS_SeqInfo_170505.rda")
	setnames(rs, 'SAMPLE_DATE', 'SEQ_DATE')
	#
	#	load demographic info on all individuals
	tmp		<- RakaiCirc.epi.get.info.170208()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#rn		<- tmp$rn
	ra		<- tmp$ra
	#set(rn, NULL, 'RID', rn[, as.character(RID)])
	#rn		<- merge(rn, subset(rd, select=c(RID, FIRSTPOSVIS, FIRSTPOSDATE)), by='RID', all.x=1)
	#tmp		<- rn[, which(is.na(FIRSTPOSDATE) & !is.na(RECENTVLDATE) & TIMESINCEVL==0)]	#this is dodgy
	#set(rn, tmp, 'FIRSTPOSDATE', rn[tmp, RECENTVLDATE])		
	#rd		<- rbind(rd, rn, use.names=TRUE, fill=TRUE)	#do not consider individuals in the neuro study that are not part of RCCS
	set(rd, NULL, c('PID','SID'), NULL)
	set(rd, NULL, 'SEX', rd[, as.character(SEX)])
	set(rd, NULL, 'RECENTVL', rd[, as.numeric(gsub('< 150','1',gsub('> ','',gsub('BD','',gsub(',','',as.character(RECENTVL))))))])
	set(rd, NULL, 'CAUSE_OF_DEATH', rd[, as.character(CAUSE_OF_DEATH)])
	#	fixup rd: 
	#	remove HIV reverters without sequence
	rd		<- subset(rd, !RID%in%c("C117824","C119303","E118889","K067249"))
	#	fixup complex serology
	set(rd, rd[, which(RID=='B106184')], 'FIRSTPOSDATE', rd[which(RID=='B106184'),DATE])
	set(rd, rd[, which(RID=='B106184')], c('LASTNEGVIS','LASTNEGDATE'), NA_real_)
	set(rd, rd[, which(RID=='B106184')], c('HIVPREV'), 1)
	set(rd, rd[, which(RID=='A008742')], 'FIRSTPOSDATE', rd[which(RID=='A008742'),DATE])
	set(rd, rd[, which(RID=='A008742')], c('HIVPREV'), 1)
	#	fixup rd: 
	#	missing first pos date
	rd		<- subset(rd, RID!='A038432')	#has missing firstposdate and not in PANGEA anyway
	rd		<- subset(rd, RID!='H013226')	#has missing firstposdate and not in PANGEA anyway
	rd		<- subset(rd, RID!='K008173')	#has missing firstposdate and not in PANGEA anyway
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
						if(any(!is.na(FIRSTPOSVIS)))
							ans[is.na(FIRSTPOSVIS) | (!is.na(FIRSTPOSVIS) & FIRSTPOSVIS!=min(FIRSTPOSVIS, na.rm=TRUE))]	<- TRUE							
						list(DUMMY=DUMMY, DELETE=ans)		
					}, by=c('RID')], by=c('RID','DUMMY'))
	tmp		<- subset(tmp, !DELETE)
	set(tmp, NULL,c('DUMMY','DELETE'), NULL)
	set(rd, NULL, c('BIRTHDATE','LASTNEGDATE','FIRSTPOSVIS','FIRSTPOSDATE','ARVSTARTDATE','EST_DATEDIED'), NULL)
	rd		<- merge(rd, tmp, by='RID')	
	tmp		<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
	stopifnot(!nrow(merge(subset(tmp[, length(BIRTHDATE), by='RID'], V1>1), tmp, by='RID')))	
	
	#
	#	from every phyloscanner run, select all individuals that had minimum data 
	#	(we need this for sequence sampling bias)
	#rsmpl<- copy(rpng.ok)
	infiles	<- data.table(F=list.files(indir, pattern='pairwise_relationships.rda', full.names=TRUE))
	infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(F)))]
	setkey(infiles, PTY_RUN)
	rsmpl <- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr197_pairwise_relationships.rda'
				load(F)
				list(ID=rplkl[, unique(c(ID1,ID2))])							
			}, by=c('PTY_RUN')]		
	rsmpl	<- unique(subset(rsmpl, select='ID'))
	rsmpl[, MIN_PNG_OUTPUT:=1L]
	setnames(rsmpl, 'ID', 'RID')
	#	add to rsmpl individuals with any sequence data
	load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75_part2.rda')
	tmp		<- unique(subset(dc, !is.na(SID), select=RID))
	tmp[, BAM_OUTPUT:=1L]
	rsmpl	<- merge(tmp, rsmpl, all=1, by='RID')
	set(rsmpl, rsmpl[, which(is.na(MIN_PNG_OUTPUT))],'MIN_PNG_OUTPUT',0L)
	#
	#rsmpl<- unique(subset(rsmpl, HIV==1, c(RID, BAM_OUTPUT, MIN_PNG_OUTPUT)))
	#	add info on all
	load('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/participation_sequencing_age_rccs.rda')
	tmp		<- as.data.table(alltable)
	setnames(tmp, c('RCCS_studyid','hivprev','visit','arvmed','pangea','COMM_NUM'), c('RID','HIV','VISIT','ART','BAM_OUTPUT_2','COMM_NUM_RAW'))
	#	add randomized communities
	tmp2	<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv'))
	tmp		<- merge(tmp, tmp2, by='COMM_NUM_RAW')	
	rsmpl	<- merge(tmp, rsmpl, all.x=1, by='RID')
	set(rsmpl, rsmpl[, which(is.na(MIN_PNG_OUTPUT))],'MIN_PNG_OUTPUT',0L)
	set(rsmpl, rsmpl[, which(is.na(BAM_OUTPUT))],'BAM_OUTPUT',0L)
	# 	treat as cross sectional study
	#	if ever sampled, take earliest
	setkey(rsmpl, RID, VISIT)
	rsmpl	<- rsmpl[, {
				z<- ifelse(length(which.min(HIV==1))==0, 1L, which.min(HIV==1)[1])
				list(	HIV=HIV[z], 
						COMM_NUM=COMM_NUM[z], 
						COMM_NUM_A=COMM_NUM_A[z],
						COMM_TYPE=COMM_TYPE[z],
						AGEYRS=AGEYRS[z], 
						SEX=SEX[z], 
						ART=ART[z], 
						BAM_OUTPUT_2=BAM_OUTPUT_2[z], BAM_OUTPUT=BAM_OUTPUT[z], MIN_PNG_OUTPUT=MIN_PNG_OUTPUT[z])
			}, by=c('RID')]
	setnames(rsmpl, 'BAM_OUTPUT_2', 'HAS_PID')
	rsmpl2	<- as.data.table(seqtable)
	setnames(rsmpl2, 'COMM_NUM', 'COMM_NUM_RAW')
	rsmpl2	<- merge(rsmpl2, tmp2, by='COMM_NUM_RAW')
	rsmpl3	<- as.data.table(partable)	
	setnames(rsmpl3, 'COMM_NUM', 'COMM_NUM_RAW')
	set(rsmpl3, NULL, 'COMM_NUM_RAW', rsmpl3[, as.integer(COMM_NUM_RAW)])
	rsmpl3	<- merge(rsmpl3, tmp2, by='COMM_NUM_RAW')
	setnames(rsmpl3, colnames(rsmpl3), toupper(colnames(rsmpl3)))
	
	save(rsmpl, rsmpl2, rsmpl3, file=outfile)		
}

Rakai.BEEHIVE.compare<- function()
{
	infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/beehive_compare/ORversion_pr.rda'
	load(infile)
	
	#
	#	redo relationships just to be sure we use the latest code version
	trmw.close.brl			<- 0.035
	trmw.distant.brl		<- 0.08
	relationship.types		<- c("TYPE_PAIR_DI","TYPE_PAIRSCORE_DI","TYPE_PAIR_TO","TYPE_PAIR_TODI2x2","TYPE_DIR_TODI3","TYPE_DIRSCORE_TODI3","TYPE_PAIR_TODI2","TYPE_PAIR_TODI","TYPE_PAIRSCORE_TODI", "TYPE_CHAIN_TODI")
	prior.keff				<- 2
	prior.neff				<- 3
	prior.calibrated.prob	<- 0.5
	
	dwin	<- subset(dwin, select=c(SUFFIX, ID1,ID2,TYPE_RAW,PATRISTIC_DISTANCE,ADJACENT,PATHS_12,PATHS_21,ID1_L,ID1_R,ID2_L,ID2_R))	
	dwin	<- phsc.get.basic.pairwise.relationships(dwin, trmw.close.brl, trmw.distant.brl)	
	setnames(dwin, 'TYPE_BASIC', 'TYPE_DIR_TODI7x3')	#for backwards compatibility
	dwin	<- phsc.get.pairwise.relationships(dwin, get.groups=relationship.types, make.pretty.labels=FALSE)
	setnames(dwin, 'TYPE_DIR_TODI7x3', 'TYPE_BASIC')
	set(dwin, NULL, 'W_FROM', dwin[, as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\1', SUFFIX))])
	set(dwin, NULL, 'W_TO', dwin[, as.integer(gsub('[^0-9]*([0-9]+)_to_([0-9]+).*','\\2', SUFFIX))])
	rplkl	<- phsc.get.pairwise.relationships.keff.and.neff(dwin, relationship.types)
	rplkl	<- phsc.get.pairwise.relationships.posterior(rplkl, n.type=prior.keff, n.obs=prior.neff, confidence.cut=prior.calibrated.prob)
	#	make TYBE_BASIC labels nice
	tmp		<- rplkl[, which(GROUP=='TYPE_BASIC')]
	set(rplkl, tmp, 'TYPE', rplkl[tmp, gsub('other_withintermediate_distant','other_distant',gsub('other_withintermediate_close','other_close',gsub('other_withintermediate$','other',gsub('other_nointermediate$','other',gsub('other_nointermediate_distant','other_distant',TYPE)))))])	
	set(rplkl, tmp, 'TYPE', rplkl[tmp, gsub('other_no','other\nno',gsub('([ho])intermediate','\\1 intermediate',gsub('intermediate_','intermediate\n',gsub('intermingled_','intermingled\n',gsub('(chain_[fm][mf])_','\\1\n',gsub('(chain_[12][21])_','\\1\n',TYPE))))))])
	set(rplkl, tmp, 'TYPE', rplkl[tmp, gsub('_',' ',TYPE)])
	set(dwin, NULL, 'TYPE_BASIC', dwin[, gsub('other_withintermediate_distant','other_distant',gsub('other_withintermediate_close','other_close',gsub('other_withintermediate$','other',gsub('other_nointermediate$','other',gsub('other_nointermediate_distant','other_distant',TYPE_BASIC)))))])	
	set(dwin, NULL, 'TYPE_BASIC', dwin[, gsub('other_no','other\nno',gsub('([ho])intermediate','\\1 intermediate',gsub('intermediate_','intermediate\n',gsub('intermingled_','intermingled\n',gsub('(chain_[fm][mf])_','\\1\n',gsub('(chain_[12][21])_','\\1\n',TYPE_BASIC))))))])
	set(dwin, NULL, 'TYPE_BASIC', dwin[, gsub('_',' ',TYPE_BASIC)])
	#	melt dwin
	rpw			<- melt(dwin, variable.name='GROUP', value.name='TYPE', measure.vars=c('TYPE_BASIC',relationship.types))
	set(rpw, NULL, 'ID_R_MAX', rpw[, pmax(ID1_R,ID2_R)])
	set(rpw, NULL, 'ID_R_MIN', rpw[, pmin(ID1_R,ID2_R)])
	
	#
	#	select likely transmitters
	#	ML likely transmission pairs by distance + topology
	rtp		<- subset(rplkl, GROUP=='TYPE_PAIR_TODI2')[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
	rtp		<- subset(rtp, TYPE_MLE=='likely pair')
	rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIR_TODI2' & TYPE=='likely pair'), by=c('ID1','ID2'), all.x=1)			
	rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
	#	save pairs with ML assignment 'likely pair' in rtp.ml
	rtp.ml	<- copy(rtp)
	#	select likely pairs above 50% threshold
	rtp		<- subset(rtp, POSTERIOR_SCORE>0.5)

	#
	#	select likely pairs with resolved direction	
	rtpd	<- merge(subset(rtp, select=c('ID1','ID2')), subset(rplkl, GROUP=='TYPE_DIR_TODI3'), by=c('ID1','ID2'))
	rtpd	<- rtpd[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
	rtpd	<- subset(rtpd, TYPE_MLE!='ambiguous')
	setnames(rtpd, 'TYPE_MLE', 'TYPE')
	rtpd	<- merge(rtpd, subset(rplkl, GROUP=='TYPE_DIRSCORE_TODI3'), by=c('ID1','ID2','TYPE'), all.x=1)			
	rtpd[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
	#	save pairs with ML assignment '12' or '21' in rtpd.ml
	rtpd.ml	<- copy(rtpd)
	#	select likely pairs with direction above 50% threshold
	rtpd	<- subset(rtpd, POSTERIOR_SCORE>0.5)
	
	#
	#	summarize
	cat('\n#likely pairs by maximum likelihood n=', nrow(rtp.ml))							#50
	cat('\n#likely pairs that meet POSTERIOR_SCORE>50% n=', nrow(rtp))					#35
	cat('\n#likely pairs with direction by maximum likelihood n=', nrow(rtpd.ml))			#26
	cat('\n#likely pairs with direction that meet POSTERIOR_SCORE>50% n=', nrow(rtpd))	#20
		
	#
	#	plot all rtp.ml pairs for comparison
	rtpa	<- copy(rtp.ml)
	rtpa[, SELECT:='patients_ambiguous_if_pair_or_not_pair']
	rtpa[, DUMMY:=1:nrow(rtpa)]
	tmp	<- merge(subset(rtpa, select=c(ID1,ID2,DUMMY)),subset(rtp, select=c(ID1,ID2)),by=c('ID1','ID2'))[, DUMMY]
	set(rtpa, tmp, 'SELECT', 'patients_most_likely_a_pair_direction_not_resolved')
	tmp	<- merge(subset(rtpa, select=c(ID1,ID2,DUMMY)),subset(rtpd.ml, select=c(ID1,ID2)),by=c('ID1','ID2'))[, DUMMY]
	set(rtpa, tmp, 'SELECT', 'patients_most_likely_a_pair_direction_not_resolved')
	tmp	<- merge(subset(rtpa, select=c(ID1,ID2,DUMMY)),subset(rtpd, select=c(ID1,ID2)),by=c('ID1','ID2'))[, DUMMY]
	set(rtpa, tmp, 'SELECT', 'patients_most_likely_a_pair_with_resolved_direction')
	
	
	#
	#	define ordered PAIRID
	set(rtpa, NULL, 'SELECT', rtpa[, factor(SELECT, levels=c('patients_ambiguous_if_pair_or_not_pair','patients_most_likely_a_pair_direction_not_resolved','patients_most_likely_a_pair_with_resolved_direction'))])
	rtpa	<- rtpa[order(SELECT,POSTERIOR_SCORE),]
	rtpa[, PAIRID:= rtpa[,factor(DUMMY, levels=rtpa$DUMMY, labels=paste(rtpa$SELECT,rtpa$ID1,rtpa$ID2))]]	
	save(rtpa, rtp.ml, rtp, rtpd.ml, rtpd, rpw, rplkl, file=gsub('\\.rda','_with_selected_pairs.rda',infile))
	
	#
	#	make birds eye view plot
	#	for each pair:
	#	get posterior prob for likely pair
	rpp		<- subset(rplkl, GROUP=='TYPE_PAIR_TODI2' & TYPE=='likely pair')
	#	merge PAIRID
	rpp		<- merge(rpp,unique(subset(rtpa, select=c(ID1,ID2,PAIRID))),by=c('ID1','ID2'))
	#	sort by posterior median highest for likely pair
	rpp		<- rpp[order(qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA)),]
	set(rpp, NULL, 'PAIRID',rpp[, factor(PAIRID, levels=rpp$PAIRID)])	
	p3	<- ggplot(rpp, aes(x=PAIRID)) +
			geom_segment(aes(xend=PAIRID, y=qbeta(0.025, POSTERIOR_ALPHA, POSTERIOR_BETA),yend=qbeta(0.975, POSTERIOR_ALPHA, POSTERIOR_BETA)), colour='grey50') +
			geom_point(aes(y=qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA)), colour='black') +			
			scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.25), labels=scales::percent) +					
			theme_bw() + 			
			#theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
			theme(axis.text.x=element_blank()) +
			theme(axis.ticks.x=element_blank()) +
			theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
			labs(x='',y='',fill='')
	#	for each pair:
	#	make pair topology assignments: ancestral, intermingled, sibling, disconnected
	rplkl2				<- subset(rplkl, GROUP=='TYPE_BASIC')
	rplkl2[, TYPE_TO:= 'disconnected/ not close']
	set(rplkl2, rplkl2[,which(grepl('close', TYPE) & grepl('no intermediate', TYPE) & grepl('chain', TYPE))], 'TYPE_TO', 'close ancestral')
	set(rplkl2, rplkl2[,which(grepl('close', TYPE) & grepl('no intermediate', TYPE) & grepl('intermingled', TYPE))], 'TYPE_TO', 'close intermingled')
	set(rplkl2, rplkl2[,which(grepl('close', TYPE) & grepl('no intermediate', TYPE) & grepl('other', TYPE))], 'TYPE_TO', 'close sibling')
	rplkl2				<- rplkl2[, list(N=N[1],NEFF=NEFF[1],K=sum(K),KEFF=sum(KEFF)), by=c('ID1','ID2','TYPE_TO')]
	set(rplkl2, NULL, 'TYPE_TO', rplkl2[, factor(TYPE_TO, levels=c('close ancestral','close intermingled','close sibling','disconnected/ not close'))])	
	rplkl2				<- merge(rplkl2, unique(subset(rpp, select=c(ID1,ID2,PAIRID))), by=c('ID1','ID2'))
	p2	<- ggplot(rplkl2, aes(x=PAIRID, y=KEFF, fill=TYPE_TO)) +
			geom_bar(stat='identity',position='stack') +
			scale_y_continuous(expand=c(0,0)) +
			scale_fill_manual(values=c("close ancestral"=brewer.pal(11, 'PuOr')[2], "close intermingled"=brewer.pal(11, 'PuOr')[4], 'close sibling'=rev(brewer.pal(11, 'PuOr'))[c(3)], "disconnected/ not close"=rev(brewer.pal(11, 'RdGy'))[4])) +
			labs(x='', y='', fill='') +
			theme_bw() + 
			theme(legend.position='bottom') +
			theme(axis.ticks.x=element_blank()) +
			theme(panel.border=element_blank()) +
			theme(axis.text.x=element_blank()) +
			theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank())
	#	for each pair:
	#	get phyloscanner distances confidence intervals and median values
	rpw2	<- subset(rpw, GROUP=='TYPE_BASIC')[, list(	PD_CL=quantile(PATRISTIC_DISTANCE,p=0.025),
					PD_IL=quantile(PATRISTIC_DISTANCE,p=0.25),
					PD_M=mean(PATRISTIC_DISTANCE),
					PD_IU=quantile(PATRISTIC_DISTANCE,p=0.75),
					PD_CU=quantile(PATRISTIC_DISTANCE,p=0.975)
			),by=c('ID1','ID2')]
	rpw2	<- merge(rpw2, unique(subset(rpp, select=c(ID1,ID2,PAIRID))), by=c('ID1','ID2'))
	set(rpw2, rpw2[, which(PD_M>0.2)], 'PD_M', 0.2)
	p1		<- ggplot(rpw2, aes(x=PAIRID)) +
			geom_segment(aes(xend=PAIRID, y=PD_CL,yend=PD_CU), colour='grey50') +
			geom_point(aes(y=PD_M), colour='black') +
			scale_y_continuous(expand=c(0,0), labels=scales::percent) +
			theme_bw() + 
			theme(axis.text.x=element_blank()) +
			theme(axis.ticks.x=element_blank()) +
			theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
			labs(x='',y='',fill='') +	
			coord_cartesian(ylim=c(0,0.201))
	#	for each pair:
	#	ambiguity measure -- proportion of windows that are not assigned the most likely relationship type
	rplkl2	<- subset(rplkl, GROUP=='TYPE_PAIR_TODI2')
	rpa		<- rplkl2[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]	
	rpa		<- merge(rpa, rplkl2, by=c('ID1','ID2'))
	rpa[, TYPE_AGREE:= factor(TYPE==TYPE_MLE, levels=c(TRUE,FALSE), labels=c('y','n'))]
	rpa		<- rpa[, list(N=N[1], NEFF=NEFF[1], K=sum(K), KEFF=sum(KEFF), PEFF=sum(KEFF)/NEFF[1]), by=c('ID1','ID2','TYPE_AGREE')]
	rpa		<- subset(rpa, TYPE_AGREE=='n')
	rpa		<- merge(rpa, unique(subset(rpp, select=c(ID1,ID2,PAIRID))), by=c('ID1','ID2'))
	p4		<- ggplot(rpa, aes(x=PAIRID)) +			
			geom_point(aes(y=PEFF), colour='black') +
			scale_y_continuous(expand=c(0,0), limits=c(-5e-3,1), labels=scales::percent, breaks=seq(0,1,0.1)) +
			theme_bw() + 
			theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
			#theme(axis.text.x=element_blank()) +
			theme(axis.ticks.x=element_blank()) +
			theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
			labs(x='',y='',fill='') +
			coord_cartesian(ylim=c(-0.002,0.51))	
	plot.file	<- gsub('\\.rda','_with_selected_pairs_birdseyeview.pdf',infile)
	pdf(file=plot.file, h=20, w=30, useDingbats=FALSE)
	grid.newpage()	
	pushViewport(viewport(layout=grid.layout(4, 1, heights=unit(c(3,1.5,1.5,4), "null"))))   	
	print(p1, vp=viewport(layout.pos.row=2, layout.pos.col=1))
	print(p2, vp=viewport(layout.pos.row=1, layout.pos.col=1))
	print(p3, vp=viewport(layout.pos.row=3, layout.pos.col=1))
	print(p4, vp=viewport(layout.pos.row=4, layout.pos.col=1))
	dev.off()	
	
		
	#
	#	plot windows of identified transmission pairs
	for(select in c('patients_ambiguous_if_pair_or_not_pair','patients_most_likely_a_pair_direction_not_resolved','patients_most_likely_a_pair_with_resolved_direction'))
	{
		#select		<- 'couple_ambiguous_if_pair_or_not_pair'
		rps		<- subset(rtpa, SELECT==select, select=c(ID1,ID2,PAIRID))
		setnames(rps, c('ID1','ID2','PAIRID'),c('FEMALE_SANGER_ID','MALE_SANGER_ID','LABEL'))
		rps[, PTY_RUN:=1L]
		for(group in c('TYPE_BASIC','TYPE_PAIR_TODI2','TYPE_DIR_TODI3'))
		{
			plot.file	<- gsub('\\.rda',paste0('_selected_',select,'_windowssummary_',group,'.pdf'),infile)
			rpw2		<- subset(rpw, GROUP==group)
			setnames(rpw2, c('ID1','ID2'),c('FEMALE_SANGER_ID','MALE_SANGER_ID'))
			rpw2[,PTY_RUN:=1L]
			rplkl2		<- subset(rplkl, GROUP==group)
			setnames(rplkl2, c('ID1','ID2'),c('FEMALE_SANGER_ID','MALE_SANGER_ID'))
			rplkl2[,PTY_RUN:=1L]	
			phsc.plot.windowsummaries.for.pairs(rps, rpw2, rplkl2, plot.file, cols=NULL, group=group)			
		}
	}

	
}

RakaiFull.preprocess.trmpairs.todi.phyloscanneroutput.170421<- function()
{
	require(data.table)
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_tb_blnormed'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170421.rda'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl5_blnormed'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170428.rda'
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170428_cl3.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_cl3.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun34'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170610_cl3.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun34d23'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170610_cl3_prior34d23.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun23'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170610_cl3_prior23.rda'	
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min10'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min10.rda'	
	
	
	use.posterior.mode		<- 1
	
	#use.direction.prior.23	<- 1
	#
	#	load couples to search for in phyloscanner output
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_info.rda")
	#
	#	load sequence data
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/circumcision/RCCS_SeqInfo_170505.rda")
	setnames(rs, 'SAMPLE_DATE', 'SEQ_DATE')
	#
	#	load demographic info on all individuals
	tmp		<- RakaiCirc.epi.get.info.170208()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#rn		<- tmp$rn
	ra		<- tmp$ra
	#set(rn, NULL, 'RID', rn[, as.character(RID)])
	#rn		<- merge(rn, subset(rd, select=c(RID, FIRSTPOSVIS, FIRSTPOSDATE)), by='RID', all.x=1)
	#tmp		<- rn[, which(is.na(FIRSTPOSDATE) & !is.na(RECENTVLDATE) & TIMESINCEVL==0)]	#this is dodgy
	#set(rn, tmp, 'FIRSTPOSDATE', rn[tmp, RECENTVLDATE])		
	#rd		<- rbind(rd, rn, use.names=TRUE, fill=TRUE)	#do not consider individuals in the neuro study that are not part of RCCS
	set(rd, NULL, c('PID','SID'), NULL)
	set(rd, NULL, 'SEX', rd[, as.character(SEX)])
	set(rd, NULL, 'RECENTVL', rd[, as.numeric(gsub('< 150','1',gsub('> ','',gsub('BD','',gsub(',','',as.character(RECENTVL))))))])
	set(rd, NULL, 'CAUSE_OF_DEATH', rd[, as.character(CAUSE_OF_DEATH)])
	#	fixup rd: 
	#	remove HIV reverters without sequence
	rd		<- subset(rd, !RID%in%c("C117824","C119303","E118889","K067249"))
	#	fixup complex serology
	set(rd, rd[, which(RID=='B106184')], 'FIRSTPOSDATE', rd[which(RID=='B106184'),DATE])
	set(rd, rd[, which(RID=='B106184')], c('LASTNEGVIS','LASTNEGDATE'), NA_real_)
	set(rd, rd[, which(RID=='B106184')], c('HIVPREV'), 1)
	set(rd, rd[, which(RID=='A008742')], 'FIRSTPOSDATE', rd[which(RID=='A008742'),DATE])
	set(rd, rd[, which(RID=='A008742')], c('HIVPREV'), 1)
	#	fixup rd: 
	#	missing first pos date
	rd		<- subset(rd, RID!='A038432')	#has missing firstposdate and not in PANGEA anyway
	rd		<- subset(rd, RID!='H013226')	#has missing firstposdate and not in PANGEA anyway
	rd		<- subset(rd, RID!='K008173')	#has missing firstposdate and not in PANGEA anyway
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
						if(any(!is.na(FIRSTPOSVIS)))
							ans[is.na(FIRSTPOSVIS) | (!is.na(FIRSTPOSVIS) & FIRSTPOSVIS!=min(FIRSTPOSVIS, na.rm=TRUE))]	<- TRUE							
						list(DUMMY=DUMMY, DELETE=ans)		
					}, by=c('RID')], by=c('RID','DUMMY'))
	tmp		<- subset(tmp, !DELETE)
	set(tmp, NULL,c('DUMMY','DELETE'), NULL)
	set(rd, NULL, c('BIRTHDATE','LASTNEGDATE','FIRSTPOSVIS','FIRSTPOSDATE','ARVSTARTDATE','EST_DATEDIED'), NULL)
	rd		<- merge(rd, tmp, by='RID')	
	tmp		<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
	stopifnot(!nrow(merge(subset(tmp[, length(BIRTHDATE), by='RID'], V1>1), tmp, by='RID')))	
		
	#
	#	from every phyloscanner run, select pairs that are closely related 
	#	
	infiles	<- data.table(F=list.files(indir, pattern='pairwise_relationships.rda', full.names=TRUE))
	infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(F)))]
	setkey(infiles, PTY_RUN)
	rtp.todi2<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun23/ptyr106_pairwise_relationships.rda'
				#cat(PTY_RUN,'\n')
				load(F)
				#	ML likely transmission pairs by distance
				rtp		<- subset(rplkl, GROUP=='TYPE_PAIR_DI')[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, TYPE_MLE=='close')
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIR_DI' & TYPE=='close'), by=c('ID1','ID2'), all.x=1)
				if(!use.posterior.mode)
					rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				if(use.posterior.mode)
					rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]								
				ans		<- copy(rtp)
				#	ML likely transmission pairs by distance + topology
				rtp		<- subset(rplkl, GROUP=='TYPE_PAIR_TODI2')[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, TYPE_MLE=='likely pair')
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIR_TODI2' & TYPE=='likely pair'), by=c('ID1','ID2'), all.x=1)
				if(!use.posterior.mode)
					rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				if(use.posterior.mode)
					rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]
				ans		<- rbind(ans, rtp)
				#	ML directed likely transmission pairs by distance + topology
				#	among those pairs that are likely transmission pairs
				rtp		<- merge(subset(rtp, select=c('ID1','ID2')), subset(rplkl, GROUP=='TYPE_DIR_TODI3'), by=c('ID1','ID2'))
				rtp		<- rtp[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, TYPE_MLE!='ambiguous')
				setnames(rtp, 'TYPE_MLE', 'TYPE')
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_DIRSCORE_TODI3'), by=c('ID1','ID2','TYPE'), all.x=1)
				#if(use.direction.prior.23)
				#{
				#	set(rtp, NULL, c('N_TYPE','PAR_PRIOR','POSTERIOR_ALPHA','POSTERIOR_BETA'), NULL)
				#	rtp	<- phsc.get.pairwise.relationships.posterior(rtp, n.type=2, n.obs=3, confidence.cut=0.66)
				#}					
				if(!use.posterior.mode)
					rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				if(use.posterior.mode)
					rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]				
				rtp[, TYPE_MLE:=TYPE]
				ans		<- rbind(ans, rtp, use.names=TRUE)
				ans				
			}, by=c('PTY_RUN')]		
	set(rtp.todi2, NULL, 'TYPE_MLE', NULL)	
	set(rtp.todi2, NULL, 'GROUP', rtp.todi2[, as.character(GROUP)])
	#
	#	re-arrange to male-female
	#
	#rtp.todi2.o	<- copy(rtp.todi2)
	tmp			<- unique(subset(rd, select=c('RID','SEX')))
	setnames(tmp, colnames(tmp), paste0('ID1_',colnames(tmp)))
	setnames(tmp, c('ID1_RID'), c('ID1'))
	stopifnot( !length(setdiff(rtp.todi2[, ID1], tmp[, ID1])) )
	rtp.todi2	<- merge(rtp.todi2, tmp, by=c('ID1'))	
	setnames(tmp, colnames(tmp), gsub('ID1','ID2',colnames(tmp)))
	stopifnot( !length(setdiff(rtp.todi2[, ID2], tmp[, ID2])) )	
	rtp.todi2	<- merge(rtp.todi2, tmp, by=c('ID2'))
	tmp			<- rtp.todi2[, which(TYPE=='12')]
	set(rtp.todi2, tmp, 'TYPE', rtp.todi2[tmp, tolower(paste0(ID1_SEX,ID2_SEX))])
	tmp			<- rtp.todi2[, which(TYPE=='21')]
	set(rtp.todi2, tmp, 'TYPE', rtp.todi2[tmp, tolower(paste0(ID2_SEX,ID1_SEX))])
	set(rtp.todi2, NULL, 'ID1_SEX', rtp.todi2[, as.character(ID1_SEX)])
	set(rtp.todi2, NULL, 'ID2_SEX', rtp.todi2[, as.character(ID2_SEX)])
	#
	#	prepare all dwin and rplkl and save separately
	#
	rplkl	<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr1_pairwise_relationships.rda'
				#cat(PTY_RUN,'\n')
				load(F)
				rplkl			
			}, by='PTY_RUN']
	rpw		<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr1_pairwise_relationships.rda'
				#cat(PTY_RUN,'\n')
				load(F)
				dwin			
			}, by='PTY_RUN']
	#	add male female
	tmp			<- unique(subset(rd, select=c('RID','SEX')))
	setnames(tmp, colnames(tmp), paste0('ID1_',colnames(tmp)))
	setnames(tmp, c('ID1_RID'), c('ID1'))	
	rplkl		<- merge(rplkl, tmp, by=c('ID1'))
	rpw			<- merge(rpw, tmp, by=c('ID1'))
	setnames(tmp, colnames(tmp), gsub('ID1','ID2',colnames(tmp)))		
	rplkl		<- merge(rplkl, tmp, by=c('ID2'))
	rpw			<- merge(rpw, tmp, by=c('ID2'))
	#	melt rpw
	rpw			<- melt(rpw, variable.name='GROUP', value.name='TYPE', measure.vars=c("TYPE_RAW","TYPE_BASIC","TYPE_DIR_TODI3","TYPE_DIRSCORE_TODI3","TYPE_DIR_TODI4","TYPE_PAIR_TODI2","TYPE_CHAIN_TODI","TYPE_PAIR_DI2","TYPE_PAIR_TO","TYPE_PAIR_TODI2x2"))
	set(rpw, NULL, 'ID_R_MAX', rpw[, pmax(ID1_R,ID2_R)])
	set(rpw, NULL, 'ID_R_MIN', rpw[, pmin(ID1_R,ID2_R)])	
	#	save
	save(rplkl, rpw, file=gsub('\\.rda','_allwindows.rda',outfile))
	#
	#	prepare just the dwin and rplkl that we need for further analysis of the pairs
	#
	tmp			<- unique( subset(rtp.todi2, select=c(ID1, ID2, PTY_RUN)) )	
	infiles		<- subset(infiles, PTY_RUN%in%tmp$PTY_RUN)
	rplkl		<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr1_pairwise_relationships.rda'
				cat(PTY_RUN,'\n')
				load(F)
				ans	<- merge(unique(subset(tmp, select=c('ID1','ID2'))), rplkl, by=c('ID1','ID2'))
				ans			
			}, by='PTY_RUN']
	rpw		<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr1_pairwise_relationships.rda'
				cat(PTY_RUN,'\n')
				load(F)
				ans	<- merge(unique(subset(tmp, select=c('ID1','ID2'))), dwin, by=c('ID1','ID2'))
				ans			
			}, by='PTY_RUN']
	rpw			<- melt(rpw, variable.name='GROUP', value.name='TYPE', measure.vars=c("TYPE_RAW","TYPE_BASIC","TYPE_DIR_TODI3","TYPE_DIRSCORE_TODI3","TYPE_DIR_TODI4","TYPE_PAIR_TODI2","TYPE_CHAIN_TODI","TYPE_PAIR_DI2","TYPE_PAIR_TO","TYPE_PAIR_TODI2x2"))
	set(rpw, NULL, 'ID_R_MAX', rpw[, pmax(ID1_R,ID2_R)])
	set(rpw, NULL, 'ID_R_MIN', rpw[, pmin(ID1_R,ID2_R)])		
	save(rp, rd, rh, ra, rs, rtp.todi2, rplkl, rpw, file=outfile)	
}

RakaiFull.preprocess.trmpairs.todi.phyloscanneroutput.170811<- function()
{
	require(data.table)
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min10'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min10.rda'
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min50'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min50.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p25_stagetwo_rerun23_min30'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl25_prior23_min30.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p45_stagetwo_rerun23_min30'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl45_prior23_min30.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s10_p35_stagetwo_rerun23_min30'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_s10.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s40_p35_stagetwo_rerun23_min30'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_s40.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_d30'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_d30.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_d100'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_d100.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_d1000'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_d1000.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_adj'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_adj.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_rg20'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_rg20.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_rg1'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_rg1.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_prt'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_prt.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_mt'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mt.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_mf6'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mf6.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_mf4'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mf4.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_mf3'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mf3.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_adj'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_adj.rda'		
	indir	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30_zbl'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_zbl.rda'		
	
	#
	#	load couples to search for in phyloscanner output
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_info.rda")
	#
	#	load sequence data
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/circumcision/RCCS_SeqInfo_170505.rda")
	setnames(rs, 'SAMPLE_DATE', 'SEQ_DATE')
	#
	#	load demographic info on all individuals
	tmp		<- RakaiCirc.epi.get.info.170208()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#rn		<- tmp$rn
	ra		<- tmp$ra
	#set(rn, NULL, 'RID', rn[, as.character(RID)])
	#rn		<- merge(rn, subset(rd, select=c(RID, FIRSTPOSVIS, FIRSTPOSDATE)), by='RID', all.x=1)
	#tmp		<- rn[, which(is.na(FIRSTPOSDATE) & !is.na(RECENTVLDATE) & TIMESINCEVL==0)]	#this is dodgy
	#set(rn, tmp, 'FIRSTPOSDATE', rn[tmp, RECENTVLDATE])		
	#rd		<- rbind(rd, rn, use.names=TRUE, fill=TRUE)	#do not consider individuals in the neuro study that are not part of RCCS
	set(rd, NULL, c('PID','SID'), NULL)
	set(rd, NULL, 'SEX', rd[, as.character(SEX)])
	set(rd, NULL, 'RECENTVL', rd[, as.numeric(gsub('< 150','1',gsub('> ','',gsub('BD','',gsub(',','',as.character(RECENTVL))))))])
	set(rd, NULL, 'CAUSE_OF_DEATH', rd[, as.character(CAUSE_OF_DEATH)])
	#	fixup rd: 
	#	remove HIV reverters without sequence
	rd		<- subset(rd, !RID%in%c("C117824","C119303","E118889","K067249"))
	#	fixup complex serology
	set(rd, rd[, which(RID=='B106184')], 'FIRSTPOSDATE', rd[which(RID=='B106184'),DATE])
	set(rd, rd[, which(RID=='B106184')], c('LASTNEGVIS','LASTNEGDATE'), NA_real_)
	set(rd, rd[, which(RID=='B106184')], c('HIVPREV'), 1)
	set(rd, rd[, which(RID=='A008742')], 'FIRSTPOSDATE', rd[which(RID=='A008742'),DATE])
	set(rd, rd[, which(RID=='A008742')], c('HIVPREV'), 1)
	#	fixup rd: 
	#	missing first pos date
	rd		<- subset(rd, RID!='A038432')	#has missing firstposdate and not in PANGEA anyway
	rd		<- subset(rd, RID!='H013226')	#has missing firstposdate and not in PANGEA anyway
	rd		<- subset(rd, RID!='K008173')	#has missing firstposdate and not in PANGEA anyway
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
						if(any(!is.na(FIRSTPOSVIS)))
							ans[is.na(FIRSTPOSVIS) | (!is.na(FIRSTPOSVIS) & FIRSTPOSVIS!=min(FIRSTPOSVIS, na.rm=TRUE))]	<- TRUE							
						list(DUMMY=DUMMY, DELETE=ans)		
					}, by=c('RID')], by=c('RID','DUMMY'))
	tmp		<- subset(tmp, !DELETE)
	set(tmp, NULL,c('DUMMY','DELETE'), NULL)
	set(rd, NULL, c('BIRTHDATE','LASTNEGDATE','FIRSTPOSVIS','FIRSTPOSDATE','ARVSTARTDATE','EST_DATEDIED'), NULL)
	rd		<- merge(rd, tmp, by='RID')	
	tmp		<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
	stopifnot(!nrow(merge(subset(tmp[, length(BIRTHDATE), by='RID'], V1>1), tmp, by='RID')))	
	
	#
	#	from every phyloscanner run, select pairs that are closely related 
	#	
	infiles	<- data.table(F=list.files(indir, pattern='pairwise_relationships.rda', full.names=TRUE))
	infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(F)))]
	setkey(infiles, PTY_RUN)
	rtp.todi2<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun23/ptyr106_pairwise_relationships.rda'
				#cat(PTY_RUN,'\n')
				load(F)
				#	ML likely transmission pairs by distance
				rtp		<- subset(rplkl, GROUP=='TYPE_PAIR_DI')[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, TYPE_MLE=='close')
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIR_DI' & TYPE=='close'), by=c('ID1','ID2'), all.x=1)
				rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]								
				ans		<- copy(rtp)
				#	ML likely transmission pairs by distance + topology
				rtp		<- subset(rplkl, GROUP=='TYPE_PAIR_TODI2')[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, TYPE_MLE=='linked')
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIR_TODI2' & TYPE=='linked'), by=c('ID1','ID2'), all.x=1)
				rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]
				ans		<- rbind(ans, rtp)
				#	ML directed likely transmission pairs by distance + topology
				#	among those pairs that are likely transmission pairs
				rtp		<- merge(subset(rtp, select=c('ID1','ID2')), subset(rplkl, GROUP=='TYPE_NETWORK_SCORES'), by=c('ID1','ID2'))
				rtp		<- rtp[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
				rtp		<- subset(rtp, !TYPE_MLE%in%c('ambiguous','not close/disconnected'))
				setnames(rtp, 'TYPE_MLE', 'TYPE')
				rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_DIR_TODI2'), by=c('ID1','ID2','TYPE'), all.x=1)
				rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]				
				rtp[, TYPE_MLE:=TYPE]
				ans		<- rbind(ans, rtp, use.names=TRUE)
				ans				
			}, by=c('PTY_RUN')]		
	set(rtp.todi2, NULL, 'TYPE_MLE', NULL)	
	set(rtp.todi2, NULL, 'GROUP', rtp.todi2[, as.character(GROUP)])
	#
	#	re-arrange to male-female
	#
	#rtp.todi2.o	<- copy(rtp.todi2)
	tmp			<- unique(subset(rd, select=c('RID','SEX')))
	setnames(tmp, colnames(tmp), paste0('ID1_',colnames(tmp)))
	setnames(tmp, c('ID1_RID'), c('ID1'))
	stopifnot( !length(setdiff(rtp.todi2[, ID1], tmp[, ID1])) )
	rtp.todi2	<- merge(rtp.todi2, tmp, by=c('ID1'))	
	setnames(tmp, colnames(tmp), gsub('ID1','ID2',colnames(tmp)))
	stopifnot( !length(setdiff(rtp.todi2[, ID2], tmp[, ID2])) )	
	rtp.todi2	<- merge(rtp.todi2, tmp, by=c('ID2'))
	tmp			<- rtp.todi2[, which(TYPE=='12')]
	set(rtp.todi2, tmp, 'TYPE', rtp.todi2[tmp, tolower(paste0(ID1_SEX,ID2_SEX))])
	tmp			<- rtp.todi2[, which(TYPE=='21')]
	set(rtp.todi2, tmp, 'TYPE', rtp.todi2[tmp, tolower(paste0(ID2_SEX,ID1_SEX))])
	set(rtp.todi2, NULL, 'ID1_SEX', rtp.todi2[, as.character(ID1_SEX)])
	set(rtp.todi2, NULL, 'ID2_SEX', rtp.todi2[, as.character(ID2_SEX)])
	#
	#	prepare all dwin and rplkl and save separately
	#
	rplkl	<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr1_pairwise_relationships.rda'
				#cat(PTY_RUN,'\n')
				load(F)
				rplkl			
			}, by='PTY_RUN']
	rpw		<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr1_pairwise_relationships.rda'
				#cat(PTY_RUN,'\n')
				load(F)
				dwin			
			}, by='PTY_RUN']
	#	add male female
	tmp			<- unique(subset(rd, select=c('RID','SEX')))
	setnames(tmp, colnames(tmp), paste0('ID1_',colnames(tmp)))
	setnames(tmp, c('ID1_RID'), c('ID1'))	
	rplkl		<- merge(rplkl, tmp, by=c('ID1'))
	rpw			<- merge(rpw, tmp, by=c('ID1'))
	setnames(tmp, colnames(tmp), gsub('ID1','ID2',colnames(tmp)))		
	rplkl		<- merge(rplkl, tmp, by=c('ID2'))
	rpw			<- merge(rpw, tmp, by=c('ID2'))
	#	melt rpw
	rpw			<- melt(rpw, variable.name='GROUP', value.name='TYPE', measure.vars=c("TYPE_RAW","TYPE_BASIC","TYPE_PAIR_DI2","TYPE_PAIR_TO","TYPE_PAIR_TODI2x2","TYPE_PAIR_TODI2","TYPE_DIR_TODI2","TYPE_NETWORK_SCORES","TYPE_CHAIN_TODI"))	
	set(rpw, NULL, 'ID_R_MAX', rpw[, pmax(ID1_R,ID2_R)])
	set(rpw, NULL, 'ID_R_MIN', rpw[, pmin(ID1_R,ID2_R)])	
	#	save
	save(rplkl, rpw, file=gsub('\\.rda','_allwindows.rda',outfile))
	#
	#	prepare just the dwin and rplkl that we need for further analysis of the pairs
	#
	tmp			<- unique( subset(rtp.todi2, select=c(ID1, ID2, PTY_RUN)) )	
	infiles		<- subset(infiles, PTY_RUN%in%tmp$PTY_RUN)
	rplkl		<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr1_pairwise_relationships.rda'
				cat(PTY_RUN,'\n')
				load(F)
				ans	<- merge(unique(subset(tmp, select=c('ID1','ID2'))), rplkl, by=c('ID1','ID2'))
				ans			
			}, by='PTY_RUN']
	rpw		<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s25_resume_sk20_cl3_blnormed/ptyr1_pairwise_relationships.rda'
				cat(PTY_RUN,'\n')
				load(F)
				ans	<- merge(unique(subset(tmp, select=c('ID1','ID2'))), dwin, by=c('ID1','ID2'))
				ans			
			}, by='PTY_RUN']
	rpw			<- melt(rpw, variable.name='GROUP', value.name='TYPE', measure.vars=c("TYPE_RAW","TYPE_BASIC","TYPE_PAIR_DI2","TYPE_PAIR_TO","TYPE_PAIR_TODI2x2","TYPE_PAIR_TODI2","TYPE_DIR_TODI2","TYPE_NETWORK_SCORES","TYPE_CHAIN_TODI"))	
	set(rpw, NULL, 'ID_R_MAX', rpw[, pmax(ID1_R,ID2_R)])
	set(rpw, NULL, 'ID_R_MIN', rpw[, pmin(ID1_R,ID2_R)])		
	save(rp, rd, rh, ra, rs, rtp.todi2, rplkl, rpw, file=outfile)	
}

RakaiFull.plot.epitimeline<- function(df, t.posneg, t.cd4, t.vl, t.seq, age.adult=14, xlim=NA)
{
	
	df[, PAIR:= paste('score',round(POSTERIOR_SCORE, d=2), 'm', MALE_RID, 'f', FEMALE_RID)]
	#	prepare t.posneg for plotting
	t.posneg[, DOA:= DOB+age.adult]
	t.posneg[, DOB:=NULL]
	setnames(t.posneg, colnames(t.posneg), paste0('MALE_',colnames(t.posneg)))
	tmp		<- merge(subset(df, select=c('MALE_RID','FEMALE_RID','PAIR','TYPE')), t.posneg, by='MALE_RID')
	setnames(t.posneg, colnames(t.posneg), gsub('MALE_','FEMALE_',colnames(t.posneg)))
	t.posneg<- merge(tmp, t.posneg, by='FEMALE_RID')
	#	prepare t.seq for plotting
	setnames(t.seq, colnames(t.seq), paste0('MALE_',colnames(t.seq)))
	t.seq.m	<- merge(unique(subset(df, select=c('MALE_RID','PAIR','TYPE'))), t.seq, by='MALE_RID')
	setnames(t.seq, colnames(t.seq), gsub('MALE_','FEMALE_',colnames(t.seq)))
	t.seq.f	<- merge(unique(subset(df, select=c('FEMALE_RID','PAIR','TYPE'))), t.seq, by='FEMALE_RID')
	#	prepare t.cd4 for plotting
	setnames(t.cd4, colnames(t.cd4), paste0('MALE_',colnames(t.cd4)))
	t.cd4.m	<- merge(unique(subset(df, select=c('MALE_RID','PAIR','TYPE'))), t.cd4, by='MALE_RID')
	setnames(t.cd4, colnames(t.cd4), gsub('MALE_','FEMALE_',colnames(t.cd4)))
	t.cd4.f	<- merge(unique(subset(df, select=c('FEMALE_RID','PAIR','TYPE'))), t.cd4, by='FEMALE_RID')
	#	prepare t.vl for plotting
	setnames(t.vl, colnames(t.vl), paste0('MALE_',colnames(t.vl)))
	t.vl.m	<- merge(unique(subset(df, select=c('MALE_RID','PAIR','TYPE'))), t.vl, by='MALE_RID')
	setnames(t.vl, colnames(t.vl), gsub('MALE_','FEMALE_',colnames(t.vl)))
	t.vl.f	<- merge(unique(subset(df, select=c('FEMALE_RID','PAIR','TYPE'))), t.vl, by='FEMALE_RID')
	#	prepare t.posneg for plotting
	tmp		<- t.posneg[, max( max(MALE_FIRSTPOSDATE, na.rm=TRUE), max(FEMALE_FIRSTPOSDATE, na.rm=TRUE))]
	tmp		<- ifelse(any(is.na(xlim)), tmp, max(tmp,max(xlim)))		
	set(t.posneg, NULL, c('MALE_DATE_MAX','FEMALE_DATE_MAX'), tmp)
	tmp		<- t.posneg[, min( min(MALE_FIRSTPOSDATE, na.rm=TRUE), min(FEMALE_FIRSTPOSDATE, na.rm=TRUE), min(MALE_LASTNEGDATE, na.rm=TRUE), min(FEMALE_LASTNEGDATE, na.rm=TRUE))]
	tmp		<- ifelse(any(is.na(xlim)), tmp, min(tmp,min(xlim)))
	set(t.posneg, NULL, c('MALE_DATE_MIN','FEMALE_DATE_MIN'), tmp)
	set(t.posneg, t.posneg[, which(is.na(MALE_LASTNEGDATE))], c('MALE_DOA','MALE_DATE_MIN'), NA_real_)
	set(t.posneg, t.posneg[, which(is.na(FEMALE_LASTNEGDATE))], c('FEMALE_DOA','FEMALE_DATE_MIN'), NA_real_)
	tmp		<- t.posneg[, which(!is.na(MALE_DOA))]
	set(t.posneg, tmp, 'MALE_DATE_MIN', t.posneg[tmp, pmax(MALE_DOA, MALE_DATE_MIN)])
	tmp		<- t.posneg[, which(!is.na(FEMALE_DOA))]
	set(t.posneg, tmp, 'FEMALE_DATE_MIN', t.posneg[tmp, pmax(FEMALE_DOA, FEMALE_DATE_MIN)])	
	tmp		<- t.posneg[, which(!is.na(MALE_DOD))]
	set(t.posneg, tmp, 'MALE_DATE_MAX', t.posneg[tmp, pmin(MALE_DOD, MALE_DATE_MAX)])
	tmp		<- t.posneg[, which(!is.na(FEMALE_DOD))]
	set(t.posneg, tmp, 'FEMALE_DATE_MAX', t.posneg[tmp, pmin(FEMALE_DOD, FEMALE_DATE_MAX)])
	t.posneg[, MALE_EVER_ART:= t.posneg[, !is.na(MALE_ARVSTARTDATE) & MALE_ARVSTARTDATE<MALE_DATE_MAX]]
	t.posneg[, FEMALE_EVER_ART:= t.posneg[, !is.na(FEMALE_ARVSTARTDATE) & FEMALE_ARVSTARTDATE<FEMALE_DATE_MAX]]
	t.posneg[, MALE_ARVSTARTDATE2:= MALE_ARVSTARTDATE]
	t.posneg[, FEMALE_ARVSTARTDATE2:= FEMALE_ARVSTARTDATE]	
	set(t.posneg, t.posneg[, which(!MALE_EVER_ART)], 'MALE_ARVSTARTDATE2', NA_real_)
	set(t.posneg, t.posneg[, which(!FEMALE_EVER_ART)], 'FEMALE_ARVSTARTDATE2', NA_real_)
	t.posneg<- melt(t.posneg, id.vars=c('PAIR','TYPE','MALE_RID','FEMALE_RID','MALE_DOA','FEMALE_DOA','MALE_DOD','FEMALE_DOD','MALE_EVER_ART','FEMALE_EVER_ART'), measure.vars=c('MALE_FIRSTPOSDATE','FEMALE_FIRSTPOSDATE','MALE_ARVSTARTDATE2','FEMALE_ARVSTARTDATE2','FEMALE_LASTNEGDATE','MALE_LASTNEGDATE','MALE_DATE_MIN','FEMALE_DATE_MIN','MALE_DATE_MAX','FEMALE_DATE_MAX','MALE_ARVSTARTDATE','FEMALE_ARVSTARTDATE'))
	t.posneg<- subset(t.posneg, !is.na(value))	
	t.posneg[, DRAW_POINT:= TRUE]
	set(t.posneg, t.posneg[, which(variable=='FEMALE_DATE_MIN' & is.na(FEMALE_DOA))], 'DRAW_POINT', FALSE)
	set(t.posneg, t.posneg[, which(variable=='FEMALE_DATE_MAX' & is.na(FEMALE_DOD))], 'DRAW_POINT', FALSE)
	set(t.posneg, t.posneg[, which(variable=='MALE_DATE_MIN' & is.na(MALE_DOA))], 'DRAW_POINT', FALSE)
	set(t.posneg, t.posneg[, which(variable=='MALE_DATE_MAX' & is.na(MALE_DOD))], 'DRAW_POINT', FALSE)
	set(t.posneg, NULL, c('FEMALE_DOA','FEMALE_DOD','MALE_DOA','MALE_DOD'), NULL)
	set(t.posneg, t.posneg[, which(!MALE_EVER_ART & variable%in%c('MALE_FIRSTPOSDATE','MALE_DATE_MAX'))], 'variable', 'male_pos')
	set(t.posneg, t.posneg[, which(MALE_EVER_ART & variable%in%c('MALE_FIRSTPOSDATE','MALE_ARVSTARTDATE'))], 'variable', 'male_pos')
	set(t.posneg, t.posneg[, which(MALE_EVER_ART & variable%in%c('MALE_ARVSTARTDATE2','MALE_DATE_MAX'))], 'variable', 'male_artstarted')	
	set(t.posneg, t.posneg[, which(!FEMALE_EVER_ART & variable%in%c('FEMALE_FIRSTPOSDATE','FEMALE_DATE_MAX'))], 'variable', 'female_pos')
	set(t.posneg, t.posneg[, which(FEMALE_EVER_ART & variable%in%c('FEMALE_FIRSTPOSDATE','FEMALE_ARVSTARTDATE'))], 'variable', 'female_pos')
	set(t.posneg, t.posneg[, which(FEMALE_EVER_ART & variable%in%c('FEMALE_ARVSTARTDATE2','FEMALE_DATE_MAX'))], 'variable', 'female_artstarted')
	set(t.posneg, t.posneg[, which(variable%in%c('MALE_LASTNEGDATE','MALE_DATE_MIN'))], 'variable', 'male_neg')
	set(t.posneg, t.posneg[, which(variable%in%c('FEMALE_LASTNEGDATE','FEMALE_DATE_MIN'))], 'variable', 'female_neg')
	t.posneg<- subset(t.posneg, grepl('male',variable))
	set(t.posneg, NULL, 'variable', t.posneg[, factor(as.character(variable), levels=c('female_neg','male_neg','male_pos','female_pos','male_artstarted','female_artstarted'))])
	#	plot	
	p	<- ggplot(df) + 
			geom_line(data=subset(t.posneg, grepl('^male',variable) & !is.na(value)), aes(x=value, y=PAIR, colour=variable), position=position_nudge(x=0, y=0.1), size=0.5) + 
			geom_line(data=subset(t.posneg, grepl('^female',variable) & !is.na(value)), aes(x=value, y=PAIR, colour=variable), position=position_nudge(x=0, y=-0.1), size=.5) +
			geom_point(data=subset(t.posneg, grepl('^male',variable) & !is.na(value) & DRAW_POINT), aes(x=value, y=PAIR, colour=variable), position=position_nudge(x=0, y=0.1), size=2, pch=17) + 
			geom_point(data=subset(t.posneg, grepl('^female',variable) & !is.na(value) & DRAW_POINT), aes(x=value, y=PAIR, colour=variable), position=position_nudge(x=0, y=-0.1), size=2, pch=17) +
			geom_point(data=t.vl.m, aes(x=MALE_RECENTVLDATE, y=PAIR, colour=MALE_RECENTVL), position=position_nudge(x=0, y=0.1), pch=1, size=3, stroke=1, fill='transparent') +
			geom_point(data=t.vl.f, aes(x=FEMALE_RECENTVLDATE, y=PAIR, colour=FEMALE_RECENTVL), position=position_nudge(x=0, y=-0.1), pch=1, size=3, stroke=1, fill='transparent') +			
			geom_point(data=t.cd4.m, aes(x=MALE_RECENTCD4DATE, y=PAIR, fill=MALE_RECENTCD4), position=position_nudge(x=0, y=0.1), size=2, pch=23, colour='transparent') +
			geom_point(data=t.cd4.f, aes(x=FEMALE_RECENTCD4DATE, y=PAIR, fill=FEMALE_RECENTCD4), position=position_nudge(x=0, y=-0.1), size=2, pch=23, colour='transparent') +
			geom_point(data=subset(t.seq.m, !is.na(MALE_SEQ_DATE)), aes(x=MALE_SEQ_DATE, y=PAIR), position=position_nudge(x=0, y=0.1), size=2, pch=115) +
			geom_point(data=subset(t.seq.f, !is.na(FEMALE_SEQ_DATE)), aes(x=FEMALE_SEQ_DATE, y=PAIR), position=position_nudge(x=0, y=-0.1), size=2, pch=115) +	# 115 for s
			scale_color_manual(values=c('male_neg'="#5AAE61",'female_neg'="#A6DBA0",'male_pos'="#9970AB",'female_pos'="#C2A5CF",'male_artstarted'="grey50",'female_artstarted'="grey70",
							'<200'="#FEB24C",'200-1,000'="#FD8D3C",'1,000-100,000'="#FC4E2A",'100,000+'="#E31A1C")) +
			scale_fill_manual(values=c('<200'="#DEEBF7", '200-349'="#9ECAE1",'350-499'="#6BAED6",'500-799'="#4292C6",'800+'="#2171B5")) +			
			scale_x_continuous(breaks=seq(1900,2200,2),minor_breaks=seq(1900,2200,1), expand=c(0,0)) +
			coord_cartesian(xlim=t.posneg[, range(value, na.rm=TRUE)]) +
			theme_bw() +
			facet_grid(TYPE~., scales='free_y', space='free_y') +
			labs(x='', y='', colour='', fill='CD4 cells/mm3')
	p
}

RakaiFull.analyze.trmpairs.todi.170421.phylogeography<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170428_withmetadata.rda"		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170428_"	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_withmetadata.rda"		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_"	
	zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	load(infile)	
	#set(rtpdm, NULL, c('MALE_FIRSTPOSVIS.y','MALE_FIRSTPOSDATE.y','FEMALE_FIRSTPOSVIS.y','FEMALE_FIRSTPOSDATE.y'), NULL)
	#setnames(rtpdm, c('MALE_FIRSTPOSVIS.x','MALE_FIRSTPOSDATE.x','FEMALE_FIRSTPOSVIS.x','FEMALE_FIRSTPOSDATE.x'), c('MALE_FIRSTPOSVIS','MALE_FIRSTPOSDATE','FEMALE_FIRSTPOSVIS','FEMALE_FIRSTPOSDATE'))
	nrow(rtpdm)
	#	stage 1: 307 transmissions with direction resolved to 307 recipients with unique transmitter
	#	stage 2: 252 transmissions with unique transmitters	
	
	rtpdm[, AGEDIFF:= rtpdm[, FEMALE_BIRTHDATE-MALE_BIRTHDATE]]
	set(rtpdm, NULL, 'PAIR_ID', rtpdm[, paste0(MALE_RID,'-',FEMALE_RID)])
	set(rtpdm, NULL, 'MALE_SEX', 'M')
	set(rtpdm, NULL, 'FEMALE_SEX', 'F')	
	
	subset(rsm, MIN_PNG_OUTPUT>0)[, quantile(MIN_PNG_OUTPUT/HIV, p=c(0,0.5,1))]
	
	#
	#	some helper data.tables
	rmf		<- subset(rtpdm, TYPE=='mf')
	rfm		<- subset(rtpdm, TYPE=='fm')
	rtr2	<- copy(rmf)
	setnames(rtr2,colnames(rtr2),gsub('FEMALE','REC',colnames(rtr2)))
	setnames(rtr2,colnames(rtr2),gsub('MALE','TR',colnames(rtr2)))
	tmp		<- copy(rfm)
	setnames(tmp,colnames(tmp),gsub('FEMALE','TR',colnames(tmp)))
	setnames(tmp,colnames(tmp),gsub('MALE','REC',colnames(tmp)))
	rtr2	<- rbind(rtr2,tmp)
	
	#
	#	Bayesian source attriution model
	#	adjust for incomplete sampling
	#
	dc	<- rtr2[, list(TR_OBS=length(PAIR_ID)), by=c('TR_COMM_NUM_A','REC_COMM_NUM_A')]	
	#	Bayesian model: add uniform prior
	if(0)
	{
		#	(which is Dirichlet 1 among all communities pairs that have a connection either way
		tmp	<- subset(dc, select=c(REC_COMM_NUM_A, TR_COMM_NUM_A))	
		setnames(tmp, c('REC_COMM_NUM_A','TR_COMM_NUM_A'), c('TR_COMM_NUM_A','REC_COMM_NUM_A'))
		tmp	<- merge(tmp, dc, all.x=1)
		tmp	<- subset(tmp, is.na(TR_OBS))
		set(tmp, NULL, 'TR_OBS', 0)
		dc	<- rbind(dc, tmp)		
	}
	if(1)
	{
		#	This is a bit non-standard, I just don t want the prior to have a large impact, so I chose a sparse one. 
		#	(which is Dirichlet 1 among all communities pairs that are closest)
		#	always add self if not present
		tmp	<- subset(dc[, list(UNOBSERVED_SELF=!any(REC_COMM_NUM_A==TR_COMM_NUM_A)), by='TR_COMM_NUM_A'],UNOBSERVED_SELF, TR_COMM_NUM_A)
		tmp[, REC_COMM_NUM_A:=TR_COMM_NUM_A]
		tmp[, TR_OBS:=0]
		dc	<- rbind(dc, tmp)
		#	ensure each community has at least 1 non-self community, if not add closest other community
		#	I really want to keep this sparse, so do not consider non-self to all communities
		tmp	<- unique(zc, by='COMM_NUM_A')
		tmp	<- as.data.table(t(sapply(seq_len(nrow(tmp)), function(i)
								{
									z<- sort( sqrt( (tmp[,longitude]-tmp[i,longitude])^2+(tmp[,latitude]-tmp[i,latitude])^2 ), index.return=TRUE)$ix
									c('TR_COMM_NUM_A'=tmp[i, COMM_NUM_A], 'REC_COMM_NUM_A'=tmp[z[2],COMM_NUM_A])
								})))	
		z	<- subset(dc[, list(REC_N_OBS=length(REC_COMM_NUM_A)), by='TR_COMM_NUM_A'], REC_N_OBS==1)
		tmp	<- merge(tmp, z, by='TR_COMM_NUM_A')
		tmp	<- merge(subset(tmp, select=c(TR_COMM_NUM_A, REC_COMM_NUM_A)), dc, all.x=1, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A'))
		tmp	<- subset(tmp, is.na(TR_OBS))
		set(tmp, NULL, 'TR_OBS', 0)
		dc	<- rbind(dc, tmp)
	}
	#rsm[, list(ELIGIBLE_AVG=sum(ELIGIBLE_AVG)), by='COMM_TYPE']
	dc[, TR_PRIOR:= 0.5]
	#
	#	Bayesian model first hierarchy: define Beta posterior for sampling probabilities (all alpha and betas)
	#
	tmp	<- subset(rsm, select=c(COMM_NUM_A, ELIGIBLE_AVG, PARTICIPATED_AVG, HIV, MIN_PNG_OUTPUT))
	tmp[, P_PART_EMP:= PARTICIPATED_AVG/ELIGIBLE_AVG]
	tmp[, P_PART_ALPHA:= round(PARTICIPATED_AVG)+1]
	tmp[, P_PART_BETA:= round(ELIGIBLE_AVG-PARTICIPATED_AVG)+1]
	tmp[, P_SEQ_EMP:= MIN_PNG_OUTPUT/HIV]
	tmp[, P_SEQ_ALPHA:= round(MIN_PNG_OUTPUT)+1]
	tmp[, P_SEQ_BETA:= round(HIV-MIN_PNG_OUTPUT)+1]	
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc	<- merge(dc, tmp, by='TR_COMM_NUM_A')
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc	<- merge(dc, tmp, by='REC_COMM_NUM_A')
	#
	#	Bayesian model second hierarchy: draw unobserved data to augment likelihood
	#
	mc.it	<- 1e4
	dcb		<- dc[, {
				tmp	<- 	rbeta(mc.it, TR_P_PART_ALPHA, TR_P_PART_BETA)*
						rbeta(mc.it, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA)*
						rbeta(mc.it, REC_P_PART_ALPHA, REC_P_PART_BETA)*
						rbeta(mc.it, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA)
				#print(tmp)
				tmp	<- rnbinom(mc.it, TR_OBS+TR_PRIOR, tmp)
				#print(tmp)
				list(MONTE_CARLO_IT=seq_len(mc.it), TR_PRIOR=TR_PRIOR, TR_OBS=TR_OBS, TR_MISS= tmp)
			}, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A')]	
	#
	#	Bayesian model second hierarchy: Dirichlet posterior for transmission from community i to j, pi_ij with pi_ij summing to 1
	#
	tmp		<- dcb[, list(	REC_COMM_NUM_A= REC_COMM_NUM_A, 
					TR_COMM_NUM_A= TR_COMM_NUM_A, 
					PI_IJ_ALPHA= TR_OBS+TR_MISS+TR_PRIOR				
			), by='MONTE_CARLO_IT']
	dcb		<- merge(dcb, tmp, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A','MONTE_CARLO_IT'))
	#
	#	this is the end of the source attribution inference on the WAIFM matrix
	#

	#
	#	transmission hubs
	#	proportion of transmitters in community i that have a recipient outside community
	#	this is a summary on the estimated posterior density of the WAIFM matrix
	#	TODO missing community 94
	#
	#	get parameters of posterior under augmented likelihood
	z		<- dcb[, {
				z<- which(REC_COMM_NUM_A==TR_COMM_NUM_A)
				list(P_RECOUTSIDE_ALPHA= sum(PI_IJ_ALPHA[-z]), P_RECOUTSIDE_BETA= sum(PI_IJ_ALPHA[z]))
			}, by=c('TR_COMM_NUM_A','MONTE_CARLO_IT')]	
	#	aggregate and get quantiles
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- quantile(rbeta(length(P_RECOUTSIDE_ALPHA)*mc.it, P_RECOUTSIDE_ALPHA, P_RECOUTSIDE_BETA), p=seq(0,1,0.01))
				list(P=seq(0,1,0.01), Q=unname(tmp))
			}, by='TR_COMM_NUM_A']
	#	subset to main quantities of interest
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_NUM_A~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('REC_P_OUTSIDE_CL','REC_P_OUTSIDE_IL','REC_P_OUTSIDE_M','REC_P_OUTSIDE_IU','REC_P_OUTSIDE_CU'))
	#	add TR_OBS from TR_COMM_NUM_A
	z		<- merge(z, dc[, list(REC_N_OBS=sum(TR_OBS)), by='TR_COMM_NUM_A'], by='TR_COMM_NUM_A')
	setnames(z, 'TR_COMM_NUM_A', 'COMM_NUM_A')	
	#	now after analysis, remove 94 with unknown long / lat
	#z		<- subset(z, COMM_NUM!='94')
	z	<- merge(z, unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE))), by='COMM_NUM_A')
	ggplot(z, aes(x=COMM_NUM_A, middle=REC_P_OUTSIDE_M, min=REC_P_OUTSIDE_CL, max=REC_P_OUTSIDE_CU, lower=REC_P_OUTSIDE_IL, upper=REC_P_OUTSIDE_IU, fill=REC_P_OUTSIDE_M)) + 
			geom_boxplot(outlier.shape=NA, stat='identity') +
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +
			#geom_point(data=subset(z, BS==0), pch=2, colour='red') +
			scale_y_continuous(labels=scales:::percent, breaks=seq(0,1,0.2)) +
			facet_grid(~COMM_TYPE, scales='free_x', space='free_x') +
			labs(x='\ncommunity', y="proportion of transmitters from the community\nthat have a recipient outside the community") +
			theme_bw() + guides(fill='none')
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_adjusted_confidenceintervals_anonmyized.pdf'), w=9, h=5)		
	
	#
	#	transmission hubs - central estimates
	#	among transmitters how many have a recipient outside community
	#	adjusted for sampling
	#	TODO missing community 94
	z[, REC_P_OUTSIDE_PREC:= 1 - REC_P_OUTSIDE_IU + REC_P_OUTSIDE_IL]	
	z	<- merge(unique(subset(zc, select=c(COMM_NUM_A, longitude, latitude, LONG_A, LAT_A)), by='COMM_NUM_A'), z, by='COMM_NUM_A')	
	ggmap(zm) +
			geom_point(data=z, aes(x=longitude, y=latitude, pch=COMM_TYPE, size=REC_N_OBS, fill=100*REC_P_OUTSIDE_M), stroke=1.5, alpha=0.8) +
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +			
			scale_size(breaks=c(5,10,20,40,80), range=c(5,20))+
			scale_shape_manual(values=c('agrarian'=21, 'trading'=22, 'fisherfolk'=23)) +
			#scale_colour_brewer(palette='Dark2') +
			theme(legend.position='bottom', legend.box = "vertical") +			
			labs(	size="transmissions with transmitters from community",
					pch="community type",
					fill="proportion of transmitters from the community\nthat have a recipient outside the community")
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_adjusted.pdf'), w=10, h=10)
	ggmap(zm) +
			geom_point(data=z, aes(x=LONG_A, y=LAT_A, pch=COMM_TYPE, size=REC_N_OBS, fill=100*REC_P_OUTSIDE_M), stroke=1.5, alpha=0.8) +
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +			
			scale_size(breaks=c(5,10,20,40,80), range=c(5,20))+
			scale_shape_manual(values=c('agrarian'=21, 'trading'=22, 'fisherfolk'=23)) +
			#scale_colour_brewer(palette='Dark2') +
			theme(legend.position='bottom', legend.box = "vertical") +			
			labs(	size="transmissions with transmitters from community",
					pch="community type",
					fill="proportion of transmitters from the community\nthat have a recipient outside the community")
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_adjusted_anonymized.pdf'), w=10, h=10)
	
	#
	#	recipient hubs - size is 1-IQR
	ggmap(zm) +
			geom_point(data=z, aes(x=longitude, y=latitude, pch=COMM_TYPE, fill=100*REC_P_OUTSIDE_M, size=REC_P_OUTSIDE_PREC), stroke=1.5, alpha=0.8) +
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +
			#scale_fill_distiller(palette = "RdBu") +
			scale_size(breaks=c(.6,.7,.8,.9,.95), range=c(1,10))+
			scale_shape_manual(values=c('agrarian'=21, 'trading'=22, 'fisherfolk'=23)) +			
			theme(legend.position='bottom', legend.box = "vertical") +
			labs(	size="1 - interquantile range",
					pch="community type",
					fill="proportion of transmitters from the community\nthat have a recipient outside the community")
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_adjusted_sizeIsUnctertainty.pdf'), w=10, h=10)
	ggmap(zm) +
			geom_point(data=z, aes(x=LONG_A, y=LAT_A, pch=COMM_TYPE, fill=100*REC_P_OUTSIDE_M, size=REC_P_OUTSIDE_PREC), stroke=1.5, alpha=0.8) +
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +
			#scale_fill_distiller(palette = "RdBu") +
			scale_size(breaks=c(.6,.7,.8,.9,.95), range=c(1,10))+
			scale_shape_manual(values=c('agrarian'=21, 'trading'=22, 'fisherfolk'=23)) +			
			theme(legend.position='bottom', legend.box = "vertical") +
			labs(	size="1 - interquantile range",
					pch="community type",
					fill="proportion of transmitters from the community\nthat have a recipient outside the community")
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_adjusted_sizeIsUnctertainty_anonymized.pdf'), w=10, h=10)
	
	#
	#	transmission hubs
	#	crude estimates		
	z	<- copy(rtr2)
	z	<- z[, list(REC_OUTSIDE_COMM=length(which(REC_COMM_NUM_A!=TR_COMM_NUM_A)), REC_IN_COMM=length(which(REC_COMM_NUM_A==TR_COMM_NUM_A)) ), by=c('TR_COMM_NUM_A','TR_COMM_TYPE')]
	z[, REC_N:=REC_OUTSIDE_COMM+REC_IN_COMM]
	setnames(z, 'TR_COMM_NUM_A', 'COMM_NUM_A')
	z	<- merge(unique(zc, by='COMM_NUM_A'), z, by='COMM_NUM_A')
	ggmap(zm) +
			geom_point(data=z, aes(x=longitude, y=latitude, pch=COMM_TYPE, colour=COMM_TYPE, size=REC_N, fill=100*REC_OUTSIDE_COMM/REC_N), stroke=1.5, alpha=0.8) + 
			scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=rescale(c(0, .2, 1)), space = "Lab") +
			scale_size(breaks=c(5,10,20,40,80), range=c(5,20))+
			scale_shape_manual(values=c('agrarian'=21, 'trading'=22, 'fisherfolk'=23)) +
			scale_colour_brewer(palette='Dark2') +			
			theme(legend.position='bottom', legend.box = "vertical") +
			labs(	size="transmissions with transmitters from community", 
					fill="proportion of transmitters from the community\nthat have a recipient outside the community")
	ggsave(file=paste0(outfile.base,'_hubs_transmitters_crude.pdf'), w=10, h=10)
	
	
	#
	#	geography transmitter flows into agrarian/trading/fisherolk recipient communities
	#	crude	
	tmp		<- rtr2[,list(N=length(unique(PAIR_ID))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]
	tmp[, P_CELL:= N/sum(N)]
	tmp		<- merge(tmp, tmp[, list(P_REC= N/sum(N), REC_COMM_TYPE=REC_COMM_TYPE), by='TR_COMM_TYPE'], by=c('TR_COMM_TYPE','REC_COMM_TYPE'))
	tmp		<- merge(tmp, tmp[, list(P_TR= N/sum(N), TR_COMM_TYPE=TR_COMM_TYPE), by='REC_COMM_TYPE'], by=c('TR_COMM_TYPE','REC_COMM_TYPE'))
	#tmp[, LABEL:= paste0(N, ' (',round(P_CELL,d=2)*100,'%)\ntransmitters: ',round(P_TR,d=2)*100,'%\nrecipients: ',round(P_REC,d=2)*100,'%')]
	tmp[, LABEL:= paste0(N, '\n(',round(P_TR,d=2)*100,'%)')]
	ggplot(tmp, aes(x=REC_COMM_TYPE, y=P_TR, fill=TR_COMM_TYPE)) + 
			geom_bar(stat='identity', position='dodge') +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + 
			scale_y_continuous(labels=scales::percent, limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.2)) +			
			labs(x='\nlocation likely recipient',y='location likely transmitter\n',fill='transmitter from') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_sources_crude_prop.pdf'), w=6, h=5)
	ggplot(tmp, aes(x=REC_COMM_TYPE, y=N, fill=TR_COMM_TYPE)) + 
			geom_bar(stat='identity', position='dodge') +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + 
			scale_y_continuous() +			
			labs(x='\nlocation likely recipient',y='location likely transmitter\n',fill='transmitter from') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_sources_crude_count.pdf'), w=6, h=5)	
	ggplot(tmp, aes(x=TR_COMM_TYPE, y=P_REC, fill=REC_COMM_TYPE)) + 
			geom_bar(stat='identity', position='dodge') +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + 
			scale_y_continuous(labels=scales::percent, limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.2)) +			
			labs(x='\nlocation likely transmitters',y='location likely recipients\n',fill='recipient in') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_destinations_crude_prop.pdf'), w=6, h=5)
	ggplot(tmp, aes(x=TR_COMM_TYPE, y=N, fill=REC_COMM_TYPE)) + 
			geom_bar(stat='identity', position='dodge') +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + 
			scale_y_continuous() +			
			labs(x='\nlocation likely transmitter',y='location likely recipient\n',fill='recipient in') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_destinations_crude_count.pdf'), w=6, h=5)	
	
	#
	#	geography transmitter flows into agrarian/trading/fisherolk recipient communities
	#	adjusted		
	groups	<- c('agrarian','trading','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- subset(zc, COMM_TYPE==group)$COMM_NUM_A		
				z		<- subset(dcb, REC_COMM_NUM_A%in%tmp)[, list(PI_ITYPE_ALPHA= sum(PI_IJ_ALPHA)), by=c('TR_COMM_NUM_A','MONTE_CARLO_IT')]
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				z		<- z[, list(PI_TYPETYPE_ALPHA=sum(PI_ITYPE_ALPHA)), by=c('TR_COMM_TYPE','MONTE_CARLO_IT')]	
				#	aggregate and get quantiles
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_TYPETYPE_ALPHA)
							colnames(tmp)	<- TR_COMM_TYPE
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				z		<- melt(z, id.vars='MONTE_CARLO_IT', variable.name='TR_COMM_TYPE', value.name='PI_TYPETYPE')
				z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(PI_TYPETYPE, p=seq(0,1,0.01)))), by='TR_COMM_TYPE']
				#	subset to main quantities of interest
				z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_TYPE~P, value.var='Q')
				z[, REC_COMM_TYPE:=group]
				z
			})
	z		<- do.call('rbind',z)
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('TR_CL','TR_IL','TR_M','TR_IU','TR_CU'))
	ggplot(z, aes(x=REC_COMM_TYPE, y=TR_M, ymin=TR_CL, lower=TR_IL, upper=TR_IU, ymax=TR_CU, fill=TR_COMM_TYPE)) + 
			#geom_boxplot(stat='identity', position=position_dodge(width=0.9)) +
			geom_bar(position='dodge', stat='identity') +
			geom_errorbar(width=0.2, position=position_dodge(width=0.9)) +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous(labels=scales::percent, limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.2)) +			
			labs(x='\nlocation likely recipient',y='sources of transmissions\n',fill='transmitter from') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_sources_adjusted.pdf'), w=6, h=5)	
	#
	#	geography transmitter flows from agrarian/trading/fisherolk transmitter communities
	#	adjusted		
	groups	<- c('agrarian','trading','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- subset(zc, COMM_TYPE==group)$COMM_NUM_A		
				z		<- subset(dcb, TR_COMM_NUM_A%in%tmp)[, list(PI_ITYPE_ALPHA= sum(PI_IJ_ALPHA)), by=c('REC_COMM_NUM_A','MONTE_CARLO_IT')]
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(z, tmp, by='REC_COMM_NUM_A')
				z		<- z[, list(PI_TYPETYPE_ALPHA=sum(PI_ITYPE_ALPHA)), by=c('REC_COMM_TYPE','MONTE_CARLO_IT')]	
				#	aggregate and get quantiles
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_TYPETYPE_ALPHA)
							colnames(tmp)	<- REC_COMM_TYPE
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				z		<- melt(z, id.vars='MONTE_CARLO_IT', variable.name='REC_COMM_TYPE', value.name='PI_TYPETYPE')
				z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(PI_TYPETYPE, p=seq(0,1,0.01)))), by='REC_COMM_TYPE']
				#	subset to main quantities of interest
				z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), REC_COMM_TYPE~P, value.var='Q')
				z[, TR_COMM_TYPE:=group]
				z
			})
	z		<- do.call('rbind',z)
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('REC_CL','REC_IL','REC_M','REC_IU','REC_CU'))
	ggplot(z, aes(x=TR_COMM_TYPE, y=REC_M, ymin=REC_CL, lower=REC_IL, upper=REC_IU, ymax=REC_CU, fill=REC_COMM_TYPE)) + 
			#geom_boxplot(stat='identity', position=position_dodge(width=0.9)) +
			geom_bar(position='dodge', stat='identity') +
			geom_errorbar(width=0.2, position=position_dodge(width=0.9)) +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous(labels=scales::percent, limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.2)) +			
			labs(x='\nlocation likely transmitter',y='destination of transmissions\n',fill='recipients in') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_destinations_barplot_adjusted.pdf'), w=6, h=5)	
	


	#
	#	geography flows out > flows in for agrarian/trading/fisherolk communities
	#	adjusted		
	groups	<- c('agrarian','trading','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
				set(tmp, tmp[, which(COMM_TYPE!=group)], 'COMM_TYPE', paste0('non-',group))
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
				setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				z		<- z[, list(PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
				z[, FLOW:=paste0('from ',TR_COMM_TYPE,' to ',REC_COMM_TYPE)]
				#	draw random variables
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
							colnames(tmp)	<- FLOW
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				#	transform as desired
				z[, STAT:= paste0('from ',group,' to non-',group,'\n/\nfrom non-',group,' to ',group)]				
				setnames(z, c(paste0('from ',group,' to non-',group),paste0('from non-',group,' to ',group)), c('A','D'))
				#	get quantiles
				z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(A/D, p=seq(0,1,0.01)))), by='STAT']
			})
	z		<- do.call('rbind',z)
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), STAT~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('FLOW_CL','FLOW_IL','FLOW_M','FLOW_IU','FLOW_CU'))
	z[, COMM_TYPE:= gsub('^from ([a-z]+).*$', '\\1', STAT)]
	ggplot(z, aes(x=COMM_TYPE, middle=FLOW_M, min=FLOW_CL, lower=FLOW_IL, upper=FLOW_IU, max=FLOW_CU, fill=STAT)) +
			geom_hline(yintercept=1, colour='grey50', size=1.5) +
			geom_boxplot(stat='identity') +
			scale_fill_brewer(palette='Dark2') +
			theme_bw() + theme(legend.position='bottom') +
			scale_y_log10(expand=c(0,0), breaks=c(1/4,1/3,1/2,2/3,1,3/2,2,3,4), labels=c('1/4','1/3','1/2','2/3','1','3/2','2','3','4'), limits=c(1/4,5)) +			
			coord_flip() +
			labs(x='community type\n', y='\nflows out / flows in', fill='flow') 
	ggsave(file=paste0(outfile.base,'_phylogeography_aft_flowsoutin_adjusted.pdf'), w=6, h=4)
	
	
	#
	#	geography who infects whom matrix  
	#	crude
	tmp		<- rtr2[,list(N=length(unique(PAIR_ID))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]
	tmp[, P_CELL:= N/sum(N)]
	tmp		<- merge(tmp, tmp[, list(P_REC= N/sum(N), REC_COMM_TYPE=REC_COMM_TYPE), by='TR_COMM_TYPE'], by=c('TR_COMM_TYPE','REC_COMM_TYPE'))
	tmp		<- merge(tmp, tmp[, list(P_TR= N/sum(N), TR_COMM_TYPE=TR_COMM_TYPE), by='REC_COMM_TYPE'], by=c('TR_COMM_TYPE','REC_COMM_TYPE'))
	#tmp[, LABEL:= paste0(N, ' (',round(P_CELL,d=2)*100,'%)\ntransmitters: ',round(P_TR,d=2)*100,'%\nrecipients: ',round(P_REC,d=2)*100,'%')]
	tmp[, LABEL:= paste0(N, '\n(',round(P_TR,d=2)*100,'%)')]
	#tmp[, LABEL:= paste0(N, '\n(',round(P_CELL,d=2)*100,'%)')]
	ggplot(tmp, aes(x=factor(REC_COMM_TYPE),y=factor(TR_COMM_TYPE))) + 
			geom_point(aes(size=N), colour='grey80') +
			geom_text(aes(label=LABEL), nudge_x=0, nudge_y=0, size=3, colour='black') +			
			theme_bw() + 
			scale_size(range = c(5, 50)) +
			labs(x='\nlocation likely recipient',y='location likely transmitter\n') +
			guides(size='none')
	ggsave(file=paste0(outfile.base,'_commtype_3x3.pdf'), w=5, h=5)
	
	set(tmp, NULL, 'REC_COMM_TYPE', tmp[,paste0('to_',REC_COMM_TYPE)])
	set(tmp, NULL, 'TR_COMM_TYPE', tmp[,paste0('from_',TR_COMM_TYPE)])
	tmp		<- suppressWarnings(melt(tmp, id.vars=c('REC_COMM_TYPE','TR_COMM_TYPE'), measure.vars=c('N','P_CELL')))
	tmp		<- dcast.data.table(tmp, variable+TR_COMM_TYPE~REC_COMM_TYPE, value.var='value')
	write.csv(tmp, row.names=FALSE, paste0(outfile.base,'_commtype_3x3_raw.csv'))
	#
	#	geography who infects whom matrix  
	#	adjusted P
	tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	z		<- z[, list(TR_OBS=sum(TR_OBS), TR_MISS=sum(TR_MISS), PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0('from_',TR_COMM_TYPE,' to_',REC_COMM_TYPE)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z[, TR_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\1',variable)]
	z[, REC_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\2',variable)]	
	z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(value, p=seq(0,1,0.01)))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]	
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('PADJ_CL','PADJ_IL','PADJ_M','PADJ_IU','PADJ_CU'))		
	ans		<- melt(z, id.vars=c('TR_COMM_TYPE','REC_COMM_TYPE'), measure.vars=c('PADJ_M','PADJ_CL','PADJ_CU'))
	ans		<- dcast.data.table(ans, variable+REC_COMM_TYPE~TR_COMM_TYPE, value.var='value')
	#	adjusted N
	tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	z		<- z[, list(TR_ADJ=sum(TR_OBS)+sum(TR_MISS)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	set(z, NULL, 'REC_COMM_TYPE', z[,paste0('to_',REC_COMM_TYPE)])
	set(z, NULL, 'TR_COMM_TYPE', z[,paste0('from_',TR_COMM_TYPE)])
	z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(TR_ADJ, p=seq(0,1,0.01), type=1))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]	
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('NADJ_CL','NADJ_IL','NADJ_M','NADJ_IU','NADJ_CU'))
	z		<- melt(z, id.vars=c('TR_COMM_TYPE','REC_COMM_TYPE'), measure.vars=c('NADJ_M','NADJ_CL','NADJ_CU'))
	z		<- dcast.data.table(z, variable+REC_COMM_TYPE~TR_COMM_TYPE, value.var='value')
	ans		<- rbind(z, ans)
	write.csv(ans, row.names=FALSE, paste0(outfile.base,'_commtype_3x3_adjusted.csv'))
	
	
	#
	#	geography transmitters from outside community
	z	<- copy(rtr2)
	z[, FROM_OUTSIDE:= factor(REC_COMM_NUM_A!=TR_COMM_NUM_A, levels=c(TRUE,FALSE), labels=c('outside community','same community'))]
	set(z, NULL, 'REC_COMM_TYPE', z[, factor(REC_COMM_TYPE)])
	ggplot(z, aes(x=REC_COMM_TYPE, fill=FROM_OUTSIDE)) + geom_bar() + 
			theme_bw() + 
			labs(x='\ncommunity type of recipient', fill='transmitter from')
	ggsave(file=paste0(outfile.base,'_extra_community.pdf'), w=5, h=5)
	ggplot(z, aes(x=REC_COMM_TYPE, fill=FROM_OUTSIDE)) + geom_bar(position='fill') + 
			theme_bw() + 
			scale_y_continuous(labels=scales::percent, expand=c(0,0), breaks=seq(0,1,0.2)) +
			labs(x='\ncommunity type of recipient', y='transmitters\n', fill='transmitter from')
	ggsave(file=paste0(outfile.base,'_extra_community_props.pdf'), w=5, h=5)
	#
	#	odds for unequal transmission flows agrarian/fisherfolk
	z	<- z[, list(N=length(REC_RID)), by=c('REC_COMM_TYPE','FROM_OUTSIDE')]
	tmp	<- dcast.data.table(subset(z, REC_COMM_TYPE!='trading'), FROM_OUTSIDE~REC_COMM_TYPE, value.var='N')
	zz	<- as.matrix(tmp[, 2:3, with=FALSE])
	rownames(zz)	<- tmp[[1]]
	fisher.test(zz)
	#STAGE 1
	#data:  zz
	#p-value = 0.0113
	#alternative hypothesis: true odds ratio is not equal to 1
	#95 percent confidence interval:
	#1.137942 3.390227
	#sample estimates:
	#odds ratio 
	#1.963628
	
	#STAGE 2
	#p-value = 0.0002936
	#alternative hypothesis: true odds ratio is not equal to 1
	#95 percent confidence interval:
	#1.698675 7.052786
	#sample estimates:
	#odds ratio 
	#3.427573 
	
	#
	#	geography recipients in different community
	z	<- copy(rtr2)
	z[, FROM_OUTSIDE:= factor(REC_COMM_NUM_A!=TR_COMM_NUM_A, levels=c(TRUE,FALSE), labels=c('outside community','same community'))]
	set(z, NULL, 'TR_COMM_TYPE', z[, factor(TR_COMM_TYPE)])
	ggplot(z, aes(x=TR_COMM_TYPE, fill=FROM_OUTSIDE)) + geom_bar() + 
			theme_bw() + 
			labs(x='\ncommunity type of transmitters', fill='recipient')
	ggsave(file=paste0(outfile.base,'_extra_community_of_recipients.pdf'), w=5, h=5)	
	ggplot(z, aes(x=TR_COMM_TYPE, fill=FROM_OUTSIDE)) + geom_bar(position='fill') + 
			theme_bw() + 
			scale_y_continuous(labels=scales::percent, expand=c(0,0), breaks=seq(0,1,0.2)) +
			labs(x='\ncommunity type of transmitters', y='recipients\n', fill='recipient')
	ggsave(file=paste0(outfile.base,'_extra_community_of_recipients_props.pdf'), w=5, h=5)	
	z	<- z[, list(N=length(TR_RID)), by=c('TR_COMM_TYPE','FROM_OUTSIDE')]
	tmp	<- dcast.data.table(subset(z, TR_COMM_TYPE!='trading'), FROM_OUTSIDE~TR_COMM_TYPE, value.var='N')
	zz	<- as.matrix(tmp[, 2:3, with=FALSE])
	rownames(zz)	<- tmp[[1]]
	fisher.test(zz)
	#STAGE 1
	#p-value = 0.05827
	#alternative hypothesis: true odds ratio is not equal to 1
	#95 percent confidence interval:
	#		0.9429887 3.1200202
	#sample estimates:
	#		odds ratio 
	#1.719531
	
	#STAGE 2
	#p-value = 0.003303
	#alternative hypothesis: true odds ratio is not equal to 1
	#95 percent confidence interval:
	# 1.346652 6.087743
	#sample estimates:
	#odds ratio 
	#  2.844896 
	
	
	#
	#	geography transmitters from outside community by gender of recipients
	z	<- copy(rtr2)
	z[, FROM_OUTSIDE:= factor(REC_COMM_NUM_A!=TR_COMM_NUM_A, levels=c(TRUE,FALSE), labels=c('outside community','same community'))]
	set(z, NULL, 'REC_COMM_TYPE', z[, factor(REC_COMM_TYPE)])
	set(z, NULL, 'REC_SEX', z[, paste0('gender recipient: ',REC_SEX)])	
	ggplot(z, aes(x=REC_COMM_TYPE, fill=FROM_OUTSIDE)) + geom_bar(position='fill') + 
			theme_bw() + 
			scale_y_continuous(labels=scales::percent, expand=c(0,0), breaks=seq(0,1,0.2)) +
			labs(x='\ncommunity type of recipient', y='transmitters\n', fill='transmitter from') +
			facet_grid(~REC_SEX)
	ggsave(file=paste0(outfile.base,'_extra_community_bygender_props.pdf'), w=7, h=5)	
	
	#
	#	community locations
	#		
	zc	<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/community_spatialLoc.csv'))
	set(zc, NULL, 'COMM_NUM_A', zc[,gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',as.character(COMM_NUM)))))])
	zc	<- merge(zc, unique(subset(rh, select=c(COMM_NUM, COMM_TYPE))), by='COMM_NUM')	
	ggmap(zm) +
			geom_point(data=unique(zc, by='COMM_NUM'), aes(x=longitude, y=latitude, pch=COMM_TYPE, colour=COMM_TYPE), size=8, alpha=0.8) +
			geom_text(data=unique(zc, by='COMM_NUM'), aes(x=longitude, y=latitude, label=COMM_NUM), nudge_x=0, nudge_y=0, size=3, colour='black')
	ggsave(file=paste0(outfile.base,'_hubs_comm_locations.pdf'), w=10, h=10)
	
		
}

RakaiFull.analyze.trmpairs.todi.170421.waifwm.community.types<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(gtools)	#rdirichlet
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170428_withmetadata.rda"		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170428_"	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_withmetadata.rda"		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_"	
	zm		<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	zc		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv', stringsAsFactors=FALSE))
	load(infile)	
	#set(rtpdm, NULL, c('MALE_FIRSTPOSVIS.y','MALE_FIRSTPOSDATE.y','FEMALE_FIRSTPOSVIS.y','FEMALE_FIRSTPOSDATE.y'), NULL)
	#setnames(rtpdm, c('MALE_FIRSTPOSVIS.x','MALE_FIRSTPOSDATE.x','FEMALE_FIRSTPOSVIS.x','FEMALE_FIRSTPOSDATE.x'), c('MALE_FIRSTPOSVIS','MALE_FIRSTPOSDATE','FEMALE_FIRSTPOSVIS','FEMALE_FIRSTPOSDATE'))
	nrow(rtpdm)
	#	stage 1: 307 transmissions with direction resolved to 307 recipients with unique transmitter
	#	stage 2: 252 transmissions with unique transmitters	
	
	rtpdm[, AGEDIFF:= rtpdm[, FEMALE_BIRTHDATE-MALE_BIRTHDATE]]
	set(rtpdm, NULL, 'PAIR_ID', rtpdm[, paste0(MALE_RID,'-',FEMALE_RID)])
	set(rtpdm, NULL, 'MALE_SEX', 'M')
	set(rtpdm, NULL, 'FEMALE_SEX', 'F')	
	
	#subset(rsm, MIN_PNG_OUTPUT>0)[, quantile(MIN_PNG_OUTPUT/HIV, p=c(0,0.5,1))]
	
	#
	#	some helper data.tables
	rmf		<- subset(rtpdm, TYPE=='mf')
	rfm		<- subset(rtpdm, TYPE=='fm')
	rtr2	<- copy(rmf)
	setnames(rtr2,colnames(rtr2),gsub('FEMALE','REC',colnames(rtr2)))
	setnames(rtr2,colnames(rtr2),gsub('MALE','TR',colnames(rtr2)))
	tmp		<- copy(rfm)
	setnames(tmp,colnames(tmp),gsub('FEMALE','TR',colnames(tmp)))
	setnames(tmp,colnames(tmp),gsub('MALE','REC',colnames(tmp)))
	rtr2	<- rbind(rtr2,tmp)
	
	#
	#	Bayesian source attriution model
	#	adjust for incomplete sampling
	#
	dc	<- rtr2[, list(TR_OBS=length(PAIR_ID)), by=c('TR_COMM_NUM_A','REC_COMM_NUM_A')]	
	#	Bayesian model: add uniform prior
	if(0)
	{
		#	(which is Dirichlet 1 among all communities pairs that have a connection either way
		tmp	<- subset(dc, select=c(REC_COMM_NUM_A, TR_COMM_NUM_A))	
		setnames(tmp, c('REC_COMM_NUM_A','TR_COMM_NUM_A'), c('TR_COMM_NUM_A','REC_COMM_NUM_A'))
		tmp	<- merge(tmp, dc, all.x=1)
		tmp	<- subset(tmp, is.na(TR_OBS))
		set(tmp, NULL, 'TR_OBS', 0)
		dc	<- rbind(dc, tmp)		
	}
	if(1)
	{
		#	This is a bit non-standard, I just don t want the prior to have a large impact, so I chose a sparse one. 
		#	(which is Dirichlet 1 among all communities pairs that are closest)
		#	always add self if not present
		tmp	<- subset(dc[, list(UNOBSERVED_SELF=!any(REC_COMM_NUM_A==TR_COMM_NUM_A)), by='TR_COMM_NUM_A'],UNOBSERVED_SELF, TR_COMM_NUM_A)
		tmp[, REC_COMM_NUM_A:=TR_COMM_NUM_A]
		tmp[, TR_OBS:=0]
		dc	<- rbind(dc, tmp)
		#	ensure each community has at least 1 non-self community, if not add closest other community
		#	I really want to keep this sparse, so do not consider non-self to all communities
		tmp	<- unique(zc, by='COMM_NUM_A')
		tmp	<- as.data.table(t(sapply(seq_len(nrow(tmp)), function(i)
								{
									z<- sort( sqrt( (tmp[,longitude]-tmp[i,longitude])^2+(tmp[,latitude]-tmp[i,latitude])^2 ), index.return=TRUE)$ix
									c('TR_COMM_NUM_A'=tmp[i, COMM_NUM_A], 'REC_COMM_NUM_A'=tmp[z[2],COMM_NUM_A])
								})))	
		z	<- subset(dc[, list(REC_N_OBS=length(REC_COMM_NUM_A)), by='TR_COMM_NUM_A'], REC_N_OBS==1)
		tmp	<- merge(tmp, z, by='TR_COMM_NUM_A')
		tmp	<- merge(subset(tmp, select=c(TR_COMM_NUM_A, REC_COMM_NUM_A)), dc, all.x=1, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A'))
		tmp	<- subset(tmp, is.na(TR_OBS))
		set(tmp, NULL, 'TR_OBS', 0)
		dc	<- rbind(dc, tmp)
	}
	#rsm[, list(ELIGIBLE_AVG=sum(ELIGIBLE_AVG)), by='COMM_TYPE']
	dc[, TR_PRIOR:= 0.5]
	#
	#	Bayesian model first hierarchy: define Beta posterior for sampling probabilities (all alpha and betas)
	#
	tmp	<- subset(rsm, select=c(COMM_NUM_A, ELIGIBLE_AVG, PARTICIPATED_AVG, HIV, MIN_PNG_OUTPUT))
	tmp[, P_PART_EMP:= PARTICIPATED_AVG/ELIGIBLE_AVG]
	tmp[, P_PART_ALPHA:= round(PARTICIPATED_AVG)+1]
	tmp[, P_PART_BETA:= round(ELIGIBLE_AVG-PARTICIPATED_AVG)+1]
	tmp[, P_SEQ_EMP:= MIN_PNG_OUTPUT/HIV]
	tmp[, P_SEQ_ALPHA:= round(MIN_PNG_OUTPUT)+1]
	tmp[, P_SEQ_BETA:= round(HIV-MIN_PNG_OUTPUT)+1]	
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc	<- merge(dc, tmp, by='TR_COMM_NUM_A')
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc	<- merge(dc, tmp, by='REC_COMM_NUM_A')
	#
	#	Bayesian model second hierarchy: draw unobserved data to augment likelihood
	#
	mc.it	<- 1e4
	dcb		<- dc[, {
				tmp	<- 	rbeta(mc.it, TR_P_PART_ALPHA, TR_P_PART_BETA)*
						rbeta(mc.it, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA)*
						rbeta(mc.it, REC_P_PART_ALPHA, REC_P_PART_BETA)*
						rbeta(mc.it, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA)
				#print(tmp)
				tmp	<- rnbinom(mc.it, TR_OBS+TR_PRIOR, tmp)
				#print(tmp)
				list(MONTE_CARLO_IT=seq_len(mc.it), TR_PRIOR=TR_PRIOR, TR_OBS=TR_OBS, TR_MISS= tmp)
			}, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A')]	
	#
	#	Bayesian model second hierarchy: Dirichlet posterior for transmission from community i to j, pi_ij with pi_ij summing to 1
	#
	tmp		<- dcb[, list(	REC_COMM_NUM_A= REC_COMM_NUM_A, 
					TR_COMM_NUM_A= TR_COMM_NUM_A, 
					PI_IJ_ALPHA= TR_OBS+TR_MISS+TR_PRIOR				
			), by='MONTE_CARLO_IT']
	dcb		<- merge(dcb, tmp, by=c('REC_COMM_NUM_A','TR_COMM_NUM_A','MONTE_CARLO_IT'))
	#
	#	this is the end of the source attribution inference on the WAIFM matrix
	#
	
	
	
	#
	#	geography transmission flows by agrarian/trading/fisherolk 
	#	crude	
	tmp		<- rtr2[,list(N=length(unique(PAIR_ID))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]
	tmp[, P_CELL:= N/sum(N)]
	tmp		<- merge(tmp, tmp[, list(P_REC= N/sum(N), REC_COMM_TYPE=REC_COMM_TYPE), by='TR_COMM_TYPE'], by=c('TR_COMM_TYPE','REC_COMM_TYPE'))
	tmp		<- merge(tmp, tmp[, list(P_TR= N/sum(N), TR_COMM_TYPE=TR_COMM_TYPE), by='REC_COMM_TYPE'], by=c('TR_COMM_TYPE','REC_COMM_TYPE'))
	#	P_REC	is cond prob Prob( recipient in A, T, F | source from x ) ie the entries in WAIFW matrix
	#	P_TR	is cond prob Prob( source from A, T, F | recipient in x )
	#	P_CELL 	is joint prob Prob( source from x, recipient in y )		
	set(tmp, NULL, 'P_TR', tmp[,paste0(round(P_TR,d=3)*100,'%')])
	set(tmp, NULL, 'P_REC', tmp[,paste0(round(P_REC,d=3)*100,'%')])
	set(tmp, NULL, 'P_CELL', tmp[,paste0(round(P_CELL,d=3)*100,'%')])
	set(tmp, NULL, 'REC_COMM_TYPE', tmp[,paste0('to_',REC_COMM_TYPE)])
	set(tmp, NULL, 'TR_COMM_TYPE', tmp[,paste0('from_',TR_COMM_TYPE)])
	setnames(tmp, c('N','P_CELL','P_REC','P_TR'), c('raw_number_phyloscanner_transmissions','raw_proportion_phyloscanner_transmissions','raw_conditionalprob_of_recipients_fixed_source','raw_conditionalprob_of_sources_fixed_recipient'))
	tmp		<- suppressWarnings(melt(tmp, id.vars=c('REC_COMM_TYPE','TR_COMM_TYPE')))
	tmp		<- dcast.data.table(tmp, variable+REC_COMM_TYPE~TR_COMM_TYPE, value.var='value')
	ans		<- copy(tmp)
	
	#
	#	geography joint prob of sources agrarian/trading/fisherolk and recipients agrarian/trading/fisherolk  
	#	adjusted for sequence sampling and participation
	tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	z		<- z[, list(TR_OBS=sum(TR_OBS), TR_MISS=sum(TR_MISS), PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0('from_',TR_COMM_TYPE,' to_',REC_COMM_TYPE)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars='MONTE_CARLO_IT')	
	z[, TR_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\1',variable)]
	z[, REC_COMM_TYPE:= gsub('(from_[a-z]+) (to_[a-z]+)','\\2',variable)]	
	z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(value, p=seq(0,1,0.01)))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]	
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('PADJ_CL','PADJ_IL','PADJ_M','PADJ_IU','PADJ_CU'))
	z[, P_CELL:= z[,paste0(round(PADJ_M,d=3)*100,'% [', round(PADJ_CL,d=3)*100,'%-', round(PADJ_CU,d=3)*100,'%]')]]
	z[, variable:='adjusted-for-par-seq_proportion_phyloscanner_transmissions']
	z		<- dcast.data.table(z, variable+REC_COMM_TYPE~TR_COMM_TYPE, value.var='P_CELL')	
	ans		<- rbind(ans, z)
	
	#	geography % sources agrarian/trading/fisherolk for given recipient communities
	#	adjusted for sequence sampling and participation
	groups	<- c('agrarian','trading','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- subset(zc, COMM_TYPE==group)$COMM_NUM_A		
				z		<- subset(dcb, REC_COMM_NUM_A%in%tmp)[, list(PI_ITYPE_ALPHA= sum(PI_IJ_ALPHA)), by=c('TR_COMM_NUM_A','MONTE_CARLO_IT')]
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				z		<- z[, list(PI_TYPETYPE_ALPHA=sum(PI_ITYPE_ALPHA)), by=c('TR_COMM_TYPE','MONTE_CARLO_IT')]	
				#	aggregate and get quantiles
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_TYPETYPE_ALPHA)
							colnames(tmp)	<- TR_COMM_TYPE
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				z		<- melt(z, id.vars='MONTE_CARLO_IT', variable.name='TR_COMM_TYPE', value.name='PI_TYPETYPE')
				z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(PI_TYPETYPE, p=seq(0,1,0.01)))), by='TR_COMM_TYPE']
				#	subset to main quantities of interest
				z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_TYPE~P, value.var='Q')
				z[, REC_COMM_TYPE:=group]
				z
			})
	z		<- do.call('rbind',z)
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('TR_CL','TR_IL','TR_M','TR_IU','TR_CU'))
	z[, P_TR:= z[,paste0(round(TR_M,d=3)*100,'% [', round(TR_CL,d=3)*100,'%-', round(TR_CU,d=3)*100,'%]')]]
	z[, variable:='adjusted-for-par-seq_conditionalprob_of_sources_fixed_recipient']
	set(z, NULL, 'REC_COMM_TYPE', z[,paste0('to_',REC_COMM_TYPE)])
	set(z, NULL, 'TR_COMM_TYPE', z[,paste0('from_',TR_COMM_TYPE)])
	z		<- dcast.data.table(z, variable+REC_COMM_TYPE~TR_COMM_TYPE, value.var='P_TR')
	ans		<- rbind(ans, z)
	
	#
	#	geography % recipients agrarian/trading/fisherolk for given source communities
	#	adjusted for sequence sampling and participation		
	groups	<- c('agrarian','trading','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- subset(zc, COMM_TYPE==group)$COMM_NUM_A		
				z		<- subset(dcb, TR_COMM_NUM_A%in%tmp)[, list(PI_ITYPE_ALPHA= sum(PI_IJ_ALPHA)), by=c('REC_COMM_NUM_A','MONTE_CARLO_IT')]
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(z, tmp, by='REC_COMM_NUM_A')
				z		<- z[, list(PI_TYPETYPE_ALPHA=sum(PI_ITYPE_ALPHA)), by=c('REC_COMM_TYPE','MONTE_CARLO_IT')]	
				#	aggregate and get quantiles
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_TYPETYPE_ALPHA)
							colnames(tmp)	<- REC_COMM_TYPE
							tmp		<- as.data.table(tmp)								
						}, by=c('MONTE_CARLO_IT')]
				z		<- melt(z, id.vars='MONTE_CARLO_IT', variable.name='REC_COMM_TYPE', value.name='PI_TYPETYPE')
				z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(PI_TYPETYPE, p=seq(0,1,0.01)))), by='REC_COMM_TYPE']
				#	subset to main quantities of interest
				z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), REC_COMM_TYPE~P, value.var='Q')
				z[, TR_COMM_TYPE:=group]
				z
			})
	z		<- do.call('rbind',z)
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('REC_CL','REC_IL','REC_M','REC_IU','REC_CU'))
	z[, P_REC:= z[,paste0(round(REC_M,d=3)*100,'% [', round(REC_CL,d=3)*100,'%-', round(REC_CU,d=3)*100,'%]')]]
	z[, variable:='adjusted-for-par-seq_conditionalprob_of_recipients_fixed_sources']
	set(z, NULL, 'REC_COMM_TYPE', z[,paste0('to_',REC_COMM_TYPE)])
	set(z, NULL, 'TR_COMM_TYPE', z[,paste0('from_',TR_COMM_TYPE)])
	z		<- dcast.data.table(z, variable+REC_COMM_TYPE~TR_COMM_TYPE, value.var='P_REC')
	ans		<- rbind(ans, z)
	
	
	suradj	<- rsm[, list(RCCS_ELIGIBLE_AVG=round(sum(ELIGIBLE_AVG))), by='COMM_TYPE']
	suradj	<- merge(suradj, data.table(COMM_TYPE=c('agrarian','trading','fisherfolk'), RAKAI_BEST_GUESS=round(c(271814.29,29425.87,19989.65))), by='COMM_TYPE')
	suradj[, RAKAI_BEST_GUESS_P:= RAKAI_BEST_GUESS/sum(RAKAI_BEST_GUESS)]  
	suradj[, RCCS_ELIGIBLE_AVG_P:= RCCS_ELIGIBLE_AVG/sum(RCCS_ELIGIBLE_AVG)]
	suradj[, BEST_GUESS_ADJ:= RAKAI_BEST_GUESS_P/RCCS_ELIGIBLE_AVG_P]
	
	#
	#	geography joint prob of sources agrarian/trading/fisherolk and recipients agrarian/trading/fisherolk  
	#	adjusted for sequence sampling and participation
	tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
	setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
	z		<- merge(dcb, tmp, by='REC_COMM_NUM_A')
	setnames(tmp, c('REC_COMM_NUM_A','REC_COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
	z		<- merge(z, tmp, by='TR_COMM_NUM_A')
	z		<- z[, list(TR_OBS=sum(TR_OBS), TR_MISS=sum(TR_MISS), PI_ST_ALPHA=sum(PI_IJ_ALPHA)), by=c('REC_COMM_TYPE','TR_COMM_TYPE','MONTE_CARLO_IT')]
	z[, FLOW:=paste0('from_',TR_COMM_TYPE,' to_',REC_COMM_TYPE)]
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_ST_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)	
				tmp[, MONTE_CARLO_IT2:= seq_len(mc.it)]
			}, by=c('MONTE_CARLO_IT')]
	z		<- melt(z, id.vars=c('MONTE_CARLO_IT','MONTE_CARLO_IT2'))	
	z[, TR_COMM_TYPE:= gsub('from_([a-z]+) to_([a-z]+)','\\1',variable)]
	z[, REC_COMM_TYPE:= gsub('from_([a-z]+) to_([a-z]+)','\\2',variable)]
	#	adjust for community selection
	tmp		<- subset(suradj, select=c(COMM_TYPE, BEST_GUESS_ADJ))
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))	
	z		<- merge(z, tmp, by='TR_COMM_TYPE')
	tmp		<- subset(suradj, select=c(COMM_TYPE, BEST_GUESS_ADJ))
	setnames(tmp, colnames(tmp), paste0('REC_',colnames(tmp)))	
	z		<- merge(z, tmp, by='REC_COMM_TYPE')	
	z		<- z[, list(TR_COMM_TYPE=TR_COMM_TYPE, REC_COMM_TYPE=REC_COMM_TYPE, value=value, value_ADJ=value*TR_BEST_GUESS_ADJ*REC_BEST_GUESS_ADJ/sum(value*TR_BEST_GUESS_ADJ*REC_BEST_GUESS_ADJ)), by=c('MONTE_CARLO_IT','MONTE_CARLO_IT2')]
	
	z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(value_ADJ, p=seq(0,1,0.01)))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]	
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_TYPE+REC_COMM_TYPE~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('PADJ_CL','PADJ_IL','PADJ_M','PADJ_IU','PADJ_CU'))
	z[, P_CELL:= z[,paste0(round(PADJ_M,d=3)*100,'%')]]
	z[, variable:='adjusted-for-bestguesssur-par-seq_proportion_phyloscanner_transmissions']
	set(z, NULL, 'REC_COMM_TYPE', z[,paste0('to_',REC_COMM_TYPE)])
	set(z, NULL, 'TR_COMM_TYPE', z[,paste0('from_',TR_COMM_TYPE)])	
	z		<- dcast.data.table(z, variable+REC_COMM_TYPE~TR_COMM_TYPE, value.var='P_CELL')	
	ans		<- rbind(ans, z)
	
	
	#	geography % sources agrarian/trading/fisherolk for given recipient communities
	#	adjusted for survey community selection and sequence sampling and participation
	groups	<- c('agrarian','trading','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- subset(zc, COMM_TYPE==group)$COMM_NUM_A		
				z		<- subset(dcb, REC_COMM_NUM_A%in%tmp)[, list(PI_ITYPE_ALPHA= sum(PI_IJ_ALPHA)), by=c('TR_COMM_NUM_A','MONTE_CARLO_IT')]
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('TR_COMM_NUM_A','TR_COMM_TYPE'))
				z		<- merge(z, tmp, by='TR_COMM_NUM_A')
				z		<- z[, list(PI_TYPETYPE_ALPHA=sum(PI_ITYPE_ALPHA)), by=c('TR_COMM_TYPE','MONTE_CARLO_IT')]	
				#	aggregate and get quantiles
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_TYPETYPE_ALPHA)
							colnames(tmp)	<- TR_COMM_TYPE
							tmp		<- as.data.table(tmp)
							tmp[, MONTE_CARLO_IT2:= seq_len(mc.it)]
						}, by=c('MONTE_CARLO_IT')]
				z		<- melt(z, id.vars=c('MONTE_CARLO_IT','MONTE_CARLO_IT2'), variable.name='COMM_TYPE', value.name='PI_TYPETYPE')
				#	adjust for community selection -- we choose not to account for uncertainty, only adjust for means
				z		<- merge(z, subset(suradj, select=c(COMM_TYPE, BEST_GUESS_ADJ)), by='COMM_TYPE')
				z		<- z[, list(COMM_TYPE=COMM_TYPE, PI_TYPETYPE=PI_TYPETYPE, PI_TYPETYPE_ADJ=PI_TYPETYPE*BEST_GUESS_ADJ/sum(PI_TYPETYPE*BEST_GUESS_ADJ)), by=c('MONTE_CARLO_IT','MONTE_CARLO_IT2')]				
				z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(PI_TYPETYPE_ADJ, p=seq(0,1,0.01)))), by='COMM_TYPE']
				setnames(z, 'COMM_TYPE', 'TR_COMM_TYPE')
				#	subset to main quantities of interest
				z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_COMM_TYPE~P, value.var='Q')
				z[, REC_COMM_TYPE:=group]
				z
			})
	z		<- do.call('rbind',z)
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('TR_CL','TR_IL','TR_M','TR_IU','TR_CU'))
	z[, P_TR:= z[,paste0(round(TR_M,d=3)*100,'%')]]
	z[, variable:='adjusted-for-bestguesssur-par-seq_conditionalprob_of_sources_fixed_recipient']
	set(z, NULL, 'REC_COMM_TYPE', z[,paste0('to_',REC_COMM_TYPE)])
	set(z, NULL, 'TR_COMM_TYPE', z[,paste0('from_',TR_COMM_TYPE)])
	z		<- dcast.data.table(z, variable+REC_COMM_TYPE~TR_COMM_TYPE, value.var='P_TR')
	ans		<- rbind(ans, z)
	
	#
	#	geography % recipients agrarian/trading/fisherolk for given source communities
	#	adjusted for survey community selection and sequence sampling and participation		
	groups	<- c('agrarian','trading','fisherfolk')
	z		<- lapply(groups, function(group)
			{				
				tmp		<- subset(zc, COMM_TYPE==group)$COMM_NUM_A		
				z		<- subset(dcb, TR_COMM_NUM_A%in%tmp)[, list(PI_ITYPE_ALPHA= sum(PI_IJ_ALPHA)), by=c('REC_COMM_NUM_A','MONTE_CARLO_IT')]
				tmp		<- unique(subset(zc, select=c(COMM_NUM_A, COMM_TYPE)))
				setnames(tmp, c('COMM_NUM_A','COMM_TYPE'), c('REC_COMM_NUM_A','REC_COMM_TYPE'))
				z		<- merge(z, tmp, by='REC_COMM_NUM_A')
				z		<- z[, list(PI_TYPETYPE_ALPHA=sum(PI_ITYPE_ALPHA)), by=c('REC_COMM_TYPE','MONTE_CARLO_IT')]	
				#	aggregate and get quantiles
				mc.it	<- 1e2
				z		<- z[, {												
							tmp		<- rdirichlet(mc.it, PI_TYPETYPE_ALPHA)
							colnames(tmp)	<- REC_COMM_TYPE
							tmp		<- as.data.table(tmp)
							tmp[, MONTE_CARLO_IT2:= seq_len(mc.it)]
						}, by=c('MONTE_CARLO_IT')]
				z		<- melt(z, id.vars=c('MONTE_CARLO_IT','MONTE_CARLO_IT2'), variable.name='COMM_TYPE', value.name='PI_TYPETYPE')
				#	adjust for community selection -- we choose not to account for uncertainty, only adjust for means
				z		<- merge(z, subset(suradj, select=c(COMM_TYPE, BEST_GUESS_ADJ)), by='COMM_TYPE')
				z		<- z[, list(COMM_TYPE=COMM_TYPE, PI_TYPETYPE=PI_TYPETYPE, PI_TYPETYPE_ADJ=PI_TYPETYPE*BEST_GUESS_ADJ/sum(PI_TYPETYPE*BEST_GUESS_ADJ)), by=c('MONTE_CARLO_IT','MONTE_CARLO_IT2')]				
				z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(PI_TYPETYPE_ADJ, p=seq(0,1,0.01)))), by='COMM_TYPE']
				setnames(z, 'COMM_TYPE', 'REC_COMM_TYPE')
				#	subset to main quantities of interest
				z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), REC_COMM_TYPE~P, value.var='Q')
				z[, TR_COMM_TYPE:=group]
				z
			})
	z		<- do.call('rbind',z)
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('REC_CL','REC_IL','REC_M','REC_IU','REC_CU'))
	z[, P_REC:= z[,paste0(round(REC_M,d=3)*100,'%')]]
	z[, variable:='adjusted-for-bestguesssur-par-seq_conditionalprob_of_recipients_fixed_sources']
	set(z, NULL, 'REC_COMM_TYPE', z[,paste0('to_',REC_COMM_TYPE)])
	set(z, NULL, 'TR_COMM_TYPE', z[,paste0('from_',TR_COMM_TYPE)])
	z		<- dcast.data.table(z, variable+REC_COMM_TYPE~TR_COMM_TYPE, value.var='P_REC')
	ans		<- rbind(ans, z)
	
	write.csv(ans, file=paste0(outfile.base, 'WAIFW_communitytypes.csv'))
}

RakaiFull.analyze.trmpairs.todi.170421.age<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ggmap)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(viridis)
	
	g_legend<-function(a.gplot)
	{
		tmp <- ggplot_gtable(ggplot_build(a.gplot))
		leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
		legend <- tmp$grobs[[leg]]
		legend
	}
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170428_withmetadata.rda"		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170428_"	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_withmetadata.rda"		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_"	
	
	load(infile)	
	
	nrow(rtpdm)
	#	stage 1: 307 transmissions with direction resolved to 307 recipients with unique transmitter
	#	stage 2: 252 transmissions with unique transmitters	
	
	rtpdm[, AGEDIFF:= rtpdm[, FEMALE_BIRTHDATE-MALE_BIRTHDATE]]
	set(rtpdm, NULL, 'PAIR_ID', rtpdm[, paste0(MALE_RID,'-',FEMALE_RID)])
	set(rtpdm, NULL, 'MALE_SEX', 'M')
	set(rtpdm, NULL, 'FEMALE_SEX', 'F')	
	rtpdm[, MALE_AGE_AT_CONCPOS:= pmax(MALE_FIRSTPOSDATE,FEMALE_FIRSTPOSDATE)-MALE_BIRTHDATE]
	rtpdm[, FEMALE_AGE_AT_CONCPOS:= pmax(MALE_FIRSTPOSDATE,FEMALE_FIRSTPOSDATE)-FEMALE_BIRTHDATE]
	rtpdm[, MALE_AGE_AT_CONCPOS_C:= rtpdm[, as.character(cut(MALE_AGE_AT_CONCPOS, right=FALSE, breaks=c(seq(15,45,5),65), labels=paste0(seq(15,50,5)[-8], '-', seq(15,50,5)[-1]-1)))]]
	rtpdm[, FEMALE_AGE_AT_CONCPOS_C:= rtpdm[, as.character(cut(FEMALE_AGE_AT_CONCPOS, right=FALSE, breaks=c(seq(15,45,5),65), labels=paste0(seq(15,50,5)[-8], '-', seq(15,50,5)[-1]-1)))]]
	cat('\nno birthdate for', paste(c(subset(rtpdm, is.na(MALE_BIRTHDATE))[, MALE_RID], subset(rtpdm, is.na(FEMALE_BIRTHDATE))[, FEMALE_RID]), collapse=', '))
	#subset(rd, RID%in%c(subset(rtpdm, is.na(MALE_BIRTHDATE))[, MALE_RID], subset(rtpdm, is.na(FEMALE_BIRTHDATE))[, FEMALE_RID]))	
	rtpdm	<- subset(rtpdm, !is.na(MALE_BIRTHDATE) & !is.na(FEMALE_BIRTHDATE))	
	rsma	<- rsma[, list(PARTICIPATED=sum(PARTICIPATED), HIV=sum(HIV), HAS_PID=sum(HAS_PID), BAM_OUTPUT=sum(BAM_OUTPUT), MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT)), by=c('SEX','AGE_C')]
	set(rsma, NULL, 'AGE_C', rsma[, as.character(AGE_C)])
	#
	#	some helper data.tables
	rmf		<- subset(rtpdm, TYPE=='mf')
	rfm		<- subset(rtpdm, TYPE=='fm')
	rtr2	<- copy(rmf)
	setnames(rtr2,colnames(rtr2),gsub('FEMALE','REC',colnames(rtr2)))
	setnames(rtr2,colnames(rtr2),gsub('MALE','TR',colnames(rtr2)))
	tmp		<- copy(rfm)
	setnames(tmp,colnames(tmp),gsub('FEMALE','TR',colnames(tmp)))
	setnames(tmp,colnames(tmp),gsub('MALE','REC',colnames(tmp)))
	rtr2	<- rbind(rtr2,tmp)
	#	TODO double check when these 9 with missing birth dates were born
	write.csv(subset(rtpdm, AGEDIFF>15), file=paste0(outfile.base,'check_age_difference_g_15_years.csv'))
	#
	#	adjust for incomplete sampling
	#	Bayesian model
	dc	<- rtr2[, list(TR_OBS=length(PAIR_ID)), by=c('TR_SEX','REC_SEX','REC_AGE_AT_CONCPOS_C','TR_AGE_AT_CONCPOS_C')]	
	#	Bayesian model: add uniform prior
	if(1)
	{
		#	which is Dirichlet among all possible pairs
		tmp	<- subset(as.data.table(expand.grid(TR_SEX=dc[, sort(unique(TR_SEX))], REC_SEX=dc[, sort(unique(REC_SEX))], REC_AGE_AT_CONCPOS_C=dc[, sort(unique(REC_AGE_AT_CONCPOS_C))], TR_AGE_AT_CONCPOS_C=dc[, sort(unique(TR_AGE_AT_CONCPOS_C))])), TR_SEX!=REC_SEX)
		dc	<- merge(tmp, dc, all.x=1)
		set(dc, dc[, which(is.na(TR_OBS))], 'TR_OBS', 0)
	}
	if(0)
	{
		#	(which is Dirichlet 1 among all communities pairs that have a connection either way
		tmp	<- subset(dc, select=c(REC_COMM_NUM, TR_COMM_NUM))	
		setnames(tmp, c('REC_COMM_NUM','TR_COMM_NUM'), c('TR_COMM_NUM','REC_COMM_NUM'))
		tmp	<- merge(tmp, dc, all.x=1)
		tmp	<- subset(tmp, is.na(TR_OBS))
		set(tmp, NULL, 'TR_OBS', 0)
		dc	<- rbind(dc, tmp)		
	}
	if(0)
	{
		#	(which is Dirichlet 1 among all communities pairs that are closest)
		#	always add self if not present
		tmp	<- subset(dc[, list(UNOBSERVED_SELF=!any(REC_COMM_NUM==TR_COMM_NUM)), by='TR_COMM_NUM'],UNOBSERVED_SELF, TR_COMM_NUM)
		tmp[, REC_COMM_NUM:=TR_COMM_NUM]
		tmp[, TR_OBS:=0]
		dc	<- rbind(dc, tmp)
		#	ensure each community has at least 1 non-self community, if not add closest other community
		#	I really want to keep this sparse, so do not consider non-self to all communities
		tmp	<- unique(zc, by='COMM_NUM')
		tmp	<- as.data.table(t(sapply(seq_len(nrow(tmp)), function(i)
								{
									z<- sort( sqrt( (tmp[,longitude]-tmp[i,longitude])^2+(tmp[,latitude]-tmp[i,latitude])^2 ), index.return=TRUE)$ix
									c('TR_COMM_NUM'=tmp[i, COMM_NUM], 'REC_COMM_NUM'=tmp[z[2],COMM_NUM])
								})))	
		z	<- subset(dc[, list(REC_N_OBS=length(REC_COMM_NUM)), by='TR_COMM_NUM'], REC_N_OBS==1)
		tmp	<- merge(tmp, z, by='TR_COMM_NUM')
		tmp	<- merge(subset(tmp, select=c(TR_COMM_NUM, REC_COMM_NUM)), dc, all.x=1, by=c('REC_COMM_NUM','TR_COMM_NUM'))
		tmp	<- subset(tmp, is.na(TR_OBS))
		set(tmp, NULL, 'TR_OBS', 0)
		dc	<- rbind(dc, tmp)
	}
	dc[, TR_PRIOR:= 0.5]
	#	Bayesian model first hierarchy: define Beta posterior for sampling probabilities (all alpha and betas)
	tmp	<- subset(rsma, select=c(SEX, AGE_C, HIV, MIN_PNG_OUTPUT))
	tmp[, P_SEQ_EMP:= MIN_PNG_OUTPUT/HIV]
	tmp[, P_SEQ_ALPHA:= round(MIN_PNG_OUTPUT)+1]
	tmp[, P_SEQ_BETA:= round(HIV-MIN_PNG_OUTPUT)+1]
	setnames(tmp, 'AGE_C', 'AGE_AT_CONCPOS_C')
	setnames(tmp, colnames(tmp), paste0('TR_',colnames(tmp)))
	dc	<- merge(dc, tmp, by=c('TR_SEX','TR_AGE_AT_CONCPOS_C'))
	setnames(tmp, colnames(tmp), gsub('TR_','REC_',colnames(tmp)))
	dc	<- merge(dc, tmp, by=c('REC_SEX','REC_AGE_AT_CONCPOS_C'))
	setkey(dc, TR_SEX, TR_AGE_AT_CONCPOS_C, REC_AGE_AT_CONCPOS_C)
	#	Bayesian model second hierarchy: draw unobserved data to augment likelihood
	mc.it	<- 1e4
	dcb		<- dc[, {
				tmp	<- 	rbeta(mc.it, TR_P_SEQ_ALPHA, TR_P_SEQ_BETA)*						
						rbeta(mc.it, REC_P_SEQ_ALPHA, REC_P_SEQ_BETA)
				#print(tmp)
				tmp	<- rnbinom(mc.it, TR_OBS+TR_PRIOR, tmp)
				#print(tmp)
				list(MONTE_CARLO_IT=seq_len(mc.it), TR_PRIOR=TR_PRIOR, TR_OBS=TR_OBS, TR_MISS=tmp)
			}, by=c('TR_SEX','REC_SEX','REC_AGE_AT_CONCPOS_C','TR_AGE_AT_CONCPOS_C')]	
	#	Bayesian model second hierarchy: Dirichlet posterior for transmission from community i to j, pi_ij with pi_ij summing to 1
	dcb[, PI_IJ_ALPHA:= TR_OBS+TR_MISS+TR_PRIOR]	
	#
	#	age-age matrix
	#	
	z		<- subset(dcb, select=c(TR_SEX, REC_AGE_AT_CONCPOS_C, TR_AGE_AT_CONCPOS_C, MONTE_CARLO_IT, PI_IJ_ALPHA))
	z[, FLOW:=paste0(TR_AGE_AT_CONCPOS_C,'_to_',REC_AGE_AT_CONCPOS_C)]
	set(z, NULL, c('REC_AGE_AT_CONCPOS_C','TR_AGE_AT_CONCPOS_C'), NULL)
	setkey(z, TR_SEX, MONTE_CARLO_IT, FLOW)
	mc.it	<- 1e2
	tmp		<- subset(z, TR_SEX=='M')[, {												
				tmp		<- rdirichlet(mc.it, PI_IJ_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('TR_SEX','MONTE_CARLO_IT')]
	tmp		<- melt(tmp, id.vars=c('TR_SEX','MONTE_CARLO_IT'), variable.name='FLOW', value.name='PI_IJ')
	set(tmp, NULL, 'FLOW', tmp[, as.character(FLOW)])
	z		<- subset(z, TR_SEX=='F')[, {												
				tmp		<- rdirichlet(mc.it, PI_IJ_ALPHA)
				colnames(tmp)	<- FLOW
				tmp		<- as.data.table(tmp)								
			}, by=c('TR_SEX','MONTE_CARLO_IT')]
	z		<- melt(z, id.vars=c('TR_SEX','MONTE_CARLO_IT'), variable.name='FLOW', value.name='PI_IJ')
	set(z, NULL, 'FLOW', z[, as.character(FLOW)])
	z		<- rbind(z, tmp)
	z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(PI_IJ, p=seq(0,1,0.01)))), by=c('TR_SEX','FLOW')]	
	z[, TR_AGE_AT_CONCPOS_C:= z[, gsub('^(.*)_to_(.*)$','\\1',FLOW)]]
	z[, REC_AGE_AT_CONCPOS_C:= z[, gsub('^(.*)_to_(.*)$','\\2',FLOW)]]		
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_SEX+TR_AGE_AT_CONCPOS_C+REC_AGE_AT_CONCPOS_C~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('CL','IL','M','IU','CU'))
	z[, REC_SEX:= gsub('F','male likely recipients',gsub('M','female likely recipients',TR_SEX))]
	set(z, NULL, 'TR_SEX', z[,gsub('M','male likely transmitters',gsub('F','female likely transmitters',TR_SEX))])
	setkey(z, REC_SEX, REC_AGE_AT_CONCPOS_C, TR_AGE_AT_CONCPOS_C)
	ggplot(z, aes(x=REC_AGE_AT_CONCPOS_C, y=M, min=CL, lower=IL, upper=IU, max=CU, fill=TR_AGE_AT_CONCPOS_C)) + 
			geom_bar(stat='identity', position='dodge') +
			geom_errorbar(width=0.2, position=position_dodge(width=0.9)) +
			scale_fill_viridis(option="magma", discrete=TRUE, begin=0.85, end=0.15) +			
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(REC_SEX~.) +
			scale_y_continuous(labels=scales::percent, limits=c(0,0.13), expand=c(0,0), breaks=seq(0,1,0.05)) +			
			labs(x='age group of recipients',y='transmissions\nfrom age group\n',fill='age group of transmitters') 
	ggsave(file=paste0(outfile.base,'_age_mf_sources_allcombinations_barplot_adjusted.pdf'), w=10, h=7)		
	ggplot(z, aes(x=TR_AGE_AT_CONCPOS_C, y=M, min=CL, lower=IL, upper=IU, max=CU, fill=REC_AGE_AT_CONCPOS_C)) + 
			geom_bar(stat='identity', position='dodge') +
			geom_errorbar(width=0.2, position=position_dodge(width=0.9)) +
			scale_fill_viridis(option="magma", discrete=TRUE, begin=0.85, end=0.15) +			
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(TR_SEX~.) +
			scale_y_continuous(labels=scales::percent, limits=c(0,0.13), expand=c(0,0), breaks=seq(0,1,0.05)) +			
			labs(x='age group of transmitters',y='transmissions\nto age group\n',fill='age group of recipients') 
	ggsave(file=paste0(outfile.base,'_age_mf_destinations_allcombinations_barplot_adjusted.pdf'), w=10, h=7)	
	p1	<- ggplot(subset(z, TR_SEX=='male likely transmitters'), aes(x=REC_AGE_AT_CONCPOS_C, y=TR_AGE_AT_CONCPOS_C, fill=100*M)) + 
			geom_tile() +			
			scale_fill_viridis(option="magma", discrete=FALSE, begin=0.15, end=0.85, limits=c(0,9)) +
			geom_abline(slope=1, intercept=0, colour='white') +
			theme_bw() + theme(legend.position='bottom') +
			scale_x_discrete(expand=c(0,0)) +
			scale_y_discrete(expand=c(0,0)) +
			labs(x='female likely recipients',y='male likely transmitters',fill='proportion of transmissions') 	
	p2	<- ggplot(subset(z, TR_SEX=='female likely transmitters'), aes(x=REC_AGE_AT_CONCPOS_C, y=TR_AGE_AT_CONCPOS_C, fill=100*M)) + 
			geom_tile() +			
			scale_fill_viridis(option="magma", discrete=FALSE, begin=0.15, end=0.85, limits=c(0,9)) +
			geom_abline(slope=1, intercept=0, colour='white') +
			theme_bw() + theme(legend.position='bottom') +
			scale_x_discrete(expand=c(0,0)) +
			scale_y_discrete(expand=c(0,0)) +
			labs(x='male likely recipients',y='female likely transmitters',fill='proportion of transmissions') 	
	p3		<- g_legend(p1)	
	p3$vp	<- viewport(layout.pos.row=2, layout.pos.col=1:2)	
	pdf(file=paste0(outfile.base,'_age_mf_sources_allcombinations_matrix_adjusted.pdf'), w=8, h=5)
	grid.newpage()	
	pushViewport(viewport(layout = grid.layout(2, 2, heights=unit(c(4,1), "null"), widths=unit(c(5, 5), "null"))))   	
	print(p1+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
	print(p2+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
	grid.draw(p3)
	dev.off()
	
	#
	#	male and female transmitters
	#
	setkey(dcb, TR_SEX, TR_AGE_AT_CONCPOS_C)
	z	<- dcb[, list(PI_I_ALPHA= sum(PI_IJ_ALPHA)), by=c('TR_SEX','TR_AGE_AT_CONCPOS_C','MONTE_CARLO_IT')]
	setkey(z, TR_SEX, MONTE_CARLO_IT, TR_AGE_AT_CONCPOS_C)
	#	aggregate and get quantiles
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_I_ALPHA)
				colnames(tmp)	<- TR_AGE_AT_CONCPOS_C
				tmp		<- as.data.table(tmp)								
			}, by=c('TR_SEX','MONTE_CARLO_IT')]
	z		<- melt(z, id.vars=c('TR_SEX','MONTE_CARLO_IT'), variable.name='TR_AGE_AT_CONCPOS_C', value.name='PI_I')
	set(z, NULL, 'TR_AGE_AT_CONCPOS_C', z[, as.character(TR_AGE_AT_CONCPOS_C)])
	setkey(z, TR_SEX, TR_AGE_AT_CONCPOS_C)
	z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(PI_I, p=seq(0,1,0.01)))), by=c('TR_SEX','TR_AGE_AT_CONCPOS_C')]
	#	subset to main quantities of interest
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), TR_SEX+TR_AGE_AT_CONCPOS_C~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('TR_CL','TR_IL','TR_M','TR_IU','TR_CU'))
	set(z, NULL, 'TR_SEX', z[, gsub('M','male likely transmitters', gsub('F','female likely transmitters', TR_SEX))])
	ggplot(z, aes(x=TR_AGE_AT_CONCPOS_C, middle=TR_M, min=TR_CL, lower=TR_IL, upper=TR_IU, max=TR_CU, fill=TR_SEX)) + 
			geom_boxplot(stat='identity', position=position_dodge(width=0.9)) +
			scale_fill_manual(values=c('male likely transmitters'='skyblue2', 'female likely transmitters'='hotpink2')) +
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(~TR_SEX) +
			scale_y_continuous(labels=scales::percent, limits=c(0,0.4), expand=c(0,0), breaks=seq(0,1,0.1)) +			
			labs(x='\nage group of likely transmitters',y='transmissions\nfrom age group\n',fill='') 
	ggsave(file=paste0(outfile.base,'_age_mf_sources_barplot_adjusted.pdf'), w=8, h=5)	
	
	
	#
	#	male and female recipients
	#
	setkey(dcb, REC_SEX, REC_AGE_AT_CONCPOS_C)
	z	<- dcb[, list(PI_J_ALPHA= sum(PI_IJ_ALPHA)), by=c('REC_SEX','REC_AGE_AT_CONCPOS_C','MONTE_CARLO_IT')]
	setkey(z, REC_SEX, MONTE_CARLO_IT, REC_AGE_AT_CONCPOS_C)
	#	aggregate and get quantiles
	mc.it	<- 1e2
	z		<- z[, {												
				tmp		<- rdirichlet(mc.it, PI_J_ALPHA)
				colnames(tmp)	<- REC_AGE_AT_CONCPOS_C
				tmp		<- as.data.table(tmp)								
			}, by=c('REC_SEX','MONTE_CARLO_IT')]
	z		<- melt(z, id.vars=c('REC_SEX','MONTE_CARLO_IT'), variable.name='REC_AGE_AT_CONCPOS_C', value.name='PI_J')
	set(z, NULL, 'REC_AGE_AT_CONCPOS_C', z[, as.character(REC_AGE_AT_CONCPOS_C)])
	setkey(z, REC_SEX, REC_AGE_AT_CONCPOS_C)
	z		<- z[, list(P=seq(0,1,0.01), Q=unname(quantile(PI_J, p=seq(0,1,0.01)))), by=c('REC_SEX','REC_AGE_AT_CONCPOS_C')]
	#	subset to main quantities of interest
	z		<- dcast.data.table(subset(z, P%in%c(0.03, 0.25, 0.5, 0.75, 0.97)), REC_SEX+REC_AGE_AT_CONCPOS_C~P, value.var='Q')
	setnames(z, c('0.03','0.25','0.5','0.75','0.97'), c('REC_CL','REC_IL','REC_M','REC_IU','REC_CU'))
	set(z, NULL, 'REC_SEX', z[, gsub('M','male likely recipients', gsub('F','female likely recipients', REC_SEX))])	
	ggplot(z, aes(x=REC_AGE_AT_CONCPOS_C, middle=REC_M, min=REC_CL, lower=REC_IL, upper=REC_IU, max=REC_CU, fill=REC_SEX)) + 
			geom_boxplot(stat='identity', position=position_dodge(width=0.9)) +
			scale_fill_manual(values=c('male likely recipients'='skyblue2', 'female likely recipients'='hotpink2')) +
			theme_bw() + theme(legend.position='bottom') +
			facet_grid(~REC_SEX) +
			scale_y_continuous(labels=scales::percent, limits=c(0,0.4), expand=c(0,0), breaks=seq(0,1,0.1)) +			
			labs(x='\nage group of likely recipients',y='transmissions\nto age group\n',fill='') 
	ggsave(file=paste0(outfile.base,'_age_mf_destinations_barplot_adjusted.pdf'), w=8, h=5)	
	
	
	#
	#	Tulio mixing matrix, stratifying by female, male transmitters
	#	crude
	z	<- subset(dc, select=c(REC_SEX, REC_AGE_AT_CONCPOS_C, TR_SEX, TR_AGE_AT_CONCPOS_C, TR_OBS))	
	z	<- merge(z, z[, list(TR_SUM=sum(TR_OBS)), by='TR_SEX'],by='TR_SEX')
	z[, TR_P:=TR_OBS/TR_SUM]
	set(z, NULL, 'REC_SEX', z[,gsub('M','male likely recipients',gsub('F','female likely recipients',REC_SEX))])
	set(z, NULL, 'TR_SEX', z[,gsub('M','male likely transmitters',gsub('F','female likely transmitters',TR_SEX))])
	setkey(z, REC_SEX, REC_AGE_AT_CONCPOS_C, TR_AGE_AT_CONCPOS_C)
	p1	<- ggplot(subset(z, TR_SEX=='male likely transmitters'), aes(x=REC_AGE_AT_CONCPOS_C, y=TR_AGE_AT_CONCPOS_C, fill=100*TR_P)) + 
			geom_tile() +			
			scale_fill_viridis(option="magma", discrete=FALSE, begin=0.15, end=0.85, limits=c(0,14)) +
			geom_abline(slope=1, intercept=0, colour='white') +
			theme_bw() + theme(legend.position='bottom') +
			scale_x_discrete(expand=c(0,0)) +
			scale_y_discrete(expand=c(0,0)) +
			labs(x='female likely recipients',y='male likely transmitters',fill='proportion of transmissions') 	
	p2	<- ggplot(subset(z, TR_SEX=='female likely transmitters'), aes(x=REC_AGE_AT_CONCPOS_C, y=TR_AGE_AT_CONCPOS_C, fill=100*TR_P)) + 
			geom_tile() +			
			scale_fill_viridis(option="magma", discrete=FALSE, begin=0.15, end=0.85, limits=c(0,14)) +
			geom_abline(slope=1, intercept=0, colour='white') +
			theme_bw() + theme(legend.position='bottom') +
			scale_x_discrete(expand=c(0,0)) +
			scale_y_discrete(expand=c(0,0)) +
			labs(x='male likely recipients',y='female likely transmitters',fill='proportion of transmissions') 	
	p3		<- g_legend(p1)	
	p3$vp	<- viewport(layout.pos.row=2, layout.pos.col=1:2)	
	pdf(file=paste0(outfile.base,'_age_mf_sources_allcombinations_matrix_crude.pdf'), w=8, h=5)
	grid.newpage()	
	pushViewport(viewport(layout = grid.layout(2, 2, heights=unit(c(4,1), "null"), widths=unit(c(5, 5), "null"))))   	
	print(p1+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
	print(p2+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
	grid.draw(p3)
	dev.off()
	
	#	--> looks like there is way more transmitters <25 into older groups than from older groups into <25 years recipients 
	#	--> no substantial difference in M->F vs F->M
	library(ggExtra)
	z	<- copy(rtpdm)
	z[, COUPLE_C:= gsub('couple F->M at enrollment|couple M->F at enrollment|couple seroinc at enrollment|couple seropos at enrollment','married couple',gsub('no couple','casual pair',COUPLE))]
	set(z, NULL, 'TYPE', z[, gsub('mf','direction of transmission m->f',gsub('fm','direction of transmission f->m',TYPE))])	
	p1	<- ggplot(subset(z, TYPE=='direction of transmission m->f'), aes(y=MALE_AGE_AT_CONCPOS, x=FEMALE_AGE_AT_CONCPOS)) +
			geom_bin2d(aes(alpha = ..count..), breaks=list(x=seq(0,55,5), y=seq(0,55,5))) +
			geom_density_2d(bins=3) +
			geom_point(colour='black') +			
			geom_abline(intercept=0, slope=1, colour='black', linetype='dotted') +
			scale_x_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +			
			theme_bw() + theme(legend.position='bottom', legend.box="vertical") +
			labs(x='\nage of female likely recipient\n(at time concordant positive)', y='age of male likely transmitter\n(at time concordant positive)\n', colour='couple status') +
			guides(alpha='none', pch='none')
	p1	<- ggExtra::ggMarginal(p1, type="histogram", binwidth=2.5)
	ggsave(p1, file=paste0(outfile.base,'_age_dots_bygender_M.pdf'), w=5, h=5)	
	p1	<- ggplot(subset(z, TYPE=='direction of transmission f->m'), aes(x=MALE_AGE_AT_CONCPOS, y=FEMALE_AGE_AT_CONCPOS)) +
			geom_bin2d(aes(alpha = ..count..), breaks=list(x=seq(0,55,5), y=seq(0,55,5))) +
			geom_density_2d(bins=3) +
			geom_point(colour='black') +			
			geom_abline(intercept=0, slope=1, colour='black', linetype='dotted') +
			scale_x_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +			
			theme_bw() + theme(legend.position='bottom', legend.box="vertical") +
			labs(x='\nage of male likely recipient\n(at time concordant positive)', y='age of female likely transmitter\n(at time concordant positive)\n', colour='couple status') +
			guides(alpha='none', pch='none')
	p1	<- ggExtra::ggMarginal(p1, type="histogram", binwidth=2.5)
	ggsave(p1, file=paste0(outfile.base,'_age_dots_bygender_F.pdf'), w=5, h=5)
	
	
	#
	#	Tulio mixing matrix, stratifying by casual pair / couple
	#
	z	<- copy(rtpdm)
	z[, COUPLE_C:= gsub('couple F->M at enrollment|couple M->F at enrollment|couple seroinc at enrollment|couple seropos at enrollment','married couple',gsub('no couple','casual pair',COUPLE))]
	set(z, NULL, 'TYPE', z[, gsub('mf','direction of transmission m->f',gsub('fm','direction of transmission f->m',TYPE))])	
	p1	<- ggplot(subset(z, TYPE=='direction of transmission m->f'), aes(y=MALE_AGE_AT_CONCPOS, x=FEMALE_AGE_AT_CONCPOS)) +
			geom_bin2d(aes(alpha = ..count..), breaks=list(x=seq(0,55,5), y=seq(0,55,5))) +
			geom_density_2d(aes(colour=COUPLE_C), bins=3) +
			geom_point(aes(pch=COUPLE_C, colour=COUPLE_C)) +			
			geom_abline(intercept=0, slope=1, colour='black', linetype='dotted') +
			scale_x_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
			scale_colour_manual(values=c('married couple'='black', 'casual pair'='DarkRed')) +
			theme_bw() + theme(legend.position='bottom', legend.box="vertical") +
			labs(x='\nage of female likely recipient\n(at time concordant positive)', y='age of male likely transmitter\n(at time concordant positive)\n', colour='couple status') +
			guides(alpha='none', pch='none')
	p1	<- ggExtra::ggMarginal(p1, type="histogram", binwidth=2.5)
	ggsave(p1, file=paste0(outfile.base,'_age_dots_bygender_M_bycouplestatus.pdf'), w=5, h=5)	
	p1	<- ggplot(subset(z, TYPE=='direction of transmission f->m'), aes(x=MALE_AGE_AT_CONCPOS, y=FEMALE_AGE_AT_CONCPOS)) +
			geom_bin2d(aes(alpha = ..count..), breaks=list(x=seq(0,55,5), y=seq(0,55,5))) +
			geom_density_2d(aes(colour=COUPLE_C), bins=3) +
			geom_point(aes(pch=COUPLE_C, colour=COUPLE_C)) +			
			geom_abline(intercept=0, slope=1, colour='black', linetype='dotted') +
			scale_x_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
			scale_colour_manual(values=c('married couple'='black', 'casual pair'='DarkRed')) +
			theme_bw() + theme(legend.position='bottom', legend.box="vertical") +
			labs(x='\nage of male likely recipient\n(at time concordant positive)', y='age of female likely transmitter\n(at time concordant positive)\n', colour='couple status') +
			guides(alpha='none', pch='none')
	p1	<- ggExtra::ggMarginal(p1, type="histogram", binwidth=2.5)
	ggsave(p1, file=paste0(outfile.base,'_age_dots_bygender_F_bycouplestatus.pdf'), w=5, h=5)
	


	#
	#	Tulio mixing matrix, stratifying by location female
	#
	z	<- copy(rtpdm)
	z[, COUPLE_C:= gsub('couple F->M at enrollment|couple M->F at enrollment|couple seroinc at enrollment|couple seropos at enrollment','married couple',gsub('no couple','casual pair',COUPLE))]
	set(z, NULL, 'TYPE', z[, gsub('mf','direction of transmission m->f',gsub('fm','direction of transmission f->m',TYPE))])	
	p1	<- ggplot(subset(z, TYPE=='direction of transmission m->f'), aes(y=MALE_AGE_AT_CONCPOS, x=FEMALE_AGE_AT_CONCPOS)) +
			geom_bin2d(aes(alpha = ..count..), breaks=list(x=seq(0,55,5), y=seq(0,55,5))) +
			geom_density_2d(aes(colour=MALE_COMM_TYPE), bins=3) +
			geom_point(aes(pch=MALE_COMM_TYPE, colour=MALE_COMM_TYPE)) +			
			geom_abline(intercept=0, slope=1, colour='black', linetype='dotted') +
			scale_x_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
			scale_colour_manual(values=c('fisherfolk'='black', 'trading'='DarkRed', 'agrarian'='DarkGreen')) +
			theme_bw() + theme(legend.position='bottom', legend.box="vertical") +
			labs(x='\nage of female likely recipient\n(at time concordant positive)', y='age of male likely transmitter\n(at time concordant positive)\n', colour='location likely transmitter') +
			guides(alpha='none', pch='none')
	p1	<- ggExtra::ggMarginal(p1, type="histogram", binwidth=2.5)
	ggsave(p1, file=paste0(outfile.base,'_age_dots_bygender_M_byloctransmitter.pdf'), w=5, h=5)	
	p1	<- ggplot(subset(z, TYPE=='direction of transmission f->m'), aes(x=MALE_AGE_AT_CONCPOS, y=FEMALE_AGE_AT_CONCPOS)) +
			geom_bin2d(aes(alpha = ..count..), breaks=list(x=seq(0,55,5), y=seq(0,55,5))) +
			geom_density_2d(aes(colour=FEMALE_COMM_TYPE), bins=3) +
			geom_point(aes(pch=FEMALE_COMM_TYPE, colour=FEMALE_COMM_TYPE)) +			
			geom_abline(intercept=0, slope=1, colour='black', linetype='dotted') +
			scale_x_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
			scale_colour_manual(values=c('fisherfolk'='black', 'trading'='DarkRed', 'agrarian'='DarkGreen')) +
			theme_bw() + theme(legend.position='bottom', legend.box="vertical") +
			labs(x='\nage of male likely recipient\n(at time concordant positive)', y='age of female likely transmitter\n(at time concordant positive)\n', colour='location likely transmitter') +
			guides(alpha='none', pch='none')
	p1	<- ggExtra::ggMarginal(p1, type="histogram", binwidth=2.5)
	ggsave(p1, file=paste0(outfile.base,'_age_dots_bygender_F_byloctransmitter.pdf'), w=5, h=5)


	#
	#	age difference
	#
	z		<- subset(rtpdm, !is.na(AGEDIFF))	
	z[, MALE_AGE_AT_CONCPOS:= pmax(MALE_FIRSTPOSDATE,FEMALE_FIRSTPOSDATE)-MALE_BIRTHDATE]
	z[, FEMALE_AGE_AT_CONCPOS:= pmax(MALE_FIRSTPOSDATE,FEMALE_FIRSTPOSDATE)-FEMALE_BIRTHDATE]		
	z[, COUPLE_C:= gsub('couple F->M at enrollment|couple M->F at enrollment|couple seroinc at enrollment|couple seropos at enrollment','married couple',gsub('no couple','casual pair',COUPLE))]	
	ggplot(z, aes(x=as.character(factor(TYPE, levels=c('fm','mf'), labels=c('female is\nlikely transmitter','female is\nlikely recipient'))), y=AGEDIFF)) + 
			geom_violin(bw=4, fill='LightBlue', trim=TRUE, draw_quantiles=0.5) + 
			geom_dotplot(binaxis='y', binwidth=1, stackdir='center', fill='DarkBlue', dotsize=1, width =0.8, stackratio = 1) +
			facet_grid(COUPLE_C~FEMALE_COMM_TYPE) +
			scale_y_continuous() +				
			theme_bw() + theme(legend.position='bottom', plot.title = element_text(hjust = 0.5)) +
			guides(fill=guide_legend(ncol=2)) +
			labs(x='', y=paste0('age difference\n(years female younger)','\n'))
	ggsave(file=paste0(outfile.base,'_agediff_by_commtype_couplestatus.pdf'), w=10, h=10)
	#	--> female recipients are much younger in casual pairs than male partners
	#		whereas female transmitters are of similar age compare to male partners
	#
	#	age difference by age of female at time conc positive
	#
	p1	<- ggplot(subset(z, TYPE=='fm'), aes(x=FEMALE_AGE_AT_CONCPOS, y=AGEDIFF)) +
		geom_bin2d(aes(alpha = ..count..), breaks=list(x=seq(0,55,5), y=seq(-20,20,5))) +
		geom_smooth(se=FALSE, span=1, method='loess') +	
		geom_point() +
		geom_abline(intercept=0, slope=0, colour='black', linetype='dotted') +
		scale_x_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
		scale_y_continuous(breaks=seq(-20,20,5), limits=c(-20, 20), expand=c(0,0)) +
		#scale_colour_manual(values=c('fisherfolk'='black', 'trading'='DarkRed', 'agrarian'='DarkGreen')) +
		theme_bw() + theme(legend.position='bottom', legend.box="vertical") +
		labs(x='\nage of female likely transmitter\n(at time concordant positive)', y='age difference\n(years female younger)\n') +
		guides(alpha='none', pch='none')
	p1	<- ggExtra::ggMarginal(p1, type="histogram", binwidth=2.5)
	ggsave(p1, file=paste0(outfile.base,'_agediff_bygender_F.pdf'), w=5, h=5)
	p1	<- ggplot(subset(z, TYPE=='mf'), aes(x=FEMALE_AGE_AT_CONCPOS, y=AGEDIFF)) +
			geom_bin2d(aes(alpha = ..count..), breaks=list(x=seq(0,55,5), y=seq(-20,20,5))) +
			geom_smooth(se=FALSE, span=1, method='loess') +	
			geom_point() +
			geom_abline(intercept=0, slope=0, colour='black', linetype='dotted') +
			scale_x_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(-20,20,5), limits=c(-20, 20), expand=c(0,0)) +
			#scale_colour_manual(values=c('fisherfolk'='black', 'trading'='DarkRed', 'agrarian'='DarkGreen')) +
			theme_bw() + theme(legend.position='bottom', legend.box="vertical") +
			labs(x='\nage of female likely recipient\n(at time concordant positive)', y='age difference\n(years female younger)\n') +
			guides(alpha='none', pch='none')
	p1	<- ggExtra::ggMarginal(p1, type="histogram", binwidth=2.5)
	ggsave(p1, file=paste0(outfile.base,'_agediff_bygender_M.pdf'), w=5, h=5)
	#
	#	age difference by age of female at time conc positive and couple status
	#
	p1	<- ggplot(subset(z, TYPE=='fm'), aes(x=FEMALE_AGE_AT_CONCPOS, y=AGEDIFF)) +
			geom_bin2d(aes(alpha = ..count..), breaks=list(x=seq(0,55,5), y=seq(-20,20,5))) +
			geom_density_2d(aes(colour=COUPLE_C), bins=3) +
			geom_smooth(se=FALSE, span=1, method='loess') +	
			geom_point(aes(pch=COUPLE_C, colour=COUPLE_C)) +
			geom_abline(intercept=0, slope=0, colour='black', linetype='dotted') +
			scale_x_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(-20,20,5), limits=c(-20, 20), expand=c(0,0)) +
			scale_colour_manual(values=c('married couple'='black', 'casual pair'='DarkRed')) +
			theme_bw() + theme(legend.position='bottom', legend.box="vertical") +
			labs(x='\nage of female likely transmitter\n(at time concordant positive)', y='age difference\n(years female younger)\n', colour='couple status') +
			guides(alpha='none', pch='none')
	p1	<- ggExtra::ggMarginal(p1, type="histogram", binwidth=2.5)
	ggsave(p1, file=paste0(outfile.base,'_agediff_bygender_F_bycouplestatus.pdf'), w=5, h=5)
	p1	<- ggplot(subset(z, TYPE=='mf'), aes(x=FEMALE_AGE_AT_CONCPOS, y=AGEDIFF)) +
			geom_bin2d(aes(alpha = ..count..), breaks=list(x=seq(0,55,5), y=seq(-20,20,5))) +
			geom_density_2d(aes(colour=COUPLE_C), bins=3) +
			geom_smooth(se=FALSE, span=1, method='loess') +	
			geom_point(aes(pch=COUPLE_C, colour=COUPLE_C)) +
			geom_abline(intercept=0, slope=0, colour='black', linetype='dotted') +
			scale_x_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(-20,20,5), limits=c(-20, 20), expand=c(0,0)) +
			scale_colour_manual(values=c('married couple'='black', 'casual pair'='DarkRed')) +
			theme_bw() + theme(legend.position='bottom', legend.box="vertical") +
			labs(x='\nage of female likely recipient\n(at time concordant positive)', y='age difference\n(years female younger)\n', colour='couple status') +
			guides(alpha='none', pch='none')
	p1	<- ggExtra::ggMarginal(p1, type="histogram", binwidth=2.5)
	ggsave(p1, file=paste0(outfile.base,'_agediff_bygender_M_bycouplestatus.pdf'), w=5, h=5)
	#
	#	age difference by age of female at time conc positive and location status of transmitter
	#
	p1	<- ggplot(subset(z, TYPE=='fm'), aes(x=FEMALE_AGE_AT_CONCPOS, y=AGEDIFF)) +
			geom_bin2d(aes(alpha = ..count..), breaks=list(x=seq(0,55,5), y=seq(-20,20,5))) +
			geom_density_2d(aes(colour=FEMALE_COMM_TYPE), bins=3) +
			geom_smooth(se=FALSE, span=1, method='loess') +	
			geom_point(aes(pch=FEMALE_COMM_TYPE, colour=FEMALE_COMM_TYPE)) +
			geom_abline(intercept=0, slope=0, colour='black', linetype='dotted') +
			scale_x_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(-20,20,5), limits=c(-20, 20), expand=c(0,0)) +
			scale_colour_manual(values=c('fisherfolk'='black', 'trading'='DarkRed', 'agrarian'='DarkGreen')) +
			theme_bw() + theme(legend.position='bottom', legend.box="vertical") +
			labs(x='\nage of female likely transmitter\n(at time concordant positive)', y='age difference\n(years female younger)\n', colour='location likely transmitter') +
			guides(alpha='none', pch='none')
	p1	<- ggExtra::ggMarginal(p1, type="histogram", binwidth=2.5)
	ggsave(p1, file=paste0(outfile.base,'_agediff_bygender_F_byloctransmitter.pdf'), w=5, h=5)
	p1	<- ggplot(subset(z, TYPE=='mf'), aes(x=FEMALE_AGE_AT_CONCPOS, y=AGEDIFF)) +
			geom_bin2d(aes(alpha = ..count..), breaks=list(x=seq(0,55,5), y=seq(-20,20,5))) +
			geom_density_2d(aes(colour=MALE_COMM_TYPE), bins=3) +
			geom_smooth(se=FALSE, span=1, method='loess') +	
			geom_point(aes(pch=MALE_COMM_TYPE, colour=MALE_COMM_TYPE)) +
			geom_abline(intercept=0, slope=0, colour='black', linetype='dotted') +
			scale_x_continuous(breaks=seq(0,55,5), limits=c(15, 55), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(-20,20,5), limits=c(-20, 20), expand=c(0,0)) +
			scale_colour_manual(values=c('fisherfolk'='black', 'trading'='DarkRed', 'agrarian'='DarkGreen')) +
			theme_bw() + theme(legend.position='bottom', legend.box="vertical") +
			labs(x='\nage of female likely recipient\n(at time concordant positive)', y='age difference\n(years female younger)\n', colour='location likely transmitter') +
			guides(alpha='none', pch='none')
	p1	<- ggExtra::ggMarginal(p1, type="histogram", binwidth=2.5)
	ggsave(p1, file=paste0(outfile.base,'_agediff_bygender_M_byloctransmitter.pdf'), w=5, h=5)
	
	
	
	
	#
	#rtr2[, table(PAIR_TYPE)]
	#
	#	did any transmitter start ART before the recipient was diagnosed?
	subset(rtr2, TR_ARVSTARTDATE<REC_FIRSTPOSDATE)	
	#	F026858:J104288 --> stable couple, rec male, first diagnosed with v high CD4 (2400), about 2 years after female started ART 
	#	C066263:K077878 --> no couple, rec female, first diagnosed with v high CD4, about 5m after male started ART
	
	#
	#	make basic epi plot: when positive, when negative, when sequenced
	#	
	t.posneg	<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
	setnames(t.posneg, c('BIRTHDATE','EST_DATEDIED'), c('DOB','DOD'))
	t.seq		<- unique(subset(rs, !is.na(PID), select=c(RID, PID, DATE)))
	setnames(t.seq, 'DATE', 'SEQ_DATE')
	t.cd4		<- unique(subset(rd, !is.na(RECENTCD4DATE) & !is.na(RECENTCD4), select=c(RID, RECENTCD4DATE, RECENTCD4)))
	set(t.cd4, NULL, 'RECENTCD4', t.cd4[, cut(RECENTCD4, breaks=c(-1,250,350,500,800,Inf), labels=c('<200','200-349','350-499','500-799','800+'))])
	t.vl		<- unique(subset(rd, !is.na(RECENTVLDATE) & !is.na(RECENTVL), select=c(RID, RECENTVLDATE, RECENTVL)))
	set(t.vl, NULL, 'RECENTVL', t.vl[, cut(RECENTVL, breaks=c(-1,200,1e3,1e5,Inf), labels=c('<200','200-1,000','1,000-100,000','100,000+'))])
	#	plot 3 timelines per page
	setkey(rtpd, TYPE, POSTERIOR_SCORE)
	rtpd[, DUMMY:= ceiling(seq_len(nrow(rtpd))/50)]	
	p			<- lapply(rtpd[, unique(DUMMY)], function(i){
				cat('\n', i)
				df		<- subset(rtpd, DUMMY==i)
				p		<- RakaiFull.plot.epitimeline(df, copy(t.posneg), copy(t.cd4), copy(t.vl), copy(t.seq), age.adult=14)				
			})
	pdf.n	<- 2
	pi	<- data.table(IDX=seq_along(p))
	pi[, PLOT:= ceiling(IDX/pdf.n)]
	pi[, PLOT_IDX:= (IDX-1)%%pdf.n+1]	
	pdf(file=paste0(outfile.base,'epilines.pdf'), w=20, h=12)	
	for(plot in pi[, unique(PLOT)])
	{
		idx			<- subset(pi, PLOT==plot)[, IDX]
		plot.idx	<- subset(pi, PLOT==plot)[, PLOT_IDX]
		grid.newpage()
		pushViewport(viewport(layout=grid.layout(1, pdf.n)))
		for(i in seq_along(idx))
			print(p[[idx[i]]], vp = viewport(layout.pos.row=1, layout.pos.col=plot.idx[i]))
	}
	dev.off()
	rtpd[, DUMMY:=NULL]
	#pdf(file=paste0(outfile.base,'epilines.pdf'), w=8, h=12)
	#print(p)
	#dev.off()
			
	

	#	plot 3 timelines per page
	setkey(rtpd, TYPE, POSTERIOR_SCORE)
	rtpd[, DUMMY:= ceiling(seq_len(nrow(rtpd))/50)]	
	p			<- lapply(rtpd[, unique(DUMMY)], function(i){
				cat('\n', i)
				df		<- subset(rtpd, DUMMY==i)
				p		<- RakaiFull.plot.epitimeline(df, copy(t.posneg), copy(t.cd4), copy(t.vl), copy(t.seq), age.adult=14)				
			})

	
	
	
	
	#
	#	how many transmitters were positive for 6m before the recipient was found positive
	subset(rtr2, (TR_FIRSTPOSDATE+.5)<REC_FIRSTPOSDATE)	
	#	26
	
	
	subset(rtr, MALE_COMM_TYPE==FEMALE_COMM_TYPE & FEMALE_COMM_TYPE!='trading')[, {
				z	<- binconf( length(which(PHSC_DIR=='m->f')), length(PHSC_DIR) )				
				list(K=length(which(PHSC_DIR=='m->f')), N=length(PHSC_DIR), P=z[1], QL=z[2], QU=z[3])
			}, by=c('MALE_COMM_TYPE')]	
	#	   MALE_COMM_TYPE  	K  N  P         QL        QU
	#1:       agrarian 		27 38 0.7105263 0.5524286 0.8299672
	#2:     fisherfolk 		55 85 0.6470588 0.5411250 0.7402751
	
	#
	#	is there a difference in male->female transmission by couple type?
	#	results: no	
	#
	tmp		<- copy(rtr)
	set(tmp, tmp[, which(PAIR_TYPE!='stable cohabiting')], 'PAIR_TYPE', 'no stable pair')
	tmp[, {
				z	<- binconf( length(which(PHSC_DIR=='m->f')), length(PHSC_DIR) )				
				list(K=length(which(PHSC_DIR=='m->f')), N=length(PHSC_DIR), P=z[1], QL=z[2], QU=z[3])
			}, by='PAIR_TYPE']	
	#			PAIR_TYPE  	 K  N         P        QL        QU
	#1:    stable cohabiting 59 86 0.6860465 0.5817960 0.7743870
	#2:    no stable pair 	 25 42 0.5952381 0.4449431 0.7295714	
	chisq_test(factor(PHSC_DIR) ~ factor(PAIR_TYPE), data=tmp, distribution="exact")
	#	chi-squared = 1.0315, p-value = 0.3278
	
	
	#	are transmitters younger in fisherfolk sites?
	#	results: yes
	tmp		<- subset(rtr2, TR_COMM_TYPE!='trading')
	ggplot(tmp, aes(x=TR_COMM_TYPE, y=TR_BIRTHDATE)) + geom_boxplot()
	independence_test(TR_BIRTHDATE~factor(TR_COMM_TYPE), data=tmp, distribution = "exact")
	#	Z = -2.3289, p-value = 0.01934
	ggplot(tmp, aes(x=REC_COMM_TYPE, y=REC_BIRTHDATE)) + geom_boxplot()
	independence_test(REC_BIRTHDATE~factor(REC_COMM_TYPE), data=tmp, distribution = "exact")
	#	Z = -2.2714, p-value = 0.02236
	#	summary(rq(TR_BIRTHDATE~TR_COMM_TYPE, tau=.5, data=tmp, method='fn'), se='nid')
	
	#
	#	is there a difference in age gap between male->female transmission / female->male transmission ?
	#	results: not significant but outside couples, men tend to be infected by much younger women
	#	
	tmp		<- subset(rtr, !grepl('cohabiting',PAIR_TYPE) & !is.na(AGEDIFF) & FEMALE_COMM_TYPE==MALE_COMM_TYPE & MALE_COMM_TYPE=='fisherfolk')
	independence_test(AGEDIFF~factor(PHSC_DIR), data=tmp, distribution = "exact")
	#	Z = 1.5902, p-value = 0.1134
	tmp		<- subset(rtr, !grepl('cohabiting',PAIR_TYPE) & !is.na(AGEDIFF) & FEMALE_COMM_TYPE==MALE_COMM_TYPE & MALE_COMM_TYPE=='agrarian')
	independence_test(AGEDIFF~factor(PHSC_DIR), data=tmp, distribution = "exact")
	#	Z = 1.7439, p-value = 0.1429
	
	
	ggplot(rtr, aes(x=PHSC_DIR, y=AGEDIFF)) + geom_boxplot()	
	tmp		<- subset(rtr, !is.na(MALE_BIRTHDATE) & !is.na(FEMALE_BIRTHDATE), select=c(PHSC_DIR,AGEDIFF))
	set(tmp, NULL, 'PHSC_DIR', tmp[, as.integer(as.character(factor(PHSC_DIR, levels=c('f->m','m->f'), labels=c('0','1'))))])
	summary(gamlss(data=tmp, PHSC_DIR~AGEDIFF, family=LO))
	summary(gamlss(data=tmp, AGEDIFF~PHSC_DIR))
	#				Estimate Std. Error t value Pr(>|t|)    
	#(Intercept)  0.626998   0.071579   8.759 1.98e-13 ***
	#AGEDIFF     -0.002135   0.008422  -0.254      0.8    
	tmp		<- subset(rtr, MALE_COMM_TYPE!='trading' & MALE_COMM_TYPE==FEMALE_COMM_TYPE)
	ggplot(tmp, aes(x=PHSC_DIR, y=AGEDIFF)) + 
			geom_boxplot() + 
			facet_grid(~MALE_COMM_TYPE)	+
			theme_bw() + labs(x='\nestimated direction of transmission', y='age difference male-female\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-directionpairs_direction-agegap-commtype.pdf',sep='')), w=4, h=6)
	
	#	AAA
	subset(rtr2, TR_COMM_TYPE!='trading' & PAIR_TYPE!='m and f not in couple' & REC_RID%in%c(rc$MALE_RID,rc$FEMALE_RID))[, {
				m<- length(which(PAIR_TYPE=='stable cohabiting'))
				n<- length(PAIR_TYPE)
				z<- unname(as.numeric(binconf(m, n)))
				list(P=z[1], PL=z[2], PU=z[3], M=m, N=n, TYPE='stable cohabiting')				
			}, by=c('TR_COMM_TYPE','REC_SEX')]
	
	#
	#	Does the primary occupation differ between transmitters / recipients? 
	#	
	tmp			<- subset(rtr, FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_occupation2.pdf',sep=''))
	tmp			<- subset(rtr, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_occupation2_stablecouples.pdf',sep=''))
	tmp			<- subset(rtr, !grepl('cohabiting',PAIR_TYPE) & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_occupation2_nocouples.pdf',sep=''))
	#	
	tmp2		<- unique(subset(ra, VISIT!=17 & SEX=='F' & !is.na(FIRSTPOSDATE) & COMM_TYPE!='trading', c(RID, OCCUP_OLLI, COMM_TYPE)))
	setnames(tmp2, colnames(tmp2), paste0('FEMALE_',colnames(tmp2)))
	tmp2[, PHSC_DIR:='denominator']
	tmp			<- rbind(tmp, tmp2, fill=TRUE)
	tmp2		<- unique(subset(ra, VISIT!=17 & SEX=='M' & !is.na(FIRSTPOSDATE) & COMM_TYPE!='trading', c(RID, OCCUP_OLLI, COMM_TYPE)))
	setnames(tmp2, colnames(tmp2), paste0('MALE_',colnames(tmp2)))
	tmp2[, PHSC_DIR:='denominator']
	tmp			<- rbind(tmp, tmp2, fill=TRUE)
	setnames(tmp, c('FEMALE_OCCUP_OLLI','MALE_OCCUP_OLLI'), c('FEMALE_FACTOR','MALE_FACTOR'))
	cols		<- tmp[, length(na.omit(unique(c(FEMALE_FACTOR, MALE_FACTOR))))]
	cols		<- colorRampPalette(brewer.pal(min(11,cols), "Set3"))( cols )		
	names(cols)	<- tmp[, na.omit(sort(unique(c(FEMALE_FACTOR, MALE_FACTOR))))]
	Rakai.plot.directed.pairs.discrete1(tmp, cols, outfile, 'occupation at diagnosis', w=10, h=7)
	
	# number female Bar/waitress that are transmitters
	ntf	<- nrow(subset(tmp, PHSC_DIR=='f->m' & FEMALE_FACTOR=='Bar/waitress' & FEMALE_COMM_TYPE=='fisherfolk'))
	# number female Bar/waitress that are recipients
	nrf	<- nrow(subset(tmp, PHSC_DIR=='m->f' & FEMALE_FACTOR=='Bar/waitress' & FEMALE_COMM_TYPE=='fisherfolk'))
	# number female Bar/waitress HIV-infected
	ndf	<- nrow(subset(tmp, PHSC_DIR=='denominator' & FEMALE_FACTOR=='Bar/waitress' & FEMALE_COMM_TYPE=='fisherfolk'))
	# number female that are transmitters
	nt	<- nrow(subset(tmp, PHSC_DIR=='f->m' & FEMALE_COMM_TYPE=='fisherfolk'))
	# number female that are recipients
	nr	<- nrow(subset(tmp, PHSC_DIR=='m->f' & FEMALE_COMM_TYPE=='fisherfolk'))
	# number female HIV-infected
	nd	<- nrow(subset(tmp, PHSC_DIR=='denominator' & FEMALE_COMM_TYPE=='fisherfolk'))
	# odds ratio transmitter / recipient
	# 'a' is exposed cases (exposed=Bar/waitress, case=transmitter)
	# I resample by taking p=ntf/nt as the best estimate of the proportion of female Bar/waitress that are transmitters
	# and then adding uncertainty around p based on p(1-p)/n
	bs	<- 1e4
	a	<- round(nt*rnorm(bs, mean=ntf/nt, sd=sqrt( (ntf/nt) * (1-ntf/nt) / nt )))
	# 'b' is exposed non-cases (exposed=Bar/waitress, non-case=recipients)
	b	<- round(nr*rnorm(bs, mean=nrf/nr, sd=sqrt( (nrf/nr) * (1-nrf/nr) / nr )))
	c	<- nt-a
	d	<- nr-b
	tmp2<- quantile( (a/c) / (b/d), p=c(0.025,0.975))
	a	<- ntf
	b	<- nrf
	c	<- nt-a
	d	<- nr-b
	tmp2<- c( (a/c) / (b/d), tmp2)
	
	bs	<- 1e4
	tmp2<- quantile( rnorm(bs, mean=ntf/nt, sd=sqrt( (ntf/nt) * (1-ntf/nt) / nt )) / rnorm(bs, mean=ndf/nd, sd=sqrt( (ndf/nd) * (1-ndf/nd) / nd )), p=c(0.025,0.975))
	tmp2<- c((ntf/nt) / (ndf/nd),tmp2)
	tmp3<- quantile( rnorm(bs, mean=nrf/nr, sd=sqrt( (nrf/nr) * (1-nrf/nr) / nr )) / rnorm(bs, mean=ndf/nd, sd=sqrt( (ndf/nd) * (1-ndf/nd) / nd )), p=c(0.025,0.975))
	tmp3<- c((nrf/nr) / (ndf/nd),tmp3)
	
	
	
	ntf	<- nrow(subset(tmp, PHSC_DIR=='m->f' & MALE_FACTOR=='Fishing' & MALE_COMM_TYPE=='fisherfolk'))
	nrf	<- nrow(subset(tmp, PHSC_DIR=='f->m' & MALE_FACTOR=='Fishing' & MALE_COMM_TYPE=='fisherfolk'))
	ndf	<- nrow(subset(tmp, PHSC_DIR=='denominator' & MALE_FACTOR=='Fishing' & MALE_COMM_TYPE=='fisherfolk'))
	nt	<- nrow(subset(tmp, PHSC_DIR=='m->f' & MALE_COMM_TYPE=='fisherfolk'))
	nr	<- nrow(subset(tmp, PHSC_DIR=='f->m' & MALE_COMM_TYPE=='fisherfolk'))
	nd	<- nrow(subset(tmp, PHSC_DIR=='denominator' & MALE_COMM_TYPE=='fisherfolk'))
	
	bs	<- 1e4
	tmp	<- quantile( rnorm(bs, mean=ntf/nt, sd=sqrt( (ntf/nt) * (1-ntf/nt) / nt )) / rnorm(bs, mean=ndf/nd, sd=sqrt( (ndf/nd) * (1-ndf/nd) / nd )), p=c(0.025,0.975))
	tmp	<- c((ntf/nt) / (ndf/nd),tmp)
	tmp	<- quantile( rnorm(bs, mean=nrf/nr, sd=sqrt( (nrf/nr) * (1-nrf/nr) / nr )) / rnorm(bs, mean=ndf/nd, sd=sqrt( (ndf/nd) * (1-ndf/nd) / nd )), p=c(0.025,0.975))
	tmp	<- c((nrf/nr) / (ndf/nd),tmp)
	#
	#	Age group
	#	
	tmp			<- subset(rtr, FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_birthdate.pdf',sep=''))
	setnames(tmp, c('FEMALE_BIRTHDATE','MALE_BIRTHDATE'), c('FEMALE_FACTOR','MALE_FACTOR'))	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'birth date', bw=4, dotsize=0.5, w=10, h=7)
	tmp			<- subset(rtr, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_birthdate_stablecouples.pdf',sep=''))
	setnames(tmp, c('FEMALE_BIRTHDATE','MALE_BIRTHDATE'), c('FEMALE_FACTOR','MALE_FACTOR'))	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'birth date', bw=4, dotsize=0.5, w=10, h=7)
	tmp			<- subset(rtr, PAIR_TYPE=='not registered as couple' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_birthdate_nocouples.pdf',sep=''))
	setnames(tmp, c('FEMALE_BIRTHDATE','MALE_BIRTHDATE'), c('FEMALE_FACTOR','MALE_FACTOR'))	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'birth date', bw=4, dotsize=0.5, w=10, h=7)
	#	
	tmp			<- subset(rtr, FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_agediff.pdf',sep=''))
	setnames(tmp, c('AGEDIFF'), c('MALE_FACTOR'))
	tmp[, FEMALE_FACTOR:= MALE_FACTOR]	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'age difference\n(years female younger)', bw=4, dotsize=0.5, w=10, h=7)
	tmp			<- subset(rtr, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_agediff_stablecouples.pdf',sep=''))
	setnames(tmp, c('AGEDIFF'), c('MALE_FACTOR'))
	tmp[, FEMALE_FACTOR:= MALE_FACTOR]	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'age difference\n(years female younger)', bw=4, dotsize=0.5, w=10, h=7)
	tmp			<- subset(rtr, PAIR_TYPE=='not registered as couple' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_agediff_nocouples.pdf',sep=''))
	setnames(tmp, c('AGEDIFF'), c('MALE_FACTOR'))
	tmp[, FEMALE_FACTOR:= MALE_FACTOR]	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'age difference\n(years female younger)', bw=4, dotsize=0.5, w=10, h=7)
	
	#
	#	Marriage Status
	#
	tmp			<- subset(rtr, FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_marriagestatus.pdf',sep=''))
	tmp			<- subset(rtr, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_marriagestatus_stablecouples.pdf',sep=''))
	tmp			<- subset(rtr, !grepl('cohabiting',PAIR_TYPE) & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_marriagestatus_nocouples.pdf',sep=''))
	
	set(tmp, NULL, 'MALE_MARSTAT', tmp[, gsub('Never Married \\+ casual partner','Never Married',gsub('Previously Married \\+ casual partner','Previously Married',MALE_MARSTAT))])
	set(tmp, NULL, 'FEMALE_MARSTAT', tmp[, gsub('Never Married \\+ casual partner','Never Married',gsub('Previously Married \\+ casual partner','Previously Married',FEMALE_MARSTAT))])
	tmp2		<- unique(subset(ra, SEX=='F' & !is.na(FIRSTPOSDATE) & COMM_TYPE!='trading', c(RID, MARSTAT, COMM_TYPE)))
	setnames(tmp2, colnames(tmp2), paste0('FEMALE_',colnames(tmp2)))
	tmp2[, PHSC_DIR:='denominator']
	tmp			<- rbind(tmp, tmp2, fill=TRUE)
	tmp2		<- unique(subset(ra, SEX=='M' & !is.na(FIRSTPOSDATE) & COMM_TYPE!='trading', c(RID, MARSTAT, COMM_TYPE)))
	setnames(tmp2, colnames(tmp2), paste0('MALE_',colnames(tmp2)))
	tmp2[, PHSC_DIR:='denominator']
	tmp			<- rbind(tmp, tmp2, fill=TRUE)	
	setnames(tmp, c('FEMALE_MARSTAT','MALE_MARSTAT'), c('FEMALE_FACTOR','MALE_FACTOR'))
	cols		<- tmp[, length(unique(c(FEMALE_FACTOR, MALE_FACTOR)))]
	cols		<- colorRampPalette(brewer.pal(min(8,cols), "Set2"))( cols )	
	#cols		<- rainbow_hcl(tmp[, length(unique(c(FEMALE_FACTOR, MALE_FACTOR)))], start = 20, end = 340, c=100, l=60)
	names(cols)	<- tmp[, sort(unique(c(FEMALE_FACTOR, MALE_FACTOR)))]
	Rakai.plot.directed.pairs.discrete1(tmp, cols, outfile, 'marital &\nself-reported\nnon-marital\nrelationships,\n', w=10, h=7)
	
	
	
	#
	#	Education
	#
	tmp			<- subset(rtr, FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	setnames(tmp, c('FEMALE_EDUCAT','MALE_EDUCAT'), c('FEMALE_FACTOR','MALE_FACTOR'))
	cols		<- tmp[, length(unique(c(FEMALE_FACTOR, MALE_FACTOR)))]
	cols		<- colorRampPalette(brewer.pal(min(11,cols), "Set1"))( cols )	
	#cols		<- rainbow_hcl(tmp[, length(unique(c(FEMALE_FACTOR, MALE_FACTOR)))], start = 20, end = 340, c=100, l=60)
	names(cols)	<- tmp[, sort(unique(c(FEMALE_FACTOR, MALE_FACTOR)))]
	Rakai.plot.directed.pairs.discrete1(tmp, cols, file.path(dir, paste(run,'-phsc-directionpairs_education.pdf',sep='')), 'Educational status', w=10, h=7)	
}

RakaiFull.preprocess.couples.todi.addingmetadata.170522<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	
	confidence.cut			<- 0.5	# this was used previously for the probability mass selection criterion that is now replaced with the posterior mode criterion
	confidence.cut			<- 0.66	# do not change, because the prior is calibrated for 0.66
	#infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_couples_170428_cl3.rda'
	#outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_couples_170428_"
	#infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_couples_170516_cl3.rda'
	#outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170522/todi_couples_170522_"
	infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170610/todi_couples_170610_cl3_prior23.rda"
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170610/todi_couples_170610_cl3_prior23_"
	infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170610/todi_couples_170610_cl3_prior34d23.rda"
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170610/todi_couples_170610_cl3_prior34d23_"
	
	outfile.save			<- paste0(outfile.base, 'withmetadata.rda')
	load(infile.trmpairs.todi)
	#	
	rtp.tpairs	<- rtp.todi2
	rtp.tpairs	<- subset(rtp.tpairs, ID1%in%rd$RID & ID2%in%rd$RID)
	
	
	if(0)	#quick overview
	{
		group 		<- 'TYPE_PAIR_TODI2'
		#
		#	get couples that are most likely unlinked 
		#	
		rex			<- subset(rtp.tpairs, SELECT== 'couple most likely not a pair' & GROUP==group)[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)]), by=c('ID1','ID2')]
		rex			<- subset(rex, ID1!=ID2)
		rex			<- merge(rex, subset(rtp.tpairs, GROUP==group), by=c('ID1','ID2','PTY_RUN'))	
		rex			<- subset(rex, POSTERIOR_SCORE>=confidence.cut)		# sep16: 83; stage 2: 91 
		rex[, length(unique(c(ID1,ID2)))]								# sep16: 166; stage 2: 182 
		
		#
		#	likely transmission pairs, using topology and distance
		#	select one patient pairing across runs: that with lowest evidence
		rtp		<- subset(rtp.tpairs, grepl('most likely a pair',SELECT) & GROUP==group)[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)]), by=c('ID1','ID2')]
		rtp		<- subset(rtp, ID1!=ID2)
		rtp		<- merge(rtp, subset(rtp.tpairs, GROUP==group), by=c('ID1','ID2','PTY_RUN'))	
		rtp		<- subset(rtp, POSTERIOR_SCORE>=confidence.cut)		# 218 on couples run; stage1 614; stage2 125
		rtp[, length(unique(c(ID1,ID2)))]							# 366 on couples run; stage1 671; stage2 248
		stopifnot(!nrow(merge(subset(rtp, select=c(ID1,ID2)), subset(rex, select=c(ID1,ID2)), by=c('ID1','ID2'))))
		
		
		#	
		#	directed likely transmission pairs, using topology and distance	
		tmp		<- unique(subset(rtp, select=c('ID1','ID2','PTY_RUN')))		
		rtpd	<- merge(tmp, subset(rtp.tpairs, GROUP=='TYPE_DIRSCORE_TODI3'), by=c('ID1','ID2','PTY_RUN'))
		rtpd	<- subset(rtpd, POSTERIOR_SCORE>=confidence.cut)		# 136 on couples run; stage1 285; stage2 75
		rtpd[, length(unique(c(ID1,ID2)))]							# 237 on couples run; stage1 394; stage2 149
		
		#
		#	select only m->f, f->m
		#	
		rtpd[, table(ID1_SEX, ID2_SEX)]	
		# on couples run:
		#            ID2_SEX
		#	ID1_SEX  F  M
		#		  F 17 95
		#		  M 92 13		30/217= 13.8%
		# now on stage 2:
		#			ID2_SEX
		#ID1_SEX  F  M
		#		F  0 47
		#		M 28  0
	}
	
	#
	#	change to Male Female throughout
	#
	rtp.tpairs	<- subset(rtp.tpairs, ID1_SEX!=ID2_SEX)
	tmp			<- subset(rtp.tpairs, ID1_SEX=='M')
	setnames(tmp, colnames(tmp), gsub('ID1','MALE',colnames(tmp)))
	setnames(tmp, colnames(tmp), gsub('ID2','FEMALE',colnames(tmp)))
	set(tmp, tmp[, which(TYPE=='12')], 'TYPE', 'mf')
	set(tmp, tmp[, which(TYPE=='21')], 'TYPE', 'fm')
	rtp.tpairs	<- subset(rtp.tpairs, ID1_SEX=='F')
	setnames(rtp.tpairs, colnames(rtp.tpairs), gsub('ID1','FEMALE',colnames(rtp.tpairs)))
	setnames(rtp.tpairs, colnames(rtp.tpairs), gsub('ID2','MALE',colnames(rtp.tpairs)))
	set(rtp.tpairs, rtp.tpairs[, which(TYPE=='12')], 'TYPE', 'fm')
	set(rtp.tpairs, rtp.tpairs[, which(TYPE=='21')], 'TYPE', 'mf')
	rtp.tpairs	<- rbind(rtp.tpairs, tmp, use.names=TRUE)	 
	setnames(rtp.tpairs, c('MALE','FEMALE'), c('MALE_RID','FEMALE_RID'))
	set(rtp.tpairs, NULL, c('FEMALE_SEX','MALE_SEX'), NULL)
	tmp			<- copy(rplkl)
	setnames(tmp, colnames(tmp), gsub('21','FM',gsub('12','MF',gsub('ID2','FEMALE_RID',gsub('ID1','MALE_RID',colnames(tmp))))))	
	set(tmp, NULL, 'TYPE', tmp[, gsub('12','mf',TYPE)])
	set(tmp, NULL, 'TYPE', tmp[, gsub('21','fm',TYPE)])	
	setnames(rplkl, colnames(rplkl), gsub('21','MF',gsub('12','FM',gsub('ID2','MALE_RID',gsub('ID1','FEMALE_RID',colnames(rplkl))))))
	set(rplkl, NULL, 'TYPE', rplkl[, gsub('12','fm',TYPE)])
	set(rplkl, NULL, 'TYPE', rplkl[, gsub('21','mf',TYPE)])	
	rplkl		<- rbind(rplkl, tmp)
	rplkl		<- merge(rplkl, unique(subset(rp, select=c(MALE_RID, FEMALE_RID))), by=c('MALE_RID','FEMALE_RID'))	
	#rpw		<- melt(rpw, id.vars=c('PTY_RUN','ID1','ID2','W_FROM','W_TO','SUFFIX','TYPE_RAW','PATRISTIC_DISTANCE','ADJACENT','CONTIGUOUS','PATHS_12','PATHS_21','ID1_L','ID1_R','ID2_L','ID2_R'), variable.name='GROUP', value.name='TYPE')	
	tmp			<- copy(rpw)
	setnames(tmp, colnames(tmp), gsub('21','FM',gsub('12','MF',gsub('ID2','FEMALE',gsub('ID1','MALE',colnames(tmp))))))	
	set(tmp, NULL, 'TYPE', tmp[, gsub('12','mf',TYPE)])
	set(tmp, NULL, 'TYPE', tmp[, gsub('21','fm',TYPE)])	
	#set(tmp, NULL, 'TYPE_RAW', tmp[, gsub('12','mf',TYPE_RAW)])
	#set(tmp, NULL, 'TYPE_RAW', tmp[, gsub('21','fm',TYPE_RAW)])	
	setnames(rpw, colnames(rpw), gsub('21','MF',gsub('12','FM',gsub('ID2','MALE',gsub('ID1','FEMALE',colnames(rpw))))))
	set(rpw, NULL, 'TYPE', rpw[, gsub('12','fm',TYPE)])
	set(rpw, NULL, 'TYPE', rpw[, gsub('21','mf',TYPE)])	
	rpw			<- rbind(rpw, tmp)
	setnames(rpw, c('MALE','FEMALE'), c('MALE_RID','FEMALE_RID'))
	rpw			<- merge(rpw, unique(subset(rp, select=c(MALE_RID, FEMALE_RID))), by=c('MALE_RID','FEMALE_RID'))
	
	#
	#	make selections
	group 		<- 'TYPE_PAIR_TODI2'
	rex			<- subset(rtp.tpairs, SELECT=='couple most likely not a pair' & GROUP==group)[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)]), by=c('MALE_RID','FEMALE_RID')]	
	rex			<- merge(rex, subset(rtp.tpairs, GROUP==group), by=c('MALE_RID','FEMALE_RID','PTY_RUN'))	
	rex			<- subset(rex, POSTERIOR_SCORE>=confidence.cut)		 
	rtp			<- subset(rtp.tpairs, grepl('most likely a pair',SELECT) & GROUP==group)[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)]), by=c('MALE_RID','FEMALE_RID')]	
	rtp			<- merge(rtp, subset(rtp.tpairs, GROUP==group), by=c('MALE_RID','FEMALE_RID','PTY_RUN'))	
	rtp			<- subset(rtp, POSTERIOR_SCORE>=confidence.cut)		# 218 on couples run; stage1 614; stage2 125
	tmp			<- unique(subset(rtp, select=c('MALE_RID','FEMALE_RID','PTY_RUN')))		
	rtpd		<- merge(tmp, subset(rtp.tpairs, GROUP=='TYPE_DIRSCORE_TODI3'), by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	rtpd		<- subset(rtpd, POSTERIOR_SCORE>=confidence.cut)		
	rtpa		<- unique(subset(rplkl, select=c(FEMALE_RID, MALE_RID)))
	rtpa		<- subset(merge(rtpa, subset(rtp, select=c(FEMALE_RID, MALE_RID, PTY_RUN)), by=c('FEMALE_RID','MALE_RID'), all.x=1), is.na(PTY_RUN))
	rtpa		<- subset(merge(subset(rtpa, select=c(FEMALE_RID, MALE_RID)), subset(rex, select=c(FEMALE_RID, MALE_RID, PTY_RUN)), by=c('FEMALE_RID','MALE_RID'), all.x=1), is.na(PTY_RUN))
	rtpa		<- merge(subset(rtpa, select=c(FEMALE_RID, MALE_RID)), subset(rplkl, GROUP==group & TYPE=='likely pair'), by=c('MALE_RID','FEMALE_RID'))
	#rtpa[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
	rtpa[, POSTERIOR_SCORE:=(POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-2)]	
	rtpa		<- merge(rtpa, rtpa[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)]), by=c('MALE_RID','FEMALE_RID')], by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	rtpa[, SELECT:= 'couple ambiguous if pair or not pair']
	#
	#	make overall data.table
	#
	tmp		<- subset(rtpd, select=c(FEMALE_RID, MALE_RID))
	tmp[, DUMMY:='Y']
	tmp		<- merge(rtp, tmp, by=c('FEMALE_RID','MALE_RID'), all.x=1)
	tmp		<- subset(tmp, is.na(DUMMY))
	set(tmp, NULL, 'SELECT', 'couple most likely a pair direction not resolved')
	tmp[, DUMMY:=NULL]
	rca		<- rbind(rtpd, tmp, use.names=TRUE, fill=TRUE)
	rca		<- rbind(rca, rex, use.names=TRUE, fill=TRUE)
	rca		<- rbind(rca, rtpa, use.names=TRUE, fill=TRUE)
	#	add couples with insufficient sequence data for phyloscanner 	
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75.rda')
	tmp		<- copy(pty.runs)	
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75_part2.rda')
	tmp		<- rbind(tmp, pty.runs)
	set(tmp, NULL, c('BATCH','PTY_RUN'), NULL)
	tmp		<- unique(subset(tmp, select=RID))
	setnames(tmp, 'RID', 'MALE_RID')
	tmp2	<- merge(rp, tmp, by='MALE_RID')
	setnames(tmp, 'MALE_RID', 'FEMALE_RID')
	tmp2	<- merge(tmp2, tmp, by='FEMALE_RID')
	tmp2	<- unique(subset(tmp2, select=c(MALE_RID, FEMALE_RID)))
	tmp		<- merge(tmp2, unique(subset(rca, select=c('FEMALE_RID','MALE_RID','PTY_RUN'))), by=c('FEMALE_RID','MALE_RID'), all.x=1)	
	tmp		<- subset(tmp, is.na(PTY_RUN), c(MALE_RID, FEMALE_RID))
	tmp[, SELECT:='insufficient deep sequence data for at least one partner of couple']
	rca		<- rbind(rca, tmp, use.names=TRUE, fill=TRUE)
	#	if couples have less than 750 nt then call insufficient data:
	#	for our 250nt genomic windows, this is NEFF=3
	tmp		<- rca[, which(!grepl('insufficient',SELECT) & NEFF<3)]
	set(rca, tmp, 'SELECT', 'insufficient deep sequence data for at least one partner of couple')
	#
	#	add couple status
	tmp		<- unique(subset(rp, select=c(FEMALE_RID, MALE_RID, COUP_SC)))
	rca		<- merge(rca, tmp, by=c('FEMALE_RID','MALE_RID'), all.x=1)
	
	
	#
	#	determine first concordant pos visit
	#	add metadata at first concordant pos visit
	tmp		<- unique(subset(rd, select=c(RID, FIRSTPOSVIS, FIRSTPOSDATE)))
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	rca		<- merge(rca, tmp, by='MALE_RID')	
	setnames(tmp, colnames(tmp), gsub('MALE','FEMALE',colnames(tmp)))
	rca		<- merge(rca, tmp, by='FEMALE_RID')			
	rca[, VISIT_FIRSTCONCPOS:= rca[, pmax(MALE_FIRSTPOSVIS,FEMALE_FIRSTPOSVIS)]]
	set(rca, NULL, c('MALE_FIRSTPOSVIS','MALE_FIRSTPOSDATE','FEMALE_FIRSTPOSVIS','FEMALE_FIRSTPOSDATE'), NULL)
	# male rd
	tmp		<- unique(subset(rd, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'MALE_RID')
	tmp		<- unique(merge(rca, tmp, by='MALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('MALE_RID','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'MALE_RID', 'RID')
	tmp2	<- unique(subset(rd, select=c('RID','VISIT', "DATE", "REGION", "COMM_NUM", "COMM_NUM_A", "HH_NUM", "RECENTCD4","RECENTCD4DATE", "RECENTVL", "RECENTVLDATE", "SELFREPORTART", "EVERSELFREPORTART", "FIRSTSELFREPORTART","RELIGION","COHORT","BIRTHDATE","LASTNEGDATE","FIRSTPOSVIS","FIRSTPOSDATE","ARVSTARTDATE","EST_DATEDIED")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	setnames(tmp, 'MALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'MALE_VISIT', NULL)	
	rca		<- merge(rca, tmp, by=c('MALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	# male rh
	tmp		<- unique(subset(rh, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'MALE_RID')
	tmp		<- unique(merge(rca, tmp, by='MALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('MALE_RID','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'MALE_RID', 'RID')
	tmp2	<- unique(subset(rh, select=c("RID","VISIT","CIRCUM","MARSTAT","EDUCAT","RELCAT","OCAT","OCCUP_OLLI","SEXP1YR","SEXP1OUT","COMM_TYPE")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	setnames(tmp, 'MALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'MALE_VISIT', NULL)	
	rca		<- merge(rca, tmp, by=c('MALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	# female rd
	tmp		<- unique(subset(rd, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'FEMALE_RID')
	tmp		<- unique(merge(rca, tmp, by='FEMALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('FEMALE_RID','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'FEMALE_RID', 'RID')
	tmp2	<- unique(subset(rd, select=c('RID','VISIT', "DATE", "REGION", "COMM_NUM", "COMM_NUM_A", "HH_NUM", "RECENTCD4","RECENTCD4DATE", "RECENTVL", "RECENTVLDATE", "SELFREPORTART", "EVERSELFREPORTART", "FIRSTSELFREPORTART","RELIGION","COHORT","BIRTHDATE","LASTNEGDATE","FIRSTPOSVIS","FIRSTPOSDATE","ARVSTARTDATE","EST_DATEDIED")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('FEMALE_',colnames(tmp)))
	setnames(tmp, 'FEMALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'FEMALE_VISIT', NULL)	
	rca		<- merge(rca, tmp, by=c('FEMALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	# female rh
	tmp		<- unique(subset(rh, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'FEMALE_RID')
	tmp		<- unique(merge(rca, tmp, by='FEMALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('FEMALE_RID','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'FEMALE_RID', 'RID')
	tmp2	<- unique(subset(rh, select=c("RID","VISIT","CIRCUM","MARSTAT","EDUCAT","RELCAT","OCAT","OCCUP_OLLI","SEXP1YR","SEXP1OUT","COMM_TYPE")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('FEMALE_',colnames(tmp)))
	setnames(tmp, 'FEMALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'FEMALE_VISIT', NULL)	
	rca		<- merge(rca, tmp, by=c('FEMALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	
	
	#
	#	define pairs with consistent sero-history
	#
	rca[, SDC_TYPE:=NA_character_]
	set(rca, rca[, which(TYPE=='fm' & FEMALE_LASTNEGDATE>MALE_FIRSTPOSDATE)], 'SDC_TYPE', 'incorrect')
	set(rca, rca[, which(TYPE=='fm' & FEMALE_FIRSTPOSDATE<=MALE_LASTNEGDATE)], 'SDC_TYPE', 'correct')	
	set(rca, rca[, which(TYPE=='mf' & MALE_LASTNEGDATE>FEMALE_FIRSTPOSDATE)], 'SDC_TYPE', 'incorrect')
	set(rca, rca[, which(TYPE=='mf' & MALE_FIRSTPOSDATE<=FEMALE_LASTNEGDATE)], 'SDC_TYPE', 'correct')
	
	#
	#	plot windows of identified transmission pairs
	if(0)
	{
		rps			<- subset(rtp, select=c(MALE_RID, FEMALE_RID, PTY_RUN, TYPE, POSTERIOR_SCORE))
		setkey(rps, TYPE, POSTERIOR_SCORE)
		write.csv(rps, file=paste0(outfile.base,'_summary_lklpairs.csv'))
		rps[, DUMMY:=seq_len(nrow(rps))]
		rps[, LABEL:=rps[, factor(DUMMY, levels=DUMMY, labels=paste0('m ',MALE_RID,' f ', FEMALE_RID,'\n',TYPE,' ',round(POSTERIOR_SCORE, d=3),'\n',PTY_RUN))]]
		rpw[, ID_R_MAX:= pmax(FEMALE_R, MALE_R)]
		rpw[, ID_R_MIN:= pmin(FEMALE_R, MALE_R)]
		
		group		<- 'TYPE_BASIC'
		group		<- 'TYPE_PAIR_TODI'					
		rpw2		<- subset(rpw, GROUP==group)
		rplkl2		<- subset(rplkl, GROUP==group)	
		plot.file	<- paste0(outfile.base,'windows_summary_',group,'.pdf')
		setnames(rps, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rpw2, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rplkl2, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		phsc.plot.windowsummaries.for.pairs(rps, rpw2, rplkl2, plot.file, cols=NULL, group=group)				
	}
	if(0)
	{
		require(colorspace)
		zz	<- rtp[, which(MALE_RID%in%c('D030388','B035048'))]
		#for(ii in seq_len(nrow(rtpdm))[-1])
		for(ii in zz)
		{		
			indir		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun'
			# load dfr and phs
			load( file.path(indir, paste0('ptyr',rtp[ii,PTY_RUN],'_trees.rda')) )
			# setup plotting
			ids			<- c(rtp[ii, MALE_RID],rtp[ii, FEMALE_RID])
			dfs			<- subset(dfr, select=c(W_FROM, W_TO, IDX))
			dfs[, MALE_RID:=ids[1]]
			dfs[, FEMALE_RID:=ids[2]]				
			dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', rtp[ii, PTY_RUN], '\nwindow ', W_FROM,'-', W_TO, sep='')]]
			plot.file	<- paste0(outfile.base, 'run_', rtp[ii, PTY_RUN],'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
			invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.blacklisted=FALSE, drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10, tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'))						
		}
	}
	save(rca, rp, rd, rh, ra, rtp.todi2, rplkl, rpw, rtp, rex, rtpd, rtpa, file=outfile.save)
}

RakaiFull.preprocess.couples.todi.addingmetadata.170811<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	confidence.cut			<- 0.66	# do not change, because the prior is calibrated for 0.66
	#confidence.cut			<- 0.5
	neff.cut				<- 3
	infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23.rda"
	infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10.rda"
	infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10_zbl.rda"
	infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30.rda"
	infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_adj.rda"
	#infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min50.rda"
	#infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl25_prior23_min30.rda"
	#infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl45_prior23_min30.rda"	
	#infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl45_prior23_min30.rda"
	#infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_d30.rda"
	#infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_d100.rda"
	#infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_d1000.rda"
	#infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_s10.rda"
	#infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_s40.rda"
	infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_mf3.rda"
	infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_mf4.rda"
	infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_mf6.rda"
	infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_mt.rda"
	infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_prt.rda"
	infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_rg1.rda"
	infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_rg20.rda"
	infile.trmpairs.todi	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_zbl.rda"
	
	outfile.base			<- gsub('\\.rda','_',infile.trmpairs.todi) 
	if(confidence.cut!=0.66 | neff.cut!=3)
		outfile.base		<- paste0(outfile.base,'_pcut',100*confidence.cut,'_ncut',neff.cut,'_')
	outfile.save			<- paste0(outfile.base, 'withmetadata.rda')
	load(infile.trmpairs.todi)	
	rtp.tpairs	<- rtp.todi2
	rtp.tpairs	<- subset(rtp.tpairs, ID1%in%rd$RID & ID2%in%rd$RID)
	if(0)	#quick overview
	{
		group 		<- 'TYPE_PAIR_TODI2'
		#
		#	get couples that are most likely unlinked 
		#	
		rex			<- subset(rtp.tpairs, SELECT== 'couple most likely not a pair' & GROUP==group)[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)]), by=c('ID1','ID2')]
		rex			<- subset(rex, ID1!=ID2)
		rex			<- merge(rex, subset(rtp.tpairs, GROUP==group), by=c('ID1','ID2','PTY_RUN'))	
		rex			<- subset(rex, POSTERIOR_SCORE>=confidence.cut)		# sep16: 83; stage 2: 91 
		rex[, length(unique(c(ID1,ID2)))]								# sep16: 166; stage 2: 182 
		
		#
		#	likely transmission pairs, using topology and distance
		#	select one patient pairing across runs: that with lowest evidence
		rtp		<- subset(rtp.tpairs, grepl('most likely a pair',SELECT) & GROUP==group)[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)]), by=c('ID1','ID2')]
		rtp		<- subset(rtp, ID1!=ID2)
		rtp		<- merge(rtp, subset(rtp.tpairs, GROUP==group), by=c('ID1','ID2','PTY_RUN'))	
		rtp		<- subset(rtp, POSTERIOR_SCORE>=confidence.cut)		# 218 on couples run; stage1 614; stage2 125
		rtp[, length(unique(c(ID1,ID2)))]							# 366 on couples run; stage1 671; stage2 248
		stopifnot(!nrow(merge(subset(rtp, select=c(ID1,ID2)), subset(rex, select=c(ID1,ID2)), by=c('ID1','ID2'))))
		
		
		#	
		#	directed likely transmission pairs, using topology and distance	
		tmp		<- unique(subset(rtp, select=c('ID1','ID2','PTY_RUN')))		
		rtpd	<- merge(tmp, subset(rtp.tpairs, GROUP=='TYPE_DIRSCORE_TODI3'), by=c('ID1','ID2','PTY_RUN'))
		rtpd	<- subset(rtpd, POSTERIOR_SCORE>=confidence.cut)		# 136 on couples run; stage1 285; stage2 75
		rtpd[, length(unique(c(ID1,ID2)))]							# 237 on couples run; stage1 394; stage2 149
		
		#
		#	select only m->f, f->m
		#	
		rtpd[, table(ID1_SEX, ID2_SEX)]	
		# on couples run:
		#            ID2_SEX
		#	ID1_SEX  F  M
		#		  F 17 95
		#		  M 92 13		30/217= 13.8%
		# now on stage 2:
		#			ID2_SEX
		#ID1_SEX  F  M
		#		F  0 47
		#		M 28  0
	}
	
	#
	#	change to Male Female throughout
	#
	rtp.tpairs	<- subset(rtp.tpairs, ID1_SEX!=ID2_SEX)
	tmp			<- subset(rtp.tpairs, ID1_SEX=='M')
	setnames(tmp, colnames(tmp), gsub('ID1','MALE',colnames(tmp)))
	setnames(tmp, colnames(tmp), gsub('ID2','FEMALE',colnames(tmp)))
	set(tmp, tmp[, which(TYPE=='12')], 'TYPE', 'mf')
	set(tmp, tmp[, which(TYPE=='21')], 'TYPE', 'fm')
	rtp.tpairs	<- subset(rtp.tpairs, ID1_SEX=='F')
	setnames(rtp.tpairs, colnames(rtp.tpairs), gsub('ID1','FEMALE',colnames(rtp.tpairs)))
	setnames(rtp.tpairs, colnames(rtp.tpairs), gsub('ID2','MALE',colnames(rtp.tpairs)))
	set(rtp.tpairs, rtp.tpairs[, which(TYPE=='12')], 'TYPE', 'fm')
	set(rtp.tpairs, rtp.tpairs[, which(TYPE=='21')], 'TYPE', 'mf')
	rtp.tpairs	<- rbind(rtp.tpairs, tmp, use.names=TRUE)	 
	setnames(rtp.tpairs, c('MALE','FEMALE'), c('MALE_RID','FEMALE_RID'))
	set(rtp.tpairs, NULL, c('FEMALE_SEX','MALE_SEX'), NULL)
	tmp			<- copy(rplkl)
	setnames(tmp, colnames(tmp), gsub('21','FM',gsub('12','MF',gsub('ID2','FEMALE_RID',gsub('ID1','MALE_RID',colnames(tmp))))))	
	set(tmp, NULL, 'TYPE', tmp[, gsub('12','mf',TYPE)])
	set(tmp, NULL, 'TYPE', tmp[, gsub('21','fm',TYPE)])	
	setnames(rplkl, colnames(rplkl), gsub('21','MF',gsub('12','FM',gsub('ID2','MALE_RID',gsub('ID1','FEMALE_RID',colnames(rplkl))))))
	set(rplkl, NULL, 'TYPE', rplkl[, gsub('12','fm',TYPE)])
	set(rplkl, NULL, 'TYPE', rplkl[, gsub('21','mf',TYPE)])	
	rplkl		<- rbind(rplkl, tmp)
	rplkl		<- merge(rplkl, unique(subset(rp, select=c(MALE_RID, FEMALE_RID))), by=c('MALE_RID','FEMALE_RID'))	
	#rpw		<- melt(rpw, id.vars=c('PTY_RUN','ID1','ID2','W_FROM','W_TO','SUFFIX','TYPE_RAW','PATRISTIC_DISTANCE','ADJACENT','CONTIGUOUS','PATHS_12','PATHS_21','ID1_L','ID1_R','ID2_L','ID2_R'), variable.name='GROUP', value.name='TYPE')	
	tmp			<- copy(rpw)
	setnames(tmp, colnames(tmp), gsub('21','FM',gsub('12','MF',gsub('ID2','FEMALE',gsub('ID1','MALE',colnames(tmp))))))	
	set(tmp, NULL, 'TYPE', tmp[, gsub('12','mf',TYPE)])
	set(tmp, NULL, 'TYPE', tmp[, gsub('21','fm',TYPE)])	
	#set(tmp, NULL, 'TYPE_RAW', tmp[, gsub('12','mf',TYPE_RAW)])
	#set(tmp, NULL, 'TYPE_RAW', tmp[, gsub('21','fm',TYPE_RAW)])	
	setnames(rpw, colnames(rpw), gsub('21','MF',gsub('12','FM',gsub('ID2','MALE',gsub('ID1','FEMALE',colnames(rpw))))))
	set(rpw, NULL, 'TYPE', rpw[, gsub('12','fm',TYPE)])
	set(rpw, NULL, 'TYPE', rpw[, gsub('21','mf',TYPE)])	
	rpw			<- rbind(rpw, tmp)
	setnames(rpw, c('MALE','FEMALE'), c('MALE_RID','FEMALE_RID'))
	rpw			<- merge(rpw, unique(subset(rp, select=c(MALE_RID, FEMALE_RID))), by=c('MALE_RID','FEMALE_RID'))
	
	#
	#	select run with max NEFF for each pair
	#
	tmp			<- unique(subset(rtp.tpairs, GROUP=='TYPE_PAIR_TODI2', c(MALE_RID, FEMALE_RID, PTY_RUN, NEFF)))
	tmp			<- tmp[, list(PTY_RUN=PTY_RUN[which.max(NEFF)]), by=c('MALE_RID','FEMALE_RID')]
	rtp.tpairs	<- merge(rtp.tpairs, tmp, by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	rpw			<- merge(rpw, tmp, by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	rplkl		<- merge(rplkl, tmp, by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	#
	#	make selections
	group 		<- 'TYPE_PAIR_TODI2'
	rex			<- subset(rtp.tpairs, SELECT=='couple most likely not a pair' & GROUP==group)[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)]), by=c('MALE_RID','FEMALE_RID')]	
	rex			<- merge(rex, subset(rtp.tpairs, GROUP==group), by=c('MALE_RID','FEMALE_RID','PTY_RUN'))	
	rex			<- subset(rex, POSTERIOR_SCORE>=confidence.cut)		 
	rtp			<- subset(rtp.tpairs, grepl('most likely a pair',SELECT) & GROUP==group)[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)]), by=c('MALE_RID','FEMALE_RID')]	
	rtp			<- merge(rtp, subset(rtp.tpairs, GROUP==group), by=c('MALE_RID','FEMALE_RID','PTY_RUN'))	
	rtp			<- subset(rtp, POSTERIOR_SCORE>=confidence.cut)		# 218 on couples run; stage1 614; stage2 125
	tmp			<- unique(subset(rtp, select=c('MALE_RID','FEMALE_RID','PTY_RUN')))		
	rtpd		<- merge(tmp, subset(rtp.tpairs, GROUP=='TYPE_DIR_TODI2'), by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	rtpd		<- subset(rtpd, NEFF>=neff.cut & POSTERIOR_SCORE>=confidence.cut)		
	rtpa		<- unique(subset(rplkl, select=c(FEMALE_RID, MALE_RID)))
	rtpa		<- subset(merge(rtpa, subset(rtp, select=c(FEMALE_RID, MALE_RID, PTY_RUN)), by=c('FEMALE_RID','MALE_RID'), all.x=1), is.na(PTY_RUN))
	rtpa		<- subset(merge(subset(rtpa, select=c(FEMALE_RID, MALE_RID)), subset(rex, select=c(FEMALE_RID, MALE_RID, PTY_RUN)), by=c('FEMALE_RID','MALE_RID'), all.x=1), is.na(PTY_RUN))
	rtpa		<- merge(subset(rtpa, select=c(FEMALE_RID, MALE_RID)), subset(rplkl, GROUP==group & TYPE=='linked'), by=c('MALE_RID','FEMALE_RID'))
	#rtpa[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
	rtpa[, POSTERIOR_SCORE:=(POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-2)]	
	rtpa		<- merge(rtpa, rtpa[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)]), by=c('MALE_RID','FEMALE_RID')], by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	rtpa[, SELECT:= 'couple ambiguous if pair or not pair']
	#
	#	make overall data.table
	#
	tmp		<- subset(rtpd, select=c(FEMALE_RID, MALE_RID))
	tmp[, DUMMY:='Y']
	tmp		<- merge(rtp, tmp, by=c('FEMALE_RID','MALE_RID'), all.x=1)
	tmp		<- subset(tmp, is.na(DUMMY))
	set(tmp, NULL, 'SELECT', 'couple most likely a pair direction not resolved')
	tmp[, DUMMY:=NULL]
	rca		<- rbind(rtpd, tmp, use.names=TRUE, fill=TRUE)
	rca		<- rbind(rca, rex, use.names=TRUE, fill=TRUE)
	rca		<- rbind(rca, rtpa, use.names=TRUE, fill=TRUE)
	#
	#	set to ambiguous if NEFF not above cut
	#
	tmp		<- unique(subset(rplkl, GROUP=='TYPE_BASIC', c(MALE_RID, FEMALE_RID, PTY_RUN, NEFF)))
	setnames(tmp, 'NEFF', 'NEFF_TYPE_BASIC')
	rca		<- merge(rca, tmp, by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	set(rca, rca[, which(NEFF_TYPE_BASIC<=neff.cut)], 'SELECT', 'insufficient deep sequence data for at least one partner of couple')
	#
	#	add couples with insufficient sequence data for phyloscanner
	#
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75.rda')
	tmp		<- copy(pty.runs)	
	load('/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170704_b75_part2.rda')
	tmp		<- rbind(tmp, pty.runs)
	set(tmp, NULL, c('BATCH','PTY_RUN'), NULL)
	tmp		<- unique(subset(tmp, select=RID))
	setnames(tmp, 'RID', 'MALE_RID')
	tmp2	<- merge(rp, tmp, by='MALE_RID')
	setnames(tmp, 'MALE_RID', 'FEMALE_RID')
	tmp2	<- merge(tmp2, tmp, by='FEMALE_RID')
	tmp2	<- unique(subset(tmp2, select=c(MALE_RID, FEMALE_RID)))
	tmp		<- merge(tmp2, unique(subset(rca, select=c('FEMALE_RID','MALE_RID','PTY_RUN'))), by=c('FEMALE_RID','MALE_RID'), all.x=1)	
	tmp		<- subset(tmp, is.na(PTY_RUN), c(MALE_RID, FEMALE_RID))
	tmp[, SELECT:='insufficient deep sequence data for at least one partner of couple']
	rca		<- rbind(rca, tmp, use.names=TRUE, fill=TRUE)
	#
	#	add couple status
	tmp		<- unique(subset(rp, select=c(FEMALE_RID, MALE_RID, COUP_SC)))
	rca		<- merge(rca, tmp, by=c('FEMALE_RID','MALE_RID'), all.x=1)
	
	
	#
	#	determine first concordant pos visit
	#	add metadata at first concordant pos visit
	tmp		<- unique(subset(rd, select=c(RID, FIRSTPOSVIS, FIRSTPOSDATE)))
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	rca		<- merge(rca, tmp, by='MALE_RID')	
	setnames(tmp, colnames(tmp), gsub('MALE','FEMALE',colnames(tmp)))
	rca		<- merge(rca, tmp, by='FEMALE_RID')			
	rca[, VISIT_FIRSTCONCPOS:= rca[, pmax(MALE_FIRSTPOSVIS,FEMALE_FIRSTPOSVIS)]]
	set(rca, NULL, c('MALE_FIRSTPOSVIS','MALE_FIRSTPOSDATE','FEMALE_FIRSTPOSVIS','FEMALE_FIRSTPOSDATE'), NULL)
	# male rd
	tmp		<- unique(subset(rd, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'MALE_RID')
	tmp		<- unique(merge(rca, tmp, by='MALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('MALE_RID','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'MALE_RID', 'RID')
	tmp2	<- unique(subset(rd, select=c('RID','VISIT', "DATE", "REGION", "COMM_NUM", "COMM_NUM_A", "HH_NUM", "RECENTCD4","RECENTCD4DATE", "RECENTVL", "RECENTVLDATE", "SELFREPORTART", "EVERSELFREPORTART", "FIRSTSELFREPORTART","RELIGION","COHORT","BIRTHDATE","LASTNEGDATE","FIRSTPOSVIS","FIRSTPOSDATE","ARVSTARTDATE","EST_DATEDIED")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	setnames(tmp, 'MALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'MALE_VISIT', NULL)	
	rca		<- merge(rca, tmp, by=c('MALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	# male rh
	tmp		<- unique(subset(rh, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'MALE_RID')
	tmp		<- unique(merge(rca, tmp, by='MALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('MALE_RID','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'MALE_RID', 'RID')
	tmp2	<- unique(subset(rh, select=c("RID","VISIT","CIRCUM","MARSTAT","EDUCAT","RELCAT","OCAT","OCCUP_OLLI","SEXP1YR","SEXP1OUT","COMM_TYPE")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	setnames(tmp, 'MALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'MALE_VISIT', NULL)	
	rca		<- merge(rca, tmp, by=c('MALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	# female rd
	tmp		<- unique(subset(rd, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'FEMALE_RID')
	tmp		<- unique(merge(rca, tmp, by='FEMALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('FEMALE_RID','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'FEMALE_RID', 'RID')
	tmp2	<- unique(subset(rd, select=c('RID','VISIT', "DATE", "REGION", "COMM_NUM", "COMM_NUM_A", "HH_NUM", "RECENTCD4","RECENTCD4DATE", "RECENTVL", "RECENTVLDATE", "SELFREPORTART", "EVERSELFREPORTART", "FIRSTSELFREPORTART","RELIGION","COHORT","BIRTHDATE","LASTNEGDATE","FIRSTPOSVIS","FIRSTPOSDATE","ARVSTARTDATE","EST_DATEDIED")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('FEMALE_',colnames(tmp)))
	setnames(tmp, 'FEMALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'FEMALE_VISIT', NULL)	
	rca		<- merge(rca, tmp, by=c('FEMALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	# female rh
	tmp		<- unique(subset(rh, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'FEMALE_RID')
	tmp		<- unique(merge(rca, tmp, by='FEMALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('FEMALE_RID','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'FEMALE_RID', 'RID')
	tmp2	<- unique(subset(rh, select=c("RID","VISIT","CIRCUM","MARSTAT","EDUCAT","RELCAT","OCAT","OCCUP_OLLI","SEXP1YR","SEXP1OUT","COMM_TYPE")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('FEMALE_',colnames(tmp)))
	setnames(tmp, 'FEMALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'FEMALE_VISIT', NULL)	
	rca		<- merge(rca, tmp, by=c('FEMALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	
	
	#
	#	define pairs with consistent sero-history
	#
	rca[, SDC_TYPE:=NA_character_]
	set(rca, rca[, which(TYPE=='fm' & FEMALE_LASTNEGDATE>MALE_FIRSTPOSDATE)], 'SDC_TYPE', 'incorrect')
	set(rca, rca[, which(TYPE=='fm' & FEMALE_FIRSTPOSDATE<=MALE_LASTNEGDATE)], 'SDC_TYPE', 'correct')	
	set(rca, rca[, which(TYPE=='mf' & MALE_LASTNEGDATE>FEMALE_FIRSTPOSDATE)], 'SDC_TYPE', 'incorrect')
	set(rca, rca[, which(TYPE=='mf' & MALE_FIRSTPOSDATE<=FEMALE_LASTNEGDATE)], 'SDC_TYPE', 'correct')
	
	#
	#	plot windows of identified transmission pairs
	if(0)
	{
		rps			<- subset(rtp, select=c(MALE_RID, FEMALE_RID, PTY_RUN, TYPE, POSTERIOR_SCORE))
		setkey(rps, TYPE, POSTERIOR_SCORE)
		write.csv(rps, file=paste0(outfile.base,'_summary_lklpairs.csv'))
		rps[, DUMMY:=seq_len(nrow(rps))]
		rps[, LABEL:=rps[, factor(DUMMY, levels=DUMMY, labels=paste0('m ',MALE_RID,' f ', FEMALE_RID,'\n',TYPE,' ',round(POSTERIOR_SCORE, d=3),'\n',PTY_RUN))]]
		
		group		<- 'TYPE_BASIC'
		group		<- 'TYPE_PAIR_TODI2'					
		rpw2		<- subset(rpw, GROUP==group)
		rplkl2		<- subset(rplkl, GROUP==group)	
		plot.file	<- paste0(outfile.base,'windows_summary_',group,'.pdf')
		setnames(rps, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rpw2, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rplkl2, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		phsc.plot.windowsummaries.for.pairs(rps, rpw2, rplkl2, plot.file, cols=NULL, group=group)				
	}
	#
	#	plot windows of other couples
	if(0)
	{
		rps			<- subset(rca, SELECT%in%c('couple ambiguous if pair or not pair','couple most likely not a pair','insufficient deep sequence data for at least one partner of couple'), select=c(MALE_RID, FEMALE_RID, PTY_RUN, TYPE, POSTERIOR_SCORE))
		rps			<- subset(rps, !is.na(PTY_RUN))
		setkey(rps, TYPE, POSTERIOR_SCORE)
		write.csv(rps, file=paste0(outfile.base,'_summary_notlklpairs.csv'))
		rps[, DUMMY:=seq_len(nrow(rps))]
		rps[, LABEL:=rps[, factor(DUMMY, levels=DUMMY, labels=paste0('m ',MALE_RID,' f ', FEMALE_RID,'\n',TYPE,' ',round(POSTERIOR_SCORE, d=3),'\n',PTY_RUN))]]
		
		group		<- 'TYPE_BASIC'
		#group		<- 'TYPE_PAIR_TODI2'					
		rpw2		<- subset(rpw, GROUP==group)
		rplkl2		<- subset(rplkl, GROUP==group)	
		plot.file	<- paste0(outfile.base,'windows_summary_notlklpairs_',group,'.pdf')
		setnames(rps, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rpw2, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rplkl2, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		phsc.plot.windowsummaries.for.pairs(rps, rpw2, rplkl2, plot.file, cols=NULL, group=group)				
	}
	if(0)
	{
		require(colorspace)
		zz	<- rtp[, which(MALE_RID%in%c('D030388','B035048'))]
		#for(ii in seq_len(nrow(rtpdm))[-1])
		for(ii in zz)
		{		
			indir		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun'
			# load dfr and phs
			load( file.path(indir, paste0('ptyr',rtp[ii,PTY_RUN],'_trees.rda')) )
			# setup plotting
			ids			<- c(rtp[ii, MALE_RID],rtp[ii, FEMALE_RID])
			dfs			<- subset(dfr, select=c(W_FROM, W_TO, IDX))
			dfs[, MALE_RID:=ids[1]]
			dfs[, FEMALE_RID:=ids[2]]				
			dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', rtp[ii, PTY_RUN], '\nwindow ', W_FROM,'-', W_TO, sep='')]]
			plot.file	<- paste0(outfile.base, 'run_', rtp[ii, PTY_RUN],'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
			invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.blacklisted=FALSE, drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10, tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'))						
		}
	}
	save(rca, rp, rd, rh, ra, rtp.todi2, rplkl, rpw, rtp, rex, rtpd, rtpa, file=outfile.save)
}


RakaiFull.preprocess.trmpairs.todi.addingmetadata.170421<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
		
	confidence.cut			<- 0.5	
	confidence.cut			<- 0.66	
	neff.cut				<- 3
	group 					<- 'TYPE_PAIR_TODI2'
	#load rsmpl data.tables
	load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/community_sampling_170522.rda')
	#outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170428_"
	#infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_cl3.rda'		
	#outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_"
	infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170610_cl3_prior23.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170610_cl3_prior23_"	
	
	outfile.save			<- paste0(outfile.base, 'withmetadata.rda')
	#	now load second stage output
	load(infile.trmpairs.todi)
			
	#
	#	likely transmission pairs, using topology and distance
	#	select one patient pairing across runs: that with lowest evidence		
	rtp		<- subset(rtp.todi2, GROUP==group)[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)]), by=c('ID1','ID2')]
	rtp		<- subset(rtp, ID1!=ID2)
	rtp		<- merge(rtp, subset(rtp.todi2, GROUP==group), by=c('ID1','ID2','PTY_RUN'))	
	rtp		<- subset(rtp, POSTERIOR_SCORE>confidence.cut & NEFF>=neff.cut)		# first stage 2723, second stage 628, p mode 530	
	rtp[, length(unique(c(ID1,ID2)))]							# first stage 1511, second stage 1098, p mode 974 
	
	#	
	#	directed likely transmission pairs, using topology and distance	
	tmp		<- unique(subset(rtp, select=c('ID1','ID2','PTY_RUN')))		
	rtpd	<- merge(tmp, subset(rtp.todi2, GROUP=='TYPE_DIRSCORE_TODI3'), by=c('ID1','ID2','PTY_RUN'))
	rtpd	<- subset(rtpd, POSTERIOR_SCORE>confidence.cut & NEFF>=neff.cut)		# first stage 1536, second stage 364, p mode 260
	rtpd[, length(unique(c(ID1,ID2)))]							# first stage 936, second stage 668, p mode 492
	
	#
	#	select only m->f, f->m
	#	
	rtpd[, table(ID1_SEX, ID2_SEX)]	
	# on stage 1:
	#            ID2_SEX
	#	ID1_SEX  F  M
	#			F 507 496
	#			M 317 216		723/1511= 47%
	
	# on stage 2:
	#            ID2_SEX
	#	ID1_SEX  F  M
	#		  F  51 140
	#		  M 117  56			107/364= 29%

	# on stage 2 using 23 and neff>=3
	#	ID1_SEX   F   M    
    #  F  55 121
    #  M  95  44 				99/315= 31%

	rtpd.all<- copy(rtpd)	
	rtpd	<- subset(rtpd, ID1_SEX!=ID2_SEX)					
	#	first stage 824 pairs, second stage 257, pmode 192
	tmp		<- subset(rtpd, ID1_SEX=='M')
	setnames(tmp, colnames(tmp), gsub('ID1','MALE',colnames(tmp)))
	setnames(tmp, colnames(tmp), gsub('ID2','FEMALE',colnames(tmp)))
	set(tmp, tmp[, which(TYPE=='12')], 'TYPE', 'mf')
	set(tmp, tmp[, which(TYPE=='21')], 'TYPE', 'fm')
	rtpd	<- subset(rtpd, ID1_SEX=='F')
	setnames(rtpd, colnames(rtpd), gsub('ID1','FEMALE',colnames(rtpd)))
	setnames(rtpd, colnames(rtpd), gsub('ID2','MALE',colnames(rtpd)))
	set(rtpd, rtpd[, which(TYPE=='12')], 'TYPE', 'fm')
	set(rtpd, rtpd[, which(TYPE=='21')], 'TYPE', 'mf')
	rtpd	<- rbind(rtpd, tmp, use.names=TRUE)	 
	setnames(rtpd, c('MALE','FEMALE'), c('MALE_RID','FEMALE_RID'))
	set(rtpd, NULL, c('FEMALE_SEX','MALE_SEX'), NULL)
	#	same for rtp
	rtp.all	<- copy(rtp)	
	rtp	<- subset(rtp, ID1_SEX!=ID2_SEX)					
	tmp		<- subset(rtp, ID1_SEX=='M')
	setnames(tmp, colnames(tmp), gsub('ID1','MALE',colnames(tmp)))
	setnames(tmp, colnames(tmp), gsub('ID2','FEMALE',colnames(tmp)))
	set(tmp, tmp[, which(TYPE=='12')], 'TYPE', 'mf')
	set(tmp, tmp[, which(TYPE=='21')], 'TYPE', 'fm')
	rtp	<- subset(rtp, ID1_SEX=='F')
	setnames(rtp, colnames(rtp), gsub('ID1','FEMALE',colnames(rtp)))
	setnames(rtp, colnames(rtp), gsub('ID2','MALE',colnames(rtp)))
	set(rtp, rtp[, which(TYPE=='12')], 'TYPE', 'fm')
	set(rtp, rtp[, which(TYPE=='21')], 'TYPE', 'mf')
	rtp	<- rbind(rtp, tmp, use.names=TRUE)	 
	setnames(rtp, c('MALE','FEMALE'), c('MALE_RID','FEMALE_RID'))
	set(rtp, NULL, c('FEMALE_SEX','MALE_SEX'), NULL)
	
	#
	#	select only pairs with consistent sero-history
	#
	rtpd.hsx<- copy(rtpd)	
	tmp		<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
	stopifnot(!nrow(merge(subset(tmp[, length(BIRTHDATE), by='RID'], V1>1), tmp, by='RID')))	
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	rtpd	<- merge(rtpd, tmp, by='MALE_RID')	
	setnames(tmp, colnames(tmp), gsub('MALE','FEMALE',colnames(tmp)))
	rtpd	<- merge(rtpd, tmp, by='FEMALE_RID')			
	rtpd[, SDC_TYPE:=NA_character_]
	set(rtpd, rtpd[, which(TYPE=='fm' & FEMALE_LASTNEGDATE>MALE_FIRSTPOSDATE)], 'SDC_TYPE', 'incorrect')
	set(rtpd, rtpd[, which(TYPE=='fm' & FEMALE_FIRSTPOSDATE<=MALE_LASTNEGDATE)], 'SDC_TYPE', 'correct')	
	set(rtpd, rtpd[, which(TYPE=='mf' & MALE_LASTNEGDATE>FEMALE_FIRSTPOSDATE)], 'SDC_TYPE', 'incorrect')
	set(rtpd, rtpd[, which(TYPE=='mf' & MALE_FIRSTPOSDATE<=FEMALE_LASTNEGDATE)], 'SDC_TYPE', 'correct')
	subset(rtpd, !is.na(SDC_TYPE))[, table(TYPE, SDC_TYPE)]
	#	stage 1
	#		 correct 	incorrect
	#fm      29         2
	#mf      20         1				3/52= 0.0576
	#--> pretty good	
	
	#	stage 2
	#		 correct 	incorrect
	#fm      16         2
	#mf      12         1				3/31= 0.096
	#--> OK-ish	

	#	on stage 2 using 23 and neff>=3
	#TYPE correct incorrect
  	#fm      14         0
  	#mf       8         1				1/23= 0.0434

	#	plot
	if(0)
	{
		t.posneg	<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
		setnames(t.posneg, c('BIRTHDATE','EST_DATEDIED'), c('DOB','DOD'))
		t.seq		<- unique(subset(rs, !is.na(PID), select=c(RID, PID, SEQ_DATE)))
		t.cd4		<- unique(subset(rd, !is.na(RECENTCD4DATE) & !is.na(RECENTCD4), select=c(RID, RECENTCD4DATE, RECENTCD4)))
		set(t.cd4, NULL, 'RECENTCD4', t.cd4[, cut(RECENTCD4, breaks=c(-1,250,350,500,800,Inf), labels=c('<200','200-349','350-499','500-799','800+'))])
		t.vl		<- unique(subset(rd, !is.na(RECENTVLDATE) & !is.na(RECENTVL), select=c(RID, RECENTVLDATE, RECENTVL)))
		set(t.vl, NULL, 'RECENTVL', t.vl[, cut(RECENTVL, breaks=c(-1,200,1e3,1e5,Inf), labels=c('<200','200-1,000','1,000-100,000','100,000+'))])
		df			<- subset(rtpd, SDC_TYPE=='incorrect')
		p			<- RakaiFull.plot.epitimeline(df, copy(t.posneg), copy(t.cd4), copy(t.vl), copy(t.seq), age.adult=14)
		pdf(file=paste0(outfile.base,'epilines_pairs_with_inconsistent_serohistory.pdf'), w=8, h=5)
		print(p)
		dev.off()
		df			<- subset(rtpd, SDC_TYPE=='correct')
		p			<- RakaiFull.plot.epitimeline(df, copy(t.posneg), copy(t.cd4), copy(t.vl), copy(t.seq), age.adult=14)
		pdf(file=paste0(outfile.base,'epilines_pairs_with_consistent_serohistory.pdf'), w=8, h=20)
		print(p)
		dev.off()
	}
	rtpd	<- subset(rtpd, is.na(SDC_TYPE) | SDC_TYPE=='correct')
	rtpd	<- subset(rtpd, select=c(FEMALE_RID, MALE_RID, PTY_RUN, GROUP, TYPE, K, KEFF, N, NEFF, N_TYPE, PAR_PRIOR, POSTERIOR_ALPHA, POSTERIOR_BETA, POSTERIOR_SCORE))
	#
	#	check for multiple donors to same recipient
	#
	tmp	<- subset(rtpd, TYPE=='mf')[, list(FEMALE_N_TR=length(MALE_RID)), by='FEMALE_RID']
	tmp[, table(FEMALE_N_TR>1)]
	#stage 1
	#FALSE  TRUE 
	#  171   97 
	#stage 2
	#FALSE  TRUE 
	#147     1
	#on stage 2 using 23 and neff>=3
	#FALSE 
	#126
	tmp2<- subset(rtpd, TYPE=='fm')[, list(MALE_N_TR=length(FEMALE_RID)), by='MALE_RID']
	tmp2[, table(MALE_N_TR>1)]
	#stage 1
	#FALSE  TRUE 
	# 117    39
	#stage 2
	#FALSE 
	#105 
	#on stage 2 using 23 and neff>=3
	#FALSE TRUE
	#   87 1 
	#--> stage 2 is awesome
	rtpd	<- merge(rtpd, tmp, by='FEMALE_RID', all.x=1)
	rtpd	<- merge(rtpd, tmp2, by='MALE_RID', all.x=1)
	rtpd	<- subset(rtpd, MALE_N_TR==1 | FEMALE_N_TR==1)
	set(rtpd, NULL, c('FEMALE_N_TR','MALE_N_TR'), NULL)
	
	#
	#	add if in reported couple yes/no 
	#
	rtpd[, MALE_IN_COUPLE:= rtpd[, as.integer(MALE_RID%in%rp$MALE_RID)]]
	rtpd[, FEMALE_IN_COUPLE:= rtpd[, as.integer(FEMALE_RID%in%rp$FEMALE_RID)]]
	tmp		<- unique(subset(rp, select=c(FEMALE_RID, MALE_RID, COUP_SC)))
	rtpd	<- merge(rtpd, tmp, by=c('FEMALE_RID','MALE_RID'), all.x=1)
	setnames(rtpd, 'COUP_SC', 'COUPLE')
	tmp		<- rtpd[, which(!is.na(COUPLE))]
	set(rtpd, tmp, 'COUPLE', rtpd[tmp, paste('couple',COUPLE,'at enrollment')])
	set(rtpd, rtpd[, which(is.na(COUPLE))], 'COUPLE', 'no couple')
	rtpd[, table(COUPLE)]
	#	stage 1: 731 not within stable couples, 79 directed pairs among couples 
	#	stage 2: 167 not within stable couples, 87 directed pairs among couples			
	#	pmode:	 125 not within stable couples, 67 directed pairs among couples 
	#
	#	determine first concordant pos visit
	#	add metadata at first concordant pos visit
	tmp		<- unique(subset(rd, select=c(RID, FIRSTPOSVIS, FIRSTPOSDATE)))
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	rtpdm	<- merge(rtpd, tmp, by='MALE_RID')	
	setnames(tmp, colnames(tmp), gsub('MALE','FEMALE',colnames(tmp)))
	rtpdm	<- merge(rtpdm, tmp, by='FEMALE_RID')			
	rtpdm[, VISIT_FIRSTCONCPOS:= rtpdm[, pmax(MALE_FIRSTPOSVIS,FEMALE_FIRSTPOSVIS)]]
	set(rtpdm, NULL, c('MALE_FIRSTPOSVIS','MALE_FIRSTPOSDATE','FEMALE_FIRSTPOSVIS','FEMALE_FIRSTPOSDATE'), NULL)
	# male rd
	tmp		<- unique(subset(rd, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'MALE_RID')
	tmp		<- unique(merge(rtpdm, tmp, by='MALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('MALE_RID','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'MALE_RID', 'RID')
	tmp2	<- unique(subset(rd, select=c('RID','VISIT', "DATE", "REGION", "COMM_NUM", "COMM_NUM_A", "HH_NUM", "RECENTCD4","RECENTCD4DATE", "RECENTVL", "RECENTVLDATE", "SELFREPORTART", "EVERSELFREPORTART", "FIRSTSELFREPORTART","RELIGION","COHORT","BIRTHDATE","LASTNEGDATE","FIRSTPOSVIS","FIRSTPOSDATE","ARVSTARTDATE","EST_DATEDIED")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	setnames(tmp, 'MALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'MALE_VISIT', NULL)	
	rtpdm	<- merge(rtpdm, tmp, by=c('MALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	# male rh
	tmp		<- unique(subset(rh, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'MALE_RID')
	tmp		<- unique(merge(rtpdm, tmp, by='MALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('MALE_RID','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'MALE_RID', 'RID')
	tmp2	<- unique(subset(rh, select=c("RID","VISIT","CIRCUM","MARSTAT","EDUCAT","RELCAT","OCAT","OCCUP_OLLI","SEXP1YR","SEXP1OUT","COMM_TYPE")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	setnames(tmp, 'MALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'MALE_VISIT', NULL)	
	rtpdm	<- merge(rtpdm, tmp, by=c('MALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	# female rd
	tmp		<- unique(subset(rd, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'FEMALE_RID')
	tmp		<- unique(merge(rtpdm, tmp, by='FEMALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('FEMALE_RID','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'FEMALE_RID', 'RID')
	tmp2	<- unique(subset(rd, select=c('RID','VISIT', "DATE", "REGION", "COMM_NUM", "COMM_NUM_A", "HH_NUM", "RECENTCD4","RECENTCD4DATE", "RECENTVL", "RECENTVLDATE", "SELFREPORTART", "EVERSELFREPORTART", "FIRSTSELFREPORTART","RELIGION","COHORT","BIRTHDATE","LASTNEGDATE","FIRSTPOSVIS","FIRSTPOSDATE","ARVSTARTDATE","EST_DATEDIED")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('FEMALE_',colnames(tmp)))
	setnames(tmp, 'FEMALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'FEMALE_VISIT', NULL)	
	rtpdm	<- merge(rtpdm, tmp, by=c('FEMALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	# female rh
	tmp		<- unique(subset(rh, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'FEMALE_RID')
	tmp		<- unique(merge(rtpdm, tmp, by='FEMALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('FEMALE_RID','VISIT_FIRSTCONCPOS')])	
	setnames(tmp, 'FEMALE_RID', 'RID')
	tmp2	<- unique(subset(rh, select=c("RID","VISIT","CIRCUM","MARSTAT","EDUCAT","RELCAT","OCAT","OCCUP_OLLI","SEXP1YR","SEXP1OUT","COMM_TYPE")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('FEMALE_',colnames(tmp)))
	setnames(tmp, 'FEMALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'FEMALE_VISIT', NULL)	
	rtpdm	<- merge(rtpdm, tmp, by=c('FEMALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	# some clean up
	set(rtpdm, rtpdm[, which(FEMALE_COMM_NUM==23)], 'FEMALE_COMM_TYPE', 'agrarian')
	set(rtpdm, rtpdm[, which(MALE_COMM_NUM==23)], 'MALE_COMM_TYPE', 'agrarian')
	set(rtpdm, rtpdm[, which(FEMALE_COMM_NUM==38)], 'FEMALE_COMM_TYPE', 'fisherfolk')
	set(rtpdm, rtpdm[, which(MALE_COMM_NUM==38)], 'MALE_COMM_TYPE', 'fisherfolk')	
	
	#
	#	make sample data.table for geography
	rsm		<- rsmpl[, list(	PARTICIPATED=length(RID),
								COMM_NUM_A= COMM_NUM_A[1],
								COMM_TYPE= COMM_TYPE[1],
								HIV=sum(HIV),
								HAS_PID=sum(HAS_PID), 
								BAM_OUTPUT=sum(BAM_OUTPUT), 
								MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT)), by='COMM_NUM']	
	#ggplot(suppressWarnings(melt(rsm, id.vars=c('COMM_NUM','HIV'), measure.vars=c('HAS_PID','BAM_OUTPUT','MIN_PNG_OUTPUT'))), aes(x=COMM_NUM, y=value/HIV, colour=variable)) + geom_point()
	tmp		<- rsmpl3[, list(ELIGIBLE_AVG=mean(ELIGIBLE), PARTICIPATED_AVG=mean(PARTICIPATED)), by='COMM_NUM']	
	rsm		<- merge(rsm, tmp, all=1, by='COMM_NUM')
	#
	#	make sample data.table by age and gender and comm type
	rsmpl[, AGE_C:= rsmpl[, cut(AGEYRS, right=FALSE, breaks=c(seq(15,45,5),65), labels=paste0(seq(15,50,5)[-8], '-', seq(15,50,5)[-1]-1))]]	
	stopifnot( !nrow(subset(rsmpl, is.na(AGE_C))) )
	rsma	<- rsmpl[, list(	PARTICIPATED=length(RID),
								HIV=sum(HIV),
								HAS_PID=sum(HAS_PID), 
								BAM_OUTPUT=sum(BAM_OUTPUT), 
								MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT)), by=c('SEX','AGE_C','COMM_TYPE')]		
	#
	#	convert rpw into MALE FEMALE format needed for plotting
	rpw			<- melt(rpw, variable.name='GROUP', value.name='TYPE', measure.vars=c("TYPE_RAW","TYPE_BASIC","TYPE_PAIR_DI","TYPE_PAIRSCORE_DI","TYPE_PAIR_TO","TYPE_PAIR_TODI2x2","TYPE_DIR_TODI3","TYPE_DIRSCORE_TODI3","TYPE_PAIR_TODI2","TYPE_PAIR_TODI","TYPE_PAIRSCORE_TODI", "TYPE_CHAIN_TODI"))
	set(rpw, NULL, 'ID_R_MAX', rpw[, pmax(ID1_R,ID2_R)])
	set(rpw, NULL, 'ID_R_MIN', rpw[, pmin(ID1_R,ID2_R)])	
	tmp			<- copy(rpw)
	setnames(tmp, c('ID1','ID2','PATHS_12','PATHS_21','ID1_L','ID1_R','ID2_L','ID2_R'), c('ID2','ID1','PATHS_21','PATHS_12','ID2_L','ID2_R','ID1_L','ID1_R'))						
	set(tmp, NULL, 'TYPE', tmp[, gsub('12','XX33XX',TYPE)])						
	set(tmp, NULL, 'TYPE', tmp[, gsub('21','12',TYPE)])
	set(tmp, NULL, 'TYPE', tmp[, gsub('XX33XX','21',TYPE)])						
	rpw			<- rbind(tmp, rpw)
	setnames(rpw, c('ID1','ID2','PATHS_12','PATHS_21','ID1_L','ID1_R','ID2_L','ID2_R'), c('MALE_RID','FEMALE_RID','PATHS_MF','PATHS_FM','MALE_L','MALE_R','FEMALE_L','FEMALE_R'))
	set(rpw, NULL, 'TYPE', rpw[, gsub('21','fm',gsub('12','mf',TYPE))])
	rpw			<- merge(rpw, unique(subset(rtpdm, select=MALE_RID)), by='MALE_RID')
	rpw			<- merge(rpw, unique(subset(rtpdm, select=FEMALE_RID)), by='FEMALE_RID')
	#
	#	convert rplkl into MALE FEMALE format needed for plotting
	tmp			<- copy(rplkl)
	setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))						
	set(tmp, NULL, 'TYPE', tmp[, gsub('12','XX33XX',TYPE)])						
	set(tmp, NULL, 'TYPE', tmp[, gsub('21','12',TYPE)])
	set(tmp, NULL, 'TYPE', tmp[, gsub('XX33XX','21',TYPE)])						
	rplkl		<- rbind(tmp, rplkl)
	setnames(rplkl, c('ID1','ID2'), c('MALE_RID','FEMALE_RID'))
	set(rplkl, NULL, 'TYPE', rplkl[, gsub('21','fm',gsub('12','mf',TYPE))])
	rplkl		<- merge(rplkl, unique(subset(rtpdm, select=MALE_RID)), by='MALE_RID')
	rplkl		<- merge(rplkl, unique(subset(rtpdm, select=FEMALE_RID)), by='FEMALE_RID')
	
	#
	#	plot windows of identified transmission pairs
	if(0)
	{
		rps			<- subset(rtpdm, select=c(MALE_RID, FEMALE_RID, PTY_RUN, TYPE, POSTERIOR_SCORE))
		setkey(rps, TYPE, POSTERIOR_SCORE)
		write.csv(rps, file=paste0(outfile.base,'_summary_lklpairswdirection.csv'))
		rps[, DUMMY:=seq_len(nrow(rps))]
		rps[, LABEL:=rps[, factor(DUMMY, levels=DUMMY, labels=paste0('m ',MALE_RID,' f ', FEMALE_RID,'\n',TYPE,' ',round(POSTERIOR_SCORE, d=3),'\n',PTY_RUN))]]
		
		group		<- 'TYPE_BASIC'
		#group		<- 'TYPE_PAIR_TODI2'					
		rpw2		<- subset(rpw, GROUP==group)
		rplkl2		<- subset(rplkl, GROUP==group)	
		plot.file	<- paste0(outfile.base,'windows_summary_',group,'.pdf')
		setnames(rps, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rpw2, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rplkl2, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		phsc.plot.windowsummaries.for.pairs(rps, rpw2, rplkl2, plot.file, cols=NULL, group=group)				
	}	
	
	require(colorspace)
	#for(ii in seq_len(nrow(rtpdm))[-1])
	for(ii in 111:252)
	{		
		if(rtpdm[ii, PTY_RUN]!=1)
		{
			indir		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun'
			# load dfr and phs
			load( file.path(indir, paste0('ptyr',rtpdm[ii,PTY_RUN],'_trees.rda')) )
			# setup plotting
			ids			<- c(rtpdm[ii, MALE_RID],rtpdm[ii, FEMALE_RID])
			dfs			<- subset(dfr, select=c(W_FROM, W_TO, IDX))
			dfs[, MALE_RID:=ids[1]]
			dfs[, FEMALE_RID:=ids[2]]
			dfs			<- merge(dfs, subset(rpw, GROUP=='TYPE_RAW'), by=c('MALE_RID','FEMALE_RID','W_FROM','W_TO'), all.x=1)	
			dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', PTY_RUN, '\nwindow ', W_FROM,'-', W_TO,'\npmf',PATHS_MF,' pfm',PATHS_FM, ' ',ADJACENT,' ',CONTIGUOUS,' ',TYPE, '\n', round(PATRISTIC_DISTANCE, d=5), sep='')]]
			plot.file	<- paste0(outfile.base, 'trees/todi_pairs_170516_run_', rtpdm[ii, PTY_RUN],'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
			invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10, tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'))			
		}				
	}
	
	save(rtpdm, rp, rs, rd, rh, ra, rtp.todi2, rplkl, rpw, rtpd.hsx, rtp.all, rtp, rtpd, rsm, rsma, rsmpl, rsmpl2, rsmpl3, file=outfile.save)
}

RakaiFull.preprocess.trmpairs.todi.addingmetadata.170811<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
		
	confidence.cut			<- 0.66	
	neff.cut				<- 3
	group 					<- 'TYPE_PAIR_TODI2'
	#load rsmpl data.tables
	load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/community_sampling_170522.rda')
	#outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170428_"
	#infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_cl3.rda'		
	#outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_"
	infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_"	
	infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min10.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min10_"	
	infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_rg20.rda'
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_rg20_"
	infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_rg1.rda'
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_rg1_"
	infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_prt.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_prt_"	
	infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mt.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mt_"
	infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mf6.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mf6_"	
	infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mf4.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mf4_"	
	infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mf3.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mf3_"	
	infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_adj.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_adj_"	
	infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_zbl.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_zbl_"	
	
	outfile.save			<- paste0(outfile.base, 'withmetadata.rda')
	#	now load second stage output
	load(infile.trmpairs.todi)
	#
	#	above NEFF cut
	tmp			<- subset(rtp.todi2, GROUP=='TYPE_PAIR_TODI2' & TYPE=='linked' & NEFF>neff.cut, select=c(ID1, ID2, PTY_RUN))
	rtp.todi2	<- merge(rtp.todi2, tmp, by=c('ID1','ID2','PTY_RUN'))
	#
	#	likely transmission pairs, using topology and distance
	#	select one patient pairing across runs: that with lowest evidence		
	rtp		<- subset(rtp.todi2, GROUP==group)[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)]), by=c('ID1','ID2')]
	rtp		<- subset(rtp, ID1!=ID2)
	rtp		<- merge(rtp, subset(rtp.todi2, GROUP==group), by=c('ID1','ID2','PTY_RUN'))	
	rtp		<- subset(rtp, POSTERIOR_SCORE>confidence.cut)		# first stage 2723, second stage 628, p mode 530	
	rtp[, length(unique(c(ID1,ID2)))]							# first stage 1511, second stage 1098, p mode 974 
	
	#	
	#	directed likely transmission pairs, using topology and distance	
	tmp		<- unique(subset(rtp, select=c('ID1','ID2','PTY_RUN')))		
	rtpd	<- merge(tmp, subset(rtp.todi2, GROUP=='TYPE_DIR_TODI2'), by=c('ID1','ID2','PTY_RUN'))
	rtpd	<- subset(rtpd, POSTERIOR_SCORE>confidence.cut)		# first stage 1536, second stage 364, p mode 260
	rtpd[, length(unique(c(ID1,ID2)))]							# first stage 936, second stage 668, p mode 492
	
	#
	#	select only m->f, f->m
	#	
	rtpd[, table(ID1_SEX, ID2_SEX)]	
	# on stage 1:
	#            ID2_SEX
	#	ID1_SEX  F  M
	#			F 507 496
	#			M 317 216		723/1511= 47%
	
	# on stage 2:
	#            ID2_SEX
	#	ID1_SEX  F  M
	#		  F  51 140
	#		  M 117  56			107/364= 29%
	
	# on stage 2 using 23 and neff>=3
	#	ID1_SEX   F   M    
	#  F  55 121
	#  M  95  44 				99/315= 31%
	
	rtpd.all<- copy(rtpd)	
	rtpd	<- subset(rtpd, ID1_SEX!=ID2_SEX)					
	#	first stage 824 pairs, second stage 257, pmode 192
	tmp		<- subset(rtpd, ID1_SEX=='M')
	setnames(tmp, colnames(tmp), gsub('ID1','MALE',colnames(tmp)))
	setnames(tmp, colnames(tmp), gsub('ID2','FEMALE',colnames(tmp)))
	set(tmp, tmp[, which(TYPE=='12')], 'TYPE', 'mf')
	set(tmp, tmp[, which(TYPE=='21')], 'TYPE', 'fm')
	rtpd	<- subset(rtpd, ID1_SEX=='F')
	setnames(rtpd, colnames(rtpd), gsub('ID1','FEMALE',colnames(rtpd)))
	setnames(rtpd, colnames(rtpd), gsub('ID2','MALE',colnames(rtpd)))
	set(rtpd, rtpd[, which(TYPE=='12')], 'TYPE', 'fm')
	set(rtpd, rtpd[, which(TYPE=='21')], 'TYPE', 'mf')
	rtpd	<- rbind(rtpd, tmp, use.names=TRUE)	 
	setnames(rtpd, c('MALE','FEMALE'), c('MALE_RID','FEMALE_RID'))
	set(rtpd, NULL, c('FEMALE_SEX','MALE_SEX'), NULL)
	#	same for rtp
	rtp.all	<- copy(rtp)	
	rtp	<- subset(rtp, ID1_SEX!=ID2_SEX)					
	tmp		<- subset(rtp, ID1_SEX=='M')
	setnames(tmp, colnames(tmp), gsub('ID1','MALE',colnames(tmp)))
	setnames(tmp, colnames(tmp), gsub('ID2','FEMALE',colnames(tmp)))
	set(tmp, tmp[, which(TYPE=='12')], 'TYPE', 'mf')
	set(tmp, tmp[, which(TYPE=='21')], 'TYPE', 'fm')
	rtp	<- subset(rtp, ID1_SEX=='F')
	setnames(rtp, colnames(rtp), gsub('ID1','FEMALE',colnames(rtp)))
	setnames(rtp, colnames(rtp), gsub('ID2','MALE',colnames(rtp)))
	set(rtp, rtp[, which(TYPE=='12')], 'TYPE', 'fm')
	set(rtp, rtp[, which(TYPE=='21')], 'TYPE', 'mf')
	rtp	<- rbind(rtp, tmp, use.names=TRUE)	 
	setnames(rtp, c('MALE','FEMALE'), c('MALE_RID','FEMALE_RID'))
	set(rtp, NULL, c('FEMALE_SEX','MALE_SEX'), NULL)
	
	#
	#	select only pairs with consistent sero-history
	#
	rtpd.hsx<- copy(rtpd)	
	tmp		<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
	stopifnot(!nrow(merge(subset(tmp[, length(BIRTHDATE), by='RID'], V1>1), tmp, by='RID')))	
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	rtpd	<- merge(rtpd, tmp, by='MALE_RID')	
	setnames(tmp, colnames(tmp), gsub('MALE','FEMALE',colnames(tmp)))
	rtpd	<- merge(rtpd, tmp, by='FEMALE_RID')			
	rtpd[, SDC_TYPE:=NA_character_]
	set(rtpd, rtpd[, which(TYPE=='fm' & FEMALE_LASTNEGDATE>MALE_FIRSTPOSDATE)], 'SDC_TYPE', 'incorrect')
	set(rtpd, rtpd[, which(TYPE=='fm' & FEMALE_FIRSTPOSDATE<=MALE_LASTNEGDATE)], 'SDC_TYPE', 'correct')	
	set(rtpd, rtpd[, which(TYPE=='mf' & MALE_LASTNEGDATE>FEMALE_FIRSTPOSDATE)], 'SDC_TYPE', 'incorrect')
	set(rtpd, rtpd[, which(TYPE=='mf' & MALE_FIRSTPOSDATE<=FEMALE_LASTNEGDATE)], 'SDC_TYPE', 'correct')
	subset(rtpd, !is.na(SDC_TYPE))[, table(TYPE, SDC_TYPE)]
	#	stage 1
	#		 correct 	incorrect
	#fm      29         2
	#mf      20         1				3/52= 0.0576
	#--> pretty good	
	
	#	stage 2
	#		 correct 	incorrect
	#fm      16         2
	#mf      12         1				3/31= 0.096
	#--> OK-ish	
	
	#	on stage 2 using 23 and neff>=3
	#TYPE correct incorrect
	#fm      14         0
	#mf       8         1				1/23= 0.0434
	
	#	plot
	if(0)
	{
		t.posneg	<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
		setnames(t.posneg, c('BIRTHDATE','EST_DATEDIED'), c('DOB','DOD'))
		t.seq		<- unique(subset(rs, !is.na(PID), select=c(RID, PID, SEQ_DATE)))
		t.cd4		<- unique(subset(rd, !is.na(RECENTCD4DATE) & !is.na(RECENTCD4), select=c(RID, RECENTCD4DATE, RECENTCD4)))
		set(t.cd4, NULL, 'RECENTCD4', t.cd4[, cut(RECENTCD4, breaks=c(-1,250,350,500,800,Inf), labels=c('<200','200-349','350-499','500-799','800+'))])
		t.vl		<- unique(subset(rd, !is.na(RECENTVLDATE) & !is.na(RECENTVL), select=c(RID, RECENTVLDATE, RECENTVL)))
		set(t.vl, NULL, 'RECENTVL', t.vl[, cut(RECENTVL, breaks=c(-1,200,1e3,1e5,Inf), labels=c('<200','200-1,000','1,000-100,000','100,000+'))])
		df			<- subset(rtpd, SDC_TYPE=='incorrect')
		p			<- RakaiFull.plot.epitimeline(df, copy(t.posneg), copy(t.cd4), copy(t.vl), copy(t.seq), age.adult=14)
		pdf(file=paste0(outfile.base,'epilines_pairs_with_inconsistent_serohistory.pdf'), w=8, h=5)
		print(p)
		dev.off()
		df			<- subset(rtpd, SDC_TYPE=='correct')
		p			<- RakaiFull.plot.epitimeline(df, copy(t.posneg), copy(t.cd4), copy(t.vl), copy(t.seq), age.adult=14)
		pdf(file=paste0(outfile.base,'epilines_pairs_with_consistent_serohistory.pdf'), w=8, h=20)
		print(p)
		dev.off()
	}
	rtpd	<- subset(rtpd, is.na(SDC_TYPE) | SDC_TYPE=='correct')
	rtpd	<- subset(rtpd, select=c(FEMALE_RID, MALE_RID, PTY_RUN, GROUP, TYPE, K, KEFF, N, NEFF, N_TYPE, PAR_PRIOR, POSTERIOR_ALPHA, POSTERIOR_BETA, POSTERIOR_SCORE))
	#
	#	check for multiple donors to same recipient
	#
	tmp	<- subset(rtpd, TYPE=='mf')[, list(FEMALE_N_TR=length(MALE_RID)), by='FEMALE_RID']
	tmp[, table(FEMALE_N_TR>1)]
	#stage 1
	#FALSE  TRUE 
	#  171   97 
	#stage 2
	#FALSE  TRUE 
	#147     1
	#on stage 2 using 23 and neff>=3
	#FALSE 
	#126
	tmp2<- subset(rtpd, TYPE=='fm')[, list(MALE_N_TR=length(FEMALE_RID)), by='MALE_RID']
	tmp2[, table(MALE_N_TR>1)]
	#stage 1
	#FALSE  TRUE 
	# 117    39
	#stage 2
	#FALSE 
	#105 
	#on stage 2 using 23 and neff>=3
	#FALSE TRUE
	#   87 1 
	#--> stage 2 is awesome
	rtpd	<- merge(rtpd, tmp, by='FEMALE_RID', all.x=1)
	rtpd	<- merge(rtpd, tmp2, by='MALE_RID', all.x=1)
	rtpd	<- subset(rtpd, MALE_N_TR==1 | FEMALE_N_TR==1)
	set(rtpd, NULL, c('FEMALE_N_TR','MALE_N_TR'), NULL)
	
	#
	#	add if in reported couple yes/no 
	#
	rtpd[, MALE_IN_COUPLE:= rtpd[, as.integer(MALE_RID%in%rp$MALE_RID)]]
	rtpd[, FEMALE_IN_COUPLE:= rtpd[, as.integer(FEMALE_RID%in%rp$FEMALE_RID)]]
	tmp		<- unique(subset(rp, select=c(FEMALE_RID, MALE_RID, COUP_SC)))
	rtpd	<- merge(rtpd, tmp, by=c('FEMALE_RID','MALE_RID'), all.x=1)
	setnames(rtpd, 'COUP_SC', 'COUPLE')
	tmp		<- rtpd[, which(!is.na(COUPLE))]
	set(rtpd, tmp, 'COUPLE', rtpd[tmp, paste('couple',COUPLE,'at enrollment')])
	set(rtpd, rtpd[, which(is.na(COUPLE))], 'COUPLE', 'no couple')
	rtpd[, table(COUPLE)]
	#	stage 1: 731 not within stable couples, 79 directed pairs among couples 
	#	stage 2: 167 not within stable couples, 87 directed pairs among couples			
	#	pmode:	 125 not within stable couples, 67 directed pairs among couples 
	#
	#	determine first concordant pos visit
	#	add metadata at first concordant pos visit
	tmp		<- unique(subset(rd, select=c(RID, FIRSTPOSVIS, FIRSTPOSDATE)))
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	rtpdm	<- merge(rtpd, tmp, by='MALE_RID')	
	setnames(tmp, colnames(tmp), gsub('MALE','FEMALE',colnames(tmp)))
	rtpdm	<- merge(rtpdm, tmp, by='FEMALE_RID')			
	rtpdm[, VISIT_FIRSTCONCPOS:= rtpdm[, pmax(MALE_FIRSTPOSVIS,FEMALE_FIRSTPOSVIS)]]
	set(rtpdm, NULL, c('MALE_FIRSTPOSVIS','MALE_FIRSTPOSDATE','FEMALE_FIRSTPOSVIS','FEMALE_FIRSTPOSDATE'), NULL)
	# male rd
	tmp		<- unique(subset(rd, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'MALE_RID')
	tmp		<- unique(merge(rtpdm, tmp, by='MALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('MALE_RID','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'MALE_RID', 'RID')
	tmp2	<- unique(subset(rd, select=c('RID','VISIT', "DATE", "REGION", "COMM_NUM", "COMM_NUM_A", "HH_NUM", "RECENTCD4","RECENTCD4DATE", "RECENTVL", "RECENTVLDATE", "SELFREPORTART", "EVERSELFREPORTART", "FIRSTSELFREPORTART","RELIGION","COHORT","BIRTHDATE","LASTNEGDATE","FIRSTPOSVIS","FIRSTPOSDATE","ARVSTARTDATE","EST_DATEDIED")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	setnames(tmp, 'MALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'MALE_VISIT', NULL)	
	rtpdm	<- merge(rtpdm, tmp, by=c('MALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	# male rh
	tmp		<- unique(subset(rh, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'MALE_RID')
	tmp		<- unique(merge(rtpdm, tmp, by='MALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('MALE_RID','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'MALE_RID', 'RID')
	tmp2	<- unique(subset(rh, select=c("RID","VISIT","CIRCUM","MARSTAT","EDUCAT","RELCAT","OCAT","OCCUP_OLLI","SEXP1YR","SEXP1OUT","COMM_TYPE")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	setnames(tmp, 'MALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'MALE_VISIT', NULL)	
	rtpdm	<- merge(rtpdm, tmp, by=c('MALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	# female rd
	tmp		<- unique(subset(rd, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'FEMALE_RID')
	tmp		<- unique(merge(rtpdm, tmp, by='FEMALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('FEMALE_RID','VISIT_FIRSTCONCPOS')])
	setnames(tmp, 'FEMALE_RID', 'RID')
	tmp2	<- unique(subset(rd, select=c('RID','VISIT', "DATE", "REGION", "COMM_NUM", "COMM_NUM_A", "HH_NUM", "RECENTCD4","RECENTCD4DATE", "RECENTVL", "RECENTVLDATE", "SELFREPORTART", "EVERSELFREPORTART", "FIRSTSELFREPORTART","RELIGION","COHORT","BIRTHDATE","LASTNEGDATE","FIRSTPOSVIS","FIRSTPOSDATE","ARVSTARTDATE","EST_DATEDIED")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('FEMALE_',colnames(tmp)))
	setnames(tmp, 'FEMALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'FEMALE_VISIT', NULL)	
	rtpdm	<- merge(rtpdm, tmp, by=c('FEMALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	# female rh
	tmp		<- unique(subset(rh, select=c(RID, VISIT)))
	setnames(tmp, 'RID', 'FEMALE_RID')
	tmp		<- unique(merge(rtpdm, tmp, by='FEMALE_RID')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('FEMALE_RID','VISIT_FIRSTCONCPOS')])	
	setnames(tmp, 'FEMALE_RID', 'RID')
	tmp2	<- unique(subset(rh, select=c("RID","VISIT","CIRCUM","MARSTAT","EDUCAT","RELCAT","OCAT","OCCUP_OLLI","SEXP1YR","SEXP1OUT","COMM_TYPE")), by=c('RID','VISIT'))
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
	setnames(tmp, colnames(tmp), paste0('FEMALE_',colnames(tmp)))
	setnames(tmp, 'FEMALE_VISIT_FIRSTCONCPOS', 'VISIT_FIRSTCONCPOS')
	set(tmp, NULL, 'FEMALE_VISIT', NULL)	
	rtpdm	<- merge(rtpdm, tmp, by=c('FEMALE_RID','VISIT_FIRSTCONCPOS'), all.x=1)
	# some clean up
	set(rtpdm, rtpdm[, which(FEMALE_COMM_NUM==23)], 'FEMALE_COMM_TYPE', 'agrarian')
	set(rtpdm, rtpdm[, which(MALE_COMM_NUM==23)], 'MALE_COMM_TYPE', 'agrarian')
	set(rtpdm, rtpdm[, which(FEMALE_COMM_NUM==38)], 'FEMALE_COMM_TYPE', 'fisherfolk')
	set(rtpdm, rtpdm[, which(MALE_COMM_NUM==38)], 'MALE_COMM_TYPE', 'fisherfolk')	
	
	#
	#	make sample data.table for geography
	rsm		<- rsmpl[, list(	PARTICIPATED=length(RID),
					COMM_NUM_A= COMM_NUM_A[1],
					COMM_TYPE= COMM_TYPE[1],
					HIV=sum(HIV),
					HAS_PID=sum(HAS_PID), 
					BAM_OUTPUT=sum(BAM_OUTPUT), 
					MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT)), by='COMM_NUM']	
	#ggplot(suppressWarnings(melt(rsm, id.vars=c('COMM_NUM','HIV'), measure.vars=c('HAS_PID','BAM_OUTPUT','MIN_PNG_OUTPUT'))), aes(x=COMM_NUM, y=value/HIV, colour=variable)) + geom_point()
	tmp		<- rsmpl3[, list(ELIGIBLE_AVG=mean(ELIGIBLE), PARTICIPATED_AVG=mean(PARTICIPATED)), by='COMM_NUM']	
	rsm		<- merge(rsm, tmp, all=1, by='COMM_NUM')
	#
	#	make sample data.table by age and gender and comm type
	rsmpl[, AGE_C:= rsmpl[, cut(AGEYRS, right=FALSE, breaks=c(seq(15,45,5),65), labels=paste0(seq(15,50,5)[-8], '-', seq(15,50,5)[-1]-1))]]	
	stopifnot( !nrow(subset(rsmpl, is.na(AGE_C))) )
	rsma	<- rsmpl[, list(	PARTICIPATED=length(RID),
					HIV=sum(HIV),
					HAS_PID=sum(HAS_PID), 
					BAM_OUTPUT=sum(BAM_OUTPUT), 
					MIN_PNG_OUTPUT=sum(MIN_PNG_OUTPUT)), by=c('SEX','AGE_C','COMM_TYPE')]		
	#
	#	convert rpw into MALE FEMALE format needed for plotting
	tmp			<- copy(rpw)
	setnames(tmp, c('ID1','ID2','PATHS_12','PATHS_21','ID1_L','ID1_R','ID2_L','ID2_R'), c('ID2','ID1','PATHS_21','PATHS_12','ID2_L','ID2_R','ID1_L','ID1_R'))						
	set(tmp, NULL, 'TYPE', tmp[, gsub('12','XX33XX',TYPE)])						
	set(tmp, NULL, 'TYPE', tmp[, gsub('21','12',TYPE)])
	set(tmp, NULL, 'TYPE', tmp[, gsub('XX33XX','21',TYPE)])						
	rpw			<- rbind(tmp, rpw)
	setnames(rpw, c('ID1','ID2','PATHS_12','PATHS_21','ID1_L','ID1_R','ID2_L','ID2_R'), c('MALE_RID','FEMALE_RID','PATHS_MF','PATHS_FM','MALE_L','MALE_R','FEMALE_L','FEMALE_R'))
	set(rpw, NULL, 'TYPE', rpw[, gsub('21','fm',gsub('12','mf',TYPE))])
	rpw			<- merge(rpw, unique(subset(rtpdm, select=MALE_RID)), by='MALE_RID')
	rpw			<- merge(rpw, unique(subset(rtpdm, select=FEMALE_RID)), by='FEMALE_RID')
	#
	#	convert rplkl into MALE FEMALE format needed for plotting
	tmp			<- copy(rplkl)
	setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))						
	set(tmp, NULL, 'TYPE', tmp[, gsub('12','XX33XX',TYPE)])						
	set(tmp, NULL, 'TYPE', tmp[, gsub('21','12',TYPE)])
	set(tmp, NULL, 'TYPE', tmp[, gsub('XX33XX','21',TYPE)])						
	rplkl		<- rbind(tmp, rplkl)
	setnames(rplkl, c('ID1','ID2'), c('MALE_RID','FEMALE_RID'))
	set(rplkl, NULL, 'TYPE', rplkl[, gsub('21','fm',gsub('12','mf',TYPE))])
	rplkl		<- merge(rplkl, unique(subset(rtpdm, select=MALE_RID)), by='MALE_RID')
	rplkl		<- merge(rplkl, unique(subset(rtpdm, select=FEMALE_RID)), by='FEMALE_RID')
	
	#
	#	plot windows of identified transmission pairs
	if(0)
	{
		rps			<- subset(rtpdm, select=c(MALE_RID, FEMALE_RID, PTY_RUN, TYPE, POSTERIOR_SCORE))
		setkey(rps, TYPE, POSTERIOR_SCORE)
		write.csv(rps, file=paste0(outfile.base,'_summary_lklpairswdirection.csv'))
		rps[, DUMMY:=seq_len(nrow(rps))]
		rps[, LABEL:=rps[, factor(DUMMY, levels=DUMMY, labels=paste0('m ',MALE_RID,' f ', FEMALE_RID,'\n',TYPE,' ',round(POSTERIOR_SCORE, d=3),'\n',PTY_RUN))]]
		
		group		<- 'TYPE_BASIC'
		#group		<- 'TYPE_PAIR_TODI2'					
		rpw2		<- subset(rpw, GROUP==group)
		rplkl2		<- subset(rplkl, GROUP==group)	
		plot.file	<- paste0(outfile.base,'windows_summary_',group,'.pdf')
		setnames(rps, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rpw2, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rplkl2, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		phsc.plot.windowsummaries.for.pairs(rps, rpw2, rplkl2, plot.file, cols=NULL, group=group)				
	}	
	if(0)
	{
		require(colorspace)
		#for(ii in seq_len(nrow(rtpdm))[-1])
		for(ii in 111:252)
		{		
			if(rtpdm[ii, PTY_RUN]!=1)
			{
				indir		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun'
				# load dfr and phs
				load( file.path(indir, paste0('ptyr',rtpdm[ii,PTY_RUN],'_trees.rda')) )
				# setup plotting
				ids			<- c(rtpdm[ii, MALE_RID],rtpdm[ii, FEMALE_RID])
				dfs			<- subset(dfr, select=c(W_FROM, W_TO, IDX))
				dfs[, MALE_RID:=ids[1]]
				dfs[, FEMALE_RID:=ids[2]]
				dfs			<- merge(dfs, subset(rpw, GROUP=='TYPE_RAW'), by=c('MALE_RID','FEMALE_RID','W_FROM','W_TO'), all.x=1)	
				dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', PTY_RUN, '\nwindow ', W_FROM,'-', W_TO,'\npmf',PATHS_MF,' pfm',PATHS_FM, ' ',ADJACENT,' ',CONTIGUOUS,' ',TYPE, '\n', round(PATRISTIC_DISTANCE, d=5), sep='')]]
				plot.file	<- paste0(outfile.base, 'trees/todi_pairs_170516_run_', rtpdm[ii, PTY_RUN],'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
				invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10, tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'))			
			}				
		}
	}	
	
	save(rtpdm, rp, rs, rd, rh, ra, rtp.todi2, rplkl, rpw, rtpd.hsx, rtpd.all, rtp, rtpd, rsm, rsma, rsmpl, rsmpl2, rsmpl3, file=outfile.save)
}

RakaiFull.analyze.trmpairs.community.anonymize.170522<- function()
{
	#
	#	anonymize IDs
	#
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv'	
	dc		<- data.table(	COMM_NUM_RAW=	c("1","2","3","4","5","6","7","8","9","14","15","16","18","19","22","23","24","25","29","32","33","34","35","36","38","40","44","45","46","51","52","53","54","55", "56","57","58","59","60","61","62","65","67","74","77","81","84","89","94","95","103","106","107","108","109","120","177", "183", "256", "370","391","401","451", "468","602", "754", "755", "760", "770","771","772","773","774","776"),
							COMM_TYPE=		c("T","A","A","T","A","A","A","A","A", "A", "A", "T", "A", "A", "T", "A", "T", "A", "A", "A", "T", "A", "A", "A", "F", "A", "A", "A", "A", "T", "A", "A", "A", "A",  "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",  "A","A",   "A",  "T",  "A",  "A",  "A",  "A",   "A",   "A",   "A",  "A", "A",  "A",    "A",  "A",  "A",    "A",   "A", "F",  "F",  "A",  "A",  "F",   "T"))
	set(dc, NULL, 'COMM_TYPE', dc[, as.character(factor(COMM_TYPE, levels=c('A','T','F'), labels=c('agrarian','trading','fisherfolk')))])		
	set(dc, NULL, 'COMM_NUM', dc[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',COMM_NUM_RAW))))])
	
	seed	<- 42L 
	set.seed(seed)
	tmp		<- unique(dc, by='COMM_NUM')
	tmp		<- tmp[, list( 	COMM_NUM=COMM_NUM,
							COMM_NUM_A=paste0(substring(COMM_TYPE,1,1), letters[round(runif(length(COMM_NUM), min=1, max=24), d=0)], letters[round(runif(length(COMM_NUM), min=1, max=24), d=0)])
							), by='COMM_TYPE']
	dc		<- merge(dc, tmp, by=c('COMM_TYPE','COMM_NUM'))
	setkey(dc, COMM_NUM_A)
	
	set(dc, dc[, which(COMM_NUM=='89')], 'COMM_NUM_A', 'aop')
	set(dc, dc[, which(COMM_NUM=='45')], 'COMM_NUM_A', 'arb')
	set(dc, dc[, which(COMM_NUM=='755')], 'COMM_NUM_A', 'ate')
	set(dc, dc[, which(COMM_NUM=='84')], 'COMM_NUM_A', 'awb')
	stopifnot( !nrow(unique(dc, by='COMM_NUM')[which(duplicated(unique(dc, by='COMM_NUM'),by='COMM_NUM_A')),]) )
	write.csv(dc, row.names=FALSE, file=outfile)	
	#
	#	anonymize locations
	#
	require(ggmap)
	zm					<- get_googlemap(center="rakai district uganda", zoom=10, maptype="hybrid")
	infile.comms		<- '~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/PANGEA_Rakai_community_anonymized_IDs.csv'
	infile.loc.a		<- '~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/community_spatialLoc.csv'
	infile.loc			<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/community_geography.rda"
	outfile.base		<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_"
	#	original locations	
	dc		<- as.data.table(read.csv(infile.comms))
	zc		<- as.data.table(read.csv(infile.loc.a))
	load(infile.loc)
	comgps	<- as.data.table(comgps)
	setnames(comgps, 'COMM_NUM', 'COMM_NUM_RAW')
	set(comgps, NULL, 'COMM_NUM_RAW', comgps[, as.integer(as.character(COMM_NUM_RAW))])
	dc		<- merge(dc, comgps, by='COMM_NUM_RAW',all.x=1)
	set.seed(42L)
	tmp	<- unique(subset(dc, select=c(COMM_NUM, longitude, latitude)), by='COMM_NUM')
	tmp[, LONG_A:= rnorm(nrow(tmp), mean=longitude, sd=0.1)]
	tmp[, LAT_A:= rnorm(nrow(tmp), mean=latitude, sd=0.1)]
	set(tmp, tmp[, which(COMM_NUM=='770')], 'LONG_A', 31.8)
	set(tmp, tmp[, which(COMM_NUM=='770')], 'LAT_A', -0.59)
	set(tmp, tmp[, which(COMM_NUM=='771')], 'LONG_A', 31.75)
	set(tmp, tmp[, which(COMM_NUM=='771')], 'LAT_A', -0.76)	
	set(tmp, tmp[, which(COMM_NUM=='38')], 'LONG_A', 31.8)
	set(tmp, tmp[, which(COMM_NUM=='38')], 'LAT_A', -1.05)
	set(tmp, tmp[, which(COMM_NUM=='772')], 'LONG_A', 31.5)
	set(tmp, tmp[, which(COMM_NUM=='772')], 'LAT_A', -1)	
	set(tmp, tmp[, which(COMM_NUM=='370')], 'LONG_A', 31.59)
	set(tmp, tmp[, which(COMM_NUM=='370')], 'LAT_A', -1.04)	
	set(tmp, tmp[, which(COMM_NUM=='23')], 'LONG_A', 31.65)
	set(tmp, tmp[, which(COMM_NUM=='23')], 'LAT_A', -0.9)
	set(tmp, tmp[, which(COMM_NUM=='16m')], 'LONG_A', 31.4)
	set(tmp, tmp[, which(COMM_NUM=='16m')], 'LAT_A', -1.1)
	set(tmp, tmp[, which(COMM_NUM=='51m')], 'LONG_A', 31.66)
	set(tmp, tmp[, which(COMM_NUM=='51m')], 'LAT_A', -0.43)
	set(tmp, tmp[, which(COMM_NUM=='40')], 'LONG_A', 31.51)
	set(tmp, tmp[, which(COMM_NUM=='40')], 'LAT_A', -0.77)
	set(tmp, tmp[, which(COMM_NUM=='56')], 'LONG_A', 31.48)
	set(tmp, tmp[, which(COMM_NUM=='56')], 'LAT_A', -0.85)
	set(tmp, tmp[, which(COMM_NUM=='108')], 'LONG_A', 31.3)
	set(tmp, tmp[, which(COMM_NUM=='108')], 'LAT_A', -0.8)
	set(tmp, tmp[, which(COMM_NUM=='34')], 'LONG_A', 31.4)
	set(tmp, tmp[, which(COMM_NUM=='34')], 'LAT_A', -0.55)
	set(tmp, tmp[, which(COMM_NUM=='57')], 'LONG_A', 31.45)
	set(tmp, tmp[, which(COMM_NUM=='57')], 'LAT_A', -0.39)
	set(tmp, tmp[, which(COMM_NUM=='58')], 'LONG_A', 31.52)
	set(tmp, tmp[, which(COMM_NUM=='58')], 'LAT_A', -0.74)
	set(tmp, tmp[, which(COMM_NUM=='74')], 'LONG_A', 31.68)
	set(tmp, tmp[, which(COMM_NUM=='74')], 'LAT_A', -0.40)
	set(tmp, tmp[, which(COMM_NUM=='754')], 'LONG_A', 31.5)
	set(tmp, tmp[, which(COMM_NUM=='754')], 'LAT_A', -0.8)
	set(tmp, tmp[, which(COMM_NUM=='5')], 'LONG_A', 31.63)
	set(tmp, tmp[, which(COMM_NUM=='5')], 'LAT_A', -0.5)
	set(tmp, tmp[, which(COMM_NUM=='6')], 'LONG_A', 31.41)
	set(tmp, tmp[, which(COMM_NUM=='6')], 'LAT_A', -0.48)
	set(tmp, tmp[, which(COMM_NUM=='62')], 'LONG_A', 31.51)
	set(tmp, tmp[, which(COMM_NUM=='62')], 'LAT_A', -0.46)
	set(tmp, tmp[, which(COMM_NUM=='89')], 'LONG_A', 31.55)
	set(tmp, tmp[, which(COMM_NUM=='89')], 'LAT_A', -0.62)
	set(tmp, tmp[, which(COMM_NUM=='391')], 'LONG_A', 31.53)
	set(tmp, tmp[, which(COMM_NUM=='391')], 'LAT_A', -0.65)	
	dc	<- merge(dc, subset(tmp, select=c(COMM_NUM, LONG_A, LAT_A)), by='COMM_NUM')
	write.csv(dc, row.names=FALSE, file=infile.comms)
	ggmap(zm) +
			geom_point(data=unique(dc, by='COMM_NUM'), aes(x=LONG_A, y=LAT_A, pch=COMM_TYPE, colour=COMM_TYPE), size=8, alpha=0.8) +
			geom_text(data=unique(dc, by='COMM_NUM'), aes(x=LONG_A, y=LAT_A, label=COMM_NUM_A), nudge_x=0, nudge_y=0, size=3, colour='black')
	ggsave(file=paste0(outfile.base,'_hubs_comm_locations_anonymized.pdf'), w=10, h=10)
	
}

RakaiFull.divergent.clades.classify.170811<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
		
	infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/divsubgraphs_170704_w250_s20_p35_stagetwo_rerun23_min30_cnt_info.rda'
	load(infile)
	dd[, F:=NULL]
	set(dd, NULL, 'INDIVIDUAL', dd[, as.character(INDIVIDUAL)])
	set(dd, NULL, 'SPLIT', dd[, as.character(SPLIT)])
	df		<- as.data.frame(dd)
	
}
	




RakaiFull.analyze.couples.todi.170421<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_couples_170428_withmetadata.rda"		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_couples_170428_"	
	load(infile)	
		
	
	#	extra couple transmission among stable couples
	#	COUPLE==1 means both male or female form stable couple
	#	All RIDs in this data.table are part of a stable couple
			
	# extra couple transmission involving stable couples
	# count: within couple transmission vs transmission to someone else vs transmission from someone else	
	#	- fisherfolk vs agrarian
	z	<- subset(rca, TRM_TYPE1!='neither likely pair nor unlinked')
	z[, DUMMY:=NA_character_]
	set(z, z[, which(TRM_TYPE1=='within couple' & TYPE%in%c('mf','fm'))], 'DUMMY', 'within couple, direction resolved')
	set(z, z[, which(TRM_TYPE1=='extra couple' & TYPE%in%c('mf','fm') & TRM_TYPE2=='couple sink')], 'DUMMY', 'extra couple, transmission into couple, direction resolved')
	set(z, z[, which(TRM_TYPE1=='extra couple' & TYPE%in%c('mf','fm') & TRM_TYPE2=='couple source')], 'DUMMY', 'extra couple, transmission from couple, direction resolved')
	set(z, z[, which(TRM_TYPE1=='extra couple' & TYPE%in%c('other'))], 'DUMMY', 'extra couple, transmission into couple, unsampled transmitter')
	subset(z, !is.na(DUMMY) & FEMALE_COMM_TYPE==MALE_COMM_TYPE)[, table(DUMMY, FEMALE_COMM_TYPE)]
	#																 COMM_TYPE
	#DUMMY                                                           agrarian fisherfolk trading
	#extra couple, transmission from couple, direction resolved           1         33       0
	#extra couple, transmission into couple, direction resolved           1         16       0
	#extra couple, transmission into couple, unsampled transmitter       10         67       2
	#within couple, direction resolved                                   17         42       1
	subset(z, !is.na(DUMMY) & FEMALE_COMM_TYPE==MALE_COMM_TYPE & !grepl('from couple',DUMMY))[, round(prop.table(table(DUMMY, FEMALE_COMM_TYPE), 2), d=3)]
	#																 COMM_TYPE
	#DUMMY                                                           agrarian fisherfolk trading
	#extra couple, transmission into couple, direction resolved       0.036      0.128   0.000
	#extra couple, transmission into couple, unsampled transmitter    0.357      0.536   0.667
	#within couple, direction resolved                                0.607      0.336   0.333
	#
	#--> 	we found much less transmission within stable couples in fisherfolk communities
	#		note that "unsampled transmitter" may mean that two infected individuals then formed a couple

	
	# exclude couples with ++ individuals
	tmp	<- subset(rp, COUP_SC!='seropos')
	z	<- subset(rca, TRM_TYPE1!='neither likely pair nor unlinked' & (MALE_RID%in%tmp$MALE_RID | FEMALE_RID%in%tmp$FEMALE_RID) )
	z[, DUMMY:=NA_character_]
	set(z, z[, which(TRM_TYPE1=='within couple' & TYPE%in%c('mf','fm'))], 'DUMMY', 'within couple, direction resolved')
	set(z, z[, which(TRM_TYPE1=='extra couple' & TYPE%in%c('mf','fm') & TRM_TYPE2=='couple sink')], 'DUMMY', 'extra couple, transmission into couple, direction resolved')
	set(z, z[, which(TRM_TYPE1=='extra couple' & TYPE%in%c('mf','fm') & TRM_TYPE2=='couple source')], 'DUMMY', 'extra couple, transmission from couple, direction resolved')
	set(z, z[, which(TRM_TYPE1=='extra couple' & TYPE%in%c('other'))], 'DUMMY', 'extra couple, transmission into couple, unsampled transmitter')
	subset(z, !is.na(DUMMY) & FEMALE_COMM_TYPE==MALE_COMM_TYPE)[, table(DUMMY, FEMALE_COMM_TYPE)]
	#																 FEMALE_COMM_TYPE
	#DUMMY                                                           agrarian fisherfolk trading
	#extra couple, transmission from couple, direction resolved           1          7       0
	#extra couple, transmission into couple, direction resolved           0          5       0
	#extra couple, transmission into couple, unsampled transmitter        2         13       1
	#within couple, direction resolved                                    4         11       0
	subset(z, !is.na(DUMMY) & FEMALE_COMM_TYPE==MALE_COMM_TYPE & !grepl('from couple',DUMMY))[, round(prop.table(table(DUMMY, FEMALE_COMM_TYPE), 2), d=3)]
	#																 FEMALE_COMM_TYPE
	#DUMMY                                                           agrarian fisherfolk trading
	#extra couple, transmission into couple, direction resolved       0.000      0.172   0.000
	#extra couple, transmission into couple, unsampled transmitter    0.333      0.448   1.000
	#within couple, direction resolved                                0.667      0.379   0.000
	#
	#--> 	same overall conclusion


	

	subset(rca, COUPLE==1 & FEMALE_COMM_TYPE=='trading')[, table(TRM_TYPE1, FEMALE_COMM_TYPE)]
	
	#
	#	make basic epi plot: when positive, when negative, when sequenced
	#	
	t.posneg	<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
	setnames(t.posneg, c('BIRTHDATE','EST_DATEDIED'), c('DOB','DOD'))
	t.seq		<- unique(subset(rs, !is.na(PID), select=c(RID, PID, DATE)))
	setnames(t.seq, 'DATE', 'SEQ_DATE')
	t.cd4		<- unique(subset(rd, !is.na(RECENTCD4DATE) & !is.na(RECENTCD4), select=c(RID, RECENTCD4DATE, RECENTCD4)))
	set(t.cd4, NULL, 'RECENTCD4', t.cd4[, cut(RECENTCD4, breaks=c(-1,250,350,500,800,Inf), labels=c('<200','200-349','350-499','500-799','800+'))])
	t.vl		<- unique(subset(rd, !is.na(RECENTVLDATE) & !is.na(RECENTVL), select=c(RID, RECENTVLDATE, RECENTVL)))
	set(t.vl, NULL, 'RECENTVL', t.vl[, cut(RECENTVL, breaks=c(-1,200,1e3,1e5,Inf), labels=c('<200','200-1,000','1,000-100,000','100,000+'))])
	#	plot 3 timelines per page
	setkey(rtpd, TYPE, POSTERIOR_SCORE)
	rtpd[, DUMMY:= ceiling(seq_len(nrow(rtpd))/50)]	
	p			<- lapply(rtpd[, unique(DUMMY)], function(i){
				cat('\n', i)
				df		<- subset(rtpd, DUMMY==i)
				p		<- RakaiFull.plot.epitimeline(df, copy(t.posneg), copy(t.cd4), copy(t.vl), copy(t.seq), age.adult=14)				
			})
	pdf.n	<- 2
	pi	<- data.table(IDX=seq_along(p))
	pi[, PLOT:= ceiling(IDX/pdf.n)]
	pi[, PLOT_IDX:= (IDX-1)%%pdf.n+1]	
	pdf(file=paste0(outfile.base,'epilines.pdf'), w=20, h=12)	
	for(plot in pi[, unique(PLOT)])
	{
		idx			<- subset(pi, PLOT==plot)[, IDX]
		plot.idx	<- subset(pi, PLOT==plot)[, PLOT_IDX]
		grid.newpage()
		pushViewport(viewport(layout=grid.layout(1, pdf.n)))
		for(i in seq_along(idx))
			print(p[[idx[i]]], vp = viewport(layout.pos.row=1, layout.pos.col=plot.idx[i]))
	}
	dev.off()
	rtpd[, DUMMY:=NULL]
	#pdf(file=paste0(outfile.base,'epilines.pdf'), w=8, h=12)
	#print(p)
	#dev.off()
	
	#
	#	select sero-discordant pairs
	#
	sdc		<- copy(rtpd)
	tmp		<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	sdc		<- merge(sdc, tmp, by='MALE_RID')	
	setnames(tmp, colnames(tmp), gsub('MALE','FEMALE',colnames(tmp)))
	sdc		<- merge(sdc, tmp, by='FEMALE_RID')			
	sdc[, SDC_TYPE:=NA_character_]
	set(sdc, sdc[, which(TYPE=='fm' & FEMALE_LASTNEGDATE>MALE_FIRSTPOSDATE)], 'SDC_TYPE', 'incorrect')
	set(sdc, sdc[, which(TYPE=='fm' & FEMALE_FIRSTPOSDATE<=MALE_LASTNEGDATE)], 'SDC_TYPE', 'correct')	
	set(sdc, sdc[, which(TYPE=='mf' & MALE_LASTNEGDATE>FEMALE_FIRSTPOSDATE)], 'SDC_TYPE', 'incorrect')
	set(sdc, sdc[, which(TYPE=='mf' & MALE_FIRSTPOSDATE<=FEMALE_LASTNEGDATE)], 'SDC_TYPE', 'correct')
	sdc		<- subset(sdc, !is.na(SDC_TYPE))
	sdc[, table(TYPE, SDC_TYPE)]
	#		 correct 	incorrect
	#  fm      20         2
	#  mf      18         1
	#--> pretty good
	t.posneg	<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
	setnames(t.posneg, c('BIRTHDATE','EST_DATEDIED'), c('DOB','DOD'))
	t.seq		<- unique(subset(rs, !is.na(PID), select=c(RID, PID, DATE)))
	setnames(t.seq, 'DATE', 'SEQ_DATE')
	t.cd4		<- unique(subset(rd, !is.na(RECENTCD4DATE) & !is.na(RECENTCD4), select=c(RID, RECENTCD4DATE, RECENTCD4)))
	set(t.cd4, NULL, 'RECENTCD4', t.cd4[, cut(RECENTCD4, breaks=c(-1,250,350,500,800,Inf), labels=c('<200','200-349','350-499','500-799','800+'))])
	t.vl		<- unique(subset(rd, !is.na(RECENTVLDATE) & !is.na(RECENTVL), select=c(RID, RECENTVLDATE, RECENTVL)))
	set(t.vl, NULL, 'RECENTVL', t.vl[, cut(RECENTVL, breaks=c(-1,200,1e3,1e5,Inf), labels=c('<200','200-1,000','1,000-100,000','100,000+'))])
	df			<- subset(sdc, SDC_TYPE=='incorrect')
	p			<- RakaiFull.plot.epitimeline(df, copy(t.posneg), copy(t.cd4), copy(t.vl), copy(t.seq), age.adult=14)
	pdf(file=paste0(outfile.base,'epilines_pairs_with_inconsistent_serohistory.pdf'), w=8, h=5)
	print(p)
	dev.off()
	df			<- subset(sdc, SDC_TYPE=='correct')
	p			<- RakaiFull.plot.epitimeline(df, copy(t.posneg), copy(t.cd4), copy(t.vl), copy(t.seq), age.adult=14)
	pdf(file=paste0(outfile.base,'epilines_pairs_with_consistent_serohistory.pdf'), w=8, h=20)
	print(p)
	dev.off()
	
	
	#
	#	multiple donors to same recipient
	#
	tmp	<- subset(rtpd, TYPE=='mf')[, list(N_TR=length(MALE_RID)), by='FEMALE_RID']
	tmp[, table(N_TR>1)]
	#FALSE  TRUE 
	# 179    83 
	#--> not so good
	tmp	<- subset(rtpd, TYPE=='fm')[, list(N_TR=length(FEMALE_RID)), by='MALE_RID']
	tmp[, table(N_TR>1)]
	#FALSE  TRUE 
	#118    33
	#--> not so good
	#are these part of chains ???
	
	
	#	plot 3 timelines per page
	setkey(rtpd, TYPE, POSTERIOR_SCORE)
	rtpd[, DUMMY:= ceiling(seq_len(nrow(rtpd))/50)]	
	p			<- lapply(rtpd[, unique(DUMMY)], function(i){
				cat('\n', i)
				df		<- subset(rtpd, DUMMY==i)
				p		<- RakaiFull.plot.epitimeline(df, copy(t.posneg), copy(t.cd4), copy(t.vl), copy(t.seq), age.adult=14)				
			})
	
	
	
	
	#load( "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )	
	#load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")	
	load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/circumcision/RCCS_SeqInfo_160816.rda')
	rs		<- subset(rs, !is.na(VISIT))		
	#	get epi info
	tmp		<- RakaiCirc.epi.get.info.170208()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#	add sequence dates to rd
	tmp		<- unique(subset(rs, !is.na(PID), select=c(PID, DATE)),by='PID')
	setnames(tmp, 'DATE','SEQDATE')
	rd		<- merge(rd, tmp, by='PID',all.x=1)
	#	focus on those with PANGEA seqs
	rd		<- subset(rd, !is.na(PID))
	#	focus on clinical data and location data closest to time of diagnosis
	tmp		<- rd[, list(VISIT= VISIT[which.min(abs(DATE-FIRSTPOSDATE))]), by='RID']
	ri		<- subset(merge(unique(rd, by=c('RID','VISIT')), tmp, by=c('RID','VISIT')), select=c(RID, VISIT, DATE, BIRTHDATE, RELIGION, REGION, COMM_NUM, HH_NUM, SEX, LASTNEGDATE, FIRSTPOSDATE, RECENTCD4, RECENTCD4DATE, RECENTVL, RECENTVLDATE, ARVSTARTDATE))
	#	add all PANGEA sequences
	ri		<- merge(ri, unique(subset(rd, select=c(RID, PID, SEQDATE))), by='RID')
	#	focus on behaviour data closest to time of diagnosis
	tmp		<- unique(subset(ri, select=c(RID, VISIT)))
	setnames(tmp, 'VISIT','VISIT_DIAG')
	tmp		<- merge(rh, tmp, by=c('RID'))
	tmp		<- merge(tmp[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_DIAG))]), by='RID'], tmp, by=c('RID','VISIT'))
	setnames(tmp, c('VISIT','VISIT_DIAG'),c('VISIT_H','VISIT'))
	set(tmp, NULL, c('SEX','COMM_NUM'), NULL)
	ri		<- merge(ri, tmp, by=c('RID','VISIT'))
	setnames(ri, 'DATE', 'VISIT_DATE')
	#	merge with rd
	tmp		<- copy(ri)            
	setnames(tmp, colnames(tmp), paste0('ID1_',colnames(tmp)))
	setnames(tmp, 'ID1_RID', 'ID1')
	
	dwin	<- merge(dwin, tmp, by='MALE_PID')
	setnames(tmp, colnames(tmp), gsub('MALE','FEMALE',colnames(tmp)))
	dwin	<- merge(dwin, tmp, by='FEMALE_PID')						
	#	reduce likely pairs to male-1 and female-2
	dwin	<- subset(dwin, MALE_SEX=='M' & FEMALE_SEX=='F')
	set(dwin, NULL, 'TYPE', dwin[, gsub('12','mf',TYPE)])						
	set(dwin, NULL, 'TYPE', dwin[, gsub('21','fm',TYPE)])
	set(dwin, NULL, 'TYPE', dwin[, as.character(TYPE)])
	setnames(dwin, colnames(dwin), gsub('ID1','MALE_SANGER_ID',colnames(dwin)))
	setnames(dwin, colnames(dwin), gsub('ID2','FEMALE_SANGER_ID',colnames(dwin)))
	set(dwin, NULL, 'RUN', RUN)
	
	
	
	
	subset(rtp.todi, GROUP=='TYPE_PAIR_TODI')
	
	'TYPE_DIRSCORE_TODI3'
	
	
	run		<- 'RCCS_170410_w250_trB_blNormedOnFly_dirlklprs_'
	dir		<- '/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410'	
	#	load denominator
	tmp		<- RakaiCirc.epi.get.info.170208()
	ra		<- tmp$ra		
	# load couples "rp"
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	rc		<- copy(rp)
	# load pty.run
	load( "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load rd, rh, rs, rp, rpw, rplkl, ptc
	load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/RCCS_170410_w250_trmp_allpairs_posteriors_cmptoprv.rda')
	#
	# select final run
	#
	tmp		<- "RCCS_170410_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl3.5_d8"
	tmp		<- "RCCS_170410_w250_d50_st20_trB_blInScriptNormed_mr20_mt1_cl3.5_d8"
	rpw		<- subset(rpw, RUN%in%tmp )
	rplkl	<- subset(rplkl, RUN%in%tmp )	
	#	add info on pair types to rplkl
	rp		<- copy(rpw)
	set(rp, NULL, c('DIR','FILE','RUN','W_FROM','W_TO','TYPE_RAW','TYPE','GROUP','PATRISTIC_DISTANCE','ADJACENT','CONTIGUOUS','PATHS_12','PATHS_21','MALE_SANGER_ID_L','MALE_SANGER_ID_R','FEMALE_SANGER_ID_L','FEMALE_SANGER_ID_R','CHUNK','CHUNK_L','CHUNK_N','ID_R_MIN','ID_R_MAX'), NULL)
	rp		<- unique(rp)
	#	make COUPID
	rp[, COUPID:= paste0(MALE_RID,':',FEMALE_RID)]	
	#	add PAIR_TYPE
	tmp		<- unique(subset(rc, select=c(COUPID, MALE_HH_NUM, FEMALE_HH_NUM, COUP_SC, PAIR_TYPE)))	
	setnames(tmp, 'COUP_SC', 'COUP_TYPE')
	set(tmp, NULL, c('MALE_HH_NUM','FEMALE_HH_NUM'), NULL)
	rp		<- merge(rp, tmp, by=c('COUPID'),all.x=1)	
	set(rp, rp[, which(!MALE_RID%in%rc[, MALE_RID] & !FEMALE_RID%in%rc[, FEMALE_RID])], 'PAIR_TYPE', 'm and f not in couple')
	set(rp, rp[, which(is.na(PAIR_TYPE))], 'PAIR_TYPE', 'f or m not in couple')	
	tmp		<- subset(rp, select=c(FEMALE_SANGER_ID, MALE_SANGER_ID, MALE_RID, FEMALE_RID, COUPID, PTY_RUN, COUP_TYPE, PAIR_TYPE))
	set(rplkl, NULL, c('MALE_RID','FEMALE_RID','COUPID','COUP_TYPE','PAIR_TYPE'), NULL)
	rplkl	<- merge(tmp, rplkl, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN'))
	set(rplkl, NULL, 'FEMALE_SANGER_ID', rplkl[, as.character(FEMALE_SANGER_ID)])
	set(rplkl, NULL, 'MALE_SANGER_ID', rplkl[, as.character(MALE_SANGER_ID)])
	rplkl	<- unique(rplkl)
	set(rpw, NULL, 'FEMALE_SANGER_ID', rpw[, as.character(FEMALE_SANGER_ID)])
	set(rpw, NULL, 'MALE_SANGER_ID', rpw[, as.character(MALE_SANGER_ID)])		
	#
	#	basic info on selection
	#
	if(0)
	{
		tmp		<- unique(subset(rc, !is.na(MALE_TAXA) & !is.na(FEMALE_TAXA) & PAIR_TYPE=='stable cohabiting'), by='COUPID')
		z		<- unique(subset(ri, select=c(COMM_NUM, COMM_TYPE)))
		setnames(z, c('COMM_NUM','COMM_TYPE'),c('MALE_COMM_NUM','MALE_COMM_TYPE'))
		tmp		<- merge(tmp, z, by= 'MALE_COMM_NUM')
		setnames(z, c('MALE_COMM_NUM','MALE_COMM_TYPE'), c('FEMALE_COMM_NUM','FEMALE_COMM_TYPE'))
		tmp		<- merge(tmp, z, by= 'FEMALE_COMM_NUM')
		
		tmp		<- unique(rp, by='COUPID')
		tmp[, c(length(unique(MALE_RID)),length(unique(FEMALE_RID)))]
		nrow(unique(subset(tmp, PAIR_TYPE=='stable cohabiting'), by=c('FEMALE_RID','MALE_RID')))
		subset(tmp, PAIR_TYPE=='stable cohabiting')[, c(length(unique(MALE_RID)),length(unique(FEMALE_RID)))]
		nrow(subset(tmp, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_NUM!=MALE_COMM_NUM))
		subset(tmp, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_NUM==MALE_COMM_NUM)[, table(FEMALE_COMM_TYPE)]
		subset(tmp, !is.na(COUP_TYPE))[, c(length(unique(MALE_RID)),length(unique(FEMALE_RID)),length(unique(COUPID)))]		
		#	2 males with new partners:  G110085 J189465
		#	4 females with new partners:  C106054 H104287 B105985 E111070
	}
	
	#	select likely transmitters (unsampled intermediate not necessarily excluded) 
	#	find pairs for whom 'likely pair' is most likely state
	#	(does not depend on prior or confidence cut)
	rex		<- subset(rplkl, GROUP=='TYPE_PAIR_TODI')[, list(TYPE_MLE=TYPE[which.max(KEFF)], KEFF=max(KEFF)), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','COUPID')]
	rex		<- subset(rex, TYPE_MLE=='distant')
	#	select one sequence pairing per couple: that with highest evidence
	rex		<- rex[, {
				z<- which.max(KEFF)
				list(MALE_SANGER_ID=MALE_SANGER_ID[z], FEMALE_SANGER_ID=FEMALE_SANGER_ID[z], PTY_RUN=PTY_RUN[z])
			}, by='COUPID']
	set(rex, NULL, 'COUPID', NULL)
	#	calculate confidence score and select
	confidence.cut	<- 0.5
	rex		<- merge(rex, subset(rplkl, GROUP=='TYPE_PAIRSCORE_TODI' & TYPE=='distant'), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all.x=1)	
	rex[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
	rex		<- subset(rex, POSTERIOR_SCORE>confidence.cut)
	
	
	#	select likely transmitters (unsampled intermediate not necessarily excluded) 
	#	find pairs for whom 'likely pair' is most likely state
	#	(does not depend on prior or confidence cut)
	rtp		<- subset(rplkl, GROUP=='TYPE_PAIR_TODI')[, list(TYPE_MLE=TYPE[which.max(KEFF)], KEFF=max(KEFF)), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','COUPID')]
	rtp		<- subset(rtp, TYPE_MLE=='likely pair')
	#	select one sequence pairing per couple: that with highest evidence
	rtp		<- rtp[, {
				z<- which.max(KEFF)
				list(MALE_SANGER_ID=MALE_SANGER_ID[z], FEMALE_SANGER_ID=FEMALE_SANGER_ID[z], PTY_RUN=PTY_RUN[z])
			}, by='COUPID']
	set(rtp, NULL, 'COUPID', NULL)
	#	calculate confidence score and select
	confidence.cut	<- 0.5
	rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIRSCORE_TODI' & TYPE=='likely pair'), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all.x=1)	
	rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
	rtp		<- subset(rtp, POSTERIOR_SCORE>confidence.cut)
	
	#	resolve direction
	#	find likely pairs for whom 'mf' or 'fm' is most likely state
	#	(does not depend on prior or confidence cut)
	rtpd	<- subset(rtp, select=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	rtpd	<- merge(rtpd, subset(rplkl, GROUP=='TYPE_DIR_TODI3'), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rtpd	<- rtpd[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	rtpd	<- subset(rtpd, TYPE_MLE!='ambiguous')	
	#	calculate confidence score and select
	confidence.cut	<- 0.5
	setnames(rtpd, 'TYPE_MLE','TYPE')
	rtpd	<- merge(rtpd, subset(rplkl, GROUP=='TYPE_DIRSCORE_TODI3'), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','TYPE'))
	rtpd[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
	rtpd	<- subset(rtpd, POSTERIOR_SCORE>confidence.cut)
	rmf		<- subset(rtpd, TYPE=='mf')
	rfm		<- subset(rtpd, TYPE=='fm')
	
	#subset(rtp, PTY_RUN==28 & MALE_SANGER_ID=='15714_1_84' & FEMALE_SANGER_ID=='15862_1_86')
	#subset(rtpd, PTY_RUN==67 & MALE_SANGER_ID=='15965_1_24' & FEMALE_SANGER_ID=='15977_1_52')
	#subset(rplkl, PTY_RUN==28 & MALE_SANGER_ID=='15714_1_84' & FEMALE_SANGER_ID=='15862_1_86' & GROUP=='TYPE_DIRSCORE_TODI3')
	#subset(rplkl, PTY_RUN==28 & MALE_SANGER_ID=='15714_1_84' & FEMALE_SANGER_ID=='15862_1_86' & GROUP=='TYPE_DIR_TODI7x3')
	
	#	info		
	cat('\ncouples with phyloscanner assessment, n=',				nrow(unique(rplkl, by='COUPID')))	
	cat('\ncouples not implicated in transmission, n=',				nrow(unique(rex, by='COUPID')))
	unique(rex,by='COUPID')[, table(PAIR_TYPE)]
	cat('\ncouples that are likely pairs, n=',						nrow(unique(rtp, by='COUPID')))
	cat('\ncouples that are likely pairs with evidence M->F, n=',	nrow(unique(rmf, by='COUPID')))
	cat('\ncouples that are likely pairswith evidence F->M, n=',	nrow(unique(rfm, by='COUPID')))
	#	pairings assessed, n= 1741
	#	couples not implicated in transmission, n= 1402
	#		    f or m not in couple m and f not in couple not always cohabiting     stable cohabiting 
	#             1311                    16                     8                    67  
	#	couples that are likely pairs, n= 209
	#	likely direction resolved, n= 127
	#	   not always cohabiting 			   not registered as couple        stable cohabiting 
	#                  2                       41                              85
	#	couples that are likely pairs with evidence M->F, n= 84
	#	couples that are likely pairswith evidence F->M, n= 43
	
	#	define two helper data.table
	rmf		<- merge(unique(subset(rmf, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID))), rp, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rmf[, PHSC_DIR:='m->f']
	rfm		<- merge(unique(subset(rfm, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID))), rp, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rfm[, PHSC_DIR:='f->m']
	rtr		<- rbind(rmf, rfm)	
	rtr[, AGEDIFF:= FEMALE_BIRTHDATE-MALE_BIRTHDATE]
	rtr[, AVGAGE:= (MALE_BIRTHDATE+FEMALE_BIRTHDATE)/2]	
	rtr2	<- copy(rmf)
	setnames(rtr2,colnames(rtr2),gsub('FEMALE','REC',colnames(rtr2)))
	setnames(rtr2,colnames(rtr2),gsub('MALE','TR',colnames(rtr2)))
	tmp		<- copy(rfm)
	setnames(tmp,colnames(tmp),gsub('FEMALE','TR',colnames(tmp)))
	setnames(tmp,colnames(tmp),gsub('MALE','REC',colnames(tmp)))
	rtr2	<- rbind(rtr2,tmp)
	
	if(0)	#for inscript
	{
		rps			<- subset(rtr, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID))
		outfile.base<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/lklpr_TODI_'
		write.csv(rps, file=paste0(outfile.base,'_summary_versionmaxscore_inscript.csv'))
		#
		#	get difference from manual check (versionstrict2) to more relaxed version
		#		
		tmp2		<- subset(as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/lklpr_TODI__summary_versionmaxscore.csv')), select=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
		tmp2[, VERSION:='on the fly']
		tmp			<- copy(rps)	
		tmp[, VERSION_NEW:='in script']
		tmp			<- merge(tmp, tmp2, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all=1)
		tmp[, table(is.na(VERSION_NEW), is.na(VERSION))]
		#		FALSE TRUE
		#FALSE   109    4
		#TRUE     18    0
		rps			<- subset(tmp, is.na(VERSION) | is.na(VERSION_NEW))	
		write.csv(rps, file=paste0(outfile.base,'_summary_diffaftermaxscore_inscript.csv'))
		set(rps, NULL, c('VERSION','VERSION_NEW'), NULL)
		group		<- 'TYPE_DIR_TODI7x3'
		#group		<- 'TYPE_PAIR_TODI'
		run			<- "RCCS_170410_w250_d50_st20_trB_blInScriptNormed_mr20_mt1_cl3.5_d8"	
		plot.select	<- unique(subset(merge(rplkl, rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP==group), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
		rpw2		<- subset(rpw, RUN==run & GROUP==group)
		rplkl2		<- subset(rplkl, RUN==run & GROUP==group)	
		plot.file	<- paste0(outfile.base,'_windows_summary_',group,'_diffaftermaxscore_inscript.pdf')	
		phsc.plot.windowsummaries.for.pairs(plot.select, rpw2, rplkl2, plot.file, cols=NULL, group=group)
	}	
	if(0)	#for onthefly
	{
		rps			<- subset(rtr, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID))
		outfile.base<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/lklpr_TODI_'
		write.csv(rps, file=paste0(outfile.base,'_summary_versionmaxscore.csv'))
		group		<- 'TYPE_DIR_TODI7x3'
		#group		<- 'TYPE_PAIR_TODI'
		run			<- "RCCS_170410_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl3.5_d8"	
		plot.select	<- unique(subset(merge(rplkl, rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP==group), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
		rpw2		<- subset(rpw, RUN==run & GROUP==group)
		rplkl2		<- subset(rplkl, RUN==run & GROUP==group)	
		plot.file	<- paste0(outfile.base,'_windows_summary_',group,'_versionmaxscore.pdf')	
		phsc.plot.windowsummaries.for.pairs(plot.select, rpw2, rplkl2, plot.file, cols=NULL, group=group)
		#
		#	get difference from manual check (versionstrict2) to more relaxed version
		#		
		tmp2		<- subset(as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/lklpr_TODI__summary_versionstrict2.csv')), select=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
		tmp2[, VERSION:='checked manually']
		tmp			<- copy(rps)	
		tmp[, VERSION_NEW:='first max then score']
		tmp			<- merge(tmp, tmp2, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all=1)
		tmp[, table(is.na(VERSION_NEW), is.na(VERSION))]
		#		FALSE TRUE
		#   FALSE    59   68
		rps			<- subset(tmp, is.na(VERSION) | is.na(VERSION_NEW))	
		write.csv(rps, file=paste0(outfile.base,'_summary_diffaftermaxscore.csv'))
		set(rps, NULL, c('VERSION','VERSION_NEW'), NULL)
		group		<- 'TYPE_DIR_TODI7x3'
		group		<- 'TYPE_PAIR_TODI'
		run			<- "RCCS_170410_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl3.5_d8"	
		plot.select	<- unique(subset(merge(rplkl, rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP==group), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
		rpw2		<- subset(rpw, RUN==run & GROUP==group)
		rplkl2		<- subset(rplkl, RUN==run & GROUP==group)	
		plot.file	<- paste0(outfile.base,'_windows_summary_',group,'_diffaftermaxscore.pdf')	
		phsc.plot.windowsummaries.for.pairs(plot.select, rpw2, rplkl2, plot.file, cols=NULL, group=group)		
	}
	
	
	
	
	#	plot phylogenies for a few examples
	tmp			<- subset(rps, PTY_RUN%in%c(113) & FEMALE_SANGER_ID=='15958_1_47')
	set(tmp, NULL, 'FEMALE_SANGER_ID', tmp[, as.character(FEMALE_SANGER_ID)])
	set(tmp, NULL, 'MALE_SANGER_ID', tmp[, as.character(MALE_SANGER_ID)])
	run			<- 'RCCS_170410_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl3.5_d8'
	rpw2		<- unique(subset(rpw, RUN==run, select=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM','W_TO','PATHS_12','PATHS_21','PATRISTIC_DISTANCE','CONTIGUOUS','TYPE_RAW')))
	set(rpw2, NULL, 'FEMALE_SANGER_ID', rpw2[, as.character(FEMALE_SANGER_ID)])
	set(rpw2, NULL, 'MALE_SANGER_ID', rpw2[, as.character(MALE_SANGER_ID)])
	set(rpw2, NULL, 'CONTIGUOUS', rpw2[, as.integer(CONTIGUOUS)])
	set(rpw2, NULL, 'PATRISTIC_DISTANCE', rpw2[, round(PATRISTIC_DISTANCE, d=4)])	
	#load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/RCCS_170410_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl3.5_d8_phscout.rda')
	invisible(sapply(seq_len(nrow(tmp)), function(ii)
					{	
						#ii<- 1
						ids			<- c(tmp[ii, MALE_SANGER_ID],tmp[ii, FEMALE_SANGER_ID])
						pty.run		<- tmp[ii, PTY_RUN]
						dfs			<- subset(dtrees, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
						dfs[, MALE_SANGER_ID:=ids[1]]
						dfs[, FEMALE_SANGER_ID:=ids[2]]
						dfs			<- merge(dfs, rpw2, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM','W_TO'), all.x=TRUE)
						dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', pty.run, '\nwindow ', W_FROM,'-', W_TO,'\n',PATHS_12,' ',PATHS_21, ' ',CONTIGUOUS,' ',TYPE_RAW, '\n', PATRISTIC_DISTANCE, sep='')]]								
						plot.file	<- paste0(outfile.base, pty.run,'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
						invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10))					
					}))	
	
	subset(rplkl, PTY_RUN==44 & MALE_SANGER_ID=='15743_1_40' & FEMALE_SANGER_ID=='15115_1_22' & GROUP=='TYPE_PAIRSCORE_TODI' & TYPE=='likely pair')[, pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
	subset(rplkl, PTY_RUN==44 & MALE_SANGER_ID=='15743_1_40' & FEMALE_SANGER_ID=='15115_1_22' & GROUP=='TYPE_DIRSCORE_TODI3' & TYPE=='fm')[, pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
	
	#
	#rtr2[, table(PAIR_TYPE)]
	#
	#	who infects whom matrix
	#
	tmp		<- rtr2[,list(N=length(unique(COUPID))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]
	ggplot(tmp, aes(x=factor(REC_COMM_TYPE),y=factor(TR_COMM_TYPE))) + 
			geom_point(aes(size=N), colour='grey80') +
			geom_text(aes(label=N), nudge_x=0, nudge_y=0, size=3, colour='black') +			
			theme_bw() + 
			scale_size(range = c(5, 50)) +
			labs(x='\nlocation likely recipient',y='location likely transmitter\n') +
			guides(size='none')
	ggsave(file=file.path(dir, paste(run,'-phsc-directionpairs_direction-numbers-commtype.pdf',sep='')), w=4, h=4)
	#
	#	did any transmitter start ART before the recipient was diagnosed?
	subset(rtr2, TR_ARVSTARTDATE<REC_FIRSTPOSDATE)	
	#	F026858:J104288 --> stable couple, rec male, first diagnosed with v high CD4 (2400), about 2 years after female started ART 
	#	C066263:K077878 --> no couple, rec female, first diagnosed with v high CD4, about 5m after male started ART
	
	#
	#	how many transmitters were positive for 6m before the recipient was found positive
	subset(rtr2, (TR_FIRSTPOSDATE+.5)<REC_FIRSTPOSDATE)	
	#	26
	
	
	subset(rtr, MALE_COMM_TYPE==FEMALE_COMM_TYPE & FEMALE_COMM_TYPE!='trading')[, {
				z	<- binconf( length(which(PHSC_DIR=='m->f')), length(PHSC_DIR) )				
				list(K=length(which(PHSC_DIR=='m->f')), N=length(PHSC_DIR), P=z[1], QL=z[2], QU=z[3])
			}, by=c('MALE_COMM_TYPE')]	
	#	   MALE_COMM_TYPE  	K  N  P         QL        QU
	#1:       agrarian 		27 38 0.7105263 0.5524286 0.8299672
	#2:     fisherfolk 		55 85 0.6470588 0.5411250 0.7402751
	
	#
	#	is there a difference in male->female transmission by couple type?
	#	results: no	
	#
	tmp		<- copy(rtr)
	set(tmp, tmp[, which(PAIR_TYPE!='stable cohabiting')], 'PAIR_TYPE', 'no stable pair')
	tmp[, {
				z	<- binconf( length(which(PHSC_DIR=='m->f')), length(PHSC_DIR) )				
				list(K=length(which(PHSC_DIR=='m->f')), N=length(PHSC_DIR), P=z[1], QL=z[2], QU=z[3])
			}, by='PAIR_TYPE']	
	#			PAIR_TYPE  	 K  N         P        QL        QU
	#1:    stable cohabiting 59 86 0.6860465 0.5817960 0.7743870
	#2:    no stable pair 	 25 42 0.5952381 0.4449431 0.7295714	
	chisq_test(factor(PHSC_DIR) ~ factor(PAIR_TYPE), data=tmp, distribution="exact")
	#	chi-squared = 1.0315, p-value = 0.3278
	
	
	#	are transmitters younger in fisherfolk sites?
	#	results: yes
	tmp		<- subset(rtr2, TR_COMM_TYPE!='trading')
	ggplot(tmp, aes(x=TR_COMM_TYPE, y=TR_BIRTHDATE)) + geom_boxplot()
	independence_test(TR_BIRTHDATE~factor(TR_COMM_TYPE), data=tmp, distribution = "exact")
	#	Z = -2.3289, p-value = 0.01934
	ggplot(tmp, aes(x=REC_COMM_TYPE, y=REC_BIRTHDATE)) + geom_boxplot()
	independence_test(REC_BIRTHDATE~factor(REC_COMM_TYPE), data=tmp, distribution = "exact")
	#	Z = -2.2714, p-value = 0.02236
	#	summary(rq(TR_BIRTHDATE~TR_COMM_TYPE, tau=.5, data=tmp, method='fn'), se='nid')
	
	#
	#	is there a difference in age gap between male->female transmission / female->male transmission ?
	#	results: not significant but outside couples, men tend to be infected by much younger women
	#	
	tmp		<- subset(rtr, !grepl('cohabiting',PAIR_TYPE) & !is.na(AGEDIFF) & FEMALE_COMM_TYPE==MALE_COMM_TYPE & MALE_COMM_TYPE=='fisherfolk')
	independence_test(AGEDIFF~factor(PHSC_DIR), data=tmp, distribution = "exact")
	#	Z = 1.5902, p-value = 0.1134
	tmp		<- subset(rtr, !grepl('cohabiting',PAIR_TYPE) & !is.na(AGEDIFF) & FEMALE_COMM_TYPE==MALE_COMM_TYPE & MALE_COMM_TYPE=='agrarian')
	independence_test(AGEDIFF~factor(PHSC_DIR), data=tmp, distribution = "exact")
	#	Z = 1.7439, p-value = 0.1429
	
	
	ggplot(rtr, aes(x=PHSC_DIR, y=AGEDIFF)) + geom_boxplot()	
	tmp		<- subset(rtr, !is.na(MALE_BIRTHDATE) & !is.na(FEMALE_BIRTHDATE), select=c(PHSC_DIR,AGEDIFF))
	set(tmp, NULL, 'PHSC_DIR', tmp[, as.integer(as.character(factor(PHSC_DIR, levels=c('f->m','m->f'), labels=c('0','1'))))])
	summary(gamlss(data=tmp, PHSC_DIR~AGEDIFF, family=LO))
	summary(gamlss(data=tmp, AGEDIFF~PHSC_DIR))
	#				Estimate Std. Error t value Pr(>|t|)    
	#(Intercept)  0.626998   0.071579   8.759 1.98e-13 ***
	#AGEDIFF     -0.002135   0.008422  -0.254      0.8    
	tmp		<- subset(rtr, MALE_COMM_TYPE!='trading' & MALE_COMM_TYPE==FEMALE_COMM_TYPE)
	ggplot(tmp, aes(x=PHSC_DIR, y=AGEDIFF)) + 
			geom_boxplot() + 
			facet_grid(~MALE_COMM_TYPE)	+
			theme_bw() + labs(x='\nestimated direction of transmission', y='age difference male-female\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-directionpairs_direction-agegap-commtype.pdf',sep='')), w=4, h=6)
	
	#	AAA
	subset(rtr2, TR_COMM_TYPE!='trading' & PAIR_TYPE!='m and f not in couple' & REC_RID%in%c(rc$MALE_RID,rc$FEMALE_RID))[, {
				m<- length(which(PAIR_TYPE=='stable cohabiting'))
				n<- length(PAIR_TYPE)
				z<- unname(as.numeric(binconf(m, n)))
				list(P=z[1], PL=z[2], PU=z[3], M=m, N=n, TYPE='stable cohabiting')				
			}, by=c('TR_COMM_TYPE','REC_SEX')]
	
	#
	#	Does the primary occupation differ between transmitters / recipients? 
	#	
	tmp			<- subset(rtr, FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_occupation2.pdf',sep=''))
	tmp			<- subset(rtr, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_occupation2_stablecouples.pdf',sep=''))
	tmp			<- subset(rtr, !grepl('cohabiting',PAIR_TYPE) & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_occupation2_nocouples.pdf',sep=''))
	#	
	tmp2		<- unique(subset(ra, VISIT!=17 & SEX=='F' & !is.na(FIRSTPOSDATE) & COMM_TYPE!='trading', c(RID, OCCUP_OLLI, COMM_TYPE)))
	setnames(tmp2, colnames(tmp2), paste0('FEMALE_',colnames(tmp2)))
	tmp2[, PHSC_DIR:='denominator']
	tmp			<- rbind(tmp, tmp2, fill=TRUE)
	tmp2		<- unique(subset(ra, VISIT!=17 & SEX=='M' & !is.na(FIRSTPOSDATE) & COMM_TYPE!='trading', c(RID, OCCUP_OLLI, COMM_TYPE)))
	setnames(tmp2, colnames(tmp2), paste0('MALE_',colnames(tmp2)))
	tmp2[, PHSC_DIR:='denominator']
	tmp			<- rbind(tmp, tmp2, fill=TRUE)
	setnames(tmp, c('FEMALE_OCCUP_OLLI','MALE_OCCUP_OLLI'), c('FEMALE_FACTOR','MALE_FACTOR'))
	cols		<- tmp[, length(na.omit(unique(c(FEMALE_FACTOR, MALE_FACTOR))))]
	cols		<- colorRampPalette(brewer.pal(min(11,cols), "Set3"))( cols )		
	names(cols)	<- tmp[, na.omit(sort(unique(c(FEMALE_FACTOR, MALE_FACTOR))))]
	Rakai.plot.directed.pairs.discrete1(tmp, cols, outfile, 'occupation at diagnosis', w=10, h=7)
	
	# number female Bar/waitress that are transmitters
	ntf	<- nrow(subset(tmp, PHSC_DIR=='f->m' & FEMALE_FACTOR=='Bar/waitress' & FEMALE_COMM_TYPE=='fisherfolk'))
	# number female Bar/waitress that are recipients
	nrf	<- nrow(subset(tmp, PHSC_DIR=='m->f' & FEMALE_FACTOR=='Bar/waitress' & FEMALE_COMM_TYPE=='fisherfolk'))
	# number female Bar/waitress HIV-infected
	ndf	<- nrow(subset(tmp, PHSC_DIR=='denominator' & FEMALE_FACTOR=='Bar/waitress' & FEMALE_COMM_TYPE=='fisherfolk'))
	# number female that are transmitters
	nt	<- nrow(subset(tmp, PHSC_DIR=='f->m' & FEMALE_COMM_TYPE=='fisherfolk'))
	# number female that are recipients
	nr	<- nrow(subset(tmp, PHSC_DIR=='m->f' & FEMALE_COMM_TYPE=='fisherfolk'))
	# number female HIV-infected
	nd	<- nrow(subset(tmp, PHSC_DIR=='denominator' & FEMALE_COMM_TYPE=='fisherfolk'))
	# odds ratio transmitter / recipient
	# 'a' is exposed cases (exposed=Bar/waitress, case=transmitter)
	# I resample by taking p=ntf/nt as the best estimate of the proportion of female Bar/waitress that are transmitters
	# and then adding uncertainty around p based on p(1-p)/n
	bs	<- 1e4
	a	<- round(nt*rnorm(bs, mean=ntf/nt, sd=sqrt( (ntf/nt) * (1-ntf/nt) / nt )))
	# 'b' is exposed non-cases (exposed=Bar/waitress, non-case=recipients)
	b	<- round(nr*rnorm(bs, mean=nrf/nr, sd=sqrt( (nrf/nr) * (1-nrf/nr) / nr )))
	c	<- nt-a
	d	<- nr-b
	tmp2<- quantile( (a/c) / (b/d), p=c(0.025,0.975))
	a	<- ntf
	b	<- nrf
	c	<- nt-a
	d	<- nr-b
	tmp2<- c( (a/c) / (b/d), tmp2)
	
	bs	<- 1e4
	tmp2<- quantile( rnorm(bs, mean=ntf/nt, sd=sqrt( (ntf/nt) * (1-ntf/nt) / nt )) / rnorm(bs, mean=ndf/nd, sd=sqrt( (ndf/nd) * (1-ndf/nd) / nd )), p=c(0.025,0.975))
	tmp2<- c((ntf/nt) / (ndf/nd),tmp2)
	tmp3<- quantile( rnorm(bs, mean=nrf/nr, sd=sqrt( (nrf/nr) * (1-nrf/nr) / nr )) / rnorm(bs, mean=ndf/nd, sd=sqrt( (ndf/nd) * (1-ndf/nd) / nd )), p=c(0.025,0.975))
	tmp3<- c((nrf/nr) / (ndf/nd),tmp3)
	
	
	
	ntf	<- nrow(subset(tmp, PHSC_DIR=='m->f' & MALE_FACTOR=='Fishing' & MALE_COMM_TYPE=='fisherfolk'))
	nrf	<- nrow(subset(tmp, PHSC_DIR=='f->m' & MALE_FACTOR=='Fishing' & MALE_COMM_TYPE=='fisherfolk'))
	ndf	<- nrow(subset(tmp, PHSC_DIR=='denominator' & MALE_FACTOR=='Fishing' & MALE_COMM_TYPE=='fisherfolk'))
	nt	<- nrow(subset(tmp, PHSC_DIR=='m->f' & MALE_COMM_TYPE=='fisherfolk'))
	nr	<- nrow(subset(tmp, PHSC_DIR=='f->m' & MALE_COMM_TYPE=='fisherfolk'))
	nd	<- nrow(subset(tmp, PHSC_DIR=='denominator' & MALE_COMM_TYPE=='fisherfolk'))
	
	bs	<- 1e4
	tmp	<- quantile( rnorm(bs, mean=ntf/nt, sd=sqrt( (ntf/nt) * (1-ntf/nt) / nt )) / rnorm(bs, mean=ndf/nd, sd=sqrt( (ndf/nd) * (1-ndf/nd) / nd )), p=c(0.025,0.975))
	tmp	<- c((ntf/nt) / (ndf/nd),tmp)
	tmp	<- quantile( rnorm(bs, mean=nrf/nr, sd=sqrt( (nrf/nr) * (1-nrf/nr) / nr )) / rnorm(bs, mean=ndf/nd, sd=sqrt( (ndf/nd) * (1-ndf/nd) / nd )), p=c(0.025,0.975))
	tmp	<- c((nrf/nr) / (ndf/nd),tmp)
	#
	#	Age group
	#	
	tmp			<- subset(rtr, FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_birthdate.pdf',sep=''))
	setnames(tmp, c('FEMALE_BIRTHDATE','MALE_BIRTHDATE'), c('FEMALE_FACTOR','MALE_FACTOR'))	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'birth date', bw=4, dotsize=0.5, w=10, h=7)
	tmp			<- subset(rtr, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_birthdate_stablecouples.pdf',sep=''))
	setnames(tmp, c('FEMALE_BIRTHDATE','MALE_BIRTHDATE'), c('FEMALE_FACTOR','MALE_FACTOR'))	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'birth date', bw=4, dotsize=0.5, w=10, h=7)
	tmp			<- subset(rtr, PAIR_TYPE=='not registered as couple' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_birthdate_nocouples.pdf',sep=''))
	setnames(tmp, c('FEMALE_BIRTHDATE','MALE_BIRTHDATE'), c('FEMALE_FACTOR','MALE_FACTOR'))	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'birth date', bw=4, dotsize=0.5, w=10, h=7)
	#	
	tmp			<- subset(rtr, FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_agediff.pdf',sep=''))
	setnames(tmp, c('AGEDIFF'), c('MALE_FACTOR'))
	tmp[, FEMALE_FACTOR:= MALE_FACTOR]	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'age difference\n(years female younger)', bw=4, dotsize=0.5, w=10, h=7)
	tmp			<- subset(rtr, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_agediff_stablecouples.pdf',sep=''))
	setnames(tmp, c('AGEDIFF'), c('MALE_FACTOR'))
	tmp[, FEMALE_FACTOR:= MALE_FACTOR]	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'age difference\n(years female younger)', bw=4, dotsize=0.5, w=10, h=7)
	tmp			<- subset(rtr, PAIR_TYPE=='not registered as couple' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_agediff_nocouples.pdf',sep=''))
	setnames(tmp, c('AGEDIFF'), c('MALE_FACTOR'))
	tmp[, FEMALE_FACTOR:= MALE_FACTOR]	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'age difference\n(years female younger)', bw=4, dotsize=0.5, w=10, h=7)
	
	#
	#	Marriage Status
	#
	tmp			<- subset(rtr, FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_marriagestatus.pdf',sep=''))
	tmp			<- subset(rtr, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_marriagestatus_stablecouples.pdf',sep=''))
	tmp			<- subset(rtr, !grepl('cohabiting',PAIR_TYPE) & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_marriagestatus_nocouples.pdf',sep=''))
	
	set(tmp, NULL, 'MALE_MARSTAT', tmp[, gsub('Never Married \\+ casual partner','Never Married',gsub('Previously Married \\+ casual partner','Previously Married',MALE_MARSTAT))])
	set(tmp, NULL, 'FEMALE_MARSTAT', tmp[, gsub('Never Married \\+ casual partner','Never Married',gsub('Previously Married \\+ casual partner','Previously Married',FEMALE_MARSTAT))])
	tmp2		<- unique(subset(ra, SEX=='F' & !is.na(FIRSTPOSDATE) & COMM_TYPE!='trading', c(RID, MARSTAT, COMM_TYPE)))
	setnames(tmp2, colnames(tmp2), paste0('FEMALE_',colnames(tmp2)))
	tmp2[, PHSC_DIR:='denominator']
	tmp			<- rbind(tmp, tmp2, fill=TRUE)
	tmp2		<- unique(subset(ra, SEX=='M' & !is.na(FIRSTPOSDATE) & COMM_TYPE!='trading', c(RID, MARSTAT, COMM_TYPE)))
	setnames(tmp2, colnames(tmp2), paste0('MALE_',colnames(tmp2)))
	tmp2[, PHSC_DIR:='denominator']
	tmp			<- rbind(tmp, tmp2, fill=TRUE)	
	setnames(tmp, c('FEMALE_MARSTAT','MALE_MARSTAT'), c('FEMALE_FACTOR','MALE_FACTOR'))
	cols		<- tmp[, length(unique(c(FEMALE_FACTOR, MALE_FACTOR)))]
	cols		<- colorRampPalette(brewer.pal(min(8,cols), "Set2"))( cols )	
	#cols		<- rainbow_hcl(tmp[, length(unique(c(FEMALE_FACTOR, MALE_FACTOR)))], start = 20, end = 340, c=100, l=60)
	names(cols)	<- tmp[, sort(unique(c(FEMALE_FACTOR, MALE_FACTOR)))]
	Rakai.plot.directed.pairs.discrete1(tmp, cols, outfile, 'marital &\nself-reported\nnon-marital\nrelationships,\n', w=10, h=7)
	
	
	
	#
	#	Education
	#
	tmp			<- subset(rtr, FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	setnames(tmp, c('FEMALE_EDUCAT','MALE_EDUCAT'), c('FEMALE_FACTOR','MALE_FACTOR'))
	cols		<- tmp[, length(unique(c(FEMALE_FACTOR, MALE_FACTOR)))]
	cols		<- colorRampPalette(brewer.pal(min(11,cols), "Set1"))( cols )	
	#cols		<- rainbow_hcl(tmp[, length(unique(c(FEMALE_FACTOR, MALE_FACTOR)))], start = 20, end = 340, c=100, l=60)
	names(cols)	<- tmp[, sort(unique(c(FEMALE_FACTOR, MALE_FACTOR)))]
	Rakai.plot.directed.pairs.discrete1(tmp, cols, file.path(dir, paste(run,'-phsc-directionpairs_education.pdf',sep='')), 'Educational status', w=10, h=7)	
}

RakaiFull.analyze.couples.todi.170811.NGS.success<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	infile		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170704_assemblystatus.rda'
	infile.bam	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/bam_stats_171018.rda'
	outfile.base<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/all_bams_'
	load(infile)
	load(infile.bam)
	setnames(bam.cov, 'FILE_ID', 'SID')
	setnames(bam.cov200, 'FILE_ID', 'SID')
	bam.cov200	<- subset(bam.cov200, !is.na(COV))
	setnames(bam.cov250, 'FILE_ID', 'SID')
	bam.cov250	<- subset(bam.cov250, !is.na(COV))
	setnames(bam.cov280, 'FILE_ID', 'SID')
	bam.cov280	<- subset(bam.cov280, !is.na(COV))
	
	subset(dc, PROC_STATUS=='ThoseWithoutFastqs')[, table(WTSI_STATUS)]
	#Assume sequencing failed 
	#					   16 
	dc		<- subset(dc, !is.na(SID))
	
	#	individuals with WTSI output
	dc[, length(unique(RID))]
	#	4074
	
	bam.covm	<- do.call('rbind',lapply(c(1,10,30,50), function(x)
					{
						bam.covm	<- subset(bam.cov, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']						
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.covm	<- merge(subset(dc, select=c(SID, RID)), bam.covm, by='SID')
	tmp			<- bam.covm[, list(SID=SID[which.max(HCOV)]), by='RID']
	bam.covm	<- merge(tmp, bam.covm, by=c('RID','SID'))	
	ans			<- bam.covm[, list(	N=length(RID), 
									XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
									HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
									), by='COV_MIN']
	ans[, READL:=1L]
		
	#	same query but only on short reads that are 
	#	at least 200 bp long
	#	NOTE: paired ends are not merged!!
	bam.cov200m	<- do.call('rbind',lapply(c(1,10,30, 50), function(x)
					{
						bam.covm	<- subset(bam.cov200, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']												
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.cov200m	<- merge(subset(dc, select=c(SID, RID)), bam.cov200m, by='SID')
	tmp			<- bam.cov200m[, list(SID=SID[which.max(HCOV)]), by='RID']
	bam.cov200m	<- merge(tmp, bam.cov200m, by=c('RID','SID'))
	tmp			<- bam.cov200m[, list(	N=length(RID), 
									XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
									HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
							), by='COV_MIN']
	tmp[, READL:=200L]
	ans			<- rbind(ans, tmp)	
	
	#	same query but only on short reads that are 
	#	at least 250 bp long
	#	NOTE: paired ends are not merged!!
	bam.cov250m	<- do.call('rbind',lapply(c(1,10,30, 50), function(x)
					{
						bam.covm	<- subset(bam.cov250, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']												
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.cov250m	<- merge(subset(dc, select=c(SID, RID)), bam.cov250m, by='SID')
	tmp			<- bam.cov250m[, list(SID=SID[which.max(HCOV)]), by='RID']
	bam.cov250m	<- merge(tmp, bam.cov250m, by=c('RID','SID'))
	tmp			<- bam.cov250m[, list(	N=length(RID), 
										XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
										HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
								), by='COV_MIN']
	tmp[, READL:=250L]
	ans			<- rbind(ans, tmp)	
	
	#	same query but only on short reads that are 
	#	at least 280 bp long
	#	NOTE: paired ends are not merged!!
	bam.cov280m	<- do.call('rbind',lapply(c(1,10,30, 50), function(x)
					{
						bam.covm	<- subset(bam.cov280, COV>=x)
						bam.covm	<- bam.covm[, list(HCOV=sum(REP)/9719, XCOV= sum(COV*REP)/9719), by='SID']												
						bam.covm[, COV_MIN:=x]
						bam.covm
					}))
	bam.cov280m	<- merge(subset(dc, select=c(SID, RID)), bam.cov280m, by='SID')
	tmp			<- bam.cov280m[, list(SID=SID[which.max(HCOV)]), by='RID']
	bam.cov280m	<- merge(tmp, bam.cov280m, by=c('RID','SID'))
	tmp			<- bam.cov280m[, list(	N=length(RID), 
										XCOV_MEAN=mean(XCOV), XCOV_MIN= min(XCOV), XCOV_QL= quantile(XCOV, p=0.025), XCOV_QU= quantile(XCOV, p=0.975), XCOV_MAX= max(XCOV),
										HCOV_MEAN=mean(HCOV), HCOV_MIN= min(HCOV), HCOV_QL= quantile(HCOV, p=0.025), HCOV_QU= quantile(HCOV, p=0.975), HCOV_MAX= max(HCOV)
								), by='COV_MIN']
	tmp[, READL:=280L]
	ans			<- rbind(ans, tmp)	
		
	#	make table
	ans[, XCOV_LABEL:= paste0(round(XCOV_MEAN,d=0),'x ( ', round(XCOV_QL, d=1),'x - ',round(XCOV_QU, d=0),'x )')]
	ans[, HCOV_LABEL:= paste0(round(100*HCOV_MEAN,d=1),'% (', round(100*HCOV_QL, d=1),'% - ',round(100*HCOV_QU, d=1),'%)')]
	ans[, P:= paste0( round(100*N / ans[COV_MIN==1 & READL==1, N], d=1), '%')]
	set(ans, NULL, 'COV_MIN', ans[, factor(COV_MIN, levels=c(1,10,30,50), labels=c('1X','10X','30X','50X'))])
	ans			<- subset(ans, select=c(COV_MIN, READL, N, P, HCOV_LABEL))
	setkey(ans, COV_MIN, READL)
	write.csv(ans, row.names=FALSE, file=paste0(outfile.base,'NGSoutput_info.csv'))
	
	
	
	#	roughly where are these 250 bp reads?
	bam.cov250.30	<- subset(bam.cov250, COV>=30)
	bam.cov250.30	<- merge(subset(dc, select=c(SID, RID)), bam.cov250.30, by='SID')
	tmp				<- bam.cov250.30[, list(SUM_REP=sum(REP)), by=c('SID','RID')]	
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
			scale_x_continuous(breaks=seq(0,10e3,1e3), expand=c(0,0)) +
			scale_y_continuous(breaks=seq(0,10e3,1e3), expand=c(0,0)) +
			coord_cartesian(ylim=c(0,3100)) +
			labs(x='\ngenomic position of mapped short reads', y='individuals\nwith minimum NGS output\n') +
			theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
	ggsave(file=paste0(outfile.base,'horizontal_coverage_histogram_minNGSoutput_reads250.pdf'), w=10, h=4)
}


RakaiFull.analyze.couples.todi.170811.computing.effort.couples<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_withmetadata.rda'	
	load(infile)
		
	load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170704_stagetwo.rda')
	tmp			<- subset(pty.runs, PTY_RUN!=1)
	set(tmp, NULL, 'PID', NULL)
	load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170704_stagethree.rda')
	pty.runs	<- rbind(tmp, pty.runs)
	
	tmp			<- unique(subset(rp, select=c(FEMALE_RID, MALE_RID)))
	ptyc		<- tmp[, {
				#FEMALE_RID<- 'A090383'; MALE_RID<- 'D110040'	
				z	<- subset(pty.runs, RID%in%c(FEMALE_RID,MALE_RID), select=c(PTY_RUN, RID))
				z	<- z[, list(HAS_COUPLE= length(unique(RID))>1), by='PTY_RUN']
				z	<- subset(z, HAS_COUPLE)[, PTY_RUN]
				list(PTY_RUN=z, N=length(z))				
			}, by=c('FEMALE_RID','MALE_RID')]
	unique(ptyc, by=c('FEMALE_RID','MALE_RID'))[, table(N)]
	#N
  	#	0   1   2   3   4 
  	#	2 117 346  15   6 
	
	#	couples are in how many runs?
	ptycr	<- unique(subset(ptyc, select=PTY_RUN))
	#	324
	
	#	how many trees in these runs
	indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30'
	tmp		<- data.table(F=list.files(indir, pattern='trees.rda$', full.names=TRUE))
	tmp[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(F)))]
	ptycr	<- merge(ptycr, tmp, by='PTY_RUN')
	ptycr	<- ptycr[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30/ptyr10_trees.rda'
				load(F)
				list(N_PHY=length(phs), MIN_TAXA=min(sapply(phs, Ntip)), MAX_TAXA=max(sapply(phs, Ntip)))
			}, by='PTY_RUN']
	
	ptycr[, c( sum(N_PHY), min(MIN_TAXA), max(MAX_TAXA))]
	#	87731    39 23438
	
}
	
RakaiFull.analyze.couples.todi.170811.demographic.table<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_withmetadata.rda'
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_"
	load(infile)
	
	setkey(rca, MALE_RID, FEMALE_RID)
	subset(rca[, list(N=length(PTY_RUN)), by=c('MALE_RID','FEMALE_RID')], N>1)
	unique(rca, by=c('MALE_RID','FEMALE_RID'))
	
	#	table
	#	location female (community type)
	group	<- 'FEMALE_COMM_TYPE'
	tmp		<- melt(rca[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='all']
	ans		<- copy(tmp)
	tmp		<- melt(subset(rca, !grepl('insufficient deep sequence data', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='data']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('couple most likely not a pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='not a pair']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('couple ambiguous if pair or not pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='ambiguous if pair']	
	ans		<- rbind(ans, tmp)	
	tmp		<- melt(subset(rca, grepl('couple most likely a pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='pair']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('with resolved direction', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='pair direction']
	ans		<- rbind(ans, tmp)
	#	year first concordant positive
	rca[, FIRSTCONCPOS:=as.character(round(pmax(MALE_FIRSTPOSDATE,FEMALE_FIRSTPOSDATE), d=0))]
	set(rca, rca[, which(FIRSTCONCPOS<'2012')], 'FIRSTCONCPOS', '<2012')
	group	<- 'FIRSTCONCPOS'
	tmp		<- melt(rca[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='all']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, !grepl('insufficient deep sequence data', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='data']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('couple most likely not a pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='not a pair']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('couple ambiguous if pair or not pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='ambiguous if pair']	
	ans		<- rbind(ans, tmp)	
	tmp		<- melt(subset(rca, grepl('couple most likely a pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='pair']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('with resolved direction', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='pair direction']
	ans		<- rbind(ans, tmp)
	rca[, COUP_SC2:= 'sero-positive']
	set(rca, rca[, which(grepl('seroinc',COUP_SC))], 'COUP_SC2', 'sero-incident')
	set(rca, rca[, which(grepl('M->F|F->M',COUP_SC))], 'COUP_SC2', 'sero-discordant')
	group	<- 'COUP_SC2'
	tmp		<- melt(rca[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='all']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, !grepl('insufficient deep sequence data', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='data']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('couple most likely not a pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='not a pair']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('couple ambiguous if pair or not pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='ambiguous if pair']	
	ans		<- rbind(ans, tmp)	
	tmp		<- melt(subset(rca, grepl('couple most likely a pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='pair']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('with resolved direction', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='pair direction']
	ans		<- rbind(ans, tmp)
	#	age of female at time first conc pos
	rca[, FEMALE_AGE:=cut(pmax(MALE_FIRSTPOSDATE,FEMALE_FIRSTPOSDATE)-FEMALE_BIRTHDATE, breaks=c(0,20,25,30,35,40,45,80),labels=c('<20','20-24','25-29','30-34','35-39','40-44','45+'))]
	group	<- 'FEMALE_AGE'
	tmp		<- melt(rca[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='all']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, !grepl('insufficient deep sequence data', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='data']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('couple most likely not a pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='not a pair']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('couple ambiguous if pair or not pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='ambiguous if pair']	
	ans		<- rbind(ans, tmp)	
	tmp		<- melt(subset(rca, grepl('couple most likely a pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='pair']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('with resolved direction', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='pair direction']
	ans		<- rbind(ans, tmp)
	#	age of male at time first conc pos
	rca[, MALE_AGE:=cut(pmax(MALE_FIRSTPOSDATE,FEMALE_FIRSTPOSDATE)-MALE_BIRTHDATE, breaks=c(0,20,25,30,35,40,45,80),labels=c('<20','20-24','25-29','30-34','35-39','40-44','45+'))]
	group	<- 'MALE_AGE'
	tmp		<- melt(rca[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='all']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, !grepl('insufficient deep sequence data', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='data']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('couple most likely not a pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='not a pair']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('couple ambiguous if pair or not pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='ambiguous if pair']	
	ans		<- rbind(ans, tmp)	
	tmp		<- melt(subset(rca, grepl('couple most likely a pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='pair']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('with resolved direction', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='pair direction']
	ans		<- rbind(ans, tmp)
	#	age difference
	rca[, AGEDIFF:= cut(FEMALE_BIRTHDATE-MALE_BIRTHDATE, breaks=c(-30,-5,-1,1,5,30), labels=c('female >5 yrs older','female 1-5 yrs older','male/female +-1 yrs','male 1-5 yrs older','male >5 yrs older'))]
	group	<- 'AGEDIFF'
	tmp		<- melt(rca[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='all']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, !grepl('insufficient deep sequence data', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='data']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('couple most likely not a pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='not a pair']	
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('couple ambiguous if pair or not pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='ambiguous if pair']	
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('couple most likely a pair', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='pair']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rca, grepl('with resolved direction', SELECT))[, list(N=length(FEMALE_RID)), by=group], id.vars='N', measure.vars=c(group), value.name='FACTOR', variable.name='GROUP')
	tmp[, WHO:='pair direction']
	ans		<- rbind(ans, tmp)
	
	ans		<- merge(ans, ans[, list(P=N/sum(N), FACTOR=FACTOR), by=c('GROUP','WHO')],by=c('GROUP','WHO','FACTOR'))
	ans[, LABEL:= paste0(N,' (',round(P*100,d=0),'%)')]
	
	tmp		<- dcast.data.table(ans, GROUP+FACTOR~WHO, value.var='LABEL')
	set(tmp, tmp[, which(is.na(FACTOR))], 'FACTOR', 'Unknown')
	set(tmp, NULL, 'FACTOR', tmp[, factor(FACTOR, levels=c("agrarian","trading","fisherfolk",                         
									"<2012","2012","2013","2014",
									"sero-positive","sero-incident","sero-discordant",     
									"<20","20-24","25-29","30-34","35-39","40-44","45+","Unknown",                 	               
									"female >5 yrs older","female 1-5 yrs older","male/female +-1 yrs","male 1-5 yrs older","male >5 yrs older"))])	
	set(tmp, NULL, 'GROUP', tmp[, factor(GROUP, levels=c('FIRSTCONCPOS','FEMALE_COMM_TYPE','COUP_SC2','FEMALE_AGE','MALE_AGE','AGEDIFF'))])
	setkey(tmp, GROUP, FACTOR)
	tmp		<- subset(tmp, select=c('GROUP','FACTOR','all','data','not a pair','ambiguous if pair','pair','pair direction'))
	write.csv(tmp, row.names=FALSE, file=paste0(outfile.base,'couples_demographics.csv'))
	
	
}

RakaiFull.analyze.trmpairs.todi.170522.demographic.table<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	#infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170522/todi_couples_170522_withmetadata.rda'		
	#outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170522/todi_couples_170522_"	
	load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/community_sampling_170522.rda')
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_"	
	rsmpl[, AGE_C:= rsmpl[, as.character(cut(AGEYRS, right=FALSE, breaks=c(seq(15,45,5),65), labels=paste0(seq(15,50,5)[-8], '-', seq(15,50,5)[-1]-1)))]]
	#	make table		cols: agrarian fisherfolk trading
	#	rows by mf
	#	eligible 
	#	participated
	#	HIV+
	#	submitted for seq 
	#	with suff data	
	
	tmp		<- melt(rsmpl[, list(N=length(RID)), by=c('COMM_TYPE','SEX')], id.vars=c('COMM_TYPE','N'), measure.vars='SEX', variable.name='GROUP', value.name='FACTOR')
	tmp[, WHO:='participated']
	ans		<- copy(tmp)
	tmp		<- melt(subset(rsmpl, HIV==1)[, list(N=length(RID)), by=c('COMM_TYPE','SEX')], id.vars=c('COMM_TYPE','N'), measure.vars='SEX', variable.name='GROUP', value.name='FACTOR')
	tmp[, WHO:='infected']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rsmpl, BAM_OUTPUT==1)[, list(N=length(RID)), by=c('COMM_TYPE','SEX')], id.vars=c('COMM_TYPE','N'), measure.vars='SEX', variable.name='GROUP', value.name='FACTOR')
	tmp[, WHO:='successfully_sequenced']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rsmpl, MIN_PNG_OUTPUT==1)[, list(N=length(RID)), by=c('COMM_TYPE','SEX')], id.vars=c('COMM_TYPE','N'), measure.vars='SEX', variable.name='GROUP', value.name='FACTOR')
	tmp[, WHO:='enough_sequence_phyloscanner']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(rsmpl[, list(N=length(RID)), by=c('COMM_TYPE','AGE_C')], id.vars=c('COMM_TYPE','N'), measure.vars='AGE_C', variable.name='GROUP', value.name='FACTOR')
	tmp[, WHO:='participated']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rsmpl, HIV==1)[, list(N=length(RID)), by=c('COMM_TYPE','AGE_C')], id.vars=c('COMM_TYPE','N'), measure.vars='AGE_C', variable.name='GROUP', value.name='FACTOR')
	tmp[, WHO:='infected']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rsmpl, BAM_OUTPUT==1)[, list(N=length(RID)), by=c('COMM_TYPE','AGE_C')], id.vars=c('COMM_TYPE','N'), measure.vars='AGE_C', variable.name='GROUP', value.name='FACTOR')
	tmp[, WHO:='successfully_sequenced']
	ans		<- rbind(ans, tmp)
	tmp		<- melt(subset(rsmpl, MIN_PNG_OUTPUT==1)[, list(N=length(RID)), by=c('COMM_TYPE','AGE_C')], id.vars=c('COMM_TYPE','N'), measure.vars='AGE_C', variable.name='GROUP', value.name='FACTOR')
	tmp[, WHO:='enough_sequence_phyloscanner']
	ans		<- rbind(ans, tmp)
	ans		<- rbind(ans, ans[, list(FACTOR=paste0('ALL_',GROUP), N=sum(N)), by=c('COMM_TYPE','WHO','GROUP')])
	set(ans, NULL, 'COMM_TYPE', ans[, factor(COMM_TYPE, levels=c("agrarian","trading","fisherfolk"))])
	set(ans, NULL, 'FACTOR', ans[, factor(FACTOR, levels=c("M","F","ALL_SEX","15-19","20-24","25-29","30-34","35-39","40-44","45-49",'Unknown',"ALL_AGE_C"), labels=c('Male','Female','Total_Sex',"15-19","20-24","25-29","30-34","35-39","40-44","45-49",'Unknown','Total_Age'))])
	set(ans, NULL, 'WHO', ans[, factor(WHO, levels=c("participated","infected","successfully_sequenced","enough_sequence_phyloscanner"))])
	
	setkey(ans, WHO, COMM_TYPE, GROUP, FACTOR)
	
	tmp		<- dcast.data.table(ans, COMM_TYPE+FACTOR~WHO, value.var='N')
	setkey(tmp, COMM_TYPE, FACTOR)
	tmp[, sequenced_p:= enough_sequence_phyloscanner/infected]
	tmp[, infected_l:= paste0(infected,' (',100*round(infected/participated, d=2),'%)')]
	tmp[, sequenced_l:= paste0(enough_sequence_phyloscanner,' (',100*round(sequenced_p, d=2),'%)')]
	
	write.csv(subset(tmp, select=c(COMM_TYPE, FACTOR, participated, infected_l, sequenced_l)), row.names=FALSE, file=paste0(outfile.base,'RCCS_demographics.csv'))
}

RakaiFull.analyze.concordance.phylo.monogamous.170522<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_concordance_170516_"
	#
	#	load couples
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_info.rda")
	rps		<- unique(subset(rp, select=c(FEMALE_RID, FEMALE_MARSTAT, FEMALE_ALC, FEMALE_SEXP1YR, FEMALE_SEXP1OUT, MALE_RID, MALE_MARSTAT, MALE_ALC, MALE_SEXP1YR, MALE_SEXP1OUT)))
	setnames(rps, gsub('MALE','ID1',gsub('FEMALE','ID2',colnames(rps))))
	tmp		<- copy(rps)
	setnames(tmp, gsub('DUMMY','ID2',gsub('ID2','ID1',gsub('ID1','DUMMY',colnames(tmp)))))
	rps		<- rbind(rps, tmp)
	setnames(rps, c('ID1_RID','ID2_RID'),c('ID1','ID2'))	
	#
	#	select likely pairs 
	confidence.cut	<- 0.5			
	indir		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run'	
	infiles		<- data.table(F=list.files(indir, pattern='close[0-9]+.rda$', full.names=TRUE))
	rtp			<- infiles[, {
				#F			<- '/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_ff_170516_close10.rda'
				load(F)
				#	select couples from those runs with most data
				pairings		<- unique(subset(rplkl, select=c(ID1,ID2,PTY_RUN,NEFF)))[, list(PTY_RUN=PTY_RUN[which.max(NEFF)]), by=c('ID2','ID1')]
				pairings		<- merge(pairings, rps, by=c('ID2','ID1'))
				rplkl2			<- merge(pairings, rplkl, by=c('ID2','ID1','PTY_RUN'))
				#	likely transmission pairs, using distance
				rtp				<- subset(rplkl2, GROUP=='TYPE_PAIR_DI')[, list(TYPE=TYPE[which.max(KEFF)]), by=c('ID1','ID2','PTY_RUN')]
				rtp				<- merge(rtp, subset(rplkl2, GROUP=='TYPE_PAIR_DI'), by=c('ID1','ID2','PTY_RUN','TYPE'))
				rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				ans				<- copy(rtp)
				#	likely transmission pairs, using distance+topoloy
				group			<- 'TYPE_PAIR_TODI2'
				#group			<- 'TYPE_PAIR_TODI'
				rtp				<- subset(rplkl2, GROUP==group)[, list(TYPE=TYPE[which.max(KEFF)]), by=c('ID1','ID2','PTY_RUN')]
				rtp				<- merge(rtp, subset(rplkl2, GROUP==group), by=c('ID1','ID2','PTY_RUN','TYPE'))
				rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				ans				<- rbind(rtp, ans)
				ans[, SELECT:=NA_character_]
				set(ans, ans[, which(grepl('close|likely pair',TYPE) & POSTERIOR_SCORE>confidence.cut)], 'SELECT', 'couple linked')
				set(ans, ans[, which(grepl('distant|other',TYPE) & POSTERIOR_SCORE>confidence.cut)], 'SELECT', 'couple unlinked')
				
				#	likely transmission pairs, using distance+topoloy, direction resolved
				#	for direct transmission, not enough to know couple unlinked, we need to have the transmitter in the data
				#	want all directed pairs  
				rtp				<- subset(rplkl, GROUP==group)[, list(TYPE=TYPE[which.max(KEFF)]), by=c('ID1','ID2','PTY_RUN')]
				rtp				<- subset(rtp, TYPE=='likely pair')				
				rtp				<- merge(rtp, subset(rplkl, GROUP==group), by=c('ID1','ID2','PTY_RUN','TYPE'))
				rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				rtp				<- subset(rtp, POSTERIOR_SCORE>confidence.cut)
				rtp				<- merge(subset(rtp, select=c('ID1','ID2','PTY_RUN')), subset(rplkl, GROUP=='TYPE_DIR_TODI3'), by=c('ID1','ID2','PTY_RUN'))
				rtp				<- rtp[, list(TYPE=TYPE[which.max(KEFF)]), by=c('ID1','ID2','PTY_RUN')]
				rtp				<- subset(rtp, TYPE!='ambiguous')
				rtp				<- merge(rtp, subset(rplkl, GROUP=='TYPE_DIRSCORE_TODI3'), by=c('ID1','ID2','PTY_RUN','TYPE'))
				rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				rtp				<- subset(rtp, POSTERIOR_SCORE>confidence.cut)
				#	now select directed pairs with recipient in couple
				tmp				<- subset(rtp, TYPE=='12' & ID2_COUPLE==1)
				tmp[, DUMMY:= 1:nrow(tmp)]
				tmp[, SELECT:= 'couple recipient infected from extra-marital partner']
				tmp2			<- merge(subset(tmp, select=c(ID1,ID2,DUMMY)), subset(rps, select=c(ID1,ID2)), by=c('ID1','ID2'))[, DUMMY]
				set(tmp, tmp2, 'SELECT', 'couple recipient infected from spouse')
				tmp[, DUMMY:= NULL]
				rtp				<- subset(rtp, TYPE=='21' & ID1_COUPLE==1)
				rtp[, DUMMY:= 1:nrow(rtp)]
				rtp[, SELECT:= 'couple recipient infected from extra-marital partner']
				tmp2			<- merge(subset(rtp, select=c(ID1,ID2,DUMMY)), subset(rps, select=c(ID1,ID2)), by=c('ID1','ID2'))[, DUMMY]
				set(rtp, tmp2, 'SELECT', 'couple recipient infected from spouse')
				rtp[, DUMMY:= NULL]
				rtp				<- rbind(rtp, tmp)
				#	merge sexual behaviour data, need time first conc pos
				tmp		<- unique(subset(rd, select=c(RID, FIRSTPOSVIS)))
				setnames(tmp, colnames(tmp), gsub('_RID','',paste0('ID1_',colnames(tmp))))
				rtp		<- merge(rtp, tmp, by='ID1')	
				setnames(tmp, colnames(tmp), gsub('ID1','ID2',colnames(tmp)))
				rtp		<- merge(rtp, tmp, by='ID2')	
				rtp[, VISIT_FIRSTCONCPOS:=pmax(ID1_FIRSTPOSVIS,ID2_FIRSTPOSVIS)]
				# 	now merge data for ID1 from closest visit
				tmp		<- unique(subset(rh, select=c(RID, VISIT)))
				setnames(tmp, 'RID', 'ID1')
				tmp		<- unique(merge(rtp, tmp, by='ID1')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('ID1','VISIT_FIRSTCONCPOS')])
				setnames(tmp, 'ID1', 'RID')
				tmp2	<- unique(subset(rh, select=c(RID,VISIT,MARSTAT,ALC,SEXP1YR,SEXP1OUT)), by=c('RID','VISIT'))
				tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
				set(tmp, NULL, 'VISIT', NULL)
				setnames(tmp, colnames(tmp), gsub('ID1_VISIT_FIRSTCONCPOS','VISIT_FIRSTCONCPOS',gsub('_RID','',paste0('ID1_',colnames(tmp)))))					
				rtp		<- merge(rtp, tmp, by=c('ID1','VISIT_FIRSTCONCPOS'), all.x=1)
				# 	now merge data for ID2 from closest visit
				tmp		<- unique(subset(rh, select=c(RID, VISIT)))
				setnames(tmp, 'RID', 'ID2')
				tmp		<- unique(merge(rtp, tmp, by='ID2')[, list(VISIT= VISIT[which.min(abs(VISIT-VISIT_FIRSTCONCPOS))]), by=c('ID2','VISIT_FIRSTCONCPOS')])
				setnames(tmp, 'ID2', 'RID')
				tmp2	<- unique(subset(rh, select=c(RID,VISIT,MARSTAT,ALC,SEXP1YR,SEXP1OUT)), by=c('RID','VISIT'))
				tmp		<- merge(tmp, tmp2, by=c('RID','VISIT'))
				set(tmp, NULL, 'VISIT', NULL)
				setnames(tmp, colnames(tmp), gsub('ID2_VISIT_FIRSTCONCPOS','VISIT_FIRSTCONCPOS',gsub('_RID','',paste0('ID2_',colnames(tmp)))))					
				rtp		<- merge(rtp, tmp, by=c('ID2','VISIT_FIRSTCONCPOS'), all.x=1)
				set(rtp, NULL, c('VISIT_FIRSTCONCPOS','ID1_FIRSTPOSVIS','ID2_FIRSTPOSVIS'), NULL)
				#	
				ans		<- rbind(rtp, ans)
				ans		<- subset(ans, !is.na(SELECT))
				ans
			}, by='F']
	#		
	#	select monogamous couples
	#	count linked/unlinked couples
	cs		<- subset(rtp, grepl('couple linked|couple unlinked',SELECT) & grepl('Monogamous',ID1_MARSTAT) & grepl('Monogamous',ID2_MARSTAT))
	cs[, EXTRA_PARTNERS:='none']
	set(cs, cs[, which(grepl('casual',ID1_MARSTAT) | grepl('casual',ID2_MARSTAT))], 'EXTRA_PARTNERS', 'at least one')
	csn		<- cs[, list(N=length(ID1)), by=c('CLOSE_BRL','GROUP','SELECT','EXTRA_PARTNERS')]
	#
	#	count couple recipients with / without external partners
	tmp		<- subset(rtp, grepl('couple recipient',SELECT) & TYPE=='12' & ID2_COUPLE==1 & grepl('Monogamous',ID2_MARSTAT))
	tmp[, EXTRA_PARTNERS:='recipient no']
	set(tmp, tmp[, which(grepl('casual',ID2_MARSTAT))], 'EXTRA_PARTNERS', 'recipient yes')
	cs		<- subset(rtp, grepl('couple recipient',SELECT) & TYPE=='21' & ID1_COUPLE==1 & grepl('Monogamous',ID1_MARSTAT))
	cs[, EXTRA_PARTNERS:='recipient no']
	set(cs, cs[, which(grepl('casual',ID1_MARSTAT))], 'EXTRA_PARTNERS', 'recipient yes')
	cs		<- rbind(cs, tmp)
	tmp		<- cs[, list(N=length(ID1)), by=c('CLOSE_BRL','GROUP','SELECT','EXTRA_PARTNERS')]
	csn		<- rbind(csn, tmp)
	#
	set(csn, NULL, 'GROUP', csn[, factor(GROUP, levels=c(	'TYPE_PAIR_DI','TYPE_PAIR_TODI2','TYPE_DIRSCORE_TODI3'),
												labels=c(	'at least one spouse reporting extra-marital partners\nif phylogenetically unlinked\ndefined by distance between within-host trees',
															'at least one spouse reporting extra-marital partners\nif phylogenetically unlinked\ndefined by distance and topology between within-host trees',
															'likely recipient reporting extra-marital partners\nif phylogenetic evidence for acquiring infection externally\ndefined by distance and topology between within-host trees'))])
	set(csn, NULL, 'SELECT', csn[, factor(SELECT,levels=c('couple unlinked','couple linked','couple recipient infected from extra-marital partner','couple recipient infected from spouse'))])
	set(csn, NULL, 'EXTRA_PARTNERS', csn[, factor(EXTRA_PARTNERS, levels=c('at least one','none','recipient yes','recipient no'))])
	setkey(csn, CLOSE_BRL, GROUP, SELECT, EXTRA_PARTNERS)
	csn[, LABEL:=paste0(SELECT, '-',EXTRA_PARTNERS)]
	#	odds ratios	
	tmp	<- csn[, {
				z<- log( (N[1]/N[2])/(N[3]/N[4]) )				
				list(LOR=z, LOR_CL=z-qnorm(0.975)*sqrt(sum(1/N)), LOR_CU=z+qnorm(0.975)*sqrt(sum(1/N)) )
			}, by=c('CLOSE_BRL','GROUP')]	
	ggplot(tmp, aes(x=factor(CLOSE_BRL), y=exp(LOR), ymin=exp(LOR_CL), ymax=exp(LOR_CU), colour=GROUP)) +
			geom_point(position=position_dodge(width=0.4), size=2) + 
			geom_errorbar(position=position_dodge(width=0.4), width=0.2) +
			scale_y_continuous() +			
			#scale_colour_manual(values=c('only distance between within-host trees'=rev(brewer.pal(11, 'PuOr'))[3],'distance and topology between within-host trees'=brewer.pal(11, 'PuOr')[2])) +
			theme_bw() + theme(legend.position='bottom') +
			labs(	x='\ndistance between within-host trees of different individuals\n(threshold for calling a likely transmission pair in subst/site)',
					y='odds ratio\nfor being a likely transmission pair and\nno self-reported extra-marital partners\n',
					colour='information used')
	ggsave(file=paste0(outfile.base,'oddsratio.pdf'), w=4, h=5)
		
	#ggplot(tmp, aes(x=GROUP, fill=LABEL, y=N)) +
	#		geom_bar(stat='identity', position='stack') +
	#		theme_bw() + theme(legend.position='bottom', panel.spacing=grid::unit(0, "lines")) +
	#		facet_grid(~CLOSE_BRL)
	
	#tmp	<- subset(cs, select=c(ID1,ID2,GROUP,CLOSE_BRL,SELECT_AS_LKL_PAIR,EXTRA_PARTNERS))
	setkey(tmp, CLOSE_BRL, GROUP, SELECT_AS_LKL_PAIR, EXTRA_PARTNERS)
	tmp2	<- tmp[, {
				pe		<- sum(N[c(1,2)]) / sum(N) * sum(N[c(1,3)]) / sum(N) + sum(N[c(3,4)]) / sum(N) * sum(N[c(2,4)]) / sum(N)
				pc		<- sum(N[c(1,4)]) / sum(N)
				list(N=sum(N), PE=pe, PC=pc, COHEN_KAPPA=(pc-pe)/(1-pe))
			}, by=c('CLOSE_BRL','GROUP')]
	tmp2[, COHEN_KAPPA_CL:= COHEN_KAPPA-qnorm(0.975)*sqrt( PC*(1-PC)/(1-PE)/(1-PE)/N )]
	tmp2[, COHEN_KAPPA_CU:= COHEN_KAPPA+qnorm(0.975)*sqrt( PC*(1-PC)/(1-PE)/(1-PE)/N )]	
	set(tmp2, NULL, 'GROUP', tmp2[, factor(GROUP, levels=c('TYPE_PAIR_DI','TYPE_PAIR_TODI2'), labels=c('only distance between within-host trees','distance and topology between within-host trees'))])
	
	ggplot(tmp2, aes(x=factor(CLOSE_BRL), y=COHEN_KAPPA, ymin=COHEN_KAPPA_CL, ymax=COHEN_KAPPA_CU, colour=GROUP)) +
			geom_point(position=position_dodge(width=0.4), size=2) + 
			geom_errorbar(position=position_dodge(width=0.4), width=0.2) +
			scale_y_continuous(labels=scales::percent) +			
			scale_colour_manual(values=c('only distance between within-host trees'=rev(brewer.pal(11, 'PuOr'))[3],'distance and topology between within-host trees'=brewer.pal(11, 'PuOr')[2])) +
			theme_bw() + theme(legend.position='bottom') +
			labs(	x='\ndistance between within-host trees of different individuals\n(threshold for calling a likely transmission pair in subst/site)',
					y='concordance\nbetween phylogenetically likely transmitters\nand self reported extra-marital partners\n(Cohens kappa)\n',
					colour='information used')
	ggsave(file=paste0(outfile.base,'CohensKappa.pdf'), w=4, h=5)
}

RakaiFull.analyze.ffpairs.todi.170522<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170516_cl3_allwindows.rda"		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_ff_170516_"
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170610_cl3_prior23_allwindows.rda"		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_ff_170610_cl3_prior23_"	
	load(infile)
	
	#	add to rpw and rplkl if ID is part of couple	
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_info.rda")
	rpw[, ID1_COUPLE:= 0L]
	set(rpw, rpw[, which(ID1%in%rp$MALE_RID)], 'ID1_COUPLE', 1L)
	set(rpw, rpw[, which(ID1%in%rp$FEMALE_RID)], 'ID1_COUPLE', 1L)	
	rpw[, ID2_COUPLE:= 0L]
	set(rpw, rpw[, which(ID2%in%rp$MALE_RID)], 'ID2_COUPLE', 1L)
	set(rpw, rpw[, which(ID2%in%rp$FEMALE_RID)], 'ID2_COUPLE', 1L)		
	rplkl[, ID1_COUPLE:= 0L]
	set(rplkl, rplkl[, which(ID1%in%rp$MALE_RID)], 'ID1_COUPLE', 1L)
	set(rplkl, rplkl[, which(ID1%in%rp$FEMALE_RID)], 'ID1_COUPLE', 1L)	
	rplkl[, ID2_COUPLE:= 0L]
	set(rplkl, rplkl[, which(ID2%in%rp$MALE_RID)], 'ID2_COUPLE', 1L)
	set(rplkl, rplkl[, which(ID2%in%rp$FEMALE_RID)], 'ID2_COUPLE', 1L)	
	
	
	#	per window:
	#	make window topology assignments: TYPE_PAIR_TO ancestral/intermingled; withintermediate; other (include adjacent other here)
	#	make boxplot
	cols.typet			<- c(brewer.pal(11, 'PuOr')[c(2,4)], rev(brewer.pal(11, 'PuOr'))[c(3,4)], rev(brewer.pal(11, 'RdGy'))[4])
	names(cols.typet)	<- c("ancestral", "intermingled", "with intermediate", 'sibling', "disconnected")
	#
	rpw2	<- subset(rpw, GROUP=='TYPE_BASIC')
	rpw2[, TYPE_TO:= 'disconnected']
	set(rpw2, rpw2[,which(grepl('chain', TYPE) & grepl('no intermediate', TYPE))], 'TYPE_TO', 'ancestral')
	set(rpw2, rpw2[,which(grepl('intermingled', TYPE) & grepl('no intermediate', TYPE))], 'TYPE_TO', 'intermingled')
	set(rpw2, rpw2[,which(grepl('with intermediate', TYPE))], 'TYPE_TO', 'with intermediate')
	set(rpw2, rpw2[,which(ADJACENT & grepl('other', TYPE))], 'TYPE_TO', 'sibling')
	rpw2[, PATRISTIC_DISTANCE_LOG:= log10(PATRISTIC_DISTANCE)]
	rpw2[, PATRISTIC_DISTANCE_LOGC:= cut(PATRISTIC_DISTANCE_LOG, breaks=log10(c(1e-12, 0.005, 0.0075, 0.01, 0.02, 0.035, 0.08, 2e3)),labels=c('<0.5%', '0.5%-0.75%', '0.75%-1%', '1%-2%', '2%-3.5%', '3.5%-8%', '>8%'))]
	rpw2[, PATRISTIC_DISTANCE_LOGC2:= cut(PATRISTIC_DISTANCE_LOG, breaks=log10(c(1e-12, 0.035, 0.08, 2e3)),labels=c('<3.5%', '3.5%-8%', '>8%'))]	
	rpw2[, TYPE_PAIR:= tolower(paste0(ID1_SEX, ID2_SEX))]
	set(rpw2, rpw2[, which(TYPE_PAIR%in%c('fm','mf'))], 'TYPE_PAIR', 'hsx')
	rpw2	<- subset(rpw2, TYPE_PAIR!='mm')	
	rpw2	<- subset(rpw2, ID1_COUPLE==1L | ID2_COUPLE==1L)
	tmp		<- rpw2[, list(N=length(W_FROM)), by=c('TYPE_PAIR','PATRISTIC_DISTANCE_LOGC','TYPE_TO')]
	#
	#	plot topology assignment by distance across all windows
	#
	set(tmp, NULL, 'TYPE_TO', tmp[, factor(TYPE_TO, levels=c("ancestral","intermingled", "with intermediate", 'sibling', "disconnected"))])
	ggplot(tmp, aes(x=PATRISTIC_DISTANCE_LOGC, y=N, fill=TYPE_TO)) +
			geom_bar(stat='identity', position='stack') +
			#scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), expand=c(0,0)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			facet_grid(~TYPE_PAIR) +
			labs(x='\npatristic distance between within-host clades of different individuals\n(subst/site)', y='genomic windows\n', fill='pairwise\ntopological relationship')
	ggsave(file=paste0(outfile.base, 'topodist_counts.pdf'), w=10, h=5)
	ggplot(tmp, aes(x=PATRISTIC_DISTANCE_LOGC, y=N, fill=TYPE_TO)) +
			geom_bar(stat='identity', position='fill') +
			scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), expand=c(0,0)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			facet_grid(~TYPE_PAIR) +
			labs(x='\npatristic distance between within-host clades of different individuals\n(subst/site)', y='genomic windows\n', fill='pairwise\ntopological relationship')
	ggsave(file=paste0(outfile.base, 'topodist_prop.pdf'), w=10, h=5)
	#
	#	wash out and look at MLE and mean distance
	#
	rplkl2	<- subset(rplkl, GROUP=='TYPE_BASIC')
	rplkl2[, TYPE_TO:= 'disconnected']
	set(rplkl2, rplkl2[,which(grepl('no intermediate', TYPE) & grepl('chain', TYPE))], 'TYPE_TO', 'ancestral')
	set(rplkl2, rplkl2[,which(grepl('no intermediate', TYPE) & grepl('intermingled', TYPE))], 'TYPE_TO', 'intermingled')
	set(rplkl2, rplkl2[,which(grepl('with intermediate', TYPE))], 'TYPE_TO', 'with intermediate')
	set(rplkl2, rplkl2[,which(grepl('no intermediate', TYPE) & grepl('other', TYPE))], 'TYPE_TO', 'sibling')
	rplkl2	<- rplkl2[, list(N=N[1], NEFF=NEFF[1], K=sum(K), KEFF=sum(KEFF)), by=c('ID2','ID1','PTY_RUN','ID1_SEX','ID2_SEX','TYPE_TO','ID1_COUPLE','ID2_COUPLE')]
	#	select pairings from those runs with most data
	pairings<- unique(subset(rplkl2, select=c(ID1,ID2,PTY_RUN,NEFF)))[, list(PTY_RUN=PTY_RUN[which.max(NEFF)]), by=c('ID2','ID1')]
	rplkl2	<- merge(pairings, rplkl2, by=c('ID2','ID1','PTY_RUN'))
	#	determine most likely state
	rplkl2	<- rplkl2[, {
				z<- which.max(KEFF)
				list(N=N[1], NEFF=NEFF[1], TYPE_TO=TYPE_TO[z], K=K[z], KEFF=KEFF[z])
			}, by=c('ID2','ID1','PTY_RUN','ID1_SEX','ID2_SEX','ID1_COUPLE','ID2_COUPLE')]
	#	select hsx and ff
	rplkl2[, TYPE_PAIR:= tolower(paste0(ID1_SEX, ID2_SEX))]
	set(rplkl2, rplkl2[, which(TYPE_PAIR%in%c('fm','mf'))], 'TYPE_PAIR', 'hsx')
	rplkl2	<- subset(rplkl2, TYPE_PAIR!='mm')
	#	select touch couple
	rplkl2	<- subset(rplkl2, ID1_COUPLE==1L | ID2_COUPLE==1L)
	#	deselect RUN==1 for now
	rplkl2	<- subset(rplkl2, PTY_RUN!=1) 
	#	add average distance
	tmp		<- subset(rpw, GROUP=='TYPE_BASIC')[, list(PD_MEAN=mean(PATRISTIC_DISTANCE), PD_IQL=quantile(PATRISTIC_DISTANCE,p=0.25), PD_IQU=quantile(PATRISTIC_DISTANCE,p=0.75)), by=c('ID2','ID1','PTY_RUN')]
	rplkl2	<- merge(rplkl2, tmp, by=c('ID2','ID1','PTY_RUN'))	
	rplkl2[, PD_MEAN_LOGC:= cut(log10(PD_MEAN), breaks=log10(c(1e-12, 0.005, 0.0075, 0.01, 0.02, 0.035, 0.08, 2e3)),labels=c('<0.5%', '0.5%-0.75%', '0.75%-1%', '1%-2%', '2%-3.5%', '3.5%-8%', '>8%'))]
	rplkl2[, PD_MEAN_LOGC:= cut(log10(PD_MEAN), breaks=log10(c(1e-12, 0.005, 0.01, 0.02, 0.035, 0.08, 2e3)),labels=c('<0.5%', '0.5%-1%', '1%-2%', '2%-3.5%', '3.5%-8%', '>8%'))]
	rplkl2[, PD_MEAN_LOGC2:= cut(log10(PD_MEAN), breaks=log10(c(1e-12, 0.035, 0.08, 2e3)),labels=c('<3.5%', '3.5%-8%', '>8%'))]	
	#	plot	
	set(rplkl2, NULL, 'TYPE_TO', rplkl2[, factor(TYPE_TO, levels=c("ancestral","intermingled", "with intermediate", 'sibling', "disconnected"))])
	ggplot(rplkl2, aes(x=PD_MEAN_LOGC, fill=TYPE_TO)) +
			geom_bar(position='stack') +			
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			scale_y_continuous(expand=c(0,0)) +
			facet_grid(~TYPE_PAIR) +
			coord_cartesian(ylim=c(0,200)) +
			labs(x='\npatristic distance between within-host clades of different individuals\n(average across genomic windows, in subst/site)', y='most likely topology classification\n(#pairings)\n', fill='pairwise\ntopological relationship')
	ggsave(file=paste0(outfile.base, 'topodist_counts_MLE.pdf'), w=10, h=5)
	ggplot(rplkl2, aes(x=PD_MEAN_LOGC, fill=TYPE_TO)) +
			geom_bar(position='stack') +			
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			scale_y_continuous(expand=c(0,0)) +
			facet_grid(~TYPE_PAIR) +
			coord_cartesian(ylim=c(200,6000)) +
			labs(x='\npatristic distance between within-host clades of different individuals\n(average across genomic windows, in subst/site)', y='most likely topology classification\n(#pairings)\n', fill='pairwise\ntopological relationship')
	ggsave(file=paste0(outfile.base, 'topodist_counts_MLE2.pdf'), w=10, h=5)
	
	ggplot(rplkl2, aes(x=PD_MEAN_LOGC, fill=TYPE_TO)) +
			geom_bar(position='fill') +
			scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), expand=c(0,0)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			facet_grid(~TYPE_PAIR) +
			labs(x='\npatristic distance between within-host clades of different individuals\n(average across genomic windows, in subst/site)', y='most likely topology classification\n(to pairs of individuals)\n', fill='pairwise\ntopological relationship')
	ggsave(file=paste0(outfile.base, 'topodist_prop_MLE.pdf'), w=10, h=5)
	
	#
	#	proportion of disconnected / with intermediate pairs among FF pairs with avg distance <3.5%
	subset(rplkl2, PD_MEAN<0.035 & TYPE_PAIR=='ff')[, table(TYPE_TO)]
	#	ancestral      intermingled with intermediate           sibling      disconnected 
	#		   42                 6                 4                 2                66 
	#	70/120=58% reduction in false positive FF pairs through topology

	subset(rplkl2, PD_MEAN<0.035 & TYPE_PAIR=='hsx')[, table(TYPE_TO)]
	#	ancestral      intermingled with intermediate           sibling      disconnected 
	#		  236                30                 7                 7               134
	#	141/414=34% reduction in hsx pairs


	#
	#	OK now do in earnest on pairs that are actually selected with
	#	TYPE_PAIR_DI and TYPE_PAIR_TODI
	#	at different distances
	prior.neff		<- 3
	prior.keff		<- 2
	confidence.cut	<- 0.66
	indir			<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun23'	
	infiles			<- data.table(F=list.files(indir, pattern='trmStatsPerWindow.rda$', full.names=TRUE))
	for(close in seq(0.01,0.04,0.005))
	{		
		ans.rplkl	<- vector('list', nrow(infiles))
		ans.rpw		<- vector('list', nrow(infiles))
		j			<- 1
		for(i in seq_len(nrow(infiles)))
		{			
			infile		<- infiles[i,F]
			cat('file',infile,'\nclose',close,'\n')
			#infile	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun/ptyr126_trmStatsPerWindow.rda'
			load(infile)	
			pty.run		<- as.numeric(gsub('^ptyr([0-9]+)_.*','\\1',basename(infile)))
			#	add male, female, IN_COUPLE
			setnames(tt, 	c('PAT.1','PAT.2','PAT.1_TIPS','PAT.2_TIPS','PAT.1_READS','PAT.2_READS','PATHS.12','PATHS.21'),
					c('ID1','ID2','ID1_L','ID2_L','ID1_R','ID2_R','PATHS_12','PATHS_21'))
			tmp			<- unique(subset(rd, select=c('RID','SEX')))
			setnames(tmp, colnames(tmp), paste0('ID1_',colnames(tmp)))
			setnames(tmp, c('ID1_RID'), c('ID1'))	
			tt			<- merge(tt, tmp, by=c('ID1'))	
			setnames(tmp, colnames(tmp), gsub('ID1','ID2',colnames(tmp)))		
			tt			<- merge(tt, tmp, by=c('ID2'))
			tt[, ID1_COUPLE:= 0L]
			set(tt, tt[, which(ID1%in%rp$MALE_RID)], 'ID1_COUPLE', 1L)
			set(tt, tt[, which(ID1%in%rp$FEMALE_RID)], 'ID1_COUPLE', 1L)	
			tt[, ID2_COUPLE:= 0L]
			set(tt, tt[, which(ID2%in%rp$MALE_RID)], 'ID2_COUPLE', 1L)
			set(tt, tt[, which(ID2%in%rp$FEMALE_RID)], 'ID2_COUPLE', 1L)		
			tt			<- subset(tt, (ID1_COUPLE==1L | ID2_COUPLE==1L) & !(ID1_SEX=='M' & ID2_SEX=='M'))
			setnames(tt, 	c('ID1','ID2','ID1_L','ID2_L','ID1_R','ID2_R','PATHS_12','PATHS_21'),
							c('PAT.1','PAT.2','PAT.1_TIPS','PAT.2_TIPS','PAT.1_READS','PAT.2_READS','PATHS.12','PATHS.21'))					
			if(nrow(tt))
			{
				#	redo relationships + likelihood at new distance threshold		
				tmp			<- phsc.get.pairwise.relationships.likelihoods(tt, trmw.min.reads=20, trmw.min.tips=1, close, trmw.distant.brl=0.08, prior.keff=prior.keff, prior.neff=prior.neff, prior.calibrated.prob=confidence.cut, relationship.types=c('TYPE_PAIR_DI','TYPE_PAIR_TODI','TYPE_PAIR_TODI2','TYPE_DIR_TODI3','TYPE_DIRSCORE_TODI3'), verbose=FALSE)
				rpw			<- tmp$dwin	
				rplkl		<- tmp$rplkl				
				#	mop up
				rplkl		<- merge(rplkl, unique(subset(rpw, select=c(ID1,ID2,ID1_SEX,ID2_SEX,ID1_COUPLE,ID2_COUPLE))), by=c('ID1','ID2'))
				rpw			<- melt(rpw, variable.name='GROUP', value.name='TYPE', measure.vars=c('TYPE_PAIR_DI','TYPE_PAIR_TODI','TYPE_PAIR_TODI2','TYPE_DIR_TODI3','TYPE_DIRSCORE_TODI3'))
				set(rpw, NULL, 'ID_R_MAX', rpw[, pmax(ID1_R,ID2_R)])
				set(rpw, NULL, 'ID_R_MIN', rpw[, pmin(ID1_R,ID2_R)])	
				rpw[, PTY_RUN:=pty.run]
				rplkl[, PTY_RUN:=pty.run]
				rpw[, CLOSE_BRL:=close]
				rplkl[, CLOSE_BRL:=close]		
				#	save
				ans.rplkl[[j]]	<- copy(rplkl)
				ans.rpw[[j]]	<- copy(rpw)
				j			<- j+1
			}			
		}
		length(ans.rplkl)	<- j
		length(ans.rpw)		<- j
		rplkl	<- do.call('rbind',ans.rplkl)
		rpw		<- do.call('rbind',ans.rpw)
		save(rpw, rplkl, file=paste0(outfile.base,'close',close*1e3,'.rda'))
		gc()
	}
	#
	#	select likely pairs 
	confidence.cut	<- 0.5
	confidence.cut	<- 0.66
	indir		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run'	
	infiles		<- data.table(F=list.files(indir, pattern='todi_ff_170610_.*close[0-9]+.rda$', full.names=TRUE))
	rtp			<- infiles[, {
				#F			<- '/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_ff_170610_cl3_prior23_close10.rda'
				load(F)
				#	select pairings from those runs with most data
				pairings		<- unique(subset(rplkl, select=c(ID1,ID2,PTY_RUN,NEFF)))[, list(PTY_RUN=PTY_RUN[which.max(NEFF)]), by=c('ID2','ID1')]
				rplkl			<- merge(pairings, rplkl, by=c('ID2','ID1','PTY_RUN'))
				#	likely transmission pairs, using distance
				rtp			<- subset(rplkl, GROUP=='TYPE_PAIR_DI')[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2','PTY_RUN')]
				rtp			<- merge(subset(rtp, TYPE_MLE=='close', c('ID1','ID2')), subset(rplkl, GROUP=='TYPE_PAIR_DI' & TYPE=='close'), by=c('ID1','ID2'), all.x=1)
				#rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
				rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-2)]						
				rtp			<- subset(rtp, POSTERIOR_SCORE>confidence.cut)
				tmp			<- copy(rtp)
				#	likely transmission pairs, using distance+topoloy
				rtp			<- subset(rplkl, GROUP=='TYPE_PAIR_TODI2')[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2','PTY_RUN')]
				rtp			<- merge(subset(rtp, TYPE_MLE=='likely pair', c('ID1','ID2')), subset(rplkl, GROUP=='TYPE_PAIR_TODI2' & TYPE=='likely pair'), by=c('ID1','ID2'), all.x=1)
				rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-2)]
				#rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]			
				rtp			<- subset(rtp, POSTERIOR_SCORE>confidence.cut)
				rtp			<- rbind(rtp, tmp)	
				rtp
			}, by='F']
	
	#	count FF, total positive by TYPE_PAIR_DI, TYPE_PAIR_TODI2, TYPE_PAIR
	rtp[, TYPE_PAIR:= tolower(paste0(ID1_SEX, ID2_SEX))]
	set(rtp, rtp[, which(TYPE_PAIR%in%c('fm','mf'))], 'TYPE_PAIR', 'hsx')
	#	exclude PTY_RUN 1 for now 
	#rtp		<- subset(rtp, PTY_RUN!=1)
	ffs		<- rtp[, list(N=length(ID1)), by=c('GROUP','TYPE_PAIR','CLOSE_BRL')]
	ffs		<- rbind(ffs, ffs[, list(N=sum(N), TYPE_PAIR='ff_hsx'), by=c('GROUP','CLOSE_BRL')])
	#	how many individuals in total did we run stage3 phyloscanner on:
	load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_stagethree.rda')
	tmp		<- unique(subset(pty.runs, select=RID))
	load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_stagetwo.rda')
	tmp		<- unique(rbind(subset(pty.runs, PTY_RUN!=1, RID),tmp))
	tmp		<- merge(tmp, unique(subset(rd, select=c(RID,SEX))),by='RID')
	ffs[, N_FEMALE:= nrow(subset(tmp, SEX=='F'))]
	ffs[, N_MALE:= nrow(subset(tmp, SEX=='M'))]
	ffs[, N_COUPLES:=nrow(unique(subset(rp, select=c(MALE_RID,FEMALE_RID))))]
	ffs		<- dcast.data.table(ffs, GROUP+CLOSE_BRL+N_FEMALE+N_COUPLES~TYPE_PAIR, value.var='N')	
	ffs[, FPR:= ff/((N_FEMALE-1)*N_COUPLES/2) ]
	ffs[, FDR:= ff/hsx]
	tmp		<- suppressWarnings(melt(ffs, id.vars=c('GROUP','CLOSE_BRL'), measure.vars=c('ff','FPR','FDR')))
	set(tmp, NULL, 'GROUP', tmp[, factor(GROUP, levels=c('TYPE_PAIR_DI','TYPE_PAIR_TODI2'), labels=c('only distance between within-host trees','distance and topology between within-host trees'))])
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('ff','FPR','FDR'), labels=c('FF pairs\nfalsely identified as\nlikely transmission pair','false positive rate\n(estimated from FF pairs)\n','false discovery rate\n(estimated from FF pairs)\n'))])
	ggplot(tmp, aes(x=CLOSE_BRL, y=value, fill=GROUP)) +
			geom_bar(stat='identity', position='dodge') +
			scale_x_continuous(labels=scales::percent) +
			scale_fill_manual(values=c('only distance between within-host trees'=rev(brewer.pal(11, 'PuOr'))[3],'distance and topology between within-host trees'=brewer.pal(11, 'PuOr')[2])) +
			theme_bw() + theme(legend.position='bottom') +
			facet_wrap(~variable, scales='free_y') +
			labs( x='\ndistance between within-host trees of different individuals\n(threshold for calling a likely transmission pair in subst/site)',
				  y='',
				  fill='information used')
  	ggsave(file=paste0(outfile.base, 'topodist_FF_fpr_fdr.pdf'), w=7, h=4)
	
}



RakaiFull.analyze.ffpairs.todi.170811<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	neff.cut<- 3
	cuts	<- data.table(	ID=   seq(1,7),
							PHSC= c(0.01,0.015,0.02,0.025,0.03,0.035,0.04),
							RWGD= c(0.0185, 0.026, 0.032, 0.0375, 0.041, 0.045, 0.048),
							FT=   c(0.037, 0.053, 0.066, 0.077, 0.088, 0.098, 0.107))
	cuts	<- melt(cuts, measure.vars=c('PHSC','RWGD','FT'))
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_ff_170811_cl3_prior23_min10_"
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_ff_170811_cl3_prior23_min30_"
	
	#	load dc
	load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170704_assemblystatus.rda')	
	# 	load rd rp
	#	load couples to search for in phyloscanner output
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_info.rda")
	#
	#	load demographic info on all individuals
	tmp		<- RakaiCirc.epi.get.info.170208()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#rn		<- tmp$rn
	ra		<- tmp$ra
	#set(rn, NULL, 'RID', rn[, as.character(RID)])
	#rn		<- merge(rn, subset(rd, select=c(RID, FIRSTPOSVIS, FIRSTPOSDATE)), by='RID', all.x=1)
	#tmp		<- rn[, which(is.na(FIRSTPOSDATE) & !is.na(RECENTVLDATE) & TIMESINCEVL==0)]	#this is dodgy
	#set(rn, tmp, 'FIRSTPOSDATE', rn[tmp, RECENTVLDATE])		
	#rd		<- rbind(rd, rn, use.names=TRUE, fill=TRUE)	#do not consider individuals in the neuro study that are not part of RCCS
	set(rd, NULL, c('PID','SID'), NULL)
	set(rd, NULL, 'SEX', rd[, as.character(SEX)])
	set(rd, NULL, 'RECENTVL', rd[, as.numeric(gsub('< 150','1',gsub('> ','',gsub('BD','',gsub(',','',as.character(RECENTVL))))))])
	set(rd, NULL, 'CAUSE_OF_DEATH', rd[, as.character(CAUSE_OF_DEATH)])
	#	fixup rd: 
	#	remove HIV reverters without sequence
	rd		<- subset(rd, !RID%in%c("C117824","C119303","E118889","K067249"))
	#	fixup complex serology
	set(rd, rd[, which(RID=='B106184')], 'FIRSTPOSDATE', rd[which(RID=='B106184'),DATE])
	set(rd, rd[, which(RID=='B106184')], c('LASTNEGVIS','LASTNEGDATE'), NA_real_)
	set(rd, rd[, which(RID=='B106184')], c('HIVPREV'), 1)
	set(rd, rd[, which(RID=='A008742')], 'FIRSTPOSDATE', rd[which(RID=='A008742'),DATE])
	set(rd, rd[, which(RID=='A008742')], c('HIVPREV'), 1)
	#	fixup rd: 
	#	missing first pos date
	rd		<- subset(rd, RID!='A038432')	#has missing firstposdate and not in PANGEA anyway
	rd		<- subset(rd, RID!='H013226')	#has missing firstposdate and not in PANGEA anyway
	rd		<- subset(rd, RID!='K008173')	#has missing firstposdate and not in PANGEA anyway
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
						if(any(!is.na(FIRSTPOSVIS)))
							ans[is.na(FIRSTPOSVIS) | (!is.na(FIRSTPOSVIS) & FIRSTPOSVIS!=min(FIRSTPOSVIS, na.rm=TRUE))]	<- TRUE							
						list(DUMMY=DUMMY, DELETE=ans)		
					}, by=c('RID')], by=c('RID','DUMMY'))
	tmp		<- subset(tmp, !DELETE)
	set(tmp, NULL,c('DUMMY','DELETE'), NULL)
	set(rd, NULL, c('BIRTHDATE','LASTNEGDATE','FIRSTPOSVIS','FIRSTPOSDATE','ARVSTARTDATE','EST_DATEDIED'), NULL)
	rd		<- merge(rd, tmp, by='RID')	
	tmp		<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
	stopifnot(!nrow(merge(subset(tmp[, length(BIRTHDATE), by='RID'], V1>1), tmp, by='RID')))	
	
	
	#	load phyloscanner results for all windows
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min10_allwindows.rda"
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_allwindows.rda"
	load(infile)	
	#	we only have raw genetic distances, 
	#	remember the FastTree patristic distances were pretty bad and ExaML never converged 
	infile	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_rawgendist.rda"
	load(infile)		
	gd		<- subset(gd, grepl('^PG',TAXA) & grepl('^PG',TAXA2))
	set(gd, NULL, 'TAXA', gd[, gsub('_WTSI.*','',TAXA)])
	set(gd, NULL, 'TAXA2', gd[, gsub('_WTSI.*','',TAXA2)])
	setnames(gd, c('TAXA'), c('PIDF'))	
	tmp		<- unique(subset(dc, !is.na(PIDF), select=c(RID, PIDF)))	
	gd		<- merge(gd, tmp, by='PIDF', all.x=TRUE)
	setnames(gd, c('PIDF','TAXA2','RID'), c('PIDF_2','PIDF','RID_2'))
	gd		<- merge(gd, tmp, by='PIDF', all.x=TRUE)
	#	there NA RIDs because we may have no overlap / deselected MRC sequences
	gd		<- subset(gd, !is.na(RID) & !is.na(RID_2) & RID_2!=RID)
	gd		<- gd[, list(CONS_GDRW=median(CONS_GDRW, na.rm=TRUE)),	by=c('RID','RID_2')]
	#	make it directly comparable to the phyloscanner study 
	#	select pairs for whom we have sufficient NGS data from both
	#	this is not quadratic because we only look at close individuals in stage2
	tmp		<- unique(subset(rplkl, GROUP=='TYPE_PAIR_DI2' & NEFF>=neff.cut, c(ID1, ID2, PTY_RUN, NEFF)))
	setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
	tmp		<- rbind( unique(subset(rplkl, GROUP=='TYPE_PAIR_DI2' & NEFF>=neff.cut, c(ID1, ID2, PTY_RUN, NEFF))), tmp)
	setnames(tmp, c('ID1','ID2'), c('RID', 'RID_2'))
	tmp[, PHSC:=1L]
	#	there are duplicates in the above, each pair can be in multiple runs
	#load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170704_assemblystatus.rda')
	#tmp2	<- subset(dc, !is.na(PIDF), select=c(RID, PIDF))
	#merge(tmp, tmp2, by="RID")	
	gd		<- merge(tmp, gd, by=c('RID','RID_2'), all.x=TRUE)
	#	add gender
	#	we are losing a few more because we killed some RIDs, most notably those in the neuro study
	tmp		<- unique(subset(rd, select=c(RID, SEX)))	
	gd		<- merge(gd, tmp, by='RID')
	setnames(tmp, c('RID','SEX'), c('RID_2','SEX_2'))
	gd		<- merge(gd, tmp, by=c('RID_2'))			
	#	we will exclude male,female in couple that are not themselves a couple
	#	as this is not part of the calculations
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_info.rda")
	tmp	<- unique(subset(rp, select=c(MALE_RID,FEMALE_RID)))
	setnames(tmp, c('MALE_RID','FEMALE_RID'), c('FEMALE_RID','MALE_RID'))
	tmp	<- rbind(subset(rp, select=c(MALE_RID,FEMALE_RID)), tmp)
	setnames(tmp, c('MALE_RID','FEMALE_RID'), c('RID','RID_2'))
	tmp[, COUPLE:=1L]
	gd	<- merge(gd, tmp, by=c('RID','RID_2'),all.x=TRUE)
	set(gd, gd[, which(is.na(COUPLE))], 'COUPLE', 0L)	
	#	make RID male or RID < RID_2 for FF pairs
	gd[, DUMMY:=NA_character_]
	tmp		<- gd[, which(COUPLE==1 & SEX=='F')]	
	set(gd, tmp, 'DUMMY', gd[tmp, RID])
	set(gd, tmp, 'RID', gd[tmp, RID_2])
	set(gd, tmp, 'RID_2', gd[tmp, DUMMY])
	set(gd, NULL, 'DUMMY', NA_character_)	
	set(gd, tmp, 'SEX', 'M')
	set(gd, tmp, 'SEX_2', 'F')
	tmp		<- gd[, which(COUPLE==0 & RID>RID_2)]	
	set(gd, tmp, 'DUMMY', gd[tmp, RID])
	set(gd, tmp, 'RID', gd[tmp, RID_2])
	set(gd, tmp, 'RID_2', gd[tmp, DUMMY])
	set(gd, NULL, 'DUMMY', NULL)
	#	select run in which pair has largest NEFF
	gd[, DUMMY:= seq_len(nrow(gd))]
	tmp		<- gd[, {
				z<- which(NEFF==max(NEFF))
				if(any(!is.na(CONS_GDRW[z])))
					z<- z[!is.na(CONS_GDRW[z])]
				list(DUMMY=DUMMY[z[1]])
			}, by=c('RID','RID_2')]
	gd		<- merge(gd, tmp, by=c('RID','RID_2','DUMMY'))
	#	we will exclude FF pairs where at least one female is not in couple
	tmp	<- unique(c(rp[, MALE_RID], rp[, FEMALE_RID]))
	gd[, RID_IN_COUPLE:=0L]
	gd[, RID_2_IN_COUPLE:=0L]
	set(gd, gd[, which(RID%in%tmp)],'RID_IN_COUPLE',1L)
	set(gd, gd[, which(RID_2%in%tmp)],'RID_2_IN_COUPLE',1L)
	
	#	select either FF pairs or couples
	gds		<- subset(gd, COUPLE==1 | (SEX=='F' & SEX_2=='F' & (RID_IN_COUPLE==1 | RID_2_IN_COUPLE==1)))
	
	#
	#	count number of FF pairs within raw genetic distance
	#	
	gdsi	<- cuts[, list(	FF_PRS_CLOSE= nrow(subset(gds, !is.na(CONS_GDRW) & SEX=='F' & SEX_2=='F' & CONS_GDRW<=value)),
							FF_PRS_NOTCLOSE= nrow(subset(gds, !is.na(CONS_GDRW) &SEX=='F' & SEX_2=='F' & CONS_GDRW>value)),
							CPLS_CLOSE= nrow(subset(gds, !is.na(CONS_GDRW) &COUPLE==1 & CONS_GDRW<=value)),
							CPLS_NOTCLOSE= nrow(subset(gds, !is.na(CONS_GDRW) & COUPLE==1 & CONS_GDRW>value)),
							FF_DENOM= nrow(subset(gds, !is.na(CONS_GDRW) & SEX=='F' & SEX_2=='F')),
							CLOSE= nrow(subset(gds, !is.na(CONS_GDRW) & CONS_GDRW<=value))
							), by=c('ID','variable','value')]
	#	this looks pretty good!!

	#
	#	process phyloscanner output at different distance cutoffs
	#
	if(0)
	{
		prior.neff		<- 3
		prior.keff		<- 2
		confidence.cut	<- 0.66
		trmw.min.reads	<- 30
		indir			<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min10'
		indir			<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30'
		infiles			<- data.table(F=list.files(indir, pattern='trmStatsPerWindow.rda$', full.names=TRUE))
		for(close in seq(0.01,0.04,0.005))
		{		
			ans.rplkl	<- vector('list', nrow(infiles))
			ans.rpw		<- vector('list', nrow(infiles))
			j			<- 1
			for(i in seq_len(nrow(infiles)))
			{			
				infile		<- infiles[i,F]
				cat('\nfile',infile,'\nclose',close,'\n')
				#infile	<- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun/ptyr126_trmStatsPerWindow.rda'
				load(infile)	
				pty.run		<- as.numeric(gsub('^ptyr([0-9]+)_.*','\\1',basename(infile)))
				#	add male, female, IN_COUPLE
				setnames(tt, 	c('PAT.1','PAT.2','PAT.1_TIPS','PAT.2_TIPS','PAT.1_READS','PAT.2_READS','PATHS.12','PATHS.21'),
						c('ID1','ID2','ID1_L','ID2_L','ID1_R','ID2_R','PATHS_12','PATHS_21'))
				tmp			<- unique(subset(rd, select=c('RID','SEX')))
				setnames(tmp, colnames(tmp), paste0('ID1_',colnames(tmp)))
				setnames(tmp, c('ID1_RID'), c('ID1'))	
				tt			<- merge(tt, tmp, by=c('ID1'))	
				setnames(tmp, colnames(tmp), gsub('ID1','ID2',colnames(tmp)))		
				tt			<- merge(tt, tmp, by=c('ID2'))
				#	define couple
				tmp	<- unique(subset(rp, select=c(MALE_RID,FEMALE_RID)))
				setnames(tmp, c('MALE_RID','FEMALE_RID'), c('FEMALE_RID','MALE_RID'))
				tmp	<- rbind(unique(subset(rp, select=c(MALE_RID,FEMALE_RID))), tmp)
				setnames(tmp, c('MALE_RID','FEMALE_RID'), c('ID1','ID2'))
				tmp[, COUPLE:=1L]
				tt	<- merge(tt, tmp, by=c('ID1','ID2'),all.x=TRUE)
				set(tt, tt[, which(is.na(COUPLE))], 'COUPLE', 0L)					
				#	make RID male or RID < RID_2 for FF pairs
				tmp	<- subset(tt, ID1_SEX=='F' & ID2_SEX=='F' & ID1>ID2)
				setnames(tmp, c('ID1','ID2','PATHS_12','PATHS_21','ID1_L','ID1_R','ID2_L','ID2_R','ID1_SEX','ID2_SEX'), c('ID2','ID1','PATHS_21','PATHS_12','ID2_L','ID2_R','ID1_L','ID1_R','ID2_SEX','ID1_SEX'))
				set(tmp, NULL, 'TYPE', tmp[,gsub('xxx','21',gsub('21','12',gsub('12','xxx',TYPE)))])
				tmp2<- rbind(subset(tt, ID1_SEX=='F' & ID2_SEX=='F' & ID1<ID2), tmp)
				tmp	<- subset(tt, ID1_SEX=='F' & ID2_SEX=='M')
				setnames(tmp, c('ID1','ID2','PATHS_12','PATHS_21','ID1_L','ID1_R','ID2_L','ID2_R','ID1_SEX','ID2_SEX'), c('ID2','ID1','PATHS_21','PATHS_12','ID2_L','ID2_R','ID1_L','ID1_R','ID2_SEX','ID1_SEX'))
				set(tmp, NULL, 'TYPE', tmp[,gsub('xxx','21',gsub('21','12',gsub('12','xxx',TYPE)))])
				tmp	<- rbind(subset(tt, ID1_SEX=='M' & ID2_SEX=='F'),tmp)
				tmp	<- rbind(tmp, tmp2)
				#	define if individual part of couple
				tt[, ID1_COUPLE:= 0L]
				set(tt, tt[, which(ID1%in%rp$MALE_RID)], 'ID1_COUPLE', 1L)
				set(tt, tt[, which(ID1%in%rp$FEMALE_RID)], 'ID1_COUPLE', 1L)	
				tt[, ID2_COUPLE:= 0L]
				set(tt, tt[, which(ID2%in%rp$MALE_RID)], 'ID2_COUPLE', 1L)
				set(tt, tt[, which(ID2%in%rp$FEMALE_RID)], 'ID2_COUPLE', 1L)		
				#	select				
				tt			<- subset(tt, COUPLE==1 | (ID1_SEX=='F' & ID2_SEX=='F' & (ID1_COUPLE==1 | ID2_COUPLE==1)))				
				setnames(tt, 	c('ID1','ID2','ID1_L','ID2_L','ID1_R','ID2_R','PATHS_12','PATHS_21'),
								c('PAT.1','PAT.2','PAT.1_TIPS','PAT.2_TIPS','PAT.1_READS','PAT.2_READS','PATHS.12','PATHS.21'))					
				if(nrow(tt))
				{
					#	redo relationships + likelihood at new distance threshold		
					tmp			<- phsc.get.pairwise.relationships.likelihoods(tt, trmw.min.reads=trmw.min.reads, trmw.min.tips=1, close, trmw.distant.brl=0.08, prior.keff=prior.keff, prior.neff=prior.neff, prior.calibrated.prob=confidence.cut, relationship.types=c('TYPE_PAIR_DI2','TYPE_PAIR_TODI2','TYPE_DIR_TODI2','TYPE_NETWORK_SCORES'), verbose=FALSE)
					rpw			<- tmp$dwin	
					rplkl		<- tmp$rplkl				
					#	mop up
					rplkl		<- merge(rplkl, unique(subset(rpw, select=c(ID1,ID2,ID1_SEX,ID2_SEX,ID1_COUPLE,ID2_COUPLE))), by=c('ID1','ID2'))
					rpw			<- melt(rpw, variable.name='GROUP', value.name='TYPE', measure.vars=c('TYPE_PAIR_DI2','TYPE_PAIR_TODI2','TYPE_DIR_TODI2','TYPE_NETWORK_SCORES'))
					set(rpw, NULL, 'ID_R_MAX', rpw[, pmax(ID1_R,ID2_R)])
					set(rpw, NULL, 'ID_R_MIN', rpw[, pmin(ID1_R,ID2_R)])	
					rpw[, PTY_RUN:=pty.run]
					rplkl[, PTY_RUN:=pty.run]
					rpw[, CLOSE_BRL:=close]
					rplkl[, CLOSE_BRL:=close]		
					#	save
					ans.rplkl[[j]]	<- copy(rplkl)
					ans.rpw[[j]]	<- copy(rpw)
					j			<- j+1
				}			
			}
			length(ans.rplkl)	<- j
			length(ans.rpw)		<- j
			rplkl	<- do.call('rbind',ans.rplkl)
			rpw		<- do.call('rbind',ans.rpw)
			save(rpw, rplkl, file=paste0(outfile.base,'close',close*1e3,'.rda'))
			gc()
		}
	}
	#
	#	select likely pairs using standard posterior mode rule
	if(0)
	{
		confidence.cut	<- 0.66
		neff.cut		<- 3
		indir		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run'	
		infiles		<- data.table(F=list.files(indir, pattern='todi_ff_170811_cl3_prior23_min30_close[0-9]+.rda$', full.names=TRUE))
		rtp			<- infiles[, {
					#F			<- '/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_ff_170811_cl3_prior23_min10_close10.rda'
					load(F)
					#	select pairings from those runs with most data
					pairings		<- unique(subset(rplkl, GROUP=='TYPE_PAIR_DI2' & NEFF>=neff.cut, select=c(ID1,ID2,PTY_RUN,NEFF)))[, list(PTY_RUN=PTY_RUN[which.max(NEFF)]), by=c('ID2','ID1')]
					rplkl			<- merge(pairings, rplkl, by=c('ID2','ID1','PTY_RUN'))
					#	likely transmission pairs, using distance
					rtp			<- subset(rplkl, GROUP=='TYPE_PAIR_DI2')[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2','PTY_RUN')]
					rtp			<- merge(subset(rtp, TYPE_MLE=='close', c('ID1','ID2')), subset(rplkl, GROUP=='TYPE_PAIR_DI2' & TYPE=='close'), by=c('ID1','ID2'), all.x=1)				
					rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-2)]						
					rtp			<- subset(rtp, POSTERIOR_SCORE>confidence.cut)
					tmp			<- copy(rtp)
					#	likely transmission pairs, using distance+topoloy
					rtp			<- subset(rplkl, GROUP=='TYPE_PAIR_TODI2')[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2','PTY_RUN')]
					rtp			<- merge(subset(rtp, TYPE_MLE=='linked', c('ID1','ID2')), subset(rplkl, GROUP=='TYPE_PAIR_TODI2' & TYPE=='linked'), by=c('ID1','ID2'), all.x=1)
					rtp[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-2)]							
					rtp			<- subset(rtp, POSTERIOR_SCORE>confidence.cut)
					rtp			<- rbind(rtp, tmp)	
					rtp
				}, by='F']
		save(rtp, file=paste0(outfile.base,'close_all.rda'))
	}
	#	select likely pairs using standard posterior mode rule
	if(0)
	{
		confidence.cut	<- 0.66
		neff.cut		<- 3
		indir		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run'	
		infiles		<- data.table(F=list.files(indir, pattern='todi_ff_170811_cl3_prior23_min30_close[0-9]+.rda$', full.names=TRUE))
		rtp			<- infiles[, {
					#F			<- '/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_ff_170811_cl3_prior23_min30_close15.rda'
					load(F)
					#	select pairings from those runs with most data
					pairings		<- unique(subset(rplkl, GROUP=='TYPE_PAIR_DI2' & NEFF>=neff.cut, select=c(ID1,ID2,PTY_RUN,NEFF)))[, list(PTY_RUN=PTY_RUN[which.max(NEFF)]), by=c('ID2','ID1')]
					rplkl			<- merge(pairings, rplkl, by=c('ID2','ID1','PTY_RUN'))
					#	likely transmission pairs, using distance
					rtp			<- subset(rplkl, GROUP=='TYPE_PAIR_DI2')[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2','PTY_RUN')]
					rtp			<- merge(subset(rtp, TYPE_MLE=='close', c('ID1','ID2')), subset(rplkl, GROUP=='TYPE_PAIR_DI2' & TYPE=='close'), by=c('ID1','ID2'), all.x=1)				
					rtp[, POSTERIOR_SCORE:= KEFF/NEFF]						
					rtp			<- subset(rtp, POSTERIOR_SCORE>confidence.cut)
					tmp			<- copy(rtp)
					#	likely transmission pairs, using distance+topoloy
					rtp			<- subset(rplkl, GROUP=='TYPE_PAIR_TODI2')[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('ID1','ID2','PTY_RUN')]
					rtp			<- merge(subset(rtp, TYPE_MLE=='linked', c('ID1','ID2')), subset(rplkl, GROUP=='TYPE_PAIR_TODI2' & TYPE=='linked'), by=c('ID1','ID2'), all.x=1)
					rtp[, POSTERIOR_SCORE:= KEFF/NEFF]							
					rtp			<- subset(rtp, POSTERIOR_SCORE>confidence.cut)
					rtp			<- rbind(rtp, tmp)	
					rtp
				}, by='F']
		save(rtp, file=paste0(outfile.base,'close_all_knratio.rda'))
	}
	load(paste0(outfile.base,'close_all.rda'))
	#load(paste0(outfile.base,'close_all_knratio.rda'))
	#	count FF, total positive by TYPE_PAIR_DI, TYPE_PAIR_TODI2, TYPE_PAIR
	rtp[, TYPE_PAIR:= tolower(paste0(ID1_SEX, ID2_SEX))]
	set(rtp, rtp[, which(TYPE_PAIR%in%c('fm','mf'))], 'TYPE_PAIR', 'hsx')
	rtp		<- subset(rtp, select=c(GROUP,TYPE_PAIR,CLOSE_BRL, ID1, ID2, ID1_SEX, ID2_SEX, NEFF))	
	#	make RID male or RID < RID_2 for FF pairs
	rtp[, DUMMY:=NA_character_]
	tmp		<- rtp[, which(ID1_SEX=='F' & ID2_SEX=='M')]	
	set(rtp, tmp, 'DUMMY', rtp[tmp, ID1])
	set(rtp, tmp, 'ID1', rtp[tmp, ID2])
	set(rtp, tmp, 'ID2', rtp[tmp, DUMMY])
	set(rtp, NULL, 'DUMMY', NA_character_)	
	set(rtp, tmp, 'ID1_SEX', 'M')
	set(rtp, tmp, 'ID2_SEX', 'F')
	tmp		<- rtp[, which(ID1_SEX=='F' & ID2_SEX=='F' & ID1>ID2)]	
	set(rtp, tmp, 'DUMMY', rtp[tmp, ID1])
	set(rtp, tmp, 'ID1', rtp[tmp, ID2])
	set(rtp, tmp, 'ID2', rtp[tmp, DUMMY])
	set(rtp, NULL, 'DUMMY', NULL)	
	rtp		<- unique(rtp, by=c('GROUP','CLOSE_BRL','ID1','ID2'))		
	tmp		<- subset(gds, select=c(RID,RID_2,CONS_GDRW))
	setnames(tmp, c('RID','RID_2'),c('ID1','ID2'))	
	rtp		<- merge(rtp, tmp, by=c('ID1','ID2'),all.x=TRUE)
	
	
	#	scale up
	ffs		<- rtp[, list(N=length(ID1)), by=c('GROUP','TYPE_PAIR','CLOSE_BRL')]
	ffs		<- rbind(ffs, ffs[, list(N=sum(N), TYPE_PAIR='ff_hsx'), by=c('GROUP','CLOSE_BRL')])
	ffs		<- dcast.data.table(ffs, GROUP+CLOSE_BRL~TYPE_PAIR, value.var='N')
	ffs[, FF_DENOM:= nrow(subset(gds, COUPLE==0))]	
	setnames(ffs, c('ff','hsx','ff_hsx'), c('FF_PRS_CLOSE', 'CPLS_CLOSE', 'CLOSE'))
	tmp		<- subset(gdsi, variable=='RWGD', select=c(ID, FF_PRS_CLOSE, CPLS_CLOSE, FF_DENOM, CLOSE))
	tmp		<- merge(subset(gdsi, variable=='PHSC', select=c(ID, variable, value)),tmp,by='ID')
	setnames(tmp, c('variable','value'), c('GROUP', 'CLOSE_BRL'))
	tmp[, FF_DENOM_SCALE:= FF_DENOM/ffs$FF_DENOM[1]]
	tmp		<- melt(tmp, id.vars=c('ID','GROUP','CLOSE_BRL','FF_DENOM_SCALE'))
	tmp[, N:=round(value/FF_DENOM_SCALE)]
	tmp		<- dcast.data.table(tmp,GROUP+CLOSE_BRL~variable, value.var='N')
	ffs		<- rbind(ffs, tmp)
	ffs[, FPR:= FF_PRS_CLOSE/FF_DENOM]
	ffs[, FDR:= (FF_PRS_CLOSE+0)/CLOSE]
	tmp		<- suppressWarnings(melt(ffs, id.vars=c('GROUP','CLOSE_BRL'), measure.vars=c('FF_PRS_CLOSE','FPR','FDR')))
	set(tmp, NULL, 'GROUP', tmp[, factor(GROUP, levels=c('PHSC','TYPE_PAIR_DI2','TYPE_PAIR_TODI2'), labels=c('consensus sequences','NGS, only subtree distance','NGS, subtree distance and topology'))])
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('FF_PRS_CLOSE','FPR','FDR'), labels=c('FF pairs\nfalsely identified as\nlikely transmission pair','false positive rate\n(estimated from FF pairs)\n','false discovery rate\n(estimated from FF pairs)\n'))])
	ggplot(tmp, aes(x=CLOSE_BRL, y=value, fill=GROUP)) +
			geom_bar(stat='identity', position='dodge') +
			scale_x_continuous(labels=scales::percent) +
			scale_fill_manual(values=c('consensus sequences'='lightpink3','NGS, only subtree distance'=rev(brewer.pal(11, 'PuOr'))[3],'NGS, subtree distance and topology'=brewer.pal(11, 'PuOr')[2])) +
			theme_bw() + theme(legend.position='bottom') +
			facet_wrap(~variable, scales='free_y') +
			labs( x='\ndistance threshold for phylogenetic linkage\n(subst/site)',
					y='',
					fill='information used')
	ggsave(file=paste0(outfile.base, 'topodist_FF_fpr_fdr.pdf'), w=7, h=4)
	save(ffs, file=paste0(outfile.base, 'topodist_FF_fpr_fdr.rda'))
	
	#	scale down
	ffs		<- subset(rtp, !is.na(CONS_GDRW))[, list(N=length(ID1)), by=c('GROUP','TYPE_PAIR','CLOSE_BRL')]
	ffs		<- rbind(ffs, ffs[, list(N=sum(N), TYPE_PAIR='ff_hsx'), by=c('GROUP','CLOSE_BRL')])
	ffs		<- dcast.data.table(ffs, GROUP+CLOSE_BRL~TYPE_PAIR, value.var='N')
	ffs[, FF_DENOM:= nrow(subset(gds, !is.na(CONS_GDRW) & COUPLE==0))]	
	setnames(ffs, c('ff','hsx','ff_hsx'), c('FF_PRS_CLOSE', 'CPLS_CLOSE', 'CLOSE'))
	tmp		<- subset(gdsi, variable=='RWGD', select=c(ID, FF_PRS_CLOSE, CPLS_CLOSE, FF_DENOM, CLOSE))
	tmp		<- merge(subset(gdsi, variable=='PHSC', select=c(ID, variable, value)),tmp,by='ID')
	setnames(tmp, c('variable','value'), c('GROUP', 'CLOSE_BRL'))
	tmp[, FF_DENOM_SCALE:= FF_DENOM/ffs$FF_DENOM[1]]
	tmp		<- melt(tmp, id.vars=c('ID','GROUP','CLOSE_BRL','FF_DENOM_SCALE'))
	tmp[, N:=round(value/FF_DENOM_SCALE)]
	tmp		<- dcast.data.table(tmp,GROUP+CLOSE_BRL~variable, value.var='N')
	ffs		<- rbind(ffs, tmp)
	ffs[, FPR:= FF_PRS_CLOSE/FF_DENOM]
	ffs[, FDR:= (FF_PRS_CLOSE+0)/CLOSE]
	tmp		<- suppressWarnings(melt(ffs, id.vars=c('GROUP','CLOSE_BRL'), measure.vars=c('FF_PRS_CLOSE','FPR','FDR')))
	set(tmp, NULL, 'GROUP', tmp[, factor(GROUP, levels=c('PHSC','TYPE_PAIR_DI2','TYPE_PAIR_TODI2'), labels=c('consensus sequences','NGS, only subtree distance','NGS, subtree distance and topology'))])
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('FF_PRS_CLOSE','FPR','FDR'), labels=c('FF pairs\nfalsely identified as\nlikely transmission pair','false positive rate\n(estimated from FF pairs)\n','false discovery rate\n(estimated from FF pairs)\n'))])
	ggplot(tmp, aes(x=CLOSE_BRL, y=value, fill=GROUP)) +
			geom_bar(stat='identity', position='dodge') +
			scale_x_continuous(labels=scales::percent) +
			scale_fill_manual(values=c('consensus sequences'='lightpink3','NGS, only subtree distance'=rev(brewer.pal(11, 'PuOr'))[3],'NGS, subtree distance and topology'=brewer.pal(11, 'PuOr')[2])) +
			theme_bw() + theme(legend.position='bottom') +
			facet_wrap(~variable, scales='free_y') +
			labs( x='\ndistance threshold for phylogenetic linkage\n(subst/site)',
					y='',
					fill='information used')
	ggsave(file=paste0(outfile.base, 'topodist_FF_fpr_fdr_scaledown.pdf'), w=7, h=4)
	save(ffs, file=paste0(outfile.base, 'topodist_FF_fpr_fdr_scaledown.rda'))
	
	#	make word table for main text
	ans	<- copy(ffs)
	set(ans, NULL, 'FPR', ans[, round(FPR, d=3)])
	set(ans, NULL, 'FDR', ans[, paste0(round(100*FDR),'%')])
	ans	<- melt(ans, id.vars=c('GROUP','CLOSE_BRL'), measure.vars=c('FF_PRS_CLOSE','FPR','FDR'))
	ans	<- dcast.data.table(ans, variable+CLOSE_BRL~GROUP, value.var='value')	
	ans	<- subset(ans, CLOSE_BRL%in%c(0.01,0.02,0.03,0.04), select=c(variable, CLOSE_BRL, PHSC, TYPE_PAIR_DI2, TYPE_PAIR_TODI2))	
	write.csv(ans, row.names=FALSE, file=paste0(outfile.base, 'topodist_FF_fpr_fdr_scaledown.csv'))
	
	#	make supplementary figure: false pos rate by NEFF groups <5, 5-9, 10-14, >15
	require(gamlss)
	require(zoo)

	dr	<- subset(gds, SEX=='F' & SEX_2=='F', c(RID, RID_2, NEFF))
	dr[, DUMMY:=1L]
	dr	<- merge(dr, data.table(DUMMY=1, CLOSE_BRL=c(0.01,0.02,0.03,0.04)), by='DUMMY', allow.cartesian=TRUE)
	setnames(dr, c('RID','RID_2'), c('ID1','ID2'))
	tmp	<- subset(rtp, TYPE_PAIR=='ff' & GROUP=='TYPE_PAIR_TODI2' & CLOSE_BRL%in%c(0.01,0.02,0.03,0.04), c(ID1, ID2, CLOSE_BRL))
	tmp[, LINKED:=1L]
	dr	<- merge(dr, tmp, by=c('ID1','ID2','CLOSE_BRL'), all.x=TRUE)
	set(dr, dr[, which(is.na(LINKED))],'LINKED',0L)
	dr[, NEFFC:= dr[,cut(NEFF, right=FALSE, breaks=c(0,5,10,15,20,25,Inf), labels=c('<5','5-9','10-14','15-19','20-24','>=25'))]]	
	setkey(dr, CLOSE_BRL, NEFF, LINKED)
	#dr	<- dr[order(CLOSE_BRL, -NEFF, -LINKED)]
	tmp	<- dr[, list(NEFF=NEFF, LINKEDR=rollmean(LINKED, 1000, fill='extend')), by='CLOSE_BRL']
	tmp	<- tmp[, list(LINKEDR=mean(LINKEDR)), by=c('CLOSE_BRL','NEFF')]
	set(tmp, NULL, 'CLOSE_BRL', tmp[, paste0('subtree distance threshold\n',100*CLOSE_BRL,'%')])
	ggplot(tmp, aes(x=NEFF)) + 
			geom_line(aes(y=LINKEDR)) + 
			facet_grid(~CLOSE_BRL) +
			coord_cartesian(xlim=c(0,20)) +
			theme_bw() +
			labs(x='\ngenomic windows', y='false positive rate\n(rolling mean over 1000 observations)\n')
	ggsave(file=paste0(outfile.base,'FPR_by_windows.pdf'), w=10, h=5)
		
	gds[, NEFFC:= gds[,cut(NEFF, right=FALSE, breaks=c(0,5,10,15,20,25,Inf), labels=c('<5','5-9','10-14','15-19','20-24','>=25'))]]
	rtp[, NEFFC:= rtp[,cut(NEFF, right=FALSE, breaks=c(0,5,10,15,20,25,Inf), labels=c('<5','5-9','10-14','15-19','20-24','>=25'))]]
	ffs		<- subset(rtp, !is.na(CONS_GDRW))[, list(N=length(ID1)), by=c('GROUP','TYPE_PAIR','CLOSE_BRL','NEFFC')]
	ffs		<- rbind(ffs, ffs[, list(N=sum(N), TYPE_PAIR='ff_hsx'), by=c('GROUP','CLOSE_BRL','NEFFC')])
	ffs		<- dcast.data.table(ffs, GROUP+CLOSE_BRL+NEFFC~TYPE_PAIR, value.var='N')
	set(ffs, ffs[, which(is.na(ff))], 'ff', 0L)	
	tmp		<- subset(gds, !is.na(CONS_GDRW) & COUPLE==0)[, list(FF_DENOM=length(RID)), by='NEFFC']
	ffs		<- merge(ffs, tmp, by='NEFFC')
	setnames(ffs, c('ff','hsx','ff_hsx'), c('FF_PRS_CLOSE', 'CPLS_CLOSE', 'CLOSE'))
	ffs[, FPR:= FF_PRS_CLOSE/FF_DENOM]
	ffs[, FDR:= (FF_PRS_CLOSE+0)/CLOSE]
	tmp		<- ffs[, {
				z	<- as.numeric(binconf(FF_PRS_CLOSE,FF_DENOM))
				zz	<- as.numeric(binconf(FF_PRS_CLOSE,CLOSE))
				list(FPRL=z[2], FPRU=z[3], FDRL=zz[2], FDRU=zz[3])	
			}, by=c('NEFFC','GROUP','CLOSE_BRL')]
	ffs		<- merge(ffs, tmp, by=c('NEFFC','GROUP','CLOSE_BRL'))
	ffs		<- subset(ffs, CLOSE_BRL%in%c(0.01,0.02,0.03,0.04))
	set(ffs, NULL, 'CLOSE_BRL', ffs[, paste0('subtree distance threshold\n',100*CLOSE_BRL,'%')])
	ggplot(subset(ffs, GROUP=='TYPE_PAIR_TODI2'), aes(x=NEFFC)) + 
			geom_point(aes(y=FPR, size=FF_DENOM)) +
			geom_errorbar(aes(ymin=FPRL, ymax=FPRU)) +
			facet_grid(~CLOSE_BRL) +
			theme_bw() + labs(x='\ngenomic windows', y='false positive rate\n',size='denominator')
	ggsave(file=paste0(outfile.base,'FPR_by_windows_neffcgroups.pdf'), w=10, h=5)
	
	
	tmp		<- dr[, {
				tmp	<- gamlss(formula=LINKED~NEFF, family=BI())				
				tmp	<- summary(tmp, save=TRUE)				
				tmp	<- tmp$coef.table				
				list(COEF=rownames(tmp), ML=tmp[,1], PVAL=tmp[,4])				
			}, by='CLOSE_BRL']
	subset(tmp, COEF=='NEFF')
	#	increase with NEFF only significant for 0.03 and 0.04 
	#	exp(0.06)/(exp(0.06)+1)	
	
}

RakaiFull.analyze.ffpairs.intermingled.170811<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	confidence.cut			<- 0.66	
	neff.cut				<- 3
	infile.trmpairs.todi	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min10.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_ff_170811_cl3_prior23_min10_"	
	#	load second stage output
	load(infile.trmpairs.todi)
	#
	#	above NEFF cut
	tmp			<- subset(rtp.todi2, GROUP=='TYPE_PAIR_TODI2' & TYPE=='linked' & NEFF>neff.cut, select=c(ID1, ID2, PTY_RUN))
	rtp.todi2	<- merge(rtp.todi2, tmp, by=c('ID1','ID2','PTY_RUN'))
	#
	#	likely transmission pairs, using topology and distance
	#	select one patient pairing across runs: that with lowest evidence		
	rtp		<- subset(rtp.todi2, GROUP=='TYPE_PAIR_TODI2')[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)]), by=c('ID1','ID2')]	
	rtp		<- subset(rtp, ID1!=ID2)
	rtp		<- merge(rtp, subset(rtp.todi2, GROUP=='TYPE_PAIR_TODI2'), by=c('ID1','ID2','PTY_RUN'))	
	rtp		<- subset(rtp, POSTERIOR_SCORE>confidence.cut)			
	#
	#	find likely FF pairs
	#	see if predominantly close intermingled
	rff		<- subset(rtp, ID1_SEX=='F' & ID2_SEX=='F')
	rpw2	<- merge(subset(rff, select=c(ID1,ID2,PTY_RUN)), subset(rpw, GROUP=='TYPE_BASIC'), by=c('ID1','ID2','PTY_RUN'))
	rpw2[, TYPE_TO:= 'disconnected/not close']
	set(rpw2, rpw2[,which(grepl('close', TYPE) & grepl('chain', TYPE) & grepl('no intermediate', TYPE))], 'TYPE_TO', 'close ancestral')
	set(rpw2, rpw2[,which(grepl('close', TYPE) & grepl('intermingled', TYPE) & grepl('no intermediate', TYPE))], 'TYPE_TO', 'close intermingled')
	#set(rpw2, rpw2[,which(grepl('close', TYPE) & grepl('with intermediate', TYPE))], 'TYPE_TO', 'close with intermediate')
	set(rpw2, rpw2[,which(ADJACENT & grepl('other', TYPE))], 'TYPE_TO', 'close sibling')
	tmp		<- rpw2[, list(K=length(W_FROM)), by=c('ID1','ID2','PTY_RUN','TYPE_TO')]
	tmp		<- tmp[, list(TYPE_MLE=TYPE_TO[which.max(K)], K_MLE=max(K), N=sum(K)), by=c('ID1','ID2','PTY_RUN')]
	rffci	<- subset(tmp, TYPE_MLE=='close intermingled')
	#
	#	plot windows and phylogenies
	#	plot windows of identified transmission pairs
	if(0)
	{
		rps			<- subset(rffci, select=c(ID1, ID2, PTY_RUN))
		write.csv(rps, file=paste0(outfile.base,'summary.csv'))
		rps[, DUMMY:=seq_len(nrow(rps))]
		rps[, LABEL:=rps[, factor(DUMMY, levels=DUMMY, labels=paste0('f ',ID1,' f ', ID2,'\n',PTY_RUN))]]
		
		group		<- 'TYPE_BASIC'
		#group		<- 'TYPE_PAIR_TODI2'					
		rpw2		<- subset(rpw, GROUP==group)
		rplkl2		<- subset(rplkl, GROUP==group)	
		plot.file	<- paste0(outfile.base,'windows_summary_',group,'.pdf')
		setnames(rps, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rpw2, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rplkl2, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		phsc.plot.windowsummaries.for.pairs(rps, rpw2, rplkl2, plot.file, cols=NULL, group=group)				
	}	
	
	require(colorspace)
	for(ii in seq_len(nrow(rffci)))
	#for(ii in 111:252)
	{		
		indir		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min10'
		# load dfr and phs
		load( file.path(indir, paste0('ptyr',rffci[ii,PTY_RUN],'_trees.rda')) )
		# setup plotting
		ids			<- c(rffci[ii, ID1],rffci[ii, ID2])
		dfs			<- subset(dfr, select=c(W_FROM, W_TO, IDX))
		dfs[, ID1:=ids[1]]
		dfs[, ID2:=ids[2]]
		dfs			<- merge(dfs, subset(rpw, GROUP=='TYPE_RAW'), by=c('ID1','ID2','W_FROM','W_TO'), all.x=1)	
		dfs[, TITLE:= dfs[, paste('female1 ', ids[1],'\nfemale2 ',ids[2],'\nrun ', PTY_RUN, '\nwindow ', W_FROM,'-', W_TO,'\np21',PATHS_12,' p12',PATHS_21, ' ',ADJACENT,' ',CONTIGUOUS,' ',TYPE, '\n', round(PATRISTIC_DISTANCE, d=5), sep='')]]
		plot.file	<- paste0(outfile.base, '_', rffci[ii, PTY_RUN],'_F_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
		invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red4','hotpink1'), drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10, tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'))										
	}
	
}

RakaiFull.analyze.couples.todi.170522.windowssummaryplots<- function()
{		
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170522/todi_couples_170522_withmetadata.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170522/todi_couples_170522_"	
	load(infile)	
	setkey(rca, MALE_RID, FEMALE_RID)
	rca[, PAIRID:=seq_len(nrow(rca))]
	set(rca, NULL, "SELECT", rca[,gsub(' ','_',SELECT)])
	rpw[, ID_R_MAX:=pmax(MALE_R,FEMALE_R)]
	rpw[, ID_R_MIN:=pmin(MALE_R,FEMALE_R)]
	#
	#	plot windows 
	for(select in c('couple_ambiguous_if_pair_or_not_pair','couple_most_likely_a_pair_direction_not_resolved','couple_most_likely_a_pair_with_resolved_direction','couple_most_likely_not_a_pair'))
	{
		#select		<- 'couple_ambiguous_if_pair_or_not_pair'
		rps		<- subset(rca, SELECT==select, c(MALE_RID,FEMALE_RID,PTY_RUN,PAIRID))
		setnames(rps, c('MALE_RID','FEMALE_RID','PAIRID'),c('FEMALE_SANGER_ID','MALE_SANGER_ID','LABEL'))		
		for(group in c('TYPE_BASIC','TYPE_PAIR_TODI2','TYPE_DIR_TODI3'))
		{
			plot.file	<- paste0(outfile.base,'_selected_',select,'_windowssummary_',group,'.pdf')
			rpw2		<- subset(rpw, GROUP==group)
			setnames(rpw2, c('MALE_RID','FEMALE_RID'),c('FEMALE_SANGER_ID','MALE_SANGER_ID'))
			rplkl2		<- subset(rplkl, GROUP==group)
			setnames(rplkl2, c('MALE_RID','FEMALE_RID'),c('FEMALE_SANGER_ID','MALE_SANGER_ID'))
			phsc.plot.windowsummaries.for.pairs(rps, rpw2, rplkl2, plot.file, cols=NULL, group=group)
		}
	}
	
	
	
}

RakaiFull.analyze.trmpairs.todi.170811.intermingledwithintermediate<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(colorspace)
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min10_allwindows.rda'
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/intermingledintermediate_170811_cl3_prior23_min10_"
	#infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10_withmetadata.rda'
	#outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10_"
	
	load(infile)	
	#setkey(rca, MALE_RID, FEMALE_RID)
	#rca[, PAIRID:=seq_len(nrow(rca))]
	rpw2	<- subset(rpw, GROUP=='TYPE_BASIC')
	#	for every couple keep PTY_RUN with longest N
	tmp		<- unique(subset(rplkl, GROUP=='TYPE_BASIC', select=c(ID1, ID2, PTY_RUN, NEFF)))
	tmp		<- tmp[, list(PTY_RUN=PTY_RUN[which.max(NEFF)]), by=c('ID1','ID2')]
	rpw2	<- merge(rpw2, tmp, by=c('ID1','ID2', 'PTY_RUN') )
	#	define new type to find pairs that have most frequently intermingled windows
	rpw2[, TYPE_TO:= 'disconnected/not close']
	set(rpw2, rpw2[,which(grepl('close', TYPE) & grepl('chain', TYPE) & grepl('no intermediate', TYPE))], 'TYPE_TO', 'close ancestral')
	set(rpw2, rpw2[,which(grepl('close', TYPE) & grepl('intermingled', TYPE) & grepl('no intermediate', TYPE))], 'TYPE_TO', 'close intermingled')		
	set(rpw2, rpw2[,which(ADJACENT & grepl('other', TYPE))], 'TYPE_TO', 'close sibling')
	tmp		<- rpw2[, list(K=length(W_FROM)), by=c('ID1','ID2','PTY_RUN','TYPE_TO')]
	tmp		<- tmp[, list(TYPE_MLE=TYPE_TO[which.max(K)], K_MLE=max(K), N=sum(K)), by=c('ID1','ID2','PTY_RUN')]
	tmp		<- subset(tmp, TYPE_MLE=='close intermingled')
	#	of those, see how many intermediates/other windows there are, that could indicate an intermediate
	rpw3	<- merge(rpw2, subset(tmp, select=c(ID1, ID2, PTY_RUN)), by=c('ID1','ID2','PTY_RUN'))
	rpw3[, TYPE_TO:= 'not other/intermediate']
	set(rpw3, rpw3[,which(grepl('other', TYPE) | grepl('with intermediate', TYPE))], 'TYPE_TO', 'other/intermediate')
	tmp		<- rpw3[, list(K=length(W_FROM)), by=c('ID1','ID2','PTY_RUN','TYPE_TO')]
	tmp		<- tmp[, list(TYPE_TO=TYPE_TO, K=K, N=sum(K), P=K/sum(K)), by=c('ID1','ID2','PTY_RUN')]
	subset(tmp, TYPE_TO=='other/intermediate')
	
	
	
	
	#	define new type to find pairs that have most frequently windows with intermediates
	rpw2[, TYPE_TO:= 'disconnected/not close']
	set(rpw2, rpw2[,which(grepl('close', TYPE) & grepl('chain', TYPE) & grepl('no intermediate', TYPE))], 'TYPE_TO', 'close ancestral')
	set(rpw2, rpw2[,which(grepl('close', TYPE) & grepl('intermingled', TYPE) & grepl('no intermediate', TYPE))], 'TYPE_TO', 'close intermingled')
	set(rpw2, rpw2[,which(grepl('close', TYPE) & grepl('with intermediate', TYPE))], 'TYPE_TO', 'close with intermediate')	
	set(rpw2, rpw2[,which(ADJACENT & grepl('other', TYPE))], 'TYPE_TO', 'close sibling')
	tmp		<- rpw2[, list(K=length(W_FROM)), by=c('ID1','ID2','PTY_RUN','TYPE_TO')]
	tmp		<- tmp[, list(TYPE_MLE=TYPE_TO[which.max(K)], K_MLE=max(K), N=sum(K)), by=c('ID1','ID2','PTY_RUN')]
	tmp		<- subset(tmp, TYPE_MLE=='close with intermediate')
	
	#	among those pairs, see how many are most frequently intermediate
	rpw3	<- merge(rpw2, subset(tmp, select=c(ID1, ID2, PTY_RUN)), by=c('ID1','ID2','PTY_RUN'))
	rpw3[, TYPE_TO:= 'disconnected/not close']
	set(rpw3, rpw3[,which(grepl('close', TYPE) & grepl('chain', TYPE) & grepl('no intermediate', TYPE))], 'TYPE_TO', 'close ancestral')
	set(rpw3, rpw3[,which(grepl('close', TYPE) & grepl('intermingled', TYPE) & grepl('no intermediate', TYPE))], 'TYPE_TO', 'close intermingled')
	set(rpw3, rpw3[,which(grepl('close', TYPE) & grepl('intermingled', TYPE) & grepl('with intermediate', TYPE))], 'TYPE_TO', 'close with intermediate intermingled')
	set(rpw3, rpw3[,which(grepl('close', TYPE) & grepl('chain 12', TYPE) & grepl('with intermediate', TYPE))], 'TYPE_TO', 'close with intermediate chain12')
	set(rpw3, rpw3[,which(grepl('close', TYPE) & grepl('chain 21', TYPE) & grepl('with intermediate', TYPE))], 'TYPE_TO', 'close with intermediate chain21')
	set(rpw3, rpw3[,which(ADJACENT & grepl('other', TYPE))], 'TYPE_TO', 'close sibling')
	tmp		<- rpw3[, list(K=length(W_FROM)), by=c('ID1','ID2','PTY_RUN','TYPE_TO')]
	tmp		<- tmp[, list(TYPE_MLE=TYPE_TO[which.max(K)], K_MLE=max(K), N=sum(K)), by=c('ID1','ID2','PTY_RUN')]
	tmp		<- subset(tmp, TYPE_MLE=='close with intermediate intermingled')
	
	tmp[, table(TYPE_MLE)]
	
	
	#	try couples
	rpw2	<- merge(subset(rtp, select=c(MALE_RID,FEMALE_RID,PTY_RUN)), subset(rpw, GROUP=='TYPE_BASIC'), by=c('MALE_RID','FEMALE_RID','PTY_RUN'))	
	rpw2[, TYPE_TO:= 'disconnected/not close']
	set(rpw2, rpw2[,which(grepl('close', TYPE) & grepl('chain', TYPE) & grepl('no intermediate', TYPE))], 'TYPE_TO', 'close ancestral')
	set(rpw2, rpw2[,which(grepl('close', TYPE) & grepl('intermingled', TYPE) & grepl('no intermediate', TYPE))], 'TYPE_TO', 'close intermingled')
	set(rpw2, rpw2[,which(grepl('close', TYPE) & grepl('with intermediate', TYPE))], 'TYPE_TO', 'close with intermediate')	
	set(rpw2, rpw2[,which(ADJACENT & grepl('other', TYPE))], 'TYPE_TO', 'close sibling')
	tmp		<- rpw2[, list(K=length(W_FROM)), by=c('MALE_RID','FEMALE_RID','PTY_RUN','TYPE_TO')]
	tmp		<- tmp[, list(TYPE_MLE=TYPE_TO[which.max(K)], K_MLE=max(K), N=sum(K)), by=c('MALE_RID','FEMALE_RID','PTY_RUN')]
	rffci	<- subset(tmp, TYPE_MLE=='close with intermediate')
	
	
	
	tmp		<- subset(rpw, GROUP=='TYPE_BASIC')
	tmp		<- subset(tmp, grepl('intermingled',TYPE) & grepl('with intermediate',TYPE))
	rtpdm	<- unique(subset(tmp, select=c(MALE_RID,FEMALE_RID,PTY_RUN)))
	for(ii in 2:7)
	{		
		indir		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun23'
		# load dfr and phs
		load( file.path(indir, paste0('ptyr',rtpdm[ii,PTY_RUN],'_trees.rda')) )
		# setup plotting
		ids			<- c(rtpdm[ii, MALE_RID],rtpdm[ii, FEMALE_RID])
		dfs			<- subset(dfr, select=c(W_FROM, W_TO, IDX))
		dfs[, MALE_RID:=ids[1]]
		dfs[, FEMALE_RID:=ids[2]]
		dfs			<- merge(dfs, subset(rpw, GROUP=='TYPE_RAW'), by=c('MALE_RID','FEMALE_RID','W_FROM','W_TO'), all.x=1)	
		dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', PTY_RUN, '\nwindow ', W_FROM,'-', W_TO,'\npmf',PATHS_MF,' pfm',PATHS_FM, ' ',ADJACENT,' ',TYPE, '\n', round(PATRISTIC_DISTANCE, d=5), sep='')]]
		plot.file	<- paste0(outfile.base, 'checkintermingled/todi_pairs_170610_run_', rtpdm[ii, PTY_RUN],'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
		invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10, tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'))					
	}
	
	tmp		<- subset(rpw, GROUP=='TYPE_BASIC')
	subset(tmp, grepl('chain mf',TYPE) &  grepl('no intermediate',TYPE) & grepl('close',TYPE))[, table(PATHS_FM, ADJACENT)]
	tmp		<- subset(tmp, PATHS_FM==0 & PATHS_MF==0 & grepl('close',TYPE) & !ADJACENT)
	rtpdm	<- unique(subset(tmp, select=c(MALE_RID,FEMALE_RID,PTY_RUN)))
	for(ii in 1:7)
	{		
		indir		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun23'
		# load dfr and phs
		load( file.path(indir, paste0('ptyr',rtpdm[ii,PTY_RUN],'_trees.rda')) )
		# setup plotting
		ids			<- c(rtpdm[ii, MALE_RID],rtpdm[ii, FEMALE_RID])
		dfs			<- subset(dfr, select=c(W_FROM, W_TO, IDX))
		dfs[, MALE_RID:=ids[1]]
		dfs[, FEMALE_RID:=ids[2]]
		dfs			<- merge(dfs, subset(rpw, GROUP=='TYPE_RAW'), by=c('MALE_RID','FEMALE_RID','W_FROM','W_TO'), all.x=1)	
		dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', PTY_RUN, '\nwindow ', W_FROM,'-', W_TO,'\npmf',PATHS_MF,' pfm',PATHS_FM, ' ',ADJACENT,' ',TYPE, '\n', round(PATRISTIC_DISTANCE, d=5), sep='')]]
		plot.file	<- paste0(outfile.base, 'checkcloseother/todi_pairs_170610_run_', rtpdm[ii, PTY_RUN],'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
		invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10, tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'))					
	}
}

RakaiFull.analyze.couples.todi.170811.shorter.than.fulllength<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
		
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_withmetadata.rda'
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_"
	load(infile)	
	setkey(rca, MALE_RID, FEMALE_RID)
	rca[, PAIRID:=seq_len(nrow(rca))]
	
	#	average number of windows per couple with enough data
	subset(rca, !grepl('insufficient', SELECT))[, mean(NEFF)]
	#	9.76
	
	#	select pairs with more than 15 windows 	
	neff.select	<- 15 
	replicate.n	<- 100
	de			<- subset(rca, NEFF>neff.select, c(FEMALE_RID, MALE_RID, PTY_RUN, SELECT, NEFF))
	setnames(de, c('FEMALE_RID','MALE_RID','PTY_RUN'), c('FEMALE_RID_Q','MALE_RID_Q','PTY_RUN_Q'))
	#	add how many non-overlapping windows to be chosen at random
	de			<- merge(de, de[, list(N_CHOOSE=c(1,3,5,10)), by=c('FEMALE_RID_Q','MALE_RID_Q','PTY_RUN_Q')],by=c('FEMALE_RID_Q','MALE_RID_Q','PTY_RUN_Q'))	
	de			<- merge(de, as.data.table(expand.grid(N_CHOOSE=c(1,3,5,10), REP=seq_len(replicate.n))), by=c('N_CHOOSE'), allow.cartesian=TRUE)
	#	
	
	#
	#	select windows at random, 
	#	then calculate most likely assignment across subsampled windows
	#	for TYPE_PAIR_TODI2 and TYPE_NETWORK_SCORES
	set.seed(42)
	trmw.close.brl			<- 0.035 
	trmw.distant.brl		<- 0.08
	prior.keff				<- prior.keff.dir <- 2
	prior.neff				<- prior.neff.dir <- 3 
	prior.calibrated.prob	<- 0.66	
	get.groups				<- c('TYPE_PAIR_TODI2','TYPE_NETWORK_SCORES')	
	de						<- de[, {
				#FEMALE_RID_Q 	<- 'K110876'; MALE_RID_Q <- 'F110995'; PTY_RUN_Q<- 318; N_CHOOSE<- 10
				#	subset to query individual
				female_rid	<- FEMALE_RID_Q; male_rid<- MALE_RID_Q; pty.run<- PTY_RUN_Q	
				tmp			<- unique(subset(rpw, GROUP=='TYPE_RAW' & FEMALE_RID==female_rid & MALE_RID==male_rid & PTY_RUN==pty.run, c(FEMALE_RID,MALE_RID,PTY_RUN,W_FROM,W_TO)))
				#	determine distinct windows
				tmp[, NEXT_NON_OVERLAPPING:= NA_integer_]
				set(tmp, 1L, 'NEXT_NON_OVERLAPPING', 1L)
				i			<- tail(which(!is.na(tmp$NEXT_NON_OVERLAPPING)))
				while( length(i) && i<nrow(tmp) )
				{
					set(tmp, tmp[, which(W_FROM<tmp$W_TO[i] & W_TO>tmp$W_TO[i])], 'NEXT_NON_OVERLAPPING', 0L)
					i	<- tmp[, which(W_FROM>=tmp$W_TO[i])]
					if(length(i))
					{
						i	<- i[1] 
						set(tmp, i, 'NEXT_NON_OVERLAPPING', 1L)		
					}		
				}
				tmp		<- tmp[sample(which(NEXT_NON_OVERLAPPING==1), N_CHOOSE),]
				tmp		<- merge(subset(rpw, GROUP=='TYPE_RAW'), tmp, by=c('FEMALE_RID','MALE_RID','PTY_RUN','W_FROM','W_TO'))
				setnames(tmp, c('MALE_RID','FEMALE_RID','TYPE','PATHS_MF','PATHS_FM'), c('ID1','ID2','TYPE_RAW','PATHS_12','PATHS_21'))
				set(tmp, NULL, 'TYPE_RAW', tmp[, gsub('fm','21',gsub('mf','12',TYPE_RAW))])
				set(tmp, NULL, 'GROUP', NULL)
				dwin	<- phsc.get.basic.pairwise.relationships(tmp, trmw.close.brl, trmw.distant.brl, verbose=FALSE)
				setnames(dwin, 'TYPE_BASIC', 'TYPE_DIR_TODI7x3')	
				dwin	<- phsc.get.pairwise.relationships(dwin, get.groups=get.groups, make.pretty.labels=FALSE)
				setnames(dwin, 'TYPE_DIR_TODI7x3', 'TYPE_BASIC')
				dl		<- phsc.get.pairwise.relationships.keff.and.neff(dwin, get.groups)
				#	get most likely assignment on subsampled windows
				dl		<- subset(dl, GROUP%in%get.groups)[, 	{
																	z<- which(KEFF==max(KEFF))
																	list(TYPE_MLE=TYPE[ifelse(length(z)>1, sample(z,1), z)])
																}, by=c('ID1','ID2','GROUP')]
				set(dl, NULL, 'TYPE_MLE', dl[, gsub('21','fm',gsub('12','mf',TYPE_MLE))])
				set(dl, NULL, c('ID1','ID2'), NULL)
			}, by=c('FEMALE_RID_Q','MALE_RID_Q','PTY_RUN_Q','N_CHOOSE','REP')]	
	setnames(de, c('FEMALE_RID_Q','MALE_RID_Q','PTY_RUN_Q'), c('FEMALE_RID','MALE_RID','PTY_RUN'))
	#	add most likely assignment on real data
	tmp	<- unique(subset(de, select=c('FEMALE_RID','MALE_RID','PTY_RUN','GROUP')))
	tmp	<- merge(tmp, rplkl, by=c('FEMALE_RID','MALE_RID','PTY_RUN','GROUP'))
	tmp	<- tmp[, 	{
						z<- which(KEFF==max(KEFF))
						list(TYPE_MLE_FULLLEN=TYPE[ifelse(length(z)>1, sample(z,1), z)])
					}, by=c('FEMALE_RID','MALE_RID','PTY_RUN','GROUP')]
	de	<- merge(de, tmp, by=c('FEMALE_RID','MALE_RID','PTY_RUN','GROUP'))	
	save(de, file=paste0(outfile.base,'fulllen_vs_shorterlen.rda'))
	#	summarise misclassification error for linked couples
	tmp	<- subset(rca, grepl('likely a pair',SELECT), select=c(FEMALE_RID,MALE_RID,PTY_RUN))
	tmp	<- merge(tmp, de, by=c('FEMALE_RID','MALE_RID','PTY_RUN'))
	tmp	<- tmp[, list(E= length(which(TYPE_MLE!=TYPE_MLE_FULLLEN))/length(TYPE_MLE) ), by=c('GROUP','N_CHOOSE','REP')]
	tmp	<- tmp[, list(P=paste0('p',c(0.025,0.25,0.5,0.75,0.975)), Q=quantile(E, prob=c(0.025,0.25,0.5,0.75,0.975))), by=c('GROUP','N_CHOOSE')]
	tmp	<- dcast.data.table(tmp, GROUP+N_CHOOSE~P, value.var='Q')
	tmp[, ANS:= paste0(round(p0.5*100, d=1), ' [',round(p0.025*100, d=1),', ',round(p0.975*100, d=1),']')]
	tmp	<- dcast.data.table(tmp, GROUP~N_CHOOSE, value.var='ANS')
	write.csv(tmp, file=paste0(outfile.base,'fulllen_vs_shorterlen_summary_among_linked_couples.csv'))
	#                 GROUP                 1                 3                 5                10
	#1: TYPE_NETWORK_SCORES 56.5 [39.1, 73.9] 47.8 [30.4, 71.8] 43.5 [30.4, 60.9] 34.8 [17.4, 47.8]
	#2:     TYPE_PAIR_TODI2     8.7 [0, 21.7]       4.3 [0, 13]        0 [0, 8.7]        0 [0, 4.3]
	
	#	summarise misclassification error for all couples
	tmp	<- de[, list(E= length(which(TYPE_MLE!=TYPE_MLE_FULLLEN))/length(TYPE_MLE) ), by=c('GROUP','N_CHOOSE','REP')]
	tmp	<- tmp[, list(P=paste0('p',c(0.025,0.25,0.5,0.75,0.975)), Q=quantile(E, prob=c(0.025,0.25,0.5,0.75,0.975))), by=c('GROUP','N_CHOOSE')]
	#	make table
	tmp	<- dcast.data.table(tmp, GROUP+N_CHOOSE~P, value.var='Q')
	tmp[, ANS:= paste0(round(p0.5*100, d=1), ' [',round(p0.025*100, d=1),', ',round(p0.975*100, d=1),']')]
	tmp	<- dcast.data.table(tmp, GROUP~N_CHOOSE, value.var='ANS')
	write.csv(tmp, file=paste0(outfile.base,'fulllen_vs_shorterlen_summary.csv'))
	#	              GROUP              1             3            5          10
	#1: TYPE_NETWORK_SCORES 	34 [22, 44] 30 [18.9, 39] 26 [18, 37] 22 [12, 29]
	#2:     TYPE_PAIR_TODI2  	 10 [4, 17]     6 [2, 12]   4 [2, 10]    2 [0, 6]
	
}

RakaiFull.analyze.couples.todi.170811.birdseyeview<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170522/todi_couples_170522_withmetadata.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170522/todi_couples_170522_"
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170610/todi_couples_170610_cl3_prior23_withmetadata.rda'
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170610/todi_couples_170610_cl3_prior23_"
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_withmetadata.rda'
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_"
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10_withmetadata.rda'
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10_"
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_withmetadata.rda'
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_"	
	load(infile)	
	setkey(rca, MALE_RID, FEMALE_RID)
	rca[, PAIRID:=seq_len(nrow(rca))]
	rca					<- subset(rca, !grepl('insufficient',SELECT))
	#	for each pair:
	#	get posterior prob for likely pair
	rpp					<- subset(rplkl, GROUP=='TYPE_PAIR_TODI2' & TYPE=='linked')
	#	merge PAIRID
	rpp					<- merge(rpp,unique(subset(rca, select=c(MALE_RID,FEMALE_RID,PTY_RUN,PAIRID))),by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	#	sort by posterior median highest for likely pair
	#rpp					<- rpp[order(qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA)),]
	#	sort by posterior mode highest for likely pair
	rpp					<- rpp[order((POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-2)),]	
	set(rpp, NULL, 'PAIRID',rpp[, factor(PAIRID, levels=rpp$PAIRID)])
		
	rpd					<- subset(rplkl, GROUP=='TYPE_DIR_TODI2')
	rpd					<- merge(rpd, rpd[, list(TYPE=TYPE[which.max(KEFF)]), by=c('MALE_RID','FEMALE_RID','PTY_RUN')], by=c('MALE_RID','FEMALE_RID','PTY_RUN','TYPE'))
	setnames(rpd, c('TYPE','POSTERIOR_ALPHA','POSTERIOR_BETA'), c('TYPE_DIR','POSTERIOR_ALPHA_DIR','POSTERIOR_BETA_DIR'))
	rpd					<- merge(subset(rpp, select=c(MALE_RID, FEMALE_RID, PTY_RUN, PAIRID, TYPE, POSTERIOR_ALPHA, POSTERIOR_BETA)), rpd, by=c('MALE_RID','FEMALE_RID','PTY_RUN'), all.x=TRUE)
	set(rpd, rpd[, which((POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-2)<0.66)],c('POSTERIOR_ALPHA_DIR','POSTERIOR_BETA_DIR'), NA_real_)	
	tmp					<- subset(rplkl, GROUP=='TYPE_NETWORK_SCORES')
	rpd					<- merge(rpd, tmp[, list(TYPE_SCORE=TYPE[which.max(KEFF)]), by=c('MALE_RID','FEMALE_RID','PTY_RUN')], by=c('MALE_RID','FEMALE_RID','PTY_RUN'))	
	set(rpd, rpd[, which(grepl('disconnected|ambiguous',TYPE_SCORE))],c('POSTERIOR_ALPHA_DIR','POSTERIOR_BETA_DIR'), NA_real_)
	set(rpd, rpd[, which(grepl('disconnected|ambiguous',TYPE_SCORE))],c('TYPE_DIR'), NA_character_)
	
	p3	<- ggplot(rpp, aes(x=PAIRID)) +
			geom_segment(aes(xend=PAIRID, y=qbeta(0.025, POSTERIOR_ALPHA, POSTERIOR_BETA),yend=qbeta(0.975, POSTERIOR_ALPHA, POSTERIOR_BETA)), colour='grey50') +
			#geom_point(aes(y=qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA)), colour='black') +
			geom_point(aes(y=(POSTERIOR_ALPHA-1)/(POSTERIOR_ALPHA+POSTERIOR_BETA-2)), colour='black') +
			scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.25), labels=scales::percent) +					
			theme_bw() + 			
			theme(axis.text.x=element_blank()) +
			theme(axis.ticks.x=element_blank()) +
			theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
			labs(x='',y='',fill='')
	p5	<- ggplot(rpd, aes(x=PAIRID)) +
			geom_segment(aes(xend=PAIRID, y=qbeta(0.025, POSTERIOR_ALPHA_DIR, POSTERIOR_BETA_DIR),yend=qbeta(0.975, POSTERIOR_ALPHA_DIR, POSTERIOR_BETA_DIR)), colour='grey50') +
			geom_point(aes(y=(POSTERIOR_ALPHA_DIR-1)/(POSTERIOR_ALPHA_DIR+POSTERIOR_BETA_DIR-2), colour=TYPE_DIR)) +
			scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.25), labels=scales::percent) +					
			scale_colour_manual(values=c('mf'='steelblue2','fm'='hotpink2')) +
			theme_bw() + 			
			theme(axis.text.x=element_blank()) +
			theme(axis.ticks.x=element_blank()) +
			theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
			labs(x='',y='',fill='') + 
			guides(colour=FALSE)
			
	
	#	for each pair:
	#	make pair topology assignments: ancestral, intermingled, sibling, disconnected
	if(1)
	{
		rplkl2				<- subset(rplkl, GROUP=='TYPE_BASIC')
		rplkl2[, TYPE_TO:= 'disconnected or with intermediate']
		set(rplkl2, rplkl2[,which(grepl('no intermediate', TYPE) & grepl('chain', TYPE))], 'TYPE_TO', 'ancestral & no intermediate')
		set(rplkl2, rplkl2[,which(grepl('no intermediate', TYPE) & grepl('intermingled', TYPE))], 'TYPE_TO', 'intermingled & no intermediate')
		set(rplkl2, rplkl2[,which(grepl('no intermediate', TYPE) & grepl('other', TYPE))], 'TYPE_TO', 'adjacent & no intermediate')
		rplkl2				<- rplkl2[, list(N=N[1],NEFF=NEFF[1],K=sum(K),KEFF=sum(KEFF)), by=c('MALE_RID','FEMALE_RID','PTY_RUN','TYPE_TO')]
		set(rplkl2, NULL, 'TYPE_TO', rplkl2[, factor(TYPE_TO, levels=c('ancestral & no intermediate','intermingled & no intermediate','adjacent & no intermediate','disconnected or with intermediate'))])
		#	copy sorted PAIRID to rplkl2
		rplkl2				<- merge(rplkl2, unique(subset(rpp, select=c(MALE_RID,FEMALE_RID,PTY_RUN,PAIRID))), by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
		p2	<- ggplot(rplkl2, aes(x=PAIRID, y=KEFF, fill=TYPE_TO)) +
				geom_bar(stat='identity',position='stack') +
				scale_y_continuous(expand=c(0,0)) +
				scale_fill_manual(values=c("ancestral & no intermediate"=brewer.pal(11, 'PuOr')[2], "intermingled & no intermediate"=brewer.pal(11, 'PuOr')[4], 'adjacent & no intermediate'=rev(brewer.pal(11, 'PuOr'))[c(3)], "disconnected or with intermediate"=rev(brewer.pal(11, 'RdGy'))[4])) +
				labs(x='', y='', fill='') +
				theme_bw() + 
				theme(legend.position='bottom') +
				theme(axis.ticks.x=element_blank()) +
				theme(panel.border=element_blank()) +
				theme(axis.text.x=element_blank()) +
				theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank())
	}
	if(0)
	{
		rplkl2				<- subset(rplkl, GROUP=='TYPE_BASIC')
		rplkl2[, TYPE_TO:= 'disconnected/ not close']
		set(rplkl2, rplkl2[,which(grepl('close', TYPE) & grepl('no intermediate', TYPE) & grepl('chain', TYPE))], 'TYPE_TO', 'close ancestral')
		set(rplkl2, rplkl2[,which(grepl('close', TYPE) & grepl('no intermediate', TYPE) & grepl('intermingled', TYPE))], 'TYPE_TO', 'close intermingled')
		set(rplkl2, rplkl2[,which(grepl('close', TYPE) & grepl('no intermediate', TYPE) & grepl('other', TYPE))], 'TYPE_TO', 'close sibling')
		rplkl2				<- rplkl2[, list(N=N[1],NEFF=NEFF[1],K=sum(K),KEFF=sum(KEFF)), by=c('MALE_RID','FEMALE_RID','PTY_RUN','TYPE_TO')]
		set(rplkl2, NULL, 'TYPE_TO', rplkl2[, factor(TYPE_TO, levels=c('close ancestral','close intermingled','close sibling','disconnected/ not close'))])
		#	copy sorted PAIRID to rplkl2
		rplkl2				<- merge(rplkl2, unique(subset(rpp, select=c(MALE_RID,FEMALE_RID,PTY_RUN,PAIRID))), by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
		#
		p2	<- ggplot(rplkl2, aes(x=PAIRID, y=KEFF, fill=TYPE_TO)) +
				geom_bar(stat='identity',position='stack') +
				scale_y_continuous(expand=c(0,0)) +
				scale_fill_manual(values=c("close ancestral"=brewer.pal(11, 'PuOr')[2], "close intermingled"=brewer.pal(11, 'PuOr')[4], 'close sibling'=rev(brewer.pal(11, 'PuOr'))[c(3)], "disconnected/ not close"=rev(brewer.pal(11, 'RdGy'))[4])) +
				labs(x='', y='', fill='') +
				theme_bw() + 
				theme(legend.position='bottom') +
				theme(axis.ticks.x=element_blank()) +
				theme(panel.border=element_blank()) +
				theme(axis.text.x=element_blank()) +
				theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank())
	}	
	#	for each pair:
	#	get phyloscanner distances confidence intervals and median values
	rpw2	<- subset(rpw, GROUP=='TYPE_BASIC')[, list(	PD_CL=quantile(PATRISTIC_DISTANCE,p=0.025),
														PD_IL=quantile(PATRISTIC_DISTANCE,p=0.25),
														PD_M=median(PATRISTIC_DISTANCE),
														PD_IU=quantile(PATRISTIC_DISTANCE,p=0.75),
														PD_CU=quantile(PATRISTIC_DISTANCE,p=0.975)
														),by=c('MALE_RID','FEMALE_RID','PTY_RUN')]
	rpw2	<- merge(rpw2, unique(subset(rpp, select=c(MALE_RID,FEMALE_RID,PTY_RUN,PAIRID))), by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	set(rpw2, rpw2[, which(PD_M>0.2)], 'PD_M', 0.2)
	p1		<- ggplot(rpw2, aes(x=PAIRID)) +
			geom_segment(aes(xend=PAIRID, y=PD_CL,yend=PD_CU), colour='grey50') +
			geom_point(aes(y=PD_M), colour='black') +
			scale_y_continuous(expand=c(0,0), labels=scales::percent) +
			theme_bw() + 
			theme(axis.text.x=element_blank()) +
			theme(axis.ticks.x=element_blank()) +
			theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
			labs(x='',y='',fill='') +	
			coord_cartesian(ylim=c(0,0.201))
	#	for each pair:
	#	ambiguity measure -- proportion of windows that are not assigned the most likely relationship type
	rplkl2	<- subset(rplkl, GROUP=='TYPE_PAIR_TODI2')
	rpa		<- rplkl2[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('MALE_RID','FEMALE_RID','PTY_RUN')]	
	rpa		<- merge(rpa, rplkl2, by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	rpa[, TYPE_AGREE:= factor(TYPE==TYPE_MLE, levels=c(TRUE,FALSE), labels=c('y','n'))]
	rpa		<- rpa[, list(N=N[1], NEFF=NEFF[1], K=sum(K), KEFF=sum(KEFF), PEFF=sum(KEFF)/NEFF[1]), by=c('MALE_RID','FEMALE_RID','PTY_RUN','TYPE_AGREE')]
	rpa		<- subset(rpa, TYPE_AGREE=='n')
	rpa		<- merge(rpa, unique(subset(rpp, select=c(MALE_RID,FEMALE_RID,PTY_RUN,PAIRID))), by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	p4		<- ggplot(rpa, aes(x=PAIRID)) +			
			geom_point(aes(y=PEFF), colour='black') +
			scale_y_continuous(expand=c(0,0), limits=c(-5e-3,1), labels=scales::percent, breaks=seq(0,1,0.1)) +
			theme_bw() + 
			theme(axis.text.x=element_blank()) +
			theme(axis.ticks.x=element_blank()) +
			theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
			labs(x='',y='',fill='') +
			coord_cartesian(ylim=c(-0.002,0.51))
	plot.file	<- paste0(outfile.base, 'birdseyeview.pdf')
	pdf(file=plot.file, h=8.27*1.25, w=11.69*1.5, useDingbats=FALSE)
	grid.newpage()	
	pushViewport(viewport(layout=grid.layout(4, 1, heights=unit(c(3.5,2,2,2), "null"))))   	
	print(p1, vp=viewport(layout.pos.row=2, layout.pos.col=1))
	print(p2, vp=viewport(layout.pos.row=1, layout.pos.col=1))
	print(p3, vp=viewport(layout.pos.row=3, layout.pos.col=1))
	print(p5, vp=viewport(layout.pos.row=4, layout.pos.col=1))
	#grid.draw(p3)
	dev.off()	
	
}

RakaiFull.analyze.couples.todi.170811.DI.vs.TODI.vs.DIR.xaxis.consensus<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_withmetadata.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_"		
	load(infile)
	
	#	load raw genetic distances 
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_rawgendist.rda"
	load(infile)		
	gd		<- subset(gd, grepl('^PG',TAXA) & grepl('^PG',TAXA2))
	set(gd, NULL, 'TAXA', gd[, gsub('_WTSI.*','',TAXA)])
	set(gd, NULL, 'TAXA2', gd[, gsub('_WTSI.*','',TAXA2)])	
	setnames(gd, c('TAXA','TAXA2'), c('MALE_TAXA','FEMALE_TAXA'))
	tmp		<- copy(gd)
	setnames(tmp, c('MALE_TAXA','FEMALE_TAXA'), c('FEMALE_TAXA','MALE_TAXA'))
	gd		<- rbind(gd, tmp, use.names=TRUE)
	tmp		<- subset(rp, !is.na(MALE_TAXA) & !is.na(FEMALE_TAXA), select=c(MALE_RID,FEMALE_RID,MALE_TAXA,FEMALE_TAXA))
	gd		<- merge(gd, tmp, by=c('FEMALE_TAXA','MALE_TAXA'))
	gd		<- gd[, list(CONS_GDRW=mean(CONS_GDRW)), by=c('MALE_RID','FEMALE_RID')]
	
	
	#	per window:
	#	make window topology assignments: TYPE_PAIR_TO ancestral/intermingled; withintermediate; other (include adjacent other here)
	#	make boxplot
	rpw2	<- subset(rpw, GROUP=='TYPE_BASIC')
	rpw2[, TYPE_TO:= 'disconnected or with intermediate']
	set(rpw2, rpw2[,which(grepl('no intermediate', TYPE) & grepl('chain', TYPE))], 'TYPE_TO', 'ancestral & no intermediate')
	set(rpw2, rpw2[,which(grepl('no intermediate', TYPE) & grepl('intermingled', TYPE))], 'TYPE_TO', 'intermingled & no intermediate')
	set(rpw2, rpw2[,which(grepl('no intermediate', TYPE) & grepl('other', TYPE))], 'TYPE_TO', 'adjacent & no intermediate')	
	set(rpw2, NULL, 'TYPE_TO', rpw2[, factor(TYPE_TO, levels=c('ancestral & no intermediate','intermingled & no intermediate','adjacent & no intermediate','disconnected or with intermediate'))])
	rpw2[, PATRISTIC_DISTANCE_LOG:= log10(PATRISTIC_DISTANCE)]
	rpw2[, PATRISTIC_DISTANCE_LOGC:= cut(PATRISTIC_DISTANCE_LOG, breaks=log10(c(1e-12, 0.005, 0.01, 0.02, 0.035, 0.08, 2e3)),labels=c('<0.5%', '0.5%-1%', '1%-2%', '2%-3.5%', '3.5%-8%', '>8%'))]
	rpw2[, PATRISTIC_DISTANCE_LOGC2:= cut(PATRISTIC_DISTANCE_LOG, breaks=log10(c(1e-12, 0.035, 0.08, 2e3)),labels=c('<3.5%', '3.5%-8%', '>8%'))]
	rpw2	<- merge(rpw2, gd, by=c('MALE_RID','FEMALE_RID'), all.x=TRUE)
	rpw2[, CONS_GDRW_LOGC:= cut(CONS_GDRW, breaks=c(1e-12, 0.005, 0.01, 0.02, 0.035, 0.08, 2e3),labels=c('<0.5%', '0.5%-1%', '1%-2%', '2%-3.5%', '3.5%-8%', '>8%'))]
	set(rpw2, rpw2[, which(is.na(CONS_GDRW_LOGC))], 'CONS_GDRW_LOGC', 'consensus\n<700nt')
	#
	rpw2[, table(TYPE_TO)]
	rpw2[, table(TYPE_TO)/nrow(rpw2)]
	#      ancestral & no intermediate    intermingled & no intermediate        adjacent & no intermediate disconnected or with intermediate 
	#                    0.34229600                        0.14258437                        0.05163276                        0.46348687
	
	#
	#	plot topology assignment by distance
	#		
	cols.typet			<- c("ancestral & no intermediate"=brewer.pal(11, 'PuOr')[2], "intermingled & no intermediate"=brewer.pal(11, 'PuOr')[4], 'adjacent & no intermediate'=rev(brewer.pal(11, 'PuOr'))[c(3)], "disconnected or with intermediate"=rev(brewer.pal(11, 'RdGy'))[4])	
	ggplot(subset(rpw2, !grepl('consensus',CONS_GDRW_LOGC)), aes(x=CONS_GDRW_LOGC, fill=TYPE_TO)) +
			geom_bar(position='fill') +
			scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), expand=c(0,0)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			labs(x='\ngenetic distance between consensus sequences\n(subst/site)', y='subtree relationships\nacross genomic windows\n', fill='pairwise\ntopological relationship') +
			guides(fill=FALSE)
	ggsave(file=paste0(outfile.base, 'topocons_pairtypes_prop.pdf'), w=5, h=2.25)
	ggplot(rpw2, aes(x=CONS_GDRW_LOGC, fill=TYPE_TO)) +
			geom_bar(position='stack') +
			scale_y_continuous(breaks=seq(0,20e3,2e3), expand=c(0,0)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			labs(x='\ngenetic distance between consensus sequences\n(subst/site)', y='subtree relationships\nacross genomic windows\n', fill='pairwise\ntopological relationship') +
			guides(fill=FALSE)
	ggsave(file=paste0(outfile.base, 'topocons_pairtypes_count.pdf'), w=6, h=2.25)	
}

RakaiFull.analyze.couples.todi.170811.DI.vs.TODI.vs.DIR<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170522/todi_couples_170522_withmetadata.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170522/todi_couples_170522_"
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170610/todi_couples_170610_cl3_prior23_withmetadata.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170610/todi_couples_170610_cl3_prior23_"	
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10_withmetadata.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10_"	
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_withmetadata.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_"	
	
	load(infile)
	#	per pair:
	#	make pair topology assignments: TYPE_PAIR_TO ancestral/intermingled; other; ambiguous as usual
	#	stratify mean patristic distance
	#	make boxplot
	#	OK if not too coarse, just a pain
	
	#	per window:
	#	make window topology assignments: TYPE_PAIR_TO ancestral/intermingled; withintermediate; other (include adjacent other here)
	#	make boxplot
	rpw2	<- subset(rpw, GROUP=='TYPE_BASIC')
	#rpw2[, TYPE_TO:= 'disconnected']
	#set(rpw2, rpw2[,which(grepl('no intermediate', TYPE))], 'TYPE_TO', 'ancestral/intermingled')
	#set(rpw2, rpw2[,which(grepl('with intermediate', TYPE))], 'TYPE_TO', 'with intermediate')
	#set(rpw2, rpw2[,which(ADJACENT & grepl('other', TYPE))], 'TYPE_TO', 'sibling')
	#set(rpw2, NULL, 'TYPE_TO', rpw2[, factor(TYPE_TO, levels=c("ancestral/intermingled", "with intermediate", 'sibling', "disconnected"))])
	rpw2[, TYPE_TO:= 'disconnected or with intermediate']
	set(rpw2, rpw2[,which(grepl('no intermediate', TYPE) & grepl('chain', TYPE))], 'TYPE_TO', 'ancestral & no intermediate')
	set(rpw2, rpw2[,which(grepl('no intermediate', TYPE) & grepl('intermingled', TYPE))], 'TYPE_TO', 'intermingled & no intermediate')
	set(rpw2, rpw2[,which(grepl('no intermediate', TYPE) & grepl('other', TYPE))], 'TYPE_TO', 'adjacent & no intermediate')	
	set(rpw2, NULL, 'TYPE_TO', rpw2[, factor(TYPE_TO, levels=c('ancestral & no intermediate','intermingled & no intermediate','adjacent & no intermediate','disconnected or with intermediate'))])
	rpw2[, PATRISTIC_DISTANCE_LOG:= log10(PATRISTIC_DISTANCE)]
	rpw2[, PATRISTIC_DISTANCE_LOGC:= cut(PATRISTIC_DISTANCE_LOG, breaks=log10(c(1e-12, 0.005, 0.01, 0.02, 0.035, 0.08, 2e3)),labels=c('<0.5%', '0.5%-1%', '1%-2%', '2%-3.5%', '3.5%-8%', '>8%'))]
	rpw2[, PATRISTIC_DISTANCE_LOGC2:= cut(PATRISTIC_DISTANCE_LOG, breaks=log10(c(1e-12, 0.035, 0.08, 2e3)),labels=c('<3.5%', '3.5%-8%', '>8%'))]	
	#
	rpw2[, table(TYPE_TO)]
	rpw2[, table(TYPE_TO)/nrow(rpw2)]
	#      ancestral & no intermediate    intermingled & no intermediate        adjacent & no intermediate disconnected or with intermediate 
    #                    0.34229600                        0.14258437                        0.05163276                        0.46348687

	#
	#	plot topology assignment by distance
	#		
	#cols.typet			<- c(brewer.pal(11, 'PuOr')[2], rev(brewer.pal(11, 'PuOr'))[c(3,4)], rev(brewer.pal(11, 'RdGy'))[4])
	#names(cols.typet)	<- c("ancestral/intermingled", "with intermediate", 'sibling', "disconnected")
	cols.typet			<- c("ancestral & no intermediate"=brewer.pal(11, 'PuOr')[2], "intermingled & no intermediate"=brewer.pal(11, 'PuOr')[4], 'adjacent & no intermediate'=rev(brewer.pal(11, 'PuOr'))[c(3)], "disconnected or with intermediate"=rev(brewer.pal(11, 'RdGy'))[4])	
	ggplot(rpw2, aes(x=PATRISTIC_DISTANCE_LOGC, fill=TYPE_TO)) +
			geom_bar(position='fill') +
			scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), expand=c(0,0)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			labs(x='\nsubtree distance\n(subst/site)', y='subtree relationships\nacross genomic windows\n', fill='pairwise\ntopological relationship') +
			guides(fill=FALSE)
	ggsave(file=paste0(outfile.base, 'topodist_pairtypes_prop.pdf'), w=5, h=2.25)
	ggplot(rpw2, aes(x=PATRISTIC_DISTANCE_LOGC, fill=TYPE_TO)) +
			geom_bar(position='stack') +
			scale_y_continuous(breaks=seq(0,20e3,2e3), expand=c(0,0)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			labs(x='\nsubtree distance\n(subst/site)', y='subtree relationships\nacross genomic windows\n', fill='pairwise\ntopological relationship') +
			guides(fill=FALSE)
	ggsave(file=paste0(outfile.base, 'topodist_pairtypes_count.pdf'), w=5, h=2.25)
	ggplot(rpw2, aes(x=PATRISTIC_DISTANCE_LOGC2, fill=TYPE_TO)) +
			geom_bar(position='fill') +
			scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), expand=c(0,0)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			labs(x='\npatristic distance between WH clades\n(subst/site)', y='genomic windows\nacross couples\n', fill='pairwise\ntopological relationship')
	ggsave(file=paste0(outfile.base, 'topodist_pairtypes_prop_3diststrat.pdf'), w=4, h=5)
	
	#
	#	table DI vs TODI
	rex			<- subset(rplkl, GROUP=='TYPE_PAIR_DI')[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('MALE_RID','FEMALE_RID','PTY_RUN')]
	rex			<- merge(subset(rex, TYPE_MLE=='distant', c('MALE_RID','FEMALE_RID')), subset(rplkl, GROUP=='TYPE_PAIR_DI' & TYPE=='distant'), by=c('MALE_RID','FEMALE_RID'), all.x=1)
	rex[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]	
	rex			<- rex[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)], POSTERIOR_SCORE=POSTERIOR_SCORE[which.max(POSTERIOR_SCORE)]), by=c('MALE_RID','FEMALE_RID')]	
	rex			<- subset(rex, POSTERIOR_SCORE>confidence.cut)
	rex[, SELECT_DI:= 'couple most likely not a pair']	
	rtp			<- subset(rplkl, GROUP=='TYPE_PAIR_DI')[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('MALE_RID','FEMALE_RID','PTY_RUN')]
	rtp			<- merge(subset(rtp, TYPE_MLE=='close', c('MALE_RID','FEMALE_RID')), subset(rplkl, GROUP=='TYPE_PAIR_DI' & TYPE=='close'), by=c('MALE_RID','FEMALE_RID'), all.x=1)
	rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]	
	rtp			<- rtp[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)], POSTERIOR_SCORE=POSTERIOR_SCORE[which.max(POSTERIOR_SCORE)]), by=c('MALE_RID','FEMALE_RID')]	
	rtp			<- subset(rtp, POSTERIOR_SCORE>confidence.cut)
	rtp[, SELECT_DI:= 'couple most likely a pair']
	tmp			<- rbind(rex, rtp)
	tmp[, POSTERIOR_SCORE:=NULL]
	rca			<- merge(rca, tmp, by=c('MALE_RID','FEMALE_RID','PTY_RUN'), all.x=1)
	set(rca, rca[, which(SELECT=='insufficient deep sequence data for at least one partner of couple')], 'SELECT_DI', 'insufficient deep sequence data for at least one partner of couple')
	set(rca, rca[, which(is.na(SELECT_DI))], 'SELECT_DI', 'couple ambiguous if pair or not pair')

	tmp	<- subset(rca, select=c(MALE_RID, FEMALE_RID, PTY_RUN, SELECT, SELECT_DI))
	tmp	<- subset(tmp, SELECT!='insufficient deep sequence data for at least one partner of couple')
	set(tmp, tmp[, which(grepl('couple most likely a pair', SELECT))], 'SELECT', 'couple most likely a pair')
	
	set(tmp, NULL, 'SELECT', tmp[, factor(SELECT, levels=c('couple most likely not a pair','couple ambiguous if pair or not pair','couple most likely a pair'))])
	set(tmp, NULL, 'SELECT_DI', tmp[, factor(SELECT_DI, levels=c('couple most likely not a pair','couple ambiguous if pair or not pair','couple most likely a pair'))])
	tmp[, table(SELECT_DI, SELECT)]
	#                                      SELECT
	#SELECT_DI                              couple most likely not a pair couple ambiguous if pair or not pair couple most likely a pair	TOTAL
  	#couple most likely not a pair                                  135                                    0                         0		135
  	#couple ambiguous if pair or not pair                             5                                    5                         0		10
  	#couple most likely a pair                                        6                                   14                       143		163
	#TOTAL															146									  19					   143

	#20/163=12.3%		143/163= 87.7%

	#	add direction
	set(rpw2, NULL, 'TYPE_TO', rpw2[, as.character(TYPE_TO)]) 
	set(rpw2, rpw2[,which(grepl('chain', TYPE) & grepl('no intermediate', TYPE))], 'TYPE_TO', 'ancestral')
	set(rpw2, rpw2[,which(grepl('intermingled', TYPE) & grepl('no intermediate', TYPE))], 'TYPE_TO', 'intermingled')
	cols.typet			<- c(brewer.pal(11, 'PuOr')[c(2,4)], rev(brewer.pal(11, 'PuOr'))[c(3,4)], rev(brewer.pal(11, 'RdGy'))[4])
	names(cols.typet)	<- c("ancestral", "intermingled", "with intermediate", 'sibling', "disconnected")
	set(rpw2, NULL, 'TYPE_TO', rpw2[, factor(TYPE_TO, levels=c("ancestral","intermingled", "with intermediate", 'sibling', "disconnected"))])
	ggplot(rpw2, aes(x=PATRISTIC_DISTANCE_LOGC, fill=TYPE_TO)) +
			geom_bar(position='fill') +
			scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), expand=c(0,0)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			labs(x='\npatristic distance between WH clades\n(subst/site)', y='genomic windows\nacross couples\n', fill='pairwise\ntopological relationship')
	ggsave(file=paste0(outfile.base, 'topodist_pairdirtypes_prop.pdf'), w=7, h=4.5)
	ggplot(rpw2, aes(x=PATRISTIC_DISTANCE_LOGC2, fill=TYPE_TO)) +
			geom_bar(position='fill') +
			scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), expand=c(0,0)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			labs(x='\npatristic distance between WH clades\n(subst/site)', y='genomic windows\nacross couples\n', fill='pairwise\ntopological relationship')
	ggsave(file=paste0(outfile.base, 'topodist_pairdirtypes_prop_3diststrat.pdf'), w=4, h=5)
	
	
	#	TODO extend definition of COUP_SC==seropos: 
	#	if one partner seropos for more than >7 years before other partner found seropos
	#
	#	validate direction on serodiscordant pairs
	tmp	<- subset(rca, COUP_SC%in%c('M->F','F->M'))
	tmp[, table(SELECT)]
	#SELECT
	#insufficient deep sequence data for at least one partner of couple 	10
	#couple most likely not a pair 							10
	#couple ambiguous if pair or not pair 					2 
	#couple most likely a pair direction not resolved 		7 
	#couple most likely a pair with resolved direction 		13 
	 
	rps	<- subset(rca, SELECT=='couple most likely a pair with resolved direction' & COUP_SC%in%c('M->F','F->M') & TYPE%in%c('mf','fm'))
	rps[, table(TYPE, COUP_SC)]
	#TYPE F->M M->F
  	#fm    5    1
  	#mf    0    7
	binconf(1,13)
	#PointEst       Lower     Upper
	#0.07692308 0.003945638 0.3331395
	subset(rca, !is.na(SDC_TYPE))
	#incorrect J056208    A078484     212
	#Argh, I mentioned this before: zero branch lengths red/blue should come out as intermingled, not ancestral.
	#talk to Matthew again
	#this error is entirely avoidable, and the current assignment is not optimal

	#	
	#	plot epilines
	if(0)
	{	
		t.posneg	<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
		setnames(t.posneg, c('BIRTHDATE','EST_DATEDIED'), c('DOB','DOD'))
		t.seq		<- unique(subset(rs, !is.na(PID), select=c(RID, PID, SEQ_DATE)))
		t.cd4		<- unique(subset(rd, !is.na(RECENTCD4DATE) & !is.na(RECENTCD4), select=c(RID, RECENTCD4DATE, RECENTCD4)))
		set(t.cd4, NULL, 'RECENTCD4', t.cd4[, cut(RECENTCD4, breaks=c(-1,250,350,500,800,Inf), labels=c('<200','200-349','350-499','500-799','800+'))])
		t.vl		<- unique(subset(rd, !is.na(RECENTVLDATE) & !is.na(RECENTVL), select=c(RID, RECENTVLDATE, RECENTVL)))
		set(t.vl, NULL, 'RECENTVL', t.vl[, cut(RECENTVL, breaks=c(-1,200,1e3,1e5,Inf), labels=c('<200','200-1,000','1,000-100,000','100,000+'))])
		#	plot 3 timelines per page
		setkey(rps, TYPE, POSTERIOR_SCORE)		
		p		<- RakaiFull.plot.epitimeline(subset(rps, SDC_TYPE=='correct'), copy(t.posneg), copy(t.cd4), copy(t.vl), copy(t.seq), age.adult=14)				
		pdf(file=paste0(outfile.base,'correctdirection_epilines.pdf'), w=8, h=8)
		print(p)
		dev.off()
		p		<- RakaiFull.plot.epitimeline(subset(rps, SDC_TYPE=='incorrect'), copy(t.posneg), copy(t.cd4), copy(t.vl), copy(t.seq), age.adult=14, xlim=c(2000,2014.5))				
		pdf(file=paste0(outfile.base,'incorrectdirection_epilines.pdf'), w=8, h=2)
		print(p)
		dev.off()
		
		
		pdf.n	<- 2
		pi	<- data.table(IDX=seq_along(p))
		pi[, PLOT:= ceiling(IDX/pdf.n)]
		pi[, PLOT_IDX:= (IDX-1)%%pdf.n+1]	
		pdf(file=paste0(outfile.base,'epilines.pdf'), w=20, h=12)	
		for(plot in pi[, unique(PLOT)])
		{
			idx			<- subset(pi, PLOT==plot)[, IDX]
			plot.idx	<- subset(pi, PLOT==plot)[, PLOT_IDX]
			grid.newpage()
			pushViewport(viewport(layout=grid.layout(1, pdf.n)))
			for(i in seq_along(idx))
				print(p[[idx[i]]], vp = viewport(layout.pos.row=1, layout.pos.col=plot.idx[i]))
		}
		dev.off()
		rtpd[, DUMMY:=NULL]
		
	}


	if(0)
	{
		#	plot incorrect J056208    A078484     212
		rps			<- subset(rtpd, MALE_RID=='J056208', select=c(MALE_RID, FEMALE_RID, PTY_RUN, TYPE, POSTERIOR_SCORE))
		rps[, DUMMY:=seq_len(nrow(rps))]
		rps[, LABEL:=rps[, factor(DUMMY, levels=DUMMY, labels=paste0('m ',MALE_RID,' f ', FEMALE_RID,'\n',TYPE,' ',round(POSTERIOR_SCORE, d=3),'\n',PTY_RUN))]]
		rpw[, ID_R_MAX:= pmax(FEMALE_R, MALE_R)]
		rpw[, ID_R_MIN:= pmin(FEMALE_R, MALE_R)]
		
		group		<- 'TYPE_BASIC'
		#group		<- 'TYPE_PAIR_TODI'					
		rpw2		<- subset(rpw, GROUP==group)
		rplkl2		<- subset(rplkl, GROUP==group)	
		plot.file	<- paste0(outfile.base,'incorrectdirection_windows_summary_',group,'.pdf')
		setnames(rps, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rpw2, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rplkl2, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		phsc.plot.windowsummaries.for.pairs(rps, rpw2, rplkl2, plot.file, cols=NULL, group=group)				
	}
	if(0)
	{
		require(colorspace)
		#	plot incorrect J056208    A078484     212
		zz		<- rtpd[, which(MALE_RID%in%c('J056208'))]
		indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170301_w250_s20_p35_stagetwo_rerun'
		#for(ii in seq_len(nrow(rtpdm))[-1])
		for(ii in zz)
		{		
			
			# load dfr and phs
			load( file.path(indir, paste0('ptyr',rtpd[ii,PTY_RUN],'_trees.rda')) )
			# setup plotting
			ids			<- c(rtpd[ii, MALE_RID],rtpd[ii, FEMALE_RID])
			dfs			<- subset(dfr, select=c(W_FROM, W_TO, IDX))
			dfs[, MALE_RID:=ids[1]]
			dfs[, FEMALE_RID:=ids[2]]				
			dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', rtpd[ii, PTY_RUN], '\nwindow ', W_FROM,'-', W_TO, sep='')]]
			plot.file	<- paste0(outfile.base, 'incorrectdirection_', rtpd[ii, PTY_RUN],'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
			invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.blacklisted=FALSE, drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10, tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'))						
		}
	}
}


RakaiFull.analyze.couples.todi.170811.compare.FF<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	cuts	<- data.table(	ID=   seq(1,7),
			PHSC= c(0.01,0.015,0.02,0.025,0.03,0.035,0.04),
			RWGD= c(0.0185, 0.026, 0.032, 0.0375, 0.041, 0.045, 0.048),
			FT=   c(0.037, 0.053, 0.066, 0.077, 0.088, 0.098, 0.107))
	cuts	<- melt(cuts, measure.vars=c('PHSC','RWGD','FT'))
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_ff_170811_cl3_prior23_min10_"
	
	#	load dc
	load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170704_assemblystatus.rda')	
	# 	load rd rp
	#	load couples to search for in phyloscanner output
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_info.rda")
	#
	#	load demographic info on all individuals
	tmp		<- RakaiCirc.epi.get.info.170208()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#rn		<- tmp$rn
	ra		<- tmp$ra
	#set(rn, NULL, 'RID', rn[, as.character(RID)])
	#rn		<- merge(rn, subset(rd, select=c(RID, FIRSTPOSVIS, FIRSTPOSDATE)), by='RID', all.x=1)
	#tmp		<- rn[, which(is.na(FIRSTPOSDATE) & !is.na(RECENTVLDATE) & TIMESINCEVL==0)]	#this is dodgy
	#set(rn, tmp, 'FIRSTPOSDATE', rn[tmp, RECENTVLDATE])		
	#rd		<- rbind(rd, rn, use.names=TRUE, fill=TRUE)	#do not consider individuals in the neuro study that are not part of RCCS
	set(rd, NULL, c('PID','SID'), NULL)
	set(rd, NULL, 'SEX', rd[, as.character(SEX)])
	set(rd, NULL, 'RECENTVL', rd[, as.numeric(gsub('< 150','1',gsub('> ','',gsub('BD','',gsub(',','',as.character(RECENTVL))))))])
	set(rd, NULL, 'CAUSE_OF_DEATH', rd[, as.character(CAUSE_OF_DEATH)])
	#	fixup rd: 
	#	remove HIV reverters without sequence
	rd		<- subset(rd, !RID%in%c("C117824","C119303","E118889","K067249"))
	#	fixup complex serology
	set(rd, rd[, which(RID=='B106184')], 'FIRSTPOSDATE', rd[which(RID=='B106184'),DATE])
	set(rd, rd[, which(RID=='B106184')], c('LASTNEGVIS','LASTNEGDATE'), NA_real_)
	set(rd, rd[, which(RID=='B106184')], c('HIVPREV'), 1)
	set(rd, rd[, which(RID=='A008742')], 'FIRSTPOSDATE', rd[which(RID=='A008742'),DATE])
	set(rd, rd[, which(RID=='A008742')], c('HIVPREV'), 1)
	#	fixup rd: 
	#	missing first pos date
	rd		<- subset(rd, RID!='A038432')	#has missing firstposdate and not in PANGEA anyway
	rd		<- subset(rd, RID!='H013226')	#has missing firstposdate and not in PANGEA anyway
	rd		<- subset(rd, RID!='K008173')	#has missing firstposdate and not in PANGEA anyway
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
						if(any(!is.na(FIRSTPOSVIS)))
							ans[is.na(FIRSTPOSVIS) | (!is.na(FIRSTPOSVIS) & FIRSTPOSVIS!=min(FIRSTPOSVIS, na.rm=TRUE))]	<- TRUE							
						list(DUMMY=DUMMY, DELETE=ans)		
					}, by=c('RID')], by=c('RID','DUMMY'))
	tmp		<- subset(tmp, !DELETE)
	set(tmp, NULL,c('DUMMY','DELETE'), NULL)
	set(rd, NULL, c('BIRTHDATE','LASTNEGDATE','FIRSTPOSVIS','FIRSTPOSDATE','ARVSTARTDATE','EST_DATEDIED'), NULL)
	rd		<- merge(rd, tmp, by='RID')	
	tmp		<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
	stopifnot(!nrow(merge(subset(tmp[, length(BIRTHDATE), by='RID'], V1>1), tmp, by='RID')))	
	
	outfile.base	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_ff_170811_'
	infiles			<- data.table(F= c(	'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min10_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min50_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl25_prior23_min30_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl45_prior23_min30_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_s10_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_s40_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_d30_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_d100_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_d1000_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mf3_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mf4_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mf6_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_mt_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_prt_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_rg1_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_rg20_allwindows.rda',
										'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_zbl_allwindows.rda'
												))
	infiles[, CONFIDENCE_CUT:=0.66]
	infiles[, NEFF_CUT:=3L]
	tmp				<- data.table(	F='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min30_allwindows.rda', 
									CONFIDENCE_CUT=seq(0.5,0.9,0.1),
									NEFF_CUT=3L)
	infiles			<- rbind(infiles, tmp)	
	rff				<- infiles[, {
				#F<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_allwindows.rda'
				load(F)	
				confidence.cut	<- CONFIDENCE_CUT
				neff.cut		<- NEFF_CUT
				cat('\n',confidence.cut)
				rff				<- subset(rplkl, GROUP=='TYPE_PAIR_TODI2' & TYPE=='linked' & NEFF>=neff.cut & ID1_SEX=='F' & ID2_SEX=='F')
				rff[, POSTERIOR_SCORE:= (POSTERIOR_ALPHA-1) / (POSTERIOR_ALPHA+POSTERIOR_BETA-2)]
				rff				<- subset(rff, POSTERIOR_SCORE>=confidence.cut)
				#	eliminate duplicates: make RID < RID_2 for FF pairs
				tmp				<- subset(rff, ID1_SEX=='F' & ID2_SEX=='F' & ID1>ID2)
				setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))	
				rff				<- rbind(subset(rff, ID1_SEX=='F' & ID2_SEX=='F' & ID1<ID2), tmp)
				#	define if individual part of couple
				rff[, ID1_COUPLE:= 0L]
				set(rff, rff[, which(ID1%in%rp$FEMALE_RID)], 'ID1_COUPLE', 1L)	
				rff[, ID2_COUPLE:= 0L]
				set(rff, rff[, which(ID2%in%rp$FEMALE_RID)], 'ID2_COUPLE', 1L)		
				rff	<- subset(rff, ID1_COUPLE==1 | ID2_COUPLE==1)
				#	eliminate duplicates: only pty.run with max NEFF
				tmp				<- rff[, list(PTY_RUN=PTY_RUN[which.max(NEFF)]), by=c('ID1','ID2')]
				rff				<- merge(rff, tmp, by=c('ID1','ID2','PTY_RUN'))
				rff
			}, by=c('F','CONFIDENCE_CUT','NEFF_CUT')]
	save(rff, file=paste0(outfile.base,'compare_runs.rda'))
}


RakaiFull.analyze.couples.todi.170811.compare.DIRext<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	outfile.base<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_'
	infiles		<- data.table(F=c(	'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min50_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl25_prior23_min30_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl45_prior23_min30_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_d30_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_d100_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_d1000_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_s10_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_s40_withmetadata.rda',											
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_mf3_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_mf4_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_mf6_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_mt_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_prt_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_rg1_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_rg20_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_zbl_withmetadata.rda',									
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_pcut50_ncut3_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_pcut60_ncut3_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_pcut70_ncut3_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_pcut80_ncut3_withmetadata.rda',
									'~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_pcut90_ncut3_withmetadata.rda'											
									))	
	#
	#	overall linked unlinked etc
	#
	rca			<- infiles[, {
				cat('\nprocess',F)	
				load(F)
				setkey(rca, MALE_RID, FEMALE_RID)
				rca[, PAIRID:=seq_len(nrow(rca))]
				#	prepare extending serodiscordant couples	
				rca[, EXT_TYPE:=NA_character_]
				set(rca, rca[, which(FEMALE_LASTNEGDATE>=MALE_FIRSTPOSDATE)], 'EXT_TYPE', 'serodisc-mf')		
				set(rca, rca[, which(MALE_LASTNEGDATE>=FEMALE_FIRSTPOSDATE)], 'EXT_TYPE', 'serodisc-fm')	
				#	add extra couples in who one has CD4<400 and the other has CD4>800
				#	for male potential recipient - evaluate CD4 around time male first positive
				tmp	<- rca[, which(is.na(EXT_TYPE) &
										MALE_RECENTCD4>800 & abs(MALE_RECENTCD4DATE-MALE_FIRSTPOSDATE)<2 &
										FEMALE_RECENTCD4<400 & abs(FEMALE_RECENTCD4DATE-MALE_FIRSTPOSDATE)<2)]
				set(rca, tmp, 'EXT_TYPE', 'CD4disc-fm')
				#	for female potential recipient - evaluate CD4 around time female first positive
				tmp	<- rca[, which(is.na(EXT_TYPE) &
										FEMALE_RECENTCD4>800 & abs(FEMALE_RECENTCD4DATE-FEMALE_FIRSTPOSDATE)<2 &
										MALE_RECENTCD4<400 & abs(MALE_RECENTCD4DATE-FEMALE_FIRSTPOSDATE)<2)]
				set(rca, tmp, 'EXT_TYPE', 'CD4disc-mf')
				#	separate into EXT_TYPE and EXT_DIR
				set(rca, NULL, 'EXT_DIR', rca[, gsub('^([a-zA-Z0-9]+)-([a-zA-Z]+)$','\\2',EXT_TYPE)])
				set(rca, NULL, 'EXT_TYPE', rca[, gsub('^([a-zA-Z0-9]+)-([a-zA-Z]+)$','\\1',EXT_TYPE)])
				#	subset to couples for who we have extended type
				df	<- subset(rca, !is.na(EXT_TYPE))
				df[, EXT_EVAL:= NA_character_]
				tmp	<- df[, which(TYPE%in%c('mf','fm'))]
				set(df, tmp, 'EXT_EVAL', df[tmp, as.character(factor(TYPE==EXT_DIR, levels=c(TRUE,FALSE),labels=c('correct','incorrect')))])
				tmp	<- df[, which(is.na(EXT_EVAL))]
				set(df, tmp, 'EXT_EVAL', df[tmp, SELECT])
				set(df, NULL, 'EXT_EVAL', df[, factor(EXT_EVAL, levels=c("insufficient deep sequence data for at least one partner of couple","couple most likely not a pair", "couple ambiguous if pair or not pair", "couple most likely a pair direction not resolved", "correct", "incorrect"))])
				rca	<- merge(rca, subset(df, select=c(PAIRID, EXT_EVAL)), by='PAIRID', all.x=TRUE)				
				rca
			}, by='F']			
	set(rca, NULL, 'F', rca[, basename(F)])
	df	<- rca[, list(N=length(PAIRID)), by=c('F','SELECT')]
	df	<- merge(as.data.table(expand.grid(F=unique(df$F), SELECT=unique(df$SELECT))), df, by=c('F', 'SELECT'), all.x=TRUE)
	set(df, df[, which(is.na(N))],'N',0L)
	set(df, NULL, 'SELECT', df[, factor(SELECT, levels=c( 	"insufficient deep sequence data for at least one partner of couple",								
															"couple most likely not a pair",                                     
															"couple ambiguous if pair or not pair",
															"couple most likely a pair direction not resolved",
															"couple most likely a pair with resolved direction"                
														),
												labels=c(	"data_insufficient",
															"unlinked",
															"linked_ambiguous",
															"linked_nodir",
															"linked_yesdir"))])	
	df	<- dcast.data.table(df, F~SELECT, value.var='N')
	df[, data_insufficient_pc:= paste0(round(100*data_insufficient/486,d=1),'%')]
	df[, linked:= linked_nodir+linked_yesdir]
	df[, linked_pc:= paste0(round(100*linked/(linked+unlinked+linked_ambiguous),d=1),'%')]
	df[, linked_yesdir_pc:= paste0(round(100*linked_yesdir/linked,d=1),'%')]
	#
	#	add sources correct incorrect
	#
	tmp	<- subset(rca, !is.na(EXT_TYPE) & EXT_EVAL%in%c('correct','incorrect'))[, list(N=length(PAIRID)), by=c('F','EXT_EVAL')]
	tmp	<- merge(as.data.table(expand.grid(F=unique(tmp$F), EXT_EVAL=unique(tmp$EXT_EVAL))), tmp, by=c('F', 'EXT_EVAL'), all.x=TRUE)
	set(tmp, tmp[, which(is.na(N))],'N',0L)
	set(tmp, NULL, 'EXT_EVAL', tmp[, factor(EXT_EVAL, levels=c('correct','incorrect'), labels=c('src_correct','src_incorrect'))])
	tmp	<- dcast.data.table(tmp, F~EXT_EVAL, value.var='N')
	#tmp[, src_FDR:= paste0(round(100*src_incorrect/(src_correct+src_incorrect),d=1),'%')]	
	tmp		<- merge(tmp, tmp[, {
						z<- round(100*as.numeric(binconf(src_incorrect, src_correct+src_incorrect)), d=1)
						list(src_FDR=paste0(z[1],'% [',z[2],'%-',z[3],'%]'))
					}, by='F'], by='F')
	df	<- merge(df, tmp, by='F')
	set(df, NULL, 'F', df[, gsub('_withmetadata.rda','',gsub('todi_couples_','',F))])
	#
	#	add FF linkages
	#
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_ff_170811_compare_runs.rda")
	tmp		<- rff[, list(linked_ff=length(ID1)), by=c('F','CONFIDENCE_CUT','NEFF_CUT')]
	set(tmp, NULL, 'F', tmp[, gsub('_allwindows.rda','',gsub('todi_pairs_','',basename(F)))])
	tmp2	<- tmp[, which(CONFIDENCE_CUT!=0.66 | NEFF_CUT!=3)]
	set(tmp, tmp2, 'F', tmp[tmp2,paste0(F,'_pcut',100*CONFIDENCE_CUT,'_ncut',NEFF_CUT)])
	set(tmp, NULL, c('CONFIDENCE_CUT','NEFF_CUT'), NULL)
	df		<- merge(df, tmp, by='F')
	#df[, linked_FDR:= paste0(round(100*linked_ff/linked, d=1),'%')]
	df		<- merge(df, df[, {
				z<- round(100*as.numeric(binconf(linked_ff, linked)), d=1)
				list(linked_FDR=paste0(z[1],'% [',z[2],'%-',z[3],'%]'))
			}, by='F'], by='F')
	#
	#	add options
	#
	df[, CONFIDENCE_CUT:= 0.66]
	set(df, df[, which(grepl('pcut',F))], 'CONFIDENCE_CUT', df[grepl('pcut',F),as.numeric(gsub('.*_pcut([0-9]+)_.*','\\1',F))/100])
	df[, NEFF_CUT:= 3]
	set(df, df[, which(grepl('ncut',F))], 'NEFF_CUT', df[grepl('ncut',F),as.numeric(gsub('.*_ncut([0-9]+).*','\\1',F))])
	df[, MIN_READ:= 20]
	set(df, df[, which(grepl('min',F))], 'MIN_READ', df[grepl('min',F),as.numeric(gsub('.*_min([0-9]+).*','\\1',F))])
	df[, DIST_CLOSE:= 3.5]
	set(df, df[, which(grepl('cl',F))], 'DIST_CLOSE', df[grepl('cl',F),as.numeric(gsub('3','35',gsub('.*_cl([0-9]+).*','\\1',F)))/10])
	df[, SANKHOFF_K:= 20]
	set(df, df[, which(grepl('s',F))], 'SANKHOFF_K', df[grepl('s',F),as.numeric(gsub('.*_s([0-9]+).*','\\1',F))])
	df[, DOWNSAMPLE:= 50]
	set(df, df[, which(grepl('d',F))], 'DOWNSAMPLE', df[grepl('d',F),as.numeric(gsub('.*_d([0-9]+).*','\\1',F))])
	df[, ZBL:=0L]
	set(df, df[, which(grepl('zbl',F))], 'ZBL', 1L)
	df[, PRT:=0L]
	set(df, df[, which(grepl('prt',F))], 'PRT', 1L)	
	df[, MT:=0L]
	set(df, df[, which(grepl('_mt',F))], 'MT', 1L)
	df[, ROGUE_MIN:= 10]
	set(df, df[, which(grepl('rg',F))], 'ROGUE_MIN', df[grepl('rg',F),as.numeric(gsub('.*_rg([0-9]+).*','\\1',F))])
	df[, MULTIFURC:= 10^(-5)]
	set(df, df[, which(grepl('mf',F))], 'MULTIFURC', df[grepl('mf',F),10^(-as.numeric(gsub('.*_mf([0-9]+).*','\\1',F)))])	
	df[, CENTRAL:=0L]
	set(df, df[, which(DOWNSAMPLE==50 & SANKHOFF_K==20 & DIST_CLOSE==3.5 & MIN_READ==30 & NEFF_CUT==3 & CONFIDENCE_CUT==0.66 & ZBL==0 & PRT==0 & MT==0 & ROGUE_MIN==10 & MULTIFURC==10^(-5))], 'CENTRAL', 1L)	
	#
	save(df, file=paste0(outfile.base,'compare_linked_srces.rda'))
	#
	#	make supp table
	#
	#	rogue
	tmp	<- subset(df, 	ROGUE_MIN!=10 | CENTRAL==1)[order(ROGUE_MIN), ]
	tmp[, ROW:= 0+seq_len(nrow(tmp))]
	ans	<- copy(tmp)
	#	downsample
	tmp	<- subset(df, 	DOWNSAMPLE!=50 | CENTRAL==1)[order(DOWNSAMPLE), ]
	tmp[, ROW:= 10+seq_len(nrow(tmp))]
	ans	<- rbind(ans, tmp)
	#	minimum reads
	tmp	<- subset(df, 	MIN_READ!=30 | CENTRAL==1)[order(MIN_READ), ]
	tmp[, ROW:= 20+seq_len(nrow(tmp))]
	ans	<- rbind(ans, tmp)	
	#	prune blacklists before anything else
	tmp	<- subset(df, 	PRT==1 | CENTRAL==1)[order(PRT), ]
	tmp[, ROW:= 30+seq_len(nrow(tmp))]
	ans	<- rbind(ans, tmp)
	#	zero branch length
	tmp	<- subset(df, 	ZBL==1 | CENTRAL==1)[order(ZBL), ]
	tmp[, ROW:= 40+seq_len(nrow(tmp))]
	ans	<- rbind(ans, tmp)
	#	collapse into multifurcations
	tmp	<- subset(df, 	MULTIFURC!=10^(-5) | CENTRAL==1)[order(MULTIFURC), ]
	tmp[, ROW:= 50+seq_len(nrow(tmp))]
	ans	<- rbind(ans, tmp)
	#	Sankoff
	tmp	<- subset(df, 	SANKHOFF_K!=20 | CENTRAL==1)[order(SANKHOFF_K), ]
	tmp[, ROW:= 60+seq_len(nrow(tmp))]
	ans	<- rbind(ans, tmp)
	#	allow multitrans
	tmp	<- subset(df, 	MT==1 | CENTRAL==1)[order(MT), ]
	tmp[, ROW:= 70+seq_len(nrow(tmp))]
	ans	<- rbind(ans, tmp)
	#	close
	tmp	<- subset(df, 	DIST_CLOSE!=3.5 | CENTRAL==1)[order(DIST_CLOSE), ]
	tmp[, ROW:= 80+seq_len(nrow(tmp))]
	ans	<- rbind(ans, tmp)
	#	confidence cut
	tmp	<- subset(df, 	CONFIDENCE_CUT!=0.66 | CENTRAL==1)[order(CONFIDENCE_CUT), ]
	tmp[, ROW:= 90+seq_len(nrow(tmp))]	
	ans	<- rbind(ans, tmp)	
	ans[, linked_label:= paste0(linked, ' (',linked_pc,')')]
	ans[, data_insufficient_label:= paste0(data_insufficient, ' (',data_insufficient_pc,')')]
	ans[, linked_yesdir_label:= paste0(linked_yesdir, ' (',linked_yesdir_pc,')')]
	ans[, type:= paste(ROGUE_MIN, DOWNSAMPLE, MIN_READ, PRT, ZBL, MULTIFURC, SANKHOFF_K, MT, DIST_CLOSE, CONFIDENCE_CUT)]
	ans	<- subset(ans, select=c(type, data_insufficient_label, linked_label, linked_FDR, linked_yesdir_label, src_FDR))
	write.csv(ans, file=paste0(outfile.base,'compare_linked_srces.csv'))
}


RakaiFull.analyze.couples.todi.170811.why.no.direction<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_withmetadata.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_"	
	
	load(infile)
	setkey(rca, MALE_RID, FEMALE_RID)
	rca[, PAIRID:=seq_len(nrow(rca))]
	rca		<- subset(rca, SELECT=='couple most likely a pair direction not resolved')
	rplkl2	<- subset(rplkl, GROUP=='TYPE_NETWORK_SCORES')
	rplkl2	<- merge(subset(rca, select=c(MALE_RID,FEMALE_RID,PTY_RUN)), rplkl2, by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	
	da		<- rplkl2[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('MALE_RID','FEMALE_RID','PTY_RUN')]
	da[, table(TYPE_MLE)]
	#	ambiguous                     fm                     mf not close/disconnected 
	#		   19                     17                     20                      3 
	subset(da, TYPE_MLE=='not close/disconnected')
	#	these have intermediates
		
	rplkl2	<- subset(rplkl, GROUP=='TYPE_BASIC')
	rplkl2	<- merge(subset(rca, select=c(MALE_RID,FEMALE_RID,PTY_RUN)), rplkl2, by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	rplkl2[, TYPE_TO:= NA_character_]
	set(rplkl2, rplkl2[,which(grepl('close', TYPE))], 'TYPE_TO', 'close other')
	set(rplkl2, rplkl2[,which(grepl('close', TYPE) & grepl('no intermediate', TYPE) & grepl('chain mf', TYPE))], 'TYPE_TO', 'mf')
	set(rplkl2, rplkl2[,which(grepl('close', TYPE) & grepl('no intermediate', TYPE) & grepl('chain fm', TYPE))], 'TYPE_TO', 'fm')
	set(rplkl2, rplkl2[,which(grepl('close', TYPE) & grepl('no intermediate', TYPE) & grepl('intermingled', TYPE))], 'TYPE_TO', 'intermingled')
	set(rplkl2, rplkl2[,which(grepl('close', TYPE) & grepl('no intermediate', TYPE) & grepl('other', TYPE))], 'TYPE_TO', 'adjacent')
	rplkl2	<- subset(rplkl2, !is.na(TYPE_TO))
	rplkl2	<- rplkl2[, list(KEFF=sum(KEFF)), by=c('MALE_RID','FEMALE_RID','PTY_RUN','TYPE_TO')]
	da		<- rplkl2[, list(TYPE_MLE=TYPE_TO[which.max(KEFF)]), by=c('MALE_RID','FEMALE_RID','PTY_RUN')]
	da[, table(TYPE_MLE)]
	# TYPE_MLE
    #	adjacent  close other           fm intermingled           mf 
    #   	   6            1           19            6           27 
	set(da, NULL, 'TYPE_MLE', da[, gsub('fm|mf', 'fm or mf',TYPE_MLE)])
	set(da, NULL, 'TYPE_MLE', da[, gsub('close other', 'intermingled',TYPE_MLE)])
	
	#	for all, compute median subtree distances
	#	compute confidence intervals
	rpw2	<- subset(rpw, GROUP=='TYPE_PAIR_TODI2')
	rpw2	<- merge(subset(rca, select=c(MALE_RID,FEMALE_RID,PTY_RUN)), rpw2, by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	da		<- merge(da, rpw2[, list(PHSC_Q50=median(PATRISTIC_DISTANCE)), by=c('MALE_RID','FEMALE_RID','PTY_RUN')], by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
	tmp		<- da[, list(N=length(MALE_RID), P=paste0('p',c(0.025,0.25,0.5,0.75,0.975)), Q=quantile(PHSC_Q50, p=c(0.025,0.25,0.5,0.75,0.975))), by='TYPE_MLE']
	tmp		<- dcast.data.table(tmp, TYPE_MLE+N~P, value.var='Q')
	tmp[, PHSC_Q5050:= paste0(round(p0.5*1e2, d=2), '% [',round(p0.25*1e2, d=2),'%-',round(p0.75*1e2, d=2),'%]')]
	write.csv( subset(tmp, select=c(TYPE_MLE, N, PHSC_Q5050)), file=paste0(outfile.base, 'unresolveddir.csv'))
	
}

RakaiFull.analyze.couples.todi.170811.sources.bias.from.read.coverage<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	indir			<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811'
	outfile.base	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_"
	
	infiles	<- data.table(F=list.files(indir, pattern='rda$', full.names=TRUE))
	infiles	<- subset(infiles, grepl('cl3_prior23_min30_d[0-9]+_withmetadata',F) | grepl('cl3_prior23_min30_withmetadata',F))
	infiles[, DOWNS:= gsub('.*_d([0-9]+)_.*','\\1',F)]
	set(infiles, infiles[, which(grepl('rda',DOWNS))],'DOWNS','50')
	set(infiles, NULL, 'DOWNS', infiles[,factor(DOWNS,levels=c('1000','100','50','30'),labels=c('none','100','50','30'))])
	
	ds		<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_d100_withmetadata.rda'
				load(F)
				ds	<- subset(rca, grepl('with resolved direction', SELECT), c(MALE_RID,FEMALE_RID,PTY_RUN,TYPE,KEFF,NEFF))
				setnames(ds, 'TYPE','TYPE_DIR')
				ds	<- merge(subset(rpw, GROUP=='TYPE_NETWORK_SCORES'), ds, by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
				ds	<- subset(ds, TYPE%in%c('fm','mf'))	
				tmp	<- ds[, which(TYPE=='mf')]
				set(ds,tmp,'READ_RATIO_R',ds[tmp,MALE_R/FEMALE_R])
				set(ds,tmp,'READ_RATIO_U',ds[tmp,MALE_L/FEMALE_L])
				tmp	<- ds[, which(TYPE=='fm')]
				set(ds,tmp,'READ_RATIO_R',ds[tmp,FEMALE_R/MALE_R])
				set(ds,tmp,'READ_RATIO_U',ds[tmp,FEMALE_L/MALE_L])
				ds
			}, by=c('DOWNS')]	
	#ds[,list(READ_RATIO_MEDIAN=median(READ_RATIO_U)), by='DOWNS']	
	ggplot(ds, aes(x=DOWNS, y=READ_RATIO_U)) + 
			geom_boxplot() + 
			theme_bw() + 
			scale_y_log10(breaks=c(0.02,0.1,0.2,0.33,0.5,0.66,1,1.5,2,3,5,10,50)) +
			coord_flip() +
			labs(	x='number of reads per individual\nin NGS phylogeny\nafter downsampling\n', 
					y='\n# unique reads in inferred source /\n# unqiue reads in inferred recipient')		
	ggsave(file=paste0(outfile.base,'downsampling_direction_among_coupleswithdirection.pdf'), w=12, h=3.5)
	
	ds		<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_d100_withmetadata.rda'
				load(F)
				ds	<- subset(rpw, GROUP=='TYPE_NETWORK_SCORES' & TYPE%in%c('fm','mf'))					
				tmp	<- ds[, which(TYPE=='mf')]
				set(ds,tmp,'READ_RATIO_R',ds[tmp,MALE_R/FEMALE_R])
				set(ds,tmp,'READ_RATIO_U',ds[tmp,MALE_L/FEMALE_L])
				tmp	<- ds[, which(TYPE=='fm')]
				set(ds,tmp,'READ_RATIO_R',ds[tmp,FEMALE_R/MALE_R])
				set(ds,tmp,'READ_RATIO_U',ds[tmp,FEMALE_L/MALE_L])
				ds
			}, by=c('DOWNS')]	
	ds[,list(READ_RATIO_MEDIAN=mean(READ_RATIO_U)), by='DOWNS']
	ggplot(ds, aes(x=DOWNS, y=READ_RATIO_U)) + 
			geom_boxplot() + 
			theme_bw() + 
			scale_y_log10(breaks=c(0.02,0.1,0.2,0.33,0.5,0.66,1,1.5,2,3,5,10,50)) +
			coord_flip() +
			labs(	x='number of reads per individual\nin NGS phylogeny\nafter downsampling\n', 
					y='\n# unique reads in inferred source /\n# unqiue reads in inferred recipient')		
	ggsave(file=paste0(outfile.base,'downsampling_direction_among_allcouples.pdf'), w=12, h=3.5)
	
}

RakaiFull.analyze.couples.todi.170811.siblings.not.contiguous<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_withmetadata.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_min30_"	
	
	#
	#	check adjacent for sibling --> OK
	#	
	load(infile)
	setkey(rca, MALE_RID, FEMALE_RID)
	rca[, PAIRID:=seq_len(nrow(rca))]	
	tmp		<- subset(rpw, PATHS_FM==0 & PATHS_MF==0 & ADJACENT & !CONTIGUOUS & GROUP=='TYPE_BASIC')
	rtpdm	<- unique(subset(tmp, select=c(MALE_RID, FEMALE_RID, PTY_RUN)))
	
	if(0)
	{
		require(colorspace)
		#for(ii in seq_len(nrow(rtpdm))[-1])
		for(ii in 1:4)
		{		
			if(rtpdm[ii, PTY_RUN]!=1)
			{
				indir		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30'
				# load dfr and phs
				load( file.path(indir, paste0('ptyr',rtpdm[ii,PTY_RUN],'_trees.rda')) )
				# setup plotting
				ids			<- c(rtpdm[ii, MALE_RID],rtpdm[ii, FEMALE_RID])
				dfs			<- subset(dfr, select=c(W_FROM, W_TO, IDX))
				dfs[, MALE_RID:=ids[1]]
				dfs[, FEMALE_RID:=ids[2]]
				dfs			<- merge(dfs, subset(rpw, GROUP=='TYPE_RAW'), by=c('MALE_RID','FEMALE_RID','W_FROM','W_TO'), all.x=1)	
				dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', PTY_RUN, '\nwindow ', W_FROM,'-', W_TO,'\npmf',PATHS_MF,' pfm',PATHS_FM, ' ',ADJACENT,' ',CONTIGUOUS,' ',TYPE, '\n', round(PATRISTIC_DISTANCE, d=5), sep='')]]
				plot.file	<- paste0(outfile.base, 'trees/adjnotcont_170811_run_', rtpdm[ii, PTY_RUN],'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
				invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10, tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'))			
			}				
		}
	}
	if(0)
	{
		require(colorspace)
		#for(ii in seq_len(nrow(rtpdm))[-1])
		for(ii in 2)
		{		
			if(rtpdm[ii, PTY_RUN]!=1)
			{
				indir		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min30'
				# load dfr and phs
				load( file.path(indir, paste0('ptyr',rtpdm[ii,PTY_RUN],'_trees.rda')) )
				# setup plotting
				ids			<- c(rtpdm[ii, MALE_RID],rtpdm[ii, FEMALE_RID])
				dfs			<- subset(dfr, select=c(W_FROM, W_TO, IDX))
				dfs[, MALE_RID:=ids[1]]
				dfs[, FEMALE_RID:=ids[2]]
				dfs			<- merge(dfs, subset(rpw, GROUP=='TYPE_RAW'), by=c('MALE_RID','FEMALE_RID','W_FROM','W_TO'), all.x=1)	
				dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', PTY_RUN, '\nwindow ', W_FROM,'-', W_TO,'\npmf',PATHS_MF,' pfm',PATHS_FM, ' ',ADJACENT,' ',CONTIGUOUS,' ',TYPE, '\n', round(PATRISTIC_DISTANCE, d=5), sep='')]]
				plot.file	<- paste0(outfile.base, 'trees/adjnotcont_nodrop_170811_run_', rtpdm[ii, PTY_RUN],'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
				invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.blacklisted=FALSE, drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10, tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'))			
			}				
		}
		subset(rpw, MALE_RID=='A106044' & PATHS_FM==0 & PATHS_MF>0 & ADJACENT & !CONTIGUOUS & GROUP=='TYPE_BASIC')
	}	 
	#
	#	check contiguous for ancestral 
	#	
	tmp		<- subset(rpw, PATHS_FM==0 & PATHS_MF>0 & ADJACENT & !CONTIGUOUS & GROUP=='TYPE_BASIC')
	subset(tmp, MALE_RID=='A106044')
	rtpdm	<- unique(subset(tmp, select=c(MALE_RID, FEMALE_RID, PTY_RUN)))

}

RakaiFull.analyze.couples.todi.170811.DIRext<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170522/todi_couples_170522_withmetadata.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170522/todi_couples_170522_"
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_withmetadata.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_"		
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10_withmetadata.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10_"	
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10_zbl_withmetadata.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10_zbl_"	
	infile					<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_withmetadata.rda'		
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_"	
	
	load(infile)
	setkey(rca, MALE_RID, FEMALE_RID)
	rca[, PAIRID:=seq_len(nrow(rca))]
	#
	#	prepare extending serodiscordant couples	
	rca[, EXT_TYPE:=NA_character_]
	set(rca, rca[, which(FEMALE_LASTNEGDATE>=MALE_FIRSTPOSDATE)], 'EXT_TYPE', 'serodisc-mf')		
	set(rca, rca[, which(MALE_LASTNEGDATE>=FEMALE_FIRSTPOSDATE)], 'EXT_TYPE', 'serodisc-fm')	
	#	add extra couples in who one has CD4<400 and the other has CD4>800
	#	for male potential recipient - evaluate CD4 around time male first positive
	tmp	<- rca[, which(is.na(EXT_TYPE) &
							MALE_RECENTCD4>800 & abs(MALE_RECENTCD4DATE-MALE_FIRSTPOSDATE)<2 &
							FEMALE_RECENTCD4<400 & abs(FEMALE_RECENTCD4DATE-MALE_FIRSTPOSDATE)<2)]
	set(rca, tmp, 'EXT_TYPE', 'CD4disc-fm')
	#	for female potential recipient - evaluate CD4 around time female first positive
	tmp	<- rca[, which(is.na(EXT_TYPE) &
							FEMALE_RECENTCD4>800 & abs(FEMALE_RECENTCD4DATE-FEMALE_FIRSTPOSDATE)<2 &
							MALE_RECENTCD4<400 & abs(MALE_RECENTCD4DATE-FEMALE_FIRSTPOSDATE)<2)]
	set(rca, tmp, 'EXT_TYPE', 'CD4disc-mf')
	#	add extra couples in who one has VL>1e5 and the other has VL<5e3
	#	for male potential recipient - evaluate VL around time male first positive
	#tmp	<- rca[, which(is.na(EXT_TYPE) & 
	#				MALE_RECENTVL>1e5 & abs(MALE_RECENTVLDATE-MALE_FIRSTPOSDATE)<1.5 & 
	#				FEMALE_RECENTVL<5e4 & abs(FEMALE_RECENTVLDATE-MALE_FIRSTPOSDATE)<1.5)]
	#set(rca, tmp, 'EXT_TYPE', 'VLdisc-fm')
	#	for female potential recipient - evaluate VL around time female first positive
	#tmp	<- rca[, which(is.na(EXT_TYPE) & 
	#						FEMALE_RECENTVL>1e5 & abs(FEMALE_RECENTVLDATE-FEMALE_FIRSTPOSDATE)<1.5 & 
	#						MALE_RECENTVL<5e4 & abs(MALE_RECENTVLDATE-FEMALE_FIRSTPOSDATE)<1.5)]
	#set(rca, tmp, 'EXT_TYPE', 'VLdisc-mf')
	#	separate into EXT_TYPE and EXT_DIR
	set(rca, NULL, 'EXT_DIR', rca[, gsub('^([a-zA-Z0-9]+)-([a-zA-Z]+)$','\\2',EXT_TYPE)])
	set(rca, NULL, 'EXT_TYPE', rca[, gsub('^([a-zA-Z0-9]+)-([a-zA-Z]+)$','\\1',EXT_TYPE)])
	#
	#	subset to couples for who we have extended type
	rca	<- subset(rca, !is.na(EXT_TYPE))
	rca[, EXT_EVAL:= NA_character_]
	tmp	<- rca[, which(TYPE%in%c('mf','fm'))]
	set(rca, tmp, 'EXT_EVAL', rca[tmp, as.character(factor(TYPE==EXT_DIR, levels=c(TRUE,FALSE),labels=c('correct','incorrect')))])
	tmp	<- rca[, which(is.na(EXT_EVAL))]
	set(rca, tmp, 'EXT_EVAL', rca[tmp, SELECT])
	set(rca, NULL, 'EXT_EVAL', rca[, factor(EXT_EVAL, levels=c("insufficient deep sequence data for at least one partner of couple","couple most likely not a pair", "couple ambiguous if pair or not pair", "couple most likely a pair direction not resolved", "correct", "incorrect"))])
	#
	rca[, table(EXT_EVAL, EXT_TYPE, useNA='if')]
	#  	insufficient deep sequence data for at least one partner of couple      21       10
  	#	couple most likely not a pair                                           14       10
  	#	couple ambiguous if pair or not pair                                     2        2
  	#	couple most likely a pair direction not resolved                         5        7
  	#	correct                                                                 10       12
  	#	incorrect                                                                1        1
	
	#	
	#	plot epilines
	if(0)
	{	
		plot.file	<- paste0(outfile.base,'testdirection_summary.pdf')
		t.posneg	<- unique(subset(rd, select=c(RID, BIRTHDATE, LASTNEGDATE, FIRSTPOSDATE, ARVSTARTDATE, EST_DATEDIED)))
		setnames(t.posneg, c('BIRTHDATE','EST_DATEDIED'), c('DOB','DOD'))
		t.seq		<- unique(subset(rs, !is.na(PID), select=c(RID, PID, SEQ_DATE)))
		t.cd4		<- unique(subset(rd, !is.na(RECENTCD4DATE) & !is.na(RECENTCD4), select=c(RID, RECENTCD4DATE, RECENTCD4)))		
		t.vl		<- unique(subset(rd, !is.na(RECENTVLDATE) & !is.na(RECENTVL), select=c(RID, RECENTVLDATE, RECENTVL)))
		set(t.vl, NULL, 'RECENTVL', t.vl[, cut(RECENTVL, breaks=c(-1,200,1e3,1e5,Inf), labels=c('<200','200-1,000','1,000-100,000','100,000+'))])
		#		
		
		df		<- subset(rca, EXT_EVAL=='correct'|EXT_EVAL=='incorrect', c(PAIRID, PTY_RUN, MALE_RID, FEMALE_RID, POSTERIOR_ALPHA, POSTERIOR_BETA, POSTERIOR_SCORE, EXT_EVAL, TYPE))
		df		<- df[order(TYPE, -POSTERIOR_SCORE), ]		
		set(df, NULL, 'PAIRID', df[, factor(PAIRID, levels=df$PAIRID)])
		rplkl2	<- merge(subset(df, select=c(PAIRID,MALE_RID,FEMALE_RID,PTY_RUN)), subset(rplkl, GROUP=='TYPE_DIR_TODI2'), by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
		rplkl3	<- merge(subset(df, select=c(PAIRID,MALE_RID,FEMALE_RID,PTY_RUN)), subset(rplkl, GROUP=='TYPE_NETWORK_SCORES'), by=c('MALE_RID','FEMALE_RID','PTY_RUN'))
		
		#	prepare t.posneg for plotting
		tpm	<- merge(t.posneg, data.table(RID=df$MALE_RID, SEX='M', PAIRID=df$PAIRID), by='RID')
		tpf	<- merge(t.posneg, data.table(RID=df$FEMALE_RID, SEX='F', PAIRID=df$PAIRID), by='RID')
		tmp	<- rbind( subset(tpm, select=c(PAIRID, FIRSTPOSDATE)), subset(tpf, select=c(PAIRID, FIRSTPOSDATE)) )
		cc	<- tmp[, list(FIRSTCONCPOS= max(FIRSTPOSDATE)), by='PAIRID']
		tpm	<- merge(tpm, cc, by='PAIRID')
		tpf	<- merge(tpf, cc, by='PAIRID')
		tpf[, DATE:= FIRSTPOSDATE-FIRSTCONCPOS]
		tpm[, DATE:= FIRSTPOSDATE-FIRSTCONCPOS]
		tpf[, NDATE:= LASTNEGDATE-FIRSTCONCPOS]
		tpm[, NDATE:= LASTNEGDATE-FIRSTCONCPOS]
		#	prepare t.seq for plotting
		tsm	<- merge(t.seq, data.table(RID=df$MALE_RID, SEX='M', PAIRID=df$PAIRID), by='RID')
		tsf	<- merge(t.seq, data.table(RID=df$FEMALE_RID, SEX='F', PAIRID=df$PAIRID), by='RID')
		tsm	<- merge(tsm, cc, by='PAIRID')
		tsf	<- merge(tsf, cc, by='PAIRID')
		tsf[, DATE:= SEQ_DATE-FIRSTCONCPOS]
		tsm[, DATE:= SEQ_DATE-FIRSTCONCPOS]
		#	prepare t.cd4 for plotting
		tcm	<- merge(t.cd4, data.table(RID=df$MALE_RID, SEX='M', PAIRID=df$PAIRID), by='RID')
		tcf	<- merge(t.cd4, data.table(RID=df$FEMALE_RID, SEX='F', PAIRID=df$PAIRID), by='RID')
		tcm	<- merge(tcm, cc, by='PAIRID')
		tcf	<- merge(tcf, cc, by='PAIRID')
		tcf[, DATE:= RECENTCD4DATE-FIRSTCONCPOS]
		tcm[, DATE:= RECENTCD4DATE-FIRSTCONCPOS]
		set(tcf, NULL, 'RECENTCD4_C', tcf[, cut(RECENTCD4, breaks=c(-1,400,800,Inf), labels=c('<400','400-799','800+'))])	
		set(tcm, NULL, 'RECENTCD4_C', tcm[, cut(RECENTCD4, breaks=c(-1,400,800,Inf), labels=c('<400','400-799','800+'))])
		#	prepare male/female
		sm	<- data.table(SEX='M', PAIRID=df$PAIRID, DATE= -1+min(min(tpm$NDATE, na.rm=TRUE), min(tpf$NDATE, na.rm=TRUE)))
		sf	<- data.table(SEX='F', PAIRID=df$PAIRID, DATE= -1+min(min(tpm$NDATE, na.rm=TRUE), min(tpf$NDATE, na.rm=TRUE)))
		
		plot.min.date	<- sm$DATE[1]
		plot.max.date	<- max(c(tcm$DATE,tcf$DATE,tsm$DATE,tsf$DATE,tpm$DATE,tpf$DATE))
		plot.nudge		<- 0.2
		plot.nudge.seq	<- 0.1
		plot.nudge.cd4	<- 0.2
		#	epilines
		p2	<- ggplot(df) + 				
				geom_text(data=sm, aes(x=DATE, y=PAIRID, label='M'), position=position_nudge(x=0, y=plot.nudge), size=2.5) +
				geom_text(data=sf, aes(x=DATE, y=PAIRID, label='F'), position=position_nudge(x=0, y=-plot.nudge), size=2.5) +				
				geom_point(data=subset(tpm, !is.na(NDATE)), aes(x=NDATE, y=PAIRID), position=position_nudge(x=0, y=plot.nudge), pch=23, size=5, fill='white') +
				geom_point(data=subset(tpf, !is.na(NDATE)), aes(x=NDATE, y=PAIRID), position=position_nudge(x=0, y=-plot.nudge), pch=23, size=5, fill='white') +
				geom_text(data=subset(tpm, !is.na(NDATE)), aes(x=NDATE, y=PAIRID, label='N'), position=position_nudge(x=0, y=plot.nudge), size=3) +
				geom_text(data=subset(tpf, !is.na(NDATE)), aes(x=NDATE, y=PAIRID, label='N'), position=position_nudge(x=0, y=-plot.nudge), size=3) +
				geom_point(data=subset(tpm, !is.na(DATE)), aes(x=DATE, y=PAIRID), position=position_nudge(x=0, y=plot.nudge), pch=23, size=5, fill='black') +
				geom_point(data=subset(tpf, !is.na(DATE)), aes(x=DATE, y=PAIRID), position=position_nudge(x=0, y=-plot.nudge), pch=23, size=5, fill='black') +
				geom_text(data=subset(tpm, !is.na(DATE)), aes(x=DATE, y=PAIRID, label='P'), position=position_nudge(x=0, y=plot.nudge), size=3, colour='white') +
				geom_text(data=subset(tpf, !is.na(DATE)), aes(x=DATE, y=PAIRID, label='P'), position=position_nudge(x=0, y=-plot.nudge), size=3, colour='white') +				
				geom_point(data=subset(tcm, !is.na(DATE)), aes(x=DATE, y=PAIRID, fill=RECENTCD4_C), position=position_nudge(x=0, y=plot.nudge.cd4), pch=23, size=3.5, colour='transparent') +
				geom_point(data=subset(tcf, !is.na(DATE)), aes(x=DATE, y=PAIRID, fill=RECENTCD4_C), position=position_nudge(x=0, y=-plot.nudge.cd4), pch=23, size=3.5, colour='transparent') +
				geom_text(data=subset(tsm, !is.na(DATE)), aes(x=DATE, y=PAIRID, label='S'), position=position_nudge(x=0, y=plot.nudge), size=2) +
				geom_text(data=subset(tsf, !is.na(DATE)), aes(x=DATE, y=PAIRID, label='S'), position=position_nudge(x=0, y=-plot.nudge), size=2) +				
				#geom_point(data=subset(tsm, !is.na(DATE)), aes(x=DATE, y=PAIRID), position=position_nudge(x=0, y=plot.nudge.seq), pch=83, size=3, colour='black') +
				#geom_point(data=subset(tsf, !is.na(DATE)), aes(x=DATE, y=PAIRID), position=position_nudge(x=0, y=-plot.nudge.seq), pch=83, size=3, colour='black') +
				scale_x_continuous(breaks=seq(-20,20,5), minor_breaks=seq(-20,20,1), limit=c(plot.min.date-.2,plot.max.date+.2), expand=c(0,0)) +
				scale_fill_manual(values=c('<400'=brewer.pal(9, 'YlOrRd')[6],'400-799'=brewer.pal(9, 'YlOrRd')[4],'800+'=brewer.pal(9, 'YlGn')[6])) +
				labs(x='', y='', fill='') +
				theme_bw() +
				theme(panel.grid.minor.x=element_line(colour="transparent", size=0.25), panel.grid.major.x=element_line(colour="white", size=0.5, linetype='dotted')) +
				theme(panel.grid.major.y=element_line(colour="grey85", size=18)) +
				theme(axis.text.y=element_blank()) +
				theme(axis.ticks.y=element_blank(), panel.border=element_blank())
		# correct yes/no
		p1	<- ggplot(df, aes(x=1, y=PAIRID, fill=EXT_EVAL)) +
				geom_tile() +
				geom_text(aes(label=factor(EXT_EVAL, levels=c('incorrect','correct'), labels=c('no','yes')), colour=EXT_EVAL)) +
				scale_x_continuous(expand=c(0,0)) +
				#scale_y_discrete(expand=c(0,0)) +
				scale_fill_manual(values=c('correct'='white','incorrect'='black')) +
				scale_colour_manual(values=c('incorrect'='white','correct'='black')) +
				labs(x='', y='') + guides(fill='none', colour='none') +
				theme_bw() + 
				theme(axis.text.x=element_text(colour="transparent"), axis.ticks.x=element_line(colour="transparent")) +				
				theme(axis.ticks.y=element_blank(), panel.border=element_blank())
		# number windows	
		p3	<- ggplot(rplkl2, aes(x=PAIRID, y=KEFF, fill=TYPE)) + 
				geom_bar(stat='identity',position='dodge',width=0.9) +
				scale_y_continuous(expand=c(0,0)) +
				scale_fill_manual(values=c('mf'='steelblue2','fm'='hotpink2','ambiguous'='grey50')) +			
				theme_bw() + 
				coord_flip() +
				labs(x='',y='',fill='') +
				theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) +
				theme(axis.text.y=element_blank()) +
				theme(axis.ticks.y=element_blank(), panel.border=element_blank())
		# posterior		
		p4	<- ggplot(rplkl2, aes(x=PAIRID,fill=TYPE,middle= qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA), lower=qbeta(0.25, POSTERIOR_ALPHA, POSTERIOR_BETA), upper=qbeta(0.75, POSTERIOR_ALPHA, POSTERIOR_BETA),ymin=qbeta(0.025, POSTERIOR_ALPHA, POSTERIOR_BETA),ymax=qbeta(0.975, POSTERIOR_ALPHA, POSTERIOR_BETA))) + 
				geom_boxplot(stat='identity',position='dodge') +
				scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.2), labels=scales::percent) +
				scale_fill_manual(values=c('mf'='steelblue2','fm'='hotpink2','ambiguous'='grey50')) +		
				theme_bw() + 
				theme(axis.text.y=element_blank()) +
				theme(axis.ticks.y=element_blank(), panel.border=element_blank()) +
				labs(x='',y='',fill='') +	
				coord_flip()
		pdf(file=plot.file, w=8.27*1.5, h=11.69*0.8)
		grid.newpage()	
		pushViewport(viewport(layout=grid.layout(1, 4, widths=unit(c(1.3,6,3.5,2.5), "null"))))   	
		print(p1, vp=viewport(layout.pos.row=1, layout.pos.col=1))
		print(p2, vp=viewport(layout.pos.row=1, layout.pos.col=2))
		print(p3, vp=viewport(layout.pos.row=1, layout.pos.col=3))
		print(p4, vp=viewport(layout.pos.row=1, layout.pos.col=4))
		#grid.draw(p3)
		dev.off()			
	}		
	if(0)
	{
		#	plot incorrect 
		rps			<- subset(rca, EXT_EVAL=='incorrect',select=c(MALE_RID, FEMALE_RID, PTY_RUN, TYPE, POSTERIOR_SCORE))
		#rps			<- subset(rtpd, MALE_RID=='J056208', select=c(MALE_RID, FEMALE_RID, PTY_RUN, TYPE, POSTERIOR_SCORE))
		rps[, DUMMY:=seq_len(nrow(rps))]
		rps[, LABEL:=rps[, factor(DUMMY, levels=DUMMY, labels=paste0('m ',MALE_RID,' f ', FEMALE_RID,'\n',TYPE,' ',round(POSTERIOR_SCORE, d=3),'\n',PTY_RUN))]]
		rpw[, ID_R_MAX:= pmax(FEMALE_R, MALE_R)]
		rpw[, ID_R_MIN:= pmin(FEMALE_R, MALE_R)]
		
		group		<- 'TYPE_BASIC'
		#group		<- 'TYPE_PAIR_TODI'					
		rpw2		<- subset(rpw, GROUP==group)
		rplkl2		<- subset(rplkl, GROUP==group)	
		plot.file	<- paste0(outfile.base,'incorrectdirection_windows_summary_',group,'.pdf')
		setnames(rps, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rpw2, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setnames(rplkl2, c('MALE_RID','FEMALE_RID'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
		phsc.plot.windowsummaries.for.pairs(rps, rpw2, rplkl2, plot.file, cols=NULL, group=group)				
	}
	if(0)
	{
		require(colorspace)
		#	plot incorrect J056208    A078484     212
		tmp		<- subset(rca, EXT_EVAL=='incorrect')[, MALE_RID]
		zz		<- rtpd[, which(MALE_RID%in%tmp)]
		indir	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min10'
		#for(ii in seq_len(nrow(rtpdm))[-1])
		for(ii in zz)
		{		
			
			# load dfr and phs
			load( file.path(indir, paste0('ptyr',rtpd[ii,PTY_RUN],'_trees.rda')) )
			# setup plotting
			ids			<- c(rtpd[ii, MALE_RID],rtpd[ii, FEMALE_RID])
			dfs			<- subset(dfr, select=c(W_FROM, W_TO, IDX))
			dfs[, MALE_RID:=ids[1]]
			dfs[, FEMALE_RID:=ids[2]]				
			dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', rtpd[ii, PTY_RUN], '\nwindow ', W_FROM,'-', W_TO, sep='')]]
			plot.file	<- paste0(outfile.base, 'incorrectdirection_', rtpd[ii, PTY_RUN],'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
			invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.blacklisted=FALSE, drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10, tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'))						
		}
	}
}

RakaiFull.analyze.couples.todi.170522.distance.histogram<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
		
	infile			<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_withmetadata.rda'		
	outfile.base	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min30_"		
	load(infile)
	dfd				<- subset(rpw, GROUP=='TYPE_BASIC')
	dfd				<- dfd[, list(	PHSC_PD_MEAN=mean(PATRISTIC_DISTANCE, na.rm=TRUE),
									PHSC_PD_Q025=quantile(PATRISTIC_DISTANCE, p=0.025, na.rm=TRUE),
									PHSC_PD_Q25=quantile(PATRISTIC_DISTANCE, p=0.25, na.rm=TRUE),
									PHSC_PD_Q75=quantile(PATRISTIC_DISTANCE, p=0.75, na.rm=TRUE),
									PHSC_PD_Q975=quantile(PATRISTIC_DISTANCE, p=0.975, na.rm=TRUE)	
									), by=c('PTY_RUN','MALE_RID','FEMALE_RID')]
	
	#
	require(mclust)
	dfds	<- subset(dfd, PHSC_PD_MEAN>1e-5 & PHSC_PD_MEAN<1)
	m1		<- densityMclust( log(dfds$PHSC_PD_MEAN), G=2)
	m2		<- densityMclust( log(dfds$PHSC_PD_MEAN))	#best fit really is 2D
	pdf(file=paste0(outfile.base,'_distances_consPatristic_lognormalmixture_pdf.pdf'), w=5, h=5)
	plot(m1, what='density', data= log(dfds$PHSC_PD_MEAN), breaks=50)
	dev.off()
	pdf(file=paste0(outfile.base,'_distances_consPatristic_lognormalmixture_cdf.pdf'), w=5, h=5)
	densityMclust.diagnostic(m1, type='cdf')	
	dev.off()
	#	mean and variance of first component	
	tmp		<- log(seq(1e-4,1,0.0001))
	th1		<- unname(c( summary(m1, parameters=TRUE)$mean, summary(m1, parameters=TRUE)$variance, summary(m1, parameters=TRUE)$pro)[c(1,3,5)])
	dens1	<- data.table(X=tmp, Y=dnorm(tmp, mean=th1[1], sd=sqrt(th1[2])))
	th2		<- unname(c( summary(m1, parameters=TRUE)$mean, summary(m1, parameters=TRUE)$variance, summary(m1, parameters=TRUE)$pro)[c(2,4,6)])
	dens2	<- data.table(X=tmp, Y=dnorm(tmp, mean=th2[1], sd=sqrt(th2[2])))
	densm	<- data.table(X=tmp, Y=predict.densityMclust(m1, tmp))
	
	#	3.5% threshold corresponds to 9.14% quantile of the LogNormal
	pnorm(log(0.035), mean=th1[1], sd=sqrt(th1[2]), lower.tail=FALSE)
	#	3.5% threshold corresponds to 0.6% quantile of the LogNormal
	pnorm(log(0.035), mean=th2[1], sd=sqrt(th2[2]), lower.tail=TRUE)	
	#	8% threshold corresponds to 6.936269e-07 quantile of the LogNormal
	pnorm(log(0.08), mean=th2[1], sd=sqrt(th2[2]), lower.tail=TRUE)
	
	#	10% quantile is 3.3% threshold
	exp(qnorm(0.9, mean=th1[1], sd=sqrt(th1[2])))
	#tmp		<- seq(1e-4,0.1,0.0001)
	#tmp		<- data.table(X=tmp, Y=dlnorm(tmp, meanlog=th1[1], sdlog=sqrt(th1[2]), log=FALSE))
	
	ggplot(dfds) +
			annotate("rect", xmin=log(2e-4), xmax=log(0.035), ymin=-Inf, ymax=Inf, fill=brewer.pal(11, 'PuOr')[2], alpha=0.5) +
			annotate("rect", xmin=log(0.08), xmax=log(1), ymin=-Inf, ymax=Inf, fill=rev(brewer.pal(11, 'RdGy'))[4], alpha=0.5) +
			#geom_vline(xintercept=log(c(0.035,0.08)), colour='grey70') +
			#geom_vline(xintercept=qnorm(c(1e-3, 0.005, 0.01), mean=th2[1], sd=sqrt(th2[2])), colour='blue') +
			#geom_vline(xintercept=qnorm(c(0.8,0.9,0.95), mean=th1[1], sd=sqrt(th1[2])), colour='red') +
			geom_histogram(aes(x=log(PHSC_PD_MEAN), y=..density..), binwidth=0.2, colour='white', fill='skyblue3') +
			geom_line(data=densm, aes(x=X, y=Y), lwd=0.8, colour='black') +			
			#geom_line(data=dens1, aes(x=X, y=Y*0.57), lwd=1.25, colour=brewer.pal(11, 'PuOr')[2]) +
			#geom_line(data=dens2, aes(x=X, y=Y*0.432), lwd=1.25, colour=rev(brewer.pal(11, 'RdGy'))[4]) +			
			#scale_colour_brewer(palette='Set1') +			
			scale_y_continuous(expand=c(0,0), limits=c(0,0.52)) + 
			scale_x_continuous(limits=log(c(2e-4, 1)), expand=c(0,0), 
					breaks=log(c(0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5)), 
					labels=paste0(100*c(0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5),'%')) +
			theme_bw() + theme(legend.position='bottom') +			
			labs(x='\naverage subtree distance between spouses\n( subst/site )', y='density', colour='')
	ggsave(file=paste0(outfile.base,'_distances_consPatristic_lognormalfitted.pdf'), w=5, h=2.25)
	
}


RakaiFull.analyze.couples.todi.170811.consensus.distances.read.from.FASTA<- function()
{	
	require(Phyloscanner.R.utilities)
	#	load latest sequence data
	load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170704_assemblystatus.rda')
	#	load couples to search for in phyloscanner output
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_info.rda")
	#	get all sequence comparisons between couples
	rps		<- unique(subset(rp, select=c(FEMALE_RID, MALE_RID, COUPID)))	
	tmp		<- unique(subset(dc, !is.na(SID), c(RID, PIDF)))
	setnames(tmp, c('RID','PIDF'), c('MALE_RID','MALE_PIDF'))
	rps		<- merge(rps, tmp, by='MALE_RID')
	setnames(tmp, c('MALE_RID','MALE_PIDF'), c('FEMALE_RID','FEMALE_PIDF'))
	rps		<- merge(rps, tmp, by='FEMALE_RID')
	
	#	get patristic distances
	if(0)
	{
		indir	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_bootstrap_trees'
		indir	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_couples_bootstrap_trees'
		infiles	<- data.table(F=list.files(indir, full.names=TRUE, pattern='^PANGEA.*newick'))
		infiles[, REP:= as.numeric(gsub('.*_ft\\.([0-9]+)\\.newick','\\1',F))]		
	}
	if(1)
	{
		indir	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_couples_ExaMLbootstrap_trees'
		infiles	<- data.table(F=list.files(indir, full.names=TRUE, pattern='^ExaML_result.PANGEA.*'))
		infiles[, REP:= as.numeric(gsub('.*finaltree\\.([0-9]+)','\\1',F))]		
	}
	
	dfd		<- infiles[, {
				#F<- '/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_bootstrap_trees/PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_ft.000.newick'
				phf			<- read.tree(F)	
				tmp			<- cophenetic.phylo(phf)
				tmp			<- as.matrix(tmp)
				tmp			<- as.data.table(melt(tmp))
				setnames(tmp, c('Var1','Var2','value'), c('MALE_PIDF','FEMALE_PIDF','PD'))
				set(tmp, NULL, 'MALE_PIDF', tmp[, as.character(MALE_PIDF)])
				set(tmp, NULL, 'FEMALE_PIDF', tmp[, as.character(FEMALE_PIDF)])
				tmp			<- merge(rps, tmp, by=c('MALE_PIDF','FEMALE_PIDF'))
				tmp
			}, by='REP']
	#	save
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_FastTree_patristicdistances.rda'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_couples_FastTree_patristicdistances.rda'
	outfile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_couples_ExaML_patristicdistances.rda'
	save(rps, dfd, file=outfile)	
}

RakaiFull.analyze.couples.todi.170811.plot.consensus.tree.for.paper<- function()
{	
	require(data.table)
	require(ape)
	require(ggtree)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(splines)
	
	infile		<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_170811_cl3_prior23_min10.rda"
	infile.tree	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_ft_bs100.newick'
	outfile.base<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/examples_'
	load(infile)
	
	#	load tree and rename tips
	ph	<- read.tree(infile.tree)
	dph	<- data.table(TAXA=ph$tip.label, TIP_ID=seq_len(Ntip(ph)))
	dph[, PIDF:= gsub('_WTSI.*','',TAXA)]
	dph	<- merge(dph, unique(subset(rs, select=c(PIDF,RID))), all.x=TRUE, by='PIDF')
	tmp	<- dph[, which(is.na(RID))]
	set(dph, tmp, 'RID', dph[tmp,PIDF])
	setkey(dph, TIP_ID)
	ph$tip.label	<- dph[, RID]
	ph$node.label[ph$node.label=='']	<- "0"
	ph$node.label	<- as.numeric(ph$node.label)/100
	
	dph[, TIP_COL:='black']
	set(dph, dph[, which(RID%in%c('K107705','G107639','J116744','K107569'))], 'TIP_COL', 'DarkRed')
	set(dph, dph[, which(RID%in%c('D105505','B104540'))], 'TIP_COL', 'DarkRed')
	set(dph, dph[, which(RID%in%c('B114802','H115007'))], 'TIP_COL', 'DarkRed')
	
	#
	#plot tree direct transmission
	tmp	<- drop.tip(ph, setdiff(dph[, RID], c('PG14-UG000568-S01881','K107705','G107639','J116744','K107569')), subtree=TRUE) 
	tmp	<- phytools:::reroot(tmp, which(tmp$tip.label=='PG14-UG000568-S01881'), position=0)
	dt	<- data.table(TAXA=tmp$tip.label, IDX=seq_along(tmp$tip.label), COUNT=sub('\\[([0-9]+)_tip.*','\\1',tmp$tip.label))
	dt[, COLOUR:= as.character(factor(grepl('^[0-9]+$',COUNT), levels=c(TRUE,FALSE),labels=c('grey50','DarkRed')))]
	set(dt, dt[, which(!grepl('^[0-9]+$',COUNT))], 'COUNT', "1")
	set(dt, NULL, 'COUNT', dt[, as.integer(COUNT)])	
	ggtree(tmp) %<+% dt +
			geom_text2(aes(subset=!isTip, label=label), hjust=-.3, size=2, colour='grey50') +
			geom_tippoint(aes(size=COUNT), shape=18) +
			scale_size(range=c(1,20)) +
			#scale_fill_hue(na.value="black") +								
			theme(legend.position="none") +
			#scale_x_continuous() +
			geom_tiplab(align=T, aes(col=I(COLOUR)), size=3, linetype='dotted', linesize=NA) +
			guides(size='none') +
			theme_tree2() +
			theme(	axis.line.x=element_line(),
					panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3),
					legend.position="bottom") + 
			ggplot2::xlim(0, max(node.depth.edgelength(tmp)[1:Ntip(tmp)])*1.15) +
			labs(x='subst/site')
	ggsave(file=paste0(outfile.base,'tree_direct_transmission.pdf'), w=4, h=7)
	
	#
	#plot tree no resolution
	tmp	<- drop.tip(ph, setdiff(dph[, RID], c('PG14-UG000568-S01881','D105505','B104540','H096194','C107300')), subtree=TRUE) 
	tmp	<- phytools:::reroot(tmp, which(tmp$tip.label=='PG14-UG000568-S01881'), position=0)
	dt	<- data.table(TAXA=tmp$tip.label, IDX=seq_along(tmp$tip.label), COUNT=sub('\\[([0-9]+)_tip.*','\\1',tmp$tip.label))
	dt[, COLOUR:= as.character(factor(grepl('^[0-9]+$',COUNT), levels=c(TRUE,FALSE),labels=c('grey50','DarkRed')))]
	set(dt, dt[, which(!grepl('^[0-9]+$',COUNT))], 'COUNT', "1")
	set(dt, NULL, 'COUNT', dt[, as.integer(COUNT)])	
	ggtree(tmp) %<+% dt +
			geom_text2(aes(subset=!isTip, label=label), hjust=-.3, size=2, colour='grey50') +
			geom_tippoint(aes(size=COUNT), shape=18) +
			scale_size(range=c(1,20)) +
			#scale_fill_hue(na.value="black") +								
			theme(legend.position="none") +
			#scale_x_continuous() +
			geom_tiplab(align=T, aes(col=I(COLOUR)), size=3, linetype='dotted', linesize=NA) +
			guides(size='none') +
			theme_tree2() +
			theme(	axis.line.x=element_line(),
					panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3),
					legend.position="bottom") + 
			ggplot2::xlim(0, max(node.depth.edgelength(tmp)[1:Ntip(tmp)])*1.15) +
			labs(x='subst/site')
	ggsave(file=paste0(outfile.base,'tree_no_resolution.pdf'), w=4, h=7)
	
	pdf(file=paste0(infile.tree,'.pdf'), w=15, h=200)
	plot(ph, cex=0.2, edge.width=0.5, align.tip.label=FALSE, show.node.label=TRUE, no.margin=FALSE, label.offset=0.001, tip.color=dph[, TIP_COL])
	axisPhylo()
	dev.off()
	
	
	#
	#plot tree dual infection
	tmp	<- drop.tip(ph, setdiff(dph[, RID], c('PG14-UG000568-S01881','B114802','H115007','C114764')), subtree=TRUE) 
	tmp	<- phytools:::reroot(tmp, which(tmp$tip.label=='PG14-UG000568-S01881'), position=0)
	dt	<- data.table(TAXA=tmp$tip.label, IDX=seq_along(tmp$tip.label), COUNT=sub('\\[([0-9]+)_tip.*','\\1',tmp$tip.label))
	dt[, COLOUR:= as.character(factor(grepl('^[0-9]+$',COUNT), levels=c(TRUE,FALSE),labels=c('grey50','DarkRed')))]
	set(dt, dt[, which(!grepl('^[0-9]+$',COUNT))], 'COUNT', "1")
	set(dt, NULL, 'COUNT', dt[, as.integer(COUNT)])	
	ggtree(tmp) %<+% dt +
			#geom_text(aes(x=branch, label=edge.length)) +
			geom_text2(aes(subset=!isTip, label=label), hjust=-.3, size=2, colour='grey50') +
			geom_tippoint(aes(size=COUNT), shape=18) +
			scale_size(range=c(1,20)) +
			#scale_fill_hue(na.value="black") +								
			theme(legend.position="none") +
			#scale_x_continuous() +
			geom_tiplab(align=T, aes(col=I(COLOUR)), size=3, linetype='dotted', linesize=NA) +
			guides(size='none') +
			theme_tree2() +
			theme(	axis.line.x=element_line(),
					panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3),
					legend.position="bottom") + 
			ggplot2::xlim(0, max(node.depth.edgelength(tmp)[1:Ntip(tmp)])*1.15) +
			labs(x='subst/site')
	ggsave(file=paste0(outfile.base,'tree_dual_infection.pdf'), w=4, h=7)
	
	#
	#	plot entire tree
	pdf(file=paste0(infile.tree,'.pdf'), w=15, h=200)
	plot(ph, cex=0.2, edge.width=0.5, align.tip.label=FALSE, show.node.label=TRUE, no.margin=FALSE, label.offset=0.001, tip.color=dph[, TIP_COL])
	axisPhylo()
	dev.off()
	
	#
	#	divergence in dual infection pair
	tmp		<- subset(dph, RID%in%c('B114802','H115007'))[, TAXA]
	infiles	<- data.table(F=list.files('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_couples_bootstrap_trees', pattern='^PANGEA.*newick$', full.names=TRUE))
	ds		<- infiles[, {
				#F	<- '/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_couples_bootstrap_trees/PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_couples_ft.099.newick'
				tr	<- read.tree(F)
				stopifnot(length(setdiff(tmp,tr$tip.label))==0)
				list(PD=cophenetic(tr)[tmp[1],tmp[2]])
			}, by='F']
	ds[, quantile(PD, p=c(0.025,0.5,0.975))]
	#      2.5%        50%      97.5% 
	#0.07329975 0.08853500 0.10583525
}

RakaiFull.analyze.couples.todi.170811.compare.to.consensus<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	require(splines)
	#outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_"	
	#outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couplesExaMLcouples_170811_"
	#outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couplescouples_170811_"
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couplescouples_170811_min10_"
	#	
	#	load patristic distances 
	#	we only have this for a subset of couples with at least 700 nt long consensus
	#infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_FastTree_patristicdistances.rda'
	#infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_couples_ExaML_patristicdistances.rda'
	infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_couples_FastTree_patristicdistances.rda'
	load(infile)
	
	#
	#	load preprocessed couples and add TAXA + SIDs
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couples_170811_cl3_prior23_min10_withmetadata.rda"	
	load(infile)
	setkey(rca, MALE_RID, FEMALE_RID)
	rca[, PAIRID:=seq_len(nrow(rca))]
	
	rpw2	<- subset(rpw, GROUP=='TYPE_PAIR_DI2')
		
	#	load raw genetic distances 
	infile	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/consensus/PANGEA_HIV_Imperial_v170704_UG_bestcov_cov700_rawgendist.rda"
	load(infile)		
	gd		<- subset(gd, grepl('^PG',TAXA) & grepl('^PG',TAXA2))
	set(gd, NULL, 'TAXA', gd[, gsub('_WTSI.*','',TAXA)])
	set(gd, NULL, 'TAXA2', gd[, gsub('_WTSI.*','',TAXA2)])	
	setnames(gd, c('TAXA','TAXA2'), c('MALE_TAXA','FEMALE_TAXA'))
	tmp		<- copy(gd)
	setnames(tmp, c('MALE_TAXA','FEMALE_TAXA'), c('FEMALE_TAXA','MALE_TAXA'))
	gd		<- rbind(gd, tmp, use.names=TRUE)
	tmp		<- subset(rp, !is.na(MALE_TAXA) & !is.na(FEMALE_TAXA), select=c(MALE_RID,FEMALE_RID,MALE_TAXA,FEMALE_TAXA))
	gd		<- merge(gd, tmp, by=c('FEMALE_TAXA','MALE_TAXA'))
	gd		<- gd[, list(CONS_GDRW=mean(CONS_GDRW)), by=c('MALE_RID','FEMALE_RID')]
	
	#
	#	get 95% confidence intervals for consensus distances and phyloscanner distances
	#	
	tmp		<- rpw2[,	list(	PHSC_PD_MEAN=mean(PATRISTIC_DISTANCE, na.rm=TRUE),
					PHSC_PD_Q025=quantile(PATRISTIC_DISTANCE, p=0.025, na.rm=TRUE),
					PHSC_PD_Q25=quantile(PATRISTIC_DISTANCE, p=0.25, na.rm=TRUE),
					PHSC_PD_Q50=quantile(PATRISTIC_DISTANCE, p=0.5, na.rm=TRUE),
					PHSC_PD_Q75=quantile(PATRISTIC_DISTANCE, p=0.75, na.rm=TRUE),
					PHSC_PD_Q975=quantile(PATRISTIC_DISTANCE, p=0.975, na.rm=TRUE)	
			), by=c('MALE_RID','FEMALE_RID')]
	dfd2	<- dfd[, list(	CONS_PD_MEAN=mean(PD, na.rm=TRUE),
					CONS_PD_Q025=quantile(PD, p=0.025, na.rm=TRUE),
					CONS_PD_Q25=quantile(PD, p=0.25, na.rm=TRUE),
					CONS_PD_Q50=quantile(PD, p=0.5, na.rm=TRUE),
					CONS_PD_Q75=quantile(PD, p=0.75, na.rm=TRUE),
					CONS_PD_Q975=quantile(PD, p=0.975, na.rm=TRUE)	
			),	by=c('COUPID','MALE_RID','FEMALE_RID')]
	dfd2	<- merge(dfd2, tmp, by=c('MALE_RID','FEMALE_RID'))
	dfd2	<- merge(dfd2, gd, by=c('MALE_RID','FEMALE_RID'), all.x=TRUE)
	set(dfd2, NULL, 'PHSC_PD_ASYM', dfd2[, as.numeric(PHSC_PD_MEAN<PHSC_PD_Q25 | PHSC_PD_MEAN>PHSC_PD_Q75)])
	set(dfd2, NULL, 'CONS_PD_ASYM', dfd2[, as.numeric(CONS_PD_MEAN<CONS_PD_Q25 | CONS_PD_MEAN>CONS_PD_Q75)])
	
	#
	#	correlation between consensus raw genetic dist and phylogenetic dist	
	tmp		<- subset(dfd2, !is.na(CONS_GDRW))
	#	determine outliers in regression curve based on Cooks Distance	
	tmp[, LOG_CONS_GDRW:= log10(CONS_GDRW)]
	tmp[, LOG_CONS_PD_MEAN:= log10(CONS_PD_MEAN)]
	tmp[, LOG_CONS_PD_Q50:= log10(CONS_PD_Q50)]
	m1		<- lm(LOG_CONS_GDRW~ns(LOG_CONS_PD_Q50, df=2), data=tmp)
	#	df=2 spline looks good
	tmp[, PR:= 10^predict(m1,type='response')]	
	#	ggplot(tmp, aes(x=CONS_PD_MEAN)) + geom_point(size=1, aes(y=CONS_GDRW)) + geom_line(aes(y=PR)) + scale_x_log10(labels=percent, expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) + scale_y_log10(labels=percent, expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) + coord_cartesian(xlim=c(0.002, 0.75), ylim=c(0.002, 0.5))  
	tmp[, COOKSD:=cooks.distance(m1)]
	tmp[, COOKSD_OUTLIER:= as.integer(COOKSD>6*mean(COOKSD))]	
	tmp[, CI_WIDTH:= (log10(CONS_PD_Q75)-log10(CONS_PD_Q25))]
	tmp[, CONS_PD_HV:= as.integer(CI_WIDTH>quantile(CI_WIDTH, p=0.95))]
	ggplot(tmp, aes(x=CONS_PD_Q50, xmin=CONS_PD_Q25, xmax=CONS_PD_Q75)) +
			geom_abline(slope=1, intercept=0, colour='black', linetype='dotted') +
			geom_line(aes(y=PR)) +
			geom_errorbarh(aes(y=CONS_GDRW), size=.5, alpha=0.5, colour='grey40', height = 0) +
			geom_point(aes(y=CONS_GDRW, colour=as.character(CONS_PD_HV)), size=1) +			
			scale_x_log10(labels=percent, expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			scale_y_log10(labels=percent, expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			coord_cartesian(xlim=c(0.002, 0.75), ylim=c(0.002, 0.5)) + 
			theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=15)) +
			guides(colour=FALSE) +
			labs(	x='\npatristic distance between consensus sequences\n(subst/site)',
					y='patristic distance between WH clades\n(avg subst/site across windows)\n')
	ggsave(file=paste0(outfile.base,'_distances_correlation_patristic_rawgenetic.pdf'), w=7, h=7)
	#	flag patristic distance outliers
	dfd2	<- merge(dfd2, subset(tmp, select=c(MALE_RID, FEMALE_RID, COOKSD_OUTLIER, CONS_PD_HV)), by=c('MALE_RID','FEMALE_RID'), all.x=TRUE)
	set(dfd2, dfd2[, which(is.na(COOKSD_OUTLIER))], 'COOKSD_OUTLIER', 2L)
	set(dfd2, dfd2[, which(is.na(CONS_PD_HV))], 'CONS_PD_HV', 2L)
	#dfd2	<- subset(dfd2, COOKSD_OUTLIER==0)
	#	
	if(1)
	{
		#Rose ARHR: raw genetic distance of 4-5.3% optimal for detecting linkage
		zz		<- data.table(CONS_PD_Q50=seq(0.07,0.14,0.001), LOG_CONS_PD_Q50=log10(seq(0.07,0.14,0.001)))
		zz[, PR:= predict(m1,type='response',newdata=zz)]
		zz[, PR_CONS_GDRW:= 10^PR]		
	}


	#
	#	plot histograms to determine cut off points
	#
	tmp		<- melt(dfd2, id.vars=c('COUPID','MALE_RID','FEMALE_RID'), measure.vars=c('CONS_PD_Q50','PHSC_PD_Q50'))
	ggplot(tmp, aes(x=log10(value), colour=variable)) +
		geom_vline(xintercept=log10(c(0.035)), colour='blue', linetype='dotted') +
		geom_vline(xintercept=log10(c(0.08)), colour='red', linetype='dotted') +
		stat_ecdf() +
		scale_colour_brewer(palette='Set1') +			
		scale_y_continuous(expand=c(0,0), limits=c(0,1)) + 
		scale_x_continuous(		limits=log10(c(0.0009, 1)), expand=c(0,0), 
				breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5)), 
				labels=paste0(100*c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5),'%')) +
		theme_bw() + theme(legend.position='bottom') +			
		labs(x='\nsubst/site', y='', colour='')
	ggsave(file=paste0(outfile.base,'_distances_ecdf.pdf'), w=5, h=5)	
	#
	#	correlation plot
	#
	tmp3	<- copy(dfd2)		
	set(tmp3, tmp3[,which(PHSC_PD_Q50< 10^(-3.1))], 'PHSC_PD_Q50', 10^(-3.1))
	ggplot(tmp3, aes(x=log10(CONS_PD_Q50))) +			
			geom_linerange(aes(ymin=log10(PHSC_PD_Q25), ymax=log10(PHSC_PD_Q75)), size=.5, alpha=0.2, colour='black') +
			geom_errorbarh(aes(y=log10(PHSC_PD_Q50), xmin=log10(CONS_PD_Q25), xmax=log10(CONS_PD_Q75)), size=.5, alpha=0.2, colour='black', height = 0) +
			geom_point(size=1, aes(y=log10(PHSC_PD_Q50))) +			
			geom_abline(slope=1, intercept=0, colour='black', linetype='dotted') +				
			scale_x_continuous(labels=paste0(c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25),'%'), expand=c(0,0), breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25))) +
			scale_y_continuous(labels=paste0(c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25),'%'), expand=c(0,0), breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25))) +
			coord_cartesian(xlim=log10(c(min(tmp3$CONS_PD_Q50), 0.5)), ylim=log10(c(0.0007, 0.5))) +
			theme_bw() +
			guides(colour=FALSE) +
			theme(axis.title=element_text(size=15), axis.text=element_text(size=15)) +
			labs(	x='\npatristic distance between consensus sequences\n(subst/site)',
					y='subtree distance across genomic windows\n(subst/site)\n')
	ggsave(file=paste0(outfile.base,'_distances_consPatristic.pdf'), w=7, h=7)
		
	dfd2[, cor(log10(CONS_PD_MEAN), log10(PHSC_PD_MEAN))]
	#	0.932
	dfd2[, cor(log10(CONS_PD_Q50), log10(PHSC_PD_Q50))]
	#	0.757
	dfd2[, cor(log10(CONS_PD_MEAN), log10(PHSC_PD_MEAN), method='spearman')]
	#	0.91
	dfd2[, cor(log10(CONS_PD_Q50), log10(PHSC_PD_Q50), method='spearman')]	
	#	0.8951977


	#
	#	define outliers
	#	determine outliers as before through regression curve based on Cooks Distance	
	#	exclude from regression fit: highly uncertain points & major outliers
	#
	tmp	<- copy(dfd2)
	tmp[, LOG_PHSC_PD_Q50:= log10(PHSC_PD_Q50)]	
	tmp[, LOG_CONS_PD_Q50:= log10(CONS_PD_Q50)]
	tmp[, NEW_Y:=LOG_PHSC_PD_Q50-LOG_CONS_PD_Q50]
	if(1)
	{
		tmp2	<- subset(tmp, CONS_PD_HV==0 & PHSC_PD_Q50>=1e-4 & !(CONS_PD_Q50>0.07 & PHSC_PD_Q50<0.01))
		m1		<- lm(LOG_PHSC_PD_Q50~ns(LOG_CONS_PD_Q50, df=3), data=tmp2)
		tmp[, PR:= predict(m1,type='response',newdata=tmp)]
		tmp[, PRSD:= predict(m1,type='response',newdata=tmp, se.fit=TRUE)$se.fit]
		tmp[, PR_U:=PR+2*PRSD]
		tmp[, PR_L:=PR-2*PRSD]
		tmp2[, COOKSD:=cooks.distance(m1)]
		tmp		<- merge(tmp, subset(tmp2, select=c(MALE_RID, FEMALE_RID, COOKSD)), by=c('MALE_RID','FEMALE_RID'), all.x=TRUE)
		set(tmp, tmp[,which(is.na(COOKSD) & CONS_PD_HV==0)], 'COOKSD', 1)
		
		zz		<- data.table(CONS_PD_Q50=seq(0.05,0.14,0.001), LOG_CONS_PD_Q50=log10(seq(0.05,0.14,0.001)))
		zz[, PR:= predict(m1,type='response',newdata=zz)]
		zz[, PR_PHSC_PD_Q50:= 10^PR]
		
		zz		<- data.table(CONS_PD_Q50=seq(0.01,0.09,0.001), LOG_CONS_PD_Q50=log10(seq(0.01,0.09,0.001)))
		zz[, PR:= predict(m1,type='response',newdata=zz)]
		zz[, PR_PHSC_PD_Q50:= 10^PR]
	}	
	if(0)
	{		
		tmp2	<- subset(tmp, PHSC_PD_Q50>=1e-4 & !(CONS_PD_Q50>0.07 & PHSC_PD_Q50<0.01))
		m1		<- lm(LOG_PHSC_PD_Q50~LOG_CONS_PD_Q50, data=tmp2)
		tmp[, PR:= predict(m1,type='response',newdata=tmp)]		
		tmp[, PRSD:= predict(m1,type='response',newdata=tmp, se.fit=TRUE)$se.fit]
		tmp[, PR_U:=PR+2*PRSD]
		tmp[, PR_L:=PR-2*PRSD]
		tmp2[, COOKSD:=cooks.distance(m1)]
		tmp		<- merge(tmp, subset(tmp2, select=c(MALE_RID, FEMALE_RID, COOKSD)), by=c('MALE_RID','FEMALE_RID'), all.x=TRUE)
		set(tmp, tmp[,which(is.na(COOKSD))], 'COOKSD', 1)
	}
	if(0)
	{
		m1		<- lm(NEW_Y~1, data=tmp)	
		tmp[, PR:= 10^(LOG_CONS_PD_Q50+predict(m1,type='response'))]			
	}	
	if(0)
	{
		tmp2	<- subset(tmp, PHSC_PD_Q50>0.035)
		m1		<- lm(LOG_PHSC_PD_Q50~LOG_CONS_PD_Q50, data=tmp2)
		tmp[, PR:= 10^(predict(m1,type='response',newdata=tmp))]
	}
	if(0)
	{ 
		m1		<- lm(NEW_Y~1, data=subset(tmp, PHSC_PD_Q50>0.035))
		tmp[, PR:= 10^(LOG_CONS_PD_Q50+predict(m1,type='response',newdata=tmp))]
	}
	if(1)
	{
		#tmp[, OUT:= CONS_PD_Q50>0.07 & PHSC_PD_Q50<0.01]
		tmp3	<- copy(tmp)		
		set(tmp3, tmp3[,which(LOG_PHSC_PD_Q50< -3.1)], 'LOG_PHSC_PD_Q50', -3.1)		
		set(tmp3, NULL, 'CLASS', tmp3[,as.character(as.numeric(COOKSD>0.027/2))])
		set(tmp3, tmp3[,which(CONS_PD_HV>=1)], 'CLASS', '2')
		ggplot(tmp3, aes(x=LOG_CONS_PD_Q50)) +
				geom_ribbon(aes(ymin=PR_L, ymax=PR_U), fill='black',alpha=0.15) +
				geom_line(aes(y=PR)) +
				geom_point(size=1, aes(y=LOG_PHSC_PD_Q50, colour=CLASS)) +
				scale_colour_manual(values=c('1'='red','0'='black','2'='blue')) +
				geom_abline(slope=1, intercept=0, colour='black', linetype='dotted') +				
				scale_x_continuous(labels=paste0(c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25),'%'), expand=c(0,0), breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25))) +
				scale_y_continuous(labels=paste0(c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25),'%'), expand=c(0,0), breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25))) +
				coord_cartesian(xlim=c(min(tmp3$LOG_CONS_PD_Q50), log10(0.5)), ylim=log10(c(0.0007, 0.5))) +
				theme_bw() +
				guides(colour=FALSE) +
				theme(axis.title=element_text(size=15), axis.text=element_text(size=15)) +
				labs(	x='\npatristic distance between consensus sequences\n(subst/site)',
						y='subtree distance across genomic windows\n(subst/site)\n')
		ggsave(file=paste0(outfile.base,'_distances_consPatristic_outliers.pdf'), w=7, h=7)
		dfd2	<- merge(dfd2, subset(tmp3, select=c(MALE_RID,FEMALE_RID,CLASS, PR, PR_U, PR_L)), by=c('MALE_RID','FEMALE_RID'))
		set(dfd2, NULL, 'PR', dfd2[, 10^PR])
		set(dfd2, NULL, 'PR_L', dfd2[, 10^PR_L])
		set(dfd2, NULL, 'PR_U', dfd2[, 10^PR_U])
	}
	if(0)
	{
		require(BLR)
		data(wheat)     #Loads the wheat dataset
		y=Y[,1]
		### Creates a testing set with 100 observations
		whichNa<-sample(1:length(y),size=100,replace=FALSE)
		yNa<-y
		yNa[whichNa]<-NA
		### Runs the Gibbs sampler
		fm<-BLR(	y=yNa,
					XL=X,
					GF=list(ID=1:nrow(A),A=A),
					prior=list(varE=list(df=3,S=0.25),
					varU=list(df=3,S=0.63),
					lambda=list(shape=0.52,rate=1e-4,type='random',value=30)),
					nIter=1000,
					burnIn=500,
					thin=1,
					saveAt="example_")
		MSE.tst<-mean((fm$yHat[whichNa]-y[whichNa])^2)
		MSE.tst
		MSE.trn<-mean((fm$yHat[-whichNa]-y[-whichNa])^2)
		MSE.trn
		COR.tst<-cor(fm$yHat[whichNa],y[whichNa])
		COR.tstCOR.trn<-cor(fm$yHat[-whichNa],y[-whichNa])
		COR.trn
		plot(fm$yHat~y,xlab="Phenotype",ylab="Pred. Gen. Value" ,cex=.8)
		points(x=y[whichNa],y=fm$yHat[whichNa],col=2,cex=.8,pch=19)
		x11()
		plot(scan('example_varE.dat'),type="o",ylab=expression(paste(sigma[epsilon]^2)))
	}
	
	
	#
	#	look at outliers
	#
	#	read existing manual comments on previous classification scheme
	dfd3	<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couplescouples_170811_min10__distances_consPatristic_outliers_EDITED_OR.csv'))
	dfd3	<- subset(dfd3, select=c(MALE_RID,FEMALE_RID,ASSIGNMENT,PAIR.CHAIN,OTHER))
	dfd3	<- merge(dfd2, dfd3, by=c('MALE_RID','FEMALE_RID'),all.x=TRUE)
	dfd3	<- subset(dfd3, CLASS>=1 & is.na(ASSIGNMENT))
	
	#	check distribution of patristic distance estimates
	tmp		<- subset(dfd3, select=c(MALE_RID,FEMALE_RID,CONS_PD_Q50,PHSC_PD_Q50,CLASS, PR, PR_U, PR_L))
	dfd4	<- merge(tmp,subset(rpw2, select=c(MALE_RID,FEMALE_RID,PTY_RUN,W_FROM,PATRISTIC_DISTANCE)),by=c('MALE_RID','FEMALE_RID'))
	tmp		<- merge(tmp, subset(dfd, select=c(MALE_RID,FEMALE_RID,PD)),by=c('MALE_RID','FEMALE_RID'))
	setnames(tmp, 'PD','PATRISTIC_DISTANCE')
	dfd4	<- rbind(dfd4,tmp,fill=TRUE)
	dfd4[, TYPE:= dfd4[,factor(is.na(PTY_RUN),levels=c(TRUE,FALSE),labels=c('CONS','PHSC'))]]
	dfd4	<- merge(dfd4,unique(subset(rca, select=c(MALE_RID,FEMALE_RID,PAIRID))),by=c('MALE_RID','FEMALE_RID'))
	ggplot(dfd4) + 
			geom_point(aes(x=PATRISTIC_DISTANCE, y=PAIRID+as.numeric(TYPE=='CONS')-0.5, colour=TYPE)) +
			geom_errorbarh(aes(y=PAIRID, x=PR, xmin=PR_L, xmax=PR_U), colour='black') +
			theme_bw() +
			scale_x_log10(labels=paste0(c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25),'%'), expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			facet_wrap(~paste(PAIRID, MALE_RID, FEMALE_RID,sep=' '), ncol=1, scales='free_y')
	ggsave(file=paste0(outfile.base,'_distances_consPatristic_outliers_v2_distances.pdf'), w=20, h=20)
	
	#	inspect manually
	#dfco<- subset(merge(rca, dfd3, by=c('MALE_RID','FEMALE_RID')), select=c(MALE_RID, FEMALE_RID, PTY_RUN, SELECT, COUP_SC, CONS_GDRW, CONS_PD_MEAN, CONS_PD_Q50, PHSC_PD_Q50, CLASS, PR))
	#dfco<- subset(merge(rca, subset(dfd5, is.na(CLASS_MANUAL) & PHSC_PD_Q50<0.035), by=c('MALE_RID','FEMALE_RID')), select=c(MALE_RID, FEMALE_RID, PTY_RUN, SELECT, COUP_SC, CONS_GDRW, CONS_PD_MEAN, CONS_PD_Q50, PHSC_PD_Q50, CLASS, PR))
	#dfco<- subset(merge(rca, subset(dfd5, is.na(CLASS_MANUAL) & PHSC_PD_Q50>=0.035), by=c('MALE_RID','FEMALE_RID')), select=c(MALE_RID, FEMALE_RID, PTY_RUN, SELECT, COUP_SC, CONS_GDRW, CONS_PD_MEAN, CONS_PD_Q50, PHSC_PD_Q50, CLASS, PR))
	dfco<- subset(merge(rca, subset(dfd5, is.na(CLASS_MANUAL)), by=c('MALE_RID','FEMALE_RID')), select=c(MALE_RID, FEMALE_RID, PTY_RUN, SELECT, COUP_SC, CONS_GDRW, CONS_PD_MEAN, CONS_PD_Q50, PHSC_PD_Q50))
	dfco<- dfco[order(CONS_PD_Q50), ]	
	#write.csv(dfco, file=paste0(outfile.base,'_distances_consPatristic_outliers_v2.csv'))
	#write.csv(dfco, file=paste0(outfile.base,'_distances_consPatristic_nooutliersclose.csv'))
	#write.csv(dfco, file=paste0(outfile.base,'_distances_consPatristic_nooutliersnotclose.csv'))
	write.csv(dfco, file=paste0(outfile.base,'_distances_consPatristic_fromPTYRUN1.csv'))
	if(0)
	{
		require(colorspace)
		#zz	<- dfco[, which(MALE_RID%in%c('H115099'))]
		for(ii in seq_len(nrow(dfco)))
		#for(ii in zz)
		{		
			indir		<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_output_170704_w250_s20_p35_stagetwo_rerun23_min10'
			# load dfr and phs
			load( file.path(indir, paste0('ptyr',dfco[ii,PTY_RUN],'_trees.rda')) )
			# setup plotting
			ids			<- c(dfco[ii, as.character(MALE_RID)],dfco[ii, as.character(FEMALE_RID)])
			dfs			<- subset(dfr, select=c(W_FROM, W_TO, IDX))
			dfs[, MALE_RID:=ids[1]]
			dfs[, FEMALE_RID:=ids[2]]				
			dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', dfco[ii, PTY_RUN], '\nwindow ', W_FROM,'-', W_TO, sep='')]]
			plot.file	<- paste0(outfile.base, 'run_', dfco[ii, PTY_RUN],'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
			invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.blacklisted=FALSE, drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10, tip.regex='^(.*)_fq[0-9]+_read_([0-9]+)_count_([0-9]+)$'))						
		}
	}
	
	#	read existing manual comments on previous classification scheme
	dfd5	<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couplescouples_170811_min10__distances_consPatristic_outliers_EDITED_OR.csv'))
	dfd5	<- subset(dfd5, select=c(MALE_RID,FEMALE_RID,COMPLEX_TOPOLOGY,PAIR_CHAIN,OTHER))
	tmp		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couplescouples_170811_min10__distances_consPatristic_outliers_v2_EDITED_OR.csv'))
	tmp		<- subset(tmp, select=c(MALE_RID,FEMALE_RID,COMPLEX_TOPOLOGY,PAIR_CHAIN,OTHER))
	dfd5	<- rbind(dfd5,tmp)
	tmp		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couplescouples_170811_min10__distances_consPatristic_nooutliersclose_EDITED_OR.csv'))
	tmp		<- subset(tmp, select=c(MALE_RID,FEMALE_RID,COMPLEX_TOPOLOGY,PAIR_CHAIN,OTHER))
	dfd5	<- rbind(dfd5,tmp)
	tmp		<- as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170811/todi_couplescouples_170811_min10__distances_consPatristic_nooutliersnotclose_EDITED_OR.csv'))
	tmp		<- subset(tmp, select=c(MALE_RID,FEMALE_RID,COMPLEX_TOPOLOGY,PAIR_CHAIN,OTHER))
	dfd5	<- rbind(dfd5,tmp)
	dfd5	<- subset(dfd5, !is.na(COMPLEX_TOPOLOGY) & COMPLEX_TOPOLOGY!='')	
	dfd5[, CLASS_MANUAL:=NA_character_]
	set(dfd5, dfd5[, which(grepl('no clear signal',tolower(COMPLEX_TOPOLOGY)))],'CLASS_MANUAL','unclear')
	set(dfd5, dfd5[, which(grepl('nothing',tolower(COMPLEX_TOPOLOGY)))],'CLASS_MANUAL','typical topological configurations')
	set(dfd5, dfd5[, which(grepl('insufficient',tolower(COMPLEX_TOPOLOGY)))],'CLASS_MANUAL','unclear')
	set(dfd5, dfd5[, which(grepl('yes',tolower(COMPLEX_TOPOLOGY)))],'CLASS_MANUAL','multiple infection/dual infection/recombinants')
	set(dfd5, dfd5[, which(grepl('yes! ish',tolower(COMPLEX_TOPOLOGY)))],'CLASS_MANUAL','multiple infection/dual infection/recombinants, weak signal')
	set(dfd5, dfd5[, which(grepl('long single tip branches|siblings|long tip branch|tree artifact',tolower(COMPLEX_TOPOLOGY)))],'CLASS_MANUAL','unusual subtree distances')
	stopifnot( !nrow(subset(dfd5, is.na(CLASS_MANUAL))) )
	dfd5	<- subset(dfd5, select=c(MALE_RID,FEMALE_RID,CLASS_MANUAL))
	dfd5	<- merge(dfd2,dfd5,by=c('MALE_RID','FEMALE_RID'),all.x=TRUE)
	
	
	
	tmp3	<- copy(dfd5)		
	set(tmp3, tmp3[,which(log10(PHSC_PD_Q50)< -3.1)], 'PHSC_PD_Q50', 10^(-3.1))
	set(tmp3, NULL, 'CLASS', tmp3[, factor(CLASS, levels=c('0','1','2'), labels=c('broadly comparable','discordant','highly variable'))])
	set(tmp3, tmp3[, which(is.na(CLASS_MANUAL))], 'CLASS_MANUAL', 'not manually checked')
	ggplot(tmp3, aes(x=log10(CONS_PD_Q50))) +
			#geom_ribbon(aes(ymin=PR_L, ymax=PR_U), fill='black',alpha=0.15) +
			#geom_line(aes(y=PR)) +
			geom_point(size=1.5, aes(y=log10(PHSC_PD_Q50), shape=CLASS, colour=CLASS_MANUAL)) +
			scale_shape_manual(values=c('discordant'=15,'broadly comparable'=16,'highly variable'=17)) +
			scale_colour_manual(values=c(	"multiple infection/dual infection/recombinants"='orangered',              
											"multiple infection/dual infection/recombinants, weak signal"='plum', 
											"insufficient data"='grey50',                              
											"not manually checked"='black',
											"typical topological configurations"='black',
											"unclear"='bisque4',                                        
											"unusual subtree distances"='purple')) +
			geom_abline(slope=1, intercept=0, colour='black', linetype='dotted') +				
			scale_x_continuous(labels=paste0(c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25),'%'), expand=c(0,0), breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25))) +
			scale_y_continuous(labels=paste0(c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25),'%'), expand=c(0,0), breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25))) +
			coord_cartesian(xlim=log10(c(min(tmp3$CONS_PD_Q50), 0.5)), ylim=log10(c(0.0007, 0.5))) +
			theme_bw() +
			theme(axis.title=element_text(size=15), axis.text=element_text(size=15)) +
			labs(	x='\npatristic distance between consensus sequences\n(subst/site)',
					y='subtree distance across genomic windows\n(subst/site)\n',
					pch='subtree distances versus\nconsensus distances',
					colour='topological configurations')
	ggsave(file=paste0(outfile.base,'_distances_consPatristic_manualassessed_v2.pdf'), w=10, h=7)
	
	
	ggplot(subset(tmp3, CLASS_MANUAL!='not manually checked'), aes(x=log10(CONS_GDRW))) +
			#geom_ribbon(aes(ymin=PR_L, ymax=PR_U), fill='black',alpha=0.15) +
			#geom_line(aes(y=PR)) +
			geom_point(size=1.5, aes(y=log10(PHSC_PD_Q50), colour=CLASS_MANUAL)) +
			scale_shape_manual(values=c('discordant'=15,'broadly comparable'=16,'highly variable'=17)) +
			scale_colour_manual(values=c(	"multiple infection/dual infection/recombinants"='orangered',              
							"multiple infection/dual infection/recombinants, weak signal"='plum', 
							"insufficient data"='grey50',                              
							"not manually checked"='black',
							"typical topological configurations"='black',
							"unclear"='bisque4',                                        
							"unusual subtree distances"='purple')) +
			geom_abline(slope=1, intercept=0, colour='black', linetype='dotted') +				
			scale_x_continuous(labels=paste0(c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25),'%'), expand=c(0,0), breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25))) +
			scale_y_continuous(labels=paste0(c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25),'%'), expand=c(0,0), breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25))) +
			coord_cartesian(xlim=log10(c(min(tmp3$CONS_GDRW, na.rm=TRUE), 0.25)), ylim=log10(c(0.0007, 0.5))) +
			theme_bw() +
			theme(axis.title=element_text(size=15), axis.text=element_text(size=15)) +
			labs(	x='\ngenetic distance between consensus sequences\n(subst/site)',
					y='subtree distance across genomic windows\n(subst/site)\n',
					pch='subtree distances versus\nconsensus distances',
					colour='topological configurations')
	ggsave(file=paste0(outfile.base,'_distances_consRawGenetic_manualassessed_v2.pdf'), w=10, h=7)

	
	#
	#	define model between phs distance and raw genetic distance
	#
	tmp		<- subset(dfd5, CLASS_MANUAL=='typical topological configurations')	
	tmp[, LOG_PHSC_PD_Q50:= log10(PHSC_PD_Q50)]	
	tmp[, LOG_CONS_GDRW:= log10(CONS_GDRW)]		
	tmp2	<- subset(tmp, !is.na(CONS_GDRW) & CONS_PD_HV==0 & PHSC_PD_Q50>=1e-4 & !(CONS_PD_Q50>0.07 & PHSC_PD_Q50<0.01))
	#tmp2	<- subset(tmp, !is.na(CONS_GDRW) & CONS_PD_HV==0 & !(CONS_PD_Q50>0.07 & PHSC_PD_Q50<0.01))
	m1		<- lm(LOG_PHSC_PD_Q50~ns(LOG_CONS_GDRW, df=4), data=tmp2)
	tmp[, PR:= predict(m1,type='response',newdata=tmp)]
	tmp[, PRSD:= predict(m1,type='response',newdata=tmp, se.fit=TRUE)$se.fit]
	tmp[, PR_U:=PR+2*PRSD]
	tmp[, PR_L:=PR-2*PRSD]
	
	tmp3	<- subset(tmp, !is.na(CONS_GDRW))		
	set(tmp3, tmp3[,which(LOG_PHSC_PD_Q50< -3.1)], 'LOG_PHSC_PD_Q50', -3.1)				
	ggplot(tmp3, aes(x=LOG_CONS_GDRW)) +
			geom_ribbon(aes(ymin=PR_L, ymax=PR_U), fill='black',alpha=0.15) +
			geom_line(aes(y=PR)) +
			geom_point(size=1, aes(y=LOG_PHSC_PD_Q50)) +				
			geom_abline(slope=1, intercept=0, colour='black', linetype='dotted') +				
			scale_x_continuous(labels=paste0(c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25),'%'), expand=c(0,0), breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25))) +
			scale_y_continuous(labels=paste0(c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25),'%'), expand=c(0,0), breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25))) +
			coord_cartesian(xlim=c(min(tmp3$LOG_CONS_GDRW), log10(0.5)), ylim=log10(c(0.0007, 0.5))) +
			theme_bw() +				
			theme(axis.title=element_text(size=15), axis.text=element_text(size=15)) +
			labs(	x='\npatristic distance between consensus sequences\n(subst/site)',
					y='subtree distance across genomic windows\n(subst/site)\n')
	ggsave(file=paste0(outfile.base,'_distances_consRawGenetic_spline.pdf'), w=7, h=7)	
	zz		<- data.table(CONS_GDRW=seq(0.015,0.07,0.001), LOG_CONS_GDRW=log10(seq(0.015,0.07,0.001)))
	zz[, PR:= predict(m1,type='response',newdata=zz)]
	zz[, PR_PHSC_PD_Q50:= 10^PR]	
		
}

RakaiFull.analyze.couples.todi.170522.compare.to.consensus<- function()
{	
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
			
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170522/todi_couples_170522_"
	outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170522/todi_couples_170611_"
	#	
	#	load patristic distance matrix
	#	we only have this at present for a subset of couples -- never mind
	infile	<- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_UG_gag_coinfinput_160219.rda'
	load(infile)
	#	prepare patristic distance data.table
	ph.gdtr	<- as.data.table(melt(ph.gdtr, varnames=c('TAXA1','TAXA2')))
	setnames(ph.gdtr, 'value', 'PD')
	ph.gdtr	<- subset(ph.gdtr, TAXA1!=TAXA2)
	set(ph.gdtr, NULL, 'TAXA1', ph.gdtr[, gsub('_','-',as.character(TAXA1))])
	set(ph.gdtr, NULL, 'TAXA2', ph.gdtr[, gsub('_','-',as.character(TAXA2))])
	#	load genetic distance matrix with overlap
	infile		<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/circumcision/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_gd.rda'
	load(infile)	#loads sq.gd
	setnames(sq.gd, 'PD', 'GD')
	#	load PANGEA sequence info
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/circumcision/RCCS_SeqInfo_170505.rda")
	setnames(rs, 'SAMPLE_DATE', 'SEQ_DATE')
	
	#
	#	load couples and add TAXA + SIDs
	infile					<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170522/todi_couples_170522_withmetadata.rda"	
	load(infile)
	load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_stagetwo.rda')
	dfd			<- subset(rca, select=c(MALE_RID, FEMALE_RID, PTY_RUN))
	tmp2		<- subset(pty.runs, select=c('RID','SID','PTY_RUN'))
	setnames(tmp2, c('RID','SID'), c('MALE_RID','MALE_SID'))	
	dfd			<- merge(dfd, tmp2, by=c('MALE_RID','PTY_RUN'))
	setnames(tmp2, c('MALE_RID','MALE_SID'), c('FEMALE_RID','FEMALE_SID'))
	dfd			<- merge(dfd, tmp2, by=c('FEMALE_RID','PTY_RUN'))
	set(dfd, NULL, 'MALE_RID', dfd[, as.character(MALE_RID)])
	set(dfd, NULL, 'FEMALE_RID', dfd[, as.character(FEMALE_RID)])
	tmp2		<- subset(rs, select=c(RID, PIDF, SID))
	setnames(tmp2, colnames(tmp2), paste0('MALE_',colnames(tmp2)))
	dfd			<- merge(dfd, tmp2, by=c('MALE_RID','MALE_SID'))
	setnames(tmp2, colnames(tmp2), gsub('MALE_','FEMALE_',colnames(tmp2)))
	dfd			<- merge(dfd, tmp2, by=c('FEMALE_RID','FEMALE_SID'))
	
	#
	#	add patristic distances and genetic distances
	setnames(dfd, c('MALE_PIDF','FEMALE_PIDF'), c('TAXA1','TAXA2'))
	dfd			<- merge(dfd, sq.gd, by=c('TAXA1','TAXA2'), all.x=1)				#in fact sq.gd is symmetric so could shortcut this	
	dfd			<- merge(dfd, ph.gdtr, by=c('TAXA1','TAXA2'), all.x=1)	#in fact ph.gdtr is symmetric
	
	#	for each couple take average consensus if there are multiple combinations
	setkey(dfd, MALE_RID, FEMALE_RID, PTY_RUN)
	stopifnot( nrow(unique(dfd,by=c('MALE_RID','FEMALE_RID','PTY_RUN')))==nrow(unique(dfd,by=c('MALE_RID','FEMALE_RID'))) )
	dfd			<- dfd[, list(PD= mean(PD[which(!is.na(PD))]), GD= mean(GD[which(!is.na(GD))]) ), by=c('MALE_RID','FEMALE_RID')]
	
	#
	#	focus on mean distances + quantiles among couples
	#
	
	tmp			<- subset(rpw, GROUP=='TYPE_PAIR_DI')
	tmp			<- merge(unique(subset(rp, select=c(MALE_RID, FEMALE_RID))), tmp, by=c('MALE_RID','FEMALE_RID'))
	tmp2		<- subset(tmp, W_FROM>=800 & W_TO<=4650)[, list(	PHSC_PD_MEAN=mean(PATRISTIC_DISTANCE, na.rm=TRUE),
					PHSC_PD_Q025=quantile(PATRISTIC_DISTANCE, p=0.025, na.rm=TRUE),
					PHSC_PD_Q25=quantile(PATRISTIC_DISTANCE, p=0.25, na.rm=TRUE),
					PHSC_PD_Q75=quantile(PATRISTIC_DISTANCE, p=0.75, na.rm=TRUE),
					PHSC_PD_Q975=quantile(PATRISTIC_DISTANCE, p=0.975, na.rm=TRUE)	
			), by=c('PTY_RUN','MALE_RID','FEMALE_RID')]
	tmp2[, PHSC_W:='gag+pol']
	tmp			<- tmp[, list(	PHSC_PD_MEAN=mean(PATRISTIC_DISTANCE, na.rm=TRUE),
					PHSC_PD_Q025=quantile(PATRISTIC_DISTANCE, p=0.025, na.rm=TRUE),
					PHSC_PD_Q25=quantile(PATRISTIC_DISTANCE, p=0.25, na.rm=TRUE),
					PHSC_PD_Q75=quantile(PATRISTIC_DISTANCE, p=0.75, na.rm=TRUE),
					PHSC_PD_Q975=quantile(PATRISTIC_DISTANCE, p=0.975, na.rm=TRUE)	
			), by=c('PTY_RUN','MALE_RID','FEMALE_RID')]
	tmp[, PHSC_W:='all']
	tmp			<- rbind(tmp, tmp2)
	dfd			<- merge(dfd, tmp, by=c('MALE_RID','FEMALE_RID'), all=TRUE)
	#
	#	add sexual risk variables to couples
	#
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_info.rda")
	tmp		<- unique(subset(rp, select=c(FEMALE_RID, FEMALE_MARSTAT, FEMALE_ALC, FEMALE_SEXP1YR, FEMALE_SEXP1OUT, FEMALE_COMM_TYPE, MALE_RID, MALE_MARSTAT, MALE_ALC, MALE_SEXP1YR, MALE_SEXP1OUT, MALE_COMM_TYPE, COUP_SC)))
	dfd		<- merge(dfd, tmp, by=c('MALE_RID','FEMALE_RID'))	
	dfd[, EXTRA_PARTNERS:=NA_character_]
	set(dfd, dfd[, which(grepl('Monogamous$',MALE_MARSTAT) | grepl('Monogamous$',FEMALE_MARSTAT))], 'EXTRA_PARTNERS', 'none')
	set(dfd, dfd[, which(grepl('Monogamous \\+ casual',MALE_MARSTAT) | grepl('Monogamous \\+ casual',FEMALE_MARSTAT))], 'EXTRA_PARTNERS', 'at least one')
	set(dfd, dfd[, which(grepl('Monogamous \\+ casual',MALE_MARSTAT) & grepl('Monogamous \\+ casual',FEMALE_MARSTAT))], 'EXTRA_PARTNERS', 'both')
	set(dfd, NULL, 'EXTRA_PARTNERS', dfd[, factor(EXTRA_PARTNERS, levels=c('both','at least one','none'))])	
	dfd[, ALC:='no']
	set(dfd, dfd[, which(MALE_ALC=='Y'|FEMALE_ALC=='Y')],'ALC','yes')
	dfd[, SEXP1OUT:='unknown']
	set(dfd, dfd[, which(MALE_SEXP1OUT=='0'& FEMALE_SEXP1OUT=='0')],'SEXP1OUT','no')
	set(dfd, dfd[, which(MALE_SEXP1OUT!='0'| FEMALE_SEXP1OUT!='0')],'SEXP1OUT','yes')			
	set(dfd, NULL, 'SEXP1OUT', dfd[, factor(SEXP1OUT, levels=c('yes','no'))])
				
	#
	#	bimodality
	#	
	tmp		<- subset(dfd, PHSC_W=='all' & !is.na(GD) & !is.na(PHSC_PD_MEAN))
	tmp		<- melt(tmp, id.vars=c('MALE_RID','FEMALE_RID','PTY_RUN'), measure.vars=c('PHSC_PD_MEAN','GD'))
	set(tmp, NULL, 'variable', tmp[, as.character(factor(as.character(variable), levels=c('PHSC_PD_MEAN','GD'), labels=c('patristic distance between WH clades','genetic distance between consensus sequences')))])
	set(tmp, tmp[, which(value<1e-3)], 'value', 1e-3)	
	#tmp		<- subset(tmp, variable=='PHSC_PD_MEAN')
	ggplot(tmp, aes(x=log10(value), colour=variable)) +
			geom_vline(xintercept=log10(c(0.035,0.08)), colour='grey70') +
			geom_density(adjust=1.2, kernel='epanechnikov') +
			scale_colour_brewer(palette='Set1') +			
			scale_y_continuous(expand=c(0,0), limits=c(0,1)) + 
			scale_x_continuous(		limits=log10(c(0.0009, 1)), expand=c(0,0), 
					breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5)), 
					labels=paste0(100*c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5),'%')) +
			theme_bw() + theme(legend.position='bottom') +			
			labs(x='\nsubst/site', y='', colour='')
	ggsave(file=paste0(outfile.base,'_distances_consGenetic_bimodal.pdf'), w=5, h=5)
	#
	tmp		<- subset(dfd, PHSC_W=='all' & !is.na(PD) & !is.na(PHSC_PD_MEAN) & !is.na(EXTRA_PARTNERS) & COUP_SC!='seropos')
	tmp		<- subset(dfd, PHSC_W=='all' & !is.na(PD) & !is.na(PHSC_PD_MEAN) & FEMALE_COMM_TYPE!='trading' & !is.na(EXTRA_PARTNERS))
	ggplot(tmp, aes(x=log10(PHSC_PD_MEAN), fill=EXTRA_PARTNERS)) +
			geom_vline(xintercept=log10(c(0.035,0.08)), colour='grey70') +
			geom_histogram(binwidth=0.1) +
			scale_fill_brewer(palette='Dark2') +
			#scale_fill_manual(values=c('none'="#80CDC1", 'at least one'="#35978F", 'both'="#01665E")) +				 
			#scale_y_continuous(expand=c(0,0), limits=c(0,32)) + 
			scale_x_continuous(		limits=log10(c(0.0009, 1)), expand=c(0,0), 
					breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5)), 
					labels=paste0(100*c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5),'%')) +
			theme_bw() + theme(legend.position='bottom') +			
			labs(x='\nsubst/site', y='', fill='partners reporting extra-marital sexual contacts in last year') +
			facet_grid(FEMALE_COMM_TYPE~., scale='free_y')
	ggsave(file=paste0(outfile.base,'_distances_consPatristic_bimodal_extraPartners.pdf'), w=5, h=8)
	#
	dfds	<- subset(dfd, PHSC_W=='all' & !is.na(PD) & !is.na(PHSC_PD_MEAN))
	dfds	<- melt(dfds, id.vars=c('MALE_RID','FEMALE_RID','PTY_RUN'), measure.vars=c('PHSC_PD_MEAN','GD'))
	set(dfds, NULL, 'variable', dfds[, as.character(factor(as.character(variable), levels=c('PHSC_PD_MEAN','GD'), labels=c('patristic distance between WH clades','patristic distance between consensus sequences')))])
	set(dfds, dfds[, which(value<1e-3)], 'value', 1e-3)	
	#dfds	<- subset(dfds, variable=='PHSC_PD_MEAN')
	ggplot(dfds, aes(x=log10(value), colour=variable)) +
			geom_vline(xintercept=log10(c(0.035,0.08)), colour='grey70') +
			geom_density(adjust=1.2, kernel='epanechnikov') +
			scale_colour_brewer(palette='Set1') +			
			scale_y_continuous(expand=c(0,0), limits=c(0,1)) + 
			scale_x_continuous(		limits=log10(c(0.0009, 1)), expand=c(0,0), 
					breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5)), 
					labels=paste0(100*c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5),'%')) +
			theme_bw() + theme(legend.position='bottom') +			
			labs(x='\nsubst/site', y='', colour='')
	ggsave(file=paste0(outfile.base,'_distances_consPatristic_bimodal.pdf'), w=5, h=5)
	
	#
	#	fit 2 component mixture model
	#	
	require(mclust)
	dfds	<- subset(dfd, PHSC_W=='all' & !is.na(PD) & !is.na(PHSC_PD_MEAN) & PHSC_PD_MEAN>1e-5 & PHSC_PD_MEAN<1)
	m1		<- densityMclust( log(dfds$PHSC_PD_MEAN), G=2)
	pdf(file=paste0(outfile.base,'_distances_consPatristic_lognormalmixture_pdf.pdf'), w=5, h=5)
	plot(m1, what='density', data= log(dfds$PHSC_PD_MEAN), breaks=50)
	dev.off()
	pdf(file=paste0(outfile.base,'_distances_consPatristic_lognormalmixture_cdf.pdf'), w=5, h=5)
	densityMclust.diagnostic(m1, type='cdf')	
	dev.off()
	#	mean and variance of first component	
	tmp		<- log(seq(1e-4,1,0.0001))
	th1		<- unname(c( summary(m1, parameters=TRUE)$mean, summary(m1, parameters=TRUE)$variance, summary(m1, parameters=TRUE)$pro)[c(1,3,5)])
	dens1	<- data.table(X=tmp, Y=dnorm(tmp, mean=th1[1], sd=sqrt(th1[2])))
	th2		<- unname(c( summary(m1, parameters=TRUE)$mean, summary(m1, parameters=TRUE)$variance, summary(m1, parameters=TRUE)$pro)[c(2,4,6)])
	dens2	<- data.table(X=tmp, Y=dnorm(tmp, mean=th2[1], sd=sqrt(th2[2])))
	densm	<- data.table(X=tmp, Y=predict.densityMclust(m1, tmp))
	
	#	3.5% threshold corresponds to 9.14% quantile of the LogNormal
	pnorm(log(0.035), mean=th1[1], sd=sqrt(th1[2]), lower.tail=FALSE)
	#	3.5% threshold corresponds to 0.6% quantile of the LogNormal
	pnorm(log(0.035), mean=th2[1], sd=sqrt(th2[2]), lower.tail=TRUE)	
	#	8% threshold corresponds to 6.936269e-07 quantile of the LogNormal
	pnorm(log(0.08), mean=th2[1], sd=sqrt(th2[2]), lower.tail=TRUE)
	
	#	10% quantile is 3.3% threshold
	exp(qnorm(0.9, mean=th1[1], sd=sqrt(th1[2])))
	#tmp		<- seq(1e-4,0.1,0.0001)
	#tmp		<- data.table(X=tmp, Y=dlnorm(tmp, meanlog=th1[1], sdlog=sqrt(th1[2]), log=FALSE))
	
	ggplot(dfds) +
			annotate("rect", xmin=log(2e-4), xmax=log(0.035), ymin=-Inf, ymax=Inf, fill=brewer.pal(11, 'PuOr')[2], alpha=0.5) +
			annotate("rect", xmin=log(0.08), xmax=log(1), ymin=-Inf, ymax=Inf, fill=rev(brewer.pal(11, 'RdGy'))[4], alpha=0.5) +
			#geom_vline(xintercept=log(c(0.035,0.08)), colour='grey70') +
			#geom_vline(xintercept=qnorm(c(1e-3, 0.005, 0.01), mean=th2[1], sd=sqrt(th2[2])), colour='blue') +
			#geom_vline(xintercept=qnorm(c(0.8,0.9,0.95), mean=th1[1], sd=sqrt(th1[2])), colour='red') +
			geom_histogram(aes(x=log(PHSC_PD_MEAN), y=..density..), binwidth=0.2, colour='white', fill='skyblue3') +
			geom_line(data=densm, aes(x=X, y=Y), lwd=0.8, colour='black') +			
			#geom_line(data=dens1, aes(x=X, y=Y*0.57), lwd=1.25, colour=brewer.pal(11, 'PuOr')[2]) +
			#geom_line(data=dens2, aes(x=X, y=Y*0.432), lwd=1.25, colour=rev(brewer.pal(11, 'RdGy'))[4]) +			
			#scale_colour_brewer(palette='Set1') +			
			scale_y_continuous(expand=c(0,0), limits=c(0,0.52)) + 
			scale_x_continuous(limits=log(c(2e-4, 1)), expand=c(0,0), 
					breaks=log(c(0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5)), 
					labels=paste0(100*c(0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5),'%')) +
			theme_bw() + theme(legend.position='bottom') +			
			labs(x='\npatristic distance between within-host subtrees of spouses\n( subst/site )', y='density', colour='')
	ggsave(file=paste0(outfile.base,'_distances_consPatristic_lognormalfitted.pdf'), w=5, h=5)
	
	
	#
	#	correlation plot genetic distance
	ggplot(subset(dfd, PHSC_W=='all' & !is.na(GD)), aes(x=GD, y=PHSC_PD_MEAN, ymin=PHSC_PD_Q25, ymax=PHSC_PD_Q75)) +
			geom_rect(xmin=log10(0.035), xmax=log10(0.08), ymin=log10(0.0001), ymax=log10(1), fill='grey85') +
			geom_rect(ymin=log10(0.035), ymax=log10(0.08), xmin=log10(0.0001), xmax=log10(1), fill='grey85') +			
			geom_abline(slope=1, intercept=0, colour='black', linetype='dotted') +
			#geom_linerange(size=.5, alpha=0.5, colour='grey40') + 
			geom_point(size=1) + 
			scale_x_log10(labels=percent, expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			scale_y_log10(labels=percent, expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			coord_cartesian(xlim=c(0.002, 0.5), ylim=c(0.002, 0.5)) + 
			theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=15)) +
			labs(	x='\ngenetic distance between consensus sequences\n(subst/site)',
					y='patristic distance between WH clades\n(avg subst/site across windows)\n')
	ggsave(file=paste0(outfile.base,'_distances_consGenetic.pdf'), w=7, h=7)
	ggplot(subset(dfd, PHSC_W=='all' & !is.na(PD)), aes(x=PD, y=PHSC_PD_MEAN, ymin=PHSC_PD_Q25, ymax=PHSC_PD_Q75)) +
			geom_rect(xmin=log10(0.05), xmax=log10(0.12), ymin=log10(0.0001), ymax=log10(1), fill='grey85') +
			geom_rect(ymin=log10(0.035), ymax=log10(0.08), xmin=log10(0.0001), xmax=log10(1), fill='grey85') +
			geom_abline(slope=1, intercept=0, colour='black', linetype='dotted') +
			#geom_linerange(size=.5, alpha=0.5, colour='grey40') + 
			geom_point(size=1) + 
			scale_x_log10(labels=percent, expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			scale_y_log10(labels=percent, expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			coord_cartesian(xlim=c(0.002, 0.5), ylim=c(0.002, 0.5)) + 
			theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=15)) +
			labs(	x='\npatristic distance between consensus sequences\n(subst/site)',
					y='patristic distance between WH clades\n(avg subst/site across windows)\n')
	ggsave(file=paste0(outfile.base,'_distances_consPatristic.pdf'), w=7, h=7)		
	ggplot(subset(dfd, PHSC_W=='all' & !is.na(PD)), aes(x=PD, y=PHSC_PD_MEAN, ymin=PHSC_PD_Q25, ymax=PHSC_PD_Q75, colour=EXTRA_PARTNERS)) +
			geom_rect(xmin=log10(0.05), xmax=log10(0.12), ymin=log10(0.0001), ymax=log10(1), fill='grey85', colour='grey85') +
			geom_rect(ymin=log10(0.035), ymax=log10(0.08), xmin=log10(0.0001), xmax=log10(1), fill='grey85', colour='grey85') +
			geom_abline(slope=1, intercept=0, colour='black', linetype='dotted') +
			#geom_linerange(size=.5, alpha=0.5, colour='grey40') + 
			geom_point(size=1) +			
			scale_x_log10(labels=percent, expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			scale_y_log10(labels=percent, expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			coord_cartesian(xlim=c(0.002, 0.5), ylim=c(0.002, 0.5)) + 
			theme_bw() + theme(axis.title=element_text(size=15), axis.text=element_text(size=15)) +
			labs(	x='\npatristic distance between consensus sequences\n(subst/site)',
					y='patristic distance between WH clades\n(avg subst/site across windows)\n')
	ggsave(file=paste0(outfile.base,'_distances_consPatristic_extraPartners.pdf'), w=7, h=7)	
	
	
	
	subset(dfd, PHSC_W=='all' & !is.na(PD))[, cor(PD, PHSC_PD_MEAN)]
	#	0.4963385
	subset(dfd, PHSC_W=='all' & !is.na(GD))[, cor(GD, PHSC_PD_MEAN)]
	#	0.4727337
	subset(dfd, PHSC_W=='gag+pol' & !is.na(GD))[, cor(GD, PHSC_PD_MEAN)]
	#	0.4942915

	
	tmp	<- subset(dfd, PHSC_W=='all' & PD>0.05 & PHSC_PD_MEAN<0.035)
	subset(merge(rtp, tmp, by=c('MALE_RID','FEMALE_RID')), select=c(MALE_RID, FEMALE_RID, PTY_RUN.x))
	#	MALE_RID FEMALE_RID PTY_RUN.x
	#1:  A106044    C106054       146
	#2:  A108832    A108688       109
	#3:  B035048    J035045       214
	#4:  D030388    J104165       131
	#5:  H104368    G104325       206
	#6:  J106848    K107014        13
	#6:  F108382    F108764        43
	#5:  D066337    B066335       116
	tmp	<- subset(dfd, PHSC_W=='all' & PD<0.01 & PHSC_PD_MEAN>0.05)
	subset(merge(rca, tmp, by=c('MALE_RID','FEMALE_RID')), select=c(MALE_RID, FEMALE_RID, PTY_RUN.x))
	#	   MALE_RID FEMALE_RID PTY_RUN.x
	#	1:  B114802    H115007       165

	tmp	<- merge(rca, unique(subset(dfd, PHSC_W=='all', c(MALE_RID, FEMALE_RID, PD, GD))), by=c('MALE_RID','FEMALE_RID'), all.x=1)
	tmp[, table(SELECT,is.finite(PD))]
	#SELECT                                                               FALSE TRUE
	#couple ambiguous if pair or not pair                                   4   15
	#couple most likely a pair direction not resolved                       3   52
	#couple most likely a pair with resolved direction                     12   76
	#couple most likely not a pair                                         35  111
	#insufficient deep sequence data for at least one partner of couple   178    0
	tmp[, table(SELECT,is.finite(GD))]
	#SELECT                                                               FALSE TRUE
	#couple ambiguous if pair or not pair                                   0   19
	#couple most likely a pair direction not resolved                       0   55
	#couple most likely a pair with resolved direction                      0   88
	#couple most likely not a pair                                          0  146
	#insufficient deep sequence data for at least one partner of couple   178    0

	
	rex			<- subset(rplkl, GROUP=='TYPE_PAIR_DI')[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('MALE_RID','FEMALE_RID','PTY_RUN')]
	rex			<- merge(subset(rex, TYPE_MLE=='distant', c('MALE_RID','FEMALE_RID')), subset(rplkl, GROUP=='TYPE_PAIR_DI' & TYPE=='distant'), by=c('MALE_RID','FEMALE_RID'), all.x=1)
	rex[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]	
	rex			<- rex[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)], POSTERIOR_SCORE=POSTERIOR_SCORE[which.max(POSTERIOR_SCORE)]), by=c('MALE_RID','FEMALE_RID')]	
	rex			<- subset(rex, POSTERIOR_SCORE>confidence.cut)
	rex[, SELECT_DI:= 'couple most likely not a pair']	
	rtp			<- subset(rplkl, GROUP=='TYPE_PAIR_DI')[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('MALE_RID','FEMALE_RID','PTY_RUN')]
	rtp			<- merge(subset(rtp, TYPE_MLE=='close', c('MALE_RID','FEMALE_RID')), subset(rplkl, GROUP=='TYPE_PAIR_DI' & TYPE=='close'), by=c('MALE_RID','FEMALE_RID'), all.x=1)
	rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]	
	rtp			<- rtp[, list(PTY_RUN=PTY_RUN[which.max(POSTERIOR_SCORE)], POSTERIOR_SCORE=POSTERIOR_SCORE[which.max(POSTERIOR_SCORE)]), by=c('MALE_RID','FEMALE_RID')]	
	rtp			<- subset(rtp, POSTERIOR_SCORE>confidence.cut)
	rtp[, SELECT_DI:= 'couple most likely a pair']
	tmp			<- rbind(rex, rtp)
	tmp[, POSTERIOR_SCORE:=NULL]
	rca			<- merge(rca, tmp, by=c('MALE_RID','FEMALE_RID','PTY_RUN'), all.x=1)
	set(rca, rca[, which(SELECT=='insufficient deep sequence data for at least one partner of couple')], 'SELECT_DI', 'insufficient deep sequence data for at least one partner of couple')
	set(rca, rca[, which(is.na(SELECT_DI))], 'SELECT_DI', 'couple ambiguous if pair or not pair')
	
	
	tmp	<- merge(rca, unique(subset(dfd, PHSC_W=='all', c(MALE_RID, FEMALE_RID, PD, GD))), by=c('MALE_RID','FEMALE_RID'), all.x=1)
	tmp	<- subset(tmp, !is.na(PD))
	tmp[, PHSC_DI:=NA_character_]
	set(tmp, tmp[, which(grepl('couple ambiguous',SELECT_DI))], 'PHSC_DI', 'ambiguous')
	set(tmp, tmp[, which(grepl('couple most likely a pair',SELECT_DI))], 'PHSC_DI', 'close')
	set(tmp, tmp[, which(grepl('couple most likely not a pair',SELECT_DI))], 'PHSC_DI', 'distant')
	tmp[, CONS_DI:= cut(PD, breaks=c(-Inf, 0.05, 0.12, Inf), labels=c('close','ambiguous','distant'))]
	set(tmp, NULL, 'PHSC_DI', tmp[, factor(PHSC_DI, levels=c('close','ambiguous','distant'))])
	set(tmp, NULL, 'CONS_DI', tmp[, factor(CONS_DI, levels=c('close','ambiguous','distant'))])
	tmp[, table(CONS_DI, PHSC_DI, useNA='if')]
	#           PHSC_DI
	#CONS_DI     close ambiguous distant
  	#close       126         3       0
  	#ambiguous     9         2       0
  	#distant       3         2     109
}


RakaiFull.preprocess.closepairs.comparetocouples.170421	<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)	
	#	load denominator
	tmp		<- RakaiCirc.epi.get.info.170208()
	ra		<- tmp$ra		
	# load couples "rp"
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	rc		<- copy(rp)
	# load pty.run
	load( "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load rd, rh, rs, rp, rpw, rplkl, ptc
	load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/RCCS_170410_w250_trmp_allpairs_posteriors_cmptoprv.rda')
	#
	# select final run
	#
	run		<- "RCCS_170410_w250_d50_st20_trB_blInScriptNormed_mr20_mt1_cl3.5_d8"
	rpw		<- subset(rpw, RUN%in%run )
	rplkl	<- subset(rplkl, RUN%in%run )	
	#	add info on pair types to rplkl
	rp		<- copy(rpw)
	set(rp, NULL, c('DIR','FILE','RUN','W_FROM','W_TO','TYPE_RAW','TYPE','GROUP','PATRISTIC_DISTANCE','ADJACENT','CONTIGUOUS','PATHS_12','PATHS_21','MALE_SANGER_ID_L','MALE_SANGER_ID_R','FEMALE_SANGER_ID_L','FEMALE_SANGER_ID_R','CHUNK','CHUNK_L','CHUNK_N','ID_R_MIN','ID_R_MAX'), NULL)
	rp		<- unique(rp)
	#	make COUPID
	rp[, COUPID:= paste0(MALE_RID,':',FEMALE_RID)]	
	#	add PAIR_TYPE
	tmp		<- unique(subset(rc, select=c(COUPID, MALE_HH_NUM, FEMALE_HH_NUM, COUP_SC, PAIR_TYPE)))	
	setnames(tmp, 'COUP_SC', 'COUP_TYPE')
	set(tmp, NULL, c('MALE_HH_NUM','FEMALE_HH_NUM'), NULL)
	rp		<- merge(rp, tmp, by=c('COUPID'),all.x=1)	
	set(rp, rp[, which(!MALE_RID%in%rc[, MALE_RID] & !FEMALE_RID%in%rc[, FEMALE_RID])], 'PAIR_TYPE', 'm and f not in couple')
	set(rp, rp[, which(is.na(PAIR_TYPE))], 'PAIR_TYPE', 'f or m not in couple')	
	tmp		<- subset(rp, select=c(FEMALE_SANGER_ID, MALE_SANGER_ID, MALE_RID, FEMALE_RID, COUPID, PTY_RUN, COUP_TYPE, PAIR_TYPE))
	set(rplkl, NULL, c('MALE_RID','FEMALE_RID','COUPID','COUP_TYPE','PAIR_TYPE'), NULL)
	rplkl	<- merge(tmp, rplkl, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN'))
	set(rplkl, NULL, 'FEMALE_SANGER_ID', rplkl[, as.character(FEMALE_SANGER_ID)])
	set(rplkl, NULL, 'MALE_SANGER_ID', rplkl[, as.character(MALE_SANGER_ID)])
	rplkl	<- unique(rplkl)
	set(rpw, NULL, 'FEMALE_SANGER_ID', rpw[, as.character(FEMALE_SANGER_ID)])
	set(rpw, NULL, 'MALE_SANGER_ID', rpw[, as.character(MALE_SANGER_ID)])		
	#	select likely transmitters (unsampled intermediate not necessarily excluded) 
	#	find pairs for whom 'likely pair' is most likely state
	#	(does not depend on prior or confidence cut)
	mle.group	<- 'TYPE_PAIR_DI'
	mle.state	<- 'close'
	conf.group	<- 'TYPE_PAIR_DI'
	conf.state	<- 'close'
	rtpc		<- subset(rplkl, GROUP==mle.group)[, list(TYPE_MLE=TYPE[which.max(KEFF)], KEFF=max(KEFF)), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','COUPID')]
	rtpc		<- subset(rtpc, TYPE_MLE==mle.state)
	#	select one sequence pairing per couple: that with highest evidence
	rtpc		<- rtpc[, {
				z<- which.max(KEFF)
				list(MALE_SANGER_ID=MALE_SANGER_ID[z], FEMALE_SANGER_ID=FEMALE_SANGER_ID[z], PTY_RUN=PTY_RUN[z])
			}, by='COUPID']
	set(rtpc, NULL, 'COUPID', NULL)
	#	calculate confidence score and select	
	rtpc		<- merge(rtpc, subset(rplkl, GROUP==conf.group & TYPE==conf.state), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all.x=1)	
	rtpc[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
	rtpc2		<- subset(rtpc, select=c(PTY_RUN, MALE_RID, FEMALE_RID, KEFF, NEFF, POSTERIOR_ALPHA, POSTERIOR_BETA, POSTERIOR_SCORE))
	tmp			<- setdiff(colnames(rtpc2),c('MALE_RID','FEMALE_RID'))
	setnames(rtpc2, tmp, paste0(tmp,'_COUPLES'))
	#
	#	compare to first batch FULL RUN
	#	this makes sense because the previous couples run is exactly those individuals that are in the first batch full run
	#	
	load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170421.rda') 
	set(rtp, NULL, c('TYPE_MLE','GROUP','TYPE','N_TYPE','PAR_PRIOR','K','N'), NULL)
	#	same couple can be in one batch and hence be listed more than once
	rtp		<- rtp[, list(PTY_RUN=PTY_RUN[1], KEFF=mean(KEFF), NEFF=mean(NEFF), POSTERIOR_ALPHA=mean(POSTERIOR_ALPHA), POSTERIOR_BETA=mean(POSTERIOR_BETA), POSTERIOR_SCORE=mean(POSTERIOR_SCORE)), by=c('ID1','ID2')]
	#	from first batch we have 1683 close pairs -- that s exciting!!

	tmp		<- copy(rtp)
	setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
	rtp		<- rbind(rtp, tmp)
	setnames(rtp, c('ID1','ID2'), c('MALE_RID','FEMALE_RID'))
	tmp		<- setdiff(colnames(rtp),c('MALE_RID','FEMALE_RID'))
	setnames(rtp, tmp, paste0(tmp,'_FULL'))
	rtp		<- merge(rtpc2, rtp, by=c('MALE_RID','FEMALE_RID'), all.x=1)
	
	rtp[, mean(!is.na(POSTERIOR_SCORE_FULL))]
	#[1] 0.9471545
	#	--> 95% of pairs in couples run that are most likely close are also most likely close in full run 
	ggplot( subset(rtp, !is.na(POSTERIOR_SCORE_FULL)), aes(x=POSTERIOR_SCORE_COUPLES, y=POSTERIOR_SCORE_FULL)) +
			geom_point() +
			theme_bw()
	#	this is not too bad
	#	we could probably run on all MLE without confidence.cut anyhow?
	#	PLUS we could probably be a bit generous on the close cut-off?

}

RakaiAll.analyze.pairs.170418.direction<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	run		<- 'RCCS_170410_w250_trB_blNormedOnFly_dirlklprs_'
	dir		<- '/Users/Oliver/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410'	
	#	load denominator
	tmp		<- RakaiCirc.epi.get.info.170208()
	ra		<- tmp$ra		
	# load couples "rp"
	load("~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	rc		<- copy(rp)
	# load pty.run
	load( "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load rd, rh, rs, rp, rpw, rplkl, ptc
	load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/RCCS_170410_w250_trmp_allpairs_posteriors_cmptoprv.rda')
	#
	# select final run
	#
	tmp		<- "RCCS_170410_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl3.5_d8"
	tmp		<- "RCCS_170410_w250_d50_st20_trB_blInScriptNormed_mr20_mt1_cl3.5_d8"
	rpw		<- subset(rpw, RUN%in%tmp )
	rplkl	<- subset(rplkl, RUN%in%tmp )	
	#	add info on pair types to rplkl
	rp		<- copy(rpw)
	set(rp, NULL, c('DIR','FILE','RUN','W_FROM','W_TO','TYPE_RAW','TYPE','GROUP','PATRISTIC_DISTANCE','ADJACENT','CONTIGUOUS','PATHS_12','PATHS_21','MALE_SANGER_ID_L','MALE_SANGER_ID_R','FEMALE_SANGER_ID_L','FEMALE_SANGER_ID_R','CHUNK','CHUNK_L','CHUNK_N','ID_R_MIN','ID_R_MAX'), NULL)
	rp		<- unique(rp)
	#	make COUPID
	rp[, COUPID:= paste0(MALE_RID,':',FEMALE_RID)]	
	#	add PAIR_TYPE
	tmp		<- unique(subset(rc, select=c(COUPID, MALE_HH_NUM, FEMALE_HH_NUM, COUP_SC, PAIR_TYPE)))	
	setnames(tmp, 'COUP_SC', 'COUP_TYPE')
	set(tmp, NULL, c('MALE_HH_NUM','FEMALE_HH_NUM'), NULL)
	rp		<- merge(rp, tmp, by=c('COUPID'),all.x=1)	
	set(rp, rp[, which(!MALE_RID%in%rc[, MALE_RID] & !FEMALE_RID%in%rc[, FEMALE_RID])], 'PAIR_TYPE', 'm and f not in couple')
	set(rp, rp[, which(is.na(PAIR_TYPE))], 'PAIR_TYPE', 'f or m not in couple')	
	tmp		<- subset(rp, select=c(FEMALE_SANGER_ID, MALE_SANGER_ID, MALE_RID, FEMALE_RID, COUPID, PTY_RUN, COUP_TYPE, PAIR_TYPE))
	set(rplkl, NULL, c('MALE_RID','FEMALE_RID','COUPID','COUP_TYPE','PAIR_TYPE'), NULL)
	rplkl	<- merge(tmp, rplkl, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN'))
	set(rplkl, NULL, 'FEMALE_SANGER_ID', rplkl[, as.character(FEMALE_SANGER_ID)])
	set(rplkl, NULL, 'MALE_SANGER_ID', rplkl[, as.character(MALE_SANGER_ID)])
	rplkl	<- unique(rplkl)
	set(rpw, NULL, 'FEMALE_SANGER_ID', rpw[, as.character(FEMALE_SANGER_ID)])
	set(rpw, NULL, 'MALE_SANGER_ID', rpw[, as.character(MALE_SANGER_ID)])		
	#
	#	basic info on selection
	#
	if(0)
	{
		tmp		<- unique(subset(rc, !is.na(MALE_TAXA) & !is.na(FEMALE_TAXA) & PAIR_TYPE=='stable cohabiting'), by='COUPID')
		z		<- unique(subset(ri, select=c(COMM_NUM, COMM_TYPE)))
		setnames(z, c('COMM_NUM','COMM_TYPE'),c('MALE_COMM_NUM','MALE_COMM_TYPE'))
		tmp		<- merge(tmp, z, by= 'MALE_COMM_NUM')
		setnames(z, c('MALE_COMM_NUM','MALE_COMM_TYPE'), c('FEMALE_COMM_NUM','FEMALE_COMM_TYPE'))
		tmp		<- merge(tmp, z, by= 'FEMALE_COMM_NUM')
		
		tmp		<- unique(rp, by='COUPID')
		tmp[, c(length(unique(MALE_RID)),length(unique(FEMALE_RID)))]
		nrow(unique(subset(tmp, PAIR_TYPE=='stable cohabiting'), by=c('FEMALE_RID','MALE_RID')))
		subset(tmp, PAIR_TYPE=='stable cohabiting')[, c(length(unique(MALE_RID)),length(unique(FEMALE_RID)))]
		nrow(subset(tmp, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_NUM!=MALE_COMM_NUM))
		subset(tmp, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_NUM==MALE_COMM_NUM)[, table(FEMALE_COMM_TYPE)]
		subset(tmp, !is.na(COUP_TYPE))[, c(length(unique(MALE_RID)),length(unique(FEMALE_RID)),length(unique(COUPID)))]		
		#	2 males with new partners:  G110085 J189465
		#	4 females with new partners:  C106054 H104287 B105985 E111070
	}
	
	#	select likely transmitters (unsampled intermediate not necessarily excluded) 
	#	find pairs for whom 'likely pair' is most likely state
	#	(does not depend on prior or confidence cut)
	rex		<- subset(rplkl, GROUP=='TYPE_PAIR_TODI')[, list(TYPE_MLE=TYPE[which.max(KEFF)], KEFF=max(KEFF)), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','COUPID')]
	rex		<- subset(rex, TYPE_MLE=='distant')
	#	select one sequence pairing per couple: that with highest evidence
	rex		<- rex[, {
				z<- which.max(KEFF)
				list(MALE_SANGER_ID=MALE_SANGER_ID[z], FEMALE_SANGER_ID=FEMALE_SANGER_ID[z], PTY_RUN=PTY_RUN[z])
			}, by='COUPID']
	set(rex, NULL, 'COUPID', NULL)
	#	calculate confidence score and select
	confidence.cut	<- 0.5
	rex		<- merge(rex, subset(rplkl, GROUP=='TYPE_PAIRSCORE_TODI' & TYPE=='distant'), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all.x=1)	
	rex[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
	rex		<- subset(rex, POSTERIOR_SCORE>confidence.cut)
	
	
	#	select likely transmitters (unsampled intermediate not necessarily excluded) 
	#	find pairs for whom 'likely pair' is most likely state
	#	(does not depend on prior or confidence cut)
	rtp		<- subset(rplkl, GROUP=='TYPE_PAIR_TODI')[, list(TYPE_MLE=TYPE[which.max(KEFF)], KEFF=max(KEFF)), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','COUPID')]
	rtp		<- subset(rtp, TYPE_MLE=='likely pair')
	#	select one sequence pairing per couple: that with highest evidence
	rtp		<- rtp[, {
				z<- which.max(KEFF)
				list(MALE_SANGER_ID=MALE_SANGER_ID[z], FEMALE_SANGER_ID=FEMALE_SANGER_ID[z], PTY_RUN=PTY_RUN[z])
			}, by='COUPID']
	set(rtp, NULL, 'COUPID', NULL)
	#	calculate confidence score and select
	confidence.cut	<- 0.5
	rtp		<- merge(rtp, subset(rplkl, GROUP=='TYPE_PAIRSCORE_TODI' & TYPE=='likely pair'), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all.x=1)	
	rtp[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
	rtp		<- subset(rtp, POSTERIOR_SCORE>confidence.cut)
	
	#	resolve direction
	#	find likely pairs for whom 'mf' or 'fm' is most likely state
	#	(does not depend on prior or confidence cut)
	rtpd	<- subset(rtp, select=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	rtpd	<- merge(rtpd, subset(rplkl, GROUP=='TYPE_DIR_TODI3'), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rtpd	<- rtpd[, list(TYPE_MLE=TYPE[which.max(KEFF)]), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	rtpd	<- subset(rtpd, TYPE_MLE!='ambiguous')	
	#	calculate confidence score and select
	confidence.cut	<- 0.5
	setnames(rtpd, 'TYPE_MLE','TYPE')
	rtpd	<- merge(rtpd, subset(rplkl, GROUP=='TYPE_DIRSCORE_TODI3'), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','TYPE'))
	rtpd[, POSTERIOR_SCORE:=pbeta(1/N_TYPE+(1-1/N_TYPE)/(N_TYPE+1), POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
	rtpd	<- subset(rtpd, POSTERIOR_SCORE>confidence.cut)
	rmf		<- subset(rtpd, TYPE=='mf')
	rfm		<- subset(rtpd, TYPE=='fm')
	
	#subset(rtp, PTY_RUN==28 & MALE_SANGER_ID=='15714_1_84' & FEMALE_SANGER_ID=='15862_1_86')
	#subset(rtpd, PTY_RUN==67 & MALE_SANGER_ID=='15965_1_24' & FEMALE_SANGER_ID=='15977_1_52')
	#subset(rplkl, PTY_RUN==28 & MALE_SANGER_ID=='15714_1_84' & FEMALE_SANGER_ID=='15862_1_86' & GROUP=='TYPE_DIRSCORE_TODI3')
	#subset(rplkl, PTY_RUN==28 & MALE_SANGER_ID=='15714_1_84' & FEMALE_SANGER_ID=='15862_1_86' & GROUP=='TYPE_DIR_TODI7x3')
	
	#	info		
	cat('\ncouples with phyloscanner assessment, n=',				nrow(unique(rplkl, by='COUPID')))	
	cat('\ncouples not implicated in transmission, n=',				nrow(unique(rex, by='COUPID')))
	unique(rex,by='COUPID')[, table(PAIR_TYPE)]
	cat('\ncouples that are likely pairs, n=',						nrow(unique(rtp, by='COUPID')))
	cat('\ncouples that are likely pairs with evidence M->F, n=',	nrow(unique(rmf, by='COUPID')))
	cat('\ncouples that are likely pairswith evidence F->M, n=',	nrow(unique(rfm, by='COUPID')))
	#	pairings assessed, n= 1741
	#	couples not implicated in transmission, n= 1402
	#		    f or m not in couple m and f not in couple not always cohabiting     stable cohabiting 
    #             1311                    16                     8                    67  
	#	couples that are likely pairs, n= 209
	#	likely direction resolved, n= 127
	#	   not always cohabiting 			   not registered as couple        stable cohabiting 
	#                  2                       41                              85
	#	couples that are likely pairs with evidence M->F, n= 84
	#	couples that are likely pairswith evidence F->M, n= 43
	
	#	define two helper data.table
	rmf		<- merge(unique(subset(rmf, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID))), rp, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rmf[, PHSC_DIR:='m->f']
	rfm		<- merge(unique(subset(rfm, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID))), rp, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rfm[, PHSC_DIR:='f->m']
	rtr		<- rbind(rmf, rfm)	
	rtr[, AGEDIFF:= FEMALE_BIRTHDATE-MALE_BIRTHDATE]
	rtr[, AVGAGE:= (MALE_BIRTHDATE+FEMALE_BIRTHDATE)/2]	
	rtr2	<- copy(rmf)
	setnames(rtr2,colnames(rtr2),gsub('FEMALE','REC',colnames(rtr2)))
	setnames(rtr2,colnames(rtr2),gsub('MALE','TR',colnames(rtr2)))
	tmp		<- copy(rfm)
	setnames(tmp,colnames(tmp),gsub('FEMALE','TR',colnames(tmp)))
	setnames(tmp,colnames(tmp),gsub('MALE','REC',colnames(tmp)))
	rtr2	<- rbind(rtr2,tmp)
	
	if(0)	#for inscript
	{
		rps			<- subset(rtr, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID))
		outfile.base<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/lklpr_TODI_'
		write.csv(rps, file=paste0(outfile.base,'_summary_versionmaxscore_inscript.csv'))
		#
		#	get difference from manual check (versionstrict2) to more relaxed version
		#		
		tmp2		<- subset(as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/lklpr_TODI__summary_versionmaxscore.csv')), select=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
		tmp2[, VERSION:='on the fly']
		tmp			<- copy(rps)	
		tmp[, VERSION_NEW:='in script']
		tmp			<- merge(tmp, tmp2, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all=1)
		tmp[, table(is.na(VERSION_NEW), is.na(VERSION))]
		#		FALSE TRUE
		#FALSE   109    4
		#TRUE     18    0
		rps			<- subset(tmp, is.na(VERSION) | is.na(VERSION_NEW))	
		write.csv(rps, file=paste0(outfile.base,'_summary_diffaftermaxscore_inscript.csv'))
		set(rps, NULL, c('VERSION','VERSION_NEW'), NULL)
		group		<- 'TYPE_DIR_TODI7x3'
		#group		<- 'TYPE_PAIR_TODI'
		run			<- "RCCS_170410_w250_d50_st20_trB_blInScriptNormed_mr20_mt1_cl3.5_d8"	
		plot.select	<- unique(subset(merge(rplkl, rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP==group), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
		rpw2		<- subset(rpw, RUN==run & GROUP==group)
		rplkl2		<- subset(rplkl, RUN==run & GROUP==group)	
		plot.file	<- paste0(outfile.base,'_windows_summary_',group,'_diffaftermaxscore_inscript.pdf')	
		phsc.plot.windowsummaries.for.pairs(plot.select, rpw2, rplkl2, plot.file, cols=NULL, group=group)
	}	
	if(0)	#for onthefly
	{
		rps			<- subset(rtr, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID))
		outfile.base<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/lklpr_TODI_'
		write.csv(rps, file=paste0(outfile.base,'_summary_versionmaxscore.csv'))
		group		<- 'TYPE_DIR_TODI7x3'
		#group		<- 'TYPE_PAIR_TODI'
		run			<- "RCCS_170410_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl3.5_d8"	
		plot.select	<- unique(subset(merge(rplkl, rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP==group), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
		rpw2		<- subset(rpw, RUN==run & GROUP==group)
		rplkl2		<- subset(rplkl, RUN==run & GROUP==group)	
		plot.file	<- paste0(outfile.base,'_windows_summary_',group,'_versionmaxscore.pdf')	
		phsc.plot.windowsummaries.for.pairs(plot.select, rpw2, rplkl2, plot.file, cols=NULL, group=group)
		#
		#	get difference from manual check (versionstrict2) to more relaxed version
		#		
		tmp2		<- subset(as.data.table(read.csv('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/lklpr_TODI__summary_versionstrict2.csv')), select=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
		tmp2[, VERSION:='checked manually']
		tmp			<- copy(rps)	
		tmp[, VERSION_NEW:='first max then score']
		tmp			<- merge(tmp, tmp2, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all=1)
		tmp[, table(is.na(VERSION_NEW), is.na(VERSION))]
		#		FALSE TRUE
		#   FALSE    59   68
		rps			<- subset(tmp, is.na(VERSION) | is.na(VERSION_NEW))	
		write.csv(rps, file=paste0(outfile.base,'_summary_diffaftermaxscore.csv'))
		set(rps, NULL, c('VERSION','VERSION_NEW'), NULL)
		group		<- 'TYPE_DIR_TODI7x3'
		group		<- 'TYPE_PAIR_TODI'
		run			<- "RCCS_170410_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl3.5_d8"	
		plot.select	<- unique(subset(merge(rplkl, rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP==group), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
		rpw2		<- subset(rpw, RUN==run & GROUP==group)
		rplkl2		<- subset(rplkl, RUN==run & GROUP==group)	
		plot.file	<- paste0(outfile.base,'_windows_summary_',group,'_diffaftermaxscore.pdf')	
		phsc.plot.windowsummaries.for.pairs(plot.select, rpw2, rplkl2, plot.file, cols=NULL, group=group)		
	}
	
	
	
	
	#	plot phylogenies for a few examples
	tmp			<- subset(rps, PTY_RUN%in%c(113) & FEMALE_SANGER_ID=='15958_1_47')
	set(tmp, NULL, 'FEMALE_SANGER_ID', tmp[, as.character(FEMALE_SANGER_ID)])
	set(tmp, NULL, 'MALE_SANGER_ID', tmp[, as.character(MALE_SANGER_ID)])
	run			<- 'RCCS_170410_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl3.5_d8'
	rpw2		<- unique(subset(rpw, RUN==run, select=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM','W_TO','PATHS_12','PATHS_21','PATRISTIC_DISTANCE','CONTIGUOUS','TYPE_RAW')))
	set(rpw2, NULL, 'FEMALE_SANGER_ID', rpw2[, as.character(FEMALE_SANGER_ID)])
	set(rpw2, NULL, 'MALE_SANGER_ID', rpw2[, as.character(MALE_SANGER_ID)])
	set(rpw2, NULL, 'CONTIGUOUS', rpw2[, as.integer(CONTIGUOUS)])
	set(rpw2, NULL, 'PATRISTIC_DISTANCE', rpw2[, round(PATRISTIC_DISTANCE, d=4)])	
	#load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/RCCS_170410_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl3.5_d8_phscout.rda')
	invisible(sapply(seq_len(nrow(tmp)), function(ii)
					{	
						#ii<- 1
						ids			<- c(tmp[ii, MALE_SANGER_ID],tmp[ii, FEMALE_SANGER_ID])
						pty.run		<- tmp[ii, PTY_RUN]
						dfs			<- subset(dtrees, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
						dfs[, MALE_SANGER_ID:=ids[1]]
						dfs[, FEMALE_SANGER_ID:=ids[2]]
						dfs			<- merge(dfs, rpw2, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM','W_TO'), all.x=TRUE)
						dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', pty.run, '\nwindow ', W_FROM,'-', W_TO,'\n',PATHS_12,' ',PATHS_21, ' ',CONTIGUOUS,' ',TYPE_RAW, '\n', PATRISTIC_DISTANCE, sep='')]]								
						plot.file	<- paste0(outfile.base, pty.run,'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
						invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10))					
					}))	
	
	subset(rplkl, PTY_RUN==44 & MALE_SANGER_ID=='15743_1_40' & FEMALE_SANGER_ID=='15115_1_22' & GROUP=='TYPE_PAIRSCORE_TODI' & TYPE=='likely pair')[, pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
	subset(rplkl, PTY_RUN==44 & MALE_SANGER_ID=='15743_1_40' & FEMALE_SANGER_ID=='15115_1_22' & GROUP=='TYPE_DIRSCORE_TODI3' & TYPE=='fm')[, pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)]
	
	#
	#rtr2[, table(PAIR_TYPE)]
	#
	#	who infects whom matrix
	#
	tmp		<- rtr2[,list(N=length(unique(COUPID))), by=c('TR_COMM_TYPE','REC_COMM_TYPE')]
	ggplot(tmp, aes(x=factor(REC_COMM_TYPE),y=factor(TR_COMM_TYPE))) + 
			geom_point(aes(size=N), colour='grey80') +
			geom_text(aes(label=N), nudge_x=0, nudge_y=0, size=3, colour='black') +			
			theme_bw() + 
			scale_size(range = c(5, 50)) +
			labs(x='\nlocation likely recipient',y='location likely transmitter\n') +
			guides(size='none')
	ggsave(file=file.path(dir, paste(run,'-phsc-directionpairs_direction-numbers-commtype.pdf',sep='')), w=4, h=4)
	#
	#	did any transmitter start ART before the recipient was diagnosed?
	subset(rtr2, TR_ARVSTARTDATE<REC_FIRSTPOSDATE)	
	#	F026858:J104288 --> stable couple, rec male, first diagnosed with v high CD4 (2400), about 2 years after female started ART 
	#	C066263:K077878 --> no couple, rec female, first diagnosed with v high CD4, about 5m after male started ART
	
	#
	#	how many transmitters were positive for 6m before the recipient was found positive
	subset(rtr2, (TR_FIRSTPOSDATE+.5)<REC_FIRSTPOSDATE)	
	#	26
	
	
	subset(rtr, MALE_COMM_TYPE==FEMALE_COMM_TYPE & FEMALE_COMM_TYPE!='trading')[, {
				z	<- binconf( length(which(PHSC_DIR=='m->f')), length(PHSC_DIR) )				
				list(K=length(which(PHSC_DIR=='m->f')), N=length(PHSC_DIR), P=z[1], QL=z[2], QU=z[3])
			}, by=c('MALE_COMM_TYPE')]	
	#	   MALE_COMM_TYPE  	K  N  P         QL        QU
	#1:       agrarian 		27 38 0.7105263 0.5524286 0.8299672
	#2:     fisherfolk 		55 85 0.6470588 0.5411250 0.7402751
	
	#
	#	is there a difference in male->female transmission by couple type?
	#	results: no	
	#
	tmp		<- copy(rtr)
	set(tmp, tmp[, which(PAIR_TYPE!='stable cohabiting')], 'PAIR_TYPE', 'no stable pair')
	tmp[, {
				z	<- binconf( length(which(PHSC_DIR=='m->f')), length(PHSC_DIR) )				
				list(K=length(which(PHSC_DIR=='m->f')), N=length(PHSC_DIR), P=z[1], QL=z[2], QU=z[3])
			}, by='PAIR_TYPE']	
	#			PAIR_TYPE  	 K  N         P        QL        QU
	#1:    stable cohabiting 59 86 0.6860465 0.5817960 0.7743870
	#2:    no stable pair 	 25 42 0.5952381 0.4449431 0.7295714	
	chisq_test(factor(PHSC_DIR) ~ factor(PAIR_TYPE), data=tmp, distribution="exact")
	#	chi-squared = 1.0315, p-value = 0.3278
	
	
	#	are transmitters younger in fisherfolk sites?
	#	results: yes
	tmp		<- subset(rtr2, TR_COMM_TYPE!='trading')
	ggplot(tmp, aes(x=TR_COMM_TYPE, y=TR_BIRTHDATE)) + geom_boxplot()
	independence_test(TR_BIRTHDATE~factor(TR_COMM_TYPE), data=tmp, distribution = "exact")
	#	Z = -2.3289, p-value = 0.01934
	ggplot(tmp, aes(x=REC_COMM_TYPE, y=REC_BIRTHDATE)) + geom_boxplot()
	independence_test(REC_BIRTHDATE~factor(REC_COMM_TYPE), data=tmp, distribution = "exact")
	#	Z = -2.2714, p-value = 0.02236
	#	summary(rq(TR_BIRTHDATE~TR_COMM_TYPE, tau=.5, data=tmp, method='fn'), se='nid')
	
	#
	#	is there a difference in age gap between male->female transmission / female->male transmission ?
	#	results: not significant but outside couples, men tend to be infected by much younger women
	#	
	tmp		<- subset(rtr, !grepl('cohabiting',PAIR_TYPE) & !is.na(AGEDIFF) & FEMALE_COMM_TYPE==MALE_COMM_TYPE & MALE_COMM_TYPE=='fisherfolk')
	independence_test(AGEDIFF~factor(PHSC_DIR), data=tmp, distribution = "exact")
	#	Z = 1.5902, p-value = 0.1134
	tmp		<- subset(rtr, !grepl('cohabiting',PAIR_TYPE) & !is.na(AGEDIFF) & FEMALE_COMM_TYPE==MALE_COMM_TYPE & MALE_COMM_TYPE=='agrarian')
	independence_test(AGEDIFF~factor(PHSC_DIR), data=tmp, distribution = "exact")
	#	Z = 1.7439, p-value = 0.1429
	
	
	ggplot(rtr, aes(x=PHSC_DIR, y=AGEDIFF)) + geom_boxplot()	
	tmp		<- subset(rtr, !is.na(MALE_BIRTHDATE) & !is.na(FEMALE_BIRTHDATE), select=c(PHSC_DIR,AGEDIFF))
	set(tmp, NULL, 'PHSC_DIR', tmp[, as.integer(as.character(factor(PHSC_DIR, levels=c('f->m','m->f'), labels=c('0','1'))))])
	summary(gamlss(data=tmp, PHSC_DIR~AGEDIFF, family=LO))
	summary(gamlss(data=tmp, AGEDIFF~PHSC_DIR))
	#				Estimate Std. Error t value Pr(>|t|)    
	#(Intercept)  0.626998   0.071579   8.759 1.98e-13 ***
	#AGEDIFF     -0.002135   0.008422  -0.254      0.8    
	tmp		<- subset(rtr, MALE_COMM_TYPE!='trading' & MALE_COMM_TYPE==FEMALE_COMM_TYPE)
	ggplot(tmp, aes(x=PHSC_DIR, y=AGEDIFF)) + 
			geom_boxplot() + 
			facet_grid(~MALE_COMM_TYPE)	+
			theme_bw() + labs(x='\nestimated direction of transmission', y='age difference male-female\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-directionpairs_direction-agegap-commtype.pdf',sep='')), w=4, h=6)
	
	#	AAA
	subset(rtr2, TR_COMM_TYPE!='trading' & PAIR_TYPE!='m and f not in couple' & REC_RID%in%c(rc$MALE_RID,rc$FEMALE_RID))[, {
				m<- length(which(PAIR_TYPE=='stable cohabiting'))
				n<- length(PAIR_TYPE)
				z<- unname(as.numeric(binconf(m, n)))
				list(P=z[1], PL=z[2], PU=z[3], M=m, N=n, TYPE='stable cohabiting')				
			}, by=c('TR_COMM_TYPE','REC_SEX')]
	
	#
	#	Does the primary occupation differ between transmitters / recipients? 
	#	
	tmp			<- subset(rtr, FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_occupation2.pdf',sep=''))
	tmp			<- subset(rtr, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_occupation2_stablecouples.pdf',sep=''))
	tmp			<- subset(rtr, !grepl('cohabiting',PAIR_TYPE) & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_occupation2_nocouples.pdf',sep=''))
	#	
	tmp2		<- unique(subset(ra, VISIT!=17 & SEX=='F' & !is.na(FIRSTPOSDATE) & COMM_TYPE!='trading', c(RID, OCCUP_OLLI, COMM_TYPE)))
	setnames(tmp2, colnames(tmp2), paste0('FEMALE_',colnames(tmp2)))
	tmp2[, PHSC_DIR:='denominator']
	tmp			<- rbind(tmp, tmp2, fill=TRUE)
	tmp2		<- unique(subset(ra, VISIT!=17 & SEX=='M' & !is.na(FIRSTPOSDATE) & COMM_TYPE!='trading', c(RID, OCCUP_OLLI, COMM_TYPE)))
	setnames(tmp2, colnames(tmp2), paste0('MALE_',colnames(tmp2)))
	tmp2[, PHSC_DIR:='denominator']
	tmp			<- rbind(tmp, tmp2, fill=TRUE)
	setnames(tmp, c('FEMALE_OCCUP_OLLI','MALE_OCCUP_OLLI'), c('FEMALE_FACTOR','MALE_FACTOR'))
	cols		<- tmp[, length(na.omit(unique(c(FEMALE_FACTOR, MALE_FACTOR))))]
	cols		<- colorRampPalette(brewer.pal(min(11,cols), "Set3"))( cols )		
	names(cols)	<- tmp[, na.omit(sort(unique(c(FEMALE_FACTOR, MALE_FACTOR))))]
	Rakai.plot.directed.pairs.discrete1(tmp, cols, outfile, 'occupation at diagnosis', w=10, h=7)
	
	# number female Bar/waitress that are transmitters
	ntf	<- nrow(subset(tmp, PHSC_DIR=='f->m' & FEMALE_FACTOR=='Bar/waitress' & FEMALE_COMM_TYPE=='fisherfolk'))
	# number female Bar/waitress that are recipients
	nrf	<- nrow(subset(tmp, PHSC_DIR=='m->f' & FEMALE_FACTOR=='Bar/waitress' & FEMALE_COMM_TYPE=='fisherfolk'))
	# number female Bar/waitress HIV-infected
	ndf	<- nrow(subset(tmp, PHSC_DIR=='denominator' & FEMALE_FACTOR=='Bar/waitress' & FEMALE_COMM_TYPE=='fisherfolk'))
	# number female that are transmitters
	nt	<- nrow(subset(tmp, PHSC_DIR=='f->m' & FEMALE_COMM_TYPE=='fisherfolk'))
	# number female that are recipients
	nr	<- nrow(subset(tmp, PHSC_DIR=='m->f' & FEMALE_COMM_TYPE=='fisherfolk'))
	# number female HIV-infected
	nd	<- nrow(subset(tmp, PHSC_DIR=='denominator' & FEMALE_COMM_TYPE=='fisherfolk'))
	# odds ratio transmitter / recipient
	# 'a' is exposed cases (exposed=Bar/waitress, case=transmitter)
	# I resample by taking p=ntf/nt as the best estimate of the proportion of female Bar/waitress that are transmitters
	# and then adding uncertainty around p based on p(1-p)/n
	bs	<- 1e4
	a	<- round(nt*rnorm(bs, mean=ntf/nt, sd=sqrt( (ntf/nt) * (1-ntf/nt) / nt )))
	# 'b' is exposed non-cases (exposed=Bar/waitress, non-case=recipients)
	b	<- round(nr*rnorm(bs, mean=nrf/nr, sd=sqrt( (nrf/nr) * (1-nrf/nr) / nr )))
	c	<- nt-a
	d	<- nr-b
	tmp2<- quantile( (a/c) / (b/d), p=c(0.025,0.975))
	a	<- ntf
	b	<- nrf
	c	<- nt-a
	d	<- nr-b
	tmp2<- c( (a/c) / (b/d), tmp2)
	
	bs	<- 1e4
	tmp2<- quantile( rnorm(bs, mean=ntf/nt, sd=sqrt( (ntf/nt) * (1-ntf/nt) / nt )) / rnorm(bs, mean=ndf/nd, sd=sqrt( (ndf/nd) * (1-ndf/nd) / nd )), p=c(0.025,0.975))
	tmp2<- c((ntf/nt) / (ndf/nd),tmp2)
	tmp3<- quantile( rnorm(bs, mean=nrf/nr, sd=sqrt( (nrf/nr) * (1-nrf/nr) / nr )) / rnorm(bs, mean=ndf/nd, sd=sqrt( (ndf/nd) * (1-ndf/nd) / nd )), p=c(0.025,0.975))
	tmp3<- c((nrf/nr) / (ndf/nd),tmp3)
	
	
	
	ntf	<- nrow(subset(tmp, PHSC_DIR=='m->f' & MALE_FACTOR=='Fishing' & MALE_COMM_TYPE=='fisherfolk'))
	nrf	<- nrow(subset(tmp, PHSC_DIR=='f->m' & MALE_FACTOR=='Fishing' & MALE_COMM_TYPE=='fisherfolk'))
	ndf	<- nrow(subset(tmp, PHSC_DIR=='denominator' & MALE_FACTOR=='Fishing' & MALE_COMM_TYPE=='fisherfolk'))
	nt	<- nrow(subset(tmp, PHSC_DIR=='m->f' & MALE_COMM_TYPE=='fisherfolk'))
	nr	<- nrow(subset(tmp, PHSC_DIR=='f->m' & MALE_COMM_TYPE=='fisherfolk'))
	nd	<- nrow(subset(tmp, PHSC_DIR=='denominator' & MALE_COMM_TYPE=='fisherfolk'))
	
	bs	<- 1e4
	tmp	<- quantile( rnorm(bs, mean=ntf/nt, sd=sqrt( (ntf/nt) * (1-ntf/nt) / nt )) / rnorm(bs, mean=ndf/nd, sd=sqrt( (ndf/nd) * (1-ndf/nd) / nd )), p=c(0.025,0.975))
	tmp	<- c((ntf/nt) / (ndf/nd),tmp)
	tmp	<- quantile( rnorm(bs, mean=nrf/nr, sd=sqrt( (nrf/nr) * (1-nrf/nr) / nr )) / rnorm(bs, mean=ndf/nd, sd=sqrt( (ndf/nd) * (1-ndf/nd) / nd )), p=c(0.025,0.975))
	tmp	<- c((nrf/nr) / (ndf/nd),tmp)
	#
	#	Age group
	#	
	tmp			<- subset(rtr, FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_birthdate.pdf',sep=''))
	setnames(tmp, c('FEMALE_BIRTHDATE','MALE_BIRTHDATE'), c('FEMALE_FACTOR','MALE_FACTOR'))	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'birth date', bw=4, dotsize=0.5, w=10, h=7)
	tmp			<- subset(rtr, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_birthdate_stablecouples.pdf',sep=''))
	setnames(tmp, c('FEMALE_BIRTHDATE','MALE_BIRTHDATE'), c('FEMALE_FACTOR','MALE_FACTOR'))	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'birth date', bw=4, dotsize=0.5, w=10, h=7)
	tmp			<- subset(rtr, PAIR_TYPE=='not registered as couple' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_birthdate_nocouples.pdf',sep=''))
	setnames(tmp, c('FEMALE_BIRTHDATE','MALE_BIRTHDATE'), c('FEMALE_FACTOR','MALE_FACTOR'))	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'birth date', bw=4, dotsize=0.5, w=10, h=7)
	#	
	tmp			<- subset(rtr, FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_agediff.pdf',sep=''))
	setnames(tmp, c('AGEDIFF'), c('MALE_FACTOR'))
	tmp[, FEMALE_FACTOR:= MALE_FACTOR]	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'age difference\n(years female younger)', bw=4, dotsize=0.5, w=10, h=7)
	tmp			<- subset(rtr, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_agediff_stablecouples.pdf',sep=''))
	setnames(tmp, c('AGEDIFF'), c('MALE_FACTOR'))
	tmp[, FEMALE_FACTOR:= MALE_FACTOR]	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'age difference\n(years female younger)', bw=4, dotsize=0.5, w=10, h=7)
	tmp			<- subset(rtr, PAIR_TYPE=='not registered as couple' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_agediff_nocouples.pdf',sep=''))
	setnames(tmp, c('AGEDIFF'), c('MALE_FACTOR'))
	tmp[, FEMALE_FACTOR:= MALE_FACTOR]	
	Rakai.plot.directed.pairs.continuous1(tmp, outfile, 'age difference\n(years female younger)', bw=4, dotsize=0.5, w=10, h=7)
	
	#
	#	Marriage Status
	#
	tmp			<- subset(rtr, FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_marriagestatus.pdf',sep=''))
	tmp			<- subset(rtr, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_marriagestatus_stablecouples.pdf',sep=''))
	tmp			<- subset(rtr, !grepl('cohabiting',PAIR_TYPE) & FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	outfile		<- file.path(dir, paste(run,'-phsc-directionpairs_marriagestatus_nocouples.pdf',sep=''))
	
	set(tmp, NULL, 'MALE_MARSTAT', tmp[, gsub('Never Married \\+ casual partner','Never Married',gsub('Previously Married \\+ casual partner','Previously Married',MALE_MARSTAT))])
	set(tmp, NULL, 'FEMALE_MARSTAT', tmp[, gsub('Never Married \\+ casual partner','Never Married',gsub('Previously Married \\+ casual partner','Previously Married',FEMALE_MARSTAT))])
	tmp2		<- unique(subset(ra, SEX=='F' & !is.na(FIRSTPOSDATE) & COMM_TYPE!='trading', c(RID, MARSTAT, COMM_TYPE)))
	setnames(tmp2, colnames(tmp2), paste0('FEMALE_',colnames(tmp2)))
	tmp2[, PHSC_DIR:='denominator']
	tmp			<- rbind(tmp, tmp2, fill=TRUE)
	tmp2		<- unique(subset(ra, SEX=='M' & !is.na(FIRSTPOSDATE) & COMM_TYPE!='trading', c(RID, MARSTAT, COMM_TYPE)))
	setnames(tmp2, colnames(tmp2), paste0('MALE_',colnames(tmp2)))
	tmp2[, PHSC_DIR:='denominator']
	tmp			<- rbind(tmp, tmp2, fill=TRUE)	
	setnames(tmp, c('FEMALE_MARSTAT','MALE_MARSTAT'), c('FEMALE_FACTOR','MALE_FACTOR'))
	cols		<- tmp[, length(unique(c(FEMALE_FACTOR, MALE_FACTOR)))]
	cols		<- colorRampPalette(brewer.pal(min(8,cols), "Set2"))( cols )	
	#cols		<- rainbow_hcl(tmp[, length(unique(c(FEMALE_FACTOR, MALE_FACTOR)))], start = 20, end = 340, c=100, l=60)
	names(cols)	<- tmp[, sort(unique(c(FEMALE_FACTOR, MALE_FACTOR)))]
	Rakai.plot.directed.pairs.discrete1(tmp, cols, outfile, 'marital &\nself-reported\nnon-marital\nrelationships,\n', w=10, h=7)
	
	
	
	#
	#	Education
	#
	tmp			<- subset(rtr, FEMALE_COMM_TYPE!='trading' & MALE_COMM_TYPE!='trading')
	setnames(tmp, c('FEMALE_EDUCAT','MALE_EDUCAT'), c('FEMALE_FACTOR','MALE_FACTOR'))
	cols		<- tmp[, length(unique(c(FEMALE_FACTOR, MALE_FACTOR)))]
	cols		<- colorRampPalette(brewer.pal(min(11,cols), "Set1"))( cols )	
	#cols		<- rainbow_hcl(tmp[, length(unique(c(FEMALE_FACTOR, MALE_FACTOR)))], start = 20, end = 340, c=100, l=60)
	names(cols)	<- tmp[, sort(unique(c(FEMALE_FACTOR, MALE_FACTOR)))]
	Rakai.plot.directed.pairs.discrete1(tmp, cols, file.path(dir, paste(run,'-phsc-directionpairs_education.pdf',sep='')), 'Educational status', w=10, h=7)	
}

RakaiAll.analyze.pairs.170405.comparetoprevious<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	#	load rd, rh, rs, rp, rpw, rplkl
	load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170405/RCCS_170405_w250_trmp_allpairs_posteriors_cmptoprv.rda')	
	#
	#	select likely pairs and compare between two runs
	#
	rps		<- unique(subset(rpw, select=c(FEMALE_SANGER_ID, MALE_SANGER_ID, MALE_RID, FEMALE_RID, PTY_RUN)))
	rps[, COUPID:= paste0(MALE_RID,':',FEMALE_RID)]
	tmp		<- unique(subset(rc, !is.na(FEMALE_SANGER_ID) & !is.na(MALE_SANGER_ID), select=c(COUPID, FEMALE_SANGER_ID, MALE_SANGER_ID, MALE_HH_NUM, FEMALE_HH_NUM, COUP_SC, PAIR_TYPE)))
	setnames(tmp, 'COUP_SC', 'COUP_TYPE')
	rps		<- merge(rps, tmp, by=c('COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID'), all.x=1)
	set(rps, rps[, which(is.na(PAIR_TYPE))], 'PAIR_TYPE', 'not registered as couple')
	rps		<- merge(rps, rplkl, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN'))
	confidence.cut	<- 0.5
	rps[, LIKELY_PAIR:= as.numeric(pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confidence.cut)]
	#	differences without including other+close as evidence for likely pair	
	rpd		<- subset(rps, GROUP%in%c('TYPE_PAIR_TODI') & TYPE=='likely pair')	
	tmp		<- rpd[, list(LIKELY_PAIR_ALWAYS=prod(LIKELY_PAIR)), by=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN','PAR_PRIOR')]
	rpd		<- merge(rpd, tmp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN','PAR_PRIOR'))	
	rpd[, list(LIKELY_PAIR_N=sum(LIKELY_PAIR)), by=c('RUN','PAR_PRIOR')]	
	#                                             RUN PAR_PRIOR LIKELY_PAIR_N
	#1:              RCCS_170320_w250_d50_st20_mr20_mt1_cl2_d5       0.1           199
	#2: RCCS_170405_w250_d50_st20_trB_blNormed_mr20_mt1_cl2_d5       0.1           194
	#2: RCCS_170405_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl2_d5  0.1           205
	#3:          RCCS_170405_w250_d50_st20_trB_mr20_mt1_cl2_d5       0.1           196
	#4:          RCCS_170405_w250_d50_st20_trC_mr20_mt1_cl2_d5       0.1           200
	#5:          RCCS_170405_w250_d50_st20_trU_mr20_mt1_cl2_d5       0.1           189
	#6:              RCCS_170320_w250_d50_st20_mr20_mt1_cl2_d5       3.0           158
	#7: RCCS_170405_w250_d50_st20_trB_blNormed_mr20_mt1_cl2_d5       3.0           167
	#8: RCCS_170405_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl2_d5  3.0           174
	#8:          RCCS_170405_w250_d50_st20_trB_mr20_mt1_cl2_d5       3.0           166
	#9:          RCCS_170405_w250_d50_st20_trC_mr20_mt1_cl2_d5       3.0           173
	#10:          RCCS_170405_w250_d50_st20_trU_mr20_mt1_cl2_d5      3.0           164
	#	--> with tipRule=c, we obtain most transmission pairs
	
	rpd		<- subset(rps, GROUP%in%c('TYPE_PAIR_TODI_NEW') & TYPE=='likely pair')	
	tmp		<- rpd[, list(LIKELY_PAIR_ALWAYS=prod(LIKELY_PAIR)), by=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN','PAR_PRIOR')]
	rpd		<- merge(rpd, tmp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN','PAR_PRIOR'))	
	rpd[, list(LIKELY_PAIR_N=sum(LIKELY_PAIR)), by=c('RUN','PAR_PRIOR')]	
	#	                                                       RUN PAR_PRIOR LIKELY_PAIR_N
	#1:              RCCS_170320_w250_d50_st20_mr20_mt1_cl2_d5       0.1           199
	#2: RCCS_170405_w250_d50_st20_trB_blNormed_mr20_mt1_cl2_d5       0.1           216
	#3:          RCCS_170405_w250_d50_st20_trB_mr20_mt1_cl2_d5       0.1           214
	#4:          RCCS_170405_w250_d50_st20_trC_mr20_mt1_cl2_d5       0.1           215
	#5:          RCCS_170405_w250_d50_st20_trU_mr20_mt1_cl2_d5       0.1           208
	#6:              RCCS_170320_w250_d50_st20_mr20_mt1_cl2_d5       3.0           158
	#7: RCCS_170405_w250_d50_st20_trB_blNormed_mr20_mt1_cl2_d5       3.0           185
	#8:          RCCS_170405_w250_d50_st20_trB_mr20_mt1_cl2_d5       3.0           185
	#9:          RCCS_170405_w250_d50_st20_trC_mr20_mt1_cl2_d5       3.0           187
	#10:          RCCS_170405_w250_d50_st20_trU_mr20_mt1_cl2_d5      3.0           181
	#	--> differences are a bit smaller; 9 cases in C that are not in U
	#									   4 cases in U that are not in C

	#	
	#	plot difference between c and u
	#	NOTE: this is using 'TYPE_PAIR_TODI'
	tmp		<- setdiff( 	subset(rpd, RUN=='RCCS_170405_w250_d50_st20_trU_mr20_mt1_cl2_d5' & LIKELY_PAIR==1 & PAR_PRIOR==3.0)[, LABEL_SH],
							subset(rpd, RUN=='RCCS_170405_w250_d50_st20_trC_mr20_mt1_cl2_d5' & LIKELY_PAIR==1 & PAR_PRIOR==3.0)[, LABEL_SH]	)

	tmp		<- setdiff( 	subset(rpd, RUN=='RCCS_170405_w250_d50_st20_trC_mr20_mt1_cl2_d5' & LIKELY_PAIR==1 & PAR_PRIOR==3.0)[, LABEL_SH],
							subset(rpd, RUN=='RCCS_170405_w250_d50_st20_trU_mr20_mt1_cl2_d5' & LIKELY_PAIR==1 & PAR_PRIOR==3.0)[, LABEL_SH]	)
	rps		<- unique(subset(merge(data.table(LABEL_SH=tmp), rpd, by='LABEL_SH'), select=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN')))
	write.csv(rps, file=paste0(outfile.base,'_summary.csv'))
	#	make detailed plots for selection
	outfile.base<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170405/with_TiesRuleC_butnot_TiesRuleU_'
	group		<- 'TYPE_DIR_TODI7x3'
	run			<- 'RCCS_170405_w250_d50_st20_trC_mr20_mt1_cl2_d5'	
	plot.select	<- unique(subset(merge(rplkl, rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP==group), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rpw2		<- subset(rpw, RUN==run & GROUP==group)
	rplkl2		<- subset(rplkl, RUN==run & GROUP==group & PAR_PRIOR==3)
	plot.file	<- paste0(outfile.base,'_plot_for_tiesRuleC_',group,'.pdf')	
	phsc.plot.windowsummaries.for.pairs(plot.select, rpw2, rplkl2, plot.file, cols=NULL, group=group)	
	run			<- 'RCCS_170405_w250_d50_st20_trB_blNormed_mr20_mt1_cl2_d5'	
	plot.select	<- unique(subset(merge(rplkl, rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP==group), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rpw2		<- subset(rpw, RUN==run & GROUP==group)
	rplkl2		<- subset(rplkl, RUN==run & GROUP==group & PAR_PRIOR==3)
	plot.file	<- paste0(outfile.base,'_plot_for_tiesRuleBNormed_',group,'.pdf')	
	phsc.plot.windowsummaries.for.pairs(plot.select, rpw2, rplkl2, plot.file, cols=NULL, group=group)
	run			<- 'RCCS_170405_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl2_d5'	
	plot.select	<- unique(subset(merge(rplkl, rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP==group), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rpw2		<- subset(rpw, RUN==run & GROUP==group)
	rplkl2		<- subset(rplkl, RUN==run & GROUP==group & PAR_PRIOR==3)
	plot.file	<- paste0(outfile.base,'_plot_for_tiesRuleBNormedOnFly_',group,'.pdf')	
	phsc.plot.windowsummaries.for.pairs(plot.select, rpw2, rplkl2, plot.file, cols=NULL, group=group)		
	run			<- 'RCCS_170405_w250_d50_st20_trU_mr20_mt1_cl2_d5'	
	plot.select	<- unique(subset(merge(rplkl, rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP==group), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rpw2		<- subset(rpw, RUN==run & GROUP==group)
	rplkl2		<- subset(rplkl, RUN==run & GROUP==group & PAR_PRIOR==3)
	plot.file	<- paste0(outfile.base,'_plot_for_tiesRuleU_',group,'.pdf')	
	#	plot distances
	tmp			<- subset(rpw, MALE_SANGER_ID=='15878_1_51' & FEMALE_SANGER_ID=='15677_1_49' & PTY_RUN==48 & GROUP=='TYPE_PAIR_TO', c(W_FROM, RUN, PATRISTIC_DISTANCE))
	ggplot(tmp, aes(x=W_FROM, y=PATRISTIC_DISTANCE, colour=RUN)) + geom_line()
	ggsave(file=paste(outfile.base,'_plot_for_distances.pdf'), w=25, h=10)
	phsc.plot.windowsummaries.for.pairs(plot.select, rpw2, rplkl2, plot.file, cols=NULL, group=group)
	#	plot phylogenies for pairs with little evidence in either of the two runs	
	tmp			<- copy(rps)
	set(tmp, NULL, 'FEMALE_SANGER_ID', tmp[, as.character(FEMALE_SANGER_ID)])
	set(tmp, NULL, 'MALE_SANGER_ID', tmp[, as.character(MALE_SANGER_ID)])
	run			<- 'RCCS_170405_w250_d50_st20_trC_mr20_mt1_cl2_d5'
	run			<- 'RCCS_170405_w250_d50_st20_trU_mr20_mt1_cl2_d5'
	rpw2		<- unique(subset(rpw, RUN==run, select=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM','W_TO','PATHS_12','PATHS_21','PATRISTIC_DISTANCE','CONTIGUOUS','TYPE_RAW')))
	set(rpw2, NULL, 'FEMALE_SANGER_ID', rpw2[, as.character(FEMALE_SANGER_ID)])
	set(rpw2, NULL, 'MALE_SANGER_ID', rpw2[, as.character(MALE_SANGER_ID)])
	set(rpw2, NULL, 'CONTIGUOUS', rpw2[, as.integer(CONTIGUOUS)])
	set(rpw2, NULL, 'PATRISTIC_DISTANCE', rpw2[, round(PATRISTIC_DISTANCE, d=4)])
	
	load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170405/RCCS_170405_w250_d50_st20_trC_mr20_mt1_cl2_d5_phscout.rda')
	load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170405/RCCS_170405_w250_d50_st20_trU_mr20_mt1_cl2_d5_phscout.rda')
	outfile.base<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170405/with_TiesRuleC_butnot_TiesRuleU_phU_'
	invisible(sapply(seq_len(nrow(tmp)), function(ii)
					{	
						#ii<- 1
						ids			<- c(tmp[ii, MALE_SANGER_ID],tmp[ii, FEMALE_SANGER_ID])
						pty.run		<- tmp[ii, PTY_RUN]
						dfs			<- subset(dtrees, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
						dfs[, MALE_SANGER_ID:=ids[1]]
						dfs[, FEMALE_SANGER_ID:=ids[2]]
						dfs			<- merge(dfs, rpw2, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM','W_TO'), all.x=TRUE)
						dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', pty.run, '\nwindow ', W_FROM,'-', W_TO,'\n',PATHS_12,' ',PATHS_21, ' ',CONTIGUOUS,' ',TYPE_RAW, '\n', PATRISTIC_DISTANCE, sep='')]]								
						plot.file	<- paste0(outfile.base, pty.run,'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
						invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10))					
					}))	
	
	
	tmp		<- setdiff( 	subset(rpd, RUN=='RCCS_170405_w250_d50_st20_trB_blNormed_mr20_mt1_cl2_d5' & LIKELY_PAIR==1 & PAR_PRIOR==3.0)[, LABEL_SH],
							subset(rpd, RUN=='RCCS_170405_w250_d50_st20_trB_mr20_mt1_cl2_d5' & LIKELY_PAIR==1 & PAR_PRIOR==3.0)[, LABEL_SH]	)
	
	
	 
}

RakaiAll.analyze.pairs.170410.comparetoprevious<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	#	load rd, rh, rs, rp, rpw, rplkl
	load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170410/RCCS_170410_w250_trmp_allpairs_posteriors_cmptoprv.rda')	
	#
	#	select likely pairs and compare between two runs
	#
	rps		<- unique(subset(rpw, select=c(FEMALE_SANGER_ID, MALE_SANGER_ID, MALE_RID, FEMALE_RID, PTY_RUN)))
	rps[, COUPID:= paste0(MALE_RID,':',FEMALE_RID)]
	tmp		<- unique(subset(rc, !is.na(FEMALE_SANGER_ID) & !is.na(MALE_SANGER_ID), select=c(COUPID, FEMALE_SANGER_ID, MALE_SANGER_ID, MALE_HH_NUM, FEMALE_HH_NUM, COUP_SC, PAIR_TYPE)))
	setnames(tmp, 'COUP_SC', 'COUP_TYPE')
	rps		<- merge(rps, tmp, by=c('COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID'), all.x=1)
	set(rps, rps[, which(is.na(PAIR_TYPE))], 'PAIR_TYPE', 'not registered as couple')
	rps		<- merge(rps, rplkl, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN'))
	confidence.cut	<- 0.5
	rps[, LIKELY_PAIR:= as.numeric(pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confidence.cut)]
	#	differences without including other+close as evidence for likely pair		
	rpd		<- subset(rps, GROUP%in%c('TYPE_PAIR_TODI_NEW') & TYPE=='likely pair')	
	tmp		<- rpd[, list(LIKELY_PAIR_ALWAYS=prod(LIKELY_PAIR)), by=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN','PAR_PRIOR')]
	rpd		<- merge(rpd, tmp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN','PAR_PRIOR'))	
	rpd[, list(LIKELY_PAIR_N=sum(LIKELY_PAIR)), by=c('RUN','PAR_PRIOR')]	
	
	#	                                                             RUN PAR_PRIOR LIKELY_PAIR_N
	#1: RCCS_170410_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl3.5_d8       0.1           255
	#2:               RCCS_170410_w250_d50_st20_trB_mr20_mt1_cl3.5_d8       0.1           245
	#3:               RCCS_170410_w250_d50_st20_trC_mr20_mt1_cl3.5_d8       0.1           242
	#4:               RCCS_170410_w250_d50_st20_trU_mr20_mt1_cl3.5_d8       0.1           244
	#5: RCCS_170410_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl3.5_d8       3.0           232
	#6:               RCCS_170410_w250_d50_st20_trB_mr20_mt1_cl3.5_d8       3.0           216
	#7:               RCCS_170410_w250_d50_st20_trC_mr20_mt1_cl3.5_d8       3.0           214
	#8:               RCCS_170410_w250_d50_st20_trU_mr20_mt1_cl3.5_d8       3.0           218
	
	#	OK this is looking good
	#	for RCCS_170410_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl3.5_d8, plot "edge cases"
	
	
	tmp		<- setdiff( 	subset(rpd, RUN=='RCCS_170410_w250_d50_st20_trC_mr20_mt1_cl3.5_d8' & LIKELY_PAIR==1 & PAR_PRIOR==3.0)[, LABEL_SH],
			subset(rpd, RUN=='RCCS_170410_w250_d50_st20_trU_mr20_mt1_cl3.5_d8' & LIKELY_PAIR==1 & PAR_PRIOR==3.0)[, LABEL_SH]	)
	
	
	#	--> differences are a bit smaller; 9 cases in C that are not in U
	#									   4 cases in U that are not in C
	
	#	
	#	plot difference between c and u
	#	NOTE: this is using 'TYPE_PAIR_TODI'
	tmp		<- setdiff( 	subset(rpd, RUN=='RCCS_170405_w250_d50_st20_trU_mr20_mt1_cl2_d5' & LIKELY_PAIR==1 & PAR_PRIOR==3.0)[, LABEL_SH],
			subset(rpd, RUN=='RCCS_170405_w250_d50_st20_trC_mr20_mt1_cl2_d5' & LIKELY_PAIR==1 & PAR_PRIOR==3.0)[, LABEL_SH]	)
	
	tmp		<- setdiff( 	subset(rpd, RUN=='RCCS_170405_w250_d50_st20_trC_mr20_mt1_cl2_d5' & LIKELY_PAIR==1 & PAR_PRIOR==3.0)[, LABEL_SH],
			subset(rpd, RUN=='RCCS_170405_w250_d50_st20_trU_mr20_mt1_cl2_d5' & LIKELY_PAIR==1 & PAR_PRIOR==3.0)[, LABEL_SH]	)
	rps		<- unique(subset(merge(data.table(LABEL_SH=tmp), rpd, by='LABEL_SH'), select=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN')))
	write.csv(rps, file=paste0(outfile.base,'_summary.csv'))
	#	make detailed plots for selection
	outfile.base<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170405/with_TiesRuleC_butnot_TiesRuleU_'
	group		<- 'TYPE_DIR_TODI7x3'
	run			<- 'RCCS_170405_w250_d50_st20_trC_mr20_mt1_cl2_d5'	
	plot.select	<- unique(subset(merge(rplkl, rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP==group), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rpw2		<- subset(rpw, RUN==run & GROUP==group)
	rplkl2		<- subset(rplkl, RUN==run & GROUP==group & PAR_PRIOR==3)
	plot.file	<- paste0(outfile.base,'_plot_for_tiesRuleC_',group,'.pdf')	
	phsc.plot.windowsummaries.for.pairs(plot.select, rpw2, rplkl2, plot.file, cols=NULL, group=group)	
	run			<- 'RCCS_170405_w250_d50_st20_trB_blNormed_mr20_mt1_cl2_d5'	
	plot.select	<- unique(subset(merge(rplkl, rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP==group), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rpw2		<- subset(rpw, RUN==run & GROUP==group)
	rplkl2		<- subset(rplkl, RUN==run & GROUP==group & PAR_PRIOR==3)
	plot.file	<- paste0(outfile.base,'_plot_for_tiesRuleBNormed_',group,'.pdf')	
	phsc.plot.windowsummaries.for.pairs(plot.select, rpw2, rplkl2, plot.file, cols=NULL, group=group)
	run			<- 'RCCS_170405_w250_d50_st20_trB_blNormedOnFly_mr20_mt1_cl2_d5'	
	plot.select	<- unique(subset(merge(rplkl, rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP==group), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rpw2		<- subset(rpw, RUN==run & GROUP==group)
	rplkl2		<- subset(rplkl, RUN==run & GROUP==group & PAR_PRIOR==3)
	plot.file	<- paste0(outfile.base,'_plot_for_tiesRuleBNormedOnFly_',group,'.pdf')	
	phsc.plot.windowsummaries.for.pairs(plot.select, rpw2, rplkl2, plot.file, cols=NULL, group=group)		
	run			<- 'RCCS_170405_w250_d50_st20_trU_mr20_mt1_cl2_d5'	
	plot.select	<- unique(subset(merge(rplkl, rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP==group), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rpw2		<- subset(rpw, RUN==run & GROUP==group)
	rplkl2		<- subset(rplkl, RUN==run & GROUP==group & PAR_PRIOR==3)
	plot.file	<- paste0(outfile.base,'_plot_for_tiesRuleU_',group,'.pdf')	
	#	plot distances
	tmp			<- subset(rpw, MALE_SANGER_ID=='15878_1_51' & FEMALE_SANGER_ID=='15677_1_49' & PTY_RUN==48 & GROUP=='TYPE_PAIR_TO', c(W_FROM, RUN, PATRISTIC_DISTANCE))
	ggplot(tmp, aes(x=W_FROM, y=PATRISTIC_DISTANCE, colour=RUN)) + geom_line()
	ggsave(file=paste(outfile.base,'_plot_for_distances.pdf'), w=25, h=10)
	phsc.plot.windowsummaries.for.pairs(plot.select, rpw2, rplkl2, plot.file, cols=NULL, group=group)
	#	plot phylogenies for pairs with little evidence in either of the two runs	
	tmp			<- copy(rps)
	set(tmp, NULL, 'FEMALE_SANGER_ID', tmp[, as.character(FEMALE_SANGER_ID)])
	set(tmp, NULL, 'MALE_SANGER_ID', tmp[, as.character(MALE_SANGER_ID)])
	run			<- 'RCCS_170405_w250_d50_st20_trC_mr20_mt1_cl2_d5'
	run			<- 'RCCS_170405_w250_d50_st20_trU_mr20_mt1_cl2_d5'
	rpw2		<- unique(subset(rpw, RUN==run, select=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM','W_TO','PATHS_12','PATHS_21','PATRISTIC_DISTANCE','CONTIGUOUS','TYPE_RAW')))
	set(rpw2, NULL, 'FEMALE_SANGER_ID', rpw2[, as.character(FEMALE_SANGER_ID)])
	set(rpw2, NULL, 'MALE_SANGER_ID', rpw2[, as.character(MALE_SANGER_ID)])
	set(rpw2, NULL, 'CONTIGUOUS', rpw2[, as.integer(CONTIGUOUS)])
	set(rpw2, NULL, 'PATRISTIC_DISTANCE', rpw2[, round(PATRISTIC_DISTANCE, d=4)])
	
	load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170405/RCCS_170405_w250_d50_st20_trC_mr20_mt1_cl2_d5_phscout.rda')
	load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170405/RCCS_170405_w250_d50_st20_trU_mr20_mt1_cl2_d5_phscout.rda')
	outfile.base<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/170405/with_TiesRuleC_butnot_TiesRuleU_phU_'
	invisible(sapply(seq_len(nrow(tmp)), function(ii)
					{	
						#ii<- 1
						ids			<- c(tmp[ii, MALE_SANGER_ID],tmp[ii, FEMALE_SANGER_ID])
						pty.run		<- tmp[ii, PTY_RUN]
						dfs			<- subset(dtrees, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
						dfs[, MALE_SANGER_ID:=ids[1]]
						dfs[, FEMALE_SANGER_ID:=ids[2]]
						dfs			<- merge(dfs, rpw2, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM','W_TO'), all.x=TRUE)
						dfs[, TITLE:= dfs[, paste('male ', ids[1],'\nfemale ',ids[2],'\nrun ', pty.run, '\nwindow ', W_FROM,'-', W_TO,'\n',PATHS_12,' ',PATHS_21, ' ',CONTIGUOUS,' ',TYPE_RAW, '\n', PATRISTIC_DISTANCE, sep='')]]								
						plot.file	<- paste0(outfile.base, pty.run,'_M_',ids[1],'_F_', ids[2],'_collapsed.pdf')					
						invisible(phsc.plot.phycollapsed.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), drop.less.than.n.ids=2, plot.file=plot.file, pdf.h=10, pdf.rw=5, pdf.ntrees=20, pdf.title.size=10))					
					}))	
	
	
	tmp		<- setdiff( 	subset(rpd, RUN=='RCCS_170405_w250_d50_st20_trB_blNormed_mr20_mt1_cl2_d5' & LIKELY_PAIR==1 & PAR_PRIOR==3.0)[, LABEL_SH],
			subset(rpd, RUN=='RCCS_170405_w250_d50_st20_trB_mr20_mt1_cl2_d5' & LIKELY_PAIR==1 & PAR_PRIOR==3.0)[, LABEL_SH]	)
	
	
	
}

RakaiCouples.save.couples.to.rda<- function()
{
	require(data.table)
	wdir				<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples"
	#	get epi info
	tmp		<- RakaiCirc.epi.get.info.170208()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#	get sequence info and add to recipients
	load('~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/circumcision/RCCS_SeqInfo_170505.rda')		
	rs		<- subset(rs, !is.na(VISIT))	
	
	#	load combination of all couples sequences
	#infile	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda'
	#	load Kates couple data
	infile	<- '~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/Pangea_Couples.csv'
	rc		<- as.data.table(read.csv(infile, stringsAsFactors=FALSE))
	setnames(rc, c('female.PANGEA.ID','male.PANGEA.ID'), c('female.TAXA','male.TAXA'))
	setnames(rc, colnames(rc), gsub('\\.','_',toupper(colnames(rc))))
	setnames(rc, c('MALE_RCCS_STUDYID','FEMALE_RCCS_STUDYID'), c('MALE_RID','FEMALE_RID'))
	set(rc, NULL, 'MALE_DATE', rc[, hivc.db.Date2numeric(as.Date(MALE_DATE))])
	set(rc, NULL, 'FEMALE_DATE', rc[, hivc.db.Date2numeric(as.Date(FEMALE_DATE))])
	set(rc, NULL, 'MALE_LASTNEGDATE', rc[, hivc.db.Date2numeric(as.Date(MALE_LASTNEGDATE))])
	set(rc, NULL, 'FEMALE_LASTNEGDATE', rc[, hivc.db.Date2numeric(as.Date(FEMALE_LASTNEGDATE))])	
	set(rc, NULL, 'MALE_FIRSTPOSDATE', rc[, hivc.db.Date2numeric(as.Date(MALE_FIRSTPOSDATE))])
	set(rc, NULL, 'FEMALE_FIRSTPOSDATE', rc[, hivc.db.Date2numeric(as.Date(FEMALE_FIRSTPOSDATE))])
	#	update sequence info on visits
	set(rc, NULL, c('FEMALE_PANGEA_ID','FEMALE_SANGER_ID','FEMALE_TAXA','MALE_PANGEA_ID','MALE_SANGER_ID','MALE_TAXA'), NULL)
	rc		<- unique(rc)
	setnames(rs, colnames(rs), paste0('MALE_',colnames(rs)))
	setnames(rs, 'MALE_VISIT', 'VISIT')
	rc		<- merge(rc, rs, by=c('MALE_RID','VISIT'), all.x=1)
	setnames(rs, colnames(rs), gsub('MALE_','FEMALE_',colnames(rs)))
	rc		<- merge(rc, rs, by=c('FEMALE_RID','VISIT'), all.x=1)
	setnames(rc, c('MALE_SID','FEMALE_SID','MALE_PIDF','FEMALE_PIDF'), c('MALE_SANGER_ID','FEMALE_SANGER_ID','MALE_TAXA','FEMALE_TAXA'))
	set(rc, NULL, c('FEMALE_RCCS_SHIP_BATCH','MALE_RCCS_SHIP_BATCH','MALE_SEQDAT','FEMALE_SEQDAT'), NULL)
	setnames(rs, colnames(rs), gsub('FEMALE_','',colnames(rs)))
	rc[, length(unique(COUPID))]	#	486 couples
	#
	#	get all possible unique pairings of SANGER_IDs from couples
	#	422
	rp		<- unique(subset(rc, !is.na(MALE_SANGER_ID), c(COUPID, MALE_RID, MALE_SANGER_ID, MALE_TAXA )))
	tmp		<- unique(subset(rc, !is.na(FEMALE_SANGER_ID), c(COUPID, FEMALE_RID, FEMALE_SANGER_ID, FEMALE_TAXA )))	
	rp		<- merge(rp, tmp, by='COUPID')
	rp		<- unique(rp, by=c('FEMALE_SANGER_ID','MALE_SANGER_ID'))
	cat('\nnumber of pairs with sequence from both partners=',rp[, length(unique(COUPID))])	
	#
	#	add couples where sequences from both are not available
	#
	tmp		<- unique(subset(rc, select=c("COUPID", "FEMALE_RID", "MALE_RID", "MALE_SANGER_ID", "MALE_TAXA", "FEMALE_SANGER_ID", "FEMALE_TAXA"))) 
	tmp2	<- tmp[, list(	BOTH_SEQUENCED= any(!is.na(FEMALE_SANGER_ID)) & any(!is.na(MALE_SANGER_ID))), by='COUPID']
	tmp		<- unique(merge(tmp, subset(tmp2, !BOTH_SEQUENCED), by='COUPID'), by='COUPID')
	set(tmp, NULL, c('BOTH_SEQUENCED'), NULL)
	rp		<- rbind(rp, tmp)
	#
	#	add visit when both first positive
	#	did couples always / never live together before then?
	#
	tmp		<- unique(subset(rd, select=c(RID, VISIT, HH_NUM)))
	setnames(tmp, colnames(tmp)[colnames(tmp)!='VISIT'], paste('MALE_',colnames(tmp)[colnames(tmp)!='VISIT'],sep=''))
	rc		<- merge(rc, tmp, by=c('MALE_RID','VISIT'), all.x=1)
	setnames(tmp, colnames(tmp), gsub('MALE','FEMALE',colnames(tmp)))
	rc		<- merge(rc, tmp, by=c('FEMALE_RID','VISIT'), all.x=1)
	tmp		<- rc[, {
				z	<- which(MALE_FIRSTPOSDATE<=(MALE_DATE+.2) & FEMALE_FIRSTPOSDATE<=(FEMALE_DATE+.2))
				if(length(z)==0)
					z	<- which.max(VISIT)
				z	<- min(z)
				tmp	<- c( sum(na.omit(MALE_HH_NUM[1:z]==FEMALE_HH_NUM[1:z])), sum(na.omit(MALE_HH_NUM[1:z]!=FEMALE_HH_NUM[1:z])))				
				list(VISIT_FIRSTCONCPOS=VISIT[z], PAIR_TYPE=ifelse(tmp[1]>0 & tmp[2]==0, 'stable cohabiting', ifelse(tmp[2]>0, 'not always cohabiting', 'no joint household data')))
			}, by='COUPID']
	rp		<- merge(rp, tmp, by='COUPID')
	#
	#	add epi info closest to first +/+ positive
	#	for males
	tmp		<- subset(rd, !is.na(PID), select=c(RID, VISIT, SEX, REGION, COMM_NUM, HH_NUM, BIRTHDATE, LASTNEGVIS, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, RELIGION))
	tmp2	<- unique(subset(rp, select=c('MALE_RID','VISIT_FIRSTCONCPOS')))
	setnames(tmp2, 'MALE_RID', 'RID')
	tmp		<- merge(tmp, tmp2, by='RID')
	tmp2	<- tmp[, list(VISIT= VISIT[which.min(abs(VISIT_FIRSTCONCPOS-VISIT))]), by=c('RID','VISIT_FIRSTCONCPOS')]
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT_FIRSTCONCPOS','VISIT'))
	set(tmp, NULL, 'VISIT', NULL)	
	tmp		<- unique(tmp, by=c('RID','VISIT_FIRSTCONCPOS'))
	setnames(tmp, colnames(tmp)[colnames(tmp)!='VISIT_FIRSTCONCPOS'], paste('MALE_',colnames(tmp)[colnames(tmp)!='VISIT_FIRSTCONCPOS'],sep=''))
	rp		<- merge(rp,tmp,by=c('MALE_RID','VISIT_FIRSTCONCPOS'))
	#	for females
	tmp		<- subset(rd, !is.na(PID), select=c(RID, VISIT, SEX, REGION, COMM_NUM, HH_NUM, BIRTHDATE, LASTNEGVIS, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, RELIGION))
	tmp2	<- unique(subset(rp, select=c('FEMALE_RID','VISIT_FIRSTCONCPOS')))
	setnames(tmp2, 'FEMALE_RID', 'RID')
	tmp		<- merge(tmp, tmp2, by='RID')
	tmp2	<- tmp[, list(VISIT= VISIT[which.min(abs(VISIT_FIRSTCONCPOS-VISIT))]), by=c('RID','VISIT_FIRSTCONCPOS')]
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT_FIRSTCONCPOS','VISIT'))
	set(tmp, NULL, 'VISIT', NULL)	
	tmp		<- unique(tmp, by=c('RID','VISIT_FIRSTCONCPOS'))
	setnames(tmp, colnames(tmp)[colnames(tmp)!='VISIT_FIRSTCONCPOS'], paste('FEMALE_',colnames(tmp)[colnames(tmp)!='VISIT_FIRSTCONCPOS'],sep=''))
	rp		<- merge(rp,tmp,by=c('FEMALE_RID','VISIT_FIRSTCONCPOS'))
	#
	#	add RH info closest to first +/+ positive
	#	for males
	tmp		<- unique(rh)
	set(tmp, NULL, c('COMM_NUM','SEX'), NULL)
	tmp2	<- unique(subset(rp, select=c('MALE_RID','VISIT_FIRSTCONCPOS')))
	setnames(tmp2, 'MALE_RID', 'RID')
	tmp		<- merge(tmp, tmp2, by='RID')
	tmp2	<- tmp[, list(VISIT= VISIT[which.min(abs(VISIT_FIRSTCONCPOS-VISIT))]), by=c('RID','VISIT_FIRSTCONCPOS')]
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT_FIRSTCONCPOS','VISIT'))
	set(tmp, NULL, 'VISIT', NULL)	
	tmp		<- unique(tmp, by=c('RID','VISIT_FIRSTCONCPOS'))
	setnames(tmp, colnames(tmp)[colnames(tmp)!='VISIT_FIRSTCONCPOS'], paste('MALE_',colnames(tmp)[colnames(tmp)!='VISIT_FIRSTCONCPOS'],sep=''))
	rp		<- merge(rp,tmp,by=c('MALE_RID','VISIT_FIRSTCONCPOS'))
	#	for females
	tmp		<- unique(rh)
	set(tmp, NULL, c('COMM_NUM','SEX'), NULL)
	tmp2	<- unique(subset(rp, select=c('FEMALE_RID','VISIT_FIRSTCONCPOS')))
	setnames(tmp2, 'FEMALE_RID', 'RID')
	tmp		<- merge(tmp, tmp2, by='RID')
	tmp2	<- tmp[, list(VISIT= VISIT[which.min(abs(VISIT_FIRSTCONCPOS-VISIT))]), by=c('RID','VISIT_FIRSTCONCPOS')]
	tmp		<- merge(tmp, tmp2, by=c('RID','VISIT_FIRSTCONCPOS','VISIT'))
	set(tmp, NULL, 'VISIT', NULL)	
	tmp		<- unique(tmp, by=c('RID','VISIT_FIRSTCONCPOS'))
	setnames(tmp, colnames(tmp)[colnames(tmp)!='VISIT_FIRSTCONCPOS'], paste('FEMALE_',colnames(tmp)[colnames(tmp)!='VISIT_FIRSTCONCPOS'],sep=''))
	rp		<- merge(rp,tmp,by=c('FEMALE_RID','VISIT_FIRSTCONCPOS'))
	#
	#	get info on sequences: date of sampling
	#	for males	
	tmp		<- unique(subset(rs, !is.na(SID), select=c(SID, SAMPLE_DATE)), by='SID')
	setnames(tmp, c('SID','SAMPLE_DATE'), c('MALE_SANGER_ID','MALE_SEQDATE'))
	rp		<- merge(rp, tmp, by='MALE_SANGER_ID', all.x=1)
	tmp		<- unique(subset(rs, !is.na(SID), select=c(SID, SAMPLE_DATE)), by='SID')
	setnames(tmp, c('SID','SAMPLE_DATE'), c('FEMALE_SANGER_ID','FEMALE_SEQDATE'))	
	rp		<- merge(rp, tmp, by='FEMALE_SANGER_ID', all.x=1)	
	#
	#	define direction based on seroconversion dates
	#
	rp[, COUP_SC:='seropos']
	set(rp, rp[, which(MALE_LASTNEGDATE>=FEMALE_FIRSTPOSDATE)],'COUP_SC','F->M')
	set(rp, rp[, which(FEMALE_LASTNEGDATE>=MALE_FIRSTPOSDATE)],'COUP_SC','M->F')
	set(rp, rp[, which(FEMALE_FIRSTPOSDATE==MALE_FIRSTPOSDATE)],'COUP_SC','seroinc')
	rp[, table(COUP_SC)]
	#	F->M    M->F seroinc seropos 
	#	20      32      79     426 

	#	link individual household data from most recent visit
	set(rp, rp[, which(PAIR_TYPE=='no joint household data' & FEMALE_HH_NUM==MALE_HH_NUM)], 'PAIR_TYPE', 'stable cohabiting')
	set(rp, rp[, which(PAIR_TYPE=='no joint household data' & FEMALE_HH_NUM!=MALE_HH_NUM)], 'PAIR_TYPE', 'not always cohabiting')	
	#	save
	save(rp, file="~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v170505_info.rda")	
}

hivc.db.Date2numeric<- function( x )
{
	if(!class(x)%in%c('Date','character'))	return( x )
	x	<- as.POSIXlt(x)
	tmp	<- x$year + 1900
	x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
	x	
}
######################################################################################
project.Rakai.aliRegion1.597<- function()
{	
	require(big.phylo)
	require(plyr)
	
	#	load 597 new files
	infile	<- '~/Dropbox (SPH Imperial College)/PANGEA_alignments/PANGEA_UG full alignment/PANGEA_HIV_n597_Imperial_v160916_RakaiPlates.fasta'
	sqn		<- read.dna(infile, format='fa')
	sqni	<- data.table(TAXA=rownames(sqn))
	sqni	<- subset(sqni, grepl('HXB2|consensus',TAXA))
	sqni[, SID:= gsub('_consensus','',TAXA)]
	set(sqni, sqni[, which(grepl('HXB2',TAXA))],'SID',NA_character_)
	#merge(sqni, unique(subset(rs, !is.na(SID), select=c(RID,SEQTYPE,SID,PID))), by='SID',all.x=1)
	sqn		<- sqn[ sqni[,TAXA], ]
	
	#	load PANGEA alignment
	infile	<- '~/Dropbox (SPH Imperial College)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_HXB2.fasta'
	sqp		<- read.dna(infile,format='fa')
	
	#	required: HXB2 in alignment
	ans		<- seq.align.based.on.common.reference(sqn, sqp, return.common.sites=TRUE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	outfile	<- '~/Dropbox (SPH Imperial College)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_597_p_HXB2.fasta'
	write.dna(ans, file=outfile, format='fasta', colsep='', nbcol=-1)
	
	#	get genetic distance matrix
	sq			<- ans
	sq			<- as.character(sq)
	sq[sq=='?']	<- '-'
	sq			<- as.DNAbin(sq)
	sq.gd		<- dist.dna(sq, pairwise.deletion=TRUE)
	sq.gd		<- as.data.table(melt(as.matrix(sq.gd), varnames=c('TAXA1','TAXA2')))
	setnames(sq.gd, 'value', 'PD')
	sq.gd		<- subset(sq.gd, TAXA1!=TAXA2)
	set(sq.gd, NULL, 'TAXA1', sq.gd[, as.character(TAXA1)])
	set(sq.gd, NULL, 'TAXA2', sq.gd[, as.character(TAXA2)])		
	#	add overlap	
	setkey(sq.gd, TAXA1, TAXA2)
	tmp			<- as.character(sq)
	tmp			<- !( tmp=='-' | tmp=='?' | tmp=='n' )
	tmp[]		<- as.integer(tmp) 
	tmp2		<- sq.gd[, list(OVERLAP= sum(bitwAnd(tmp[TAXA1,], tmp[TAXA2,])) ), by=c('TAXA1','TAXA2')]	
	sq.gd		<- merge(sq.gd, tmp2, by=c('TAXA1','TAXA2'))
	save(sq, sq.gd, file=gsub('\\.fasta','\\.rda',outfile))
	
	#
	#	resolve SANGER IDs
	#
	load('~/Dropbox (SPH Imperial College)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_597_p_HXB2.rda')
	sqi		<- data.table(TAXA=rownames(sq), ID=seq_len(nrow(sq)))
	
	infile	<- '~/Dropbox (SPH Imperial College)/PANGEA_data/2016-09-18_PAN_SANGER_IDs.txt'
	pni		<- as.data.table(read.table(infile, header=TRUE, sep='\t', comment.char='%'))
	setnames(pni, colnames(pni), toupper(colnames(pni)))
	set(pni, NULL, 'SANGERID', pni[, gsub('#','_',SANGERID)])
	pni[, TAXA:= paste(SANGERID,'_consensus',sep='')]
	sqi		<- merge(sqi, pni, by='TAXA',all.x=1)
	stopifnot( nrow(subset(sqi, grepl('consensus',TAXA) & is.na(SANGERID)))==0 )
	#	set new taxa names
	sqi[, TAXA_NEW:=TAXA]
	tmp		<- sqi[, which(!is.na(PANGEA_ID))]
	set(sqi, tmp, 'TAXA_NEW', sqi[tmp, PANGEA_ID])
	setkey(sqi, ID)
	rownames(sq)	<- sqi[, TAXA_NEW]
	setnames(sqi, c('TAXA','TAXA_NEW'), c('TAXA1','TAXA1_NEW'))
	sq.gd			<- merge(sq.gd, subset(sqi, select=c(TAXA1,TAXA1_NEW)), by='TAXA1')
	setnames(sqi, c('TAXA1','TAXA1_NEW'), c('TAXA2','TAXA2_NEW'))
	sq.gd			<- merge(sq.gd, subset(sqi, select=c(TAXA2,TAXA2_NEW)), by='TAXA2')
	set(sq.gd, NULL, c('TAXA1','TAXA2'), NULL)
	setnames(sq.gd, c('TAXA1_NEW','TAXA2_NEW'), c('TAXA1','TAXA2'))
	#	write to file
	outfile	<- '~/Dropbox (SPH Imperial College)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_597_p_HXB2.fasta'
	write.dna(sq, file=outfile, format='fasta', colsep='', nbcol=-1)	
	#
	#	read COMETv0.5
	#
	infile	<- '~/Dropbox (SPH Imperial College)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_597_p_HXB2_COMETv0.5.txt'
	sqi		<- as.data.table(read.table(infile, header=TRUE, sep='\t',stringsAsFactors=FALSE))
	setnames(sqi, c('name','subtype','length'), c('TAXA','COMET_ST','COMET_N'))
	sqi[, COMET_V:='0.5']
	#
	#	read COMETv2.1
	#
	infile	<- '~/Dropbox (SPH Imperial College)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_597_p_HXB2_COMETv2.1.txt'
	tmp		<- as.data.table(read.table(infile, header=TRUE, sep='\t',stringsAsFactors=FALSE))
	setnames(tmp, c('name','subtype','bootstrap.support'), c('TAXA','COMET_ST','COMET_BOOTSTRAP'))
	tmp[, COMET_V:='2.1']
	sqi		<- rbind(sqi, tmp, use.names=TRUE, fill=TRUE)
	#	set to factors to save a bit of space
	set(sq.gd, NULL, 'TAXA1', sq.gd[,factor(TAXA1)])
	set(sq.gd, NULL, 'TAXA2', sq.gd[,factor(TAXA2)])
	save(sq, sqi, sq.gd, file=gsub('\\.fasta','\\.rda',outfile))
}
######################################################################################
######################################################################################
project.Rakai.aliRegion1.add.HXB2.RCCSmissing<- function()
{
	wdir	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/circumcision"
	#	load Susanna s codon alignment
	susa.f	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/Susannah3/Region1_codon aligned.fasta'
	susa.s	<- read.dna(susa.f, format='fasta')
	susa.d	<- data.table(	ID=gsub('\\*.*','',rownames(susa.s)),
							DATA= factor(grepl('^PG', rownames(susa.s)), levels=c(TRUE,FALSE), label=c('PNG','LNL')))
	infile.relabel	<- "~/Dropbox (SPH Imperial College)/PANGEA_alignments/Regional phylogenetic analyses/Region 1 gag analysis/SummaryofGagSequenceData.rda"				
	load(infile.relabel)
	s		<- as.data.table(summaryData)
	s[, TAXA:= gsub('a\\||b\\|','',seqid)]		#taxa names have been tinkered with. does not merge.
	s[, ID:=gsub('\\*.*','',TAXA)]	
	stopifnot( nrow(subset(s, grepl('\\|',TAXA)))==0 )
	stopifnot( nrow(subset(s, is.na(ID)))==0 )
	#
	#	see what s missing
	#
	ssu		<- merge(s, susa.d, by='ID', all=1)
	cat('\nseqs in codon alignment that are not in summaryData', nrow(subset(ssu, is.na(seqid))))
	#	0, yay
	cat('\nseqs in summaryData that are not in codon alignment', nrow(subset(ssu, is.na(DATA))))
	#	3094 
	#	subset(ssu, is.na(DATA))[, table(dataSource, useNA='if')]
	#	select missing RCCS historical and HXB2 
	tmp		<- subset(ssu, is.na(DATA) & dataSource=='RCCS historical' )
	tmp		<- rbind(tmp,subset(ssu, is.na(DATA) & dataSource=='reference' & grepl('HXB2',TAXA)))
	tmp[, relabel2:= gsub('^[^\\*]+\\*','',relabel)]
	# 
	#	get missing sequences
	#
	infile.gag	<- '~/Dropbox (SPH Imperial College)/PANGEA_alignments/Rakai Data for IqTree/gag.sqn.fasta'
	seq.m		<- read.dna(infile.gag, format='fasta')
	seq.m		<- seq.m[ tmp[, relabel2], ]
	file.seq.m	<- file.path('/Users/Oliver/duke/tmp', 'tmp_missingseq_150825.fasta')
	write.dna(seq.m, file=file.seq.m, format='fasta')
	#
	#	add missing sequences to codon alignment, potentially clipping missing seqs
	#
	outfile		<- file.path('/Users/Oliver/duke/tmp', '150825_Region1_codonaligned_p_missing_p_HXB2.fasta')
	cmd			<- paste('mafft --anysymbol',' --keeplength --mapout',' --add "',file.seq.m,'" --auto "',file.path('/Users/Oliver/duke/tmp',gsub(' ','_',basename(susa.f))),'" > ', outfile, '', sep='')
	system(cmd)
	#	inspecting mapout: HXB2 no internal sites removed, only end clipped 
	#
	#	move to '~/Dropbox (SPH Imperial College)/PANGEA_alignments/Regional Alignments'
	#
	file.rename(outfile, file.path('~/Dropbox (SPH Imperial College)/PANGEA_alignments/Regional Alignments',basename(outfile)))
}
######################################################################################
#
######################################################################################
project.Rakai.aliRegion1.merge.with.PANGEA.160825<- function()
{	
	#	load Susanna s codon alignment and remove any PANGEA seqs
	susa.f	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151129.fasta'
	susa.f	<- '~/Dropbox (SPH Imperial College)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_missing_p_HXB2.fasta'
	susa.s	<- read.dna(susa.f, format='fasta')
	susa.d	<- data.table(TAXA=rownames(susa.s), DATA= factor(grepl('^PG', rownames(susa.s)), levels=c(TRUE,FALSE), label=c('PNG','LNL')))	
	in.s	<- susa.s[subset(susa.d, DATA!='PNG')[, TAXA],]	
	#	load PANGEA alignment	
	png.f	<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/PANGEA_orig/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
	load(png.f)	#sq, sqi, si
	sq		<- sq[grepl('^PG[0-9]+|HXB2',rownames(sq)),]
	#	required: HXB2 in alignment
	sqn		<- seq.align.based.on.common.reference(in.s, sq, return.common.sites=TRUE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	outfile	<- '~/Dropbox (SPH Imperial College)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_HXB2.fasta'
	write.dna(sqn, file=outfile, format='fasta', colsep='', nbcol=-1)
	
	#
	#	subset to Rakai and relabel
	#
	
	#	load summaryData
	infile.relabel	<- "~/Dropbox (SPH Imperial College)/PANGEA_alignments/Regional phylogenetic analyses/Region 1 gag analysis/SummaryofGagSequenceData.rda"				
	load(infile.relabel)
	s		<- as.data.table(summaryData)
	s[, TAXA:= gsub('a\\||b\\|','',seqid)]		#taxa names have been tinkered with. does not merge.
	s[, ID:=gsub('\\*.*','',TAXA)]	
	stopifnot( nrow(subset(s, grepl('\\|',TAXA)))==0 )
	stopifnot( nrow(subset(s, is.na(ID)))==0 )	
	#	get info on new alignment
	sqni		<- data.table(SEQID=rownames(sqn), IDX=seq_len(nrow(sqn)))		
	sqni[, SRC:= factor(grepl('^PG',SEQID),levels=c(TRUE,FALSE),labels=c('PNG','LANL'))]
	set(sqni, NULL, "TAXA_ID", sqni[, gsub('RCCS\\*|Ref\\*', '', SEQID)])
	
	sqni[, ID:= gsub('\\*.*','',TAXA_ID)]	
	sqni		<- merge(sqni, s, by='ID', all.x=1)
	#
	#	collect info on taxa missing in summaryData
	#
	stopifnot( nrow(subset(sqni, is.na(IDX)))==0 )
	not.in.summaryData	<- subset(sqni, is.na(relabel) & SRC!='PNG' & !grepl('HXB2', SEQID))
	#	deselect PANGEA ZA and BW
	sqni		<- subset(sqni, !(is.na(relabel) & SRC=='PNG' & !grepl('UG', SEQID)))
	#	
	set(sqni, sqni[, which(is.na(relabel) & SRC=='PNG' & grepl('UG', SEQID))], 'dataSource', 'UG-MRC')
	set(sqni, sqni[, which(grepl('HXB2', SEQID))], 'dataSource', 'reference')
	set(sqni, sqni[, which(is.na(relabel) & SRC!='PNG' & grepl('RCCS', SEQID))], 'dataSource', 'RCCS historical')
	set(sqni, sqni[, which(is.na(relabel) & SRC!='PNG' & !grepl('RCCS', SEQID))], 'dataSource', 'reference')
	#
	setkey(sqni, IDX)
	set(sqni, NULL, c('IDX','TAXA_ID','ID','SRC','TAXA'), NULL)
	set(sqni, NULL, 'index', 1:nrow(sqni))
	sq			<- sqn[sqni[, SEQID],]
	outfile	<- '~/Dropbox (SPH Imperial College)/PANGEA_alignments/Regional Alignments/150825_Region1UG_codonaligned_p_PANGEA151113_p_HXB2.fasta'
	write.dna(sqn, file=outfile, format='fasta', colsep='', nbcol=-1)
	summaryData	<- copy(sqni)
	save(sq, summaryData, file=gsub('\\.fasta','.rda',outfile))
	
	
	#
	#	old code perhaps useful at some point
	#
	
	tmp			<- sqni[, which(SRC=='LANL' & IDX<300)]
	set(sqni, tmp, 'SRC', 'COMPENDIUM')	
	tmp			<- as.character(sqn)
	tmp			<- data.table(	TAXA=rownames(tmp), 
								FIRST= apply( tmp, 1, function(x) which(!x%in%c('-','?'))[1] ),
								LAST= ncol(tmp)-apply( tmp, 1, function(x) which(!rev(x)%in%c('-','?'))[1] ) + 1L,						
								COV=apply(tmp, 1, function(x) length(which(!x%in%c('-','?')))	))
	sqni		<- merge(sqni,tmp,by='TAXA')
	tmp			<- subset(sqni, SRC=='LANL')[, max(LAST)]
	tmp			<- as.character(sqn[,seq_len(tmp)])
	tmp			<- data.table(	TAXA=rownames(tmp), 
								COV_REGION=apply(tmp, 1, function(x) length(which(!x%in%c('-','?')))	))
	sqni		<- merge(sqni,tmp,by='TAXA')
	sqni[, SUBT:=NA_character_]
	tmp			<- sqni[, which(SRC=='COMPENDIUM')]
	set(sqni, tmp, 'SUBT', sqni[tmp, sapply(strsplit(TAXA,'.',fixed=1),'[[',1)])
	tmp			<- sqni[, which(SRC=='LANL')]
	set(sqni, tmp, 'SUBT', sqni[tmp, sapply(strsplit(TAXA,'*',fixed=1),'[[',3)])
	tmp			<- sqni[, which(SUBT%in%c('-','U'))]
	set(sqni, tmp, 'SUBT', NA_character_)
	sqni[, ACCN:=NA_character_]
	tmp			<- sqni[, which(SRC=='COMPENDIUM')]
	set(sqni, tmp, 'ACCN', sqni[tmp, sapply(strsplit(TAXA,'.',fixed=1),function(x) rev(x)[1])])
	tmp			<- sqni[, which(SRC=='LANL')]
	set(sqni, tmp, 'ACCN', sqni[tmp, sapply(strsplit(TAXA,'*',fixed=1),'[[',1)])
	write.dna( sqn, format='fasta', file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201.fasta', colsep='', nbcol=-1)
	save(sqni, sqn, file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201.rda')
	
	sqn				<- read.dna(file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201_rm99gps.fasta',format='fasta')
	rownames(sqn)	<- gsub(' (stripped)','',rownames(sqn),fixed=1)
	#	get info
	sqni		<- data.table(TAXA=rownames(sqn), IDX=seq_len(nrow(sqn)))		
	sqni[, SRC:= factor(grepl('^PG',TAXA),levels=c(TRUE,FALSE),labels=c('PNG','LANL'))]
	tmp			<- sqni[, which(SRC=='LANL' & IDX<300)]
	set(sqni, tmp, 'SRC', 'COMPENDIUM')	
	tmp			<- as.character(sqn)
	tmp			<- data.table(	TAXA=rownames(tmp), 
			FIRST= apply( tmp, 1, function(x) which(!x%in%c('-','?'))[1] ),
			LAST= ncol(tmp)-apply( tmp, 1, function(x) which(!rev(x)%in%c('-','?'))[1] ) + 1L,						
			COV=apply(tmp, 1, function(x) length(which(!x%in%c('-','?')))	))
	sqni		<- merge(sqni,tmp,by='TAXA')
	tmp			<- subset(sqni, SRC=='LANL')[, max(LAST)]
	tmp			<- as.character(sqn[,seq_len(tmp)])
	tmp			<- data.table(	TAXA=rownames(tmp), 
			COV_REGION=apply(tmp, 1, function(x) length(which(!x%in%c('-','?')))	))
	sqni		<- merge(sqni,tmp,by='TAXA')
	sqni[, SUBT:=NA_character_]
	tmp			<- sqni[, which(SRC=='COMPENDIUM')]
	set(sqni, tmp, 'SUBT', sqni[tmp, sapply(strsplit(TAXA,'.',fixed=1),'[[',1)])
	tmp			<- sqni[, which(SRC=='LANL')]
	set(sqni, tmp, 'SUBT', sqni[tmp, sapply(strsplit(TAXA,'*',fixed=1),'[[',3)])
	tmp			<- sqni[, which(SUBT%in%c('-','U'))]
	set(sqni, tmp, 'SUBT', NA_character_)
	sqni[, ACCN:=NA_character_]
	tmp			<- sqni[, which(SRC=='COMPENDIUM')]
	set(sqni, tmp, 'ACCN', sqni[tmp, sapply(strsplit(TAXA,'.',fixed=1),function(x) rev(x)[1])])
	tmp			<- sqni[, which(SRC=='LANL')]
	set(sqni, tmp, 'ACCN', sqni[tmp, sapply(strsplit(TAXA,'*',fixed=1),'[[',1)])
	save(sqni, sqn, file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201_rm99gps.rda')	
	write.dna( sqn, format='fasta', file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201_rm99gps.fasta', colsep='', nbcol=-1)
}
######################################################################################
project.Rakai.checkforARVs.150910<- function()
{
	require(data.table)
	#	Kate, I generally use data.table rather than data.frame since it s faster and more versatile. 
	#	warning: it s a bit of a learning curve though to use it
	
	require(big.phylo)
	#	library(devtools)
	#	install_github("olli0601/big.phylo")
	
	f.arv	<- '~/Dropbox (SPH Imperial College)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_150909.rda'
	f.rccsid<- '~/Dropbox (SPH Imperial College)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/150625_PangeaBoxes.csv'
	f.sid	<- '~/Dropbox (SPH Imperial College)/PANGEA_data/SangerUpdates/2015-07-20_PANGEA_3.csv'
	f.seq	<- '~/Dropbox (SPH Imperial College)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908.fasta'
	#
	#	get IDs into OK format
	#
	load(f.arv)	#	loads arvdat
	arvdat	<- as.data.table(arvdat)
	d.rccsid<- as.data.table(read.csv(f.rccsid, stringsAsFactors=0))
	d.rccsid<- subset(d.rccsid, select=c(Pangea.ID, Rakai.Study.ID))
	setnames(d.rccsid, c('Rakai.Study.ID','Pangea.ID'), c('RCCS_ID','PNG_ID'))
	set(d.rccsid, NULL, 'RCCS_ID', d.rccsid[, substr(RCCS_ID, 1, nchar(RCCS_ID)-3)])
	set(d.rccsid, NULL, 'PNG_ID', d.rccsid[, gsub('-Sy*','',PNG_ID)])
	setnames(arvdat, 'RCCS_studyid', 'RCCS_ID')
	d.sid	<- as.data.table(read.csv(f.sid, stringsAsFactors=0, header=0))
	setnames(d.sid, c('V1','V2','V3'), c('PNG_ID_FULL', 'SNG_ID', 'x'))
	d.sid[, PNG_ID:= gsub('-S[0-9]+','',PNG_ID_FULL)]
	set(d.sid, NULL, 'SNG_ID', d.sid[, gsub('#','_',SNG_ID)])
	d.sid[, x:=NULL]
	#
	#	merge
	#	
	subset(arvdat, any.pangea==1)	
	#	126		
	arvdat	<- merge( arvdat, d.rccsid, by='RCCS_ID' )
	#	85
	arvdat	<- merge( arvdat, d.sid, by='PNG_ID')
	#	76
	
	#
	#	check for DRMs
	#
	seq				<- read.dna(file=f.seq, format='fasta')	
	seq				<- seq[c("B.FR.83.HXB2_LAI_IIIB_BRU.K03455",arvdat[, PNG_ID_FULL]),]
	rownames(seq)[1]	<- 'HXB2'
	outfile			<- '~/Dropbox (SPH Imperial College)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_PotentialDRMs_OR_150910.R'
	tmp				<- big.phylo:::seq.rm.drugresistance(seq, outfile=outfile)
	nodr.info		<- tmp$nodr.info
	nodr.seq		<- tmp$nodr.seq
}
######################################################################################
project.Rakai.checkforARVs.150911<- function()
{
	require(data.table)
	#	Kate, I generally use data.table rather than data.frame since it s faster and more versatile. 
	#	warning: it s a bit of a learning curve though to use it
	
	require(big.phylo)
	#	library(devtools)
	#	install_github("olli0601/big.phylo")
	
	f.arv			<- '~/Dropbox (SPH Imperial College)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_150911.rda'
	f.rccsid		<- '~/Dropbox (SPH Imperial College)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/150625_PangeaBoxes.csv'
	f.sid			<- '~/Dropbox (SPH Imperial College)/PANGEA_data/SangerUpdates/2015-07-20_PANGEA_3.csv'
	f.seq			<- '~/Dropbox (SPH Imperial College)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908.fasta'
	#
	load(f.arv)	#	loads rakaidat
	rakaidat		<- as.data.table(rakaidat)
	setnames(rakaidat, 'ProjectID', 'PNG_ID')
	#	482
	seq				<- read.dna(file=f.seq, format='fasta')	
	rownames(seq)	<- gsub('-S[0-9]+','',rownames(seq))
	d.seq			<- merge( data.table(PNG_ID=rownames(seq)), rakaidat, by='PNG_ID' )
	#	459 with PANGEA sequence
	
	
	#
	#	check for DRMs
	#
	seq				<- seq[c("B.FR.83.HXB2_LAI_IIIB_BRU.K03455",d.seq[, PNG_ID]),]		
	rownames(seq)[1]	<- 'HXB2'
	outfile			<- '~/Dropbox (SPH Imperial College)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_PotentialDRMs_OR_150911.R'
	tmp				<- big.phylo:::seq.rm.drugresistance(seq, outfile=outfile)
	nodr.info		<- tmp$nodr.info
	nodr.seq		<- tmp$nodr.seq
}