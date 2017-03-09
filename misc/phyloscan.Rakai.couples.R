RakaiCouples.setup.phyloscan.runs<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ape)
	
	wdir	<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples"
	#
	#	load table of sanger ids and pangea ids
	load("~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda")
	#	we need from load: dm,  sqi
	#
	#	load couples data
	infile	<- '~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/Pangea_Couples.csv'
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
	infile	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_UG_gag_coinfinput_160219.rda'
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
	infile		<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_gd.rda'
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
	#	fixup inconsistent phylotype run assignments
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
	save(pty.runs, file=file.path('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples','Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda'))
	#	add second dummy individual
	ans			<- rbind(ans, as.data.table(expand.grid(TAXA='PG14-UG501310-S03314', FILE_ID='15777_1_5', PARTNER='Other', COUPID='Other', PTY_RUN=ans[, sort(unique(PTY_RUN))])))
	setkey(ans, PTY_RUN, TAXA)
	pty.runs	<- copy(ans)
	save(pty.runs, file=file.path('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples','Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns_2dummy.rda'))
}

RakaiCouples.setup.phyloscan.runs.info<- function()
{
	
}
######################################################################################
######################################################################################
project.Rakai.checkMissingRakai.150307<- function()
{
	png.f	<- '~/Dropbox (Infectious Disease)/PANGEA_data/2016-01-20_PANGEA_stats_by_sample.csv'	
	#infile	<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/CheckPangeaId.csv'
	#infile	<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/CheckPangeaId2.csv'
	infile	<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/SequenceDataNOTAvailable_Documented.csv'
	infile	<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/SequenceDataAvailable_NotDocumented.csv'
	
	#
	#
	#
	dpng	<- as.data.table(read.csv(png.f))
	setnames(dpng, 'ProjectID', 'PANGEA_ID')
	df		<- as.data.table(read.csv(infile))
	setnames(df, 'x', 'PANGEA_ID')
	df		<- merge(dpng, df, by='PANGEA_ID')	
	write.csv(df, row.names=FALSE, file=gsub('.csv','_info.csv',infile))
	
	png.f	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/PANGEA_orig/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
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
#
#	add first sequences
#	
RakaiCirc.circ.timelines.addfirstseq.160816<- function(rt, rs)
{
	tmp		<- rs[, {
						z	<- which.min(DATE)
						ps	<- NA_character_						
						if(any(SEQTYPE=='Sanger not started by Jul2016'))
							ps	<- 'Sanger not started by Jul2016'												
						if(any(SEQTYPE=='Sanger failed'))
							ps	<- 'Sanger failed'												
						if(any(SEQTYPE=='Sanger completed without IVA'))
							ps	<- 'Sanger completed without IVA'												
						if(any(SEQTYPE=='Sanger completed with IVA'))
							ps	<- 'Sanger completed with IVA'						
						if(any(SEQTYPE=='PANGEA'))
						{
							ps	<- 'PANGEA'
							z	<- which.min(DATE)
						}													
						list(	VISIT=VISIT[z], 
								DATE=DATE[z], 
								SEQ_HIST= any(SEQTYPE=='HISTORIC'),
								SEQ_HIST_GAG= ifelse(any(SEQTYPE=='HISTORIC'), any(GENE_REGION[SEQTYPE=='HISTORIC']=='Region1'), FALSE),
								SEQ_PANGEA_NOW= any(SEQTYPE=='PANGEA'),
								SEQ_PANGEA_NOW_GAG= ifelse(any(SEQTYPE=='PANGEA'), any(GENE_REGION[SEQTYPE=='PANGEA']=='Region1'), FALSE),
								SEQ_N=length(DATE),
								PANGEA_SEQTYPE=ps,								 
								SEQ_NOW_GAG= any(GENE_REGION=='Region1'))
					}, by=RID]
	#	complete to all available RIDs
	tmp		<- merge(unique(subset(rt, select=RID)), tmp, all.x=1, by='RID')			
	#	define SEQ_TYPE of interest
	tmp[, SEQ_TYPE:= 'no_sequence']
	tmp2	<- tmp[, which(grepl('Sanger',PANGEA_SEQTYPE))]
	set(tmp, tmp2, 'SEQ_TYPE', gsub(' ','_',tmp[tmp2, PANGEA_SEQTYPE]))	
	set(tmp, tmp[, which(!SEQ_HIST_GAG & SEQ_HIST & !SEQ_PANGEA_NOW)], 'SEQ_TYPE', "other_gene_historic_only")
	set(tmp, tmp[, which(SEQ_HIST_GAG & SEQ_HIST & !SEQ_PANGEA_NOW)], 'SEQ_TYPE', "gag_or_partial_gag_historic_only")
	set(tmp, tmp[, which(!SEQ_PANGEA_NOW_GAG & SEQ_PANGEA_NOW)], 'SEQ_TYPE', "other_gene_PANGEA")
	set(tmp, tmp[, which(SEQ_PANGEA_NOW_GAG & SEQ_PANGEA_NOW)], 'SEQ_TYPE', "gag_or_partial_gag_PANGEA")	
	#	
	tmp[, SEQ_NOW:= factor(SEQ_TYPE%in%c('gag_or_partial_gag_historic_only','gag_or_partial_gag_PANGEA','other_gene_historic_only'), levels=c(TRUE,FALSE), labels=c('Y','N'))]
	set(tmp, NULL, 'SEQ_NOW_GAG', tmp[, factor(SEQ_NOW=='Y' & !is.na(SEQ_NOW_GAG) & SEQ_NOW_GAG, levels=c(TRUE,FALSE), labels=c('Y','N'))])
	set(tmp, NULL, c('VISIT','SEQ_HIST','SEQ_PANGEA_NOW','SEQ_N', 'PANGEA_SEQTYPE','SEQ_HIST_GAG','SEQ_PANGEA_NOW_GAG'), NULL)
	setnames(tmp, 'DATE', 'FIRST_SEQ_DATE')
	rt		<- merge(rt, tmp, by='RID')
	rt
}
#
#	make individual timeline over visits
#
RakaiCirc.circ.timelines.init.160816<- function(rh, rd)
{	
	rt		<- rbind( subset(rh, select=c(RID, VISIT)), subset(rd, select=c(RID, VISIT)) )
	setkey(rt, RID, VISIT)
	rt		<- unique(rt)		
	cat('\nall visit data', nrow(rt), '\tall data in history',nrow(rh))
	rt[, ROUND:=VISIT]
	set(rt, rt[, which(ROUND==15.1)], 'ROUND', 15)
	#	expand to all possible rounds after entry into cohort
	tmp		<- as.data.table(expand.grid(RID=rt[, unique(RID)], ROUND=rt[, seq.int(1,max(ROUND))]))
	rt		<- merge(rt, tmp, by=c('RID','ROUND'),all=1)	
	#	add basic demographic data
	tmp		<- subset(rd, select=c(RID, SEX, RELIGION, CAUSE_OF_DEATH, BIRTHDATE, LASTNEGDATE, FIRSTPOSDATE, EST_DATEDIED))
	setkey(tmp, RID)
	tmp		<- unique(tmp)
	set(tmp, tmp[, which(is.na(RELIGION))], 'RELIGION', 'missing')
	set(tmp, NULL, 'SEX', tmp[, factor(SEX, levels=c('M','F'), labels=c('male','female'))])
	rt		<- merge(tmp, rt, by='RID')
	#	add start and end of visit rounds	
	tmp		<- rd[, list(ROUND_ST= round(min(DATE),d=1), ROUND_E= round(max(DATE),d=1) ), by='VISIT']
	setnames(tmp, 'VISIT','ROUND')
	#	TODO: rounds are overlapping? Here, I cut the end times of each period so that the periods do not overlap (hack).
	setkey(tmp, ROUND)
	set(tmp, NULL, 'ROUND_E2', c(tmp[, ROUND_ST-.1][-1], Inf) )
	tmp2	<- tmp[, which(ROUND_E2<ROUND_E)]
	set(tmp, tmp2, 'ROUND_E', tmp[tmp2,ROUND_E2])
	set(tmp, NULL, 'ROUND_E2', NULL)
	#	TODO: confirm no data from round 5?
	tmp		<- rbind(tmp, data.table(ROUND=5, ROUND_ST=1998.4, ROUND_E=1999.2))	
	rt		<- merge(rt, tmp, by='ROUND')
	setkey(rt, RID, ROUND)
	#	delete rounds before diagnosis or known visit
	rt		<- subset(rt, !is.na(VISIT) | (is.na(VISIT) & FIRSTPOSDATE<ROUND_E))
	stopifnot( nrow(subset(rt[, all(is.na(VISIT)), by=RID], V1))==0 )
	#	delete rounds after recorded death event
	rt		<- subset(rt, is.na(EST_DATEDIED) | (!is.na(EST_DATEDIED) & ROUND_ST<EST_DATEDIED))
	stopifnot( nrow(subset(rt[, all(is.na(VISIT)), by=RID], V1))==0 )
	#	define missing rounds as 
	#		'no_follow_up'			temporarily not in follow up, visits before and after the current round
	#		'before_first_visit'	already estimated to be HIV+, first visit after the current round
	#		'no_obs_est_died'		after last visit, estimated to have died in the current round
	#		'not_seen_again'		after last visit and no record of death
	tmp		<- rt[, list(ROUND_FIRST_ST=min(ROUND_ST[!is.na(VISIT)]), ROUND_LAST_E=max(ROUND_E[!is.na(VISIT)])), by='RID']	
	rt		<- merge(rt, tmp, by='RID')
	set(rt, NULL, 'VISIT', rt[, as.character(VISIT)])
	set(rt, rt[, which(is.na(VISIT) & ROUND_ST>ROUND_FIRST_ST & ROUND_E<ROUND_LAST_E)], 'VISIT', 'no_follow_up')	
	set(rt, rt[, which(is.na(VISIT) & ROUND_ST<ROUND_FIRST_ST)], 'VISIT', 'before_first_visit')
	set(rt, rt[, which(is.na(VISIT) & ROUND_E==ROUND_LAST_E)], 'VISIT', 'no_obs_est_died')
	set(rt, rt[, which(is.na(VISIT) & ROUND_E>ROUND_LAST_E)], 'VISIT', 'not_seen_again')
	stopifnot( nrow(subset(rt, is.na(VISIT)))==0 )
	#	define first circumcision 
	tmp		<- subset(rh, CIRC=='Y')[, list(FIRSTCIRC_V= min(VISIT) ), by='RID']
	tmp2	<- subset(rh, CIRC=='N')[, list(LASTNOTCIRC_V= max(VISIT) ), by='RID']
	tmp		<- merge(tmp, tmp2, by='RID', all=1)	
	#	TODO individuals with conflicting circumcision data. For now ignore conflicts (hack).
	tmp2	<- subset(tmp, FIRSTCIRC_V<=LASTNOTCIRC_V)
	write.csv(tmp2, row.names=FALSE, file=file.path(wdir, 'check_conflictingCircumcisionStatus.csv'))
	cat('\nWarning: Found conflicting circumcision data for n=',nrow(tmp2))
	set(tmp, tmp[, which(FIRSTCIRC_V<=LASTNOTCIRC_V)], 'LASTNOTCIRC_V', NA_real_)
	#	define ever circumcised and all missing
	tmp2	<- rh[, list(EVERCIRC=as.integer(any(CIRC=='Y', na.rm=1)), NODATA=all(is.na(CIRC))), by='RID']
	tmp		<- merge(tmp, tmp2, by='RID', all=1)	
	set(tmp, tmp[, which(NODATA)], 'EVERCIRC', 2L)
	#	tmp[, table(EVERCIRC, NODATA, useNA='if')]
	set(tmp, NULL, 'EVERCIRC', tmp[, factor(EVERCIRC, levels=c(0,1,2), labels=c('N','Y','missing_all'))])
	set(tmp, NULL, 'NODATA', NULL)
	rt		<- merge(rt, tmp, by='RID')
	#	define circumcision timeline
	rt[, CIRC:=NA_character_]	
	set(rt, rt[, which(FIRSTCIRC_V<=ROUND)], 'CIRC', 'Y')
	set(rt, rt[, which(ROUND<=LASTNOTCIRC_V)], 'CIRC', 'N')
	set(rt, rt[, which(EVERCIRC=='missing_all')], 'CIRC', 'missing_all')
	set(rt, rt[, which(EVERCIRC=='N')], 'CIRC', 'N')
	set(rt, rt[, which(is.na(CIRC) & FIRSTCIRC_V>ROUND)], 'CIRC', 'missing_before_circ')
	set(rt, rt[, which(is.na(CIRC) & LASTNOTCIRC_V<ROUND)], 'CIRC', 'missing_after_notcirc')
	stopifnot( nrow(subset(rt, is.na(CIRC)))==0 )
	set(rt, NULL, c('FIRSTCIRC_V','LASTNOTCIRC_V'), NULL)
	rt
}

RakaiCirc.recipient.female.get.info<- function(wdir=NA)
{
	infile				<- "~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda"
	indir.historicseq	<- "~/Dropbox (Infectious Disease)/PANGEA_alignments/Rakai Data for IqTree"
	wdir				<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
	#
	#	prepare RCCS data
	#	
	load(infile)
	#	a bit of clean up 
	rd		<- as.data.table(rccsData)
	setnames(rd, colnames(rd), gsub('\\.','_',toupper(colnames(rd))))	
	for(x in colnames(rd))
		if(class(rd[[x]])=='Date')
			set(rd, NULL, x, hivc.db.Date2numeric(rd[[x]]))
	rh		<- as.data.table(rccsHistory)
	setnames(rh, colnames(rh), gsub('\\.','_',toupper(colnames(rh))))
	#rd[, table(VISIT)]
	#	make shorter
	setnames(rd, 'RCCS_STUDYID', 'RID')
	setnames(rd, 'PANGEA_ID', 'PID')
	setnames(rd, 'STUDYID', 'SID')
	setnames(rh, 'RCCS_STUDYID', 'RID')	
	
	rrec	<- subset(rd, SEX=='F') 
	stopifnot( nrow(subset(rrec, is.na(LASTNEGVIS) & !is.na(LASTNEGDATE)))==0 )
	stopifnot( nrow(subset(rrec, !is.na(LASTNEGVIS) & is.na(LASTNEGDATE)))==0 )
	rrec	<- subset(rrec, !is.na(LASTNEGVIS), select=c(RID, SEX, VISIT, DATE, REGION, COMM_NUM, HH_NUM, BIRTHDATE, LASTNEGVIS, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE))
	stopifnot( nrow(subset(rrec, FIRSTPOSDATE<LASTNEGDATE))==0 )
	rrec[, SC_WINDOW:= FIRSTPOSDATE-LASTNEGDATE]		
	#	add first sequences. load sequence info. expect "rs".
	#	(if not exist, run RakaiCirc.seq.get.info() )
	load(file.path(wdir,'RCCS_SeqInfo_160816.rda'))	
	rs		<- subset(rs, !is.na(VISIT))	
	rrec	<- RakaiCirc.circ.timelines.addfirstseq.160816(rrec, rs)	
	#	rm B009659 (check) has missing first pos date, was sequenced, but Sanger_completed_without_IVA
	rrec	<- subset(rrec, RID!='B009659')
	#	rm B106184 (check) sequence date before last neg date 
	rrec	<- subset(rrec, RID!='B106184')
	stopifnot( nrow(subset(rrec, is.na(FIRSTPOSDATE)))==0 )
	stopifnot( nrow(subset(rrec, FIRSTPOSDATE<=LASTNEGDATE))==0 )
	#	make unique by RID
	set(rrec, NULL, c('VISIT','DATE'), NULL)
	rrec	<- unique(rrec)
	if(!is.na(wdir))
	{
		#	plot recipient sequence coverage
		tmp		<- subset(rrec, SC_WINDOW<3.4)
		tmp[, SC_WINDOWc:= cut(SC_WINDOW, breaks=seq(0,3.4,.2))]		
		set(tmp, NULL, 'SEQ_TYPE', tmp[, factor(SEQ_TYPE, levels=c('gag_or_partial_gag_PANGEA','gag_or_partial_gag_historic_only','Sanger_completed_with_IVA','Sanger_completed_without_IVA','Sanger_not_started_by_Jul2016','other_gene_historic_only','Sanger_failed','no_sequence'))])                           
		tmp		<- tmp[, list(N=length(RID)), by=c('SC_WINDOWc','SEQ_TYPE')]
		tmp		<- merge( tmp, data.table(expand.grid(SC_WINDOWc= tmp[, unique(SC_WINDOWc)], SEQ_TYPE= tmp[, unique(SEQ_TYPE)])), by=c('SC_WINDOWc','SEQ_TYPE'), all=1)
		set(tmp, tmp[, which(is.na(N))], 'N', 0)	
		setkey(tmp, SEQ_TYPE, SC_WINDOWc)
		tmp		<- tmp[, list(CN= cumsum(N), SC_WINDOWc=SC_WINDOWc), by='SEQ_TYPE']	
		ggplot(tmp, aes(x=SC_WINDOWc, y=CN, fill=SEQ_TYPE)) + geom_bar(position='stack', stat='identity') + 			
				scale_y_continuous(breaks=seq(0,5e3,2e2)) +
				scale_fill_brewer(palette='Set3') + 
				theme_bw() + theme(legend.position='bottom') + guides(fill=guide_legend(ncol=2)) +
				labs(x='\nseroconversion window\n(years)', y='female RCCS participants\nwith last negative test\n(cumulated counts)\n')
		ggsave(file.path(wdir, 'RCCSfemalerecipient_by_sequencecoverage.pdf'), w=14, h=8)
	}	
	rrec
}

RakaiCirc.epi.get.info<- function()
{
	wdir				<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"	
	infile				<- "~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda"
	load(infile)
	#	a bit of clean up 
	rd		<- as.data.table(rccsData)
	setnames(rd, colnames(rd), gsub('\\.','_',toupper(colnames(rd))))	
	for(x in colnames(rd))
		if(class(rd[[x]])=='Date')
			set(rd, NULL, x, hivc.db.Date2numeric(rd[[x]]))
	rh		<- as.data.table(rccsHistory)
	setnames(rh, colnames(rh), gsub('\\.','_',toupper(colnames(rh))))
	#rd[, table(VISIT)]
	#	make shorter
	setnames(rd, 'RCCS_STUDYID', 'RID')
	setnames(rd, 'PANGEA_ID', 'PID')
	setnames(rd, 'STUDYID', 'SID')
	setnames(rh, 'RCCS_STUDYID', 'RID')	
	#	data checks
	setkey(rh, VISIT, RID)
	stopifnot(nrow(rh)==nrow(unique(rh)))
	#	define circumcision	
	set(rh, rh[, which(!CIRC%in%c(1,2))], 'CIRC', NA_integer_)
	set(rh, NULL, 'CIRC', rh[, factor(CIRC, levels=c(1,2), labels=c('Y','N'))])
	#	TODO Warning: found 1 female circumcised A006734
	tmp		<- rh[, which(CIRC=='Y' & SEX=='F')]
	cat('\nWarning: found female circumcised',rh[tmp, paste(RID, collapse=' ')])
	set(rh, tmp, 'CIRC', NA_integer_)
	list(rd=rd, rh=rh)
}

RakaiCirc.seq.get.info<- function()
{
	infile				<- "~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData.rda"	
	infile.sangerstats	<- "~/Dropbox (Infectious Disease)/PANGEA_data/2016-07-07_PANGEA_stats_by_sample.csv"
	infile.relabel		<- "~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional phylogenetic analyses/Region 1 gag analysis/SummaryofGagSequenceData.rda"
	infile.assembly		<- "~/Dropbox (Infectious Disease)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908_Summary.csv"
	wdir				<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
	#
	#	read all processed RCCS sequences 
	#	
	infile.region1		<- '/Users/Oliver/Dropbox (Infectious Disease)/PANGEA_alignments/Rakai Data for IqTree/Prelim_RakaiPangeaSqnAndMetaData_IqTree_Region1.rda'
	infile.region2		<- '/Users/Oliver/Dropbox (Infectious Disease)/PANGEA_alignments/Rakai Data for IqTree/Prelim_RakaiPangeaSqnAndMetaData_IqTree_Region2.rda'
	infile.region3		<- '/Users/Oliver/Dropbox (Infectious Disease)/PANGEA_alignments/Rakai Data for IqTree/Prelim_RakaiPangeaSqnAndMetaData_IqTree_Region3.rda'
	infile.region4		<- '/Users/Oliver/Dropbox (Infectious Disease)/PANGEA_alignments/Rakai Data for IqTree/Prelim_RakaiPangeaSqnAndMetaData_IqTree_Region4.rda'
	
	load(infile.region1)
	tmp		<- as.character(gag.sqn)
	tmp		<- data.table(	relabel= rownames(tmp),
							BASE_N= apply(tmp, 1, function(x) sum(as.numeric(x%in%c('a','t','g','c')))) )
	tmp		<- merge(tmp, as.data.table(subset(summaryData, select=c('seqid', 'relabel'))), by='relabel', all.x=1)
	#	get Pangea data from region 1
	z2		<- subset(as.data.table(pangeaMetaData), select=c(RCCS_studyid, visit, date, seqid, CometSubtype, Pangea.id))	
	z2		<- merge(z2, tmp, by='seqid', all.x=1)
	stopifnot( nrow(subset(z2, is.na(Pangea.id)))==0 )
	z2		<- subset(z2, BASE_N>0)
	z2[, SEQTYPE:='PANGEA']
	z2[, GENE_REGION:= 'Region1']
	rs		<- copy(z2)
	#	add historical data from region 1
	z2		<- subset(as.data.table(rakhistMetaData), select=c(RCCS_studyid, visit, date, seqid, CometSubtype))	
	z2		<- merge(z2, tmp, by='seqid', all.x=1)
	stopifnot( nrow(subset(z2, is.na(RCCS_studyid)))==0 )
	z2		<- subset(z2, BASE_N>0)
	z2[, SEQTYPE:='HISTORIC']
	z2[, GENE_REGION:= 'Region1']
	rs		<- rbind(rs, z2, use.names=TRUE, fill=TRUE)
	#	
	load(infile.region2)
	tmp		<- as.character(pol.sqn)
	tmp		<- data.table(	relabel= rownames(tmp),
			BASE_N= apply(tmp, 1, function(x) sum(as.numeric(x%in%c('a','t','g','c')))) )
	tmp		<- merge(tmp, as.data.table(subset(summaryData, select=c('seqid', 'relabel'))), by='relabel', all.x=1)
	#	get Pangea data from region 2
	z2		<- subset(as.data.table(pangeaMetaData), select=c(RCCS_studyid, visit, date, seqid, CometSubtype, Pangea.id))	
	z2		<- merge(z2, tmp, by='seqid', all.x=1)
	stopifnot( nrow(subset(z2, is.na(Pangea.id)))==0 )
	z2		<- subset(z2, BASE_N>0)
	z2[, SEQTYPE:='PANGEA']
	z2[, GENE_REGION:= 'Region2']
	rs		<- rbind(rs, z2, use.names=TRUE, fill=TRUE)
	#	no historical data from region 2
	#	
	load(infile.region3)
	#	TODO: there are 729 duplicates in env.sqn ie RCCS*A108646*vis15*reg13*com774*hh75*F*17*prev*2012.37
	tmp		<- as.character(env.sqn)
	tmp		<- data.table(	relabel= rownames(tmp),
							BASE_N= apply(tmp, 1, function(x) sum(as.numeric(x%in%c('a','t','g','c')))) )
	tmp		<- subset(tmp, BASE_N>0)	
	setkey(tmp, relabel)	
	tmp		<- merge(unique(tmp), as.data.table(subset(summaryData, select=c('seqid', 'relabel'))), by='relabel', all.x=1)
	#	get Pangea data from region 3
	z2		<- subset(as.data.table(pangeaMetaData), select=c(RCCS_studyid, visit, date, seqid, CometSubtype, Pangea.id))	
	z2		<- merge(z2, tmp, by='seqid', all.x=1)
	stopifnot( nrow(subset(z2, is.na(Pangea.id)))==0 )
	z2		<- subset(z2, BASE_N>0)
	z2[, SEQTYPE:='PANGEA']
	z2[, GENE_REGION:= 'Region3']
	rs		<- rbind(rs, z2, use.names=TRUE, fill=TRUE)
	#	no historical data from region 2
	#
	load(infile.region4)
	tmp		<- as.character(gp.sqn)
	tmp		<- data.table(	relabel= rownames(tmp),
							BASE_N= apply(tmp, 1, function(x) sum(as.numeric(x%in%c('a','t','g','c')))) )
	tmp		<- merge(tmp, as.data.table(subset(summaryData, select=c('seqid', 'relabel'))), by='relabel', all.x=1)
	#
	#	TODO: no PANGEA data from region 4
	#
	#	add historical data from region 4
	z2		<- subset(as.data.table(rakhistMetaData), select=c(RCCS_studyid, visit, date, seqid, CometSubtype))	
	z2		<- merge(z2, tmp, by='seqid', all.x=1)
	stopifnot( nrow(subset(z2, is.na(RCCS_studyid)))==0 )
	z2		<- subset(z2, BASE_N>0)
	z2[, SEQTYPE:='HISTORIC']
	z2[, GENE_REGION:= 'Region4']
	rs		<- rbind(rs, z2, use.names=TRUE, fill=TRUE)
	#
	#	reading done 
	#	clean up names etc
	#
	setnames(rs, colnames(rs), gsub('\\.','_',toupper(colnames(rs))))
	setnames(rs, c('RCCS_STUDYID','PANGEA_ID'), c('RID','PID'))
	set(rs, NULL, 'DATE', hivc.db.Date2numeric(rs[['DATE']]))
	set(rs, NULL, 'RELABEL', NULL)
	rs		<- subset(rs, !is.na(RID))		
	#
	#	add unprocessed but shipped PANGEA sequences
	#	
	load(infile)
	#	a bit of clean up 
	rd		<- as.data.table(rccsData)
	setnames(rd, colnames(rd), gsub('\\.','_',toupper(colnames(rd))))	
	for(x in colnames(rd))
		if(class(rd[[x]])=='Date')
			set(rd, NULL, x, hivc.db.Date2numeric(rd[[x]]))
	#	make shorter
	setnames(rd, 'RCCS_STUDYID', 'RID')
	setnames(rd, 'PANGEA_ID', 'PID')	
	#	
	tmp		<- subset(rd, !is.na(PID), select=c(RID, VISIT, DATE, PID))
	tmp2	<- setdiff(tmp[, PID], subset(rs, !is.na(PID))[, unique(PID)])
	cat('\nWarning: Found PANGEA sequences not in rd (error)? n=', length(setdiff(subset(rs, !is.na(PID))[, unique(PID)], tmp[, PID] )))	
	cat('\nWarning: PANGEA sequences not yet processed n=', length(tmp2))	
	setkey(tmp, RID)	
	#	TODO out of curiosity why are these 349 duplicated?
	cat('\nFound participants with more than one PANGEA sequence. n=',tmp[duplicated(tmp),][, length(unique(RID))],'\nids=',tmp[duplicated(tmp),][, paste(unique(RID), collapse=' ')])
	#	TODO there are 2674 unprocessed PANGEA sequences
	tmp		<- merge(tmp, data.table(PID=tmp2), by='PID')	
	tmp[, SEQTYPE:='PANGEA_not_yet_processed']	
	rs		<- rbind(rs, tmp, fill=TRUE, use.names=TRUE)
	#
	#	add 160120 statistics from Sanger
	#
	tmp		<- as.data.table(read.csv(infile.sangerstats, stringsAsFactors=FALSE))	
	setnames(tmp, colnames(tmp), gsub('\\.','_',toupper(colnames(tmp))))
	setnames(tmp, c('PROJECTID','STATUS'), c('PID','SANGER_STATUS'))
	tmp		<- subset(tmp, select=c(PID, SANGER_STATUS, ASSEMBLED, HIVCONTIG))
	tmp		<- tmp[-nrow(tmp),]
	set(tmp, tmp[, which(HIVCONTIG=='HIVcontig')],  'SANGER_STATUS',  'Sanger completed with IVA')
	set(tmp, tmp[, which(SANGER_STATUS=='Assume no HIV')], 'SANGER_STATUS', 'Sanger completed without IVA')
	set(tmp, tmp[, which(SANGER_STATUS=='Assume assembly failed')], 'SANGER_STATUS', 'Sanger failed')
	set(tmp, tmp[, which(SANGER_STATUS=='Assume sequencing failed')], 'SANGER_STATUS', 'Sanger failed')
	set(tmp, tmp[, which(SANGER_STATUS=='No sample')], 'SANGER_STATUS', 'Sanger failed')
	set(tmp, tmp[, which(SANGER_STATUS=='ID not in system yet')], 'SANGER_STATUS', 'Sanger not started by Jul2016')
	rs		<- merge(rs, subset(tmp, select=c(PID, SANGER_STATUS)), by='PID', all.x=1)
	tmp		<- rs[, which(SEQTYPE=='PANGEA_not_yet_processed')]
	set(rs, tmp, 'SEQTYPE', rs[tmp, SANGER_STATUS])
	set(rs, rs[, which(is.na(SEQTYPE))], 'SEQTYPE', 'Sanger not started by Jul2016')
	set(rs, NULL, 'SANGER_STATUS', NULL)
	#
	#	add Sanger IDs
	#
	tmp		<- as.data.table(read.csv(infile.assembly, stringsAsFactors=FALSE))
	setnames(tmp, c('Sanger.ID','PANGEA.ID'), c('SID','PID'))
	set(tmp, NULL, 'PID', tmp[, gsub('-S[0-9]+','',PID)])
	set(tmp, NULL, 'PID', tmp[, gsub('^\\s','',PID)])
	rs		<- merge(rs, subset(tmp, select=c(PID, SID)), by='PID', all.x=1)
	#	NOTE: there are a few PID with multiple SIDs !
	#
	#	add Phylo ID for gag tree
	#
	load(infile.relabel)
	tmp		<- subset(as.data.table(summaryData), select=c(seqid,relabel))
	setnames(tmp, c('seqid','relabel'),c('SEQID','SEQIDb'))
	rs		<- merge(rs, tmp, by='SEQID',all.x=1)
	#	> rs[, table(SEQTYPE, GENE_REGION)]
	#                               GENE_REGION
	#	SEQTYPE                         Region1 Region2 Region3 Region4 <NA>
  	#	HISTORIC                         2381       0       0    2608    0
  	#	PANGEA                           2905    2401    1302       0    0
  	#	Sanger completed with IVA           0       0       0       0  598
  	#	Sanger completed without IVA        0       0       0       0  286
  	#	Sanger failed                       0       0       0       0  292
  	#	Sanger not started by Jul2016       0       0       0       0 1501
	#
	#	save to file
	#
	save(rs, file=file.path(wdir,'RCCS_SeqInfo_160816.rda'))
}

RakaiCirc.various<- function()
{
	require(ape)
	require(data.table)
	if(0)
	{
		wdir	<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
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
		
		wdir	<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
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
	wdir	<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
	#
	#	start with Kate s GP24 FastTree
	#
	infile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional phylogenetic analyses/Region 1 gag analysis/gagAllSmall.nwk'
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
	infile	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_UG_gag_fasttree.rda'
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
	infile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1UG_codonaligned_p_PANGEA151113_p_HXB2.rda'
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

RakaiCouples.process.couples.161007<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	setkey(rp, COUPID)
	unique(rp)[, table(COUP_SC)]
	#	total of couples with assigned SANGER_ID in RCCS
	#   F->M    M->F seroinc seropos 
	# 	10      19      52     266 	
	merge(unique(rp), unique(subset(pty.runs, COUPID!='Other', COUPID)), by='COUPID')[, table(COUP_SC)]
	#	total of couples with SANGER_ID for which I have data
	#	F->M    M->F seroinc seropos 
	#	7       16      45     235 
	#
	#	collect runs
	#
	infiles	<- data.table(	FILE= c(	'~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161007/RCCS_161007_w270_phscout.rda' ))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_phscout.rda','',basename(FILE))]				
	#
	#	for each run: get list of pairs
	#	
	rps		<- infiles[, {
				#FILE	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161007/RCCS_161007_w270_phscout.rda'
				load(FILE)	#loads phs dtrms dtrees
				#	select likely pairs -- these are ordered
				dtrms[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]				
				tmp		<- dcast.data.table(dtrms, PAIR_ID~TYPE, value.var='WIN_OF_TYPE_P')
				set(tmp, tmp[, which(is.na(disconnected))], 'disconnected', 0)
				set(tmp, tmp[, which(is.na(int))], 'int', 0)
				set(tmp, tmp[, which(is.na(unint))], 'unint', 0)
				set(tmp, tmp[, which(is.na(cher))], 'cher', 0)
				set(tmp, tmp[, which(is.na(trans_12))], 'trans_12', 0)
				set(tmp, tmp[, which(is.na(trans_21))], 'trans_21', 0)
				tmp		<- melt.data.table(tmp, id.vars='PAIR_ID', value.name='WIN_OF_TYPE_P', variable.name='TYPE')				
				stopifnot( tmp[, list(CHECK=sum(WIN_OF_TYPE_P)),by='PAIR_ID'][, all(abs(CHECK-1)<=2*.Machine$double.eps)] )				
				dtrms	<- merge(unique(subset(dtrms, select=c(PAIR_ID, ID1, ID2, PTY_RUN, WIN_TOTAL, SCORE))), tmp, by='PAIR_ID')
				#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
				tmp		<- copy(dtrms)
				setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
				set(tmp, tmp[, which(TYPE=='trans_12')], 'TYPE', 'trans')
				set(tmp, tmp[, which(TYPE=='trans_21')], 'TYPE', 'trans_12')
				set(tmp, tmp[, which(TYPE=='trans')], 'TYPE', 'trans_21')
				dtrms	<- rbind(tmp, dtrms)
				#	first individual is always male	
				setnames(dtrms, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
				set(dtrms, dtrms[, which(TYPE=='trans_12')], 'TYPE', 'trans_mf')
				set(dtrms, dtrms[, which(TYPE=='trans_21')], 'TYPE', 'trans_fm')
				set(dtrms, NULL, 'TYPE', dtrms[, as.character(TYPE)])
				dtrms	<- merge(dtrms, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))				
				dtrms
			}, by=c('RUN','DIR','FILE')]
	save(rps, file= '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161007/RCCS_161007_w270_assignments.rda')
	#
	cat('\nNumber of couples',rps[, length(unique(COUPID))])
	setkey(rps, RUN, COUPID, MALE_SANGER_ID, FEMALE_SANGER_ID)
	cat('\nNumber of sequence pairs from couples',nrow(unique(rps)))
	setkey(rps, RUN, COUPID, PAIR_ID)
	cat('\nNumber of pairings (including repeated sequence pairs)',nrow(unique(rps)))	
	#
	#	for each run: plot pairs	
	#	
	#for( run in rps[, unique(RUN)] )
	#{			
	run		<- 'RCCS_161007_w270'
	dir		<- subset(rps, RUN==run)[1,DIR]
	cat('\ndir is',dir,'\trun is',run)
	df		<- subset(rps, RUN==run)
	setkey(df, MALE_SANGER_ID, FEMALE_SANGER_ID, PAIR_ID)					
	#
	#	plot evidence
	#		
	tmp		<- unique(df)
	tmp[, PLOT_ID:= as.numeric(gsub('-','\\.',PAIR_ID))]
	tmp		<- tmp[order(-PLOT_ID),]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair', PAIR_ID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- merge(subset(tmp, select=c(PAIR_ID, LABEL)), df, by='PAIR_ID')
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('trans_mf','trans_fm','unint','int','cher','disconnected'), labels=c('M transmit to F','F transmit to M','M, F are unint','M, F are intermingled','M, F are a cherry','M, F are disconnected'))])
	tmp		<- subset(tmp, select=c(PAIR_ID, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, WIN_OF_TYPE_P, WIN_TOTAL, TYPE, LABEL))
	tmp[, WIN_OF_TYPE_N:=WIN_OF_TYPE_P*WIN_TOTAL]
	#tmp		<- melt(tmp, measure.vars=c('WIN_OF_TYPE_N', 'WIN_OF_TYPE_P'))
	#set(tmp, tmp[, which(variable=='WIN_OF_TYPE_N')],'variable','number of read windows')
	#set(tmp, tmp[, which(variable=='WIN_OF_TYPE_P')],'variable','proportion of read windows')	
	tmp2	<- subset(tmp, COUP_SC=='F->M')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC, ncol=2) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_F2M.pdf',sep='')), w=25, h=max(4,0.23*nrow(tmp2)), limitsize = FALSE)
	tmp2	<- subset(tmp, COUP_SC=='M->F')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_grid(~COUP_SC) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_M2F.pdf',sep='')), w=25, h=max(3,0.23*nrow(tmp2)), limitsize = FALSE)
	tmp2	<- subset(tmp, COUP_SC=='seroinc')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_grid(~COUP_SC) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_SeroInc.pdf',sep='')), w=25, h=max(3,0.23*nrow(tmp2)), limitsize = FALSE)
	tmp2	<- subset(tmp, COUP_SC=='seropos')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_grid(~COUP_SC) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_SeroPos.pdf',sep='')), w=25, h=max(3,0.23*nrow(tmp2)), limitsize = FALSE)
	#
	#	tabulate correct classifications with phyloscanner for serodiscordant couples
	#
	#	correct: trm between M and F.
	rpa		<- subset(df, COUP_SC=='F->M' | COUP_SC=='M->F')[, list(CLASS='ancestral in either direction\nor intermingled', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='trans_mf'|TYPE=='trans_fm'|TYPE=='int'|TYPE=='cher'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(df, COUP_SC=='F->M')[, list(CLASS='ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='trans_fm'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp)
	tmp		<- subset(df, COUP_SC=='M->F')[, list(CLASS='ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='trans_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp)
	#	rank each couple
	tmp		<- rpa[order(CLASS, -CLASS_PROP),][, list(PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)), by='CLASS']
	rpa		<- merge(rpa, tmp, by=c('CLASS','PAIR_ID'))	
	setkey(rpa, CLASS, CLASS_RANK)
	#	plot by rank
	ggplot(rpa, aes(x=CLASS_PROP, y=CLASS_RANK, colour=CLASS)) + 
			geom_point() + geom_step() +
			geom_text(aes(label=PAIR_ID), size=2, nudge_x=-.05, nudge_y=.5) +
			scale_y_continuous(breaks=seq(0,50,5)) +
			scale_x_reverse(labels = scales::percent, breaks=seq(0,1,0.1)) +
			scale_colour_brewer(palette='Set1') +
			facet_grid(~CLASS) +
			labs(	x= '\nminimum proportion\n(proportion of ancestral windows out of all windows\nthat have reads from both individuals is at least x%)', 
					y='sequence pairs\n(#)\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-correctancestral.pdf',sep='')), w=12, h=7)
	#	write to file
	tmp			<- subset(rpa, CLASS=='ancestral in either direction\nor intermingled', select=c(COUPID, PAIR_ID, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, CLASS, CLASS_PROP, CLASS_RANK))
	write.csv(tmp, row.names=FALSE, file=file.path(dir, paste(run,'-phsc-serodiscpairs-assignments.csv',sep='')) )
	#	numbers
	tmp			<- subset(rpa, CLASS=='ancestral in either direction\nor intermingled')
	cat('\nNumber of couples',tmp[, length(unique(COUPID))])
	setkey(tmp, COUPID, MALE_SANGER_ID, FEMALE_SANGER_ID)
	cat('\nNumber of sequence pairs from couples',nrow(unique(tmp)))
	setkey(tmp, RUN, COUPID, PAIR_ID)
	cat('\nNumber of pairings (including repeated sequence pairs)',nrow(unique(tmp)))	
	#
	#	plot on proportion of assignments in epidemiologically possible direction
	#
	rpb		<- subset(df, COUP_SC=='F->M')[, list(CLASS='prop ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='trans_fm'])/sum(WIN_OF_TYPE_P[TYPE=='trans_fm'|TYPE=='trans_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]	
	tmp		<- subset(df, COUP_SC=='M->F')[, list(CLASS='prop ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='trans_mf'])/sum(WIN_OF_TYPE_P[TYPE=='trans_fm'|TYPE=='trans_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpb		<- rbind(rpb, tmp)
	tmp		<- rpb[order(CLASS, -CLASS_PROP),][, list(PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)), by='CLASS']
	rpb		<- merge(rpb, tmp, by=c('CLASS','PAIR_ID'))	
	setkey(rpb, CLASS, CLASS_RANK)
	ggplot(rpb, aes(x=CLASS_RANK, y=cumsum(CLASS_PROP))) + geom_line() + geom_point() +
			coord_cartesian(xlim=c(0, max(rpb[,CLASS_RANK])), ylim=c(0,max(rpb[,CLASS_RANK]))) +
			geom_abline(intercept=0, slope=0.5, colour='blue') +
			geom_abline(intercept=0, slope=1, colour='blue') +
			geom_text(aes(label=PAIR_ID), size=2, nudge_x=-.2, nudge_y=.8) +
			scale_y_continuous(expand=c(0,0)) +
			scale_x_continuous(expand=c(0,0)) +
			theme_bw() +
			labs(	x='\nsequence pairs\nof serodiscordant couples whose uninfected partners turns positive\n(cumulated)',
					y='# ancestral assignments in direction that is epidemiologically possible\nout of all ancestral assignments\n(cumulated)')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-only_possible_direction_assigned.pdf',sep='')), w=7, h=7)
	#
	#	plot on number of ancestral windows 
	#
	df[, WIN_OF_TYPE_N:= WIN_OF_TYPE_P*WIN_TOTAL]
	rpa		<- subset(df, COUP_SC=='F->M')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='trans_fm']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='trans_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(df, COUP_SC=='M->F')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='trans_mf']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='trans_fm'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= IN_DIR+AGAINST_DIR]
	rpa		<- melt(rpa, measure.vars=c('IN_DIR', 'AGAINST_DIR'))
	set(rpa, rpa[, which(variable=='IN_DIR')], 'variable', 'trm assignment in the only epidemiologically possible direction')
	set(rpa, rpa[, which(variable=='AGAINST_DIR')], 'variable', 'trm assignment against the only epidemiologically possible direction')
	setkey(rpa, PAIR_ID)	
	tmp		<- unique(rpa)[order(-WIN_TRM),][, list(COUPID=COUPID, PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(PAIR_ID, ' (', COUPID, ')', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PAIR_ID','COUPID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set1') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			#scale_x_discrete(labels=rpa[,PAIR_ID]) +
			labs(	x= '\nsequence pairs of couples', 
					y='number of transmission windows with direction\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trm_windows.pdf',sep='')), w=10, h=7)
	#
	#	plot on all windows 
	#	
	rpa		<- subset(df, COUP_SC=='F->M')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int'|TYPE=='cher']), NO_ANC=sum(WIN_OF_TYPE_N[!(TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int'|TYPE=='cher')])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(df, COUP_SC=='M->F')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='trans_mf'|TYPE=='trans_mf'|TYPE=='int'|TYPE=='cher']), NO_ANC=sum(WIN_OF_TYPE_N[!(TYPE=='trans_mf'|TYPE=='trans_mf'|TYPE=='int'|TYPE=='cher')])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= ANY_ANC]
	rpa		<- melt(rpa, measure.vars=c('ANY_ANC', 'NO_ANC'))
	set(rpa, rpa[, which(variable=='ANY_ANC')], 'variable', 'transmission assignments or intermingled or cherry')
	set(rpa, rpa[, which(variable=='NO_ANC')], 'variable', 'unint or disconnected')
	setkey(rpa, PAIR_ID)	
	tmp		<- unique(rpa)[order(-WIN_TRM),][, list(COUPID=COUPID, PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(PAIR_ID, ' (', COUPID, ')', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PAIR_ID','COUPID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set2') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			#scale_x_discrete(labels=rpa[,PAIR_ID]) +
			labs(	x= '\nsequence pairs of couples', 
					y='number of windows\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_other_windows.pdf',sep='')), w=10, h=7)
	
	
	
	
	#
	#	inconsistent pairings of same sequences
	#	inconsistent pairings of same couples
	
	#}
	#
	#	re-examine phylogenies for all sero-discordant couples
	#	
	tmp		<- subset(df, COUP_SC=='F->M' | COUP_SC=='M->F')
	setkey(tmp, RUN, PAIR_ID, MALE_SANGER_ID, FEMALE_SANGER_ID)
	dir		<- tmp[1,DIR]
	cat('\ndir is',dir,'\trun is',run)			
	setkey(tmp, MALE_SANGER_ID, FEMALE_SANGER_ID, PAIR_ID)				
	rpoku	<- unique(tmp)
	load( tmp[1, FILE] )	#loads phs dtrms dtrees
	invisible(sapply(seq_len(nrow(rpoku)), function(ii)
					{								
						pair.id		<- rpoku[ii, PAIR_ID]
						pty.run		<- rpoku[ii, PTY_RUN]
						dfs			<- subset(dtrees, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
						dfs[, TITLE:= dfs[, paste('pair', pair.id, '\n', rpoku[ii, COUP_SC], '\nid M: ', rpoku[ii, MALE_RID], ' (', rpoku[ii, MALE_SANGER_ID], ')\nid F: ', rpoku[ii, FEMALE_RID], ' (', rpoku[ii, FEMALE_SANGER_ID], ')\nrun ', pty.run, '\nwindow ', W_FROM,'-', W_TO,sep='')]]			
						plot.file	<- file.path(dir, paste(run,'-phsc-serodiscpairs-',rpoku[ii, COUP_SC],'-M-', rpoku[ii, MALE_RID],'-F-',rpoku[ii, FEMALE_RID],'-', pair.id,'.pdf',sep=''))			
						invisible(phsc.plot.selected.pairs(phs, dfs, rpoku[ii, MALE_SANGER_ID], rpoku[ii, FEMALE_SANGER_ID], plot.file=plot.file, pdf.h=150, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40))
					}))				
	
}

RakaiCouples.process.couples.161027<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )	
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	setkey(rp, COUPID)
	unique(rp)[, table(COUP_SC)]
	#	total of couples with assigned SANGER_ID in RCCS
	#   F->M    M->F seroinc seropos 
    # 	10      19      52     266 	
	tmp		<- merge(unique(rp), unique(subset(pty.runs, COUPID!='Other', COUPID)), by='COUPID')	
	tmp[, table(COUP_SC)]
	#	total of couples with SANGER_ID for which I have data
	#	F->M    M->F seroinc seropos 
	#	7       19   50      251 			total: 327
	#
	#	collect runs
	#
	infiles	<- data.table(	FILE= c(	'~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161027/RCCS_161027_w270_d20_phscout.rda',
										'~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161027/RCCS_161027_w270_d50_phscout.rda',
										'~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161027/RCCS_161027_w270_d200_phscout.rda'))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_phscout.rda','',basename(FILE))]				
	#
	#	for each run: get list of pairs
	#	
	rps		<- infiles[, {
				#FILE	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161027/RCCS_161027_w270_d20_phscout.rda'
				cat('\n',FILE)
				load(FILE)	#loads phs dtrms dtrees
				#	select likely pairs -- these are ordered
				dtrms[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]			
				tmp		<- dtrms[, list(CH=WIN_TOTAL-sum(WIN_OF_TYPE)), by=c('PTY_RUN','ID1','ID2')]
				stopifnot( tmp[, all(CH==0)] )
				tmp		<- dcast.data.table(dtrms, PAIR_ID~TYPE, value.var='WIN_OF_TYPE_P')
				set(tmp, tmp[, which(is.na(disconnected))], 'disconnected', 0)
				set(tmp, tmp[, which(is.na(int))], 'int', 0)
				set(tmp, tmp[, which(is.na(unint))], 'unint', 0)
				set(tmp, tmp[, which(is.na(cher))], 'cher', 0)
				set(tmp, tmp[, which(is.na(trans_12))], 'trans_12', 0)
				set(tmp, tmp[, which(is.na(trans_21))], 'trans_21', 0)
				tmp		<- melt.data.table(tmp, id.vars='PAIR_ID', value.name='WIN_OF_TYPE_P', variable.name='TYPE')				
				#
				dtrms	<- merge(unique(subset(dtrms, select=c(PAIR_ID, ID1, ID2, PTY_RUN, WIN_TOTAL, SCORE))), tmp, by='PAIR_ID')
				#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
				tmp		<- copy(dtrms)
				setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
				set(tmp, tmp[, which(TYPE=='trans_12')], 'TYPE', 'trans')
				set(tmp, tmp[, which(TYPE=='trans_21')], 'TYPE', 'trans_12')
				set(tmp, tmp[, which(TYPE=='trans')], 'TYPE', 'trans_21')
				dtrms	<- rbind(tmp, dtrms)
				#	first individual is always male	
				setnames(dtrms, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
				set(dtrms, dtrms[, which(TYPE=='trans_12')], 'TYPE', 'trans_mf')
				set(dtrms, dtrms[, which(TYPE=='trans_21')], 'TYPE', 'trans_fm')
				set(dtrms, NULL, 'TYPE', dtrms[, as.character(TYPE)])
				dtrms	<- merge(dtrms, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))				
				dtrms
			}, by=c('RUN','DIR','FILE')]
	save(rps, file= '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161027/RCCS_161027_w270_assignments.rda')
	rpsn	<- copy(rps)
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161007/RCCS_161007_w270_assignments.rda')
	rps		<- rbind(rpsn, rps)
	#
	cat('\nNumber of couples',paste(rps[, length(unique(COUPID)), by='RUN'], collapse=''))
	setkey(rps, RUN, COUPID, MALE_SANGER_ID, FEMALE_SANGER_ID)
	cat('\nNumber of sequence pairs from couples',nrow(unique(rps)))
	setkey(rps, RUN, COUPID, PAIR_ID)
	cat('\nNumber of pairings (including repeated sequence pairs)',nrow(unique(rps)))
	#
	#	check inconsistent
	#
	if(0)
	{
		rpsn	<- copy(rps)
		rps.new<- unique(subset(rpsn, RUN=='RCCS_161027_w270_d200', select=c(COUPID,MALE_TAXA,FEMALE_TAXA)))
		rps.new[, TYPE:='NEW']	
		rps.old<- unique(subset(rps, RUN=='RCCS_161007_w270', select=c(COUPID,MALE_TAXA,FEMALE_TAXA)))
		rps.old[, TYPE:='OLD']
		tmp		<- merge(rps.new, rps.old, by=c('COUPID','MALE_TAXA','FEMALE_TAXA'),all=1)
		subset(tmp, is.na(TYPE.x))
		subset(pty.runs, TAXA=='PG14-UG503118-S05122' | TAXA=='PG14-UG503081-S05085')
		subset(pty.runs, TAXA=='PG14-UG500525-S02529' | TAXA=='PG14-UG500526-S02530')
		
		subset(rps, MALE_TAXA=='PG14-UG503118-S05122' | FEMALE_TAXA=='PG14-UG503081-S05085')
	}
	set(rps, NULL, 'RUN', rps[, factor(RUN, levels=c("RCCS_161007_w270","RCCS_161027_w270_d20","RCCS_161027_w270_d50","RCCS_161027_w270_d200"))])
	#
	#	for each run: plot pairs	
	#	
	run		<- 'RCCS_161027_w270_dxxx'
	dir		<- rps$DIR[1]
	rpp		<- subset(rps, RUN==rps$RUN[1])
	setkey(rpp, MALE_SANGER_ID, FEMALE_SANGER_ID, PAIR_ID)
	tmp		<- unique(rpp)	
	tmp[, PLOT_ID:= as.numeric(gsub('-','\\.',PAIR_ID))]
	tmp		<- tmp[order(-PLOT_ID),]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair', PAIR_ID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- merge(subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL)), rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('trans_mf','trans_fm','unint','int','cher','disconnected'), labels=c('M transmit to F','F transmit to M','M, F are unint','M, F are intermingled','M, F are a cherry','M, F are disconnected'))])
	tmp		<- subset(tmp, select=c(RUN, PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, WIN_OF_TYPE_P, WIN_TOTAL, TYPE, LABEL))
	tmp[, WIN_OF_TYPE_N:=WIN_OF_TYPE_P*WIN_TOTAL]
	setkey(tmp, COUP_SC, LABEL, RUN, TYPE)
	#	F->M
	tmp2	<- subset(tmp, COUP_SC=='F->M')
	ggplot(tmp2, aes(x=RUN, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC+LABEL, ncol=1) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_F2M.pdf',sep='')), w=25, h=max(4,0.23*nrow(tmp2)), limitsize = FALSE)
	#	M->F
	tmp2	<- subset(tmp, COUP_SC=='M->F')
	ggplot(tmp2, aes(x=RUN, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC+LABEL, ncol=1) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_M2F.pdf',sep='')), w=25, h=max(4,0.23*nrow(tmp2)), limitsize = FALSE)
	#	seroinc
	tmp2	<- subset(tmp, COUP_SC=='seroinc')
	ggplot(tmp2, aes(x=RUN, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC+LABEL, ncol=1) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_SeroInc.pdf',sep='')), w=25, h=max(4,0.23*nrow(tmp2)), limitsize = FALSE)
	#	seropos
	tmp2	<- subset(tmp, COUP_SC=='seropos')
	ggplot(tmp2, aes(x=RUN, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC+LABEL, ncol=1) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_SeroPos.pdf',sep='')), w=25, h=max(4,0.23*nrow(tmp2)), limitsize = FALSE)
	#
	#	plot number of directional trm assignments in only possible direction 
	#
	rps[, WIN_OF_TYPE_N:= WIN_OF_TYPE_P*WIN_TOTAL]
	rpa		<- subset(rps, COUP_SC=='F->M')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='trans_fm']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='trans_mf'])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(rps, COUP_SC=='M->F')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='trans_mf']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='trans_fm'])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= IN_DIR+AGAINST_DIR]
	rpa		<- melt(rpa, measure.vars=c('IN_DIR', 'AGAINST_DIR'))
	set(rpa, rpa[, which(variable=='IN_DIR')], 'variable', 'trm assignment in the only epidemiologically possible direction')
	set(rpa, rpa[, which(variable=='AGAINST_DIR')], 'variable', 'trm assignment against the only epidemiologically possible direction')	
	setkey(rpa, PAIR_ID)		
	tmp		<- unique(subset(rpa, RUN==rpa$RUN[1]))[order(-WIN_TRM),][, list(COUPID=COUPID, MALE_SANGER_ID=MALE_SANGER_ID, FEMALE_SANGER_ID=FEMALE_SANGER_ID, PTY_RUN=PTY_RUN, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PTY_RUN','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set1') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~RUN) +
			coord_flip() +
			labs(	x= '\nsequence pairs of couples', 
					y='number of transmission windows with direction\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trm_windows.pdf',sep='')), w=20, h=10)
	#
	#
	#
	rpa		<- subset(rps, COUP_SC=='F->M')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int']), CHER=WIN_OF_TYPE_N[TYPE=='cher'], NO_ANC=sum(WIN_OF_TYPE_N[!(TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int'|TYPE=='cher')])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(rps, COUP_SC=='M->F')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='trans_mf'|TYPE=='trans_mf'|TYPE=='int']), CHER=WIN_OF_TYPE_N[TYPE=='cher'], NO_ANC=sum(WIN_OF_TYPE_N[!(TYPE=='trans_mf'|TYPE=='trans_mf'|TYPE=='int'|TYPE=='cher')])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= ANY_ANC]
	rpa		<- melt(rpa, measure.vars=c('ANY_ANC', 'NO_ANC','CHER'))
	set(rpa, rpa[, which(variable=='ANY_ANC')], 'variable', 'transmission assignments or intermingled')
	set(rpa, rpa[, which(variable=='NO_ANC')], 'variable', 'unint or disconnected')
	set(rpa, rpa[, which(variable=='CHER')], 'variable', 'cherry')
	setkey(rpa, PAIR_ID)	
	tmp		<- unique(subset(rpa, RUN==rpa$RUN[1]))[order(-WIN_TRM),][, list(COUPID=COUPID, MALE_SANGER_ID=MALE_SANGER_ID, FEMALE_SANGER_ID=FEMALE_SANGER_ID, PTY_RUN=PTY_RUN, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PTY_RUN','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set2') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~RUN) +
			coord_flip() +
			labs(	x= '\nsequence pairs of couples', 
					y='number of windows\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_other_windows.pdf',sep='')), w=20, h=10)
	
	#
	#	inspect rogue case by window
	#
	dr		<- as.data.table(read.csv('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161027/couple_374_rogues_OR.csv'))
	setnames(dr, c('window.start','J110207','G105508'), c('W_FROM','15778_1_82','15758_1_76'))
	tmp		<- melt(subset(dr, select=c('W_FROM','15778_1_82','15758_1_76')), id.vars='W_FROM',value.name='BRL',variable.name='ID')  
	dr		<- melt(subset(dr, select=c('W_FROM','rogue_15778_1_82','rogue_15758_1_76')), id.vars='W_FROM',value.name='ROGUE',variable.name='ID')
	set(dr, NULL, 'ID', dr[, gsub('rogue_','',ID)])
	dr		<- merge(tmp, dr, by=c('W_FROM','ID'))
	
	#infiles<- c(	'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d20_rerun/ptyr66_trmStatsPerWindow.rda',
	#				'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d50_rerun/ptyr66_trmStatsPerWindow.rda',
	#				'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d200_rerun/ptyr66_trmStatsPerWindow.rda'	)
	infiles<- c(	'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160915_couples_w270/ptyr66_trmStatsPerWindow.rda')
	id.m	<- '15778_1_82'
	id.f	<- '15758_1_76'	
	df		<- phsc.get.assignments.by.window.for.couple(id1=paste(id.m,'.bam',sep=''), id2=paste(id.f,'.bam',sep=''), infiles)
	set(df, NULL, 'ID1', df[, gsub('\\.bam','',ID1)])
	set(df, NULL, 'ID2', df[, gsub('\\.bam','',ID2)])
	setnames(df, colnames(df), gsub('_rerun','',gsub('couples_','',gsub('Rakai_ptoutput_','',colnames(df)))))
	
	setnames(dr, c('ID','BRL','ROGUE'), c('ID1','BRL1','ROGUE1'))
	df		<- merge(df,dr,by=c('ID1','W_FROM'))
	setnames(dr, c('ID1','BRL1','ROGUE1'), c('ID2','BRL2','ROGUE2'))
	df		<- merge(df,dr,by=c('ID2','W_FROM'))
	df		<- subset(df, !is.na(BRL1) & !is.na(BRL2) & !is.na(ROGUE1) & !is.na(ROGUE2))
	#
	#	branches among rogues and non rogues
	#
	tmp		<- subset(df, select=c(BRL1,ROGUE1))
	setnames(tmp, c('BRL1','ROGUE1'), c('BRL','ROGUE'))
	tmp2	<- subset(df, select=c(BRL2,ROGUE2))
	setnames(tmp2, c('BRL2','ROGUE2'), c('BRL','ROGUE'))	
	tmp		<- rbind(tmp, tmp2)
	ggplot(tmp, aes(x=BRL, fill=factor(ROGUE))) + geom_histogram() + facet_grid(~ROGUE)
	
	subset(df, (BRL1>0.05 & !ROGUE1) | (BRL2>0.05 & !ROGUE2))
	subset(df, (BRL1<0.07 & ROGUE1) | (BRL2<0.07 & ROGUE2))
	#	calculate patristic distance..
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161027/RCCS_161027_w270_d200_phscout.rda')
	#ph		<- phs[[ subset(dtrees, PTY_RUN==66 & W_FROM==1225)[, IDX] ]]
	#tmp		<- as.matrix(cophenetic.phylo(ph))
	#max( tmp['15758_1_76.bam_read_7_count_1', grepl('15758_1_76', colnames(tmp))] )
	#	collect branch lengths from individuals that I classified manually
	tmp		<- merge(subset(dtrees, PTY_RUN==66, c(W_FROM, IDX)), subset(df, select=c(ID1, W_FROM, ROGUE1)),by='W_FROM')
	setnames(tmp, c('ID1','ROGUE1'), c('ID','ROGUE'))
	dr		<- tmp[, {
				#	IDX<- 26801; ID1<- '15778_1_82'
				ph	<- phs[[ IDX ]]
				z	<- ph$tip.label[!grepl(ID,ph$tip.label)]
				ph	<- drop.tip(ph,z)
				list(BRL=ph$edge.length) 
			}, by=c('W_FROM','ID','ROGUE')]	
	tmp		<- merge(subset(dtrees, PTY_RUN==66, c(W_FROM, IDX)), subset(df, select=c(ID2, W_FROM, ROGUE2)),by='W_FROM')
	setnames(tmp, c('ID2','ROGUE2'), c('ID','ROGUE'))
	tmp		<- tmp[, {
				#	IDX<- 26801; ID1<- '15778_1_82'
				ph	<- phs[[ IDX ]]
				z	<- ph$tip.label[!grepl(ID,ph$tip.label)]
				ph	<- drop.tip(ph,z)
				list(BRL=ph$edge.length) 
			}, by=c('W_FROM','ID','ROGUE')]
	dr		<- rbind(dr, tmp)	
	#	calculate prob that the max BRL is > x under the "null" that all distances are from the same distribution
	#	using Weibull as "null" model because (1) for x positive and (2) it is easy to calculate 1-CDF(max X_i))	
	require(gamlss)
	dp		<- dr[, {
				x			<- max(BRL)
				p			<- NA_real_
				z			<- BRL[ BRL>1e-5 ]
				if(length(z)>3)	# don't think makes sense to fit Weibull to 3 data points
				{	
					cat('\n', W_FROM, ID) # for debugging					 
					w		<- gamlss(formula=z~1, family=WEI, trace=1)
					w.l		<- exp(coef(w, what='mu'))
					w.k		<- exp(coef(w, what='sigma'))
					p			<- 1 - ( 1 - exp( -( max(z)/w.l )^w.k ) )^length(z)					
				}
				list(P=p, BRL.mx=max(BRL))
			}, by=c('W_FROM','ID','ROGUE')]	
	#
	#	define rogue as all reads that are not in largest subtree
	#
	require(phangorn)
	tmp		<- merge(subset(dtrees, PTY_RUN==66, c(W_FROM, IDX)), subset(df, select=c(ID1, W_FROM, ROGUE1)),by='W_FROM')
	setnames(tmp, c('ID1','ROGUE1'), c('ID','ROGUE'))
	dst		<- tmp[, {
				#	IDX<- 26817; ID<- '15778_1_82'
				ph		<- phs[[ IDX ]]
				tips.id	<- which(grepl(ID,ph$tip.label))
				rogue	<- rep(FALSE, length(tips.id))
				z		<- Ancestors(ph, tips.id)
				#	along each path find last with taxa==ID only
				z		<- sapply(seq_along(z), function(i){
							zz	<- Descendants(ph, z[[i]], type="tips")							 
							zzz	<- sapply(zz, function(x) all(grepl(ID,ph$tip.label[x])) )
							zzz	<- which(zzz)
							ifelse(length(zzz), z[[i]][zzz[1]], tips.id[i])							
						})
				z		<- unique(z)
				#	determine number of reads in each tip
				zz		<- sapply(Descendants(ph, z, type="tips"), function(x) sum(as.numeric(gsub('.*_count_([0-9]+)','\\1',ph$tip.label[x])))	)
				if(length(zz)>1 && any(zz/sum(zz)<.01))
				{
					rogue.cl<- z[which( zz/sum(zz)<.01 )]
					rogue.id<- unlist(Descendants(ph, rogue.cl, type="tips"))
					rogue	<- sapply(tips.id, function(x) any(x==rogue.id))					
				}
				list(TAXA=ph$tip.label[tips.id], ROGUE=rogue)				
			}, by=c('W_FROM','ID','LARGEST_SUBTREE_ROGUE')]	
	#	
	#

	dp[, BRL.mx.c:= cut(BRL.mx, breaks=c(0,0.04,Inf), labels=c('<4%','>=4%'))]		
	subset(dp, !is.na(P))[, table(ROGUE, BRL.mx.c)]
	subset(dp, !is.na(P))[, table(ROGUE, P<.05)]
	subset(dp, !is.na(P))[, table(ROGUE, P<.0001)]	
	ggplot(subset(dp, !is.na(P)), aes(x=factor(ROGUE), y=P)) + 
			geom_boxplot() + 
			facet_grid(.~BRL.mx.c) +
			labs(x='manual classification rogues',y='prob that max BRL is > observed max\nunder Weibull model')
	
	#
	#	inspect couple K061956:J061939 ( M:15103_1_74 F:15861_1_22 run:91 )
	#
	infiles<- c(	'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d20_rerun/ptyr91_trmStatsPerWindow.rda',
					'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d50_rerun/ptyr91_trmStatsPerWindow.rda',
					'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d200_rerun/ptyr91_trmStatsPerWindow.rda'	)
	id.m	<- '15103_1_74'
	id.f	<- '15861_1_22'	
	df		<- phsc.get.assignments.by.window.for.couple(id1=paste(id.m,'.bam',sep=''), id2=paste(id.f,'.bam',sep=''), infiles)
	set(df, NULL, 'ID1', df[, gsub('\\.bam','',ID1)])
	set(df, NULL, 'ID2', df[, gsub('\\.bam','',ID2)])
	setnames(df, colnames(df), gsub('_rerun','',gsub('couples_','',gsub('Rakai_ptoutput_','',colnames(df)))))

	infiles<- c(	'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d20_rerun/ptyr62_trmStatsPerWindow.rda',
					'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d50_rerun/ptyr62_trmStatsPerWindow.rda',
					'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161027_couples_w270_d200_rerun/ptyr62_trmStatsPerWindow.rda'	)
	id.m	<- '15861_1_26'
	id.f	<- '15861_1_22'	
	df		<- phsc.get.assignments.by.window.for.couple(id1=paste(id.m,'.bam',sep=''), id2=paste(id.f,'.bam',sep=''), infiles)
	set(df, NULL, 'ID1', df[, gsub('\\.bam','',ID1)])
	set(df, NULL, 'ID2', df[, gsub('\\.bam','',ID2)])
	setnames(df, colnames(df), gsub('_rerun','',gsub('couples_','',gsub('Rakai_ptoutput_','',colnames(df)))))
}

RakaiCouples.process.couples.161107<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )	
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	setkey(rp, COUPID)
	unique(rp)[, table(COUP_SC)]
	#	total of couples with assigned SANGER_ID in RCCS
	#   F->M    M->F seroinc seropos 
	# 	10      19      52     266 	
	tmp		<- merge(unique(rp), unique(subset(pty.runs, COUPID!='Other', COUPID)), by='COUPID')	
	tmp[, table(COUP_SC)]
	#	total of couples with SANGER_ID for which I have data
	#	F->M    M->F seroinc seropos 
	#	7       19   50      251 			total: 327
	
	load('~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161107/RCCS_161027_w270d200_r004_mr1_phscout.rda')
	dwin[, RUN:= 'mr1']
	dwin[, min(ID1_R)]
	#
	#	for each run: save trm assignments for couples
	#
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161107', pattern='_phscout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_phscout.rda','',basename(FILE))]
	infiles	<- subset(infiles, grepl('d50',RUN))
	invisible(infiles[, {
				#FILE	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161107/RCCS_161027_w270d20_r004_mr1_phscout.rda'
				cat('\n',FILE)
				load(FILE)	#loads phs dtrms dtrees dwin
				#
				#	extract couples transmissions summary
				#
				#	check
				tmp		<- dtrms[, list(CH=WIN_TOTAL-sum(WIN_OF_TYPE)), by=c('PTY_RUN','ID1','ID2')]
				stopifnot( tmp[, all(CH==0)] )
				#	add zeros
				tmp		<- dcast.data.table(dtrms, PAIR_ID~TYPE, value.var='WIN_OF_TYPE')
				set(tmp, tmp[, which(is.na(disconnected))], 'disconnected', 0)
				set(tmp, tmp[, which(is.na(int))], 'int', 0)
				set(tmp, tmp[, which(is.na(unint))], 'unint', 0)
				set(tmp, tmp[, which(is.na(cher))], 'cher', 0)
				set(tmp, tmp[, which(is.na(trans_12))], 'trans_12', 0)
				set(tmp, tmp[, which(is.na(trans_21))], 'trans_21', 0)
				tmp		<- melt.data.table(tmp, id.vars='PAIR_ID', value.name='WIN_OF_TYPE', variable.name='TYPE')								
				dtrms	<- merge(unique(subset(dtrms, select=c(PAIR_ID, ID1, ID2, ID1_R, ID1_L, ID2_R, ID2_L, PTY_RUN, WIN_TOTAL, SCORE))), tmp, by='PAIR_ID')
				#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
				tmp		<- copy(dtrms)
				setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
				set(tmp, tmp[, which(TYPE=='trans_12')], 'TYPE', 'trans')
				set(tmp, tmp[, which(TYPE=='trans_21')], 'TYPE', 'trans_12')
				set(tmp, tmp[, which(TYPE=='trans')], 'TYPE', 'trans_21')
				dtrms	<- rbind(tmp, dtrms)
				#	select couples; first individual is always male	
				setnames(dtrms, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
				set(dtrms, dtrms[, which(TYPE=='trans_12')], 'TYPE', 'trans_mf')
				set(dtrms, dtrms[, which(TYPE=='trans_21')], 'TYPE', 'trans_fm')
				set(dtrms, NULL, 'TYPE', dtrms[, as.character(TYPE)])
				dtrms	<- merge(dtrms, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
				#	
				dtrms[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]
				set(dtrms, NULL, 'RUN', RUN)
				save(dtrms, file=gsub('phscout.rda','trmsout.rda',FILE))
				#
				#	extract couples transmissions assignments per window
				#
				#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
				tmp		<- copy(dwin)
				setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
				set(tmp, tmp[, which(TYPE=='trans_12')], 'TYPE', 'trans')
				set(tmp, tmp[, which(TYPE=='trans_21')], 'TYPE', 'trans_12')
				set(tmp, tmp[, which(TYPE=='trans')], 'TYPE', 'trans_21')
				dwin	<- rbind(tmp, dwin)
				#	select couples; first individual is always male	
				setnames(dwin, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
				set(dwin, dwin[, which(TYPE=='trans_12')], 'TYPE', 'trans_mf')
				set(dwin, dwin[, which(TYPE=='trans_21')], 'TYPE', 'trans_fm')
				set(dwin, NULL, 'TYPE', dwin[, as.character(TYPE)])
				dwin	<- merge(dwin, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
				set(dwin, NULL, 'RUN', RUN)
				save(dwin, file=gsub('phscout.rda','trmwout.rda',FILE))				
			}, by=c('RUN','DIR','FILE')])
	#		
	#	load transmission summary assignments		
	#		
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161107', pattern='_trmsout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_trmsout.rda','',basename(FILE))]					
	rps		<- infiles[, {
				load(FILE)
				dtrms
			}, by=c('DIR','FILE')]
	save(rps, file= '~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161107/RCCS_161107_w270_trms_assignments.rda')
	#		
	#	load transmission window assignments		
	#		
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161107', pattern='_trmwout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_trmwout.rda','',basename(FILE))]					
	rpw		<- infiles[, {
				load(FILE)
				dwin
			}, by=c('DIR','FILE')]
	save(rpw, file= '~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161107/RCCS_161107_w270_trmw_assignments.rda')	
}

RakaiCouples.process.couples.161110<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )	
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	setkey(rp, COUPID)
	unique(rp)[, table(COUP_SC)]
	#	total of couples with assigned SANGER_ID in RCCS
	#   F->M    M->F seroinc seropos 
	# 	10      19      52     266 	
	tmp		<- merge(unique(rp), unique(subset(pty.runs, COUPID!='Other', COUPID)), by='COUPID')	
	tmp[, table(COUP_SC)]
	#	total of couples with SANGER_ID for which I have data
	#	F->M    M->F seroinc seropos 
	#	7       19   50      251 			total: 327
	
	#
	#	for each run: save trm assignments for couples
	#
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161110', pattern='_phscout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_phscout.rda','',basename(FILE))]
	#infiles	<- subset(infiles, grepl('d50',RUN))
	invisible(infiles[, {
						#FILE	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161107/RCCS_161027_w270d20_r004_mr1_phscout.rda'
						cat('\n',FILE)
						load(FILE)	#loads phs dtrms dtrees dwin
						#
						#	extract couples transmissions summary
						#
						#	check
						tmp		<- dtrms[, list(CH=WIN_TOTAL-sum(WIN_OF_TYPE)), by=c('PTY_RUN','ID1','ID2')]
						stopifnot( tmp[, all(CH==0)] )
						#	add zeros
						tmp		<- dcast.data.table(dtrms, PAIR_ID~TYPE, value.var='WIN_OF_TYPE')
						set(tmp, tmp[, which(is.na(disconnected))], 'disconnected', 0)
						set(tmp, tmp[, which(is.na(int))], 'int', 0)
						set(tmp, tmp[, which(is.na(unint))], 'unint', 0)
						set(tmp, tmp[, which(is.na(cher))], 'cher', 0)
						set(tmp, tmp[, which(is.na(trans_12))], 'trans_12', 0)
						set(tmp, tmp[, which(is.na(trans_21))], 'trans_21', 0)
						tmp		<- melt.data.table(tmp, id.vars='PAIR_ID', value.name='WIN_OF_TYPE', variable.name='TYPE')								
						dtrms	<- merge(unique(subset(dtrms, select=c(PAIR_ID, ID1, ID2, ID1_R, ID1_L, ID2_R, ID2_L, PTY_RUN, WIN_TOTAL, SCORE))), tmp, by='PAIR_ID')
						#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
						tmp		<- copy(dtrms)
						setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
						set(tmp, tmp[, which(TYPE=='trans_12')], 'TYPE', 'trans')
						set(tmp, tmp[, which(TYPE=='trans_21')], 'TYPE', 'trans_12')
						set(tmp, tmp[, which(TYPE=='trans')], 'TYPE', 'trans_21')
						dtrms	<- rbind(tmp, dtrms)
						#	select couples; first individual is always male	
						setnames(dtrms, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
						set(dtrms, dtrms[, which(TYPE=='trans_12')], 'TYPE', 'trans_mf')
						set(dtrms, dtrms[, which(TYPE=='trans_21')], 'TYPE', 'trans_fm')
						set(dtrms, NULL, 'TYPE', dtrms[, as.character(TYPE)])
						dtrms	<- merge(dtrms, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
						#	
						dtrms[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]
						set(dtrms, NULL, 'RUN', RUN)
						save(dtrms, file=gsub('phscout.rda','trmsout.rda',FILE))
						#
						#	extract couples transmissions assignments per window
						#
						#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
						tmp		<- copy(dwin)
						setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
						set(tmp, tmp[, which(TYPE=='trans_12')], 'TYPE', 'trans')
						set(tmp, tmp[, which(TYPE=='trans_21')], 'TYPE', 'trans_12')
						set(tmp, tmp[, which(TYPE=='trans')], 'TYPE', 'trans_21')
						dwin	<- rbind(tmp, dwin)
						#	select couples; first individual is always male	
						setnames(dwin, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
						set(dwin, dwin[, which(TYPE=='trans_12')], 'TYPE', 'trans_mf')
						set(dwin, dwin[, which(TYPE=='trans_21')], 'TYPE', 'trans_fm')
						set(dwin, NULL, 'TYPE', dwin[, as.character(TYPE)])
						dwin	<- merge(dwin, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
						set(dwin, NULL, 'RUN', RUN)
						save(dwin, file=gsub('phscout.rda','trmwout.rda',FILE))				
					}, by=c('RUN','DIR','FILE')])
	#		
	#	load transmission summary assignments		
	#		
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161110', pattern='_trmsout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_trmsout.rda','',basename(FILE))]					
	rps		<- infiles[, {
				load(FILE)
				dtrms
			}, by=c('DIR','FILE')]
	save(rps, file= '~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161110/RCCS_161110_w270_trms_assignments.rda')
	#		
	#	load transmission window assignments		
	#		
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161110', pattern='_trmwout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_trmwout.rda','',basename(FILE))]					
	rpw		<- infiles[, {
				load(FILE)
				dwin
			}, by=c('DIR','FILE')]
	save(rpw, file= '~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161110/RCCS_161110_w270_trmw_assignments.rda')	
}

RakaiCouples.process.couples.161213<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )	
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	setkey(rp, COUPID)
	unique(rp)[, table(COUP_SC)]
	#	total of couples with assigned SANGER_ID in RCCS
	#   F->M    M->F seroinc seropos 
	# 	10      19      52     266 	
	tmp		<- merge(unique(rp), unique(subset(pty.runs, COUPID!='Other', COUPID)), by='COUPID')	
	tmp[, table(COUP_SC)]
	#	total of couples with SANGER_ID for which I have data
	#	F->M    M->F seroinc seropos 
	#	7       19   50      251 			total: 327
	
	#
	#	for each run: save trm assignments for couples
	#
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161213', pattern='_phscout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_phscout.rda','',basename(FILE))]
	#infiles	<- subset(infiles, grepl('d50',RUN))
	invisible(infiles[, {
						#FILE	<- '/Users/Oliver/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161213/RCCS_161213_w270_d50_p001_mr20_mt2_cl1_phscout.rda'
						cat('\n',FILE)
						load(FILE)	#loads phs dtrms dtrees dwin
						#
						#	extract couples transmissions summary
						#
						#	check
						tmp		<- dtrms[, list(CH=WIN_TOTAL-sum(WIN_OF_TYPE)), by=c('PTY_RUN','ID1','ID2')]
						stopifnot( tmp[, all(CH==0)] )
						#	add zeros
						tmp		<- dcast.data.table(dtrms, PAIR_ID~TYPE, value.var='WIN_OF_TYPE')
						for(x in setdiff(colnames(tmp),'PAIR_ID'))
							set(tmp, which(is.na(tmp[[x]])), x, 0)						
						tmp		<- melt.data.table(tmp, id.vars='PAIR_ID', value.name='WIN_OF_TYPE', variable.name='TYPE')								
						dtrms	<- merge(unique(subset(dtrms, select=c(PAIR_ID, ID1, ID2, ID1_R, ID1_L, ID2_R, ID2_L, PTY_RUN, WIN_TOTAL, SCORE))), tmp, by='PAIR_ID')
						#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
						tmp		<- copy(dtrms)
						setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
						set(tmp, tmp[, which(TYPE=='anc_12')], 'TYPE', 'anc')
						set(tmp, tmp[, which(TYPE=='anc_21')], 'TYPE', 'anc_12')
						set(tmp, tmp[, which(TYPE=='anc')], 'TYPE', 'anc_21')
						set(tmp, tmp[, which(TYPE=='close_anc_12')], 'TYPE', 'anc')
						set(tmp, tmp[, which(TYPE=='close_anc_21')], 'TYPE', 'close_anc_12')
						set(tmp, tmp[, which(TYPE=='anc')], 'TYPE', 'close_anc_21')						
						dtrms	<- rbind(tmp, dtrms)
						#	select couples; first individual is always male	
						setnames(dtrms, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
						set(dtrms, dtrms[, which(TYPE=='anc_12')], 'TYPE', 'anc_mf')
						set(dtrms, dtrms[, which(TYPE=='anc_21')], 'TYPE', 'anc_fm')
						set(dtrms, dtrms[, which(TYPE=='close_anc_12')], 'TYPE', 'close_anc_mf')
						set(dtrms, dtrms[, which(TYPE=='close_anc_21')], 'TYPE', 'close_anc_fm')						
						set(dtrms, NULL, 'TYPE', dtrms[, as.character(TYPE)])
						dtrms	<- merge(dtrms, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
						#	
						dtrms[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]
						set(dtrms, NULL, 'RUN', RUN)
						save(dtrms, file=gsub('phscout.rda','trmsout.rda',FILE))
						#
						#	extract couples transmissions assignments per window
						#
						#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
						tmp		<- copy(dwin)
						setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
						set(tmp, tmp[, which(TYPE=='anc_12')], 'TYPE', 'anc')
						set(tmp, tmp[, which(TYPE=='anc_21')], 'TYPE', 'anc_12')
						set(tmp, tmp[, which(TYPE=='anc')], 'TYPE', 'anc_21')
						set(tmp, tmp[, which(TYPE=='close_anc_12')], 'TYPE', 'anc')
						set(tmp, tmp[, which(TYPE=='close_anc_21')], 'TYPE', 'close_anc_12')
						set(tmp, tmp[, which(TYPE=='anc')], 'TYPE', 'close_anc_21')												
						dwin	<- rbind(tmp, dwin)
						#	select couples; first individual is always male	
						setnames(dwin, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
						set(dwin, dwin[, which(TYPE=='anc_12')], 'TYPE', 'anc_mf')
						set(dwin, dwin[, which(TYPE=='anc_21')], 'TYPE', 'anc_fm')
						set(dwin, dwin[, which(TYPE=='close_anc_12')], 'TYPE', 'close_anc_mf')
						set(dwin, dwin[, which(TYPE=='close_anc_21')], 'TYPE', 'close_anc_fm')												
						set(dwin, NULL, 'TYPE', dwin[, as.character(TYPE)])
						dwin	<- merge(dwin, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
						set(dwin, NULL, 'RUN', RUN)
						save(dwin, file=gsub('phscout.rda','trmwout.rda',FILE))				
					}, by=c('RUN','DIR','FILE')])
	#		
	#	load transmission summary assignments		
	#		
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161213', pattern='_trmsout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_trmsout.rda','',basename(FILE))]					
	rps		<- infiles[, {
				load(FILE)
				dtrms
			}, by=c('DIR','FILE')]
	save(rps, file= '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161213/RCCS_161213_w270_trms_assignments.rda')
	#		
	#	load transmission window assignments		
	#		
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161213', pattern='_trmwout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_trmwout.rda','',basename(FILE))]					
	rpw		<- infiles[, {
				load(FILE)
				dwin
			}, by=c('DIR','FILE')]
	save(rpw, file= '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161213/RCCS_161213_w270_trmw_assignments.rda')	
}

RakaiCouples.process.couples.161219<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )	
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	setkey(rp, COUPID)
	unique(rp)[, table(COUP_SC)]
	#	total of couples with assigned SANGER_ID in RCCS
	#   F->M    M->F seroinc seropos 
	# 	10      19      52     266 	
	tmp		<- merge(unique(rp), unique(subset(pty.runs, COUPID!='Other', COUPID)), by='COUPID')	
	tmp[, table(COUP_SC)]
	#	total of couples with SANGER_ID for which I have data
	#	F->M    M->F seroinc seropos 
	#	7       19   50      251 			total: 327
	
	#
	#	for each run: save trm assignments for couples
	#
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219', pattern='_phscout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_phscout.rda','',basename(FILE))]
	#infiles	<- subset(infiles, grepl('d50',RUN))
	invisible(infiles[, {
						#FILE	<- '/Users/Oliver/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_d50_p001_mr20_mt1_cl1_d7_phscout.rda'
						cat('\n',FILE)
						load(FILE)	#loads phs dtrms dtrees dwin
						#
						#	extract couples transmissions summary
						#
						#	check
						tmp		<- dtrms[, list(CH=WIN_TOTAL-sum(WIN_OF_TYPE)), by=c('PTY_RUN','ID1','ID2')]
						stopifnot( tmp[, all(CH==0)] )
						#	keep association between TYPE and TYPE_PAIR
						dtypes	<- unique(subset(dtrms, select=c(TYPE, TYPE_PAIR)))
						set(dtypes, NULL, 'TYPE', dtypes[, gsub('12','mf',TYPE)])						
						set(dtypes, NULL, 'TYPE', dtypes[, gsub('21','fm',TYPE)])						
						#	add zeros
						tmp		<- dcast.data.table(dtrms, PAIR_ID~TYPE, value.var='WIN_OF_TYPE')
						for(x in setdiff(colnames(tmp),'PAIR_ID'))
							set(tmp, which(is.na(tmp[[x]])), x, 0)						
						tmp		<- melt.data.table(tmp, id.vars='PAIR_ID', value.name='WIN_OF_TYPE', variable.name='TYPE')								
						dtrms	<- merge(unique(subset(dtrms, select=c(PAIR_ID, ID1, ID2, ID1_R, ID1_L, ID2_R, ID2_L, PTY_RUN, WIN_TOTAL, SCORE)), by=c('ID1','ID2','PTY_RUN')), tmp, by='PAIR_ID')
						#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
						tmp		<- copy(dtrms)
						setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))						
						set(tmp, NULL, 'TYPE', tmp[, gsub('12','XX33XX',TYPE)])						
						set(tmp, NULL, 'TYPE', tmp[, gsub('21','12',TYPE)])
						set(tmp, NULL, 'TYPE', tmp[, gsub('XX33XX','21',TYPE)])						
						dtrms	<- rbind(tmp, dtrms)
						#	select couples; first individual is always male	
						setnames(dtrms, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
						set(dtrms, NULL, 'TYPE', dtrms[, gsub('12','mf',TYPE)])						
						set(dtrms, NULL, 'TYPE', dtrms[, gsub('21','fm',TYPE)])
						set(dtrms, NULL, 'TYPE', dtrms[, as.character(TYPE)])
						dtrms	<- merge(dtrms, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
						#	
						dtrms[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]
						set(dtrms, NULL, 'RUN', RUN)
						dtrms	<- merge(dtrms, dtypes, by='TYPE')
						save(dtrms, file=gsub('phscout.rda','trmsout.rda',FILE))
						#
						#	extract couples transmissions assignments per window
						#
						#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
						tmp		<- copy(dwin)
						setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
						set(tmp, NULL, 'TYPE', tmp[, gsub('12','XX33XX',TYPE)])						
						set(tmp, NULL, 'TYPE', tmp[, gsub('21','12',TYPE)])
						set(tmp, NULL, 'TYPE', tmp[, gsub('XX33XX','21',TYPE)])						
						dwin	<- rbind(tmp, dwin)
						#	select couples; first individual is always male	
						setnames(dwin, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
						set(dwin, NULL, 'TYPE', dwin[, gsub('12','mf',TYPE)])						
						set(dwin, NULL, 'TYPE', dwin[, gsub('21','fm',TYPE)])
						set(dwin, NULL, 'TYPE', dwin[, as.character(TYPE)])
						dwin	<- merge(dwin, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
						set(dwin, NULL, 'RUN', RUN)
						save(dwin, file=gsub('phscout.rda','trmwout.rda',FILE))				
					}, by=c('RUN','DIR','FILE')])
	#		
	#	save transmission summary assignments		
	#		
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219', pattern='_trmsout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_trmsout.rda','',basename(FILE))]					
	rps		<- infiles[, {
				load(FILE)
				dtrms
			}, by=c('DIR','FILE')]
	save(rps, file= '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trms_assignments.rda')
	#		
	#	save transmission window assignments		
	#		
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219', pattern='_trmwout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_trmwout.rda','',basename(FILE))]					
	rpw		<- infiles[, {
				load(FILE)
				dwin
			}, by=c('DIR','FILE')]
	save(rpw, file= '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments.rda')	
}

RakaiAll.preprocess.pairs.170120<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )	
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision/RCCS_SeqInfo_160816.rda')
	#	get epi info
	tmp		<- RakaiCirc.epi.get.info()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#	add sequence dates to rd
	tmp		<- unique(subset(rs, !is.na(PID), select=c(PID, DATE)),by='PID')
	setnames(tmp, 'DATE','SEQDATE')
	rd		<- merge(rd, tmp, by='PID',all.x=1)
	#	focus on those with PANGEA seqs
	rd		<- subset(rd, !is.na(PID))
	#	focus on clinical data and location data closest to time of diagnosis
	tmp		<- rd[, list(VISIT= VISIT[which.min(DATE-FIRSTPOSDATE)]), by='RID']
	ri		<- subset(merge(unique(rd, by=c('RID','VISIT')), tmp, by=c('RID','VISIT')), select=c(RID, VISIT, DATE, BIRTHDATE, RELIGION, REGION, COMM_NUM, HH_NUM, SEX, LASTNEGDATE, FIRSTPOSDATE, RECENTCD4, RECENTCD4DATE, RECENTVL, RECENTVLDATE, ARVSTARTDATE))
	#	add all PANGEA sequences
	ri		<- merge(ri, unique(subset(rd, select=c(RID, PID, SEQDATE))), by='RID')
	#	focus on behaviour data closest to time of diagnosis
	tmp		<- unique(subset(ri, select=c(RID, VISIT)))
	tmp		<- subset(merge(rh, tmp, by=c('RID','VISIT')), select=c(RID, VISIT, SEXYEAR, EVERSEX, CIRC, OCCUP1, OCCUP2, SEXP1YR, SEXP1OUT, SEXPEVER, SEXC))
	ri		<- merge(ri, tmp, by=c('RID','VISIT'))
	setnames(ri, 'DATE', 'VISIT_DATE')
	#	add community types to rd	
	# 	from Kate:
	#	Trading communities: 1   4   16  22  24  33  51  107 776
	#	Fishing communities: 38, 770, 771, 774		
	tmp		<- data.table(	COMM_NUM=	c("1","2","3","4","5","6","7","8","16","18","19","22","23","24","25","29","33","34","36","38","40","51","54","55", "56","57","58","59","62","74","77","81","89","94","95","103","106","107","108","109","120","177", "370","391","401","451", "602", "754", "760", "770","771","772","773","774","776"),
							COMM_TYPE=	c("T","A","A","T","A","A","A","A", "T", "A", "A", "T", "A", "T", "A", "A", "T", "A", "A", "F", "A", "T", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",  "A","A",   "A",  "T",  "A",  "A",  "A",  "A",   "A",  "A", "A",  "A",    "A",  "A",    "A", "F",  "F",  "A",  "A",  "F",   "T"))
	set(tmp, NULL, 'COMM_TYPE', tmp[, as.character(factor(COMM_TYPE, levels=c('A','T','F'), labels=c('agrarian','trading','fisherfolk')))])
	set(tmp, NULL, 'COMM_NUM', tmp[, as.integer(COMM_NUM)])
	ri		<- merge(ri, tmp, by='COMM_NUM')
	
	pty.runs[, PID:= gsub('-S[0-9]+$','',TAXA)]
	stopifnot( !length(setdiff(pty.runs[, sort(unique(PID))], ri[, sort(unique(PID))])) )
	
	#
	#	for each run: save trm assignments for couples
	#
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219', pattern='_phscout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_phscout.rda','',basename(FILE))]
	infiles	<- subset(infiles, grepl('d50_p001_mr20_mt1_cl2_d5',RUN))
	invisible(infiles[, {
						#FILE	<- '/Users/Oliver/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5_phscout.rda'
						cat('\n',FILE)
						load(FILE)	#loads phs dtrms dtrees dwin
						#
						#	extract couples transmissions summary
						#
						#	check
						tmp		<- dtrms[, list(CH=WIN_TOTAL-sum(WIN_OF_TYPE)), by=c('PTY_RUN','ID1','ID2')]
						stopifnot( tmp[, all(CH==0)] )
						#	keep association between TYPE and TYPE_PAIR
						dtypes	<- unique(subset(dtrms, select=c(TYPE, TYPE_PAIR)))
						set(dtypes, NULL, 'TYPE', dtypes[, gsub('12','mf',TYPE)])						
						set(dtypes, NULL, 'TYPE', dtypes[, gsub('21','fm',TYPE)])						
						#	add zeros
						tmp		<- dcast.data.table(dtrms, PAIR_ID~TYPE, value.var='WIN_OF_TYPE')
						for(x in setdiff(colnames(tmp),'PAIR_ID'))
							set(tmp, which(is.na(tmp[[x]])), x, 0)						
						tmp		<- melt.data.table(tmp, id.vars='PAIR_ID', value.name='WIN_OF_TYPE', variable.name='TYPE')								
						dtrms	<- merge(unique(subset(dtrms, select=c(PAIR_ID, ID1, ID2, ID1_R, ID1_L, ID2_R, ID2_L, PTY_RUN, WIN_TOTAL, SCORE)), by=c('ID1','ID2','PTY_RUN')), tmp, by='PAIR_ID')
						#	double the likely pairs (preserving the inferred direction), 
						#	making it easier to select male / female covariates below
						tmp		<- copy(dtrms)
						setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))						
						set(tmp, NULL, 'TYPE', tmp[, gsub('12','XX33XX',TYPE)])						
						set(tmp, NULL, 'TYPE', tmp[, gsub('21','12',TYPE)])
						set(tmp, NULL, 'TYPE', tmp[, gsub('XX33XX','21',TYPE)])						
						dtrms	<- rbind(tmp, dtrms)
						#	add PANGEAIDs
						tmp		<- unique(subset(pty.runs, select=c(FILE_ID,TAXA)))
						setnames(tmp, c('FILE_ID','TAXA'), c('ID1','MALE_TAXA'))
						dtrms	<- merge(dtrms, tmp, by='ID1')
						setnames(tmp, c('ID1','MALE_TAXA'), c('ID2','FEMALE_TAXA'))
						dtrms	<- merge(dtrms, tmp, by='ID2')
						dtrms[, MALE_PID:= gsub('-S[0-9]+$','',MALE_TAXA)]
						dtrms[, FEMALE_PID:= gsub('-S[0-9]+$','',FEMALE_TAXA)]
						#	merge with rd
						tmp		<- copy(ri)            
						setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))	
						dtrms	<- merge(dtrms, tmp, by='MALE_PID')
						setnames(tmp, colnames(tmp), gsub('MALE','FEMALE',colnames(tmp)))
						dtrms	<- merge(dtrms, tmp, by='FEMALE_PID')						
						#	reduce likely pairs to male-1 and female-2
						dtrms	<- subset(dtrms, MALE_SEX=='M' & FEMALE_SEX=='F')
						set(dtrms, NULL, 'TYPE', dtrms[, gsub('12','mf',TYPE)])						
						set(dtrms, NULL, 'TYPE', dtrms[, gsub('21','fm',TYPE)])
						set(dtrms, NULL, 'TYPE', dtrms[, as.character(TYPE)])
						setnames(dtrms, colnames(dtrms), gsub('ID1','MALE_SANGER_ID',colnames(dtrms)))
						setnames(dtrms, colnames(dtrms), gsub('ID2','FEMALE_SANGER_ID',colnames(dtrms)))
						#	
						dtrms[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]
						set(dtrms, NULL, 'RUN', RUN)
						dtrms	<- merge(dtrms, dtypes, by='TYPE')
						save(dtrms, file=gsub('phscout.rda','allpairs_trmsout.rda',FILE))
						#
						#	extract couples transmissions assignments per window
						#
						#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
						tmp		<- copy(dwin)
						setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
						set(tmp, NULL, 'TYPE', tmp[, gsub('12','XX33XX',TYPE)])						
						set(tmp, NULL, 'TYPE', tmp[, gsub('21','12',TYPE)])
						set(tmp, NULL, 'TYPE', tmp[, gsub('XX33XX','21',TYPE)])						
						dwin	<- rbind(tmp, dwin)
						#	add PANGEAIDs
						tmp		<- unique(subset(pty.runs, select=c(FILE_ID,TAXA)))
						setnames(tmp, c('FILE_ID','TAXA'), c('ID1','MALE_TAXA'))
						dwin	<- merge(dwin, tmp, by='ID1')
						setnames(tmp, c('ID1','MALE_TAXA'), c('ID2','FEMALE_TAXA'))
						dwin	<- merge(dwin, tmp, by='ID2')
						dwin[, MALE_PID:= gsub('-S[0-9]+$','',MALE_TAXA)]
						dwin[, FEMALE_PID:= gsub('-S[0-9]+$','',FEMALE_TAXA)]
						#	merge with rd
						tmp		<- copy(ri)            
						setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))	
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
						save(dwin, file=gsub('phscout.rda','allpairs_trmwout.rda',FILE))				
					}, by=c('RUN','DIR','FILE')])
	#		
	#	save transmission summary assignments		
	#		
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219', pattern='_allpairs_trmsout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_allpairs_trmsout.rda','',basename(FILE))]					
	rps		<- infiles[, {
				load(FILE)
				dtrms
			}, by=c('DIR','FILE')]
	save(rps, file= '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trms_assignments_allpairs.rda')
	#		
	#	save transmission window assignments		
	#		
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219', pattern='_allpairs_trmwout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_allpairs_trmwout.rda','',basename(FILE))]					
	rpw		<- infiles[, {
				load(FILE)
				dwin
			}, by=c('DIR','FILE')]
	save(rpw, file= '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments_allpairs.rda')	
}

RakaiAll.preprocess.pairs.170227<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )	
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision/RCCS_SeqInfo_160816.rda')
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
	#	add community types to rd	
	# 	from Kate:
	#	Trading communities: 1   4   16  22  24  33  51  107 776
	#	Fishing communities: 38, 770, 771, 774		
	tmp		<- data.table(	COMM_NUM=	c("1","2","3","4","5","6","7","8","16","18","19","22","23","24","25","29","33","34","36","38","40","51","54","55", "56","57","58","59","62","74","77","81","89","94","95","103","106","107","108","109","120","177", "370","391","401","451", "602", "754", "760", "770","771","772","773","774","776"),
			COMM_TYPE=	c("T","A","A","T","A","A","A","A", "T", "A", "A", "T", "A", "T", "A", "A", "T", "A", "A", "F", "A", "T", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",  "A","A",   "A",  "T",  "A",  "A",  "A",  "A",   "A",  "A", "A",  "A",    "A",  "A",    "A", "F",  "F",  "A",  "A",  "F",   "T"))
	set(tmp, NULL, 'COMM_TYPE', tmp[, as.character(factor(COMM_TYPE, levels=c('A','T','F'), labels=c('agrarian','trading','fisherfolk')))])
	set(tmp, NULL, 'COMM_NUM', tmp[, as.integer(COMM_NUM)])
	ri		<- merge(ri, tmp, by='COMM_NUM')
	
	pty.runs[, PID:= gsub('-S[0-9]+$','',TAXA)]
	stopifnot( !length(setdiff(pty.runs[, sort(unique(PID))], ri[, sort(unique(PID))])) )	
	#
	#	for each run: save trm assignments for couples
	#
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/170227', pattern='_phscout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_phscout.rda','',basename(FILE))]
	invisible(infiles[, {
						#FILE	<- '/Users/Oliver/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5_phscout.rda'
						cat('\n',FILE)
						load(FILE)	#loads phs dtrms dtrees dwin
						#
						#	extract couples transmissions summary
						#
						#	check
						tmp		<- dtrms[, list(CH=WIN_TOTAL-sum(WIN_OF_TYPE)), by=c('PTY_RUN','ID1','ID2')]
						stopifnot( tmp[, all(CH==0)] )
						#	keep association between TYPE and TYPE_PAIR
						dtypes	<- unique(subset(dtrms, select=c(TYPE, TYPE_PAIR)))
						set(dtypes, NULL, 'TYPE', dtypes[, gsub('12','mf',TYPE)])						
						set(dtypes, NULL, 'TYPE', dtypes[, gsub('21','fm',TYPE)])						
						#	add zeros
						tmp		<- dcast.data.table(dtrms, PAIR_ID~TYPE, value.var='WIN_OF_TYPE')
						for(x in setdiff(colnames(tmp),'PAIR_ID'))
							set(tmp, which(is.na(tmp[[x]])), x, 0)						
						tmp		<- melt.data.table(tmp, id.vars='PAIR_ID', value.name='WIN_OF_TYPE', variable.name='TYPE')								
						dtrms	<- merge(unique(subset(dtrms, select=c(PAIR_ID, ID1, ID2, ID1_R, ID1_L, ID2_R, ID2_L, PTY_RUN, WIN_TOTAL, SCORE)), by=c('ID1','ID2','PTY_RUN')), tmp, by='PAIR_ID')
						#	double the likely pairs (preserving the inferred direction), 
						#	making it easier to select male / female covariates below
						tmp		<- copy(dtrms)
						setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))						
						set(tmp, NULL, 'TYPE', tmp[, gsub('12','XX33XX',TYPE)])						
						set(tmp, NULL, 'TYPE', tmp[, gsub('21','12',TYPE)])
						set(tmp, NULL, 'TYPE', tmp[, gsub('XX33XX','21',TYPE)])						
						dtrms	<- rbind(tmp, dtrms)
						#	add PANGEAIDs
						tmp		<- unique(subset(pty.runs, select=c(FILE_ID,TAXA)))
						setnames(tmp, c('FILE_ID','TAXA'), c('ID1','MALE_TAXA'))
						dtrms	<- merge(dtrms, tmp, by='ID1')
						setnames(tmp, c('ID1','MALE_TAXA'), c('ID2','FEMALE_TAXA'))
						dtrms	<- merge(dtrms, tmp, by='ID2')
						dtrms[, MALE_PID:= gsub('-S[0-9]+$','',MALE_TAXA)]
						dtrms[, FEMALE_PID:= gsub('-S[0-9]+$','',FEMALE_TAXA)]
						#	merge with rd
						tmp		<- copy(ri)            
						setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))	
						dtrms	<- merge(dtrms, tmp, by='MALE_PID')
						setnames(tmp, colnames(tmp), gsub('MALE','FEMALE',colnames(tmp)))
						dtrms	<- merge(dtrms, tmp, by='FEMALE_PID')						
						#	reduce likely pairs to male-1 and female-2
						dtrms	<- subset(dtrms, MALE_SEX=='M' & FEMALE_SEX=='F')
						set(dtrms, NULL, 'TYPE', dtrms[, gsub('12','mf',TYPE)])						
						set(dtrms, NULL, 'TYPE', dtrms[, gsub('21','fm',TYPE)])
						set(dtrms, NULL, 'TYPE', dtrms[, as.character(TYPE)])
						setnames(dtrms, colnames(dtrms), gsub('ID1','MALE_SANGER_ID',colnames(dtrms)))
						setnames(dtrms, colnames(dtrms), gsub('ID2','FEMALE_SANGER_ID',colnames(dtrms)))
						#	
						dtrms[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]
						set(dtrms, NULL, 'RUN', RUN)
						dtrms	<- merge(dtrms, dtypes, by='TYPE')
						save(dtrms, file=gsub('phscout.rda','allpairs_trmsout.rda',FILE))
						#
						#	extract couples transmissions assignments per window
						#
						#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
						tmp		<- copy(dwin)
						setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
						set(tmp, NULL, 'TYPE', tmp[, gsub('12','XX33XX',TYPE)])						
						set(tmp, NULL, 'TYPE', tmp[, gsub('21','12',TYPE)])
						set(tmp, NULL, 'TYPE', tmp[, gsub('XX33XX','21',TYPE)])						
						dwin	<- rbind(tmp, dwin)
						#	add PANGEAIDs
						tmp		<- unique(subset(pty.runs, select=c(FILE_ID,TAXA)))
						setnames(tmp, c('FILE_ID','TAXA'), c('ID1','MALE_TAXA'))
						dwin	<- merge(dwin, tmp, by='ID1')
						setnames(tmp, c('ID1','MALE_TAXA'), c('ID2','FEMALE_TAXA'))
						dwin	<- merge(dwin, tmp, by='ID2')
						dwin[, MALE_PID:= gsub('-S[0-9]+$','',MALE_TAXA)]
						dwin[, FEMALE_PID:= gsub('-S[0-9]+$','',FEMALE_TAXA)]
						#	merge with rd
						tmp		<- copy(ri)            
						setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))	
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
						save(dwin, file=gsub('phscout.rda','allpairs_trmwout.rda',FILE))				
					}, by=c('RUN','DIR','FILE')])
	#		
	#	save transmission window assignments		
	#		
	infiles	<- data.table(FILE=list.files('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/170227', pattern='_allpairs_trmwout.rda', full.names=TRUE))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_allpairs_trmwout.rda','',basename(FILE))]					
	rpw		<- infiles[, {
				load(FILE)
				dwin
			}, by=c('DIR','FILE')]
	save(rpw, file= '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/170227/RCCS_170227_w270_trmw_assignments_allpairs.rda')	 	
}

RakaiCouples.analyze.couples.161107.trmw<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161107/RCCS_161107_w270_trmw_assignments.rda')
	rpw		<- subset(rpw, RUN%in%c("RCCS_161107_w270_d200_r004_mr1", "RCCS_161107_w270_d50_r004_mr1", "RCCS_161107_w270_d20_r004_mr1", "RCCS_161107_w270_d50_r004_mr20", "RCCS_161107_w270_d50_r004_mr30", "RCCS_161107_w270_d50_r004_mr40", "RCCS_161107_w270_d50_r004_mr50", "RCCS_161107_w270_d50_r004_mr100") )
	#
	set(rpw, NULL, 'RUN', rpw[, factor(RUN, levels=c( 	"RCCS_161107_w270_d20_r004_mr1","RCCS_161107_w270_d50_r004_mr1","RCCS_161107_w270_d200_r004_mr1",
														"RCCS_161107_w270_d50_r004_mr20","RCCS_161107_w270_d50_r004_mr30","RCCS_161107_w270_d50_r004_mr40","RCCS_161107_w270_d50_r004_mr50","RCCS_161107_w270_d50_r004_mr100"))])
	rpw[, table(RUN, useNA='if')]
	#
	#	for each run: plot pairs	
	#	
	#	define plotting order: largest number of trm assignments
	tmp		<- rpw[, list( WIN_TR=length(which(grepl('trans',TYPE))) ), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	tmp		<- tmp[order(-WIN_TR),]
	tmp[, PLOT_ID:=seq_len(nrow(tmp))]	
	#	define label
	tmp		<- merge(tmp, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair ', COUPID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- merge(subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL)), rpw, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('trans_mf','trans_fm','unint','int','cher','disconnected'), labels=c('M transmit to F','F transmit to M','M, F are unint','M, F are intermingled','M, F are a cherry','M, F are disconnected'))])
	rpwp	<- subset(tmp, select=c(RUN, PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, LABEL, W_FROM,  W_TO, TYPE, ID1_R, ID1_L, ID2_R, ID2_L ))
	tmp		<- rpwp[, list(	ID_R_MIN=min(ID1_R, ID2_R),
							ID_R_MAX=max(ID1_R, ID2_R)), by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM')]
	rpwp	<- merge(rpwp, tmp, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM'))
	set(rpwp, NULL, 'RUN', rpwp[, gsub('_','\n',gsub('w270_','',gsub('RCCS_','',as.character(RUN))))])
	set(rpwp, NULL, 'RUN', rpwp[, factor(RUN, levels=c(	"161107\nd200\nr004\nmr1", "161107\nd50\nr004\nmr1", "161107\nd20\nr004\nmr1",   
														"161107\nd50\nr004\nmr20","161107\nd50\nr004\nmr30","161107\nd50\nr004\nmr40","161107\nd50\nr004\nmr50","161107\nd50\nr004\nmr100"))])
	
	run		<- 'RCCS_161107_w270_dxxx'
	dir		<- rpw$DIR[1]
	#	F->M
	tmp		<- subset(rpwp, COUP_SC=='F->M')
	p		<- lapply(tmp[, unique(LABEL)], function(label)
				{
					p	<- ggplot(subset(tmp, LABEL==label), aes(x=W_FROM)) +			
							geom_bar(aes(y=ID_R_MAX, colour=TYPE), stat='identity', fill='transparent') +
							geom_bar(aes(y=ID_R_MIN, fill=TYPE), stat='identity', colour='transparent') +
							labs(x='window start', y='number of reads', fill='topology of clades\nbetween patient pairs', colour='topology of clades\nbetween patient pairs') +
							scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +
							scale_colour_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +
							scale_x_continuous(breaks=seq(0,1e4,250), limits=c(tmp[, min(W_FROM)], tmp[, max(W_FROM)])) +
							scale_y_log10(breaks=c(10,100,1000,1e4,1e5)) +
							theme_bw() + theme(legend.position='top') +
							facet_grid(RUN~LABEL) +
							guides(fill=guide_legend(ncol=2)) 
					p
				})
	pdf(file=file.path(dir, paste(run,'-phsc-windowassignments_F2M.pdf',sep='')), w=25, h=1*tmp[, length(unique(RUN))])
	for(i in seq_along(p))	print(p[[i]])
	dev.off()
	#	M->F
	tmp		<- subset(rpwp, COUP_SC=='M->F')
	p		<- lapply(tmp[, unique(LABEL)], function(label)
			{
				p	<- ggplot(subset(tmp, LABEL==label), aes(x=W_FROM)) +			
						geom_bar(aes(y=ID_R_MAX, colour=TYPE), stat='identity', fill='transparent') +
						geom_bar(aes(y=ID_R_MIN, fill=TYPE), stat='identity', colour='transparent') +
						labs(x='window start', y='number of reads', fill='topology of clades\nbetween patient pairs', colour='topology of clades\nbetween patient pairs') +
						scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +
						scale_colour_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +
						scale_x_continuous(breaks=seq(0,1e4,250), limits=c(tmp[, min(W_FROM)], tmp[, max(W_FROM)])) +
						scale_y_log10(breaks=c(10,100,1000,1e4,1e5)) +
						theme_bw() + theme(legend.position='top') +
						facet_grid(RUN~LABEL) +
						guides(fill=guide_legend(ncol=2)) 
				p
			})
	pdf(file=file.path(dir, paste(run,'-phsc-windowassignments_M2F.pdf',sep='')), w=25, h=1*tmp[, length(unique(RUN))])
	for(i in seq_along(p))	print(p[[i]])
	dev.off()
	
}

RakaiCouples.analyze.couples.161110.bbm.model<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	#
	#	plot how the variance changes with n and pi
	#
	
	#	alpha, beta --> Inf
	db	<- data.table(expand.grid(ptrue=seq(0.1, 0.9, 0.2), k=1, m=c(0.1, 1, 10, 100, 1000)))
	db[, n:=1]
	
	#	plot out posterior
	tmp	<- db[, list(	x=seq(0.01,0.99,0.01), 
					y=dbeta(seq(0.01,0.99,0.01), k+m*ptrue, n-k+m*(1-ptrue))), by=c('ptrue','m')]
	ggplot(tmp, aes(x=x, y=y)) + geom_line() + facet_grid(m~ptrue, scales='free_y')		
}
	

RakaiCouples.consensus.to.couples.161205.model<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)	
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	#	load PANGEA alignment to create raw genetic distances
	load('~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda')	
		
	#
	#	raw genetic distance
	#
	
	#	do we have replicates in the couples? yes ...
	tmp	<- subset(dm, grepl('-R2',TAXA) & !is.na(STUDY_ID))[, as.character(STUDY_ID)]
	intersect(tmp,rp[, FEMALE_RID])
	intersect(tmp,rp[, MALE_RID])
	subset(rp, FEMALE_RID=="B115996")	# and these are not tracked in 'rp'
	#	need to update rp
	tmp	<- subset(dm, grepl('-R2',TAXA) & !is.na(STUDY_ID), select=c(STUDY_ID, TAXA, SANGER_ID))
	setnames(tmp, c('STUDY_ID','TAXA','SANGER_ID'), c('MALE_RID','MALE_TAXA_NEW','MALE_SANGER_ID'))
	rp	<- merge(rp, tmp, by=c('MALE_RID','MALE_SANGER_ID'), all.x=1)
	setnames(tmp, c('MALE_RID','MALE_TAXA_NEW','MALE_SANGER_ID'), c('FEMALE_RID','FEMALE_TAXA_NEW','FEMALE_SANGER_ID'))
	rp	<- merge(rp, tmp, by=c('FEMALE_RID','FEMALE_SANGER_ID'), all.x=1)
	tmp	<- rp[, which(!is.na(MALE_TAXA_NEW))]
	set(rp, tmp, 'MALE_TAXA', rp[tmp, MALE_TAXA_NEW])
	tmp	<- rp[, which(!is.na(FEMALE_TAXA_NEW))]
	set(rp, tmp, 'FEMALE_TAXA', rp[tmp, FEMALE_TAXA_NEW])
	set(rp, NULL, c('FEMALE_TAXA_NEW','MALE_TAXA_NEW'), NULL)
	#
	tmp		<- unique(data.table(TAXA= c(rp[, MALE_TAXA], rp[, FEMALE_TAXA])))
	cat('no PANGEA sequences for',setdiff( tmp[,TAXA] , rownames(sq) ))	
	rp.sq	<- sq[intersect( tmp[,TAXA] , rownames(sq) ),]
	rp.cd	<- as.matrix(dist.dna(rp.sq, model='raw', pairwise.deletion=TRUE))	
	rp.cd	<- as.data.table(melt(rp.cd, value.name="CONS_GD"))
	setnames(rp.cd, c('Var1','Var2'), c('MALE_TAXA','FEMALE_TAXA'))
	rp		<- merge(rp, rp.cd, by=c('MALE_TAXA','FEMALE_TAXA'))
	
	#
	#	patristic distance in gag fasttree
	#		

	#	load 'dgd', 'sqi', 'sq'
	load('~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/explaingaps/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda')
	tmp			<- dcast.data.table(subset(dgd, GENE%in%c('GAG','POL')), TAXA~GENE, value.var='ACTG')	
	phsc.input	<- unique(subset(pty.runs, select=c('TAXA','FILE_ID')))
	phsc.input	<- merge(phsc.input, tmp, by='TAXA')
	phsc.input[, GAGPOL_ACTG:= GAG*1801+POL*(4650-1762+1)]
	#	keep taxa that rouhghly meet the phsc requirement of a contiguous 270 bp read
	phsc.input	<- subset(phsc.input, GAGPOL_ACTG>270)		
	#	select gag+pol for these
	tmp			<- sq[ c(phsc.input[, TAXA], subset(sqi, PNG=='N')[, TAXA]), 1:4650]
	#	replace ? with -
	tmp			<- as.character(tmp)
	tmp[tmp=='?']	<- '-'
	tmp			<- as.DNAbin(tmp)
	file.ftin	<- "/Users/Oliver/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_gagpol.fasta"
	file.ftout	<- "/Users/Oliver/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_gagpol_fasttree.newick"
	write.dna(tmp, file.ftin, format='fa', colsep='', nbcol=-1)
	#	get FastTree and ExaML tree
	require(big.phylo)
	tmp			<- cmd.fasttree(file.ftin, file.ftout)
	cat(tmp)	
	file.ftin	<- "/Users/Oliver/duke/tmp/Couples_PANGEA_HIV_n4562_Imperial_v151113_gagpol.fasta"
	tmp2		<- cmd.examl.bootstrap(dirname(file.ftin), basename(file.ftin), bs.from=0, bs.to=0, outdir=dirname(file.ftin))	
	cat(tmp2)
	#	read trees
	phf			<- read.tree(file.ftout)
	file.ftout	<- "/Users/Oliver/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_gagpol_examl.newick"
	phe			<- read.tree(file.ftout)
	#	get patristic distances
	phd			<- as.data.table(melt(as.matrix(cophenetic.phylo(phf)), varnames=c('TAXA1','TAXA2'), value.name='PD_FT'))
	set(phd, NULL, 'TAXA1', phd[, as.character(TAXA1)])
	set(phd, NULL, 'TAXA2', phd[, as.character(TAXA2)])
	tmp			<- as.data.table(melt(as.matrix(cophenetic.phylo(phe)), varnames=c('TAXA1','TAXA2'), value.name='PD_EX'))
	set(tmp, NULL, 'TAXA1', tmp[, as.character(TAXA1)])
	set(tmp, NULL, 'TAXA2', tmp[, as.character(TAXA2)])
	phd			<- merge(phd, tmp, by=c('TAXA1','TAXA2'))
	setnames(phd, c('TAXA1','TAXA2'), c('MALE_TAXA','FEMALE_TAXA'))	
	rp			<- merge(rp, phd, by=c('MALE_TAXA','FEMALE_TAXA'), all.x=1)
	setnames(rp, c('PD_FT','PD_EX'), c('CONS_PD_FT','CONS_PD_EX'))
	save(rp, file="~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
}



RakaiCouples.analyze.couples.161110.bb.model<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)	
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161110/RCCS_161110_w270_trmw_assignments.rda')
	rpw		<- subset(rpw, RUN%in%c("RCCS_161110_w270_d200_r004_mr1_mt1", "RCCS_161110_w270_d50_r004_mr1_mt1", "RCCS_161110_w270_d20_r004_mr1_mt1", "RCCS_161110_w270_d50_r004_mr20_mt2", "RCCS_161110_w270_d50_r004_mr30_mt2", "RCCS_161110_w270_d50_r004_mr40_mt2", "RCCS_161110_w270_d50_r004_mr50_mt2", "RCCS_161110_w270_d50_r004_mr100_mt2") )
	#
	set(rpw, NULL, 'RUN', rpw[, factor(RUN, levels=c( 	"RCCS_161110_w270_d20_r004_mr1_mt1","RCCS_161110_w270_d50_r004_mr1_mt1","RCCS_161110_w270_d200_r004_mr1_mt1",
									"RCCS_161110_w270_d50_r004_mr20_mt2","RCCS_161110_w270_d50_r004_mr30_mt2","RCCS_161110_w270_d50_r004_mr40_mt2","RCCS_161110_w270_d50_r004_mr50_mt2","RCCS_161110_w270_d50_r004_mr100_mt2"))])
	rpw[, table(RUN, useNA='if')]
	rpw[, table(TYPE, useNA='if')]
	rpw		<- subset(rpw, COUP_SC=='seroinc')
	#
	#	define TPAIR_TYPE
	#
	set(rpw, NULL, 'TPAIR_TYPE', NA_character_) 
	set(rpw, rpw[, which(TYPE%in%c('trans_fm', 'trans_mf'))], 'TPAIR_TYPE', 'linked')		
	set(rpw, rpw[, which(TYPE%in%c('disconnected', 'unint'))], 'TPAIR_TYPE', 'unlinked')
	set(rpw, rpw[, which(TYPE%in%c('cher', 'int'))], 'TPAIR_TYPE', 'ambiguous')
	stopifnot( rpw[, !any(is.na(TPAIR_TYPE))])
	
	#	define plotting order: largest number of trm assignments
	tmp		<- rpw[, list( WIN_TR=length(which(grepl('trans',TYPE))) ), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	tmp		<- tmp[order(-WIN_TR),]
	tmp[, PLOT_ID:=seq_len(nrow(tmp))]	
	#	define legend for couples
	tmp		<- merge(tmp, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
	setkey(tmp, PLOT_ID)
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair ', COUPID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp[, LABEL_SH:=tmp[, factor(PLOT_ID, levels=PLOT_ID, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))]]	
	rpleg	<- subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL, LABEL_SH, PLOT_ID))
	setkey(rpleg, LABEL)
	#	define effectively independent number of windows
	#	select only W_FROM that are > W_TO	
	tmp		<- rpw[, {
				z		<- 1
				repeat
				{
					zz		<- which(W_FROM>max(W_TO[z]))[1]
					if(length(zz)==0 | is.na(zz))	break
					z		<- c(z, zz)			
				}
				list(W_FROM=W_FROM[z], W_TO=W_TO[z])
			}, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN')]
	tmp[, OVERLAP:= 0L]
	rpw		<- merge(rpw, tmp, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','W_FROM','W_TO','RUN'), all.x=1 )
	set(rpw, rpw[, which(is.na(OVERLAP))], 'OVERLAP', 1L)	
	#	for each pair count windows by type (k) 
	rplkl	<- rpw[, list(K=length(W_FROM)), by=c('TPAIR_TYPE','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN')]
	rplkl	<- dcast.data.table(rplkl, RUN+PTY_RUN+FEMALE_SANGER_ID+MALE_SANGER_ID~TPAIR_TYPE, value.var='K')
	for(x in setdiff(colnames(rplkl),c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID')))
		set(rplkl, which(is.na(rplkl[[x]])), x, 0L)
	rplkl	<- melt(rplkl, id.vars=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID'), variable.name='TPAIR_TYPE', value.name='K')
	#	for each pair count total windows (n) and effective windows (neff)	
	tmp		<- rpw[, list(N=length(W_FROM), NEFF=sum(1-OVERLAP)), by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN')]
	rplkl	<- merge(rplkl, tmp, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN'))
	#	define effective number of windows of type
	rplkl[, KEFF:= K/N*NEFF]	
	#	add prior (equal probabilities to all types)
	rplkl[, DIR_PRIOR:= 1L]
	#	add posterior
	rplkl[, DIR_PO:= DIR_PRIOR+KEFF]
	#	add marginal posterior mean and 95%
	tmp		<- rplkl[, {
				alpha	<- DIR_PO
				beta	<- sum(DIR_PO)-DIR_PO				
				list(	TPAIR_TYPE=TPAIR_TYPE, BE_ALPHA=alpha, BE_BETA=beta, BE_MEAN=alpha/(alpha+beta), 
						BE_QL=qbeta(0.025, alpha, beta), BE_QU=qbeta(0.975, alpha, beta),
						BE_IL=qbeta(0.25, alpha, beta), BE_IU=qbeta(0.75, alpha, beta))	
			}, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN')]
	rplkl	<- merge(rplkl, tmp, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','TPAIR_TYPE','RUN'))
	#
	#	plot posterior densities
	#	
	run		<- 'RCCS_161110_w270_dxxx'
	dir		<- rpw$DIR[1]

	rpwp	<- merge(rplkl, rpleg, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID'))
	set(rpwp, NULL, 'RUN', rpwp[, gsub('_','\n',gsub('w270_','',gsub('RCCS_','',as.character(RUN))))])
	set(rpwp, NULL, 'RUN', rpwp[, factor(RUN, levels=c(	"161110\nd200\nr004\nmr1\nmt1", "161110\nd50\nr004\nmr1\nmt1", "161110\nd20\nr004\nmr1\nmt1",   
														"161110\nd50\nr004\nmr20\nmt2","161110\nd50\nr004\nmr30\nmt2","161110\nd50\nr004\nmr40\nmt2","161110\nd50\nr004\nmr50\nmt2","161110\nd50\nr004\nmr100\nmt2"))])
	set(rpwp, NULL, 'TPAIR_TYPE', rplkl[, factor(as.character(TPAIR_TYPE), levels=c('linked','unlinked','ambiguous'))])
	rpwp	<- merge(rpwp, subset(rp, select=c(FEMALE_SANGER_ID, MALE_SANGER_ID, CONS_GD, CONS_PD_FTUGGAG)), by=c('FEMALE_SANGER_ID','MALE_SANGER_ID'))
	#
	#	long plot
	#
	p		<- lapply(rpwp[, unique(LABEL)], function(label)
			{
				p	<- ggplot(subset(rpwp, LABEL==label), aes(x=TPAIR_TYPE, fill=TPAIR_TYPE)) + facet_grid(RUN~LABEL) +
						geom_boxplot(aes(lower=BE_IL, middle=BE_MEAN, upper=BE_IU, ymin=BE_QL, ymax=BE_QU), stat="identity", show.legend=FALSE) +
						scale_y_continuous(labels=percent, expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.1)) +
						#scale_fill_manual(values=c('linked'="#9E0142",'ambiguous'="#3288BD",'unlinked'='grey50')) +
						scale_fill_manual(values=c('linked'="#F46D43",'ambiguous'="#3288BD",'unlinked'='grey50')) +			
						theme_bw() + coord_flip() +
						labs(x='phyloscan classification (v16-11-16)\n(transmission pair)\n', y='\nposterior probability')
				p
			})
	pdf(file=file.path(dir, paste(run,'-phsc-windowassignments_SEROINC_trmpairposterior.pdf',sep='')), w=10, h=0.8*rpwp[, length(unique(RUN))])
	for(i in seq_along(p))	print(p[[i]])
	dev.off()
	#
	#	small plot
	#
	tmp		<- subset(rpwp, RUN%in%c('161110\nd200\nr004\nmr1\nmt1','161110\nd50\nr004\nmr1\nmt1','161110\nd20\nr004\nmr1\nmt1','161110\nd50\nr004\nmr20\nmt2'))
	tmp		<- subset(rpwp, RUN%in%c('161110\nd50\nr004\nmr20\nmt2'))
	ggplot(tmp, aes(x=TPAIR_TYPE, fill=TPAIR_TYPE)) + 
			facet_grid(LABEL_SH~RUN) +
			geom_boxplot(aes(lower=BE_IL, middle=BE_MEAN, upper=BE_IU, ymin=BE_QL, ymax=BE_QU), stat="identity") +
			scale_y_continuous(labels=percent, expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.1)) +
			#scale_fill_manual(values=c('linked'="#9E0142",'ambiguous'="#3288BD",'unlinked'='grey50')) +
			scale_fill_manual(values=c('linked'="#F46D43",'ambiguous'="#3288BD",'unlinked'='grey50')) +			
			theme_bw() + coord_flip() + theme(legend.position='top', strip.text.y = element_text(angle=0)) +			
			labs(x='sequence pairs of couples\n', fill='phyloscanner classification v16-11-16\n(transmission pair)\n', y='\nposterior probability')
	ggsave(file=file.path(dir, paste(run,'-phsc-windowassignments_SEROINC_trmpairposterior_small.pdf',sep='')), w=10, h=30, limitsize = FALSE)
	#
	#	for each run: plot pairs	
	#	
	#	define plotting order: largest number of trm assignments
	tmp		<- rpw[, list( WIN_TR=length(which(grepl('trans',TYPE))) ), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	tmp		<- tmp[order(-WIN_TR),]
	tmp[, PLOT_ID:=seq_len(nrow(tmp))]	
	#	define label
	tmp		<- merge(tmp, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair ', COUPID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- merge(subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL)), rpw, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('trans_mf','trans_fm','unint','int','cher','disconnected'), labels=c('M transmit to F','F transmit to M','M, F are unint','M, F are intermingled','M, F are a cherry','M, F are disconnected'))])
	rpwp	<- subset(tmp, select=c(RUN, PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, LABEL, W_FROM,  W_TO, TYPE, ID1_R, ID1_L, ID2_R, ID2_L ))
	tmp		<- rpwp[, list(	ID_R_MIN=min(ID1_R, ID2_R),
					ID_R_MAX=max(ID1_R, ID2_R)), by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM')]
	rpwp	<- merge(rpwp, tmp, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM'))
	set(rpwp, NULL, 'RUN', rpwp[, gsub('_','\n',gsub('w270_','',gsub('RCCS_','',as.character(RUN))))])
	set(rpwp, NULL, 'RUN', rpwp[, factor(RUN, levels=c(	"161110\nd200\nr004\nmr1\nmt1", "161110\nd50\nr004\nmr1\nmt1", "161110\nd20\nr004\nmr1\nmt1",   
									"161110\nd50\nr004\nmr20\nmt2","161110\nd50\nr004\nmr30\nmt2","161110\nd50\nr004\nmr40\nmt2","161110\nd50\nr004\nmr50\nmt2","161110\nd50\nr004\nmr100\nmt2"))])
	#
	#	compare to consensus approach
	#		
	tmp		<- subset(rpwp, RUN%in%c('161110\nd50\nr004\nmr20\nmt2'))
	ggplot(tmp, aes(x=CONS_GD, colour=TPAIR_TYPE)) +
			geom_point(aes(y=BE_MEAN)) +
			geom_errorbar(aes(ymin=BE_IL, ymax=BE_IU)) +
			#geom_boxplot(aes(lower=BE_IL, middle=BE_MEAN, upper=BE_IU, ymin=BE_QL, ymax=BE_QU), stat="identity") +
			scale_y_continuous(labels=percent, expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.1)) +			
			scale_colour_manual(values=c('linked'="#F46D43",'ambiguous'="#3288BD",'unlinked'='grey50')) +			
			facet_grid(~TPAIR_TYPE) +
			theme_bw() + theme(legend.position='bottom') +
			labs(x='raw genetic distance between sequences\nof seroincident pairs', colour='phyloscanner classification v16-11-16\n(transmission pair)\n', y='\nposterior probability')
	ggsave(file=file.path(dir, paste(run,'-phsc-windowassignments_SEROINC_trmpairposterior_small_vs_rawgenetic.pdf',sep='')), w=10, h=6, limitsize = FALSE)
	tmp		<- subset(rpwp, RUN%in%c('161110\nd50\nr004\nmr20\nmt2'))
	ggplot(tmp, aes(x=CONS_PD_FTUGGAG, colour=TPAIR_TYPE)) +
			geom_point(aes(y=BE_MEAN)) +
			geom_errorbar(aes(ymin=BE_IL, ymax=BE_IU)) +
			#geom_boxplot(aes(lower=BE_IL, middle=BE_MEAN, upper=BE_IU, ymin=BE_QL, ymax=BE_QU), stat="identity") +
			scale_y_continuous(labels=percent, expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.1)) +			
			scale_colour_manual(values=c('linked'="#F46D43",'ambiguous'="#3288BD",'unlinked'='grey50')) +			
			facet_grid(~TPAIR_TYPE) +
			theme_bw() + theme(legend.position='bottom') +
			labs(x='patristic distance between sequences\nof seroincident pairs from FastTree on gag', colour='phyloscanner classification v16-11-16\n(transmission pair)\n', y='\nposterior probability')
	ggsave(file=file.path(dir, paste(run,'-phsc-windowassignments_SEROINC_trmpairposterior_small_vs_patrFTgag.pdf',sep='')), w=10, h=6, limitsize = FALSE)
	#
	#	look at pairs with <2% patristic distance and <50% trm prob
	#
	tmp	<- subset(rpwp, CONS_GD<0.02 & TPAIR_TYPE=='linked' & BE_MEAN<.5 & RUN%in%c('161110\nd50\nr004\nmr20\nmt2')) 
}


RakaiCouples.analyze.couples.161110.trmw<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161110/RCCS_161110_w270_trmw_assignments.rda')
	rpw		<- subset(rpw, RUN%in%c("RCCS_161110_w270_d200_r004_mr1_mt1", "RCCS_161110_w270_d50_r004_mr1_mt1", "RCCS_161110_w270_d20_r004_mr1_mt1", "RCCS_161110_w270_d50_r004_mr20_mt2", "RCCS_161110_w270_d50_r004_mr30_mt2", "RCCS_161110_w270_d50_r004_mr40_mt2", "RCCS_161110_w270_d50_r004_mr50_mt2", "RCCS_161110_w270_d50_r004_mr100_mt2") )
	#
	set(rpw, NULL, 'RUN', rpw[, factor(RUN, levels=c( 	"RCCS_161110_w270_d20_r004_mr1_mt1","RCCS_161110_w270_d50_r004_mr1_mt1","RCCS_161110_w270_d200_r004_mr1_mt1",
														"RCCS_161110_w270_d50_r004_mr20_mt2","RCCS_161110_w270_d50_r004_mr30_mt2","RCCS_161110_w270_d50_r004_mr40_mt2","RCCS_161110_w270_d50_r004_mr50_mt2","RCCS_161110_w270_d50_r004_mr100_mt2"))])
	rpw[, table(RUN, useNA='if')]
	#
	#	for each run: plot pairs	
	#	
	#	define plotting order: largest number of trm assignments
	tmp		<- rpw[, list( WIN_TR=length(which(grepl('trans',TYPE))) ), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	tmp		<- tmp[order(-WIN_TR),]
	tmp[, PLOT_ID:=seq_len(nrow(tmp))]	
	#	define label
	tmp		<- merge(tmp, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair ', COUPID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- merge(subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL)), rpw, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('trans_mf','trans_fm','unint','int','cher','disconnected'), labels=c('M transmit to F','F transmit to M','M, F are unint','M, F are intermingled','M, F are a cherry','M, F are disconnected'))])
	rpwp	<- subset(tmp, select=c(RUN, PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, LABEL, W_FROM,  W_TO, TYPE, ID1_R, ID1_L, ID2_R, ID2_L ))
	tmp		<- rpwp[, list(	ID_R_MIN=min(ID1_R, ID2_R),
					ID_R_MAX=max(ID1_R, ID2_R)), by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM')]
	rpwp	<- merge(rpwp, tmp, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM'))
	set(rpwp, NULL, 'RUN', rpwp[, gsub('_','\n',gsub('w270_','',gsub('RCCS_','',as.character(RUN))))])
	set(rpwp, NULL, 'RUN', rpwp[, factor(RUN, levels=c(	"161110\nd200\nr004\nmr1\nmt1", "161110\nd50\nr004\nmr1\nmt1", "161110\nd20\nr004\nmr1\nmt1",   
									"161110\nd50\nr004\nmr20\nmt2","161110\nd50\nr004\nmr30\nmt2","161110\nd50\nr004\nmr40\nmt2","161110\nd50\nr004\nmr50\nmt2","161110\nd50\nr004\nmr100\nmt2"))])
	
	run		<- 'RCCS_161110_w270_dxxx'
	dir		<- rpw$DIR[1]
	#	SEROINC
	tmp		<- subset(rpwp, COUP_SC=='seroinc')
	p		<- lapply(tmp[, unique(LABEL)], function(label)
			{
				p	<- ggplot(subset(tmp, LABEL==label), aes(x=W_FROM)) +			
						geom_bar(aes(y=ID_R_MAX, colour=TYPE), stat='identity', fill='transparent') +
						geom_bar(aes(y=ID_R_MIN, fill=TYPE), stat='identity', colour='transparent') +
						labs(x='window start', y='number of reads', fill='topology of clades\nbetween patient pairs', colour='topology of clades\nbetween patient pairs') +
						scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +
						scale_colour_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +
						scale_x_continuous(breaks=seq(0,1e4,250), limits=c(tmp[, min(W_FROM)], tmp[, max(W_FROM)])) +
						scale_y_log10(breaks=c(10,100,1000,1e4,1e5)) +
						theme_bw() + theme(legend.position='top') +
						facet_grid(RUN~LABEL) +
						guides(fill=guide_legend(ncol=2)) 
				p
			})
	pdf(file=file.path(dir, paste(run,'-phsc-windowassignments_SEROINC.pdf',sep='')), w=25, h=1*tmp[, length(unique(RUN))])
	for(i in seq_along(p))	print(p[[i]])
	dev.off()
	#
	#	get distribution of cherries on manually identified unlinked and linked
	#
	tmp	<- as.data.table(matrix( c('83','15715_1_17','15715_1_18',
							'36','15715_1_17','15715_1_18',
							'85','15699_1_38','15699_1_40',
							'104','15699_1_3','15699_1_4',
							'106','15736_1_28','15102_1_16',
							'121','16059_1_66','16059_1_52',
							'31','15981_1_31','15981_1_33',
							'98','15916_1_7','15916_1_21',
							'64','15916_1_7','15916_1_21',
							'122','15910_1_76','15910_1_82',
							'24','15981_1_17','15835_1_32',
							'70','15833_1_84','15833_1_83',
							'104','16059_1_66','16059_1_52',
							'24','16016_1_11','16017_1_44',
							'24','15915_1_51','15776_1_32',
							'15','16033_1_76','16033_1_60',
							'7','16033_1_76','16033_1_60',
							'18','16019_1_45','16018_1_44',
							'67','15892_1_3','15892_1_1',
							'91','16021_1_66','16021_1_58',
							'29','16021_1_66','16021_1_58',
							'116','15915_1_84','15915_1_83',
							'64','15915_1_84','15915_1_83',
							'50','15115_1_4','15115_1_9',
							'86','15978_1_12','15777_1_42',
							'68','15978_1_12','15777_1_42',
							'44','16056_1_83','15099_1_59'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp[, CLASS_OR:= 'OR_disconnected']	
	tmp2	<- as.data.table(matrix( c('111','15950_1_27','15950_1_25',
							'122','16056_1_58','16056_1_57',
							'38','15105_1_35','16016_1_4',
							'113','15950_1_5','15958_1_47',				
							'99','15915_1_61','15915_1_74',
							'117','15915_1_61','15915_1_74',
							'60','16056_1_83','16056_1_85',
							'39','16034_1_86','16034_1_82',
							'120','15981_1_19','15981_1_23',				
							'79','15978_1_36','15978_1_41'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp2[, CLASS_OR:= 'OR_ambiguous']
	tmp		<- rbind(tmp, tmp2)
	tmp2	<- as.data.table(matrix( c('101','16003_1_76','16003_1_79',
							'37','15114_1_64','15758_1_73',
							'37','16016_1_5','16016_1_4',
							'20','16057_1_15','16057_1_16',
							'59','15699_1_5','15867_1_21',
							'99','15915_1_78','15915_1_64',
							'109','15777_1_29','15777_1_21',
							'63','16003_1_86','16003_1_85',
							'109','15776_1_19','15776_1_28',
							'27','16219_1_78','16219_1_80',
							'120','16059_1_73','16059_1_53',
							'25','15675_1_87','15675_1_77',
							'70','15835_1_21','15835_1_22'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp2[, CLASS_OR:= 'OR_transmission']
	tmp		<- rbind(tmp, tmp2)
	set(tmp, NULL, 'PTY_RUN', tmp[, as.integer(PTY_RUN)])
	rpi		<- copy(tmp)
	#
	#	try genetic distance
	#
	rpd		<- merge(rpi, subset(rpw, TYPE=='cher' & RUN=='RCCS_161110_w270_d50_r004_mr20_mt2'), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rpd		<- subset(rpd, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, W_FROM, W_TO, PATRISTIC_DIST, COUPID, CLASS_OR))
	tmp		<- unique(subset(rpd, select=COUPID))
	setkey(tmp, COUPID)
	tmp[, COUPID_RNK:= seq_len(nrow(tmp))]
	rpd		<- merge(rpd, tmp, by='COUPID')
	rpd[, LEGEND:= paste(COUPID_RNK, ' ',COUPID, '(M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID,')',sep='')]
	ggplot(rpd, aes(x=W_FROM, y=PATRISTIC_DIST, colour=LEGEND)) + 
			geom_point() + 
			geom_hline(aes(yintercept=0.01)) +
			geom_hline(aes(yintercept=0.02)) +
			scale_y_log10(breaks=c(0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2)) +
			facet_grid(~CLASS_OR) + theme_bw()
	ggsave(file=file.path(dir, paste(run,'-phsc-windowassignments_SEROINC_cherries_patristic.pdf',sep='')), w=10, h=5)
	#
	#	try Weibull..
	#	
	tmp		<- subset(rpd, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, W_FROM, W_TO))
	tmp		<- melt(tmp, measure.vars=c('MALE_SANGER_ID','FEMALE_SANGER_ID'), value.name='ID')
	load('~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161110/RCCS_161110_w270_d50_r004_mr20_mt2_phscout.rda')	
	require(gamlss)
	dr		<- tmp[, {
				pty_run	<- PTY_RUN
				w_from	<- W_FROM
				id		<- ID
				z		<- subset(dtrees, PTY_RUN==pty_run & W_FROM==w_from)[, IDX]
				cat('\nRun',pty_run,'Window',w_from,'ID',id,'tree index',z)
				ph		<- phs[[ z ]]				
				#	collect branch lengths
				z	<- ph$tip.label[!grepl(id, ph$tip.label)]
				z	<- drop.tip(ph,z,root.edge=0)	
				brls<- z$edge.length
				#	get length of root edge
				z	<- getMRCA(ph, which(grepl(id, ph$tip.label)))
				rl	<- ph$edge.length[ which(ph$edge[,2]==z) ]
				#
				p	<- NA_real_		# this will be the tail area probability Pr(Wei>rl)
				brls<- brls[ brls>1e-5 ]		
				cat('\nNumber branch lengths >1e-5',length(brls))
				if(length(brls)>3)	# don't think makes sense to fit Weibull to 3 data points
				{	
					#cat('\n', W_FROM, ID) # for debugging		
					#z		<- data.frame(BRL=brls)
					#print(z)
					w		<- gamlss(formula=brls~1, family=WEI, trace=FALSE)
					w.l		<- exp(coef(w, what='mu'))
					w.k		<- exp(coef(w, what='sigma'))
					p		<- pweibull(rl, w.k, scale=w.l, lower.tail=FALSE)					
				}				
				list(ROOTL=rl, ROOTL_P=p) 
			}, by=c('PTY_RUN','W_FROM','W_TO','variable','ID')]
	dr[, variable:=NULL]
	setnames(dr, c('ID','ROOTL','ROOTL_P'), c('MALE_SANGER_ID','MALE_ROOTL','MALE_ROOTL_P'))	
	rpd		<- merge(rpd, dr, by=c('PTY_RUN','W_FROM','W_TO','MALE_SANGER_ID'))
	setnames(dr, c('MALE_SANGER_ID','MALE_ROOTL','MALE_ROOTL_P'), c('FEMALE_SANGER_ID','FEMALE_ROOTL','FEMALE_ROOTL_P'))
	rpd		<- merge(rpd, dr, by=c('PTY_RUN','W_FROM','W_TO','FEMALE_SANGER_ID'))
	tmp		<- rpd[, list(CHERRY_P= max(MALE_ROOTL_P, FEMALE_ROOTL_P)), by=c('PTY_RUN','W_FROM','W_TO','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	rpd		<- merge(rpd, tmp, by=c('PTY_RUN','W_FROM','W_TO','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	#
	#
	#
	tmp		<- subset(rpd, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, W_FROM, W_TO))
	tmp		<- melt(tmp, measure.vars=c('MALE_SANGER_ID','FEMALE_SANGER_ID'), value.name='ID')
	dr		<- tmp[, {
				pty_run	<- PTY_RUN #<- 7
				w_from	<- W_FROM #<- 1750
				id		<- ID #<- '16033_1_76'
				z		<- subset(dtrees, PTY_RUN==pty_run & W_FROM==w_from)[, IDX]
				cat('\nRun',pty_run,'Window',w_from,'ID',id,'tree index',z)
				ph		<- phs[[ z ]]				
				#	collect branch lengths
				z	<- ph$tip.label[!grepl(id, ph$tip.label)]
				z	<- drop.tip(ph,z,root.edge=0)	
				brls<- z$edge.length
				#	get length of root edge
				z	<- getMRCA(ph, which(grepl(id, ph$tip.label)))
				rl	<- ph$edge.length[ which(ph$edge[,2]==z) ]
				#	get empirical p-value
				z	<- brls[ brls>1e-5 ]
				p	<- length(which(z>=rl))	 / length(z)
				list(ROOTL_EP=p, BRL_MAX=max(brls)) 
			}, by=c('PTY_RUN','W_FROM','W_TO','variable','ID')]
	dr[, variable:=NULL]
	setnames(dr, c('ID','ROOTL_EP','BRL_MAX'), c('MALE_SANGER_ID','MALE_ROOTL_EP','MALE_BRL_MAX'))	
	rpd		<- merge(rpd, dr, by=c('PTY_RUN','W_FROM','W_TO','MALE_SANGER_ID'))
	setnames(dr, c('MALE_SANGER_ID','MALE_ROOTL_EP','MALE_BRL_MAX'), c('FEMALE_SANGER_ID','FEMALE_ROOTL_EP','FEMALE_BRL_MAX'))
	rpd		<- merge(rpd, dr, by=c('PTY_RUN','W_FROM','W_TO','FEMALE_SANGER_ID'))	
	tmp		<- rpd[, list(	CHERRY_EP= max(MALE_ROOTL_EP, FEMALE_ROOTL_EP),
							CHERRY_BRL_MAX=max(MALE_BRL_MAX,FEMALE_BRL_MAX)), by=c('PTY_RUN','W_FROM','W_TO','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	rpd		<- merge(rpd, tmp, by=c('PTY_RUN','W_FROM','W_TO','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	
	tmp		<- subset(rpd, CHERRY_BRL_MAX<=.05) #exclude couples with super long internal branch lengths
	tmp		<- melt(tmp, measure.vars=c('PATRISTIC_DIST','CHERRY_P','CHERRY_EP'), id.vars=c('W_FROM','LEGEND','MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN','CLASS_OR','COUPID_RNK'))
	tmp		<- subset(tmp, variable!='CHERRY_P')
	tmp		<- subset(tmp, !is.na(value))
	set(tmp, NULL, 'variable', tmp[, factor(variable, 	levels=c('PATRISTIC_DIST','CHERRY_P','CHERRY_EP'),
														labels=c('patristic distance\n(subst/site)','length of root edge < internal branch lengths\n(probability under fitted Weibull)','length of root edge < internal branch lengths\n(empirical probability)'))])				
	ggplot(tmp, aes(x=W_FROM, y=value)) +
			scale_y_continuous(breaks=seq(0,1,0.1), minor_breaks=seq(0,1,0.02)) +
			geom_point(aes(colour=LEGEND), size=2) + geom_text(aes(label=COUPID_RNK), size=1) +			
			#geom_hline(aes(yintercept=0.01)) +
			#geom_hline(aes(yintercept=0.02)) +
			#scale_y_log10(breaks=c(0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2)) +
			facet_grid(variable~CLASS_OR, scales='free_y') + theme_bw()
	ggsave(file=file.path(dir, paste(run,'-phsc-windowassignments_SEROINC_cherries_metrics.pdf',sep='')), w=17, h=10)
	#
	#	some tables
	#
	subset(rpd, CHERRY_BRL_MAX<=.05)[, table(CLASS_OR, PATRISTIC_DIST<.02)]	
	subset(rpd, CHERRY_BRL_MAX<=.05)[, table(CLASS_OR, CHERRY_EP>.5)]
}

RakaiCouples.analyze.couples.161213.trmw<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161213/RCCS_161213_w270_trmw_assignments.rda')
	rpw		<- subset(rpw, RUN%in%c("RCCS_161213_w270_d50_p001_mr20_mt2_cl1_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_cl2_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_cl4_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_clInf_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_clInf_dInf") )
	#
	set(rpw, NULL, 'RUN', rpw[, factor(RUN, levels=c( 	"RCCS_161213_w270_d50_p001_mr20_mt2_cl1_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_cl2_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_cl4_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_clInf_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_clInf_dInf"))])
	rpw[, table(RUN, useNA='if')]
	#
	#	for each run: plot pairs	
	#	
	#	define plotting order: largest number of trm assignments
	tmp		<- rpw[, list( WIN_TR=length(which(grepl('close|anc',TYPE))) ), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	tmp		<- tmp[order(-WIN_TR),]
	tmp[, PLOT_ID:=seq_len(nrow(tmp))]	
	#	define label
	tmp		<- merge(tmp, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair ', COUPID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- merge(subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL)), rpw, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, 	levels=c("close_anc_mf", "close_anc_fm", "close_cher_unint",'anc_mf','anc_fm','unint','int','cher','disconnected'), 
												labels=c('M close & ancestral to F','F close & ancestral to M','M, F are close & unint/cherry',
														 'M ancestral to F','F ancestral to M','M, F are unint','M, F are intermingled','M, F are a cherry','M, F are disconnected'))])	
	rpwp	<- subset(tmp, select=c(RUN, PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, LABEL, W_FROM,  W_TO, TYPE, ID1_R, ID1_L, ID2_R, ID2_L ))
	tmp		<- rpwp[, list(	ID_R_MIN=min(ID1_R, ID2_R),
							ID_R_MAX=max(ID1_R, ID2_R)), by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM')]
	rpwp	<- merge(rpwp, tmp, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM'))
	set(rpwp, NULL, 'RUN', rpwp[, gsub('_','\n',gsub('w270_','',gsub('RCCS_','',as.character(RUN))))])
	set(rpwp, NULL, 'RUN', rpwp[, factor(RUN, levels=c(	"161213\nd50\np001\nmr20\nmt2\ncl1\nd7", "161213\nd50\np001\nmr20\nmt2\ncl2\nd7",    
														"161213\nd50\np001\nmr20\nmt2\ncl4\nd7","161213\nd50\np001\nmr20\nmt2\nclInf\nd7","161213\nd50\np001\nmr20\nmt2\nclInf\ndInf"))])
	
	run		<- 'RCCS_161213_w270_dxxx'
	dir		<- rpw$DIR[1]
	#	SEROINC
	tmp		<- subset(rpwp, COUP_SC=='seroinc')
	cols	<- c(	'M close & ancestral to F'="#542788",
					'M ancestral to F'="#8073AC",
					'F close & ancestral to M'="#8C510A",
					'F ancestral to M'="#BF812D",
					'M, F are close & unint/cherry'="#B2182B",
					'M, F are a cherry'="#D6604D",
					'M, F are unint'="#F4A582",
					'M, F are intermingled'="#2166AC",						
					'M, F are disconnected'='grey50')
	p		<- lapply(tmp[, unique(LABEL)], function(label)
			{
				p	<- ggplot(subset(tmp, LABEL==label), aes(x=W_FROM)) +			
						geom_bar(aes(y=ID_R_MAX, colour=TYPE), stat='identity', fill='transparent') +
						geom_bar(aes(y=ID_R_MIN, fill=TYPE), stat='identity', colour='transparent') +
						labs(x='window start', y='number of reads', fill='topology of clades\nbetween patient pairs', colour='topology of clades\nbetween patient pairs') +
						scale_fill_manual(values=cols) +
						scale_colour_manual(values=cols) +
						scale_x_continuous(breaks=seq(0,1e4,250), limits=c(tmp[, min(W_FROM)], tmp[, max(W_FROM)])) +
						scale_y_log10(breaks=c(10,100,1000,1e4,1e5)) +
						theme_bw() + theme(legend.position='top') +
						facet_grid(RUN~LABEL) +
						guides(fill=guide_legend(ncol=2)) 
				p
			})
	pdf(file=file.path(dir, paste(run,'-phsc-windowassignments_SEROINC.pdf',sep='')), w=25, h=2*tmp[, length(unique(RUN))])
	for(i in seq_along(p))	print(p[[i]])
	dev.off()
	#
	#	get distribution of cherries on manually identified unlinked and linked
	#
	tmp	<- as.data.table(matrix( c('83','15715_1_17','15715_1_18',
							'36','15715_1_17','15715_1_18',
							'85','15699_1_38','15699_1_40',
							'104','15699_1_3','15699_1_4',
							'106','15736_1_28','15102_1_16',
							'121','16059_1_66','16059_1_52',
							'31','15981_1_31','15981_1_33',
							'98','15916_1_7','15916_1_21',
							'64','15916_1_7','15916_1_21',
							'122','15910_1_76','15910_1_82',
							'24','15981_1_17','15835_1_32',
							'70','15833_1_84','15833_1_83',
							'104','16059_1_66','16059_1_52',
							'24','16016_1_11','16017_1_44',
							'24','15915_1_51','15776_1_32',
							'15','16033_1_76','16033_1_60',
							'7','16033_1_76','16033_1_60',
							'18','16019_1_45','16018_1_44',
							'67','15892_1_3','15892_1_1',
							'91','16021_1_66','16021_1_58',
							'29','16021_1_66','16021_1_58',
							'116','15915_1_84','15915_1_83',
							'64','15915_1_84','15915_1_83',
							'50','15115_1_4','15115_1_9',
							'86','15978_1_12','15777_1_42',
							'68','15978_1_12','15777_1_42',
							'44','16056_1_83','15099_1_59'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp[, CLASS_OR:= 'OR_disconnected']	
	tmp2	<- as.data.table(matrix( c('111','15950_1_27','15950_1_25',
							'122','16056_1_58','16056_1_57',
							'38','15105_1_35','16016_1_4',
							'113','15950_1_5','15958_1_47',				
							'99','15915_1_61','15915_1_74',
							'117','15915_1_61','15915_1_74',
							'60','16056_1_83','16056_1_85',
							'39','16034_1_86','16034_1_82',
							'120','15981_1_19','15981_1_23',				
							'79','15978_1_36','15978_1_41'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp2[, CLASS_OR:= 'OR_ambiguous']
	tmp		<- rbind(tmp, tmp2)
	tmp2	<- as.data.table(matrix( c('101','16003_1_76','16003_1_79',
							'37','15114_1_64','15758_1_73',
							'37','16016_1_5','16016_1_4',
							'20','16057_1_15','16057_1_16',
							'59','15699_1_5','15867_1_21',
							'99','15915_1_78','15915_1_64',
							'109','15777_1_29','15777_1_21',
							'63','16003_1_86','16003_1_85',
							'109','15776_1_19','15776_1_28',
							'27','16219_1_78','16219_1_80',
							'120','16059_1_73','16059_1_53',
							'25','15675_1_87','15675_1_77',
							'70','15835_1_21','15835_1_22'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp2[, CLASS_OR:= 'OR_transmission']
	tmp		<- rbind(tmp, tmp2)
	set(tmp, NULL, 'PTY_RUN', tmp[, as.integer(PTY_RUN)])
	rpi		<- copy(tmp)
	ubset(rps, COUP_SC=='F->M')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int']), CHER=WIN_OF_TYPE_N[TYPE=='cher'], NO_ANC=sum(WIN_OF_TYPE_N[TYPE=='disconnected'|TYPE=='unint'])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(rps, COUP_SC=='M->F')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='trans_mf'|TYPE=='trans_mf'|TYPE=='int']), CHER=WIN_OF_TYPE_N[TYPE=='cher'], NO_ANC=sum(WIN_OF_TYPE_N[TYPE=='disconnected'|TYPE=='unint'])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= ANY_ANC]
	rpa		<- melt(rpa, measure.vars=c('ANY_ANC', 'NO_ANC','CHER'))
	set(rpa, rpa[, which(variable=='ANY_ANC')], 'variable', 'transmission assignments or intermingled')
	set(rpa, rpa[, which(variable=='NO_ANC')], 'variable', 'unint or disconnected')
	set(rpa, rpa[, which(variable=='CHER')], 'variable', 'cherry')
	setkey(rpa, PAIR_ID)	
	tmp		<- unique(subset(rpa, RUN=='RCCS_161107_w270_d200_r004_mr1'))[order(-WIN_TRM),][, list(COUPID=COUPID, MALE_SANGER_ID=MALE_SANGER_ID, FEMALE_SANGER_ID=FEMALE_SANGER_ID, PTY_RUN=PTY_RUN, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PTY_RUN','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set2') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~RUN) +
			coord_flip() +
			labs(	x= '\nsequence pairs of couples', 
					y='number of windows\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trmpair_windows.pdf',sep='')), w=35, h=10)
	#
	#	define three transmission evidence groups manually
	#	
	rpi	<- subset(rpa, RUN=='RCCS_161107_w270_d200_r004_mr1', c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID))
	setkey(rpi, CLASS_RANK)
	tmp	<- as.data.table(matrix( c('5','15861_1_13','15861_1_18',
				'118','15097_1_88','15714_1_60',
				'56','15743_1_11','16002_1_5',
				'62','15777_1_41','15977_1_69',
				'7','16057_1_2','15834_1_58',
				'1','15861_1_13','15861_1_18',
				'111','16057_1_2','15834_1_58',
				'7','15861_1_13','15103_1_68',
				'8','15949_1_52','15714_1_60',
				'37','15950_1_6','15950_1_16',
				'10','16002_1_17','15758_1_66',
				'81','16035_1_3','15758_1_83',
				'105','16035_1_3','15758_1_83',
				'32','15861_1_13','15103_1_68'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp[, CLASS_OR:= 'OR_disconnected']	
	tmp2	<- as.data.table(matrix( c('66','15778_1_82','15758_1_76',
				'71','16005_1_55','16005_1_54',
				'11','16056_1_74','15736_1_34',				
				'10','15861_1_26','15080_1_12',
				'12','15861_1_20','15861_1_24',
				'64','16005_1_55','15097_1_82',
				'97','16005_1_55','15097_1_82',				
				'80','15080_1_3','15758_1_66'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp2[, CLASS_OR:= 'OR_ambiguous']
	tmp		<- rbind(tmp, tmp2)
	tmp2	<- as.data.table(matrix( c('110','15806_1_32','15806_1_27',
				'78','16016_1_43','16017_1_42',
				'62','15861_1_26','15861_1_22',
				'42','16056_1_73','15736_1_7',
				'102','15736_1_20','15880_1_7',
				'26','15805_1_22','15805_1_23',
				'98','15964_1_3','15758_1_84',
				'91','15103_1_74','15861_1_22',
				'93','15804_1_81','15804_1_60',
				'77','15804_1_53','15804_1_57',
				'16','15097_1_51','15805_1_23',
				'80','15103_1_74','15080_1_12'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp2[, CLASS_OR:= 'OR_transmission']
	tmp		<- rbind(tmp, tmp2)
	set(tmp, NULL, 'PTY_RUN', tmp[, as.integer(PTY_RUN)])
	rpi		<- unique(merge(rpi, tmp, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')))
	rps		<- merge(rps, rpi, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all.x=1)
	#
	#	get false pos and true pos for the two unambiguous groups	
	#
	stopifnot( nrow(subset(rps, (COUP_SC=='F->M'|COUP_SC=='M->F') & is.na(CLASS_OR)))==0 )
	tmp		<- subset(rps, (COUP_SC=='F->M'|COUP_SC=='M->F') & CLASS_OR=='OR_disconnected')
	rpb		<- tmp[, list(	POT_FP= sum(WIN_OF_TYPE_N[TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int']), 
							POT_TN=sum(WIN_OF_TYPE_N[TYPE=='disconnected'|TYPE=='unint']),
							POT_TP=0,
							POT_FN=0), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(rps, (COUP_SC=='F->M'|COUP_SC=='M->F') & CLASS_OR=='OR_transmission')
	tmp		<- tmp[, list(	POT_TP= sum(WIN_OF_TYPE_N[TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int']), 
							POT_FN=sum(WIN_OF_TYPE_N[TYPE=='disconnected'|TYPE=='unint']),
							POT_FP=0,
							POT_TN=0), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpb		<- rbind(rpb, tmp, fill=TRUE, use.names=TRUE)
	rpb[, TOTAL:= POT_FP+POT_FN+POT_TP+POT_TN]
	#	by choice of the three groups, there are only TN and TPs.. ..so this is not interesting here.
	rpb[, Maj_Class:= NA_character_]
	set(rpb, rpb[, which(POT_FP<POT_TN)], 'Maj_Class', 'TN')
	set(rpb, rpb[, which(POT_FP>=POT_TN)], 'Maj_Class', 'FP')
	set(rpb, rpb[, which(POT_TP>POT_FN)], 'Maj_Class', 'TP')
	set(rpb, rpb[, which(POT_TP<=POT_FN)], 'Maj_Class', 'FN')					
	#
	#
	tmp		<- melt(rpb, id.vars=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC','TOTAL'), measure.vars=c('POT_FP','POT_FN'))	
	tmp		<- tmp[, list(Potential_False_Class= sum(value), Potential_True_Class= sum(TOTAL)-sum(value) ), by=c('RUN','variable')]	
	set(tmp, tmp[,which(variable=='POT_FP')], 'variable', 'Manually identified unlinked pair')
	set(tmp, tmp[,which(variable=='POT_FN')], 'variable', 'Manually identified trm pair')	
	tmp[, Total_Class:= Potential_False_Class+Potential_True_Class]
	
	ggplot( tmp, aes(x=RUN, y=Potential_False_Class/Total_Class, fill=variable) ) + 
			geom_bar(stat='identity') + 
			scale_fill_brewer(palette='Set1') +
			scale_y_continuous(breaks=seq(0,1,0.02), labels=percent) +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~variable) +
			coord_flip() +
			labs(	y= '\nproportion of "likely false" window assignments\ne.g. % of windows from manually identified trm pairs that are falsely classified as unlinked', 
					x='run\n',
					colour='',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')	
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trmpair_windows_falsepos.pdf',sep='')), w=10, h=5)
	#	number of couples with less than 5 trm assignments
	tmp		<- unique(subset(rpb, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID)))
	tmp2	<- unique(subset(rpb, select=RUN))
	tmp[, D:=1]
	tmp2[, D:=1]
	tmp		<- merge(tmp, tmp2, by='D', allow.cartesian=TRUE)
	tmp[, D:=NULL]
	tmp		<- merge(tmp, rpb, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all.x=1)
	set(tmp, tmp[, which(is.na(TOTAL))], 'TOTAL', 0)	
	tmp		<- subset(tmp, TOTAL<5)[, list(WIN_L_FIVE=length(PAIR_ID)), by='RUN']	
	setkey(tmp, RUN)
	#	also track pairs that cannot be any longer classified
	
	#
	#	plot number of directional trm assignments in only possible direction 
	#	
	rpa		<- subset(rps, COUP_SC=='F->M')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='trans_fm']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='trans_mf'])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(rps, COUP_SC=='M->F')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='trans_mf']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='trans_fm'])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= IN_DIR+AGAINST_DIR]
	rpa		<- melt(rpa, measure.vars=c('IN_DIR', 'AGAINST_DIR'))
	set(rpa, rpa[, which(variable=='IN_DIR')], 'variable', 'trm assignment in the only epidemiologically possible direction')
	set(rpa, rpa[, which(variable=='AGAINST_DIR')], 'variable', 'trm assignment against the only epidemiologically possible direction')	
	setkey(rpa, PAIR_ID)		
	tmp		<- unique(subset(rpa, RUN==rpa$RUN[1]))[order(-WIN_TRM),][, list(COUPID=COUPID, MALE_SANGER_ID=MALE_SANGER_ID, FEMALE_SANGER_ID=FEMALE_SANGER_ID, PTY_RUN=PTY_RUN, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PTY_RUN','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set1') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~RUN) +
			coord_flip() +
			labs(	x= '\nsequence pairs of couples', 
					y='number of transmission windows with direction\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trmdir_windows.pdf',sep='')), w=35, h=10)
	#
	#	get false pos and true pos for the manual classification group
	#
	tmp		<- unique(subset(rps, RUN=='RCCS_161107_w270_d20_r004_mr100' & CLASS_OR=='OR_transmission',c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')))
	rpa		<- merge(rpa, tmp, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	tmp		<- rpa[, {
							list( TRM_OK_MAJ= as.integer(sum(value[variable=='trm assignment in the only epidemiologically possible direction'])>sum(value[variable!='trm assignment in the only epidemiologically possible direction'])), 
								  TRM_OK_66pc= as.integer(sum(value[variable=='trm assignment in the only epidemiologically possible direction']) / sum(value) >= 2/3))
						}, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	rpa		<- merge(rpa, tmp, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	tmp		<- rpa[, {
				tmp	<- sum(value[variable=='trm assignment in the only epidemiologically possible direction'])				
				list(	V=	c( tmp / (sum(WIN_TRM)/2), 1 - tmp / (sum(WIN_TRM)/2), (sum(TRM_OK_MAJ)/2) / (length(TRM_OK_MAJ)/2), 1 - (sum(TRM_OK_MAJ)/2) / (length(TRM_OK_MAJ)/2), (sum(TRM_OK_66pc)/2) / (length(TRM_OK_66pc)/2), 1-(sum(TRM_OK_66pc)/2) / (length(TRM_OK_66pc)/2)),
						TYPE= c('TP', 'FP', 'TP', 'FP', 'TP', 'FP'),
						RULE= c('by window', 'by window', 'by majority rule', 'by majority rule', 'by 66% majority', 'by 66% majority')
						)
			}, by=c('RUN')]	
	ggplot( tmp, aes(x=RUN, y=V, fill= TYPE) ) + 
			geom_bar(stat='identity', position='stack') + 
			scale_fill_brewer(palette='Set1') +
			scale_y_continuous(breaks=seq(0,1,0.1), labels=percent) +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			coord_flip() +
			labs(	y= '\nproportion of "false positives" and "false negatives"', 
					x='run\n',
					colour='',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n') +
			facet_grid(~RULE)
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trm_windows_falsepos.pdf',sep='')), w=14, h=6)
	
	tmp		<- unique(subset(rpa, !TRM_OK_66pc, c('RUN','PTY_RUN','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID')))
	setkey(tmp, RUN)
	#
	#
	#
	#
}

RakaiCouples.analyze.couples.161219.pair.bb.model<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	#
	rpw[, table(RUN, useNA='if')]
	#	define plotting order: largest number of trm assignments	
	tmp		<- rpw[, list( WIN_TR=length(which(grepl('close|anc',TYPE_DIR_TODI7x3))) ), by=c('PTY_RUN','COUP_SC','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	tmp		<- tmp[order(-WIN_TR),]
	tmp[, PLOT_ID:=seq_len(nrow(tmp))]	
	#	define label
	tmp		<- merge(tmp, rp, by=c('COUP_SC','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	setkey(tmp, PLOT_ID)
	tmp[, LABEL_SH:= factor(PLOT_ID, levels=PLOT_ID, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair ', COUPID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL, LABEL_SH))
	rpw		<- merge(tmp, rpw, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	#	define min/max reads
	tmp		<- rpw[, list(	ID_R_MIN=min(ID1_R, ID2_R),
					ID_R_MAX=max(ID1_R, ID2_R)), by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM')]
	rpw		<- merge(rpw, tmp, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM'))
	#
	#	define topo/distance types of phylogenetic relationship between pairs
	#
	#	given the bimodal distribution of patristic distances from phyloscanner
	#	I ususally consider 3 distance states: close, neither close nor distant, distant
	#
	#	TYPE_DIR_TODI7x3		7 topology states (this is the "basic topological characterization" on Slack that I derive from Matthew's output) 
	#							chain_12_nointermediate, chain_12_withintermediate, chain_21_nointermediate, chain_21_withintermediate, intermingled_nointermediate, intermingled_withintermediate, other  
	#							times 
	#							3 distance states
	#
	#	TYPE_PAIR_TODI3x3		3 topology states
	#							pair (chain_XX_nointermediate+intermingled_nointermediate), withintermediate (XX_withintermediate), other
	#							times 
	#							3 distance states
	#
	#	TYPE_PAIR_TODI2x2		2 topology states
	#							ancestral/intermingled, other
	#							times 
	#							2 distance states (close, not close)
	#
	#	TYPE_PAIR_TODI3			non-factorial design that combines distant or withintermediate, and close or ancestral/intermingled
	#							the idea here is that 
	#							evidence for extra-couple transmission comes from large patristic distance OR intermediates	
	#							evidence for extra-couple transmission comes from ancestral/intermingled OR short patristic distance	
	#
	#	TYPE_PAIR_DI			3 distance states
	#
	#
	#	TYPE_DIR_TO				3 topological direction states
	#							chain_fm	(chain_fm_nointermediate regardless of distance); chain_mf	(chain_mf_nointermediate regardless of distance); other
	#							
	rpw[, TYPE_PAIR_TO:= 'other']
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair_close','pair','pair_distant'))], 'TYPE_PAIR_TO', 'ancestral/\nintermingled')
	rpw[, TYPE_PAIR_TODI3:= 'with intermediate\nor distant']
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair_close','other_close'))], 'TYPE_PAIR_TODI3', 'no intermediate\n and close')
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair','other'))], 'TYPE_PAIR_TODI3', 'no intermediate\n but not close')
	rpw[, TYPE_PAIR_DI:= cut(PATRISTIC_DISTANCE, breaks=c(1e-12,0.02,0.05,2), labels=c('close','intermediate\ndistance','distant'))] 
	tmp		<- 	rpw[, list(PD_MEAN=mean(PATRISTIC_DISTANCE)), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	rpw[, TYPE_PAIR_TODI2x2:= 'not close other']
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair_close'))], 'TYPE_PAIR_TODI2x2', 'close ancestral/\nintermingled')
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair','pair_distant'))], 'TYPE_PAIR_TODI2x2', 'not close ancestral/\nintermingled')
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('withintermediate_close','other_close'))], 'TYPE_PAIR_TODI2x2', 'close other')	
	tmp[, TYPE_PD_MEAN:= cut(PD_MEAN, breaks=c(1e-12,0.02,0.05,2), labels=c('close','ambiguous','distant'))]				
	rpw		<- merge(rpw, tmp, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rpw[, TYPE_DIR_TO3:= NA_character_]
	set(rpw, rpw[, which(grepl('chain_fm',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TO3', 'fm')
	set(rpw, rpw[, which(grepl('chain_mf',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TO3', 'mf')
	set(rpw, rpw[, which(grepl('intermingled',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TO3', 'intermingled')		
	rpw[, TYPE_DIR_TODI3:= NA_character_]
	set(rpw, rpw[, which(grepl('chain_fm',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'fm')
	set(rpw, rpw[, which(grepl('chain_mf',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'mf')
	set(rpw, rpw[, which(grepl('intermingled',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'intermingled')				
	set(rpw, NULL, 'TYPE_DIR_TODI7x3', rpw[, gsub('intermediate',' intermediate',gsub('intermediate_','intermediate\n',gsub('intermingled_','intermingled\n',gsub('(chain_[fm][mf])_','\\1\n',TYPE_DIR_TODI7x3))))])
	#	define effectively independent number of windows
	#	select only W_FROM that are > W_TO	
	tmp		<- rpw[, {
				z		<- 1
				repeat
				{
					zz		<- which(W_FROM>max(W_TO[z]))[1]
					if(length(zz)==0 | is.na(zz))	break
					z		<- c(z, zz)			
				}
				list(W_FROM=W_FROM[z], W_TO=W_TO[z])
			}, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN')]
	tmp[, OVERLAP:= 0L]
	rpw		<- merge(rpw, tmp, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','W_FROM','W_TO','RUN'), all.x=1 )
	set(rpw, rpw[, which(is.na(OVERLAP))], 'OVERLAP', 1L)		
	rpw		<- melt(rpw, measure.vars=c('TYPE_PD_MEAN','TYPE_PAIR_TO','TYPE_PAIR_TODI3x3','TYPE_PAIR_TODI3','TYPE_PAIR_DI','TYPE_PAIR_TODI2x2','TYPE_DIR_TODI7x3','TYPE_DIR_TO3','TYPE_DIR_TODI3'), variable.name='GROUP', value.name='TYPE')
	set(rpw, NULL, 'TYPE', rpw[,gsub('_',' ', TYPE)])	
	rpw		<- subset(rpw, !is.na(TYPE))
	#
	#	for each pair count windows by TYPE_PAIR_DI (k)
	#
	rplkl	<- rpw[, list(K=length(W_FROM)), by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','COUP_SC','LABEL','LABEL_SH','GROUP','TYPE')]	
	rplkl	<- dcast.data.table(rplkl, RUN+PTY_RUN+FEMALE_SANGER_ID+MALE_SANGER_ID+COUP_SC+LABEL+LABEL_SH~GROUP+TYPE, value.var='K')
	for(x in setdiff(colnames(rplkl),c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','COUP_SC','LABEL','LABEL_SH','GROUP')))
		set(rplkl, which(is.na(rplkl[[x]])), x, 0L)	
	for(x in colnames(rplkl)[grepl('TYPE_PD_MEAN',colnames(rplkl))])
		set(rplkl, which(rplkl[[x]]>0), x, 1L)	
	rplkl	<- melt(rplkl, id.vars=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','COUP_SC','LABEL','LABEL_SH'), variable.name='GROUP', value.name='K')
	rplkl[, TYPE:= gsub('.*_([^_]+)$','\\1',GROUP)]
	set(rplkl, NULL, 'GROUP', rplkl[, gsub('(.*)_[^_]+$','\\1',GROUP)])	
	#	for each pair count total windows (n) and effective windows (neff)	
	tmp		<- rpw[, list(N=length(W_FROM), NEFF=sum(1-OVERLAP)), by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN','GROUP')]
	rplkl	<- merge(rplkl, tmp, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN','GROUP'))
	#	define effective number of windows of type
	rplkl[, KEFF:= K/N*NEFF]	
	#
	#	add prior (equal probabilities to all types), posterior, and marginal probabilities
	#
	rplkl[, DIR_PRIOR:= 0.1]
	rplkl[, DIR_PO:= DIR_PRIOR+KEFF]
	tmp		<- rplkl[, {
				alpha	<- DIR_PO
				beta	<- sum(DIR_PO)-DIR_PO				
				list(	TYPE=TYPE, 
						POSTERIOR_ALPHA=alpha, 
						POSTERIOR_BETA=beta, 
						POSTERIOR_MEAN=alpha/(alpha+beta), 
						POSTERIOR_QL=qbeta(0.025, alpha, beta), POSTERIOR_QU=qbeta(0.975, alpha, beta),
						POSTERIOR_IL=qbeta(0.25, alpha, beta), POSTERIOR_IU=qbeta(0.75, alpha, beta))	
			}, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN','GROUP')]
	rplkl	<- merge(rplkl, tmp, by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','GROUP','TYPE'))
	#	fixup TYPE_PD_MEAN: there are no posteriors here
	set(rplkl, rplkl[, which(GROUP=='TYPE_PD_MEAN' & K==0)], c('POSTERIOR_MEAN'), 0)
	set(rplkl, rplkl[, which(GROUP=='TYPE_PD_MEAN' & K>0)], c('POSTERIOR_MEAN'), 1)
	set(rplkl, rplkl[, which(GROUP=='TYPE_PD_MEAN')], c('KEFF','DIR_PO','DIR_PRIOR','POSTERIOR_ALPHA','POSTERIOR_BETA','POSTERIOR_QL','POSTERIOR_QU','POSTERIOR_IL','POSTERIOR_IU'), NA_real_)
	
	#	save
	save(rplkl, file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmp_bbmodels.rda')
	
	#
	#	make table by TYPE of effective counts,  posterior probs,	alpha/beta
	#	as on blackboard
	#
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]
	for(group in c('TYPE_PAIR_TODI2x2','TYPE_PAIR_TODI3'))
	{
		rpt	<- subset(rplkl, GROUP==group)
		rpt[, POSTERIOR:= paste0( round(100*POSTERIOR_MEAN, d=1), '% (', round(100*POSTERIOR_IL, d=1), '%, ',round(100*POSTERIOR_IU, d=1),'%)')]
		suppressWarnings({rpt	<- melt(rpt, measure.vars=c('KEFF','POSTERIOR'))})
		rpt	<- dcast.data.table(rpt, RUN+PTY_RUN+FEMALE_SANGER_ID+MALE_SANGER_ID+LABEL+LABEL_SH+COUP_SC~variable+TYPE, value.var='value')
		setkey(rpt, COUP_SC, LABEL_SH)
		set(rpt, NULL, c('LABEL_SH'), NULL)
		write.csv(rpt, row.names=FALSE, file=file.path(dir, paste0(run,'_table_',group,'.csv')))		
	}	 
	#
	#	prepare colours
	#	
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
	#
	#	make detailed plots for each pair
	#		topology assignments across genome
	#		distances across genome
	#		number of windows of certain type
	#		estimated posterior probabilities on unknown phylogenetic relationship
	#
	groups		<- c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3','TYPE_PAIR_TODI2x2','TYPE_PAIR_TODI3','TYPE_PAIR_DI')
	#groups		<- c('TYPE_DIR_TODI7x3')
	coup_sc		<- 'seroinc'
	coup_sc		<- c('M->F','F->M')
	for(group in groups)
	{
		widths	<- unit(c(4, 6), "null")
		heights	<- unit(c(2, 3.5, 4, 5), "null")
		height	<- 9
		if(group%in%c('TYPE_DIR_TODI7x3'))
		{
			widths	<- unit(c(4, 6), "null")
			heights	<- unit(c(2, 3.5, 4, 15), "null")
			height	<- 17
		}		
		if(group%in%c('TYPE_PAIR_TODI2x2'))
		{
			heights	<- unit(c(2, 3.5, 4, 3.75), "null")
			height	<- 8
		}
		if(group%in%c('TYPE_PAIR_TODI3','TYPE_PAIR_DI'))
		{
			heights	<- unit(c(2, 3.5, 4, 3.5), "null")
			height	<- 7
		}
		pdf(file=file.path(dir, paste(run,'-phsc-relationships_',ifelse(any(coup_sc=='M->F'), 'serodisc', coup_sc),'_',group,'.pdf',sep='')), w=10, h=height)	
		plot.tmp	<- unique(subset(rplkl, COUP_SC%in%coup_sc & GROUP==group, c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','LABEL')), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setkey(plot.tmp, LABEL)		
		for(i in seq_len(nrow(plot.tmp)))
		{
			#pty_run	<- 38; id1		<- '16016_1_4'; id2		<- '15105_1_35'
			pty_run	<- plot.tmp[i, PTY_RUN]; id1		<- plot.tmp[i, FEMALE_SANGER_ID]; id2		<- plot.tmp[i, MALE_SANGER_ID]		
			tmp		<- subset(rpw, PTY_RUN==pty_run & GROUP==group & FEMALE_SANGER_ID==id1 & MALE_SANGER_ID==id2)
			p1		<- ggplot(tmp, aes(x=W_FROM)) +			
					geom_bar(aes(y=ID_R_MAX, colour=TYPE), stat='identity', fill='transparent') +
					geom_bar(aes(y=ID_R_MIN, fill=TYPE), stat='identity', colour='transparent') +
					labs(x='', y='number of reads', fill='phylogenetic\nrelationship\n', colour='phylogenetic\nrelationship\n') +
					scale_fill_manual(values=cols.type[[group]]) +
					scale_colour_manual(values=cols.type[[group]]) +
					scale_x_continuous(breaks=seq(0,1e4,500), minor_breaks=seq(0,1e4,100), limits=c(rpw[, min(W_FROM)], rpw[, max(W_FROM)])) +
					scale_y_log10(breaks=c(10,100,1000,1e4,1e5)) +
					theme_bw() + theme(legend.position='left') +			
					guides(fill=FALSE, colour=FALSE)
			p2		<- ggplot(tmp, aes(x=W_FROM, y=PATRISTIC_DISTANCE)) +
					geom_point(size=1) +					
					labs(x='window start\n\n', y='patristic distance') +
					scale_x_continuous(breaks=seq(0,1e4,500), minor_breaks=seq(0,1e4,100), limits=c(rpw[, min(W_FROM)], rpw[, max(W_FROM)])) +
					scale_y_log10(labels=percent, limits=c(0.001, 0.7), expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
					theme_bw() + theme(legend.position='left')
			tmp		<- subset(rplkl, GROUP==group & PTY_RUN==pty_run & FEMALE_SANGER_ID==id1 & MALE_SANGER_ID==id2)
			p3		<- ggplot(tmp, aes(x=TYPE, y=KEFF, fill=TYPE)) + geom_bar(stat='identity') +
					scale_fill_manual(values=cols.type[[group]]) +
					theme_bw() + theme(legend.position='bottom') +
					coord_flip() + guides(fill=FALSE) +			
					labs(x='', y='\nnon-overlapping windows\n(number)', fill='phylogenetic\nrelationship\n')
			p4		<- ggplot(tmp, aes(x=TYPE, middle=qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA), ymin=POSTERIOR_QL, ymax=POSTERIOR_QU, lower=POSTERIOR_IL, upper=POSTERIOR_IU, fill=TYPE)) + 
					geom_boxplot(stat='identity') +
					scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), limits=c(0,1), expand=c(0,0)) +
					scale_fill_manual(values=cols.type[[group]]) +
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
	#
	#	compare classification 'TYPE_PAIR_DI' to 'TYPE_PAIR_TODI3' 
	#
	g_legend<-function(a.gplot)
	{
		tmp <- ggplot_gtable(ggplot_build(a.gplot))
		leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
		legend <- tmp$grobs[[leg]]
		legend
	}
	
	tmp		<- subset(rplkl, COUP_SC=='seroinc' & GROUP%in%c('TYPE_PAIR_DI','TYPE_PAIR_TODI3'))
	#tmp		<- subset(rplkl, COUP_SC%in%c('M->F','F->M') & GROUP%in%c('TYPE_PAIR_DI','TYPE_PAIR_TODI3'))
	setkey(tmp, LABEL_SH)
	set(tmp, NULL, 'GROUP', tmp[, factor(GROUP, levels=rev(c('TYPE_PAIR_DI','TYPE_PAIR_TODI3','TYPE_PAIR_TODI2x2')))])
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=rev(c("close","intermediate\ndistance","distant","no intermediate\n and close","no intermediate\n but not close","with intermediate\nor distant","close ancestral/\nintermingled","close other","not close ancestral/\nintermingled","not close other")))])
	tmp2	<- c(cols.type[['TYPE_PAIR_DI']],cols.type[['TYPE_PAIR_TODI3']],cols.type[['TYPE_PAIR_TODI2x2']])	
	p1		<- ggplot(tmp, aes(x=GROUP, y=KEFF, fill=TYPE)) + 
			geom_bar(stat='identity',position='stack') +
			scale_y_continuous(expand=c(0,0)) +
			scale_fill_manual(values=tmp2) +
			theme_bw() + 
			theme(	legend.position='bottom', axis.text.y=element_blank(),
					axis.ticks.y=element_blank(),
					panel.spacing=unit(0.4, "lines"), strip.text.y = element_text(angle=180),
					strip.background=element_rect(fill="transparent", colour="transparent"),
					panel.border=element_rect(color="transparent")) +
			facet_grid(LABEL_SH~., switch='y') +
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE)) +
			labs(	x='', 
					y='non-overlapping windows\n(number)\n',
					fill='phylogenetic\nrelationship')	
	p2		<- ggplot(subset(tmp, TYPE%in%c("no intermediate\n and close","close")), aes(x=GROUP, middle=qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA), lower=POSTERIOR_IL, upper=POSTERIOR_IU, ymin=POSTERIOR_QL, ymax=POSTERIOR_QU, fill=TYPE)) + 
			geom_boxplot(stat='identity') +
			scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.2), labels=percent) +
			scale_fill_manual(values=tmp2) +
			theme_bw() + 
			theme(	legend.position='bottom', axis.text.y=element_blank(),
					axis.ticks.y=element_blank(),
					panel.spacing=unit(0.4, "lines"), strip.text.y=element_blank(),
					strip.background=element_blank(),
					panel.border=element_rect(color="transparent")) +
			facet_grid(LABEL_SH~., switch='y') +
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE)) +
			labs(	x='', 
					y='posterior probability\n\n',
					fill='phylogenetic\nrelationship')
	p3		<- g_legend(p1)	
	p3$vp	<- viewport(layout.pos.row=2, layout.pos.col=1:2)	
	pdf(file=file.path(dir, paste(run,'-phsc-relationships_seroinc_compare_DI_to_TODI3.pdf',sep='')), w=12, h=15)
	#pdf(file=file.path(dir, paste(run,'-phsc-relationships_serodisc_compare_DI_to_TODI3.pdf',sep='')), w=12, h=15)
	grid.newpage()	
	pushViewport(viewport(layout = grid.layout(2, 2, heights=unit(c(10,1), "null"), widths=unit(c(7, 3), "null"))))   	
	print(p1+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
	print(p2+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
	grid.draw(p3)
	dev.off()	
	#
	#	compare classification 'TYPE_PAIR_DI' 'TYPE_PAIR_TO' 'TYPE_PAIR_TODI3' next to each other
	#		
	tmp		<- subset(rplkl, COUP_SC=='seroinc' & GROUP%in%c('TYPE_PAIR_DI','TYPE_PAIR_TO','TYPE_PAIR_TODI3'))
	setkey(tmp, LABEL_SH)
	set(tmp, NULL, 'GROUP', tmp[, factor(GROUP, levels=rev(c('TYPE_PAIR_DI','TYPE_PAIR_TO','TYPE_PAIR_TODI3')))])
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=rev(c("close","intermediate\ndistance","distant","no intermediate\n and close","no intermediate\n but not close","with intermediate\nor distant","ancestral/\nintermingled","other")))])
	tmp2	<- c(cols.type[['TYPE_PAIR_DI']],cols.type[['TYPE_PAIR_TO']],cols.type[['TYPE_PAIR_TODI3']])	
	p1		<- ggplot(tmp, aes(x=GROUP, y=KEFF, fill=TYPE)) + 
			geom_bar(stat='identity',position='stack') +
			scale_y_continuous(expand=c(0,0)) +
			scale_fill_manual(values=tmp2) +
			theme_bw() + 
			theme(	legend.position='bottom', axis.text.y=element_blank(),
					axis.ticks.y=element_blank(),
					panel.spacing=unit(0.4, "lines"), strip.text.y = element_text(angle=180),
					strip.background=element_rect(fill="transparent", colour="transparent"),
					panel.border=element_rect(color="transparent")) +
			facet_grid(LABEL_SH~., switch='y') +
			coord_flip() +
			guides(fill=guide_legend(ncol=4, byrow = TRUE)) +
			labs(	x='', 
					y='non-overlapping windows\n(number)\n',
					fill='phylogenetic\nrelationship')	
	p2		<- ggplot(subset(tmp, TYPE%in%c("no intermediate\n and close","ancestral/\nintermingled","close")), aes(x=GROUP, middle=qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA), lower=POSTERIOR_IL, upper=POSTERIOR_IU, ymin=POSTERIOR_QL, ymax=POSTERIOR_QU, fill=TYPE)) + 
			geom_boxplot(stat='identity') +
			scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.2), labels=percent) +
			scale_fill_manual(values=tmp2) +
			theme_bw() + 
			theme(	legend.position='bottom', axis.text.y=element_blank(),
					axis.ticks.y=element_blank(),
					panel.spacing=unit(0.4, "lines"), strip.text.y=element_blank(),
					strip.background=element_blank(),
					panel.border=element_rect(color="transparent")) +
			facet_grid(LABEL_SH~., switch='y') +
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE)) +
			labs(	x='', 
					y='posterior probability\n\n',
					fill='phylogenetic\nrelationship')
	p3		<- g_legend(p1)	
	p3$vp	<- viewport(layout.pos.row=2, layout.pos.col=1:2)	
	pdf(file=file.path(dir, paste(run,'-phsc-relationships_compare_DI_to_TO_to_TODI3.pdf',sep='')), w=12, h=20)
	grid.newpage()	
	pushViewport(viewport(layout = grid.layout(2, 2, heights=unit(c(10,1), "null"), widths=unit(c(7, 3), "null"))))   	
	print(p1+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
	print(p2+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
	grid.draw(p3)
	dev.off()
	#
	#	compare classifications  'TYPE_PAIR_DI' 'TYPE_PAIR_TODI2x2' next to each other
	#
	tmp		<- subset(rplkl, COUP_SC=='seroinc' & GROUP%in%c('TYPE_PAIR_DI','TYPE_PAIR_TODI2x2'))
	#tmp		<- subset(rplkl, COUP_SC%in%c('M->F','F->M') & GROUP%in%c('TYPE_PAIR_DI','TYPE_PAIR_TODI3'))
	setkey(tmp, LABEL_SH)
	set(tmp, NULL, 'GROUP', tmp[, factor(GROUP, levels=rev(c('TYPE_PAIR_DI','TYPE_PAIR_TODI2x2')))])
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=rev(c("close","intermediate\ndistance","distant","close ancestral/\nintermingled","close other","not close ancestral/\nintermingled","not close other")))])
	tmp2	<- c(cols.type[['TYPE_PAIR_DI']],cols.type[['TYPE_PAIR_TODI2x2']])	
	p1		<- ggplot(tmp, aes(x=GROUP, y=KEFF, fill=TYPE)) + 
			geom_bar(stat='identity',position='stack') +
			scale_y_continuous(expand=c(0,0)) +
			scale_fill_manual(values=tmp2) +
			theme_bw() + 
			theme(	legend.position='bottom', axis.text.y=element_blank(),
					axis.ticks.y=element_blank(),
					panel.spacing=unit(0.4, "lines"), strip.text.y = element_text(angle=180),
					strip.background=element_rect(fill="transparent", colour="transparent"),
					panel.border=element_rect(color="transparent")) +
			facet_grid(LABEL_SH~., switch='y') +
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE)) +
			labs(	x='', 
					y='non-overlapping windows\n(number)\n',
					fill='phylogenetic\nrelationship')	
	p2		<- ggplot(subset(tmp, TYPE%in%c("close ancestral/\nintermingled","close")), aes(x=GROUP, middle=qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA), lower=POSTERIOR_IL, upper=POSTERIOR_IU, ymin=POSTERIOR_QL, ymax=POSTERIOR_QU, fill=TYPE)) + 
			geom_boxplot(stat='identity') +
			scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.2), labels=percent) +
			scale_fill_manual(values=tmp2) +
			theme_bw() + 
			theme(	legend.position='bottom', axis.text.y=element_blank(),
					axis.ticks.y=element_blank(),
					panel.spacing=unit(0.4, "lines"), strip.text.y=element_blank(),
					strip.background=element_blank(),
					panel.border=element_rect(color="transparent")) +
			facet_grid(LABEL_SH~., switch='y') +
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE)) +
			labs(	x='', 
					y='posterior probability\n\n',
					fill='phylogenetic\nrelationship')
	p3		<- g_legend(p1)	
	p3$vp	<- viewport(layout.pos.row=2, layout.pos.col=1:2)	
	pdf(file=file.path(dir, paste(run,'-phsc-relationships_seroinc_compare_DI_to_TODI2x2.pdf',sep='')), w=12, h=15)
	#pdf(file=file.path(dir, paste(run,'-phsc-relationships_serodisc_compare_DI_to_TODI3.pdf',sep='')), w=12, h=15)
	grid.newpage()	
	pushViewport(viewport(layout = grid.layout(2, 2, heights=unit(c(10,1), "null"), widths=unit(c(7, 3), "null"))))   	
	print(p1+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
	print(p2+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
	grid.draw(p3)
	dev.off()	
		
}

RakaiAll.analyze.pairs.170120.distances<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ape)
	
	#
	#	quick info on how many couples are close in consensus tree
	#
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]	
	
	
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision/RCCS_SeqInfo_160816.rda')
	#	get epi info
	tmp		<- RakaiCirc.epi.get.info()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#	add sequence dates to rd
	rd		<- unique(subset(rd, !is.na(PID), select=c("RID","PID","SEX","REGION", "COMM_NUM", "HH_NUM","BIRTHDATE","LASTNEGVIS","LASTNEGDATE","FIRSTPOSVIS","FIRSTPOSDATE","RELIGION")), by='PID')
	tmp		<- unique(subset(rs, !is.na(PID), select=c(PID, DATE)),by='PID')
	setnames(tmp, 'DATE','SEQDATE')
	rd		<- merge(rd, tmp, by='PID',all.x=1)
	#	add community types to rd
	tmp		<- data.table(	COMM_NUM=	c("1","2","3","4","5","6","7","8","16","18","19","22","23","24","25","29","33","34","36","38","40","51","54","55", "56","57","58","59","62","74","77","81","89","94","95","103","106","107","108","109","120","177", "370","391","401","451", "602", "754", "760", "770","771","772","773","774","776"),
							COMM_TYPE=	c("T","A","A","T","A","A","A","A", "T", "A", "A", "T", "A", "T", "A", "A", "T", "A", "A", "F", "A", "T", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",  "A","A",   "A",  "T",  "A",  "A",  "A",  "A",   "A",  "A", "A",  "A",    "A",  "A",    "A", "F",  "F",  "A",  "A",  "F",   "T"))
	set(tmp, NULL, 'COMM_TYPE', tmp[, as.character(factor(COMM_TYPE, levels=c('A','T','F'), labels=c('agrarian','trading','fisherfolk')))])
	set(tmp, NULL, 'COMM_NUM', tmp[, as.integer(COMM_NUM)])
	rd		<- merge(rd, tmp, by='COMM_NUM')
	
	#	load patristic distance matrix
	infile	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_UG_gag_coinfinput_160219.rda'
	load(infile)
	#	prepare patristic distance data.table
	ph.gdtr	<- as.data.table(melt(ph.gdtr, varnames=c('TAXA1','TAXA2')))
	setnames(ph.gdtr, 'value', 'PD')
	ph.gdtr	<- subset(ph.gdtr, TAXA1!=TAXA2)
	set(ph.gdtr, NULL, 'TAXA1', ph.gdtr[, gsub('_','-',as.character(TAXA1))])
	set(ph.gdtr, NULL, 'TAXA2', ph.gdtr[, gsub('_','-',as.character(TAXA2))])
	#	load genetic distance matrix with overlap
	infile		<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment_gd.rda'
	load(infile)	#loads sq.gd
	setnames(sq.gd, 'PD', 'GD')
	#	load selected phyloscanner runs
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	#	for each run construct all pairings that were evaluated
	ptc			<- lapply(pty.runs[, unique(PTY_RUN)], function(pty_run)
			{
				#pty_run<- 1
				z			<- subset(pty.runs, PTY_RUN==pty_run)[, TAXA]
				tmp			<- subset(sq.gd, TAXA1%in%z & TAXA2%in%z)	# this is symmetric
				tmp[, PTY_RUN:=pty_run]				
				tmp2		<- subset(ph.gdtr, TAXA1%in%z & TAXA2%in%z)
				tmp			<- merge(tmp, tmp2, by=c('TAXA1','TAXA2'), all.x=1)
				tmp				
			})
	ptc			<- do.call('rbind',ptc)
	ptc[, MALE_PID:= gsub('-S[0-9]+','',TAXA1)]
	ptc[, FEMALE_PID:= gsub('-S[0-9]+','',TAXA2)]
	ptcc		<- copy(ptc)
	#	merge with Rakai IDs and make sure male is first
	tmp			<- subset(rd, select=c(PID, RID, SEX))
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	ptc			<- merge(ptc, tmp, by='MALE_PID')
	setnames(tmp, colnames(tmp), gsub('MALE','FEMALE',colnames(tmp)))
	ptc			<- merge(ptc, tmp, by='FEMALE_PID')
	ptc			<- subset(ptc, MALE_SEX=='M' & FEMALE_SEX=='F')
	ptc			<- unique(ptc, by=c('TAXA1','TAXA2'))
	#	merge with couples	
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	tmp			<- unique(subset(rp, select=c(MALE_RID,FEMALE_RID,COUPID,COUP_SC)))
	ptc			<- merge(ptc, tmp, by=c('FEMALE_RID','MALE_RID'), all.x=1)
	ptc[, GD_DI:= cut(GD, breaks=c(-Inf,0.035,0.07,Inf), labels=c('close','unclear','distant'))]
	ptc			<- subset(ptc, !is.na(GD_DI))
	
	unique(ptc,by=c('FEMALE_RID','MALE_RID'))[, table(!is.na(COUPID))]
	unique(subset(ptc, GD<0.035),by=c('FEMALE_RID','MALE_RID'))[, table(!is.na(COUPID))]	# 4052
	unique(subset(ptc, GD>0.07),by=c('FEMALE_RID','MALE_RID'))[, table(!is.na(COUPID))]		# 2369
	
	#subset(ptc, !is.na())
	
	#
	#	compare to those with sufficient PANGEA data for phyloscanner
	#
	
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	rc		<- copy(rp)
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments_allpairs.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	# load pairwise probabilities
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmp_allpairs_posteriors.rda')
	
	#	select: distant, close, unclear from phyloscanner distances
	rpd		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_DI') & TYPE=='distant' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>0.5)		
	rpd[, PHSC_DI:='distant']
	tmp		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_DI') & TYPE=='close' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>0.5)		
	tmp[, PHSC_DI:='close']
	rpd		<- rbind(rpd, tmp)
	tmp		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_DI') & (TYPE=='distant' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)<=0.5) | (TYPE=='close' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)<=0.5))		
	tmp[, PHSC_DI:='unclear']
	tmp		<- unique(tmp, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID'))
	rpd		<- rbind(rpd, tmp)
	rpd		<- merge(rpd, unique(subset(rpw, select=c('FEMALE_RID','MALE_RID','FEMALE_PID','MALE_PID','MALE_SANGER_ID','FEMALE_SANGER_ID'))), by=c('FEMALE_SANGER_ID','MALE_SANGER_ID'))	
	rpd		<- unique(subset(rpd, select=c('FEMALE_RID','MALE_RID','FEMALE_PID','MALE_PID','PHSC_DI')))
	
	#	merge distances 
	rpd		<- merge(ptc, rpd, by=c('FEMALE_RID','MALE_RID','FEMALE_PID','MALE_PID'), all.x=1)	
	rpd		<- unique(rpd,by=c('FEMALE_RID','MALE_RID'))
	#	merge communities
	tmp		<- unique(subset(rd, select=c(RID, COMM_NUM, COMM_TYPE)), by='RID')
	setnames(tmp, colnames(tmp), paste0('MALE_',colnames(tmp)))
	rpd		<- merge(rpd, tmp, by='MALE_RID')
	setnames(tmp, colnames(tmp), gsub('MALE_','FEMALE_',colnames(tmp)))
	rpd		<- merge(rpd, tmp, by='FEMALE_RID')
	
	subset(rpd, !is.na(COUP_SC))[, table(MALE_COMM_TYPE, FEMALE_COMM_TYPE)]
	#MALE_COMM_TYPE agrarian fisherfolk trading
    #	agrarian         66          0       0
    #	fisherfolk        3        231       0
    #	trading           0          0       6

	subset(rpd, !is.na(COUP_SC) & !is.na(PHSC_DI))[, table(MALE_COMM_TYPE, FEMALE_COMM_TYPE)]
	#MALE_COMM_TYPE agrarian fisherfolk trading
    #agrarian         54          0       0
    #fisherfolk        3        155       0
    #trading           0          0       4
	
	rpd[, list(INCOUPLE=any(!is.na(COUP_SC))), by='FEMALE_RID'][,table(INCOUPLE)]
	#	FALSE  TRUE 
  	#	123   301 
	rpd[, list(INCOUPLE=any(!is.na(COUP_SC))), by='MALE_RID'][,table(INCOUPLE)]
	#	FALSE  TRUE 
	#	   63   298
	rpd[, list(INCOUPLE=any(!is.na(COUP_SC))), by='FEMALE_RID'][,table(INCOUPLE)]
	#	FALSE  TRUE 
	#	123   301 
	rpd[, list(INCOUPLE=any(!is.na(COUP_SC))), by='MALE_RID'][,table(INCOUPLE)]
	#	FALSE  TRUE 
	#	   63   298

	rpd[, table(GD_DI, PHSC_DI, useNA='if')]
	#         PHSC_DI
	#GD_DI     	  close distant unclear <NA>
  	#close     		194       2      17 3830
  	#unclear    	 18      31      57 1969
  	#distant     	  5    1293      87  982
	subset(rpd, !is.na(COUP_SC))[, table(GD_DI, PHSC_DI, useNA='if')]
	#			PHSC_DI
	#GD_DI     close distant unclear   <NA>
	#close       122       1       6   49
	#unclear       8       2       1   21
	#distant       0      72       4   20

	#
	#	double check NA's
	#
	tmp		<- phsc.combine.phyloscanner.output('~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_161213_couples_w270_d50_p001_rerun')
	dwin	<- tmp$dwin
	tmp		<- subset(pty.runs, select=c('TAXA','FILE_ID','PTY_RUN'))
	setnames(tmp, c('TAXA','FILE_ID'), c('TAXA1','ID1'))
	dwin	<- merge(dwin, unique(tmp), by=c('ID1','PTY_RUN'))	#something does not quite match here
	setnames(tmp, c('TAXA1','ID1'), c('TAXA2','ID2') )
	dwin	<- merge(dwin, unique(tmp), by=c('ID2','PTY_RUN'))	#something does not quite match here
	
	rpch	<- subset(rpd, is.na(PHSC_DI), select=c('FEMALE_RID','MALE_RID','TAXA1','TAXA2','PTY_RUN'))
	tmp		<- copy(rpch)
	setnames(tmp, c('FEMALE_RID','MALE_RID','TAXA1','TAXA2','PTY_RUN'), c('MALE_RID','FEMALE_RID','TAXA2','TAXA1','PTY_RUN'))
	rpch	<- rbind(rpch, tmp)	
	rpch	<- merge(rpch, dwin, by=c('TAXA1','TAXA2','PTY_RUN'))
	
	tmp		<- rpch[, list(MXMIN_R= max(pmin(ID1_R, ID2_R))), by=c('FEMALE_RID','MALE_RID','TAXA1','TAXA2','PTY_RUN')]
	
	unique(merge(rpch, dwin, by=c('TAXA1','TAXA2','PTY_RUN')), by=c('TAXA1','TAXA2','PTY_RUN'))
	unique(rpch, by=c('TAXA1','TAXA2','PTY_RUN'))
	

	subset(rpd, !is.na(COUP_SC) & is.na(PHSC_DI) & GD_DI=='close')

	
	
}

RakaiAll.analyze.pairs.170120.comparedirection<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]	
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	rc		<- copy(rp)
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments_allpairs.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	# load pairwise probabilities
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmp_allpairs_posteriors.rda')
	
	#	make pair data table
	rp		<- copy(rpw)
	set(rp, NULL, c('DIR','FILE','PTY_RUN','RUN','W_FROM','W_TO','TYPE_RAW','TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3','PATRISTIC_DISTANCE','CONTIGUOUS','PATHS_12','PATHS_21','MALE_SANGER_ID_L','MALE_SANGER_ID_R','FEMALE_SANGER_ID_L','FEMALE_SANGER_ID_R'), NULL)
	rp		<- unique(rp)
	#	make COUPID
	rp[, COUPID:= paste0(MALE_RID,':',FEMALE_RID)]	
	#	add PAIR_TYPE
	tmp		<- unique(subset(rc, select=c(COUPID, MALE_HH_NUM, FEMALE_HH_NUM, COUP_SC, PAIR_TYPE)))	
	setnames(tmp, 'COUP_SC', 'COUP_TYPE')
	rp	<- merge(rp, tmp, by=c('COUPID','MALE_HH_NUM','FEMALE_HH_NUM'),all.x=1)
	set(rp, rp[, which(is.na(PAIR_TYPE))], 'PAIR_TYPE', 'not registered as couple')
	
	#	add RCCS IDs COUPLE_TYPE PAIR_TYPE COUPID
	tmp		<- subset(rp, select=c(FEMALE_SANGER_ID, MALE_SANGER_ID, MALE_RID, FEMALE_RID, COUPID, COUP_TYPE, PAIR_TYPE))
	rplkl	<- merge(tmp, rplkl, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
	
	#
	#	compare distance / distance + topology
	#	results: 	with TYPE_PAIR_TODI3 pairs in close transmission chains end up in unlinked 
	#				contiguous cannot be interpreted as 'A->C->B' it can also be 'A->B/C' or 'A->B and A->C'
	confint	<- 0.5
	unique(subset(rplkl, GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='no intermediate\n and close' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confint), by='COUPID')
	#	193 -- not sure why I get 191 below???
	unique(subset(rplkl, GROUP%in%c('TYPE_PAIR_DI') & TYPE=='close' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confint), by='COUPID')
	
	rpd		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_DI') & TYPE=='close')
	tmp		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='no intermediate\n and close')
	rpd		<- rbind(rpd, tmp)
	rpd[, SELECT_P50:= as.character(factor(pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confint, levels=c(TRUE,FALSE), labels=c('Y','N')))]
	set(rpd, NULL, c('COUPID','COUP_TYPE','RUN','PAIR_TYPE'), NULL)
	rpd		<- merge(rp, rpd, by=c('FEMALE_SANGER_ID','MALE_SANGER_ID','MALE_RID','FEMALE_RID'))
	rpd		<- dcast.data.table(rpd, COUPID+COUP_TYPE+LABEL+PTY_RUN+FEMALE_SANGER_ID+MALE_SANGER_ID+MALE_RID+FEMALE_RID~GROUP, value.var='SELECT_P50')
	
	unique(subset(rpd, !is.na(COUP_TYPE)), by='COUPID')[, table(TYPE_PAIR_DI, TYPE_PAIR_TODI3)]
	#            TYPE_PAIR_TODI3
	#TYPE_PAIR_DI   N   Y
    #       N  		85   0
    #       Y   	5    126
	unique(rpd, by='COUPID')[, table(TYPE_PAIR_DI, TYPE_PAIR_TODI3)]
	#            TYPE_PAIR_TODI3
	#TYPE_PAIR_DI    N    Y
    #       N 	  1519    0
    #       Y       32    191

	tmp		<- subset(rpd, TYPE_PAIR_DI=='Y' & TYPE_PAIR_TODI3=='N')
	unique(tmp, by=c('MALE_RID','FEMALE_RID'))
}
	
RakaiAll.analyze.pairs.170227.comparetoprevious<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	#	get epi info
	tmp		<- RakaiCirc.epi.get.info.170208()
	rh		<- tmp$rh
	rd		<- tmp$rd
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision/RCCS_SeqInfo_160816.rda')		
	rs		<- subset(rs, !is.na(VISIT))	
	
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	rc		<- copy(rp)
	# load transmission summaries	
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments_allpairs.rda')	
	tmp		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/170227/RCCS_170227_w270_trmw_assignments_allpairs.rda')
	set(tmp, NULL, setdiff(colnames(tmp),colnames(rpw)), NULL)
	rpw		<- rbind(rpw, tmp, fill=TRUE)		
	run		<- 'RCCS_170227_w270_dxxx'
	dir		<- '/Users/Oliver/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/170227'		
	# load pairwise probabilities
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/170227/RCCS_170227_w270_trmp_allpairs_posteriors.rda')
	set(rplkl, NULL, 'MALE_SANGER_ID', rplkl[, as.character(MALE_SANGER_ID)])
	set(rplkl, NULL, 'FEMALE_SANGER_ID', rplkl[, as.character(FEMALE_SANGER_ID)])
	
	#	make pair data table
	rp		<- copy(rpw)
	set(rp, NULL, c('DIR','FILE','PTY_RUN','RUN','W_FROM','W_TO','TYPE_RAW','TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3','PATRISTIC_DISTANCE','CONTIGUOUS','PATHS_12','PATHS_21','MALE_SANGER_ID_L','MALE_SANGER_ID_R','FEMALE_SANGER_ID_L','FEMALE_SANGER_ID_R'), NULL)
	rp		<- unique(rp)
	#	make COUPID
	rp[, COUPID:= paste0(MALE_RID,':',FEMALE_RID)]	
	#	add PAIR_TYPE
	tmp		<- unique(subset(rc, select=c(COUPID, MALE_HH_NUM, FEMALE_HH_NUM, COUP_SC, PAIR_TYPE)))	
	setnames(tmp, 'COUP_SC', 'COUP_TYPE')
	rp	<- merge(rp, tmp, by=c('COUPID','MALE_HH_NUM','FEMALE_HH_NUM'),all.x=1)
	set(rp, rp[, which(is.na(PAIR_TYPE))], 'PAIR_TYPE', 'not registered as couple')
	#	add RCCS IDs COUPLE_TYPE PAIR_TYPE COUPID
	tmp		<- subset(rp, select=c(FEMALE_SANGER_ID, MALE_SANGER_ID, MALE_RID, FEMALE_RID, COUPID, COUP_TYPE, PAIR_TYPE))
	tmp		<- unique(tmp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rplkl	<- merge(tmp, rplkl, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
	#	select first couples for whom transmission cannot be excluded 
	#	no intermediate\n and close > 50%
	confidence.cut	<- 0.5
	rpd		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='no intermediate\n and close' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confidence.cut)
	#
	#	plot phylogenies for all inconsistent pairs in either of the two runs
	#
	tmp		<- dcast.data.table(rpd, LABEL_SH~RUN, value.var='KEFF')
	tmp		<- subset(tmp, is.na(RCCS_170227_w270_d50_st20_mr20_mt1_cl2_d5) | is.na(RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5))
	dch		<- merge(tmp, unique(subset(rpd, select=c('LABEL_SH','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','MALE_RID','FEMALE_RID','COUPID'))), by='LABEL_SH')	
	#	check RCCS_170227_w270_d50_st20_mr20_mt1_cl2_d5
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/170227/RCCS_170227_w270_d50_st20_mr20_mt1_cl2_d5_phscout.rda')
	outfile	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/170227/noTP_RCCS_170227_w270_d50_st20_mr20_mt1_cl2_d5_'
	tmp		<- subset(dch, is.na(RCCS_170227_w270_d50_st20_mr20_mt1_cl2_d5))
	invisible(sapply(seq_len(nrow(tmp)), function(ii)
				{	
					#ii<- 14
					ids			<- c(tmp[ii, MALE_SANGER_ID],tmp[ii, FEMALE_SANGER_ID])
					pty.run		<- tmp[ii, PTY_RUN]
					dfs			<- subset(dtrees, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
					dfs[, TITLE:= dfs[, paste('couple', tmp[ii,COUPID], '\nmale ', ids[1],'\nfemale ',ids[2],'\nrun ', pty.run, '\nwindow ', W_FROM,'-', W_TO,sep='')]]			
					plot.file	<- paste0(outfile, pty.run,'_',ids[1],'_', ids[2],'.pdf')
					invisible(phsc.plot.selected.individuals(phs, dfs, ids, plot.cols=c('red','blue'), group.redo=TRUE, plot.file=plot.file, pdf.h=150, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40))					
				}))	
	#	are the missing couples in our dataset before subset confidence.cut is taken?
	#	yes, none of them removed with new blacklister
	subset(merge(subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID)), rplkl, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='no intermediate\n and close' & RUN=='RCCS_170227_w270_d50_st20_mr20_mt1_cl2_d5')

	tmp		<- subset(dch, is.na(RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5))
	subset(merge(subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID)), rplkl, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')), GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='no intermediate\n and close' & RUN=='RCCS_170227_w270_d50_st20_mr20_mt1_cl2_d5')
	
	tmp		<- dcast.data.table(rpd, LABEL_SH~RUN, value.var='KEFF')
	tmp[, table(is.na(RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5), is.na(RCCS_170227_w270_d50_st20_mr20_mt1_cl2_d5))]
	subset(tmp, is.na(RCCS_170227_w270_d50_st20_mr20_mt1_cl2_d5))	
}

	
RakaiAll.analyze.pairs.170120.direction<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	rc		<- copy(rp)
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments_allpairs.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]		
	# load pairwise probabilities
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmp_allpairs_posteriors.rda')
	
	#	make pair data table
	rp		<- copy(rpw)
	set(rp, NULL, c('DIR','FILE','PTY_RUN','RUN','W_FROM','W_TO','TYPE_RAW','TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3','PATRISTIC_DISTANCE','CONTIGUOUS','PATHS_12','PATHS_21','MALE_SANGER_ID_L','MALE_SANGER_ID_R','FEMALE_SANGER_ID_L','FEMALE_SANGER_ID_R'), NULL)
	rp		<- unique(rp)
	#	make COUPID
	rp[, COUPID:= paste0(MALE_RID,':',FEMALE_RID)]	
	#	add PAIR_TYPE
	tmp		<- unique(subset(rc, select=c(COUPID, MALE_HH_NUM, FEMALE_HH_NUM, COUP_SC, PAIR_TYPE)))	
	setnames(tmp, 'COUP_SC', 'COUP_TYPE')
	rp	<- merge(rp, tmp, by=c('COUPID','MALE_HH_NUM','FEMALE_HH_NUM'),all.x=1)
	set(rp, rp[, which(is.na(PAIR_TYPE))], 'PAIR_TYPE', 'not registered as couple')

	#
	#	basic info on selection
	#
	if(0)
	{
		tmp		<- unique(rp, by='COUPID')
		tmp[, c(length(unique(MALE_RID)),length(unique(FEMALE_RID)))]
		subset(tmp, PAIR_TYPE=='cohabiting couple')[, c(length(unique(MALE_RID)),length(unique(FEMALE_RID)))]
		subset(tmp, !is.na(COUP_TYPE))[, c(length(unique(MALE_RID)),length(unique(FEMALE_RID)),length(unique(COUPID)))]		
	}

	#	add RCCS IDs COUPLE_TYPE PAIR_TYPE COUPID
	tmp		<- subset(rp, select=c(FEMALE_SANGER_ID, MALE_SANGER_ID, MALE_RID, FEMALE_RID, COUPID, COUP_TYPE, PAIR_TYPE))
	rplkl	<- merge(tmp, rplkl, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
	#	select first couples for whom transmission cannot be excluded 
	#	no intermediate\n and close > 50%
	confidence.cut	<- 0.5
	rpd		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='no intermediate\n and close' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confidence.cut)		
	tmp		<- subset(rpd, select=c(PTY_RUN, FEMALE_SANGER_ID, MALE_SANGER_ID))
	tmp		<- merge(tmp, rplkl, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	rmf		<- subset(tmp, GROUP%in%c('TYPE_DIR_TODI3') & TYPE=='mf' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confidence.cut)
	rfm		<- subset(tmp, GROUP%in%c('TYPE_DIR_TODI3') & TYPE=='fm' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confidence.cut)
	rex		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='with intermediate\nor distant' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confidence.cut)
	
	cat('\ncouples with phyloscanner assessment, n=',nrow(unique(rplkl, by='COUPID')))
	cat('\ncouples for whom transmission cannot be excluded, n=',nrow(unique(rpd, by='COUPID')))
	cat('\ncouples with evidence M->F, n=',nrow(unique(rmf, by='COUPID')))
	cat('\ncouples with evidence F->M, n=',nrow(unique(rfm, by='COUPID')))
	cat('\ncouples for whom transmission can be excluded, n=',nrow(unique(rex, by='COUPID')))
	#	couples with phyloscanner assessment					1742	
	#	for whom transmission can be excluded 					1414	(81.2%)
	#	ambiguous												135		(7.7%)
	#	no intermediate and close								193		(11%)	
	#	of those we cannot determine likely direction for		53		(31%)
	#	with evidence for M->F 									87		(60%)
	#	with evidence for F->M									53 		(39%)
	
	unique(rpd, by='COUPID')[, table(PAIR_TYPE)]
	#	   not always cohabiting 	not registered as couple        stable cohabiting 
    #                          3                       	  71        119  

	rmf		<- merge(unique(subset(rmf, select=COUPID)), unique(rp, by='COUPID'),by='COUPID')
	rmf[, PHSC_DIR:='m->f']
	rfm		<- merge(unique(subset(rfm, select=COUPID)), unique(rp, by='COUPID'),by='COUPID')
	rfm[, PHSC_DIR:='f->m']
	rtr		<- rbind(rmf, rfm)
	
	rtr[, AGEDIFF:= MALE_BIRTHDATE-FEMALE_BIRTHDATE]
	rtr[, AVGAGE:= (MALE_BIRTHDATE+FEMALE_BIRTHDATE)/2]
	
	rtr2	<- copy(rmf)
	setnames(rtr2,colnames(rtr2),gsub('FEMALE','REC',colnames(rtr2)))
	setnames(rtr2,colnames(rtr2),gsub('MALE','TR',colnames(rtr2)))
	tmp		<- copy(rfm)
	setnames(tmp,colnames(tmp),gsub('FEMALE','TR',colnames(tmp)))
	setnames(tmp,colnames(tmp),gsub('MALE','REC',colnames(tmp)))
	rtr2	<- rbind(rtr2,tmp)
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
	
	rtr2[, table(TR_COMM_TYPE,REC_COMM_TYPE)]
	#				 REC_COMM_TYPE
	#			     agrarian fisherfolk trading
	#TR_COMM_TYPE 
	#agrarian         36          3       2
  	#fisherfolk        2         95       0
  	#trading           0          0       2
	tmp		<- rtr2[,list(N=length(unique(COUPID))), by=c('TR_COMM_NUM','REC_COMM_NUM')]
	ggplot(tmp, aes(x=factor(REC_COMM_NUM),y=factor(TR_COMM_NUM))) + 
			geom_point(aes(size=N), colour='grey80') +
			geom_text(aes(label=N), size=3, colour='black') +
			scale_size(range = c(1, 20)) +
			theme_bw() + 
			labs(x='\nlocation likely recipient',y='location likely transmitter\n') +
			guides(size='none')
	ggsave(file=file.path(dir, paste(run,'-phsc-directionpairs_direction-numbers-commnum.pdf',sep='')), w=6, h=6)
	
	#
	#	did any transmitter start ART before the recipient was diagnosed?
	#
	subset(rtr2, TR_ARVSTARTDATE<REC_FIRSTPOSDATE)	
	#	F026858:J104288 --> 

	subset(rtr2, (TR_FIRSTPOSDATE+.5)<REC_FIRSTPOSDATE)	

	#
	#	Does the primary occupation differ between transmitters / recipients? 
	#
	ggplot(subset(rtr, FEMALE_COMM_TYPE!='trading'), aes(x=PHSC_DIR, fill=FEMALE_OCCUP1)) + 
			geom_bar(position='fill') + 
			facet_grid(~FEMALE_COMM_TYPE) +
			scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), expand=c(0,0)) +
			scale_fill_brewer(palette='Set3') +
			theme_bw() + theme() +
			labs(x='\ninferred direction of transmission', y='female partner\n', fill='occupation\nfemale partner\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-directionpairs_occupationfemales.pdf',sep='')), w=6, h=8)

	colourCount	<- length(unique(rtr$MALE_OCCUP1))
	getPalette	<- colorRampPalette(brewer.pal(11, "Set3"))	
	ggplot(subset(rtr, MALE_COMM_TYPE!='trading'), aes(x=PHSC_DIR, fill=MALE_OCCUP1)) + 
			geom_bar(position='fill') + 
			facet_grid(~MALE_COMM_TYPE) +
			scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), expand=c(0,0)) +
			scale_fill_manual(values=getPalette(colourCount)) +
			#scale_fill_brewer(palette='Set3') +
			#scale_fill_distiller(palette='Set3') +
			theme_bw() + theme() +
			labs(x='\ninferred direction of transmission', y='male partner\n', fill='occupation\nmale partner\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-directionpairs_occupationmales.pdf',sep='')), w=8, h=8)
	
	subset(rtr, MALE_COMM_TYPE=='fisherfolk' & PHSC_DIR=='m->f')[, table(MALE_OCCUP1)]
	binconf(41,52)
	#
	#	housework significantly associated with being a female recipient?
	#
	tmp	<- subset(rtr, FEMALE_COMM_TYPE=='fisherfolk')
	tmp[, HOUSEWORK:= factor(grepl('Housework',FEMALE_OCCUP1), levels=c(TRUE,FALSE), labels=c('Y','N'))]
	chisq_test(factor(PHSC_DIR) ~ HOUSEWORK, data=tmp, distribution="exact")
	#	chi-squared = 1.2911, p-value = 0.2774

	
	#
	#	do recipients / transmitters report other partners?
	#	

	ggplot(subset(rtr, PAIR_TYPE=='stable cohabiting' & FEMALE_COMM_TYPE!='trading'), aes(x=PHSC_DIR, fill=FEMALE_SEXC)) + 
			geom_bar(position='fill') + 
			facet_grid(~FEMALE_COMM_TYPE) +
			scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), expand=c(0,0)) +
			scale_fill_brewer(palette='Set2') +
			theme_bw() + theme() +
			labs(x='\ninferred direction of transmission', y='female partner\n', fill='self-reported relationships\nfemale partner\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-directionpairs_relationshiptypefemales.pdf',sep='')), w=6, h=8)

	ggplot(subset(rtr, PAIR_TYPE=='stable cohabiting' & MALE_COMM_TYPE!='trading'), aes(x=PHSC_DIR, fill=MALE_SEXC)) + 
			geom_bar(position='fill') + 
			facet_grid(~MALE_COMM_TYPE) +
			scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), expand=c(0,0)) +
			scale_fill_brewer(palette='Set2') +
			theme_bw() + theme() +
			labs(x='\ninferred direction of transmission', y='male partner\n', fill='self-reported relationships\nmale partner\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-directionpairs_relationshiptypemales.pdf',sep='')), w=6, h=8)
	
	#
	#	is there a difference in male->female transmission by community type?
	#	results: no	- but it seems there could be a difference between stable/casual couples
	#
	subset(rtr, MALE_COMM_TYPE==FEMALE_COMM_TYPE)[, {
				z	<- binconf( length(which(PHSC_DIR=='m->f')), length(PHSC_DIR) )				
				list(K=length(which(PHSC_DIR=='m->f')), N=length(PHSC_DIR), P=z[1], QL=z[2], QU=z[3])
			}, by='MALE_COMM_TYPE']	
	#	MALE_COMM_TYPE  K  N         P        QL        QU
	#1:       agrarian 24 36 0.6666667 0.5033400 0.7978538
	#2:     fisherfolk 58 95 0.6105263 0.5100024 0.7024590
	#3:        trading  2  2 1.0000000 0.3423802 1.0000000
	
	#
	#	is there a difference in male->female transmission by couple type?
	#	results: no	
	#			 but I still find this interesting!
	rtr[, {
				z	<- binconf( length(which(PHSC_DIR=='m->f')), length(PHSC_DIR) )				
				list(K=length(which(PHSC_DIR=='m->f')), N=length(PHSC_DIR), P=z[1], QL=z[2], QU=z[3])
			}, by='PAIR_TYPE']	
	#			PAIR_TYPE  K  N         P        QL        QU
	#1: stable cohabiting 			57 85 0.6705882 0.56520191 0.7612223
	#2: not registered as couple 	29 53 0.5471698 0.41453975 0.6734242
	#3:    not always cohabiting  	 1  2 0.5000000 0.02564665 0.9743534	
	chisq_test(factor(PHSC_DIR) ~ factor(PAIR_TYPE), data=rtr, distribution="exact")
	#	chi-squared = 1.8666, p-value = 0.1953


	#	are transmitters younger in fisherfolk sites?
	#	results: no
	tmp		<- subset(rtr2, TR_COMM_TYPE!='trading')
	ggplot(tmp, aes(x=TR_COMM_TYPE, y=TR_BIRTHDATE)) + geom_boxplot()
	independence_test(TR_BIRTHDATE~factor(TR_COMM_TYPE), data=tmp, distribution = "exact")
	#Z = -1.4964, p-value = 0.1351
	summary(rq(TR_BIRTHDATE~TR_COMM_TYPE, tau=.5, data=tmp, method='fn'), se='nid')

	#
	#	is there a difference in age gap between male->female transmission / female->male transmission ?
	#	results: no
	#
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



	#
	#	are couples overall younger in male->female transmission ?
	#	results: no
	ggplot(rtr, aes(x=PHSC_DIR, y=AVGAGE)) + geom_boxplot()
	tmp		<- subset(rtr, !is.na(AVGAGE), select=c(PHSC_DIR,AVGAGE))
	summary(rq(AVGAGE~PHSC_DIR, tau=.5, data=tmp, method='fn'), se='nid')
	
	
	#
	#	are female transmitters younger than female recipients ?
	#	results: yes but not significant if all couples are considered
	#			 yes and signficant if fairly old couples are excluded (female age < 40)
	#			 yes, especially in fisherfolk sites!
	ggplot(subset(rtr, FEMALE_COMM_TYPE!='trading'), aes(x=PHSC_DIR, y=FEMALE_BIRTHDATE)) + geom_boxplot() + facet_grid(~FEMALE_COMM_TYPE)
	ggplot(subset(rtr,FEMALE_BIRTHDATE>1975), aes(x=PHSC_DIR, y=FEMALE_BIRTHDATE)) + geom_boxplot()	
	#	the female birthdate distribution among m->f is asymmetric 
	#	the female birthdate distribution has outliers in f->m
	#	try median
	tmp		<- subset(rtr, !is.na(FEMALE_BIRTHDATE), select=c(PHSC_DIR,FEMALE_BIRTHDATE))
	summary(rq(FEMALE_BIRTHDATE~PHSC_DIR, tau=.5, data=tmp, method='fn'), se='nid')
	#	             Value      Std. Error t value    Pr(>|t|)  
	#             Value      Std. Error t value    Pr(>|t|)  
	#(Intercept)  1986.37500    0.84098 2361.97122    0.00000
	#PHSC_DIRm->f   -1.74100    1.20906   -1.43996    0.15228	
	tmp		<- subset(rtr, !is.na(FEMALE_BIRTHDATE) & FEMALE_BIRTHDATE>1975, select=c(PHSC_DIR,FEMALE_BIRTHDATE))
	summary(rq(FEMALE_BIRTHDATE~PHSC_DIR, tau=.5, data=tmp, method='fn'), se='nid')
	#             	Value      Std. Error t value    Pr(>|t|)  
	#(Intercept)  1986.91500    0.70542 2816.64130    0.00000
	#PHSC_DIRm->f   -2.19252    1.03707   -2.11415    0.03652
	
	
	#
	#	are male transmitters younger than male recipients ?
	#	results: overall, male transmitter tend to be older than male recipients but this is not significant 
	#			 in agrarian, male transmitters tend to be younger than male recipients, but this is not significant
	tmp		<- subset(rtr, MALE_COMM_TYPE!='trading' & MALE_COMM_TYPE==FEMALE_COMM_TYPE)
	tmp		<- melt(tmp, measure.vars=c('MALE_BIRTHDATE','FEMALE_BIRTHDATE'))
	set(tmp, NULL, 'variable', tmp[, as.character(factor(variable, levels=c('MALE_BIRTHDATE','FEMALE_BIRTHDATE'), labels=c('male partner','female partner')))])
	ggplot(tmp, aes(x=PHSC_DIR, y=value)) + 
			geom_boxplot() + 
			facet_grid(~MALE_COMM_TYPE+variable) +
			theme_bw() + labs(x='\nestimated direction of transmission', y='birth date\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-directionpairs_direction-age-commtype.pdf',sep='')), w=6, h=6)
	
	
	tmp		<- subset(rtr, !is.na(MALE_BIRTHDATE), select=c(PHSC_DIR,MALE_COMM_TYPE,MALE_BIRTHDATE))
	summary(rq(MALE_BIRTHDATE~PHSC_DIR, tau=.5, data=tmp, method='fn'), se='nid')
	#				Value      Std. Error t value    Pr(>|t|)  
	#(Intercept)  1981.35600    1.61670 1225.55813    0.00000
	#PHSC_DIRm->f   -0.60556    1.84886   -0.32753    0.74381
	summary(rq(MALE_BIRTHDATE~PHSC_DIR, tau=.5, data=subset(tmp,MALE_COMM_TYPE=='agrarian'), method='fn'), se='nid')
	#				Value      Std. Error t value    Pr(>|t|)  
	#(Intercept)  1977.30100    3.97448  497.49959    0.00000
	#PHSC_DIRm->f    3.03691    4.39768    0.69057    0.49414

	tmp		<- subset(rtr, !is.na(MALE_BIRTHDATE) & FEMALE_BIRTHDATE>1975 & MALE_BIRTHDATE>1975, select=c(PHSC_DIR,MALE_BIRTHDATE))
	summary(rq(MALE_BIRTHDATE~PHSC_DIR, tau=.5, data=tmp, method='fn'), se='nid')
	#	

}
	
RakaiAll.addposteriors.pairs.170120<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments_allpairs.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	#
	rpw[, table(RUN, useNA='if')]
	#	define plotting order: largest number of trm assignments	
	tmp		<- rpw[, list( WIN_TR=length(which(grepl('close|anc',TYPE_DIR_TODI7x3))) ), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	tmp		<- tmp[order(-WIN_TR),]
	tmp[, PLOT_ID:=seq_len(nrow(tmp))]	
	#	make pair data table
	rp		<- copy(rpw)
	set(rp, NULL, c('PTY_RUN','RUN','W_FROM','W_TO','TYPE_RAW','TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3','PATRISTIC_DISTANCE','CONTIGUOUS','PATHS_12','PATHS_21','MALE_SANGER_ID_L','MALE_SANGER_ID_R','FEMALE_SANGER_ID_L','FEMALE_SANGER_ID_R'), NULL)
	rp		<- unique(rp)
	#	make COUPID
	rp[, COUPID:= paste0(MALE_RID,':',FEMALE_RID)]
	#	define label
	tmp		<- merge(tmp, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
	setkey(tmp, PLOT_ID)
	tmp[, LABEL_SH:= factor(PLOT_ID, levels=PLOT_ID, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair ', COUPID,' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL, LABEL_SH))
	rpw		<- merge(tmp, rpw, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	#	define min/max reads
	tmp		<- rpw[, list(	ID_R_MIN=min(MALE_SANGER_ID_R, FEMALE_SANGER_ID_R),
							ID_R_MAX=max(MALE_SANGER_ID_R, FEMALE_SANGER_ID_R)), by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM')]
	rpw		<- merge(rpw, tmp, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM'))
	#
	#	define topo/distance types of phylogenetic relationship between pairs
	#
	#	given the bimodal distribution of patristic distances from phyloscanner
	#	I ususally consider 3 distance states: close, neither close nor distant, distant
	#
	#	TYPE_DIR_TODI7x3		7 topology states (this is the "basic topological characterization" on Slack that I derive from Matthew's output) 
	#							chain_12_nointermediate, chain_12_withintermediate, chain_21_nointermediate, chain_21_withintermediate, intermingled_nointermediate, intermingled_withintermediate, other  
	#							times 
	#							3 distance states
	#
	#	TYPE_PAIR_TODI3x3		3 topology states
	#							pair (chain_XX_nointermediate+intermingled_nointermediate), withintermediate (XX_withintermediate), other
	#							times 
	#							3 distance states
	#
	#	TYPE_PAIR_TODI2x2		2 topology states
	#							ancestral/intermingled, other
	#							times 
	#							2 distance states (close, not close)
	#
	#	TYPE_PAIR_TODI3			non-factorial design that combines distant or withintermediate, and close or ancestral/intermingled
	#							the idea here is that 
	#							evidence for extra-couple transmission comes from large patristic distance OR intermediates	
	#							evidence for extra-couple transmission comes from ancestral/intermingled OR short patristic distance	
	#
	#	TYPE_PAIR_DI			3 distance states
	#
	#
	#	TYPE_DIR_TO3			3 topological direction states
	#							chain_fm	(chain_fm_nointermediate regardless of distance); chain_mf	(chain_mf_nointermediate regardless of distance); other
	#							Rest: intermingled -> NA.
	#
	#	TYPE_DIR_TODI3			3 topological direction states if close.  
	#							chain_fm	(chain_fm_nointermediate regardless of distance); chain_mf	(chain_mf_nointermediate regardless of distance); other
	#							Rest: intermingled / distant / intermediate distance -> NA.
	#
	#	TYPE_CHAIN_TODI2		2 topological states if close.  
	#							chain	(xxx_withintermediate_close); pair	(xxx_nointermediate_close, other)
	#							Rest: distant / intermediate distance -> NA.
	rpw[, TYPE_PAIR_TO:= 'other']
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair_close','pair','pair_distant'))], 'TYPE_PAIR_TO', 'ancestral/\nintermingled')
	rpw[, TYPE_PAIR_TODI3:= 'with intermediate\nor distant']
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair_close','other_close'))], 'TYPE_PAIR_TODI3', 'no intermediate\n and close')
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair','other'))], 'TYPE_PAIR_TODI3', 'no intermediate\n but not close')
	rpw[, TYPE_PAIR_DI:= cut(PATRISTIC_DISTANCE, breaks=c(1e-12,0.02,0.05,2), labels=c('close','intermediate\ndistance','distant'))] 	
	rpw[, TYPE_PAIR_TODI2x2:= 'not close other']
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair_close'))], 'TYPE_PAIR_TODI2x2', 'close ancestral/\nintermingled')
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair','pair_distant'))], 'TYPE_PAIR_TODI2x2', 'not close ancestral/\nintermingled')
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('withintermediate_close','other_close'))], 'TYPE_PAIR_TODI2x2', 'close other')	
	rpw[, TYPE_DIR_TO3:= NA_character_]
	set(rpw, rpw[, which(grepl('chain_fm',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TO3', 'fm')
	set(rpw, rpw[, which(grepl('chain_mf',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TO3', 'mf')
	set(rpw, rpw[, which(grepl('intermingled',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TO3', 'intermingled')		
	rpw[, TYPE_DIR_TODI3:= NA_character_]
	set(rpw, rpw[, which(grepl('chain_fm',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'fm')
	set(rpw, rpw[, which(grepl('chain_mf',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'mf')
	set(rpw, rpw[, which(grepl('intermingled',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'intermingled')			
	rpw[, TYPE_CHAIN_TODI2:= NA_character_]
	set(rpw, rpw[, which(grepl('withintermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_CHAIN_TODI2', 'chain')
	set(rpw, rpw[, which(grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_CHAIN_TODI2', 'pair')
	set(rpw, rpw[, which(grepl('other',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_CHAIN_TODI2', 'pair')
	set(rpw, rpw[, which(grepl('intermingled',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'intermingled')			
	set(rpw, NULL, 'TYPE_DIR_TODI7x3', rpw[, gsub('intermediate',' intermediate',gsub('intermediate_','intermediate\n',gsub('intermingled_','intermingled\n',gsub('(chain_[fm][mf])_','\\1\n',TYPE_DIR_TODI7x3))))])
	#	define effectively independent number of windows
	#	select only W_FROM that are > W_TO	
	tmp		<- rpw[, {
				z		<- 1
				repeat
				{
					zz		<- which(W_FROM>max(W_TO[z]))[1]
					if(length(zz)==0 | is.na(zz))	break
					z		<- c(z, zz)			
				}
				list(W_FROM=W_FROM[z], W_TO=W_TO[z])
			}, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN')]
	tmp[, OVERLAP:= 0L]
	rpw		<- merge(rpw, tmp, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','W_FROM','W_TO','RUN'), all.x=1 )
	set(rpw, rpw[, which(is.na(OVERLAP))], 'OVERLAP', 1L)		
	rpw		<- melt(rpw, measure.vars=c('TYPE_PAIR_TO','TYPE_CHAIN_TODI2','TYPE_PAIR_TODI3x3','TYPE_PAIR_TODI3','TYPE_PAIR_DI','TYPE_PAIR_TODI2x2','TYPE_DIR_TODI7x3','TYPE_DIR_TO3','TYPE_DIR_TODI3'), variable.name='GROUP', value.name='TYPE')
	set(rpw, NULL, 'TYPE', rpw[,gsub('_',' ', TYPE)])	
	rpw		<- subset(rpw, !is.na(TYPE))
	#
	#	for each pair count windows by TYPE_PAIR_DI (k)
	#
	rplkl	<- rpw[, list(K=length(W_FROM)), by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','LABEL','LABEL_SH','GROUP','TYPE')]	
	rplkl	<- dcast.data.table(rplkl, RUN+PTY_RUN+FEMALE_SANGER_ID+MALE_SANGER_ID+LABEL+LABEL_SH~GROUP+TYPE, value.var='K')
	for(x in setdiff(colnames(rplkl),c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','LABEL','LABEL_SH','GROUP')))
		set(rplkl, which(is.na(rplkl[[x]])), x, 0L)	
	for(x in colnames(rplkl)[grepl('TYPE_PD_MEAN',colnames(rplkl))])
		set(rplkl, which(rplkl[[x]]>0), x, 1L)	
	rplkl	<- melt(rplkl, id.vars=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','LABEL','LABEL_SH'), variable.name='GROUP', value.name='K')
	rplkl[, TYPE:= gsub('.*_([^_]+)$','\\1',GROUP)]
	set(rplkl, NULL, 'GROUP', rplkl[, gsub('(.*)_[^_]+$','\\1',GROUP)])	
	#	for each pair count total windows (n) and effective windows (neff)	
	tmp		<- rpw[, list(N=length(W_FROM), NEFF=sum(1-OVERLAP)), by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN','GROUP')]
	rplkl	<- merge(rplkl, tmp, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN','GROUP'))
	#	define effective number of windows of type
	rplkl[, KEFF:= K/N*NEFF]	
	#
	#	add prior (equal probabilities to all types), posterior, and marginal probabilities
	#
	rplkl[, DIR_PRIOR:= 0.1]
	rplkl[, DIR_PO:= DIR_PRIOR+KEFF]
	tmp		<- rplkl[, {
				alpha	<- DIR_PO
				beta	<- sum(DIR_PO)-DIR_PO				
				list(	TYPE=TYPE, 
						POSTERIOR_ALPHA=alpha, 
						POSTERIOR_BETA=beta)	
			}, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN','GROUP')]
	rplkl	<- merge(rplkl, tmp, by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','GROUP','TYPE'))
	#	save
	save(rplkl, file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmp_allpairs_posteriors.rda')
	
	#
	#	prepare colours
	#	
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
	#
	#	make detailed plots for each pair
	#		topology assignments across genome
	#		distances across genome
	#		number of windows of certain type
	#		estimated posterior probabilities on unknown phylogenetic relationship
	#
	groups		<- c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3','TYPE_PAIR_TODI2x2','TYPE_PAIR_TODI3','TYPE_PAIR_DI')
	group		<- c('TYPE_DIR_TODI7x3')
	
	widths	<- unit(c(4, 6), "null")
	heights	<- unit(c(2, 3.5, 4, 5), "null")
	height	<- 9
	if(group%in%c('TYPE_DIR_TODI7x3'))
	{
		widths	<- unit(c(4, 6), "null")
		heights	<- unit(c(2, 3.5, 4, 15), "null")
		height	<- 17
	}		
	if(group%in%c('TYPE_PAIR_TODI2x2'))
	{
		heights	<- unit(c(2, 3.5, 4, 3.75), "null")
		height	<- 8
	}
	if(group%in%c('TYPE_PAIR_TODI3','TYPE_PAIR_DI'))
	{
		heights	<- unit(c(2, 3.5, 4, 3.5), "null")
		height	<- 7
	}
	pdf(file=file.path(dir, paste(run,'-phsc-relationships_allpairs','_',group,'.pdf',sep='')), w=10, h=height)	
	plot.tmp	<- unique(subset(rplkl, GROUP==group, c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','LABEL')), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	setkey(plot.tmp, LABEL)		
	for(i in seq_len(nrow(plot.tmp)))
	{
		#pty_run	<- 38; id1		<- '16016_1_4'; id2		<- '15105_1_35'
		pty_run	<- plot.tmp[i, PTY_RUN]; id1		<- plot.tmp[i, FEMALE_SANGER_ID]; id2		<- plot.tmp[i, MALE_SANGER_ID]		
		tmp		<- subset(rpw, PTY_RUN==pty_run & GROUP==group & FEMALE_SANGER_ID==id1 & MALE_SANGER_ID==id2)
		p1		<- ggplot(tmp, aes(x=W_FROM)) +			
				geom_bar(aes(y=ID_R_MAX, colour=TYPE), stat='identity', fill='transparent') +
				geom_bar(aes(y=ID_R_MIN, fill=TYPE), stat='identity', colour='transparent') +
				labs(x='', y='number of reads', fill='phylogenetic\nrelationship\n', colour='phylogenetic\nrelationship\n') +
				scale_fill_manual(values=cols.type[[group]]) +
				scale_colour_manual(values=cols.type[[group]]) +
				scale_x_continuous(breaks=seq(0,1e4,500), minor_breaks=seq(0,1e4,100), limits=c(rpw[, min(W_FROM)], rpw[, max(W_FROM)])) +
				scale_y_log10(breaks=c(10,100,1000,1e4,1e5)) +
				theme_bw() + theme(legend.position='left') +			
				guides(fill=FALSE, colour=FALSE)
		p2		<- ggplot(tmp, aes(x=W_FROM, y=PATRISTIC_DISTANCE)) +
				geom_point(size=1) +					
				labs(x='window start\n\n', y='patristic distance') +
				scale_x_continuous(breaks=seq(0,1e4,500), minor_breaks=seq(0,1e4,100), limits=c(rpw[, min(W_FROM)], rpw[, max(W_FROM)])) +
				scale_y_log10(labels=percent, limits=c(0.001, 0.7), expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
				theme_bw() + theme(legend.position='left')
		tmp		<- subset(rplkl, GROUP==group & PTY_RUN==pty_run & FEMALE_SANGER_ID==id1 & MALE_SANGER_ID==id2)
		p3		<- ggplot(tmp, aes(x=TYPE, y=KEFF, fill=TYPE)) + geom_bar(stat='identity') +
				scale_fill_manual(values=cols.type[[group]]) +
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
				scale_fill_manual(values=cols.type[[group]]) +
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

RakaiAll.addposteriors.pairs.170201<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments_allpairs.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	#
	rpw[, table(RUN, useNA='if')]
	#	define plotting order: largest number of trm assignments	
	tmp		<- rpw[, list( WIN_TR=length(which(grepl('close|anc',TYPE_DIR_TODI7x3))) ), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	tmp		<- tmp[order(-WIN_TR),]
	tmp[, PLOT_ID:=seq_len(nrow(tmp))]	
	#	make pair data table
	rp		<- copy(rpw)
	set(rp, NULL, c('PTY_RUN','RUN','W_FROM','W_TO','TYPE_RAW','TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3','PATRISTIC_DISTANCE','CONTIGUOUS','PATHS_12','PATHS_21','MALE_SANGER_ID_L','MALE_SANGER_ID_R','FEMALE_SANGER_ID_L','FEMALE_SANGER_ID_R'), NULL)
	rp		<- unique(rp)
	#	make COUPID
	rp[, COUPID:= paste0(MALE_RID,':',FEMALE_RID)]
	#	define label
	tmp		<- merge(tmp, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
	setkey(tmp, PLOT_ID)
	tmp[, LABEL_SH:= factor(PLOT_ID, levels=PLOT_ID, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair ', COUPID,' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL, LABEL_SH))
	rpw		<- merge(tmp, rpw, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	#	define min/max reads
	tmp		<- rpw[, list(	ID_R_MIN=min(MALE_SANGER_ID_R, FEMALE_SANGER_ID_R),
					ID_R_MAX=max(MALE_SANGER_ID_R, FEMALE_SANGER_ID_R)), by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM')]
	rpw		<- merge(rpw, tmp, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM'))
	#
	#	define topo/distance types of phylogenetic relationship between pairs
	#
	#	given the bimodal distribution of patristic distances from phyloscanner
	#	I ususally consider 3 distance states: close, neither close nor distant, distant
	#
	#	TYPE_DIR_TODI7x3		7 topology states (this is the "basic topological characterization" on Slack that I derive from Matthew's output) 
	#							chain_12_nointermediate, chain_12_withintermediate, chain_21_nointermediate, chain_21_withintermediate, intermingled_nointermediate, intermingled_withintermediate, other  
	#							times 
	#							3 distance states
	#
	#	TYPE_PAIR_TODI3x3		3 topology states
	#							pair (chain_XX_nointermediate+intermingled_nointermediate), withintermediate (XX_withintermediate), other
	#							times 
	#							3 distance states
	#
	#	TYPE_PAIR_TODI2x2		2 topology states
	#							ancestral/intermingled, other
	#							times 
	#							2 distance states (close, not close)
	#
	#	TYPE_PAIR_TODI3			non-factorial design that combines distant or withintermediate, and close or ancestral/intermingled
	#							the idea here is that 
	#							evidence for extra-couple transmission comes from large patristic distance OR intermediates	
	#							evidence for extra-couple transmission comes from ancestral/intermingled OR short patristic distance	
	#
	#	TYPE_PAIR_DI			3 distance states
	#
	#
	#	TYPE_DIR_TO3			3 topological direction states
	#							chain_fm	(chain_fm_nointermediate regardless of distance); chain_mf	(chain_mf_nointermediate regardless of distance); other
	#							Rest: intermingled -> NA.
	#
	#	TYPE_DIR_TODI3			3 topological direction states if close.  
	#							chain_fm	(chain_fm_nointermediate regardless of distance); chain_mf	(chain_mf_nointermediate regardless of distance); other
	#							Rest: intermingled / distant / intermediate distance -> NA.
	#
	#	TYPE_CHAIN_TODI2		2 topological states if close.  
	#							chain	(xxx_withintermediate_close); pair	(xxx_nointermediate_close, other)
	#							Rest: distant / intermediate distance -> NA.
	rpw[, TYPE_PAIR_TO:= 'other']
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair_close','pair','pair_distant'))], 'TYPE_PAIR_TO', 'ancestral/\nintermingled')
	rpw[, TYPE_PAIR_TODI3:= 'with intermediate\nor distant']
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair_close','other_close'))], 'TYPE_PAIR_TODI3', 'no intermediate\n and close')
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair','other'))], 'TYPE_PAIR_TODI3', 'no intermediate\n but not close')
	rpw[, TYPE_PAIR_DI:= cut(PATRISTIC_DISTANCE, breaks=c(1e-12,0.02,0.05,2), labels=c('close','intermediate\ndistance','distant'))] 	
	rpw[, TYPE_PAIR_TODI2x2:= 'not close other']
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair_close'))], 'TYPE_PAIR_TODI2x2', 'close ancestral/\nintermingled')
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair','pair_distant'))], 'TYPE_PAIR_TODI2x2', 'not close ancestral/\nintermingled')
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('withintermediate_close','other_close'))], 'TYPE_PAIR_TODI2x2', 'close other')	
	rpw[, TYPE_DIR_TO3:= NA_character_]
	set(rpw, rpw[, which(grepl('chain_fm',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TO3', 'fm')
	set(rpw, rpw[, which(grepl('chain_mf',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TO3', 'mf')
	set(rpw, rpw[, which(grepl('intermingled',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TO3', 'intermingled')		
	rpw[, TYPE_DIR_TODI3:= NA_character_]
	set(rpw, rpw[, which(grepl('chain_fm',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'fm')
	set(rpw, rpw[, which(grepl('chain_mf',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'mf')
	set(rpw, rpw[, which(grepl('intermingled',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'intermingled')			
	rpw[, TYPE_CHAIN_TODI2:= NA_character_]
	set(rpw, rpw[, which(grepl('withintermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_CHAIN_TODI2', 'chain')
	set(rpw, rpw[, which(grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_CHAIN_TODI2', 'pair')
	set(rpw, rpw[, which(grepl('other',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_CHAIN_TODI2', 'pair')
	set(rpw, rpw[, which(grepl('intermingled',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'intermingled')			
	set(rpw, NULL, 'TYPE_DIR_TODI7x3', rpw[, gsub('intermediate',' intermediate',gsub('intermediate_','intermediate\n',gsub('intermingled_','intermingled\n',gsub('(chain_[fm][mf])_','\\1\n',TYPE_DIR_TODI7x3))))])
	#
	#	identify chunks of contiguous windows
	#	
	setkey(rpw, RUN, PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, W_FROM)
	rpw.slide	<- rpw[, {
				ans	<- NA_integer_
				tmp	<- diff(W_FROM)
				if(length(tmp))
					ans	<- min(tmp)
				list(W_SLIDE=ans)
			}, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','RUN')][, min(na.omit(W_SLIDE))]
	#	melt because chunks are dependent on topology types: if there are NAs, then the chunks may change
	rpw		<- melt(rpw, measure.vars=c('TYPE_PAIR_TO','TYPE_CHAIN_TODI2','TYPE_PAIR_TODI3x3','TYPE_PAIR_TODI3','TYPE_PAIR_DI','TYPE_PAIR_TODI2x2','TYPE_DIR_TODI7x3','TYPE_DIR_TO3','TYPE_DIR_TODI3'), variable.name='GROUP', value.name='TYPE')
	set(rpw, NULL, 'TYPE', rpw[,gsub('_',' ', TYPE)])	
	rpw		<- subset(rpw, !is.na(TYPE))			
	#	define chunks
	setkey(rpw, RUN, PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, GROUP, W_FROM)
	tmp		<- rpw[, {
				tmp<- as.integer( c(TRUE,(W_FROM[-length(W_FROM)]+rpw.slide)!=W_FROM[-1]) )
				list(W_FROM=W_FROM, W_TO=W_TO, CHUNK=cumsum(tmp))
			}, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','RUN','GROUP')]
	rpw		<- merge(rpw,tmp,by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','RUN','W_FROM','W_TO','GROUP'))
	#	define chunk length in terms of non-overlapping windows	& number of windows in chunk
	tmp		<- rpw[, {
				list(W_FROM=W_FROM, W_TO=W_TO, CHUNK_L=(max(W_TO+1L)-min(W_FROM))/(W_TO[1]+1L-W_FROM[1]), CHUNK_N=length(W_FROM))
			}, by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','GROUP','CHUNK')]
	rpw		<- merge(rpw,tmp,by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','RUN','GROUP','CHUNK','W_FROM','W_TO'))	
	#	for each chunk, count: windows by type and effective length of chunk
	#	then sum chunks
	rplkl	<- rpw[, list(	K= length(W_FROM), KEFF= length(W_FROM)/CHUNK_N[1] * CHUNK_L[1]), by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','CHUNK','LABEL','LABEL_SH','GROUP','TYPE')]	
	rplkl	<- rplkl[, list(STAT=c('K','KEFF'), V=c(sum(K),sum(KEFF))), by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','LABEL','LABEL_SH','GROUP','TYPE')]	
	#	add zeros
	rplkl	<- dcast.data.table(rplkl, RUN+PTY_RUN+FEMALE_SANGER_ID+MALE_SANGER_ID+LABEL+LABEL_SH~GROUP+TYPE+STAT, value.var='V')
	for(x in setdiff(colnames(rplkl),c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','LABEL','LABEL_SH','GROUP')))
		set(rplkl, which(is.na(rplkl[[x]])), x, 0L)	
	for(x in colnames(rplkl)[grepl('TYPE_PD_MEAN',colnames(rplkl))])
		set(rplkl, which(rplkl[[x]]>0), x, 1L)	
	rplkl	<- melt(rplkl, id.vars=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','LABEL','LABEL_SH'), variable.name='GROUP', value.name='V')
	rplkl[, STAT:= gsub('.*_([^_]+)$','\\1',GROUP)]
	set(rplkl, NULL, 'GROUP', rplkl[, gsub('(.*)_[^_]+$','\\1',GROUP)])	
	rplkl[, TYPE:= gsub('.*_([^_]+)$','\\1',GROUP)]
	set(rplkl, NULL, 'GROUP', rplkl[, gsub('(.*)_[^_]+$','\\1',GROUP)])	
	rplkl	<- dcast.data.table(rplkl, RUN+PTY_RUN+FEMALE_SANGER_ID+MALE_SANGER_ID+LABEL+LABEL_SH+GROUP+TYPE~STAT, value.var='V')	
	#	calculate N and NEFF
	tmp		<- rplkl[, list(N= sum(K), NEFF= sum(KEFF)), by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','GROUP')]	
	rplkl	<- merge(rplkl, tmp, by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','GROUP'))
	#	add parameters for marginal posterior probabilities (alpha, beta)
	par.prior	<- 0.1	
	tmp		<- rplkl[, list(	TYPE=TYPE, POSTERIOR_ALPHA=KEFF+par.prior, POSTERIOR_BETA=sum(KEFF+par.prior)-(KEFF+par.prior))	, by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','GROUP')]
	rplkl	<- merge(rplkl, tmp, by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','GROUP','TYPE'))
	#	save
	save(rplkl, file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmp_allpairs_posteriors.rda')
	#
	#	prepare colours
	#	
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
	#
	#	make detailed plots for each pair
	#		topology assignments across genome
	#		distances across genome
	#		number of windows of certain type
	#		estimated posterior probabilities on unknown phylogenetic relationship
	#
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]	

	groups		<- c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3','TYPE_PAIR_TODI2x2','TYPE_PAIR_TODI3','TYPE_PAIR_DI')
	group		<- c('TYPE_DIR_TODI7x3')
	
	widths	<- unit(c(4, 6), "null")
	heights	<- unit(c(2, 3.5, 4, 5), "null")
	height	<- 9
	if(group%in%c('TYPE_DIR_TODI7x3'))
	{
		widths	<- unit(c(4, 6), "null")
		heights	<- unit(c(2, 3.5, 4, 15), "null")
		height	<- 17
	}		
	if(group%in%c('TYPE_PAIR_TODI2x2'))
	{
		heights	<- unit(c(2, 3.5, 4, 3.75), "null")
		height	<- 8
	}
	if(group%in%c('TYPE_PAIR_TODI3','TYPE_PAIR_DI'))
	{
		heights	<- unit(c(2, 3.5, 4, 3.5), "null")
		height	<- 7
	}
	pdf(file=file.path(dir, paste(run,'-phsc-relationships_allpairs','_',group,'.pdf',sep='')), w=10, h=height)	
	plot.tmp	<- unique(subset(rplkl, GROUP==group, c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','LABEL')), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	setkey(plot.tmp, LABEL)		
	for(i in seq_len(nrow(plot.tmp)))
	{
		#pty_run	<- 38; id1		<- '16016_1_4'; id2		<- '15105_1_35'
		pty_run	<- plot.tmp[i, PTY_RUN]; id1		<- plot.tmp[i, FEMALE_SANGER_ID]; id2		<- plot.tmp[i, MALE_SANGER_ID]		
		tmp		<- subset(rpw, PTY_RUN==pty_run & GROUP==group & FEMALE_SANGER_ID==id1 & MALE_SANGER_ID==id2)
		p1		<- ggplot(tmp, aes(x=W_FROM)) +			
				geom_bar(aes(y=ID_R_MAX, colour=TYPE), stat='identity', fill='transparent') +
				geom_bar(aes(y=ID_R_MIN, fill=TYPE), stat='identity', colour='transparent') +
				labs(x='', y='number of reads', fill='phylogenetic\nrelationship\n', colour='phylogenetic\nrelationship\n') +
				scale_fill_manual(values=cols.type[[group]]) +
				scale_colour_manual(values=cols.type[[group]]) +
				scale_x_continuous(breaks=seq(0,1e4,500), minor_breaks=seq(0,1e4,100), limits=c(rpw[, min(W_FROM)], rpw[, max(W_FROM)])) +
				scale_y_log10(breaks=c(10,100,1000,1e4,1e5)) +
				theme_bw() + theme(legend.position='left') +			
				guides(fill=FALSE, colour=FALSE)
		p2		<- ggplot(tmp, aes(x=W_FROM, y=PATRISTIC_DISTANCE)) +
				geom_point(size=1) +					
				labs(x='window start\n\n', y='patristic distance') +
				scale_x_continuous(breaks=seq(0,1e4,500), minor_breaks=seq(0,1e4,100), limits=c(rpw[, min(W_FROM)], rpw[, max(W_FROM)])) +
				scale_y_log10(labels=percent, limits=c(0.001, 0.7), expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
				theme_bw() + theme(legend.position='left')
		tmp		<- subset(rplkl, GROUP==group & PTY_RUN==pty_run & FEMALE_SANGER_ID==id1 & MALE_SANGER_ID==id2)
		p3		<- ggplot(tmp, aes(x=TYPE, y=KEFF, fill=TYPE)) + geom_bar(stat='identity') +
				scale_fill_manual(values=cols.type[[group]]) +
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
				scale_fill_manual(values=cols.type[[group]]) +
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

RakaiAll.addposteriors.pairs.170227<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments_allpairs.rda')	
	tmp		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/170227/RCCS_170227_w270_trmw_assignments_allpairs.rda')
	set(tmp, NULL, setdiff(colnames(tmp),colnames(rpw)), NULL)
	rpw		<- rbind(rpw, tmp, fill=TRUE)	
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	#
	rpw[, table(RUN, useNA='if')]
	#	check if we have all pty.runs
	stopifnot(!length(setdiff( 	subset(rpw, RUN=='RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5')[, sort(unique(PTY_RUN))],
								subset(rpw, RUN=='RCCS_170227_w270_d50_st20_mr20_mt1_cl2_d5')[, sort(unique(PTY_RUN))]	)))	
	#	define plotting order: largest number of trm assignments	
	tmp		<- rpw[, list( WIN_TR=length(which(grepl('close|anc',TYPE_DIR_TODI7x3))) ), by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	tmp		<- tmp[, list(WIN_TR=max(WIN_TR)), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	tmp		<- tmp[order(-WIN_TR),]
	tmp[, PLOT_ID:=seq_len(nrow(tmp))]	
	#	make pair data table
	rp		<- copy(rpw)
	set(rp, NULL, c('RUN','FILE','DIR','W_FROM','W_TO','TYPE_RAW','TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3','PATRISTIC_DISTANCE','CONTIGUOUS','PATHS_12','PATHS_21','MALE_SANGER_ID_L','MALE_SANGER_ID_R','FEMALE_SANGER_ID_L','FEMALE_SANGER_ID_R'), NULL)
	rp		<- unique(rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN'))
	#	make COUPID
	rp[, COUPID:= paste0(MALE_RID,':',FEMALE_RID)]
	#	define label
	tmp		<- merge(tmp, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN'))
	setkey(tmp, PLOT_ID)
	tmp[, LABEL_SH:= factor(PLOT_ID, levels=PLOT_ID, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair ', COUPID,' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL, LABEL_SH))
	rpw		<- merge(tmp, rpw, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	#	define min/max reads
	tmp		<- rpw[, list(	ID_R_MIN=min(MALE_SANGER_ID_R, FEMALE_SANGER_ID_R),
					ID_R_MAX=max(MALE_SANGER_ID_R, FEMALE_SANGER_ID_R)), by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM')]
	rpw		<- merge(rpw, tmp, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM'))
	#
	#	define topo/distance types of phylogenetic relationship between pairs
	#
	#	given the bimodal distribution of patristic distances from phyloscanner
	#	I ususally consider 3 distance states: close, neither close nor distant, distant
	#
	#	TYPE_DIR_TODI7x3		7 topology states (this is the "basic topological characterization" on Slack that I derive from Matthew's output) 
	#							chain_12_nointermediate, chain_12_withintermediate, chain_21_nointermediate, chain_21_withintermediate, intermingled_nointermediate, intermingled_withintermediate, other  
	#							times 
	#							3 distance states
	#
	#	TYPE_PAIR_TODI3x3		3 topology states
	#							pair (chain_XX_nointermediate+intermingled_nointermediate), withintermediate (XX_withintermediate), other
	#							times 
	#							3 distance states
	#
	#	TYPE_PAIR_TODI2x2		2 topology states
	#							ancestral/intermingled, other
	#							times 
	#							2 distance states (close, not close)
	#
	#	TYPE_PAIR_TODI3			non-factorial design that combines distant or withintermediate, and close or ancestral/intermingled
	#							the idea here is that 
	#							evidence for extra-couple transmission comes from large patristic distance OR intermediates	
	#							evidence for extra-couple transmission comes from ancestral/intermingled OR short patristic distance	
	#
	#	TYPE_PAIR_DI			3 distance states
	#
	#
	#	TYPE_DIR_TO3			3 topological direction states
	#							chain_fm	(chain_fm_nointermediate regardless of distance); chain_mf	(chain_mf_nointermediate regardless of distance); other
	#							Rest: intermingled -> NA.
	#
	#	TYPE_DIR_TODI3			3 topological direction states if close.  
	#							chain_fm	(chain_fm_nointermediate regardless of distance); chain_mf	(chain_mf_nointermediate regardless of distance); other
	#							Rest: intermingled / distant / intermediate distance -> NA.
	#
	#	TYPE_CHAIN_TODI2		2 topological states if close.  
	#							chain	(xxx_withintermediate_close); pair	(xxx_nointermediate_close, other)
	#							Rest: distant / intermediate distance -> NA.
	rpw[, TYPE_PAIR_TO:= 'other']
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair_close','pair','pair_distant'))], 'TYPE_PAIR_TO', 'ancestral/\nintermingled')
	rpw[, TYPE_PAIR_TODI3:= 'with intermediate\nor distant']
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair_close','other_close'))], 'TYPE_PAIR_TODI3', 'no intermediate\n and close')
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair','other'))], 'TYPE_PAIR_TODI3', 'no intermediate\n but not close')
	rpw[, TYPE_PAIR_DI:= cut(PATRISTIC_DISTANCE, breaks=c(1e-12,0.02,0.05,2), labels=c('close','intermediate\ndistance','distant'))] 	
	rpw[, TYPE_PAIR_TODI2x2:= 'not close other']
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair_close'))], 'TYPE_PAIR_TODI2x2', 'close ancestral/\nintermingled')
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair','pair_distant'))], 'TYPE_PAIR_TODI2x2', 'not close ancestral/\nintermingled')
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('withintermediate_close','other_close'))], 'TYPE_PAIR_TODI2x2', 'close other')	
	rpw[, TYPE_DIR_TO3:= NA_character_]
	set(rpw, rpw[, which(grepl('chain_fm',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TO3', 'fm')
	set(rpw, rpw[, which(grepl('chain_mf',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TO3', 'mf')
	set(rpw, rpw[, which(grepl('intermingled',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TO3', 'intermingled')		
	rpw[, TYPE_DIR_TODI3:= NA_character_]
	set(rpw, rpw[, which(grepl('chain_fm',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'fm')
	set(rpw, rpw[, which(grepl('chain_mf',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'mf')
	set(rpw, rpw[, which(grepl('intermingled',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'intermingled')			
	rpw[, TYPE_CHAIN_TODI2:= NA_character_]
	set(rpw, rpw[, which(grepl('withintermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_CHAIN_TODI2', 'chain')
	set(rpw, rpw[, which(grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_CHAIN_TODI2', 'pair')
	set(rpw, rpw[, which(grepl('other',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_CHAIN_TODI2', 'pair')
	set(rpw, rpw[, which(grepl('intermingled',TYPE_DIR_TODI7x3) & grepl('nointermediate',TYPE_DIR_TODI7x3) & grepl('close',TYPE_DIR_TODI7x3))], 'TYPE_DIR_TODI3', 'intermingled')			
	set(rpw, NULL, 'TYPE_DIR_TODI7x3', rpw[, gsub('intermediate',' intermediate',gsub('intermediate_','intermediate\n',gsub('intermingled_','intermingled\n',gsub('(chain_[fm][mf])_','\\1\n',TYPE_DIR_TODI7x3))))])
	#
	#	identify chunks of contiguous windows
	#	
	setkey(rpw, RUN, PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, W_FROM)
	rpw.slide	<- rpw[, {
				ans	<- NA_integer_
				tmp	<- diff(W_FROM)
				if(length(tmp))
					ans	<- min(tmp)
				list(W_SLIDE=ans)
			}, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','RUN')][, min(na.omit(W_SLIDE))]
	#	melt because chunks are dependent on topology types: if there are NAs, then the chunks may change
	rpw		<- melt(rpw, measure.vars=c('TYPE_PAIR_TO','TYPE_CHAIN_TODI2','TYPE_PAIR_TODI3x3','TYPE_PAIR_TODI3','TYPE_PAIR_DI','TYPE_PAIR_TODI2x2','TYPE_DIR_TODI7x3','TYPE_DIR_TO3','TYPE_DIR_TODI3'), variable.name='GROUP', value.name='TYPE')
	set(rpw, NULL, 'TYPE', rpw[,gsub('_',' ', TYPE)])	
	rpw		<- subset(rpw, !is.na(TYPE))			
	#	define chunks
	setkey(rpw, RUN, PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, GROUP, W_FROM)
	tmp		<- rpw[, {
				tmp<- as.integer( c(TRUE,(W_FROM[-length(W_FROM)]+rpw.slide)!=W_FROM[-1]) )
				list(W_FROM=W_FROM, W_TO=W_TO, CHUNK=cumsum(tmp))
			}, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','GROUP')]
	rpw		<- merge(rpw,tmp,by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM','W_TO','GROUP'))
	#	define chunk length in terms of non-overlapping windows	& number of windows in chunk
	tmp		<- rpw[, {
				list(W_FROM=W_FROM, W_TO=W_TO, CHUNK_L=(max(W_TO+1L)-min(W_FROM))/(W_TO[1]+1L-W_FROM[1]), CHUNK_N=length(W_FROM))
			}, by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','GROUP','CHUNK')]
	rpw		<- merge(rpw,tmp,by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','GROUP','CHUNK','W_FROM','W_TO'))	
	#	for each chunk, count: windows by type and effective length of chunk
	#	then sum chunks
	rplkl	<- rpw[, list(	K= length(W_FROM), KEFF= length(W_FROM)/CHUNK_N[1] * CHUNK_L[1]), by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','CHUNK','LABEL','LABEL_SH','GROUP','TYPE')]	
	rplkl	<- rplkl[, list(STAT=c('K','KEFF'), V=c(sum(K),sum(KEFF))), by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','LABEL','LABEL_SH','GROUP','TYPE')]	
	#	add zeros
	rplkl	<- dcast.data.table(rplkl, RUN+PTY_RUN+FEMALE_SANGER_ID+MALE_SANGER_ID+LABEL+LABEL_SH~GROUP+TYPE+STAT, value.var='V')
	for(x in setdiff(colnames(rplkl),c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','LABEL','LABEL_SH','GROUP')))
		set(rplkl, which(is.na(rplkl[[x]])), x, 0L)	
	for(x in colnames(rplkl)[grepl('TYPE_PD_MEAN',colnames(rplkl))])
		set(rplkl, which(rplkl[[x]]>0), x, 1L)	
	rplkl	<- melt(rplkl, id.vars=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','LABEL','LABEL_SH'), variable.name='GROUP', value.name='V')
	rplkl[, STAT:= gsub('.*_([^_]+)$','\\1',GROUP)]
	set(rplkl, NULL, 'GROUP', rplkl[, gsub('(.*)_[^_]+$','\\1',GROUP)])	
	rplkl[, TYPE:= gsub('.*_([^_]+)$','\\1',GROUP)]
	set(rplkl, NULL, 'GROUP', rplkl[, gsub('(.*)_[^_]+$','\\1',GROUP)])	
	rplkl	<- dcast.data.table(rplkl, RUN+PTY_RUN+FEMALE_SANGER_ID+MALE_SANGER_ID+LABEL+LABEL_SH+GROUP+TYPE~STAT, value.var='V')	
	#	calculate N and NEFF
	tmp		<- rplkl[, list(N= sum(K), NEFF= sum(KEFF)), by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','GROUP')]	
	rplkl	<- merge(rplkl, tmp, by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','GROUP'))
	#	add parameters for marginal posterior probabilities (alpha, beta)
	par.prior	<- 0.1	
	tmp		<- rplkl[, list(	TYPE=TYPE, POSTERIOR_ALPHA=KEFF+par.prior, POSTERIOR_BETA=sum(KEFF+par.prior)-(KEFF+par.prior))	, by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','GROUP')]
	rplkl	<- merge(rplkl, tmp, by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','GROUP','TYPE'))
	#	save
	save(rplkl, file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/170227/RCCS_170227_w270_trmp_allpairs_posteriors.rda')
	#
	#	check runs missing?
	#
	tmp		<- unique(subset(rplkl, select=c(RUN, PTY_RUN)))
	tmp		<- merge(tmp,tmp[,list(N_RUNS=length(RUN)),by='PTY_RUN'],by='PTY_RUN')
	dcast.data.table(subset(tmp, N_RUNS<4), PTY_RUN~RUN, value.var='N_RUNS')
	
	#	
	#	check presence/absence by dual
	#
	tmp		<- unique(subset(rplkl, select=c(RUN, PTY_RUN, FEMALE_SANGER_ID, MALE_SANGER_ID)))
	tmp		<- merge(tmp, tmp[, list(N_RUNS=length(RUN)),by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID')],by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID'))
	
	subset(tmp, N_RUNS<4)[, table(PTY_RUN)]
	#
	#	prepare colours
	#	
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
	#
	#	make detailed plots for each pair
	#		topology assignments across genome
	#		distances across genome
	#		number of windows of certain type
	#		estimated posterior probabilities on unknown phylogenetic relationship
	#
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]	
	
	groups		<- c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3','TYPE_PAIR_TODI2x2','TYPE_PAIR_TODI3','TYPE_PAIR_DI')
	group		<- c('TYPE_DIR_TODI7x3')
	
	widths	<- unit(c(4, 6), "null")
	heights	<- unit(c(2, 3.5, 4, 5), "null")
	height	<- 9
	if(group%in%c('TYPE_DIR_TODI7x3'))
	{
		widths	<- unit(c(4, 6), "null")
		heights	<- unit(c(2, 3.5, 4, 15), "null")
		height	<- 17
	}		
	if(group%in%c('TYPE_PAIR_TODI2x2'))
	{
		heights	<- unit(c(2, 3.5, 4, 3.75), "null")
		height	<- 8
	}
	if(group%in%c('TYPE_PAIR_TODI3','TYPE_PAIR_DI'))
	{
		heights	<- unit(c(2, 3.5, 4, 3.5), "null")
		height	<- 7
	}
	pdf(file=file.path(dir, paste(run,'-phsc-relationships_allpairs','_',group,'.pdf',sep='')), w=10, h=height)	
	plot.tmp	<- unique(subset(rplkl, GROUP==group, c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','LABEL')), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	setkey(plot.tmp, LABEL)		
	for(i in seq_len(nrow(plot.tmp)))
	{
		#pty_run	<- 38; id1		<- '16016_1_4'; id2		<- '15105_1_35'
		pty_run	<- plot.tmp[i, PTY_RUN]; id1		<- plot.tmp[i, FEMALE_SANGER_ID]; id2		<- plot.tmp[i, MALE_SANGER_ID]		
		tmp		<- subset(rpw, PTY_RUN==pty_run & GROUP==group & FEMALE_SANGER_ID==id1 & MALE_SANGER_ID==id2)
		p1		<- ggplot(tmp, aes(x=W_FROM)) +			
				geom_bar(aes(y=ID_R_MAX, colour=TYPE), stat='identity', fill='transparent') +
				geom_bar(aes(y=ID_R_MIN, fill=TYPE), stat='identity', colour='transparent') +
				labs(x='', y='number of reads', fill='phylogenetic\nrelationship\n', colour='phylogenetic\nrelationship\n') +
				scale_fill_manual(values=cols.type[[group]]) +
				scale_colour_manual(values=cols.type[[group]]) +
				scale_x_continuous(breaks=seq(0,1e4,500), minor_breaks=seq(0,1e4,100), limits=c(rpw[, min(W_FROM)], rpw[, max(W_FROM)])) +
				scale_y_log10(breaks=c(10,100,1000,1e4,1e5)) +
				theme_bw() + theme(legend.position='left') +			
				guides(fill=FALSE, colour=FALSE)
		p2		<- ggplot(tmp, aes(x=W_FROM, y=PATRISTIC_DISTANCE)) +
				geom_point(size=1) +					
				labs(x='window start\n\n', y='patristic distance') +
				scale_x_continuous(breaks=seq(0,1e4,500), minor_breaks=seq(0,1e4,100), limits=c(rpw[, min(W_FROM)], rpw[, max(W_FROM)])) +
				scale_y_log10(labels=percent, limits=c(0.001, 0.7), expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
				theme_bw() + theme(legend.position='left')
		tmp		<- subset(rplkl, GROUP==group & PTY_RUN==pty_run & FEMALE_SANGER_ID==id1 & MALE_SANGER_ID==id2)
		p3		<- ggplot(tmp, aes(x=TYPE, y=KEFF, fill=TYPE)) + geom_bar(stat='identity') +
				scale_fill_manual(values=cols.type[[group]]) +
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
				scale_fill_manual(values=cols.type[[group]]) +
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

RakaiCouples.extracoupletrm.170106.analyze.extratrm.couples.age<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)	
	require(coin)
	require(Hmisc)
	require(gamlss)
	# load pty.runs
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]	
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	# load posterior relationships
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_extracoupletrm_bbmodels.rda')
	# add agrarian etc
	# 	from Kate:
	#	Trading communities: 1   4   16  22  24  33  51  107 776
	#	Fishing communities: 38, 770, 771, 774	
	# 	tmp		<- rp[, sort(unique(c(MALE_COMM_NUM,FEMALE_COMM_NUM)))]; cat(paste('"', tmp, '",', collapse='', sep=''))
	tmp		<- data.table(	COMM_NUM=	c("1","2","4","7","8","16","22","23","24","33","34","38","40","51","56","57","58","89","94","106","107","108","370","391","770","771","772","773","774","776"),
			COMM_TYPE=	c("T","A","T","A","A", "T", "T", "A", "T", "T", "A", "F", "A", "T", "A", "A", "A", "A", "A",  "A",  "T",  "A",  "A",  "A", "F",  "F",  "A",  "A",  "F",   "T"))
	set(tmp, NULL, 'COMM_TYPE', tmp[, as.character(factor(COMM_TYPE, levels=c('A','T','F'), labels=c('agrarian','trading','fisherfolk')))])
	set(tmp, NULL, 'COMM_NUM', tmp[, as.integer(COMM_NUM)])
	setnames(tmp, c('COMM_NUM','COMM_TYPE'), c('MALE_COMM_NUM','MALE_COMM_TYPE'))
	rp		<- merge(rp, tmp, by='MALE_COMM_NUM')
	setnames(tmp, c('MALE_COMM_NUM','MALE_COMM_TYPE'), c('FEMALE_COMM_NUM','FEMALE_COMM_TYPE'))
	rp		<- merge(rp, tmp, by='FEMALE_COMM_NUM')
	#
	#	select couples	
	#
	rpd		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='with intermediate\nor distant')	
	#	make central selection: post prob of extra-couple trm > 50%
	rpd[, SELECT_P50:= as.character(factor(pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>0.8, levels=c(TRUE,FALSE), labels=c('Y','N')))]
	rpd[, SELECT_P80:= as.character(factor(pbeta(0.8, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>0.8, levels=c(TRUE,FALSE), labels=c('Y','N')))]
	rpd		<- subset(rpd, select=c(PTY_RUN, FEMALE_SANGER_ID, MALE_SANGER_ID, NEFF, KEFF, POSTERIOR_ALPHA, POSTERIOR_BETA, LABEL_SH, LABEL, SELECT_P50, SELECT_P80))
	rpd		<- merge(rp, rpd, by=c('FEMALE_SANGER_ID','MALE_SANGER_ID'), all.x=1)
	set(rpd, rpd[, which(is.na(SELECT_P50))], 'SELECT_P50', 'No_PANGEA_data')
	set(rpd, rpd[, which(is.na(SELECT_P80))], 'SELECT_P80', 'No_PANGEA_data')
	#
	#	among cohabiting couples: age at diagnosis vs extra-couple transmission
	#
	#	results:
	#	the men are substantially younger in extra-transmission couples		
	#	this is not the case for women (if anything they are older in extra-transmission couples)
	tmp		<- unique(subset(rpd, MALE_HH_NUM==FEMALE_HH_NUM & SELECT_P50%in%c('Y','N')), by='COUPID')
	tmp[, MALE_FIRSTPOSDATE_AGE:= MALE_FIRSTPOSDATE-MALE_BIRTHDATE]
	tmp[, FEMALE_FIRSTPOSDATE_AGE:= FEMALE_FIRSTPOSDATE-FEMALE_BIRTHDATE]
	tmp		<- merge(tmp, tmp[, list(	MALE_BOTHPOS_AGE= max(MALE_FIRSTPOSDATE, FEMALE_FIRSTPOSDATE)-MALE_BIRTHDATE, 
							FEMALE_BOTHPOS_AGE= max(MALE_FIRSTPOSDATE, FEMALE_FIRSTPOSDATE)-FEMALE_BIRTHDATE,
							BOTHPOS_AVGAGE= max(MALE_FIRSTPOSDATE, FEMALE_FIRSTPOSDATE)-(FEMALE_BIRTHDATE+FEMALE_BIRTHDATE)/2), by='COUPID'], by='COUPID')	
	tmp[, FIRSTPOSDATE_AGE_MFDIFF:= MALE_BIRTHDATE-FEMALE_BIRTHDATE]
	tmp2	<- melt(tmp, measure.vars=c('MALE_FIRSTPOSDATE_AGE','FEMALE_FIRSTPOSDATE_AGE'))
	set(tmp2, NULL,'variable', tmp2[,factor(variable, 	levels=c('MALE_FIRSTPOSDATE_AGE','FEMALE_FIRSTPOSDATE_AGE'),labels=c('male partner','female partner'))])
	set(tmp2, NULL,'SELECT_P50', tmp2[, factor(SELECT_P50, levels=c('Y','N'), labels=c('Yes','No'))])
	ggplot(tmp2, aes(x=SELECT_P50, y=value, fill=variable, alpha=SELECT_P50)) + geom_boxplot(outlier.shape=NA) +			
			theme_bw() + labs(x='\nphylogenetically unlinked virus among +/+ couples\n(posterior probability >50% with confidence >80%)', y='age at diagnosis\n') +
			scale_y_continuous(breaks=seq(10,80,5)) +
			scale_fill_brewer(palette='Set1') +
			scale_alpha_manual(values=c('Yes'=1, 'No'=0.5)) +
			facet_grid(~variable) +
			guides(fill='none', alpha='none')
	ggsave(file=file.path(dir, paste0(run,'-phsc-extracoupletrm_couples_all_by_agepartners.pdf')), w=4, h=6)
	#
	ggplot(subset(tmp2, COUP_SC%in%c('seroinc','M->F','F->M')), aes(x=SELECT_P50, y=value, fill=variable, alpha=SELECT_P50)) + geom_boxplot(outlier.shape=NA) +			
			theme_bw() + labs(x='\nextra-couple trm among incident / initially discordant couples\n(posterior probability >50% with confidence >80%)', y='age at diagnosis\n') +
			scale_y_continuous(breaks=seq(10,80,5)) +
			scale_fill_brewer(palette='Set1') +
			scale_alpha_manual(values=c('Yes'=1, 'No'=0.5)) +
			facet_grid(~variable) +
			guides(fill='none', alpha='none')
	ggsave(file=file.path(dir, paste0(run,'-phsc-extracoupletrm_couples_incdisc_by_agepartners.pdf')), w=4, h=6)	
	#
	#	age difference versus unlinked virus
	#
	#	results: no difference in age differences between male/female in linked/unlinked couples
	tmp2	<- subset(tmp, !is.na(FEMALE_FIRSTPOSDATE_AGE) & !is.na(MALE_FIRSTPOSDATE_AGE))
	tmp2[, COUP_SC2:= '+/+ at enrollment']
	set(tmp2, tmp2[, which(COUP_SC%in%c('seroinc','M->F','F->M'))], 'COUP_SC2', '-/-, -/+, +/- at enrollment')
	set(tmp2, NULL,'SELECT_P50', tmp2[, factor(SELECT_P50, levels=c('Y','N'), labels=c('Yes','No'))])
	ggplot(tmp2, aes(x=SELECT_P50, y=FIRSTPOSDATE_AGE_MFDIFF, fill=SELECT_P50)) + 
			geom_boxplot(outlier.shape=NA) +
			scale_fill_brewer(palette='PRGn') +
			scale_alpha_manual(values=c('Yes'=1, 'No'=0.5)) +			
			scale_y_continuous(breaks=seq(-80,80,5),limits=c(-10,25), expand=c(0,0)) +
			labs(x='\nphylogenetically unlinked virus among male and female partner\n(posterior probability >50% with confidence >80%)', y='age difference\n(years)\n') +
			facet_grid(~COUP_SC2) +
			theme_bw() + guides(fill='none')
	ggsave(file=file.path(dir, paste0(run,'-phsc-extracoupletrm_couples_incdisc_by_agedifferencepartners.pdf')), w=4, h=6)
	#	
	#	among -/- -/+ +/- couples, BOTH the male and female partners are younger in unlinked couples
	tmp2	<- subset(tmp, !is.na(FEMALE_FIRSTPOSDATE_AGE) & !is.na(MALE_FIRSTPOSDATE_AGE))
	tmp2[, COUP_SC2:= '+/+ at enrollment']
	set(tmp2, tmp2[, which(COUP_SC%in%c('seroinc','M->F','F->M'))], 'COUP_SC2', '-/-, -/+, +/- at enrollment')
	set(tmp2, NULL,'SELECT_P50', tmp2[, factor(SELECT_P50, levels=c('Y','N'), labels=c('Yes','No'))])		
	tmp2	<- melt(tmp2, measure.vars=c('MALE_BIRTHDATE','FEMALE_BIRTHDATE'))
	set(tmp2, NULL, 'variable', tmp2[, factor(variable, levels=c('MALE_BIRTHDATE','FEMALE_BIRTHDATE'), labels=c('male','female'))])
	ggplot(tmp2, aes(x=SELECT_P50, y=value, fill=SELECT_P50)) +
			geom_boxplot(outlier.shape=NA) +
			scale_fill_brewer(palette='PRGn') +
			facet_grid(COUP_SC2~variable) +
			theme_bw() + guides(fill='none') +
			labs(x='\nphylogenetically unlinked virus among male and female partner\n(posterior probability >50% with confidence >80%)', y='birth date')
	ggsave(file=file.path(dir, paste0(run,'-phsc-extracoupletrm_couples_incdisc_by_birthdate.pdf')), w=6, h=6)
	#
	#	age when both partners are positive
	#
	tmp2	<- subset(tmp, !is.na(FEMALE_FIRSTPOSDATE_AGE) & !is.na(MALE_FIRSTPOSDATE_AGE))
	tmp2[, COUP_SC2:= '+/+ at enrollment']
	set(tmp2, tmp2[, which(COUP_SC%in%c('seroinc','M->F','F->M'))], 'COUP_SC2', '-/-, -/+, +/- at enrollment')
	set(tmp2, NULL,'SELECT_P50', tmp2[, factor(SELECT_P50, levels=c('Y','N'), labels=c('Yes','No'))])	
	tmp2	<- melt(tmp2, measure.vars=c('MALE_BOTHPOS_AGE','FEMALE_BOTHPOS_AGE'))
	set(tmp2, NULL, 'variable', tmp2[, factor(variable, levels=c('MALE_BOTHPOS_AGE','FEMALE_BOTHPOS_AGE'), labels=c('male','female'))])
	ggplot(tmp2, aes(x=SELECT_P50, y=value, fill=SELECT_P50)) +
			geom_boxplot(outlier.shape=NA) +
			scale_fill_brewer(palette='PRGn') +
			facet_grid(COUP_SC2~variable) +
			theme_bw() + guides(fill='none') +
			labs(	x='\nphylogenetically unlinked virus among male and female partner\n(posterior probability >50% with confidence >80%)', 
					y='age at time both partners were first positive')
	ggsave(file=file.path(dir, paste0(run,'-phsc-extracoupletrm_couples_incdisc_by_agebothpos.pdf')), w=6, h=6)
	#
	#	average age when both partners are positive
	#
	tmp2	<- subset(tmp, !is.na(FEMALE_FIRSTPOSDATE_AGE) & !is.na(MALE_FIRSTPOSDATE_AGE))
	tmp2[, COUP_SC2:= '+/+ at enrollment']
	set(tmp2, tmp2[, which(COUP_SC%in%c('seroinc','M->F','F->M'))], 'COUP_SC2', '-/-, -/+, +/- at enrollment')
	set(tmp2, NULL,'SELECT_P50', tmp2[, factor(SELECT_P50, levels=c('Y','N'), labels=c('Yes','No'))])	
	ggplot(tmp2, aes(x=SELECT_P50, y=BOTHPOS_AVGAGE, fill=SELECT_P50)) +
			geom_boxplot(outlier.shape=NA) +
			scale_fill_brewer(palette='PRGn') +
			facet_grid(~COUP_SC2) +
			theme_bw() + guides(fill='none') +
			labs(	x='\nphylogenetically unlinked virus among male and female partner\n(posterior probability >50% with confidence >80%)', 
					y='average age of both partners at time both partners were first positive')
	ggsave(file=file.path(dir, paste0(run,'-phsc-extracoupletrm_couples_incdisc_by_avgagebothpos.pdf')), w=4, h=6)
	#
	#	regression model
	#
	#	none significant in multivariable analysis
	tmp		<- unique(subset(rpd, 	MALE_HH_NUM==FEMALE_HH_NUM & 
							COUP_SC%in%c('seroinc','M->F','F->M') &
							!is.na(FEMALE_BIRTHDATE) & !is.na(MALE_BIRTHDATE) & 
							SELECT_P50%in%c('Y','N')), by='COUPID')
	tmp[, MALE_FIRSTPOSDATE_AGE:= MALE_FIRSTPOSDATE-MALE_BIRTHDATE]
	tmp[, FEMALE_FIRSTPOSDATE_AGE:= FEMALE_FIRSTPOSDATE-FEMALE_BIRTHDATE]
	tmp		<- merge(tmp, tmp[, list(	MALE_BOTHPOS_AGE= max(MALE_FIRSTPOSDATE, FEMALE_FIRSTPOSDATE)-MALE_BIRTHDATE, 
							FEMALE_BOTHPOS_AGE= max(MALE_FIRSTPOSDATE, FEMALE_FIRSTPOSDATE)-FEMALE_BIRTHDATE,
							BOTHPOS_AVGAGE= max(MALE_FIRSTPOSDATE, FEMALE_FIRSTPOSDATE)-(FEMALE_BIRTHDATE+FEMALE_BIRTHDATE)/2), by='COUPID'], by='COUPID')	
	tmp[, FIRSTPOSDATE_AGE_MFDIFF:= MALE_BIRTHDATE-FEMALE_BIRTHDATE]
	set(tmp, NULL,'SELECT_P50', tmp[, as.integer(as.character(factor(SELECT_P50, levels=c('Y','N'), labels=c('1','0'))))])	
	tmp2	<- subset(tmp, select=c(SELECT_P50, MALE_COMM_TYPE, BOTHPOS_AVGAGE, FEMALE_BOTHPOS_AGE, MALE_BOTHPOS_AGE, FIRSTPOSDATE_AGE_MFDIFF))	
	summary(gamlss(data=tmp2, SELECT_P50~MALE_COMM_TYPE, family=LO))
	summary(gamlss(data=tmp2, SELECT_P50~MALE_COMM_TYPE+BOTHPOS_AVGAGE, family=LO))
	summary(gamlss(data=tmp2, SELECT_P50~MALE_COMM_TYPE+BOTHPOS_AVGAGE+FIRSTPOSDATE_AGE_MFDIFF, family=LO))	
}

RakaiCouples.extracoupletrm.170106.analyze.extratrm.couples.communities.auxsampling<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)	
	require(coin)
	require(Hmisc)
	require(gamlss)
	# load pty.runs
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]	
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	# load posterior relationships
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_extracoupletrm_bbmodels.rda')
	# add agrarian etc
	# 	from Kate:
	#	Trading communities: 1   4   16  22  24  33  51  107 776
	#	Fishing communities: 38, 770, 771, 774	
	# 	tmp		<- rp[, sort(unique(c(MALE_COMM_NUM,FEMALE_COMM_NUM)))]; cat(paste('"', tmp, '",', collapse='', sep=''))
	tmp		<- data.table(	COMM_NUM=	c("1","2","4","7","8","16","22","23","24","33","34","38","40","51","56","57","58","89","94","106","107","108","370","391","770","771","772","773","774","776"),
							COMM_TYPE=	c("T","A","T","A","A", "T", "T", "A", "T", "T", "A", "F", "A", "T", "A", "A", "A", "A", "A",  "A",  "T",  "A",  "A",  "A", "F",  "F",  "A",  "A",  "F",   "T"))
	set(tmp, NULL, 'COMM_TYPE', tmp[, as.character(factor(COMM_TYPE, levels=c('A','T','F'), labels=c('agrarian','trading','fisherfolk')))])
	set(tmp, NULL, 'COMM_NUM', tmp[, as.integer(COMM_NUM)])
	setnames(tmp, c('COMM_NUM','COMM_TYPE'), c('MALE_COMM_NUM','MALE_COMM_TYPE'))
	rp		<- merge(rp, tmp, by='MALE_COMM_NUM')
	setnames(tmp, c('MALE_COMM_NUM','MALE_COMM_TYPE'), c('FEMALE_COMM_NUM','FEMALE_COMM_TYPE'))
	rp		<- merge(rp, tmp, by='FEMALE_COMM_NUM')
	#
	#	select cohabiting couples	
	#
	rpd		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='with intermediate\nor distant')
	rpd		<- merge(rp, rpd, by=c('FEMALE_SANGER_ID','MALE_SANGER_ID','COUP_SC'))
	rpd		<- unique(subset(rpd, MALE_HH_NUM==FEMALE_HH_NUM), by='COUPID')
	#
	#	Monte Carlo procedure that defines an auxiliary variable Z:= truely unlinked 
	#		where Z has as pdf the posterior distribution from phyloscanner  
	#		this fully accounts for uncertainty in the phylogenetic relationships 
	#		and avoids cut-offs
	
	#	select couples
	rpds	<- subset(rpd, select=c(COUPID, COUP_SC, MALE_COMM_TYPE, POSTERIOR_ALPHA, POSTERIOR_BETA))
	rpds[, COUP_SC:='all +/+ couples']
	tmp2	<- subset(rpd, COUP_SC%in%c('seroinc','M->F','F->M'), select=c(COUPID, COUP_SC, MALE_COMM_TYPE, POSTERIOR_ALPHA, POSTERIOR_BETA))
	tmp2[, COUP_SC:='-/-, -/+, +/- at enrollment']
	rpds	<- rbind(rpds, tmp2)
	#
	tmp2	<- tmp2[, list(N=length(COUPID)), by=c('MALE_COMM_TYPE','COUP_SC')]
	setkey(tmp2, MALE_COMM_TYPE, COUP_SC)
	tmp3	<- data.table(	COUP_SC			= c('seroinc', 'seroinc',   'seroinc','M->F',    'M->F',      'M->F',   'F->M',    'F->M',      'F->M'),
							MALE_COMM_TYPE	= c('agrarian','fisherfolk','trading','agrarian','fisherfolk','trading','agrarian','fisherfolk','trading'),
							EXTRACOUPLE		= c((.05+.29)/2,(.16+.56)/2,(.08+.42)/2,.05,      .16,         .08,      .29,       .56,         .42))
	tmp2	<- merge(tmp2,tmp3,by=c('MALE_COMM_TYPE','COUP_SC'))
	tmp2[, list(EXTRACOUPLE= sum( N/sum(N)*EXTRACOUPLE) ), by='MALE_COMM_TYPE']
	#   MALE_COMM_TYPE EXTRACOUPLE
	#1:       agrarian   0.1300000
	#2:     fisherfolk   0.3535484
	#3:        trading   0.1650000
	
	#	
	mc.it	<- 1e4	
	set.seed(42)
	# 1e7 iterations seem reasonable to get accurate mean to tolerance 1e-3
	# 	mean(rbinom(1e7,1,rbeta(1e7, 0.1, 1.2))) - 0.1/(0.1+1.2)
	rpds	<- do.call('rbind',lapply(seq_len(mc.it), function(i)
					{
						#	one Monte Carlo iteration:
						#	1- draw Z from posterior probabilities
						#	2- count successes Z=1 (unlinked couples) and evaluate binomial prob and confidence intervals
						#i		<- 1
						if(i%%1e2==0) cat('\n',i)
						mc		<- rpds[, list(	MC_IT=i, 
												MC_Z= rbinom(1,1,rbeta(1, POSTERIOR_ALPHA, POSTERIOR_BETA))), by=c('COUPID','COUP_SC','MALE_COMM_TYPE')]
						tmp		<- mc[, 	{
									z	<- binconf( length(which(MC_Z==1)), length(MC_Z) )				
									list(K=length(which(MC_Z==1)), N=length(MC_Z), P=z[1], QL=z[2], QU=z[3])
								}, by=c('MALE_COMM_TYPE','COUP_SC','MC_IT')]	
						tmp						
					}))
	save(rpds, file=file.path(dir, paste0(run,'-phsc-extracoupletrm_couples_all_by_commtype_MonteCarlo.rda')))
	# 
	tmp		<- rpds[, list(N=mean(N), P=mean(P), QL=mean(QL), QU=mean(QU)), by=c('COUP_SC','MALE_COMM_TYPE')]
	set(tmp, NULL, 'MALE_COMM_TYPE', tmp[, factor(MALE_COMM_TYPE, levels=c('agrarian','trading','fisherfolk'))])
	set(tmp, tmp[, which(COUP_SC=='-/-, -/+, +/- at enrollment')],'COUP_SC','+/+ couples,\n-/-, -/+, +/- at enrollment')
	#                                     COUP_SC MALE_COMM_TYPE   N         P         QL        QU
	#1:                           all +/+ couples       agrarian  53 0.2336755 0.14093054 0.3624179
	#2:                           all +/+ couples     fisherfolk 146 0.4148418 0.33822165 0.4958283
	#3:                           all +/+ couples        trading   4 0.4676250 0.13353785 0.8294649
	#4: +/+ couples,\n-/-, -/+, +/- at enrollment       agrarian  18 0.2157944 0.08888292 0.4425456
	#5: +/+ couples,\n-/-, -/+, +/- at enrollment     fisherfolk  31 0.4611935 0.30106451 0.6298798
	#6: +/+ couples,\n-/-, -/+, +/- at enrollment        trading   2 0.4442000 0.03338005 0.9284103
	ggplot(tmp, aes(x=MALE_COMM_TYPE)) +
			geom_point(aes(y=P)) + geom_errorbar(aes(ymax=QU, ymin=QL), width=0.5) +
			geom_text(aes(y=QU+.05, label=paste0('n=',N))) + 
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous(limits=c(0,1), labels=percent, expand=c(0,0), breaks=seq(0,1,0.2)) +
			facet_grid(~COUP_SC) +
			labs(x='\ncommunity type', y='phylogenetically unlinked couples\n(proportion)\n')
	ggsave(file=file.path(dir, paste0(run,'-phsc-extracoupletrm_couples_all_by_commtype_MonteCarlo.pdf')), w=5, h=6)
	
	#
	#	repeat for comparison using distance only, ie no topological information
	#
	rpd		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_DI') & TYPE=='distant')
	rpd		<- merge(rp, rpd, by=c('FEMALE_SANGER_ID','MALE_SANGER_ID','COUP_SC'))
	rpd		<- unique(subset(rpd, MALE_HH_NUM==FEMALE_HH_NUM), by='COUPID')
	rpds	<- subset(rpd, select=c(COUPID, COUP_SC, MALE_COMM_TYPE, POSTERIOR_ALPHA, POSTERIOR_BETA))
	rpds[, COUP_SC:='all +/+ couples']
	tmp2	<- subset(rpd, COUP_SC%in%c('seroinc','M->F','F->M'), select=c(COUPID, COUP_SC, MALE_COMM_TYPE, POSTERIOR_ALPHA, POSTERIOR_BETA))
	tmp2[, COUP_SC:='-/-, -/+, +/- at enrollment']
	rpds	<- rbind(rpds, tmp2)
	mc.it	<- 1e4	
	set.seed(42)
	rpds	<- do.call('rbind',lapply(seq_len(mc.it), function(i)
					{
						#	one Monte Carlo iteration:
						#	1- draw Z from posterior probabilities
						#	2- count successes Z=1 (unlinked couples) and evaluate binomial prob and confidence intervals
						#i		<- 1
						if(i%%1e2==0) cat('\n',i)
						mc		<- rpds[, list(	MC_IT=i, 
										MC_Z= rbinom(1,1,rbeta(1, POSTERIOR_ALPHA, POSTERIOR_BETA))), by=c('COUPID','COUP_SC','MALE_COMM_TYPE')]
						tmp		<- mc[, 	{
									z	<- binconf( length(which(MC_Z==1)), length(MC_Z) )				
									list(K=length(which(MC_Z==1)), N=length(MC_Z), P=z[1], QL=z[2], QU=z[3])
								}, by=c('MALE_COMM_TYPE','COUP_SC','MC_IT')]	
						tmp						
					}))
	save(rpds, file=file.path(dir, paste0(run,'-phsc-extracoupletrm_couples_all_by_commtype_MonteCarlo_usingdistanceonly.rda')))
	tmp		<- rpds[, list(N=mean(N), P=mean(P), QL=mean(QL), QU=mean(QU)), by=c('COUP_SC','MALE_COMM_TYPE')]
	set(tmp, NULL, 'MALE_COMM_TYPE', tmp[, factor(MALE_COMM_TYPE, levels=c('agrarian','trading','fisherfolk'))])
	set(tmp, tmp[, which(COUP_SC=='-/-, -/+, +/- at enrollment')],'COUP_SC','+/+ couples,\n-/-, -/+, +/- at enrollment')	
	#                                     COUP_SC MALE_COMM_TYPE   N         P         QL        QU
	#1:                           all +/+ couples       agrarian  53 0.1810358 0.10082490 0.3043593
	#2:                           all +/+ couples     fisherfolk 146 0.3590110 0.28582459 0.4394263
	#3:                           all +/+ couples        trading   4 0.4664500 0.13280196 0.8288349
	#4: +/+ couples,\n-/-, -/+, +/- at enrollment       agrarian  18 0.1656944 0.06000133 0.3885805
	#5: +/+ couples,\n-/-, -/+, +/- at enrollment     fisherfolk  31 0.3395161 0.19990158 0.5145191
	#6: +/+ couples,\n-/-, -/+, +/- at enrollment        trading   2 0.4430500 0.03326284 0.9277401
	
}

RakaiCouples.extracoupletrm.170106.analyze.extratrm.couples.communities<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)	
	require(coin)
	require(Hmisc)
	require(gamlss)
	
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]	
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	rc		<- copy(rp)
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments_allpairs.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	# load pairwise probabilities
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmp_allpairs_posteriors.rda')	
	#	make pair data table
	rp		<- copy(rpw)
	set(rp, NULL, c('DIR','FILE','PTY_RUN','RUN','W_FROM','W_TO','TYPE_RAW','TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3','PATRISTIC_DISTANCE','CONTIGUOUS','PATHS_12','PATHS_21','MALE_SANGER_ID_L','MALE_SANGER_ID_R','FEMALE_SANGER_ID_L','FEMALE_SANGER_ID_R'), NULL)
	rp		<- unique(rp)
	#	make COUPID
	rp[, COUPID:= paste0(MALE_RID,':',FEMALE_RID)]	
	#	add PAIR_TYPE
	tmp		<- unique(subset(rc, select=c(COUPID, MALE_HH_NUM, FEMALE_HH_NUM, COUP_SC, PAIR_TYPE)))	
	setnames(tmp, 'COUP_SC', 'COUP_TYPE')
	rp	<- merge(rp, tmp, by=c('COUPID','MALE_HH_NUM','FEMALE_HH_NUM'),all.x=1)
	set(rp, rp[, which(is.na(PAIR_TYPE))], 'PAIR_TYPE', 'not registered as couple')
	
	#
	#	select stable couples
	#	
	rp	<- subset(rp, PAIR_TYPE!='not registered as couple')
	
	#
	#	compare distance / distance + topology
	#	results: 	with TYPE_PAIR_TODI3 pairs in close transmission chains end up in unlinked 
	#				contiguous cannot be interpreted as 'A->C->B' it can also be 'A->B/C' or 'A->B and A->C'
	confint	<- 0.5
	rpd		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_DI') & TYPE=='distant')
	tmp		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='with intermediate\nor distant')
	rpd		<- rbind(rpd, tmp)
	rpd[, SELECT_P50:= as.character(factor(pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confint, levels=c(TRUE,FALSE), labels=c('Y','N')))]
	rpd		<- merge(rp, rpd, by=c('FEMALE_SANGER_ID','MALE_SANGER_ID'))
	rpd		<- dcast.data.table(rpd, COUPID+COUP_TYPE+LABEL+PTY_RUN+FEMALE_SANGER_ID+MALE_SANGER_ID+MALE_RID+FEMALE_RID~GROUP, value.var='SELECT_P50')
		
	tmp		<- subset(rpd, TYPE_PAIR_DI=='N' & TYPE_PAIR_TODI3=='Y')
	unique(tmp, by=c('MALE_RID','FEMALE_RID'))
	
	#
	#	select unlinked couples	by distance + topology
	#

	rpd		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_DI') & TYPE=='distant')		
	confint	<- 0.5
	rpd[, SELECT_P50:= as.character(factor(pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confint, levels=c(TRUE,FALSE), labels=c('Y','N')))]
	rpd[, SELECT_P80:= as.character(factor(pbeta(0.8, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confint, levels=c(TRUE,FALSE), labels=c('Y','N')))]
	rpd		<- subset(rpd, select=c(PTY_RUN, FEMALE_SANGER_ID, MALE_SANGER_ID, NEFF, KEFF, POSTERIOR_ALPHA, POSTERIOR_BETA, LABEL_SH, LABEL, SELECT_P50, SELECT_P80))
	rpd		<- merge(rp, rpd, by=c('FEMALE_SANGER_ID','MALE_SANGER_ID'), all.x=1)
	
	#
	#	region among cohabiting couples and extra-couple transmission
	#
	tmp	<- unique(subset(rpd, MALE_HH_NUM==FEMALE_HH_NUM & SELECT_P50%in%c('Y','N')), by='COUPID')
	tmp[, {
				z	<- binconf( length(which(SELECT_P50=='Y')), length(which(SELECT_P50%in%c('Y','N'))) )				
				list(K=length(which(SELECT_P50=='Y')), N=length(which(SELECT_P50%in%c('Y','N'))), P=z[1], QL=z[2], QU=z[3])
			}, by='MALE_REGION']	
	#REGION          K  N         P         QL        QU
	#1:           4  2 14 0.1428571 0.04009392 0.3994138
	#2:          12 29 93 0.3118280 0.22672869 0.4118559
	#3:           9  3 17 0.1764706 0.06191127 0.4102946
	#4:           5  3  8 0.3750000 0.13684429 0.6942576
	#5:          13 30 72 0.4166667 0.30985217 0.5319230
	#
	#	extra-couple trm seems higher in 5, 12, 13
	#	if I square the maps in Chang 2016 with the IDs in Collinson-Streng 2009
	#	4, 9 seem agragrian
	
	#
	#	community type among cohabiting couples and extra-couple transmission
	#
	#	results:
	#	significantly more phylog discordance among couples in fisherfolk communities
	tmp	<- unique(subset(rpd, PAIR_TYPE=='stable cohabiting' & SELECT_P50%in%c('Y','N')), by='COUPID')
	tmp2	<- tmp[, {
				z	<- binconf( length(which(SELECT_P50=='Y')), length(which(SELECT_P50%in%c('Y','N'))) )				
				list(K=length(which(SELECT_P50=='Y')), N=length(which(SELECT_P50%in%c('Y','N'))), P=z[1], QL=z[2], QU=z[3])
			}, by=c('MALE_COMM_TYPE')]
	setkey(tmp2, MALE_COMM_TYPE)	
	set(tmp2, NULL, 'MALE_COMM_TYPE', tmp2[, factor(MALE_COMM_TYPE, levels=c('agrarian','trading','fisherfolk'))])
	#   MALE_COMM_TYPE  K   N         P        QL        QU
	#1:       agrarian  9  53 0.1698113 0.09199945 0.2922528
	#2:     fisherfolk 56 147 0.3809524 0.30642783 0.4615405
	#3:        trading  2   4 0.5000000 0.15003899 0.8499610
	ggplot(subset(tmp2, MALE_COMM_TYPE!='trading'), aes(x=MALE_COMM_TYPE)) +
			geom_point(aes(y=P)) + geom_errorbar(aes(ymax=QU, ymin=QL), width=0.5) +
			geom_text(aes(y=QU+.05, label=paste0('n=',N))) + 
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous(limits=c(0,1), labels=percent, expand=c(0,0), breaks=seq(0,1,0.2)) +
			labs(x='\ncommunity type', y='phylogenetically unlinked couples among +/+ couples\n(proportion)\n')
	ggsave(file=file.path(dir, paste0(run,'-phsc-extracoupletrm_couples_all_by_commtype.pdf')), w=3, h=6)
	set(tmp, NULL,'SELECT_P50', tmp[, as.integer(as.character(factor(SELECT_P50, levels=c('Y','N'), labels=c('1','0'))))])	
	summary(gamlss(data=subset(tmp, select=c(SELECT_P50, MALE_COMM_TYPE)), SELECT_P50~MALE_COMM_TYPE, family=LO))
	#(Intercept)               0.10632    0.06625   1.605  0.11008   
	#MALE_COMM_TYPEfisherfolk  0.21257    0.07727   2.751  0.00649 **
	#MALE_COMM_TYPEtrading     0.39368    0.25008   1.574  0.11700 
	#
	#	repeat community type among cohabiting seroinc, M->F, F->M couples and extra-couple transmission
	#
	#	results:
	#	significantly more extra-couple trms among couples in fisherfolk communities
	tmp	<- unique(subset(rpd, COUP_TYPE%in%c('seroinc','M->F','F->M') & PAIR_TYPE=='stable cohabiting' & SELECT_P50%in%c('Y','N')), by='COUPID')
	tmp2	<- tmp[, {
				z	<- binconf( length(which(SELECT_P50=='Y')), length(which(SELECT_P50%in%c('Y','N'))) )				
				list(K=length(which(SELECT_P50=='Y')), N=length(which(SELECT_P50%in%c('Y','N'))), P=z[1], QL=z[2], QU=z[3])
			}, by=c('MALE_COMM_TYPE')]
	setkey(tmp2, MALE_COMM_TYPE)	
	set(tmp2, NULL, 'MALE_COMM_TYPE', tmp2[, factor(MALE_COMM_TYPE, levels=c('agrarian','trading','fisherfolk'))])
	#      MALE_COMM_TYPE  	K  N         P         QL        QU
	#1:       agrarian  2 18 0.1111111 0.03101952 0.3279977
	#2:     fisherfolk 10 31 0.3225806 0.18569427 0.4985899
	#3:        trading  1  2 0.5000000 0.02564665 0.9743534
	ggplot(subset(tmp2, MALE_COMM_TYPE!='trading'), aes(x=MALE_COMM_TYPE)) +
			geom_point(aes(y=P)) + geom_errorbar(aes(ymax=QU, ymin=QL), width=0.5) +
			geom_text(aes(y=QU+.05, label=paste0('n=',N))) + 
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous(limits=c(0,1), labels=percent, expand=c(0,0), breaks=seq(0,1,0.2)) +
			labs(x='\ncommunity type', y='extra-couple transmission among incident/initially discordant couples\n(posterior probability)\n')
	ggsave(file=file.path(dir, paste0(run,'-phsc-extracoupletrm_couples_incdisc_by_commtype.pdf')), w=3, h=6)
	set(tmp, NULL,'SELECT_P50', tmp[, as.integer(as.character(factor(SELECT_P50, levels=c('Y','N'), labels=c('1','0'))))])	
	summary(gamlss(data=subset(tmp, select=c(SELECT_P50, MALE_COMM_TYPE)), SELECT_P50~MALE_COMM_TYPE, family=LO))
	#                         Estimate Std. Error t value Pr(>|t|)  
	#(Intercept)               0.05694    0.09610   0.592    0.556
	#MALE_COMM_TYPEfisherfolk  0.16674    0.12082   1.380    0.174
	#MALE_COMM_TYPEtrading     0.44306    0.30391   1.458    0.151
	

	#
	#	How many individuals report extra partner among unlinked couples?
	#

	tmp	<- unique(subset(rpd, PAIR_TYPE=='stable cohabiting' & SELECT_P50%in%c('Y')), by='COUPID')
	tmp	<- melt(tmp, measure.vars=c('MALE_SEXC','FEMALE_SEXC'))
	set(tmp, NULL, 'variable', tmp[, factor(variable, levels=c('MALE_SEXC','FEMALE_SEXC'), labels=c('male','female') )])
	tmp[, value2:= 'No']
	set(tmp, tmp[, which(grepl('other partner',value))], 'value2', 'Yes')
	
	tmp2	<- tmp[, {
				z	<- binconf( length(which(value2=='Yes')), length(value2) )				
				list(K=length(which(value2=='Yes')), N=length(value2), P=z[1], QL=z[2], QU=z[3])
			}, by=c('MALE_COMM_TYPE','variable')]
	setkey(tmp2, MALE_COMM_TYPE)	
	#   MALE_COMM_TYPE variable  K  N         P          QL        QU
	#1:       agrarian     male  1  9 0.1111111 0.005699255 0.4349997
	#2:       agrarian   female  1  9 0.1111111 0.005699255 0.4349997
	#3:     fisherfolk     male 28 56 0.5000000 0.373317388 0.6266826
	#4:     fisherfolk   female 15 56 0.2678571 0.169573054 0.3959456
	#5:        trading     male  0  2 0.0000000 0.000000000 0.6576198
	#6:        trading   female  1  2 0.5000000 0.025646647 0.9743534
	ggplot(subset(tmp2, MALE_COMM_TYPE!='trading'), aes(x= variable, y=P, ymin=QL, ymax=QU)) + 
		geom_point() + geom_errorbar() + 
		facet_grid(~MALE_COMM_TYPE) +
		scale_y_continuous(expand=c(0,0), limits=c(0,1), labels=percent, breaks=seq(0,1,0.2)) +
		theme_bw() + 
		labs(x='\nindividuals in phylogenetically unlinked\nstable relationships', y='extra-marital partner reported\n(proportion)\n')
	ggsave(file=file.path(dir, paste0(run,'-phsc-extracoupletrm_couples_all_extrapartnerreported_conf.pdf')), w=4, h=6)
	ggplot(subset(tmp, MALE_COMM_TYPE!='trading'), aes(x= variable, fill=value2)) + 
			geom_bar(position='stack') + 
			facet_grid(~MALE_COMM_TYPE) +
			scale_y_continuous(expand=c(0,0), limits=c(0,65)) +
			theme_bw() + 
			labs(x='\npartner', y='phylogenetically unlinked +/+ couples\n', fill='extra-marital partner\nreported')
	ggsave(file=file.path(dir, paste0(run,'-phsc-extracoupletrm_couples_all_extrapartnerreported.pdf')), w=6, h=8)
	
	#
	#	repeat community type among couples that were found +/+ 
	#
	#	results:
	#	still significantly more extra-couple trms among couples in fisherfolk communities
	#	but discrepancy slightly lower, suggesting that both individuals were already infected before couple formation
	tmp	<- unique(subset(rpd, COUP_SC%in%c('seropos') & MALE_HH_NUM==FEMALE_HH_NUM & SELECT_P50%in%c('Y','N')), by='COUPID')
	tmp2	<- tmp[, {
				z	<- binconf( length(which(SELECT_P50=='Y')), length(which(SELECT_P50%in%c('Y','N'))) )				
				list(K=length(which(SELECT_P50=='Y')), N=length(which(SELECT_P50%in%c('Y','N'))), P=z[1], QL=z[2], QU=z[3])
			}, by=c('MALE_COMM_TYPE')]
	setkey(tmp2, MALE_COMM_TYPE)	
	set(tmp2, NULL, 'MALE_COMM_TYPE', tmp2[, factor(MALE_COMM_TYPE, levels=c('agrarian','trading','fisherfolk'))])
	#      MALE_COMM_TYPE  K   N         P         QL        QU
	#1:       agrarian  6  35 0.1714286 0.08102640 0.3268228
	#2:     fisherfolk 42 115 0.3652174 0.28289758 0.4562507
	#3:        trading  1   2 0.5000000 0.02564665 0.9743534
	ggplot(tmp2, aes(x=MALE_COMM_TYPE)) +
			geom_point(aes(y=P)) + geom_errorbar(aes(ymax=QU, ymin=QL), width=0.5) +
			geom_text(aes(y=QU+.05, label=paste0('n=',N))) + 
			theme_bw() + theme(legend.position='bottom') +
			scale_y_continuous(limits=c(0,1), labels=percent, expand=c(0,0), breaks=seq(0,1,0.2)) +
			labs(x='\ncommunity type', y='phylogenetically unlinked couples among couples that were +/+ at enrollment\n(proportion)\n')
	ggsave(file=file.path(dir, paste0(run,'-phsc-extracoupletrm_couples_seropos_by_commtype.pdf')), w=4, h=6)
	#
	#	community among cohabiting couples and extra-couple transmission
	#
	tmp	<- unique(subset(rpd, MALE_HH_NUM==FEMALE_HH_NUM & SELECT_P50%in%c('Y','N')), by='COUPID')	
	tmp	<- tmp[, {
				z	<- binconf( length(which(SELECT_P50=='Y')), length(which(SELECT_P50%in%c('Y','N'))) )				
				list(K=length(which(SELECT_P50=='Y')), N=length(which(SELECT_P50%in%c('Y','N'))), P=z[1], QL=z[2], QU=z[3])
			}, by=c('MALE_REGION','MALE_COMM_NUM')]
	setkey(tmp, MALE_REGION, MALE_COMM_NUM)
	#	numbers are small but 771 sticks out!
	#	38 is also fairly high and the largest fishing community in Rakai
	#	770 771 are fishing villages in Masaka district
	#	on the other hand, 770 seems fairly low, hmmm
	
	#	MALE_REGION MALE_COMM_NUM  K  N         P         QL        QU
	#1:           4             7  0  4 0.0000000 0.00000000 0.4898908
	#2:           4             8  0  1 0.0000000 0.00000000 0.9487067
	#3:           4           106  2  9 0.2222222 0.06322511 0.5474110
	#4:           5             2  1  4 0.2500000 0.01282332 0.6993582
	#5:           5            16  0  1 0.0000000 0.00000000 0.9487067
	#6:           5            33  1  2 0.5000000 0.02564665 0.9743534
	#7:           5           107  1  1 1.0000000 0.05129329 1.0000000
	#8:           9            34  1  5 0.2000000 0.01025866 0.6244654
	#9:           9            40  0  6 0.0000000 0.00000000 0.3903343
	#10:           9            56  0  1 0.0000000 0.00000000 0.9487067
	#11:           9            89  1  1 1.0000000 0.05129329 1.0000000
	#12:           9            94  0  1 0.0000000 0.00000000 0.9487067
	#13:           9           108  1  3 0.3333333 0.01709776 0.7923404
	#14:          12            23  0  9 0.0000000 0.00000000 0.2991450
	#15:          12            38 28 79 0.3544304 0.25795398 0.4644073
	#16:          12           370  2  5 0.4000000 0.11762077 0.7692757
	#17:          13           770  3 19 0.1578947 0.05520472 0.3756548
	#18:          13           771 20 37 0.5405405 0.38384047 0.6896143
	#19:          13           772  0  3 0.0000000 0.00000000 0.5614970
	#20:          13           773  0  1 0.0000000 0.00000000 0.9487067
	#21:          13           774  4 11 0.3636364 0.15166471 0.6461988
	
	
	tmp		<- unique(subset(rpd, 	MALE_HH_NUM==FEMALE_HH_NUM & 
							COUP_SC%in%c('seroinc','M->F','F->M') &
							!is.na(FEMALE_BIRTHDATE) & !is.na(MALE_BIRTHDATE) & 
							SELECT_P50%in%c('Y','N')), by='COUPID')
	tmp[, MALE_FIRSTPOSDATE_AGE:= MALE_FIRSTPOSDATE-MALE_BIRTHDATE]
	tmp[, FEMALE_FIRSTPOSDATE_AGE:= FEMALE_FIRSTPOSDATE-FEMALE_BIRTHDATE]
	tmp		<- merge(tmp, tmp[, list(	MALE_BOTHPOS_AGE= max(MALE_FIRSTPOSDATE, FEMALE_FIRSTPOSDATE)-MALE_BIRTHDATE, 
							FEMALE_BOTHPOS_AGE= max(MALE_FIRSTPOSDATE, FEMALE_FIRSTPOSDATE)-FEMALE_BIRTHDATE,
							BOTHPOS_AVGAGE= max(MALE_FIRSTPOSDATE, FEMALE_FIRSTPOSDATE)-(FEMALE_BIRTHDATE+FEMALE_BIRTHDATE)/2), by='COUPID'], by='COUPID')	
	tmp[, FIRSTPOSDATE_AGE_MFDIFF:= MALE_BIRTHDATE-FEMALE_BIRTHDATE]
	set(tmp, NULL,'SELECT_P50', tmp[, as.integer(as.character(factor(SELECT_P50, levels=c('Y','N'), labels=c('1','0'))))])	
	tmp2	<- subset(tmp, select=c(SELECT_P50, MALE_COMM_TYPE, BOTHPOS_AVGAGE, FEMALE_BOTHPOS_AGE, MALE_BOTHPOS_AGE, FIRSTPOSDATE_AGE_MFDIFF))
	summary(gamlss(data=tmp2, SELECT_P50~MALE_COMM_TYPE, family=LO))
	summary(gamlss(data=tmp2, SELECT_P50~MALE_COMM_TYPE+BOTHPOS_AVGAGE, family=LO))
	summary(gamlss(data=tmp2, SELECT_P50~MALE_COMM_TYPE+BOTHPOS_AVGAGE+FIRSTPOSDATE_AGE_MFDIFF, family=LO))	
}

RakaiCouples.extracoupletrm.170106.analyze.extratrm.couples.noncohabiting<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)	
	require(coin)
	require(Hmisc)
	require(gamlss)
	# load pty.runs
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]	
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	# load posterior relationships
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_extracoupletrm_bbmodels.rda')
	# add agrarian etc
	# 	from Kate:
	#	Trading communities: 1   4   16  22  24  33  51  107 776
	#	Fishing communities: 38, 770, 771, 774	
	# 	tmp		<- rp[, sort(unique(c(MALE_COMM_NUM,FEMALE_COMM_NUM)))]; cat(paste('"', tmp, '",', collapse='', sep=''))
	tmp		<- data.table(	COMM_NUM=	c("1","2","4","7","8","16","22","23","24","33","34","38","40","51","56","57","58","89","94","106","107","108","370","391","770","771","772","773","774","776"),
			COMM_TYPE=	c("T","A","T","A","A", "T", "T", "A", "T", "T", "A", "F", "A", "T", "A", "A", "A", "A", "A",  "A",  "T",  "A",  "A",  "A", "F",  "F",  "A",  "A",  "F",   "T"))
	set(tmp, NULL, 'COMM_TYPE', tmp[, as.character(factor(COMM_TYPE, levels=c('A','T','F'), labels=c('agrarian','trading','fisherfolk')))])
	set(tmp, NULL, 'COMM_NUM', tmp[, as.integer(COMM_NUM)])
	setnames(tmp, c('COMM_NUM','COMM_TYPE'), c('MALE_COMM_NUM','MALE_COMM_TYPE'))
	rp		<- merge(rp, tmp, by='MALE_COMM_NUM')
	setnames(tmp, c('MALE_COMM_NUM','MALE_COMM_TYPE'), c('FEMALE_COMM_NUM','FEMALE_COMM_TYPE'))
	rp		<- merge(rp, tmp, by='FEMALE_COMM_NUM')
	#
	#	select couples	
	#
	rpd		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='with intermediate\nor distant')	
	#	make central selection: post prob of extra-couple trm > 50%
	rpd[, SELECT_P50:= as.character(factor(pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>0.8, levels=c(TRUE,FALSE), labels=c('Y','N')))]
	rpd[, SELECT_P80:= as.character(factor(pbeta(0.8, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>0.8, levels=c(TRUE,FALSE), labels=c('Y','N')))]
	rpd		<- subset(rpd, select=c(PTY_RUN, FEMALE_SANGER_ID, MALE_SANGER_ID, NEFF, KEFF, POSTERIOR_ALPHA, POSTERIOR_BETA, LABEL_SH, LABEL, SELECT_P50, SELECT_P80))
	rpd		<- merge(rp, rpd, by=c('FEMALE_SANGER_ID','MALE_SANGER_ID'), all.x=1)
	set(rpd, rpd[, which(is.na(SELECT_P50))], 'SELECT_P50', 'No_PANGEA_data')
	set(rpd, rpd[, which(is.na(SELECT_P80))], 'SELECT_P80', 'No_PANGEA_data')
	#
	unique(rpd, by='COUPID')[, table(SELECT_P50)]
	#	308 couples, 211 with PANGEA data, 73 with strong evidence for extra-couple trm 
	#
	#	cohabitation and extra-couple transmission
	#
	unique(subset(rpd, MALE_HH_NUM==FEMALE_HH_NUM), by='COUPID')[, table(SELECT_P50)]
	binconf(64,199)
	# 291 couples cohabitate
	# 64 / 199 (32%) with PANGEA data extra-couple trms (at least)
	# conf int:	0.2606181-0.3893548
	unique(subset(rpd, MALE_HH_NUM!=FEMALE_HH_NUM), by='COUPID')[, table(SELECT_P50)]
	unique(subset(rpd, MALE_HH_NUM!=FEMALE_HH_NUM & !is.na(FEMALE_TAXA) & !is.na(MALE_TAXA) & !is.na(NEFF)), by='COUPID')[, COUPID]
	#	"A114712:H104287" "H115099:C106054" "F105662:B096376" "C114447:J114639" "F049477:F115144" "H061328:B104087" "B110153:K105200" "E104488:K105535" "H087552:A032173" "G107639:C107790" "F108382:A107756" "A104696:D107343"
	
	# 17 couples do not cohabitate in same household
	# 9 / 12 (75%) with PANGEA data extra couple trms
	# conf int: 0.4676947-0.9110583
	unique(subset(rpd, MALE_COMM_NUM!=FEMALE_COMM_NUM), by='COUPID')	
	# 3 couples do not live in same community
	# 2/3 extra couple trms
	unique(subset(rpd, MALE_REGION!=FEMALE_REGION), by='COUPID')		
	# 1 couples does not live in same region
	# 1/1 extra couple trms
	#
	#	cohabitation / extra-couple transmission significant ?
	#	
	rpd[, COHABITATE:= as.character(factor(MALE_HH_NUM==FEMALE_HH_NUM, levels=c(TRUE,FALSE),labels=c('Y','N')))]
	chisq_test(factor(SELECT_P50) ~ factor(COHABITATE), data=subset(rpd, SELECT_P50%in%c('Y','N')), distribution="exact")
	# p-value = 0.001212
}

RakaiCouples.extracoupletrm.170106.analyze.extratrm.couples.indexpartner<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)	
	require(coin)
	require(Hmisc)
	require(gamlss)
	# load pty.runs
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]	
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	# load posterior relationships
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_extracoupletrm_bbmodels.rda')
	# add agrarian etc
	# 	from Kate:
	#	Trading communities: 1   4   16  22  24  33  51  107 776
	#	Fishing communities: 38, 770, 771, 774	
	# 	tmp		<- rp[, sort(unique(c(MALE_COMM_NUM,FEMALE_COMM_NUM)))]; cat(paste('"', tmp, '",', collapse='', sep=''))
	tmp		<- data.table(	COMM_NUM=	c("1","2","4","7","8","16","22","23","24","33","34","38","40","51","56","57","58","89","94","106","107","108","370","391","770","771","772","773","774","776"),
			COMM_TYPE=	c("T","A","T","A","A", "T", "T", "A", "T", "T", "A", "F", "A", "T", "A", "A", "A", "A", "A",  "A",  "T",  "A",  "A",  "A", "F",  "F",  "A",  "A",  "F",   "T"))
	set(tmp, NULL, 'COMM_TYPE', tmp[, as.character(factor(COMM_TYPE, levels=c('A','T','F'), labels=c('agrarian','trading','fisherfolk')))])
	set(tmp, NULL, 'COMM_NUM', tmp[, as.integer(COMM_NUM)])
	setnames(tmp, c('COMM_NUM','COMM_TYPE'), c('MALE_COMM_NUM','MALE_COMM_TYPE'))
	rp		<- merge(rp, tmp, by='MALE_COMM_NUM')
	setnames(tmp, c('MALE_COMM_NUM','MALE_COMM_TYPE'), c('FEMALE_COMM_NUM','FEMALE_COMM_TYPE'))
	rp		<- merge(rp, tmp, by='FEMALE_COMM_NUM')
	#
	#	select couples	
	#
	rpd		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='with intermediate\nor distant')	
	#	make central selection: post prob of extra-couple trm > 50%
	rpd[, SELECT_P50:= as.character(factor(pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>0.8, levels=c(TRUE,FALSE), labels=c('Y','N')))]
	rpd[, SELECT_P80:= as.character(factor(pbeta(0.8, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>0.8, levels=c(TRUE,FALSE), labels=c('Y','N')))]
	rpd		<- subset(rpd, select=c(PTY_RUN, FEMALE_SANGER_ID, MALE_SANGER_ID, NEFF, KEFF, POSTERIOR_ALPHA, POSTERIOR_BETA, LABEL_SH, LABEL, SELECT_P50, SELECT_P80))
	rpd		<- merge(rp, rpd, by=c('FEMALE_SANGER_ID','MALE_SANGER_ID'), all.x=1)
	set(rpd, rpd[, which(is.na(SELECT_P50))], 'SELECT_P50', 'No_PANGEA_data')
	set(rpd, rpd[, which(is.na(SELECT_P80))], 'SELECT_P80', 'No_PANGEA_data')
	#
	#	by index partner
	#
	tmp	<- unique(subset(rpd, COUP_SC%in%c('M->F','F->M') & MALE_HH_NUM==FEMALE_HH_NUM & SELECT_P50%in%c('Y','N')), by='COUPID')
	tmp2	<- tmp[, {
				z	<- binconf( length(which(SELECT_P50=='Y')), length(which(SELECT_P50%in%c('Y','N'))) )				
				list(K=length(which(SELECT_P50=='Y')), N=length(which(SELECT_P50%in%c('Y','N'))), P=z[1], QL=z[2], QU=z[3])
			}, by=c('COUP_SC')]
	#   COUP_SC K  N         P        QL        QU
	#1:    M->F 4 14 0.2857143 0.1172138 0.5464908
	#2:    F->M 3  6 0.5000000 0.1876163 0.8123837
	set(tmp, NULL,'SELECT_P50', tmp[, as.integer(as.character(factor(SELECT_P50, levels=c('Y','N'), labels=c('1','0'))))])
	summary(gamlss(data=subset(tmp, select=c(SELECT_P50, COUP_SC)), SELECT_P50~COUP_SC, family=LO))	
}

RakaiCouples.extracoupletrm.170106.analyze.extratrm.couples.samplingbias<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)	
	require(coin)
	require(Hmisc)
	require(gamlss)
	# load pty.runs
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]	
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	# load posterior relationships
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_extracoupletrm_bbmodels.rda')
	# add agrarian etc
	# 	from Kate:
	#	Trading communities: 1   4   16  22  24  33  51  107 776
	#	Fishing communities: 38, 770, 771, 774	
	# 	tmp		<- rp[, sort(unique(c(MALE_COMM_NUM,FEMALE_COMM_NUM)))]; cat(paste('"', tmp, '",', collapse='', sep=''))
	tmp		<- data.table(	COMM_NUM=	c("1","2","4","7","8","16","22","23","24","33","34","38","40","51","56","57","58","89","94","106","107","108","370","391","770","771","772","773","774","776"),
			COMM_TYPE=	c("T","A","T","A","A", "T", "T", "A", "T", "T", "A", "F", "A", "T", "A", "A", "A", "A", "A",  "A",  "T",  "A",  "A",  "A", "F",  "F",  "A",  "A",  "F",   "T"))
	set(tmp, NULL, 'COMM_TYPE', tmp[, as.character(factor(COMM_TYPE, levels=c('A','T','F'), labels=c('agrarian','trading','fisherfolk')))])
	set(tmp, NULL, 'COMM_NUM', tmp[, as.integer(COMM_NUM)])
	setnames(tmp, c('COMM_NUM','COMM_TYPE'), c('MALE_COMM_NUM','MALE_COMM_TYPE'))
	rp		<- merge(rp, tmp, by='MALE_COMM_NUM')
	setnames(tmp, c('MALE_COMM_NUM','MALE_COMM_TYPE'), c('FEMALE_COMM_NUM','FEMALE_COMM_TYPE'))
	rp		<- merge(rp, tmp, by='FEMALE_COMM_NUM')
	#
	#	sampling bias of couples 
	#	-- currently not meaningful
	sb		<- unique(rp, by='COUPID')
	sb[, COUPLE_SEQUENCED:= factor(!is.na(MALE_SANGER_ID)&!is.na(FEMALE_SANGER_ID), levels=c(TRUE,FALSE),labels=c('Y','N'))]	
	sb[, MALE_FIRSTPOSDATE_AGE:= MALE_FIRSTPOSDATE-MALE_BIRTHDATE]
	sb[, FEMALE_FIRSTPOSDATE_AGE:= FEMALE_FIRSTPOSDATE-FEMALE_BIRTHDATE]
	sb[, FIRSTPOSDATE_AGE_MFDIFF:= MALE_FIRSTPOSDATE_AGE-FEMALE_FIRSTPOSDATE_AGE]	
	tmp		<- subset(sb, MALE_HH_NUM==FEMALE_HH_NUM & !is.na(MALE_FIRSTPOSDATE_AGE) & !is.na(FEMALE_FIRSTPOSDATE_AGE))
	#	by agrarian/trading/fisherfolk	
	tmp[, {
				z	<- binconf( length(which(COUPLE_SEQUENCED=='Y')), length(which(COUPLE_SEQUENCED%in%c('Y','N'))) )				
				list(K=length(which(COUPLE_SEQUENCED=='Y')), N=length(which(COUPLE_SEQUENCED%in%c('Y','N'))), P=z[1], QL=z[2], QU=z[3])
			}, by='MALE_COMM_TYPE']
	#	by age
	tmp2	<- melt(tmp, measure.vars=c('MALE_FIRSTPOSDATE_AGE','FEMALE_FIRSTPOSDATE_AGE'))
	set(tmp2, NULL,'variable', tmp2[,factor(variable, 	levels=c('MALE_FIRSTPOSDATE_AGE','FEMALE_FIRSTPOSDATE_AGE'),labels=c('male partner','female partner'))])
	set(tmp2, NULL,'COUPLE_SEQUENCED', tmp2[, factor(COUPLE_SEQUENCED, levels=c('Y','N'), labels=c('Yes','No'))])
	ggplot(tmp2, aes(x=COUPLE_SEQUENCED, y=value)) +
			geom_boxplot() +
			facet_grid(~MALE_COMM_TYPE) 
}

RakaiCouples.extracoupletrm.170106.addposteriors<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)	
	# load pty.runs
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]	
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	#
	rpw[, table(RUN, useNA='if')]
	#	define plotting order: largest number of trm assignments	
	tmp		<- rpw[, list( WIN_TR=length(which(grepl('close|anc',TYPE_DIR_TODI7x3))) ), by=c('PTY_RUN','COUP_SC','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	tmp		<- tmp[order(-WIN_TR),]
	tmp[, PLOT_ID:=seq_len(nrow(tmp))]	
	#	define label
	tmp		<- merge(tmp, rp, by=c('COUP_SC','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	setkey(tmp, PLOT_ID)
	tmp[, LABEL_SH:= factor(PLOT_ID, levels=PLOT_ID, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair ', COUPID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL, LABEL_SH))
	rpw		<- merge(tmp, rpw, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	#	define min/max reads
	tmp		<- rpw[, list(	ID_R_MIN=min(ID1_R, ID2_R),
					ID_R_MAX=max(ID1_R, ID2_R)), by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM')]
	rpw		<- merge(rpw, tmp, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM'))
	#
	#	prepare colours
	#	
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
	#
	#	define topo types
	#
	rpw[, TYPE_PAIR_TO:= 'other']
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair_close','pair','pair_distant'))], 'TYPE_PAIR_TO', 'ancestral/\nintermingled')
	rpw[, TYPE_PAIR_TODI3:= 'with intermediate\nor distant']
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair_close','other_close'))], 'TYPE_PAIR_TODI3', 'no intermediate\n and close')
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair','other'))], 'TYPE_PAIR_TODI3', 'no intermediate\n but not close')
	rpw[, TYPE_PAIR_DI:= cut(PATRISTIC_DISTANCE, breaks=c(1e-12,0.02,0.05,2), labels=c('close','intermediate\ndistance','distant'))] 
	tmp		<- 	rpw[, list(PD_MEAN=mean(PATRISTIC_DISTANCE)), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	rpw[, TYPE_PAIR_TODI2x2:= 'not close other']
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair_close'))], 'TYPE_PAIR_TODI2x2', 'close ancestral/\nintermingled')
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('pair','pair_distant'))], 'TYPE_PAIR_TODI2x2', 'not close ancestral/\nintermingled')
	set(rpw, rpw[, which(TYPE_PAIR_TODI3x3%in%c('withintermediate_close','other_close'))], 'TYPE_PAIR_TODI2x2', 'close other')	
	set(rpw, NULL, 'TYPE_DIR_TODI7x3', rpw[, gsub('intermediate',' intermediate',gsub('intermediate_','intermediate\n',gsub('intermingled_','intermingled\n',gsub('(chain_[fm][mf])_','\\1\n',TYPE_DIR_TODI7x3))))])	
	#	define effectively independent number of windows
	#	select only W_FROM that are > W_TO	
	tmp		<- rpw[, {
				z		<- 1
				repeat
				{
					zz		<- which(W_FROM>max(W_TO[z]))[1]
					if(length(zz)==0 | is.na(zz))	break
					z		<- c(z, zz)			
				}
				list(W_FROM=W_FROM[z], W_TO=W_TO[z])
			}, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN')]
	tmp[, OVERLAP:= 0L]
	rpw		<- merge(rpw, tmp, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','W_FROM','W_TO','RUN'), all.x=1 )
	set(rpw, rpw[, which(is.na(OVERLAP))], 'OVERLAP', 1L)		
	rpw		<- melt(rpw, measure.vars=c('TYPE_PAIR_TO','TYPE_PAIR_TODI3x3','TYPE_PAIR_TODI3','TYPE_PAIR_DI','TYPE_PAIR_TODI2x2','TYPE_DIR_TODI7x3'), variable.name='GROUP', value.name='TYPE')
	set(rpw, NULL, 'TYPE', rpw[,gsub('_',' ', TYPE)])	
	#
	#	for each pair count windows by TYPE_PAIR_DI (k)
	#
	rplkl	<- rpw[, list(K=length(W_FROM)), by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','COUP_SC','LABEL','LABEL_SH','GROUP','TYPE')]	
	rplkl	<- dcast.data.table(rplkl, RUN+PTY_RUN+FEMALE_SANGER_ID+MALE_SANGER_ID+COUP_SC+LABEL+LABEL_SH~GROUP+TYPE, value.var='K')
	for(x in setdiff(colnames(rplkl),c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','COUP_SC','LABEL','LABEL_SH','GROUP')))
		set(rplkl, which(is.na(rplkl[[x]])), x, 0L)	
	for(x in colnames(rplkl)[grepl('TYPE_PD_MEAN',colnames(rplkl))])
		set(rplkl, which(rplkl[[x]]>0), x, 1L)	
	rplkl	<- melt(rplkl, id.vars=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','COUP_SC','LABEL','LABEL_SH'), variable.name='GROUP', value.name='K')
	rplkl[, TYPE:= gsub('.*_([^_]+)$','\\1',GROUP)]
	set(rplkl, NULL, 'GROUP', rplkl[, gsub('(.*)_[^_]+$','\\1',GROUP)])	
	#	for each pair count total windows (n) and effective windows (neff)	
	tmp		<- rpw[, list(N=length(W_FROM), NEFF=sum(1-OVERLAP)), by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN','GROUP')]
	rplkl	<- merge(rplkl, tmp, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN','GROUP'))
	#	define effective number of windows of type
	rplkl[, KEFF:= K/N*NEFF]	
	#
	#	add prior (equal probabilities to all types), posterior, and marginal probabilities
	#
	rplkl[, DIR_PRIOR:= 0.1]
	rplkl[, DIR_PO:= DIR_PRIOR+KEFF]
	tmp		<- rplkl[, {
							alpha	<- DIR_PO
							beta	<- sum(DIR_PO)-DIR_PO				
							list(	TYPE=TYPE, 
									POSTERIOR_ALPHA=alpha, 
									POSTERIOR_BETA=beta)	
						}, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN','GROUP')]
	rplkl	<- merge(rplkl, tmp, by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','GROUP','TYPE'))
	#
	#	plot detailed summary plots per couple
	#	
	groups		<- c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3','TYPE_PAIR_TODI2x2','TYPE_PAIR_TODI3','TYPE_PAIR_DI')
	groups		<- c('TYPE_PAIR_TODI3')
	for(group in groups)
	{
		widths	<- unit(c(4, 6), "null")
		heights	<- unit(c(2, 3.5, 4, 5), "null")
		height	<- 9
		if(group%in%c('TYPE_DIR_TODI7x3'))
		{
			widths	<- unit(c(4, 6), "null")
			heights	<- unit(c(2, 3.5, 4, 15), "null")
			height	<- 17
		}		
		if(group%in%c('TYPE_PAIR_TODI2x2'))
		{
			heights	<- unit(c(2, 3.5, 4, 3.75), "null")
			height	<- 8
		}
		if(group%in%c('TYPE_PAIR_TODI3','TYPE_PAIR_DI'))
		{
			heights	<- unit(c(2, 3.5, 4, 3.5), "null")
			height	<- 7
		}
		pdf(file=file.path(dir, paste(run,'-phsc-extracoupletrm_',group,'.pdf',sep='')), w=10, h=height)	
		plot.tmp	<- unique(subset(rplkl, GROUP==group, c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','LABEL')), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
		setkey(plot.tmp, LABEL)		
		for(i in seq_len(nrow(plot.tmp)))
		{
			#pty_run	<- 38; id1		<- '16016_1_4'; id2		<- '15105_1_35'
			pty_run	<- plot.tmp[i, PTY_RUN]; id1		<- plot.tmp[i, FEMALE_SANGER_ID]; id2		<- plot.tmp[i, MALE_SANGER_ID]		
			tmp		<- subset(rpw, PTY_RUN==pty_run & GROUP==group & FEMALE_SANGER_ID==id1 & MALE_SANGER_ID==id2)
			p1		<- ggplot(tmp, aes(x=W_FROM)) +			
					geom_bar(aes(y=ID_R_MAX, colour=TYPE), stat='identity', fill='transparent') +
					geom_bar(aes(y=ID_R_MIN, fill=TYPE), stat='identity', colour='transparent') +
					labs(x='', y='number of reads', fill='phylogenetic\nrelationship\n', colour='phylogenetic\nrelationship\n') +
					scale_fill_manual(values=cols.type[[group]]) +
					scale_colour_manual(values=cols.type[[group]]) +
					scale_x_continuous(breaks=seq(0,1e4,500), minor_breaks=seq(0,1e4,100), limits=c(rpw[, min(W_FROM)], rpw[, max(W_FROM)])) +
					scale_y_log10(breaks=c(10,100,1000,1e4,1e5)) +
					theme_bw() + theme(legend.position='left') +			
					guides(fill=FALSE, colour=FALSE)
			p2		<- ggplot(tmp, aes(x=W_FROM, y=PATRISTIC_DISTANCE)) +
					geom_point(size=1) +					
					labs(x='window start\n\n', y='patristic distance') +
					scale_x_continuous(breaks=seq(0,1e4,500), minor_breaks=seq(0,1e4,100), limits=c(rpw[, min(W_FROM)], rpw[, max(W_FROM)])) +
					scale_y_log10(labels=percent, limits=c(0.001, 0.7), expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
					theme_bw() + theme(legend.position='left')
			tmp		<- subset(rplkl, GROUP==group & PTY_RUN==pty_run & FEMALE_SANGER_ID==id1 & MALE_SANGER_ID==id2)
			p3		<- ggplot(tmp, aes(x=TYPE, y=KEFF, fill=TYPE)) + geom_bar(stat='identity') +
					scale_fill_manual(values=cols.type[[group]]) +
					theme_bw() + theme(legend.position='bottom') +
					coord_flip() + guides(fill=FALSE) +			
					labs(x='', y='\nnon-overlapping windows\n(number)', fill='phylogenetic\nrelationship\n')
			p4		<- ggplot(tmp, aes(x=TYPE, 	middle=qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA), 
												ymin=qbeta(0.025, POSTERIOR_ALPHA, POSTERIOR_BETA), 
												ymax=qbeta(0.975, POSTERIOR_ALPHA, POSTERIOR_BETA), 
												lower=qbeta(0.25, POSTERIOR_ALPHA, POSTERIOR_BETA), 
												upper=qbeta(0.75, POSTERIOR_ALPHA, POSTERIOR_BETA), fill=TYPE)) + 
					geom_boxplot(stat='identity') +
					scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), limits=c(0,1), expand=c(0,0)) +
					scale_fill_manual(values=cols.type[[group]]) +
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
	#
	#	plot 'TYPE_PAIR_TODI3': keff and posterior prob of extra couple trm
	#
	g_legend<-function(a.gplot)
	{
		tmp <- ggplot_gtable(ggplot_build(a.gplot))
		leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
		legend <- tmp$grobs[[leg]]
		legend
	}
	
	tmp		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_TODI3'))
	setkey(tmp, LABEL_SH)	
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=rev(c("no intermediate\n and close","no intermediate\n but not close","with intermediate\nor distant")))])	
	#tmp2	<- c(cols.type[['TYPE_PAIR_DI']],cols.type[['TYPE_PAIR_TODI3']],cols.type[['TYPE_PAIR_TODI2x2']])	
	p1		<- ggplot(tmp, aes(x=LABEL_SH, y=KEFF, fill=TYPE)) + 
			geom_bar(stat='identity',position='stack') +
			scale_y_continuous(expand=c(0,0)) +
			scale_fill_manual(values=cols.type[['TYPE_PAIR_TODI3']]) +			
			theme_bw() + 
			theme(	legend.position='bottom', 
					axis.ticks.y=element_blank()) +			
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE)) +
			labs(	x='', 
					y='non-overlapping windows with ancestral assignments\n(number)\n',
					fill='phylogenetic\nrelationship')
	p2		<- ggplot(	subset(tmp, TYPE%in%c("with intermediate\nor distant")), 
						aes(x=LABEL_SH, 	middle=qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA), 
										lower=qbeta(0.25, POSTERIOR_ALPHA, POSTERIOR_BETA), 
										upper=qbeta(0.75, POSTERIOR_ALPHA, POSTERIOR_BETA), 
										ymin=qbeta(0.025, POSTERIOR_ALPHA, POSTERIOR_BETA), 
										ymax=qbeta(0.975, POSTERIOR_ALPHA, POSTERIOR_BETA), fill=TYPE)) + 
			geom_boxplot(stat='identity') +
			scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.2), labels=percent) +
			scale_fill_manual(values=cols.type[['TYPE_PAIR_TODI3']]) +
			theme_bw() + 
			theme(	legend.position='bottom', axis.text.y=element_blank(),
					axis.ticks.y=element_blank()) +
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE)) +
			labs(	x='', 
					y='posterior probability\n\n',
					fill='phylogenetic\nrelationship')
	p3		<- g_legend(p1)	
	p3$vp	<- viewport(layout.pos.row=2, layout.pos.col=1:2)	
	pdf(file=file.path(dir, paste(run,'-phsc-extracoupletrms_relationships_TODI3.pdf',sep='')), w=15, h=60)
	#pdf(file=file.path(dir, paste(run,'-phsc-relationships_serodisc_compare_DI_to_TODI3.pdf',sep='')), w=12, h=15)
	grid.newpage()	
	pushViewport(viewport(layout = grid.layout(2, 2, heights=unit(c(40,1), "null"), widths=unit(c(7, 3), "null"))))   	
	print(p1+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
	print(p2+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
	grid.draw(p3)
	dev.off()	
	#
	#	plot 'TYPE_PAIR_TODI3': keff and posterior prob of no intermediate and close
	#	
	tmp		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_TODI3'))
	setkey(tmp, LABEL_SH)	
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=rev(c("no intermediate\n and close","no intermediate\n but not close","with intermediate\nor distant")))])	
	#tmp2	<- c(cols.type[['TYPE_PAIR_DI']],cols.type[['TYPE_PAIR_TODI3']],cols.type[['TYPE_PAIR_TODI2x2']])	
	p1		<- ggplot(tmp, aes(x=LABEL_SH, y=KEFF, fill=TYPE)) + 
			geom_bar(stat='identity',position='stack') +
			scale_y_continuous(expand=c(0,0)) +
			scale_fill_manual(values=cols.type[['TYPE_PAIR_TODI3']]) +			
			theme_bw() + 
			theme(	legend.position='bottom', 
					axis.ticks.y=element_blank()) +			
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE)) +
			labs(	x='', 
					y='non-overlapping windows with ancestral assignments\n(number)\n',
					fill='phylogenetic\nrelationship')
	p2		<- ggplot(	subset(tmp, TYPE%in%c("no intermediate\n and close")), 
					aes(x=LABEL_SH, 	middle=qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA), 
							lower=qbeta(0.25, POSTERIOR_ALPHA, POSTERIOR_BETA), 
							upper=qbeta(0.75, POSTERIOR_ALPHA, POSTERIOR_BETA), 
							ymin=qbeta(0.025, POSTERIOR_ALPHA, POSTERIOR_BETA), 
							ymax=qbeta(0.975, POSTERIOR_ALPHA, POSTERIOR_BETA), fill=TYPE)) + 
			geom_boxplot(stat='identity') +
			scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.2), labels=percent) +
			scale_fill_manual(values=cols.type[['TYPE_PAIR_TODI3']]) +
			theme_bw() + 
			theme(	legend.position='bottom', axis.text.y=element_blank(),
					axis.ticks.y=element_blank()) +
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE)) +
			labs(	x='', 
					y='posterior probability\n\n',
					fill='phylogenetic\nrelationship')
	p3		<- g_legend(p1)	
	p3$vp	<- viewport(layout.pos.row=2, layout.pos.col=1:2)	
	pdf(file=file.path(dir, paste(run,'-phsc-likelypair_relationships_TODI3.pdf',sep='')), w=15, h=60)	
	grid.newpage()	
	pushViewport(viewport(layout = grid.layout(2, 2, heights=unit(c(40,1), "null"), widths=unit(c(7, 3), "null"))))   	
	print(p1+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
	print(p2+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
	grid.draw(p3)
	dev.off()	
	#
	#	save
	save(rplkl, file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_extracoupletrm_bbmodels.rda')	
}

RakaiCouples.analyze.couples.161219.direction.serodisc<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	# load pairwise probabilities
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmp_bbmodels.rda')
	
	#
	#	select close ancestral/intermingled > 50%
	#
	
	rpd		<- subset(rplkl, COUP_SC%in%c('M->F','F->M') & GROUP%in%c('TYPE_PAIR_TODI2x2') & TYPE=='close ancestral/\nintermingled')
	rpd		<- subset(rpd, POSTERIOR_IL>0.5, c(RUN, PTY_RUN, FEMALE_SANGER_ID, MALE_SANGER_ID))
	rpd		<- merge(rpd, subset(rpw, COUP_SC%in%c('M->F','F->M')), by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID'))
	
	#	define plotting order: largest number of trm assignments	
	tmp		<- rpd[, list( WIN_TR=length(which(grepl('close|anc',TYPE_DIR_TODI7x3))) ), by=c('PTY_RUN','COUP_SC','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	tmp		<- tmp[order(-WIN_TR),]
	tmp[, PLOT_ID:=seq_len(nrow(tmp))]	
	#	define label
	tmp		<- merge(tmp, rp, by=c('COUP_SC','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	setkey(tmp, PLOT_ID)
	tmp[, LABEL_SH:= factor(PLOT_ID, levels=PLOT_ID, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair ', COUPID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL, LABEL_SH))
	rpd		<- merge(tmp, rpd, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	#	define min/max reads
	tmp		<- rpd[, list(	ID_R_MIN=min(ID1_R, ID2_R),
					ID_R_MAX=max(ID1_R, ID2_R)), by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM')]
	rpd		<- merge(rpd, tmp, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM'))	
	
	#
	#	define topo types
	#
	rpd[, TYPE_DIR_ANC:=NA_character_]
	set(rpd, rpd[, which(COUP_SC=='M->F' & grepl('chain_mf',TYPE_DIR_TODI7x3))], 'TYPE_DIR_ANC', 'chain_correct_direction')
	set(rpd, rpd[, which(COUP_SC=='F->M' & grepl('chain_fm',TYPE_DIR_TODI7x3))], 'TYPE_DIR_ANC', 'chain_correct_direction')
	set(rpd, rpd[, which(COUP_SC=='M->F' & grepl('chain_fm',TYPE_DIR_TODI7x3))], 'TYPE_DIR_ANC', 'chain_incorrect_direction')
	set(rpd, rpd[, which(COUP_SC=='F->M' & grepl('chain_mf',TYPE_DIR_TODI7x3))], 'TYPE_DIR_ANC', 'chain_incorrect_direction')
	rpd[, TYPE_DIR_NUMBER:='no chain']
	set(rpd, rpd[, which(grepl('chain',TYPE_DIR_TODI7x3))], 'TYPE_DIR_NUMBER', 'chain')
	set(rpw, NULL, 'TYPE_DIR_TODI7x3', rpw[, gsub('intermediate',' intermediate',gsub('intermediate_','intermediate\n',gsub('intermingled_','intermingled\n',gsub('(chain_[fm][mf])_','\\1\n',TYPE_DIR_TODI7x3))))])
	#
	#	define effectively independent number of windows
	#	select only W_FROM that are > W_TO
	#
	tmp		<- rpd[, {
				z		<- 1
				repeat
				{
					zz		<- which(W_FROM>max(W_TO[z]))[1]
					if(length(zz)==0 | is.na(zz))	break
					z		<- c(z, zz)			
				}
				list(W_FROM=W_FROM[z], W_TO=W_TO[z])
			}, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN')]
	tmp[, OVERLAP:= 0L]
	rpd		<- merge(rpd, tmp, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','W_FROM','W_TO','RUN'), all.x=1 )
	set(rpd, rpd[, which(is.na(OVERLAP))], 'OVERLAP', 1L)
	tmp		<- subset(rpd, !is.na(TYPE_DIR_ANC))[, {
				z		<- 1
				repeat
				{
					zz		<- which(W_FROM>max(W_TO[z]))[1]
					if(length(zz)==0 | is.na(zz))	break
					z		<- c(z, zz)			
				}
				list(W_FROM=W_FROM[z], W_TO=W_TO[z])
			}, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN')]
	tmp[, OVERLAP_DIR_ANC:= 0L]
	rpd		<- merge(rpd, tmp, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','W_FROM','W_TO','RUN'), all.x=1 )
	set(rpd, rpd[, which(!is.na(TYPE_DIR_ANC) & is.na(OVERLAP_DIR_ANC))], 'OVERLAP_DIR_ANC', 1L)
	#
	rpd		<- melt(rpd, measure.vars=c('TYPE_PAIR_TODI3x3','TYPE_DIR_ANC','TYPE_DIR_NUMBER'), variable.name='GROUP', value.name='TYPE')
	set(rpd, NULL, 'TYPE', rpd[,gsub('_',' ', TYPE)])
	#	remove windows with no chain assignments for type TYPE_DIR_ANC
	rpd		<- subset(rpd, !is.na(TYPE))	
	#	reset OVERLAP
	tmp		<- rpd[, which(GROUP=='TYPE_DIR_ANC')]
	set(rpd, tmp, 'OVERLAP', rpd[tmp, OVERLAP_DIR_ANC])
	rpd[, OVERLAP_DIR_ANC:=NULL]
	#
	#	for each pair count windows
	#
	rplkl	<- rpd[, list(K=length(W_FROM)), by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','COUP_SC','LABEL','LABEL_SH','GROUP','TYPE')]	
	rplkl	<- dcast.data.table(rplkl, RUN+PTY_RUN+FEMALE_SANGER_ID+MALE_SANGER_ID+COUP_SC+LABEL+LABEL_SH~GROUP+TYPE, value.var='K')
	for(x in setdiff(colnames(rplkl),c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','COUP_SC','LABEL','LABEL_SH','GROUP')))
		set(rplkl, which(is.na(rplkl[[x]])), x, 0L)	
	for(x in colnames(rplkl)[grepl('TYPE_PD_MEAN',colnames(rplkl))])
		set(rplkl, which(rplkl[[x]]>0), x, 1L)	
	rplkl	<- melt(rplkl, id.vars=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','COUP_SC','LABEL','LABEL_SH'), variable.name='GROUP', value.name='K')
	rplkl[, TYPE:= gsub('.*_([^_]+)$','\\1',GROUP)]
	set(rplkl, NULL, 'GROUP', rplkl[, gsub('(.*)_[^_]+$','\\1',GROUP)])
	#
	#	for each pair count total windows (n) and effective windows (neff)
	#
	tmp		<- rpd[, list(N=length(W_FROM), NEFF=sum(1-OVERLAP)), by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN','GROUP')]	
	rplkl	<- merge(rplkl, tmp, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN','GROUP'))
	#	define effective number of windows of type
	rplkl[, KEFF:= K/N*NEFF]	
	
	
	#
	#	prepare colours
	#	
	cols.type	<- list()
	tmp2		<- data.table(	TYPE= c("chain","no chain"),
			COLS= c(rev(brewer.pal(11, 'PuOr'))[2], rev(brewer.pal(11, 'RdBu'))[3]))					
	tmp2		<- { tmp<- tmp2[, COLS]; names(tmp) <- tmp2[, TYPE]; tmp }
	cols.type[['TYPE_DIR_NUMBER']]	<- tmp2	
	tmp2		<- data.table(	TYPE= c("chain correct direction","chain incorrect direction"),
			COLS= c(rev(brewer.pal(9, 'BuGn'))[6], rev(brewer.pal(11, 'PuOr'))[2]))					
	tmp2		<- { tmp<- tmp2[, COLS]; names(tmp) <- tmp2[, TYPE]; tmp }
	cols.type[['TYPE_DIR_ANC']]	<- tmp2
	
	#
	#	plot ancestral assignments out of all assignments
	#
	
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]
	tmp		<- subset(rplkl, GROUP%in%c('TYPE_DIR_NUMBER'))
	setkey(tmp, LABEL_SH)	
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=rev(c("chain","no chain")))])
	#tmp2	<- c(cols.type[['TYPE_PAIR_DI']],cols.type[['TYPE_PAIR_TODI3']],cols.type[['TYPE_PAIR_TODI2x2']])	
	ggplot(tmp, aes(x=LABEL, y=KEFF, fill=TYPE)) + 
			geom_bar(stat='identity',position='stack') +
			scale_y_continuous(expand=c(0,0)) +
			scale_fill_manual(values=cols.type[['TYPE_DIR_NUMBER']]) +			
			theme_bw() + 
			theme(	legend.position='bottom', 
					axis.ticks.y=element_blank()) +			
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE)) +
			labs(	x='', 
					y='non-overlapping windows\n(number)\n',
					fill='phylogenetic\nrelationship')	
	ggsave(file=file.path(dir, paste0(run,'-phsc-relationships_direction_vs_other.pdf')), w=12, h=12)
	#	average proportion of ancestral assignments out of all assignments
	tmp		<- subset(rplkl, GROUP=='TYPE_DIR_NUMBER' & TYPE=='chain')
	tmp[, list(P_ANCESTRAL_MEAN=mean(KEFF/NEFF), P_ANCESTRAL_MEDIAN=median(KEFF/NEFF))]
	#	   		P_ANCESTRAL_MEAN 	P_ANCESTRAL_MEDIAN
	#        	0.6706306           0.65625
	
	#
	#	plot ancestral assignments in correct direction
	#
	
	tmp		<- subset(rplkl, GROUP%in%c('TYPE_DIR_ANC'))
	setkey(tmp, LABEL_SH)	
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=rev(c("chain correct direction","chain incorrect direction")))])
	#tmp2	<- c(cols.type[['TYPE_PAIR_DI']],cols.type[['TYPE_PAIR_TODI3']],cols.type[['TYPE_PAIR_TODI2x2']])	
	ggplot(tmp, aes(x=LABEL, y=KEFF, fill=TYPE)) + 
			geom_bar(stat='identity',position='stack') +
			scale_y_continuous(expand=c(0,0)) +
			scale_fill_manual(values=cols.type[['TYPE_DIR_ANC']]) +			
			theme_bw() + 
			theme(	legend.position='bottom', 
					axis.ticks.y=element_blank()) +			
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE)) +
			labs(	x='', 
					y='non-overlapping windows with ancestral assignments\n(number)\n',
					fill='phylogenetic\nrelationship')	
	ggsave(file=file.path(dir, paste0(run,'-phsc-relationships_directioncorrect.pdf')), w=12, h=12)
	
	tmp		<- subset(rplkl, GROUP%in%c('TYPE_DIR_ANC') & TYPE=='chain correct direction')	
	tmp[, list(P_DIR_MEAN=mean(KEFF/NEFF), P_DIR_MEDIAN=median(KEFF/NEFF))]
	#		P_DIR_MEAN 		P_DIR_MEDIAN
	#1:  	0.7705332    	0.8113208
	
	#
	#	add prior (equal probabilities to all types), posterior, and marginal probabilities
	#
	rplkl[, DIR_PRIOR:= 0.1]
	rplkl[, DIR_PO:= DIR_PRIOR+KEFF]
	tmp		<- rplkl[, {
				alpha	<- DIR_PO
				beta	<- sum(DIR_PO)-DIR_PO				
				list(	TYPE=TYPE, 
						POSTERIOR_ALPHA=alpha, 
						POSTERIOR_BETA=beta, 
						POSTERIOR_MEAN=alpha/(alpha+beta), 
						POSTERIOR_QL=qbeta(0.025, alpha, beta), POSTERIOR_QU=qbeta(0.975, alpha, beta),
						POSTERIOR_IL=qbeta(0.25, alpha, beta), POSTERIOR_IU=qbeta(0.75, alpha, beta))	
			}, by=c('PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','RUN','GROUP')]
	rplkl	<- merge(rplkl, tmp, by=c('RUN','PTY_RUN','FEMALE_SANGER_ID','MALE_SANGER_ID','GROUP','TYPE'))
	#	fixup TYPE_PD_MEAN: there are no posteriors here
	set(rplkl, rplkl[, which(GROUP=='TYPE_PD_MEAN' & K==0)], c('POSTERIOR_MEAN'), 0)
	set(rplkl, rplkl[, which(GROUP=='TYPE_PD_MEAN' & K>0)], c('POSTERIOR_MEAN'), 1)
	set(rplkl, rplkl[, which(GROUP=='TYPE_PD_MEAN')], c('KEFF','DIR_PO','DIR_PRIOR','POSTERIOR_ALPHA','POSTERIOR_BETA','POSTERIOR_QL','POSTERIOR_QU','POSTERIOR_IL','POSTERIOR_IU'), NA_real_)
	#
	#	make table by TYPE of effective counts,  posterior probs,	alpha/beta
	#	as on blackboard
	#
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]
	for(group in c('TYPE_PAIR_TODI2x2','TYPE_PAIR_TODI3'))
	{
		rpt	<- subset(rplkl, GROUP==group)
		rpt[, POSTERIOR:= paste0( round(100*POSTERIOR_MEAN, d=1), '% (', round(100*POSTERIOR_IL, d=1), '%, ',round(100*POSTERIOR_IU, d=1),'%)')]
		suppressWarnings({rpt	<- melt(rpt, measure.vars=c('KEFF','POSTERIOR'))})
		rpt	<- dcast.data.table(rpt, RUN+PTY_RUN+FEMALE_SANGER_ID+MALE_SANGER_ID+LABEL+LABEL_SH+COUP_SC~variable+TYPE, value.var='value')
		setkey(rpt, COUP_SC, LABEL_SH)
		set(rpt, NULL, c('LABEL_SH'), NULL)
		write.csv(rpt, row.names=FALSE, file=file.path(dir, paste0(run,'_table_',group,'.csv')))		
	}	 
	
	
	#
	#	plot 'TYPE_PAIR_DI' 'TYPE_PAIR_TODI3' next to each other
	#
	g_legend<-function(a.gplot)
	{
		tmp <- ggplot_gtable(ggplot_build(a.gplot))
		leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
		legend <- tmp$grobs[[leg]]
		legend
	}
	
	tmp		<- subset(rplkl, COUP_SC=='seroinc' & GROUP%in%c('TYPE_PAIR_DI','TYPE_PAIR_TODI3'))
	setkey(tmp, LABEL_SH)
	set(tmp, NULL, 'GROUP', tmp[, factor(GROUP, levels=rev(c('TYPE_PAIR_DI','TYPE_PAIR_TODI3','TYPE_PAIR_TODI2x2')))])
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=rev(c("close","intermediate\ndistance","distant","no intermediate\n and close","no intermediate\n but not close","with intermediate\nor distant","close ancestral/\nintermingled","close other","not close ancestral/\nintermingled","not close other")))])
	tmp2	<- c(cols.type[['TYPE_PAIR_DI']],cols.type[['TYPE_PAIR_TODI3']],cols.type[['TYPE_PAIR_TODI2x2']])	
	p1		<- ggplot(tmp, aes(x=GROUP, y=KEFF, fill=TYPE)) + 
			geom_bar(stat='identity',position='stack') +
			scale_y_continuous(expand=c(0,0)) +
			scale_fill_manual(values=tmp2) +
			theme_bw() + 
			theme(	legend.position='bottom', axis.text.y=element_blank(),
					axis.ticks.y=element_blank(),
					panel.spacing=unit(0.4, "lines"), strip.text.y = element_text(angle=180),
					strip.background=element_rect(fill="transparent", colour="transparent"),
					panel.border=element_rect(color="transparent")) +
			facet_grid(LABEL_SH~., switch='y') +
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE)) +
			labs(	x='', 
					y='non-overlapping windows\n(number)\n',
					fill='phylogenetic\nrelationship')	
	p2		<- ggplot(subset(tmp, TYPE%in%c("no intermediate\n and close","close")), aes(x=GROUP, middle=qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA), lower=POSTERIOR_IL, upper=POSTERIOR_IU, ymin=POSTERIOR_QL, ymax=POSTERIOR_QU, fill=TYPE)) + 
			geom_boxplot(stat='identity') +
			scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.2), labels=percent) +
			scale_fill_manual(values=tmp2) +
			theme_bw() + 
			theme(	legend.position='bottom', axis.text.y=element_blank(),
					axis.ticks.y=element_blank(),
					panel.spacing=unit(0.4, "lines"), strip.text.y=element_blank(),
					strip.background=element_blank(),
					panel.border=element_rect(color="transparent")) +
			facet_grid(LABEL_SH~., switch='y') +
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE)) +
			labs(	x='', 
					y='posterior probability\n\n',
					fill='phylogenetic\nrelationship')
	p3		<- g_legend(p1)	
	p3$vp	<- viewport(layout.pos.row=2, layout.pos.col=1:2)	
	pdf(file=file.path(dir, paste(run,'-phsc-relationships_compare_DI_to_TODI3.pdf',sep='')), w=12, h=15)
	grid.newpage()	
	pushViewport(viewport(layout = grid.layout(2, 2, heights=unit(c(10,1), "null"), widths=unit(c(7, 3), "null"))))   	
	print(p1+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
	print(p2+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
	grid.draw(p3)
	dev.off()	
	#
	#	plot 'TYPE_PAIR_DI' 'TYPE_PAIR_TO' 'TYPE_PAIR_TODI3' next to each other
	#		
	tmp		<- subset(rplkl, COUP_SC=='seroinc' & GROUP%in%c('TYPE_PAIR_DI','TYPE_PAIR_TO','TYPE_PAIR_TODI3'))
	setkey(tmp, LABEL_SH)
	set(tmp, NULL, 'GROUP', tmp[, factor(GROUP, levels=rev(c('TYPE_PAIR_DI','TYPE_PAIR_TO','TYPE_PAIR_TODI3')))])
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=rev(c("close","intermediate\ndistance","distant","no intermediate\n and close","no intermediate\n but not close","with intermediate\nor distant","ancestral/\nintermingled","other")))])
	tmp2	<- c(cols.type[['TYPE_PAIR_DI']],cols.type[['TYPE_PAIR_TO']],cols.type[['TYPE_PAIR_TODI3']])	
	p1		<- ggplot(tmp, aes(x=GROUP, y=KEFF, fill=TYPE)) + 
			geom_bar(stat='identity',position='stack') +
			scale_y_continuous(expand=c(0,0)) +
			scale_fill_manual(values=tmp2) +
			theme_bw() + 
			theme(	legend.position='bottom', axis.text.y=element_blank(),
					axis.ticks.y=element_blank(),
					panel.spacing=unit(0.4, "lines"), strip.text.y = element_text(angle=180),
					strip.background=element_rect(fill="transparent", colour="transparent"),
					panel.border=element_rect(color="transparent")) +
			facet_grid(LABEL_SH~., switch='y') +
			coord_flip() +
			guides(fill=guide_legend(ncol=4, byrow = TRUE)) +
			labs(	x='', 
					y='non-overlapping windows\n(number)\n',
					fill='phylogenetic\nrelationship')	
	p2		<- ggplot(subset(tmp, TYPE%in%c("no intermediate\n and close","ancestral/\nintermingled","close")), aes(x=GROUP, middle=qbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA), lower=POSTERIOR_IL, upper=POSTERIOR_IU, ymin=POSTERIOR_QL, ymax=POSTERIOR_QU, fill=TYPE)) + 
			geom_boxplot(stat='identity') +
			scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.2), labels=percent) +
			scale_fill_manual(values=tmp2) +
			theme_bw() + 
			theme(	legend.position='bottom', axis.text.y=element_blank(),
					axis.ticks.y=element_blank(),
					panel.spacing=unit(0.4, "lines"), strip.text.y=element_blank(),
					strip.background=element_blank(),
					panel.border=element_rect(color="transparent")) +
			facet_grid(LABEL_SH~., switch='y') +
			coord_flip() +
			guides(fill=guide_legend(ncol=3, byrow = TRUE)) +
			labs(	x='', 
					y='posterior probability\n\n',
					fill='phylogenetic\nrelationship')
	p3		<- g_legend(p1)	
	p3$vp	<- viewport(layout.pos.row=2, layout.pos.col=1:2)	
	pdf(file=file.path(dir, paste(run,'-phsc-relationships_compare_DI_to_TO_to_TODI3.pdf',sep='')), w=12, h=20)
	grid.newpage()	
	pushViewport(viewport(layout = grid.layout(2, 2, heights=unit(c(10,1), "null"), widths=unit(c(7, 3), "null"))))   	
	print(p1+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
	print(p2+theme(legend.position = 'none'), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
	grid.draw(p3)
	dev.off()	
}

RakaiCouples.analyze.couples.170120.direction.serodisc<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	# load pairwise probabilities
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmp_bbmodels.rda')
	
	#	add MALE_RID, FEMALE_RID, COUPID
	tmp		<- subset(rp, select=c(FEMALE_SANGER_ID, MALE_SANGER_ID, MALE_RID, FEMALE_RID, COUPID))
	rplkl	<- merge(tmp, rplkl, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	#
	#	select lkl transmission pairs based on "no intermediate and close"
	#	evaluate proportion of initially serodisc couples with correct transmission estimated
	#
	confidence.cut	<- 0.5
	rpd				<- subset(rplkl, COUP_SC%in%c('M->F','F->M') & GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='no intermediate\n and close' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confidence.cut)		
	tmp				<- subset(rpd, select=c(PTY_RUN, FEMALE_SANGER_ID, MALE_SANGER_ID))
	tmp				<- merge(tmp, rplkl, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	rmf				<- subset(tmp, GROUP%in%c('TYPE_DIR_TO3') & COUP_SC=='M->F')
	rfm				<- subset(tmp, GROUP%in%c('TYPE_DIR_TO3') & COUP_SC=='F->M')
	#
	unique(subset(rplkl, COUP_SC%in%c('M->F','F->M') & GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='no intermediate\n and close'), by='COUPID')
	#	initially +/- -/+ couples: 20
	cat('\ncouples with phyloscanner evidence M->F, n=',nrow(unique(rmf, by='COUPID')))
	cat('\ncouples with phyloscanner evidence F->M, n=',nrow(unique(rfm, by='COUPID')))
	#	initially +/- couples with phyloscanner evidence for lkl trm pair, n= 9	
	#	initially -/+ couples with phyloscanner evidence for lkl trm pair, n= 3	
	rtr		<- rbind(rmf, rfm)
	rtr[, DIR_SELECT:= pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confidence.cut]
	unique(subset(rtr, (COUP_SC=='M->F' & TYPE=='mf' | COUP_SC=='F->M' & TYPE=='fm')), by='COUPID')[, mean(DIR_SELECT)]
	#	0.75
	unique(subset(rtr, (COUP_SC=='M->F' & TYPE=='fm' | COUP_SC=='F->M' & TYPE=='mf')), by='COUPID')[, mean(DIR_SELECT)]
	#	0.08333333 (1/12)
		
}

RakaiCouples.analyze.couples.170120.direction.all.analyze<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(grid)
	require(gridExtra)
	require(RColorBrewer)
	require(Hmisc)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5") )
	setnames(rpw, c('TYPE','TYPE_PAIR'), c('TYPE_DIR_TODI7x3','TYPE_PAIR_TODI3x3'))
	# load pairwise probabilities
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmp_bbmodels.rda')
	
	# 	add agrarian etc
	# 	from Kate:
	#	Trading communities: 1   4   16  22  24  33  51  107 776
	#	Fishing communities: 38, 770, 771, 774	
	# 	tmp		<- rp[, sort(unique(c(MALE_COMM_NUM,FEMALE_COMM_NUM)))]; cat(paste('"', tmp, '",', collapse='', sep=''))
	tmp		<- data.table(	COMM_NUM=	c("1","2","4","7","8","16","22","23","24","33","34","38","40","51","56","57","58","89","94","106","107","108","370","391","770","771","772","773","774","776"),
			COMM_TYPE=	c("T","A","T","A","A", "T", "T", "A", "T", "T", "A", "F", "A", "T", "A", "A", "A", "A", "A",  "A",  "T",  "A",  "A",  "A", "F",  "F",  "A",  "A",  "F",   "T"))
	set(tmp, NULL, 'COMM_TYPE', tmp[, as.character(factor(COMM_TYPE, levels=c('A','T','F'), labels=c('agrarian','trading','fisherfolk')))])
	set(tmp, NULL, 'COMM_NUM', tmp[, as.integer(COMM_NUM)])
	setnames(tmp, c('COMM_NUM','COMM_TYPE'), c('MALE_COMM_NUM','MALE_COMM_TYPE'))
	rp		<- merge(rp, tmp, by='MALE_COMM_NUM')
	setnames(tmp, c('MALE_COMM_NUM','MALE_COMM_TYPE'), c('FEMALE_COMM_NUM','FEMALE_COMM_TYPE'))
	rp		<- merge(rp, tmp, by='FEMALE_COMM_NUM')
	
	#	add RCCS IDs
	tmp		<- subset(rp, select=c(FEMALE_SANGER_ID, MALE_SANGER_ID, MALE_RID, FEMALE_RID, COUPID))
	rplkl	<- merge(tmp, rplkl, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
	#	select first couples for whom transmission cannot be excluded 
	#	no intermediate\n and close > 50%
	confidence.cut	<- 0.5
	rpd		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='no intermediate\n and close' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confidence.cut)		
	tmp		<- subset(rpd, select=c(PTY_RUN, FEMALE_SANGER_ID, MALE_SANGER_ID))
	tmp		<- merge(tmp, rplkl, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	rmf		<- subset(tmp, GROUP%in%c('TYPE_DIR_TO') & TYPE=='mf' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confidence.cut)
	rfm		<- subset(tmp, GROUP%in%c('TYPE_DIR_TO') & TYPE=='fm' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confidence.cut)
	rex		<- subset(rplkl, GROUP%in%c('TYPE_PAIR_TODI3') & TYPE=='with intermediate\nor distant' & pbeta(0.5, POSTERIOR_ALPHA, POSTERIOR_BETA, lower.tail=FALSE)>confidence.cut)
	
	cat('\ncouples with phyloscanner assessment, n=',nrow(unique(rplkl, by='COUPID')))
	cat('\ncouples for whom transmission cannot be excluded, n=',nrow(unique(rpd, by='COUPID')))
	cat('\ncouples with evidence M->F, n=',nrow(unique(rmf, by='COUPID')))
	cat('\ncouples with evidence F->M, n=',nrow(unique(rfm, by='COUPID')))
	cat('\ncouples for whom transmission can be excluded, n=',nrow(unique(rex, by='COUPID')))
	#	couples with phyloscanner assessment					215	
	#	for whom transmission can be excluded 					83		(38.6%)
	#	ambiguous												8		(3.7%)
	#	no intermediate and close								124		(57.7%)	
	#	of those we cannot determine likely direction for		37		(37/124=30%)
	#	with evidence for M->F 									56		(56/87 =64%)
	#	with evidence for F->M									31 		(31/87 =36%)
	
	rmf		<- merge(unique(subset(rmf, select=COUPID)), unique(rp, by='COUPID'),by='COUPID')
	rmf[, PHSC_DIR:='m->f']
	rfm		<- merge(unique(subset(rfm, select=COUPID)), unique(rp, by='COUPID'),by='COUPID')
	rfm[, PHSC_DIR:='f->m']
	rtr		<- rbind(rmf, rfm)
	
	rtr[, AGEDIFF:= MALE_BIRTHDATE-FEMALE_BIRTHDATE]
	rtr[, AVGAGE:= (MALE_BIRTHDATE+FEMALE_BIRTHDATE)/2]
	#
	#	is there a difference in male->female transmission by community type?
	#	results: no	
	rtr[, {
				z	<- binconf( length(which(PHSC_DIR=='m->f')), length(PHSC_DIR) )				
				list(K=length(which(PHSC_DIR=='m->f')), N=length(PHSC_DIR), P=z[1], QL=z[2], QU=z[3])
			}, by='MALE_COMM_TYPE']	
	#	MALE_COMM_TYPE  K  N         P        QL        QU
	#1:       agrarian 21 32 0.6562500 0.4831109 0.7958956
	#2:     fisherfolk 33 53 0.6226415 0.4880690 0.7406373
	#3:        trading  2  2 1.0000000 0.3423802 1.0000000
	
	#
	#	is there a difference in age gap between male->female transmission / female->male transmission ?
	#	results: no
	ggplot(rtr, aes(x=PHSC_DIR, y=AGEDIFF)) + geom_boxplot()
	tmp		<- subset(rtr, !is.na(MALE_BIRTHDATE) & !is.na(FEMALE_BIRTHDATE), select=c(PHSC_DIR,AGEDIFF))
	set(tmp, NULL, 'PHSC_DIR', tmp[, as.integer(as.character(factor(PHSC_DIR, levels=c('f->m','m->f'), labels=c('0','1'))))])
	summary(gamlss(data=tmp, PHSC_DIR~AGEDIFF, family=LO))
	summary(gamlss(data=tmp, AGEDIFF~PHSC_DIR))
	#				Estimate Std. Error t value Pr(>|t|)    
	#(Intercept)  0.626998   0.071579   8.759 1.98e-13 ***
	#AGEDIFF     -0.002135   0.008422  -0.254      0.8    
	
	#
	#	are couples overall younger in male->female transmission ?
	#	results: no
	ggplot(rtr, aes(x=PHSC_DIR, y=AVGAGE)) + geom_boxplot()
	tmp		<- subset(rtr, !is.na(AVGAGE), select=c(PHSC_DIR,AVGAGE))
	summary(rq(AVGAGE~PHSC_DIR, tau=.5, data=tmp, method='fn'), se='nid')
	
	
	#
	#	are female transmitters younger than female recipients ?
	#	results: yes but not significant if all couples are considered
	#			 yes and signficant if fairly old couples are excluded (female age < 40)
	ggplot(rtr, aes(x=PHSC_DIR, y=FEMALE_BIRTHDATE)) + geom_boxplot()	
	ggplot(subset(rtr,FEMALE_BIRTHDATE>1975), aes(x=PHSC_DIR, y=FEMALE_BIRTHDATE)) + geom_boxplot()	
	#	the female birthdate distribution among m->f is asymmetric 
	#	the female birthdate distribution has outliers in f->m
	#	try median
	tmp		<- subset(rtr, !is.na(FEMALE_BIRTHDATE), select=c(PHSC_DIR,FEMALE_BIRTHDATE))
	summary(rq(FEMALE_BIRTHDATE~PHSC_DIR, tau=.5, data=tmp, method='fn'), se='nid')
	#	             Value      Std. Error t value    Pr(>|t|)  
	#(Intercept)  1986.31800    1.16809 1700.48527    0.00000
	#PHSC_DIRm->f   -1.80805    1.54660   -1.16904    0.24565

	tmp		<- subset(rtr, !is.na(FEMALE_BIRTHDATE) & FEMALE_BIRTHDATE>1975, select=c(PHSC_DIR,FEMALE_BIRTHDATE))
	summary(rq(FEMALE_BIRTHDATE~PHSC_DIR, tau=.5, data=tmp, method='fn'), se='nid')
	#             	Value      Std. Error t value    Pr(>|t|)  
	#(Intercept)  1987.95300    0.78721 2525.31962    0.00000
	#PHSC_DIRm->f   -3.27014    1.33152   -2.45595    0.01625
	

	#
	#	are male transmitters younger than male recipients ?
	#	results: no, male transmitter tend to be older than male recipients 
	#			 but this is not significant
	ggplot(rtr, aes(x=PHSC_DIR, y=MALE_BIRTHDATE)) + geom_boxplot()	
	tmp		<- subset(rtr, !is.na(MALE_BIRTHDATE), select=c(PHSC_DIR,MALE_BIRTHDATE))
	summary(rq(MALE_BIRTHDATE~PHSC_DIR, tau=.5, data=tmp, method='fn'), se='nid')
	#				Value      Std. Error t value    Pr(>|t|)  
	#(Intercept)  1981.17500    2.16624  914.56890    0.00000
	#PHSC_DIRm->f   -2.29300    2.55394   -0.89783    0.37184
	tmp		<- subset(rtr, !is.na(MALE_BIRTHDATE) & FEMALE_BIRTHDATE>1975 & MALE_BIRTHDATE>1975, select=c(PHSC_DIR,MALE_BIRTHDATE))
	summary(rq(MALE_BIRTHDATE~PHSC_DIR, tau=.5, data=tmp, method='fn'), se='nid')
	#	
}

RakaiCouples.analyze.couples.161219.distances<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	require(ape)
	require(phyloscan)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments.rda')
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl1_d5") )
	#	
	rpw[, table(RUN, useNA='if')]
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]
	
	#	define order: largest number of trm assignments
	tmp		<- rpw[, list( WIN_TR=length(which(grepl('close|anc',TYPE))) ), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	tmp		<- tmp[order(-WIN_TR),]
	tmp[, PLOT_ID:=seq_len(nrow(tmp))]	
	#	define short label
	tmp		<- merge(tmp, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
	tmp[, LABEL_SH:= factor(PLOT_ID, levels=PLOT_ID, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair ', COUPID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
								'\n<->', 
								'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
								'\n',sep=''))]
	tmp		<- subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL, LABEL_SH))	
	rpw		<- merge(rpw, tmp, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))		
	
	
	#
	#	patristic distances across the genome and smooth
	#	
	ggplot(rpw, aes(x=W_FROM, y=PATRISTIC_DISTANCE)) + geom_point(alpha=0.2) +
			scale_y_log10(labels=percent, limits=c(0.001, 0.7), expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			scale_x_continuous(breaks=seq(1,1e4,500)) +
			geom_smooth(method='loess', span=.5) +
			theme_bw() + theme(legend.position='top') +
			labs(x='phyloscanner window', y='patristic distance in read tree\n(subst/site)')			
	ggsave(file=paste0(dir,'/',run,'_distances_across_genome.pdf'), w=15, h=6, limitsize=FALSE)
	
	#	select windows from gag+pol
	rpw		<- subset(rpw, W_FROM>=800 & W_TO<=4650)
	
	#
	#	how good is the mean?
	#
	
	#	plot couples boxplots with repetitions in colours by couple
	tmp		<- subset(rpw, COUP_SC=='seroinc')
	tmp		<- melt(tmp, id.vars=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN','LABEL_SH','W_FROM'), measure.vars=c('PATRISTIC_DISTANCE','CONS_GD','CONS_PD_FT','CONS_PD_EX'))
	ggplot(tmp, aes(x=W_FROM, y=value, colour=variable)) + geom_point() +
			theme_bw() + theme(legend.position='top') +
			labs(x='phyloscanner window', y='subst/site', colour='distance measure') +
			facet_wrap(~LABEL_SH, scales='free_y', ncol=1)
	ggsave(file=paste0(dir,'/',run,'_distances_seroinc.pdf'), w=8, h=50, limitsize=FALSE)
	
	#
	#	calculate mean and 50% CIs 95% CIs
	#
	
	rpd		<- rpw[, list(	PHSC_PD_MEAN=mean(PATRISTIC_DISTANCE),
				PHSC_PD_Q025=quantile(PATRISTIC_DISTANCE, p=0.025),
				PHSC_PD_Q25=quantile(PATRISTIC_DISTANCE, p=0.25),
				PHSC_PD_Q75=quantile(PATRISTIC_DISTANCE, p=0.75),
				PHSC_PD_Q975=quantile(PATRISTIC_DISTANCE, p=0.975),	
				CONS_GD=CONS_GD[1], CONS_PD_FT=CONS_PD_FT[1], CONS_PD_EX=CONS_PD_EX[1]
				), by=c('RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN','LABEL_SH','COUP_SC')]

	#	
	#	correlation plots
	#

	#	correlation plot ExaML
	ggplot(rpd, aes(x=CONS_PD_EX, y=PHSC_PD_MEAN, ymin=PHSC_PD_Q25, ymax=PHSC_PD_Q75)) + 
			geom_abline(slope=1, intercept=0, colour='black', linetype='dotted') +
			geom_linerange(size=.5, alpha=0.5, colour='grey40') + geom_point(size=1) + 
			scale_x_log10(labels=percent, expand=c(0,0), limits=c(0.002, 0.5), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			scale_y_log10(labels=percent, expand=c(0,0), limits=c(0.002, 0.5), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			theme_bw() +
			labs(	x='\npatristic distance in consensus tree\n(subst/site)',
					y='patristic distance in read tree\nmean and interquartile range\n(subst/site)\n')
	ggsave(file=paste0(dir,'/',run,'_distances_consExaML.pdf'), w=7, h=7)
	#	correlation plot FastTree
	ggplot(rpd, aes(x=CONS_PD_FT, y=PHSC_PD_MEAN, ymin=PHSC_PD_Q25, ymax=PHSC_PD_Q75)) + 
			geom_abline(slope=1, intercept=0, colour='black', linetype='dotted') +
			geom_linerange(size=.5, alpha=0.5, colour='grey40') + geom_point(size=1) + 
			scale_x_log10(labels=percent, expand=c(0,0), limits=c(0.002, 0.5), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			scale_y_log10(labels=percent, expand=c(0,0), limits=c(0.002, 0.5), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			theme_bw() +
			labs(	x='\npatristic distance in consensus tree\n(subst/site)',
					y='patristic distance in read tree\nmean and interquartile range\n(subst/site)\n')
	ggsave(file=paste0(dir,'/',run,'_distances_consFastTree.pdf'), w=7, h=7)	
	#	correlation plot genetic distance
	ggplot(rpd, aes(x=CONS_GD, y=PHSC_PD_MEAN, ymin=PHSC_PD_Q25, ymax=PHSC_PD_Q75)) + 
			geom_abline(slope=1, intercept=0, colour='black', linetype='dotted') +
			geom_linerange(size=.5, alpha=0.5, colour='grey40') + geom_point(size=1) + 
			scale_x_log10(labels=percent, expand=c(0,0), limits=c(0.002, 0.5), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			scale_y_log10(labels=percent, expand=c(0,0), limits=c(0.002, 0.5), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			theme_bw() +
			labs(	x='\ngenetic distance between consensus sequences\n(subst/site)',
					y='patristic distance in read tree\nmean and interquartile range\n(subst/site)\n')
	ggsave(file=paste0(dir,'/',run,'_distances_consGenetic.pdf'), w=7, h=7)
	
	#
	#	correlations	
	#	
	
	rpd[, cor( PHSC_PD_MEAN, CONS_PD_EX)]				# 0.910854
	rpd[, cor( log10(PHSC_PD_MEAN), log10(CONS_PD_EX))]	# 0.9211325
	rpd[, cor( PHSC_PD_MEAN, CONS_PD_FT)]				# 0.9191956
	rpd[, cor( log10(PHSC_PD_MEAN), log10(CONS_PD_FT))]	# 0.927517
	
	#
	#	bimodality
	#

	tmp		<- melt(rpd, id.vars=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN','LABEL_SH'), measure.vars=c('PHSC_PD_MEAN','CONS_PD_EX'))
	set(tmp, NULL, 'variable', tmp[, as.character(factor(as.character(variable), levels=c('PHSC_PD_MEAN','CONS_PD_EX'), labels=c('average patristic distance in read trees','patristic distance in consensus tree')))])
	#tmp		<- subset(tmp, variable=='PHSC_PD_MEAN')
	ggplot(tmp, aes(x=log10(value), colour=variable)) +
			geom_density(adjust=0.5) +
			scale_colour_brewer(palette='Set1') +
			scale_y_continuous(expand=c(0,0), limits=c(0,1.1)) + 
 			scale_x_continuous(		limits=log10(c(0.002, 0.5)), expand=c(0,0), 
									breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)), 
									labels=paste0(100*c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),'%')) +
			theme_bw() + theme(legend.position='bottom') + 
			labs(x='\npatristic distance\n(subst/site)', y='', colour='')
	ggsave(file=paste0(dir,'/',run,'_distances_bimodal.pdf'), w=6, h=6)	
	
	tmp		<- melt(rpd, id.vars=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PTY_RUN','LABEL_SH'), measure.vars=c('PHSC_PD_MEAN','CONS_GD'))
	set(tmp, NULL, 'variable', tmp[, as.character(factor(as.character(variable), levels=c('PHSC_PD_MEAN','CONS_GD'), labels=c('average patristic distance in read trees','genetic distance among consensus sequences')))])
	#tmp		<- subset(tmp, variable=='PHSC_PD_MEAN')
	ggplot(tmp, aes(x=log10(value), colour=variable)) +
			geom_density(adjust=0.5) +
			scale_colour_brewer(palette='Set1') +
			scale_y_continuous(expand=c(0,0), limits=c(0,1.7)) + 
			scale_x_continuous(		limits=log10(c(0.002, 0.5)), expand=c(0,0), 
					breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)), 
					labels=paste0(100*c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),'%')) +
			geom_vline(xintercept=log10(c(0.035,0.07))) +
			theme_bw() + theme(legend.position='bottom') + 
			labs(x='\npatristic distance\n(subst/site)', y='', colour='')
	ggsave(file=paste0(dir,'/',run,'_distances_bimodal_gd.pdf'), w=6, h=6)	
}
	
RakaiCouples.analyze.couples.161219.distances.window.topology<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl1_d5", "RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5", "RCCS_161219_w270_d50_p001_mr20_mt1_cl1_dInf", "RCCS_161219_w270_d50_p001_mr20_mt1_cl2_dInf", "RCCS_161219_w270_d50_p001_mr20_mt1_clInf_dInf") )
	#
	set(rpw, NULL, 'RUN', rpw[, factor(RUN, levels=c( 	"RCCS_161219_w270_d50_p001_mr20_mt1_cl1_d5", "RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5", "RCCS_161219_w270_d50_p001_mr20_mt1_cl1_dInf", "RCCS_161219_w270_d50_p001_mr20_mt1_cl2_dInf", "RCCS_161219_w270_d50_p001_mr20_mt1_clInf_dInf"))])
	rpw[, table(RUN, useNA='if')]
	#
	#	for each run: plot pairs	
	#	
	#	define plotting order: largest number of trm assignments
	tmp		<- rpw[, list( WIN_TR=length(which(grepl('close|anc',TYPE))) ), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	tmp		<- tmp[order(-WIN_TR),]
	tmp[, PLOT_ID:=seq_len(nrow(tmp))]	
	#	define label
	tmp		<- merge(tmp, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair ', COUPID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL))
	rpw		<- merge(tmp, rpw, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	#	define min/max reads
	tmp		<- rpw[, list(	ID_R_MIN=min(ID1_R, ID2_R),
					ID_R_MAX=max(ID1_R, ID2_R)), by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM')]
	rpw	<- merge(rpw, tmp, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM'))
	#	define topo types
	rpw[, TYPE_PAIR_TOPO:= gsub('_.*','',TYPE_PAIR)]
	#	make pretty
	set(rpw, NULL, 'RUN', rpw[, gsub('_','\n',gsub('w270_','',gsub('RCCS_','',as.character(RUN))))])
	set(rpw, NULL, 'RUN', rpw[, factor(RUN, levels=c(	"161219\nd50\np001\nmr20\nmt1\ncl1\nd5", "161219\nd50\np001\nmr20\nmt1\ncl2\nd5",    
									"161219\nd50\np001\nmr20\nmt1\ncl1\ndInf", "161219\nd50\np001\nmr20\nmt1\ncl2\ndInf",
									"161219\nd50\np001\nmr20\nmt1\nclInf\ndInf"))])
	
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]
	#
	#	select SEROINC
	#
	rpw	<- subset(rpw, RUN=="161219\nd50\np001\nmr20\nmt1\ncl2\nd5")
	
	#
	#	define colours for TYPE_PAIR_TOPO	
	#

	cols.typet	<- data.table(	TYPE= c("pair", "withintermediate", "other"),
								COLS= c(brewer.pal(11, 'PuOr')[2], rev(brewer.pal(11, 'RdBu'))[4], rev(brewer.pal(11, 'RdGy'))[4]))					
	cols.typet	<- { tmp<- cols.typet[, COLS]; names(tmp) <- cols.typet[, TYPE]; tmp }
	
	#
	#	plot distances by topology assignment
	#

	set(rpw, NULL, 'TYPE_PAIR_TOPO', rpw[, factor(TYPE_PAIR_TOPO, levels=c('pair','withintermediate','other'))])
	ggplot(rpw, aes(x=TYPE_PAIR_TOPO, fill=TYPE_PAIR_TOPO, y=log10(PATRISTIC_DISTANCE))) +
			geom_violin(adjust=1.5, scale='area') +
			scale_y_continuous(		limits=log10(c(0.0009, 0.8)), expand=c(0,0), 
									breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)), 
									labels=paste0(100*c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),'%')) +
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			facet_grid(~COUP_SC) +
			labs(x='', y='patristic distance in read tree\n(subst/site)\n', fill='topology assignment')
	ggsave(file=file.path(dir, paste(run,'-topodist_pairtypes_wide.pdf',sep='')), w=15, h=6)
	ggplot(rpw, aes(x=TYPE_PAIR_TOPO, fill=TYPE_PAIR_TOPO, y=log10(PATRISTIC_DISTANCE))) +
			geom_violin(adjust=2, scale='area') +
			scale_y_continuous(		limits=log10(c(0.0009, 0.8)), expand=c(0,0), 
					breaks=log10(c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)), 
					labels=paste0(100*c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),'%')) +
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			labs(x='', y='patristic distance in read tree\n(subst/site)\n', fill='topology assignment')
	ggsave(file=file.path(dir, paste(run,'-topodist_pairtypes.pdf',sep='')), w=6, h=6)
	
	#
	#	plot topology assignment by distance
	#
	rpw[, PATRISTIC_DISTANCE_LOG:= log10(PATRISTIC_DISTANCE)]
	rpw[, PATRISTIC_DISTANCE_LOGC:= cut(PATRISTIC_DISTANCE_LOG, breaks=log10(c(1e-12, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 2)),
																labels=c('<0.2%', '0.2%-0.4%', '0.4%-0.6%', '0.6%-0.8%', '0.8%-1%', '1%-2%', '2%-4%', '4%-6%', '6%-8%', '8%-10%', '10%-20%', '>20%'))]
	ggplot(rpw, aes(x=PATRISTIC_DISTANCE_LOGC, fill=TYPE_PAIR_TOPO)) +
			geom_bar(position='fill') +
			scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), expand=c(0,0)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			labs(x='\npatristic distance in read tree\n(subst/site)', y='', fill='topology assignment')
	ggsave(file=file.path(dir, paste(run,'-topodist_pairtypes2.pdf',sep='')), w=10, h=5)
	ggplot(rpw, aes(x=PATRISTIC_DISTANCE_LOGC, fill=TYPE_PAIR_TOPO)) +
			geom_bar(position='fill') +
			scale_y_continuous(labels=percent, breaks=seq(0,1,0.2), expand=c(0,0)) +
			theme_bw() + theme(legend.position='bottom') +
			scale_fill_manual(values=cols.typet) +
			facet_grid(COUP_SC~.) +
			labs(x='\npatristic distance in read tree\n(subst/site)', y='', fill='topology assignment')
	ggsave(file=file.path(dir, paste(run,'-topodist_pairtypes2_long.pdf',sep='')), w=10, h=10)
}
	
RakaiCouples.analyze.couples.161219.trmw<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trmw_assignments.rda')	
	rpw		<- subset(rpw, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl1_d5", "RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5", "RCCS_161219_w270_d50_p001_mr20_mt1_cl1_dInf", "RCCS_161219_w270_d50_p001_mr20_mt1_cl2_dInf", "RCCS_161219_w270_d50_p001_mr20_mt1_clInf_dInf") )
	#
	set(rpw, NULL, 'RUN', rpw[, factor(RUN, levels=c( 	"RCCS_161219_w270_d50_p001_mr20_mt1_cl1_d5", "RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5", "RCCS_161219_w270_d50_p001_mr20_mt1_cl1_dInf", "RCCS_161219_w270_d50_p001_mr20_mt1_cl2_dInf", "RCCS_161219_w270_d50_p001_mr20_mt1_clInf_dInf"))])
	rpw[, table(RUN, useNA='if')]
	#
	#	for each run: plot pairs	
	#	
	#	define plotting order: largest number of trm assignments
	tmp		<- rpw[, list( WIN_TR=length(which(grepl('close|anc',TYPE))) ), by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	tmp		<- tmp[order(-WIN_TR),]
	tmp[, PLOT_ID:=seq_len(nrow(tmp))]	
	#	define label
	tmp		<- merge(tmp, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair ', COUPID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL))
	rpw		<- merge(tmp, rpw, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))		
	tmp		<- rpw[, list(	ID_R_MIN=min(ID1_R, ID2_R),
							ID_R_MAX=max(ID1_R, ID2_R)), by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM')]
	rpw	<- merge(rpw, tmp, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID','W_FROM'))
	set(rpw, NULL, 'RUN', rpw[, gsub('_','\n',gsub('w270_','',gsub('RCCS_','',as.character(RUN))))])
	set(rpw, NULL, 'RUN', rpw[, factor(RUN, levels=c(	"161219\nd50\np001\nmr20\nmt1\ncl1\nd5", "161219\nd50\np001\nmr20\nmt1\ncl2\nd5",    
														"161219\nd50\np001\nmr20\nmt1\ncl1\ndInf", "161219\nd50\np001\nmr20\nmt1\ncl2\ndInf",
														"161219\nd50\np001\nmr20\nmt1\nclInf\ndInf"))])
	
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rpw$DIR[1]
	#
	#	select SEROINC
	#
	rpw	<- subset(rpw, COUP_SC=='seroinc')
	#
	#	define colours for TYPE and TYPE_PAIR	
	#		
	cols.type	<- do.call('rbind',list(
					data.table(	TYPE= c("chain_fm_nointermediate_close","chain_fm_nointermediate","chain_fm_nointermediate_distant"),
							COLS= brewer.pal(11, 'PiYG')[c(1,2,4)]),
					data.table(	TYPE= c("chain_mf_nointermediate_close","chain_mf_nointermediate","chain_mf_nointermediate_distant"),
							COLS= brewer.pal(11, 'PuOr')[c(1,2,4)]),
					data.table(	TYPE= c("intermingled_nointermediate_close","intermingled_nointermediate","intermingled_nointermediate_distant"),
							COLS= brewer.pal(11, 'PRGn')[c(1,2,4)]),
					data.table(	TYPE= c("chain_fm_withintermediate_close","chain_fm_withintermediate","chain_fm_withintermediate_distant"),
							COLS= rev(brewer.pal(11, 'BrBG'))[c(3,4,5)]),
					data.table(	TYPE= c("chain_mf_withintermediate_close","chain_mf_withintermediate","chain_mf_withintermediate_distant"),
							COLS= rev(brewer.pal(11, 'PRGn'))[c(3,4,5)]),
					data.table(	TYPE= c("intermingled_withintermediate_close","intermingled_withintermediate","intermingled_withintermediate_distant"),
							COLS= rev(brewer.pal(11, 'RdBu'))[c(3,4,5)]),
					data.table(	TYPE= c("other_close","other","other_distant"),
							COLS= rev(brewer.pal(11, 'RdGy'))[c(3,4,5)])))
	cols.type	<- { tmp<- cols.type[, COLS]; names(tmp) <- cols.type[, TYPE]; tmp }
	cols.typep	<- do.call('rbind',list(
					data.table(	TYPE= c("pair_close","pair","pair_distant"),
							COLS= brewer.pal(11, 'PuOr')[c(1,2,4)]),
					data.table(	TYPE= c("withintermediate_close","withintermediate","withintermediate_distant"),
							COLS= rev(brewer.pal(11, 'RdBu'))[c(3,4,5)]),
					data.table(	TYPE= c("other_close","other","other_distant"),
							COLS= rev(brewer.pal(11, 'RdGy'))[c(3,4,5)])))
	cols.typep	<- { tmp<- cols.typep[, COLS]; names(tmp) <- cols.typep[, TYPE]; tmp }
	#
	#	plot
	#
	p		<- lapply(rpw[, unique(LABEL)], function(label)
			{
				p	<- ggplot(subset(rpw, LABEL==label), aes(x=W_FROM)) +			
						geom_bar(aes(y=ID_R_MAX, colour=TYPE), stat='identity', fill='transparent') +
						geom_bar(aes(y=ID_R_MIN, fill=TYPE), stat='identity', colour='transparent') +
						labs(x='window start', y='number of reads', fill='topology of clades\nbetween patient pairs', colour='topology of clades\nbetween patient pairs') +
						scale_fill_manual(values=cols.type) +
						scale_colour_manual(values=cols.type) +
						scale_x_continuous(breaks=seq(0,1e4,250), limits=c(rpw[, min(W_FROM)], tmp[, max(W_FROM)])) +
						scale_y_log10(breaks=c(10,100,1000,1e4,1e5)) +
						theme_bw() + theme(legend.position='top') +
						facet_grid(RUN~LABEL) +
						guides(fill=guide_legend(ncol=2)) 
				p
			})
	pdf(file=file.path(dir, paste(run,'-phsc-seroincpairs-windowassignments_alltypes.pdf',sep='')), w=25, h=2*tmp[, length(unique(RUN))])
	for(i in seq_along(p))	print(p[[i]])
	dev.off()
	p		<- lapply(rpw[, unique(LABEL)], function(label)
			{
				p	<- ggplot(subset(rpw, LABEL==label), aes(x=W_FROM)) +			
						geom_bar(aes(y=ID_R_MAX, colour=TYPE_PAIR), stat='identity', fill='transparent') +
						geom_bar(aes(y=ID_R_MIN, fill=TYPE_PAIR), stat='identity', colour='transparent') +
						labs(x='window start', y='number of reads', fill='topology of clades\nbetween patient pairs', colour='topology of clades\nbetween patient pairs') +
						scale_fill_manual(values=cols.typep) +
						scale_colour_manual(values=cols.typep) +
						scale_x_continuous(breaks=seq(0,1e4,250), limits=c(rpw[, min(W_FROM)], tmp[, max(W_FROM)])) +
						scale_y_log10(breaks=c(10,100,1000,1e4,1e5)) +
						theme_bw() + theme(legend.position='top') +
						facet_grid(RUN~LABEL) +
						guides(fill=guide_legend(ncol=2)) 
				p
			})
	pdf(file=file.path(dir, paste(run,'-phsc-seroincpairs-windowassignments_pairtypes.pdf',sep='')), w=25, h=2*tmp[, length(unique(RUN))])
	for(i in seq_along(p))	print(p[[i]])
	dev.off()
	#
	#	get distribution of cherries on manually identified unlinked and linked
	#
	tmp	<- as.data.table(matrix( c('83','15715_1_17','15715_1_18',
							'36','15715_1_17','15715_1_18',
							'85','15699_1_38','15699_1_40',
							'104','15699_1_3','15699_1_4',
							'106','15736_1_28','15102_1_16',
							'121','16059_1_66','16059_1_52',
							'31','15981_1_31','15981_1_33',
							'98','15916_1_7','15916_1_21',
							'64','15916_1_7','15916_1_21',
							'122','15910_1_76','15910_1_82',
							'24','15981_1_17','15835_1_32',
							'70','15833_1_84','15833_1_83',
							'104','16059_1_66','16059_1_52',
							'24','16016_1_11','16017_1_44',
							'24','15915_1_51','15776_1_32',
							'15','16033_1_76','16033_1_60',
							'7','16033_1_76','16033_1_60',
							'18','16019_1_45','16018_1_44',
							'67','15892_1_3','15892_1_1',
							'91','16021_1_66','16021_1_58',
							'29','16021_1_66','16021_1_58',
							'116','15915_1_84','15915_1_83',
							'64','15915_1_84','15915_1_83',
							'50','15115_1_4','15115_1_9',
							'86','15978_1_12','15777_1_42',
							'68','15978_1_12','15777_1_42',
							'44','16056_1_83','15099_1_59'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp[, CLASS_OR:= 'OR_disconnected']	
	tmp2	<- as.data.table(matrix( c('111','15950_1_27','15950_1_25',
							'122','16056_1_58','16056_1_57',
							'38','15105_1_35','16016_1_4',
							'113','15950_1_5','15958_1_47',				
							'99','15915_1_61','15915_1_74',
							'117','15915_1_61','15915_1_74',
							'60','16056_1_83','16056_1_85',
							'39','16034_1_86','16034_1_82',
							'120','15981_1_19','15981_1_23',				
							'79','15978_1_36','15978_1_41'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp2[, CLASS_OR:= 'OR_ambiguous']
	tmp		<- rbind(tmp, tmp2)
	tmp2	<- as.data.table(matrix( c('101','16003_1_76','16003_1_79',
							'37','15114_1_64','15758_1_73',
							'37','16016_1_5','16016_1_4',
							'20','16057_1_15','16057_1_16',
							'59','15699_1_5','15867_1_21',
							'99','15915_1_78','15915_1_64',
							'109','15777_1_29','15777_1_21',
							'63','16003_1_86','16003_1_85',
							'109','15776_1_19','15776_1_28',
							'27','16219_1_78','16219_1_80',
							'120','16059_1_73','16059_1_53',
							'25','15675_1_87','15675_1_77',
							'70','15835_1_21','15835_1_22'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp2[, CLASS_OR:= 'OR_transmission']
	tmp		<- rbind(tmp, tmp2)
	set(tmp, NULL, 'PTY_RUN', tmp[, as.integer(PTY_RUN)])
	rpi		<- copy(tmp)
	ubset(rps, COUP_SC=='F->M')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int']), CHER=WIN_OF_TYPE_N[TYPE=='cher'], NO_ANC=sum(WIN_OF_TYPE_N[TYPE=='disconnected'|TYPE=='unint'])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(rps, COUP_SC=='M->F')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='trans_mf'|TYPE=='trans_mf'|TYPE=='int']), CHER=WIN_OF_TYPE_N[TYPE=='cher'], NO_ANC=sum(WIN_OF_TYPE_N[TYPE=='disconnected'|TYPE=='unint'])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= ANY_ANC]
	rpa		<- melt(rpa, measure.vars=c('ANY_ANC', 'NO_ANC','CHER'))
	set(rpa, rpa[, which(variable=='ANY_ANC')], 'variable', 'transmission assignments or intermingled')
	set(rpa, rpa[, which(variable=='NO_ANC')], 'variable', 'unint or disconnected')
	set(rpa, rpa[, which(variable=='CHER')], 'variable', 'cherry')
	setkey(rpa, PAIR_ID)	
	tmp		<- unique(subset(rpa, RUN=='RCCS_161107_w270_d200_r004_mr1'))[order(-WIN_TRM),][, list(COUPID=COUPID, MALE_SANGER_ID=MALE_SANGER_ID, FEMALE_SANGER_ID=FEMALE_SANGER_ID, PTY_RUN=PTY_RUN, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PTY_RUN','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set2') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~RUN) +
			coord_flip() +
			labs(	x= '\nsequence pairs of couples', 
					y='number of windows\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trmpair_windows.pdf',sep='')), w=35, h=10)
	#
	#	define three transmission evidence groups manually
	#	
	rpi	<- subset(rpa, RUN=='RCCS_161107_w270_d200_r004_mr1', c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID))
	setkey(rpi, CLASS_RANK)
	tmp	<- as.data.table(matrix( c('5','15861_1_13','15861_1_18',
							'118','15097_1_88','15714_1_60',
							'56','15743_1_11','16002_1_5',
							'62','15777_1_41','15977_1_69',
							'7','16057_1_2','15834_1_58',
							'1','15861_1_13','15861_1_18',
							'111','16057_1_2','15834_1_58',
							'7','15861_1_13','15103_1_68',
							'8','15949_1_52','15714_1_60',
							'37','15950_1_6','15950_1_16',
							'10','16002_1_17','15758_1_66',
							'81','16035_1_3','15758_1_83',
							'105','16035_1_3','15758_1_83',
							'32','15861_1_13','15103_1_68'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp[, CLASS_OR:= 'OR_disconnected']	
	tmp2	<- as.data.table(matrix( c('66','15778_1_82','15758_1_76',
							'71','16005_1_55','16005_1_54',
							'11','16056_1_74','15736_1_34',				
							'10','15861_1_26','15080_1_12',
							'12','15861_1_20','15861_1_24',
							'64','16005_1_55','15097_1_82',
							'97','16005_1_55','15097_1_82',				
							'80','15080_1_3','15758_1_66'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp2[, CLASS_OR:= 'OR_ambiguous']
	tmp		<- rbind(tmp, tmp2)
	tmp2	<- as.data.table(matrix( c('110','15806_1_32','15806_1_27',
							'78','16016_1_43','16017_1_42',
							'62','15861_1_26','15861_1_22',
							'42','16056_1_73','15736_1_7',
							'102','15736_1_20','15880_1_7',
							'26','15805_1_22','15805_1_23',
							'98','15964_1_3','15758_1_84',
							'91','15103_1_74','15861_1_22',
							'93','15804_1_81','15804_1_60',
							'77','15804_1_53','15804_1_57',
							'16','15097_1_51','15805_1_23',
							'80','15103_1_74','15080_1_12'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp2[, CLASS_OR:= 'OR_transmission']
	tmp		<- rbind(tmp, tmp2)
	set(tmp, NULL, 'PTY_RUN', tmp[, as.integer(PTY_RUN)])
	rpi		<- unique(merge(rpi, tmp, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')))
	rps		<- merge(rps, rpi, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all.x=1)
	#
	#	get false pos and true pos for the two unambiguous groups	
	#
	stopifnot( nrow(subset(rps, (COUP_SC=='F->M'|COUP_SC=='M->F') & is.na(CLASS_OR)))==0 )
	tmp		<- subset(rps, (COUP_SC=='F->M'|COUP_SC=='M->F') & CLASS_OR=='OR_disconnected')
	rpb		<- tmp[, list(	POT_FP= sum(WIN_OF_TYPE_N[TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int']), 
					POT_TN=sum(WIN_OF_TYPE_N[TYPE=='disconnected'|TYPE=='unint']),
					POT_TP=0,
					POT_FN=0), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(rps, (COUP_SC=='F->M'|COUP_SC=='M->F') & CLASS_OR=='OR_transmission')
	tmp		<- tmp[, list(	POT_TP= sum(WIN_OF_TYPE_N[TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int']), 
					POT_FN=sum(WIN_OF_TYPE_N[TYPE=='disconnected'|TYPE=='unint']),
					POT_FP=0,
					POT_TN=0), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpb		<- rbind(rpb, tmp, fill=TRUE, use.names=TRUE)
	rpb[, TOTAL:= POT_FP+POT_FN+POT_TP+POT_TN]
	#	by choice of the three groups, there are only TN and TPs.. ..so this is not interesting here.
	rpb[, Maj_Class:= NA_character_]
	set(rpb, rpb[, which(POT_FP<POT_TN)], 'Maj_Class', 'TN')
	set(rpb, rpb[, which(POT_FP>=POT_TN)], 'Maj_Class', 'FP')
	set(rpb, rpb[, which(POT_TP>POT_FN)], 'Maj_Class', 'TP')
	set(rpb, rpb[, which(POT_TP<=POT_FN)], 'Maj_Class', 'FN')					
	#
	#
	tmp		<- melt(rpb, id.vars=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC','TOTAL'), measure.vars=c('POT_FP','POT_FN'))	
	tmp		<- tmp[, list(Potential_False_Class= sum(value), Potential_True_Class= sum(TOTAL)-sum(value) ), by=c('RUN','variable')]	
	set(tmp, tmp[,which(variable=='POT_FP')], 'variable', 'Manually identified unlinked pair')
	set(tmp, tmp[,which(variable=='POT_FN')], 'variable', 'Manually identified trm pair')	
	tmp[, Total_Class:= Potential_False_Class+Potential_True_Class]
	
	ggplot( tmp, aes(x=RUN, y=Potential_False_Class/Total_Class, fill=variable) ) + 
			geom_bar(stat='identity') + 
			scale_fill_brewer(palette='Set1') +
			scale_y_continuous(breaks=seq(0,1,0.02), labels=percent) +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~variable) +
			coord_flip() +
			labs(	y= '\nproportion of "likely false" window assignments\ne.g. % of windows from manually identified trm pairs that are falsely classified as unlinked', 
					x='run\n',
					colour='',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')	
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trmpair_windows_falsepos.pdf',sep='')), w=10, h=5)
	#	number of couples with less than 5 trm assignments
	tmp		<- unique(subset(rpb, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID)))
	tmp2	<- unique(subset(rpb, select=RUN))
	tmp[, D:=1]
	tmp2[, D:=1]
	tmp		<- merge(tmp, tmp2, by='D', allow.cartesian=TRUE)
	tmp[, D:=NULL]
	tmp		<- merge(tmp, rpb, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all.x=1)
	set(tmp, tmp[, which(is.na(TOTAL))], 'TOTAL', 0)	
	tmp		<- subset(tmp, TOTAL<5)[, list(WIN_L_FIVE=length(PAIR_ID)), by='RUN']	
	setkey(tmp, RUN)
	#	also track pairs that cannot be any longer classified
	
	#
	#	plot number of directional trm assignments in only possible direction 
	#	
	rpa		<- subset(rps, COUP_SC=='F->M')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='trans_fm']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='trans_mf'])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(rps, COUP_SC=='M->F')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='trans_mf']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='trans_fm'])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= IN_DIR+AGAINST_DIR]
	rpa		<- melt(rpa, measure.vars=c('IN_DIR', 'AGAINST_DIR'))
	set(rpa, rpa[, which(variable=='IN_DIR')], 'variable', 'trm assignment in the only epidemiologically possible direction')
	set(rpa, rpa[, which(variable=='AGAINST_DIR')], 'variable', 'trm assignment against the only epidemiologically possible direction')	
	setkey(rpa, PAIR_ID)		
	tmp		<- unique(subset(rpa, RUN==rpa$RUN[1]))[order(-WIN_TRM),][, list(COUPID=COUPID, MALE_SANGER_ID=MALE_SANGER_ID, FEMALE_SANGER_ID=FEMALE_SANGER_ID, PTY_RUN=PTY_RUN, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PTY_RUN','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set1') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~RUN) +
			coord_flip() +
			labs(	x= '\nsequence pairs of couples', 
					y='number of transmission windows with direction\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trmdir_windows.pdf',sep='')), w=35, h=10)
	#
	#	get false pos and true pos for the manual classification group
	#
	tmp		<- unique(subset(rps, RUN=='RCCS_161107_w270_d20_r004_mr100' & CLASS_OR=='OR_transmission',c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')))
	rpa		<- merge(rpa, tmp, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	tmp		<- rpa[, {
				list( TRM_OK_MAJ= as.integer(sum(value[variable=='trm assignment in the only epidemiologically possible direction'])>sum(value[variable!='trm assignment in the only epidemiologically possible direction'])), 
						TRM_OK_66pc= as.integer(sum(value[variable=='trm assignment in the only epidemiologically possible direction']) / sum(value) >= 2/3))
			}, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')]
	rpa		<- merge(rpa, tmp, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	tmp		<- rpa[, {
				tmp	<- sum(value[variable=='trm assignment in the only epidemiologically possible direction'])				
				list(	V=	c( tmp / (sum(WIN_TRM)/2), 1 - tmp / (sum(WIN_TRM)/2), (sum(TRM_OK_MAJ)/2) / (length(TRM_OK_MAJ)/2), 1 - (sum(TRM_OK_MAJ)/2) / (length(TRM_OK_MAJ)/2), (sum(TRM_OK_66pc)/2) / (length(TRM_OK_66pc)/2), 1-(sum(TRM_OK_66pc)/2) / (length(TRM_OK_66pc)/2)),
						TYPE= c('TP', 'FP', 'TP', 'FP', 'TP', 'FP'),
						RULE= c('by window', 'by window', 'by majority rule', 'by majority rule', 'by 66% majority', 'by 66% majority')
				)
			}, by=c('RUN')]	
	ggplot( tmp, aes(x=RUN, y=V, fill= TYPE) ) + 
			geom_bar(stat='identity', position='stack') + 
			scale_fill_brewer(palette='Set1') +
			scale_y_continuous(breaks=seq(0,1,0.1), labels=percent) +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			coord_flip() +
			labs(	y= '\nproportion of "false positives" and "false negatives"', 
					x='run\n',
					colour='',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n') +
			facet_grid(~RULE)
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trm_windows_falsepos.pdf',sep='')), w=14, h=6)
	
	tmp		<- unique(subset(rpa, !TRM_OK_66pc, c('RUN','PTY_RUN','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID')))
	setkey(tmp, RUN)
	#
	#
	#
	#
}

RakaiCouples.analyze.couples.161110.trms<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/OR_Work/2016/2016_Rakai_Couples/161110/RCCS_161110_w270_trms_assignments.rda')
	rps		<- subset(rps, RUN%in%c("RCCS_161110_w270_d200_r004_mr1_mt1", "RCCS_161110_w270_d50_r004_mr1_mt1", "RCCS_161110_w270_d20_r004_mr1_mt1", "RCCS_161110_w270_d50_r004_mr20_mt2", "RCCS_161110_w270_d50_r004_mr30_mt2", "RCCS_161110_w270_d50_r004_mr40_mt2", "RCCS_161110_w270_d50_r004_mr50_mt2", "RCCS_161110_w270_d50_r004_mr100_mt2") )
	#
	rpsn	<- copy(rps)
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161007/RCCS_161007_w270_assignments.rda')	
	rps		<- rbind(rpsn, rps, fill=TRUE, use.names=TRUE)	
	
	set(rps, NULL, 'RUN', rps[, factor(RUN, levels=c( 	"RCCS_161007_w270", "RCCS_161110_w270_d20_r004_mr1_mt1","RCCS_161110_w270_d50_r004_mr1_mt1","RCCS_161110_w270_d200_r004_mr1_mt1",
														"RCCS_161110_w270_d50_r004_mr20_mt2","RCCS_161110_w270_d50_r004_mr30_mt2","RCCS_161110_w270_d50_r004_mr40_mt2","RCCS_161110_w270_d50_r004_mr50_mt2","RCCS_161110_w270_d50_r004_mr100_mt2"))])
	rps[, table(RUN)]	
	rps		<- subset(rps, COUPID!='F108560:F108560')
	rps[, WIN_OF_TYPE_N:= WIN_OF_TYPE_P*WIN_TOTAL]
	
	
	#
	#
	#	for each run: plot pairs	
	#	
	run		<- 'RCCS_161110_w270_dxxx'
	dir		<- rps$DIR[1]
	rpp		<- subset(rps, RUN==rps$RUN[1])
	setkey(rpp, MALE_SANGER_ID, FEMALE_SANGER_ID, PAIR_ID)
	tmp		<- unique(rpp)	
	tmp[, PLOT_ID:= as.numeric(gsub('-','\\.',PAIR_ID))]
	tmp		<- tmp[order(-PLOT_ID),]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair', PAIR_ID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- merge(subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL)), rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('trans_mf','trans_fm','unint','int','cher','disconnected'), labels=c('M transmit to F','F transmit to M','M, F are unint','M, F are intermingled','M, F are a cherry','M, F are disconnected'))])
	tmp		<- subset(tmp, select=c(RUN, PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, WIN_OF_TYPE_P, WIN_TOTAL, TYPE, LABEL))
	tmp[, WIN_OF_TYPE_N:=WIN_OF_TYPE_P*WIN_TOTAL]
	setkey(tmp, COUP_SC, LABEL, RUN, TYPE)
	#	seroinc
	tmp2	<- subset(tmp, COUP_SC=='seroinc')
	ggplot(tmp2, aes(x=RUN, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M transmit to F'="#9E0142",'F transmit to M'="#F46D43",'M, F are a cherry'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are unint'="blue",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC+LABEL, ncol=1) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_SeroInc.pdf',sep='')), w=25, h=max(4,0.23*nrow(tmp2)), limitsize = FALSE)
	#
	#	re-examine phylogenies for all sero-discordant couples
	#	
	tmp		<- subset(rps, COUP_SC=='seroinc' & RUN=='RCCS_161110_w270_d50_r004_mr1_mt1')
	setkey(tmp, RUN, PAIR_ID, MALE_SANGER_ID, FEMALE_SANGER_ID)
	dir		<- tmp[1,DIR]
	cat('\ndir is',dir,'\trun is',run)			
	setkey(tmp, MALE_SANGER_ID, FEMALE_SANGER_ID, PAIR_ID)				
	rpoku	<- unique(tmp)
	load( tmp[1, gsub('trmsout','phscout',FILE)] )	#loads phs dtrms dtrees
	invisible(sapply(seq_len(nrow(rpoku)), function(ii)
					{								
						pair.id		<- rpoku[ii, PAIR_ID]
						pty.run		<- rpoku[ii, PTY_RUN]
						dfs			<- subset(dtrees, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
						dfs[, TITLE:= dfs[, paste('pair', pair.id, '\n', rpoku[ii, COUP_SC], '\nid M: ', rpoku[ii, MALE_RID], ' (', rpoku[ii, MALE_SANGER_ID], ')\nid F: ', rpoku[ii, FEMALE_RID], ' (', rpoku[ii, FEMALE_SANGER_ID], ')\nrun ', pty.run, '\nwindow ', W_FROM,'-', W_TO,sep='')]]			
						plot.file	<- file.path(dir, paste(run,'-phsc-serodiscpairs-',rpoku[ii, COUP_SC],'-M-', rpoku[ii, MALE_RID],'-F-',rpoku[ii, FEMALE_RID],'-', pair.id,'.pdf',sep=''))			
						invisible(phsc.plot.selected.pairs(phs, dfs, rpoku[ii, MALE_SANGER_ID], rpoku[ii, FEMALE_SANGER_ID], plot.file=plot.file, pdf.h=150, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40))
					}))	
	#
	#	plot evidence for or against transmission
	#	
	rpa		<- subset(rps, COUP_SC=='seroinc')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int']), CHER=WIN_OF_TYPE_N[TYPE=='cher'], NO_ANC=sum(WIN_OF_TYPE_N[TYPE=='disconnected'|TYPE=='unint'])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa[, WIN_TRM:= ANY_ANC]
	rpa		<- melt(rpa, measure.vars=c('ANY_ANC', 'NO_ANC','CHER'))
	set(rpa, rpa[, which(variable=='ANY_ANC')], 'variable', 'transmission assignments or intermingled')
	set(rpa, rpa[, which(variable=='NO_ANC')], 'variable', 'unint or disconnected')
	set(rpa, rpa[, which(variable=='CHER')], 'variable', 'cherry')
	setkey(rpa, PAIR_ID)	
	tmp		<- unique(subset(rpa, RUN=='RCCS_161110_w270_d50_r004_mr1_mt1'))[order(-WIN_TRM),][, list(COUPID=COUPID, MALE_SANGER_ID=MALE_SANGER_ID, FEMALE_SANGER_ID=FEMALE_SANGER_ID, PTY_RUN=PTY_RUN, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PTY_RUN','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set2') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~RUN) +
			coord_flip() +
			labs(	x= '\nsequence pairs of couples', 
					y='number of windows\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS seroincident couples\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-seroincpairs-number_trmpair_windows.pdf',sep='')), w=35, h=10)
	#
	#	define three transmission evidence groups manually
	#	
	rpi	<- subset(rpa, RUN=='RCCS_161110_w270_d50_r004_mr1_mt1', c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID))
	setkey(rpi, CLASS_RANK)
	tmp	<- as.data.table(matrix( c('83','15715_1_17','15715_1_18',
							'36','15715_1_17','15715_1_18',
							'85','15699_1_38','15699_1_40',
							'104','15699_1_3','15699_1_4',
							'106','15736_1_28','15102_1_16',
							'121','16059_1_66','16059_1_52',
							'31','15981_1_31','15981_1_33',
							'98','15916_1_7','15916_1_21',
							'64','15916_1_7','15916_1_21',
							'122','15910_1_76','15910_1_82',
							'24','15981_1_17','15835_1_32',
							'70','15833_1_84','15833_1_83',
							'104','16059_1_66','16059_1_52',
							'24','16016_1_11','16017_1_44',
							'24','15915_1_51','15776_1_32',
							'15','16033_1_76','16033_1_60',
							'7','16033_1_76','16033_1_60',
							'18','16019_1_45','16018_1_44',
							'67','15892_1_3','15892_1_1',
							'91','16021_1_66','16021_1_58',
							'29','16021_1_66','16021_1_58',
							'116','15915_1_84','15915_1_83',
							'64','15915_1_84','15915_1_83',
							'50','15115_1_4','15115_1_9',
							'86','15978_1_12','15777_1_42',
							'68','15978_1_12','15777_1_42',
							'44','16056_1_83','15099_1_59'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp[, CLASS_OR:= 'OR_disconnected']	
	tmp2	<- as.data.table(matrix( c('111','15950_1_27','15950_1_25',
							'122','16056_1_58','16056_1_57',
							'38','15105_1_35','16016_1_4',
							'113','15950_1_5','15958_1_47',				
							'99','15915_1_61','15915_1_74',
							'117','15915_1_61','15915_1_74',
							'60','16056_1_83','16056_1_85',
							'39','16034_1_86','16034_1_82',
							'120','15981_1_19','15981_1_23',				
							'79','15978_1_36','15978_1_41'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp2[, CLASS_OR:= 'OR_ambiguous']
	tmp		<- rbind(tmp, tmp2)
	tmp2	<- as.data.table(matrix( c('101','16003_1_76','16003_1_79',
							'37','15114_1_64','15758_1_73',
							'37','16016_1_5','16016_1_4',
							'20','16057_1_15','16057_1_16',
							'59','15699_1_5','15867_1_21',
							'99','15915_1_78','15915_1_64',
							'109','15777_1_29','15777_1_21',
							'63','16003_1_86','16003_1_85',
							'109','15776_1_19','15776_1_28',
							'27','16219_1_78','16219_1_80',
							'120','16059_1_73','16059_1_53',
							'25','15675_1_87','15675_1_77',
							'70','15835_1_21','15835_1_22'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp2[, CLASS_OR:= 'OR_transmission']
	tmp		<- rbind(tmp, tmp2)
	set(tmp, NULL, 'PTY_RUN', tmp[, as.integer(PTY_RUN)])
	rpi		<- unique(merge(rpi, tmp, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')))
	rps		<- merge(rps, rpi, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all.x=1)
	#
	#	get false pos and true pos for the two unambiguous groups	
	#
	stopifnot( nrow(subset(rps, (COUP_SC=='seroinc') & is.na(CLASS_OR)))==0 )
	tmp		<- subset(rps, COUP_SC=='seroinc' & CLASS_OR=='OR_disconnected')
	rpb		<- tmp[, list(	POT_FP= sum(WIN_OF_TYPE_N[TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int']), 
					POT_TN=sum(WIN_OF_TYPE_N[TYPE=='disconnected'|TYPE=='unint']),
					POT_TP=0,
					POT_FN=0), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(rps, COUP_SC=='seroinc' & CLASS_OR=='OR_transmission')
	tmp		<- tmp[, list(	POT_TP= sum(WIN_OF_TYPE_N[TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int']), 
					POT_FN=sum(WIN_OF_TYPE_N[TYPE=='disconnected'|TYPE=='unint']),
					POT_FP=0,
					POT_TN=0), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpb		<- rbind(rpb, tmp, fill=TRUE, use.names=TRUE)
	rpb[, TOTAL:= POT_FP+POT_FN+POT_TP+POT_TN]
	#	by choice of the three groups, there are only TN and TPs.. ..so this is not interesting here.
	rpb[, Maj_Class:= NA_character_]
	set(rpb, rpb[, which(POT_FP<POT_TN)], 'Maj_Class', 'TN')
	set(rpb, rpb[, which(POT_FP>=POT_TN)], 'Maj_Class', 'FP')
	set(rpb, rpb[, which(POT_TP>POT_FN)], 'Maj_Class', 'TP')
	set(rpb, rpb[, which(POT_TP<=POT_FN)], 'Maj_Class', 'FN')					
	#
	#
	tmp		<- melt(rpb, id.vars=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC','TOTAL'), measure.vars=c('POT_FP','POT_FN'))	
	tmp		<- tmp[, list(Potential_False_Class= sum(value), Potential_True_Class= sum(TOTAL)-sum(value) ), by=c('RUN','variable')]	
	set(tmp, tmp[,which(variable=='POT_FP')], 'variable', 'Manually identified unlinked pair')
	set(tmp, tmp[,which(variable=='POT_FN')], 'variable', 'Manually identified trm pair')	
	tmp[, Total_Class:= Potential_False_Class+Potential_True_Class]
	
	ggplot( tmp, aes(x=RUN, y=Potential_False_Class/Total_Class, fill=variable) ) + 
			geom_bar(stat='identity') + 
			scale_fill_brewer(palette='Set1') +
			scale_y_continuous(breaks=seq(0,1,0.02), labels=percent) +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~variable) +
			coord_flip() +
			labs(	y= '\nproportion of "likely false" window assignments\ne.g. % of windows from manually identified trm pairs that are falsely classified as unlinked', 
					x='run\n',
					colour='',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')	
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trmpair_windows_falsepos.pdf',sep='')), w=10, h=5)
	#	number of couples with less than 5 trm assignments
	tmp		<- unique(subset(rpb, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID)))
	tmp2	<- unique(subset(rpb, select=RUN))
	tmp[, D:=1]
	tmp2[, D:=1]
	tmp		<- merge(tmp, tmp2, by='D', allow.cartesian=TRUE)
	tmp[, D:=NULL]
	tmp		<- merge(tmp, rpb, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all.x=1)
	set(tmp, tmp[, which(is.na(TOTAL))], 'TOTAL', 0)	
	tmp		<- subset(tmp, TOTAL<5)[, list(WIN_L_FIVE=length(PAIR_ID)), by='RUN']	
	setkey(tmp, RUN)
	#	also track pairs that cannot be any longer classified

	
}

RakaiCouples.analyze.couples.161213.trms<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161213/RCCS_161213_w270_trms_assignments.rda')
	rps		<- subset(rps, RUN%in%c("RCCS_161213_w270_d50_p001_mr20_mt2_cl1_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_cl2_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_cl4_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_clInf_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_clInf_dInf") )	
	set(rps, NULL, 'RUN', rps[, factor(RUN, levels=c( 	"RCCS_161213_w270_d50_p001_mr20_mt2_cl1_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_cl2_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_cl4_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_clInf_d7", "RCCS_161213_w270_d50_p001_mr20_mt2_clInf_dInf"))])
	rps[, table(RUN)]	
	rps		<- subset(rps, COUPID!='F108560:F108560')
	rps[, WIN_OF_TYPE_N:= WIN_OF_TYPE_P*WIN_TOTAL]
	#
	#
	#	for each run: plot pairs	
	#	
	run		<- 'RCCS_161213_w270_dxxx'
	dir		<- rps$DIR[1]
	rpp		<- subset(rps, RUN==rps$RUN[1])
	setkey(rpp, MALE_SANGER_ID, FEMALE_SANGER_ID, PAIR_ID)
	tmp		<- unique(rpp)	
	tmp[, PLOT_ID:= as.numeric(gsub('-','\\.',PAIR_ID))]
	tmp		<- tmp[order(-PLOT_ID),]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair', PAIR_ID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- merge(subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL)), rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, 	levels=c("close_anc_mf", "close_anc_fm", "close_cher_unint",'anc_mf','anc_fm','unint','int','cher','disconnected'), 
												labels=c('M close & ancestral to F','F close & ancestral to M','M, F are close & unint/cherry',
														'M ancestral to F','F ancestral to M','M, F are unint','M, F are intermingled','M, F are a cherry','M, F are disconnected'))])	
	tmp		<- subset(tmp, select=c(RUN, PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, WIN_OF_TYPE_P, WIN_TOTAL, TYPE, LABEL))
	tmp[, WIN_OF_TYPE_N:=WIN_OF_TYPE_P*WIN_TOTAL]
	setkey(tmp, COUP_SC, LABEL, RUN, TYPE)
	#	seroinc
	tmp2	<- subset(tmp, COUP_SC=='seroinc')
	cols	<- c(	'M close & ancestral to F'="#542788",
					'M ancestral to F'="#8073AC",
					'F close & ancestral to M'="#8C510A",
					'F ancestral to M'="#BF812D",
					'M, F are close & unint/cherry'="#B2182B",
					'M, F are a cherry'="#D6604D",
					'M, F are unint'="#F4A582",
					'M, F are intermingled'="#2166AC",						
					'M, F are disconnected'='grey50')
	ggplot(tmp2, aes(x=RUN, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=cols) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC+LABEL, ncol=1) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_SeroInc.pdf',sep='')), w=25, h=max(2.5,0.15*nrow(tmp2)), limitsize = FALSE)
	#
	#	plot evidence for or against transmission
	#	
	rpa		<- subset(rps, COUP_SC=='seroinc')[, list(ANY_CLOSE= sum(WIN_OF_TYPE_N[grepl('close',TYPE)]), ANY_ANC=sum(WIN_OF_TYPE_N[grepl('anc',TYPE)]), CHERUNINT=WIN_OF_TYPE_N[TYPE=='cher'|TYPE=='unint'], NO_ANC=sum(WIN_OF_TYPE_N[TYPE=='disconnected'])), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa[, WIN_TRM:= ANY_ANC]
	rpa		<- melt(rpa, measure.vars=c('ANY_CLOSE','ANY_ANC', 'NO_ANC','CHERUNINT'))
	set(rpa, rpa[, which(variable=='ANY_CLOSE')], 'variable', 'close')
	set(rpa, rpa[, which(variable=='ANY_ANC')], 'variable', 'ancestral')
	set(rpa, rpa[, which(variable=='NO_ANC')], 'variable', 'disconnected')
	set(rpa, rpa[, which(variable=='CHERUNINT')], 'variable', 'cherry/unint')
	setkey(rpa, PAIR_ID)	
	tmp		<- unique(subset(rpa, RUN=='RCCS_161213_w270_d50_p001_mr20_mt2_cl1_d7'))[order(-WIN_TRM),][, list(COUPID=COUPID, MALE_SANGER_ID=MALE_SANGER_ID, FEMALE_SANGER_ID=FEMALE_SANGER_ID, PTY_RUN=PTY_RUN, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PTY_RUN','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID'))	
	setkey(rpa, variable, CLASS_RANK)
	cols	<- c(	'close'="#542788",
					'ancestral'="#8073AC",
					'cherry/unint'="#F4A582",									
					'disconnected'='grey50')	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_manual(values=cols) +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~RUN) +
			coord_flip() +
			labs(	x= '\nsequence pairs of couples', 
					y='number of windows\n',
					colour='phyloscanner\nassignments',
					title='\nphyloscanner assignments\nto RCCS seroincident couples\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-seroincpairs-number_trmpair_windows.pdf',sep='')), w=35, h=10)
	#
	#	define three transmission evidence groups manually
	#	
	rpi	<- subset(rpa, RUN=='RCCS_161110_w270_d50_r004_mr1_mt1', c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID))
	setkey(rpi, CLASS_RANK)
	tmp	<- as.data.table(matrix( c('83','15715_1_17','15715_1_18',
							'36','15715_1_17','15715_1_18',
							'85','15699_1_38','15699_1_40',
							'104','15699_1_3','15699_1_4',
							'106','15736_1_28','15102_1_16',
							'121','16059_1_66','16059_1_52',
							'31','15981_1_31','15981_1_33',
							'98','15916_1_7','15916_1_21',
							'64','15916_1_7','15916_1_21',
							'122','15910_1_76','15910_1_82',
							'24','15981_1_17','15835_1_32',
							'70','15833_1_84','15833_1_83',
							'104','16059_1_66','16059_1_52',
							'24','16016_1_11','16017_1_44',
							'24','15915_1_51','15776_1_32',
							'15','16033_1_76','16033_1_60',
							'7','16033_1_76','16033_1_60',
							'18','16019_1_45','16018_1_44',
							'67','15892_1_3','15892_1_1',
							'91','16021_1_66','16021_1_58',
							'29','16021_1_66','16021_1_58',
							'116','15915_1_84','15915_1_83',
							'64','15915_1_84','15915_1_83',
							'50','15115_1_4','15115_1_9',
							'86','15978_1_12','15777_1_42',
							'68','15978_1_12','15777_1_42',
							'44','16056_1_83','15099_1_59'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp[, CLASS_OR:= 'OR_disconnected']	
	tmp2	<- as.data.table(matrix( c('111','15950_1_27','15950_1_25',
							'122','16056_1_58','16056_1_57',
							'38','15105_1_35','16016_1_4',
							'113','15950_1_5','15958_1_47',				
							'99','15915_1_61','15915_1_74',
							'117','15915_1_61','15915_1_74',
							'60','16056_1_83','16056_1_85',
							'39','16034_1_86','16034_1_82',
							'120','15981_1_19','15981_1_23',				
							'79','15978_1_36','15978_1_41'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp2[, CLASS_OR:= 'OR_ambiguous']
	tmp		<- rbind(tmp, tmp2)
	tmp2	<- as.data.table(matrix( c('101','16003_1_76','16003_1_79',
							'37','15114_1_64','15758_1_73',
							'37','16016_1_5','16016_1_4',
							'20','16057_1_15','16057_1_16',
							'59','15699_1_5','15867_1_21',
							'99','15915_1_78','15915_1_64',
							'109','15777_1_29','15777_1_21',
							'63','16003_1_86','16003_1_85',
							'109','15776_1_19','15776_1_28',
							'27','16219_1_78','16219_1_80',
							'120','16059_1_73','16059_1_53',
							'25','15675_1_87','15675_1_77',
							'70','15835_1_21','15835_1_22'), ncol=3, byrow=TRUE, dimnames=list(c(),c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))))
	tmp2[, CLASS_OR:= 'OR_transmission']
	tmp		<- rbind(tmp, tmp2)
	set(tmp, NULL, 'PTY_RUN', tmp[, as.integer(PTY_RUN)])
	rpi		<- unique(merge(rpi, tmp, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID')))
	rps		<- merge(rps, rpi, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all.x=1)
	#
	#	get false pos and true pos for the two unambiguous groups	
	#
	stopifnot( nrow(subset(rps, (COUP_SC=='seroinc') & is.na(CLASS_OR)))==0 )
	tmp		<- subset(rps, COUP_SC=='seroinc' & CLASS_OR=='OR_disconnected')
	rpb		<- tmp[, list(	POT_FP= sum(WIN_OF_TYPE_N[TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int']), 
					POT_TN=sum(WIN_OF_TYPE_N[TYPE=='disconnected'|TYPE=='unint']),
					POT_TP=0,
					POT_FN=0), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(rps, COUP_SC=='seroinc' & CLASS_OR=='OR_transmission')
	tmp		<- tmp[, list(	POT_TP= sum(WIN_OF_TYPE_N[TYPE=='trans_fm'|TYPE=='trans_mf'|TYPE=='int']), 
					POT_FN=sum(WIN_OF_TYPE_N[TYPE=='disconnected'|TYPE=='unint']),
					POT_FP=0,
					POT_TN=0), by=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpb		<- rbind(rpb, tmp, fill=TRUE, use.names=TRUE)
	rpb[, TOTAL:= POT_FP+POT_FN+POT_TP+POT_TN]
	#	by choice of the three groups, there are only TN and TPs.. ..so this is not interesting here.
	rpb[, Maj_Class:= NA_character_]
	set(rpb, rpb[, which(POT_FP<POT_TN)], 'Maj_Class', 'TN')
	set(rpb, rpb[, which(POT_FP>=POT_TN)], 'Maj_Class', 'FP')
	set(rpb, rpb[, which(POT_TP>POT_FN)], 'Maj_Class', 'TP')
	set(rpb, rpb[, which(POT_TP<=POT_FN)], 'Maj_Class', 'FN')					
	#
	#
	tmp		<- melt(rpb, id.vars=c('RUN','PTY_RUN','COUPID','PAIR_ID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC','TOTAL'), measure.vars=c('POT_FP','POT_FN'))	
	tmp		<- tmp[, list(Potential_False_Class= sum(value), Potential_True_Class= sum(TOTAL)-sum(value) ), by=c('RUN','variable')]	
	set(tmp, tmp[,which(variable=='POT_FP')], 'variable', 'Manually identified unlinked pair')
	set(tmp, tmp[,which(variable=='POT_FN')], 'variable', 'Manually identified trm pair')	
	tmp[, Total_Class:= Potential_False_Class+Potential_True_Class]
	
	ggplot( tmp, aes(x=RUN, y=Potential_False_Class/Total_Class, fill=variable) ) + 
			geom_bar(stat='identity') + 
			scale_fill_brewer(palette='Set1') +
			scale_y_continuous(breaks=seq(0,1,0.02), labels=percent) +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~variable) +
			coord_flip() +
			labs(	y= '\nproportion of "likely false" window assignments\ne.g. % of windows from manually identified trm pairs that are falsely classified as unlinked', 
					x='run\n',
					colour='',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')	
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trmpair_windows_falsepos.pdf',sep='')), w=10, h=5)
	#	number of couples with less than 5 trm assignments
	tmp		<- unique(subset(rpb, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID)))
	tmp2	<- unique(subset(rpb, select=RUN))
	tmp[, D:=1]
	tmp2[, D:=1]
	tmp		<- merge(tmp, tmp2, by='D', allow.cartesian=TRUE)
	tmp[, D:=NULL]
	tmp		<- merge(tmp, rpb, by=c('RUN','PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'), all.x=1)
	set(tmp, tmp[, which(is.na(TOTAL))], 'TOTAL', 0)	
	tmp		<- subset(tmp, TOTAL<5)[, list(WIN_L_FIVE=length(PAIR_ID)), by='RUN']	
	setkey(tmp, RUN)
	#	also track pairs that cannot be any longer classified
	
	
}

RakaiCouples.analyze.couples.161219.trms<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	# load pty.run
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	# load couples "rp"
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	# load transmission summaries
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/161219/RCCS_161219_w270_trms_assignments.rda')
	rps		<- subset(rps, RUN%in%c("RCCS_161219_w270_d50_p001_mr20_mt1_cl1_d5", "RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5", "RCCS_161219_w270_d50_p001_mr20_mt1_cl1_dInf", "RCCS_161219_w270_d50_p001_mr20_mt1_cl2_dInf", "RCCS_161219_w270_d50_p001_mr20_mt1_clInf_dInf") )	
	set(rps, NULL, 'RUN', rps[, factor(RUN, levels=c( 	"RCCS_161219_w270_d50_p001_mr20_mt1_cl1_d5", "RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5", "RCCS_161219_w270_d50_p001_mr20_mt1_cl1_dInf", "RCCS_161219_w270_d50_p001_mr20_mt1_cl2_dInf", "RCCS_161219_w270_d50_p001_mr20_mt1_clInf_dInf"))])
	rps[, table(RUN)]	
	rps		<- subset(rps, COUPID!='F108560:F108560')
	rps[, WIN_OF_TYPE_N:= WIN_OF_TYPE_P*WIN_TOTAL]
	#
	#	make labels	
	#	
	run		<- 'RCCS_161219_w270_dxxx'
	dir		<- rps$DIR[1]
	rpp		<- subset(rps, RUN==rps$RUN[1])	
	tmp		<- unique(rpp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID','PAIR_ID'))	
	tmp[, PLOT_ID:= as.numeric(gsub('-','\\.',PAIR_ID))]
	tmp		<- tmp[order(-PLOT_ID),]
	tmp[, LABEL_SH:= factor(PLOT_ID, levels=PLOT_ID, labels=paste(COUPID, ' ( M:', MALE_SANGER_ID,' F:',FEMALE_SANGER_ID, ' run:', PTY_RUN, ' )', sep=''))]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair', PAIR_ID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	rps		<- merge(subset(tmp, select=c(PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, LABEL, LABEL_SH)), rps, by=c('PTY_RUN','MALE_SANGER_ID','FEMALE_SANGER_ID'))
	rps		<- subset(rps, select=c(RUN, PTY_RUN, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, WIN_OF_TYPE_P, WIN_TOTAL, TYPE, TYPE_PAIR, LABEL, LABEL_SH))
	rps[, WIN_OF_TYPE_N:=WIN_OF_TYPE_P*WIN_TOTAL]
	setkey(rps, COUP_SC, LABEL, RUN, TYPE)
	#
	#	select seroincident couples	
	#		
	rps		<- subset(rps, COUP_SC=='seroinc')
	#
	#	define colours for TYPE and TYPE_PAIR	
	#		
	cols.type	<- do.call('rbind',list(
					data.table(	TYPE= c("chain_fm_nointermediate_close","chain_fm_nointermediate","chain_fm_nointermediate_distant"),
							COLS= brewer.pal(11, 'PiYG')[c(1,2,4)]),
					data.table(	TYPE= c("chain_mf_nointermediate_close","chain_mf_nointermediate","chain_mf_nointermediate_distant"),
							COLS= brewer.pal(11, 'PuOr')[c(1,2,4)]),
					data.table(	TYPE= c("intermingled_nointermediate_close","intermingled_nointermediate","intermingled_nointermediate_distant"),
							COLS= brewer.pal(11, 'PRGn')[c(1,2,4)]),
					data.table(	TYPE= c("chain_fm_withintermediate_close","chain_fm_withintermediate","chain_fm_withintermediate_distant"),
							COLS= rev(brewer.pal(11, 'BrBG'))[c(3,4,5)]),
					data.table(	TYPE= c("chain_mf_withintermediate_close","chain_mf_withintermediate","chain_mf_withintermediate_distant"),
							COLS= rev(brewer.pal(11, 'PRGn'))[c(3,4,5)]),
					data.table(	TYPE= c("intermingled_withintermediate_close","intermingled_withintermediate","intermingled_withintermediate_distant"),
							COLS= rev(brewer.pal(11, 'RdBu'))[c(3,4,5)]),
					data.table(	TYPE= c("other_close","other","other_distant"),
							COLS= rev(brewer.pal(11, 'RdGy'))[c(3,4,5)])))
	cols.type	<- { tmp<- cols.type[, COLS]; names(tmp) <- cols.type[, TYPE]; tmp }
	cols.typep	<- do.call('rbind',list(
					data.table(	TYPE= c("pair_close","pair","pair_distant"),
							COLS= brewer.pal(11, 'PuOr')[c(1,2,4)]),
					data.table(	TYPE= c("withintermediate_close","withintermediate","withintermediate_distant"),
							COLS= rev(brewer.pal(11, 'RdBu'))[c(3,4,5)]),
					data.table(	TYPE= c("other_close","other","other_distant"),
							COLS= rev(brewer.pal(11, 'RdGy'))[c(3,4,5)])))
	cols.typep	<- { tmp<- cols.typep[, COLS]; names(tmp) <- cols.typep[, TYPE]; tmp }
	#
	#	plot all types
	#
	ggplot(rps, aes(x=RUN, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=cols.type) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC+LABEL, ncol=1) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-seroincpairs-numberwindows_alltypes_long.pdf',sep='')), w=20, h=max(2.5,0.05*nrow(rps)), limitsize = FALSE)
	ggplot(rps, aes(x=LABEL_SH, y=WIN_OF_TYPE_N, fill=TYPE)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_manual(values=cols.type) +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~RUN) +
			coord_flip() +
			guides(fill=guide_legend(ncol=7)) +
			labs(	x= '\nsequence pairs of couples', 
					y='number of windows\n',
					colour='phyloscanner\nassignments',
					title='\nphyloscanner assignments\nto RCCS seroincident couples\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-seroincpairs-numberwindows_alltypes.pdf',sep='')), w=35, h=10)
	#
	#	plot pair types
	#	
	tmp	<- rps[, list(WIN_OF_TYPE_N=sum(WIN_OF_TYPE_N)), by=c('RUN','TYPE_PAIR','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC','LABEL','LABEL_SH')]
	ggplot(tmp, aes(x=RUN, y=WIN_OF_TYPE_N, fill=TYPE_PAIR)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology_distance') +
			scale_fill_manual(values=cols.typep) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC+LABEL, ncol=1) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-seroincpairs-numberwindows_pairtypes_long.pdf',sep='')), w=20, h=max(2.5,0.07*nrow(tmp)), limitsize = FALSE)	
	ggplot(tmp, aes(x=LABEL_SH, y=WIN_OF_TYPE_N, fill=TYPE_PAIR)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_manual(values=cols.typep) +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			facet_grid(~RUN) +
			coord_flip() +
			guides(fill=guide_legend(ncol=3)) +
			labs(	x= '\nsequence pairs of couples', 
					y='number of windows\n',
					fill='topology_distance',
					title='\nphyloscanner assignments\nto RCCS seroincident couples\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-seroincpairs-numberwindows_pairtypes.pdf',sep='')), w=35, h=10)
	#
	#	the contingency table is weird because different patients have different numbers of windows
	#
	tmp	<- subset(rps, RUN=='RCCS_161219_w270_d50_p001_mr20_mt1_cl2_d5')
}


RakaiCouples.process.couples.160930<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	load( "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda" )
	
	load("~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")
	setkey(rp, COUPID)
	unique(rp)[, table(COUP_SC)]
	#	total of couples with assigned SANGER_ID in RCCS
	#   F->M    M->F seroinc seropos 
	# 	10      19      52     266 	
	merge(unique(rp), unique(subset(pty.runs, COUPID!='Other', COUPID)), by='COUPID')[, table(COUP_SC)]
	#	total of couples with SANGER_ID for which I have data
	#	F->M    M->F seroinc seropos 
	#	7       16      45     235 
	#
	#	collect runs
	#
	infiles	<- data.table(	FILE= c(	'~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/160930/RCCS_160930_w270_phscout.rda' ))	
	infiles[, DIR:= dirname(FILE)]
	infiles[, RUN:= gsub('_phscout.rda','',basename(FILE))]				
	#
	#	for each run: get list of pairs
	#	
	rpso		<- infiles[, {
				#F_TRM	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/RCCS_160919_w270_trmStats.rda'; F_PH	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/RCCS_160919_w270_trees.rda'
				load(FILE)	#loads phs dtrms dtrees
				#	select likely pairs -- these are ordered
				dtrms[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]				
				tmp		<- dcast.data.table(dtrms, PAIR_ID~TYPE, value.var='WIN_OF_TYPE_P')
				set(tmp, tmp[, which(is.na(disconnected))], 'disconnected', 0)
				set(tmp, tmp[, which(is.na(sib))], 'sib', 0)
				set(tmp, tmp[, which(is.na(int))], 'int', 0)
				set(tmp, tmp[, which(is.na(anc_12))], 'anc_12', 0)
				set(tmp, tmp[, which(is.na(anc_21))], 'anc_21', 0)
				tmp		<- melt.data.table(tmp, id.vars='PAIR_ID', value.name='WIN_OF_TYPE_P', variable.name='TYPE')
				dtrms	<- merge(unique(subset(dtrms, select=c(PAIR_ID, ID1, ID2, PTY_RUN, WIN_TOTAL, SCORE))), tmp, by='PAIR_ID')
				#	double the likely pairs (preserving the inferred direction), so that all get matched with the pairs in rp
				tmp		<- copy(dtrms)
				setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
				set(tmp, tmp[, which(TYPE=='anc_12')], 'TYPE', 'anc')
				set(tmp, tmp[, which(TYPE=='anc_21')], 'TYPE', 'anc_12')
				set(tmp, tmp[, which(TYPE=='anc')], 'TYPE', 'anc_21')
				dtrms	<- rbind(tmp, dtrms)
				#	first individual is always male	
				setnames(dtrms, c('ID1','ID2'), c('MALE_SANGER_ID','FEMALE_SANGER_ID'))
				set(dtrms, dtrms[, which(TYPE=='anc_12')], 'TYPE', 'anc_mf')
				set(dtrms, dtrms[, which(TYPE=='anc_21')], 'TYPE', 'anc_fm')
				set(dtrms, NULL, 'TYPE', dtrms[, as.character(TYPE)])
				dtrms	<- merge(dtrms, rp, by=c('MALE_SANGER_ID','FEMALE_SANGER_ID'))				
				dtrms
			}, by=c('RUN','DIR','FILE')]
	#
	cat('\nNumber of couples',rps[, length(unique(COUPID))])
	setkey(rps, RUN, COUPID, MALE_SANGER_ID, FEMALE_SANGER_ID)
	cat('\nNumber of sequence pairs from couples',nrow(unique(rps)))
	setkey(rps, RUN, COUPID, PAIR_ID)
	cat('\nNumber of pairings (including repeated sequence pairs)',nrow(unique(rps)))	
	#
	#	for each run: plot pairs	
	#	
	#for( run in rps[, unique(RUN)] )
	#{		
	run		<- 'RCCS_160919_w270'
	run		<- 'RCCS_160930_w270'
	dir		<- subset(rps, RUN==run)[1,DIR]
	cat('\ndir is',dir,'\trun is',run)
	df		<- subset(rps, RUN==run)
	setkey(df, MALE_SANGER_ID, FEMALE_SANGER_ID, PAIR_ID)					
	#
	#	plot evidence
	#		
	tmp		<- unique(df)
	tmp[, PLOT_ID:= as.numeric(gsub('-','\\.',PAIR_ID))]
	tmp		<- tmp[order(-PLOT_ID),]
	tmp[, LABEL:= factor(PLOT_ID, levels=PLOT_ID, labels=paste('Pair', PAIR_ID,' -type=', COUP_SC, ' -phsc.run=',PTY_RUN, '\nPerson M ', MALE_RID, ' ', MALE_SANGER_ID,' -loc:',MALE_REGION,',',MALE_COMM_NUM,',',MALE_HH_NUM,' -birth:',MALE_BIRTHDATE,' -neg:',MALE_LASTNEGDATE,' -pos:',MALE_FIRSTPOSDATE,' -seq:',MALE_SEQDATE,
							'\n<->', 
							'\nPerson F ', FEMALE_RID, ' ', FEMALE_SANGER_ID,' -loc:',FEMALE_REGION,',',FEMALE_COMM_NUM,',',FEMALE_HH_NUM,' -birth:',FEMALE_BIRTHDATE,' -neg:',FEMALE_LASTNEGDATE,' -pos:',FEMALE_FIRSTPOSDATE,' -seq:',FEMALE_SEQDATE,																				
							'\n',sep=''))]
	tmp		<- merge(subset(tmp, select=c(PAIR_ID, LABEL)), df, by='PAIR_ID')
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('anc_mf','anc_fm','sib','int','disconnected'), labels=c('M ancestral to F','F ancestral to M','M, F are siblings','M, F are intermingled','M, F are disconnected'))])
	tmp		<- subset(tmp, select=c(PAIR_ID, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, WIN_OF_TYPE_P, WIN_TOTAL, TYPE, LABEL))
	tmp[, WIN_OF_TYPE_N:=WIN_OF_TYPE_P*WIN_TOTAL]
	#tmp		<- melt(tmp, measure.vars=c('WIN_OF_TYPE_N', 'WIN_OF_TYPE_P'))
	#set(tmp, tmp[, which(variable=='WIN_OF_TYPE_N')],'variable','number of read windows')
	#set(tmp, tmp[, which(variable=='WIN_OF_TYPE_P')],'variable','proportion of read windows')	
	tmp2	<- subset(tmp, COUP_SC=='F->M')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M ancestral to F'="#9E0142",'F ancestral to M'="#F46D43",'M, F are siblings'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_wrap(~COUP_SC, ncol=2) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_F2M.pdf',sep='')), w=25, h=max(4,0.23*nrow(tmp2)), limitsize = FALSE)
	tmp2	<- subset(tmp, COUP_SC=='M->F')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M ancestral to F'="#9E0142",'F ancestral to M'="#F46D43",'M, F are siblings'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_grid(~COUP_SC) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_M2F.pdf',sep='')), w=25, h=max(3,0.23*nrow(tmp2)), limitsize = FALSE)
	tmp2	<- subset(tmp, COUP_SC=='seroinc')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M ancestral to F'="#9E0142",'F ancestral to M'="#F46D43",'M, F are siblings'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_grid(~COUP_SC) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_SeroInc.pdf',sep='')), w=25, h=max(3,0.23*nrow(tmp2)), limitsize = FALSE)
	tmp2	<- subset(tmp, COUP_SC=='seropos')
	ggplot(tmp2, aes(x=LABEL, y=WIN_OF_TYPE_N, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('M ancestral to F'="#9E0142",'F ancestral to M'="#F46D43",'M, F are siblings'="#ABDDA4",'M, F are intermingled'="#3288BD",'M, F are disconnected'='grey50')) +				
			theme_bw() + theme(legend.position='top') +
			facet_grid(~COUP_SC) +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs_SeroPos.pdf',sep='')), w=25, h=max(3,0.23*nrow(tmp2)), limitsize = FALSE)
	#
	#	tabulate correct classifications with phyloscanner for serodiscordant couples
	#
	#	correct: trm between M and F.
	rpa		<- subset(df, COUP_SC=='F->M' | COUP_SC=='M->F')[, list(CLASS='ancestral in either direction\nor intermingled', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='anc_mf'|TYPE=='anc_fm'|TYPE=='int'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(df, COUP_SC=='F->M')[, list(CLASS='ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='anc_fm'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp)
	tmp		<- subset(df, COUP_SC=='M->F')[, list(CLASS='ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='anc_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp)
	#	rank each couple
	tmp		<- rpa[order(CLASS, -CLASS_PROP),][, list(PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)), by='CLASS']
	rpa		<- merge(rpa, tmp, by=c('CLASS','PAIR_ID'))	
	setkey(rpa, CLASS, CLASS_RANK)
	#	plot by rank
	ggplot(rpa, aes(x=CLASS_PROP, y=CLASS_RANK, colour=CLASS)) + 
			geom_point() + geom_step() +
			geom_text(aes(label=PAIR_ID), size=2, nudge_x=-.05, nudge_y=.5) +
			scale_y_continuous(breaks=seq(0,50,5)) +
			scale_x_reverse(labels = scales::percent, breaks=seq(0,1,0.1)) +
			scale_colour_brewer(palette='Set1') +
			facet_grid(~CLASS) +
			labs(	x= '\nminimum proportion\n(proportion of ancestral windows out of all windows\nthat have reads from both individuals is at least x%)', 
					y='sequence pairs\n(#)\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-correctancestral.pdf',sep='')), w=12, h=7)
	#	write to file
	tmp			<- subset(rpa, CLASS=='ancestral in either direction\nor intermingled', select=c(COUPID, PAIR_ID, MALE_SANGER_ID, FEMALE_SANGER_ID, COUP_SC, CLASS, CLASS_PROP, CLASS_RANK))
	write.csv(tmp, row.names=FALSE, file=file.path(dir, paste(run,'-phsc-serodiscpairs-assignments.csv',sep='')) )
	#	numbers
	tmp			<- subset(rpa, CLASS=='ancestral in either direction\nor intermingled')
	cat('\nNumber of couples',tmp[, length(unique(COUPID))])
	setkey(tmp, COUPID, MALE_SANGER_ID, FEMALE_SANGER_ID)
	cat('\nNumber of sequence pairs from couples',nrow(unique(tmp)))
	setkey(tmp, RUN, COUPID, PAIR_ID)
	cat('\nNumber of pairings (including repeated sequence pairs)',nrow(unique(tmp)))	
	#
	#	plot on proportion of assignments in epidemiologically possible direction
	#
	rpb		<- subset(df, COUP_SC=='F->M')[, list(CLASS='prop ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='anc_fm'])/sum(WIN_OF_TYPE_P[TYPE=='anc_fm'|TYPE=='anc_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]	
	tmp		<- subset(df, COUP_SC=='M->F')[, list(CLASS='prop ancestral in correct direction', CLASS_PROP= sum(WIN_OF_TYPE_P[TYPE=='anc_mf'])/sum(WIN_OF_TYPE_P[TYPE=='anc_fm'|TYPE=='anc_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpb		<- rbind(rpb, tmp)
	tmp		<- rpb[order(CLASS, -CLASS_PROP),][, list(PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)), by='CLASS']
	rpb		<- merge(rpb, tmp, by=c('CLASS','PAIR_ID'))	
	setkey(rpb, CLASS, CLASS_RANK)
	ggplot(rpb, aes(x=CLASS_RANK, y=cumsum(CLASS_PROP))) + geom_line() + geom_point() +
			coord_cartesian(xlim=c(0, max(rpb[,CLASS_RANK])), ylim=c(0,max(rpb[,CLASS_RANK]))) +
			geom_abline(intercept=0, slope=0.5, colour='blue') +
			geom_abline(intercept=0, slope=1, colour='blue') +
			geom_text(aes(label=PAIR_ID), size=2, nudge_x=-.2, nudge_y=.8) +
			scale_y_continuous(expand=c(0,0)) +
			scale_x_continuous(expand=c(0,0)) +
			theme_bw() +
			labs(	x='\nsequence pairs\nof serodiscordant couples whose uninfected partners turns positive\n(cumulated)',
					y='# ancestral assignments in direction that is epidemiologically possible\nout of all ancestral assignments\n(cumulated)')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-only_possible_direction_assigned.pdf',sep='')), w=7, h=7)
	#
	#	plot on number of ancestral windows 
	#
	df[, WIN_OF_TYPE_N:= WIN_OF_TYPE_P*WIN_TOTAL]
	rpa		<- subset(df, COUP_SC=='F->M')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='anc_fm']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='anc_mf'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(df, COUP_SC=='M->F')[, list(IN_DIR= sum(WIN_OF_TYPE_N[TYPE=='anc_mf']), AGAINST_DIR=sum(WIN_OF_TYPE_N[TYPE=='anc_fm'])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= IN_DIR+AGAINST_DIR]
	rpa		<- melt(rpa, measure.vars=c('IN_DIR', 'AGAINST_DIR'))
	set(rpa, rpa[, which(variable=='IN_DIR')], 'variable', 'trm assignment in the only epidemiologically possible direction')
	set(rpa, rpa[, which(variable=='AGAINST_DIR')], 'variable', 'trm assignment against the only epidemiologically possible direction')
	setkey(rpa, PAIR_ID)	
	tmp		<- unique(rpa)[order(-WIN_TRM),][, list(COUPID=COUPID, PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(PAIR_ID, ' (', COUPID, ')', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PAIR_ID','COUPID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set1') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			#scale_x_discrete(labels=rpa[,PAIR_ID]) +
			labs(	x= '\nsequence pairs of couples', 
					y='number of ancestral windows\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_trm_windows.pdf',sep='')), w=10, h=7)
	#
	#	plot on all windows 
	#	
	rpa		<- subset(df, COUP_SC=='F->M')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='anc_fm'|TYPE=='anc_mf'|TYPE=='int']), NO_ANC=sum(WIN_OF_TYPE_N[!(TYPE=='anc_fm'|TYPE=='anc_mf'|TYPE=='int')])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	tmp		<- subset(df, COUP_SC=='M->F')[, list(ANY_ANC= sum(WIN_OF_TYPE_N[TYPE=='anc_mf'|TYPE=='anc_mf'|TYPE=='int']), NO_ANC=sum(WIN_OF_TYPE_N[!(TYPE=='anc_fm'|TYPE=='anc_mf'|TYPE=='int')])), by=c('PAIR_ID','COUPID','MALE_SANGER_ID','FEMALE_SANGER_ID','COUP_SC')]
	rpa		<- rbind(rpa, tmp, use.names=TRUE)
	rpa[, WIN_TRM:= ANY_ANC]
	rpa		<- melt(rpa, measure.vars=c('ANY_ANC', 'NO_ANC'))
	set(rpa, rpa[, which(variable=='ANY_ANC')], 'variable', 'ancestral assignments or intermingled')
	set(rpa, rpa[, which(variable=='NO_ANC')], 'variable', 'sibling or disconnected')
	setkey(rpa, PAIR_ID)	
	tmp		<- unique(rpa)[order(-WIN_TRM),][, list(COUPID=COUPID, PAIR_ID=PAIR_ID, CLASS_RANK=seq_along(PAIR_ID)) ]
	set(tmp, NULL, 'CLASS_RANK', tmp[, factor(CLASS_RANK, levels=CLASS_RANK, labels=paste(PAIR_ID, ' (', COUPID, ')', sep=''))])
	rpa		<- merge(rpa, tmp, by=c('PAIR_ID','COUPID'))	
	setkey(rpa, variable, CLASS_RANK)	
	ggplot(rpa, aes(x=CLASS_RANK, y=value, fill=variable)) + 
			geom_bar(stat='identity',position='stack') +
			scale_fill_brewer(palette='Set2') +
			theme_bw() + theme(legend.position='bottom',axis.text.x = element_text(angle = 90, hjust = 1)) +
			#scale_x_discrete(labels=rpa[,PAIR_ID]) +
			labs(	x= '\nsequence pairs of couples', 
					y='number of ancestral windows\n',
					colour='phyloscanner\ntransmission assignments',
					title='\nphyloscanner transmission assignments\nto RCCS serodiscordant couples\nin which the uninfected partner is infected during follow-up\n')
	ggsave(file=file.path(dir, paste(run,'-phsc-serodiscpairs-number_other_windows.pdf',sep='')), w=10, h=7)
	
	
	
	
	#
	#	inconsistent pairings of same sequences
	#	inconsistent pairings of same couples
	
	#}
	#
	#	re-examine phylogenies for all sero-discordant couples
	#	
	tmp		<- subset(df, COUP_SC=='F->M' | COUP_SC=='M->F')
	setkey(tmp, RUN, PAIR_ID, MALE_SANGER_ID, FEMALE_SANGER_ID)
	dir		<- tmp[1,DIR]
	cat('\ndir is',dir,'\trun is',run)			
	setkey(tmp, MALE_SANGER_ID, FEMALE_SANGER_ID, PAIR_ID)				
	rpoku	<- unique(tmp)
	load( tmp[1, FILE] )	#loads phs dtrms dtrees
	invisible(sapply(seq_len(nrow(rpoku)), function(ii)
					{								
						pair.id		<- rpoku[ii, PAIR_ID]
						pty.run		<- rpoku[ii, PTY_RUN]
						dfs			<- subset(dtrees, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
						dfs[, TITLE:= dfs[, paste('pair', pair.id, '\n', rpoku[ii, COUP_SC], '\nid M: ', rpoku[ii, MALE_RID], ' (', rpoku[ii, MALE_SANGER_ID], ')\nid F: ', rpoku[ii, FEMALE_RID], ' (', rpoku[ii, FEMALE_SANGER_ID], ')\nrun ', pty.run, '\nwindow ', W_FROM,'-', W_TO,sep='')]]			
						plot.file	<- file.path(dir, paste(run,'-phsc-serodiscpairs-',rpoku[ii, COUP_SC],'-M-', rpoku[ii, MALE_RID],'-F-',rpoku[ii, FEMALE_RID],'-', pair.id,'.pdf',sep=''))			
						invisible(phsc.plot.selected.pairs(phs, dfs, rpoku[ii, MALE_SANGER_ID], rpoku[ii, FEMALE_SANGER_ID], plot.file=plot.file, pdf.h=150, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40))
					}))				
	
}

RakaiCouples.save.couples.to.rda<- function()
{
	require(data.table)
	wdir				<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples"
	#	get epi info
	tmp		<- RakaiCirc.epi.get.info()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#	get sequence info and add to recipients
	load('~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision/RCCS_SeqInfo_160816.rda')		
	rs		<- subset(rs, !is.na(VISIT))	
	
	#	load combination of all couples sequences
	infile	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_phscruns.rda'
	#	load Kates couple data
	infile	<- '~/Dropbox (Infectious Disease)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/Pangea_Couples.csv'
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
	#	486 couples
	#
	#	get all possible unique pairings of SANGER_IDs from couples
	#	347
	rp		<- subset(rc, !is.na(MALE_SANGER_ID), c(COUPID, MALE_RID, MALE_SANGER_ID, MALE_TAXA ))
	tmp		<- subset(rc, !is.na(FEMALE_SANGER_ID), c(COUPID, FEMALE_RID, FEMALE_SANGER_ID, FEMALE_TAXA ))	
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
	#	get info on sequences: date of sampling
	#	for males
	tmp		<- unique(subset(rs, !is.na(SID), select=c(SID, DATE)), by='SID')
	setnames(tmp, c('SID','DATE'), c('MALE_SANGER_ID','MALE_SEQDATE'))
	rp		<- merge(rp, tmp, by='MALE_SANGER_ID', all.x=1)
	tmp		<- unique(subset(rs, !is.na(SID), select=c(SID, DATE)), by='SID')
	setnames(tmp, c('SID','DATE'), c('FEMALE_SANGER_ID','FEMALE_SEQDATE'))	
	rp		<- merge(rp, tmp, by='FEMALE_SANGER_ID', all.x=1)				
	#
	#	define direction based on seroconversion dates
	#
	rp[, COUP_SC:='seropos']
	set(rp, rp[, which(MALE_LASTNEGDATE>=FEMALE_FIRSTPOSDATE)],'COUP_SC','F->M')
	set(rp, rp[, which(FEMALE_LASTNEGDATE>=MALE_FIRSTPOSDATE)],'COUP_SC','M->F')
	set(rp, rp[, which(FEMALE_FIRSTPOSDATE==MALE_FIRSTPOSDATE)],'COUP_SC','seroinc')
	#	F->M    M->F seroinc seropos 
	#	20      31      69     417 

	#	link individual household data from most recent visit
	set(rp, rp[, which(PAIR_TYPE=='no joint household data' & FEMALE_HH_NUM==MALE_HH_NUM)], 'PAIR_TYPE', 'stable cohabiting')
	set(rp, rp[, which(PAIR_TYPE=='no joint household data' & FEMALE_HH_NUM!=MALE_HH_NUM)], 'PAIR_TYPE', 'not always cohabiting')
	
	#	save
	save(rp, file="~/Dropbox (Infectious Disease)/Rakai Fish Analysis/couples/Couples_PANGEA_HIV_n4562_Imperial_v151113_info.rda")	
}

RakaiCirc.circ.dev160907<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	wdir				<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
	
	#	get epi info
	tmp		<- RakaiCirc.epi.get.info()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#	get sequence info and add to recipients
	load(file.path(wdir,'RCCS_SeqInfo_160816.rda'))		
	rs		<- subset(rs, !is.na(VISIT))	
	#	load female recipients
	rrec	<- RakaiCirc.recipient.female.get.info()
	#	add sequence info to recipients
	rrec	<- merge(rrec, unique(subset(rs, select=c(RID, PID, SID))), by='RID', all.x=1)
	#	select female recipients with at least one PANGEA sequence
	rrec	<- subset(rrec, !is.na(SID))
	
	
	#	select likely pairs, use same selection as before just for consistency
	select.discsib	<- 0.65	
	#
	#	collect runs
	#
	infiles	<- data.table(	F_TRM= c(	'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w200/RCCS_160902_w200_trmStats.rda',
										'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w220/RCCS_160902_w220_trmStats.rda',
										'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w250/RCCS_160902_w250_trmStats.rda',
										'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w270/RCCS_160902_w270_trmStats.rda',
										'~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w280/RCCS_160902_w280_trmStats.rda'
							))
	infiles[, F_PH:= gsub('trmStats.rda','trees.rda', F_TRM)]
	infiles[, DIR:= dirname(F_TRM)]
	infiles[, RUN:= gsub('_trmStats.rda','',basename(F_TRM))]				
	#
	#	for each run: get list of pairs
	#	
	rp		<- infiles[, {
				#F_TRM	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w270/RCCS_160902_w270_trmStats.rda'
				#F_PH	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/Rakai_ptoutput_160902_w270/RCCS_160902_w270_trees.rda'
				load(F_TRM)	#loads df
				dlkl	<- copy(df)
				load(F_PH)	#loads phs and dfr
				#	select likely pairs
				dlkl[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]
				tmp		<- dcast.data.table(dlkl, PAIR_ID~TYPE, value.var='WIN_OF_TYPE_P')
				set(tmp, tmp[, which(is.na(disconnected))], 'disconnected', 0)
				set(tmp, tmp[, which(is.na(sib))], 'sib', 0)
				tmp		<- subset(tmp, disconnected+sib < select.discsib, PAIR_ID)
				dlkl	<- merge(dlkl, tmp, by='PAIR_ID')
				#	get likely pairs that involve female recipients
				tmp		<- copy(dlkl)
				setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
				set(tmp, tmp[, which(TYPE=='anc_12')], 'TYPE', 'anc')
				set(tmp, tmp[, which(TYPE=='anc_21')], 'TYPE', 'anc_12')
				set(tmp, tmp[, which(TYPE=='anc')], 'TYPE', 'anc_21')
				tmp		<- rbind(tmp, dlkl)
				setnames(tmp, 'ID1', 'SID')
				rp		<- merge(tmp, rrec, by='SID')
				setnames(rp, 'ID2', 'TR_SID')
				cat('\nFound female recipients among likely pairs, n=', rp[, length(unique(RID))])
				#	get info on transmitters: RID and PID
				tmp		<- subset(rs, !is.na(SID), select=c(RID, PID, SID))
				setkey(tmp, SID)
				tmp		<- unique(tmp)
				setnames(tmp, colnames(tmp), paste('TR_',colnames(tmp),sep=''))	
				#
				#	TODO check whereabouts of 12559_1_9 15065_1_32 15034_1_79 15081_1_71 15430_1_73
				#
				tmp2	<- setdiff( rp[, unique(TR_SID)], tmp[, TR_SID] )
				if(length(tmp2))
					cat('\nWarning: Sanger ID in rp that is not in dictionary', tmp2)
				rp		<- merge(rp, tmp, by='TR_SID')
				#	get info on transmitters: demographic stuff	
				tmp		<- subset(rd, !is.na(PID), select=c(RID, SEX, REGION, COMM_NUM, HH_NUM, BIRTHDATE, LASTNEGVIS, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, RELIGION))
				setkey(tmp, RID)
				tmp		<- unique(tmp)
				setnames(tmp, colnames(tmp), paste('TR_',colnames(tmp),sep=''))
				rp		<- merge(rp,tmp,by='TR_RID')
				#	get info on sequences: date of sampling
				tmp		<- unique(subset(rs, !is.na(SID), select=c(SID, DATE)))
				setnames(tmp, 'DATE', 'SEQ_DATE')
				rp		<- merge(rp, tmp, by='SID')
				setnames(tmp, c('SID','SEQ_DATE'), c('TR_SID','TR_SEQ_DATE'))
				rp		<- merge(rp, tmp, by='TR_SID')				
				
				rp
			}, by=c('RUN','DIR','F_TRM','F_PH')]
	#
	#	for each run: plot pairs & trees	
	#
	setkey(rp, RUN, PAIR_ID, SID, TR_SID)
	for( run in rp[, unique(RUN)] )
	{		
		#run		<- 'RCCS_160902_w270'
		dir		<- subset(rp, RUN==run)[1,DIR]
		cat('\ndir is',dir,'\trun is',run)
		df		<- subset(rp, RUN==run)
		rpok	<- subset(df, TR_SEX=='M')	
		rpff	<- subset(df, TR_SEX=='F')
		rpoku	<- unique(rpok)
		rpffu	<- unique(rpff)
		
		#	for every F-F pair, keep only one instance
		#	and select different individuals
		rpffu		<- merge(rpffu, rpffu[, list(SID=SID[1],TR_SID=TR_SID[1]), by='PAIR_ID'], by=c('PAIR_ID','SID','TR_SID'))	
		rpffu		<- subset(rpffu, RID!=TR_RID)		
		#
		#	plot evidence
		#		
		tmp		<- copy(rpoku)
		tmp		<- tmp[order(-PAIR_ID),]
		tmp[, LABEL:= factor(PAIR_ID, levels=PAIR_ID, labels=paste('Pair',PAIR_ID, 	'\nPerson1 ', RID, ' ', SID,' -sex:',SEX,' -loc:',REGION,',',COMM_NUM,',',HH_NUM,' -birth:',BIRTHDATE,' -neg:',LASTNEGDATE,' -pos:',FIRSTPOSDATE,' -seq:',SEQ_DATE,
								'\n<->', 
								'\nPerson2 ', TR_RID, ' ', TR_SID,' -sex:',TR_SEX,' -loc:',TR_REGION,',',TR_COMM_NUM,',',TR_HH_NUM,' -birth:',TR_BIRTHDATE,' -neg:',TR_LASTNEGDATE,' -pos:',TR_FIRSTPOSDATE,' -seq:',TR_SEQ_DATE,																				
								'\n',sep=''))]
		tmp		<- merge(subset(tmp, select=c(PAIR_ID, LABEL)), df, by='PAIR_ID')
		set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('anc_12','anc_21','sib','int','disconnected'), labels=c('1 ancestral to 2','2 ancestral to 1','1, 2 are siblings','1, 2 are intermingled','1, 2 are disconnected'))])
		ggplot(tmp, aes(x=LABEL, y=WIN_OF_TYPE, fill=TYPE)) +
				geom_bar(stat='identity', position='stack') +
				coord_flip() +
				labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
				scale_fill_manual(values=c('1 ancestral to 2'="#9E0142",'2 ancestral to 1'="#F46D43",'1, 2 are siblings'="#ABDDA4",'1, 2 are intermingled'="#3288BD",'1, 2 are disconnected'='grey50')) +
				theme_bw() + theme(legend.position='top') +
				guides(fill=guide_legend(ncol=2))
		ggsave(file=file.path(dir, paste(run,'_RCCS_lkltransmissionpairs_withfemaleseroconverter_discsiblt',100*select.discsib,'.pdf',sep='')), w=15, h=max(3,0.23*nrow(tmp)), limitsize = FALSE)	
		#	look at F-F pairs
		tmp			<- copy(rpffu)
		tmp			<- tmp[order(-PAIR_ID),]
		tmp[, LABEL:= factor(PAIR_ID, levels=PAIR_ID, labels=paste('Pair',PAIR_ID, 	'\nPerson1 ', RID, ' ', SID,' -sex:',SEX,' -loc:',REGION,',',COMM_NUM,',',HH_NUM,' -birth:',BIRTHDATE,' -neg:',LASTNEGDATE,' -pos:',FIRSTPOSDATE,' -seq:',SEQ_DATE,
								'\n<->', 
								'\nPerson2 ', TR_RID, ' ', TR_SID,' -sex:',TR_SEX,' -loc:',TR_REGION,',',TR_COMM_NUM,',',TR_HH_NUM,' -birth:',TR_BIRTHDATE,' -neg:',TR_LASTNEGDATE,' -pos:',TR_FIRSTPOSDATE,' -seq:',TR_SEQ_DATE,																				
								'\n',sep=''))]	
		tmp		<- merge(subset(tmp, select=c(PAIR_ID, LABEL)), df, by='PAIR_ID')
		set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('anc_12','anc_21','sib','int','disconnected'), labels=c('1 ancestral to 2','2 ancestral to 1','1, 2 are siblings','1, 2 are intermingled','1, 2 are disconnected'))])
		ggplot(tmp, aes(x=LABEL, y=WIN_OF_TYPE, fill=TYPE)) +
				geom_bar(stat='identity', position='stack') +
				coord_flip() +
				labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
				scale_fill_manual(values=c('1 ancestral to 2'="#9E0142",'2 ancestral to 1'="#F46D43",'1, 2 are siblings'="#ABDDA4",'1, 2 are intermingled'="#3288BD",'1, 2 are disconnected'='grey50')) +
				theme_bw() + theme(legend.position='top') +
				guides(fill=guide_legend(ncol=2))
		ggsave(file=file.path(dir, paste(run,'_RCCS_lkltransmissionpairs_femalefemale_discsiblt',100*select.discsib,'.pdf',sep='')), w=15, h=max(3,0.2*nrow(tmp)), limitsize = FALSE)				
	}
	#
	#	re-examine phylogenies
	#
	setkey(rp, RUN, PAIR_ID, SID, TR_SID)
	for( run in rp[, unique(RUN)] )
	{
		dir		<- subset(rp, RUN==run)[1,DIR]
		cat('\ndir is',dir,'\trun is',run)
		df		<- subset(rp, RUN==run)
		rpok	<- subset(df, TR_SEX=='M')	
		rpff	<- subset(df, TR_SEX=='F')
		rpoku	<- unique(rpok)
		rpffu	<- unique(rpff)
		rpffu	<- merge(rpffu, rpffu[, list(SID=SID[1],TR_SID=TR_SID[1]), by='PAIR_ID'], by=c('PAIR_ID','SID','TR_SID'))	
		rpffu	<- subset(rpffu, RID!=TR_RID)		
		
		load( df[1, F_PH] )	#loads phs and dfr
		setkey(rpok, PAIR_ID, SID, TR_SID)
		rpoku	<- unique(rpok)
		invisible(sapply(seq_len(nrow(rpoku)), function(ii)
						{								
							pair.id		<- rpoku[ii, PAIR_ID]
							pty.run		<- rpoku[ii, PTY_RUN]
							dfs			<- subset(dfr, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
							dfs[, TITLE:= dfs[, paste('pair', pair.id, 'ids', rpoku[ii, TR_SID], rpoku[ii, SID], '\nrun', pty.run, '\nwindow', W_FROM,'-', W_TO)]]			
							plot.file	<- file.path(dir, paste(run,'_RCCS_lkltransmissionpairs_withfemaleseroconverter_discsiblt65_pair',pair.id,'.pdf',sep=''))			
							invisible(phsc.plot.selected.pairs(phs, dfs, rpoku[ii, TR_SID], rpoku[ii, SID], plot.file=plot.file, pdf.h=50, pdf.rw=10))
						}))		
		tmp		<- copy(rpffu)
		invisible(sapply(seq_len(nrow(tmp)), function(ii)
						{
							pair.id		<- tmp[ii, PAIR_ID]
							pty.run		<- tmp[ii, PTY_RUN]
							dfs			<- subset(dfr, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
							dfs[, TITLE:= dfs[, paste('pair', pair.id, 'ids', tmp[ii, TR_SID], tmp[ii, SID], '\nrun', pty.run, '\nwindow', W_FROM,'-', W_TO)]]			
							plot.file	<- file.path(dir, paste(run,'_RCCS_lkltransmissionpairs_femalefemale_discsiblt65_pair',pair.id,'.pdf',sep=''))			
							invisible(phsc.plot.selected.pairs(phs, dfs, tmp[ii, TR_SID], tmp[ii, SID], plot.file=plot.file, pdf.h=50, pdf.rw=10))
						}))
	}

	
	#	the 200 + 270 phylogenies look identical!?!?!
	subset(rp, grepl('w270',RUN) & PAIR_ID==3)
	subset(dfr, grepl('w270',RUN), select=c(PTY_RUN, W_FROM, W_TO, IDX))
	
	subset(rp, grepl('w200',RUN) & PAIR_ID==5)
	
}	

RakaiCirc.circ.dev160901<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	wdir				<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
	
	#	get epi info
	tmp		<- RakaiCirc.epi.get.info()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#	get sequence info and add to recipients
	load(file.path(wdir,'RCCS_SeqInfo_160816.rda'))		
	rs		<- subset(rs, !is.na(VISIT))	
	
	
	#	load likely transmissions summary from phyloscanner
	infile.phsc.trms	<- "~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/pty_Rakai_160825/RCCS_run160825_lkltrms_summary.rda"
	load(infile.phsc.trms)
	dlkl				<- copy(df)
	#	load trees from phyloscanner
	infile.phsc.trees	<- '~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/pty_Rakai_160825/RCCS_run160825_all_trees.rda'	
	load(infile.phsc.trees)
	phs					<- tmp$phs
	stat.phs			<- tmp$dfr
	
	
	
	#	select likely pairs
	select.discsib	<- 0.65
	dlkl[, WIN_OF_TYPE_P:=WIN_OF_TYPE/WIN_TOTAL]
	tmp		<- dcast.data.table(dlkl, PAIR_ID~TYPE, value.var='WIN_OF_TYPE_P')
	set(tmp, tmp[, which(is.na(disconnected))], 'disconnected', 0)
	set(tmp, tmp[, which(is.na(sib))], 'sib', 0)
	tmp		<- subset(tmp, disconnected+sib < select.discsib, PAIR_ID)
	dlkl	<- merge(dlkl, tmp, by='PAIR_ID')
	#	load female recipients
	rrec	<- RakaiCirc.recipient.female.get.info()
	#	add sequence info to recipients
	rrec	<- merge(rrec, unique(subset(rs, select=c(RID, PID, SID))), by='RID', all.x=1)
	#	select female recipients with at least one PANGEA sequence
	rrec	<- subset(rrec, !is.na(SID))
	#	get likely pairs that involve female recipients
	tmp		<- copy(dlkl)
	setnames(tmp, c('ID1','ID2'), c('ID2','ID1'))
	set(tmp, tmp[, which(TYPE=='anc_12')], 'TYPE', 'anc')
	set(tmp, tmp[, which(TYPE=='anc_21')], 'TYPE', 'anc_12')
	set(tmp, tmp[, which(TYPE=='anc')], 'TYPE', 'anc_21')
	tmp		<- rbind(tmp, dlkl)
	setnames(tmp, 'ID1', 'SID')
	rp		<- merge(tmp, rrec, by='SID')
	setnames(rp, 'ID2', 'TR_SID')
	cat('\nFound female recipients among likely pairs, n=', rp[, length(unique(RID))])
	#	get info on transmitters: RID and PID
	tmp		<- subset(rs, !is.na(SID), select=c(RID, PID, SID))
	setkey(tmp, SID)
	tmp		<- unique(tmp)
	setnames(tmp, colnames(tmp), paste('TR_',colnames(tmp),sep=''))	
	#
	#	TODO check whereabouts of 12559_1_9 15065_1_32 15034_1_79 15081_1_71 15430_1_73
	#
	tmp2	<- setdiff( rp[, unique(TR_SID)], tmp[, TR_SID] )
	if(length(tmp2))
		cat('\nWarning: Sanger ID in rp that is not in dictionary', tmp2)
	rp		<- merge(rp, tmp, by='TR_SID')
	#	get info on transmitters: demographic stuff	
	tmp		<- subset(rd, !is.na(PID), select=c(RID, SEX, REGION, COMM_NUM, HH_NUM, BIRTHDATE, LASTNEGVIS, LASTNEGDATE, FIRSTPOSVIS, FIRSTPOSDATE, RELIGION))
	setkey(tmp, RID)
	tmp		<- unique(tmp)
	setnames(tmp, colnames(tmp), paste('TR_',colnames(tmp),sep=''))
	rp		<- merge(rp,tmp,by='TR_RID')
	#	get info on sequences: date of sampling
	tmp		<- unique(subset(rs, !is.na(SID), select=c(SID, DATE)))
	setnames(tmp, 'DATE', 'SEQ_DATE')
	rp		<- merge(rp, tmp, by='SID')
	setnames(tmp, c('SID','SEQ_DATE'), c('TR_SID','TR_SEQ_DATE'))
	rp		<- merge(rp, tmp, by='TR_SID')
	setkey(rp, SID, TR_SID)
	rpok	<- subset(rp, TR_SEX=='M')	
	rpff	<- subset(rp, TR_SEX=='F')
	rpoku	<- unique(rpok)
	rpffu	<- unique(rpff)
	#	for every F-F pair, keep only one instance
	#	and select different individuals
	rpffu		<- merge(rpffu, rpffu[, list(SID=SID[1],TR_SID=TR_SID[1]), by='PAIR_ID'], by=c('PAIR_ID','SID','TR_SID'))	
	rpffu		<- subset(rpffu, RID!=TR_RID)		
	#
	#	plot evidence
	#		
	tmp		<- copy(rpoku)
	tmp		<- tmp[order(-PAIR_ID),]
	tmp[, LABEL:= factor(PAIR_ID, levels=PAIR_ID, labels=paste('Pair',PAIR_ID, 	'\nPerson1 ', RID, ' ', SID,' -sex:',SEX,' -loc:',REGION,',',COMM_NUM,',',HH_NUM,' -birth:',BIRTHDATE,' -neg:',LASTNEGDATE,' -pos:',FIRSTPOSDATE,' -seq:',SEQ_DATE,
							'\n<->', 
							'\nPerson2 ', TR_RID, ' ', TR_SID,' -sex:',TR_SEX,' -loc:',TR_REGION,',',TR_COMM_NUM,',',TR_HH_NUM,' -birth:',TR_BIRTHDATE,' -neg:',TR_LASTNEGDATE,' -pos:',TR_FIRSTPOSDATE,' -seq:',TR_SEQ_DATE,																				
							'\n',sep=''))]
	tmp		<- merge(subset(tmp, select=c(PAIR_ID, LABEL)), df, by='PAIR_ID')
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('anc_12','anc_21','sib','int','disconnected'), labels=c('1 ancestral to 2','2 ancestral to 1','1, 2 are siblings','1, 2 are intermingled','1, 2 are disconnected'))])
	ggplot(tmp, aes(x=LABEL, y=WIN_OF_TYPE, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('1 ancestral to 2'="#9E0142",'2 ancestral to 1'="#F46D43",'1, 2 are siblings'="#ABDDA4",'1, 2 are intermingled'="#3288BD",'1, 2 are disconnected'='grey50')) +
			theme_bw() + theme(legend.position='top') +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(wdir, paste('150901_RCCS_lkltransmissionpairs_withfemaleseroconverter_discsiblt',100*select.discsib,'.pdf',sep='')), w=15, h=0.23*nrow(tmp), limitsize = FALSE)	
	#	look at F-F pairs
	tmp			<- copy(rpffu)
	tmp			<- tmp[order(-PAIR_ID),]
	tmp[, LABEL:= factor(PAIR_ID, levels=PAIR_ID, labels=paste('Pair',PAIR_ID, 	'\nPerson1 ', RID, ' ', SID,' -sex:',SEX,' -loc:',REGION,',',COMM_NUM,',',HH_NUM,' -birth:',BIRTHDATE,' -neg:',LASTNEGDATE,' -pos:',FIRSTPOSDATE,' -seq:',SEQ_DATE,
							'\n<->', 
							'\nPerson2 ', TR_RID, ' ', TR_SID,' -sex:',TR_SEX,' -loc:',TR_REGION,',',TR_COMM_NUM,',',TR_HH_NUM,' -birth:',TR_BIRTHDATE,' -neg:',TR_LASTNEGDATE,' -pos:',TR_FIRSTPOSDATE,' -seq:',TR_SEQ_DATE,																				
							'\n',sep=''))]	
	tmp		<- merge(subset(tmp, select=c(PAIR_ID, LABEL)), df, by='PAIR_ID')
	set(tmp, NULL, 'TYPE', tmp[, factor(TYPE, levels=c('anc_12','anc_21','sib','int','disconnected'), labels=c('1 ancestral to 2','2 ancestral to 1','1, 2 are siblings','1, 2 are intermingled','1, 2 are disconnected'))])
	ggplot(tmp, aes(x=LABEL, y=WIN_OF_TYPE, fill=TYPE)) +
			geom_bar(stat='identity', position='stack') +
			coord_flip() +
			labs(x='', y='number of read windows', fill='topology of clades\nbetween patient pairs') +
			scale_fill_manual(values=c('1 ancestral to 2'="#9E0142",'2 ancestral to 1'="#F46D43",'1, 2 are siblings'="#ABDDA4",'1, 2 are intermingled'="#3288BD",'1, 2 are disconnected'='grey50')) +
			theme_bw() + theme(legend.position='top') +
			guides(fill=guide_legend(ncol=2))
	ggsave(file=file.path(wdir, paste('150901_RCCS_lkltransmissionpairs_femalefemale_discsiblt',100*select.discsib,'.pdf',sep='')), w=15, h=0.2*nrow(tmp), limitsize = FALSE)	
	#
	#	re-examine phylogenies
	#	
	setkey(rpok, PAIR_ID, SID, TR_SID)
	rpoku	<- unique(rpok)
	invisible(sapply(seq_len(nrow(rpoku)), function(ii)
					{
						pair.id		<- rpoku[ii, PAIR_ID]
						pty.run		<- rpoku[ii, PTY_RUN]
						dfs			<- subset(dfr, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
						dfs[, TITLE:= dfs[, paste('pair', pair.id, 'ids', rpoku[ii, TR_SID], rpoku[ii, SID], '\nrun', pty.run, '\nwindow', W_FROM,'-', W_TO)]]			
						plot.file	<- file.path(wdir, paste('150901_RCCS_lkltransmissionpairs_withfemaleseroconverter_discsiblt65_pair',pair.id,'.pdf',sep=''))			
						invisible(phsc.plot.selected.pairs(phs, dfs, rpoku[ii, TR_SID], rpoku[ii, SID], plot.file=plot.file, pdf.h=50, pdf.rw=10))
					}))
	
	tmp		<- copy(rpffu)
	invisible(sapply(seq_len(nrow(tmp)), function(ii)
					{
						pair.id		<- tmp[ii, PAIR_ID]
						pty.run		<- tmp[ii, PTY_RUN]
						dfs			<- subset(dfr, PTY_RUN==pty.run, select=c(PTY_RUN, W_FROM, W_TO, IDX))
						dfs[, TITLE:= dfs[, paste('pair', pair.id, 'ids', tmp[ii, TR_SID], tmp[ii, SID], '\nrun', pty.run, '\nwindow', W_FROM,'-', W_TO)]]			
						plot.file	<- file.path(wdir, paste('150901_RCCS_lkltransmissionpairs_femalefemale_discsiblt65_pair',pair.id,'.pdf',sep=''))			
						invisible(phsc.plot.selected.pairs(phs, dfs, tmp[ii, TR_SID], tmp[ii, SID], plot.file=plot.file, pdf.h=50, pdf.rw=10))
					}))
}	
	
RakaiCirc.circ.dev160815<- function()
{
	require(data.table)
	require(scales)
	require(ggplot2)
	wdir				<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"	
	tmp		<- RakaiCirc.epi.get.info()
	rh		<- tmp$rh
	rd		<- tmp$rd
	#
	#	make individual timeline over visits
	#
	rt	<- RakaiCirc.circ.timelines.init.160816(rh, rd)
	#
	#	add first sequences
	#	load sequence info. expect "rs".
	#	(if not exist, run RakaiCirc.seq.get.info() )
	#
	load(file.path(wdir,'RCCS_SeqInfo_160816.rda'))	
	#	TODO: has not visit: E21593M2
	rs		<- subset(rs, !is.na(VISIT))
	rt		<- RakaiCirc.circ.timelines.addfirstseq.160816(rt, rs)
	setkey(rt, RID, ROUND)
	#RakaiCirc.circ.timelines.plots(rt, wdir)
	
	rrec	<- RakaiCirc.recipient.female.get.info(wdir=NA)		
	#
	#	load info from p24 tree, PANGEA tree, genetic distances
	#
	if(0)
	{
		tmp		<- RakaiCirc.seq.get.phylogenies()
		phf		<- tmp$phf
		php		<- tmp$php
		phfi	<- tmp$phfi
		phpi	<- tmp$phpi		
	}
	load(file=file.path(wdir,'RCCS_PhInfo_160825.rda'))
	#
	#	load info on phylotype runs
	#
	load("~/Dropbox (Infectious Disease)/2015_PANGEA_DualPairsFromFastQIVA/data/PANGEA_HIV_n5003_Imperial_v160110_UG_gag_coinfinput_160219.rda")	
	ptyi	<- subset(pty.runs, select=c(TAXA, FILE_ID, PTY_RUN))
	set(ptyi,NULL,'TAXA', ptyi[,gsub('_','-',gsub('_S[0-9]+','',TAXA))])
	setnames(ptyi, c('TAXA','FILE_ID'), c('PID','SID'))
	#
	#	check how many rec in tree, how many prob transmitters, and in which sequ sample
	#
	#pd.max	<- 0.02
	pd.max	<- 0.04
	#	may have multiple seqs per female recipient, keep all
	rrecs	<- merge(rrec, unique(subset(rs, select=c(RID, PID, SEQID, SEQIDb))), by='RID', all.x=1)
	rrecs	<- merge(rrecs, ptyi, by='PID', all.x=1)
	#	add taxa in p24 tree that are close
	#	phfi is only lower triangular, complete 
	tmp		<- subset(phfi, PD<=pd.max)	
	tmp2	<- copy(tmp)	
	setnames(tmp2, c('SEQIDb','SEQIDb2'), c('SEQIDb2','SEQIDb'))	
	tmp		<- rbind(tmp, tmp2, use.names=TRUE)
	set(tmp, NULL, 'SEQIDb', tmp[, as.character(SEQIDb)])
	set(tmp, NULL, 'SEQIDb2', tmp[, as.character(SEQIDb2)])
	#	now add
	rrecs	<- merge(rrecs, tmp, all=1, by='SEQIDb')
	rrecs	<- subset(rrecs, !is.na(RID))
	#	add basic epi info and IDs for transmitters
	setnames(rrecs, 'SEQIDb2', 'TR_SEQIDb') 
	tmp		<- unique(subset(rs, !is.na(SEQIDb), select=c(RID, PID, SEQIDb, SEQTYPE)))
	tmp2	<- subset(rd, select=c(RID, SEX))
	setkey(tmp2, RID)
	tmp		<- merge(tmp, unique(tmp2), by='RID')
	#	add basic info about phylotype runs for transmitters	
	tmp		<- merge(tmp, ptyi, by='PID', all.x=1)
	setnames(tmp, c('RID', 'PID', 'SEQIDb', 'SID', 'PTY_RUN','SEQTYPE', 'SEX'), c('TR_RID', 'TR_PID', 'TR_SEQIDb', 'TR_SID', 'TR_PTY_RUN', 'TR_SEQTYPE', 'TR_SEX'))	
	rp		<- merge(rrecs, tmp, by='TR_SEQIDb', all.x=1)	
	tmp		<- subset(rp, PTY_RUN!=TR_PTY_RUN & TR_SEX=='M')	
	if(nrow(tmp))
	{
		cat('\nWarning: potential M->F transmission pair not in phylotype groupings\n', paste( tmp[, paste(TR_PID, '->', PID, sep='')], collapse='\n'))
		save(tmp, file=file.path(wdir, paste('150830_FastTree_gagAllSmall_pairs_patristic',pd.max*100,'_missedpotentialpairs.rda',sep='')))
	}
	#	plot tree with pairs highlighted
	if(0)
	{
		tmp		<- data.table(SEQIDb=phf$tip.label, IDX=seq_len(Ntip(phf)), TIP_COL='black')
		tmp2	<- data.table(SEQIDb=subset(rp, !is.na(TR_SEQIDb))[, unique(na.omit(c(SEQIDb, TR_SEQIDb))) ], PAIR='Y')	
		tmp		<- merge(tmp, tmp2, by='SEQIDb', all.x=1)
		set(tmp, tmp[, which(PAIR=='Y')], 'TIP_COL', 'red')
		setkey(tmp, IDX)
		file	<- file.path(wdir, paste('150825_FastTree_gagAllSmall_pairs_patristic',pd.max*100,'pc.pdf',sep=''))
		invisible(hivc.clu.plot(phf, tip.color=tmp[,TIP_COL], file=file, pdf.scaley=100, show.tip.label=TRUE, pdf.width=30))		
	}
	#
	#	get sample counts
	#
	#	for every female recipient, count prob transmitters by SEQ_TYPE
	tmp		<- subset(rp, SC_WINDOW<=2)
	rpi		<- tmp[, list(	FIRSTPOSDATE=			FIRSTPOSDATE[1],
							SEQ_R_N= 				ifelse(any(!is.na(SEQID)), length(unique(SEQID)), 0L),
							SEQ_R_PANGEA= 			any(grepl('PANGEA', SEQ_TYPE)),
							PHP24_R_N= 				ifelse(any(!is.na(SEQIDb)), length(unique(SEQIDb)), 0L),
							PHP24_TR_N= 			ifelse(any(!is.na(TR_SEQIDb)), length(unique(TR_SEQIDb)), 0L),
							PHP24_MTR_N= 			ifelse(any(!is.na(TR_SEQIDb[TR_SEX=='M'])), length(unique(TR_SEQIDb[TR_SEX=='M'])), 0L),
							PHP24_MTR_MINPD=		ifelse(any(!is.na(PD[TR_SEX=='M'])), min(PD[TR_SEX=='M'], na.rm=TRUE), NA_real_),
							PHP24_MTR_PANGEA= 		length(which(grepl('PANGEA', TR_SEQTYPE[TR_SEX=='M']))),
							PHP24_MTR_PANGEA_MINPD=	ifelse(any(!is.na(PD[TR_SEX=='M' & TR_SEQTYPE=='PANGEA'])), min(PD[TR_SEX=='M' & TR_SEQTYPE=='PANGEA'], na.rm=TRUE), NA_real_)
							), by='RID']
	#	
	set(rpi, NULL, 'FIRSTPOSDATEc', rpi[, cut(FIRSTPOSDATE, breaks=c(1994,1999, 2011, 2016), labels=c('<1999','<2011','2011-2015'))])
	rpi[, REC_S:= as.integer(SEQ_R_N>0)]
	set(rpi, rpi[, which(PHP24_R_N>0)], 'REC_S', 2L)
	set(rpi, rpi[, which(SEQ_R_PANGEA & REC_S==1L)], 'REC_S', 3L)
	set(rpi, rpi[, which(SEQ_R_PANGEA & REC_S==2L)], 'REC_S', 4L)	
	set(rpi, NULL, 'REC_S', rpi[, factor(REC_S, levels=c(0L,1L,2L,3L,4L), labels=c('recipient has no sequence','recipient sequenced HISTORIC only','recipient in p24 tree HISTORIC only','recipient sequenced PANGEA','recipient in p24 tree PANGEA'))])
	#
	rpi[, TR_S:= as.integer(PHP24_R_N>0)]
	set(rpi, rpi[, which(PHP24_TR_N>0)], 'TR_S', 1L)	
	set(rpi, rpi[, which(PHP24_MTR_N>0)], 'TR_S', 2L)
	set(rpi, rpi[, which(PHP24_MTR_PANGEA>0)], 'TR_S', 3L)
	set(rpi, rpi[, which(!is.na(PHP24_MTR_MINPD) & !is.na(PHP24_MTR_PANGEA_MINPD) & PHP24_MTR_MINPD==PHP24_MTR_PANGEA_MINPD)], 'TR_S', 4L)
	set(rpi, NULL, 'TR_S', rpi[, factor(TR_S, levels=c(0L,1L,2L,3L,4L), labels=c('recipient not in p24 tree','recipient in p24 tree has only female pr tr','recipient in p24 tree has male pr tr in HISTORIC only','recipient in p24 tree with male pr tr in PANGEA that is not closest','recipient in p24 tree with male pr tr in PANGEA that is closest'))])	
	#	save
	save(rp, rpi, file=file.path(wdir, paste('150825_FastTree_gagAllSmall_pairs_patristic',pd.max*100,'pc.rda',sep='')))	

	 #	plot
	ggplot(rpi, aes(x=REC_S, fill=as.factor(FIRSTPOSDATEc))) + geom_bar() + coord_flip() +
			theme_bw() + theme(legend.position='bottom') +
			labs(y='count',x='',fill='date first positive', title='RCCS female seroconverter\n(seroconversion window <= 2yrs)\n')
	ggsave(file=file.path(wdir,paste('150825_FastTree_Recipients_SCw2_patristic',pd.max*100,'.pdf',sep='')), h=7, w=7)
	
	ggplot(rpi, aes(x=TR_S, fill=as.factor(FIRSTPOSDATEc))) + geom_bar() + coord_flip() +
			theme_bw() + theme(legend.position='bottom') + facet_grid(.~REC_S) +
			scale_y_continuous(minor_breaks=seq(0,1e3,10)) +
			labs(y='count',x='',fill='date first positive', title=paste('pairs with RCCS female seroconverter\n(seroconversion window <= 2yrs)\n(patristic distance <',100*pd.max,'%)\n'))
	ggsave(file=file.path(wdir,paste('150825_FastTree_Pairs_SCw2_patristic',pd.max*100,'.pdf',sep='')), h=7, w=14)
	
	
	rpi[, table(REC_S, TR_S)]	
	subset(rpi, REC_S=='recipient in p24 tree PANGEA' & TR_S=='recipient in p24 tree with male pr tr in PANGEA that is closest')
	#	34		at 4%
	#	24	 	at 2%
	subset(rpi, REC_S%in%c('recipient in p24 tree HISTORIC only','recipient in p24 tree PANGEA') &
				TR_S%in%c('recipient in p24 tree with male pr tr in PANGEA that is not closest','recipient in p24 tree with male pr tr in PANGEA that is closest','recipient in p24 tree has male pr tr in HISTORIC only'))
	#	191		at 4%
	#	121		at 2%

	#	extract phylotype candidates
	pty.rp	<- subset(rp, !is.na(PID) & !is.na(TR_PID) & !is.na(TR_SID) & !is.na(SID) & TR_SEX=='M')
	#	add geographic info on transmitters
	#	TODO have these been moving in the RCCS?
	tmp2	<- subset(rd, !is.na(PID) & SEX=='M', select=c(RID, REGION, COMM_NUM, HH_NUM, BIRTHDATE))
	setkey(tmp2, RID)
	tmp2	<- unique(tmp2)
	setnames(tmp2, c('RID','REGION','COMM_NUM','HH_NUM','BIRTHDATE'), c('TR_RID','TR_REGION','TR_COMM_NUM','TR_HH_NUM','TR_BIRTHDATE'))
	pty.rp	<- merge(pty.rp, tmp2, by='TR_RID')	
	setkey(pty.rp, PTY_RUN, RID)
	write.csv(pty.rp, row.names=FALSE, file=file.path(wdir, paste('150830_FastTree_Pairs_SCw2_patristic',pd.max*100,'_PANGEA_M2F_hypothetical_pairs.csv',sep='')))
	save(pty.rp, file=file.path(wdir, paste('150830_FastTree_Pairs_SCw2_patristic',pd.max*100,'_PANGEA_M2F_hypothetical_pairs.rda',sep='')))
	#
	#	TODO add circumcision status on transmitters (in rt --> now move to timelines)
	#	TODO add ART trajectories
	#
	
	#	check whereabouts of those not yet processed
	tmp		<- unique(subset(rt, SEQ_TYPE=='Sanger not started by Jul2016', select=c(RID, ROUND, SEX)))
	tmp		<- merge(tmp, subset(rd, select=c(RID, PID)), by='RID')
	write.table(tmp, file.path(wdir, 'RCCSparticipants_sequencingunclear.csv'))
	infile.sangerstats	<- "~/Dropbox (Infectious Disease)/PANGEA_data/2016-07-07_PANGEA_stats_by_sample.csv"
	tmp2	<- as.data.table(read.csv(infile.sangerstats, stringsAsFactors=FALSE))	
	setnames(tmp2, colnames(tmp2), gsub('\\.','_',toupper(colnames(tmp2))))
	setnames(tmp2, c('PROJECTID','STATUS'), c('PID','SANGER_STATUS'))
	merge(tmp, tmp2, by='PID')	# they are not there!
}
######################################################################################
######################################################################################
hivc.db.Date2numeric<- function( x )
{
	if(!class(x)%in%c('Date','character'))	return( x )
	x	<- as.POSIXlt(x)
	tmp	<- x$year + 1900
	x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
	x	
}
######################################################################################
######################################################################################
project.Rakai.aliRegion1.597<- function()
{	
	require(big.phylo)
	require(plyr)
	
	#	load 597 new files
	infile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/PANGEA_UG full alignment/PANGEA_HIV_n597_Imperial_v160916_RakaiPlates.fasta'
	sqn		<- read.dna(infile, format='fa')
	sqni	<- data.table(TAXA=rownames(sqn))
	sqni	<- subset(sqni, grepl('HXB2|consensus',TAXA))
	sqni[, SID:= gsub('_consensus','',TAXA)]
	set(sqni, sqni[, which(grepl('HXB2',TAXA))],'SID',NA_character_)
	#merge(sqni, unique(subset(rs, !is.na(SID), select=c(RID,SEQTYPE,SID,PID))), by='SID',all.x=1)
	sqn		<- sqn[ sqni[,TAXA], ]
	
	#	load PANGEA alignment
	infile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_HXB2.fasta'
	sqp		<- read.dna(infile,format='fa')
	
	#	required: HXB2 in alignment
	ans		<- seq.align.based.on.common.reference(sqn, sqp, return.common.sites=TRUE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_597_p_HXB2.fasta'
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
	load('~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_597_p_HXB2.rda')
	sqi		<- data.table(TAXA=rownames(sq), ID=seq_len(nrow(sq)))
	
	infile	<- '~/Dropbox (Infectious Disease)/PANGEA_data/2016-09-18_PAN_SANGER_IDs.txt'
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
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_597_p_HXB2.fasta'
	write.dna(sq, file=outfile, format='fasta', colsep='', nbcol=-1)	
	#
	#	read COMETv0.5
	#
	infile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_597_p_HXB2_COMETv0.5.txt'
	sqi		<- as.data.table(read.table(infile, header=TRUE, sep='\t',stringsAsFactors=FALSE))
	setnames(sqi, c('name','subtype','length'), c('TAXA','COMET_ST','COMET_N'))
	sqi[, COMET_V:='0.5']
	#
	#	read COMETv2.1
	#
	infile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_597_p_HXB2_COMETv2.1.txt'
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
	wdir	<- "~/Dropbox (Infectious Disease)/Rakai Fish Analysis/circumcision"
	#	load Susanna s codon alignment
	susa.f	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah3/Region1_codon aligned.fasta'
	susa.s	<- read.dna(susa.f, format='fasta')
	susa.d	<- data.table(	ID=gsub('\\*.*','',rownames(susa.s)),
							DATA= factor(grepl('^PG', rownames(susa.s)), levels=c(TRUE,FALSE), label=c('PNG','LNL')))
	infile.relabel	<- "~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional phylogenetic analyses/Region 1 gag analysis/SummaryofGagSequenceData.rda"				
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
	infile.gag	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Rakai Data for IqTree/gag.sqn.fasta'
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
	#	move to '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments'
	#
	file.rename(outfile, file.path('~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments',basename(outfile)))
}
######################################################################################
#
######################################################################################
project.Rakai.aliRegion1.merge.with.PANGEA.160825<- function()
{	
	#	load Susanna s codon alignment and remove any PANGEA seqs
	susa.f	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151129.fasta'
	susa.f	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_missing_p_HXB2.fasta'
	susa.s	<- read.dna(susa.f, format='fasta')
	susa.d	<- data.table(TAXA=rownames(susa.s), DATA= factor(grepl('^PG', rownames(susa.s)), levels=c(TRUE,FALSE), label=c('PNG','LNL')))	
	in.s	<- susa.s[subset(susa.d, DATA!='PNG')[, TAXA],]	
	#	load PANGEA alignment	
	png.f	<- '~/Dropbox (Infectious Disease)/Rakai Fish Analysis/PANGEA_orig/PANGEA_HIV_n4562_Imperial_v151113_GlobalAlignment.rda'
	load(png.f)	#sq, sqi, si
	sq		<- sq[grepl('^PG[0-9]+|HXB2',rownames(sq)),]
	#	required: HXB2 in alignment
	sqn		<- seq.align.based.on.common.reference(in.s, sq, return.common.sites=TRUE, regexpr.reference='HXB2', regexpr.nomatch='-|\\?')
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1_codonaligned_p_PANGEA151113_p_HXB2.fasta'
	write.dna(sqn, file=outfile, format='fasta', colsep='', nbcol=-1)
	
	#
	#	subset to Rakai and relabel
	#
	
	#	load summaryData
	infile.relabel	<- "~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional phylogenetic analyses/Region 1 gag analysis/SummaryofGagSequenceData.rda"				
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
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEA_alignments/Regional Alignments/150825_Region1UG_codonaligned_p_PANGEA151113_p_HXB2.fasta'
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
	write.dna( sqn, format='fasta', file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201.fasta', colsep='', nbcol=-1)
	save(sqni, sqn, file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201.rda')
	
	sqn				<- read.dna(file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201_rm99gps.fasta',format='fasta')
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
	save(sqni, sqn, file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201_rm99gps.rda')	
	write.dna( sqn, format='fasta', file='~/Dropbox (Infectious Disease)/Rakai Fish Analysis/Susannah/PANGEA_Region1_Final2_151201_rm99gps.fasta', colsep='', nbcol=-1)
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
	
	f.arv	<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_150909.rda'
	f.rccsid<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/150625_PangeaBoxes.csv'
	f.sid	<- '~/Dropbox (Infectious Disease)/PANGEA_data/SangerUpdates/2015-07-20_PANGEA_3.csv'
	f.seq	<- '~/Dropbox (Infectious Disease)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908.fasta'
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
	outfile			<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_PotentialDRMs_OR_150910.R'
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
	
	f.arv			<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_150911.rda'
	f.rccsid		<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/data/150625_PangeaBoxes.csv'
	f.sid			<- '~/Dropbox (Infectious Disease)/PANGEA_data/SangerUpdates/2015-07-20_PANGEA_3.csv'
	f.seq			<- '~/Dropbox (Infectious Disease)/PANGEA_data/PANGEAconsensuses_2015-09_Imperial/PANGEA_HIV_n4562_Imperial_v150908.fasta'
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
	outfile			<- '~/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_Fisherfolk/PANGEA_ARV/RakaiARVData_PotentialDRMs_OR_150911.R'
	tmp				<- big.phylo:::seq.rm.drugresistance(seq, outfile=outfile)
	nodr.info		<- tmp$nodr.info
	nodr.seq		<- tmp$nodr.seq
}