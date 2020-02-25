femalenetworks.get.outdegree <- function()
{
	#	load networks
	load('~/sandbox/DeepSeqProjects/RakaiPopSample_phyloscanner_analysis_190706/Rakai_phscnetworks_190706.rda')
	dnet <- as.data.table(dnet)
	
	#	calculate out-degree of all individuals
	dg <- subset(dnet, CATEGORISATION=='close.and.adjacent.cat' & TYPE%in%c('12','21'))
	tmp <- subset(dnet, CATEGORISATION=='close.and.adjacent.cat' &  TYPE=='complex.or.no.ancestry')
	set(tmp, NULL, 'SCORE', tmp[,SCORE/2])
	set(tmp, NULL, 'TYPE', '12')
	dg <- rbind(dg, tmp)
	set(tmp, NULL, 'TYPE', '21')
	dg <- rbind(dg, tmp)
	dg <- dg[, list(SCORE=sum(SCORE)), by=c('H1','H2','TYPE')]
	subset(dg, TYPE=='12')
	tmp <- dg[TYPE=='12', list(OUTDEGREE_S=sum(SCORE), OUTDEGREE_N=length(unique(H2))), by='H1']
	dg <- dg[TYPE=='21', list(OUTDEGREE_S=sum(SCORE), OUTDEGREE_N=length(unique(H1))), by='H2']
	setnames(dg, 'H2', 'AID')
	setnames(tmp, 'H1', 'AID')
	dg <- rbind(dg, tmp)
	dg <- dg[, list(OUTDEGREE_S=sum(OUTDEGREE_S), OUTDEGREE_N=sum(OUTDEGREE_N)), by='AID']
	
	#	load selected individuals
	df <- as.data.table(read.csv('~/Box Sync/OR_Work/2019/2019_RCCS_high_risk_women/191011_femaleonly_primary_selection.csv', stringsAsFactors=FALSE))
	df[, SELECT:='primary']
	tmp <- as.data.table(read.csv('~/Box Sync/OR_Work/2019/2019_RCCS_high_risk_women/191011_femaleonly_backup_selection.csv', stringsAsFactors=FALSE))
	tmp[, SELECT:='backup']
	df <- rbind(df,tmp)
	set(df, NULL, c('X','SEX'), NULL)	
	dg <- merge(df, dg, by='AID')
	
	outfile <- '~/Box Sync/OR_Work/2019/2019_RCCS_high_risk_women/191011_femaleonly_outdegrees.csv'
	write.csv(dg, file=outfile)
	
	dfa	<- as.data.table(read.csv("~/Box/OR_Work/PANGEA/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_anonymised_RIDs.csv", stringsAsFactors=FALSE))
	dfa	<- unique(subset(dfa, select=c(ID,AID)))		
	setnames(dfa, 'ID', 'RID')
	dg <- merge(dfa, dg, by='AID')
	
	rd <- merge(rd, dfa, by='RID')
	rd <- subset(rd, select=c(RID, AID, SEX, COMM_NUM_A))
	rd[, COMM_TYPE:= as.character(factor(substr(COMM_NUM_A,1,1)=='f', levels=c(TRUE,FALSE), labels=c('fishing','inland')))]
	rd <- subset(rd, select=c(AID, COMM_TYPE))
	
}


RakaiFull.analyze.trmpairs.femaleoutdegrees<- function()
{
	infile				<- '~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/todi_pairs_171122_cl25_d50_prior23_min30_networksallpairs.rda'		
	outfile.base		<- gsub('_networksallpairs.rda','',infile)		
	
	confidence.cut		<- 0.6
	load(infile)
	
	#
	#	calculate females with high out degree
	#
	if(1)
	{
		outfile	<- paste0(outfile.base,'_transmission_networks_high_out_degree_females.csv')
		#	of all edges in the transmission network, keep 
		#	the more likely direction
		rtm	<- rtn[, list(TYPE= TYPE[ which.max(KEFF) ]), by=c('ID1','ID2','PTY_RUN','ID1_SEX','ID2_SEX')]
		rtm	<- subset(rtm, !grepl('disconnected',TYPE))
		tmp	<- subset(rtm, TYPE=='ambiguous')
		set(tmp, NULL, 'TYPE', '12')
		rtm	<- rbind(rtm, tmp)
		set(tmp, NULL, 'TYPE', '21')
		rtm	<- rbind(rtm, tmp)
		rtm	<- subset(rtm, TYPE!='ambiguous')
		tmp	<- subset(rtm, TYPE=='21')
		setnames(tmp, c('ID1','ID2','ID1_SEX','ID2_SEX'), c('ID2','ID1','ID2_SEX','ID1_SEX'))
		set(tmp, NULL, 'TYPE', '12')
		rtm	<- rbind(subset(rtm, TYPE=='12'), tmp)
		dg	<- rtm[, list(OUT_DEGREE=length(ID2)), by=c('ID1','ID1_SEX')]
		dg	<- dg[order(ID1_SEX,-OUT_DEGREE),]
		dg	<- subset(dg, ID1_SEX=='F' & OUT_DEGREE>1)	
		setnames(dg, c('ID1','ID1_SEX','OUT_DEGREE'), c('ID','SEX','DEGREE'))
		# 30 females
	}
	#
	#	calculate females with high in+out degree
	#
	if(0)
	{
		outfile	<- paste0(outfile.base,'_transmission_networks_high_degree_females.csv')
		#	get igraph from all edges in the transmission network
		#	count linkages regardless of direction
		tmp		<- unique(subset(rtn, select=c(ID1, ID2)))			
		dg		<- graph.data.frame(tmp, directed=FALSE, vertices=NULL)
		#	calculate degrees
		dg		<- degree(dg)
		dg		<- data.table(ID=names(dg), DEGREE=as.numeric(dg))
		#	get gender	
		tmp		<- unique(subset(rd, select=c(RID,SEX)))
		setnames(tmp, 'RID', 'ID')
		dg		<- merge(dg, tmp, by='ID')
		#	sort by gender, degree
		dg		<- dg[order(SEX,-DEGREE),]
		dg		<- subset(dg, SEX=='F' & DEGREE>3)		
		# 41 females		
	}
	#
	#	write to csv file
	#
	write.csv(dg, row.names=FALSE, file=outfile)	
	
	#
	#	plot transmission networks with high degree females shown as triangles
	#
	idclus	<- rtn[ ID1%in%dg$ID | ID2%in%dg$ID, sort(unique(IDCLU)) ]	# 25 networks
	pns		<- lapply(seq_along(idclus), function(i)
			{
				idclu	<- idclus[i]
				di		<- unique(subset(rd, select=c(RID,SEX)))
				setnames(di, 'RID', 'ID')
				di		<- merge(di, subset(dg, select=c(ID,DEGREE)), by='ID', all.x=TRUE)
				di[, HIGH_DEGREE:= as.character(as.integer(!is.na(DEGREE)))]				
				df		<- subset(rtn, IDCLU==idclu)
				set(df, NULL, c('ID1_SEX','ID2_SEX'), NULL)
				p		<- phsc.plot.transmission.network(df, di, point.size=10, 
						edge.gap=0.04, 
						edge.size=0.4, 
						curvature= -0.2, 
						arrow=arrow(length=unit(0.04, "npc"), type="open"), 
						curv.shift=0.06, 
						label.size=3, 
						node.label='ID', 			 
						node.fill='SEX', 
						node.shape='HIGH_DEGREE',
						node.shape.values=c('0'=16, '1'=17), 
						node.fill.values=c('F'='hotpink2', 'M'='steelblue2'),
						threshold.linked=0.6)	
				p	
			})
	pdf(file=gsub('csv','pdf',outfile), w=8, h=8)
	for(i in seq_along(pns))	
		print(pns[[i]])
	dev.off()	
}