date2numeric<- function( x )
{
	if(!class(x)%in%c('Date','character'))	return( x )
	x	<- as.POSIXlt(x)
	tmp	<- x$year + 1900
	x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
	x	
}


vl.get.data.round17<- function()
{
	require(data.table)
	prjdir	<- '~/Box Sync/OR_Work/2018/2018_RakaiViralLoad'
	infile	<- 'data_raw/ViralLoad_Data_Pangea_Ratmann.rda'
	load( file.path(prjdir, infile) )
	
	# get individuals surveyed in round 17
	ds		<- subset(as.data.table(survey_data), visit==17)
	# reset dates from Date format to numeric
	for(x in c('visit_date','lastNegDate','firstPosDate'))
	{
		set(ds, NULL, x, date2numeric(ds[[x]]))
	}
	# make all column names upper case
	setnames(ds, colnames(ds), toupper(colnames(ds)))
	
	
	# merge GPS coordinates to surveyed individuals
	dg	<- as.data.table(gpsdat)	
	# bring dates into common format
	setnames(dg, colnames(dg), gsub('\\.','_',toupper(colnames(dg))))
	tmp	<- which(dg[, grepl('([0-9]+)/([0-9]+)/([0-9]+)',GPS_DATE)])
	set(dg, tmp, 'GPS_DATE', dg[tmp,gsub('([0-9]+)/([0-9]+)/([0-9]+)','\\3-\\1-\\2',GPS_DATE)])
	tmp	<- which(dg[, grepl('([0-9]+)-([A-Za-z]+)-([0-9]+)',GPS_DATE)])
	set(dg, tmp, 'GPS_DATE', dg[tmp,gsub('([0-9]+)-([A-Za-z]+)-([0-9]+)','20\\3-\\2-\\1',GPS_DATE)])	
	set(dg, NULL, 'GPS_DATE', dg[,gsub('Nov','11',gsub('Oct','10',gsub('Sep','09',gsub('Aug','08',gsub('July','07',gsub('Jun','06',gsub('May','05',GPS_DATE)))))))])
	# reset dates from character format to numeric
	set(dg, NULL, 'GPS_DATE', date2numeric(dg[,GPS_DATE]))
	# make households per date unique
	dg	<- unique(dg, by=c('HHID','GPS_DATE'))
	
	# merge surveyed individuals from round 17. a bit tricky
	tmp	<- unique(subset(ds, select=c(RCCS_STUDYID, VISIT_DATE, HHID)))
	tmp	<- merge(tmp, dg, by='HHID', all.x=TRUE)
	# TODO check households with missing geocodes
	ch	<- subset(tmp, is.na(LATITUDE_JITTER) | is.na(LATITUDE_JITTER))
	write.csv(ch, file=file.path(prjdir,'data/check_missing_coordinates.csv'))
	# for every individual, extract house closest in time
	tmp		<- subset(tmp, !is.na(LATITUDE_JITTER) & !is.na(LATITUDE_JITTER))
	tmp2	<- tmp[, list(GPS_DATE= GPS_DATE[which.min(abs(GPS_DATE-VISIT_DATE))[1]]), by=c('RCCS_STUDYID','VISIT_DATE')]
	tmp		<- merge(tmp, tmp2, by=c('RCCS_STUDYID','VISIT_DATE','GPS_DATE'))
	stopifnot(nrow(tmp)==nrow(tmp2))
	set(tmp, NULL, c('COMM','HOUSE'), NULL)
	#	TODO we are losing some individuals here, likely due to missing GPS. double check when above fixed.
	ds		<- merge(ds, tmp, by=c('RCCS_STUDYID','VISIT_DATE','HHID'))
	
	# extract viral loads from round 17
	dvl		<- subset(as.data.table(viralLoads), visit==17)
	setnames(dvl, colnames(dvl), toupper(colnames(dvl)))
	setnames(dvl, c('DATE','COPIES','DONEBY'), c('VL_DATE','VL_COPIES','VL_DONEBY'))
	set(dvl, NULL, 'VL_DATE', date2numeric(dvl[,VL_DATE]))
	set(dvl, NULL, 'VL_UNDETECTABLE', dvl[, as.integer(VL_COPIES<1000)])
	# merge surveyed individuals from round 17. a bit tricky
	tmp	<- unique(subset(ds, select=c(RCCS_STUDYID, VISIT_DATE)))
	# TODO we are losing some individuals here, could be because of error above. double check when fixed.
	dvl	<- merge(dvl, tmp, by='RCCS_STUDYID')
	stopifnot(nrow(dvl)==nrow(unique(dvl, by=c('RCCS_STUDYID','VISIT','VL_DATE'))))
	set(dvl, NULL, 'VISIT', dvl[, as.numeric(VISIT)])
	# ok, one viral load per individual. 
	# merge viral load with surveillance data.
	ds	<- merge(ds, dvl, by=c('RCCS_STUDYID','VISIT','VISIT_DATE'), all.x=TRUE)
	
	# check if viral load for all infected
	ch	<- subset(ds, HIV_STATUS==1 & is.na(VL_COPIES))
	write.csv(ch, file=file.path(prjdir,'data/check_missing_viralloads.csv'))
	
	ds	<- subset(ds, HIV_STATUS==0 | (HIV_STATUS==1 & !is.na(VL_COPIES)))
	save(ds, file=file.path(prjdir,'data/merged_round17_vl_gps.rda'))
}


