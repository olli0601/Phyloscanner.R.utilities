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

make.map.190129	<- function()
{
	require(data.table)
	require(rgdal)
	require(rgeos)
	library(raster)
	
	# load data
	infile	<- '~/Box Sync/OR_Work/2018/2018_RakaiViralLoad/data/merged_round17_vl_gps.rda'
	load(infile)
		
	#convert the data into a data table
	dt<- as.data.table(ds)
	dt<- dt[,.(RCCS_STUDYID, SEX, AGEYRS, HIV_STATUS, LATITUDE_JITTER, LONGITUDE_JITTER, VL_COPIES, VL_UNDETECTABLE)]
	#set the NA VL to 0
	dt[is.na(VL_COPIES), VL_COPIES:=0]
	dt[,VL_DETECTABLE := as.numeric(VL_COPIES>=1000)]
	dt[,RCCS_STUDYID2:= seq_len(nrow(dt)) ]
		
	#################################################### load in maps
	# Load in Uganda Shape files 
	uganda1<-raster::getData('GADM',country="UGA",level=1)# Admin unit 1
	uganda3<- raster::getData('GADM', country='UGA', level=3)
	rakai1<-subset(uganda1, NAME_1=="Rakai")
	rakai3<- subset(uganda3, NAME_1=="Rakai")
	masaka1<-subset(uganda1, NAME_1=="Masaka")
	# Create a smaller Rakai for plotting (not current Rakai region no longer includes kabula subdistrict 3)
	#minirak<-rakai3[which(rakai3$NAME_2!="Kabula" | rakai3$NAME_3=="Lyantonde Tc" | rakai3$NAME_3=="Lyantonde"),]
	minirak<-rakai3[which(rakai3$NAME_2!="Kabula" | rakai3$NAME_3=="Lyantonde Tc"),]
	
	####################################################### Convert the data to meters
	#set the coordinates of the data
	coordinates(dt)<- ~ LONGITUDE_JITTER+LATITUDE_JITTER
	#set coordinate system to match uganda files
	proj4string(dt) <- proj4string(uganda1)
	
	#convert to m in order to build a 30x30m grid
	newcrs <- CRS("+proj=robin +datum=WGS84")
	dtnew<- spTransform(dt, newcrs)
	rakai1trans<- spTransform(rakai1, newcrs)
	miniraktrans<- spTransform(minirak, newcrs)
	masaka1trans<- spTransform(masaka1, newcrs)
	
	###################################################### Build Grid
	#Combine rakai1trans and masaka1trans
	outline<- union(rakai1trans, masaka1trans)
	#find the extent of the data
	exnew<- extent(dtnew)
	#extent of the maps
	exmap<- extent(outline)
	
	#chose extent to cover all the data and rakai district
	
	#With a 30m grid, I think the same individuals are usually entering calculations for a large number of grid points
	#Do we really need a 30m grid? Why not 100m?

	grid<- raster(xmn=min(exnew[1], exmap[1]), xmx= exnew[2], ymn=exmap[3], ymx=exnew[4], res=100 )
	grid[]<- 1:ncell(grid)
	
	# set the coordinate reference system to match
	proj4string(grid)<- proj4string(dtnew) 
	
	#restrict grid to map
	gridmask<- mask(grid, outline)
	
	#plot(gridmask)
	
	#consider the grid points in a data frame
	id<- as.data.table(1:ncell(gridmask))
	setnames(id, "V1", "ID")
	griddf<- as.data.table(SpatialPoints(grid))
	griddf<- data.table(id, griddf)
	setnames(griddf, gsub('y','LAT_GRID',gsub('x','LONG_GRID',colnames(griddf))))
	
	bw			<- 3000
	#require(mvtnorm)
	#dmvnorm( c(3.84,0) )	# ~ 9.996634e-05 
	threshold	<- bw*3.84 	# cut if density is < 1e-4
	threshold	<- threshold*threshold	# square the threshold, to avoid sqrt calculations in loop 		
	
	tmp			<- griddf[1:1e4,]
	anst<- system.time({
		ans	<- tmp[, {
					z1	<- LONG_GRID - dtnew@coords[,'LONGITUDE_JITTER']
					z2	<- LAT_GRID - dtnew@coords[,'LATITUDE_JITTER']
					z1	<- z1*z1 + z2*z2 		# square distance
					z2	<- which(z1<threshold)	# avoid sqrt on 2e4 entries
					# to avoid very large output data, calculate directly all smooths here
					z1	<- sqrt(z1[z2])			# sqrt on few entries					
					w	<- dnorm(z1, mu=0, sd=bw)
					# code assumes @coords and @data has same order. 
					list( 	HIV_STATUS_MEAN=mean( dtnew@data$HIV_STATUS[z2] ),				#no weighting by distance
							HIV_STATUS_KERNEL=sum( dtnew@data$HIV_STATUS[z2]*w )/sum(w),	#Gaussian kernel
							)
				}, by=c('ID','LONG_GRID','LAT_GRID')]
	})
}


