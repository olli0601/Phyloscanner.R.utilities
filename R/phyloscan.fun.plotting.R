#' @export
#' @import data.table grid ggtree ggnet
#' @title Plot probability network with most likely edges
#' @description This function plots the network showing the most likely edge types.  
#' @param df data.table with the following columns  "IDCLU","ID1", "ID2", "TYPE","KEFF","LKL_MAX","POSTERIOR_SCORE" 
#' @param di data.table with meta-data to customize the plot with columns  "ID", node.shape, node.label, node.fill 
#' @param point.size size of the individual points
#' @param point.sizec.couple size of the outer ring around individuals in couples
#' @param edge.gap value to adjust start / end points of edges
#' @param edge.size multiplier by which the size of edges is shrunk/magnified
#' @param curvature curvature of directed edges  
#' @param arrow type of arrow to be plotted
#' @param curv.shift offset to place the label for directed edges
#' @param label.size size of label
#' @param node.label Text displayed on top of each node 
#' @param node.shape column name in di by which the shape of each node is drawn 
#' @param node.fill column name in di by which each node is coloured
#' @param node.shape.values named vector associating shapes to the values in the node.shape column
#' @param node.fill.values named vector associating colours to the values in the node.fill column
#' @param threshold.linked treshold value between 0 and 1. Edges with weight above this treshold are shown in black.
#' @return ggplot object
phsc.plot.maxedge.network<- function(df, di, point.size=10, point.size.couple=point.size*1.4, edge.gap=0.04, edge.size=0.4, curvature= -0.2, arrow=arrow(length=unit(0.04, "npc"), type="open"), curv.shift=0.08, label.size=3, node.label='ID', node.shape=NA_character_, node.fill=NA_character_, node.shape.values=NA_integer_, node.fill.values=NA_character_, edge.label=NA_character_, edge.label.values=NA_character_, threshold.linked=NA_real_)
{	
	#point.size=10; point.size.couple=14; edge.gap=0.04; edge.size=0.4; curvature= -0.2; arrow=arrow(length=unit(0.04, "npc"), type="open"); curv.shift=0.08; label.size=3
	#node.label='ID'; threshold.linked=0.6; node.shape=NA_character_; node.fill=NA_character_; node.shape.values=NA_integer_; node.fill.values=NA_character_
	#node.shape='IN_COUPLE'; node.fill='SEX'
	#node.fill.values=c('F'='hotpink2', 'M'='steelblue2')
	#node.shape.values=c('not in long-term\nrelationship'=18,'in long-term\nrelationship'=16)
	if(is.na(node.label))
	{
		node.label<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.label, NA_character_)
	}
	if(is.na(edge.label))
	{
		edge.label<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(df, NULL, edge.label, df[, TYPE])
		edge.label.values	<- c('12'='red','21'='red','ambiguous'='grey50')
	}
	if(is.na(node.shape))
	{
		node.shape<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.shape, 'NA')
	}
	if(is.na(node.fill))
	{
		node.fill<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.fill, 'NA')
	}
	if(any(is.na(node.fill.values)))
	{
		z						<- unique(di[[node.fill]])
		node.fill.values		<- heat.colors(length(z))
		names(node.fill.values)	<- z
	}
	if(any(is.na(node.shape.values)))
	{
		z						<- unique(di[[node.shape]])
		node.shape.values		<- seq_along(z)
		names(node.shape.values)<- z
	}
	setnames(di, c(node.label, node.shape, node.fill), c('NODE_LABEL','NODE_SHAPE','NODE_FILL'))
	tmp	<- c('NODE_LABEL','NODE_SHAPE','NODE_FILL')[which(c(node.label, node.shape, node.fill)=='ID')]
	if(length(tmp))
		set(di, NULL, 'ID', di[[tmp]])
	di	<- subset(di, select=c(ID, NODE_LABEL, NODE_SHAPE, NODE_FILL))
	setnames(df, c(edge.label), c('EDGE_LABEL'))
	
	layout	<- as.data.table(ggnet2(network(unique(subset(df, select=c(ID1,ID2))), directed=FALSE, matrix.type="edgelist"))$data[,c("label", "x", "y")])
	setnames(layout, c('label','x','y'), c('ID1','ID1_X','ID1_Y'))
	df		<- merge(df, layout, by='ID1')
	setnames(layout, c('ID1','ID1_X','ID1_Y'), c('ID2','ID2_X','ID2_Y'))
	df		<- merge(df, layout, by='ID2')
	setnames(layout, c('ID2','ID2_X','ID2_Y'),  c('ID','X','Y'))	
	layout	<- merge(layout,di, by='ID')	
	
	df[, EDGETEXT_X:= (ID1_X+ID2_X)/2]
	df[, EDGETEXT_Y:= (ID1_Y+ID2_Y)/2]
	#
	#	calculate score for linked
	if(is.na(threshold.linked))
	{
		df	<- merge(df,df[, 	{
							z<- rep('alpha_1', length(TYPE))
							z[which.max(POSTERIOR_SCORE)]	<- 'alpha_2'
							list(ALPHA=z, TYPE=TYPE)	
						}, by=c('ID1','ID2')], by=c('ID1','ID2','TYPE'))		
	}
	if(!is.na(threshold.linked))
	{
		tmp	<- subset(df, TYPE!='not close/disconnected')[, list( ALPHA=as.character(factor(sum(POSTERIOR_SCORE)>=threshold.linked, levels=c(TRUE, FALSE), labels=c('alpha_2','alpha_1'))) ), by=c('ID1','ID2')]
		df	<- merge(df, tmp, by=c('ID1','ID2'))		
	}	
	#	for edges, move the start and end points on the line between X and Y
	#	define unit gradient
	df[, MX:= (ID2_X - ID1_X)]	
	df[, MY:= (ID2_Y - ID1_Y)]
	tmp		<- df[, sqrt(MX*MX+MY*MY)]
	set(df, NULL, 'MX', df[, MX/tmp])
	set(df, NULL, 'MY', df[, MY/tmp])	
	set(df, NULL, 'ID1_X', df[, ID1_X + MX*edge.gap])
	set(df, NULL, 'ID1_Y', df[, ID1_Y + MY*edge.gap])
	set(df, NULL, 'ID2_X', df[, ID2_X - MX*edge.gap])
	set(df, NULL, 'ID2_Y', df[, ID2_Y - MY*edge.gap])	
	#	label could just be move on the tangent vector to the line
	#	define unit tangent
	df[, TX:= -MY]
	df[, TY:= MX]
	tmp		<- df[, which(TYPE=='12')]
	set(df, tmp, 'EDGETEXT_X', df[tmp, EDGETEXT_X + TX*curv.shift])
	set(df, tmp, 'EDGETEXT_Y', df[tmp, EDGETEXT_Y + TY*curv.shift])
	tmp		<- df[, which(TYPE=='21')]
	set(df, tmp, 'EDGETEXT_X', df[tmp, EDGETEXT_X - TX*curv.shift])
	set(df, tmp, 'EDGETEXT_Y', df[tmp, EDGETEXT_Y - TY*curv.shift])
	
	tmp		<- df[, list(TYPE=TYPE[which.max(KEFF)]), by=c('ID1','ID2')]
	df		<- merge(df, tmp, by=c('ID1','ID2','TYPE'))
	#
	p		<- ggplot() +			
			geom_point(data=layout, aes(x=X, y=Y, fill=NODE_FILL, pch=NODE_SHAPE), size=point.size) +
			geom_segment(data=subset(df, TYPE=='ambiguous' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, alpha=ALPHA, colour=EDGE_LABEL), lineend="butt") +
			geom_curve(data=subset(df, TYPE=='12' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, alpha=ALPHA, colour=EDGE_LABEL), curvature=curvature, arrow=arrow, lineend="butt") +
			geom_curve(data=subset(df, TYPE=='21' & KEFF>0), aes(x=ID2_X, xend=ID1_X, y=ID2_Y, yend=ID1_Y, size=edge.size*KEFF, alpha=ALPHA, colour=EDGE_LABEL), curvature=curvature, arrow=arrow, lineend="butt") +
			scale_shape_manual(values=c(node.shape.values, 'NA'=16)) +
			scale_fill_manual(values=c(node.fill.values, 'NA'='grey50')) +
			scale_alpha_manual(values=c('alpha_1'=0.5,'alpha_2'=1, 'NA'=0)) +
			scale_colour_manual(values=c(edge.label.values, 'NA'='grey50')) +
			scale_size_identity() +
			geom_text(data=subset(df, TYPE!='not close/disconnected' & KEFF>0), aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=paste0(round(100*POSTERIOR_SCORE,d=1),'%')), size=label.size) +
			geom_text(data=layout, aes(x=X, y=Y, label=NODE_LABEL)) +
			theme_void() +
			guides(colour='none', fill='none',size='none', pch='none') 
	layout		<- subset(layout, select=c(ID,X,Y))
	setnames(layout, c('ID','X','Y'), c('label','x','y'))	
	p$layout	<- layout
	p
}

#' @export
#' @import data.table grid ggtree
#' @title Plot short read phylogenies and highlight individuals
#' @description This function plots short read phylogenies and highlights the clades of two individuals in red and blue.  
#' @param phs List of trees in ape format
#' @param dfs data.table with mandatory column 'IDX' and optional column 'TITLE'. IDX is the index of all phylogenies in 'phs' that are to be plotted. TITLE is a title for each sub-plot, for example specifying a window.
#' @param ids Vector of regular expressions that identify individuals to be highlighted in colour.
#' @param plot.cols Vector of colours for each individual
#' @param group.redo Logical, indicating if the colour groups should be recalculated from ids.
#' @param drop.blacklisted Logical, indicating if all blacklisted taxa should be dropped prior to plotting.  
#' @param pdf.h	Height of the pdf file in inches.
#' @param pdf.rw Relative width of the pdf file, internally multiplied by the number of phylogenies to give the total width in inches.
#' @param pdf.ntrees Number of trees per pdf.
#' @param pdf.title.size Size of pdf title in inches.
#' @param plot.file If not missing, the phylogenies will be printed to file.	
#' @return List of ggtree objects, ready for printing.
phsc.plot.phy.selected.individuals<- function(phs, dfs, ids, plot.cols=rainbow(length(ids)), plot.file=NA, group.redo=FALSE, drop.less.than.n.ids=2, drop.blacklisted=FALSE, pdf.h=50, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40)
{	
	#	determine which phylogenies contain at least one of the requested individuals
	tmp		<- copy(dfs)
	tmp		<- merge(tmp, tmp[, {
						ph	<- phs[[ IDX ]]
						z	<- ph$tip.label
						if(drop.blacklisted)
							z	<- as.character(attr(ph, "INDIVIDUAL"))
						list(HAS_N_IND=length(which(sapply(ids, function(id) any(grepl(id, z)) ))))						  				
					}, by='IDX'], by='IDX')
	tmp		<- subset(tmp, HAS_N_IND>=drop.less.than.n.ids)
	phps	<- lapply(seq_len(nrow(tmp)), function(i){
				ph.title	<- NULL
				if('TITLE'%in%colnames(tmp))
					ph.title	<- tmp[i, TITLE]										
				ph			<- phs[[ tmp[i, IDX] ]]
				#
				#
				if(drop.blacklisted)
				{
					z						<- data.table(FROM=ph$edge[,1],TO=ph$edge[,2])
					z						<- merge(z, data.table(TO= seq_len(length(attr(ph,'INDIVIDUAL'))), GROUP= attr(ph,'INDIVIDUAL')), by='TO')
					z						<- subset(z, is.na(GROUP) & TO<=Ntip(ph))[, TO]
					ph						<- drop.tip(ph, z)
					attr(ph,'NODE_SHAPES')	<- rep(FALSE, Nnode(ph, internal.only=FALSE))
					attr(ph,'INDIVIDUAL')	<- NULL
					attr(ph,'SPLIT')		<- NULL
				}
				#
				#
				if(group.redo||drop.blacklisted)
				{
					phb			<- data.table(IDX=seq_along(ph$tip.label), TAXA=ph$tip.label, ID='none')
					for(id in ids)
						set(phb, phb[, which(grepl(id,TAXA))], 'ID', id)
					tmp				<- lapply( phb[, unique(ID)], function(x)	subset(phb, ID==x)[, TAXA]	)
					names(tmp)		<- phb[, unique(ID)]
					ph				<- groupOTU(ph, tmp, group='INDIVIDUAL')
					#z	<- merge(data.table(FROM=ph$edge[,1],IDX=ph$edge[,2]), phb, by='IDX', all=1)
					#z[, GROUP:= attr(ph,'INDIVIDUAL')[1:nrow(ph$edge)]]
					#z	<- unique(subset(z, !is.na(ID), select=c(ID, GROUP)))
					#attr(ph,'INDIVIDUAL')	<- factor(attr(ph,'INDIVIDUAL'), levels=c(0,z[,as.character(GROUP)]), labels=c('not characterized',z[,FILE_ID]))
					#ph									
				}				
				#
				#
				cols		<- rep('grey50', length(attr(ph, "INDIVIDUAL"))) 
				for(i in seq_along(ids))
					cols[ grepl(ids[i], attr(ph, "INDIVIDUAL")) ]<- plot.cols[i]
				attr(ph, 'COLOUR')	<- cols								
				#						
				p 			<- ggtree(ph, aes(color=I(COLOUR))) +
						#geom_point2(shape = 16, size=3, aes(subset=NODE_SHAPES)) +
						scale_fill_hue(na.value="black") +								
						theme(legend.position="none") +
						geom_tiplab(aes(col=I(COLOUR))) +
						theme_tree2() +
						theme(legend.position="bottom", plot.title = element_text(size=pdf.title.size)) + 
						ggplot2::xlim(0, max(node.depth.edgelength(ph)[1:Ntip(ph)])*1.3) +
						labs(x='subst/site', title=ph.title)
				#ggsave(plot.file, w=10, h=200, limitsize = FALSE)
				p
			})
	#
	#	single page plot
	#		
	if(!is.na(plot.file))					
	{
		if(length(phps)<=pdf.ntrees)
		{
			cat('Plotting to file', plot.file,'...\n')
			pdf(file=plot.file, w=pdf.rw*length(phps), h=pdf.h)
			grid.newpage()
			pushViewport(viewport(layout=grid.layout(1, length(phps))))
			for(i in seq_along(phps))
				print(phps[[i]], vp = viewport(layout.pos.row=1, layout.pos.col=i))
			dev.off()
		}
		if(length(phps)>pdf.ntrees)
		{
			pi	<- data.table(IDX=seq_along(phps))
			pi[, PLOT:= ceiling(IDX/pdf.ntrees)]
			pi[, PLOT_IDX:= (IDX-1)%%pdf.ntrees+1]
			pi[,{
						cat('Plotting to file', gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file),'...\n')
						pdf(file=gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file), w=pdf.rw*pdf.ntrees, h=pdf.h)
						grid.newpage()
						pushViewport(viewport(layout=grid.layout(1, pdf.ntrees)))
						for(i in seq_along(IDX))
							print(phps[[IDX[i]]], vp = viewport(layout.pos.row=1, layout.pos.col=PLOT_IDX[i]))
						dev.off()
					}, by='PLOT']
		}
	}
	phps	
}

#' @export
#' @import data.table grid ggtree ggnet
#' @title Plot most likely transmission chain
#' @description This function plots a most likely transmission chain, showing one edges with two labels: 
#' L (linked) indicates the proportion of deep-sequence phylogenies in whom the two individuals are phylogenetically close and adjacent.
#' D (direction) indicates the proportion of deep-sequence phylogenies that support the indicated direction of transmission, out of all deep-sequence phylogenies that support either direction of transmission.    
#' @param df data.table with the following columns  "IDCLU","ID1", "ID2", "TYPE","KEFF","LKL_MAX","POSTERIOR_SCORE" 
#' @param di data.table with meta-data to customize the plot with columns  "ID", node.shape, node.label, node.fill 
#' @param point.size size of the individual points
#' @param point.sizec.couple size of the outer ring around individuals in couples
#' @param edge.gap value to adjust start / end points of edges
#' @param edge.size multiplier by which the size of edges is shrunk/magnified
#' @param curvature curvature of directed edges  
#' @param arrow type of arrow to be plotted
#' @param curv.shift offset to place the label for directed edges
#' @param label.size size of label
#' @param node.label Text displayed on top of each node 
#' @param node.shape column name in di by which the shape of each node is drawn 
#' @param node.fill column name in di by which each node is coloured
#' @param node.shape.values named vector associating shapes to the values in the node.shape column
#' @param node.fill.values named vector associating colours to the values in the node.fill column
#' @param threshold.linked treshold value between 0 and 1. Edges with weight above this treshold are shown in black.
#' @return ggplot object
phsc.plot.most.likely.transmission.chain<- function(df, di, point.size=10, edge.gap=0.04, edge.size=0.4, curvature= -0.2, arrow=arrow(length=unit(0.04, "npc"), type="open"), curv.shift=0.08, label.size=3, node.label='ID', node.shape=NA_character_, node.fill=NA_character_, node.shape.values=NA_integer_, node.fill.values=NA_character_, threshold.linked=NA_real_, layout=NULL)
{
	#point.size=10; point.size.couple=14; edge.gap=0.04; edge.size=0.4; curvature= -0.2; arrow=arrow(length=unit(0.04, "npc"), type="open"); curv.shift=0.08; label.size=3
	#node.label='ID'; node.shape='IN_COUPLE'; node.fill='SEX'
	#node.fill.values=c('F'='hotpink2', 'M'='steelblue2')
	#node.shape.values=c('not in long-term\nrelationship'=18,'in long-term\nrelationship'=16)	
	if(is.na(node.label))
	{
		node.label<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.label, NA_character_)
	}
	if(is.na(node.shape))
	{
		node.shape<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.shape, 'NA')
	}
	if(is.na(node.fill))
	{
		node.fill<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.fill, 'NA')
	}
	if(any(is.na(node.fill.values)))
	{
		z						<- unique(di[[node.fill]])
		node.fill.values		<- heat.colors(length(z))
		names(node.fill.values)	<- z
	}
	if(any(is.na(node.shape.values)))
	{
		z						<- unique(di[[node.shape]])
		node.shape.values		<- seq_along(z)
		names(node.shape.values)<- z
	}
	setnames(di, c(node.label, node.shape, node.fill), c('NODE_LABEL','NODE_SHAPE','NODE_FILL'))
	tmp	<- c('NODE_LABEL','NODE_SHAPE','NODE_FILL')[which(c(node.label, node.shape, node.fill)=='ID')]
	if(length(tmp))
		set(di, NULL, 'ID', di[[tmp]])
	di	<- subset(di, select=c(ID, NODE_LABEL, NODE_SHAPE, NODE_FILL))
	if(is.null(layout))
	{
		layout	<- as.data.table(ggnet2(network(unique(subset(df, select=c(ID1,ID2))), directed=FALSE, matrix.type="edgelist"))$data[,c("label", "x", "y")])		
	}		
	if(any(grepl('label', colnames(layout))))
	{		
		setnames(layout, c('label','x','y'), c('ID1','ID1_X','ID1_Y'))
	}
	if(any(colnames(layout)=='ID'))
	{
		setnames(layout, c('ID','X','Y'), c('ID1','ID1_X','ID1_Y'))		
	}	
	df		<- merge(df, layout, by='ID1')
	setnames(layout, c('ID1','ID1_X','ID1_Y'), c('ID2','ID2_X','ID2_Y'))
	df		<- merge(df, layout, by='ID2')
	setnames(layout, c('ID2','ID2_X','ID2_Y'),  c('ID','X','Y'))	
	layout	<- merge(layout,di, by='ID')		
	df[, EDGETEXT_X:= (ID1_X+ID2_X)/2]
	df[, EDGETEXT_Y:= (ID1_Y+ID2_Y)/2]
	df[, EDGE_LABEL:= paste0('L',round(100*POSTERIOR_SCORE_LINKED,d=1),'%',' // ','D',round(100*pmax(POSTERIOR_SCORE_12,POSTERIOR_SCORE_21),d=1),'%') ]
	df[, EDGE_COL:= as.character(factor(POSTERIOR_SCORE_LINKED>threshold.linked, levels=c(TRUE,FALSE),labels=c('edge_col_2','edge_col_1')))]	
	#	for edges, move the start and end points on the line between X and Y
	#	define unit gradient
	df[, MX:= (ID2_X - ID1_X)]	
	df[, MY:= (ID2_Y - ID1_Y)]
	tmp		<- df[, sqrt(MX*MX+MY*MY)]
	set(df, NULL, 'MX', df[, MX/tmp])
	set(df, NULL, 'MY', df[, MY/tmp])	
	set(df, NULL, 'ID1_X', df[, ID1_X + MX*edge.gap])
	set(df, NULL, 'ID1_Y', df[, ID1_Y + MY*edge.gap])
	set(df, NULL, 'ID2_X', df[, ID2_X - MX*edge.gap])
	set(df, NULL, 'ID2_Y', df[, ID2_Y - MY*edge.gap])	
	p		<- ggplot() +			
			geom_point(data=layout, aes(x=X, y=Y, colour=NODE_FILL, pch=NODE_SHAPE), size=point.size) +
			geom_segment(data=subset(df, LINK_12==1), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*MX_KEFF_12, colour=EDGE_COL), arrow=arrow, lineend="butt") +			
			geom_segment(data=subset(df, LINK_21==1), aes(x=ID2_X, xend=ID1_X, y=ID2_Y, yend=ID1_Y, size=edge.size*MX_KEFF_21, colour=EDGE_COL), arrow=arrow, lineend="butt") +						
			scale_colour_manual(values=c(node.fill.values, 'edge_col_1'='grey80', 'edge_col_2'='grey40','NA'='grey50')) +
			scale_shape_manual(values=c(node.shape.values, 'NA'=16)) +
			scale_fill_manual(values=c(node.fill.values, 'NA'='grey50')) +
			scale_size_identity() +
			geom_text(data=df, aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=EDGE_LABEL), size=label.size) +
			geom_text(data=layout, aes(x=X, y=Y, label=NODE_LABEL)) +
			theme_void() +
			guides(colour='none', fill='none',size='none', pch='none')
	layout		<- subset(layout, select=c(ID,X,Y))
	setnames(layout, c('ID','X','Y'), c('label','x','y'))
	p$layout	<- layout
	p	
}

#' @export
#' @import data.table grid ggtree ggnet
#' @title Plot transmission network
#' @description This function plots a phylogenetic transmission network, showing three types of edges: 
#' two directed edges respectively in the 1->2 and 2->1 direction, and an undirected edge that represents phylogenetic support of close and adjacent individuals without evidence into the direction transmission.  
#' @param df data.table with the following columns  "IDCLU","ID1", "ID2", "TYPE","KEFF","LKL_MAX","POSTERIOR_SCORE" 
#' @param di data.table with meta-data to customize the plot with columns  "ID", node.shape, node.label, node.fill 
#' @param point.size size of the individual points
#' @param point.sizec.couple size of the outer ring around individuals in couples
#' @param edge.gap value to adjust start / end points of edges
#' @param edge.size multiplier by which the size of edges is shrunk/magnified
#' @param curvature curvature of directed edges  
#' @param arrow type of arrow to be plotted
#' @param curv.shift offset to place the label for directed edges
#' @param label.size size of label
#' @param node.label Text displayed on top of each node 
#' @param node.shape column name in di by which the shape of each node is drawn 
#' @param node.fill column name in di by which each node is coloured
#' @param node.shape.values named vector associating shapes to the values in the node.shape column
#' @param node.fill.values named vector associating colours to the values in the node.fill column
#' @param threshold.linked treshold value between 0 and 1. Edges with weight above this treshold are shown in black.
#' @return ggplot object
phsc.plot.transmission.network<- function(df, di, point.size=10, point.size.couple=point.size*1.4, edge.gap=0.04, edge.size=0.4, curvature= -0.2, arrow=arrow(length=unit(0.04, "npc"), type="open"), curv.shift=0.08, label.size=3, node.label='ID', node.shape=NA_character_, node.fill=NA_character_, node.shape.values=NA_integer_, node.fill.values=NA_character_, threshold.linked=NA_real_)
{	
	#point.size=10; point.size.couple=14; edge.gap=0.04; edge.size=0.4; curvature= -0.2; arrow=arrow(length=unit(0.04, "npc"), type="open"); curv.shift=0.08; label.size=3
	#node.label='ID'; threshold.linked=0.6; node.shape=NA_character_; node.fill=NA_character_; node.shape.values=NA_integer_; node.fill.values=NA_character_
	#node.shape='IN_COUPLE'; node.fill='SEX'
	#node.fill.values=c('F'='hotpink2', 'M'='steelblue2')
	#node.shape.values=c('not in long-term\nrelationship'=18,'in long-term\nrelationship'=16)
	if(is.na(node.label))
	{
		node.label<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.label, NA_character_)
	}
	if(is.na(node.shape))
	{
		node.shape<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.shape, 'NA')
	}
	if(is.na(node.fill))
	{
		node.fill<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
		set(di, NULL, node.fill, 'NA')
	}
	if(any(is.na(node.fill.values)))
	{
		z						<- unique(di[[node.fill]])
		node.fill.values		<- heat.colors(length(z))
		names(node.fill.values)	<- z
	}
	if(any(is.na(node.shape.values)))
	{
		z						<- unique(di[[node.shape]])
		node.shape.values		<- seq_along(z)
		names(node.shape.values)<- z
	}
	setnames(di, c(node.label, node.shape, node.fill), c('NODE_LABEL','NODE_SHAPE','NODE_FILL'))
	tmp	<- c('NODE_LABEL','NODE_SHAPE','NODE_FILL')[which(c(node.label, node.shape, node.fill)=='ID')]
	if(length(tmp))
		set(di, NULL, 'ID', di[[tmp]])
	di	<- subset(di, select=c(ID, NODE_LABEL, NODE_SHAPE, NODE_FILL))
	
	layout	<- as.data.table(ggnet2(network(unique(subset(df, select=c(ID1,ID2))), directed=FALSE, matrix.type="edgelist"))$data[,c("label", "x", "y")])
	setnames(layout, c('label','x','y'), c('ID1','ID1_X','ID1_Y'))
	df		<- merge(df, layout, by='ID1')
	setnames(layout, c('ID1','ID1_X','ID1_Y'), c('ID2','ID2_X','ID2_Y'))
	df		<- merge(df, layout, by='ID2')
	setnames(layout, c('ID2','ID2_X','ID2_Y'),  c('ID','X','Y'))	
	layout	<- merge(layout,di, by='ID')	
	
	df[, EDGETEXT_X:= (ID1_X+ID2_X)/2]
	df[, EDGETEXT_Y:= (ID1_Y+ID2_Y)/2]
	#
	#	calculate score for linked
	if(is.na(threshold.linked))
	{
		df	<- merge(df,df[, 	{
							z<- rep('edge_col_1', length(TYPE))
							z[which.max(POSTERIOR_SCORE)]	<- 'edge_col_2'
							list(EDGE_COL=z, TYPE=TYPE)	
						}, by=c('ID1','ID2')], by=c('ID1','ID2','TYPE'))		
	}
	if(!is.na(threshold.linked))
	{
		tmp	<- subset(df, TYPE!='not close/disconnected')[, list( EDGE_COL=as.character(factor(sum(POSTERIOR_SCORE)>=threshold.linked, levels=c(TRUE, FALSE), labels=c('edge_col_2','edge_col_1'))) ), by=c('ID1','ID2')]
		df	<- merge(df, tmp, by=c('ID1','ID2'))		
	}	
	#	for edges, move the start and end points on the line between X and Y
	#	define unit gradient
	df[, MX:= (ID2_X - ID1_X)]	
	df[, MY:= (ID2_Y - ID1_Y)]
	tmp		<- df[, sqrt(MX*MX+MY*MY)]
	set(df, NULL, 'MX', df[, MX/tmp])
	set(df, NULL, 'MY', df[, MY/tmp])	
	set(df, NULL, 'ID1_X', df[, ID1_X + MX*edge.gap])
	set(df, NULL, 'ID1_Y', df[, ID1_Y + MY*edge.gap])
	set(df, NULL, 'ID2_X', df[, ID2_X - MX*edge.gap])
	set(df, NULL, 'ID2_Y', df[, ID2_Y - MY*edge.gap])	
	#	label could just be move on the tangent vector to the line
	#	define unit tangent
	df[, TX:= -MY]
	df[, TY:= MX]
	tmp		<- df[, which(TYPE=='12')]
	set(df, tmp, 'EDGETEXT_X', df[tmp, EDGETEXT_X + TX*curv.shift])
	set(df, tmp, 'EDGETEXT_Y', df[tmp, EDGETEXT_Y + TY*curv.shift])
	tmp		<- df[, which(TYPE=='21')]
	set(df, tmp, 'EDGETEXT_X', df[tmp, EDGETEXT_X - TX*curv.shift])
	set(df, tmp, 'EDGETEXT_Y', df[tmp, EDGETEXT_Y - TY*curv.shift])
	#
	print(layout)
	print(node.fill.values)
	print(node.shape.values)
	
	p		<- ggplot() +			
			geom_point(data=layout, aes(x=X, y=Y, colour=NODE_FILL, pch=NODE_SHAPE), size=point.size) +
			geom_segment(data=subset(df, EDGE_COL=='edge_col_1' & TYPE=='ambiguous' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), lineend="butt") +
			geom_curve(data=subset(df, EDGE_COL=='edge_col_1' & TYPE=='12' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +
			geom_curve(data=subset(df, EDGE_COL=='edge_col_1' & TYPE=='21' & KEFF>0), aes(x=ID2_X, xend=ID1_X, y=ID2_Y, yend=ID1_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +
			geom_segment(data=subset(df, EDGE_COL=='edge_col_2' & TYPE=='ambiguous' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), lineend="butt") +
			geom_curve(data=subset(df, EDGE_COL=='edge_col_2' & TYPE=='12' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +
			geom_curve(data=subset(df, EDGE_COL=='edge_col_2' & TYPE=='21' & KEFF>0), aes(x=ID2_X, xend=ID1_X, y=ID2_Y, yend=ID1_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +									
			scale_colour_manual(values=c(node.fill.values, 'edge_col_1'='grey80', 'edge_col_2'='grey40','NA'='grey50')) +
			scale_shape_manual(values=c(node.shape.values, 'NA'=21)) +
			scale_fill_manual(values=c(node.fill.values, 'NA'='grey50')) +
			scale_size_identity() +
			geom_text(data=subset(df, TYPE!='not close/disconnected' & KEFF>0), aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=paste0(round(100*POSTERIOR_SCORE,d=1),'%')), size=label.size) +
			geom_text(data=layout, aes(x=X, y=Y, label=NODE_LABEL)) +
			theme_void() +
			guides(colour='none', fill='none',size='none', pch='none') 
	layout		<- subset(layout, select=c(ID,X,Y))
	setnames(layout, c('ID','X','Y'), c('label','x','y'))	
	p$layout	<- layout
	p
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
#' @param pdf.ntrees Number of trees per pdf.
#' @param pdf.title.size Size of pdf title in inches.
#' @param plot.file If not missing, the phylogenies will be printed to file.	
#' @return List of ggtree objects, ready for printing.
phsc.plot.phy.selected.pairs<- function(phs, dfs, id1, id2, plot.file=NA, pdf.h=50, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40)
{
	#	determine which phylogenies contain both individuals
	tmp		<- copy(dfs)
	tmp		<- merge(tmp, tmp[, {
						ph	<- phs[[ IDX ]]
						list(HAS_BOTH_IND= any(grepl(id1, attr(ph, "INDIVIDUAL"))) & any(grepl(id2, attr(ph, "INDIVIDUAL"))))  				
					}, by='IDX'], by='IDX')
	tmp		<- subset(tmp, HAS_BOTH_IND)
	phps	<- lapply(seq_len(nrow(tmp)), function(i){
				ph.title	<- NULL
				if('TITLE'%in%colnames(tmp))
					ph.title	<- tmp[i, TITLE]										
				ph			<- phs[[ tmp[i, IDX] ]]
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
						theme(legend.position="bottom", plot.title = element_text(size=pdf.title.size)) + 
						ggplot2::xlim(0, max(node.depth.edgelength(ph)[1:Ntip(ph)])*1.3) +
						labs(x='subst/site', title=ph.title)						
				p
			})
	#
	#	single page plot
	#		
	if(!is.na(plot.file))					
	{
		if(length(phps)<=pdf.ntrees)
		{
			cat('Plotting to file', plot.file,'...\n')
			pdf(file=plot.file, w=pdf.rw*length(phps), h=pdf.h)
			grid.newpage()
			pushViewport(viewport(layout=grid.layout(1, length(phps))))
			for(i in seq_along(phps))
				print(phps[[i]], vp = viewport(layout.pos.row=1, layout.pos.col=i))
			dev.off()
		}
		if(length(phps)>pdf.ntrees)
		{
			pi	<- data.table(IDX=seq_along(phps))
			pi[, PLOT:= ceiling(IDX/pdf.ntrees)]
			pi[, PLOT_IDX:= (IDX-1)%%pdf.ntrees+1]
			pi[,{
						cat('Plotting to file', gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file),'...\n')
						pdf(file=gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file), w=pdf.rw*pdf.ntrees, h=pdf.h)
						grid.newpage()
						pushViewport(viewport(layout=grid.layout(1, pdf.ntrees)))
						for(i in seq_along(IDX))
							print(phps[[IDX[i]]], vp = viewport(layout.pos.row=1, layout.pos.col=PLOT_IDX[i]))
						dev.off()
					}, by='PLOT']
		}
	}
	phps	
}

#' @export
#' @import data.table ggplot2 
#' @title Plot phyloscan
#' @description 
#' This function generates scan plots that summarize reconstructed viral phylogenetic relationships of two individuals.
#' Several pairs of individuals can be processed simultaneously. For each pair of individuals, the scan plot shows the 
#' phylogenetic distance (y-axis) and topological relationship (colours) between subgraphs from both individuals in each 
#' deep-sequence phylogeny across the genome. The genomic position on the x-axis indicates the start of each 250bp read alignment. 
#' @param rpw2 data.table containing the basic phyloscanner statistics for each genomic window. 
#' @param id.cols name of columns in rpw2 that identify the two individuals 
#' @param ylim limits of y-axis. 
#' @param cols.typet colour for each phylogenetic relationship type 
#' @return ggplot object
phsc.plot.phyloscan<- function(rpw2, id.cols=c('ID1','ID2'), ylim=NULL, cols.typet=NULL)
{		
	#	make manual plot to show intermingled
	if(is.null(ylim))
		ylim	<- c(1e-3,0.4)
	if(is.null(cols.typet))
		cols.typet			<- c(	"ancestral 1->2"=brewer.pal(11, 'PiYG')[2],
				'ancestral m->f'='steelblue2',
				"ancestral 2->1"=brewer.pal(11, 'PiYG')[10],
				'ancestral f->m'='hotpink2',
				"intermingled"=brewer.pal(11, 'PuOr')[4], 
				'sibling'=rev(brewer.pal(11, 'PuOr'))[c(3)], 
				"disconnected"=rev(brewer.pal(11, 'RdGy'))[4])		
	group		<- 'TYPE_BASIC'
	stopifnot(group%in%unique(rpw2$GROUP))
	stopifnot(id.cols%in%colnames(rpw2))
	#
	rpw3		<- subset(rpw2, GROUP==group)
	rpw3[, TYPE_TO:= 'disconnected']
	rpw3[, Y:=1e-3]
	set(rpw3, rpw3[, which(PATRISTIC_DISTANCE<1e-3)],'PATRISTIC_DISTANCE',1.1e-3)
	set(rpw3, rpw3[,which(grepl('chain 12', TYPE))], 'TYPE_TO', 'ancestral 1->2')
	set(rpw3, rpw3[,which(grepl('chain 21', TYPE))], 'TYPE_TO', 'ancestral 2->1')
	set(rpw3, rpw3[,which(grepl('chain mf', TYPE))], 'TYPE_TO', 'ancestral m->f')
	set(rpw3, rpw3[,which(grepl('chain fm', TYPE))], 'TYPE_TO', 'ancestral f->m')	
	set(rpw3, rpw3[,which(grepl('intermingled', TYPE))], 'TYPE_TO', 'intermingled')
	set(rpw3, rpw3[,which(grepl('sibling', TYPE))], 'TYPE_TO', 'sibling')	
	set(rpw3, NULL, 'TYPE_TO', rpw3[, factor(TYPE_TO, levels=c('ancestral 1->2','ancestral m->f','ancestral 2->1','ancestral f->m','intermingled','sibling','disconnected'))])	
	setnames(rpw3, id.cols, c('PRIVATECOL_ID1','PRIVATECOL_ID2'))	
	#	
	p <- ggplot(rpw3, aes(x=W_FROM)) +
			geom_hline(yintercept=0.025, colour='grey50') +
			geom_bar(aes(y=Y, fill=TYPE_TO), colour='transparent', stat='identity', width=25) +
			geom_point(aes(y=PATRISTIC_DISTANCE), size=1) +				
			labs(x='\ngenomic position\n(relative to HXB2)', y='subgraph distance\n(subst/site)\n',fill='topological subgraph\nrelationship') +
			scale_x_continuous(breaks=seq(0,1e4,500), minor_breaks=seq(0,1e4,100), limits=c(rpw3[, min(W_FROM)]-diff(rpw3[1:2, W_FROM]), rpw3[, max(W_FROM)]+diff(rpw3[1:2, W_FROM]))) +
			scale_y_log10(labels=percent, expand=c(0,0), breaks=c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25)) +
			coord_cartesian(ylim=ylim) +
			scale_fill_manual(values=cols.typet) +
			theme_bw() + 
			theme(legend.position='bottom', panel.spacing = unit(1, "lines")) +
			facet_grid(PRIVATECOL_ID1+PRIVATECOL_ID2~.)
	p	
}	

#' @export
#' @import data.table grid ggtree colorspace
#' @title Plot short read phylogenies and highlight individuals
#' @description This function plots short read phylogenies and highlights the clades of two individuals in red and blue.  
#' @param phs List of trees in ape format
#' @param dfs data.table with mandatory column 'IDX' and optional column 'TITLE'. IDX is the index of all phylogenies in 'phs' that are to be plotted. TITLE is a title for each sub-plot, for example specifying a window.
#' @param ids Vector of regular expressions that identify individuals to be highlighted in colour.
#' @param plot.cols Vector of colours for each individual
#' @param group.redo Logical, indicating if the colour groups should be recalculated from ids.
#' @param drop.blacklisted Logical, indicating if all blacklisted taxa should be dropped prior to plotting.  
#' @param pdf.h	Height of the pdf file in inches.
#' @param pdf.rw Relative width of the pdf file, internally multiplied by the number of phylogenies to give the total width in inches.
#' @param pdf.ntrees Number of trees per pdf.
#' @param pdf.title.size Size of pdf title in inches.
#' @param plot.file If not missing, the phylogenies will be printed to file.	
#' @return List of ggtree objects, ready for printing.
phsc.plot.phycollapsed.selected.individuals<- function(phs, dfs, ids, plot.cols=rainbow(length(ids)), plot.cols.background=function(n){ rainbow_hcl(n, start = 60, end = 240, c=30, l=80) }, drop.blacklisted=TRUE, tip.regex='(.*)_read_([0-9]+)_count_([0-9]+)', plot.file=NA, drop.less.than.n.ids=2, pdf.h=50, pdf.rw=10, pdf.ntrees=20, pdf.title.size=40)
{	
	#	determine which phylogenies contain at least one of the requested individuals
	tmp		<- copy(dfs)
	tmp		<- merge(tmp, tmp[, {
						ph	<- phs[[ IDX ]]
						z	<- ph$tip.label
						if(drop.blacklisted)
							z	<- as.character(attr(ph, "INDIVIDUAL"))
						list(HAS_N_IND=length(which(sapply(ids, function(id) any(grepl(id, z)) ))))						  				
					}, by='IDX'], by='IDX')
	tmp		<- subset(tmp, HAS_N_IND>=drop.less.than.n.ids)
	#	define colours so that they don t change across windows
	cols	<- tmp[, {
				ph	<- phs[[ IDX ]]
				phb			<- data.table(	TAXA=ph$tip.label,
						IDX=seq_along(ph$tip.label), 											 
						ID= gsub(tip.regex,'\\1',ph$tip.label),
						COUNT=sub(tip.regex,'\\3',ph$tip.label)		)
				set(phb, phb[, which(ID==COUNT)], 'COUNT', NA_character_)
				set(phb, NULL, 'COUNT', phb[, as.integer(COUNT)])
				set(phb, phb[, which(is.na(COUNT))], 'ID','REFERENCE')	
				list(ID= phb[, unique(ID)])				
			}, by='IDX']
	cols	<- unique(subset(cols, select=ID))
	cols[,COL:= plot.cols.background(nrow(cols))]	
	set(cols, cols[,which(grepl('REFERENCE',ID))], 'COL', 'grey50')
	for(i in seq_along(ids))
		set(cols, cols[, which(grepl(ids[i], ID))], 'COL', plot.cols[i])
	#	
	phps	<- lapply(seq_len(nrow(tmp)), function(i){
				#cat(i,'\n')
				ph.title	<- NULL
				if('TITLE'%in%colnames(tmp))
					ph.title	<- tmp[i, TITLE]										
				ph			<- phs[[ tmp[i, IDX] ]]
				ph			<- phsc.phy.collapse.monophyletic.clades(ph, drop.blacklisted=TRUE, tip.regex=tip.regex)
				#
				#	define IDs, find COUNTS, find references
				phb			<- data.table(	TAXA=ph$tip.label,
						IDX=seq_along(ph$tip.label), 											 
						ID= gsub(tip.regex,'\\1',ph$tip.label),
						COUNT=sub(tip.regex,'\\3',ph$tip.label)		)
				set(phb, phb[, which(ID==COUNT)], 'COUNT', NA_character_)
				set(phb, NULL, 'COUNT', phb[, as.integer(COUNT)])
				set(phb, phb[, which(is.na(COUNT))], 'ID','REFERENCE')	
				set(phb, phb[, which(is.na(COUNT))], 'COUNT', 1L)
				#	define cols for this tree
				tmp					<- data.table(IDX= seq_along(attr(ph, "INDIVIDUAL")), ID= attr(ph, "INDIVIDUAL"))
				set(tmp, tmp[, which(is.na(ID))], 'ID', 'REFERENCE')
				tmp					<- merge(tmp, cols, by='ID')
				attr(ph, 'COLOUR')	<- tmp[order(IDX),][,COL]	
				#						
				p 			<- ggtree(ph, aes(color=I(COLOUR))) %<+% phb +
						geom_tippoint(aes(size=COUNT), shape=18) +
						#geom_point2(shape=16, aes(size=NODE_S, subset=NODE_S)) +
						scale_fill_hue(na.value="black") +								
						theme(legend.position="none") +
						geom_tiplab(align=T, aes(col=I(COLOUR)), size=1, linetype=NA, linesize=NA) +
						guides(size='none') +
						theme_tree2() +
						theme(	axis.line.x=element_line(),
								panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3),
								legend.position="bottom",								
								plot.title = element_text(hjust = 0.5, size=pdf.title.size)) + 
						ggplot2::xlim(0, max(node.depth.edgelength(ph)[1:Ntip(ph)])*1.15) +
						labs(x='subst/site', title=ph.title)
				#ggsave(plot.file, w=10, h=200, limitsize = FALSE)
				p
			})
	#
	#	single page plot
	#		
	if(!is.na(plot.file))					
	{
		if(length(phps)<=pdf.ntrees)
		{
			cat('Plotting to file', plot.file,'...\n')
			pdf(file=plot.file, w=pdf.rw*length(phps), h=pdf.h)
			grid.newpage()
			pushViewport(viewport(layout=grid.layout(1, length(phps))))
			for(i in seq_along(phps))
				print(phps[[i]], vp = viewport(layout.pos.row=1, layout.pos.col=i))
			dev.off()
		}
		if(length(phps)>pdf.ntrees)
		{
			pi	<- data.table(IDX=seq_along(phps))
			pi[, PLOT:= ceiling(IDX/pdf.ntrees)]
			pi[, PLOT_IDX:= (IDX-1)%%pdf.ntrees+1]
			pi[,{
						cat('Plotting to file', gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file),'...\n')
						pdf(file=gsub('\\.pdf',paste('_plot',PLOT,'\\.pdf',sep=''),plot.file), w=pdf.rw*pdf.ntrees, h=pdf.h)
						grid.newpage()
						pushViewport(viewport(layout=grid.layout(1, pdf.ntrees)))
						for(i in seq_along(IDX))
							print(phps[[IDX[i]]], vp = viewport(layout.pos.row=1, layout.pos.col=PLOT_IDX[i]))
						dev.off()
					}, by='PLOT']
		}
	}
	phps	
}
