% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloscan.fun.private.R
\name{phsc.get.most.likely.transmission.chains.RBGLedmonds}
\alias{phsc.get.most.likely.transmission.chains.RBGLedmonds}
\title{Construct maximum probability transmission network}
\usage{
phsc.get.most.likely.transmission.chains.RBGLedmonds(rtnn, verbose = 0)
}
\arguments{
\item{rtn}{data.table with network scores for all individuals that could form a network. Must contain columns 'ID1','ID2','IDCLU','GROUP','TYPE','POSTERIOR_SCORE','KEFF'.}
}
\value{
new data.table with added columns LINK_12 LINK_21 (either 1 or 0), and MX_PROB_12 MX_PROB_21 (associated posterior probabilities)
}
\description{
This function reconstructs a maximum probility transmission 
network from the scores associated with directed and undirected edges.
The algorithm starts by keeping the edge with highest score.
It then removes the competitor in the opposite direction, and any conflicting edges that would result in indegrees larger than one.
By construction, all removed edges have lower probability.
The algorithm proceeds until all edges have been processed.
}
\keyword{internal}
