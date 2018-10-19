# Reconstructing transmission networks

## Introduction
[In the previous tutorial](Rakai.01.run_phyloscanner.md), we performed *phyloscanner* analyses to reconstruct phylogenetic relationships between individuals that together form a potential transmission network. The analysis proceeded in batches, in which each transmission network was given a batch number, and small networks were combined into a single batch. Each batch was identified by the prefix `ptyrX_` where `X` is the batch number. 

**In this tutorial, we will take a look at the output.** We will:
1. Collect pairs of individuals between whom phylogenetic linkage cannot be excluded based on the distance and topological relationship of viral reads from both individuals. 
2. Reconstruct transmission networks from these pairs. Unlike typical phylogenetic clusters, the transmission networks that we reconstruct from deep-sequence data contain information into the direction of transmission. 
3. Calculate the proportion of male-female, male-male, and female-female links in the most likely transmission chains. 

## Setting up the analysis
[In the previous tutorial](Rakai.01.run_phyloscanner.md), we generated a large number of files, including `ptyrX_pairwise_relationships.rda`. We will need only these files to reconstruct HIV-1 transmission networks. 

To be safe, it might be a good idea to copy them into a new directory. I copied them to the following analysis directory, and all further output will be stored in the same directory:
```r
HOME <- '/Users/Oliver/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA'
indir <- file.path(HOME, 'RakaiPopSample_phyloscanner_analysis')
outdir <- indir
outfile.base <- file.path(outdir,'phsc_analysis_of_dataset_S1')
neff.cut <- 3
conf.cut <- 0.6		 
```
 
## Find all pairs of individuals between whom linkage cannot be excluded
The following code snippet processes all pairwise phylogenetic relationships that were reconstructed with *phyloscanner* [in the previous tutorial](Rakai.01.run_phyloscanner.md), and returns a data.table that contains all pairs of individuals between whom phylogenetic linkage is not excluded based on distance and adjacency. There are three input arguments: 
1. `batch.regex` identifies the batch number from the file names of *phyloscanner* output; 
2. `neff.cut`  specifies the minimum number of deep-sequence phylogenies with sufficient reads from two individuals in order to make any phylogenetic inferences (default is 3); 
3. `conf.cut` specifies the proportion of deep-sequence phylogenies with distant/disconnected subgraphs above which pairs are considered phylogenetically unlinked (default is 60%).

At this stage, it is also easy to add any further individual-level meta-data to the analysis output. Simply specify a data.table to the `dmeta` input variable. 

For the analysis of the Rakai population-based sample:
```r
infile <- "~/Dropbox (SPH Imperial College)/2017_phyloscanner_validation/Supp_Data/Data_Set_S2.csv"
dmeta <- as.data.table(read.csv(infile, stringsAsFactors=FALSE))
tmp <- phsc.find.linked.pairs(indir, batch.regex='^ptyr([0-9]+)_.*', conf.cut=0.6, neff.cut=3, verbose=TRUE, dmeta=dmeta)
rtp <- copy(tmp$linked.pairs)
rpw <- copy(tmp$windows)
rplkl <- copy(tmp$relationship.counts)
```

This will produce the following output (when `verbose=TRUE`):
```text
Found phylogenetic relationship files, n= 345
Processing...
Found (potentially duplicate) pairs between whom linkage is not excluded phylogenetically, n= 1705
Collect phylogenetic relationship counts for each pair...
Collect basic phyloscanner statistics (distance, adjacency, paths between subgraphs) for each pair...
Re-arrange pairs so that ID1<ID2...
If pairs are in several batches, select batch with most deep-sequence phylogenies...
Left with pairs between whom linkage is not excluded phylogenetically, n= 1326
Add meta-data...
Done. Found pairs, n= 1251 . Found relationship counts, n= 60570 . Found phyloscanner statistics, n= 1076768 .
```

Let us save the output, and then take a closer look:
```r
save(rtp, rplkl, rpw, file=paste0(outfile.base,'_allpairs.rda'))	 
```

## *phyloscanner* statistics for each window
The data.table `rpw` describes for each pair (`ID1`, `ID2`) the basic *phyloscanner* statistics (patristic distance, adjacency, contiguity, paths from subgraphs of individual 1 to subgraphs of individual 2, and vice versa) across the genome, plus any meta-data that was potentially added: 
 
<p align="center"><img src="figures/Rakai.02.reconstruct_transmission_networks.rpw.png" alt="Output of phyloscanner statistics for each window."/></p>

In addition, `rpw` describes how the phylogenetic relationship of the two individuals is classified in each genomic window. There are several classifications: 
1. The most important one is probably `GROUP=='TYPE_ADJ_NETWORK_SCORES'`. In each genomic window, pairs are classified as either '12' (subgraphs of individual 1 are close, adjacent and ancestral to those from individual 2), '21' (subgraphs of individual 2 are close, adjacent and ancestral to those from individual 1), 'ambiguous' (subgraphs are close, adjacent and either intermingled or sibling), or 'not close/disconnected' (subgraphs are either not close, or not adjacent).
2. There are two more important classification types. The relationship types that support phylogenetic linkage are '12', '21', and 'ambiguous'. In the classification `GROUP=='TYPE_CHAIN_TODI'` these three states are already summarised to state 'chain'. 
3. Finally, the classification `GROUP=='TYPE_ADJ_DIR_TODI2'` only considers genomic windows with phylogenetic support into the direction of transmission. The possible states are either '12', '21', and all other genomic windows are assigned NA.     

We can visualise the phylogenetic relationships between two individuals across the genome with the function `phsc.plot.phyloscan`:	
```r
rpw2 <- subset(rpw, ID1=='RkA04565F' & ID2=='RkA05315F')		
p	<- phsc.plot.phyloscan(rpw2)
p
ggsave(file=paste0(outfile.base,'_phyloscan_RkA04565F_RkA05315F.png'), width=6, height=2.8, units='in', dpi=400)
```	
<p align="center"><img src="figures/phsc_analysis_of_dataset_S1_phyloscan_RkA04565F_RkA05315F.png" alt="Output of phyloscanner statistics for each window."/></p>

The plot shows the phylogenetic relationship between the two females *RkA04565F* *RkA05315F* across 55 overlapping deep-sequence trees on reads that cover the *gag* gene. The start position of each 250bp genomic window is plotted on the x-axis. The genetic distance between the subgraphs of both females are shown on the y-axis: the two females have nearly identical virus in all trees. The topological relationship of the subgraphs of the two females is indicated in colours: the two females have intermingled virus across nearly all deep-sequence phylogenies. Since HIV is extremely rarely transmitter between women, the important conclusion is that even when virus is heavily intermingled and nearly identical, it is not possible to prove that transmission occurred between the corresponding two individuals.    

## Relationship counts for each pair
The data.table `rplkl` is a summary of the information in `rpw` for each genomic window, and gives for each pair the counts of how often certain phylogenetic relationships were seen across the genome. 
1. For each pair, column `N` gives the total number of deep-sequence phylogenies in which both individuals had sufficient reads for phylogenetic analysis.
2. Column `TYPE` gives a particular phylogenetic relationship type, for example '12', and column `K` gives the number of deep-sequence phylogenies in whom the subgraphs of both individuals were of that type.
3. Colums `NEFF` and `KEFF` are similar to `N` and `K`, but adjust for potential overlap in read alignments.     
Here is a screenshot of data.table `rplkl` for the *phyloscanner* analysis of the Rakai population-based sample:   
<p align="center"><img src="figures/Rakai.02.reconstruct_transmission_networks.rplkl.png" alt="Phylogenetic relationship counts."/></p>


## Finally reconstruct transmission networks
plot phyloscans for some pairs
get networks
plot networks
plot ML chains


