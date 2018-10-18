# Reconstructing transmission networks

## Introduction
[In the previous tutorial](Rakai.01.run_phyloscanner.md), we performed *phyloscanner* analyses to reconstruct phylogenetic relationships between individuals that together form a potential transmission network. The analysis proceeded in batches, in which each transmission network was given a batch number, and small networks were combined into a single batch. Each batch was identified by the prefix `ptyrX_` where `X` is the batch number. 

**In this tutorial, we will take a look at the output.** We will:
1. Collect pairs of individuals between whom phylogenetic linkage cannot be excluded based on the distance and topological relationship of viral reads from both individuals. 
2. Reconstruct transmission networks from these pairs. Unlike typical phylogenetic clusters, the transmission networks that we reconstruct from deep-sequence data contain information into the direction of transmission. 
3. Calculate the proportion of male-female, male-male, and female-female links in the most likely transmission chains. 

## Setting up the analysis
[In the previous tutorial](Rakai.01.run_phyloscanner.md), we generated a large number of files, including `ptyrX_pairwise_relationships.rda`. We will need only these files to reconstruct HIV-1 transmission networks. 

To be safe, it might be a good idea to copy them into a new directory. I copied them to the analysis directory:
```r
indir <- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiPopSample_phyloscanner_analysis' 
```
 
## Find all pairs of individuals between whom linkage cannot be excluded, then reconstruct transmission networks
The following code snippet processes all pairwise phylogenetic relationships that were reconstructed with *phyloscanner* [in the previous tutorial](Rakai.01.run_phyloscanner.md), and returns a data.table that contains all pairs of individuals between whom phylogenetic linkage is not excluded based on distance and adjacency. There are three input arguments: 
1. `batch.regex` identifies the batch number from the file names of *phyloscanner* output; 
2. `neff.cut`  specifies the minimum number of deep-sequence phylogenies with sufficient reads from two individuals in order to make any phylogenetic inferences (default is 3); 
3. `conf.cut` specifies the proportion of deep-sequence phylogenies with distant/disconnected subgraphs above which pairs are considered phylogenetically unlinked (default is 60%).

At this stage, it is also easy to add any further individual-level meta-data to the analysis output. Simply specify a data.table to the `dmeta` input variable. 

For the analysis of the Rakai population-based sample:
```r
infile <- "~/Dropbox (SPH Imperial College)/2017_phyloscanner_validation/Supp_Data/Data_Set_S2.csv"
dmeta <- as.data.table(read.csv(infile, stringsAsFactors=FALSE))
	
indir <- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiPopSample_phyloscanner_analysis'
outfile <- '~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiPopSample_phyloscanner_analysis/todi_pairs_171122_cl25_d50_prior23_min30.rda'
tmp <- phsc.find.linked.pairs(indir, batch.regex='^ptyr([0-9]+)_.*', conf.cut=0.6, neff.cut=3, verbose=TRUE, dmeta=dmeta)
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
rtp <- copy(tmp$linked.pairs)
rplkl <- copy(tmp$relationship.counts)
rpw <- copy(tmp$windows)
save(rtp, rplkl, rpw, file=outfile)	 
```

<p align="center"><img src="Rakai.02.reconstruct_transmission_networks.rpw.png" alt="Output of phyloscanner statistics for each window."/></p>

<p align="center"><img src="Rakai.02.reconstruct_transmission_networks.rplkl.png" alt="Phylogenetic relationship counts."/></p>

<p align="center"><img src="Rakai.02.reconstruct_transmission_networks.rtp.png" alt="Pairs between whom transmission cannot be excluded phylogenetically."/></p>


plot phyloscans for some pairs
get networks
plot networks
plot ML chains


