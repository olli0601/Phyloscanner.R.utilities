# Phyloscanner.R.utilities

# Utility functions for deep-sequence phylogenetic analyses with *phyloscanner*

## Overview
*Phyloscanner.R.utilities* is a set of functions used to help apply [the phyloscanner software](https://github.com/BDI-pathogens/phyloscanner) to large datasets of a population-based sample, and reconstruct transmission networks from multiple *phyloscanner* runs. Please go to the [*phyloscanner* github page](https://github.com/BDI-pathogens/phyloscanner) for access to the main software package.

The software package comprises
1. R functions for setting up and running a large number of *phyloscanner* analyses on batches of deep-sequence data in parallel. 
2. R functions for collecting *phyloscanner* output of multiple batches, and for reconstructing transmission networks.
3. R functions to visualise reconstructed transmission networks and primary *phyloscanner* output. 

## Installation
Both *phyloscanner* and *Phyloscanner.R.utilities* are supported on *Linux* and *MacOS*. 
1. Instructions for installing *phyloscanner* on either system [are available here](https://github.com/BDI-pathogens/phyloscanner/blob/master/InfoAndInputs/InstallationNotesForMakingTrees.sh). In addition, *phyloscanner* requires *R* version >= 3.1, and the following *R* packages: argparse, ape, data.table, dplyr, dtplyr, ff, ggplot2, grid, gridExtra, gtable, ggtree, kimisc, pegas, phangorn, phytools, prodlim, RColorBrewer, reshape, reshape2, scales. These should be installed as part of the `devtools:::install_github` command below. If you have issues with installation of *phyloscanner*, [please report it here and we will get back to you](https://github.com/BDI-pathogens/phyloscanner/issues). 
2. *Phyloscanner.R.utilities* requires *R* version >= 3.1, and the following additional *R* packages: colorspace, devtools, ggnet, igraph, Rsamtools, RBGL, sna. These should be installed as part of the `devtools:::install_github` command below. If you have issues with installation/running of *Phyloscanner.R.utilities*, [please report it here and we will get back to you](https://github.com/olli0601/Phyloscanner.R.utilities/issues).
```r
devtools:::install_github("olli0601/Phyloscanner.R.utilities", dependencies=TRUE, build_vignettes=FALSE)
require(Phyloscanner.R.utilities)
```

## General protocol for analyses at the population-level
It is computationally challenging to reconstruct viral trees from 
deep-sequence reads of hundreds or thousands of individuals. To
address this challenge, we generally proceed in two stages. 

**In a first stage**, the large population-based sample is divided into groups of 50
to 75 individuals, and then *phyloscanner* is run on all possible pairs of groups
to generate read alignments and deep-sequence phylogenies. From these trees, potentially phylogenetically close pairs are identifed and from those, networks of
individuals that were connected through at least one common, phylogenetically close
individual. 

**In a second stage**, *phyloscanner* analyses are repeated on all individuals that together form a potential transmission network. In addition to reads from all individuals in a potential transmission network, each analysis also includes reads from closely related individuals that act as control sequences. This step allows us to resolve the ordering of
transmission events within transmission networks, and to confirm potential
transmission pairs within a network. 

**Both stages 1 and 2 consist of very similar steps.** In either stage, phylogenetic trees are first generated from deep-sequence reads of batches of individuals (either from two randomly paired groups of individuals in stage 1, or phylogenetically closely connected individuals in a potential transmission network in stage 2). Thereafter, phylogenetic relationships between these individuals are estimated from the reconstructed deep-sequence trees. 

**The similarity of the two stages allows us to use the same high-throughput R scripts for analysis**  

## Tutorials for *phyloscanner* analyses of large population-based samples
To demonstrate analysis of a deep-sequence reads from of large population-based sample of individuals, we here provide several tutorials on data from Rakai District, Uganda. [The data are described here.](articles/Rakai.01.data_description.html)


### [Running *phyloscanner* in high-throughput](articles/Rakai.01.run_phyloscanner.html)
[This tutorial](articles/Rakai.01.run_phyloscanner.html) starts with *phyloscanner* analyses of deep-sequence phylogenetic trees for each potential transmission network identified on this population-based sample of 2,652 infected individuals. The main objective is to illustrate how large numbers of phyloscanner runs can be generated and run in parallel with the utility functions in this software package, without too much computational overhead.

Expected runtime: *phyloscanner* analysis of each batch of transmission networks takes about 2 hours, and the 345 batches can be processed in parallel.

### [Reconstructing transmission networks](articles/Rakai.02.reconstruct_transmission_networks.html)
[This tutorial](articles/Rakai.02.reconstruct_transmission_networks.html) uses the generated *phyloscanner* output to reconstruct partially sampled HIV-1 transmission networks, of individuals who are in the population-based sample. We will illustrate a range of analytical functions to create these networks, and plotting functions to visualise phylogenetic inferences.  

Expected runtime: 10 minutes.

## Testing
The software has been run and tested on MacOS X 10.10, 10.11, 10.12 (Yosemite, El Capitan, Sierra) and CentOS Linux release 7.3.1611.  

