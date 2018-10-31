# Phyloscanner.R.utilities

# Utility functions for deep-sequence phylogenetic analyses with *phyloscanner*

## Overview
*Phyloscanner.R.utilities* is a set of functions used to help apply [the phyloscanner software](https://github.com/BDI-pathogens/phyloscanner) to large datasets of a population-based sample, and reconstruct transmission networks from multiple *phyloscanner* runs. Please go to the [*phyloscanner* github page](https://github.com/BDI-pathogens/phyloscanner) for access to the main software package.

The software package comprises
1. R functions for setting up and running a large number of *phyloscanner* analyses on batches of deep-sequence data in parallel. 
2. R functions for collecting *phyloscanner* output of multiple batches, and for reconstructing transmission networks.
3. R functions to visualise reconstructed transmission networks and primary *phyloscanner* output. 

## System requirements
Both *phyloscanner* and *Phyloscanner.R.utilities* are supported on *Linux* and *MacOS*.
1. *phyloscanner* builds on standard tools for deep-sequence data analysis and phylogeny reconstruction, such as *samtools* and *RAxML*. Please install these first; [instructions are available here for *Linux* and *MacOS*](https://github.com/BDI-pathogens/phyloscanner/blob/master/InfoAndInputs/InstallationNotesForMakingTrees.sh). If you have installation issues, [please report it here and we will get back to you](https://github.com/BDI-pathogens/phyloscanner/issues). This step may take up to 60 minutes. 
2. *phyloscanner* and *Phyloscanner.R.utilities* depend on several *R* packages. We find it easiest to install them as follows. First run the following `install_github` command, and then install any packages that could not be installed manually (these packages not on CRAN and need to be installed from Bioconductor/ github) That is:    
```r
devtools:::install_github("olli0601/Phyloscanner.R.utilities", dependencies=TRUE, build_vignettes=FALSE)
``` 
On a fresh R build, installation will fail with
```text
ERROR: dependencies ‘ggtree’, ‘ggnet’, ‘Rsamtools’, ‘RBGL’ are not available
```
You will then need
```r
source("https://bioconductor.org/biocLite.R")
biocLite(c("Rsamtools","RBGL","ggtree"),dependencies=TRUE, build_vignettes=FALSE)
devtools:::install_github("briatte/ggnet", dependencies=TRUE)
``` 
This step may take up to 60 minutes. If you have issues with installation/running of *Phyloscanner.R.utilities*, [please report it here and we will get back to you](https://github.com/olli0601/Phyloscanner.R.utilities/issues).

## Installation
1. [Download *phyloscanner* version 1.1.2](https://github.com/olli0601/Phyloscanner.R.utilities/tree/master/misc/phyloscanner_v1.1.2.tar.gz) and unzip to a directory of your choice.
2. Install *Phyloscanner.R.utilities* in R:
```r
devtools:::install_github("olli0601/Phyloscanner.R.utilities", dependencies=TRUE, build_vignettes=FALSE)
require(Phyloscanner.R.utilities)
``` 
This step takes less than 5 minutes. If you have issues with installation/running of *Phyloscanner.R.utilities*, [please report it here and we will get back to you](https://github.com/olli0601/Phyloscanner.R.utilities/issues). 

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


### [Running *phyloscanner* in high-throughput](articles/Rakai.02.run_phyloscanner.html)
[This tutorial](articles/Rakai.02.run_phyloscanner.html) starts with *phyloscanner* analyses of deep-sequence phylogenetic trees for each potential transmission network identified on this population-based sample of 2,652 infected individuals. The main objective is to illustrate how large numbers of phyloscanner runs can be generated and run in parallel with the utility functions in this software package, without too much computational overhead.

Expected runtime: *phyloscanner* analysis of each batch of transmission networks takes about 2 hours, and the 345 batches can be processed in parallel.

### [Reconstructing transmission networks](articles/Rakai.03.reconstruct_transmission_networks.html)
[This tutorial](articles/Rakai.03.reconstruct_transmission_networks.html) uses the generated *phyloscanner* output to reconstruct partially sampled HIV-1 transmission networks, of individuals who are in the population-based sample. We will illustrate a range of analytical functions to create these networks, and plotting functions to visualise phylogenetic inferences.  

Expected runtime: 10 minutes.

### [Inferring the direction of transmission](articles/Rakai.04.direction_of_transmission.html)
[This tutorial](articles/Rakai.04.direction_of_transmission.html) builds on the reconstructed HIV-1 transmission networks, and identifies phylogenetically highly supported source-recipient pairs. We will validate our phylogenetic inference into the direction of transmission against available epidemiologic and clinical data. 

Expected runtime: 10 minutes.

## More stuff
1. [Tutorial to generate read alignments for stage 1 analysis](articles/Stage1.create.read.alignments.html) .

## Testing
The software has been run and tested on MacOS X 10.10, 10.11, 10.12 (Yosemite, El Capitan, Sierra) and CentOS Linux release 7.3.1611; and RAxML (8.2.9), 
MAFFT (7.212), samtools (1.2), anaconda (2.3.0), pysam (0.8.1), R (3.1), argparse (1.0.4), ape (4.1), data.table (1.10.5), dplyr (0.7.4), dtplyr (0.0.1), ff (2.2-13), ggplot2 (2.2.1.9000), grid (3.3.3), gridExtra (2.2.1), gtable (0.2.0), ggtree (1.6.9), kimisc (0.3), pegas (0.9), phangorn (2.1.1), phytools (0.5-64), prodlim (1.5.7), RColorBrewer (1.1-2), reshape (0.8.6), reshape2 (1.4.3), scales (0.5.0.9000), colorspace (1.3-2), devtools (1.13.4), ggnet (0.1.0), igraph (1.0.1), Rsamtools (1.26.1), RBGL (1.55.1), sna (2.4), knitr (1.17), rmarkdown (1.7).