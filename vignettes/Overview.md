# Utility functions for deep-sequence phylogenetic analyses at the population level with *phyloscanner*

## Overview
*Phyloscanner.R.utilities* is a set of functions used to help apply [the phyloscanner software](https://github.com/BDI-pathogens/phyloscanner) to large datasets of a population-based sample, and reconstruct transmission networks from multiple *phyloscanner* runs. 

The software package comprises
1. R functions for setting up and running a large number of *phyloscanner* analyses on batches of deep-sequence data in parallel. 
2. R functions for collecting *phyloscanner* output of multiple batches, and for reconstructing transmission networks.
3. R functions to visualise reconstructed transmission networks and primary *phyloscanner* output. 

## Installation
Both *phyloscanner* and *Phyloscanner.R.utilities* are supported on *Linux* and *MacOS*. 
1. Instructions for installing *phyloscanner* on either system [are available here](https://github.com/BDI-pathogens/phyloscanner/blob/master/InfoAndInputs/InstallationNotesForMakingTrees.sh). If you have issues with installation of *phyloscanner*, [please report it here and we will get back to you](https://github.com/BDI-pathogens/phyloscanner/issues). 
2. *Phyloscanner.R.utilities* requires *R* version >= 3.1, and the following *R* packages: ape, argparse, phytools, phangorn, reshape2, data.table, RColorBrewer, colorspace, devtools, grid, gridExtra, ggplot2, ggtree, zoo. Please install these first. When done, the software package can be installed via
```r
devtools:::install_github("olli0601/Phyloscanner.R.utilities")
require(Phyloscanner.R.utilities)
``` 
If you have issues with installation/running of *Phyloscanner.R.utilities*, [please report it here and we will get back to you](https://github.com/olli0601/Phyloscanner.R.utilities/issues). 

## Testing
The software has been run and tested on MacOS X 10.10, 10.11, 10.12 (Yosemite, El Capitan, Sierra) and CentOS Linux release 7.3.1611.  

## Vignettes for *phyloscanner* analyses of large population-based samples
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

**Both stages 1 and 2 consist of very similar steps.** In either stage, phylogenetic trees
are first generated from deep-sequence reads of batches of individuals (either from two randomly paired groups of individuals in stage 1, or phylogenetically closely connected individuals in a potential transmission network in stage 2). Thereafter, phylogenetic relationships between these individuals are estimated from the reconstructed deep-sequence trees. 

**This tutorial describes the steps to infer phylogenetic
relationships from a large number of deep-sequence phylogenies of individuals in
the same potential transmission network (stage 2).**  

The code below assumes that *phyloscanner_make_trees.py* was already run to
generate read alignments and deep-sequence phylogenies for individuals in the
same potential transmission network (see here); and that output from this step
is available in the following file structure.
1. Each analysis of a potential transmission network is identified with the
   prefix `ptyrX_` where `X` is an integer. Note that, to minimise computations,
   individuals in small transmission networks can be grouped into a single
   *phyloscanner* analysis, which we call a *batch* in this tutorial. 
2. Three files are available for each batch. First, a text file listing all
   individuals in the batch, in file `ptyrX_patients.txt`. 
3. Second, the read alignments generated with *phyloscanner*, zipped into file
   `ptyrX_trees_fasta.zip`. 
4. Third, the deep-sequence phylogenies generated with *phyloscanner*, zipped
   into file `ptyrX_trees_newick.zip`. 

Data Set S1 contains these files for 345 batches of individuals of the Rakai
population-based sample, that comprise 1,426 potential transmission pairs and
closely related control sequences. 

