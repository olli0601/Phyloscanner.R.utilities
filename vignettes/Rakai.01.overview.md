# Utility functions for deep-sequence phylogenetic analyses at the population level with *phyloscanner*

## Overview
*Phyloscanner.R.utilities* is a set of functions used to help apply [the phyloscanner software](https://github.com/BDI-pathogens/phyloscanner) to large datasets of a population-based sample, and reconstruct transmission networks from multiple *phyloscanner* runs. 

The utility software package comprises
1. R functions for setting up and running a large number of *phyloscanner* analyses on batches of deep-sequence data in parallel. 
2. R functions for collecting *phyloscanner* output of multiple batches, and for reconstructing transmission networks.
3. R functions to visualise reconstructed transmission networks and primary *phyloscanner* output. 

## Installation
Both *phyloscanner* and *Phyloscanner.R.utilities* are supported on *Linux* and *MacOS*. 
1. Instructions for installing *phyloscanner* on either system [are available here](https://github.com/BDI-pathogens/phyloscanner/blob/master/InfoAndInputs/InstallationNotesForMakingTrees.sh). If you have an issue with installation of *phyloscanner*, [please report it here and we will get back to you](https://github.com/BDI-pathogens/phyloscanner/issues). 
2. *Phyloscanner.R.utilities* requires *R* version >= 3.1, and the following *R* packages: ape, argparse, phytools, phangorn, reshape2, data.table, RColorBrewer, colorspace, devtools, grid, gridExtra, ggplot2, ggtree, zoo. Please install these first. When done, the software package can be installed via
```r
devtools:::install_github("olli0601/Phyloscanner.R.utilities")
require(Phyloscanner.R.utilities)
``` 
The software was run on MacOS X El Capitan, and CentOS Linux release 7.3.1611, and tested on both systems.  

## Vignettes for *phyloscanner* analyses of large population-based samples
 

