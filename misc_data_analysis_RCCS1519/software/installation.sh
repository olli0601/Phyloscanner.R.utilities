echo -e "\nCreating Conda environment for R packages: phylostan"
conda create -n phylostan -y python=2.7
conda activate phylostan
conda install -c conda-forge r-base=4.0 -y
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
echo -e "\nInstalling Dependencies: R packages via conda"
conda install cython matplotlib pysam biopython mafft raxml curl -y
conda install -c conda-forge r-devtools r-ggplot2 r-network r-scales r-sna r-igraph r-intergraph r-knitr r-rcolorbrewer r-testthat r-argparse r-ape r-extradistr r-ff r-ggally r-gtable r-kimisc r-pegas r-phangorn r-phytools r-prodlim r-reshape2 r-tidyverse bioconductor-treeio r-viridis r-digest r-gtable r-lazyeval r-mass r-mgcv r-matrix r-lattice r-nlme r-rlang r-r6 r-rcpp r-cli=3.2.0 r-tibble r-assertthat r-fansi r-pillar r-vctrs r-backports r-ellipsis r-glue r-zeallot r-vctrs r-findpython r-jsonlite r-colorspace r-data.table r-callr r-processx r-git2r r-httr r-mime r-pkgbuild r-pkgload r-rstudioapi r-rcmdcheck r-roxygen2 r-brew r-commonmark r-purrr r-stringi r-stringr r-xml2 r-evaluate r-dplyr r-bh r-tidyselect r-dtplyr r-highr r-markdown r-xfun r-quadprog r-numderiv r-reshape r-rmarkdown r-tinytex r-units r-sf -y
conda install bioconductor-rsamtools bioconductor-rbgl bioconductor-genomeinfodb bioconductor-genomicranges bioconductor-genomeinfodbdata -y
conda install -c bioconda bioconductor-treeio bioconductor-ggtree -y
conda install -c bioconda iqtree -y
conda install -c conda-forge libiconv -y 
conda install -c conda-forge r-usethis r-rstan -y
# Install additional dependencies not available via anaconda

R -e 'if(!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")};BiocManager::install("RBGL")'

echo -e "\nInstalling Additional Dependencies not available via anaconda"
R -e 'options(unzip = "internal");library(devtools);devtools::install_github("briatte/ggnet");install.packages("https://cran.r-project.org/src/contrib/Archive/rvcheck/rvcheck_0.1.8.tar.gz", repos = NULL);install_github("BDI-pathogens/phyloscanner/phyloscannerR", dependencies=TRUE, INSTALL_opts="--no-staged-install")'
conda install r-seqinr r-optparse r-here -y
