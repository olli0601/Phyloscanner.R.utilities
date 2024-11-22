# Install miniforge

conda env list
conda remove -n phylo --all

if [ -d $HOME/miniforge3 ]; then
    echo -e "\nminiforge3 already installed\n\n";
else
    echo -e "\ninstalling miniforge3\n\n";
    module load miniforge/3
    miniforge-setup
    eval "$(~/miniforge3/bin/conda shell.bash hook)"
    conda config --set auto_activate_base false
fi

# Create new conda environment called "phylo"
if [ -d $HOME/anaconda3/envs/phylo ]; then
    echo "###############################################"
    echo -e "\nphylo conda environment is already present"
    echo -e "\nIf you wish you re-install please remove the conda environment first with:"
    echo -e "\tconda remove -n phylo --all -y"
    echo -e "\n\n###############################################"
    exit 1
else
    echo -e "\nCreating Conda environment for R packages: phylo"
    eval "$(~/miniforge3/bin/conda shell.bash hook)"
    conda env list
    conda create -n phylo python=2.7 -y
    conda activate phylo
    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
fi


conda env list
# Install initial dependencies
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda


echo -e "\nInstalling Dependencies: R packages via conda"
conda install cython matplotlib pysam biopython mafft raxml curl gcc_linux-64 gxx_linux-64 -y
conda install -c conda-forge r>=4.0.0
conda install -c r r-base r-devtools r-ggplot2 r-network r-scales -y
conda install -c r r-sna r-igraph r-intergraph -y
conda install -c r r-argparse r-ape r-extradistr r-ff r-ggally r-gtable r-kimisc -y
conda install -c r r-pegas r-phangorn r-phytools r-prodlim r-reshape2 -y
conda install -c r r-tidyverse r-viridis r-digest r-gtable -y
conda install -c r r-lazyeval r-mass r-mgcv r-matrix r-lattice -y
conda install -c r r-zeallot r-git2r r-bh r-markdown r-reshape -y

# Install additional dependencies not available via anaconda
echo -e "\nInstalling Additional Dependencies not available via anaconda"
R -e 'options(unzip = "internal");library(devtools);devtools::install_github("briatte/ggnet")'
R -e 'options(unzip = "internal");install.packages("BiocManager", repos="http://cran.us.r-project.org");BiocManager::install(c("GenomeInfoDb", "GenomeInfoDbData", "treeio", "Rsamtools", "RBGL", "GenomicRanges"));'

conda install -c bioconda bioconductor-ggtree
conda install -c bioconda iqtree
# Installing phyloscanner
echo -e "\nInstalling phyloscanner"
cd github
git clone https://github.com/BDI-pathogens/phyloscanner.git
cd phyloscanner/phyloscannerR
R CMD INSTALL .

# Testing
echo -e "\nTesting phyloscanner library loads correctly"
R -e 'library(phyloscannerR)'

conda install -c conda-forge r-ragg
conda install -c conda-forge r-pkgdown
conda install -c conda-forge r-devtools

# Installing phyloscanner.R.utilities
echo -e "\nInstalling Phyloscanner.R.utilities."
R -e 'options(unzip = "internal");library(devtools);devtools::install_github("olli0601/Phyloscanner.R.utilities")'

# Testingla
echo -e "\nTesting phyloscanner.R.utilities"
R -e 'library(Phyloscanner.R.utilities)'

echo "========================================="
echo "INSTALLATION COMPLETE"


echo -e "#############################################################"
echo -e "#############################################################\n\n"

# INSTRUCTIONS OF USE:
echo -e "\t\tINSTRUCTIONS FOR USE"
echo -e "\n#############################################################\n\n"
echo -e "Ensure you add the following lines in your jobscript:\n\n"
echo -e "\tmodule load miniforge/3"
echo -e "\teval ~/miniforge3/bin/conda shell.bash hook..."
echo -e "\tconda activate phylo"
echo -e "\t"'export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH'
echo -e "\n\n"
echo -e "\n#############################################################\n\n"
