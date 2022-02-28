# The phyloscanner run - find transmission networks 

# Preamble
# This script aims to find transmission networks and the most likely transmission chains 
# using the phyloscanner outputs.

# Load the required packages
library(data.table)
library(tidyverse)
library(dplyr)
library(glue)
library(igraph)
library(RBGL)
library(phyloscannerR)

#
# Define input arguments that can be changed by users
#
option_list <- list(
  optparse::make_option(
    c("-v", "--verbose"),
    action = "store_true",
    default = FALSE,
    help = "Print extra output [default]",
    dest = "verbose"
  ),
  optparse::make_option(
    "--seed",
    type = "integer",
    default = 42L,
    help = "Random number seed [default %default]",
    dest = "seed"
  ),
  optparse::make_option(
    "--save_data",
    action = "store_true",
    default = TRUE,
    help = "Save data [default]",
    dest = 'if_save_data'
  ),
  optparse::make_option(
    "--env_name",
    type = "character",
    default = 'phyloenv',
    help = "Conda environment name [default]",
    dest = 'env_name'
  ),
  optparse::make_option(
    "--classification_rule",
    type = "character",
    default = 'o',
    help = "Rules for classifying linked and directed pairs. It takes values o or m. 
    o: a pair is linked if the linkage score > 0.6, and directed if the ancestral / (ancestral + descedant) > 0.6.
    m: a pair is linked if the linkage score > 0.5, and directed if the ancestral / all > 0.33. 
    b: both [default]",
    dest = 'classif_rule'
  ),
  optparse::make_option(
    "--pkg_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to package directory, used as long we don t build an R package [default]",
    dest = 'prj.dir'
  ),
  optparse::make_option(
    "--out_dir_base",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where all the tree outputs are stored [default]",
    dest = 'out.dir'
  ),
  optparse::make_option(
    "--phylo_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where the phyloscanner outputs are stored [default]",
    dest = 'phylo.dir'
  ),
  optparse::make_option(
    "--date",
    type = 'character',
    default = '2022-02-04',
    metavar = '"YYYY-MM-DD"',
    help = 'As of date to extract data from.  Defaults to today.',
    dest = 'date'
  )
)

args <-
  optparse::parse_args(optparse::OptionParser(option_list = option_list))

#
# test
#
if(0)
{
  args <- list(
    verbose = T,
    seed = 42,
    if_save_data = T,
    date = '2022-02-04',
    env_name = 'phylo',
    classif_rule = 'o',
    out.dir = NA,
    prj.dir = NA,
    phylo.dir = NA
  )
}


#
# use manually specified directories when args$out.dir is NA
#
tmp <- Sys.info()
if (tmp["user"] == "xx4515")
{
  if (is.na(args$out.dir))
  {
    args$out.dir <-
      "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/"
  }
  if (is.na(args$prj.dir))
  {
    args$prj.dir <-
      "~/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/"
  }
  if (is.na(args$phylo.dir))
  {
    args$phylo.dir <-
      "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd"
  }
}

# if prj.dir and out.dir are not manually set, default to here()
if (is.na(args$prj.dir))
{
  args$prj.dir <- here::here()
}
if (is.na(args$out.dir))
{
  args$out.dir <- here::here()
}
if (is.na(args$phylo.dir))
{
  args$phylo.dir <- here::here()
}

#
# Add constants that should not be changed by the user
#
max.per.run <- 4900

# Set default output directories relative to out.dir
args$out.dir.data <-
  file.path(args$out.dir, paste0(args$date, "_phsc_input"))
args$out.dir.work <-
  file.path(args$out.dir, paste0(args$date, "_phsc_work"))
args$out.dir.output <-
  file.path(args$out.dir, paste0(args$date, "_phsc_output"))


cat('Load phyloscanner outputs... \n')
infiles	<- data.table(F=list.files(args$phylo.dir, pattern='*workspace.rda$', full.names=TRUE))
infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(F)))]
setkey(infiles, PTY_RUN)
dca	<- infiles[, {
  cat(PTY_RUN,'\n')
  load(F)
  dc
}, by='PTY_RUN']
dwina <- infiles[, {
  cat(PTY_RUN,'\n')
  load(F)
  dwin
}, by='PTY_RUN']


# Change the format
dca[,CNTRL1:=FALSE]
dca[,CNTRL2:=FALSE]
dca[grepl('CNTRL-',host.1),CNTRL1:=TRUE]
dca[grepl('CNTRL-',host.2),CNTRL2:=TRUE]
dca[grepl('CNTRL-',host.1),host.1:=gsub('CNTRL-','',host.1)]
dca[grepl('CNTRL-',host.2),host.2:=gsub('CNTRL-','',host.2)]
dwina[,CNTRL1:=FALSE]
dwina[,CNTRL2:=FALSE]
dwina[grepl('CNTRL-',host.1),CNTRL1:=TRUE]
dwina[grepl('CNTRL-',host.2),CNTRL2:=TRUE]
dwina[grepl('CNTRL-',host.1),host.1:=gsub('CNTRL-','',host.1)]
dwina[grepl('CNTRL-',host.2),host.2:=gsub('CNTRL-','',host.2)]

# Sort
tmp			<- subset(dwina, host.1>host.2)
setnames(tmp, c('host.1','host.2','paths12','paths21','nodes1','nodes2','CNTRL1','CNTRL2'),
         c('host.2','host.1','paths21','paths12','nodes2','nodes1','CNTRL2','CNTRL1'))
set(tmp, NULL, 'close.and.contiguous.and.directed.cat',
    tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',close.and.contiguous.and.directed.cat)))])
set(tmp, NULL, 'close.and.adjacent.and.directed.cat',
    tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',close.and.adjacent.and.directed.cat)))])
set(tmp, NULL, 'close.and.contiguous.and.ancestry.cat',
    tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',close.and.contiguous.and.ancestry.cat)))])
set(tmp, NULL, 'close.and.adjacent.and.ancestry.cat',
    tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',close.and.adjacent.and.ancestry.cat)))])
dwina		<- rbind(subset(dwina, !(host.1>host.2)), tmp)

# Sort
tmp			<- subset(dca, host.1>host.2)
setnames(tmp, c('host.1','host.2','CNTRL1','CNTRL2'),
         c('host.2','host.1','CNTRL2','CNTRL1'))
set(tmp, NULL, 'type',
    tmp[,gsub('xx','21',gsub('21','12',gsub('12','xx',type)))])
dca		<- rbind(subset(dca, !(host.1>host.2)), tmp)  
tmp <- unique(subset(dwina,select=c('PTY_RUN','host.1','host.2')))
tmp <- tmp[,list(PTY_RUN=PTY_RUN[1]),by=c('host.1','host.2')]
dwina <- merge(dwina,tmp, by=c('host.1','host.2'))
dca <- merge(dca,tmp, by=c('host.1','host.2','PTY_RUN'))

if('PTY_RUN.y' %in% colnames(dwina)){
  dwina$PTY_RUN.y=NULL
  setnames(dwina,'PTY_RUN.x','PTY_RUN',skip_absent=T)
}
if('PTY_RUN.y' %in% colnames(dca)){
  dca$PTY_RUN.y=NULL
  setnames(dca,'PTY_RUN.x','PTY_RUN',skip_absent=T)
}




if(args$classif_rule=='o'|args$classif_rule='b'){
  control <- list(linked.group='close.and.adjacent.cat',
                  linked.no='not.close.or.nonadjacent',
                  linked.yes='close.and.adjacent', 
                  dir.group = "close.and.adjacent.and.directed.cat",
                  conf.cut=0.6, 
                  neff.cut=3,
                  weight.complex.or.no.ancestry=0.5)
  # Find pairs
  tmp <- find.pairs.in.networks(dwina, dca, control=control, verbose=TRUE)
  dpl <- copy(tmp$network.pairs)
  dc <- copy(tmp$relationship.counts)
  dw <- copy(tmp$windows)
  save(dpl, dc, dw, file=file.path(args$phylo.dir,'Rakai_phscnetworks_allpairs_ruleo.rda'))
  # Find chains
  tmp <- find.networks(dc, control=control, verbose=TRUE)
  dnet <- copy(tmp$transmission.networks)
  dchain <- copy(tmp$most.likely.transmission.chains)
  save(dpl, dc, dw, dnet, dchain, file=file.path(args$phylo.dir,'Rakai_phscnetworks_ruleo.rda'))
}

if(args$classif_rule=='m'|args$classif_rule='b'){
  # Find pairs
  control <- list(linked.group='close.and.adjacent.cat',
                  linked.no='not.close.or.nonadjacent',
                  linked.yes='close.and.adjacent',
                  dir.group="close.and.adjacent.and.ancestry.cat",
                  conf.cut=0.5, 
                  neff.cut=3,
                  weight.complex.or.no.ancestry=0.5)
  tmp <- find.pairs.in.networks(dwina, dca, control=control)
  dpl <- copy(tmp$network.pairs)
  dc <- copy(tmp$relationship.counts)
  dw <- copy(tmp$windows)
  save(dpl, dc, dw, file=file.path(args$phylo.dir,'Rakai_phscnetworks_allpairs_rulem.rda'))
  # Find chains
  control<- list(linked.group='close.and.adjacent.cat',
                 linked.no='not.close.or.nonadjacent',
                 linked.yes='close.and.adjacent',
                 dir.group="close.and.adjacent.and.ancestry.cat", 
                 neff.cut=3, 
                 weight.complex.or.no.ancestry=0.5)
  tmp <- find.networks(dc, control=control, verbose=TRUE)
  dnet <- copy(tmp$transmission.networks)
  dchain <- copy(tmp$most.likely.transmission.chains)
  save(dpl, dc, dw, dnet, dchain, file=file.path(args$phylo.dir,'Rakai_phscnetworks_rulem.rda'))
}

  
if(!args$classif_rule %in% c('o','m','b')){
  stop('Please input --classification_rule as o, m or b')
}



