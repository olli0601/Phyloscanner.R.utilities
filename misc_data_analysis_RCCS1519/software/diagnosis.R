# The phyloscanner run - compare phyloscans in two directories

# Preamble
# This script aims to compare phyloscans in two directories for any sets of individuals. 

# Load the required packages
library(phyloscannerR)
library(tibble)
library(ggpubr)
library(ggplot2)
library(data.table)

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
    "--individuals",
    type = "character",
    default = NA_character_,
    help = " Individuals of which phyloscans will be diagnosed. Enter the IDs of individuals separated by commas, e.g. ID1,ID2,ID3 [default]",
    dest = 'ind'
  ),
  optparse::make_option(
    "--pkg_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to package directory, used as long we don t build an R package [default]",
    dest = 'prj.dir'
  ),
  optparse::make_option(
    "--out_dir_base1",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where all output is stored [default]",
    dest = 'out.dir1'
  ),
  optparse::make_option(
    "--out_dir_base2",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where all output is stored for comparison [default]",
    dest = 'out.dir2'
  ),
  optparse::make_option(
    "--name_out_dir_base1",
    type = "character",
    default = NA_character_,
    help = "Name of the first directory for plots [default]",
    dest = 'name.out.dir1'
  ),
  optparse::make_option(
    "--name_out_dir_base2",
    type = "character",
    default = NA_character_,
    help = "Name of the second directory for plots [default]",
    dest = 'name.out.dir2'
  )
  
)

args <-
  optparse::parse_args(optparse::OptionParser(option_list = option_list))

#
# test
#
if (0) {
  args <- list(
    verbose = T,
    seed = 42,
    if_save_data = T,
    ind = 'AID0039,AID8583,AID6594',
    out.dir1 = NA,
    out.dir2 = NA,
    name.out.dir1 = '80% threshold',
    name.out.dir2 = '50% threshold'
  )
}


#
# use manually specified directories when args$out.dir is NA
#
tmp <- Sys.info()
if (tmp["user"] == "xx4515")
{
  if (is.na(args$out.dir1))
  {
    args$out.dir1 <-
      "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd"
  }
  if (is.na(args$out.dir2))
  {
    args$out.dir2 <-
      "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd"
  }
  if (is.na(args$name.out.dir1)){
    args$name.out.dir1 = '80% threshold'
    args$name.out.dir2 = '50% threshold'
  }
}

# if out.dir1 and out.dir2 are not manually set, default to here()
if (is.na(args$out.dir1))
{
  args$out.dir1 <- here::here()
}
if (is.na(args$out.dir2))
{
  args$out.dir2 <- here::here()
}

#
# Add constants that should not be changed by the user
#
control <- list(
  yintercept_close = 0.025,
  yintercept_dist = 1.0,
  breaks_x = seq(0, 1e4, 500),
  minor_breaks_x = seq(0, 1e4, 100),
  breaks_y = c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25),
  limits_y = c(1e-3, 0.4),
  fill.topology = c(
    "12" = "deepskyblue1",
    "21" = "deepskyblue4",
    "complex.or.no.ancestralestry" = "#FDB863",
    "not.close.or.nonadjacent" = "grey80"
  )
)
out.dir <- file.path(args$out.dir1,'plots')
dir.create(out.dir)
#
aids <- strsplit(args$ind, split = ',')[[1]]

cat(' ---------- Load phyloscanner outputs from directory ', 
    args$out.dir1,' ---------- \n ')
infile <- list.files(args$out.dir1, pattern = '*_phscnetworks.rda', full.names = T)
stopifnot(length(infile)==1)
load(infile)

# format
dw1 <- copy(dw)
setnames(
  dw1,
  c(
    'H1',
    'H2',
    'TREE_ID',
    'CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT',
    'PATRISTIC_DISTANCE'
  ),
  c(
    'host.1',
    'host.2',
    'tree.id',
    'basic.classification',
    'patristic.distance'
  )
)
# 
tmp <- dw1[, list(PTY_RUN = min(PTY_RUN)), by = c('host.1', 'host.2')]
dw1 <- merge(dw1, tmp, by = c('host.1', 'host.2', 'PTY_RUN'))

dc1 <- copy(dc)

cat(' ---------- Load phyloscanner outputs from directory ', 
    args$out.dir2,'----------')
infile <- list.files(args$out.dir2, pattern = '*_phscnetworks.rda', full.names = T)
stopifnot(length(infile)==1)
load(infile)

# format
dw2 <- copy(dw)
setnames(
  dw2,
  c(
    'H1',
    'H2',
    'TREE_ID',
    'CLOSE_AND_ADJACENT_AND_ANCESTRY_CAT',
    'PATRISTIC_DISTANCE'
  ),
  c(
    'host.1',
    'host.2',
    'tree.id',
    'basic.classification',
    'patristic.distance'
  )
)
# 
tmp <- dw2[, list(PTY_RUN = min(PTY_RUN)), by = c('host.1', 'host.2')]
dw2 <- merge(dw2, tmp, by = c('host.1', 'host.2', 'PTY_RUN'))

dc2 <- copy(dc)

aidps <- dw1[host.1 %in% aids & host.2 %in% aids,]
aidps <- unique(subset(aidps, select = c('host.1','host.2')))
tmp <- dw2[host.1 %in% aids & host.2 %in% aids,]
tmp <- unique(subset(tmp, select = c('host.1','host.2')))
aidps <- merge(aidps, tmp, by = c('host.1','host.2'))

dc1 <- dc1[H1%in% aidps$host.1 & H2 %in% aidps$host.2 & CATEGORISATION=='close.and.adjacent.and.ancestry.cat']
dc1 <- dcast(dc1, H1 + H2 ~ TYPE, value.var = 'SCORE')
dc1[`12`+`21`+`complex.or.no.ancestry`>=0.6, linked := T]
dc1[linked==T & `12`/(`12`+`21`)>=0.6, direction:= 'ancestral']
dc1[linked==T & `21`/(`12`+`21`)>=0.6, direction:= 'descendants']
dc1[linked==T & is.na(direction), direction:='complex']
dc1[is.na(linked), direction:='unlinked']

dc2 <- dc2[H1%in% aidps$host.1 & H2 %in% aidps$host.2 & CATEGORISATION=='close.and.adjacent.and.ancestry.cat']
dc2 <- dcast(dc2, H1 + H2 ~ TYPE, value.var = 'SCORE')
dc2[`12`+`21`+`complex.or.no.ancestry`>=0.6, linked := T]
dc2[linked==T & `12`/(`12`+`21`)>=0.6, direction:= 'ancestral']
dc2[linked==T & `21`/(`12`+`21`)>=0.6, direction:= 'descendants']
dc2[linked==T & is.na(direction), direction:='complex']
dc2[is.na(linked), direction:='unlinked']

cat(' ---------- Make plots ----------')
for (i in 1:nrow(aidps)) {
  cat('---------- ',aidps$host.1[i],'-',aidps$host.2,'---------- \n')
  hosts <- c(aidps$host.1[i],aidps$host.2[i])
  tmpc1 <- dc1[(H1==hosts[1] & H2==hosts[2])|(H1==hosts[2] & H2==hosts[1])]$direction
  tmpc2 <- dc2[(H1==hosts[1] & H2==hosts[2])|(H1==hosts[2] & H2==hosts[1])]$direction  
  tmp <- dw1[(host.1==hosts[1] & host.2==hosts[2])|(host.1==hosts[2] & host.2==hosts[1]),'tree.id']
  tmp <- rbind(tmp,dw2[(host.1==hosts[1] & host.2==hosts[2])|(host.1==hosts[2] & host.2==hosts[1]),'tree.id'])
  tmp <- as.integer(gsub('([0-9]+)_to_([0-9]+)','\\1',tmp$tree.id))
  control$limits_x <- range(tmp) + c(-25,25) # bin width
  tmpp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw1), inclusion = "both",control=control)
  g1 <- tmpp$graph +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  tmpp	<- produce.pairwise.graphs2(ptrees=NULL, hosts=hosts,dwin=as_tibble(dw2), inclusion = "both",control=control)
  g2 <- tmpp$graph +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  arrange <- ggarrange(g1 + ggtitle(paste0(args$name.out.dir1,'\n', 'classified as ',tmpc1)),
                       g2 + ggtitle(paste0(args$name.out.dir2,'\n', 'classified as ',tmpc2)), 
                       ncol = 1, nrow = 2, common.legend = T, legend='bottom')
  ggsave(file.path(out.dir,paste0('phyloscan_',hosts[1],'_',hosts[2],'.pdf')),arrange,width = 8, height = 8)
}