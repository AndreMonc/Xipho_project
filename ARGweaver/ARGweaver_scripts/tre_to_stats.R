### R script to run on Linux  ###
### Thus is a script first developed by Hejase et al (2020) for their study on Sporophila
### Andre Moncrieff updates to script primarily involve the ability to handle different number of individuals per species/population
### as well as compute a few new ARG-based statistics
### Original tre_to_stats.R script available here:
### https://github.com/CshlSiepelLab/bird_capuchino_analysis/tree/master/ARG_analysis/smc_to_stats

### Update Contig to Pseudochr ...actually, needs to be a little more nuanced, since I have scaffold and Chrom starts to files
############################################################################
##### read tree file and compute stats for every tree                   ####
####  Stats computed for every tree:                                    ####
####  - clade60age - age of youngest clade of size >= 60                ####
####                 (proxy for general depth of clades in tree)        ####
####  - clade12age - age of youngest clade of size >= 12                ####
####                 (proxy for partial sweeps)                         ####
####  - <pop>_clade12age - pop-specific clade12age                      ####
####                       (proxy for pop-secific sweep                 ####
####  - <pop>_enrich - highest enrichment score for pop of any clade    ####
####  -                ( - log10(pval) under hypergeometric null)       ####      
#### Notes:                                                             ####
####  * all clade12age values are normalized by clade60age              ####
####  * hypergeometric p-values are for p(count >= observed count)      ####
############################################################################
library("ape", "phytools")
library("castor")
############################################################################
# get dir of this script (and also the treeStatFunctions.R script)
# - different code for interactive or batch (Rscript)
script_dir <- "./scripts/"
if(interactive()) {
  script_dir <- dirname(parent.frame(2)$ofile)
} else {
  library(tibble)
  library(tidyr)
  library(plyr)
  library(dplyr)
  get_script_dir <- function(){
    commandArgs() %>% 
       tibble::enframe(name=NULL) %>%
       tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
       dplyr::filter(key == "--file") %>%
       dplyr::pull(value) %>%
       dirname()
  }
  script_dir <- get_script_dir()
}
# print(script_dir)
############################################################################

############################################################################
# script arguments and i/o
# load functions
source(paste(script_dir,"treeStatFunctions.R",sep="/"))
args         <- commandArgs(trailingOnly = TRUE) # TRUE
# example: args         <- c("argTreeFiles/Contig1.50001-1650000.1000.tre.gz","argStats")
treeFile     <- args[1]
statsDir     <- args[2]
############################################################################

############################################################################
# set these based on your specific data set
total_n    <- 66      # total number of haploid samples in data set
#species_n  <- 24      # number of haploid samples per species
tapajos_n  <- 20      # number of haploid samples in tapajos population (RTH_half = 10)
xingu_n    <- 26      # number of haploid samples in xingu population (RTH_half = 13)
belem_n    <- 20      # number of haploid samples in belem population (RTH_half = 10)
# average of half the haploid number is 11 (33/3=11)
  
############################################################################


############################################################################
# output stats file
if (grepl("/Chromosome", treeFile)) {
  statsFile <- sub(".*/Chromosome", "Chromosome", treeFile)
} else if (grepl("/scaffold", treeFile)) {
  statsFile <- sub(".*/scaffold", "scaffold", treeFile)
} else {
  stop("Invalid file name: does not start with Chromosome or scaffold")
}
statsFile <- sub("\\.tre(\\.gz)?$", ".stats", statsFile)
statsFile <- paste(statsDir,statsFile,sep="/")
dir.create(statsDir,showWarnings=FALSE)
res <- file.create(statsFile)
############################################################################

############################################################################
# other setup
# pop-individual key: note that pop key is already modified and pruned
pop_ind_key  <- read.csv("infoTables/individual-species-key-modified_xipho.txt",sep="\t")
popNames     <- unique(as.vector(pop_ind_key[,"Species"]))
if (grepl("/Chromosome", treeFile)) {
  scaffold <- sub(".*/Chromosome", "Chromosome", treeFile)
  scaffold <- sub("\\..*", "", scaffold)
} else if (grepl("/scaffold", treeFile)) {
  scaffold <- sub(".*/scaffold", "scaffold", treeFile)
  scaffold <- sub("\\..*", "", scaffold)
} else {
  stop("Invalid file name: does not start with Chromosome or scaffold")
}
trees        <- ape::read.tree(treeFile)
coordinates  <- names(trees)
# write header
# OLD SPOROPHILA NAMES: treeStats    <- c("chrom", "pos", "clade60age", "clade12age_norm", paste(popNames,"clade12age_norm",sep="_"), paste(popNames, "enrich", sep="_"))
# SPOROPHILA NAMES: treeStats    <- c("chrom", "pos", "TMRCAH_all", "RT12", paste(popNames,"RTH",sep="_"), paste(popNames, "enrich", sep="_"))
treeStats    <- c("chrom", "pos", "TMRCA_all", "TMRCAH_all", "RT11", "TMRCAH_Tap_Xin", "TMRCAH_Tap_Xin_corr1_divbyTMRCAall", "TMRCAH_Tap_Xin_corr2_divbyTMRCAHall", "TMRCAH_Tap_Xin_joint_count", 
                  "TMRCAH_Tap_Xin_tap_count", "TMRCAH_Tap_Xin_xin_count", "test", "TMRCAH_Xin_Bel", "TMRCAH_Xin_Bel_corr1_divbyTMRCAall", "TMRCAH_Xin_Bel_corr2_divbyTMRCAHall", "TMRCAH_Xin_Bel_joint_count",
                  "TMRCAH_Xin_Bel_xin_count", "TMRCAH_Xin_Bel_bel_count", paste(popNames,"RTH",sep="_"), paste(popNames,"RT3Q",sep="_"), paste(popNames,"RTall",sep="_"), paste(popNames, "enrich", sep="_"))
write.table(t(treeStats), file=statsFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)
############################################################################


############################################################################
# main loop - read all trees
counter=1
for(coord in coordinates) {
  # replace individual labels with pop labels and prune individuals not in key
  tree <- trees[[coord]]
  tree <- switchToPopLabels(tree , pop_ind_key)
  totalTreeCounts <- getSubtreeStats(tree, popNames)

  tree_stats_df <- data.frame(matrix(nrow=1,ncol=length(treeStats)))
  colnames(tree_stats_df) <- treeStats
  tree_stats_df$chrom  <- scaffold
  tree_stats_df$pos    <- coord

  # secondary loop - get subtree stats
  subtrees <- ape::subtrees(tree)
  subtree_stats_df <- getSubtreeStats(NULL,popNames)
  for(subtree in subtrees) {
    subtree_stats_df <- rbind( subtree_stats_df , getSubtreeStats(subtree,popNames) )
  }
  # replace NAs with 0s
  subtree_stats_df[is.na(subtree_stats_df)] <- 0
  # add enrichment scores
  subtree_stats_df <- computePopEnrichment(subtree_stats_df, popNames, totalTreeCounts) 
  
  # TMRCA_all, TMRCAH_all, and RT11
  tree_stats_df$TMRCA_all <- getCladeNage(subtree_stats_df,"total",n=total_n)
  tree_stats_df$TMRCAH_all <- getCladeNage(subtree_stats_df,"total",n=total_n/2)
  tree_stats_df$RT11 <- round(getCladeNage(subtree_stats_df,"total",n=11) / tree_stats_df$TMRCAH_all , digits=4)
  
  # Add joint counts to subtree stats dataframe
  subtree_stats_df$tap_xin_joint_count <- subtree_stats_df$xipho_tapajos_count + subtree_stats_df$xipho_xingu_count
  subtree_stats_df$xin_bel_joint_count <- subtree_stats_df$xipho_xingu_count + subtree_stats_df$xipho_belem_count
  
  # Calculate TMRCAH_Tap-Xin (half of Tap + Xin = 23+ haploid samples)--This is the youngest subtree with 23+ Tap + Xin samples
  # Also required that this tree includes 5+ Tap individuls AND 7+ Xin individuals (i.e., 1/4 of samples in these pops)
  # If multiple subtrees have same age and 23+ Tap + Xin samples, then the subtree with fewest individuals (closest to 23) chosen
  # Outputs node age
  # Outputs joint count of Tap + Xin samples in the subtree
  # Outputs the separate counts of Tap and Xin samples in the subtree 
  subset1_df <- subtree_stats_df[subtree_stats_df$tap_xin_joint_count >=23, ]
  subset2_df <- subset1_df[subset1_df$xipho_tapajos_count >=5, ]
  subset3_df <- subset2_df[subset2_df$xipho_xingu_count >=7, ]
  lowest_age <- min(subset3_df$age)
  subset4_df <- subset3_df[subset3_df$age == lowest_age, ]
  lowest_joint_count <- min(subset4_df$tap_xin_joint_count)
  subset5_df <- subset4_df[subset4_df$tap_xin_joint_count == lowest_joint_count, ]
  tree_stats_df$TMRCAH_Tap_Xin <- subset5_df[1, "age"]
  tree_stats_df$TMRCAH_Tap_Xin_corr1_divbyTMRCAall <- round(subset5_df[1, "age"] / tree_stats_df$TMRCA_all , digits=4)
  tree_stats_df$TMRCAH_Tap_Xin_corr2_divbyTMRCAHall <- round(subset5_df[1, "age"] / tree_stats_df$TMRCAH_all , digits=4)
  tree_stats_df$TMRCAH_Tap_Xin_joint_count <- subset5_df[1, "tap_xin_joint_count"]
  tree_stats_df$TMRCAH_Tap_Xin_tap_count <- subset5_df[1, "xipho_tapajos_count"]
  tree_stats_df$TMRCAH_Tap_Xin_xin_count <- subset5_df[1, "xipho_xingu_count"]
  
  # Calculate TMRCAH_Xin-Bel (half of xin + bel = 23+ haploid samples)--This is the youngest subtree with 23+ xin + bel samples
  # Also required that this tree includes 7+ Xin individuls AND 5+ Bel individuals (i.e., 1/4 of samples in these pops)
  # If multiple subtrees have same age and 23+ xin + bel samples, then the subtree with fewest individuals (closest to 23) chosen
  # Outputs node age
  # Outputs joint count of xin + bel samples in the subtree
  # Outputs the separate counts of xin and bel samples in the subtree 
  subset1_df <- subtree_stats_df[subtree_stats_df$xin_bel_joint_count >=23, ]
  subset2_df <- subset1_df[subset1_df$xipho_xingu_count >=7, ]
  subset3_df <- subset2_df[subset2_df$xipho_belem_count >=5, ]
  lowest_age <- min(subset3_df$age)
  subset4_df <- subset3_df[subset3_df$age == lowest_age, ]
  lowest_joint_count <- min(subset4_df$xin_bel_joint_count)
  subset5_df <- subset4_df[subset4_df$xin_bel_joint_count == lowest_joint_count, ]
  tree_stats_df$TMRCAH_Xin_Bel <- subset5_df[1, "age"]
  tree_stats_df$TMRCAH_Xin_Bel_corr1_divbyTMRCAall <- round(subset5_df[1, "age"] / tree_stats_df$TMRCA_all , digits=4)
  tree_stats_df$TMRCAH_Xin_Bel_corr2_divbyTMRCAHall <- round(subset5_df[1, "age"] / tree_stats_df$TMRCAH_all , digits=4)
  tree_stats_df$TMRCAH_Xin_Bel_joint_count <- subset5_df[1, "xin_bel_joint_count"]
  tree_stats_df$TMRCAH_Xin_Bel_xin_count <- subset5_df[1, "xipho_xingu_count"]
  tree_stats_df$TMRCAH_Xin_Bel_bel_count <- subset5_df[1, "xipho_belem_count"]
  
  # RTH', RT3Q, RTall, and enrichment scores
  for(pop in popNames) {
    if (pop=='xipho_tapajos') {
      pop_RTH <- paste(pop,"RTH",sep="_")
      pop_RT3Q <- paste(pop,"RT3Q",sep="_")
      pop_RTall <- paste(pop,"RTall",sep="_")
      pop_enrich     <- paste(pop,"enrich",sep="_")
      tree_stats_df[1,pop_RTH] <- round(getCladeNage(subtree_stats_df,pop,n=tapajos_n/2) / tree_stats_df$TMRCAH_all , digits=4)
      tree_stats_df[1,pop_RT3Q] <- round(getCladeNage(subtree_stats_df,pop,n=tapajos_n*0.75) / tree_stats_df$TMRCAH_all , digits=4)
      tree_stats_df[1,pop_RTall] <- round(getCladeNage(subtree_stats_df,pop,n=tapajos_n) / tree_stats_df$TMRCAH_all , digits=4)
      tree_stats_df[1,pop_enrich]     <- max(subtree_stats_df[,pop_enrich])
    }
    if (pop=='xipho_xingu') {
      pop_RTH <- paste(pop,"RTH",sep="_")
      pop_RT3Q <- paste(pop,"RT3Q",sep="_")
      pop_RTall <- paste(pop,"RTall",sep="_")
      pop_enrich     <- paste(pop,"enrich",sep="_")
      tree_stats_df[1,pop_RTH] <- round(getCladeNage(subtree_stats_df,pop,n=xingu_n/2) / tree_stats_df$TMRCAH_all , digits=4)
      tree_stats_df[1,pop_RT3Q] <- round(getCladeNage(subtree_stats_df,pop,n=20) / tree_stats_df$TMRCAH_all , digits=4)
      tree_stats_df[1,pop_RTall] <- round(getCladeNage(subtree_stats_df,pop,n=xingu_n) / tree_stats_df$TMRCAH_all , digits=4)
      tree_stats_df[1,pop_enrich]     <- max(subtree_stats_df[,pop_enrich])
    }
    if (pop=='xipho_belem') {
      pop_RTH <- paste(pop,"RTH",sep="_")
      pop_RT3Q <- paste(pop,"RT3Q",sep="_")
      pop_RTall <- paste(pop,"RTall",sep="_")
      pop_enrich     <- paste(pop,"enrich",sep="_")
      tree_stats_df[1,pop_RTH] <- round(getCladeNage(subtree_stats_df,pop,n=belem_n/2) / tree_stats_df$TMRCAH_all , digits=4)
      tree_stats_df[1,pop_RT3Q] <- round(getCladeNage(subtree_stats_df,pop,n=belem_n*0.75) / tree_stats_df$TMRCAH_all , digits=4)
      tree_stats_df[1,pop_RTall] <- round(getCladeNage(subtree_stats_df,pop,n=belem_n) / tree_stats_df$TMRCAH_all , digits=4)
      tree_stats_df[1,pop_enrich]     <- max(subtree_stats_df[,pop_enrich])
    }
  }
  

  # write stats line
  write.table(tree_stats_df[1,treeStats], file=statsFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  counter = counter+1
} # end of for(coord)
############################################################################

# gzip file
system(paste("gzip",statsFile,"-f"))

############################################################################
