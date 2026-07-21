### tre_to_stats_xipho_parallel_part2.R
###
### Convert one midpoint TRE file into one gzipped STAT file.
###
### Input TRE format:
###   position<TAB>Newick_tree
###
### Output STAT format:
###   one row per local tree / midpoint position
###
### This script is designed for the parallel ARGweaver pipeline:
###   argTreeFiles_midpoint/<region>/<region>.<start>-<end>.<iter>.tre.gz
###   argStats_midpoint/<region>/<region>.<start>-<end>.<iter>.stat.gz

library("ape")
library("phytools")
library("plyr")

### -----------------------------
### 1. Parse command-line arguments
### -----------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript tre_to_stats_xipho_parallel_CC.R <treeFile> <statsDir>")
}

treeFile <- args[1]
statsDir <- args[2]

if (!file.exists(treeFile)) {
  stop(sprintf("Missing tree file: %s", treeFile))
}

### -----------------------------
### 2. Load helper functions
### -----------------------------

functionFile1 <- "./treeStatFunctions.R"
functionFile2 <- "./scripts/treeStatFunctions.R"

if (file.exists(functionFile1)) {
  source(functionFile1)
} else if (file.exists(functionFile2)) {
  source(functionFile2)
} else {
  stop("Could not find treeStatFunctions.R in ./ or ./scripts/")
}

### -----------------------------
### 3. Load individual/population key
### -----------------------------

keyFile <- "infoTables/individual-species-key-xipho.txt"

if (!file.exists(keyFile)) {
  stop(sprintf("Missing individual species key: %s", keyFile))
}

pop_ind_key <- read.csv(keyFile, sep = "\t", stringsAsFactors = FALSE)

popNames <- unique(pop_ind_key$Species)

tap_pop <- popNames[grepl("tap", popNames, ignore.case = TRUE)]
xin_pop <- popNames[grepl("xin", popNames, ignore.case = TRUE)]
bel_pop <- popNames[grepl("bel", popNames, ignore.case = TRUE)]

if (length(tap_pop) != 1 || length(xin_pop) != 1 || length(bel_pop) != 1) {
  stop("Could not uniquely identify Tapajos, Xingu, and Belem populations.")
}

pop_counts <- table(pop_ind_key$Species)

total_n      <- sum(pop_counts)
total_half_n <- total_n / 2

tapajos_n <- as.numeric(pop_counts[tap_pop])
xingu_n   <- as.numeric(pop_counts[xin_pop])

### -----------------------------
### 4. Define output file
### -----------------------------

baseName <- basename(treeFile)

statsFile <- sub("\\.tre\\.gz$", ".stat", baseName)
statsFile <- sub("\\.tre$", ".stat", statsFile)
statsFile <- file.path(statsDir, statsFile)

dir.create(statsDir, showWarnings = FALSE, recursive = TRUE)

if (file.exists(statsFile)) file.remove(statsFile)
if (file.exists(paste0(statsFile, ".gz"))) file.remove(paste0(statsFile, ".gz"))

### Region/block name, e.g. scaffold100-553
argBlock <- sub("\\..*", "", baseName)

### -----------------------------
### 5. Read trees
### -----------------------------

trees <- ape::read.tree(treeFile)
coordinates <- names(trees)

### -----------------------------
### 6. Define statistics to output
### -----------------------------

treeStats <- c(
  "chrom",
  "pos",

  ### Age of the full tree/root.
  "TMRCA_all",

  ### Age of the youngest clade containing at least half of all sampled haplotypes.
  "TMRCAH_all",

  ### Tapajos-Xingu half-half clade statistics.
  ### Tap_Xin_TMRCAH:
  ###   age of youngest clade containing >= half Tapajos and >= half Xingu.
  ### Tap_Xin_RTH_original:
  ###   Tap_Xin_TMRCAH / TMRCA_all.
  ### Tap_Xin_RTH_prime:
  ###   Tap_Xin_TMRCAH / TMRCAH_all.
  "Tap_Xin_TMRCAH",
  "Tap_Xin_RTH_original",
  "Tap_Xin_RTH_prime",

  ### Xingu-Belem cross-coalescence statistics.
  ### Xin_Bel_CC:
  ###   age of youngest clade containing >= 1 Xingu and >= 1 Belem.
  ### Xin_Bel_CC_original:
  ###   Xin_Bel_CC / TMRCA_all.
  ### Xin_Bel_CC_prime:
  ###   Xin_Bel_CC / TMRCAH_all.
  "Xin_Bel_CC",
  "Xin_Bel_CC_original",
  "Xin_Bel_CC_prime",

  ### Joint coalescence ratio v2:
  ###   Xin_Bel_CC / Tap_Xin_TMRCAH.
  "Tap_Xin_JCR_v2"
)

write.table(
  t(treeStats),
  file = statsFile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

### -----------------------------
### 7. Helper for selecting joint clades
### -----------------------------
###
### For a pairwise comparison, multiple subtrees may satisfy the
### requested population-count rule.
###
### We choose:
###   1. The youngest such clade.
###   2. If tied, the clade with the smallest joint count.
###      This favors the most specific pairwise clade.

choose_youngest_clade <- function(df, joint_col) {
  if (nrow(df) == 0) return(NULL)

  lowest_age <- min(df$age, na.rm = TRUE)
  df <- df[df$age == lowest_age, ]

  lowest_joint_count <- min(df[[joint_col]], na.rm = TRUE)
  df <- df[df[[joint_col]] == lowest_joint_count, ]

  df[1, ]
}

safe_ratio <- function(numerator, denominator, digits = 4) {
  if (is.na(numerator) || is.na(denominator) || denominator == 0) return(NA)
  round(numerator / denominator, digits = digits)
}

### -----------------------------
### 8. Main loop over local trees
### -----------------------------

for (coord in coordinates) {

  ### Convert tip labels from individual IDs to population labels.
  tree <- trees[[coord]]
  tree <- switchToPopLabels(tree, pop_ind_key)

  tree_stats_df <- data.frame(matrix(nrow = 1, ncol = length(treeStats)))
  colnames(tree_stats_df) <- treeStats

  tree_stats_df$chrom <- argBlock
  tree_stats_df$pos   <- coord

  ### Generate all subtrees/clades and summarize their population composition.
  subtrees <- ape::subtrees(tree)
  subtree_stats_df <- getSubtreeStats(NULL, popNames)

  for (subtree in subtrees) {
    subtree_stats_df <- rbind(
      subtree_stats_df,
      getSubtreeStats(subtree, popNames)
    )
  }

  subtree_stats_df[is.na(subtree_stats_df)] <- 0

  ### Overall tree-depth statistics.
  tree_stats_df$TMRCA_all <- getCladeNage(
    subtree_stats_df,
    "total",
    n = total_n
  )

  tree_stats_df$TMRCAH_all <- getCladeNage(
    subtree_stats_df,
    "total",
    n = total_half_n
  )

  ### Population count columns.
  tap_count_col <- paste(tap_pop, "count", sep = "_")
  xin_count_col <- paste(xin_pop, "count", sep = "_")
  bel_count_col <- paste(bel_pop, "count", sep = "_")

  ### Pairwise joint sizes are used only for tie-breaking.
  subtree_stats_df$tap_xin_joint_count <-
    subtree_stats_df[[tap_count_col]] + subtree_stats_df[[xin_count_col]]

  subtree_stats_df$xin_bel_joint_count <-
    subtree_stats_df[[xin_count_col]] + subtree_stats_df[[bel_count_col]]

  ### -----------------------------
  ### 8a. Tapajos-Xingu half-half clade
  ### -----------------------------

  tap_xin_df <- subtree_stats_df[
    subtree_stats_df[[tap_count_col]] >= tapajos_n / 2 &
      subtree_stats_df[[xin_count_col]] >= xingu_n / 2,
  ]

  tap_xin_clade <- choose_youngest_clade(
    tap_xin_df,
    "tap_xin_joint_count"
  )

  if (!is.null(tap_xin_clade)) {
    tree_stats_df$Tap_Xin_TMRCAH <- tap_xin_clade$age

    tree_stats_df$Tap_Xin_RTH_original <- safe_ratio(
      tap_xin_clade$age,
      tree_stats_df$TMRCA_all
    )

    tree_stats_df$Tap_Xin_RTH_prime <- safe_ratio(
      tap_xin_clade$age,
      tree_stats_df$TMRCAH_all
    )
  }

  ### -----------------------------
  ### 8b. Xingu-Belem cross-coalescence clade
  ### -----------------------------

  xin_bel_cc_df <- subtree_stats_df[
    subtree_stats_df[[xin_count_col]] >= 1 &
      subtree_stats_df[[bel_count_col]] >= 1,
  ]

  xin_bel_cc_clade <- choose_youngest_clade(
    xin_bel_cc_df,
    "xin_bel_joint_count"
  )

  if (!is.null(xin_bel_cc_clade)) {
    tree_stats_df$Xin_Bel_CC <- xin_bel_cc_clade$age

    tree_stats_df$Xin_Bel_CC_original <- safe_ratio(
      xin_bel_cc_clade$age,
      tree_stats_df$TMRCA_all
    )

    tree_stats_df$Xin_Bel_CC_prime <- safe_ratio(
      xin_bel_cc_clade$age,
      tree_stats_df$TMRCAH_all
    )
  }

  ### Ratio comparing Xingu-Belem cross-coalescence age
  ### to Tapajos-Xingu half-half clade age.
  tree_stats_df$Tap_Xin_JCR_v2 <- safe_ratio(
    tree_stats_df$Xin_Bel_CC,
    tree_stats_df$Tap_Xin_TMRCAH
  )

  ### Append one row for this local tree/midpoint position.
  write.table(
    tree_stats_df[1, treeStats],
    file = statsFile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
  )
}

### -----------------------------
### 9. Compress final stat file
### -----------------------------

system(sprintf("gzip -f %s", shQuote(statsFile)))
