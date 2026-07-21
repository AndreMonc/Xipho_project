### bed_to_tre_midpoint_xipho_parallel.R
### Convert one ARG BED file into one gzipped midpoint TREE file.

args <- commandArgs(trailingOnly = TRUE)

bedFile  <- args[1]
outDir   <- args[2]
scaffold <- args[3]
startInd <- as.numeric(args[4])
endInd   <- as.numeric(args[5])

mcmcIter <- sub(".*_out\\.([0-9]+)\\.bed\\.gz$", "\\1", basename(bedFile))
argBlock <- sub("_out\\.[0-9]+\\.bed\\.gz$", "", basename(bedFile))

region_out_dir <- file.path(outDir, argBlock)
dir.create(region_out_dir, showWarnings = FALSE, recursive = TRUE)

treeFile <- sprintf(
  "%s/%s.%s-%s.%s.tre",
  region_out_dir,
  argBlock,
  format(startInd, scientific = FALSE),
  format(endInd, scientific = FALSE),
  mcmcIter
)

command <- sprintf(
  "tabix %s %s:%.0f-%.0f",
  shQuote(bedFile),
  scaffold,
  startInd,
  endInd
)

treeTable <- read.table(
  pipe(command),
  header = FALSE,
  stringsAsFactors = FALSE
)

names(treeTable) <- c("chrom", "chromStart", "chromEnd", "MCMC", "tree")

treeTable$trimStart <- pmax(treeTable$chromStart, startInd - 1)
treeTable$trimEnd   <- pmin(treeTable$chromEnd, endInd)

treeTable <- treeTable[treeTable$trimEnd > treeTable$trimStart, ]

treeTable$midpoint <- floor((treeTable$trimStart + treeTable$trimEnd) / 2) + 1

out <- data.frame(
  pos = format(treeTable$midpoint, scientific = FALSE),
  tree = treeTable$tree
)

write.table(
  out,
  file = treeFile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

system(sprintf("gzip -f %s", shQuote(treeFile)))