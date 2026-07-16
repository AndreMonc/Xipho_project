inputFile <- "./infoTables/ARGblock-coordinates-trimmed.txt"
outFile   <- "./infoTables/region_tree_job_table_midpoint.txt"

lines <- read.table(
  inputFile,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(lines) <- c("region", "scaffold", "startInd", "endInd")

lines$startInd <- as.numeric(lines$startInd)
lines$endInd   <- as.numeric(lines$endInd)
lines$blockLen <- lines$endInd - lines$startInd + 1

lines <- lines[lines$blockLen > 0, ]
lines <- lines[order(lines$region), ]

lines$job_id <- seq_len(nrow(lines))

output <- lines[, c("job_id", "region", "scaffold", "startInd", "endInd")]

write.table(
  output,
  file = outFile,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

cat("Wrote:", outFile, "\n")
cat("Number of region jobs:", nrow(output), "\n")