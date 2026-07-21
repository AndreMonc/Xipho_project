workDir   <- "./infoTables"
inputFile <- file.path(workDir, "ARGblock-coordinates.txt")
outFile   <- file.path(workDir, "ARGblock-coordinates-trimmed.txt")
infoFile  <- file.path(workDir, "ARGblock-coordinates-info.txt")

min_block_len  <- 101000
segment_len    <- 1000
trim_len       <- 50000

lines <- read.table(
  inputFile,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Expected columns:
# region chromosome start end
colnames(lines) <- c("blockName", "scaffold", "startInd", "endInd")

# Extract numeric block index after hyphen:
# scaffold10-135 -> 135
lines$blockIndex <- as.integer(sub(".*-", "", lines$blockName))

total_length_contigs <- function(df) {
  sum(df$endInd - df$startInd + 1)
}

divideBlocksToRegionsBySegments <- function(args_df, trim_len, segment_len) {
  ordered_args_df <- args_df[order(args_df$blockIndex), ]

  for (i in seq_len(nrow(ordered_args_df))) {

    if (i == 1 ||
        ordered_args_df[i, "blockIndex"] != ordered_args_df[i - 1, "blockIndex"] + 1) {
      ordered_args_df[i, "startInd"] <- ordered_args_df[i, "startInd"] + trim_len
    } else {
      ordered_args_df[i, "startInd"] <- ordered_args_df[i - 1, "endInd"] + 1
    }

    num_segments <- floor(
      (ordered_args_df[i, "endInd"] - trim_len - ordered_args_df[i, "startInd"] + 1) /
        segment_len
    )

    ordered_args_df[i, "endInd"] <-
      ordered_args_df[i, "startInd"] + num_segments * segment_len - 1
  }

  ordered_args_df
}

num_scaffolds <- length(unique(lines$scaffold))

lines_long  <- subset(lines, endInd - startInd + 1 >= min_block_len)
lines_short <- subset(lines, endInd - startInd + 1 <  min_block_len)

num_long_arg_blocks    <- nrow(lines_long)
num_short_arg_blocks   <- nrow(lines_short)

total_len_long_blocks  <- total_length_contigs(lines_long)
total_len_short_blocks <- total_length_contigs(lines_short)

write(paste("number of analyzed scaffolds is", num_scaffolds), infoFile)
write(paste("number of ARG blocks of size at least", min_block_len,
            "bp is", num_long_arg_blocks), infoFile, append = TRUE)
write(paste("total length of these ARG blocks is:", total_len_long_blocks),
      infoFile, append = TRUE)
write(paste("number of ARG blocks shorter than", min_block_len,
            "bp is", num_short_arg_blocks), infoFile, append = TRUE)
write(paste("total length of these ARG blocks is:", total_len_short_blocks),
      infoFile, append = TRUE)

output <- data.frame()

for (scaffold in unique(lines_long$scaffold)) {
  arg_blocks_per_scaffold <- lines_long[
    lines_long$scaffold == scaffold,
  ]

  output <- rbind(
    output,
    divideBlocksToRegionsBySegments(arg_blocks_per_scaffold, trim_len, segment_len)
  )
}

write.table(
  output[, c("blockName", "scaffold", "startInd", "endInd")],
  file = outFile,
  row.names = FALSE,
  quote = FALSE,
  col.names = FALSE,
  sep = "\t"
)