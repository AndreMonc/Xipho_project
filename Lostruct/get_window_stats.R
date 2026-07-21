#!/usr/bin/env Rscript

library(tidyverse)

input_files <- tibble(
  window_set = c("win500", "win1000"),
  snps_per_window = c(500, 1000),
  file = c(
    "/Users/moncrieff/Library/CloudStorage/Dropbox/Work/Postdoc/Manuscript--Xiphorhynchus/Xipho_revision/loStruct/allele_stats.lostruct.win500.windows.rds",
    "/Users/moncrieff/Library/CloudStorage/Dropbox/Work/Postdoc/Manuscript--Xiphorhynchus/Xipho_revision/loStruct/allele_stats.lostruct.win1000.windows.rds"
  )
)

summarize_windows <- function(window_set, snps_per_window, file) {
  x <- readRDS(file)
  
  if (!all(c("start", "end") %in% names(x))) {
    stop("File is missing start/end columns: ", file)
  }
  
  x %>%
    mutate(window_bp = end - start) %>%
    summarize(
      window_set = window_set,
      snps_per_window = snps_per_window,
      n_windows = n(),
      mean_bp = mean(window_bp, na.rm = TRUE),
      median_bp = median(window_bp, na.rm = TRUE),
      sd_bp = sd(window_bp, na.rm = TRUE),
      min_bp = min(window_bp, na.rm = TRUE),
      max_bp = max(window_bp, na.rm = TRUE)
    )
}

summary_tbl <- pmap_dfr(input_files, summarize_windows)

write_tsv(
  summary_tbl,
  "/Users/moncrieff/Library/CloudStorage/Dropbox/Work/Postdoc/Manuscript--Xiphorhynchus/Xipho_revision/loStruct/lostruct_window_size_summary.tsv"
)

print(summary_tbl)