#!/usr/bin/env Rscript

library(data.table)

# ----------------
# Input/output files
# ----------------
base_dir <- "/Users/moncrieff/Library/CloudStorage/Dropbox/Work/Postdoc/Manuscript--Xiphorhynchus/Xipho_revision/Fst_Dxy_Pi"

fst_file <- file.path(base_dir, "output/pixy_output_fst.txt")
dxy_file <- file.path(base_dir, "output/pixy_output_dxy.txt")
pi_file  <- file.path(base_dir, "output/pixy_output_pi.txt")

windows_file <- file.path(base_dir, "windows.bed")
out_file <- file.path(base_dir, "pixy_merged_fst_dxy_pi.tsv")

dup_file <- file.path(base_dir, "duplicate_windows_after_merge.tsv")

# ----------------
# Read master windows
# ----------------
windows <- fread(
  windows_file,
  col.names = c("chromosome", "start", "end")
)

windows[, window_order := .I]

# ----------------
# Helper for pair labels
# ----------------
pair_label <- function(pop1, pop2) {
  x <- paste(pop1, pop2, sep = "_")
  
  fifelse(x %in% c("Tap_Xin", "Xin_Tap"), "Tap_Xin",
          fifelse(x %in% c("Tap_Bel", "Bel_Tap"), "Tap_Bel",
                  fifelse(x %in% c("Xin_Bel", "Bel_Xin"), "Xin_Bel", NA_character_)))
}

# ----------------
# Fst
# ----------------
fst <- fread(fst_file)

fst[, pair := pair_label(pop1, pop2)]
fst <- fst[!is.na(pair)]

fst_wide <- dcast(
  fst,
  chromosome + window_pos_1 + window_pos_2 ~ pair,
  value.var = "avg_wc_fst"
)

setnames(
  fst_wide,
  old = c("window_pos_1", "window_pos_2", "Tap_Xin", "Tap_Bel", "Xin_Bel"),
  new = c("start", "end", "Tap_Xin_Fst", "Tap_Bel_Fst", "Xin_Bel_Fst"),
  skip_absent = TRUE
)

fst_snps <- fst[, .(
  no_snps_fst = min(no_snps, na.rm = TRUE)
), by = .(chromosome, start = window_pos_1, end = window_pos_2)]

fst_wide <- merge(
  fst_wide,
  fst_snps,
  by = c("chromosome", "start", "end"),
  all.x = TRUE
)

# ----------------
# Dxy
# ----------------
dxy <- fread(dxy_file)

dxy[, pair := pair_label(pop1, pop2)]
dxy <- dxy[!is.na(pair)]

dxy_wide <- dcast(
  dxy,
  chromosome + window_pos_1 + window_pos_2 ~ pair,
  value.var = "avg_dxy"
)

setnames(
  dxy_wide,
  old = c("window_pos_1", "window_pos_2", "Tap_Xin", "Tap_Bel", "Xin_Bel"),
  new = c("start", "end", "Tap_Xin_Dxy", "Tap_Bel_Dxy", "Xin_Bel_Dxy"),
  skip_absent = TRUE
)

dxy_sites <- dxy[, .(
  no_sites_dxy = min(no_sites, na.rm = TRUE)
), by = .(chromosome, start = window_pos_1, end = window_pos_2)]

dxy_wide <- merge(
  dxy_wide,
  dxy_sites,
  by = c("chromosome", "start", "end"),
  all.x = TRUE
)

# ----------------
# Pi
# ----------------
pi <- fread(pi_file)

pi_wide <- dcast(
  pi,
  chromosome + window_pos_1 + window_pos_2 ~ pop,
  value.var = "avg_pi"
)

setnames(
  pi_wide,
  old = c("window_pos_1", "window_pos_2", "Tap", "Xin", "Bel"),
  new = c("start", "end", "Tap_pi", "Xin_pi", "Bel_pi"),
  skip_absent = TRUE
)

pi_sites <- pi[, .(
  no_sites_pi = min(no_sites, na.rm = TRUE)
), by = .(chromosome, start = window_pos_1, end = window_pos_2)]

pi_wide <- merge(
  pi_wide,
  pi_sites,
  by = c("chromosome", "start", "end"),
  all.x = TRUE
)

# ----------------
# Merge everything to master windows
# ----------------
merged <- merge(
  windows,
  fst_wide,
  by = c("chromosome", "start", "end"),
  all.x = TRUE
)

merged <- merge(
  merged,
  dxy_wide,
  by = c("chromosome", "start", "end"),
  all.x = TRUE
)

merged <- merge(
  merged,
  pi_wide,
  by = c("chromosome", "start", "end"),
  all.x = TRUE
)

# ----------------
# Check for accidental duplicate windows
# ----------------
dup_check <- merged[, .N, by = .(chromosome, start, end)][N > 1]

if (nrow(dup_check) > 0) {
  fwrite(dup_check, dup_file, sep = "\t", quote = FALSE)
  stop(paste("Duplicate windows created during merge. See:", dup_file))
}

# Restore original windows.bed order
setorder(merged, window_order)

# ----------------
# Final column order
# ----------------
final_cols <- c(
  "chromosome", "start", "end",
  "Tap_Xin_Fst", "Tap_Bel_Fst", "Xin_Bel_Fst",
  "Tap_Xin_Dxy", "Tap_Bel_Dxy", "Xin_Bel_Dxy",
  "Tap_pi", "Xin_pi", "Bel_pi",
  "no_snps_fst", "no_sites_dxy", "no_sites_pi"
)

for (col in final_cols) {
  if (!col %in% names(merged)) {
    merged[, (col) := NA]
  }
}

merged <- merged[, ..final_cols]

# ----------------
# Write output
# ----------------
fwrite(merged, out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Wrote:", out_file, "\n")
cat("Rows written:", nrow(merged), "\n")
cat("Expected rows from windows.bed:", nrow(windows), "\n")