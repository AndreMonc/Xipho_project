#!/usr/bin/env Rscript

# ============================================================
# rCNV pipeline for one chromosome/scaffold VCF
#
# Usage:
#   Rscript run_rCNV_one_chrom.R input.vcf.gz output_dir
#
# Outputs per scaffold:
#   *.allele_info.full.tsv       traditional allele.info() output
#   *.deviants.full.tsv          traditional dupGet() output
#   *.cnv.full.tsv               traditional cnv() output
#   *.allele_info_WGS.full.tsv   WGS likelihood-based allele.info.WGS() output
#   *.summary.tsv                scaffold-level run summary
# ============================================================

suppressPackageStartupMessages(library(rCNV))

# Disable automatic plotting
pdf(NULL)

# -----------------------------
# Parse command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript run_rCNV_one_chrom.R input.vcf.gz output_dir")
}

vcf.file.path <- args[1]
outdir <- args[2]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

chrom <- sub("\\.vcf\\.gz$", "", basename(vcf.file.path))

message("rCNV version: ", packageVersion("rCNV"))
message("========================================")
message("Running rCNV on: ", chrom)
message("Input VCF: ", vcf.file.path)
message("Output dir: ", outdir)
message("========================================")

# -----------------------------
# Import VCF
# -----------------------------
xipho <- readVCF(vcf.file.path, verbose = FALSE)

# -----------------------------
# Filter SNPs by missingness
# -----------------------------
mss_snp_pre <- get.miss(xipho, verbose = FALSE, plot = FALSE)

write.table(
  mss_snp_pre$perSNP,
  file.path(outdir, paste0(chrom, ".pre_filter.snp_missingness.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

bad_snps <- which(mss_snp_pre$perSNP$f_miss > 0.2)

if (length(bad_snps) > 0) {
  message("Removing ", length(bad_snps), " SNPs with >20% missing data")
  xipho$vcf <- xipho$vcf[-bad_snps, , drop = FALSE]
} else {
  message("No SNPs removed for >20% missingness")
}

min_snps <- 1000

message(
  "Dimensions after SNP filtering: ",
  nrow(xipho$vcf), " rows x ", ncol(xipho$vcf), " columns"
)

if (nrow(xipho$vcf) < min_snps) {
  warning("Skipping ", chrom, ": fewer than ", min_snps,
          " SNPs remain after SNP filtering")
  
  summary_tab <- data.frame(
    chromosome_file = chrom,
    status = "skipped_too_few_snps_after_snp_filter",
    n_sites_after_snp_filter = nrow(xipho$vcf),
    n_samples_after_snp_filter = ncol(xipho$vcf) - 9
  )
  
  write.table(
    summary_tab,
    file.path(outdir, paste0(chrom, ".summary.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  quit(save = "no", status = 0)
}

# -----------------------------
# Filter samples by missingness
# -----------------------------
mss_sample_post_snp <- get.miss(xipho, verbose = FALSE, plot = FALSE)

write.table(
  mss_sample_post_snp$perSample,
  file.path(outdir, paste0(chrom, ".post_snp_filter.sample_missingness.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

bad_samples <- which(mss_sample_post_snp$perSample$f_miss > 0.2) + 9

if (length(bad_samples) > 0) {
  message("Removing ", length(bad_samples),
          " samples with >20% missing data after SNP filtering")
  
  keep_cols <- setdiff(seq_len(ncol(xipho$vcf)), bad_samples)
  xipho$vcf <- xipho$vcf[, keep_cols, with = FALSE]
} else {
  message("No samples removed for >20% missingness after SNP filtering")
}

message(
  "Dimensions after sample filtering: ",
  nrow(xipho$vcf), " rows x ", ncol(xipho$vcf), " columns"
)

n_samples <- ncol(xipho$vcf) - 9

if (n_samples < 2) {
  warning("Skipping ", chrom, ": fewer than 2 samples remain after filtering")
  
  summary_tab <- data.frame(
    chromosome_file = chrom,
    status = "skipped_too_few_samples_after_filter",
    n_sites_after_filter = nrow(xipho$vcf),
    n_samples_after_filter = n_samples
  )
  
  write.table(
    summary_tab,
    file.path(outdir, paste0(chrom, ".summary.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  quit(save = "no", status = 0)
}

# -----------------------------
# Prepare filtered VCF object
# -----------------------------
xipho <- data.frame(xipho$vcf)

message(
  "Filtered VCF dimensions: ",
  nrow(xipho), " SNPs x ", ncol(xipho), " columns"
)

# -----------------------------
# Calculate heterozygosity and global Fis
# -----------------------------
hz <- h.zygosity(xipho, verbose = FALSE)
fis <- mean(hz$Fis, na.rm = TRUE)

write.table(
  hz,
  file.path(outdir, paste0(chrom, ".heterozygosity.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# Generate allele-depth and genotype tables
# -----------------------------
ad.tab <- hetTgen(xipho, info.type = "AD", verbose = FALSE)
gt <- hetTgen(xipho, info.type = "GT", verbose = FALSE)

# Replace blank allele-depth cells, if present
ad.tab[ad.tab == ""] <- "0,0"

# -----------------------------
# Traditional rCNV path:
# correction -> normalization -> allele.info -> dupGet -> cnv
# -----------------------------
ad.corrected <- ad.correct(
  ad.tab,
  gt.table = gt
)

ad.nor <- cpm.normal(
  ad.corrected,
  method = "MedR",
  verbose = FALSE,
  plot = FALSE
)

A.info <- allele.info(
  X = ad.corrected,
  x.norm = ad.nor,
  Fis = fis,
  plot.allele.cov = FALSE,
  verbose = FALSE
)

deviants <- dupGet(
  A.info,
  Fis = fis,
  test = c("z.05", "chi.05"),
  plot = FALSE,
  verbose = FALSE
)

CV <- cnv(
  A.info,
  test = c("z.05", "chi.05"),
  filter = "intersection",
  WGS = TRUE,
  plot = FALSE,
  verbose = FALSE
)

# -----------------------------
# WGS likelihood-based rCNV path:
# allele.info.WGS output saved separately
# -----------------------------
out <- tryCatch(
  {
    allele.info.WGS(
      ad = ad.tab,
      gt = gt,
      fis = fis,
      parallel = FALSE
    )
  },
  error = function(e) {
    warning("allele.info.WGS() failed for ", chrom, ": ", conditionMessage(e))
    
    data.frame(
      chromosome_file = chrom,
      status = "allele_info_WGS_failed",
      error_message = conditionMessage(e)
    )
  }
)

# -----------------------------
# Add scaffold provenance
# -----------------------------
A.info$chromosome_file <- chrom
deviants$chromosome_file <- chrom
CV$chromosome_file <- chrom
out$chromosome_file <- chrom

# -----------------------------
# Create scaffold-level summary
# -----------------------------
summary_tab <- data.frame(
  chromosome_file = chrom,
  status = "completed",
  n_sites = nrow(CV),
  fis = fis,
  n_Ainfo_rows = nrow(A.info),
  n_deviants_rows = nrow(deviants),
  n_CV_rows = nrow(CV),
  n_WGS_rows = nrow(out),
  n_cnv = sum(CV$dup.stat == "cnv", na.rm = TRUE),
  n_non_cnv = sum(CV$dup.stat == "non-cnv", na.rm = TRUE)
)

# -----------------------------
# Save outputs
# -----------------------------
write.table(
  A.info,
  file.path(outdir, paste0(chrom, ".allele_info.full.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  deviants,
  file.path(outdir, paste0(chrom, ".deviants.full.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  CV,
  file.path(outdir, paste0(chrom, ".cnv.full.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  out,
  file.path(outdir, paste0(chrom, ".allele_info_WGS.full.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  summary_tab,
  file.path(outdir, paste0(chrom, ".summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Finished rCNV for: ", chrom)