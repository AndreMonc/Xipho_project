args <- commandArgs(trailingOnly = TRUE)

arg_block <- args[1]
iter      <- args[2]
smc_dir   <- args[3]
log_dir   <- args[4]
out_dir   <- args[5]

smc2bed_path <- "/home/amonc/ARGweaver/bin/smc2bed"

log_file <- sprintf("%s/%s_out.log", log_dir, arg_block)

smc_file <- sprintf(
  "%s/%s/%s_out.%s.smc.gz",
  smc_dir,
  arg_block,
  arg_block,
  iter
)

region_out_dir <- sprintf("%s/%s", out_dir, arg_block)
dir.create(region_out_dir, showWarnings = FALSE, recursive = TRUE)

out_bed_file <- sprintf(
  "%s/%s_out.%s.bed.gz",
  region_out_dir,
  arg_block,
  iter
)

if (!file.exists(smc_file)) {
  stop(sprintf("Missing SMC file: %s", smc_file))
}

if (!file.exists(log_file)) {
  stop(sprintf("Missing log file: %s", log_file))
}

cat("Creating bed file:", out_bed_file, "\n")

cmd <- sprintf(
  "%s --log-file %s %s | bgzip > %s",
  smc2bed_path,
  shQuote(log_file),
  shQuote(smc_file),
  shQuote(out_bed_file)
)

status <- system(cmd)
if (status != 0) {
  stop(sprintf("smc2bed failed for %s", smc_file))
}

status <- system(sprintf("tabix -f %s", shQuote(out_bed_file)))
if (status != 0) {
  stop(sprintf("tabix failed for %s", out_bed_file))
}