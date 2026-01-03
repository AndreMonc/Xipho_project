#!/usr/bin/env Rscript

## ============================================================================
## FULL PIPELINE FOR SUMMARIZING FST PEAKS (FOR XIPHORHYNCHUS MANUSCRIPT)
##
## Includes:
## 1) Filter master dataframe
## 2) Fst outlier windows (5 SD above mean), autosomes vs Z separately
## 3) Cluster outliers into peaks (<=50 kb gaps)
## 4) ARG significance counts (windows + peaks) + patterned barplots
## 5) RAiSD bed import + 99th percentile thresholds (auto vs Z separately)
## 6) Peak model assignment:
##      - Selection-bottleneck requires XOR signals (one pop only) for enrichment/RTH
##      - Genomic architecture reflects SHARED signals (both pops) for enrichment/RTH
##      - If a peak has support for >1 model => "Overlapping models"
## 7) Heatmap (percentile-based; auto vs Z percentiles separately)
##      - Y-axis labeled by MODEL blocks only + BLANK SPACE between blocks
## 8) Per-peak Manhattan panel PDFs (15 stats) saved into folders by model block.
##      - 15 Mb span; peak centered in window; rolling smooth line included
## ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(rlang)
  library(ggpattern)
  library(slider)
  library(data.table)   # NEW (fast overlap joins)
})

## OPTIONAL: set working directory ------------------------------------------
setwd("/Users/andremoncrieff/Dropbox/Work/Postdoc/Manuscript--Xiphorhynchus/new_paradigm")

## USER INPUTS --------------------------------------------------------------
infile_master      <- "master_dataframe.txt"
infile_ARG_thresh  <- "ARG_thresholds_p0.001.txt"
infile_raisd_belem <- "combined_RAiSD_output.sorted.belem.bed"
infile_raisd_xingu <- "combined_RAiSD_output.sorted.xingu.bed"

outfile_filtered        <- "master_dataframe_filtered.txt"
outfile_auto_out        <- "Fst_xingu_belem_autosome_outliers_5SD.txt"
outfile_Z_out           <- "Fst_xingu_belem_Z_outliers_5SD.txt"
outfile_peaks           <- "Fst_Fst_outlier_peaks_50kb.txt"
outfile_barplot_windows <- "ARG_significant_Fst_outliers_barplot.pdf"
outfile_barplot_peaks   <- "ARG_significant_Fst_peaks_barplot.pdf"

outfile_heat_table      <- "Fst_peak_heatmap_percentiles_table.txt"
outfile_heatmap         <- "Fst_peak_heatmap_percentiles.pdf"

outfile_model_table     <- "Fst_peak_model_assignments.txt"

## NEW OUTPUT ROOT for per-peak Manhattan panels ----------------------------
outdir_peak_panels_root <- "Fst_peak_manhattan_panels"  # folders created under this

## PARAMETERS ---------------------------------------------------------------
chrom_Z_name <- "Chromosome_Z_RagTag"
peak_gap_bp  <- 50000   # 50 kb

## ARG statistics to consider (order on x-axis) -----------------------------
arg_stats_order <- c(
  "Xingu_enrichment",
  "Belem_enrichment",
  "Xingu_RTH",
  "Belem_RTH",
  "Xingu_Belem_RTH",
  "Tapajos_Xingu_RTH",
  "Tapajos_Xingu_JCR"
)

## STATS TO USE IN HEAT MAP (from master dataframe) -------------------------
heat_stats_master <- c(
  "Fst_tapajos_belem",
  "Fst_xingu_belem",
  "Fst_tapajos_xingu",
  "dxy_tapajos_belem",
  "dxy_xingu_belem",
  "dxy_tapajos_xingu",
  "pi_tapajos",
  "pi_xingu",
  "pi_belem",
  "Tapajos_enrichment",
  "Xingu_enrichment",
  "Belem_enrichment",
  "Tapajos_RTH",
  "Xingu_RTH",
  "Belem_RTH",
  "Xingu_Belem_RTH",
  "Tapajos_Xingu_RTH",
  "Tapajos_Xingu_JCR",
  "recombination_rate_50site_minimum"
)

## RAiSD stats to add to the heat map ---------------------------------------
heat_stats_raisd <- c("RAiSD_u_belem", "RAiSD_u_xingu")

## Combined heat-map stat order on x-axis -----------------------------------
heat_stats_order <- c(heat_stats_master, heat_stats_raisd)

## Model order TOP -> BOTTOM for heatmap blocks -----------------------------
model_block_order <- c(
  "Selection-bottleneck: Selection",
  "Selection-bottleneck: Bottleneck",
  "Selection-bottleneck: Other",
  "Genomic architecture",
  "Introgression",
  "Deep lineage sorting",
  "Overlapping models",
  "Unassigned to model"
)

## =============================================================================
## HELPERS
## =============================================================================
stop_if_missing_cols <- function(df, cols, df_name = "dataframe") {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop("Missing required columns in ", df_name, ": ", paste(missing, collapse = ", "))
  }
}

compute_is_sig_arg_vec <- function(chrom_vec, x_vec, thr_auto, thr_Z, thr_type, chrom_Z_name) {
  isZ <- chrom_vec == chrom_Z_name
  out <- rep(FALSE, length(x_vec))
  ok  <- !is.na(x_vec)
  
  if (thr_type == "high") {
    out[ ok & !isZ ] <- x_vec[ ok & !isZ ] > thr_auto
    out[ ok &  isZ ] <- x_vec[ ok &  isZ ] > thr_Z
  } else if (thr_type == "low") {
    out[ ok & !isZ ] <- x_vec[ ok & !isZ ] < thr_auto
    out[ ok &  isZ ] <- x_vec[ ok &  isZ ] < thr_Z
  } else {
    warning("Unknown threshold_type = ", thr_type, " (expected 'high' or 'low'). Returning FALSE.")
  }
  out
}

safe_percent_rank <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  if (sd(x, na.rm = TRUE) == 0) return(rep(0.5, length(x)))  # constant -> middle
  percent_rank(x)
}

## Chrom ordering: numbered chromosomes, then Z, then scaffolds
chrom_order_df <- function(chrom_vec, chrom_Z_name) {
  tibble(chrom = chrom_vec) %>%
    mutate(
      chrom_group = case_when(
        chrom == chrom_Z_name ~ "Z",
        str_detect(chrom, "^Chromosome_[0-9]+_RagTag$") ~ "num",
        TRUE ~ "scaf"
      ),
      chrom_num = if_else(
        chrom_group == "num",
        suppressWarnings(as.numeric(str_extract(chrom, "[0-9]+"))),
        NA_real_
      )
    )
}

## Safe directory/file name
sanitize_name <- function(x) {
  x %>%
    str_replace_all("[/\\\\:;\\*\\?\\\"\\<\\>\\|\\s]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

## Prefer an existing column; otherwise compute inverse = 1/x (useful fallback)
get_or_inverse <- function(df, col_name, fallback_raw = NULL, new_name = col_name) {
  if (col_name %in% names(df)) {
    df[[new_name]] <- as.numeric(df[[col_name]])
    return(df)
  }
  if (!is.null(fallback_raw) && fallback_raw %in% names(df)) {
    x <- as.numeric(df[[fallback_raw]])
    df[[new_name]] <- ifelse(is.na(x) | x == 0, NA_real_, 1 / x)
    return(df)
  }
  stop("Missing required column '", col_name, "' (and no fallback '", fallback_raw, "') in master dataframe.")
}

## =============================================================================
## STEP 1: READ & FILTER MASTER DATA
## =============================================================================
dat_raw <- readr::read_tsv(infile_master, show_col_types = FALSE)

stop_if_missing_cols(
  dat_raw,
  c("chrom", "start", "end", "Fst_xingu_belem", "sites_summ_stats", "Belem_RTH"),
  "master_dataframe"
)

## Also ensure heatmap stats exist (fail early)
stop_if_missing_cols(dat_raw, heat_stats_master, "master_dataframe (heat_stats_master)")

dat <- dat_raw %>%
  filter(!is.na(Fst_xingu_belem)) %>%
  mutate(Fst_xingu_belem = if_else(Fst_xingu_belem < 0, 0, Fst_xingu_belem)) %>%
  filter(sites_summ_stats >= 50) %>%
  filter(!is.na(Belem_RTH), is.finite(Belem_RTH)) %>%
  mutate(row_id = row_number())

write_tsv(dat, outfile_filtered)
cat("Number of rows after filtering:", nrow(dat), "\n\n")

## =============================================================================
## STEP 2: OUTLIER WINDOWS (5 SD ABOVE MEAN) AUTOSOMES VS Z
## =============================================================================
dat_auto <- dat %>% filter(chrom != chrom_Z_name)
dat_Z    <- dat %>% filter(chrom == chrom_Z_name)

auto_mean <- mean(dat_auto$Fst_xingu_belem)
auto_sd   <- sd(dat_auto$Fst_xingu_belem)
auto_thr  <- auto_mean + 5 * auto_sd
auto_outliers <- dat_auto %>% filter(Fst_xingu_belem >= auto_thr)

cat("AUTOSOMES (chrom != ", chrom_Z_name, ")\n", sep = "")
cat("  Mean:", auto_mean, "\n")
cat("  SD  :", auto_sd, "\n")
cat("  Thr :", auto_thr, "\n")
cat("  Outliers:", nrow(auto_outliers), "\n\n")

if (nrow(dat_Z) > 0) {
  Z_mean <- mean(dat_Z$Fst_xingu_belem)
  Z_sd   <- sd(dat_Z$Fst_xingu_belem)
  Z_thr  <- Z_mean + 5 * Z_sd
  Z_outliers <- dat_Z %>% filter(Fst_xingu_belem >= Z_thr)
  
  cat("Z CHROMOSOME (chrom == ", chrom_Z_name, ")\n", sep = "")
  cat("  Mean:", Z_mean, "\n")
  cat("  SD  :", Z_sd, "\n")
  cat("  Thr :", Z_thr, "\n")
  cat("  Outliers:", nrow(Z_outliers), "\n\n")
} else {
  Z_outliers <- dat_Z[0, ]
  cat("No rows found for ", chrom_Z_name, ".\n\n", sep = "")
}

write_tsv(auto_outliers, outfile_auto_out)
write_tsv(Z_outliers,    outfile_Z_out)

fst_outliers <- bind_rows(auto_outliers, Z_outliers)
total_fst_outliers <- nrow(fst_outliers)
cat("Total Fst outlier windows:", total_fst_outliers, "\n\n")
if (total_fst_outliers == 0) stop("No Fst outlier windows found; stopping.")

## =============================================================================
## STEP 3: CLUSTER OUTLIERS INTO PEAKS (<=50 kb gaps)
## =============================================================================
fst_outliers_with_peaks <- fst_outliers %>%
  arrange(chrom, start) %>%
  group_by(chrom) %>%
  mutate(
    gap_from_prev = start - lag(end),
    new_peak_flag = if_else(is.na(gap_from_prev) | gap_from_prev > peak_gap_bp, 1L, 0L),
    peak_id       = cumsum(new_peak_flag)
  ) %>%
  ungroup()

fst_peaks <- fst_outliers_with_peaks %>%
  group_by(chrom, peak_id) %>%
  summarise(
    peak_start = min(start),
    peak_end   = max(end),
    n_windows  = n(),
    .groups    = "drop"
  ) %>%
  arrange(chrom, peak_start)

n_peaks <- nrow(fst_peaks)
cat("Number of Fst peaks:", n_peaks, "\n\n")
write_tsv(fst_peaks, outfile_peaks)

## =============================================================================
## STEP 4: READ ARG THRESHOLDS
## =============================================================================
arg_thresh <- readr::read_tsv(infile_ARG_thresh, show_col_types = FALSE)
stop_if_missing_cols(arg_thresh, c("ARG_stats","autosome_thresholds","Z_chromosome_thresholds","threshold_type"), "ARG_thresholds")

arg_thresh_sub <- arg_thresh %>%
  filter(ARG_stats %in% arg_stats_order) %>%
  mutate(ARG_stats = factor(ARG_stats, levels = arg_stats_order)) %>%
  arrange(ARG_stats)

## =============================================================================
## STEP 5: WINDOW-LEVEL SIGNIFICANCE FLAGS FOR 7 ARG STATS
## =============================================================================
out_df <- fst_outliers_with_peaks
for (stat_name in arg_stats_order) {
  thr_row <- arg_thresh_sub %>% filter(ARG_stats == stat_name)
  if (nrow(thr_row) != 1) stop("Threshold row missing/duplicated for ARG stat: ", stat_name)
  
  thr_auto <- thr_row$autosome_thresholds[[1]]
  thr_Z    <- thr_row$Z_chromosome_thresholds[[1]]
  thr_type <- thr_row$threshold_type[[1]]
  
  sig_col <- paste0("sig_", stat_name)
  out_df[[sig_col]] <- compute_is_sig_arg_vec(
    chrom_vec = out_df$chrom,
    x_vec     = out_df[[stat_name]],
    thr_auto  = thr_auto,
    thr_Z     = thr_Z,
    thr_type  = thr_type,
    chrom_Z_name = chrom_Z_name
  )
}
fst_outliers_with_peaks <- out_df

## =============================================================================
## STEP 6: SIGNIFICANCE COUNTS (WINDOWS + PEAKS) + BARPLOTS
## =============================================================================
arg_label_levels <- paste0(
  as.character(arg_thresh_sub$ARG_stats),
  " (",
  arg_thresh_sub$threshold_type,
  ")"
)

arg_results <- map_dfr(arg_stats_order, function(stat_name) {
  sig_col <- paste0("sig_", stat_name)
  n_sig <- sum(fst_outliers_with_peaks[[sig_col]], na.rm = TRUE)
  tibble(ARG_stat = stat_name, n_sig = n_sig, percent = 100 * n_sig / total_fst_outliers)
}) %>%
  left_join(
    arg_thresh_sub %>% mutate(ARG_stat = as.character(ARG_stats)) %>% select(ARG_stat, threshold_type),
    by = "ARG_stat"
  ) %>%
  mutate(
    ARG_stat  = factor(ARG_stat, levels = arg_stats_order),
    ARG_label = paste0(ARG_stat, " (", threshold_type, ")"),
    ARG_label = factor(ARG_label, levels = arg_label_levels)
  )

arg_peak_results <- map_dfr(arg_stats_order, function(stat_name) {
  sig_col <- paste0("sig_", stat_name)
  n_sig_peaks <- fst_outliers_with_peaks %>%
    filter(.data[[sig_col]]) %>%
    distinct(chrom, peak_id) %>%
    nrow()
  tibble(ARG_stat = stat_name, n_sig = n_sig_peaks, percent = 100 * n_sig_peaks / n_peaks)
}) %>%
  left_join(arg_results %>% select(ARG_stat, threshold_type, ARG_label), by = "ARG_stat") %>%
  mutate(
    ARG_stat  = factor(ARG_stat, levels = arg_stats_order),
    ARG_label = factor(ARG_label, levels = arg_label_levels)
  )

prep_bar_data <- function(df) {
  df %>%
    mutate(
      fill_color = case_when(
        ARG_stat %in% c("Xingu_enrichment", "Xingu_RTH") ~ "#6ea8db",
        ARG_stat %in% c("Belem_enrichment", "Belem_RTH") ~ "#fdbb4a",
        ARG_stat == "Xingu_Belem_RTH" ~ "#6ea8db",
        ARG_stat %in% c("Tapajos_Xingu_RTH", "Tapajos_Xingu_JCR") ~ "#89ca8e",
        TRUE ~ "grey70"
      ),
      pattern_fill_color = case_when(
        ARG_stat == "Xingu_Belem_RTH" ~ "#fdbb4a",
        ARG_stat %in% c("Tapajos_Xingu_RTH", "Tapajos_Xingu_JCR") ~ "#6ea8db",
        TRUE ~ NA_character_
      ),
      pattern = if_else(is.na(pattern_fill_color), "none", "stripe"),
      label_text = paste0(n_sig, "\n(", sprintf("%.1f", percent), "%)")
    )
}

bar_data_windows <- prep_bar_data(arg_results)
bar_data_peaks   <- prep_bar_data(arg_peak_results)

plot_bars <- function(bar_df, ylab) {
  ggplot(
    bar_df,
    aes(
      x = ARG_label,
      y = n_sig,
      fill = fill_color,
      pattern = pattern,
      pattern_fill = pattern_fill_color
    )
  ) +
    geom_col_pattern(
      color = "black",
      pattern_angle = 45,
      pattern_density = 0.5,
      pattern_spacing = 0.02,
      pattern_color = NA,
      na.rm = TRUE
    ) +
    scale_fill_identity() +
    scale_pattern_identity() +
    scale_pattern_fill_identity() +
    geom_text(aes(label = label_text), vjust = -0.3, lineheight = 0.9, size = 3) +
    labs(x = NULL, y = ylab) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    expand_limits(y = max(bar_df$n_sig, na.rm = TRUE) * 1.2)
}

p_windows <- plot_bars(bar_data_windows, "Fst outlier windows with significant ARG statistics")
pdf(outfile_barplot_windows, width = 7, height = 5); print(p_windows); dev.off()

p_peaks <- plot_bars(bar_data_peaks, "Fst peaks with significant ARG statistics")
pdf(outfile_barplot_peaks, width = 7, height = 5); print(p_peaks); dev.off()

cat("Window-based barplot saved to:", outfile_barplot_windows, "\n")
cat("Peak-based barplot saved to:", outfile_barplot_peaks, "\n")
cat("Fst peaks saved to:", outfile_peaks, "\n\n")

## =============================================================================
## STEP 7: RAiSD IMPORT + 99th PERCENTILE THRESHOLDS (AUTO vs Z) + PEAK OVERLAP
## =============================================================================
read_raisd <- function(infile) {
  read_tsv(infile, col_names = c("chrom","start","end","u"), show_col_types = FALSE) %>%
    mutate(chrom_type = if_else(chrom == chrom_Z_name, "Z", "auto"))
}

raisd_belem <- read_raisd(infile_raisd_belem)
raisd_xingu <- read_raisd(infile_raisd_xingu)

raisd_thr <- function(df) {
  df %>%
    group_by(chrom_type) %>%
    summarise(u_p99 = quantile(u, probs = 0.99, na.rm = TRUE), .groups = "drop")
}

thr_belem <- raisd_thr(raisd_belem)
thr_xingu <- raisd_thr(raisd_xingu)

peak_has_raisd_p99 <- function(fst_peaks, raisd_df, thr_df) {
  if (nrow(fst_peaks) == 0) return(tibble(chrom = character(), peak_id = integer(), has_p99 = logical()))
  
  map_dfr(seq_len(nrow(fst_peaks)), function(i) {
    p_chr <- fst_peaks$chrom[i]
    p_id  <- fst_peaks$peak_id[i]
    p_s   <- fst_peaks$peak_start[i]
    p_e   <- fst_peaks$peak_end[i]
    
    chrom_type <- if_else(p_chr == chrom_Z_name, "Z", "auto")
    u_thr <- thr_df %>% filter(chrom_type == !!chrom_type) %>% pull(u_p99)
    if (length(u_thr) == 0 || is.na(u_thr)) u_thr <- Inf
    
    hits <- raisd_df %>% filter(chrom == p_chr, end >= p_s, start <= p_e)
    
    tibble(
      chrom = p_chr,
      peak_id = p_id,
      has_p99 = if (nrow(hits) == 0) FALSE else any(hits$u > u_thr, na.rm = TRUE)
    )
  })
}

peak_p99_belem <- peak_has_raisd_p99(fst_peaks, raisd_belem, thr_belem) %>% rename(has_p99_belem = has_p99)
peak_p99_xingu <- peak_has_raisd_p99(fst_peaks, raisd_xingu, thr_xingu) %>% rename(has_p99_xingu = has_p99)

peak_raisd_p99 <- fst_peaks %>%
  select(chrom, peak_id) %>%
  left_join(peak_p99_belem, by = c("chrom","peak_id")) %>%
  left_join(peak_p99_xingu, by = c("chrom","peak_id")) %>%
  mutate(
    has_p99_belem = replace_na(has_p99_belem, FALSE),
    has_p99_xingu = replace_na(has_p99_xingu, FALSE),
    has_p99_any   = has_p99_belem | has_p99_xingu
  )

## =============================================================================
## STEP 8: MODEL ASSIGNMENT
## =============================================================================
w <- fst_outliers_with_peaks

sig_Xen  <- w[["sig_Xingu_enrichment"]]
sig_Ben  <- w[["sig_Belem_enrichment"]]
sig_Xrth <- w[["sig_Xingu_RTH"]]
sig_Brth <- w[["sig_Belem_RTH"]]
sig_XB   <- w[["sig_Xingu_Belem_RTH"]]
sig_TXR  <- w[["sig_Tapajos_Xingu_RTH"]]
sig_TXJ  <- w[["sig_Tapajos_Xingu_JCR"]]

## Selection-bottleneck: EXACTLY ONE population significant (XOR), not both.
enrichment_xor <- (sig_Xen & !sig_Ben) | (!sig_Xen & sig_Ben)
rth_xor        <- (sig_Xrth & !sig_Brth) | (!sig_Xrth & sig_Brth)

## Genomic architecture: SHARED signals (both pops)
enrichment_both <- sig_Xen & sig_Ben
rth_both        <- sig_Xrth & sig_Brth

w <- w %>%
  mutate(
    support_introgression = sig_TXR | sig_TXJ,
    support_dls           = sig_XB,
    support_selbot        = enrichment_xor | rth_xor,
    support_genoarch      = enrichment_both | rth_both
  )

peak_support <- w %>%
  group_by(chrom, peak_id) %>%
  summarise(
    outlier_n_windows = n(),
    any_intro  = any(support_introgression, na.rm = TRUE),
    any_dls    = any(support_dls, na.rm = TRUE),
    any_selbot = any(support_selbot, na.rm = TRUE),
    any_geno   = any(support_genoarch, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    fst_peaks %>% select(chrom, peak_id, peak_start, peak_end, peak_n_windows = n_windows),
    by = c("chrom","peak_id")
  ) %>%
  left_join(peak_raisd_p99, by = c("chrom","peak_id")) %>%
  mutate(
    n_models_supported = (any_intro + any_dls + any_selbot + any_geno),
    primary_model = case_when(
      n_models_supported >= 2 ~ "Overlapping models",
      n_models_supported == 0 ~ "Unassigned to model",
      any_selbot ~ "Selection-bottleneck",
      any_geno   ~ "Genomic architecture",
      any_intro  ~ "Introgression",
      any_dls    ~ "Deep lineage sorting",
      TRUE ~ "Unassigned to model"
    )
  )

## Selection-bottleneck subcategory:
## - "Selection": >1 outlier window AND any RAiSD p99 (either pop)
## - "Bottleneck": exactly 1 outlier window AND no RAiSD p99
## - otherwise "Other"
peak_support <- peak_support %>%
  mutate(
    model_block = case_when(
      primary_model != "Selection-bottleneck" ~ primary_model,
      outlier_n_windows > 1 & has_p99_any ~ "Selection-bottleneck: Selection",
      outlier_n_windows == 1 & !has_p99_any ~ "Selection-bottleneck: Bottleneck",
      TRUE ~ "Selection-bottleneck: Other"
    ),
    model_block = factor(model_block, levels = model_block_order)
  )

## NEW: which model blocks are overlapping (peak-level)
peak_support <- peak_support %>%
  mutate(
    overlap_models = pmap_chr(
      list(any_selbot, any_geno, any_intro, any_dls),
      function(selbot, geno, intro, dls) {
        mods <- c(
          if (isTRUE(selbot)) "Selection-bottleneck" else NULL,
          if (isTRUE(geno))   "Genomic architecture" else NULL,
          if (isTRUE(intro))  "Introgression" else NULL,
          if (isTRUE(dls))    "Deep lineage sorting" else NULL
        )
        if (length(mods) == 0) NA_character_ else paste(mods, collapse = " + ")
      }
    )
  )


model_table <- peak_support %>%
  select(
    chrom, peak_id, peak_start, peak_end,
    outlier_n_windows, peak_n_windows,
    has_p99_belem, has_p99_xingu, has_p99_any,
    any_intro, any_dls, any_selbot, any_geno,
    overlap_models,
    primary_model, model_block
  ) %>%
  arrange(model_block, chrom, peak_start)

write_tsv(model_table, outfile_model_table)
cat("Model assignment table saved to:", outfile_model_table, "\n\n")
cat("MODEL BLOCK COUNTS:\n")
print(model_table %>% count(model_block, name = "n_peaks") %>% arrange(model_block))
cat("\n")

## =============================================================================
## STEP 8.5: WHY "OVERLAPPING MODELS"?  (peak-level list of sig ARG stats)
## =============================================================================

## Make a lookup for threshold direction (high/low) for each ARG stat
arg_thr_lookup <- arg_thresh_sub %>%
  mutate(
    stat = as.character(ARG_stats),
    thr_type = as.character(threshold_type)
  ) %>%
  select(stat, thr_type)

## Summarize which sig_* flags occur anywhere within each peak
peak_sig_long <- fst_outliers_with_peaks %>%
  select(chrom, peak_id, starts_with("sig_")) %>%
  group_by(chrom, peak_id) %>%
  summarise(across(starts_with("sig_"), ~ any(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(
    cols = starts_with("sig_"),
    names_to = "sig_col",
    values_to = "is_sig"
  ) %>%
  filter(is_sig) %>%
  mutate(stat = str_replace(sig_col, "^sig_", "")) %>%
  left_join(arg_thr_lookup, by = "stat") %>%
  mutate(
    thr_type = replace_na(thr_type, "NA"),
    stat_label = paste0(stat, " (", thr_type, ")")
  )

## Collapse to one reason string per peak
peak_overlap_reasons <- peak_sig_long %>%
  group_by(chrom, peak_id) %>%
  summarise(
    overlap_reason = paste(stat_label, collapse = "; "),
    .groups = "drop"
  )

## Add it onto model_table (NA for non-overlap peaks is fine)
model_table <- model_table %>%
  left_join(peak_overlap_reasons, by = c("chrom", "peak_id"))


## =============================================================================
## STEP 9: HEAT MAP PREP – MASTER STATS PERCENTILES (AUTO vs Z separately)
## =============================================================================
master_long_all <- dat %>%
  mutate(chrom_type = if_else(chrom == chrom_Z_name, "Z", "auto")) %>%
  select(row_id, chrom, start, end, chrom_type, all_of(heat_stats_master)) %>%
  pivot_longer(cols = all_of(heat_stats_master), names_to = "stat", values_to = "value") %>%
  group_by(chrom_type, stat) %>%
  mutate(percentile = safe_percent_rank(value)) %>%
  ungroup()

master_long_peaks <- master_long_all %>%
  inner_join(
    fst_outliers_with_peaks %>% select(row_id, chrom, peak_id),
    by = c("row_id", "chrom")
  )

peak_master_long <- master_long_peaks %>%
  group_by(chrom, peak_id, stat) %>%
  summarise(
    stat_percentile = {
      vals <- percentile[!is.na(percentile)]
      if (length(vals) == 0) NA_real_ else vals[which.max(abs(vals - 0.5))]
    },
    .groups = "drop"
  )

## =============================================================================
## STEP 10: RAiSD WINDOW-MEAN PERCENTILES (AUTO vs Z separately)
##   - compute mean(u) per filtered window in dat
##   - compute percentiles over window means
##   - for peaks: choose the outlier window percentile farthest from 0.5
## =============================================================================

## Helper: mean RAiSD u per dat window using fast overlap joins (data.table::foverlaps)
window_mean_raisd <- function(dat_windows, raisd_df, chrom_Z_name) {
  # dat_windows must have: row_id, chrom, start, end
  wdt <- as.data.table(dat_windows %>%
                         transmute(row_id, chrom, start, end,
                                   chrom_type = if_else(chrom == chrom_Z_name, "Z", "auto"))
  )
  rdt <- as.data.table(raisd_df %>%
                         transmute(chrom, start, end, u = as.numeric(u),
                                   chrom_type = if_else(chrom == chrom_Z_name, "Z", "auto"))
  )
  
  setkey(wdt, chrom, start, end)
  setkey(rdt, chrom, start, end)
  
  # overlap RAiSD intervals onto windows
  ov <- foverlaps(rdt, wdt, by.x = c("chrom","start","end"),
                  by.y = c("chrom","start","end"),
                  type = "any", nomatch = NA)
  
  # mean u per window (row_id)
  out <- ov[!is.na(row_id),
            .(raisd_mean_u = mean(u, na.rm = TRUE)),
            by = .(row_id, chrom, chrom_type)]
  
  # ensure every window exists (NA if no overlaps)
  out_full <- merge(
    wdt[, .(row_id, chrom, chrom_type)],
    out,
    by = c("row_id","chrom","chrom_type"),
    all.x = TRUE
  )
  
  as_tibble(out_full)
}

## Compute window means for each population
raisd_belem_win <- window_mean_raisd(dat, raisd_belem, chrom_Z_name) %>%
  rename(raisd_mean_u_belem = raisd_mean_u)

raisd_xingu_win <- window_mean_raisd(dat, raisd_xingu, chrom_Z_name) %>%
  rename(raisd_mean_u_xingu = raisd_mean_u)

## Compute percentiles over WINDOW MEANS (auto vs Z separately)
raisd_belem_win_p <- raisd_belem_win %>%
  group_by(chrom_type) %>%
  mutate(raisd_pctl_belem = safe_percent_rank(raisd_mean_u_belem)) %>%
  ungroup()

raisd_xingu_win_p <- raisd_xingu_win %>%
  group_by(chrom_type) %>%
  mutate(raisd_pctl_xingu = safe_percent_rank(raisd_mean_u_xingu)) %>%
  ungroup()

## Attach peak_id ONLY for Fst outlier windows (like the master stats)
raisd_belem_outlier <- fst_outliers_with_peaks %>%
  select(row_id, chrom, peak_id) %>%
  left_join(raisd_belem_win_p %>% select(row_id, chrom, raisd_pctl_belem),
            by = c("row_id","chrom"))

raisd_xingu_outlier <- fst_outliers_with_peaks %>%
  select(row_id, chrom, peak_id) %>%
  left_join(raisd_xingu_win_p %>% select(row_id, chrom, raisd_pctl_xingu),
            by = c("row_id","chrom"))

## For each peak: choose percentile farthest from 0.5 (same rule as master stats)
raisd_belem_peaks <- raisd_belem_outlier %>%
  group_by(chrom, peak_id) %>%
  summarise(
    stat = "RAiSD_u_belem",
    stat_percentile = {
      vals <- raisd_pctl_belem[!is.na(raisd_pctl_belem)]
      if (length(vals) == 0) NA_real_ else vals[which.max(abs(vals - 0.5))]
    },
    .groups = "drop"
  )

raisd_xingu_peaks <- raisd_xingu_outlier %>%
  group_by(chrom, peak_id) %>%
  summarise(
    stat = "RAiSD_u_xingu",
    stat_percentile = {
      vals <- raisd_pctl_xingu[!is.na(raisd_pctl_xingu)]
      if (length(vals) == 0) NA_real_ else vals[which.max(abs(vals - 0.5))]
    },
    .groups = "drop"
  )


## =============================================================================
## STEP 11: COMBINE HEAT DATA + SAVE WIDE TABLE
## =============================================================================
heat_peak_long <- bind_rows(
  peak_master_long,
  raisd_belem_peaks,
  raisd_xingu_peaks
) %>%
  left_join(fst_peaks, by = c("chrom", "peak_id")) %>%
  left_join(model_table %>% select(chrom, peak_id, model_block), by = c("chrom","peak_id"))

heat_table_wide <- heat_peak_long %>%
  select(chrom, peak_id, peak_start, peak_end, model_block, stat, stat_percentile) %>%
  distinct() %>%
  pivot_wider(names_from = stat, values_from = stat_percentile) %>%
  left_join(chrom_order_df(unique(.$chrom), chrom_Z_name), by = "chrom") %>%
  arrange(model_block, factor(chrom_group, levels = c("num","Z","scaf")), chrom_num, chrom, peak_start) %>%
  select(-chrom_group, -chrom_num)

write_tsv(heat_table_wide, outfile_heat_table)
cat("Heat-map percentile table saved to:", outfile_heat_table, "\n\n")

## =============================================================================
## STEP 12: HEATMAP PLOT (MODEL BLOCK LABELS ONLY + SPACER ROWS)
## =============================================================================
peaks_for_plot <- fst_peaks %>%
  left_join(model_table %>% select(chrom, peak_id, model_block), by = c("chrom","peak_id")) %>%
  left_join(chrom_order_df(unique(fst_peaks$chrom), chrom_Z_name), by = "chrom") %>%
  mutate(
    model_block = factor(model_block, levels = model_block_order),
    peak_label  = paste0(chrom, ":", peak_start, "-", peak_end)
  ) %>%
  arrange(
    model_block,
    factor(chrom_group, levels = c("num","Z","scaf")),
    chrom_num,
    chrom,
    peak_start
  ) %>%
  select(chrom, peak_id, peak_start, peak_end, model_block, peak_label)

y_levels <- c()
y_label_map <- c()

for (mb in model_block_order) {
  block_peaks <- peaks_for_plot %>% filter(model_block == mb)
  
  if (nrow(block_peaks) > 0) {
    y_levels <- c(y_levels, block_peaks$peak_label)
    
    y_label_map[block_peaks$peak_label[1]] <- as.character(mb)
    if (nrow(block_peaks) > 1) y_label_map[block_peaks$peak_label[-1]] <- ""
    
    spacer_key <- paste0("SPACER__", sanitize_name(mb))
    y_levels <- c(y_levels, spacer_key)
    y_label_map[spacer_key] <- ""
  } else {
    spacer_key <- paste0("SPACER__", sanitize_name(mb))
    y_levels <- c(y_levels, spacer_key)
    y_label_map[spacer_key] <- ""
  }
}

heat_peak_long2 <- heat_peak_long %>%
  mutate(
    peak_label = paste0(chrom, ":", peak_start, "-", peak_end),
    stat       = factor(stat, levels = heat_stats_order)
  ) %>%
  filter(peak_label %in% peaks_for_plot$peak_label)

## FIXED: spacer rows built with expand_grid() to match sizes
spacer_labels <- y_levels[str_detect(y_levels, "^SPACER__")]
spacer_rows <- tidyr::expand_grid(
  peak_label = spacer_labels,
  stat       = factor(heat_stats_order, levels = heat_stats_order)
) %>%
  mutate(stat_percentile = NA_real_)

heat_plot_df <- heat_peak_long2 %>%
  select(peak_label, stat, stat_percentile) %>%
  bind_rows(spacer_rows) %>%
  mutate(
    peak_label = factor(peak_label, levels = y_levels),
    stat = factor(stat, levels = heat_stats_order)
  )

p_heat <- ggplot(heat_plot_df, aes(x = stat, y = peak_label, fill = stat_percentile)) +
  geom_tile(color = NA) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0.5, limits = c(0, 1), na.value = "white"
  ) +
  scale_y_discrete(
    limits = rev(y_levels),                    # top-to-bottom in requested block order
    labels = function(x) unname(y_label_map[x])
  ) +
  labs(x = NULL, y = NULL, fill = "Percentile") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(outfile_heatmap, p_heat, width = 12, height = 10, units = "in", dpi = 300)
cat("Heatmap saved to:", outfile_heatmap, "\n\n")

## =============================================================================
## STEP 13 (NEW): PER-PEAK MANHATTAN PANEL PDFs (15 STATS) BY MODEL FOLDER
## =============================================================================

## Create output root + category folders
dir.create(outdir_peak_panels_root, showWarnings = FALSE, recursive = TRUE)
for (mb in model_block_order) {
  dir.create(file.path(outdir_peak_panels_root, sanitize_name(mb)), showWarnings = FALSE, recursive = TRUE)
}

## Prepare master data columns for Manhattan panels:
## We'll use pos = midpoint of [start,end], in Mb for plotting.
dat2 <- dat %>%
  mutate(
    pos_bp = (start + end) / 2,
    pos_Mb = pos_bp / 1e6
  )

## Ensure inverse RTH columns exist for plotting:
## If your master already has Xingu_RTH_inverse / Belem_RTH_inverse, it uses them;
## otherwise it computes inverse from Xingu_RTH / Belem_RTH.
dat2 <- get_or_inverse(dat2, "Xingu_RTH_inverse", fallback_raw = "Xingu_RTH", new_name = "Xingu_RTH_inverse")
dat2 <- get_or_inverse(dat2, "Belem_RTH_inverse", fallback_raw = "Belem_RTH", new_name = "Belem_RTH_inverse")

## Ensure additional Manhattan columns exist
manhattan_master_cols <- c(
  "Fst_xingu_belem",
  "dxy_xingu_belem",
  "pi_xingu",
  "pi_belem",
  "Xingu_enrichment",
  "Belem_enrichment",
  "Xingu_RTH_inverse",
  "Belem_RTH_inverse",
  "Xingu_Belem_RTH",
  "Tapajos_Xingu_RTH",
  "Tapajos_Xingu_JCR",
  "recombination_rate_50site_minimum",
  "avg_depth_from_bam_files"
)
stop_if_missing_cols(dat2, manhattan_master_cols, "master_dataframe (for Manhattan panels)")

## Read RAiSD and make midpoints
raisd_to_plot <- function(df, stat_label) {
  df %>%
    transmute(
      chrom = chrom,
      pos_bp = (start + end) / 2,
      pos_Mb = pos_bp / 1e6,
      stat = stat_label,
      value = as.numeric(u)
    )
}
raisd_xingu_plot <- raisd_to_plot(raisd_xingu, "U statistic Xingu (RAISD)")
raisd_belem_plot <- raisd_to_plot(raisd_belem, "U statistic Belem (RAISD)")

## Rolling smooth parameters (match your guide style) :contentReference[oaicite:1]{index=1}
k_before <- 5
k_after  <- 5

## Manhattan stat mapping + order (15 panels)
panel_levels <- c(
  "Fst_xingu_belem",
  "dxy_xingu_belem",
  "pi_xingu",
  "pi_belem",
  "Xingu_enrichment",
  "Belem_enrichment",
  "Xingu_RTH_inverse",
  "Belem_RTH_inverse",
  "Xingu_Belem_RTH",
  "Tapajos_Xingu_RTH",
  "Tapajos_Xingu_JCR",
  "recombination_rate_50site_minimum",
  "avg_depth_from_bam_files",
  "U statistic Xingu (RAISD)",
  "U statistic Belem (RAISD)"
)

panel_labels <- c(
  Fst_xingu_belem = "Xingu-Belem Fst",
  dxy_xingu_belem = "Xingu-Belem Dxy",
  pi_xingu = "Xingu Pi",
  pi_belem = "Belem Pi",
  Xingu_enrichment = "Xingu enrichment",
  Belem_enrichment = "Belem enrichment",
  Xingu_RTH_inverse = "Xingu RTH (inverse)",
  Belem_RTH_inverse = "Belem RTH (inverse)",
  Xingu_Belem_RTH = "Xingu-Belem RTH",
  Tapajos_Xingu_RTH = "Tapajos-Xingu RTH",
  Tapajos_Xingu_JCR = "Tapajos-Xingu JCR",
  recombination_rate_50site_minimum = "Recombination rate",
  avg_depth_from_bam_files = "Average depth",
  `U statistic Xingu (RAISD)` = "U statistic Xingu (RAiSD)",
  `U statistic Belem (RAISD)` = "U statistic Belem (RAiSD)"
)

## Function to build and save one peak PDF
save_peak_panel <- function(peak_row) {
  chr <- peak_row$chrom
  pstart <- peak_row$peak_start
  pend   <- peak_row$peak_end
  mb     <- as.character(peak_row$model_block)
  
  center_bp <- (pstart + pend) / 2
  half_window <- 7.5e6
  win_start <- max(0, center_bp - half_window)
  win_end   <- center_bp + half_window
  
  ## Master long for this region
  region_master <- dat2 %>%
    filter(chrom == chr, pos_bp >= win_start, pos_bp <= win_end) %>%
    transmute(
      chrom = chrom,
      pos_bp = pos_bp,
      pos_Mb = pos_Mb,
      `Fst_xingu_belem` = as.numeric(Fst_xingu_belem),
      `dxy_xingu_belem` = as.numeric(dxy_xingu_belem),
      `pi_xingu` = as.numeric(pi_xingu),
      `pi_belem` = as.numeric(pi_belem),
      `Xingu_enrichment` = as.numeric(Xingu_enrichment),
      `Belem_enrichment` = as.numeric(Belem_enrichment),
      `Xingu_RTH_inverse` = as.numeric(Xingu_RTH_inverse),
      `Belem_RTH_inverse` = as.numeric(Belem_RTH_inverse),
      `Xingu_Belem_RTH` = as.numeric(Xingu_Belem_RTH),
      `Tapajos_Xingu_RTH` = as.numeric(Tapajos_Xingu_RTH),
      `Tapajos_Xingu_JCR` = as.numeric(Tapajos_Xingu_JCR),
      `recombination_rate_50site_minimum` = as.numeric(recombination_rate_50site_minimum),
      `avg_depth_from_bam_files` = as.numeric(avg_depth_from_bam_files)
    ) %>%
    pivot_longer(cols = -c(chrom, pos_bp, pos_Mb), names_to = "stat", values_to = "value")
  
  ## RAiSD long for this region
  region_raisd <- bind_rows(
    raisd_xingu_plot,
    raisd_belem_plot
  ) %>%
    filter(chrom == chr, pos_bp >= win_start, pos_bp <= win_end)
  
  plot_df <- bind_rows(
    region_master,
    region_raisd
  ) %>%
    mutate(
      stat = factor(stat, levels = panel_levels)
    )
  
  if (nrow(plot_df) == 0) {
    message("No data for peak panel: ", chr, ":", pstart, "-", pend)
    return(invisible(NULL))
  }
  
  smooth_df <- plot_df %>%
    group_by(stat) %>%
    arrange(pos_bp, .by_group = TRUE) %>%
    mutate(
      smooth = slide_dbl(
        value,
        .f = ~ mean(.x, na.rm = TRUE),
        .before = k_before,
        .after = k_after,
        .complete = FALSE
      )
    ) %>%
    ungroup()
  
  ## Peak center line (Mb)
  center_Mb <- center_bp / 1e6
  
  ## --- Title/subtitle ---
  title_line1 <- paste0("Peak: ", chr, ":", pstart, "-", pend)
  title_line2 <- paste0("Model: ", mb)
  
  overlap_models <- peak_row$overlap_models
  overlap_reason <- peak_row$overlap_reason
  
  subtitle_text <- title_line2
  
  if (mb == "Overlapping models") {
    
    if (!is.null(overlap_models) && !is.na(overlap_models) && overlap_models != "") {
      subtitle_text <- paste0(
        subtitle_text,
        "\nOverlapping models: ",
        stringr::str_wrap(overlap_models, width = 110)
      )
    }
    
    if (!is.null(overlap_reason) && !is.na(overlap_reason) && overlap_reason != "") {
      subtitle_text <- paste0(
        subtitle_text,
        "\nOverlap signals: ",
        stringr::str_wrap(overlap_reason, width = 110)
      )
    }
  }
  
  p <- ggplot(smooth_df, aes(x = pos_Mb, y = value)) +
    geom_point(alpha = 0.6, size = 0.35) +
    geom_line(aes(y = smooth), linewidth = 0.5, colour = "red") +
    geom_vline(xintercept = center_Mb, linetype = "dashed", linewidth = 0.4) +
    scale_x_continuous(
      name = paste0("Position on ", chr, " (Mb)"),
      breaks = scales::pretty_breaks(n = 10),
      expand = c(0, 0)
    ) +
    facet_grid(stat ~ ., scales = "free_y", switch = "y",
               labeller = as_labeller(panel_labels)) +
    labs(
      y = NULL,
      title = title_line1,
      subtitle = subtitle_text
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1),
      strip.placement    = "outside",
      strip.background   = element_blank(),
      strip.text.y.left  = element_text(size = 9, angle = 0),
      strip.text.y.right = element_text(size = 9, angle = 0),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      
      ## key additions:
      plot.title.position = "plot",                # title uses full plot width
      plot.title = element_text(hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(hjust = 0, margin = margin(b = 6)),
      plot.margin = margin(t = 10, r = 30, b = 10, l = 10)  # more right/top breathing room
    )
  
  outdir <- file.path(outdir_peak_panels_root, sanitize_name(mb))
  fname  <- paste0(sanitize_name(chr), "__", pstart, "_", pend, ".pdf")
  outfile <- file.path(outdir, fname)
  
  ggsave(outfile, p, width = 8.5, height = 11, units = "in", dpi = 300)
  invisible(outfile)
}

## Iterate peaks and save
peaks_with_models <- fst_peaks %>%
  left_join(
    model_table %>% select(chrom, peak_id, model_block, overlap_models, overlap_reason),
    by = c("chrom","peak_id")
  ) %>%
  mutate(model_block = factor(model_block, levels = model_block_order)) %>%
  arrange(model_block, chrom, peak_start)

cat("Saving per-peak Manhattan panels into: ", outdir_peak_panels_root, "\n")
outfiles <- vector("character", nrow(peaks_with_models))
for (i in seq_len(nrow(peaks_with_models))) {
  outfiles[i] <- save_peak_panel(peaks_with_models[i, ])
}
cat("Done. Manhattan panels written (non-empty): ", sum(!is.na(outfiles) & outfiles != ""), "\n")
