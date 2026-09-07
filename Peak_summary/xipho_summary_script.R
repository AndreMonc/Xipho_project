#!/usr/bin/env Rscript

library(data.table)

if (!requireNamespace("ape", quietly = TRUE)) {
  warning("Package 'ape' is not installed; Manhattan panels will be generated without tree plots.")
}

if (!requireNamespace("pdftools", quietly = TRUE) &&
    !requireNamespace("magick", quietly = TRUE)) {
  warning("Install either 'pdftools' or 'magick' to embed manhattan_footnotes.pdf in Manhattan panels.")
}

# -----------------------------
# Paths
# -----------------------------
# The script is expected to live in:
#   Xipho_revision/Peak_summary
#
# Paths below are defined relative to the script location, so the script
# should work across computers as long as the Dropbox folder structure is
# the same. Run with, for example:
#   Rscript xipho_summary_script.R

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grep("^--file=", args)]

if (length(script_arg) > 0) {
  script_path <- sub("^--file=", "", script_arg[1])
  work_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
} else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  if (nzchar(script_path)) {
    work_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
  } else {
    work_dir <- normalizePath(getwd(), mustWork = TRUE)
    warning("Could not identify active script path; using current working directory as work_dir.")
  }
} else {
  work_dir <- normalizePath(getwd(), mustWork = TRUE)
  warning("Could not identify script path from commandArgs(); using current working directory as work_dir.")
}

setwd(work_dir)

# -----------------------------
# Plot creation switches
# -----------------------------
# Set these switches to TRUE only for the plot groups you want to regenerate.
CREATE_BARPLOT_PDFS <- FALSE

# All one-off PDFs written to output_files/figures, except the three peak-size
# distribution PDFs, which have their own grouped switch below.
CREATE_FIGURES_PDFS <- FALSE

# Controls all three Xingu-Belem Fst peak-size distribution PDFs:
#   1) uncolored all-peaks plot,
#   2) bars colored by the at-least-one-tree model assignments, and
#   3) bars colored by the >=25% model assignments.
CREATE_PEAK_SIZE_DISTRIBUTION_PDFS <- FALSE

# Peak-by-peak Manhattan PDFs: independent controls for each model-assignment
# method and for the 2000th-iteration versus multi-iteration plot sets.
CREATE_MANHATTAN_ANY_SINGLE_ITER_PDFS <- FALSE
CREATE_MANHATTAN_MIN25_SINGLE_ITER_PDFS <- FALSE
CREATE_MANHATTAN_ANY_MULTI_ITER_PDFS <- FALSE
CREATE_MANHATTAN_MIN25_MULTI_ITER_PDFS <- FALSE

# Controls the slow genome-wide summaries across all 50 saved ARG MCMC
# iterations (1510-2000 every 10 iterations), including their output tables
# and the 50-iteration model/statistic summaries. This does NOT control the
# independently switched multi-iteration Manhattan PDF generation.
RUN_MULTI_ITER_ARG_SUMMARIES <- FALSE


# -----------------------------
# Fst-peak vs. control-region significance subsampling
# -----------------------------
# Run matched-control subsampling for every genomic statistic that already has
# a peak-level TRUE/FALSE significance designation in the comprehensive tables.
RUN_FST_PEAK_CONTROL_SUBSAMPLING <- TRUE

# Number of independent null datasets. Increase this value later for more
# precise tail probabilities/critical values.
CONTROL_SUBSAMPLING_N <- 10000L

# Fixed seed makes the subsampling exactly reproducible.
CONTROL_SUBSAMPLING_SEED <- 24680L

pixy_file <- file.path(
  "..",
  "Fst_Dxy_Pi",
  "pixy_merged_fst_dxy_pi.chrom_type.tsv"
)

arg_stat_dir <- file.path(
  "..",
  "ARGweaver",
  "stat_files"
)

arg_stat_file_2000 <- file.path(
  arg_stat_dir,
  "argStats_midpoint.2000.stat.gz"
)

arg_cc_stat_dir <- file.path(
  "..",
  "ARGweaver",
  "stat_files_CC"
)

arg_cc_stat_file_2000 <- file.path(
  arg_cc_stat_dir,
  "argStats_midpoint_CC.2000.stat.gz"
)

output_dir <- file.path(work_dir, "output_files")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Purpose-specific output subfolders keep output_files organized.
fst_outlier_dir <- file.path(output_dir, "Fst_outlier_definitions")
significance_thresholds_dir <- file.path(output_dir, "Significance_thresholds")
arg_barplots_significance_dir <- file.path(output_dir, "ARG_barplots_and_significance")
arg_coverage_dir <- file.path(output_dir, "ARG_coverage")
arg_diagnostics_dir <- file.path(output_dir, "ARG_diagnostics")
heatmap_dir <- file.path(output_dir, "Heatmap")
figures_dir <- file.path(output_dir, "figures")
tables_dir <- file.path(output_dir, "tables")

for (dir_to_make in c(
  fst_outlier_dir,
  significance_thresholds_dir,
  arg_barplots_significance_dir,
  arg_coverage_dir,
  arg_diagnostics_dir,
  heatmap_dir,
  figures_dir,
  tables_dir
)) {
  dir.create(dir_to_make, showWarnings = FALSE, recursive = TRUE)
}

outlier_windows_bed <- file.path(
  fst_outlier_dir,
  "Xin_Bel_Fst.outlier_windows.autosomes.min50snps.bed"
)

outlier_peaks_bed <- file.path(
  fst_outlier_dir,
  "Xin_Bel_Fst.outlier_peaks.autosomes.min50snps.bed"
)

comparison_windows_file <- file.path(
  fst_outlier_dir,
  "Xin_Bel_Fst.outlier_vs_control_window_coordinates.tsv"
)

window_stat_summary_file <- file.path(
  fst_outlier_dir,
  "Xin_Bel_Fst.outlier_vs_control_stat_summary.tsv"
)

region_set_pi_recombination_summary_file <- file.path(
  fst_outlier_dir,
  "Xin_Bel_Fst.region_set_pi_recombination_rate_summary.tsv"
)

pixy_threshold_summary_file <- file.path(
  significance_thresholds_dir,
  "pixy_control_window_empirical_significance_thresholds.tsv"
)

arg_threshold_summary_file <- file.path(
  significance_thresholds_dir,
  "ARGweaver_control_window_empirical_significance_thresholds.iter2000.tsv"
)

arg_cc_threshold_summary_file <- file.path(
  significance_thresholds_dir,
  "ARGweaver_CC_control_window_empirical_significance_thresholds.iter2000.tsv"
)

arg_cc_stat_summary_file <- file.path(
  arg_diagnostics_dir,
  "ARGweaver_CC_outlier_vs_control_stat_summary.iter2000.tsv"
)

arg_cc_value_one_summary_file <- file.path(
  arg_diagnostics_dir,
  "ARGweaver_CC_Xin_Bel_CC_value_one_summary.iter2000.tsv"
)

arg_cc_jcr_rcc_scatter_pdf <- file.path(
  arg_diagnostics_dir,
  "ARGweaver_CC_Tap_Xin_JCR_v2_vs_Xin_Bel_CC_original.control_sample.pdf"
)

fd_rth_correlation_pdf <- file.path(
  arg_diagnostics_dir,
  "Tapajos_Xingu_fd_vs_Tap_Xin_RTH_original.10kb_window_correlation.pdf"
)

fd_rth_correlation_data_file <- file.path(
  arg_diagnostics_dir,
  "Tapajos_Xingu_fd_vs_Tap_Xin_RTH_original.10kb_window_matched_values.tsv"
)

fd_rth_correlation_summary_file <- file.path(
  arg_diagnostics_dir,
  "Tapajos_Xingu_fd_vs_Tap_Xin_RTH_original.10kb_window_correlation_summary.tsv"
)

arg_significant_window_counts_file <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_significant_ARG_stat_window_counts.iter2000.tsv"
)

arg_significant_window_barplot_pdf <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_significant_ARG_stat_window_counts.Fst_outlier_windows.iter2000.pdf"
)

arg_significant_window_counts_min5_file <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_significant_ARG_stat_window_counts.min5percent.iter2000.tsv"
)

arg_significant_window_barplot_min5_pdf <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_significant_ARG_stat_window_counts.min5percent.Fst_outlier_windows.iter2000.pdf"
)

arg_min_prop_values <- c(0.05, 0.10, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00)

arg_min_prop_suffixes <- paste0(
  "min",
  as.integer(arg_min_prop_values * 100),
  "percent"
)

arg_significant_window_counts_minprop_files <- setNames(
  file.path(
    arg_barplots_significance_dir,
    paste0(
      "ARGweaver_significant_ARG_stat_window_counts.",
      arg_min_prop_suffixes,
      ".iter2000.tsv"
    )
  ),
  arg_min_prop_suffixes
)

arg_significant_window_barplot_minprop_pdfs <- setNames(
  file.path(
    arg_barplots_significance_dir,
    paste0(
      "ARGweaver_significant_ARG_stat_window_counts.",
      arg_min_prop_suffixes,
      ".Fst_outlier_windows.iter2000.pdf"
    )
  ),
  arg_min_prop_suffixes
)

arg_significant_peak_counts_file <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_significant_ARG_stat_peak_counts.iter2000.tsv"
)

arg_significant_peak_barplot_pdf <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_significant_ARG_stat_peak_counts.Xin_Bel_Fst_peaks.iter2000.pdf"
)

arg_significant_peak_counts_minprop_files <- setNames(
  file.path(
    arg_barplots_significance_dir,
    paste0(
      "ARGweaver_significant_ARG_stat_peak_counts.",
      arg_min_prop_suffixes,
      ".iter2000.tsv"
    )
  ),
  arg_min_prop_suffixes
)

arg_significant_peak_barplot_minprop_pdfs <- setNames(
  file.path(
    arg_barplots_significance_dir,
    paste0(
      "ARGweaver_significant_ARG_stat_peak_counts.",
      arg_min_prop_suffixes,
      ".Xin_Bel_Fst_peaks.iter2000.pdf"
    )
  ),
  arg_min_prop_suffixes
)

arg_window_significance_status_file <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_ARG_stat_significance_status_by_Fst_outlier_window.iter2000.tsv"
)

arg_peak_significance_status_file <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_ARG_stat_significance_status_by_Xin_Bel_Fst_peak.iter2000.tsv"
)

arg_window_significance_status_summary_file <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_ARG_stat_significance_status_by_Fst_outlier_window.summary.iter2000.tsv"
)

arg_peak_significance_status_summary_file <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_ARG_stat_significance_status_by_Xin_Bel_Fst_peak.summary.iter2000.tsv"
)

arg_window_venn_counts_file <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_ARG_stat_venn_counts_by_Fst_outlier_window.iter2000.tsv"
)

arg_peak_venn_counts_file <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_ARG_stat_venn_counts_by_Xin_Bel_Fst_peak.iter2000.tsv"
)

arg_window_ARG_coverage_file <- file.path(
  arg_coverage_dir,
  "ARGweaver_ARG_based_data_coverage_by_Fst_outlier_window.iter2000.tsv"
)

arg_peak_ARG_coverage_file <- file.path(
  arg_coverage_dir,
  "ARGweaver_ARG_based_data_coverage_by_Xin_Bel_Fst_peak.iter2000.tsv"
)

arg_ARG_coverage_summary_file <- file.path(
  arg_coverage_dir,
  "ARGweaver_ARG_based_data_coverage_summary.iter2000.tsv"
)

tap_xin_rth_value_one_summary_file <- file.path(
  arg_diagnostics_dir,
  "ARGweaver_Tap_Xin_RTH_value_one_summary.iter2000.tsv"
)

tap_xin_rth_decile_distribution_file <- file.path(
  arg_diagnostics_dir,
  "ARGweaver_Tap_Xin_RTH_decile_distribution.iter2000.tsv"
)

tap_xin_rth_decile_barplot_pdf <- file.path(
  arg_diagnostics_dir,
  "ARGweaver_Tap_Xin_RTH_decile_distribution.iter2000.pdf"
)

# -----------------------------
# Additional genome-wide statistics
# -----------------------------
d_stats_file <- file.path(
  "..",
  "D_statistics_windows",
  "xiph.Dstats.csv"
)

relernn_file <- file.path(
  "..",
  "ReLERNN",
  "ReLERNN_clean_data.bed"
)

xingu_raisd_file <- file.path(
  "..",
  "RAiSD",
  "combined_RAiSD_output_xingu.sorted.cleaned.bed"
)

belem_raisd_file <- file.path(
  "..",
  "RAiSD",
  "combined_RAiSD_output_belem.sorted.cleaned.bed"
)

xingu_ihs_file <- file.path(
  "..",
  "selscan",
  "Xin_selscan_results",
  "merged_files",
  "Xin_all_ihs.norm.tsv"
)

xingu_nsl_file <- file.path(
  "..",
  "selscan",
  "Xin_selscan_results",
  "merged_files",
  "Xin_all_nsl.norm.tsv"
)

belem_ihs_file <- file.path(
  "..",
  "selscan",
  "Bel_selscan_results",
  "merged_files",
  "Bel_all_ihs.norm.tsv"
)

belem_nsl_file <- file.path(
  "..",
  "selscan",
  "Bel_selscan_results",
  "merged_files",
  "Bel_all_nsl.norm.tsv"
)

rcnv_file <- file.path(
  "..",
  "rCNV",
  "genomewide.rCNV.allele_info_WGS.10kb.tsv"
)

arg_based_fst_file <- file.path(
  "..",
  "ARGweaver",
  "ARG_based_Fst",
  "xipho_arg_based_branch_fst.iter2000.cleaned.tsv"
)

argweaver_mask_file <- file.path(
  "..",
  "ARGweaver",
  "bed_files",
  "ARGweaver_mask.bed"
)

repeat_regions_file <- file.path(
  "..",
  "ARGweaver",
  "bed_files",
  "repeat_regions.bed"
)


bel_popA_allele_stats_file <- file.path(
  "..",
  "allele_stats",
  "bel_popA",
  "allele_stats_by_window_belpopA.csv"
)

xin_popA_allele_stats_file <- file.path(
  "..",
  "allele_stats",
  "xin_popA",
  "allele_stats_by_window_xinpopA.csv"
)

tree_folder_coordinate_map_file <- file.path(
  work_dir,
  "tree_folder_coordinate_map.tsv"
)

manhattan_footnotes_pdf <- file.path(
  work_dir,
  "manhattan_footnotes.pdf"
)

arg_tree_files_midpoint_dir <- file.path(
  "..",
  "ARGweaver",
  "argTreeFiles_midpoint"
)

tree_sampling_seed <- 12345L

additional_threshold_summary_file <- file.path(
  significance_thresholds_dir,
  "additional_control_window_empirical_significance_thresholds.tsv"
)

additional_control_value_summary_file <- file.path(
  significance_thresholds_dir,
  "additional_control_window_value_summary.tsv"
)

heatmap_window_values_file <- file.path(
  heatmap_dir,
  "Xin_Bel_Fst_peak_heatmap_10kb_window_values.tsv"
)

heatmap_window_percentiles_file <- file.path(
  heatmap_dir,
  "Xin_Bel_Fst_peak_heatmap_10kb_window_percentiles.tsv"
)

heatmap_peak_representative_file <- file.path(
  heatmap_dir,
  "Xin_Bel_Fst_peak_heatmap_representative_percentiles.tsv"
)

heatmap_pdf <- file.path(
  heatmap_dir,
  "Xin_Bel_Fst_peak_statistic_percentile_heatmap.pdf"
)

heatmap_pca_any_pdf <- file.path(
  heatmap_dir,
  "Xin_Bel_Fst_peak_statistic_percentile_PCA.at_least_one_tree_signif_in_peak.pdf"
)

heatmap_pca_min25_pdf <- file.path(
  heatmap_dir,
  "Xin_Bel_Fst_peak_statistic_percentile_PCA.25perc_trees_signif_in_peak.pdf"
)

heatmap_pca_any_scores_file <- file.path(
  heatmap_dir,
  "Xin_Bel_Fst_peak_statistic_percentile_PCA_scores.at_least_one_tree_signif_in_peak.tsv"
)

heatmap_pca_min25_scores_file <- file.path(
  heatmap_dir,
  "Xin_Bel_Fst_peak_statistic_percentile_PCA_scores.25perc_trees_signif_in_peak.tsv"
)

figure_3_1_pdf <- file.path(
  figures_dir,
  "Figure_3_1_heatmap.at_least_one_tree_signif_in_peak.pdf"
)

figure_3_2_pdf <- file.path(
  figures_dir,
  "Figure_3_2_heatmap.25perc_trees_signif_in_peak.pdf"
)

figure_4_pdf <- file.path(
  figures_dir,
  "Figure_4_peak7_multi_iter_manhattan_tree_panel.pdf"
)

figure_5_pdf <- file.path(
  figures_dir,
  "Figure_5_peak51_2000th_iter_manhattan_tree_panel.pdf"
)

figure_peak_size_distribution_pdf <- file.path(
  figures_dir,
  "Xin_Bel_Fst_peak_size_distribution.all_Fst_peaks.pdf"
)

figure_peak_size_distribution_any_model_pdf <- file.path(
  figures_dir,
  "Xin_Bel_Fst_peak_size_distribution.model_colored.at_least_one_tree_signif_in_peak.pdf"
)

figure_peak_size_distribution_min25_model_pdf <- file.path(
  figures_dir,
  "Xin_Bel_Fst_peak_size_distribution.model_colored.25perc_trees_signif_in_peak.pdf"
)

figure_ARG_barplot_p_sensitivity_any_pdf <- file.path(
  figures_dir,
  "ARGweaver_significant_ARG_stat_barplots.P_sensitivity.any_significant_local_tree.pdf"
)

figure_ARG_barplot_p_sensitivity_min25_pdf <- file.path(
  figures_dir,
  "ARGweaver_significant_ARG_stat_barplots.P_sensitivity.min25percent_significant_local_trees.pdf"
)

figure_large_repeat_adjacent_Fst_manhattan_pdf <- file.path(
  figures_dir,
  "Xin_Bel_Fst_peaks_near_large_repeat_regions.with_ARG_data.Manhattan_panels.pdf"
)

large_repeat_peak_proximity_file <- file.path(
  fst_outlier_dir,
  "Xin_Bel_Fst_peaks_with_large_repeat_regions_within_50kb.tsv"
)

model_assignment_any_dir <- file.path(
  output_dir,
  "Model_assignments_at_least_one_tree_signif_in_peak"
)

model_assignment_min25_dir <- file.path(
  output_dir,
  "Model_assignments_25perc_trees_signif_in_peak"
)

dir.create(model_assignment_any_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(model_assignment_min25_dir, showWarnings = FALSE, recursive = TRUE)

model_assignment_any_file <- file.path(
  model_assignment_any_dir,
  "Xin_Bel_Fst_peak_model_assignments.at_least_one_tree_signif_in_peak.tsv"
)

model_assignment_any_summary_file <- file.path(
  model_assignment_any_dir,
  "Xin_Bel_Fst_peak_model_assignment_summary.at_least_one_tree_signif_in_peak.tsv"
)

model_assignment_min25_file <- file.path(
  model_assignment_min25_dir,
  "Xin_Bel_Fst_peak_model_assignments.25perc_trees_signif_in_peak.tsv"
)

model_assignment_min25_summary_file <- file.path(
  model_assignment_min25_dir,
  "Xin_Bel_Fst_peak_model_assignment_summary.25perc_trees_signif_in_peak.tsv"
)

manhattan_any_dir <- file.path(
  model_assignment_any_dir,
  "Manhattan_panels_2000th_iter"
)

manhattan_min25_dir <- file.path(
  model_assignment_min25_dir,
  "Manhattan_panels_2000th_iter"
)

manhattan_any_multi_iter_dir <- file.path(
  model_assignment_any_dir,
  "Manhattan_panels_multi_iter"
)

manhattan_min25_multi_iter_dir <- file.path(
  model_assignment_min25_dir,
  "Manhattan_panels_multi_iter"
)

dir.create(manhattan_any_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(manhattan_min25_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(manhattan_any_multi_iter_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(manhattan_min25_multi_iter_dir, showWarnings = FALSE, recursive = TRUE)

manhattan_any_summary_file <- file.path(
  manhattan_any_dir,
  "Manhattan_panel_summary.at_least_one_tree_signif_in_peak.tsv"
)

manhattan_min25_summary_file <- file.path(
  manhattan_min25_dir,
  "Manhattan_panel_summary.25perc_trees_signif_in_peak.tsv"
)

manhattan_any_multi_iter_summary_file <- file.path(
  manhattan_any_multi_iter_dir,
  "Manhattan_panel_summary.multi_iter.at_least_one_tree_signif_in_peak.tsv"
)

manhattan_min25_multi_iter_summary_file <- file.path(
  manhattan_min25_multi_iter_dir,
  "Manhattan_panel_summary.multi_iter.25perc_trees_signif_in_peak.tsv"
)

arg_multi_iter_window_counts_file <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_significant_ARG_stat_window_counts.multi_iter.tsv"
)

arg_multi_iter_peak_counts_file <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_significant_ARG_stat_peak_counts.multi_iter.tsv"
)

arg_multi_iter_model_assignments_file <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_peak_model_assignments.multi_iter.tsv"
)

arg_multi_iter_model_assignment_summary_file <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_peak_model_assignment_summary.multi_iter.tsv"
)

arg_50iter_stat_summary_any_file <- file.path(
  arg_barplots_significance_dir,
  "Summary_ARG_statistics_50_MCMC_iterations_at_least_one_tree_signif_in_peak.tsv"
)

arg_50iter_stat_summary_min25_file <- file.path(
  arg_barplots_significance_dir,
  "Summary_ARG_statistics_50_MCMC_iterations_at_25perc_trees_signif_in_peak.tsv"
)

model_category_definitions_xlsx <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_peak_model_category_definitions.xlsx"
)

model_category_definitions_tsv <- file.path(
  arg_barplots_significance_dir,
  "ARGweaver_peak_model_category_definitions.tsv"
)


fst_peak_control_subsampling_any_file <- file.path(
  tables_dir,
  "Fst_peak_significance_overrepresentation_control_subsampling.at_least_one_tree_signif_in_peak.tsv"
)

fst_peak_control_subsampling_min25_file <- file.path(
  tables_dir,
  "Fst_peak_significance_overrepresentation_control_subsampling.25perc_trees_signif_in_peak.tsv"
)

# -----------------------------
# Helper functions
# -----------------------------

# Summarize one ARG statistic in every Xin-Bel Fst control window. This uses the
# same local-tree significance rule used for the empirical Fst peaks.
make_ARG_control_window_status <- function(
  arg_joined_dt,
  control_windows_dt,
  stat,
  threshold,
  tail,
  min_prop = NA_real_
) {
  if (!stat %in% names(arg_joined_dt)) {
    stop("Statistic is absent from arg_joined_dt: ", stat)
  }
  if (length(threshold) == 0 || is.na(threshold)) {
    stop("Missing significance threshold for control subsampling statistic: ", stat)
  }

  control_ref <- control_windows_dt[, .(
    chromosome,
    start = as.integer(start),
    end = as.integer(end),
    window_id = paste(chromosome, start, end, sep = ":")
  )]

  tmp <- arg_joined_dt[
    window_class == "control" &
      !is.na(get(stat)),
    .(
      n_local_trees = .N,
      n_significant_local_trees = if (tail == "upper") {
        sum(get(stat) >= threshold)
      } else if (tail == "lower") {
        sum(get(stat) <= threshold)
      } else {
        stop("Unexpected tail value: ", tail)
      }
    ),
    by = .(chromosome, window_start, window_end, window_id)
  ]

  if (nrow(tmp) > 0) {
    setnames(
      tmp,
      old = c("window_start", "window_end"),
      new = c("start", "end")
    )
  }

  out <- merge(
    control_ref,
    tmp,
    by = c("chromosome", "start", "end", "window_id"),
    all.x = TRUE,
    sort = FALSE
  )

  out[is.na(n_local_trees), n_local_trees := 0L]
  out[is.na(n_significant_local_trees), n_significant_local_trees := 0L]

  out[, informative := n_local_trees > 0L]
  out[, prop_significant_local_trees := fifelse(
    informative,
    n_significant_local_trees / n_local_trees,
    NA_real_
  )]

  if (is.na(min_prop)) {
    out[, is_significant := informative & n_significant_local_trees >= 1L]
  } else {
    out[, is_significant :=
          informative &
          n_significant_local_trees >= 1L &
          !is.na(prop_significant_local_trees) &
          prop_significant_local_trees >= min_prop]
  }

  setorder(out, chromosome, start, end)
  out[]
}



# Build significance/informativeness flags for statistics that are already
# represented by one value per pixy control window.
make_aligned_control_window_status <- function(
  control_windows_dt,
  stat,
  threshold,
  tail,
  use_abs = FALSE
) {
  out <- control_windows_dt[, .(
    chromosome,
    start = as.integer(start),
    end = as.integer(end),
    window_id = paste(chromosome, start, end, sep = ":"),
    value = as.numeric(get(stat))
  )]

  out[, informative := !is.na(value)]

  if (use_abs) {
    out[, is_significant := informative & abs(value) > threshold]
  } else if (tail == "upper") {
    out[, is_significant := informative & value > threshold]
  } else if (tail == "lower") {
    out[, is_significant := informative & value < threshold]
  } else {
    stop("Unexpected tail value: ", tail)
  }

  out[, `:=`(
    n_values = as.integer(informative),
    n_significant_values = as.integer(is_significant)
  )]

  out[, .(
    chromosome, start, end, window_id,
    n_values, n_significant_values,
    informative, is_significant
  )]
}


# Build control-window significance flags from raw interval records. A control
# window is informative if at least one retained raw value overlaps it, and is
# significant if at least one overlapping value passes the same rule used for
# Fst-peak significance.
make_raw_interval_control_window_status <- function(
  dt,
  stat,
  control_windows_dt,
  threshold = NA_real_,
  tail = "upper",
  use_abs = FALSE,
  nonzero_is_significant = FALSE
) {
  control_ref <- control_windows_dt[, .(
    chromosome = standardize_scaffold_names(chromosome),
    start = as.integer(start),
    end = as.integer(end),
    window_end_inclusive = as.integer(end - 1L),
    window_id = paste(chromosome, start, end, sep = ":")
  )]

  tmp <- copy(dt)
  tmp[, chromosome := standardize_scaffold_names(chromosome)]
  tmp <- tmp[!is.na(get(stat))]
  tmp[, interval_start := as.integer(interval_start)]
  tmp[, interval_end_inclusive := as.integer(interval_end - 1L)]
  tmp <- tmp[interval_end_inclusive >= interval_start]

  if (nrow(tmp) == 0) {
    return(control_ref[, .(
      chromosome, start, end, window_id,
      n_values = 0L,
      n_significant_values = 0L,
      informative = FALSE,
      is_significant = FALSE
    )])
  }

  if (nonzero_is_significant) {
    tmp[, raw_is_significant := get(stat) != 0]
  } else if (use_abs) {
    tmp[, raw_is_significant := abs(get(stat)) > threshold]
  } else if (tail == "upper") {
    tmp[, raw_is_significant := get(stat) > threshold]
  } else if (tail == "lower") {
    tmp[, raw_is_significant := get(stat) < threshold]
  } else {
    stop("Unexpected tail value: ", tail)
  }

  setkey(tmp, chromosome, interval_start, interval_end_inclusive)
  setkey(control_ref, chromosome, start, window_end_inclusive)

  hits <- foverlaps(
    tmp[, .(
      chromosome,
      interval_start,
      interval_end_inclusive,
      raw_is_significant
    )],
    control_ref,
    by.x = c("chromosome", "interval_start", "interval_end_inclusive"),
    by.y = c("chromosome", "start", "window_end_inclusive"),
    nomatch = 0
  )

  status <- hits[, .(
    n_values = .N,
    n_significant_values = sum(raw_is_significant, na.rm = TRUE)
  ), by = window_id]

  out <- merge(
    control_ref[, .(chromosome, start, end, window_id)],
    status,
    by = "window_id",
    all.x = TRUE,
    sort = FALSE
  )

  out[is.na(n_values), n_values := 0L]
  out[is.na(n_significant_values), n_significant_values := 0L]
  out[, informative := n_values > 0L]
  out[, is_significant := informative & n_significant_values > 0L]
  out[]
}


make_raw_point_control_window_status <- function(
  dt,
  stat,
  control_windows_dt,
  threshold = NA_real_,
  tail = "upper",
  use_abs = FALSE,
  nonzero_is_significant = FALSE
) {
  tmp <- copy(dt)
  tmp[, interval_start := as.integer(pos)]
  tmp[, interval_end := as.integer(pos + 1L)]
  make_raw_interval_control_window_status(
    dt = tmp,
    stat = stat,
    control_windows_dt = control_windows_dt,
    threshold = threshold,
    tail = tail,
    use_abs = use_abs,
    nonzero_is_significant = nonzero_is_significant
  )
}


# ARG-based branch-Fst values are raw point records, but the >=25% assignment
# rule is based on the proportion of significant point/local-tree values within
# each control window. This helper mirrors that rule.
make_point_proportion_control_window_status <- function(
  dt,
  stat,
  control_windows_dt,
  threshold,
  tail,
  min_prop = NA_real_,
  inclusive_threshold = FALSE
) {
  control_ref <- control_windows_dt[, .(
    chromosome = standardize_scaffold_names(chromosome),
    start = as.integer(start),
    end = as.integer(end),
    window_end_inclusive = as.integer(end - 1L),
    window_id = paste(chromosome, start, end, sep = ":")
  )]

  tmp <- copy(dt)
  tmp[, chromosome := standardize_scaffold_names(chromosome)]
  tmp <- tmp[!is.na(get(stat))]
  tmp[, point_start := as.integer(pos)]
  tmp[, point_end := as.integer(pos)]

  if (tail == "upper") {
    if (inclusive_threshold) {
      tmp[, point_is_significant := get(stat) >= threshold]
    } else {
      tmp[, point_is_significant := get(stat) > threshold]
    }
  } else if (tail == "lower") {
    if (inclusive_threshold) {
      tmp[, point_is_significant := get(stat) <= threshold]
    } else {
      tmp[, point_is_significant := get(stat) < threshold]
    }
  } else {
    stop("Unexpected tail value: ", tail)
  }

  if (nrow(tmp) == 0) {
    return(control_ref[, .(
      chromosome, start, end, window_id,
      n_local_trees = 0L,
      n_significant_local_trees = 0L,
      prop_significant_local_trees = NA_real_,
      informative = FALSE,
      is_significant = FALSE
    )])
  }

  setkey(tmp, chromosome, point_start, point_end)
  setkey(control_ref, chromosome, start, window_end_inclusive)

  hits <- foverlaps(
    tmp[, .(chromosome, point_start, point_end, point_is_significant)],
    control_ref,
    by.x = c("chromosome", "point_start", "point_end"),
    by.y = c("chromosome", "start", "window_end_inclusive"),
    nomatch = 0
  )

  status <- hits[, .(
    n_local_trees = .N,
    n_significant_local_trees = sum(point_is_significant, na.rm = TRUE)
  ), by = window_id]

  out <- merge(
    control_ref[, .(chromosome, start, end, window_id)],
    status,
    by = "window_id",
    all.x = TRUE,
    sort = FALSE
  )

  out[is.na(n_local_trees), n_local_trees := 0L]
  out[is.na(n_significant_local_trees), n_significant_local_trees := 0L]
  out[, informative := n_local_trees > 0L]
  out[, prop_significant_local_trees := fifelse(
    informative,
    n_significant_local_trees / n_local_trees,
    NA_real_
  )]

  if (is.na(min_prop)) {
    out[, is_significant := informative & n_significant_local_trees >= 1L]
  } else {
    out[, is_significant :=
          informative &
          n_significant_local_trees >= 1L &
          !is.na(prop_significant_local_trees) &
          prop_significant_local_trees >= min_prop]
  }

  out[]
}


# Append one statistic column to a wide subsampling output table while retaining
# the single Descriptor column requested for the final TSV files.
append_control_subsampling_column <- function(
  output_dt,
  statistic_column_name,
  statistic_display_name,
  empirical_count,
  subsample_counts
) {
  col_dt <- make_control_subsampling_output_column(
    statistic_name = statistic_display_name,
    empirical_count = empirical_count,
    subsample_counts = subsample_counts
  )
  setnames(col_dt, "value", statistic_column_name)

  if (is.null(output_dt)) {
    return(col_dt)
  }

  merge(
    output_dt,
    col_dt,
    by = "Descriptor",
    all = TRUE,
    sort = FALSE
  )
}


# Enumerate every contiguous region made exclusively from Fst control windows
# that exactly matches a requested empirical peak length.
#
# A candidate region is retained only if it contains at least one informative
# control window for the statistic being evaluated. Individual windows within
# the region are allowed to lack values; this mirrors the peak analysis, where
# a peak can remain informative when ARG data exist in only part of the peak.
enumerate_exact_length_control_regions <- function(
  control_status_dt,
  region_length_bp
) {
  if (nrow(control_status_dt) == 0) {
    return(data.table())
  }

  dt <- copy(control_status_dt)
  setorder(dt, chromosome, start, end)

  # A new run begins whenever adjacent control windows are not directly
  # contiguous. Therefore every candidate below consists only of uninterrupted
  # control windows; intermediate/outlier windows cannot occur inside it.
  dt[, previous_end := shift(end), by = chromosome]
  dt[, new_control_run :=
       is.na(previous_end) | start != previous_end,
     by = chromosome]
  dt[, control_run_id := cumsum(new_control_run)]
  dt[, previous_end := NULL]

  dt[, run_index := seq_len(.N), by = .(chromosome, control_run_id)]
  dt[, cum_informative := cumsum(as.integer(informative)),
     by = .(chromosome, control_run_id)]
  dt[, cum_significant := cumsum(as.integer(is_significant)),
     by = .(chromosome, control_run_id)]

  starts <- dt[, .(
    start_index = run_index,
    region_start = start,
    start_cum_informative_before =
      shift(cum_informative, fill = 0L),
    start_cum_significant_before =
      shift(cum_significant, fill = 0L)
  ), by = .(chromosome, control_run_id)]

  starts[, target_end := as.integer(region_start + region_length_bp)]

  ends <- dt[, .(
    chromosome,
    control_run_id,
    end_index = run_index,
    target_end = end,
    end_cum_informative = cum_informative,
    end_cum_significant = cum_significant
  )]

  candidates <- merge(
    starts,
    ends,
    by = c("chromosome", "control_run_id", "target_end"),
    all = FALSE,
    sort = FALSE
  )

  candidates <- candidates[end_index >= start_index]
  if (nrow(candidates) == 0) {
    return(data.table())
  }

  candidates[, n_informative_windows :=
               end_cum_informative - start_cum_informative_before]
  candidates[, n_significant_windows :=
               end_cum_significant - start_cum_significant_before]

  candidates <- candidates[n_informative_windows > 0L]

  candidates[, `:=`(
    region_end = target_end,
    region_length_bp = as.integer(region_length_bp),
    region_is_significant = n_significant_windows > 0L
  )]

  candidates[, .(
    chromosome,
    region_start,
    region_end,
    region_length_bp,
    n_informative_windows,
    n_significant_windows,
    region_is_significant
  )]
}


# Generate the null distribution for one statistic/assignment method.
#
# Each null replicate contains exactly one sampled control region for every
# empirical Fst peak, preserving the complete empirical peak-length
# distribution. Sampling is with replacement, as requested.
subsample_control_regions_matching_Fst_peak_sizes <- function(
  peaks_dt,
  control_status_dt,
  n_subsamples = 10000L,
  seed = 1L
) {
  peak_lengths <- as.integer(peaks_dt$end - peaks_dt$start)

  # The caller supplies the empirical peak set to be matched. In the current
  # overrepresentation analysis this is the 141 Fst peaks with ARG-based data.
  # Every null replicate therefore inherits exactly the number and length
  # distribution of the supplied empirical peaks.
  if (length(peak_lengths) == 0L) {
    stop("No empirical Fst peaks were supplied for matched-control subsampling.")
  }

  unique_lengths <- sort(unique(peak_lengths))

  candidate_by_length <- setNames(
    lapply(unique_lengths, function(len) {
      enumerate_exact_length_control_regions(
        control_status_dt = control_status_dt,
        region_length_bp = len
      )
    }),
    as.character(unique_lengths)
  )

  n_candidates_by_length <- vapply(
    candidate_by_length,
    nrow,
    integer(1)
  )

  missing_lengths <- unique_lengths[n_candidates_by_length == 0L]
  if (length(missing_lengths) > 0L) {
    stop(
      "No informative contiguous control region could be found for the ",
      "following empirical Fst peak length(s): ",
      paste(missing_lengths, collapse = ", "),
      " bp. Exact size matching therefore cannot be completed."
    )
  }

  set.seed(seed)

  subsample_counts <- integer(n_subsamples)

  for (subsample_i in seq_len(n_subsamples)) {
    sampled_significance <- logical(length(peak_lengths))

    for (peak_i in seq_along(peak_lengths)) {
      candidate_dt <- candidate_by_length[[as.character(peak_lengths[peak_i])]]
      sampled_row <- sample.int(nrow(candidate_dt), size = 1L)
      sampled_significance[peak_i] <-
        candidate_dt$region_is_significant[sampled_row]
    }

    # Every replicate contains the same number and exact size distribution of
    # regions as the empirical Fst-peak dataset.
    if (length(sampled_significance) != nrow(peaks_dt)) {
      stop("Internal error: sampled control-region count does not match Fst peak count.")
    }

    subsample_counts[subsample_i] <- sum(sampled_significance)
  }

  list(
    subsample_counts = subsample_counts,
    candidate_counts_by_length = data.table(
      region_length_bp = unique_lengths,
      n_candidate_control_regions = n_candidates_by_length
    )
  )
}


# Convert one null distribution into the requested column format.
# The P thresholds are upper-tail critical values because the focal hypothesis
# is overrepresentation of significant regions inside Fst peaks.
make_control_subsampling_output_column <- function(
  statistic_name,
  empirical_count,
  subsample_counts
) {
  critical_values <- as.integer(quantile(
    subsample_counts,
    probs = c(0.95, 0.99, 0.999),
    type = 1,
    names = FALSE
  ))

  data.table(
    Descriptor = c(
      "Statistic_name",
      "Empirical_count",
      "P_0.05_threshold",
      "P_0.01_threshold",
      "P_0.001_threshold",
      paste0("Subsample ", seq_along(subsample_counts))
    ),
    value = c(
      statistic_name,
      as.character(empirical_count),
      as.character(critical_values[1]),
      as.character(critical_values[2]),
      as.character(critical_values[3]),
      as.character(subsample_counts)
    )
  )
}


safe_inverse <- function(x) {
  y <- rep(NA_real_, length(x))
  y[!is.na(x) & x != 0] <- 1 / x[!is.na(x) & x != 0]
  y
}

summarize_stats <- function(dt, class_name, stats) {
  rbindlist(lapply(stats, function(stat) {
    values <- dt[[stat]]
    non_na_values <- values[!is.na(values)]

    data.table(
      window_class = class_name,
      statistic = stat,
      n_windows = length(non_na_values),
      mean = mean(non_na_values),
      sd = sd(non_na_values),
      min = min(non_na_values),
      max = max(non_na_values)
    )
  }))
}

summarize_pi_recombination_by_region_set <- function(
  all_windows_dt,
  outlier_windows_dt,
  control_windows_dt,
  recombination_dt
) {
  region_sets <- list(
    all_regions = copy(all_windows_dt),
    Xin_Bel_Fst_outlier = copy(outlier_windows_dt),
    control = copy(control_windows_dt)
  )

  rbindlist(lapply(names(region_sets), function(region_set_name) {
    win <- region_sets[[region_set_name]]

    if (nrow(win) == 0) {
      return(data.table(
        region_set = region_set_name,
        n_windows = 0L,
        mean_Tap_pi = NA_real_,
        mean_Xin_pi = NA_real_,
        mean_Bel_pi = NA_real_,
        mean_recombination_rate = NA_real_
      ))
    }

    win[, window_id := paste(chromosome, start, end, sep = ":")]
    win_intervals <- win[, .(
      window_id,
      chromosome,
      window_start = as.integer(start),
      window_end_inclusive = as.integer(end - 1L)
    )]

    recomb_by_window <- mean_by_window_from_intervals(
      recombination_dt,
      "recombination_rate",
      win_intervals
    )

    data.table(
      region_set = region_set_name,
      n_windows = nrow(win),
      mean_Tap_pi = mean(win$Tap_pi, na.rm = TRUE),
      mean_Xin_pi = mean(win$Xin_pi, na.rm = TRUE),
      mean_Bel_pi = mean(win$Bel_pi, na.rm = TRUE),
      mean_recombination_rate = mean(recomb_by_window$value, na.rm = TRUE)
    )
  }), fill = TRUE)
}


make_outlier_window_selection_flags_from_intervals <- function(
  dt,
  stat,
  outlier_windows_dt,
  threshold,
  use_abs = FALSE
) {
  if (nrow(dt) == 0 || nrow(outlier_windows_dt) == 0 || is.na(threshold)) {
    return(outlier_windows_dt[, .(
      window_id = paste(chromosome, start, end, sep = ":"),
      statistic = stat,
      n_significant_values_in_Fst_outlier_window = 0L,
      significant_selection_stat_in_Fst_outlier_window = FALSE
    )])
  }

  win <- outlier_windows_dt[, .(
    window_id = paste(chromosome, start, end, sep = ":"),
    chromosome,
    window_start = as.integer(start),
    window_end_inclusive = as.integer(end - 1L)
  )]
  win[, chromosome := standardize_scaffold_names(chromosome)]

  tmp <- copy(dt)
  tmp[, chromosome := standardize_scaffold_names(chromosome)]
  tmp <- tmp[!is.na(get(stat))]
  if (nrow(tmp) == 0) {
    return(win[, .(
      window_id,
      statistic = stat,
      n_significant_values_in_Fst_outlier_window = 0L,
      significant_selection_stat_in_Fst_outlier_window = FALSE
    )])
  }

  tmp[, interval_start := as.integer(interval_start)]
  tmp[, interval_end_inclusive := as.integer(interval_end - 1L)]
  tmp <- tmp[interval_end_inclusive >= interval_start]
  if (use_abs) {
    tmp <- tmp[abs(get(stat)) > threshold]
  } else {
    tmp <- tmp[get(stat) > threshold]
  }

  if (nrow(tmp) == 0) {
    return(win[, .(
      window_id,
      statistic = stat,
      n_significant_values_in_Fst_outlier_window = 0L,
      significant_selection_stat_in_Fst_outlier_window = FALSE
    )])
  }

  setkey(tmp, chromosome, interval_start, interval_end_inclusive)
  setkey(win, chromosome, window_start, window_end_inclusive)

  joined <- foverlaps(
    tmp[, .(chromosome, interval_start, interval_end_inclusive)],
    win,
    by.x = c("chromosome", "interval_start", "interval_end_inclusive"),
    by.y = c("chromosome", "window_start", "window_end_inclusive"),
    nomatch = 0
  )

  counts <- joined[, .(
    n_significant_values_in_Fst_outlier_window = .N
  ), by = window_id]

  out <- merge(
    win[, .(window_id)],
    counts,
    by = "window_id",
    all.x = TRUE,
    sort = FALSE
  )
  out[is.na(n_significant_values_in_Fst_outlier_window),
      n_significant_values_in_Fst_outlier_window := 0L]
  out[, statistic := stat]
  out[, significant_selection_stat_in_Fst_outlier_window :=
        n_significant_values_in_Fst_outlier_window > 0]
  out[]
}

make_outlier_window_selection_flags_from_points <- function(
  dt,
  stat,
  outlier_windows_dt,
  threshold,
  use_abs = FALSE
) {
  if (nrow(dt) == 0 || nrow(outlier_windows_dt) == 0 || is.na(threshold)) {
    return(outlier_windows_dt[, .(
      window_id = paste(chromosome, start, end, sep = ":"),
      statistic = stat,
      n_significant_values_in_Fst_outlier_window = 0L,
      significant_selection_stat_in_Fst_outlier_window = FALSE
    )])
  }

  win <- outlier_windows_dt[, .(
    window_id = paste(chromosome, start, end, sep = ":"),
    chromosome,
    window_start = as.integer(start),
    window_end_inclusive = as.integer(end - 1L)
  )]
  win[, chromosome := standardize_scaffold_names(chromosome)]

  tmp <- copy(dt)
  tmp[, chromosome := standardize_scaffold_names(chromosome)]
  tmp <- tmp[!is.na(get(stat))]
  if (nrow(tmp) == 0) {
    return(win[, .(
      window_id,
      statistic = stat,
      n_significant_values_in_Fst_outlier_window = 0L,
      significant_selection_stat_in_Fst_outlier_window = FALSE
    )])
  }

  if (use_abs) {
    tmp <- tmp[abs(get(stat)) > threshold]
  } else {
    tmp <- tmp[get(stat) > threshold]
  }

  if (nrow(tmp) == 0) {
    return(win[, .(
      window_id,
      statistic = stat,
      n_significant_values_in_Fst_outlier_window = 0L,
      significant_selection_stat_in_Fst_outlier_window = FALSE
    )])
  }

  tmp[, point_start := as.integer(pos)]
  tmp[, point_end := as.integer(pos)]

  setkey(tmp, chromosome, point_start, point_end)
  setkey(win, chromosome, window_start, window_end_inclusive)

  joined <- foverlaps(
    tmp[, .(chromosome, point_start, point_end)],
    win,
    by.x = c("chromosome", "point_start", "point_end"),
    by.y = c("chromosome", "window_start", "window_end_inclusive"),
    nomatch = 0
  )

  counts <- joined[, .(
    n_significant_values_in_Fst_outlier_window = .N
  ), by = window_id]

  out <- merge(
    win[, .(window_id)],
    counts,
    by = "window_id",
    all.x = TRUE,
    sort = FALSE
  )
  out[is.na(n_significant_values_in_Fst_outlier_window),
      n_significant_values_in_Fst_outlier_window := 0L]
  out[, statistic := stat]
  out[, significant_selection_stat_in_Fst_outlier_window :=
        n_significant_values_in_Fst_outlier_window > 0]
  out[]
}

make_peak_selection_significance_flags <- function(
  outlier_windows_with_peak_dt,
  selection_window_flag_list
) {
  selection_stats <- c(
    "Xingu_RAiSD_u",
    "xingu_norm_ihs",
    "xingu_norm_nsl",
    "Belem_RAiSD_u",
    "belem_norm_ihs",
    "belem_norm_nsl"
  )

  peak_ref <- unique(outlier_windows_with_peak_dt[, .(peak_id)])
  if (nrow(peak_ref) == 0) {
    return(data.table())
  }

  win_peak <- outlier_windows_with_peak_dt[, .(
    peak_id,
    window_id = paste(chromosome, start, end, sep = ":")
  )]

  out <- copy(peak_ref)

  for (stat in selection_stats) {
    flag_dt <- selection_window_flag_list[[stat]]
    flag_col_name <- paste0("significant_", stat, "_in_Fst_outlier_windows")
    count_col_name <- paste0("n_significant_", stat, "_values_in_Fst_outlier_windows")

    if (is.null(flag_dt) || nrow(flag_dt) == 0) {
      out[, (flag_col_name) := FALSE]
      out[, (count_col_name) := 0L]
      next
    }

    # These peak-level flags are based on raw selection-statistic values, not
    # 10-kb averages. A peak is flagged if any raw value exceeding the relevant
    # threshold overlaps any 10-kb Fst-outlier window in that peak. For iHS and
    # nSL, the raw absolute value is evaluated against threshold 2 upstream.
    tmp <- merge(
      win_peak,
      flag_dt[
        statistic == stat,
        .(
          window_id,
          n_significant_values_in_Fst_outlier_window,
          significant_selection_stat_in_Fst_outlier_window
        )
      ],
      by = "window_id",
      all.x = TRUE,
      sort = FALSE
    )

    tmp[is.na(n_significant_values_in_Fst_outlier_window),
        n_significant_values_in_Fst_outlier_window := 0L]
    tmp[is.na(significant_selection_stat_in_Fst_outlier_window),
        significant_selection_stat_in_Fst_outlier_window := FALSE]

    peak_flags <- tmp[, .(
      flag = any(
        significant_selection_stat_in_Fst_outlier_window %in% TRUE,
        na.rm = TRUE
      ),
      n_significant_values = sum(
        n_significant_values_in_Fst_outlier_window,
        na.rm = TRUE
      )
    ), by = peak_id]

    out <- merge(
      out,
      peak_flags,
      by = "peak_id",
      all.x = TRUE,
      sort = FALSE
    )
    setnames(out, "flag", flag_col_name)
    setnames(out, "n_significant_values", count_col_name)
    out[is.na(get(flag_col_name)), (flag_col_name) := FALSE]
    out[is.na(get(count_col_name)), (count_col_name) := 0L]
  }

  out[]
}



# Build a TRUE/FALSE peak-level flag from values already aligned to the pixy
# 10-kb windows. Only the Fst-outlier windows belonging to each peak are tested.
make_peak_flag_from_pixy_windows <- function(
  outlier_windows_dt,
  stat,
  threshold,
  tail,
  output_col
) {
  out <- unique(outlier_windows_dt[, .(peak_id)])
  out[, (output_col) := FALSE]
  if (!stat %in% names(outlier_windows_dt) || is.na(threshold)) return(out[])

  tmp <- outlier_windows_dt[!is.na(get(stat))]
  if (tail == "upper") {
    tmp[, is_sig := get(stat) > threshold]
  } else if (tail == "lower") {
    tmp[, is_sig := get(stat) < threshold]
  } else {
    stop("Unexpected tail: ", tail)
  }

  flags <- tmp[, .(flag = any(is_sig, na.rm = TRUE)), by = peak_id]
  out[flags, on = "peak_id", (output_col) := i.flag]
  out[]
}

# Build a peak-level flag from raw interval records. A peak is TRUE if at least
# one significant raw interval overlaps any Fst-outlier window in that peak.
make_peak_flag_from_raw_intervals <- function(
  dt,
  stat,
  outlier_windows_dt,
  threshold = NA_real_,
  tail = "upper",
  output_col,
  use_abs = FALSE,
  nonzero_is_significant = FALSE
) {
  out <- unique(outlier_windows_dt[, .(peak_id)])
  out[, (output_col) := FALSE]
  if (nrow(dt) == 0 || !stat %in% names(dt)) return(out[])

  tmp <- copy(dt)
  tmp[, chromosome := standardize_scaffold_names(chromosome)]
  tmp <- tmp[!is.na(get(stat))]
  if (nonzero_is_significant) {
    tmp <- tmp[get(stat) != 0]
  } else if (use_abs) {
    tmp <- tmp[abs(get(stat)) > threshold]
  } else if (tail == "upper") {
    tmp <- tmp[get(stat) > threshold]
  } else if (tail == "lower") {
    tmp <- tmp[get(stat) < threshold]
  } else {
    stop("Unexpected tail: ", tail)
  }
  if (nrow(tmp) == 0) return(out[])

  tmp[, interval_start := as.integer(interval_start)]
  tmp[, interval_end_inclusive := as.integer(interval_end - 1L)]
  tmp <- tmp[interval_end_inclusive >= interval_start]

  wins <- outlier_windows_dt[, .(
    peak_id,
    chromosome = standardize_scaffold_names(chromosome),
    window_start = as.integer(start),
    window_end_inclusive = as.integer(end - 1L)
  )]

  setkey(tmp, chromosome, interval_start, interval_end_inclusive)
  setkey(wins, chromosome, window_start, window_end_inclusive)
  hits <- foverlaps(
    tmp[, .(chromosome, interval_start, interval_end_inclusive)],
    wins,
    by.x = c("chromosome", "interval_start", "interval_end_inclusive"),
    by.y = c("chromosome", "window_start", "window_end_inclusive"),
    nomatch = 0
  )
  if (nrow(hits) > 0) {
    flagged <- unique(hits[, .(peak_id, flag = TRUE)])
    out[flagged, on = "peak_id", (output_col) := i.flag]
  }
  out[]
}

make_peak_flag_from_raw_points <- function(
  dt,
  stat,
  outlier_windows_dt,
  threshold,
  tail = "upper",
  output_col,
  use_abs = FALSE
) {
  tmp <- copy(dt)
  tmp[, interval_start := as.integer(pos)]
  tmp[, interval_end := as.integer(pos + 1L)]
  make_peak_flag_from_raw_intervals(
    dt = tmp,
    stat = stat,
    outlier_windows_dt = outlier_windows_dt,
    threshold = threshold,
    tail = tail,
    output_col = output_col,
    use_abs = use_abs
  )
}

# ARG point statistics follow the assignment method represented by each table.
# For the >=25% approach, a peak is TRUE when at least one constituent Fst
# outlier window has >=25% significant local-tree rows for the statistic.
make_ARG_peak_flag_by_assignment_method <- function(
  dt,
  stat,
  outlier_windows_dt,
  threshold,
  tail,
  output_col,
  min_prop = NA_real_
) {
  out <- unique(outlier_windows_dt[, .(peak_id)])
  out[, (output_col) := FALSE]
  if (nrow(dt) == 0 || !stat %in% names(dt) || is.na(threshold)) return(out[])

  pts <- copy(dt)
  pts[, chromosome := standardize_scaffold_names(chromosome)]
  pts <- pts[!is.na(get(stat))]
  pts[, point_start := as.integer(pos)]
  pts[, point_end := as.integer(pos)]
  if (tail == "upper") pts[, is_sig := get(stat) > threshold]
  else if (tail == "lower") pts[, is_sig := get(stat) < threshold]
  else stop("Unexpected tail: ", tail)

  wins <- outlier_windows_dt[, .(
    peak_id,
    window_id = paste(chromosome, start, end, sep = ":"),
    chromosome = standardize_scaffold_names(chromosome),
    window_start = as.integer(start),
    window_end_inclusive = as.integer(end - 1L)
  )]
  setkey(pts, chromosome, point_start, point_end)
  setkey(wins, chromosome, window_start, window_end_inclusive)
  joined <- foverlaps(
    pts[, .(chromosome, point_start, point_end, is_sig)],
    wins,
    by.x = c("chromosome", "point_start", "point_end"),
    by.y = c("chromosome", "window_start", "window_end_inclusive"),
    nomatch = 0
  )
  if (nrow(joined) == 0) return(out[])

  window_status <- joined[, .(
    n_values = .N,
    n_significant = sum(is_sig, na.rm = TRUE),
    prop_significant = mean(is_sig, na.rm = TRUE)
  ), by = .(peak_id, window_id)]
  if (is.na(min_prop)) {
    window_status[, window_is_significant := n_significant > 0L]
  } else {
    window_status[, window_is_significant := prop_significant >= min_prop]
  }
  peak_flags <- window_status[, .(
    flag = any(window_is_significant, na.rm = TRUE)
  ), by = peak_id]
  out[peak_flags, on = "peak_id", (output_col) := i.flag]
  out[]
}

merge_peak_flag_tables <- function(flag_tables, all_peak_ids) {
  out <- data.table(peak_id = all_peak_ids)
  for (x in flag_tables) {
    if (!is.null(x) && nrow(x) > 0) out <- merge(out, x, by = "peak_id", all.x = TRUE, sort = FALSE)
  }
  flag_cols <- setdiff(names(out), "peak_id")
  for (col in flag_cols) out[is.na(get(col)), (col) := FALSE]
  out[]
}

expand_model_assignment_to_all_peaks <- function(
  model_dt,
  peaks_dt,
  comprehensive_flags_dt,
  repeat_peak_ids
) {
  peak_ref <- peaks_dt[, .(
    peak_id,
    chromosome,
    start,
    end,
    n_outlier_windows,
    max_Xin_Bel_Fst,
    mean_Xin_Bel_Fst
  )]

  model_extra_cols <- setdiff(names(model_dt), names(peak_ref))
  model_extra <- model_dt[, c("peak_id", model_extra_cols), with = FALSE]
  out <- merge(peak_ref, model_extra, by = "peak_id", all.x = TRUE, sort = FALSE)
  out <- merge(out, comprehensive_flags_dt, by = "peak_id", all.x = TRUE, sort = FALSE)
  out[, near_large_repeat_region_within_50kb := peak_id %in% repeat_peak_ids]

  # Requested significance columns are explicit TRUE/FALSE for every peak.
  requested_flag_cols <- c(
    grep("_(high|low)$", names(out), value = TRUE),
    "near_large_repeat_region_within_50kb"
  )
  for (col in unique(requested_flag_cols)) out[is.na(get(col)), (col) := FALSE]
  setorder(out, chromosome, start, end)
  out[]
}

make_fd_rth_correlation_outputs <- function(
  arg_points_dt,
  d_stats_file,
  all_windows_dt,
  outfile_pdf,
  outfile_data,
  outfile_summary
) {
  if (nrow(arg_points_dt) == 0 || nrow(all_windows_dt) == 0) {
    warning("No ARG points or pixy windows available for Tapajos_Xingu_fd vs Tap_Xin_RTH_original correlation.")
    return(invisible(data.table()))
  }

  # Use the pixy 10-kb windows as the reference. ARG local-tree values are
  # averaged within each pixy window, while D-statistic fd values are retained
  # only when their window coordinates exactly match a pixy window after
  # converting D-statistic coordinates to BED-style [start, end).
  win <- copy(all_windows_dt)
  win[, window_id := paste(chromosome, start, end, sep = ":")]
  win_intervals <- win[, .(
    window_id,
    chromosome,
    window_start = as.integer(start),
    window_end = as.integer(end),
    window_end_inclusive = as.integer(end - 1L)
  )]

  arg_tmp <- copy(arg_points_dt)
  arg_tmp <- arg_tmp[!is.na(Tap_Xin_RTH_original)]
  if (nrow(arg_tmp) == 0) {
    warning("No non-missing Tap_Xin_RTH_original values available for correlation plot.")
    return(invisible(data.table()))
  }

  setkey(arg_tmp, chromosome, pos_start, pos_end)
  setkey(win_intervals, chromosome, window_start, window_end_inclusive)

  arg_by_window <- foverlaps(
    arg_tmp[, .(chromosome, pos_start, pos_end, Tap_Xin_RTH_original)],
    win_intervals,
    by.x = c("chromosome", "pos_start", "pos_end"),
    by.y = c("chromosome", "window_start", "window_end_inclusive"),
    nomatch = 0
  )[, .(
    mean_Tap_Xin_RTH_original = mean(Tap_Xin_RTH_original, na.rm = TRUE),
    n_ARG_local_trees = .N
  ), by = .(window_id, chromosome, window_start, window_end)]

  d_fd <- fread(
    d_stats_file,
    select = c("scaffold", "start", "end", "sitesUsed", "D", "fd")
  )
  setnames(d_fd, old = c("scaffold", "fd"), new = c("chromosome", "Tapajos_Xingu_fd"))
  d_fd[, chromosome := standardize_scaffold_names(chromosome)]
  d_fd[!is.na(D) & D < 0, Tapajos_Xingu_fd := 0]
  # Exclude biologically invalid fd values before downstream matching/correlation.
  d_fd <- d_fd[
    is.na(Tapajos_Xingu_fd) |
      (Tapajos_Xingu_fd >= 0 & Tapajos_Xingu_fd <= 1)
  ]
  d_fd <- d_fd[sitesUsed >= 100]
  d_fd[, window_start := as.integer(start - 1L)]
  d_fd[, window_end := as.integer(end)]

  # Exact coordinate match to pixy windows. This intentionally skips any
  # D-statistic windows that differ from the pixy reference window boundaries.
  fd_by_window <- merge(
    win_intervals[, .(window_id, chromosome, window_start, window_end)],
    d_fd[, .(chromosome, window_start, window_end, Tapajos_Xingu_fd, D, sitesUsed)],
    by = c("chromosome", "window_start", "window_end"),
    all.x = FALSE,
    all.y = FALSE,
    sort = FALSE
  )

  plot_dt <- merge(
    fd_by_window,
    arg_by_window,
    by = c("window_id", "chromosome", "window_start", "window_end"),
    all = FALSE,
    sort = FALSE
  )

  plot_dt <- plot_dt[
    !is.na(Tapajos_Xingu_fd) &
      !is.na(mean_Tap_Xin_RTH_original) &
      is.finite(Tapajos_Xingu_fd) &
      is.finite(mean_Tap_Xin_RTH_original)
  ]

  if (nrow(plot_dt) < 3) {
    warning("Fewer than three matched windows available for Tapajos_Xingu_fd vs Tap_Xin_RTH_original correlation.")
    fwrite(plot_dt, outfile_data, sep = "\t")
    summary_dt <- data.table(
      comparison = "Tapajos_Xingu_fd_vs_mean_Tap_Xin_RTH_original",
      n_matched_windows = nrow(plot_dt),
      pearson_r = NA_real_,
      pearson_p_value = NA_real_,
      r_squared = NA_real_,
      slope = NA_real_,
      intercept = NA_real_,
      mean_Tapajos_Xingu_fd = mean(plot_dt$Tapajos_Xingu_fd, na.rm = TRUE),
      mean_Tap_Xin_RTH_original = mean(plot_dt$mean_Tap_Xin_RTH_original, na.rm = TRUE)
    )
    fwrite(summary_dt, outfile_summary, sep = "\t")
    return(invisible(summary_dt))
  }

  cor_test <- cor.test(
    plot_dt$Tapajos_Xingu_fd,
    plot_dt$mean_Tap_Xin_RTH_original,
    method = "pearson"
  )

  lm_fit <- lm(mean_Tap_Xin_RTH_original ~ Tapajos_Xingu_fd, data = plot_dt)
  lm_coef <- coef(lm_fit)

  summary_dt <- data.table(
    comparison = "Tapajos_Xingu_fd_vs_mean_Tap_Xin_RTH_original",
    n_matched_windows = nrow(plot_dt),
    pearson_r = unname(cor_test$estimate),
    pearson_p_value = cor_test$p.value,
    r_squared = unname(cor_test$estimate)^2,
    slope = unname(lm_coef[["Tapajos_Xingu_fd"]]),
    intercept = unname(lm_coef[["(Intercept)"]]),
    mean_Tapajos_Xingu_fd = mean(plot_dt$Tapajos_Xingu_fd, na.rm = TRUE),
    mean_Tap_Xin_RTH_original = mean(plot_dt$mean_Tap_Xin_RTH_original, na.rm = TRUE),
    median_Tapajos_Xingu_fd = median(plot_dt$Tapajos_Xingu_fd, na.rm = TRUE),
    median_Tap_Xin_RTH_original = median(plot_dt$mean_Tap_Xin_RTH_original, na.rm = TRUE)
  )

  fwrite(plot_dt, outfile_data, sep = "\t")
  fwrite(summary_dt, outfile_summary, sep = "\t")

  pdf(outfile_pdf, width = 7.4, height = 6.4, useDingbats = FALSE)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  par(mar = c(5.1, 5.2, 4.2, 1.6))
  plot(
    plot_dt$Tapajos_Xingu_fd,
    plot_dt$mean_Tap_Xin_RTH_original,
    pch = 16,
    cex = 0.45,
    col = adjustcolor("black", alpha.f = 0.28),
    xlab = "Tapajos-Xingu fd",
    ylab = "Mean Tapajos-Xingu RTH across local trees per 10-kb window",
    main = "Tapajos-Xingu fd vs. ARG-based Tapajos-Xingu RTH"
  )
  abline(lm_fit, col = "red", lwd = 1.4)
  grid(col = "grey88", lty = "dotted")
  box()

  usr <- par("usr")
  annotation <- c(
    paste0("Matched 10-kb windows = ", format(nrow(plot_dt), big.mark = ",")),
    paste0("Pearson r = ", signif(summary_dt$pearson_r, 4)),
    paste0("P = ", format.pval(summary_dt$pearson_p_value, digits = 3, eps = 1e-300)),
    paste0("R^2 = ", signif(summary_dt$r_squared, 4))
  )
  text(
    x = usr[1] + 0.03 * diff(usr[1:2]),
    y = usr[4] - 0.05 * diff(usr[3:4]),
    labels = paste(annotation, collapse = "\n"),
    adj = c(0, 1),
    cex = 0.82
  )

  invisible(summary_dt)
}

make_threshold_rows <- function(values, stat, threshold_basis, tails) {
  values <- values[!is.na(values)]

  rbindlist(lapply(tails, function(tail) {
    if (tail == "upper") {
      probs <- 0.999
      percentile <- 99.9
    } else if (tail == "lower") {
      probs <- 0.001
      percentile <- 0.1
    } else {
      stop("Unexpected tail value: ", tail)
    }

    data.table(
      statistic = stat,
      threshold_basis = threshold_basis,
      tail = tail,
      p_value = 0.001,
      percentile = percentile,
      threshold = if (length(values) > 0) {
        as.numeric(quantile(values, probs = probs, na.rm = TRUE))
      } else {
        NA_real_
      },
      n_values = length(values)
    )
  }))
}

make_value_summary <- function(values, stat, threshold_basis) {
  values <- values[!is.na(values)]

  data.table(
    statistic = stat,
    threshold_basis = threshold_basis,
    n_values = length(values),
    mean = if (length(values) > 0) mean(values) else NA_real_,
    sd = if (length(values) > 1) sd(values) else NA_real_,
    min = if (length(values) > 0) min(values) else NA_real_,
    max = if (length(values) > 0) max(values) else NA_real_
  )
}

extract_control_values_from_intervals <- function(dt, stat, control_interval_dt) {
  tmp <- copy(dt)
  tmp <- tmp[!is.na(get(stat))]
  if (nrow(tmp) == 0) {
    return(numeric(0))
  }

  tmp[, interval_start := as.integer(interval_start)]
  tmp[, interval_end_inclusive := as.integer(interval_end - 1L)]
  tmp <- tmp[interval_end_inclusive >= interval_start]

  setkey(tmp, chromosome, interval_start, interval_end_inclusive)

  ctrl <- copy(control_interval_dt)
  setkey(ctrl, chromosome, window_start, window_end_inclusive)

  joined <- foverlaps(
    tmp,
    ctrl,
    by.x = c("chromosome", "interval_start", "interval_end_inclusive"),
    by.y = c("chromosome", "window_start", "window_end_inclusive"),
    nomatch = 0
  )

  joined[[stat]]
}

extract_control_values_from_points <- function(dt, stat, control_interval_dt) {
  tmp <- copy(dt)
  tmp <- tmp[!is.na(get(stat))]
  if (nrow(tmp) == 0) {
    return(numeric(0))
  }

  tmp[, point_start := as.integer(pos)]
  tmp[, point_end := as.integer(pos)]

  setkey(tmp, chromosome, point_start, point_end)

  ctrl <- copy(control_interval_dt)
  setkey(ctrl, chromosome, window_start, window_end_inclusive)

  joined <- foverlaps(
    tmp,
    ctrl,
    by.x = c("chromosome", "point_start", "point_end"),
    by.y = c("chromosome", "window_start", "window_end_inclusive"),
    nomatch = 0
  )

  joined[[stat]]
}

make_fixed_threshold_rows <- function(stats, threshold_basis, n_values_lookup) {
  rbindlist(lapply(stats, function(stat) {
    n_values <- n_values_lookup[statistic == stat, n_values]
    if (length(n_values) == 0) n_values <- NA_integer_

    rbind(
      data.table(
        statistic = stat,
        threshold_basis = threshold_basis,
        tail = "lower",
        p_value = NA_real_,
        percentile = NA_real_,
        threshold = -2,
        n_values = n_values[1]
      ),
      data.table(
        statistic = stat,
        threshold_basis = threshold_basis,
        tail = "upper",
        p_value = NA_real_,
        percentile = NA_real_,
        threshold = 2,
        n_values = n_values[1]
      )
    )
  }), fill = TRUE)
}

mean_by_window_from_intervals <- function(dt, stat, window_interval_dt) {
  tmp <- copy(dt)
  tmp <- tmp[!is.na(get(stat))]
  if (nrow(tmp) == 0) {
    return(window_interval_dt[, .(window_id, value = NA_real_)])
  }

  tmp[, interval_start := as.integer(interval_start)]
  tmp[, interval_end_inclusive := as.integer(interval_end - 1L)]
  tmp <- tmp[interval_end_inclusive >= interval_start]

  setkey(tmp, chromosome, interval_start, interval_end_inclusive)

  win <- copy(window_interval_dt)
  setkey(win, chromosome, window_start, window_end_inclusive)

  joined <- foverlaps(
    tmp,
    win,
    by.x = c("chromosome", "interval_start", "interval_end_inclusive"),
    by.y = c("chromosome", "window_start", "window_end_inclusive"),
    nomatch = 0
  )

  out <- joined[, .(value = mean(get(stat), na.rm = TRUE)), by = window_id]
  merge(window_interval_dt[, .(window_id)], out, by = "window_id", all.x = TRUE, sort = FALSE)
}

mean_by_window_from_points <- function(dt, stat, window_interval_dt, use_abs = FALSE) {
  tmp <- copy(dt)
  tmp <- tmp[!is.na(get(stat))]
  if (nrow(tmp) == 0) {
    return(window_interval_dt[, .(window_id, value = NA_real_)])
  }

  tmp[, point_start := as.integer(pos)]
  tmp[, point_end := as.integer(pos)]

  setkey(tmp, chromosome, point_start, point_end)

  win <- copy(window_interval_dt)
  setkey(win, chromosome, window_start, window_end_inclusive)

  joined <- foverlaps(
    tmp,
    win,
    by.x = c("chromosome", "point_start", "point_end"),
    by.y = c("chromosome", "window_start", "window_end_inclusive"),
    nomatch = 0
  )

  if (use_abs) {
    out <- joined[, .(value = mean(abs(get(stat)), na.rm = TRUE)), by = window_id]
  } else {
    out <- joined[, .(value = mean(get(stat), na.rm = TRUE)), by = window_id]
  }

  merge(window_interval_dt[, .(window_id)], out, by = "window_id", all.x = TRUE, sort = FALSE)
}

add_window_stat_column <- function(window_values, stat, values_dt) {
  tmp <- values_dt[, .(window_id, value)]
  setnames(tmp, "value", stat)
  merge(window_values, tmp, by = "window_id", all.x = TRUE, sort = FALSE)
}

percentile_rank_0_1 <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x)
  n <- sum(ok)
  if (n == 0) return(out)
  if (n == 1) {
    out[ok] <- 0.5
  } else {
    out[ok] <- (rank(x[ok], ties.method = "average") - 1) / (n - 1)
  }
  out
}

make_peak_heatmap_pdf <- function(heatmap_dt, stat_order, plot_labels, outfile) {
  mat <- as.matrix(heatmap_dt[, ..stat_order])
  rownames(mat) <- paste0(
    heatmap_dt$peak_id, ": ",
    heatmap_dt$chromosome, ":",
    heatmap_dt$start, "-",
    heatmap_dt$end
  )

  n_peaks <- nrow(mat)
  n_stats <- ncol(mat)

  # Wider PDF and larger margins prevent clipping of peak labels,
  # diagonal statistic labels, and the percentile legend.
  pdf(
    outfile,
    width = max(14, n_stats * 0.45 + 3),
    height = max(7, n_peaks * 0.18 + 4.5)
  )

  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  })

  par(
    mar = c(11.5, 12.5, 4.8, 5.2),
    xpd = NA
  )

  plot(
    NA,
    xlim = c(0.5, n_stats + 0.5),
    ylim = c(0.5, n_peaks + 0.5),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    main = paste0(
      "Xingu-Belem Fst peak statistic heat map\n",
      "Total number of peaks included = ",
      n_peaks
    )
  )

  cols <- colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(101)

  for (i in seq_len(n_peaks)) {
    y <- n_peaks - i + 1
    for (j in seq_len(n_stats)) {
      val <- mat[i, j]
      col <- if (is.na(val)) {
        "grey85"
      } else {
        cols[pmin(pmax(floor(val * 100) + 1L, 1L), 101L)]
      }
      rect(j - 0.5, y - 0.5, j + 0.5, y + 0.5, col = col, border = "grey80")
    }
  }

  # Draw x-axis statistic labels manually at a diagonal angle to avoid clipping
  # and improve readability.
  axis(
    1,
    at = seq_len(n_stats),
    labels = FALSE
  )

  text(
    x = seq_len(n_stats),
    y = par("usr")[3] - 0.035 * diff(par("usr")[3:4]),
    labels = unname(plot_labels[stat_order]),
    srt = 45,
    adj = 1,
    xpd = TRUE,
    cex = 0.65
  )

  # Large left margin prevents peak labels from being clipped.
  axis(
    2,
    at = seq_len(n_peaks),
    labels = rev(rownames(mat)),
    las = 2,
    cex.axis = 0.55
  )

  box()

  # Legend is placed outside the plotting panel, with extra right margin.
  legend_x0 <- n_stats + 0.95
  legend_x1 <- n_stats + 1.25
  legend_y <- seq(0.8, n_peaks + 0.2, length.out = 101)

  if (n_peaks >= 5) {
    for (k in seq_len(100)) {
      rect(
        legend_x0,
        legend_y[k],
        legend_x1,
        legend_y[k + 1],
        col = cols[k],
        border = NA,
        xpd = TRUE
      )
    }

    rect(
      legend_x0,
      min(legend_y),
      legend_x1,
      max(legend_y),
      col = NA,
      border = "black",
      xpd = TRUE
    )

    text(legend_x1 + 0.12, min(legend_y), "0", adj = 0, cex = 0.7, xpd = TRUE)
    text(legend_x1 + 0.12, mean(range(legend_y)), "0.5", adj = 0, cex = 0.7, xpd = TRUE)
    text(legend_x1 + 0.12, max(legend_y), "1", adj = 0, cex = 0.7, xpd = TRUE)

    text(
      x = (legend_x0 + legend_x1) / 2,
      y = max(legend_y) + 0.9,
      labels = "percentile",
      adj = 0.5,
      cex = 0.75,
      xpd = TRUE
    )
  }
}



make_peak_size_distribution_barplot <- function(peaks_dt, outfile) {
  if (nrow(peaks_dt) == 0) {
    warning("No Xingu-Belem Fst peaks available for peak-size barplot: ", outfile)
    return(invisible(data.table()))
  }

  plot_dt <- copy(peaks_dt)
  plot_dt[, peak_size_bp := as.numeric(end - start)]
  plot_dt[, peak_size_kb := peak_size_bp / 1000]
  setorder(plot_dt, peak_size_bp, chromosome, start, end)
  plot_dt[, peak_rank := seq_len(.N)]

  pdf(
    outfile,
    width = 8.2,
    height = 5.6,
    useDingbats = FALSE
  )
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  par(mar = c(5.0, 5.1, 4.0, 1.2))
  barplot(
    height = plot_dt$peak_size_kb,
    names.arg = rep("", nrow(plot_dt)),
    border = NA,
    xlab = "Xingu-Belem Fst peaks ordered by size",
    ylab = "Peak size (kb)",
    main = paste0(
      "Distribution of Xingu-Belem Fst peak sizes\n",
      "All Fst peaks included = ",
      format(nrow(plot_dt), big.mark = ",")
    )
  )
  box()

  invisible(plot_dt)
}


make_model_colored_peak_size_distribution_barplot <- function(
  peaks_dt,
  model_dt,
  outfile,
  assignment_label
) {
  if (nrow(peaks_dt) == 0) {
    warning("No Xingu-Belem Fst peaks available for model-colored peak-size barplot: ", outfile)
    return(invisible(data.table()))
  }

  plot_dt <- copy(peaks_dt)
  plot_dt[, peak_size_bp := as.numeric(end - start)]
  plot_dt[, peak_size_kb := peak_size_bp / 1000]

  model_lookup <- unique(model_dt[, .(peak_id, model_assignment)])
  plot_dt <- merge(
    plot_dt,
    model_lookup,
    by = "peak_id",
    all.x = TRUE,
    sort = FALSE
  )

  # Peaks lacking an ARG-based model assignment are retained and shown in black.
  plot_dt[
    is.na(model_assignment) | !nzchar(model_assignment),
    model_assignment := "No ARG data"
  ]

  model_colors <- c(
    "Selection-bottleneck model" = "#0072B2",
    "Selection-recombination model" = "#E69F00",
    "Deep lineage sorting model" = "#009E73",
    "Introgression model" = "#CC79A7",
    "Overlapping models" = "#FFFFFF",
    "Unassigned to model" = "#808080",
    "No ARG data" = "#000000"
  )

  plot_dt[, bar_color := unname(model_colors[model_assignment])]
  plot_dt[is.na(bar_color), bar_color := "#808080"]

  setorder(plot_dt, peak_size_bp, chromosome, start, end)
  plot_dt[, peak_rank := seq_len(.N)]

  pdf(
    outfile,
    width = 9.0,
    height = 6.0,
    useDingbats = FALSE
  )
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  par(mar = c(5.0, 5.1, 4.5, 1.2))
  barplot(
    height = plot_dt$peak_size_kb,
    names.arg = rep("", nrow(plot_dt)),
    col = plot_dt$bar_color,
    border = "black",
    lwd = 0.25,
    xlab = "Xingu-Belem Fst peaks ordered by size",
    ylab = "Peak size (kb)",
    main = paste0(
      "Distribution of Xingu-Belem Fst peak sizes\n",
      "Bars colored by model assignment: ",
      assignment_label,
      "\nAll Fst peaks included = ",
      format(nrow(plot_dt), big.mark = ",")
    )
  )
  box()

  legend_order <- names(model_colors)
  legend(
    "topleft",
    legend = legend_order,
    fill = unname(model_colors[legend_order]),
    border = "black",
    bty = "n",
    cex = 0.72,
    ncol = 2
  )

  invisible(plot_dt)
}

make_model_filtered_heatmap_pdf <- function(
  heatmap_dt,
  model_dt,
  stat_order,
  plot_labels,
  outfile,
  title_text
) {
  model_order <- c(
    "Introgression model",
    "Deep lineage sorting model",
    "Selection-bottleneck model",
    "Selection-recombination model"
  )

  if (nrow(heatmap_dt) == 0 || nrow(model_dt) == 0) {
    warning("No heatmap/model data available for ", outfile)
    return(invisible(data.table()))
  }

  plot_dt <- merge(
    heatmap_dt,
    model_dt[, .(peak_id, model_assignment)],
    by = "peak_id",
    all.x = FALSE,
    all.y = FALSE,
    sort = FALSE
  )

  plot_dt <- plot_dt[model_assignment %in% model_order]
  if (nrow(plot_dt) == 0) {
    warning("No specifically assigned model peaks available for ", outfile)
    return(invisible(data.table()))
  }

  plot_dt[, model_order_index := match(model_assignment, model_order)]
  setorder(plot_dt, model_order_index, chromosome, start, end)

  rows <- list()
  row_meta <- list()
  y_labels <- data.table()

  row_i <- 0L
  for (model_i in model_order) {
    group_dt <- plot_dt[model_assignment == model_i]
    if (nrow(group_dt) == 0) next

    if (row_i > 0L) {
      row_i <- row_i + 1L
      rows[[length(rows) + 1L]] <- rep(NA_real_, length(stat_order))
      row_meta[[length(row_meta) + 1L]] <- data.table(
        row_index = row_i,
        model_assignment = NA_character_,
        is_spacer = TRUE
      )
    }

    group_start <- row_i + 1L
    for (j in seq_len(nrow(group_dt))) {
      row_i <- row_i + 1L
      rows[[length(rows) + 1L]] <- as.numeric(unlist(group_dt[j, ..stat_order], use.names = FALSE))
      row_meta[[length(row_meta) + 1L]] <- data.table(
        row_index = row_i,
        model_assignment = model_i,
        is_spacer = FALSE
      )
    }
    group_end <- row_i

    label <- paste0(gsub(" model$", "", model_i), " (", nrow(group_dt), ")")
    y_labels <- rbind(
      y_labels,
      data.table(
        model_assignment = model_i,
        label = label,
        row_midpoint = mean(c(group_start, group_end))
      ),
      fill = TRUE
    )
  }

  mat <- do.call(rbind, rows)
  n_rows <- nrow(mat)
  n_stats <- ncol(mat)

  pdf(
    outfile,
    width = 9.2,
    height = max(4.6, min(7.2, 2.3 + 0.035 * n_rows + 0.10 * n_stats)),
    useDingbats = FALSE
  )
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  })

  par(mar = c(8.0, 6.8, 3.1, 1.0), xpd = NA)

  plot(
    NA,
    xlim = c(0.5, n_stats + 0.5),
    ylim = c(0.5, n_rows + 0.5),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    main = paste0(title_text, "\nPeaks included = ", nrow(plot_dt))
  )

  cols <- colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(101)

  for (i in seq_len(n_rows)) {
    y <- n_rows - i + 1
    if (all(is.na(mat[i, ]))) {
      rect(0.5, y - 0.5, n_stats + 0.5, y + 0.5, col = "white", border = NA)
      next
    }

    for (j in seq_len(n_stats)) {
      val <- mat[i, j]
      col <- if (is.na(val)) {
        "grey88"
      } else {
        cols[pmin(pmax(floor(val * 100) + 1L, 1L), 101L)]
      }
      rect(j - 0.5, y - 0.5, j + 0.5, y + 0.5, col = col, border = NA)
    }
  }

  box()
  axis(1, at = seq_len(n_stats), labels = FALSE)
  text(
    x = seq_len(n_stats),
    y = par("usr")[3] - 0.04 * diff(par("usr")[3:4]),
    labels = unname(plot_labels[stat_order]),
    srt = 45,
    adj = 1,
    cex = 0.52,
    xpd = TRUE
  )

  if (nrow(y_labels) > 0) {
    axis(
      2,
      at = n_rows - y_labels$row_midpoint + 1,
      labels = y_labels$label,
      las = 2,
      tick = FALSE,
      cex.axis = 0.62
    )
  }

  # Compact horizontal percentile legend.
  legend_x0 <- 0.60
  legend_x1 <- min(n_stats + 0.4, 7.6)
  legend_y0 <- par("usr")[3] - 1.28
  legend_y1 <- par("usr")[3] - 0.78
  legend_x <- seq(legend_x0, legend_x1, length.out = 102)
  for (k in seq_len(101)) {
    rect(legend_x[k], legend_y0, legend_x[k + 1], legend_y1, col = cols[k], border = NA, xpd = TRUE)
  }
  rect(legend_x0, legend_y0, legend_x1, legend_y1, border = "black", col = NA, xpd = TRUE)
  text(legend_x0, legend_y0 - 0.18, "0", cex = 0.62, adj = c(0.5, 1), xpd = TRUE)
  text(mean(c(legend_x0, legend_x1)), legend_y0 - 0.18, "0.5", cex = 0.62, adj = c(0.5, 1), xpd = TRUE)
  text(legend_x1, legend_y0 - 0.18, "1", cex = 0.62, adj = c(0.5, 1), xpd = TRUE)
  text(legend_x0, legend_y1 + 0.14, "percentile", cex = 0.65, adj = c(0, 0), xpd = TRUE)

  invisible(plot_dt)
}


make_peak_percentile_pca_pdf <- function(
  heatmap_peak_representative,
  stat_order,
  plot_labels,
  model_dt,
  outfile_pdf,
  outfile_scores
) {
  if (nrow(heatmap_peak_representative) == 0) {
    warning("No heatmap peaks available for PCA; skipping: ", outfile_pdf)
    return(invisible(data.table()))
  }

  missing_stats <- setdiff(stat_order, names(heatmap_peak_representative))
  if (length(missing_stats) > 0) {
    warning(
      "The following heatmap statistics are missing from the PCA input and will be skipped: ",
      paste(missing_stats, collapse = ", ")
    )
  }

  pca_stats <- intersect(stat_order, names(heatmap_peak_representative))
  if (length(pca_stats) < 2) {
    warning("Fewer than two statistics available for PCA; skipping: ", outfile_pdf)
    return(invisible(data.table()))
  }

  mat_dt <- copy(heatmap_peak_representative[, ..pca_stats])

  for (stat in pca_stats) {
    mat_dt[[stat]] <- as.numeric(mat_dt[[stat]])

    # PCA must use a complete matrix. Because percentile values are on a
    # 0-to-1 scale, missing statistic values are imputed to 0.5, the neutral
    # midpoint of the percentile distribution. This keeps all heat-map peaks
    # in the PCA without making missing values artificially extreme.
    mat_dt[is.na(get(stat)), (stat) := 0.5]
  }

  # Remove zero-variance columns because they cannot contribute to a scaled PCA.
  stat_sds <- vapply(mat_dt, sd, numeric(1), na.rm = TRUE)
  pca_stats_kept <- names(stat_sds)[is.finite(stat_sds) & stat_sds > 0]

  if (length(pca_stats_kept) < 2) {
    warning("Fewer than two non-constant statistics available for PCA; skipping: ", outfile_pdf)
    return(invisible(data.table()))
  }

  pca_mat <- as.matrix(mat_dt[, ..pca_stats_kept])
  rownames(pca_mat) <- heatmap_peak_representative$peak_id

  pca <- prcomp(pca_mat, center = TRUE, scale. = TRUE)
  variance_explained <- (pca$sdev^2) / sum(pca$sdev^2)

  pca_scores <- data.table(
    peak_id = heatmap_peak_representative$peak_id,
    chromosome = heatmap_peak_representative$chromosome,
    start = heatmap_peak_representative$start,
    end = heatmap_peak_representative$end,
    n_outlier_windows = heatmap_peak_representative$n_outlier_windows,
    n_missing_percentile_scores = rowSums(is.na(heatmap_peak_representative[, ..pca_stats])),
    n_statistics_used_for_PCA = length(pca_stats_kept)
  )

  n_pcs_to_save <- min(4L, ncol(pca$x))
  for (pc_i in seq_len(n_pcs_to_save)) {
    pca_scores[, paste0("PC", pc_i) := pca$x[, pc_i]]
  }

  # Keep expected PC columns in the output table even if fewer components are
  # available, so downstream file structure is predictable.
  for (pc_i in seq_len(4L)) {
    pc_col <- paste0("PC", pc_i)
    if (!pc_col %in% names(pca_scores)) {
      pca_scores[, (pc_col) := NA_real_]
    }
  }

  model_cols <- c("peak_id", "model_assignment")
  if (all(model_cols %in% names(model_dt))) {
    pca_scores <- merge(
      pca_scores,
      model_dt[, ..model_cols],
      by = "peak_id",
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    pca_scores[, model_assignment := NA_character_]
  }

  pca_scores[is.na(model_assignment), model_assignment := "Unassigned to model"]

  model_colors <- c(
    "Introgression model" = "#A8DADC",
    "Deep lineage sorting model" = "#F4A261",
    "Selection-bottleneck model" = "#A3B18B",
    "Selection-recombination model" = "#D8B4E2",
    "Overlapping models" = "#808080",
    "Unassigned to model" = "#808080"
  )

  pca_scores[, point_color := model_colors[model_assignment]]
  pca_scores[is.na(point_color), point_color := "#808080"]

  fwrite(pca_scores, outfile_scores, sep = "\t")

  pdf(outfile_pdf, width = 7.2, height = 6.6, useDingbats = FALSE)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  plot_pca_pair <- function(x_pc, y_pc) {
    x_col <- paste0("PC", x_pc)
    y_col <- paste0("PC", y_pc)

    if (!all(c(x_col, y_col) %in% names(pca_scores)) ||
        all(is.na(pca_scores[[x_col]])) ||
        all(is.na(pca_scores[[y_col]]))) {
      plot.new()
      text(
        0.5,
        0.5,
        paste0(x_col, " vs ", y_col, " unavailable\nfewer than ", y_pc, " PCA axes available"),
        cex = 0.9
      )
      return(invisible(NULL))
    }

    par(mar = c(5, 5, 4.2, 2))

    plot(
      pca_scores[[x_col]],
      pca_scores[[y_col]],
      pch = 16,
      cex = 1.05,
      col = pca_scores$point_color,
      xlab = paste0(
        x_col,
        " (",
        round(variance_explained[x_pc] * 100, 1),
        "% variance explained)"
      ),
      ylab = paste0(
        y_col,
        " (",
        round(variance_explained[y_pc] * 100, 1),
        "% variance explained)"
      ),
      main = paste0(
        "PCA of Xingu-Belem Fst peak statistic percentile profiles\n",
        x_col,
        " vs ",
        y_col,
        "; peaks included = ",
        nrow(pca_scores)
      )
    )

    legend_categories <- c(
      "Introgression model",
      "Deep lineage sorting model",
      "Selection-bottleneck model",
      "Selection-recombination model"
    )

    legend(
      "topright",
      legend = legend_categories,
      col = model_colors[legend_categories],
      pch = 16,
      pt.cex = 1.05,
      bty = "n",
      cex = 0.78
    )

    box()
    invisible(NULL)
  }

  # Each PCA PDF contains three pages:
  #   1) PC1 vs PC2
  #   2) PC1 vs PC3
  #   3) PC1 vs PC4
  plot_pca_pair(1L, 2L)
  plot_pca_pair(1L, 3L)
  plot_pca_pair(1L, 4L)

  invisible(pca_scores)
}

rolling_mean_centered <- function(x, half_window = 5L) {
  n <- length(x)
  out <- rep(NA_real_, n)

  if (n == 0) {
    return(out)
  }

  for (i in seq_len(n)) {
    lo <- max(1L, i - half_window)
    hi <- min(n, i + half_window)
    vals <- x[lo:hi]
    if (all(is.na(vals))) {
      out[i] <- NA_real_
    } else {
      out[i] <- mean(vals, na.rm = TRUE)
    }
  }

  out
}


# Deprecated helper retained for compatibility but no longer used.
# Manhattan plots now draw all points as vector black circles.
make_point_raster <- function(
  x,
  y,
  xlim,
  ylim,
  width_px = 2200L,
  height_px = 320L,
  point_radius_px = 2.4,
  point_col = "#000000"
) {
  ok <- !is.na(x) & !is.na(y) &
    x >= xlim[1] & x <= xlim[2] &
    y >= ylim[1] & y <= ylim[2]

  x <- x[ok]
  y <- y[ok]

  # Store alpha values separately so overlapping points become fully black.
  alpha_mat <- matrix(0, nrow = height_px, ncol = width_px)

  if (length(x) == 0) {
    return(as.raster(matrix("#FFFFFF00", nrow = height_px, ncol = width_px)))
  }

  x_pixel <- as.integer(round((x - xlim[1]) / diff(xlim) * (width_px - 1L))) + 1L
  y_pixel <- as.integer(round((y - ylim[1]) / diff(ylim) * (height_px - 1L))) + 1L

  x_pixel <- pmin(pmax(x_pixel, 1L), width_px)
  y_pixel <- pmin(pmax(y_pixel, 1L), height_px)

  # Raster row 1 is the top of the image, so invert y.
  row_pixel <- height_px - y_pixel + 1L

  max_offset <- ceiling(point_radius_px + 1)
  offsets <- expand.grid(
    dx = seq.int(-max_offset, max_offset),
    dy = seq.int(-max_offset, max_offset)
  )
  offsets$distance <- sqrt(offsets$dx^2 + offsets$dy^2)

  # Anti-aliased filled circle. Pixels inside the radius are opaque;
  # pixels at the edge are partially transparent. This avoids the
  # diamond/cross-like appearance of coarse binary raster points.
  offsets$alpha <- pmin(
    1,
    pmax(0, point_radius_px + 0.5 - offsets$distance)
  )
  offsets <- offsets[offsets$alpha > 0, ]

  for (offset_i in seq_len(nrow(offsets))) {
    xx <- x_pixel + offsets$dx[offset_i]
    yy <- row_pixel + offsets$dy[offset_i]

    keep <- xx >= 1L & xx <= width_px & yy >= 1L & yy <= height_px
    if (!any(keep)) next

    alpha_mat[cbind(yy[keep], xx[keep])] <- pmax(
      alpha_mat[cbind(yy[keep], xx[keep])],
      offsets$alpha[offset_i]
    )
  }

  alpha_hex <- sprintf("%02X", pmin(255L, pmax(0L, round(alpha_mat * 255))))
  raster_mat <- paste0(point_col, alpha_hex)

  as.raster(raster_mat)
}


read_tree_folder_coordinate_map <- function(tree_map_file) {
  if (!file.exists(tree_map_file)) {
    warning("Tree folder coordinate map not found: ", tree_map_file)
    return(data.table())
  }

  tree_map <- fread(tree_map_file)

  expected_cols <- c(
    "scaffold",
    "tree_subfolder",
    "start_1based_inclusive",
    "end_1based_inclusive"
  )
  missing_cols <- setdiff(expected_cols, names(tree_map))
  if (length(missing_cols) > 0) {
    stop(
      "tree_folder_coordinate_map.tsv is missing expected column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  tree_map[, scaffold := standardize_scaffold_names(scaffold)]
  tree_map[, start_1based_inclusive := as.integer(start_1based_inclusive)]
  tree_map[, end_1based_inclusive := as.integer(end_1based_inclusive)]
  tree_map[]
}

find_tree_files_for_window <- function(
  chromosome,
  window_start,
  window_end,
  tree_map,
  tree_base_dir,
  iteration = 2000L
) {
  if (nrow(tree_map) == 0 || is.na(window_start) || is.na(window_end)) {
    return(data.table())
  }

  # Pixy/Fst windows are handled as BED-style [start, end).
  # Tree-coordinate map entries are 1-based inclusive, so overlap is:
  #   map_start <= window_end - 1 and map_end >= window_start
  #
  # Important: a single 10-kb window can overlap more than one ARG tree block.
  # In that case, return all matching .tre.gz files so random sampling can use
  # all local trees available in the source window.
  candidates <- tree_map[
    scaffold == chromosome &
      start_1based_inclusive <= as.integer(window_end - 1L) &
      end_1based_inclusive >= as.integer(window_start)
  ]

  if (nrow(candidates) == 0) {
    return(data.table())
  }

  setorder(candidates, start_1based_inclusive)

  out <- rbindlist(lapply(seq_len(nrow(candidates)), function(i) {
    row <- candidates[i]

    expected_file <- file.path(
      tree_base_dir,
      row$tree_subfolder,
      paste0(
        row$tree_subfolder,
        ".",
        row$start_1based_inclusive,
        "-",
        row$end_1based_inclusive,
        ".",
        iteration,
        ".tre.gz"
      )
    )

    if (file.exists(expected_file)) {
      tree_file <- expected_file
    } else {
      fallback <- list.files(
        file.path(tree_base_dir, row$tree_subfolder),
        pattern = paste0("\\.", iteration, "\\.tre\\.gz$"),
        full.names = TRUE
      )

      if (length(fallback) == 0) {
        warning(
          "No iteration ", iteration, " .tre.gz file found for ",
          chromosome, ":", window_start, "-", window_end,
          " in ", row$tree_subfolder
        )
        tree_file <- NA_character_
      } else {
        tree_file <- fallback[1]
      }
    }

    data.table(
      chromosome = chromosome,
      tree_subfolder = row$tree_subfolder,
      tree_block_start = row$start_1based_inclusive,
      tree_block_end = row$end_1based_inclusive,
      tree_file = tree_file
    )
  }), fill = TRUE)

  out[!is.na(tree_file) & file.exists(tree_file)]
}

simplify_tree_tip_labels <- function(phy) {
  if (is.null(phy) || is.null(phy$tip.label)) {
    return(phy)
  }

  phy$tip.label <- ifelse(
    grepl("^Tap", phy$tip.label),
    "Tapajos",
    ifelse(
      grepl("^Xin", phy$tip.label),
      "Xingu",
      ifelse(grepl("^Bel", phy$tip.label), "Belem", phy$tip.label)
    )
  )

  phy
}

sample_trees_from_window <- function(
  chromosome,
  window_start,
  window_end,
  window_label,
  n_trees,
  tree_map,
  tree_base_dir,
  iteration = 2000L
) {
  tree_files <- find_tree_files_for_window(
    chromosome = chromosome,
    window_start = window_start,
    window_end = window_end,
    tree_map = tree_map,
    tree_base_dir = tree_base_dir,
    iteration = iteration
  )

  if (nrow(tree_files) == 0) {
    return(data.table(
      tree_role = window_label,
      chromosome = chromosome,
      tree_position = NA_integer_,
      newick = NA_character_,
      tree_file = NA_character_,
      n_local_trees_in_window = 0L
    )[0])
  }

  tre_window_list <- lapply(seq_len(nrow(tree_files)), function(i) {
    tree_file_i <- tree_files$tree_file[i]

    tre <- fread(
      tree_file_i,
      header = FALSE,
      sep = "\t",
      col.names = c("pos", "newick")
    )

    tre[, pos := as.integer(pos)]

    # Local trees are selected from the BED-style window interval [start, end).
    tre <- tre[
      !is.na(pos) &
        pos >= as.integer(window_start) &
        pos < as.integer(window_end)
    ]

    if (nrow(tre) == 0) {
      return(data.table())
    }

    tre[, tree_file := tree_file_i]
    tre[]
  })

  tre_window <- rbindlist(tre_window_list, fill = TRUE)

  if (nrow(tre_window) == 0) {
    warning(
      "No local trees found in ", window_label, " for ",
      chromosome, ":", window_start, "-", window_end,
      " across ", nrow(tree_files), " overlapping tree block(s)."
    )
    return(data.table(
      tree_role = window_label,
      chromosome = chromosome,
      tree_position = NA_integer_,
      newick = NA_character_,
      tree_file = NA_character_,
      n_local_trees_in_window = 0L
    )[0])
  }

  n_local_trees <- nrow(tre_window)
  n_to_sample <- min(n_trees, n_local_trees)
  sampled_idx <- sample(seq_len(n_local_trees), size = n_to_sample, replace = FALSE)
  sampled <- tre_window[sampled_idx]

  data.table(
    tree_role = window_label,
    chromosome = chromosome,
    tree_position = sampled$pos,
    newick = sampled$newick,
    tree_file = sampled$tree_file,
    n_local_trees_in_window = n_local_trees
  )
}

sample_manhattan_panel_trees <- function(
  chromosome,
  focal_window_midpoint,
  nonfocal_window_1_midpoint,
  nonfocal_window_2_midpoint,
  manhattan_window_reference,
  tree_map,
  tree_base_dir,
  iteration = 2000L
) {
  if (nrow(tree_map) == 0) {
    return(data.table())
  }

  chromosome_i <- chromosome

  lookup_window <- function(midpoint) {
    if (is.na(midpoint)) {
      return(data.table())
    }

    manhattan_window_reference[
      chromosome == chromosome_i &
        abs(window_midpoint - midpoint) < 1e-6
    ][1]
  }

  focal_window <- lookup_window(focal_window_midpoint)
  nonfocal_window_1 <- lookup_window(nonfocal_window_1_midpoint)
  nonfocal_window_2 <- lookup_window(nonfocal_window_2_midpoint)

  out <- list()

  if (nrow(focal_window) > 0) {
    out[[length(out) + 1L]] <- sample_trees_from_window(
      chromosome = chromosome_i,
      window_start = focal_window$window_start,
      window_end = focal_window$window_end,
      window_label = "focal window",
      n_trees = 2L,
      tree_map = tree_map,
      tree_base_dir = tree_base_dir,
      iteration = iteration
    )
  }

  if (nrow(nonfocal_window_1) > 0) {
    out[[length(out) + 1L]] <- sample_trees_from_window(
      chromosome = chromosome_i,
      window_start = nonfocal_window_1$window_start,
      window_end = nonfocal_window_1$window_end,
      window_label = "non-focal window 1",
      n_trees = 1L,
      tree_map = tree_map,
      tree_base_dir = tree_base_dir,
      iteration = iteration
    )
  }

  if (nrow(nonfocal_window_2) > 0) {
    out[[length(out) + 1L]] <- sample_trees_from_window(
      chromosome = chromosome_i,
      window_start = nonfocal_window_2$window_start,
      window_end = nonfocal_window_2$window_end,
      window_label = "non-focal window 2",
      n_trees = 1L,
      tree_map = tree_map,
      tree_base_dir = tree_base_dir,
      iteration = iteration
    )
  }

  tree_dt <- rbindlist(out, fill = TRUE)
  if (nrow(tree_dt) == 0) {
    return(tree_dt)
  }

  tree_dt[, tree_panel_label := paste0(
    tree_role,
    "\nrandom tree sampled from ",
    n_local_trees_in_window,
    " local coalescence trees in window",
    "\n",
    chromosome,
    ":",
    tree_position
  )]

  tree_dt[]
}

plot_single_tree_panel <- function(tree_row) {
  if (nrow(tree_row) == 0 ||
      is.na(tree_row$newick[1]) ||
      !nzchar(tree_row$newick[1])) {
    plot.new()
    text(0.5, 0.5, "No tree sampled", cex = 0.8)
    return(invisible(NULL))
  }

  title_lines <- strsplit(tree_row$tree_panel_label[1], "\n", fixed = TRUE)[[1]]
  if (length(title_lines) < 3) {
    title_lines <- c(title_lines, rep("", 3 - length(title_lines)))
  }

  if (!requireNamespace("ape", quietly = TRUE)) {
    plot.new()
    text(0.5, 0.68, paste(title_lines, collapse = "\n"), cex = 0.78)
    text(0.5, 0.34, "Install ape to plot tree", cex = 0.75)
    return(invisible(NULL))
  }

  phy <- tryCatch(
    ape::read.tree(text = tree_row$newick[1]),
    error = function(e) NULL
  )

  if (is.null(phy)) {
    plot.new()
    text(0.5, 0.68, paste(title_lines, collapse = "\n"), cex = 0.78)
    text(0.5, 0.34, "Tree parse failed", cex = 0.75)
    return(invisible(NULL))
  }

  phy <- simplify_tree_tip_labels(phy)

  # Keep labels and tree in the same plotting panel. This avoids the fragile
  # par(fig = par("fig"), new = TRUE) overlay that can fail inside layout()
  # and can clip title text. The enlarged top margin gives the three-line
  # tree label enough room, while xpd = NA lets labels extend slightly beyond
  # the plotting region if needed.
  old_xpd <- par(xpd = NA)
  old_mar <- par(mar = c(0.8, 0.8, 5.8, 0.8))
  on.exit({
    par(xpd = old_xpd)
    par(mar = old_mar)
  }, add = TRUE)

  ape::plot.phylo(
    phy,
    type = "phylogram",
    show.tip.label = TRUE,
    cex = 0.55,
    edge.width = 0.55,
    no.margin = FALSE
  )

  mtext(title_lines[1], side = 3, line = 3.7, cex = 0.55, font = 1)
  mtext(title_lines[2], side = 3, line = 2.7, cex = 0.48, font = 1)
  mtext(title_lines[3], side = 3, line = 1.7, cex = 0.50, font = 1)

  invisible(NULL)
}


wrap_axis_label <- function(label, width = 20L) {
  if (is.na(label) || !nzchar(label)) {
    return(label)
  }

  paste(strwrap(label, width = width), collapse = "\n")
}

format_scaffold_axis_label <- function(chromosome) {
  label <- gsub("^scaffold_", "Scaffold ", chromosome)
  paste0(label, " position (Mb)")
}

standardize_scaffold_names <- function(x) {
  x <- sub("-.*$", "", x)
  sub("^scaffold(?!_)", "scaffold_", x, perl = TRUE)
}

# Read repeat regions for visual annotation in peak-by-peak Manhattan plots.
# repeat_regions.bed has a header and standard BED-style coordinates:
#   chrom  chromStart  chromEnd
# BED intervals are treated as [chromStart, chromEnd).
read_repeat_regions_for_plotting <- function(repeat_regions_file) {
  if (!file.exists(repeat_regions_file)) {
    warning("Repeat-region BED file not found: ", repeat_regions_file)
    return(data.table())
  }

  repeat_regions <- fread(repeat_regions_file)

  expected_cols <- c("chrom", "chromStart", "chromEnd")
  missing_cols <- setdiff(expected_cols, names(repeat_regions))
  if (length(missing_cols) > 0) {
    stop(
      "repeat_regions.bed is missing expected column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  repeat_regions <- repeat_regions[, .(
    chromosome = standardize_scaffold_names(chrom),
    repeat_start = as.integer(chromStart),
    repeat_end = as.integer(chromEnd)
  )]

  repeat_regions <- repeat_regions[
    !is.na(chromosome) &
      !is.na(repeat_start) &
      !is.na(repeat_end) &
      repeat_end > repeat_start
  ]

  setkey(repeat_regions, chromosome, repeat_start, repeat_end)
  repeat_regions[]
}

merge_contiguous_repeat_regions <- function(
  repeat_regions,
  min_size_bp = 50000L,
  max_gap_bp = 1000L
) {
  empty_out <- data.table(
    chromosome = character(),
    repeat_block_start = integer(),
    repeat_block_end = integer(),
    repeat_block_size_bp = integer(),
    n_repeat_intervals_in_block = integer(),
    n_gaps_merged_in_block = integer(),
    max_gap_merged_bp = integer(),
    total_gap_merged_bp = integer()
  )

  if (nrow(repeat_regions) == 0) {
    return(empty_out)
  }

  repeat_blocks <- copy(repeat_regions)
  repeat_blocks <- repeat_blocks[
    !is.na(chromosome) &
      !is.na(repeat_start) &
      !is.na(repeat_end) &
      repeat_end > repeat_start
  ]

  if (nrow(repeat_blocks) == 0) {
    return(empty_out)
  }

  setorder(repeat_blocks, chromosome, repeat_start, repeat_end)

  # Merge repeat intervals into larger repeat blocks when they overlap,
  # directly abut, or are separated by a small non-repeat gap. BED coordinates
  # are [start, end), so gap_bp = next_start - current_end. A gap of 0 means
  # intervals are directly adjacent. Setting max_gap_bp = 1000 allows very
  # nearly continuous repeat regions to be treated as one repeat block before
  # applying the >=50-kb size filter.
  merged_list <- list()
  out_i <- 0L

  for (chr_i in unique(repeat_blocks$chromosome)) {
    chr_dt <- repeat_blocks[chromosome == chr_i]
    current_start <- chr_dt$repeat_start[1]
    current_end <- chr_dt$repeat_end[1]
    current_n <- 1L
    current_n_gaps <- 0L
    current_max_gap <- 0L
    current_total_gap <- 0L

    if (nrow(chr_dt) > 1L) {
      for (row_i in 2:nrow(chr_dt)) {
        next_start <- chr_dt$repeat_start[row_i]
        next_end <- chr_dt$repeat_end[row_i]
        gap_bp <- max(0L, as.integer(next_start - current_end))

        if (gap_bp <= max_gap_bp) {
          if (gap_bp > 0L) {
            current_n_gaps <- current_n_gaps + 1L
            current_max_gap <- max(current_max_gap, gap_bp)
            current_total_gap <- current_total_gap + gap_bp
          }
          current_end <- max(current_end, next_end)
          current_n <- current_n + 1L
        } else {
          out_i <- out_i + 1L
          merged_list[[out_i]] <- data.table(
            chromosome = chr_i,
            repeat_block_start = current_start,
            repeat_block_end = current_end,
            repeat_block_size_bp = current_end - current_start,
            n_repeat_intervals_in_block = current_n,
            n_gaps_merged_in_block = current_n_gaps,
            max_gap_merged_bp = current_max_gap,
            total_gap_merged_bp = current_total_gap
          )

          current_start <- next_start
          current_end <- next_end
          current_n <- 1L
          current_n_gaps <- 0L
          current_max_gap <- 0L
          current_total_gap <- 0L
        }
      }
    }

    out_i <- out_i + 1L
    merged_list[[out_i]] <- data.table(
      chromosome = chr_i,
      repeat_block_start = current_start,
      repeat_block_end = current_end,
      repeat_block_size_bp = current_end - current_start,
      n_repeat_intervals_in_block = current_n,
      n_gaps_merged_in_block = current_n_gaps,
      max_gap_merged_bp = current_max_gap,
      total_gap_merged_bp = current_total_gap
    )
  }

  merged <- rbindlist(merged_list, fill = TRUE)
  merged <- merged[repeat_block_size_bp >= min_size_bp]
  setorder(merged, chromosome, repeat_block_start, repeat_block_end)
  merged[]
}

make_large_repeat_peak_proximity_table <- function(
  peaks_dt,
  repeat_regions,
  outfile,
  min_repeat_size_bp = 50000L,
  max_distance_bp = 50000L,
  max_repeat_gap_bp = 1000L
) {
  if (nrow(peaks_dt) == 0) {
    warning("No Xingu-Belem Fst peaks available for repeat-proximity table: ", outfile)
    empty <- data.table()
    fwrite(empty, outfile, sep = "\t")
    return(invisible(empty))
  }

  large_repeats <- merge_contiguous_repeat_regions(
    repeat_regions,
    min_size_bp = min_repeat_size_bp,
    max_gap_bp = max_repeat_gap_bp
  )

  if (nrow(large_repeats) == 0) {
    warning(
      "No contiguous repeat blocks >= ",
      min_repeat_size_bp,
      " bp were found after merging repeat intervals separated by <= ",
      max_repeat_gap_bp,
      " bp; writing an empty repeat-proximity table."
    )
    empty <- data.table()
    fwrite(empty, outfile, sep = "\t")
    return(invisible(empty))
  }

  peak_ref <- copy(peaks_dt)
  peak_ref[, chromosome := standardize_scaffold_names(chromosome)]
  setorder(peak_ref, chromosome, start, end)

  proximity_list <- lapply(seq_len(nrow(peak_ref)), function(i) {
    peak_i <- peak_ref[i]
    chr_repeats <- large_repeats[chromosome == peak_i$chromosome]

    upstream <- chr_repeats[
      repeat_block_end <= peak_i$start &
        (peak_i$start - repeat_block_end) <= max_distance_bp
    ]
    if (nrow(upstream) > 0) {
      upstream[, distance_to_peak_bp := peak_i$start - repeat_block_end]
      setorder(upstream, distance_to_peak_bp, -repeat_block_size_bp)
      upstream <- upstream[1]
    }

    downstream <- chr_repeats[
      repeat_block_start >= peak_i$end &
        (repeat_block_start - peak_i$end) <= max_distance_bp
    ]
    if (nrow(downstream) > 0) {
      downstream[, distance_to_peak_bp := repeat_block_start - peak_i$end]
      setorder(downstream, distance_to_peak_bp, -repeat_block_size_bp)
      downstream <- downstream[1]
    }

    overlapping <- chr_repeats[
      repeat_block_start < peak_i$end &
        repeat_block_end > peak_i$start
    ]
    if (nrow(overlapping) > 0) {
      overlapping[, overlap_bp := pmin(repeat_block_end, peak_i$end) -
                    pmax(repeat_block_start, peak_i$start)]
      setorder(overlapping, -overlap_bp, -repeat_block_size_bp)
      overlapping <- overlapping[1]
    }

    if (nrow(upstream) == 0 && nrow(downstream) == 0 && nrow(overlapping) == 0) {
      return(NULL)
    }

    data.table(
      peak_id = peak_i$peak_id,
      chromosome = peak_i$chromosome,
      peak_start = peak_i$start,
      peak_end = peak_i$end,
      peak_size_bp = peak_i$end - peak_i$start,
      n_outlier_windows = peak_i$n_outlier_windows,
      max_Xin_Bel_Fst = peak_i$max_Xin_Bel_Fst,
      mean_Xin_Bel_Fst = peak_i$mean_Xin_Bel_Fst,

      upstream_repeat_start = if (nrow(upstream) > 0) upstream$repeat_block_start else NA_integer_,
      upstream_repeat_end = if (nrow(upstream) > 0) upstream$repeat_block_end else NA_integer_,
      upstream_repeat_size_bp = if (nrow(upstream) > 0) upstream$repeat_block_size_bp else NA_integer_,
      upstream_repeat_distance_to_peak_bp = if (nrow(upstream) > 0) upstream$distance_to_peak_bp else NA_integer_,
      upstream_repeat_n_intervals_merged = if (nrow(upstream) > 0) upstream$n_repeat_intervals_in_block else NA_integer_,
      upstream_repeat_n_gaps_merged = if (nrow(upstream) > 0) upstream$n_gaps_merged_in_block else NA_integer_,
      upstream_repeat_max_gap_merged_bp = if (nrow(upstream) > 0) upstream$max_gap_merged_bp else NA_integer_,
      upstream_repeat_total_gap_merged_bp = if (nrow(upstream) > 0) upstream$total_gap_merged_bp else NA_integer_,

      downstream_repeat_start = if (nrow(downstream) > 0) downstream$repeat_block_start else NA_integer_,
      downstream_repeat_end = if (nrow(downstream) > 0) downstream$repeat_block_end else NA_integer_,
      downstream_repeat_size_bp = if (nrow(downstream) > 0) downstream$repeat_block_size_bp else NA_integer_,
      downstream_repeat_distance_to_peak_bp = if (nrow(downstream) > 0) downstream$distance_to_peak_bp else NA_integer_,
      downstream_repeat_n_intervals_merged = if (nrow(downstream) > 0) downstream$n_repeat_intervals_in_block else NA_integer_,
      downstream_repeat_n_gaps_merged = if (nrow(downstream) > 0) downstream$n_gaps_merged_in_block else NA_integer_,
      downstream_repeat_max_gap_merged_bp = if (nrow(downstream) > 0) downstream$max_gap_merged_bp else NA_integer_,
      downstream_repeat_total_gap_merged_bp = if (nrow(downstream) > 0) downstream$total_gap_merged_bp else NA_integer_,

      overlapping_repeat_start = if (nrow(overlapping) > 0) overlapping$repeat_block_start else NA_integer_,
      overlapping_repeat_end = if (nrow(overlapping) > 0) overlapping$repeat_block_end else NA_integer_,
      overlapping_repeat_size_bp = if (nrow(overlapping) > 0) overlapping$repeat_block_size_bp else NA_integer_,
      overlapping_repeat_overlap_with_peak_bp = if (nrow(overlapping) > 0) overlapping$overlap_bp else NA_integer_,
      overlapping_repeat_n_intervals_merged = if (nrow(overlapping) > 0) overlapping$n_repeat_intervals_in_block else NA_integer_,
      overlapping_repeat_n_gaps_merged = if (nrow(overlapping) > 0) overlapping$n_gaps_merged_in_block else NA_integer_,
      overlapping_repeat_max_gap_merged_bp = if (nrow(overlapping) > 0) overlapping$max_gap_merged_bp else NA_integer_,
      overlapping_repeat_total_gap_merged_bp = if (nrow(overlapping) > 0) overlapping$total_gap_merged_bp else NA_integer_,

      has_upstream_large_repeat_within_50kb = nrow(upstream) > 0,
      has_downstream_large_repeat_within_50kb = nrow(downstream) > 0,
      has_overlapping_large_repeat = nrow(overlapping) > 0
    )
  })

  out <- rbindlist(proximity_list, fill = TRUE)

  if (nrow(out) > 0) {
    out[, any_large_repeat_within_or_near_peak :=
          has_upstream_large_repeat_within_50kb |
          has_downstream_large_repeat_within_50kb |
          has_overlapping_large_repeat]
    setorder(out, chromosome, peak_start, peak_end)
  }

  fwrite(out, outfile, sep = "\t")
  invisible(out)
}


make_large_repeat_adjacent_Fst_manhattan_figure <- function(
  repeat_peak_dt,
  arg_peak_coverage_dt,
  fst_windows_dt,
  repeat_regions,
  scaffold_lengths_dt,
  fst_threshold,
  outfile,
  panel_width_bp = 2000000L
) {
  if (nrow(repeat_peak_dt) == 0 || nrow(arg_peak_coverage_dt) == 0) {
    warning("No repeat-adjacent peaks or ARG coverage data available for: ", outfile)
    return(invisible(data.table()))
  }

  plot_peaks <- copy(repeat_peak_dt)

  if (!"has_any_ARG_based_data" %in% names(plot_peaks)) {
    plot_peaks <- merge(
      plot_peaks,
      arg_peak_coverage_dt[, .(
        peak_id,
        n_ARG_local_trees,
        n_ARG_CC_local_trees,
        has_ARG_data,
        has_ARG_CC_data,
        has_any_ARG_based_data
      )],
      by = "peak_id",
      all.x = TRUE,
      sort = FALSE
    )
  }

  plot_peaks[is.na(n_ARG_local_trees), n_ARG_local_trees := 0L]
  plot_peaks[is.na(n_ARG_CC_local_trees), n_ARG_CC_local_trees := 0L]
  plot_peaks[is.na(has_ARG_data), has_ARG_data := FALSE]
  plot_peaks[is.na(has_ARG_CC_data), has_ARG_CC_data := FALSE]
  plot_peaks[is.na(has_any_ARG_based_data), has_any_ARG_based_data := FALSE]

  plot_peaks <- plot_peaks[has_any_ARG_based_data == TRUE]
  if (nrow(plot_peaks) == 0) {
    warning("No repeat-adjacent Xingu-Belem Fst peaks have ARG-based data; skipping: ", outfile)
    return(invisible(plot_peaks))
  }

  plot_peaks[, chromosome := standardize_scaffold_names(chromosome)]
  setorder(plot_peaks, chromosome, peak_start, peak_end)

  fst_dt <- copy(fst_windows_dt)
  fst_dt[, chromosome := standardize_scaffold_names(chromosome)]
  fst_dt[, x := (as.numeric(start) + as.numeric(end)) / 2]
  fst_dt <- fst_dt[!is.na(Xin_Bel_Fst) & is.finite(Xin_Bel_Fst)]

  panel_width_bp <- as.integer(panel_width_bp)
  half_width <- panel_width_bp / 2

  panel_bounds <- lapply(seq_len(nrow(plot_peaks)), function(i) {
    peak_i <- plot_peaks[i]
    scaffold_end <- scaffold_lengths_dt[
      chromosome == peak_i$chromosome,
      scaffold_max_end
    ]
    if (length(scaffold_end) == 0 || !is.finite(scaffold_end[1])) {
      scaffold_end <- max(fst_dt[chromosome == peak_i$chromosome]$end, na.rm = TRUE)
    } else {
      scaffold_end <- scaffold_end[1]
    }

    peak_mid <- (peak_i$peak_start + peak_i$peak_end) / 2
    panel_start <- max(0, floor(peak_mid - half_width))
    panel_end <- panel_start + panel_width_bp

    if (panel_end > scaffold_end) {
      panel_end <- scaffold_end
      panel_start <- max(0, panel_end - panel_width_bp)
    }

    data.table(
      peak_id = peak_i$peak_id,
      panel_start = as.numeric(panel_start),
      panel_end = as.numeric(panel_end)
    )
  })
  panel_bounds <- rbindlist(panel_bounds)
  plot_peaks <- merge(plot_peaks, panel_bounds, by = "peak_id", all.x = TRUE, sort = FALSE)
  setorder(plot_peaks, chromosome, peak_start, peak_end)

  n_panels <- nrow(plot_peaks)
  pdf(
    outfile,
    width = 8.3,
    height = max(6.0, 2.35 * n_panels),
    useDingbats = FALSE
  )
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  layout(matrix(seq_len(n_panels), ncol = 1L), heights = rep(1, n_panels))

  for (i in seq_len(n_panels)) {
    peak_i <- plot_peaks[i]
    panel_dt <- fst_dt[
      chromosome == peak_i$chromosome &
        x >= peak_i$panel_start &
        x <= peak_i$panel_end
    ]
    setorder(panel_dt, x)

    y_values <- panel_dt$Xin_Bel_Fst
    y_all <- c(y_values, fst_threshold)
    y_all <- y_all[is.finite(y_all)]
    if (length(y_all) == 0) y_all <- c(0, 1)
    y_min <- min(0, min(y_all))
    y_max <- max(y_all)
    if (y_max <= y_min) y_max <- y_min + 1
    y_pad <- 0.08 * (y_max - y_min)
    y_limits <- c(y_min, y_max + y_pad)

    par(mar = c(4.4, 5.5, 3.0, 1.1), xpd = FALSE)

    plot(
      NA,
      xlim = c(peak_i$panel_start, peak_i$panel_end),
      ylim = y_limits,
      xaxt = "n",
      xlab = "",
      ylab = "Xingu-Belem Fst",
      cex.lab = 0.82,
      cex.axis = 0.72,
      main = paste0(
        "Peak ", peak_i$peak_id, ": ",
        gsub("^scaffold_", "Scaffold ", peak_i$chromosome),
        " ",
        format(peak_i$peak_start, big.mark = ","),
        "-",
        format(peak_i$peak_end, big.mark = ",")
      ),
      cex.main = 0.82
    )

    draw_repeat_region_bars(
      repeat_regions = repeat_regions,
      chromosome_i = peak_i$chromosome,
      x_limits_bp = c(peak_i$panel_start, peak_i$panel_end),
      y_limits = y_limits
    )

    if (nrow(panel_dt) > 0) {
      points(
        panel_dt$x,
        panel_dt$Xin_Bel_Fst,
        pch = 16,
        cex = 0.72,
        col = "black"
      )
      if (nrow(panel_dt) >= 2L) {
        lines(
          panel_dt$x,
          rolling_mean_centered(panel_dt$Xin_Bel_Fst, half_window = 5L),
          col = "red",
          lwd = 1.0
        )
      }
    }

    if (is.finite(fst_threshold)) {
      abline(h = fst_threshold, lty = 2, lwd = 0.9)
    }

    # Mark the focal peak boundaries without introducing focal/non-focal windows.
    abline(v = c(peak_i$peak_start, peak_i$peak_end), lty = 3, lwd = 0.75)

    ticks_bp <- pretty(c(peak_i$panel_start, peak_i$panel_end), n = 5)
    ticks_bp <- ticks_bp[
      ticks_bp >= peak_i$panel_start &
        ticks_bp <= peak_i$panel_end
    ]
    axis(
      1,
      at = ticks_bp,
      labels = format(round(ticks_bp / 1e6, 2), trim = TRUE, nsmall = 2),
      cex.axis = 0.70
    )
    mtext(
      paste0(gsub("^scaffold_", "Scaffold ", peak_i$chromosome), " position (MB)"),
      side = 1,
      line = 2.7,
      cex = 0.78
    )
    box()
  }

  invisible(plot_peaks)
}

draw_repeat_region_bars <- function(
  repeat_regions,
  chromosome_i,
  x_limits_bp,
  y_limits,
  bar_height_fraction = 0.055,
  bar_offset_fraction = 0.018,
  fill_col = adjustcolor("gray30", alpha.f = 0.24)
) {
  if (is.null(repeat_regions) || nrow(repeat_regions) == 0) {
    return(invisible(NULL))
  }

  reps <- repeat_regions[
    chromosome == chromosome_i &
      repeat_end > x_limits_bp[1] &
      repeat_start < x_limits_bp[2]
  ]

  if (nrow(reps) == 0) {
    return(invisible(NULL))
  }

  y_range <- diff(y_limits)
  y_bottom <- y_limits[1] + bar_offset_fraction * y_range
  y_top <- y_bottom + bar_height_fraction * y_range

  # Draw first so repeat annotations stay behind data points, rolling means,
  # thresholds, and focal/non-focal guide lines.
  for (rep_i in seq_len(nrow(reps))) {
    rect(
      xleft = max(reps$repeat_start[rep_i], x_limits_bp[1]),
      ybottom = y_bottom,
      xright = min(reps$repeat_end[rep_i], x_limits_bp[2]),
      ytop = y_top,
      col = fill_col,
      border = NA
    )
  }

  invisible(NULL)
}

# Read ARGweaver mask/problem regions once and filter them from ARG-based
# datasets before downstream thresholding, model assignment, heat maps, and
# Manhattan plots. Intervals are treated as BED-style [start, end).
read_argweaver_mask <- function(argweaver_mask_file) {
  # ARGweaver_mask.bed is expected to be a BED file with either:
  #   chrom  chromStart  chromEnd
  # or no header and at least three columns.
  #
  # BED convention is [chromStart, chromEnd), so chromEnd is converted to
  # chromEnd - 1 for foverlaps(), which expects inclusive interval ends.
  mask_raw <- fread(argweaver_mask_file)

  if (all(c("chrom", "chromStart", "chromEnd") %in% names(mask_raw))) {
    mask_regions <- mask_raw[, .(
      chromosome = chrom,
      mask_start = as.integer(chromStart),
      mask_end = as.integer(chromEnd)
    )]
  } else {
    mask_raw <- fread(
      argweaver_mask_file,
      header = FALSE,
      select = 1:3,
      col.names = c("chromosome", "mask_start", "mask_end")
    )

    # If the file actually had a header but not the canonical BED names,
    # drop that header-like first row after reading without a header.
    mask_regions <- mask_raw[
      suppressWarnings(!is.na(as.integer(mask_start))) &
        suppressWarnings(!is.na(as.integer(mask_end)))
    ]
    mask_regions[, `:=`(
      mask_start = as.integer(mask_start),
      mask_end = as.integer(mask_end)
    )]
  }

  mask_regions[, chromosome := standardize_scaffold_names(chromosome)]
  mask_regions[, mask_end_inclusive := as.integer(mask_end - 1L)]
  mask_regions <- mask_regions[mask_end_inclusive >= mask_start]

  setkey(mask_regions, chromosome, mask_start, mask_end_inclusive)
  mask_regions[]
}

filter_points_outside_mask <- function(dt, mask_regions, pos_col = "pos") {
  if (nrow(dt) == 0 || nrow(mask_regions) == 0) {
    return(dt)
  }

  tmp <- copy(dt)
  tmp[, repeat_filter_row_id := .I]
  tmp[, point_start := as.integer(get(pos_col))]
  tmp[, point_end := as.integer(get(pos_col))]

  setkey(tmp, chromosome, point_start, point_end)

  overlaps <- foverlaps(
    tmp[, .(repeat_filter_row_id, chromosome, point_start, point_end)],
    mask_regions,
    by.x = c("chromosome", "point_start", "point_end"),
    by.y = c("chromosome", "mask_start", "mask_end_inclusive"),
    nomatch = 0
  )

  repeat_rows <- unique(overlaps$repeat_filter_row_id)
  tmp <- tmp[!repeat_filter_row_id %in% repeat_rows]
  tmp[, c("repeat_filter_row_id", "point_start", "point_end") := NULL]
  tmp[]
}

filter_intervals_outside_mask <- function(
  dt,
  mask_regions,
  start_col = "interval_start",
  end_col = "interval_end"
) {
  if (nrow(dt) == 0 || nrow(mask_regions) == 0) {
    return(dt)
  }

  tmp <- copy(dt)
  tmp[, repeat_filter_row_id := .I]
  tmp[, interval_start_for_repeat_filter := as.integer(get(start_col))]
  tmp[, interval_end_inclusive_for_repeat_filter := as.integer(get(end_col) - 1L)]
  tmp <- tmp[interval_end_inclusive_for_repeat_filter >= interval_start_for_repeat_filter]

  setkey(tmp, chromosome, interval_start_for_repeat_filter, interval_end_inclusive_for_repeat_filter)

  overlaps <- foverlaps(
    tmp[, .(
      repeat_filter_row_id,
      chromosome,
      interval_start_for_repeat_filter,
      interval_end_inclusive_for_repeat_filter
    )],
    mask_regions,
    by.x = c(
      "chromosome",
      "interval_start_for_repeat_filter",
      "interval_end_inclusive_for_repeat_filter"
    ),
    by.y = c("chromosome", "mask_start", "mask_end_inclusive"),
    nomatch = 0
  )

  repeat_rows <- unique(overlaps$repeat_filter_row_id)
  tmp <- tmp[!repeat_filter_row_id %in% repeat_rows]
  tmp[, c(
    "repeat_filter_row_id",
    "interval_start_for_repeat_filter",
    "interval_end_inclusive_for_repeat_filter"
  ) := NULL]
  tmp[]
}

sanitize_path_component <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

model_label_without_redundant_model <- function(x) {
  gsub(" model$", "", x, ignore.case = FALSE)
}

get_threshold_for_manhattan <- function(stat, thresholds_dt) {
  out <- thresholds_dt[statistic == stat]

  if (nrow(out) == 0) {
    return(data.table(
      stat = stat,
      tail = NA_character_,
      threshold = NA_real_
    ))
  }

  data.table(
    stat = stat,
    tail = out$tail,
    threshold = out$threshold
  )
}

make_model_justification <- function(model_row) {
  sigs <- character()

  if (isTRUE(model_row$significant_Tap_Xin_RTH_original)) {
    sigs <- c(sigs, "Tapajos-Xingu RTH")
  }
  if (isTRUE(model_row$significant_Xin_Bel_CC_original)) {
    sigs <- c(sigs, "Xingu-Belem RCC")
  }
  if (isTRUE(model_row$significant_xingu_enrich)) {
    sigs <- c(sigs, "Xingu enrichment")
  }
  if (isTRUE(model_row$significant_belem_enrich)) {
    sigs <- c(sigs, "Belem enrichment")
  }
  if (isTRUE(model_row$significant_xingu_RTH)) {
    sigs <- c(sigs, "Xingu RTH'")
  }
  if (isTRUE(model_row$significant_belem_RTH)) {
    sigs <- c(sigs, "Belem RTH'")
  }

  if (length(sigs) == 0) {
    return("no modeled ARG statistic significant")
  }

  paste(sigs, collapse = "; ")
}

make_supported_model_string <- function(model_row) {
  models <- character()

  if (isTRUE(model_row$introgression_model)) {
    models <- c(models, "Introgression")
  }
  if (isTRUE(model_row$deep_lineage_sorting_model)) {
    models <- c(models, "Deep lineage sorting")
  }
  if (isTRUE(model_row$selection_bottleneck_model)) {
    models <- c(models, "Selection-bottleneck")
  }
  if (isTRUE(model_row$selection_recombination_model)) {
    models <- c(models, "Selection-recombination")
  }

  paste(models, collapse = "; ")
}

make_overlapping_signal_string <- function(model_row) {
  make_model_justification(model_row)
}

native_values_in_window <- function(native_dt, chromosome_i, start_i, end_i, stats = NULL) {
  tmp <- native_dt[
    chromosome == chromosome_i &
      x >= start_i &
      x < end_i &
      !is.na(value)
  ]

  if (!is.null(stats)) {
    tmp <- tmp[stat %in% stats]
  }

  nrow(tmp) > 0
}

choose_nonfocal_window <- function(
  candidates,
  target_pos,
  all_peaks,
  focal_peak_id,
  native_dt,
  arg_data_stats,
  exclude_window_midpoints = numeric(0)
) {
  if (nrow(candidates) == 0) {
    return(data.table())
  }

  tmp <- copy(candidates)
  tmp[, target_distance := abs(window_midpoint - target_pos)]
  setorder(tmp, target_distance)

  for (i in seq_len(nrow(tmp))) {
    pos <- tmp$window_midpoint[i]

    # Avoid choosing the focal window or reusing the same non-focal window.
    if (length(exclude_window_midpoints) > 0 &&
        any(abs(pos - exclude_window_midpoints) < 1e-6, na.rm = TRUE)) {
      next
    }

    overlaps_other_peak <- nrow(all_peaks[
      chromosome == tmp$chromosome[i] &
        start <= pos &
        end >= pos
    ]) > 0

    has_arg_data <- native_values_in_window(
      native_dt = native_dt,
      chromosome_i = tmp$chromosome[i],
      start_i = tmp$window_start[i],
      end_i = tmp$window_end[i],
      stats = arg_data_stats
    )

    if (!overlaps_other_peak && has_arg_data) {
      return(tmp[i])
    }
  }

  data.table()
}

plot_manhattan_footnote_pdf_panel <- function(footnote_pdf_file) {
  old_mar <- par(mar = c(0, 0, 0, 0), xpd = NA)
  on.exit(par(old_mar), add = TRUE)

  if (is.null(footnote_pdf_file) ||
      is.na(footnote_pdf_file) ||
      !nzchar(footnote_pdf_file) ||
      !file.exists(footnote_pdf_file)) {
    plot.new()
    text(
      0.5,
      0.5,
      "Footnote PDF not found: manhattan_footnotes.pdf",
      cex = 0.8
    )
    return(invisible(NULL))
  }

  footnote_raster <- NULL

  if (requireNamespace("pdftools", quietly = TRUE)) {
    footnote_raster <- tryCatch(
      as.raster(pdftools::pdf_render_page(
        pdf = footnote_pdf_file,
        page = 1,
        dpi = 200
      )),
      error = function(e) NULL
    )
  }

  if (is.null(footnote_raster) && requireNamespace("magick", quietly = TRUE)) {
    footnote_raster <- tryCatch(
      as.raster(magick::image_read_pdf(
        path = footnote_pdf_file,
        pages = 1,
        density = 200
      )),
      error = function(e) NULL
    )
  }

  if (is.null(footnote_raster)) {
    plot.new()
    text(
      0.5,
      0.5,
      paste0(
        "Could not render manhattan_footnotes.pdf.\n",
        "Install either the R package 'pdftools' or 'magick'."
      ),
      cex = 0.8
    )
    return(invisible(NULL))
  }

  # Important: do not call plot.new() before this plot() call. In a layout(),
  # each new plot advances to the next layout cell; calling plot.new() here and
  # then plot() again would consume the footer cell and push the footnote PDF
  # onto a second page.
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    bty = "n"
  )

  # Preserve the original width:height proportions of manhattan_footnotes.pdf
  # when fitting it into the bottom footer panel.
  raster_dim <- dim(footnote_raster)
  raster_height_px <- raster_dim[1]
  raster_width_px <- raster_dim[2]
  raster_aspect <- raster_width_px / raster_height_px

  panel_inches <- par("pin")
  panel_aspect <- panel_inches[1] / panel_inches[2]

  if (is.finite(raster_aspect) && is.finite(panel_aspect) &&
      raster_aspect > 0 && panel_aspect > 0) {
    if (raster_aspect >= panel_aspect) {
      # PDF is wider than the panel: use full panel width and center vertically.
      image_width <- 1
      image_height <- panel_aspect / raster_aspect
      xleft <- 0
      xright <- 1
      ybottom <- (1 - image_height) / 2
      ytop <- ybottom + image_height
    } else {
      # PDF is taller than the panel: use full panel height and center horizontally.
      image_height <- 1
      image_width <- raster_aspect / panel_aspect
      xleft <- (1 - image_width) / 2
      xright <- xleft + image_width
      ybottom <- 0
      ytop <- 1
    }
  } else {
    xleft <- 0
    xright <- 1
    ybottom <- 0
    ytop <- 1
  }

  rasterImage(
    footnote_raster,
    xleft = xleft,
    ybottom = ybottom,
    xright = xright,
    ytop = ytop,
    interpolate = TRUE
  )

  invisible(NULL)
}


make_manhattan_panel_pdf <- function(
  panel_dt,
  peak_row,
  model_row,
  thresholds_dt,
  manhattan_stats,
  plot_labels,
  focal_window_midpoint,
  nonfocal_window_1_midpoint,
  nonfocal_window_2_midpoint,
  outfile,
  tree_plot_dt = NULL,
  footnote_pdf_file = NULL,
  multi_iter_panel_dt = NULL,
  multi_iter_stats = character(),
  multi_iter_colors = NULL,
  repeat_regions = data.table()
) {
  n_stats <- length(manhattan_stats)

  if (is.null(multi_iter_colors)) {
    multi_iter_colors <- c(
      "1600" = "#FE6100",
      "1700" = "#FFB000",
      "1800" = "#648FFF",
      "1900" = "#785EF0",
      "2000" = "#DC267F"
    )
  }

  pdf(
    outfile,
    width = 15.2,
    height = max(14.8, 1.02 * n_stats + 7.6),
    useDingbats = FALSE,
    compress = TRUE
  )
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  })

  # Equal layout heights plus shared outer margins prevents the first and last
  # panels from being squeezed by the title and x-axis labels.
  # The right column is reserved for four sampled ARG local trees.
  right_ids <- rep(seq.int(n_stats + 1L, n_stats + 4L), length.out = n_stats)
  right_ids <- right_ids[order(rep(seq_len(4L), length.out = n_stats))]
  tree_breaks <- floor(seq(0, n_stats, length.out = 5L))
  for (tree_i in seq_len(4L)) {
    rows_i <- seq.int(tree_breaks[tree_i] + 1L, tree_breaks[tree_i + 1L])
    right_ids[rows_i] <- n_stats + tree_i
  }

  spacer_layout_id <- n_stats + 5L
  footnote_layout_id <- n_stats + 6L

  layout(
    rbind(
      cbind(seq_len(n_stats), right_ids),
      c(spacer_layout_id, spacer_layout_id),
      c(footnote_layout_id, footnote_layout_id)
    ),
    widths = c(4.05, 1.95),
    heights = c(rep(1, n_stats), 0.65, 3.4)
  )
  # Footnotes are drawn from manhattan_footnotes.pdf as a dedicated bottom
  # panel. The spacer row separates the footer from the lowest Manhattan
  # subplot so the x-axis title is not clipped or crowded.
  par(oma = c(1.1, 0, 7.4, 0))

  make_model_signal_text <- function(model_name, signals) {
    paste0(model_name, " (due to ", paste(signals, collapse = "; "), ")")
  }

  introgression_signals <- if (isTRUE(model_row$introgression_model)) {
    "Tapajos-Xingu RTH"
  } else {
    character()
  }

  deep_lineage_sorting_signals <- if (isTRUE(model_row$deep_lineage_sorting_model)) {
    c("Xingu-Belem RCC", "Xingu enrichment", "Belem enrichment")
  } else {
    character()
  }

  selection_bottleneck_signals <- character()
  if (xor(isTRUE(model_row$significant_xingu_enrich), isTRUE(model_row$significant_belem_enrich))) {
    if (isTRUE(model_row$significant_xingu_enrich)) {
      selection_bottleneck_signals <- c(selection_bottleneck_signals, "Xingu enrichment")
    }
    if (isTRUE(model_row$significant_belem_enrich)) {
      selection_bottleneck_signals <- c(selection_bottleneck_signals, "Belem enrichment")
    }
  }
  if (xor(isTRUE(model_row$significant_xingu_RTH), isTRUE(model_row$significant_belem_RTH))) {
    if (isTRUE(model_row$significant_xingu_RTH)) {
      selection_bottleneck_signals <- c(selection_bottleneck_signals, "Xingu RTH'")
    }
    if (isTRUE(model_row$significant_belem_RTH)) {
      selection_bottleneck_signals <- c(selection_bottleneck_signals, "Belem RTH'")
    }
  }

  selection_recombination_signals <- character()
  if (!isTRUE(model_row$significant_Xin_Bel_CC_original) &
      isTRUE(model_row$significant_xingu_enrich) &
      isTRUE(model_row$significant_belem_enrich)) {
    selection_recombination_signals <- c(
      selection_recombination_signals,
      "Xingu enrichment",
      "Belem enrichment"
    )
  }
  if (isTRUE(model_row$significant_xingu_RTH) & isTRUE(model_row$significant_belem_RTH)) {
    selection_recombination_signals <- c(
      selection_recombination_signals,
      "Xingu RTH'",
      "Belem RTH'"
    )
  }

  supported_model_text <- character()
  if (isTRUE(model_row$introgression_model)) {
    supported_model_text <- c(
      supported_model_text,
      make_model_signal_text("Introgression", introgression_signals)
    )
  }
  if (isTRUE(model_row$deep_lineage_sorting_model)) {
    supported_model_text <- c(
      supported_model_text,
      make_model_signal_text("Deep lineage sorting", deep_lineage_sorting_signals)
    )
  }
  if (isTRUE(model_row$selection_bottleneck_model)) {
    supported_model_text <- c(
      supported_model_text,
      make_model_signal_text("Selection-bottleneck", selection_bottleneck_signals)
    )
  }
  if (isTRUE(model_row$selection_recombination_model)) {
    supported_model_text <- c(
      supported_model_text,
      make_model_signal_text("Selection-recombination", selection_recombination_signals)
    )
  }

  if ("rcc_plus_selection_bottleneck_overlap" %in% names(model_row) &&
      isTRUE(model_row$rcc_plus_selection_bottleneck_overlap)) {
    supported_model_text <- c(
      supported_model_text,
      make_model_signal_text("Xingu-Belem RCC signal", "Xingu-Belem RCC")
    )
  }

  if (length(supported_model_text) == 0) {
    supported_model_text <- "Unassigned (due to no modeled ARG statistic significant)"
  }

  assignment_label <- model_label_without_redundant_model(model_row$model_assignment)

  unassigned_reason <- if ("assignment_threshold_type" %in% names(model_row) &&
                             identical(model_row$assignment_threshold_type, "25perc_trees_signif_in_peak")) {
    "Model: Unassigned (due to no ARG statistic with >=25% of local trees significant within any Fst outlier window in the focal peak)"
  } else {
    "Model: Unassigned (due to no ARG statistic significant within Fst outlier windows of focal peak)"
  }

  model_text <- if (identical(model_row$model_assignment, "Overlapping models")) {
    paste0("Overlapping models: ", paste(supported_model_text, collapse = "; "))
  } else if (identical(model_row$model_assignment, "Unassigned to model")) {
    unassigned_reason
  } else {
    paste0("Model: ", assignment_label, " (due to ", make_model_justification(model_row), ")")
  }

  peak_text <- paste0(
    "Focal Xingu-Belem peak: ",
    peak_row$chromosome,
    ":",
    peak_row$start,
    "-",
    peak_row$end
  )

  # Footnotes are wrapped without hanging indents so "Footnotes:",
  # "1.", "2.", "3.", and all wrapped continuation lines share one fixed
  # left edge. This avoids inconsistent alignment caused by leading spaces in
  # mtext() with proportional PDF fonts.
  footnote_text_1 <- paste0(
    "1. We defined the focal Xingu-Belem Fst window as the window with the highest Fst value in the focal peak. ",
    "This window may not be centered in the Manhattan plot due to the peak occurring close to the beginning or end of a scaffold."
  )

  footnote_text_2 <- paste0(
    "2. All empirical significance thresholds in the Manhattan plots (shown as horizontal dashed lines) indicate upper-tail ",
    "thresholds (values above line are significant), except for Belem pi, Xingu pi, and recombination rate which show ",
    "lower-tail thresholds (values below line are significant)."
  )

  footnote_text_3 <- paste0(
    "3. Empirical significance thresholds for the traditional and ARG-based Xingu-Belem Fst are set at five standard ",
    "deviations above the mean. Empirical significance thresholds for all remaining statistics, except recombination rate, ",
    "are set at the 99.9th or 0.1th percentile based on value distributions across control window regions. The empirical ",
    "significance threshold for recombination rate is set at the 0.5th percentile because the 0.1th percentile results in ",
    "a threshold of zero (due to a floor effect caused by the presence of some recombination rate values of zero)."
  )

  footnote_wrapped_lines <- c(
    "Footnotes:",
    strwrap(footnote_text_1, width = 180),
    strwrap(footnote_text_2, width = 180),
    strwrap(footnote_text_3, width = 180)
  )

  # Peak/model text is added after the first subplot is initialized.
  # Calling mtext(..., outer = TRUE) before plot.new() can trigger:
  #   "plot.new has not been called yet"

  vertical_positions <- c(
    focal_window_midpoint,
    nonfocal_window_1_midpoint,
    nonfocal_window_2_midpoint
  )

  vertical_labels <- c(
    "focal window",
    "non-focal window 1",
    "non-focal window 2"
  )

  for (i in seq_along(manhattan_stats)) {
    current_stat <- manhattan_stats[i]
    tmp <- panel_dt[panel_dt[["stat"]] == current_stat & !is.na(value)]
    setorder(tmp, x)

    par(
      mar = if (i == 1) {
        c(0.85, 5.2, 0.85, 1.4)
      } else {
        c(0.85, 5.2, 0.85, 1.4)
      }
    )

    # Strictly subset to the threshold row(s) for the statistic currently
    # being plotted. Use base indexing to avoid data.table name-shadowing
    # between the column called "stat" and the loop variable current_stat.
    threshold_rows <- unique(
      thresholds_dt[
        which(thresholds_dt[["stat"]] == current_stat &
                !is.na(thresholds_dt[["threshold"]]))
      ][, .(stat, tail, threshold)]
    )
    threshold_values <- threshold_rows$threshold

    y_all <- c(tmp$value, threshold_values)
    if (current_stat %in% multi_iter_stats &&
        !is.null(multi_iter_panel_dt) &&
        nrow(multi_iter_panel_dt) > 0) {
      y_all <- c(
        y_all,
        multi_iter_panel_dt[
          stat == current_stat,
          value
        ]
      )
    }
    if (length(y_all) == 0 || all(is.na(y_all))) {
      y_all <- c(0, 1)
    }

    y_min <- min(y_all, na.rm = TRUE)
    y_max <- max(y_all, na.rm = TRUE)
    if (!is.finite(y_min) || !is.finite(y_max) || y_min == y_max) {
      y_min <- y_min - 0.5
      y_max <- y_max + 0.5
    }
    pad <- 0.08 * (y_max - y_min)
    y_limits <- c(y_min - pad, y_max + pad)

    x_limits_bp <- c(peak_row$panel_start, peak_row$panel_end)

    y_label <- wrap_axis_label(unname(plot_labels[current_stat]), width = 20L)

    plot(
      NA,
      xlim = x_limits_bp,
      ylim = y_limits,
      xaxt = "n",
      xlab = "",
      ylab = y_label,
      cex.lab = 0.74,
      cex.axis = 0.622,
      main = ""
    )

    draw_repeat_region_bars(
      repeat_regions = repeat_regions,
      chromosome_i = peak_row$chromosome,
      x_limits_bp = x_limits_bp,
      y_limits = y_limits
    )

    if (i == 1) {
      mtext(peak_text, side = 3, outer = TRUE, line = 5.8, cex = 0.86, font = 1)
      mtext(
        paste(strwrap(model_text, width = 135), collapse = "\n"),
        side = 3,
        outer = TRUE,
        line = 4.2,
        cex = 0.86,
        font = 1
      )
    }

    if (current_stat %in% multi_iter_stats &&
        !is.null(multi_iter_panel_dt) &&
        nrow(multi_iter_panel_dt) > 0) {
      multi_tmp <- multi_iter_panel_dt[
        stat == current_stat &
          !is.na(value)
      ]
      if (nrow(multi_tmp) > 0) {
        for (iter_i in sort(unique(multi_tmp$iteration))) {
          iter_tmp <- multi_tmp[iteration == iter_i]
          setorder(iter_tmp, x)
          iter_tmp[, rolling_mean := rolling_mean_centered(value, half_window = 5L)]

          iter_col <- multi_iter_colors[as.character(iter_i)]
          if (is.na(iter_col)) iter_col <- "red"

          lines(
            iter_tmp$x,
            iter_tmp$rolling_mean,
            col = iter_col,
            lwd = 1.15
          )
        }

        legend(
          "topleft",
          legend = names(multi_iter_colors),
          col = unname(multi_iter_colors),
          lwd = 1.4,
          bty = "n",
          cex = 0.58,
          title = "Iteration"
        )
      }
    } else {
      # Plot all native-resolution points as vector black circles.
      # This intentionally avoids rasterization because the raster point layer
      # produced visual artifacts in some Manhattan panels.
      if (nrow(tmp) > 0) {
        points(
          tmp$x,
          tmp$value,
          pch = 16,
          cex = 0.75,
          col = "black"
        )
      }

      if (nrow(tmp) > 0) {
        tmp[, rolling_mean := rolling_mean_centered(value, half_window = 5L)]
        lines(
          tmp$x,
          tmp$rolling_mean,
          col = "red",
          lwd = 1.0
        )
      }
    }

    # Only draw the threshold(s) for the statistic currently being plotted.
    # iHS/nSL have two statistic-specific thresholds (-2 and 2); most other
    # statistics have one statistic-specific empirical threshold.
    if (nrow(threshold_rows) > 0) {
      for (thr_i in seq_len(nrow(threshold_rows))) {
        if (!is.na(threshold_rows$threshold[thr_i])) {
          abline(h = threshold_rows$threshold[thr_i], col = "grey35", lty = 2, lwd = 1)
        }
      }
    }

    for (v in vertical_positions) {
      if (!is.na(v)) {
        abline(v = v, lty = 2, lwd = 1, col = "grey35")
      }
    }

    if (i == 1) {
      usr <- par("usr")
      label_y <- usr[4] + 0.09 * diff(usr[3:4])
      old_xpd <- par(xpd = NA)

      # Draw these labels after the rasterized point layer, rolling mean,
      # thresholds, and vertical guide lines so they remain visible and are
      # not clipped by the panel plotting region.
      for (lab_i in seq_along(vertical_positions)) {
        if (!is.na(vertical_positions[lab_i])) {
          text(
            x = vertical_positions[lab_i],
            y = label_y,
            labels = vertical_labels[lab_i],
            srt = 45,
            adj = c(0, 0),
            cex = 0.62,
            col = "black"
          )
        }
      }

      par(xpd = old_xpd)
    }

    if (i == n_stats) {
      axis_tick_step_bp <- 250000L
      axis_ticks_bp <- seq(
        ceiling(x_limits_bp[1] / axis_tick_step_bp) * axis_tick_step_bp,
        floor(x_limits_bp[2] / axis_tick_step_bp) * axis_tick_step_bp,
        by = axis_tick_step_bp
      )

      axis(
        1,
        at = axis_ticks_bp,
        labels = format(round(axis_ticks_bp / 1e6, 2), trim = TRUE),
        cex.axis = 0.622
      )

      mtext(
        format_scaffold_axis_label(peak_row$chromosome),
        side = 1,
        line = 3.0,
        cex = 0.82
      )
    }

    box()
  }
  # Plot sampled local trees in the right-side panels.
  for (tree_i in seq_len(4L)) {
    par(mar = c(1.0, 0.8, 3.4, 0.8), xpd = NA)
    if (!is.null(tree_plot_dt) && nrow(tree_plot_dt) >= tree_i) {
      plot_single_tree_panel(tree_plot_dt[tree_i])
    } else {
      plot.new()
      text(0.5, 0.5, "No tree sampled", cex = 0.7)
    }
  }

  # Add a small blank spacer row before the externally prepared footnote PDF.
  # This keeps the footer from clipping the x-axis labels/title of the lowest
  # Manhattan subplot.
  par(mar = c(0, 0, 0, 0))
  plot.new()

  # Add the externally prepared footnote PDF as a dedicated bottom
  # panel. This preserves the footnote formatting exactly as prepared in
  # manhattan_footnotes.pdf and avoids fragile mtext() alignment issues.
  plot_manhattan_footnote_pdf_panel(footnote_pdf_file)

}

create_manhattan_panels_for_assignments <- function(
  model_dt,
  assignment_dir,
  manhattan_native_values,
  manhattan_window_reference,
  thresholds_dt,
  peaks_dt,
  plot_labels,
  manhattan_stats,
  scaffold_lengths,
  peak_window_lookup_dt,
  tree_map = data.table(),
  tree_base_dir = NULL,
  tree_iteration = 2000L,
  footnote_pdf_file = NULL,
  panel_subdir = "Manhattan_panels_2000th_iter",
  multi_iter_native_values = NULL,
  multi_iter_stats = character(),
  multi_iter_colors = NULL,
  multi_iter_iterations = NULL,
  multi_iter_arg_stat_dir = NULL,
  multi_iter_arg_cc_stat_dir = NULL,
  multi_iter_argweaver_mask = NULL,
  repeat_regions = data.table()
) {
  if (nrow(model_dt) == 0) {
    return(data.table())
  }

  panel_base_dir <- file.path(assignment_dir, panel_subdir)
  dir.create(panel_base_dir, showWarnings = FALSE, recursive = TRUE)

  arg_data_stats <- c(
    "Xin_Bel_Fst_ARG_based",
    "xingu_enrich",
    "belem_enrich",
    "xingu_RTH_inverse",
    "belem_RTH_inverse",
    "Xin_Bel_CC_original",
    "Tap_Xin_RTH_original_inverse"
  )

  panel_summaries <- list()

  for (row_i in seq_len(nrow(model_dt))) {
    model_row <- model_dt[row_i]
    peak_id_i <- model_row$peak_id
    chrom_i <- model_row$chromosome

    scaffold_max <- scaffold_lengths[
      chromosome == chrom_i,
      scaffold_max_end
    ]

    if (length(scaffold_max) == 0 || is.na(scaffold_max)) {
      next
    }

    peak_center <- floor((model_row$start + model_row$end) / 2)
    panel_start <- max(0L, as.integer(peak_center - 1000000L))
    panel_end <- panel_start + 2000000L

    if (panel_end > scaffold_max) {
      panel_end <- as.integer(scaffold_max)
      panel_start <- max(0L, as.integer(panel_end - 2000000L))
    }

    panel_dt <- manhattan_native_values[
      chromosome == chrom_i &
        x >= panel_start &
        x <= panel_end &
        stat %in% manhattan_stats
    ]

    if (nrow(panel_dt) == 0) {
      next
    }

    panel_windows <- manhattan_window_reference[
      chromosome == chrom_i &
        window_start >= panel_start &
        window_end <= panel_end
    ]

    if (nrow(panel_windows) == 0) {
      next
    }

    focal_candidates <- merge(
      peak_window_lookup_dt[peak_id == peak_id_i],
      panel_windows,
      by = c("chromosome", "window_start", "window_end"),
      all.x = FALSE,
      all.y = FALSE,
      sort = FALSE
    )
    focal_candidates <- focal_candidates[!is.na(Xin_Bel_Fst)]

    if (nrow(focal_candidates) == 0) {
      focal_candidates <- panel_windows[!is.na(Xin_Bel_Fst)]
    }

    focal_window <- focal_candidates[which.max(Xin_Bel_Fst)]
    focal_mid <- focal_window$window_midpoint[1]

    left_edge_close <- (model_row$start - panel_start) <= 500000
    right_edge_close <- (panel_end - model_row$end) <= 500000

    if (left_edge_close) {
      target_1 <- focal_mid + 500000
      target_2 <- focal_mid + 1000000
    } else if (right_edge_close) {
      target_1 <- focal_mid - 1000000
      target_2 <- focal_mid - 500000
    } else {
      target_1 <- focal_mid - 750000
      target_2 <- focal_mid + 750000
    }

    candidate_windows <- panel_windows[!is.na(Xin_Bel_Fst) & is_Xin_Bel_Fst_outlier == FALSE]

    nonfocal_1 <- choose_nonfocal_window(
      candidate_windows,
      target_1,
      peaks_dt,
      peak_id_i,
      native_dt = panel_dt,
      arg_data_stats = arg_data_stats,
      exclude_window_midpoints = focal_mid
    )

    nonfocal_2 <- choose_nonfocal_window(
      candidate_windows,
      target_2,
      peaks_dt,
      peak_id_i,
      native_dt = panel_dt,
      arg_data_stats = arg_data_stats,
      exclude_window_midpoints = c(
        focal_mid,
        if (nrow(nonfocal_1) > 0) nonfocal_1$window_midpoint[1] else NA_real_
      )
    )

    nonfocal_positions <- sort(c(
      if (nrow(nonfocal_1) > 0) nonfocal_1$window_midpoint[1] else NA_real_,
      if (nrow(nonfocal_2) > 0) nonfocal_2$window_midpoint[1] else NA_real_
    ), na.last = TRUE)

    model_folder_name <- if (
      is.na(model_row$model_assignment) ||
        !nzchar(as.character(model_row$model_assignment))
    ) {
      "ARG_data_missing"
    } else {
      sanitize_path_component(model_row$model_assignment)
    }

    peak_model_dir <- file.path(
      panel_base_dir,
      model_folder_name
    )
    dir.create(peak_model_dir, showWarnings = FALSE, recursive = TRUE)

    outfile <- file.path(
      peak_model_dir,
      paste0(
        "peak_",
        peak_id_i,
        "_",
        chrom_i,
        "_",
        panel_start,
        "_",
        panel_end,
        ".pdf"
      )
    )

    peak_row_for_plot <- copy(model_row)
    peak_row_for_plot[, `:=`(
      panel_start = panel_start,
      panel_end = panel_end
    )]

    tree_plot_dt <- sample_manhattan_panel_trees(
      chromosome = chrom_i,
      focal_window_midpoint = focal_mid,
      nonfocal_window_1_midpoint = nonfocal_positions[1],
      nonfocal_window_2_midpoint = nonfocal_positions[2],
      manhattan_window_reference = manhattan_window_reference,
      tree_map = tree_map,
      tree_base_dir = tree_base_dir,
      iteration = tree_iteration
    )

    multi_iter_panel_dt <- NULL
    if (length(multi_iter_stats) > 0) {
      if (!is.null(multi_iter_native_values) &&
          nrow(multi_iter_native_values) > 0) {
        # Backward-compatible path: use a precomputed multi-iteration object
        # if one was explicitly provided.
        multi_iter_panel_dt <- multi_iter_native_values[
          chromosome == chrom_i &
            x >= panel_start &
            x <= panel_end &
            stat %in% multi_iter_stats
        ]
      } else if (!is.null(multi_iter_iterations) &&
                 !is.null(multi_iter_arg_stat_dir) &&
                 !is.null(multi_iter_arg_cc_stat_dir) &&
                 !is.null(multi_iter_argweaver_mask)) {
        # Memory-efficient path: read and retain only values inside the
        # current 2-Mb Manhattan panel. This avoids storing all iterations
        # genome-wide at once, which can exceed the Mac R vector memory limit.
        multi_iter_panel_dt <- read_multi_iter_values_for_panel(
          chromosome_i = chrom_i,
          panel_start = panel_start,
          panel_end = panel_end,
          iterations = multi_iter_iterations,
          arg_stat_dir = multi_iter_arg_stat_dir,
          arg_cc_stat_dir = multi_iter_arg_cc_stat_dir,
          argweaver_mask = multi_iter_argweaver_mask
        )
      }
    }

    make_manhattan_panel_pdf(
      panel_dt = panel_dt,
      peak_row = peak_row_for_plot,
      model_row = model_row,
      thresholds_dt = thresholds_dt,
      manhattan_stats = manhattan_stats,
      plot_labels = plot_labels,
      focal_window_midpoint = focal_mid,
      nonfocal_window_1_midpoint = nonfocal_positions[1],
      nonfocal_window_2_midpoint = nonfocal_positions[2],
      outfile = outfile,
      tree_plot_dt = tree_plot_dt,
      footnote_pdf_file = footnote_pdf_file,
      multi_iter_panel_dt = multi_iter_panel_dt,
      multi_iter_stats = multi_iter_stats,
      multi_iter_colors = multi_iter_colors,
      repeat_regions = repeat_regions
    )

    panel_summaries[[length(panel_summaries) + 1L]] <- data.table(
      peak_id = peak_id_i,
      chromosome = chrom_i,
      peak_start = model_row$start,
      peak_end = model_row$end,
      panel_start = panel_start,
      panel_end = panel_end,
      model_assignment = model_row$model_assignment,
      n_supported_models = model_row$n_supported_models,
      focal_window_midpoint = focal_mid,
      nonfocal_window_1_midpoint = nonfocal_positions[1],
      nonfocal_window_2_midpoint = nonfocal_positions[2],
      n_trees_sampled_for_plot = nrow(tree_plot_dt),
      sampled_tree_positions = if (nrow(tree_plot_dt) > 0) {
        paste(tree_plot_dt$tree_position, collapse = ";")
      } else {
        NA_character_
      },
      sampled_tree_files = if (nrow(tree_plot_dt) > 0) {
        paste(unique(tree_plot_dt$tree_file), collapse = ";")
      } else {
        NA_character_
      },
      output_pdf = outfile
    )
  }

  rbindlist(panel_summaries, fill = TRUE)
}


read_arg_iteration_for_windows <- function(iteration, arg_stat_dir, argweaver_mask, comparison_windows) {
  arg_file <- file.path(
    arg_stat_dir,
    paste0("argStats_midpoint.", iteration, ".stat.gz")
  )

  if (!file.exists(arg_file)) {
    warning("ARG stat file not found for iteration ", iteration, ": ", arg_file)
    return(data.table())
  }

  arg_iter <- fread(
    arg_file,
    select = c(
      "chrom",
      "pos",
      "Tap_Xin_RTH_original",
      "belem_RTH",
      "xingu_RTH",
      "belem_enrich",
      "xingu_enrich"
    )
  )

  arg_iter[, chromosome := standardize_scaffold_names(chrom)]
  arg_iter[, chrom := NULL]
  arg_iter <- filter_points_outside_mask(arg_iter, argweaver_mask, pos_col = "pos")

  arg_iter[, xingu_RTH_inverse := safe_inverse(xingu_RTH)]
  arg_iter[, belem_RTH_inverse := safe_inverse(belem_RTH)]
  arg_iter[, Tap_Xin_RTH_original_inverse := safe_inverse(Tap_Xin_RTH_original)]

  arg_iter[, pos_start := as.integer(pos)]
  arg_iter[, pos_end := as.integer(pos)]
  arg_iter[, pos := NULL]

  iter_windows <- comparison_windows[, .(
    chromosome,
    window_start = as.integer(start),
    window_end = as.integer(end),
    window_end_inclusive = as.integer(end - 1L),
    window_id,
    window_class
  )]

  setkey(arg_iter, chromosome, pos_start, pos_end)
  setkey(iter_windows, chromosome, window_start, window_end_inclusive)

  joined <- foverlaps(
    arg_iter,
    iter_windows,
    by.x = c("chromosome", "pos_start", "pos_end"),
    by.y = c("chromosome", "window_start", "window_end_inclusive"),
    nomatch = 0
  )

  joined[, pos := pos_start]
  joined[, c("pos_start", "pos_end", "window_end_inclusive") := NULL]
  joined[, iteration := as.integer(iteration)]

  joined[]
}

read_arg_cc_iteration_for_windows <- function(iteration, arg_cc_stat_dir, argweaver_mask, comparison_windows) {
  arg_cc_file <- file.path(
    arg_cc_stat_dir,
    paste0("argStats_midpoint_CC.", iteration, ".stat.gz")
  )

  if (!file.exists(arg_cc_file)) {
    warning("ARG CC stat file not found for iteration ", iteration, ": ", arg_cc_file)
    return(data.table())
  }

  arg_cc_iter <- fread(
    arg_cc_file,
    select = c("chrom", "pos", "Xin_Bel_CC_original")
  )

  arg_cc_iter[, chromosome := standardize_scaffold_names(chrom)]
  arg_cc_iter[, chrom := NULL]
  arg_cc_iter <- filter_points_outside_mask(arg_cc_iter, argweaver_mask, pos_col = "pos")

  arg_cc_iter[, pos_start := as.integer(pos)]
  arg_cc_iter[, pos_end := as.integer(pos)]
  arg_cc_iter[, pos := NULL]

  iter_windows <- comparison_windows[, .(
    chromosome,
    window_start = as.integer(start),
    window_end = as.integer(end),
    window_end_inclusive = as.integer(end - 1L),
    window_id,
    window_class
  )]

  setkey(arg_cc_iter, chromosome, pos_start, pos_end)
  setkey(iter_windows, chromosome, window_start, window_end_inclusive)

  joined <- foverlaps(
    arg_cc_iter,
    iter_windows,
    by.x = c("chromosome", "pos_start", "pos_end"),
    by.y = c("chromosome", "window_start", "window_end_inclusive"),
    nomatch = 0
  )

  joined[, pos := pos_start]
  joined[, c("pos_start", "pos_end", "window_end_inclusive") := NULL]
  joined[, iteration := as.integer(iteration)]

  joined[]
}


read_multi_iter_values_for_panel <- function(
  chromosome_i,
  panel_start,
  panel_end,
  iterations,
  arg_stat_dir,
  arg_cc_stat_dir,
  argweaver_mask
) {
  out <- list()

  for (iter in iterations) {
    arg_file <- file.path(
      arg_stat_dir,
      paste0("argStats_midpoint.", iter, ".stat.gz")
    )

    if (file.exists(arg_file)) {
      arg_iter <- fread(
        arg_file,
        select = c(
          "chrom",
          "pos",
          "Tap_Xin_RTH_original",
          "belem_RTH",
          "xingu_RTH",
          "belem_enrich",
          "xingu_enrich"
        )
      )

      arg_iter[, chromosome := standardize_scaffold_names(chrom)]
      arg_iter[, chrom := NULL]
      arg_iter[, pos := as.integer(pos)]

      # Keep only the current Manhattan panel before masking/reshaping.
      # This is the key memory-saving step.
      arg_iter <- arg_iter[
        chromosome == chromosome_i &
          pos >= as.integer(panel_start) &
          pos <= as.integer(panel_end)
      ]

      if (nrow(arg_iter) > 0) {
        arg_iter <- filter_points_outside_mask(arg_iter, argweaver_mask, pos_col = "pos")
        arg_iter[, xingu_RTH_inverse := safe_inverse(xingu_RTH)]
        arg_iter[, belem_RTH_inverse := safe_inverse(belem_RTH)]
        arg_iter[, Tap_Xin_RTH_original_inverse := safe_inverse(Tap_Xin_RTH_original)]
        arg_iter[, x := as.numeric(pos)]

        for (stat in c(
          "xingu_enrich",
          "belem_enrich",
          "xingu_RTH_inverse",
          "belem_RTH_inverse",
          "Tap_Xin_RTH_original_inverse"
        )) {
          out[[length(out) + 1L]] <- data.table(
            iteration = as.integer(iter),
            chromosome = arg_iter$chromosome,
            x = arg_iter$x,
            stat = stat,
            value = as.numeric(arg_iter[[stat]])
          )
        }
      }

      rm(arg_iter)
      gc()
    } else {
      warning("ARG stat file not found for multi-iteration Manhattan plotting: ", arg_file)
    }

    arg_cc_file <- file.path(
      arg_cc_stat_dir,
      paste0("argStats_midpoint_CC.", iter, ".stat.gz")
    )

    if (file.exists(arg_cc_file)) {
      arg_cc_iter <- fread(
        arg_cc_file,
        select = c("chrom", "pos", "Xin_Bel_CC_original")
      )

      arg_cc_iter[, chromosome := standardize_scaffold_names(chrom)]
      arg_cc_iter[, chrom := NULL]
      arg_cc_iter[, pos := as.integer(pos)]

      # Keep only the current Manhattan panel before masking/reshaping.
      arg_cc_iter <- arg_cc_iter[
        chromosome == chromosome_i &
          pos >= as.integer(panel_start) &
          pos <= as.integer(panel_end)
      ]

      if (nrow(arg_cc_iter) > 0) {
        arg_cc_iter <- filter_points_outside_mask(arg_cc_iter, argweaver_mask, pos_col = "pos")
        arg_cc_iter[, x := as.numeric(pos)]

        out[[length(out) + 1L]] <- data.table(
          iteration = as.integer(iter),
          chromosome = arg_cc_iter$chromosome,
          x = arg_cc_iter$x,
          stat = "Xin_Bel_CC_original",
          value = as.numeric(arg_cc_iter$Xin_Bel_CC_original)
        )
      }

      rm(arg_cc_iter)
      gc()
    } else {
      warning("ARG CC stat file not found for multi-iteration Manhattan plotting: ", arg_cc_file)
    }
  }

  if (length(out) == 0) {
    return(data.table(
      iteration = integer(),
      chromosome = character(),
      x = numeric(),
      stat = character(),
      value = numeric()
    ))
  }

  rbindlist(out, fill = TRUE)
}

make_multi_iter_native_values <- function(
  iterations,
  arg_stat_dir,
  arg_cc_stat_dir,
  argweaver_mask
) {
  out <- list()

  for (iter in iterations) {
    arg_file <- file.path(
      arg_stat_dir,
      paste0("argStats_midpoint.", iter, ".stat.gz")
    )

    if (file.exists(arg_file)) {
      arg_iter <- fread(
        arg_file,
        select = c(
          "chrom",
          "pos",
          "Tap_Xin_RTH_original",
          "belem_RTH",
          "xingu_RTH",
          "belem_enrich",
          "xingu_enrich"
        )
      )

      arg_iter[, chromosome := standardize_scaffold_names(chrom)]
      arg_iter[, chrom := NULL]
      arg_iter <- filter_points_outside_mask(arg_iter, argweaver_mask, pos_col = "pos")

      arg_iter[, xingu_RTH_inverse := safe_inverse(xingu_RTH)]
      arg_iter[, belem_RTH_inverse := safe_inverse(belem_RTH)]
      arg_iter[, Tap_Xin_RTH_original_inverse := safe_inverse(Tap_Xin_RTH_original)]
      arg_iter[, x := as.numeric(pos)]

      for (stat in c(
        "xingu_enrich",
        "belem_enrich",
        "xingu_RTH_inverse",
        "belem_RTH_inverse",
        "Tap_Xin_RTH_original_inverse"
      )) {
        out[[length(out) + 1L]] <- data.table(
          iteration = as.integer(iter),
          chromosome = arg_iter$chromosome,
          x = arg_iter$x,
          stat = stat,
          value = as.numeric(arg_iter[[stat]])
        )
      }

      rm(arg_iter)
      gc()
    } else {
      warning("ARG stat file not found for multi-iteration Manhattan plotting: ", arg_file)
    }

    arg_cc_file <- file.path(
      arg_cc_stat_dir,
      paste0("argStats_midpoint_CC.", iter, ".stat.gz")
    )

    if (file.exists(arg_cc_file)) {
      arg_cc_iter <- fread(
        arg_cc_file,
        select = c("chrom", "pos", "Xin_Bel_CC_original")
      )

      arg_cc_iter[, chromosome := standardize_scaffold_names(chrom)]
      arg_cc_iter[, chrom := NULL]
      arg_cc_iter <- filter_points_outside_mask(arg_cc_iter, argweaver_mask, pos_col = "pos")
      arg_cc_iter[, x := as.numeric(pos)]

      out[[length(out) + 1L]] <- data.table(
        iteration = as.integer(iter),
        chromosome = arg_cc_iter$chromosome,
        x = arg_cc_iter$x,
        stat = "Xin_Bel_CC_original",
        value = as.numeric(arg_cc_iter$Xin_Bel_CC_original)
      )

      rm(arg_cc_iter)
      gc()
    } else {
      warning("ARG CC stat file not found for multi-iteration Manhattan plotting: ", arg_cc_file)
    }
  }

  rbindlist(out, fill = TRUE)
}

summarize_multi_iter_ARG_significance <- function(
  iterations,
  barplot_stats,
  arg_stat_dir,
  arg_cc_stat_dir,
  argweaver_mask,
  comparison_windows,
  peak_window_lookup,
  window_reference,
  peak_reference,
  peak_coverage_dt,
  arg_threshold_summary,
  arg_cc_threshold_summary,
  plot_labels
) {
  window_counts_list <- list()
  peak_counts_list <- list()
  model_assignment_list <- list()
  model_summary_list <- list()

  for (iter in iterations) {
    cat("Processing multi-iteration ARG summaries for iteration ", iter, "\n", sep = "")

    arg_iter_joined <- read_arg_iteration_for_windows(
      iteration = iter,
      arg_stat_dir = arg_stat_dir,
      argweaver_mask = argweaver_mask,
      comparison_windows = comparison_windows
    )

    arg_cc_iter_joined <- read_arg_cc_iteration_for_windows(
      iteration = iter,
      arg_cc_stat_dir = arg_cc_stat_dir,
      argweaver_mask = argweaver_mask,
      comparison_windows = comparison_windows
    )

    arg_iter_joined_peaks <- if (nrow(arg_iter_joined) > 0) {
      merge(
        arg_iter_joined[window_class == "Xin_Bel_Fst_outlier"],
        peak_window_lookup,
        by = c("chromosome", "window_start", "window_end"),
        all.x = FALSE,
        all.y = FALSE
      )
    } else {
      data.table()
    }

    arg_cc_iter_joined_peaks <- if (nrow(arg_cc_iter_joined) > 0) {
      merge(
        arg_cc_iter_joined[window_class == "Xin_Bel_Fst_outlier"],
        peak_window_lookup,
        by = c("chromosome", "window_start", "window_end"),
        all.x = FALSE,
        all.y = FALSE
      )
    } else {
      data.table()
    }

    iter_window_status <- rbindlist(lapply(seq_len(nrow(barplot_stats)), function(i) {
      row <- barplot_stats[i]

      if (row$source == "ARG") {
        threshold_row <- arg_threshold_summary[
          statistic == row$statistic &
            tail == row$tail
        ]

        make_window_significance_status(
          dt = arg_iter_joined,
          stat = row$statistic,
          threshold = threshold_row$threshold[1],
          tail = row$tail,
          source_name = row$source,
          plot_label = row$plot_label,
          window_reference = window_reference,
          min_props = c(NA_real_, 0.25)
        )
      } else {
        threshold_row <- arg_cc_threshold_summary[
          statistic == row$statistic &
            tail == row$tail
        ]

        make_window_significance_status(
          dt = arg_cc_iter_joined,
          stat = row$statistic,
          threshold = threshold_row$threshold[1],
          tail = row$tail,
          source_name = row$source,
          plot_label = row$plot_label,
          window_reference = window_reference,
          min_props = c(NA_real_, 0.25)
        )
      }
    }), fill = TRUE)

    iter_peak_status <- rbindlist(lapply(seq_len(nrow(barplot_stats)), function(i) {
      row <- barplot_stats[i]

      if (row$source == "ARG") {
        threshold_row <- arg_threshold_summary[
          statistic == row$statistic &
            tail == row$tail
        ]

        make_peak_significance_status(
          dt = arg_iter_joined_peaks,
          stat = row$statistic,
          threshold = threshold_row$threshold[1],
          tail = row$tail,
          source_name = row$source,
          plot_label = row$plot_label,
          peak_reference = peak_reference,
          min_props = c(NA_real_, 0.25)
        )
      } else {
        threshold_row <- arg_cc_threshold_summary[
          statistic == row$statistic &
            tail == row$tail
        ]

        make_peak_significance_status(
          dt = arg_cc_iter_joined_peaks,
          stat = row$statistic,
          threshold = threshold_row$threshold[1],
          tail = row$tail,
          source_name = row$source,
          plot_label = row$plot_label,
          peak_reference = peak_reference,
          min_props = c(NA_real_, 0.25)
        )
      }
    }), fill = TRUE)

    iter_window_counts <- iter_window_status[, .(
      n_windows_with_data = sum(n_local_trees > 0),
      n_significant_windows = sum(is_significant == TRUE, na.rm = TRUE)
    ), by = .(
      statistic,
      source,
      plot_label,
      tail,
      threshold,
      threshold_type,
      min_prop_significant_local_trees
    )]

    iter_window_counts[, iteration := as.integer(iter)]

    iter_peak_counts <- iter_peak_status[, .(
      n_peaks_with_data = sum(n_local_trees > 0),
      n_significant_peaks = sum(is_significant == TRUE, na.rm = TRUE)
    ), by = .(
      statistic,
      source,
      plot_label,
      tail,
      threshold,
      threshold_type,
      min_prop_significant_local_trees
    )]

    iter_peak_counts[, iteration := as.integer(iter)]

    iter_model_any <- make_peak_model_assignments(
      iter_peak_status,
      peak_coverage_dt,
      min_prop = NA_real_
    )
    iter_model_any[, iteration := as.integer(iter)]

    iter_model_min25 <- make_peak_model_assignments(
      iter_peak_status,
      peak_coverage_dt,
      min_prop = 0.25
    )
    iter_model_min25[, iteration := as.integer(iter)]

    iter_models <- rbind(iter_model_any, iter_model_min25, fill = TRUE)

    iter_model_summary <- iter_models[, .(
      n_peaks = .N
    ), by = .(
      iteration,
      assignment_threshold_type,
      min_prop_significant_local_trees,
      model_assignment
    )]

    window_counts_list[[length(window_counts_list) + 1L]] <- iter_window_counts
    peak_counts_list[[length(peak_counts_list) + 1L]] <- iter_peak_counts
    model_assignment_list[[length(model_assignment_list) + 1L]] <- iter_models
    model_summary_list[[length(model_summary_list) + 1L]] <- iter_model_summary

    rm(
      arg_iter_joined,
      arg_cc_iter_joined,
      arg_iter_joined_peaks,
      arg_cc_iter_joined_peaks,
      iter_window_status,
      iter_peak_status,
      iter_window_counts,
      iter_peak_counts,
      iter_model_any,
      iter_model_min25,
      iter_models,
      iter_model_summary
    )
    gc()
  }

  list(
    window_counts = rbindlist(window_counts_list, fill = TRUE),
    peak_counts = rbindlist(peak_counts_list, fill = TRUE),
    model_assignments = rbindlist(model_assignment_list, fill = TRUE),
    model_summary = rbindlist(model_summary_list, fill = TRUE)
  )
}


make_multi_iter_model_assignment_wide <- function(
  multi_iter_model_dt,
  assignment_threshold_type_filter,
  iterations
) {
  model_categories <- c(
    "Introgression model",
    "Deep lineage sorting model",
    "Selection-bottleneck model",
    "Selection-recombination model",
    "Overlapping models",
    "Unassigned to model"
  )

  percentage_cols <- c(
    "percent_iterations_Introgression_model",
    "percent_iterations_Deep_lineage_sorting_model",
    "percent_iterations_Selection_bottleneck_model",
    "percent_iterations_Selection_recombination_model",
    "percent_iterations_Overlapping_models",
    "percent_iterations_Unassigned_to_model"
  )

  if (nrow(multi_iter_model_dt) == 0) {
    return(data.table())
  }

  tmp <- multi_iter_model_dt[
    assignment_threshold_type == assignment_threshold_type_filter &
      iteration %in% iterations,
    .(peak_id, iteration, model_assignment)
  ]

  if (nrow(tmp) == 0) {
    return(data.table())
  }

  # Protect against accidental duplicate rows for a peak/iteration.
  tmp <- unique(tmp, by = c("peak_id", "iteration"))
  tmp[, iteration_model_col := paste0("model_assignment_iter", iteration)]

  wide <- dcast(
    tmp,
    peak_id ~ iteration_model_col,
    value.var = "model_assignment"
  )

  expected_cols <- paste0("model_assignment_iter", iterations)
  for (col_i in expected_cols) {
    if (!col_i %in% names(wide)) {
      wide[, (col_i) := NA_character_]
    }
  }
  setcolorder(wide, c("peak_id", expected_cols))

  # Percentages are calculated across the full requested set of MCMC
  # iterations. Peaks with no model assignment in any iteration retain NA
  # percentages and an NA top model.
  assignment_matrix <- as.matrix(wide[, ..expected_cols])
  n_requested_iterations <- length(iterations)
  n_nonmissing <- rowSums(!is.na(assignment_matrix))

  for (i in seq_along(model_categories)) {
    category_i <- model_categories[i]
    col_i <- percentage_cols[i]
    counts_i <- rowSums(assignment_matrix == category_i, na.rm = TRUE)
    values_i <- 100 * counts_i / n_requested_iterations
    values_i[n_nonmissing == 0L] <- NA_real_
    wide[, (col_i) := values_i]
  }

  percent_matrix <- as.matrix(wide[, ..percentage_cols])
  wide[, top_model := vapply(seq_len(.N), function(i) {
    vals <- percent_matrix[i, ]
    if (all(is.na(vals))) return(NA_character_)
    max_val <- max(vals, na.rm = TRUE)
    # Ties are reported explicitly rather than resolved arbitrarily.
    winners <- model_categories[which(vals == max_val)]
    paste(winners, collapse = "; ")
  }, character(1))]

  wide[, percent_iterations_top_model := vapply(seq_len(.N), function(i) {
    vals <- percent_matrix[i, ]
    if (all(is.na(vals))) return(NA_real_)
    max(vals, na.rm = TRUE)
  }, numeric(1))]

  wide[]
}

make_multi_iter_statistic_significance_percentages <- function(
  multi_iter_model_dt,
  assignment_threshold_type_filter,
  iterations
) {
  model_stats <- c(
    "xingu_enrich",
    "belem_enrich",
    "xingu_RTH",
    "belem_RTH",
    "Xin_Bel_CC_original",
    "Tap_Xin_RTH_original"
  )
  significance_cols <- paste0("significant_", model_stats)
  output_cols <- paste0("percent_iterations_significant_", model_stats)

  required_cols <- c(
    "peak_id",
    "iteration",
    "assignment_threshold_type",
    significance_cols
  )
  if (nrow(multi_iter_model_dt) == 0 ||
      length(setdiff(required_cols, names(multi_iter_model_dt))) > 0) {
    return(data.table())
  }

  tmp <- multi_iter_model_dt[
    assignment_threshold_type == assignment_threshold_type_filter &
      iteration %in% iterations,
    c("peak_id", "iteration", significance_cols),
    with = FALSE
  ]
  if (nrow(tmp) == 0) return(data.table())

  # Each peak/iteration should occur once for a given assignment definition.
  tmp <- unique(tmp, by = c("peak_id", "iteration"))
  n_requested_iterations <- length(iterations)

  out <- tmp[, lapply(.SD, function(x) {
    value <- 100 * sum(x %in% TRUE, na.rm = TRUE) / n_requested_iterations
    if (all(is.na(x))) NA_real_ else value
  }), by = peak_id, .SDcols = significance_cols]

  setnames(out, significance_cols, output_cols)
  out[]
}

add_multi_iter_model_assignments_to_manhattan_summary <- function(
  manhattan_summary_dt,
  model_dt_2000,
  multi_iter_model_dt,
  assignment_threshold_type,
  iterations,
  all_peaks_dt
) {
  iter_wide <- make_multi_iter_model_assignment_wide(
    multi_iter_model_dt = multi_iter_model_dt,
    assignment_threshold_type_filter = assignment_threshold_type,
    iterations = iterations
  )

  stat_percentage_wide <- make_multi_iter_statistic_significance_percentages(
    multi_iter_model_dt = multi_iter_model_dt,
    assignment_threshold_type_filter = assignment_threshold_type,
    iterations = iterations
  )

  # Begin with all Fst peaks so peaks lacking ARG coverage remain represented.
  out <- all_peaks_dt[, .(
    peak_id,
    chromosome,
    peak_start = start,
    peak_end = end,
    n_outlier_windows,
    max_Xin_Bel_Fst,
    mean_Xin_Bel_Fst
  )]

  if (nrow(manhattan_summary_dt) > 0) {
    manhattan_extra_cols <- setdiff(
      names(manhattan_summary_dt),
      c("peak_id", "chromosome", "peak_start", "peak_end")
    )
    if (length(manhattan_extra_cols) > 0) {
      out <- merge(
        out,
        manhattan_summary_dt[, c("peak_id", manhattan_extra_cols), with = FALSE],
        by = "peak_id",
        all.x = TRUE,
        sort = FALSE
      )
    }
  }

  if (nrow(model_dt_2000) > 0) {
    model_2000_cols <- intersect(
      c(
        "peak_id",
        "model_assignment",
        "assignment_threshold_type",
        "min_prop_significant_local_trees"
      ),
      names(model_dt_2000)
    )
    # Use the model table as the authoritative source for the 2000th-iteration
    # assignment, avoiding duplicate model_assignment columns from a Manhattan
    # summary that may already contain the same information.
    replace_cols <- intersect(setdiff(model_2000_cols, "peak_id"), names(out))
    if (length(replace_cols) > 0) out[, (replace_cols) := NULL]
    out <- merge(
      out,
      unique(model_dt_2000[, ..model_2000_cols], by = "peak_id"),
      by = "peak_id",
      all.x = TRUE,
      sort = FALSE
    )
  }

  expected_cols <- paste0("model_assignment_iter", iterations)
  percentage_cols <- c(
    "percent_iterations_Introgression_model",
    "percent_iterations_Deep_lineage_sorting_model",
    "percent_iterations_Selection_bottleneck_model",
    "percent_iterations_Selection_recombination_model",
    "percent_iterations_Overlapping_models",
    "percent_iterations_Unassigned_to_model"
  )
  statistic_percentage_cols <- paste0(
    "percent_iterations_significant_",
    c(
      "xingu_enrich",
      "belem_enrich",
      "xingu_RTH",
      "belem_RTH",
      "Xin_Bel_CC_original",
      "Tap_Xin_RTH_original"
    )
  )
  summary_cols <- c(
    percentage_cols,
    statistic_percentage_cols,
    "top_model",
    "percent_iterations_top_model"
  )

  # Remove stale versions of the regenerated columns before merging.
  stale_cols <- intersect(c(expected_cols, summary_cols), names(out))
  if (length(stale_cols) > 0) out[, (stale_cols) := NULL]

  if (nrow(iter_wide) > 0) {
    out <- merge(out, iter_wide, by = "peak_id", all.x = TRUE, sort = FALSE)
  } else {
    for (col_i in expected_cols) out[, (col_i) := NA_character_]
    for (col_i in percentage_cols) out[, (col_i) := NA_real_]
    out[, `:=`(top_model = NA_character_, percent_iterations_top_model = NA_real_)]
  }

  if (nrow(stat_percentage_wide) > 0) {
    out <- merge(
      out,
      stat_percentage_wide,
      by = "peak_id",
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    for (col_i in statistic_percentage_cols) out[, (col_i) := NA_real_]
  }

  for (col_i in expected_cols) {
    if (!col_i %in% names(out)) out[, (col_i) := NA_character_]
  }
  for (col_i in percentage_cols) {
    if (!col_i %in% names(out)) out[, (col_i) := NA_real_]
  }
  for (col_i in statistic_percentage_cols) {
    if (!col_i %in% names(out)) out[, (col_i) := NA_real_]
  }
  if (!"top_model" %in% names(out)) out[, top_model := NA_character_]
  if (!"percent_iterations_top_model" %in% names(out)) {
    out[, percent_iterations_top_model := NA_real_]
  }

  # Explicitly remove the three columns the updated table no longer uses.
  remove_cols <- intersect(
    c(
      "model_assignments_agree_all_iterations",
      "n_iterations_with_model_assignment",
      "n_supported_models",
      "model_assignments"
    ),
    names(out)
  )
  if (length(remove_cols) > 0) out[, (remove_cols) := NULL]

  preferred_order <- c(
    "peak_id",
    "chromosome",
    "peak_start",
    "peak_end",
    "n_outlier_windows",
    "max_Xin_Bel_Fst",
    "mean_Xin_Bel_Fst",
    "panel_start",
    "panel_end",
    "model_assignment",
    expected_cols,
    percentage_cols,
    statistic_percentage_cols,
    "top_model",
    "percent_iterations_top_model"
  )
  setcolorder(out, c(intersect(preferred_order, names(out)), setdiff(names(out), preferred_order)))
  setorder(out, chromosome, peak_start, peak_end)
  out[]
}

make_50_iteration_ARG_statistic_summary <- function(
  window_counts_dt,
  peak_counts_dt,
  threshold_type_filter,
  iterations,
  statistics,
  min_prop_filter = NA_real_
) {
  filter_counts <- function(dt) {
    out <- copy(dt)[
      threshold_type == threshold_type_filter & iteration %in% iterations
    ]
    if (!is.na(min_prop_filter)) {
      out <- out[
        !is.na(min_prop_significant_local_trees) &
          abs(min_prop_significant_local_trees - min_prop_filter) < 1e-9
      ]
    }
    out
  }

  win <- filter_counts(window_counts_dt)
  peak <- filter_counts(peak_counts_dt)

  win <- win[, .(
    n_significant_windows = sum(n_significant_windows, na.rm = TRUE)
  ), by = .(statistic, iteration)]

  peak <- peak[, .(
    n_significant_peaks = sum(n_significant_peaks, na.rm = TRUE)
  ), by = .(statistic, iteration)]

  template <- CJ(statistic = statistics, iteration = as.integer(iterations), unique = TRUE)
  win <- merge(template, win, by = c("statistic", "iteration"), all.x = TRUE, sort = FALSE)
  peak <- merge(template, peak, by = c("statistic", "iteration"), all.x = TRUE, sort = FALSE)

  win[, window_col := paste0("n_significant_windows_iter", iteration)]
  peak[, peak_col := paste0("n_significant_peaks_iter", iteration)]

  win_wide <- dcast(win, statistic ~ window_col, value.var = "n_significant_windows")
  peak_wide <- dcast(peak, statistic ~ peak_col, value.var = "n_significant_peaks")

  expected_window_cols <- paste0("n_significant_windows_iter", iterations)
  expected_peak_cols <- paste0("n_significant_peaks_iter", iterations)
  for (col_i in expected_window_cols) if (!col_i %in% names(win_wide)) win_wide[, (col_i) := NA_integer_]
  for (col_i in expected_peak_cols) if (!col_i %in% names(peak_wide)) peak_wide[, (col_i) := NA_integer_]

  setcolorder(win_wide, c("statistic", expected_window_cols))
  setcolorder(peak_wide, c("statistic", expected_peak_cols))

  win_matrix <- as.matrix(win_wide[, ..expected_window_cols])
  peak_matrix <- as.matrix(peak_wide[, ..expected_peak_cols])

  safe_row_min <- function(m) apply(m, 1, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
  safe_row_max <- function(m) apply(m, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  safe_row_mean <- function(m) apply(m, 1, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
  safe_row_sd <- function(m) apply(m, 1, function(x) if (sum(!is.na(x)) < 2L) NA_real_ else sd(x, na.rm = TRUE))

  win_wide[, `:=`(
    min_n_significant_windows = safe_row_min(win_matrix),
    max_n_significant_windows = safe_row_max(win_matrix),
    average_n_significant_windows = safe_row_mean(win_matrix),
    standard_deviation_n_significant_windows = safe_row_sd(win_matrix)
  )]

  peak_wide[, `:=`(
    min_n_significant_peaks = safe_row_min(peak_matrix),
    max_n_significant_peaks = safe_row_max(peak_matrix),
    average_n_significant_peaks = safe_row_mean(peak_matrix),
    standard_deviation_n_significant_peaks = safe_row_sd(peak_matrix)
  )]

  out <- merge(win_wide, peak_wide, by = "statistic", all = TRUE, sort = FALSE)
  out[, statistic_order := match(statistic, statistics)]
  setorder(out, statistic_order)
  out[, statistic_order := NULL]

  setcolorder(out, c(
    "statistic",
    expected_window_cols,
    "min_n_significant_windows",
    "max_n_significant_windows",
    "average_n_significant_windows",
    "standard_deviation_n_significant_windows",
    expected_peak_cols,
    "min_n_significant_peaks",
    "max_n_significant_peaks",
    "average_n_significant_peaks",
    "standard_deviation_n_significant_peaks"
  ))
  out[]
}

count_significant_windows <- function(dt, stat, threshold, tail, source_name) {
  if (length(threshold) == 0 || is.na(threshold)) {
    stop("Missing threshold for statistic: ", stat)
  }

  if (tail == "upper") {
    out <- dt[
      !is.na(get(stat)),
      .(
        n_windows_with_data = uniqueN(window_id),
        n_significant_windows = uniqueN(window_id[get(stat) >= threshold])
      ),
      by = window_class
    ]
  } else if (tail == "lower") {
    out <- dt[
      !is.na(get(stat)),
      .(
        n_windows_with_data = uniqueN(window_id),
        n_significant_windows = uniqueN(window_id[get(stat) <= threshold])
      ),
      by = window_class
    ]
  } else {
    stop("Unexpected tail value: ", tail)
  }

  out[, `:=`(
    statistic = stat,
    source = source_name,
    tail = tail,
    threshold = threshold
  )]

  setcolorder(out, c(
    "statistic",
    "source",
    "window_class",
    "tail",
    "threshold",
    "n_windows_with_data",
    "n_significant_windows"
  ))

  out[]
}

count_significant_windows_min_prop <- function(
  dt,
  stat,
  threshold,
  tail,
  source_name,
  min_prop = 0.05
) {
  if (length(threshold) == 0 || is.na(threshold)) {
    stop("Missing threshold for statistic: ", stat)
  }

  if (tail == "upper") {
    out <- dt[
      !is.na(get(stat)),
      .(
        n_windows_with_data = uniqueN(window_id),
        n_local_trees = .N,
        n_significant_local_trees = sum(get(stat) >= threshold)
      ),
      by = .(window_class, window_id)
    ]
  } else if (tail == "lower") {
    out <- dt[
      !is.na(get(stat)),
      .(
        n_windows_with_data = uniqueN(window_id),
        n_local_trees = .N,
        n_significant_local_trees = sum(get(stat) <= threshold)
      ),
      by = .(window_class, window_id)
    ]
  } else {
    stop("Unexpected tail value: ", tail)
  }

  out[, prop_significant_local_trees :=
        n_significant_local_trees / n_local_trees]

  out_sig <- out[
    n_significant_local_trees >= 1 &
      prop_significant_local_trees >= min_prop
  ]

  summary_out <- out[, .(
    n_windows_with_data = .N
  ), by = window_class]

  sig_summary <- out_sig[, .(
    n_significant_windows = .N,
    mean_prop_significant_among_significant_windows = mean(prop_significant_local_trees),
    min_prop_significant_among_significant_windows = min(prop_significant_local_trees),
    max_prop_significant_among_significant_windows = max(prop_significant_local_trees)
  ), by = window_class]

  summary_out <- merge(
    summary_out,
    sig_summary,
    by = "window_class",
    all.x = TRUE
  )

  summary_out[is.na(n_significant_windows), n_significant_windows := 0L]
  summary_out[is.na(mean_prop_significant_among_significant_windows),
              mean_prop_significant_among_significant_windows := NA_real_]
  summary_out[is.na(min_prop_significant_among_significant_windows),
              min_prop_significant_among_significant_windows := NA_real_]
  summary_out[is.na(max_prop_significant_among_significant_windows),
              max_prop_significant_among_significant_windows := NA_real_]

  summary_out[, `:=`(
    statistic = stat,
    source = source_name,
    tail = tail,
    threshold = threshold,
    min_prop_significant_local_trees = min_prop
  )]

  setcolorder(summary_out, c(
    "statistic",
    "source",
    "window_class",
    "tail",
    "threshold",
    "min_prop_significant_local_trees",
    "n_windows_with_data",
    "n_significant_windows",
    "mean_prop_significant_among_significant_windows",
    "min_prop_significant_among_significant_windows",
    "max_prop_significant_among_significant_windows"
  ))

  summary_out[]
}

count_significant_peaks <- function(dt, stat, threshold, tail, source_name) {
  if (length(threshold) == 0 || is.na(threshold)) {
    stop("Missing threshold for statistic: ", stat)
  }

  if (!"peak_id" %in% names(dt)) {
    stop("Input data table must contain peak_id")
  }

  if (tail == "upper") {
    out <- dt[
      !is.na(get(stat)),
      .(
        n_peaks_with_data = uniqueN(peak_id),
        n_significant_peaks = uniqueN(peak_id[get(stat) >= threshold])
      )
    ]
  } else if (tail == "lower") {
    out <- dt[
      !is.na(get(stat)),
      .(
        n_peaks_with_data = uniqueN(peak_id),
        n_significant_peaks = uniqueN(peak_id[get(stat) <= threshold])
      )
    ]
  } else {
    stop("Unexpected tail value: ", tail)
  }

  out[, `:=`(
    statistic = stat,
    source = source_name,
    peak_class = "Xin_Bel_Fst_peak",
    tail = tail,
    threshold = threshold
  )]

  setcolorder(out, c(
    "statistic",
    "source",
    "peak_class",
    "tail",
    "threshold",
    "n_peaks_with_data",
    "n_significant_peaks"
  ))

  out[]
}

count_significant_peaks_min_prop <- function(
  dt,
  stat,
  threshold,
  tail,
  source_name,
  min_prop = 0.05
) {
  if (length(threshold) == 0 || is.na(threshold)) {
    stop("Missing threshold for statistic: ", stat)
  }

  if (!"peak_id" %in% names(dt) || !"window_id" %in% names(dt)) {
    stop("Input data table must contain peak_id and window_id")
  }

  if (tail == "upper") {
    out_window <- dt[
      !is.na(get(stat)),
      .(
        n_local_trees = .N,
        n_significant_local_trees = sum(get(stat) >= threshold)
      ),
      by = .(peak_id, window_id)
    ]
  } else if (tail == "lower") {
    out_window <- dt[
      !is.na(get(stat)),
      .(
        n_local_trees = .N,
        n_significant_local_trees = sum(get(stat) <= threshold)
      ),
      by = .(peak_id, window_id)
    ]
  } else {
    stop("Unexpected tail value: ", tail)
  }

  out_window[, prop_significant_local_trees :=
               n_significant_local_trees / n_local_trees]

  out_sig_window <- out_window[
    n_significant_local_trees >= 1 &
      prop_significant_local_trees >= min_prop
  ]

  out_peak <- out_window[, .(
    n_local_trees = sum(n_local_trees),
    n_significant_local_trees = sum(n_significant_local_trees),
    max_window_prop_significant_local_trees =
      max(prop_significant_local_trees, na.rm = TRUE)
  ), by = peak_id]

  out_sig_peak <- out_sig_window[, .(
    max_prop_significant_among_significant_windows =
      max(prop_significant_local_trees, na.rm = TRUE)
  ), by = peak_id]

  summary_out <- data.table(
    n_peaks_with_data = uniqueN(out_peak$peak_id),
    n_significant_peaks = uniqueN(out_sig_peak$peak_id),
    mean_prop_significant_among_significant_peaks =
      if (nrow(out_sig_peak) > 0) mean(out_sig_peak$max_prop_significant_among_significant_windows) else NA_real_,
    min_prop_significant_among_significant_peaks =
      if (nrow(out_sig_peak) > 0) min(out_sig_peak$max_prop_significant_among_significant_windows) else NA_real_,
    max_prop_significant_among_significant_peaks =
      if (nrow(out_sig_peak) > 0) max(out_sig_peak$max_prop_significant_among_significant_windows) else NA_real_
  )

  summary_out[, `:=`(
    statistic = stat,
    source = source_name,
    peak_class = "Xin_Bel_Fst_peak",
    tail = tail,
    threshold = threshold,
    min_prop_significant_local_trees = min_prop
  )]

  setcolorder(summary_out, c(
    "statistic",
    "source",
    "peak_class",
    "tail",
    "threshold",
    "min_prop_significant_local_trees",
    "n_peaks_with_data",
    "n_significant_peaks",
    "mean_prop_significant_among_significant_peaks",
    "min_prop_significant_among_significant_peaks",
    "max_prop_significant_among_significant_peaks"
  ))

  summary_out[]
}

make_window_significance_status <- function(
  dt,
  stat,
  threshold,
  tail,
  source_name,
  plot_label,
  window_reference,
  min_props = c(NA_real_)
) {
  if (length(threshold) == 0 || is.na(threshold)) {
    stop("Missing threshold for statistic: ", stat)
  }

  tmp <- dt[
    window_class == "Xin_Bel_Fst_outlier" &
      !is.na(get(stat)),
    .(
      n_local_trees = .N,
      n_significant_local_trees = if (tail == "upper") {
        sum(get(stat) >= threshold)
      } else if (tail == "lower") {
        sum(get(stat) <= threshold)
      } else {
        stop("Unexpected tail value: ", tail)
      }
    ),
    by = .(chromosome, window_start, window_end, window_id)
  ]

  out_base <- merge(
    window_reference,
    tmp,
    by = c("chromosome", "window_start", "window_end", "window_id"),
    all.x = TRUE,
    sort = FALSE
  )

  out_base[is.na(n_local_trees), n_local_trees := 0L]
  out_base[is.na(n_significant_local_trees), n_significant_local_trees := 0L]
  out_base[, prop_significant_local_trees := fifelse(
    n_local_trees > 0,
    n_significant_local_trees / n_local_trees,
    NA_real_
  )]

  rbindlist(lapply(min_props, function(min_prop) {
    out <- copy(out_base)

    if (is.na(min_prop)) {
      out[, threshold_type := "any_significant_local_tree"]
      out[, min_prop_significant_local_trees := NA_real_]
      out[, is_significant := n_significant_local_trees >= 1]
    } else {
      out[, threshold_type := "minimum_proportion_significant_local_trees"]
      out[, min_prop_significant_local_trees := min_prop]
      out[, is_significant :=
            n_significant_local_trees >= 1 &
            !is.na(prop_significant_local_trees) &
            prop_significant_local_trees >= min_prop]
    }

    out[, `:=`(
      statistic = stat,
      source = source_name,
      plot_label = plot_label,
      tail = tail,
      threshold = threshold
    )]

    setcolorder(out, c(
      "statistic",
      "source",
      "plot_label",
      "tail",
      "threshold",
      "threshold_type",
      "min_prop_significant_local_trees",
      "chromosome",
      "window_start",
      "window_end",
      "window_id",
      "n_local_trees",
      "n_significant_local_trees",
      "prop_significant_local_trees",
      "is_significant"
    ))

    out[]
  }))
}

make_peak_significance_status <- function(
  dt,
  stat,
  threshold,
  tail,
  source_name,
  plot_label,
  peak_reference,
  min_props = c(NA_real_)
) {
  if (length(threshold) == 0 || is.na(threshold)) {
    stop("Missing threshold for statistic: ", stat)
  }

  if (!"window_id" %in% names(dt)) {
    stop("Input data table must contain window_id for peak-level significance summaries")
  }

  # Peak-level summaries are based only on ARG local trees that fall within
  # the Fst-outlier windows that define each peak. For minimum-proportion
  # thresholds, the test is performed within each Fst-outlier window first,
  # and a peak is significant if any one of its outlier windows passes the
  # minimum-proportion rule.
  tmp_window <- dt[
    !is.na(get(stat)),
    .(
      n_local_trees = .N,
      n_significant_local_trees = if (tail == "upper") {
        sum(get(stat) >= threshold)
      } else if (tail == "lower") {
        sum(get(stat) <= threshold)
      } else {
        stop("Unexpected tail value: ", tail)
      }
    ),
    by = .(peak_id, window_id)
  ]

  tmp_window[, prop_significant_local_trees :=
               n_significant_local_trees / n_local_trees]

  tmp_peak <- tmp_window[, .(
    n_local_trees = sum(n_local_trees),
    n_significant_local_trees = sum(n_significant_local_trees),
    n_windows_with_ARG_data = .N,
    n_windows_with_any_significant_local_tree =
      sum(n_significant_local_trees >= 1),
    max_window_prop_significant_local_trees =
      max(prop_significant_local_trees, na.rm = TRUE)
  ), by = peak_id]

  out_base <- merge(
    peak_reference,
    tmp_peak,
    by = "peak_id",
    all.x = TRUE,
    sort = FALSE
  )

  out_base[is.na(n_local_trees), n_local_trees := 0L]
  out_base[is.na(n_significant_local_trees), n_significant_local_trees := 0L]
  out_base[is.na(n_windows_with_ARG_data), n_windows_with_ARG_data := 0L]
  out_base[is.na(n_windows_with_any_significant_local_tree),
           n_windows_with_any_significant_local_tree := 0L]
  out_base[, prop_significant_local_trees := fifelse(
    n_local_trees > 0,
    n_significant_local_trees / n_local_trees,
    NA_real_
  )]
  out_base[is.infinite(max_window_prop_significant_local_trees),
           max_window_prop_significant_local_trees := NA_real_]

  rbindlist(lapply(min_props, function(min_prop) {
    out <- copy(out_base)

    if (is.na(min_prop)) {
      out[, threshold_type := "any_significant_local_tree"]
      out[, min_prop_significant_local_trees := NA_real_]
      out[, is_significant := n_significant_local_trees >= 1]
    } else {
      sig_windows <- tmp_window[
        n_significant_local_trees >= 1 &
          !is.na(prop_significant_local_trees) &
          prop_significant_local_trees >= min_prop,
        .(
          n_windows_passing_min_prop = .N,
          max_prop_among_windows_passing_min_prop =
            max(prop_significant_local_trees, na.rm = TRUE)
        ),
        by = peak_id
      ]

      out <- merge(out, sig_windows, by = "peak_id", all.x = TRUE, sort = FALSE)
      out[is.na(n_windows_passing_min_prop), n_windows_passing_min_prop := 0L]
      out[is.infinite(max_prop_among_windows_passing_min_prop),
          max_prop_among_windows_passing_min_prop := NA_real_]

      out[, threshold_type := "minimum_proportion_significant_local_trees"]
      out[, min_prop_significant_local_trees := min_prop]
      out[, is_significant := n_windows_passing_min_prop >= 1]
    }

    if (!"n_windows_passing_min_prop" %in% names(out)) {
      out[, n_windows_passing_min_prop := NA_integer_]
    }
    if (!"max_prop_among_windows_passing_min_prop" %in% names(out)) {
      out[, max_prop_among_windows_passing_min_prop := NA_real_]
    }

    out[, `:=`(
      statistic = stat,
      source = source_name,
      plot_label = plot_label,
      tail = tail,
      threshold = threshold
    )]

    setcolorder(out, c(
      "statistic",
      "source",
      "plot_label",
      "tail",
      "threshold",
      "threshold_type",
      "min_prop_significant_local_trees",
      "peak_id",
      "chromosome",
      "start",
      "end",
      "n_outlier_windows",
      "n_local_trees",
      "n_significant_local_trees",
      "prop_significant_local_trees",
      "n_windows_with_ARG_data",
      "n_windows_with_any_significant_local_tree",
      "max_window_prop_significant_local_trees",
      "n_windows_passing_min_prop",
      "max_prop_among_windows_passing_min_prop",
      "is_significant"
    ))

    out[]
  }))
}

count_unique_significant_windows <- function(status_dt, min_prop = NA_real_) {
  if (is.na(min_prop)) {
    uniqueN(status_dt[
      threshold_type == "any_significant_local_tree" & is_significant == TRUE,
      window_id
    ])
  } else {
    uniqueN(status_dt[
      threshold_type == "minimum_proportion_significant_local_trees" &
        min_prop_significant_local_trees == min_prop &
        is_significant == TRUE,
      window_id
    ])
  }
}

count_unique_significant_peaks <- function(status_dt, min_prop = NA_real_) {
  if (is.na(min_prop)) {
    uniqueN(status_dt[
      threshold_type == "any_significant_local_tree" & is_significant == TRUE,
      peak_id
    ])
  } else {
    uniqueN(status_dt[
      threshold_type == "minimum_proportion_significant_local_trees" &
        min_prop_significant_local_trees == min_prop &
        is_significant == TRUE,
      peak_id
    ])
  }
}

count_entities_with_ARG_based_data <- function(coverage_dt) {
  sum(coverage_dt$has_any_ARG_based_data == TRUE, na.rm = TRUE)
}


make_peak_model_assignments <- function(
  peak_status_dt,
  peak_coverage_dt,
  min_prop = NA_real_
) {
  # Restrict to peaks with ARG-based data, because only these peaks can be
  # evaluated for ARG-based model support.
  peaks_with_arg <- peak_coverage_dt[
    has_any_ARG_based_data == TRUE,
    .(
      peak_id,
      chromosome,
      start,
      end,
      n_outlier_windows,
      n_ARG_local_trees,
      n_ARG_CC_local_trees,
      n_total_ARG_based_local_trees,
      has_ARG_data,
      has_ARG_CC_data,
      has_any_ARG_based_data,
      has_both_ARG_sources
    )
  ]

  if (nrow(peaks_with_arg) == 0) {
    return(data.table())
  }

  if (is.na(min_prop)) {
    status_subset <- peak_status_dt[
      threshold_type == "any_significant_local_tree"
    ]
    threshold_label <- "at_least_one_tree_signif_in_peak"
    min_prop_value <- NA_real_
  } else {
    status_subset <- peak_status_dt[
      threshold_type == "minimum_proportion_significant_local_trees" &
        min_prop_significant_local_trees == min_prop
    ]
    threshold_label <- paste0(
      as.integer(min_prop * 100),
      "perc_trees_signif_in_peak"
    )
    min_prop_value <- min_prop
  }

  model_stats <- c(
    "Tap_Xin_RTH_original",
    "Xin_Bel_CC_original",
    "xingu_enrich",
    "belem_enrich",
    "xingu_RTH",
    "belem_RTH"
  )

  status_wide <- dcast(
    status_subset[statistic %in% model_stats],
    peak_id ~ statistic,
    value.var = "is_significant",
    fun.aggregate = function(x) any(x == TRUE, na.rm = TRUE),
    fill = FALSE
  )

  out <- merge(
    peaks_with_arg,
    status_wide,
    by = "peak_id",
    all.x = TRUE,
    sort = FALSE
  )

  for (stat in model_stats) {
    if (!stat %in% names(out)) {
      out[, (stat) := FALSE]
    }
    out[is.na(get(stat)), (stat) := FALSE]
  }

  setnames(
    out,
    old = model_stats,
    new = paste0("significant_", model_stats)
  )

  out[, introgression_model := significant_Tap_Xin_RTH_original]

  # Deep lineage sorting requires all three peak-level signals within the
  # Fst-outlier portions of the peak: Xingu-Belem RCC plus Xingu and Belem
  # enrichment. The peak-level significance calls have already been made
  # according to the selected assignment approach.
  out[, deep_lineage_sorting_model :=
        significant_Xin_Bel_CC_original &
        significant_xingu_enrich &
        significant_belem_enrich]

  out[, selection_bottleneck_model :=
        xor(significant_xingu_enrich, significant_belem_enrich) |
        xor(significant_xingu_RTH, significant_belem_RTH)]

  # If a peak has a population-specific bottleneck-like signal plus significant
  # Xingu-Belem RCC, it should not be assigned as bottleneck-only. This is an
  # ambiguous combination of signals and is therefore forced into the
  # Overlapping models category, unless it already satisfies the stricter Deep
  # lineage sorting definition (RCC + both Xingu and Belem enrichment).
  out[, rcc_plus_selection_bottleneck_overlap :=
        significant_Xin_Bel_CC_original &
        selection_bottleneck_model &
        !deep_lineage_sorting_model]

  # Shared Xingu + Belem enrichment supports Selection-recombination only
  # in the absence of significant Xingu-Belem RCC. If RCC is significant
  # together with both enrichment signals, the peak is assigned to Deep
  # lineage sorting rather than being treated as overlapping with
  # Selection-recombination due only to shared enrichment.
  out[, selection_recombination_model :=
        ((!significant_Xin_Bel_CC_original) &
           significant_xingu_enrich & significant_belem_enrich) |
        (significant_xingu_RTH & significant_belem_RTH)]

  model_cols <- c(
    "introgression_model",
    "deep_lineage_sorting_model",
    "selection_bottleneck_model",
    "selection_recombination_model"
  )

  out[, n_supported_models := rowSums(.SD), .SDcols = model_cols]

  # Count the RCC + bottleneck conflict as an additional overlapping signal for
  # assignment purposes. This prevents peaks with significant RCC plus
  # population-specific enrichment/RTH from being labeled as bottleneck-only.
  out[
    rcc_plus_selection_bottleneck_overlap == TRUE &
      n_supported_models == 1L,
    n_supported_models := 2L
  ]

  out[, model_assignment := fifelse(
    n_supported_models == 0,
    "Unassigned to model",
    fifelse(
      n_supported_models > 1,
      "Overlapping models",
      fifelse(
        introgression_model,
        "Introgression model",
        fifelse(
          deep_lineage_sorting_model,
          "Deep lineage sorting model",
          fifelse(
            selection_bottleneck_model,
            "Selection-bottleneck model",
            "Selection-recombination model"
          )
        )
      )
    )
  )]

  out[, `:=`(
    assignment_threshold_type = threshold_label,
    min_prop_significant_local_trees = min_prop_value
  )]

  setcolorder(out, c(
    "assignment_threshold_type",
    "min_prop_significant_local_trees",
    "model_assignment",
    "n_supported_models",
    "introgression_model",
    "deep_lineage_sorting_model",
    "selection_bottleneck_model",
    "selection_recombination_model",
    "rcc_plus_selection_bottleneck_overlap",
    "significant_Tap_Xin_RTH_original",
    "significant_Xin_Bel_CC_original",
    "significant_xingu_enrich",
    "significant_belem_enrich",
    "significant_xingu_RTH",
    "significant_belem_RTH",
    "peak_id",
    "chromosome",
    "start",
    "end",
    "n_outlier_windows",
    "n_ARG_local_trees",
    "n_ARG_CC_local_trees",
    "n_total_ARG_based_local_trees",
    "has_ARG_data",
    "has_ARG_CC_data",
    "has_any_ARG_based_data",
    "has_both_ARG_sources"
  ))

  out[]
}

summarize_peak_model_assignments <- function(model_dt) {
  if (nrow(model_dt) == 0) {
    return(data.table())
  }

  model_dt[, .(
    n_peaks = .N
  ), by = .(
    assignment_threshold_type,
    min_prop_significant_local_trees,
    model_assignment
  )][, prop_peaks := n_peaks / sum(n_peaks),
     by = .(assignment_threshold_type, min_prop_significant_local_trees)][]
}



make_model_category_definitions <- function() {
  data.table(
    model_category = c(
      "Introgression model",
      "Deep lineage sorting model",
      "Selection-bottleneck model",
      "Selection-recombination model",
      "Overlapping models",
      "Unassigned to model"
    ),
    at_least_one_tree_signif_in_peak_definition = c(
      paste0(
        "Peak has ARG-based data and contains at least one ARG local tree within an Fst-outlier 10-kb window ",
        "of the focal peak with significant Tap_Xin_RTH_original."
      ),
      paste0(
        "Peak has ARG-based data and the Fst-outlier windows of the focal peak contain support for all three ",
        "signals: significant Xin_Bel_CC_original, significant xingu_enrich, and significant belem_enrich. ",
        "Each peak-level signal requires at least one significant local tree in the Fst-outlier portions of the peak."
      ),
      paste0(
        "Peak has ARG-based data and contains a population-specific signal in the Fst-outlier portions of the peak: ",
        "significant xingu_enrich or belem_enrich, but not both, or significant xingu_RTH or belem_RTH, but not both. ",
        "If significant Xin_Bel_CC_original is also present, the peak is assigned to Overlapping models rather than bottleneck-only."
      ),
      paste0(
        "Peak has ARG-based data and contains a shared Xingu-Belem signal in the Fst-outlier portions of the peak: ",
        "significant xingu_enrich and belem_enrich together only if Xin_Bel_CC_original is not significant, or ",
        "significant xingu_RTH and belem_RTH together."
      ),
      paste0(
        "Peak has ARG-based data and fulfills criteria for more than one biological model after applying the Deep ",
        "lineage sorting exception. Peaks with significant Xin_Bel_CC_original plus both enrichment signals are Deep ",
        "lineage sorting and are not Selection-recombination due only to shared enrichment. Peaks with significant ",
        "Xin_Bel_CC_original plus population-specific bottleneck-like signals are assigned to Overlapping models."
      ),
      paste0(
        "Peak has ARG-based data but does not fulfill criteria for Introgression, Deep lineage sorting, ",
        "Selection-bottleneck, or Selection-recombination."
      )
    ),
    twentyfive_percent_trees_signif_in_peak_definition = c(
      paste0(
        "Peak has ARG-based data and at least one Fst-outlier 10-kb window in the focal peak has >=25% of local ",
        "ARG trees significant for Tap_Xin_RTH_original."
      ),
      paste0(
        "Peak has ARG-based data and the Fst-outlier windows of the focal peak contain window-level support for all ",
        "three signals: Xin_Bel_CC_original, xingu_enrich, and belem_enrich. Window-level support means at least one ",
        "Fst-outlier 10-kb window has >=25% of local ARG trees significant for that statistic."
      ),
      paste0(
        "Peak has ARG-based data and at least one Fst-outlier 10-kb window in the focal peak has >=25% significant ",
        "local ARG trees for a population-specific signal: xingu_enrich or belem_enrich, but not both, or xingu_RTH ",
        "or belem_RTH, but not both. If significant Xin_Bel_CC_original is also present under the same 25% rule, ",
        "the peak is assigned to Overlapping models rather than bottleneck-only."
      ),
      paste0(
        "Peak has ARG-based data and at least one Fst-outlier 10-kb window in the focal peak has >=25% significant ",
        "local ARG trees for a shared Xingu-Belem signal: xingu_enrich and belem_enrich together only if ",
        "Xin_Bel_CC_original is not significant, or xingu_RTH and belem_RTH together."
      ),
      paste0(
        "Peak has ARG-based data and fulfills >=25% window-level criteria for more than one biological model after ",
        "applying the Deep lineage sorting exception. Peaks with significant Xin_Bel_CC_original plus both enrichment ",
        "signals are Deep lineage sorting and are not Selection-recombination due only to shared enrichment. Peaks with ",
        "significant Xin_Bel_CC_original plus population-specific bottleneck-like signals are assigned to Overlapping models."
      ),
      paste0(
        "Peak has ARG-based data but no ARG statistic fulfills the >=25% local-tree criterion within any Fst-outlier ",
        "10-kb window in the focal peak."
      )
    ),
    role_of_local_trees_windows_peaks = c(
      "Local trees are tested for significant statistic values; peak support is evaluated only within Fst-outlier windows of the focal peak.",
      "Local trees provide statistic-level evidence; windows are the 10-kb Fst-outlier windows; peaks inherit support from their Fst-outlier windows.",
      "The xingu/belem signal must be asymmetric at the peak level under the selected significance approach.",
      "The xingu/belem signal must be shared at the peak level; shared enrichment is ignored for this model when Xin_Bel_CC_original is significant.",
      "Overlapping is decided after all model booleans are computed for a peak.",
      "Only peaks with ARG-based data are assigned; peaks without ARG-based data are not evaluated for model support."
    )
  )
}

write_model_category_definitions <- function(definitions_dt, xlsx_file, tsv_file) {
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    openxlsx::write.xlsx(
      definitions_dt,
      file = xlsx_file,
      overwrite = TRUE
    )
  } else {
    warning(
      "Package 'openxlsx' is not installed; writing model-category definitions as TSV instead of XLSX: ",
      tsv_file
    )
  }

  fwrite(definitions_dt, tsv_file, sep = "\t")
}



get_status_for_min_prop <- function(status_dt, min_prop = NA_real_) {
  if (is.na(min_prop)) {
    status_dt[threshold_type == "any_significant_local_tree"]
  } else {
    status_dt[
      threshold_type == "minimum_proportion_significant_local_trees" &
        min_prop_significant_local_trees == min_prop
    ]
  }
}

get_venn_counts <- function(
  status_dt,
  entity_col,
  stat_a,
  stat_b,
  min_prop = NA_real_
) {
  tmp <- get_status_for_min_prop(status_dt, min_prop)

  set_a <- unique(tmp[
    statistic == stat_a & is_significant == TRUE,
    get(entity_col)
  ])

  set_b <- unique(tmp[
    statistic == stat_b & is_significant == TRUE,
    get(entity_col)
  ])

  both <- intersect(set_a, set_b)
  only_a <- setdiff(set_a, set_b)
  only_b <- setdiff(set_b, set_a)

  data.table(
    stat_a = stat_a,
    stat_b = stat_b,
    threshold_type = ifelse(
      is.na(min_prop),
      "any_significant_local_tree",
      "minimum_proportion_significant_local_trees"
    ),
    min_prop_significant_local_trees = min_prop,
    n_stat_a = length(set_a),
    n_stat_b = length(set_b),
    n_overlap = length(both),
    n_stat_a_only = length(only_a),
    n_stat_b_only = length(only_b),
    n_union = length(union(set_a, set_b))
  )
}

circle_intersection_area <- function(r1, r2, d) {
  if (r1 <= 0 || r2 <= 0) return(0)
  if (d >= r1 + r2) return(0)
  if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)

  part1 <- r1^2 * acos((d^2 + r1^2 - r2^2) / (2 * d * r1))
  part2 <- r2^2 * acos((d^2 + r2^2 - r1^2) / (2 * d * r2))
  part3 <- 0.5 * sqrt(
    max(
      0,
      (-d + r1 + r2) *
        (d + r1 - r2) *
        (d - r1 + r2) *
        (d + r1 + r2)
    )
  )

  part1 + part2 - part3
}

find_circle_distance <- function(r1, r2, overlap_area) {
  if (r1 <= 0 || r2 <= 0 || overlap_area <= 0) {
    return(r1 + r2 + 0.05 * max(r1, r2, 1))
  }

  max_overlap <- pi * min(r1, r2)^2
  if (overlap_area >= max_overlap) {
    return(abs(r1 - r2))
  }

  uniroot(
    function(d) circle_intersection_area(r1, r2, d) - overlap_area,
    lower = abs(r1 - r2),
    upper = r1 + r2
  )$root
}

draw_two_set_venn <- function(
  venn_counts,
  label_a,
  label_b,
  color_a,
  color_b,
  main_title,
  subtitle = NULL,
  entity_label = "features"
) {
  n_a <- venn_counts$n_stat_a[1]
  n_b <- venn_counts$n_stat_b[1]
  n_overlap <- venn_counts$n_overlap[1]
  n_a_only <- venn_counts$n_stat_a_only[1]
  n_b_only <- venn_counts$n_stat_b_only[1]
  n_union <- venn_counts$n_union[1]

  # Use a generous bottom margin so the Venn summary text does not overlap.
  par(mar = c(7.2, 4.5, 4.2, 1.5))

  plot.new()

  if ((n_a + n_b) == 0) {
    title(main = main_title, line = 1)
    if (!is.null(subtitle)) {
      mtext(subtitle, side = 3, line = 0.2, cex = 0.8)
    }
    text(0.5, 0.5, paste0("No significant ", entity_label, " for either statistic"))

    mtext(paste0(label_a, " = 0"), side = 1, line = 2.0, cex = 0.75)
    mtext(paste0(label_b, " = 0"), side = 1, line = 2.9, cex = 0.75)
    mtext("overlap = 0", side = 1, line = 3.8, cex = 0.75)
    mtext("union = 0", side = 1, line = 4.7, cex = 0.75)

    return(invisible(NULL))
  }

  r1 <- sqrt(max(n_a, 0) / pi)
  r2 <- sqrt(max(n_b, 0) / pi)
  overlap_area <- max(n_overlap, 0)
  d <- find_circle_distance(r1, r2, overlap_area)

  x1 <- -d / 2
  x2 <- d / 2
  y1 <- 0
  y2 <- 0

  xmin <- min(x1 - r1, x2 - r2)
  xmax <- max(x1 + r1, x2 + r2)
  ymin <- -max(r1, r2)
  ymax <- max(r1, r2)

  pad_x <- 0.18 * max(xmax - xmin, 1)
  pad_y <- 0.35 * max(ymax - ymin, 1)

  plot.window(
    xlim = c(xmin - pad_x, xmax + pad_x),
    ylim = c(ymin - pad_y, ymax + pad_y),
    asp = 1
  )

  theta <- seq(0, 2 * pi, length.out = 400)

  polygon(
    x1 + r1 * cos(theta),
    y1 + r1 * sin(theta),
    col = adjustcolor(color_a, alpha.f = 0.45),
    border = color_a,
    lwd = 2
  )

  polygon(
    x2 + r2 * cos(theta),
    y2 + r2 * sin(theta),
    col = adjustcolor(color_b, alpha.f = 0.45),
    border = color_b,
    lwd = 2
  )

  text(x1 - 0.35 * r1, y1, labels = n_a_only, cex = 1.25, font = 2)
  text(x2 + 0.35 * r2, y2, labels = n_b_only, cex = 1.25, font = 2)
  text((x1 + x2) / 2, 0, labels = n_overlap, cex = 1.25, font = 2)

  text(x1, y1 + r1 + 0.08 * max(r1, r2), labels = label_a, cex = 0.9, font = 2)
  text(x2, y2 + r2 + 0.08 * max(r1, r2), labels = label_b, cex = 0.9, font = 2)

  title(main = main_title, line = 1)

  # Put the cutoff/subtitle near the title instead of under the plot, so it
  # does not collide with the count summary.
  if (!is.null(subtitle)) {
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8)
  }

  # Put the summary counts on separate lines with enough spacing to avoid
  # overlap in the PDF output.
  mtext(paste0(label_a, " = ", n_a), side = 1, line = 2.0, cex = 0.75)
  mtext(paste0(label_b, " = ", n_b), side = 1, line = 2.9, cex = 0.75)
  mtext(paste0("overlap = ", n_overlap), side = 1, line = 3.8, cex = 0.75)
  mtext(paste0("union = ", n_union), side = 1, line = 4.7, cex = 0.75)

  box()
}

add_venn_pages_to_pdf <- function(
  status_dt,
  entity_col,
  min_prop = NA_real_,
  entity_label = "features",
  title_context = "Fst outlier windows"
) {
  threshold_label <- if (is.na(min_prop)) {
    "any significant local tree"
  } else {
    paste0(">=", as.integer(min_prop * 100), "% significant local trees")
  }

  enrich_counts <- get_venn_counts(
    status_dt = status_dt,
    entity_col = entity_col,
    stat_a = "xingu_enrich",
    stat_b = "belem_enrich",
    min_prop = min_prop
  )

  draw_two_set_venn(
    venn_counts = enrich_counts,
    label_a = "Xingu enrichment",
    label_b = "Belem enrichment",
    color_a = "#6ea8db",
    color_b = "#fdbb4a",
    main_title = paste0(title_context, ": Xingu vs Belem enrichment"),
    subtitle = threshold_label,
    entity_label = entity_label
  )

  rth_counts <- get_venn_counts(
    status_dt = status_dt,
    entity_col = entity_col,
    stat_a = "xingu_RTH",
    stat_b = "belem_RTH",
    min_prop = min_prop
  )

  draw_two_set_venn(
    venn_counts = rth_counts,
    label_a = "Xingu RTH'",
    label_b = "Belem RTH'",
    color_a = "#6ea8db",
    color_b = "#fdbb4a",
    main_title = paste0(title_context, ": Xingu vs Belem RTH'"),
    subtitle = threshold_label,
    entity_label = entity_label
  )
}

summarize_venn_counts_for_thresholds <- function(
  status_dt,
  entity_col,
  entity_level,
  min_props
) {
  rbindlist(lapply(c(NA_real_, min_props), function(min_prop) {
    rbind(
      cbind(
        data.table(entity_level = entity_level, comparison = "xingu_enrich_vs_belem_enrich"),
        get_venn_counts(status_dt, entity_col, "xingu_enrich", "belem_enrich", min_prop)
      ),
      cbind(
        data.table(entity_level = entity_level, comparison = "xingu_RTH_vs_belem_RTH"),
        get_venn_counts(status_dt, entity_col, "xingu_RTH", "belem_RTH", min_prop)
      ),
      fill = TRUE
    )
  }), fill = TRUE)
}

make_rth_decile_distribution_barplot <- function(plot_dt, outfile) {
  if (exists("CREATE_BARPLOT_PDFS") && !isTRUE(CREATE_BARPLOT_PDFS)) {
    return(invisible(NULL))
  }

  pdf(outfile, width = 8.5, height = 5.5)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  })

  stats <- unique(plot_dt$statistic)

  for (stat in stats) {
    tmp <- plot_dt[statistic == stat]
    classes <- unique(tmp$window_class)

    ymax <- max(tmp$percent_local_trees, na.rm = TRUE) * 1.15
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1

    class_colors <- c(
      control = "grey70",
      Xin_Bel_Fst_outlier = "grey30"
    )

    mat <- dcast(
      tmp,
      bin_index + bin_label ~ window_class,
      value.var = "percent_local_trees",
      fill = 0
    )

    setorder(mat, bin_index)

    ymat <- as.matrix(mat[, ..classes])
    ymat <- t(ymat)

    par(mar = c(7.5, 4.5, 3, 1))

    bp <- barplot(
      ymat,
      beside = TRUE,
      col = class_colors[classes],
      border = "black",
      ylim = c(0, ymax),
      ylab = "Percent of local trees",
      xlab = "",
      main = paste0(unname(plot_labels[stat]), " distribution"),
      names.arg = mat$bin_label,
      las = 2,
      cex.names = 0.8
    )

    legend(
      "topright",
      legend = classes,
      fill = class_colors[classes],
      border = "black",
      bty = "n",
      cex = 0.85
    )
  }
}

make_significant_window_barplot <- function(
  plot_dt,
  outfile,
  y_col = "n_significant_windows",
  ylab = "Number of Xin-Bel Fst outlier windows",
  main_title = "Fst outlier windows containing significant ARG local trees",
  annotation_lines = NULL,
  venn_status_dt = NULL,
  venn_entity_col = NULL,
  venn_min_prop = NA_real_,
  venn_entity_label = "features",
  venn_title_context = main_title,
  interval_dt = NULL,
  interval_low_col = "interval_low",
  interval_high_col = "interval_high",
  interval_label = "MCMC iteration range"
) {
  if (exists("CREATE_BARPLOT_PDFS") && !isTRUE(CREATE_BARPLOT_PDFS)) {
    return(invisible(NULL))
  }

  pdf(outfile, width = 7.5, height = 5.5)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  })

  par(mar = c(8.5, 4.5, 2, 1))

  plot_dt <- copy(plot_dt)
  plot_dt[, .barplot_row_order := .I]

  if (!is.null(interval_dt) && nrow(interval_dt) > 0) {
    interval_cols_required <- c("statistic", interval_low_col, interval_high_col)
    missing_interval_cols <- setdiff(interval_cols_required, names(interval_dt))
    if (length(missing_interval_cols) > 0) {
      warning(
        "Interval data are missing required column(s): ",
        paste(missing_interval_cols, collapse = ", "),
        ". Plotting bars without intervals for ",
        outfile
      )
      interval_dt <- NULL
    } else {
      interval_subset <- copy(interval_dt[, ..interval_cols_required])
      setnames(
        interval_subset,
        old = c(interval_low_col, interval_high_col),
        new = c("bar_interval_low", "bar_interval_high")
      )
      plot_dt <- merge(
        plot_dt,
        interval_subset,
        by = "statistic",
        all.x = TRUE,
        sort = FALSE
      )
      setorder(plot_dt, .barplot_row_order)
    }
  }

  if (!"bar_interval_low" %in% names(plot_dt)) {
    plot_dt[, bar_interval_low := NA_real_]
  }
  if (!"bar_interval_high" %in% names(plot_dt)) {
    plot_dt[, bar_interval_high := NA_real_]
  }

  y <- plot_dt[[y_col]]
  x <- seq_along(y)
  ymax <- max(c(y, plot_dt$bar_interval_high), na.rm = TRUE) * 1.18
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1

  plot(
    x,
    y,
    type = "n",
    xaxt = "n",
    xlab = "",
    ylab = ylab,
    ylim = c(0, ymax),
    xlim = c(0.4, length(x) + 0.6),
    main = main_title
  )

  bar_width <- 0.72

  draw_striped_bar <- function(xmid, height, col1, col2) {
    left <- xmid - bar_width / 2
    right <- xmid + bar_width / 2

    rect(left, 0, right, height, col = "white", border = "black")

    old_clip <- par("usr")
    clip(left, right, 0, height)

    spacing <- max(height / 14, 1)
    starts <- seq(-height, height + spacing * 4, by = spacing)

    for (j in seq_along(starts)) {
      col_use <- if (j %% 2 == 1) col1 else col2
      segments(
        x0 = left,
        y0 = starts[j],
        x1 = right,
        y1 = starts[j] + height,
        col = col_use,
        lwd = 5,
        lend = "butt"
      )
    }

    do.call(clip, as.list(old_clip))
    rect(left, 0, right, height, col = NA, border = "black")
  }

  for (i in seq_along(x)) {
    stat <- plot_dt$statistic[i]
    if (stat == "Xin_Bel_CC_original") {
      draw_striped_bar(x[i], y[i], "#6ea8db", "#fdbb4a")
    } else if (stat == "Tap_Xin_RTH_original") {
      draw_striped_bar(x[i], y[i], "#89ca8e", "#6ea8db")
    } else {
      rect(
        x[i] - bar_width / 2,
        0,
        x[i] + bar_width / 2,
        y[i],
        col = plot_dt$fill_color[i],
        border = "black"
      )
    }
  }

  has_intervals <- any(
    !is.na(plot_dt$bar_interval_low) &
      !is.na(plot_dt$bar_interval_high)
  )

  if (has_intervals) {
    interval_cap_width <- bar_width * 0.42
    for (i in seq_along(x)) {
      low_i <- plot_dt$bar_interval_low[i]
      high_i <- plot_dt$bar_interval_high[i]
      if (is.na(low_i) || is.na(high_i)) next

      segments(
        x0 = x[i],
        y0 = low_i,
        x1 = x[i],
        y1 = high_i,
        lwd = 1.4,
        col = "black"
      )
      segments(
        x0 = x[i] - interval_cap_width / 2,
        y0 = low_i,
        x1 = x[i] + interval_cap_width / 2,
        y1 = low_i,
        lwd = 1.4,
        col = "black"
      )
      segments(
        x0 = x[i] - interval_cap_width / 2,
        y0 = high_i,
        x1 = x[i] + interval_cap_width / 2,
        y1 = high_i,
        lwd = 1.4,
        col = "black"
      )
    }

    legend(
      "topleft",
      legend = interval_label,
      lty = 1,
      lwd = 1.4,
      bty = "n",
      cex = 0.72
    )
  }

  axis(
    1,
    at = x,
    labels = FALSE
  )

  text(
    x = x,
    y = par("usr")[3] - 0.055 * diff(par("usr")[3:4]),
    labels = plot_dt$plot_label,
    srt = 45,
    adj = 1,
    xpd = TRUE,
    cex = 0.8
  )

  label_y <- ifelse(
    !is.na(plot_dt$bar_interval_high),
    pmax(y, plot_dt$bar_interval_high, na.rm = TRUE),
    y
  )

  text(
    x = x,
    y = label_y,
    labels = y,
    pos = 3,
    cex = 0.75
  )

  if (!is.null(annotation_lines) && length(annotation_lines) > 0) {
    usr <- par("usr")
    text(
      x = usr[2] - 0.02 * diff(usr[1:2]),
      y = usr[4] - 0.08 * diff(usr[3:4]),
      labels = paste(annotation_lines, collapse = "\n"),
      adj = c(1, 1),
      cex = 0.72,
      xpd = TRUE
    )
  }

  if (!is.null(venn_status_dt) && !is.null(venn_entity_col)) {
    add_venn_pages_to_pdf(
      status_dt = venn_status_dt,
      entity_col = venn_entity_col,
      min_prop = venn_min_prop,
      entity_label = venn_entity_label,
      title_context = venn_title_context
    )
  }
}


make_multi_iter_barplot_interval_dt <- function(
  multi_iter_counts_dt,
  y_col,
  threshold_type_value,
  min_prop_value = NA_real_,
  iterations = NULL
) {
  if (is.null(multi_iter_counts_dt) || nrow(multi_iter_counts_dt) == 0) {
    return(data.table())
  }

  required_cols <- c("statistic", "threshold_type", y_col, "iteration")
  missing_cols <- setdiff(required_cols, names(multi_iter_counts_dt))
  if (length(missing_cols) > 0) {
    warning(
      "Cannot make multi-iteration barplot intervals; missing column(s): ",
      paste(missing_cols, collapse = ", ")
    )
    return(data.table())
  }

  tmp <- copy(multi_iter_counts_dt)
  tmp <- tmp[threshold_type == threshold_type_value]
  if (!is.null(iterations)) {
    tmp <- tmp[iteration %in% iterations]
  }

  if (!is.na(min_prop_value)) {
    if (!"min_prop_significant_local_trees" %in% names(tmp)) {
      warning(
        "Cannot filter multi-iteration intervals by minimum proportion; ",
        "column min_prop_significant_local_trees is missing."
      )
      return(data.table())
    }
    tmp <- tmp[
      !is.na(min_prop_significant_local_trees) &
        abs(min_prop_significant_local_trees - min_prop_value) < 1e-9
    ]
  }

  if (nrow(tmp) == 0) {
    return(data.table())
  }

  out <- tmp[, .(
    interval_low = min(get(y_col), na.rm = TRUE),
    interval_high = max(get(y_col), na.rm = TRUE),
    n_iterations_in_interval = uniqueN(iteration)
  ), by = statistic]

  out[!is.finite(interval_low), interval_low := NA_real_]
  out[!is.finite(interval_high), interval_high := NA_real_]
  out[]
}


make_ARG_threshold_summary_for_p_value <- function(
  p_value,
  arg_joined_dt,
  arg_cc_joined_dt,
  barplot_stats_dt
) {
  if (is.na(p_value) || p_value <= 0 || p_value >= 1) {
    stop("p_value must be between 0 and 1")
  }

  rbindlist(lapply(seq_len(nrow(barplot_stats_dt)), function(i) {
    row <- barplot_stats_dt[i]

    control_dt <- if (row$source == "ARG") {
      arg_joined_dt[window_class == "control"]
    } else if (row$source == "ARG_CC") {
      arg_cc_joined_dt[window_class == "control"]
    } else {
      stop("Unexpected source in barplot_stats: ", row$source)
    }

    values <- control_dt[[row$statistic]]
    values <- values[!is.na(values)]

    if (row$tail == "upper") {
      prob <- 1 - p_value
      percentile <- 100 * (1 - p_value)
    } else if (row$tail == "lower") {
      prob <- p_value
      percentile <- 100 * p_value
    } else {
      stop("Unexpected tail value: ", row$tail)
    }

    data.table(
      statistic = row$statistic,
      source = row$source,
      tail = row$tail,
      p_value = p_value,
      quantile_probability = prob,
      percentile = percentile,
      threshold = if (length(values) > 0) {
        as.numeric(quantile(values, probs = prob, na.rm = TRUE))
      } else {
        NA_real_
      },
      n_values = length(values)
    )
  }), fill = TRUE)
}

make_ARG_barplot_counts_for_p_value <- function(
  p_value,
  entity_level = c("window", "peak"),
  significance_definition = c("any_significant_local_tree", "minimum_proportion_significant_local_trees"),
  min_prop = NA_real_,
  arg_joined_dt,
  arg_cc_joined_dt,
  arg_joined_peaks_dt,
  arg_cc_joined_peaks_dt,
  barplot_stats_dt
) {
  entity_level <- match.arg(entity_level)
  significance_definition <- match.arg(significance_definition)

  thresholds_dt <- make_ARG_threshold_summary_for_p_value(
    p_value = p_value,
    arg_joined_dt = arg_joined_dt,
    arg_cc_joined_dt = arg_cc_joined_dt,
    barplot_stats_dt = barplot_stats_dt
  )

  counts_dt <- rbindlist(lapply(seq_len(nrow(barplot_stats_dt)), function(i) {
    row <- barplot_stats_dt[i]
    threshold_row <- thresholds_dt[
      statistic == row$statistic &
        source == row$source &
        tail == row$tail
    ]

    if (nrow(threshold_row) == 0 || is.na(threshold_row$threshold[1])) {
      stop("Missing threshold for statistic: ", row$statistic, " at P = ", p_value)
    }

    if (entity_level == "window") {
      input_dt <- if (row$source == "ARG") arg_joined_dt else arg_cc_joined_dt

      if (significance_definition == "any_significant_local_tree") {
        count_significant_windows(
          dt = input_dt,
          stat = row$statistic,
          threshold = threshold_row$threshold[1],
          tail = row$tail,
          source_name = row$source
        )[window_class == "Xin_Bel_Fst_outlier"]
      } else {
        count_significant_windows_min_prop(
          dt = input_dt,
          stat = row$statistic,
          threshold = threshold_row$threshold[1],
          tail = row$tail,
          source_name = row$source,
          min_prop = min_prop
        )[window_class == "Xin_Bel_Fst_outlier"]
      }
    } else {
      input_dt <- if (row$source == "ARG") arg_joined_peaks_dt else arg_cc_joined_peaks_dt

      if (significance_definition == "any_significant_local_tree") {
        count_significant_peaks(
          dt = input_dt,
          stat = row$statistic,
          threshold = threshold_row$threshold[1],
          tail = row$tail,
          source_name = row$source
        )
      } else {
        count_significant_peaks_min_prop(
          dt = input_dt,
          stat = row$statistic,
          threshold = threshold_row$threshold[1],
          tail = row$tail,
          source_name = row$source,
          min_prop = min_prop
        )
      }
    }
  }), fill = TRUE)

  counts_dt <- merge(
    counts_dt,
    barplot_stats_dt[, .(
      statistic,
      direction_label,
      plot_label,
      fill_color
    )],
    by = "statistic",
    all.x = TRUE,
    sort = FALSE
  )

  counts_dt <- merge(
    counts_dt,
    thresholds_dt[, .(
      statistic,
      p_value_for_threshold = p_value,
      threshold_percentile = percentile,
      quantile_probability,
      threshold
    )],
    by = c("statistic", "threshold"),
    all.x = TRUE,
    sort = FALSE
  )

  counts_dt[, statistic_order := match(statistic, barplot_stats_dt$statistic)]
  setorder(counts_dt, statistic_order)
  counts_dt[]
}

format_p_sensitivity_title <- function(p_value) {
  if (abs(p_value - 0.0001) < 1e-12) {
    "P = 0.0001 (0.01st / 99.99th percentiles)"
  } else if (abs(p_value - 0.001) < 1e-12) {
    "P = 0.001 (0.1st / 99.9th percentiles)"
  } else if (abs(p_value - 0.01) < 1e-12) {
    "P = 0.01 (1st / 99th percentiles)"
  } else {
    paste0("P = ", signif(p_value, 4))
  }
}

draw_ARG_p_sensitivity_barplot_panel <- function(
  plot_dt,
  y_col,
  ylab,
  main_title
) {
  plot_dt <- copy(plot_dt)
  setorder(plot_dt, statistic_order)

  y <- plot_dt[[y_col]]
  x <- seq_along(y)
  ymax <- max(y, na.rm = TRUE) * 1.22
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1

  plot(
    x,
    y,
    type = "n",
    xaxt = "n",
    xlab = "",
    ylab = ylab,
    ylim = c(0, ymax),
    xlim = c(0.4, length(x) + 0.6),
    main = main_title,
    cex.main = 0.9,
    cex.lab = 0.85,
    cex.axis = 0.8
  )

  bar_width <- 0.72

  draw_striped_bar <- function(xmid, height, col1, col2) {
    left <- xmid - bar_width / 2
    right <- xmid + bar_width / 2

    rect(left, 0, right, height, col = "white", border = "black")

    old_clip <- par("usr")
    clip(left, right, 0, height)

    spacing <- max(height / 14, 1)
    starts <- seq(-height, height + spacing * 4, by = spacing)

    for (j in seq_along(starts)) {
      col_use <- if (j %% 2 == 1) col1 else col2
      segments(
        x0 = left,
        y0 = starts[j],
        x1 = right,
        y1 = starts[j] + height,
        col = col_use,
        lwd = 4.2,
        lend = "butt"
      )
    }

    do.call(clip, as.list(old_clip))
    rect(left, 0, right, height, col = NA, border = "black")
  }

  for (i in seq_along(x)) {
    stat <- plot_dt$statistic[i]
    if (stat == "Xin_Bel_CC_original") {
      draw_striped_bar(x[i], y[i], "#6ea8db", "#fdbb4a")
    } else if (stat == "Tap_Xin_RTH_original") {
      draw_striped_bar(x[i], y[i], "#89ca8e", "#6ea8db")
    } else {
      rect(
        x[i] - bar_width / 2,
        0,
        x[i] + bar_width / 2,
        y[i],
        col = plot_dt$fill_color[i],
        border = "black"
      )
    }
  }

  axis(1, at = x, labels = FALSE)
  text(
    x = x,
    y = par("usr")[3] - 0.06 * diff(par("usr")[3:4]),
    labels = plot_dt$plot_label,
    srt = 45,
    adj = 1,
    xpd = TRUE,
    cex = 0.68
  )

  text(
    x = x,
    y = y,
    labels = y,
    pos = 3,
    cex = 0.66
  )

  box()
}

make_ARG_p_sensitivity_combined_barplot_pdf <- function(
  outfile,
  significance_definition = c("any_significant_local_tree", "minimum_proportion_significant_local_trees"),
  p_values = c(0.0001, 0.001, 0.01),
  min_prop = NA_real_,
  arg_joined_dt,
  arg_cc_joined_dt,
  arg_joined_peaks_dt,
  arg_cc_joined_peaks_dt,
  barplot_stats_dt
) {
  significance_definition <- match.arg(significance_definition)

  pdf(outfile, width = 12, height = 13.2, useDingbats = FALSE)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  layout(matrix(seq_len(length(p_values) * 2), nrow = length(p_values), ncol = 2, byrow = TRUE))

  for (p_value in p_values) {
    window_counts <- make_ARG_barplot_counts_for_p_value(
      p_value = p_value,
      entity_level = "window",
      significance_definition = significance_definition,
      min_prop = min_prop,
      arg_joined_dt = arg_joined_dt,
      arg_cc_joined_dt = arg_cc_joined_dt,
      arg_joined_peaks_dt = arg_joined_peaks_dt,
      arg_cc_joined_peaks_dt = arg_cc_joined_peaks_dt,
      barplot_stats_dt = barplot_stats_dt
    )

    peak_counts <- make_ARG_barplot_counts_for_p_value(
      p_value = p_value,
      entity_level = "peak",
      significance_definition = significance_definition,
      min_prop = min_prop,
      arg_joined_dt = arg_joined_dt,
      arg_cc_joined_dt = arg_cc_joined_dt,
      arg_joined_peaks_dt = arg_joined_peaks_dt,
      arg_cc_joined_peaks_dt = arg_cc_joined_peaks_dt,
      barplot_stats_dt = barplot_stats_dt
    )

    if (significance_definition == "any_significant_local_tree") {
      definition_label <- "any significant ARG local tree"
    } else {
      definition_label <- paste0(
        ">=",
        as.integer(min_prop * 100),
        "% significant ARG local trees"
      )
    }

    par(mar = c(7.8, 4.7, 3.0, 1.1))
    draw_ARG_p_sensitivity_barplot_panel(
      plot_dt = window_counts,
      y_col = "n_significant_windows",
      ylab = "Number of Xin-Bel Fst outlier windows",
      main_title = paste0("Windows: ", definition_label, "\n", format_p_sensitivity_title(p_value))
    )

    par(mar = c(7.8, 4.7, 3.0, 1.1))
    draw_ARG_p_sensitivity_barplot_panel(
      plot_dt = peak_counts,
      y_col = "n_significant_peaks",
      ylab = "Number of Xin-Bel Fst peaks",
      main_title = paste0("Peaks: ", definition_label, "\n", format_p_sensitivity_title(p_value))
    )
  }

  invisible(NULL)
}


# -----------------------------
# Plotting labels for later use
# -----------------------------
plot_labels <- c(
  Tap_Xin_Fst = "Tapajos-Xingu Fst",
  Tap_Bel_Fst = "Tapajos-Belem Fst",
  Xin_Bel_Fst = "Xingu-Belem Fst",
  Tap_Xin_Dxy = "Tapajos-Xingu Dxy",
  Tap_Bel_Dxy = "Tapajos-Belem Dxy",
  Xin_Bel_Dxy = "Xingu-Belem Dxy",
  Tap_pi = "Tapajos pi",
  Xin_pi = "Xingu pi",
  Bel_pi = "Belem pi",

  tapajos_enrich = "Tapajos enrichment",
  xingu_enrich = "Xingu enrichment",
  belem_enrich = "Belem enrichment",

  tapajos_RTH = "Tapajos RTH'",
  xingu_RTH = "Xingu RTH'",
  belem_RTH = "Belem RTH'",
  Tap_Xin_RTH_prime = "Tapajos-Xingu RTH'",
  Tap_Xin_RTH_original = "Tapajos-Xingu RTH",

  tapajos_RTH_inverse = "Tapajos RTH' (inverse)",
  xingu_RTH_inverse = "Xingu RTH' (inverse)",
  belem_RTH_inverse = "Belem RTH' (inverse)",
  Tap_Xin_RTH_prime_inverse = "Tapajos-Xingu RTH' (inverse)",
  Tap_Xin_RTH_original_inverse = "Tapajos-Xingu RTH (inverse)",

  Tapajos_Xingu_fd = "Tapajos-Xingu fd",
  recombination_rate = "recombination rate",
  Xingu_RAiSD_u = "Xingu RAiSD u",
  xingu_norm_ihs = "Xingu iHS",
  xingu_norm_nsl = "Xingu nSL",
  Belem_RAiSD_u = "Belem RAiSD u",
  belem_norm_ihs = "Belem iHS",
  belem_norm_nsl = "Belem nSL",
  dup_sites = "rCNV duplicate sites",
  perc_dup_sites = "rCNV % duplicate sites",
  Tap_Xin_Fst_ARG_based = "Tapajos-Xingu Fst (ARG-based)",
  Tap_Bel_Fst_ARG_based = "Tapajos-Belem Fst (ARG-based)",
  Xin_Bel_Fst_ARG_based = "Xingu-Belem Fst (ARG-based)",

  Xin_Bel_CC_original = "Xingu-Belem RCC",
  Xin_Bel_CC_prime = "Xingu-Belem RCC'",
  Tap_Xin_JCR_v2 = "Tapajos-Xingu JCR v2"
)

# -----------------------------
# Read ARGweaver mask/problem regions
# -----------------------------
argweaver_mask <- read_argweaver_mask(argweaver_mask_file)

# -----------------------------
# Read repeat regions for Manhattan plot annotation
# -----------------------------
repeat_regions_for_manhattan <- read_repeat_regions_for_plotting(repeat_regions_file)

# -----------------------------
# Read pixy data
# -----------------------------
df <- fread(pixy_file)

# -----------------------------
# Convert negative Fst values to zero
# -----------------------------
fst_cols <- c("Tap_Xin_Fst", "Tap_Bel_Fst", "Xin_Bel_Fst")

for (col in fst_cols) {
  df[get(col) < 0, (col) := 0]
}

# -----------------------------
# Filter to autosomes and windows with >=50 SNPs
# -----------------------------
df_clean <- df[
  chrom_type != "sex_chromosome" &
    no_snps_fst >= 50
]

# -----------------------------
# Calculate Xin-Bel Fst thresholds
# -----------------------------
mean_fst <- mean(df_clean$Xin_Bel_Fst, na.rm = TRUE)
sd_fst   <- sd(df_clean$Xin_Bel_Fst, na.rm = TRUE)

threshold <- mean_fst + 5 * sd_fst
control_threshold <- mean_fst + sd_fst

# -----------------------------
# Define outlier, control, and intermediate windows
# -----------------------------
outliers <- df_clean[
  !is.na(Xin_Bel_Fst) &
    Xin_Bel_Fst > threshold
]

control_windows <- df_clean[
  !is.na(Xin_Bel_Fst) &
    Xin_Bel_Fst <= control_threshold
]

intermediate_windows <- df_clean[
  !is.na(Xin_Bel_Fst) &
    Xin_Bel_Fst > control_threshold &
    Xin_Bel_Fst <= threshold
]

# -----------------------------
# Write outlier windows
# -----------------------------
outlier_bed <- outliers[, .(
  chromosome,
  start,
  end,
  Xin_Bel_Fst
)]

fwrite(outlier_bed, outlier_windows_bed, sep = "\t", col.names = TRUE)

# -----------------------------
# Store coordinates for outlier and control windows
# -----------------------------
outlier_coords <- outliers[, .(
  chromosome,
  start,
  end,
  window_class = "Xin_Bel_Fst_outlier"
)]

control_coords <- control_windows[, .(
  chromosome,
  start,
  end,
  window_class = "control"
)]

comparison_windows <- rbind(outlier_coords, control_coords)
comparison_windows[, window_id := paste(chromosome, start, end, sep = ":")]

fwrite(comparison_windows, comparison_windows_file, sep = "\t")

# -----------------------------
# Summarize Fst, Dxy, and pi for outlier and control windows
# -----------------------------
stats_to_summarize <- c(
  "Xin_Bel_Fst",
  "Tap_Xin_Dxy",
  "Tap_Bel_Dxy",
  "Xin_Bel_Dxy",
  "Tap_pi",
  "Xin_pi",
  "Bel_pi"
)

window_stat_summary <- rbind(
  summarize_stats(outliers, "Xin_Bel_Fst_outlier", stats_to_summarize),
  summarize_stats(control_windows, "control", stats_to_summarize)
)

fwrite(window_stat_summary, window_stat_summary_file, sep = "\t")

# -----------------------------
# Empirical pixy significance thresholds from control windows
# -----------------------------
pixy_high_tail_stats <- c(
  "Tap_Xin_Fst",
  "Tap_Bel_Fst",
  "Xin_Bel_Fst",
  "Tap_Xin_Dxy",
  "Tap_Bel_Dxy",
  "Xin_Bel_Dxy"
)

pixy_low_tail_stats <- c(
  "Tap_pi",
  "Xin_pi",
  "Bel_pi"
)

pixy_threshold_summary <- rbindlist(c(
  lapply(pixy_high_tail_stats, function(stat) {
    data.table(
      statistic = stat,
      threshold_basis = "pixy_control_windows",
      tail = "upper",
      p_value = 0.001,
      percentile = 99.9,
      threshold = as.numeric(quantile(
        control_windows[[stat]],
        probs = 0.999,
        na.rm = TRUE
      )),
      n_values = sum(!is.na(control_windows[[stat]]))
    )
  }),

  lapply(pixy_low_tail_stats, function(stat) {
    data.table(
      statistic = stat,
      threshold_basis = "pixy_control_windows",
      tail = "lower",
      p_value = 0.001,
      percentile = 0.1,
      threshold = as.numeric(quantile(
        control_windows[[stat]],
        probs = 0.001,
        na.rm = TRUE
      )),
      n_values = sum(!is.na(control_windows[[stat]]))
    )
  })
))

fwrite(pixy_threshold_summary, pixy_threshold_summary_file, sep = "\t")

# =============================================================
# Additional genome-wide statistics: control-window thresholds
# =============================================================
# The thresholds below are calculated from the distribution of values
# overlapping the same Xin-Bel Fst control windows defined above.
# For interval/window files, intervals are treated as BED-style intervals
# [start, end). For D-statistic windows, which are 1-based, start is
# converted to BED-style coordinates by subtracting 1.

control_interval_dt <- control_windows[, .(
  chromosome,
  window_start = as.integer(start),
  window_end_inclusive = as.integer(end - 1L),
  window_id = paste(chromosome, start, end, sep = ":")
)]

additional_threshold_list <- list()
additional_value_summary_list <- list()

add_additional_thresholds <- function(values, stat, threshold_basis, tails) {
  additional_threshold_list[[length(additional_threshold_list) + 1L]] <<-
    make_threshold_rows(values, stat, threshold_basis, tails)

  additional_value_summary_list[[length(additional_value_summary_list) + 1L]] <<-
    make_value_summary(values, stat, threshold_basis)
}

add_recombination_rate_threshold <- function(values) {
  values <- values[!is.na(values)]

  additional_threshold_list[[length(additional_threshold_list) + 1L]] <<-
    data.table(
      statistic = "recombination_rate",
      threshold_basis = "ReLERNN_control_windows",
      tail = "lower",
      p_value = 0.005,
      percentile = 0.5,
      threshold = if (length(values) > 0) {
        as.numeric(quantile(values, probs = 0.005, na.rm = TRUE))
      } else {
        NA_real_
      },
      n_values = length(values)
    )

  additional_value_summary_list[[length(additional_value_summary_list) + 1L]] <<-
    make_value_summary(values, "recombination_rate", "ReLERNN_control_windows")
}


# Store selection-statistic flags while the raw selection datasets are in memory.
# These flags are NOT based on 10-kb averaged selection values. Instead, they
# record whether any raw selection-statistic value exceeds its threshold inside
# each 10-kb Fst-outlier window; later they are collapsed to peak-level flags.
selection_window_significance_list <- list()

# -----------------------------
# D statistics: Tapajos-Xingu fd
# -----------------------------
d_stats <- fread(
  d_stats_file,
  select = c("scaffold", "start", "end", "sitesUsed", "D", "fd")
)

setnames(
  d_stats,
  old = c("scaffold", "fd"),
  new = c("chromosome", "Tapajos_Xingu_fd")
)

# fd is not meaningful when D is negative, so set fd to zero for those windows.
d_stats[!is.na(D) & D < 0, Tapajos_Xingu_fd := 0]

# Exclude biologically invalid fd values before all downstream summaries,
# thresholds, heat maps, Manhattan plots, and correlation analyses.
d_stats <- d_stats[
  is.na(Tapajos_Xingu_fd) |
    (Tapajos_Xingu_fd >= 0 & Tapajos_Xingu_fd <= 1)
]

d_stats <- d_stats[sitesUsed >= 100]
d_stats[, interval_start := as.integer(start - 1L)]
d_stats[, interval_end := as.integer(end)]

values <- extract_control_values_from_intervals(
  d_stats[, .(chromosome, interval_start, interval_end, Tapajos_Xingu_fd)],
  "Tapajos_Xingu_fd",
  control_interval_dt
)
add_additional_thresholds(
  values,
  "Tapajos_Xingu_fd",
  "D_statistics_control_windows_sitesUsed_ge_100",
  "upper"
)
rm(d_stats, values)
gc()

# -----------------------------
# ReLERNN recombination rate
# -----------------------------
relernn <- fread(
  relernn_file,
  header = FALSE,
  col.names = c("chromosome", "interval_start", "interval_end", "recombination_rate")
)

region_set_pi_recombination_summary <- summarize_pi_recombination_by_region_set(
  all_windows_dt = df_clean,
  outlier_windows_dt = outliers,
  control_windows_dt = control_windows,
  recombination_dt = relernn
)

fwrite(
  region_set_pi_recombination_summary,
  region_set_pi_recombination_summary_file,
  sep = "\t"
)

values <- extract_control_values_from_intervals(
  relernn,
  "recombination_rate",
  control_interval_dt
)
add_recombination_rate_threshold(values)
rm(relernn, values)
gc()

# -----------------------------
# RAiSD U statistics
# -----------------------------
xingu_raisd <- fread(
  xingu_raisd_file,
  header = FALSE,
  col.names = c("chromosome", "interval_start", "interval_end", "Xingu_RAiSD_u")
)

values <- extract_control_values_from_intervals(
  xingu_raisd,
  "Xingu_RAiSD_u",
  control_interval_dt
)
xingu_raisd_threshold <- as.numeric(quantile(values, probs = 0.999, na.rm = TRUE))
selection_window_significance_list[["Xingu_RAiSD_u"]] <-
  make_outlier_window_selection_flags_from_intervals(
    xingu_raisd,
    "Xingu_RAiSD_u",
    outliers,
    threshold = xingu_raisd_threshold,
    use_abs = FALSE
  )
add_additional_thresholds(
  values,
  "Xingu_RAiSD_u",
  "Xingu_RAiSD_control_windows",
  "upper"
)
rm(xingu_raisd, values, xingu_raisd_threshold)
gc()

belem_raisd <- fread(
  belem_raisd_file,
  header = FALSE,
  col.names = c("chromosome", "interval_start", "interval_end", "Belem_RAiSD_u")
)

values <- extract_control_values_from_intervals(
  belem_raisd,
  "Belem_RAiSD_u",
  control_interval_dt
)
belem_raisd_threshold <- as.numeric(quantile(values, probs = 0.999, na.rm = TRUE))
selection_window_significance_list[["Belem_RAiSD_u"]] <-
  make_outlier_window_selection_flags_from_intervals(
    belem_raisd,
    "Belem_RAiSD_u",
    outliers,
    threshold = belem_raisd_threshold,
    use_abs = FALSE
  )
add_additional_thresholds(
  values,
  "Belem_RAiSD_u",
  "Belem_RAiSD_control_windows",
  "upper"
)
rm(belem_raisd, values, belem_raisd_threshold)
gc()

# -----------------------------
# selscan iHS and nSL statistics
# -----------------------------
xingu_ihs <- fread(
  xingu_ihs_file,
  select = c("chr", "pos", "norm_ihs")
)
setnames(
  xingu_ihs,
  old = c("chr", "norm_ihs"),
  new = c("chromosome", "xingu_norm_ihs")
)
values <- extract_control_values_from_points(
  xingu_ihs,
  "xingu_norm_ihs",
  control_interval_dt
)
selection_window_significance_list[["xingu_norm_ihs"]] <-
  make_outlier_window_selection_flags_from_points(
    xingu_ihs,
    "xingu_norm_ihs",
    outliers,
    threshold = 2,
    use_abs = TRUE
  )
add_additional_thresholds(
  values,
  "xingu_norm_ihs",
  "Xingu_iHS_control_windows",
  c("lower", "upper")
)
rm(xingu_ihs, values)
gc()

xingu_nsl <- fread(
  xingu_nsl_file,
  select = c("chr", "pos", "norm_nsl")
)
setnames(
  xingu_nsl,
  old = c("chr", "norm_nsl"),
  new = c("chromosome", "xingu_norm_nsl")
)
values <- extract_control_values_from_points(
  xingu_nsl,
  "xingu_norm_nsl",
  control_interval_dt
)
selection_window_significance_list[["xingu_norm_nsl"]] <-
  make_outlier_window_selection_flags_from_points(
    xingu_nsl,
    "xingu_norm_nsl",
    outliers,
    threshold = 2,
    use_abs = TRUE
  )
add_additional_thresholds(
  values,
  "xingu_norm_nsl",
  "Xingu_nSL_control_windows",
  c("lower", "upper")
)
rm(xingu_nsl, values)
gc()

belem_ihs <- fread(
  belem_ihs_file,
  select = c("chr", "pos", "norm_ihs")
)
setnames(
  belem_ihs,
  old = c("chr", "norm_ihs"),
  new = c("chromosome", "belem_norm_ihs")
)
values <- extract_control_values_from_points(
  belem_ihs,
  "belem_norm_ihs",
  control_interval_dt
)
selection_window_significance_list[["belem_norm_ihs"]] <-
  make_outlier_window_selection_flags_from_points(
    belem_ihs,
    "belem_norm_ihs",
    outliers,
    threshold = 2,
    use_abs = TRUE
  )
add_additional_thresholds(
  values,
  "belem_norm_ihs",
  "Belem_iHS_control_windows",
  c("lower", "upper")
)
rm(belem_ihs, values)
gc()

belem_nsl <- fread(
  belem_nsl_file,
  select = c("chr", "pos", "norm_nsl")
)
setnames(
  belem_nsl,
  old = c("chr", "norm_nsl"),
  new = c("chromosome", "belem_norm_nsl")
)
values <- extract_control_values_from_points(
  belem_nsl,
  "belem_norm_nsl",
  control_interval_dt
)
selection_window_significance_list[["belem_norm_nsl"]] <-
  make_outlier_window_selection_flags_from_points(
    belem_nsl,
    "belem_norm_nsl",
    outliers,
    threshold = 2,
    use_abs = TRUE
  )
add_additional_thresholds(
  values,
  "belem_norm_nsl",
  "Belem_nSL_control_windows",
  c("lower", "upper")
)
rm(belem_nsl, values)
gc()

# -----------------------------
# rCNV statistics
# -----------------------------
rcnv <- fread(
  rcnv_file,
  select = c("scaffold", "start", "end", "total.sites", "dup.sites", "perc.dup.sites"),
  na.strings = c("NA", "")
)
setnames(
  rcnv,
  old = c("scaffold", "start", "end", "total.sites", "dup.sites", "perc.dup.sites"),
  new = c("chromosome", "interval_start", "interval_end", "total_sites", "dup_sites", "perc_dup_sites")
)
rcnv[, chromosome := standardize_scaffold_names(chromosome)]
rcnv <- rcnv[!is.na(total_sites) & total_sites >= 50]
values <- extract_control_values_from_intervals(
  rcnv[, .(chromosome, interval_start, interval_end, dup_sites)],
  "dup_sites",
  control_interval_dt
)
add_additional_thresholds(
  values,
  "dup_sites",
  "rCNV_control_windows",
  "upper"
)
rm(values)

values <- extract_control_values_from_intervals(
  rcnv[, .(chromosome, interval_start, interval_end, perc_dup_sites)],
  "perc_dup_sites",
  control_interval_dt
)
add_additional_thresholds(
  values,
  "perc_dup_sites",
  "rCNV_control_windows",
  "upper"
)
rm(rcnv, values)
gc()

# -----------------------------
# ARG-based branch Fst statistics
# -----------------------------
arg_based_fst <- fread(
  arg_based_fst_file,
  select = c(
    "region",
    "coordinate",
    "Fst_belem_tapajos",
    "Fst_belem_xingu",
    "Fst_tapajos_xingu"
  )
)

setnames(
  arg_based_fst,
  old = c(
    "region",
    "coordinate",
    "Fst_tapajos_xingu",
    "Fst_belem_tapajos",
    "Fst_belem_xingu"
  ),
  new = c(
    "chromosome",
    "pos",
    "Tap_Xin_Fst_ARG_based",
    "Tap_Bel_Fst_ARG_based",
    "Xin_Bel_Fst_ARG_based"
  )
)
arg_based_fst[, chromosome := standardize_scaffold_names(chromosome)]
arg_based_fst <- filter_points_outside_mask(arg_based_fst, argweaver_mask, pos_col = "pos")

for (stat in c(
  "Tap_Xin_Fst_ARG_based",
  "Tap_Bel_Fst_ARG_based",
  "Xin_Bel_Fst_ARG_based"
)) {
  values <- extract_control_values_from_points(
    arg_based_fst[, c("chromosome", "pos", stat), with = FALSE],
    stat,
    control_interval_dt
  )
  add_additional_thresholds(
    values,
    stat,
    "ARG_based_branch_Fst_control_windows",
    "upper"
  )
  rm(values)
}
rm(arg_based_fst)
gc()

additional_threshold_summary <- rbindlist(additional_threshold_list, fill = TRUE)
additional_control_value_summary <- rbindlist(additional_value_summary_list, fill = TRUE)

# selscan recommends fixed thresholds of -2 and 2 for normalized iHS/nSL.
# Remove the empirical percentile thresholds for these statistics and replace
# them with fixed lower/upper thresholds.
selscan_threshold_stats <- c(
  "xingu_norm_ihs",
  "xingu_norm_nsl",
  "belem_norm_ihs",
  "belem_norm_nsl"
)

additional_threshold_summary <- additional_threshold_summary[
  !statistic %in% selscan_threshold_stats
]

additional_threshold_summary <- rbind(
  additional_threshold_summary,
  make_fixed_threshold_rows(
    selscan_threshold_stats,
    "selscan_recommended_fixed_threshold",
    additional_control_value_summary
  ),
  fill = TRUE
)

fwrite(
  additional_threshold_summary,
  additional_threshold_summary_file,
  sep = "\t"
)

fwrite(
  additional_control_value_summary,
  additional_control_value_summary_file,
  sep = "\t"
)

# -----------------------------
# Identify peaks:
# outlier windows on same chromosome within 50 kb
# -----------------------------
setorder(outliers, chromosome, start, end)

outliers[, gap_from_previous := start - shift(end), by = chromosome]
outliers[, new_peak := is.na(gap_from_previous) | gap_from_previous > 50000, by = chromosome]
outliers[, peak_id := cumsum(new_peak)]

peaks <- outliers[, .(
  start = min(start),
  end   = max(end),
  n_outlier_windows = .N,
  max_Xin_Bel_Fst = max(Xin_Bel_Fst, na.rm = TRUE),
  mean_Xin_Bel_Fst = mean(Xin_Bel_Fst, na.rm = TRUE)
), by = .(chromosome, peak_id)]

setorder(peaks, chromosome, start)

peak_bed <- peaks[, .(
  chromosome,
  start,
  end,
  peak_id,
  n_outlier_windows,
  max_Xin_Bel_Fst,
  mean_Xin_Bel_Fst
)]

fwrite(peak_bed, outlier_peaks_bed, sep = "\t", col.names = TRUE)

large_repeat_peak_proximity <- make_large_repeat_peak_proximity_table(
  peaks_dt = peaks,
  repeat_regions = repeat_regions_for_manhattan,
  outfile = large_repeat_peak_proximity_file,
  min_repeat_size_bp = 50000L,
  max_distance_bp = 50000L,
  max_repeat_gap_bp = 1000L
)

# Map each outlier Fst window to its merged Xin-Bel Fst peak.
# This allows ARG local trees assigned to outlier windows to be pooled by peak.
peak_window_lookup <- outliers[, .(
  chromosome,
  window_start = start,
  window_end = end,
  peak_id
)]

selection_peak_significance_flags <- make_peak_selection_significance_flags(
  outlier_windows_with_peak_dt = outliers,
  selection_window_flag_list = selection_window_significance_list
)

rm(df, outlier_bed, peak_bed)
gc()

# =============================================================
# ARGweaver section: iteration 2000
# =============================================================

# -----------------------------
# Read only needed ARGweaver columns
# -----------------------------
arg_cols_needed <- c(
  "chrom",
  "pos",
  "TMRCAH_all",
  "Tap_Xin_TMRCAH",
  "Tap_Xin_RTH_original",
  "Tap_Xin_RTH_prime",
  "belem_RTH",
  "tapajos_RTH",
  "xingu_RTH",
  "belem_enrich",
  "tapajos_enrich",
  "xingu_enrich"
)

arg <- fread(
  arg_stat_file_2000,
  select = arg_cols_needed
)

# -----------------------------
# Standardize ARG chromosome names
#
# scaffold100-553 -> scaffold_100
# scaffold17-202  -> scaffold_17
# -----------------------------
arg[, chromosome := standardize_scaffold_names(chrom)]
arg[, chrom := NULL]
arg <- filter_points_outside_mask(arg, argweaver_mask, pos_col = "pos")

# -----------------------------
# Add inverse RTH statistics
# -----------------------------
arg[, tapajos_RTH_inverse := safe_inverse(tapajos_RTH)]
arg[, xingu_RTH_inverse := safe_inverse(xingu_RTH)]
arg[, belem_RTH_inverse := safe_inverse(belem_RTH)]
arg[, Tap_Xin_RTH_prime_inverse := safe_inverse(Tap_Xin_RTH_prime)]

# -----------------------------
# Assign ARG local trees to outlier/control pixy windows
# using BED-style intervals: start <= pos < end
# -----------------------------
arg[, pos_start := pos]
arg[, pos_end := pos]
arg[, pos := NULL]

fd_rth_correlation_summary <- make_fd_rth_correlation_outputs(
  arg_points_dt = arg,
  d_stats_file = d_stats_file,
  all_windows_dt = df_clean,
  outfile_pdf = fd_rth_correlation_pdf,
  outfile_data = fd_rth_correlation_data_file,
  outfile_summary = fd_rth_correlation_summary_file
)

gc()

arg_windows <- comparison_windows[, .(
  chromosome,
  window_start = start,
  window_end = end,
  window_end_inclusive = end - 1,
  window_id,
  window_class
)]

setkey(arg, chromosome, pos_start, pos_end)
setkey(arg_windows, chromosome, window_start, window_end_inclusive)

arg_joined <- foverlaps(
  arg,
  arg_windows,
  by.x = c("chromosome", "pos_start", "pos_end"),
  by.y = c("chromosome", "window_start", "window_end_inclusive"),
  nomatch = 0
)

rm(arg)
gc()

arg_joined[, pos := pos_start]
arg_joined[, c("pos_start", "pos_end", "window_end_inclusive") := NULL]

arg_joined_peaks <- merge(
  arg_joined[window_class == "Xin_Bel_Fst_outlier"],
  peak_window_lookup,
  by = c("chromosome", "window_start", "window_end"),
  all.x = FALSE,
  all.y = FALSE
)

# -----------------------------
# Tap_Xin_RTH_original and Tap_Xin_RTH_prime diagnostics
# -----------------------------
tap_xin_rth_stats <- c(
  "Tap_Xin_RTH_original",
  "Tap_Xin_RTH_prime"
)

tap_xin_rth_value_one_summary <- rbindlist(lapply(
  tap_xin_rth_stats,
  function(stat) {
    rbindlist(lapply(unique(arg_joined$window_class), function(class_name) {
      tmp <- arg_joined[window_class == class_name]
      values <- tmp[[stat]]
      non_na_values <- values[!is.na(values)]

      data.table(
        statistic = stat,
        window_class = class_name,
        n_local_trees = length(non_na_values),
        n_local_trees_equal_1 = sum(non_na_values == 1),
        prop_local_trees_equal_1 = mean(non_na_values == 1),
        percent_local_trees_equal_1 = 100 * mean(non_na_values == 1),
        n_windows = uniqueN(tmp[!is.na(get(stat)), window_id]),
        n_windows_with_any_equal_1 = uniqueN(tmp[get(stat) == 1, window_id])
      )
    }))
  }
))

fwrite(
  tap_xin_rth_value_one_summary,
  tap_xin_rth_value_one_summary_file,
  sep = "\t"
)

tap_xin_rth_decile_distribution <- rbindlist(lapply(
  tap_xin_rth_stats,
  function(stat) {
    rbindlist(lapply(unique(arg_joined$window_class), function(class_name) {
      tmp <- arg_joined[
        window_class == class_name &
          !is.na(get(stat))
      ]

      # RTH values are expected to be between 0 and 1.
      # Values are assigned to 10 fixed bins:
      # [0.0,0.1), [0.1,0.2), ..., [0.9,1.0].
      values <- tmp[[stat]]
      bin_index <- pmin(pmax(floor(values * 10) + 1L, 1L), 10L)

      out <- data.table(
        bin_index = bin_index
      )[, .(
        n_local_trees = .N
      ), by = bin_index]

      all_bins <- data.table(bin_index = 1:10)
      out <- merge(all_bins, out, by = "bin_index", all.x = TRUE)
      out[is.na(n_local_trees), n_local_trees := 0L]

      out[, `:=`(
        statistic = stat,
        window_class = class_name,
        bin_lower = (bin_index - 1) / 10,
        bin_upper = bin_index / 10,
        total_local_trees = sum(n_local_trees)
      )]

      out[, prop_local_trees := n_local_trees / total_local_trees]
      out[, percent_local_trees := 100 * prop_local_trees]

      out[, bin_label := paste0(
        sprintf("%.1f", bin_lower),
        "-",
        sprintf("%.1f", bin_upper)
      )]

      setcolorder(out, c(
        "statistic",
        "window_class",
        "bin_index",
        "bin_label",
        "bin_lower",
        "bin_upper",
        "n_local_trees",
        "total_local_trees",
        "prop_local_trees",
        "percent_local_trees"
      ))

      out[]
    }))
  }
))

fwrite(
  tap_xin_rth_decile_distribution,
  tap_xin_rth_decile_distribution_file,
  sep = "\t"
)

make_rth_decile_distribution_barplot(
  tap_xin_rth_decile_distribution,
  tap_xin_rth_decile_barplot_pdf
)

# -----------------------------
# ARG empirical significance thresholds
# based on control windows/local trees
# -----------------------------
arg_control_trees <- arg_joined[window_class == "control"]

arg_high_tail_tree_stats <- c(
  "tapajos_enrich",
  "xingu_enrich",
  "belem_enrich",
  "tapajos_RTH_inverse",
  "xingu_RTH_inverse",
  "belem_RTH_inverse",
  "Tap_Xin_RTH_prime_inverse"
)

arg_low_tail_tree_stats <- c(
  "xingu_RTH",
  "belem_RTH",
  "Tap_Xin_RTH_original"
)

arg_threshold_high_tree_summary <- rbindlist(lapply(arg_high_tail_tree_stats, function(stat) {
  data.table(
    statistic = stat,
    threshold_basis = "ARG_local_trees_in_control_windows",
    tail = "upper",
    p_value = 0.001,
    percentile = 99.9,
    threshold = as.numeric(quantile(
      arg_control_trees[[stat]],
      probs = 0.999,
      na.rm = TRUE
    )),
    n_values = sum(!is.na(arg_control_trees[[stat]]))
  )
}))

arg_threshold_low_tree_summary <- rbindlist(lapply(arg_low_tail_tree_stats, function(stat) {
  data.table(
    statistic = stat,
    threshold_basis = "ARG_local_trees_in_control_windows",
    tail = "lower",
    p_value = 0.001,
    percentile = 0.1,
    threshold = as.numeric(quantile(
      arg_control_trees[[stat]],
      probs = 0.001,
      na.rm = TRUE
    )),
    n_values = sum(!is.na(arg_control_trees[[stat]]))
  )
}))

arg_threshold_summary <- rbind(
  arg_threshold_high_tree_summary,
  arg_threshold_low_tree_summary,
  fill = TRUE
)

fwrite(arg_threshold_summary, arg_threshold_summary_file, sep = "\t")

rm(arg_control_trees)
gc()

# =============================================================
# ARGweaver cross-coalescence section: iteration 2000
# =============================================================

# -----------------------------
# Read only needed ARGweaver CC columns
# -----------------------------
arg_cc_cols_needed <- c(
  "chrom",
  "pos",
  "Xin_Bel_CC_original",
  "Xin_Bel_CC_prime",
  "Tap_Xin_JCR_v2"
)

arg_cc <- fread(
  arg_cc_stat_file_2000,
  select = arg_cc_cols_needed
)

# -----------------------------
# Standardize ARG CC chromosome names
#
# scaffold100-553 -> scaffold_100
# scaffold17-202  -> scaffold_17
# -----------------------------
arg_cc[, chromosome := standardize_scaffold_names(chrom)]
arg_cc[, chrom := NULL]
arg_cc <- filter_points_outside_mask(arg_cc, argweaver_mask, pos_col = "pos")

# -----------------------------
# Assign ARG CC local trees to outlier/control pixy windows
# using BED-style intervals: start <= pos < end
# -----------------------------
arg_cc[, pos_start := pos]
arg_cc[, pos_end := pos]
arg_cc[, pos := NULL]

gc()

arg_cc_windows <- comparison_windows[, .(
  chromosome,
  window_start = start,
  window_end = end,
  window_end_inclusive = end - 1,
  window_id,
  window_class
)]

setkey(arg_cc, chromosome, pos_start, pos_end)
setkey(arg_cc_windows, chromosome, window_start, window_end_inclusive)

arg_cc_joined <- foverlaps(
  arg_cc,
  arg_cc_windows,
  by.x = c("chromosome", "pos_start", "pos_end"),
  by.y = c("chromosome", "window_start", "window_end_inclusive"),
  nomatch = 0
)

rm(arg_cc)
gc()

arg_cc_joined[, pos := pos_start]
arg_cc_joined[, c("pos_start", "pos_end", "window_end_inclusive") := NULL]

arg_cc_joined_peaks <- merge(
  arg_cc_joined[window_class == "Xin_Bel_Fst_outlier"],
  peak_window_lookup,
  by = c("chromosome", "window_start", "window_end"),
  all.x = FALSE,
  all.y = FALSE
)

# =============================================================
# Heat map: representative percentile scores for Xin-Bel Fst peaks
# =============================================================
# This section creates 10-kb window-level averages for all requested
# statistics using the pixy windows as the common coordinate reference.
# It then converts each statistic to genome-wide percentile scores and
# represents each Fst peak by the outlier window whose percentile score is
# farthest from 0.5 for that statistic.

heatmap_window_values <- df_clean[, .(
  chromosome,
  window_start = start,
  window_end = end,
  window_id = paste(chromosome, start, end, sep = ":"),
  Tap_Xin_Fst,
  Tap_Bel_Fst,
  Xin_Bel_Fst,
  Tap_Xin_Dxy,
  Tap_Bel_Dxy,
  Xin_Bel_Dxy,
  Tap_pi,
  Xin_pi,
  Bel_pi
)]

heatmap_window_intervals <- heatmap_window_values[, .(
  chromosome,
  window_start = as.integer(window_start),
  window_end_inclusive = as.integer(window_end - 1L),
  window_id
)]

# D statistics: already 10-kb but 1-based; convert to BED-style.
d_stats_heat <- fread(
  d_stats_file,
  select = c("scaffold", "start", "end", "sitesUsed", "D", "fd")
)
setnames(
  d_stats_heat,
  old = c("scaffold", "fd"),
  new = c("chromosome", "Tapajos_Xingu_fd")
)
# fd is not meaningful when D is negative, so set fd to zero for those windows.
d_stats_heat[!is.na(D) & D < 0, Tapajos_Xingu_fd := 0]
# Exclude biologically invalid fd values before heat-map window averaging.
d_stats_heat <- d_stats_heat[
  is.na(Tapajos_Xingu_fd) |
    (Tapajos_Xingu_fd >= 0 & Tapajos_Xingu_fd <= 1)
]
d_stats_heat <- d_stats_heat[sitesUsed >= 100]
d_stats_heat[, interval_start := as.integer(start - 1L)]
d_stats_heat[, interval_end := as.integer(end)]
heatmap_window_values <- add_window_stat_column(
  heatmap_window_values,
  "Tapajos_Xingu_fd",
  mean_by_window_from_intervals(
    d_stats_heat[, .(chromosome, interval_start, interval_end, Tapajos_Xingu_fd)],
    "Tapajos_Xingu_fd",
    heatmap_window_intervals
  )
)
rm(d_stats_heat)
gc()

# ReLERNN recombination-rate intervals.
relernn_heat <- fread(
  relernn_file,
  header = FALSE,
  col.names = c("chromosome", "interval_start", "interval_end", "recombination_rate")
)
heatmap_window_values <- add_window_stat_column(
  heatmap_window_values,
  "recombination_rate",
  mean_by_window_from_intervals(relernn_heat, "recombination_rate", heatmap_window_intervals)
)
rm(relernn_heat)
gc()

# RAiSD point/short-interval statistics.
xingu_raisd_heat <- fread(
  xingu_raisd_file,
  header = FALSE,
  col.names = c("chromosome", "interval_start", "interval_end", "Xingu_RAiSD_u")
)
heatmap_window_values <- add_window_stat_column(
  heatmap_window_values,
  "Xingu_RAiSD_u",
  mean_by_window_from_intervals(xingu_raisd_heat, "Xingu_RAiSD_u", heatmap_window_intervals)
)
rm(xingu_raisd_heat)
gc()

belem_raisd_heat <- fread(
  belem_raisd_file,
  header = FALSE,
  col.names = c("chromosome", "interval_start", "interval_end", "Belem_RAiSD_u")
)
heatmap_window_values <- add_window_stat_column(
  heatmap_window_values,
  "Belem_RAiSD_u",
  mean_by_window_from_intervals(belem_raisd_heat, "Belem_RAiSD_u", heatmap_window_intervals)
)
rm(belem_raisd_heat)
gc()

# rCNV 10-kb windows.
rcnv_heat <- fread(
  rcnv_file,
  select = c("scaffold", "start", "end", "total.sites", "dup.sites", "perc.dup.sites"),
  na.strings = c("NA", "")
)
setnames(
  rcnv_heat,
  old = c("scaffold", "start", "end", "total.sites", "dup.sites", "perc.dup.sites"),
  new = c("chromosome", "interval_start", "interval_end", "total_sites", "dup_sites", "perc_dup_sites")
)
rcnv_heat[, chromosome := standardize_scaffold_names(chromosome)]
rcnv_heat <- rcnv_heat[!is.na(total_sites) & total_sites >= 50]
heatmap_window_values <- add_window_stat_column(
  heatmap_window_values,
  "dup_sites",
  mean_by_window_from_intervals(
    rcnv_heat[, .(chromosome, interval_start, interval_end, dup_sites)],
    "dup_sites",
    heatmap_window_intervals
  )
)
heatmap_window_values <- add_window_stat_column(
  heatmap_window_values,
  "perc_dup_sites",
  mean_by_window_from_intervals(
    rcnv_heat[, .(chromosome, interval_start, interval_end, perc_dup_sites)],
    "perc_dup_sites",
    heatmap_window_intervals
  )
)
rm(rcnv_heat)
gc()

# ARG-based branch Fst point statistics.
arg_based_fst_heat <- fread(
  arg_based_fst_file,
  select = c(
    "region",
    "coordinate",
    "Fst_belem_tapajos",
    "Fst_belem_xingu",
    "Fst_tapajos_xingu"
  )
)
setnames(
  arg_based_fst_heat,
  old = c(
    "region",
    "coordinate",
    "Fst_tapajos_xingu",
    "Fst_belem_tapajos",
    "Fst_belem_xingu"
  ),
  new = c(
    "chromosome",
    "pos",
    "Tap_Xin_Fst_ARG_based",
    "Tap_Bel_Fst_ARG_based",
    "Xin_Bel_Fst_ARG_based"
  )
)
arg_based_fst_heat[, chromosome := standardize_scaffold_names(chromosome)]
arg_based_fst_heat <- filter_points_outside_mask(arg_based_fst_heat, argweaver_mask, pos_col = "pos")
for (stat in c("Tap_Xin_Fst_ARG_based", "Tap_Bel_Fst_ARG_based", "Xin_Bel_Fst_ARG_based")) {
  heatmap_window_values <- add_window_stat_column(
    heatmap_window_values,
    stat,
    mean_by_window_from_points(
      arg_based_fst_heat[, c("chromosome", "pos", stat), with = FALSE],
      stat,
      heatmap_window_intervals
    )
  )
}
rm(arg_based_fst_heat)
gc()

# ARG local-tree point statistics.
arg_heat <- fread(
  arg_stat_file_2000,
  select = c(
    "chrom",
    "pos",
    "Tap_Xin_RTH_original",
    "belem_RTH",
    "tapajos_RTH",
    "xingu_RTH",
    "belem_enrich",
    "tapajos_enrich",
    "xingu_enrich"
  )
)
arg_heat[, chromosome := standardize_scaffold_names(chrom)]
arg_heat[, chrom := NULL]
arg_heat <- filter_points_outside_mask(arg_heat, argweaver_mask, pos_col = "pos")
for (stat in c(
  "tapajos_enrich",
  "xingu_enrich",
  "belem_enrich",
  "tapajos_RTH",
  "xingu_RTH",
  "belem_RTH",
  "Tap_Xin_RTH_original"
)) {
  heatmap_window_values <- add_window_stat_column(
    heatmap_window_values,
    stat,
    mean_by_window_from_points(
      arg_heat[, c("chromosome", "pos", stat), with = FALSE],
      stat,
      heatmap_window_intervals
    )
  )
}
rm(arg_heat)
gc()

# ARG cross-coalescence point statistic.
arg_cc_heat <- fread(
  arg_cc_stat_file_2000,
  select = c("chrom", "pos", "Xin_Bel_CC_original")
)
arg_cc_heat[, chromosome := standardize_scaffold_names(chrom)]
arg_cc_heat[, chrom := NULL]
arg_cc_heat <- filter_points_outside_mask(arg_cc_heat, argweaver_mask, pos_col = "pos")
heatmap_window_values <- add_window_stat_column(
  heatmap_window_values,
  "Xin_Bel_CC_original",
  mean_by_window_from_points(arg_cc_heat, "Xin_Bel_CC_original", heatmap_window_intervals)
)
rm(arg_cc_heat)
gc()

# selscan point statistics. For heat map purposes only, use absolute values
# before averaging within 10-kb windows.
xingu_ihs_heat <- fread(xingu_ihs_file, select = c("chr", "pos", "norm_ihs"))
setnames(xingu_ihs_heat, c("chr", "norm_ihs"), c("chromosome", "xingu_norm_ihs"))
heatmap_window_values <- add_window_stat_column(
  heatmap_window_values,
  "xingu_norm_ihs",
  mean_by_window_from_points(xingu_ihs_heat, "xingu_norm_ihs", heatmap_window_intervals, use_abs = TRUE)
)
rm(xingu_ihs_heat)
gc()

xingu_nsl_heat <- fread(xingu_nsl_file, select = c("chr", "pos", "norm_nsl"))
setnames(xingu_nsl_heat, c("chr", "norm_nsl"), c("chromosome", "xingu_norm_nsl"))
heatmap_window_values <- add_window_stat_column(
  heatmap_window_values,
  "xingu_norm_nsl",
  mean_by_window_from_points(xingu_nsl_heat, "xingu_norm_nsl", heatmap_window_intervals, use_abs = TRUE)
)
rm(xingu_nsl_heat)
gc()

belem_ihs_heat <- fread(belem_ihs_file, select = c("chr", "pos", "norm_ihs"))
setnames(belem_ihs_heat, c("chr", "norm_ihs"), c("chromosome", "belem_norm_ihs"))
heatmap_window_values <- add_window_stat_column(
  heatmap_window_values,
  "belem_norm_ihs",
  mean_by_window_from_points(belem_ihs_heat, "belem_norm_ihs", heatmap_window_intervals, use_abs = TRUE)
)
rm(belem_ihs_heat)
gc()

belem_nsl_heat <- fread(belem_nsl_file, select = c("chr", "pos", "norm_nsl"))
setnames(belem_nsl_heat, c("chr", "norm_nsl"), c("chromosome", "belem_norm_nsl"))
heatmap_window_values <- add_window_stat_column(
  heatmap_window_values,
  "belem_norm_nsl",
  mean_by_window_from_points(belem_nsl_heat, "belem_norm_nsl", heatmap_window_intervals, use_abs = TRUE)
)
rm(belem_nsl_heat)
gc()

heatmap_stat_order <- c(
  "Tap_Xin_Fst",
  "Tap_Xin_Fst_ARG_based",
  "Tap_Bel_Fst",
  "Tap_Bel_Fst_ARG_based",
  "Xin_Bel_Fst",
  "Xin_Bel_Fst_ARG_based",
  "Tap_Xin_Dxy",
  "Tap_Bel_Dxy",
  "Xin_Bel_Dxy",
  "Tap_pi",
  "Xin_pi",
  "Bel_pi",
  "tapajos_enrich",
  "xingu_enrich",
  "belem_enrich",
  "tapajos_RTH",
  "xingu_RTH",
  "belem_RTH",
  "Xin_Bel_CC_original",
  "Tap_Xin_RTH_original",
  "Tapajos_Xingu_fd",
  "recombination_rate",
  "Xingu_RAiSD_u",
  "xingu_norm_ihs",
  "xingu_norm_nsl",
  "Belem_RAiSD_u",
  "belem_norm_ihs",
  "belem_norm_nsl",
  "dup_sites",
  "perc_dup_sites"
)

fwrite(heatmap_window_values, heatmap_window_values_file, sep = "\t")

heatmap_window_percentiles <- copy(heatmap_window_values)
for (stat in heatmap_stat_order) {
  heatmap_window_percentiles[, (stat) := percentile_rank_0_1(get(stat))]
}

fwrite(heatmap_window_percentiles, heatmap_window_percentiles_file, sep = "\t")

heatmap_peak_window_percentiles <- merge(
  peak_window_lookup,
  heatmap_window_percentiles,
  by = c("chromosome", "window_start", "window_end"),
  all.x = TRUE,
  sort = FALSE
)

heatmap_arg_stats <- c(
  "Tap_Xin_Fst_ARG_based",
  "Tap_Bel_Fst_ARG_based",
  "Xin_Bel_Fst_ARG_based",
  "tapajos_enrich",
  "xingu_enrich",
  "belem_enrich",
  "tapajos_RTH",
  "xingu_RTH",
  "belem_RTH",
  "Xin_Bel_CC_original",
  "Tap_Xin_RTH_original"
)

peaks_with_ARG_heatmap_data <- heatmap_peak_window_percentiles[, .(
  has_ARG_heatmap_data = any(!is.na(unlist(.SD)))
), by = peak_id, .SDcols = heatmap_arg_stats][
  has_ARG_heatmap_data == TRUE, peak_id
]

heatmap_peak_representative <- peaks[
  peak_id %in% peaks_with_ARG_heatmap_data,
  .(peak_id, chromosome, start, end, n_outlier_windows)
]

for (stat in heatmap_stat_order) {
  representative_values <- heatmap_peak_window_percentiles[
    peak_id %in% peaks_with_ARG_heatmap_data & !is.na(get(stat)),
    .SD[which.max(abs(get(stat) - 0.5))][1],
    by = peak_id,
    .SDcols = stat
  ]
  setnames(representative_values, stat, paste0(stat, "_representative_percentile"))
  heatmap_peak_representative <- merge(
    heatmap_peak_representative,
    representative_values,
    by = "peak_id",
    all.x = TRUE,
    sort = FALSE
  )
}

# Use clean statistic names as heat-map matrix columns.
for (stat in heatmap_stat_order) {
  old_name <- paste0(stat, "_representative_percentile")
  if (old_name %in% names(heatmap_peak_representative)) {
    setnames(heatmap_peak_representative, old_name, stat)
  }
}

setorder(heatmap_peak_representative, chromosome, start, end)

fwrite(heatmap_peak_representative, heatmap_peak_representative_file, sep = "\t")

make_peak_heatmap_pdf(
  heatmap_peak_representative,
  heatmap_stat_order,
  plot_labels,
  heatmap_pdf
)

rm(
  heatmap_window_values,
  heatmap_window_intervals,
  heatmap_window_percentiles,
  heatmap_peak_window_percentiles
)
gc()

# -----------------------------
# ARG CC empirical significance thresholds
# based on local trees in control windows
# -----------------------------
arg_cc_control_trees <- arg_cc_joined[window_class == "control"]

arg_cc_high_tail_tree_stats <- c(
  "Xin_Bel_CC_original",
  "Xin_Bel_CC_prime",
  "Tap_Xin_JCR_v2"
)

arg_cc_threshold_summary <- rbindlist(lapply(arg_cc_high_tail_tree_stats, function(stat) {
  data.table(
    statistic = stat,
    threshold_basis = "ARG_CC_local_trees_in_control_windows",
    tail = "upper",
    p_value = 0.001,
    percentile = 99.9,
    threshold = as.numeric(quantile(
      arg_cc_control_trees[[stat]],
      probs = 0.999,
      na.rm = TRUE
    )),
    n_values = sum(!is.na(arg_cc_control_trees[[stat]]))
  )
}))

fwrite(arg_cc_threshold_summary, arg_cc_threshold_summary_file, sep = "\t")

# -----------------------------
# Summarize CC statistics in outlier and control windows
# -----------------------------
arg_cc_stats_to_summarize <- c(
  "Xin_Bel_CC_original",
  "Xin_Bel_CC_prime",
  "Tap_Xin_JCR_v2"
)

arg_cc_stat_summary <- rbind(
  summarize_stats(
    arg_cc_joined[window_class == "Xin_Bel_Fst_outlier"],
    "Xin_Bel_Fst_outlier",
    arg_cc_stats_to_summarize
  ),
  summarize_stats(
    arg_cc_joined[window_class == "control"],
    "control",
    arg_cc_stats_to_summarize
  )
)

fwrite(arg_cc_stat_summary, arg_cc_stat_summary_file, sep = "\t")

# -----------------------------
# Compare ceiling effects:
# Xin_Bel_CC_original vs Xin_Bel_CC_prime
# -----------------------------
arg_cc_value_one_summary <- rbindlist(lapply(
  c("Xin_Bel_CC_original", "Xin_Bel_CC_prime"),
  function(stat) {
    rbindlist(lapply(unique(arg_cc_joined$window_class), function(class_name) {
      tmp <- arg_cc_joined[window_class == class_name]

      data.table(
        statistic = stat,
        window_class = class_name,
        n_local_trees = sum(!is.na(tmp[[stat]])),
        n_local_trees_equal_1 = sum(tmp[[stat]] == 1, na.rm = TRUE),
        prop_local_trees_equal_1 = mean(tmp[[stat]] == 1, na.rm = TRUE),
        percent_local_trees_equal_1 = 100 * mean(tmp[[stat]] == 1, na.rm = TRUE),
        n_windows = uniqueN(tmp$window_id),
        n_windows_with_any_equal_1 = uniqueN(tmp[get(stat) == 1, window_id])
      )
    }))
  }
))

fwrite(arg_cc_value_one_summary, arg_cc_value_one_summary_file, sep = "\t")

# -----------------------------
# Scatterplot:
# Tap_Xin_JCR_v2 vs Xin_Bel_CC_original
# Control windows only, one ARG CC local tree sampled per 10-kb window
# -----------------------------
set.seed(123)

arg_cc_control_sample <- arg_cc_joined[
  window_class == "control" &
    !is.na(Tap_Xin_JCR_v2) &
    !is.na(Xin_Bel_CC_original),
  .SD[sample(.N, 1)],
  by = window_id
]

if (nrow(arg_cc_control_sample) >= 2) {
  jcr_cc_model <- lm(
    Xin_Bel_CC_original ~ Tap_Xin_JCR_v2,
    data = arg_cc_control_sample
  )

  jcr_cc_r2 <- summary(jcr_cc_model)$r.squared

  pdf(arg_cc_jcr_rcc_scatter_pdf, width = 6, height = 6)

  plot(
    arg_cc_control_sample$Tap_Xin_JCR_v2,
    arg_cc_control_sample$Xin_Bel_CC_original,
    pch = 16,
    cex = 0.55,
    xlab = "Tapajos-Xingu JCR v2",
    ylab = "Xingu-Belem RCC",
    main = paste0(
      "Control windows only\n",
      "One local tree sampled per 10-kb window; R^2 = ",
      round(jcr_cc_r2, 4)
    )
  )

  abline(jcr_cc_model, lwd = 2)

  dev.off()
} else {
  jcr_cc_r2 <- NA_real_

  pdf(arg_cc_jcr_rcc_scatter_pdf, width = 6, height = 6)
  plot.new()
  text(
    0.5,
    0.5,
    "Not enough non-missing control-window data\nto plot Tap_Xin_JCR_v2 vs Xin_Bel_CC_original"
  )
  dev.off()
}

rm(arg_cc_control_trees, arg_cc_control_sample)
gc()

# -----------------------------
# Count Fst outlier/control windows containing
# any significant local tree for selected ARG statistics
# -----------------------------
barplot_stats <- data.table(
  statistic = c(
    "xingu_enrich",
    "belem_enrich",
    "xingu_RTH",
    "belem_RTH",
    "Xin_Bel_CC_original",
    "Tap_Xin_RTH_original"
  ),
  source = c(
    "ARG",
    "ARG",
    "ARG",
    "ARG",
    "ARG_CC",
    "ARG"
  ),
  tail = c(
    "upper",
    "upper",
    "lower",
    "lower",
    "upper",
    "lower"
  ),
  direction_label = c(
    "high",
    "high",
    "low",
    "low",
    "high",
    "low"
  ),
  fill_color = c(
    "#6ea8db",
    "#fdbb4a",
    "#6ea8db",
    "#fdbb4a",
    "#6ea8db",
    "#89ca8e"
  )
)

barplot_stats[, plot_label := paste0(
  unname(plot_labels[statistic]),
  " (",
  direction_label,
  ")"
)]

# -----------------------------
# Window- and peak-level significance status tables
# -----------------------------
outlier_window_reference <- outliers[, .(
  chromosome,
  window_start = start,
  window_end = end
)]
outlier_window_reference[, window_id := paste(chromosome, window_start, window_end, sep = ":")]

peak_reference <- peaks[, .(
  peak_id,
  chromosome,
  start,
  end,
  n_outlier_windows
)]

# -----------------------------
# ARG-based data coverage for Fst outlier windows and peaks
# -----------------------------
arg_window_ARG_counts <- arg_joined[
  window_class == "Xin_Bel_Fst_outlier",
  .(n_ARG_local_trees = .N),
  by = .(chromosome, window_start, window_end, window_id)
]

arg_cc_window_ARG_counts <- arg_cc_joined[
  window_class == "Xin_Bel_Fst_outlier",
  .(n_ARG_CC_local_trees = .N),
  by = .(chromosome, window_start, window_end, window_id)
]

arg_window_ARG_coverage <- merge(
  outlier_window_reference,
  arg_window_ARG_counts,
  by = c("chromosome", "window_start", "window_end", "window_id"),
  all.x = TRUE,
  sort = FALSE
)

arg_window_ARG_coverage <- merge(
  arg_window_ARG_coverage,
  arg_cc_window_ARG_counts,
  by = c("chromosome", "window_start", "window_end", "window_id"),
  all.x = TRUE,
  sort = FALSE
)

arg_window_ARG_coverage[is.na(n_ARG_local_trees), n_ARG_local_trees := 0L]
arg_window_ARG_coverage[is.na(n_ARG_CC_local_trees), n_ARG_CC_local_trees := 0L]
arg_window_ARG_coverage[, has_ARG_data := n_ARG_local_trees > 0]
arg_window_ARG_coverage[, has_ARG_CC_data := n_ARG_CC_local_trees > 0]
arg_window_ARG_coverage[, has_any_ARG_based_data := has_ARG_data | has_ARG_CC_data]
arg_window_ARG_coverage[, has_both_ARG_sources := has_ARG_data & has_ARG_CC_data]
arg_window_ARG_coverage[, n_total_ARG_based_local_trees :=
                          n_ARG_local_trees + n_ARG_CC_local_trees]

fwrite(
  arg_window_ARG_coverage,
  arg_window_ARG_coverage_file,
  sep = "\t"
)

arg_peak_ARG_counts <- arg_joined_peaks[, .(
  n_ARG_local_trees = .N
), by = peak_id]

arg_cc_peak_ARG_counts <- arg_cc_joined_peaks[, .(
  n_ARG_CC_local_trees = .N
), by = peak_id]

arg_peak_ARG_coverage <- merge(
  peak_reference,
  arg_peak_ARG_counts,
  by = "peak_id",
  all.x = TRUE,
  sort = FALSE
)

arg_peak_ARG_coverage <- merge(
  arg_peak_ARG_coverage,
  arg_cc_peak_ARG_counts,
  by = "peak_id",
  all.x = TRUE,
  sort = FALSE
)

arg_peak_ARG_coverage[is.na(n_ARG_local_trees), n_ARG_local_trees := 0L]
arg_peak_ARG_coverage[is.na(n_ARG_CC_local_trees), n_ARG_CC_local_trees := 0L]
arg_peak_ARG_coverage[, has_ARG_data := n_ARG_local_trees > 0]
arg_peak_ARG_coverage[, has_ARG_CC_data := n_ARG_CC_local_trees > 0]
arg_peak_ARG_coverage[, has_any_ARG_based_data := has_ARG_data | has_ARG_CC_data]
arg_peak_ARG_coverage[, has_both_ARG_sources := has_ARG_data & has_ARG_CC_data]
arg_peak_ARG_coverage[, n_total_ARG_based_local_trees :=
                        n_ARG_local_trees + n_ARG_CC_local_trees]

fwrite(
  arg_peak_ARG_coverage,
  arg_peak_ARG_coverage_file,
  sep = "\t"
)

# Add ARG-based-data coverage to the repeat-proximity supplementary table.
# The table retains all repeat-adjacent peaks; the one-off Manhattan figure
# below uses only rows with has_any_ARG_based_data == TRUE.
if (exists("large_repeat_peak_proximity") &&
    is.data.table(large_repeat_peak_proximity) &&
    nrow(large_repeat_peak_proximity) > 0) {
  repeat_arg_cols <- arg_peak_ARG_coverage[, .(
    peak_id,
    n_ARG_local_trees,
    n_ARG_CC_local_trees,
    has_ARG_data,
    has_ARG_CC_data,
    has_any_ARG_based_data,
    has_both_ARG_sources,
    n_total_ARG_based_local_trees
  )]

  large_repeat_peak_proximity <- merge(
    large_repeat_peak_proximity,
    repeat_arg_cols,
    by = "peak_id",
    all.x = TRUE,
    sort = FALSE
  )

  large_repeat_peak_proximity[is.na(n_ARG_local_trees), n_ARG_local_trees := 0L]
  large_repeat_peak_proximity[is.na(n_ARG_CC_local_trees), n_ARG_CC_local_trees := 0L]
  large_repeat_peak_proximity[is.na(has_ARG_data), has_ARG_data := FALSE]
  large_repeat_peak_proximity[is.na(has_ARG_CC_data), has_ARG_CC_data := FALSE]
  large_repeat_peak_proximity[is.na(has_any_ARG_based_data), has_any_ARG_based_data := FALSE]
  large_repeat_peak_proximity[is.na(has_both_ARG_sources), has_both_ARG_sources := FALSE]
  large_repeat_peak_proximity[
    is.na(n_total_ARG_based_local_trees),
    n_total_ARG_based_local_trees := 0L
  ]
  setorder(large_repeat_peak_proximity, chromosome, peak_start, peak_end)
  fwrite(large_repeat_peak_proximity, large_repeat_peak_proximity_file, sep = "\t")
}

n_Fst_outlier_windows_with_ARG_based_data <- count_entities_with_ARG_based_data(
  arg_window_ARG_coverage
)

n_Fst_peaks_with_ARG_based_data <- count_entities_with_ARG_based_data(
  arg_peak_ARG_coverage
)

arg_ARG_coverage_summary <- rbind(
  data.table(
    entity_level = "Fst_outlier_window",
    total_entities = nrow(arg_window_ARG_coverage),
    n_with_ARG_data = sum(arg_window_ARG_coverage$has_ARG_data),
    n_with_ARG_CC_data = sum(arg_window_ARG_coverage$has_ARG_CC_data),
    n_with_any_ARG_based_data = n_Fst_outlier_windows_with_ARG_based_data,
    n_with_both_ARG_sources = sum(arg_window_ARG_coverage$has_both_ARG_sources),
    n_without_any_ARG_based_data =
      nrow(arg_window_ARG_coverage) - n_Fst_outlier_windows_with_ARG_based_data,
    prop_with_any_ARG_based_data =
      n_Fst_outlier_windows_with_ARG_based_data / nrow(arg_window_ARG_coverage)
  ),
  data.table(
    entity_level = "Xin_Bel_Fst_peak",
    total_entities = nrow(arg_peak_ARG_coverage),
    n_with_ARG_data = sum(arg_peak_ARG_coverage$has_ARG_data),
    n_with_ARG_CC_data = sum(arg_peak_ARG_coverage$has_ARG_CC_data),
    n_with_any_ARG_based_data = n_Fst_peaks_with_ARG_based_data,
    n_with_both_ARG_sources = sum(arg_peak_ARG_coverage$has_both_ARG_sources),
    n_without_any_ARG_based_data =
      nrow(arg_peak_ARG_coverage) - n_Fst_peaks_with_ARG_based_data,
    prop_with_any_ARG_based_data =
      n_Fst_peaks_with_ARG_based_data / nrow(arg_peak_ARG_coverage)
  )
)

fwrite(
  arg_ARG_coverage_summary,
  arg_ARG_coverage_summary_file,
  sep = "\t"
)

status_min_props <- c(NA_real_, arg_min_prop_values)

arg_window_significance_status <- rbindlist(lapply(seq_len(nrow(barplot_stats)), function(i) {
  row <- barplot_stats[i]

  if (row$source == "ARG") {
    threshold_row <- arg_threshold_summary[
      statistic == row$statistic &
        tail == row$tail
    ]

    make_window_significance_status(
      dt = arg_joined,
      stat = row$statistic,
      threshold = threshold_row$threshold[1],
      tail = row$tail,
      source_name = row$source,
      plot_label = row$plot_label,
      window_reference = outlier_window_reference,
      min_props = status_min_props
    )
  } else if (row$source == "ARG_CC") {
    threshold_row <- arg_cc_threshold_summary[
      statistic == row$statistic &
        tail == row$tail
    ]

    make_window_significance_status(
      dt = arg_cc_joined,
      stat = row$statistic,
      threshold = threshold_row$threshold[1],
      tail = row$tail,
      source_name = row$source,
      plot_label = row$plot_label,
      window_reference = outlier_window_reference,
      min_props = status_min_props
    )
  } else {
    stop("Unexpected source in barplot_stats")
  }
}))

arg_window_significance_status[, statistic_order := match(
  statistic,
  barplot_stats$statistic
)]
setorder(
  arg_window_significance_status,
  threshold_type,
  min_prop_significant_local_trees,
  chromosome,
  window_start,
  window_end,
  statistic_order
)

fwrite(
  arg_window_significance_status,
  arg_window_significance_status_file,
  sep = "\t"
)

arg_window_significance_status_summary <- rbind(
  data.table(
    threshold_type = "any_significant_local_tree",
    min_prop_significant_local_trees = NA_real_,
    total_Fst_outlier_windows = nrow(outlier_window_reference),
    n_Fst_outlier_windows_with_ARG_based_data =
      n_Fst_outlier_windows_with_ARG_based_data,
    n_Fst_outlier_windows_with_any_significant_ARG_stat =
      count_unique_significant_windows(arg_window_significance_status, NA_real_)
  ),
  rbindlist(lapply(arg_min_prop_values, function(min_prop) {
    data.table(
      threshold_type = "minimum_proportion_significant_local_trees",
      min_prop_significant_local_trees = min_prop,
      total_Fst_outlier_windows = nrow(outlier_window_reference),
      n_Fst_outlier_windows_with_ARG_based_data =
        n_Fst_outlier_windows_with_ARG_based_data,
      n_Fst_outlier_windows_with_any_significant_ARG_stat =
        count_unique_significant_windows(arg_window_significance_status, min_prop)
    )
  }))
)

fwrite(
  arg_window_significance_status_summary,
  arg_window_significance_status_summary_file,
  sep = "\t"
)

arg_peak_significance_status <- rbindlist(lapply(seq_len(nrow(barplot_stats)), function(i) {
  row <- barplot_stats[i]

  if (row$source == "ARG") {
    threshold_row <- arg_threshold_summary[
      statistic == row$statistic &
        tail == row$tail
    ]

    make_peak_significance_status(
      dt = arg_joined_peaks,
      stat = row$statistic,
      threshold = threshold_row$threshold[1],
      tail = row$tail,
      source_name = row$source,
      plot_label = row$plot_label,
      peak_reference = peak_reference,
      min_props = status_min_props
    )
  } else if (row$source == "ARG_CC") {
    threshold_row <- arg_cc_threshold_summary[
      statistic == row$statistic &
        tail == row$tail
    ]

    make_peak_significance_status(
      dt = arg_cc_joined_peaks,
      stat = row$statistic,
      threshold = threshold_row$threshold[1],
      tail = row$tail,
      source_name = row$source,
      plot_label = row$plot_label,
      peak_reference = peak_reference,
      min_props = status_min_props
    )
  } else {
    stop("Unexpected source in barplot_stats")
  }
}))

arg_peak_significance_status[, statistic_order := match(
  statistic,
  barplot_stats$statistic
)]
setorder(
  arg_peak_significance_status,
  threshold_type,
  min_prop_significant_local_trees,
  peak_id,
  statistic_order
)

fwrite(
  arg_peak_significance_status,
  arg_peak_significance_status_file,
  sep = "\t"
)

arg_peak_significance_status_summary <- rbind(
  data.table(
    threshold_type = "any_significant_local_tree",
    min_prop_significant_local_trees = NA_real_,
    total_Fst_peaks = nrow(peak_reference),
    n_Fst_peaks_with_ARG_based_data =
      n_Fst_peaks_with_ARG_based_data,
    n_Fst_peaks_with_any_significant_ARG_stat =
      count_unique_significant_peaks(arg_peak_significance_status, NA_real_)
  ),
  rbindlist(lapply(arg_min_prop_values, function(min_prop) {
    data.table(
      threshold_type = "minimum_proportion_significant_local_trees",
      min_prop_significant_local_trees = min_prop,
      total_Fst_peaks = nrow(peak_reference),
      n_Fst_peaks_with_ARG_based_data =
        n_Fst_peaks_with_ARG_based_data,
      n_Fst_peaks_with_any_significant_ARG_stat =
        count_unique_significant_peaks(arg_peak_significance_status, min_prop)
    )
  }))
)

fwrite(
  arg_peak_significance_status_summary,
  arg_peak_significance_status_summary_file,
  sep = "\t"
)

arg_window_venn_counts <- summarize_venn_counts_for_thresholds(
  status_dt = arg_window_significance_status,
  entity_col = "window_id",
  entity_level = "Fst_outlier_window",
  min_props = arg_min_prop_values
)

fwrite(
  arg_window_venn_counts,
  arg_window_venn_counts_file,
  sep = "\t"
)

arg_peak_venn_counts <- summarize_venn_counts_for_thresholds(
  status_dt = arg_peak_significance_status,
  entity_col = "peak_id",
  entity_level = "Xin_Bel_Fst_peak",
  min_props = arg_min_prop_values
)

fwrite(
  arg_peak_venn_counts,
  arg_peak_venn_counts_file,
  sep = "\t"
)

# -----------------------------
# Model assignments for Xin-Bel Fst peaks with ARG-based data
# -----------------------------
# =============================================================
# Comprehensive peak-level significance flags for supplementary tables
# =============================================================
# Non-ARG statistics use the same call in both assignment tables: a peak is
# significant when any significant raw/window value occurs in an Fst-outlier
# window belonging to that peak. ARG point statistics added here follow the
# table-specific any-tree or >=25%-within-an-outlier-window rule.

get_threshold_value <- function(summary_dt, statistic_name, tail_name) {
  value <- summary_dt[statistic == statistic_name & tail == tail_name, threshold]
  if (length(value) == 0) NA_real_ else as.numeric(value[1])
}

all_peak_ids <- peaks$peak_id
non_arg_flag_tables <- list()

# Traditional pixy Fst, Dxy, and pi values are already aligned to the outlier windows.
for (stat in pixy_high_tail_stats) {
  # Xin_Bel_Fst uses the same mean + 5 SD threshold that defines the focal
  # outlier windows and its Manhattan track. Other traditional Fst and Dxy
  # statistics use their existing upper 99.9th-percentile control thresholds.
  threshold_i <- if (stat == "Xin_Bel_Fst") {
    threshold
  } else {
    get_threshold_value(pixy_threshold_summary, stat, "upper")
  }
  non_arg_flag_tables[[length(non_arg_flag_tables) + 1L]] <- make_peak_flag_from_pixy_windows(
    outliers, stat, threshold_i,
    "upper", paste0(stat, "_high")
  )
}
for (stat in pixy_low_tail_stats) {
  non_arg_flag_tables[[length(non_arg_flag_tables) + 1L]] <- make_peak_flag_from_pixy_windows(
    outliers, stat,
    get_threshold_value(pixy_threshold_summary, stat, "lower"),
    "lower", paste0(stat, "_low")
  )
}

# Tapajos-Xingu fd.
dt_tmp <- fread(d_stats_file, select = c("scaffold", "start", "end", "sitesUsed", "D", "fd"))
setnames(dt_tmp, c("scaffold", "fd"), c("chromosome", "Tapajos_Xingu_fd"))
dt_tmp[!is.na(D) & D < 0, Tapajos_Xingu_fd := 0]
dt_tmp <- dt_tmp[sitesUsed >= 100 & !is.na(Tapajos_Xingu_fd) & Tapajos_Xingu_fd >= 0 & Tapajos_Xingu_fd <= 1]
dt_tmp[, `:=`(interval_start = as.integer(start - 1L), interval_end = as.integer(end))]
non_arg_flag_tables[[length(non_arg_flag_tables) + 1L]] <- make_peak_flag_from_raw_intervals(
  dt_tmp, "Tapajos_Xingu_fd", outliers,
  get_threshold_value(additional_threshold_summary, "Tapajos_Xingu_fd", "upper"),
  "upper", "Tapajos_Xingu_fd_high"
)
rm(dt_tmp); gc()

# Recombination rate: lower 0.5th-percentile threshold.
dt_tmp <- fread(relernn_file, header = FALSE,
  col.names = c("chromosome", "interval_start", "interval_end", "recombination_rate"))
non_arg_flag_tables[[length(non_arg_flag_tables) + 1L]] <- make_peak_flag_from_raw_intervals(
  dt_tmp, "recombination_rate", outliers,
  get_threshold_value(additional_threshold_summary, "recombination_rate", "lower"),
  "lower", "recombination_rate_low"
)
rm(dt_tmp); gc()

# RAiSD U.
for (spec in list(
  list(file = xingu_raisd_file, stat = "Xingu_RAiSD_u"),
  list(file = belem_raisd_file, stat = "Belem_RAiSD_u")
)) {
  dt_tmp <- fread(spec$file, header = FALSE,
    col.names = c("chromosome", "interval_start", "interval_end", spec$stat))
  non_arg_flag_tables[[length(non_arg_flag_tables) + 1L]] <- make_peak_flag_from_raw_intervals(
    dt_tmp, spec$stat, outliers,
    get_threshold_value(additional_threshold_summary, spec$stat, "upper"),
    "upper", paste0(spec$stat, "_high")
  )
  rm(dt_tmp); gc()
}

# Absolute normalized iHS/nSL, fixed threshold >2.
for (spec in list(
  list(file = xingu_ihs_file, stat = "xingu_norm_ihs", input = "norm_ihs"),
  list(file = xingu_nsl_file, stat = "xingu_norm_nsl", input = "norm_nsl"),
  list(file = belem_ihs_file, stat = "belem_norm_ihs", input = "norm_ihs"),
  list(file = belem_nsl_file, stat = "belem_norm_nsl", input = "norm_nsl")
)) {
  dt_tmp <- fread(spec$file, select = c("chr", "pos", spec$input))
  setnames(dt_tmp, c("chr", spec$input), c("chromosome", spec$stat))
  non_arg_flag_tables[[length(non_arg_flag_tables) + 1L]] <- make_peak_flag_from_raw_points(
    dt_tmp, spec$stat, outliers, 2, "upper", paste0(spec$stat, "_high"), use_abs = TRUE
  )
  rm(dt_tmp); gc()
}

# rCNV values after total.sites >=50 filtering.
dt_tmp <- fread(rcnv_file,
  select = c("scaffold", "start", "end", "total.sites", "dup.sites", "perc.dup.sites"),
  na.strings = c("NA", ""))
setnames(dt_tmp,
  c("scaffold", "start", "end", "total.sites", "dup.sites", "perc.dup.sites"),
  c("chromosome", "interval_start", "interval_end", "total_sites", "dup_sites", "perc_dup_sites"))
dt_tmp <- dt_tmp[!is.na(total_sites) & total_sites >= 50]
for (stat in c("dup_sites", "perc_dup_sites")) {
  non_arg_flag_tables[[length(non_arg_flag_tables) + 1L]] <- make_peak_flag_from_raw_intervals(
    dt_tmp, stat, outliers,
    get_threshold_value(additional_threshold_summary, stat, "upper"),
    "upper", paste0(stat, "_high")
  )
}
rm(dt_tmp); gc()

# U20/U50/Q95 allele statistics: any nonzero value is significant.
for (spec in list(
  list(file = bel_popA_allele_stats_file, suffix = "belem_popA"),
  list(file = xin_popA_allele_stats_file, suffix = "xingu_popA")
)) {
  if (!file.exists(spec$file)) stop("Allele-statistics file not found: ", spec$file)
  dt_tmp <- fread(spec$file)
  required <- c("chromosome", "start", "end", "U20", "U50", "Q95", "informative_sites")
  missing <- setdiff(required, names(dt_tmp))
  if (length(missing) > 0) stop("Missing allele-statistics column(s): ", paste(missing, collapse = ", "))
  dt_tmp[, `:=`(interval_start = as.integer(start), interval_end = as.integer(end))]
  for (stat in c("U20", "U50", "Q95")) {
    col_name <- paste0(stat, "_", spec$suffix, "_high")
    non_arg_flag_tables[[length(non_arg_flag_tables) + 1L]] <- make_peak_flag_from_raw_intervals(
      dt_tmp, stat, outliers, output_col = col_name,
      nonzero_is_significant = TRUE
    )
  }
  rm(dt_tmp); gc()
}

non_arg_peak_flags <- merge_peak_flag_tables(non_arg_flag_tables, all_peak_ids)
rm(non_arg_flag_tables); gc()

# ARG-based Fst thresholds:
#   * Tapajos-Xingu and Tapajos-Belem use upper-tail empirical thresholds
#     calculated from local-tree values assigned to control windows
#     (99.9th percentile; P = 0.001).
#   * Xingu-Belem retains the full filtered-distribution mean + 5 SD rule
#     because it is the focal comparison used to define the Fst peaks.
arg_fst_points <- fread(arg_based_fst_file, select = c(
  "region", "coordinate", "Fst_belem_tapajos", "Fst_belem_xingu", "Fst_tapajos_xingu"
))
setnames(arg_fst_points,
  c("region", "coordinate", "Fst_tapajos_xingu", "Fst_belem_tapajos", "Fst_belem_xingu"),
  c("chromosome", "pos", "Tap_Xin_Fst_ARG_based", "Tap_Bel_Fst_ARG_based", "Xin_Bel_Fst_ARG_based"))
arg_fst_points[, chromosome := standardize_scaffold_names(chromosome)]
arg_fst_points <- filter_points_outside_mask(arg_fst_points, argweaver_mask, pos_col = "pos")
for (stat in c("Tap_Xin_Fst_ARG_based", "Tap_Bel_Fst_ARG_based", "Xin_Bel_Fst_ARG_based")) {
  arg_fst_points[!is.na(get(stat)) & get(stat) < 0, (stat) := 0]
}

arg_fst_thresholds <- rbindlist(list(
  data.table(
    statistic = "Tap_Xin_Fst_ARG_based",
    threshold = get_threshold_value(
      additional_threshold_summary,
      "Tap_Xin_Fst_ARG_based",
      "upper"
    )
  ),
  data.table(
    statistic = "Tap_Bel_Fst_ARG_based",
    threshold = get_threshold_value(
      additional_threshold_summary,
      "Tap_Bel_Fst_ARG_based",
      "upper"
    )
  ),
  data.table(
    statistic = "Xin_Bel_Fst_ARG_based",
    threshold = mean(arg_fst_points$Xin_Bel_Fst_ARG_based, na.rm = TRUE) +
      5 * sd(arg_fst_points$Xin_Bel_Fst_ARG_based, na.rm = TRUE)
  )
))

# Approach 1 model table.
model_assignment_any_core <- make_peak_model_assignments(
  peak_status_dt = arg_peak_significance_status,
  peak_coverage_dt = arg_peak_ARG_coverage,
  min_prop = NA_real_
)
model_assignment_any_core <- merge(
  model_assignment_any_core, selection_peak_significance_flags,
  by = "peak_id", all.x = TRUE, sort = FALSE
)
for (col in grep("^significant_.*_in_Fst_outlier_windows$", names(model_assignment_any_core), value = TRUE)) {
  model_assignment_any_core[is.na(get(col)), (col) := FALSE]
}

arg_extra_any <- list()
for (stat in c("tapajos_enrich", "tapajos_RTH")) {
  tail_i <- if (stat == "tapajos_enrich") "upper" else "lower"
  arg_extra_any[[length(arg_extra_any) + 1L]] <- make_ARG_peak_flag_by_assignment_method(
    arg_joined[, .(chromosome, pos, value = get(stat))][, (stat) := value][, value := NULL],
    stat, outliers, get_threshold_value(arg_threshold_summary, stat, tail_i), tail_i,
    paste0(stat, "_", ifelse(tail_i == "upper", "high", "low")), NA_real_
  )
}
for (stat in arg_fst_thresholds$statistic) {
  arg_extra_any[[length(arg_extra_any) + 1L]] <- make_ARG_peak_flag_by_assignment_method(
    arg_fst_points, stat, outliers,
    arg_fst_thresholds[statistic == stat, threshold][1], "upper", paste0(stat, "_high"), NA_real_
  )
}
arg_extra_any <- merge_peak_flag_tables(arg_extra_any, all_peak_ids)
comprehensive_any_flags <- merge(non_arg_peak_flags, arg_extra_any, by = "peak_id", all = TRUE, sort = FALSE)

# Approach 2 model table.
model_assignment_min25_core <- make_peak_model_assignments(
  peak_status_dt = arg_peak_significance_status,
  peak_coverage_dt = arg_peak_ARG_coverage,
  min_prop = 0.25
)
model_assignment_min25_core <- merge(
  model_assignment_min25_core, selection_peak_significance_flags,
  by = "peak_id", all.x = TRUE, sort = FALSE
)
for (col in grep("^significant_.*_in_Fst_outlier_windows$", names(model_assignment_min25_core), value = TRUE)) {
  model_assignment_min25_core[is.na(get(col)), (col) := FALSE]
}

arg_extra_min25 <- list()
for (stat in c("tapajos_enrich", "tapajos_RTH")) {
  tail_i <- if (stat == "tapajos_enrich") "upper" else "lower"
  arg_extra_min25[[length(arg_extra_min25) + 1L]] <- make_ARG_peak_flag_by_assignment_method(
    arg_joined[, .(chromosome, pos, value = get(stat))][, (stat) := value][, value := NULL],
    stat, outliers, get_threshold_value(arg_threshold_summary, stat, tail_i), tail_i,
    paste0(stat, "_", ifelse(tail_i == "upper", "high", "low")), 0.25
  )
}
for (stat in arg_fst_thresholds$statistic) {
  arg_extra_min25[[length(arg_extra_min25) + 1L]] <- make_ARG_peak_flag_by_assignment_method(
    arg_fst_points, stat, outliers,
    arg_fst_thresholds[statistic == stat, threshold][1], "upper", paste0(stat, "_high"), 0.25
  )
}
arg_extra_min25 <- merge_peak_flag_tables(arg_extra_min25, all_peak_ids)
comprehensive_min25_flags <- merge(non_arg_peak_flags, arg_extra_min25, by = "peak_id", all = TRUE, sort = FALSE)

repeat_peak_ids <- large_repeat_peak_proximity$peak_id
model_assignment_any_summary <- summarize_peak_model_assignments(model_assignment_any_core)
model_assignment_min25_summary <- summarize_peak_model_assignments(model_assignment_min25_core)

model_assignment_any <- expand_model_assignment_to_all_peaks(
  model_assignment_any_core, peaks, comprehensive_any_flags, repeat_peak_ids
)
model_assignment_min25 <- expand_model_assignment_to_all_peaks(
  model_assignment_min25_core, peaks, comprehensive_min25_flags, repeat_peak_ids
)

fwrite(model_assignment_any, model_assignment_any_file, sep = "\t")
fwrite(model_assignment_any_summary, model_assignment_any_summary_file, sep = "\t")
fwrite(model_assignment_min25, model_assignment_min25_file, sep = "\t")
fwrite(model_assignment_min25_summary, model_assignment_min25_summary_file, sep = "\t")


# =============================================================
# Fst-peak significance overrepresentation vs. matched control regions
# All peak-level genomic significance statistics
# =============================================================
if (isTRUE(RUN_FST_PEAK_CONTROL_SUBSAMPLING)) {

  # -----------------------------------------------------------
  # Empirical peak set used for ALL overrepresentation tests
  # -----------------------------------------------------------
  #
  # Restrict the empirical comparison to Fst peaks that contain any ARG-based
  # data. This is the same ARG-coverage definition used elsewhere in the script
  # for deciding whether a peak can be evaluated for ARG-based model support.
  #
  # Importantly, this restriction is applied to EVERY statistic in the
  # overrepresentation analysis, including non-ARG statistics that may have data
  # for all 159 Fst peaks. This keeps the empirical denominator and the peak-size
  # distribution identical across ARG and non-ARG statistics.
  subsampling_peak_ids <- arg_peak_ARG_coverage[
    has_any_ARG_based_data == TRUE,
    peak_id
  ]

  subsampling_peaks <- peaks[
    peak_id %in% subsampling_peak_ids
  ]
  setorder(subsampling_peaks, chromosome, start, end)

  # The current dataset is expected to contain 141 Fst peaks with ARG-based
  # data. Stop if that number changes so that an upstream coverage change does
  # not silently alter the null sampling design.
  if (nrow(subsampling_peaks) != 141L) {
    stop(
      "Expected 141 Fst peaks with ARG-based data for overrepresentation ",
      "subsampling, but found ", nrow(subsampling_peaks), "."
    )
  }

  # Filter both comprehensive empirical significance tables to exactly these
  # same 141 peak IDs. Thus empirical counts for every statistic are based on
  # the identical set of peaks used to define the null peak-length distribution.
  model_assignment_any_subsampling <- model_assignment_any[
    peak_id %in% subsampling_peak_ids
  ]
  model_assignment_min25_subsampling <- model_assignment_min25[
    peak_id %in% subsampling_peak_ids
  ]

  if (nrow(model_assignment_any_subsampling) != 141L ||
      nrow(model_assignment_min25_subsampling) != 141L) {
    stop(
      "The comprehensive empirical significance tables do not contain exactly ",
      "the expected 141 ARG-covered Fst peaks."
    )
  }

  cat(
    "Matched-control overrepresentation analysis restricted to ",
    nrow(subsampling_peaks),
    " Fst peaks with ARG-based data.\n",
    sep = ""
  )

  # Each statistic is represented once. Bookkeeping/model-category booleans and
  # duplicate selection-statistic flags are intentionally excluded.
  #
  # NOTE: Xin_Bel_Fst is included because the comprehensive peak table contains
  # an explicit TRUE/FALSE significance designation for it. Its enrichment test
  # is circular by construction because Xin_Bel_Fst itself defines the peaks;
  # interpret that column only as a pipeline sanity check, not as an independent
  # biological enrichment test.

  subsampling_table_any <- NULL
  subsampling_table_min25 <- NULL

  stat_index <- 0L

  run_one_subsampling_stat <- function(
    statistic_column_name,
    statistic_display_name,
    empirical_col_any,
    empirical_col_min25 = empirical_col_any,
    control_status_any,
    control_status_min25 = control_status_any
  ) {
    stat_index <<- stat_index + 1L
    seed_i <- CONTROL_SUBSAMPLING_SEED + stat_index - 1L

    if (!empirical_col_any %in% names(model_assignment_any_subsampling)) {
      stop("Missing empirical peak significance column: ", empirical_col_any)
    }
    if (!empirical_col_min25 %in% names(model_assignment_min25_subsampling)) {
      stop("Missing empirical peak significance column: ", empirical_col_min25)
    }

    # Observed significant-peak counts are calculated only across the common
    # set of 141 Fst peaks with ARG-based data.
    empirical_any <- sum(
      model_assignment_any_subsampling[[empirical_col_any]] %in% TRUE
    )
    empirical_min25 <- sum(
      model_assignment_min25_subsampling[[empirical_col_min25]] %in% TRUE
    )

    # Add an explicit tail designation to every statistic column in the final
    # matched-control subsampling tables. For most non-ARG statistics this can
    # be inferred directly from the existing empirical TRUE/FALSE column name.
    # The primary ARG statistics use generic "significant_*" empirical columns,
    # so their tested tail is supplied here explicitly.
    tail_designation <- if (grepl("_high$", empirical_col_any)) {
      "high"
    } else if (grepl("_low$", empirical_col_any)) {
      "low"
    } else {
      primary_ARG_tail_lookup <- c(
        "Tap_Xin_RTH_original" = "low",
        "xingu_enrich" = "high",
        "belem_enrich" = "high",
        "xingu_RTH" = "low",
        "belem_RTH" = "low",
        "Xin_Bel_CC_original" = "high"
      )

      unname(primary_ARG_tail_lookup[statistic_column_name])
    }

    if (length(tail_designation) == 0L ||
        is.na(tail_designation) ||
        !tail_designation %in% c("high", "low")) {
      stop(
        "Could not determine high/low tail designation for subsampling statistic: ",
        statistic_column_name
      )
    }

    statistic_column_name_with_tail <- paste0(
      statistic_column_name,
      "_",
      tail_designation
    )

    null_any <- subsample_control_regions_matching_Fst_peak_sizes(
      peaks_dt = subsampling_peaks,
      control_status_dt = control_status_any,
      n_subsamples = CONTROL_SUBSAMPLING_N,
      seed = seed_i
    )

    # The same statistic-specific seed is used for both methods. Whenever the
    # informative candidate sets are identical, this causes the exact same
    # genomic control regions to be sampled, so only the significance rule
    # differs between the two tables.
    null_min25 <- subsample_control_regions_matching_Fst_peak_sizes(
      peaks_dt = subsampling_peaks,
      control_status_dt = control_status_min25,
      n_subsamples = CONTROL_SUBSAMPLING_N,
      seed = seed_i
    )

    subsampling_table_any <<- append_control_subsampling_column(
      output_dt = subsampling_table_any,
      statistic_column_name = statistic_column_name_with_tail,
      statistic_display_name = paste0(statistic_display_name, "_", tail_designation),
      empirical_count = empirical_any,
      subsample_counts = null_any$subsample_counts
    )

    subsampling_table_min25 <<- append_control_subsampling_column(
      output_dt = subsampling_table_min25,
      statistic_column_name = statistic_column_name_with_tail,
      statistic_display_name = paste0(statistic_display_name, "_", tail_designation),
      empirical_count = empirical_min25,
      subsample_counts = null_min25$subsample_counts
    )

    cat(
      "  completed matched-control subsampling for ",
      statistic_display_name,
      "\n",
      sep = ""
    )
  }


  # -----------------------------------------------------------
  # 1. Traditional pixy statistics
  # -----------------------------------------------------------
  for (stat in pixy_high_tail_stats) {
    threshold_i <- if (stat == "Xin_Bel_Fst") {
      threshold
    } else {
      get_threshold_value(pixy_threshold_summary, stat, "upper")
    }

    status_i <- make_aligned_control_window_status(
      control_windows_dt = control_windows,
      stat = stat,
      threshold = threshold_i,
      tail = "upper"
    )

    run_one_subsampling_stat(
      statistic_column_name = stat,
      statistic_display_name = unname(plot_labels[stat]),
      empirical_col_any = paste0(stat, "_high"),
      control_status_any = status_i
    )
  }

  for (stat in pixy_low_tail_stats) {
    status_i <- make_aligned_control_window_status(
      control_windows_dt = control_windows,
      stat = stat,
      threshold = get_threshold_value(pixy_threshold_summary, stat, "lower"),
      tail = "lower"
    )

    run_one_subsampling_stat(
      statistic_column_name = stat,
      statistic_display_name = unname(plot_labels[stat]),
      empirical_col_any = paste0(stat, "_low"),
      control_status_any = status_i
    )
  }


  # -----------------------------------------------------------
  # 2. Non-ARG raw/interval statistics
  # -----------------------------------------------------------

  # Tapajos-Xingu fd.
  dt_tmp <- fread(
    d_stats_file,
    select = c("scaffold", "start", "end", "sitesUsed", "D", "fd")
  )
  setnames(dt_tmp, c("scaffold", "fd"), c("chromosome", "Tapajos_Xingu_fd"))
  dt_tmp[!is.na(D) & D < 0, Tapajos_Xingu_fd := 0]
  dt_tmp <- dt_tmp[
    sitesUsed >= 100 &
      !is.na(Tapajos_Xingu_fd) &
      Tapajos_Xingu_fd >= 0 &
      Tapajos_Xingu_fd <= 1
  ]
  dt_tmp[, `:=`(
    interval_start = as.integer(start - 1L),
    interval_end = as.integer(end)
  )]
  status_i <- make_raw_interval_control_window_status(
    dt = dt_tmp,
    stat = "Tapajos_Xingu_fd",
    control_windows_dt = control_windows,
    threshold = get_threshold_value(
      additional_threshold_summary,
      "Tapajos_Xingu_fd",
      "upper"
    ),
    tail = "upper"
  )
  run_one_subsampling_stat(
    "Tapajos_Xingu_fd",
    unname(plot_labels["Tapajos_Xingu_fd"]),
    "Tapajos_Xingu_fd_high",
    control_status_any = status_i
  )
  rm(dt_tmp, status_i); gc()

  # Recombination rate.
  dt_tmp <- fread(
    relernn_file,
    header = FALSE,
    col.names = c(
      "chromosome", "interval_start", "interval_end", "recombination_rate"
    )
  )
  status_i <- make_raw_interval_control_window_status(
    dt = dt_tmp,
    stat = "recombination_rate",
    control_windows_dt = control_windows,
    threshold = get_threshold_value(
      additional_threshold_summary,
      "recombination_rate",
      "lower"
    ),
    tail = "lower"
  )
  run_one_subsampling_stat(
    "recombination_rate",
    unname(plot_labels["recombination_rate"]),
    "recombination_rate_low",
    control_status_any = status_i
  )
  rm(dt_tmp, status_i); gc()

  # RAiSD U.
  for (spec in list(
    list(file = xingu_raisd_file, stat = "Xingu_RAiSD_u"),
    list(file = belem_raisd_file, stat = "Belem_RAiSD_u")
  )) {
    dt_tmp <- fread(
      spec$file,
      header = FALSE,
      col.names = c("chromosome", "interval_start", "interval_end", spec$stat)
    )
    status_i <- make_raw_interval_control_window_status(
      dt = dt_tmp,
      stat = spec$stat,
      control_windows_dt = control_windows,
      threshold = get_threshold_value(
        additional_threshold_summary,
        spec$stat,
        "upper"
      ),
      tail = "upper"
    )
    run_one_subsampling_stat(
      spec$stat,
      unname(plot_labels[spec$stat]),
      paste0(spec$stat, "_high"),
      control_status_any = status_i
    )
    rm(dt_tmp, status_i); gc()
  }

  # Absolute normalized iHS and nSL; fixed |value| > 2.
  for (spec in list(
    list(file = xingu_ihs_file, stat = "xingu_norm_ihs", input = "norm_ihs"),
    list(file = xingu_nsl_file, stat = "xingu_norm_nsl", input = "norm_nsl"),
    list(file = belem_ihs_file, stat = "belem_norm_ihs", input = "norm_ihs"),
    list(file = belem_nsl_file, stat = "belem_norm_nsl", input = "norm_nsl")
  )) {
    dt_tmp <- fread(spec$file, select = c("chr", "pos", spec$input))
    setnames(dt_tmp, c("chr", spec$input), c("chromosome", spec$stat))
    status_i <- make_raw_point_control_window_status(
      dt = dt_tmp,
      stat = spec$stat,
      control_windows_dt = control_windows,
      threshold = 2,
      tail = "upper",
      use_abs = TRUE
    )
    run_one_subsampling_stat(
      spec$stat,
      unname(plot_labels[spec$stat]),
      paste0(spec$stat, "_high"),
      control_status_any = status_i
    )
    rm(dt_tmp, status_i); gc()
  }

  # rCNV duplicate-site statistics.
  dt_tmp <- fread(
    rcnv_file,
    select = c(
      "scaffold", "start", "end",
      "total.sites", "dup.sites", "perc.dup.sites"
    ),
    na.strings = c("NA", "")
  )
  setnames(
    dt_tmp,
    c(
      "scaffold", "start", "end",
      "total.sites", "dup.sites", "perc.dup.sites"
    ),
    c(
      "chromosome", "interval_start", "interval_end",
      "total_sites", "dup_sites", "perc_dup_sites"
    )
  )
  dt_tmp <- dt_tmp[!is.na(total_sites) & total_sites >= 50]

  for (stat in c("dup_sites", "perc_dup_sites")) {
    status_i <- make_raw_interval_control_window_status(
      dt = dt_tmp,
      stat = stat,
      control_windows_dt = control_windows,
      threshold = get_threshold_value(
        additional_threshold_summary,
        stat,
        "upper"
      ),
      tail = "upper"
    )
    run_one_subsampling_stat(
      stat,
      unname(plot_labels[stat]),
      paste0(stat, "_high"),
      control_status_any = status_i
    )
  }
  rm(dt_tmp, status_i); gc()

  # Allele statistics: any nonzero U20/U50/Q95 value is significant.
  for (spec in list(
    list(
      file = bel_popA_allele_stats_file,
      suffix = "belem_popA",
      display_suffix = "Belem popA"
    ),
    list(
      file = xin_popA_allele_stats_file,
      suffix = "xingu_popA",
      display_suffix = "Xingu popA"
    )
  )) {
    dt_tmp <- fread(spec$file)
    dt_tmp[, `:=`(
      chromosome = standardize_scaffold_names(chromosome),
      interval_start = as.integer(start),
      interval_end = as.integer(end)
    )]

    for (stat in c("U20", "U50", "Q95")) {
      output_name <- paste0(stat, "_", spec$suffix)
      status_i <- make_raw_interval_control_window_status(
        dt = dt_tmp,
        stat = stat,
        control_windows_dt = control_windows,
        nonzero_is_significant = TRUE
      )
      run_one_subsampling_stat(
        statistic_column_name = output_name,
        statistic_display_name = paste(stat, spec$display_suffix),
        empirical_col_any = paste0(output_name, "_high"),
        control_status_any = status_i
      )
    }
    rm(dt_tmp, status_i); gc()
  }


  # -----------------------------------------------------------
  # 3. ARG statistics used for model assignment
  # -----------------------------------------------------------

  # Standard ARG statistics. The same threshold used for empirical peak
  # significance is applied to the control local-tree values.
  arg_subsampling_specs <- data.table(
    stat = c(
      "Tap_Xin_RTH_original",
      "xingu_enrich",
      "belem_enrich",
      "xingu_RTH",
      "belem_RTH",
      "tapajos_enrich",
      "tapajos_RTH"
    ),
    tail = c(
      "lower",
      "upper",
      "upper",
      "lower",
      "lower",
      "upper",
      "lower"
    ),
    empirical_col = c(
      "significant_Tap_Xin_RTH_original",
      "significant_xingu_enrich",
      "significant_belem_enrich",
      "significant_xingu_RTH",
      "significant_belem_RTH",
      "tapajos_enrich_high",
      "tapajos_RTH_low"
    )
  )

  for (spec_i in seq_len(nrow(arg_subsampling_specs))) {
    spec <- arg_subsampling_specs[spec_i]
    threshold_i <- get_threshold_value(
      arg_threshold_summary,
      spec$stat,
      spec$tail
    )

    # The six primary model-assignment statistics use the inclusive >= / <=
    # rule implemented by make_peak_significance_status(). The two additional
    # Tapajos statistics use make_ARG_peak_flag_by_assignment_method(), which is
    # strict > / <. Handle those separately to mirror the empirical calls.
    if (spec$stat %in% c("tapajos_enrich", "tapajos_RTH")) {
      arg_points_i <- arg_joined[, .(
        chromosome,
        pos,
        value = get(spec$stat)
      )]
      setnames(arg_points_i, "value", spec$stat)

      status_any <- make_point_proportion_control_window_status(
        dt = arg_points_i,
        stat = spec$stat,
        control_windows_dt = control_windows,
        threshold = threshold_i,
        tail = spec$tail,
        min_prop = NA_real_,
        inclusive_threshold = FALSE
      )
      status_min25 <- make_point_proportion_control_window_status(
        dt = arg_points_i,
        stat = spec$stat,
        control_windows_dt = control_windows,
        threshold = threshold_i,
        tail = spec$tail,
        min_prop = 0.25,
        inclusive_threshold = FALSE
      )
      rm(arg_points_i)
    } else {
      status_any <- make_ARG_control_window_status(
        arg_joined_dt = arg_joined,
        control_windows_dt = control_windows,
        stat = spec$stat,
        threshold = threshold_i,
        tail = spec$tail,
        min_prop = NA_real_
      )
      status_min25 <- make_ARG_control_window_status(
        arg_joined_dt = arg_joined,
        control_windows_dt = control_windows,
        stat = spec$stat,
        threshold = threshold_i,
        tail = spec$tail,
        min_prop = 0.25
      )
    }

    run_one_subsampling_stat(
      statistic_column_name = spec$stat,
      statistic_display_name = unname(plot_labels[spec$stat]),
      empirical_col_any = spec$empirical_col,
      control_status_any = status_any,
      control_status_min25 = status_min25
    )
    rm(status_any, status_min25); gc()
  }

  # Xingu-Belem RCC is stored in the separate ARG_CC table.
  cc_threshold_i <- get_threshold_value(
    arg_cc_threshold_summary,
    "Xin_Bel_CC_original",
    "upper"
  )
  status_any <- make_ARG_control_window_status(
    arg_joined_dt = arg_cc_joined,
    control_windows_dt = control_windows,
    stat = "Xin_Bel_CC_original",
    threshold = cc_threshold_i,
    tail = "upper",
    min_prop = NA_real_
  )
  status_min25 <- make_ARG_control_window_status(
    arg_joined_dt = arg_cc_joined,
    control_windows_dt = control_windows,
    stat = "Xin_Bel_CC_original",
    threshold = cc_threshold_i,
    tail = "upper",
    min_prop = 0.25
  )
  run_one_subsampling_stat(
    "Xin_Bel_CC_original",
    unname(plot_labels["Xin_Bel_CC_original"]),
    "significant_Xin_Bel_CC_original",
    control_status_any = status_any,
    control_status_min25 = status_min25
  )
  rm(status_any, status_min25); gc()


  # -----------------------------------------------------------
  # 4. ARG-based branch Fst statistics
  # -----------------------------------------------------------
  for (stat in arg_fst_thresholds$statistic) {
    threshold_i <- arg_fst_thresholds[
      statistic == stat,
      threshold
    ][1]

    status_any <- make_point_proportion_control_window_status(
      dt = arg_fst_points,
      stat = stat,
      control_windows_dt = control_windows,
      threshold = threshold_i,
      tail = "upper",
      min_prop = NA_real_,
      inclusive_threshold = FALSE
    )
    status_min25 <- make_point_proportion_control_window_status(
      dt = arg_fst_points,
      stat = stat,
      control_windows_dt = control_windows,
      threshold = threshold_i,
      tail = "upper",
      min_prop = 0.25,
      inclusive_threshold = FALSE
    )

    run_one_subsampling_stat(
      statistic_column_name = stat,
      statistic_display_name = unname(plot_labels[stat]),
      empirical_col_any = paste0(stat, "_high"),
      control_status_any = status_any,
      control_status_min25 = status_min25
    )
    rm(status_any, status_min25); gc()
  }


  # Restore the requested row order explicitly after repeated merges.
  descriptor_order <- c(
    "Statistic_name",
    "Empirical_count",
    "P_0.05_threshold",
    "P_0.01_threshold",
    "P_0.001_threshold",
    paste0("Subsample ", seq_len(CONTROL_SUBSAMPLING_N))
  )

  subsampling_table_any[
    , descriptor_order_index := match(Descriptor, descriptor_order)
  ]
  setorder(subsampling_table_any, descriptor_order_index)
  subsampling_table_any[, descriptor_order_index := NULL]
  setcolorder(
    subsampling_table_any,
    c("Descriptor", setdiff(names(subsampling_table_any), "Descriptor"))
  )

  subsampling_table_min25[
    , descriptor_order_index := match(Descriptor, descriptor_order)
  ]
  setorder(subsampling_table_min25, descriptor_order_index)
  subsampling_table_min25[, descriptor_order_index := NULL]
  setcolorder(
    subsampling_table_min25,
    c("Descriptor", setdiff(names(subsampling_table_min25), "Descriptor"))
  )

  fwrite(
    subsampling_table_any,
    fst_peak_control_subsampling_any_file,
    sep = "\t"
  )
  fwrite(
    subsampling_table_min25,
    fst_peak_control_subsampling_min25_file,
    sep = "\t"
  )

  cat(
    "Completed Fst-peak vs. matched-control significance subsampling for ",
    ncol(subsampling_table_any) - 1L,
    " statistics using ",
    CONTROL_SUBSAMPLING_N,
    " replicates per statistic.\n",
    sep = ""
  )
}


# Create the three peak-size distribution PDFs together under one dedicated
# switch. These are generated after model assignments are available so the two
# colored versions can associate each peak with its assignment category.
if (isTRUE(CREATE_PEAK_SIZE_DISTRIBUTION_PDFS)) {
  make_peak_size_distribution_barplot(
    peaks_dt = peaks,
    outfile = figure_peak_size_distribution_pdf
  )

  make_model_colored_peak_size_distribution_barplot(
    peaks_dt = peaks,
    model_dt = model_assignment_any,
    outfile = figure_peak_size_distribution_any_model_pdf,
    assignment_label = "at least one significant local tree in a peak"
  )

  make_model_colored_peak_size_distribution_barplot(
    peaks_dt = peaks,
    model_dt = model_assignment_min25,
    outfile = figure_peak_size_distribution_min25_model_pdf,
    assignment_label = ">=25% significant local trees in an Fst-outlier window"
  )
}

rm(
  arg_fst_points, arg_fst_thresholds,
  arg_extra_any, arg_extra_min25,
  comprehensive_any_flags, comprehensive_min25_flags,
  non_arg_peak_flags
)
gc()

model_category_definitions <- make_model_category_definitions()

write_model_category_definitions(
  definitions_dt = model_category_definitions,
  xlsx_file = model_category_definitions_xlsx,
  tsv_file = model_category_definitions_tsv
)

arg_significant_window_counts <- rbindlist(lapply(seq_len(nrow(barplot_stats)), function(i) {
  row <- barplot_stats[i]

  if (row$source == "ARG") {
    threshold_row <- arg_threshold_summary[
      statistic == row$statistic &
        tail == row$tail
    ]

    count_significant_windows(
      dt = arg_joined,
      stat = row$statistic,
      threshold = threshold_row$threshold[1],
      tail = row$tail,
      source_name = row$source
    )
  } else if (row$source == "ARG_CC") {
    threshold_row <- arg_cc_threshold_summary[
      statistic == row$statistic &
        tail == row$tail
    ]

    count_significant_windows(
      dt = arg_cc_joined,
      stat = row$statistic,
      threshold = threshold_row$threshold[1],
      tail = row$tail,
      source_name = row$source
    )
  } else {
    stop("Unexpected source in barplot_stats")
  }
}))

arg_significant_window_counts <- merge(
  arg_significant_window_counts,
  barplot_stats[, .(
    statistic,
    direction_label,
    plot_label,
    fill_color
  )],
  by = "statistic",
  all.x = TRUE,
  sort = FALSE
)

arg_significant_window_counts[, statistic_order := match(
  statistic,
  barplot_stats$statistic
)]

setorder(
  arg_significant_window_counts,
  window_class,
  statistic_order
)

fwrite(
  arg_significant_window_counts,
  arg_significant_window_counts_file,
  sep = "\t"
)

arg_significant_window_barplot_data <- arg_significant_window_counts[
  window_class == "Xin_Bel_Fst_outlier"
]

setorder(
  arg_significant_window_barplot_data,
  statistic_order
)

make_significant_window_barplot(
  arg_significant_window_barplot_data,
  arg_significant_window_barplot_pdf,
  annotation_lines = c(
    paste0("Total number of Fst outlier windows = ", nrow(outlier_window_reference)),
    paste0(
      "Total number of Fst outlier windows with ARG-based data = ",
      n_Fst_outlier_windows_with_ARG_based_data
    ),
    paste0(
      "Total number of Fst outlier windows with a significant ARG-based statistic = ",
      count_unique_significant_windows(arg_window_significance_status, NA_real_)
    )
  ),
  venn_status_dt = arg_window_significance_status,
  venn_entity_col = "window_id",
  venn_min_prop = NA_real_,
  venn_entity_label = "Fst outlier windows",
  venn_title_context = "Fst outlier windows"
)

# -----------------------------
# Count Fst outlier/control windows containing
# significant local trees for selected ARG statistics,
# requiring minimum proportions of significant local trees
# in the window: 5%, 10%, 20%, 25%, 30%, 40%, 50%, 75%, and 100%
# -----------------------------
arg_significant_window_counts_minprop_list <- list()

for (min_prop in arg_min_prop_values) {
  suffix <- paste0("min", as.integer(min_prop * 100), "percent")

  tmp_counts <- rbindlist(lapply(seq_len(nrow(barplot_stats)), function(i) {
    row <- barplot_stats[i]

    if (row$source == "ARG") {
      threshold_row <- arg_threshold_summary[
        statistic == row$statistic &
          tail == row$tail
      ]

      count_significant_windows_min_prop(
        dt = arg_joined,
        stat = row$statistic,
        threshold = threshold_row$threshold[1],
        tail = row$tail,
        source_name = row$source,
        min_prop = min_prop
      )
    } else if (row$source == "ARG_CC") {
      threshold_row <- arg_cc_threshold_summary[
        statistic == row$statistic &
          tail == row$tail
      ]

      count_significant_windows_min_prop(
        dt = arg_cc_joined,
        stat = row$statistic,
        threshold = threshold_row$threshold[1],
        tail = row$tail,
        source_name = row$source,
        min_prop = min_prop
      )
    } else {
      stop("Unexpected source in barplot_stats")
    }
  }))

  tmp_counts <- merge(
    tmp_counts,
    barplot_stats[, .(
      statistic,
      direction_label,
      plot_label,
      fill_color
    )],
    by = "statistic",
    all.x = TRUE,
    sort = FALSE
  )

  tmp_counts[, statistic_order := match(
    statistic,
    barplot_stats$statistic
  )]

  setorder(
    tmp_counts,
    window_class,
    statistic_order
  )

  fwrite(
    tmp_counts,
    arg_significant_window_counts_minprop_files[[suffix]],
    sep = "	"
  )

  tmp_barplot_data <- tmp_counts[
    window_class == "Xin_Bel_Fst_outlier"
  ]

  setorder(
    tmp_barplot_data,
    statistic_order
  )

  make_significant_window_barplot(
    tmp_barplot_data,
    arg_significant_window_barplot_minprop_pdfs[[suffix]],
    annotation_lines = c(
      paste0("Total number of Fst outlier windows = ", nrow(outlier_window_reference)),
      paste0(
        "Total number of Fst outlier windows with ARG-based data = ",
        n_Fst_outlier_windows_with_ARG_based_data
      ),
      paste0(
        "Total number of Fst outlier windows with a significant ARG-based statistic at >= ",
        as.integer(min_prop * 100),
        "% frequency = ",
        count_unique_significant_windows(arg_window_significance_status, min_prop)
      )
    ),
    venn_status_dt = arg_window_significance_status,
    venn_entity_col = "window_id",
    venn_min_prop = min_prop,
    venn_entity_label = "Fst outlier windows",
    venn_title_context = paste0("Fst outlier windows with >=", as.integer(min_prop * 100), "% significant local trees")
  )

  arg_significant_window_counts_minprop_list[[suffix]] <- tmp_counts
}

# Preserve the original 5% object names for backwards compatibility
arg_significant_window_counts_min5 <- arg_significant_window_counts_minprop_list[["min5percent"]]
arg_significant_window_counts_min5_file <- arg_significant_window_counts_minprop_files[["min5percent"]]
arg_significant_window_barplot_min5_pdf <- arg_significant_window_barplot_minprop_pdfs[["min5percent"]]

# -----------------------------
# Count Xin-Bel Fst peaks containing any significant local tree
# for selected ARG statistics
# -----------------------------
arg_significant_peak_counts <- rbindlist(lapply(seq_len(nrow(barplot_stats)), function(i) {
  row <- barplot_stats[i]

  if (row$source == "ARG") {
    threshold_row <- arg_threshold_summary[
      statistic == row$statistic &
        tail == row$tail
    ]

    count_significant_peaks(
      dt = arg_joined_peaks,
      stat = row$statistic,
      threshold = threshold_row$threshold[1],
      tail = row$tail,
      source_name = row$source
    )
  } else if (row$source == "ARG_CC") {
    threshold_row <- arg_cc_threshold_summary[
      statistic == row$statistic &
        tail == row$tail
    ]

    count_significant_peaks(
      dt = arg_cc_joined_peaks,
      stat = row$statistic,
      threshold = threshold_row$threshold[1],
      tail = row$tail,
      source_name = row$source
    )
  } else {
    stop("Unexpected source in barplot_stats")
  }
}))

arg_significant_peak_counts <- merge(
  arg_significant_peak_counts,
  barplot_stats[, .(
    statistic,
    direction_label,
    plot_label,
    fill_color
  )],
  by = "statistic",
  all.x = TRUE,
  sort = FALSE
)

arg_significant_peak_counts[, statistic_order := match(
  statistic,
  barplot_stats$statistic
)]

setorder(
  arg_significant_peak_counts,
  statistic_order
)

fwrite(
  arg_significant_peak_counts,
  arg_significant_peak_counts_file,
  sep = "\t"
)

arg_significant_peak_barplot_data <- copy(arg_significant_peak_counts)

setorder(
  arg_significant_peak_barplot_data,
  statistic_order
)

make_significant_window_barplot(
  arg_significant_peak_barplot_data,
  arg_significant_peak_barplot_pdf,
  y_col = "n_significant_peaks",
  ylab = "Number of Xin-Bel Fst peaks",
  main_title = "Xin-Bel Fst peaks containing significant ARG local trees",
  annotation_lines = c(
    paste0("Total number of Fst peaks = ", nrow(peak_reference)),
    paste0(
      "Total number of Fst peaks with ARG-based data = ",
      n_Fst_peaks_with_ARG_based_data
    ),
    paste0(
      "Total number of Fst peaks with a significant ARG-based statistic = ",
      count_unique_significant_peaks(arg_peak_significance_status, NA_real_)
    )
  ),
  venn_status_dt = arg_peak_significance_status,
  venn_entity_col = "peak_id",
  venn_min_prop = NA_real_,
  venn_entity_label = "Fst peaks",
  venn_title_context = "Xin-Bel Fst peaks"
)

# -----------------------------
# Count Xin-Bel Fst peaks containing significant local trees,
# requiring minimum proportions of significant local trees
# within each peak: 5%, 10%, 20%, 25%, 30%, 40%, 50%, 75%, and 100%
# -----------------------------
arg_significant_peak_counts_minprop_list <- list()

for (min_prop in arg_min_prop_values) {
  suffix <- paste0("min", as.integer(min_prop * 100), "percent")

  tmp_peak_counts <- rbindlist(lapply(seq_len(nrow(barplot_stats)), function(i) {
    row <- barplot_stats[i]

    if (row$source == "ARG") {
      threshold_row <- arg_threshold_summary[
        statistic == row$statistic &
          tail == row$tail
      ]

      count_significant_peaks_min_prop(
        dt = arg_joined_peaks,
        stat = row$statistic,
        threshold = threshold_row$threshold[1],
        tail = row$tail,
        source_name = row$source,
        min_prop = min_prop
      )
    } else if (row$source == "ARG_CC") {
      threshold_row <- arg_cc_threshold_summary[
        statistic == row$statistic &
          tail == row$tail
      ]

      count_significant_peaks_min_prop(
        dt = arg_cc_joined_peaks,
        stat = row$statistic,
        threshold = threshold_row$threshold[1],
        tail = row$tail,
        source_name = row$source,
        min_prop = min_prop
      )
    } else {
      stop("Unexpected source in barplot_stats")
    }
  }))

  tmp_peak_counts <- merge(
    tmp_peak_counts,
    barplot_stats[, .(
      statistic,
      direction_label,
      plot_label,
      fill_color
    )],
    by = "statistic",
    all.x = TRUE,
    sort = FALSE
  )

  tmp_peak_counts[, statistic_order := match(
    statistic,
    barplot_stats$statistic
  )]

  setorder(
    tmp_peak_counts,
    statistic_order
  )

  fwrite(
    tmp_peak_counts,
    arg_significant_peak_counts_minprop_files[[suffix]],
    sep = "\t"
  )

  tmp_peak_barplot_data <- copy(tmp_peak_counts)

  setorder(
    tmp_peak_barplot_data,
    statistic_order
  )

  make_significant_window_barplot(
    tmp_peak_barplot_data,
    arg_significant_peak_barplot_minprop_pdfs[[suffix]],
    y_col = "n_significant_peaks",
    ylab = "Number of Xin-Bel Fst peaks",
    main_title = paste0(
      "Xin-Bel Fst peaks with >=",
      as.integer(min_prop * 100),
      "% significant ARG local trees"
    ),
    annotation_lines = c(
      paste0("Total number of Fst peaks = ", nrow(peak_reference)),
      paste0(
        "Total number of Fst peaks with ARG-based data = ",
        n_Fst_peaks_with_ARG_based_data
      ),
      paste0(
        "Total number of Fst peaks with a significant ARG-based statistic at >= ",
        as.integer(min_prop * 100),
        "% frequency = ",
        count_unique_significant_peaks(arg_peak_significance_status, min_prop)
      )
    ),
    venn_status_dt = arg_peak_significance_status,
    venn_entity_col = "peak_id",
    venn_min_prop = min_prop,
    venn_entity_label = "Fst peaks",
    venn_title_context = paste0("Xin-Bel Fst peaks with >=", as.integer(min_prop * 100), "% significant local trees")
  )

  arg_significant_peak_counts_minprop_list[[suffix]] <- tmp_peak_counts
}

arg_significant_peak_counts_min5 <- arg_significant_peak_counts_minprop_list[["min5percent"]]
arg_significant_peak_counts_min5_file <- arg_significant_peak_counts_minprop_files[["min5percent"]]
arg_significant_peak_barplot_min5_pdf <- arg_significant_peak_barplot_minprop_pdfs[["min5percent"]]



# -----------------------------
# One-off ARG barplot P-value sensitivity figures
# -----------------------------
# These combined figures show how ARG-significant Fst outlier window and
# peak counts change when empirical control-window thresholds are defined at
# P = 0.0001, P = 0.001, and P = 0.01. For one-sided tests this corresponds
# to lower-tail percentiles of 0.01st, 0.1st, and 1st, or upper-tail
# percentiles of 99.99th, 99.9th, and 99th, respectively.
if (isTRUE(CREATE_FIGURES_PDFS)) {
make_ARG_p_sensitivity_combined_barplot_pdf(
  outfile = figure_ARG_barplot_p_sensitivity_any_pdf,
  significance_definition = "any_significant_local_tree",
  p_values = c(0.0001, 0.001, 0.01),
  min_prop = NA_real_,
  arg_joined_dt = arg_joined,
  arg_cc_joined_dt = arg_cc_joined,
  arg_joined_peaks_dt = arg_joined_peaks,
  arg_cc_joined_peaks_dt = arg_cc_joined_peaks,
  barplot_stats_dt = barplot_stats
)

make_ARG_p_sensitivity_combined_barplot_pdf(
  outfile = figure_ARG_barplot_p_sensitivity_min25_pdf,
  significance_definition = "minimum_proportion_significant_local_trees",
  p_values = c(0.0001, 0.001, 0.01),
  min_prop = 0.25,
  arg_joined_dt = arg_joined,
  arg_cc_joined_dt = arg_cc_joined,
  arg_joined_peaks_dt = arg_joined_peaks,
  arg_cc_joined_peaks_dt = arg_cc_joined_peaks,
  barplot_stats_dt = barplot_stats
)
}



# =============================================================
# Multi-iteration ARG significance summaries
# =============================================================
# These tables summarize the bar-plot-style counts and model assignments
# for all 50 saved MCMC iterations, every 10 iterations from 1510 through 2000.
# Significance thresholds remain based on the 2000th iteration.

all_ARG_MCMC_iterations <- seq.int(1510L, 2000L, by = 10L)

# This is intentionally OFF by default because reading/summarizing all 50
# genome-wide ARG iterations is one of the slowest parts of the pipeline.
if (isTRUE(RUN_MULTI_ITER_ARG_SUMMARIES)) {

multi_iter_ARG_summary <- summarize_multi_iter_ARG_significance(
  iterations = all_ARG_MCMC_iterations,
  barplot_stats = barplot_stats,
  arg_stat_dir = arg_stat_dir,
  arg_cc_stat_dir = arg_cc_stat_dir,
  argweaver_mask = argweaver_mask,
  comparison_windows = comparison_windows,
  peak_window_lookup = peak_window_lookup,
  window_reference = outlier_window_reference,
  peak_reference = peak_reference,
  peak_coverage_dt = arg_peak_ARG_coverage,
  arg_threshold_summary = arg_threshold_summary,
  arg_cc_threshold_summary = arg_cc_threshold_summary,
  plot_labels = plot_labels
)

fwrite(
  multi_iter_ARG_summary$window_counts,
  arg_multi_iter_window_counts_file,
  sep = "\t"
)

fwrite(
  multi_iter_ARG_summary$peak_counts,
  arg_multi_iter_peak_counts_file,
  sep = "\t"
)

fwrite(
  multi_iter_ARG_summary$model_assignments,
  arg_multi_iter_model_assignments_file,
  sep = "\t"
)

fwrite(
  multi_iter_ARG_summary$model_summary,
  arg_multi_iter_model_assignment_summary_file,
  sep = "\t"
)


model_assignment_ARG_statistics <- c(
  "xingu_enrich",
  "belem_enrich",
  "xingu_RTH",
  "belem_RTH",
  "Xin_Bel_CC_original",
  "Tap_Xin_RTH_original"
)

ARG_stat_summary_50iter_any <- make_50_iteration_ARG_statistic_summary(
  window_counts_dt = multi_iter_ARG_summary$window_counts,
  peak_counts_dt = multi_iter_ARG_summary$peak_counts,
  threshold_type_filter = "any_significant_local_tree",
  iterations = all_ARG_MCMC_iterations,
  statistics = model_assignment_ARG_statistics
)

ARG_stat_summary_50iter_min25 <- make_50_iteration_ARG_statistic_summary(
  window_counts_dt = multi_iter_ARG_summary$window_counts,
  peak_counts_dt = multi_iter_ARG_summary$peak_counts,
  threshold_type_filter = "minimum_proportion_significant_local_trees",
  iterations = all_ARG_MCMC_iterations,
  statistics = model_assignment_ARG_statistics,
  min_prop_filter = 0.25
)

fwrite(ARG_stat_summary_50iter_any, arg_50iter_stat_summary_any_file, sep = "\t")
fwrite(ARG_stat_summary_50iter_min25, arg_50iter_stat_summary_min25_file, sep = "\t")

# -----------------------------
# Add multi-iteration range intervals to the four primary ARG barplots
# -----------------------------
# The bars themselves remain the 2000th-iteration values already calculated
# above. The intervals show the minimum-to-maximum number of significant
# windows or peaks observed across all 50 MCMC iterations from 1510 through 2000.
if (exists("CREATE_BARPLOT_PDFS") && isTRUE(CREATE_BARPLOT_PDFS)) {
  window_any_intervals <- make_multi_iter_barplot_interval_dt(
    multi_iter_counts_dt = multi_iter_ARG_summary$window_counts,
    y_col = "n_significant_windows",
    threshold_type_value = "any_significant_local_tree",
    iterations = all_ARG_MCMC_iterations
  )

  window_min25_intervals <- make_multi_iter_barplot_interval_dt(
    multi_iter_counts_dt = multi_iter_ARG_summary$window_counts,
    y_col = "n_significant_windows",
    threshold_type_value = "minimum_proportion_significant_local_trees",
    min_prop_value = 0.25,
    iterations = all_ARG_MCMC_iterations
  )

  peak_any_intervals <- make_multi_iter_barplot_interval_dt(
    multi_iter_counts_dt = multi_iter_ARG_summary$peak_counts,
    y_col = "n_significant_peaks",
    threshold_type_value = "any_significant_local_tree",
    iterations = all_ARG_MCMC_iterations
  )

  peak_min25_intervals <- make_multi_iter_barplot_interval_dt(
    multi_iter_counts_dt = multi_iter_ARG_summary$peak_counts,
    y_col = "n_significant_peaks",
    threshold_type_value = "minimum_proportion_significant_local_trees",
    min_prop_value = 0.25,
    iterations = all_ARG_MCMC_iterations
  )

  make_significant_window_barplot(
    arg_significant_window_barplot_data,
    arg_significant_window_barplot_pdf,
    annotation_lines = c(
      paste0("Total number of Fst outlier windows = ", nrow(outlier_window_reference)),
      paste0(
        "Total number of Fst outlier windows with ARG-based data = ",
        n_Fst_outlier_windows_with_ARG_based_data
      ),
      paste0(
        "Total number of Fst outlier windows with a significant ARG-based statistic = ",
        count_unique_significant_windows(arg_window_significance_status, NA_real_)
      ),
      "Intervals show min-max counts across 50 MCMC iterations (1510-2000)"
    ),
    venn_status_dt = arg_window_significance_status,
    venn_entity_col = "window_id",
    venn_min_prop = NA_real_,
    venn_entity_label = "Fst outlier windows",
    venn_title_context = "Fst outlier windows",
    interval_dt = window_any_intervals,
    interval_label = "range across MCMC iterations"
  )

  make_significant_window_barplot(
    arg_significant_window_counts_minprop_list[["min25percent"]][window_class == "Xin_Bel_Fst_outlier"],
    arg_significant_window_barplot_minprop_pdfs[["min25percent"]],
    annotation_lines = c(
      paste0("Total number of Fst outlier windows = ", nrow(outlier_window_reference)),
      paste0(
        "Total number of Fst outlier windows with ARG-based data = ",
        n_Fst_outlier_windows_with_ARG_based_data
      ),
      paste0(
        "Total number of Fst outlier windows with a significant ARG-based statistic at >= 25% frequency = ",
        count_unique_significant_windows(arg_window_significance_status, 0.25)
      ),
      "Intervals show min-max counts across 50 MCMC iterations (1510-2000)"
    ),
    venn_status_dt = arg_window_significance_status,
    venn_entity_col = "window_id",
    venn_min_prop = 0.25,
    venn_entity_label = "Fst outlier windows",
    venn_title_context = "Fst outlier windows with >=25% significant local trees",
    interval_dt = window_min25_intervals,
    interval_label = "range across MCMC iterations"
  )

  make_significant_window_barplot(
    arg_significant_peak_barplot_data,
    arg_significant_peak_barplot_pdf,
    y_col = "n_significant_peaks",
    ylab = "Number of Xin-Bel Fst peaks",
    main_title = "Xin-Bel Fst peaks containing significant ARG local trees",
    annotation_lines = c(
      paste0("Total number of Fst peaks = ", nrow(peak_reference)),
      paste0(
        "Total number of Fst peaks with ARG-based data = ",
        n_Fst_peaks_with_ARG_based_data
      ),
      paste0(
        "Total number of Fst peaks with a significant ARG-based statistic = ",
        count_unique_significant_peaks(arg_peak_significance_status, NA_real_)
      ),
      "Intervals show min-max counts across 50 MCMC iterations (1510-2000)"
    ),
    venn_status_dt = arg_peak_significance_status,
    venn_entity_col = "peak_id",
    venn_min_prop = NA_real_,
    venn_entity_label = "Fst peaks",
    venn_title_context = "Xin-Bel Fst peaks",
    interval_dt = peak_any_intervals,
    interval_label = "range across MCMC iterations"
  )

  make_significant_window_barplot(
    arg_significant_peak_counts_minprop_list[["min25percent"]],
    arg_significant_peak_barplot_minprop_pdfs[["min25percent"]],
    y_col = "n_significant_peaks",
    ylab = "Number of Xin-Bel Fst peaks",
    main_title = "Xin-Bel Fst peaks with >=25% significant ARG local trees",
    annotation_lines = c(
      paste0("Total number of Fst peaks = ", nrow(peak_reference)),
      paste0(
        "Total number of Fst peaks with ARG-based data = ",
        n_Fst_peaks_with_ARG_based_data
      ),
      paste0(
        "Total number of Fst peaks with a significant ARG-based statistic at >= 25% frequency = ",
        count_unique_significant_peaks(arg_peak_significance_status, 0.25)
      ),
      "Intervals show min-max counts across 50 MCMC iterations (1510-2000)"
    ),
    venn_status_dt = arg_peak_significance_status,
    venn_entity_col = "peak_id",
    venn_min_prop = 0.25,
    venn_entity_label = "Fst peaks",
    venn_title_context = "Xin-Bel Fst peaks with >=25% significant local trees",
    interval_dt = peak_min25_intervals,
    interval_label = "range across MCMC iterations"
  )
}

multi_iter_model_assignments_all_iterations <- copy(
  multi_iter_ARG_summary$model_assignments
)

rm(multi_iter_ARG_summary)
gc()

} else {
  # Keep a defined empty object so downstream code can safely test whether
  # multi-iteration model assignments are available without triggering the
  # expensive summary calculation.
  multi_iter_model_assignments_all_iterations <- data.table()
  cat(
    "Skipping multi-iteration ARG summary tables because ",
    "RUN_MULTI_ITER_ARG_SUMMARIES = FALSE.\n",
    sep = ""
  )
}


# =============================================================
# PCA of heat-map percentile profiles for Xin-Bel Fst peaks
# =============================================================
# PCA is based on the same peak-level representative percentile scores used
# in the heat map. Points are colored by the 2000th-iteration model assignment,
# with separate PCA PDFs for the two model-categorization strategies. These PCA
# plots are generated before Manhattan plots.
heatmap_pca_any_scores <- make_peak_percentile_pca_pdf(
  heatmap_peak_representative = heatmap_peak_representative,
  stat_order = heatmap_stat_order,
  plot_labels = plot_labels,
  model_dt = model_assignment_any,
  outfile_pdf = heatmap_pca_any_pdf,
  outfile_scores = heatmap_pca_any_scores_file
)

heatmap_pca_min25_scores <- make_peak_percentile_pca_pdf(
  heatmap_peak_representative = heatmap_peak_representative,
  stat_order = heatmap_stat_order,
  plot_labels = plot_labels,
  model_dt = model_assignment_min25,
  outfile_pdf = heatmap_pca_min25_pdf,
  outfile_scores = heatmap_pca_min25_scores_file
)

rm(heatmap_pca_any_scores, heatmap_pca_min25_scores)
gc()


# =============================================================
# Manuscript figure helpers
# =============================================================
make_oneoff_manhattan_tree_figure <- function(
  peak_id_to_plot,
  outfile,
  figure_stats,
  model_dt,
  manhattan_native_values,
  manhattan_window_reference,
  thresholds_dt,
  peaks_dt,
  plot_labels,
  scaffold_lengths,
  peak_window_lookup_dt,
  tree_map,
  tree_base_dir,
  multi_iter = FALSE,
  multi_iter_stats = character(),
  multi_iter_colors = NULL,
  multi_iter_iterations = NULL,
  multi_iter_arg_stat_dir = NULL,
  multi_iter_arg_cc_stat_dir = NULL,
  multi_iter_argweaver_mask = NULL,
  focal_rule = c("max_fst", "closest_to_peak_center"),
  figure_title = NULL,
  rolling_half_window_overrides = NULL
) {
  focal_rule <- match.arg(focal_rule)
  model_row <- model_dt[as.character(peak_id) == as.character(peak_id_to_plot)][1]
  if (nrow(model_row) == 0 || is.na(model_row$peak_id[1])) {
    warning("Peak ", peak_id_to_plot, " not found in model table; skipping ", outfile)
    return(invisible(data.table()))
  }

  chrom_i <- model_row$chromosome[1]
  scaffold_max <- scaffold_lengths[chromosome == chrom_i, scaffold_max_end]
  if (length(scaffold_max) == 0 || is.na(scaffold_max)) {
    warning("No scaffold length found for ", chrom_i, "; skipping ", outfile)
    return(invisible(data.table()))
  }

  peak_center <- floor((model_row$start[1] + model_row$end[1]) / 2)
  panel_start <- max(0L, as.integer(peak_center - 1000000L))
  panel_end <- panel_start + 2000000L
  if (panel_end > scaffold_max) {
    panel_end <- as.integer(scaffold_max)
    panel_start <- max(0L, as.integer(panel_end - 2000000L))
  }

  panel_dt <- manhattan_native_values[
    chromosome == chrom_i &
      x >= panel_start &
      x <= panel_end &
      stat %in% figure_stats
  ]

  panel_windows <- manhattan_window_reference[
    chromosome == chrom_i &
      window_start >= panel_start &
      window_end <= panel_end
  ]

  focal_candidates <- merge(
    peak_window_lookup_dt[as.character(peak_id) == as.character(peak_id_to_plot)],
    panel_windows,
    by = c("chromosome", "window_start", "window_end"),
    all.x = FALSE,
    all.y = FALSE,
    sort = FALSE
  )
  focal_candidates <- focal_candidates[!is.na(Xin_Bel_Fst)]

  if (nrow(focal_candidates) == 0) {
    warning("No focal Fst outlier window found for peak ", peak_id_to_plot, "; skipping ", outfile)
    return(invisible(data.table()))
  }

  if (identical(focal_rule, "closest_to_peak_center")) {
    focal_candidates[, distance_to_peak_center := abs(window_midpoint - peak_center)]
    setorder(focal_candidates, distance_to_peak_center)
    focal_window <- focal_candidates[1]
  } else {
    focal_window <- focal_candidates[which.max(Xin_Bel_Fst)]
  }
  focal_mid <- focal_window$window_midpoint[1]

  arg_data_stats <- c(
    "Xin_Bel_Fst_ARG_based",
    "xingu_enrich",
    "belem_enrich",
    "xingu_RTH_inverse",
    "belem_RTH_inverse",
    "Xin_Bel_CC_original",
    "Tap_Xin_RTH_original_inverse"
  )

  candidate_windows <- panel_windows[
    !is.na(Xin_Bel_Fst) &
      is_Xin_Bel_Fst_outlier == FALSE &
      window_midpoint > model_row$end[1]
  ]

  target_right <- min(panel_end - 1, model_row$end[1] + 750000)
  nonfocal_window <- choose_nonfocal_window(
    candidate_windows,
    target_right,
    peaks_dt,
    peak_id_to_plot,
    native_dt = panel_dt,
    arg_data_stats = arg_data_stats,
    exclude_window_midpoints = focal_mid
  )

  if (nrow(nonfocal_window) == 0) {
    candidate_windows <- panel_windows[
      !is.na(Xin_Bel_Fst) &
        is_Xin_Bel_Fst_outlier == FALSE &
        window_midpoint > focal_mid
    ]
    nonfocal_window <- choose_nonfocal_window(
      candidate_windows,
      target_right,
      peaks_dt,
      peak_id_to_plot,
      native_dt = panel_dt,
      arg_data_stats = arg_data_stats,
      exclude_window_midpoints = focal_mid
    )
  }

  nonfocal_mid <- if (nrow(nonfocal_window) > 0) nonfocal_window$window_midpoint[1] else NA_real_

  tree_plot_dt <- rbindlist(list(
    sample_trees_from_window(
      chromosome = chrom_i,
      window_start = focal_window$window_start[1],
      window_end = focal_window$window_end[1],
      window_label = "focal window",
      n_trees = 1L,
      tree_map = tree_map,
      tree_base_dir = tree_base_dir,
      iteration = 2000L
    ),
    if (nrow(nonfocal_window) > 0) {
      sample_trees_from_window(
        chromosome = chrom_i,
        window_start = nonfocal_window$window_start[1],
        window_end = nonfocal_window$window_end[1],
        window_label = "non-focal window",
        n_trees = 1L,
        tree_map = tree_map,
        tree_base_dir = tree_base_dir,
        iteration = 2000L
      )
    } else {
      data.table()
    }
  ), fill = TRUE)

  if (nrow(tree_plot_dt) > 0) {
    tree_plot_dt[, tree_panel_label := paste0(
      tree_role,
      "\nrandom tree sampled from ",
      n_local_trees_in_window,
      " local coalescence trees in window",
      "\n",
      chromosome,
      ":",
      tree_position
    )]
  }

  multi_iter_panel_dt <- NULL
  if (isTRUE(multi_iter)) {
    multi_iter_panel_dt <- read_multi_iter_values_for_panel(
      chromosome_i = chrom_i,
      panel_start = panel_start,
      panel_end = panel_end,
      iterations = multi_iter_iterations,
      arg_stat_dir = multi_iter_arg_stat_dir,
      arg_cc_stat_dir = multi_iter_arg_cc_stat_dir,
      argweaver_mask = multi_iter_argweaver_mask
    )
  }

  if (is.null(multi_iter_colors)) {
    multi_iter_colors <- c(
      "1600" = "#FE6100",
      "1700" = "#FFB000",
      "1800" = "#648FFF",
      "1900" = "#785EF0",
      "2000" = "#DC267F"
    )
  }

  n_stats <- length(figure_stats)
  right_ids <- rep(n_stats + 1L, n_stats)
  right_ids[seq.int(floor(n_stats / 2) + 1L, n_stats)] <- n_stats + 2L

  pdf(outfile, width = 7.3, height = max(8.3, 0.62 * n_stats + 4.1), useDingbats = FALSE, compress = TRUE)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  })

  layout(
    cbind(seq_len(n_stats), right_ids),
    widths = c(3.15, 1.75),
    heights = rep(1, n_stats)
  )
  # Keep all Manhattan plotting rows exactly the same height. The x-axis is
  # drawn directly on the lowest Manhattan subplot to keep tick marks/labels
  # aligned with that panel, while the shared bottom outer margin gives the
  # labels enough room.
  par(oma = c(2.4, 0, 4.7, 0))

  vertical_positions <- c(focal_mid, nonfocal_mid)
  vertical_labels <- c("focal window", "non-focal window")
  xlim_mb <- c(panel_start, panel_end) / 1e6

  for (i in seq_along(figure_stats)) {
    current_stat <- figure_stats[i]

    rolling_half_window_i <- 5L
    if (!is.null(rolling_half_window_overrides) &&
        current_stat %in% names(rolling_half_window_overrides)) {
      rolling_half_window_i <- as.integer(rolling_half_window_overrides[[current_stat]])
    }

    tmp <- panel_dt[stat == current_stat & !is.na(value)]
    thresh <- thresholds_dt[stat == current_stat & !is.na(threshold)]

    y_values <- tmp$value
    if (isTRUE(multi_iter) && current_stat %in% multi_iter_stats && !is.null(multi_iter_panel_dt)) {
      y_values <- multi_iter_panel_dt[stat == current_stat & !is.na(value), value]
    }
    y_values <- c(y_values, thresh$threshold)
    ylim <- range(y_values, na.rm = TRUE)
    if (!all(is.finite(ylim)) || diff(ylim) == 0) {
      ylim <- c(0, 1)
    } else {
      pad <- diff(ylim) * 0.08
      ylim <- c(ylim[1] - pad, ylim[2] + pad)
    }

    par(mar = c(1.35, 6.25, 0.25, 0.5), xpd = FALSE)
    plot(
      NA,
      xlim = xlim_mb,
      ylim = ylim,
      xaxt = "n",
      xlab = "",
      ylab = wrap_axis_label(unname(plot_labels[current_stat]), width = 18),
      cex.lab = 0.68,
      cex.axis = 0.62
    )

    if (nrow(thresh) > 0) {
      abline(h = thresh$threshold, lty = 2, lwd = 0.8, col = "grey35")
    }

    for (v_i in seq_along(vertical_positions)) {
      if (!is.na(vertical_positions[v_i])) {
        abline(v = vertical_positions[v_i] / 1e6, lty = 2, col = "grey40", lwd = 0.8)
      }
    }

    if (isTRUE(multi_iter) && current_stat %in% multi_iter_stats && !is.null(multi_iter_panel_dt)) {
      multi_tmp <- multi_iter_panel_dt[stat == current_stat & !is.na(value)]
      if (nrow(multi_tmp) > 0) {
        for (iter_i in sort(unique(multi_tmp$iteration))) {
          iter_dt <- multi_tmp[iteration == iter_i]
          setorder(iter_dt, x)
          lines(
            iter_dt$x / 1e6,
            rolling_mean_centered(iter_dt$value, half_window = rolling_half_window_i),
            col = multi_iter_colors[as.character(iter_i)],
            lwd = 1.0
          )
        }
        legend(
          "topleft",
          legend = names(multi_iter_colors),
          col = unname(multi_iter_colors),
          lwd = 1.0,
          bty = "n",
          cex = 0.45
        )
      }
    } else {
      if (nrow(tmp) > 0) {
        setorder(tmp, x)
        point_cex <- if (identical(current_stat, "Xin_Bel_Fst")) 0.78 else 0.42
        points(tmp$x / 1e6, tmp$value, pch = 16, cex = point_cex, col = "black")
        lines(
          tmp$x / 1e6,
          rolling_mean_centered(tmp$value, half_window = rolling_half_window_i),
          col = "red",
          lwd = 0.75
        )
      }
    }

    if (i == n_stats) {
      x_ticks <- seq(
        floor(xlim_mb[1] * 4) / 4,
        ceiling(xlim_mb[2] * 4) / 4,
        by = 0.25
      )
      x_ticks <- x_ticks[x_ticks >= xlim_mb[1] & x_ticks <= xlim_mb[2]]
      axis(
        1,
        at = x_ticks,
        labels = format(x_ticks, trim = TRUE),
        cex.axis = 0.62,
        line = 0.05
      )
      mtext(
        format_scaffold_axis_label(chrom_i),
        side = 1,
        line = 1.55,
        cex = 0.68
      )
    }

    if (i == 1L) {
      usr <- par("usr")
      par(xpd = NA)
      text(vertical_positions / 1e6, usr[4] + 0.04 * diff(usr[3:4]), labels = vertical_labels, cex = 0.50, xpd = NA)
      par(xpd = FALSE)
      title_txt <- if (is.null(figure_title)) paste0("Peak ", peak_id_to_plot) else figure_title
      mtext(
        paste0(title_txt, ": ", chrom_i, ":", model_row$start[1], "-", model_row$end[1]),
        side = 3,
        outer = TRUE,
        line = 2.7,
        cex = 0.85,
        font = 1
      )
      mtext(
        paste0("Model: ", model_label_without_redundant_model(model_row$model_assignment[1])),
        side = 3,
        outer = TRUE,
        line = 1.6,
        cex = 0.75,
        font = 1
      )
    }

    box()
  }

  for (tree_i in seq_len(2L)) {
    par(mar = c(1.0, 0.7, 3.2, 0.7), xpd = NA)
    if (nrow(tree_plot_dt) >= tree_i) {
      plot_single_tree_panel(tree_plot_dt[tree_i])
    } else {
      plot.new()
      text(0.5, 0.5, "No tree sampled", cex = 0.7)
    }
  }

  invisible(data.table(
    peak_id = peak_id_to_plot,
    chromosome = chrom_i,
    peak_start = model_row$start[1],
    peak_end = model_row$end[1],
    panel_start = panel_start,
    panel_end = panel_end,
    focal_window_midpoint = focal_mid,
    nonfocal_window_midpoint = nonfocal_mid,
    output_pdf = outfile
  ))
}

# =============================================================
# Manhattan panels for Xin-Bel Fst peaks with ARG-based data
# =============================================================
# Each panel spans a 2-Mb region where possible. For Manhattan plots, values
# are plotted at their imported/native resolution after the same quality
# filtering already used above. RTH statistics are plotted as inverses, and
# selscan iHS/nSL values are kept signed.

manhattan_stats <- c(
  "Xin_Bel_Fst",
  "Xin_Bel_Fst_ARG_based",
  "Xin_Bel_Dxy",
  "Xin_pi",
  "Bel_pi",
  "xingu_enrich",
  "belem_enrich",
  "xingu_RTH_inverse",
  "belem_RTH_inverse",
  "Xin_Bel_CC_original",
  "Tap_Xin_RTH_original_inverse",
  "Tapajos_Xingu_fd",
  "recombination_rate",
  "Xingu_RAiSD_u",
  "xingu_norm_ihs",
  "xingu_norm_nsl",
  "Belem_RAiSD_u",
  "belem_norm_ihs",
  "belem_norm_nsl",
  "dup_sites",
  "perc_dup_sites"
)

make_native_stat_dt <- function(dt, stat, x_col = "x") {
  data.table(
    chromosome = dt$chromosome,
    x = as.numeric(dt[[x_col]]),
    stat = stat,
    value = as.numeric(dt[[stat]])
  )
}

manhattan_native_list <- list()

# 10-kb pixy-window statistics are plotted at their window midpoints.
manhattan_window_reference <- df_clean[, .(
  chromosome,
  window_start = as.integer(start),
  window_end = as.integer(end),
  window_id = paste(chromosome, start, end, sep = ":"),
  window_midpoint = (start + end) / 2,
  Xin_Bel_Fst,
  is_Xin_Bel_Fst_outlier = !is.na(Xin_Bel_Fst) & Xin_Bel_Fst > threshold,
  Xin_Bel_Dxy,
  Xin_pi,
  Bel_pi
)]

for (stat in c("Xin_Bel_Fst", "Xin_Bel_Dxy", "Xin_pi", "Bel_pi")) {
  tmp <- manhattan_window_reference[, .(
    chromosome,
    x = window_midpoint,
    stat = stat,
    value = as.numeric(get(stat))
  )]
  manhattan_native_list[[length(manhattan_native_list) + 1L]] <- tmp
}

# D statistics: already window-based but 1-based; plot imported windows at
# their midpoint after sitesUsed filtering and D<0 -> fd=0 conversion.
d_stats_manhattan <- fread(
  d_stats_file,
  select = c("scaffold", "start", "end", "sitesUsed", "D", "fd")
)
setnames(
  d_stats_manhattan,
  old = c("scaffold", "fd"),
  new = c("chromosome", "Tapajos_Xingu_fd")
)
d_stats_manhattan[!is.na(D) & D < 0, Tapajos_Xingu_fd := 0]
# Exclude biologically invalid fd values before Manhattan/native-value plotting.
d_stats_manhattan <- d_stats_manhattan[
  is.na(Tapajos_Xingu_fd) |
    (Tapajos_Xingu_fd >= 0 & Tapajos_Xingu_fd <= 1)
]
d_stats_manhattan <- d_stats_manhattan[sitesUsed >= 100]
d_stats_manhattan[, x := (as.numeric(start - 1L) + as.numeric(end)) / 2]
manhattan_native_list[[length(manhattan_native_list) + 1L]] <-
  make_native_stat_dt(d_stats_manhattan, "Tapajos_Xingu_fd")
rm(d_stats_manhattan)
gc()

# ReLERNN recombination rate: plot imported interval midpoints.
relernn_manhattan <- fread(
  relernn_file,
  header = FALSE,
  col.names = c("chromosome", "interval_start", "interval_end", "recombination_rate")
)
relernn_manhattan[, x := (as.numeric(interval_start) + as.numeric(interval_end)) / 2]
manhattan_native_list[[length(manhattan_native_list) + 1L]] <-
  make_native_stat_dt(relernn_manhattan, "recombination_rate")
rm(relernn_manhattan)
gc()

# RAiSD U statistics: plot imported interval midpoints.
xingu_raisd_manhattan <- fread(
  xingu_raisd_file,
  header = FALSE,
  col.names = c("chromosome", "interval_start", "interval_end", "Xingu_RAiSD_u")
)
xingu_raisd_manhattan[, x := (as.numeric(interval_start) + as.numeric(interval_end)) / 2]
manhattan_native_list[[length(manhattan_native_list) + 1L]] <-
  make_native_stat_dt(xingu_raisd_manhattan, "Xingu_RAiSD_u")
rm(xingu_raisd_manhattan)
gc()

belem_raisd_manhattan <- fread(
  belem_raisd_file,
  header = FALSE,
  col.names = c("chromosome", "interval_start", "interval_end", "Belem_RAiSD_u")
)
belem_raisd_manhattan[, x := (as.numeric(interval_start) + as.numeric(interval_end)) / 2]
manhattan_native_list[[length(manhattan_native_list) + 1L]] <-
  make_native_stat_dt(belem_raisd_manhattan, "Belem_RAiSD_u")
rm(belem_raisd_manhattan)
gc()

# rCNV 10-kb windows: plot imported windows at their midpoint.
rcnv_manhattan <- fread(
  rcnv_file,
  select = c("scaffold", "start", "end", "total.sites", "dup.sites", "perc.dup.sites"),
  na.strings = c("NA", "")
)
setnames(
  rcnv_manhattan,
  old = c("scaffold", "start", "end", "total.sites", "dup.sites", "perc.dup.sites"),
  new = c("chromosome", "interval_start", "interval_end", "total_sites", "dup_sites", "perc_dup_sites")
)
rcnv_manhattan[, chromosome := standardize_scaffold_names(chromosome)]
rcnv_manhattan <- rcnv_manhattan[!is.na(total_sites) & total_sites >= 50]
rcnv_manhattan[, x := (as.numeric(interval_start) + as.numeric(interval_end)) / 2]
manhattan_native_list[[length(manhattan_native_list) + 1L]] <-
  make_native_stat_dt(rcnv_manhattan, "dup_sites")
manhattan_native_list[[length(manhattan_native_list) + 1L]] <-
  make_native_stat_dt(rcnv_manhattan, "perc_dup_sites")
rm(rcnv_manhattan)
gc()

# ARG-based branch Fst point statistic: plot imported point values.
arg_based_fst_manhattan <- fread(
  arg_based_fst_file,
  select = c(
    "region",
    "coordinate",
    "Fst_belem_tapajos",
    "Fst_belem_xingu",
    "Fst_tapajos_xingu"
  )
)
setnames(
  arg_based_fst_manhattan,
  old = c(
    "region",
    "coordinate",
    "Fst_tapajos_xingu",
    "Fst_belem_tapajos",
    "Fst_belem_xingu"
  ),
  new = c(
    "chromosome",
    "pos",
    "Tap_Xin_Fst_ARG_based",
    "Tap_Bel_Fst_ARG_based",
    "Xin_Bel_Fst_ARG_based"
  )
)
arg_based_fst_manhattan[, chromosome := standardize_scaffold_names(chromosome)]
arg_based_fst_manhattan <- filter_points_outside_mask(arg_based_fst_manhattan, argweaver_mask, pos_col = "pos")
arg_based_fst_manhattan[, x := as.numeric(pos)]
manhattan_native_list[[length(manhattan_native_list) + 1L]] <-
  make_native_stat_dt(arg_based_fst_manhattan, "Xin_Bel_Fst_ARG_based")
rm(arg_based_fst_manhattan)
gc()

# ARG local-tree point statistics. RTH statistics are inverted for plotting.
arg_manhattan <- fread(
  arg_stat_file_2000,
  select = c(
    "chrom",
    "pos",
    "Tap_Xin_RTH_original",
    "belem_RTH",
    "xingu_RTH",
    "belem_enrich",
    "xingu_enrich"
  )
)
arg_manhattan[, chromosome := standardize_scaffold_names(chrom)]
arg_manhattan[, chrom := NULL]
arg_manhattan <- filter_points_outside_mask(arg_manhattan, argweaver_mask, pos_col = "pos")
arg_manhattan[, xingu_RTH_inverse := safe_inverse(xingu_RTH)]
arg_manhattan[, belem_RTH_inverse := safe_inverse(belem_RTH)]
arg_manhattan[, Tap_Xin_RTH_original_inverse := safe_inverse(Tap_Xin_RTH_original)]
arg_manhattan[, x := as.numeric(pos)]

for (stat in c(
  "xingu_enrich",
  "belem_enrich",
  "xingu_RTH_inverse",
  "belem_RTH_inverse",
  "Tap_Xin_RTH_original_inverse"
)) {
  manhattan_native_list[[length(manhattan_native_list) + 1L]] <-
    make_native_stat_dt(arg_manhattan, stat)
}
rm(arg_manhattan)
gc()

# ARG cross-coalescence point statistic.
arg_cc_manhattan <- fread(
  arg_cc_stat_file_2000,
  select = c("chrom", "pos", "Xin_Bel_CC_original")
)
arg_cc_manhattan[, chromosome := standardize_scaffold_names(chrom)]
arg_cc_manhattan[, chrom := NULL]
arg_cc_manhattan <- filter_points_outside_mask(arg_cc_manhattan, argweaver_mask, pos_col = "pos")
arg_cc_manhattan[, x := as.numeric(pos)]
manhattan_native_list[[length(manhattan_native_list) + 1L]] <-
  make_native_stat_dt(arg_cc_manhattan, "Xin_Bel_CC_original")
rm(arg_cc_manhattan)
gc()

# selscan point statistics. Use absolute values for Manhattan plots.
xingu_ihs_manhattan <- fread(xingu_ihs_file, select = c("chr", "pos", "norm_ihs"))
setnames(xingu_ihs_manhattan, c("chr", "norm_ihs"), c("chromosome", "xingu_norm_ihs"))
xingu_ihs_manhattan[, xingu_norm_ihs := abs(xingu_norm_ihs)]
xingu_ihs_manhattan[, x := as.numeric(pos)]
manhattan_native_list[[length(manhattan_native_list) + 1L]] <-
  make_native_stat_dt(xingu_ihs_manhattan, "xingu_norm_ihs")
rm(xingu_ihs_manhattan)
gc()

xingu_nsl_manhattan <- fread(xingu_nsl_file, select = c("chr", "pos", "norm_nsl"))
setnames(xingu_nsl_manhattan, c("chr", "norm_nsl"), c("chromosome", "xingu_norm_nsl"))
xingu_nsl_manhattan[, xingu_norm_nsl := abs(xingu_norm_nsl)]
xingu_nsl_manhattan[, x := as.numeric(pos)]
manhattan_native_list[[length(manhattan_native_list) + 1L]] <-
  make_native_stat_dt(xingu_nsl_manhattan, "xingu_norm_nsl")
rm(xingu_nsl_manhattan)
gc()

belem_ihs_manhattan <- fread(belem_ihs_file, select = c("chr", "pos", "norm_ihs"))
setnames(belem_ihs_manhattan, c("chr", "norm_ihs"), c("chromosome", "belem_norm_ihs"))
belem_ihs_manhattan[, belem_norm_ihs := abs(belem_norm_ihs)]
belem_ihs_manhattan[, x := as.numeric(pos)]
manhattan_native_list[[length(manhattan_native_list) + 1L]] <-
  make_native_stat_dt(belem_ihs_manhattan, "belem_norm_ihs")
rm(belem_ihs_manhattan)
gc()

belem_nsl_manhattan <- fread(belem_nsl_file, select = c("chr", "pos", "norm_nsl"))
setnames(belem_nsl_manhattan, c("chr", "norm_nsl"), c("chromosome", "belem_norm_nsl"))
belem_nsl_manhattan[, belem_norm_nsl := abs(belem_norm_nsl)]
belem_nsl_manhattan[, x := as.numeric(pos)]
manhattan_native_list[[length(manhattan_native_list) + 1L]] <-
  make_native_stat_dt(belem_nsl_manhattan, "belem_norm_nsl")
rm(belem_nsl_manhattan)
gc()

manhattan_native_values <- rbindlist(manhattan_native_list, fill = TRUE)
rm(manhattan_native_list)
gc()

multi_iter_ARG_iterations <- c(1600L, 1700L, 1800L, 1900L, 2000L)

multi_iter_ARG_colors <- c(
  "1600" = "#FE6100",
  "1700" = "#FFB000",
  "1800" = "#648FFF",
  "1900" = "#785EF0",
  "2000" = "#DC267F"
)

multi_iter_manhattan_stats <- c(
  "xingu_enrich",
  "belem_enrich",
  "xingu_RTH_inverse",
  "belem_RTH_inverse",
  "Xin_Bel_CC_original",
  "Tap_Xin_RTH_original_inverse"
)

# Do not build a genome-wide multi-iteration Manhattan object here.
# Multi-iteration ARG values are read lazily inside each 2-Mb panel to keep
# memory use low on machines with a 24-GB R vector memory limit.

# Thresholds to draw on Manhattan panels. Each subplot receives only the
# threshold row(s) for its own statistic. Xin_Bel_Fst and ARG-based Xin_Bel_Fst
# use mean + 5 SD rather than the empirical control-window thresholds.
arg_based_fst_threshold_values <- manhattan_native_values[
  stat == "Xin_Bel_Fst_ARG_based" & !is.na(value),
  value
]

manhattan_threshold_rows <- rbindlist(list(
  data.table(
    stat = "Xin_Bel_Fst",
    tail = "upper",
    threshold = mean(
      df_clean$Xin_Bel_Fst,
      na.rm = TRUE
    ) + 5 * sd(df_clean$Xin_Bel_Fst, na.rm = TRUE)
  ),
  data.table(
    stat = "Xin_Bel_Fst_ARG_based",
    tail = "upper",
    threshold = mean(arg_based_fst_threshold_values, na.rm = TRUE) +
      5 * sd(arg_based_fst_threshold_values, na.rm = TRUE)
  ),
  get_threshold_for_manhattan("Xin_Bel_Dxy", pixy_threshold_summary),
  get_threshold_for_manhattan("Xin_pi", pixy_threshold_summary),
  get_threshold_for_manhattan("Bel_pi", pixy_threshold_summary),
  get_threshold_for_manhattan("xingu_enrich", arg_threshold_summary),
  get_threshold_for_manhattan("belem_enrich", arg_threshold_summary),
  get_threshold_for_manhattan("xingu_RTH_inverse", arg_threshold_summary),
  get_threshold_for_manhattan("belem_RTH_inverse", arg_threshold_summary),
  get_threshold_for_manhattan("Xin_Bel_CC_original", arg_cc_threshold_summary),
  data.table(
    stat = "Tap_Xin_RTH_original_inverse",
    tail = "upper",
    threshold = as.numeric(quantile(
      safe_inverse(arg_joined[window_class == "control"]$Tap_Xin_RTH_original),
      probs = 0.999,
      na.rm = TRUE
    ))
  ),
  get_threshold_for_manhattan("Tapajos_Xingu_fd", additional_threshold_summary),
  get_threshold_for_manhattan("recombination_rate", additional_threshold_summary),
  get_threshold_for_manhattan("Xingu_RAiSD_u", additional_threshold_summary),
  data.table(stat = "xingu_norm_ihs", tail = "upper", threshold = 2),
  data.table(stat = "xingu_norm_nsl", tail = "upper", threshold = 2),
  get_threshold_for_manhattan("Belem_RAiSD_u", additional_threshold_summary),
  data.table(stat = "belem_norm_ihs", tail = "upper", threshold = 2),
  data.table(stat = "belem_norm_nsl", tail = "upper", threshold = 2),
  get_threshold_for_manhattan("dup_sites", additional_threshold_summary),
  get_threshold_for_manhattan("perc_dup_sites", additional_threshold_summary)
), fill = TRUE)

scaffold_lengths <- df_clean[, .(
  scaffold_max_end = max(end, na.rm = TRUE)
), by = chromosome]

# Always produce a supplementary stacked Manhattan figure for repeat-adjacent
# peaks that also have ARG-based data. Only Xin_Bel_Fst is plotted; each panel
# has its own scaffold-specific x axis, and no trees, model text, or footnotes.
large_repeat_adjacent_Fst_manhattan_summary <- data.table()
if (isTRUE(CREATE_FIGURES_PDFS)) {
  large_repeat_adjacent_Fst_manhattan_summary <-
    make_large_repeat_adjacent_Fst_manhattan_figure(
      repeat_peak_dt = large_repeat_peak_proximity,
      arg_peak_coverage_dt = arg_peak_ARG_coverage,
      fst_windows_dt = df_clean,
      repeat_regions = repeat_regions_for_manhattan,
      scaffold_lengths_dt = scaffold_lengths,
      fst_threshold = threshold,
      outfile = figure_large_repeat_adjacent_Fst_manhattan_pdf,
      panel_width_bp = 2000000L
    )
}

tree_folder_coordinate_map <- read_tree_folder_coordinate_map(
  tree_folder_coordinate_map_file
)

set.seed(tree_sampling_seed)

# =============================================================
# Manuscript figures
# =============================================================
figure_4_summary <- data.table()
figure_5_summary <- data.table()

if (isTRUE(CREATE_FIGURES_PDFS)) {
make_model_filtered_heatmap_pdf(
  heatmap_dt = heatmap_peak_representative,
  model_dt = model_assignment_any,
  stat_order = heatmap_stat_order,
  plot_labels = plot_labels,
  outfile = figure_3_1_pdf,
  title_text = "Figure 3.1. Peak percentile profiles by model assignment"
)

make_model_filtered_heatmap_pdf(
  heatmap_dt = heatmap_peak_representative,
  model_dt = model_assignment_min25,
  stat_order = heatmap_stat_order,
  plot_labels = plot_labels,
  outfile = figure_3_2_pdf,
  title_text = "Figure 3.2. Peak percentile profiles by model assignment"
)

figure_4_stats <- c(
  "Xin_Bel_Fst",
  "xingu_enrich",
  "belem_enrich",
  "xingu_RTH_inverse",
  "belem_RTH_inverse",
  "Xin_Bel_CC_original",
  "Tap_Xin_RTH_original_inverse"
)

figure_5_stats <- c(
  "Xin_Bel_Fst",
  "Xin_Bel_Fst_ARG_based",
  "xingu_enrich",
  "belem_enrich",
  "xingu_RTH_inverse",
  "belem_RTH_inverse",
  "Xin_Bel_CC_original",
  "Tap_Xin_RTH_original_inverse",
  "Belem_RAiSD_u",
  "belem_norm_nsl"
)

figure_4_summary <- make_oneoff_manhattan_tree_figure(
  peak_id_to_plot = "7",
  outfile = figure_4_pdf,
  figure_stats = figure_4_stats,
  model_dt = model_assignment_any,
  manhattan_native_values = manhattan_native_values,
  manhattan_window_reference = manhattan_window_reference,
  thresholds_dt = manhattan_threshold_rows,
  peaks_dt = peaks,
  plot_labels = plot_labels,
  scaffold_lengths = scaffold_lengths,
  peak_window_lookup_dt = peak_window_lookup,
  tree_map = tree_folder_coordinate_map,
  tree_base_dir = arg_tree_files_midpoint_dir,
  multi_iter = TRUE,
  multi_iter_stats = multi_iter_manhattan_stats,
  multi_iter_colors = multi_iter_ARG_colors,
  multi_iter_iterations = multi_iter_ARG_iterations,
  multi_iter_arg_stat_dir = arg_stat_dir,
  multi_iter_arg_cc_stat_dir = arg_cc_stat_dir,
  multi_iter_argweaver_mask = argweaver_mask,
  focal_rule = "closest_to_peak_center",
  figure_title = "Figure 4"
)

figure_5_summary <- make_oneoff_manhattan_tree_figure(
  peak_id_to_plot = "51",
  outfile = figure_5_pdf,
  figure_stats = figure_5_stats,
  model_dt = model_assignment_any,
  manhattan_native_values = manhattan_native_values,
  manhattan_window_reference = manhattan_window_reference,
  thresholds_dt = manhattan_threshold_rows,
  peaks_dt = peaks,
  plot_labels = plot_labels,
  scaffold_lengths = scaffold_lengths,
  peak_window_lookup_dt = peak_window_lookup,
  tree_map = tree_folder_coordinate_map,
  tree_base_dir = arg_tree_files_midpoint_dir,
  multi_iter = FALSE,
  focal_rule = "max_fst",
  figure_title = "Figure 5",
  rolling_half_window_overrides = c(
    "Xin_Bel_Fst_ARG_based" = 20L,
    "xingu_enrich" = 20L,
    "belem_enrich" = 20L,
    "xingu_RTH_inverse" = 20L,
    "belem_RTH_inverse" = 20L,
    "Xin_Bel_CC_original" = 20L,
    "Tap_Xin_RTH_original_inverse" = 20L
  )
)

fwrite(figure_4_summary, file.path(figures_dir, "Figure_4_summary.tsv"), sep = "\t")
fwrite(figure_5_summary, file.path(figures_dir, "Figure_5_summary.tsv"), sep = "\t")
}


manhattan_any_summary <- data.table()
if (isTRUE(CREATE_MANHATTAN_ANY_SINGLE_ITER_PDFS)) {
manhattan_any_summary <- create_manhattan_panels_for_assignments(
  model_dt = model_assignment_any,
  assignment_dir = model_assignment_any_dir,
  manhattan_native_values = manhattan_native_values,
  manhattan_window_reference = manhattan_window_reference,
  thresholds_dt = manhattan_threshold_rows,
  peaks_dt = peaks,
  plot_labels = plot_labels,
  manhattan_stats = manhattan_stats,
  scaffold_lengths = scaffold_lengths,
  peak_window_lookup_dt = peak_window_lookup,
  tree_map = tree_folder_coordinate_map,
  tree_base_dir = arg_tree_files_midpoint_dir,
  tree_iteration = 2000L,
  footnote_pdf_file = manhattan_footnotes_pdf,
  repeat_regions = repeat_regions_for_manhattan
)

fwrite(
  manhattan_any_summary,
  manhattan_any_summary_file,
  sep = "\t"
)

}

manhattan_min25_summary <- data.table()
if (isTRUE(CREATE_MANHATTAN_MIN25_SINGLE_ITER_PDFS)) {
manhattan_min25_summary <- create_manhattan_panels_for_assignments(
  model_dt = model_assignment_min25,
  assignment_dir = model_assignment_min25_dir,
  manhattan_native_values = manhattan_native_values,
  manhattan_window_reference = manhattan_window_reference,
  thresholds_dt = manhattan_threshold_rows,
  peaks_dt = peaks,
  plot_labels = plot_labels,
  manhattan_stats = manhattan_stats,
  scaffold_lengths = scaffold_lengths,
  peak_window_lookup_dt = peak_window_lookup,
  tree_map = tree_folder_coordinate_map,
  tree_base_dir = arg_tree_files_midpoint_dir,
  tree_iteration = 2000L,
  footnote_pdf_file = manhattan_footnotes_pdf,
  repeat_regions = repeat_regions_for_manhattan
)

fwrite(
  manhattan_min25_summary,
  manhattan_min25_summary_file,
  sep = "	"
)

}

manhattan_any_multi_iter_summary <- data.table()
if (isTRUE(CREATE_MANHATTAN_ANY_MULTI_ITER_PDFS)) {
manhattan_any_multi_iter_summary <- create_manhattan_panels_for_assignments(
  model_dt = model_assignment_any,
  assignment_dir = model_assignment_any_dir,
  manhattan_native_values = manhattan_native_values,
  manhattan_window_reference = manhattan_window_reference,
  thresholds_dt = manhattan_threshold_rows,
  peaks_dt = peaks,
  plot_labels = plot_labels,
  manhattan_stats = manhattan_stats,
  scaffold_lengths = scaffold_lengths,
  peak_window_lookup_dt = peak_window_lookup,
  tree_map = tree_folder_coordinate_map,
  tree_base_dir = arg_tree_files_midpoint_dir,
  tree_iteration = 2000L,
  footnote_pdf_file = manhattan_footnotes_pdf,
  panel_subdir = "Manhattan_panels_multi_iter",
  multi_iter_stats = multi_iter_manhattan_stats,
  multi_iter_colors = multi_iter_ARG_colors,
  multi_iter_iterations = multi_iter_ARG_iterations,
  multi_iter_arg_stat_dir = arg_stat_dir,
  multi_iter_arg_cc_stat_dir = arg_cc_stat_dir,
  multi_iter_argweaver_mask = argweaver_mask,
  repeat_regions = repeat_regions_for_manhattan
)

fwrite(
  manhattan_any_multi_iter_summary,
  manhattan_any_multi_iter_summary_file,
  sep = "	"
)

}

manhattan_min25_multi_iter_summary <- data.table()
if (isTRUE(CREATE_MANHATTAN_MIN25_MULTI_ITER_PDFS)) {
manhattan_min25_multi_iter_summary <- create_manhattan_panels_for_assignments(
  model_dt = model_assignment_min25,
  assignment_dir = model_assignment_min25_dir,
  manhattan_native_values = manhattan_native_values,
  manhattan_window_reference = manhattan_window_reference,
  thresholds_dt = manhattan_threshold_rows,
  peaks_dt = peaks,
  plot_labels = plot_labels,
  manhattan_stats = manhattan_stats,
  scaffold_lengths = scaffold_lengths,
  peak_window_lookup_dt = peak_window_lookup,
  tree_map = tree_folder_coordinate_map,
  tree_base_dir = arg_tree_files_midpoint_dir,
  tree_iteration = 2000L,
  footnote_pdf_file = manhattan_footnotes_pdf,
  panel_subdir = "Manhattan_panels_multi_iter",
  multi_iter_stats = multi_iter_manhattan_stats,
  multi_iter_colors = multi_iter_ARG_colors,
  multi_iter_iterations = multi_iter_ARG_iterations,
  multi_iter_arg_stat_dir = arg_stat_dir,
  multi_iter_arg_cc_stat_dir = arg_cc_stat_dir,
  multi_iter_argweaver_mask = argweaver_mask,
  repeat_regions = repeat_regions_for_manhattan
)

fwrite(
  manhattan_min25_multi_iter_summary,
  manhattan_min25_multi_iter_summary_file,
  sep = "	"
)

}

# Add model assignments from all 50 saved MCMC iterations to the two
# multi-iteration Manhattan summary tables only when the slow genome-wide
# multi-iteration ARG summaries were explicitly requested. The Manhattan PDFs
# themselves remain independently controlled by their own switches.
if (isTRUE(RUN_MULTI_ITER_ARG_SUMMARIES)) {

manhattan_any_multi_iter_summary <- add_multi_iter_model_assignments_to_manhattan_summary(
  manhattan_summary_dt = manhattan_any_multi_iter_summary,
  model_dt_2000 = model_assignment_any,
  multi_iter_model_dt = multi_iter_model_assignments_all_iterations,
  assignment_threshold_type = "at_least_one_tree_signif_in_peak",
  iterations = all_ARG_MCMC_iterations,
  all_peaks_dt = peaks
)

fwrite(
  manhattan_any_multi_iter_summary,
  manhattan_any_multi_iter_summary_file,
  sep = "\t"
)

manhattan_min25_multi_iter_summary <- add_multi_iter_model_assignments_to_manhattan_summary(
  manhattan_summary_dt = manhattan_min25_multi_iter_summary,
  model_dt_2000 = model_assignment_min25,
  multi_iter_model_dt = multi_iter_model_assignments_all_iterations,
  assignment_threshold_type = "25perc_trees_signif_in_peak",
  iterations = all_ARG_MCMC_iterations,
  all_peaks_dt = peaks
)

fwrite(
  manhattan_min25_multi_iter_summary,
  manhattan_min25_multi_iter_summary_file,
  sep = "\t"
)

}

rm(
  manhattan_native_values,
  manhattan_window_reference,
  manhattan_threshold_rows,
  scaffold_lengths,
  arg_based_fst_threshold_values
)
gc()

# Print summary
# -----------------------------
cat("\nXin_Bel_Fst outlier summary\n")
cat("===========================\n")
cat("Input file: ", pixy_file, "\n", sep = "")
cat("Autosomal windows with >=50 SNPs: ", nrow(df_clean), "\n", sep = "")
cat("Mean Xin_Bel_Fst: ", mean_fst, "\n", sep = "")
cat("SD Xin_Bel_Fst: ", sd_fst, "\n", sep = "")
cat("Control threshold (mean + 1xSD): ", control_threshold, "\n", sep = "")
cat("Outlier threshold (mean + 5xSD): ", threshold, "\n", sep = "")
cat("Control windows: ", nrow(control_windows), "\n", sep = "")
cat("Intermediate windows (>1 SD and <=5 SD): ", nrow(intermediate_windows), "\n", sep = "")
cat("Outlier windows: ", nrow(outliers), "\n", sep = "")
cat("Number of outlier peaks: ", nrow(peaks), "\n", sep = "")

cat("\nWindow accounting\n")
cat("=================\n")
cat("Total filtered windows: ", nrow(df_clean), "\n", sep = "")
cat("Outlier windows: ", nrow(outliers), "\n", sep = "")
cat("Control windows: ", nrow(control_windows), "\n", sep = "")
cat("Intermediate windows: ", nrow(intermediate_windows), "\n", sep = "")
cat("NA Xin_Bel_Fst windows: ",
    sum(is.na(df_clean$Xin_Bel_Fst)), "\n", sep = "")

cat("\nOutlier vs control pixy statistic summary\n")
cat("=========================================\n")
print(window_stat_summary)

cat("\nPixy empirical thresholds from control windows\n")
cat("==============================================\n")
print(pixy_threshold_summary)

cat("\nAdditional empirical thresholds from control windows\n")
cat("====================================================\n")
print(additional_threshold_summary)

cat("\nARGweaver iteration 2000 summary\n")
cat("================================\n")
cat("ARG stat file: ", arg_stat_file_2000, "\n", sep = "")
cat("ARG local trees assigned to control/outlier windows: ",
    nrow(arg_joined), "\n", sep = "")
cat("ARG local trees assigned to control windows: ",
    nrow(arg_joined[window_class == "control"]), "\n", sep = "")
cat("ARG local trees assigned to outlier windows: ",
    nrow(arg_joined[window_class == "Xin_Bel_Fst_outlier"]), "\n", sep = "")

cat("\nARG empirical thresholds from control windows\n")
cat("=============================================\n")
print(arg_threshold_summary)

cat("\nTap_Xin_RTH values equal to 1\n")
cat("==============================\n")
print(tap_xin_rth_value_one_summary)

cat("\nTap_Xin_RTH decile distribution\n")
cat("===============================\n")
print(tap_xin_rth_decile_distribution)
cat("Bar plot PDF: ", tap_xin_rth_decile_barplot_pdf, "\n", sep = "")

cat("\nARGweaver CC iteration 2000 summary\n")
cat("===================================\n")
cat("ARG CC stat file: ", arg_cc_stat_file_2000, "\n", sep = "")
cat("ARG CC local trees assigned to control/outlier windows: ",
    nrow(arg_cc_joined), "\n", sep = "")
cat("ARG CC local trees assigned to control windows: ",
    nrow(arg_cc_joined[window_class == "control"]), "\n", sep = "")
cat("ARG CC local trees assigned to outlier windows: ",
    nrow(arg_cc_joined[window_class == "Xin_Bel_Fst_outlier"]), "\n", sep = "")

cat("\nARG CC empirical thresholds from control windows\n")
cat("================================================\n")
print(arg_cc_threshold_summary)

cat("\nARG CC outlier vs control statistic summary\n")
cat("============================================\n")
print(arg_cc_stat_summary)

cat("\nXin_Bel_CC values equal to 1\n")
cat("============================\n")
print(arg_cc_value_one_summary)

cat("\nTap_Xin_JCR_v2 vs Xin_Bel_CC_original control scatterplot\n")
cat("=========================================================\n")
cat("R^2: ", jcr_cc_r2, "\n", sep = "")
cat("PDF: ", arg_cc_jcr_rcc_scatter_pdf, "\n", sep = "")

cat("\nARG-based data coverage for Fst outlier windows and peaks\n")
cat("===========================================================\n")
print(arg_ARG_coverage_summary)

cat("\nSignificant ARG statistic window counts\n")
cat("=======================================\n")
print(arg_significant_window_counts)
cat("Bar plot PDF: ", arg_significant_window_barplot_pdf, "\n", sep = "")

cat("\nSignificant ARG statistic window counts, minimum 5% local trees\n")
cat("===============================================================\n")
print(arg_significant_window_counts_min5)
cat("Bar plot PDF: ", arg_significant_window_barplot_min5_pdf, "\n", sep = "")

cat("\nSignificant ARG statistic peak counts\n")
cat("=====================================\n")
print(arg_significant_peak_counts)
cat("Bar plot PDF: ", arg_significant_peak_barplot_pdf, "\n", sep = "")

cat("\nSignificant ARG statistic peak counts, minimum 5% local trees\n")
cat("=============================================================\n")
print(arg_significant_peak_counts_min5)
cat("Bar plot PDF: ", arg_significant_peak_barplot_min5_pdf, "\n", sep = "")

cat("\nXin-Bel Fst peak model assignments: at least one significant local tree\n")
cat("=====================================================================\n")
print(model_assignment_any_summary)

cat("\nXin-Bel Fst peak model assignments: >=25% significant local trees\n")
cat("===============================================================\n")
print(model_assignment_min25_summary)

cat("\nOutput files:\n")
cat("  ", outlier_windows_bed, "\n", sep = "")
cat("  ", outlier_peaks_bed, "\n", sep = "")
cat("  ", comparison_windows_file, "\n", sep = "")
cat("  ", window_stat_summary_file, "\n", sep = "")
cat("  ", region_set_pi_recombination_summary_file, "\n", sep = "")
cat("  ", large_repeat_peak_proximity_file, "\n", sep = "")
cat("  ", figure_large_repeat_adjacent_Fst_manhattan_pdf, "\n", sep = "")
cat("  ", pixy_threshold_summary_file, "\n", sep = "")
cat("  ", additional_threshold_summary_file, "\n", sep = "")
cat("  ", additional_control_value_summary_file, "\n", sep = "")
cat("  ", arg_threshold_summary_file, "\n", sep = "")
cat("  ", tap_xin_rth_value_one_summary_file, "\n", sep = "")
cat("  ", tap_xin_rth_decile_distribution_file, "\n", sep = "")
cat("  ", tap_xin_rth_decile_barplot_pdf, "\n", sep = "")
cat("  ", arg_cc_threshold_summary_file, "\n", sep = "")
cat("  ", arg_cc_stat_summary_file, "\n", sep = "")
cat("  ", arg_cc_value_one_summary_file, "\n", sep = "")
cat("  ", arg_cc_jcr_rcc_scatter_pdf, "\n", sep = "")
cat("  ", arg_significant_window_counts_file, "\n", sep = "")
cat("  ", arg_significant_window_barplot_pdf, "\n", sep = "")
cat("  ", arg_significant_window_counts_min5_file, "\n", sep = "")
cat("  ", arg_significant_window_barplot_min5_pdf, "\n", sep = "")
cat("  ", arg_significant_peak_counts_file, "\n", sep = "")
cat("  ", arg_significant_peak_barplot_pdf, "\n", sep = "")
cat("  ", arg_window_significance_status_file, "\n", sep = "")
cat("  ", arg_peak_significance_status_file, "\n", sep = "")
cat("  ", arg_window_significance_status_summary_file, "\n", sep = "")
cat("  ", arg_peak_significance_status_summary_file, "\n", sep = "")
cat("  ", arg_window_venn_counts_file, "\n", sep = "")
cat("  ", arg_peak_venn_counts_file, "\n", sep = "")
cat("  ", arg_window_ARG_coverage_file, "\n", sep = "")
cat("  ", arg_peak_ARG_coverage_file, "\n", sep = "")
cat("  ", arg_ARG_coverage_summary_file, "\n", sep = "")
cat("  ", model_assignment_any_file, "\n", sep = "")
cat("  ", model_assignment_any_summary_file, "\n", sep = "")
cat("  ", model_assignment_min25_file, "\n", sep = "")
cat("  ", model_assignment_min25_summary_file, "\n", sep = "")
cat("  ", manhattan_any_summary_file, "\n", sep = "")
cat("  ", manhattan_min25_summary_file, "\n", sep = "")
for (suffix in names(arg_significant_peak_counts_minprop_files)) {
  cat("  ", arg_significant_peak_counts_minprop_files[[suffix]], "\n", sep = "")
  cat("  ", arg_significant_peak_barplot_minprop_pdfs[[suffix]], "\n", sep = "")
}
