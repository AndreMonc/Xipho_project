# Xipho_project
Repository for data, scripts, notes, and ARG tutorial associated with manuscript on *Xiphorhynchus spixii*:

*Ancestral recombination graphs highlight selection and demography as drivers of genomic islands of differentiation in an Amazonian bird species*

# Notes

- For quick accessibility, I have added relatively small files (scripts, some input files, step-by-step notes, some output files, etc.) to this repository.
- The Manhattan plots and coalescence trees for all *F*st peaks are found [here](https://github.com/AndreMonc/Xipho_project/tree/main/Peak_summary/Manhattan_plots_and_trees/). Some helpful notes for interpretation:
    * Repeat regions are illustrated in the Manhattan plots as gray bars above the x axis
    * We defined the focal Xingu-Belem *F*st window as the window with the highest *F*st value in the focal peak. This window may not be centered in the Manhattan plot due to the peak occurring close to the beginning or end of a scaffold.
    * All empirical outlier thresholds in the Manhattan plots (shown as horizontal dashed lines) indicate upper-tail thresholds (values above line are outliers), except for Belem pi, Xingu pi, and recombination rate which show lower-tail thresholds (values below line are outliers).
    * Empirical outlier thresholds for the traditional and ARG-based Xingu-Belem *F*st are set at five standard deviations above the mean. Empirical outlier thresholds for all remaining statistics, except recombination rate, iHS, and nSL are set at the 99.9th or 0.1st percentile based on value distributions across control window regions. The empirical outlier threshold for recombination rate is set at the 0.5th percentile because the 0.1st percentile results in a threshold of zero (due to a floor effect caused by the presence of some recombination rate values of zero). The outlier thresholds for iHS and nSL are set at the default thresholds of 2 for absolute normalized values.
    * Model assignments shown at the top of each figure use ARG-based statistics with a peak-level significance defined as having an outlier value (empirical *P* value = 0.001) for at least 25% of the local trees within an *F*st outlier window of the peak. The ARG-based statistics are from the 2000th MCMC iteration of our ARGweaver analysis. 
   
- This set of files is a subset of all the files available on Dryad.
- Files only found on Dryad include but are not limited to: all VCF datasets, phylip input file for IQ-TREE, full PopCluster output, full GADMA output, Pixy output, reference genome and annotation, ARGweaver output smc.gz files, ARG-based *F*st output, RAiSD output, and selscan output.

# Repository files
```
.
├── ARG_based_Fst
│   ├── clean_fst_file.sh
│   ├── individual-species-key-xipho.txt
│   ├── run_branch_fst_iter2000.sbatch
│   └── tre_to_branch_fst_xipho.py
├── ARG_tutorial
│   ├── ARGjobs_xipho.py
│   ├── ARGweaver_jobinfo.txt
│   ├── README.md
│   ├── recomb_map.bed
│   ├── scaffolds_under_110kb.bed
│   └── xipho_ARGweaver_tutorial.vcf
├── ARGweaver_processing
│   ├── ARG_parallel_processing_files
│   │   ├── bed_to_tre_midpoint_xipho_parallel.R
│   │   ├── create_arg_block_region_file_xipho_parallel.sh
│   │   ├── create_bed_files_parallel.sbatch
│   │   ├── create_stat_files_midpoint_parallel.sbatch
│   │   ├── individual-species-key-xipho.txt
│   │   ├── make_region_job_table.sh
│   │   ├── make_region_tree_job_table_midpoint_xipho_parallel.R
│   │   ├── region_job_table.txt
│   │   ├── smc_to_bed_xipho_parallel.R
│   │   ├── tre_to_stats_xipho_parallel.R
│   │   ├── tre_to_stats_xipho_parallel_part2.R
│   │   ├── treeStatFunctions.R
│   │   └── trim_arg_blocks_xipho_parallel.R
│   ├── ARG_run_updated.sbatch
│   ├── ARGweaver_parallel_processing.md
│   ├── ARGweaver_processing.md
│   ├── ARGweaver_running.md
│   ├── Create_repeat_mask_bed_file.md
│   ├── argweaver_parallel_pipeline_README.md
│   ├── bed_files
│   │   ├── ARGweaver_mask.bed
│   │   ├── ARGweaver_windows.bed
│   │   ├── GATK_poor_sites.bed
│   │   ├── ReLERNN_clean_data.bed
│   │   ├── ReLERNN_ft50sites.bed
│   │   ├── bad_scaffolds.bed
│   │   ├── lowQUAL_variant_regions.bed
│   │   ├── non_ReLERNN_scaffs.bed
│   │   ├── repeat_regions.bed
│   │   ├── scaff_ft250SNPs.bed
│   │   ├── sex_chrom_scaffolds.bed
│   │   └── small_scaffolds_110kb.bed
│   └── masked_percentages.txt
├── ARGweaver_trace_plots # contains 573 trace plots across MCMC iterations
│   ├── scaffold1-10_out.png
│   ├── scaffold1-11_out.png
│   ├── scaffold1-12_out.png
│   ├── scaffold1-13_out.png
│   ├── scaffold1-14_out.png
│   ├── scaffold1-15_out.png
│   ├── scaffold1-16_out.png
│   ├── scaffold1-17_out.png
│   ├── scaffold1-18_out.png
│   ├── scaffold1-19_out.png
│   └── etc.
├── Allele_stats
│   ├── VCF_filtering_B_allele_stats.md
│   ├── allele_stats.md
│   ├── bel_popA
│   │   ├── all_sites_map.csv
│   │   ├── allele_stats.py
│   │   ├── allele_stats_Bel_popA.sbatch
│   │   ├── allele_stats_by_site.csv
│   │   ├── allele_stats_by_window_belpopA.csv
│   │   ├── alternate_sites_map.csv
│   │   ├── genome_file.txt
│   │   ├── popKey.txt
│   │   ├── windows.bed
│   │   └── xiph_elegans_ref.fa.fna.fai
│   └── xin_popA
│       ├── all_sites_map.csv
│       ├── allele_stats.py
│       ├── allele_stats_Xin_popA.sbatch
│       ├── allele_stats_by_site.csv
│       ├── allele_stats_by_window_xinpopA.csv
│       ├── alternate_sites_map.csv
│       ├── genome_file.txt
│       ├── popKey.txt
│       ├── windows.bed
│       └── xiph_elegans_ref.fa.fna.fai
├── D_statistics_windows
│   ├── D_stats.md
│   └── xiph.Dstats.csv
├── Dsuite
│   ├── Dsuite.md
│   ├── Dsuite_output
│   │   ├── Dtrios_31098527_stderr.txt
│   │   ├── Dtrios_31098527_stdout.txt
│   │   ├── SETS_BBAA.txt
│   │   ├── SETS_Dmin.txt
│   │   ├── SETS_combine.txt
│   │   └── SETS_combine_stderr.txt
│   ├── Dtrios.sbatch
│   ├── SETS.txt
│   ├── VCF_filtering_new_for_Dstats.md
│   └── vcf.samples.txt
├── Fst_Dxy_Pi
│   ├── Format_pixy_output.R
│   ├── Pixy.md
│   ├── VCF_filtering_for_Fst_Dxy_Pi.md
│   ├── add_chrom_type.sh
│   ├── sex_chrom_scaffolds.bed
│   └── windows.bed
├── GADMA
│   ├── Divergence_bounds.xlsx
│   ├── GADMA.md
│   ├── VCF_filtering_GADMA.md
│   ├── easySFS.md
│   ├── param_file_xipho.yaml
│   └── xiph_pops.txt
├── IQTREE2
│   ├── IQTREE_running.md
│   ├── VCF_filtering_for_IQTREE.md
│   ├── output
│   │   ├── xipho.min4.phy.varsites.phy.bionj
│   │   ├── xipho.min4.phy.varsites.phy.ckp.gz
│   │   ├── xipho.min4.phy.varsites.phy.contree
│   │   ├── xipho.min4.phy.varsites.phy.iqtree
│   │   ├── xipho.min4.phy.varsites.phy.log
│   │   ├── xipho.min4.phy.varsites.phy.mldist
│   │   ├── xipho.min4.phy.varsites.phy.model.gz
│   │   ├── xipho.min4.phy.varsites.phy.splits.nex
│   │   └── xipho.min4.phy.varsites.phy.treefile
│   ├── rename_collapse_tips.R
│   └── tip_name_update.txt
├── Lostruct
│   ├── LoStruct_inversions_HPC_win1000.R
│   ├── LoStruct_inversions_HPCwin500.R
│   ├── VCF_filtering_lostruct.md
│   ├── allele_stats.lostruct.win1000.pdf
│   ├── allele_stats.lostruct.win1000.windows.rds
│   ├── allele_stats.lostruct.win500.pdf
│   ├── allele_stats.lostruct.win500.windows.rds
│   ├── cluster_genotypes.win1000.tsv
│   ├── cluster_genotypes.win500.tsv
│   ├── get_window_stats.R
│   ├── loStruct.md
│   ├── lostruct_window_size_summary.tsv
│   ├── mds_info.win1000.tsv
│   ├── mds_info.win500.tsv
│   ├── sample_info.tsv
│   └── samplelist.txt
├── PCA
│   ├── Genetic_variance_explained_by_eigenvectors.pdf
│   ├── PCA_Eigenvalues.pdf
│   ├── PCA_adegenet_33ind_maxm75.R
│   ├── Raw_PCA.pdf
│   ├── Raw_PCA_with_individual_labels.pdf
│   └── pop_assignments.txt
├── Peak_summary
│   ├── Manhattan_plots_and_trees
│   │   ├── ARG_data_missing
│   │   │   ├── peak_119_scaffold_53_2525000_4525000.pdf
│   │   │   ├── peak_121_scaffold_58_0_2000000.pdf
│   │   │   ├── peak_132_scaffold_73_115000_2115000.pdf
│   │   │   ├── peak_17_scaffold_116_0_1600000.pdf
│   │   │   ├── peak_19_scaffold_12_14330000_16330000.pdf
│   │   │   ├── peak_20_scaffold_121_0_1320000.pdf
│   │   │   ├── peak_21_scaffold_121_0_1320000.pdf
│   │   │   ├── peak_24_scaffold_139_0_1086750.pdf
│   │   │   ├── peak_28_scaffold_14_14780000_16780000.pdf
│   │   │   ├── peak_29_scaffold_147_0_230000.pdf
│   │   │   ├── peak_35_scaffold_161_0_160000.pdf
│   │   │   ├── peak_43_scaffold_192_0_250000.pdf
│   │   │   ├── peak_78_scaffold_271_0_40000.pdf
│   │   │   ├── peak_8_scaffold_109_0_1807944.pdf
│   │   │   ├── peak_93_scaffold_30_0_2000000.pdf
│   │   │   ├── peak_94_scaffold_30_0_2000000.pdf
│   │   │   ├── peak_98_scaffold_332_0_40000.pdf
│   │   │   └── peak_9_scaffold_109_0_1807944.pdf
│   │   ├── Deep_lineage_sorting_model
│   │   │   ├── peak_65_scaffold_23_3510000_5510000.pdf
│   │   │   ├── peak_81_scaffold_28_4485000_6485000.pdf
│   │   │   └── peak_95_scaffold_30_10880000_12880000.pdf
│   │   ├── Overlapping_models
│   │   │   ├── peak_118_scaffold_52_0_2000000.pdf
│   │   │   ├── peak_143_scaffold_8_6845000_8845000.pdf
│   │   │   ├── peak_5_scaffold_1_32860000_34860000.pdf
│   │   │   ├── peak_71_scaffold_25_6890000_8890000.pdf
│   │   │   └── peak_84_scaffold_29_1125000_3125000.pdf
│   │   ├── Selection-bottleneck_model
│   │   │   ├── peak_101_scaffold_38_7500000_9500000.pdf
│   │   │   ├── peak_102_scaffold_39_6620000_8620000.pdf
│   │   │   ├── peak_104_scaffold_40_2075000_4075000.pdf
│   │   │   ├── peak_105_scaffold_40_2855000_4855000.pdf
│   │   │   ├── peak_108_scaffold_41_1035000_3035000.pdf
│   │   │   ├── peak_10_scaffold_11_5035000_7035000.pdf
│   │   │   ├── peak_114_scaffold_5_3655000_5655000.pdf
│   │   │   ├── peak_115_scaffold_5_5475000_7475000.pdf
│   │   │   ├── peak_117_scaffold_50_5980000_7980000.pdf
│   │   │   ├── peak_11_scaffold_11_6060000_8060000.pdf
│   │   │   ├── peak_120_scaffold_53_3695000_5695000.pdf
│   │   │   ├── peak_125_scaffold_6_13905000_15905000.pdf
│   │   │   ├── peak_126_scaffold_64_0_2000000.pdf
│   │   │   ├── peak_12_scaffold_11_6185000_8185000.pdf
│   │   │   ├── peak_133_scaffold_75_0_2000000.pdf
│   │   │   ├── peak_134_scaffold_76_150000_2150000.pdf
│   │   │   ├── peak_136_scaffold_77_955000_2955000.pdf
│   │   │   ├── peak_137_scaffold_77_1135000_3135000.pdf
│   │   │   ├── peak_138_scaffold_77_1425000_3425000.pdf
│   │   │   ├── peak_13_scaffold_11_8095000_10095000.pdf
│   │   │   ├── peak_142_scaffold_8_3280000_5280000.pdf
│   │   │   ├── peak_146_scaffold_8_20450000_22450000.pdf
│   │   │   ├── peak_147_scaffold_80_0_2000000.pdf
│   │   │   ├── peak_148_scaffold_86_625000_2625000.pdf
│   │   │   ├── peak_149_scaffold_86_725000_2725000.pdf
│   │   │   ├── peak_152_scaffold_9_4795000_6795000.pdf
│   │   │   ├── peak_153_scaffold_9_5235000_7235000.pdf
│   │   │   ├── peak_156_scaffold_9_18380000_20380000.pdf
│   │   │   ├── peak_159_scaffold_99_0_2000000.pdf
│   │   │   ├── peak_15_scaffold_11_16340000_18340000.pdf
│   │   │   ├── peak_16_scaffold_115_0_1240000.pdf
│   │   │   ├── peak_18_scaffold_12_0_2000000.pdf
│   │   │   ├── peak_1_scaffold_1_6515000_8515000.pdf
│   │   │   ├── peak_22_scaffold_13_15200000_17200000.pdf
│   │   │   ├── peak_27_scaffold_14_13585000_15585000.pdf
│   │   │   ├── peak_2_scaffold_1_24170000_26170000.pdf
│   │   │   ├── peak_31_scaffold_15_10435000_12435000.pdf
│   │   │   ├── peak_33_scaffold_15_14590000_16590000.pdf
│   │   │   ├── peak_34_scaffold_16_6315000_8315000.pdf
│   │   │   ├── peak_36_scaffold_17_0_2000000.pdf
│   │   │   ├── peak_37_scaffold_17_5435000_7435000.pdf
│   │   │   ├── peak_38_scaffold_17_6095000_8095000.pdf
│   │   │   ├── peak_3_scaffold_1_30800000_32800000.pdf
│   │   │   ├── peak_40_scaffold_19_0_2000000.pdf
│   │   │   ├── peak_42_scaffold_19_12830000_14830000.pdf
│   │   │   ├── peak_44_scaffold_2_835000_2835000.pdf
│   │   │   ├── peak_45_scaffold_2_6465000_8465000.pdf
│   │   │   ├── peak_46_scaffold_2_17165000_19165000.pdf
│   │   │   ├── peak_47_scaffold_2_18620000_20620000.pdf
│   │   │   ├── peak_4_scaffold_1_32365000_34365000.pdf
│   │   │   ├── peak_51_scaffold_20_2015000_4015000.pdf
│   │   │   ├── peak_52_scaffold_20_3825000_5825000.pdf
│   │   │   ├── peak_53_scaffold_20_4065000_6065000.pdf
│   │   │   ├── peak_54_scaffold_20_4820000_6820000.pdf
│   │   │   ├── peak_56_scaffold_20_13400000_15400000.pdf
│   │   │   ├── peak_59_scaffold_22_3285000_5285000.pdf
│   │   │   ├── peak_60_scaffold_22_7845000_9845000.pdf
│   │   │   ├── peak_64_scaffold_23_3175000_5175000.pdf
│   │   │   ├── peak_67_scaffold_23_8235000_10235000.pdf
│   │   │   ├── peak_68_scaffold_24_0_2000000.pdf
│   │   │   ├── peak_70_scaffold_25_3715000_5715000.pdf
│   │   │   ├── peak_72_scaffold_25_7015000_9015000.pdf
│   │   │   ├── peak_79_scaffold_28_4030000_6030000.pdf
│   │   │   ├── peak_7_scaffold_102_185000_2185000.pdf
│   │   │   ├── peak_80_scaffold_28_4140000_6140000.pdf
│   │   │   ├── peak_82_scaffold_28_5100000_7100000.pdf
│   │   │   ├── peak_83_scaffold_29_540000_2540000.pdf
│   │   │   ├── peak_88_scaffold_3_13325000_15325000.pdf
│   │   │   ├── peak_89_scaffold_3_16905000_18905000.pdf
│   │   │   ├── peak_90_scaffold_3_18295000_20295000.pdf
│   │   │   ├── peak_92_scaffold_3_21545000_23545000.pdf
│   │   │   ├── peak_97_scaffold_32_10515000_12515000.pdf
│   │   │   └── peak_99_scaffold_34_10480000_12480000.pdf
│   │   ├── Selection-recombination_model
│   │   │   ├── peak_130_scaffold_7_11085000_13085000.pdf
│   │   │   ├── peak_135_scaffold_76_2260000_4260000.pdf
│   │   │   ├── peak_155_scaffold_9_15055000_17055000.pdf
│   │   │   ├── peak_23_scaffold_13_15200000_17200000.pdf
│   │   │   ├── peak_48_scaffold_2_18775000_20775000.pdf
│   │   │   ├── peak_50_scaffold_20_85000_2085000.pdf
│   │   │   ├── peak_61_scaffold_22_8315000_10315000.pdf
│   │   │   ├── peak_62_scaffold_22_11115000_13115000.pdf
│   │   │   ├── peak_63_scaffold_23_1885000_3885000.pdf
│   │   │   └── peak_77_scaffold_27_4025000_6025000.pdf
│   │   └── Unassigned_to_model
│   │       ├── peak_100_scaffold_38_5315000_7315000.pdf
│   │       ├── peak_103_scaffold_4_25015000_27015000.pdf
│   │       ├── peak_106_scaffold_40_3675000_5675000.pdf
│   │       ├── peak_107_scaffold_40_6875000_8875000.pdf
│   │       ├── peak_109_scaffold_41_2565000_4565000.pdf
│   │       ├── peak_110_scaffold_41_5075000_7075000.pdf
│   │       ├── peak_111_scaffold_42_7983172_9983172.pdf
│   │       ├── peak_112_scaffold_44_4605000_6605000.pdf
│   │       ├── peak_113_scaffold_49_4905000_6905000.pdf
│   │       ├── peak_116_scaffold_50_2845000_4845000.pdf
│   │       ├── peak_122_scaffold_59_2415000_4415000.pdf
│   │       ├── peak_123_scaffold_6_3875000_5875000.pdf
│   │       ├── peak_124_scaffold_6_4755000_6755000.pdf
│   │       ├── peak_127_scaffold_64_775000_2775000.pdf
│   │       ├── peak_128_scaffold_65_0_2000000.pdf
│   │       ├── peak_129_scaffold_7_1525000_3525000.pdf
│   │       ├── peak_131_scaffold_7_12805000_14805000.pdf
│   │       ├── peak_139_scaffold_79_205000_2205000.pdf
│   │       ├── peak_140_scaffold_79_2011640_4011640.pdf
│   │       ├── peak_141_scaffold_79_2011640_4011640.pdf
│   │       ├── peak_144_scaffold_8_9555000_11555000.pdf
│   │       ├── peak_145_scaffold_8_20450000_22450000.pdf
│   │       ├── peak_14_scaffold_11_16225000_18225000.pdf
│   │       ├── peak_150_scaffold_9_2095000_4095000.pdf
│   │       ├── peak_151_scaffold_9_4515000_6515000.pdf
│   │       ├── peak_154_scaffold_9_14505000_16505000.pdf
│   │       ├── peak_157_scaffold_9_18380000_20380000.pdf
│   │       ├── peak_158_scaffold_93_0_2000000.pdf
│   │       ├── peak_25_scaffold_14_5395000_7395000.pdf
│   │       ├── peak_26_scaffold_14_11415000_13415000.pdf
│   │       ├── peak_30_scaffold_15_5725000_7725000.pdf
│   │       ├── peak_32_scaffold_15_13505000_15505000.pdf
│   │       ├── peak_39_scaffold_18_3615000_5615000.pdf
│   │       ├── peak_41_scaffold_19_11495000_13495000.pdf
│   │       ├── peak_49_scaffold_2_24995000_26995000.pdf
│   │       ├── peak_55_scaffold_20_12345000_14345000.pdf
│   │       ├── peak_57_scaffold_21_10975000_12975000.pdf
│   │       ├── peak_58_scaffold_22_1935000_3935000.pdf
│   │       ├── peak_66_scaffold_23_4135000_6135000.pdf
│   │       ├── peak_69_scaffold_24_12010000_14010000.pdf
│   │       ├── peak_6_scaffold_10_3645000_5645000.pdf
│   │       ├── peak_73_scaffold_25_9775000_11775000.pdf
│   │       ├── peak_74_scaffold_25_11790000_13790000.pdf
│   │       ├── peak_75_scaffold_26_5945000_7945000.pdf
│   │       ├── peak_76_scaffold_26_10675000_12675000.pdf
│   │       ├── peak_85_scaffold_29_6455000_8455000.pdf
│   │       ├── peak_86_scaffold_3_4735000_6735000.pdf
│   │       ├── peak_87_scaffold_3_8855000_10855000.pdf
│   │       ├── peak_91_scaffold_3_21015000_23015000.pdf
│   │       └── peak_96_scaffold_32_265000_2265000.pdf
│   └── xipho_summary_script.R
├── PopCluster
│   ├── PopCluster.dat
│   ├── PopCluster.md
│   ├── VCF_filtering_PopCluster.md
│   └── xipho.PcPjt
├── RAiSD
│   ├── Belem_RAiSD.sbatch
│   ├── RAiSD_Belem.md
│   ├── RAiSD_Xingu.md
│   ├── RAiSD_to_bed.py
│   ├── VCF_filtering_RAiSD.md
│   ├── Xingu_RAiSD.sbatch
│   └── windows.bed
├── RCNV
│   ├── genomewide.rCNV.allele_info_WGS.10kb.tsv
│   ├── rCNV.sbatch.sh
│   ├── rCNV_HPC.md
│   ├── rCNV_vcf_split.md
│   └── run_rCNV_one_chrom.R
├── RagTag
│   ├── RagTag_Xiphorhynchus_elegans_to_Chiroxiphia_lanceolata.md
│   ├── ragtag.scaffold.agp
│   ├── sex_chrom_scaffolds.txt
│   └── xiph_elegans_ref_scaffolds.txt
├── ReLERNN
│   ├── ReLERNN_clean_data.bed
│   ├── ReLERNN_notes.md
│   ├── ReLERNN_scaffolds.bed
│   ├── VCF_filtering_new_for_ReLERNN.md
│   ├── raw_output
│   │   └── RelERNN_biallelic_snps_tapajos.PREDICT.txt
│   ├── relernn_30707823_stderr.txt
│   └── relernn_30707823_stdout.txt
├── Selscan
│   ├── VCF_filtering_for_selscan.md
│   ├── selscan.md
│   ├── selscan_Bel.sbatch
│   └── selscan_Xin.sbatch
├── SnpArcher
│   ├── Downloading_X_spixii.md
│   ├── SNParcher_notes_OU_HPC.md
│   ├── Xipho_sample_sheet.csv
│   ├── bam2vcf_gatk_intervals.smk
│   ├── config
│   │   └── config.yaml
│   ├── output
│   │   ├── FINAL_XIPHO.king
│   │   ├── FINAL_XIPHO.king.id
│   │   ├── FINAL_XIPHO_callable_sites.bed
│   │   └── FINAL_XIPHO_qc.html
│   ├── rename_fasta.py
│   ├── scaff_len.py
│   ├── scaffold_list.txt
│   ├── scaffold_rename.csv
│   └── slurm
│       └── config.yaml
├── full_file_tree.txt
└── github_full_file_tree.txt

38 directories, 904 files
```