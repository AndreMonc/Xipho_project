# xipho_ARGs_project
Repository for data, scripts, notes, and ARG tutorial associated with manuscript on *Xiphorhynchus spixii*

# Notes

- For quick accessibility, I have added relatively small files (scripts, some input files, step-by-step notes, some output files, etc.) to this repository.
- The Manhattan plot panels for all *F*st peaks are found [here](https://github.com/AndreMonc/Xipho_project/tree/main/peak_summary/output_files/Manhattan_panels).
- This set of files is a subset of all the files available on Dryad.
- Files only found on Dryad include: all VCF files, phylip input file for IQ-TREE, full PopCluster output, full GADMA output, reference genomes, ARGweaver input .bed files, ARGweaver output SMC.gz and log files, and RAiSD output.

# Repository files
```
.
├── ARG_tutorial
│   ├── ARGjobs_xipho.py
│   ├── ARGweaver_jobinfo.txt
│   ├── README.md
│   ├── recomb_map.bed
│   ├── scaffolds_under_110kb.bed
│   └── xipho_ARGweaver_tutorial.vcf
├── ARGweaver
│   ├── ARGweaver_notes_downstream_processing_notes.md
│   ├── ARGweaver_notes_running.md
│   ├── ARGweaver_scripts
│   │   └── tre_to_stats.R
│   ├── Create_500bp_Bam_depth_masks.md
│   ├── Create_repeat_mask_bedfile.md
│   └── ind_mask_file.txt
├── GADMA
│   ├── Divergence_bounds.xlsx
│   ├── GADMA_VCF_filtering_notes.md
│   ├── GADMA_runs
│   │   ├── GADMA.log
│   │   ├── best_aic_model.png
│   │   ├── best_aic_model_dadi_code.py
│   │   ├── best_aic_model_demes_code.py.yml
│   │   ├── best_aic_model_moments_code.py
│   │   ├── best_logLL_model.png
│   │   ├── best_logLL_model_dadi_code.py
│   │   ├── best_logLL_model_demes_code.py.yml
│   │   ├── best_logLL_model_moments_code.py
│   │   ├── extra_params_file
│   │   └── params_file
│   ├── Running_GADMA.md
│   ├── Running_easySFS.md
│   ├── param_file_xipho_run5.yaml
│   └── xiph_pops.txt
├── IQTREE
│   ├── VCF_filtering_for_IQTREE_and_running_IQTREE.md
│   ├── rename_collapse_tips.R
│   ├── xipho.min4.phy.varsites.phy.treefile
│   └── xipho_tree_polytomies.tre
├── PCA
│   ├── PCA_adegenet_33ind_maxm75.R
│   └── pop_assignments.txt
├── PopCluster
│   ├── PopCluster_xipho.qsub
│   ├── VCF_filtering_cluster_analysis_notes.md
│   ├── output_files
│   │   └── xipho_popc.K
│   └── xipho.PcPjt
├── RAiSD
│   ├── RAiSD_Belem_notes.md
│   ├── RAiSD_Xingu_notes.md
│   ├── RAiSD_results_Belem
│   │   ├── RAiSD_to_bed.py
│   │   ├── xipho_10kbwindows_112136.bed
│   │   └── xipho_master_stats_with_inverseRTHs_with_mu.txt
│   └── RAiSD_results_Xingu
│       ├── RAiSD_to_bed.py
│       ├── xipho_10kbwindows_112136.bed
│       └── xipho_master_stats_with_inverseRTHs_with_Belmu_and_Xinmu.txt
├── README.md
├── RagTag
│   └── RagTag_Xiphorhynchus_elegans_pseudochromosome.md
├── ReLERNN
│   ├── VCF_filtering_new_for_ReLERNN_notes.md
│   ├── raw_output
│   │   ├── ReLERNN_Tapajos_Z_hiqual_final.recode.PREDICT.txt
│   │   └── ReLERNN_Tapajos_autosomes_hiqual_final.recode.PREDICT.txt
│   └── xipho_ReLERNN_notes.md
├── allele_stats_calcs
│   ├── allele_stats.md
│   ├── bel_popA
│   │   ├── all_sites_map.csv
│   │   ├── allele_stats.py
│   │   ├── allele_stats_Bel_popA.qsub
│   │   ├── allele_stats_Bel_popA.qsub.e414298
│   │   ├── allele_stats_Bel_popA.qsub.o414298
│   │   ├── allele_stats_by_site.csv
│   │   ├── allele_stats_by_window.csv
│   │   ├── alternate_sites_map.csv
│   │   ├── genome_file.txt
│   │   ├── popKey.txt
│   │   └── windows.bed
│   ├── bel_popB
│   │   ├── all_sites_map.csv
│   │   ├── allele_stats.py
│   │   ├── allele_stats_Bel_popB.qsub
│   │   ├── allele_stats_Bel_popB.qsub.e414338
│   │   ├── allele_stats_Bel_popB.qsub.o414338
│   │   ├── allele_stats_by_site.csv
│   │   ├── allele_stats_by_window.csv
│   │   ├── alternate_sites_map.csv
│   │   ├── genome_file.txt
│   │   ├── popKey.txt
│   │   └── windows.bed
│   └── popKey.xlsx
├── dataframes
│   ├── control_dataframe_autosomes_and_Z.csv
│   └── master_dataframe.txt
├── peak_summary
│   ├── input_files
│   │   ├── ARG_thresholds_p0.0001.txt
│   │   ├── ARG_thresholds_p0.001.txt
│   │   ├── ARG_thresholds_p0.01.txt
│   │   └── master_dataframe.txt
│   ├── output_files
│   │   ├── ARG_significant_Fst_outliers_barplot.pdf
│   │   ├── ARG_significant_Fst_peaks_barplot.pdf
│   │   ├── Fst_Fst_outlier_peaks_50kb.txt
│   │   ├── Fst_peak_heatmap_percentiles.pdf
│   │   ├── Fst_peak_heatmap_percentiles_table.txt
│   │   ├── Fst_peak_model_assignments.txt
│   │   ├── Fst_xingu_belem_Z_outliers_5SD.txt
│   │   ├── Fst_xingu_belem_autosome_outliers_5SD.txt
│   │   ├── Manhattan_panels
│   │   │   ├── Deep_lineage_sorting
│   │   │   │   ├── Chromosome_14_RagTag__13110000_13120000.pdf
│   │   │   │   └── Chromosome_5_RagTag__15900000_15920000.pdf
│   │   │   ├── Genomic_architecture
│   │   │   │   ├── Chromosome_15_RagTag__19100000_19110000.pdf
│   │   │   │   ├── Chromosome_19_RagTag__870000_880000.pdf
│   │   │   │   ├── Chromosome_7_RagTag__2670000_2690000.pdf
│   │   │   │   └── Chromosome_Z_RagTag__38400000_38460000.pdf
│   │   │   ├── Introgression
│   │   │   │   ├── Chromosome_1_RagTag__65190000_65200000.pdf
│   │   │   │   ├── Chromosome_3_RagTag__23690000_23730000.pdf
│   │   │   │   ├── Chromosome_4_RagTag__43440000_43450000.pdf
│   │   │   │   ├── Chromosome_5_RagTag__28880000_28890000.pdf
│   │   │   │   ├── Chromosome_5_RagTag__60080000_60090000.pdf
│   │   │   │   └── Chromosome_Z_RagTag__47680000_47690000.pdf
│   │   │   ├── Overlapping_models
│   │   │   │   ├── Chromosome_10_RagTag__18250000_18270000.pdf
│   │   │   │   ├── Chromosome_19_RagTag__720000_8e+05.pdf
│   │   │   │   ├── Chromosome_1_RagTag__119880000_119910000.pdf
│   │   │   │   ├── Chromosome_1_RagTag__71330000_71350000.pdf
│   │   │   │   ├── Chromosome_1_RagTag__86830000_86880000.pdf
│   │   │   │   ├── Chromosome_25_RagTag__6650000_6840000.pdf
│   │   │   │   ├── Chromosome_3_RagTag__36680000_36730000.pdf
│   │   │   │   ├── Chromosome_3_RagTag__47030000_47040000.pdf
│   │   │   │   ├── Chromosome_4_RagTag__40180000_40200000.pdf
│   │   │   │   ├── Chromosome_4_RagTag__71450000_71470000.pdf
│   │   │   │   └── Chromosome_7_RagTag__22120000_22150000.pdf
│   │   │   ├── Selection-bottleneck_Bottleneck
│   │   │   │   ├── Chromosome_11_RagTag__22430000_22440000.pdf
│   │   │   │   ├── Chromosome_12_RagTag__9300000_9310000.pdf
│   │   │   │   ├── Chromosome_13_RagTag__15150000_15160000.pdf
│   │   │   │   ├── Chromosome_14_RagTag__3570000_3580000.pdf
│   │   │   │   ├── Chromosome_16_RagTag__13030000_13040000.pdf
│   │   │   │   ├── Chromosome_1_RagTag__128630000_128640000.pdf
│   │   │   │   ├── Chromosome_1_RagTag__28570000_28580000.pdf
│   │   │   │   ├── Chromosome_1_RagTag__37440000_37450000.pdf
│   │   │   │   ├── Chromosome_1_RagTag__66770000_66780000.pdf
│   │   │   │   ├── Chromosome_1_RagTag__76060000_76070000.pdf
│   │   │   │   ├── Chromosome_1_RagTag__77560000_77570000.pdf
│   │   │   │   ├── Chromosome_3_RagTag__116360000_116370000.pdf
│   │   │   │   ├── Chromosome_3_RagTag__52380000_52390000.pdf
│   │   │   │   ├── Chromosome_4_RagTag__31490000_31500000.pdf
│   │   │   │   ├── Chromosome_4_RagTag__51430000_51440000.pdf
│   │   │   │   ├── Chromosome_4_RagTag__70390000_70400000.pdf
│   │   │   │   ├── Chromosome_5_RagTag__15320000_15330000.pdf
│   │   │   │   ├── Chromosome_5_RagTag__33690000_33700000.pdf
│   │   │   │   ├── Chromosome_6_RagTag__34460000_34470000.pdf
│   │   │   │   ├── Chromosome_8_RagTag__33920000_33930000.pdf
│   │   │   │   ├── Chromosome_9_RagTag__3200000_3210000.pdf
│   │   │   │   ├── Chromosome_Z_RagTag__14200000_14210000.pdf
│   │   │   │   └── scaffold_113__780000_790000.pdf
│   │   │   ├── Selection-bottleneck_Other
│   │   │   │   ├── Chromosome_10_RagTag__17880000_17890000.pdf
│   │   │   │   ├── Chromosome_12_RagTag__6e+06_6020000.pdf
│   │   │   │   ├── Chromosome_12_RagTag__9170000_9190000.pdf
│   │   │   │   ├── Chromosome_14_RagTag__1540000_1550000.pdf
│   │   │   │   ├── Chromosome_1_RagTag__42600000_42640000.pdf
│   │   │   │   ├── Chromosome_1_RagTag__53850000_53880000.pdf
│   │   │   │   ├── Chromosome_1_RagTag__99790000_99810000.pdf
│   │   │   │   ├── Chromosome_21_RagTag__360000_370000.pdf
│   │   │   │   ├── Chromosome_26_RagTag__2010000_2020000.pdf
│   │   │   │   ├── Chromosome_29_RagTag__2480000_2530000.pdf
│   │   │   │   ├── Chromosome_2_RagTag__29800000_29810000.pdf
│   │   │   │   ├── Chromosome_2_RagTag__97340000_97350000.pdf
│   │   │   │   ├── Chromosome_3_RagTag__12650000_12660000.pdf
│   │   │   │   ├── Chromosome_3_RagTag__39500000_39520000.pdf
│   │   │   │   ├── Chromosome_4_RagTag__64820000_64840000.pdf
│   │   │   │   ├── Chromosome_4_RagTag__70520000_70530000.pdf
│   │   │   │   ├── Chromosome_5_RagTag__22810000_22820000.pdf
│   │   │   │   ├── Chromosome_5_RagTag__32900000_32920000.pdf
│   │   │   │   ├── Chromosome_5_RagTag__49820000_49830000.pdf
│   │   │   │   ├── Chromosome_6_RagTag__18640000_18660000.pdf
│   │   │   │   ├── Chromosome_6_RagTag__34330000_34340000.pdf
│   │   │   │   ├── Chromosome_6_RagTag__49360000_49380000.pdf
│   │   │   │   └── Chromosome_6_RagTag__51160000_51170000.pdf
│   │   │   ├── Selection-bottleneck_Selection
│   │   │   │   ├── Chromosome_10_RagTag__14660000_14720000.pdf
│   │   │   │   ├── Chromosome_10_RagTag__24600000_24620000.pdf
│   │   │   │   ├── Chromosome_13_RagTag__12890000_1.3e+07.pdf
│   │   │   │   ├── Chromosome_14_RagTag__19930000_19950000.pdf
│   │   │   │   ├── Chromosome_14_RagTag__20070000_20260000.pdf
│   │   │   │   ├── Chromosome_14_RagTag__2640000_2700000.pdf
│   │   │   │   ├── Chromosome_16_RagTag__17490000_17510000.pdf
│   │   │   │   ├── Chromosome_17_RagTag__13860000_13890000.pdf
│   │   │   │   ├── Chromosome_19_RagTag__12720000_12840000.pdf
│   │   │   │   ├── Chromosome_1_RagTag__119380000_119450000.pdf
│   │   │   │   ├── Chromosome_1_RagTag__122660000_122700000.pdf
│   │   │   │   ├── Chromosome_1_RagTag__137790000_137930000.pdf
│   │   │   │   ├── Chromosome_2_RagTag__39620000_39640000.pdf
│   │   │   │   ├── Chromosome_4_RagTag__70750000_70940000.pdf
│   │   │   │   ├── Chromosome_5_RagTag__42100000_42250000.pdf
│   │   │   │   ├── Chromosome_6_RagTag__34570000_34590000.pdf
│   │   │   │   ├── Chromosome_6_RagTag__37560000_37580000.pdf
│   │   │   │   ├── Chromosome_6_RagTag__40130000_40190000.pdf
│   │   │   │   ├── Chromosome_7_RagTag__2780000_2870000.pdf
│   │   │   │   └── Chromosome_9_RagTag__20980000_21010000.pdf
│   │   │   └── Unassigned_to_model
│   │   │       ├── Chromosome_10_RagTag__24310000_24360000.pdf
│   │   │       ├── Chromosome_10_RagTag__24780000_24790000.pdf
│   │   │       ├── Chromosome_10_RagTag__2e+05_210000.pdf
│   │   │       ├── Chromosome_10_RagTag__440000_450000.pdf
│   │   │       ├── Chromosome_10_RagTag__8080000_8090000.pdf
│   │   │       ├── Chromosome_11_RagTag__14040000_14050000.pdf
│   │   │       ├── Chromosome_12_RagTag__12060000_12070000.pdf
│   │   │       ├── Chromosome_12_RagTag__15030000_15050000.pdf
│   │   │       ├── Chromosome_14_RagTag__3650000_3660000.pdf
│   │   │       ├── Chromosome_15_RagTag__11120000_11130000.pdf
│   │   │       ├── Chromosome_15_RagTag__14880000_14890000.pdf
│   │   │       ├── Chromosome_16_RagTag__1780000_1790000.pdf
│   │   │       ├── Chromosome_17_RagTag__12530000_12540000.pdf
│   │   │       ├── Chromosome_17_RagTag__130000_140000.pdf
│   │   │       ├── Chromosome_18_RagTag__13520000_13530000.pdf
│   │   │       ├── Chromosome_19_RagTag__5550000_5560000.pdf
│   │   │       ├── Chromosome_1_RagTag__104820000_104830000.pdf
│   │   │       ├── Chromosome_1_RagTag__113500000_113510000.pdf
│   │   │       ├── Chromosome_1_RagTag__114860000_114870000.pdf
│   │   │       ├── Chromosome_1_RagTag__116150000_116160000.pdf
│   │   │       ├── Chromosome_1_RagTag__128290000_128300000.pdf
│   │   │       ├── Chromosome_1_RagTag__139570000_139580000.pdf
│   │   │       ├── Chromosome_1_RagTag__19870000_19880000.pdf
│   │   │       ├── Chromosome_1_RagTag__30730000_30740000.pdf
│   │   │       ├── Chromosome_1_RagTag__32690000_32700000.pdf
│   │   │       ├── Chromosome_1_RagTag__3770000_3780000.pdf
│   │   │       ├── Chromosome_1_RagTag__38620000_38630000.pdf
│   │   │       ├── Chromosome_1_RagTag__54060000_54070000.pdf
│   │   │       ├── Chromosome_1_RagTag__69480000_69490000.pdf
│   │   │       ├── Chromosome_1_RagTag__69680000_69750000.pdf
│   │   │       ├── Chromosome_1_RagTag__71960000_71970000.pdf
│   │   │       ├── Chromosome_1_RagTag__7700000_7710000.pdf
│   │   │       ├── Chromosome_21_RagTag__480000_5e+05.pdf
│   │   │       ├── Chromosome_23_RagTag__7160000_7170000.pdf
│   │   │       ├── Chromosome_27_RagTag__6210000_6220000.pdf
│   │   │       ├── Chromosome_29_RagTag__1570000_1580000.pdf
│   │   │       ├── Chromosome_29_RagTag__2620000_2630000.pdf
│   │   │       ├── Chromosome_2_RagTag__108440000_108450000.pdf
│   │   │       ├── Chromosome_2_RagTag__115640000_115650000.pdf
│   │   │       ├── Chromosome_2_RagTag__14230000_14240000.pdf
│   │   │       ├── Chromosome_2_RagTag__15120000_15130000.pdf
│   │   │       ├── Chromosome_2_RagTag__15410000_15420000.pdf
│   │   │       ├── Chromosome_2_RagTag__25960000_25980000.pdf
│   │   │       ├── Chromosome_2_RagTag__26280000_26290000.pdf
│   │   │       ├── Chromosome_2_RagTag__27460000_27470000.pdf
│   │   │       ├── Chromosome_2_RagTag__30350000_30360000.pdf
│   │   │       ├── Chromosome_2_RagTag__40060000_40080000.pdf
│   │   │       ├── Chromosome_2_RagTag__40340000_40350000.pdf
│   │   │       ├── Chromosome_2_RagTag__42760000_42770000.pdf
│   │   │       ├── Chromosome_2_RagTag__46790000_46800000.pdf
│   │   │       ├── Chromosome_2_RagTag__55110000_55120000.pdf
│   │   │       ├── Chromosome_2_RagTag__6920000_6930000.pdf
│   │   │       ├── Chromosome_2_RagTag__70580000_70590000.pdf
│   │   │       ├── Chromosome_2_RagTag__77980000_77990000.pdf
│   │   │       ├── Chromosome_2_RagTag__79800000_79810000.pdf
│   │   │       ├── Chromosome_2_RagTag__80650000_80660000.pdf
│   │   │       ├── Chromosome_2_RagTag__9270000_9280000.pdf
│   │   │       ├── Chromosome_30_RagTag__820000_830000.pdf
│   │   │       ├── Chromosome_3_RagTag__103390000_103400000.pdf
│   │   │       ├── Chromosome_3_RagTag__12580000_12590000.pdf
│   │   │       ├── Chromosome_3_RagTag__17120000_17130000.pdf
│   │   │       ├── Chromosome_3_RagTag__19480000_19490000.pdf
│   │   │       ├── Chromosome_3_RagTag__21670000_21680000.pdf
│   │   │       ├── Chromosome_3_RagTag__23580000_23590000.pdf
│   │   │       ├── Chromosome_3_RagTag__24720000_24740000.pdf
│   │   │       ├── Chromosome_3_RagTag__2770000_2780000.pdf
│   │   │       ├── Chromosome_3_RagTag__29880000_29890000.pdf
│   │   │       ├── Chromosome_3_RagTag__30770000_30780000.pdf
│   │   │       ├── Chromosome_3_RagTag__34780000_34790000.pdf
│   │   │       ├── Chromosome_3_RagTag__38510000_38520000.pdf
│   │   │       ├── Chromosome_3_RagTag__38750000_38760000.pdf
│   │   │       ├── Chromosome_3_RagTag__39390000_39400000.pdf
│   │   │       ├── Chromosome_3_RagTag__48710000_48720000.pdf
│   │   │       ├── Chromosome_3_RagTag__52620000_52630000.pdf
│   │   │       ├── Chromosome_3_RagTag__52760000_52770000.pdf
│   │   │       ├── Chromosome_3_RagTag__55200000_55210000.pdf
│   │   │       ├── Chromosome_3_RagTag__58970000_58980000.pdf
│   │   │       ├── Chromosome_3_RagTag__66260000_66270000.pdf
│   │   │       ├── Chromosome_3_RagTag__66650000_66660000.pdf
│   │   │       ├── Chromosome_3_RagTag__72380000_72390000.pdf
│   │   │       ├── Chromosome_3_RagTag__76400000_76410000.pdf
│   │   │       ├── Chromosome_3_RagTag__83190000_83210000.pdf
│   │   │       ├── Chromosome_3_RagTag__87770000_87780000.pdf
│   │   │       ├── Chromosome_3_RagTag__98020000_98030000.pdf
│   │   │       ├── Chromosome_4_RagTag__3.2e+07_32010000.pdf
│   │   │       ├── Chromosome_4_RagTag__33570000_33590000.pdf
│   │   │       ├── Chromosome_4_RagTag__37030000_37040000.pdf
│   │   │       ├── Chromosome_4_RagTag__38100000_38110000.pdf
│   │   │       ├── Chromosome_4_RagTag__50910000_50980000.pdf
│   │   │       ├── Chromosome_4_RagTag__55310000_55320000.pdf
│   │   │       ├── Chromosome_4_RagTag__57820000_57860000.pdf
│   │   │       ├── Chromosome_4_RagTag__64480000_64490000.pdf
│   │   │       ├── Chromosome_5_RagTag__15510000_15520000.pdf
│   │   │       ├── Chromosome_5_RagTag__31390000_31400000.pdf
│   │   │       ├── Chromosome_5_RagTag__32090000_32100000.pdf
│   │   │       ├── Chromosome_5_RagTag__37630000_37640000.pdf
│   │   │       ├── Chromosome_5_RagTag__41670000_41680000.pdf
│   │   │       ├── Chromosome_5_RagTag__43450000_43460000.pdf
│   │   │       ├── Chromosome_5_RagTag__50420000_50430000.pdf
│   │   │       ├── Chromosome_5_RagTag__52960000_52970000.pdf
│   │   │       ├── Chromosome_5_RagTag__61970000_61980000.pdf
│   │   │       ├── Chromosome_5_RagTag__65040000_65050000.pdf
│   │   │       ├── Chromosome_6_RagTag__10960000_10970000.pdf
│   │   │       ├── Chromosome_6_RagTag__23120000_23130000.pdf
│   │   │       ├── Chromosome_6_RagTag__24890000_24900000.pdf
│   │   │       ├── Chromosome_6_RagTag__27240000_27250000.pdf
│   │   │       ├── Chromosome_6_RagTag__42770000_42780000.pdf
│   │   │       ├── Chromosome_6_RagTag__56430000_56440000.pdf
│   │   │       ├── Chromosome_6_RagTag__57960000_57970000.pdf
│   │   │       ├── Chromosome_7_RagTag__12980000_12990000.pdf
│   │   │       ├── Chromosome_7_RagTag__21990000_2.2e+07.pdf
│   │   │       ├── Chromosome_7_RagTag__24320000_24330000.pdf
│   │   │       ├── Chromosome_7_RagTag__26490000_26500000.pdf
│   │   │       ├── Chromosome_7_RagTag__32510000_32520000.pdf
│   │   │       ├── Chromosome_8_RagTag__20980000_20990000.pdf
│   │   │       ├── Chromosome_8_RagTag__26010000_26020000.pdf
│   │   │       ├── Chromosome_8_RagTag__32820000_32830000.pdf
│   │   │       ├── Chromosome_8_RagTag__34120000_34130000.pdf
│   │   │       ├── Chromosome_9_RagTag__19540000_19550000.pdf
│   │   │       ├── Chromosome_9_RagTag__21150000_21160000.pdf
│   │   │       ├── Chromosome_9_RagTag__27370000_27380000.pdf
│   │   │       ├── Chromosome_9_RagTag__8840000_8850000.pdf
│   │   │       ├── Chromosome_Z_RagTag__18500000_18510000.pdf
│   │   │       ├── Chromosome_Z_RagTag__18750000_18760000.pdf
│   │   │       ├── Chromosome_Z_RagTag__52490000_52500000.pdf
│   │   │       └── scaffold_113__860000_870000.pdf
│   │   └── master_dataframe_filtered.txt
│   └── xipho_summary_script.R
├── snpArcher
│   ├── FINAL_pseudo_XIPHO_qc.html
│   ├── SNP_archer_xiphorhynchus_no_W.md
│   ├── Xiphorhynchus_sample_sheet.csv
│   └── config.yaml
└── traditional_summ_stats
    ├── Fst_final_output.txt
    └── VCF_filtering_for_Fst_and_runningFst_notes.md
```