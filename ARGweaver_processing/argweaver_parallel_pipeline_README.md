# Parallel ARGweaver Processing Pipeline Across MCMC Iterations

This pipeline processes ARGweaver `.smc.gz` files from multiple MCMC iterations across many genomic regions.

The workflow is designed for ARGweaver files from iterations:

```text
1510, 1520, 1530, ..., 2000
```

That is:

```text
50 MCMC iterations per region
```

The input files are stored in scaffold-region folders:

```text
/scratch/amonc/xipho_revision/arg_processing/verified_smcs/
└── scaffold100-553/
    ├── scaffold100-553_out.1510.smc.gz
    ├── scaffold100-553_out.1520.smc.gz
    ├── ...
    └── scaffold100-553_out.2000.smc.gz
```

The working directory for this pipeline is:

```bash
/scratch/amonc/xipho_revision/arg_processing/args_parallel
```

---

## Directory setup

From:

```bash
cd /scratch/amonc/xipho_revision/arg_processing
```

create the clean working directory:

```bash
mkdir -p args_parallel/{scripts,logs,infoTables,argBedFiles,argTreeFiles_midpoint,argStats_midpoint,tmp}
cd args_parallel
```

The final directory structure should look like:

```text
args_parallel/
├── argBedFiles/
├── argStats_midpoint/
├── argTreeFiles_midpoint/
├── infoTables/
├── logs/
├── scripts/
└── tmp/
```

---

## Required metadata files

The population/species key should be located here:

```text
infoTables/individual-species-key-xipho.txt
```

The tree-stat helper functions should be accessible as either:

```text
treeStatFunctions.R
```

or:

```text
scripts/treeStatFunctions.R
```

---

# Step 1. Create and trim ARG block coordinate files

## 1A. Create untrimmed ARG block coordinate table

Script:

```text
scripts/create_arg_block_region_file_xipho_parallel.sh
```

Purpose:

This script scans the `1510` SMC file for each region and extracts the `REGION` line. We only need one iteration per region because genomic coordinates are constant across MCMC iterations.

Input:

```text
/scratch/amonc/xipho_revision/arg_processing/verified_smcs/<region>/<region>_out.1510.smc.gz
```

Output:

```text
infoTables/ARGblock-coordinates.txt
```

Run:

```bash
cd /scratch/amonc/xipho_revision/arg_processing/args_parallel

chmod +x scripts/create_arg_block_region_file_xipho_parallel.sh
./scripts/create_arg_block_region_file_xipho_parallel.sh
```

Check:

```bash
head infoTables/ARGblock-coordinates.txt
wc -l infoTables/ARGblock-coordinates.txt
```

Expected:

```text
576 lines
```

This is:

```text
575 input regions + 1 header
```

---

## 1B. Trim ARG blocks

Script:

```text
scripts/trim_arg_blocks_xipho_parallel.R
```

Purpose:

This script trims ARG blocks to remove edge effects and drops blocks shorter than the minimum analyzable length.

Main parameters:

```r
min_block_len  <- 101000
segment_len    <- 1000
trim_len       <- 50000
```

Input:

```text
infoTables/ARGblock-coordinates.txt
```

Outputs:

```text
infoTables/ARGblock-coordinates-trimmed.txt
infoTables/ARGblock-coordinates-info.txt
```

Run:

```bash
module load R/4.2.1-foss-2022a

cd /scratch/amonc/xipho_revision/arg_processing/args_parallel
Rscript scripts/trim_arg_blocks_xipho_parallel.R
```

Check:

```bash
wc -l infoTables/ARGblock-coordinates-trimmed.txt
cat infoTables/ARGblock-coordinates-info.txt
```

Observed result:

```text
575 input regions
2 too short after filtering
573 retained trimmed regions
```

Expected trimmed table line count:

```text
573 lines
```

Note: `ARGblock-coordinates-trimmed.txt` has no header.

---

# Step 2. Generate ARG BED files

This step converts SMC files to indexed BED files using `smc2bed`.

This step uses the **untrimmed** coordinate table:

```text
infoTables/ARGblock-coordinates.txt
```

because BED files should represent the full ARGweaver SMC output. Trimming is applied later when extracting midpoint trees.

---

## 2A. Create per-region job table

Script:

```text
scripts/make_region_job_table.sh
```

Purpose:

Creates one job-table row per region.

Input:

```text
infoTables/ARGblock-coordinates.txt
```

Output:

```text
infoTables/region_job_table.txt
```

Run:

```bash
cd /scratch/amonc/xipho_revision/arg_processing/args_parallel

chmod +x scripts/make_region_job_table.sh
./scripts/make_region_job_table.sh
```

Check:

```bash
head infoTables/region_job_table.txt
wc -l infoTables/region_job_table.txt
```

Expected:

```text
576 lines
```

This is:

```text
575 regions + 1 header
```

---

## 2B. Convert SMC files to BED files

Scripts:

```text
scripts/create_bed_files_parallel.sbatch
scripts/smc_to_bed_xipho_parallel.R
```

Strategy:

One Slurm array task is submitted per region:

```text
575 array tasks
```

Each region task loops over all 50 MCMC iterations:

```text
1510, 1520, ..., 2000
```

Input SMC files:

```text
/scratch/amonc/xipho_revision/arg_processing/verified_smcs/<region>/<region>_out.<iter>.smc.gz
```

Input ARGweaver log files:

```text
/scratch/amonc/xipho_revision/arg_processing/verified_logs/<region>_out.log
```

Output BED files:

```text
argBedFiles/<region>/<region>_out.<iter>.bed.gz
argBedFiles/<region>/<region>_out.<iter>.bed.gz.tbi
```

Example:

```text
argBedFiles/scaffold100-553/scaffold100-553_out.1510.bed.gz
argBedFiles/scaffold100-553/scaffold100-553_out.1510.bed.gz.tbi
```

Submit:

```bash
sbatch scripts/create_bed_files_parallel.sbatch
```

Check progress:

```bash
squeue -u $USER
```

Check final counts:

```bash
find argBedFiles -name "*.bed.gz" | wc -l
find argBedFiles -name "*.bed.gz.tbi" | wc -l
find argBedFiles -name "*.bed.gz" -size 0 | wc -l
```

Expected:

```text
28750 .bed.gz files
28750 .bed.gz.tbi files
0 empty .bed.gz files
```

This is:

```text
575 regions × 50 iterations = 28,750 BED files
```

Important runtime note:

If `smc2bed` gives a `GLIBCXX` or `CXXABI` error, the job is using an old system C++ runtime. Load a newer GCC module or set `LD_LIBRARY_PATH` before running `smc2bed`, for example in the sbatch script:

```bash
module load GCC/11.3.0
```

or, if using conda:

```bash
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
```

---

# Step 3. Generate midpoint TRE files

This step converts indexed BED files into midpoint TRE files.

Each output TRE file contains one row per local-tree span:

```text
midpoint_position    Newick_tree
```

The midpoint is calculated after clipping each local-tree span to the trimmed analyzed region.

This step uses the **trimmed** coordinate table:

```text
infoTables/ARGblock-coordinates-trimmed.txt
```

---

## 3A. Create midpoint tree job table

Script:

```text
scripts/make_region_tree_job_table_midpoint_xipho_parallel.R
```

Input:

```text
infoTables/ARGblock-coordinates-trimmed.txt
```

Output:

```text
infoTables/region_tree_job_table_midpoint.txt
```

Run:

```bash
module load R/4.2.1-foss-2022a

cd /scratch/amonc/xipho_revision/arg_processing/args_parallel
Rscript scripts/make_region_tree_job_table_midpoint_xipho_parallel.R
```

Check:

```bash
head infoTables/region_tree_job_table_midpoint.txt
wc -l infoTables/region_tree_job_table_midpoint.txt
```

Expected:

```text
574 lines
```

This is:

```text
573 trimmed regions + 1 header
```

---

## 3B. Convert BED files to midpoint TRE files

Scripts:

```text
scripts/create_tree_files_midpoint_parallel.sbatch
scripts/bed_to_tre_midpoint_xipho_parallel.R
```

Strategy:

One Slurm array task is submitted per trimmed region:

```text
573 array tasks
```

Each region task loops over all 50 MCMC iterations.

Input BED files:

```text
argBedFiles/<region>/<region>_out.<iter>.bed.gz
```

Output TRE files:

```text
argTreeFiles_midpoint/<region>/<region>.<start>-<end>.<iter>.tre.gz
```

Example:

```text
argTreeFiles_midpoint/scaffold100-553/scaffold100-553.50001-1950000.1510.tre.gz
```

Submit:

```bash
sbatch scripts/create_tree_files_midpoint_parallel.sbatch
```

Check final counts:

```bash
find argTreeFiles_midpoint -name "*.tre.gz" | wc -l
find argTreeFiles_midpoint -name "*.tre.gz" -size 0 | wc -l
```

Expected:

```text
28650 .tre.gz files
0 empty .tre.gz files
```

This is:

```text
573 trimmed regions × 50 iterations = 28,650 TRE files
```

Spot-check:

```bash
find argTreeFiles_midpoint -name "*.tre.gz" | head -1
zcat $(find argTreeFiles_midpoint -name "*.tre.gz" | head -1) | head
```

Expected format:

```text
position    Newick_tree
```

---

# Step 4. Generate STAT files

This step converts midpoint TRE files into per-tree statistic files.

Each output STAT file contains one row per local tree / midpoint position.

---

## 4A. Inputs for stats calculation

Required files:

```text
argTreeFiles_midpoint/<region>/<region>.<start>-<end>.<iter>.tre.gz
infoTables/individual-species-key-xipho.txt
treeStatFunctions.R
```

The individual-species key maps sample IDs to populations, for example:

```text
Tap
Xin
Bel
```

The helper function file provides functions used by the statistics script, including:

```r
switchToPopLabels()
getSubtreeStats()
computePopEnrichment()
getCladeNage()
```

---

## 4B. Convert midpoint TRE files to STAT files

Scripts:

```text
scripts/create_stat_files_midpoint_parallel.sbatch
scripts/tre_to_stats_xipho_parallel.R
```

Job table reused from Step 3:

```text
infoTables/region_tree_job_table_midpoint.txt
```

Strategy:

One Slurm array task is submitted per trimmed region:

```text
573 array tasks
```

Each region task loops over all 50 MCMC iterations.

Input TRE files:

```text
argTreeFiles_midpoint/<region>/<region>.<start>-<end>.<iter>.tre.gz
```

Output STAT files:

```text
argStats_midpoint/<region>/<region>.<start>-<end>.<iter>.stat.gz
```

Example:

```text
argStats_midpoint/scaffold100-553/scaffold100-553.50001-1950000.1510.stat.gz
```

Submit:

```bash
module load R/4.2.1-foss-2022a

cd /scratch/amonc/xipho_revision/arg_processing/args_parallel
sbatch scripts/create_stat_files_midpoint_parallel.sbatch
```

Check final counts:

```bash
find argStats_midpoint -name "*.stat.gz" | wc -l
find argStats_midpoint -name "*.stat.gz" -size 0 | wc -l
```

Expected:

```text
28650 .stat.gz files
0 empty .stat.gz files
```

Spot-check:

```bash
find argStats_midpoint -name "*.stat.gz" | head -1
zcat $(find argStats_midpoint -name "*.stat.gz" | head -1) | head
```

---

# Statistics produced in Step 4

Each `.stat.gz` file includes:

```text
chrom
pos
TMRCA_all
TMRCAH_all
Tap_Xin_TMRCAH
Tap_Xin_RTH_original
Tap_Xin_RTH_prime
Xin_Bel_TMRCAH
Xin_Bel_RTH_original
Xin_Bel_RTH_prime
Tap_Xin_JCR
<pop>_RTH
<pop>_enrich
```

## Statistic definitions

### `TMRCA_all`

Age of the full tree/root.

### `TMRCAH_all`

Age of the youngest clade containing at least half of all sampled haplotypes.

### `Tap_Xin_TMRCAH`

Age of the youngest clade containing:

```text
>= half Tapajos haplotypes
>= half Xingu haplotypes
```

### `Tap_Xin_RTH_original`

```text
Tap_Xin_TMRCAH / TMRCA_all
```

### `Tap_Xin_RTH_prime`

```text
Tap_Xin_TMRCAH / TMRCAH_all
```

### `Xin_Bel_TMRCAH`

Age of the youngest clade containing:

```text
>= half Xingu haplotypes
>= half Belem haplotypes
```

### `Xin_Bel_RTH_original`

```text
Xin_Bel_TMRCAH / TMRCA_all
```

### `Xin_Bel_RTH_prime`

```text
Xin_Bel_TMRCAH / TMRCAH_all
```

### `Tap_Xin_JCR`

Joint Coalescence Ratio:

```text
Xin_Bel_TMRCAH / Tap_Xin_TMRCAH
```

Interpretation:

```text
> 1    Xingu-Belem joint coalescence is older than Tapajos-Xingu
< 1    Tapajos-Xingu joint coalescence is older than Xingu-Belem
≈ 1    Similar joint coalescence depths
```

### `<pop>_RTH`

For each population, the age of the youngest clade containing at least half of that population, divided by `TMRCAH_all`.

### `<pop>_enrich`

For each population, the maximum population enrichment observed in any subtree.

---

# Complete script list

## Shell scripts

```text
scripts/create_arg_block_region_file_xipho_parallel.sh
scripts/make_region_job_table.sh
```

## R scripts

```text
scripts/trim_arg_blocks_xipho_parallel.R
scripts/smc_to_bed_xipho_parallel.R
scripts/make_region_tree_job_table_midpoint_xipho_parallel.R
scripts/bed_to_tre_midpoint_xipho_parallel.R
scripts/tre_to_stats_xipho_parallel.R
```

## Slurm sbatch scripts

```text
scripts/create_bed_files_parallel.sbatch
scripts/create_tree_files_midpoint_parallel.sbatch
scripts/create_stat_files_midpoint_parallel.sbatch
```

## Helper scripts/files

```text
treeStatFunctions.R
```

or:

```text
scripts/treeStatFunctions.R
```

---

# Complete info table list

```text
infoTables/individual-species-key-xipho.txt
infoTables/ARGblock-coordinates.txt
infoTables/ARGblock-coordinates-trimmed.txt
infoTables/ARGblock-coordinates-info.txt
infoTables/region_job_table.txt
infoTables/region_tree_job_table_midpoint.txt
```

Optional/diagnostic table from earlier per-iteration testing:

```text
infoTables/region_iteration_job_table.txt
infoTables/region_iteration_job_table_untrimmed.txt
```

These are not required for the final per-region pipeline.

---

# Expected final output counts

## Input SMC files

```bash
find /scratch/amonc/xipho_revision/arg_processing/verified_smcs \
  -mindepth 2 -maxdepth 2 -type f -name "*_out.*.smc.gz" | wc -l
```

Expected:

```text
28750
```

## BED files

```bash
find argBedFiles -name "*.bed.gz" | wc -l
find argBedFiles -name "*.bed.gz.tbi" | wc -l
```

Expected:

```text
28750
28750
```

## TRE files

```bash
find argTreeFiles_midpoint -name "*.tre.gz" | wc -l
```

Expected:

```text
28650
```

## STAT files

```bash
find argStats_midpoint -name "*.stat.gz" | wc -l
```

Expected:

```text
28650
```

---

# Pipeline summary

```text
575 regions with complete SMC input
× 50 MCMC iterations
= 28,750 SMC files

After trimming:
573 retained regions
× 50 MCMC iterations
= 28,650 TRE files
= 28,650 STAT files
```

The BED step uses all 575 untrimmed regions, while the TRE and STAT steps use only the 573 retained trimmed regions.
