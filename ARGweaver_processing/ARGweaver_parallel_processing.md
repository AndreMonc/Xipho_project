# Processing .smc file outputs from ARGweaver
# Across 50 MCMC iterations - 1510 through 2000

Starting on 12 June 2026
Summary at beginning of processing: 631 windows out of an original 634 windows successfully ran in ARGweaver, so 99.5%%
Then, I removed 56 windows that had >50% of sites masked in ARGweaver. This information is provided in the output .log files from ARGweaver.

This left me with 575 windows for calculating ARG statistics.

## Basing my processing of the ARGs on this page (smc to stats pipeline):
https://github.com/CshlSiepelLab/bird_capuchino_analysis/tree/master/ARG_analysis

## All of my smc files have this naming format
## 209 is the window number out of 634
## 1630 is the MCMC iteration number
scaffold17-209_out.1630.smc.gz

# Rather than organizing all my smc files into one folder (too many files once I start processing different MCMC iterations)
# I have all regions organized in their own folders here:

/scratch/amonc/xipho_revision/arg_processing/verified_smcs

# Basic folder structure to start:

cd /scratch/amonc/xipho_revision/arg_processing

mkdir -p args_parallel/{scripts,logs,infoTables,argBedFiles,argTreeFiles_midpoint,argStats_midpoint,tmp}

# Created species key file
individual-species-key-xipho.txt
/scratch/amonc/xipho_revision/arg_processing/args_parallel/infoTables/individual-species-key-xipho.txt

# What I'm working with:
575 regions × 50 iterations = 28,750 smc.gz files

# Step 1 - Create block region files

# Created this script as create_arg_block_region_file_xipho.sh
```
chmod +x scripts/create_arg_block_region_file_xipho_parallel.sh
./scripts/create_arg_block_region_file_xipho_parallel.sh

module load R/4.2.1-foss-2022a

cd /scratch/amonc/xipho_revision/arg_processing/args_parallel
Rscript scripts/trim_arg_blocks_xipho_parallel.R
```
```
575 input regions
2 too short after filtering
573 retained trimmed regions
```

## Step 2 - Generate ARG bed files
cd /scratch/amonc/xipho_revision/arg_processing/args_parallel

chmod +x scripts/make_region_job_table.sh
./scripts/make_region_job_table.sh

sbatch scripts/create_bed_files_parallel.sbatch

## Step 3 -- Generate ARG tree files
module load R/4.2.1-foss-2022a
Rscript scripts/make_region_tree_job_table_midpoint_xipho_parallel.R

sbatch scripts/create_tree_files_midpoint_parallel.sbatch


## Step 4 - Generate stat files
module load R/4.2.1-foss-2022a
Rscript make_command_lines_for_stats_midpoint_xipho.R 
sbatch scripts/create_stat_files_midpoint_parallel.sbatch


# Testing fst calculation from genealogies

Working directory:
/scratch/amonc/xipho_revision/arg_processing/bird_capuchino_analysis/ARG_analysis/testing_fst


conda create -n arg_fst_env python=3.11 -y
conda activate arg_fst_env

conda install -c conda-forge numpy tskit -y

# Versions in use
numpy v2.4.6
tskit v1.0.3

pip install git+https://github.com/tskit-dev/tsconvert.git

pip install newick # Successfully installed newick-1.11.0
conda install -c conda-forge pandas -y


python - <<'PY'
import tsconvert

newick = "((A:1,B:1):2,C:3);"

ts = tsconvert.from_newick(
    newick,
    span=1,
    min_edge_length=1e-6
)

print(ts)
print("samples:", ts.samples())

for n in ts.samples():
    node = ts.node(n)
    print(n, node.metadata)
PY