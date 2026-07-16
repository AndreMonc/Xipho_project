# Processing .smc file output from ARGweaver
Starting on 9 June 2026
Summary at beginning of processing: 631 windows out of an original 634 windows successfully ran in ARGweaver, so 99.5%%
Then, I removed 56 windows that had >50% of sites masked in ARGweaver. This information is provided in the output .log files from ARGweaver.

This left me with 575 windows for calculating ARG statistics.

## Basing my processing of the ARGs on this page (smc to stats pipeline):
https://github.com/CshlSiepelLab/bird_capuchino_analysis/tree/master/ARG_analysis

## All of my smc files have this naming format
## 209 is the window number out 634
## 1630 is the MCMC iteration number
scaffold17-209_out.1630.smc.gz

# Rather than organizing all my smc files into one folder (too many files once I start processing different MCMC iterations)
# I have all regions organized in their own folders.

# Once I trim the 50-kb overlaps off of each region, I will calculate number of local trees and tree spans per region. I can then create a summary across all 575 regions.

cd /scratch/amonc/xipho_revision/arg_processing
git clone https://github.com/CshlSiepelLab/bird_capuchino_analysis.git
cd bird_capuchino_analysis/ARG_analysis/smc_to_stats
ls -lh
cat README.md

# Created species key file
individual-species-key-xipho.txt


# Step 1.1 - Create block region file

# Created this script as create_arg_block_region_file_xipho.sh

bash create_arg_block_region_file_xipho.sh

# Step 1.2 - Filter and trim ARG blocks
I used the following parameters for trimming coordinates:
mcmc_iter      <- 2000
min_block_len  <- 101000
segment_len    <- 1000
trim_len       <- 50000

```
module load R/4.2.1-foss-2022a
Rscript trim_arg_blocks_xipho.R
```
# Output from step 1.2:
number of analyzed scaffolds is 141
number of ARG blocks with at least 2000 MCMC iterations, of size at least 101000 bp is 573
total length of these ARG blocks is: 998057441
number of ARG blocks with at least 2000 MCMC iterations, but shorter than 101000 bp is 2
total length of these ARG blocks is: 164815


## Step 2 - Generate ARG bed files
number of jobs in parallel = 20
mcmc iterations = 2000
min ARG block length (101,000)

module load R/4.2.1-foss-2022a
Rscript make_command_lines_for_beds_xipho.R # using inputFile   <- "./infoTables/ARGblock-coordinates.txt"

## Step 3 - Generate tree files
module load R/4.2.1-foss-2022a
Rscript make_command_lines_for_trees_midpoint_xipho.R

## Step 4 - Generate stat files
module load R/4.2.1-foss-2022a
Rscript make_command_lines_for_stats_midpoint_xipho.R 
sbatch run_tre2stat_midpoint.sbatch


