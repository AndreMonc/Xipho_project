# Time to process the smc files output from ARGweaver
# Starting on 20 March 2025
# Summary at beginning of processing: 592 windows out of an original 604 windows successfully ran in ARGweaver, so 98.0%

# Basing my processing of the ARGs on this page (smc to stats pipeline):
https://github.com/CshlSiepelLab/bird_capuchino_analysis/tree/master/ARG_analysis

# All of my smc files are already in the correct format this time around:
Chromosome6-482_out.2000.smc.gz

# Put all ARGweaver output files (.SMC.GZ ) in a single smc source dir:
/scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/output/all_smc_files

# Put all ARGweaver output log files (.LOG) in a single log source dir:
/scratch/a_monc/postdoc/xipho_project/ARGweaver/Runs/output/log_files

# Create directory infoTables/ in current directory:
/scratch/a_monc/postdoc/xipho_project/ARGweaver/ARG_processing/infoTables

# smc2bed program
/ddnA/work/a_monc/postdoc/ARGweaver/bin/smc2bed

individual-species-key-xipho.txt

# Set bash variable for directory

scriptDir='/scratch/a_monc/postdoc/xipho_project/ARGweaver/ARG_processing'

# Step 1 - Filter and trim ARG blocks
# step 1.1
--------
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -q single
#PBS -N create_arg_block_region_file

cd /scratch/a_monc/postdoc/xipho_project/ARGweaver/ARG_processing

bash /ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ARG_processing/create_arg_block_region_file.sh


-----
number of analyzed scaffolds is 48
number of ARG blocks with at least 2000 MCMC iterations, of size at least 110000  bp is 589
total length of these ARG blocks is: 1131205430
number of ARG blocks with at least 2000 MCMC iterations, but shorter than 110000  bp is 3
total length of these ARG blocks is: 47932
-----


# step 1.2
# Settings for filtering and trimming windows:
(1) MCMC sampling iteration for filtering = 2000; (2) minimum ARG block length for filtering = 110,000; (3) window length = 10,000; (4) trim size in each end of an ARG block (default=50,0000).

module load r/4.2.0/gcc-11.2.0

--------
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -q single
#PBS -N trim_arg_blocks

cd /scratch/a_monc/postdoc/xipho_project/ARGweaver/ARG_processing

Rscript /ddnA/work/a_monc/postdoc/xipho_project/ARGweaver/ARG_processing/trim_arg_blocks.R

-------
# Got an error when running above
Error in if (i == 1 || ordered_args_df[i, "blockIndex"] != ordered_args_df[i -  : 
  missing value where TRUE/FALSE needed
Calls: rbind -> divideBlocksToRegionsBySegments
In addition: Warning message:
In divideBlocksToRegionsBySegments(arg_blocks_per_scaffold, trim_len,  :
  NAs introduced by coercion
Execution halted

---------

# Updating the trim_arg_blocks.R file

#### IMPORTANT ### Ok, I had to update line 37 in the trim_arg_blocks.R script, changing sub(".*\\.","",args_df$blockName)) to sub(".*[-]","",args_df$blockName))


# Step 2 - Generate ARG bed files
number of jobs in parallel = 1 (had issues with more jobs when testing on chromosome 1)
mcmc iterations = 2000
min ARG block length (120,000)


--------


--------
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -q single
#PBS -N create_bed_files

module load r/4.2.0/gcc-11.2.0

cd /scratch/a_monc/postdoc/xipho_project/ARGweaver/ARG_processing

bash ./create_bed_files.sh


# Step 3 - Generate Tree Files
make command line for trees file:
    set interval between sampled trees = 100
    
----
# Basically, had an issue with the formatting of the scaffold/chromosome names in my original ARG files, which got carried down the line. 
# The program bed_to_tre.R cannot handle the variety of scaffold/chromosome names in my files, so I need to update bed_to_tre.R
# I made the following prompt for chatgpt, which helped me to find the solution

Here is some code of mine for processing raw .bed.gz bedFile names to give a scaffold name:

scaffold <- sub(".*/","",bedFile)
scaffold <- sub("[.]bed[.]gz","",scaffold)
scaffold <- sub("_out\\..*","",scaffold)
scaffold <- sub("\\..*","",scaffold)

I need to update this code to handle somewhat variable bedFile names. For example here are some bedFile names: 
Chromosome14-141_out.2000.bed.gz
scaffold225-570_out.2000.bed.gz
Chromosome3-339_out.2000.bed.gz
scaffold219-569_out.2000.bed.gz
Chromosome2-207_out.2000.bed.gz.tbi
ChromosomeZ-600_out.2000.bed.gz

For the bedFile names that start with the word “Chromosome”, I need the code to produce scaffold names in the format Chromosome_[num or text]_RagTag (e.g., Chromosome_14_RagTag or Chromosome_Z_RagTag)

For bedFile names that start with the word “scaffold”, I need the code to produce scaffold names in the format scaffold_num (e.g., scaffold_219)

Can you help me update my code to do this?

# Code from chatGPT to add to bed_to_tre.R for the scaffold name processing (which will then allow tabix to properly query regions)

# Extract the filename from the full path
scaffold <- sub(".*/", "", bedFile)

# Check if the filename starts with "Chromosome" or "scaffold"
if (grepl("^Chromosome", scaffold)) {
    # For files starting with "Chromosome", format as "Chromosome_[num or text]_RagTag"
    scaffold <- sub("Chromosome([A-Za-z0-9]+)-.*", "Chromosome_\\1_RagTag", scaffold)
} else if (grepl("^scaffold", scaffold)) {
    # For files starting with "scaffold", format as "scaffold_num"
    scaffold <- sub("scaffold([0-9]+)-.*", "scaffold_\\1", scaffold)
}

# Now 'scaffold' will be in the desired format based on the bedFile name
print(scaffold)

--------------
# Also a little update to remove the after hyphen part of the chromosome names
# Process the filename to extract the scaffold name
----
# chat gpt prompt
Here is some code for processing .bed.gz filenames:
scaffold <- sub(".*/","",bedFile)
scaffold <- sub("[.]bed[.]gz","",scaffold)
scaffold <- sub("_out\\..*","",scaffold)
scaffold <- sub("\\..*","",scaffold)

How does the code process the filenames below?
scaffold219-569_out.2000.bed.gz
Chromosome8-540_out.2000.bed.gz
------
scaffold <- sub(".*/","",bedFile)  # Remove directory path
scaffold <- sub("[.]bed[.]gz","",scaffold)  # Remove .bed.gz extension
scaffold <- sub("_out\\..*","",scaffold)  # Remove _out.* part
scaffold <- sub("-.*","",scaffold)  # Remove everything after the hyphen

# The resulting scaffold will be in the desired format
print(scaffold)

----------------------

#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -q single
#PBS -N create_tree_files

module load r/4.2.0/gcc-11.2.0

cd /scratch/a_monc/postdoc/xipho_project/ARGweaver/ARG_processing

bash ./create_tree_files.sh

# Step 4 - Generate Stat files

# Needed to update for the variable files names I have (both chrom and scaff)
# ChatGPT prompt:
Here are some paths, each stored as a variable named treefile:
./argTreeFiles/Chromosome19.1950001-3850000.2000.tre.gz
./argTreeFiles/Chromosome5.79850001-79880000.2000.tre.gz
./argTreeFiles/scaffold225.50001-80000.2000.tre.gz
./argTreeFiles/scaffold216.50001-120000.2000.tre.gz

I need code to process these treeFile variables that takes into account whether the .tre.gz file starts with Chromosome or scaffold. For example, if it starts with Chromosome, then:
statsFile <- sub(“.*/Chromsome, “Chromosome”, treeFile)

Or, if it starts with scaffold:
statsFile <- sub(“.*/scaffold, “scaffold”, treeFile)

How do I do that?

----------------------
# Answer, adding to tre_to_stats.R file

if (grepl("/Chromosome", treeFile)) {
  statsFile <- sub(".*/Chromosome", "Chromosome", treeFile)
} else if (grepl("/scaffold", treeFile)) {
  statsFile <- sub(".*/scaffold", "scaffold", treeFile)
} else {
  stop("Invalid file name: does not start with Chromosome or scaffold")
}
--------------------

#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=48:00:00
#PBS -q single
#PBS -N create_stat_files

module load r/4.2.0/gcc-11.2.0

cd /scratch/a_monc/postdoc/xipho_project/ARGweaver/ARG_processing

bash ./create_stat_files.sh


-----------------------
# the Create State Files step takes a LONG time!!
Ok, I am going to use a tmux session to run this stat-generating step

tmux new -s argstats #give whatever session name you want, created on mike1 head node

qsub -I -l walltime=72:00:00,nodes=1:ppn=20 -A hpc_argweaver3

module load r/4.2.0/gcc-11.2.0

bash ./create_stat_files.sh

# Ctrl-b d # use to detach from tmux window/session while program within continues to run

tmux a -t argstats #attaches to tmux session with specified name
Ctrl-b c #makes new window in tmux session
Ctrl-b n #switches between windows in tmux session

---------------------
tmux new -s argstats_pt2 # started on mike2 head node

qsub -I -l walltime=72:00:00,nodes=1:ppn=20 -A hpc_argweaver3

module load r/4.2.0/gcc-11.2.0

bash ./create_stat_files_pt2.sh
---------------------
tmux new -s argstats_pt3

qsub -I -l walltime=72:00:00,nodes=1:ppn=20 -A hpc_argweaver3

module load r/4.2.0/gcc-11.2.0

bash ./create_stat_files_pt3.sh
---------------------
tmux new -s argstats_pt4

qsub -I -l walltime=48:00:00,nodes=1:ppn=20 -A hpc_argweaver3

module load r/4.2.0/gcc-11.2.0

bash ./create_stat_files_pt4.sh
---------------------
tmux new -s argstats_pt5

qsub -I -l walltime=48:00:00,nodes=1:ppn=20 -A hpc_argweaver3

module load r/4.2.0/gcc-11.2.0

bash ./create_stat_files_pt5.sh
---------------------
tmux new -s argstats_pt6

qsub -I -l walltime=48:00:00,nodes=1:ppn=20 -A hpc_argweaver3

module load r/4.2.0/gcc-11.2.0

bash ./create_stat_files_pt6.sh

---------------------
kill all sessions
tmux kill-session -t argstats_pt2
tmux kill-session -t argstats_pt3
tmux kill-session -t argstats_pt4
tmux kill-session -t argstats_pt5
tmux kill-session -t argstats_pt6



# Step 5 - Generate Stat BED files

Just made a couple of minor changes to window_stats.R related to rounding (and remove the age and enrich stat variables)
and "pop" to "pos" in the line below:
message <- paste(scaffold," gap in end of ARG block",i,"( pos",arg_block_table$endPos[i],"-",pop,")")

-----------
tmux new -s argwinstats

qsub -I -l walltime=72:00:00,nodes=1:ppn=20 -A hpc_argweaver3

module load r/4.2.0/gcc-11.2.0

Rscript make_all_window_stats.R


tmux a -t argwinstats #attaches to tmux session with specified name
Ctrl-b c #makes new window in tmux session
Ctrl-b n #switches between windows in tmux session

---------------

# Sweet, above worked. Now have 48 scaffolds.
# Why have 49 scaffolds with recomb data and only 48 with ARG data? Check
scaffold_227 failed in the ARGweaver run. There was only one VCF 2mb window for that scaffold, so it did not continue to later steps.

# Next, combine all beds into one:
head -n 1 Chromosome1.stat.bed > xipho_2025_argstats.txt; tail -n +2 -q *.bed >> xipho_2025_argstats.txt