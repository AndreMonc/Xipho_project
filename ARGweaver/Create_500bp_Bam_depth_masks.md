# Goal is to get BED depth maps for xipho individuals, for input to ARGweaver
# Set up conda environment
conda create -n mosdepth_env -c bioconda mosdepth
conda activate mosdepth_env

# Need SAMtools
wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
tar -xvjf samtools-1.21.tar.bz2
cd samtools-1.21
./configure
make
make install # doesn't work on cluster--error, but apparently that doesn't matter because samtools worked

# Create genome file .fai
/ddnA/work/a_monc/postdoc/samtools-1.21/samtools faidx xipho_elegans_ragtagRef_no_W.fa

# Create a bed file with 500 bp windows that will serve as the map for all calculations
bedtools makewindows -g /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_elegans_ragtagRef_no_W.fa.fai -w 500 > xipho_depth_windows.bed

# Calculate depth for files
https://github.com/brentp/mosdepth

# make directory for output
mkdir Depth_from_BAMs_500bp_windows

# HPC prompt example
#!/bin/bash
#PBS -A hpc_argweaver3
#PBS -l nodes=1:ppn=64
#PBS -l walltime=24:00:00
#PBS -q checkpt
#PBS -N bam_depth

cd /scratch/a_monc/postdoc/xipho_project/Depth_from_BAMs_500bp_windows

source activate mosdepth_env

mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Bel_xsFRC041_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Bel_xsFRC041_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Bel_xsGUR156_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Bel_xsGUR156_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Bel_xsLCA23_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Bel_xsLCA23_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Bel_xsMLV157_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Bel_xsMLV157_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Bel_xsRDP017_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Bel_xsRDP017_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Bel_xsREBIO002_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Bel_xsREBIO002_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Bel_xsTP32151_Toe /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Bel_xsTP32151_Toe_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Bel_xsTP37350_Toe /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Bel_xsTP37350_Toe_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Bel_xsTP38597_Toe /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Bel_xsTP38597_Toe_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Bel_xsTP51965_Toe /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Bel_xsTP51965_Toe_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Tap_xsA08267_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Tap_xsA08267_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Tap_xsBR163-028_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Tap_xsBR163-028_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Tap_xsBR163-145_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Tap_xsBR163-145_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Tap_xsBR163-212_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Tap_xsBR163-212_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Tap_xsMPDS1217_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Tap_xsMPDS1217_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Tap_xsMPDS1294_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Tap_xsMPDS1294_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Tap_xsMSF111_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Tap_xsMSF111_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Tap_xsPIME217_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Tap_xsPIME217_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Tap_xsSER013_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Tap_xsSER013_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Tap_xsTM005_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Tap_xsTM005_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Xin_xsCOUFT0424_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Xin_xsCOUFT0424_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Xin_xsFTA012_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Xin_xsFTA012_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Xin_xsGAPTO037_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Xin_xsGAPTO037_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Xin_xsGAPTO271_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Xin_xsGAPTO271_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Xin_xsGAPX047_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Xin_xsGAPX047_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Xin_xsMAYA066_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Xin_xsMAYA066_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Xin_xsMOP011_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Xin_xsMOP011_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Xin_xsMRJ498_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Xin_xsMRJ498_Mus_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Xin_xsTP36025_Toe /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Xin_xsTP36025_Toe_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Xin_xsTP36276_Toe /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Xin_xsTP36276_Toe_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Xin_xsTP48649_Toe /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Xin_xsTP48649_Toe_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Xin_xsTP81164_Toe /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Xin_xsTP81164_Toe_final.bam
mosdepth --no-per-base --fast-mode --by /ddnA/work/a_monc/postdoc/refs/Xiphorhynchus_elegans/xipho_depth_windows.bed --threads 3 Xin_xsUHE455_Mus /scratch/a_monc/postdoc/xipho_project/snpArcher_wd_no_Wchrom/results/xipho_elegans_ragtagRef_no_W.fa/bams/Xin_xsUHE455_Mus_final.bam

# The above outputs bed.gz files
# counting the line numbers on the bed.gz does not give actual line numbers of bed files!!!!!!! Confusion and lost time here.
# Need to gunzip all the bed.gz files

for file in *bed.gz; do
    gunzip $file
done

for file in *bed; do
    echo $file
done

# Then count lines for each bed file and make sure the same. Good to go. 2240367 lines.
for file in *.bed; do
    wc -l $file
done

### Ok, for ARGweaver I'm going to use individual depth masks based on these bam-informed depth-interval bed files
# I filtered windows to depth <5 and >50

awk '($4 < 5 || $4 > 50) {print $0}' Bel_xsFRC041_Mus.regions.bed > Bel_xsFRC041_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Bel_xsGUR156_Mus.regions.bed > Bel_xsGUR156_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Bel_xsLCA23_Mus.regions.bed > Bel_xsLCA23_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Bel_xsMLV157_Mus.regions.bed > Bel_xsMLV157_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Bel_xsRDP017_Mus.regions.bed > Bel_xsRDP017_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Bel_xsREBIO002_Mus.regions.bed > Bel_xsREBIO002_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Bel_xsTP32151_Toe.regions.bed > Bel_xsTP32151_Toe_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Bel_xsTP37350_Toe.regions.bed > Bel_xsTP37350_Toe_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Bel_xsTP38597_Toe.regions.bed > Bel_xsTP38597_Toe_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Bel_xsTP51965_Toe.regions.bed > Bel_xsTP51965_Toe_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Tap_xsA08267_Mus.regions.bed > Tap_xsA08267_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Tap_xsBR163-028_Mus.regions.bed > Tap_xsBR163-028_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Tap_xsBR163-145_Mus.regions.bed > Tap_xsBR163-145_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Tap_xsBR163-212_Mus.regions.bed > Tap_xsBR163-212_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Tap_xsMPDS1217_Mus.regions.bed > Tap_xsMPDS1217_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Tap_xsMPDS1294_Mus.regions.bed > Tap_xsMPDS1294_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Tap_xsMSF111_Mus.regions.bed > Tap_xsMSF111_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Tap_xsPIME217_Mus.regions.bed > Tap_xsPIME217_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Tap_xsSER013_Mus.regions.bed > Tap_xsSER013_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Tap_xsTM005_Mus.regions.bed > Tap_xsTM005_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Xin_xsCOUFT0424_Mus.regions.bed > Xin_xsCOUFT0424_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Xin_xsFTA012_Mus.regions.bed > Xin_xsFTA012_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Xin_xsGAPTO037_Mus.regions.bed > Xin_xsGAPTO037_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Xin_xsGAPTO271_Mus.regions.bed > Xin_xsGAPTO271_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Xin_xsGAPX047_Mus.regions.bed > Xin_xsGAPX047_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Xin_xsMAYA066_Mus.regions.bed > Xin_xsMAYA066_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Xin_xsMOP011_Mus.regions.bed > Xin_xsMOP011_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Xin_xsMRJ498_Mus.regions.bed > Xin_xsMRJ498_Mus_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Xin_xsTP36025_Toe.regions.bed > Xin_xsTP36025_Toe_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Xin_xsTP36276_Toe.regions.bed > Xin_xsTP36276_Toe_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Xin_xsTP48649_Toe.regions.bed > Xin_xsTP48649_Toe_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Xin_xsTP81164_Toe.regions.bed > Xin_xsTP81164_Toe_depthMask.bed
awk '($4 < 5 || $4 > 50) {print $0}' Xin_xsUHE455_Mus.regions.bed > Xin_xsUHE455_Mus_depthMask.bed


# Merge book-ended intervals
bedtools merge -d 0 -i Bel_xsFRC041_Mus_depthMask.bed > Bel_xsFRC041_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Bel_xsGUR156_Mus_depthMask.bed > Bel_xsGUR156_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Bel_xsLCA23_Mus_depthMask.bed > Bel_xsLCA23_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Bel_xsMLV157_Mus_depthMask.bed > Bel_xsMLV157_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Bel_xsRDP017_Mus_depthMask.bed > Bel_xsRDP017_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Bel_xsREBIO002_Mus_depthMask.bed > Bel_xsREBIO002_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Bel_xsTP32151_Toe_depthMask.bed > Bel_xsTP32151_Toe_depthMask_merged.bed
bedtools merge -d 0 -i Bel_xsTP37350_Toe_depthMask.bed > Bel_xsTP37350_Toe_depthMask_merged.bed
bedtools merge -d 0 -i Bel_xsTP38597_Toe_depthMask.bed > Bel_xsTP38597_Toe_depthMask_merged.bed
bedtools merge -d 0 -i Bel_xsTP51965_Toe_depthMask.bed > Bel_xsTP51965_Toe_depthMask_merged.bed
bedtools merge -d 0 -i Tap_xsA08267_Mus_depthMask.bed > Tap_xsA08267_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Tap_xsBR163-028_Mus_depthMask.bed > Tap_xsBR163-028_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Tap_xsBR163-145_Mus_depthMask.bed > Tap_xsBR163-145_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Tap_xsBR163-212_Mus_depthMask.bed > Tap_xsBR163-212_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Tap_xsMPDS1217_Mus_depthMask.bed > Tap_xsMPDS1217_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Tap_xsMPDS1294_Mus_depthMask.bed > Tap_xsMPDS1294_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Tap_xsMSF111_Mus_depthMask.bed > Tap_xsMSF111_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Tap_xsPIME217_Mus_depthMask.bed > Tap_xsPIME217_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Tap_xsSER013_Mus_depthMask.bed > Tap_xsSER013_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Tap_xsTM005_Mus_depthMask.bed > Tap_xsTM005_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Xin_xsCOUFT0424_Mus_depthMask.bed > Xin_xsCOUFT0424_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Xin_xsFTA012_Mus_depthMask.bed > Xin_xsFTA012_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Xin_xsGAPTO037_Mus_depthMask.bed > Xin_xsGAPTO037_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Xin_xsGAPTO271_Mus_depthMask.bed > Xin_xsGAPTO271_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Xin_xsGAPX047_Mus_depthMask.bed > Xin_xsGAPX047_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Xin_xsMAYA066_Mus_depthMask.bed > Xin_xsMAYA066_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Xin_xsMOP011_Mus_depthMask.bed > Xin_xsMOP011_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Xin_xsMRJ498_Mus_depthMask.bed > Xin_xsMRJ498_Mus_depthMask_merged.bed
bedtools merge -d 0 -i Xin_xsTP36025_Toe_depthMask.bed > Xin_xsTP36025_Toe_depthMask_merged.bed
bedtools merge -d 0 -i Xin_xsTP36276_Toe_depthMask.bed > Xin_xsTP36276_Toe_depthMask_merged.bed
bedtools merge -d 0 -i Xin_xsTP48649_Toe_depthMask.bed > Xin_xsTP48649_Toe_depthMask_merged.bed
bedtools merge -d 0 -i Xin_xsTP81164_Toe_depthMask.bed > Xin_xsTP81164_Toe_depthMask_merged.bed
bedtools merge -d 0 -i Xin_xsUHE455_Mus_depthMask.bed > Xin_xsUHE455_Mus_depthMask_merged.bed

# Final individual mask file. Did some simple excel work to create this file. Input for ARGweaver.
/ddnA/work/a_monc/postdoc/xipho_project/Depth_from_BAMs_500bp_windows/ind_mask_file.txt