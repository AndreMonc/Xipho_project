# May 21 2026, for Xiphorhynchus revision
# Running rCNV on the OSCER HPC

# Installed recommend newer GitHub version of rCNV package ‘1.4.900’
(not from R)


#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=rCNV
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=rjob_%j.out
#SBATCH --error=rjob_%j.err

module load R/4.3.1

Rscript myscript.R

# For running on cluster, see:
run_rCNV_array.sh
run_rCNV_one_chrom.R

# Results written here:
/scratch/amonc/xipho_revision/rCNV/results

# Merge the summary files for each scaffold 
- This will provide info on which scaffolds failed and which had usable data (bigger scaffolds)
```
Rscript merge_rCNV_summaries.R
```
This produced a genome-wide summary file for 396 scaffolds that produced summary files.

# Now, use "n_WGS_rows" column to remove scaffolds with values of NA or 1
awk -F'\t' '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i=="n_WGS_rows") col=i
    }
    print
    next
}
$col != "" && $col != "NA" && $col != 1
' genomewide.rCNV.summary.tsv > genomewide.rCNV.summary.filtered.tsv

# Merging all the allele.info.WGS() function output (separate scaffold-specific files):
```
awk -F'\t' '
NR==FNR {
    if (FNR==1) next
    keep[$1]=1
    next
}
FNR==1 {
    file=FILENAME
    sub(/^.*\//, "", file)
    sub(/\.allele_info_WGS\.full\.tsv$/, "", file)

    if (!(file in keep)) nextfile

    if (!printed_header) {
        print
        printed_header=1
    }
    next
}
{
    print
}
' genomewide.rCNV.summary.filtered.tsv *.allele_info_WGS.full.tsv \
> genomewide.rCNV.allele_info_WGS.full.filtered.tsv
```
# Counting the number of unique scaffolds in genomewide.rCNV.allele_info_WGS.full.filtered.tsv
tail -n +2 genomewide.rCNV.allele_info_WGS.full.filtered.tsv | cut -f1 | sort -u | wc -l

result: 88 (this checks out, yay)


# Remove last column (unecessary duplicate of first column--this will reduce file size)
awk 'BEGIN{FS=OFS="\t"} {
    NF--
    print
}' genomewide.rCNV.allele_info_WGS.full.filtered.tsv > genomewide.rCNV.allele_info_WGS.full.filtered.nochrom.tsv


# Sort the genomewide file by chrom and position
{
    head -n 1 genomewide.rCNV.allele_info_WGS.full.filtered.nochrom.tsv
    tail -n +2 genomewide.rCNV.allele_info_WGS.full.filtered.nochrom.tsv | \
        sort -t $'\t' -k1,1 -k2,2n
} > genomewide.rCNV.allele_info_WGS.full.filtered.nochrom.sorted.tsv

# Sort the 10-kb windows bed file by chrom and position
sort -t $'\t' -k1,1 -k2,2n windows.bed > windows.sorted.bed


# Quick script to calculate the percentage of duplicated sites per 10-kb window across the genome
# Correctly handles partial windows at end of scaffolds
```
awk -F'\t' '
BEGIN{OFS="\t"}

FNR==NR {
    key=$1 FS $2
    scaffold[key]=$1
    start[key]=$2
    end[key]=$3
    order[++n]=key
    dup[key]=0
    nondup[key]=0
    next
}

FNR==1 {
    for(i=1;i<=NF;i++) {
        if($i=="CHROM") chrom_col=i
        if($i=="POS") pos_col=i
        if($i=="dup.stat") dup_col=i
    }
    next
}

{
    chrom=$chrom_col
    pos=$pos_col
    win_start=int((pos-1)/10000)*10000
    key=chrom FS win_start

    if(key in dup) {
        if($dup_col=="duplicated") dup[key]++
        else if($dup_col=="non-duplicated") nondup[key]++
    }
}

END{
    print "scaffold","start","end","dup.sites","nondup.sites","total.sites","perc.dup.sites"

    for(i=1;i<=n;i++) {
        key=order[i]
        d=dup[key]
        nd=nondup[key]
        total=d+nd

        if(total>0) perc=100*d/total
        else perc="NA"

        print scaffold[key],start[key],end[key],d,nd,total,perc
    }
}
' windows.sorted.bed genomewide.rCNV.allele_info_WGS.full.filtered.nochrom.sorted.tsv \
> genomewide.rCNV.allele_info_WGS.10kb.tsv
```
