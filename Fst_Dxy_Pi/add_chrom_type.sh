#!/bin/bash

set -euo pipefail

merged="pixy_merged_fst_dxy_pi.tsv"
sex_bed="sex_chrom_scaffolds.bed"
out="pixy_merged_fst_dxy_pi.chrom_type.tsv"

awk -F'\t' -v OFS='\t' '

# Read sex chromosome scaffold names from BED file
FNR==NR {
    if (FNR == 1) next   # skip BED header
    sex[$1] = 1
    next
}

# Add new column name to TSV header
FNR==1 {
    print $0, "chrom_type"
    next
}

# Annotate each row
{
    if ($1 in sex)
        print $0, "sex_chromosome"
    else
        print $0, "autosome"
}

' "$sex_bed" "$merged" > "$out"

echo "Wrote: $out"
echo ""

echo "Unique chromosome counts:"
awk -F'\t' '
NR==1 { next }

{
    all[$1] = 1

    if ($NF == "sex_chromosome")
        sex[$1] = 1
    else if ($NF == "autosome")
        auto[$1] = 1
}

END {
    n_all = n_sex = n_auto = 0

    for (x in all)  n_all++
    for (x in sex)  n_sex++
    for (x in auto) n_auto++

    print "all_chromosomes:\t" n_all
    print "sex_chromosomes:\t" n_sex
    print "autosomes:\t" n_auto
}
' "$out"

echo ""

echo "Window counts:"
awk -F'\t' '
NR==1 { next }

{
    total++

    if ($NF == "sex_chromosome")
        sex++
    else if ($NF == "autosome")
        auto++
}

END {
    print "total_windows:\t" total
    print "sex_windows:\t" sex
    print "autosome_windows:\t" auto
}
' "$out"
