#!/bin/bash
set -euo pipefail

IN="/scratch/amonc/xipho_revision/arg_processing/arg_based_fst/xipho_arg_based_branch_fst.iter2000.tsv"
OUT="/scratch/amonc/xipho_revision/arg_processing/arg_based_fst/xipho_arg_based_branch_fst.iter2000.cleaned.tsv"

awk -v OFS='\t' '
NR == 1 {
    print "region", "coordinate", "Fst_belem_tapajos", "Fst_belem_xingu", "Fst_tapajos_xingu"
    next
}
{
    region = $1

    # Convert scaffold100-553 -> scaffold_100
    sub(/-.*/, "", region)
    sub(/^scaffold/, "scaffold_", region)

    fst_bt = $4
    fst_bx = $5
    fst_tx = $6

    # Convert negative Fst values to zero
    if (fst_bt < 0) fst_bt = 0
    if (fst_bx < 0) fst_bx = 0
    if (fst_tx < 0) fst_tx = 0

    print region, $3, \
          sprintf("%.5f", fst_bt), \
          sprintf("%.5f", fst_bx), \
          sprintf("%.5f", fst_tx)
}
' "$IN" > "$OUT"

echo "Cleaned file written to:"
echo "$OUT"
