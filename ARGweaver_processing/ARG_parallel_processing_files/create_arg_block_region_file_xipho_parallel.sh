#!/bin/bash
set -euo pipefail

SMC_DIR="/scratch/amonc/xipho_revision/arg_processing/verified_smcs"
OUT="infoTables/ARGblock-coordinates.txt"

echo -e "region\tchromosome\tstart\tend" > "$OUT"

find "$SMC_DIR" -mindepth 2 -maxdepth 2 -type f -name "*_out.1510.smc.gz" | sort | while read -r f
do
    base=$(basename "$f")
    region=${base%_out.1510.smc.gz}

    set +o pipefail
    region_line=$(gzip -dc "$f" | awk '$1=="REGION" {print $2"\t"$3"\t"$4; exit}')
    set -o pipefail

    if [[ -z "$region_line" ]]; then
        echo "WARNING: no REGION line found in $f" >&2
        continue
    fi

    echo -e "${region}\t${region_line}" >> "$OUT"
done

echo "Wrote: $OUT"
wc -l "$OUT"