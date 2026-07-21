#!/bin/bash
set -euo pipefail

INPUT="infoTables/ARGblock-coordinates.txt"
OUT="infoTables/region_job_table.txt"

echo -e "job_id\tregion\tchromosome\tstart\tend" > "$OUT"

job_id=0

tail -n +2 "$INPUT" | while read -r region chromosome start end
do
    job_id=$((job_id + 1))
    echo -e "${job_id}\t${region}\t${chromosome}\t${start}\t${end}" >> "$OUT"
done

echo "Wrote: $OUT"
wc -l "$OUT"