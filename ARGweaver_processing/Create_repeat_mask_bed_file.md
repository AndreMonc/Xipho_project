# Background:
Greg Thom, in Fall 2024, used RepeatMasker to soft-mask the new Xiphorhynchus elegans reference genome.
I can use the soft masking to create a bed file 

# First, check if there are any Ns or ns in the reference
```
grep -v "^>" xiph_elegans_ref.fa  | tr -d '\n' | grep -o "[Nn]" | wc -l
```

Ok, there aren't any

# Next how many lowercase bases are in the reference? Percentage of reference?
```
awk '
/^>/ {next}
{
  total += length($0)
  lower += gsub(/[acgt]/, "&")
}
END {
  printf "Lowercase bases: %d\nTotal bases: %d\nPercent lowercase: %.4f%%\n", lower, total, (lower/total)*100
}' xiph_elegans_ref.fa
```
Lowercase bases: 170758451
Total bases: 1120096857
Percent lowercase: 15.2450%

# No non actg or ACTG characters, Great!
```
grep -v "^>" xiph_elegans_ref.fa | grep -q '[^ACTGactg]' && echo "Non-ACTG found" || echo "Clean"
```
Clean


# ## Goal is to get bed file for all soft masked regions (acgt)
# AWK-based method to create bed file. This behaves as expected after testing.
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --output=repeat_bed_%J_stdout.txt
#SBATCH --error=repeat_bed_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=repeat_bed
#SBATCH --chdir=/scratch/amonc/xipho_revision/ref
#
#################################################

awk '
/^>/ {
    seqname = substr($0, 2)
    pos = 0
    in_lower = 0
    next
}
{
    for (i = 1; i <= length($0); i++) {
        base = substr($0, i, 1)

        if (base ~ /[acgt]/) {
            if (!in_lower) {
                start = pos
                in_lower = 1
            }
        } else {
            if (in_lower) {
                print seqname "\t" start "\t" pos
                in_lower = 0
            }
        }
        pos++
    }
}
END {
    if (in_lower) {
        print seqname "\t" start "\t" pos
    }
}
' xiph_elegans_ref.fa > repeat_regions.bed