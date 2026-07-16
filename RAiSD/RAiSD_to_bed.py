# !/usr/bin/env python
# encoding: utf-8

"""
Goal: Convert RAiSD output to bed file format

To run:
python RAiSD_to_bed.py --input_file RAiSD_Report.vcf_run1.Chromosome_3_RagTag --output_file RAiSD_Report.vcf_run1.Chromosome_3_RagTag.bed

_______________________________________________________________________________
Copyright 2024 Andre E. Moncrieff. All rights reserved.

"""

import argparse

def process_custom_file(input_path, output_path):
    with open(input_path, 'r') as infile:
        lines = infile.readlines()

    # Get chromosome name from first line
    if not lines or not lines[0].startswith("//"):
        raise ValueError("First line must start with '// ' and contain the chromosome name.")
    chrom = lines[0].strip().lstrip('// ').strip()

    processed_lines = []

    # Process all remaining lines as data
    for line in lines[1:]:
        if not line.strip():
            continue  # Skip empty lines

        fields = line.strip().split('\t')
        try:
            pos = int(fields[0])  # First column is the position
            bed_start = pos - 1   # Convert to 0-based BED start
            bed_end = pos         # BED end is exclusive
        except ValueError:
            continue  # Skip malformed lines

        new_line = [chrom, str(bed_start), str(bed_end)] + fields[1:]
        processed_lines.append('\t'.join(new_line))

    # Write output
    with open(output_path, 'w') as outfile:
        for line in processed_lines:
            outfile.write(line + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert RAiSD tab-delimited report file to 0-based BED format.")
    parser.add_argument("--input_file", help="Input file path (custom format)")
    parser.add_argument("--output_file", help="Output file path (BED format)")

    args = parser.parse_args()
    process_custom_file(args.input_file, args.output_file)
