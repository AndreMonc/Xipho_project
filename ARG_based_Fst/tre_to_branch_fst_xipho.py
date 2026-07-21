#!/usr/bin/env python

import argparse
import gzip
import re
from pathlib import Path

import pandas as pd
import tsconvert

BASE_DIR = Path("/scratch/amonc/xipho_revision/arg_processing")
WORK_DIR = BASE_DIR / "arg_based_fst"

MAP_FILE = WORK_DIR / "individual-species-key-xipho.txt"
OUT_DIR = WORK_DIR / "per_region_branch_fst"

OUT_DIR.mkdir(parents=True, exist_ok=True)

header = [
    "region",
    "tree_file",
    "coordinate",
    "Fst_belem_tapajos",
    "Fst_belem_xingu",
    "Fst_tapajos_xingu",
]

popmap = pd.read_csv(MAP_FILE, sep="\t")


def count_trees(tree_file):
    with gzip.open(tree_file, "rt") as f:
        return sum(1 for line in f if line.strip())


def output_is_complete(out_file, tree_file):
    if not out_file.exists() or out_file.stat().st_size == 0:
        return False

    try:
        with open(out_file, "r") as f:
            out_header = f.readline().rstrip("\n").split("\t")
            if out_header != header:
                return False
            n_rows = sum(1 for line in f if line.strip())
    except Exception:
        return False

    return n_rows == count_trees(tree_file)


def branch_fst_for_line(line, region, tree_file):
    coord, newick = line.strip().split(maxsplit=1)

    # Remove ARGweaver NHX comments
    newick = re.sub(r"\[&&NHX:[^\]]+\]", "", newick)

    ts = tsconvert.from_newick(
        newick,
        span=1,
        min_edge_length=1e-6,
    )

    samples = {ts.node(n).metadata["name"]: int(n) for n in ts.samples()}

    missing = sorted(set(popmap["Individual_hap"]) - set(samples))
    if missing:
        raise ValueError(
            f"Missing samples in {tree_file} at coordinate {coord}: "
            + ", ".join(missing)
        )

    belem = [
        samples[x]
        for x in popmap.loc[popmap["Species"] == "belem", "Individual_hap"]
    ]
    tapajos = [
        samples[x]
        for x in popmap.loc[popmap["Species"] == "tapajos", "Individual_hap"]
    ]
    xingu = [
        samples[x]
        for x in popmap.loc[popmap["Species"] == "xingu", "Individual_hap"]
    ]

    return [
        region,
        tree_file.name,
        coord,
        ts.Fst([belem, tapajos], indexes=[(0, 1)], mode="branch")[0],
        ts.Fst([belem, xingu], indexes=[(0, 1)], mode="branch")[0],
        ts.Fst([tapajos, xingu], indexes=[(0, 1)], mode="branch")[0],
    ]


def process_tree_file(tree_file):
    tree_file = Path(tree_file)

    if not tree_file.exists():
        raise FileNotFoundError(f"Tree file does not exist: {tree_file}")

    region = tree_file.parent.name
    out_file = OUT_DIR / f"{region}.branch_fst.iter2000.tsv"
    tmp_file = OUT_DIR / f"{region}.branch_fst.iter2000.tsv.tmp"

    if output_is_complete(out_file, tree_file):
        print(f"Skipping completed region: {region}", flush=True)
        return

    print(f"Processing region {region}: {tree_file}", flush=True)

    with gzip.open(tree_file, "rt") as f, open(tmp_file, "w") as out:
        out.write("\t".join(header) + "\n")

        for i, line in enumerate(f, start=1):
            if not line.strip():
                continue

            if i % 1000 == 0:
                print(f"  {region}: processed {i} trees", flush=True)

            row = branch_fst_for_line(line, region, tree_file)
            out.write("\t".join(map(str, row)) + "\n")

    tmp_file.replace(out_file)
    print(f"Finished region {region}: {out_file}", flush=True)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate tskit branch-mode Fst for one iteration-2000 ARGweaver .tre.gz file."
    )
    parser.add_argument(
        "--tree-file",
        required=True,
        help="Path to one *.2000.tre.gz file.",
    )
    args = parser.parse_args()

    process_tree_file(args.tree_file)


if __name__ == "__main__":
    main()