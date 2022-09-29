#!/usr/bin/env python
# coding: utf-8

import argparse
import pickle
from collections import defaultdict
from glob import glob
from pathlib import Path
import pandas as pd

from ehdn_lrdn import Locus, Sample


def main():
    parser = argparse.ArgumentParser(
        description="Merge ehdn_lrdn results into a single TSV file."
    )
    parser.add_argument(
        "dir",
        type=Path,
    )
    parser.add_argument(
        "output",
        type=Path,
    )

    args = parser.parse_args()
    loci_counts_merged = defaultdict(list)

    for path in glob(str(args.dir) + "/*"):
        with open(path, "rb") as f:
            sample = pickle.load(f)
            for locus, value in sample.loci_counts.items():
                entry = f"{sample.name}:{value}"
                loci_counts_merged[locus].append(entry)

    # Build dataframe for final result

    contig = []
    start = []
    end = []
    motif = []
    counts = []

    for locus, sample_counts in loci_counts_merged.items():
        contig.append(locus.contig)
        start.append(locus.start)
        end.append(locus.end)
        motif.append(locus.motif)
        counts.append(",".join(sample_counts))

    df = pd.DataFrame({"contig": contig, "start": start, "end": end, "counts": counts})

    output = str(args.output)
    if not output.endswith(".tsv"):
        output += ".tsv"

    df.to_csv(output, index=False, sep="\t")


if __name__ == "__main__":
    main()
