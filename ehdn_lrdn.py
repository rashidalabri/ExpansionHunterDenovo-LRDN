#!/usr/bin/env python
# coding: utf-8
import argparse
import pandas as pd
from pathlib import Path
import pysam
import pysamstats
import numpy as np
import json
from typing import Tuple
import pickle


class Locus(object):
    def __init__(self, contig: str, start: int, end: int, motif: str):
        self.contig = contig
        self.start = start
        self.end = end
        self.motif = motif

    @property
    def length(self):
        return self.end - self.start

    def __eq__(self, other: object) -> bool:
        return (
            hasattr(other, "chrom")
            and hasattr(other, "start")
            and hasattr(other, "end")
            and hasattr(other, "motif")
            and self.contig == other.chrom
            and self.start == other.start
            and self.end == other.end
            and self.motif == other.motif
        )

    def __hash__(self) -> int:
        return hash((self.contig, self.start, self.end, self.motif))

    def __str__(self) -> str:
        return f"Locus<{self.contig}:{self.start}-{self.end}, {self.motif}>"

    def __repr__(self) -> str:
        return self.__str__()


class Sample(object):
    def __init__(self, name: str):
        self.name = name
        self.loci_counts = dict()

    def build_loci_from_calls(self, calls: pd.DataFrame) -> None:
        for _, row in calls.iterrows():
            locus = Locus(row["contig"], row["start"], row["end"], row["motif"])
            count = self.find_count_in_row(row)
            self.loci_counts[locus] = count

    def find_count_in_row(self, row: pd.Series) -> int:
        counts = [entry.split(":") for entry in row["counts"].split(",")]
        for sample, count in counts:
            if sample == self.name:
                return float(count)
        return 0

    def save_to_file(self, path: Path) -> None:
        with open(path, "wb") as f:
            pickle.dump(self, f)


class LocalReadDepthNormalizer(object):
    def __init__(
        self,
        sample: Sample,
        reads: pysam.AlignmentFile,
        global_depth: float,
        padding: int,
    ):
        self.sample = sample
        self.padding = padding
        self.reads = reads
        self.global_depth = global_depth

    def normalize_sample(self) -> None:
        for locus, count in self.sample.loci_counts.items():
            local_normalized_count = self.normalize_locus_count(locus, count)
            self.sample.loci_counts[locus] = local_normalized_count

    def pad_locus(self, locus: Locus) -> Tuple[int, int]:
        return locus.start - self.padding, locus.end + self.padding

    def get_locus_local_depth(self, locus: Locus) -> float:
        padded_start, padded_end = self.pad_locus(locus)
        local_coverage = pysamstats.load_coverage(
            self.reads, chrom=locus.contig, start=padded_start, end=padded_end
        )
        return np.mean(local_coverage.reads_all)

    def normalize_locus_count(self, locus: Locus, count: float) -> float:
        """Normalizes a locus' Anchored IRR count by the local read depth.

        The Anchored IRR count supplied by ExpansionHunter denovo is normalized
        by the global depth by default. Global normalization is achieved by
        dividing the raw count by the global depth found in the str_profile.json file.

        To retrieve the raw count, we multiply the global normalized count by
        the global depth. Then, we divide the raw count by the local depth to
        find the locally normalized count.

        We define the local depth as the average number of reads overlapping
        each basepair in the locus in addition to a buffer window on either side.

        Args:
            locus (Locus): The locus to normalize.
            count (float): The globally normalized Anchored IRR count.

        Returns:
            float: The locally normalized Anchored IRR count.
        """
        local_depth = self.get_locus_local_depth(locus)
        return count * self.global_depth / local_depth


class Parameters(object):
    def __init__(self, args: argparse.Namespace):
        self.args = args

    @property
    def sample_name(self) -> str:
        return self.args.sample_name

    @property
    def ehdn_calls_df(self) -> pd.DataFrame:
        return pd.read_table(self.args.calls)

    @property
    def alignment_file(self) -> pysam.AlignmentFile:
        with open(self.args.reads, "rb") as f:
            return pysam.AlignmentFile(f)

    @property
    def global_depth(self) -> float:
        with open(self.args.profile, "r") as f:
            profile = json.load(f)
        return profile["Depth"]

    @property
    def padding(self) -> int:
        return self.args.padding


def main():
    parser = argparse.ArgumentParser(
        description="Normalize EHDN calls to local read depth."
    )
    parser.add_argument(
        "-r",
        "--reads",
        required=True,
        type=Path,
        help="Alignment file to use for normalization (BAM, CRAM, SAM).",
    )
    parser.add_argument(
        "-s",
        "--sample_name",
        type=str,
        required=True,
        help="The name of the sample to normalize.",
    )
    parser.add_argument(
        "-c",
        "--calls",
        type=Path,
        required=True,
        help="The path to the EHDN calls file.",
    )
    parser.add_argument(
        "-p",
        "--profile",
        type=Path,
        required=True,
        help="Path to the STR profile associated with the sample.",
    )
    parser.add_argument(
        "-w",
        "--padding",
        type=int,
        default=1000,
        help="The number of basepairs to pad around each locus.",
    )

    args = parser.parse_args()
    parameters = Parameters(args)

    sample = Sample(parameters.sample_name)
    sample.build_loci_from_calls(parameters.ehdn_calls_df)
    normalizer = LocalReadDepthNormalizer(
        sample, parameters.alignment_file, parameters.global_depth, parameters.padding
    )
    normalizer.normalize_sample()
    sample.save_to_file(parameters.sample_name + ".pkl")


if __name__ == "__main__":
    main()
