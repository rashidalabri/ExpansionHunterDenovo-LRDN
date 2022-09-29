#!/usr/bin/env python
# coding: utf-8

import argparse
import json
import pickle
import subprocess
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd


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
            hasattr(other, "contig")
            and hasattr(other, "start")
            and hasattr(other, "end")
            and hasattr(other, "motif")
            and self.contig == other.contig
            and self.start == other.start
            and self.end == other.end
            and self.motif == other.motif
        )

    def __hash__(self) -> int:
        return hash((self.contig, self.start, self.end, self.motif))

    def __str__(self) -> str:
        return f"{self.contig}:{self.start}-{self.end}"

    def __repr__(self) -> str:
        return f"Locus({self.contig}:{self.start}-{self.end};{self.motif})"


class Sample(object):
    def __init__(self, name: str, alignment: Path, profile: Path):
        self.name = name
        self.alignment = alignment
        self.profile = profile
        self.loci_counts = dict()

    @property
    def global_depth(self) -> float:
        with open(self.profile, "r") as f:
            profile = json.load(f)
        return profile["Depth"]

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
    def __init__(self, sample: Sample, padding: int):
        self.sample = sample
        self.padding = padding

    def normalize_sample(self) -> None:
        for locus, count in self.sample.loci_counts.items():
            local_normalized_count = self.normalize_locus_count(locus, count)
            self.sample.loci_counts[locus] = local_normalized_count

    def pad_locus(self, locus: Locus) -> Tuple[int, int]:
        return max(0, locus.start - self.padding), locus.end + self.padding

    def get_locus_local_depth(self, locus: Locus) -> float:
        start_padded, stop_padded = self.pad_locus(locus)
        region = f"{locus.contig}:{start_padded}-{stop_padded}"
        samtools_cmd = ["samtools", "depth", "-r", region, str(self.sample.alignment)]
        awk_cmd = ["awk", "{d+=$3}END{print d}"]

        ps = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE)
        result = subprocess.check_output(awk_cmd, stdin=ps.stdout)
        ps.wait()
        total_read_depth = int(result.strip())

        return total_read_depth / (stop_padded - start_padded)

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
        return count * self.sample.global_depth / local_depth


class Parameters(object):
    def __init__(self, args: argparse.Namespace):
        self.args = args

    @property
    def sample_name(self) -> str:
        return self.args.name

    @property
    def ehdn_calls_df(self) -> pd.DataFrame:
        return pd.read_table(self.args.calls)

    @property
    def alignment(self) -> Path:
        return self.args.alignment

    @property
    def profile(self) -> str:
        return self.args.profile

    @property
    def global_depth(self) -> float:
        with open(self.args.profile, "r") as f:
            profile = json.load(f)
        return profile["Depth"]

    @property
    def padding(self) -> int:
        return self.args.padding

    @property
    def output_dir(self) -> int:
        return self.args.output_dir


def main():
    parser = argparse.ArgumentParser(
        description="Normalize EHDN repeat expansion calls to sample local read depth."
    )
    parser.add_argument(
        "alignment",
        type=Path,
        help="Alignment file to use for normalization (BAM, CRAM, SAM).",
    )
    parser.add_argument(
        "name",
        type=str,
        help="The name of the sample to normalize.",
    )
    parser.add_argument(
        "calls",
        type=Path,
        help="The path to the EHDN repeat expansion calls file.",
    )
    parser.add_argument(
        "profile",
        type=Path,
        help="Path to the STR profile associated with the sample.",
    )
    parser.add_argument(
        "-p",
        "--padding",
        type=int,
        default=1000,
        help="The number of basepairs to pad around each locus.",
    )
    parser.add_argument(
        "output_dir", type=Path, help="The directory to save the output"
    )

    params = Parameters(parser.parse_args())

    sample = Sample(params.sample_name, params.alignment, params.profile)
    sample.build_loci_from_calls(params.ehdn_calls_df)
    normalizer = LocalReadDepthNormalizer(sample, params.padding)
    normalizer.normalize_sample()
    sample.save_to_file(params.output_dir / (params.sample_name + ".pkl"))


if __name__ == "__main__":
    main()
