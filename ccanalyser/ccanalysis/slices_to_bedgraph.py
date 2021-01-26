#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Converts the tsv output of ccanalyser.py to a bedgraph by intersecting
the reporter slices with a bed file. The pipeline uses the generated restriction
map of the genome  (ccanalysis/restriction_enzyme_map/genome.digest.bed.gz).

"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
from pybedtools import BedTool
import tempfile


class CCBedgraph(object):
    def __init__(
        self,
        path=None,
        df=None,
        capture_name="",
        capture_chrom="",
        capture_start="",
        capture_end="",
    ):

        self.fn = path
        self.df = df

        if self.fn:
            self.df = pd.read_csv(
                self.fn, sep="\t", header=None, names=["chrom", "start", "end", "score"]
            )

        self.capture_name = capture_name
        self.capture_chrom = capture_chrom
        self.capture_start = capture_start
        self.capture_end = capture_end

    @property
    def score(self):
        return self.df.rename(columns={"score": self.fn})[self.fn]

    @property
    def coordinates(self):
        return self.df.loc[:, "chrom":"end"]

    def to_bedtool(self):
        return self.df.pipe(BedTool.from_dataframe)

    def normalise_by_n_cis(self, reporter_distance=1e5, n_read_scale=1e6):

        n_reporters = (
            self.df.loc[lambda df: (df["chrom"] == self.capture_chrom)]
            .loc[lambda df: df["start"] >= (self.capture_start - reporter_distance)]
            .loc[lambda df: df["end"] <= (self.capture_end + reporter_distance)][
                "score"
            ]
            .sum()
        )

        return (self / (n_reporters / n_read_scale))

    def to_file(self, path):
        self.df.to_csv(path, sep="\t", header=None, index=None)

    def __add__(self, other):
        if isinstance(other, CCBedgraph):
            self.df["score"] = self.df["score"] + other.df["score"]
            return self

        elif isinstance(other, (np.ndarray, pd.Series, int, float)):
            self.df["score"] = self.df["score"] + other
            return self

        else:
            return NotImplementedError()

    def __sub__(self, other):
        if isinstance(other, CCBedgraph):
            self.df["score"] = self.df["score"] - other.df["score"]
            return self

        elif isinstance(other, (np.ndarray, pd.Series, int, float)):
            self.df["score"] = self.df["score"] - other
            return self

        else:
            return NotImplementedError()

    def __mul__(self, other):
        if isinstance(other, CCBedgraph):
            self.df["score"] = self.df["score"] * other.df["score"]
            return self

        elif isinstance(other, (np.ndarray, pd.Series, int, float)):
            self.df["score"] = self.df["score"] * other
            return self

        else:
            return NotImplementedError()

    def __truediv__(self, other):
        if isinstance(other, CCBedgraph):
            self.df["score"] = self.df["score"] / other.df["score"]
            return self

        elif isinstance(other, (np.ndarray, pd.Series, int, float)):
            self.df["score"] = self.df["score"] / other
            return self

        else:
            return NotImplementedError()


class CCBedgraphCollection(object):
    def __init__(self, bedgraphs):
        self.bdg_fnames = bedgraphs
        self.bedgraphs = [CCBedgraph(bg) for bg in self.bdg_fnames]

    def normalise_bedgraphs(self, how="region", **kwargs):

        if how == "region":
            self.bedgraphs = [bg.normalise_by_region(**kwargs) for bg in self.bedgraphs]
        else:
            raise NotImplementedError("No other methods")
        return self

    def get_average_bedgraph(self):
        coords = self.bedgraphs[0].coordinates
        scores = [bg.score for bg in self.bedgraphs]
        df_scores = pd.concat(scores, axis=1)
        df_scores_mean = df_scores.mean(axis=1)
        bdg = (
            pd.concat([coords, df_scores_mean], axis=1).pipe(BedTool.from_dataframe).fn
        )
        return CCBedgraph(bdg)


def make_bedgraph(reporters, bed):

    df_reporters = pd.read_csv(reporters, sep="\t").query('capture == "."')

    df_reporters[["start", "end"]] = df_reporters[["start", "end"]].astype(int)

    bt_reporters = BedTool.from_dataframe(
        df_reporters[["chrom", "start", "end", "slice_name"]]
    )

    bt_bed = BedTool(bed)
    bedgraph = (
        bt_bed.intersect(bt_reporters, c=True)
        .to_dataframe()
        .sort_values(["chrom", "start"])
        .iloc[:, [0, 1, 2, 4]]
    )

    return bedgraph


def main(
    slices,
    bed,
    output="out.bedgraph",
    normalise=None,
    normalise_reporter_distance=1e5,
    normalise_scale=1e6,
):

    """
    Converts reporter tsv to bedgraph

    Args:
        slices - slices tsv file name
        bed - bed file to intersect with reporters
        output - file name for output file
    """

    df_rf_map = pd.read_csv(
        bed, sep="\t", header=None, names=["chrom", "start", "end", "name"]
    )

    for ii, df_slices in enumerate(pd.read_csv(slices, sep="\t", chunksize=2e6)):

        print("=" * 4)
        print(f"Iteration: {ii}")
        print("=" * 4)

        df_reporters = df_slices.query("capture_count == 0")
        df_captures = df_slices.query("capture_count > 0")

        # Slight hack to enable bedgraphs to be generated for multicapture fragments
        df_primary_captures = df_captures.groupby("parent_read").first()
        df_not_primary_capture = df_captures.loc[
            ~(df_captures["slice_name"].isin(df_primary_captures["slice_name"]))
        ]
        df_reporters = pd.concat([df_reporters, df_not_primary_capture])

        # Count restriction fragments and merge with rf bed file
        bedgraph = (
            df_reporters["restriction_fragment"]
            .value_counts()
            .reset_index()
            .rename(columns={"index": "rf", "restriction_fragment": "score"})
            .merge(df_rf_map, left_on="rf", right_on="name", how="outer")[
                ["chrom", "start", "end", "score"]
            ]
            .sort_values(["chrom", "start"])
            .fillna(0)
        )

        # Get capture site coords (will obtain the best guess of coords from the captures dataframe)
        capture_site = {
            "chrom": df_captures["chrom"].mode()[0],
            "start": df_captures["start"].min(),
            "end": df_captures["end"].max(),
            "capture": df_captures["capture"].mode()[0],
        }

        # Create CCBedgraph
        ccbdg = CCBedgraph(
            df=bedgraph,
            capture_chrom=capture_site["chrom"],
            capture_start=capture_site["start"],
            capture_end=capture_site["end"],
            capture_name=capture_site["capture"],
        )

        print("Merging bedgraphs")
        if ii == 0:
            bedgraph = ccbdg
        else:
            bedgraph = bedgraph + ccbdg

    # Run normalisation if required
    if normalise == "n_cis":
        bedgraph = bedgraph.normalise_by_n_cis(
            reporter_distance=normalise_reporter_distance, n_read_scale=normalise_scale
        )

    print("Outputting to file")
    bedgraph.to_file(output)

