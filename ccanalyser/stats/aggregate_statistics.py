import glob
import itertools
import os
import sys
import argparse
import numpy as np
import pandas as pd
import re


class DeduplicationStats:
    def __init__(self, fnames):
        self.fnames = sorted(fnames)

    def _get_sample_name(self, fn):
        sample_name_re = re.compile(r".*\/(.*)\.")
        matches = sample_name_re.match(fn)
        return matches.group(1)

    def _read_file(self, fn, sample_name):
        df = pd.read_csv(fn, sep="\t", header=None, names=["count"], index_col=0)
        return df.transpose().assign(sample=sample_name).set_index("sample")

    @property
    def processed_dataframe(self):
        return pd.concat(
            [self._read_file(fn, self._get_sample_name(fn)) for fn in self.fnames]
        )

    @property
    def total_reads(self):
        return (
            self.processed_dataframe["Read_pairs_processed"]
            .to_frame()
            .assign(read_type="flashed")
            .rename(columns={"Read_pairs_processed": "total_reads"})
            .reset_index()
            .set_index(["sample", "read_type"])
        )

    @property
    def unique_reads(self):
        return (
            self.processed_dataframe["Read_pairs_unique"]
            .reset_index()
            .assign(read_type="flashed")
            .rename(
                columns={
                    "index": "sample",
                    "Read_pairs_unique": "read_pairs_with_unique_sequence",
                }
            )
            .reset_index(drop=True)
            .set_index(["sample", "read_type"])
        )


class DigestionStats:
    def __init__(self, fnames):
        self.fnames = fnames

    def _read_file(self, fn, sample_name):
        return pd.read_csv(fn).assign(sample=sample_name)

    def _get_sample_name(self, fn):
        sample_name_re = re.compile(r".*\/(.*).(flashed|pe).*")
        matches = sample_name_re.match(fn)
        return matches.group(1)

    @property
    def processed_dataframe(self):
        dframes = [self._read_file(fn, self._get_sample_name(fn)) for fn in self.fnames]
        df = pd.concat(dframes)
        return df.groupby(["sample", "stat", "read_type", "bin"]).sum().reset_index()

    @property
    def flashed_reads(self):
        return (
            self.processed_dataframe.loc[
                lambda df: (df["read_type"].isin(["flashed", "r1"]))
                & (df["stat"] == "total")
            ]
            .groupby(["sample", "read_type"])["frequency"]
            .sum()
            .reset_index()
            .assign(read_type=lambda df: df["read_type"].str.replace("r1", "pe"))
            .rename(columns={"frequency": "flashed or unflashed"})
            .set_index(["sample", "read_type"])
        )

    @property
    def digested_reads(self):
        return (
            self.processed_dataframe.loc[
                lambda df: (df["stat"] == "valid")
                & (df["read_type"] != "r2")
                & (df["bin"] != 0)
            ]
            .drop(columns="bin")
            .groupby(["sample", "read_type"])
            .sum()
            .reset_index()[["sample", "read_type", "frequency"]]
            .assign(read_type=lambda df: df["read_type"].str.replace("r1", "pe"))
            .rename(columns={"frequency": "read_pairs_with_restriction_site(s)"})
            .set_index(["sample", "read_type"])
        )


class SliceStats:
    def __init__(self, fnames):
        self.fnames = fnames

    def _read_file(self, fn, sample_name, read_type):
        return (
            pd.read_csv(fn, index_col=0, sep="\t")
            .transpose()
            .reset_index()
            .rename(columns={"index": "filter_stage"})
            .assign(sample=sample_name, read_type=read_type)
        )

    def _get_sample_name(self, fn):
        sample_name_re = re.compile(r".*\/(.*).(flashed|pe).*")
        matches = sample_name_re.match(fn)
        return matches.group(1)

    def _get_read_type(self, fn):
        read_type_re = re.compile(r".*\/(.*).(flashed|pe).*")
        matches = read_type_re.match(fn)
        return matches.group(2)

    @property
    def processed_dataframe(self):
        dframes = [
            self._read_file(fn, self._get_sample_name(fn), self._get_read_type(fn))
            for fn in self.fnames
        ]

        df = pd.concat(dframes)

        # Need to account for nunique samples (need to max these)
        agg_columns = ["sample", "filter_stage", "read_type"]
        agg_dict = {
            stat: np.sum if not "unique" in stat else np.max for stat in df.columns
        }

        for col in agg_columns:
            del agg_dict[col]

        return df.groupby(agg_columns).agg(agg_dict).reset_index()

    @property
    def filtered_reads(self):
        return self.processed_dataframe.pivot(
            index=["sample", "read_type"], columns="filter_stage", values="mapped"
        )


class ReporterStats:
    def __init__(self, fnames):
        self.fnames = fnames

    def _read_file(self, fn, sample_name, read_type):
        return pd.read_csv(
            fn, sep="\t", names=["capture", "cis_or_trans", "count"], header=0
        ).assign(sample=sample_name, read_type=read_type)

    def _get_sample_name(self, fn):
        sample_name_re = re.compile(r".*\/(.*).(flashed|pe).*")
        matches = sample_name_re.match(fn)
        return matches.group(1)

    def _get_read_type(self, fn):
        read_type_re = re.compile(r".*\/(.*).(flashed|pe).*")
        matches = read_type_re.match(fn)
        return matches.group(2)

    @property
    def processed_dataframe(self):
        dframes = [
            self._read_file(fn, self._get_sample_name(fn), self._get_read_type(fn))
            for fn in self.fnames
        ]

        df = pd.concat(dframes)

        return (
            df.groupby(["sample", "capture", "read_type", "cis_or_trans"])
            .sum()
            .reset_index()
        )


def main(
    deduplication_stats,
    digestion_stats,
    ccanalyser_stats,
    reporter_stats,
    output_dir=".",
):

    stats_dict = dict(
        zip(
            [
                "deduplication_stats",
                "digestion_stats",
                "ccanalyser_stats",
                "reporter_stats",
            ],
            [
                DeduplicationStats(deduplication_stats),
                DigestionStats(digestion_stats),
                SliceStats(ccanalyser_stats),
                ReporterStats(reporter_stats),
            ],
        )
    )

    # Output individual aggregated stats files
    out_dir = output_dir.rstrip("/")
    for name, stats in stats_dict.items():
        stats.processed_dataframe.to_csv(f"{out_dir}/{name}.tsv", sep="\t")

    # Combined stats
    readpair_stats = [stats_dict['deduplication_stats'].total_reads,
                      stats_dict['deduplication_stats'].unique_reads,
                      stats_dict['digestion_stats'].flashed_reads,
                      stats_dict['digestion_stats'].digested_reads,
                      stats_dict['ccanalyser_stats'].filtered_reads,]
    
    readpair_stats = (readpair_stats[0]
                                    .join(readpair_stats[1:], how='outer')
                                    .fillna(0)
                                    .reset_index()
                                    .melt(id_vars=['sample', 'read_type'], var_name='stat_type', value_name='read_pairs')
                     )

    readpair_stats.to_csv(f"{out_dir}/combined_stats.tsv", sep="\t")

