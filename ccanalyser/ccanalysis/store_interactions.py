import argparse
import os
import re
import sys
from typing import Union

import cooler
import pandas as pd
from ccanalyser.utils.helpers import get_human_readable_number_of_bp
from pandas.core.frame import DataFrame
from pybedtools import BedTool
from natsort import natsort_key
import numpy as np


def get_bins(genome="mm9", bin_size=1000) -> pd.DataFrame:

    bt = BedTool()

    if os.path.isfile(genome):
        windows = bt.makewindows(w=bin_size, g=genome, i="srcwinnum")
    else:
        windows = bt.makewindows(w=bin_size, genome=genome, i="srcwinnum")

    return (
        windows.to_dataframe()
        .reset_index()
        .rename(columns={"index": "bin"})[["chrom", "start", "end", "bin"]]
    )


def get_midpoints(bed: Union[BedTool, pd.DataFrame]) -> pd.DataFrame:

    if isinstance(bed, BedTool):
        bed = bed.to_dataframe()

    return bed.assign(
        midpoint_start=lambda df: df["start"] + ((df["end"] - df["start"]) // 2),
        midpoint_end=lambda df: df["midpoint_start"] + 1,
    )[["chrom", "midpoint_start", "midpoint_end", "name"]]


def format_restriction_fragment(
    restriction_fragment_map: Union[pd.DataFrame, BedTool], bin_method="overlap"
) -> pd.DataFrame:
    if bin_method == "midpoints":
        restriction_fragment = get_midpoints(bed=restriction_fragment_map)
    elif bin_method == "overlap":
        restriction_fragment = restriction_fragment_map
    else:
        raise ValueError(
            'Incorrect bin method provided. Use either "midpoints" or "overlap'
        )

    return restriction_fragment



def bin_restriction_fragment(
    bins: Union[pd.DataFrame, BedTool],
    restriction_fragment_map: Union[pd.DataFrame, BedTool],
    overlap_fraction: int = 0.5,
) -> pd.DataFrame:

    if isinstance(bins, pd.DataFrame):
        bins = BedTool.from_dataframe(bins)

    if isinstance(restriction_fragment_map, pd.DataFrame):
        restriction_fragment_map = BedTool.from_dataframe(restriction_fragment_map)

    # Intersect bins with restriction fragments
    return (
        bins.intersect(restriction_fragment_map, loj=True, F=overlap_fraction)
        .to_dataframe()[["name", "thickEnd"]]
        .rename(columns={"name": "bin", "thickEnd": "restriction_fragment"})
    )


def aggregate_binned_counts(
    counts: pd.DataFrame, binned_restriction_fragments: pd.DataFrame
) -> pd.DataFrame:

    return (
        counts.merge(
            binned_restriction_fragments,
            left_on="rf1",
            right_on="restriction_fragment",
        )
        .merge(
            binned_restriction_fragments,
            left_on="rf2",
            right_on="restriction_fragment",
        )
        .rename(columns={"bin_x": "bin1_id", "bin_y": "bin2_id"})
        .assign(bin1_id_corrected=lambda df: np.where(df['bin1_id'] > df['bin2_id'], df['bin2_id'], df['bin1_id']),
                bin2_id_corrected=lambda df: np.where(df['bin1_id'] == df['bin1_id_corrected'], df['bin2_id'], df['bin1_id'])
        )
        .drop(columns=['bin1_id', 'bin2_id'])
        .rename(columns=lambda col: col.replace('_corrected', ''))
        .sort_values(["bin1_id", "bin2_id"])
        .groupby(["bin1_id", "bin2_id"])["count"]
        .sum()
        .reset_index()
    )


def store_counts_restriction_fragment(
    counts: pd.DataFrame,
    restriction_fragment_map: pd.DataFrame,
    outfile: str = "restriction_fragment_counts.hdf5",
) -> str:

    # Generate a mapping of restriction_fragment name to order in restriction_fragment dataframe (index)
    # e.g. {chr1_DpnII_0: 0}
    df_restriction_fragment_map = restriction_fragment_map.reset_index(drop=True)
    restriction_fragment_mapping = {
        v: k for k, v in restriction_fragment_map["name"].iteritems()
    }

    # Generate "pixels" dataframe. Effectively coordinate based sparse matrix

    #TODO: Look further at bug bin1_id > bin2_id. For now just swaping the bins over as they are equivalent

    pixels = (counts.assign(
        bin1_id=counts["rf1"].map(restriction_fragment_mapping).astype(int),
        bin2_id=counts["rf2"].map(restriction_fragment_mapping).astype(int),
        bin1_id_corrected=lambda df: np.where(df['bin1_id'] > df['bin2_id'], df['bin2_id'], df['bin1_id']),
        bin2_id_corrected=lambda df: np.where(df['bin1_id'] == df['bin1_id_corrected'], df['bin2_id'], df['bin1_id'])
        )
        .drop(columns=['bin1_id', 'bin2_id'])
        .rename(columns=lambda col: col.replace('_corrected', ''))
        .sort_values(["bin1_id", "bin2_id"])
        [["bin1_id", "bin2_id", "count"]])


    cooler.create_cooler(
        outfile,
        bins=df_restriction_fragment_map.drop(columns="name"),
        pixels=pixels[["bin1_id", "bin2_id", "count"]],
        ordered=False,
    )

    return outfile

def scale_binned_counts(counts: pd.DataFrame,
                        scale_factor: int = 1e6, 
                        cis_only=False):

    counts = counts.copy()
    if cis_only:
        raise NotImplementedError('Need to identify which bins belong to each chrom')
    
    counts['count_scaled'] = (counts['counts'] / scale_factor) * counts['counts'].sum()

    return counts


def main(
    restriction_fragment_counts,
    restriction_fragment_map,
    genome,
    binsizes=None,
    output_prefix="",
    bin_method="overlap",
    overlap_fraction=0.5,
    store_restriction_fragment_counts=True,

):
    
    # Format output prefix
    output_prefix = output_prefix.replace(".hdf5", "")
    
    # Load counts
    df_restriction_fragment_counts = pd.read_csv(
        restriction_fragment_counts, sep="\t"
    ).sort_values(["rf1", "rf2"], key=natsort_key)

    # Load restriction_fragment map
    df_restriction_fragment_map = pd.read_csv(
        restriction_fragment_map,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name"],
    )
    
    # Store interactions at the restriction_fragment level
    if store_restriction_fragment_counts:
        store_counts_restriction_fragment(
            counts=df_restriction_fragment_counts,
            restriction_fragment_map=df_restriction_fragment_map,
            outfile=f"{output_prefix}_restriction_fragments.hdf5",
        )

    # Bin counts and store
    for bin_size in binsizes:

        bins = get_bins(genome=genome, bin_size=bin_size)
        restriction_fragment = format_restriction_fragment(
            restriction_fragment_map=df_restriction_fragment_map, bin_method=bin_method
        )

        restriction_fragment_binned = bin_restriction_fragment(
            bins=bins,
            restriction_fragment_map=restriction_fragment,
            overlap_fraction=overlap_fraction if not bin_method == "midpoint" else 1,
        )
        binned_counts = aggregate_binned_counts(
            counts=df_restriction_fragment_counts,
            binned_restriction_fragments=restriction_fragment_binned,
        )

        # Scale binned counts
        binned_counts = scale_binned_counts(binned_counts)

        # Store binned counts
        cooler.create_cooler(
            f"{output_prefix}_{get_human_readable_number_of_bp(bin_size)}_binned.hdf5",
            bins=bins[["chrom", "start", "end"]],
            pixels=binned_counts[["bin1_id", "bin2_id", "count", "count_scaled"]],
            columns=['count_scaled'],
            dtypes=[np.float32],
            ordered=True,
        )



if __name__ == "__main__":
    pass

    # parser = argparse.ArgumentParser()
    # parser.add_argument("-c", "--restriction_fragment_counts", required=True)
    # parser.add_argument("-m", "--restriction_fragment_map", required=True)
    # parser.add_argument("-g", "--genome", default="mm9", required=True, type=str)
    # parser.add_argument("-b", "--binsizes", default=1000, nargs='+', type=int)
    # parser.add_argument("--bin_method", default='overlap')
    # parser.add_argument("-f", '--overlap_fraction', default=0.5, type=float)
    # parser.add_argument("-o", "--output_prefix", default="counts.hdf5")
    # args = parser.parse_args()

    # main(**vars(args))
