import os
import re
import sys
import argparse
from collections import Counter
from typing import Union

import cooler
from h5py._hl.selections2 import ScalarReadSelection
from iced.filter import filter_high_counts, filter_low_counts
import pandas as pd
from pybedtools import BedTool
import numpy as np
from scipy import sparse
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import iced
from seaborn.matrix import _matrix_mask

from ccanalyser.utils.helpers import (
    bed_has_duplicate_names,
    bed_has_name,
    get_human_readable_number_of_bp,
)
from ccanalyser.ccanalysis.store_interactions import (
    format_restriction_fragment,
    bin_restriction_fragment,
)


def format_coordinates(coordinates):

    """Returns a bedtool object with coordinates"""

    coordinates = str(coordinates)
    pattern_genomic_coord = re.compile(r"chr[0-2xXyYmM][0-9]*:\d+-\d+(\s\w)*$")
    pattern_bed_file = re.compile(r"(.*).bed")

    # breakpoint()
    if pattern_genomic_coord.match(coordinates):

        coordinates_split = re.split(":|-", coordinates)
        if len(coordinates_split) < 4:
            coordinates_split.append("region_0")

        bt = BedTool(" ".join(coordinates_split), from_string=True)

    elif pattern_bed_file.match(coordinates) and bed_has_name(coordinates):
        bt = BedTool(coordinates)

    elif pattern_bed_file.match(coordinates) and not bed_has_name(coordinates):
        bt = (
            BedTool(coordinates)
            .to_dataframe()
            .reset_index()
            .assign(name=lambda df: "region_" + df["index"].astype("string"))[
                ["chrom", "start", "end", "name"]
            ]
            .pipe(BedTool.from_dataframe)
        )
    else:
        raise ValueError(
            """Coordinates not provided in the correct format. Provide coordinates in the form chr[NUMBER]:[START]-[END] or a .bed file"""
        )

    return bt


def convert_interval_to_coords(interval):
    return f'{interval["chrom"]}:{interval["start"]}-{interval["end"]}'


def plot_matrix(matrix, figsize=(10, 10), axis_labels=None, cmap=None, thresh=1):

    if not thresh:
        thresh = np.percentile(matrix, 98)

    fig, ax = plt.subplots()
    ax.imshow(matrix, vmax=thresh, vmin=0.0001, cmap=cmap)
    ax.axis("off")

    return fig


def get_region_raw(
    cooler_binned: cooler.Cooler,
    cooler_restriction_fragments: cooler.Cooler,
    region: str,
    **normalisation_kwargs,
):

    return cooler_binned.matrix(balance=False).fetch(region).astype(float)


def normalise_region_ice(
    cooler_binned: cooler.Cooler,
    cooler_restriction_fragments: cooler.Cooler,
    region: str,
    filter_low_counts: float = 0,
    filter_high_counts: float = 0,
    **normalisation_kwargs,
):

    if normalisation_kwargs['scale_counts']:
        matrix = cooler_binned.matrix(balance=False, field='count_scaled').fetch(region).astype(float)
    else:
        matrix = cooler_binned.matrix(balance=False).fetch(region).astype(float)


    if filter_low_counts:
        matrix = iced.filter.filter_low_counts(matrix, percentage=filter_low_counts)

    if filter_high_counts:
        matrix = iced.filter.filter_high_counts(matrix)

    # Need to remove unwanted keyword args
    del normalisation_kwargs["binning_method"]
    del normalisation_kwargs["overlap_fraction"]
    del normalisation_kwargs['scale_counts']

    matrix_normalised = iced.normalization.ICE_normalization(matrix, **normalisation_kwargs)

    return matrix_normalised

def extract_bins_from_cooler(cooler_obj, region=None):

    if region:
        return (
            cooler_obj.bins()
            .fetch(region)
            .reset_index()
            .rename(columns={"index": "bin"})[["chrom", "start", "end", "bin"]]
        )
    else:
        return cooler_obj.bins()[:]


def normalise_region_scaled(
    cooler_binned: cooler.Cooler,
    cooler_restriction_fragments: cooler.Cooler,
    region: str,
    binning_method="overlap",
    overlap_fraction=0.5,
    **normalisation_kwargs,
):
    def normalise_tric_counts(df, scaling_factor=1000000):
        n_restriction_fragments = df["bin1_n_rf"] * df["bin2_n_rf"]
        total = df["count"].sum()

        return (df["count"] / n_restriction_fragments) * (scaling_factor / total)

    bins_binned = extract_bins_from_cooler(cooler_binned, region=region)
    bins_restriction_fragments = extract_bins_from_cooler(
        cooler_restriction_fragments, region=region
    )
    bins_restriction_fragments = format_restriction_fragment(
        bin_restriction_fragment, bin_method=binning_method
    )

    binned_restriction_fragments = bin_restriction_fragment(
        bins_binned, bins_restriction_fragments, overlap_fraction=overlap_fraction
    )

    n_rf_per_bin = (
        binned_restriction_fragments.groupby("bin")
        .size()
        .reset_index()
        .rename(columns={0: "rf_per_bin"})["rf_per_bin"]
    )

    # Determine the number of rf per bin  for the normalisation
    df_sparse_counts = cooler_binned.pixels().fetch(region)
    df_sparse_counts = (
        df_sparse_counts.merge(n_rf_per_bin, left_on="bin1_id", right_index=True)
        .rename(columns={"rf_per_bin": "bin1_n_rf"})
        .merge(n_rf_per_bin, left_on="bin1_id", right_index=True)
        .rename(columns={"rf_per_bin": "bin2_n_rf"})
    )

    # Normalise count data
    df_sparse_counts["normalised_count"] = normalise_tric_counts(df_sparse_counts)

    # Generate a new sparse matrix
    offset = cooler_binned.offset(region)
    x, y = (
        (df_sparse_counts["bin1_id"] - offset),
        (df_sparse_counts["bin2_id"] - offset),
    )
    data = df_sparse_counts["normalised_count"]
    sm = sparse.coo_matrix((data, (x, y)))

    return sm.todense()


def get_normalisation(normalisation: str, method: str):
    method_default_dict = {
        "tri": normalise_region_scaled,
        "tiled": normalise_region_ice,
    }
    normalisations_dict = {
        "infer": method_default_dict[method],
        "raw": get_region_raw,
        "ice": normalise_region_ice,
        "scaling_factor": normalise_region_scaled,
    }

    return normalisations_dict[normalisation]


def main(
    method: str,
    coordinates: str,
    cooler_binned: str,
    cooler_restriction_fragments: str,
    scale_counts: bool = True,
    normalisation=None,
    output_prefix: str = None,
    output_format: str = "png",
    scaling_factor=1e5,
    binning_method="overlap",
    overlap_fraction=0.5,
    filter_low_counts=0.02,
    filter_high_counts=0.02,
    cmap="viridis",
    thresh=0,
    output_matrix=None,
):

    # Count data
    cooler_binned = cooler.Cooler(cooler_binned)
    cooler_restriction_fragments = cooler.Cooler(cooler_restriction_fragments)
    resolution = get_human_readable_number_of_bp(cooler_binned.binsize)



    # Get normalisation function to use
    normalisation_func = get_normalisation(normalisation=normalisation, method=method)
    normalisation_kwargs = {'filter_low_counts': filter_low_counts,
                            'filter_high_counts': filter_high_counts,
                            'binning_method': binning_method,
                            'overlap_fraction': overlap_fraction,
                            'scale_counts': scale_counts}


    for interval in format_coordinates(coordinates):

        interval_coords = convert_interval_to_coords(interval)
        figure_path = f'{output_prefix}_{resolution}_{interval["name"]}.{output_format}'

        matrix_normalised = normalisation_func(
            cooler_binned,
            cooler_restriction_fragments,
            region=interval_coords,
            **normalisation_kwargs
        )


        fig = plot_matrix(matrix=matrix_normalised, cmap=cmap, thresh=thresh)
        fig.savefig(figure_path)

        if output_matrix:
            np.savetxt(output_matrix)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "cooler_restriction_fragments",
        help="HDF5 file in Cooler format containing counts of restriction fragment combination counts",
    )
    parser.add_argument(
        "cooler_binned",
        help="HDF5 file in Cooler format containing binned counts of restriction fragment combination counts",
    )
    parser.add_argument("-c", "--coordinates", required=True, type=str)
    parser.add_argument("-m", "--method", choices=["tri", "tiled"])
    parser.add_argument(
        "-n",
        "--normalisation",
        default=None,
        choices=["infer", "raw", "ice", "scaling_factor"],
    )
    parser.add_argument("--scaling_factor", default=1e5, type=int)
    parser.add_argument("--binning_method", default="overlap")
    parser.add_argument("--overlap_fraction", default=0.5, type=float)
    parser.add_argument("--filter_low_counts", default=0.02, type=float)
    parser.add_argument("--filter_high_counts", default=0, type=float)
    parser.add_argument("-o", "--output_prefix", default="")
    parser.add_argument(
        "-f", "--output_format", default="png", choices=["png", "svg", "jpeg"]
    )
    parser.add_argument("--cmap", help="Colour map to use", default="viridis")
    parser.add_argument('--thresh', type=int, default=0)
    args = parser.parse_args()

    main(**vars(args))
