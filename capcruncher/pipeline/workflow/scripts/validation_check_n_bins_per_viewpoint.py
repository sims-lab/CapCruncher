"""
Aim: Check that there is only one restriction fragment per viewpoint.
"""

import pandas as pd
import numpy as np
import pyranges as pr
import polars as pl
import tabulate
import pathlib
from loguru import logger


with logger.catch():
    logger.info("Checking that there is only one restriction fragment per viewpoint")

    bins = snakemake.input.bins
    viewpoints = snakemake.input.viewpoints

    df_bins = pl.read_csv(
        bins,
        separator="\t",
        has_header=False,
        new_columns=["Chromosome", "Start", "End", "Name"],
    )
    gr_bins = pr.PyRanges(df_bins.to_pandas())

    df_viewpoints = pl.read_csv(
        viewpoints,
        separator="\t",
        has_header=False,
        new_columns=["Chromosome", "Start", "End", "Name"],
    )
    gr_viewpoints = pr.PyRanges(df_viewpoints.to_pandas())

    # Generate a table with the number of restriction fragments overlapped by each viewpoint
    gr_bin_vp_overlap = gr_bins.join(gr_viewpoints, suffix="_viewpoints")
    bin_vp_counts = gr_bin_vp_overlap.df["Name_viewpoints"].value_counts()
    has_multiple_bins = bin_vp_counts > 1

    df_viewpoints = df_viewpoints.to_pandas()
    df_viewpoints["n_restriction_fragments_overlapped"] = 1
    df_viewpoints = df_viewpoints.set_index("Name")
    df_viewpoints.loc[
        has_multiple_bins, "n_restriction_fragments_overlapped"
    ] = bin_vp_counts[has_multiple_bins].values
    df_viewpoints = df_viewpoints.reset_index().rename(columns={"index": "Viewpoint"})
    df_viewpoints.to_csv(snakemake.output.n_bins_per_viewpoint, sep="\t", index=False)

    # df_rf_counts = (
    #     df_viewpoints.to_pandas()
    #     .set_index("Name")
    #     .loc[has_multiple_bins]
    #     .assign(
    #         n_restriction_fragments_overlapped=bin_vp_counts[has_multiple_bins].values
    #     )
    # )

    # df_rf_counts = pd.concat(
    #     [
    #         df_rf_counts,
    #         df_viewpoints.to_pandas()
    #         .set_index("Name")
    #         .loc[~has_multiple_bins]
    #         .assign(n_restriction_fragments_overlapped=1),
    #     ]
    # )

    # df_rf_counts = (
    #     df_rf_counts.reset_index()
    #     .rename(columns={"index": "Viewpoint"})
    #     .to_csv(snakemake.output.n_bins_per_viewpoint, sep="\t", index=False)
    # )

    if (
        has_multiple_bins.any()
        and not snakemake.params.ignore_multiple_bins_per_viewpoint
    ):
        tbl = tabulate.tabulate(df_rf_counts, headers="keys", tablefmt="psql")

        raise ValueError(
            f"""The following viewpoints overlap multiple restriction fragments:\n{df_rf_counts}\n"""
        )

    else:
        pathlib.Path(snakemake.output.sentinel).touch()
