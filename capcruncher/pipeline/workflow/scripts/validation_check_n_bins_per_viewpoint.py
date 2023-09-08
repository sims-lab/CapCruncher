"""
Aim: Check that there is only one restriction fragment per viewpoint.
"""

import pandas as pd
import numpy as np
import pyranges as pr
import polars as pl
import tabulate

bins = snakemake.input.bins
viewpoints = snakemake.input.viewpoints


df_bins = pl.read_csv(bins, sep="\t", columns=["Chromosome", "Start", "End", "Name"])
gr_bins = pr.PyRanges(df_bins.to_pandas())

gr_viewpoints = pr.read_bed(viewpoints)

# Generate a table with the number of restriction fragments overlapped by each viewpoint
gr_bin_vp_overlap = gr_bins.join(gr_viewpoints)
bin_vp_counts = gr_bin_vp_overlap.Name.value_counts()
has_multiple_bins = bin_vp_counts > 1
df_rf_counts = (
    gr_viewpoints.df.set_index("Name")
    .loc[has_multiple_bins]
    .assign(n_restriction_fragments_overlapped=bin_vp_counts[has_multiple_bins])
)
df_rf_counts = pd.concat(
    [
        df_rf_counts,
        gr_viewpoints.df.set_index("Name")
        .loc[~has_multiple_bins]
        .assign(n_restriction_fragments_overlapped=1),
    ]
)


df_rf_counts = (
    df_rf_counts.reset_index()
    .rename(columns={"index": "Viewpoint"})
    .to_csv(snakemake.output.n_bins_per_viewpoint, sep="\t", index=False)
)

if has_multiple_bins.any() and not snakemake.params.ignore_multiple_bins_per_viewpoint:
    tbl = tabulate.tabulate(df_rf_counts, headers="keys", tablefmt="psql")

    raise ValueError(
        f"""The following viewpoints overlap multiple restriction fragments:\n{df_rf_counts}\n"""
    )
