"""
Aim: Ensure that all viewpoints are found in the annotated slices.
"""

import pandas as pd
import numpy as np
import pyranges as pr
import polars as pl
import tabulate
import pathlib

slices = snakemake.input.slices
viewpoints = snakemake.input.viewpoints

ignore_missing_viewpoints = snakemake.params.ignore_missing_viewpoints


gr_viewpoints = pr.read_bed(viewpoints)

with pl.StringCache():
    vp_counts = []
    for pq in slices:
        df = pl.read_parquet(pq, columns=["capture"])
        vp_counts.append(df["capture"].value_counts())

    df_counts = (
        pl.concat(vp_counts).groupby("capture").agg(pl.sum("counts")).to_pandas()
    )


df_counts.rename(columns={"counts": "n_slices"}).to_csv(
    snakemake.output.viewpoints_present, sep="\t", index=True
)


if not gr_viewpoints.df.Name.isin(df_counts.capture).all():
    # check which viewpoints are missing
    missing = gr_viewpoints.df.Name[~viewpoints.df.Name.isin(df_counts.capture)]
    raise ValueError(f"Not all viewpoints are present in the annotation: {missing}")

else:
    pathlib.Path(snakemake.output.sentinel).touch()
