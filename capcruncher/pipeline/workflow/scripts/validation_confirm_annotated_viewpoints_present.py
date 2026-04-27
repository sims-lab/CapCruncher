"""
Aim: Ensure that all viewpoints are found in the annotated slices.
"""

import pyranges1 as pr
import polars as pl
import pathlib


def count_annotated_viewpoints(slices):
    with pl.StringCache():
        vp_counts = []
        for pq in slices:
            df = pl.read_parquet(pq, columns=["capture"])
            vp_counts.append(df["capture"].value_counts())

        return (
            pl.concat(vp_counts)
            .group_by("capture")
            .agg(pl.col("count").sum().alias("n_slices"))
            .to_pandas()
        )


def validate_viewpoints_present(slices, viewpoints, output_counts, output_sentinel):
    df_counts = count_annotated_viewpoints(slices)
    df_counts.to_csv(output_counts, sep="\t", index=True)

    gr_viewpoints = pr.read_bed(viewpoints)
    if not gr_viewpoints.Name.isin(df_counts.capture).all():
        missing = gr_viewpoints.Name[~gr_viewpoints.Name.isin(df_counts.capture)]
        raise ValueError(f"Not all viewpoints are present in the annotation: {missing}")

    pathlib.Path(output_sentinel).touch()


if "snakemake" in globals():
    snakemake = globals()["snakemake"]
    validate_viewpoints_present(
        snakemake.input.slices,
        snakemake.input.viewpoints,
        snakemake.output.viewpoints_present,
        snakemake.output.sentinel,
    )
