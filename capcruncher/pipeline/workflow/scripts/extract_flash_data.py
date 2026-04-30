import json

import polars as pl

from capcruncher.api.statistics import FlashStats


def extract_flash_stats(flash_summary_path):
    df_stats = pl.read_csv(flash_summary_path, separator="\t")
    if "combopairs" not in df_stats.columns:
        df_stats = df_stats.with_columns(
            (
                pl.col("totalpairs")
                - pl.col("discardpairs").fill_null(0)
                - pl.col("uncombopairs")
            ).alias("combopairs")
        )

    df_stats = (
        df_stats.with_columns(
            pl.col("Sample").str.split("_part").list.first().alias("sample")
        )
        .group_by("sample")
        .agg(
            pl.col("combopairs").sum(),
            pl.col("uncombopairs").sum(),
        )
        .sort("sample")
    )

    return [
        FlashStats(
            sample=row["sample"],
            n_combined=row["combopairs"],
            n_uncombined=row["uncombopairs"],
        )
        for row in df_stats.iter_rows(named=True)
    ]


def write_flash_stats(flash_summary_path, output_path):
    stats_json = [
        stat.model_dump_json() for stat in extract_flash_stats(flash_summary_path)
    ]
    with open(output_path, "w") as f:
        f.write(json.dumps(stats_json, indent=4))


if "snakemake" in globals():
    write_flash_stats(
        globals()["snakemake"].input[0],
        globals()["snakemake"].output[0],
    )
