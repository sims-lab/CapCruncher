import json

import polars as pl

from capcruncher.api.statistics import FastqTrimmingStatistics


def extract_trimming_stats(trimming_summary_path):
    df_stats = pl.read_csv(trimming_summary_path, separator="\t")
    df_stats_agg = (
        df_stats.with_columns(
            pl.col("Sample").str.split("_").list.last().cast(pl.Int64).alias("read_number"),
            pl.col("Sample").str.extract(r"(.+)_part\d+_\d+", 1).alias("sample"),
        )
        .group_by("sample", "read_number")
        .agg(
            pl.col("r_processed").sum(),
            pl.col("r_written").sum(),
            pl.col("r_with_adapters").sum(),
        )
        .sort("sample", "read_number")
    )

    return [
        FastqTrimmingStatistics(
            sample=row["sample"],
            read_number=row["read_number"],
            reads_input=row["r_processed"],
            reads_output=row["r_written"],
            reads_with_adapter_identified=row["r_with_adapters"],
        )
        for row in df_stats_agg.iter_rows(named=True)
    ]


def write_trimming_stats(trimming_summary_path, output_path):
    stats_json = [
        stat.model_dump_json() for stat in extract_trimming_stats(trimming_summary_path)
    ]
    with open(output_path, "w") as f:
        f.write(json.dumps(stats_json, indent=4))


if "snakemake" in globals():
    write_trimming_stats(
        globals()["snakemake"].input[0],
        globals()["snakemake"].output[0],
    )
