import os
import polars as pl


def count_identified_viewpoints(parquet_path):
    if os.path.isdir(parquet_path):
        parquet_path = os.path.join(parquet_path, "*.parquet")

    return (
        pl.scan_parquet(parquet_path)
        .select(["viewpoint", "pe"])
        .unique()
        .filter((pl.col("viewpoint") != "") & pl.col("viewpoint").is_not_null())
        .filter((pl.col("pe") != "") & pl.col("pe").is_not_null())
        .collect()
    )


def write_identified_viewpoints(parquet_path, output_path):
    count_identified_viewpoints(parquet_path).write_csv(output_path, separator="\t")


if "snakemake" in globals():
    write_identified_viewpoints(snakemake.params.slices_dir, snakemake.output[0])
