import os
import polars as pl

parquet_path = snakemake.params.slices_dir
if os.path.isdir(parquet_path):
    parquet_path = os.path.join(parquet_path, "*.parquet")

unique_viewpoints = (
    pl.scan_parquet(parquet_path)
    .select(["viewpoint", "pe"])
    .unique()
    .filter((pl.col("viewpoint") != "") & pl.col("viewpoint").is_not_null())
    .filter((pl.col("pe") != "") & pl.col("pe").is_not_null())
    .collect()
)

unique_viewpoints.write_csv(snakemake.output[0], separator="\t")
