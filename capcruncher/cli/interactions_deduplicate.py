import pandas as pd
import os
import ibis
from ibis import _
import pyarrow.dataset as ds
import shutil
from loguru import logger

ibis.options.interactive = False


def deduplicate(
    slices: os.PathLike,
    output: os.PathLike,
    read_type: str = "flashed",
    sample_name: str = "",
    stats_prefix: os.PathLike = "",
):
    logger.info("Connecting to DuckDB")
    con = ibis.duckdb.connect()

    if not os.path.isdir(slices):
        slices_tbl_raw = con.register(f"parquet://{slices}", table_name="slices_tbl")
    else:
        slices_tbl_raw = con.register(
            f"parquet://{slices}/*.parquet", table_name="slices_tbl"
        )

    if read_type == "pe":
        logger.info("Read type is PE")
        logger.info("Identifying unique fragment IDs")
        query = (
            slices_tbl_raw[["chrom", "start", "end", "parent_id"]]
            .order_by(["chrom", "start", "end", "parent_id"])
            .group_by(by="parent_id", order_by=["chrom", "start", "end"])
            .order_by(["chrom", "start", "end", "parent_id"])
            .mutate(
                slice_f_chrom=_.chrom.first(),
                slice_f_start=_.start.first(),
                slice_l_end=_.end.last(),
            )
            .group_by(["slice_f_chrom", "slice_f_start", "slice_l_end"])
            .mutate(pid=_.parent_id.first())[["pid"]]
            .distinct()["pid"]
        )
    elif read_type == "flashed":
        logger.info("Read type is Flashed")
        logger.info("Identifying unique fragment IDs")

        query = (
            slices_tbl_raw[["coordinates", "parent_id"]]
            .group_by(by="parent_id", order_by=["coordinates"])
            .aggregate(coordinates=lambda t: t.coordinates.group_concat(","))
            .group_by("coordinates")
            .mutate(parent_id_unique=_.parent_id.first())[["parent_id_unique"]]
            .distinct()["parent_id_unique"]
        )

    parent_ids_unique = query.execute(limit=None)

    logger.info("Writing deduplicated slices to disk")
    slices_unfiltered_ds = ds.dataset(slices, format="parquet")
    scanner = slices_unfiltered_ds.scanner(
        filter=ds.field("parent_id").isin(parent_ids_unique)
    )

    if os.path.exists(output):
        shutil.rmtree(output)

    ds.write_dataset(
        scanner,
        output,
        format="parquet",
        partitioning_flavor="hive",
        min_rows_per_group=0,
    )

    # If the output directory is empty, create a dummy file to prevent downstream errors
    if not os.path.exists(output):
        os.makedirs(output)
        df_dummy = scanner.to_table().to_pandas()
        df_dummy.to_parquet(f"{output}/dummy.parquet")

    logger.info("Calculating deduplication stats")
    # Calculate the number of slices in the input

    n_reads_total = (
        slices_tbl_raw.group_by("parent_id")
        .agg([_.count().name("count")])["count"]
        .sum()
        .execute(limit=None)
    )

    # Calculate the number of slices in the output
    n_reads_unique = parent_ids_unique.shape[0]

    # Prepare stats
    df_stats = pd.DataFrame()
    df_stats["stat_type"] = ["not-deduplicated", "deduplicated"]
    df_stats["stat"] = [n_reads_total, n_reads_unique]
    df_stats["sample"] = sample_name
    df_stats["read_type"] = read_type
    df_stats["read_number"] = 0
    df_stats["stage"] = "deduplicate_slices"
    df_stats["stat"] = df_stats["stat"].fillna(0)

    logger.info("Deduplication stats:")
    logger.info(f"\n{df_stats.to_string()}")

    logger.info("Writing deduplication stats to disk")
    df_stats.to_csv(f"{stats_prefix}.read.stats.csv", index=False)

    return df_stats
