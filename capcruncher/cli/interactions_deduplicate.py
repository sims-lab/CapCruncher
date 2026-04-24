import os
import ibis
from ibis import _
import pyarrow.dataset as ds
import shutil
from loguru import logger

from capcruncher.api.statistics import AlignmentDeduplicationStats

ibis.options.interactive = False


def read_parquet_with_ibis(con, path: os.PathLike):
    parquet_path = str(path)
    if os.path.isdir(parquet_path):
        parquet_path = os.path.join(parquet_path, "*.parquet")

    return con.read_parquet(parquet_path)


def deduplicate(
    slices: os.PathLike,
    output: os.PathLike,
    read_type: str = "flashed",
    sample_name: str = "sampleX",
    statistics: os.PathLike = "deduplication_stats.json",
):
    logger.info("Connecting to DuckDB")
    con = ibis.duckdb.connect()

    slices_tbl_raw = read_parquet_with_ibis(con, slices)
    
    n_slices_raw = slices_tbl_raw[["slice_id"]].distinct().count().execute()
    n_reads_raw = slices_tbl_raw[["parent_id"]].distinct().count().execute()

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

    parent_ids_unique = query.execute()

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
        max_rows_per_file=int(2e6),
    )

    # If the output directory is empty, create a dummy file to prevent downstream errors
    if not os.path.exists(output):
        os.makedirs(output)
        df_dummy = scanner.to_table().to_pandas()
        df_dummy.to_parquet(f"{output}/dummy.parquet")

    logger.info("Calculating deduplication stats")

    # Calculate the number of reads in the output
    n_reads_unique = parent_ids_unique.shape[0]
    
    # Calculate the number of slices in the output
    tbl_dedup = read_parquet_with_ibis(con, output)
    n_slices_unique = tbl_dedup[["slice_id"]].distinct().count().execute()
    

    stats = AlignmentDeduplicationStats(
        sample=sample_name,
        read_type=read_type,
        n_total_reads=n_reads_raw,
        n_unique_reads=n_reads_unique,
        n_total_slices=n_slices_raw,
        n_unique_slices=n_slices_unique,
    )
    
    with open(statistics, "w") as f:
        f.write(stats.model_dump_json())
    
