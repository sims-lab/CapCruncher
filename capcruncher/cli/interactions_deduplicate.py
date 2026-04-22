import os
import shutil

import duckdb
import pyarrow.dataset as ds
from loguru import logger

from capcruncher.api.statistics import AlignmentDeduplicationStats


def parquet_source(path: os.PathLike) -> str:
    parquet_path = f"{path}/*.parquet" if os.path.isdir(path) else str(path)
    return parquet_path.replace("'", "''")


def deduplicate(
    slices: os.PathLike,
    output: os.PathLike,
    read_type: str = "flashed",
    sample_name: str = "sampleX",
    statistics: os.PathLike = "deduplication_stats.json",
):
    logger.info("Connecting to DuckDB")
    con = duckdb.connect()
    slices_source = parquet_source(slices)

    n_slices_raw = con.sql(
        f"SELECT COUNT(DISTINCT slice_id) AS n_slices FROM read_parquet('{slices_source}')"
    ).fetchone()[0]
    n_reads_raw = con.sql(
        f"SELECT COUNT(DISTINCT parent_id) AS n_reads FROM read_parquet('{slices_source}')"
    ).fetchone()[0]

    if read_type == "pe":
        logger.info("Read type is PE")
        logger.info("Identifying unique fragment IDs")
        parent_ids_unique = con.sql(
            f"""
            WITH ranked AS (
                SELECT
                    chrom,
                    start,
                    "end",
                    parent_id,
                    row_number() OVER (
                        PARTITION BY parent_id
                        ORDER BY chrom, start, "end", parent_id
                    ) AS rn_first,
                    row_number() OVER (
                        PARTITION BY parent_id
                        ORDER BY chrom DESC, start DESC, "end" DESC, parent_id DESC
                    ) AS rn_last
                FROM read_parquet('{slices_source}')
            ),
            parent_ranges AS (
                SELECT
                    parent_id,
                    max(CASE WHEN rn_first = 1 THEN chrom END) AS slice_f_chrom,
                    max(CASE WHEN rn_first = 1 THEN start END) AS slice_f_start,
                    max(CASE WHEN rn_last = 1 THEN "end" END) AS slice_l_end
                FROM ranked
                GROUP BY parent_id
            )
            SELECT min(parent_id) AS parent_id_unique
            FROM parent_ranges
            GROUP BY slice_f_chrom, slice_f_start, slice_l_end
            """
        ).df()["parent_id_unique"].to_numpy()
    elif read_type == "flashed":
        logger.info("Read type is Flashed")
        logger.info("Identifying unique fragment IDs")
        parent_ids_unique = con.sql(
            f"""
            WITH parent_coordinates AS (
                SELECT
                    parent_id,
                    string_agg(coordinates, ',' ORDER BY coordinates) AS coordinates
                FROM read_parquet('{slices_source}')
                GROUP BY parent_id
            )
            SELECT min(parent_id) AS parent_id_unique
            FROM parent_coordinates
            GROUP BY coordinates
            """
        ).df()["parent_id_unique"].to_numpy()

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
    output_source = parquet_source(output)
    n_slices_unique = con.sql(
        f"SELECT COUNT(DISTINCT slice_id) AS n_slices FROM read_parquet('{output_source}')"
    ).fetchone()[0]
    

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
    
