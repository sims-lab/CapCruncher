import os
import shutil
from pathlib import Path

import polars as pl
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
from loguru import logger

from capcruncher.api.statistics import AlignmentDeduplicationStats
from capcruncher.types import ReadType, VALID_READ_TYPES, validate_choice


def read_parquet(path: Path | str) -> pl.LazyFrame:
    parquet_path = str(path)
    if os.path.isdir(parquet_path):
        parquet_path = os.path.join(parquet_path, "*.parquet")

    return pl.scan_parquet(parquet_path)


def remove_unused_dictionary_values(table: pa.Table) -> pa.Table:
    for index, field in enumerate(table.schema):
        if pa.types.is_dictionary(field.type):
            column = pc.dictionary_encode(pc.cast(table.column(field.name), pa.string()))
            table = table.set_column(index, field.name, column)
    return table


def deduplicate(
    slices: Path | str,
    output: Path | str,
    read_type: ReadType | str = ReadType.FLASHED,
    sample_name: str = "sampleX",
    statistics: Path | str = "deduplication_stats.json",
) -> None:
    logger.info("Loading parquet input")
    read_type = validate_choice(read_type, VALID_READ_TYPES, "read_type")
    slices_tbl_raw = read_parquet(slices)

    n_slices_raw = (
        slices_tbl_raw.select(pl.col("slice_id").n_unique().alias("count"))
        .collect()
        .item()
    )
    n_reads_raw = (
        slices_tbl_raw.select(pl.col("parent_id").n_unique().alias("count"))
        .collect()
        .item()
    )

    if read_type == ReadType.PE:
        logger.info("Read type is PE")
        logger.info("Identifying unique fragment IDs")
        query = (
            slices_tbl_raw.select(["chrom", "start", "end", "parent_id"])
            .sort(["parent_id", "chrom", "start", "end"])
            .group_by("parent_id")
            .agg(
                slice_f_chrom=pl.col("chrom").first(),
                slice_f_start=pl.col("start").first(),
                slice_l_end=pl.col("end").last(),
            )
            .group_by(["slice_f_chrom", "slice_f_start", "slice_l_end"])
            .agg(pl.col("parent_id").first().alias("pid"))
            .select(pl.col("pid").unique())
        )
    elif read_type == ReadType.FLASHED:
        logger.info("Read type is Flashed")
        logger.info("Identifying unique fragment IDs")

        query = (
            slices_tbl_raw.select(["coordinates", "parent_id"])
            .with_columns(pl.col("coordinates").cast(pl.Utf8))
            .sort(["parent_id", "coordinates"])
            .group_by("parent_id")
            .agg(
                pl.col("coordinates")
                .sort()
                .implode()
                .list.join(",")
                .alias("coordinates")
            )
            .group_by("coordinates")
            .agg(pl.col("parent_id").first().alias("parent_id_unique"))
            .select(pl.col("parent_id_unique").unique())
        )
    else:
        raise ValueError(f"Unsupported read_type: {read_type}")

    parent_ids_unique = query.collect().to_series(0).to_list()

    logger.info("Writing deduplicated slices to disk")
    slices_unfiltered_ds = ds.dataset(slices, format="parquet")
    parent_id_type = slices_unfiltered_ds.schema.field("parent_id").type
    scanner = slices_unfiltered_ds.scanner(
        filter=ds.field("parent_id").isin(
            pa.array(parent_ids_unique, type=parent_id_type)
        )
    )

    if os.path.exists(output):
        shutil.rmtree(output)

    deduplicated_slices = remove_unused_dictionary_values(scanner.to_table())

    ds.write_dataset(
        deduplicated_slices,
        output,
        format="parquet",
        partitioning_flavor="hive",
        min_rows_per_group=0,
        max_rows_per_file=int(2e6),
    )

    # If the output directory is empty, create a dummy file to prevent downstream errors
    if not os.path.exists(output):
        os.makedirs(output)
        df_dummy = deduplicated_slices.to_pandas()
        df_dummy.to_parquet(f"{output}/dummy.parquet")

    logger.info("Calculating deduplication stats")

    # Calculate the number of reads in the output
    n_reads_unique = len(parent_ids_unique)

    # Calculate the number of slices in the output
    tbl_dedup = read_parquet(output)
    n_slices_unique = (
        tbl_dedup.select(pl.col("slice_id").n_unique().alias("count"))
        .collect()
        .item()
    )

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
