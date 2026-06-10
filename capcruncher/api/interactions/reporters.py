from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl
from loguru import logger


@dataclass(frozen=True)
class ReporterViewpointSummary:
    viewpoints: list[str]
    viewpoint_sizes: dict[str, int]
    viewpoint_sizes_table: pl.DataFrame
    low_memory: bool
    partitions: list[str] | None


def valid_viewpoint_names(viewpoint_path: Path | str) -> list[str]:
    """Return unique viewpoint names from a BED-like viewpoint file."""
    viewpoints = pl.read_csv(
        viewpoint_path,
        separator="\t",
        has_header=False,
        columns=[3],
        new_columns=["name"],
    )
    return (
        viewpoints.select(pl.col("name").drop_nulls().cast(pl.Utf8).unique())
        .get_column("name")
        .to_list()
    )


def parquet_files(path: Path | str) -> list[Path]:
    """Return parquet files represented by a file path or directory path."""
    path = Path(path)
    if path.is_dir():
        return sorted(path.glob("*.parquet"))
    return [path]


def _normalise_nullable_viewpoints(reporters_df: pl.DataFrame) -> pl.DataFrame:
    viewpoint = pl.col("viewpoint").cast(pl.Utf8)
    return reporters_df.with_columns(
        pl.when(viewpoint.is_in(["", "None", "nan", "<NA>"]))
        .then(None)
        .otherwise(viewpoint)
        .alias("viewpoint")
    )


def _validate_reporter_columns(reporters_df: pl.DataFrame, parquet_file: Path) -> None:
    required_columns = {"viewpoint"}
    missing_columns = required_columns - set(reporters_df.columns)
    if missing_columns:
        raise ValueError(
            f"Reporter file {parquet_file} is missing required column(s): "
            f"{', '.join(sorted(missing_columns))}"
        )


def _read_reporter_columns(reporters: Path | str, columns: list[str]) -> pl.DataFrame:
    frames = [
        pl.read_parquet(path, columns=columns) for path in parquet_files(reporters)
    ]
    if not frames:
        return pl.DataFrame()
    return pl.concat(frames, how="diagonal_relaxed")


def write_countable_reporters(
    reporters: Path | str, viewpoint_path: Path | str, output_dir: Path | str
) -> Path:
    """Write reporter parquet files with viewpoint categories limited to real baits.

    ``capcruncher-tools`` expects the reporter ``viewpoint`` category set to contain
    only viewpoints from the bait BED. Older CapCruncher reporter files can carry
    unused synthetic categories, so this normalises categories while still rejecting
    actual non-viewpoint values.
    """
    valid_viewpoints = valid_viewpoint_names(viewpoint_path)
    if not valid_viewpoints:
        raise ValueError(f"No viewpoints found in {viewpoint_path}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for index, parquet_file in enumerate(parquet_files(reporters)):
        reporters_df = pl.read_parquet(parquet_file)
        _validate_reporter_columns(reporters_df, parquet_file)
        reporters_df = _normalise_nullable_viewpoints(reporters_df)
        invalid_viewpoints = sorted(
            set(
                reporters_df.select(
                    pl.col("viewpoint").drop_nulls().cast(pl.Utf8).unique()
                )
                .get_column("viewpoint")
                .to_list()
            )
            - set(valid_viewpoints)
        )
        if invalid_viewpoints:
            raise ValueError(
                "Reporter file contains viewpoint values not present in "
                f"{viewpoint_path}: {invalid_viewpoints}"
            )

        reporters_df = reporters_df.with_columns(
            pl.col("viewpoint").cast(pl.Enum(valid_viewpoints))
        )
        reporters_df.write_parquet(output_dir / f"part-{index}.parquet")

    return output_dir


def summarise_reporter_viewpoints(reporters: Path | str) -> ReporterViewpointSummary:
    logger.info("Extracting viewpoint names and sizes")

    viewpoint_df = _read_reporter_columns(reporters, columns=["viewpoint"])
    viewpoint_dtype = viewpoint_df.schema["viewpoint"]
    if hasattr(viewpoint_dtype, "categories"):
        viewpoints = viewpoint_dtype.categories.to_list()
    else:
        viewpoints = (
            viewpoint_df.select(pl.col("viewpoint").drop_nulls().cast(pl.Utf8).unique())
            .get_column("viewpoint")
            .to_list()
        )

    viewpoint_sizes_table = (
        viewpoint_df.drop_nulls("viewpoint")
        .group_by("viewpoint")
        .len(name="n_slices")
        .sort("viewpoint")
    )
    viewpoint_sizes_dict = {
        row["viewpoint"]: row["n_slices"] for row in viewpoint_sizes_table.to_dicts()
    }

    logger.info(f"Number of viewpoints: {len(viewpoints)}")
    logger.info(f"Number of slices per viewpoint:\n{viewpoint_sizes_table}")

    if any(size > 2e6 for size in viewpoint_sizes_dict.values()):
        logger.warning(
            "High number of slices per viewpoint detected. Switching to low memory mode"
        )
        partition_df = _read_reporter_columns(reporters, columns=["bam"])
        partitions = (
            partition_df.select(pl.col("bam").drop_nulls().cast(pl.Utf8).unique())
            .get_column("bam")
            .to_list()
        )
        low_memory = True
    else:
        partitions = None
        low_memory = False

    return ReporterViewpointSummary(
        viewpoints=viewpoints,
        viewpoint_sizes=viewpoint_sizes_dict,
        viewpoint_sizes_table=viewpoint_sizes_table,
        low_memory=low_memory,
        partitions=partitions,
    )
