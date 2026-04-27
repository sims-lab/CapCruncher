from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from loguru import logger


@dataclass(frozen=True)
class ReporterViewpointSummary:
    viewpoints: list[str]
    viewpoint_sizes: dict[str, int]
    viewpoint_sizes_table: pd.DataFrame
    low_memory: bool
    partitions: list[str] | None


def valid_viewpoint_names(viewpoint_path: Path | str) -> list[str]:
    """Return unique viewpoint names from a BED-like viewpoint file."""
    viewpoints = pd.read_csv(
        viewpoint_path,
        sep="\t",
        header=None,
        usecols=[3],
        names=["name"],
    )
    return viewpoints["name"].dropna().astype(str).drop_duplicates().tolist()


def parquet_files(path: Path | str) -> list[Path]:
    """Return parquet files represented by a file path or directory path."""
    path = Path(path)
    if path.is_dir():
        return sorted(path.glob("*.parquet"))
    return [path]


def _normalise_nullable_viewpoints(reporters_df: pd.DataFrame) -> pd.Series:
    return reporters_df["viewpoint"].astype("string").replace(
        {"": pd.NA, "None": pd.NA, "nan": pd.NA, "<NA>": pd.NA}
    )


def _validate_reporter_columns(reporters_df: pd.DataFrame, parquet_file: Path) -> None:
    required_columns = {"viewpoint"}
    missing_columns = required_columns - set(reporters_df.columns)
    if missing_columns:
        raise ValueError(
            f"Reporter file {parquet_file} is missing required column(s): "
            f"{', '.join(sorted(missing_columns))}"
        )


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
        reporters_df = pd.read_parquet(parquet_file)
        _validate_reporter_columns(reporters_df, parquet_file)
        reporters_df["viewpoint"] = _normalise_nullable_viewpoints(reporters_df)
        invalid_viewpoints = sorted(
            set(reporters_df["viewpoint"].dropna().astype(str)) - set(valid_viewpoints)
        )
        if invalid_viewpoints:
            raise ValueError(
                "Reporter file contains viewpoint values not present in "
                f"{viewpoint_path}: {invalid_viewpoints}"
            )

        for column in reporters_df.select_dtypes(include="category").columns:
            reporters_df[column] = reporters_df[column].cat.remove_unused_categories()

        reporters_df["viewpoint"] = pd.Categorical(
            reporters_df["viewpoint"],
            categories=valid_viewpoints,
        )
        reporters_df.to_parquet(output_dir / f"part-{index}.parquet", index=False)

    return output_dir


def summarise_reporter_viewpoints(reporters: Path | str) -> ReporterViewpointSummary:
    logger.info("Extracting viewpoint names and sizes")

    viewpoint_df = pd.read_parquet(reporters, engine="pyarrow", columns=["viewpoint"])
    if hasattr(viewpoint_df["viewpoint"], "cat") and hasattr(
        viewpoint_df["viewpoint"].cat, "categories"
    ):
        viewpoints = viewpoint_df["viewpoint"].cat.categories.to_list()
    else:
        viewpoints = (
            viewpoint_df["viewpoint"].dropna().drop_duplicates().astype(str).to_list()
        )

    viewpoint_sizes = viewpoint_df["viewpoint"].value_counts()
    viewpoint_sizes_dict = viewpoint_sizes.to_dict()
    viewpoint_sizes_table = pd.DataFrame.from_dict(
        viewpoint_sizes_dict, orient="index", columns=["n_slices"]
    )

    logger.info(f"Number of viewpoints: {len(viewpoints)}")
    logger.info(f"Number of slices per viewpoint:\n{viewpoint_sizes_table}")

    if any(size > 2e6 for size in viewpoint_sizes_dict.values()):
        logger.warning(
            "High number of slices per viewpoint detected. Switching to low memory mode"
        )
        partition_df = pd.read_parquet(reporters, engine="pyarrow", columns=["bam"])
        partitions = partition_df["bam"].cat.categories.to_list()
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
